use std::{path::PathBuf, str::FromStr};

use anyhow::Context;
use anyhow::Result;
use clap::Parser;
use env_logger::Env;
use fgoxide::io::Io;

use crate::utils::built_info;
use noodles::bed;
use noodles::fasta;

use noodles::core::Position;

/// Find repeats in your reference FASTA
#[derive(Parser, Debug)]
#[clap(name = "repeat-finder", verbatim_doc_comment, version = built_info::VERSION.as_str())]
pub struct Opts {
    /// The input reference genome FASTA
    #[clap(short = 'i', long, display_order = 1)]
    pub input: PathBuf,

    /// The output BED file with MSI loci
    #[clap(short = 'o', long, display_order = 2)]
    pub output: PathBuf,

    /// The minimum number of bases for a repeat region must span to be called a microsatellite
    #[clap(short = 'm', default_value = "10", long, display_order = 3)]
    pub min_bases: usize,

    /// The maximum number of bases for a repeat region must span to be called a microsatellite
    #[clap(short = 'M', default_value = "100", long, display_order = 4)]
    pub max_bases: usize,

    /// The minimum number of repeats for a repeat region must span to be called a microsatellite
    #[clap(short = 'r', default_value = "3", long, display_order = 5)]
    pub min_repeats: usize,

    /// The minimum k-mer length
    #[clap(short = 'l', default_value = "1", long, display_order = 6)]
    pub min_repeat_length: usize,

    /// The maximum k-mer length
    #[clap(short = 'L', default_value = "5", long, display_order = 7)]
    pub max_repeat_length: usize,
}

// Run extract
#[allow(clippy::too_many_lines)]
pub fn run(opts: &Opts) -> Result<(), anyhow::Error> {
    // Build the finders
    let mut finders: Vec<KmerFinder> =
        Vec::with_capacity(opts.max_repeat_length - opts.min_repeat_length + 1);
    let mut i = opts.min_repeat_length;
    while i <= opts.max_repeat_length {
        finders.push(KmerFinder::new(i));
        i += 1;
    }

    // Open the input and output
    let io = Io::default();

    let mut reader = io
        .new_reader(&opts.input)
        .map(fasta::Reader::new)
        .with_context(|| format!("Could not open FASTA for reading: {:?}", opts.input))?;
    let mut writer = io
        .new_writer(&opts.output)
        .map(bed::Writer::new)
        .with_context(|| format!("Could not open BED for writing: {:?}", opts.output))?;

    // Go through each contig, one at a time
    for (index, result) in reader.records().enumerate() {
        let record =
            result.with_context(|| format!("Could not parse the {}th FASTA record", index + 1))?;
        let contig = record.name();

        for finder in &mut finders {
            finder.reset();
        }

        // Go through each base, one at a time
        let mut i = 1;
        while i <= record.sequence().len() {
            let position = Position::try_from(i)
                .with_context(|| format!("Could not create a Position from {}", i))?;
            let base: u8 = record
                .sequence()
                .get(position)
                .with_context(|| format!("Could not retrieve base at {}:{}", record.name(), i))
                .unwrap()
                & 0xdf; // to upper case

            // Output the smallest repeat found ending at this position.
            let mut found = false;
            for finder in &mut finders {
                if let Some(repeat) = finder.add_maybe_emit(base, !found) {
                    if found {
                        continue;
                    } else if let Some(rec) = to_bed_record(opts, contig, i, finder, &repeat) {
                        writer
                            .write_record(&rec)
                            .with_context(|| format!("Could not write BED record {:?}", rec))
                            .unwrap();
                        found = true;
                    }
                }
            }
            i += 1;
        }

        // Emit any repeat that goes to the end of the contig
        for finder in &mut finders {
            if let Some(repeat) = finder.emit() {
                // only output the first repeat found at this position
                if let Some(rec) = to_bed_record(opts, contig, i, finder, &repeat) {
                    writer
                        .write_record(&rec)
                        .with_context(|| format!("Could not write BED record {:?}", rec))
                        .unwrap();
                    break;
                }
            }
        }
    }

    Ok(())
}

/// Converts the given repeat to a [`bed::Record`] if the repeat passes all the filters, otherwise
/// None.
fn to_bed_record(
    opts: &Opts,
    contig: &str,
    position: usize,
    finder: &KmerFinder,
    repeat: &Repeat,
) -> Option<bed::Record<4>> {
    if repeat.span < opts.min_bases
        || opts.max_bases < repeat.span
        || repeat.span / finder.len() < opts.min_repeats
        || is_repeat(&repeat.kmer)
    {
        None
    } else {
        let unit = String::from_utf8_lossy(&repeat.kmer);
        let name = format!("({}){:?}", unit, repeat.span / finder.len());
        let start_position = Position::try_from(position - repeat.span)
            .with_context(|| {
                format!("Could not create BED start Position: {}", position - repeat.span)
            })
            .unwrap();
        let end_position = Position::try_from(position - 1)
            .with_context(|| format!("Could not create BED end Position: {}", position - 1))
            .unwrap();
        let name = bed::record::Name::from_str(&name)
            .with_context(|| format!("Could not create BED name: {}", name))
            .unwrap();
        let bed_record = bed::Record::<4>::builder()
            .set_reference_sequence_name(contig)
            .set_start_position(start_position)
            .set_end_position(end_position)
            .set_name(name)
            .build()
            .with_context(|| {
                format!(
                    "Could not build a BED record at {}:{} for repeat {:?}",
                    contig, position, repeat
                )
            })
            .unwrap();
        Some(bed_record)
    }
}

// TODO: caching this may help speed things up
/// Returns true if the given sequence is composed of a repeated smaller sequence.
fn is_repeat(unit: &[u8]) -> bool {
    let mut kmer_size = unit.len() / 2;
    if kmer_size == unit.len() {
        kmer_size = unit.len() - 1;
    }
    while 0 < kmer_size {
        if unit.len() % kmer_size == 0 {
            let mut i = 0;
            let mut j = kmer_size;
            while j < unit.len() && unit[i] == unit[j] {
                i += 1;
                j += 1;
            }
            if j == unit.len() {
                return true;
            }
        }
        kmer_size -= 1;
    }

    false
}

/// Parse args and set up logging / tracing
pub fn setup() -> Opts {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    Opts::parse()
}

/// Stores information about a repeat that has been found.  T
#[derive(Debug, Clone, PartialEq, Hash, Eq)]
pub struct Repeat {
    /// The bases of the repeat
    pub kmer: Vec<u8>,
    /// The number of times the repeat is repeated
    pub num_repeats: usize,
    /// The span of the repeat.  The repeat may span a few extra bases if the kmer is not
    /// repeated exactly.  Eg. a kmer of ACG could be seen twice ACGACGA so has span 7.
    pub span: usize,
}

/// Struct to help find repeats of a given size when a contiguous sequence is provided one base
/// at a time.
#[derive(Debug, Clone, PartialEq, Hash, Eq)]
pub struct KmerFinder {
    // The current set of bases seen of a given size
    pub kmer: Vec<u8>,
    // The index of the next base
    index: usize,
    // The span of the current repeat
    span: usize,
}

impl KmerFinder {
    fn new(size: usize) -> KmerFinder {
        KmerFinder { kmer: vec![b'n'; size], index: 0, span: 0 }
    }

    fn len(&self) -> usize {
        self.kmer.len()
    }

    /// Adds a new base, and optionally returns a repeat if found.
    ///
    /// Conceptually for a repeat with unit length N, stores up to the last N bases.  For each new
    /// base provided to [`KmerFinder::add`], if the span is smaller than the unit length the next
    /// base is updated, otherwise the new base is compared to the previous base at the appropriate
    /// offset in the repeat.  In the latter case, if the base matches, the span is incremented and
    /// None is returned, otherwise the repeat up to this point is returned and the kmer finder
    /// reset to a new kmer that includes the new base as the last base.
    fn add_maybe_emit(&mut self, base: u8, emit: bool) -> Option<Repeat> {
        let retval = if self.span < self.len() {
            // Not enough bases added yet, so update the kmer
            self.kmer[self.index] = base;
            None
        } else if self.kmer[self.index] == base {
            // The given base matches the expected base
            None
        } else {
            // The given base does not match the expected base, so return the current repeat (if
            // long enough).  Reset the span to the unit length to start a new repeat.
            let repeat = if emit { self.emit() } else { None };
            self.span = self.len() - 1; // span is incremented below, so subtract one here
            self.kmer[self.index] = base;
            repeat
        };
        // increment index and span
        self.index = if self.index == self.len() - 1 { 0 } else { self.index + 1 };
        self.span += 1;
        retval
    }

    #[cfg(test)]
    fn add(&mut self, base: u8) -> Option<Repeat> {
        self.add_maybe_emit(base, true)
    }

    /// Emit the current repeat seen so far.  Returns None if not enough bases have been added.
    fn emit(&self) -> Option<Repeat> {
        if self.span < self.len() {
            None
        } else {
            let mut kmer = vec![b'n'; self.len()];
            let mut i = 0;
            let mut j = {
                let remainder = self.span % self.len();
                if remainder <= self.index {
                    self.index - remainder
                } else {
                    self.index + self.len() - remainder
                }
            };
            while i < self.len() {
                kmer[i] = self.kmer[j];
                i += 1;
                j = if j == self.len() - 1 { 0 } else { j + 1 };
            }
            let num_repeats = self.span / self.len();

            Some(Repeat { kmer, num_repeats, span: self.span })
        }
    }

    /// Resets the repeat finder, for example when examining a new contig
    fn reset(&mut self) {
        self.kmer = vec![b'n'; self.len()];
        self.index = 0;
        self.span = 0;
    }
}

#[cfg(test)]
mod test {
    use super::is_repeat;
    use super::run;
    use super::KmerFinder;
    use super::Opts;
    use super::Repeat;
    use noodles::bed;
    use noodles::fasta;
    use rstest::rstest;
    use std::{
        fs::File,
        io::{BufReader, BufWriter},
    };
    use tempfile::tempdir;

    /// Helper function to verify values in [`Repeat`]
    fn check_repeat(repeat: Repeat, kmer: &str, num_repeats: usize, span: usize) {
        assert_eq!(String::from_utf8(repeat.kmer).unwrap(), kmer);
        assert_eq!(repeat.num_repeats, num_repeats);
        assert_eq!(repeat.span, span);
    }

    #[rstest]
    #[case("A", false)]
    #[case("AA", true)]
    #[case("CGCG", true)]
    #[case("CAGCAGCA", false)]
    #[case("CAGCAGCAG", true)]
    #[case("CAGCTGCAG", false)]
    #[case("CAGCAGCAGCAGCAG", true)]
    #[case("CAGCAGCAGCAGCAGC", false)]
    #[case("CAGCAGCAGCAGCAGCA", false)]
    fn test_is_repeat(#[case] unit: &str, #[case] expected_value: bool) {
        assert_eq!(is_repeat(unit.as_bytes()), expected_value);
    }

    #[test]
    fn test_emit() {
        // too few bases seen to emit a repeat
        assert_eq!(KmerFinder::new(3).emit(), None);
        assert_eq!(
            KmerFinder { kmer: [b'a', b'c', b'g'].to_vec(), index: 1, span: 1 }.emit(),
            None
        );
        assert_eq!(
            KmerFinder { kmer: [b'a', b'c', b'g'].to_vec(), index: 2, span: 2 }.emit(),
            None
        );

        // emit repeats
        check_repeat(
            KmerFinder { kmer: [b'a', b'c', b'g'].to_vec(), index: 0, span: 3 }.emit().unwrap(),
            "acg",
            1,
            3,
        );
        check_repeat(
            KmerFinder { kmer: [b'a', b'c', b'g'].to_vec(), index: 1, span: 4 }.emit().unwrap(),
            "acg",
            1,
            4,
        );
        check_repeat(
            KmerFinder { kmer: [b'a', b'c', b'g'].to_vec(), index: 2, span: 5 }.emit().unwrap(),
            "acg",
            1,
            5,
        );
        check_repeat(
            KmerFinder { kmer: [b'a', b'c', b'g'].to_vec(), index: 0, span: 6 }.emit().unwrap(),
            "acg",
            2,
            6,
        );

        // emit repeats repeats with various indexes, but the same span
        check_repeat(
            KmerFinder { kmer: [b'a', b'c', b'g'].to_vec(), index: 0, span: 5 }.emit().unwrap(),
            "cga",
            1,
            5,
        );
        check_repeat(
            KmerFinder { kmer: [b'a', b'c', b'g'].to_vec(), index: 1, span: 5 }.emit().unwrap(),
            "gac",
            1,
            5,
        );
        check_repeat(
            KmerFinder { kmer: [b'a', b'c', b'g'].to_vec(), index: 2, span: 5 }.emit().unwrap(),
            "acg",
            1,
            5,
        );
    }

    #[rstest]
    #[case(1)]
    #[case(2)]
    #[case(3)]
    #[case(4)]
    #[case(5)]
    #[case(6)]
    #[case(7)]
    #[case(8)]
    #[case(9)]
    fn test_simple_add(#[case] kmer_len: usize) {
        let mut finder = KmerFinder::new(kmer_len);

        let mut i = 0;
        while i < kmer_len {
            assert_eq!(finder.add(b'a'), None);
            if i < kmer_len - 1 {
                assert_eq!(finder.emit(), None);
            }
            i += 1;
        }
        check_repeat(finder.emit().unwrap(), &"a".repeat(kmer_len), 1, kmer_len);
    }

    /// Tests sequential calls to add looking for a trinucleotide repeat.  Also checks emit to
    /// verify the in progress repeat.
    #[test]
    fn test_trinuc() {
        let mut finder = KmerFinder::new(3);

        // add the tri-nuc
        assert_eq!(finder.add(b'a'), None);
        assert_eq!(finder.emit(), None);
        assert_eq!(finder.add(b'c'), None);
        assert_eq!(finder.emit(), None);
        assert_eq!(finder.add(b'g'), None);
        check_repeat(finder.emit().unwrap(), "acg", 1, 3);

        // add a few more copies
        let mut num_repeats = 1;
        let mut span = 3;
        while num_repeats <= 10 {
            assert_eq!(finder.add(b'a'), None);
            span += 1;
            check_repeat(finder.emit().unwrap(), "acg", num_repeats, span);
            assert_eq!(finder.add(b'c'), None);
            span += 1;
            check_repeat(finder.emit().unwrap(), "acg", num_repeats, span);
            assert_eq!(finder.add(b'g'), None);
            num_repeats += 1;
            span += 1;
            check_repeat(finder.emit().unwrap(), "acg", num_repeats, span);
        }

        // add one more base, so not a fully new trinuc repeat
        assert_eq!(finder.add(b'a'), None);
        span += 1;
        check_repeat(finder.emit().unwrap(), "acg", num_repeats, span);

        // add a mismatching base, which should yields the previous repeat
        check_repeat(finder.add(b't').unwrap(), "acg", num_repeats, span);
        // emitting the current repeat, yields a new one, with one copy (atg tga gat)
        check_repeat(finder.emit().unwrap(), "gat", 1, 3);
    }

    /// Helper method to create a new [`fasta::Record`]
    fn to_contig(name: &str, sequence: &str, writer: &mut fasta::Writer<BufWriter<File>>) {
        let definition = fasta::record::Definition::new(name, None);
        let sequence = fasta::record::Sequence::from(sequence.as_bytes().to_vec());
        let record = fasta::Record::new(definition, sequence);
        writer.write_record(&record).unwrap();
    }

    /// Helper method to check values in a [`bed::Record`]
    fn check_bed(record: &bed::Record<4>, contig: &str, start: usize, end: usize, name: &str) {
        assert_eq!(record.reference_sequence_name(), contig);
        assert_eq!(usize::from(record.start_position()), start);
        assert_eq!(usize::from(record.end_position()), end);
        assert_eq!(record.name().unwrap().to_string(), name);
    }

    #[test]
    fn test_repeat_finder() {
        let tempdir = tempdir().unwrap();
        let in_fasta = tempdir.path().join("input.fasta");
        let out_bed = tempdir.path().join("output.bed");
        let opts = Opts {
            input: in_fasta.clone(),
            output: out_bed.clone(),
            min_bases: 10,
            max_bases: 20,
            min_repeats: 3,
            min_repeat_length: 2,
            max_repeat_length: 5,
        };

        // Write a fasta
        {
            let mut fasta_writer: fasta::Writer<BufWriter<File>> =
                File::create(in_fasta).map(BufWriter::new).map(fasta::Writer::new).unwrap();

            // add contigs that will have no repeats returned
            to_contig("too_few_bases", &"ACG".repeat(3), &mut fasta_writer);
            to_contig("too_many_bases", &"ACG".repeat(34), &mut fasta_writer);
            to_contig("too_few_repeats", &"AGGAT".repeat(2), &mut fasta_writer);
            to_contig("too_small_repeat_length", &"A".repeat(15), &mut fasta_writer);
            to_contig("too_large_repeat_length", &"ACGTGA".repeat(2), &mut fasta_writer);

            // // add contigs with repeats on the parameter boundaries
            to_contig("(CG)5", &"CG".repeat(5), &mut fasta_writer); // min_bases and min_repeat_length
            to_contig("(CG)10", &"CG".repeat(10), &mut fasta_writer); // max_bases
            to_contig("(AGGAT)3", &"AGGAT".repeat(3), &mut fasta_writer); // min_repeats and max_repeat_length

            // add a complicate contig
            let mut sequence = "CG".repeat(5); // repeat (CG)5
            sequence.push('C'); // extra base, but still only 5 repeats of CG
            sequence.push_str(&"CG".repeat(10)); // repeat (CG)10
            sequence.push_str("TCGTT"); // non-sense
            sequence.push_str(&"ACG".repeat(3)); // too few bases
            sequence.push_str(&"TGGAT".repeat(3)); // repeat (TGGAT)3
            sequence.push_str("AGG"); // extra bases, but still only 3 repeats of AGGAT
            to_contig("complicated", &sequence, &mut fasta_writer);
        }

        run(&opts).unwrap();

        // Read in the output bED
        let mut bed_reader = File::open(out_bed).map(BufReader::new).map(bed::Reader::new).unwrap();
        let records: Vec<bed::Record<4>> =
            bed_reader.records::<4>().map(std::result::Result::unwrap).into_iter().collect();

        // Check the records
        assert_eq!(records.len(), 6);
        check_bed(&records[0], "(CG)5", 1, 10, "(CG)5");
        check_bed(&records[1], "(CG)10", 1, 20, "(CG)10");
        check_bed(&records[2], "(AGGAT)3", 1, 15, "(AGGAT)3");
        check_bed(&records[3], "complicated", 1, 11, "(CG)5");
        check_bed(&records[4], "complicated", 12, 31, "(CG)10");
        check_bed(&records[5], "complicated", 46, 60, "(TGGAT)3");
    }
}
