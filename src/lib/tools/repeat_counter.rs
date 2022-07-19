use std::fs::File;
use std::{path::PathBuf, str::FromStr};

use anyhow::Context;
use anyhow::Result;
use clap::Parser;
use env_logger::Env;
use fgoxide::io::Io;

use crate::utils::built_info;
use noodles::bam;
use noodles::bed;
use noodles::bgzf;
use noodles::core as noodles_core;
use noodles::fasta;
use noodles::sam;
use noodles::sam::record::sequence::Base;

use flume::{bounded, unbounded, Receiver, Sender};
use noodles::core::Position;

/// Counts repeats in your BAM given a BED of repeats
#[derive(Parser, Debug)]
#[clap(name = "repeat-finder", verbatim_doc_comment, version = built_info::VERSION.as_str())]
pub struct Opts {
    /// The input BAM for the normal sample
    #[clap(short = 'n', long, display_order = 1)]
    pub normal: PathBuf,

    // /// The input BAM for the tumor sample
    // #[clap(short = 't', long, display_order = 2)]
    // pub tumor: PathBuf,
    /// The input BED file with MSI loci
    #[clap(short = 'b', long, display_order = 3)]
    pub bedfile: PathBuf,

    // /// The input reference genome FASTA
    // #[clap(short = 'g', long, display_order = 4)]
    // pub genome: PathBuf,

    // // /// The path to the output
    // // #[clap(short = 'o', long, display_order = 5)]
    // // pub output: PathBuf,

    // // /// The minimum average per-base read quality for a read to pass the quality control filters
    // // #[clap(long, default_value = "25", display_order = 6)]
    // // pub min_read_mean_base_quality: u64,

    // // /// Minimum average per-base quality for the bases contained within the microsatellite locus.
    // // /// Reads that pass the read quality filter (above) will still fail quality control if the
    // // /// locus quality scores are too low.
    // // #[clap(long, default_value = "30", display_order = 7)]
    // // pub min_locus_mean_base_quality: u64,

    // // /// Minimum read length for a read to pass quality control. Only bases that are not clipped
    // // /// will be considered; in other words, soft-clipped or hard-clipped parts of the read do not
    // // /// count towards the length.
    // // #[clap(long, default_value = "35", display_order = 7)]
    // // pub min_read_length: u64,
    /// The number of threads to use
    #[clap(long, default_value = "1", display_order = 3)]
    pub threads: u64,
}

struct BamRecordJob {
    locus: bed::Record<4>,
    record: sam::alignment::Record,
    is_normal: bool,
}

fn get_index(bam: &PathBuf) -> Result<bam::bai::Index> {
    let bam_bai = bam.with_extension("bam.bai");
    let index = if bam_bai.exists() {
        bam::bai::read(bam_bai)
    } else {
        bam::bai::read(bam.with_extension("bai"))
    };
    index.with_context(|| format!("Could not open BAM index for BAM: {:?}", bam))
}

fn query_reads(
    bam_header: &sam::Header,
    bam_reader: &mut bam::Reader<bgzf::Reader<File>>,
    bam_index: &bam::bai::Index,
    is_normal: bool,
    bed_record: &bed::Record<4>,
    region: &noodles_core::Region,
    bam_record_job_tx: &Sender<BamRecordJob>,
) {
    let name = if is_normal { "normal" } else { "tumor" };

    let query = bam_reader
        .query(bam_header.reference_sequences(), bam_index, &region)
        .with_context(|| format!("Could not query {} BAM for region: {:?}", name, region))
        .unwrap();
    for bam_result in query {
        let bam_record = bam_result
            .with_context(|| format!("Could not parse {} reads in region: {:?}", name, region))
            .unwrap();
        if !bam_record.flags().is_unmapped()
            && !bam_record.cigar().is_empty()
            && bam_record.alignment_start().unwrap() <= bed_record.start_position()
            && bed_record.end_position() <= bam_record.alignment_end().unwrap()
        {
            let job = BamRecordJob { locus: bed_record.clone(), record: bam_record, is_normal };
            bam_record_job_tx.send(job).with_context(|| "Could not send BamRecordJob").unwrap();
        }
    }
}

/// Convenience type for functions that return [`PoolError`].
type JobResult<T> = Result<T, anyhow::Error>;

// Run extract
#[allow(clippy::too_many_lines)]
pub fn run(opts: &Opts) -> Result<(), anyhow::Error> {
    // Open the input and output
    let io = Io::default();
    let mut bed_reader = io
        .new_reader(&opts.bedfile)
        .map(bed::Reader::new)
        .with_context(|| format!("Could not open BED for reading: {:?}", opts.bedfile))?;
    let mut normal_bam_reader = File::open(&opts.normal)
        .map(bam::Reader::new)
        .with_context(|| format!("Could not open normal BAM for reading: {:?}", opts.normal))?;
    // let mut tumor_bam_reader = File::open(&opts.tumor)
    //     .map(bam::Reader::new)
    //     .with_context(|| format!("Could not open tumor BAM for reading: {:?}", opts.tumor))?;

    // Read the BAM indexes
    let normal_bai = get_index(&opts.normal)?;
    // let tumor_bai = get_index(&opts.tumor)?;

    // Read the SAM headers
    let normal_bam_header: sam::Header = normal_bam_reader.read_header()?.parse()?;
    // let tumor_bam_header: sam::Header = tumor_bam_reader.read_header()?.parse()?;

    let (bam_record_job_tx, bam_record_job_rx): (Sender<BamRecordJob>, Receiver<BamRecordJob>) =
        flume::unbounded();

    let workers: Vec<std::thread::JoinHandle<JobResult<()>>> = (0..opts.threads)
        .map(|_i| {
            let rx = bam_record_job_rx.clone();
            std::thread::spawn(move || {
                while let Ok(job) = rx.recv() {
                    // TODO
                    println!(
                        "record: {:?}",
                        String::from_utf8_lossy(job.record.read_name().unwrap().as_ref())
                    );
                }
                Ok(())
            })
        })
        .collect();

    // Iterate loci by loci
    for (index, bed_result) in bed_reader.records::<4>().enumerate() {
        let bed_record = bed_result
            .with_context(|| format!("Could not parse the {}th BED record", index + 1))?;

        let region = noodles_core::Region::new(
            bed_record.reference_sequence_name(),
            bed_record.start_position()..=bed_record.end_position(),
        );

        query_reads(
            &normal_bam_header,
            &mut normal_bam_reader,
            &normal_bai,
            true,
            &bed_record,
            &region,
            &bam_record_job_tx,
        );
        // query_reads(
        //     &tumor_bam_header,
        //     &mut tumor_bam_reader,
        //     &tumor_bai,
        //     true,
        //     &bed_record,
        //     &region,
        //     &bam_record_job_tx,
        // );
    }

    Ok(())
}
