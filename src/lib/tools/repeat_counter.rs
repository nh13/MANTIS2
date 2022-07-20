use log::{debug, info};
use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::File;
use std::hash::Hash;
use std::thread::JoinHandle;
use std::{path::PathBuf, str::FromStr};

use anyhow::Context;
use anyhow::Result;
use clap::Parser;
use env_logger::Env;
use fgoxide::io::Io;
use noodles::bam::bai::Index;

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
#[derive(Parser, Debug, Clone)]
#[clap(name = "repeat-finder", verbatim_doc_comment, version = built_info::VERSION.as_str())]
pub struct Opts {
    /// The input BAM for the normal sample
    #[clap(short = 'n', long, display_order = 1)]
    pub normal: PathBuf,

    /// The input BAM for the tumor sample
    #[clap(short = 't', long, display_order = 2)]
    pub tumor: PathBuf,

    /// The input BED file with MSI loci
    #[clap(short = 'b', long, display_order = 3)]
    pub bedfile: PathBuf,

    /// The input reference genome FASTA
    #[clap(short = 'g', long, display_order = 4)]
    pub genome: PathBuf,

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
    pub threads: usize,
}

struct WorkerJob {
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
    bam_record_job_tx: &Sender<WorkerJob>,
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
            let job = WorkerJob { locus: bed_record.clone(), record: bam_record, is_normal };
            bam_record_job_tx.send(job).with_context(|| "Could not send BamRecordJob").unwrap();
        }
    }
}

type AnyhowResult<T> = Result<T, anyhow::Error>;

struct CounterJob {
    locus: bed::Record<4>,
    count: usize,
    is_normal: bool,
}

struct BamQueryJob {
    locus: bed::Record<4>,
    bam: PathBuf,
    bai: bam::bai::Index,
    header: sam::Header,
    is_normal: bool,
}

struct BaiAndHeader {
    bai: bam::bai::Index,
    header: sam::Header,
}

impl BaiAndHeader {
    fn new(io: &Io, path: &PathBuf, is_normal: bool) -> BaiAndHeader {
        let name = if is_normal { "normal" } else { "tumor" };
        let mut reader = File::open(path)
            .map(bam::Reader::new)
            .with_context(|| format!("Could not open {} BAM for reading: {:?}", name, path))
            .unwrap();
        let bai = get_index(&path).unwrap();
        let header: sam::Header = reader.read_header().unwrap().parse().unwrap();
        BaiAndHeader { bai, header }
    }
}

fn build_contig_to_offset(fasta: &PathBuf) -> HashMap<String, u64> {
    let fai_path: PathBuf = {
        let mut os_string = OsString::from(fasta);
        os_string.push(".");
        os_string.push("fai");
        PathBuf::from(os_string)
    };
    let fai = fasta::fai::read(fai_path.clone())
        .with_context(|| format!("Could not open FASTA index: {:?}", fai_path))
        .unwrap();
    fai.iter().map(|rec| (rec.name().to_string(), rec.offset())).collect()
}

// Run extract
#[allow(clippy::too_many_lines)]
pub fn run(opts: &Opts) -> Result<(), anyhow::Error> {
    let (bam_query_job_tx, bam_query_job_rx): (Sender<BamQueryJob>, Receiver<BamQueryJob>) =
        flume::bounded(opts.threads * 8);

    let (worker_job_tx, worker_job_rx): (Sender<WorkerJob>, Receiver<WorkerJob>) =
        flume::bounded(opts.threads * 1024 * 1024);

    let (counter_job_tx, counter_job_rx): (Sender<CounterJob>, Receiver<CounterJob>) =
        flume::bounded(opts.threads * 1024 * 1024);

    let contig_to_offset = build_contig_to_offset(&opts.genome);
    let counter_handle = std::thread::spawn(move || {
        info!("Computing counts by locus");
        let mut key_to_locus: HashMap<String, bed::Record<4>> = HashMap::new();
        let mut counter: HashMap<String, HashMap<usize, HashMap<bool, usize>>> = HashMap::new();
        while let Ok(count_job) = counter_job_rx.recv() {
            // Add to the key -> locus map
            let locus_key = count_job.locus.to_string();
            key_to_locus.entry(locus_key.clone()).or_insert(count_job.locus);

            // Get the map from repeat count to T/N count map
            let repeat = counter.entry(locus_key).or_insert_with(HashMap::new);

            // Get the map from is_normal to count
            let repeat_count = repeat.entry(count_job.count).or_insert_with(HashMap::new);

            // Increment the count for the repeat count at the given locus
            *repeat_count.entry(count_job.is_normal).or_insert(0) += 1;
        }

        info!("Collecting counts by locus");
        let mut loci: Vec<(&String, &bed::Record<4>)> = counter
            .keys()
            .map(|key| {
                let locus = match key_to_locus.get(key) {
                    Some(l) => l,
                    None => panic!("Bug: locus not found {}", key),
                };
                (key, locus)
            })
            .collect();

        info!("Sorting counts by locus");
        loci.sort_by_cached_key(|(key, locus)| {
            let offset = match contig_to_offset.get(locus.reference_sequence_name()) {
                Some(offset) => offset,
                None => panic!(
                    "Cannot find BED contig in FASTA index: {}",
                    locus.reference_sequence_name()
                ),
            };
            (offset, locus.start_position(), locus.end_position())
        });

        info!("Outputting counts by locus");
        for (key, locus) in loci {
            let repeat = match counter.get(key) {
                Some(r) => r,
                None => panic!("Bug: locus not found {}", key),
            };

            let mut lengths: Vec<&usize> = repeat.keys().collect();
            lengths.sort();

            for length in lengths {
                let repeat_count = match repeat.get(length) {
                    Some(rc) => rc,
                    None => panic!("Bug: Cannot find repeat length: {}", length),
                };
                let normal_count = repeat_count.get(&true).unwrap_or(&0);
                let tumor_count = repeat_count.get(&false).unwrap_or(&0);
                println!(
                    "{}:{}-{}\t{}\t{}\t{}",
                    locus.reference_sequence_name(),
                    locus.start_position(),
                    locus.end_position(),
                    length,
                    normal_count,
                    tumor_count
                );
            }
        }
    });

    let worker_handles: Vec<std::thread::JoinHandle<AnyhowResult<()>>> = (0..opts.threads)
        .map(|_i| {
            let rx = worker_job_rx.clone();
            let tx = counter_job_tx.clone();
            std::thread::spawn(move || {
                while let Ok(worker_job) = rx.recv() {
                    let count_job = CounterJob {
                        locus: worker_job.locus,
                        count: 0,
                        is_normal: worker_job.is_normal,
                    };
                    tx.send(count_job).unwrap();
                }
                drop(tx);
                Ok(())
            })
        })
        // Collect is needed to force the evaluation of the closure and start the loops
        .collect();

    let bam_query_handles: Vec<std::thread::JoinHandle<AnyhowResult<()>>> = (0..opts.threads)
        .map(|_i| {
            let rx = bam_query_job_rx.clone();
            let tx = worker_job_tx.clone();
            std::thread::spawn(move || {
                while let Ok(bam_query_job) = rx.recv() {
                    let name = if bam_query_job.is_normal { "normal" } else { "tumor" };
                    let locus = bam_query_job.locus;
                    let region = noodles_core::Region::new(
                        locus.reference_sequence_name(),
                        locus.start_position()..=locus.end_position(),
                    );

                    let mut reader = File::open(&bam_query_job.bam)
                        .map(bam::Reader::new)
                        .with_context(|| {
                            format!(
                                "Could not open {} BAM for reading: {:?}",
                                name, bam_query_job.bam
                            )
                        })
                        .unwrap();

                    query_reads(
                        &bam_query_job.header,
                        &mut reader,
                        &bam_query_job.bai,
                        bam_query_job.is_normal,
                        &locus,
                        &region,
                        &tx,
                    );
                }
                drop(tx);
                Ok(())
            })
        })
        // Collect is needed to force the evaluation of the closure and start the loops
        .collect();

    // Read in the BED, and send jobs to query the BAM
    let io = Io::default();
    let mut bed_reader = io
        .new_reader(&opts.bedfile)
        .map(bed::Reader::new)
        .with_context(|| format!("Could not open BED for reading: {:?}", opts.bedfile))
        .unwrap();
    let normal_bai_and_header = BaiAndHeader::new(&io, &opts.normal, true);
    let tumor_bai_and_header = BaiAndHeader::new(&io, &opts.tumor, false);

    let contig_to_offset = build_contig_to_offset(&opts.genome);
    for (index, bed_result) in bed_reader.records::<4>().enumerate() {
        let locus = bed_result
            .with_context(|| format!("Could not parse the {}th BED record", index + 1))
            .unwrap();

        assert!(
            contig_to_offset.get(locus.reference_sequence_name()) != None,
            "Cannot find BED contig in FASTA index: {}",
            locus.reference_sequence_name()
        );

        if (index + 1) % 1000 == 0 {
            info!(
                "Processed {} loci; last: {}:{}-{}",
                index + 1,
                locus.reference_sequence_name(),
                locus.start_position(),
                locus.end_position()
            );
        }

        let normal_job = BamQueryJob {
            locus: locus.clone(),
            bam: opts.normal.clone(),
            bai: normal_bai_and_header.bai.clone(),
            header: normal_bai_and_header.header.clone(),
            is_normal: true,
        };
        let tumor_job = BamQueryJob {
            locus: locus.clone(),
            bam: opts.tumor.clone(),
            bai: tumor_bai_and_header.bai.clone(),
            header: tumor_bai_and_header.header.clone(),
            is_normal: false,
        };

        bam_query_job_tx.send(normal_job)?;
        bam_query_job_tx.send(tumor_job)?;
    }
    drop(bam_query_job_tx);

    // Close the worker handles
    bam_query_handles.into_iter().try_for_each(|handle| match handle.join() {
        Ok(result) => result,
        Err(e) => std::panic::resume_unwind(e),
    })?;
    drop(worker_job_tx);

    // Close the worker handles
    worker_handles.into_iter().try_for_each(|handle| match handle.join() {
        Ok(result) => result,
        Err(e) => std::panic::resume_unwind(e),
    })?;
    drop(counter_job_tx);

    // Close the counter handle
    match counter_handle.join() {
        Ok(result) => result,
        Err(e) => std::panic::resume_unwind(e),
    };

    Ok(())
}
