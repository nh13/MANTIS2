use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;
use env_logger::Env;

use crate::utils::built_info;

/// Detect microsatellite instability from a matched tumor-normal pair
#[derive(Parser, Debug)]
#[clap(name = "detect", verbatim_doc_comment, version = built_info::VERSION.as_str())]
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

    /// The path to the output
    #[clap(short = 'o', long, display_order = 5)]
    pub output: PathBuf,

    /// The minimum average per-base read quality for a read to pass the quality control filters
    #[clap(long, default_value = "25", display_order = 6)]
    pub min_read_mean_base_quality: u64,

    /// Minimum average per-base quality for the bases contained within the microsatellite locus.
    /// Reads that pass the read quality filter (above) will still fail quality control if the
    /// locus quality scores are too low.
    #[clap(long, default_value = "30", display_order = 7)]
    pub min_locus_mean_base_quality: u64,

    /// Minimum read length for a read to pass quality control. Only bases that are not clipped
    /// will be considered; in other words, soft-clipped or hard-clipped parts of the read do not
    /// count towards the length.
    #[clap(long, default_value = "35", display_order = 7)]
    pub min_read_length: u64,

    /// Minimum coverage (after QC filters) required for each of the normal and tumor samples for a
    ///  locus to be considered in the calculations.
    #[clap(long, default_value = "30", display_order = 7)]
    pub min_locus_coverage: u64,

    /// Minimum reads supporting a specific repeat count. Repeat counts that have less than this
    /// value will be discarded as part of outlier filtering
    #[clap(long, default_value = "3", display_order = 7)]
    pub min_repeat_reads: u64,

    /// Standard deviations from the mean before a repeat count is considered an outlier and
    /// discarded.
    #[clap(long, default_value = "3.0", display_order = 7)]
    pub outlier_standard_deviation: f64,

    /// The number of threads to use
    #[clap(long, default_value = "1", display_order = 3)]
    pub threads: u64,
}

// Run index
#[allow(clippy::too_many_lines)]
pub fn run(opts: &Opts) -> Result<(), anyhow::Error> {
    println!("{:?}", opts);
    Ok(())
}

/// Parse args and set up logging / tracing
pub fn setup() -> Opts {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    Opts::parse()
}
