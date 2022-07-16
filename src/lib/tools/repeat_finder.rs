use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;
use env_logger::Env;

use crate::utils::built_info;

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
    pub min_bases: u64,

    /// The maximum number of bases for a repeat region must span to be called a microsatellite
    #[clap(short = 'M', default_value = "100", long, display_order = 4)]
    pub max_bases: u64,

    /// The minimum number of repeats for a repeat region must span to be called a microsatellite
    #[clap(short = 'r', default_value = "3", long, display_order = 5)]
    pub min_repeats: u64,

    /// The minimum k-mer length
    #[clap(short = 'l', default_value = "1", long, display_order = 6)]
    pub min_repeat_length: u64,

    /// The maximum k-mer length
    #[clap(short = 'L', default_value = "5", long, display_order = 7)]
    pub max_repeat_length: u64,
}

// Run extract
#[allow(clippy::too_many_lines)]
pub fn run(opts: &Opts) -> Result<(), anyhow::Error> {
    println!("{:?}", opts);
    // TODO
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
