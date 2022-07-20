#![allow(clippy::missing_errors_doc, clippy::missing_panics_doc)]
use std::process::exit;

use clap::{Parser, Subcommand};
use env_logger::Env;
use log::error;
use mantis_msi2_lib::tools::detect::{run as extract, Opts as DetectOps};
use mantis_msi2_lib::tools::repeat_counter::{run as count, Opts as RepeatCounterOpts};
use mantis_msi2_lib::tools::repeat_finder::{run as index, Opts as RepeatFinderOpts};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
#[clap(propagate_version = true)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

// This is where sub-commands are added.  The value of each enum should be the corresponding option
// struct
#[derive(Subcommand)]
enum Commands {
    /// Detect microsatellite instability from a matched tumor-normal pair
    Detect(DetectOps),
    /// Counts repeats in your BAM given a BED of repeats
    RepeatCounter(RepeatCounterOpts),
    /// Find repeats in your reference FASTA
    RepeatFinder(RepeatFinderOpts),
}

#[cfg(not(tarpaulin_include))]
#[allow(clippy::too_many_lines)]
fn main() {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();

    let result = match &cli.command {
        Commands::Detect(opts) => extract(opts),
        Commands::RepeatCounter(opts) => count(opts),
        Commands::RepeatFinder(opts) => index(opts),
    };

    if let Err(err) = result {
        error!("{:#}", err);
        exit(1);
    }
}
