
[![Check](https://github.com/nh13/mantis2/actions/workflows/build_and_test.yml/badge.svg?branch=develop-rs)](https://github.com/nh13/mantis2/actions/workflows/build_and_test.yml)
[![License](http://img.shields.io/badge/license-GPLv3-blue.svg)](https://github.com/nh13/mantis2/blob/develop-rs/LICENSE)
[![Language](http://img.shields.io/badge/language-rust-brightgreen.svg)](http://www.https://www.rust-lang.org/)

MANTIS 2
====

***Please do not use this for any serious work until a release has been made***


**This is a fork of https://github.com/OSU-SRLab/MANTIS as the latter is no longer being develop-rstained**

MANTIS2 (Microsatellite Analysis for Normal-Tumor InStability) is a program developed for detecting microsatellite instability from paired-end tumor-normal BAM files.
To perform analysis, the program needs a tumor BAM and a matched normal BAM file (produced using the same pipeline) to determine the instability score between the two samples within the pair.
Longer reads (ideally, 100 bp or longer) are recommended, as shorter reads are unlikely to entirely cover the microsatellite loci, and will be discarded after failing the quality control filters.

<!---toc start-->
* [MANTIS2](#mantis2)
* [Quickstart](#quickstart)
   * [Help](#help)
   * [To Build](#to-build)
   * [To Test](#to-test)
* [Microsatellite Loci BED File Format](#microsatellite-loci-bed-file-format)
* [Citation](#citation)

<!---toc end-->

# Quickstart

Find repeats in your reference FASTA (only needs to be run once per genome):

```bash
mantis-msi2 repeat-finder -i /path/to/genome.fasta -o /path/to/msi.loci.bed
```

Detect microsatellite instability (run for each tumor-normal matched pair):
```bash
mantis-msi2 detect -b /path/to/msi.loci.bed -g /path/to/genome.bed -n /path/to/normal.bam -t /path/to/tumor.bam -o /path/to/output/file.txt
```

## Help

List the available commands:

```bash
mantis-msi2 --help
```

List the help for a given command:

```bash
mantis-msi2 <command> --help
```

## To Build

```bash
cargo build --release
```

The executable is located in:

```bash
target/release/mantis-msi2
```

## To Test

```bash
cargo test
```

# Microsatellite Loci BED File Format

The tool requires input in a 6-column BED format, and can be produced via the `mantis-msi2 repeat-finder` tool.

In particular, the fourth (name) column of the BED file must contain the targeted repeating k-mer (e.g. AC) and the reference repeat count for it, e.g. “(AC)12”. 
A sample entryin your BED file might look like:

```
chr15 33256217 33256249 (AC)16 0 +
```

With the format being:


|   Column    |      1      |      2      |     3     |       4       |   5    |   6 |
| --- | --- | --- | --- | --- | --- | --- |
| Description |  Chromosome | Locus Start | Locus End | K-Mer Feature | Unused | Unused |


# Citation

For the original MANTIS software, please cite the [MANTIS publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5352334/).
Kautto, E. A., Bonneville, R., Miya, J., Yu, L., Krook, M. A., Reeser, J. W., & Roychowdhury, S. (2017). Performance evaluation for rapid detection of pan-cancer microsatellite instability with MANTIS. Oncotarget, 8(5), 7452–7463. http://doi.org/10.18632/oncotarget.13918. PMID: 27980218.

For tools, please cite this repository.
