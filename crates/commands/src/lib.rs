//! One module per ported RSeQC command. Each module exposes a pure
//! computation entry point consumed by both `rseqc-cli` and `rseqc-python`
//! — no argument parsing or plot rendering lives here.

pub mod bam_stat;
pub mod split_paired_bam;
