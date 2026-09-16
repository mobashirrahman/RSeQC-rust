//! One module per ported RSeQC command. Each module exposes a pure
//! computation entry point consumed by both `rseqc-cli` and `rseqc-python`
//! — no argument parsing or plot rendering lives here.

pub mod bam2fq;
pub mod bam_stat;
pub mod clipping_profile;
pub mod divide_bam;
pub mod read_duplication;
pub mod read_gc;
pub mod read_nvc;
pub mod read_quality;
pub mod split_paired_bam;
