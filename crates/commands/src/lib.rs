//! One module per ported RSeQC command. Each module exposes a pure
//! computation entry point consumed by both `rseqc-cli` and `rseqc-python`
//! — no argument parsing or plot rendering lives here.

pub mod bam2fq;
pub mod bam_stat;
pub mod fpkm_count;
pub mod fpkm_uq;
pub mod clipping_profile;
pub mod deletion_profile;
pub mod mismatch_profile;
pub mod python_fmt;
pub mod infer_experiment;
pub mod inner_distance;
pub mod junction_annotation;
pub mod junction_saturation;
pub mod insertion_profile;
pub mod divide_bam;
pub mod read_duplication;
pub mod read_gc;
pub mod read_nvc;
pub mod read_distribution;
pub mod read_quality;
pub mod rna_fragment_size;
pub mod rpkm_saturation;
pub mod split_bam;
pub mod split_paired_bam;
pub mod tin;
