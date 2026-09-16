//! Port of `bam_stat.py`: summarize mapping statistics for a BAM/SAM file.
//! Contract: see compatibility/commands.yaml entry `bam_stat.py`.

#[derive(Debug, Default, Clone, PartialEq)]
pub struct BamStatCounts {
    pub total_records: u64,
    pub qc_failed: u64,
    pub duplicates: u64,
}
