//! Shared CIGAR-derived coordinate helpers, ported from
//! `oracle/upstream-src/src/qcmodule/bam_cigar.py`. Each upstream function
//! has its own, sometimes inconsistent, per-op-kind branch table (e.g.
//! `fetch_exon`'s soft-clip handling differs from `fetch_intron`'s); port
//! each one explicitly rather than sharing a single generic walker, so
//! those upstream quirks are reproduced deliberately, not accidentally
//! unified away.

use noodles_sam::alignment::record::cigar::{Op, op::Kind};

/// Exon (aligned-match) blocks for a CIGAR string, as 0-based
/// `[start, end)` reference coordinates. Ported from `bam_cigar.fetch_exon`.
///
/// Upstream only advances/emits on BAM op code 0 (`M`); codes 7/8 (`=`/`X`,
/// [`Kind::SequenceMatch`]/[`Kind::SequenceMismatch`]) fall into upstream's
/// `else: continue` and are silently ignored (no coordinate advance, no
/// block emitted) even though they also represent aligned bases per the SAM
/// spec. This is very likely an upstream limitation predating widespread
/// `=`/`X` usage, but it is preserved here rather than "fixed" per
/// PORTING_PLAN.md's compatibility-mode policy — a CIGAR built entirely
/// from `=`/`X` would silently yield zero exon blocks from this function,
/// same as upstream.
pub fn fetch_exon_blocks(start: usize, cigar: impl IntoIterator<Item = Op>) -> Vec<(usize, usize)> {
    let mut chrom_st = start;
    let mut blocks = Vec::new();

    for op in cigar {
        match op.kind() {
            Kind::Match => {
                blocks.push((chrom_st, chrom_st + op.len()));
                chrom_st += op.len();
            }
            Kind::Insertion => {
                // Insertion relative to the reference: consumes read, not
                // reference, so no coordinate advance.
            }
            Kind::Deletion | Kind::Skip | Kind::SoftClip => {
                chrom_st += op.len();
            }
            Kind::HardClip | Kind::Pad | Kind::SequenceMatch | Kind::SequenceMismatch => {
                // Matches upstream's `else: continue` -- no advance, no block.
            }
        }
    }

    blocks
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_and_gaps_advance_and_split_blocks() {
        // 10M 5N 8M 3S: two match blocks separated by a skip, trailing
        // soft-clip advances the position but doesn't extend/start a block.
        let cigar = vec![
            Op::new(Kind::Match, 10),
            Op::new(Kind::Skip, 5),
            Op::new(Kind::Match, 8),
            Op::new(Kind::SoftClip, 3),
        ];
        let blocks = fetch_exon_blocks(100, cigar);
        assert_eq!(blocks, vec![(100, 110), (115, 123)]);
    }

    #[test]
    fn insertion_does_not_advance_reference_position() {
        // 4M 2I 4M: the insertion doesn't consume reference, so the second
        // match block starts immediately after the first, not offset by 2.
        let cigar = vec![
            Op::new(Kind::Match, 4),
            Op::new(Kind::Insertion, 2),
            Op::new(Kind::Match, 4),
        ];
        let blocks = fetch_exon_blocks(0, cigar);
        assert_eq!(blocks, vec![(0, 4), (4, 8)]);
    }

    #[test]
    fn sequence_match_mismatch_ops_are_silently_ignored() {
        // Preserves the upstream quirk: '=' and 'X' ops produce no blocks
        // and do not advance the coordinate, unlike plain 'M'.
        let cigar = vec![
            Op::new(Kind::SequenceMatch, 5),
            Op::new(Kind::Match, 3),
            Op::new(Kind::SequenceMismatch, 2),
        ];
        let blocks = fetch_exon_blocks(50, cigar);
        assert_eq!(blocks, vec![(50, 53)]);
    }
}
