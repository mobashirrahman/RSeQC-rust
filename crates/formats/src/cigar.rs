//! Shared CIGAR-derived coordinate helpers, ported from
//! `oracle/upstream-src/src/qcmodule/bam_cigar.py`. Each upstream function
//! has its own, sometimes inconsistent, per-op-kind branch table (e.g.
//! `fetch_exon`'s soft-clip handling differs from `fetch_intron`'s); port
//! each one explicitly rather than sharing a single generic walker, so
//! those upstream quirks are reproduced deliberately, not accidentally
//! unified away.

use noodles_sam::alignment::record::cigar::{Op, op::Kind};

/// Reference-genome `[start, end)` span covered by a CIGAR (0-based).
/// Ported from `bam_cigar.map_bounds`: M/D/N ops advance the span size,
/// I/S/H/P/=/X do not (matches upstream's `else: continue`, including the
/// same "no advance for '='/'X'" limitation noted on
/// [`fetch_exon_blocks`]). Equivalent to pysam's `reference_end`.
pub fn reference_span(start: usize, cigar: impl IntoIterator<Item = Op>) -> (usize, usize) {
    let mut span_size = 0usize;
    for op in cigar {
        if matches!(op.kind(), Kind::Match | Kind::Deletion | Kind::Skip) {
            span_size += op.len();
        }
    }
    (start, start + span_size)
}

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

/// Expands a CIGAR into a "long string" of per-read-base operation codes,
/// one byte per read base consumed. Ported from `bam_cigar.list2longstr`:
/// only M/I/S/=/X (op codes 0/1/4/7/8, all of which consume the read)
/// contribute bytes; D/N/H/P (codes 2/3/5/6, which don't consume the read)
/// contribute nothing and are skipped entirely -- not even a placeholder.
/// Sum of lengths of the M/I/S/=/X operations equals the returned length,
/// matching upstream's own docstring.
pub fn expand_cigar_to_read_ops(cigar: impl IntoIterator<Item = Op>) -> Vec<u8> {
    let mut out = Vec::new();

    for op in cigar {
        let byte = match op.kind() {
            Kind::Match => Some(b'M'),
            Kind::Insertion => Some(b'I'),
            Kind::SoftClip => Some(b'S'),
            Kind::SequenceMatch => Some(b'='),
            Kind::SequenceMismatch => Some(b'X'),
            Kind::Deletion | Kind::Skip | Kind::HardClip | Kind::Pad => None,
        };
        if let Some(byte) = byte {
            out.extend(std::iter::repeat_n(byte, op.len()));
        }
    }

    out
}

/// Deletion positions within a read, as `(read_position, deletion_size)`
/// pairs, read-position 0-based. Ported from `bam_cigar.fetch_deletion_range`.
///
/// Only M/S/I (op codes 0/4/1) advance the read-position counter (they
/// consume the read); a `D` op records `(current_position, its_length)`
/// WITHOUT advancing the counter further (a deletion doesn't consume read
/// bases); N/H/P and, matching the same upstream limitation as
/// [`fetch_exon_blocks`] and [`expand_cigar_to_read_ops`], `=`/`X` do
/// neither (upstream's `else: continue` -- codes 7/8 predate wide `=`/`X`
/// usage and are silently not treated as read-consuming here either).
pub fn fetch_deletion_range(cigar: impl IntoIterator<Item = Op>) -> Vec<(usize, usize)> {
    let mut read_pos = 0usize;
    let mut bounds = Vec::new();

    for op in cigar {
        match op.kind() {
            Kind::Match | Kind::SoftClip | Kind::Insertion => {
                read_pos += op.len();
            }
            Kind::Deletion => {
                bounds.push((read_pos, op.len()));
            }
            Kind::Skip | Kind::HardClip | Kind::Pad | Kind::SequenceMatch | Kind::SequenceMismatch => {}
        }
    }

    bounds
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

    #[test]
    fn expand_cigar_matches_docstring_example() {
        // [(0, 9), (4, 1)] ==> "MMMMMMMMMS", from bam_cigar.list2longstr's
        // own docstring.
        let cigar = vec![Op::new(Kind::Match, 9), Op::new(Kind::SoftClip, 1)];
        assert_eq!(expand_cigar_to_read_ops(cigar), b"MMMMMMMMMS".to_vec());
    }

    #[test]
    fn expand_cigar_skips_deletion_skip_hardclip_pad_entirely() {
        // 3M 2D 4N 2H 1P 2M: only the M ops (and I/S/=/X, not present here)
        // contribute bytes; D/N/H/P vanish completely, not even a placeholder.
        let cigar = vec![
            Op::new(Kind::Match, 3),
            Op::new(Kind::Deletion, 2),
            Op::new(Kind::Skip, 4),
            Op::new(Kind::HardClip, 2),
            Op::new(Kind::Pad, 1),
            Op::new(Kind::Match, 2),
        ];
        assert_eq!(expand_cigar_to_read_ops(cigar), b"MMMMM".to_vec());
    }

    #[test]
    fn expand_cigar_includes_insertion_and_sequence_match_mismatch() {
        let cigar = vec![
            Op::new(Kind::Insertion, 2),
            Op::new(Kind::SequenceMatch, 2),
            Op::new(Kind::SequenceMismatch, 1),
        ];
        assert_eq!(expand_cigar_to_read_ops(cigar), b"II==X".to_vec());
    }

    #[test]
    fn deletion_range_records_position_without_advancing() {
        // 5M2D3M: deletion recorded at read-position 5 (after the first
        // 5M), and the following 3M starts at position 5 too (the 2D does
        // not consume read bases).
        let cigar = vec![
            Op::new(Kind::Match, 5),
            Op::new(Kind::Deletion, 2),
            Op::new(Kind::Match, 3),
        ];
        assert_eq!(fetch_deletion_range(cigar), vec![(5, 2)]);
    }

    #[test]
    fn deletion_range_soft_clip_and_insertion_also_advance() {
        // 2S3M1I1D4M: read-position advances through S/M/I (2+3+1=6),
        // then the 1D is recorded at position 6.
        let cigar = vec![
            Op::new(Kind::SoftClip, 2),
            Op::new(Kind::Match, 3),
            Op::new(Kind::Insertion, 1),
            Op::new(Kind::Deletion, 1),
            Op::new(Kind::Match, 4),
        ];
        assert_eq!(fetch_deletion_range(cigar), vec![(6, 1)]);
    }

    #[test]
    fn deletion_range_multiple_deletions_and_skip_does_not_advance() {
        // 2M1D2M3N1D2M: skip (N) doesn't advance read-position, so the
        // second deletion is recorded at position 4 (2+2), not 7.
        let cigar = vec![
            Op::new(Kind::Match, 2),
            Op::new(Kind::Deletion, 1),
            Op::new(Kind::Match, 2),
            Op::new(Kind::Skip, 3),
            Op::new(Kind::Deletion, 1),
            Op::new(Kind::Match, 2),
        ];
        assert_eq!(fetch_deletion_range(cigar), vec![(2, 1), (4, 1)]);
    }

    #[test]
    fn reference_span_sums_match_deletion_skip_only() {
        // 10M5N8M2I3S: M/N/M advance (10+5+8=23), I/S do not.
        let cigar = vec![
            Op::new(Kind::Match, 10),
            Op::new(Kind::Skip, 5),
            Op::new(Kind::Match, 8),
            Op::new(Kind::Insertion, 2),
            Op::new(Kind::SoftClip, 3),
        ];
        assert_eq!(reference_span(100, cigar), (100, 123));
    }
}
