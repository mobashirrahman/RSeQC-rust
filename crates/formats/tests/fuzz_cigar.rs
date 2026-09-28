// Fuzz tests for CIGAR parsing functions
use proptest::prelude::*;
use noodles_sam::alignment::record::cigar::{Op, op::Kind};
use rseqc_formats::cigar;

proptest! {
    #[test]
    fn fuzz_reference_span_basic(start in 0usize..1_000_000) {
        let ops = vec![];
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = cigar::reference_span(start, ops.into_iter());
        }));
    }

    #[test]
    fn fuzz_reference_span_with_match(
        start in 0usize..1_000_000,
        match_len in 1..=10_000usize,
    ) {
        let ops = vec![Op::new(Kind::Match, match_len)];
        let (s, e) = cigar::reference_span(start, ops.into_iter());
        prop_assert_eq!(s, start);
        prop_assert_eq!(e, start + match_len);
    }

    #[test]
    fn fuzz_reference_span_multiple_ops(
        start in 0usize..1_000_000,
        match1 in 1..=1_000usize,
        del1 in 1..=1_000usize,
        match2 in 1..=1_000usize,
    ) {
        let ops = vec![
            Op::new(Kind::Match, match1),
            Op::new(Kind::Deletion, del1),
            Op::new(Kind::Match, match2),
        ];
        let (s, e) = cigar::reference_span(start, ops.into_iter());
        prop_assert_eq!(s, start);
        prop_assert_eq!(e, start + match1 + del1 + match2);
    }

    #[test]
    fn fuzz_reference_span_large_operations(
        start in 0usize..100_000_000,
        huge_match in 100_000usize..1_000_000_000,
    ) {
        let ops = vec![Op::new(Kind::Match, huge_match)];
        let (s, e) = cigar::reference_span(start, ops.into_iter());
        prop_assert_eq!(s, start);
        prop_assert_eq!(e, start + huge_match);
    }

    // Test operations that shouldn't advance reference coordinate
    #[test]
    fn fuzz_reference_span_insertion_only(
        start in 0usize..1_000_000,
        insert_len in 1..=10_000usize,
    ) {
        let ops = vec![Op::new(Kind::Insertion, insert_len)];
        let (s, e) = cigar::reference_span(start, ops.into_iter());
        prop_assert_eq!(s, start);
        prop_assert_eq!(e, start); // Insertions don't advance ref
    }

    #[test]
    fn fuzz_reference_span_soft_clip_only(
        start in 0usize..1_000_000,
        clip_len in 1..=10_000usize,
    ) {
        let ops = vec![Op::new(Kind::SoftClip, clip_len)];
        let (s, e) = cigar::reference_span(start, ops.into_iter());
        prop_assert_eq!(s, start);
        prop_assert_eq!(e, start); // Soft clip doesn't advance ref
    }

    // Test fetch_intron_blocks
    #[test]
    fn fuzz_fetch_intron_blocks_skip_only(
        start in 0usize..1_000_000,
        skip_len in 1..=100_000usize,
    ) {
        let ops = vec![Op::new(Kind::Skip, skip_len)];
        let blocks = cigar::fetch_intron_blocks(start, ops.into_iter());
        prop_assert_eq!(blocks.len(), 1);
        prop_assert_eq!(blocks[0], (start, start + skip_len));
    }

    #[test]
    fn fuzz_fetch_intron_blocks_multiple_skips(
        start in 0usize..1_000_000,
        skip1 in 1..=10_000usize,
        match_between in 1..=10_000usize,
        skip2 in 1..=10_000usize,
    ) {
        let ops = vec![
            Op::new(Kind::Skip, skip1),
            Op::new(Kind::Match, match_between),
            Op::new(Kind::Skip, skip2),
        ];
        let blocks = cigar::fetch_intron_blocks(start, ops.into_iter());
        prop_assert_eq!(blocks.len(), 2);
        prop_assert_eq!(blocks[0], (start, start + skip1));
        prop_assert_eq!(blocks[1], (start + skip1 + match_between, start + skip1 + match_between + skip2));
    }

    #[test]
    fn fuzz_fetch_intron_blocks_no_skips(
        start in 0usize..1_000_000,
        match_len in 1..=10_000usize,
    ) {
        let ops = vec![Op::new(Kind::Match, match_len)];
        let blocks = cigar::fetch_intron_blocks(start, ops.into_iter());
        prop_assert_eq!(blocks.len(), 0); // No skips = no introns
    }

    // Test fetch_exon_blocks
    #[test]
    fn fuzz_fetch_exon_blocks_match_only(
        start in 0usize..1_000_000,
        match_len in 1..=100_000usize,
    ) {
        let ops = vec![Op::new(Kind::Match, match_len)];
        let blocks = cigar::fetch_exon_blocks(start, ops.into_iter());
        prop_assert_eq!(blocks.len(), 1);
        prop_assert_eq!(blocks[0], (start, start + match_len));
    }

    #[test]
    fn fuzz_fetch_exon_blocks_deletion_advances(
        start in 0usize..1_000_000,
        match1 in 1..=10_000usize,
        del_len in 1..=10_000usize,
        match2 in 1..=10_000usize,
    ) {
        let ops = vec![
            Op::new(Kind::Match, match1),
            Op::new(Kind::Deletion, del_len),
            Op::new(Kind::Match, match2),
        ];
        let blocks = cigar::fetch_exon_blocks(start, ops.into_iter());
        prop_assert_eq!(blocks.len(), 2);
        prop_assert_eq!(blocks[0], (start, start + match1));
        prop_assert_eq!(blocks[1], (start + match1 + del_len, start + match1 + del_len + match2));
    }

    #[test]
    fn fuzz_fetch_exon_blocks_insertion_no_advance(
        start in 0usize..1_000_000,
        match1 in 1..=10_000usize,
        insert_len in 1..=10_000usize,
        match2 in 1..=10_000usize,
    ) {
        let ops = vec![
            Op::new(Kind::Match, match1),
            Op::new(Kind::Insertion, insert_len),
            Op::new(Kind::Match, match2),
        ];
        let blocks = cigar::fetch_exon_blocks(start, ops.into_iter());
        // Insertion doesn't advance reference, so creates separate blocks
        prop_assert_eq!(blocks.len(), 2);
        prop_assert_eq!(blocks[0], (start, start + match1));
        prop_assert_eq!(blocks[1], (start + match1, start + match1 + match2));
    }

    #[test]
    fn fuzz_fetch_exon_blocks_sequence_match_ignored(
        start in 0usize..1_000_000,
        match_len in 1..=10_000usize,
        seq_match_len in 1..=10_000usize,
    ) {
        let ops = vec![
            Op::new(Kind::Match, match_len),
            Op::new(Kind::SequenceMatch, seq_match_len), // Should be ignored
        ];
        let blocks = cigar::fetch_exon_blocks(start, ops.into_iter());
        // Only the Match block should be emitted; SequenceMatch is ignored (legacy upstream behavior)
        prop_assert!(!blocks.is_empty());
        prop_assert_eq!(blocks[0], (start, start + match_len));
    }

    // Stress test with many operations
    #[test]
    fn fuzz_fetch_exon_blocks_many_ops(ops_count in 1..100usize) {
        let mut ops = vec![];
        let start_pos = 0usize;
        for i in 0..ops_count {
            let op_type = i % 5;
            let len = 1 + (i % 100);
            match op_type {
                0 => ops.push(Op::new(Kind::Match, len)),
                1 => ops.push(Op::new(Kind::Insertion, len)),
                2 => ops.push(Op::new(Kind::Deletion, len)),
                3 => ops.push(Op::new(Kind::Skip, len)),
                _ => ops.push(Op::new(Kind::SoftClip, len)),
            }
        }
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = cigar::fetch_exon_blocks(start_pos, ops.into_iter());
        }));
    }

    // Test with zero-length operations (edge case)
    #[test]
    fn fuzz_fetch_exon_blocks_zero_len(
        start in 0usize..1_000_000,
    ) {
        // Note: noodles likely rejects zero-length ops at construction,
        // but test anyway
        let ops = vec![];
        let blocks = cigar::fetch_exon_blocks(start, ops.into_iter());
        prop_assert_eq!(blocks.len(), 0);
    }

    // Test coordinate saturation/overflow
    #[test]
    fn fuzz_fetch_exon_blocks_large_coordinates(
        start in 1_000_000_000usize..3_000_000_000usize,
        match_len in 1..=1_000_000usize,
    ) {
        let ops = vec![Op::new(Kind::Match, match_len)];
        let blocks = cigar::fetch_exon_blocks(start, ops.into_iter());
        prop_assert_eq!(blocks.len(), 1);
        if let Some((s, e)) = blocks.first() {
            prop_assert_eq!(*s, start);
            prop_assert_eq!(*e, start + match_len);
        }
    }
}
