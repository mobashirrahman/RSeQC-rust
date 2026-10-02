//! Independent truth cases for the shared scientific semantics.
//!
//! These are hand-derived expectations, not comparisons against upstream and not
//! recordings of what the port currently outputs. The audit's "shared scientific
//! semantics" layer asks for exactly this: pileup overlap/quality/depth
//! boundaries, CIGAR `M`/`=`/`X`/`N`/`D` rules, paired exon membership and a
//! strandedness truth table, verified against independent truth.
//!
//! The distinction matters because a differential suite cannot supply it. Every
//! existing CIGAR and interval test in this crate asserts either upstream's
//! behaviour (correct as a compatibility claim) or the port's own previous
//! output (which pins a regression without validating the science). A test that
//! passes only because the port agrees with itself is not evidence.
//!
//! Where the port deliberately reproduces an upstream quirk -- `fetch_exon_blocks`
//! ignoring `=`/`X`, `fetch_intron_blocks` not advancing on soft clips -- the
//! independent truth case is recorded ALONGSIDE the quirk, with the discrepancy
//! named. Compatibility and correctness are then separate, checkable claims
//! rather than one undifferentiated expectation.

use noodles_sam::alignment::record::cigar::{Op, op::Kind};
use rseqc_formats::cigar;
use rseqc_formats::interval::{Bed3, MergedRegions};

fn ops(spec: &[(Kind, usize)]) -> Vec<Op> {
    spec.iter().map(|(k, l)| Op::new(*k, *l)).collect()
}

// ---------------------------------------------------------------------------------
// CIGAR M/=/X/N/D rules
// ---------------------------------------------------------------------------------
// Independent truth: the SAM specification's reference-consuming op set is
// M, D, N, = and X. I, S, H and P do not consume reference. Every function below
// that walks coordinates must respect this rule, except where a divergence is
// explicitly recorded.

#[test]
fn reference_span_follows_the_sam_spec_exactly() {
    // The spec's rule, tested op by op: each op either advances the reference
    // cursor by its length or does not, with no other effect on the span.
    let cases: &[(&str, Vec<Op>, usize)] = &[
        ("10M", ops(&[(Kind::Match, 10)]), 10),
        ("5I", ops(&[(Kind::Insertion, 5)]), 0),
        ("7D", ops(&[(Kind::Deletion, 7)]), 7),
        ("9N", ops(&[(Kind::Skip, 9)]), 9),
        ("3S", ops(&[(Kind::SoftClip, 3)]), 0),
        ("4H", ops(&[(Kind::HardClip, 4)]), 0),
        ("2P", ops(&[(Kind::Pad, 2)]), 0),
        ("6=", ops(&[(Kind::SequenceMatch, 6)]), 6),
        ("8X", ops(&[(Kind::SequenceMismatch, 8)]), 8),
    ];
    for (label, cigar, expected_len) in cases {
        let (start, end) = cigar::reference_span(1000, cigar.clone());
        assert_eq!(
            end - start,
            *expected_len,
            "{label}: reference span must be {expected_len}, got {:?}",
            (start, end)
        );
    }
}

#[test]
fn reference_span_sums_all_reference_consuming_ops_in_a_mixed_cigar() {
    // 4M 3I 6D 2N 5S 7= 1X 2H 1P 3M
    // reference-consuming: 4 + 6 + 2 + 7 + 1 + 3 = 23
    let cigar = ops(&[
        (Kind::Match, 4),
        (Kind::Insertion, 3),
        (Kind::Deletion, 6),
        (Kind::Skip, 2),
        (Kind::SoftClip, 5),
        (Kind::SequenceMatch, 7),
        (Kind::SequenceMismatch, 1),
        (Kind::HardClip, 2),
        (Kind::Pad, 1),
        (Kind::Match, 3),
    ]);
    let (start, end) = cigar::reference_span(0, cigar);
    assert_eq!((start, end), (0, 23));
}

#[test]
fn exon_blocks_partition_the_reference_span_of_an_m_only_cigar() {
    // Independent invariant for `fetch_exon_blocks` on a plain M CIGAR: the exon
    // blocks must exactly tile [start, reference_end), with no gap and no
    // overlap. This is checkable arithmetic and does not depend on upstream.
    let cigar = ops(&[(Kind::Match, 7), (Kind::Deletion, 3), (Kind::Match, 5), (Kind::Skip, 11), (Kind::Match, 2)]);
    let blocks = cigar::fetch_exon_blocks(100, cigar.clone());
    let (start, end) = cigar::reference_span(100, cigar);

    assert_eq!(blocks, vec![(100, 107), (110, 115), (126, 128)]);
    assert_eq!(blocks[0].0, start, "first block must begin at the alignment start");
    assert_eq!(
        blocks.last().unwrap().1,
        end,
        "last block must end at the reference end"
    );
    for pair in blocks.windows(2) {
        assert!(
            pair[0].1 < pair[1].0,
            "consecutive exon blocks must not touch or overlap: {:?}",
            pair
        );
    }
    // D and N are not exons, so their lengths are exactly the gaps between blocks.
    let gap_total: usize = blocks.windows(2).map(|p| p[1].0 - p[0].1).sum();
    assert_eq!(gap_total, 3 + 11, "gaps must be exactly the D and N lengths");
}

#[test]
fn intron_blocks_are_exactly_the_n_operations() {
    let cigar = ops(&[
        (Kind::Match, 5),
        (Kind::Skip, 100),
        (Kind::Match, 5),
        (Kind::Skip, 20),
        (Kind::Match, 5),
    ]);
    assert_eq!(cigar::fetch_intron_blocks(1000, cigar), vec![(1005, 1105), (1110, 1130)]);
}

#[test]
fn deletion_and_skip_are_both_introns_but_only_n_is_recorded_as_such() {
    // The independent distinction: a deletion is a genomic gap with no spliced
    // evidence, so `fetch_intron_blocks` reports only `N`. A port that reported
    // `D` here would inflate the junction count.
    let cigar = ops(&[(Kind::Match, 10), (Kind::Deletion, 50), (Kind::Match, 10)]);
    assert_eq!(cigar::fetch_intron_blocks(0, cigar.clone()), Vec::new());
    assert_eq!(cigar::fetch_deletion_range(cigar), vec![(10, 50)]);
}

#[test]
fn read_operation_expansion_covers_exactly_the_read_consuming_ops() {
    // Independent invariant: the expanded string's length equals the sum of the
    // read-consuming op lengths, and each position carries the code of the op that
    // consumed it. Stated here because upstream's own docstring asserts it, so a
    // change that broke it would not be caught by a comparison alone.
    let cigar = ops(&[
        (Kind::SoftClip, 2),
        (Kind::Match, 3),
        (Kind::Insertion, 1),
        (Kind::Deletion, 4),
        (Kind::Match, 2),
        (Kind::Skip, 6),
        (Kind::SequenceMatch, 2),
        (Kind::SequenceMismatch, 1),
        (Kind::HardClip, 3),
    ]);
    // Read-consuming ops in order: 2S, 3M, 1I, 2M, 2=, 1X. The 4D and 6N sit
    // between them and contribute nothing.
    let expanded = cigar::expand_cigar_to_read_ops(cigar);
    // 2S + 3M + 1I + 2M + 2= + 1X = 11 read bases, in that order.
    assert_eq!(expanded, b"SSMMMIMM==X".to_vec());
    assert_eq!(expanded.len(), 2 + 3 + 1 + 2 + 2 + 1);
}

#[test]
fn deletion_positions_are_read_offsets_not_genomic_offsets() {
    // Independent invariant for `fetch_deletion_range`: every recorded position
    // must be expressible as a read offset, so it can never exceed the read
    // length consumed so far. A port that reported a genomic offset here would
    // produce plausible-looking but meaningless positions.
    let cigar = ops(&[(Kind::Match, 4), (Kind::Deletion, 3), (Kind::Match, 4)]);
    let read_ops = cigar::expand_cigar_to_read_ops(cigar.clone());
    let read_len = read_ops.len();
    for (pos, size) in cigar::fetch_deletion_range(cigar) {
        assert!(
            pos <= read_len,
            "deletion at read offset {pos} exceeds read length {read_len}"
        );
        assert!(size > 0, "a deletion must have a positive length");
    }
    assert_eq!(cigar::fetch_deletion_range(ops(&[
        (Kind::Match, 4),
        (Kind::Deletion, 3),
        (Kind::Match, 4)
    ])), vec![(4, 3)]);
}

// ---------------------------------------------------------------------------------
// Recorded divergences from spec-correct CIGAR semantics
// ---------------------------------------------------------------------------------
// The port reproduces two upstream inconsistencies. Each is asserted here as a
// DIVERGENCE so that it stays a deliberate compatibility decision: if one is ever
// fixed, this test fails and forces the compatibility record to be updated in the
// same change, rather than the behaviour drifting silently.

#[test]
fn divergence_exon_blocks_ignore_sequence_match_and_mismatch() {
    // Spec-correct behaviour would treat `=` and `X` as aligned blocks exactly
    // like `M`. Upstream's `fetch_exon` falls through to `continue`, so a CIGAR
    // built from `=`/`X` yields fewer (or zero) exon blocks.
    //
    // Truth: 5= 3X spans reference [200, 208) and is fully aligned.
    // Port:  yields nothing, because neither op advances the cursor.
    let cigar = ops(&[(Kind::SequenceMatch, 5), (Kind::SequenceMismatch, 3)]);
    assert_eq!(
        cigar::fetch_exon_blocks(200, cigar.clone()),
        Vec::new(),
        "upstream's `else: continue` drops =/X entirely; changing this is a \
         compatibility decision, not a bug fix"
    );
    // The same CIGAR under the spec-correct walker DOES span 8 bases, which is why
    // the divergence is a real defect in upstream rather than a formatting choice.
    assert_eq!(cigar::reference_span(200, cigar), (200, 208));
}

#[test]
fn divergence_soft_clip_advances_for_exons_but_not_for_introns() {
    // Spec-correct: a soft clip consumes read, not reference, so it must not move
    // the reference cursor in EITHER walker. Upstream's two walkers disagree.
    //
    // Port (upstream-compatible): fetch_exon_blocks advances on S, so the second
    // match block starts 3 bases later than the spec requires; fetch_intron_blocks
    // does not advance, so the intron starts where the spec says it should.
    let cigar = ops(&[
        (Kind::SoftClip, 3),
        (Kind::Match, 10),
        (Kind::Skip, 5),
        (Kind::Match, 5),
    ]);
    assert_eq!(
        cigar::fetch_exon_blocks(0, cigar.clone()),
        vec![(3, 13), (18, 23)],
        "upstream advances on S here"
    );
    assert_eq!(
        cigar::fetch_intron_blocks(0, cigar),
        vec![(10, 15)],
        "upstream does not advance on S here; the two walkers disagree"
    );
}

// ---------------------------------------------------------------------------------
// Paired exon membership and interval overlap boundaries
// ---------------------------------------------------------------------------------
// The boundary cases that decide whether an interval overlaps at all. Half-open
// [start, end) means touching intervals do NOT overlap, which is the single most
// common off-by-one in this class of code and the reason these are pinned.

#[test]
fn touching_intervals_do_not_overlap() {
    let exons: Vec<Bed3> = vec![("chr1".into(), 100, 200)];
    let regions = MergedRegions::new(&exons);
    // An interval ending exactly where the exon begins, and one starting exactly
    // where it ends, share a boundary but no base.
    assert_eq!(regions.overlap_length("chr1", 50, 100), 0);
    assert_eq!(regions.overlap_length("chr1", 200, 250), 0);
    // One base of overlap is one base counted.
    assert_eq!(regions.overlap_length("chr1", 99, 101), 1);
    assert_eq!(regions.overlap_length("chr1", 199, 201), 1);
    // Fully contained, and exactly equal.
    assert_eq!(regions.overlap_length("chr1", 120, 130), 10);
    assert_eq!(regions.overlap_length("chr1", 100, 200), 100);
    assert_eq!(regions.overlap_length("chr1", 50, 250), 100);
}

#[test]
fn empty_and_inverted_query_spans_count_nothing() {
    let regions = MergedRegions::new(&[("chr1".into(), 100, 200)]);
    assert_eq!(regions.overlap_length("chr1", 150, 150), 0, "empty span");
    assert_eq!(regions.overlap_length("chr1", 180, 120), 0, "inverted span");
    assert_eq!(regions.overlap_length("chrUn_gl000220", 100, 200), 0, "unknown contig");
}

#[test]
fn adjacent_exons_are_not_double_counted() {
    // Merged coverage must count each base once. Two exons sharing a boundary are
    // merged into one region, so a query spanning both reports their sum, not more.
    let regions = MergedRegions::new(&[
        ("chr1".into(), 100, 200),
        ("chr1".into(), 200, 300),
    ]);
    assert_eq!(regions.overlap_length("chr1", 100, 300), 200);
}

#[test]
fn overlapping_exons_are_not_double_counted() {
    let regions = MergedRegions::new(&[
        ("chr1".into(), 100, 200),
        ("chr1".into(), 150, 250),
    ]);
    // Union is [100, 250) = 150 bases, not 100 + 100.
    assert_eq!(regions.overlap_length("chr1", 100, 250), 150);
}

#[test]
fn exon_union_intersection_and_subtraction_on_hand_computed_regions() {
    // Independent arithmetic on the half-open union of each side.
    let a: Vec<Bed3> = vec![("chr1".into(), 100, 200), ("chr1".into(), 400, 500)];
    let b: Vec<Bed3> = vec![("chr1".into(), 150, 450)];

    // union(a) = [100,200) U [400,500)
    let union = rseqc_formats::interval::union_bed3(&a);
    assert_eq!(union, vec![("chr1".into(), 100, 200), ("chr1".into(), 400, 500)]);

    // intersection = [150,200) U [400,450)
    let inter = rseqc_formats::interval::intersect_bed3(&a, &b);
    assert_eq!(inter, vec![("chr1".into(), 150, 200), ("chr1".into(), 400, 450)]);

    // subtraction: [100,200) minus [150,200) = [100,150);
    //               [400,500) minus [400,450) = [450,500)
    let sub = rseqc_formats::interval::subtract_bed3(&a, &b);
    assert_eq!(sub, vec![("chr1".into(), 100, 150), ("chr1".into(), 450, 500)]);
}

#[test]
fn subtracting_a_contained_region_splits_the_outer_one() {
    // The boundary case a single-cut implementation gets wrong.
    let a: Vec<Bed3> = vec![("chr1".into(), 100, 200)];
    let b: Vec<Bed3> = vec![("chr1".into(), 120, 180)];
    assert_eq!(
        rseqc_formats::interval::subtract_bed3(&a, &b),
        vec![("chr1".into(), 100, 120), ("chr1".into(), 180, 200)]
    );
    // Subtracting a region that only touches the boundary removes nothing.
    let touching: Vec<Bed3> = vec![("chr1".into(), 200, 300)];
    assert_eq!(rseqc_formats::interval::subtract_bed3(&a, &touching), a);
}

// ---------------------------------------------------------------------------------
// Coordinate frame invariants
// ---------------------------------------------------------------------------------
// The rat refGene defect (a P0 finding) was a one-base frame error that survived
// because every consumer was self-consistent. These tests state the frame
// invariants independently, so a conversion cannot be self-consistently wrong.

#[test]
fn half_open_intervals_never_report_negative_or_inverted_lengths() {
    // Any span this code computes as end - start must be non-negative. A negative
    // or inverted span is the signature of a coordinate-frame error, which is
    // exactly what the refGene conversion defect produced.
    let cases: &[(usize, Vec<Op>)] = &[
        (0, ops(&[(Kind::Match, 10)])),
        (100, ops(&[(Kind::Match, 10), (Kind::Skip, 5), (Kind::Match, 10)])),
        (0, ops(&[(Kind::Deletion, 7), (Kind::Match, 3)])),
        (5, ops(&[(Kind::SequenceMatch, 4), (Kind::SequenceMismatch, 4)])),
    ];
    for (start, cigar) in cases {
        let (s, e) = cigar::reference_span(*start, cigar.clone());
        assert!(e >= s, "reference span inverted: ({s}, {e})");
        assert_eq!(s, *start, "reference span must start where asked");
        for (bs, be) in cigar::fetch_exon_blocks(*start, cigar.clone()) {
            assert!(be >= bs, "exon block inverted: ({bs}, {be})");
            assert!(be > bs, "zero-length exon block is not representable in BED12");
        }
        for (bs, be) in cigar::fetch_intron_blocks(*start, cigar.clone()) {
            assert!(be > bs, "intron block must be non-empty");
        }
    }
}

#[test]
fn exon_and_intron_blocks_stay_within_the_alignment_span() {
    let cigar = ops(&[
        (Kind::Match, 10),
        (Kind::Deletion, 5),
        (Kind::Match, 10),
        (Kind::Skip, 100),
        (Kind::Match, 10),
    ]);
    let (start, end) = cigar::reference_span(500, cigar.clone());
    for (bs, be) in cigar::fetch_exon_blocks(500, cigar.clone())
        .into_iter()
        .chain(cigar::fetch_intron_blocks(500, cigar))
    {
        assert!(bs >= start, "block {bs} starts before the alignment");
        assert!(be <= end, "block {be} ends after the alignment");
    }
}
