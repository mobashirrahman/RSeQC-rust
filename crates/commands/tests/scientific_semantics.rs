//! Independent truth cases for strandedness inference and paired exon membership.
//!
//! Two claims are separated here that are usually conflated:
//!
//!   * `infer_experiment`'s **arithmetic** -- which (read-end strand, gene strand)
//!     pairs belong to which fraction. That is pure combinatorics, checkable
//!     against a hand-built truth table with no reference implementation
//!     involved.
//!   * the **protocol mapping** from those fractions to a library-preparation
//!     conclusion. That is a biological interpretation, and the audit's P1
//!     scientific-interpretation finding is precisely that E1's strand
//!     expectation was taken from PCR selection rather than from the protocol.
//!     `spec1`/`spec2` do not name a protocol on their own, so this file asserts
//!     only what the numbers mean, and never asserts a protocol label.
//!
//! `split_bam`'s paired exon membership is checked the same way: the membership
//! rule is stated as a truth table over (read start in exon, mate start in exon,
//! mate mapped), which is checkable arithmetic, while the documented upstream
//! quirks (secondary/duplicate reads not filtered, an unmapped mate's start
//! ignored) are asserted as deliberate divergences.

use noodles_sam::{
    self as sam,
    alignment::{
        record::{Flags, MappingQuality, cigar::Op, cigar::op::Kind},
        record_buf::{Cigar, RecordBuf},
    },
    header::record::value::{Map, map::ReferenceSequence},
};
use rseqc_commands::infer_experiment::{compute_experiment, GeneRanges, Protocol};
use rseqc_commands::split_bam::{build_exon_ranges, classify_alignment, Category};
use std::collections::BTreeSet;
use std::num::NonZeroUsize;

fn header_of(names: &[&str]) -> sam::Header {
    let mut builder = sam::Header::builder();
    for name in names {
        builder = builder.add_reference_sequence(
            *name,
            Map::<ReferenceSequence>::new(NonZeroUsize::new(10_000).unwrap()),
        );
    }
    builder.build()
}

fn header() -> sam::Header {
    header_of(&["chr1"])
}

fn bam_records(header: &sam::Header, records: &[RecordBuf]) -> Vec<noodles_bam::Record> {
    use noodles_sam::alignment::io::Write as _;
    let mut buf = Vec::new();
    {
        let mut writer = noodles_bam::io::Writer::new(&mut buf);
        writer.write_header(header).unwrap();
        for r in records {
            writer.write_alignment_record(header, r).unwrap();
        }
    }
    let mut reader = noodles_bam::io::Reader::new(buf.as_slice());
    reader.read_header().unwrap();
    reader.records().map(|r| r.unwrap()).collect()
}

/// A read built through the builder, so every field is set explicitly. The builder
/// takes ownership, which is why the mutable "modify one flag" style used in
/// other test modules is not available here.
fn read(
    flags: Flags,
    ref_id: usize,
    start: usize,
    mapq: u8,
    mate: Option<(usize, usize)>,
) -> RecordBuf {
    let mut builder = RecordBuf::builder()
        .set_flags(flags)
        .set_reference_sequence_id(ref_id)
        .set_alignment_start(noodles_core::Position::new(start).unwrap())
        .set_mapping_quality(MappingQuality::new(mapq).unwrap())
        .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]));
    if let Some((mate_ref, mate_start)) = mate {
        builder = builder
            .set_mate_reference_sequence_id(mate_ref)
            .set_mate_alignment_start(noodles_core::Position::new(mate_start).unwrap());
    }
    builder.build()
}

/// One single-end read at `start` (1-based SAM position), forward or reverse.
fn se_read(start: usize, reverse: bool) -> RecordBuf {
    let flags = if reverse {
        Flags::REVERSE_COMPLEMENTED
    } else {
        Flags::empty()
    };
    read(flags, 0, start, 40, None)
}

/// One paired-end read: `first` selects read 1; the mate start is on chr1.
fn pe_read(start: usize, reverse: bool, first: bool, mate_start: usize) -> RecordBuf {
    let mut flags = Flags::SEGMENTED | Flags::MATE_REVERSE_COMPLEMENTED;
    flags |= if first {
        Flags::FIRST_SEGMENT
    } else {
        Flags::LAST_SEGMENT
    };
    if reverse {
        flags |= Flags::REVERSE_COMPLEMENTED;
    }
    read(flags, 0, start, 40, Some((0, mate_start)))
}

// ---------------------------------------------------------------------------------
// Strandedness: the (read strand, gene strand) truth table
// ---------------------------------------------------------------------------------
// The buckets, from arithmetic alone. A read overlapping a plus-strand gene and
// aligned forward is "++"; overlapping a minus-strand gene and aligned forward is
// "+-"; and so on. A read overlapping both strands joins them with ':', which
// matches neither bucket and so falls out of both fractions.

#[test]
fn single_end_strand_truth_table() {
    // Disjoint exons: [100,200) plus-strand, [5000,5100) minus-strand.
    let bed = "chr1\t100\t200\tplusGene\t0\t+\nchr1\t5000\t5100\tminusGene\t0\t-\n";
    let (ranges, _) = GeneRanges::parse(bed.as_bytes()).unwrap();
    let header = header();

    let cases: &[(&str, usize, bool, f64, f64)] = &[
        ("forward read over plus gene is ++", 150, false, 1.0, 0.0),
        ("reverse read over plus gene is -+", 150, true, 0.0, 1.0),
        ("forward read over minus gene is +-", 5050, false, 0.0, 1.0),
        ("reverse read over minus gene is --", 5050, true, 1.0, 0.0),
    ];

    for (label, start, reverse, spec1, spec2) in cases {
        let records = bam_records(&header, &[se_read(*start, *reverse)]);
        let r =
            compute_experiment(records.into_iter().map(Ok), &header, &ranges, 100_000, 30).unwrap();
        assert_eq!(r.protocol, Protocol::SingleEnd, "{label}");
        assert!(
            (r.spec1 - spec1).abs() < 1e-12,
            "{label}: spec1 (++,--) expected {spec1}, got {}",
            r.spec1
        );
        assert!(
            (r.spec2 - spec2).abs() < 1e-12,
            "{label}: spec2 (+-,-+) expected {spec2}, got {}",
            r.spec2
        );
    }
}

#[test]
fn single_end_fractions_are_complementary_and_sum_to_one() {
    // Independent invariant: spec1, spec2 and undetermined are three disjoint
    // buckets over the sampled reads, so they sum to exactly 1.
    let bed = "chr1\t100\t200\tplusGene\t0\t+\nchr1\t5000\t5100\tminusGene\t0\t-\n";
    let (ranges, _) = GeneRanges::parse(bed.as_bytes()).unwrap();
    let header = header();

    for reverse in [false, true] {
        let records = bam_records(
            &header,
            &[
                se_read(150, reverse),
                se_read(5050, !reverse),
                se_read(160, !reverse),
                se_read(5060, reverse),
            ],
        );
        let r =
            compute_experiment(records.into_iter().map(Ok), &header, &ranges, 100_000, 30).unwrap();
        let total = r.spec1 + r.spec2 + r.undetermined;
        assert!(
            (total - 1.0).abs() < 1e-12,
            "fractions must sum to 1, got {total} (spec1={} spec2={} undet={})",
            r.spec1,
            r.spec2,
            r.undetermined
        );
        assert!((0.0..=1.0).contains(&r.spec1));
        assert!((0.0..=1.0).contains(&r.spec2));
    }
}

#[test]
fn reads_overlapping_both_strands_are_undetermined_not_assigned() {
    // Independent rule: a read overlapping a plus-strand AND a minus-strand gene
    // cannot be assigned to either fraction. Getting this wrong would silently
    // inflate whichever fraction the overlapping reads favoured.
    let bed = "chr1\t100\t300\tplusGene\t0\t+\nchr1\t200\t400\tminusGene\t0\t-\n";
    let (ranges, _) = GeneRanges::parse(bed.as_bytes()).unwrap();
    let header = header();

    let records = bam_records(&header, &[se_read(250, false)]);
    let r =
        compute_experiment(records.into_iter().map(Ok), &header, &ranges, 100_000, 30).unwrap();
    assert_eq!(r.sampled_count, 1, "the read is still sampled");
    assert_eq!(r.spec1, 0.0, "an ambiguous read is not counted in ++/--");
    assert_eq!(r.spec2, 0.0, "an ambiguous read is not counted in +-/-+");
    assert!(
        (r.undetermined - 1.0).abs() < 1e-12,
        "an ambiguous read must be undetermined, got {}",
        r.undetermined
    );
}

#[test]
fn paired_end_buckets_group_by_read_number() {
    // Independent truth for the paired-end grouping, read off the bucket names the
    // report itself prints: "1++,1--,2+-,2-+" and "1+-,1-+,2++,2--".
    let bed = "chr1\t100\t200\tplusGene\t0\t+\nchr1\t5000\t5100\tminusGene\t0\t-\n";
    let (ranges, _) = GeneRanges::parse(bed.as_bytes()).unwrap();
    let header = header();

    let cases: &[(&str, usize, bool, bool, f64)] = &[
        // label, start, reverse, first_segment, expected spec1
        ("read1 forward over plus gene -> 1++ -> spec1", 150, false, true, 1.0),
        ("read1 forward over minus gene -> 1+- -> spec2", 5050, false, true, 0.0),
        ("read2 forward over plus gene -> 2++ -> spec2", 150, false, false, 0.0),
        ("read2 forward over minus gene -> 2+- -> spec1", 5050, false, false, 1.0),
    ];

    for (label, start, reverse, first, spec1) in cases {
        let records = bam_records(&header, &[pe_read(*start, *reverse, *first, *start)]);
        let r =
            compute_experiment(records.into_iter().map(Ok), &header, &ranges, 100_000, 30).unwrap();
        assert_eq!(r.protocol, Protocol::PairEnd, "{label}");
        assert!(
            (r.spec1 - spec1).abs() < 1e-12,
            "{label}: expected spec1 {spec1}, got {}",
            r.spec1
        );
    }
}

#[test]
fn paired_and_single_evidence_together_report_a_mixture() {
    // Independent rule: a sample with both paired and single-end usable reads is a
    // mixture with no fractions attributed, because neither model fits. A port
    // that picked one would make a mixed library look like a clean stranded or
    // unstranded preparation.
    let bed = "chr1\t100\t200\tplusGene\t0\t+\n";
    let (ranges, _) = GeneRanges::parse(bed.as_bytes()).unwrap();
    let header = header();

    let records = bam_records(
        &header,
        &[se_read(150, false), pe_read(150, false, true, 160)],
    );
    let r =
        compute_experiment(records.into_iter().map(Ok), &header, &ranges, 100_000, 30).unwrap();
    assert_eq!(r.protocol, Protocol::Mixture);
    assert_eq!(r.spec1, 0.0);
    assert_eq!(r.spec2, 0.0);
    assert_eq!(r.undetermined, 0.0);
}

#[test]
fn low_quality_duplicate_and_unmapped_reads_are_not_sampled() {
    // Independent rule about what counts as evidence: a read that would be
    // discarded must not contribute to any fraction. This boundary decides
    // whether a fraction is computed from real reads at all.
    let bed = "chr1\t100\t200\tplusGene\t0\t+\n";
    let (ranges, _) = GeneRanges::parse(bed.as_bytes()).unwrap();
    let header = header();

    let records = bam_records(
        &header,
        &[
            read(Flags::empty(), 0, 150, 5, None),                 // below the cutoff
            read(Flags::DUPLICATE, 0, 150, 40, None),             // duplicate
            read(Flags::SECONDARY, 0, 150, 40, None),             // secondary
            read(Flags::QC_FAIL, 0, 150, 40, None),               // QC fail
        ],
    );
    let r =
        compute_experiment(records.into_iter().map(Ok), &header, &ranges, 100_000, 30).unwrap();
    assert_eq!(r.sampled_count, 0, "none of these reads is usable evidence");
    assert_eq!(r.protocol, Protocol::Mixture);
    assert_eq!(r.spec1, 0.0);
    assert_eq!(r.spec2, 0.0);
}

#[test]
fn no_fraction_is_reported_without_strandedness_evidence() {
    // The audit's P1 finding, pinned structurally: with no usable read overlapping
    // an annotated strand, nothing is attributed to any fraction. `spec1`/`spec2`
    // are the measured quantities; mapping them to a library-preparation protocol
    // is the caller's decision and is deliberately not encoded here, because
    // deriving the expected protocol from PCR selection metadata is what the
    // audit found unsupported.
    let (ranges, _) = GeneRanges::parse("chr1\t0\t100\tg\t0\t+\n".as_bytes()).unwrap();
    let header = header();
    let records = bam_records(&header, &[se_read(5000, false)]);
    let r =
        compute_experiment(records.into_iter().map(Ok), &header, &ranges, 100_000, 30).unwrap();
    assert_eq!(r.sampled_count, 0);
    assert_eq!(r.protocol, Protocol::Mixture);
    assert_eq!(r.spec1, 0.0);
    assert_eq!(r.spec2, 0.0);
}

#[test]
fn the_report_states_buckets_not_a_protocol_label_for_a_stratum() {
    // The rendered text names which (read, gene) strand combinations each
    // fraction covers. That is the auditable claim; it does not name a library
    // preparation, which is what makes it safe to report without protocol
    // metadata.
    use rseqc_commands::infer_experiment::{render_results, ExperimentResult};
    let r = ExperimentResult {
        protocol: Protocol::SingleEnd,
        spec1: 0.9876,
        spec2: 0.0100,
        undetermined: 0.0024,
        sampled_count: 0,
        stopped_at_eof: false,
    };
    let out = render_results(&r);
    assert!(out.contains(r#""++,--""#), "must name its first bucket");
    assert!(out.contains(r#""+-,-+""#), "must name its second bucket");
    assert!(!out.to_lowercase().contains("unstranded"));
    assert!(!out.to_lowercase().contains("stranded protocol"));
}

// ---------------------------------------------------------------------------------
// Paired exon membership: the split_bam truth table
// ---------------------------------------------------------------------------------
// "In" iff the read's start or (mapped) mate's start lands on a base covered by an
// exon; "Ex" otherwise; "Junk" for QC-fail or unmapped. Exons are half-open
// [start, end), so a start exactly at `end` is outside.

fn exon_ranges() -> rseqc_formats::interval::MergedRegions {
    build_exon_ranges(&[("chr1".into(), 100, 200), ("chr1".into(), 1000, 1100)])
}

fn classify(start: usize, mate_start: Option<usize>) -> Category {
    let header = header();
    let ranges = exon_ranges();
    let record = match mate_start {
        Some(ms) => read(
            Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::MATE_REVERSE_COMPLEMENTED,
            0,
            start,
            60,
            Some((0, ms)),
        ),
        None => read(Flags::MATE_UNMAPPED, 0, start, 60, None),
    };
    let records = bam_records(&header, &[record]);
    classify_alignment(&records[0], &header, &ranges).unwrap()
}

#[test]
fn paired_exon_membership_truth_table() {
    // Positions are 1-based SAM coordinates, so the 0-based coordinate the
    // classifier works in is `position - 1`. Exons are half-open [100, 200), so
    // 0-based 100..=199 are inside and 0-based 200 is the first base outside --
    // which is SAM position 201. Getting this conversion wrong in either
    // direction is precisely the class of defect the rat refGene finding was, so
    // the boundary cases are stated on both frames explicitly.
    let cases: &[(&str, usize, usize, Category)] = &[
        ("read start inside exon -> In", 150, 5000, Category::In),
        ("mate start inside exon -> In", 5000, 150, Category::In),
        ("both inside -> In", 150, 1050, Category::In),
        ("neither inside -> Ex", 5000, 6000, Category::Ex),
        // 0-based 100, the exon's first base.
        ("read start at exon start -> In", 101, 6000, Category::In),
        // 0-based 99, one base before the exon.
        ("read start one base before exon -> Ex", 100, 6000, Category::Ex),
        // 0-based 199, the exon's last base.
        ("read start at exon last base -> In", 200, 6000, Category::In),
        // 0-based 200, one past the exon's half-open end.
        ("read start one base past exon end -> Ex", 201, 6000, Category::Ex),
        // Second exon [1000, 1100): 0-based 1000 and 1099 are in, 1100 is not.
        ("second exon, at its first base -> In", 1001, 6000, Category::In),
        ("second exon, at its last base -> In", 1100, 6000, Category::In),
        ("second exon, one past its end -> Ex", 1101, 6000, Category::Ex),
    ];

    for (label, start, mate, expected) in cases {
        assert_eq!(classify(*start, Some(*mate)), *expected, "{label}");
    }
}

#[test]
fn unmapped_mate_start_is_not_tested() {
    // Documented upstream behaviour: when the mate is unmapped its start field is
    // meaningless (typically 0), so only the read's own start is tested. A port
    // consulting the zero position would classify reads by an artefact.
    let header = header();
    let ranges = exon_ranges();
    // Mate start field deliberately set to 1 (which would be inside nothing, but
    // set to a value inside an exon instead to make the point). Position 150 is
    // inside the first exon.
    let record = read(
        Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::MATE_UNMAPPED,
        0,
        5000,
        60,
        None,
    );
    let records = bam_records(&header, &[record]);
    // Read start 5000 is outside every exon, and the mate is ignored.
    assert_eq!(
        classify_alignment(&records[0], &header, &ranges).unwrap(),
        Category::Ex
    );
}

#[test]
fn qc_fail_and_unmapped_reads_are_junk() {
    let header = header();
    let ranges = exon_ranges();

    // Even though its start IS inside an exon, a QC-fail read is junk, not "In".
    let records = bam_records(&header, &[read(Flags::QC_FAIL, 0, 150, 60, None)]);
    assert_eq!(
        classify_alignment(&records[0], &header, &ranges).unwrap(),
        Category::Junk
    );

    // An unmapped read, which carries no reference position at all.
    let records = bam_records(
        &header,
        &[RecordBuf::builder()
            .set_flags(Flags::UNMAPPED)
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .build()],
    );
    assert_eq!(
        classify_alignment(&records[0], &header, &ranges).unwrap(),
        Category::Junk
    );
}

#[test]
fn divergence_secondary_and_duplicates_are_not_filtered() {
    // Documented upstream divergence: split_bam routes ONLY QC-fail and unmapped
    // reads to junk. Secondary and duplicate alignments are classified normally,
    // unlike most other commands in this port. Asserted so a future "fix" has to
    // update the compatibility record rather than change behaviour silently.
    let header = header();
    let ranges = exon_ranges();

    for (label, flags) in [
        ("secondary", Flags::SECONDARY),
        ("duplicate", Flags::DUPLICATE),
    ] {
        let records = bam_records(&header, &[read(flags, 0, 150, 60, None)]);
        assert_eq!(
            classify_alignment(&records[0], &header, &ranges).unwrap(),
            Category::In,
            "{label} reads are classified, not discarded"
        );
    }
}

#[test]
fn exon_ranges_are_uppercased_on_both_sides() {
    // Upstream uppercases chromosome names when building the exon tree and the
    // classifier uppercases the query too, so a lowercase annotation still
    // matches. Pinned because a one-sided change would silently classify
    // everything as "Ex" -- plausible output, wrong in a way no parity test would
    // flag if both arms shared the bug.
    let ranges = build_exon_ranges(&[("chr1".into(), 100, 200)]);
    let header = header();
    let records = bam_records(&header, &[se_read(150, false)]);
    assert_eq!(
        classify_alignment(&records[0], &header, &ranges).unwrap(),
        Category::In
    );
}

#[test]
fn a_read_on_an_unannotated_contig_is_ex_not_junk() {
    // Independent rule: "not in an exon" and "not usable" are different verdicts.
    // A mapped read on a contig with no annotation is a legitimate non-exonic
    // read, so calling it junk would corrupt the in/ex/junk partition.
    let header = header_of(&["chr1", "chrUn_ungl"]);
    let ref_id = (0..header.reference_sequences().len())
        .find(|i| {
            header
                .reference_sequences()
                .get_index(*i)
                .map(|(n, _)| n == "chrUn_ungl")
                .unwrap_or(false)
        })
        .expect("chrUn_ungl must be in the header");

    let ranges = exon_ranges();
    let records = bam_records(&header, &[read(Flags::empty(), ref_id, 10, 60, None)]);
    assert_eq!(
        classify_alignment(&records[0], &header, &ranges).unwrap(),
        Category::Ex
    );
}

#[test]
fn gene_ranges_strand_set_is_ordered_deterministically() {
    // Independent determinism requirement: `find_strands` returns a BTreeSet, so
    // the joined key is stable regardless of insertion order. Without this, the
    // reported fraction would depend on iteration order and could differ between
    // runs on identical input.
    let bed = "chr1\t100\t300\tplusGene\t0\t+\nchr1\t200\t400\tminusGene\t0\t-\n";
    let (ranges, _) = GeneRanges::parse(bed.as_bytes()).unwrap();
    let first = ranges.find_strands("chr1", 250, 260);
    assert_eq!(first, BTreeSet::from(["+".to_string(), "-".to_string()]));
    for _ in 0..8 {
        assert_eq!(ranges.find_strands("chr1", 250, 260), first);
    }
}
