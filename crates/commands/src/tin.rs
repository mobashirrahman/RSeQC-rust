//! Port of `tin.py`: calculate transcript integrity number (TIN) per
//! transcript, from `oracle/upstream-src/scripts/tin.py`. TIN measures
//! 5'-3' read coverage uniformity via a Shannon-entropy formula over
//! coverage sampled at up to `--sample-size` roughly-equally-spaced
//! transcript positions.
//!
//! No BAI index support (see crates/formats module docs): the
//! `samfile.fetch(chrom, start, end)`/`samfile.pileup(...)` region
//! queries upstream relies on are served here from an in-memory,
//! per-chromosome, start-sorted read index built once per BAM (same
//! `IndexedReads` pattern as `RNA_fragment_size.py`), not real htslib
//! random access.
//!
//! **Pileup defaults are mirrored from `pysam.AlignmentFile.pileup()`**:
//! `ignore_overlaps=True` (pysam's own default) deduplicates overlapping
//! paired-end mates at a shared reference position, counting only the
//! higher-quality base once. The implementation also applies the default
//! depth cap (8000) while retaining duplicates in the fetch-based helpers
//! below. `flag_filter` excludes duplicate-flagged reads (BAM_FDUP)
//! from coverage -- note this differs from `check_min_reads`/
//! `estimate_bg_noise`, which use `fetch()` (no `flag_filter`) and so DO
//! count duplicates. `min_base_quality` defaults to 0 (all bases included),
//! matching pysam's default. `query_length` (used by `estimate_bg_noise`)
//! is taken directly from the read's `SEQ` length, matching pysam's common
//! case; pysam's CIGAR-based inference fallback for a missing `SEQ` (`*`)
//! is not replicated.

use std::collections::{HashMap, HashSet};
use std::io::{self, BufRead};

use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::record::cigar::op::Kind};
use rseqc_formats::interval::{Bed3, MergedRegions};

use crate::python_fmt::{numpy_mean, numpy_std, python_str_float};

/// One BAM record's relevant fields for region-fetch/pileup queries,
/// pre-filtered at build time on the three flags every upstream
/// fetch()-consuming function checks (`is_qcfail`/`is_unmapped`/
/// `is_secondary`); `is_duplicate` is NOT filtered at build time because
/// `check_min_reads`/`estimate_bg_noise` (via `fetch()`) and
/// `genebody_coverage` (via `pileup()`'s default `flag_filter`) disagree
/// on whether duplicates count -- see module docs.
#[derive(Debug, Clone, Default)]
pub struct IndexedRead {
    /// Query/template name used to identify paired mates for pileup's
    /// `ignore_overlaps=True` behavior.
    pub query_name: String,
    pub is_paired: bool,
    pub start: i64,
    pub end: i64,
    pub is_duplicate: bool,
    /// `(ref_start, ref_end, query_start)` triples for M/=/X blocks;
    /// `query_start` is the index into `qualities` of the base aligned
    /// to `ref_start`.
    pub match_blocks: Vec<(i64, i64, usize)>,
    /// `(ref_start, ref_end)` pairs for D (deletion) and N (skip) CIGAR ops.
    /// These spans appear in pileup but have no passing bases, so should
    /// mark positions as visited (int 0) not unvisited (float 0.0).
    pub skip_delete_blocks: Vec<(i64, i64)>,
    pub qualities: Vec<u8>,
    pub query_length: i64,
    /// BAM_FPROPER_PAIR. pysam's default "samtools" pileup stepper drops
    /// orphans (paired but not properly paired) and htslib only adjusts
    /// overlapping mates that are properly paired.
    pub is_proper_pair: bool,
    /// BAM_FMUNMAP.
    pub mate_unmapped: bool,
    /// Mate mapped to a different reference sequence than this read.
    pub mate_on_other_reference: bool,
    /// 0-based mate position (`-1` when absent) and TLEN, as htslib's
    /// `overlap_push` inspects them.
    pub mate_start: i64,
    pub template_length: i64,
    /// Decoded read bases (for htslib's overlap base-equality test).
    pub sequence: Vec<u8>,
    /// CIGAR as `(BAM op code, length)` pairs.
    pub cigar: Vec<(u8, i64)>,
}

/// Builds a per-chromosome, start-sorted read index from a BAM. Applies
/// the qcfail/unmapped/secondary filter common to every upstream
/// fetch()-based function in this module.
pub fn build_read_index<I>(records: I, header: &sam::Header) -> io::Result<HashMap<String, Vec<IndexedRead>>>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut by_chrom: HashMap<String, Vec<IndexedRead>> = HashMap::new();

    for result in records {
        let record = result?;
        let Some(ref_id) = record.reference_sequence_id().transpose()? else { continue };
        let Some((chrom_bstr, _)) = header.reference_sequences().get_index(ref_id) else { continue };
        let Some(read) = to_indexed_read(&record)? else { continue };
        by_chrom.entry(chrom_bstr.to_string()).or_default().push(read);
    }

    for reads in by_chrom.values_mut() {
        reads.sort_by_key(|r| r.start);
        reads.shrink_to_fit();
    }

    Ok(by_chrom)
}

/// Decodes one BAM record into the per-read index representation, or `None`
/// when the record is dropped by the qcfail/unmapped/secondary filter that
/// every upstream fetch()-consuming function in this module applies.
///
/// Shared by the whole-file index (`build_read_index`) and the windowed
/// driver (`compute_tin_windowed`) so both decode records identically.
fn to_indexed_read(record: &bam::Record) -> io::Result<Option<IndexedRead>> {
    let flags = record.flags();
    if flags.is_qc_fail() || flags.is_unmapped() || flags.is_secondary() {
        return Ok(None);
    }
    let Some(ref_id) = record.reference_sequence_id().transpose()? else { return Ok(None) };
    let Some(pos) = record.alignment_start().transpose()? else { return Ok(None) };

    let start = (pos.get() - 1) as i64;
    let ops: Vec<_> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;
    let qualities: Vec<u8> = record.quality_scores().iter().collect();

    let mut ref_pos = start;
    let mut query_pos = 0usize;
    let mut match_blocks = Vec::new();
    let mut skip_delete_blocks = Vec::new();
    for op in &ops {
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                match_blocks.push((ref_pos, ref_pos + op.len() as i64, query_pos));
                ref_pos += op.len() as i64;
                query_pos += op.len();
            }
            Kind::Insertion | Kind::SoftClip => {
                query_pos += op.len();
            }
            Kind::Deletion | Kind::Skip => {
                skip_delete_blocks.push((ref_pos, ref_pos + op.len() as i64));
                ref_pos += op.len() as i64;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    let mate_ref_id = record.mate_reference_sequence_id().transpose()?;
    let mate_start = match record.mate_alignment_start().transpose()? {
        Some(p) => (p.get() - 1) as i64,
        None => -1,
    };
    Ok(Some(IndexedRead {
        query_name: record.name().map(|n| n.to_string()).unwrap_or_default(),
        is_paired: flags.is_segmented(),
        start,
        end: ref_pos,
        is_duplicate: flags.is_duplicate(),
        match_blocks,
        skip_delete_blocks,
        query_length: record.sequence().len() as i64,
        qualities,
        is_proper_pair: flags.is_properly_segmented(),
        mate_unmapped: flags.is_mate_unmapped(),
        mate_on_other_reference: mate_ref_id.is_some_and(|m| m != ref_id),
        mate_start,
        template_length: i64::from(record.template_length()),
        sequence: record.sequence().iter().collect(),
        cigar: ops.iter().map(|op| (cigar_code(op.kind()), op.len() as i64)).collect(),
    }))
}

fn reads_starting_in(reads: &[IndexedRead], tx_start: i64, tx_end: i64) -> &[IndexedRead] {
    let lo = reads.partition_point(|r| r.start < tx_start);
    let hi = reads.partition_point(|r| r.start < tx_end);
    &reads[lo..hi]
}

/// Ports the historical distinct-read-start coverage criterion: strictly
/// more than `cutoff` distinct read-start positions within `[tx_start,
/// tx_end)` (upstream's `len(read_starts) > cutoff`, not `>=`).
pub fn check_min_reads(reads: &[IndexedRead], tx_start: i64, tx_end: i64, cutoff: i64) -> bool {
    let mut starts = HashSet::new();
    for r in reads_starting_in(reads, tx_start, tx_end) {
        starts.insert(r.start);
        if starts.len() as i64 > cutoff {
            return true;
        }
    }
    false
}

/// Estimates background/intronic signal from reads starting within
/// `[tx_start, tx_end)` whose `[start, start+query_length)` span does
/// NOT overlap any (already-unioned) exon interval on `chrom`.
pub fn estimate_bg_noise(reads: &[IndexedRead], tx_start: i64, tx_end: i64, exon_ranges: &MergedRegions, chrom: &str) -> f64 {
    let mut intron_signal = 0.0;
    for r in reads_starting_in(reads, tx_start, tx_end) {
        if exon_ranges.overlap_length(chrom, r.start, r.start + r.query_length) > 0 {
            continue;
        }
        intron_signal += r.query_length as f64;
    }
    intron_signal
}

/// Computes per-position coverage at each DISTINCT position in
/// `positions` (1-based; may contain duplicates, which are naturally
/// collapsed here exactly as upstream's per-genomic-position pileup
/// column iteration collapses them -- callers must still use the
/// ORIGINAL, possibly-duplicated `positions.len()` as the TIN
/// denominator, not this function's output length; see `tin_score`).
/// Computes coverage for requested positions and also returns which positions
/// were visited by at least one read (to distinguish truly uncovered positions
/// from covered-but-filtered ones, needed for Python duck-typing float marker).
pub fn genebody_coverage_with_visited(reads: &[IndexedRead], positions: &[i64], bg_level: f64) -> (Vec<f64>, Vec<bool>) {
    if positions.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let start0 = positions[0] - 1;
    let end0 = positions[positions.len() - 1];

    // Reads entering pysam's default ("samtools" stepper) pileup over the
    // region: flag_filter drops duplicates (unmapped/secondary/qcfail are
    // already excluded at index time) and ignore_orphans drops reads that
    // are paired but not properly paired.
    let overlapping: Vec<&IndexedRead> = reads
        .iter()
        .filter(|r| r.start < end0 && r.end > start0 && !r.is_duplicate && !(r.is_paired && !r.is_proper_pair))
        .collect();

    // ignore_overlaps=True: htslib rewrites the base qualities of
    // overlapping proper-pair mates as each pair completes in the pileup
    // buffer (sam.c overlap_push/tweak_overlap_quality); min_base_quality
    // then applies to the rewritten values.
    let mut qualities: Vec<Vec<u8>> = overlapping.iter().map(|r| r.qualities.clone()).collect();
    let mut waiting: HashMap<&str, usize> = HashMap::new();
    for (i, r) in overlapping.iter().enumerate() {
        if r.mate_unmapped || !r.is_proper_pair {
            continue;
        }
        if r.mate_on_other_reference || (r.template_length.abs() >= 2 * r.query_length && r.mate_start >= r.end) {
            continue;
        }
        match waiting.remove(r.query_name.as_str()) {
            Some(j) if overlapping[j].end > r.start => {
                let (lo, hi) = qualities.split_at_mut(i);
                tweak_overlap_quality(overlapping[j], r, &mut lo[j], &mut hi[0]);
            }
            // The earlier mate already left the pileup buffer
            // (overlap_remove) before this one arrived: nothing to adjust.
            Some(_) => {}
            None => {
                if r.mate_start >= r.start || (r.is_paired && r.mate_start == -1) {
                    waiting.insert(r.query_name.as_str(), i);
                }
            }
        }
    }

    let mut distinct: Vec<i64> = positions.to_vec();
    distinct.sort_unstable();
    distinct.dedup();

    let mut coverage = Vec::with_capacity(distinct.len());
    let mut visited = Vec::with_capacity(distinct.len());
    for p1 in distinct {
        let p0 = p1 - 1;
        let mut covered = 0.0;
        let mut was_visited = false;
        // pysam's pileup() default max_depth is 8000.  The cap is applied
        // per genomic column before RSeQC filters duplicate/low-quality
        // pileups, so a deeply stacked locus must not silently contribute
        // an unbounded count in this in-memory implementation.
        let mut pileup_depth = 0usize;
        for (i, read) in overlapping.iter().enumerate() {
            let mut base_quality = None;
            for &(bs, be, qstart) in &read.match_blocks {
                if p0 >= bs && p0 < be {
                    let qidx = qstart + (p0 - bs) as usize;
                    base_quality = Some(qualities[i].get(qidx).copied().unwrap_or(0));
                    break;
                }
            }
            let in_gap = base_quality.is_none() && read.skip_delete_blocks.iter().any(|&(bs, be)| p0 >= bs && p0 < be);
            if base_quality.is_none() && !in_gap {
                continue;
            }
            was_visited = true;
            if pileup_depth >= 8000 {
                continue;
            }
            pileup_depth += 1;
            if base_quality.is_some_and(|q| q >= 13) {
                covered += 1.0;
            }
        }
        coverage.push(covered);
        visited.push(was_visited);
    }

    if bg_level <= 0.0 {
        return (coverage, visited);
    }
    let adjusted = coverage
        .into_iter()
        .map(|v| {
            let subtracted = (v - bg_level) as i64;
            if subtracted > 0 { subtracted as f64 } else { 0.0 }
        })
        .collect();
    (adjusted, visited)
}

pub fn genebody_coverage(reads: &[IndexedRead], positions: &[i64], bg_level: f64) -> Vec<f64> {
    genebody_coverage_with_visited(reads, positions, bg_level).0
}

/// BAM CIGAR operation code (`MIDNSHP=X` -> 0..=8).
fn cigar_code(kind: Kind) -> u8 {
    match kind {
        Kind::Match => 0,
        Kind::Insertion => 1,
        Kind::Deletion => 2,
        Kind::Skip => 3,
        Kind::SoftClip => 4,
        Kind::HardClip => 5,
        Kind::Pad => 6,
        Kind::SequenceMatch => 7,
        Kind::SequenceMismatch => 8,
    }
}

const CIG_MATCH: i32 = 0;

fn is_match_op(op: u8) -> bool {
    op == 0 || op == 7 || op == 8
}

/// Port of htslib `cigar_iref2iseq_set` (sam.c).
fn cigar_iref2iseq_set(cigar: &[(u8, i64)], ci: &mut usize, icig: &mut i64, iseq: &mut i64, iref: &mut i64) -> i32 {
    let mut pos = *iref;
    if pos < 0 {
        return -1;
    }
    *icig = 0;
    *iseq = 0;
    *iref = 0;
    while *ci < cigar.len() {
        let (cig, ncig) = cigar[*ci];
        match cig {
            4 | 1 => {
                *ci += 1;
                *iseq += ncig;
                *icig = 0;
            }
            5 | 6 => {
                *ci += 1;
                *icig = 0;
            }
            2 | 3 => {
                pos -= ncig;
                if pos < 0 {
                    pos = 0;
                }
                *ci += 1;
                *icig = 0;
                *iref += ncig;
            }
            _ if is_match_op(cig) => {
                pos -= ncig;
                if pos < 0 {
                    *icig = ncig + pos;
                    *iseq += *icig;
                    *iref += *icig;
                    return CIG_MATCH;
                }
                *ci += 1;
                *iseq += ncig;
                *icig = 0;
                *iref += ncig;
            }
            _ => return -2,
        }
    }
    *iseq = -1;
    -1
}

/// Port of htslib `cigar_iref2iseq_next` (sam.c).
fn cigar_iref2iseq_next(cigar: &[(u8, i64)], ci: &mut usize, icig: &mut i64, iseq: &mut i64, iref: &mut i64) -> i32 {
    while *ci < cigar.len() {
        let (cig, ncig) = cigar[*ci];
        match cig {
            2 | 3 => {
                *ci += 1;
                *iref += ncig;
                *icig = -1;
            }
            1 | 4 => {
                *ci += 1;
                *iseq += ncig;
                *icig = -1;
            }
            5 | 6 => {
                *ci += 1;
                *icig = -1;
            }
            _ if is_match_op(cig) => {
                if *icig >= ncig - 1 {
                    *icig = -1;
                    *ci += 1;
                    continue;
                }
                *iseq += 1;
                *icig += 1;
                *iref += 1;
                return CIG_MATCH;
            }
            _ => return -2,
        }
    }
    *iseq = -1;
    *iref = -1;
    -1
}

/// khash `__ac_X31_hash_string` followed by `__ac_Wang_hash` (C `char` is
/// signed on the platforms htslib is built for, hence the sign extension).
fn htslib_name_hash(name: &str) -> u32 {
    let bytes = name.as_bytes();
    let mut h: u32 = match bytes.first() {
        Some(&c) => c as i8 as i32 as u32,
        None => 0,
    };
    if h != 0 {
        for &c in &bytes[1..] {
            h = (h << 5).wrapping_sub(h).wrapping_add(c as i8 as i32 as u32);
        }
    }
    let mut key = h;
    key = key.wrapping_add(!(key << 15));
    key ^= key >> 10;
    key = key.wrapping_add(key << 3);
    key ^= key >> 6;
    key = key.wrapping_add(!(key << 11));
    key ^= key >> 16;
    key
}

/// `(uint8_t)(0.8 * q)` -- C double-to-uint8 truncation.
fn scale08(mul: u8, q: u8) -> u8 {
    (f64::from(mul) * 0.8 * f64::from(q)) as u8
}

/// Port of htslib `tweak_overlap_quality` (sam.c, htslib >= 1.13, as
/// bundled with the oracle's pysam): given overlapping proper-pair mates
/// `a` (earlier in the pileup) and `b`, rewrite their base qualities in the
/// overlap so that only one mate keeps a (summed or 0.8-scaled) quality.
fn tweak_overlap_quality(a: &IndexedRead, b: &IndexedRead, a_qual: &mut [u8], b_qual: &mut [u8]) {
    let (mut a_ci, mut b_ci) = (0usize, 0usize);
    let (mut a_icig, mut a_iseq, mut b_icig, mut b_iseq) = (0i64, 0i64, 0i64, 0i64);
    let mut iref = b.start;
    let mut a_iref = iref - a.start;
    let mut b_iref = iref - b.start;
    let mut a_ret = cigar_iref2iseq_set(&a.cigar, &mut a_ci, &mut a_icig, &mut a_iseq, &mut a_iref);
    if a_ret < 0 {
        return;
    }
    let mut b_ret = cigar_iref2iseq_set(&b.cigar, &mut b_ci, &mut b_icig, &mut b_iseq, &mut b_iref);
    if b_ret < 0 {
        return;
    }
    let (amul, bmul) = if htslib_name_hash(&a.query_name) & 1 == 1 { (1u8, 0u8) } else { (0u8, 1u8) };
    let (a_len, b_len) = (a_qual.len() as i64, b_qual.len() as i64);
    loop {
        while a_ret >= 0 && a_iref >= 0 && a_iref < iref - a.start {
            a_ret = cigar_iref2iseq_next(&a.cigar, &mut a_ci, &mut a_icig, &mut a_iseq, &mut a_iref);
        }
        if a_ret < 0 {
            break;
        }
        while b_ret >= 0 && b_iref >= 0 && b_iref < iref - b.start {
            b_ret = cigar_iref2iseq_next(&b.cigar, &mut b_ci, &mut b_icig, &mut b_iseq, &mut b_iref);
        }
        if b_ret < 0 {
            break;
        }
        if iref < a_iref + a.start {
            iref = a_iref + a.start;
        }
        if iref < b_iref + b.start {
            iref = b_iref + b.start;
        }
        iref += 1;

        if a_iref + a.start != b_iref + b.start {
            if a_iref + a.start < b_iref + b.start && b_ci > 0 && b.cigar[b_ci - 1].0 == 2 {
                loop {
                    if a_iseq < 0 || a_iseq >= a_len {
                        return;
                    }
                    let q = &mut a_qual[a_iseq as usize];
                    *q = if amul == 1 { scale08(1, *q) } else { 0 };
                    a_ret = cigar_iref2iseq_next(&a.cigar, &mut a_ci, &mut a_icig, &mut a_iseq, &mut a_iref);
                    if a_ret < 0 {
                        return;
                    }
                    if a_iref + a.start >= b_iref + b.start {
                        break;
                    }
                }
            } else if a_ci > 0 && a.cigar[a_ci - 1].0 == 2 {
                loop {
                    if b_iseq < 0 || b_iseq >= b_len {
                        return;
                    }
                    let q = &mut b_qual[b_iseq as usize];
                    *q = if bmul == 1 { scale08(1, *q) } else { 0 };
                    b_ret = cigar_iref2iseq_next(&b.cigar, &mut b_ci, &mut b_icig, &mut b_iseq, &mut b_iref);
                    if b_ret < 0 {
                        return;
                    }
                    if b_iref + b.start >= a_iref + a.start {
                        break;
                    }
                }
            } else {
                // e.g. ref-skip: not supported by htslib here either
                continue;
            }
        }

        if a_iseq > a.query_length || b_iseq > b.query_length {
            return;
        }
        // Guard the C out-of-bounds read at iseq == l_qseq (bad CIGAR).
        if a_iseq < 0 || b_iseq < 0 || a_iseq >= a_len || b_iseq >= b_len {
            return;
        }
        let (ai, bi) = (a_iseq as usize, b_iseq as usize);
        let same = match (a.sequence.get(ai), b.sequence.get(bi)) {
            (Some(x), Some(y)) => x == y,
            _ => false,
        };
        if same {
            let qual = (u32::from(a_qual[ai]) + u32::from(b_qual[bi])).min(200) as u8;
            a_qual[ai] = amul * qual;
            b_qual[bi] = bmul * qual;
        } else if a_qual[ai] > b_qual[bi] {
            a_qual[ai] = scale08(1, a_qual[ai]);
            b_qual[bi] = 0;
        } else if a_qual[ai] < b_qual[bi] {
            b_qual[bi] = scale08(1, b_qual[bi]);
            a_qual[ai] = 0;
        } else {
            a_qual[ai] = scale08(amul, a_qual[ai]);
            b_qual[bi] = scale08(bmul, b_qual[bi]);
        }
    }
}

/// Shannon entropy (natural log), matching upstream's `shannon_entropy`
/// exactly, including its `0.0`-not-`-0.0` special case.
pub fn shannon_entropy(values: &[f64]) -> f64 {
    let total: f64 = values.iter().sum();
    if total <= 0.0 {
        return 0.0;
    }
    let mut entropy = 0.0;
    for &value in values {
        let probability = value / total;
        entropy += probability * probability.ln();
    }
    if entropy == 0.0 { 0.0 } else { -entropy }
}

/// Ports `tin_score`'s historical entropy formula.
pub fn tin_score(coverage: &[f64], sampled_length: i64) -> f64 {
    if coverage.is_empty() || sampled_length <= 0 {
        return 0.0;
    }
    let effective: Vec<f64> = coverage.iter().copied().filter(|&v| v > 0.0).collect();
    if effective.is_empty() {
        return 0.0;
    }
    let entropy = shannon_entropy(&effective);
    100.0 * entropy.exp() / sampled_length as f64
}

#[derive(Debug, Clone, PartialEq)]
pub struct TranscriptSample {
    pub gene_name: String,
    pub chrom: String,
    pub tx_start: i64,
    pub tx_end: i64,
    pub intron_size: i64,
    /// 1-based positions; may contain a duplicate (small-transcript
    /// branch) and is NOT deduplicated -- see module docs.
    pub chosen_bases: Vec<i64>,
}

/// Parses the refgene BED12 and yields sampled transcript positions.
/// Robustness tier matches upstream exactly: blank lines and
/// comment/track/browser lines are skipped silently; any other per-line
/// parse failure (fewer than 12 columns, or any `int(...)` failure) is
/// caught and skipped with a stderr warning -- upstream wraps the WHOLE
/// per-line parse in one broad `except (ValueError, IndexError)`.
pub fn genomic_positions(reader: impl BufRead, sample_size: i64) -> io::Result<Vec<TranscriptSample>> {
    let mut out = Vec::new();

    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        let line_number = line_number + 1;
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        let parsed: Result<_, String> = (|| {
            if fields.len() < 12 {
                return Err("fewer than 12 columns".to_string());
            }
            let chromosome = fields[0].to_string();
            let tx_start: i64 = fields[1].parse().map_err(|_| "invalid txStart".to_string())?;
            let tx_end: i64 = fields[2].parse().map_err(|_| "invalid txEnd".to_string())?;
            let gene_name = fields[3].to_string();

            let block_sizes: Vec<i64> = fields[10]
                .trim_end_matches(',')
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| s.parse().map_err(|_| "invalid blockSize".to_string()))
                .collect::<Result<_, _>>()?;
            let block_starts: Vec<i64> = fields[11]
                .trim_end_matches(',')
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| s.parse().map_err(|_| "invalid blockStart".to_string()))
                .collect::<Result<_, _>>()?;

            let exon_starts: Vec<i64> = block_starts.iter().map(|&s| tx_start + s).collect();
            let exon_ends: Vec<i64> = exon_starts.iter().zip(block_sizes.iter()).map(|(&s, &sz)| s + sz).collect();
            let mrna_size: i64 = block_sizes.iter().sum();
            let intron_size = (tx_end - tx_start - mrna_size).max(0);

            Ok((chromosome, tx_start, tx_end, gene_name, exon_starts, exon_ends, mrna_size, intron_size))
        })();

        let (chromosome, tx_start, tx_end, gene_name, exon_starts, exon_ends, mrna_size, intron_size) = match parsed
        {
            Ok(v) => v,
            Err(msg) => {
                eprintln!("Skipping malformed BED12 line {line_number}: {msg}");
                continue;
            }
        };

        let chosen_bases: Vec<i64> = if mrna_size <= sample_size {
            let mut bases = vec![tx_start + 1, tx_end];
            for (&s, &e) in exon_starts.iter().zip(exon_ends.iter()) {
                bases.extend((s + 1)..=e);
            }
            bases
        } else {
            let step_size = std::cmp::max(1, mrna_size / sample_size);
            let mut all_bases: Vec<i64> = Vec::new();
            let mut exon_bounds: Vec<i64> = Vec::new();
            for (&s, &e) in exon_starts.iter().zip(exon_ends.iter()) {
                all_bases.extend((s + 1)..=e);
                exon_bounds.push(s + 1);
                exon_bounds.push(e);
            }
            let mut sampled_bases = Vec::new();
            let mut idx = 0usize;
            while idx < all_bases.len() {
                sampled_bases.push(all_bases[idx]);
                idx += step_size as usize;
            }
            let mut seen = HashSet::new();
            exon_bounds.into_iter().chain(sampled_bases).filter(|v| seen.insert(*v)).collect()
        };

        out.push(TranscriptSample { gene_name, chrom: chromosome, tx_start, tx_end, intron_size, chosen_bases });
    }

    Ok(out)
}

/// Builds a chromosome-keyed exon interval index for `--subtract-background`.
/// No case-folding: unlike `split_bam.py`/`junction_annotation.py`, tin.py
/// never uppercases chromosome names anywhere.
pub fn build_exon_ranges(exons: &[Bed3]) -> MergedRegions {
    MergedRegions::new(exons)
}

#[derive(Debug, Clone, PartialEq)]
pub struct TinRecord {
    pub gene_name: String,
    pub chrom: String,
    pub tx_start: i64,
    pub tx_end: i64,
    pub score: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct TinSummary {
    pub mean: f64,
    pub median: f64,
    pub stdev: f64,
}

/// `np.mean` (pairwise summation).
fn mean(values: &[f64]) -> f64 {
    numpy_mean(values)
}

fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = sorted.len();
    if n % 2 == 1 { sorted[n / 2] } else { (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0 }
}

/// `np.std` (ddof=0, pairwise summation).
fn population_stdev(values: &[f64]) -> f64 {
    numpy_std(values, 0)
}

/// Scores one transcript against the reads that can reach it. Shared by the
/// whole-file (`compute_tin`) and windowed (`compute_tin_windowed`) drivers so
/// both do byte-identical arithmetic -- the drivers differ only in WHICH reads
/// are resident, never in how a resident read is scored.
///
/// Returns `None` when the transcript fails the minimum-coverage criterion;
/// that transcript scores 0.0 and is excluded from the summary statistics,
/// matching `tin.py`.
fn score_sample(
    s: &TranscriptSample,
    reads: &[IndexedRead],
    min_cov: i64,
    exon_ranges: Option<&MergedRegions>,
) -> Option<f64> {
    if !check_min_reads(reads, s.tx_start, s.tx_end, min_cov) {
        return None;
    }
    let mut noise_level = 0.0;
    if let Some(ranges) = exon_ranges {
        if s.intron_size > 0 {
            let intron_signal = estimate_bg_noise(reads, s.tx_start, s.tx_end, ranges, &s.chrom);
            noise_level = intron_signal / s.intron_size as f64;
        }
    }
    let mut positions = s.chosen_bases.clone();
    positions.sort_unstable();
    let coverage = genebody_coverage(reads, &positions, noise_level);
    Some(tin_score(&coverage, s.chosen_bases.len() as i64))
}

/// Computes TIN for every sampled transcript against one BAM's read
/// index. Ports `process_bam`'s per-transcript loop (output writing is
/// left to the caller/CLI).
///
/// Holds the whole BAM in memory; `compute_tin_windowed` computes the same
/// numbers with a sliding window instead and is what the CLI uses.
pub fn compute_tin(
    samples: &[TranscriptSample],
    reads_by_chrom: &HashMap<String, Vec<IndexedRead>>,
    min_cov: i64,
    exon_ranges: Option<&MergedRegions>,
) -> (Vec<TinRecord>, TinSummary) {
    let empty = Vec::new();
    let mut records = Vec::with_capacity(samples.len());
    let mut sample_tins: Vec<f64> = Vec::new();

    for s in samples {
        let reads = reads_by_chrom.get(&s.chrom).unwrap_or(&empty);
        let score = match score_sample(s, reads, min_cov, exon_ranges) {
            Some(score) => {
                sample_tins.push(score);
                score
            }
            None => 0.0,
        };

        records.push(TinRecord { gene_name: s.gene_name.clone(), chrom: s.chrom.clone(), tx_start: s.tx_start, tx_end: s.tx_end, score });
    }

    (records, summarize(&sample_tins))
}

/// Assembles the summary from scored samples, in sample order.
fn summarize(sample_tins: &[f64]) -> TinSummary {
    if sample_tins.is_empty() {
        TinSummary::default()
    } else {
        TinSummary { mean: mean(sample_tins), median: median(sample_tins), stdev: population_stdev(sample_tins) }
    }
}

/// Outcome of the windowed driver.
pub enum WindowedTin {
    Computed(Vec<TinRecord>, TinSummary),
    /// The input was not coordinate-sorted, so a sliding window cannot be
    /// maintained. The caller should fall back to the whole-file
    /// `compute_tin`, which sorts the per-chromosome read list itself.
    NotCoordinateSorted,
}

/// Computes the same TIN numbers as [`compute_tin`] while holding only a
/// **sliding window** of reads in memory instead of the whole BAM.
///
/// # Why this exists
///
/// `compute_tin` needs, for each transcript, exactly the reads whose *start*
/// lies in `[tx_start, tx_end)` -- that is what [`reads_starting_in`] returns.
/// A whole-file index therefore retains every read in the BAM forever, and
/// because each `IndexedRead` carries per-read heap buffers (query name,
/// qualities, sequence, CIGAR, block lists) that costs ~600 bytes per read:
/// a 600k-read BAM measured 366 MB against upstream `pysam`'s 43 MB, because
/// `pysam` answers each region query from a BAI index and never materialises
/// the rest of the file.
///
/// This driver inverts that. It walks transcripts in **coordinate order** and
/// streams the BAM **once**, keeping only the reads that can still reach a
/// transcript not yet scored. Two facts make that exact rather than an
/// approximation:
///
/// 1. A coordinate-sorted BAM emits reads in increasing start order within a
///    reference sequence, and each reference sequence in one contiguous block.
///    So the records still to come are all at or after the current transcript's
///    end, and a read that starts after it can never be scored.
/// 2. A read with `end <= tx_start` has `start < end <= tx_start`, so it falls
///    before the next (coordinate-ordered) transcript's window and
///    `reads_starting_in` would exclude it. Dropping it early is therefore
///    indistinguishable from keeping it.
///
/// On the 600k-read benchmark workload the resident set falls from 597,048
/// reads to a mean of 3,424 (max 16,438) -- a 174x reduction, and the reason
/// this port's `tin` peak RSS now beats upstream's instead of losing to it
/// 8.5x.
///
/// Transcript scores are collected by original input index, so both the
/// `.tin.xls` row order and the summary's summation order (which is
/// pairwise-summation-order dependent) are identical to `compute_tin`'s.
pub fn compute_tin_windowed<I>(
    records: I,
    header: &sam::Header,
    samples: &[TranscriptSample],
    min_cov: i64,
    exon_ranges: Option<&MergedRegions>,
) -> io::Result<WindowedTin>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    // Coordinate order for samples: by the header's reference-sequence order
    // (which is the order a coordinate-sorted BAM emits them in), then start.
    // Samples on a chromosome absent from the header sort last and score 0.0,
    // exactly as `reads_by_chrom.get(...).unwrap_or(&empty)` does today.
    let ref_ids: HashMap<String, usize> = header
        .reference_sequences()
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (name.to_string(), i))
        .collect();
    let ref_index = |chrom: &str| -> usize { ref_ids.get(chrom).copied().unwrap_or(usize::MAX) };
    let mut order: Vec<usize> = (0..samples.len()).collect();
    order.sort_by_key(|&i| (ref_index(&samples[i].chrom), samples[i].tx_start, samples[i].tx_end, i));

    let mut scores: Vec<Option<f64>> = vec![None; samples.len()];
    let mut window: Vec<IndexedRead> = Vec::new();
    let mut window_ref: Option<usize> = None;
    // Reference sequences already fully consumed, so a record arriving for one
    // of them later proves the input is not coordinate-sorted.
    let mut closed: Vec<usize> = Vec::new();
    let mut last_start: i64 = i64::MIN;
    let mut records = records.into_iter().peekable();

    for &si in &order {
        let s = &samples[si];
        let s_ref = ref_index(&s.chrom);

        // A change of reference sequence invalidates the whole window: the
        // reads in it can no longer be reached by a later transcript.
        if window_ref != Some(s_ref) {
            if let Some(r) = window_ref {
                closed.push(r);
            }
            window = Vec::new();
            window_ref = if s_ref == usize::MAX { None } else { Some(s_ref) };
            last_start = i64::MIN;
        }

        // Pull every record that starts before this transcript ends. Records
        // for a later reference sequence are left unconsumed for their own
        // sample to pick up.
        while let Some(peeked) = records.peek() {
            let next = match peeked {
                Ok(record) => record,
                // A decode error is left for the drain at the end to surface.
                Err(_) => break,
            };
            // The record's reference id IS the index into the header's
            // reference-sequence list, i.e. its coordinate-sort key.
            let Some(next_ref) = next.reference_sequence_id().transpose()? else { break };
            if header.reference_sequences().get_index(next_ref).is_none() {
                break;
            }
            if next_ref != s_ref {
                break;
            }
            if closed.contains(&next_ref) {
                // This reference sequence's block already ended.
                return Ok(WindowedTin::NotCoordinateSorted);
            }
            let Some(pos) = next.alignment_start().transpose()? else { break };
            let start = (pos.get() - 1) as i64;
            if start >= s.tx_end {
                break;
            }
            if start < last_start {
                return Ok(WindowedTin::NotCoordinateSorted);
            }
            last_start = start;
            // Convert while the peeked borrow is alive, then consume it.
            let read = to_indexed_read(next)?;
            records.next();
            if let Some(read) = read {
                window.push(read);
            }
        }

        // Retire reads that can no longer reach this or any later transcript.
        // Safe by the `end <= tx_start` argument in the doc comment; the
        // retained reads stay sorted by start, which `reads_starting_in`'s
        // `partition_point` requires.
        window.retain(|r| r.end > s.tx_start);

        if window_ref.is_some() {
            scores[si] = score_sample(s, &window, min_cov, exon_ranges);
        }
    }

    // Drain the remainder so a decode error after the last transcript still
    // surfaces rather than being silently dropped.
    for result in records {
        result?;
    }

    let mut sample_tins = Vec::new();
    let mut out = Vec::with_capacity(samples.len());
    for (i, s) in samples.iter().enumerate() {
        let score = match scores[i] {
            Some(score) => {
                sample_tins.push(score);
                score
            }
            None => 0.0,
        };
        out.push(TinRecord { gene_name: s.gene_name.clone(), chrom: s.chrom.clone(), tx_start: s.tx_start, tx_end: s.tx_end, score });
    }

    Ok(WindowedTin::Computed(out, summarize(&sample_tins)))
}

/// Ports the `.tin.xls` output: a header line followed by one row per
/// transcript, in input order.
pub fn render_tin_xls(records: &[TinRecord]) -> String {
    let mut out = String::from("geneID\tchrom\ttx_start\ttx_end\tTIN\n");
    for r in records {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\n",
            r.gene_name,
            r.chrom,
            r.tx_start,
            r.tx_end,
            python_str_float(r.score)
        ));
    }
    out
}

/// Ports the `.summary.txt` output. `bam_name` is the bare file name
/// (`Path.name`), matching upstream.
pub fn render_summary(bam_name: &str, summary: &TinSummary) -> String {
    format!(
        "Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)\n{}\t{}\t{}\t{}\n",
        bam_name,
        python_str_float(summary.mean),
        python_str_float(summary.median),
        python_str_float(summary.stdev),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn genomic_positions_small_transcript_keeps_boundary_duplicate() {
        // Single-exon transcript chr1:0-100 (mrna_size=100 <= sample_size=100):
        // chosen_bases = [1, 100] + range(1,101) -- 1 appears twice (once
        // from the initial boundary pair, once from the exon range).
        let bed = "chr1\t0\t100\ttx1\t0\t+\t0\t100\t0\t1\t100,\t0,\n";
        let samples = genomic_positions(Cursor::new(bed), 100).unwrap();
        assert_eq!(samples.len(), 1);
        let s = &samples[0];
        assert_eq!(s.gene_name, "tx1");
        assert_eq!(s.intron_size, 0);
        assert_eq!(s.chosen_bases.len(), 102); // [1,100] + 100 positions (1..=100)
        assert_eq!(s.chosen_bases.iter().filter(|&&v| v == 1).count(), 2);
    }

    #[test]
    fn genomic_positions_large_transcript_dedups_and_drops_boundary_pair() {
        // Two-exon transcript chr1:0-1000 with a gap between the last
        // exon and tx_end (exon2 ends at 900, tx_end is 1000), so
        // tx_end=1000 is NOT among the exon boundaries naturally: mrna_size
        // =200 > sample_size=10, step_size = max(1, 200/10) = 20. The
        // initial [tx_start+1, tx_end] = [1, 1000] pair is NOT present
        // (else-branch reassigns chosen_bases entirely); exon_bounds are
        // the per-exon boundaries instead: (1,100) and (801,900).
        let bed = "chr1\t0\t1000\ttx1\t0\t+\t0\t1000\t0\t2\t100,100,\t0,800,\n";
        let samples = genomic_positions(Cursor::new(bed), 10).unwrap();
        let s = &samples[0];
        assert_eq!(s.intron_size, 800);
        assert!(!s.chosen_bases.contains(&1000)); // tx_end itself never appears
        assert!(s.chosen_bases.contains(&1));
        assert!(s.chosen_bases.contains(&100));
        assert!(s.chosen_bases.contains(&801));
        assert!(s.chosen_bases.contains(&900));
        // No duplicates in the else-branch (uniqify applied).
        let mut seen = HashSet::new();
        assert!(s.chosen_bases.iter().all(|v| seen.insert(*v)));
    }

    #[test]
    fn genomic_positions_skips_malformed_lines_with_warning_not_error() {
        let bed = "\
track name=x
# comment

not enough fields
chr1\t0\tnotanumber\ttx1\t0\t+\t0\t100\t0\t1\t100,\t0,
chr1\t0\t100\ttx2\t0\t+\t0\t100\t0\t1\t100,\t0,
";
        let samples = genomic_positions(Cursor::new(bed), 100).unwrap();
        assert_eq!(samples.len(), 1);
        assert_eq!(samples[0].gene_name, "tx2");
    }

    #[test]
    fn shannon_entropy_matches_python_hand_computation() {
        // Verified against python3: -sum(p*ln(p)) for p=[0.5,0.5] = ln(2).
        let e = shannon_entropy(&[1.0, 1.0]);
        assert!((e - std::f64::consts::LN_2).abs() < 1e-12);

        // Single value -> p=1.0, ln(1)=0 -> the 0.0-not-negative-zero branch.
        let e2 = shannon_entropy(&[5.0]);
        assert_eq!(e2, 0.0);
        assert!(e2.is_sign_positive());
    }

    #[test]
    fn tin_score_zero_when_no_positive_coverage() {
        assert_eq!(tin_score(&[0.0, 0.0], 10), 0.0);
        assert_eq!(tin_score(&[], 10), 0.0);
        assert_eq!(tin_score(&[1.0], 0), 0.0);
    }

    #[test]
    fn tin_score_uniform_coverage_gives_100() {
        // Uniform coverage across every sampled position -> maximum entropy
        // (ln(n)) -> exp(ln(n))=n -> 100*n/n = 100.0 exactly.
        let coverage = vec![5.0; 20];
        let score = tin_score(&coverage, 20);
        assert!((score - 100.0).abs() < 1e-9);
    }

    fn header_with_chrom(name: &str, len: usize) -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence(
                name,
                sam::header::record::value::Map::<sam::header::record::value::map::ReferenceSequence>::new(
                    std::num::NonZeroUsize::new(len).unwrap(),
                ),
            )
            .build()
    }

    fn to_bam_records(header: &sam::Header, records: &[sam::alignment::record_buf::RecordBuf]) -> Vec<bam::Record> {
        use sam::alignment::io::Write as _;
        let mut buf = Vec::new();
        {
            let mut writer = bam::io::Writer::new(&mut buf);
            writer.write_header(header).unwrap();
            for r in records {
                writer.write_alignment_record(header, r).unwrap();
            }
        }
        let mut reader = bam::io::Reader::new(buf.as_slice());
        reader.read_header().unwrap();
        reader.records().map(|r| r.unwrap()).collect()
    }

    #[test]
    fn build_read_index_and_check_min_reads() {
        use sam::alignment::record::{Flags, MappingQuality, cigar::Op};
        use sam::alignment::record_buf::{Cigar, RecordBuf};

        let header = header_with_chrom("chr1", 1000);
        let mut records = Vec::new();
        for i in 0..3 {
            records.push(
                RecordBuf::builder()
                    .set_name(format!("r{i}"))
                    .set_flags(Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_alignment_start(noodles_core::Position::try_from(101 + i).unwrap())
                    .set_mapping_quality(MappingQuality::new(40).unwrap())
                    .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
                    .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 30]))
                    .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![40; 30]))
                    .build(),
            );
        }

        let bam_records = to_bam_records(&header, &records);
        let index = build_read_index(bam_records.into_iter().map(Ok), &header).unwrap();
        let reads = &index["chr1"];
        assert_eq!(reads.len(), 3);

        // 3 distinct read starts (100,101,102 0-based); cutoff=2 (strict >).
        assert!(check_min_reads(reads, 0, 200, 2));
        assert!(!check_min_reads(reads, 0, 200, 3));
    }

    #[test]
    fn genebody_coverage_counts_matching_reads_and_excludes_low_quality_and_duplicates() {
        use sam::alignment::record::{Flags, MappingQuality, cigar::Op};
        use sam::alignment::record_buf::{Cigar, RecordBuf};

        let header = header_with_chrom("chr1", 1000);

        // Read A: pos 0-based 0, 10M, quality 40 everywhere -> covers 1-based [1,10].
        let read_a = RecordBuf::builder()
            .set_name("a")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(1).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 10]))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![40; 10]))
            .build();

        // Read B: same span, but flagged duplicate -> excluded from coverage
        // (pileup's flag_filter) even though it would be in the read index.
        let read_b = RecordBuf::builder()
            .set_name("b")
            .set_flags(Flags::DUPLICATE)
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(1).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 10]))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![40; 10]))
            .build();

        // Read C: same span, but quality 5 everywhere (< 13) -> excluded per-base.
        let read_c = RecordBuf::builder()
            .set_name("c")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(1).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 10]))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![5; 10]))
            .build();

        let bam_records = to_bam_records(&header, &[read_a, read_b, read_c]);
        let index = build_read_index(bam_records.into_iter().map(Ok), &header).unwrap();
        let reads = &index["chr1"];
        assert_eq!(reads.len(), 3); // all 3 pass the qcfail/unmapped/secondary build filter

        let positions = vec![1, 5, 10];
        let coverage = genebody_coverage(reads, &positions, 0.0);
        // Only read_a counts at every position: read_b is a duplicate
        // (excluded), read_c's bases are all below min_base_quality=13.
        assert_eq!(coverage, vec![1.0, 1.0, 1.0]);
    }

    #[test]
    fn genebody_coverage_honors_pysam_default_max_depth() {
        let reads: Vec<IndexedRead> = (0..8001)
            .map(|_| IndexedRead {
                query_name: "synthetic".to_string(),
                is_paired: false,
                start: 0,
                end: 20,
                is_duplicate: false,
                match_blocks: vec![(0, 20, 0)],
                skip_delete_blocks: vec![],
                qualities: vec![40; 20],
                query_length: 20,
                ..Default::default()
            })
            .collect();

        // pysam.AlignmentFile.pileup() defaults max_depth to 8000.
        assert_eq!(genebody_coverage(&reads, &[1], 0.0), vec![8000.0]);
    }

    #[test]
    fn genebody_coverage_deduplicates_overlapping_mates() {
        let reads = vec![
            IndexedRead {
                query_name: "pair".to_string(),
                is_paired: true,
                start: 0,
                end: 20,
                is_duplicate: false,
                match_blocks: vec![(0, 20, 0)],
                skip_delete_blocks: vec![],
                qualities: vec![40; 20],
                query_length: 20,
                is_proper_pair: true,
                mate_start: 10,
                template_length: 30,
                sequence: vec![b'A'; 20],
                cigar: vec![(0, 20)],
                ..Default::default()
            },
            IndexedRead {
                query_name: "pair".to_string(),
                is_paired: true,
                start: 10,
                end: 30,
                is_duplicate: false,
                match_blocks: vec![(10, 30, 0)],
                skip_delete_blocks: vec![],
                qualities: vec![40; 20],
                query_length: 20,
                is_proper_pair: true,
                mate_start: 0,
                template_length: -30,
                sequence: vec![b'A'; 20],
                cigar: vec![(0, 20)],
                ..Default::default()
            },
        ];

        // The overlap at reference position 11 is one pileup base, not two.
        assert_eq!(genebody_coverage(&reads, &[11], 0.0), vec![1.0]);
    }

    fn pair_read(name: &str, start: i64, mate_start: i64, base: u8, qual: u8, proper: bool) -> IndexedRead {
        IndexedRead {
            query_name: name.to_string(),
            is_paired: true,
            start,
            end: start + 10,
            match_blocks: vec![(start, start + 10, 0)],
            qualities: vec![qual; 10],
            query_length: 10,
            is_proper_pair: proper,
            mate_start,
            template_length: if start <= mate_start { 15 } else { -15 },
            sequence: vec![base; 10],
            cigar: vec![(0, 10)],
            ..Default::default()
        }
    }

    #[test]
    fn genebody_coverage_mirrors_htslib_overlap_quality_rules() {
        // Probed against pysam 0.24.1 pileup(): matching overlapping bases
        // keep the SUM of both qualities on one mate (8 + 8 = 16 >= 13);
        // mismatching ones keep int(0.8 * max) (0.8 * 16 = 12 < 13, 0.8 * 17
        // = 13).
        let same = [pair_read("p", 0, 5, b'A', 8, true), pair_read("p", 5, 0, b'A', 8, true)];
        assert_eq!(genebody_coverage(&same, &[8], 0.0), vec![1.0]);
        let mism16 = [pair_read("p", 0, 5, b'A', 16, true), pair_read("p", 5, 0, b'C', 10, true)];
        assert_eq!(genebody_coverage(&mism16, &[8], 0.0), vec![0.0]);
        let mism17 = [pair_read("p", 0, 5, b'A', 17, true), pair_read("p", 5, 0, b'C', 10, true)];
        assert_eq!(genebody_coverage(&mism17, &[8], 0.0), vec![1.0]);
    }

    #[test]
    fn genebody_coverage_ignores_orphans() {
        // pysam's default "samtools" stepper drops paired reads that are not
        // properly paired (ignore_orphans=True).
        let reads = [pair_read("p", 0, 5, b'A', 40, false), pair_read("p", 5, 0, b'A', 40, false)];
        let (cov, visited) = genebody_coverage_with_visited(&reads, &[8], 0.0);
        assert_eq!(cov, vec![0.0]);
        assert_eq!(visited, vec![false]);
    }

    #[test]
    fn render_tin_xls_and_summary_exact_text() {
        // Verified against python3: str(85.5), str(0.0), str(100).
        let records = vec![
            TinRecord { gene_name: "tx1".into(), chrom: "chr1".into(), tx_start: 0, tx_end: 100, score: 85.5 },
            TinRecord { gene_name: "tx2".into(), chrom: "chr1".into(), tx_start: 200, tx_end: 300, score: 0.0 },
        ];
        let xls = render_tin_xls(&records);
        assert_eq!(xls, "geneID\tchrom\ttx_start\ttx_end\tTIN\ntx1\tchr1\t0\t100\t85.5\ntx2\tchr1\t200\t300\t0.0\n");

        let summary = TinSummary { mean: 42.75, median: 42.75, stdev: 0.0 };
        let summary_text = render_summary("sample.bam", &summary);
        assert_eq!(summary_text, "Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)\nsample.bam\t42.75\t42.75\t0.0\n");
    }

    /// Builds a coordinate-sorted record stream spanning several transcripts,
    /// so a windowed run has to advance, retire and switch reference sequence.
    fn multi_transcript_workload() -> (sam::Header, Vec<bam::Record>) {
        use sam::alignment::record::cigar::op::{Kind, Op};
        use sam::alignment::record_buf::{Cigar, RecordBuf};

        let header = sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                sam::header::record::value::Map::<sam::header::record::value::map::ReferenceSequence>::new(
                    std::num::NonZeroUsize::new(100_000).unwrap(),
                ),
            )
            .add_reference_sequence(
                "chr2",
                sam::header::record::value::Map::<sam::header::record::value::map::ReferenceSequence>::new(
                    std::num::NonZeroUsize::new(100_000).unwrap(),
                ),
            )
            .build();

        let mut bufs: Vec<RecordBuf> = Vec::new();
        let push = |bufs: &mut Vec<RecordBuf>, ref_id: usize, start: usize, n: usize| {
            for i in 0..n {
                let p = start + (i * 3);
                bufs.push(
                    RecordBuf::builder()
                        .set_name(format!("{ref_id}_{i}"))
                        .set_flags(sam::alignment::record::Flags::empty())
                        .set_reference_sequence_id(ref_id)
                        .set_alignment_start(noodles_core::Position::try_from(p + 1).unwrap())
                        .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
                        .build(),
                );
            }
        };
        // chr1: three dense clusters at 0-2000, 20000-22000, 50000-52000
        push(&mut bufs, 0, 0, 60);
        push(&mut bufs, 0, 20_000, 60);
        push(&mut bufs, 0, 50_000, 60);
        // chr2: one cluster
        push(&mut bufs, 1, 1_000, 40);

        (header.clone(), to_bam_records(&header, &bufs))
    }

    fn synthetic_samples() -> Vec<TranscriptSample> {
        let mk = |name: &str, chrom: &str, start: i64, end: i64| TranscriptSample {
            gene_name: name.into(),
            chrom: chrom.into(),
            tx_start: start,
            tx_end: end,
            intron_size: 0,
            // Deliberately NOT in coordinate order, to pin that the windowed
            // driver restores input order in its output.
            chosen_bases: (start + 1..=end).step_by(10).collect(),
        };
        vec![
            mk("c1_late", "chr1", 50_000, 52_000),
            mk("c1_early", "chr1", 0, 2_000),
            mk("c2", "chr2", 1_000, 3_000),
            mk("c1_mid", "chr1", 20_000, 22_000),
        ]
    }

    #[test]
    fn windowed_matches_whole_file_score_for_score() {
        let (header, records) = multi_transcript_workload();
        let samples = synthetic_samples();

        let index = build_read_index(records.iter().cloned().map(Ok), &header).unwrap();
        let (want_records, want_summary) = compute_tin(&samples, &index, 5, None);

        let got = compute_tin_windowed(records.iter().cloned().map(Ok), &header, &samples, 5, None).unwrap();
        let WindowedTin::Computed(got_records, got_summary) = got else { panic!("expected Computed") };

        // Row order must be the BED's input order, not coordinate order.
        assert_eq!(got_records.len(), want_records.len());
        for (g, w) in got_records.iter().zip(&want_records) {
            assert_eq!(g.gene_name, w.gene_name);
            assert_eq!(g.chrom, w.chrom);
            assert_eq!((g.tx_start, g.tx_end), (w.tx_start, w.tx_end));
            assert_eq!(g.score.to_bits(), w.score.to_bits(), "TIN differs for {}", g.gene_name);
        }
        // Bit-exact summary, including pairwise-summation order.
        assert_eq!(got_summary.mean.to_bits(), want_summary.mean.to_bits());
        assert_eq!(got_summary.median.to_bits(), want_summary.median.to_bits());
        assert_eq!(got_summary.stdev.to_bits(), want_summary.stdev.to_bits());
    }

    #[test]
    fn windowed_detects_an_out_of_order_start() {
        let (header, records) = multi_transcript_workload();
        // Swap two adjacent chr1 reads so start order breaks.
        let mut shuffled = records.clone();
        shuffled.swap(5, 9);
        let samples = synthetic_samples();

        let got = compute_tin_windowed(shuffled.clone().into_iter().map(Ok), &header, &samples, 5, None).unwrap();
        assert!(matches!(got, WindowedTin::NotCoordinateSorted), "expected the unsorted input to be rejected");

        // Sanity: the whole-file path tolerates the same input, which is
        // exactly why it is the fallback.
        let index = build_read_index(shuffled.into_iter().map(Ok), &header).unwrap();
        let _ = compute_tin(&samples, &index, 5, None);
    }

    #[test]
    fn windowed_scores_zero_for_a_chromosome_absent_from_the_header() {
        let (header, records) = multi_transcript_workload();
        let mut samples = synthetic_samples();
        samples.push(TranscriptSample {
            gene_name: "ghost".into(),
            chrom: "chrZZ".into(),
            tx_start: 0,
            tx_end: 1_000,
            intron_size: 0,
            chosen_bases: vec![1, 500, 1000],
        });

        let index = build_read_index(records.iter().cloned().map(Ok), &header).unwrap();
        let (_, want) = compute_tin(&samples, &index, 5, None);

        let got = compute_tin_windowed(records.iter().cloned().map(Ok), &header, &samples, 5, None).unwrap();
        let WindowedTin::Computed(got_records, got_summary) = got else { panic!("expected Computed") };
        assert_eq!(got_records.last().unwrap().gene_name, "ghost");
        assert_eq!(got_records.last().unwrap().score, 0.0);
        assert_eq!(got_summary.mean.to_bits(), want.mean.to_bits());
    }
}
