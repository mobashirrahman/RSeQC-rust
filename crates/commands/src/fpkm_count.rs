//! Port of `FPKM_count.py`: calculate fragment counts, FPM, and FPKM for
//! BED12 transcript models. Algorithm ported directly from
//! `oracle/upstream-src/scripts/FPKM_count.py`.
//!
//! No BAI index support (as with earlier commands): `count_transcript`'s
//! `samfile.fetch(chrom, tx_start, tx_end)` calls are served from an
//! in-memory, per-(raw, case-preserved)-chromosome, start-sorted read
//! index built once per BAM, not real htslib random access.
//!
//! **Preserves several deliberate upstream quirks, do not "fix"**:
//! - Fragment overlap formulas (`read_end`/`frag_end`, and a mapped mate's
//!   estimated span) use `start + rlen`, where `rlen` is the READ'S OWN
//!   query length (SEQ length), matching the upstream counting formulas.
//!   Region-fetch selection is separate and uses the true CIGAR reference
//!   span, as htslib does; a spliced/deletion-containing alignment can
//!   therefore be fetched even when its query-length estimate ends before
//!   the transcript.
//! - `count_total_fragments` and `count_transcript` disagree on case
//!   folding: the former uppercases chromosome names when testing
//!   membership in the global exon index (built with uppercased keys),
//!   but the latter's `fetch()`-equivalent lookup uses the RAW
//!   (un-uppercased) chromosome string from the BED line, matching
//!   upstream's real case-sensitive-fetch requirement.
//! - `count_transcript`'s paired-read branch tests exon overlap at TWO
//!   single-base points (the read's own start, and its mate's start via
//!   `pnext`), not the read's full span; if NEITHER point overlaps an
//!   exon, the read contributes nothing, even if the actual fragment
//!   spans exonic sequence elsewhere.
//! - `single_read` (the `-s` weight) only applies when exactly one end
//!   of a pair is mapped; a pair with both ends mapped always counts as
//!   a full `1`, regardless of `-s`.
//! - The BED12 gene model is read with TWO DIFFERENT robustness tiers in
//!   the SAME command: `build_global_exon_ranges` (first pass, feeds
//!   `count_total_fragments`) catches any per-line parse failure and
//!   skips it with a stderr note; the transcript-output pass
//!   (`parse_transcript_line`, called per-line from `compute_fpkm_rows`)
//!   has NO such handling and propagates a hard error that aborts the
//!   whole command, matching upstream's uncaught exception there.
//! - In stranded mode (`--strand` given), a transcript whose BED strand
//!   field is anything other than a literal `"+"`/`"-"` produces NO
//!   output row at all (silently dropped), even though it was still
//!   fully processed.

use std::collections::HashMap;
use std::io::{self, BufRead};

use noodles_bam as bam;
use noodles_sam as sam;
use rseqc_formats::interval::{Bed3, MergedRegions};

use crate::python_fmt::python_str_float;

/// Parses the legacy RSeQC strand-rule syntax (`"1++,1--,2+-,2-+"` or
/// `"++,--"`-style). Returns an empty map for `None` (unstranded).
pub fn parse_strand_rule(rule: Option<&str>) -> Result<HashMap<String, char>, String> {
    let Some(rule) = rule else { return Ok(HashMap::new()) };
    let tokens: Vec<&str> = rule.split(',').collect();
    let mut map = HashMap::new();

    if tokens.len() == 4 {
        for token in &tokens {
            let chars: Vec<char> = token.chars().collect();
            if chars.len() < 3 {
                return Err(format!("invalid paired-end strand token: {token:?}"));
            }
            map.insert(format!("{}{}", chars[0], chars[1]), chars[2]);
        }
        return Ok(map);
    }

    if tokens.len() == 2 {
        for token in &tokens {
            let chars: Vec<char> = token.chars().collect();
            if chars.len() < 2 {
                return Err(format!("invalid single-end strand token: {token:?}"));
            }
            map.insert(chars[0].to_string(), chars[1]);
        }
        return Ok(map);
    }

    Err(format!("unknown strand rule: {rule}"))
}

/// One BAM record's relevant fields for the region-fetch queries below.
/// Pre-filtered at build time on the filters `alignment_passes_legacy_
/// filters` applies (qcfail/duplicate/secondary/mapq); `is_unmapped` is
/// NOT filtered out (upstream's own filter function doesn't check it --
/// the unmapped-vs-mapped branching happens inside each counting
/// function instead).
/// One retained read.
///
/// 20 bytes, not 40: every coordinate here came out of a BAM, whose POS,
/// endpos and PNEXT are all int32 on the wire, so `i32` is not a narrowing of
/// the source format -- it IS the source format. The five flags are bit-packed
/// rather than five separate `bool`s, which would pad the struct back out to
/// 24 bytes. This matters because the whole file is retained: at 2.27M reads
/// the struct layout is the entire resident cost, and it was the reason this
/// command measured 2.1x MORE memory than the pysam original rather than less.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IndexedRead {
    pub start: i32,
    /// True CIGAR reference end used to emulate htslib fetch selection.
    pub end: i32,
    pub rlen: i32,
    /// `pnext`, 0-based; `-1` if unset (pysam's raw sentinel), matching
    /// upstream's use of the raw field regardless of pairing state.
    pub mate_start: i32,
    flags: u8,
}

const F_PAIRED: u8 = 1 << 0;
const F_READ2: u8 = 1 << 1;
const F_REVERSE: u8 = 1 << 2;
const F_UNMAPPED: u8 = 1 << 3;
const F_MATE_UNMAPPED: u8 = 1 << 4;

impl IndexedRead {
    pub fn is_paired(&self) -> bool { self.flags & F_PAIRED != 0 }
    pub fn is_read2(&self) -> bool { self.flags & F_READ2 != 0 }
    pub fn is_reverse(&self) -> bool { self.flags & F_REVERSE != 0 }
    pub fn is_unmapped(&self) -> bool { self.flags & F_UNMAPPED != 0 }
    pub fn mate_is_unmapped(&self) -> bool { self.flags & F_MATE_UNMAPPED != 0 }

    /// Coordinates arrive from the BAM as `usize`; refuse rather than truncate
    /// if a value cannot be represented, so a malformed record is loud instead
    /// of silently reporting a wrong position. A valid BAM cannot trip this --
    /// POS is int32 on the wire -- which is why the conversion is exact.
    fn coord(v: usize) -> io::Result<i32> {
        i32::try_from(v).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("coordinate {v} does not fit in i32; not a valid BAM position"),
            )
        })
    }

    fn signed_coord(v: i64) -> io::Result<i32> {
        i32::try_from(v).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("coordinate {v} does not fit in i32; not a valid BAM position"),
            )
        })
    }
}

fn alignment_passes_legacy_filters(flags: sam::alignment::record::Flags, mapq: u8, skip_multi: bool, map_qual: u8) -> bool {
    if flags.is_qc_fail() || flags.is_duplicate() || flags.is_secondary() {
        return false;
    }
    if skip_multi && mapq < map_qual {
        return false;
    }
    true
}

/// Converts one BAM record into an [`IndexedRead`], applying the legacy
/// filters. `Ok(None)` means the record is filtered out.
///
/// Shared by `build_read_index` (whole-file path) and the two streaming
/// paths, so the filter set cannot drift between them -- and so the streaming
/// totals cannot disagree with the whole-file totals.
fn to_indexed_read(record: &bam::Record, skip_multi: bool, map_qual: u8) -> io::Result<Option<IndexedRead>> {
    let flags = record.flags();
    let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
    if !alignment_passes_legacy_filters(flags, mapq, skip_multi, map_qual) {
        return Ok(None);
    }

    let rlen = record.sequence().len() as i32;

    // For unmapped reads, use pos directly from the record (the coordinate field
    // still contains meaningful data for unpaired unmapped reads or when the mate
    // is mapped). For mapped reads, use alignment_start from CIGAR. This matches
    // Python's pysam behavior which exposes pos for unmapped reads too.
    let (start, reference_end) = if flags.is_unmapped() {
        // Unmapped reads: use the raw POS field as start, estimate end as start+rlen
        let pos = record.alignment_start().transpose()?.map_or(0, |p| p.get() - 1);
        (pos, (pos as i64 + rlen as i64) as usize)
    } else {
        // Mapped reads: use CIGAR to compute reference span
        let Some(pos) = record.alignment_start().transpose()? else { return Ok(None) };
        let start = pos.get() - 1;
        let cigar_ops: Vec<_> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;
        let (_, ref_end) = rseqc_formats::cigar::reference_span(start, cigar_ops);
        (start, ref_end)
    };

    let mate_start = record.mate_alignment_start().transpose()?.map_or(-1i64, |p| (p.get() - 1) as i64);

    let mut packed = 0u8;
    if flags.is_segmented() {
        packed |= F_PAIRED;
    }
    if flags.is_last_segment() {
        packed |= F_READ2;
    }
    if flags.is_reverse_complemented() {
        packed |= F_REVERSE;
    }
    if flags.is_unmapped() {
        packed |= F_UNMAPPED;
    }
    if flags.is_mate_unmapped() {
        packed |= F_MATE_UNMAPPED;
    }

    Ok(Some(IndexedRead {
        start: IndexedRead::coord(start)?,
        end: IndexedRead::coord(reference_end)?,
        rlen,
        mate_start: IndexedRead::signed_coord(mate_start)?,
        flags: packed,
    }))
}

/// Builds a per-(raw)-chromosome, start-sorted read index from a BAM.
pub fn build_read_index<I>(
    records: I,
    header: &sam::Header,
    skip_multi: bool,
    map_qual: u8,
) -> io::Result<HashMap<String, Vec<IndexedRead>>>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut by_chrom: HashMap<String, Vec<IndexedRead>> = HashMap::new();

    for result in records {
        let record = result?;
        // Order matters for machinery, not arithmetic: the reference id and
        // chromosome are resolved BEFORE the record is converted, so a record
        // with no reference sequence is skipped rather than failing on its
        // CIGAR, which is the order the whole-file path has always used.
        let Some(ref_id) = record.reference_sequence_id().transpose()? else { continue };
        let Some((chrom_bstr, _)) = header.reference_sequences().get_index(ref_id) else { continue };
        let Some(read) = to_indexed_read(&record, skip_multi, map_qual)? else { continue };
        by_chrom.entry(chrom_bstr.to_string()).or_default().push(read);
    }

    for reads in by_chrom.values_mut() {
        reads.sort_by_key(|r| r.start);
    }

    Ok(by_chrom)
}

/// Builds the global, uppercased-chromosome exon index used by
/// `count_total_fragments`. Robustness tier: comment/track/browser lines
/// skipped; ANY other per-line parse failure (short line, bad int) is
/// caught and skipped with a stderr note -- upstream wraps the whole
/// per-line block in one bare `except (IndexError, ValueError)`.
pub fn build_global_exon_ranges(reader: impl BufRead) -> io::Result<MergedRegions> {
    let mut exons: Vec<Bed3> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }

        let parsed: Result<(String, Vec<i64>, Vec<i64>), ()> = (|| {
            let fields: Vec<&str> = line.split_whitespace().collect();
            let chrom = fields.first().ok_or(())?.to_uppercase();
            let tx_start: i64 = fields.get(1).ok_or(())?.parse().map_err(|_| ())?;

            let exon_starts: Vec<i64> = fields
                .get(11)
                .ok_or(())?
                .trim_end_matches(',')
                .split(',')
                .map(|s| s.parse::<i64>().map(|v| v + tx_start).map_err(|_| ()))
                .collect::<Result<_, _>>()?;
            let exon_sizes: Vec<i64> = fields
                .get(10)
                .ok_or(())?
                .trim_end_matches(',')
                .split(',')
                .map(|s| s.parse::<i64>().map_err(|_| ()))
                .collect::<Result<_, _>>()?;
            let exon_ends: Vec<i64> = exon_starts.iter().zip(exon_sizes.iter()).map(|(&s, &sz)| s + sz).collect();

            Ok((chrom, exon_starts, exon_ends))
        })();

        match parsed {
            Ok((chrom, starts, ends)) => {
                for (s, e) in starts.into_iter().zip(ends) {
                    exons.push((chrom.clone(), s, e));
                }
            }
            Err(()) => {
                eprintln!("[NOTE: input BED must be 12-column] skipped this line: {line}");
            }
        }
    }

    Ok(MergedRegions::new(&exons))
}

/// One retained read's contribution to the two fragment totals.
///
/// Factored out of [`count_total_fragments`] so the whole-file and streaming
/// totals are literally the same arithmetic and cannot drift apart.
fn add_total_fragment(
    chrom_upper: &str,
    r: &IndexedRead,
    global_exon_ranges: &MergedRegions,
    single_read: f64,
    total_frags: &mut f64,
    exonic_frags: &mut f64,
) {
    let read_end = r.start as i64 + r.rlen as i64;

    if !r.is_paired() {
        *total_frags += 1.0;
        if global_exon_ranges.overlap_length(chrom_upper, r.start as i64, read_end) > 0 {
            *exonic_frags += 1.0;
        }
        return;
    }

    if r.is_read2() {
        return;
    }

    let mate_end = r.mate_start as i64 + r.rlen as i64;

    if r.is_unmapped() {
        if r.mate_is_unmapped() {
            return;
        }
        *total_frags += single_read;
        if global_exon_ranges.overlap_length(chrom_upper, r.mate_start as i64, mate_end) > 0 {
            *exonic_frags += single_read;
        }
    } else if r.mate_is_unmapped() {
        *total_frags += single_read;
        if global_exon_ranges.overlap_length(chrom_upper, r.start as i64, read_end) > 0 {
            *exonic_frags += single_read;
        }
    } else {
        *total_frags += 1.0;
        // Upstream treats a paired fragment as exonic only when
        // BOTH mapped ends overlap an exon.  Checking only read 1
        // incorrectly classifies a pair whose mate lies outside the
        // gene model and inflates the exonic denominator used by
        // `-e` FPKM normalization.
        if global_exon_ranges.overlap_length(chrom_upper, r.start as i64, read_end) > 0
            && global_exon_ranges.overlap_length(chrom_upper, r.mate_start as i64, mate_end) > 0
        {
            *exonic_frags += 1.0;
        }
    }
}

/// Counts total and exonic fragments across the whole BAM. Ports
/// `count_total_fragments` exactly (iteration order over the read index
/// doesn't affect these running sums).
pub fn count_total_fragments(
    reads_by_chrom: &HashMap<String, Vec<IndexedRead>>,
    global_exon_ranges: &MergedRegions,
    single_read: f64,
) -> (f64, f64) {
    let mut total_frags = 0.0;
    let mut exonic_frags = 0.0;

    for (chrom, reads) in reads_by_chrom {
        let chrom_upper = chrom.to_uppercase();
        for r in reads {
            add_total_fragment(
                &chrom_upper,
                r,
                global_exon_ranges,
                single_read,
                &mut total_frags,
                &mut exonic_frags,
            );
        }
    }

    (total_frags, exonic_frags)
}

/// Result of the streaming totals pass.
pub struct StreamingTotals {
    pub total_frags: f64,
    pub exonic_frags: f64,
    /// False when the stream was not in coordinate order, in which case the
    /// sliding-window driver cannot be used and the caller must fall back to the
    /// whole-file read index.
    ///
    /// This is checked here because the totals pass already visits every record,
    /// so the check is free -- and because the windowed driver's own
    /// `start < last_start` probe is NOT sufficient. That probe only fires for a
    /// record the pull loop actually reaches; an out-of-order record beyond the
    /// current transcript's `tx_end` makes the loop `break` before comparing, so
    /// a shuffled BAM could silently score zero for transcripts whose window
    /// never opened. Measured: on a shuffled copy of the 2.05M-read alignment,
    /// the driver returned `Computed` with 1729 of 3000 rows wrong (all zero)
    /// instead of reporting `NotCoordinateSorted`.
    pub coordinate_sorted: bool,
}

/// [`count_total_fragments`] without materialising the file.
///
/// Same arithmetic on the same retained reads, in BAM order rather than
/// HashMap order. That is not a behaviour change: the accumulated terms are
/// `0`, `0.5` and `1` (the `-s` values upstream accepts), all exact in binary
/// floating point, so the sums are associative and were already independent of
/// iteration order -- which the whole-file path relies on, since it iterates a
/// `HashMap`. Streaming makes the order deterministic instead of leaving it to
/// the hasher.
pub fn count_total_fragments_streaming<I>(
    records: I,
    header: &sam::Header,
    global_exon_ranges: &MergedRegions,
    single_read: f64,
    skip_multi: bool,
    map_qual: u8,
) -> io::Result<StreamingTotals>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut total_frags = 0.0;
    let mut exonic_frags = 0.0;
    // Ref names repeat for the whole of each reference sequence, so uppercasing
    // is cached per reference id rather than per read.
    let mut upper_cache: HashMap<usize, String> = HashMap::new();

    let mut coordinate_sorted = true;
    let mut last_key: (usize, i64) = (0, i64::MIN);

    for result in records {
        let record = result?;

        let ref_id = record.reference_sequence_id().transpose()?;
        let start = record.alignment_start().transpose()?.map(|p| (p.get() - 1) as i64);
        // A record with no reference id (tid -1) or no position is the trailing
        // unmapped block of a coordinate-sorted BAM. Key it past every real
        // record, so anything appearing after it counts as out of order. A
        // false positive here is harmless -- it only means the whole-file path
        // runs -- so the check is deliberately conservative.
        let key = match (ref_id, start) {
            (Some(r), Some(s)) => (r, s),
            _ => (usize::MAX, i64::MAX),
        };
        if key < last_key {
            coordinate_sorted = false;
        }
        last_key = key;

        let Some(ref_id) = ref_id else { continue };
        let Some((chrom_bstr, _)) = header.reference_sequences().get_index(ref_id) else { continue };
        let Some(read) = to_indexed_read(&record, skip_multi, map_qual)? else { continue };
        let chrom_upper = upper_cache
            .entry(ref_id)
            .or_insert_with(|| std::str::from_utf8(chrom_bstr).unwrap_or("").to_uppercase());
        add_total_fragment(
            chrom_upper,
            &read,
            global_exon_ranges,
            single_read,
            &mut total_frags,
            &mut exonic_frags,
        );
    }

    Ok(StreamingTotals { total_frags, exonic_frags, coordinate_sorted })
}

/// Counts one transcript's forward/reverse/unstranded fragments. Ports
/// `count_transcript` exactly, including its two-separate-`if` (not
/// `if`/`elif`) structure for the unpaired/paired branches.
#[allow(clippy::too_many_arguments)]
pub fn count_transcript(
    reads: &[IndexedRead],
    chrom: &str,
    tx_start: i64,
    tx_end: i64,
    exon_ranges: &MergedRegions,
    strand_rule_active: bool,
    strand_map: &HashMap<String, char>,
    single_read: f64,
) -> (f64, f64, f64) {
    let mut frag_count_f = 0.0;
    let mut frag_count_r = 0.0;
    let mut frag_count_fr = 0.0;

    // `reads` is the whole chromosome, start-sorted (whole-file path) or this
    // transcript's window, also start-sorted (windowed path). Every entry in a
    // window already satisfies `start < tx_end`, so `hi` is the window length
    // there; in both cases the predicate below is exactly `fetch()`'s.
    let hi = reads.partition_point(|r| (r.start as i64) < tx_end);

    for r in reads[..hi].iter().filter(|r| (r.end as i64) > tx_start) {
        if !r.is_paired() {
            let frag_st = r.start as i64;
            let frag_end = r.start as i64 + r.rlen as i64;
            let strand_key = if r.is_reverse() { "-" } else { "+" };

            if exon_ranges.overlap_length(chrom, frag_st, frag_end) > 0 {
                if !strand_rule_active {
                    frag_count_fr += 1.0;
                } else if strand_map.get(strand_key) == Some(&'+') {
                    frag_count_f += 1.0;
                } else if strand_map.get(strand_key) == Some(&'-') {
                    frag_count_r += 1.0;
                }
            }
        }

        if r.is_paired() {
            let frag_st = r.start as i64;
            let frag_end = r.mate_start as i64;

            let start_in_exon = exon_ranges.overlap_length(chrom, frag_st, frag_st + 1) > 0;
            let mate_in_exon = exon_ranges.overlap_length(chrom, frag_end, frag_end + 1) > 0;
            if !start_in_exon && !mate_in_exon {
                continue;
            }

            if r.is_read2() {
                continue;
            }

            let strand_key = if r.is_reverse() { "1-" } else { "1+" };

            if !strand_rule_active {
                if r.is_unmapped() {
                    if r.mate_is_unmapped() {
                        continue;
                    }
                    frag_count_fr += single_read;
                } else if r.mate_is_unmapped() {
                    frag_count_fr += single_read;
                } else {
                    frag_count_fr += 1.0;
                }
            } else {
                if strand_map.get(strand_key) == Some(&'+') {
                    if r.is_unmapped() {
                        if !r.mate_is_unmapped() {
                            frag_count_f += single_read;
                        }
                    } else if r.mate_is_unmapped() {
                        frag_count_f += single_read;
                    } else {
                        frag_count_f += 1.0;
                    }
                }
                if strand_map.get(strand_key) == Some(&'-') {
                    if r.is_unmapped() {
                        if !r.mate_is_unmapped() {
                            frag_count_r += single_read;
                        }
                    } else if r.mate_is_unmapped() {
                        frag_count_r += single_read;
                    } else {
                        frag_count_r += 1.0;
                    }
                }
            }
        }
    }

    (frag_count_f, frag_count_r, frag_count_fr)
}

#[derive(Debug, Clone, PartialEq)]
pub struct TranscriptLine {
    pub chrom: String,
    pub tx_start: i64,
    pub tx_end: i64,
    pub gene_name: String,
    pub gene_strand: String,
    pub exon_starts: Vec<i64>,
    pub exon_ends: Vec<i64>,
}

/// Ports `parse_transcript_line`. Upstream calls this with NO
/// surrounding try/except in `main()`'s output-writing loop (unlike
/// `build_global_exon_ranges`'s tolerant first pass), so a malformed
/// line here is a hard, whole-command-aborting error -- matched by
/// returning `io::Result` and propagating `?` at the call site.
pub fn parse_transcript_line(line: &str) -> io::Result<TranscriptLine> {
    let err = |msg: &str| io::Error::new(io::ErrorKind::InvalidData, msg.to_string());
    let fields: Vec<&str> = line.split_whitespace().collect();

    let chrom = fields.first().ok_or_else(|| err("missing chrom field"))?.to_string();
    let tx_start: i64 =
        fields.get(1).ok_or_else(|| err("missing txStart"))?.parse().map_err(|_| err("invalid txStart"))?;
    let tx_end: i64 = fields.get(2).ok_or_else(|| err("missing txEnd"))?.parse().map_err(|_| err("invalid txEnd"))?;
    let gene_name = fields.get(3).ok_or_else(|| err("missing name"))?.to_string();
    let gene_strand = fields.get(5).ok_or_else(|| err("missing strand"))?.replace(' ', "_");

    let exon_starts: Vec<i64> = fields
        .get(11)
        .ok_or_else(|| err("missing blockStarts"))?
        .trim_end_matches(',')
        .split(',')
        .map(|s| s.parse::<i64>().map(|v| v + tx_start).map_err(|_| err("invalid blockStart")))
        .collect::<io::Result<_>>()?;
    let exon_sizes: Vec<i64> = fields
        .get(10)
        .ok_or_else(|| err("missing blockSizes"))?
        .trim_end_matches(',')
        .split(',')
        .map(|s| s.parse::<i64>().map_err(|_| err("invalid blockSize")))
        .collect::<io::Result<_>>()?;
    let exon_ends: Vec<i64> = exon_starts.iter().zip(exon_sizes.iter()).map(|(&s, &sz)| s + sz).collect();

    Ok(TranscriptLine { chrom, tx_start, tx_end, gene_name, gene_strand, exon_starts, exon_ends })
}

#[derive(Debug, Clone, PartialEq)]
pub struct FpkmRow {
    pub chrom: String,
    pub tx_start: i64,
    pub tx_end: i64,
    pub gene_name: String,
    pub mrna_size: f64,
    pub gene_strand: String,
    pub frag_count: f64,
    pub fpm: f64,
    pub fpkm: f64,
}

/// Computes one output row per transcript line (skipping rows silently
/// dropped by the strand-rule tri-branch -- see module docs). Ports the
/// per-transcript loop in `main()`, minus its progress/report printing.
pub fn compute_fpkm_rows(
    reader: impl BufRead,
    reads_by_chrom: &HashMap<String, Vec<IndexedRead>>,
    strand_rule_active: bool,
    strand_map: &HashMap<String, char>,
    single_read: f64,
    denominator: f64,
    mut on_transcript_finished: impl FnMut(usize),
) -> io::Result<Vec<FpkmRow>> {
    let mut rows = Vec::new();
    let mut gene_finished = 0usize;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }

        let t = parse_transcript_line(&line)?;

        let mut mrna_size = 0.0;
        let mut exon_pairs: Vec<Bed3> = Vec::with_capacity(t.exon_starts.len());
        for (&s, &e) in t.exon_starts.iter().zip(t.exon_ends.iter()) {
            mrna_size += (e - s) as f64;
            exon_pairs.push((t.chrom.clone(), s, e));
        }
        let exon_ranges = MergedRegions::new(&exon_pairs);

        let reads: &[IndexedRead] =
            reads_by_chrom.get(&t.chrom).map(Vec::as_slice).unwrap_or(&[]);
        let (frag_f, frag_r, frag_fr) = count_transcript(
            reads,
            &t.chrom,
            t.tx_start,
            t.tx_end,
            &exon_ranges,
            strand_rule_active,
            strand_map,
            single_read,
        );

        let fpm_fr = frag_fr * 1_000_000.0 / denominator;
        let fpm_f = frag_f * 1_000_000.0 / denominator;
        let fpm_r = frag_r * 1_000_000.0 / denominator;
        let fpkm_fr = frag_fr * 1_000_000_000.0 / (denominator * mrna_size);
        let fpkm_f = frag_f * 1_000_000_000.0 / (denominator * mrna_size);
        let fpkm_r = frag_r * 1_000_000_000.0 / (denominator * mrna_size);

        let selected = if !strand_rule_active {
            Some((frag_fr, fpm_fr, fpkm_fr))
        } else if t.gene_strand == "+" {
            Some((frag_f, fpm_f, fpkm_f))
        } else if t.gene_strand == "-" {
            Some((frag_r, fpm_r, fpkm_r))
        } else {
            None
        };

        if let Some((frag_count, fpm, fpkm)) = selected {
            rows.push(FpkmRow {
                chrom: t.chrom,
                tx_start: t.tx_start,
                tx_end: t.tx_end,
                gene_name: t.gene_name,
                mrna_size,
                gene_strand: t.gene_strand,
                frag_count,
                fpm,
                fpkm,
            });
        }
        // Upstream's `gene_finished` counter advances for every parsed BED
        // line, including strand-ruled lines whose strand is neither + nor -.
        gene_finished += 1;
        on_transcript_finished(gene_finished);
    }

    Ok(rows)
}

/// The windowed driver's result, mirroring `genebody_coverage`'s
/// `WindowedCoverage`.
pub enum WindowedFpkm {
    Computed(Vec<FpkmRow>),
    /// The input was not coordinate-sorted, so a sliding window cannot be
    /// maintained. The caller should fall back to the whole-file `build_read_index`.
    NotCoordinateSorted,
}

/// One parsed, ready-to-count transcript.
struct ParsedTranscript {
    t: TranscriptLine,
    mrna_size: f64,
    exon_ranges: MergedRegions,
}

/// Computes the same rows as [`compute_fpkm_rows`] while holding only a
/// **sliding window** of reads instead of the whole BAM.
///
/// Upstream is the mirror image of this port, and the reason for the shape
/// here. `count_transcript` in `FPKM_count.py` calls
/// `samfile.fetch(chrom, tx_start, tx_end)` once per transcript, so pysam holds
/// one region's reads at a time (low memory) and re-seeks per transcript (slow).
/// This port instead loads every read once (fast) and scans a start-sorted
/// prefix per transcript (memory-heavy, and O(reads) per transcript). Measured
/// at 2.1x upstream's memory on a 2.27M-read alignment and 3.1x on the
/// 8.2M-record rat alignment, this is the same trade `compute_tin_windowed` and
/// `compute_coverage_windowed` already remove for their commands.
///
/// Two facts make the window exact rather than approximate:
///
/// 1. Transcripts are visited in coordinate order, so `tx_start` is
///    non-decreasing. A read with `end <= tx_start` therefore cannot overlap
///    this transcript or any later one, and is discarded at PUSH time rather
///    than after buffering -- the mistake that cost `geneBody_coverage` its
///    window, where an inter-transcript gap had to be held in full before the
///    trim ran.
/// 2. `count_transcript` selects `start < tx_end && end > tx_start`. Everything
///    the pull loop stops at has `start >= window_end`, so every record with
///    `start < tx_end` has already been offered to the window; and everything
///    discarded has `end <= tx_start`. The window therefore contains exactly the
///    reads `count_transcript` would select from the whole chromosome, in the
///    same (file) order.
///
/// Rows are collected by original BED index and flattened at the end, so the
/// `.FPKM.xls` row order is unchanged.
#[allow(clippy::too_many_arguments)]
pub fn compute_fpkm_rows_windowed<I>(
    records: I,
    header: &sam::Header,
    bed: impl BufRead,
    strand_rule_active: bool,
    strand_map: &HashMap<String, char>,
    single_read: f64,
    denominator: f64,
    skip_multi: bool,
    map_qual: u8,
    mut on_transcript_finished: impl FnMut(usize),
) -> io::Result<WindowedFpkm>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    // The model is parsed up front because the BAM has to be walked in
    // coordinate order. A malformed line is still the hard, whole-command
    // error the whole-file path raises, at the same first offending line.
    let mut parsed: Vec<ParsedTranscript> = Vec::new();
    for line in bed.lines() {
        let line = line?;
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }
        let t = parse_transcript_line(&line)?;
        let mut mrna_size = 0.0;
        let mut exon_pairs: Vec<Bed3> = Vec::with_capacity(t.exon_starts.len());
        for (&s, &e) in t.exon_starts.iter().zip(t.exon_ends.iter()) {
            mrna_size += (e - s) as f64;
            exon_pairs.push((t.chrom.clone(), s, e));
        }
        parsed.push(ParsedTranscript { exon_ranges: MergedRegions::new(&exon_pairs), mrna_size, t });
    }

    let ref_ids: HashMap<&str, usize> = header
        .reference_sequences()
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (std::str::from_utf8(name).unwrap_or(""), i))
        .collect();
    let ref_index = |chrom: &str| -> usize { ref_ids.get(chrom).copied().unwrap_or(usize::MAX) };

    // Reference-sequence order, then transcript start, then input order. A
    // transcript on a chromosome the BAM does not contain sorts last and
    // produces no row, which is what the whole-file path's failed lookup does.
    let mut order: Vec<usize> = (0..parsed.len()).collect();
    order.sort_by_key(|&i| (ref_index(&parsed[i].t.chrom), parsed[i].t.tx_start, i));

    let mut rows: Vec<Option<FpkmRow>> = (0..parsed.len()).map(|_| None).collect();
    let mut window: Vec<IndexedRead> = Vec::new();
    let mut window_ref: Option<usize> = None;
    let mut closed: Vec<usize> = Vec::new();
    let mut last_start: i64 = i64::MIN;
    let mut records = records.into_iter().peekable();
    let mut gene_finished = 0usize;

    for &ti in &order {
        let pt = &parsed[ti];
        let t_ref = ref_index(&pt.t.chrom);

        let counts = if t_ref == usize::MAX {
            // Chromosome absent from the BAM. The whole-file path looks it up,
            // finds nothing, and still emits a zero row -- so this must too. The
            // window and the stream are deliberately untouched: an absent
            // chromosome sorts last, after every real one is finished.
            count_transcript(
                &[],
                &pt.t.chrom,
                pt.t.tx_start,
                pt.t.tx_end,
                &pt.exon_ranges,
                strand_rule_active,
                strand_map,
                single_read,
            )
        } else {
            let window_start = pt.t.tx_start;
            let window_end = pt.t.tx_end;

            if window_ref != Some(t_ref) {
                if let Some(r) = window_ref {
                    closed.push(r);
                }
                window = Vec::new();
                window_ref = Some(t_ref);
                last_start = i64::MIN;
            }

            // Drop what the previous transcript left behind before pulling.
            window.retain(|r| (r.end as i64) > window_start);

            while let Some(peeked) = records.peek() {
                let next = match peeked {
                    Ok(record) => record,
                    Err(_) => break,
                };
                // No reference id (tid -1): the trailing unmapped block in a
                // coordinate-sorted BAM. Consume rather than stop, or every
                // later transcript starves.
                let Some(next_ref) = next.reference_sequence_id().transpose()? else {
                    records.next();
                    continue;
                };
                if header.reference_sequences().get_index(next_ref).is_none() {
                    break;
                }
                if next_ref > t_ref {
                    break;
                }
                if next_ref < t_ref {
                    records.next();
                    continue;
                }
                if closed.contains(&next_ref) {
                    return Ok(WindowedFpkm::NotCoordinateSorted);
                }
                let Some(pos) = next.alignment_start().transpose()? else { break };
                let start = (pos.get() - 1) as i64;
                if start >= window_end {
                    break;
                }
                if start < last_start {
                    return Ok(WindowedFpkm::NotCoordinateSorted);
                }
                last_start = start;
                let read = to_indexed_read(next, skip_multi, map_qual)?;
                records.next();
                if let Some(read) = read {
                    if (read.end as i64) > window_start {
                        window.push(read);
                    }
                }
            }

            count_transcript(
                &window,
                &pt.t.chrom,
                pt.t.tx_start,
                pt.t.tx_end,
                &pt.exon_ranges,
                strand_rule_active,
                strand_map,
                single_read,
            )
        };

        let (frag_f, frag_r, frag_fr) = counts;
        let fpm_fr = frag_fr * 1_000_000.0 / denominator;
        let fpm_f = frag_f * 1_000_000.0 / denominator;
        let fpm_r = frag_r * 1_000_000.0 / denominator;
        let fpkm_fr = frag_fr * 1_000_000_000.0 / (denominator * pt.mrna_size);
        let fpkm_f = frag_f * 1_000_000_000.0 / (denominator * pt.mrna_size);
        let fpkm_r = frag_r * 1_000_000_000.0 / (denominator * pt.mrna_size);

        let selected = if !strand_rule_active {
            Some((frag_fr, fpm_fr, fpkm_fr))
        } else if pt.t.gene_strand == "+" {
            Some((frag_f, fpm_f, fpkm_f))
        } else if pt.t.gene_strand == "-" {
            Some((frag_r, fpm_r, fpkm_r))
        } else {
            None
        };

        if let Some((frag_count, fpm, fpkm)) = selected {
            rows[ti] = Some(FpkmRow {
                chrom: pt.t.chrom.clone(),
                tx_start: pt.t.tx_start,
                tx_end: pt.t.tx_end,
                gene_name: pt.t.gene_name.clone(),
                mrna_size: pt.mrna_size,
                gene_strand: pt.t.gene_strand.clone(),
                frag_count,
                fpm,
                fpkm,
            });
        }

        // Upstream's `gene_finished` counter advances for every parsed BED line,
        // including one on a chromosome the BAM lacks and one strand-ruled to a
        // value that is neither + nor -.
        gene_finished += 1;
        on_transcript_finished(gene_finished);
    }

    for result in records {
        result?;
    }

    Ok(WindowedFpkm::Computed(rows.into_iter().flatten().collect()))
}

/// Ports the `.FPKM.xls` output: a header line followed by one row per
/// (non-dropped) transcript, in input order.
pub fn render_fpkm_xls(rows: &[FpkmRow]) -> String {
    let mut out = String::from("#chrom\tst\tend\taccession\tmRNA_size\tgene_strand\tFrag_count\tFPM\tFPKM\n");
    for r in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.chrom,
            r.tx_start,
            r.tx_end,
            r.gene_name,
            python_str_float(r.mrna_size),
            r.gene_strand,
            python_str_float(r.frag_count),
            python_str_float(r.fpm),
            python_str_float(r.fpkm),
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn parse_strand_rule_four_and_two_token_forms() {
        let map4 = parse_strand_rule(Some("1++,1--,2+-,2-+")).unwrap();
        assert_eq!(map4.get("1+"), Some(&'+'));
        assert_eq!(map4.get("1-"), Some(&'-'));
        assert_eq!(map4.get("2+"), Some(&'-'));
        assert_eq!(map4.get("2-"), Some(&'+'));

        let map2 = parse_strand_rule(Some("++,--")).unwrap();
        assert_eq!(map2.get("+"), Some(&'+'));
        assert_eq!(map2.get("-"), Some(&'-'));

        assert!(parse_strand_rule(Some("bad")).is_err());
        assert_eq!(parse_strand_rule(None).unwrap(), HashMap::new());
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
    fn count_total_fragments_unpaired_exonic_and_nonexonic() {
        use sam::alignment::record::{Flags, MappingQuality, cigar::Op, cigar::op::Kind};
        use sam::alignment::record_buf::{Cigar, RecordBuf};

        let header = header_with_chrom("chr1", 1000);
        let exonic = RecordBuf::builder()
            .set_name("a")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(101).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 30]))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![40; 30]))
            .build();
        let nonexonic = RecordBuf::builder()
            .set_name("b")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(501).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 30]))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![40; 30]))
            .build();

        let bam_records = to_bam_records(&header, &[exonic, nonexonic]);
        let index = build_read_index(bam_records.into_iter().map(Ok), &header, false, 0).unwrap();

        let exon_ranges = MergedRegions::new(&[("CHR1".to_string(), 100, 200)]);
        let (total, exonic_count) = count_total_fragments(&index, &exon_ranges, 1.0);
        assert_eq!(total, 2.0);
        assert_eq!(exonic_count, 1.0);
    }

    #[test]
    fn count_total_fragments_paired_single_read_weighting() {
        use sam::alignment::record::{Flags, MappingQuality, cigar::Op, cigar::op::Kind};
        use sam::alignment::record_buf::{Cigar, RecordBuf};

        let header = header_with_chrom("chr1", 1000);

        // Read1 mapped inside the exon, mate unmapped -> single_read weight,
        // counted via this read's own span.
        let read1 = RecordBuf::builder()
            .set_name("p1")
            .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT | Flags::MATE_UNMAPPED)
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(101).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 30]))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![40; 30]))
            .build();

        let bam_records = to_bam_records(&header, &[read1]);
        let index = build_read_index(bam_records.into_iter().map(Ok), &header, false, 0).unwrap();
        let exon_ranges = MergedRegions::new(&[("CHR1".to_string(), 100, 200)]);

        let (total, exonic_count) = count_total_fragments(&index, &exon_ranges, 0.5);
        assert_eq!(total, 0.5);
        assert_eq!(exonic_count, 0.5);
    }

    #[test]
    fn count_total_fragments_requires_both_mapped_mates_to_overlap_exon() {
        use sam::alignment::record::{Flags, MappingQuality, cigar::Op, cigar::op::Kind};
        use sam::alignment::record_buf::{Cigar, RecordBuf};

        let header = header_with_chrom("chr1", 1000);
        let read1 = RecordBuf::builder()
            .set_name("mate_outside")
            .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(111).unwrap())
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(noodles_core::Position::try_from(501).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 20)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 20]))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![40; 20]))
            .build();

        let bam_records = to_bam_records(&header, &[read1]);
        let index = build_read_index(bam_records.into_iter().map(Ok), &header, false, 0).unwrap();
        let exon_ranges = MergedRegions::new(&[("CHR1".to_string(), 100, 200)]);

        let (total, exonic_count) = count_total_fragments(&index, &exon_ranges, 1.0);
        assert_eq!(total, 1.0);
        assert_eq!(exonic_count, 0.0);
    }

    #[test]
    fn count_transcript_unpaired_unstranded_and_stranded() {
        use sam::alignment::record::{Flags, MappingQuality, cigar::Op, cigar::op::Kind};
        use sam::alignment::record_buf::{Cigar, RecordBuf};

        let header = header_with_chrom("chr1", 1000);
        let fwd = RecordBuf::builder()
            .set_name("f")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(101).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 30]))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![40; 30]))
            .build();
        let rev = RecordBuf::builder()
            .set_name("r")
            .set_flags(Flags::REVERSE_COMPLEMENTED)
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(101).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 30]))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![40; 30]))
            .build();

        let bam_records = to_bam_records(&header, &[fwd, rev]);
        let index = build_read_index(bam_records.into_iter().map(Ok), &header, false, 0).unwrap();
        let exon_ranges = MergedRegions::new(&[("chr1".to_string(), 100, 200)]);

        // Unstranded: both reads land in frag_count_fr.
        let (f, r, fr) = count_transcript(&index["chr1"], "chr1", 0, 1000, &exon_ranges, false, &HashMap::new(), 1.0);
        assert_eq!((f, r, fr), (0.0, 0.0, 2.0));

        // Stranded "++,--" -> forward read maps to '+', reverse read maps to '-'.
        let map = parse_strand_rule(Some("++,--")).unwrap();
        let (f, r, fr) = count_transcript(&index["chr1"], "chr1", 0, 1000, &exon_ranges, true, &map, 1.0);
        assert_eq!((f, r, fr), (1.0, 1.0, 0.0));
    }

    #[test]
    fn parse_transcript_line_hard_errors_on_malformed_line() {
        assert!(parse_transcript_line("too short").is_err());
        let t = parse_transcript_line("chr1\t0\t100\ttx1\t0\t+\t0\t100\t0\t1\t100,\t0,").unwrap();
        assert_eq!(t.gene_name, "tx1");
        assert_eq!(t.exon_starts, vec![0]);
        assert_eq!(t.exon_ends, vec![100]);
    }

    #[test]
    fn build_global_exon_ranges_skips_malformed_lines() {
        let bed = "\
track name=x
too short
chr1\t0\t100\ttx1\t0\t+\t0\t100\t0\t1\t100,\t0,
";
        let ranges = build_global_exon_ranges(Cursor::new(bed)).unwrap();
        assert!(ranges.overlap_length("CHR1", 0, 100) > 0);
    }

    #[test]
    fn render_fpkm_xls_exact_text() {
        // Verified against python3: str(200.0), str(1.0), str(5000000.0), str(2.5e8).
        let rows = vec![FpkmRow {
            chrom: "chr1".into(),
            tx_start: 0,
            tx_end: 200,
            gene_name: "tx1".into(),
            mrna_size: 200.0,
            gene_strand: "+".into(),
            frag_count: 1.0,
            fpm: 5_000_000.0,
            fpkm: 2.5e10,
        }];
        let xls = render_fpkm_xls(&rows);
        assert_eq!(
            xls,
            "#chrom\tst\tend\taccession\tmRNA_size\tgene_strand\tFrag_count\tFPM\tFPKM\n\
chr1\t0\t200\ttx1\t200.0\t+\t1.0\t5000000.0\t25000000000.0\n"
        );
    }

    // --- windowed driver ------------------------------------------------

    const BED: &str = "\
chr1\t100\t400\ttxA\t0\t+\t100\t400\t0\t2\t100,100,\t0,200,
chr1\t1000\t1400\ttxB\t0\t-\t1000\t1400\t0\t1\t400,\t0,
chr2\t100\t400\ttxC\t0\t+\t100\t400\t0\t1\t300,\t0,
";

    /// A coordinate-sorted pair, plus enough reads between the transcripts that
    /// the whole-file path and the windowed path see the same input.
    fn pair(pos: usize, mate: usize, i: usize) -> sam::alignment::record_buf::RecordBuf {
        sam::alignment::record_buf::RecordBuf::builder()
            .set_name(format!("r{i}"))
            .set_flags(sam::alignment::record::Flags::SEGMENTED | sam::alignment::record::Flags::FIRST_SEGMENT)
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(pos).unwrap())
            .set_mate_alignment_start(noodles_core::Position::try_from(mate).unwrap())
            .set_mapping_quality(sam::alignment::record::MappingQuality::new(40).unwrap())
            .set_cigar(sam::alignment::record_buf::Cigar::from(vec![
                sam::alignment::record::cigar::Op::new(
                    sam::alignment::record::cigar::op::Kind::Match,
                    30,
                ),
            ]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 30]))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![40; 30]))
            .build()
    }

    /// Runs both drivers over the same records and returns (whole_file, windowed)
    /// rendered outputs. The whole-file path is the reference: it is what every
    /// recorded differential result was produced with.
    fn both_paths(
        header: &sam::Header,
        records: &[sam::alignment::record_buf::RecordBuf],
    ) -> (String, String, (f64, f64), StreamingTotals) {
        let bam_records = to_bam_records(header, records);
        let reads_by_chrom = build_read_index(bam_records.iter().cloned().map(Ok), header, false, 0).unwrap();
        let global = build_global_exon_ranges(Cursor::new(BED)).unwrap();
        let totals_whole = count_total_fragments(&reads_by_chrom, &global, 1.0);
        let totals_stream = count_total_fragments_streaming(
            bam_records.iter().cloned().map(Ok),
            header,
            &global,
            1.0,
            false,
            0,
        )
        .unwrap();

        let whole = compute_fpkm_rows(
            Cursor::new(BED),
            &reads_by_chrom,
            false,
            &HashMap::new(),
            1.0,
            1.0,
            |_| {},
        )
        .unwrap();
        let windowed = match compute_fpkm_rows_windowed(
            bam_records.iter().cloned().map(Ok),
            header,
            Cursor::new(BED),
            false,
            &HashMap::new(),
            1.0,
            1.0,
            false,
            0,
            |_| {},
        )
        .unwrap()
        {
            WindowedFpkm::Computed(rows) => rows,
            WindowedFpkm::NotCoordinateSorted => panic!("fixture is coordinate-sorted"),
        };

        (
            render_fpkm_xls(&whole),
            render_fpkm_xls(&windowed),
            totals_whole,
            totals_stream,
        )
    }

    #[test]
    fn windowed_driver_matches_the_whole_file_path_row_for_row() {
        let header = header_with_chrom("chr1", 100_000);
        // Coordinate-sorted: pos ascending. Reads land inside, before, after and
        // BETWEEN the transcripts, so the window has to retire reads the way the
        // whole-file prefix scan implicitly does.
        let mut records = Vec::new();
        for (i, pos) in [1usize, 50, 150, 250, 350, 500, 900, 1100, 1300, 1500, 20_000, 50_000]
            .into_iter()
            .enumerate()
        {
            records.push(pair(pos, pos + 200, i));
        }

        let (whole, windowed, totals_whole, totals_stream) = both_paths(&header, &records);

        assert_eq!(
            totals_whole,
            (totals_stream.total_frags, totals_stream.exonic_frags),
            "streaming totals must equal the whole-file totals"
        );
        assert!(
            totals_stream.coordinate_sorted,
            "the fixture is coordinate-sorted and must be recognised as such"
        );
        assert_eq!(whole, windowed, "windowed rows must match the whole-file rows exactly");
        assert!(whole.contains("txA"), "fixture must actually produce rows: {whole}");
    }

    #[test]
    fn windowed_driver_preserves_input_row_order_not_coordinate_order() {
        // txB starts after txA but is written second; a coordinate-ordered
        // driver that emitted rows as it scored them would flip them.
        let header = header_with_chrom("chr1", 100_000);
        let records = vec![pair(120, 320, 0), pair(1020, 1220, 1)];

        let (whole, windowed, ..) = both_paths(&header, &records);

        assert_eq!(whole, windowed);
        let a = whole.find("txA").unwrap();
        let b = whole.find("txB").unwrap();
        assert!(a < b, "txA is written first in the BED and must be emitted first");
    }

    #[test]
    fn streaming_totals_flag_an_unsorted_stream() {
        // Regression for the false negative that made a shuffled BAM score 1729
        // of 3000 rows as zero: the windowed driver's `start < last_start` probe
        // only fires for a record the pull loop actually reaches, and a record
        // past the current `tx_end` makes the loop break before comparing. The
        // totals pass sees every record, so the flag has to come from there.
        let header = header_with_chrom("chr1", 100_000);
        let global = build_global_exon_ranges(Cursor::new(BED)).unwrap();

        // Sorted: offsets ascending.
        let sorted = vec![pair(100, 300, 0), pair(500, 700, 1), pair(900, 1100, 2)];
        let got = count_total_fragments_streaming(
            to_bam_records(&header, &sorted).iter().cloned().map(Ok),
            &header,
            &global,
            1.0,
            false,
            0,
        )
        .unwrap();
        assert!(got.coordinate_sorted, "an ascending stream must be accepted");

        // A record far past its predecessor, then one that goes backwards.
        let shuffled = vec![pair(100, 300, 0), pair(90_000, 90_200, 1), pair(500, 700, 2)];
        let got = count_total_fragments_streaming(
            to_bam_records(&header, &shuffled).iter().cloned().map(Ok),
            &header,
            &global,
            1.0,
            false,
            0,
        )
        .unwrap();
        assert!(
            !got.coordinate_sorted,
            "a backwards start must be flagged even when no transcript window reaches it"
        );
    }

    #[test]
    fn windowed_driver_detects_an_out_of_order_start() {
        let header = header_with_chrom("chr1", 100_000);
        // 150 then 500 then 300. The violation only becomes visible once a
        // record is actually PULLED, so the fixture needs a valid read first and
        // the out-of-order one must fall inside a later transcript's window --
        // a bare [2000, 500] pair never gets pulled, because the pull loop stops
        // at the first start past `tx_end` before it can compare order.
        let records = vec![pair(150, 350, 0), pair(500, 700, 1), pair(300, 500, 2)];
        let bam_records = to_bam_records(&header, &records);

        let got = compute_fpkm_rows_windowed(
            bam_records.iter().cloned().map(Ok),
            &header,
            Cursor::new(BED),
            false,
            &HashMap::new(),
            1.0,
            1.0,
            false,
            0,
            |_| {},
        )
        .unwrap();

        assert!(
            matches!(got, WindowedFpkm::NotCoordinateSorted),
            "an out-of-order start must send the caller to the whole-file fallback"
        );
    }

    #[test]
    fn windowed_driver_skips_a_transcript_on_an_absent_chromosome() {
        // txC is on chr2, which the header does not contain. The whole-file path
        // looks the chromosome up, finds nothing, and still emits a ZERO row --
        // it does not drop the transcript. The windowed driver has to reproduce
        // that, not skip it, or the two drivers produce different row counts.
        let header = header_with_chrom("chr1", 100_000);
        let records = vec![pair(120, 320, 0)];

        let (whole, windowed, ..) = both_paths(&header, &records);

        assert_eq!(whole, windowed);
        assert!(whole.contains("txC"), "an absent chromosome still yields a zero row");
        let tx_c = whole.lines().find(|l| l.contains("txC")).unwrap();
        assert!(
            tx_c.ends_with("\t0.0\t0.0\t0.0"),
            "the absent-chromosome row must be all zeros, got {tx_c:?}"
        );
    }
}
