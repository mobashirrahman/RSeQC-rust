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
#[derive(Debug, Clone)]
pub struct IndexedRead {
    pub start: i64,
    /// True CIGAR reference end used to emulate htslib fetch selection.
    pub end: i64,
    pub rlen: i64,
    /// `pnext`, 0-based; `-1` if unset (pysam's raw sentinel), matching
    /// upstream's use of the raw field regardless of pairing state.
    pub mate_start: i64,
    pub is_paired: bool,
    pub is_read2: bool,
    pub is_reverse: bool,
    pub is_unmapped: bool,
    pub mate_is_unmapped: bool,
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
        let flags = record.flags();
        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        if !alignment_passes_legacy_filters(flags, mapq, skip_multi, map_qual) {
            continue;
        }

        let Some(ref_id) = record.reference_sequence_id().transpose()? else { continue };
        let Some((chrom_bstr, _)) = header.reference_sequences().get_index(ref_id) else { continue };
        let chrom = chrom_bstr.to_string();

        let rlen = record.sequence().len() as i64;

        // For unmapped reads, use pos directly from the record (the coordinate field
        // still contains meaningful data for unpaired unmapped reads or when the mate
        // is mapped). For mapped reads, use alignment_start from CIGAR. This matches
        // Python's pysam behavior which exposes pos for unmapped reads too.
        let (start, reference_end) = if flags.is_unmapped() {
            // Unmapped reads: use the raw POS field as start, estimate end as start+rlen
            let pos = record.alignment_start().transpose()?.map(|p| (p.get() - 1) as i64).unwrap_or(0);
            (pos, pos + rlen)
        } else {
            // Mapped reads: use CIGAR to compute reference span
            let Some(pos) = record.alignment_start().transpose()? else { continue };
            let start = (pos.get() - 1) as i64;
            let cigar_ops: Vec<_> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;
            let (_, ref_end) = rseqc_formats::cigar::reference_span(start as usize, cigar_ops);
            (start, ref_end as i64)
        };

        let mate_start = record.mate_alignment_start().transpose()?.map(|p| (p.get() - 1) as i64).unwrap_or(-1);

        by_chrom.entry(chrom).or_default().push(IndexedRead {
            start,
            end: reference_end,
            rlen,
            mate_start,
            is_paired: flags.is_segmented(),
            is_read2: flags.is_last_segment(),
            is_reverse: flags.is_reverse_complemented(),
            is_unmapped: flags.is_unmapped(),
            mate_is_unmapped: flags.is_mate_unmapped(),
        });
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
            let read_end = r.start + r.rlen;

            if !r.is_paired {
                total_frags += 1.0;
                if global_exon_ranges.overlap_length(&chrom_upper, r.start, read_end) > 0 {
                    exonic_frags += 1.0;
                }
                continue;
            }

            if r.is_read2 {
                continue;
            }

            let mate_end = r.mate_start + r.rlen;

            if r.is_unmapped {
                if r.mate_is_unmapped {
                    continue;
                }
                total_frags += single_read;
                if global_exon_ranges.overlap_length(&chrom_upper, r.mate_start, mate_end) > 0 {
                    exonic_frags += single_read;
                }
            } else if r.mate_is_unmapped {
                total_frags += single_read;
                if global_exon_ranges.overlap_length(&chrom_upper, r.start, read_end) > 0 {
                    exonic_frags += single_read;
                }
            } else {
                total_frags += 1.0;
                // Upstream treats a paired fragment as exonic only when
                // BOTH mapped ends overlap an exon.  Checking only read 1
                // incorrectly classifies a pair whose mate lies outside the
                // gene model and inflates the exonic denominator used by
                // `-e` FPKM normalization.
                if global_exon_ranges.overlap_length(&chrom_upper, r.start, read_end) > 0
                    && global_exon_ranges.overlap_length(&chrom_upper, r.mate_start, mate_end) > 0
                {
                    exonic_frags += 1.0;
                }
            }
        }
    }

    (total_frags, exonic_frags)
}

/// Counts one transcript's forward/reverse/unstranded fragments. Ports
/// `count_transcript` exactly, including its two-separate-`if` (not
/// `if`/`elif`) structure for the unpaired/paired branches.
#[allow(clippy::too_many_arguments)]
pub fn count_transcript(
    reads_by_chrom: &HashMap<String, Vec<IndexedRead>>,
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

    let Some(reads) = reads_by_chrom.get(chrom) else { return (0.0, 0.0, 0.0) };
    let hi = reads.partition_point(|r| r.start < tx_end);

    for r in reads[..hi].iter().filter(|r| r.end > tx_start) {
        if !r.is_paired {
            let frag_st = r.start;
            let frag_end = r.start + r.rlen;
            let strand_key = if r.is_reverse { "-" } else { "+" };

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

        if r.is_paired {
            let frag_st = r.start;
            let frag_end = r.mate_start;

            let start_in_exon = exon_ranges.overlap_length(chrom, frag_st, frag_st + 1) > 0;
            let mate_in_exon = exon_ranges.overlap_length(chrom, frag_end, frag_end + 1) > 0;
            if !start_in_exon && !mate_in_exon {
                continue;
            }

            if r.is_read2 {
                continue;
            }

            let strand_key = if r.is_reverse { "1-" } else { "1+" };

            if !strand_rule_active {
                if r.is_unmapped {
                    if r.mate_is_unmapped {
                        continue;
                    }
                    frag_count_fr += single_read;
                } else if r.mate_is_unmapped {
                    frag_count_fr += single_read;
                } else {
                    frag_count_fr += 1.0;
                }
            } else {
                if strand_map.get(strand_key) == Some(&'+') {
                    if r.is_unmapped {
                        if !r.mate_is_unmapped {
                            frag_count_f += single_read;
                        }
                    } else if r.mate_is_unmapped {
                        frag_count_f += single_read;
                    } else {
                        frag_count_f += 1.0;
                    }
                }
                if strand_map.get(strand_key) == Some(&'-') {
                    if r.is_unmapped {
                        if !r.mate_is_unmapped {
                            frag_count_r += single_read;
                        }
                    } else if r.mate_is_unmapped {
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

        let (frag_f, frag_r, frag_fr) = count_transcript(
            reads_by_chrom,
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
        let (f, r, fr) = count_transcript(&index, "chr1", 0, 1000, &exon_ranges, false, &HashMap::new(), 1.0);
        assert_eq!((f, r, fr), (0.0, 0.0, 2.0));

        // Stranded "++,--" -> forward read maps to '+', reverse read maps to '-'.
        let map = parse_strand_rule(Some("++,--")).unwrap();
        let (f, r, fr) = count_transcript(&index, "chr1", 0, 1000, &exon_ranges, true, &map, 1.0);
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
}
