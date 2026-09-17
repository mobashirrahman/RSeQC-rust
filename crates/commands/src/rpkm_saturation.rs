//! Port of `RPKM_saturation.py`: assess whether transcript RPKM
//! estimates have reached sequencing saturation, plus the quartile
//! boxplot R script that summarizes the result. Algorithm ported from
//! `ParseBAM.saturation_RPKM` in `oracle/upstream-src/src/qcmodule/
//! SAM.py` (lines 4061-4236, the `q_cut`-taking variant used by the CLI
//! script) and the CLI script's own `show_saturation`/`parse_rpkm_table`/
//! `square_error`.
//!
//! No BAI index support needed here (unlike `tin.py`/`FPKM_count.py`):
//! the read scan is purely sequential, matching upstream's own
//! `next(self.samfile)` loop.
//!
//! SAM-text input (DIV-0002/0004) is supported at the CLI layer via
//! `rseqc_formats::open_alignments`; this module's own functions were
//! unaffected (already generic over
//! `IntoIterator<Item = io::Result<bam::Record>>`).
//!
//! **Preserves several deliberate upstream quirks, do not "fix"**:
//! - Exon coordinates contribute a SINGLE MIDPOINT per exon block
//!   (`exon_start + (exon_len // 2)`), not the full exon span, to the
//!   sampling population -- a spliced read with N exon blocks
//!   contributes N independent midpoints.
//! - `sample_size` (the RPKM denominator) always uses the OVERALL total
//!   exon-block count (`cUR_num`), even in strand-specific mode where
//!   the actual sampled ranges are built from the strand-specific counts
//!   (`cUR_plus`/`cUR_minus`) -- normalization never becomes
//!   strand-specific, only the numerator does.
//! - The per-transcript "how many sampled points fall in this exon"
//!   query counts RAW OVERLAPPING POINTS (duplicates included), unlike
//!   most other ports' presence-only exon-overlap checks -- implemented
//!   here with a small local `PointCounts` (sorted-coordinate binary
//!   search per chromosome), not the shared `MergedRegions` (which would
//!   silently collapse overlapping points and undercount).
//! - If the computed strand key (e.g. `"1+"`) has no entry in a
//!   strand-specific rule's map, upstream crashes with an uncaught
//!   `KeyError` (a user error: rule/data pairing-mode mismatch) rather
//!   than silently treating the read as unassigned -- replicated here as
//!   a propagated `io::Error`, not a silent skip.
//! - `.eRPKM.xls`/`.rawCount.xls` row text has a literal SPACE between
//!   the key's trailing tab and the first value (`print(key+'\t',
//!   end=' ')` immediately followed by another `print(...)` call) --
//!   same class of `print(..., end=' ')` quirk seen in
//!   `junction_annotation.py`.
//! - The gene-model BED12 parser here is a SIXTH distinct robustness
//!   tier across this port: a bare, unqualified `except:` around the
//!   whole per-line block (broader than the `(IndexError, ValueError)`-
//!   scoped tiers used elsewhere).
//! - A transcript whose computed `mRNA_len` is `0`, or a percentile
//!   whose `sample_size` is `0`, is an unconditional hard abort
//!   (`sys.exit(1)`, which is NOT caught by the CLI's own
//!   `except (OSError, ValueError, RuntimeError, IndexError)` clause) --
//!   replicated as a propagated `io::Error`.

use std::collections::{HashMap, HashSet};
use std::io::{self, BufRead};

use noodles_bam as bam;
use noodles_sam as sam;
use rand::seq::SliceRandom;
use rseqc_formats::cigar::fetch_exon_blocks;

use crate::python_fmt::python_str_float;

#[derive(Debug, Default, Clone)]
pub struct BlockLists {
    pub block_list: Vec<(String, i64)>,
    pub block_list_plus: Vec<(String, i64)>,
    pub block_list_minus: Vec<(String, i64)>,
    pub cur_num: i64,
    pub cur_plus: i64,
    pub cur_minus: i64,
}

/// Builds the exon-midpoint sampling population from a BAM. Ports the
/// `saturation_RPKM` read-scan loop exactly (see module docs for the
/// strand-key/KeyError-crash quirk).
pub fn build_block_lists<I>(
    records: I,
    header: &sam::Header,
    skip_multi: bool,
    map_qual: u8,
    strand_rule_active: bool,
    strand_map: &HashMap<String, char>,
) -> io::Result<BlockLists>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut lists = BlockLists::default();

    for result in records {
        let record = result?;
        let flags = record.flags();
        if flags.is_qc_fail() || flags.is_duplicate() || flags.is_secondary() || flags.is_unmapped() {
            continue;
        }
        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        if skip_multi && mapq < map_qual {
            continue;
        }

        let Some(ref_id) = record.reference_sequence_id().transpose()? else { continue };
        let Some((chrom_bstr, _)) = header.reference_sequences().get_index(ref_id) else { continue };
        let chrom = chrom_bstr.to_string().to_uppercase();

        let Some(pos) = record.alignment_start().transpose()? else { continue };
        let hit_st = (pos.get() - 1) as i64;

        let read_id = if flags.is_segmented() {
            if flags.is_first_segment() {
                "1"
            } else if flags.is_last_segment() {
                "2"
            } else {
                ""
            }
        } else {
            ""
        };
        let map_strand = if flags.is_reverse_complemented() { "-" } else { "+" };
        let strand_key = format!("{read_id}{map_strand}");

        let ops: Vec<_> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;
        let exon_blocks = fetch_exon_blocks(hit_st as usize, ops.iter().copied());
        lists.cur_num += exon_blocks.len() as i64;

        if strand_rule_active {
            let assigned = *strand_map.get(&strand_key).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("strand_rule has no mapping for computed key {strand_key:?} (rule/pairing-mode mismatch)"),
                )
            })?;
            if assigned == '+' {
                lists.cur_plus += exon_blocks.len() as i64;
            }
            if assigned == '-' {
                lists.cur_minus += exon_blocks.len() as i64;
            }
            for &(s, e) in &exon_blocks {
                let (s, e) = (s as i64, e as i64);
                let midpoint = s + (e - s) / 2;
                if assigned == '+' {
                    lists.block_list_plus.push((chrom.clone(), midpoint));
                }
                if assigned == '-' {
                    lists.block_list_minus.push((chrom.clone(), midpoint));
                }
            }
        } else {
            for &(s, e) in &exon_blocks {
                let (s, e) = (s as i64, e as i64);
                let midpoint = s + (e - s) / 2;
                lists.block_list.push((chrom.clone(), midpoint));
            }
        }
    }

    Ok(lists)
}

pub fn shuffle_block_lists(lists: &mut BlockLists, rng: &mut impl rand::Rng) {
    lists.block_list_plus.shuffle(rng);
    lists.block_list_minus.shuffle(rng);
    lists.block_list.shuffle(rng);
}

/// Counts of raw, possibly-overlapping single-base points per
/// chromosome, queryable by range -- NOT a merged/deduplicated interval
/// set (see module docs for why `MergedRegions` would be wrong here).
struct PointCounts {
    by_chrom: HashMap<String, Vec<i64>>,
}

impl PointCounts {
    fn new(points: &[(String, i64)]) -> Self {
        let mut by_chrom: HashMap<String, Vec<i64>> = HashMap::new();
        for (chrom, coord) in points {
            by_chrom.entry(chrom.clone()).or_default().push(*coord);
        }
        for v in by_chrom.values_mut() {
            v.sort_unstable();
        }
        Self { by_chrom }
    }

    fn count_in_range(&self, chrom: &str, start: i64, end: i64) -> i64 {
        if start >= end {
            return 0;
        }
        let Some(coords) = self.by_chrom.get(chrom) else { return 0 };
        let lo = coords.partition_point(|&c| c < start);
        let hi = coords.partition_point(|&c| c < end);
        (hi - lo) as i64
    }
}

fn slice_range(list: &[(String, i64)], lo: i64, hi: i64) -> &[(String, i64)] {
    let len = list.len() as i64;
    let lo = lo.clamp(0, len) as usize;
    let hi = hi.clamp(0, len) as usize;
    if lo >= hi { &[] } else { &list[lo..hi] }
}

#[derive(Debug, Clone)]
pub struct TranscriptExons {
    pub key: String,
    pub gene_name: String,
    pub chrom_upper: String,
    pub strand: String,
    pub exon_starts: Vec<i64>,
    pub exon_ends: Vec<i64>,
    pub mrna_len: i64,
}

/// Parses the reference BED12 gene model. Robustness tier: comment/
/// track/browser lines skipped; ANY other per-line parse failure is
/// caught and skipped with a stderr note -- upstream's bare `except:`.
pub fn parse_gene_model(reader: impl BufRead) -> io::Result<Vec<TranscriptExons>> {
    let mut out = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }

        let parsed: Option<TranscriptExons> = (|| {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 12 {
                return None;
            }
            let chrom = fields[0].to_string();
            let tx_start: i64 = fields[1].parse().ok()?;
            let tx_end: i64 = fields[2].parse().ok()?;
            let gene_name = fields[3].to_string();
            let strand = fields[5].to_string();

            let exon_starts: Vec<i64> = fields[11]
                .trim_end_matches(',')
                .split(',')
                .map(|s| s.parse::<i64>().ok().map(|v| v + tx_start))
                .collect::<Option<_>>()?;
            let exon_sizes: Vec<i64> =
                fields[10].trim_end_matches(',').split(',').map(|s| s.parse::<i64>().ok()).collect::<Option<_>>()?;
            let exon_ends: Vec<i64> = exon_starts.iter().zip(exon_sizes.iter()).map(|(&s, &sz)| s + sz).collect();
            let mrna_len: i64 = exon_sizes.iter().sum();

            let key = format!("{}\t{}\t{}\t{}\t0\t{}", chrom.to_lowercase(), tx_start, tx_end, gene_name, strand);

            Some(TranscriptExons { key, gene_name, chrom_upper: chrom.to_uppercase(), strand, exon_starts, exon_ends, mrna_len })
        })();

        match parsed {
            Some(t) => out.push(t),
            None => eprintln!("[NOTE:input bed must be 12-column] skipped this line: {line}"),
        }
    }

    Ok(out)
}

#[derive(Debug, Default, Clone)]
pub struct SaturationResult {
    /// Percentile column labels, e.g. `"5%"`, in sampling order.
    pub header: Vec<String>,
    /// Transcript keys in first-seen order (stable output row order).
    pub keys_in_order: Vec<String>,
    pub rpkm_table: HashMap<String, Vec<f64>>,
    pub raw_table: HashMap<String, Vec<i64>>,
}

/// Runs the cumulative percentile resampling loop, computing RPKM for
/// every transcript at every percentile. Ports the main loop of
/// `saturation_RPKM` (lines 4157-4227) exactly.
///
/// **Critical: the sampled point population is CUMULATIVE across
/// percentile iterations, not rebuilt fresh each time.** Upstream's
/// `ranges`/`ranges_plus`/`ranges_minus` dicts are declared ONCE before
/// the `for pertl in tmp:` loop and are never cleared inside it -- each
/// iteration's `[int(cUR_num*percent_st):int(cUR_num*percent_end)]`
/// slice is ADDED on top of whatever earlier iterations already
/// inserted (same accumulation pattern as `junction_saturation.py`'s
/// `uniqSpliceSites`). Confirmed by reading the actual loop body, not
/// inferred from the earlier percentile-list docs. Getting this wrong
/// (rebuilding an independent per-percentile point set, as an earlier
/// version of this port did) produces a completely different, WRONG,
/// non-monotonic-looking saturation curve -- verified by an actual live
/// diff against the real upstream CLI, not just reasoning about the
/// source.
pub fn compute_saturation(
    transcripts: &[TranscriptExons],
    lists: &BlockLists,
    strand_rule_active: bool,
    sample_start: i64,
    sample_end: i64,
    sample_step: i64,
    refbed_path: &str,
) -> io::Result<SaturationResult> {
    let mut percentiles = Vec::new();
    let mut p = sample_start;
    while p < sample_end {
        percentiles.push(p);
        p += sample_step;
    }
    percentiles.push(100);

    let mut result = SaturationResult::default();
    let mut seen_keys: HashSet<String> = HashSet::new();

    // Cumulative point buffers, grown (never cleared) across iterations,
    // matching upstream's persistent `ranges`/`ranges_plus`/`ranges_minus`.
    let mut cumulative_plus: Vec<(String, i64)> = Vec::new();
    let mut cumulative_minus: Vec<(String, i64)> = Vec::new();
    let mut cumulative_all: Vec<(String, i64)> = Vec::new();

    for &pertl in &percentiles {
        let percent_st = (((pertl - sample_step) as f64) / 100.0).max(0.0);
        let percent_end = (pertl as f64) / 100.0;
        let sample_size = lists.cur_num as f64 * percent_end;
        result.header.push(format!("{pertl}%"));

        let (ranges_plus, ranges_minus, ranges) = if strand_rule_active {
            let lo_p = (lists.cur_plus as f64 * percent_st) as i64;
            let hi_p = (lists.cur_plus as f64 * percent_end) as i64;
            eprintln!("sampling {pertl}% ({}) forward strand fragments ...", (lists.cur_plus as f64 * percent_end) as i64);
            cumulative_plus.extend_from_slice(slice_range(&lists.block_list_plus, lo_p, hi_p));

            let lo_m = (lists.cur_minus as f64 * percent_st) as i64;
            let hi_m = (lists.cur_minus as f64 * percent_end) as i64;
            eprintln!("sampling {pertl}% ({}) reverse strand fragments ...", (lists.cur_minus as f64 * percent_end) as i64);
            cumulative_minus.extend_from_slice(slice_range(&lists.block_list_minus, lo_m, hi_m));

            let plus = PointCounts::new(&cumulative_plus);
            let minus = PointCounts::new(&cumulative_minus);
            (Some(plus), Some(minus), None)
        } else {
            let lo = (lists.cur_num as f64 * percent_st) as i64;
            let hi = (lists.cur_num as f64 * percent_end) as i64;
            eprintln!("sampling {pertl}% ({}) fragments ...", sample_size as i64);
            cumulative_all.extend_from_slice(slice_range(&lists.block_list, lo, hi));
            let all = PointCounts::new(&cumulative_all);
            (None, None, Some(all))
        };

        eprintln!("assign reads to transcripts in {refbed_path} ...");

        for t in transcripts {
            let mut mrna_count = 0i64;
            for (&s, &e) in t.exon_starts.iter().zip(t.exon_ends.iter()) {
                if strand_rule_active {
                    if t.strand == "+" {
                        if let Some(pr) = &ranges_plus {
                            mrna_count += pr.count_in_range(&t.chrom_upper, s, e);
                        }
                    }
                    if t.strand == "-" {
                        if let Some(mr) = &ranges_minus {
                            mrna_count += mr.count_in_range(&t.chrom_upper, s, e);
                        }
                    }
                } else if let Some(r) = &ranges {
                    mrna_count += r.count_in_range(&t.chrom_upper, s, e);
                }
            }

            if t.mrna_len == 0 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, format!("{} has 0 nucleotides. Exit!", t.gene_name)));
            }
            if sample_size == 0.0 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "Too few reads to sample. Exit!"));
            }

            let mrna_rpkm = (mrna_count as f64 * 1_000_000_000.0) / (t.mrna_len as f64 * sample_size);

            if seen_keys.insert(t.key.clone()) {
                result.keys_in_order.push(t.key.clone());
            }
            result.rpkm_table.entry(t.key.clone()).or_default().push(mrna_rpkm);
            result.raw_table.entry(t.key.clone()).or_default().push(mrna_count);
        }
        // Upstream: `print("", file=sys.stderr)` -- a bare blank line
        // after each percentile's transcript-assignment pass.
        eprintln!();
    }

    Ok(result)
}

/// Renders the `.eRPKM.xls` table text, including the literal
/// space-after-key quirk described in the module docs.
pub fn render_rpkm_xls(result: &SaturationResult) -> String {
    render_table(result, python_str_float)
}

/// Renders the `.rawCount.xls` table text (plain integer counts).
pub fn render_raw_xls(result: &SaturationResult) -> String {
    let mut out = String::new();
    let head: Vec<&str> = ["#chr", "start", "end", "name", "score", "strand"].into_iter().collect();
    out.push_str(&head.join("\t"));
    for h in &result.header {
        out.push('\t');
        out.push_str(h);
    }
    out.push('\n');
    for key in &result.keys_in_order {
        let values = &result.raw_table[key];
        let joined: Vec<String> = values.iter().map(|v| v.to_string()).collect();
        out.push_str(key);
        out.push_str("\t ");
        out.push_str(&joined.join("\t"));
        out.push('\n');
    }
    out
}

fn render_table(result: &SaturationResult, fmt: impl Fn(f64) -> String) -> String {
    let mut out = String::new();
    let head: Vec<&str> = ["#chr", "start", "end", "name", "score", "strand"].into_iter().collect();
    out.push_str(&head.join("\t"));
    for h in &result.header {
        out.push('\t');
        out.push_str(h);
    }
    out.push('\n');
    for key in &result.keys_in_order {
        let values = &result.rpkm_table[key];
        let joined: Vec<String> = values.iter().map(|&v| fmt(v)).collect();
        out.push_str(key);
        out.push_str("\t ");
        out.push_str(&joined.join("\t"));
        out.push('\n');
    }
    out
}

/// Ports `square_error`: relative error of each sampled value against
/// the final (100%, full-depth) value, or `None` if the true value or
/// the value range is zero.
pub fn square_error(values: &[f64]) -> Option<Vec<f64>> {
    let true_rpkm = *values.last()?;
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
    if true_rpkm == 0.0 || (max - min) == 0.0 {
        return None;
    }
    Some(values.iter().map(|v| (v - true_rpkm).abs() / true_rpkm).collect())
}

#[derive(Debug, Default, Clone)]
pub struct QuartilePlotData {
    /// Percentile labels WITHOUT the trailing `%`, in header order.
    pub header: Vec<String>,
    pub keys_in_order: Vec<String>,
    pub rpkm_errors: HashMap<String, Vec<f64>>,
    pub rpkm_means: HashMap<String, f64>,
}

/// Ports `parse_rpkm_table`'s filtering, operating directly on the
/// in-memory `SaturationResult` rather than re-reading the `.eRPKM.xls`
/// file upstream writes and re-parses -- a lossless simplification (the
/// file is an exact serialization of this same data), not a behavioral
/// change.
pub fn build_quartile_plot_data(result: &SaturationResult, rpkm_cutoff: f64) -> QuartilePlotData {
    let header: Vec<String> = result.header.iter().map(|h| h.trim_end_matches('%').to_string()).collect();
    let mut data = QuartilePlotData { header, ..Default::default() };

    for key in &result.keys_in_order {
        let values = &result.rpkm_table[key];
        let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
        if max == 0.0 {
            continue;
        }
        if max - min == 0.0 {
            continue;
        }
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        if mean < rpkm_cutoff {
            continue;
        }
        let Some(errors) = square_error(values) else { continue };

        data.keys_in_order.push(key.clone());
        data.rpkm_errors.insert(key.clone(), errors);
        data.rpkm_means.insert(key.clone(), mean);
    }

    data
}

/// Ports `show_saturation`'s quartile boxplot R script generation.
/// `out_prefix` is the CLI's `-o` value (matching upstream's
/// `outfile.with_suffix(".pdf")` on the `{out_prefix}.saturation.r`
/// path, i.e. the PDF is `{out_prefix}.saturation.pdf`).
pub fn render_saturation_r_script(data: &QuartilePlotData, out_prefix: &str) -> io::Result<String> {
    let gene_count = data.keys_in_order.len();
    if gene_count == 0 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "no transcripts passed the saturation-plot filtering criteria"));
    }

    let mut sorted_genes: Vec<(&String, f64)> = data.keys_in_order.iter().map(|k| (k, data.rpkm_means[k])).collect();
    sorted_genes.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let plot_headers: &[String] = if !data.header.is_empty() { &data.header[..data.header.len() - 1] } else { &[] };
    if plot_headers.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "at least two sampling percentages are required"));
    }

    let mut out = String::new();
    out.push_str(&format!("pdf('{out_prefix}.saturation.pdf')\n"));
    out.push_str("par(mfrow=c(2,2))\n");

    for (label, lo, hi) in [("Q1", 0.0, 0.25), ("Q2", 0.25, 0.50), ("Q3", 0.50, 0.75), ("Q4", 0.75, 1.00)] {
        let mut normalized_rpkm: HashMap<&str, Vec<String>> = HashMap::new();

        for (i, (gene_key, _)) in sorted_genes.iter().enumerate() {
            let line_count = (i + 1) as f64;
            if line_count > gene_count as f64 * lo && line_count <= gene_count as f64 * hi {
                let errors = &data.rpkm_errors[*gene_key];
                for (idx, e) in errors.iter().enumerate() {
                    normalized_rpkm.entry(data.header[idx].as_str()).or_default().push(python_str_float(*e));
                }
            }
        }

        let name_vec: Vec<String> = plot_headers.iter().map(|h| format!("'{h}'")).collect();
        out.push_str(&format!("name=c({})\n", name_vec.join(",")));

        for sample in plot_headers {
            let values = normalized_rpkm.get(sample.as_str()).cloned().unwrap_or_default();
            out.push_str(&format!("S{sample}=c({})\n", values.join(",")));
        }

        let series: Vec<String> = plot_headers.iter().map(|s| format!("100*S{s}")).collect();
        out.push_str(&format!(
            "boxplot({},names=name,outline=FALSE,ylab='Percent Relative Error',main='{label}',xlab='Resampling percentage')\n",
            series.join(",")
        ));
    }

    out.push_str("dev.off()\n");
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn parse_gene_model_lowercases_key_chrom_and_skips_short_lines() {
        let bed = "\
track name=x
too short
chr1\t0\t300\ttx1\t0\t+\t0\t300\t0\t2\t100,100,\t0,200,
";
        let transcripts = parse_gene_model(Cursor::new(bed)).unwrap();
        assert_eq!(transcripts.len(), 1);
        let t = &transcripts[0];
        assert_eq!(t.key, "chr1\t0\t300\ttx1\t0\t+");
        assert_eq!(t.chrom_upper, "CHR1");
        assert_eq!(t.exon_starts, vec![0, 200]);
        assert_eq!(t.exon_ends, vec![100, 300]);
        assert_eq!(t.mrna_len, 200);
    }

    #[test]
    fn compute_saturation_percentile_list_always_ends_with_literal_100() {
        let lists = BlockLists::default();
        let result = compute_saturation(&[], &lists, false, 5, 50, 10, "refgene.bed").unwrap();
        assert_eq!(result.header, vec!["5%", "15%", "25%", "35%", "45%", "100%"]);
    }

    #[test]
    fn compute_saturation_counts_raw_overlapping_points_not_merged() {
        // Two midpoints at the SAME coordinate (100) inside one exon
        // [50,200) -- a merged/presence-only check would undercount this
        // as 1; the real upstream Intersecter.find() counts both.
        let lists = BlockLists { block_list: vec![("CHR1".to_string(), 100), ("CHR1".to_string(), 100)], cur_num: 2, ..Default::default() };
        let t = TranscriptExons {
            key: "chr1\t50\t200\ttx1\t0\t+".to_string(),
            gene_name: "tx1".to_string(),
            chrom_upper: "CHR1".to_string(),
            strand: "+".to_string(),
            exon_starts: vec![50],
            exon_ends: vec![200],
            mrna_len: 150,
        };
        // percentiles: range(100,100,100) is empty, tmp=[100] only ->
        // single 100% sample. sample_step=100 makes percent_st =
        // (100-100)/100.0 = 0.0, so the slice covers the FULL population
        // (lo=0,hi=cur_num) -- matches upstream's `percent_st =
        // (pertl-sample_step)/100.0` exactly.
        let result = compute_saturation(&[t], &lists, false, 100, 100, 100, "refgene.bed").unwrap();
        assert_eq!(result.header, vec!["100%"]);
        let raw = &result.raw_table["chr1\t50\t200\ttx1\t0\t+"];
        assert_eq!(raw, &vec![2]);
        // RPKM = 2 * 1e9 / (150 * (2*1.0)) = 2e9/300 = 6666666.666...
        let rpkm = result.rpkm_table["chr1\t50\t200\ttx1\t0\t+"][0];
        assert!((rpkm - (2_000_000_000.0 / 300.0)).abs() < 1e-6);
    }

    #[test]
    fn compute_saturation_point_population_is_cumulative_across_percentiles() {
        // Critical regression guard: upstream's `ranges` dict is built
        // ONCE before the percentile loop and never cleared, so each
        // iteration's mRNA_count reflects ALL points sampled so far, not
        // just that iteration's own slice. With 4 points and 2 evenly
        // spaced percentiles (50%, 100%), the FIRST iteration's slice is
        // points[0..2] and the SECOND iteration's slice is points[2..4]
        // -- but because the population is cumulative, the second
        // iteration's mRNA_count must be 4 (all points), not 2 (just its
        // own slice). Getting this wrong was a real, previously-shipped
        // bug found via a live diff against the real upstream CLI.
        let lists = BlockLists {
            block_list: vec![
                ("CHR1".to_string(), 60),
                ("CHR1".to_string(), 70),
                ("CHR1".to_string(), 80),
                ("CHR1".to_string(), 90),
            ],
            cur_num: 4,
            ..Default::default()
        };
        let t = TranscriptExons {
            key: "chr1\t50\t200\ttx1\t0\t+".to_string(),
            gene_name: "tx1".to_string(),
            chrom_upper: "CHR1".to_string(),
            strand: "+".to_string(),
            exon_starts: vec![50],
            exon_ends: vec![200],
            mrna_len: 150,
        };
        // percentiles: range(50,100,50) -> [50], plus literal 100 -> [50,100].
        let result = compute_saturation(&[t], &lists, false, 50, 100, 50, "refgene.bed").unwrap();
        assert_eq!(result.header, vec!["50%", "100%"]);
        let raw = &result.raw_table["chr1\t50\t200\ttx1\t0\t+"];
        // 50%: slice [0,2) -> 2 points. 100%: CUMULATIVE union of
        // [0,2) and [2,4) -> all 4 points, not just the second slice's 2.
        assert_eq!(raw, &vec![2, 4]);
    }

    #[test]
    fn compute_saturation_errors_on_zero_mrna_len() {
        let lists = BlockLists { cur_num: 10, ..Default::default() };
        let t = TranscriptExons {
            key: "k".to_string(),
            gene_name: "tx0".to_string(),
            chrom_upper: "CHR1".to_string(),
            strand: "+".to_string(),
            exon_starts: vec![],
            exon_ends: vec![],
            mrna_len: 0,
        };
        let err = compute_saturation(&[t], &lists, false, 100, 100, 1, "refgene.bed").unwrap_err();
        assert!(err.to_string().contains("tx0"));
    }

    #[test]
    fn square_error_none_when_true_rpkm_or_range_is_zero() {
        assert_eq!(square_error(&[1.0, 2.0, 0.0]), None); // true_rpkm==0
        assert_eq!(square_error(&[5.0, 5.0, 5.0]), None); // range==0
        let errors = square_error(&[2.0, 8.0, 4.0]).unwrap();
        // true_rpkm=4.0: |2-4|/4=0.5, |8-4|/4=1.0, |4-4|/4=0.0
        assert_eq!(errors, vec![0.5, 1.0, 0.0]);
    }

    #[test]
    fn render_rpkm_xls_has_space_after_key_tab() {
        let mut result = SaturationResult { header: vec!["5%".to_string(), "100%".to_string()], ..Default::default() };
        result.keys_in_order.push("chr1\t0\t100\ttx1\t0\t+".to_string());
        result.rpkm_table.insert("chr1\t0\t100\ttx1\t0\t+".to_string(), vec![1.0, 2.5]);
        let xls = render_rpkm_xls(&result);
        assert_eq!(
            xls,
            "#chr\tstart\tend\tname\tscore\tstrand\t5%\t100%\nchr1\t0\t100\ttx1\t0\t+\t 1.0\t2.5\n"
        );
    }

    #[test]
    fn render_saturation_r_script_exact_text() {
        // Two genes, two percentile columns (last one dropped from
        // plot_headers). Verified against a hand-computed rendering of
        // upstream's literal print()/%-format lines with these inputs.
        let mut data = QuartilePlotData { header: vec!["5".to_string(), "100".to_string()], ..Default::default() };
        data.keys_in_order.push("g1".to_string());
        data.keys_in_order.push("g2".to_string());
        data.rpkm_errors.insert("g1".to_string(), vec![0.1, 0.0]);
        data.rpkm_errors.insert("g2".to_string(), vec![0.5, 0.0]);
        data.rpkm_means.insert("g1".to_string(), 10.0);
        data.rpkm_means.insert("g2".to_string(), 20.0);

        let script = render_saturation_r_script(&data, "out").unwrap();
        // Cross-verified byte-for-byte against a `python3 -c` run of the
        // literal upstream show_saturation() expressions with these same
        // inputs (gene_count=2, g1 mean=10 sorted first, g2 mean=20).
        let expected = "\
pdf('out.saturation.pdf')
par(mfrow=c(2,2))
name=c('5')
S5=c()
boxplot(100*S5,names=name,outline=FALSE,ylab='Percent Relative Error',main='Q1',xlab='Resampling percentage')
name=c('5')
S5=c(0.1)
boxplot(100*S5,names=name,outline=FALSE,ylab='Percent Relative Error',main='Q2',xlab='Resampling percentage')
name=c('5')
S5=c()
boxplot(100*S5,names=name,outline=FALSE,ylab='Percent Relative Error',main='Q3',xlab='Resampling percentage')
name=c('5')
S5=c(0.5)
boxplot(100*S5,names=name,outline=FALSE,ylab='Percent Relative Error',main='Q4',xlab='Resampling percentage')
dev.off()
";
        assert_eq!(script, expected);
    }
}
