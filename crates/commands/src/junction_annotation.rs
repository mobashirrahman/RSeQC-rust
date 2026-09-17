//! Port of `junction_annotation.py`: annotate splice junctions against a
//! reference gene model. Contract: see `compatibility/commands.yaml` entry
//! `junction_annotation.py`. Core classification ported from
//! `annotate_junction()` in `oracle/upstream-src/src/qcmodule/SAM.py`
//! (~lines 3752-3886); `.bed`/`.Interact.bed` generation ported directly
//! from `oracle/upstream-src/scripts/junction_annotation.py`'s
//! `generate_bed12`/`generate_interact` (modern-style helpers that
//! re-derive their output from the `.junction.xls` rows -- implemented
//! here as operating on the in-memory rows directly, not a text re-parse).
//!
//! **Preserves upstream's known-junction test being two INDEPENDENT
//! membership checks, not a paired lookup**: `refIntronStarts`/
//! `refIntronEnds` are separate per-chromosome sets built from ALL
//! transcripts' intron boundaries; a junction is "known" iff its start is
//! in the starts-set AND its end is in the ends-set, even if those two
//! boundaries never co-occur on the same reference transcript. This is
//! upstream's actual (approximate) definition of "known", not a bug to
//! fix.
//!
//! **Two-pass design preserved**: the read-processing pass counts every
//! splice *event* (i.e. with duplicates across reads) into `event_*`
//! counters for the first pie chart; a second pass over the *distinct*
//! junctions (in first-seen order) recomputes known/novel counts fresh
//! for the `.xls` rows and the second pie chart. These are genuinely
//! different numbers, not a bug.
//!
//! SAM-text and CRAM input (DIV-0002/0004) are supported at the CLI layer via
//! `rseqc_formats::open_alignments`; this module's own functions were
//! unaffected (already generic over
//! `IntoIterator<Item = io::Result<bam::Record>>`).

use std::collections::{HashMap, HashSet};
use std::io::{self, BufRead};

use noodles_bam as bam;
use noodles_sam as sam;
use rseqc_formats::cigar::fetch_intron_blocks;

pub struct JunctionModel {
    ref_intron_starts: HashMap<String, HashSet<i64>>,
    ref_intron_ends: HashMap<String, HashSet<i64>>,
}

/// Parses the reference BED12 gene model, matching upstream's specific
/// (partial) defensiveness: comment/track/browser lines are skipped;
/// lines with fewer than 12 whitespace-separated fields are skipped with
/// a stderr note; anything else that fails to parse as an integer is an
/// uncaught error (propagated as `io::Error`), matching upstream's bare
/// `int(...)` calls with no surrounding try/except. Upstream's `if
/// int(fields[9] == 1): continue` is DEAD CODE (compares the string
/// `fields[9]` to the int `1`, always False) and is correctly never
/// replicated here -- there is no "skip single-exon transcripts" behavior
/// to port, despite the comment implying one was intended.
pub fn build_model(reader: impl BufRead) -> io::Result<JunctionModel> {
    let mut ref_intron_starts: HashMap<String, HashSet<i64>> = HashMap::new();
    let mut ref_intron_ends: HashMap<String, HashSet<i64>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 12 {
            continue;
        }

        let bad = |msg: &str| io::Error::new(io::ErrorKind::InvalidData, msg.to_string());
        let chrom = fields[0].to_uppercase();
        let tx_start: i64 = fields[1].parse().map_err(|_| bad("invalid txStart"))?;
        let block_sizes: Vec<i64> = fields[10]
            .trim_end_matches(',')
            .split(',')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse().map_err(|_| bad("invalid blockSizes")))
            .collect::<io::Result<_>>()?;
        let relative_starts: Vec<i64> = fields[11]
            .trim_end_matches(',')
            .split(',')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse().map_err(|_| bad("invalid blockStarts")))
            .collect::<io::Result<_>>()?;

        let exon_starts: Vec<i64> = relative_starts.iter().map(|&rs| tx_start + rs).collect();
        let exon_ends: Vec<i64> = exon_starts.iter().zip(block_sizes.iter()).map(|(&s, &sz)| s + sz).collect();

        if exon_starts.len() < 2 {
            continue; // no introns possible with 0 or 1 exon
        }
        let intron_starts = &exon_ends[..exon_ends.len() - 1];
        let intron_ends = &exon_starts[1..];
        for (&s, &e) in intron_starts.iter().zip(intron_ends.iter()) {
            ref_intron_starts.entry(chrom.clone()).or_default().insert(s);
            ref_intron_ends.entry(chrom.clone()).or_default().insert(e);
        }
    }

    Ok(JunctionModel { ref_intron_starts, ref_intron_ends })
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JunctionRow {
    /// Upstream's literal `chrom.replace("CHR","chr")` (substring
    /// replace, not a case-insensitive prefix strip).
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub count: u64,
    pub annotation: &'static str,
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct JunctionResult {
    pub event_total: u64,
    pub event_known: u64,
    pub event_novel35: u64,
    pub event_novel3or5: u64,
    pub filtered: u64,
    pub rows: Vec<JunctionRow>,
    pub junction_known: u64,
    pub junction_novel35: u64,
    pub junction_novel3or5: u64,
}

fn is_known(model: &JunctionModel, chrom: &str, start: i64, end: i64) -> (bool, bool) {
    let is_start = model.ref_intron_starts.get(chrom).is_some_and(|s| s.contains(&start));
    let is_end = model.ref_intron_ends.get(chrom).is_some_and(|s| s.contains(&end));
    (is_start, is_end)
}

pub fn compute_junction_annotation<I>(
    records: I,
    header: &sam::Header,
    model: &JunctionModel,
    q_cut: u8,
    min_intron: i64,
) -> io::Result<JunctionResult>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut event_total = 0u64;
    let mut event_known = 0u64;
    let mut event_novel35 = 0u64;
    let mut event_novel3or5 = 0u64;
    let mut filtered = 0u64;

    let mut order: Vec<(String, i64, i64)> = Vec::new();
    let mut counts: HashMap<(String, i64, i64), u64> = HashMap::new();

    for result in records {
        let record = result?;
        let flags = record.flags();
        if flags.is_qc_fail() || flags.is_duplicate() || flags.is_secondary() || flags.is_unmapped() {
            continue;
        }
        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        if mapq < q_cut {
            continue;
        }

        let Some(ref_id) = record.reference_sequence_id().transpose()? else { continue };
        let Some((chrom_bstr, _)) = header.reference_sequences().get_index(ref_id) else { continue };
        let chrom = chrom_bstr.to_string().to_uppercase();

        let Some(pos) = record.alignment_start().transpose()? else { continue };
        let hit_st = (pos.get() - 1) as i64;

        let ops: Vec<_> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;
        let intron_blocks = fetch_intron_blocks(hit_st as usize, ops.iter().copied());
        if intron_blocks.is_empty() {
            continue;
        }

        for (s, e) in intron_blocks {
            let (s, e) = (s as i64, e as i64);
            event_total += 1;
            if e - s < min_intron {
                filtered += 1;
                continue;
            }

            let key = (chrom.clone(), s, e);
            if !counts.contains_key(&key) {
                order.push(key.clone());
            }
            *counts.entry(key).or_insert(0) += 1;

            let (is_start, is_end) = is_known(model, &chrom, s, e);
            if is_start && is_end {
                event_known += 1;
            } else if !is_start && !is_end {
                event_novel35 += 1;
            } else {
                event_novel3or5 += 1;
            }
        }
    }

    let mut rows = Vec::with_capacity(order.len());
    let mut junction_known = 0u64;
    let mut junction_novel35 = 0u64;
    let mut junction_novel3or5 = 0u64;

    for key in &order {
        let (chrom, s, e) = key;
        let count = counts[key];
        let (is_start, is_end) = is_known(model, chrom, *s, *e);
        let annotation = if is_start && is_end {
            junction_known += 1;
            "annotated"
        } else if !is_start && !is_end {
            junction_novel35 += 1;
            "complete_novel"
        } else {
            junction_novel3or5 += 1;
            "partial_novel"
        };
        rows.push(JunctionRow {
            chrom: chrom.replace("CHR", "chr"),
            start: *s,
            end: *e,
            count,
            annotation,
        });
    }

    Ok(JunctionResult {
        event_total,
        event_known,
        event_novel35,
        event_novel3or5,
        filtered,
        rows,
        junction_known,
        junction_novel35,
        junction_novel3or5,
    })
}

pub fn render_xls(result: &JunctionResult) -> String {
    let mut lines = vec!["chrom\tintron_st(0-based)\tintron_end(1-based)\tread_count\tannotation".to_string()];
    for row in &result.rows {
        // Exact upstream quirk: a literal space between the trailing tab
        // and the annotation word (two separate print() calls upstream,
        // the first ending with `end=' '` instead of a newline).
        lines.push(format!("{}\t{}\t{}\t{}\t {}", row.chrom, row.start, row.end, row.count, row.annotation));
    }
    // Trailing newline: the second print() per row (the annotation
    // word) uses the default end='\n', including for the last row.
    format!("{}\n", lines.join("\n"))
}

/// Python 3's `round()`: round-half-to-even, NOT round-half-away-from-zero.
/// Only used for the pie-chart label percentages (cosmetic `%d%%` text);
/// the `events=c(...)`/`junction=c(...)` data lines use full float
/// precision via ordinary `to_string()`, unaffected by this.
fn python_round(x: f64) -> i64 {
    let floor = x.floor();
    let diff = x - floor;
    if diff < 0.5 {
        floor as i64
    } else if diff > 0.5 {
        floor as i64 + 1
    } else if (floor as i64) % 2 == 0 {
        floor as i64
    } else {
        floor as i64 + 1
    }
}

pub fn render_r_script(result: &JunctionResult, out_prefix: &str) -> String {
    let mut lines = Vec::new();

    // Upstream's pie-1 denominator is `total_junc`, the RAW event total
    // that was incremented before the min_intron filter check runs --
    // i.e. it includes filtered introns too, never subtracting
    // `filtered_junc`. Use `event_total` as-is, unmodified.
    let denom1 = result.event_total as f64;
    let events = [result.event_novel3or5, result.event_novel35, result.event_known];
    let events_csv: Vec<String> = events.iter().map(|&v| (v as f64 * 100.0 / denom1).to_string()).collect();
    let events_pct: Vec<i64> = events.iter().map(|&v| python_round(v as f64 * 100.0 / denom1)).collect();

    lines.push(format!("pdf(\"{out_prefix}.splice_events.pdf\")"));
    lines.push(format!("events=c({})", events_csv.join(",")));
    lines.push(format!(
        "pie(events,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main=\"splicing events\",labels=c(\"partial_novel {}%\",\"complete_novel {}%\",\"known {}%\"))",
        events_pct[0], events_pct[1], events_pct[2]
    ));
    lines.push("dev.off()".to_string());
    lines.push(String::new());

    let denom2 = (result.junction_known + result.junction_novel35 + result.junction_novel3or5) as f64;
    let junctions = [result.junction_novel3or5, result.junction_novel35, result.junction_known];
    let junction_csv: Vec<String> = junctions.iter().map(|&v| (v as f64 * 100.0 / denom2).to_string()).collect();
    let junction_pct: Vec<i64> = junctions.iter().map(|&v| python_round(v as f64 * 100.0 / denom2)).collect();

    lines.push(format!("pdf(\"{out_prefix}.splice_junction.pdf\")"));
    lines.push(format!("junction=c({})", junction_csv.join(",")));
    lines.push(format!(
        "pie(junction,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main=\"splicing junctions\",labels=c(\"partial_novel {}%\",\"complete_novel {}%\",\"known {}%\"))",
        junction_pct[0], junction_pct[1], junction_pct[2]
    ));
    lines.push("dev.off()".to_string());

    // Trailing newline: upstream's plain `print(...)` calls each add
    // their own trailing newline, including the final `dev.off()`.
    format!("{}\n", lines.join("\n"))
}

fn junction_color(annotation: &str) -> &'static str {
    match annotation {
        "annotated" => "205,0,0",
        "partial_novel" => "0,205,0",
        "complete_novel" => "0,0,205",
        _ => "0,0,0",
    }
}

/// Ports `generate_bed12` (default `size=1`), operating directly on the
/// already-computed rows rather than re-parsing the `.xls` text.
pub fn render_bed12(result: &JunctionResult, size: i64) -> String {
    let mut lines = Vec::new();
    for row in &result.rows {
        let start = row.start - size;
        let end = row.end + size;
        let block_sizes = format!("{size},{size}");
        let block_starts = format!("0,{}", end - size - start);
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}\t2\t{}\t{}",
            row.chrom,
            start,
            end,
            row.annotation,
            row.count,
            start,
            end,
            junction_color(row.annotation),
            block_sizes,
            block_starts,
        ));
    }
    // Trailing newline: upstream's plain `print(...)` call per row adds
    // its own trailing newline, including for the last row.
    format!("{}\n", lines.join("\n"))
}

/// Ports `generate_interact` (default `size=1`). `bam_file` is the input
/// BAM path as displayed in the track description line.
pub fn render_interact(result: &JunctionResult, bam_file: &str, size: i64) -> String {
    let mut lines = Vec::new();
    lines.push(format!(
        "track type=interact name=\"Splice junctions\" description=\"Splice junctions detected from {bam_file}\" maxHeightPixels=200:200:50 visibility=full"
    ));
    for row in &result.rows {
        let chrom_start = row.start - size;
        let chrom_end = row.end + size;
        let name = format!("{}:{}-{}_{}", row.chrom, chrom_start, chrom_end, row.annotation);
        // Upstream: `value = float(score)`, printed via Python's
        // `str(float)` which always shows a decimal point (e.g. "5.0"),
        // unlike Rust's default f64 Display ("5"). row.count is always a
        // whole number here, so a fixed ".0" suffix matches exactly.
        let value = format!("{}.0", row.count);
        let source_start = chrom_start;
        let source_end = chrom_start + size;
        let source_name = format!("{}:{}-{}", row.chrom, source_start, source_end);
        let target_start = chrom_end - size;
        let target_end = chrom_end;
        let target_name = format!("{}:{}-{}", row.chrom, target_start, target_end);
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}\t{}\t.",
            row.chrom,
            chrom_start,
            chrom_end,
            name,
            row.count,
            value,
            "RNAseq_junction",
            junction_color(row.annotation),
            row.chrom,
            source_start,
            source_end,
            source_name,
            row.chrom,
            target_start,
            target_end,
            target_name,
        ));
    }
    // Trailing newline: upstream's plain `print(...)` call per row adds
    // its own trailing newline, including for the last row (and the
    // leading `track type=interact ...` line above it).
    format!("{}\n", lines.join("\n"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_model_skips_short_lines_and_dead_branch_never_skips() {
        // blockCount field (index 9) is "1" here -- if the dead
        // `int(fields[9]==1)` branch were (incorrectly) implemented as a
        // real check, this single-exon transcript would be skipped and
        // contribute no intron boundaries. It must NOT be skipped, but a
        // single exon has no introns anyway, so verify via a 2-exon
        // record instead, and a too-short line is skipped.
        let text = "chr1\t0\t500\tgeneA\t0\t+\t0\t500\t0\t2\t100,100,\t0,300,\nshort line\n";
        let model = build_model(text.as_bytes()).unwrap();
        assert!(model.ref_intron_starts.get("CHR1").unwrap().contains(&100));
        // exon_starts=[0,300], exon_ends=[100,400]; intron = (exon_ends[0],
        // exon_starts[1]) = (100, 300).
        assert!(model.ref_intron_ends.get("CHR1").unwrap().contains(&300));
    }

    #[test]
    fn render_xls_exact_row_format() {
        let result = JunctionResult {
            rows: vec![JunctionRow { chrom: "chr1".into(), start: 100, end: 200, count: 5, annotation: "annotated" }],
            ..Default::default()
        };
        let text = render_xls(&result);
        assert_eq!(
            text,
            "chrom\tintron_st(0-based)\tintron_end(1-based)\tread_count\tannotation\nchr1\t100\t200\t5\t annotated\n"
        );
    }

    #[test]
    fn python_round_half_to_even() {
        assert_eq!(python_round(12.5), 12);
        assert_eq!(python_round(13.5), 14);
        assert_eq!(python_round(20.4), 20);
        assert_eq!(python_round(20.6), 21);
    }

    #[test]
    fn render_r_script_exact_text() {
        // Independently derived via python3 -c from oracle/upstream-src/
        // src/qcmodule/SAM.py lines 3839-3886.
        let result = JunctionResult {
            event_total: 10,
            event_novel3or5: 2,
            event_novel35: 3,
            event_known: 5,
            junction_novel3or5: 1,
            junction_novel35: 2,
            junction_known: 5,
            ..Default::default()
        };
        let output = render_r_script(&result, "testprefix");
        let expected = "pdf(\"testprefix.splice_events.pdf\")\n\
events=c(20,30,50)\n\
pie(events,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main=\"splicing events\",labels=c(\"partial_novel 20%\",\"complete_novel 30%\",\"known 50%\"))\n\
dev.off()\n\
\n\
pdf(\"testprefix.splice_junction.pdf\")\n\
junction=c(12.5,25,62.5)\n\
pie(junction,col=c(2,3,4),init.angle=30,angle=c(60,120,150),density=c(70,70,70),main=\"splicing junctions\",labels=c(\"partial_novel 12%\",\"complete_novel 25%\",\"known 62%\"))\n\
dev.off()\n";
        assert_eq!(output, expected);
    }

    #[test]
    fn render_bed12_and_interact_smoke() {
        let result = JunctionResult {
            rows: vec![JunctionRow { chrom: "chr1".into(), start: 100, end: 200, count: 5, annotation: "annotated" }],
            ..Default::default()
        };
        let bed = render_bed12(&result, 1);
        // start=100-1=99, end=200+1=201; block_starts second value =
        // end-size-start = 201-1-99 = 101.
        assert_eq!(bed, "chr1\t99\t201\tannotated\t5\t.\t99\t201\t205,0,0\t2\t1,1\t0,101\n");

        // Exact text, independently verified against a real Python
        // execution of generate_interact's literal expressions.
        let interact = render_interact(&result, "sample.bam", 1);
        let expected = "track type=interact name=\"Splice junctions\" description=\"Splice junctions detected from sample.bam\" maxHeightPixels=200:200:50 visibility=full\n\
chr1\t99\t201\tchr1:99-201_annotated\t5\t5.0\tRNAseq_junction\t205,0,0\tchr1\t99\t100\tchr1:99-100\t.\tchr1\t200\t201\tchr1:200-201\t.\n";
        assert_eq!(interact, expected);
    }
}
