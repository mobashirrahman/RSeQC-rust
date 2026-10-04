//! Port of `geneBody_coverage2.py`: gene-body coverage from a single
//! BigWig signal file (the "deprecated" BigWig-input sibling of
//! `geneBody_coverage.py`, which takes BAM input). Ported from
//! `oracle/upstream-src/scripts/geneBody_coverage2.py`.
//!
//! Despite the similar name and purpose, this is NOT a thin variant of
//! `genebody_coverage.rs` -- its transcript-filtering and strand
//! handling use genuinely different mechanisms (see below), so nothing
//! is shared between the two modules except `percentile_list` (reused
//! directly, same `mystat.percentile_list` function upstream).
//!
//! **Preserves a real upstream bug, do not "fix"**: the "skip
//! transcripts shorter than 100bp" check is misplaced INSIDE the
//! per-exon accumulation loop (`for start,end in zip(...): gene_all_bases
//! .extend(...); if len(gene_all_bases) < 100: skip_gene = True; break`)
//! rather than after it (contrast `geneBody_coverage.py`'s sibling
//! function, which correctly checks the TOTAL after the full loop).
//! Because the running total only ever grows, this check can only ever
//! be false-then-stay-false or true-and-break -- meaning a multi-exon
//! transcript is skipped entirely unless its FIRST exon alone already
//! has >=100 bases, even if its full spliced length is well over 100bp.
//! Most real multi-exon transcripts (small individual exons) fail this
//! filter regardless of total mRNA length. Documented as DIV-0014.
//!
//! **Preserves a second, subtler quirk**: strand handling is NOT "compute
//! ascending percentiles then reverse the 100-value output" (as
//! `geneBody_coverage.py` does) -- it's "sort the RAW genomic positions
//! descending BEFORE calling `percentile_list`" for `strand == "-"`.
//! These are NOT equivalent due to `percentile_list`'s asymmetric
//! floor/ceil interpolation formula; replicated via the exact same
//! sort-then-percentile mechanism, not the reverse-the-output shortcut.
use std::collections::{HashMap, HashSet};
use std::io::{self, BufRead};

use rseqc_formats::bigwig::BigWigReader;

use crate::genebody_coverage::percentile_list;
use crate::python_fmt::python_str_float;

#[derive(Debug, Clone)]
struct Bed12Record {
    chrom: String,
    strand: String,
    exon_starts: Vec<i64>,
    exon_ends: Vec<i64>,
}

/// Ports `parse_bed12_line`: hard-fails (returns `None`, caller prints
/// a warning and skips the line) on fewer than 12 columns or any
/// non-integer coordinate field.
fn parse_bed12_line(line: &str) -> Option<Bed12Record> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() < 12 {
        return None;
    }
    let chrom = fields[0].to_string();
    let tx_start: i64 = fields[1].parse().ok()?;
    let strand = fields[5].to_string();

    let exon_starts: Vec<i64> = fields[11].trim_end_matches(',').split(',').filter(|s| !s.is_empty()).map(|s| s.parse::<i64>().ok().map(|v| v + tx_start)).collect::<Option<_>>()?;
    let exon_sizes: Vec<i64> = fields[10].trim_end_matches(',').split(',').filter(|s| !s.is_empty()).map(|s| s.parse::<i64>().ok()).collect::<Option<_>>()?;
    let exon_ends: Vec<i64> = exon_starts.iter().zip(exon_sizes.iter()).map(|(&s, &sz)| s + sz).collect();

    Some(Bed12Record { chrom, strand, exon_starts, exon_ends })
}

/// Runs the full per-transcript coverage-sampling pass. Ports
/// `coverage_gene_body_bigwig`. Returns `(coverage_by_index,
/// gene_count)`; `coverage_by_index` is empty if no transcript ever
/// passed the (buggy, see module docs) length filter, matching
/// upstream's `defaultdict(float)` staying empty in that case (NOT a
/// 100-length all-zero vector).
pub fn coverage_gene_body_bigwig(bw: &mut BigWigReader, refbed: impl BufRead) -> io::Result<(Vec<f64>, i64)> {
    let chrom_lengths: HashMap<String, u64> =
        bw.chroms().into_iter().map(|(n, l)| (n, l as u64)).collect();
    let chrom_set: HashSet<String> = chrom_lengths.keys().cloned().collect();

    let mut coverage: Vec<f64> = Vec::new();
    let mut gene_count = 0i64;

    for (line_number, line) in refbed.lines().enumerate() {
        let line = line?;
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }

        let Some(rec) = parse_bed12_line(&line) else {
            eprintln!("[NOTE: input BED must be 12-column] skipped line {}: {}", line_number + 1, line.trim_end());
            continue;
        };
        if !chrom_set.contains(&rec.chrom) {
            continue;
        }
        gene_count += 1;

        let mut gene_all_bases: Vec<i64> = Vec::new();
        let mut skip_gene = false;
        for (&s, &e) in rec.exon_starts.iter().zip(rec.exon_ends.iter()) {
            for p in (s + 1)..=e {
                gene_all_bases.push(p);
            }
            if (gene_all_bases.len() as i64) < 100 {
                skip_gene = true;
                break;
            }
        }
        if skip_gene {
            continue;
        }

        if rec.strand == "-" {
            gene_all_bases.sort_unstable_by(|a, b| b.cmp(a));
        } else {
            gene_all_bases.sort_unstable();
        }
        let percentile_bases = percentile_list(&gene_all_bases);

        if coverage.len() < percentile_bases.len() {
            coverage.resize(percentile_bases.len(), 0.0);
        }
        for (index, &genomic_position) in percentile_bases.iter().enumerate() {
            // A sampled position past the end of its chromosome is REFUSED, not
            // silently treated as zero.
            //
            // Upstream calls pyBigWig's `values(chrom, pos-1, pos)`, which raises
            // "Invalid interval bounds!" and exits 1 when that interval lies outside
            // the chromosome. This port's BigWig reader returns NaN for such an
            // interval instead of failing, and the NaN was mapped to 0.0 -- so a gene
            // model extending past the last base of a chromosome produced a
            // confidently wrong coverage curve, with the uncovered tail silently
            // contributing nothing rather than the run refusing. That is the audit's
            // "no silent metric loss and no successful corrupt output" failure: the
            // answer looks like an answer and understates coverage for every affected
            // transcript.
            //
            // It was found by verification/verify_stream_declarations.py, which runs
            // each declared command once and noticed that this one SUCCEEDED where
            // upstream refused -- the mirror image of the three declaration bugs it
            // found, and the reason that tool reports per-command outcomes at all.
            if let Some(&chrom_len) = chrom_lengths.get(&rec.chrom) {
                if genomic_position < 1 || genomic_position as u64 > chrom_len {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "gene model position {} on {} lies outside that \
                             chromosome in the BigWig ({} bp). Upstream refuses this \
                             input; it is refused here too rather than counted as zero \
                             coverage, which would silently understate this transcript.",
                            genomic_position, rec.chrom, chrom_len
                        ),
                    ));
                }
            }
            let signal = bw.values(&rec.chrom, (genomic_position - 1) as u32, genomic_position as u32)?;
            let v = signal[0];
            coverage[index] += if v.is_nan() { 0.0 } else { v as f64 };
        }

        // Upstream: `print(f"\t{gene_count} genes finished\r", end=' ')`
        // -- the literal's own `\r` plus a trailing space from `end`.
        eprint!("\t{gene_count} genes finished\r ");
    }
    eprintln!();

    Ok((coverage, gene_count))
}

/// Renders the `.geneBodyCoverage.txt` table: `percentile\tcount`
/// header, then one `index\tvalue` row per populated index (NOT
/// necessarily 100 rows -- see module docs).
pub fn render_coverage_txt(coverage: &[f64]) -> String {
    let mut out = String::from("percentile\tcount\n");
    for (index, value) in coverage.iter().enumerate() {
        out.push_str(&format!("{index}\t{}\n", python_str_float(*value)));
    }
    out
}

/// Renders the single-vector step-plot R script. Ports the R-writing
/// half of `coverage_gene_body_bigwig`.
pub fn render_r_script(coverage: &[f64], gene_count: i64, output_prefix: &str, graph_type: &str) -> String {
    let plot_file = format!("{output_prefix}.geneBodyCoverage.{graph_type}");
    let y_values: Vec<String> = coverage.iter().map(|&v| python_str_float(v)).collect();
    format!(
        "{graph_type}('{plot_file}')\nx=1:100\ny=c({})\nplot(x, y/{gene_count}, xlab=\"percentile of gene body (5'->3')\", ylab='average wigsum', type='s')\ndev.off()\n",
        y_values.join(",")
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn parse_bed12_line_extracts_fields() {
        let line = "chr1\t0\t300\ttx1\t0\t+\t0\t300\t0\t2\t100,100,\t0,200,";
        let rec = parse_bed12_line(line).unwrap();
        assert_eq!(rec.chrom, "chr1");
        assert_eq!(rec.strand, "+");
        assert_eq!(rec.exon_starts, vec![0, 200]);
        assert_eq!(rec.exon_ends, vec![100, 300]);
    }

    #[test]
    fn parse_bed12_line_rejects_short_lines() {
        assert!(parse_bed12_line("too short").is_none());
    }

    #[test]
    fn buggy_early_break_skips_multiexon_transcript_with_small_first_exon() {
        // Two-exon transcript, 50bp + 50bp = 100bp total (would clear a
        // POST-loop check), but the first exon alone is only 50bp: the
        // upstream bug breaks out after exon 1 before ever reaching
        // exon 2, so this transcript must be skipped entirely.
        let bed = "chr1\t0\t200\ttx1\t0\t+\t0\t200\t0\t2\t50,50,\t0,100,\n";
        // Single-exon transcript, exactly 150bp: passes on exon 1 alone.
        let bed_ok = "chr1\t0\t150\ttx2\t0\t+\t0\t150\t0\t1\t150,\t0,\n";

        // We can't easily construct a real BigWigReader in a unit test
        // (no write support pulled in -- see crates/formats/src/
        // bigwig.rs module docs), so this test exercises the pure
        // gene_all_bases/skip_gene logic directly via a hand-copy of the
        // loop body instead of the full coverage_gene_body_bigwig path
        // (that path is exercised by the CLI smoke test instead).
        fn would_skip(bed_line: &str) -> bool {
            let rec = parse_bed12_line(bed_line).unwrap();
            let mut gene_all_bases: Vec<i64> = Vec::new();
            let mut skip_gene = false;
            for (&s, &e) in rec.exon_starts.iter().zip(rec.exon_ends.iter()) {
                for p in (s + 1)..=e {
                    gene_all_bases.push(p);
                }
                if (gene_all_bases.len() as i64) < 100 {
                    skip_gene = true;
                    break;
                }
            }
            skip_gene
        }

        assert!(would_skip(bed.trim_end()));
        assert!(!would_skip(bed_ok.trim_end()));
    }

    #[test]
    fn render_coverage_txt_matches_python_str_float_formatting() {
        let coverage = vec![0.0, 1.5, 3.0];
        let txt = render_coverage_txt(&coverage);
        assert_eq!(txt, "percentile\tcount\n0\t0.0\n1\t1.5\n2\t3.0\n");
    }

    #[test]
    fn render_coverage_txt_empty_when_no_transcript_passed() {
        assert_eq!(render_coverage_txt(&[]), "percentile\tcount\n");
    }

    #[test]
    fn render_r_script_exact_text() {
        let coverage = vec![1.0, 2.0];
        let script = render_r_script(&coverage, 5, "out", "png");
        assert_eq!(script, "png('out.geneBodyCoverage.png')\nx=1:100\ny=c(1.0,2.0)\nplot(x, y/5, xlab=\"percentile of gene body (5'->3')\", ylab='average wigsum', type='s')\ndev.off()\n");
    }

    #[test]
    fn a_model_position_past_the_chromosome_end_is_refused_not_zero_filled() {
        // The silent-metric-loss case this guards.
        //
        // Upstream calls pyBigWig's values(chrom, pos-1, pos), which raises
        // "Invalid interval bounds!" and exits 1 when the interval is outside the
        // chromosome. The port's BigWig reader returns NaN there instead of failing,
        // and NaN was mapped to 0.0 -- so a model extending past a chromosome's last
        // base produced a confident coverage curve whose uncovered tail contributed
        // nothing instead of the run refusing. Found by verify_stream_declarations.py
        // noticing this command SUCCEEDING where upstream refused.
        //
        // (A `"""` docstring cannot be used here: Rust reads it as `""` followed by a
        // string, so the embedded quotes below it would be a syntax error.)
        let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../formats/tests/fixtures/pybigwig_test.bw");
        let mut bw = BigWigReader::open(&fixture).unwrap();
        // Chromosome "1" is 195,471,971 bp in the fixture, so a 100 bp single-exon
        // transcript at 195,471,950 puts sampled positions past the end.
        let bed = "1\t195471950\t195472050\tOUTOFRANGE\t0\t+\t195471950\t195472050\t0\t1\t100,\t0,\n";
        let err = coverage_gene_body_bigwig(&mut bw, Cursor::new(bed))
            .expect_err("a position past the chromosome end must be refused, not zero-filled");
        let text = err.to_string();
        assert!(text.contains("outside that chromosome"), "message: {text}");
        assert!(text.contains("195471971"), "message must name the length: {text}");
        // The reported position is the first SAMPLED position past the end, not the
        // transcript start, because the check runs per percentile sample. So the
        // invariant is "the position it names is beyond the chromosome", not "it names
        // the position I happened to write in the fixture".
        let reported: u64 = text
            .split_whitespace()
            .nth(3)
            .and_then(|w| w.parse().ok())
            .unwrap_or_else(|| panic!("message must name a numeric position: {text}"));
        assert!(
            reported > 195_471_971,
            "the position reported as out of range ({reported}) must exceed the \
             chromosome length the message also states"
        );
    }

    #[test]
    fn a_model_inside_the_chromosome_is_still_accepted() {
        // The refusal must not turn into a blanket failure: the same fixture, a
        // transcript well inside chromosome "1", has to keep working.
        let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../formats/tests/fixtures/pybigwig_test.bw");
        let mut bw = BigWigReader::open(&fixture).unwrap();
        let bed = "1\t1000\t1200\tINSIDE\t0\t+\t1000\t1200\t0\t1\t200,\t0,\n";
        let (coverage, gene_count) =
            coverage_gene_body_bigwig(&mut bw, Cursor::new(bed)).unwrap();
        assert_eq!(gene_count, 1);
        assert_eq!(coverage.len(), 100);
    }

    #[test]
    fn coverage_gene_body_bigwig_returns_empty_for_no_matching_chrom() {
        // Exercise the actual function with a real (fixture) BigWig
        // reader for the "no BED lines match any BigWig chromosome"
        // path, which needs no real signal data to verify.
        let fixture = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../formats/tests/fixtures/pybigwig_test.bw");
        let mut bw = BigWigReader::open(&fixture).unwrap();
        let bed = "chrNoSuchChrom\t0\t300\ttx1\t0\t+\t0\t300\t0\t1\t300,\t0,\n";
        let (coverage, gene_count) = coverage_gene_body_bigwig(&mut bw, Cursor::new(bed)).unwrap();
        assert_eq!(gene_count, 0);
        assert!(coverage.is_empty());
    }
}
