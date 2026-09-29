//! Port of `normalize_bigwig.py`: rescale a BigWig signal to a target
//! total WIG sum, writing variableStep WIG or bedGraph text. Ported
//! from `oracle/upstream-src/scripts/normalize_bigwig.py`.
//!
//! **Preserves a genuinely output-affecting upstream quirk, do not
//! "fix"**: `write_bedgraph_chromosome`'s contiguous-same-value run
//! detection (`value_to_positions`/`range_to_value`, built via
//! `itertools.groupby`) is scoped PER CHUNK, not per chromosome --
//! both the grouping state and the final sorted-print step live INSIDE
//! the per-chunk loop in upstream. A run of identical values that
//! happens to straddle a chunk boundary is therefore reported as TWO
//! separate bedGraph rows instead of one merged row. This makes
//! bedGraph output genuinely dependent on `--chunk`, unlike WIG output
//! (each position is printed independently there, so chunking never
//! changes WIG content). Replicated by chunking
//! `write_bedgraph_chromosome` for real (not just for memory,
//! cosmetically) rather than processing a whole chromosome in one pass.
//!
//! Simplifications with NO effect on output bytes: `signal_sum` always
//! calls `.values()` and sums non-NaN entries directly, skipping
//! upstream's separate `.stats()` pre-check call (mathematically
//! identical: an empty/no-data region sums to 0.0 either way). The
//! "does this chromosome have ANY data" check
//! (`bigwig.stats(chrom,0,size)[0] is None`) is replicated via "does
//! `intervals()` return anything", cheaper than materializing a
//! whole-chromosome values array just for an existence test.
use std::collections::HashMap;
use std::io::{self, BufRead};

use rseqc_formats::bigwig::BigWigReader;

use crate::python_fmt::python_str_float;

/// Ports `BED.tillingBed`: tiles `[0, chrom_size)` into `chunk_size`
/// windows, clamping the last one to `chrom_size`.
pub fn chromosome_chunks(chrom_size: i64, chunk_size: i64) -> Vec<(i64, i64)> {
    let mut out = Vec::new();
    let mut start = 0;
    while start < chrom_size {
        let end = (start + chunk_size).min(chrom_size);
        out.push((start, end));
        start += chunk_size;
    }
    out
}

/// Sums non-NaN signal over `[start, end)`. Ports `signal_sum`.
fn signal_sum(bw: &mut BigWigReader, chrom: &str, start: i64, end: i64) -> io::Result<f64> {
    let values = bw.values(chrom, start as u32, end as u32)?;
    Ok(values.iter().filter(|v| !v.is_nan()).map(|&v| v as f64).sum())
}

/// Sums signal over merged exon regions only. Ports
/// `calculate_exonic_wigsum`, including its 3 stderr progress lines
/// ("Extract exons from...", printed by the CLI beforehand since it
/// alone has the refgene path string; "Merge overlapping exons ...";
/// "Calculate WIG sum covered by {refgene_bed} only"). A merged-exon
/// chromosome name absent from the BigWig is skipped, matching
/// upstream's `except RuntimeError: continue` around pyBigWig's
/// invalid-chromosome error.
pub fn calculate_exonic_wigsum(bw: &mut BigWigReader, refbed: impl BufRead, refgene_path: &str) -> io::Result<f64> {
    let exons = rseqc_formats::bed::get_exon(refbed)?;
    eprintln!("Merge overlapping exons ...");
    let merged = rseqc_formats::interval::union_bed3(&exons);
    eprintln!("Calculate WIG sum covered by {refgene_path} only");
    let chrom_set: std::collections::HashSet<String> = bw.chroms().into_iter().map(|(n, _)| n).collect();

    let mut wig_sum = 0.0;
    for (chrom, start, end) in merged {
        if !chrom_set.contains(&chrom) {
            continue;
        }
        wig_sum += signal_sum(bw, &chrom, start, end)?;
    }
    Ok(wig_sum)
}

/// Sums signal across the whole genome, chunked. Ports
/// `calculate_genome_wigsum`.
pub fn calculate_genome_wigsum(bw: &mut BigWigReader, chrom_sizes: &[(String, i64)], chunk_size: i64) -> io::Result<f64> {
    let mut wig_sum = 0.0;
    for (chrom, size) in chrom_sizes {
        if bw.intervals(chrom, 0, *size as u32)?.is_empty() {
            eprintln!("Skip {chrom}!");
            continue;
        }
        eprintln!("Processing {chrom} ...");
        for (start, end) in chromosome_chunks(*size, chunk_size) {
            wig_sum += signal_sum(bw, chrom, start, end)?;
        }
    }
    Ok(wig_sum)
}

/// Writes one chromosome's `variableStep` WIG body: NaN treated as 0,
/// scaled by `weight`, zero-valued positions omitted. Ports
/// `write_wig_chromosome`.
pub fn write_wig_chromosome(out: &mut String, bw: &mut BigWigReader, chrom: &str, chrom_size: i64, chunk_size: i64, weight: f64) -> io::Result<()> {
    eprintln!("Writing {chrom} ...");
    out.push_str(&format!("variableStep chrom={chrom}\n"));
    for (start, end) in chromosome_chunks(chrom_size, chunk_size) {
        let values = bw.values(chrom, start as u32, end as u32)?;
        let mut coordinate = start;
        for &v in &values {
            coordinate += 1;
            let scaled = if v.is_nan() { 0.0 } else { v as f64 * weight };
            if scaled != 0.0 {
                out.push_str(&format!("{coordinate}\t{scaled:.2}\n"));
            }
        }
    }
    Ok(())
}

/// Groups one chunk's non-zero scaled values into contiguous
/// same-value runs and appends their bedGraph rows. Ports the per-chunk
/// body of `write_bedgraph_chromosome` (see module docs for why this
/// must NOT be merged across chunk boundaries).
fn write_bedgraph_chunk(out: &mut String, chrom: &str, chunk_start: i64, scaled_values: &[f64]) {
    let mut value_to_positions: HashMap<u64, (f64, Vec<i64>)> = HashMap::new();
    let mut coordinate = chunk_start;
    for &v in scaled_values {
        coordinate += 1;
        if v != 0.0 {
            value_to_positions.entry(v.to_bits()).or_insert_with(|| (v, Vec::new())).1.push(coordinate);
        }
    }

    let mut range_to_value: HashMap<i64, (i64, f64)> = HashMap::new();
    for (value, positions) in value_to_positions.values() {
        let mut i = 0;
        while i < positions.len() {
            let mut j = i;
            while j + 1 < positions.len() && positions[j + 1] == positions[j] + 1 {
                j += 1;
            }
            let run_start = positions[i];
            let length = (j - i + 1) as i64;
            range_to_value.insert(run_start - 1, (length, *value));
            i = j + 1;
        }
    }

    let mut starts: Vec<i64> = range_to_value.keys().copied().collect();
    starts.sort_unstable();
    for range_start in starts {
        let (length, value) = range_to_value[&range_start];
        out.push_str(&format!("{chrom}\t{range_start}\t{}\t{}\n", range_start + length, python_str_float(value)));
    }
}

/// Writes one chromosome's bedGraph body, chunk by chunk. Ports
/// `write_bedgraph_chromosome`.
pub fn write_bedgraph_chromosome(out: &mut String, bw: &mut BigWigReader, chrom: &str, chrom_size: i64, chunk_size: i64, weight: f64) -> io::Result<()> {
    eprintln!("Writing {chrom} ...");
    for (start, end) in chromosome_chunks(chrom_size, chunk_size) {
        let values = bw.values(chrom, start as u32, end as u32)?;
        let scaled: Vec<f64> = values.iter().map(|&v| if v.is_nan() { 0.0 } else { v as f64 * weight }).collect();
        write_bedgraph_chunk(out, chrom, start, &scaled);
    }
    Ok(())
}

#[derive(Debug)]
pub struct WigsumResult {
    pub chrom_sizes: Vec<(String, i64)>,
    pub observed_wigsum: f64,
    pub weight: f64,
}

/// Computes the observed WIG sum and normalization weight -- the FIRST
/// half of upstream's `normalize_bigwig`, up through the point it
/// prints "Total WIG sum is .../Normalization factor: ...". Split out
/// from the write phase (`render_normalized_body`) so the CLI can print
/// "Normalization factor: .../Normalizing BigWig file ..." at the
/// correct point in the sequence -- upstream's `main()` calls
/// `normalize_bigwig` as ONE function that does both phases back to
/// back with no opportunity to interleave caller-side prints; this
/// port's own `main`-equivalent (crates/cli/src/bin/normalize_bigwig.rs)
/// needs to emit the same lines, in the same order, without a String
/// round-trip through both phases first.
///
/// The "\nTotal WIG sum is {:.2}\n" diagnostic is printed here, not by
/// the caller, and unconditionally before the zero/negative-sum check --
/// matching upstream's `print(...)` call sitting textually BEFORE its
/// `if observed_wigsum <= 0: raise ValueError(...)`. Deferring that print
/// to the caller (only on `Ok`) used to swallow it whenever this
/// function errored, so the port's stderr silently dropped the "Total
/// WIG sum is <negative-or-zero-value>" line upstream always emits.
pub fn calculate_wigsum(bw: &mut BigWigReader, refgene: Option<impl BufRead>, total_wigsum: f64, chunk_size: i64, refgene_path: &str) -> io::Result<WigsumResult> {
    let chrom_sizes: Vec<(String, i64)> = bw.chroms().into_iter().map(|(n, l)| (n, l as i64)).collect();

    let observed_wigsum = match refgene {
        Some(r) => calculate_exonic_wigsum(bw, r, refgene_path)?,
        None => calculate_genome_wigsum(bw, &chrom_sizes, chunk_size)?,
    };

    eprintln!("\nTotal WIG sum is {observed_wigsum:.2}\n");

    if observed_wigsum <= 0.0 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "observed WIG sum is zero; normalization cannot be calculated"));
    }
    let weight = total_wigsum / observed_wigsum;

    Ok(WigsumResult { chrom_sizes, observed_wigsum, weight })
}

/// Writes every chromosome's normalized body. Ports the write loop at
/// the end of upstream's `normalize_bigwig`, called AFTER the CLI has
/// already printed "Total WIG sum/Normalization factor/Normalizing
/// BigWig file ..." (see `calculate_wigsum`'s doc comment).
pub fn render_normalized_body(bw: &mut BigWigReader, chrom_sizes: &[(String, i64)], chunk_size: i64, weight: f64, out_format: &str) -> io::Result<String> {
    let mut body = String::new();
    for (chrom, size) in chrom_sizes {
        if bw.intervals(chrom, 0, *size as u32)?.is_empty() {
            eprintln!("Skip {chrom}!");
            continue;
        }
        if out_format == "wig" {
            write_wig_chromosome(&mut body, bw, chrom, *size, chunk_size, weight)?;
        } else {
            write_bedgraph_chromosome(&mut body, bw, chrom, *size, chunk_size, weight)?;
        }
    }
    Ok(body)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chromosome_chunks_tiles_and_clamps_last_chunk() {
        assert_eq!(chromosome_chunks(1000, 300), vec![(0, 300), (300, 600), (600, 900), (900, 1000)]);
        assert_eq!(chromosome_chunks(300, 300), vec![(0, 300)]);
    }

    #[test]
    fn write_wig_chromosome_skips_zero_and_scales() {
        let mut out = String::new();
        let mut bw = fixture_reader();
        // chrom "1": values at [0,3)=0.1,0.2,0.3; gap (NaN) until 100;
        // [100,150)=1.4; [150,151)=1.5.
        write_wig_chromosome(&mut out, &mut bw, "1", 151, 1000, 2.0).unwrap();
        assert!(out.starts_with("variableStep chrom=1\n"));
        assert!(out.contains("1\t0.20\n")); // position 1 (1-based) = value[0]=0.1*2=0.2
        assert!(out.contains("3\t0.60\n")); // position 3 = value[2]=0.3*2=0.6
        assert!(!out.contains("\n4\t")); // NaN position must be omitted (scaled==0)
        assert!(out.contains("101\t2.80\n")); // position 101 = value[100]=1.4*2=2.8
        assert!(out.contains("151\t3.00\n")); // position 151 = value[150]=1.5*2=3.0
    }

    #[test]
    fn write_bedgraph_chunk_merges_contiguous_runs_of_one_value() {
        let mut out = String::new();
        // positions 1..=5 all value 2.0 (0-indexed values[0..5]).
        let values = vec![2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 3.0];
        write_bedgraph_chunk(&mut out, "chr1", 0, &values);
        // run of 2.0 at 1-based positions 1..5 -> 0-based bedGraph [0,5).
        assert!(out.contains("chr1\t0\t5\t2.0\n"));
        // single 3.0 at position 7 (1-based) -> 0-based [6,7).
        assert!(out.contains("chr1\t6\t7\t3.0\n"));
    }

    #[test]
    fn write_bedgraph_chunk_does_not_merge_across_a_chunk_boundary() {
        // Same logical run of value 5.0 split into two chunks by
        // calling write_bedgraph_chunk twice, matching upstream's
        // per-chunk-scoped grouping (see module docs).
        let mut out = String::new();
        write_bedgraph_chunk(&mut out, "chr1", 0, &[5.0, 5.0]); // positions 1,2
        write_bedgraph_chunk(&mut out, "chr1", 2, &[5.0, 5.0]); // positions 3,4 (contiguous with the above!)
        // Two separate rows, NOT one merged [0,4) row.
        assert!(out.contains("chr1\t0\t2\t5.0\n"));
        assert!(out.contains("chr1\t2\t4\t5.0\n"));
        assert!(!out.contains("chr1\t0\t4\t"));
    }

    #[test]
    fn calculate_genome_wigsum_matches_hand_summed_fixture_signal() {
        let mut bw = fixture_reader();
        let chrom_sizes: Vec<(String, i64)> = bw.chroms().into_iter().map(|(n, l)| (n, l as i64)).collect();
        // Only check chrom "10" (small enough to reason about exactly):
        // single interval (200,300,2.0) -> sum = 100*2.0 = 200.0.
        let sum10 = signal_sum(&mut bw, "10", 0, 130694993).unwrap();
        assert_eq!(sum10, 200.0);
        // Full genome sum should be >= chrom10's contribution alone.
        let total = calculate_genome_wigsum(&mut bw, &chrom_sizes, 500_000).unwrap();
        assert!(total >= 200.0);
    }

    #[test]
    fn calculate_wigsum_errors_on_non_positive_observed_sum() {
        // Regression for `normalize_bigwig_refgene_bgr_sig2`: upstream's
        // `if observed_wigsum <= 0: raise ValueError(...)`
        // (normalize_bigwig.py, in `normalize_bigwig()`) sits textually
        // AFTER its unconditional
        // `print(f"\nTotal WIG sum is {observed_wigsum:.2f}\n", file=sys.stderr)`.
        // A refgene region with no overlapping signal sums to exactly
        // 0.0, which must still hit the same `<= 0` branch (not just
        // `== 0`, and not silently succeed with weight = +inf).
        let mut bw = fixture_reader();
        // chrom "1" has intervals only within [0,151); [500,600) has no
        // data at all, so calculate_exonic_wigsum's signal_sum is 0.0.
        let refbed = "1\t500\t600\tno_signal\t0\t+\t500\t600\t0\t1\t100,\t0,\n";
        let err = calculate_wigsum(&mut bw, Some(refbed.as_bytes()), 100_000_000.0, 500_000, "test.bed12")
            .unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert_eq!(err.to_string(), "observed WIG sum is zero; normalization cannot be calculated");
    }

    fn fixture_reader() -> BigWigReader {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../formats/tests/fixtures/pybigwig_test.bw");
        BigWigReader::open(&path).unwrap()
    }
}
