//! Port of `read_hexamer.py`: calculate normalized hexamer frequencies
//! from FASTA/FASTQ files. Ported from `oracle/upstream-src/scripts/
//! read_hexamer.py` and the `seq_generator`/`word_generator`/
//! `all_possible_kmer`/`kmer_freq_file` functions of `oracle/
//! upstream-src/src/qcmodule/FrameKmer.py` (`kmer_ratio` is dead code
//! for this command, not ported).
//!
//! No BAM/noodles APIs and no new crate dependencies are needed: the
//! upstream sequence parser is a hand-rolled, deliberately permissive
//! line scanner, not a spec-compliant FASTA/FASTQ parser.
//!
//! **Preserves several deliberate upstream quirks, do not "fix"**:
//! - Sequence lines are matched against `^[ACGTN]+$` (after
//!   uppercasing); any other line -- a FASTQ `+` separator, a quality
//!   string, garbage -- is silently dropped, not an error. This means a
//!   FASTQ file is effectively read as if it were FASTA (quality lines
//!   just never match and vanish).
//! - The WHOLE file's sequences collapse into ONE cumulative k-mer count
//!   table (`kmer_freq_file`), not one table per sequence record.
//! - Kmers containing `N` are entirely excluded from the returned table
//!   (upstream enumerates all 5^6 `A,C,G,T,N` combinations then filters
//!   `N`-containing ones out; removing a digit from a lexicographic
//!   enumeration preserves the relative order of what remains, so this
//!   is implemented directly as a 4-symbol `A,C,G,T` enumeration).
use std::collections::{HashMap, HashSet};
use std::io::{self, BufRead};
use std::path::Path;

use crate::python_fmt::python_g12;

/// Yields `(name, sequence)` pairs from one file, replicating
/// `seq_generator`. Each line is trimmed and uppercased -- and, unlike a
/// typical "just normalize the sequence" reading, upstream REASSIGNS
/// `line` to that uppercased copy before doing anything else with it, so
/// the header name is extracted from the uppercased text too (`name =
/// line.split()[0][1:]` runs against the already-`.upper()`'d `line`,
/// there is no separately preserved original-case copy anywhere in
/// upstream). Confirmed via a live probe of the real installed
/// `qcmodule.FrameKmer.seq_generator` against a mixed-case header: it
/// returns `'MIXEDCASEHEADER'`, not `'MixedCaseHeader'`. Not currently
/// observable through `read_hexamer.py`'s own output (its only caller,
/// `kmer_freq_file`, discards the name and keeps only the sequence), but
/// ported exactly anyway since this is a shared, documented-as-exact
/// utility function, not a one-off formatter. A `#`-prefixed or empty
/// line is skipped. A `>`/`@`-prefixed line starts a new record: if a
/// sequence was already accumulating, it is flushed first; the new
/// record's `name` is the (uppercased) header line's first whitespace
/// token with its leading `>`/`@` stripped. Any other line is appended
/// to the current sequence only if every character is `A`/`C`/`G`/`T`/`N`.
/// The final accumulated `(name, sequence)` is always yielded once at
/// EOF, even if empty.
pub fn seq_generator(reader: impl BufRead) -> io::Result<Vec<(String, String)>> {
    let mut out = Vec::new();
    let mut name = String::new();
    let mut tmpseq = String::new();

    for line in reader.lines() {
        let raw = line?;
        let upper = raw.trim().to_uppercase();
        if upper.is_empty() || upper.starts_with('#') {
            continue;
        }
        if upper.starts_with('>') || upper.starts_with('@') {
            if !tmpseq.is_empty() {
                out.push((std::mem::take(&mut name), std::mem::take(&mut tmpseq)));
            }
            let first_token = upper.split_whitespace().next().unwrap_or(&upper);
            name = first_token.chars().skip(1).collect();
        } else if upper.chars().all(|c| matches!(c, 'A' | 'C' | 'G' | 'T' | 'N')) {
            tmpseq.push_str(&upper);
        }
    }
    out.push((name, tmpseq));

    Ok(out)
}

/// Yields `word_size`-length windows of `seq` starting at
/// `frame, frame+step_size, frame+2*step_size, ...`, dropping any
/// truncated trailing window. Ports `word_generator` exactly (Python:
/// `for i in range(frame, len(seq), step_size): word = seq[i:i+word_size];
/// if len(word) == word_size: yield word`).
pub fn word_generator(seq: &str, word_size: usize, step_size: usize, frame: usize) -> Vec<String> {
    let bytes = seq.as_bytes();
    let mut out = Vec::new();
    let mut i = frame;
    while i < bytes.len() {
        if i + word_size <= bytes.len() {
            out.push(String::from_utf8_lossy(&bytes[i..i + word_size]).into_owned());
        }
        i += step_size;
    }
    out
}

/// Enumerates all `4^word_size` combinations of exactly `A,C,G,T` (no
/// `N`) in the same relative order upstream's `itertools.product(['A',
/// 'C','G','T','N'], repeat=word_size)` produces once `N`-containing
/// tuples are filtered out (see module docs).
pub fn all_possible_kmer(word_size: usize) -> Vec<String> {
    let bases = *b"ACGT";
    let mut out = vec![String::new()];
    for _ in 0..word_size {
        let mut next = Vec::with_capacity(out.len() * 4);
        for prefix in &out {
            for &b in &bases {
                let mut s = prefix.clone();
                s.push(b as char);
                next.push(s);
            }
        }
        out = next;
    }
    out
}

/// Builds one cumulative 6-mer count table across every sequence in the
/// file, then returns exactly the 4096 `A,C,G,T`-only 6-mers (unseen
/// ones default to 0; any `N`-containing 6-mer is excluded entirely).
/// Ports `kmer_freq_file(fastafile, word_size=6, step_size=1, frame=0)`.
pub fn kmer_freq_file(reader: impl BufRead) -> io::Result<HashMap<String, i64>> {
    let mut count_table: HashMap<String, i64> = HashMap::new();
    for (_name, seq) in seq_generator(reader)? {
        for word in word_generator(&seq, 6, 1, 0) {
            *count_table.entry(word).or_insert(0) += 1;
        }
    }

    let mut result = HashMap::new();
    for kmer in all_possible_kmer(6) {
        result.insert(kmer.clone(), *count_table.get(&kmer).unwrap_or(&0));
    }
    Ok(result)
}

/// Creates a stable, unique display name for `path`, mutating `existing`
/// to record it. Ports `unique_display_name`: prefer the bare file name;
/// if taken, fall back to the full path string; if THAT is also taken,
/// append `#2`, `#3`, ... until unique.
pub fn unique_display_name(path: &Path, existing: &mut HashSet<String>) -> String {
    let file_name = path.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
    if !existing.contains(&file_name) {
        existing.insert(file_name.clone());
        return file_name;
    }

    let full = path.to_string_lossy().into_owned();
    if !existing.contains(&full) {
        existing.insert(full.clone());
        return full;
    }

    let mut suffix = 2u32;
    loop {
        let candidate = format!("{full}#{suffix}");
        if !existing.contains(&candidate) {
            existing.insert(candidate.clone());
            return candidate;
        }
        suffix += 1;
    }
}

/// Returns `counts.get(kmer, 0.0) / total`, or `0.0` if `total <= 0`.
/// Ports `normalized_frequency` exactly.
pub fn normalized_frequency(counts: &HashMap<String, i64>, total: f64, kmer: &str) -> f64 {
    if total <= 0.0 {
        return 0.0;
    }
    *counts.get(kmer).unwrap_or(&0) as f64 / total
}

/// Renders the `Hexamer\t<name>...` report: one row per non-`N` 6-mer
/// (in `all_possible_kmer` order), each column normalized against that
/// file's own total. Ports `write_report`.
pub fn render_report(names: &[String], tables: &HashMap<String, HashMap<String, i64>>, totals: &HashMap<String, f64>) -> String {
    let mut out = String::new();
    out.push_str("Hexamer\t");
    out.push_str(&names.join("\t"));
    out.push('\n');

    for kmer in all_possible_kmer(6) {
        out.push_str(&kmer);
        for name in names {
            out.push('\t');
            let value = normalized_frequency(&tables[name], totals[name], &kmer);
            out.push_str(&python_g12(value));
        }
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn seq_generator_parses_fasta_and_yields_final_record() {
        let text = ">seq1\nACGT\nACGT\n>seq2\nGGG\nNNN\n>seq3\n";
        let sequences = seq_generator(Cursor::new(text)).unwrap();
        assert_eq!(
            sequences,
            vec![("SEQ1".to_string(), "ACGTACGT".to_string()), ("SEQ2".to_string(), "GGGNNN".to_string()), ("SEQ3".to_string(), "".to_string())]
        );
    }

    #[test]
    fn seq_generator_uppercases_header_name_and_skips_comments() {
        // Upstream reassigns `line = line.strip().upper()` before
        // extracting the name -- confirmed via a live probe of the real
        // qcmodule.FrameKmer.seq_generator: a ">Seq1 description" header
        // comes back as "SEQ1", not "Seq1".
        let text = "# comment\n>Seq1 description\nacgt\n\n# another\n>seq2\nTTT\n";
        let sequences = seq_generator(Cursor::new(text)).unwrap();
        assert_eq!(sequences, vec![("SEQ1".to_string(), "ACGT".to_string()), ("SEQ2".to_string(), "TTT".to_string())]);
    }

    #[test]
    fn seq_generator_drops_non_acgtn_lines() {
        let text = ">seq1\nACGT\n+\n!@#$\nIIIIIII\n>seq2\nCCCC\n";
        let sequences = seq_generator(Cursor::new(text)).unwrap();
        assert_eq!(sequences, vec![("SEQ1".to_string(), "ACGT".to_string()), ("SEQ2".to_string(), "CCCC".to_string())]);
    }

    #[test]
    fn word_generator_matches_python_range_semantics() {
        // Cross-checked: for seq of length 12, word_size 6, step 1,
        // frame 0, positions 0..=6 (7 windows) are valid; i=7 onward
        // truncates below word_size and is dropped.
        let windows = word_generator("ACGTACGTACGT", 6, 1, 0);
        assert_eq!(windows, vec!["ACGTAC", "CGTACG", "GTACGT", "TACGTA", "ACGTAC", "CGTACG", "GTACGT"]);
    }

    #[test]
    fn word_generator_drops_truncated_trailing_window() {
        let windows = word_generator("ACGTAC", 6, 1, 0);
        assert_eq!(windows, vec!["ACGTAC"]);
    }

    #[test]
    fn all_possible_kmer_order_matches_itertools_product_after_n_filter() {
        let kmers = all_possible_kmer(2);
        assert_eq!(kmers, vec!["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]);
    }

    #[test]
    fn all_possible_kmer_is_4096_acgt_only_combinations() {
        let kmers = all_possible_kmer(6);
        assert_eq!(kmers.len(), 4096);
        assert!(kmers.iter().all(|k| k.len() == 6 && k.chars().all(|c| matches!(c, 'A' | 'C' | 'G' | 'T'))));
        assert_eq!(kmers[0], "AAAAAA");
        assert_eq!(*kmers.last().unwrap(), "TTTTTT");
    }

    #[test]
    fn kmer_freq_file_counts_across_whole_file_cumulatively() {
        // seq1="ACGTACGT" (8 chars) -> windows at i=0,1,2 -> ACGTAC, CGTACG, GTACGT
        // seq2="ACGTACGTACGT" (12 chars) -> 7 windows (see word_generator test)
        // Combined counts: ACGTAC: 1(seq1) + 2(seq2, i=0,4) = 3
        //                  CGTACG: 1(seq1) + 2(seq2, i=1,5) = 3
        //                  GTACGT: 1(seq1) + 2(seq2, i=2,6) = 3
        //                  TACGTA: 0(seq1) + 1(seq2, i=3)   = 1
        let text = ">seq1\nACGTACGT\n>seq2\nACGTACGTACGT\n";
        let counts = kmer_freq_file(Cursor::new(text)).unwrap();
        assert_eq!(counts["ACGTAC"], 3);
        assert_eq!(counts["CGTACG"], 3);
        assert_eq!(counts["GTACGT"], 3);
        assert_eq!(counts["TACGTA"], 1);
        assert_eq!(counts["AAAAAA"], 0);
    }

    #[test]
    fn unique_display_name_dedup_sequence() {
        let mut existing = HashSet::new();
        let n1 = unique_display_name(Path::new("dir1/file1.fa"), &mut existing);
        assert_eq!(n1, "file1.fa");
        let n2 = unique_display_name(Path::new("dir2/file1.fa"), &mut existing);
        assert_eq!(n2, "dir2/file1.fa");
        let n3 = unique_display_name(Path::new("dir2/file1.fa"), &mut existing);
        assert_eq!(n3, "dir2/file1.fa#2");
    }

    #[test]
    fn python_g12_reexport_smoke_check() {
        // Full python3-cross-checked coverage lives with the
        // implementation now, in python_fmt.rs::tests -- this just
        // confirms the re-export wires up correctly.
        assert_eq!(python_g12(1.0), "1");
        assert_eq!(python_g12(0.5), "0.5");
    }

    #[test]
    fn normalized_frequency_zero_total_is_zero() {
        let counts = HashMap::from([("AAAAAA".to_string(), 5)]);
        assert_eq!(normalized_frequency(&counts, 0.0, "AAAAAA"), 0.0);
        assert_eq!(normalized_frequency(&counts, 10.0, "AAAAAA"), 0.5);
        assert_eq!(normalized_frequency(&counts, 10.0, "CCCCCC"), 0.0);
    }

    #[test]
    fn render_report_header_and_row_shape() {
        let names = vec!["a.fa".to_string(), "b.fa".to_string()];
        let mut tables = HashMap::new();
        tables.insert("a.fa".to_string(), HashMap::from([("AAAAAA".to_string(), 1i64)]));
        tables.insert("b.fa".to_string(), HashMap::from([("AAAAAA".to_string(), 3i64)]));
        let totals = HashMap::from([("a.fa".to_string(), 4.0), ("b.fa".to_string(), 4.0)]);

        let report = render_report(&names, &tables, &totals);
        let mut lines = report.lines();
        assert_eq!(lines.next().unwrap(), "Hexamer\ta.fa\tb.fa");
        assert_eq!(lines.next().unwrap(), "AAAAAA\t0.25\t0.75");
        assert_eq!(report.lines().count(), 1 + 4096);
    }
}
