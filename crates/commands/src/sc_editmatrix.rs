//! Port of `sc_editMatrix.py`: compare raw (sequencer-reported) vs
//! error-corrected cellular-barcode/UMI BAM tags and render nucleotide
//! edit-count matrices, plus the `pheatmap` R script generator shared
//! with `sc_seqQual.py`. Ports `qcmodule.scbam.barcode_edits`/
//! `diff_str` and `qcmodule.heatmap.make_heatmap`.
//!
//! **Preserves upstream quirks, do not "fix"**:
//! - `diff_str` only reports per-position edits when the raw and
//!   corrected barcodes are the SAME LENGTH; a length mismatch still
//!   counts toward `*_diff` (plain string inequality), but contributes
//!   ZERO edit-matrix entries (no positional comparison is possible).
//! - Frequency tables are sorted by descending count using a STABLE
//!   sort (Python's `sorted(..., reverse=True)`), so ties keep their
//!   first-seen (BAM record) order -- replicated with an explicit
//!   insertion-order side list, not a `HashMap`'s arbitrary iteration
//!   order.
//! - The edit-count CSV replicates `pandas.DataFrame.from_dict(...).
//!   fillna(0)`'s behavior of a mixed sparse dict becoming an
//!   ALL-FLOAT table (every cell prints with a trailing `.0`, even
//!   though every value is a whole count) -- confirmed via a real
//!   `pandas` run in `oracle/venv`, not just inferred from source.
//! - `heatmap.make_heatmap`'s actual Rscript subprocess invocation is
//!   HARDCODED to the literal string `"Rscript"`, ignoring whatever
//!   `--rscript` executable the CLI's own (separate) dependency-check
//!   step resolved -- a real upstream inconsistency between the
//!   dependency-check step and the actual render step, replicated
//!   exactly (not unified to always use the user's `--rscript` value).
use std::collections::{BTreeMap, HashMap};
use std::io;

use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::record::data::field::{Tag, Value}};

fn tag(name: &str) -> Tag {
    let b = name.as_bytes();
    Tag::new(b[0], b[1])
}

fn tag_string<'d, D: sam::alignment::record::Data<'d>>(data: &D, t: Tag) -> Option<String> {
    match data.get(&t)?.ok()? {
        Value::String(s) => Some(String::from_utf8_lossy(s).into_owned()),
        _ => None,
    }
}

#[derive(Debug, Clone)]
pub struct BarcodeTagNames {
    pub cr: String,
    pub cb: String,
    pub ur: String,
    pub ub: String,
}

impl Default for BarcodeTagNames {
    fn default() -> Self {
        Self { cr: "CR".to_string(), cb: "CB".to_string(), ur: "UR".to_string(), ub: "UB".to_string() }
    }
}

/// Ports `diff_str`: per-position `(index, original_char, corrected_char)`
/// triples where `s1`/`s2` differ, or empty if their lengths differ.
pub fn diff_str(s1: &str, s2: &str) -> Vec<(usize, char, char)> {
    let c1: Vec<char> = s1.chars().collect();
    let c2: Vec<char> = s2.chars().collect();
    if c1.len() != c2.len() {
        return Vec::new();
    }
    c1.iter().zip(c2.iter()).enumerate().filter(|(_, (a, b))| a != b).map(|(i, (&a, &b))| (i, a, b)).collect()
}

/// One barcode kind's (CB or UMI) accumulated statistics.
#[derive(Debug, Default, Clone)]
pub struct EditCounts {
    pub miss: i64,
    pub same: i64,
    pub diff: i64,
    /// corrected-value -> occurrence count.
    pub freq: HashMap<String, i64>,
    /// corrected-value first-seen order, for stable descending-count sort.
    pub freq_order: Vec<String>,
    /// position -> ("from:to" substitution -> count).
    pub corrected_bases: BTreeMap<i64, BTreeMap<String, i64>>,
}

impl EditCounts {
    fn record(&mut self, original: &str, corrected: &str) {
        if !self.freq.contains_key(corrected) {
            self.freq_order.push(corrected.to_string());
        }
        *self.freq.entry(corrected.to_string()).or_insert(0) += 1;

        if original != corrected {
            self.diff += 1;
            for (pos, from, to) in diff_str(original, corrected) {
                *self.corrected_bases.entry(pos as i64).or_default().entry(format!("{from}:{to}")).or_insert(0) += 1;
            }
        } else {
            self.same += 1;
        }
    }
}

#[derive(Debug, Default, Clone)]
pub struct BarcodeEditStats {
    pub total_alignments: i64,
    pub cb: EditCounts,
    pub umi: EditCounts,
}

/// Runs the full per-alignment barcode/UMI edit-counting pass. Ports
/// `barcode_edits`. `limit` stops after that many alignments have been
/// processed (upstream's `--limit`).
pub fn barcode_edits<I>(records: I, tags: &BarcodeTagNames, limit: Option<i64>) -> io::Result<BarcodeEditStats>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let cr_tag = tag(&tags.cr);
    let cb_tag = tag(&tags.cb);
    let ur_tag = tag(&tags.ur);
    let ub_tag = tag(&tags.ub);

    let mut s = BarcodeEditStats::default();
    let mut iter = records.into_iter();

    loop {
        // Upstream: `total_alignments += 1` runs BEFORE `next(samfile)`
        // inside the loop body, so when the record stream is exhausted
        // the increment has ALREADY happened and is never rolled back
        // (`StopIteration` is caught OUTSIDE the loop) -- a genuine
        // upstream off-by-one (`total_alignments` == real record count
        // + 1 whenever the stream runs out before `limit` is reached).
        // Preserved here rather than "fixed": reported via
        // `BarcodeEditStats::total_alignments` and printed by the CLI's
        // "Total alignments processed: N" line.
        let Some(item) = iter.next() else {
            s.total_alignments += 1;
            break;
        };
        let record = item?;
        s.total_alignments += 1;
        let data = record.data();

        match (tag_string(&data, cr_tag), tag_string(&data, cb_tag)) {
            (Some(o), Some(c)) => s.cb.record(&o.replace("-1", ""), &c.replace("-1", "")),
            _ => s.cb.miss += 1,
        }
        match (tag_string(&data, ur_tag), tag_string(&data, ub_tag)) {
            (Some(o), Some(c)) => s.umi.record(&o.replace("-1", ""), &c.replace("-1", "")),
            _ => s.umi.miss += 1,
        }

        if let Some(l) = limit {
            if s.total_alignments >= l {
                break;
            }
        }
    }

    Ok(s)
}

/// Renders a `<barcode>\t<count>\n` frequency table, sorted by
/// descending count with a STABLE tie-break (first-seen order).
pub fn render_freq_tsv(counts: &EditCounts) -> String {
    let mut order = counts.freq_order.clone();
    order.sort_by(|a, b| counts.freq[b].cmp(&counts.freq[a]));
    let mut out = String::new();
    for bc in order {
        out.push_str(&format!("{bc}\t{}\n", counts.freq[&bc]));
    }
    out
}

/// Renders the edit-count CSV, replicating `pandas.DataFrame.
/// from_dict(...).fillna(0)`'s dtype behavior: `Index,<pos>,<pos>,...`
/// header (positions ascending), one row per substitution key
/// (alphabetically ascending). Dtype is decided PER COLUMN, not
/// globally: a position's column stays plain-integer ONLY if that
/// position's inner dict already had an entry for every substitution
/// key in the union (fully dense -- `fillna(0)` never touches it, so
/// pandas keeps whatever int dtype the raw counts had); if the position
/// was missing even one substitution key present at some OTHER
/// position, the whole column becomes float64 (every cell in it prints
/// with a trailing `.0`, including cells that had a real nonzero count
/// -- NaN-filling promotes the entire column, not just the filled
/// cells). Confirmed via live `pandas.DataFrame.from_dict` probes
/// against several dense/sparse shapes, not inferred from source alone.
/// Empty input renders just `Index\n`, matching pandas' empty-DataFrame
/// `to_csv` output.
pub fn render_edit_matrix_csv(corrected_bases: &BTreeMap<i64, BTreeMap<String, i64>>) -> String {
    let positions: Vec<i64> = corrected_bases.keys().copied().collect();

    let mut substitutions: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    for by_sub in corrected_bases.values() {
        substitutions.extend(by_sub.keys().cloned());
    }

    let mut out = String::from("Index");
    for pos in &positions {
        out.push(',');
        out.push_str(&pos.to_string());
    }
    out.push('\n');

    let column_is_dense: HashMap<i64, bool> = positions
        .iter()
        .map(|&pos| {
            let dense = corrected_bases.get(&pos).is_some_and(|m| substitutions.iter().all(|s| m.contains_key(s)));
            (pos, dense)
        })
        .collect();

    for sub in &substitutions {
        out.push_str(sub);
        for pos in &positions {
            let count = corrected_bases.get(pos).and_then(|m| m.get(sub)).copied().unwrap_or(0);
            if column_is_dense[pos] {
                out.push_str(&format!(",{count}"));
            } else {
                out.push_str(&format!(",{count}.0"));
            }
        }
        out.push('\n');
    }

    out
}

/// Renders the `pheatmap` R script. Ports `heatmap.make_heatmap`'s
/// R-code-generation half (the `Rscript`-invocation half is a CLI
/// concern, not part of this pure function).
#[allow(clippy::too_many_arguments)]
pub fn render_heatmap_r_script(infile: &str, outfile_prefix: &str, filetype: &str, cell_width: i64, cell_height: i64, col_angle: i64, font_size: i64, text_color: &str, no_numbers: bool, log2_scale: bool) -> String {
    let mut out = String::new();
    out.push_str("if(!require(pheatmap)){install.packages(\"pheatmap\")}\n");
    out.push_str("library(pheatmap)\n");
    out.push_str(&format!("dat = read.table(file = '{infile}',sep=',',header=T,row.names=1,check.names=FALSE)\n"));

    let plot_file = format!("{outfile_prefix}.{filetype}");
    let matrix_expr = if log2_scale { "log2(as.matrix(dat)+1)".to_string() } else { "as.matrix(dat)".to_string() };

    let line = if no_numbers {
        format!(
            "pheatmap({matrix_expr}, filename='{plot_file}', cellwidth = {cell_width}, cellheight = {cell_height}, display_numbers = FALSE, angle_col={col_angle}, fontsize={font_size},cluster_rows=F, cluster_cols=F, scale='none',color = colorRampPalette(c('#f1eef6','#d7b5d8','#df65b0', '#ce1256'))(50))\n"
        )
    } else {
        format!(
            "pheatmap({matrix_expr}, filename='{plot_file}', cellwidth = {cell_width}, cellheight = {cell_height}, display_numbers = TRUE, angle_col={col_angle}, fontsize={font_size}, number_color='{text_color}',cluster_rows=F, cluster_cols=F, scale='none',color = colorRampPalette(c('#f1eef6','#d7b5d8','#df65b0', '#ce1256'))(50))\n"
        )
    };
    out.push_str(&line);

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use sam::alignment::record::{Flags, cigar::Op, cigar::op::Kind};
    use sam::alignment::record_buf::{Cigar, RecordBuf, data::Data, data::field::Value as BufValue};

    fn header_with_chrom(name: &str, len: usize) -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence(name, sam::header::record::value::Map::<sam::header::record::value::map::ReferenceSequence>::new(std::num::NonZeroUsize::new(len).unwrap()))
            .build()
    }

    fn to_bam_records(header: &sam::Header, records: &[RecordBuf]) -> Vec<bam::Record> {
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

    fn record_with_tags(name: &str, cr: Option<&str>, cb: Option<&str>) -> RecordBuf {
        let mut data = Data::default();
        if let Some(v) = cr {
            data.insert(tag("CR"), BufValue::String(v.into()));
        }
        if let Some(v) = cb {
            data.insert(tag("CB"), BufValue::String(v.into()));
        }
        RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(1).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 4)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 4]))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![40; 4]))
            .set_data(data)
            .build()
    }

    #[test]
    fn diff_str_reports_positions_only_for_equal_length_strings() {
        assert_eq!(diff_str("AAAA", "ATAA"), vec![(1, 'A', 'T')]);
        assert_eq!(diff_str("AAAA", "AAAA"), vec![]);
        assert_eq!(diff_str("AAA", "AAAA"), vec![]); // length mismatch -> empty
    }

    #[test]
    fn barcode_edits_counts_miss_same_diff_and_matrix() {
        let header = header_with_chrom("chr1", 1000);
        let recs = vec![
            record_with_tags("r1", Some("AAAA-1"), Some("ATAA-1")), // edited at pos1 A->T
            record_with_tags("r2", Some("CCCC"), Some("CCCC")),     // same
            record_with_tags("r3", None, None),                    // missing both tags
        ];
        let bam_records = to_bam_records(&header, &recs);
        let tags = BarcodeTagNames::default();
        let stats = barcode_edits(bam_records.into_iter().map(Ok), &tags, None).unwrap();

        // Regression coverage for a real upstream off-by-one, confirmed
        // live (3-record BAM through the actual sc_editMatrix.py CLI
        // prints "Total alignments processed: 4"): `total_alignments`
        // is incremented BEFORE the loop's `next(samfile)` call, and
        // that increment is never rolled back when the stream is
        // exhausted (`StopIteration` is caught outside the loop) --
        // so with no `--limit`, `total_alignments` == real record count
        // + 1, not the real count.
        assert_eq!(stats.total_alignments, 4);
        assert_eq!(stats.cb.diff, 1);
        assert_eq!(stats.cb.same, 1);
        assert_eq!(stats.cb.miss, 1);
        assert_eq!(stats.cb.corrected_bases[&1]["A:T"], 1);
        // "-1" suffix stripped before comparison and frequency counting.
        assert_eq!(stats.cb.freq["ATAA"], 1);
        assert_eq!(stats.cb.freq["CCCC"], 1);
    }

    #[test]
    fn barcode_edits_respects_limit() {
        let header = header_with_chrom("chr1", 1000);
        let recs = vec![record_with_tags("r1", Some("AAAA"), Some("AAAA")), record_with_tags("r2", Some("CCCC"), Some("CCCC")), record_with_tags("r3", Some("GGGG"), Some("GGGG"))];
        let bam_records = to_bam_records(&header, &recs);
        let tags = BarcodeTagNames::default();
        let stats = barcode_edits(bam_records.into_iter().map(Ok), &tags, Some(2)).unwrap();
        assert_eq!(stats.total_alignments, 2);
    }

    #[test]
    fn render_freq_tsv_sorts_by_descending_count_stable_on_ties() {
        let mut counts = EditCounts::default();
        // Insert in this order: B(1), A(2), C(1) -- ties (B,C both count 1)
        // must keep first-seen order (B before C) per Python's stable sort.
        for bc in ["B", "A", "A", "C"] {
            counts.record(bc, bc); // same == same, just builds freq/freq_order
        }
        let tsv = render_freq_tsv(&counts);
        let lines: Vec<&str> = tsv.lines().collect();
        assert_eq!(lines, vec!["A\t2", "B\t1", "C\t1"]);
    }

    #[test]
    fn render_edit_matrix_csv_matches_real_pandas_output() {
        // Cross-checked byte-for-byte against a real pandas
        // DataFrame.from_dict(...).fillna(0).to_csv(...) run in
        // oracle/venv with this exact input shape.
        let mut matrix: BTreeMap<i64, BTreeMap<String, i64>> = BTreeMap::new();
        matrix.entry(1).or_default().insert("A:T".to_string(), 1);
        matrix.entry(1).or_default().insert("A:C".to_string(), 1);
        matrix.entry(2).or_default().insert("C:G".to_string(), 1);
        matrix.entry(3).or_default().insert("G:T".to_string(), 1);

        let csv = render_edit_matrix_csv(&matrix);
        let expected = "Index,1,2,3\nA:C,1.0,0.0,0.0\nA:T,1.0,0.0,0.0\nC:G,0.0,1.0,0.0\nG:T,0.0,0.0,1.0\n";
        assert_eq!(csv, expected);
    }

    #[test]
    fn render_edit_matrix_csv_dtype_is_decided_per_column_not_globally() {
        // Cross-checked byte-for-byte against a real
        // pandas.DataFrame.from_dict(...).fillna(0).to_csv(...) run:
        // column 2 only has an "A:T" entry (missing "C:G", which
        // appears at columns 1 and 3), so ONLY column 2 gets promoted
        // to float by fillna(0) -- columns 1 and 3 are fully dense
        // (every row present) and stay plain integers, even though the
        // matrix as a whole is NOT fully dense.
        let mut matrix: BTreeMap<i64, BTreeMap<String, i64>> = BTreeMap::new();
        matrix.entry(1).or_default().insert("A:T".to_string(), 2);
        matrix.entry(1).or_default().insert("C:G".to_string(), 1);
        matrix.entry(2).or_default().insert("A:T".to_string(), 5);
        matrix.entry(3).or_default().insert("A:T".to_string(), 1);
        matrix.entry(3).or_default().insert("C:G".to_string(), 3);

        let csv = render_edit_matrix_csv(&matrix);
        let expected = "Index,1,2,3\nA:T,2,5.0,1\nC:G,1,0.0,3\n";
        assert_eq!(csv, expected);
    }

    #[test]
    fn render_edit_matrix_csv_single_dense_column_stays_integer() {
        // A single (position, substitution) pair is trivially dense --
        // no NaN is ever introduced, so pandas keeps the int64 dtype and
        // `to_csv` prints "1", not "1.0". This is the exact shape that
        // was previously mis-rendered as "1.0" (the bug this test
        // guards against).
        let mut matrix: BTreeMap<i64, BTreeMap<String, i64>> = BTreeMap::new();
        matrix.entry(3).or_default().insert("A:T".to_string(), 1);
        assert_eq!(render_edit_matrix_csv(&matrix), "Index,3\nA:T,1\n");
    }

    #[test]
    fn render_edit_matrix_csv_empty_is_just_the_header() {
        assert_eq!(render_edit_matrix_csv(&BTreeMap::new()), "Index\n");
    }

    #[test]
    fn render_heatmap_r_script_exact_text_all_four_branches() {
        let base = ("in.csv", "out.CB_edits_heatmap", "pdf", 15, 10, 45, 8, "black");
        let with_numbers_linear = render_heatmap_r_script(base.0, base.1, base.2, base.3, base.4, base.5, base.6, base.7, false, false);
        assert!(with_numbers_linear.contains("pheatmap(as.matrix(dat), filename='out.CB_edits_heatmap.pdf', cellwidth = 15, cellheight = 10, display_numbers = TRUE, angle_col=45, fontsize=8, number_color='black',cluster_rows=F, cluster_cols=F, scale='none',color = colorRampPalette(c('#f1eef6','#d7b5d8','#df65b0', '#ce1256'))(50))\n"));

        let no_numbers_log2 = render_heatmap_r_script(base.0, base.1, base.2, base.3, base.4, base.5, base.6, base.7, true, true);
        assert!(no_numbers_log2.contains("pheatmap(log2(as.matrix(dat)+1), filename='out.CB_edits_heatmap.pdf', cellwidth = 15, cellheight = 10, display_numbers = FALSE, angle_col=45, fontsize=8,cluster_rows=F, cluster_cols=F, scale='none',color = colorRampPalette(c('#f1eef6','#d7b5d8','#df65b0', '#ce1256'))(50))\n"));
    }
}
