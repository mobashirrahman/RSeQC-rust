//! Port of `read_quality.py`: calculate per-cycle Phred quality-score distributions.
//! Contract: see `compatibility/commands.yaml` entry `read_quality.py`; algorithm
//! ported from `readsQual_boxplot()` in `oracle/upstream-src/src/qcmodule/SAM.py`
//! (lines 3051-3117).
//!
//! **Notable upstream quirk, confirmed by reading the literal source (do
//! NOT "fix")**: `readsQual_boxplot`'s `is_unmapped`/`is_qcfail` checks
//! are commented out in upstream (`#if aligned_read.is_unmapped:continue`).
//! Only the MAPQ filter is actually live -- unmapped, QC-fail, duplicate,
//! and secondary-alignment reads are ALL included as long as their MAPQ
//! clears `q_cut`. (An earlier version of this port's doc comment and
//! code incorrectly filtered unmapped/QC-fail reads anyway; found and
//! fixed via `verification/run_diff.py`'s `read_quality_basic` case,
//! which caught a real read-count mismatch against the actual upstream
//! CLI.) Quality scores are aggregated per position; the output table's
//! row count equals the last processed read's length (not the max
//! length seen across all records) -- same quirk as read_NVC.py.
//!
//! SAM-text and CRAM input (DIV-0002/0004) are supported at the CLI layer via
//! `rseqc_formats::open_alignments`; this module's own functions were
//! unaffected (already generic over
//! `IntoIterator<Item = io::Result<bam::Record>>`).
//!
//! **Second upstream quirk, found via a real per-read diff against a kept
//! synthetic-sweep failure**: `aligned_read.qqual` is pysam's *deprecated
//! alias for `query_alignment_qualities`*, NOT `query_qualities` -- despite
//! the docstring ("calculate phred quality score for each base in read")
//! suggesting the whole read. `query_alignment_qualities` is
//! `query_qualities[query_alignment_start:query_alignment_end]`: it
//! excludes SOFT-CLIPPED bases at either end of the query (confirmed live:
//! `pysam.AlignedSegment.qqual` on a read with a leading/trailing `S` CIGAR
//! op returns exactly `query_alignment_qualities`, byte-for-byte, while
//! `query_qualities` includes the clipped bases). Hard clips never appear
//! in SEQ/QUAL at all (already excluded by BAM's own encoding), so only
//! `S` ops matter here. A read with no CIGAR (typically: unmapped) has no
//! clip information, so its "alignment" quality scores are the full
//! sequence, same as pysam's fallback. `trim_soft_clips` reproduces this
//! by finding the leading/trailing contiguous soft-clip length (skipping
//! over any flanking hard clip, which doesn't consume query bases) and
//! slicing the full-length quality vector to match, BEFORE the
//! is_reverse reversal below (same order upstream applies: `qqual` is
//! already clip-trimmed storage-order data, then optionally reversed).

use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::Write as _;
use std::process::{Command, Stdio};

use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::op::Kind as CigarOpKind;

use crate::exec_resolve;

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct QualityHistogram {
    pub quality: HashMap<usize, HashMap<u8, u64>>, // position -> quality_score -> count
    pub q_min: u8,
    pub q_max: u8,
    pub read_len: usize, // Last processed record's length (the quirk)
}

/// Slices `full_qual_scores` down to pysam's `query_alignment_qualities`
/// range: drops a leading and/or trailing soft-clip run (any flanking
/// hard clip is skipped over, since it doesn't consume query bases and so
/// was never present in `full_qual_scores` to begin with). A CIGAR with
/// no soft clips at that end (including an empty/absent CIGAR, e.g. for
/// an unmapped read) leaves that end untouched, matching pysam's
/// full-sequence fallback when there is no clip information.
fn trim_soft_clips(full_qual_scores: &[u8], ops: &[noodles_sam::alignment::record::cigar::Op]) -> Vec<u8> {
    let start = ops
        .iter()
        .find(|op| op.kind() != CigarOpKind::HardClip)
        .filter(|op| op.kind() == CigarOpKind::SoftClip)
        .map(|op| op.len())
        .unwrap_or(0);

    let end_clip = ops
        .iter()
        .rev()
        .find(|op| op.kind() != CigarOpKind::HardClip)
        .filter(|op| op.kind() == CigarOpKind::SoftClip)
        .map(|op| op.len())
        .unwrap_or(0);

    let end = full_qual_scores.len().saturating_sub(end_clip);
    let start = start.min(end);
    full_qual_scores[start..end].to_vec()
}

/// Computes quality score histogram for a sequence of BAM records,
/// matching the upstream algorithm exactly.
/// Pure computation: no I/O, no printing.
pub fn compute_quality<I>(records: I, q_cut: u8) -> io::Result<QualityHistogram>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut quality = HashMap::new();
    let mut q_max = 1; // Start at 1 since upstream uses -1 initially
    let mut q_min = 93; // Start at 93 since upstream uses 10000 initially
    let mut read_len = 0;

    for result in records {
        let record = result?;

        // Only the MAPQ filter is actually live upstream -- see the
        // module doc comment for why unmapped/QC-fail reads are NOT
        // filtered here, unlike almost every other command in this port.
        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        if mapq < q_cut {
            continue;
        }

        // Get quality scores (raw Phred scores, not ASCII), then trim
        // soft-clipped bases to match upstream's `qqual`
        // (== `query_alignment_qualities`) -- see the module doc comment.
        let full_qual_scores: Vec<u8> = record.quality_scores().iter().collect();
        let ops: Vec<_> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;
        let mut qual_scores = trim_soft_clips(&full_qual_scores, &ops);
        // Upstream's table length is `aligned_read.rlen` (full l_qseq,
        // soft clips included), not the trimmed `qqual` length.
        read_len = record.sequence().len();

        // Reverse quality scores if the read is reverse complemented
        if record.flags().is_reverse_complemented() {
            qual_scores.reverse();
        }
        
        // Process each position's quality score
        for (i, q) in qual_scores.iter().enumerate() {
            q_max = q_max.max(*q);
            q_min = q_min.min(*q);
            quality.entry(i).or_insert_with(HashMap::new);
            *quality.get_mut(&i).unwrap().entry(*q).or_insert(0) += 1;
        }
    }

    Ok(QualityHistogram {
        quality,
        q_min,
        q_max,
        read_len,
    })
}

/// Renders the quality histogram as an R script, matching upstream's .qual.r format
pub fn render_qual_r_script(hist: &QualityHistogram, shrink: u64, output_prefix: &str) -> String {
    let mut lines = Vec::new();
    
    // Build per-position R vector-construction expressions
    let mut i_box = HashMap::new(); // position -> R expression string
    let mut q_list = Vec::new(); // flat list of counts, position-major, q-minor
    
    for p in 0..hist.read_len {
        let mut val = Vec::new();
        let mut occurrence = Vec::new();
        
        for q in hist.q_min..=hist.q_max {
            if let Some(pos_scores) = hist.quality.get(&p) {
                if let Some(count) = pos_scores.get(&q) {
                    val.push(q.to_string());
                    occurrence.push(count.to_string());
                    q_list.push(count.to_string());
                } else {
                    q_list.push("0".to_string());
                }
            } else {
                q_list.push("0".to_string());
            }
        }
        
        if !val.is_empty() {
            i_box.insert(p, format!(
                "rep(c({}),times=c({})/{})",
                val.join(","),
                occurrence.join(","),
                shrink
            ));
        } else {
            i_box.insert(p, "rep(c(),times=c())".to_string());
        }
    }
    
    // Output: single file "<prefix>.qual.r", written as (each on its own line, in order):
    lines.push(format!("pdf('{}.qual.boxplot.pdf')", output_prefix));
    
    // Print position vectors in ascending order
    for i in 0..hist.read_len {
        if let Some(expr) = i_box.get(&i) {
            lines.push(format!("p{}<-{}", i, expr));
        }
    }
    
    // Upstream always prints this line, even with no data: an empty
    // `i_box` makes `','.join([])` == "", so the Python output is the
    // (oddly-shaped but real) `boxplot(,xlab=...)`. Reproduce that exactly
    // rather than substituting a blank line.
    let plot_vars: Vec<String> = (0..hist.read_len).map(|i| format!("p{}", i)).collect();
    lines.push(format!(
        "boxplot({},xlab=\"Position of Read(5'->3')\",ylab=\"Phred Quality Score\",outline=F)",
        plot_vars.join(",")
    ));
    
    lines.push("dev.off()".to_string());
    // Upstream's `print('\n', file=FO)` writes the literal string "\n"
    // PLUS print's own trailing newline -- that's TWO blank lines
    // between the sections, not one.
    lines.push("".to_string());
    lines.push("".to_string());

    // Heatmap section
    lines.push(format!("pdf('{}.qual.heatmap.pdf')", output_prefix));
    lines.push(format!("qual=c({})", q_list.join(",")));
    lines.push(format!("mat=matrix(qual,ncol={},byrow=F)", hist.read_len));
    lines.push("Lab.palette <- colorRampPalette(c(\"blue\", \"orange\", \"red3\",\"red2\",\"red1\",\"red\"), space = \"rgb\",interpolate=c('spline'))".to_string());
    lines.push(format!(
        "heatmap(mat,Rowv=NA,Colv=NA,xlab=\"Position of Read\",ylab=\"Phred Quality Score\",labRow=seq(from={},to={}),col = Lab.palette(256),scale=\"none\" )",
        hist.q_min, hist.q_max
    ));
    lines.push("dev.off()".to_string());

    // Trailing newline: upstream's plain `print(...)` calls each add
    // their own trailing newline, including the final `dev.off()`.
    format!("{}\n", lines.join("\n"))
}

/// Refuses an output prefix whose parent directory does not exist, before
/// any work. Mirrors `rseqc_cli::require_existing_output_parent` exactly
/// (same `Path::parent` mapping, same message); kept local because the
/// commands crate cannot depend on the CLI crate.
fn require_output_parent(output_prefix: &str) -> io::Result<()> {
    use std::path::Path;
    let prefix = Path::new(output_prefix);
    let parent = match prefix.parent() {
        Some(p) if p.as_os_str().is_empty() => Path::new("."),
        Some(p) => p,
        None => Path::new("."),
    };
    if parent.is_dir() {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("output directory does not exist: {}", parent.display()),
        ))
    }
}

/// Runs the `read_quality.py` CLI body over an already-opened record stream
/// (C1 multi-driver pattern).
///
/// This is exactly what the standalone binary's `run()` does -- same
/// progress lines, same files with the same bytes, same Rscript contract,
/// same error propagation -- except the record source is a
/// caller-supplied iterator and stdout/stderr are caller-supplied sinks.
/// The standalone binary delegates to this (passing the process streams);
/// `rseqc_multi` passes one record broadcast plus per-command stream
/// files. `compute_quality` and the renderer are untouched.
///
/// Transport note (same as `run_read_gc`/`run_read_nvc`): the standalone
/// binary used to run Rscript with inherited stdio (`Command::status`);
/// here the child is spawned piped and its captured stdout/stderr are
/// copied to the command's sinks, so the bytes stay attributable under
/// multi while remaining identical standalone.
///
/// The argument list is long on purpose and carries a targeted lint
/// allowance: this is the C1 multi-driver pattern (one callable per
/// command carrying its full CLI surface plus the two sinks), and
/// bundling the flags into a struct would only hide them from the C2
/// cards that copy this signature command by command.
#[allow(clippy::too_many_arguments)]
pub fn run_read_quality<I>(
    records: I,
    q_cut: u8,
    out_prefix: &str,
    reduce: u64,
    skip_plot: bool,
    rscript: &str,
    stdout: &mut dyn io::Write,
    stderr: &mut dyn io::Write,
) -> io::Result<()>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    // Same rule as the binary's `require_existing_output_parent_or_exit`
    // (and upstream's `validate_args`), in `Result` form; see above.
    require_output_parent(out_prefix)?;

    // Upstream: `if self.bam_format: print("Read BAM file ... ", end=' ')`
    // -- always the BAM branch in practice (htslib's `'rb'` open is
    // lenient about actual content and succeeds for genuine plain-text
    // SAM too). The literal's own trailing space plus `end=' '` gives
    // two spaces before "Done".
    write!(stderr, "Read BAM file ...  ")?;
    let hist = compute_quality(records, q_cut)?;
    writeln!(stderr, "Done")?;

    let r_path = format!("{out_prefix}.qual.r");
    File::create(&r_path)?.write_all(render_qual_r_script(&hist, reduce, out_prefix).as_bytes())?;

    if !skip_plot {
        let rscript_path = exec_resolve::which(rscript).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("Rscript executable not found: {rscript}"),
            )
        })?;
        let out = Command::new(&rscript_path)
            .arg(&r_path)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()?
            .wait_with_output()?;
        stdout.write_all(&out.stdout)?;
        stderr.write_all(&out.stderr)?;
        if !out.status.success() {
            return Err(io::Error::other(format!("R plotting failed for {r_path}")));
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::{
        self as sam,
        alignment::{
            io::Write as _,
            record::{Flags, MappingQuality, cigar::Op, cigar::op::Kind as CigarOpKind},
            record_buf::{Cigar, RecordBuf},
        },
        header::record::value::{Map, map::ReferenceSequence},
    };
    use std::num::NonZeroUsize;

    fn test_header() -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
            )
            .build()
    }

    /// Round-trips hand-built records through a real BAM byte stream so the
    /// test exercises the same decode path `compute_quality` runs against in
    /// production, not just in-memory builder state.
    fn to_bam_records(header: &sam::Header, records: &[RecordBuf]) -> Vec<bam::Record> {
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
    fn unmapped_and_qcfail_reads_are_not_filtered_when_mapq_clears_the_cutoff() {
        // Regression test for a real bug found via verification/
        // run_diff.py's read_quality_basic case: upstream's
        // is_unmapped/is_qcfail checks in readsQual_boxplot are
        // commented out (only the MAPQ filter is live), but an earlier
        // version of this port filtered them anyway. Both reads here
        // have MAPQ 40 (clears the default cutoff) despite being
        // flagged unmapped/QC-fail, and must still be counted.
        let header = test_header();

        let unmapped = RecordBuf::builder()
            .set_flags(Flags::UNMAPPED)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 2)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"AA".to_vec()))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![10, 20]))
            .build();

        let qcfail = RecordBuf::builder()
            .set_flags(Flags::QC_FAIL)
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 2)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"AA".to_vec()))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![15, 25]))
            .build();

        let bam_records = to_bam_records(&header, &[unmapped, qcfail]);
        let hist = compute_quality(bam_records.into_iter().map(Ok), 30).unwrap();

        assert_eq!(hist.read_len, 2);
        assert_eq!(hist.quality[&0][&10], 1);
        assert_eq!(hist.quality[&0][&15], 1);
        assert_eq!(hist.quality[&1][&20], 1);
        assert_eq!(hist.quality[&1][&25], 1);
    }

    #[test]
    fn soft_clipped_bases_are_excluded_like_pysam_qqual() {
        // Regression test for a real per-read divergence found via a kept
        // synthetic-sweep failure (read_quality_pe): upstream's
        // `aligned_read.qqual` is pysam's deprecated alias for
        // `query_alignment_qualities`, which excludes soft-clipped bases
        // at either end of the query -- NOT `query_qualities` (the full
        // sequence). A read with a 2bp leading soft clip and a 1bp
        // trailing soft clip must only contribute its 3 aligned-region
        // quality scores, not all 6.
        let header = test_header();

        let leading_and_trailing_clip = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![
                Op::new(CigarOpKind::SoftClip, 2),
                Op::new(CigarOpKind::Match, 3),
                Op::new(CigarOpKind::SoftClip, 1),
            ]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"AAAAAA".to_vec()))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![
                90, 91, 10, 20, 30, 92,
            ]))
            .build();

        // A hard clip flanking a soft clip must be skipped over (it
        // consumes no query bases, so it was never in the quality vector
        // to begin with) rather than blocking the soft-clip detection.
        let hard_then_soft_clip = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![
                Op::new(CigarOpKind::HardClip, 5),
                Op::new(CigarOpKind::SoftClip, 1),
                Op::new(CigarOpKind::Match, 2),
            ]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"AAA".to_vec()))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![93, 40, 50]))
            .build();

        let bam_records = to_bam_records(&header, &[leading_and_trailing_clip, hard_then_soft_clip]);
        let hist = compute_quality(bam_records.into_iter().map(Ok), 30).unwrap();

        // Table length is the last record's `rlen` (full SEQ length, soft
        // clip included, hard clip excluded) -- 3, not its 2 aligned bases.
        // Confirmed with pysam: 5H1S2M over "AAA" gives rlen 3, qqual len 2.
        assert_eq!(hist.read_len, 3);
        // leading_and_trailing_clip: aligned-region quals are [10, 20, 30].
        assert_eq!(hist.quality[&0][&10], 1);
        assert_eq!(hist.quality[&1][&20], 1);
        // hard_then_soft_clip: aligned-region quals are [40, 50] (the
        // leading hard clip's 5bp and soft clip's 1bp are both excluded).
        assert_eq!(hist.quality[&0][&40], 1);
        assert_eq!(hist.quality[&1][&50], 1);
        // None of the clipped sentinel values (90/91/92/93) were counted.
        for pos in hist.quality.values() {
            for sentinel in [90u8, 91, 92, 93] {
                assert!(!pos.contains_key(&sentinel), "clipped base leaked into histogram");
            }
        }
    }

    #[test]
    fn matches_hand_computed_counts() {
        let header = test_header();

        // Forward read with quality scores: 30, 40, 50 (length 3)
        let forward = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 3)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"AAA".to_vec()))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![30, 40, 50]))
            .build();

        // Reverse read with quality scores: 20, 30, 40 (length 3) - should be reversed to 40, 30, 20
        let reverse = RecordBuf::builder()
            .set_flags(Flags::REVERSE_COMPLEMENTED)
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 3)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"AAA".to_vec()))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![20, 30, 40]))
            .build();

        // Short read with quality scores: 10 (length 1) - this should determine the table length
        let short = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 1)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"A".to_vec()))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![10]))
            .build();

        // Low MAPQ read (should be filtered out)
        let low_mapq = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(10).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 2)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"AA".to_vec()))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![25, 35]))
            .build();

        let records = vec![forward, reverse, short, low_mapq];
        let bam_records = to_bam_records(&header, &records);
        let hist = compute_quality(bam_records.into_iter().map(Ok), 30).unwrap();

        // Hand-computed from the fixture above:
        // - forward: positions 0:30, 1:40, 2:50 (filtered because short read determines length)
        // - reverse: positions 0:40, 1:30 (reversed from 20,30,40 -> 40,30,20, but only first 2 kept)
        // - short: position 0:10
        // - low_mapq: filtered out
        // Table length is 1 (the LAST record processed is `short`, length 1)
        
        assert_eq!(hist.read_len, 1, "Table length should match last read's length");
        assert_eq!(hist.q_min, 10, "Min quality should be 10");
        assert_eq!(hist.q_max, 50, "Max quality should be 50");

        // Position 0: 30 (forward) + 40 (reverse) + 10 (short) = 80 total
        let pos0_counts = hist.quality.get(&0).unwrap();
        assert_eq!(pos0_counts.get(&10), Some(&1), "Quality 10 count should be 1");
        assert_eq!(pos0_counts.get(&30), Some(&1), "Quality 30 count should be 1");
        assert_eq!(pos0_counts.get(&40), Some(&1), "Quality 40 count should be 1");
    }

    #[test]
    fn render_qual_r_script_format() {
        let hist = QualityHistogram {
            quality: {
                let mut map = HashMap::new();
                let mut pos0 = HashMap::new();
                pos0.insert(10, 5);
                pos0.insert(20, 3);
                map.insert(0, pos0);
                let mut pos1 = HashMap::new();
                pos1.insert(15, 2);
                pos1.insert(25, 4);
                map.insert(1, pos1);
                map
            },
            q_min: 10,
            q_max: 25,
            read_len: 2,
        };

        let output = render_qual_r_script(&hist, 1000, "test_output");
        // Cross-checked byte-for-byte against a `python3 -c` run of the
        // literal upstream print() sequence (single-quoted pdf(), and
        // print('\n', file=FO)'s content-newline PLUS its own trailing
        // newline giving two blank lines between the sections).
        let expected = r#"pdf('test_output.qual.boxplot.pdf')
p0<-rep(c(10,20),times=c(5,3)/1000)
p1<-rep(c(15,25),times=c(2,4)/1000)
boxplot(p0,p1,xlab="Position of Read(5'->3')",ylab="Phred Quality Score",outline=F)
dev.off()


pdf('test_output.qual.heatmap.pdf')
qual=c(5,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,4)
mat=matrix(qual,ncol=2,byrow=F)
Lab.palette <- colorRampPalette(c("blue", "orange", "red3","red2","red1","red"), space = "rgb",interpolate=c('spline'))
heatmap(mat,Rowv=NA,Colv=NA,xlab="Position of Read",ylab="Phred Quality Score",labRow=seq(from=10,to=25),col = Lab.palette(256),scale="none" )
dev.off()
"#;

        assert_eq!(output, expected);
    }

    #[test]
    fn empty_case_no_records_pass_filter() {
        let header = test_header();
        let records = Vec::new();
        let bam_records = to_bam_records(&header, &records);
        let hist = compute_quality(bam_records.into_iter().map(Ok), 30).unwrap();

        // Edge case: if NO records pass the mapq filter, read_len stays 0,
        // q_min/q_max stay at their initial sentinel values (93/1)
        assert_eq!(hist.read_len, 0);
        assert_eq!(hist.q_min, 93);
        assert_eq!(hist.q_max, 1);
        assert_eq!(hist.quality.len(), 0);

        // Should produce empty R script
        let output = render_qual_r_script(&hist, 1000, "test_output");
        let expected = r#"pdf('test_output.qual.boxplot.pdf')
boxplot(,xlab="Position of Read(5'->3')",ylab="Phred Quality Score",outline=F)
dev.off()


pdf('test_output.qual.heatmap.pdf')
qual=c()
mat=matrix(qual,ncol=0,byrow=F)
Lab.palette <- colorRampPalette(c("blue", "orange", "red3","red2","red1","red"), space = "rgb",interpolate=c('spline'))
heatmap(mat,Rowv=NA,Colv=NA,xlab="Position of Read",ylab="Phred Quality Score",labRow=seq(from=93,to=1),col = Lab.palette(256),scale="none" )
dev.off()
"#;

        assert_eq!(output, expected);
    }
}