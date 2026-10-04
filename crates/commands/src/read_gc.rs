//! Port of `read_GC.py`: calculate GC content distribution of aligned reads.
//! Contract: see `compatibility/commands.yaml` entry `read_GC.py`; algorithm
//! ported from `readGC()` in `oracle/upstream-src/src/qcmodule/SAM.py`
//! (lines 3119-3158).
//!
//! Notable upstream behavior preserved exactly: unmapped reads, QC-fail reads,
//! and low MAPQ reads are filtered out; GC percent is calculated as a string
//! with 2 decimal places; insertion order of distinct GC percentages is preserved.
//!
//! SAM-text and CRAM input (DIV-0002/0004) are supported at the CLI layer via
//! `rseqc_formats::open_alignments`; this module's own functions were
//! unaffected (already generic over
//! `IntoIterator<Item = io::Result<bam::Record>>`).

use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::Write as _;
use std::process::{Command, Stdio};

use noodles_bam as bam;

use crate::exec_resolve;

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

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct GcHistogram {
    pub buckets: Vec<(String, u64)>, // (gc_percent_string, count) in first-seen order
}

/// Computes GC content histogram for a sequence of BAM records,
/// matching the upstream algorithm exactly.
/// Pure computation: no I/O, no printing.
pub fn compute_gc<I>(records: I, q_cut: u8) -> io::Result<GcHistogram>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut gc_hist: HashMap<String, usize> = HashMap::new(); // Maps gc_percent_string -> index in buckets vec
    let mut buckets: Vec<(String, u64)> = Vec::new(); // Maintains insertion order

    for result in records {
        let record = result?;
        
        // Filter out unmapped reads
        if record.flags().is_unmapped() {
            continue;
        }
        
        // Filter out QC fail reads
        if record.flags().is_qc_fail() {
            continue;
        }
        
        // Filter out low MAPQ reads
        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        if mapq < q_cut {
            continue;
        }

        // Get uppercase sequence
        let sequence = record.sequence();
        let rna_read: String = sequence.iter().map(|b| b as char).collect();
        
        // Calculate GC percent: (C + G) / total_length * 100, formatted as "%.2f"
        let c_count = rna_read.chars().filter(|&c| c == 'C').count();
        let g_count = rna_read.chars().filter(|&c| c == 'G').count();
        let total_length = rna_read.len();
        
        // Avoid division by zero (though upstream doesn't handle this case)
        if total_length == 0 {
            continue;
        }
        
        let gc_percent = (c_count + g_count) as f64 / total_length as f64 * 100.0;
        let gc_percent_str = format!("{:.2}", gc_percent);
        
        // Update histogram while preserving insertion order
        if gc_hist.contains_key(&gc_percent_str) {
            // Already seen this GC percentage, increment count
            let index = gc_hist[&gc_percent_str];
            buckets[index].1 += 1;
        } else {
            // New GC percentage, add to both HashMap and Vec
            let index = buckets.len();
            gc_hist.insert(gc_percent_str.clone(), index);
            buckets.push((gc_percent_str, 1));
        }
    }

    Ok(GcHistogram { buckets })
}

/// Renders the GC histogram as tab-separated text, matching upstream's
/// `.GC.xls` format. Every line, including the last, ends with `\n`
/// (upstream's plain `print(...)` calls each add their own trailing
/// newline).
pub fn render_gc_table(hist: &GcHistogram) -> String {
    let mut out = String::from("GC%\tread_count\n");
    for (gc_percent, count) in &hist.buckets {
        out.push_str(&format!("{gc_percent}\t{count}\n"));
    }
    out
}

/// Renders the GC histogram as an R script, matching upstream's
/// `.GC_plot.r` format (see `render_gc_table` for the trailing-newline
/// note -- applies here too, including after the final `dev.off()`).
pub fn render_gc_r_script(hist: &GcHistogram, output_prefix: &str) -> String {
    let gc_values: Vec<String> = hist.buckets.iter().map(|(gc, _)| gc.clone()).collect();
    let count_values: Vec<String> = hist.buckets.iter().map(|(_, count)| count.to_string()).collect();

    format!(
        "pdf(\"{output_prefix}.GC_plot.pdf\")\ngc=rep(c({}),times=c({}))\nhist(gc,probability=T,breaks=100,xlab=\"GC content (%)\",ylab=\"Density of Reads\",border=\"blue\",main=\"\")\ndev.off()\n",
        gc_values.join(","),
        count_values.join(","),
    )
}

/// Runs the `read_GC.py` CLI body over an already-opened record stream
/// (C1 multi-driver pattern).
///
/// This is exactly what the standalone binary's `run()` does -- same
/// progress lines, same files with the same bytes, same Rscript contract,
/// same error propagation -- except the record source is a
/// caller-supplied iterator and stdout/stderr are caller-supplied sinks.
/// The standalone binary delegates to this (passing the process streams);
/// `rseqc_multi` passes one record broadcast plus per-command stream
/// files. `compute_gc` and both renderers are untouched.
///
/// Two deliberate transport notes for the multi driver. First, the output
/// parent check uses the `Result` form (`require_existing_output_parent`)
/// rather than the binary's `_or_exit` form: the binary keeps its own
/// `_or_exit` call first (so standalone exit 2 is unchanged) while the
/// multi worker maps this `Err` to the same `prog: error:` line in its
/// stream file. Second, the Rscript child is spawned piped (not
/// inherited) and its captured stdout/stderr are copied to these sinks:
/// per-stream bytes are identical to inheritance (the parent emits nothing
/// between spawn and wait), and piping is what lets a multi worker
/// attribute the child's output to its own stream files.
pub fn run_read_gc<I>(
    records: I,
    q_cut: u8,
    out_prefix: &str,
    skip_plot: bool,
    rscript: &str,
    stdout: &mut dyn io::Write,
    stderr: &mut dyn io::Write,
) -> io::Result<()>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    // Same rule as the binary's `require_existing_output_parent_or_exit`
    // (and upstream's `validate_args`), in `Result` form so a multi worker
    // can attribute the failure to its own stream file instead of exiting
    // the whole driver process. Message text matches exactly.
    require_output_parent(out_prefix)?;

    // Upstream: `if self.bam_format: print("Read BAM file ... ", end=' ')
    // else: print("Read SAM file ... ", end=' ')` -- `self.bam_format`
    // comes from `pysam.Samfile(path, 'rb')` succeeding, which it does
    // even for genuine plain-text SAM content (htslib auto-detects,
    // ignoring the 'b' mode hint; confirmed via a live diff for
    // bam_stat.py/read_NVC.py, same underlying pysam.Samfile call here).
    // The "Read SAM file" branch is practically dead code for any valid
    // input. The literal's own trailing space plus `end=' '` gives two
    // spaces before "Done".
    write!(stderr, "Read BAM file ...  ")?;
    let hist = compute_gc(records, q_cut)?;
    writeln!(stderr, "Done")?;

    writeln!(stderr, "writing GC content ...")?;
    let xls_path = format!("{out_prefix}.GC.xls");
    File::create(&xls_path)?.write_all(render_gc_table(&hist).as_bytes())?;

    writeln!(stderr, "writing R script ...")?;
    let r_path = format!("{out_prefix}.GC_plot.r");
    File::create(&r_path)?.write_all(render_gc_r_script(&hist, out_prefix).as_bytes())?;

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
    /// test exercises the same decode path `compute_gc` runs against in
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
    fn matches_hand_computed_counts() {
        let header = test_header();

        // Read with 50% GC: "ATGC" (2 C+G out of 4 = 50.00%)
        let read_50_percent = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"ATGC".to_vec()))
            .build();

        // Read with 25% GC: "ATCG" (1 C+G out of 4 = 25.00%)
        let read_25_percent = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"ATCG".to_vec()))
            .build();

        // Read with 75% GC: "GGCT" (3 C+G out of 4 = 75.00%)
        let read_75_percent = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"GGCT".to_vec()))
            .build();

        // Another read with 50% GC: "CGTA" (2 C+G out of 4 = 50.00%)
        let read_50_percent_2 = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"CGTA".to_vec()))
            .build();

        let records = vec![read_50_percent, read_25_percent, read_75_percent, read_50_percent_2];
        let bam_records = to_bam_records(&header, &records);
        let hist = compute_gc(bam_records.into_iter().map(Ok), 30).unwrap();

        // Expected: 50.00% (3 reads), 75.00% (1 read)
        // Order should be: 50.00, 75.00 (first-seen order from the test data)
        assert_eq!(hist.buckets.len(), 2);
        
        // Check first-seen order: 50% appears first, then 75%
        assert_eq!(hist.buckets[0], ("50.00".to_string(), 3));
        assert_eq!(hist.buckets[1], ("75.00".to_string(), 1));
    }

    #[test]
    fn render_gc_table_format() {
        let hist = GcHistogram {
            buckets: vec![
                ("25.00".to_string(), 10),
                ("50.00".to_string(), 25),
                ("75.00".to_string(), 5),
            ],
        };

        let output = render_gc_table(&hist);
        let expected = "GC%\tread_count\n25.00\t10\n50.00\t25\n75.00\t5\n";

        assert_eq!(output, expected);
    }

    #[test]
    fn render_gc_r_script_format() {
        let hist = GcHistogram {
            buckets: vec![
                ("25.00".to_string(), 10),
                ("50.00".to_string(), 25),
                ("75.00".to_string(), 5),
            ],
        };

        let output = render_gc_r_script(&hist, "test_output");
        let expected = "pdf(\"test_output.GC_plot.pdf\")\ngc=rep(c(25.00,50.00,75.00),times=c(10,25,5))\nhist(gc,probability=T,breaks=100,xlab=\"GC content (%)\",ylab=\"Density of Reads\",border=\"blue\",main=\"\")\ndev.off()\n";

        assert_eq!(output, expected);
    }

    #[test]
    fn gc_percent_formatting() {
        // Test edge cases in GC percent formatting
        let header = test_header();

        // Read with 0% GC: "AAAA" (0 C+G out of 4 = 0.00%)
        let read_0_percent = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"AAAA".to_vec()))
            .build();

        // Read with 100% GC: "CCCC" (4 C+G out of 4 = 100.00%)
        let read_100_percent = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"CCCC".to_vec()))
            .build();

        // Read with 33.333...% GC: "ACT" (1 C+G out of 3 ≈ 33.33%)
        let read_33_percent = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 3)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"ACT".to_vec()))
            .build();

        let records = vec![read_0_percent, read_100_percent, read_33_percent];
        let bam_records = to_bam_records(&header, &records);
        let hist = compute_gc(bam_records.into_iter().map(Ok), 30).unwrap();

        // Check exact formatting matches Python's "%.2f" behavior
        // Order should be: 0.00, 100.00, 33.33 (first-seen order from the test data)
        assert_eq!(hist.buckets[0], ("0.00".to_string(), 1));
        assert_eq!(hist.buckets[1], ("100.00".to_string(), 1));
        assert_eq!(hist.buckets[2], ("33.33".to_string(), 1));
    }
}