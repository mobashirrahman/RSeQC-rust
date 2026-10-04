//! Port of `read_NVC.py`: calculate nucleotide frequency at each read cycle.
//! Contract: see `compatibility/commands.yaml` entry `read_NVC.py`; algorithm
//! ported from `readsNVC()` in `oracle/upstream-src/src/qcmodule/SAM.py`
//! (lines 2954-3009).
//!
//! Notable upstream behavior preserved exactly: ONLY mapq filtering is applied
//! (no QC-fail/duplicate/secondary/unmapped filtering); reverse-strand reads
//! are reverse-complemented; the output table's row count equals the last
//! processed read's length (not the max length seen across all records).
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

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct NvcTable {
    pub rows: Vec<NvcRow>,
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct NvcRow {
    pub position: usize,
    pub a: u64,
    pub c: u64,
    pub g: u64,
    pub t: u64,
    pub n: u64,
    pub x: u64,
}

fn complement(base: u8) -> u8 {
    // Reuse the same complement table as bam2fq.py
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' => b'N',
        b'X' => b'X',
        other => other,
    }
}

/// Computes nucleotide composition table for a sequence of BAM records,
/// matching the branch order and category semantics of upstream's `readsNVC()`.
/// Pure computation: no I/O, no printing.
pub fn compute_nvc<I>(records: I, q_cut: u8) -> io::Result<NvcTable>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut base_freq = HashMap::new();
    let mut last_read_length = 0;

    for result in records {
        let record = result?;
        
        // NOTE: is_unmapped and is_qcfail checks are commented out in
        // upstream (dead code) -- do NOT filter on them. Only the mapq
        // cutoff applies.
        
        // pysam's raw `.mapq` reports a missing MAPQ as byte value 255;
        // noodles represents "missing" as `None`. 255 always compares >=
        // any realistic q_cut, so treating `None` as 255 here reproduces
        // upstream's comparison behavior without a special case.
        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        
        if mapq < q_cut {
            continue;
        }

        let sequence = record.sequence();
        last_read_length = sequence.len();
        
        let mut seq: Vec<u8> = sequence.iter().collect();
        
        if record.flags().is_reverse_complemented() {
            seq.reverse();
            for b in seq.iter_mut() {
                *b = complement(*b);
            }
        }
        
        for (i, base) in seq.iter().enumerate() {
            // Only count ACGTNX bases as upstream does
            match base {
                b'A' | b'C' | b'G' | b'T' | b'N' | b'X' => {
                    *base_freq.entry((i, *base)).or_insert(0) += 1;
                }
                _ => {
                    // Other bases are not counted in the output, matching upstream
                }
            }
        }
    }

    // Build the table with rows equal to the last read's length
    let mut rows = Vec::with_capacity(last_read_length);
    for position in 0..last_read_length {
        let row = NvcRow {
            position,
            a: *base_freq.get(&(position, b'A')).unwrap_or(&0),
            c: *base_freq.get(&(position, b'C')).unwrap_or(&0),
            g: *base_freq.get(&(position, b'G')).unwrap_or(&0),
            t: *base_freq.get(&(position, b'T')).unwrap_or(&0),
            n: *base_freq.get(&(position, b'N')).unwrap_or(&0),
            x: *base_freq.get(&(position, b'X')).unwrap_or(&0),
        };
        rows.push(row);
    }

    Ok(NvcTable { rows })
}

/// Renders the NVC table as tab-separated text, matching upstream's output format
/// Renders the `.NVC.xls` table. Data rows replicate a literal upstream
/// quirk: the position/A/C/G/T/N fields are each written via
/// `print(value+'\t', end=' ')`, immediately followed by the next
/// `print()` call -- this inserts a literal SPACE after every one of
/// those tabs (same class of `print(..., end=' ')` quirk seen in
/// junction_annotation.py/RPKM_saturation.py). The final field (X) is
/// `print(value+'\t', file=FO)` -- still has its OWN trailing tab (via
/// string concatenation), but no `end=' '`, so the line ends
/// `...\t<X>\t\n` (tab-X-tab-newline), not a bare `<X>\n`. The header
/// row has no such quirk (it's a single literal `print()` call).
pub fn render_nvc_table(table: &NvcTable) -> String {
    let mut out = String::from("Position\tA\tC\tG\tT\tN\tX\n");

    for row in &table.rows {
        out.push_str(&format!("{}\t {}\t {}\t {}\t {}\t {}\t {}\t\n", row.position, row.a, row.c, row.g, row.t, row.n, row.x));
    }

    out
}

/// Renders the `.NVC_plot.r` script text. Ports the R-generation half
/// of `readsNVC()` (lines 3011-3047): each column vector is joined from
/// plain integer counts (no floating-point formatting anywhere -- the
/// `total`/`ym`/`yn`/division expressions are literal R code text, not
/// values computed in Python), and `nx` selects between two nearly-
/// identical branches that differ in which count vectors are included
/// (N/X or not) and the legend's label/color lists.
pub fn render_nvc_r_script(table: &NvcTable, out_prefix: &str, nx: bool) -> String {
    let n = table.rows.len();
    let position: Vec<String> = (0..n).map(|i| i.to_string()).collect();
    let a: Vec<String> = table.rows.iter().map(|r| r.a.to_string()).collect();
    let c: Vec<String> = table.rows.iter().map(|r| r.c.to_string()).collect();
    let g: Vec<String> = table.rows.iter().map(|r| r.g.to_string()).collect();
    let t: Vec<String> = table.rows.iter().map(|r| r.t.to_string()).collect();
    let nn: Vec<String> = table.rows.iter().map(|r| r.n.to_string()).collect();
    let x: Vec<String> = table.rows.iter().map(|r| r.x.to_string()).collect();

    let mut out = String::new();
    out.push_str(&format!("position=c({})\n", position.join(",")));
    out.push_str(&format!("A_count=c({})\n", a.join(",")));
    out.push_str(&format!("C_count=c({})\n", c.join(",")));
    out.push_str(&format!("G_count=c({})\n", g.join(",")));
    out.push_str(&format!("T_count=c({})\n", t.join(",")));
    out.push_str(&format!("N_count=c({})\n", nn.join(",")));
    out.push_str(&format!("X_count=c({})\n", x.join(",")));

    let legend_x = n as i64 - 10;
    let pdf_path = format!("{out_prefix}.NVC_plot.pdf");

    if nx {
        out.push_str("total= A_count + C_count + G_count + T_count + N_count + X_count\n");
        out.push_str("ym=max(A_count/total,C_count/total,G_count/total,T_count/total,N_count/total,X_count/total) + 0.05\n");
        out.push_str("yn=min(A_count/total,C_count/total,G_count/total,T_count/total,N_count/total,X_count/total)\n");
        out.push_str(&format!("pdf(\"{pdf_path}\")\n"));
        out.push_str("plot(position,A_count/total,type=\"o\",pch=20,ylim=c(yn,ym),col=\"dark green\",xlab=\"Position of Read\",ylab=\"Nucleotide Frequency\")\n");
        out.push_str("lines(position,T_count/total,type=\"o\",pch=20,col=\"red\")\n");
        out.push_str("lines(position,G_count/total,type=\"o\",pch=20,col=\"blue\")\n");
        out.push_str("lines(position,C_count/total,type=\"o\",pch=20,col=\"cyan\")\n");
        out.push_str("lines(position,N_count/total,type=\"o\",pch=20,col=\"black\")\n");
        out.push_str("lines(position,X_count/total,type=\"o\",pch=20,col=\"grey\")\n");
        out.push_str(&format!(
            "legend({legend_x},ym,legend=c(\"A\",\"T\",\"G\",\"C\",\"N\",\"X\"),col=c(\"dark green\",\"red\",\"blue\",\"cyan\",\"black\",\"grey\"),lwd=2,pch=20,text.col=c(\"dark green\",\"red\",\"blue\",\"cyan\",\"black\",\"grey\"))\n"
        ));
        out.push_str("dev.off()\n");
    } else {
        out.push_str("total= A_count + C_count + G_count + T_count\n");
        out.push_str("ym=max(A_count/total,C_count/total,G_count/total,T_count/total) + 0.05\n");
        out.push_str("yn=min(A_count/total,C_count/total,G_count/total,T_count/total)\n");
        out.push_str(&format!("pdf(\"{pdf_path}\")\n"));
        out.push_str("plot(position,A_count/total,type=\"o\",pch=20,ylim=c(yn,ym),col=\"dark green\",xlab=\"Position of Read\",ylab=\"Nucleotide Frequency\")\n");
        out.push_str("lines(position,T_count/total,type=\"o\",pch=20,col=\"red\")\n");
        out.push_str("lines(position,G_count/total,type=\"o\",pch=20,col=\"blue\")\n");
        out.push_str("lines(position,C_count/total,type=\"o\",pch=20,col=\"cyan\")\n");
        out.push_str(&format!(
            "legend({legend_x},ym,legend=c(\"A\",\"T\",\"G\",\"C\"),col=c(\"dark green\",\"red\",\"blue\",\"cyan\"),lwd=2,pch=20,text.col=c(\"dark green\",\"red\",\"blue\",\"cyan\"))\n"
        ));
        out.push_str("dev.off()\n");
    }

    out
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

/// Runs the `read_NVC.py` CLI body over an already-opened record stream
/// (C1 multi-driver pattern).
///
/// This is exactly what the standalone binary's `run()` does -- same
/// progress lines, same files with the same bytes, same Rscript contract,
/// same error propagation -- except the record source is a
/// caller-supplied iterator and stdout/stderr are caller-supplied sinks.
/// The standalone binary delegates to this (passing the process streams);
/// `rseqc_multi` passes one record broadcast plus per-command stream
/// files. `compute_nvc` and both renderers are untouched. See
/// `run_read_gc` for the transport notes (parent-check `Result` form,
/// piped Rscript child), which apply here unchanged.
///
/// The argument list is long on purpose and carries a targeted lint
/// allowance: this is the C1 multi-driver pattern (one callable per
/// command carrying its full CLI surface plus the two sinks), and
/// bundling the flags into a struct would only hide them from the C2
/// cards that copy this signature command by command.
#[allow(clippy::too_many_arguments)]
pub fn run_read_nvc<I>(
    records: I,
    q_cut: u8,
    out_prefix: &str,
    nx: bool,
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
    // SAM too; the "Read SAM file" branch is practically dead code for
    // any valid input). The literal's own trailing space plus `end=' '`
    // gives two spaces before "Done".
    write!(stderr, "Read BAM file ...  ")?;
    let table = compute_nvc(records, q_cut)?;
    writeln!(stderr, "Done")?;

    writeln!(stderr, "generating data matrix ...")?;
    let nvc_path = format!("{out_prefix}.NVC.xls");
    File::create(&nvc_path)?.write_all(render_nvc_table(&table).as_bytes())?;

    // Upstream: `print("generating R script  ...", ...)` -- literal has
    // two spaces between "script" and "...".
    writeln!(stderr, "generating R script  ...")?;
    let r_path = format!("{out_prefix}.NVC_plot.r");
    File::create(&r_path)?.write_all(render_nvc_r_script(&table, out_prefix, nx).as_bytes())?;

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
    /// test exercises the same decode path `compute_nvc` runs against in
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

        // Forward read: ACGT (length 4)
        let forward = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"ACGT".to_vec()))
            .build();

        // Reverse read: ACGT -> reverse complement is ACGT (palindromic)
        let reverse = RecordBuf::builder()
            .set_flags(Flags::REVERSE_COMPLEMENTED)
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"ACGT".to_vec()))
            .build();

        // Short read: AT (length 2) - this should determine the table length
        let short = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 2)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"AT".to_vec()))
            .build();

        // Low MAPQ read (should be filtered out)
        let low_mapq = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(10).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 3)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"ACG".to_vec()))
            .build();

        let records = vec![forward, reverse, short, low_mapq];
        let bam_records = to_bam_records(&header, &records);
        let table = compute_nvc(bam_records.into_iter().map(Ok), 30).unwrap();

        // Hand-computed from the fixture above (all three mapq>=30 reads
        // contribute at positions 0 and 1; low_mapq is filtered out):
        // - forward "ACGT": 0:A, 1:C, 2:G, 3:T
        // - reverse "ACGT" -> reverse-complemented -> "ACGT" (palindromic): 0:A, 1:C, 2:G, 3:T
        // - short "AT": 0:A, 1:T
        // Position 0: A = forward's A + reverse's A + short's A = 3
        // Position 1: C = forward's C + reverse's C = 2; T = short's T = 1
        // Table length is 2 (the LAST record processed is `short`, length 2),
        // matching the preserved upstream quirk (DIV-0008) even though the
        // forward/reverse reads were longer.

        assert_eq!(table.rows.len(), 2, "Table length should match last read's length");

        assert_eq!(table.rows[0].position, 0);
        assert_eq!(table.rows[0].a, 3);
        assert_eq!(table.rows[0].c, 0);
        assert_eq!(table.rows[0].g, 0);
        assert_eq!(table.rows[0].t, 0);
        assert_eq!(table.rows[0].n, 0);
        assert_eq!(table.rows[0].x, 0);

        assert_eq!(table.rows[1].position, 1);
        assert_eq!(table.rows[1].a, 0);
        assert_eq!(table.rows[1].c, 2);
        assert_eq!(table.rows[1].g, 0);
        assert_eq!(table.rows[1].t, 1);
        assert_eq!(table.rows[1].n, 0);
        assert_eq!(table.rows[1].x, 0);
    }

    #[test]
    fn render_output_format() {
        let table = NvcTable {
            rows: vec![
                NvcRow { position: 0, a: 10, c: 5, g: 3, t: 8, n: 1, x: 0 },
                NvcRow { position: 1, a: 8, c: 12, g: 6, t: 4, n: 0, x: 1 },
            ],
        };

        let output = render_nvc_table(&table);
        // Cross-checked byte-for-byte against a `python3 -c` run of the
        // literal upstream print(...)/end=' ' sequence with this data.
        let expected = "Position\tA\tC\tG\tT\tN\tX\n0\t 10\t 5\t 3\t 8\t 1\t 0\t\n1\t 8\t 12\t 6\t 4\t 0\t 1\t\n";

        assert_eq!(output, expected);
    }

    #[test]
    fn render_nvc_r_script_without_nx_exact_text() {
        // Cross-checked byte-for-byte against a real python3 -c run of
        // the literal upstream print()/string-join sequence with this
        // exact data (readsNVC lines 3011-3047, nx=False branch).
        let table = NvcTable {
            rows: vec![
                NvcRow { position: 0, a: 3, c: 0, g: 0, t: 0, n: 0, x: 0 },
                NvcRow { position: 1, a: 0, c: 2, g: 0, t: 1, n: 0, x: 0 },
            ],
        };
        let script = render_nvc_r_script(&table, "out", false);
        let expected = "\
position=c(0,1)
A_count=c(3,0)
C_count=c(0,2)
G_count=c(0,0)
T_count=c(0,1)
N_count=c(0,0)
X_count=c(0,0)
total= A_count + C_count + G_count + T_count
ym=max(A_count/total,C_count/total,G_count/total,T_count/total) + 0.05
yn=min(A_count/total,C_count/total,G_count/total,T_count/total)
pdf(\"out.NVC_plot.pdf\")
plot(position,A_count/total,type=\"o\",pch=20,ylim=c(yn,ym),col=\"dark green\",xlab=\"Position of Read\",ylab=\"Nucleotide Frequency\")
lines(position,T_count/total,type=\"o\",pch=20,col=\"red\")
lines(position,G_count/total,type=\"o\",pch=20,col=\"blue\")
lines(position,C_count/total,type=\"o\",pch=20,col=\"cyan\")
legend(-8,ym,legend=c(\"A\",\"T\",\"G\",\"C\"),col=c(\"dark green\",\"red\",\"blue\",\"cyan\"),lwd=2,pch=20,text.col=c(\"dark green\",\"red\",\"blue\",\"cyan\"))
dev.off()
";
        assert_eq!(script, expected);
    }

    #[test]
    fn render_nvc_r_script_with_nx_includes_n_and_x_series() {
        let table = NvcTable { rows: vec![NvcRow { position: 0, a: 1, c: 2, g: 3, t: 4, n: 5, x: 6 }] };
        let script = render_nvc_r_script(&table, "out", true);
        assert!(script.contains("N_count=c(5)\n"));
        assert!(script.contains("X_count=c(6)\n"));
        assert!(script.contains("total= A_count + C_count + G_count + T_count + N_count + X_count\n"));
        assert!(script.contains("lines(position,N_count/total,type=\"o\",pch=20,col=\"black\")\n"));
        assert!(script.contains("lines(position,X_count/total,type=\"o\",pch=20,col=\"grey\")\n"));
        assert!(script.contains("legend=c(\"A\",\"T\",\"G\",\"C\",\"N\",\"X\")"));
    }
}