//! Port of `insertion_profile.py`: calculate the distribution of inserted
//! nucleotides across reads. Contract: see `compatibility/commands.yaml`
//! entry `insertion_profile.py`; algorithm ported from
//! `insertion_profile()` in `oracle/upstream-src/src/qcmodule/SAM.py`
//! (lines 3362-3481+), which is structurally identical to
//! `clipping_profile()` with `type="I"` instead of `type="S"` and
//! different labels/file names.
//!
//! Reuses [`crate::clipping_profile::compute_single_end`] and
//! [`crate::clipping_profile::compute_paired_end`] directly (called with
//! `clip_char = b'I'`) rather than duplicating the aggregation logic —
//! only the rendering (column headers, R variable names, plot titles) and
//! file names differ from clipping_profile.py. Same DIV-0010-pattern
//! last-record-length quirk applies here too (inherited from the shared
//! compute functions).

use std::fs::File;
use std::io;
use std::io::Write as _;
use std::process::{Command, Stdio};

use noodles_bam as bam;

use crate::clipping_profile::{
    PairedEndProfile, SingleEndProfile, compute_paired_end, compute_single_end,
};
use crate::exec_resolve;

fn fmt_float(n: u64) -> String {
    format!("{n}.0")
}

/// Same `defaultdict(int)` duck-typing quirk as
/// `clipping_profile::fmt_clip_count`: an untouched position prints bare
/// `"0"`, a touched one prints `"N.0"`. Applies to Insert_nt/insert_count,
/// NOT to Non_insert_nt/noninsert_count (always float, see that function's
/// doc comment for the full explanation).
fn fmt_clip_count(n: u64) -> String {
    if n == 0 { "0".to_string() } else { format!("{n}.0") }
}

/// Every line, including the last, ends with `\n` (upstream's plain
/// `print(...)` calls each add their own trailing newline) -- applies
/// to all four render functions in this module, same as
/// clipping_profile.rs (structurally identical upstream function).
pub fn render_single_table(p: &SingleEndProfile) -> String {
    let mut out = String::from("Position\tInsert_nt\tNon_insert_nt\n");
    for (i, &c) in p.clip_count.iter().enumerate() {
        let non_insert = p.total_read - c;
        out.push_str(&format!("{i}\t{}\t{}\n", fmt_clip_count(c), fmt_float(non_insert)));
    }
    out
}

pub fn render_single_r_script(p: &SingleEndProfile, out_prefix: &str) -> String {
    let read_pos: Vec<String> = (0..p.clip_count.len()).map(|i| i.to_string()).collect();
    let insert_strs: Vec<String> = p.clip_count.iter().map(|&c| fmt_clip_count(c)).collect();

    format!(
        "pdf(\"{out_prefix}.insertion_profile.pdf\")\nread_pos=c({})\ninsert_count=c({})\nnoninsert_count= {} - insert_count\nplot(read_pos, insert_count*100/(insert_count+noninsert_count),col=\"blue\",main=\"Insertion profile\",xlab=\"Position of read\",ylab=\"Insertion %\",type=\"b\")\ndev.off()\n",
        read_pos.join(","),
        insert_strs.join(","),
        p.total_read,
    )
}

pub fn render_paired_table(p: &PairedEndProfile) -> String {
    let mut out = String::from("Position\tInsert_nt\tNon_insert_nt\nRead-1:\n");
    for (i, &c) in p.r1_clip_count.iter().enumerate() {
        let non_insert = p.total_read1 - c;
        out.push_str(&format!("{i}\t{}\t{}\n", fmt_clip_count(c), fmt_float(non_insert)));
    }
    out.push_str("Read-2:\n");
    for (i, &c) in p.r2_clip_count.iter().enumerate() {
        let non_insert = p.total_read2 - c;
        out.push_str(&format!("{i}\t{}\t{}\n", fmt_clip_count(c), fmt_float(non_insert)));
    }
    out
}

pub fn render_paired_r_script(p: &PairedEndProfile, out_prefix: &str) -> String {
    let read_pos: Vec<String> = (0..p.r1_clip_count.len()).map(|i| i.to_string()).collect();
    let r1_strs: Vec<String> = p.r1_clip_count.iter().map(|&c| fmt_clip_count(c)).collect();
    let r2_strs: Vec<String> = p.r2_clip_count.iter().map(|&c| fmt_clip_count(c)).collect();
    let read_pos_csv = read_pos.join(",");

    format!(
        "pdf(\"{out_prefix}.insertion_profile.R1.pdf\")\nread_pos=c({read_pos_csv})\nr1_insert_count=c({})\nr1_noninsert_count = {} - r1_insert_count\nplot(read_pos, r1_insert_count*100/(r1_insert_count + r1_noninsert_count),col=\"blue\",main=\"Insertion profile\",xlab=\"Position of read (read-1)\",ylab=\"Insertion %\",type=\"b\")\ndev.off()\n\
pdf(\"{out_prefix}.insertion_profile.R2.pdf\")\nread_pos=c({read_pos_csv})\nr2_insert_count=c({})\nr2_noninsert_count = {} - r2_insert_count\nplot(read_pos, r2_insert_count*100/(r2_insert_count + r2_noninsert_count),col=\"blue\",main=\"Insertion profile\",xlab=\"Position of read (read-2)\",ylab=\"Insertion %\",type=\"b\")\ndev.off()\n",
        r1_strs.join(","),
        p.total_read1,
        r2_strs.join(","),
        p.total_read2,
    )
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

/// Runs the `insertion_profile.py` CLI body over an already-opened record
/// stream (C1 multi-driver pattern).
///
/// This is exactly what the standalone binary's `run()` does -- same
/// progress lines (including the upstream `Totoal` typo), same files with
/// the same bytes, same Rscript contract, same error propagation --
/// except the record source is a caller-supplied iterator and
/// stdout/stderr are caller-supplied sinks. The standalone binary
/// delegates to this (passing the process streams); `rseqc_multi` passes
/// one record broadcast plus per-command stream files. The shared compute
/// functions and all four renderers are untouched.
///
/// `sequencing` is `"SE"` or `"PE"`, the only values the standalone
/// binary's `--sequencing` accepts (clap enforces it there); anything else
/// is an error, which keeps the driver's validation honest rather than
/// silently picking a branch.
///
/// The clip character is `b'I'` (upstream's `type="I"`), distinguishing
/// this from `clipping_profile`'s `b'S'` while sharing the aggregation.
///
/// The argument list is long on purpose and carries a targeted lint
/// allowance: this is the C1 multi-driver pattern (one callable per
/// command carrying its full CLI surface plus the two sinks), and
/// bundling the flags into a struct would only hide them from the C2
/// cards that copy this signature command by command.
#[allow(clippy::too_many_arguments)]
pub fn run_insertion_profile<I>(
    records: I,
    q_cut: u8,
    out_prefix: &str,
    sequencing: &str,
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

    // Upstream: `if self.bam_format: print("Load BAM file ... ", end=' ')`
    // -- always the BAM branch in practice; the literal's own trailing
    // space plus `end=' '` gives two spaces before "Done".
    write!(stderr, "Load BAM file ...  ")?;

    let (table_text, r_script_text) = match sequencing {
        "SE" => {
            let profile = compute_single_end(records, q_cut, b'I')?;
            writeln!(stderr, "Done")?;
            // Upstream: `print("Totoal reads used: %d" % ...)` -- a
            // literal upstream typo ("Totoal"), preserved exactly.
            writeln!(stderr, "Totoal reads used: {}", profile.total_read)?;
            (
                render_single_table(&profile),
                render_single_r_script(&profile, out_prefix),
            )
        }
        "PE" => {
            let profile = compute_paired_end(records, q_cut, b'I')?;
            writeln!(stderr, "Done")?;
            // Upstream prints these as TWO SEPARATE lines (also with
            // the same "Totoal" typo), not one combined line.
            writeln!(stderr, "Totoal read-1 used: {}", profile.total_read1)?;
            writeln!(stderr, "Totoal read-2 used: {}", profile.total_read2)?;
            (
                render_paired_table(&profile),
                render_paired_r_script(&profile, out_prefix),
            )
        }
        other => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("invalid sequencing layout: {other} (expected SE or PE)"),
            ));
        }
    };

    let xls_path = format!("{out_prefix}.insertion_profile.xls");
    File::create(&xls_path)?.write_all(table_text.as_bytes())?;

    let r_path = format!("{out_prefix}.insertion_profile.r");
    File::create(&r_path)?.write_all(r_script_text.as_bytes())?;

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

    #[test]
    fn render_single_table_uses_insert_labels() {
        let profile = SingleEndProfile {
            total_read: 4,
            clip_count: vec![1, 3],
        };
        assert_eq!(
            render_single_table(&profile),
            "Position\tInsert_nt\tNon_insert_nt\n0\t1.0\t3.0\n1\t3.0\t1.0\n"
        );
    }

    #[test]
    fn render_single_table_untouched_position_prints_bare_zero() {
        // Same defaultdict(int)-vs-float duck-typing quirk as
        // clipping_profile.rs: Insert_nt at an untouched position is "0",
        // not "0.0"; Non_insert_nt is always float. Verified via python3 -c.
        let profile = SingleEndProfile {
            total_read: 3,
            clip_count: vec![0, 3, 0],
        };
        assert_eq!(
            render_single_table(&profile),
            "Position\tInsert_nt\tNon_insert_nt\n0\t0\t3.0\n1\t3.0\t0.0\n2\t0\t3.0\n"
        );
    }

    #[test]
    fn render_single_r_script_exact_text() {
        // Independently derived by running the equivalent Python
        // print()/%-format expressions from oracle/upstream-src/src/
        // qcmodule/SAM.py lines 3410-3415 via `python3 -c`.
        let profile = SingleEndProfile {
            total_read: 4,
            clip_count: vec![1, 3],
        };
        let output = render_single_r_script(&profile, "test_output");
        let expected = "pdf(\"test_output.insertion_profile.pdf\")\n\
read_pos=c(0,1)\n\
insert_count=c(1.0,3.0)\n\
noninsert_count= 4 - insert_count\n\
plot(read_pos, insert_count*100/(insert_count+noninsert_count),col=\"blue\",main=\"Insertion profile\",xlab=\"Position of read\",ylab=\"Insertion %\",type=\"b\")\n\
dev.off()\n";
        assert_eq!(output, expected);
    }

    #[test]
    fn render_paired_r_script_exact_text() {
        // Independently derived from oracle/upstream-src/src/qcmodule/
        // SAM.py lines 3471-3481 via `python3 -c`.
        let profile = PairedEndProfile {
            total_read1: 5,
            total_read2: 4,
            r1_clip_count: vec![2, 0],
            r2_clip_count: vec![0, 1],
        };
        let output = render_paired_r_script(&profile, "test_output");
        let expected = "pdf(\"test_output.insertion_profile.R1.pdf\")\n\
read_pos=c(0,1)\n\
r1_insert_count=c(2.0,0)\n\
r1_noninsert_count = 5 - r1_insert_count\n\
plot(read_pos, r1_insert_count*100/(r1_insert_count + r1_noninsert_count),col=\"blue\",main=\"Insertion profile\",xlab=\"Position of read (read-1)\",ylab=\"Insertion %\",type=\"b\")\n\
dev.off()\n\
pdf(\"test_output.insertion_profile.R2.pdf\")\n\
read_pos=c(0,1)\n\
r2_insert_count=c(0,1.0)\n\
r2_noninsert_count = 4 - r2_insert_count\n\
plot(read_pos, r2_insert_count*100/(r2_insert_count + r2_noninsert_count),col=\"blue\",main=\"Insertion profile\",xlab=\"Position of read (read-2)\",ylab=\"Insertion %\",type=\"b\")\n\
dev.off()\n";
        assert_eq!(output, expected);
    }
}
