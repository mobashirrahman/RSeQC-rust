//! Port of `deletion_profile.py`: calculate the distribution of deleted
//! nucleotides across aligned reads. Contract: see
//! `compatibility/commands.yaml` entry `deletion_profile.py`; algorithm
//! ported from `deletionProfile()` in `oracle/upstream-src/src/qcmodule/
//! SAM.py` (lines 4554-4634).
//!
//! Unlike the other `*_profile.py` commands ported so far, `read_length`
//! is a required, fixed CLI parameter here (`-l/--read-align-length`), not
//! derived from the last-processed record -- so there is no DIV-0008-style
//! quirk in this command. Reads whose actual sequence length or "matched
//! portion" (sum of M/S/I op lengths) doesn't exactly equal `read_length`
//! are skipped entirely, as are reads with no deletion at all. `read_num`
//! caps the count of *qualifying* reads processed (checked before pulling
//! the next record), not the total records scanned.

use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::Write as _;
use std::process::{Command, Stdio};

use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::op::Kind;
use rseqc_formats::cigar::fetch_deletion_range;

use crate::exec_resolve;

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct DeletionProfile {
    pub count: u64,
    pub del_counts: Vec<u64>,
}

pub fn compute_deletion_profile<I>(
    records: I,
    q_cut: u8,
    read_length: usize,
    read_num: u64,
) -> io::Result<DeletionProfile>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut count = 0u64;
    let mut del_postns: HashMap<usize, u64> = HashMap::new();

    for result in records {
        if count >= read_num {
            break;
        }
        let record = result?;

        let flags = record.flags();
        if flags.is_qc_fail() || flags.is_duplicate() || flags.is_secondary() || flags.is_unmapped() {
            continue;
        }
        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        if mapq < q_cut {
            continue;
        }

        let ops: Vec<_> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;

        let has_deletion = ops.iter().any(|op| op.kind() == Kind::Deletion);
        if !has_deletion {
            continue;
        }

        if record.sequence().len() != read_length {
            continue;
        }

        let matched_portion: usize = ops
            .iter()
            .filter(|op| matches!(op.kind(), Kind::Match | Kind::SoftClip | Kind::Insertion))
            .map(|op| op.len())
            .sum();
        if matched_portion != read_length {
            continue;
        }

        count += 1;
        for (p, _size) in fetch_deletion_range(ops.iter().copied()) {
            let p = if flags.is_reverse_complemented() {
                read_length - p
            } else {
                p
            };
            // A flipped position can land exactly on `read_length` (when
            // the original position was 0), which is out of the
            // `0..read_length` output range and so is silently absent
            // from the rendered table -- matches upstream, which never
            // guards against this either.
            *del_postns.entry(p).or_insert(0) += 1;
        }
    }

    let del_counts = (0..read_length).map(|k| *del_postns.get(&k).unwrap_or(&0)).collect();
    Ok(DeletionProfile { count, del_counts })
}

/// Every line, including the last, ends with `\n` (upstream's plain
/// `print(...)` calls each add their own trailing newline).
pub fn render_deletion_table(p: &DeletionProfile) -> String {
    let mut out = String::from("read_position\tdeletion_count\n");
    for (k, &c) in p.del_counts.iter().enumerate() {
        out.push_str(&format!("{k}\t{c}\n"));
    }
    out
}

/// See `render_deletion_table` for the trailing-newline note (applies
/// here too, including after the final `dev.off()`).
pub fn render_deletion_r_script(p: &DeletionProfile, out_prefix: &str) -> String {
    let pos: Vec<String> = (0..p.del_counts.len()).map(|i| i.to_string()).collect();
    let vals: Vec<String> = p.del_counts.iter().map(|c| c.to_string()).collect();

    format!(
        "pdf(\"{out_prefix}.deletion_profile.pdf\")\npos=c({})\nvalue=c({})\nplot(pos,value,type='b', col='blue',xlab=\"Read position (5'->3')\", ylab='Deletion count')\ndev.off()\n",
        pos.join(","),
        vals.join(","),
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

/// Runs the `deletion_profile.py` CLI body over an already-opened record
/// stream (C1 multi-driver pattern).
///
/// This is exactly what the standalone binary's `run()` does -- same
/// progress lines, same two stdout blank lines (upstream's
/// `print('\n')` twice: the literal `\n` plus print's own newline), same
/// files with the same bytes, same Rscript contract, same error
/// propagation -- except the record source is a caller-supplied iterator
/// and stdout/stderr are caller-supplied sinks. The standalone binary
/// delegates to this (passing the process streams); `rseqc_multi` passes
/// one record broadcast plus per-command stream files. The compute
/// function and both renderers are untouched.
///
/// `read_align_length` and `read_num` are the command's own `-l`/`-n`
/// flags; both are forwarded verbatim because the filtering they drive is
/// the whole point of the command, and a driver default would silently
/// change which reads qualify (the driver makes `-l` required for this
/// command rather than inventing a length).
///
/// The argument list is long on purpose and carries a targeted lint
/// allowance: this is the C1 multi-driver pattern (one callable per
/// command carrying its full CLI surface plus the two sinks), and
/// bundling the flags into a struct would only hide them from the C2
/// cards that copy this signature command by command.
#[allow(clippy::too_many_arguments)]
pub fn run_deletion_profile<I>(
    records: I,
    q_cut: u8,
    out_prefix: &str,
    read_align_length: usize,
    read_num: u64,
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

    // Upstream: `print("Process BAM file ... ", end=' ', file=sys.stderr)`
    // -- the string literal's own trailing space plus `end=' '` gives two
    // spaces before "Total reads used" on the same stderr line.
    write!(stderr, "Process BAM file ...  ")?;
    let profile = compute_deletion_profile(records, q_cut, read_align_length, read_num)?;
    writeln!(stderr, "Total reads used: {}", profile.count)?;

    // Upstream's unconditional `print('\n')` twice: the literal "\n" plus
    // print's own trailing newline is two bytes each. These go to STDOUT
    // (no `file=` argument), which is why the driver routes them to the
    // per-command stdout stream file and they still compare byte-wise.
    writeln!(stdout)?;
    writeln!(stdout)?;

    let r_path = format!("{out_prefix}.deletion_profile.r");
    File::create(&r_path)?.write_all(render_deletion_r_script(&profile, out_prefix).as_bytes())?;
    // Upstream writes the table first, then the R script; both paths are
    // created here rather than in the caller so a worker writes exactly
    // what the standalone binary writes, in the same order.
    File::create(format!("{out_prefix}.deletion_profile.txt"))?
        .write_all(render_deletion_table(&profile).as_bytes())?;

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
            record::{Flags, MappingQuality, cigar::Op},
            record_buf::{Cigar, RecordBuf, Sequence},
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
    fn forward_strand_deletion_position_and_length_filters() {
        let header = test_header();

        // Sequence length 10, CIGAR 5M2D5M: matched_portion = 5+5 = 10 ==
        // read_length; qualifies. Deletion recorded at read-position 5.
        let with_deletion = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 5), Op::new(Kind::Deletion, 2), Op::new(Kind::Match, 5)]))
            .set_sequence(Sequence::from(vec![b'A'; 10]))
            .build();

        // No deletion at all -- skipped entirely.
        let no_deletion = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .set_sequence(Sequence::from(vec![b'A'; 10]))
            .build();

        // Has a deletion but wrong sequence length (5, not 10) -- skipped.
        let wrong_length = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 3), Op::new(Kind::Deletion, 1), Op::new(Kind::Match, 2)]))
            .set_sequence(Sequence::from(vec![b'A'; 5]))
            .build();

        let records = vec![with_deletion, no_deletion, wrong_length];
        let bam_records = to_bam_records(&header, &records);
        let profile = compute_deletion_profile(bam_records.into_iter().map(Ok), 30, 10, 1_000_000).unwrap();

        assert_eq!(profile.count, 1);
        assert_eq!(profile.del_counts.len(), 10);
        assert_eq!(profile.del_counts[5], 1);
        assert_eq!(profile.del_counts.iter().sum::<u64>(), 1);
    }

    #[test]
    fn reverse_strand_flips_deletion_position() {
        let header = test_header();

        // Sequence length 10, CIGAR 3M2D7M: deletion at read-position 3
        // (forward). Reverse strand -> flipped to read_length - 3 = 7.
        let rev = RecordBuf::builder()
            .set_flags(Flags::REVERSE_COMPLEMENTED)
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 3), Op::new(Kind::Deletion, 2), Op::new(Kind::Match, 7)]))
            .set_sequence(Sequence::from(vec![b'A'; 10]))
            .build();

        let bam_records = to_bam_records(&header, &[rev]);
        let profile = compute_deletion_profile(bam_records.into_iter().map(Ok), 30, 10, 1_000_000).unwrap();

        assert_eq!(profile.count, 1);
        assert_eq!(profile.del_counts[7], 1);
        assert_eq!(profile.del_counts.iter().sum::<u64>(), 1);
    }

    #[test]
    fn read_num_caps_qualifying_reads_not_total_records() {
        let header = test_header();

        let mk = || {
            RecordBuf::builder()
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_mapping_quality(MappingQuality::new(40).unwrap())
                .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 5), Op::new(Kind::Deletion, 1), Op::new(Kind::Match, 5)]))
                .set_sequence(Sequence::from(vec![b'A'; 10]))
                .build()
        };
        let records: Vec<_> = (0..5).map(|_| mk()).collect();
        let bam_records = to_bam_records(&header, &records);
        let profile = compute_deletion_profile(bam_records.into_iter().map(Ok), 30, 10, 3).unwrap();

        assert_eq!(profile.count, 3);
        assert_eq!(profile.del_counts[5], 3);
    }

    #[test]
    fn render_table_and_r_script_exact_text() {
        let profile = DeletionProfile {
            count: 7,
            del_counts: vec![0, 2, 0, 1],
        };
        assert_eq!(
            render_deletion_table(&profile),
            "read_position\tdeletion_count\n0\t0\n1\t2\n2\t0\n3\t1\n"
        );

        // Independently derived by running the equivalent Python
        // print()/%-format expressions from oracle/upstream-src/src/
        // qcmodule/SAM.py lines 4630-4634 via `python3 -c`.
        let expected = "pdf(\"test_output.deletion_profile.pdf\")\n\
pos=c(0,1,2,3)\n\
value=c(0,2,0,1)\n\
plot(pos,value,type='b', col='blue',xlab=\"Read position (5'->3')\", ylab='Deletion count')\n\
dev.off()\n";
        assert_eq!(render_deletion_r_script(&profile, "test_output"), expected);
    }
}
