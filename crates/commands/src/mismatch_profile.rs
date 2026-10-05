//! Port of `mismatch_profile.py`: calculate the distribution of mismatches
//! across aligned reads. Contract: see `compatibility/commands.yaml` entry
//! `mismatch_profile.py`; algorithm ported from `mismatchProfile()` in
//! `oracle/upstream-src/src/qcmodule/SAM.py` (lines 4415-4552).
//!
//! Requires the BAM "MD" optional tag (parsed with the same
//! `(\d+)([A-Z]+)` pattern upstream uses via `re.findall`) and the "NM"
//! tag (reads with `NM == 0` are skipped). Like `deletion_profile.py`,
//! `read_length` is a required fixed CLI parameter, and `read_num` caps
//! the count of *qualifying* reads, checked before pulling the next
//! record.
//!
//! Preserved upstream quirk: the "Total reads used: N" summary line is
//! written to the `.xls` DATA FILE itself (as its first line, before the
//! real header) if the input is exhausted naturally, but to stderr if the
//! loop instead stops early because `count` reached `read_num` -- almost
//! certainly a copy-paste slip in upstream (the early-stop path correctly
//! targets stderr), reproduced here via `loop_exhausted_naturally` on the
//! result rather than silently unified to one destination.
//!
//! If no mismatches are found across any processed read, upstream exits
//! without writing either output file (`sys.exit()` after printing "No
//! mismatches found" to stderr) -- the CLI layer checks
//! `MismatchProfile::data.is_empty()` for this, not this module.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io;
use std::io::Write as _;
use std::process::{Command, Stdio};

use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::{Tag, Value};
use regex::Regex;

pub const ALL_GENOTYPES: [&str; 12] = [
    "A2C", "A2G", "A2T", "C2A", "C2G", "C2T", "G2A", "G2C", "G2T", "T2A", "T2C", "T2G",
];

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct MismatchProfile {
    pub count: u64,
    /// read_position -> genotype ("REF2READ") -> count, sorted by position.
    pub data: BTreeMap<usize, HashMap<String, u64>>,
    pub loop_exhausted_naturally: bool,
}

pub fn compute_mismatch_profile<I>(
    records: I,
    q_cut: u8,
    read_length: usize,
    read_num: u64,
) -> io::Result<MismatchProfile>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let md_pattern = Regex::new(r"(\d+)([A-Z]+)").expect("static regex");

    let mut count = 0u64;
    let mut data: BTreeMap<usize, HashMap<String, u64>> = BTreeMap::new();
    let mut loop_exhausted_naturally = true;

    for result in records {
        if count >= read_num {
            loop_exhausted_naturally = false;
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

        let fields = record.data();

        let nm_is_zero = matches!(
            fields.get(&Tag::EDIT_DISTANCE),
            Some(Ok(v)) if v.as_int() == Some(0)
        );

        let md_string: Option<String> = match fields.get(&Tag::MISMATCHED_POSITIONS) {
            Some(Ok(Value::String(s))) => Some(s.to_string()),
            _ => None,
        };
        let has_deletion_marker = md_string.as_deref().is_some_and(|s| s.contains('^'));

        if nm_is_zero || has_deletion_marker {
            continue;
        }

        let seq_bytes: Vec<u8> = record.sequence().iter().collect();
        if seq_bytes.len() != read_length {
            continue;
        }
        if seq_bytes.contains(&b'N') {
            continue;
        }

        let ops: Vec<_> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;
        let matched_portion: usize = ops
            .iter()
            .filter(|op| op.kind() == Kind::Match)
            .map(|op| op.len())
            .sum();
        if matched_portion != read_length {
            continue;
        }

        count += 1;

        let Some(md) = md_string else { continue };
        let reverse = flags.is_reverse_complemented();
        let mut read_coord = 0usize;

        for caps in md_pattern.captures_iter(&md) {
            let match_number: usize = caps[1].parse().expect("digits");
            let ref_base = &caps[2];

            read_coord += match_number;
            if read_coord >= seq_bytes.len() {
                // Malformed/inconsistent MD relative to the sequence
                // length; upstream would raise an uncaught IndexError and
                // crash. Rather than panic, stop processing this read's
                // remaining MD matches -- already-recorded mismatches for
                // earlier matches in this read stay recorded.
                break;
            }
            let read_base = seq_bytes[read_coord];

            if ref_base.as_bytes() == [read_base] {
                read_coord += 1;
                continue;
            }

            let genotype = format!("{ref_base}2{}", read_base as char);
            let coord = if reverse {
                read_length - read_coord - 1
            } else {
                read_coord
            };
            *data.entry(coord).or_default().entry(genotype).or_insert(0) += 1;
            read_coord += 1;
        }
    }

    Ok(MismatchProfile {
        count,
        data,
        loop_exhausted_naturally,
    })
}

/// Renders the `.mismatch_profile.xls` table. Callers should check
/// `profile.data.is_empty()` first and skip writing entirely in that case
/// (matching upstream's `sys.exit()` before either file is opened for
/// real content).
pub fn render_mismatch_table(p: &MismatchProfile) -> String {
    let mut lines = Vec::new();
    if p.loop_exhausted_naturally {
        lines.push(format!("Total reads used: {}", p.count));
    }
    lines.push(format!("read_pos\tsum\t{}", ALL_GENOTYPES.join("\t")));
    for (&pos, genotypes) in &p.data {
        let sum: u64 = genotypes.values().sum();
        let mut row = vec![pos.to_string(), sum.to_string()];
        for gt in ALL_GENOTYPES {
            row.push(genotypes.get(gt).copied().unwrap_or(0).to_string());
        }
        lines.push(row.join("\t"));
    }
    // Trailing newline: upstream's plain `print(...)` calls each add
    // their own trailing newline.
    format!("{}\n", lines.join("\n"))
}

pub fn render_mismatch_r_script(p: &MismatchProfile, out_prefix: &str) -> String {
    let positions: Vec<usize> = p.data.keys().copied().collect();

    let mut lines = Vec::new();

    for gt in ALL_GENOTYPES {
        let counts: Vec<String> = positions
            .iter()
            .map(|pos| p.data[pos].get(gt).copied().unwrap_or(0).to_string())
            .collect();
        lines.push(format!("{gt}=c({})", counts.join(",")));
    }

    lines.push(
        "color_code = c(\"green\",\"powderblue\",\"lightseagreen\",\"red\",\"violetred4\",\"mediumorchid1\",\"blue\",\"royalblue\",\"steelblue1\",\"orange\",\"gold\",\"black\")"
            .to_string(),
    );

    let log_terms: Vec<String> = ALL_GENOTYPES.iter().map(|gt| format!("log10({gt}+1)")).collect();
    lines.push(format!("y_up_bound = max(c({}))", log_terms.join(",")));
    lines.push(format!("y_low_bound = min(c({}))", log_terms.join(",")));

    lines.push(format!("pdf(\"{out_prefix}.mismatch_profile.pdf\")"));

    for (i, gt) in ALL_GENOTYPES.iter().enumerate() {
        let n = i + 1;
        if n == 1 {
            lines.push(format!(
                "plot(log10({gt}+1),type=\"l\",col=color_code[{n}],ylim=c(y_low_bound,y_up_bound),ylab=\"log10(# of mismatch)\",xlab=\"Read position (5'->3')\")"
            ));
        } else {
            lines.push(format!("lines(log10({gt}+1), col=color_code[{n}])"));
        }
    }

    let legend_genotypes: Vec<String> = ALL_GENOTYPES.iter().map(|gt| format!("\"{gt}\"")).collect();
    lines.push(format!(
        "legend(13,y_up_bound,legend=c({}), fill=color_code, border=color_code, ncol=4)",
        legend_genotypes.join(",")
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

/// Runs the `mismatch_profile.py` CLI body over an already-opened record
/// stream (C1 multi-driver pattern).
///
/// This is exactly what the standalone binary's `run()` does -- same
/// progress lines, the same two stdout blank lines, the same
/// no-mismatches early return (including which files exist and what is
/// inside them in that branch), same files with the same bytes, same
/// Rscript contract, same error propagation -- except the record source
/// is a caller-supplied iterator and stdout/stderr are caller-supplied
/// sinks. The standalone binary delegates to this (passing the process
/// streams); `rseqc_multi` passes one record broadcast plus per-command
/// stream files. The compute function and both renderers are untouched.
///
/// The empty-data branch is reproduced exactly rather than "cleaned up":
/// upstream opens both output files unconditionally before checking
/// `len(data) == 0` and then `sys.exit()`s, so both files exist; only on
/// natural iterator exhaustion has the "Total reads used" line already
/// been written to the xls. That is an output contract the differential
/// suite pins, so it is preserved byte for byte.
///
/// The argument list is long on purpose and carries a targeted lint
/// allowance: this is the C1 multi-driver pattern (one callable per
/// command carrying its full CLI surface plus the two sinks), and
/// bundling the flags into a struct would only hide them from the C2
/// cards that copy this signature command by command.
#[allow(clippy::too_many_arguments)]
pub fn run_mismatch_profile<I>(
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
    // -- the literal's own trailing space plus `end=' '` gives two spaces
    // before whatever prints next on the same stderr line.
    write!(stderr, "Process BAM file ...  ")?;
    let profile = compute_mismatch_profile(records, q_cut, read_align_length, read_num)?;

    if !profile.loop_exhausted_naturally {
        writeln!(stderr, "Total reads used: {}", profile.count)?;
    }
    // Upstream's unconditional `print('\n')` twice: the literal "\n" plus
    // print's own trailing newline is two bytes each, always, regardless
    // of whether any mismatches were found. These go to STDOUT (no
    // `file=` argument), so a worker routes them to its stdout stream
    // file and they still compare byte-wise.
    writeln!(stdout)?;
    writeln!(stdout)?;

    if profile.data.is_empty() {
        // Upstream opens both output files unconditionally before this
        // check, so they exist even though `sys.exit()` fires before the
        // table header or any data rows are written. On natural iterator
        // exhaustion the "Total reads used" line was already written
        // before this check, so that single line is present in an
        // otherwise-empty xls.
        if profile.loop_exhausted_naturally {
            File::create(format!("{out_prefix}.mismatch_profile.xls"))?
                .write_all(format!("Total reads used: {}\n", profile.count).as_bytes())?;
        } else {
            File::create(format!("{out_prefix}.mismatch_profile.xls"))?;
        }
        File::create(format!("{out_prefix}.mismatch_profile.r"))?;
        writeln!(stderr, "No mismatches found")?;
        return Ok(());
    }

    File::create(format!("{out_prefix}.mismatch_profile.xls"))?
        .write_all(render_mismatch_table(&profile).as_bytes())?;

    let r_path = format!("{out_prefix}.mismatch_profile.r");
    File::create(&r_path)?.write_all(render_mismatch_r_script(&profile, out_prefix).as_bytes())?;

    if !skip_plot {
        let rscript_path = crate::exec_resolve::which(rscript).ok_or_else(|| {
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
            record_buf::{Cigar, RecordBuf, Sequence, data::Data, data::field::Value as BufValue},
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

    fn data_with_md_nm(md: &str, nm: i32) -> Data {
        let mut data = Data::default();
        data.insert(Tag::MISMATCHED_POSITIONS, BufValue::from(md));
        data.insert(Tag::EDIT_DISTANCE, BufValue::from(nm));
        data
    }

    #[test]
    fn parses_md_tag_and_records_mismatches_forward_strand() {
        let header = test_header();

        // 10M, seq "AAAAGAAAAA", MD "4A5" means: 4 matches, ref base 'A'
        // mismatched (read has 'G' at position 4), then 5 more matches.
        let record = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .set_sequence(Sequence::from(b"AAAAGAAAAA".to_vec()))
            .set_data(data_with_md_nm("4A5", 1))
            .build();

        let bam_records = to_bam_records(&header, &[record]);
        let profile = compute_mismatch_profile(bam_records.into_iter().map(Ok), 30, 10, 1_000_000).unwrap();

        assert_eq!(profile.count, 1);
        assert_eq!(profile.data.len(), 1);
        assert_eq!(profile.data[&4].get("A2G"), Some(&1));
    }

    #[test]
    fn reverse_strand_flips_position() {
        let header = test_header();

        // Same as above but reverse-complemented: position 4 (forward) ->
        // read_length - 4 - 1 = 5.
        let record = RecordBuf::builder()
            .set_flags(Flags::REVERSE_COMPLEMENTED)
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .set_sequence(Sequence::from(b"AAAAGAAAAA".to_vec()))
            .set_data(data_with_md_nm("4A5", 1))
            .build();

        let bam_records = to_bam_records(&header, &[record]);
        let profile = compute_mismatch_profile(bam_records.into_iter().map(Ok), 30, 10, 1_000_000).unwrap();

        assert_eq!(profile.data[&5].get("A2G"), Some(&1));
    }

    #[test]
    fn skips_nm_zero_and_deletion_marker_and_n_base() {
        let header = test_header();

        let nm_zero = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .set_sequence(Sequence::from(b"AAAAAAAAAA".to_vec()))
            .set_data(data_with_md_nm("10", 0))
            .build();

        let has_deletion = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 5), Op::new(Kind::Deletion, 2), Op::new(Kind::Match, 5)]))
            .set_sequence(Sequence::from(b"AAAAAAAAAA".to_vec()))
            .set_data(data_with_md_nm("5^AA5", 2))
            .build();

        let has_n = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .set_sequence(Sequence::from(b"AAAANAAAAA".to_vec()))
            .set_data(data_with_md_nm("4A5", 1))
            .build();

        let records = vec![nm_zero, has_deletion, has_n];
        let bam_records = to_bam_records(&header, &records);
        let profile = compute_mismatch_profile(bam_records.into_iter().map(Ok), 30, 10, 1_000_000).unwrap();

        assert_eq!(profile.count, 0);
        assert!(profile.data.is_empty());
    }

    #[test]
    fn render_table_includes_total_line_only_when_exhausted_naturally() {
        let mut data = BTreeMap::new();
        let mut pos4 = HashMap::new();
        pos4.insert("A2G".to_string(), 2u64);
        data.insert(4usize, pos4);

        let exhausted = MismatchProfile { count: 5, data: data.clone(), loop_exhausted_naturally: true };
        let capped = MismatchProfile { count: 5, data, loop_exhausted_naturally: false };

        assert!(render_mismatch_table(&exhausted).starts_with("Total reads used: 5\n"));
        assert!(!render_mismatch_table(&capped).contains("Total reads used"));
    }

    #[test]
    fn render_table_exact_text() {
        let mut data = BTreeMap::new();
        let mut pos4 = HashMap::new();
        pos4.insert("A2G".to_string(), 2u64);
        data.insert(4usize, pos4);

        let profile = MismatchProfile { count: 3, data, loop_exhausted_naturally: false };
        let output = render_mismatch_table(&profile);
        let expected = "read_pos\tsum\tA2C\tA2G\tA2T\tC2A\tC2G\tC2T\tG2A\tG2C\tG2T\tT2A\tT2C\tT2G\n4\t2\t0\t2\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
        assert_eq!(output, expected);
    }

    #[test]
    fn render_r_script_exact_text() {
        // Independently derived by running the equivalent Python
        // print()/%-format expressions from oracle/upstream-src/src/
        // qcmodule/SAM.py lines 4526-4551 via `python3 -c`.
        let mut data = BTreeMap::new();
        let mut pos4 = HashMap::new();
        pos4.insert("A2G".to_string(), 2u64);
        data.insert(4usize, pos4);

        let profile = MismatchProfile { count: 3, data, loop_exhausted_naturally: false };
        let output = render_mismatch_r_script(&profile, "test_output");
        let expected = "A2C=c(0)\n\
A2G=c(2)\n\
A2T=c(0)\n\
C2A=c(0)\n\
C2G=c(0)\n\
C2T=c(0)\n\
G2A=c(0)\n\
G2C=c(0)\n\
G2T=c(0)\n\
T2A=c(0)\n\
T2C=c(0)\n\
T2G=c(0)\n\
color_code = c(\"green\",\"powderblue\",\"lightseagreen\",\"red\",\"violetred4\",\"mediumorchid1\",\"blue\",\"royalblue\",\"steelblue1\",\"orange\",\"gold\",\"black\")\n\
y_up_bound = max(c(log10(A2C+1),log10(A2G+1),log10(A2T+1),log10(C2A+1),log10(C2G+1),log10(C2T+1),log10(G2A+1),log10(G2C+1),log10(G2T+1),log10(T2A+1),log10(T2C+1),log10(T2G+1)))\n\
y_low_bound = min(c(log10(A2C+1),log10(A2G+1),log10(A2T+1),log10(C2A+1),log10(C2G+1),log10(C2T+1),log10(G2A+1),log10(G2C+1),log10(G2T+1),log10(T2A+1),log10(T2C+1),log10(T2G+1)))\n\
pdf(\"test_output.mismatch_profile.pdf\")\n\
plot(log10(A2C+1),type=\"l\",col=color_code[1],ylim=c(y_low_bound,y_up_bound),ylab=\"log10(# of mismatch)\",xlab=\"Read position (5'->3')\")\n\
lines(log10(A2G+1), col=color_code[2])\n\
lines(log10(A2T+1), col=color_code[3])\n\
lines(log10(C2A+1), col=color_code[4])\n\
lines(log10(C2G+1), col=color_code[5])\n\
lines(log10(C2T+1), col=color_code[6])\n\
lines(log10(G2A+1), col=color_code[7])\n\
lines(log10(G2C+1), col=color_code[8])\n\
lines(log10(G2T+1), col=color_code[9])\n\
lines(log10(T2A+1), col=color_code[10])\n\
lines(log10(T2C+1), col=color_code[11])\n\
lines(log10(T2G+1), col=color_code[12])\n\
legend(13,y_up_bound,legend=c(\"A2C\",\"A2G\",\"A2T\",\"C2A\",\"C2G\",\"C2T\",\"G2A\",\"G2C\",\"G2T\",\"T2A\",\"T2C\",\"T2G\"), fill=color_code, border=color_code, ncol=4)\n\
dev.off()\n";
        assert_eq!(output, expected);
    }
}
