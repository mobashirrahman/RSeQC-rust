//! Port of `clipping_profile.py`: estimate the soft-clipping profile of
//! reads. Contract: see `compatibility/commands.yaml` entry
//! `clipping_profile.py`; algorithm ported from `clipping_profile()` in
//! `oracle/upstream-src/src/qcmodule/SAM.py` (lines 3239-3360).
//!
//! Shared with a future `insertion_profile.py` port: the compute functions
//! here take a `clip_char` parameter (upstream's `type="S"` vs `type="I"`)
//! so the same aggregation logic serves both; only the CLI/render layer
//! differs (file names, axis labels).
//!
//! Notable upstream behavior preserved exactly: mapq/unmapped/QC-fail
//! filtering; `total_read`(s) count every filter-passing record, including
//! ones with no clip operation at all; the output row count follows the
//! LAST filter-passing record's expanded-CIGAR length, same quirk pattern
//! as DIV-0008/0009 (see DIV-0010). Upstream's per-position counts are
//! Python floats (`+= 1.0`) and print with a trailing `.0`; reproduced here
//! via an explicit formatter rather than using floats internally, since the
//! values are always whole numbers.

use std::collections::HashMap;
use std::io;

use noodles_bam as bam;
use rseqc_formats::cigar::expand_cigar_to_read_ops;

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct SingleEndProfile {
    pub total_read: u64,
    pub clip_count: Vec<u64>,
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct PairedEndProfile {
    pub total_read1: u64,
    pub total_read2: u64,
    pub r1_clip_count: Vec<u64>,
    pub r2_clip_count: Vec<u64>,
}

fn passes_common_filters(record: &bam::Record, q_cut: u8) -> bool {
    let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
    if mapq < q_cut {
        return false;
    }
    if record.flags().is_unmapped() {
        return false;
    }
    if record.flags().is_qc_fail() {
        return false;
    }
    true
}

fn cigar_ops(record: &bam::Record) -> io::Result<Vec<noodles_sam::alignment::record::cigar::Op>> {
    record.cigar().iter().collect()
}

/// Single-end clipping profile: matches the `PE is False` branch of
/// upstream's `clipping_profile()`.
pub fn compute_single_end<I>(records: I, q_cut: u8, clip_char: u8) -> io::Result<SingleEndProfile>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut total_read = 0u64;
    let mut profile: HashMap<usize, u64> = HashMap::new();
    let mut last_len = 0usize;

    for result in records {
        let record = result?;
        if !passes_common_filters(&record, q_cut) {
            continue;
        }

        total_read += 1;
        let mut expanded = expand_cigar_to_read_ops(cigar_ops(&record)?);
        last_len = expanded.len();

        if !expanded.contains(&clip_char) {
            continue;
        }
        if record.flags().is_reverse_complemented() {
            expanded.reverse();
        }
        for (i, &b) in expanded.iter().enumerate() {
            if b == clip_char {
                *profile.entry(i).or_insert(0) += 1;
            }
        }
    }

    let clip_count = (0..last_len).map(|i| *profile.get(&i).unwrap_or(&0)).collect();
    Ok(SingleEndProfile { total_read, clip_count })
}

/// Paired-end clipping profile: matches the `PE is True` branch. Records
/// that aren't segmented (`is_paired` in pysam / SAM flag 0x1) are skipped
/// entirely, same as upstream.
pub fn compute_paired_end<I>(records: I, q_cut: u8, clip_char: u8) -> io::Result<PairedEndProfile>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut total_read1 = 0u64;
    let mut total_read2 = 0u64;
    let mut r1_profile: HashMap<usize, u64> = HashMap::new();
    let mut r2_profile: HashMap<usize, u64> = HashMap::new();
    let mut last_len = 0usize;

    for result in records {
        let record = result?;
        if !passes_common_filters(&record, q_cut) {
            continue;
        }

        let flags = record.flags();
        if !flags.is_segmented() {
            continue;
        }
        if flags.is_first_segment() {
            total_read1 += 1;
        }
        if flags.is_last_segment() {
            total_read2 += 1;
        }

        let mut expanded = expand_cigar_to_read_ops(cigar_ops(&record)?);
        if flags.is_reverse_complemented() {
            expanded.reverse();
        }
        last_len = expanded.len();

        if !expanded.contains(&clip_char) {
            continue;
        }

        if flags.is_first_segment() {
            for (i, &b) in expanded.iter().enumerate() {
                if b == clip_char {
                    *r1_profile.entry(i).or_insert(0) += 1;
                }
            }
        }
        if flags.is_last_segment() {
            for (i, &b) in expanded.iter().enumerate() {
                if b == clip_char {
                    *r2_profile.entry(i).or_insert(0) += 1;
                }
            }
        }
    }

    let r1_clip_count = (0..last_len).map(|i| *r1_profile.get(&i).unwrap_or(&0)).collect();
    let r2_clip_count = (0..last_len).map(|i| *r2_profile.get(&i).unwrap_or(&0)).collect();
    Ok(PairedEndProfile {
        total_read1,
        total_read2,
        r1_clip_count,
        r2_clip_count,
    })
}

/// Matches Python's `str(float)` for the whole-number floats this module
/// accumulates (upstream increments counts by `1.0`, never fractionally).
fn fmt_float(n: u64) -> String {
    format!("{n}.0")
}

/// Every line, including the last, ends with `\n` (upstream's plain
/// `print(...)` calls each add their own trailing newline) -- applies
/// to all four render functions in this module.
pub fn render_single_table(p: &SingleEndProfile) -> String {
    let mut out = String::from("Position\tClipped_nt\tNon_clipped_nt\n");
    for (i, &c) in p.clip_count.iter().enumerate() {
        let non_clip = p.total_read - c;
        out.push_str(&format!("{i}\t{}\t{}\n", fmt_float(c), fmt_float(non_clip)));
    }
    out
}

pub fn render_single_r_script(p: &SingleEndProfile, out_prefix: &str) -> String {
    let read_pos: Vec<String> = (0..p.clip_count.len()).map(|i| i.to_string()).collect();
    let clip_strs: Vec<String> = p.clip_count.iter().map(|&c| fmt_float(c)).collect();

    format!(
        "pdf(\"{out_prefix}.clipping_profile.pdf\")\nread_pos=c({})\nclip_count=c({})\nnonclip_count= {} - clip_count\nplot(read_pos, nonclip_count*100/(clip_count+nonclip_count),col=\"blue\",main=\"clipping profile\",xlab=\"Position of read\",ylab=\"Non-clipped %\",type=\"b\")\ndev.off()\n",
        read_pos.join(","),
        clip_strs.join(","),
        p.total_read,
    )
}

pub fn render_paired_table(p: &PairedEndProfile) -> String {
    let mut out = String::from("Position\tClipped_nt\tNon_clipped_nt\nRead-1:\n");
    for (i, &c) in p.r1_clip_count.iter().enumerate() {
        let non_clip = p.total_read1 - c;
        out.push_str(&format!("{i}\t{}\t{}\n", fmt_float(c), fmt_float(non_clip)));
    }
    out.push_str("Read-2:\n");
    for (i, &c) in p.r2_clip_count.iter().enumerate() {
        let non_clip = p.total_read2 - c;
        out.push_str(&format!("{i}\t{}\t{}\n", fmt_float(c), fmt_float(non_clip)));
    }
    out
}

pub fn render_paired_r_script(p: &PairedEndProfile, out_prefix: &str) -> String {
    // Upstream computes `read_pos` once (from the shared last cigar_str)
    // and reuses it for both the R1 and R2 sections; r1/r2_clip_count have
    // the same length by construction (both built over 0..last_len).
    let read_pos: Vec<String> = (0..p.r1_clip_count.len()).map(|i| i.to_string()).collect();
    let r1_strs: Vec<String> = p.r1_clip_count.iter().map(|&c| fmt_float(c)).collect();
    let r2_strs: Vec<String> = p.r2_clip_count.iter().map(|&c| fmt_float(c)).collect();
    let read_pos_csv = read_pos.join(",");

    format!(
        "pdf(\"{out_prefix}.clipping_profile.R1.pdf\")\nread_pos=c({read_pos_csv})\nr1_clip_count=c({})\nr1_nonclip_count = {} - r1_clip_count\nplot(read_pos, r1_nonclip_count*100/(r1_clip_count + r1_nonclip_count),col=\"blue\",main=\"clipping profile\",xlab=\"Position of read (read-1)\",ylab=\"Non-clipped %\",type=\"b\")\ndev.off()\n\
pdf(\"{out_prefix}.clipping_profile.R2.pdf\")\nread_pos=c({read_pos_csv})\nr2_clip_count=c({})\nr2_nonclip_count = {} - r2_clip_count\nplot(read_pos, r2_nonclip_count*100/(r2_clip_count + r2_nonclip_count),col=\"blue\",main=\"clipping profile\",xlab=\"Position of read (read-2)\",ylab=\"Non-clipped %\",type=\"b\")\ndev.off()\n",
        r1_strs.join(","),
        p.total_read1,
        r2_strs.join(","),
        p.total_read2,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::{
        self as sam,
        alignment::{
            io::Write as _,
            record::{Flags, MappingQuality, cigar::Op, cigar::op::Kind},
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
    fn single_end_counts_and_last_record_length_quirk() {
        let header = test_header();

        // 9M1S (length 10, one soft-clip at the end -> index 9).
        let clipped = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 9), Op::new(Kind::SoftClip, 1)]))
            .build();

        // 4M (length 4, no clipping at all) -- this is the LAST record, so
        // the output table/vectors should only have 4 rows, not 10.
        let unclipped_short = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 4)]))
            .build();

        let records = vec![clipped, unclipped_short];
        let bam_records = to_bam_records(&header, &records);
        let profile = compute_single_end(bam_records.into_iter().map(Ok), 30, b'S').unwrap();

        assert_eq!(profile.total_read, 2);
        assert_eq!(profile.clip_count.len(), 4, "last record's length (4), not the clipped record's (10)");
        assert_eq!(profile.clip_count, vec![0, 0, 0, 0]);
    }

    #[test]
    fn single_end_reverse_strand_reverses_clip_position() {
        let header = test_header();

        // 1S9M reversed (is_reverse) becomes "MMMMMMMMMS" -> clip at index 9.
        let rev_clipped = RecordBuf::builder()
            .set_flags(Flags::REVERSE_COMPLEMENTED)
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::SoftClip, 1), Op::new(Kind::Match, 9)]))
            .build();

        let records = vec![rev_clipped];
        let bam_records = to_bam_records(&header, &records);
        let profile = compute_single_end(bam_records.into_iter().map(Ok), 30, b'S').unwrap();

        assert_eq!(profile.total_read, 1);
        assert_eq!(profile.clip_count.len(), 10);
        assert_eq!(profile.clip_count[9], 1);
        assert_eq!(profile.clip_count[0..9], vec![0; 9]);
    }

    #[test]
    fn render_single_table_uses_python_float_style() {
        let profile = SingleEndProfile {
            total_read: 5,
            clip_count: vec![2, 0],
        };
        let output = render_single_table(&profile);
        assert_eq!(
            output,
            "Position\tClipped_nt\tNon_clipped_nt\n0\t2.0\t3.0\n1\t0.0\t5.0\n"
        );
    }

    #[test]
    fn render_single_r_script_exact_text() {
        let profile = SingleEndProfile {
            total_read: 5,
            clip_count: vec![2, 0],
        };
        let output = render_single_r_script(&profile, "test_output");
        let expected = "pdf(\"test_output.clipping_profile.pdf\")\n\
read_pos=c(0,1)\n\
clip_count=c(2.0,0.0)\n\
nonclip_count= 5 - clip_count\n\
plot(read_pos, nonclip_count*100/(clip_count+nonclip_count),col=\"blue\",main=\"clipping profile\",xlab=\"Position of read\",ylab=\"Non-clipped %\",type=\"b\")\n\
dev.off()\n";
        assert_eq!(output, expected);
    }

    #[test]
    fn paired_end_routes_by_segment_and_shares_read_pos_length() {
        let header = test_header();

        // Read1, 3M1S (clip at index 3).
        let r1 = RecordBuf::builder()
            .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 3), Op::new(Kind::SoftClip, 1)]))
            .build();

        // Read2, 4M (no clip) -- last record, sets read_pos length to 4.
        let r2 = RecordBuf::builder()
            .set_flags(Flags::SEGMENTED | Flags::LAST_SEGMENT)
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 4)]))
            .build();

        let records = vec![r1, r2];
        let bam_records = to_bam_records(&header, &records);
        let profile = compute_paired_end(bam_records.into_iter().map(Ok), 30, b'S').unwrap();

        assert_eq!(profile.total_read1, 1);
        assert_eq!(profile.total_read2, 1);
        assert_eq!(profile.r1_clip_count.len(), 4);
        assert_eq!(profile.r2_clip_count.len(), 4);
        assert_eq!(profile.r1_clip_count, vec![0, 0, 0, 1]);
        assert_eq!(profile.r2_clip_count, vec![0, 0, 0, 0]);
    }

    #[test]
    fn render_paired_r_script_exact_text() {
        // Expected text independently derived by running the equivalent
        // Python print()/%-format expressions from
        // oracle/upstream-src/src/qcmodule/SAM.py lines 3348-3360 via
        // `python3 -c`, not by copying this function's own output back in.
        let profile = PairedEndProfile {
            total_read1: 3,
            total_read2: 2,
            r1_clip_count: vec![1, 0],
            r2_clip_count: vec![0, 2],
        };
        let output = render_paired_r_script(&profile, "test_output");
        let expected = "pdf(\"test_output.clipping_profile.R1.pdf\")\n\
read_pos=c(0,1)\n\
r1_clip_count=c(1.0,0.0)\n\
r1_nonclip_count = 3 - r1_clip_count\n\
plot(read_pos, r1_nonclip_count*100/(r1_clip_count + r1_nonclip_count),col=\"blue\",main=\"clipping profile\",xlab=\"Position of read (read-1)\",ylab=\"Non-clipped %\",type=\"b\")\n\
dev.off()\n\
pdf(\"test_output.clipping_profile.R2.pdf\")\n\
read_pos=c(0,1)\n\
r2_clip_count=c(0.0,2.0)\n\
r2_nonclip_count = 2 - r2_clip_count\n\
plot(read_pos, r2_nonclip_count*100/(r2_clip_count + r2_nonclip_count),col=\"blue\",main=\"clipping profile\",xlab=\"Position of read (read-2)\",ylab=\"Non-clipped %\",type=\"b\")\n\
dev.off()\n";
        assert_eq!(output, expected);
    }

    #[test]
    fn unpaired_records_are_skipped_in_paired_mode() {
        let header = test_header();
        let unpaired = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 4), Op::new(Kind::SoftClip, 1)]))
            .build();
        let bam_records = to_bam_records(&header, &[unpaired]);
        let profile = compute_paired_end(bam_records.into_iter().map(Ok), 30, b'S').unwrap();
        assert_eq!(profile.total_read1, 0);
        assert_eq!(profile.total_read2, 0);
        assert_eq!(profile.r1_clip_count.len(), 0);
    }
}
