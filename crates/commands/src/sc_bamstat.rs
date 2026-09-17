//! Port of `sc_bamStat.py`: mapping-statistics summary for a Cell
//! Ranger-style tagged single-cell BAM. Ports `qcmodule.scbam.
//! mapping_stat` (`oracle/upstream-src/src/qcmodule/scbam.py`,
//! `read_match_type`/`list2str`).
//!
//! No REAL BAI index is opened or queried: upstream's per-chromosome
//! `samfile.fetch(chr_id)` loop is replicated as a SINGLE streaming
//! pass over the input, grouped by consecutive `reference_sequence_id`
//! runs assumed to already be coordinate-sorted (see `mapping_stat`'s
//! own doc comment for the full reasoning, including DIV-0019's
//! now-closed per-chromosome progress-message requirement). The
//! intermediate `<chrom>.all_reads_id.txt`/`<chrom>.confident_reads_id
//! .txt` files plus the `awk '!a[$0]++' ... | wc -l` subprocess
//! dedup-count are skipped entirely: they exist only to count DISTINCT
//! read (QNAME) values, which a plain `HashSet<String>` computes
//! directly and losslessly (same final `total_reads_n`/`confi_reads_n`,
//! no behavior change) -- see DIV-0018 for the CWD-pollution-file
//! aspect specifically (still accepted, not replicated).
//!
//! **Preserves a real, confirmed upstream bug, do not "fix"**: the
//! region-type tally is
//! ```python
//! if RE_tag in tag_dict:
//!     if tag_dict[RE_tag] == "E": exon_reads += 1
//!     elif tag_dict[RE_tag] == "I": intron_reads += 1
//! elif tag_dict[RE_tag] == "N":
//!     intergenic_reads += 1
//! else:
//!     other_reads1 += 1
//! ```
//! The `elif`/`else` are attached to the OUTER `if RE_tag in tag_dict`,
//! not the inner `E`/`I` check. So: (1) a confidently-mapped read whose
//! RE tag is present but is neither `"E"` nor `"I"` (e.g. `"N"`) hits
//! NEITHER inner branch NOR `other_reads1` -- a silent no-op, not even
//! counted as "other". (2) A confidently-mapped read whose RE tag is
//! ABSENT unconditionally evaluates `tag_dict[RE_tag]` on a missing key
//! -- an uncaught `KeyError` that aborts the whole command (caught by the
//! CLI's own `except (..., KeyError)`, so it's a clean exit-1 failure,
//! not a raw traceback, but a SINGLE RE-tag-less confident read still
//! aborts the entire analysis). Replicated exactly: present-but-other
//! value is a silent no-op; absent is a propagated `io::Error`.
use std::collections::{HashMap, HashSet};
use std::io;

use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::record::data::field::Tag};

fn tag(name: &str) -> Tag {
    let bytes = name.as_bytes();
    Tag::new(bytes[0], bytes[1])
}

#[derive(Debug, Clone)]
pub struct TagNames {
    pub cb: String,
    pub umi: String,
    pub re: String,
    pub tx: String,
    pub an: String,
    pub xf: String,
}

impl Default for TagNames {
    fn default() -> Self {
        Self { cb: "CB".to_string(), umi: "UB".to_string(), re: "RE".to_string(), tx: "TX".to_string(), an: "AN".to_string(), xf: "xf".to_string() }
    }
}

/// Ports `read_match_type`, matching directly on the parsed CIGAR
/// operation sequence rather than round-tripping through a string and
/// regex (behaviorally identical: each upstream regex requires the
/// WHOLE cigar string to consist of exactly that op sequence).
fn read_match_type(ops: &[(sam::alignment::record::cigar::op::Kind, usize)]) -> &'static str {
    use sam::alignment::record::cigar::op::Kind::*;
    match ops {
        [(Match, _)] => "Map_consecutively",
        [(Match, _), (Skip, _), (Match, _)] => "Map_with_splicing",
        [(SoftClip, _), (Match, _)] | [(Match, _), (SoftClip, _)] => "Map_with_clipping",
        [(Match, _), (Skip, _), (Match, _), (SoftClip, _)] | [(SoftClip, _), (Match, _), (Skip, _), (Match, _)] => "Map_with_splicing_and_clipping",
        _ => "Others",
    }
}

#[derive(Debug, Default, Clone)]
pub struct MappingStats {
    pub total_alignments: i64,
    pub confi_alignments: i64,
    pub total_reads_n: i64,
    pub confi_reads_n: i64,
    pub confi_reads_dup: i64,
    pub confi_reads_nondup: i64,
    pub confi_reads_rev: i64,
    pub confi_reads_fwd: i64,
    pub confi_cb: i64,
    pub confi_ub: i64,
    pub exon_reads: i64,
    pub intron_reads: i64,
    pub intergenic_reads: i64,
    pub other_reads1: i64,
    pub sense_reads: i64,
    pub anti_reads: i64,
    pub other_reads2: i64,
    pub chrm_reads: i64,
    pub read_type: HashMap<String, i64>,
}

/// Runs the full per-alignment counting pass. Ports `mapping_stat`,
/// INCLUDING its per-chromosome iteration shape (DIV-0019, now
/// closed): upstream visits every chromosome in `samfile.references`
/// order via `samfile.fetch(chr_id)` (no start/end -- the whole
/// chromosome), printing `logging.info` "Processing"/"Processed"
/// lines around each one, EVEN chromosomes with zero alignments (the
/// header-chromosome loop is unconditional). A read whose
/// `reference_sequence_id` doesn't resolve to any header chromosome at
/// all (a truly unmapped read with no RNAME) is therefore never
/// visited by any per-chromosome fetch and silently excluded from
/// EVERY count upstream computes -- replicated here by skipping such
/// records entirely, not counting them toward `total_alignments`.
/// Assumes coordinate-sorted input (same assumption `require_index`
/// already implies upstream), so a single streaming pass grouped by
/// consecutive `reference_sequence_id` matches per-chromosome fetch
/// order without needing real BAI-indexed queries (no BAI support
/// anywhere in this port, see crates/cli/src/bin/sc_bamstat.rs docs).
pub fn mapping_stat<I>(records: I, header: &sam::Header, tags: &TagNames, chrm_id: &str) -> io::Result<MappingStats>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let cb_tag = tag(&tags.cb);
    let umi_tag = tag(&tags.umi);
    let re_tag = tag(&tags.re);
    let tx_tag = tag(&tags.tx);
    let an_tag = tag(&tags.an);
    let xf_tag = tag(&tags.xf);

    let mut s = MappingStats::default();
    let mut all_reads: HashSet<String> = HashSet::new();
    let mut confi_reads: HashSet<String> = HashSet::new();

    let ref_names: Vec<String> = header.reference_sequences().keys().map(|k| k.to_string()).collect();
    let mut iter = records.into_iter().peekable();

    for (chrom_idx, chrom_name_ref) in ref_names.iter().enumerate() {
        eprintln!("Processing \"{chrom_name_ref}\" ...");
        let mut chrom_count: i64 = 0;

        loop {
            let is_current_chrom = match iter.peek() {
                Some(Ok(record)) => record.reference_sequence_id().transpose()?.map(|id| id == chrom_idx).unwrap_or(false),
                Some(Err(_)) => true, // surface the error by consuming it below
                None => false,
            };
            if !is_current_chrom {
                break;
            }

            let record = iter.next().unwrap()?;
            chrom_count += 1;
            s.total_alignments += 1;

            let read_id = record.name().map(|n| String::from_utf8_lossy(n).into_owned()).unwrap_or_default();
            all_reads.insert(read_id.clone());

            let data = record.data();

            let is_confident = data.get(&xf_tag).and_then(|r| r.ok()).and_then(|v| v.as_int()).map(|n| n & 1 != 0).unwrap_or(false);

            if is_confident {
                s.confi_alignments += 1;
                confi_reads.insert(read_id);

                if chrom_name_ref == chrm_id {
                    s.chrm_reads += 1;
                }

                if data.get(&cb_tag).is_some() {
                    s.confi_cb += 1;
                }
                if data.get(&umi_tag).is_some() {
                    s.confi_ub += 1;
                }

                let flags = record.flags();
                if flags.is_duplicate() {
                    s.confi_reads_dup += 1;
                } else {
                    s.confi_reads_nondup += 1;
                }
                if flags.is_reverse_complemented() {
                    s.confi_reads_rev += 1;
                } else {
                    s.confi_reads_fwd += 1;
                }

                match data.get(&re_tag) {
                    Some(v) => {
                        let value = v?;
                        let ch = match value {
                            sam::alignment::record::data::field::Value::Character(c) => Some(c as char),
                            sam::alignment::record::data::field::Value::String(s) => s.first().map(|&b| b as char),
                            _ => None,
                        };
                        match ch {
                            Some('E') => s.exon_reads += 1,
                            Some('I') => s.intron_reads += 1,
                            _ => {} // present but neither E nor I: silent no-op (see module docs)
                        }
                    }
                    None => {
                        return Err(io::Error::new(io::ErrorKind::NotFound, format!("'{}'", tags.re)));
                    }
                }

                if data.get(&tx_tag).is_some() {
                    s.sense_reads += 1;
                } else if data.get(&an_tag).is_some() {
                    s.anti_reads += 1;
                } else {
                    s.other_reads2 += 1;
                }

                let ops: Vec<(sam::alignment::record::cigar::op::Kind, usize)> = record.cigar().iter().map(|r| r.map(|op| (op.kind(), op.len()))).collect::<Result<_, _>>()?;
                let mt = read_match_type(&ops);
                *s.read_type.entry(mt.to_string()).or_insert(0) += 1;
            }
        }

        eprintln!("Processed {chrom_count} alignments from \"{chrom_name_ref}\"");
    }

    // Records whose reference doesn't match any (remaining) header
    // chromosome in order -- e.g. a truly unmapped read with no RNAME,
    // or one appearing out of coordinate-sorted order -- are drained
    // and skipped, matching upstream's silent per-chromosome-fetch
    // exclusion. Errors are still surfaced.
    for item in iter {
        item?;
    }

    eprintln!("Processing total {} alignments mapped to all chromosomes.", s.total_alignments);
    eprintln!("Count total mapped reads ...");
    eprintln!("Count confidently mapped reads ...");
    eprintln!("Removing intermediate files ...");

    s.total_reads_n = all_reads.len() as i64;
    s.confi_reads_n = confi_reads.len() as i64;

    Ok(s)
}

/// Renders the stdout mapping-statistics report. Ports the print-
/// statement block at the end of `mapping_stat`. Returns an error if
/// `total_reads_n` or `confi_reads_n` is zero (matching upstream's
/// uncaught `ZeroDivisionError` on the percentage calculations).
pub fn render_report(s: &MappingStats) -> io::Result<String> {
    if s.total_reads_n == 0 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "division by zero: total_reads_n is 0"));
    }
    if s.confi_reads_n == 0 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "division by zero: confi_reads_n is 0"));
    }

    let non_confi = s.total_reads_n - s.confi_reads_n;
    let t = s.total_reads_n as f64;
    let c = s.confi_reads_n as f64;
    let pct_t = |n: i64| n as f64 * 100.0 / t;
    let pct_c = |n: i64| n as f64 * 100.0 / c;

    let mut out = String::new();
    out.push('\n');
    out.push_str(&format!("\nTotal_alignments: {}\n", s.total_alignments));
    out.push_str(&format!("└--Confident_alignments: {}\n", s.confi_alignments));
    out.push('\n');
    out.push_str(&format!("Total_mapped_reads:\t{}\n", s.total_reads_n));
    out.push_str(&format!("|--Non_confidently_mapped_reads:\t{}\t({:.2}%)\n", non_confi, pct_t(non_confi)));
    out.push_str(&format!("└--Confidently_mapped_reads:\t{}\t({:.2}%)\n", s.confi_reads_n, pct_t(s.confi_reads_n)));
    out.push_str(&format!("   |--Reads_with_PCR_duplicates:\t{}\t({:.2}%)\n", s.confi_reads_dup, pct_c(s.confi_reads_dup)));
    out.push_str(&format!("   └--Reads_no_PCR_duplicates:\t{}\t({:.2}%)\n", s.confi_reads_nondup, pct_c(s.confi_reads_nondup)));
    out.push('\n');
    out.push_str(&format!("   |--Reads_map_to_forward(Waston)_strand:\t{}\t({:.2}%)\n", s.confi_reads_fwd, pct_c(s.confi_reads_fwd)));
    out.push_str(&format!("   └--Reads_map_to_Reverse(Crick)_strand:\t{}\t({:.2}%)\n", s.confi_reads_rev, pct_c(s.confi_reads_rev)));
    out.push('\n');
    out.push_str(&format!("   |--Reads_map_to_sense_strand:\t{}\t({:.2}%)\n", s.sense_reads, pct_c(s.sense_reads)));
    out.push_str(&format!("   └--Reads_map_to_antisense_strand:\t{}\t({:.2}%)\n", s.anti_reads, pct_c(s.anti_reads)));
    out.push_str(&format!("   └--Other:\t{}\t({:.2}%)\n", s.other_reads2, pct_c(s.other_reads2)));
    out.push('\n');
    out.push_str(&format!("   |--Reads_map_to_exons:\t{}\t({:.2}%)\n", s.exon_reads, pct_c(s.exon_reads)));
    out.push_str(&format!("   └--Reads_map_to_introns:\t{}\t({:.2}%)\n", s.intron_reads, pct_c(s.intron_reads)));
    out.push_str(&format!("   └--Reads_map_to_intergenic:\t{}\t({:.2}%)\n", s.intergenic_reads, pct_c(s.intergenic_reads)));
    out.push_str(&format!("   └--Other:\t{}\t({:.2}%)\n", s.other_reads1, pct_c(s.other_reads1)));
    out.push('\n');
    out.push_str(&format!("   |--Reads_with_error-corrected_barcode:\t{}\t({:.2}%)\n", s.confi_cb, pct_c(s.confi_cb)));
    out.push_str(&format!("   └--Reads_no_error-corrected_barcode:\t{}\t({:.2}%)\n", s.confi_reads_n - s.confi_cb, pct_c(s.confi_reads_n - s.confi_cb)));
    out.push('\n');
    out.push_str(&format!("   |--Reads_with_error-corrected_UMI:\t{}\t({:.2}%)\n", s.confi_ub, pct_c(s.confi_ub)));
    out.push_str(&format!("   └--Reads_no_error-corrected_UMI:\t{}\t({:.2}%)\n", s.confi_reads_n - s.confi_ub, pct_c(s.confi_reads_n - s.confi_ub)));
    out.push('\n');
    out.push_str(&format!("   |--Reads_map_to_mitochonrial_genome:\t{}\t({:.2}%)\n", s.chrm_reads, pct_c(s.chrm_reads)));
    out.push_str(&format!("   └--Reads_map_to_nuclear_genome:\t{}\t({:.2}%)\n", s.confi_reads_n - s.chrm_reads, pct_c(s.confi_reads_n - s.chrm_reads)));
    out.push('\n');

    let mut kinds: Vec<&String> = s.read_type.keys().filter(|k| k.as_str() != "Others").collect();
    kinds.sort();
    for kind in kinds {
        let n = s.read_type[kind];
        out.push_str(&format!("   |--{kind}:\t{n}\t({:.2}%)\n", pct_c(n)));
    }
    let others = s.read_type.get("Others").copied().unwrap_or(0);
    out.push_str(&format!("   └--Others:\t{others}\t({:.2}%)\n", pct_c(others)));
    out.push('\n');

    Ok(out)
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

    #[test]
    fn read_match_type_classifies_known_shapes() {
        use Kind::*;
        assert_eq!(read_match_type(&[(Match, 50)]), "Map_consecutively");
        assert_eq!(read_match_type(&[(Match, 20), (Skip, 100), (Match, 30)]), "Map_with_splicing");
        assert_eq!(read_match_type(&[(SoftClip, 5), (Match, 45)]), "Map_with_clipping");
        assert_eq!(read_match_type(&[(Match, 45), (SoftClip, 5)]), "Map_with_clipping");
        assert_eq!(read_match_type(&[(Match, 20), (Skip, 100), (Match, 20), (SoftClip, 10)]), "Map_with_splicing_and_clipping");
        assert_eq!(read_match_type(&[(Insertion, 3), (Match, 50)]), "Others");
    }

    #[allow(clippy::too_many_arguments)]
    fn confident_record(name: &str, chrom_id: usize, reverse: bool, dup: bool, cb: bool, umi: bool, re: Option<char>, tx: bool, an: bool) -> RecordBuf {
        let mut data = Data::default();
        data.insert(tag("xf"), BufValue::from(1i32));
        if cb {
            data.insert(tag("CB"), BufValue::from("AAAA"));
        }
        if umi {
            data.insert(tag("UB"), BufValue::from("TTTT"));
        }
        if let Some(c) = re {
            data.insert(tag("RE"), BufValue::Character(c as u8));
        }
        if tx {
            data.insert(tag("TX"), BufValue::from(1i32));
        }
        if an {
            data.insert(tag("AN"), BufValue::from(1i32));
        }
        let mut flags = Flags::empty();
        if reverse {
            flags |= Flags::REVERSE_COMPLEMENTED;
        }
        if dup {
            flags |= Flags::DUPLICATE;
        }
        RecordBuf::builder()
            .set_name(name)
            .set_flags(flags)
            .set_reference_sequence_id(chrom_id)
            .set_alignment_start(noodles_core::Position::try_from(1).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 10]))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![40; 10]))
            .set_data(data)
            .build()
    }

    #[test]
    fn mapping_stat_counts_confident_read_correctly() {
        let header = header_with_chrom("chr1", 1000);
        let rec = confident_record("r1", 0, false, false, true, true, Some('E'), true, false);
        let bam_records = to_bam_records(&header, &[rec]);
        let tags = TagNames::default();
        let stats = mapping_stat(bam_records.into_iter().map(Ok), &header, &tags, "chrM").unwrap();

        assert_eq!(stats.total_alignments, 1);
        assert_eq!(stats.confi_alignments, 1);
        assert_eq!(stats.total_reads_n, 1);
        assert_eq!(stats.confi_reads_n, 1);
        assert_eq!(stats.confi_cb, 1);
        assert_eq!(stats.confi_ub, 1);
        assert_eq!(stats.confi_reads_fwd, 1);
        assert_eq!(stats.confi_reads_nondup, 1);
        assert_eq!(stats.exon_reads, 1);
        assert_eq!(stats.sense_reads, 1);
        assert_eq!(stats.read_type["Map_consecutively"], 1);
    }

    #[test]
    fn mapping_stat_missing_re_tag_on_confident_read_is_a_hard_error() {
        // Replicates the upstream KeyError-on-missing-RE_tag bug.
        let header = header_with_chrom("chr1", 1000);
        let rec = confident_record("r1", 0, false, false, false, false, None, false, false);
        let bam_records = to_bam_records(&header, &[rec]);
        let tags = TagNames::default();
        let err = mapping_stat(bam_records.into_iter().map(Ok), &header, &tags, "chrM").unwrap_err();
        assert!(err.to_string().contains("RE"));
    }

    #[test]
    fn mapping_stat_re_tag_present_but_other_value_is_silent_noop() {
        // RE present but value is neither E nor I ("N") -- upstream's
        // broken if/elif/else structure means this counts as NEITHER
        // exon/intron NOR "other" -- a true no-op, not even other_reads1.
        let header = header_with_chrom("chr1", 1000);
        let rec = confident_record("r1", 0, false, false, false, false, Some('N'), false, false);
        let bam_records = to_bam_records(&header, &[rec]);
        let tags = TagNames::default();
        let stats = mapping_stat(bam_records.into_iter().map(Ok), &header, &tags, "chrM").unwrap();
        assert_eq!(stats.exon_reads, 0);
        assert_eq!(stats.intron_reads, 0);
        assert_eq!(stats.intergenic_reads, 0);
        assert_eq!(stats.other_reads1, 0); // never incremented anywhere -- confirms it stays dead code
    }

    #[test]
    fn mapping_stat_chrm_detection_and_non_confident_reads_excluded() {
        let header = sam::Header::builder()
            .add_reference_sequence("chrM", sam::header::record::value::Map::<sam::header::record::value::map::ReferenceSequence>::new(std::num::NonZeroUsize::new(1000).unwrap()))
            .build();
        let confident = confident_record("r1", 0, false, false, false, false, Some('E'), true, false);
        let mut non_confident_data = Data::default();
        // no xf tag at all -> not confident, skipped entirely.
        let non_confident = RecordBuf::builder()
            .set_name("r2")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(1).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(vec![b'A'; 10]))
            .set_quality_scores(sam::alignment::record_buf::QualityScores::from(vec![40; 10]))
            .set_data(std::mem::take(&mut non_confident_data))
            .build();

        let bam_records = to_bam_records(&header, &[confident, non_confident]);
        let tags = TagNames::default();
        let stats = mapping_stat(bam_records.into_iter().map(Ok), &header, &tags, "chrM").unwrap();

        assert_eq!(stats.total_alignments, 2);
        assert_eq!(stats.confi_alignments, 1);
        assert_eq!(stats.total_reads_n, 2); // both distinct QNAMEs counted
        assert_eq!(stats.confi_reads_n, 1);
        assert_eq!(stats.chrm_reads, 1);
    }

    #[test]
    fn mapping_stat_excludes_reads_with_no_header_chromosome() {
        // DIV-0019 restructuring fix: upstream's per-chromosome
        // `samfile.fetch(chr_id)` never visits a read whose reference
        // doesn't resolve to any header chromosome (a truly unmapped
        // read with no RNAME) -- it's silently excluded from EVERY
        // count, not just skipped as "unmapped but tallied". Confirmed
        // by grouping the scan by header-chromosome order instead of a
        // flat sequential pass.
        let header = header_with_chrom("chr1", 1000);
        let mapped = confident_record("r1", 0, false, false, false, false, Some('E'), true, false);
        // No reference_sequence_id set at all (stays unset/None) -- a
        // truly unmapped read with no RNAME, unlike a positioned-but-
        // flagged-unmapped mate.
        let mut data = Data::default();
        data.insert(tag("xf"), BufValue::from(1i32));
        let no_rname = RecordBuf::builder()
            .set_name("r2")
            .set_flags(Flags::UNMAPPED)
            .set_data(data)
            .build();

        let bam_records = to_bam_records(&header, &[mapped, no_rname]);
        let tags = TagNames::default();
        let stats = mapping_stat(bam_records.into_iter().map(Ok), &header, &tags, "chrM").unwrap();

        // Only "r1" (on chr1, a real header chromosome) is counted;
        // "r2" (no RNAME at all) is invisible to every per-chromosome
        // fetch and contributes nothing.
        assert_eq!(stats.total_alignments, 1);
        assert_eq!(stats.total_reads_n, 1);
        assert_eq!(stats.confi_alignments, 1);
    }

    #[test]
    fn render_report_percentages_use_correct_denominators() {
        let mut s = MappingStats { total_alignments: 10, confi_alignments: 8, total_reads_n: 10, confi_reads_n: 8, confi_reads_dup: 2, confi_reads_nondup: 6, ..Default::default() };
        s.read_type.insert("Map_consecutively".to_string(), 8);
        let report = render_report(&s).unwrap();
        assert!(report.contains("Total_mapped_reads:\t10\n"));
        // non_confi=2, pct against total(10) -> 20.00%
        assert!(report.contains("Non_confidently_mapped_reads:\t2\t(20.00%)"));
        // confi_reads_dup=2, pct against confi(8) -> 25.00%
        assert!(report.contains("Reads_with_PCR_duplicates:\t2\t(25.00%)"));
        assert!(report.contains("Map_consecutively:\t8\t(100.00%)"));
    }

    #[test]
    fn render_report_errors_on_zero_denominators() {
        let s = MappingStats::default();
        assert!(render_report(&s).is_err());
    }
}
