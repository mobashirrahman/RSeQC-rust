//! Port of `bam2wig.py`: convert a BAM file into WIG coverage file(s).
//! Ports `ParseBAM.bamTowig`/`ParseBAM.calWigSum` from
//! `oracle/upstream-src/src/qcmodule/SAM.py` (lines 2533-2670).
//!
//! No BAI index support (project-wide): upstream iterates
//! `self.samfile.fetch(chrom, 0, chrom_size)` once per chromosome listed
//! in the chrom-size file; since each such fetch only ever returns reads
//! truly mapped to that chromosome, a single sequential pass over the
//! whole BAM grouped by chromosome name produces an identical result
//! (not a behavioral simplification -- the final per-chromosome Fwig/
//! Rwig tables are the same either way).
//!
//! **Preserves a genuine upstream inconsistency, do not "fix"**:
//! `bamTowig`'s `skip_multi` filters on MAPQ (`mapq < q_cut`), exactly
//! like every other command in this port -- but `calWigSum`'s
//! `skip_multi` is a COMPLETELY DIFFERENT mechanism: it inspects the
//! BAM optional tags `H0`, `H1`, `H2`, `IH`, `NH` for any value `> 1`
//! (`ParseBAM.multi_hit_tags`), ignoring MAPQ entirely. These two
//! functions are called back-to-back from the same CLI script with the
//! same `--skip-multi-hits` flag, yet disagree on what "multi-hit"
//! means. Both are replicated exactly, independently.
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{self, BufRead};

use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::record::data::field::Tag};

use rseqc_formats::cigar::fetch_exon_blocks;

pub use crate::fpkm_count::parse_strand_rule;

/// Reads a two-column chromosome-size file, preserving file order
/// (needed: output chromosome order follows this order exactly) and
/// rejecting duplicates / non-positive sizes / malformed lines with a
/// hard error, matching `load_chrom_sizes`.
pub fn load_chrom_sizes(reader: impl BufRead) -> io::Result<Vec<(String, i64)>> {
    let mut out = Vec::new();
    let mut seen = HashSet::new();

    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        let stripped = line.trim();
        if stripped.is_empty() || stripped.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = stripped.split_whitespace().collect();
        if fields.len() < 2 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, format!("line {}: expected at least two columns", line_number + 1)));
        }
        let chrom = fields[0].to_string();
        let size: i64 = fields[1]
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, format!("line {}: chromosome size must be an integer", line_number + 1)))?;
        if size <= 0 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, format!("line {}: chromosome size must be positive", line_number + 1)));
        }
        if !seen.insert(chrom.clone()) {
            return Err(io::Error::new(io::ErrorKind::InvalidData, format!("line {}: duplicate chromosome {chrom:?}", line_number + 1)));
        }
        out.push((chrom, size));
    }

    if out.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "no chromosome sizes were found"));
    }
    Ok(out)
}

#[derive(Debug, Default, Clone)]
pub struct ChromWig {
    /// 1-based position -> accumulated signal. Unstranded and forward
    /// strand signal share this field (always non-negative).
    pub forward: BTreeMap<i64, f64>,
    /// 1-based position -> accumulated signal, always non-positive
    /// (upstream subtracts, not adds, for the reverse strand).
    pub reverse: BTreeMap<i64, f64>,
}

fn read_id_and_strand_key(flags: sam::alignment::record::Flags) -> String {
    let read_id = if flags.is_segmented() {
        if flags.is_first_segment() {
            "1"
        } else if flags.is_last_segment() {
            "2"
        } else {
            ""
        }
    } else {
        ""
    };
    let map_strand = if flags.is_reverse_complemented() { "-" } else { "+" };
    format!("{read_id}{map_strand}")
}

/// Scans every alignment once, accumulating per-base coverage into a
/// `ChromWig` keyed by the read's own reference-sequence name. Ports the
/// per-read body of `bamTowig` (MAPQ-based `skip_multi`, see module
/// docs). A strand-specific run with a `key` missing from `strand_map`
/// replicates upstream's uncaught `KeyError` as a propagated `io::Error`.
pub fn build_wig_signal<I>(
    records: I,
    header: &sam::Header,
    strand_rule_active: bool,
    strand_map: &HashMap<String, char>,
    skip_multi: bool,
    map_qual: u8,
) -> io::Result<HashMap<String, ChromWig>>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut result: HashMap<String, ChromWig> = HashMap::new();

    for item in records {
        let record = item?;
        let flags = record.flags();
        if flags.is_qc_fail() || flags.is_duplicate() || flags.is_secondary() || flags.is_unmapped() {
            continue;
        }
        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        if skip_multi && mapq < map_qual {
            continue;
        }

        let Some(ref_id) = record.reference_sequence_id().transpose()? else { continue };
        let Some((chrom_bstr, _)) = header.reference_sequences().get_index(ref_id) else { continue };
        let chrom = chrom_bstr.to_string();

        let Some(pos) = record.alignment_start().transpose()? else { continue };
        let hit_st = (pos.get() - 1) as i64;

        let key = read_id_and_strand_key(flags);

        let ops: Vec<_> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;
        let blocks = fetch_exon_blocks(hit_st as usize, ops.iter().copied());
        let entry = result.entry(chrom).or_default();

        for (s, e) in blocks {
            for pos1 in (s as i64 + 1)..=(e as i64) {
                if !strand_rule_active {
                    *entry.forward.entry(pos1).or_insert(0.0) += 1.0;
                } else {
                    let assigned = *strand_map.get(&key).ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, format!("strand_rule has no mapping for computed key {key:?} (rule/pairing-mode mismatch)"))
                    })?;
                    if assigned == '+' {
                        *entry.forward.entry(pos1).or_insert(0.0) += 1.0;
                    }
                    if assigned == '-' {
                        *entry.reverse.entry(pos1).or_insert(0.0) -= 1.0;
                    }
                }
            }
        }
    }

    Ok(result)
}

const MULTI_HIT_TAGS: [Tag; 5] = [Tag::new(b'H', b'0'), Tag::new(b'H', b'1'), Tag::new(b'H', b'2'), Tag::TOTAL_HIT_COUNT, Tag::ALIGNMENT_HIT_COUNT];

/// Replicates `calWigSum`'s TAG-based multi-hit check (see module
/// docs): any of `H0`/`H1`/`H2`/`IH`/`NH` present with an integer value
/// `> 1` marks the read as multi-mapped.
fn is_tagged_multi_hit(record: &bam::Record) -> io::Result<bool> {
    let data = record.data();
    for tag in MULTI_HIT_TAGS {
        if let Some(result) = data.get(&tag) {
            if let Some(n) = result?.as_int() {
                if n > 1 {
                    return Ok(true);
                }
            }
        }
    }
    Ok(false)
}

/// Sums exon-block lengths over every alignment whose reference name is
/// in `chrom_names`, applying `calWigSum`'s (tag-based, not MAPQ-based)
/// `skip_multi` filter. Ports `calWigSum`.
pub fn cal_wig_sum<I>(records: I, header: &sam::Header, chrom_names: &HashSet<String>, skip_multi: bool) -> io::Result<f64>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut wigsum = 0.0f64;

    for item in records {
        let record = item?;
        let flags = record.flags();
        if flags.is_qc_fail() || flags.is_duplicate() || flags.is_secondary() || flags.is_unmapped() {
            continue;
        }
        if skip_multi && is_tagged_multi_hit(&record)? {
            continue;
        }

        let Some(ref_id) = record.reference_sequence_id().transpose()? else { continue };
        let Some((chrom_bstr, _)) = header.reference_sequences().get_index(ref_id) else { continue };
        let chrom = chrom_bstr.to_string();
        if !chrom_names.contains(&chrom) {
            continue;
        }

        let Some(pos) = record.alignment_start().transpose()? else { continue };
        let hit_st = (pos.get() - 1) as i64;
        let ops: Vec<_> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;
        for (s, e) in fetch_exon_blocks(hit_st as usize, ops.iter().copied()) {
            wigsum += (e - s) as f64;
        }
    }

    Ok(wigsum)
}

/// Renders one chromosome's `variableStep` block plus its sorted
/// position/value lines (`"%d\t%.2f"`, optionally scaled by
/// `normalization_factor`), or just the header line if `chrom` has no
/// entry in `signal` (a listed chromosome with zero coverage).
fn render_chrom_block(out: &mut String, chrom: &str, values: Option<&BTreeMap<i64, f64>>, normalization_factor: Option<f64>) {
    out.push_str("variableStep chrom=");
    out.push_str(chrom);
    out.push('\n');
    let Some(values) = values else { return };
    for (&pos, &value) in values {
        let scaled = normalization_factor.map(|f| value * f).unwrap_or(value);
        out.push_str(&format!("{pos}\t{scaled:.2}\n"));
    }
}

/// Renders the unstranded `.wig` file body. Ports the `strandRule`-empty
/// branch of `bamTowig`'s output loop: a chromosome absent from the BAM
/// header entirely is skipped (with a caller-visible warning, not
/// rendered here); every chromosome present in the header always gets a
/// `variableStep` header line, even with zero accumulated positions.
pub fn render_unstranded_wig(chrom_sizes: &[(String, i64)], valid_chroms: &HashSet<String>, signal: &HashMap<String, ChromWig>, normalization_factor: Option<f64>) -> String {
    let mut out = String::new();
    for (chrom, _) in chrom_sizes {
        if !valid_chroms.contains(chrom) {
            continue;
        }
        let values = signal.get(chrom).map(|c| &c.forward);
        render_chrom_block(&mut out, chrom, values, normalization_factor);
    }
    out
}

/// Renders the `(forward_wig, reverse_wig)` file bodies for a
/// strand-specific run.
pub fn render_stranded_wig(chrom_sizes: &[(String, i64)], valid_chroms: &HashSet<String>, signal: &HashMap<String, ChromWig>, normalization_factor: Option<f64>) -> (String, String) {
    let mut fwd = String::new();
    let mut rev = String::new();
    for (chrom, _) in chrom_sizes {
        if !valid_chroms.contains(chrom) {
            continue;
        }
        let chrom_wig = signal.get(chrom);
        render_chrom_block(&mut fwd, chrom, chrom_wig.map(|c| &c.forward), normalization_factor);
        render_chrom_block(&mut rev, chrom, chrom_wig.map(|c| &c.reverse), normalization_factor);
    }
    (fwd, rev)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn header_with_chrom(name: &str, len: usize) -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence(name, sam::header::record::value::Map::<sam::header::record::value::map::ReferenceSequence>::new(std::num::NonZeroUsize::new(len).unwrap()))
            .build()
    }

    fn to_bam_records(header: &sam::Header, records: &[sam::alignment::record_buf::RecordBuf]) -> Vec<bam::Record> {
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
    fn build_wig_signal_unstranded_counts_exon_coverage() {
        use sam::alignment::record::{Flags, cigar::Op, cigar::op::Kind};
        use sam::alignment::record_buf::{Cigar, QualityScores, RecordBuf, Sequence};

        let header = header_with_chrom("chr1", 1000);
        let rec = RecordBuf::builder()
            .set_name("r1")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(101).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 5)]))
            .set_sequence(Sequence::from(vec![b'A'; 5]))
            .set_quality_scores(QualityScores::from(vec![40; 5]))
            .build();

        let bam_records = to_bam_records(&header, &[rec]);
        let signal = build_wig_signal(bam_records.into_iter().map(Ok), &header, false, &HashMap::new(), false, 0).unwrap();

        let chr1 = &signal["chr1"];
        // 1-based positions 101..105 inclusive (5 bases), each covered once.
        assert_eq!(chr1.forward.len(), 5);
        assert_eq!(chr1.forward[&101], 1.0);
        assert_eq!(chr1.forward[&105], 1.0);
        assert!(chr1.reverse.is_empty());
    }

    #[test]
    fn build_wig_signal_strand_specific_routes_to_forward_or_reverse() {
        use sam::alignment::record::{Flags, cigar::Op, cigar::op::Kind};
        use sam::alignment::record_buf::{Cigar, QualityScores, RecordBuf, Sequence};

        let header = header_with_chrom("chr1", 1000);
        // Single-end reverse-strand read: read_id="" + map_strand="-" -> key "-".
        let rec = RecordBuf::builder()
            .set_name("r1")
            .set_flags(Flags::REVERSE_COMPLEMENTED)
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(1).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 3)]))
            .set_sequence(Sequence::from(vec![b'A'; 3]))
            .set_quality_scores(QualityScores::from(vec![40; 3]))
            .build();

        let bam_records = to_bam_records(&header, &[rec]);
        let strand_map = HashMap::from([("-".to_string(), '-'), ("+".to_string(), '+')]);
        let signal = build_wig_signal(bam_records.into_iter().map(Ok), &header, true, &strand_map, false, 0).unwrap();

        let chr1 = &signal["chr1"];
        assert!(chr1.forward.is_empty());
        assert_eq!(chr1.reverse.len(), 3);
        assert_eq!(chr1.reverse[&1], -1.0);
    }

    #[test]
    fn cal_wig_sum_excludes_chroms_not_in_the_given_list_and_applies_tag_based_multi_hit_filter() {
        use sam::alignment::record::{Flags, cigar::Op, cigar::op::Kind};
        use sam::alignment::record_buf::{Cigar, QualityScores, RecordBuf, Sequence, data::Data, data::field::Value as BufValue};

        let header = sam::Header::builder()
            .add_reference_sequence("chr1", sam::header::record::value::Map::<sam::header::record::value::map::ReferenceSequence>::new(std::num::NonZeroUsize::new(1000).unwrap()))
            .add_reference_sequence("chr2", sam::header::record::value::Map::<sam::header::record::value::map::ReferenceSequence>::new(std::num::NonZeroUsize::new(1000).unwrap()))
            .build();

        let on_chr1 = RecordBuf::builder()
            .set_name("a")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(1).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .set_sequence(Sequence::from(vec![b'A'; 10]))
            .set_quality_scores(QualityScores::from(vec![40; 10]))
            .build();

        let on_chr2_excluded = RecordBuf::builder()
            .set_name("b")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(1)
            .set_alignment_start(noodles_core::Position::try_from(1).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .set_sequence(Sequence::from(vec![b'A'; 10]))
            .set_quality_scores(QualityScores::from(vec![40; 10]))
            .build();

        let mut multi_hit_data = Data::default();
        multi_hit_data.insert(Tag::ALIGNMENT_HIT_COUNT, BufValue::from(2));
        let multi_hit_on_chr1 = RecordBuf::builder()
            .set_name("c")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::try_from(1).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .set_sequence(Sequence::from(vec![b'A'; 10]))
            .set_quality_scores(QualityScores::from(vec![40; 10]))
            .set_data(multi_hit_data)
            .build();

        let bam_records = to_bam_records(&header, &[on_chr1, on_chr2_excluded, multi_hit_on_chr1]);
        let chrom_names = HashSet::from(["chr1".to_string()]);

        let sum_no_skip = cal_wig_sum(bam_records.clone().into_iter().map(Ok), &header, &chrom_names, false).unwrap();
        // chr2 read always excluded (not in chrom_names); both chr1 reads count when not skipping multi-hits.
        assert_eq!(sum_no_skip, 20.0);

        let sum_skip_multi = cal_wig_sum(bam_records.into_iter().map(Ok), &header, &chrom_names, true).unwrap();
        // NH=2 read excluded by the TAG-based filter -> only the plain chr1 read counts.
        assert_eq!(sum_skip_multi, 10.0);
    }

    #[test]
    fn load_chrom_sizes_preserves_order_and_rejects_duplicates() {
        let text = "chr2\t2000\nchr1\t1000\n";
        let sizes = load_chrom_sizes(Cursor::new(text)).unwrap();
        assert_eq!(sizes, vec![("chr2".to_string(), 2000), ("chr1".to_string(), 1000)]);

        let dup = "chr1\t100\nchr1\t200\n";
        assert!(load_chrom_sizes(Cursor::new(dup)).is_err());

        let non_positive = "chr1\t0\n";
        assert!(load_chrom_sizes(Cursor::new(non_positive)).is_err());
    }

    #[test]
    fn render_unstranded_wig_skips_unlisted_chrom_but_keeps_zero_coverage_header() {
        let chrom_sizes = vec![("chr1".to_string(), 1000), ("chrX".to_string(), 500)];
        let mut valid = HashSet::new();
        valid.insert("chr1".to_string());
        valid.insert("chrX".to_string());

        let mut signal = HashMap::new();
        let mut chr1 = ChromWig::default();
        chr1.forward.insert(5, 2.0);
        chr1.forward.insert(3, 1.0);
        signal.insert("chr1".to_string(), chr1);
        // chrX has zero coverage: absent from `signal` entirely.

        let wig = render_unstranded_wig(&chrom_sizes, &valid, &signal, None);
        assert_eq!(wig, "variableStep chrom=chr1\n3\t1.00\n5\t2.00\nvariableStep chrom=chrX\n");
    }

    #[test]
    fn render_unstranded_wig_drops_chrom_not_in_bam_header() {
        let chrom_sizes = vec![("chr1".to_string(), 1000), ("chrUnknown".to_string(), 500)];
        let mut valid = HashSet::new();
        valid.insert("chr1".to_string()); // chrUnknown NOT a valid BAM reference

        let wig = render_unstranded_wig(&chrom_sizes, &valid, &HashMap::new(), None);
        assert_eq!(wig, "variableStep chrom=chr1\n");
    }

    #[test]
    fn render_applies_normalization_factor_to_both_strands() {
        let chrom_sizes = vec![("chr1".to_string(), 1000)];
        let mut valid = HashSet::new();
        valid.insert("chr1".to_string());

        let mut signal = HashMap::new();
        let mut chr1 = ChromWig::default();
        chr1.forward.insert(1, 4.0);
        chr1.reverse.insert(1, -4.0);
        signal.insert("chr1".to_string(), chr1);

        let (fwd, rev) = render_stranded_wig(&chrom_sizes, &valid, &signal, Some(0.5));
        assert_eq!(fwd, "variableStep chrom=chr1\n1\t2.00\n");
        assert_eq!(rev, "variableStep chrom=chr1\n1\t-2.00\n");
    }
}
