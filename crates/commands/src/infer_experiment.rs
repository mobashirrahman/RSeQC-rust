//! Port of `infer_experiment.py`: infer RNA-seq library layout and
//! strandedness. Contract: see `compatibility/commands.yaml` entry
//! `infer_experiment.py`; algorithm ported from `configure_experiment()`
//! in `oracle/upstream-src/src/qcmodule/SAM.py` (lines 2418-2531).
//!
//! Note: upstream's `configure_experiment()` does its own ad-hoc,
//! non-exon-aware gene-range parsing (whole transcript span per line, only
//! columns 0/1/2/5) rather than using `qcmodule.BED.ParseBED`'s full BED12
//! exon parser -- this module deliberately mirrors that (a simple
//! `GeneRanges` type here, not a general BED12 reader). Other Queue C
//! commands that need real exon-level BED12 parsing (read_distribution.py,
//! junction_annotation.py, ...) will need a separate, more complete
//! `crates/formats::bed` module -- don't conflate the two.
//!
//! `readEnd = readStart + qlen` approximates the reference span from the
//! read's aligned length (pysam's `qlen` = sum of M/I/=/X op lengths,
//! excluding soft/hard clips) rather than a true CIGAR-walked reference
//! span (which would also account for N/D ops) -- an upstream
//! simplification, preserved as-is rather than "corrected".
//!
//! SAM-text and CRAM input (DIV-0002/0004) are supported at the CLI layer via
//! `rseqc_formats::open_alignments`; this module's own functions were
//! unaffected (already generic over
//! `IntoIterator<Item = io::Result<bam::Record>>`).

use std::collections::{BTreeSet, HashMap};
use std::io::{self, BufRead};

use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::record::cigar::op::Kind};

#[derive(Debug, Default, Clone)]
pub struct GeneRanges {
    by_chrom: HashMap<String, Vec<(i64, i64, String)>>,
}

impl GeneRanges {
    /// Parses a simple 6+-column (chrom, start, end, name, score, strand)
    /// gene-range file, matching upstream's ad-hoc parsing exactly: lines
    /// starting with `#`, `track`, or `browser` are skipped; any line that
    /// doesn't split into enough whitespace-separated fields, or whose
    /// start/end aren't valid integers, is skipped (upstream: bare
    /// `except:` around the whole parse). Returns the parsed ranges and a
    /// count of skipped lines (upstream logs a stderr note per skip; the
    /// CLI layer can render that from this count).
    pub fn parse(reader: impl BufRead) -> io::Result<(Self, u64)> {
        let mut by_chrom: HashMap<String, Vec<(i64, i64, String)>> = HashMap::new();
        let mut skipped = 0u64;

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
                continue;
            }
            let fields: Vec<&str> = line.split_whitespace().collect();
            let parsed = (|| -> Option<(String, i64, i64, String)> {
                let chrom = fields.first()?.to_string();
                let start: i64 = fields.get(1)?.parse().ok()?;
                let end: i64 = fields.get(2)?.parse().ok()?;
                let strand = fields.get(5)?.to_string();
                Some((chrom, start, end, strand))
            })();

            match parsed {
                Some((chrom, start, end, strand)) => {
                    by_chrom.entry(chrom).or_default().push((start, end, strand));
                }
                None => skipped += 1,
            }
        }

        Ok((Self { by_chrom }, skipped))
    }

    /// Distinct strand values of every stored range overlapping
    /// `[start, end)` on `chrom` (half-open interval overlap:
    /// `range_start < end && start < range_end`).
    pub fn find_strands(&self, chrom: &str, start: i64, end: i64) -> BTreeSet<String> {
        let mut out = BTreeSet::new();
        if let Some(ranges) = self.by_chrom.get(chrom) {
            for (rs, re, strand) in ranges {
                if *rs < end && start < *re {
                    out.insert(strand.clone());
                }
            }
        }
        out
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Protocol {
    PairEnd,
    SingleEnd,
    Mixture,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ExperimentResult {
    pub protocol: Protocol,
    pub spec1: f64,
    pub spec2: f64,
    pub undetermined: f64,
}

pub fn compute_experiment<I>(
    records: I,
    header: &sam::Header,
    gene_ranges: &GeneRanges,
    sample_size: u64,
    q_cut: u8,
) -> io::Result<ExperimentResult>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut count = 0u64;
    let mut p_strandness: HashMap<String, u64> = HashMap::new();
    let mut s_strandness: HashMap<String, u64> = HashMap::new();

    for result in records {
        if count >= sample_size {
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

        let Some(ref_id) = record.reference_sequence_id().transpose()? else {
            continue;
        };
        let Some((chrom, _)) = header.reference_sequences().get_index(ref_id) else {
            continue;
        };
        let chrom = chrom.to_string();

        let Some(pos) = record.alignment_start().transpose()? else {
            continue;
        };
        let read_start = (pos.get() - 1) as i64;

        let ops: Vec<_> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;
        let qlen: i64 = ops
            .iter()
            .filter(|op| matches!(op.kind(), Kind::Match | Kind::Insertion | Kind::SequenceMatch | Kind::SequenceMismatch))
            .map(|op| op.len() as i64)
            .sum();
        let read_end = read_start + qlen;

        let map_strand = if flags.is_reverse_complemented() { '-' } else { '+' };

        if flags.is_segmented() {
            let read_id = if flags.is_first_segment() {
                '1'
            } else if flags.is_last_segment() {
                '2'
            } else {
                continue;
            };
            let strands = gene_ranges.find_strands(&chrom, read_start, read_end);
            if strands.is_empty() {
                continue;
            }
            let strand_from_gene = strands.into_iter().collect::<Vec<_>>().join(":");
            let key = format!("{read_id}{map_strand}{strand_from_gene}");
            *p_strandness.entry(key).or_insert(0) += 1;
            count += 1;
        } else {
            let strands = gene_ranges.find_strands(&chrom, read_start, read_end);
            if strands.is_empty() {
                continue;
            }
            let strand_from_gene = strands.into_iter().collect::<Vec<_>>().join(":");
            let key = format!("{map_strand}{strand_from_gene}");
            *s_strandness.entry(key).or_insert(0) += 1;
            count += 1;
        }
    }

    let (protocol, spec1, spec2, other) = if !p_strandness.is_empty() && s_strandness.is_empty() {
        let total: u64 = p_strandness.values().sum();
        let n1 = get(&p_strandness, "1++") + get(&p_strandness, "1--") + get(&p_strandness, "2+-") + get(&p_strandness, "2-+");
        let n2 = get(&p_strandness, "1+-") + get(&p_strandness, "1-+") + get(&p_strandness, "2++") + get(&p_strandness, "2--");
        let spec1 = n1 as f64 / total as f64;
        let spec2 = n2 as f64 / total as f64;
        (Protocol::PairEnd, spec1, spec2, 1.0 - spec1 - spec2)
    } else if !s_strandness.is_empty() && p_strandness.is_empty() {
        let total: u64 = s_strandness.values().sum();
        let n1 = get(&s_strandness, "++") + get(&s_strandness, "--");
        let n2 = get(&s_strandness, "+-") + get(&s_strandness, "-+");
        let spec1 = n1 as f64 / total as f64;
        let spec2 = n2 as f64 / total as f64;
        (Protocol::SingleEnd, spec1, spec2, 1.0 - spec1 - spec2)
    } else {
        (Protocol::Mixture, 0.0, 0.0, 0.0)
    };

    Ok(ExperimentResult {
        protocol,
        spec1,
        spec2,
        undetermined: other,
    })
}

fn get(map: &HashMap<String, u64>, key: &str) -> u64 {
    map.get(key).copied().unwrap_or(0)
}

/// Matches `print_results()` in oracle/upstream-src/scripts/infer_experiment.py exactly.
pub fn render_results(r: &ExperimentResult) -> String {
    let undetermined = r.undetermined.max(0.0);
    match r.protocol {
        Protocol::PairEnd => format!(
            "\nThis is PairEnd Data\nFraction of reads failed to determine: {:.4}\nFraction of reads explained by \"1++,1--,2+-,2-+\": {:.4}\nFraction of reads explained by \"1+-,1-+,2++,2--\": {:.4}",
            undetermined, r.spec1, r.spec2
        ),
        Protocol::SingleEnd => format!(
            "\nThis is SingleEnd Data\nFraction of reads failed to determine: {:.4}\nFraction of reads explained by \"++,--\": {:.4}\nFraction of reads explained by \"+-,-+\": {:.4}",
            undetermined, r.spec1, r.spec2
        ),
        Protocol::Mixture => "Unknown data type: Mixture".to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::{
        alignment::{
            io::Write as _,
            record::{Flags, MappingQuality, cigar::Op},
            record_buf::{Cigar, RecordBuf},
        },
        header::record::value::{Map, map::ReferenceSequence},
    };
    use std::num::NonZeroUsize;

    fn test_header() -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(10_000).unwrap()),
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
    fn gene_ranges_parses_and_skips_malformed_lines() {
        let text = "#comment\ntrack name=x\nchr1\t100\t200\tgeneA\t0\t+\nchr1\tbad\t200\tgeneB\t0\t+\nchr1\t300\t400\tgeneC\t0\t-\n";
        let (ranges, skipped) = GeneRanges::parse(text.as_bytes()).unwrap();
        assert_eq!(skipped, 1);
        assert_eq!(ranges.find_strands("chr1", 150, 160), BTreeSet::from(["+".to_string()]));
        assert_eq!(ranges.find_strands("chr1", 350, 360), BTreeSet::from(["-".to_string()]));
        assert!(ranges.find_strands("chr1", 250, 260).is_empty());
    }

    #[test]
    fn single_end_forward_strand_read_over_plus_gene_is_plusplus() {
        let header = test_header();
        let (ranges, _) = GeneRanges::parse("chr1\t0\t1000\tgeneA\t0\t+\n".as_bytes()).unwrap();

        let record = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(noodles_core::Position::new(1).unwrap())
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
            .build();

        let bam_records = to_bam_records(&header, &[record]);
        let result = compute_experiment(bam_records.into_iter().map(Ok), &header, &ranges, 200_000, 30).unwrap();

        assert_eq!(result.protocol, Protocol::SingleEnd);
        // Only bucket populated is "++": spec1 should be 1.0 (all reads explained by ++/--).
        assert_eq!(result.spec1, 1.0);
        assert_eq!(result.spec2, 0.0);
    }

    #[test]
    fn render_results_exact_text() {
        let r = ExperimentResult {
            protocol: Protocol::SingleEnd,
            spec1: 0.9876,
            spec2: 0.0100,
            undetermined: 0.0024,
        };
        let output = render_results(&r);
        let expected = "\nThis is SingleEnd Data\nFraction of reads failed to determine: 0.0024\nFraction of reads explained by \"++,--\": 0.9876\nFraction of reads explained by \"+-,-+\": 0.0100";
        assert_eq!(output, expected);
    }

    #[test]
    fn render_results_pairend_exact_text() {
        let r = ExperimentResult {
            protocol: Protocol::PairEnd,
            spec1: 0.5,
            spec2: 0.25,
            undetermined: 0.25,
        };
        let output = render_results(&r);
        let expected = "\nThis is PairEnd Data\nFraction of reads failed to determine: 0.2500\nFraction of reads explained by \"1++,1--,2+-,2-+\": 0.5000\nFraction of reads explained by \"1+-,1-+,2++,2--\": 0.2500";
        assert_eq!(output, expected);
    }
}
