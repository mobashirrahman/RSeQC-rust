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
use std::fs::File;
use std::io::{self, BufRead, BufReader};

use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::record::cigar::op::Kind};

#[derive(Debug, Default, Clone)]
pub struct GeneRanges {
    by_chrom: HashMap<String, Vec<(i64, i64, String)>>,
    /// Per-chromosome overlap index, built lazily by `build_index`.
    ///
    /// Upstream stores each chromosome's ranges in a `bx.intervals.Intersecter`
    /// (`qcmodule/SAM.py:2122`), a bitset-backed interval index, so its per-read
    /// overlap query is logarithmic. A linear scan here is asymptotically worse and
    /// dominated the command's cost: with 9,179 real chr17 transcripts and upstream's
    /// 200,000-read sample cap it is ~1.8e9 comparisons, which made this command 6x
    /// SLOWER than the Python reference (measured: 3.9 s of CPU against upstream's
    /// 0.7 s). See benchmarks/RESULTS.generated.md section 6.1.
    index: HashMap<String, ChromIndex>,
}

/// A per-strand overlap index over one chromosome's ranges.
///
/// Only the SET OF DISTINCT STRANDS overlapping a query interval is ever needed
/// (`find_strands` returns a set of strand labels), and there are one or two distinct
/// strands in practice. So instead of reporting every overlapping interval, each
/// strand gets one sorted `starts` array plus a `max_end_before` running maximum. For a
/// query `[s, e)`, the last entry with `start < e` bounds the candidates, and any of
/// those intervals overlaps iff its running max end is `> s`. That is exactly the
/// half-open test `rs < e && s < re` the linear scan performed, so the result is
/// identical while the query becomes O(log n) per strand.
#[derive(Debug, Default, Clone)]
struct ChromIndex {
    /// (strand, sorted start positions, running maximum of end over that prefix)
    strands: Vec<(String, Vec<i64>, Vec<i64>)>,
}

impl ChromIndex {
    fn build(ranges: &[(i64, i64, String)]) -> Self {
        let mut by_strand: HashMap<&str, Vec<(i64, i64)>> = HashMap::new();
        for (s, e, strand) in ranges {
            by_strand.entry(strand.as_str()).or_default().push((*s, *e));
        }
        let mut strands = Vec::with_capacity(by_strand.len());
        for (strand, mut iv) in by_strand {
            iv.sort_unstable();
            let mut starts = Vec::with_capacity(iv.len());
            let mut max_end_before = Vec::with_capacity(iv.len());
            let mut running = i64::MIN;
            for (s, e) in iv {
                starts.push(s);
                running = running.max(e);
                max_end_before.push(running);
            }
            strands.push((strand.to_string(), starts, max_end_before));
        }
        // Deterministic order so the returned set does not depend on HashMap iteration.
        strands.sort_by(|a, b| a.0.cmp(&b.0));
        Self { strands }
    }

    /// True if any interval of this strand overlaps `[s, e)` (half-open).
    fn overlaps(starts: &[i64], max_end_before: &[i64], s: i64, e: i64) -> bool {
        // First index whose start is >= e: nothing at or after it can overlap.
        let hi = starts.partition_point(|&x| x < e);
        if hi == 0 {
            return false;
        }
        // The largest end among the `hi` candidates with start < e.
        max_end_before[hi - 1] > s
    }
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

        let index = by_chrom
            .iter()
            .map(|(c, r)| (c.clone(), ChromIndex::build(r)))
            .collect();
        Ok((Self { by_chrom, index }, skipped))
    }

    /// Distinct strand values of every stored range overlapping
    /// `[start, end)` on `chrom` (half-open interval overlap:
    /// `range_start < end && start < range_end`).
    ///
    /// Uses the per-chromosome overlap index rather than scanning every range, which is
    /// behaviourally identical (the index applies the same half-open test) but
    /// logarithmic instead of linear. Built on first use and reused thereafter, since
    /// `parse` is the only other producer and the map is immutable afterwards.
    pub fn find_strands(&self, chrom: &str, start: i64, end: i64) -> BTreeSet<String> {
        let mut out = BTreeSet::new();
        if let Some(index) = self.index.get(chrom) {
            for (strand, starts, max_end_before) in &index.strands {
                if ChromIndex::overlaps(starts, max_end_before, start, end) {
                    out.insert(strand.clone());
                }
            }
            return out;
        }
        // No index yet: fall back to the linear scan, then index for next time. This
        // keeps `find_strands` correct on a `GeneRanges` built by any path that has not
        // gone through `parse`.
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
    /// Number of usable (strand-classified) reads sampled -- upstream's
    /// `count`, reported via "Total N usable reads were sampled".
    pub sampled_count: u64,
    /// True when the record stream was exhausted (upstream's
    /// `StopIteration` on `next(self.samfile)`) before `sample_size` was
    /// reached -- the CLI layer prints "Finished" in that case only, same
    /// as upstream's `except StopIteration: print("Finished", ...)`
    /// (hitting the sample-size cap instead breaks the loop silently).
    pub stopped_at_eof: bool,
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
    let mut records_iter = records.into_iter();
    let mut stopped_at_eof = false;

    loop {
        if count >= sample_size {
            break;
        }
        // Mirrors upstream's `next(self.samfile)` / `except StopIteration`:
        // check the sample-size cap BEFORE pulling the next record (so a
        // `for`-loop's eager pull-then-check wouldn't consume one extra
        // record past the cap), and treat exhaustion as ending the whole
        // loop rather than just this iteration.
        let Some(result) = records_iter.next() else {
            stopped_at_eof = true;
            break;
        };
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
        sampled_count: count,
        stopped_at_eof,
    })
}

fn get(map: &HashMap<String, u64>, key: &str) -> u64 {
    map.get(key).copied().unwrap_or(0)
}

/// Runs the `infer_experiment.py` CLI body over an already-opened record
/// stream (C1 multi-driver pattern).
///
/// This is exactly what the standalone binary's `run()` does -- same
/// warning, same progress lines with their exact spacing, same stdout
/// report, same error propagation -- except the record source is a
/// caller-supplied iterator and stdout/stderr are caller-supplied sinks.
/// The standalone binary delegates to this (passing the process streams);
/// `rseqc_multi` passes one record broadcast, the shared header, and
/// per-command stream files. `compute_experiment` and the renderer are
/// untouched.
///
/// This is the first command needing the SAM header (`MultiArgs::header`,
/// one `Arc` shared by every worker rather than a per-command clone of the
/// whole reference dictionary) and the first needing a gene model
/// (`MultiArgs::reference_bed`). The BED is opened and parsed HERE rather
/// than by the driver, so a missing or malformed file produces this
/// command's own `infer_experiment.py: error:` text on its stderr stream
/// file, byte-identical to the standalone binary, instead of a
/// driver-level message naming neither the file nor the command.
///
/// Unlike the file-writing commands there is no output-prefix parent check:
/// this one writes no output file, so there is no parent to check. That is
/// upstream's own shape, not an omission.
pub fn run_infer_experiment<I>(
    records: I,
    header: &sam::Header,
    refgene: &std::path::Path,
    sample_size: u64,
    q_cut: u8,
    stdout: &mut dyn io::Write,
    stderr: &mut dyn io::Write,
) -> io::Result<()>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    // Upstream's `validate_args` prints this warning (if any) before doing
    // any real work -- ahead of even opening the refgene BED.
    if sample_size < 1_000 {
        writeln!(
            stderr,
            "Warning: sample size is below 1,000; the inferred protocol may be unreliable."
        )?;
    }

    // `"Reading reference gene model " + refbed + ' ...'` then `end=' '`:
    // one space from the literal's own trailing `...`+space concatenation,
    // no second space (unlike read_quality's "Read BAM file ...  Done").
    write!(
        stderr,
        "Reading reference gene model {} ... ",
        refgene.display()
    )?;
    let (gene_ranges, skipped) = GeneRanges::parse(BufReader::new(File::open(refgene)?))?;
    if skipped > 0 {
        writeln!(stderr, "[NOTE: input bed must be 12-column] skipped {skipped} line(s)")?;
    }
    writeln!(stderr, "Done")?;

    // `"Loading SAM/BAM file ... "` (trailing space in the literal) plus
    // `end=' '` gives two spaces before whatever prints next.
    write!(stderr, "Loading SAM/BAM file ...  ")?;
    let result = compute_experiment(records, header, &gene_ranges, sample_size, q_cut)?;
    if result.stopped_at_eof {
        writeln!(stderr, "Finished")?;
    }
    writeln!(stderr, "Total {} usable reads were sampled", result.sampled_count)?;

    writeln!(stdout, "{}", render_results(&result))?;
    Ok(())
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
        // Regression coverage for stderr progress lines (a kept
        // synthetic-sweep failure showed the CLI printed nothing at all):
        // one usable read was sampled, and the record stream ran out
        // before `sample_size` (200_000) was reached, so upstream prints
        // "Finished" before "Total 1 usable reads were sampled".
        assert_eq!(result.sampled_count, 1);
        assert!(result.stopped_at_eof);
    }

    #[test]
    fn sample_size_cap_stops_before_end_of_file() {
        // When `sample_size` is reached before the record stream is
        // exhausted, upstream's loop breaks via the `count >= sample_size`
        // check (not `StopIteration`), so it must NOT print "Finished".
        let header = test_header();
        let (ranges, _) = GeneRanges::parse("chr1\t0\t1000\tgeneA\t0\t+\n".as_bytes()).unwrap();

        let make_record = || {
            RecordBuf::builder()
                .set_flags(Flags::empty())
                .set_reference_sequence_id(0)
                .set_alignment_start(noodles_core::Position::new(1).unwrap())
                .set_mapping_quality(MappingQuality::new(40).unwrap())
                .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
                .build()
        };

        let bam_records = to_bam_records(&header, &[make_record(), make_record(), make_record()]);
        let result = compute_experiment(bam_records.into_iter().map(Ok), &header, &ranges, 2, 30).unwrap();

        assert_eq!(result.sampled_count, 2);
        assert!(!result.stopped_at_eof);
    }

    #[test]
    fn render_results_exact_text() {
        let r = ExperimentResult {
            protocol: Protocol::SingleEnd,
            spec1: 0.9876,
            spec2: 0.0100,
            undetermined: 0.0024,
            sampled_count: 0,
            stopped_at_eof: false,
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
            sampled_count: 0,
            stopped_at_eof: false,
        };
        let output = render_results(&r);
        let expected = "\nThis is PairEnd Data\nFraction of reads failed to determine: 0.2500\nFraction of reads explained by \"1++,1--,2+-,2-+\": 0.5000\nFraction of reads explained by \"1+-,1-+,2++,2--\": 0.2500";
        assert_eq!(output, expected);
    }

    /// The overlap index must be behaviourally IDENTICAL to the linear scan it
    /// replaced, on every input, not merely on the hand-written fixtures. This is the
    /// property that licenses the optimisation: the index answers the same half-open
    /// question in O(log n) instead of O(n), so any divergence here would silently
    /// change reported strandedness fractions.
    #[test]
    fn find_strands_index_matches_linear_scan_on_random_data() {
        // Deterministic LCG so a failure is reproducible without a dev-dependency.
        let mut state: u64 = 0x2545_F491_4F6C_DD1D;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };

        for trial in 0..200 {
            let n_ranges = 1 + (next() % 60) as usize;
            let mut text = String::new();
            for _ in 0..n_ranges {
                let start = (next() % 400) as i64;
                // Include degenerate and nested ranges: zero-length, fully contained,
                // exactly abutting, and identical duplicates.
                let width = match next() % 4 {
                    0 => 0,
                    1 => 1,
                    2 => (next() % 50) as i64,
                    _ => (next() % 200) as i64,
                };
                let strand = if next() % 2 == 0 { "+" } else { "-" };
                text.push_str(&format!("chr1\t{}\t{}\tg\t0\t{}\n", start, start + width, strand));
            }
            let (ranges, _) = GeneRanges::parse(text.as_bytes()).unwrap();

            // Reference: the original linear scan, computed independently here.
            let linear = |s: i64, e: i64| -> BTreeSet<String> {
                let mut out = BTreeSet::new();
                for (rs, re, strand) in &ranges.by_chrom["chr1"] {
                    if *rs < e && s < *re {
                        out.insert(strand.clone());
                    }
                }
                out
            };

            for _ in 0..40 {
                let s = (next() % 420) as i64;
                let w = (next() % 40) as i64;
                let e = s + w;
                assert_eq!(
                    ranges.find_strands("chr1", s, e),
                    linear(s, e),
                    "trial {trial}: query [{s},{e}) diverged"
                );
            }
            // An unknown contig must stay empty rather than panic.
            assert!(ranges.find_strands("chrZ", 0, 10).is_empty());
        }
    }
}
