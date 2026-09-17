//! Port of `junction_saturation.py`: assess whether splice-junction
//! discovery has reached sequencing saturation, by repeatedly
//! subsampling shuffled per-read introns and counting cumulative unique
//! known/novel junctions. Algorithm ported from `ParseBAM.
//! saturation_junction` in `oracle/upstream-src/src/qcmodule/SAM.py`
//! (lines 3928-4059) -- the `q_cut`-taking `ParseBAM` variant used by the
//! CLI script, not the `ParseSAM` variant at line 1702.
//!
//! **Preserves several deliberate upstream quirks, do not "fix"**:
//! - The refgene BED's `if int(fields[9] == 1): continue` "skip
//!   single-exon transcripts" check compares a STRING to an int and is
//!   always `False` (dead code) -- never replicated as an actual skip,
//!   same bug as `junction_annotation.py`'s `build_model`.
//! - Only introns on chromosomes that appear in the refgene BED are
//!   considered at all (`if chrom not in chrom_list: continue`) --
//!   reads on chromosomes absent from the gene model contribute nothing,
//!   not even to the "novel junction" counts.
//! - The percentile list is `range(sample_start, sample_end, sample_step)`
//!   with a LITERAL `100` always appended at the end, regardless of what
//!   `sample_end` actually is -- if `sample_end < 100`, the final sampled
//!   point is still exactly 100%, not `sample_end`.
//! - `uniqSpliceSites` accumulates across percentile steps (never reset
//!   between iterations); the "known junction" count requires `count >=
//!   recur` at each step, but the "novel junction" count only requires
//!   the junction NOT being in the known set -- no `recur` threshold
//!   applies there. This asymmetry is preserved exactly, not unified.
//! - Each percentile's bin boundaries are computed from `pertl` and the
//!   fixed `sample_step`, NOT from the previous entry in the percentile
//!   list -- if `sample_start != sample_step` the bins can gap or
//!   overlap. Preserved as a literal per-iteration formula.
//! - `random.shuffle` has no seed anywhere upstream, so the ORIGINAL
//!   command's output already varies run-to-run; this port shuffles with
//!   Rust's `rand` crate (same "uniformly random subsample" semantics,
//!   not bit-identical to Python's Mersenne Twister -- matching that
//!   would require reimplementing CPython's RNG, out of scope; see
//!   compatibility/divergences.yaml).

use std::collections::{HashMap, HashSet};
use std::io::{self, BufRead};

use noodles_bam as bam;
use noodles_sam as sam;
use rand::seq::SliceRandom;
use rseqc_formats::cigar::fetch_intron_blocks;

pub type Junction = (String, i64, i64);

/// Parses the reference BED12 gene model into a known-splice-junction set
/// plus the set of chromosomes it mentions. Robustness tier matches
/// upstream exactly: comment/track/browser lines are skipped; lines with
/// fewer than 12 whitespace-separated fields are skipped with a stderr
/// note; anything else that fails to parse as an integer is an uncaught
/// error, matching upstream's bare `int(...)` calls with no surrounding
/// try/except.
pub fn build_known_splice_sites(reader: impl BufRead) -> io::Result<(HashSet<Junction>, HashSet<String>)> {
    let mut known = HashSet::new();
    let mut chrom_list = HashSet::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 12 {
            eprintln!("Invalid bed line (skipped): {line}");
            continue;
        }

        let chrom = fields[0].to_uppercase();
        chrom_list.insert(chrom.clone());

        let tx_start: i64 = fields[1]
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "invalid txStart"))?;

        let block_sizes: Vec<i64> = fields[10]
            .trim_end_matches(',')
            .split(',')
            .map(|s| s.parse().map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "invalid blockSize")))
            .collect::<io::Result<_>>()?;
        let block_starts: Vec<i64> = fields[11]
            .trim_end_matches(',')
            .split(',')
            .map(|s| s.parse().map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "invalid blockStart")))
            .collect::<io::Result<_>>()?;

        let exon_starts: Vec<i64> = block_starts.iter().map(|&s| s + tx_start).collect();
        let exon_ends: Vec<i64> = exon_starts.iter().zip(block_sizes.iter()).map(|(&s, &sz)| s + sz).collect();

        if exon_starts.len() >= 2 {
            for i in 0..exon_starts.len() - 1 {
                known.insert((chrom.clone(), exon_ends[i], exon_starts[i + 1]));
            }
        }
    }

    Ok((known, chrom_list))
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct SaturationCounts {
    pub percentiles: Vec<i64>,
    pub known: Vec<i64>,
    pub all: Vec<i64>,
    pub novel: Vec<i64>,
}

/// Collects filtered per-read intron junctions (mirrors the upstream
/// read-scan loop) and shuffles them. Split out from the resampling loop
/// so the shuffle step is independently testable/swappable.
pub fn collect_splice_sites<I>(
    records: I,
    header: &sam::Header,
    chrom_list: &HashSet<String>,
    q_cut: u8,
    min_intron: i64,
) -> io::Result<Vec<Junction>>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut sites = Vec::new();

    for result in records {
        let record = result?;
        let flags = record.flags();
        if flags.is_qc_fail() || flags.is_duplicate() || flags.is_secondary() || flags.is_unmapped() {
            continue;
        }
        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        if mapq < q_cut {
            continue;
        }

        let Some(ref_id) = record.reference_sequence_id().transpose()? else { continue };
        let Some((chrom_bstr, _)) = header.reference_sequences().get_index(ref_id) else { continue };
        let chrom = chrom_bstr.to_string().to_uppercase();
        if !chrom_list.contains(&chrom) {
            continue;
        }

        let Some(pos) = record.alignment_start().transpose()? else { continue };
        let hit_st = (pos.get() - 1) as i64;

        let ops: Vec<_> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;
        let introns = fetch_intron_blocks(hit_st as usize, ops.iter().copied());

        for (s, e) in introns {
            let (s, e) = (s as i64, e as i64);
            if e - s < min_intron {
                continue;
            }
            sites.push((chrom.clone(), s, e));
        }
    }

    Ok(sites)
}

/// Runs the cumulative percentile resampling loop over an already-
/// shuffled site list. `sample_end` is used to build the percentile
/// list, but the list always ends with a literal `100` regardless (see
/// module docs).
pub fn compute_saturation(
    shuffled_sites: &[Junction],
    known_sites: &HashSet<Junction>,
    recur: u32,
    sample_start: i64,
    sample_end: i64,
    sample_step: i64,
) -> SaturationCounts {
    let mut percentiles = Vec::new();
    let mut p = sample_start;
    while p < sample_end {
        percentiles.push(p);
        p += sample_step;
    }
    percentiles.push(100);

    let sr_num = shuffled_sites.len() as i64;
    let mut uniq_splice_sites: HashMap<&Junction, u32> = HashMap::new();
    let mut known_junc = Vec::new();
    let mut all_junc = Vec::new();
    let mut unknown_junc = Vec::new();

    for &pertl in &percentiles {
        let mut index_st = (sr_num as f64 * ((pertl - sample_step) as f64 / 100.0)) as i64;
        let index_end = (sr_num as f64 * (pertl as f64 / 100.0)) as i64;
        if index_st < 0 {
            index_st = 0;
        }

        for i in index_st..index_end {
            if let Some(site) = shuffled_sites.get(i as usize) {
                *uniq_splice_sites.entry(site).or_insert(0) += 1;
            }
        }

        all_junc.push(uniq_splice_sites.len() as i64);

        let known_count = uniq_splice_sites
            .iter()
            .filter(|(sj, &count)| known_sites.contains(*sj) && count >= recur)
            .count() as i64;
        known_junc.push(known_count);

        let unknown_count = uniq_splice_sites.keys().filter(|sj| !known_sites.contains(**sj)).count() as i64;
        unknown_junc.push(unknown_count);
    }

    SaturationCounts { percentiles, known: known_junc, all: all_junc, novel: unknown_junc }
}

pub fn shuffle_sites(sites: &mut [Junction], rng: &mut impl rand::Rng) {
    sites.shuffle(rng);
}

/// Ports the literal R-script text written by `saturation_junction`
/// (verified byte-for-byte against a real `python3 -c` run of the
/// upstream `print(...)`/`%`-format lines).
pub fn render_r_script(counts: &SaturationCounts, out_prefix: &str) -> String {
    let join = |v: &[i64]| v.iter().map(i64::to_string).collect::<Vec<_>>().join(",");

    let known_last = counts.known.last().copied().unwrap_or(0) / 1000;
    let all_last = counts.all.last().copied().unwrap_or(0) / 1000;
    let novel_last = counts.novel.last().copied().unwrap_or(0) / 1000;
    let known_first = counts.known.first().copied().unwrap_or(0) / 1000;
    let all_first = counts.all.first().copied().unwrap_or(0) / 1000;
    let novel_first = counts.novel.first().copied().unwrap_or(0) / 1000;

    format!(
        "pdf('{out_prefix}.junctionSaturation_plot.pdf')\n\
x=c({x})\n\
y=c({y})\n\
z=c({z})\n\
w=c({w})\n\
m=max({known_last},{all_last},{novel_last})\n\
n=min({known_first},{all_first},{novel_first})\n\
plot(x,z/1000,xlab='percent of total reads',ylab='Number of splicing junctions (x1000)',type='o',col='blue',ylim=c(n,m))\n\
points(x,y/1000,type='o',col='red')\n\
points(x,w/1000,type='o',col='green')\n\
legend(5,{all_last}, legend=c(\"All junctions\",\"known junctions\", \"novel junctions\"),col=c(\"blue\",\"red\",\"green\"),lwd=1,pch=1)\n\
dev.off()\n",
        x = join(&counts.percentiles),
        y = join(&counts.known),
        z = join(&counts.all),
        w = join(&counts.novel),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn build_known_splice_sites_skips_comments_and_short_lines_dead_branch_never_skips() {
        let bed = "\
track name=x
# comment
chr1\t0\t400\ttx1\t0\t+\t0\t400\t0\t2\t100,100,\t0,300,
too short line
";
        let (known, chrom_list) = build_known_splice_sites(Cursor::new(bed)).unwrap();
        // exon_starts=[0,300], block_sizes=[100,100] -> exon_ends=[100,400];
        // one intron: (exon_ends[0]=100, exon_starts[1]=300).
        assert_eq!(known, HashSet::from([("CHR1".to_string(), 100, 300)]));
        assert_eq!(chrom_list, HashSet::from(["CHR1".to_string()]));
    }

    #[test]
    fn build_known_splice_sites_single_exon_transcript_is_not_skipped() {
        // blockCount field (fields[9]) is "1", which upstream's dead-code
        // check compares as a string to the int 1 -- always False, so a
        // single-exon transcript is still processed (and yields zero
        // introns, not a skip).
        let bed = "chr1\t0\t100\ttx1\t0\t+\t0\t100\t0\t1\t100,\t0,\n";
        let (known, chrom_list) = build_known_splice_sites(Cursor::new(bed)).unwrap();
        assert!(known.is_empty());
        assert_eq!(chrom_list, HashSet::from(["CHR1".to_string()]));
    }

    #[test]
    fn compute_saturation_percentile_list_always_ends_with_literal_100() {
        // sample_end=50 (less than 100) -- the final percentile is still
        // literally 100, not 50.
        let counts = compute_saturation(&[], &HashSet::new(), 1, 5, 50, 10);
        assert_eq!(counts.percentiles, vec![5, 15, 25, 35, 45, 100]);
    }

    #[test]
    fn compute_saturation_cumulative_known_vs_unthresholded_novel() {
        let known_sites = HashSet::from([("CHR1".to_string(), 100, 200)]);
        // Ten sites total: 5 copies of a known junction (recur=2 needs >=2
        // copies within the sampled prefix) and 5 copies of a distinct
        // novel junction, in this fixed (already "shuffled") order so the
        // percentile bins are deterministic for the test.
        let known_j = ("CHR1".to_string(), 100, 200);
        let novel_j = ("CHR1".to_string(), 500, 600);
        let sites: Vec<Junction> =
            std::iter::repeat_n(known_j.clone(), 5).chain(std::iter::repeat_n(novel_j.clone(), 5)).collect();

        // sample_start=50, sample_end=100, sample_step=50 -> percentiles [50, 100].
        let counts = compute_saturation(&sites, &known_sites, 2, 50, 100, 50);
        assert_eq!(counts.percentiles, vec![50, 100]);

        // At 50%: index_st=int(10*0/100)=0, index_end=int(10*0.5)=5 ->
        // covers the 5 known-junction copies only. known count>=2 recur
        // -> 1 known junction; novel set is empty so unknown=0; all=1.
        assert_eq!(counts.all[0], 1);
        assert_eq!(counts.known[0], 1);
        assert_eq!(counts.novel[0], 0);

        // At 100%: index_st=int(10*0.5)=5, index_end=int(10*1.0)=10 ->
        // adds the 5 novel-junction copies (uniqSpliceSites accumulates).
        // all=2 (both distinct junctions seen so far); known stays 1;
        // novel (not thresholded by recur) becomes 1.
        assert_eq!(counts.all[1], 2);
        assert_eq!(counts.known[1], 1);
        assert_eq!(counts.novel[1], 1);
    }

    #[test]
    fn render_r_script_exact_text() {
        // Verified byte-for-byte against a real `python3 -c` run of the
        // upstream print()/%-format lines with the same inputs.
        let counts = SaturationCounts {
            percentiles: vec![5, 10, 100],
            known: vec![1, 2, 3],
            all: vec![2, 4, 6],
            novel: vec![1, 2, 3],
        };
        let script = render_r_script(&counts, "out");
        let expected = "pdf('out.junctionSaturation_plot.pdf')\n\
x=c(5,10,100)\n\
y=c(1,2,3)\n\
z=c(2,4,6)\n\
w=c(1,2,3)\n\
m=max(0,0,0)\n\
n=min(0,0,0)\n\
plot(x,z/1000,xlab='percent of total reads',ylab='Number of splicing junctions (x1000)',type='o',col='blue',ylim=c(n,m))\n\
points(x,y/1000,type='o',col='red')\n\
points(x,w/1000,type='o',col='green')\n\
legend(5,0, legend=c(\"All junctions\",\"known junctions\", \"novel junctions\"),col=c(\"blue\",\"red\",\"green\"),lwd=1,pch=1)\n\
dev.off()\n";
        assert_eq!(script, expected);
    }
}
