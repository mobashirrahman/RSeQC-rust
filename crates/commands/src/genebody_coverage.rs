//! Port of `geneBody_coverage.py`: aggregate RNA-seq read coverage
//! across the gene body (5'->3', percentile-binned) for one or more BAM
//! files, plus the skewness-ranked R heatmap/curve plot script.
//! Ported from `oracle/upstream-src/scripts/geneBody_coverage.py` and
//! `qcmodule/mystat.percentile_list`.
//!
//! Reuses `crate::tin::{IndexedRead, build_read_index, genebody_coverage}`
//! for the actual per-position coverage counting: upstream's own
//! `genebody_coverage()` here calls `pysam.AlignmentFile.pileup()` with
//! the SAME unstated defaults documented for `tin.py` (DIV-0011) --
//! `flag_filter` excluding duplicates and `min_base_quality=13` are
//! implemented (via the shared primitive); `ignore_overlaps=True`
//! (mate-pair dedup) is NOT, for the same reason as tin.py.
//!
//! **Preserves several deliberate upstream quirks, do not "fix"**:
//! - `mystat.percentile_list` builds exactly 100 percentile positions
//!   via `(n-1)*i/100.0` linear interpolation, ROUNDED with Python's
//!   round-half-to-even -- NOT the same formula as numpy's percentile
//!   (see `RPKM_saturation.py`'s `percentile_st`) or this port's other
//!   percentile helpers. Only the always->=100-elements branch is
//!   implemented: the CLI enforces `--minimum-length >= 100`, and
//!   `genebody_percentile` only keeps transcripts whose base count is
//!   >= that cutoff, so `percentile_list`'s `len(N) < 100` early return
//!   > is unreachable from this command and intentionally not ported.
//! - Distinct genomic positions among the 100 percentile picks are
//!   DEDUPLICATED (`{position: 0.0 for position in positions}` is a
//!   Python dict keyed by position) before coverage is computed and
//!   aggregated by bin INDEX (not by percentile rank) -- a transcript
//!   short enough that several percentile ranks collapse onto the same
//!   base contributes FEWER than 100 values to `aggregated_coverage`,
//!   silently shifting which bin index later values land in relative to
//!   a transcript with no collisions. This is a genuine upstream
//!   behavior (not obviously intentional), replicated exactly via the
//!   same "dedup to distinct sorted positions" logic already used by
//!   `tin::genebody_coverage`.
//! - `gene_percentiles` is keyed by a composite `gene_id` string
//!   (`chrom_txStart_txEnd_name_strand`); a BED12 file with two lines
//!   producing an IDENTICAL `gene_id` has the SECOND one silently
//!   overwrite the first (Python dict semantics) -- replicated with a
//!   `HashMap` keyed the same way, not a `Vec` that would keep both.
//! - `pearson_moment_coefficient`'s "skewness" centers each value on the
//!   sample's MIDDLE-INDEXED raw value (`values[len(values)/2]`), not
//!   the mean -- an upstream design choice (possibly unconventional for
//!   a "skewness" statistic), preserved exactly, and computed from the
//!   RAW (pre-normalization) coverage values, not the min-max-normalized
//!   ones used for the R plot data.
use std::collections::HashMap;
use std::io::{self, BufRead};
use std::path::{Path, PathBuf};

use noodles_bam as bam;
use noodles_sam as sam;

use crate::python_fmt::{numpy_mean, numpy_std, python_round, python_str_float};
use crate::tin::{self, IndexedRead, build_read_index, genebody_coverage_with_visited};

/// Converts a string into a valid, safe R variable name. Ports
/// `valid_name`.
pub fn valid_name(value: &str) -> String {
    let joined: String = value.split_whitespace().collect::<Vec<_>>().join("_");
    if joined.is_empty() {
        return "sample".to_string();
    }
    let mut chars: Vec<char> = joined.chars().collect();
    if chars[0].is_ascii_digit() {
        chars.insert(0, 'V');
    }
    chars
        .into_iter()
        .map(|c| if c.is_ascii_alphanumeric() || c == '_' || c == '.' { c } else { '_' })
        .collect()
}

/// Ports `mystat.percentile_list` for the always->=100-elements case
/// (see module docs for why the `<100` branch is unreachable here).
/// `sorted_positions` must already be sorted ascending.
pub fn percentile_list(sorted_positions: &[i64]) -> Vec<i64> {
    let n = sorted_positions.len();
    let mut out = Vec::with_capacity(100);
    for i in 1..=100 {
        let k = (n - 1) as f64 * i as f64 / 100.0;
        let f = k.floor();
        let c = k.ceil();
        if f == c {
            out.push(sorted_positions[k as usize]);
        } else {
            let d0 = sorted_positions[f as usize] as f64 * (c - k);
            let d1 = sorted_positions[c as usize] as f64 * (k - f);
            out.push(python_round(d0 + d1));
        }
    }
    out
}

#[derive(Debug, Clone)]
pub struct TranscriptPercentiles {
    pub chrom: String,
    pub strand: String,
    pub positions: Vec<i64>,
}

/// Parses the reference BED12 gene model, keeping transcripts whose
/// total exon base count is `>= mrna_length_cutoff`. Ports
/// `genebody_percentile`. Returns `(transcripts, total_parsed_count)`;
/// `total_parsed_count` includes transcripts later dropped for being
/// too short (matches upstream's `transcript_count` increment point).
pub fn genebody_percentile(reader: impl BufRead, mrna_length_cutoff: i64) -> io::Result<(Vec<TranscriptPercentiles>, i64)> {
    let mut by_gene_id: HashMap<String, TranscriptPercentiles> = HashMap::new();
    let mut order: Vec<String> = Vec::new();
    let mut transcript_count = 0i64;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }

        #[allow(clippy::type_complexity)]
        let parsed: Option<(String, String, String, Vec<i64>, Vec<i64>)> = (|| {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 12 {
                return None;
            }
            let chrom = fields[0].to_string();
            let tx_start: i64 = fields[1].parse().ok()?;
            let tx_end: i64 = fields[2].parse().ok()?;
            let gene_name = fields[3].to_string();
            let strand = fields[5].to_string();
            let gene_id = format!("{chrom}_{tx_start}_{tx_end}_{gene_name}_{strand}");

            let exon_starts: Vec<i64> = fields[11].trim_end_matches(',').split(',').map(|s| s.parse::<i64>().ok().map(|v| v + tx_start)).collect::<Option<_>>()?;
            let exon_sizes: Vec<i64> = fields[10].trim_end_matches(',').split(',').map(|s| s.parse::<i64>().ok()).collect::<Option<_>>()?;
            let exon_ends: Vec<i64> = exon_starts.iter().zip(exon_sizes.iter()).map(|(&s, &sz)| s + sz).collect();

            Some((gene_id, chrom, strand, exon_starts, exon_ends))
        })();

        let Some((gene_id, chrom, strand, exon_starts, exon_ends)) = parsed else {
            eprint!("[NOTE: input BED must be 12-column] skipped this line: {line} ");
            continue;
        };
        transcript_count += 1;

        let mut gene_all_bases = Vec::new();
        for (&s, &e) in exon_starts.iter().zip(exon_ends.iter()) {
            for p in (s + 1)..=e {
                gene_all_bases.push(p);
            }
        }
        if (gene_all_bases.len() as i64) < mrna_length_cutoff {
            continue;
        }

        let positions = percentile_list(&gene_all_bases);
        if by_gene_id.insert(gene_id.clone(), TranscriptPercentiles { chrom, strand, positions }).is_none() {
            order.push(gene_id);
        }
    }

    let out = order.into_iter().filter_map(|id| by_gene_id.remove(&id)).collect();
    Ok((out, transcript_count))
}

/// Aggregates per-bin-index coverage across every transcript for one
/// BAM's already-built read index. Ports the per-BAM-file
/// `genebody_coverage()` orchestration (the per-position counting
/// itself is `tin::genebody_coverage_with_visited`, which returns both
/// coverage values and position-visit markers).
/// A transcript whose `chrom` is not a valid reference in `header` is
/// silently skipped (upstream's `next(samfile.pileup(chrom,1,2),None)`
/// existence probe).
///
/// Returns `(aggregated_coverage, has_float_marker)` where `has_float_marker`
/// indicates which indices received a value from a position that was NEVER
/// visited by any read across ALL transcripts. These positions stay as 0.0
/// (float) in Python due to dict initialization, matching the duck-typing:
/// once 0.0 is added to an int, the aggregate becomes float permanently.
/// Positions that were visited but filtered to 0 (deletions, duplicates, etc.)
/// remain int 0, correctly distinguishing true no-coverage from zero-coverage.
pub fn compute_coverage_for_bam(reads_by_chrom: &HashMap<String, Vec<IndexedRead>>, header: &sam::Header, transcripts: &[TranscriptPercentiles]) -> (Vec<i64>, Vec<bool>) {
    let valid_chroms: std::collections::HashSet<&str> = header.reference_sequences().keys().map(|k| std::str::from_utf8(k).unwrap_or("")).collect();
    let mut aggregated: Vec<i64> = Vec::new();
    let mut ever_visited: Vec<bool> = Vec::new();

    for t in transcripts {
        if t.positions.is_empty() || !valid_chroms.contains(t.chrom.as_str()) {
            continue;
        }
        let empty = Vec::new();
        let reads = reads_by_chrom.get(&t.chrom).unwrap_or(&empty);
        let (mut coverage, mut visited) = genebody_coverage_with_visited(reads, &t.positions, 0.0);
        if t.strand == "-" {
            coverage.reverse();
            visited.reverse();
        }
        if coverage.len() > aggregated.len() {
            aggregated.resize(coverage.len(), 0);
            ever_visited.resize(coverage.len(), false);
        }
        for (i, (v, was_visited)) in coverage.into_iter().zip(visited).enumerate() {
            // An index becomes float if ANY transcript had it unvisited
            // (because dict init with 0.0 causes duck-typing: 0.0 + int → float)
            if !was_visited {
                ever_visited[i] = true;  // Mark as "should be float"
            }
            aggregated[i] += v as i64;
        }
    }

    // Indices that were unvisited in ANY transcript should be marked as float
    let has_float = ever_visited;
    (aggregated, has_float)
}

/// Builds a per-chromosome read index for one BAM, reusing
/// `tin::build_read_index`.
pub fn build_index<I>(records: I, header: &sam::Header) -> io::Result<HashMap<String, Vec<IndexedRead>>>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    build_read_index(records, header)
}

/// Outcome of the windowed coverage driver.
pub enum WindowedCoverage {
    Computed(Vec<i64>, Vec<bool>),
    /// The input was not coordinate-sorted, so a sliding window cannot be
    /// maintained. The caller should fall back to [`compute_coverage_for_bam`]
    /// over a whole-file index.
    NotCoordinateSorted,
}

/// Computes the same aggregate as [`compute_coverage_for_bam`] while holding
/// only a **sliding window** of reads in memory instead of the whole BAM.
///
/// # Why this is exact
///
/// [`compute_coverage_for_bam`] accumulates two things per transcript: a sum
/// (`aggregated[i] += v`) and an OR (`ever_visited[i] |= !was_visited`). Both
/// are commutative and neither reads any other transcript's contribution, so
/// visiting transcripts in a different ORDER cannot change the result. That is
/// what licenses sorting them into coordinate order here, which the whole-file
/// driver is free not to do.
///
/// A transcript only needs reads that overlap `[positions[0]-1,
/// positions[last])`, so a read whose `end` is before the next transcript's
/// window start can be dropped: every remaining transcript starts later, and a
/// read with `start < end` and `end <= window_start` cannot overlap any of
/// them. Reads are pushed in increasing start order and retired from the front,
/// so the retained slice stays sorted by start.
///
/// The per-position counting itself is the same shared
/// `tin::genebody_coverage_with_visited` the whole-file path uses, so this
/// differs only in which reads are resident, never in how a resident read is
/// counted.
pub fn compute_coverage_windowed<I>(
    records: I,
    header: &sam::Header,
    transcripts: &[TranscriptPercentiles],
) -> io::Result<WindowedCoverage>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let valid_chroms: std::collections::HashSet<&str> =
        header.reference_sequences().keys().map(|k| std::str::from_utf8(k).unwrap_or("")).collect();
    let ref_ids: HashMap<&str, usize> = header
        .reference_sequences()
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (std::str::from_utf8(name).unwrap_or(""), i))
        .collect();
    let ref_index = |chrom: &str| -> usize { ref_ids.get(chrom).copied().unwrap_or(usize::MAX) };

    // Coordinate order: by the header's reference-sequence order, then by the
    // transcript's first sampled position. Transcripts on an unknown or
    // invalid chromosome sort last and are skipped, exactly as the whole-file
    // path skips them.
    let mut order: Vec<usize> = (0..transcripts.len())
        .filter(|&i| {
            !transcripts[i].positions.is_empty() && valid_chroms.contains(transcripts[i].chrom.as_str())
        })
        .collect();
    order.sort_by_key(|&i| (ref_index(&transcripts[i].chrom), transcripts[i].positions[0], i));

    let mut aggregated: Vec<i64> = Vec::new();
    let mut ever_visited: Vec<bool> = Vec::new();
    let mut window: Vec<IndexedRead> = Vec::new();
    let mut window_ref: Option<usize> = None;
    let mut closed: Vec<usize> = Vec::new();
    let mut last_start: i64 = i64::MIN;
    let mut records = records.into_iter().peekable();

    for &ti in &order {
        let t = &transcripts[ti];
        let positions = &t.positions;
        let window_start = positions[0] - 1;
        let window_end = positions[positions.len() - 1];
        let t_ref = ref_index(&t.chrom);

        if window_ref != Some(t_ref) {
            if let Some(r) = window_ref {
                closed.push(r);
            }
            window = Vec::new();
            window_ref = Some(t_ref);
            last_start = i64::MIN;
        }

        // Pull every record that starts before this transcript's last sampled
        // position.
        //
        // Records for a LATER reference sequence stay unconsumed for their own
        // transcript to pick up. Records for an EARLIER reference sequence are
        // leftovers from a finished chromosome and must be CONSUMED and
        // discarded: coordinate-sorted BAMs emit whole reference sequences in
        // order, so breaking on them strands the stream on the first leftover
        // and starves every remaining transcript.
        while let Some(peeked) = records.peek() {
            let next = match peeked {
                Ok(record) => record,
                Err(_) => break,
            };
            // A record with NO reference id (BAM tid -1) sorts after every real
            // reference sequence and belongs to no transcript's window, so it
            // must be CONSUMED rather than treated as a stop signal: in a
            // coordinate-sorted BAM these are the trailing unmapped block, and
            // breaking on them would abandon the rest of the stream and leave
            // every later transcript with an empty window. They are dropped by
            // `to_indexed_read`'s unmapped filter anyway.
            let Some(next_ref) = next.reference_sequence_id().transpose()? else {
                records.next();
                continue;
            };
            if header.reference_sequences().get_index(next_ref).is_none() {
                break;
            }
            // Coordinate-sorted order is header-index order, so a record for a
            // strictly later reference sequence means this one is done.
            if next_ref > t_ref {
                break;
            }
            // Strictly earlier: a leftover from a finished chromosome.
            if next_ref < t_ref {
                records.next();
                continue;
            }
            if closed.contains(&next_ref) {
                return Ok(WindowedCoverage::NotCoordinateSorted);
            }
            let Some(pos) = next.alignment_start().transpose()? else { break };
            let start = (pos.get() - 1) as i64;
            if start >= window_end {
                break;
            }
            if start < last_start {
                return Ok(WindowedCoverage::NotCoordinateSorted);
            }
            last_start = start;
            let read = tin::to_indexed_read(next)?;
            records.next();
            if let Some(read) = read {
                window.push(read);
            }
        }

        window.retain(|r| r.end > window_start);

        let (mut coverage, mut visited) = genebody_coverage_with_visited(&window, positions, 0.0);
        if t.strand == "-" {
            coverage.reverse();
            visited.reverse();
        }
        if coverage.len() > aggregated.len() {
            aggregated.resize(coverage.len(), 0);
            ever_visited.resize(coverage.len(), false);
        }
        for (i, (v, was_visited)) in coverage.into_iter().zip(visited).enumerate() {
            if !was_visited {
                ever_visited[i] = true;
            }
            aggregated[i] += v as i64;
        }
    }

    for result in records {
        result?;
    }

    Ok(WindowedCoverage::Computed(aggregated, ever_visited))
}

/// Returns a unique R-safe sample name: `name` itself the first time it
/// occurs in `seen`, else `name.N` where `N` is how many times it has
/// already occurred. Ports `make_unique_sample_name`.
pub fn make_unique_sample_name(name: &str, seen: &[String]) -> String {
    let count = seen.iter().filter(|s| s.as_str() == name).count();
    if count == 0 { name.to_string() } else { format!("{name}.{count}") }
}

/// Checks whether `path` is usable as an upstream `getBamFiles.isbamfile`
/// would: exists, `.bam` extension (case-insensitive), non-zero size,
/// and a `<path>.bai` sidecar present. This is a pure filesystem
/// PRESENCE check -- no BAI content is ever parsed, consistent with the
/// rest of this port's no-BAI-support policy.
pub fn is_bam_file(path: &Path) -> bool {
    if !path.is_file() {
        return false;
    }
    let is_bam_ext = path.extension().map(|e| e.eq_ignore_ascii_case("bam")).unwrap_or(false);
    if !is_bam_ext {
        return false;
    }
    let Ok(meta) = path.metadata() else { return false };
    if meta.len() == 0 {
        eprintln!("The size of {} is 0! Skip it.", path.display());
        return false;
    }
    let mut bai = path.as_os_str().to_owned();
    bai.push(".bai");
    if Path::new(&bai).is_file() {
        true
    } else {
        eprintln!("Warning: {}.bai does not exist! Skip it.", path.display());
        false
    }
}

/// Resolves `-i/--input` into a list of usable BAM paths: a directory
/// (recursively walked), a single BAM file, a newline-separated text
/// file of paths (`#`-prefixed lines skipped), or a comma-separated
/// list. Ports `getBamFiles.get_bam_files`.
pub fn get_bam_files(input: &str) -> Vec<PathBuf> {
    let path = Path::new(input);

    if path.is_dir() {
        let mut out = Vec::new();
        let mut stack = vec![path.to_path_buf()];
        while let Some(dir) = stack.pop() {
            let Ok(entries) = std::fs::read_dir(&dir) else { continue };
            for entry in entries.flatten() {
                let p = entry.path();
                if p.is_dir() {
                    stack.push(p);
                } else if is_bam_file(&p) {
                    out.push(p);
                }
            }
        }
        out.sort();
        return out;
    }

    if is_bam_file(path) {
        return vec![path.to_path_buf()];
    }

    if path.is_file() {
        let mut out = Vec::new();
        if let Ok(text) = std::fs::read_to_string(path) {
            for line in text.lines() {
                let line = line.trim();
                if line.is_empty() || line.starts_with('#') {
                    continue;
                }
                if is_bam_file(Path::new(line)) {
                    out.push(PathBuf::from(line));
                }
            }
        }
        return out;
    }

    input.split(',').filter(|s| is_bam_file(Path::new(s))).map(PathBuf::from).collect()
}

/// Renders the `.geneBodyCoverage.txt` table: header row
/// `Percentile\t1\t2\t...\t100`, then one row per sample.
///
/// The `float_markers` parameter (same length as each sample's values)
/// indicates which indices received 0.0 values during aggregation,
/// matching Python's duck-typing: any value that had 0.0 added to it
/// becomes float, even if the final sum is non-zero (e.g., 15.0).
pub fn render_coverage_txt(samples: &[(String, Vec<i64>, Vec<bool>)]) -> String {
    let mut out = String::new();
    out.push_str("Percentile\t");
    out.push_str(&(1..=100).map(|i| i.to_string()).collect::<Vec<_>>().join("\t"));
    out.push('\n');
    for (name, values, float_markers) in samples {
        out.push_str(name);
        for (v, &is_float) in values.iter().zip(float_markers.iter()) {
            out.push('\t');
            // Preserve Python's duck-typing behavior: any index that received
            // a 0.0 value (marked in float_markers) becomes float, even if
            // subsequent additions make it non-zero.
            if is_float {
                out.push_str(&format!("{v}.0"));
            } else {
                out.push_str(&v.to_string());
            }
        }
        out.push('\n');
    }
    out
}

/// Sample standard deviation (Bessel-corrected, `ddof=1`), matching
/// `numpy.std(values, ddof=1)`.
/// Ports `pearson_moment_coefficient`: standardizes each value against
/// the sample's MIDDLE-INDEXED raw value (not the mean -- see module
/// docs), cubes, and averages.
pub fn pearson_moment_coefficient(values: &[f64]) -> f64 {
    let middle_value = values[values.len() / 2];
    // np.std(ddof=1) and np.mean use numpy's pairwise summation, and
    // `np.float64 ** 3` is C pow(), not repeated multiplication.
    let sigma = numpy_std(values, 1);
    let cubes: Vec<f64> = values.iter().map(|&v| ((v - middle_value) / sigma).powf(3.0)).collect();
    numpy_mean(&cubes)
}

#[derive(Debug, Clone)]
pub struct DatasetEntry {
    pub name: String,
    pub normalized: Vec<f64>,
    pub skewness: f64,
}

/// Min-max normalizes each sample's raw coverage values, computes
/// skewness from the RAW values, and sorts by descending skewness.
/// Ports `load_dataset`, operating directly on the in-memory
/// `(name, raw_values)` pairs already produced by `render_coverage_txt`
/// rather than re-parsing that text back off disk -- a lossless
/// simplification (same class as `RPKM_saturation.py`'s
/// `build_quartile_plot_data`), not a behavioral change.
pub fn load_dataset(samples: &[(String, Vec<i64>, Vec<bool>)]) -> Vec<DatasetEntry> {
    let mut out: Vec<DatasetEntry> = samples
        .iter()
        .map(|(name, raw, _float_markers)| {
            let raw_f: Vec<f64> = raw.iter().map(|&v| v as f64).collect();
            let skewness = pearson_moment_coefficient(&raw_f);
            let min = raw_f.iter().cloned().fold(f64::INFINITY, f64::min);
            let max = raw_f.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let normalized = raw_f.iter().map(|&v| (v - min) / (max - min)).collect();
            DatasetEntry { name: name.clone(), normalized, skewness }
        })
        .collect();
    out.sort_by(|a, b| b.skewness.partial_cmp(&a.skewness).unwrap());
    out
}

/// Renders the R heatmap+curves plotting script. Ports `write_r_code`.
pub fn write_r_code(dataset: &[DatasetEntry], file_prefix: &str, output_format: &str) -> String {
    let fmt = output_format.to_lowercase();
    let names: Vec<&str> = dataset.iter().map(|d| d.name.as_str()).collect();
    let mut out = String::new();

    for d in dataset {
        let values: Vec<String> = d.normalized.iter().map(|&v| python_str_float(v)).collect();
        out.push_str(&format!("{} <- c({})\n", d.name, values.join(",")));
    }

    let tick = "1,10,20,30,40,50,60,70,80,90,100";
    let tick_labels = "\"1\",\"10\",\"20\",\"30\",\"40\",\"50\",\"60\",\"70\",\"80\",\"90\",\"100\"";

    if names.len() >= 3 {
        out.push_str(&format!("data_matrix <- matrix(c({}), byrow=T, ncol=100)\n", names.join(",")));
        out.push_str(&format!("rowLabel <- c({})\n", names.iter().map(|n| format!("\"{n}\"")).collect::<Vec<_>>().join(",")));
        out.push('\n');
        out.push_str(&format!("{fmt}(\"{file_prefix}.heatMap.{fmt}\")\n"));
        out.push_str("rc <- cm.colors(ncol(data_matrix))\n");
        out.push_str(&format!(
            "heatmap(data_matrix, scale=c(\"none\"),keep.dendro=F, labRow=rowLabel,Colv=NA,Rowv=NA,labCol=NA,col=cm.colors(256),margins=c(6,8),ColSideColors=rc,cexRow=1,cexCol=1,xlab=\"Gene body percentile (5'->3')\",add.expr=x_axis_expr <- axis(side=1,at=c({tick}),labels=c({tick_labels})))\n"
        ));
        out.push_str("dev.off()\n");
    }

    out.push('\n');
    out.push_str(&format!("{fmt}(\"{file_prefix}.curves.{fmt}\")\n"));
    out.push_str("x=1:100\n");
    out.push_str(&format!(
        "icolor = colorRampPalette(c(\"#7fc97f\",\"#beaed4\",\"#fdc086\",\"#ffff99\",\"#386cb0\",\"#f0027f\"))({})\n",
        names.len()
    ));

    if names.len() == 1 {
        out.push_str(&format!("plot(x,{},type='l',xlab=\"Gene body percentile (5'->3')\",ylab=\"Coverage\",lwd=0.8,col=icolor[1])\n", names[0]));
    } else if names.len() >= 2 && names.len() <= 6 {
        out.push_str(&format!("plot(x,{},type='l',xlab=\"Gene body percentile (5'->3')\",ylab=\"Coverage\",lwd=0.8,col=icolor[1])\n", names[0]));
        for (i, name) in names.iter().enumerate().skip(1) {
            out.push_str(&format!("lines(x,{name},type='l',col=icolor[{}])\n", i + 1));
        }
        out.push_str(&format!("legend(0,1,fill=icolor[1:{}],legend=c({}))\n", names.len(), names.iter().map(|n| format!("'{n}'")).collect::<Vec<_>>().join(",")));
    } else if names.len() > 6 {
        out.push_str("layout(matrix(c(1,1,1,2,1,1,1,2,1,1,1,2),4,4,byrow=TRUE))\n");
        out.push_str(&format!("plot(x,{},type='l',xlab=\"Gene body percentile (5'->3')\",ylab=\"Coverage\",lwd=0.8,col=icolor[1])\n", names[0]));
        for (i, name) in names.iter().enumerate().skip(1) {
            out.push_str(&format!("lines(x,{name},type='l',col=icolor[{}])\n", i + 1));
        }
        out.push_str("par(mar=c(1,0,2,1))\n");
        out.push_str("plot.new()\n");
        out.push_str(&format!("legend(0,1,fill=icolor[1:{}],legend=c({}))\n", names.len(), names.iter().map(|n| format!("'{n}'")).collect::<Vec<_>>().join(",")));
    }

    out.push_str("dev.off()\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    /// Encodes `RecordBuf`s as BAM records, mirroring the helper each module's
    /// tests carry.
    fn to_bam_records(header: &sam::Header, records: &[sam::alignment::record_buf::RecordBuf]) -> Vec<bam::Record> {
        use sam::alignment::io::Write as _;
        let mut buf = Vec::new();
        {
            let mut writer = bam::io::Writer::new(&mut buf);
            writer.write_header(header).unwrap();
            for record in records {
                writer.write_alignment_record(header, record).unwrap();
            }
        }
        let mut reader = bam::io::Reader::new(buf.as_slice());
        reader.read_header().unwrap();
        reader.records().collect::<Result<Vec<_>, _>>().unwrap()
    }

    #[test]
    fn valid_name_replaces_spaces_and_prefixes_leading_digit() {
        assert_eq!(valid_name("my sample-1"), "my_sample_1");
        assert_eq!(valid_name("1sample"), "V1sample");
        assert_eq!(valid_name(""), "sample");
        assert_eq!(valid_name("plain.name_ok"), "plain.name_ok");
    }

    #[test]
    fn percentile_list_matches_hand_computed_values() {
        // N = [1..=200] (1-indexed values, 0-indexed array), n=200.
        // i=1: k=199*1/100.0=1.99, f=1,c=2 (interpolate).
        //   d0=N[1]*(2-1.99)=2*0.01=0.02, d1=N[2]*(1.99-1)=3*0.99=2.97,
        //   round(2.99)=3.
        // i=50: k=199*50/100.0=99.5, f=99,c=100 (interpolate).
        //   d0=N[99]*(100-99.5)=100*0.5=50, d1=N[100]*(99.5-99)=101*0.5=50.5,
        //   round(100.5)=100 (banker's rounding: rounds to even).
        // i=100: k=199 exactly (f==c) -> N[199]=200.
        let n: Vec<i64> = (1..=200).collect();
        let pl = percentile_list(&n);
        assert_eq!(pl.len(), 100);
        assert_eq!(pl[0], 3);
        assert_eq!(pl[49], 100);
        assert_eq!(pl[99], 200);
    }

    #[test]
    fn genebody_percentile_dedups_by_gene_id_last_wins_and_counts_all_parsed() {
        let bed = "\
chr1\t0\t300\ttx1\t0\t+\t0\t300\t0\t1\t300,\t0,
chr1\t0\t300\ttx1\t0\t+\t0\t300\t0\t1\t300,\t0,
";
        let (transcripts, count) = genebody_percentile(Cursor::new(bed), 100).unwrap();
        assert_eq!(count, 2); // both lines parsed successfully
        assert_eq!(transcripts.len(), 1); // but same gene_id -> dedup to 1
    }

    #[test]
    fn genebody_percentile_skips_too_short_transcripts() {
        let bed = "chr1\t0\t50\ttx1\t0\t+\t0\t50\t0\t1\t50,\t0,\n";
        let (transcripts, count) = genebody_percentile(Cursor::new(bed), 100).unwrap();
        assert_eq!(count, 1);
        assert_eq!(transcripts.len(), 0); // 50bp < 100bp cutoff
    }

    #[test]
    fn make_unique_sample_name_appends_count_suffix() {
        let seen = vec!["a".to_string(), "a".to_string(), "b".to_string()];
        assert_eq!(make_unique_sample_name("a", &seen), "a.2");
        assert_eq!(make_unique_sample_name("b", &seen), "b.1");
        assert_eq!(make_unique_sample_name("c", &seen), "c");
    }

    #[test]
    fn pearson_moment_coefficient_centers_on_middle_value_not_mean() {
        // values = [1,2,3,4,5]; middle_value = values[5/2] = values[2] = 3
        // (0-indexed). sigma = sample std of [1,2,3,4,5] = sqrt(2.5) (ddof=1).
        let values = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let result = pearson_moment_coefficient(&values);
        let sigma = (2.5f64).sqrt();
        let expected: f64 = values.iter().map(|&v| ((v - 3.0) / sigma).powi(3)).sum::<f64>() / 5.0;
        assert!((result - expected).abs() < 1e-12);
    }

    #[test]
    fn load_dataset_normalizes_from_raw_and_sorts_by_descending_skewness() {
        let samples = vec![
            ("low_skew".to_string(), vec![1i64, 2, 3, 4, 5], vec![false; 5]),
            ("high_skew".to_string(), vec![1i64, 1, 1, 1, 100], vec![false; 5]),
        ];
        let dataset = load_dataset(&samples);
        // high_skew should sort first (larger positive skewness from the outlier).
        assert_eq!(dataset[0].name, "high_skew");
        assert_eq!(dataset[1].name, "low_skew");
        // normalized values are min-max in [0,1] with min->0.0 and max->1.0.
        assert_eq!(dataset[1].normalized[0], 0.0);
        assert_eq!(dataset[1].normalized[4], 1.0);
    }

    #[test]
    fn render_coverage_txt_header_and_row_shape() {
        let values = (1..=100).collect::<Vec<i64>>();
        let float_markers = vec![false; 100];
        let samples = vec![("s1".to_string(), values, float_markers)];
        let txt = render_coverage_txt(&samples);
        let mut lines = txt.lines();
        assert!(lines.next().unwrap().starts_with("Percentile\t1\t2\t3"));
        assert!(lines.next().unwrap().starts_with("s1\t1\t2\t3"));
        assert_eq!(txt.lines().count(), 2);
    }

    #[test]
    fn write_r_code_two_sample_curve_branch_exact_text() {
        // Cross-verified byte-for-byte against a `python3 -c` run of the
        // literal upstream write_r_code() expressions with these same
        // inputs (<3 samples -> no heatmap section at all).
        let dataset = vec![
            DatasetEntry { name: "s1".to_string(), normalized: vec![0.0, 0.5, 1.0], skewness: 1.0 },
            DatasetEntry { name: "s2".to_string(), normalized: vec![0.2, 0.4, 0.6], skewness: 0.5 },
        ];
        let script = write_r_code(&dataset, "out", "pdf");
        let expected = "\
s1 <- c(0.0,0.5,1.0)
s2 <- c(0.2,0.4,0.6)

pdf(\"out.curves.pdf\")
x=1:100
icolor = colorRampPalette(c(\"#7fc97f\",\"#beaed4\",\"#fdc086\",\"#ffff99\",\"#386cb0\",\"#f0027f\"))(2)
plot(x,s1,type='l',xlab=\"Gene body percentile (5'->3')\",ylab=\"Coverage\",lwd=0.8,col=icolor[1])
lines(x,s2,type='l',col=icolor[2])
legend(0,1,fill=icolor[1:2],legend=c('s1','s2'))
dev.off()
";
        assert_eq!(script, expected);
    }

    #[test]
    fn write_r_code_three_sample_heatmap_branch_exact_text() {
        // Cross-verified byte-for-byte against a `python3 -c` run of the
        // literal upstream write_r_code() expressions with these same
        // inputs (names a,b,c; single-value normalized vectors).
        let dataset = vec![
            DatasetEntry { name: "a".to_string(), normalized: vec![0.1], skewness: 3.0 },
            DatasetEntry { name: "b".to_string(), normalized: vec![0.2], skewness: 2.0 },
            DatasetEntry { name: "c".to_string(), normalized: vec![0.3], skewness: 1.0 },
        ];
        let script = write_r_code(&dataset, "out", "png");
        let expected_prefix = "\
a <- c(0.1)
b <- c(0.2)
c <- c(0.3)
data_matrix <- matrix(c(a,b,c), byrow=T, ncol=100)
rowLabel <- c(\"a\",\"b\",\"c\")

png(\"out.heatMap.png\")
rc <- cm.colors(ncol(data_matrix))
heatmap(data_matrix, scale=c(\"none\"),keep.dendro=F, labRow=rowLabel,Colv=NA,Rowv=NA,labCol=NA,col=cm.colors(256),margins=c(6,8),ColSideColors=rc,cexRow=1,cexCol=1,xlab=\"Gene body percentile (5'->3')\",add.expr=x_axis_expr <- axis(side=1,at=c(1,10,20,30,40,50,60,70,80,90,100),labels=c(\"1\",\"10\",\"20\",\"30\",\"40\",\"50\",\"60\",\"70\",\"80\",\"90\",\"100\")))
dev.off()

png(\"out.curves.png\")
";
        assert!(script.starts_with(expected_prefix), "actual:\n{script}");
    }

    #[test]
    fn is_bam_file_requires_bai_sidecar() {
        let dir = std::env::temp_dir().join(format!("genebody_cov_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("x.bam");
        std::fs::write(&bam, b"not really bam but nonzero").unwrap();

        assert!(!is_bam_file(&bam)); // no .bai yet
        std::fs::write(dir.join("x.bam.bai"), b"idx").unwrap();
        assert!(is_bam_file(&bam));

        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// Builds a coordinate-sorted record stream plus transcripts spanning
    /// several reference sequences and both strands, so a windowed run has to
    /// advance, retire and switch reference sequence.
    fn windowed_workload() -> (sam::Header, Vec<bam::Record>, Vec<TranscriptPercentiles>) {
        use sam::alignment::record::cigar::op::{Kind, Op};
        use sam::alignment::record_buf::{Cigar, RecordBuf};

        let header = sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                sam::header::record::value::Map::<sam::header::record::value::map::ReferenceSequence>::new(
                    std::num::NonZeroUsize::new(100_000).unwrap(),
                ),
            )
            .add_reference_sequence(
                "chr2",
                sam::header::record::value::Map::<sam::header::record::value::map::ReferenceSequence>::new(
                    std::num::NonZeroUsize::new(100_000).unwrap(),
                ),
            )
            .build();

        let mut bufs: Vec<RecordBuf> = Vec::new();
        let push = |bufs: &mut Vec<RecordBuf>, ref_id: usize, start: usize, n: usize| {
            for i in 0..n {
                let p = start + (i * 4);
                bufs.push(
                    RecordBuf::builder()
                        .set_name(format!("{ref_id}_{i}"))
                        .set_flags(sam::alignment::record::Flags::empty())
                        .set_reference_sequence_id(ref_id)
                        .set_alignment_start(noodles_core::Position::try_from(p + 1).unwrap())
                        .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 10)]))
                        .build(),
                );
            }
        };
        push(&mut bufs, 0, 0, 60);
        push(&mut bufs, 0, 20_000, 60);
        push(&mut bufs, 0, 50_000, 60);
        push(&mut bufs, 1, 1_000, 60);

        let mk = |chrom: &str, strand: &str, from: i64, to: i64| TranscriptPercentiles {
            chrom: chrom.into(),
            strand: strand.into(),
            positions: (from..=to).step_by(10).collect(),
        };
        // Deliberately NOT in coordinate order, and with a mixed strand, to pin
        // that visiting order does not change the aggregate.
        let transcripts = vec![
            mk("chr1", "-", 50_000, 52_000),
            mk("chr2", "+", 1_000, 3_000),
            mk("chr1", "+", 0, 2_000),
            mk("chr1", "-", 20_000, 22_000),
        ];

        (header.clone(), to_bam_records(&header, &bufs), transcripts)
    }

    #[test]
    fn windowed_coverage_matches_whole_file_aggregate() {
        let (header, records, transcripts) = windowed_workload();

        let index = build_index(records.iter().cloned().map(Ok), &header).unwrap();
        let (want_cov, want_float) = compute_coverage_for_bam(&index, &header, &transcripts);

        let got = compute_coverage_windowed(records.iter().cloned().map(Ok), &header, &transcripts).unwrap();
        let WindowedCoverage::Computed(got_cov, got_float) = got else { panic!("expected Computed") };

        assert_eq!(got_cov, want_cov, "aggregate coverage differs");
        assert_eq!(got_float, want_float, "float markers differ");
        assert!(!want_cov.is_empty(), "the fixture must actually produce coverage");
    }

    #[test]
    fn windowed_coverage_detects_an_out_of_order_start() {
        let (header, records, transcripts) = windowed_workload();
        let mut shuffled = records.clone();
        shuffled.swap(5, 9);

        let got = compute_coverage_windowed(shuffled.clone().into_iter().map(Ok), &header, &transcripts).unwrap();
        assert!(matches!(got, WindowedCoverage::NotCoordinateSorted), "expected the unsorted input to be rejected");

        // The whole-file path tolerates it, which is why it is the fallback.
        let index = build_index(shuffled.into_iter().map(Ok), &header).unwrap();
        let _ = compute_coverage_for_bam(&index, &header, &transcripts);
    }

    #[test]
    fn windowed_coverage_skips_transcripts_on_unknown_chromosomes() {
        let (header, records, mut transcripts) = windowed_workload();
        transcripts.push(TranscriptPercentiles {
            chrom: "chrZZ".into(),
            strand: "+".into(),
            positions: vec![1, 500, 1000],
        });

        let index = build_index(records.iter().cloned().map(Ok), &header).unwrap();
        let (want_cov, want_float) = compute_coverage_for_bam(&index, &header, &transcripts);

        let got = compute_coverage_windowed(records.iter().cloned().map(Ok), &header, &transcripts).unwrap();
        let WindowedCoverage::Computed(got_cov, got_float) = got else { panic!("expected Computed") };
        assert_eq!(got_cov, want_cov);
        assert_eq!(got_float, want_float);
    }
}
