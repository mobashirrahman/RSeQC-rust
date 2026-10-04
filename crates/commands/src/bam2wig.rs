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
    /// Accumulated signal per 1-based position, chunked (see `StrandSignal`).
    /// Unstranded and forward strand signal share this field (always
    /// non-negative).
    pub forward: StrandSignal,
    /// Accumulated signal per 1-based position, chunked, always non-positive
    /// (upstream subtracts, not adds, for the reverse strand).
    pub reverse: StrandSignal,
}

/// Positions per dense chunk (B2 memory envelope).
const CHUNK_SIZE: usize = 4096;
/// Words of the per-chunk touched bitmap (`CHUNK_SIZE / 64`).
const CHUNK_WORDS: usize = CHUNK_SIZE / 64;

/// One dense chunk of per-position signal: `CHUNK_SIZE` `f64` values plus a
/// bitmap of which positions were actually touched. A chunk is allocated on
/// first touch; untouched slots (always 0.0) are never iterated or rendered,
/// so a touched position holding 0.0 still renders while an untouched one is
/// skipped, exactly like the per-position map this replaces.
#[derive(Debug, Clone)]
struct Chunk {
    values: Box<[f64; CHUNK_SIZE]>,
    touched: [u64; CHUNK_WORDS],
}

/// A chunk starts sparse and promotes to dense past this many DISTINCT
/// touched positions. Sparse entries cost ~16 bytes each against a 32 KiB
/// dense array, so chunks in thinly covered regions never pay for positions
/// they do not hold; promotion keeps densely covered regions at 8 bytes per
/// position. The trigger is deliberately the distinct-position count, not
/// the touch count: deep RNA-seq coverage re-touches the same positions
/// thousands of times, and promoting on touches would densify chunks whose
/// distinct positions still fit in a few hundred bytes.
///
/// Sparse touches merge lazily: `add` only pushes, and the list is
/// folded (sorted, repeated offsets summed) once `SPARSE_FOLD_EVERY`
/// pushes have arrived since the last fold, promoting to dense when the
/// folded distinct count exceeds `SPARSE_PROMOTE_AT`. Eager per-add
/// merging scans the whole list per touch and costs minutes on real
/// inputs; batching keeps adds amortised O(1). The fold cadence counts
/// pushes (not list length): a chunk with mostly distinct positions must
/// not re-sort on every push once past the threshold.
const SPARSE_PROMOTE_AT: usize = 2048;
/// Fold a sparse touch list every this many pushes since the last fold.
const SPARSE_FOLD_EVERY: usize = 512;

/// Per-chunk storage: a sparse `(offset, value)` touch list until the chunk
/// proves dense, then a dense array plus touched bitmap. Offsets fit `u16`
/// (`CHUNK_SIZE` is 4096); sparse entries are recorded in arrival order with
/// one entry per `add` and merged (summed) on promotion or iteration, so
/// repeated touches of one position accumulate exactly as the old
/// `entry(pos).or_insert(0.0) += delta` did.
#[derive(Debug, Clone)]
enum ChunkData {
    /// `(offset, value)` touches in arrival order (`folded` = list length
    /// right after the last fold) plus the push count logic in `add`.
    Sparse { touches: Vec<(u16, f64)>, since_fold: usize },
    Dense(Box<Chunk>),
}

impl ChunkData {
    /// Fold a sparse touch list into a dense chunk, summing repeated
    /// touches of the same offset.
    fn promote(touches: &[(u16, f64)]) -> Chunk {
        let mut chunk = Chunk {
            values: Box::new([0.0; CHUNK_SIZE]),
            touched: [0; CHUNK_WORDS],
        };
        for &(off, val) in touches {
            let off = off as usize;
            chunk.values[off] += val;
            chunk.touched[off / 64] |= 1u64 << (off % 64);
        }
        chunk
    }

    /// Sort a touch list by offset and sum repeated touches in place,
    /// returning the folded (distinct-offset) list.
    fn fold(touches: &mut Vec<(u16, f64)>) {
        touches.sort_by_key(|&(off, _)| off);
        let mut distinct = 0;
        for i in 0..touches.len() {
            if distinct > 0 && touches[distinct - 1].0 == touches[i].0 {
                let v = touches[i].1;
                touches[distinct - 1].1 += v;
            } else {
                touches[distinct] = touches[i];
                distinct += 1;
            }
        }
        touches.truncate(distinct);
    }
}

/// Per-position coverage signal stored as fixed-size chunks (B2).
///
/// Replaces `BTreeMap<i64, f64>` (one heap node per covered position, ~44
/// bytes per covered base). Each chunk covers `CHUNK_SIZE` consecutive
/// positions and is keyed by chunk index in a `BTreeMap`; chunks start as a
/// sparse touch list and promote to a dense 32 KiB array past
/// `SPARSE_PROMOTE_AT` touches. Iteration yields exactly the touched
/// positions in ascending order with the accumulated values, so every
/// rendered value is unchanged.
#[derive(Debug, Default, Clone)]
pub struct StrandSignal {
    chunks: BTreeMap<i64, ChunkData>,
}

impl StrandSignal {
    /// Chunk index holding 1-based `pos` (positions start at 1; chunk 0
    /// holds positions 1 through `CHUNK_SIZE`).
    fn chunk_index(pos: i64) -> i64 {
        (pos - 1) / CHUNK_SIZE as i64
    }

    /// Offset of 1-based `pos` within its chunk.
    fn offset(pos: i64) -> usize {
        ((pos - 1) % CHUNK_SIZE as i64) as usize
    }

    /// Accumulate `delta` at 1-based `pos`, allocating the chunk on first
    /// touch and promoting it to dense past `SPARSE_PROMOTE_AT` distinct
    /// positions. Sparse touches merge lazily (see the threshold docs), so
    /// an add is an amortised-O(1) push; repeated touches of one offset are
    /// summed at fold time.
    pub fn add(&mut self, pos: i64, delta: f64) {
        let (idx, off) = (Self::chunk_index(pos), Self::offset(pos));
        let data = self.chunks.entry(idx).or_insert_with(|| ChunkData::Sparse {
            touches: Vec::new(),
            since_fold: 0,
        });
        if let ChunkData::Sparse { touches, since_fold } = data {
            touches.push((off as u16, delta));
            *since_fold += 1;
            if *since_fold >= SPARSE_FOLD_EVERY {
                *since_fold = 0;
                ChunkData::fold(touches);
                if touches.len() > SPARSE_PROMOTE_AT {
                    *data = ChunkData::Dense(Box::new(ChunkData::promote(touches)));
                } else if touches.capacity() > 1024 && touches.len() * 4 < touches.capacity() {
                    // Bound retained capacity: repeated folds of a mostly
                    // duplicate list would otherwise pin a large buffer.
                    touches.shrink_to_fit();
                }
            }
            return;
        }
        if let ChunkData::Dense(chunk) = data {
            chunk.values[off] += delta;
            chunk.touched[off / 64] |= 1u64 << (off % 64);
        }
    }

    /// Set an absolute value at 1-based `pos` (used by tests to build
    /// fixtures; production only accumulates via `add`). Recorded as a
    /// touch, so the position renders.
    pub fn insert(&mut self, pos: i64, value: f64) {
        self.add(pos, value - self.get(pos).unwrap_or(0.0));
    }

    /// Value at 1-based `pos` (repeated touches summed), or `None` if
    /// untouched.
    pub fn get(&self, pos: i64) -> Option<f64> {
        let data = self.chunks.get(&Self::chunk_index(pos))?;
        match data {
            ChunkData::Dense(chunk) => {
                let off = Self::offset(pos);
                if chunk.touched[off / 64] & (1u64 << (off % 64)) != 0 {
                    Some(chunk.values[off])
                } else {
                    None
                }
            }
            ChunkData::Sparse { touches, .. } => {
                let off = Self::offset(pos) as u16;
                let mut sum = 0.0;
                let mut found = false;
                for &(o, v) in touches {
                    if o == off {
                        sum += v;
                        found = true;
                    }
                }
                found.then_some(sum)
            }
        }
    }

    /// Number of touched positions. Sparse lists may hold unfolded
    /// repeats, so this folds a copy to count distinct offsets.
    pub fn len(&self) -> usize {
        self.chunks
            .values()
            .map(|data| match data {
                ChunkData::Dense(chunk) => chunk
                    .touched
                    .iter()
                    .map(|w| w.count_ones() as usize)
                    .sum::<usize>(),
                ChunkData::Sparse { touches, .. } => {
                    let mut offs: Vec<u16> =
                        touches.iter().map(|&(o, _)| o).collect();
                    offs.sort_unstable();
                    offs.dedup();
                    offs.len()
                }
            })
            .sum()
    }

    /// True when no position has been touched (chunks are only ever created
    /// by a touching `add`/`insert`, so no chunk means no positions).
    pub fn is_empty(&self) -> bool {
        self.chunks.is_empty()
    }

    /// Touched `(position, value)` pairs in ascending position order
    /// (chunks iterate by ascending index; sparse touches are folded by
    /// offset, repeated touches summed).
    pub fn iter(&self) -> impl Iterator<Item = (i64, f64)> + '_ {
        self.chunks.iter().flat_map(|(&idx, data)| {
            let base = idx * CHUNK_SIZE as i64;
            let folded: Vec<(u16, f64)> = match data {
                ChunkData::Dense(chunk) => (0..CHUNK_SIZE)
                    .filter_map(|off| {
                        if chunk.touched[off / 64] & (1u64 << (off % 64)) != 0 {
                            Some((off as u16, chunk.values[off]))
                        } else {
                            None
                        }
                    })
                    .collect(),
                ChunkData::Sparse { touches, .. } => {
                    let mut sorted = touches.clone();
                    sorted.sort_by_key(|&(off, _)| off);
                    let mut folded: Vec<(u16, f64)> =
                        Vec::with_capacity(sorted.len());
                    for (off, val) in sorted {
                        match folded.last_mut() {
                            Some(last) if last.0 == off => last.1 += val,
                            _ => folded.push((off, val)),
                        }
                    }
                    folded
                }
            };
            folded
                .into_iter()
                .map(move |(off, val)| (base + off as i64 + 1, val))
        })
    }
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
                    entry.forward.add(pos1, 1.0);
                } else {
                    let assigned = *strand_map.get(&key).ok_or_else(|| {
                        io::Error::new(io::ErrorKind::InvalidData, format!("strand_rule has no mapping for computed key {key:?} (rule/pairing-mode mismatch)"))
                    })?;
                    if assigned == '+' {
                        entry.forward.add(pos1, 1.0);
                    }
                    if assigned == '-' {
                        entry.reverse.add(pos1, -1.0);
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
/// entry in `signal` (a listed chromosome with zero coverage). Bytes are
/// written incrementally (B2: the full body is ~0.5 GB on real inputs and
/// must never sit in memory as one `String`); formatting per line is
/// unchanged.
fn render_chrom_block(out: &mut impl io::Write, chrom: &str, values: Option<&StrandSignal>, normalization_factor: Option<f64>) -> io::Result<()> {
    out.write_all(format!("variableStep chrom={chrom}\n").as_bytes())?;
    let Some(values) = values else { return Ok(()) };
    for (pos, value) in values.iter() {
        let scaled = normalization_factor.map(|f| value * f).unwrap_or(value);
        out.write_all(format!("{pos}\t{scaled:.2}\n").as_bytes())?;
    }
    Ok(())
}

/// Renders the unstranded `.wig` file body, streamed into `out` (B2: never
/// materialised as one `String`). Ports the `strandRule`-empty
/// branch of `bamTowig`'s output loop: a chromosome absent from the BAM
/// header entirely is skipped (with a caller-visible warning, not
/// rendered here); every chromosome present in the header always gets a
/// `variableStep` header line, even with zero accumulated positions.
///
/// Takes `signal` by `&mut` and drops each chromosome's entry once its block
/// is rendered (B2): this frees the signal map progressively while output
/// streams to disk, so peak RSS is the signal map alone rather than the map
/// plus a ~0.5 GB output string. Iteration follows `chrom_sizes` order
/// exactly as before, so output bytes are unchanged.
pub fn render_unstranded_wig(out: &mut impl io::Write, chrom_sizes: &[(String, i64)], valid_chroms: &HashSet<String>, signal: &mut HashMap<String, ChromWig>, normalization_factor: Option<f64>) -> io::Result<()> {
    for (chrom, _) in chrom_sizes {
        if !valid_chroms.contains(chrom) {
            continue;
        }
        let values = signal.get(chrom).map(|c| &c.forward);
        render_chrom_block(out, chrom, values, normalization_factor)?;
        signal.remove(chrom);
    }
    Ok(())
}

/// Renders the `(forward_wig, reverse_wig)` file bodies for a
/// strand-specific run, streamed into `fwd`/`rev`. Drains `signal` per
/// chromosome like `render_unstranded_wig` (B2; see its doc comment).
pub fn render_stranded_wig(fwd: &mut impl io::Write, rev: &mut impl io::Write, chrom_sizes: &[(String, i64)], valid_chroms: &HashSet<String>, signal: &mut HashMap<String, ChromWig>, normalization_factor: Option<f64>) -> io::Result<()> {
    for (chrom, _) in chrom_sizes {
        if !valid_chroms.contains(chrom) {
            continue;
        }
        let chrom_wig = signal.get(chrom);
        render_chrom_block(fwd, chrom, chrom_wig.map(|c| &c.forward), normalization_factor)?;
        render_chrom_block(rev, chrom, chrom_wig.map(|c| &c.reverse), normalization_factor)?;
        signal.remove(chrom);
    }
    Ok(())
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
        assert_eq!(chr1.forward.get(101), Some(1.0));
        assert_eq!(chr1.forward.get(105), Some(1.0));
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
        assert_eq!(chr1.reverse.get(1), Some(-1.0));
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

        let mut buf = Vec::new();
        render_unstranded_wig(&mut buf, &chrom_sizes, &valid, &mut signal, None).unwrap();
        assert_eq!(String::from_utf8(buf).unwrap(), "variableStep chrom=chr1\n3\t1.00\n5\t2.00\nvariableStep chrom=chrX\n");
    }

    #[test]
    fn render_unstranded_wig_drops_chrom_not_in_bam_header() {
        let chrom_sizes = vec![("chr1".to_string(), 1000), ("chrUnknown".to_string(), 500)];
        let mut valid = HashSet::new();
        valid.insert("chr1".to_string()); // chrUnknown NOT a valid BAM reference

        let mut buf = Vec::new();
        render_unstranded_wig(&mut buf, &chrom_sizes, &valid, &mut HashMap::new(), None).unwrap();
        assert_eq!(String::from_utf8(buf).unwrap(), "variableStep chrom=chr1\n");
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

        let mut fwd_buf = Vec::new();
        let mut rev_buf = Vec::new();
        render_stranded_wig(&mut fwd_buf, &mut rev_buf, &chrom_sizes, &valid, &mut signal, Some(0.5)).unwrap();
        assert_eq!(String::from_utf8(fwd_buf).unwrap(), "variableStep chrom=chr1\n1\t2.00\n");
        assert_eq!(String::from_utf8(rev_buf).unwrap(), "variableStep chrom=chr1\n1\t-2.00\n");
    }
}
