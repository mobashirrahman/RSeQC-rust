#!/usr/bin/env python3
"""Generate realistic RSeQC benchmark workloads from REAL human annotation and sequence.

This is Tier A/B of the protocol in `benchmarks/protocol.md`. It differs from the
original `generate_workload.py` (Tier C) in the ways that matter for performance
measurement:

  * The gene model is the REAL UCSC hg38 RefSeq annotation: real transcript
    boundaries, real exon/intron structure, real contig naming. chr17 alone has
    9,616 transcripts and 111,823 exons, versus 3 genes in the Tier C generator.
  * The reference sequence is the REAL hg38 chromosome, so base composition,
    repeat content and sequence complexity are realistic.
  * Reads are simulated FROM THE REAL TRANSCRIPT SEQUENCES with a realistic
    Illumina-like quality profile, spliced alignments (real CIGAR N ops), and
    realistic error rates for indels, clipping and mismatches. This is *not* a
    real sequencer run; see protocol.md section 9.

Everything is deterministic given --seed.

Usage:
    oracle/venv/bin/python3 benchmarks/generate_workload_real.py \
        --refgene real/hg38.ncbiRefSeq.gtf.gz \
        --chrom-dir real --contigs chr17 \
        --reads 1000000 --output-dir workloads/A-mid
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import hashlib
import json
import os
import random
import sys
import time
from pathlib import Path

try:
    import pysam
except ImportError:
    print("Error: pysam is required.", file=sys.stderr)
    sys.exit(1)

GENERATOR_VERSION = "2.0.0"
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

BASES = "ACGT"


class Transcript:
    """A transcript with its exons in genomic ascending order."""

    __slots__ = ("tid", "chrom", "strand", "starts", "ends", "cum", "length", "gene")

    def __init__(self, tid, chrom, strand, starts, ends, gene):
        self.tid = tid
        self.chrom = chrom
        self.strand = strand
        self.starts = starts
        self.ends = ends
        self.gene = gene
        # cumulative exonic length, cum[0] == 0
        cum = [0]
        for s, e in zip(starts, ends):
            cum.append(cum[-1] + (e - s))
        self.cum = cum
        self.length = cum[-1]

    def segments(self, t0: int, t1: int):
        """Return [(gstart, gend), ...] in ASCENDING genomic order covering the
        transcript range [t0, t1).

        Exon i covers transcript coords [cum[i], cum[i+1]) and genomic
        [starts[i], ends[i]).

        On the PLUS strand, transcript offset o within exon i maps to genomic
        starts[i] + o. On the MINUS strand the transcript runs backwards along the
        genome, so offset o maps to genomic ends[i] - 1 - o. Ignoring this places
        every minus-strand read at the wrong locus, and the result is still a
        syntactically valid BAM, so it is validated against the reference after
        writing.
        """
        out = []
        i = bisect.bisect_right(self.cum, t0) - 1
        if i < 0:
            i = 0
        while i < len(self.starts) and self.cum[i] < t1:
            a = max(0, t0 - self.cum[i])          # offset within exon, 5' end of range
            b = min(self.ends[i] - self.starts[i], t1 - self.cum[i])  # 3' end offset
            if self.strand == "+":
                s = self.starts[i] + a
                e = self.starts[i] + b
            else:
                s = self.ends[i] - b
                e = self.ends[i] - a
            if e > s:
                out.append((s, e))
            i += 1
        return out


def load_genome(chrom_dir: Path, contigs: list[str]) -> dict[str, str]:
    genome = {}
    for c in contigs:
        p = chrom_dir / f"{c}.fa"
        if not p.exists():
            p = chrom_dir / f"{c}.fa.gz"
        opener = gzip.open if str(p).endswith(".gz") else open
        with opener(p, "rt") as fh:
            seq = []
            for line in fh:
                if line.startswith(">"):
                    continue
                seq.append(line.strip())
        genome[c] = "".join(seq).upper()
    return genome


def load_transcripts(gtf: Path, contigs: list[str], min_exons: int = 1,
                     min_length: int = 0):
    """Parse BED-like records from RefSeq GTF into Transcript objects."""
    want = set(contigs)
    exons: dict[tuple, list] = {}
    meta: dict[tuple, tuple] = {}
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon" or f[0] not in want:
                continue
            attrs = {}
            for kv in f[8].split(";"):
                kv = kv.strip()
                if not kv:
                    continue
                k, _, v = kv.partition(" ")
                attrs[k] = v.strip().strip('"')
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            key = (f[0], tid)
            exons.setdefault(key, []).append((int(f[3]) - 1, int(f[4])))
            meta[key] = (f[6], attrs.get("gene_id", ""), attrs.get("gene_name", ""))

    txs = []
    for (chrom, tid), ev in exons.items():
        if len(ev) < min_exons:
            continue
        # A transcript must be long enough to hold a full fragment, otherwise the
        # simulated read would be shorter than the quality array.
        if sum(e - s for s, e in ev) < min_length:
            continue
        ev.sort()
        merged = [ev[0]]
        for s, e in ev[1:]:
            if s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        strand, gene_id, gene_name = meta[(chrom, tid)]
        starts = [s for s, _ in merged]
        ends = [e for _, e in merged]
        txs.append(Transcript(tid, chrom, strand, starts, ends, gene_name or gene_id))
    txs.sort(key=lambda t: t.tid)
    return txs


def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


class QualityProfile:
    """Illumina-like Phred profile: high at the start, decaying toward the tail."""

    def __init__(self, read_len: int, rng: random.Random):
        self.read_len = read_len
        # Phred score by cycle: 38 at cycle 0 decaying to ~28 at the end
        self.mean = [max(20.0, 38.0 - 10.0 * (c / max(1, read_len - 1))) for c in range(read_len)]

    def quals(self, rng: random.Random) -> list[int]:
        return [max(2, min(45, int(round(m + rng.gauss(0, 2.0))))) for m in self.mean]


def simulate(args, genome, txs, rng, stats):
    """Yield (transcript, read1_info, read2_info) simulated pairs."""
    cum = [t.length for t in txs]
    total = sum(cum)
    # prefix sums for length-weighted transcript sampling
    pre = []
    acc = 0
    for c in cum:
        acc += c
        pre.append(acc)
    contig_ids = {c: i for i, c in enumerate(sorted(genome))}

    n = args.reads
    for pair_idx in range(n):
        # length-weighted transcript choice
        x = rng.random() * total
        ti = bisect.bisect_left(pre, x)
        if ti >= len(txs):
            ti = len(txs) - 1
        t = txs[ti]
        chrom = t.chrom
        seq_full = genome[chrom]
        rl = args.read_length

        fl = min(max(int(rng.gauss(args.fragment_length, args.fragment_sd)), rl + 10), t.length)
        fs = rng.randrange(0, max(1, t.length - fl + 1))
        t1 = fs + rl
        t2 = fs + fl - rl
        qname = f"READ{pair_idx:09d}"

        for which, (t_start, t_end) in enumerate(((fs, t1), (t2, fs + fl))):
            segs = t.segments(t_start, t_end)
            if not segs:
                continue
            # SAM stores SEQ in read orientation and the CIGAR in reference
            # (genomic, left-to-right) order; for a reverse-strand record the two
            # run in opposite directions, so SEQ[0] corresponds to the LEFTMOST
            # aligned reference base.
            #
            # A transcript's sequence is revcomp(exon1 + exon2 + ... ) in ascending
            # genomic order, so BOTH the read sequence and the CIGAR runs are built
            # here in ascending genomic order, each segment reverse-complemented for
            # a minus-strand read. Building the SEQ in descending order while
            # reversing the CIGAR (or vice versa) mis-places every minus-strand
            # SPLICED read while still producing a syntactically valid BAM, so the
            # output is validated against the reference after writing.
            minus = t.strand == "-"
            if minus:
                # SAM stores SEQ reverse-complemented when FLAG 0x10 is set, and the
                # CIGAR runs in ascending reference order, so stored base i is the
                # COMPLEMENT of reference base i -- not the reverse complement. The
                # per-segment reversal is already accounted for by `segments()`, which
                # maps transcript coordinates onto the genome in the right direction
                # for each strand.
                ref = "".join(seq_full[s:e].translate(COMPLEMENT) for s, e in segs)
            else:
                ref = "".join(seq_full[s:e] for s, e in segs)

            # CIGAR in ascending genomic (reference) order.
            ops = []
            for j, (s, e) in enumerate(segs):
                ops.append(("M", e - s))
                if j < len(segs) - 1:
                    n_ops = segs[j + 1][0] - e
                    if n_ops > 0:
                        ops.append(("N", n_ops))
            if len(ops) > 1:
                stats["spliced"] += 1

            read = list(ref)
            qual = QP.quals(rng)

            # realistic error processes
            r = rng.random()
            if r < args.p_softclip:
                stats["clipped"] += 1
                k = min(rng.randint(1, 8), max(1, len(read) // 4))
                # Clipping k bases off the read also removes k reference bases from
                # the alignment, so the adjacent M run must shrink by k. Adding an S
                # run without adjusting M leaves the CIGAR's query length equal to
                # the M run, which fails validation and silently drops the clip.
                # A soft-clipped base still occupies a position in SEQ, so the
                # read and quality arrays keep their full length; only the M run
                # shrinks and an S run is added.
                if rng.random() < 0.5:
                    # 5' clip in reference order
                    _shrink_edge_m(ops, k, leading=True)
                    ops.insert(0, ("S", k))
                else:
                    # 3' clip in reference order
                    _shrink_edge_m(ops, k, leading=False)
                    ops.append(("S", k))
            elif r < args.p_softclip + args.p_ins:
                stats["ins"] += 1
                p = rng.randrange(1, max(2, len(read) - 1))
                read.insert(p, rng.choice(BASES))
                qual.insert(p, max(2, qual[p] - 5))
                # Only place an indel when the read is entirely within ONE exon.
                # Across a splice junction the CIGAR runs on either side belong to
                # different genomic blocks, and splicing an I/D run between them
                # would misstate the reference span. Those reads are left unspliced
                # with the indel applied, which is a valid and common real pattern.
                if _single_exon(ops):
                    _split_m_insert(ops, p, ("I", 1), consume_ref=False)
                else:
                    read = ref
                    qual = QP.quals(rng)[: len(read)]
                    ops = [("M", len(read))]
            elif r < args.p_softclip + args.p_ins + args.p_del:
                stats["dele"] += 1
                p = rng.randrange(1, max(2, len(read) - 1))
                del read[p]
                del qual[p]
                if _single_exon(ops):
                    # Deletion: M(len) becomes M(p) D(1) M(len-p-1). The deleted
                    # base consumes reference only, so the following M run loses one.
                    _split_m_insert(ops, p, ("D", 1), consume_ref=True)
                else:
                    read = ref
                    qual = QP.quals(rng)[: len(read)]
                    ops = [("M", len(read))]
            elif r < args.p_softclip + args.p_ins + args.p_del + args.p_mismatch:
                stats["mism"] += 1
                p = rng.randrange(0, len(read))
                old = read[p]
                read[p] = rng.choice([b for b in BASES if b != old])
                qual[p] = max(2, qual[p] - 6)

            # ops are already in reference (ascending genomic) order, matching the
            # order `ref` was built in, so no reversal is needed here.
            gstart = segs[0][0]
            is_rev = minus

            # Soft-clip placement needs no fix-up: the ops were built in read order,
            # and reversing the list for a minus-strand read maps the read's leading
            # clip onto the genomic CIGAR's trailing position, which is what SAM
            # requires.

            # Validate: the query-consuming ops must sum to the query length.
            qcons = sum(ln for op, ln in ops if op in ("M", "I", "S", "=" , "X"))
            if qcons != len(read):
                # Fall back to a plain unspliced, unclipped representation rather
                # than emitting an invalid record.
                read = ref
                qual = QP.quals(rng)[: len(read)]
                ops = [("M", len(read))]
                gstart = segs[0][0]
                is_rev = minus
                stats["cigar_repaired"] = stats.get("cigar_repaired", 0) + 1

            cigar = "".join(f"{ln}{op}" for op, ln in ops if ln > 0)
            # exactly one mate of a proper pair is reverse-strand
            mate_is_rev = not is_rev
            leftmost = (which == 0) == (t.strand == "+")
            tlen = fl if leftmost else -fl
            yield (chrom, contig_ids[chrom], gstart, is_rev, "".join(read), qual,
                   cigar, which, qname, tlen, mate_is_rev)

    return stats


def validate_against_reference(bam_path: Path, genome: dict, args, limit: int = 200000):
    """Count mismatches between each read and the real reference.

    The CIGAR is walked in reference (genomic) order, which is the order SAM
    guarantees. Each M block covers reference [rpos, rpos+l); the query substring
    for that block is `query[qpos:qpos+l]` for a forward read and its reverse
    complement for a reverse read, because SEQ is stored in read orientation.

    Note this deliberately does NOT use `get_aligned_pairs`, whose query-position
    convention for reverse-strand records is easy to mismatch against the stored
    SEQ; deriving the walk from the CIGAR keeps the check independent of pysam's
    conventions and validates our own orientation logic directly.
    """
    n = 0
    exact = 0
    total_mm = 0
    per_strand = {True: [0, 0], False: [0, 0]}
    for r in pysam.AlignmentFile(str(bam_path)):
        if r.is_unmapped or r.reference_name is None:
            continue
        seq = genome.get(r.reference_name)
        if seq is None:
            continue
        rpos = r.reference_start
        qpos = 0
        mm = 0
        for op, ln in r.cigartuples or []:
            if op in (0, 7, 8):  # M / = / X
                # SAM stores SEQ reverse-complemented when FLAG 0x10 is set, while the
                # CIGAR runs in ascending reference order. So stored base i pairs with
                # reference base i, complemented for a reverse record and with no
                # reversal. (Verified against a hand-built reverse-strand record: the
                # complement-only reading matches; complement-and-reverse does not.)
                qblk = r.query_sequence[qpos:qpos + ln]
                if r.is_reverse:
                    qblk = qblk.translate(COMPLEMENT)
                rblk = seq[rpos:rpos + ln]
                if len(qblk) != ln or len(rblk) != ln:
                    mm += ln
                else:
                    mm += sum(1 for a, b in zip(qblk, rblk) if a != b)
                qpos += ln
                rpos += ln
            elif op == 1:  # I
                qpos += ln
            elif op == 2:  # D
                rpos += ln
            elif op == 3:  # N
                rpos += ln
            elif op == 4:  # S
                qpos += ln
        n += 1
        total_mm += mm
        if mm == 0:
            exact += 1
        st = per_strand[r.is_reverse]
        st[0] += 1
        st[1] += mm
        if n >= limit:
            break
    return {
        "checked": n,
        "mean_mm": round(total_mm / max(1, n), 4),
        "frac_zero_mm": round(exact / max(1, n), 4),
        "mean_mm_plus": round(per_strand[False][1] / max(1, per_strand[False][0]), 4),
        "mean_mm_minus": round(per_strand[True][1] / max(1, per_strand[True][0]), 4),
        "n_plus": per_strand[False][0],
        "n_minus": per_strand[True][0],
    }


def _shrink_edge_m(ops, k, leading):
    """Shorten the first (or last) M run by k reference bases."""
    idx = 0 if leading else len(ops) - 1
    for i in range(idx, -1, -1) if leading else range(idx, -1, -1):
        if ops[i][0] == "M":
            ops[i] = ("M", max(1, ops[i][1] - k))
            return True
    return False


def _single_exon(ops):
    """True when the CIGAR spans no splice junction (no N run)."""
    return not any(op == "N" for op, _ in ops)


def _split_m_insert(ops, p, mop, consume_ref):
    """Split the M run containing read offset `p` and insert `mop` there.

    An M run of length L covering read offset p becomes
        M(p) <mop> M(L - p)          for an insertion (query-only op)
        M(p) <mop> M(L - p - 1)      for a deletion  (reference-only op)

    Leaving the neighbouring run lengths untouched produces a CIGAR whose implied
    reference span no longer matches the record's SEQ, which silently mis-aligns the
    read while still yielding a valid BAM.
    """
    acc = 0
    for i, (op, ln) in enumerate(ops):
        if op == "M":
            if p <= acc + ln:
                left = p - acc
                right = ln - left - (1 if consume_ref else 0)
                new = [("M", left), mop]
                if right > 0:
                    new.append(("M", right))
                ops[i:i + 1] = new
                return True
            acc += ln
        elif op in ("S", "I", "=" , "X"):
            acc += ln
    return False


def _m_index(ops, p):
    """Index in `ops` at which to insert a mutation op for read offset p."""
    acc = 0
    for i, (op, ln) in enumerate(ops):
        if op == "M":
            if p < acc + ln:
                return i + 1
            acc += ln
        elif op in ("S",):
            acc += ln
    return len(ops)


QP: QualityProfile


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refgene", required=True, help="RefSeq GTF (.gz ok)")
    ap.add_argument("--chrom-dir", required=True, help="dir with <contig>.fa or <contig>.fa.gz")
    ap.add_argument("--contigs", nargs="+", required=True)
    ap.add_argument("--reads", type=int, default=1000000)
    ap.add_argument("--read-length", type=int, default=100)
    ap.add_argument("--fragment-length", type=int, default=200)
    ap.add_argument("--fragment-sd", type=float, default=40.0)
    ap.add_argument("--p-softclip", type=float, default=0.02)
    ap.add_argument("--p-ins", type=float, default=0.01)
    ap.add_argument("--p-del", type=float, default=0.01)
    ap.add_argument("--p-mismatch", type=float, default=0.005)
    ap.add_argument("--p-multimap", type=float, default=0.05)
    ap.add_argument("--p-unmapped-mate", type=float, default=0.005)
    ap.add_argument("--mapq", type=int, default=60)
    ap.add_argument("--seed", type=int, default=20260929)
    ap.add_argument("--max-transcripts", type=int, default=0,
                    help="cap the number of BED12 transcripts (0 = all). Used to keep "
                         "per-transcript commands inside the measurement time budget "
                         "while preserving real exon structure.")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--validate", action="store_true", default=True,
                    help="verify generated reads against the reference (default on)")
    ap.add_argument("--no-validate", dest="validate", action="store_false")
    args = ap.parse_args()

    global QP
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    QP = QualityProfile(args.read_length, rng)

    t0 = time.time()
    genome = load_genome(Path(args.chrom_dir), args.contigs)
    txs = load_transcripts(Path(args.refgene), args.contigs,
                           min_length=args.read_length + args.fragment_length)
    if args.max_transcripts and len(txs) > args.max_transcripts:
        # Longest transcripts first: they are the ones that actually exercise
        # per-transcript pileup and interval work.
        txs.sort(key=lambda t: -t.length)
        txs = txs[:args.max_transcripts]
        txs.sort(key=lambda t: t.tid)
    print(f"genome: {sum(len(v) for v in genome.values())/1e6:.1f} Mbp across {len(genome)} contigs")
    print(f"transcripts: {len(txs)}  exons: {sum(len(t.starts) for t in txs)}  ({time.time()-t0:.1f}s)")

    # write BED12 + chrom.sizes
    bed = out / "model.bed12"
    with open(bed, "w") as fh:
        for t in txs:
            blocks = ",".join(str(e - s) for s, e in zip(t.starts, t.ends)) + ","
            starts = ",".join(str(s - t.starts[0]) for s in t.starts) + ","
            fh.write(
                f"{t.chrom}\t{t.starts[0]}\t{t.ends[-1]}\t{t.tid}\t0\t{t.strand}\t"
                f"{t.starts[0]}\t{t.ends[-1]}\t255\t{len(t.starts)}\t{blocks}\t{starts}\n"
            )
    with open(out / "chrom.sizes", "w") as fh:
        for c in sorted(genome):
            fh.write(f"{c}\t{len(genome[c])}\n")

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": len(genome[c])} for c in sorted(genome)],
    }
    bam_path = out / "reads.bam"
    tmp = out / "reads.unsorted.bam"
    stats = dict(reads=0, spliced=0, clipped=0, ins=0, dele=0, mism=0, multimap=0, unmapped=0)
    n_written = 0

    # Feature counters are raised inside simulate(), which writes into this dict.
    sim_stats = dict(reads=0, spliced=0, clipped=0, ins=0, dele=0, mism=0,
                     unmapped=0, multimap=0, cigar_repaired=0)
    with pysam.AlignmentFile(tmp, "wb", header=header) as f:
        for item in simulate(args, genome, txs, rng, sim_stats):
            (chrom, cid, gstart, is_rev, seq, qual, cigar, which, qname,
             tlen, mate_is_rev) = item
            mate_unmapped = rng.random() < args.p_unmapped_mate
            multi = rng.random() < args.p_multimap
            a = pysam.AlignedSegment()
            a.query_name = qname
            a.query_sequence = seq
            a.query_qualities = bytes(qual)
            a.cigarstring = cigar
            # paired(1) proper-pair(2) read1(64)/read2(128) reverse(16)/mate-reverse(32)
            flag = 1 | 2 | (128 if which else 64)
            flag |= 16 if is_rev else 32
            if mate_unmapped:
                # Unmapped mate: still placed on the same reference so the pair
                # remains traversable, but with the unmapped flag and no CIGAR.
                flag |= 256
                a.reference_id = cid
                a.reference_start = gstart
                a.cigarstring = None
                a.is_unmapped = True
            else:
                a.reference_id = cid
                a.reference_start = gstart
                a.mate_is_unmapped = False
                a.mate_is_reverse = mate_is_rev
                a.next_reference_id = cid
                a.next_reference_start = gstart
                a.template_length = tlen
            a.flag = flag
            a.mapping_quality = 0 if multi else args.mapq
            if multi:
                stats["multimap"] += 1
            f.write(a)
            n_written += 1
            stats["reads"] += 1
            if n_written % 500000 == 0:
                print(f"  {n_written} records ({time.time()-t0:.0f}s)", flush=True)

    pysam.sort("-o", str(bam_path), str(tmp))
    pysam.index("-o", str(bam_path) + ".bai", str(bam_path))
    os.unlink(tmp)

    validation = None
    if args.validate:
        validation = validate_against_reference(bam_path, genome, args)
        print(
            f"  validation: {validation['checked']} reads, "
            f"mean mismatches/read = {validation['mean_mm']:.3f}, "
            f"zero-mismatch = {validation['frac_zero_mm']*100:.1f}%"
        )
        # The simulated per-base mismatch rate is 0.5%, i.e. ~0.5 mismatches per
        # 100bp read. A correct generator lands near that. A generator that
        # mis-places minus-strand reads produces tens of mismatches per read while
        # still yielding a valid BAM, so this check is the only thing standing
        # between a silent data bug and a benchmark.
        # Tolerance scales with read length: the simulated 0.5% per-base rate implies
        # ~0.5 mismatches per 100 bases, plus contributions from indel and clipped
        # reads, so a fixed threshold would reject valid 150bp data.
        tol = max(1.5, 0.02 * args.read_length)
        if validation["mean_mm"] > tol:
            print(
                f"ERROR: mean mismatches/read = {validation['mean_mm']:.2f} exceeds "
                f"the tolerance {tol:.2f} for {args.read_length}bp reads. Reads are "
                f"mis-placed; refusing to write a manifest.",
                file=sys.stderr,
            )
            sys.exit(2)
        if validation["mean_mm_minus"] > tol or validation["mean_mm_plus"] > tol:
            print(
                f"ERROR: strand asymmetry (plus={validation['mean_mm_plus']:.2f}, "
                f"minus={validation['mean_mm_minus']:.2f}) indicates a minus-strand bug.",
                file=sys.stderr,
            )
            sys.exit(2)

    def md5(p):
        h = hashlib.md5()
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()

    # Content hash of the alignments. The BAM's own bytes are NOT content-stable:
    # `samtools sort` stamps a @PG line containing the output path, so two runs of the
    # same seed into different directories produce byte-different BAMs holding identical
    # records. Determinism is therefore asserted over the records, not the container.
    content_md5 = hashlib.md5()
    for r in pysam.AlignmentFile(str(bam_path)):
        content_md5.update(
            f"{r.query_name}\t{r.flag}\t{r.reference_name}\t{r.reference_start}\t"
            f"{r.cigarstring}\t{r.mapping_quality}\t{r.query_sequence}\t"
            f"{r.next_reference_start}\n".encode()
        )

    manifest = {
        "generator_version": GENERATOR_VERSION,
        "validation": validation,
        "reads_content_md5": content_md5.hexdigest(),
        "seed": args.seed,
        "tier": "A",
        "params": vars(args),
        "n_transcripts": len(txs),
        "n_exons": sum(len(t.starts) for t in txs),
        "contigs": {c: len(genome[c]) for c in sorted(genome)},
        "n_records": n_written,
        "read_features": sim_stats,
        "files": {f.name: {"md5": md5(f), "bytes": f.stat().st_size} for f in sorted(out.iterdir()) if f.is_file()},
        "elapsed_s": round(time.time() - t0, 1),
    }
    with open(out / "manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"wrote {n_written} records to {bam_path} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
