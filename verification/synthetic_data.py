#!/usr/bin/env python3
"""Rich, deterministic synthetic RNA-seq data for verification/synthetic_sweep.py.

Everything is derived from one ``random.Random(seed)`` stream, so the same
``--seed``/``--size`` always produce byte-identical inputs.  Unlike
benchmarks/generate_workload.py (one chromosome, 100M reads, three genes),
this exercises the edge cases RSeQC's commands actually branch on:

* several chromosomes, in a deliberately NON-lexicographic header order
  (chr2, chr1, chr10, chrM, chrX), one of them (chrX) with no genes, and a
  BED chromosome (chr4) absent from the BAM header;
* multi-exon and single-exon, coding and non-coding transcripts on both
  strands (overlapping antisense genes and alternative isoforms included);
* reads that are spliced across annotated introns and novel ones (some
  shorter than junction tools' 50 bp minimum), soft-clipped, hard-clipped,
  with insertions, deletions, mismatches and read ``N`` bases -- MD/NM
  tags are computed from the real reference;
* duplicates (flagged and positional-only), NH>1 multi-mappers with
  secondary alignments, supplementary alignments, low/boundary MAPQs,
  QC-fail, unmapped reads and mates;
* paired-end (proper, improper: wrong orientation / far / cross-chromosome
  / proper-flagged cross-chromosome, overlapping mates, mate unmapped) and
  single-end variants;
* a single-cell BAM with CB/CR/UB/UR/xf/RE/TX/AN/GN tags (some corrected
  barcodes, some missing) including chrM reads;
* two saturation BAMs whose qualifying reads all yield ONE distinct
  junction key / exon-block key, so junction_saturation.py and
  RPKM_saturation.py output is independent of their random.shuffle();
* FASTA/FASTQ (plain, .gz, .bz2) read files, a genome FASTA with
  soft-masked lowercase runs, a transcript FASTA, fixed-length barcode
  FASTA/FASTQ;
* two BigWig tracks with positive and negative values (bedGraph,
  variableStep and fixedStep sections) and a chrom.sizes file;
* a gene-information table plus a mock htseq-count for FPKM-UQ.py.

Usage:
    oracle/venv/bin/python3 verification/synthetic_data.py --seed 1 --size 500 --output-dir DIR
"""
from __future__ import annotations

import argparse
import bz2
import gzip
import os
import random
import stat
from pathlib import Path

import pysam
import pyBigWig

CHROMS = [("chr2", 24000), ("chr1", 30000), ("chr10", 12000), ("chrM", 3000), ("chrX", 6000)]
GENE_PLAN = {"chr2": 6, "chr1": 8, "chr10": 3, "chrM": 1}
READ_LEN = 50
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


class Gen:
    def __init__(self, seed: int, size: int, out: Path):
        self.rng = random.Random(seed)
        self.size = size
        self.out = out
        self.chrom_len = dict(CHROMS)
        self.tid = {c: i for i, (c, _) in enumerate(CHROMS)}
        self.genome = {c: "".join(self.rng.choice("ACGT") for _ in range(n)) for c, n in CHROMS}
        self.transcripts: list[dict] = []
        self.name_n = 0

    # ------------------------------------------------------------------ model
    def make_model(self) -> None:
        rng = self.rng
        for chrom, n in GENE_PLAN.items():
            clen = self.chrom_len[chrom]
            pos = rng.randint(150, 600)
            for g in range(n):
                nex = 1 if rng.random() < 0.25 else rng.randint(2, 5)
                exons = []
                p = pos
                for _ in range(nex):
                    el = rng.randint(60, 400)
                    exons.append((p, p + el))
                    p += el + rng.randint(80, 900)
                if exons[-1][1] >= clen - 200:
                    break
                strand = rng.choice("+-")
                name = f"{chrom}_g{g}"
                self.add_tx(chrom, name, strand, exons)
                if nex >= 3 and rng.random() < 0.5:  # alternative isoform skipping exon 2
                    self.add_tx(chrom, name + "_iso", strand, [exons[0]] + exons[2:])
                if rng.random() < 0.3:  # antisense single-exon overlapping gene
                    a = exons[0][0] + 20
                    self.add_tx(chrom, name + "_as", "-" if strand == "+" else "+",
                                [(a, min(a + rng.randint(120, 300), exons[-1][1]))])
                pos = exons[-1][1] + rng.randint(200, 1500)
        # a BED-only chromosome, absent from every BAM header
        self.add_tx("chr4", "chr4_g0", "+", [(100, 300), (500, 700)])

    def add_tx(self, chrom, name, strand, exons):
        rng = self.rng
        start, end = exons[0][0], exons[-1][1]
        if rng.random() < 0.3:
            cds = (end, end)  # non-coding
        else:
            mrna = [b for s, e in exons for b in range(s, e)]
            i = rng.randint(0, max(0, len(mrna) // 4))
            j = rng.randint(len(mrna) * 3 // 4, len(mrna) - 1)
            cds = (mrna[i], mrna[j] + 1)
        self.transcripts.append(dict(chrom=chrom, name=name, strand=strand, exons=exons,
                                     start=start, end=end, cds=cds))

    def write_model(self) -> None:
        with open(self.out / "model.bed12", "w") as fh:
            lines = []
            for t in sorted(self.transcripts, key=lambda t: (t["chrom"], t["start"])):
                ex = t["exons"]
                lines.append("\t".join(map(str, [
                    t["chrom"], t["start"], t["end"], t["name"], 0, t["strand"], t["cds"][0], t["cds"][1],
                    "0,0,0", len(ex), ",".join(str(e - s) for s, e in ex) + ",",
                    ",".join(str(s - t["start"]) for s, e in ex) + ","])) + "\n")
            fh.write("".join(lines))
        # the same model behind UCSC-style header lines (not every upstream
        # command tolerates them -- exercised by dedicated sweep cases)
        (self.out / "model_hdr.bed12").write_text(
            "browser position chr1:1-1000\ntrack name=synthetic\n# synthetic gene model\n" + "".join(lines))
        with open(self.out / "model.gtf", "w") as fh:
            for t in self.transcripts:
                for k, (s, e) in enumerate(t["exons"]):
                    fh.write(f'{t["chrom"]}\tsyn\texon\t{s + 1}\t{e}\t.\t{t["strand"]}\t.\t'
                             f'gene_id "{t["name"]}"; transcript_id "{t["name"]}.1"; exon_number "{k + 1}";\n')
        # FPKM-UQ gene information + a workload-specific mock htseq-count
        rng = self.rng
        with open(self.out / "genes.info.txt", "w") as fh:
            fh.write("gene_id\tsymbol\tchrom\tstart\tend\tstrand\tgene_type\tx\ty\tz\texon_length\n")
            for i, t in enumerate(self.transcripts):
                gtype = "protein_coding" if i % 3 else rng.choice(["lincRNA", "protein_coding", "miRNA"])
                length = sum(e - s for s, e in t["exons"])
                fh.write(f'G{i}\t{t["name"]}\t{t["chrom"]}\t{t["start"]}\t{t["end"]}\t{t["strand"]}\t'
                         f'{gtype}\ta\tb\tc\t{length}\n')
        with open(self.out / "htseq_counts.txt", "w") as fh:
            for i, _ in enumerate(self.transcripts):
                fh.write(f"G{i}\t{rng.choice([0, 0, 1, 3, 17, 250, rng.randint(0, 5000)])}\n")
            fh.write("GMISSING\t12\n__no_feature\t33\n__ambiguous\t4\n__alignment_not_unique\t9\n")
        mock = self.out / "mock_htseq_count.sh"
        mock.write_text('#!/bin/bash\n# mock htseq-count: ignores its arguments\ncat "$(dirname "$0")/htseq_counts.txt"\n')
        mock.chmod(mock.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    # ------------------------------------------------------------ sequences
    def write_sequences(self) -> None:
        rng = self.rng
        with open(self.out / "genome.fa", "w") as fh:
            for c, _ in CHROMS:
                seq = list(self.genome[c])
                for _ in range(3):  # soft-masked runs
                    a = rng.randint(0, len(seq) - 300)
                    for k in range(a, a + rng.randint(20, 300)):
                        seq[k] = seq[k].lower()
                a = rng.randint(0, len(seq) - 50)
                for k in range(a, a + rng.randint(1, 40)):
                    seq[k] = "N"
                s = "".join(seq)
                fh.write(f">{c} synthetic\n")
                for k in range(0, len(s), 60):
                    fh.write(s[k:k + 60] + "\n")
        with open(self.out / "chrom.sizes", "w") as fh:
            for c, n in CHROMS:
                fh.write(f"{c}\t{n}\n")
        with open(self.out / "mrna.fa", "w") as fh:
            for t in self.transcripts:
                if t["chrom"] not in self.genome:
                    continue
                s = "".join(self.genome[t["chrom"]][a:b] for a, b in t["exons"])
                if t["strand"] == "-":
                    s = revcomp(s)
                fh.write(f">{t['name']}\n{s}\n")

    # ------------------------------------------------------------- alignment
    def qual(self, n: int) -> list[int]:
        rng = self.rng
        base = rng.choice([38, 30, 20])
        return [max(2, min(41, base + rng.randint(-12, 3) - (k // 10))) for k in range(n)]

    def align(self, chrom: str, blocks: list[tuple[int, int]], start: int, qlen: int, *,
              clip5=0, clip3=0, hard5=0, hard3=0, p_mis=0.01, p_ins=0.004, p_del=0.004, p_n=0.002):
        """Walk the reference from ``start`` following ``blocks`` (introns between
        consecutive blocks become N ops) consuming ``qlen`` query bases.
        Returns (pos, cigar, seq, md, nm) or None if it would run off the chromosome."""
        rng = self.rng
        ref = self.genome[chrom]
        clen = len(ref)
        ops: list[list] = []
        seq: list[str] = []

        def emit(op, n=1):
            if ops and ops[-1][0] == op:
                ops[-1][1] += n
            else:
                ops.append([op, n])

        md: list = []  # tokens: int (match run) | str (mismatch base) | '^XYZ'
        nm = 0

        def md_match():
            if md and isinstance(md[-1], int):
                md[-1] += 1
            else:
                md.append(1)

        if hard5:
            emit(5, hard5)
        if clip5:
            emit(4, clip5)
            seq.extend(rng.choice("ACGT") for _ in range(clip5))
        remaining = qlen - clip5 - clip3
        if remaining <= 0:
            return None
        # locate the block containing start
        bi = 0
        while bi < len(blocks) and blocks[bi][1] <= start:
            bi += 1
        if bi >= len(blocks):
            blocks = [(start, clen)]
            bi = 0
        pos = max(start, blocks[bi][0])
        aln_start = pos
        last_op = None
        while remaining > 0:
            if pos >= clen:
                return None
            if bi < len(blocks) - 1 and pos >= blocks[bi][1]:
                emit(3, blocks[bi + 1][0] - pos)
                pos = blocks[bi + 1][0]
                bi += 1
                last_op = 3
                continue
            r = rng.random()
            if last_op == 0 and remaining > 2 and r < p_ins:
                k = min(rng.randint(1, 3), remaining - 1)
                emit(1, k)
                seq.extend(rng.choice("ACGT") for _ in range(k))
                remaining -= k
                nm += k
                last_op = 1
                continue
            if last_op == 0 and remaining > 2 and r < p_ins + p_del:
                k = rng.randint(1, 3)
                if (bi < len(blocks) - 1 and pos + k >= blocks[bi][1]) or pos + k >= clen:
                    pass
                else:
                    emit(2, k)
                    md.append("^" + ref[pos:pos + k])
                    pos += k
                    nm += k
                    last_op = 2
                    continue
            b = ref[pos]
            r2 = rng.random()
            if r2 < p_n:
                seq.append("N")
                nm += 1
                md.append(b)
            elif r2 < p_n + p_mis:
                seq.append(rng.choice([x for x in "ACGT" if x != b]))
                nm += 1
                md.append(b)
            else:
                seq.append(b)
                md_match()
            emit(0, 1)
            pos += 1
            remaining -= 1
            last_op = 0
        if clip3:
            emit(4, clip3)
            seq.extend(rng.choice("ACGT") for _ in range(clip3))
        if hard3:
            emit(5, hard3)
        # render MD with mandatory numbers between non-number tokens
        parts = []
        prev_num = False
        for t in md:
            if isinstance(t, int):
                parts.append(str(t))
                prev_num = True
            else:
                if not prev_num:
                    parts.append("0")
                parts.append(t)
                prev_num = False
        if not prev_num:
            parts.append("0")
        return aln_start, [(o, n) for o, n in ops], "".join(seq), "".join(parts), nm

    def tx_blocks_from(self, t, mrna_offset: int):
        """Genomic start of transcript-coordinate ``mrna_offset`` (5'->3' in genome order)."""
        off = mrna_offset
        for s, e in t["exons"]:
            if off < e - s:
                return s + off
            off -= e - s
        return None

    def new_name(self, prefix="r") -> str:
        self.name_n += 1
        return f"{prefix}{self.name_n:07d}"

    # -------------------------------------------------------- fragment model
    def pick_location(self):
        """Returns (chrom, blocks, left_start, frag_len, gene_strand)."""
        rng = self.rng
        r = rng.random()
        tx = [t for t in self.transcripts if t["chrom"] in self.genome]
        if r < 0.72:
            t = rng.choice(tx)
            mlen = sum(e - s for s, e in t["exons"])
            frag = rng.randint(60, 420)
            if mlen <= READ_LEN + 5:
                return None
            frag = min(frag, mlen)
            a = rng.randint(0, mlen - frag)
            left = self.tx_blocks_from(t, a)
            right_off = a + frag - READ_LEN
            right = self.tx_blocks_from(t, max(a, right_off))
            return t["chrom"], t["exons"], left, right, t["strand"]
        chrom = rng.choice(["chr2", "chr1", "chr10", "chrX", "chrM"])
        clen = self.chrom_len[chrom]
        left = rng.randint(0, clen - 800)
        if r < 0.86:  # unspliced genomic (intronic/intergenic)
            frag = rng.randint(60, 500)
            return chrom, [(0, clen)], left, left + frag - READ_LEN, None
        # novel junction (sometimes shorter than the 50bp minimum intron)
        x = rng.randint(10, READ_LEN - 10)
        gap = rng.choice([rng.randint(20, 49), rng.randint(50, 600)])
        blocks = [(left, left + x), (left + x + gap, clen)]
        frag = rng.randint(80, 400)
        return chrom, blocks, left, left + x + gap + frag - READ_LEN, None

    def read_mods(self):
        rng = self.rng
        m = {}
        if rng.random() < 0.08:
            m["clip5"] = rng.randint(1, 10)
        if rng.random() < 0.08:
            m["clip3"] = rng.randint(1, 10)
        if rng.random() < 0.02:
            m["hard5"] = rng.randint(1, 5)
        if rng.random() < 0.3:
            m.update(p_ins=0.03, p_del=0.03)
        return m

    def mapq(self) -> int:
        return self.rng.choice([255, 255, 255, 60, 60, 50, 40, 30, 30, 29, 10, 3, 1, 0])

    def seg(self, name, chrom, aln, flag, mapq, qual=None, tags=()):
        pos, cigar, seq, md, nm = aln
        a = pysam.AlignedSegment(self.header_obj)
        a.query_name = name
        a.flag = flag
        a.reference_id = self.tid[chrom]
        a.reference_start = pos
        a.mapping_quality = mapq
        a.cigartuples = cigar
        if flag & 0x10:
            # stored SEQ is reverse-complemented relative to the read, but we
            # simply keep the reference-strand sequence as aligners do.
            pass
        a.query_sequence = seq
        a.query_qualities = qual if qual is not None else self.qual(len(seq))
        a.set_tags([("NM", nm), ("MD", md)] + list(tags))
        return a

    def ref_end(self, aln):
        pos, cigar = aln[0], aln[1]
        return pos + sum(n for o, n in cigar if o in (0, 2, 3, 7, 8))

    def fragment(self, paired: bool):
        rng = self.rng
        loc = None
        while loc is None:
            loc = self.pick_location()
        chrom, blocks, left, right, gstrand = loc
        qlen = READ_LEN if rng.random() > 0.05 else rng.randint(35, READ_LEN - 1)
        a1 = self.align(chrom, blocks, left, qlen, **self.read_mods())
        if a1 is None:
            return []
        name = self.new_name()
        mq = self.mapq()
        extra = []
        nh = 1
        if rng.random() < 0.08:
            nh = rng.randint(2, 4)
            mq = rng.choice([0, 1, 3])
        tags = [("NH", nh), ("HI", 1)]
        base = 0
        if rng.random() < 0.02:
            base |= 0x200
        # strandedness: read1 follows the gene strand 85% of the time
        if gstrand is not None:
            r1_rev = (gstrand == "-") if rng.random() < 0.85 else (gstrand == "+")
        else:
            r1_rev = rng.random() < 0.5
        if not paired:
            flag = base | (0x10 if r1_rev else 0)
            recs = [self.seg(name, chrom, a1, flag, mq, tags=tags)]
            for h in range(2, nh + 1):
                c2 = rng.choice(["chr2", "chr1", "chr10"])
                p2 = rng.randint(0, self.chrom_len[c2] - 200)
                a2 = self.align(c2, [(0, self.chrom_len[c2])], p2, len(a1[2]))
                if a2:
                    recs.append(self.seg(name, c2, a2, 0x100 | (flag & 0x210), mq, tags=[("NH", nh), ("HI", h)]))
            if rng.random() < 0.01:  # supplementary piece
                a3 = self.align(chrom, [(0, self.chrom_len[chrom])], min(left + 3000, self.chrom_len[chrom] - 100), 20)
                if a3:
                    recs.append(self.seg(name, chrom, a3, 0x800 | (flag & 0x10), mq, tags=tags))
            return recs
        # ---- paired
        kind = rng.random()
        mate_chrom, mate_blocks, mate_start = chrom, blocks, right
        orient_ok = True
        proper = True
        if kind < 0.03:  # cross-chromosome mate
            mate_chrom = rng.choice([c for c in ("chr2", "chr1", "chr10") if c != chrom])
            mate_start = rng.randint(0, self.chrom_len[mate_chrom] - 200)
            mate_blocks = [(0, self.chrom_len[mate_chrom])]
            proper = rng.random() < 0.3  # sometimes (wrongly) flagged proper
        elif kind < 0.06:  # far apart
            mate_start = min(left + rng.randint(3000, 8000), self.chrom_len[chrom] - 200)
            mate_blocks = [(0, self.chrom_len[chrom])]
            proper = False
        elif kind < 0.09:  # same orientation
            orient_ok = False
            proper = False
        a2 = self.align(mate_chrom, mate_blocks, max(0, mate_start), qlen, **self.read_mods())
        if a2 is None:
            return []
        unm1 = rng.random() < 0.02
        unm2 = (not unm1) and rng.random() < 0.03
        if unm1 or unm2:
            proper = False
        # FR pair: the leftmost read is forward, its mate reverse; read1 is
        # the reverse (right) read when read1 should map to the minus strand.
        l_rev, r_rev = False, orient_ok
        if not orient_ok and rng.random() < 0.5:
            l_rev = r_rev = True
        left_is_r1 = not r1_rev
        segs = []
        for which, (c, aln, rev, is_r1) in enumerate(((chrom, a1, l_rev, left_is_r1),
                                                      (mate_chrom, a2, r_rev, not left_is_r1))):
            flag = base | 0x1 | (0x2 if proper else 0) | (0x10 if rev else 0) | (0x40 if is_r1 else 0x80)
            segs.append([c, aln, flag, mq])
        (c1, al1, f1, _), (c2, al2, f2, _) = segs
        # mate fields
        recs = []
        for me, mate, unm_me, unm_mate in ((segs[0], segs[1], unm1, unm2), (segs[1], segs[0], unm2, unm1)):
            c, aln, flag, q = me
            if unm_mate:
                flag |= 0x8
            if mate[2] & 0x10:
                flag |= 0x20
            if unm_me:
                flag |= 0x4
                flag &= ~0x2
                # placed at the mate's position, no CIGAR
                a = pysam.AlignedSegment(self.header_obj)
                a.query_name = name
                a.flag = flag & ~0x10
                a.reference_id = self.tid[mate[0]]
                a.reference_start = mate[1][0]
                a.mapping_quality = 0
                a.query_sequence = aln[2]
                a.query_qualities = self.qual(len(aln[2]))
                a.next_reference_id = self.tid[mate[0]]
                a.next_reference_start = mate[1][0]
                recs.append(a)
                continue
            a = self.seg(name, c, aln, flag, q, tags=tags)
            if unm_mate:
                a.next_reference_id = self.tid[c]
                a.next_reference_start = aln[0]
            else:
                a.next_reference_id = self.tid[mate[0]]
                a.next_reference_start = mate[1][0]
                if mate[0] == c:
                    lo = min(aln[0], mate[1][0])
                    hi = max(self.ref_end(aln), self.ref_end(mate[1]))
                    a.template_length = (hi - lo) if aln[0] <= mate[1][0] else -(hi - lo)
                a.set_tag("MC", "".join(f"{n}{'MIDNSHP=X'[o]}" for o, n in mate[1][1]))
            recs.append(a)
        # secondary alignments for multi-mappers (read1 only)
        for h in range(2, nh + 1):
            c3 = rng.choice(["chr2", "chr1", "chr10"])
            p3 = rng.randint(0, self.chrom_len[c3] - 200)
            a3 = self.align(c3, [(0, self.chrom_len[c3])], p3, len(a1[2]))
            if a3:
                s = self.seg(name, c3, a3, 0x100 | 0x1 | 0x40 | (base & 0x200), mq, tags=[("NH", nh), ("HI", h)])
                s.next_reference_id = self.tid[c1]
                s.next_reference_start = al1[0]
                recs.append(s)
        return recs

    def make_reads(self, paired: bool) -> list:
        rng = self.rng
        recs = []
        for _ in range(self.size):
            fr = self.fragment(paired)
            recs.extend(fr)
            if fr and rng.random() < 0.08:  # PCR duplicate copies
                marked = rng.random() < 0.6
                name = self.new_name("d")
                for r in fr:
                    if r.flag & 0x900:
                        continue
                    c = pysam.AlignedSegment.fromstring(r.to_string(), self.header_obj)
                    c.query_name = name
                    if marked:
                        c.flag |= 0x400
                    recs.append(c)
        # fully unmapped reads / pairs (no coordinates)
        for _ in range(max(2, self.size // 50)):
            name = self.new_name("u")
            n = 2 if paired else 1
            for k in range(n):
                a = pysam.AlignedSegment(self.header_obj)
                a.query_name = name
                a.flag = 0x4 | ((0x1 | 0x8 | (0x40 if k == 0 else 0x80)) if paired else 0)
                a.reference_id = -1
                a.reference_start = -1
                a.next_reference_id = -1
                a.next_reference_start = -1
                a.mapping_quality = 0
                a.query_sequence = "".join(rng.choice("ACGT") for _ in range(READ_LEN))
                a.query_qualities = self.qual(READ_LEN)
                recs.append(a)
        return recs

    @property
    def header_obj(self):
        if not hasattr(self, "_hdr"):
            self._hdr = pysam.AlignmentHeader.from_text(self.header_text())
        return self._hdr

    def header_text(self):
        return ("@HD\tVN:1.6\tSO:coordinate\n" + "".join(f"@SQ\tSN:{c}\tLN:{n}\n" for c, n in CHROMS)
                + "@PG\tID:synthetic\tPN:synthetic_data.py\tVN:1\n")

    def write_bam(self, name: str, recs: list) -> None:
        path = self.out / name

        def key(r):
            return (r.reference_id if r.reference_id >= 0 else 1 << 30, r.reference_start, r.flag, r.query_name)

        recs = sorted(recs, key=key)
        with pysam.AlignmentFile(str(path), "wb", header=self.header_obj) as fh:
            for r in recs:
                fh.write(pysam.AlignedSegment.fromstring(r.to_string(), self.header_obj))
        pysam.index(str(path))

    # --------------------------------------------------------- single cell
    def sc_reads(self) -> list:
        rng = self.rng
        cells = ["".join(rng.choice("ACGT") for _ in range(16)) for _ in range(12)]
        recs = []
        for i in range(max(40, self.size)):
            loc = None
            while loc is None:
                loc = self.pick_location()
            chrom, blocks, left, _, gstrand = loc
            if rng.random() < 0.06:
                chrom, blocks, left = "chrM", [(0, 3000)], rng.randint(0, 2900)
            kind = rng.random()
            mods = {}
            if kind < 0.1:
                mods["clip5"] = rng.randint(2, 8)
            elif kind < 0.2:
                mods["clip3"] = rng.randint(2, 8)
            aln = self.align(chrom, blocks, left, 60, p_ins=0.002, p_del=0.002, **mods)
            if aln is None:
                continue
            flag = 0x10 if rng.random() < 0.5 else 0
            if rng.random() < 0.1:
                flag |= 0x400
            cb = rng.choice(cells)
            cr = list(cb)
            if rng.random() < 0.2:
                k = rng.randrange(16)
                cr[k] = rng.choice([x for x in "ACGTN" if x != cr[k]])
            umi = "".join(rng.choice("ACGT") for _ in range(10))
            ur = list(umi)
            if rng.random() < 0.15:
                k = rng.randrange(10)
                ur[k] = rng.choice([x for x in "ACGT" if x != ur[k]])
            tags = [("NH", 1)]
            if rng.random() < 0.9:
                tags.append(("CR", "".join(cr)))
            if rng.random() < 0.85:
                tags.append(("CB", cb + "-1"))
            if rng.random() < 0.95:
                tags.append(("UR", "".join(ur)))
            if rng.random() < 0.85:
                tags.append(("UB", umi))
            xf = rng.choice([0, 1, 17, 25, 8, 1, 17])
            tags.append(("xf", xf))
            tags.append(("RE", rng.choice("EEEINI")))  # upstream KeyErrors on confident reads lacking RE
            r = rng.random()
            if r < 0.6:
                tags.append(("TX", "T1,+100,60M"))
                tags.append(("GN", "G1"))
            elif r < 0.8:
                tags.append(("AN", "T2,-50,60M"))
            name = self.new_name("sc")
            recs.append(self.seg(name, chrom, aln, flag, rng.choice([255, 255, 3, 1, 0]), tags=tags))
            if rng.random() < 0.05:  # same read name aligned twice (secondary)
                recs.append(self.seg(name, chrom, aln, flag | 0x100, 1, tags=tags))
        return recs

    # --------------------------------------------------------- saturation
    def saturation_bams(self) -> None:
        rng = self.rng
        tx = [t for t in self.transcripts if t["chrom"] in self.genome and len(t["exons"]) >= 2]
        t = tx[0]
        (s0, e0), (s1, e1) = t["exons"][0], t["exons"][1]
        noise = []

        def noisy_flags():
            return rng.choice([0x200, 0x400, 0x100])

        # junction: every qualifying spliced read has the same intron (e0, s1)
        recs = []
        for i in range(max(3, self.size // 20)):
            x = rng.randint(10, 40)
            start = e0 - x
            if start < s0:
                start = s0
            aln = self.align(t["chrom"], [(s0, e0), (s1, e1)], start, READ_LEN, p_ins=0, p_del=0)
            recs.append(self.seg(self.new_name("js"), t["chrom"], aln, rng.choice([0, 0x10]), 255))
        for i in range(max(5, self.size // 5)):
            loc = None
            while loc is None:
                loc = self.pick_location()
            chrom, blocks, left, _, _ = loc
            aln = self.align(chrom, blocks, left, READ_LEN)
            if aln is None:
                continue
            spliced = any(o == 3 for o, _ in aln[1])
            if spliced:
                # only non-qualifying spliced reads, or ones on a chromosome absent from the BED
                q = rng.random()
                if q < 0.5:
                    recs.append(self.seg(self.new_name("jn"), chrom, aln, noisy_flags(), 255))
                else:
                    recs.append(self.seg(self.new_name("jn"), chrom, aln, 0, rng.choice([0, 10, 29])))
            else:
                recs.append(self.seg(self.new_name("jn"), chrom, aln, rng.choice([0, 0x10]), rng.choice([255, 0])))
        self.write_bam("sat_junction.bam", recs)

        # RPKM: every qualifying read is the same unspliced alignment
        recs = []
        p = s0 + 5 if e0 - s0 > READ_LEN + 10 else s0
        blocks = [(0, self.chrom_len[t["chrom"]])]
        for i in range(max(3, self.size // 10)):
            aln = self.align(t["chrom"], blocks, p, READ_LEN, p_ins=0, p_del=0)
            flag = 0x1 | rng.choice([0x40, 0x80]) | rng.choice([0, 0x10]) | 0x8
            recs.append(self.seg(self.new_name("rs"), t["chrom"], aln, flag, 255))
        for i in range(max(5, self.size // 5)):
            loc = None
            while loc is None:
                loc = self.pick_location()
            chrom, blocks2, left, _, _ = loc
            aln = self.align(chrom, blocks2, left, READ_LEN)
            if aln is None:
                continue
            flag = 0x1 | 0x40 | 0x8
            if rng.random() < 0.5:
                recs.append(self.seg(self.new_name("rn"), chrom, aln, flag | (noisy_flags()), 255))
            else:
                recs.append(self.seg(self.new_name("rn"), chrom, aln, flag, rng.choice([0, 10, 29])))
        self.write_bam("sat_rpkm.bam", recs)

    # ------------------------------------------------------------ fastx
    def write_fastx(self, recs: list) -> None:
        rng = self.rng
        seen = set()
        fq_lines, fa_lines = [], []
        for r in recs:
            if r.flag & 0x900 or r.query_name in seen or not r.query_sequence:
                continue
            seen.add(r.query_name)
            s = r.query_sequence
            q = "".join(chr(33 + x) for x in r.query_qualities)
            if len(s) != READ_LEN:
                continue
            fq_lines.append(f"@{r.query_name}\n{s}\n+\n{q}\n")
            fa_lines.append(f">{r.query_name}\n{s}\n")
        (self.out / "reads.fq").write_text("".join(fq_lines))
        (self.out / "reads.fa").write_text("".join(fa_lines))
        with gzip.GzipFile(self.out / "reads.fq.gz", "wb", mtime=0) as fh:
            fh.write("".join(fq_lines).encode())
        (self.out / "reads.fq.bz2").write_bytes(bz2.compress("".join(fq_lines).encode()))
        # fixed-length barcode+UMI reads for the sc_* sequence tools
        cells = ["".join(rng.choice("ACGT") for _ in range(16)) for _ in range(8)]
        bfq, bfa = [], []
        for i in range(max(30, self.size)):
            s = list(rng.choice(cells) + "".join(rng.choice("ACGT") for _ in range(12)))
            if rng.random() < 0.05:
                s[rng.randrange(len(s))] = "N"
            s = "".join(s)
            q = "".join(chr(33 + rng.choice([2, 11, 25, 32, 37, 37, 41])) for _ in s)
            bfq.append(f"@bc{i}\n{s}\n+\n{q}\n")
            bfa.append(f">bc{i}\n{s}\n")
        (self.out / "barcodes.fq").write_text("".join(bfq))
        (self.out / "barcodes.fa").write_text("".join(bfa))
        with gzip.GzipFile(self.out / "barcodes.fq.gz", "wb", mtime=0) as fh:
            fh.write("".join(bfq).encode())

    # ------------------------------------------------------------ bigwig
    def write_bigwigs(self) -> None:
        rng = self.rng
        for name, neg in (("sig1.bw", 0.1), ("sig2.bw", 0.4)):
            bw = pyBigWig.open(str(self.out / name), "w")
            chroms = [(c, n) for c, n in CHROMS if not (name == "sig2.bw" and c == "chr10")]
            bw.addHeader(chroms)
            for c, n in chroms:
                if c == "chrX":
                    continue
                pos = rng.randint(0, 50)
                mode = 0
                while pos < n - 200:
                    mode = rng.choice([0, 0, 1, 2])
                    if mode == 0:  # bedGraph run
                        starts, ends, vals = [], [], []
                        for _ in range(rng.randint(3, 12)):
                            ln = rng.randint(1, 60)
                            if pos + ln >= n:
                                break
                            v = round(rng.uniform(0.5, 30), rng.choice([0, 1, 2]))
                            if rng.random() < neg:
                                v = -v
                            starts.append(pos)
                            ends.append(pos + ln)
                            vals.append(float(v))
                            pos += ln + rng.choice([0, 0, rng.randint(1, 40)])
                        if starts:
                            bw.addEntries([c] * len(starts), starts, ends=ends, values=vals)
                    elif mode == 1:  # variableStep
                        span = rng.randint(1, 10)
                        starts, vals = [], []
                        for _ in range(rng.randint(3, 10)):
                            if pos + span >= n:
                                break
                            starts.append(pos)
                            v = float(round(rng.uniform(0.25, 12), 2))
                            vals.append(-v if rng.random() < neg else v)
                            pos += span + rng.randint(0, 20)
                        if starts:
                            bw.addEntries(c, starts, values=vals, span=span)
                    else:  # fixedStep
                        span = rng.randint(1, 8)
                        step = span + rng.randint(0, 5)
                        k = rng.randint(3, 15)
                        if pos + step * k >= n:
                            break
                        vals = [float(rng.randint(-3 if neg > 0.2 else 0, 20)) for _ in range(k)]
                        bw.addEntries(c, pos, values=vals, span=span, step=step)
                        pos += step * k
                    pos += rng.randint(0, 400)
            bw.close()


def make_depth_cap_fixture(out: Path) -> None:
    """A pileup deep enough to bind pysam's default max_depth of 8000.

    Written because geneBody_coverage's real-data divergence (DIV-0024) lives in
    max_depth semantics, and the main synthetic fixture has depth around 40 -- far
    below any cap, so the differential suite passed 90/90 while the command was wrong
    on real data. A fixture that cannot reach the threshold cannot test the
    threshold.

    8,100 copies of one 20M20D20M read, so the cap binds at the 20M, the reads are
    in a deletion (is_del) through the middle, and return matched at the second 20M.
    A handful of clean single reads are included so a port that reports 0 at a
    position that clean reads demonstrably cover cannot pass.
    """
    contig = "chrCap"
    start, length = 20_000, 100
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": contig, "LN": 40_000}]}

    def segment(pos, name, cigar, seq, qual="I" * 40):
        a = pysam.AlignedSegment()
        a.query_name = name
        a.query_sequence = seq
        a.query_qualities = pysam.qualitystring_to_array(qual)
        a.flag = 0
        a.reference_id = 0
        a.reference_start = pos
        a.cigarstring = cigar
        a.mapping_quality = 60
        a.next_reference_id = 0
        a.next_reference_start = pos
        a.template_length = 0
        return a

    records = []
    for i in range(8_100):
        records.append(segment(start + 10, f"cap{i}", "20M20D20M", "A" * 40))
    for i in range(40):
        records.append(segment(start + 20 + i, f"clean{i}", "20M", "C" * 20, "I" * 20))
    records.sort(key=lambda r: r.reference_start)

    bam = out / "depth_cap.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as handle:
        for record in records:
            handle.write(record)
    pysam.index(str(bam))
    (out / "depth_cap.bed12").write_text(
        f"{contig}\t{start}\t{start + length}\tCAP1\t0\t+\t{start}\t"
        f"{start + length}\t0\t1\t{length},\t0,\n"
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--size", type=int, default=500, help="number of fragments per BAM")
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    g = Gen(args.seed, args.size, out)
    g.make_model()
    g.write_model()
    g.write_sequences()
    pe = g.make_reads(paired=True)
    g.write_bam("pe.bam", pe)
    se = g.make_reads(paired=False)
    g.write_bam("se.bam", se)
    # the same reads minus the coordinate-less unmapped records (FPKM_count.py
    # crashes on those unless -u is given: DIV-0023)
    g.write_bam("pe_placed.bam", [r for r in pe if r.reference_id >= 0])
    g.write_bam("se_placed.bam", [r for r in se if r.reference_id >= 0])
    g.write_bam("sc.bam", g.sc_reads())
    g.saturation_bams()
    g.write_fastx(pe)
    g.write_bigwigs()
    make_depth_cap_fixture(out)
    # a BAM-list file for geneBody_coverage.py/tin.py
    (out / "bams.txt").write_text(f"{out / 'pe.bam'}\n{out / 'se.bam'}\n")


if __name__ == "__main__":
    main()
