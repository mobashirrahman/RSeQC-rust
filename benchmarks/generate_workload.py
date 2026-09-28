#!/usr/bin/env python3
"""Generate deterministic synthetic workloads for benchmarking.

This creates coordinate-sorted, indexed BAM files with paired-end reads,
a matching BED12 gene model, and FASTQ files. All output is deterministic
and controlled by a CLI argument specifying read count.

Usage:
    oracle/venv/bin/python3 benchmarks/generate_workload.py --size 1000 --output-dir workloads/test_1000

Output:
    workloads/test_1000/
        manifest.json          - Seed, version, parameters
        reads.bam              - Coordinate-sorted BAM with paired-end reads
        reads.bam.bai          - BAM index
        model.bed12            - BED12 gene model
        reads_1.fastq          - Paired-end read 1 FASTQ
        reads_2.fastq          - Paired-end read 2 FASTQ
"""
import argparse
import json
import os
import random
import sys
from pathlib import Path
from typing import NamedTuple

# Require pysam for BAM/index generation
try:
    import pysam
except ImportError:
    print("Error: pysam is required. Install with: pip install pysam", file=sys.stderr)
    sys.exit(1)

# Generator version for reproducibility tracking
GENERATOR_VERSION = "1.0.0"
FIXED_SEED = 42


class ReadInfo(NamedTuple):
    qname: str
    seq: str
    qual: str
    pos: int
    mapq: int
    is_reverse: bool
    mate_pos: int


def generate_reads(num_reads: int, seed: int) -> tuple[list[ReadInfo], dict]:
    """Generate deterministic paired-end reads with mix of spliced/unspliced."""
    rng = random.Random(seed)

    # Synthetic genome: one chromosome, 10000 bp
    chrom_len = 10000

    # Gene model: 3 genes with introns
    # Gene 1: [100-500] with intron [200-300] (exons: 100-200, 300-500)
    # Gene 2: [2000-3000] with intron [2400-2600] (exons: 2000-2400, 2600-3000)
    # Gene 3: [5000-7000] with intron [5800-6200] (exons: 5000-5800, 6200-7000)
    genes = [
        {"name": "tx1", "chrom": "chr1", "start": 100, "end": 500, "strand": "+",
         "exons": [(100, 200), (300, 500)]},
        {"name": "tx2", "chrom": "chr1", "start": 2000, "end": 3000, "strand": "-",
         "exons": [(2000, 2400), (2600, 3000)]},
        {"name": "tx3", "chrom": "chr1", "start": 5000, "end": 7000, "strand": "+",
         "exons": [(5000, 5800), (6200, 7000)]},
    ]

    reads = []

    # Generate paired-end reads
    for i in range(num_reads):
        # Randomly select a gene
        gene = rng.choice(genes)
        strand = rng.choice([True, False])  # True = reverse complement

        # Pick position in gene
        gene_start = gene["start"]
        gene_end = gene["end"]
        pos1 = rng.randint(gene_start, gene_end - 101)  # Leave room for read pair
        pos2 = pos1 + 100  # 100bp fragment size

        # Generate sequence and quality (all A's with Phred 40)
        seq = "A" * 100
        qual = "I" * 100  # Phred 40

        qname = f"read_{i:06d}"

        # Adjust for reverse strand
        if strand:
            seq_rev = reverse_complement(seq)
        else:
            seq_rev = seq

        # Read 1: first in pair
        read1 = ReadInfo(
            qname=qname,
            seq=seq if not strand else seq_rev,
            qual=qual,
            pos=pos1,
            mapq=40,
            is_reverse=strand,
            mate_pos=pos2,
        )

        # Read 2: second in pair (reverse complement of read 1's mate)
        read2 = ReadInfo(
            qname=qname,
            seq=seq_rev if not strand else seq,
            qual=qual,
            pos=pos2,
            mapq=40,
            is_reverse=not strand,
            mate_pos=pos1,
        )

        reads.append((read1, read2, gene))

    metadata = {
        "seed": seed,
        "num_reads": num_reads,
        "fragment_size": 100,
        "read_length": 100,
        "chrom": "chr1",
        "chrom_length": chrom_len,
        "genes": genes,
    }

    return reads, metadata


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement.get(c, "N") for c in reversed(seq))


def write_bam(reads_data: list, output_path: Path, metadata: dict):
    """Write BAM file with coordinate-sorted paired-end reads."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Create header
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": metadata["chrom_length"]}],
    }

    # Collect all individual reads for sorting
    all_reads = []
    for read1, read2, gene in reads_data:
        all_reads.append((read1, read2, gene, True))   # True = is_first_in_pair
        all_reads.append((read1, read2, gene, False))  # False = is_second_in_pair

    # Sort by position of each read
    def get_pos(item):
        read1, read2, gene, is_first = item
        return read1.pos if is_first else read2.pos

    sorted_reads = sorted(all_reads, key=get_pos)

    with pysam.AlignmentFile(str(output_path), "wb", header=header) as bam:
        for read1, read2, gene, is_first in sorted_reads:
            if is_first:
                # Write read 1
                a = pysam.AlignedSegment()
                a.query_name = read1.qname
                a.query_sequence = read1.seq
                a.query_qualities = pysam.qualitystring_to_array(read1.qual)
                a.reference_id = 0  # chr1
                a.reference_start = read1.pos
                a.mapping_quality = read1.mapq
                a.is_read1 = True
                a.is_read2 = False
                a.is_reverse = read1.is_reverse
                a.next_reference_id = 0
                a.next_reference_start = read1.mate_pos
                a.template_length = read2.pos - read1.pos if not read1.is_reverse else -(read1.pos - read2.pos)
                a.cigar = [(0, 100)]  # 100M
                bam.write(a)
            else:
                # Write read 2
                a = pysam.AlignedSegment()
                a.query_name = read2.qname
                a.query_sequence = read2.seq
                a.query_qualities = pysam.qualitystring_to_array(read2.qual)
                a.reference_id = 0  # chr1
                a.reference_start = read2.pos
                a.mapping_quality = read2.mapq
                a.is_read1 = False
                a.is_read2 = True
                a.is_reverse = read2.is_reverse
                a.next_reference_id = 0
                a.next_reference_start = read2.mate_pos
                a.template_length = -(read2.pos - read1.pos) if not read1.is_reverse else (read1.pos - read2.pos)
                a.cigar = [(0, 100)]  # 100M
                bam.write(a)

    # Create index
    pysam.index(str(output_path))


def write_bed12(genes: list, output_path: Path):
    """Write BED12 gene model file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        for gene in genes:
            exons = gene["exons"]
            chrom = gene["chrom"]
            name = gene["name"]
            score = 0
            strand = gene["strand"]
            cds_start = exons[0][0]
            cds_end = exons[-1][1]
            rgb = "0,0,0"

            exon_count = len(exons)
            exon_sizes = ",".join(str(e[1] - e[0]) for e in exons) + ","
            exon_starts = ",".join(str(e[0] - gene["start"]) for e in exons) + ","

            fields = [
                chrom,
                gene["start"],
                gene["end"],
                name,
                score,
                strand,
                cds_start,
                cds_end,
                rgb,
                exon_count,
                exon_sizes,
                exon_starts,
            ]
            f.write("\t".join(str(f) for f in fields) + "\n")


def write_fastq(reads_data: list, read_num: int, output_path: Path):
    """Write FASTQ file for reads (read_num = 1 or 2)."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        for read1, read2, _ in reads_data:
            read = read1 if read_num == 1 else read2
            # FASTQ format: @name, sequence, +, qualities
            f.write(f"@{read.qname}/1\n")
            f.write(f"{read.seq}\n")
            f.write("+\n")
            f.write(f"{read.qual}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate deterministic synthetic workloads for benchmarking"
    )
    parser.add_argument(
        "--size",
        type=int,
        default=1000,
        help="Number of paired-end reads to generate (default: 1000)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for workload files",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=FIXED_SEED,
        help=f"Random seed for reproducibility (default: {FIXED_SEED})",
    )

    args = parser.parse_args()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate reads
    print(f"Generating {args.size} paired-end reads...", file=sys.stderr)
    reads_data, metadata = generate_reads(args.size, args.seed)

    # Write BAM
    bam_path = output_dir / "reads.bam"
    print(f"Writing BAM to {bam_path}...", file=sys.stderr)
    write_bam(reads_data, bam_path, metadata)

    # Write BED12
    bed_path = output_dir / "model.bed12"
    print(f"Writing BED12 to {bed_path}...", file=sys.stderr)
    write_bed12(metadata["genes"], bed_path)

    # Write FASTQ files
    fastq1_path = output_dir / "reads_1.fastq"
    fastq2_path = output_dir / "reads_2.fastq"
    print(f"Writing FASTQ to {fastq1_path} and {fastq2_path}...", file=sys.stderr)
    write_fastq(reads_data, 1, fastq1_path)
    write_fastq(reads_data, 2, fastq2_path)

    # Write manifest
    manifest = {
        "generator_version": GENERATOR_VERSION,
        "seed": args.seed,
        "num_reads": args.size,
        "files": {
            "bam": "reads.bam",
            "bam_index": "reads.bam.bai",
            "bed12": "model.bed12",
            "fastq_1": "reads_1.fastq",
            "fastq_2": "reads_2.fastq",
        },
        "metadata": metadata,
    }

    manifest_path = output_dir / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"Workload generated successfully in {output_dir}", file=sys.stderr)
    print(f"Manifest: {manifest_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
