#!/usr/bin/env bash
# Align one ENA run with the pinned STAR index and produce a coordinate-sorted,
# indexed BAM.
#
# Outputs land in datasets/aligned/<run>/:
#   <run>.Aligned.sortedByCoord.out.bam   what STAR wrote
#   <run>.bam                              the deliverable (sorted, indexed)
#   <run>.bam.bai
#   <run>.Log.final.out                    STAR's own report: mapping rates,
#                                          input read count, splice counts
#   <run>.align.json                       machine-readable record of what was run
#
# `Log.final.out` is kept deliberately. It is the independent record of how
# many reads STAR actually placed, which is what makes a downstream "this BAM
# has 12M reads" claim checkable rather than assumed.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF_DIR="${REF_DIR:-$HERE/reference}"
IDX_DIR="${IDX_DIR:-$HERE/star_index}"
RAW_DIR="${RAW_DIR:-$HERE/raw}"
OUT_ROOT="${OUT_ROOT:-$HERE/aligned}"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-/scratch/mdra00001/tmp/mamba}"
export MAMBA_ROOT_PREFIX
MM="${MICROMAMBA:-/scratch/mdra00001/tmp/bin/micromamba}"
ENV_NAME="${ENV_NAME:-t4star}"
THREADS="${THREADS:-$(nproc)}"

[[ $# -eq 1 ]] || { echo "usage: $0 <run_accession>" >&2; exit 2; }
RUN="$1"
OUT="$OUT_ROOT/$RUN"
mkdir -p "$OUT"

[[ -f "$IDX_DIR/SA" ]] || { echo "error: no STAR index in $IDX_DIR (run datasets/build_star_index.sh)" >&2; exit 1; }

# ENA names the two mates <run>_1.fastq.gz / _2.fastq.gz. Missing either is a
# hard error rather than a silent single-end alignment, which would quietly
# change what every downstream command sees.
R1="$RAW_DIR/${RUN}_1.fastq.gz"
R2="$RAW_DIR/${RUN}_2.fastq.gz"
for f in "$R1" "$R2"; do
  [[ -f "$f" ]] || { echo "error: missing $f (run datasets/fetch_fastq.py $RUN)" >&2; exit 1; }
done

if [[ -f "$OUT/${RUN}.bam.bai" ]]; then
  echo "ok       $RUN already aligned"
  exit 0
fi

STAR="$MM run -n $ENV_NAME STAR"
SAMBAM="$MM run -n $ENV_NAME samtools"

echo "align    $RUN  ($THREADS threads)"
$STAR --genomeDir "$IDX_DIR" \
  --readFilesIn "$R1" "$R2" \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --runThreadN "$THREADS" \
  --outFileNamePrefix "$OUT/${RUN}." \
  --outTmpDir "$OUT/tmp" \
  > "$OUT/${RUN}.align.stdout" 2>&1

# STAR writes SortedByCoordinate itself when asked, but samtools sort + index
# is kept as an explicit step so the BAM is provably well-formed and
# coordinate-ordered rather than trusted from a flag. `sort -c` then verifies it.
$SAMBAM sort -@ "$THREADS" -o "$OUT/${RUN}.bam" "$OUT/${RUN}.Aligned.sortedByCoord.out.bam"
$SAMBAM index -@ "$THREADS" "$OUT/${RUN}.bam"
$SAMBAM quickcheck -v "$OUT/${RUN}.bam"

# Record what was actually run, so an aligned BAM can be traced back to inputs
# and parameters without re-deriving them.
{
  echo "{"
  echo "  \"run\": \"$RUN\","
  echo "  \"star_version\": \"$($STAR --version 2>/dev/null | tail -1)\","
  echo "  \"threads\": $THREADS,"
  echo "  \"read1\": \"$(basename "$R1")\","
  echo "  \"read2\": \"$(basename "$R2")\","
  echo "  \"read1_md5\": \"$(md5sum "$R1" | awk '{print $1}')\","
  echo "  \"read2_md5\": \"$(md5sum "$R2" | awk '{print $1}')\","
  echo "  \"index_stamp\": \"$(tr '\n' ' ' < "$IDX_DIR/.built-"* 2>/dev/null | tr -s ' ')\","
  echo "  \"input_read_pairs\": $($SAMBAM view -c -f 1 "$OUT/${RUN}.bam" 2>/dev/null || echo null)"
  echo "}"
} > "$OUT/${RUN}.align.json"

rm -rf "$OUT/tmp"
echo "ok       $RUN -> $OUT/${RUN}.bam"
grep -E "input reads|uniquely mapped reads %|Number of splices" "$OUT/${RUN}.Log.final.out" 2>/dev/null || true
