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
IDX_ROOT="${IDX_ROOT:-$HERE/star_index}"
# Same per-assembly directory naming as datasets/build_star_index.sh, so a build and
# an alignment agree on which index is meant without either having to be told. The
# human and rat panels differ in genome, annotation, sjdbOverhang and contig set, and
# a mismatched overhang degrades splice-junction detection silently -- which is
# exactly what endpoint E5 measures.
GENOME="${GENOME:-$REF_DIR/hg38.fa.gz}"
GTF="${GTF:-$REF_DIR/gencode.v47.annotation.gtf.gz}"
CONTIGS="${CONTIGS:-chr1 chr17 chrM}"
SJDB_OVERHANG="${SJDB_OVERHANG:-149}"
INDEX_KEY="$(basename "$GENOME" .gz)-$(basename "$GTF" .gz)-o$SJDB_OVERHANG-$(echo "$CONTIGS" | tr ' ' '_')"
IDX_DIR="${IDX_DIR:-$IDX_ROOT/$INDEX_KEY}"
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

# A STAR index is only valid for the genome, annotation, overhang and contig set it
# was built from, and its parts are overwritten in place. So the check is not just
# that SA exists but that the index directory carries a stamp naming THIS
# annotation's digest: otherwise an index left over from a different build would be
# used silently and the alignment would record a provenance nobody can check.
if [[ ! -f "$IDX_DIR/SA" ]]; then
  echo "error: no STAR index in $IDX_DIR" >&2
  echo "       build it with datasets/build_star_index.sh using:" >&2
  echo "         GENOME=$GENOME GTF=$GTF CONTIGS='$CONTIGS' SJDB_OVERHANG=$SJDB_OVERHANG" >&2
  echo "       or point IDX_DIR at the right one under $IDX_ROOT:" >&2
  ls -1 "$IDX_ROOT" 2>/dev/null | sed 's/^/         /' >&2 || true
  exit 1
fi
STAMP_FILE="$(ls "$IDX_DIR"/.built-* 2>/dev/null | head -1 || true)"
if [[ -z "$STAMP_FILE" ]]; then
  echo "error: $IDX_DIR has no .built-* stamp; cannot verify which annotation this" >&2
  echo "       index was built from. Refusing to align against an unrecorded index." >&2
  exit 1
fi
WANT_GTF_SHA="$(sha256sum "$IDX_DIR/annotation.gtf" | awk '{print $1}')"
HAVE_GTF_SHA="$(grep -o 'gtf_sha256=[0-9a-f]*' "$STAMP_FILE" | head -1 | cut -d= -f2)"
if [[ -n "$HAVE_GTF_SHA" ]] && [[ "$WANT_GTF_SHA" != "$HAVE_GTF_SHA" ]]; then
  echo "error: $IDX_DIR/annotation.gtf does not match the digest its stamp records." >&2
  echo "       stamp: gtf_sha256=$HAVE_GTF_SHA" >&2
  echo "       file : $WANT_GTF_SHA" >&2
  echo "       The index on disk is not the one the stamp describes. Rebuild it." >&2
  exit 1
fi

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
  # The whole stamp, including gtf_sha256, not just its parameters. This is what makes
  # 'was this BAM aligned against that index?' answerable later: the earlier stamp
  # recorded genome/contigs/overhang but no annotation digest, so a re-alignment
  # against a rebuilt index was undetectable from the alignment's own record.
  echo "  \"index_stamp\": \"$(tr '\n' ' ' < "$STAMP_FILE" | tr -s ' ')\","
  echo "  \"index_dir\": \"$IDX_DIR\","
  echo "  \"input_read_pairs\": $($SAMBAM view -c -f 1 "$OUT/${RUN}.bam" 2>/dev/null || echo null)"
  echo "}"
} > "$OUT/${RUN}.align.json"

rm -rf "$OUT/tmp"
echo "ok       $RUN -> $OUT/${RUN}.bam"
grep -E "input reads|uniquely mapped reads %|Number of splices" "$OUT/${RUN}.Log.final.out" 2>/dev/null || true
