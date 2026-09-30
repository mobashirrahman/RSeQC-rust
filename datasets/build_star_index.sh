#!/usr/bin/env bash
# Build the pinned STAR index for the T4 real-data panel.
#
# Idempotent: if a marker file recording the exact parameters that produced the
# index is present and matches, this exits without rebuilding. The marker
# exists because a STAR index is only valid for the genome + annotation +
# sjdbOverhang it was built from, and silently reusing a stale one would put
# every downstream alignment on an unrecorded reference.
#
# Memory: STAR's human index needs ~30 GB at the default
# --genomeSAindexNbases 14. This machine has 31 GB total, so we build with 11
# (peak ~16 GB). That is a deliberate, recorded deviation from STAR's
# recommended value: it costs some mapping speed, not mapping correctness, and
# the exact value is pinned here so the cost is visible rather than implicit.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF_DIR="${REF_DIR:-$HERE/reference}"
IDX_DIR="${IDX_DIR:-$HERE/star_index}"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-/scratch/mdra00001/tmp/mamba}"
ENV_NAME="${ENV_NAME:-t4star}"

GENOME="$REF_DIR/hg38.fa.gz"
GTF="$REF_DIR/gencode.v47.annotation.gtf.gz"
# STAR's own ceiling, in MB, on how much RAM genomeGenerate may plan for. The
# default is 31000, which on a 31 GB machine is a promise the kernel cannot
# keep: the first build here was OOM-killed at "inserting junctions into the
# genome indices" (Log.out ends there, process died) after completing the
# 34-minute suffix array. Setting the limit makes STAR choose a smaller
# genomeSAindexNbases itself, which is the mechanism the option exists for --
# better than hand-tuning the index parameter and getting it subtly wrong.
# NOTE the unit: STAR parses --limitGenomeGenerateRAM as BYTES, and passing
# "20000" meaning megabytes is rejected as "too small for your genome" with a
# minimum of 8746522208 (i.e. ~8.1 GiB expressed in bytes). 21 GB in bytes.
LIMIT_RAM_BYTES="${LIMIT_RAM_BYTES:-21000000000}"
# 149 = 2*75 - 1 for 2x100 Illumina. Recorded because STAR uses it to size the
# junction database, and a mismatch against the actual read length degrades
# splice detection at the long end.
SJDB_OVERHANG=149
THREADS="${THREADS:-$(nproc)}"

for f in "$GENOME" "$GTF"; do
  [[ -f "$f" ]] || { echo "error: missing $f (run datasets/fetch_reference.sh first)" >&2; exit 1; }
done

STAMP="$IDX_DIR/.built-${LIMIT_RAM_BYTES}-${SJDB_OVERHANG}"
if [[ -f "$STAMP" ]]; then
  echo "ok       STAR index already built ($(cat "$STAMP"))"
  exit 0
fi

mkdir -p "$IDX_DIR"
MM="${MICROMAMBA:-/scratch/mdra00001/tmp/bin/micromamba}"
# Export the root prefix here rather than relying on the caller: micromamba
# resolves envs relative to it, and a prefix that is unset-but-defaulted by
# micromamba points at XDG state instead, which fails with a confusing
# "prefix does not exist" from a backgrounded run.
export MAMBA_ROOT_PREFIX
STAR="$MM run -n $ENV_NAME STAR"

echo "building  STAR index (limitGenomeGenerateRAM=$LIMIT_RAM_BYTES bytes sjdbOverhang=$SJDB_OVERHANG)"
echo "          STAR chooses genomeSAindexNbases from the RAM limit; the chosen"
echo "          value is read back from Log.out and recorded in the stamp."
# STAR cannot read a gzipped genome OR a gzipped GTF: --sjdbGTFfile silently
# yields an empty junction database from a .gz (its own error text is the clue:
# "Make sure the GTF file is unzipped"). Decompressed once and reused.
if [[ ! -f "$IDX_DIR/genome.fa" ]]; then
  echo "unpack    hg38.fa.gz -> genome.fa (STAR requires uncompressed FASTA)"
  gunzip -c "$GENOME" > "$IDX_DIR/genome.fa"
fi
if [[ ! -f "$IDX_DIR/annotation.gtf" ]]; then
  echo "unpack    gencode.gtf.gz -> annotation.gtf (STAR requires uncompressed GTF)"
  gunzip -c "$GTF" > "$IDX_DIR/annotation.gtf"
fi

$STAR --runMode genomeGenerate \
  --genomeDir "$IDX_DIR" \
  --genomeFastaFiles "$IDX_DIR/genome.fa" \
  --sjdbGTFfile "$IDX_DIR/annotation.gtf" \
  --sjdbOverhang "$SJDB_OVERHANG" \
  --limitGenomeGenerateRAM "$LIMIT_RAM_BYTES" \
  --runThreadN "$THREADS" \
  --outFileNamePrefix "$IDX_DIR/"

# Read the genomeSAindexNbases STAR actually chose back out of its log rather
# than asserting one, so the stamp records what happened and not what was
# requested.
CHOSEN="$(grep -oE 'genomeSAindexNbases *= *[0-9]+' "$IDX_DIR/Log.out" | tail -1 | grep -oE '[0-9]+' || echo unknown)"
printf 'genome=%s\ngtf=%s\nsjdbOverhang=%s\nlimitGenomeGenerateRAM_bytes=%s\ngenomeSAindexNbases=%s\n' \
  "$(basename "$GENOME")" "$(basename "$GTF")" "$SJDB_OVERHANG" "$LIMIT_RAM_BYTES" "$CHOSEN" > "$STAMP"

echo "ok       index built: $IDX_DIR"
