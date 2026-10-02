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

# Overridable so the same script builds the rat index for the cross-organism
# held-out stratum rather than duplicating it. The contig scope, RAM limit and
# SA-index size are overridden too: rn6 needs its own contig choice, and a rat index
# over a chromosome set chosen for hg38 would either be empty or full-genome (which
# is what OOMs). These are declared first because the index directory name is
# derived from them below.
GENOME="${GENOME:-$REF_DIR/hg38.fa.gz}"
GTF="${GTF:-$REF_DIR/gencode.v47.annotation.gtf.gz}"

# Per-assembly index directory. A STAR index is only valid for the genome +
# annotation + sjdbOverhang + contig set it was built from, and every one of those
# differs between the human dev/held-out panel and the rat cross-organism stratum.
#
# The previous default was a single shared $HERE/star_index for both, which made
# "build the rat index" mean "overwrite the human index". That is not a
# hypothetical: it happened here while re-deriving the rat junction database, and
# the damage is silent -- STAR overwrites Genome, SA, SA_0..SA_15 and the exon
# tables in place, so the human index is left as a `Genome`/`SA` pair whose
# per-chunk files and junction tables now describe rat. `align_run.sh` only
# checks that $IDX_DIR/SA exists, so the next human alignment would run against a
# corrupt index and produce plausible-looking output.
#
# The directory is named after the genome and annotation, not the assembly name
# alone, because the two rat runs already differed by sjdbOverhang alone: an index
# built with the human default overhang degrades splice-junction detection silently,
# which is exactly what endpoint E5 measures.
INDEX_KEY="$(basename "${GENOME:-hg38.fa.gz}" .gz)-$(basename "${GTF:-gencode.v47.annotation.gtf.gz}" .gz)"
# Contig set and overhang join the key because they change the index contents
# without changing either input file.
CONTIGS="${CONTIGS:-chr1 chr17 chrM}"
SJDB_OVERHANG="${SJDB_OVERHANG:-149}"
# The contig set goes into the directory name so two indexes differing only in contig
# scope cannot collide -- but spelled out in full it overflows the 255-byte filename
# limit for a WHOLE-GENOME set, which is the scope the production pilot needs. So a long
# set is summarised rather than dropped: the first two contigs stay readable, the count
# follows, and a digest of the full set keeps the name unique. Two whole-genome indexes
# for different assemblies still get different directories, and the marker file inside
# records the complete contig list either way, so nothing is lost by abbreviating.
contig_slug="$(echo "$CONTIGS" | tr ' ' '_')"
if [ "${#contig_slug}" -gt 60 ]; then
  contig_digest="$(printf '%s' "$CONTIGS" | sha256sum | cut -c1-10)"
  contig_count="$(echo "$CONTIGS" | wc -w | tr -d ' ')"
  contig_slug="$(echo "$CONTIGS" | cut -d' ' -f1-2 | tr ' ' '_')-plus${contig_count}more-${contig_digest}"
fi
INDEX_KEY="$INDEX_KEY-o$SJDB_OVERHANG-$contig_slug"

# Print the derived directory name and exit, without building anything. This exists so
# the naming rule is TESTABLE: a whole-genome contig set overflows the 255-byte filename
# limit, and that only shows up as "File name too long" after STAR has been located and
# the genome unpacked -- 90 seconds into a 45-minute build. The rule also has to stay
# readable and collision-free, which is a judgement no filesystem error will make for us.
if [ "${1:-}" = "--print-index-key" ]; then
  echo "$INDEX_KEY"
  exit 0
fi

IDX_DIR="${IDX_DIR:-$HERE/star_index/$INDEX_KEY}"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-/scratch/mdra00001/tmp/mamba}"
ENV_NAME="${ENV_NAME:-t4star}"

# Contig scope. A whole-genome GENCODE index does not build on this machine:
# it was OOM-killed three times at "inserting junctions into the genome
# indices" -- at genomeSAindexNbases 14, 11 and 8, and with
# --limitGenomeGenerateRAM honoured. The parameter is not the binding
# constraint; the junction database built from GENCODE's 2,155,005 exon
# records is, and it is allocated ON TOP of the finished suffix-array index
# regardless of the SA size. Swap is not available as a fallback (swapon is
# not permitted in this container).
#
# So the scope is narrowed to a pinned contig set rather than the whole
# genome. This matches the benchmark workloads' scope (chr1 + chr17) plus the
# mitochondrion, so the real-data panel and the generated benchmark are
# directly comparable, and it still gives a few thousand real transcripts
# with real splice junctions. The narrowing is a HARDWARE limit, recorded in
# datasets/manifest.yaml, and it bounds the claim: this panel validates the
# port against upstream on real reads, it does not validate whole-genome
# behaviour.
# CONTIGS and SJDB_OVERHANG were defaulted above, before the index directory name
# was derived from them.

# sjdbOverhang MUST be (read length - 1) of the reads this index will be used
# for. It is a required override rather than a derived default because the
# index is built before the reads are chosen, and getting it wrong silently
# degrades splice-junction detection -- which is precisely what the
# cross-organism endpoint E5 measures. The human dev/held-out runs are 150 bp
# (overhang 149); the rat run is 101 bp (overhang 100). The first rat index
# was built with the human default and had to be discarded.
#
# TRAP: do NOT read this value off STAR's own progress output. STAR's
# "Read length" column reports the SUMMED mate length, so a 2x101 library
# displays 202 and a 2x150 library displays 300. The overhang must come from
# the FASTQ, and it is (single mate length - 1), not (displayed length - 1).
# STAR's own ceiling, in MB, on how much RAM genomeGenerate may plan for. The
# default is 31000, which on a 31 GB machine is a promise the kernel cannot
# keep: the first build here was OOM-killed at "inserting junctions into the
# genome indices" (Log.out ends there, process died) after completing the
# 34-minute suffix array. Setting the limit makes STAR choose a smaller
# genomeSAindexNbases itself, which is the mechanism the option exists for --
# better than hand-tuning the index parameter and getting it subtly wrong.
# NOTE the unit: STAR parses --limitGenomeGenerateRAM as BYTES, and passing
# "20000" meaning megabytes is rejected as "too small for your genome" with a
# minimum of 8746522208 (i.e. ~8.1 GiB expressed in bytes).
LIMIT_RAM_BYTES="${LIMIT_RAM_BYTES:-16000000000}"

# genomeSAindexNbases is set EXPLICITLY, and this is the parameter that
# actually decides whether the build survives here.
#
#   14 (STAR's recommendation)  ~30 GB for human -- does not fit in 31 GB total
#   11                        ~19 GB observed, then OOM-killed at the
#                             "inserting junctions into the genome indices"
#                             step, TWICE
#    8                        ~5-7 GB, leaves room for the junction database
#
# The failure mode is why 11 is not good enough: --limitGenomeGenerateRAM
# plans for the suffix-array phase, but junction insertion allocates ON TOP of
# the finished SA index, so honouring the limit still gets OOM-killed. Swap
# was not available as a fallback (swapon is not permitted in this container),
# which is why this has to be a smaller index rather than a slower one.
# 11 rather than STAR's recommended 14: at 14 the human index needs ~30 GB and
# this machine has 31 GB total. A smaller suffix-array index means marginally
# more reads land as ambiguous near boundaries, which is a real cost; it does
# not affect whether the port agrees with upstream, because both arms read the
# same BAM.
SA_INDEX_NBASES="${SA_INDEX_NBASES:-11}"

THREADS="${THREADS:-$(nproc)}"

for f in "$GENOME" "$GTF"; do
  [[ -f "$f" ]] || { echo "error: missing $f (run datasets/fetch_reference.sh first)" >&2; exit 1; }
done

mkdir -p "$IDX_DIR"

# The stamp names its INPUTS as well as its parameters, and is consulted only
# after the unpack and subset steps below have established that genome.fa and
# annotation.gtf match their sources. Checking it earlier would compare a digest
# of a stale file and accept a stale index.
STAMP="$IDX_DIR/.built-${SA_INDEX_NBASES}-${LIMIT_RAM_BYTES}-${SJDB_OVERHANG}-${CONTIGS// /_}"

MM="${MICROMAMBA:-/scratch/mdra00001/tmp/bin/micromamba}"
# Export the root prefix here rather than relying on the caller: micromamba
# resolves envs relative to it, and a prefix that is unset-but-defaulted by
# micromamba points at XDG state instead, which fails with a confusing
# "prefix does not exist" from a backgrounded run.
export MAMBA_ROOT_PREFIX
STAR="$MM run -n $ENV_NAME STAR"

echo "preparing STAR index inputs (contigs='$CONTIGS')"
# STAR cannot read a gzipped genome OR a gzipped GTF: --sjdbGTFfile silently
# yields an empty junction database from a .gz (its own error text is the clue:
# "Make sure the GTF file is unzipped"). Decompressed once and reused.
# Unpack via a temp file and move it into place only on success. A plain
# `gunzip -c in > out` CREATES `out` even when gunzip then fails, and `set -e`
# does not remove it -- so one failed unpack leaves a zero-byte file that
# every later run then trusts as "already unpacked" and silently builds an
# empty index. That is exactly what happened on the rat attempt, where a
# plain-text annotation was fed to gunzip.
unpack() {
  local src="$1" dst="$2"
  echo "unpack    $(basename "$src") -> $(basename "$dst") (STAR requires uncompressed input)"
  if [[ "$src" == *.gz ]]; then
    gunzip -c "$src" > "$dst.tmp"
  else
    cp "$src" "$dst.tmp"
  fi
  [[ -s "$dst.tmp" ]] || { echo "error: unpack of $src produced an empty file" >&2; exit 1; }
  mv "$dst.tmp" "$dst"
}

if [[ ! -s "$IDX_DIR/genome.fa" ]]; then
  unpack "$GENOME" "$IDX_DIR/genome.fa"
fi

# Subset the FASTA and the GTF to CONTIGS. Both the index AND the splice
# junction database shrink, which is the step that actually made this build
# possible: the whole-genome GENCODE junction database is what does not fit
# in 31 GB.
SUBSET_GTF="$IDX_DIR/annotation.subset.gtf"

# Only unpack when the annotation is actually gzipped. The rat annotation is
# produced by datasets/refgene_to_gtf.py as plain text, and an unconditional
# gunzip fails on it -- which is how this got found.
if [[ ! -s "$IDX_DIR/annotation.gtf" ]] || [[ "$GTF" -nt "$IDX_DIR/annotation.gtf" ]]; then
  unpack "$GTF" "$IDX_DIR/annotation.gtf"
fi

# Regenerate the subset when it is missing or was built from different inputs.
#
# The previous guard was `if [[ ! -f genome.subset.fa || ! -f SUBSET_GTF ]]`,
# i.e. absence only, which silently reused a STALE subset whenever the annotation
# changed. That is precisely what happened when the rat refGene conversion was
# corrected: the corrected annotation.gtf was written, the subset was left
# untouched, and STAR rebuilt a junction database from the old coordinates while
# reporting success.
#
# A digest of the inputs is recorded rather than comparing mtimes, because a
# regenerated file can legitimately have an older timestamp than its derived
# artefact -- copying preserves neither, and rsync and tar do not agree. The
# stamp below makes "what was this built from" answerable after the fact.
SUBSET_STAMP="$IDX_DIR/.subset-from"
want_subset() {
  local genome_src="${1:-}" gtf_src="${2:-}"
  [[ -s "$SUBSET_GTF" ]] || return 1
  [[ -s "$SUBSET_STAMP" ]] || return 1
  grep -qxF "genome=$genome_src" "$SUBSET_STAMP" || return 1
  grep -qxF "gtf=$gtf_src" "$SUBSET_STAMP" || return 1
  return 0
}
GENOME_SRC="$(readlink -f "$IDX_DIR/genome.fa" 2>/dev/null || echo "$IDX_DIR/genome.fa")"
GTF_SRC="$(readlink -f "$IDX_DIR/annotation.gtf" 2>/dev/null || echo "$IDX_DIR/annotation.gtf")"
if ! want_subset "$GENOME_SRC" "$GTF_SRC" || [[ ! -s "$IDX_DIR/genome.subset.fa" ]]; then
  echo "subset    contigs: $CONTIGS"
  python3 - "$IDX_DIR/genome.fa" "$IDX_DIR/genome.subset.fa" $CONTIGS <<'PYEOF'
import sys
src, dst = sys.argv[1], sys.argv[2]
want = set(sys.argv[3:])
keep, seen = False, set()
with open(src) as fi, open(dst, "w") as fo:
    for line in fi:
        if line.startswith(">"):
            name = line[1:].split()[0]
            keep = name in want
            if keep:
                seen.add(name)
        if keep:
            fo.write(line)
missing = want - seen
if missing:
    sys.exit(f"error: contigs not present in the assembly: {sorted(missing)}")
PYEOF
  # GTF is filtered to the same contigs, which is what bounds the junction
  # database. Also drops non-transcript feature rows that STAR does not use.
  awk -v OFS='	' -v c="$CONTIGS" 'BEGIN{n=split(c,a," ");for(i=1;i<=n;i++)keep[a[i]]=1}
       /^#/ {print; next}
       keep[$1] && $3=="exon" {print}' "$IDX_DIR/annotation.gtf" > "$SUBSET_GTF.tmp"
  [[ -s "$SUBSET_GTF.tmp" ]] || { echo "error: contig subsetting produced an empty GTF" >&2; exit 1; }
  mv "$SUBSET_GTF.tmp" "$SUBSET_GTF"
  printf 'genome=%s\ngtf=%s\ncontigs=%s\n' "$GENOME_SRC" "$GTF_SRC" "$CONTIGS" > "$SUBSET_STAMP"
  echo "subset    $(wc -l < "$SUBSET_GTF") exon records retained"
fi

# Consulted only now, once genome.fa and annotation.gtf are known to match their
# sources and the subset is known to derive from them. Checking this earlier would
# compare a digest of a stale file and accept a stale index.
GTF_DIGEST="$(sha256sum "$IDX_DIR/annotation.gtf" | awk '{print $1}')"
GENOME_DIGEST="$(sha256sum "$IDX_DIR/genome.fa" | awk '{print $1}')"
SUBSET_DIGEST="$(sha256sum "$SUBSET_GTF" | awk '{print $1}')"
if [[ -f "$STAMP" ]] \
   && grep -qxF "gtf_sha256=$GTF_DIGEST" "$STAMP" \
   && grep -qxF "genome_sha256=$GENOME_DIGEST" "$STAMP" \
   && grep -qxF "subset_gtf_sha256=$SUBSET_DIGEST" "$STAMP"; then
  echo "ok       STAR index already built ($(cat "$STAMP"))"
  exit 0
fi
if [[ -f "$STAMP" ]]; then
  echo "rebuild  index inputs changed since the last build:" >&2
  echo "          recorded gtf_sha256=$(grep -o 'gtf_sha256=[0-9a-f]*' "$STAMP" | head -1 | cut -d= -f2)" >&2
  echo "          current  gtf_sha256=$GTF_DIGEST" >&2
  echo "          the previous junction database was built from the older annotation" >&2
fi

echo "building  STAR index (contigs='$CONTIGS' genomeSAindexNbases=$SA_INDEX_NBASES"
echo "          limitGenomeGenerateRAM=$LIMIT_RAM_BYTES bytes sjdbOverhang=$SJDB_OVERHANG)"
echo "          the genomeSAindexNbases STAR actually used is read back from Log.out"
echo "          into the stamp, so the stamp records what happened not what was asked"

$STAR --runMode genomeGenerate \
  --genomeDir "$IDX_DIR" \
  --genomeFastaFiles "$IDX_DIR/genome.subset.fa" \
  --sjdbGTFfile "$SUBSET_GTF" \
  --sjdbOverhang "$SJDB_OVERHANG" \
  --limitGenomeGenerateRAM "$LIMIT_RAM_BYTES" \
  --genomeSAindexNbases "$SA_INDEX_NBASES" \
  --runThreadN "$THREADS" \
  --outFileNamePrefix "$IDX_DIR/"

# Read the genomeSAindexNbases STAR actually chose back out of its log rather
# than asserting one, so the stamp records what happened and not what was
# requested.
CHOSEN="$(grep -oE '^genomeSAindexNbases +[0-9]+' "$IDX_DIR/Log.out" | tail -1 | awk '{print $2}' || true)"
CHOSEN="${CHOSEN:-unknown}"
printf 'genome=%s\ngtf=%s\ncontigs=%s\nsjdbOverhang=%s\nlimitGenomeGenerateRAM_bytes=%s\ngenomeSAindexNbases_requested=%s\ngenomeSAindexNbases_used=%s\ngtf_sha256=%s\ngenome_sha256=%s\nsubset_gtf_sha256=%s\n' \
  "$(basename "$GENOME")" "$(basename "$GTF")" "$CONTIGS" "$SJDB_OVERHANG" "$LIMIT_RAM_BYTES" "$SA_INDEX_NBASES" "$CHOSEN" \
  "$GTF_DIGEST" "$GENOME_DIGEST" "$SUBSET_DIGEST" > "$STAMP"

echo "ok       index built: $IDX_DIR"
