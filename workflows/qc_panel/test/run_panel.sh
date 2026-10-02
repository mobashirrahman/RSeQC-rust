#!/usr/bin/env bash
# Execute the QC panel through Snakemake and prove it agrees with direct invocation.
#
# A workflow example that is never run is documentation. This runs the whole DAG on a
# real alignment, then re-runs the same five commands directly and diffs every
# artifact. That catches the failures a Snakemake file can have while every command
# still works perfectly on its own:
#
#   * a wrong flag (a -o passed to read_distribution, which accepts none),
#   * an output declared in `output:` that the shell line never writes, which
#     Snakemake treats as a transient error and retries until it gives up,
#   * stdout redirected to a file that the command does not actually print to,
#   * an Rscript dependency in a QC panel that must run headless,
#   * a completion marker that disagrees with the reports it summarises.
#
# Usage:  test/run_panel.sh [BAM]
# Needs:  snakemake, the release binaries on PATH as `rseqc`, and an indexed BAM.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW="$(dirname "$HERE")"
BAM="${1:-/scratch/mdra00001/RSeQC-rust/datasets/heldout/aligned/SRR1177982/SRR1177982.bam}"

if [[ ! -f "$BAM" ]]; then
  echo "skip: no alignment at $BAM" >&2
  echo "      fetch one with: datasets/fetch_fastq.py SRR1177982" >&2
  exit 0
fi
if [[ ! -f "$BAM.bai" ]]; then
  echo "FAIL: $BAM has no .bai sidecar; this workflow requires indexed alignments" >&2
  exit 1
fi
if ! command -v snakemake >/dev/null 2>&1; then
  echo "FAIL: snakemake is not on PATH. This example is not 'tested' without it:" >&2
  echo "      a workflow that is never executed is documentation, not a test." >&2
  exit 1
fi
if ! command -v rseqc >/dev/null 2>&1 && ! command -v bam_stat >/dev/null 2>&1; then
  echo "FAIL: neither 'rseqc' nor the release binaries are on PATH" >&2
  exit 1
fi

RUNNER=(rseqc)
command -v rseqc >/dev/null 2>&1 || RUNNER=()

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

echo "== 1. run the workflow through Snakemake"
cd "$WORKFLOW"
REFGENE="${REFGENE:-/scratch/mdra00001/RSeQC-rust/datasets/heldout/reference/rn6.bed12}"
if [[ ! -f "$REFGENE" ]]; then
  echo "FAIL: no gene model at $REFGENE; infer_experiment and read_distribution require one" >&2
  exit 1
fi

cat > config/test.yaml <<EOF
samples:
  panel_sample: $BAM
refgene: $REFGENE
outdir: $WORK/snakemake
mapq: 30
EOF
snakemake --cores 2 --configfile config/test.yaml \
          --directory "$WORK/snakemake" \
          --printshellcmds > "$WORK/snakemake.log" 2>&1 || {
  echo "FAIL: snakemake did not complete" >&2
  tail -40 "$WORK/snakemake.log" >&2
  exit 1
}

echo "== 2. run the same five commands directly"
mkdir -p "$WORK/direct"
"${RUNNER[@]}" bam_stat.py       -i "$BAM" -q 30          > "$WORK/direct/bam_stat.txt"       2>/dev/null
"${RUNNER[@]}" infer_experiment.py -i "$BAM" -r "$REFGENE"  > "$WORK/direct/infer_experiment.txt" 2>/dev/null
"${RUNNER[@]}" read_distribution.py -i "$BAM" -r "$REFGENE" > "$WORK/direct/read_distribution.txt" 2>/dev/null
"${RUNNER[@]}" read_NVC.py       -i "$BAM" -o "$WORK/direct/read_NVC"     --skip-plot 2>/dev/null
"${RUNNER[@]}" read_quality.py   -i "$BAM" -o "$WORK/direct/read_quality" --skip-plot 2>/dev/null

echo "== 3. diff every artifact the workflow claims to produce"
status=0
for f in bam_stat.txt infer_experiment.txt read_distribution.txt \
         read_NVC.NVC.xls read_quality.qual.r; do
  a="$WORK/snakemake/panel_sample/$f"
  b="$WORK/direct/$f"
  if [[ ! -f "$a" ]]; then
    echo "FAIL  $f: the workflow did not produce the file it declared" >&2
    status=1
  elif [[ "$f" == *.r ]]; then
    # A generated R script embeds its OWN output directory in pdf('...') calls, so
    # two runs writing to two directories cannot be byte-identical. The data is
    # compared instead: strip only the run directory, which is the named,
    # per-artifact normalisation this project's differential suite also uses, and
    # everything else -- the counts, the matrix, the layout -- stays byte-exact.
    if sed "s|$WORK/snakemake/panel_sample|<OUT>|g" "$a" \
       | diff -q - <(sed "s|$WORK/direct|<OUT>|g" "$b") >/dev/null; then
      echo "  ok  $f ($(wc -c < "$a") bytes, identical after normalising the run directory)"
    else
      echo "FAIL  $f: differs from direct invocation even with run-directory paths normalised" >&2
      diff <(sed "s|$WORK/snakemake/panel_sample|<OUT>|g" "$a") \
           <(sed "s|$WORK/direct|<OUT>|g" "$b") | head -10 >&2
      status=1
    fi
  elif ! diff -q "$a" "$b" >/dev/null; then
    echo "FAIL  $f: the workflow's output differs from direct invocation" >&2
    diff "$a" "$b" | head -10 >&2
    status=1
  else
    echo "  ok  $f ($(wc -c < "$a") bytes identical)"
  fi
done

echo "== 4. check the completion marker agrees with the reports"
marker="$WORK/snakemake/panel_sample/panel.complete"
if [[ ! -f "$marker" ]]; then
  echo "FAIL  panel.complete is missing, so nothing downstream can depend on the panel" >&2
  status=1
elif [[ ! -s "$marker" ]]; then
  # An EMPTY marker is not a passing marker. A loop over a zero-line file checks
  # nothing and reports success, which is how an awk whose prints lived in the main
  # block -- never executed without input files -- passed this check while writing
  # nothing at all.
  echo "FAIL  panel.complete is empty; a summary with no lines summarises nothing" >&2
  status=1
else
  # Each row is a metric name and a count. A zero or non-numeric value means the
  # marker is reporting nothing about a stage that ran, so both are refused: a
  # summary that silently records zero is worse than no summary, because it looks
  # like a measurement.
  rows=$(grep -c . "$marker" || true)
  while IFS=$'\t' read -r key value; do
    [[ -z "$key" || "$key" == "sample" ]] && continue
    if [[ ! "$value" =~ ^[0-9]+$ ]] || [[ "$value" -le 0 ]]; then
      echo "FAIL  panel.complete records '$key' as '$value', which is not a positive count" >&2
      status=1
    fi
  done < "$marker"
  if [[ "$rows" -lt 6 ]]; then
    echo "FAIL  panel.complete has $rows lines; five metrics plus a sample header is 6" >&2
    status=1
  else
    echo "  ok  panel.complete ($rows rows, every metric a positive count)"
  fi
fi

echo "== 5. check no stage silently required Rscript"
# A QC panel that needs R on every node fails on the first headless sample. The
# rules pass --skip-plot precisely so that it does not.
if grep -qE 'Rscript|rpy2' "$WORK/snakemake/panel_sample/logs/"*.log 2>/dev/null; then
  echo "FAIL  a stage log mentions Rscript despite --skip-plot" >&2
  status=1
else
  echo "  ok  no stage invoked Rscript"
fi

if [[ $status -eq 0 ]]; then
  echo "PASS: the Snakemake panel reproduces direct invocation on all five artifacts"
else
  echo "FAILED: see above"
fi
exit $status