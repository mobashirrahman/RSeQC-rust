#!/usr/bin/env bash
# Re-check DIV-0024: does geneBody_coverage still disagree with upstream on real data?
#
# The exit status is INVERTED on purpose. This script exits 0 when the divergence is
# still present, because the divergence entry claims it is; it exits 1 when the
# divergence has gone, because then `compatibility/divergences.yaml` and the CI step
# that calls this are both stale and must be updated in the same change. A recorded
# divergence nobody re-checks becomes a stale claim, which is how DIV-0005 was found
# to be out of date.
#
# Needs the pinned oracle and a provisioned real-data panel:
#   oracle/venv/bin/python3 datasets/fetch_fastq.py SRR1177982
#   datasets/build_star_index.sh && datasets/align_run.sh SRR1177982
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

BAM="${RSEQC_RAT_BAM:-datasets/heldout/aligned/SRR1177982/SRR1177982.bam}"
BED="verification/fixtures/genebody_divergence_minimal.bed12"
[[ -f "$BAM" ]] || { echo "skip: $BAM not provisioned"; exit 0; }
[[ -x oracle/venv/bin/python3 ]] || { echo "skip: pinned oracle absent"; exit 0; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

PYTHONPATH=oracle/upstream-src/src oracle/venv/bin/python3 \
  oracle/upstream-src/scripts/geneBody_coverage.py \
  -i "$BAM" -r "$BED" --out-prefix "$WORK/py" --skip-plot >/dev/null 2>&1
./target/release/geneBody_coverage \
  -i "$BAM" -r "$BED" --out-prefix "$WORK/rs" --skip-plot >/dev/null 2>&1

python3 - "$WORK" <<'PYEOF'
import sys
from pathlib import Path

work = Path(sys.argv[1])
py = [float(x) for x in (work / "py.geneBodyCoverage.txt").read_text().splitlines()[1].split("\t")[1:]]
rs = [float(x) for x in (work / "rs.geneBodyCoverage.txt").read_text().splitlines()[1].split("\t")[1:]]
differing = sum(1 for a, b in zip(py, rs) if a != b)
worst = max(abs(a - b) for a, b in zip(py, rs))
print(f"DIV-0024: {differing}/100 bins differ, largest difference {worst:.0f} reads")
if differing >= 50:
    print("  the divergence is still present and the ledger entry is accurate")
    sys.exit(0)
print("  the divergence has CHANGED or gone. Update compatibility/divergences.yaml")
print("  (DIV-0024) and the CI step that calls this script in the same change.")
sys.exit(1)
PYEOF
