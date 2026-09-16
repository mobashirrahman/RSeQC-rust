# Oracle `--help` capture — summary

Date: 2026-09-16
Venv: `oracle/venv/` (Python 3.13.13). Activate with `source oracle/venv/bin/activate`.
Result: **33/33 scripts produced a successful `--help` capture (exit code 0).**
Zero failures.

## Method

Each script was invoked directly as
`oracle/venv/bin/python oracle/upstream-src/scripts/<name>.py --help`
(stdout+stderr merged), with the capture file prefixed by a `# exit=<code>`
header. The package was also installed editable (`pip install -e
oracle/upstream-src`), so `qcmodule` resolves on import and console scripts are
on the venv PATH as an alternative invocation route.

## Installed dependencies

Resolved by grepping `import`/`from` statements across
`scripts/*.py` and `src/qcmodule/*.py` (setup.py declares no install_requires —
it only lists `scripts=`):

| package     | version | note |
|-------------|---------|------|
| numpy       | 2.5.3   | wheel |
| pysam       | 0.24.1  | wheel |
| pandas      | 3.0.5   | wheel |
| matplotlib  | 3.11.2  | wheel |
| pyBigWig    | 0.3.26  | wheel |
| bx-python   | 0.15.1  | wheel |
| logomaker   | 0.8.7   | wheel (dep of sc_seqLogo) |

All installed from prebuilt manylinux wheels — **no package failed to
install**; no build-from-source, no system-library issues, no pip failures to
report.

## Per-script status

All 33 scripts listed in `oracle/upstream-src/setup.py` succeeded:

FPKM-UQ, FPKM_count, RNA_fragment_size, RPKM_saturation, bam2fq, bam2wig,
bam_stat, clipping_profile, deletion_profile, divide_bam, geneBody_coverage,
geneBody_coverage2, infer_experiment, inner_distance, insertion_profile,
junction_annotation, junction_saturation, mismatch_profile, normalize_bigwig,
overlay_bigwig, read_GC, read_NVC, read_distribution, read_duplication,
read_hexamer, read_quality, sc_bamStat, sc_editMatrix, sc_seqLogo, sc_seqQual,
split_bam, split_paired_bam, tin — each exit=0.

Captures live in `oracle/reference-runs/help/<script-name>.txt`. Note:
`scripts/update.pl` is a perl utility, not one of the 33 Python entry points,
and was intentionally skipped.

## Caveats for future oracle runs

- `qcmodule` was installed editable, so imports resolve via
  `oracle/upstream-src/src/qcmodule/` (reference code is executed, never
  imported by Rust production code).
- Full valid-workload runs are out of scope for this task.
