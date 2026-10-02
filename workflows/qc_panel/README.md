# Bulk RNA-seq QC panel, as a Snakemake workflow

A five-command QC panel wired into a real workflow engine, executed as a test.

## What it is

One rule per command, for the panel this project measures end to end:

```
bam_stat  ->  infer_experiment  ->  read_distribution  ->  read_NVC  ->  read_quality
                                                        \___________________________/
                                                                 panel.complete
```

The five stages are independent; only the summary marker depends on all of them,
which is what makes it usable as a dependency for a downstream rule.

| file | purpose |
|---|---|
| `Snakefile` | the workflow |
| `config/config.yaml` | sample panel, gene model, output directory |
| `test/run_panel.sh` | runs the workflow and diffs it against direct invocation |

## Run it

```bash
snakemake --cores 2 --configfile config/config.yaml
```

Every output lands under `{outdir}/{sample}/`. A QC report whose path does not name
its sample is a report that will eventually be compared against the wrong sample.

## Test it

```bash
test/run_panel.sh [BAM]
```

This is the part that matters. It runs the whole DAG through Snakemake on a real
alignment, re-runs the same five commands directly, and diffs every artifact — so a
wrapper bug fails the test instead of quietly changing somebody's QC numbers:

- `bam_stat.txt`, `infer_experiment.txt`, `read_distribution.txt` and
  `read_NVC.NVC.xls` must be **byte-identical** to direct invocation.
- `read_quality.qual.r` is compared after normalising the run directory, because a
  generated R script embeds its own output paths in `pdf('...')` calls. Everything
  else in it — the counts — stays byte-exact.
- `panel.complete` must be non-empty and every metric in it a positive count.
- No stage may invoke `Rscript`: a QC panel that needs R on every node fails on the
  first headless sample.

It needs `snakemake`, the release binaries on `PATH` reachable as `rseqc <name>.py`,
and an indexed BAM with a `.bai` sidecar. It **skips** rather than fails when no
alignment is present, so it is safe in a fresh checkout — but a skip is not a pass,
and CI provisions the alignment before relying on it.

## What it caught while being written

Every one of these was a real error in this file, found by running it:

- `infer_experiment` and `read_distribution` need `-r <gene model>`. The first draft
  omitted it and both rules failed with an argument error.
- `-o X` produces `X.NVC.xls`, not `X`. A rule whose `output:` does not name the file
  the shell line creates fails with `MissingOutputException` *after* the command has
  already succeeded.
- With `--skip-plot`, `read_quality`'s only deliverable is `X.qual.r`, which embeds
  the counts inline. It writes no `.xls` and no separate data table, in the port or in
  upstream. A pipeline waiting for a table will wait forever.
- awk runs its main block once per **input line**, and this marker command is given no
  input files — so a `print` outside `BEGIN` never executes and the marker comes out
  empty. The checker then reported success over a zero-line file, because a loop over
  no lines checks nothing. Both halves are fixed: the prints are in `BEGIN`, and the
  test now refuses an empty marker.
- awk's `getline var < file` returns the **line**, not a count. A `wc()` helper built on
  it stored file *content* in every "count" field. Counts are measured by `wc` and
  passed to awk with `-v`.

## Scope

This is a QC panel, not a complete RNA-seq pipeline: it consumes alignments that
already exist. Any end-to-end saving attributable to this project is therefore bounded
by the QC fraction of a real pipeline — see `docs/ENVELOPE.md` and
`benchmarks/pipeline-impact-rat-8M.json` for the measured figure.

The commands `read_NVC`, `read_quality` and the two profile commands emit an R script
for their heat maps. `--skip-plot` is used here so the panel runs headless; a
pipeline that wants the plots needs `Rscript` available, which is a real dependency
and not something this example hides.