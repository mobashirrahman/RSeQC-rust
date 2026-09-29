#!/usr/bin/env python3
"""Scaling sweeps along the implementation's actual cost drivers (protocol section 6).

One driver is varied at a time, all others held fixed, and each point is measured with
the same paired/randomised machinery as the main matrix. The point is the SHAPE of
the curve and where the knee is, not a single ratio: preliminary probing showed the
bam_stat speedup falling from 69x at 1k reads to 4.9x at 400k, so any single-point
number is uninformative.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PY = REPO / "oracle" / "venv" / "bin" / "python3"
GEN = REPO / "benchmarks" / "generate_workload_real.py"
BENCH = REPO / "benchmarks" / "bench.py"


def workload_complete(out: Path) -> bool:
    """A workload is usable only if the generator actually finished it.

    Checking only for reads.bam is not enough: a generator that fails its own
    post-write reference validation leaves a truncated, partial workload behind, and
    reusing it makes both arms run over missing input and agree perfectly on empty
    output -- a silent false pass.
    """
    need = ["reads.bam", "reads.bam.bai", "manifest.json", "model.bed12",
            "chrom.sizes", "reads_1.fastq", "reads_1.fa"]
    return all((out / f).exists() and (out / f).stat().st_size > 0 for f in need)


def gen(out: Path, real: Path, contigs, reads, seed, max_tx=0, read_len=100):
    if workload_complete(out):
        return
    if out.exists():
        import shutil as _sh
        _sh.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    cmd = [str(PY), str(GEN), "--refgene", str(real / "hg38.ncbiRefSeq.gtf.gz"),
           "--chrom-dir", str(real), "--contigs", *contigs,
           "--reads", str(reads), "--output-dir", str(out), "--seed", str(seed),
           "--read-length", str(read_len)]
    if max_tx:
        cmd += ["--max-transcripts", str(max_tx)]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    # Companion inputs the command matrix expects.
    import pysam
    fq, fa = out / "reads_1.fastq", out / "reads_1.fa"
    if not fq.exists():
        with pysam.AlignmentFile(str(out / "reads.bam")) as f, open(fq, "w") as a, open(fa, "w") as b:
            n = 0
            for r in f:
                if r.is_read1 and not r.is_unmapped:
                    a.write(f"@{r.query_name}\n{r.query_sequence}\n+\n"
                            + "".join(chr(33 + q) for q in r.query_qualities) + "\n")
                    b.write(f">{r.query_name}\n{r.query_sequence}\n")
                    n += 1
                if n >= 20000:
                    break
    if not (out / "sig1.bw").exists():
        import numpy as np, pyBigWig
        sizes = {}
        with open(out / "chrom.sizes") as fh:
            for line in fh:
                c, L = line.split()
                sizes[c] = int(L)
        for name in ("sig1.bw", "sig2.bw"):
            bw = pyBigWig.open(str(out / name), "w")
            bw.addHeader([(c, sizes[c]) for c in sorted(sizes)])
            for c in sorted(sizes):
                s = np.arange(0, 400000, 50)
                bw.addEntries([c] * len(s), s.tolist(), ends=(s + 20).tolist(),
                              values=[float(1 + i % 7) for i in range(len(s))])
            bw.close()
    # Only now is the workload complete: the BAM, its index, the manifest written by the
    # generator, and the companion inputs all exist and are non-empty.
    if not workload_complete(out):
        raise RuntimeError(f"workload generation incomplete: {out}")


def sweep(label, driver, points, real, workroot, commands, reps, base):
    """points: list of (value, kwargs-for-gen, extra-args-for-bench)."""
    out_rows = []
    for value, gkw, extra in points:
        wl = workroot / f"{label}_{driver}_{value}"
        print(f"\n### {label} / {driver} = {value}", flush=True)
        gen(wl, real, **gkw)
        od = workroot / f"res_{label}_{driver}_{value}"
        cmd = [str(PY), str(BENCH), "--workload", str(wl), "--commands", *commands,
               "--reps", str(reps), "--output-dir", str(od), "--label", f"{label}:{driver}={value}"]
        if extra:
            cmd += extra
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        d = json.loads((od / "results.json").read_text())
        for r in d["results"]:
            if "error" in r:
                continue
            w = r.get("wall_e2e", {})
            m = r.get("peak_rss_mb", {})
            out_rows.append({
                "sweep": label, "driver": driver, "value": value,
                "command": r["command"],
                "median_python": w.get("median_python"),
                "median_rust": w.get("median_rust"),
                "ratio": w.get("ratio_median"),
                "ci95": w.get("ratio_ci95"),
                "rss_py": m.get("median_python"),
                "rss_rs": m.get("median_rust"),
                "gate": r.get("equivalence", {}).get("pass"),
                "n_transcripts": (json.loads((wl / "manifest.json").read_text())
                                  .get("n_transcripts") if (wl / "manifest.json").exists() else None),
            })
            print(f"  {r['command']:20} py {w.get('median_python', 0):8.3f} "
                  f"rs {w.get('median_rust', 0):8.4f} x{w.get('ratio_median') or 0:7.2f} "
                  f"gate={r.get('equivalence', {}).get('pass')}", flush=True)
    return out_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", required=True)
    ap.add_argument("--workroot", required=True)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--which", nargs="+", default=["reads", "transcripts"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    real = Path(args.real)
    wr = Path(args.workroot)
    wr.mkdir(parents=True, exist_ok=True)
    rows = []

    if "reads" in args.which:
        pts = []
        for n in (2000, 10000, 50000, 200000, 800000):
            pts.append((n, dict(contigs=["chr17"], reads=n, seed=20260929), []))
        rows += sweep("A", "reads", pts, real, wr, ["bam_stat", "read_quality", "read_NVC"],
                      args.reps, None)

    if "transcripts" in args.which:
        pts = []
        for t in (200, 1000, 5000, 20000):
            pts.append((t, dict(contigs=["chr17"], reads=200000, seed=20260929,
                                max_tx=t), []))
        # read_distribution and FPKM_count are transcript-driven but stay inside the
        # time budget. tin is excluded here: upstream's per-transcript cost makes a
        # 20k-transcript point take minutes, and it is already covered by the main
        # matrix at a sized workload.
        rows += sweep("A", "transcripts", pts, real, wr,
                      ["read_distribution", "FPKM_count"], args.reps, None)

    if "readlen" in args.which:
        pts = []
        for L in (50, 75, 100, 150):
            pts.append((L, dict(contigs=["chr17"], reads=100000, seed=20260929,
                                read_len=L), []))
        rows += sweep("A", "readlen", pts, real, wr, ["read_quality", "read_hexamer"],
                      args.reps, None)

    Path(args.out).write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {len(rows)} rows -> {args.out}")


if __name__ == "__main__":
    main()
