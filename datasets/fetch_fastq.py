#!/usr/bin/env python3
"""Fetch FASTQ runs for the T4 real-data panel from ENA, with checksums.

Run accessions and their expected MD5s come from the ENA portal API, not from
this file, so a run's digest is always the archive's own. The *selection* of
runs -- which studies, which lanes, which split -- is recorded in
`datasets/manifest.yaml`; this script only fetches what the manifest names.

Design notes that matter for reproducibility:

  - ENA publishes one FASTQ per read for older submissions and a single
    interleaved-or-single file for newer ones. We handle both by reading the
    per-run `fastq_ftp`/`fastq_md5` fields, which are semicolon-separated and
    positionally aligned, rather than assuming a naming convention.
  - Downloads resume (`curl -C -`) and are verified against the archive's MD5
    before being accepted. A partial file is left in place for the next
    attempt rather than deleted, since a 3.4 GB run may be interrupted.
  - Nothing is deleted. `datasets/raw/` is gitignored, and a run already
    present with a matching MD5 is skipped, so re-running is cheap.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ENA_API = "https://www.ebi.ac.uk/ena/portal/api"
HERE = Path(__file__).resolve().parent

# Fields needed to (a) resolve each run's FASTQ URLs and (b) record the
# metadata testing.md section 11.1 requires in the manifest.
RUN_FIELDS = [
    "run_accession", "experiment_accession", "study_accession",
    "scientific_name", "library_layout", "library_strategy", "library_source",
    "library_selection", "instrument_model", "experiment_title",
    "read_count", "base_count", "fastq_ftp", "fastq_md5", "fastq_bytes",
    "submitted_ftp", "first_public", "last_updated",
]


def api(path: str, params: dict[str, str | int]) -> list[dict]:
    url = f"{ENA_API}/{path}?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=120) as fh:
        return json.load(fh)


def run_metadata(accession: str) -> dict:
    rows = api("filereport", {
        "accession": accession, "result": "read_run",
        "fields": ",".join(RUN_FIELDS), "format": "json", "limit": 1,
    })
    if not rows:
        raise SystemExit(f"error: ENA has no read_run for {accession}")
    return rows[0]


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path, expect_md5: str) -> None:
    if dest.is_file():
        got = md5(dest)
        if got == expect_md5:
            print(f"  ok       {dest.name} (present, md5 matches)")
            return
        print(f"  refetch  {dest.name} (present, md5 differs -> truncated?)")
    dest.parent.mkdir(parents=True, exist_ok=True)
    # --fail so a 404 cannot be written out as a 200-byte "fastq"; -C - to resume.
    subprocess.run(
        ["curl", "--fail", "--location", "--retry", "5", "--retry-delay", "5",
         "-C", "-", "-o", str(dest), url],
        check=True,
    )
    got = md5(dest)
    if got != expect_md5:
        raise SystemExit(f"error: {dest.name} md5 {got} != archive's {expect_md5}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+", help="ENA run accessions, e.g. SRR1177968")
    ap.add_argument("--out", default=str(HERE / "raw"))
    ap.add_argument("--meta-out", default=str(HERE / "logs" / "run_metadata.json"),
                    help="where to record the archive metadata testing.md 11.1 requires")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    meta_path = Path(args.meta_out)
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    metadata: dict = {}
    if meta_path.is_file():
        metadata = json.loads(meta_path.read_text())

    for acc in args.runs:
        print(f"=== {acc} ===")
        rec = run_metadata(acc)
        metadata[acc] = rec

        layout = rec.get("library_layout", "SINGLE")
        if layout != "PAIRED":
            print(f"  note     {acc} is {layout}; the T4 panel is paired-end only")
            continue
        if rec.get("submitted_ftp"):
            print("  note     archive also supplies aligned files; ignoring them on purpose (see manifest)")

        urls = (rec.get("fastq_ftp") or "").rstrip(";").split(";")
        digests = (rec.get("fastq_md5") or "").rstrip(";").split(";")
        if not urls or urls == [""]:
            raise SystemExit(f"error: {acc} has no public FASTQ")
        if len(urls) != len(digests):
            raise SystemExit(f"error: {acc} has {len(urls)} FASTQs but {len(digests)} md5s")

        for url, want in zip(urls, digests):
            dest = out / url.rsplit("/", 1)[-1]
            fetch(url, dest, want)

    meta_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(f"\nmetadata -> {meta_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
