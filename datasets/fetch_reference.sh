#!/usr/bin/env bash
# Fetch and checksum-verify the reference inputs for the T4 real-data panel.
#
# Pinned by URL + SHA256. Re-running is idempotent: a file that already matches
# its recorded digest is left alone, so this is safe to call after an
# interrupted download. A file that exists but does NOT match is treated as a
# truncated download and refetched rather than trusted.
#
# Why GENCODE rather than UCSC's RefSeq: UCSC has renamed and moved its RefSeq
# GTF several times (the `genes/hg38.ncbiRefSeq.gtf.gz` path used by the
# benchmark generator's own docstring now 404s), and GENCODE publishes
# versioned releases under stable, immutable FTP paths. testing.md section 11.1
# requires the annotation version to be recorded; a release number is a stronger
# pin than a filename.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF_DIR="${REF_DIR:-$HERE/reference}"
mkdir -p "$REF_DIR"

# name|URL|sha256
ENTRIES=(
  "hg38.fa.gz|https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz|c1dd87068c254eb53d944f71e51d1311964fce8de24d6fc0effc9c61c01527d4"
  "gencode.v47.annotation.gtf.gz|https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz|df11938c66d2b39f8ebbdb5f1720919321f9cf9d58f02b019753f6a24ca3db24"
)

have() { command -v "$1" >/dev/null 2>&1; }
if have sha256sum; then SHA=sha256sum; elif have shasum; then SHA="shasum -a 256"; else
  echo "error: need sha256sum or shasum" >&2; exit 1
fi

for entry in "${ENTRIES[@]}"; do
  IFS='|' read -r name url want <<<"$entry"
  dest="$REF_DIR/$name"

  if [[ -f "$dest" ]]; then
    got="$($SHA "$dest" | awk '{print $1}')"
    if [[ "$got" == "$want" && "$want" != __* ]]; then
      echo "ok       $name (already present, digest matches)"
      continue
    fi
    if [[ "$want" == __* ]]; then
      echo "digest   $name -> $got   (no expected digest recorded yet; accepting)"
      continue
    fi
    echo "refetch  $name (present but digest differs -- treating as truncated)"
    rm -f "$dest"
  fi

  echo "fetch    $name"
  # -C - resumes a partial transfer; --fail turns an HTTP error into a non-zero
  # exit so a 404 cannot be silently written as a 40-byte "genome".
  curl --fail --location --retry 5 --retry-delay 5 -C - -o "$dest" "$url"

  got="$($SHA "$dest" | awk '{print $1}')"
  if [[ "$want" != __* && "$got" != "$want" ]]; then
    echo "error    $name digest $got != expected $want" >&2
    exit 1
  fi
  echo "ok       $name ($got)"
done

echo
echo "reference ready in $REF_DIR"
