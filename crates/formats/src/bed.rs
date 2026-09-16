//! BED12 gene-model region extraction, ported from `qcmodule.BED.ParseBED`
//! (`oracle/upstream-src/src/qcmodule/BED.py`). Each upstream method has
//! its OWN robustness level -- do not unify them into one shared parser:
//!
//! - [`get_cds_exon`] (ports `getCDSExon`, line 513): NO comment-line skip,
//!   NO exception handling at all. A `#`/`track`/`browser` line, or any
//!   malformed line, crashes upstream outright. Returns `io::Error` here
//!   for the same reason (fail fast, don't silently skip), matching the
//!   policy already used for `RNA_fragment_size.py`'s malformed-BED case.
//! - [`get_utr`] (ports `getUTR`, line 417) and [`get_intergenic`] (ports
//!   `getIntergenic`, line 718): skip `#`/`track`/`browser` lines
//!   explicitly, but have no exception handling for anything else
//!   malformed -- also effectively fail-fast beyond the comment skip.
//! - [`get_intron`] (ports `getIntron`, line 620): skips comments AND
//!   wraps the rest in upstream's `try/except: print(note); continue` --
//!   the only one of these four that gracefully skips a bad data line
//!   rather than erroring the whole run. Returns `(Vec<Bed3>, u64)` (value,
//!   skipped-line count) instead of a bare `Result`.

use std::io::{self, BufRead};

use crate::interval::Bed3;

fn is_comment_or_header(line: &str) -> bool {
    line.starts_with('#') || line.starts_with("track") || line.starts_with("browser")
}

struct Bed12Fields {
    chrom: String,
    // Only consumed while building `exon_starts` below; kept as named
    // locals rather than struct fields to avoid dead-code warnings, since
    // no current extraction method needs the transcript span itself.
    strand: String,
    cds_start: i64,
    cds_end: i64,
    block_sizes: Vec<i64>,
    /// Absolute (tx_start + relative) exon starts.
    exon_starts: Vec<i64>,
}

fn parse_bed12_fields(fields: &[&str]) -> Result<Bed12Fields, String> {
    if fields.len() < 12 {
        return Err(format!("has {} columns; expected 12", fields.len()));
    }
    let chrom = fields[0].to_string();
    let tx_start: i64 = fields[1].parse().map_err(|_| "invalid txStart".to_string())?;
    let _tx_end: i64 = fields[2].parse().map_err(|_| "invalid txEnd".to_string())?;
    let strand = fields[5].to_string();
    let cds_start: i64 = fields[6].parse().map_err(|_| "invalid cdsStart".to_string())?;
    let cds_end: i64 = fields[7].parse().map_err(|_| "invalid cdsEnd".to_string())?;
    let block_sizes: Vec<i64> = fields[10]
        .trim_end_matches(',')
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().map_err(|_| "invalid blockSizes".to_string()))
        .collect::<Result<_, _>>()?;
    let relative_starts: Vec<i64> = fields[11]
        .trim_end_matches(',')
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().map_err(|_| "invalid blockStarts".to_string()))
        .collect::<Result<_, _>>()?;
    let exon_starts: Vec<i64> = relative_starts.iter().map(|&rs| tx_start + rs).collect();

    Ok(Bed12Fields { chrom, strand, cds_start, cds_end, block_sizes, exon_starts })
}

/// Ports `ParseBED.getCDSExon`. Fail-fast: no comment skip, no recovery.
pub fn get_cds_exon(reader: impl BufRead) -> io::Result<Vec<Bed3>> {
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split_whitespace().collect();
        let f = parse_bed12_fields(&fields)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let exon_ends: Vec<i64> = f.exon_starts.iter().zip(f.block_sizes.iter()).map(|(&s, &sz)| s + sz).collect();
        for (&base, &end) in f.exon_starts.iter().zip(exon_ends.iter()) {
            if end < f.cds_start {
                continue;
            }
            if base > f.cds_end {
                continue;
            }
            let exon_start = base.max(f.cds_start);
            let exon_end = end.min(f.cds_end);
            out.push((f.chrom.clone(), exon_start, exon_end));
        }
    }
    Ok(out)
}

/// Ports `ParseBED.getUTR(utr=3|5)`. `utr` selects which end: `3` or `5`
/// (upstream's default `35` extracting both isn't used by any ported
/// command yet, so it's not implemented here -- add it if/when needed).
pub fn get_utr(reader: impl BufRead, utr: u8) -> io::Result<Vec<Bed3>> {
    assert!(utr == 3 || utr == 5, "get_utr only supports 3 or 5 (not upstream's combined 35 mode)");
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if is_comment_or_header(&line) {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        let f = parse_bed12_fields(&fields)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let exon_ends: Vec<i64> = f.exon_starts.iter().zip(f.block_sizes.iter()).map(|(&s, &sz)| s + sz).collect();

        let want_5 = (f.strand == "+" && utr == 5) || (f.strand == "-" && utr == 3);
        let want_3 = (f.strand == "+" && utr == 3) || (f.strand == "-" && utr == 5);

        if want_5 {
            for (&st, &end) in f.exon_starts.iter().zip(exon_ends.iter()) {
                if st < f.cds_start {
                    out.push((f.chrom.clone(), st, end.min(f.cds_start)));
                }
            }
        }
        if want_3 {
            for (&st, &end) in f.exon_starts.iter().zip(exon_ends.iter()) {
                if end > f.cds_end {
                    out.push((f.chrom.clone(), st.max(f.cds_end), end));
                }
            }
        }
    }
    Ok(out)
}

/// Ports `ParseBED.getIntron`. Returns `(regions, skipped_line_count)`;
/// the intron "label" (UTR5/UTR3/CDS) upstream also computes isn't used by
/// any ported command yet (they all go through `unionBed3`, which drops
/// it) -- add it back if a future command needs it.
pub fn get_intron(reader: impl BufRead) -> io::Result<(Vec<Bed3>, u64)> {
    let mut out = Vec::new();
    let mut skipped = 0u64;
    for line in reader.lines() {
        let line = line?;
        if is_comment_or_header(&line) {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        match parse_bed12_fields(&fields) {
            Ok(f) => {
                let exon_ends: Vec<i64> = f.exon_starts.iter().zip(f.block_sizes.iter()).map(|(&s, &sz)| s + sz).collect();
                if f.exon_starts.len() < 2 {
                    continue; // no introns possible with 0 or 1 exon
                }
                let intron_starts = &exon_ends[..exon_ends.len() - 1];
                let intron_ends = &f.exon_starts[1..];
                for (&st, &end) in intron_starts.iter().zip(intron_ends.iter()) {
                    out.push((f.chrom.clone(), st, end));
                }
            }
            Err(_) => skipped += 1,
        }
    }
    Ok((out, skipped))
}

/// Ports `ParseBED.getIntergenic`. `direction` is `"up"` or `"down"`
/// (upstream's `"both"` mode isn't used by any ported command yet).
pub fn get_intergenic(reader: impl BufRead, direction: &str, size: i64) -> io::Result<Vec<Bed3>> {
    assert!(direction == "up" || direction == "down", "get_intergenic only supports up or down (not upstream's combined both mode)");
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if is_comment_or_header(&line) {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 6 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, format!("has {} columns; expected at least 6", fields.len())));
        }
        let chrom = fields[0].to_string();
        let tx_start: i64 = fields[1].parse().map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "invalid txStart"))?;
        let tx_end: i64 = fields[2].parse().map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "invalid txEnd"))?;
        let strand = fields[5];

        let upstream_is_minus_shaped = (direction == "up") == (strand == "-");
        let (region_st, region_end) = if upstream_is_minus_shaped {
            (tx_end, tx_end + size)
        } else {
            ((tx_start - size).max(0), tx_start)
        };
        out.push((chrom, region_st, region_end));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "chr1\t100\t500\tgeneA\t0\t+\t150\t450\t0\t2\t100,100,\t0,300,\n";

    #[test]
    fn cds_exon_clips_to_thick_start_end() {
        let out = get_cds_exon(SAMPLE.as_bytes()).unwrap();
        // exon1: [100,200) clipped to cdsStart 150 -> [150,200)
        // exon2: [400,500) clipped to cdsEnd 450 -> [400,450)
        assert_eq!(out, vec![("chr1".to_string(), 150, 200), ("chr1".to_string(), 400, 450)]);
    }

    #[test]
    fn cds_exon_errors_on_comment_line() {
        let err = get_cds_exon("# comment\n".as_bytes()).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn utr5_and_utr3_forward_strand() {
        let utr5 = get_utr(SAMPLE.as_bytes(), 5).unwrap();
        let utr3 = get_utr(SAMPLE.as_bytes(), 3).unwrap();
        // + strand: 5'UTR from exon1 start<cdsStart -> [100,150); 3'UTR from exon2 end>cdsEnd -> [450,500)
        assert_eq!(utr5, vec![("chr1".to_string(), 100, 150)]);
        assert_eq!(utr3, vec![("chr1".to_string(), 450, 500)]);
    }

    #[test]
    fn utr_strand_minus_swaps_5_and_3() {
        let minus = "chr1\t100\t500\tgeneA\t0\t-\t150\t450\t0\t2\t100,100,\t0,300,\n";
        let utr5 = get_utr(minus.as_bytes(), 5).unwrap();
        let utr3 = get_utr(minus.as_bytes(), 3).unwrap();
        assert_eq!(utr5, vec![("chr1".to_string(), 450, 500)]);
        assert_eq!(utr3, vec![("chr1".to_string(), 100, 150)]);
    }

    #[test]
    fn intron_between_exons_and_skips_malformed_line() {
        let text = format!("{SAMPLE}bad line\n");
        let (introns, skipped) = get_intron(text.as_bytes()).unwrap();
        assert_eq!(introns, vec![("chr1".to_string(), 200, 400)]);
        assert_eq!(skipped, 1);
    }

    #[test]
    fn intergenic_upstream_downstream_plus_strand() {
        let up = get_intergenic(SAMPLE.as_bytes(), "up", 1000).unwrap();
        let down = get_intergenic(SAMPLE.as_bytes(), "down", 1000).unwrap();
        // + strand: up = [max(100-1000,0),100) = [0,100); down = [500,1500)
        assert_eq!(up, vec![("chr1".to_string(), 0, 100)]);
        assert_eq!(down, vec![("chr1".to_string(), 500, 1500)]);
    }

    #[test]
    fn intergenic_minus_strand_swaps_up_down() {
        let minus = "chr1\t100\t500\tgeneA\t0\t-\t150\t450\t0\t2\t100,100,\t0,300,\n";
        let up = get_intergenic(minus.as_bytes(), "up", 1000).unwrap();
        let down = get_intergenic(minus.as_bytes(), "down", 1000).unwrap();
        assert_eq!(up, vec![("chr1".to_string(), 500, 1500)]);
        assert_eq!(down, vec![("chr1".to_string(), 0, 100)]);
    }
}
