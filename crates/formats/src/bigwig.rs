//! Thin wrapper over the `bigtools` crate's BigWig reader, exposing
//! only what `geneBody_coverage2.py`/`normalize_bigwig.py`/
//! `overlay_bigwig.py` need: chromosome list, per-base values over a
//! range (matching `pyBigWig.values()`), and raw intervals (matching
//! `pyBigWig.intervals()`). No BigWig WRITING support -- none of the
//! ported commands write BigWig, only WIG/bedGraph text or plot images.
//!
//! `bigtools` was chosen over hand-rolling a BigWig parser because,
//! unlike every other "format" this port has replicated so far (whose
//! upstream Python readers turned out to be simple hand-rolled
//! scanners, e.g. `read_hexamer.py`'s FASTA/FASTQ parsing), BigWig is a
//! genuinely complex bgzf-compressed, B+tree/R-tree-indexed binary
//! format that upstream itself reads via a mature C library
//! (`pyBigWig` wraps `libBigWig`) -- the same class of justified
//! dependency as this project's own use of `noodles-bam`/`noodles-sam`
//! for BAM/SAM. Verified byte-for-byte against real `pyBigWig` output
//! on a real BigWig file (see the test module).

use std::io;
use std::path::Path;

use bigtools::BigWigRead;

pub struct BigWigReader {
    inner: BigWigRead<bigtools::utils::reopen::ReopenableFile>,
}

impl BigWigReader {
    pub fn open(path: &Path) -> io::Result<Self> {
        let inner = BigWigRead::open_file(path).map_err(io::Error::other)?;
        Ok(Self { inner })
    }

    /// Returns `(chrom_name, length)` pairs in the file's own order.
    pub fn chroms(&self) -> Vec<(String, u32)> {
        self.inner.chroms().iter().map(|c| (c.name.clone(), c.length)).collect()
    }

    /// Returns per-base values for `[start, end)` (0-based, half-open),
    /// with `f32::NAN` at any position with no data. Matches
    /// `pyBigWig.values(chrom, start, end)` exactly.
    pub fn values(&mut self, chrom: &str, start: u32, end: u32) -> io::Result<Vec<f32>> {
        self.inner.values(chrom, start, end).map_err(io::Error::other)
    }

    /// Returns the raw `(start, end, value)` intervals intersecting
    /// `[start, end)`. Matches `pyBigWig.intervals(chrom, start, end)`.
    pub fn intervals(&mut self, chrom: &str, start: u32, end: u32) -> io::Result<Vec<(u32, u32, f32)>> {
        let iter = self.inner.get_interval(chrom, start, end).map_err(io::Error::other)?;
        iter.map(|r| r.map(|v| (v.start, v.end, v.value)).map_err(io::Error::other)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_path() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/pybigwig_test.bw")
    }

    /// Ground truth captured directly via `pyBigWig.open(...).chroms()`/
    /// `.values()`/`.intervals()` on this same fixture file (see
    /// `tests/fixtures/README.md`).
    #[test]
    fn chroms_match_pybigwig_ground_truth() {
        let bw = BigWigReader::open(&fixture_path()).unwrap();
        let chroms = bw.chroms();
        assert_eq!(chroms, vec![("1".to_string(), 195471971), ("10".to_string(), 130694993)]);
    }

    #[test]
    fn values_match_pybigwig_ground_truth() {
        let mut bw = BigWigReader::open(&fixture_path()).unwrap();
        let vals = bw.values("1", 0, 151).unwrap();
        assert_eq!(vals[0], 0.1);
        assert_eq!(vals[1], 0.2);
        assert_eq!(vals[2], 0.3);
        assert!(vals[3].is_nan());
        assert!(vals[99].is_nan());
        for v in &vals[100..150] {
            assert_eq!(*v, 1.4);
        }
        assert_eq!(vals[150], 1.5);
    }

    #[test]
    fn intervals_match_pybigwig_ground_truth() {
        let mut bw = BigWigReader::open(&fixture_path()).unwrap();
        let ivs = bw.intervals("10", 0, 400).unwrap();
        assert_eq!(ivs, vec![(200, 300, 2.0)]);

        let ivs1 = bw.intervals("1", 0, 4).unwrap();
        assert_eq!(ivs1, vec![(0, 1, 0.1), (1, 2, 0.2), (2, 3, 0.3)]);
    }

    #[test]
    fn open_nonexistent_file_is_an_io_error() {
        assert!(BigWigReader::open(std::path::Path::new("/nonexistent/path.bw")).is_err());
    }
}
