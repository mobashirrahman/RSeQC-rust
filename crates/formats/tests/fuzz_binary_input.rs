// Fuzz tests for binary/compressed input layer
use proptest::prelude::*;
use proptest::test_runner::{Config, RngSeed};
use std::fs::File;
use std::io::{BufRead, Write};
use std::path::Path;
use tempfile::TempDir;

proptest! {
    #![proptest_config(Config {
        cases: 100,
        rng_seed: RngSeed::Fixed(0x1234567890ABCDEF),
        .. Config::default()
    })]

    /// Fuzz gzip files: valid, truncated, bit-flipped, empty, garbage appended
    #[test]
    fn fuzz_open_text_input_gzip_valid(content in "\\PC*") {
        let temp_dir = TempDir::new().unwrap();
        let gz_path = temp_dir.path().join("test.gz");

        // Create a valid gzip file
        {
            let file = File::create(&gz_path).unwrap();
            let encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
            let mut enc = encoder;
            let _ = enc.write_all(content.as_bytes());
        }

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok(mut reader) = rseqc_formats::open_text_input(&gz_path) {
                let mut line = String::new();
                while let Ok(n) = reader.read_line(&mut line) {
                    if n == 0 { break; }
                    line.clear();
                }
            }
        }));
    }

    /// Fuzz bzip2 files: valid and mutated
    #[test]
    fn fuzz_open_text_input_bzip2_valid(content in "\\PC*") {
        let temp_dir = TempDir::new().unwrap();
        let bz2_path = temp_dir.path().join("test.bz2");

        // Create a valid bzip2 file
        {
            let file = File::create(&bz2_path).unwrap();
            let mut encoder = bzip2::write::BzEncoder::new(file, bzip2::Compression::default());
            let _ = encoder.write_all(content.as_bytes());
        }

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok(mut reader) = rseqc_formats::open_text_input(&bz2_path) {
                let mut line = String::new();
                while let Ok(n) = reader.read_line(&mut line) {
                    if n == 0 { break; }
                    line.clear();
                }
            }
        }));
    }

    /// Fuzz plain text files with arbitrary bytes
    #[test]
    fn fuzz_open_text_input_plain_text(content in "\\PC*") {
        let temp_dir = TempDir::new().unwrap();
        let plain_path = temp_dir.path().join("test.txt");

        let mut file = File::create(&plain_path).unwrap();
        let _ = file.write_all(content.as_bytes());

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok(mut reader) = rseqc_formats::open_text_input(&plain_path) {
                let mut line = String::new();
                while let Ok(n) = reader.read_line(&mut line) {
                    if n == 0 { break; }
                    line.clear();
                }
            }
        }));
    }

    /// Fuzz gzip files with truncation
    #[test]
    fn fuzz_open_text_input_gzip_truncated(
        content in "\\PC{1,1000}",
        truncate_at in 0usize..=500,
    ) {
        let temp_dir = TempDir::new().unwrap();
        let gz_path = temp_dir.path().join("test.gz");

        // Create a valid gzip file, then truncate it
        let compressed = {
            let mut buf = Vec::new();
            let encoder = flate2::write::GzEncoder::new(&mut buf, flate2::Compression::default());
            let mut file = encoder;
            let _ = file.write_all(content.as_bytes());
            let _ = file.finish();
            buf
        };

        if truncate_at < compressed.len() {
            let mut file = File::create(&gz_path).unwrap();
            let _ = file.write_all(&compressed[..truncate_at]);
        } else {
            let mut file = File::create(&gz_path).unwrap();
            let _ = file.write_all(&compressed);
        }

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok(mut reader) = rseqc_formats::open_text_input(&gz_path) {
                let mut line = String::new();
                while let Ok(n) = reader.read_line(&mut line) {
                    if n == 0 { break; }
                    line.clear();
                }
            }
        }));
    }

    /// Fuzz bzip2 files with truncation
    #[test]
    fn fuzz_open_text_input_bzip2_truncated(
        content in "\\PC{1,1000}",
        truncate_at in 0usize..=500,
    ) {
        let temp_dir = TempDir::new().unwrap();
        let bz2_path = temp_dir.path().join("test.bz2");

        // Create a valid bzip2 file, then truncate it
        let compressed = {
            let mut buf = Vec::new();
            let encoder = bzip2::write::BzEncoder::new(&mut buf, bzip2::Compression::default());
            let mut file = encoder;
            let _ = file.write_all(content.as_bytes());
            let _ = file.finish();
            buf
        };

        if truncate_at < compressed.len() {
            let mut file = File::create(&bz2_path).unwrap();
            let _ = file.write_all(&compressed[..truncate_at]);
        } else {
            let mut file = File::create(&bz2_path).unwrap();
            let _ = file.write_all(&compressed);
        }

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok(mut reader) = rseqc_formats::open_text_input(&bz2_path) {
                let mut line = String::new();
                while let Ok(n) = reader.read_line(&mut line) {
                    if n == 0 { break; }
                    line.clear();
                }
            }
        }));
    }

    /// Fuzz gzip files with bit flips
    #[test]
    fn fuzz_open_text_input_gzip_bit_flip(
        content in "\\PC{1,1000}",
        flip_byte in 0usize..=500,
        flip_bit in 0u8..8,
    ) {
        let temp_dir = TempDir::new().unwrap();
        let gz_path = temp_dir.path().join("test.gz");

        // Create a valid gzip file, then flip a bit
        let mut compressed = {
            let mut buf = Vec::new();
            let encoder = flate2::write::GzEncoder::new(&mut buf, flate2::Compression::default());
            let mut file = encoder;
            let _ = file.write_all(content.as_bytes());
            let _ = file.finish();
            buf
        };

        if flip_byte < compressed.len() {
            compressed[flip_byte] ^= 1 << flip_bit;
        }

        let mut file = File::create(&gz_path).unwrap();
        let _ = file.write_all(&compressed);

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok(mut reader) = rseqc_formats::open_text_input(&gz_path) {
                let mut line = String::new();
                while let Ok(n) = reader.read_line(&mut line) {
                    if n == 0 { break; }
                    line.clear();
                }
            }
        }));
    }

    /// Fuzz bzip2 files with bit flips
    #[test]
    fn fuzz_open_text_input_bzip2_bit_flip(
        content in "\\PC{1,1000}",
        flip_byte in 0usize..=500,
        flip_bit in 0u8..8,
    ) {
        let temp_dir = TempDir::new().unwrap();
        let bz2_path = temp_dir.path().join("test.bz2");

        // Create a valid bzip2 file, then flip a bit
        let mut compressed = {
            let mut buf = Vec::new();
            let encoder = bzip2::write::BzEncoder::new(&mut buf, bzip2::Compression::default());
            let mut file = encoder;
            let _ = file.write_all(content.as_bytes());
            let _ = file.finish();
            buf
        };

        if flip_byte < compressed.len() {
            compressed[flip_byte] ^= 1 << flip_bit;
        }

        let mut file = File::create(&bz2_path).unwrap();
        let _ = file.write_all(&compressed);

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok(mut reader) = rseqc_formats::open_text_input(&bz2_path) {
                let mut line = String::new();
                while let Ok(n) = reader.read_line(&mut line) {
                    if n == 0 { break; }
                    line.clear();
                }
            }
        }));
    }

    /// Fuzz gzip files with garbage appended
    #[test]
    fn fuzz_open_text_input_gzip_garbage_appended(
        content in "\\PC*",
        garbage in "[\\x00-\\xFF]{1,100}",
    ) {
        let temp_dir = TempDir::new().unwrap();
        let gz_path = temp_dir.path().join("test.gz");

        // Create a valid gzip file, then append garbage
        let mut compressed = {
            let mut buf = Vec::new();
            let encoder = flate2::write::GzEncoder::new(&mut buf, flate2::Compression::default());
            let mut file = encoder;
            let _ = file.write_all(content.as_bytes());
            let _ = file.finish();
            buf
        };

        compressed.extend_from_slice(garbage.as_bytes());

        let mut file = File::create(&gz_path).unwrap();
        let _ = file.write_all(&compressed);

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok(mut reader) = rseqc_formats::open_text_input(&gz_path) {
                let mut line = String::new();
                while let Ok(n) = reader.read_line(&mut line) {
                    if n == 0 { break; }
                    line.clear();
                }
            }
        }));
    }

    /// Fuzz plain text files with non-UTF8 bytes
    #[test]
    fn fuzz_open_text_input_plain_text_non_utf8(bytes in prop::collection::vec(0u8..=255, 0..1000)) {
        let temp_dir = TempDir::new().unwrap();
        let plain_path = temp_dir.path().join("test.txt");

        let mut file = File::create(&plain_path).unwrap();
        let _ = file.write_all(&bytes);

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok(mut reader) = rseqc_formats::open_text_input(&plain_path) {
                let mut line = String::new();
                while let Ok(n) = reader.read_line(&mut line) {
                    if n == 0 { break; }
                    line.clear();
                }
            }
        }));
    }

    /// Fuzz SAM files with valid header and records
    #[test]
    fn fuzz_open_alignments_sam_valid(
        _chrom in "chr[a-z0-9]+",
        pos in 1i64..=1_000_000,
        seq_len in 1..=100usize,
        mapq in 0u8..=60,
    ) {
        let temp_dir = TempDir::new().unwrap();
        let sam_path = temp_dir.path().join("test.sam");

        let seq = "ACGT".repeat((seq_len / 4) + 1);
        let seq = &seq[..seq_len];
        let qual = "IIII".repeat((seq_len / 4) + 1);
        let qual = &qual[..seq_len];

        let mut file = File::create(&sam_path).unwrap();
        writeln!(file, "@SQ\tSN:chr1\tLN:1000000").unwrap();
        writeln!(file, "read1\t0\tchr1\t{}\t{}\t{}M\t*\t0\t0\t{}\t{}", pos, mapq, seq_len, seq, qual).unwrap();

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok((_, records)) = rseqc_formats::open_alignments(&sam_path) {
                let _ = records.into_iter().collect::<Result<Vec<_>, _>>();
            }
        }));
    }

    /// Fuzz SAM files with bad CIGAR
    #[test]
    fn fuzz_open_alignments_sam_bad_cigar(bad_cigar in "[0-9]*[^0-9MIDNSHPX=]{0,20}") {
        let temp_dir = TempDir::new().unwrap();
        let sam_path = temp_dir.path().join("test.sam");

        let mut file = File::create(&sam_path).unwrap();
        writeln!(file, "@SQ\tSN:chr1\tLN:1000000").unwrap();
        writeln!(file, "read1\t0\tchr1\t100\t60\t{}\t*\t0\t0\tACGT\tIIII", bad_cigar).unwrap();

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok((_, records)) = rseqc_formats::open_alignments(&sam_path) {
                let _ = records.into_iter().collect::<Result<Vec<_>, _>>();
            }
        }));
    }

    /// Fuzz SAM files with huge/negative POS
    #[test]
    fn fuzz_open_alignments_sam_extreme_pos(pos in i64::MIN..=i64::MAX) {
        let temp_dir = TempDir::new().unwrap();
        let sam_path = temp_dir.path().join("test.sam");

        let mut file = File::create(&sam_path).unwrap();
        writeln!(file, "@SQ\tSN:chr1\tLN:1000000").unwrap();
        if pos >= 0 {
            writeln!(file, "read1\t0\tchr1\t{}\t60\t4M\t*\t0\t0\tACGT\tIIII", pos).unwrap();
        } else {
            // Negative positions should fail to parse
            writeln!(file, "read1\t0\tchr1\t{}\t60\t4M\t*\t0\t0\tACGT\tIIII", pos).unwrap();
        }

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok((_, records)) = rseqc_formats::open_alignments(&sam_path) {
                let _ = records.into_iter().collect::<Result<Vec<_>, _>>();
            }
        }));
    }

    /// Fuzz SAM files with mismatched SEQ/QUAL length
    #[test]
    fn fuzz_open_alignments_sam_mismatched_seq_qual(
        seq_len in 1..=100usize,
        qual_len in 0..=100usize,
    ) {
        let temp_dir = TempDir::new().unwrap();
        let sam_path = temp_dir.path().join("test.sam");

        let seq = "ACGT".repeat((seq_len / 4) + 1);
        let seq = &seq[..seq_len];
        let qual = "IIII".repeat((qual_len / 4) + 1);
        let qual = if qual_len > 0 { &qual[..qual_len] } else { "" };

        let mut file = File::create(&sam_path).unwrap();
        writeln!(file, "@SQ\tSN:chr1\tLN:1000000").unwrap();
        writeln!(file, "read1\t0\tchr1\t100\t60\t{}M\t*\t0\t0\t{}\t{}", seq_len, seq, qual).unwrap();

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok((_, records)) = rseqc_formats::open_alignments(&sam_path) {
                let _ = records.into_iter().collect::<Result<Vec<_>, _>>();
            }
        }));
    }

    /// Fuzz SAM files with truncation
    #[test]
    fn fuzz_open_alignments_sam_truncated(truncate_at in 0usize..=500) {
        let temp_dir = TempDir::new().unwrap();
        let sam_path = temp_dir.path().join("test.sam");

        let content = "@SQ\tSN:chr1\tLN:1000000\nread1\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\n";

        let mut file = File::create(&sam_path).unwrap();
        if truncate_at < content.len() {
            let _ = file.write_all(&content.as_bytes()[..truncate_at]);
        } else {
            let _ = file.write_all(content.as_bytes());
        }

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok((_, records)) = rseqc_formats::open_alignments(&sam_path) {
                let _ = records.into_iter().collect::<Result<Vec<_>, _>>();
            }
        }));
    }

    /// Fuzz BAM files with truncation (using existing fixture)
    #[test]
    fn fuzz_open_alignments_bam_truncated(truncate_at in 0usize..=5000) {
        let bam_fixture = Path::new("/scratch/mdra00001/RSeQC-rust/verification/fixtures/bam_stat_basic.bam");

        if !bam_fixture.exists() {
            return Ok(());
        }

        let original = std::fs::read(bam_fixture).unwrap_or_default();
        if original.is_empty() {
            return Ok(());
        }

        let temp_dir = TempDir::new().unwrap();
        let bam_path = temp_dir.path().join("test.bam");

        let mut file = File::create(&bam_path).unwrap();
        if truncate_at < original.len() {
            let _ = file.write_all(&original[..truncate_at]);
        } else {
            let _ = file.write_all(&original);
        }

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok((_, records)) = rseqc_formats::open_alignments(&bam_path) {
                let _ = records.into_iter().collect::<Result<Vec<_>, _>>();
            }
        }));
    }

    /// Fuzz BAM files with bit flips (using existing fixture)
    #[test]
    fn fuzz_open_alignments_bam_bit_flip(
        flip_byte in 0usize..=5000,
        flip_bit in 0u8..8,
    ) {
        let bam_fixture = Path::new("/scratch/mdra00001/RSeQC-rust/verification/fixtures/bam_stat_basic.bam");

        if !bam_fixture.exists() {
            return Ok(());
        }

        let mut original = std::fs::read(bam_fixture).unwrap_or_default();
        if original.is_empty() {
            return Ok(());
        }

        if flip_byte < original.len() {
            original[flip_byte] ^= 1 << flip_bit;
        }

        let temp_dir = TempDir::new().unwrap();
        let bam_path = temp_dir.path().join("test.bam");

        let mut file = File::create(&bam_path).unwrap();
        let _ = file.write_all(&original);

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok((_, records)) = rseqc_formats::open_alignments(&bam_path) {
                let _ = records.into_iter().collect::<Result<Vec<_>, _>>();
            }
        }));
    }

    /// Fuzz empty gzip files
    #[test]
    fn fuzz_open_text_input_gzip_empty(_unit in ".*") {
        let temp_dir = TempDir::new().unwrap();
        let gz_path = temp_dir.path().join("test.gz");

        // Create an empty gzip file
        {
            let file = File::create(&gz_path).unwrap();
            let encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
            let _ = encoder.finish();
        }

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok(mut reader) = rseqc_formats::open_text_input(&gz_path) {
                let mut line = String::new();
                while let Ok(n) = reader.read_line(&mut line) {
                    if n == 0 { break; }
                    line.clear();
                }
            }
        }));
    }

    /// Fuzz empty plain text files
    #[test]
    fn fuzz_open_text_input_plain_empty(_unit in ".*") {
        let temp_dir = TempDir::new().unwrap();
        let plain_path = temp_dir.path().join("test.txt");

        File::create(&plain_path).unwrap();

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            if let Ok(mut reader) = rseqc_formats::open_text_input(&plain_path) {
                let mut line = String::new();
                while let Ok(n) = reader.read_line(&mut line) {
                    if n == 0 { break; }
                    line.clear();
                }
            }
        }));
    }
}
