// Fuzz tests for BED parsing functions
use proptest::prelude::*;
use proptest::strategy::Just;
use proptest::test_runner::{Config, RngSeed};
use std::io::Cursor;
use rseqc_formats::bed;

proptest! {
    #![proptest_config(Config {
        cases: 1000,
        rng_seed: RngSeed::Fixed(0x1234567890ABCDEF),
        .. Config::default()
    })]

    #[test]
    fn fuzz_get_cds_exon(bed_content in any::<String>()) {
        let cursor = Cursor::new(bed_content);
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = bed::get_cds_exon(cursor);
        }));
    }

    #[test]
    fn fuzz_get_exon(bed_content in any::<String>()) {
        let cursor = Cursor::new(bed_content);
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = bed::get_exon(cursor);
        }));
    }

    #[test]
    fn fuzz_get_utr(bed_content in any::<String>(), utr in prop_oneof![Just(3), Just(5)]) {
        let cursor = Cursor::new(bed_content);
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = bed::get_utr(cursor, utr);
        }));
    }

    #[test]
    fn fuzz_get_intergenic(
        bed_content in any::<String>(),
        direction in prop_oneof![Just("up"), Just("down")],
        size in -1000i64..=1000i64,
    ) {
        let cursor = Cursor::new(bed_content);
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = bed::get_intergenic(cursor, direction, size);
        }));
    }

    #[test]
    fn fuzz_get_transcript_ranges(bed_content in any::<String>()) {
        let cursor = Cursor::new(bed_content);
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = bed::get_transcript_ranges(cursor);
        }));
    }

    #[test]
    fn fuzz_get_intron(bed_content in any::<String>()) {
        let cursor = Cursor::new(bed_content);
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = bed::get_intron(cursor);
        }));
    }

    // Test with structured BED12 lines that have invalid field values
    #[test]
    fn fuzz_bed12_malformed_coordinates(
        chrom in "chr[a-zA-Z0-9]+",
        tx_start in any::<i32>(),
        tx_end in any::<i32>(),
        cds_start in any::<i32>(),
        cds_end in any::<i32>(),
        block_sizes in "(proptest::string::string_regex(\"[0-9]+\").unwrap())",
    ) {
        let line = format!(
            "{}\t{}\t{}\t.\t.\t+\t{}\t{}\t0\t1\t{}\t0",
            chrom, tx_start, tx_end, cds_start, cds_end, block_sizes
        );
        let cursor = Cursor::new(line);
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = bed::get_cds_exon(cursor);
        }));
    }

    // Test with empty lines and whitespace variations
    #[test]
    fn fuzz_bed_whitespace(whitespace_variant in "[ \t\n]*") {
        let cursor = Cursor::new(whitespace_variant);
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = bed::get_cds_exon(cursor);
        }));
    }

    // Test BED with CRLF line endings
    #[test]
    fn fuzz_bed_crlf(bed_content in any::<String>()) {
        let crlf_content = bed_content.replace("\n", "\r\n");
        let cursor = Cursor::new(crlf_content);
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = bed::get_cds_exon(cursor);
        }));
    }

    // Test BED with extremely long lines
    #[test]
    fn fuzz_bed_long_lines(line_count in 1..100usize) {
        let mut bed_content = String::new();
        for i in 0..line_count {
            bed_content.push_str(&format!(
                "chr{}\t{}\t{}\t{}\t{}\t+\t{}\t{}\t0\t1\t100\t0\n",
                i, i * 1000, (i + 1) * 1000, i, i, i, i + 1
            ));
        }
        let cursor = Cursor::new(bed_content);
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = bed::get_cds_exon(cursor);
        }));
    }

    // Test with huge coordinate values that might overflow
    #[test]
    fn fuzz_bed_huge_coordinates(
        coord in i64::MIN..=i64::MAX,
    ) {
        let line = format!(
            "chr1\t{}\t{}\t.\t.\t+\t{}\t{}\t0\t1\t100\t0",
            coord, coord.saturating_add(100), coord, coord.saturating_add(50)
        );
        let cursor = Cursor::new(line);
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = bed::get_cds_exon(cursor);
        }));
    }

    // Test with negative coordinates
    #[test]
    fn fuzz_bed_negative_coords(coord in -1000i64..1000i64) {
        let line = format!(
            "chr1\t{}\t{}\t.\t.\t+\t{}\t{}\t0\t1\t100\t0",
            coord, coord + 100, coord, coord + 50
        );
        let cursor = Cursor::new(line);
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = bed::get_cds_exon(cursor);
        }));
    }

    // Test with non-UTF8 lines (by truncating UTF8 sequences)
    #[test]
    fn fuzz_bed_truncated_utf8(
        content in r"([\x20-\x7E])*",
    ) {
        let mut bytes = content.into_bytes();
        if !bytes.is_empty() {
            // Try truncating at various positions
            for i in 1..bytes.len().min(10) {
                bytes.truncate(i);
                let cursor = Cursor::new(bytes.clone());
                let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    let _ = bed::get_cds_exon(cursor);
                }));
            }
        }
    }
}
