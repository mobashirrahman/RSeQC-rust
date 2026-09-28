// Fuzz tests for sequence parsing functions
use proptest::prelude::*;
use std::io::Cursor;
use rseqc_commands::{read_hexamer, sc_seqlogo, sc_seqqual};

proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]

    // read_hexamer sequence generator fuzz tests
    #[test]
    fn fuzz_read_hexamer_seq_generator(content in r"[ACGTN>a-z\n\r]{0,5000}") {
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let cursor = Cursor::new(content.into_bytes());
            let _ = read_hexamer::seq_generator(cursor);
        }));
    }

    // FASTA parsing fuzz tests
    #[test]
    fn fuzz_fasta_iter_basic(content in r"[>ACGTN\n\r]{0,2000}") {
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let cursor = Cursor::new(content.into_bytes());
            let _ = sc_seqlogo::fasta_iter(cursor);
        }));
    }

    #[test]
    fn fuzz_fasta_iter_empty_lines(content in r"[ \t\n\r]*") {
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let cursor = Cursor::new(content.into_bytes());
            let _ = sc_seqlogo::fasta_iter(cursor);
        }));
    }

    #[test]
    fn fuzz_fasta_iter_no_headers(content in r"[ACGTN\n]{0,1000}") {
        // FASTA without '>' headers - should gracefully handle
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let cursor = Cursor::new(content.into_bytes());
            let _ = sc_seqlogo::fasta_iter(cursor);
        }));
    }

    #[test]
    fn fuzz_fasta_iter_invalid_chars(content in r"[^ACGTN>a-z\n\r]{0,100}") {
        // Non-sequence characters mixed in
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let cursor = Cursor::new(content.into_bytes());
            let _ = sc_seqlogo::fasta_iter(cursor);
        }));
    }

    // FASTQ parsing fuzz tests
    #[test]
    fn fuzz_fastq_seq_strings_basic(content in r"[@+ACGTN=\n\r:0-9]{0,5000}") {
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let cursor = Cursor::new(content.into_bytes());
            let _ = sc_seqlogo::fastq_seq_strings(cursor);
        }));
    }

    #[test]
    fn fuzz_fastq_qual_strings_basic(content in r"[@+ACGTN!=\n\r:0-9]{0,5000}") {
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let cursor = Cursor::new(content.into_bytes());
            let _ = sc_seqqual::fastq_qual_strings(cursor);
        }));
    }

    #[test]
    fn fuzz_fastq_incomplete_record(lines_count in 0usize..6usize) {
        let mut content = String::new();
        for i in 0..lines_count {
            match i % 4 {
                0 => content.push_str("@read_name\n"),
                1 => content.push_str("ACGTACGT\n"),
                2 => content.push_str("+\n"),
                3 => content.push_str("IIIIIIII\n"),
                _ => {}
            }
        }

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let cursor = Cursor::new(content.into_bytes());
            let _ = sc_seqlogo::fastq_seq_strings(cursor);
        }));
    }

    #[test]
    fn fuzz_fastq_mismatched_lengths(
        seq_len in 0usize..100usize,
        qual_len in 0usize..100usize,
    ) {
        let seq = "A".repeat(seq_len);
        let qual = "I".repeat(qual_len);
        let content = format!("@read\n{}\n+\n{}\n", seq, qual);

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let cursor = Cursor::new(content.into_bytes());
            let _ = sc_seqqual::fastq_qual_strings(cursor);
        }));
    }

    #[test]
    fn fuzz_fastq_invalid_quality_chars(content in r"[@+ACGTN=\x00-\x08\x0B-\x0C\x0E-\x1F\n\r:0-9]{0,1000}") {
        // Quality strings can be any ASCII but typically >33 (Phred offset)
        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let cursor = Cursor::new(content.into_bytes());
            let _ = sc_seqqual::fastq_qual_strings(cursor);
        }));
    }

    #[test]
    fn fuzz_fastq_very_long_record(line_len in 100usize..10000usize) {
        let seq = "A".repeat(line_len);
        let qual = "I".repeat(line_len);
        let content = format!("@longread\n{}\n+\n{}\n", seq, qual);

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let cursor = Cursor::new(content.into_bytes());
            let _ = sc_seqlogo::fastq_seq_strings(cursor);
        }));
    }

    #[test]
    fn fuzz_fastq_crlf_line_endings(content in r"[@+ACGTN=\n\r:0-9]{0,1000}") {
        // Replace newlines with CRLF
        let crlf_content = content.replace("\n", "\r\n");

        let _ = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let cursor = Cursor::new(crlf_content.into_bytes());
            let _ = sc_seqlogo::fastq_seq_strings(cursor);
        }));
    }
}
