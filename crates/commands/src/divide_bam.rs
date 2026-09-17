//! Port of `divide_bam.py`: randomly divide a BAM file into approximately equal subsets.
//! Contract: see `compatibility/commands.yaml` entry `divide_bam.py`; algorithm
//! ported from `divide_bam()` in `oracle/upstream-src/scripts/divide_bam.py`
//! (lines 110-174).
//!
//! `.bai` index generation (DIV-0006, closed) is implemented at the CLI
//! layer via `rseqc_formats::write_bai_index`; this module's own
//! `divide_bam` function is unaffected -- it only ever writes the plain
//! BAM records.

use std::io;

use noodles_bam as bam;
use noodles_sam::alignment::io::Write as _;


#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct DivideCounts {
    pub retained: u64,
    pub skipped: u64,
    pub per_file: Vec<u64>,
}

/// Mirrors upstream `main()`'s final report:
/// `print(f"{path}\t{count}")` (stdout, one line per output file) followed
/// by `print(f"Total alignments written: {n}", file=sys.stderr)` and,
/// only when `--skip-unmap` was given, `print(f"Unmapped alignments
/// skipped: {n}", file=sys.stderr)`.
pub fn render_report(paths: &[String], counts: &DivideCounts, skip_unmap: bool) -> (String, String) {
    let mut stdout = String::new();
    for (path, count) in paths.iter().zip(counts.per_file.iter()) {
        stdout.push_str(&format!("{path}\t{count}\n"));
    }

    let mut stderr = String::new();
    stderr.push_str(&format!("Total alignments written: {}\n", counts.retained));
    if skip_unmap {
        stderr.push_str(&format!("Unmapped alignments skipped: {}\n", counts.skipped));
    }

    (stdout, stderr)
}

/// Divides `records` (in file order) across multiple BAM writers using the same
/// header/template as the input. Pure logic over generic writers so it's
/// testable with in-memory sinks; no file handling here.
pub fn divide_bam<I, W>(
    records: I,
    header: &noodles_sam::Header,
    outputs: &mut [bam::io::Writer<W>],
    skip_unmapped: bool,
    rng: &mut impl rand::Rng,
) -> io::Result<DivideCounts>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
    W: io::Write,
{
    let mut counts = DivideCounts {
        per_file: vec![0; outputs.len()],
        ..Default::default()
    };

    let mut query_assignments: std::collections::HashMap<Vec<u8>, usize> = std::collections::HashMap::new();

    for result in records {
        let record = result?;
        
        if skip_unmapped && record.flags().is_unmapped() {
            counts.skipped += 1;
            continue;
        }

        let query_name = record.name().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "Alignment record missing query name",
            )
        })?;

        let output_index = query_assignments.entry(query_name.to_string().into_bytes()).or_insert_with(|| {
            rng.gen_range(0..outputs.len())
        });

        outputs[*output_index].write_alignment_record(header, &record)?;
        counts.per_file[*output_index] += 1;
        counts.retained += 1;
    }

    Ok(counts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::{
        self as sam,
        alignment::record::{Flags, MappingQuality, cigar::Op, cigar::op::Kind},
        alignment::record_buf::{Cigar, RecordBuf},
        header::record::value::{Map, map::ReferenceSequence},
    };
    use std::num::NonZeroUsize;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn test_header() -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
            )
            .build()
    }

    fn to_bam_records(header: &sam::Header, records: &[RecordBuf]) -> Vec<bam::Record> {
        let mut buf = Vec::new();
        {
            let mut writer = bam::io::Writer::new(&mut buf);
            writer.write_header(header).unwrap();
            for r in records {
                writer.write_alignment_record(header, r).unwrap();
            }
        }
        let mut reader = bam::io::Reader::new(buf.as_slice());
        reader.read_header().unwrap();
        reader.records().map(|r| r.unwrap()).collect()
    }

    #[test]
    fn divides_records_with_fixed_seed() {
        let header = test_header();

        // Create records with the same query name (paired-end mates)
        let mate1 = RecordBuf::builder()
            .set_name("read1")
            .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .build();

        let mate2 = RecordBuf::builder()
            .set_name("read1")
            .set_flags(Flags::SEGMENTED | Flags::LAST_SEGMENT)
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .build();

        // Create records with different query names
        let read2 = RecordBuf::builder()
            .set_name("read2")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .build();

        let read3 = RecordBuf::builder()
            .set_name("read3")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .build();

        let input = vec![mate1, mate2, read2, read3];
        let bam_records = to_bam_records(&header, &input);

        // Create output buffers
        let mut out1_buf = Vec::new();
        let mut out2_buf = Vec::new();
        let mut out1_writer = bam::io::Writer::new(&mut out1_buf);
        let mut out2_writer = bam::io::Writer::new(&mut out2_buf);
        out1_writer.write_header(&header).unwrap();
        out2_writer.write_header(&header).unwrap();

        let mut outputs = vec![out1_writer, out2_writer];

        // Use a fixed seed for reproducible results
        let mut rng = StdRng::seed_from_u64(42);
        let counts = divide_bam(
            bam_records.into_iter().map(Ok),
            &header,
            &mut outputs,
            false,
            &mut rng,
        ).unwrap();

        // Both mates with the same query name should go to the same file.
        // Drop the writers (releasing their mutable borrows on the output
        // buffers) before reading the buffers back.
        drop(outputs);

        let mut r1_reader = bam::io::Reader::new(out1_buf.as_slice());
        r1_reader.read_header().unwrap();
        let r1_records: Vec<_> = r1_reader.records().map(|r| r.unwrap()).collect();

        let mut r2_reader = bam::io::Reader::new(out2_buf.as_slice());
        r2_reader.read_header().unwrap();
        let r2_records: Vec<_> = r2_reader.records().map(|r| r.unwrap()).collect();

        // Check that we have 4 records total, with some distribution
        assert_eq!(counts.retained, 4);
        assert_eq!(counts.skipped, 0);
        assert_eq!(counts.per_file.len(), 2);
        assert_eq!(counts.per_file.iter().sum::<u64>(), 4);

        // Check that both mates from read1 are in the same output
        let read1_records_in_out1 = r1_records
            .iter()
            .filter(|r| r.name().map(|n| n.to_string()) == Some("read1".to_string()))
            .count();
        let read1_records_in_out2 = r2_records
            .iter()
            .filter(|r| r.name().map(|n| n.to_string()) == Some("read1".to_string()))
            .count();

        // Both mates should be in the same output file. (Note: read2/read3
        // may legitimately also land in whichever file holds read1's mates
        // — 3 distinct query names split across 2 outputs makes that likely
        // — so we only assert on read1's own co-location, not file purity.)
        assert_eq!(read1_records_in_out1 + read1_records_in_out2, 2);
        assert!(read1_records_in_out1 == 2 || read1_records_in_out2 == 2);
    }

    #[test]
    fn skips_unmapped_records() {
        let header = test_header();

        let mapped = RecordBuf::builder()
            .set_name("read1")
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .build();

        let unmapped = RecordBuf::builder()
            .set_name("read2")
            .set_flags(Flags::UNMAPPED)
            .build();

        let input = vec![mapped, unmapped];
        let bam_records = to_bam_records(&header, &input);

        let mut out_buf = Vec::new();
        let mut out_writer = bam::io::Writer::new(&mut out_buf);
        out_writer.write_header(&header).unwrap();

        let mut outputs = vec
![out_writer];

        let mut rng = StdRng::seed_from_u64(42);
        let counts = divide_bam(
            bam_records.into_iter().map(Ok),
            &header,
            &mut outputs,
            true, // skip unmapped
            &mut rng,
        ).unwrap();

        assert_eq!(counts.retained, 1);
        assert_eq!(counts.skipped, 1);
        assert_eq!(counts.per_file[0], 1);
    }

    #[test]
    fn handles_missing_query_name_error() {
        let header = test_header();

        // Create a record without a query name
        let record = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 30)]))
            .build();

        let input = vec![record];
        let bam_records = to_bam_records(&header, &input);

        let mut out_buf = Vec::new();
        let mut out_writer = bam::io::Writer::new(&mut out_buf);
        out_writer.write_header(&header).unwrap();

        let mut outputs = vec
![out_writer];

        let mut rng = StdRng::seed_from_u64(42);
        
        // This should return an error because the record has no query name
        let result = divide_bam(
            bam_records.into_iter().map(Ok),
            &header,
            &mut outputs,
            false,
            &mut rng,
        );

        assert!(result.is_err());
        if let Err(ref e) = result {
            assert_eq!(e.kind(), io::ErrorKind::InvalidData);
        }
    }

    #[test]
    fn render_report_matches_upstream_format() {
        let counts = DivideCounts {
            retained: 7,
            skipped: 2,
            per_file: vec![3, 4],
        };
        let paths = vec!["out_0.bam".to_string(), "out_1.bam".to_string()];

        let (stdout, stderr) = render_report(&paths, &counts, true);
        assert_eq!(stdout, "out_0.bam\t3\nout_1.bam\t4\n");
        assert_eq!(stderr, "Total alignments written: 7\nUnmapped alignments skipped: 2\n");

        let (_, stderr_no_skip) = render_report(&paths, &counts, false);
        assert_eq!(stderr_no_skip, "Total alignments written: 7\n");
    }
}