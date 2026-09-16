//! Port of `read_NVC.py`: calculate nucleotide frequency at each read cycle.
//! Contract: see `compatibility/commands.yaml` entry `read_NVC.py`; algorithm
//! ported from `readsNVC()` in `oracle/upstream-src/src/qcmodule/SAM.py`
//! (lines 2954-3009).
//!
//! Notable upstream behavior preserved exactly: ONLY mapq filtering is applied
//! (no QC-fail/duplicate/secondary/unmapped filtering); reverse-strand reads
//! are reverse-complemented; the output table's row count equals the last
//! processed read's length (not the max length seen across all records).
//!
//! Known gap: `-x/--nx` and plot generation are not implemented here yet —
//! disclosed follow-up. BAM-only for now (see DIV-0002 pattern).

use std::collections::HashMap;
use std::io;

use noodles_bam as bam;

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct NvcTable {
    pub rows: Vec<NvcRow>,
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct NvcRow {
    pub position: usize,
    pub a: u64,
    pub c: u64,
    pub g: u64,
    pub t: u64,
    pub n: u64,
    pub x: u64,
}

fn complement(base: u8) -> u8 {
    // Reuse the same complement table as bam2fq.py
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' => b'N',
        b'X' => b'X',
        other => other,
    }
}

/// Computes nucleotide composition table for a sequence of BAM records,
/// matching the branch order and category semantics of upstream's `readsNVC()`.
/// Pure computation: no I/O, no printing.
pub fn compute_nvc<I>(records: I, q_cut: u8) -> io::Result<NvcTable>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut base_freq = HashMap::new();
    let mut last_read_length = 0;

    for result in records {
        let record = result?;
        
        // NOTE: is_unmapped and is_qcfail checks are commented out in
        // upstream (dead code) -- do NOT filter on them. Only the mapq
        // cutoff applies.
        
        // pysam's raw `.mapq` reports a missing MAPQ as byte value 255;
        // noodles represents "missing" as `None`. 255 always compares >=
        // any realistic q_cut, so treating `None` as 255 here reproduces
        // upstream's comparison behavior without a special case.
        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        
        if mapq < q_cut {
            continue;
        }

        let sequence = record.sequence();
        last_read_length = sequence.len();
        
        let mut seq: Vec<u8> = sequence.iter().collect();
        
        if record.flags().is_reverse_complemented() {
            seq.reverse();
            for b in seq.iter_mut() {
                *b = complement(*b);
            }
        }
        
        for (i, base) in seq.iter().enumerate() {
            // Only count ACGTNX bases as upstream does
            match base {
                b'A' | b'C' | b'G' | b'T' | b'N' | b'X' => {
                    *base_freq.entry((i, *base)).or_insert(0) += 1;
                }
                _ => {
                    // Other bases are not counted in the output, matching upstream
                }
            }
        }
    }

    // Build the table with rows equal to the last read's length
    let mut rows = Vec::with_capacity(last_read_length);
    for position in 0..last_read_length {
        let row = NvcRow {
            position,
            a: *base_freq.get(&(position, b'A')).unwrap_or(&0),
            c: *base_freq.get(&(position, b'C')).unwrap_or(&0),
            g: *base_freq.get(&(position, b'G')).unwrap_or(&0),
            t: *base_freq.get(&(position, b'T')).unwrap_or(&0),
            n: *base_freq.get(&(position, b'N')).unwrap_or(&0),
            x: *base_freq.get(&(position, b'X')).unwrap_or(&0),
        };
        rows.push(row);
    }

    Ok(NvcTable { rows })
}

/// Renders the NVC table as tab-separated text, matching upstream's output format
pub fn render_nvc_table(table: &NvcTable) -> String {
    let mut lines = Vec::new();
    
    // Header row
    lines.push("Position\tA\tC\tG\tT\tN\tX".to_string());
    
    // Data rows
    for row in &table.rows {
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.position, row.a, row.c, row.g, row.t, row.n, row.x
        ));
    }
    
    lines.join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::{
        self as sam,
        alignment::{
            io::Write as _,
            record::{Flags, MappingQuality, cigar::Op, cigar::op::Kind as CigarOpKind},
            record_buf::{Cigar, RecordBuf},
        },
        header::record::value::{Map, map::ReferenceSequence},
    };
    use std::num::NonZeroUsize;

    fn test_header() -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZeroUsize::new(1000).unwrap()),
            )
            .build()
    }

    /// Round-trips hand-built records through a real BAM byte stream so the
    /// test exercises the same decode path `compute_nvc` runs against in
    /// production, not just in-memory builder state.
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
    fn matches_hand_computed_counts() {
        let header = test_header();

        // Forward read: ACGT (length 4)
        let forward = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"ACGT".to_vec()))
            .build();

        // Reverse read: ACGT -> reverse complement is ACGT (palindromic)
        let reverse = RecordBuf::builder()
            .set_flags(Flags::REVERSE_COMPLEMENTED)
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"ACGT".to_vec()))
            .build();

        // Short read: AT (length 2) - this should determine the table length
        let short = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 2)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"AT".to_vec()))
            .build();

        // Low MAPQ read (should be filtered out)
        let low_mapq = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(10).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 3)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"ACG".to_vec()))
            .build();

        let records = vec![forward, reverse, short, low_mapq];
        let bam_records = to_bam_records(&header, &records);
        let table = compute_nvc(bam_records.into_iter().map(Ok), 30).unwrap();

        // Hand-computed from the fixture above (all three mapq>=30 reads
        // contribute at positions 0 and 1; low_mapq is filtered out):
        // - forward "ACGT": 0:A, 1:C, 2:G, 3:T
        // - reverse "ACGT" -> reverse-complemented -> "ACGT" (palindromic): 0:A, 1:C, 2:G, 3:T
        // - short "AT": 0:A, 1:T
        // Position 0: A = forward's A + reverse's A + short's A = 3
        // Position 1: C = forward's C + reverse's C = 2; T = short's T = 1
        // Table length is 2 (the LAST record processed is `short`, length 2),
        // matching the preserved upstream quirk (DIV-0008) even though the
        // forward/reverse reads were longer.

        assert_eq!(table.rows.len(), 2, "Table length should match last read's length");

        assert_eq!(table.rows[0].position, 0);
        assert_eq!(table.rows[0].a, 3);
        assert_eq!(table.rows[0].c, 0);
        assert_eq!(table.rows[0].g, 0);
        assert_eq!(table.rows[0].t, 0);
        assert_eq!(table.rows[0].n, 0);
        assert_eq!(table.rows[0].x, 0);

        assert_eq!(table.rows[1].position, 1);
        assert_eq!(table.rows[1].a, 0);
        assert_eq!(table.rows[1].c, 2);
        assert_eq!(table.rows[1].g, 0);
        assert_eq!(table.rows[1].t, 1);
        assert_eq!(table.rows[1].n, 0);
        assert_eq!(table.rows[1].x, 0);
    }

    #[test]
    fn render_output_format() {
        let table = NvcTable {
            rows: vec![
                NvcRow { position: 0, a: 10, c: 5, g: 3, t: 8, n: 1, x: 0 },
                NvcRow { position: 1, a: 8, c: 12, g: 6, t: 4, n: 0, x: 1 },
            ],
        };

        let output = render_nvc_table(&table);
        let expected = "Position\tA\tC\tG\tT\tN\tX\n0\t10\t5\t3\t8\t1\t0\n1\t8\t12\t6\t4\t0\t1";
        
        assert_eq!(output, expected);
    }
}