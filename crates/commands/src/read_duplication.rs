//! Port of `read_duplication.py`: calculate sequence-based and mapping-based read duplication rates.
//! Contract: see `compatibility/commands.yaml` entry `read_duplication.py`; algorithm
//! ported from `readDupRate()` in `oracle/upstream-src/src/qcmodule/SAM.py`
//! (lines 3161-3237).
//!
//! Notable upstream behavior preserved exactly: unmapped reads, QC-fail reads,
//! and low MAPQ reads are filtered out; three outputs are generated: sequence-based
//! and position-based duplication rate tables, plus an R script for plotting.

use std::collections::HashMap;
use std::io;

use noodles_sam as sam;
use noodles_sam::alignment::record::cigar::Op;
use rseqc_formats::cigar::fetch_exon_blocks;


#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct DuplicationHistograms {
    pub seq_occurrence_counts: Vec<(u64, u64)>, // (occurrence, distinct_sequence_count), sorted by occurrence
    pub pos_occurrence_counts: Vec<(u64, u64)>, // (occurrence, distinct_position_count), sorted by occurrence
}

/// Computes sequence-based and position-based duplication histograms for a sequence of BAM records,
/// matching the upstream algorithm exactly.
/// Pure computation: no I/O, no printing.
pub fn compute_duplication<I>(
    records: I,
    header: &sam::Header,
    q_cut: u8,
) -> io::Result<DuplicationHistograms>
where
    I: IntoIterator<Item = io::Result<noodles_bam::Record>>,
{
    let mut seq_dup = HashMap::new(); // exact uppercase sequence string -> occurrence count
    let mut pos_dup = HashMap::new(); // "chrom:start:exon_boundary_string" -> occurrence count

    for result in records {
        let record = result?;
        
        // Filter out unmapped reads
        if record.flags().is_unmapped() {
            continue;
        }
        
        // Filter out QC fail reads
        if record.flags().is_qc_fail() {
            continue;
        }
        
        // Filter out low MAPQ reads
        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(255);
        if mapq < q_cut {
            continue;
        }

        // Get uppercase sequence
        let sequence = record.sequence();
        let rna_read: String = sequence.iter().map(|b| b as char).collect();
        seq_dup.insert(rna_read.clone(), seq_dup.get(&rna_read).unwrap_or(&0) + 1);

        // Get reference sequence name
        let ref_id = match record.reference_sequence_id().transpose() {
            Ok(id) => id,
            Err(e) => return Err(e),
        };
        let ref_seqs = header.reference_sequences();
        let ref_index = match ref_id {
            Some(id) => id,
            None => return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid reference sequence ID")),
        };
        let (ref_name, _) = match ref_seqs.get_index(ref_index) {
            Some(name) => name,
            None => return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid reference sequence ID")),
        };
        let chrom = ref_name.to_string();

        // Get alignment start (convert from 1-based to 0-based)
        let start = match record.alignment_start().transpose() {
            Ok(Some(pos)) => pos.get() - 1,
            Ok(None) => return Err(io::Error::new(io::ErrorKind::InvalidData, "Missing alignment start")),
            Err(e) => return Err(e),
        };

        // Get exon blocks using CIGAR
        let cigar_ops: Vec<Op> = record.cigar().iter().collect::<Result<Vec<_>, _>>()?;
        let exon_blocks = fetch_exon_blocks(start, cigar_ops);

        // Build exon boundary string
        let mut exon_boundary = String::new();
        for (ex_start, ex_end) in exon_blocks {
            exon_boundary.push_str(&format!("{}-{}:", ex_start, ex_end));
        }

        // Create position duplication key
        let key = format!("{}:{}:{}", chrom, start, exon_boundary);
        pos_dup.insert(key.clone(), pos_dup.get(&key).unwrap_or(&0) + 1);
    }

    // Convert seq_dup to occurrence counts: occurrence -> distinct_sequence_count
    let mut seq_count = HashMap::new();
    for v in seq_dup.values() {
        let count = v;
        seq_count.insert(*count, seq_count.get(count).unwrap_or(&0) + 1);
    }

    // Convert pos_dup to occurrence counts: occurrence -> distinct_position_count
    let mut pos_count = HashMap::new();
    for v in pos_dup.values() {
        let count = v;
        pos_count.insert(*count, pos_count.get(count).unwrap_or(&0) + 1);
    }

    // Sort both by occurrence (ascending). NOTE: upstream's `up_bound` is
    // NEVER used to filter seqDup_count/posDup_count -- it only bounds the
    // R plot's x-axis (`xlim=c(1,up_bound)`) and the legend x-position, both
    // purely in render_dup_r_script. Filtering the data itself here would
    // be a real behavioral divergence, not a cosmetic one.
    let mut seq_occurrence_counts: Vec<(u64, u64)> = seq_count.into_iter().collect();
    seq_occurrence_counts.sort_by_key(|&(occ, _)| occ);

    let mut pos_occurrence_counts: Vec<(u64, u64)> = pos_count.into_iter().collect();
    pos_occurrence_counts.sort_by_key(|&(occ, _)| occ);

    Ok(DuplicationHistograms {
        seq_occurrence_counts,
        pos_occurrence_counts,
    })
}

/// Renders the sequence duplication histogram as tab-separated text, matching upstream's .seq.DupRate.xls format
pub fn render_seq_dup_table(hist: &DuplicationHistograms) -> String {
    let mut lines = Vec::new();
    
    // Header row
    lines.push("Occurrence\tUniqReadNumber".to_string());
    
    // Data rows sorted by occurrence (ascending)
    for (occurrence, distinct_count) in &hist.seq_occurrence_counts {
        lines.push(format!("{}\t{}", occurrence, distinct_count));
    }
    
    lines.join("\n")
}

/// Renders the position duplication histogram as tab-separated text, matching upstream's .pos.DupRate.xls format
pub fn render_pos_dup_table(hist: &DuplicationHistograms) -> String {
    let mut lines = Vec::new();
    
    // Header row
    lines.push("Occurrence\tUniqReadNumber".to_string());
    
    // Data rows sorted by occurrence (ascending)
    for (occurrence, distinct_count) in &hist.pos_occurrence_counts {
        lines.push(format!("{}\t{}", occurrence, distinct_count));
    }
    
    lines.join("\n")
}

/// Renders the duplication histograms as an R script, matching upstream's .DupRate_plot.r format
pub fn render_dup_r_script(
    hist: &DuplicationHistograms,
    output_prefix: &str,
    up_bound: u64,
) -> String {
    let mut lines = Vec::new();
    
    // PDF output line
    lines.push(format!("pdf(\"{}.DupRate_plot.pdf\")", output_prefix));
    
    // Set plotting parameters
    lines.push("par(mar=c(5,4,4,5),las=0)".to_string());
    
    // Build comma-separated occurrence and count vectors for sequence data
    let seq_occurrences: Vec<String> = hist
        .seq_occurrence_counts
        .iter()
        .map(|&(occ, _)| occ.to_string())
        .collect();
    let seq_counts: Vec<String> = hist
        .seq_occurrence_counts
        .iter()
        .map(|&(_, count)| count.to_string())
        .collect();
    
    // Build comma-separated occurrence and count vectors for position data
    let pos_occurrences: Vec<String> = hist
        .pos_occurrence_counts
        .iter()
        .map(|&(occ, _)| occ.to_string())
        .collect();
    let pos_counts: Vec<String> = hist
        .pos_occurrence_counts
        .iter()
        .map(|&(_, count)| count.to_string())
        .collect();
    
    // R script data lines
    lines.push(format!("seq_occ=c({})", seq_occurrences.join(",")));
    lines.push(format!("seq_uniqRead=c({})", seq_counts.join(",")));
    lines.push(format!("pos_occ=c({})", pos_occurrences.join(",")));
    lines.push(format!("pos_uniqRead=c({})", pos_counts.join(",")));
    
    // `log10(pos_uniqRead)`/`log10(seq_uniqRead)` reference the R vectors
    // already defined above by name -- NOT one log10() call per value
    // (that would pass N extra positional arguments to plot()/points(),
    // breaking the call). Only `up_bound` is a real Rust-side substitution
    // here (upstream's `%d % up_bound`).
    let max_legend_pos = if up_bound > 200 { up_bound - 200 } else { 1 };

    lines.push(format!(
        "plot(pos_occ,log10(pos_uniqRead),ylab='Number of Reads (log10)',xlab='Occurrence of read',pch=4,cex=0.8,col='blue',xlim=c(1,{}),yaxt='n')",
        up_bound
    ));
    lines.push("points(seq_occ,log10(seq_uniqRead),pch=20,cex=0.8,col='red')".to_string());
    
    // `ym` is an R-side computation upstream (`ym=floor(max(log10(pos_uniqRead)))`),
    // printed as literal text -- NOT something computed in Rust and
    // substituted numerically. The `legend(...)` and `axis(side=2,...)`
    // lines both reference this R variable by name, same as upstream.
    lines.push("ym=floor(max(log10(pos_uniqRead)))".to_string());
    lines.push(format!("legend({},ym,legend=c('Sequence-based','Mapping-based'),col=c('blue','red'),pch=c(4,20))", max_legend_pos));
    lines.push("axis(side=2,at=0:ym,labels=0:ym)".to_string());

    // Right axis for percentage calculation (exact upstream text)
    lines.push("axis(side=4,at=c(log10(pos_uniqRead[1]),log10(pos_uniqRead[2]),log10(pos_uniqRead[3]),log10(pos_uniqRead[4])), labels=c(round(pos_uniqRead[1]*100/sum(pos_uniqRead*pos_occ)),round(pos_uniqRead[2]*100/sum(pos_uniqRead*pos_occ)),round(pos_uniqRead[3]*100/sum(pos_uniqRead*pos_occ)),round(pos_uniqRead[4]*100/sum(pos_uniqRead*pos_occ))))".to_string());
    lines.push("mtext(4, text = \"Reads %\", line = 2)".to_string());
    
    // End PDF device
    lines.push("dev.off()".to_string());
    
    lines.join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_core::Position;
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
            .add_reference_sequence(
                "chr2", 
                Map::<ReferenceSequence>::new(NonZeroUsize::new(2000).unwrap()),
            )
            .build()
    }

    /// Round-trips hand-built records through a real BAM byte stream so the
    /// test exercises the same decode path `compute_duplication` runs against in
    /// production, not just in-memory builder state.
    fn to_bam_records(header: &sam::Header, records: &[RecordBuf]) -> Vec<noodles_bam::Record> {
        let mut buf = Vec::new();
        {
            let mut writer = noodles_bam::io::Writer::new(&mut buf);
            writer.write_header(header).unwrap();
            for r in records {
                writer.write_alignment_record(header, r).unwrap();
            }
        }
        let mut reader = noodles_bam::io::Reader::new(buf.as_slice());
        reader.read_header().unwrap();
        reader.records().map(|r| r.unwrap()).collect()
    }

    #[test]
    fn matches_hand_computed_counts() {
        let header = test_header();

        // Read with identical sequence "ATCG" (should create seq duplication)
        let read_1 = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"ATCG".to_vec()))
            .set_alignment_start(Position::new(1).unwrap())
            .build();

        // Read with identical sequence "ATCG" (should create seq duplication)
        let read_2 = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"ATCG".to_vec()))
            .set_alignment_start(Position::new(1).unwrap())
            .build();

        // Read with different sequence "GCTA" (no duplication)
        let read_3 = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"GCTA".to_vec()))
            .set_alignment_start(Position::new(10).unwrap())
            .build();

        // Read with identical position and CIGAR but different sequence (should create pos duplication)
        let read_4 = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_mapping_quality(MappingQuality::new(40).unwrap())
            .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
            .set_sequence(sam::alignment::record_buf::Sequence::from(b"TTTT".to_vec()))
            .set_alignment_start(Position::new(1).unwrap())
            .build();

        let records = vec![read_1, read_2, read_3, read_4];
        let bam_records = to_bam_records(&header, &records);
        let hist = compute_duplication(bam_records.into_iter().map(Ok), &header, 30).unwrap();

        // Expected: 
        // seq: occurrence 2 -> 2 distinct sequences ("ATCG" appears twice, "GCTA" and "TTTT" appear once each)
        // pos: occurrence 3 -> 1 distinct position (position 1 has 3 reads)
        assert_eq!(hist.seq_occurrence_counts, vec![(1, 2), (2, 1)]); // 1 appears 2 times, 2 appears 1 time
        assert_eq!(hist.pos_occurrence_counts, vec![(1, 1), (3, 1)]); // 1 appears 1 time, 3 appears 1 time
    }

    #[test]
    fn render_seq_dup_table_format() {
        let hist = DuplicationHistograms {
            seq_occurrence_counts: vec![(1, 10), (2, 5), (3, 2)],
            pos_occurrence_counts: vec![(1, 15), (2, 3), (4, 1)],
        };

        let output = render_seq_dup_table(&hist);
        let expected = "Occurrence\tUniqReadNumber\n1\t10\n2\t5\n3\t2";
        
        assert_eq!(output, expected);
    }

    #[test]
    fn render_pos_dup_table_format() {
        let hist = DuplicationHistograms {
            seq_occurrence_counts: vec![(1, 10), (2, 5), (3, 2)],
            pos_occurrence_counts: vec![(1, 15), (2, 3), (4, 1)],
        };

        let output = render_pos_dup_table(&hist);
        let expected = "Occurrence\tUniqReadNumber\n1\t15\n2\t3\n4\t1";
        
        assert_eq!(output, expected);
    }

    #[test]
    fn render_dup_r_script_format() {
        let hist = DuplicationHistograms {
            seq_occurrence_counts: vec![(1, 10), (2, 5), (3, 2)],
            pos_occurrence_counts: vec![(1, 15), (2, 3), (4, 1)],
        };

        let output = render_dup_r_script(&hist, "test_output", 500);

        // Exact text, derived independently by running the equivalent
        // Python print()/%-format expressions from
        // oracle/upstream-src/src/qcmodule/SAM.py lines 3222-3236 via
        // `python3 -c`, not by copying this function's own output back in.
        let expected = "pdf(\"test_output.DupRate_plot.pdf\")\n\
par(mar=c(5,4,4,5),las=0)\n\
seq_occ=c(1,2,3)\n\
seq_uniqRead=c(10,5,2)\n\
pos_occ=c(1,2,4)\n\
pos_uniqRead=c(15,3,1)\n\
plot(pos_occ,log10(pos_uniqRead),ylab='Number of Reads (log10)',xlab='Occurrence of read',pch=4,cex=0.8,col='blue',xlim=c(1,500),yaxt='n')\n\
points(seq_occ,log10(seq_uniqRead),pch=20,cex=0.8,col='red')\n\
ym=floor(max(log10(pos_uniqRead)))\n\
legend(300,ym,legend=c('Sequence-based','Mapping-based'),col=c('blue','red'),pch=c(4,20))\n\
axis(side=2,at=0:ym,labels=0:ym)\n\
axis(side=4,at=c(log10(pos_uniqRead[1]),log10(pos_uniqRead[2]),log10(pos_uniqRead[3]),log10(pos_uniqRead[4])), labels=c(round(pos_uniqRead[1]*100/sum(pos_uniqRead*pos_occ)),round(pos_uniqRead[2]*100/sum(pos_uniqRead*pos_occ)),round(pos_uniqRead[3]*100/sum(pos_uniqRead*pos_occ)),round(pos_uniqRead[4]*100/sum(pos_uniqRead*pos_occ))))\n\
mtext(4, text = \"Reads %\", line = 2)\n\
dev.off()";

        assert_eq!(output, expected);
    }

    #[test]
    fn up_bound_does_not_filter_aggregated_data() {
        // Regression test: upstream's `up_bound` parameter is used ONLY in
        // render_dup_r_script (the plot's xlim and legend position) -- it
        // never filters seqDup_count/posDup_count in readDupRate. An
        // earlier draft of this port incorrectly filtered occurrence
        // counts above up_bound out of the aggregation; this asserts a
        // high-occurrence count survives intact regardless of what
        // up_bound is later rendered with.
        let header = test_header();

        let records: Vec<_> = (0..10)
            .map(|_| {
                RecordBuf::builder()
                    .set_flags(Flags::empty())
                    .set_reference_sequence_id(0)
                    .set_mapping_quality(MappingQuality::new(40).unwrap())
                    .set_cigar(Cigar::from(vec![Op::new(CigarOpKind::Match, 4)]))
                    .set_sequence(sam::alignment::record_buf::Sequence::from(b"AAAA".to_vec()))
                    .set_alignment_start(Position::new(1).unwrap())
                    .build()
            })
            .collect();

        let bam_records = to_bam_records(&header, &records);
        let hist = compute_duplication(bam_records.into_iter().map(Ok), &header, 30).unwrap();

        // All 10 reads share the same sequence: one distinct sequence with
        // occurrence count 10 -- not filtered, regardless of up_bound.
        assert_eq!(hist.seq_occurrence_counts, vec![(10, 1)]);
    }
}