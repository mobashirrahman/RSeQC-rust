//! Port of `bam2fq.py`: convert BAM alignments to FASTQ. Contract: see
//! `compatibility/commands.yaml` entry `bam2fq.py`; algorithm ported from
//! `ParseBAM.bam2fq()` in `oracle/upstream-src/src/qcmodule/SAM.py`
//! (lines 2678-2767).
//!
//! Notable upstream behavior preserved exactly: EVERY record in the file is
//! converted (no QC-fail/duplicate/secondary/unmapped filtering); a missing
//! sequence or quality string is a hard error, not a skip; reverse-strand
//! reads are reverse-complemented; in paired mode, a record that is neither
//! read1 nor read2 is silently dropped (no `else` branch upstream).
//!
//! Known gap: `-c/--compress` (gzip output) is not implemented here yet —
//! disclosed follow-up (DIV-0007). SAM-text input (DIV-0002/0004) is
//! supported at the CLI layer via `rseqc_formats::open_alignments`; this
//! module's own functions were unaffected (already generic over
//! `IntoIterator<Item = io::Result<bam::Record>>`).

use std::io::{self, Write};

use noodles_bam as bam;

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Bam2FqCounts {
    pub read1: u64,
    pub read2: u64,
    pub single: u64,
}

fn complement(base: u8) -> u8 {
    // Matches upstream's `str.maketrans("ACGTNX", "TGCANX")`: only these six
    // letters are translated, everything else passes through unchanged.
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

/// Extracts (name, sequence, quality) as FASTQ-ready bytes, applying the
/// reverse-complement upstream applies for reverse-strand reads. Quality
/// bytes are ASCII (Phred+33); a BAM quality array that is empty or entirely
/// `0xff` (BAM's "no quality string" sentinel, matching pysam's
/// `query_qualities is None`) is treated as missing, same as upstream.
fn fastq_values(record: &bam::Record) -> io::Result<(String, Vec<u8>, Vec<u8>)> {
    let name = record
        .name()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "record has no query name"))?
        .to_string();

    let sequence = record.sequence();
    if sequence.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Read {name:?} does not contain sequence information"),
        ));
    }

    let quality_scores = record.quality_scores();
    let missing_quality = quality_scores.as_bytes().iter().all(|&q| q == 0xff);
    if missing_quality {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Read {name:?} does not contain quality information"),
        ));
    }

    let mut seq: Vec<u8> = sequence.iter().collect();
    let mut qual: Vec<u8> = quality_scores.iter().map(|q| q + 33).collect();

    if record.flags().is_reverse_complemented() {
        seq.reverse();
        for b in seq.iter_mut() {
            *b = complement(*b);
        }
        qual.reverse();
    }

    Ok((name, seq, qual))
}

fn write_fastq_record<W: Write>(w: &mut W, name: &str, seq: &[u8], qual: &[u8]) -> io::Result<()> {
    writeln!(w, "@{name}")?;
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    w.write_all(qual)?;
    w.write_all(b"\n")?;
    Ok(())
}

fn suffixed(name: String, suffix: &str) -> String {
    if name.ends_with(suffix) {
        name
    } else {
        format!("{name}{suffix}")
    }
}

/// Paired mode: read1 records go to `out1` (name suffixed `/1` unless
/// already present), read2 records go to `out2` (suffixed `/2`); anything
/// that's neither is silently dropped, matching upstream's `if`/`elif` with
/// no `else`.
pub fn write_paired<I, W1, W2>(
    records: I,
    out1: &mut W1,
    out2: &mut W2,
) -> io::Result<Bam2FqCounts>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
    W1: Write,
    W2: Write,
{
    let mut counts = Bam2FqCounts::default();

    for result in records {
        let record = result?;
        let flags = record.flags();
        let (name, seq, qual) = fastq_values(&record)?;

        if flags.is_first_segment() {
            write_fastq_record(out1, &suffixed(name, "/1"), &seq, &qual)?;
            counts.read1 += 1;
        } else if flags.is_last_segment() {
            write_fastq_record(out2, &suffixed(name, "/2"), &seq, &qual)?;
            counts.read2 += 1;
        }
    }

    Ok(counts)
}

/// Single-end mode: every record is written, unsuffixed, to `out`.
pub fn write_single<I, W>(records: I, out: &mut W) -> io::Result<Bam2FqCounts>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
    W: Write,
{
    let mut counts = Bam2FqCounts::default();

    for result in records {
        let record = result?;
        let (name, seq, qual) = fastq_values(&record)?;
        write_fastq_record(out, &name, &seq, &qual)?;
        counts.single += 1;
    }

    Ok(counts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::{
        self as sam,
        alignment::io::Write as _,
        alignment::record::Flags,
        alignment::record_buf::{QualityScores, RecordBuf, Sequence},
    };

    fn header() -> sam::Header {
        sam::Header::default()
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
    fn paired_mode_splits_and_reverse_complements() {
        let h = header();

        // read1, forward strand, already-suffixed name.
        let r1 = RecordBuf::builder()
            .set_name("read.a/1")
            .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .set_sequence(Sequence::from(b"ACGT".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30, 30, 30, 30]))
            .build();

        // read2, reverse strand, unsuffixed name.
        let r2 = RecordBuf::builder()
            .set_name("read.b")
            .set_flags(Flags::SEGMENTED | Flags::LAST_SEGMENT | Flags::REVERSE_COMPLEMENTED)
            .set_sequence(Sequence::from(b"ACGT".to_vec()))
            .set_quality_scores(QualityScores::from(vec![0, 10, 20, 30]))
            .build();

        let bam_records = to_bam_records(&h, &[r1, r2]);

        let mut out1 = Vec::new();
        let mut out2 = Vec::new();
        let counts = write_paired(bam_records.into_iter().map(Ok), &mut out1, &mut out2).unwrap();

        assert_eq!(
            counts,
            Bam2FqCounts {
                read1: 1,
                read2: 1,
                single: 0
            }
        );

        let out1_text = String::from_utf8(out1).unwrap();
        assert_eq!(out1_text, "@read.a/1\nACGT\n+\n????\n");

        // ACGT reverse-complemented is ACGT (palindromic here on purpose is
        // avoided): complement(A,C,G,T) = T,G,C,A then reversed -> A,C,G,T.
        // Use a non-palindromic check on quality instead, which has no
        // ambiguity: quality [0,10,20,30] -> ASCII "!+5?" reversed -> "?5+!".
        let out2_text = String::from_utf8(out2).unwrap();
        assert_eq!(out2_text, "@read.b/2\nACGT\n+\n?5+!\n");
    }

    #[test]
    fn single_end_mode_writes_every_record() {
        let h = header();

        let r = RecordBuf::builder()
            .set_name("only")
            .set_sequence(Sequence::from(b"TTTT".to_vec()))
            .set_quality_scores(QualityScores::from(vec![40, 40, 40, 40]))
            .build();

        let bam_records = to_bam_records(&h, &[r]);
        let mut out = Vec::new();
        let counts = write_single(bam_records.into_iter().map(Ok), &mut out).unwrap();

        assert_eq!(
            counts,
            Bam2FqCounts {
                read1: 0,
                read2: 0,
                single: 1
            }
        );
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "@only\nTTTT\n+\nIIII\n"
        );
    }

    #[test]
    fn missing_quality_is_an_error() {
        let h = header();
        // An empty quality_scores buffer with a non-empty sequence is BAM's
        // encoding of "no quality string": the writer auto-fills it with
        // the 0xff sentinel per base, which is exactly what real
        // pysam/samtools-produced "missing quality" BAM records look like.
        let r = RecordBuf::builder()
            .set_name("noqual")
            .set_sequence(Sequence::from(b"ACGT".to_vec()))
            .set_quality_scores(QualityScores::from(Vec::new()))
            .build();
        let bam_records = to_bam_records(&h, &[r]);
        let mut out = Vec::new();
        let err = write_single(bam_records.into_iter().map(Ok), &mut out).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }
}
