//! Port of `insertion_profile.py`: calculate the distribution of inserted
//! nucleotides across reads. Contract: see `compatibility/commands.yaml`
//! entry `insertion_profile.py`; algorithm ported from
//! `insertion_profile()` in `oracle/upstream-src/src/qcmodule/SAM.py`
//! (lines 3362-3481+), which is structurally identical to
//! `clipping_profile()` with `type="I"` instead of `type="S"` and
//! different labels/file names.
//!
//! Reuses [`crate::clipping_profile::compute_single_end`] and
//! [`crate::clipping_profile::compute_paired_end`] directly (called with
//! `clip_char = b'I'`) rather than duplicating the aggregation logic —
//! only the rendering (column headers, R variable names, plot titles) and
//! file names differ from clipping_profile.py. Same DIV-0010-pattern
//! last-record-length quirk applies here too (inherited from the shared
//! compute functions).

use crate::clipping_profile::{PairedEndProfile, SingleEndProfile};

fn fmt_float(n: u64) -> String {
    format!("{n}.0")
}

/// Every line, including the last, ends with `\n` (upstream's plain
/// `print(...)` calls each add their own trailing newline) -- applies
/// to all four render functions in this module, same as
/// clipping_profile.rs (structurally identical upstream function).
pub fn render_single_table(p: &SingleEndProfile) -> String {
    let mut out = String::from("Position\tInsert_nt\tNon_insert_nt\n");
    for (i, &c) in p.clip_count.iter().enumerate() {
        let non_insert = p.total_read - c;
        out.push_str(&format!("{i}\t{}\t{}\n", fmt_float(c), fmt_float(non_insert)));
    }
    out
}

pub fn render_single_r_script(p: &SingleEndProfile, out_prefix: &str) -> String {
    let read_pos: Vec<String> = (0..p.clip_count.len()).map(|i| i.to_string()).collect();
    let insert_strs: Vec<String> = p.clip_count.iter().map(|&c| fmt_float(c)).collect();

    format!(
        "pdf(\"{out_prefix}.insertion_profile.pdf\")\nread_pos=c({})\ninsert_count=c({})\nnoninsert_count= {} - insert_count\nplot(read_pos, insert_count*100/(insert_count+noninsert_count),col=\"blue\",main=\"Insertion profile\",xlab=\"Position of read\",ylab=\"Insertion %\",type=\"b\")\ndev.off()\n",
        read_pos.join(","),
        insert_strs.join(","),
        p.total_read,
    )
}

pub fn render_paired_table(p: &PairedEndProfile) -> String {
    let mut out = String::from("Position\tInsert_nt\tNon_insert_nt\nRead-1:\n");
    for (i, &c) in p.r1_clip_count.iter().enumerate() {
        let non_insert = p.total_read1 - c;
        out.push_str(&format!("{i}\t{}\t{}\n", fmt_float(c), fmt_float(non_insert)));
    }
    out.push_str("Read-2:\n");
    for (i, &c) in p.r2_clip_count.iter().enumerate() {
        let non_insert = p.total_read2 - c;
        out.push_str(&format!("{i}\t{}\t{}\n", fmt_float(c), fmt_float(non_insert)));
    }
    out
}

pub fn render_paired_r_script(p: &PairedEndProfile, out_prefix: &str) -> String {
    let read_pos: Vec<String> = (0..p.r1_clip_count.len()).map(|i| i.to_string()).collect();
    let r1_strs: Vec<String> = p.r1_clip_count.iter().map(|&c| fmt_float(c)).collect();
    let r2_strs: Vec<String> = p.r2_clip_count.iter().map(|&c| fmt_float(c)).collect();
    let read_pos_csv = read_pos.join(",");

    format!(
        "pdf(\"{out_prefix}.insertion_profile.R1.pdf\")\nread_pos=c({read_pos_csv})\nr1_insert_count=c({})\nr1_noninsert_count = {} - r1_insert_count\nplot(read_pos, r1_insert_count*100/(r1_insert_count + r1_noninsert_count),col=\"blue\",main=\"Insertion profile\",xlab=\"Position of read (read-1)\",ylab=\"Insertion %\",type=\"b\")\ndev.off()\n\
pdf(\"{out_prefix}.insertion_profile.R2.pdf\")\nread_pos=c({read_pos_csv})\nr2_insert_count=c({})\nr2_noninsert_count = {} - r2_insert_count\nplot(read_pos, r2_insert_count*100/(r2_insert_count + r2_noninsert_count),col=\"blue\",main=\"Insertion profile\",xlab=\"Position of read (read-2)\",ylab=\"Insertion %\",type=\"b\")\ndev.off()\n",
        r1_strs.join(","),
        p.total_read1,
        r2_strs.join(","),
        p.total_read2,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn render_single_table_uses_insert_labels() {
        let profile = SingleEndProfile {
            total_read: 4,
            clip_count: vec![1, 3],
        };
        assert_eq!(
            render_single_table(&profile),
            "Position\tInsert_nt\tNon_insert_nt\n0\t1.0\t3.0\n1\t3.0\t1.0\n"
        );
    }

    #[test]
    fn render_single_r_script_exact_text() {
        // Independently derived by running the equivalent Python
        // print()/%-format expressions from oracle/upstream-src/src/
        // qcmodule/SAM.py lines 3410-3415 via `python3 -c`.
        let profile = SingleEndProfile {
            total_read: 4,
            clip_count: vec![1, 3],
        };
        let output = render_single_r_script(&profile, "test_output");
        let expected = "pdf(\"test_output.insertion_profile.pdf\")\n\
read_pos=c(0,1)\n\
insert_count=c(1.0,3.0)\n\
noninsert_count= 4 - insert_count\n\
plot(read_pos, insert_count*100/(insert_count+noninsert_count),col=\"blue\",main=\"Insertion profile\",xlab=\"Position of read\",ylab=\"Insertion %\",type=\"b\")\n\
dev.off()\n";
        assert_eq!(output, expected);
    }

    #[test]
    fn render_paired_r_script_exact_text() {
        // Independently derived from oracle/upstream-src/src/qcmodule/
        // SAM.py lines 3471-3481 via `python3 -c`.
        let profile = PairedEndProfile {
            total_read1: 5,
            total_read2: 4,
            r1_clip_count: vec![2, 0],
            r2_clip_count: vec![0, 1],
        };
        let output = render_paired_r_script(&profile, "test_output");
        let expected = "pdf(\"test_output.insertion_profile.R1.pdf\")\n\
read_pos=c(0,1)\n\
r1_insert_count=c(2.0,0.0)\n\
r1_noninsert_count = 5 - r1_insert_count\n\
plot(read_pos, r1_insert_count*100/(r1_insert_count + r1_noninsert_count),col=\"blue\",main=\"Insertion profile\",xlab=\"Position of read (read-1)\",ylab=\"Insertion %\",type=\"b\")\n\
dev.off()\n\
pdf(\"test_output.insertion_profile.R2.pdf\")\n\
read_pos=c(0,1)\n\
r2_insert_count=c(0.0,1.0)\n\
r2_noninsert_count = 4 - r2_insert_count\n\
plot(read_pos, r2_insert_count*100/(r2_insert_count + r2_noninsert_count),col=\"blue\",main=\"Insertion profile\",xlab=\"Position of read (read-2)\",ylab=\"Insertion %\",type=\"b\")\n\
dev.off()\n";
        assert_eq!(output, expected);
    }
}
