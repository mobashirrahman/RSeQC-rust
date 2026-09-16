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

pub fn render_single_table(p: &SingleEndProfile) -> String {
    let mut lines = vec!["Position\tInsert_nt\tNon_insert_nt".to_string()];
    for (i, &c) in p.clip_count.iter().enumerate() {
        let non_insert = p.total_read - c;
        lines.push(format!("{i}\t{}\t{}", fmt_float(c), fmt_float(non_insert)));
    }
    lines.join("\n")
}

pub fn render_single_r_script(p: &SingleEndProfile, out_prefix: &str) -> String {
    let read_pos: Vec<String> = (0..p.clip_count.len()).map(|i| i.to_string()).collect();
    let insert_strs: Vec<String> = p.clip_count.iter().map(|&c| fmt_float(c)).collect();

    let mut lines = Vec::new();
    lines.push(format!("pdf(\"{out_prefix}.insertion_profile.pdf\")"));
    lines.push(format!("read_pos=c({})", read_pos.join(",")));
    lines.push(format!("insert_count=c({})", insert_strs.join(",")));
    lines.push(format!("noninsert_count= {} - insert_count", p.total_read));
    lines.push(
        "plot(read_pos, insert_count*100/(insert_count+noninsert_count),col=\"blue\",main=\"Insertion profile\",xlab=\"Position of read\",ylab=\"Insertion %\",type=\"b\")"
            .to_string(),
    );
    lines.push("dev.off()".to_string());
    lines.join("\n")
}

pub fn render_paired_table(p: &PairedEndProfile) -> String {
    let mut lines = vec!["Position\tInsert_nt\tNon_insert_nt".to_string(), "Read-1:".to_string()];
    for (i, &c) in p.r1_clip_count.iter().enumerate() {
        let non_insert = p.total_read1 - c;
        lines.push(format!("{i}\t{}\t{}", fmt_float(c), fmt_float(non_insert)));
    }
    lines.push("Read-2:".to_string());
    for (i, &c) in p.r2_clip_count.iter().enumerate() {
        let non_insert = p.total_read2 - c;
        lines.push(format!("{i}\t{}\t{}", fmt_float(c), fmt_float(non_insert)));
    }
    lines.join("\n")
}

pub fn render_paired_r_script(p: &PairedEndProfile, out_prefix: &str) -> String {
    let read_pos: Vec<String> = (0..p.r1_clip_count.len()).map(|i| i.to_string()).collect();
    let r1_strs: Vec<String> = p.r1_clip_count.iter().map(|&c| fmt_float(c)).collect();
    let r2_strs: Vec<String> = p.r2_clip_count.iter().map(|&c| fmt_float(c)).collect();

    let mut lines = Vec::new();
    lines.push(format!("pdf(\"{out_prefix}.insertion_profile.R1.pdf\")"));
    lines.push(format!("read_pos=c({})", read_pos.join(",")));
    lines.push(format!("r1_insert_count=c({})", r1_strs.join(",")));
    lines.push(format!("r1_noninsert_count = {} - r1_insert_count", p.total_read1));
    lines.push(
        "plot(read_pos, r1_insert_count*100/(r1_insert_count + r1_noninsert_count),col=\"blue\",main=\"Insertion profile\",xlab=\"Position of read (read-1)\",ylab=\"Insertion %\",type=\"b\")"
            .to_string(),
    );
    lines.push("dev.off()".to_string());

    lines.push(format!("pdf(\"{out_prefix}.insertion_profile.R2.pdf\")"));
    lines.push(format!("read_pos=c({})", read_pos.join(",")));
    lines.push(format!("r2_insert_count=c({})", r2_strs.join(",")));
    lines.push(format!("r2_noninsert_count = {} - r2_insert_count", p.total_read2));
    lines.push(
        "plot(read_pos, r2_insert_count*100/(r2_insert_count + r2_noninsert_count),col=\"blue\",main=\"Insertion profile\",xlab=\"Position of read (read-2)\",ylab=\"Insertion %\",type=\"b\")"
            .to_string(),
    );
    lines.push("dev.off()".to_string());

    lines.join("\n")
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
            "Position\tInsert_nt\tNon_insert_nt\n0\t1.0\t3.0\n1\t3.0\t1.0"
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
dev.off()";
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
dev.off()";
        assert_eq!(output, expected);
    }
}
