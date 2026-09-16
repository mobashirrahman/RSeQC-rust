//! Per-chromosome interval-set algebra (union/intersect/subtract),
//! replicating the *observable behavior* of `qcmodule.BED`'s
//! `unionBed3`/`intersectBed3`/`subtractBed3` (which are themselves thin
//! wrappers around bx-python's per-chromosome bitset union/AND/AND-NOT).
//!
//! Upstream converts every interval list to a bitset (one bit per genomic
//! position covered by ANY input interval on that chromosome), combines
//! bitsets, then re-extracts contiguous runs of set bits as merged
//! `[start, end)` intervals -- which discards all per-interval identity
//! (which original interval a position came from) and any interval-level
//! metadata, merging touching/overlapping intervals into single runs.
//! Rather than materializing a real per-base bitset (memory-heavy for
//! chromosome-scale coordinates), this module gets the same *result* via
//! direct interval-merge/intersect/subtract algorithms on sorted interval
//! lists -- behaviorally equivalent, not a literal bitset port.

use std::collections::HashMap;

/// A 0-based half-open `[start, end)` interval on one chromosome.
pub type Bed3 = (String, i64, i64);

/// Merges (per chromosome) all input intervals into maximal non-
/// overlapping, non-touching `[start, end)` runs, matching
/// `unionBed3`. Input order and any interval shorter than 1bp (start>=end)
/// are not meaningful upstream either (an empty/inverted interval
/// contributes nothing to the bitset); such intervals are dropped here too.
pub fn union_bed3(intervals: &[Bed3]) -> Vec<Bed3> {
    let by_chrom = group_by_chrom(intervals);
    let mut out = Vec::new();
    for (chrom, mut ivs) in by_chrom {
        merge_sorted(&mut ivs);
        for (s, e) in ivs {
            out.push((chrom.clone(), s, e));
        }
    }
    out.sort();
    out
}

/// Per-chromosome intersection of the union of `a` with the union of `b`,
/// matching `intersectBed3` (which bitset-ANDs, so it's inherently an
/// intersection of the *merged* coverage of each side, not pairwise
/// per-original-interval intersection).
pub fn intersect_bed3(a: &[Bed3], b: &[Bed3]) -> Vec<Bed3> {
    let a_by_chrom = merged_by_chrom(a);
    let b_by_chrom = merged_by_chrom(b);

    let mut out = Vec::new();
    for (chrom, a_ivs) in &a_by_chrom {
        let Some(b_ivs) = b_by_chrom.get(chrom) else { continue };
        for &(a_s, a_e) in a_ivs {
            for &(b_s, b_e) in b_ivs {
                let s = a_s.max(b_s);
                let e = a_e.min(b_e);
                if s < e {
                    out.push((chrom.clone(), s, e));
                }
            }
        }
    }
    out.sort();
    out
}

/// Per-chromosome subtraction of the union of `b` from the union of `a`,
/// matching `subtractBed3` (bitset AND-NOT of merged coverage).
pub fn subtract_bed3(a: &[Bed3], b: &[Bed3]) -> Vec<Bed3> {
    let a_by_chrom = merged_by_chrom(a);
    let b_by_chrom = merged_by_chrom(b);

    let mut out = Vec::new();
    for (chrom, a_ivs) in &a_by_chrom {
        let empty = Vec::new();
        let b_ivs = b_by_chrom.get(chrom).unwrap_or(&empty);
        for &(a_s, a_e) in a_ivs {
            out.extend(subtract_one(a_s, a_e, b_ivs).into_iter().map(|(s, e)| (chrom.clone(), s, e)));
        }
    }
    out.sort();
    out
}

fn group_by_chrom(intervals: &[Bed3]) -> HashMap<String, Vec<(i64, i64)>> {
    let mut by_chrom: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
    for (chrom, s, e) in intervals {
        if s < e {
            by_chrom.entry(chrom.clone()).or_default().push((*s, *e));
        }
    }
    by_chrom
}

fn merged_by_chrom(intervals: &[Bed3]) -> HashMap<String, Vec<(i64, i64)>> {
    let mut by_chrom = group_by_chrom(intervals);
    for ivs in by_chrom.values_mut() {
        merge_sorted(ivs);
    }
    by_chrom
}

/// Sorts and merges overlapping/touching (half-open-adjacent) intervals in place.
fn merge_sorted(ivs: &mut Vec<(i64, i64)>) {
    ivs.sort();
    let mut merged: Vec<(i64, i64)> = Vec::with_capacity(ivs.len());
    for &(s, e) in ivs.iter() {
        if let Some(last) = merged.last_mut() {
            if s <= last.1 {
                last.1 = last.1.max(e);
                continue;
            }
        }
        merged.push((s, e));
    }
    *ivs = merged;
}

/// Subtracts a set of already-merged, sorted, non-overlapping `b_ivs` from
/// one `[a_s, a_e)` interval, returning the remaining piece(s).
fn subtract_one(a_s: i64, a_e: i64, b_ivs: &[(i64, i64)]) -> Vec<(i64, i64)> {
    let mut out = Vec::new();
    let mut cursor = a_s;
    for &(b_s, b_e) in b_ivs {
        if b_e <= cursor {
            continue;
        }
        if b_s >= a_e {
            break;
        }
        if b_s > cursor {
            out.push((cursor, b_s.min(a_e)));
        }
        cursor = cursor.max(b_e);
        if cursor >= a_e {
            break;
        }
    }
    if cursor < a_e {
        out.push((cursor, a_e));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bed(chrom: &str, s: i64, e: i64) -> Bed3 {
        (chrom.to_string(), s, e)
    }

    #[test]
    fn union_merges_overlapping_and_touching_intervals() {
        let input = vec![bed("chr1", 10, 20), bed("chr1", 15, 25), bed("chr1", 25, 30), bed("chr1", 100, 110)];
        let result = union_bed3(&input);
        assert_eq!(result, vec![bed("chr1", 10, 30), bed("chr1", 100, 110)]);
    }

    #[test]
    fn union_keeps_chromosomes_separate() {
        let input = vec![bed("chr1", 0, 10), bed("chr2", 0, 10)];
        let result = union_bed3(&input);
        assert_eq!(result, vec![bed("chr1", 0, 10), bed("chr2", 0, 10)]);
    }

    #[test]
    fn intersect_finds_overlap_of_merged_sides() {
        let a = vec![bed("chr1", 0, 10), bed("chr1", 20, 30)];
        let b = vec![bed("chr1", 5, 25)];
        let result = intersect_bed3(&a, &b);
        assert_eq!(result, vec![bed("chr1", 5, 10), bed("chr1", 20, 25)]);
    }

    #[test]
    fn intersect_empty_when_no_overlap_or_chrom_missing() {
        assert!(intersect_bed3(&[bed("chr1", 0, 10)], &[bed("chr1", 10, 20)]).is_empty());
        assert!(intersect_bed3(&[bed("chr1", 0, 10)], &[bed("chr2", 0, 10)]).is_empty());
    }

    #[test]
    fn subtract_removes_overlapping_portion_leaving_remainder() {
        let a = vec![bed("chr1", 0, 100)];
        let b = vec![bed("chr1", 20, 40), bed("chr1", 60, 70)];
        let result = subtract_bed3(&a, &b);
        assert_eq!(result, vec![bed("chr1", 0, 20), bed("chr1", 40, 60), bed("chr1", 70, 100)]);
    }

    #[test]
    fn subtract_full_coverage_leaves_nothing() {
        let a = vec![bed("chr1", 10, 20)];
        let b = vec![bed("chr1", 0, 100)];
        assert!(subtract_bed3(&a, &b).is_empty());
    }

    #[test]
    fn subtract_with_no_matching_b_chrom_returns_a_unchanged() {
        let a = vec![bed("chr1", 10, 20)];
        let b = vec![bed("chr2", 0, 100)];
        assert_eq!(subtract_bed3(&a, &b), vec![bed("chr1", 10, 20)]);
    }
}
