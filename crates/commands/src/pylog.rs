//! Minimal replication of Python's
//! `logging.basicConfig(format="%(asctime)s [%(levelname)s] %(message)s",
//! datefmt="%Y-%m-%d %H:%M:%S")`, used by all four `sc_*` single-cell
//! commands (`sc_bamStat.py`, `sc_editMatrix.py`, `sc_seqQual.py`,
//! `sc_seqLogo.py`; see each script's `configure_logging()`).
//!
//! Python's default `asctime` formatter uses **local** time
//! (`time.localtime`), not UTC. That differs from this port's other
//! timestamped-progress helper (`FPKM-UQ.py`/`geneBody_coverage.py`'s
//! `@ <TS>: msg` `printlog()`, in `crates/cli/src/bin/fpkm_uq.rs`), which
//! deliberately uses UTC because that line is never compared byte-for-byte
//! against upstream. The `sc_*` commands' `<TS> [LEVEL] msg` lines ARE
//! compared, via `verification/synthetic_sweep.py`'s `LOG_PREFIX_RE`,
//! which strips the whole prefix (both sides) before comparing -- so the
//! *wall-clock correctness* of the timestamp doesn't affect the diff
//! either (any digit string matching `TIMESTAMP_RE` works), but getting
//! local time right here costs nothing and matches upstream's own
//! behavior rather than an unrelated simplification.
//!
//! Uses a direct `localtime_r(3)` FFI call rather than pulling in a
//! chrono/time dependency for four call sites.

use std::ffi::c_char;
use std::time::{SystemTime, UNIX_EPOCH};

#[repr(C)]
struct Tm {
    tm_sec: i32,
    tm_min: i32,
    tm_hour: i32,
    tm_mday: i32,
    tm_mon: i32,
    tm_year: i32,
    tm_wday: i32,
    tm_yday: i32,
    tm_isdst: i32,
    tm_gmtoff: i64,
    tm_zone: *const c_char,
}

extern "C" {
    fn localtime_r(timep: *const i64, result: *mut Tm) -> *mut Tm;
}

/// `%Y-%m-%d %H:%M:%S` local-time timestamp, matching Python
/// `logging.Formatter`'s default `asctime` converter.
fn local_timestamp() -> String {
    let secs = SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_secs() as i64).unwrap_or(0);
    // SAFETY: `tm` is a plain-old-data out-parameter fully written by
    // `localtime_r` before any field is read; `secs` is a valid `time_t`.
    let tm: Tm = unsafe {
        let mut tm: Tm = std::mem::zeroed();
        localtime_r(&secs, &mut tm);
        tm
    };
    format!(
        "{:04}-{:02}-{:02} {:02}:{:02}:{:02}",
        tm.tm_year + 1900,
        tm.tm_mon + 1,
        tm.tm_mday,
        tm.tm_hour,
        tm.tm_min,
        tm.tm_sec
    )
}

/// One `logging.<level>(...)` line as upstream's `basicConfig` format
/// renders it: `<local-ts> [<LEVEL>] <msg>`.
pub fn log_line(level: &str, msg: &str) -> String {
    format!("{} [{level}] {msg}", local_timestamp())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_line_matches_python_logging_format() {
        let line = log_line("INFO", "hello world");
        // "YYYY-MM-DD HH:MM:SS [INFO] hello world"
        let re_parts: Vec<&str> = line.splitn(2, " [INFO] ").collect();
        assert_eq!(re_parts.len(), 2);
        assert_eq!(re_parts[1], "hello world");
        let ts = re_parts[0];
        assert_eq!(ts.len(), 19, "timestamp {ts:?} should be YYYY-MM-DD HH:MM:SS");
        assert_eq!(ts.as_bytes()[4], b'-');
        assert_eq!(ts.as_bytes()[7], b'-');
        assert_eq!(ts.as_bytes()[10], b' ');
        assert_eq!(ts.as_bytes()[13], b':');
        assert_eq!(ts.as_bytes()[16], b':');
    }
}
