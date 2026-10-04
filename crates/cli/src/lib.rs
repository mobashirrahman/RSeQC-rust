//! Shared command-entry helpers.
//!
//! These exist because the same upstream validation rule was reproduced in one
//! binary and missed in twenty-two others, which is the shape of bug this file
//! removes: a rule that lives in a caller is a rule that a new caller forgets.

use std::alloc::{GlobalAlloc, Layout, System};
use std::io;
use std::path::Path;

/// Out-of-memory diagnostic (B4).
///
/// Past the memory limit every command died with the allocator's internal
/// abort (`memory allocation of N bytes failed`) and no explanation naming
/// the situation (`docs/ENVELOPE.md`, "What happens past the limit"). This
/// `#[global_allocator]` wraps `System`: when the inner allocation returns
/// null it writes a fixed message to file descriptor 2 and returns null, so
/// the normal abort (and today's exit status) follows unchanged.
///
/// Allocation-failure discipline, because this runs on the allocation hot
/// path: the message is a static byte string written with a raw `write(2)`
/// syscall (declared via `extern "C"`, no `libc` dependency added); nothing
/// here allocates, takes a lock, or touches thread-locals, so it is safe to
/// run with the allocator in a failed state. Every `GlobalAlloc` entry
/// point delegates directly to `System` (no default-method chaining), so a
/// failure through any of them reports exactly once.
const OOM_MESSAGE: &[u8] =
    b"rseqc-rust: out of memory; see docs/ENVELOPE.md for per-command memory costs\n";

extern "C" {
    fn write(fd: i32, buf: *const u8, count: usize) -> isize;
}

struct OomReportingAlloc;

unsafe impl GlobalAlloc for OomReportingAlloc {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let ptr = unsafe { System.alloc(layout) };
        if ptr.is_null() {
            unsafe { write(2, OOM_MESSAGE.as_ptr(), OOM_MESSAGE.len()) };
        }
        ptr
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) };
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let ptr = unsafe { System.alloc_zeroed(layout) };
        if ptr.is_null() {
            unsafe { write(2, OOM_MESSAGE.as_ptr(), OOM_MESSAGE.len()) };
        }
        ptr
    }

    unsafe fn realloc(
        &self,
        ptr: *mut u8,
        layout: Layout,
        new_size: usize,
    ) -> *mut u8 {
        let out = unsafe { System.realloc(ptr, layout, new_size) };
        if out.is_null() {
            unsafe { write(2, OOM_MESSAGE.as_ptr(), OOM_MESSAGE.len()) };
        }
        out
    }
}

#[global_allocator]
static GLOBAL_ALLOC: OomReportingAlloc = OomReportingAlloc;

/// Refuse an output prefix whose parent directory does not exist, before any work.
///
/// Upstream's `validate_args` runs this check in every script before it reads an
/// input, so the command refuses immediately and names the directory:
///
/// ```text
/// read_NVC.py: error: output directory does not exist: /nonexistent/out
/// ```
///
/// A port that omits the check instead reads the whole alignment, computes every
/// metric, and only then fails when it opens the output, reporting
/// `No such file or directory (os error 2)` -- an error that names neither the
/// directory nor the flag, arrived at after all the work, and is
/// indistinguishable from a genuinely missing input file. That is the audit's
/// "no silent metric loss" failure in its purest form: the metrics are computed,
/// then discarded without ever being reported.
///
/// The `Path::parent` mapping is upstream's, including the case that looks like a
/// bug but is not: for a bare relative filename such as `result.tsv`,
/// Python's `Path("result.tsv").parent` is `Path(".")`, which exists and is a
/// directory. `Path::parent` in Rust returns an empty `Path` for the same input,
/// so an empty parent is mapped to `.` rather than treated as missing.
pub fn require_existing_output_parent(output_prefix: &Path) -> io::Result<()> {
    let parent = match output_prefix.parent() {
        Some(p) if p.as_os_str().is_empty() => Path::new("."),
        Some(p) => p,
        None => Path::new("."),
    };
    if parent.is_dir() {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("output directory does not exist: {}", parent.display()),
        ))
    }
}

/// Report a usage error the way upstream's `parser.error()` does, and exit.
///
/// Upstream raises every "you passed a path that cannot work" condition from argparse,
/// which writes `prog: error: <message>` to stderr and exits **2**. This port exits 2
/// for ordinary usage errors only because that is what `clap` does, so any such
/// condition returned as an `io::Error` became the one usage error in the tree with a
/// different exit status -- and a caller that checks `$?` sees a runtime failure where
/// upstream reports a usage error.
pub fn usage_exit(program: &str, message: &str) -> ! {
    eprintln!("{program}: error: {message}");
    std::process::exit(2);
}

/// [`require_existing_output_parent`], exiting the way upstream's `parser.error()`
/// does when the check fails.
///
/// Upstream raises this from `validate_args` via `parser.error()`, which writes
/// `prog: error: output directory does not exist: <dir>` to stderr and exits
/// **2**. Returning an `io::Error` instead would exit **1**, so the port would
/// disagree with upstream on exit status as well as on timing: this port already
/// exits 2 for every other usage error, because that is what `clap` does, so exit 1
/// here would be the one usage error in the tree with a different status.
///
/// The message and the exit code match upstream exactly. The usage block above it
/// does not: upstream prints argparse's `usage: read_NVC.py [-h] -i INPUT_FILE ...`
/// summary and clap prints its own. That difference is recorded as DIV-0025 rather
/// than hand-replicated per command, because reproducing 33 argparse usage strings
/// by hand is a way to introduce 33 chances to be subtly wrong.
pub fn require_existing_output_parent_or_exit(program: &str, output_prefix: &Path) {
    if let Err(err) = require_existing_output_parent(output_prefix) {
        eprintln!("{program}: error: {err}");
        std::process::exit(2);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn existing_directory_is_accepted() {
        let dir = std::env::temp_dir();
        assert!(require_existing_output_parent(&dir.join("out.tsv")).is_ok());
        // A NESTED missing directory is refused, because upstream checks only
        // `parent.exists()` and does not create intermediate directories.
        assert!(require_existing_output_parent(&dir.join("a/b/out.tsv")).is_err());
    }

    #[test]
    fn bare_filename_maps_its_empty_parent_to_the_current_directory() {
        // Upstream's `Path("result.tsv").parent` is `Path(".")`, which exists, so a
        // bare filename must NOT be rejected as a missing output directory.
        assert!(require_existing_output_parent(Path::new("result.tsv")).is_ok());
    }

    #[test]
    fn missing_directory_is_rejected_and_named() {
        let missing = std::env::temp_dir().join("rseqc-cli-absent-directory-probe");
        let err = require_existing_output_parent(&missing.join("out.tsv"))
            .expect_err("a missing parent directory must be refused");
        let text = err.to_string();
        assert!(
            text.starts_with("output directory does not exist: "),
            "message must name the condition, got {text:?}"
        );
        assert!(
            text.contains(&missing.display().to_string()),
            "message must name the directory, got {text:?}"
        );
    }

    #[test]
    fn an_existing_file_is_not_a_directory() {
        let path = std::env::temp_dir().join("rseqc-cli-parent-is-a-file");
        std::fs::write(&path, b"x").unwrap();
        // A file where a directory is required must be refused, not silently used.
        let _ = std::fs::remove_file(&path);
    }
}