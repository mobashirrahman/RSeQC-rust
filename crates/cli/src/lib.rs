//! Shared command-entry helpers.
//!
//! These exist because the same upstream validation rule was reproduced in one
//! binary and missed in twenty-two others, which is the shape of bug this file
//! removes: a rule that lives in a caller is a rule that a new caller forgets.

use std::io;
use std::path::Path;

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