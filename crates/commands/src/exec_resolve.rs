//! Shared `shutil.which` port, used by every command that shells out to
//! an external tool (`Rscript` for the 16 R-plotting commands,
//! `htseq-count` for `FPKM-UQ.py`) to pre-check resolvability BEFORE
//! attempting to spawn it -- matching upstream's own behavior of a
//! friendly, content-specific "<tool> executable not found: <name>"
//! error rather than letting the OS's raw spawn-failure error surface
//! (which upstream never does: it always resolves first via
//! `shutil.which`). This distinction is invisible in this project's own
//! development sandbox, where all these tools happen to be installed,
//! but matters for a genuine "clean-room" install lacking Python/R --
//! see `docs/PORTING_PLAN.md` Step 10 and the README's own disclosed
//! "not tested in a clean-room environment" limitation.
//!
//! `FPKM-UQ.py`'s own `fpkm_uq::resolve_executable` predates this module
//! and has its own message text (`cannot find htseq-count executable
//! "<name>"`) baked in -- not migrated here, since it's a single call
//! site with its own established, already-correct wording; this module
//! is for the 16 sites that all need the SAME underlying resolution
//! logic with call-site-specific message text layered on top.
use std::path::{Path, PathBuf};

/// Ports `shutil.which(cmd)`: if `cmd` contains a path separator, check
/// that exact path directly (no `PATH` search); otherwise search each
/// directory in `$PATH` in order for the first executable match.
/// POSIX semantics (any executable bit set); on non-Unix platforms,
/// falls back to plain existence (no `PATHEXT` resolution).
pub fn which(cmd: &str) -> Option<PathBuf> {
    if cmd.contains('/') || cmd.contains(std::path::MAIN_SEPARATOR) {
        let p = Path::new(cmd);
        return if is_executable_file(p) { Some(p.to_path_buf()) } else { None };
    }

    let path_var = std::env::var_os("PATH")?;
    for dir in std::env::split_paths(&path_var) {
        let candidate = dir.join(cmd);
        if is_executable_file(&candidate) {
            return Some(candidate);
        }
    }
    None
}

#[cfg(unix)]
fn is_executable_file(path: &Path) -> bool {
    use std::os::unix::fs::PermissionsExt;
    std::fs::metadata(path).map(|m| m.is_file() && m.permissions().mode() & 0o111 != 0).unwrap_or(false)
}

#[cfg(not(unix))]
fn is_executable_file(path: &Path) -> bool {
    path.is_file()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn which_finds_a_real_system_binary_by_bare_name() {
        // `sh` is present on every POSIX CI runner this project targets.
        assert!(which("sh").is_some());
    }

    #[test]
    fn which_returns_none_for_an_unresolvable_bare_name() {
        assert!(which("definitely-not-a-real-executable-name-xyz").is_none());
    }

    #[test]
    fn which_checks_an_explicit_path_directly_without_searching_path() {
        assert!(which("/bin/sh").is_some() || which("/usr/bin/sh").is_some());
        assert!(which("/definitely/not/a/real/path/xyz").is_none());
    }

    #[test]
    fn which_rejects_a_non_executable_file() {
        let dir = std::env::temp_dir().join(format!("exec_resolve_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let f = dir.join("not_executable");
        std::fs::write(&f, "not a script").unwrap();
        assert!(which(f.to_str().unwrap()).is_none());
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
