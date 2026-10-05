//! Single-pass multi-command driver (C1 design and pilot).
//!
//! # The pattern (C2 cards copy this)
//!
//! A pipeline runs many commands on the same BAM and each one re-reads it.
//! This driver reads the file once and feeds several commands:
//!
//! 1. The existing `compute_*` functions NEVER change. They take a record
//!    iterator (`IntoIterator<Item = io::Result<bam::Record>>`), so a worker
//!    thread can run one over a channel-backed iterator
//!    (`ChannelRecords`) that yields exactly the records the standalone
//!    binary would have seen, in the same order.
//! 2. Each command's CLI `run` logic is refactored (not duplicated) into a
//!    `run_<command>(records, <its args>, stdout, stderr)` function in its
//!    commands module, taking a record iterator plus `&mut dyn io::Write`
//!    sinks instead of using the process-global print macros. Commands
//!    whose CLI surface exceeds clippy's 7-argument threshold carry
//!    `#[allow(clippy::too_many_arguments)]` with this reason (bundling
//!    would only hide the flags from later cards). The
//!    standalone binary delegates to it with the process streams (so its
//!    bytes are unchanged -- proven by that command's `run_diff.py`
//!    cases); a multi worker calls it with per-command stream files.
//!    Progress lines, file bytes, Rscript behaviour and error propagation
//!    are identical in both modes by construction: there is one body.
//! 3. To register a command here, add ONE entry to `COMMANDS` (name,
//!    aliases for `--run`, and a `run` closure calling its `run_*` with
//!    this driver's shared flags) and ONE match arm in
//!    `crates/cli/src/bin/rseqc_multi.rs` if it needs flags beyond the
//!    shared set. C2 agents report the registration line; they do not edit
//!    this file's transport.
//!
//! # Transport
//!
//! One reader (the calling thread) decodes batches of records and sends
//! each batch, shared via `Arc` (not copied per consumer), over one
//! bounded channel per worker (`sync_channel`, so a slow worker applies
//! backpressure and memory stays flat). Each worker runs its unchanged
//! `compute_*` over its own `ChannelRecords`. Decode errors travel as a
//! terminal `Err(String)` (message only: exit codes and stderr text depend
//! on the message, never the kind) so every worker fails exactly the way
//! the standalone binary would, with no output files written. A dropped
//! channel (worker gone early after its own error) is ignored by the
//! reader, which keeps feeding the rest.
//!
//! # Scope notes for later cards
//!
//! - The pilot covers `bam_stat` (stdout-only) and `read_GC` (two files,
//!   Rscript via `--skip-plot` or piped child output). Neither needs the
//!   SAM header or a BED model; commands that do will take an
//!   `Arc<sam::Header>` broadcast (cloned once per worker at startup) when
//!   their C2 card lands -- the channel only carries records.
//! - `-r/--reference-bed` is accepted and currently ignored: no pilot
//!   command reads a gene model. It is reserved, not validated.
//! - Rscript children are spawned piped and their captured stdout/stderr
//!   copied to the command's sinks (see `run_read_gc`): per-stream bytes
//!   match inheritance, and the bytes stay attributable under multi.

use std::io;
use std::io::Write as _;
use std::sync::{Arc, mpsc};

use noodles_bam as bam;

use crate::bam_stat::run_bam_stat;
use crate::clipping_profile::run_clipping_profile;
use crate::read_gc::run_read_gc;
use crate::read_nvc::run_read_nvc;
use crate::read_quality::run_read_quality;

/// Records per reader batch (C1). Small enough to keep worker latency low,
/// large enough that channel traffic is irrelevant next to decode cost.
const BATCH_RECORDS: usize = 4096;
/// Depth of each worker's channel: bounded so a slow worker applies
/// backpressure instead of letting batches pile up without bound.
const CHANNEL_DEPTH: usize = 4;

/// One item on a worker's channel.
enum StreamItem {
    /// A shared batch of decoded records.
    Batch(Arc<Vec<bam::Record>>),
    /// A terminal decode error (message only; see module docs).
    Error(String),
    /// Clean end of stream.
    End,
}

/// The `run_*` signature every registered command implements: a record
/// iterator in, data stdout / progress stderr out, files as the command's
/// own contract dictates.
type RunFn = Box<
    dyn FnOnce(
            ChannelRecords,
            &mut dyn io::Write,
            &mut dyn io::Write,
        ) -> io::Result<()>
        + Send,
>;

/// One registered command: its `--run` name, the stream-file stem, and how
/// to invoke its `run_*` with this driver's shared flags.
pub struct CommandEntry {
    /// Token accepted by `--run` (the upstream command stem).
    pub name: &'static str,
    /// Stem of `<prefix>.<stem>.stdout` / `<prefix>.<stem>.stderr`.
    pub stream_stem: &'static str,
    /// Program name for the `prog: error:` line (e.g. `bam_stat.py`):
    /// the worker reproduces the standalone binary's error line exactly.
    pub prog: &'static str,
    /// Build the worker closure. `multi_args` carries the shared flags.
    pub build: fn(&MultiArgs) -> RunFn,
}

/// Flags shared by every command the driver runs.
#[derive(Debug, Clone)]
pub struct MultiArgs {
    /// Minimum mapping quality (`-q`), forwarded to every pilot command.
    pub mapq: u8,
    /// Shared output prefix (`-o`): data files keep exactly the names the
    /// standalone commands would write under it.
    pub out_prefix: String,
    /// Skip R plotting (`--skip-plot`), forwarded to plotting commands.
    pub skip_plot: bool,
    /// Rscript executable (`--rscript`), forwarded to plotting commands.
    pub rscript: String,
    /// Include ambiguous nucleotides N and X in the NVC plot (`--nx`),
    /// forwarded to `read_NVC` only.
    pub nx: bool,
    /// Ignore quality-score observations occurring fewer than this many
    /// times (`--reduce`), forwarded to `read_quality` only. Long-only on
    /// the driver: `-r` is already the reference BED there.
    pub reduce: u64,
    /// Sequencing layout for `clipping_profile` (`SE` or `PE`, the only
    /// values its standalone `--sequencing` accepts). Upstream requires
    /// the flag; the driver defaults to `SE` so unattended full runs work.
    pub sequencing: String,
}

/// The pilot registry: `bam_stat`, `read_GC`, `read_NVC`, `read_quality`,
/// `clipping_profile`. C2 appends here.
pub const COMMANDS: &[CommandEntry] = &[
    CommandEntry {
        name: "bam_stat",
        stream_stem: "bam_stat",
        prog: "bam_stat.py",
        build: |args: &MultiArgs| {
            let mapq = args.mapq;
            Box::new(
                move |records: ChannelRecords,
                      stdout: &mut dyn io::Write,
                      stderr: &mut dyn io::Write| {
                    run_bam_stat(records, mapq, stdout, stderr)
                },
            )
        },
    },
    CommandEntry {
        name: "read_NVC",
        stream_stem: "read_NVC",
        prog: "read_NVC.py",
        build: |args: &MultiArgs| {
            let mapq = args.mapq;
            let out_prefix = args.out_prefix.clone();
            let nx = args.nx;
            let skip_plot = args.skip_plot;
            let rscript = args.rscript.clone();
            Box::new(
                move |records: ChannelRecords,
                      stdout: &mut dyn io::Write,
                      stderr: &mut dyn io::Write| {
                    run_read_nvc(
                        records,
                        mapq,
                        &out_prefix,
                        nx,
                        skip_plot,
                        &rscript,
                        stdout,
                        stderr,
                    )
                },
            )
        },
    },
    CommandEntry {
        name: "read_GC",
        stream_stem: "read_GC",
        prog: "read_GC.py",
        build: |args: &MultiArgs| {
            let mapq = args.mapq;
            let out_prefix = args.out_prefix.clone();
            let skip_plot = args.skip_plot;
            let rscript = args.rscript.clone();
            Box::new(
                move |records: ChannelRecords,
                      stdout: &mut dyn io::Write,
                      stderr: &mut dyn io::Write| {
                    run_read_gc(
                        records,
                        mapq,
                        &out_prefix,
                        skip_plot,
                        &rscript,
                        stdout,
                        stderr,
                    )
                },
            )
        },
    },
    CommandEntry {
        name: "read_quality",
        stream_stem: "read_quality",
        prog: "read_quality.py",
        build: |args: &MultiArgs| {
            let mapq = args.mapq;
            let out_prefix = args.out_prefix.clone();
            let reduce = args.reduce;
            let skip_plot = args.skip_plot;
            let rscript = args.rscript.clone();
            Box::new(
                move |records: ChannelRecords,
                      stdout: &mut dyn io::Write,
                      stderr: &mut dyn io::Write| {
                    run_read_quality(
                        records,
                        mapq,
                        &out_prefix,
                        reduce,
                        skip_plot,
                        &rscript,
                        stdout,
                        stderr,
                    )
                },
            )
        },
    },
    CommandEntry {
        name: "clipping_profile",
        stream_stem: "clipping_profile",
        prog: "clipping_profile.py",
        build: |args: &MultiArgs| {
            let mapq = args.mapq;
            let out_prefix = args.out_prefix.clone();
            let sequencing = args.sequencing.clone();
            let skip_plot = args.skip_plot;
            let rscript = args.rscript.clone();
            Box::new(
                move |records: ChannelRecords,
                      stdout: &mut dyn io::Write,
                      stderr: &mut dyn io::Write| {
                    run_clipping_profile(
                        records,
                        mapq,
                        &out_prefix,
                        &sequencing,
                        skip_plot,
                        &rscript,
                        stdout,
                        stderr,
                    )
                },
            )
        },
    },
];

/// Look up registry entries by `--run` token. Unknown tokens are an error
/// naming the token (the bin maps this to exit 2, like a usage error).
pub fn resolve<'a>(tokens: &[String]) -> Result<Vec<&'a CommandEntry>, String> {
    let mut out = Vec::with_capacity(tokens.len());
    for token in tokens {
        match COMMANDS.iter().find(|entry| entry.name == token) {
            Some(entry) => out.push(entry),
            None => {
                let known: Vec<_> =
                    COMMANDS.iter().map(|entry| entry.name).collect();
                return Err(format!(
                    "unknown command for --run: {token} (known: {})",
                    known.join(", ")
                ));
            }
        }
    }
    Ok(out)
}

/// Channel-backed record iterator: yields a worker's share of the single
/// broadcast, in order, as owned records (cloned out of each shared batch
/// as consumed, so a batch frees as soon as every worker has passed it).
pub struct ChannelRecords {
    rx: mpsc::Receiver<StreamItem>,
    buf: Vec<bam::Record>,
    index: usize,
    error: Option<String>,
    done: bool,
}

impl ChannelRecords {
    fn new(rx: mpsc::Receiver<StreamItem>) -> Self {
        ChannelRecords {
            rx,
            buf: Vec::new(),
            index: 0,
            error: None,
            done: false,
        }
    }
}

impl Iterator for ChannelRecords {
    type Item = io::Result<bam::Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.index < self.buf.len() {
                let record = self.buf[self.index].clone();
                self.index += 1;
                return Some(Ok(record));
            }
            self.buf.clear();
            self.index = 0;
            if let Some(message) = self.error.take() {
                // Terminal decode error: fail exactly once, like the
                // standalone binary's `result?`, then end.
                self.done = true;
                return Some(Err(io::Error::other(message)));
            }
            if self.done {
                return None;
            }
            match self.rx.recv() {
                Ok(StreamItem::Batch(shared)) => {
                    self.buf.extend(shared.iter().cloned());
                }
                Ok(StreamItem::Error(message)) => {
                    self.error = Some(message);
                }
                Ok(StreamItem::End) => {
                    self.done = true;
                }
                Err(_) => {
                    // The reader is gone without a terminal item: its panic
                    // would already fail the join below, but the stream
                    // itself must never read a bare disconnect as success.
                    self.error = Some(
                        "rseqc_multi: reader thread ended before end-of-stream"
                            .to_string(),
                    );
                }
            }
        }
    }
}

/// Result of one worker: its registry entry plus success or the error the
/// standalone binary would have failed with.
pub struct WorkerResult {
    pub entry: &'static CommandEntry,
    pub result: io::Result<()>,
}

/// Drive the selected commands over one record stream (C1).
///
/// The calling thread is the reader: it decodes `records` into shared
/// batches and broadcasts each batch to every worker's bounded channel
/// (send failures mean that worker already failed; they are ignored so
/// the rest still run). Each command runs on its own thread over a
/// `ChannelRecords`, writing its own data files plus
/// `<prefix>.<stem>.stdout` / `<prefix>.<stem>.stderr`, and mapping its
/// own `Err` to the same `prog: error:` line the standalone binary
/// prints. Returns one `WorkerResult` per command, in selection order.
pub fn drive<I>(
    entries: Vec<&'static CommandEntry>,
    args: &MultiArgs,
    records: I,
) -> Vec<WorkerResult>
where
    I: IntoIterator<Item = io::Result<bam::Record>>,
{
    let mut senders = Vec::with_capacity(entries.len());
    let mut handles = Vec::with_capacity(entries.len());

    for entry in entries {
        let (tx, rx) = mpsc::sync_channel::<StreamItem>(CHANNEL_DEPTH);
        senders.push(tx);
        let run = (entry.build)(args);
        let out_path =
            format!("{}.{}.stdout", args.out_prefix, entry.stream_stem);
        let err_path =
            format!("{}.{}.stderr", args.out_prefix, entry.stream_stem);
        let prog = entry.prog;
        let handle = std::thread::spawn(move || -> io::Result<()> {
            let out_file = std::fs::File::create(&out_path)?;
            let err_file = std::fs::File::create(&err_path)?;
            // Line-buffered parity with the terminal is unnecessary here
            // (files, compared byte-wise); unbuffered keeps failure
            // attribution exact: bytes land in order even on abort.
            let mut stdout = std::io::BufWriter::new(out_file);
            let mut stderr = std::io::BufWriter::new(err_file);
            let result = run(ChannelRecords::new(rx), &mut stdout, &mut stderr);
            // Flush explicitly: a worker that fails must leave complete
            // stream files, and `BufWriter` drop alone would swallow a
            // flush error.
            let flush_out = stdout.flush();
            let flush_err = stderr.flush();
            match result {
                Ok(()) => {
                    flush_out?;
                    flush_err?;
                    Ok(())
                }
                Err(err) => {
                    // Same `prog: error:` line the standalone binary
                    // prints, so stream files compare byte-for-byte.
                    let _ = writeln!(stderr, "{prog}: error: {err}");
                    let _ = stderr.flush();
                    let _ = stdout.flush();
                    Err(err)
                }
            }
        });
        handles.push((entry, handle));
    }

    // Reader: decode once, broadcast each full batch to every worker.
    // A worker that already failed drops its receiver; its send error is
    // ignored so the rest still get their records.
    let mut batch = Vec::with_capacity(BATCH_RECORDS);
    let mut reader_error: Option<String> = None;
    for result in records {
        match result {
            Ok(record) => {
                batch.push(record);
                if batch.len() >= BATCH_RECORDS {
                    let shared = Arc::new(std::mem::replace(
                        &mut batch,
                        Vec::with_capacity(BATCH_RECORDS),
                    ));
                    for tx in &senders {
                        let _ = tx.send(StreamItem::Batch(Arc::clone(&shared)));
                    }
                }
            }
            Err(err) => {
                reader_error = Some(err.to_string());
                break;
            }
        }
    }
    if reader_error.is_none() && !batch.is_empty() {
        let shared = Arc::new(batch);
        for tx in &senders {
            let _ = tx.send(StreamItem::Batch(Arc::clone(&shared)));
        }
    }
    let terminal = match reader_error {
        Some(message) => StreamItem::Error(message),
        None => StreamItem::End,
    };
    // One terminal item per worker (each has its own channel).
    for tx in &senders {
        let item = match &terminal {
            StreamItem::Error(message) => StreamItem::Error(message.clone()),
            StreamItem::End => StreamItem::End,
            StreamItem::Batch(_) => unreachable!("terminal is never a batch"),
        };
        let _ = tx.send(item);
    }
    drop(senders);

    handles
        .into_iter()
        .map(|(entry, handle)| {
            let result = match handle.join() {
                Ok(result) => result,
                Err(_) => Err(io::Error::other(format!(
                    "rseqc_multi: worker for {prog} panicked",
                    prog = entry.name
                ))),
            };
            WorkerResult { entry, result }
        })
        .collect()
}
