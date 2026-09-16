# CLI workers — dispatch reference

CLI workers are a second delegate tier alongside the OpenCode backends. Same
principle (work runs on someone else's quota; Claude reads a bounded report),
different transport: a subscription/gateway CLI invoked via `Bash` in a git
worktree instead of `opencode run` against an API key.

Use these **in addition to** `references/backends.md` / `references/commands.md`,
not instead. Read `references/guardrails.md` and `references/quota.md` first —
every rule there applies here too.

## Roster

| worker | pool | wired roles | never |
|---|---|---|---|
| `kilo` | Kilo Gateway `:free` catalogue (balance $0) | coder, investigator | config edits, commits, paid model slugs |
| `cline` | Cline promo free models (`meta/muse-spark-1.3`) | reviewer | the gating review; bulk coding; concurrent runs with mixed roles¹ |
| `cursor` | Cursor Hobby Auto | critic (short reviews) | coder (wrapper refuses it) |

¹ cline persists plan/act mode globally; `run-cline.sh` sets it per-role before
each run, so two cline workers with different roles must not run at once. The
swarm only dispatches cline as `reviewer`, so this doesn't bite in practice.

Status and per-worker auth: `workers/<w>/AUTH.md`. A worker not marked WIRED
there is not available — fall through to OpenCode backends or direct-API
delegates.

## Wrapper commands

Run from the repo root (the wrapper `cd`s into the worktree itself):

```bash
workers/kilo/run-kilo.sh     coder        .worktrees/kilo-parser   /abs/task.md   900
workers/kilo/run-kilo.sh     investigator .worktrees/kilo-probe     /abs/q.md      300
workers/cline/run-cline.sh   reviewer     .worktrees/cline-review   /abs/review.md 600
workers/cursor/run-cursor.sh reviewer     .worktrees/cursor-critic  /abs/critic.md 300
```

Under the hood (all verified 2026-09-09):
- **kilo** (OpenCode fork v7.5.16): `kilo run --pure --auto --format json --dir
  <cwd> --agent kilo-<role> --model $KILO_MODEL -f <task>`; evidence via `kilo
  export --sanitize`. `KILO_MODEL` must be a `kilo/…:free` slug (balance $0).
- **cline** v3.0.61: `cline --json --auto-approve true -m $CLINE_MODEL <task>`
  (reviewer adds `-p`). NDJSON → normalized `{events,final,exit,degraded}`.
- **cursor** 2026.09.08: `cursor-agent -p --output-format json [--plan] <task>`.

Each prints one JSON line to stdout:

```json
{"worker":"kilo","role":"coder","cwd":"...","exit":0,"status":"ok",
 "log_file":"~/.cache/swarm-workers/kilo-coder-....log",
 "result_file":"...json","final":"<=8k of the worker's report"}
```

Parse `status` + `final`; open `log_file` / `result_file` only when you need
detail. Never paste a whole log into context.

## Worktree protocol

The coordinator owns all git operations. Workers never touch git.

1. **Create** one worktree per task, off the integration branch:
   `git worktree add .worktrees/<worker>-<task> <base-branch>`
2. **Dispatch** exactly one worker into that worktree (`cwd` arg). One worker
   per worktree, disjoint file scopes — never shard one task across workers,
   never point two workers at the same worktree.
3. **Workers never commit.** Extends the delegates-don't-commit rule. The
   wrapper's edit gate proves the worker leaves the tree uncommitted.
4. **Validate** yourself: `cd` into the worktree, run the tests/lint that gate
   the merge. A worker attesting to its own work is not evidence.
5. **Integrate**: `git -C .worktrees/<...> diff` → judge → apply to the
   integration branch (`git merge`, or cherry-pick / `git -C ... diff | git apply`
   for a partial take). Competing solutions: diff both, keep one, discard the
   rest.
6. **Remove**: `git worktree remove .worktrees/<worker>-<task>` (add `--force`
   if the worker left it dirty, which is expected).

`.worktrees/` should be gitignored in the target repo (setup adds it).

## Role split (model diversity beats same-model repetition)

- **Claude** — architect + integrator. Owns the plan, the worktree merges, the
  gating review, all commits.
- **kilo / OpenCode backends** — coders. Bounded implementation tasks with a
  machine-checkable acceptance check.
- **cline** — reviewer. Feed it the scoped diff + your test evidence as an
  NDJSON run; its structured findings are an *extra* pass.
- **cursor Auto** — critic. One short "what's the worst assumption here"
  review per integration, pool permitting. Skip when near-empty.

Same-model review shares the coder's blind spots — the gating review stays on
Claude regardless of how many workers looked at it.

## Failure handling

| signal | wrapper status | coordinator action |
|---|---|---|
| exit 124 | `timeout` | worker already killed; re-dispatch the task to a different worker/backend |
| exit 1 | `rejected` | shrink the scope, retry **once** on the same worker; still failing → reroute |
| exit 2 | `error` | setup/wiring bug — fix wrapper or auth, not a task retry |
| parse-degraded (cline) | `ok` + `degraded:true` | treat the report as low-confidence; don't let it gate anything |
| pool exhausted (see ledger) | — | don't dispatch; fall back to OpenCode free backends, then direct-API delegates |

Reroute order for a coding task: kilo → OpenCode free backends → OpenCode
metered → direct-API delegate. Never spin up a paid path silently.

## Quota ledger (check before dispatch, never after failing)

| pool | how to check | where recorded |
|---|---|---|
| Kilo Gateway | `kilo stats --days N` / `kilo profile` | `workers/kilo/AUTH.md` table |
| Cline promo | Cline account balance / `cline` status view | `workers/cline/AUTH.md` |
| Cursor Hobby | remaining shown in `cursor-agent` output / dashboard | `workers/cursor/AUTH.md` |
| OpenRouter free (shared 50/day) | unchanged — `references/quota.md` | — |

Before a batch: glance at the relevant `AUTH.md` baseline + rough spend since.
If a pool looks tight, route elsewhere up front. Update the `AUTH.md` table
when you observe a real number.

## Guardrails (delta from `references/guardrails.md`)

- Worker auth (Kilo Gateway, Cline account, Cursor Hobby) lives in each CLI's
  own credential store. Never in this repo, `.env`, a config file, or chat.
- Per-provider API keys used by `kilo.jsonc` follow the existing `{env:VAR}` +
  global-scope rule. Project-level `kilo.json` ignores `{env:}` — don't use it.
- Free CLI pools may log prompts: workers get bounded technical packets only,
  never credentials or bulk data — identical to the API backends.
- Never add a model to `kilo.jsonc` without confirming its price. The pinned
  slugs match `opencode/opencode.json`.
- `cursor` is critic-only. The wrapper enforces it; don't override.
