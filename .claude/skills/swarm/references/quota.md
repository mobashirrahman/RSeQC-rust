# Quota discipline (mandatory)

~98% of Claude tokens in a session go to re-reading history: every
coordinator turn replays full context. Codex cost ≈ reasoning time. These
rules are not optional.

## Claude

- **Batch per turn.** Each turn replays everything — three questions in one
  turn costs ~⅓ of three separate turns. Never drip-feed a delegate's
  follow-ups one per turn.
- **`/compact` before context balloons**, `/clear` between tasks. Neither
  restores quota; both stop per-turn replay cost compounding. Write a 5-line
  handoff note (objective, files touched, pending, session IDs) before
  clearing.
- **Delegate-first.** Anything machine-checkable goes to a free backend on
  the first pass, not after you've already read the code yourself.
- **Right model.** Coordinator on Sonnet; architect override to Opus only
  when judgment genuinely needs it; utility subagents on Haiku
  (`CLAUDE_CODE_SUBAGENT_MODEL=haiku`). Never Fable for routine work.
- **Window discipline.** The 5h clock starts at your first prompt and is one
  pool shared with web/desktop ChatGPT. Start heavy sessions right after
  reset; don't run parallel Claude sessions in one window; keep casual chat
  off the pool on heavy days.
- **Measure.** `/usage` before long sessions; JSONL-based spend-by-model
  reports when diagnosing a drain.

## Codex

- Global `model_reasoning_effort = "medium"`; per-thread high/max only for
  genuinely hard problems. Check `/status` before big runs.
- Near quota: finish-and-validate what exists, then stop with a resume
  handoff — never start a mutation you may not be able to finish and verify.
