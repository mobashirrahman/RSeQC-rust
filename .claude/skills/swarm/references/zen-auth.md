# OpenCode Zen auth setup

Zen backs three families: `muse-*`, `nemo-*`, `mimo-*`. Zen is optional —
everything else is key-based and unaffected.

## Walkthrough (read to the user, step by step)

1. Get the key **once** (it works on all servers): sign in at
   `https://opencode.ai/auth`, add billing details (required even for $0
   use), copy the API key.
2. On the server: run `opencode auth login` and paste the key (headless-safe,
   device-code flow), or in the TUI run `/connect` → select OpenCode Zen →
   paste.
3. Verify with a 10-second probe:
   `opencode run --pure --model opencode/mimo-v2.5-free "say OK"`.

## Failure signature

401/auth errors on any `opencode/*` model = Zen login missing, not a bug.
Redo step 2. Never go hunting credential files — token location varies by
install and is none of the coordinator's business.

## Caveats

- Zen free models are limited-time and may change; re-check
  `https://opencode.ai/zen/v1/models` and the pricing page periodically.
- Free-tier data terms apply (prompts may be used for training/logging).
  Delegates already receive bounded technical packets only — keep it that way.
