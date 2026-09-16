# Guardrails

- **Routing.** Delegates get anything with a machine-checkable acceptance
  criterion: investigation, tests, lint, mechanical refactors, draft tests.
  Claude keeps judgment: architecture, security, final pre-commit validation,
  and review of a delegate's own output.
- **One worker per backend**, disjoint file scopes. No competing solutions,
  no recursive delegation (delegates never spawn teams).
- **Never silently substitute models.** Never add a model to a provider block
  or pass one on the command line without confirming its price first. Free
  today ≠ free tomorrow — re-verify periodically.
- **Paid bans.** The paid OpenRouter route is banned; only the two pinned
  `:free` slugs on `OPENROUTER_FREE_API_KEY`. Only `glm-5.3-flash` on the
  `bai` key (metered). Only the pinned slug on each router block — never
  `orcarouter/auto`. Only the two pinned `mistralai/` slugs on the `xkiro`
  key — never a bare model name (xKiro 404s it) and never a paid/premium
  slug on the $0 allowance.
- **Metered backends last.** `bai`, `tium`, `sail`, `above`, `tiyuvta`, `nv`
  are funded, not free — prefer the twelve free backends first, check balances
  before large batches, never describe metered as free. `sail` has zero data
  retention: prefer it for anything sensitive.
- **Built-in-only rule.** Gemini and Groq run through OpenCode's built-in
  `google`/`groq` providers. Routing either through a generic
  OpenAI-compatible block breaks multi-step tool use (verified: Gemini
  rejects calls missing `thought_signature`; Groq rejects returned
  `reasoning_content`).
- **Same-model review shares blind spots.** A backend's reviewer is an extra
  pass only, never the gate. The gating review runs on Claude.
- **Validate before commit.** A free model attesting to its own work is not
  evidence. Coordinator re-runs or inspects the checks that gate a commit.
- **Data terms.** Free endpoints may log prompts. Delegates receive bounded
  technical packets only — never credentials, keys, or bulk project data.
