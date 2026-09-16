# Backend roster

`coder` = edit+bash. `investigator` = bash only, read-only. `reviewer` =
neither edit nor bash. Full permission blocks live in `opencode/agents/`.

| Agents | Model | Credential | Cost basis |
|---|---|---|---|
| `muse-*` | `opencode/muse-spark-1.3-contributor-free` | OpenCode Zen login | free, limited time |
| `glm-*` | `tokenrouter/z-ai/glm-5.3-free` | `TOKENROUTER_API_KEY` | $0, this slug only |
| `orca-*` | `orcarouter/z-ai/glm-5.3-flash-free` | `ORCAROUTER_API_KEY` | $0, never `orcarouter/auto` |
| `nemo-*` | `opencode/nemotron-3-ultra-free` | OpenCode Zen login | free, limited time |
| `mimo-*` | `opencode/mimo-v2.5-free` | OpenCode Zen login | free, limited time |
| `gem-*` | `google/gemini-3.5-flash` | `GOOGLE_GENERATIVE_AI_API_KEY` (built-in `google` provider) | free tier |
| `groq-*` | `groq/openai/gpt-oss-120b` | `GROQ_API_KEY` (built-in `groq` provider) | free dev tier |
| `zai-*` | `zai/glm-4.5-flash` | `ZAI_API_KEY` (custom `zai` block) | $0/$0, this slug only |
| `lag-*` | `openrouter-free/poolside/laguna-s-2.1:free` | `OPENROUTER_FREE_API_KEY` | coding agent; ONE shared 50/day pool with `north-*` |
| `north-*` | `openrouter-free/cohere/north-mini-code:free` | `OPENROUTER_FREE_API_KEY` | coding agent; ONE shared 50/day pool with `lag-*` |
| `bai-*` | `bai/glm-5.3-flash` | `BAI_API_KEY` | METERED ~$0.075/$0.25 per M — use last |
| `tium-*` | `tium/glm-5.3-flash` | `TIUM_API_KEY` | METERED on a weighted-token balance (per-call cost headers) |
| `sail-*` | `sail/zai-org/GLM-5.3-Flash` | `SAIL_API_KEY` | METERED paygo — zero data retention, prefer for anything sensitive |
| `above-*` | `above/glm-5.3-flash-modal` | `ABOVE_API_KEY` | $10 credit, $2/day cap on the free key |
| `tiyuvta-*` | `tiyuvta/ornith-ai/ornith-1.5-35b-a3b` | `TIYUVTA_API_KEY` | METERED — only non-GLM-family model in the swarm |
| `nv-*` | `nvidia/nvidia/nemotron-3-ultra-550b-a55b` | `NVIDIA_API_KEY` (custom `nvidia` block) | METERED paygo (stutter intended: block id + NVIDIA org prefix) |
| `xcs-*` | `xkiro/mistralai/codestral-2508` | `XKIRO_API_KEY` (custom `xkiro` block) | $0 free-tier, this slug only; shares 5M/day pool with `xdv-*` |
| `xdv-*` | `xkiro/mistralai/devstral-medium` | `XKIRO_API_KEY` (custom `xkiro` block) | $0 free-tier, this slug only; shares 5M/day pool with `xcs-*` |

Per-backend notes: `glm-*` has a 1M-token window — prefer it for tasks that
read broadly in one pass. `groq-*` is fastest per token. `gem-*` has the
deepest context on the free tier. `lag-*`/`north-*` are coding-specialized;
spend their shared 50/day pool on coding tasks only.
