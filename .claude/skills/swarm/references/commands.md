# Dispatch commands

All runs use `--pure`. The message positional comes FIRST; `--file` (a yargs
array) goes after, or it swallows the message as a filename. Save returned
session IDs. `--session ID` only for coder repairs; never reuse a coder
session for review. Env keys must be exported first
(`set -a; source .env; set +a` — OpenCode does not load the project `.env`).

```bash
# Template (muse family shown; swap agent + model per table below)
opencode run --pure --agent muse-coder --model opencode/muse-spark-1.3-contributor-free "Implement the attached task" --file /absolute/path/to/handoff.md
opencode run --pure --agent muse-investigator --model opencode/muse-spark-1.3-contributor-free "Answer the attached question" --file /absolute/path/to/question.md
opencode run --pure --agent muse-reviewer --model opencode/muse-spark-1.3-contributor-free "Review the attached changes" --file /absolute/path/to/review.md
```

| Family | `--agent` prefix | `--model` |
|---|---|---|
| muse | `muse-` | `opencode/muse-spark-1.3-contributor-free` |
| glm | `glm-` | `tokenrouter/z-ai/glm-5.3-free` |
| orca | `orca-` | `orcarouter/z-ai/glm-5.3-flash-free` |
| nemo | `nemo-` | `opencode/nemotron-3-ultra-free` |
| mimo | `mimo-` | `opencode/mimo-v2.5-free` |
| gem | `gem-` | `google/gemini-3.5-flash` |
| groq | `groq-` | `groq/openai/gpt-oss-120b` |
| zai | `zai-` | `zai/glm-4.5-flash` |
| lag | `lag-` | `openrouter-free/poolside/laguna-s-2.1:free` |
| north | `north-` | `openrouter-free/cohere/north-mini-code:free` |
| bai | `bai-` | `bai/glm-5.3-flash` |
| tium | `tium-` | `tium/glm-5.3-flash` |
| sail | `sail-` | `sail/zai-org/GLM-5.3-Flash` |
| above | `above-` | `above/glm-5.3-flash-modal` |
| tiyuvta | `tiyuvta-` | `tiyuvta/ornith-ai/ornith-1.5-35b-a3b` |
| nv | `nv-` | `nvidia/nvidia/nemotron-3-ultra-550b-a55b` |
| xcs | `xcs-` | `xkiro/mistralai/codestral-2508` |
| xdv | `xdv-` | `xkiro/mistralai/devstral-medium` |

Suffix is `-coder`, `-investigator`, or `-reviewer`.
