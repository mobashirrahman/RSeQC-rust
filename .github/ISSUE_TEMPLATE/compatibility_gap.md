---
name: Compatibility gap
about: This port's output differs from real upstream RSeQC
title: "[command.py] short description"
labels: compatibility
---

**Command and flags**

```
<the exact command line you ran, e.g. bam_stat.py -i sample.bam -q 30>
```

**Input file shape**

Describe the input (or attach/link a minimal reproduction — a tiny synthetic BAM/BED is usually
enough; please don't attach large real patient/sample data).

**Real upstream RSeQC output**

```
<paste here>
```

**This port's actual output**

```
<paste here>
```

**Have you checked `compatibility/divergences.yaml`?**

If this is a *known, disclosed* difference (status `accepted`), it's likely intentional — please
say so and explain why you think it should be revisited. If it's not listed there, it's likely a
genuine bug.

**Upstream source location (optional but very helpful)**

If you can point to the exact line(s) in RSeQC's own source responsible for the behavior, that
speeds up triage a lot.
