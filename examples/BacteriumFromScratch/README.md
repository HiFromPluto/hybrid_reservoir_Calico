# Bacterium from scratch

Modular *E. coli* chassis: one cell, then a colony, then a population.
Not a reservoir, not a chip layout.

- Job 1: [`CITATION_DOSSIER.md`](CITATION_DOSSIER.md)
- Job 2: [`PROTOCOL.md`](PROTOCOL.md) — isolated rod **PASS**
  ([`JOB2_STANDING.md`](JOB2_STANDING.md))
- Job 3: [`PROTOCOL.md`](PROTOCOL.md) — Hertzian packing **PASS**
  ([`JOB3_STANDING.md`](JOB3_STANDING.md))
- Job 3b: [`PROTOCOL.md`](PROTOCOL.md) — spatial nutrient field **PASS**
  ([`JOB3B_STANDING.md`](JOB3B_STANDING.md))
- Job 3c: [`PROTOCOL.md`](PROTOCOL.md) — colony dish, radial nutrient
  limitation **PASS** ([`JOB3C_STANDING.md`](JOB3C_STANDING.md))
- Job 5: [`PROTOCOL.md`](PROTOCOL.md) — Weber well-mixed QS
  **PARTIAL / gates 1–2 FAIL** ([`JOB5_STANDING.md`](JOB5_STANDING.md))
- Job 5b: [`PROTOCOL.md`](PROTOCOL.md) — Dilanji mm-scale spatial AHL
  transport **PASS** ([`JOB5B_STANDING.md`](JOB5B_STANDING.md))
- Job 6: [`PROTOCOL.md`](PROTOCOL.md) — zero-carbon starvation viability
  **PASS** ([`JOB6_STANDING.md`](JOB6_STANDING.md))
- Job 7: [`PROTOCOL.md`](PROTOCOL.md) — basal *E. coli* run-and-tumble
  motility **PASS** ([`JOB7_STANDING.md`](JOB7_STANDING.md))

Gap-fill campaign status and remaining coupling boundaries:
[`GAPFILL_STANDING.md`](GAPFILL_STANDING.md).

Same-day integration into a **new** PocketDish-class garage (do not
edit HybridDish claim-dish Java): copy
[`CHAT_INTEGRATION_BOOTSTRAP.md`](CHAT_INTEGRATION_BOOTSTRAP.md)
into a **new** chat.

```
compile_and_run.cmd
python check_job2.py
python check_job3.py
python check_job3b.py
python check_job3c.py
python check_job5b.py
python check_job6.py
python check_job7.py
```

Optional preview: `compile_and_run.cmd preview`  
Job 2 only: `compile_and_run.cmd job2`  
Job 3 only: `compile_and_run.cmd job3`  
Job 3b only: `compile_and_run.cmd job3b`
Job 5b only: `compile_and_run.cmd job5b`
Job 6 only: `compile_and_run.cmd job6`
Job 7 only: `compile_and_run.cmd job7`
