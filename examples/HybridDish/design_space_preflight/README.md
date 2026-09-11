# HybridDish design-space preflight

Analysis-only package for Track E0.1–E0.3. It does not run Java or BSim,
does not change historical gates, and does not authorize Waveform2,
Lorenz, IPC, AC implementation, biomarker BSim, independent-input NARMA,
or Stage99.

## Rerun from one command

From the repository root:

```
python examples/HybridDish/design_space_preflight/run_transport_screen.py
```

`PROTOCOL.md` freezes the rules before that script writes results.
Unit tests A-E run first inside the same command.

Failed 2026-08-17 mass-gate rows are retained in
`transport_validation_failed_2026-08-17.csv`.

## Package map

| File | Role |
|---|---|
| `PROTOCOL.md` | Frozen analysis rules |
| `POST_GATE_DISPOSITION.md` | Module 1; `C1_DISPOSITION: DEFER` |
| `PARAMETER_PROVENANCE.csv` / `.md` | Implemented HybridDish inventory |
| `run_transport_screen.py` | Reduced transport/kinetic surrogate |
| `TRANSPORT_MODEL_VALIDATION.md` | Frozen-gate comparison |
| `transport_validation.csv` | Per-seed validation numbers |
| `transport_screen_conditions.csv` | Predeclared screen table |
| `transport_screen_results.csv` | Raw physical metrics |
| `reset_screen_results.csv` | Residual-state curves |
| `test_mass_budget.py` | Unit tests A-E; also run from the screen command |
| `transport_validation_failed_2026-08-17.csv` | Archived mixed-horizon mass-gate failure |
| `E0_PREFLIGHT_REPORT.md` | Integrated answers |
| `E0_3_SHORTLIST.csv` | At most 12 E0.3 conditions after validation |
| `E0_3_FROZEN_BUILDER_PROMPT.md` | Next authorized builder task |

Current status: `TRANSPORT_MODEL_STATUS: VALIDATED_FOR_SCREENING`.
Screen CSVs are `DECISION_GRADE_PHYSICAL_SCREEN`. Living E0.3 BSim is
requested in `E0_3_FROZEN_BUILDER_PROMPT.md` but is not executed here.
