# PocketFill-T1 — occupancy (job protocol)

Frozen 2026-08-20. Architecture:
[`examples/PocketDish/NEXT_DESIGN.md`](../PocketDish/NEXT_DESIGN.md),
[`examples/PocketDish/POCKETFILL_T1_FROZEN_BUILDER_PROMPT.md`](../PocketDish/POCKETFILL_T1_FROZEN_BUILDER_PROMPT.md).
Geometry clone: T0 `OPEN_BUS_3`.
QS_*: [`examples/BSimDaninoD1g/BSimDaninoD1g.java`](../BSimDaninoD1g/BSimDaninoD1g.java)
verbatim. D1g remains open-loop PASS (box). This job **writes** the
field. D1/D1g/T0/Narma10b Overall lines are not rewritten.

No NARMA, ridge, living Java, Grober, two-way, vesicles, glucose.

## Question

Same Danino-class bay as T0. Does the **Danino recipe** occupy?

| Arm | Expectation |
|---|---|
| `EMPTY_POINT` | DEAD (T0 replay, HybridDish Hill \(K=1.6\)) |
| `PACKED_PLUG` | unknown; \(D_{\mathrm{eff}}=D/4\) is `HYPOTHETICAL_DESIGN_ENVELOPE` |
| `VOLUMETRIC_LUXI` | candidate; optical `H(LA)` vs `QS_KMLA=0.01` |

If volumetric DEAD: stop. Do not retune `QS_KMLA`, \(N_{\mathrm{pack}}\),
or T0 \(J_{\max}\).

## Frozen numbers

- Pocket 100×100×10 µm, bus 400×80, \(v_{\mathrm{bus}}=3\) µm/s, dx=5 µm
- \(D=159\) µm²/s
- `EMPTY_POINT` / `PACKED_PLUG`: \(k=0.0033\) s⁻¹, \(J_{\max}=6.36\times10^5\), \(u=0.5\), 9000 s, Hill \(K=1.6\), \(n=2\), \(\tau_R=15\) s
- `PACKED_PLUG`: pocket \(D/4\), bus \(D\), harmonic faces
- `VOLUMETRIC_LUXI`: \(k=2.76\times10^{-3}/60\) s⁻¹, no point source, 14400 s, last-7200 s means
- \(N_{\mathrm{pack}}=5000\), `CELL_VOL=1` µm³ (ENGINEERING)
- QS_* and `µ=ln2/1800` as D1g; ICs zero; `CONV=602.2`
- Mean-field 4-ODE; membrane flux writes the field uniformly in the pocket

Occupancy ALIVE: HybridDish arms `mean_R≥0.05`; volumetric last-7200 s mean `H(LA)≥0.05`.
