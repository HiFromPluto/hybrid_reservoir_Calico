# PocketFill-T1b — seed occupancy (job protocol)

Frozen 2026-08-20. Architecture:
[`examples/PocketDish/POCKETFILL_T1B_FROZEN_BUILDER_PROMPT.md`](../PocketDish/POCKETFILL_T1B_FROZEN_BUILDER_PROMPT.md),
[`examples/PocketDish/T1_STANDING.md`](../PocketDish/T1_STANDING.md).
Geometry / QS_* clone: T1 reduced `OPEN_BUS_3` + D1g verbatim.
T1 `VOLUMETRIC_LUXI` is **DEAD**. This job does **not** edit T1 results.

No NARMA, ridge, living Java, membrane LuxI, HybridDish \(k\), `W20`,
`QS_KMLA` hunt, extra seed hunt.

## Question

Same T1 volumetric recipe (Danino \(k\), D1g QS_*, \(N_{\mathrm{pack}}=5000\),
cells write, T0 `OPEN_BUS_3` garage). Does a **0.05 µM seed** occupy
optical LA, or does the leak still wipe AHL before autoinduction?

| Arm | Seed | Expectation |
|---|---|---|
| `VOLUMETRIC_ZERO` | none (\(C(0)=0\)) | T1 replay; **DEAD**. If ALIVE, the port broke |
| `SEED_PULSE` | pocket \(C(0)=0.05\) µM, then **free** | **primary** QS-cure; unknown |
| `SEED_BATH` | pocket extracellular **held** at 0.05 µM | diagnostic / `ENGINEERING`; D1g-class forcing |

`SEED_BATH` ALIVE with `SEED_PULSE` DEAD means: nucleation needs a
maintained bath; leak still wins. That is **not** a PocketFill occupancy
win.

If `SEED_PULSE` is DEAD: stop. Do not retune `QS_KMLA` or the seed.
Do not start NARMA. Do not move the garage to `W20`.

## Frozen numbers (copy T1 / D1g)

- Pocket 100×100×10 µm, full-width open +y, bus 400×80, \(v=3\) µm/s, dx=5 µm
- \(D=159\) µm²/s, `CONV=602.2`
- Volumetric decay \(k=2.76\times10^{-3}/60\) s⁻¹ (Danino / T1, **not** HybridDish 0.0033)
- \(N_{\mathrm{pack}}=5000\), `CELL_VOL=1` µm³ (ENGINEERING)
- QS_* and `µ=ln2/1800` as D1g; 4-ODE ICs `{0,0,0,0}`
- Mean-field write: total membrane flux
  \(N_{\mathrm{pack}} V_{\mathrm{cell}}\times 602.2\times\) `CELL_WALL_DIFF` \((A_{\mathrm{in}}-C)\)
  spread uniformly over pocket voxels
- `SEED_PULSE`: pocket voxels only at 0.05 µM; bus starts at 0; no replenish
- `SEED_BATH`: hold pocket voxels at D1g `BATH_UM=0.05` µM for the whole run
- All arms **14400 s**, sample ≤60 s, last-**7200 s** means
- No AC point source

## Gates (not NRMSE)

- \(H(x)=x^2/(K_{\mathrm{MLA}}^2+x^2)\), \(K_{\mathrm{MLA}}=0.01\)
- Optical (primary): last-7200 s pocket-mean \(H(\mathrm{LA})\)
- Field diagnostic: last-7200 s pocket-mean \(H(C)\)
- ALIVE if optical \(H(\mathrm{LA})\ge 0.05\)
- SATURATED if \(H(\mathrm{LA})>0.95\)
- Also print mean \(C\), mean LA, mean LuxI, whether LA is still rising
- Mass residual class on the field budget; \(C\ge 0\)

`TRANSPORT_MODEL_STATUS`: residual ≤1% of dominant budget term; \(C\ge 0\).
