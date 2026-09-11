# P0 protocol — packed Danino-class device (chemistry OFF)

**Gate:** P0 packed-device mechanics.  
**Chemistry:** OFF. No AHL field, no Object B circuit, no Hill R/L, no AC,
no NARMA, no LuxI.  
**Organism (mechanics only):** Scratch-class *E. coli* rods.  
**Object B:** OFF. **NOT_FIG4B.** Object A remains FAIL / FAIL_NO_IDENTITY
and is not ported.

D0 remains FAIL. D0b remains FAIL_NO_IDENTITY. D1 remains PASS on
circuit parity only; D1 does not authorize packing claims. This protocol
does not reopen Fig. 4b, does not shrink \(\mu=0.32\)–\(0.40\), and does
not retune Object B parameters.

Claim freeze:
[`examples/PocketDish/CLAIM_FREEZE_DANINO_SI.md`](../PocketDish/CLAIM_FREEZE_DANINO_SI.md).

N0 scheduler / RNG:
[`examples/PocketDish/N0_NUMERICAL_KERNEL_STANDING.md`](../PocketDish/N0_NUMERICAL_KERNEL_STANDING.md).

Geometry class:
[`examples/PocketDish/GEOMETRY_FREEZE.md`](../PocketDish/GEOMETRY_FREEZE.md).

**frozen_before_traces:** true  
**Frozen:** 2026-08-27, after N0 PASS and D1 PASS, before any P0
production run.

I0c CSVs are **not** this gate. ChassisPocket I0c is useful rod /
Hertzian / door evidence on a 100×100×1 µm garage and the legacy
`BSim.export()` clock. P0 is a new job.

---

## Scientific question

Does a dense growing monolayer of Scratch rods, in a Danino-class
square pocket of SI bulk height 1.65 µm, remain occupied and turn over
**only** by growth-driven extrusion at the open edge, with chemistry
OFF and motility OFF?

**MAY:** packed-device mechanics; extrusion-only population control;
monolayer retention after the colony reaches the mouth.

**MAY NOT:** Fig. 4b; synchronized GFP; Object B occupancy; AHL /
Hill / NARMA / AC; a physical chemostat; HybridDish living-layer PASS;
permission to start C0 if this gate FAILs.

C0 (Object B well-mixed coupling, NOT_FIG4B) may start **only** if
P0 PASSes.

---

## Device freeze (declared once, before traces)

```
                    BUS (particle sink only; chemistry OFF)
        ╔══════════════════════════════════════════════════╗
        ║         400 × 80 × 1.65 µm mask voxels           ║
        ╚══════════════╤══════════════════════╤════════════╝
                       │   OPEN EDGE  +y      │
                       │                      │
                       │  SQUARE POCKET       │  FLOW_pocket = 0
                       │  100 × 100 × 1.65 µm │  three walls: Hertzian
                       │  founder (50, 80)    │  particles: mirror
                       │  motility OFF        │
                       └──────────────────────┘
```

| Item | Frozen value | Class |
|---|---|---|
| Pocket \(L_x\) | **100.0 µm** | LITERATURE_LAYOUT + ENGINEERING square (PocketDish-A) |
| Pocket \(L_y\) | **100.0 µm** | same |
| Pocket \(L_z\) | **1.65 µm** | TAKEN Danino SI bulk trap. **Not** ChassisPocket 1 µm. **Not** HybridDish 10 µm. |
| Open-edge orientation | entire **+y** face of the pocket | GEOMETRY_FREEZE handedness |
| Pocket flow | **0** | frozen |
| Bus | **present as particle sink only** | GEOMETRY_FREEZE bus class; no advection of a chemical field |
| Bus mask | 400 µm along \(x\), 80 µm wide in \(+y\), height 1.65 µm, centred on the pocket in \(x\) | ENGINEERING layout, chemistry not advanced |
| Particle domain | pocket box only: \(100\times100\times1.65\) µm | particles that enter the bus are removed |
| Particle BC | Hertzian / reflect on \(x=0\), \(x=L_x\), \(y=0\), \(z=0\), \(z=L_z\); **no** Hertzian on \(+y\); remove if centre \(y > L_y + 10^{-3}\) µm | GEOMETRY_FREEZE particle rule |
| Chemical field | **OFF** (not constructed, not stepped) | P0 scope |
| Voxel / mask | shared N0-class mask, declared below, written at \(t=0\) | for C1 later; chemistry still OFF |
| Founder count | **1** | ENGINEERING |
| Founder centre | **(50, 80, 0.825) µm**, long axis along \(x\) | ENGINEERING: same \(xy\) door-reach placement as I0c, mid-plane of the 1.65 µm trap. **Not** a hunt after P0 traces. Centre-founded \(y=50\) does not reach \(y=100\) on this horizon (I0c \(R\approx 30\) µm at \(N=512\)). |
| Seed | **101** | ENGINEERING; injected `BSimRandom` |
| Scheduler | `BSimStepScheduler` | N0 contract. **Not** `BSim.export()` inclusive ticker |
| \(dt\) | **1.0 s** | `ChassisParameters.DT_S` |
| Mechanics period | **1.0 s** (every base step) | integer multiple of \(dt\) |
| \(T\) | **25950 s** = \(N_{\mathrm{steps}}\times dt\) with \(N_{\mathrm{steps}}=25950\) | ENGINEERING: I0c horizon; 10 × \(T_{\mathrm{div}}\) at bath |
| \(\log dt\) | **10.0 s** | ENGINEERING |
| Nutrient | uniform bath \(C_s=0.5\) mM, **no PDE** | Warren / I0c; not a glucose claim |
| Division | `SYMMETRY_BROKEN` | Job 3c / I0c |
| Hertzian | Job 3 `k_cc`, `k_ac`; `WallSpec.openYPlus()` | TAKEN springs; **do not retune**. Additive API; closed-box path unchanged |
| Motility | **OFF** | |
| Death / removal clamp | **OFF** | no `p_removal`, no global random death, no hidden overlap deletion |
| Object B / AHL / Hill / LuxI / AC / NARMA | **OFF** | |

If this height were changed to 1 µm or 10 µm, that would be ENGINEERING
and **must not** be called the Danino SI trap. P0 uses 1.65 µm.

### Voxel mask (chemistry OFF; shared with N0 later)

World origin of the **particle** pocket: southwest-bottom corner.

Mask lattice origin: \((x,y,z)=(-150, 0, 0)\) µm so the 400 µm bus is
centred on the 100 µm pocket in \(x\).

| Item | Frozen value |
|---|---|
| \(dx=dy\) | 2.0 µm |
| \(dz\) | 1.65 µm (\(n_z=1\); monolayer slab) |
| \(n_x,n_y,n_z\) | 200 × 90 × 1 |
| Pocket fluid | voxel centres inside \([0,100]\times[0,100]\times[0,1.65]\) |
| Bus fluid | voxel centres inside \([-150,250]\times[100,180]\times[0,1.65]\) |
| Solid | all other voxels (shoulders beside the pocket) |
| Pocket opening | shared face at \(y=100\) µm (fluid–fluid pocket/bus) |

The mask is built by `bsim.p0.DaninoPocketGeometry`, hashed, and written
before traces. P0 does **not** call `BSimTransportField` transport.

### Clock (N0)

For duration \(T=N\,dt\):

1. observe once at \(t=0\);
2. perform exactly \(N=25950\) updates;
3. boundary / transport / sample / integrate / deposit are no-ops
   (chemistry OFF);
4. mechanics (elongate, Hertzian, spill, divide) every base step;
5. observe the resulting state at \((step+1)\,dt\).

CSV rows are written at \(t=0\) and every \(\log dt=10\) s (integer
multiple of the base step).

---

## Cells

Scratch rods via `EcoliRodCell` / `ValdezHertzian` / `ChassisParameters`:

- Valdez-style elongation at uniform bath;
- symmetry-broken division;
- Hertzian cell–cell and cell–wall contacts;
- motility OFF;
- all stochastic draws from the injected `BSimRandom` only
  (`sim.setRandomSeed(101)`; `divisionRng = sim.getRandom().asJavaRandom()`);
- no `Math.random()`.

Population control is growth + packing + extrusion. If \(N\) is capped
by random removal anywhere in the pocket, P0 = FAIL.

---

## Pre-registered gates

Output: `examples/BSimReservoirPlanDaninoPocketP0/results/p0_seed101/`.

CSV columns (semicolon):

`t_s;seed;N;R_um;R_disc_um;R_Rdisc;offplane;z_exc_max_um;delta_cc_max_um;d_centers_min_um;spill_tick;spill_cum;N_ever;door_contact;y_max_um;wall_leak;n_drop_unexplained`

Definitions (I0c lineage, height-adjusted monolayer):

- `N` = cells still in the pocket.
- `spill_tick` = cells removed this log interval because
  `centre.y > L_y + 1e-3` µm (left through \(+y\) into the bus sink).
- `spill_cum` = running total of those removals. **May be > 0.**
- `N_ever` = `N + spill_cum` (every cell born, still in or spilled).
  Binary-fission checks apply to `N_ever`.
- `door_contact` = 1 if any remaining pole has \(y \ge L_y - r\)
  (\(r=0.5\) µm ⇒ 99.5 µm), else 0.
- `y_max_um` = max remaining pole \(y\).
- `wall_leak` = centres that left through a face that is **not** \(+y\).
  Must stay 0.
- `n_drop_unexplained` = decrease in `N` during a log interval that is
  **not** accounted for by `spill_tick` (and wall-leak removals). Must
  stay 0. This is the death-clamp / hidden-deletion detector.
- `R`, `R_disc` = Job 3c / I0c: max pole \(xy\) distance from COM;
  `R_disc = sqrt(N * footprint / 0.85 / π)` with
  `footprint = w0 * (0.5*(ℓ0+ℓ_div) + 2r)` = 3.0 µm².
- `z_exc_max_um` = \(\max_i |z_i - L_z/2|\).
- `offplane` = fraction of centres with \(|z - L_z/2| > 0.40\) µm.
  At SI height 1.65 µm a monolayer sitting in the slab has geometric
  excursion \(\le L_z/2 - r = 0.325\) µm. The I0c threshold
  \(0.5 r = 0.25\) µm is **not** reused: it would flag every cell at
  this height. 0.40 µm is the P0 monolayer cap (cannot stack a second
  1 µm rod in 1.65 µm).

| Gate | Pass if | Fail if |
|---|---|---|
| **P0.1 growth** | `N_ever ≥ 256` at \(t=25950\) s; `N_ever` never decreases; `N_ever` does not more than double between log rows | growth stalls; death clamp |
| **P0.2 packing** | `delta_cc_max ≤ 0.25` µm and `d_centres_min ≥ 0.75` µm at every logged row after `N≥8` | Hertzian overlap blow-up |
| **P0.3 monolayer** | first row with `N≥256`: `offplane ≤ 0.05`, `z_exc_max ≤ 0.40` µm, `R/R_disc ∈ [0.7, 1.3]`, `d_centres_min ≤ 1.2` µm | filament (`R/R_disc>2`), stacking, or motile gas |
| **P0.4 retention** | some row has `door_contact=1`; `N_end ≥ 256`; `wall_leak=0` on every row | never reached the mouth; `N_end=0` EMPTY (T4/T5 class); drained below 256; closed-face leak |
| **P0.5 extrusion-only** | `spill_cum=0` on every row before the first `door_contact=1`; `spill_cum>0` by \(t_{\mathrm{end}}\); `n_drop_unexplained=0` on every row | spill before door contact; no extrusion after contact; interior deletion |
| **P0.6 repeatability** | two smoke runs, same seed, **byte-identical** `colony_timeseries.csv` | RNG split or hidden `Math.random()` |
| **P0.7 honesty** | stdout and standing print `chemistry=OFF`, `motility=OFF`, `Object B=OFF`, `NOT_FIG4B`; height logged as 1.65 µm; scheduler is `BSimStepScheduler`; I0c CSVs not cited as this gate | Fig. 4b claim; silent 10 µm or 1 µm swap; `BSim.export()` used as the P0 clock |

`spill_cum>0` is **not** a fail. T4/T5 failed because `N_end=0`.

Do **not** retune `k_cc`, growth rate, \(L_z\), founder, or seed after
seeing holes.

### Smoke (before the production run)

Stop expensive mechanics at `N=32`. Scheduler duration for smoke is
**15570 s** = \(6\times T_{\mathrm{div}}\) (predeclared; long enough to
reach \(N=32\)). Confirm `offplane=0`, `R/R_disc < 2`,
`door_contact=0`, `spill_cum=0`. If filament, **stop**. Do not start
the 25950 s run.

Run smoke **twice** (same seed) for P0.6 before production.

---

## Integrity fences

Do **not**:

- port Object A / Fig. 4b into Java, logs, or the standing;
- enable Object B, AHL, Hill, LuxI, AC, or NARMA;
- rewrite D0, D0b, or D1 standings;
- use `BSim.export()` inclusive `0..N` as the P0 clock;
- enable motility to refill holes;
- add a global death clamp or hidden overlap deletion;
- retune `k_cc` or growth after seeing holes;
- silently swap HybridDish 10 µm or ChassisPocket 1 µm for 1.65 µm;
- treat I0c CSVs as this gate;
- start C0 unless P0 PASSes;
- commit unless asked.

Chassis API used: `EcoliRodCell`, `ValdezHertzian`, `ChassisParameters`.
Additive `WallSpec.openYPlus()` is allowed if Job 2/3/3b/3c CSV bytes
do not move.

---

## Commands

Smoke twice, then production. Quote the Ant `-D` property on PowerShell.

```
ant p0-packed "-Dp0.args=smokeA"
ant p0-packed "-Dp0.args=smokeB"
python examples/BSimReservoirPlanDaninoPocketP0/check_p0.py --smoke
ant p0-packed
python examples/BSimReservoirPlanDaninoPocketP0/check_p0.py
```

Standing: `examples/PocketDish/P0_PACKED_DEVICE_STANDING.md`.
Must repeat chemistry=OFF, motility=OFF, Object B=OFF, NOT_FIG4B,
and whether C0 may start.
