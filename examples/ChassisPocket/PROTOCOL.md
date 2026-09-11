# ChassisPocket I0 protocol — packed monolayer fill

**FROZEN BEFORE IMPLEMENTATION (2026-08-24).** Every threshold below is
pre-registered. Do not retune a threshold, `k_cc`, `λ_S`, `K_S`, or growth
after seeing `N`. If a gate fails, the finding is the failure.

I0 is **locked PASS** (`CHASSISPOCKET_I0_STANDING.md`). I0b is **locked
FAIL** (door not reached; `CHASSISPOCKET_I0b_STANDING.md`). This file now
also contains the **I0c freeze**. I1 / I2 each need their own freeze
before any code. Do not stack I1 or I2 into the I0c Java file.

## Naming

This is **Job I0** of `ChassisPocket`, a **new named dish**. It is not
HybridDish, not PocketDish Java, and not a chassis gap-fill job.

PocketDish T4/T5 emptied because `BSimBacterium` spheres do not jam a pad.
Jobs 2 / 3 / 3c already have the missing mechanism (Valdez capsule, Hertzian
contacts, `b_z = 1 µm` monolayer, `DivisionMode.SYMMETRY_BROKEN`). I0 asks
whether that locked colony **stays in a PocketDish-A-class garage**.

## Scientific claim I0 may make

**MAY:** a Hertzian *E. coli* rod colony, grown from one founder in a
100×100×1 µm closed garage with Job 3c division and wall contacts, remains a
monolayer disc and does not empty.

**MAY NOT:** occupancy of Hill `R`, NARMA skill, door-as-chemostat, spatial
GFP domains, death, swimming, Weber hysteresis.

I0 is **not** an exact Danino blueprint, **not** a physical chemostat, and
**not** HybridDish living-layer PASS.

## Physical-scenario freeze

I0 is a **surface colony in a square garage**. Not a millimetre HybridDish,
not a swimming chamber, not a mm AHL lane.

```
PocketDish-A class, CHASSIS_MONOLAYER amendment:

                    (neck is a WALL for particles in I0)
        ╔══════════════════════════════════════╗
        ║  three no-flux / Hertzian walls      ║
        ║                                      ║
        ║     SQUARE POCKET  100 × 100 × 1 µm  ║
        ║     founder at centre                ║
        ║     motility OFF                     ║
        ║                                      ║
        └──────────────────────────────────────┘
```

| Item | Value | Class | Provenance |
|---|---|---|---|
| footprint | `100 × 100 µm` | TAKEN layout | PocketDish-A square |
| height `b_z` | `1.0 µm` | ENGINEERING | Job 3c monolayer (`2r`); **not** PocketDish 10 µm |
| neck | **closed to particles** | ENGINEERING | I0 asks “does Hertzian jam?”, not “does the door leak cells?” |
| AHL field | **OFF** | — | AC perfume is T3 class; not this question |
| nutrient | uniform bath `C_s = 0.5 mM` | TAKEN | Warren; **no PDE in I0** |
| occupancy | unused (no PDE) | — | `CONSERVATIVE_MASS` is I1 |
| `dt` | `1.0 s` | ENGINEERING | `ChassisParameters.DT_S` |
| `t_end` | `23400 s` | ENGINEERING | Job 3c horizon; enough for N~256–512 |
| `log_dt` | `10.0 s` | ENGINEERING | `ChassisParameters.LOG_DT_S` |
| division | `SYMMETRY_BROKEN`, seed 101 | ENGINEERING | Job 3c |
| Hertzian | Job 3 `k_cc`, `k_ac`, frictions | TAKEN | do not retune after overlap |
| motility | OFF | — | Job 7 stays standalone |
| death | OFF | — | Job 6 stays standalone |
| QS / Hill / LuxI | OFF | — | |

**Why close the neck in I0:** T4 emptied through `W20` because spheres do not
jam **and** the door is 20 cell-widths. Mixing those two failures in one CSV
repeats T5. First prove the garage holds a monolayer with walls. Door-as-weir
is **I0b**, only if I0 PASSes.

**Why `b_z = 1 µm` not 10 µm:** Job 3c at 4 µm stacked (`offplane=0.71`).
PocketDish 10 µm is a swimming tank. Image-2 / Danino is a 1–2 µm monolayer.
I0 uses the locked monolayer, not the 10 µm garage that T5 already showed
does not pack.

**Why no nutrient PDE:** PocketDish protocol forbids mixing glucose into
PocketHill; I0 is a packing/fill job. Uniform Warren bath. I1 later.

Cells will **not** fill the 100 µm square at N=256 (`R≈17 µm`). That is fine.
I0 tests **non-emptying + monolayer**, not “pad is full.” Filling the square
needs N ~ 10^4 and is not this run.

Compare against T5 EMPTY (`N_end=0`, `N_max=50`) and Job 3c PASS (`N=256`,
`offplane=0`, `R/R_disc≈0.99`).

## Integrity fences

Do **not**:

- edit `examples/HybridDish/bsim/BSimHybridDish.java`
- edit any `GATE_EVIDENCE.md` or Narma10b / C1 / E5.2 Overall lines
- edit locked chassis classes in a way that changes Job 2/3/3b/3c CSV bytes
- retune `k_cc`, `λ_S`, `K_S`, `W`, `J_max`, `K`, `n` after seeing `N`
- turn motility on in the packed garage
- turn Biselli death on without a zero-carbon transfer protocol
- add a spatial AHL PDE in this 100 µm pocket and claim a spatial QS result
- put LuxI, NARMA, and Hertzian in one Java file
- call I0 an “exact Danino blueprint” or a “physical chemostat”
- call I0 HybridDish living-layer PASS
- mix glucose PDE into PocketHill Java

Chassis API used by I0: `EcoliRodCell`, `ValdezHertzian`, `ChassisParameters`.
Do not copy-paste those files. A visibility-only `public` export of
`ValdezHertzian.relaxContacts` / `stats` is allowed if it does not change
Job 2/3/3b/3c CSV bytes.

## Pre-registered gates

Output: `examples/ChassisPocket/results/i0_seed101/`.

CSV columns (semicolon):

`t_s;seed;N;R_um;R_disc_um;R_Rdisc;offplane;delta_cc_max_um;d_centers_min_um;spill`

`spill` is the count of cells with any pole strictly outside the closed box
`[0,100]×[0,100]×[0,1]` µm (tolerance `1e-3` µm). I0 must keep `spill = 0`.

Do not add a “mean_R > 0.05” gate. There is no AHL.

| Gate | Pass if | Fail if |
|---|---|---|
| **I0.1 growth** | `N_final ≥ 256` at `t=23400 s` (binary fission; N never decreases; N does not more than double between log rows) | N stalls or collapses |
| **I0.2 packing** | `delta_cc_max ≤ 0.25 µm` and `d_centres_min ≥ 0.75 µm` at every logged row after N≥8 | Hertzian overlap blow-up |
| **I0.3 morphology** | at first row with N≥256: `offplane ≤ 0.05`, `R/R_disc ∈ [0.7, 1.3]`, `d_centres_min ≤ 1.2 µm` | filament (`R/R_disc>2`) or stacking |
| **I0.4 fill (the T4/T5 question)** | last-row `N = N_final` (no emptying); `spill = 0` on every row | last cell time ≪ t_end; N→0; poles leave the box |
| **I0.5 chassis regression** | `check_job2/3/3b/3c` exit 0; Job 3b CSVs SHA256-identical to the freeze below | any lock moves |
| **I0.6 honesty** | standing memo states AHL/Hill/NARMA were OFF | claim occupancy or task score |

If I0.3 fails as a filament: that is the Job 3c `COLLINEAR` bug; check
`SYMMETRY_BROKEN` + seed 101. Do not add noise beyond the frozen
`DIVISION_*` magnitudes.

If I0.4 fails with cells escaping through a “closed” wall: the wall Hertzian
is wrong; fix the boundary, do not lower growth.

### Morphology definitions (same as Job 3c)

- `R` = max pole *xy* distance from the colony centre of mass.
- `R_disc` = `sqrt(N * footprint / 0.85 / π)` with
  `footprint = w0 * (0.5*(ℓ0+ℓ_div) + 2r)` (= 3.0 µm²).
- `offplane` = fraction of cell centres more than `0.5 r` off the mid-plane
  `z = b_z/2`. At `b_z = 2r` this is a regression guard.

### Smoke (before the science run)

Stop at `N = 32`. Confirm `offplane = 0` and `R/R_disc` falling toward 1, not
a filament (`R/R_disc > 2`). If filament, **stop**. Do not start the 23400 s
run.

Job 3c smoke at N=32 in the 60 µm monolayer: `R/R_disc = 1.44`, `offplane = 0`.
A 100 µm closed garage should look the same at N=32 (`R ≈ 8 µm`, walls far).

## I0.5 freeze hashes (byte-identical chassis CSVs)

Recorded 2026-08-24 **before** any I0 code, after `check_job2/3/3b/3c` all
exited 0. Paths relative to `examples/BacteriumFromScratch/`.

```
e5ad8c7809d6b1712b47d4f231de59d71698aab888980e644c7604392df746ae  results/job2_seed101/job2_timeseries.csv
b3a41af567834cc3076056027f9608f46bd541f8d191d88c8119c690b90eb5c9  results/job3_seed101/twobody_timeseries.csv
bceef8bc7a417a069552a82ebeb54efbd3f2ec8df66e9e67760292e512941774  results/job3_seed101/cluster_timeseries.csv
cb21d5ecce63b59813f33ef54af8388ce70fdaf26608356afb37d72a06c6cf9a  results/job3_seed101/isolated_packing_timeseries.csv
351f2f2ffa5d8dc23b71dd0697eab7dc3056bdd3fcde017925537c6ff311d829  results/job3b_seed101/depletion_timeseries.csv
e6fe563fa6c1b02db5c72c0a97d960c825e2a6658f96348050d48b2ece577a1c  results/job3b_seed101/packing_cluster_timeseries.csv
ebf19ffa47e980f8b04b827ce9bcfd02bb190c500958c2ea717ec55369e74bbc  results/job3b_seed101/uniform_timeseries.csv
f1d83f2e87ee36d3217568a78eb3ec9438c94886b712f0fa76e8c66640dc20cc  results/job3c_seed101/colony_timeseries.csv
```

Job 3b CSVs are the load-bearing byte-identity lock. The others are
regression guards for the same rule.

## OFF in I0

Nutrient PDE, AHL PDE, Hill `R,L`, LuxI, Danino 4-ODE, NARMA / ridge / 408-D,
motility, Biselli death, acid clamp-death, second species, chemotaxis,
metabolism ODE, Brownian, `F_s`, stochastic `P(ℓ)`.

## Later jobs after I0 (I0b is the next freeze, below)

I1 NutrientPocket, I2 HillOnRods, I3 ChassisHybridmm, I4 WeberHost,
I5 StarvationArm, I6 SwimChamber. Each is a new PROTOCOL section. Do not
stack. Do not start them unless I0b is scored.

---

# ChassisPocket I0b protocol — DoorWeir at T4/T5 width

**FROZEN BEFORE IMPLEMENTATION (2026-08-24), after I0 PASS.** Every
threshold below is pre-registered. Do not retune a threshold, `k_cc`,
`λ_S`, `K_S`, `W`, or founder placement after seeing `N` or spill. If a
gate fails, the finding is the failure.

## Naming

This is **Job I0b** of `ChassisPocket`. It is not HybridDish, not
PocketDish Java, and not a chassis gap-fill job. It is **not** I0: the
neck is no longer a wall.

## Scientific claim I0b may make

**MAY:** a Hertzian *E. coli* rod colony, grown from one founder in the
I0 garage with the T4/T5 particle door open at frozen `W = 20 µm`,
remains a monolayer in the garage and does not empty through that door.

**MAY NOT:** occupancy of Hill `R`, NARMA skill, door-as-chemostat,
spatial GFP, AHL leak (AHL is OFF), death, swimming, Weber hysteresis,
“1–3 µm mouth keeps cells while AHL leaks.”

I0b is **not** an exact Danino blueprint, **not** a physical chemostat,
**not** HybridDish living-layer PASS, and **not** PocketHill occupancy.

## Why `W = 20 µm`, not 1–3 µm

T4/T5 emptied through `W20_L20`. I0 closed that door so packing and
door-leak were not mixed. I0b reopens **the same particle width T4 used**.

A 1–3 µm mouth would confound **cell body vs door width**. If I0b at
`W = 20` FAILS, a **separate later job** may freeze `W = 3 µm`. That is
not a retune of a failed I0b gate. Do not add `W = 10` or `W = 50` in
this run.

`L_n = 0` for particles: a gap in the `+y` wall, not the T4 bus and not
the 20 µm corridor. Cells whose **centre** leaves through `+y` are
removed (spillover). The bus is not simulated. Mouth-only is the
harsher particle test (no neck walls to jam in). AHL is OFF, so corridor
length is not a chemical leak question.

## Why the founder is at `(50, 80)`, not I0’s centre

I0 locked the centre-founded closed garage: at `t_end = 23400 s`,
`N = 512`, `R ≈ 21.75 µm`. The `+y` wall is 50 µm from the centre, so
that colony **never meets the door** on the I0 horizon.

I0b’s question is the door. Keep the I0 clocks (`t_end = 23400 s`,
`dt = 1 s`) and place the founder so the disc can geometrically contact
`y = 100` inside that horizon:

- founder centre `(BX/2, 80, BZ/2)` µm, long axis along `x` (same as I0)
- distance to `+y` = 20 µm
- I0’s `R(N=256) ≈ 16.8 µm` and `R(N=512) ≈ 21.8 µm` ⇒ contact during
  the 8th–9th doubling

This offset is ENGINEERING, frozen from the **locked I0 standing**, not
from I0b spill. Do not move the founder after seeing `N`. Do not extend
`t_end` after seeing spill. Do not switch to a 1–3 µm mouth because the
colony aims at the hole.

## Physical-scenario freeze

```
PocketDish-A class, CHASSIS_MONOLAYER amendment, T4 W20 mouth:

                    (neck is OPEN, W = 20 µm, centred on +y)
        ╔═════════════════╤════╤═════════════════╗
        ║                 │    │                 ║
        ║     SQUARE POCKET  100 × 100 × 1 µm    ║
        ║     founder at (50, 80), motility OFF  ║
        ║                                        ║
        └────────────────────────────────────────┘
```

| Item | Value | Class | Provenance |
|---|---|---|---|
| footprint | `100 × 100 µm` | TAKEN layout | PocketDish-A / I0 |
| height `b_z` | `1.0 µm` | ENGINEERING | Job 3c / I0 monolayer |
| neck width `W` | **20 µm** | TAKEN layout | T4/T5 `W20`; not a hunt |
| neck centre | `x = 50 µm` on `+y` | ENGINEERING | PocketDish-A handedness |
| neck length `L_n` | **0** (gap in the wall) | ENGINEERING | particle weir; no bus |
| particle BC | Hertzian walls on five faces + `+y` except the gap; **remove** if centre `y > 100` | ENGINEERING | T4 absorb-on-open-edge, no motility |
| AHL field | **OFF** | — | chemical leak is not this question |
| nutrient | uniform bath `C_s = 0.5 mM` | TAKEN | Warren; **no PDE** |
| `dt` | `1.0 s` | ENGINEERING | `ChassisParameters.DT_S` |
| `t_end` | `23400 s` | ENGINEERING | I0 / Job 3c horizon |
| `log_dt` | `10.0 s` | ENGINEERING | `ChassisParameters.LOG_DT_S` |
| founder | centre `(50, 80, 0.5)` µm, axis along `x` | ENGINEERING | door contact on I0 clocks; see above |
| division | `SYMMETRY_BROKEN`, seed 101 | ENGINEERING | Job 3c / I0 |
| Hertzian | Job 3 `k_cc`, `k_ac`, frictions | TAKEN | do not retune |
| motility | OFF | — | Job 7 stays standalone |
| death | OFF | — | Job 6 stays standalone |
| QS / Hill / LuxI | OFF | — | |

Wall Hertzian on `+y` applies only for poles with
`x ∉ (50 − W/2, 50 + W/2)` = `x ∉ (40, 60)` µm. That API is an
**additive** `ValdezHertzian.WallSpec`. The no-argument
`relaxContacts(cells, bound)` path stays the closed box. Job 2/3/3b/3c
CSV bytes must not move.

Compare against T5 EMPTY (`N_end = 0`, `N_max = 50`) and I0 PASS
(`N_end = 512`, `spill = 0`, neck closed).

## Integrity fences

Same list as I0, plus:

- do not narrow `W` after seeing spill or `N_end`
- do not move the founder after seeing spill
- do not turn the bus or `L_n = 20` on in this file
- do not claim AHL occupancy; AHL is OFF
- do not edit I0 Java to open the neck; I0b is a new class

## Pre-registered gates

Output: `examples/ChassisPocket/results/i0b_seed101/`.

CSV columns (semicolon):

`t_s;seed;N;R_um;R_disc_um;R_Rdisc;offplane;delta_cc_max_um;d_centers_min_um;spill_tick;spill_cum;N_ever;door_contact;y_max_um;wall_leak`

Definitions:

- `N` = cells still in the garage (the live list).
- `spill_tick` = cells **removed this log interval** because
  `centre.y > BY + 1e-3` µm (left through `+y`).
- `spill_cum` = running total of those removals. **May be > 0.** Danino
  push-out through the mouth is allowed. T5-class failure is emptying,
  not a non-zero spill count.
- `N_ever` = `N + spill_cum` (every cell that was born, still in or
  spilled). Must not decrease. Binary-fission check applies to `N_ever`,
  not to `N` (garage `N` may drop when cells leave).
- `door_contact` = `1` if any remaining cell has a pole with
  `y ≥ BY − r` (`r = 0.5 µm`), else `0`. I0b is not a door test unless
  this is 1 by `t_end`.
- `y_max_um` = max pole `y` among remaining cells.
- `wall_leak` = count of remaining or just-removed cells whose centre
  left through a face that is **not** `+y` (`x<0`, `x>BX`, `y<0`, `z`
  out). Must stay 0. If it fires, the closed-face Hertzian is wrong;
  do not lower growth.

Morphology `R`, `R_disc`, `offplane` are the I0 / Job 3c definitions,
computed on remaining cells.

| Gate | Pass if | Fail if |
|---|---|---|
| **I0b.1 growth** | `N_ever ≥ 256` at `t = 23400 s`; `N_ever` never decreases; `N_ever` does not more than double between log rows | growth stalls; `N_ever` collapses |
| **I0b.2 packing** | `delta_cc_max ≤ 0.25 µm` and `d_centres_min ≥ 0.75 µm` at every logged row after `N ≥ 8` | Hertzian overlap blow-up |
| **I0b.3 morphology** | at first row with `N ≥ 256`: `offplane ≤ 0.05`, `R/R_disc ∈ [0.7, 1.3]`, `d_centres_min ≤ 1.2 µm` | filament (`R/R_disc>2`) or stacking |
| **I0b.4 fill (the T4/T5 door)** | `N_end ≥ 256`; `wall_leak = 0` on every row; some row has `door_contact = 1` | `N_end = 0` (EMPTY); `N_end < 256` after a colony existed (drained through `W20`); never touched the door; poles left through a closed face |
| **I0b.5 chassis regression** | same as I0.5: checkers exit 0; Job 3b CSVs SHA256-identical to the I0.5 freeze | any lock moves |
| **I0b.6 honesty** | standing memo states AHL/Hill/NARMA were OFF; `W` was not narrowed; founder was not moved after looking | claim occupancy or task score; silent `W` hunt |

`spill_cum > 0` is **not** a fail. T5 failed because `N_end = 0`.

If I0b.4 fails EMPTY at `W = 20`: write that standing. Do not open a
`W = 3` arm in this Java file.

If I0b.4 fails because `door_contact` never fired: the founder offset
was too timid; that is a protocol bug, not permission to retune `k_cc`.

If I0b.3 fails as a filament: same as I0 — check `SYMMETRY_BROKEN` +
seed 101.

### Smoke (before the science run)

Stop at `N = 32`. Confirm `offplane = 0` and `R/R_disc` falling toward
1, not a filament. At N=32, `R ≈ 8 µm`; founder at `y = 80` so the
front is near `y ≈ 88`, still short of the door. Smoke must **not**
require spill. If filament, **stop**.

## OFF in I0b

Same as I0: nutrient PDE, AHL PDE, Hill `R,L`, LuxI, Danino 4-ODE,
NARMA / ridge / 408-D, motility, Biselli death, acid clamp-death, second
species, chemotaxis, metabolism ODE, Brownian, `F_s`, stochastic `P(ℓ)`,
the T4 bus, `L_n = 20` corridor, HybridDish `p_removal`.

## Later jobs (not this freeze)

I0c is the next freeze, below. I1 NutrientPocket, I2 HillOnRods, I3
ChassisHybridmm, I4 WeberHost, I5 StarvationArm, I6 SwimChamber. A
failed I0b/I0c at `W = 20` after the door is actually tested may later
license a named `W = 3 µm` job; that is a new PROTOCOL section.

---

# ChassisPocket I0c protocol — DoorWeir, one extra doubling

**FROZEN BEFORE IMPLEMENTATION (2026-08-24), after I0b FAIL (door not
reached).** Every threshold below is pre-registered. Do not retune a
threshold, `k_cc`, `λ_S`, `K_S`, `W`, or founder placement after seeing
`N` or spill. If a gate fails, the finding is the failure.

## Naming

This is **Job I0c** of `ChassisPocket`. It is not HybridDish, not
PocketDish Java, and not a chassis gap-fill job. It is **not** I0b
edited in place: I0b Java and CSVs stay locked. The only allowed change
versus I0b is `t_end`.

## Scientific claim I0c may make

**MAY:** a Hertzian *E. coli* rod colony, grown from the I0b founder in
the I0b garage with the T4/T5 particle door open at frozen `W = 20 µm`,
run to `t_end = 25950 s` (one extra doubling), contacts that door and
does not empty through it.

**MAY NOT:** occupancy of Hill `R`, NARMA skill, door-as-chemostat,
spatial GFP, AHL leak (AHL is OFF), death, swimming, Weber hysteresis,
“1–3 µm mouth keeps cells while AHL leaks.”

I0c is **not** an exact Danino blueprint, **not** a physical chemostat,
**not** HybridDish living-layer PASS, and **not** PocketHill occupancy.

## Why `t_end = 25950`, not founder `y = 70`

I0b already placed the founder at `(50, 80)`. At `t_end = 23400 s`,
`y_max = 99.16 µm` against the door-contact threshold `99.5 µm` — short
by **0.34 µm**. Moving the founder after seeing `y_max` would look like
a placement hunt. One more doubling is the smallest change that can
clear that miss.

I0b `y_max` went 92.15 (N=256) → 99.16 (N=512). Another `√2` in radius
puts `+y` through 99.5 with margin. Hertzian cost grows; I0b was ~59 s
to N=512; expect a few minutes to N~1024, not a 35 min PDE job.

The founder-`y = 70` arm is **retired** unless I0c still misses the
door, in which case it is a **new** named job, not a quiet edit.

**Why still `W = 20`:** I0b never tested the door. Narrowing `W` now
would confound “did not reach” with “mouth too wide.”

Do not edit `ChassisPocketI0.java` or `ChassisPocketI0b.java` to change
clocks. New class `ChassisPocketI0c`.

## Physical-scenario freeze

Pick is **already made**. Do not take the founder-`y = 70` arm.

```
Same as I0b, one extra doubling:
                    (neck OPEN, W = 20 µm, centred on +y)
        ╔═════════════════╤════╤═════════════════╗
        ║                 │    │                 ║
        ║     SQUARE POCKET  100 × 100 × 1 µm    ║
        ║     founder at (50, 80), motility OFF  ║
        ║     t_end = 25950 s  (10 × T_div)      ║
        └────────────────────────────────────────┘
```

| Item | Value | Class |
|---|---|---|
| footprint | `100 × 100 × 1 µm` | I0 / I0b |
| `W` | **20 µm** on `+y` at `x = 50` | T4/T5; not a hunt |
| `L_n` | **0** (gap in the wall; no bus) | I0b |
| founder | `(50, 80, 0.5)` µm, axis along `x` | I0b; **not moved** |
| `dt` | `1.0 s` | `ChassisParameters.DT_S` |
| `t_end` | **`25950 s`** | ENGINEERING: `10 × 2595 s`; I0b was 9 doublings and 0.34 µm short |
| `log_dt` | `10.0 s` | I0b |
| nutrient | uniform `C_s = 0.5 mM`, no PDE | I0b |
| division | `SYMMETRY_BROKEN`, seed 101 | I0b |
| Hertzian | Job 3 `k_cc`, `k_ac`; `WallSpec.openNeckYPlus(20, 50)` | do not retune |
| particle BC | remove if `centre.y > 100 + 1e-3` µm | I0b |
| motility, death, AHL, Hill, NARMA | **OFF** | |

## Integrity fences

Do **not**:

- edit HybridDish or PocketDish claim dishes, or any `GATE_EVIDENCE.md`
- retune `k_cc`, `λ_S`, `K_S`, `W`, founder `(50, 80)`, or growth after seeing `N` or spill
- edit `ChassisPocketI0.java` or `ChassisPocketI0b.java` to change clocks
- turn motility, death, AHL, Hill, or nutrient PDE on
- add `L_n = 20` or the T4 bus
- call I0c an exact Danino blueprint, a chemostat, or HybridDish living-layer PASS
- claim PocketHill occupancy (`mean_R`) — there is no AHL
- start I1 / I2 from this standing

`spill_cum > 0` is **allowed** (Danino push-out). T5 failed because
`N_end = 0`.

## Pre-registered gates

Output: `examples/ChassisPocket/results/i0c_seed101/`.

CSV columns (semicolon), same as I0b:

`t_s;seed;N;R_um;R_disc_um;R_Rdisc;offplane;delta_cc_max_um;d_centers_min_um;spill_tick;spill_cum;N_ever;door_contact;y_max_um;wall_leak`

Definitions: copy I0b (this file’s I0b section). `door_contact = 1` if
any remaining pole has `y ≥ 99.5 µm`. Morphology definitions: Job 3c / I0.

| Gate | Pass if | Fail if |
|---|---|---|
| **I0c.1 growth** | `N_ever ≥ 256` at `t = 25950 s`; `N_ever` never decreases; `N_ever` does not more than double between log rows | growth stalls |
| **I0c.2 packing** | `delta_cc_max ≤ 0.25 µm` and `d_centres_min ≥ 0.75 µm` after `N ≥ 8` | overlap blow-up |
| **I0c.3 morphology** | first row with `N ≥ 256`: `offplane ≤ 0.05`, `R/R_disc ∈ [0.7, 1.3]`, `d_centres_min ≤ 1.2 µm` | filament or stacking |
| **I0c.4 fill (the T4/T5 door)** | some row has `door_contact = 1`; `N_end ≥ 256`; `wall_leak = 0` on every row | never touched the door; `N_end = 0` EMPTY; drained below 256; closed-face leak |
| **I0c.5 chassis regression** | `check_job2/3/3b/3c` exit 0; Job 3b CSVs SHA256-identical to I0.5 freeze | any lock moves |
| **I0c.6 honesty** | standing states AHL/Hill/NARMA OFF; `W = 20`; founder still `(50, 80)` | occupancy claim; silent `W` or founder hunt |

I0.5 freeze hashes: see the I0 section of this file.

**Smoke:** stop at `N = 32`. Expect I0b smoke numbers: `offplane = 0`,
`R/R_disc = 1.44`, `y_max ≈ 83`, `door_contact = 0`, no spill. If
filament (`R/R_disc > 2`), **stop**.

If I0c.4 FAILS EMPTY after `door_contact = 1`: that is the T4/T5
question answered. Standing, stop, do not hunt `W`.

If I0c.4 FAILS because `door_contact` still never fired: standing, stop,
do not move the founder in this file.

## OFF in I0c

Same as I0b: nutrient PDE, AHL PDE, Hill `R,L`, LuxI, Danino 4-ODE,
NARMA / ridge / 408-D, motility, Biselli death, acid clamp-death, second
species, chemotaxis, metabolism ODE, Brownian, `F_s`, stochastic `P(ℓ)`,
the T4 bus, `L_n = 20` corridor, HybridDish `p_removal`.

## Later jobs (not this freeze)

I1 NutrientPocket, I2 HillOnRods, I3 ChassisHybridmm. A `W = 3 µm` job
is a new PROTOCOL section, only if I0c scores the door at `W = 20` as
EMPTY. Do not start I1 / I2 from I0c.
