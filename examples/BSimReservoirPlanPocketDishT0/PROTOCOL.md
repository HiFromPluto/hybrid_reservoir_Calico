# PocketDish-T0 — transport / occupancy (job protocol)

Frozen 2026-08-19. Job-level only. Architecture freeze:
[`examples/PocketDish/PROTOCOL.md`](../PocketDish/PROTOCOL.md).
Geometry: [`examples/PocketDish/GEOMETRY_FREEZE.md`](../PocketDish/GEOMETRY_FREEZE.md).
Citations: [`examples/PocketDish/CITATION_DOSSIER.md`](../PocketDish/CITATION_DOSSIER.md).
Boundary names and “no labels in the screen”:
[`examples/HybridDish/design_space_preflight/PROTOCOL.md`](../HybridDish/design_space_preflight/PROTOCOL.md)
§4–§6 (E0.2). Through-flow 8 µm/s occupancy **DEAD**: SweepS5 / E0.3.
Chemical ≠ reporter wipe: WashoutReset.

This file does **not** authorize NARMA, Mackey–Glass, waveform AUC,
ridge, a λ grid, week-3 population Java, PocketOsc, PocketRect, vesicles,
TX–TL, Danino 4-ODE, glucose, acid bias, a second strain, Grober pucks,
Stokes rotlets, or two-way \(h(P)\).

Do not edit any `GATE_EVIDENCE.md`. Do not retune HybridDish \(K,n,\tau_R,\tau_L\)
or HybridDish claim \(J_{\max}\). Do not copy HybridDish \(J_{\max}=1.28\times10^8\)
into this pocket. Do not promote a new claim dish.

## Scientific question

On HybridDish clocks (\(D=159\) µm²/s, \(k=0.0033\) s⁻¹,
\(\tau\approx 303\) s), a 100 µm pocket has \(L^2/D\approx 63\) s, so
**diffusion to the open edge is faster than decay**. Does side-channel
exchange (bus 3 µm/s *past* the open edge) leave the Hill occupied, or
does the interface drain \(C\) below occupancy?

Hypothesis (predeclared, before occupancy):

- `CLOSED_NOFLUX` is **ALIVE** at \(u=0.5\) (closed-pocket well-mixed
  \(C_{\mathrm{ss}}=K\)).
- `THROUGH_8` is **DEAD** (negative control; SweepS5 / E0.3 class).
- `OPEN_BUS_3` is the **unknown**.
- Corner AHL hotspots on no-flux walls are **artifacts to report**,
  not to cure by retuning \(K\).

## Path

Preferred: **reduced 2-D field + Hill/\(L\) surrogate** (E0.2 style),
this package. Living BSim occupancy, seed **101** only, **only if**
reduced `OPEN_BUS_3` is not obviously `SATURATED` from the closed-form
mass check below. If reduced `OPEN_BUS_3` is `DEAD`: that is the result;
do not retune \(K\), \(k\), or \(J_{\max}\); do not start NARMA; still
run living `CLOSED_NOFLUX` seed 101 **if cheap**. Living pocket+bus
Java is **not cheap** (new mask / bus field). Reduced `CLOSED_NOFLUX`
is the \(J_{\max}\) unit sanity. If that arm is `DEAD`, fix
**units/conversion only**, not \(K\).

Do not start week-3 population Java from this job.

## Frozen dish (PocketDish-A + PocketHill)

Copied clocks: HybridDish \(D,k,K,n,\tau_R,\tau_L\). **New** geometry
and PocketDish \(J_{\max}\) from the architecture PROTOCOL.

| Quantity | Value | Class |
|---|---|---|
| Pocket | \(100\times 100\times 10\) µm | Danino-class square + HybridDish height ENGINEERING |
| Open edge | entire **+y** face | ENGINEERING handedness |
| Bus | along \(x\), width **80** µm, length **400** µm | ≥400 µm freeze; 400 chosen |
| Pocket attachment | pocket \(x\in[150,250]\) µm, \(y\in[0,100]\); bus \(y\in[100,180]\), \(x\in[0,400]\) | ENGINEERING, pocket centred on bus |
| Pocket flow | 0 | frozen |
| Grid | \(dx=dy=5\) µm, one \(z\) box (depth 10 µm) | numerics |
| `dt` (declared living) | 0.05 s | HybridDish |
| \(D\) | 159 µm²/s | HybridDish |
| \(k\) | 0.0033 s⁻¹ | HybridDish ENGINEERING |
| \(K,n,\tau_R\) | 1.6 µM, 2, 15 s | HybridDish; do not retune |
| \(\tau_L\) | 1500 s | HybridDish; do not retune |
| Conversion | 1 µM = 602.2 molecules/µm³ | DERIVED |
| \(J_{\max}\) | \(6.36\times 10^5\) molecules/s at pocket centre | predeclared; `ENGINEERING_DESIGN` |
| Command | \(J=J_{\max} u\) | one-way `addQuantity`; no vesicles |
| AC | pocket centre (50, 50, 5) µm in pocket coords | ENGINEERING |

Closed-pocket well-mixed algebra (must hold before any occupancy
interpretation):

\[
C_{\mathrm{ss}}(u=1)=2K=3.2~\mu\mathrm{M},\qquad
J_{\max}=k\,C_{\mathrm{ss}}\,V,\qquad
V=10^5~\mu\mathrm{m}^3.
\]

At \(u=0.5\), leak-free \(C_{\mathrm{ss}}=K=1.6\) µM. Open-pocket
occupancy is then a measurement of **interface leak**. Do not raise
\(J_{\max}\) if open arms are `DEAD`. Do not lower \(K\) if they are
`SATURATED`.

Reduced model: cell-centred finite volume on the 5 µm grid; centred
diffusive fluxes; first-order upwind advection; SciPy `solve_ivp`
BDF with `rtol=1e-7`, `atol=1e-10`; source discontinuities are hard
segment edges. Mass identity (molecules):

`residual = initial + injected − remaining − decay_loss − boundary_loss`

Do **not** mix STEP_ON leftover into SINGLE_PULSE accounting (E0.2
warmup/NARMA lesson). Each drive starts from \(C=R=L=0\). Unresolved
\(|\mathrm{residual}| > 1\%\) of `initial+injected`, or \(C < -10^{-9}\) µM,
is a numerical failure. Particle rules (mirror / open-edge remove) are
**not** chemical sinks.

This surrogate has **no cells**. \(N(t)\) is `NA` on reduced rows.

## Predeclared arms (do not add after occupancy)

| Arm | Chemical BC | Flow | Role |
|---|---|---|---|
| `CLOSED_NOFLUX` | four no-flux walls; no bus | 0 | occupancy sanity |
| `OPEN_ABSORBING` | no-flux on \(x=0,x=100,y=0\); Dirichlet \(C=0\) on pocket \(+y\) (E0.2 ghost, face \(C=0\)) | 0 | harsh leak diagnostic |
| `OPEN_BUS_3` | three pocket no-flux walls; pocket \(+y\) coupled by diffusion to an explicit bus; bus inlet \(C=0\); conservative bus outlet; bus far walls no-flux | bus \(v_x=3.0\) µm/s; pocket \(v=0\) | **backbone geometry** |
| `THROUGH_8` | no-flux on pocket \(x\) walls; advective inlet \(C=0\) at \(y=0\); conservative outflow at \(y=100\) | **8 µm/s through the pocket in \(+y\)** | negative control |

No extra Robin \(h\), extra \(J_{\max}\), extra \(K\), extra flow, or
extra geometry after the first occupancy number. Millimetre-dish
\(h=0.159\) µm/s is not imported.

`THROUGH_8` is through-pocket, not bus-past-the-edge. If it is `ALIVE`,
do **not** celebrate: check that mean \(v_y\) in the pocket is 8 µm/s
and that outlet flux is positive. Transit \(L/v=12.5\) s.

Occupancy means are **pocket voxels only** (bus AHL is reported
separately on `OPEN_BUS_3` and does not enter `mean_R`).

## Predeclared drives (no task labels)

Command range \([0,1]\). No NARMA file. No sequence labels.

1. **`STEP_ON`:** \(u=0.5\) held for **\(T=9000\) s** (6 \(\tau_L\);
   HybridDish warmup 18 000 s is allowed; this job declares the shorter
   plateau). Plateau: terminal pocket `mean_R` within 5% of the mean of
   `mean_R` over the last 10% of time (8100–9000 s). Occupancy flag is
   evaluated on the terminal pocket field. \(r(R,u)\) is **NA** (\(u\)
   does not vary).
2. **`SINGLE_PULSE`:** \(u=1\) for **75 s**, then zeros for **1500 s**
   (5 \(\tau_{\mathrm{AHL}}\), 1 \(\tau_L\)). Report AHL and \(L\)
   relative to their own peaks at pulse-end+1500 s. Expect chemical ≠
   reporter (WashoutReset). \(r(R,u)\) is Pearson of pocket `mean_R(t)`
   against \(u(t)\) on this drive only.

## Occupancy flags (not NRMSE)

| Flag | Rule |
|---|---|
| `ALIVE` | pocket `mean_R ≥ 0.05` and not `SATURATED` |
| `DEAD` | pocket `mean_R < 0.05` |
| `SATURATED` | pocket `mean_R > 0.95` |

Also print at plateau / pulse: pocket-mean AHL, \(R\), \(L\);
corner/interior AHL ratio; \(L^2/D\), Damköhler \(k L^2/D\), Pe;
mass budget; `THROUGH_8` flow check.

**Corners:** pocket SW and SE (two-wall no-flux) plus the two open-edge
side-wall voxels NW and NE. Report closed-end ratio
\(\langle C_{\mathrm{SW}},C_{\mathrm{SE}}\rangle / C_{\mathrm{interior}}\)
and all-four-corner ratio. Interior: pocket cells at least two cells
from every pocket wall. Do not retune \(K\) if corners are hot.

## Closed-form mass check (before living)

Well-mixed closed pocket at \(u=0.5\): \(C=1.6\) µM, \(H=0.5\),
so `SATURATED` from algebra would mean **wrong \(J_{\max}\) units**.
Leak can only lower \(C\). `OPEN_BUS_3` cannot be `SATURATED` unless
the conversion is wrong. If reduced `OPEN_BUS_3` is `DEAD`: **stop**.

## Gates that are not this job

No NRMSE, AUC, NARMA target, ridge, or λ. The checker **must refuse**
NARMA paths and must not import `BSimReservoirPlanNarma10b`.
Do not write Overall PASS/FAIL on a task. Overall lines:

- `OCCUPANCY_<ARM>=ALIVE|DEAD|SATURATED` on `STEP_ON`
- `TRANSPORT_MODEL_STATUS=VALIDATED_FOR_SCREENING|NOT_VALIDATED`
  from **this geometry** (mass, non-negativity, time convergence,
  closed-form \(J_{\max}\)). Do **not** block the screen on
  millimetre-dish NARMA voxel files.

Convergence (optional E0.2 port): `CLOSED_NOFLUX` + `SINGLE_PULSE`,
tighter BDF (`rtol=1e-9`, `atol=1e-12`) vs primary. Peak/mean \(C,R,L\)
and integrated mass terms within 1%. Failure → `NOT_VALIDATED`; still
print occupancy as `NON_DECISION_GRADE` if mass also fails.

## Outputs

- `results/TRANSPORT_SCREEN.md`
- `results/transport_screen.csv`
- occupancy maps (AHL, \(R\), \(L\)) as images
- SHA-256 of configs and checker

## Commands

```
python examples/BSimReservoirPlanPocketDishT0/check_pocketdish_t0.py --theory
python examples/BSimReservoirPlanPocketDishT0/test_mass_budget.py
python examples/BSimReservoirPlanPocketDishT0/run_transport_screen.py
python examples/BSimReservoirPlanPocketDishT0/check_pocketdish_t0.py
```

`--smoke` on the runner integrates 50 s only and must not be scored as
occupancy evidence.

## Stop list

- Do not add arms, \(h\), \(J_{\max}\), or \(K\) after occupancy
- Do not retune HybridDish clocks to occupy the pocket
- Do not start NARMA if `OPEN_BUS_3` is `DEAD`
- Do not call this an exact Danino/Prindle blueprint
- Do not rewrite WashoutReset / E0.3 / Narma10b results
- Do not start week-3 population Java
