# A0 protocol — ideal AC source on C1c_FILLED_POCKET (NOT_FIG4B)

**Gate:** A0 ideal prescribed molecular current.  
**Status label:** `A0_IDEAL_SOURCE`. **NOT_FIG4B.**  
**Device:** `C1c_FILLED_POCKET` only.  
**Organism (bacteria-ON arms):** `DaninoSI_OccupiedDF` (Object B). **NOT_FIG4B.**  
**Not ported:** Object A `Danino2010_Fig4b_bulk_twin`.  
**Not this gate:** Fig. 4b; Lentini TX–TL; A1 vesicle/Hill-gate; A2 NARMA/ridge/CHARC; W0; HybridDish L.

C1 packed-spatial standing remains **FAIL**. Do not rewrite it as PASS.  
Do **not** place this AC on `C1_OPEN_DILUTE`.  
C1b island standing remains **PASS** on the C0-equivalent voxel only. Island PASS is not the A0 device.  
C1c filled-pocket standing remains **PASS**. A0 starts from that PASS, still **NOT_FIG4B**.

Claim freeze:
[`examples/PocketDish/CLAIM_FREEZE_DANINO_SI.md`](../PocketDish/CLAIM_FREEZE_DANINO_SI.md).  
Pre-build design §5.1 / §8–9 / §14 A0:
[`examples/PocketDish/AC_DANINOPOCKET_PREBUILD_DESIGN.md`](../PocketDish/AC_DANINOPOCKET_PREBUILD_DESIGN.md).  
C1c standing:
[`examples/PocketDish/C1C_FILLED_POCKET_STANDING.md`](../PocketDish/C1C_FILLED_POCKET_STANDING.md).

**frozen_before_traces:** true  
**Frozen:** 2026-08-27, after C1c filled-pocket PASS, before any A0 integration.

This protocol answers, on an identity frozen **before traces**:

> Can a named, immobilized, one-way releasing module deposit a prescribed
> molecular current \(J(t)\) conservatively into the C1c filled-pocket AHL
> field, with commanded mass matching the N0 source ledger, a bacteria-OFF
> field predicted by N0/C1c transport, and Object B still occupying its
> C1c class when \(J=0\)?

Not a digital twin. Not Lentini TX–TL. Not a Fig. 4b repair.

Do not shrink \(\mu=0.32\)–\(0.40\). Do not `TIME_ADJ=60`. Do not retune
\(\alpha,\tau,\gamma^\ast,k_1,d,D,D_{1,\mathrm{spatial}}\). Do not start
A1/A2/W0 from a FAIL. Do not implement A0 on `C1_OPEN_DILUTE`.

---

## Named objects (frozen now)

Keep (unchanged):

| Object | Role | Status |
|---|---|---|
| `C1c_FILLED_POCKET` | occupancy identity; **A0 device** | PASS. NOT_FIG4B |
| `C1b_ISLAND` | C0-equivalent voxel occupancy | PASS (unchanged). Not the A0 device. NOT_FIG4B |
| `C1_OPEN_DILUTE` | empty chip+bus | FAIL / `NO_PERIOD`. **Forbidden A0 device.** Do not re-run. Do not place the AC here. NOT_FIG4B |
| `C1_PACKED_SPATIAL` | C1 standing | FAIL (unchanged). NOT_FIG4B |

New:

| Object | Role |
|---|---|
| `A0_IDEAL_SOURCE` | named immobilized one-way module; \(J(t)=J_{\max}u(t)\); no internal AC ODE |

Payload species: **extracellular AHL** (same \(H_e\) as Object B).  
Not IPTG. Not an unnamed “signal.”

---

## Device (frozen; copied from C1c, not retuned)

| Symbol | Frozen value |
|---|---|
| Device | `C1c_FILLED_POCKET` |
| Pocket | \(100\times100\times1.65\) µm |
| Voxel | \(dx=dy=2\) µm, \(dz=1.65\) µm, \(V_{\mathrm{vox}}=6.6\) µm³ |
| Pocket grid | \(50\times50\times1\) (all fluid) |
| \(N_{\mathrm{vox}}\) | 2500 |
| \(V_{e,\mathrm{pocket}}\) | **16500** µm³ |
| \(N\) | **11000** (bacteria-ON arms) |
| \(d\) | **0.5** exactly |
| `D1_spatial` | **800**, `SI_TAKEN_SPATIAL` |
| Open-edge BC | **no-flux** (bus not constructed; mouth is a wall) |
| Bus | chemical solids (not in \(V_e\)) |
| \(\mu\) | **0.40** on pocket fluid only |
| Mechanics | OFF |
| Clock | SI-scaled. \(\tau=10\). No `TIME_ADJ` |
| Kernel | N0 `BSimTransportField`. Signed `transferQuantity` for membrane. `addQuantity` for AC sources |
| Object B | unchanged. No \(\alpha,\tau,\gamma^\ast,k_1,d,D,D_1\) retune |

This is **not** `C1_OPEN_DILUTE`. This is **not** Fig. 4b.

---

## AC placement (frozen now)

One immobilized point in the pocket interior. Deposit into the
**containing voxel only**. Not in the bus (bus is not constructed).
Not on the mouth (\(j=49\), \(y\in[98,100]\)).

| Item | Frozen value |
|---|---|
| World \((x,y,z)\) | \((50,\ 50,\ 0.825)\) µm |
| Containing voxel \((i,j,k)\) | \((25,\ 25,\ 0)\) |
| Voxel world box | \([50,52)\times[50,52)\times[0,1.65]\) µm |
| Count | 1 |
| Mobility | immobilized |

Far-voxel probe (field delay; frozen now): \((i,j,k)=(0,0,0)\),
world centre \((1,\ 1,\ 0.825)\) µm. Interior corner, not the mouth.

---

## Ideal source law (frozen now; no internal AC ODE)

For \(t\ge 0\):

\[
J(t)=J_{\max}\,u(t)
\quad\text{(quantity per SI-time; field “molecule” units)}
\]

\(u(t)\in\{0,1\}\) is a predeclared rectangular command. No TX–TL, no
Hill gate, no leakage ODE, no payload ODE, no vesicle.

### Mass-budget \(J_{\max}\) (ENGINEERING; frozen before traces)

C1c kick mass on the pocket:

\[
\Delta M_{\mathrm{kick}}=H_e(0)\,V_{e,\mathrm{pocket}}=0.05\times 16500=\mathbf{825}.
\]

Unit pulse duration \(\tau_{\mathrm{pulse}}=\mathbf{1}\) SI-time
(short vs Object B \(T\sim 61\); 1000 scheduler steps at `rk_dt=0.001`).

\[
J_{\max}=\frac{\Delta M_{\mathrm{kick}}}{\tau_{\mathrm{pulse}}}=\mathbf{825}.
\]

A unit pulse therefore adds \(\Delta M=825\), comparable to the C1c
kick mass, deposited into the AC voxel only. **Not** fitted to later
\(I\) traces. **Not** fitted to NARMA. Do not retune \(J_{\max}\) after
seeing \(I\).

Quantity units are the SI concentration units already used by Object B
(\(H_e=0.05\) kick) times µm³. They are the N0 field ledger units.

### Rectangular commands (frozen now)

Exact overlap on each scheduler interval \([t,t+\Delta t)\):

\[
\Delta M_{\mathrm{step}}=J_{\max}\cdot\big|[t,t+\Delta t)\cap[t_{\mathrm{on}},t_{\mathrm{off}})\big|.
\]

Deposit that \(\Delta M_{\mathrm{step}}\) with `addQuantity` on voxel
\((25,25,0)\) in the scheduler **flux-deposition phase**. Do not deposit
inside `bacterium.action`. Do not use a lagged half-step.

| Command | \(u(t)\) | \(t_{\mathrm{on}}\) | \(t_{\mathrm{off}}\) | \(\Delta M_{\mathrm{command}}\) |
|---|---|---|---|---|
| `U_PULSE` | 1 on the half-open window | 0 | 1 | \(825\) |
| `U_STEP` | 1 on the half-open window | 0 | 8 | \(6600\) |
| `U_ZERO` | 0 | — | — | \(0\) |

No pulse train. No other \(u(t)\).

---

## Ledger (frozen now)

\[
\Delta M_{\mathrm{command}}=\int_0^{t_{\mathrm{end}}} J(t)\,dt
\]

\[
\Delta M_{\mathrm{field\_sources}}=\text{N0 } \texttt{sourceAdded}
\]

Gate:

\[
\frac{\bigl|\Delta M_{\mathrm{command}}-\Delta M_{\mathrm{field\_sources}}\bigr|}{M_\ast}\le 10^{-3}.
\]

\(M_\ast=\max(|\Delta M_{\mathrm{command}}|,|\Delta M_{\mathrm{field\_sources}}|,825,10^{-15})\).
Tighter is allowed; do not loosen after traces.

Initial `AHL_KICK_005` (bacteria-ON) is initial mass, not a runtime
source (`ledger.reset()` after `setConc`). Membrane uses signed
`transferQuantity` (not a source). AC `addQuantity` **is** a source
and **must** enter the AHL residual on bacteria-ON arms:

\[
R=M(0)+\mathrm{synth}+\Delta M_{\mathrm{field\_sources}}
-\mu_{\mathrm{loss}}-\gamma_{H,\mathrm{loss}}-\mathrm{outlet}-\mathrm{boundary}-M(t).
\]

Bacteria-OFF residual is the N0 transport residual (initial 0).

---

## N0 / C1c field oracle (frozen now)

Bacteria-OFF field response must match a **separate N0-only run** on
the same C1c pocket mask, \(D_1=800\), \(\mu=0.40\), no-flux faces,
same rectangular pulse deposited in `depositFluxes` by the **inline
frozen numbers** \(J_{\max}=825\), \((i,j,k)=(25,25,0)\), **not** by
calling `A0IdealSource`. Same scheduler, same `rk_dt`.

That oracle is the predeclared transport prediction. It is not an
analytic 2-D cosine series (optional extra; not required).

Compare, at every recorded sample:

- \(H_e\) at the AC voxel
- \(H_e\) at the far voxel
- planar centroid \((c_x,c_y)\) of remaining field mass

**PASS** if

\[
\max_t\frac{\bigl|H_e^{\mathrm{A0}}-H_e^{\mathrm{N0}}\bigr|}{\max(H_e^{\mathrm{N0}},10^{-15})}\le 10^{-9}
\]

at both voxels, and centroid distance \(\le 10^{-9}\) µm, for each
field arm vs its oracle.

---

## Predeclared arms

Do **not** add `A0` on `C1_OPEN_DILUTE` as an identity arm.

| id | bacteria | \(u\) | \(t_{\mathrm{end}}\) | `sample_dt` | required |
|---|---|---|---|---|---|
| `A0_FIELD_IMPULSE` | OFF (no Object B integration) | `U_PULSE` | 10 | 0.1 | mass replay; He/centroid vs N0 impulse oracle |
| `A0_FIELD_STEP` | OFF | `U_STEP` | 10 | 0.1 | mass replay; He vs N0 step oracle; delay/monotonicity below |
| `A0_J0_BACTERIA` | Object B ON, `AHL_KICK_005`, \(\mu=0.40\) | `U_ZERO` | 1000 | 0.5 | OSC in the C1c class: \(\|T-63.5833\|/63.5833\le 20\%\). Report \(T\) vs 61.25 (C1c) and 63.58 (Object B). Proves the AC is off |
| `A0_J_PULSE_BACTERIA` | Object B ON, same kick and \(\mu\) | `U_PULSE` | 1000 | 0.5 | He near AC rises vs `A0_J0_BACTERIA`; mass ledger closes. Report colony-mean \(I\) and local \(I\) at the AC voxel. **Not** a NARMA pass/fail. Do not retune \(J_{\max}\) |

`rk_dt=0.001` on every arm (C1c clock). Mechanics OFF.

### A0_FIELD_IMPULSE extra (frozen)

After \(t_{\mathrm{off}}=1\), \(H_e\) at the AC voxel must fall (diffusion
plus \(\mu=0.40\)). Far-voxel \(H_e\) must rise above its \(t=0\) value
at some recorded \(t\in(0,10]\).

### A0_FIELD_STEP extra (frozen)

While the source is on (\(0<t\le 8\)): \(H_e\) at the AC voxel exceeds
\(H_e\) at the far voxel. Time to 10% of each voxel’s recorded max
satisfies \(t_{10\%,\mathrm{AC}} < t_{10\%,\mathrm{far}}\) (diffusion
delay, not a fudge).

### A0_J_PULSE_BACTERIA extra (frozen)

At the first sample with \(t\ge 1\), \(H_e\) at the AC voxel is strictly
greater than the same sample on `A0_J0_BACTERIA`. That is the only
\(I\)-adjacent predeclared claim. Colony-mean and local \(I\) are
**reported**. Do not fit \(J_{\max}\) to them.

---

## PASS / FAIL

A0 **PASS** only if all four arms hold with **no** Object B / `D1_spatial`
/ \(d\) / \(\mu\) / \(J_{\max}\) / placement change after traces, every
log/plot/CSV/standing line says **NOT_FIG4B** and **A0_IDEAL_SOURCE**,
and the AC was not run on `C1_OPEN_DILUTE`.

FAIL: mass mismatch, field not predicted by the N0 oracle, Object B
retuned, J=0 not OSC in the 20% band, or the AC placed on the open bus
chip.

If A0 PASSes: A1 may start only on `C1c_FILLED_POCKET`, still
**NOT_FIG4B**, still not Fig. 4b, still not `C1_OPEN_DILUTE`. A1 is
**not** calibrated. A2/W0 may not start.

If A0 FAILs: stop. Do not start A1/A2/W0. Do not retune Object B. Do
not move the AC onto `C1_OPEN_DILUTE`.

---

## Commands

```
ant a0-ideal
python examples/BSimReservoirPlanDaninoPocketA0/check_a0.py
```

Standing: `examples/PocketDish/A0_IDEAL_SOURCE_STANDING.md`.  
Must say `A0_IDEAL_SOURCE`, **NOT_FIG4B**, device = `C1c_FILLED_POCKET`,
that `C1_PACKED_SPATIAL_STANDING.md` remains FAIL, that
`C1B_ISLAND_STANDING.md` remains island PASS, and that
`C1C_FILLED_POCKET_STANDING.md` remains filled-pocket PASS.
