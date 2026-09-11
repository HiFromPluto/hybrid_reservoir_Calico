# Job 2 protocol — isolated *E. coli* rod

Single cell. Cited growth + shape + division. No colony packing, no QS,
no chemotaxis ODE, no genome-scale dump. Not a reservoir task.

## Equations (TAKEN)

Valdez et al. 2025, eqs. (4)–(5), DOI 10.1038/s42005-025-02078-1:

```
ℓ(t) = ℓ₀ exp( σ · N/(N+κ) · t )     while ℓ < ℓ_div
σ    := ν · log₂(ℓ_div / ℓ₀)
```

Discrete exact step at constant bath `N` (same closed form):

```
ℓ ← ℓ · exp( σ · N/(N+κ) · Δt )
```

Warren et al. 2019, DOI 10.7554/eLife.41093, Materials and methods:

```
λ = λ_S · C / (C + K_S)
```

Job 2 identifies `ν = λ_S` and `κ = K_S`, and `C = N` as a uniform
saturating bath (no nutrient PDE).

Geometry (Warren Methods, fixed in all their simulations):

- diameter `w₀ = 1 µm` → hemisphere radius `0.5 µm`
- `ℓ_div = 3 µm`
- `ℓ_div = 2 ℓ₀ + w₀` ⇒ `ℓ₀ = 1 µm`
- aspect `ℓ_div : w₀ = 3 : 1`

Division placement (Warren Fig. 12B): the mother’s hemisphere centres
become the hemisphere centres of the two daughters; each daughter has
cylindrical length `ℓ₀`.

## ENGINEERING (named, frozen)

- Division is a **deterministic threshold** `ℓ ≥ ℓ_div`, not Valdez’s
  stochastic arctan `P(ℓ)`. Objective: seed-reproducible isolated-cell
  test, not “look alive.”
- Daughter length fluctuation is **zero** (Warren allows a small
  `ℓ_ran`; we freeze 0).
- Isolated cell: no Hertzian forces, no Brownian motion, no motility.
- Time step `Δt = 1 s` (growth timescale is hours).

## Frozen numbers

Printed at every run as `TAKEN` or `ENGINEERING`. See
`ChassisParameters.java`.

| Symbol | Value | Class | Source |
|---|---|---|---|
| `λ_S` = `ν` | `1.0 h⁻¹` | TAKEN | Warren glucose MM |
| `K_S` = `κ` | `0.02 mM` (`20 µM`) | TAKEN | Warren (Monod 1949 via Warren) |
| `N` (bath) | `0.5 mM` | TAKEN | Warren `C_s` |
| `w₀` | `1 µm` | TAKEN | Warren |
| `ℓ_div` | `3 µm` | TAKEN | Warren |
| `ℓ₀` | `1 µm` | TAKEN | Warren `ℓ_div = 2ℓ₀ + w₀` |
| `Δt` | `1 s` | ENGINEERING | integrator |
| division | `ℓ ≥ ℓ_div` | ENGINEERING | deterministic sizer |
| `ℓ_ran` | `0` | ENGINEERING | no length noise |

Expected first division time at this bath:

```
Monod = 0.5 / (0.5 + 0.02) = 0.961538…
T_div = ln(2) / (λ_S · Monod) ≈ 2595.1 s ≈ 43.25 min
```

Warren’s `λ_S = 1.0 h⁻¹` ballpark without Monod is `ln 2 / λ_S ≈ 41.6 min`.
The gate uses the Monod-corrected analytic `T_div`, not a fitted movie.

## Box

- Domain 100 × 100 × 10 µm (kinematics only; one cell does not fill it)
- `N = 1` at `t = 0`, replication on, death off
- Seed 101 (unused by the ODE; recorded)
- Export 2.5 h (`t = 9000 s`), log every 10 s
- Completeness: last CSV row `t = 9000`

CSV `results/job2_seed101/job2_timeseries.csv`:

```
t_s;N;founder_L_um;founder_L_analytic_um;L_min_um;L_max_um;radius_um;nutrient_mM;n_negative
```

## Gates

1. `n_negative = 0` on every row (length, radius, nutrient).
2. Before first division, `|founder_L − analytic| / ℓ₀ ≤ 1e-6`.
3. First time `N = 2` is within one log step (10 s) of analytic `T_div`.
4. On that row, `L_min` and `L_max` are within `0.02 µm` of `ℓ₀`.
5. Radius is `0.5 µm` on every row.
6. `N` is in `{1,2,4,8}` (binary fission, no death). Completeness `t = 9000`.

Do not retune `λ_S` or `K_S` if a gate fails.

Job 2 standing: Hertzian **off**. Job 3 does not retune these numbers.

# Job 3 protocol — Hertzian excluded-volume rods

Same capsule and Valdez (4)–(5) elongation as Job 2. Adds Valdez (6)–(10)
packing. Uniform saturating bath still `N = 0.5 mM` (no nutrient PDE).
Closed box so growth must push. Deterministic seed 101.

## Force law (non-negotiable)

`RelaxationMover` is **not** Hertzian. It integrates
`BSimCapsuleBacterium.computeNeighbourForce`:

```
repulsion = 0.4 k_cell (2R − d)^{2.5}
```

with uncited `k_cell`. Job 3 does not call it.

Implemented on the capsules, Valdez et al. 2025 DOI
10.1038/s42005-025-02078-1 eqs. (8)–(10):

```
δ_cc = d0 − d    if d0 > d else 0     (d0 = w0)

F_cc n = ( (2/3) k_cc √d0 δ_cc^{3/2} − γ_n M_eff δ_cc (v_cc · n) ) n

F_cc τ = −min( γ_t M_eff √δ_cc (v_cc · τ) ;
               (2/3) μ_cc k_cc √d0 δ_cc^{3/2} ) τ
```

Cell–wall uses (10) with `k_ac`, `μ_ac`, and the `2√2/3` prefactor.
Overdamped (6)–(7): `F_i = η_t ú_i`, `T_i = η_t θ̇_i` (RodCellVSPhage
`Integrate.cpp` uses `η_t` for both). Tirado drag with
`eta = 1e6 * μ_liq` as in that BUILD file (μm, hours) — **TAKEN**, not
raised after overlap. Each growth step solves contacts **quasi-statically**:
elastic Valdez (8)/(10) only (`v = 0`), CFL `0.02 µm`/substep, until
overlap ≤ named residual. Job 2 `Δt = 1 s` is longer than contact time;
this is not a 1 s ballistic step and does not raise `k_cc`.

`F_s` OFF. Brownian OFF. Phage OFF. Nutrient field `N(x,t)` OFF.
Stochastic `P(ℓ)` OFF. Uncited `k_ov` OFF.

## Frozen numbers (printLedger)

| Symbol | Value | Class | Source |
|---|---|---|---|
| `k_cc` | `3.0e4` | TAKEN | Warren Table S4 via RodCellVSPhage `Constants.cpp` |
| `k_ac` | `3.0e4` | TAKEN | same table `k_wc` |
| `μ_cc` | `0.1` | TAKEN | Warren Fig. 9 / Table S4 |
| `μ_ac` | `0.8` | TAKEN | Warren Fig. 9 `μ_ca` |
| `γ_t` | `1.0e4 µm^{-1} h^{-1}` | TAKEN | Warren Fig. 9 `γ_cc,t` |
| `γ_n` | `5.0e2 µm^{-1} h^{-1}` | TAKEN | RodCellVSPhage Table S4 |
| `μ_liq` | `1.00160e-3 Pa s` | TAKEN | water; Valdez (12) |
| `eta` factor | `1e6 * μ_liq` | TAKEN | RodCellVSPhage `Integrate.cpp` |
| Job 2 `λ_S`, `K_S`, `w_0`, `ℓ_0`, `ℓ_div` | locked | TAKEN | do not retune |
| contact CFL | `0.02 µm`/substep | ENGINEERING | quasi-static residual solve |
| `CLUSTER_BOX` | `8 × 6 × 4 µm` closed | ENGINEERING | named; growth must push |
| `δ_cc(0)` two-body | `0.20 µm` | ENGINEERING | known overlap |
| two-body `t_end` | `7200 s` | ENGINEERING | ~2 h at paper mobility |
| gap residual | `0.02 µm` | ENGINEERING | leftover overlap allowed |
| cluster `t_end` | `6000 s` | ENGINEERING | ≥2 divisions at Job 2 bath |
| overlap cap | `0.25 µm` | ENGINEERING | frozen before seeing overlap |

Do not raise `k_cc` after seeing overlap. Do not use `T7Plaque.txt`
(`k_cc = 1e6`). Do not enlarge drawn radius. Do not add “motility off
when overlapping.”

## Tests

CSV directory `results/job3_seed101/`. Seed column on every row.

### 1. Two-body

`twobody_timeseries.csv`

Seed two parallel capsules with `δ_cc = 0.20 µm`. No growth. After
relaxation, min surface gap `≥ −0.02 µm` and each centre displaces
`≥ 0.05 µm`. Completeness last row `t = 7200`.

### 2. Growing cluster

`cluster_timeseries.csv`

Seed two rods in `CLUSTER_BOX` (`8 × 6 × 4 µm`, closed) with seed
overlap `0.15 µm`. Job 2 bath. Grow through ≥2 divisions (`N: 2 → 4 → 8`).
Max pairwise `δ_cc` stays `≤ 0.25 µm`. No two centres closer than
`w_0 − 0.25 µm`. Completeness last row `t = 6000`.

### 3. Job 2 regression

`isolated_packing_timeseries.csv`: `N = 1` in the Job 2 box, packing
ON, `t_end = 100 s`. Packing force `~ 0`. Length still matches Valdez
analytic.

`python check_job2.py` must still PASS on the Job 2 CSV (packing off
on that run).

If (1) or (2) fail: STOP.

```
compile_and_run.cmd
python check_job3.py
python check_job2.py
```

# Job 3b protocol — spatial nutrient field + local Monod elongation

Same capsule, Valdez (4)–(5) elongation, deterministic division, and Job 3
Hertzian packing. Adds Valdez (16)–(17) / Warren two-domain nutrient PDE
with local `N(x)` at each cell centre. DOI 10.1038/s42005-025-02078-1 /
10.7554/eLife.41093.

## Equations (TAKEN)

Colony Ω_C:

```
∂N/∂t = div(D_C ∇N) − (λ_S/Y) ρ₀ N/(N+K_S)
```

Agar Ω_A:

```
∂N/∂t = div(D_A ∇N)
```

Interface: `D_C ∇N·n_C = D_A ∇N·n_A`, `N` continuous.

BC: **ENGINEERING** — low-index agar faces `i=j=k=0` Dirichlet `N = C_s`;
high faces Neumann; agar is a corner shell (not Warren's full `z<0` slab).

Steady-state acceptance (3b.3): **two** consecutive GS sweeps with
`maxΔN < 1e-9 mM` **and** `flux_rel_err ≤ 0.005` (Dirichlet inflow vs Monod
consumption). Depletion gate: **every** row with `t > 0` must satisfy the flux
tolerance — no warm-up exemption — and `gs_converged = 1` on every row.

These tolerances are **absolute and frozen at `dx = 0.5 µm`**; they are not
grid-invariant. See `JOB3B_STANDING.md` before changing `dx` or the agar pad.

Growth: sample `N(x)` at cell centre → `elongateCited(Δt, N_local)` with
`κ = K_S`, Monod factor `N/(N+K_S)` (Valdez (4)–(5), `ν = λ_S`).

Warren Appendix A1.2.7: quasi-steady solve (`∂N/∂t = 0`) each integrator tick
via Gauss–Seidel (RodCellVSPhage `Nutrients.cpp` pattern). Consumption is
Monod, not a generic decay constant.

## Frozen numbers (printLedger)

| Symbol | Value | Class | Source |
|---|---|---|---|
| `D_C` | `90 µm²/s` | TAKEN | Warren Table S2 |
| `D_A` | `600 µm²/s` | TAKEN | Warren Table S2 |
| `Y` | `0.5 gCDW/g glucose` | TAKEN | Warren Table S2 |
| `ρ_cell` | `0.137×10⁻¹² gCDW/µm³` | TAKEN | Warren Table S3 |
| `ρ₀` | `0.68 ρ_cell` | TAKEN | Warren colony density (v1 constant/voxel) |
| `C_s` | `0.5 mM` | TAKEN | Warren lateral agar |
| `MM→g/µm³` | `1.8×10⁻¹⁶` | TAKEN | `180 g/mol`, `1 L = 10¹⁵ µm³` |
| flux gate | `≤ 0.05` rel. err. | ENGINEERING | inflow vs Monod sink |
| BC geometry | low faces Dirichlet | ENGINEERING | corner agar shell |
| `K_S`, `λ_S`, `ℓ₀`, `ℓ_div`, `k_cc` | Job 2/3 locked | TAKEN | do not retune |
| grid `Δx` | `0.5 µm` | ENGINEERING | resolve CLUSTER_BOX |
| agar pad | `2` voxels = `1.0 µm` | ENGINEERING | VOXEL COUNT: shell = pad x dx; rescale pad on any dx change |
| GS increment tol | `1e-9 mM` | ENGINEERING | 3b.3; absolute, frozen at dx=0.5 |
| flux tol | `0.005` | ENGINEERING | 3b.3; inflow vs consumption, every row t>0 |

OFF: phage, QS, chemotaxis, metabolism ODE, Brownian, `F_s`, stochastic `P(ℓ)`.
Do not edit `BSimChemicalField.java`. Do not use `BSimCapsuleBacterium.k_growth`.

## Tests

CSV directory `results/job3b_seed101/`. Seed 101 on every row.

### 1. Uniform-field regression

`uniform_timeseries.csv` — Job 2 box (100×100×10 µm), single cell,
`ρ = 0` (uniform bath mode). Job 2 analytic length and `T_div` unchanged.
Completeness last row `t = 9000`.

### 2. Packing regression

`packing_cluster_timeseries.csv` — `CLUSTER_BOX` with uniform `N ≡ C_s`
(`ρ = 0`) and Job 3 Hertzian. Same gates as Job 3 cluster.
Completeness last row `t = 6000`.

`python check_job3.py` must still PASS on the unchanged Job 3 CSVs.

### 3. Depletion + coupled growth

`depletion_timeseries.csv` — closed `CLUSTER_BOX`, Monod consumption on,
two seed rods. Time series shows `N_local < C_s`, local Monod `<` bath
Monod, founder length below uniform-bath analytic; no negative `N`.
At least one division (`N ≥ 4`); max Hertzian overlap `≤ 0.25 µm`.
Completeness last row `t = 6000`.

CSV columns include `N_local_mM`, `field_min_mM`, `field_mean_mM`.

```
compile_and_run.cmd job3b
python check_job3b.py
python check_job3.py
python check_job2.py
```

# Job 3c protocol — colony dish, radial nutrient limitation

**FROZEN BEFORE IMPLEMENTATION.** Every threshold below is pre-registered.
Do not retune a threshold after seeing a run. If a gate fails, the finding is
the failure.

## Naming

This is **Job 3c**, a continuation of the colony track (3 → 3b → 3c). It is
**not** bootstrap Job 4. `CHAT_BOOTSTRAP.md` reserves Job 4 for
population/communication, and `JOB3B_STANDING.md` says "do not start Job 4 (QS
bootstrap) from this file." QS remains Job 5. Nothing here touches either.

## What 3c is for

Job 3b validated the nutrient PDE **plumbing** but produced only 0.31%
depletion, because the CLUSTER_BOX colony was `0.16` penetration depths across.
3c puts a colony in a geometry where nutrient limitation is a **first-class
biological effect**, and measures the observable 3b structurally could not
produce: a **radial growth-rate gradient**.

## The number that designs the dish (TAKEN + DERIVED)

```
sink at C_s          q = (λ_S/Y) ρ₀ · C_s/(C_s+K_S) / MM_TO_G   = 0.2765 mM/s
penetration depth    δ = √(D_C · C_s / q) = √(90 × 0.5 / 0.2765) = 12.76 µm
```

A colony of radius `R` starves at its centre once `R ≳ δ`. Job 3b sat at
`R/δ = 0.16`. Job 3c targets `R/δ ≥ 1.3`.

Monolayer footprint per cell = `w₀ · (L + 2r)`, mean `3.0 µm²` over
`L ∈ [ℓ₀, ℓ_div]`. At packing fraction `φ = 0.85`:

| N | R (µm) | R/δ | diameter | wall margin |
|---|---|---|---|---|
| 64 | 8.5 | 0.66 | 17.0 | 21.5 |
| 128 | 12.0 | 0.94 | 24.0 | 18.0 |
| **256** | **17.0** | **1.33** | **33.9** | **13.0** |
| 512 | 24.0 | 1.88 | 48.0 | 6.0 |

`N = 256` is the gate target: past one penetration depth with a real wall
margin. `N = 512` is an optional extension and needs an 80 µm chamber — at
6 µm margin the colony feels the walls mechanically and through the BC.

**The table is now MEASURED, not assumed (2026-08-23).** Under
`SYMMETRY_BROKEN` in the `b_z = 1.0 µm` monolayer chamber, the smoke test gives:

| N | R predicted (φ=0.85) | R measured | R/R_disc | offplane |
|---|---|---|---|---|
| 32 | 6.00 | 8.62 | 1.44 | 0.00 |
| 64 | 8.48 | 10.61 | 1.25 | 0.00 |
| 128 | 11.99 | 12.84 | 1.07 | 0.00 |
| **256** | **17.0** | **16.79** | **0.99** | **0.00** |

At the gate target the prediction is good to **1.2%**, `R/R_disc → 0.99`, and
`offplane` is **0.00 at every N** — the monolayer is genuinely enforced. The
φ=0.85 assumption stands; it was simply unreachable until symmetry breaking and
the `b_z` fix were both in place. Measured `R/δ = 1.32` against the assumed
1.33.

## Geometry — ENGINEERING, and different from 3b on purpose

**Chamber `60 × 60 × 1 µm`.** `b_z` equals ONE cell diameter (`2r = 1.0 µm`),
so the monolayer is **enforced by the wall Hertzian**, not hoped for.

**AMENDED 2026-08-23 (was `60 × 60 × 4 µm`).** At `b_z = 4 µm` there is room for
four cell diameters and the colony used them: the morphology smoke test measured
**71% of cells more than one radius off the mid-plane** at N=256. The colony was
compact (`R/R_disc = 0.73`) because it had left the plane, not because monolayer
packing is tight — which would have made Gate 4's radial bins compare in-plane
edge cells against buried interior cells, and made `R/δ` stop meaning what it
was derived to mean. Vertical structure is the Warren Fig. 5 problem that 3c
explicitly scopes out. `b_z = 1.0` is also an integer multiple of `dx = 1.0`.
Do not use `b_z ≈ 1.5`: non-integer in `z` at the working grid, and it still
admits a second layer.

**BC mode `DISH_LATERAL`** — new, additive:

- **Dirichlet `N = C_s` on all four lateral faces** (`x=0, x=b_x, y=0, y=b_y`),
  via an agar rim of `pad` voxels at `D_A`.
- **Flux-free on `z = 0` and `z = b_z`** (chamber floor and ceiling).

### Why not agar-below

An agar-below BC (Warren plate geometry) supplies nutrient through ~1 µm of
agar directly beneath **every** cell, interior and edge alike. The radial
supply path vanishes and **no radial gradient can form** — the gate would
return null for a geometric reason, not a biological one. Lateral supply is
required for the observable to exist.

This is a microfluidic-chamber geometry, with precedent in Melke et al. 2010
(thin rectangular layer, boundary held at fixed `A_e`; CITE-ONLY, geometry
precedent only — no numbers taken) and in the channel devices of Leaman &
Behkam 2018.

**Warren Fig. 5 vertical starvation remains explicitly OUT OF SCOPE.** It needs
a tall stack with agar below, which is a different job.

## Division mode — REQUIRED, not optional

Job 3c runs `EcoliRodCell.DivisionMode.SYMMETRY_BROKEN`, seeded from
`RNG_SEED = 101`.

`COLLINEAR` (the Jobs 2/3/3b default) places daughters exactly on the mother's
axis with no perturbation. In an OPEN chamber nothing ever breaks the founder's
axis, and the colony grows as a **1-D filament**. Measured:

| N | R (µm) COLLINEAR | R if disc | R if chain | d_centres_min |
|---|---|---|---|---|
| 4 | 3.49 | 2.12 | 3.00 | 1.981 |
| 8 | 7.48 | 3.00 | 6.00 | 1.981 |
| 16 | 15.47 | 4.24 | 12.00 | 1.982 |
| 32 | 29.60 | 6.00 | 24.00 | 1.857 |

`R ∝ N`, and `d_centres_min` is pinned at `2·ℓ₀` — lateral Hertzian contact
**never engages**. Past N=32 the filament simply jams against the chamber wall.

Job 3's CLUSTER_BOX (`8 × 6 × 4 µm`, two seeds already in lateral contact) hid
this: a chain hits the wall after ~4 cells, so Job 3's packing gate was testing
**wall-confined** morphology, not open-dish self-organisation.

Mechanism for the fix is cited — Melke et al. 2010 Methods, citing Cho et al.
2007: *"At each cell division we introduce some randomness in order to break the
axial symmetry of the system, giving two daughter cells with slightly different
sizes and imperfect alignment."* Melke's model is non-dimensional, so the
MAGNITUDES (`DIVISION_ANGLE_SD_RAD = 0.05`, `DIVISION_LENGTH_ASYM = 0.05`) are
**ENGINEERING** and remain open to a sweep.

## Grid and cost

`dx = 1.0 µm`, `pad = 1` (rim thickness `pad × dx = 1.0 µm`).

```
grid 62 × 62 × 1 = 3,844 voxels     (3b.3 CLUSTER: 2,520; was 15,376 at b_z=4)
cold acceptance solve  ~340–510 ms at frozen ω/tol
science run            ~35 min wall clock = 23,400 ticks (DT_S = 1.0 s)
```

`DT_S = 1.0 s` is the tick; `LOG_DT_S = 10.0 s` is the CSV row. Do not
confuse row count with tick count. Acceptance costs are **cold** solves
(fresh field); run ticks are **warm** — do not predict one from the other.

Occupancy is `CONSERVATIVE_MASS` (volume-fraction sub-sampling); SOR uses
`ω = 1.96`. Both are ENGINEERING, measured — not ported from 3b.

## Acceptance tolerances — DERIVED AND FROZEN

`JOB3B_STANDING.md` states the 3b.3 acceptance (`NUTRIENT_CONV_TOL_MM = 1e-9`,
`NUTRIENT_FLUX_TOL_REL = 0.005`) is **absolute and frozen at `dx = 0.5`**.
3c uses `dx = 1.0`, so the pair was re-derived from
`results/job3c_acceptance/acceptance_sweep.csv` under
`DISH_LATERAL` + `CONSERVATIVE_MASS` + `ω = 1.96`, with `flux_tol`
non-binding during the study so the residual stayed an independent
measurement.

**Frozen pair (ENGINEERING, `dx = 1.0` only — not transferable):**

| Symbol | Value | Binding case |
|---|---|---|
| `conv_tol` | `1e-10` mM | R=0 (single cell) |
| `flux_tol` | `1e-3` | residual ≤ tol/20 on every acceptance row |

R=0 is the hard case: at `conv_tol = 1e-8` its flux residual is `4.07e-3`
while R=16 sits at `1.50e-5`. Choosing off the large colony would have
repeated the 3b.2 failure. At `1e-10`, margins vs `flux_tol = 1e-3` are
24.7× (R=0), 390× (R=4), 3723× (R=12), 6725× (R=16) — all clear the 20×
rule. `1e-12` adds ~40% iterations for no gate value.

**Grid convergence (R=16, `conv_tol = 1e-8`):** marked volume stable to ~2%
across `dx ∈ {1.0, 0.5, 0.25}`; `N_centre` moves 0.23% from 1.0→0.5 and
3 ppm from 0.5→0.25. The field is converged; `dx = 1.0` is the working
grid.

Record the same pair in `JOB3C_STANDING.md` with the non-transferability
warning. The 3b.3 lock is untouched: `check_job3b.py` must still PASS on
the existing CSVs, byte-identical, after every 3c change.

## Gates — pre-registered thresholds

CSV directory `results/job3c_seed101/`. Seed 101 on every row.

### 1. Colony growth
`N` reaches **≥ 256** by `t_end`. Every division is **binary** — one mother
yields the mother plus exactly one daughter, `ΔN = +1` per dividing cell, each
at `ℓ_div`.

**No power-of-two requirement on logged `N`.** Under a radial Monod gradient
divisions necessarily desynchronize, so `N` will pass through 3, 5, 7, …
That is the *signal* Gate 4 exists to measure, not a defect. A power-of-two
lineage gate here (as in Jobs 2 and 3, which had a uniform bath) would fail for
exactly the physical reason 3c is built to detect. Gates 1 and 4 must not
contradict each other.

### 2. Packing holds at scale
Max pairwise `δ_cc ≤ 0.25 µm` (the Job 3 cap) at **every** logged row, and no
two cell centres closer than `w₀ − 0.25 = 0.75 µm`. This is the Job 3 gate
re-run two orders of magnitude larger in N.

### 3. Radial depletion — the headline
At the tick where `N` first reaches 256:

- **Centre depletion `ΔN_centre ≥ 1.55×10⁻² mM`** (absolute ENGINEERING
  threshold, pre-registered before the science run). Under
  `CONSERVATIVE_MASS` at R=16 / `dx=1.0` the acceptance field gives
  `ΔN_centre = 5.468×10⁻² mM` (3.53× margin). The earlier "10× the 3b.3
  value" framing was legacy-occupancy arithmetic and is **retired**; the
  numeric threshold itself is unchanged.
- `N_local` at the colony centre-of-mass is **strictly below** `N_local`
  averaged over cells in the outer 25% annulus.
- Field remains finite and non-negative; `n_clamp = 0`.

### 4. Radial growth-rate gradient — the Warren-type observable
Bin cells by distance from colony centre-of-mass into an **inner disc**
(`r < R/3`) and an **outer annulus** (`r > 2R/3`). At matched time-since-birth,

**mean interior elongation rate < mean edge elongation rate**, with the
separation exceeding the within-bin standard error.

`R` is defined as the **maximum pole distance from the colony centre-of-mass**
(frozen definition; not convex hull, not radius of gyration).

**Acceptance prediction (not a retune):** under TAKEN `K_S = 0.02 mM ≪ C_s`,
Monod stays near saturation across the realizable nutrient range in this
chamber. At R=16 the centre-to-edge Monod separation is **0.204%** — below
what elongation-rate binning can resolve against division-phase scatter.
Enlarging the chamber to 80–100 µm does **not** rescue Gate 4: reaching
`N ~ K_S` needs ~96% depletion (several penetration depths), outside the
feasible cost envelope. **Gate 4 stays frozen.** A miss is the pre-registered
`K_S ≪ C_s` finding. Do not rewrite this gate to the nutrient gradient
(that is Gate 3) and do not retune `K_S`.

### 4b. Morphology (new, 2026-08-23)

The colony must be a monolayer disc, not a filament and not a stack:

- **`offplane ≈ 0`** — no cell centre more than one radius off the mid-plane.
  Enforced by `b_z = 2r`; a nonzero value means the geometry is not doing its
  job.
- **`R/R_disc → ~1`** by N=256 (measured 0.99). A value growing past ~2 means
  the run has reverted to filament growth — check `DivisionMode`.
- **`d_centres_min ≈ w₀`**, not `2·ℓ₀`. Lateral Hertzian must actually engage.

### 5. Solver quality
`gs_converged = 1` on every row `t > 0`; `flux_rel_err ≤ 1e-3` (frozen
acceptance `flux_tol`), on **every** row (no warm-up exemption, per 3b.3).

### 6. Regression
`check_job2.py`, `check_job3.py`, `check_job3b.py` all still PASS, and the
3b.3 CSVs regenerate **byte-identical**.

## Frozen run parameters

| Symbol | Value | Class |
|---|---|---|
| chamber | `60 × 60 × 1 µm` (monolayer enforced) | ENGINEERING |
| division | `SYMMETRY_BROKEN`, seed 101 | ENGINEERING (mechanism TAKEN: Melke/Cho) |
| BC | `DISH_LATERAL` (4 lateral Dirichlet `C_s`, z faces flux-free) | ENGINEERING |
| occupancy | `CONSERVATIVE_MASS` | ENGINEERING |
| `ω` | `1.96` | ENGINEERING |
| `dx` | `1.0 µm` | ENGINEERING |
| `pad` | `1` voxel = `1.0 µm` rim | ENGINEERING |
| `conv_tol` | `1e-10` mM | ENGINEERING |
| `flux_tol` | `1e-3` | ENGINEERING |
| seeds | 1 founder at chamber centre | ENGINEERING |
| `t_end` | `23400 s` = **23,400 ticks** (`DT_S = 1.0 s`) = **9** `T_div` at `C_s` | ENGINEERING |

`t_end` carries deliberate margin: 8 doublings (1 → 256 cells) would take
`8 × T_div = 20760 s` under an unlimited bath, but growth **slows** as the
colony passes one penetration depth — that slowdown is the whole point of
3c. The 9th division time is headroom for it. **Do not extend `t_end`
after a near-miss on Gate 1**; a colony that cannot reach 256 in 9 `T_div`
is a finding about nutrient limitation, not a run that needs more time.

| log `dt` | `10 s` | ENGINEERING |
| seed | `101` | ENGINEERING |
| `λ_S`, `K_S`, `C_s`, `ℓ₀`, `ℓ_div`, `k_cc`, `D_C`, `D_A`, `Y`, `ρ₀` | Job 2/3/3b locked | TAKEN — do not retune |

## Non-goals — OFF in 3c

Quorum sensing, phage, chemotaxis, metabolism ODE, Brownian motion, `F_s`,
stochastic `P(ℓ)`, cell death, a second species, any reservoir readout, any
task score, ridge regression, voxel feature extraction.

Do not edit `BSimChemicalField.java`. Do not use
`BSimCapsuleBacterium.k_growth`. Do not start Job 5 (QS) from this file.

## Order of work

1. Freeze this protocol. ✓
2. Add `DISH_LATERAL` BC mode to `NutrientField` (additive; 3b path
   bit-identical, verified by diff + `check_job3b.py`). ✓
3. Run the `dx` acceptance study; derive and freeze the tolerance pair. ✓
4. Implement `Job3cSims` + `check_job3c.py` (incl. Gate 4b morphology). ✓
5. Science run; lock `JOB3C_STANDING.md`. ✓

**Status 2026-08-23: all six gates PASS.** See `JOB3C_STANDING.md`.

An earlier draft of this section predicted Gate 4 FAIL from a 0.2% Monod
separation on a hand-packed acceptance disc. The science run measured
**1.196%** separation and Gate 4 **PASS** as a sign test (34× combined SEM).
Prescribed configurations underestimate self-organised density; the `K_S ≪ C_s`
finding remains (effect is still modest) but is not a gate failure.

Stop at the first failing gate. Do not retune to make a gate pass.

# Job 6 protocol — carbon-starvation viability, standalone

**FROZEN BEFORE IMPLEMENTATION (2026-08-24).** This job fills the original
chassis mandate's missing death/viability layer without inventing a switch that
would alter the locked Job 3c dish.

Sources:

- Biselli, Schink & Gerland 2020, *Molecular Systems Biology* 16:e9478,
  DOI `10.15252/msb.20209478`
- Mandatory erratum, DOI `10.15252/msb.202010070`
- Schink, Biselli, Ammar & Gerland 2019, *Cell Systems* 9:64–73.e3,
  DOI `10.1016/j.cels.2019.06.003`

## Scope and boundary

Job 6 is a deterministic, standalone population module for *E. coli* K-12
NCM3722 after transfer to **zero carbon**. It uses days for starvation, hours
for pre-starvation growth rate, and fmol/CFU for the recycling mechanism.

It does **not** import BSim, `EcoliRodCell`, `NutrientField`, or Job 5. Death
remains OFF in Jobs 2/3/3b/3c. In particular, Job 3c's centre nutrient is
`0.276 mM = 13.8 K_S` and its Monod factor is about 0.93: that colony is not in
the zero-carbon regime measured by Biselli. Coupling this law to Job 3c would
require an uncited starvation-onset/interpolation rule and is **forbidden in
this job**.

## Published model

The central growth–death fit used for the feast–famine calculation is

```
gamma(mu) = 0.21 day^-1 exp[(1.0 h) mu]
```

where `mu` is the growth rate before starvation in `h^-1`. Viability obeys

```
dN/dt = -gamma(mu) N
N(t)/N(0) = exp[-gamma(mu)t].
```

Schink's supply/demand mechanism gives `gamma = beta/alpha`, where `beta` is
maintenance demand (`fmol CFU^-1 day^-1`) and `alpha` is recycled nutrient per
death (`fmol CFU^-1`). `alpha` is **not** Biselli's dimensionless recycling
yield `alpha_G`; substituting one for the other is dimensionally wrong.

The parameter-free mixed-live/UV-killed lag law is

```
T = (N_UV/N_viable) / gamma.
```

For feast time `T_plus` and starvation time `T_minus`,

```
f(mu) = mu T_plus - gamma(mu) T_minus
mu* = (1/a) ln[T_plus / (a gamma_0 T_minus)]
```

with all times expressed consistently and `a = 1.0 h`.

## Source-chain corrections carried explicitly

1. A surviving body-text fit prints the death prefactor with `h^-1`. The
   supplementary Table EV1 header and Fig. EV1 caption give `day^-1`, as do
   the model equation and maintenance units. **Job 6 uses `day^-1`; `h^-1`
   would make death 24× too fast.**
2. The original worked point (`T_plus=3 h`, `T_minus=3 d`) is inconsistent
   with its own equation. The erratum changes `T_minus` to **6 d**, yielding
   `mu*=0.86 h^-1` and `gamma=0.50 day^-1`.

## Frozen parameters

| Symbol | Value | Class |
|---|---|---|
| `gamma_0` | `0.21 day^-1` | TAKEN, Biselli model Eq. (2) |
| death slope `a` | `1.0 h` | TAKEN, central rounded fit used in Eq. (1) |
| WT glycerol `gamma` | `0.43 day^-1` | TAKEN, Schink |
| WT glycerol `beta` | `0.49 fmol CFU^-1 day^-1` | TAKEN, Schink |
| glycerol pulse `E_0` | `0.40 fmol CFU^-1` | TAKEN, Schink Fig. 3 |
| viability integration `dt` | `0.01 day` | ENGINEERING |
| viability horizon | `10 day` | ENGINEERING, covers EV1 survival traces |

## Gates — pre-registered

Output directory: `results/job6_starvation/`.

### 1. Units and experimental-range guard

The central law must fall within one reported SD for the WT glycerol
chemostat rows:

| `mu` (`h^-1`) | observed `gamma` (`day^-1`) |
|---|---|
| 0.10 | `0.24 ± 0.02` |
| 0.30 | `0.29 ± 0.03` |
| 0.50 | `0.33 ± 0.03` |
| 0.70 | `0.40 ± 0.04` |

The same numerical prefactor interpreted as `h^-1` must fail this gate by
more than 10×, serving as an explicit regression against the source typo.

### 2. Viability integration

RK4 integration of `dN/dt=-gamma N` must match `exp(-gamma t)` over 10 days
with maximum absolute error `< 1e-9`, remain finite/non-negative, and be
monotone non-increasing.

### 3. Parameter-free recycling lag

Using Schink's WT glycerol `gamma=0.43 day^-1`:

- 50:50 UV-killed:viable predicts `2.33 d`, within `0.10 d` of `2.3 d`;
- 30:70 predicts `1.00 d`, within `0.25 d` of `1.2 d`.

### 4. Maintenance-lag mechanism

`T=E_0/beta` with `E_0=0.40 fmol CFU^-1` and
`beta=0.49 fmol CFU^-1 day^-1` must lie within 2% of the reported `0.825 d`.

### 5. Feast–famine optimum and erratum guard

For at least three `T_plus/T_minus` ratios, a numerical argmax on a
`1e-4 h^-1` grid must agree with the closed form to `< 2e-4 h^-1`.
For the erratum point (`3 h`, `6 d`), require:

- `mu*` in `[0.84, 0.89] h^-1`;
- `gamma(mu*)` in `[0.48, 0.52] day^-1`.

The pre-erratum `3 h`, `3 d` point must **not** satisfy those intervals.

### 6. Mechanism dimensional consistency

`alpha = beta/gamma = 1.14 fmol CFU^-1` for Schink WT glycerol. The checker
must reject substitution of dimensionless `alpha_G=0.18` into
`gamma=beta/alpha`; it would predict `2.72 day^-1` instead of `0.43 day^-1`.

### 7. Chassis regression

`check_job2.py`, `check_job3.py`, `check_job3b.py`, and `check_job3c.py` all
still exit 0. Their locked CSVs remain byte-identical.

## Claims and non-goals

**MAY claim:** reproduction of the published starvation viability law, the
Schink supply/demand lag predictions, and the Biselli feast–famine optimum.

**MAY NOT claim:** death in a nutrient-limited but non-starved colony, single
cell death timing, lysis mechanics, or a Job 3c mortality prediction.

OFF: stochastic death events, lysis/recycling PDE, Monod coupling, motility,
QS, metabolism, phage, PocketDish, reservoir/task scores.

## Order of work

1. Freeze this protocol. ✓
2. Implement `StarvationViability.java` and standalone `Job6Sims.java`.
3. Implement `check_job6.py`; stop on the first failed gate.
4. Run chassis regressions and write `JOB6_STANDING.md`.

# Job 5b protocol — spatial AHL transport precondition, mm-scale lane

**FROZEN BEFORE IMPLEMENTATION (2026-08-24).** This job validates the spatial
transport layer required before any QS–chassis coupling. It does not claim to
repair Job 5's unexplained upper saddle-node and does not put a spatial field
into the 60 µm Job 3c dish.

Primary source: Dilanji et al. 2012, *JACS* 134:5618–5626,
DOI `10.1021/ja211593q`, including SI. Independent parameter checks:

- Leaman & Behkam 2018, DOI `10.1021/acssynbio.7b00406`;
- Grant et al. 2016 Appendix S8, DOI `10.15252/msb.20156590`.

## Why this is a separate geometry

Dilanji's measured lane is 32 mm long; reporter activation spreads about 1 cm
in about 10 h. With their `D=1.98 mm²/h = 550 µm²/s`:

```
10 mm: L²/(4D) = 12.63 h
60 µm: L²/(4D) = 1.64 s
```

Thus AHL is effectively well-mixed in the locked Job 3c dish relative to
Weber's 25 min autoinduction and hour-scale growth. A spatial field there
would pass a “matches well-mixed” gate by construction. Job 5b instead
validates Dilanji's mass-conserving transport problem at its published
millimetre scale.

## Frozen transport model

```
partial C / partial t = D partial²C / partial x²
C(x,0) = C_inf L/nu    for 0 <= x < nu
         0             otherwise
partial_x C(0,t) = partial_x C(L,t) = 0
```

Parameters:

| Symbol | Value | Class |
|---|---|---|
| `L` | `32 mm` | TAKEN, Dilanji |
| `nu` | `2 mm` | TAKEN, Dilanji |
| `D` | `1.98 mm²/h` (`550 µm²/s`) | TAKEN, Dilanji |
| `C_inf` | `4 nM` | TAKEN, Dilanji sensor lane |
| reporter half-activation `a` | `1.5 nM` | TAKEN, Dilanji |
| working `dx` | `0.10 mm` | ENGINEERING |
| explicit finite-volume CFL | `0.40` | ENGINEERING |
| horizon | `30 h` | ENGINEERING |

This job implements the AHL PDE only. Dilanji's GFP model also requires the
experiment-specific measured initial OD `n_0`; it is not a fixed Table 1
number. Reporter kinetics are therefore retained as a validated future layer,
not silently supplied with an invented `n_0`.

## Gates — pre-registered

Output directory: `results/job5b_spatial_ahl/`.

### 1. Unit and literature-range gate

`1.98 mm²/h` must convert to `550 µm²/s` and lie within 15% of Leaman's
`490 µm²/s`. It must also lie inside Grant's two agar estimates
`379–991 µm²/s`.

### 2. Mass, positivity, and boundary conservation

At `dx={0.20,0.10,0.05} mm`, total AHL mass must remain constant to relative
error `<1e-11`; all concentrations must remain finite and non-negative.
The conserved domain mean must equal `C_inf=4 nM`.

### 3. Grid convergence

At `x=10 mm`, the first time `C>=a=1.5 nM` must move by `<0.05 h` between
`dx=0.10` and `0.05 mm`. The 0.10 mm grid is the working point.

### 4. Diffusion identity

The numerical concentration profile at 10 h on the 0.10 mm grid must agree
with a 0.05 mm reference profile (restricted by averaging pairs) with
normalized RMSE `<0.01`.

### 5. Geometry-separation guard

The analytic diffusion scale must put 10 mm in `[10,15] h` and 60 µm below
`2 s`. This gate prevents accidental reuse of the mm-scale spatial claim in
the 60 µm Job 3c chamber.

### 6. Regression

Jobs 2/3/3b/3c and standalone Job 6 remain PASS. Job 5's partial verdict and
frozen failed upper-bound gates remain unchanged.

## Claim boundary

**MAY claim:** conservative, grid-converged 1-D 3OC6-HSL transport at the
Dilanji lane scale, with a diffusivity independently corroborated by Leaman
and bracketed by Grant.

**MAY NOT claim:** reproduction of Dilanji GFP activation without its measured
`n_0`; spatial variation in the 60 µm Job 3c dish; repair/full reproduction of
Weber Fig. 3; a coupled communicating colony.

## Order of work

1. Freeze this protocol. ✓
2. Implement standalone `SpatialAhlLane.java` and `Job5bSims.java`.
3. Implement `check_job5b.py`; run refinement and regressions.
4. Write `JOB5B_STANDING.md`.

# Job 7 protocol — unbiased *E. coli* run-and-tumble motility

**FROZEN BEFORE IMPLEMENTATION (2026-08-24).** Job 7 fills basal motility
without claiming chemotaxis. It is a standalone, seeded ensemble validation;
the locked colony remains non-motile.

Primary quantitative source: Kurzthaler et al. 2024,
*Physical Review Letters* 132:038302,
DOI `10.1103/PhysRevLett.132.038302`. Mechanism corroboration:
Saragosti, Silberzan & Buguin 2012, DOI
`10.1371/journal.pone.0035412`, and Berg & Brown 1972,
DOI `10.1038/239500a0`.

## Model and units

In a homogeneous medium, a bacterium alternates:

- RUN: straight swimming at sampled speed `v`;
- TUMBLE: zero translational speed, followed by an isotropically sampled new
  3-D direction.

Run and tumble durations are exponential. Speeds are sampled from the
measured normal population distribution, rejecting the negligible negative
tail. This is the renewal model used for the published long-time identity:

```
D_eff = (mean(v)^2 + sigma_v^2) tau_R^2
        / [3 (tau_R + tau_T)].
```

Parameters:

| Symbol | Value | Class |
|---|---|---|
| mean speed | `16.0 µm/s` | TAKEN |
| speed SD | `5.78 µm/s` | TAKEN |
| mean run `tau_R` | `2.39 s` | TAKEN |
| mean tumble `tau_T` | `0.38 s` | TAKEN |
| measured `D_eff` | `185 ± 7 µm²/s` | TAKEN observable |
| renewal prediction | about `199 ± 11 µm²/s` | DERIVED/published |
| trajectories | `20,000` | ENGINEERING |
| horizon | `1000 s` | ENGINEERING |
| seed | `101` | ENGINEERING, reproducibility |

The small fitted passive diffusivity (`0.24 µm²/s`) is OFF: it is three orders
below active `D_eff` and is not part of the stated renewal identity.

## Gates — pre-registered

Output directory: `results/job7_motility/`.

### 1. Parameter/distribution recovery

Across generated events, sampled mean run time, tumble time, speed mean, and
speed SD must each be within 2% of their frozen source values.

### 2. Run fraction and run length

Require `p_R=tau_R/(tau_R+tau_T)` within `0.01` of the measured `0.86`, and
mean run length within 2% of `16*2.39=38.24 µm`.

### 3. Effective diffusion

At 1000 s, `MSD/(6t)` must be within 7% of the renewal identity and inside
`[160,220] µm²/s`, a window containing both the measured `185±7` and the
published `198±11` prediction.

### 4. Ballistic-to-diffusive crossover

The log-slope of ensemble MSD from 0.1 to 0.5 s must be `>1.7`; from 100 to
1000 s it must be in `[0.90,1.10]`.

### 5. Isotropy and reproducibility

At 1000 s, each coordinate mean displacement must be less than 2% of RMS
displacement. Re-running with seed 101 must produce byte-identical CSVs.

### 6. Regression

Jobs 2/3/3b/3c/5b/6 remain unchanged; Job 5 remains the same partial port.

## Boundary

**MAY claim:** reproducible basal 3-D run-and-tumble motility with published
single-cell statistics and long-time dispersion.

**MAY NOT claim:** chemotaxis, Tar/Che signalling, motility inside the crowded
Job 3c monolayer, swimming through walls/cells, or nutrient-directed drift.
Those require a separate coupling job and a geometry compatible with swimming.

## Order of work

1. Freeze this protocol. ✓
2. Implement `RunTumbleMotility.java` and `Job7Sims.java`.
3. Implement `check_job7.py`, including deterministic rerun.
4. Run regressions and write `JOB7_STANDING.md`.
A Gate 4 miss under frozen TAKEN Monod parameters is a finding, not a
reason to enlarge the chamber or rewrite the observable.

# Job 5 protocol — Weber & Buceta LuxI/LuxR quorum sensing, well-mixed

**FROZEN BEFORE IMPLEMENTATION.** Every threshold below is pre-registered.
Do not retune after seeing output. If a gate fails, the failure is the finding.

Source: Weber M, Buceta J (2013) "Dynamics of the quorum sensing switch:
stochastic and non-stationary effects." *BMC Systems Biology* 7:6.
DOI **10.1186/1752-0509-7-6**. Open access. The only **BUILD** row in
`CITATION_DOSSIER.md`.

## Two declarations that are first-class, not footnotes

**1. SPECIES CROSSING — ENGINEERING.** Weber models *Vibrio fischeri*. Our
chassis freeze (`CHAT_BOOTSTRAP.md`) is one species, non-motile *E. coli*. The
lux circuit is a standard part routinely expressed in *E. coli*, so the port is
defensible — but it is a **declared engineering decision, not a TAKEN result**.
Every figure produced by Job 5 must carry this label. Leaman & Behkam 2018
(`10.1021/acssynbio.7b00406`) is the *E. coli*-native alternative and remains
CITE-ONLY because its observable is entangled with chemotaxis.

**2. `τ = 45 min` IS WEBER'S, NOT OURS.** Weber Table 1 gives a doubling time of
45 min in RM/succinate at 30 °C. Our chassis `T_div = 43.25 min` at `C_s` is a
different organism in a different medium and is **not** used here. The 4%
coincidence is a coincidence. **Do not cite it as corroboration and do not
substitute one for the other.** Job 5 growth dilution is `ln(2)/τ` with Weber's
`τ` and nothing else.

## Scope — well-mixed only

Job 5 is **step 1 of the two-step rule** in `CITATION_DOSSIER.md`: port Weber
exactly as published, gate on Weber's own observable, freeze. Spatial AHL is a
separate later job with `D_AHL` from Leaman Table 2 (`490 µm²/s`, cited to
Stewart 2003).

**Weber's `D = 10 min⁻¹` is a first-order transmembrane transport rate, not a
diffusivity.** Wrong dimensions, wrong quantity. It must never be fed into
`NutrientField` or any PDE.

Job 5 does **not** touch the Job 2/3/3b/3c chassis. No `EcoliRodCell`, no
`ValdezHertzian`, no `NutrientField`, no Monod growth, no `λ_S`/`K_S`. It is a
standalone well-mixed ODE module. Coupling is a later question.

## Units

Job 5 runs in **Weber's native units: nM and min**. The chassis runs in mM, µm,
s. There is no conversion anywhere in Job 5 because there is no coupling. Any
future coupling must convert explicitly at the boundary and declare it.

## Model — Weber Eq. (1), eleven species

Species: `luxR`, `luxI::gfp`, `mRNA_luxR`, `mRNA_luxI::gfp`, autoinducer `A`,
`luxR·A`, `(luxR·A)₂`, `DNA`, `DNA·(luxR·A)₂`, external `A_ext`.

Reactions, verbatim from Weber Eq. (1): basal transcription at `α_R k_R` and
`α_I k_I`; activated transcription at `k_R`, `k_I` from `DNA·(luxR·A)₂`;
translation at `p_R`, `p_I`; autoinducer synthesis `luxI::gfp → A + luxI::gfp`
at `k_A`; binding `luxR + A ⇌ luxR·A`; dimerisation `2(luxR·A) ⇌ (luxR·A)₂`;
promoter binding `(luxR·A)₂ + DNA ⇌ DNA·(luxR·A)₂`; membrane transport
`A ⇌ A_ext` at `D`/`rD`; degradation of every species; DNA duplication at
`ln(2)/τ`.

Growth dilution `−c_X ln(2)/τ` on all species **except `A_ext`**. DNA
duplication adds `+ (ln2/τ)(c_DNA + c_DNA·(luxR·A)₂)` so total DNA is conserved.

**Constructs:** `lux01` has `k_A = 0` (no LuxI, exogenous induction only);
`lux02` has `k_A = 0.04 min⁻¹` (autoinduction restored).

## Frozen parameters — Weber Table 1, TAKEN unless noted

| Symbol | Value | Weber's own label |
|---|---|---|
| `K_d1` LuxR–A dissociation | 100 nM | Urbanowski 2004 |
| `k⁻_1` LuxR–A unbinding | 10 min⁻¹ | estimated |
| `K_d2` dimerisation | 20 nM | fitted |
| `K⁻_2` dimer dissociation | 1 min⁻¹ | estimated |
| `k_A` A synthesis by LuxI | 0.04 min⁻¹ (`lux02`); 0 (`lux01`) | fitted |
| `K_dlux` complex → promoter | 200 nM | fitted |
| `k⁻_lux` promoter dissociation | 10 min⁻¹ | estimated |
| `b` burst size | 20 | Cai/Friedman/Xie 2006 |
| `k_R` luxR transcription | 200/b = 10 min⁻¹ | fitted |
| `k_I` luxI transcription | 50/b = 2.5 min⁻¹ | fitted |
| `p_R`, `p_I` translation | `b·d_mR` = `b·d_mI` = 6.94 min⁻¹ | derived |
| `α_R`, `α_I` basal ratio | 0.001, 0.01 | fitted |
| `d_A` autoinducer degradation | 0.001 min⁻¹ | Kaufmann 2005 |
| `d_C2`, `d_C`, `d_R` | 0.002 min⁻¹ each | estimated |
| `d_I` LuxI degradation | 0.01 min⁻¹ | estimated |
| `d_mR`, `d_mI` mRNA degradation | 0.347 min⁻¹ each | Roberts 2006 |
| `D` membrane transport | 10 min⁻¹ | Kaplan & Greenberg 1985 |
| `τ` doubling time | **45 min** | Williams 2008 |
| `V_0` cell volume | 1.5 µm³ | Trueba & Koppes 1998 |
| `V_tot` culture volume | 2×10⁻⁴ µL | — |
| `c_N` cell density | 5×10⁸ cells/mL, `N = 100` | estimated from OD 0.5 |

Derived: `V_ext = V_tot − N·V_0`, `r = V_c,tot/V_ext`.

`k_R = 200/b` and `k_I = 50/b` are **TAKEN (paper-fitted)** — verbatim Weber
Table 1, fitted by the authors to Williams 2008 response curves. They are not
inherited guesses; see the dossier correction.

## Dilution protocol

Cell density is held constant, as in the experiments Weber models. Continuous
efflux removes culture medium at the rate that compensates growth, so `V_c,tot`
stays constant; a matching influx supplies fresh medium at exogenous
concentration `c_A*`. Weber's reaction form:

`A_ext ⇌ ∅` at rate `γ` out and `γ c_A* V_tot` in, with `γ = ln(2)/τ`.

`c_A*` is **the control parameter** — the swept quantity in every gate below.

## Gates — pre-registered

Output directory `results/job5_seed101/`. Deterministic ODE; seed 101 applies
only to the optional Gillespie layer.

### 1. lux01 hysteresis window (Weber Fig. 3A)

Sweep `c_A*` up 0 → 100 nM then down, by continuation (each point starts from
the previous steady state), integrating 100 h at each point. The bistable window
is where up- and down-sweeps differ.

- lower bound in **[1.0, 3.0] nM** (Weber: 2 nM)
- upper bound in **[12.0, 18.0] nM** (Weber: 15 nM)

### 2. lux02 hysteresis window (Weber Fig. 3B)

Same sweep with `k_A = 0.04 min⁻¹`.

- lower bound **≤ 0.5 nM** — LuxI autoinduction must pull it to zero. This is
  the discriminating result between the two constructs.
- upper bound in **[12.0, 18.0] nM** (Weber: 15 nM)

### 3. GFP dynamic range

`GFP_high / GFP_low ≥ 50` at steady state. Weber states "two orders of
magnitude"; the gate allows half an order of slack below that.

### 4. Settling time

From an uninduced start at `c_A* = 100 nM`, GFP reaches **95% of its 100 h
value within 360 min** (Weber: "shorter than 6 hours").

### 5. Integrator correctness — DNA conservation

Weber states DNA duplication "compensates exactly for the cell growth dilution
such that `c_DNA,tot = c_DNA + c_DNA·(luxR·A)₂` is kept constant."

- relative drift of `c_DNA,tot` over a full 100 h run **< 1×10⁻⁶**

This is an exact invariant of the published model and therefore an independent
check on the integrator — the analogue of the `ρ=0 → N ≡ C_s` identity test
that caught discretisation errors in Job 3c.

### 6. Non-negativity and finiteness

All eleven species remain finite and `≥ 0` at every logged step, on every sweep
point, for both constructs.

### 7. Integrator numerical convergence

The integrator's exposed ENGINEERING parameter — **fixed step size, or
adaptive error tolerance, whichever the method uses** — must be measured,
not assumed. Halving that parameter (step → step/2, or tol → tol/2) must
change both hysteresis bounds by **< 0.1 nM**.

**AMENDED 2026-08-23 (before any sweep numbers).** Weber's ON state is stiff
under Table 1 parameters (promoter/binding modes ~70–280 /min at
`c_A* = 100 nM`); fixed-step RK4 at a viable margin is not practical for the
continuation sweep. Adaptive step-doubling RK4 with error control is the
allowed fix. Quasi-steady-state reduction of the fast binding equilibria is
**forbidden** — that would change Weber's model; Job 5 is a port.

Intent unchanged from the original Gate 7: an unverified numerical parameter
must not silently set the answer (Job 3c lesson).

### 8. Chassis regression

`check_job2.py`, `check_job3.py`, `check_job3b.py`, `check_job3c.py` all still
PASS, and their CSVs are byte-identical. Job 5 must not touch the chassis.

## What Job 5 may and may not claim

**MAY:** "reproduces the published model output of Weber & Buceta 2013 Fig. 3
from their Table 1 parameters."

**MAY NOT:** "reproduces experimental lux hysteresis", or any claim of agreement
with Williams 2008. See the RESOLVED section of `CITATION_DOSSIER.md`: Weber
states Williams' window as 0–15 nM, a units-recovery pass read it as 0–50 nM,
and the conflict is unresolved because we do not have the Williams PDF. It does
not affect this gate, which tests our implementation of Weber, not Weber's fit
to data.

## Non-goals — OFF in Job 5

Spatial AHL, Weber's `D` as a diffusivity, Melke's dimensionless numbers, any
Hill threshold not in Weber Table 1, LasR/3OC12HSL crosstalk (Grant 2016 is
CITE-ONLY), chemotaxis, phage, nutrient coupling, Monod growth, reservoir
readout, ridge regression, task score.

Do not couple Job 5 to the chassis in this job. Do not open Job 3d.

## Order of work

1. Freeze this protocol. ← *you are here*
2. `Job5Sims.java` — standalone well-mixed integrator, Weber units, both
   constructs. No chassis imports.
3. Integrator convergence study (Gate 7) → freeze the exposed numerical
   parameter (step or adaptive error tol) as ENGINEERING.
4. `check_job5.py` against the frozen table.
5. Sweep run → `JOB5_STANDING.md`.
6. *Optional, separate:* Gillespie layer for Weber Figs. 4/8/9 burst-size
   results. Not required for the lock.

Stop at the first failing gate. Do not retune to make a gate pass.
