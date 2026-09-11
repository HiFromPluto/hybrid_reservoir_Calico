# Post-gate E0 preflight protocol — frozen before calculation

Frozen: 2026-08-17. This package is analysis-only. It does not authorize
Java, BSim, Waveform2, Lorenz, IPC, AC implementation, biomarker BSim,
independent-input NARMA, or Stage99. Historical gate decisions remain
unchanged.

## 1. Decision scope

The package will:

1. record the completed zero-simulation disposition;
2. inventory implemented HybridDish parameters and provenance;
3. validate a reduced AHL transport/receiver/reporter model against
   existing claim-dish traces;
4. run a target-free physical screen only if the calculations are
   numerically sound;
5. generate an E0.3 living-BSim shortlist only if transport validation
   passes every frozen gate.

No NRMSE, AUC, NARMA target, waveform label, Lorenz target, or patient
target may select a transport condition.

## 2. Frozen implemented reference

- Domain: `1000 x 500 x 10` micrometres.
- Chemical grid: `50 x 25 x 1`; the reduced model is two-dimensional
  because the implemented z grid has one box. `dx=dy=20 um`; depth is
  `10 um`; cell volume is `4000 um^3`.
- AHL: `D=159 um^2/s`, first-order decay `k=0.0033 1/s`.
- Source: `(500,250,5) um`, `1.28e8 molecules/s` at command 1 for the
  first `75 s` of each `300 s` window.
- Command map: `J(t)=1.28e8*u_window molecules/s`; warmup command is
  `0.5`.
- Warmup: `18000 s`; screen sampling: every `20 s`, including the
  window start and a final sample at `299.95 s` when comparison with
  stored BSim is required.
- Receiver/reporter:
  `H(C)=C^2/(1.6^2+C^2)`,
  `dR/dt=(H-R)/15`,
  `dL/dt=(R-L)/1500`.
- Conversion: `1 uM = 602.2 molecules/um^3`.

The source-cell forcing in concentration units is
`J/(602.2*4000) uM/s`. Total mass is always reported in molecules.

## 3. Reduced-model discretization and accounting

The screen uses a cell-centred finite-volume method on the implemented
chemical grid. Diffusive face fluxes use centred differences. Positive
x advection uses a conservative first-order upwind flux. The method of
lines is integrated with SciPy `solve_ivp` using DOP853, relative
tolerance `1e-7`, absolute tolerance `1e-10`, and maximum step `2 s`.
Source discontinuities are hard segment boundaries; a solver step never
straddles a pulse edge or window edge.

The state includes spatial `C`, `R`, and `L`. Decay and boundary loss
are integrated as explicit budget states. Every result reports initial
mass, injected mass, final mass, decay loss, boundary loss, and

`residual = initial + injected - final - decay_loss - boundary_loss`.

This identity is evaluated on two windows. The 1% residual ceiling is
unchanged on both.

1. **NARMA/screen window.** `remaining(t0)` is leftover mass after
   warmup. Injected mass, decay, and boundary loss are counted only on
   `(t0, t1]`. Budget counters are zeroed at `t0`. This is the frozen
   validation mass gate and every screen-row budget.
2. **Full horizon.** `remaining(t0)=0` before warmup. Injected mass,
   decay, and boundary loss include warmup and the subsequent drive.

Mixing leftover warmup mass with warmup-plus-drive injection, while
omitting warmup decay, is an accounting error. It is not a reason to
change D, decay, source strength, or boundaries.

An unresolved absolute residual above 1% of
`initial + injected` is a numerical failure. Concentrations below
`-1e-9 uM` are also a numerical failure.

Time convergence is checked on the claim condition and single-pulse
reset with a tighter DOP853 run (`rtol=1e-9`, `atol=1e-12`,
`max_step=0.5 s`). Peak/mean C, H, R and L, integrated mass terms, and
the 37/14/5/1% reset times must differ by no more than 1% (reset times
also pass when they differ by no more than 2 s). Failure blocks
validation and shortlisting.

This is a screening surrogate, not BSim evidence. It does not represent
cells, motility, growth, population equilibrium, clamp deaths, acid
mortality, or sub-voxel cell trajectories.

## 4. Boundary definitions

Only these mathematically explicit variants are allowed:

- `NO_FLUX`: zero normal diffusive flux on all faces and zero velocity.
  This is the active HybridDish chemical boundary.
- `OUTFLOW`: zero-concentration advective inflow at x=0, conservative
  upwind advective outflow at x=1000, and zero diffusive flux on all
  faces. This is used for positive x flow.
- `ABSORBING`: zero exterior concentration for diffusive flux on every
  wall, plus the same advective inlet/outlet rule when flow is positive.
- `ROBIN_LEAKY`: outward diffusive flux `h*C` on every wall, with
  `h=D_reference/Lx=0.159 um/s`; this coefficient is a derived,
  hypothetical sensitivity value, not calibration. Positive flow also
  uses the `OUTFLOW` advective rule.

`NO_FLUX` is not combined with wall-normal positive flow. Bacterial
x-boundary removal and y reflection are particle rules and are never
counted as chemical boundary loss.

## 5. Frozen validation

Validation uses every available driven claim-dish NARMA voxel trace:
seeds `111`, `222`, and `333` under
`examples/BSimReservoirPlanNarma10b/results/`. Their AHL columns are
checked for identity; identical traces remain three provenance records
but only one independent deterministic field realization.

The implemented source, D, decay, grid, no-flux boundary, warmup, pulse,
and input sequence are used without fitting. Model values are sampled
at the exact 20x10 readout centres with the same piecewise-constant
field lookup as `VoxelAnalyzer`.

Frozen gates, evaluated before inspection:

- correlation of sampled domain-mean AHL: at least `0.98`;
- domain-mean AHL NRMSE: at most `0.10`, where NRMSE is RMSE divided by
  the observed peak-to-peak range;
- pulse-integrated domain-mean AHL relative error: at most `10%`;
- mean-R relative error: at most `10%` for every seed. Predicted R is
  weighted by each stored sample's `Den_*`; observed R is
  `sum(Den_i*Receiver_R_i)/sum(Den_i)`. Empty voxels do not become
  chemical sinks;
- unresolved mass-balance error: at most `1%`;
- numerical convergence check: pass.

Source-near `(525,275)`, middle `(775,275)`, and far `(975,275) um`
voxel traces are compared descriptively by correlation and NRMSE; they
are not additional fitted gates.

All gates must pass. Otherwise:

`TRANSPORT_MODEL_STATUS: NOT_VALIDATED`

and ranking/shortlisting stops. Raw sensitivity calculations may remain
only as `NON_DECISION_GRADE`. If all pass:

`TRANSPORT_MODEL_STATUS: VALIDATED_FOR_SCREENING`.

## 6. Frozen standardized drives

All drives use command range `[0,1]`, the frozen `75 s` pulse, and no
task labels:

1. `SINGLE_PULSE`: one command-1 window followed by eleven zero windows.
2. `STEP_ON_OFF`: two zero, four command-1, then six zero windows.
3. `ALTERNATING_BINARY`: twelve windows `1,0,1,0,...`.
4. `FIXED_RANDOM`: 64 amplitudes generated before screening by Python
   `random.Random(20260817)`, rounded to six decimal places.

Canonical random-drive serialization is one six-decimal amplitude plus
newline. Its SHA-256 is
`5919f8fe318308163b1d798967fd8f89e99ca5eaa268af0b8a62e30e0b67307f`.
The exact sequence is embedded in `run_transport_screen.py` and its hash
is checked before any simulation.

## 7. Frozen screen conditions

Primary source positions (all at z=5 um):

- `CENTER=(500,250)`;
- `UPSTREAM_CENTER=(200,250)`;
- `DOWNSTREAM_CENTER=(800,250)`;
- `UPSTREAM_OFFAXIS=(200,150)`;
- `CROSSSTREAM_OFFAXIS=(500,350)`.

Every point is at least 100 um (five chemical cells) from a wall. These
are positions of one AHL source, not a multi-AC layout.

Positive x flow levels are fixed before results:

`0, 0.25, 0.5, 0.666667, 1.0, 2.0, 3.333333, 8.0 um/s`.

They correspond to domain transit times infinity, 4000, 2000, 1500,
1000, 500, 300, and 125 s. Thus the list brackets the 300 s window,
the AHL decay time `1/k=303.03 s`, and `tau_L=1500 s`. `8 um/s` is the
known washout anchor and is ineligible for promotion.

The primary screen is the complete five-position by eight-flow table:
zero flow uses `NO_FLUX`; positive flow uses `OUTFLOW`. All use the
implemented D, decay, source strength, and identical commanded payload.

Boundary sensitivity adds `ABSORBING` and `ROBIN_LEAKY` for CENTER at
flows `0` and `0.666667 um/s`.

Because no verified HybridDish uncertainty range was found before this
freeze, four explicitly hypothetical, one-factor CENTER/zero-flow arms
are added: `D=80`, `D=240 um^2/s`, `k=0.00165`, and `k=0.0066 1/s`.
These are design envelopes, not biological ranges and cannot be called
calibration.

## 8. Frozen physical metrics

For every condition and drive report:

- initial, injected, remaining, decay-loss, and boundary-loss mass;
- mass-balance residual;
- peak and mean AHL;
- source-relative near, middle, and far probe traces;
- plume coverage, defined as the fraction of voxel-samples with
  `H(C)>=0.10`;
- wall/interior concentration ratio, using the one-cell outer ring;
- mean and variance of H, R, and L;
- weak fraction `H<0.05`;
- saturated fraction `H>0.90`;
- lag-1 temporal autocorrelation of window-end mean C, R, and L.

For `FIXED_RANDOM`, the linear memory diagnostic uses window-end
features `[meanC,meanH,meanR,meanL,nearC,middleC,farC]`. For each delay
1..10, ordinary least squares with an intercept is fit on windows
16..47 and tested on windows 48..63 after train-only standardization.
Report every held-out R2 and `sum(max(0,R2_delay))`. No regularization,
lambda selection, target task, or use in condition selection is allowed.

Reset is a separate command-1 pulse followed by 40 zero-input windows.
At every zero-window end report L2-state residuals for C, R, and L,
normalized to each state's maximum norm after source shutoff. For each
state, the 37%, 14%, 5%, and 1% time is the first downward crossing
after that post-shutoff maximum; missing crossings are reported, not
extrapolated.

## 9. Frozen exclusions and shortlist rule

A condition is excluded from promotion for any of:

- numerical/convergence failure or mass residual above 1%;
- receiver effectively dead: random-drive mean R `<0.02` or weak
  fraction `>0.90`;
- receiver nearly always saturated: random-drive mean H `>0.90` or
  saturated fraction `>0.90`;
- severe wall hotspot: wall/interior AHL ratio `>3` and more than half
  of field mass lies in the outer ring;
- reset incompatible with the proposed sequential protocol: any of C,
  R, or L remains above 1% after 30 zero windows;
- flow `8 um/s`, which is retained only as the known dead-flow anchor.

Population clamp and acid mortality are neither simulated nor optimized.
They require later living-model robustness and calibration.

Only after validation passes, E0.3 may contain at most 12 conditions.
The shortlist must include the frozen claim condition, retain the
8-um/s result as a historical anchor rather than requiring a new run,
and use eligible primary-table conditions only. It will preserve a
factorial block of at least three positions by at least three low-flow
levels when nine such cells survive. Remaining slots, if any, retain
weak/intermediate/strong transport regimes by tertiles of random-drive
mean R. Deterministic tie order is position order above, then ascending
flow. No task score enters this rule. Matched-occupancy arms, if later
requested, are separate and are not part of this primary payload-matched
screen.

Every promoted living condition requires driven living cells, silent
cells, Brownian density, field-only, unmasked kinetic surrogate,
occupancy-masked surrogate, and direct-input baseline controls.

