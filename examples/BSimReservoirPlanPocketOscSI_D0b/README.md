# PocketOsc-SI D0b — predeclared AHL-kick bulk oracle

D0-reopen of the Danino 2010 SI bulk delay-DDE. Not D1.

The first protocol
[`examples/BSimReservoirPlanPocketOscSI_D0/`](../BSimReservoirPlanPocketOscSI_D0/)
is **FAIL** and is left in place. This package is a **new** predeclared
IVP and \(\mu\) line. Diagnosis scripts in D0 are cited, not used as a
pass artifact.

`examples/BSimReservoirPlanPocketOscSI/` is the **failed prior twin**.
Do not patch it.

## What D0b is

- Same SI Hill production as D0:
  \(P=(\delta+\alpha H_\tau^2)/(1+k_1 H_\tau^2)\),
  \(H_\tau=H_i(t-\tau)\).
- Same TAKEN table, \(d=0.5\), bulk \(D_1=0\), SI-scaled time.
- Primary IC: `AHL_KICK_005` (\(A=I=0\), \(H_i=H_e=0.05\), history 0.05).
- Identity \(\mu=0.32,0.36,\ldots,0.60\). Extras include the failed D0
  line \(\mu=1\)–\(2\).
- Period from LuxI \(I\) with the D0 peak checker (unchanged).

## What D0b is not

Java circuits, packed rods, spatial AHL, channel flow, ACs, NARMA,
GFP maps, HybridDish, D1g 4-ODE, \(\mu_{\mathrm{bus}}\), `QS_KMLA`,
a retune of \(\alpha,\tau,\gamma_\ast,k_1,d\), or a rewrite of
[`D0_BULK_ORACLE_STANDING.md`](../PocketDish/D0_BULK_ORACLE_STANDING.md).

## Run

Protocol must exist with `frozen_before_traces=true` **before** the
identity grid. It does.

```
python examples/BSimReservoirPlanPocketOscSI_D0b/run_d0b.py
python examples/BSimReservoirPlanPocketOscSI_D0b/check_d0b.py
```

Standing record:
[`examples/PocketDish/D0B_BULK_ORACLE_STANDING.md`](../PocketDish/D0B_BULK_ORACLE_STANDING.md).

Equations: [`EQUATIONS.md`](EQUATIONS.md).
Frozen protocol: [`PROTOCOL.md`](PROTOCOL.md).
