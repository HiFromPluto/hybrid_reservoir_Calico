# PocketOsc-SI D0 — published Danino 2010 bulk delay-DDE oracle

Independent transcription of the Danino 2010 SI bulk delay-DDE, plus a
Figure 4 reproduction protocol. This is Gate D0 of AC–DaninoPocket.

`examples/BSimReservoirPlanPocketOscSI/` is the **failed prior twin**.
Do not patch it. Do not copy `DANINO_SI_MODEL.md` as the oracle.

## What D0 is

- SI Hill production
  \(P=(\delta+\alpha H_\tau^2)/(1+k_1 H_\tau^2)\),
  \(H_\tau=H_i(t-\tau)\).
- Bulk \(D_1=0\), \(d=0.5\), TAKEN table, SI-scaled time.
- Predeclared IC/history ensemble, because the SI does not state bulk ICs.
- Period from LuxI \(I\), aimed at Fig. 4b/c, not at the PocketOsc-SI
  occupancy screen.

## What D0 is not

Java circuits, packed rods, spatial AHL, channel flow, ACs, NARMA,
GFP maps, HybridDish, D1g 4-ODE, \(\mu_{\mathrm{bus}}\), `QS_KMLA`,
or a retune of \(\alpha,\tau,\gamma_\ast,k_1,d\).

## Run

```
python examples/BSimReservoirPlanPocketOscSI_D0/run_d0.py
python examples/BSimReservoirPlanPocketOscSI_D0/check_d0.py
```

Standing record:
[`examples/PocketDish/D0_BULK_ORACLE_STANDING.md`](../PocketDish/D0_BULK_ORACLE_STANDING.md).

Equations and SI quotes: [`EQUATIONS.md`](EQUATIONS.md).
Frozen protocol: [`PROTOCOL.md`](PROTOCOL.md).
