# Sweep Step 5 — +x flow scout

OFAT from the Narma10b claim dish. Only `flow.speed.um_s` changes.
Claim box stays `1000×500×10`. AHL stays `(500, 250, 5)`, rate `1.28e8`,
`K=1.6`, `tau_L=1500` (`alpha=delta=1/1500`), `tau_R=15`, warmup `18000`.
Readout stays `20×10` + `4×2`. Step 0 FLAT still stands. Steps 1–4 are
closed: do not use wall, offcentre, a different rate, a different
`tau_L`, or a different box.

This is not cell-drift-only. Each tick upwind-advects AHL, acid,
attractant, and repellent with the Stage 6 scheme (left inlet 0, no
recycling). Stokes +x is applied to `ReservoirBacterium` and
`BrownianParticle`. Do not copy Stage 6's 10 µm/s. Frozen speed is 8.

Frozen `u` SHA-256:
`d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`.

Claim `FLOW_SPEED=0` uses existing Narma10b seed 111 voxels. Do not
rerun it.

Production BSim (seed 111 only unless promoted):

```
javac -encoding UTF-8 -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanSweepS5.java VoxelAnalyzer.java
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepS5.BSimReservoirPlanSweepS5 sim_config_s5_flow_driven_seed111.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepS5.BSimReservoirPlanSweepS5 sim_config_s5_flow_brownian_seed111.properties
```

`sim_config_s5_smoke.properties` is HybridDish smoke timing with
`flow.speed.um_s=8`. Do not score it as NARMA.
