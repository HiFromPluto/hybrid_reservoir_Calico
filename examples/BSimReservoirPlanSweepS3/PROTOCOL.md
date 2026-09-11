# Sweep Step 3 — luminescence lag tau_L scout

OFAT from the Narma10b claim dish. Only `tau_L` changes, via
`luminescence.alpha = luminescence.delta = 1/tau_L`. Changing alpha
alone is a gain change and is a fail of this package. AHL stays
`(500, 250, 5)`, rate `1.28e8`, `K=1.6`, `tau_R=15`, warmup `18000`.
Step 0 FLAT still stands. Steps 1–2 are closed: do not use wall,
offcentre, or a different rate.

Frozen `u` SHA-256:
`d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`.

Claim `tau_L=1500` uses existing Narma10b seed 111 voxels. Do not rerun it.
Brownian is not rerun (`BrownianParticle` has no L).

Production BSim (driven seed 111 only unless promoted):

```
javac -encoding UTF-8 -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanSweepS3.java VoxelAnalyzer.java
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepS3.BSimReservoirPlanSweepS3 sim_config_s3_fast_driven_seed111.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepS3.BSimReservoirPlanSweepS3 sim_config_s3_slow_driven_seed111.properties
```

`sim_config_s3_smoke.properties` is HybridDish smoke timing only. Do not
score it as NARMA.
