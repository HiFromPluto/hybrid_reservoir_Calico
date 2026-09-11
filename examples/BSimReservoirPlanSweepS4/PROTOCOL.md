# Sweep Step 4 — box vs plume scout

OFAT from the Narma10b claim dish. Dish bounds change. Field spacing
stays ~20 µm. Pop and clamp scale with area so density is not
confounded. Relative source layout stays the claim fractions. AHL rate
stays `1.28e8` — do not scale the source with area. Acid rate stays
`2e11`. `K=1.6`, `tau_L=1500` (`alpha=delta=1/1500`), `tau_R=15`,
warmup `18000`. Readout stays `20×10` + `4×2`. Step 0 FLAT still
stands. Steps 1–3 are closed: do not use wall, offcentre, a different
rate, or a different `tau_L`.

Frozen `u` SHA-256:
`d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`.

Claim `1000×500×10` uses existing Narma10b seed 111 voxels. Do not
rerun it.

Production BSim (seed 111 only unless promoted):

```
javac -encoding UTF-8 -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanSweepS4.java VoxelAnalyzer.java
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepS4.BSimReservoirPlanSweepS4 sim_config_s4_small_driven_seed111.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepS4.BSimReservoirPlanSweepS4 sim_config_s4_small_brownian_seed111.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepS4.BSimReservoirPlanSweepS4 sim_config_s4_large_driven_seed111.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepS4.BSimReservoirPlanSweepS4 sim_config_s4_large_brownian_seed111.properties
```

`sim_config_s4_smoke.properties` is HybridDish smoke timing on the small
box only. Do not score it as NARMA.
