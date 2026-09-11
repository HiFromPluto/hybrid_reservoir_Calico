# Sweep Step 1 — one-AHL-source layout scout

OFAT from the Narma10b claim dish. Only `field.ahl.source.{x,y,z}`
changes. Acid stays `(300, 375, 5)`. Attractants stay silent. Kinetics,
`u`, windows, and ridge are not retuned. Step 0 FLAT still stands: do
not run higher voxel numbers. This is not Stage 7.

Frozen `u` SHA-256:
`d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`.

Claim `(500, 250, 5)` uses existing Narma10b seed 111 voxels. Do not
rerun it.

Production BSim (seed 111 only unless promoted):

```
javac -encoding UTF-8 -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanSweepS1.java VoxelAnalyzer.java
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepS1.BSimReservoirPlanSweepS1 sim_config_s1_offcentre_driven_seed111.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepS1.BSimReservoirPlanSweepS1 sim_config_s1_offcentre_brownian_seed111.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepS1.BSimReservoirPlanSweepS1 sim_config_s1_wall_driven_seed111.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepS1.BSimReservoirPlanSweepS1 sim_config_s1_wall_brownian_seed111.properties
```

`sim_config_s1_smoke.properties` is HybridDish smoke timing only. Do not
score it as NARMA.
