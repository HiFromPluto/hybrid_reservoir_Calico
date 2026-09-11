# Sweep Step 2 — AHL source-rate occupancy scout

OFAT from the Narma10b claim dish. Only `field.ahl.source.rate` changes.
AHL position stays `(500, 250, 5)`. Acid stays `(300, 375, 5)`.
Attractants stay silent. `K` stays 1.6. Kinetics, `u`, windows, and
ridge are not retuned. Step 0 FLAT still stands. Step 1 is closed:
do not use wall or offcentre.

Frozen `u` SHA-256:
`d6c0cdfbe4dc713bf6b6041695c55aad500ddf2d68979cdfc7e3f5ee870e4c1e`.

Claim rate `1.28e8` uses existing Narma10b seed 111 voxels. Do not rerun it.

Rates: low `64000000` (0.5×), high `256000000` (2×). Do not try 0.25× or 4×.

Production BSim (seed 111 only unless promoted):

```
javac -encoding UTF-8 -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanSweepS2.java VoxelAnalyzer.java
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepS2.BSimReservoirPlanSweepS2 sim_config_s2_low_driven_seed111.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepS2.BSimReservoirPlanSweepS2 sim_config_s2_low_brownian_seed111.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepS2.BSimReservoirPlanSweepS2 sim_config_s2_high_driven_seed111.properties
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanSweepS2.BSimReservoirPlanSweepS2 sim_config_s2_high_brownian_seed111.properties
```

`sim_config_s2_smoke.properties` is HybridDish smoke timing only. Do not
score it as NARMA.
