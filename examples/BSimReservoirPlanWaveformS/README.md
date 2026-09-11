# WaveformS — sine / square / triangle on the claim dish

Track E1. New named experiment. 16 samples per period so the sine is
an actual sine, not the five-point polygon from Waveform1. Classes
are `sine`, `square`, `triangle`. Phase is identically 0.

This does not convert Waveform1 Overall FAIL into PASS. It is not
Waveform2 / 2b / 2c / 2d. Waveform `GATE_EVIDENCE.md` is not edited.
C1 remains DEFER. E0.3 remains NO_STORY_MOVE.

Variance and power differ across classes. START_U is 0.25 / 0.50 /
0.00. Those leaks are declared. Do not rewrite templates after AUC.

## Status

**u-only construction PASS. Occupancy ALIVE.** System PASS on seeds
101 / 202 / 303: driven block AUC beats Brownian and silent on every
seed (mean `0.7639 ± 0.1187` vs Brownian `0.5139 ± 0.1211`, silent
`0.5000`). MOMENTS `0.8333` and START_U `0.8333` are declared skill;
the paper must say so. Field-only `0.9583` is ≥ driven on seeds 202
and 303: the shapes are already in the plume. Waveform1 Overall
FAIL is unchanged.

See `PROTOCOL.md` and `results/WAVEFORM_S_SCOUT.md`.

## Generate and check

From the repo root:

```
python examples/BSimReservoirPlanWaveformS/generate_waveforms_inputs.py
python examples/BSimReservoirPlanWaveformS/check_waveforms.py --u-only
```

Hashes in `PROTOCOL.md` were recorded before ridge. Do not regenerate
after seeing reservoir AUC. Do not rewrite templates.

## Compile and smoke (authorized after u-only construction PASS)

From this directory:

```
javac -encoding UTF-8 -cp "../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" -d .. BSimReservoirPlanWaveformS.java VoxelAnalyzer.java
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanWaveformS.BSimReservoirPlanWaveformS sim_config_waveforms_smoke.properties
```

Scout seed 101 driven (Brownian / silent reuse Waveform2c CSVs):

```
java -cp "..;../../dist/build;../../lib/core.jar;../../lib/vecmath.jar;../../lib/objimport.jar" BSimReservoirPlanWaveformS.BSimReservoirPlanWaveformS sim_config_waveforms_driven_seed101.properties
```

Do not retune K, n, tau_R, tau_L, source rate, clamp, mortality, flow,
or layout. Claim dish stays CENTER / FLOW=0.
