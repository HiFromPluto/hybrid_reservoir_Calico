# WaveformS protocol (frozen before ridge)

Track E1. Recognizable sine / square / triangle classification on the
frozen HybridDish claim dish, 16 samples per period. This is a **new
named experiment**. It does not convert Waveform1 Overall FAIL into
PASS. It is not Waveform2 / 2b / 2c / 2d.

Waveform `GATE_EVIDENCE.md` is not edited. Waveform2c is read/reused
for silent and Brownian CSVs only. Kinetics, layout, flow, source
rate, clamp, and mortality are not retuned. Claim dish stays CENTER /
`FLOW=0`. C1 remains DEFER. E0.3 remains NO_STORY_MOVE.

## Scientific question

Can an occupied one-way AHL → R → L dish classify recognizable
sine / square / triangle blocks, above Brownian, silent, and scalar
moment baselines, when each period is sampled 16 times?

Field-only AHL is a carrier diagnostic, not a kill switch. Square has
higher power than sine or triangle. START_U is 0.25 / 0.50 / 0.00.
Those leaks are expected. Report them. Do not rewrite templates after
AUC.

## Dish, timing, seeds (copied — do not change)

- `1000 × 500 × 10` µm, `dt=0.05`, `FLOW_SPEED=0`, `NO_FLUX`
- Field `50×25×1`, readout `20×10` + `4×2`
- `K=1.6`, `n=2`, `tau_R=15`, `tau_L=1500`
- AHL `1.28e8` at `(500, 250, 5)`. Acid `(300, 375, 5)` held 0.5 in
  analysis windows, off in warmup
- Attractants silent. Clamp `K=2000`, `INITIAL_POP=1800`
- Warmup `18000 s`, window `300 s`, pulse `75 s`, sample every `20 s`,
  **400 windows**
- Official 408. Brownian Den 20×10. Field-only driven AHL_uM 20×10
- Seeds `101 / 202 / 303`. Scout seed **101** first. 202/303 only if
  101 is ALIVE and driven block AUC beats Brownian and silent
- Last sample `399;15;299.95`. CSV 400 / 6400 / 6400

## Frozen task geometry

16 windows per block = one period. 25 blocks = 400 windows.
Chance accuracy = 1/3. Chance macro OVR AUC = 0.5.
Classes: `0=sine`, `1=square`, `2=triangle`.

| Split | Blocks | Windows |
|---|---|---|
| Washout | 0..2 | 0..47 |
| Train | 3..18 | 48..303 |
| Test | 19..24 | 304..399 |

Inner validation = last 4 train blocks (15..18), windows 240..303.
Lambda and standardization never see test. Inner standardization uses
inner-train blocks only (3..14). Final fit uses all 16 train blocks.

Class counts, placed then shuffled with `random.Random(20260818)`
inside each split (already done; vector frozen, do not reshuffle):

- washout 3: 1,1,1
- train 16: 6 sine, 5 square, 5 triangle
- test 6: 2,2,2

Length-25 class vector (comma-separated, no spaces):

```
2,0,1,1,2,0,0,1,0,2,0,1,0,1,1,2,0,2,2,0,2,2,1,1,0
```

SHA-256 of that ASCII string:
`40adec9e1e69314a1250a723173d5d877bc4540e15ad17ae918fac735f081a46`

Phase is identically 0 on every block. No circular shift.

Length-25 phase vector:

```
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
```

SHA-256 of that ASCII string:
`ecd9d5de84e1338de30087c5f516b0d872152955d3ca059810cd03aa7b85fb6c`

## Frozen templates (12-decimal strings, w = 0..15)

Sine: `u = 0.25 + 0.25 * sin(2 * pi * w / 16)`

Square: eight samples at 0.5, then eight at 0.0

Triangle: `u = 0.5 * (1 - abs(2 * (w/16) - 1))`

Clip is unnecessary: every value is already in `[0, 0.5]`.
No extra noise. Do not change templates after seeing AUC.

### Block moments (from the 12-decimal strings)

| class | mean | variance | pop SD | power | min | max |
|---|---|---|---|---|---|---|
| sine | 0.250000000000 | 0.031250000000 | 0.176776695297 | 0.093750000000 | 0.000000000000 | 0.500000000000 |
| square | 0.250000000000 | 0.062500000000 | 0.250000000000 | 0.125000000000 | 0.000000000000 | 0.500000000000 |
| triangle | 0.250000000000 | 0.021484375000 | 0.146575492494 | 0.083984375000 | 0.000000000000 | 0.500000000000 |

Means, min, and max match. Variance / power do **not**. That is
declared. START_U is 0.25 / 0.50 / 0.00.

## Input files

Generated, then hashed, **before** any ridge:

- `input_ahl_waveforms_400.txt` — 400 lines of 12-decimal `u`
- `input_acid_held05_400.txt` — byte copy from Waveform2c
- `waveforms_labels.csv` — `n;block;class;phase;u` with `phase=0`
- `class_phase_schedule.txt` — class vector and all-zero phase vector

SHA-256 of the 400-line 12-decimal `u` body (newline after every
value, including the last):

`12941fff63a4b02641803b3c11fdb89cd812832b54a3e06c6fee5f9a4b9d7c92`

First 16 values are the triangle template (class vector starts with 2).

## File hashes (recorded before ridge)

SHA-256 of each file as stored on disk, recorded after input
generation and before `--u-only`. The 400-line `u` body hash above
is the construction gate; the AHL file hash below includes the
comment header.

| File | SHA-256 |
|---|---|
| BSimReservoirPlanWaveformS.java | `4c1732d8fe85ef24f8c07c2f3d0363b1ce876213bc1b82c29138d6a2560b6cba` |
| VoxelAnalyzer.java | `23734a2c7cbf2f60df97bcb2898b00aa7c35367d6a95383fcda2e3476eede5db` |
| check_waveforms.py | `feb1c9695ee1c03c349ff4195f3d085d9e1237af3cb3b72a7805d4e3e8490836` |
| generate_waveforms_inputs.py | `734d4798238413e7280a2779a69a9d111f0306ae6f12a4176b87853adf61f0c5` |
| sim_config_waveforms.properties | `a372396f3a8e9eaae0dc20314b2699450033a8f937bb0e2d5b70149a6dea93d9` |
| sim_config_waveforms_smoke.properties | `244bf42368f81efbd76770ccd1dae766c372458b9dd5e61871a4bdccfcba83c3` |
| sim_config_waveforms_driven_seed101.properties | `b2fa19f918f803c6d436caf8ee7c6f82e8c40c5f61db085c5178440988ee910c` |
| sim_config_waveforms_driven_seed202.properties | `4558e6937a2934a51e4a32a91bbb13afcf008710623d0df03b3cd573f0145c73` |
| sim_config_waveforms_driven_seed303.properties | `c24ab93b39edc59a030a87a8703cfdef66e55caa71b61970e95d9faaf0082420` |
| input_ahl_waveforms_400.txt | `fab0af071e453a58130c4d8f04fb6848e95e9ac23503015e72a6840408843466` |
| input_acid_held05_400.txt | `9223d8ad6add51b023d2b495d5a6ebdf3703b86133b3d2531705126af279e3e1` |
| waveforms_labels.csv | `1f73e5f00cc571c9a75ffb59800a273ceeec1e88480019e888254110cb11a3b3` |
| class_phase_schedule.txt | `28f679d19a014f68bbc1b8ddec74273c9c7c3a9aea9cf4c647842cc0328bc97f` |

## u-only audit (mandatory stop, before BSim)

Block-level one-vs-rest ridge, same lambda rules as Waveform2c
(16 windows/block, 25 blocks).

Baselines: MEAN_ONLY, POWER_ONLY, VARIANCE_ONLY, MOMENTS
(mean + variance + power), START_U, RAW_U_16.

**STOP (construction bug), no BSim:**

- MEAN_ONLY test macro AUC ≥ 0.70 → means are not matched
- RAW_U_16 test macro AUC < 0.90 → the three 16-vectors are not
  linearly distinct

**Not a stop (expected skill, must print):**

- POWER, VARIANCE, MOMENTS, START_U may be high. Record them. If
  MOMENTS or START_U ≥ 0.70, the paper must say so. Still run BSim if
  the two STOP rules pass.

### u-only audit result

**PASS.** Block-level test macro OVR AUC:

- MEAN_ONLY 0.5000 (stop if >= 0.70)
- POWER_ONLY 0.8333 (record)
- VARIANCE_ONLY 0.8333 (record)
- MOMENTS 0.8333 (record)
- START_U 0.8333 (record)
- RAW_U_16 1.0000 (stop if < 0.90)

MOMENTS and/or START_U are >= 0.70. That is declared skill. The paper must say so. Not a construction stop.

See `results/u_only_baselines.md`.

BSim is authorized to the smoke, then seed-101 scout.


## Ridge

Preserve Waveform2c OVR behaviour: three binary heads, lambda grid
`{1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6}`, unregularized intercept,
inner-train standardization, lambda by highest macro OVR validation
AUC, ties → larger lambda.

**Primary:** mean the 16 window OVR scores inside each block, then
macro OVR AUC and accuracy on the 6 test blocks.

**Secondary:** window-level macro OVR AUC on test windows; per-class
OVR AUC; confusion counts; accuracy vs 1/3 (not a gate).

Biology is **not** required to beat RAW_U_16 or field.

Occupancy DEAD if `mean_R < 0.05` or `|r(mean_R, u)| < 0.5` on all
400 windows. If DEAD: NOT_SCORED. Do not raise `Jmax` or `K`.

## Arms

1. Driven. New BSim. AHL file above. Acid held 0.5.
2. Brownian. Reuse Waveform2c Brownian CSVs after 400 / 6400 / 6400
   and last sample `399;15;299.95`. Particles do not sense AHL
   (`BrownianParticle` has no AHL field; voxels are `Den_` only).
3. Silent. Reuse Waveform2c silent CSVs under the same check.
   200-window Waveform1 / Stage 6 silent is forbidden.

## Standing after scores

Write `results/WAVEFORM_S_SCOUT.md` and `results/waveforms_scout.csv`.
Do not invent a new Overall for Waveform1. Do not edit Waveform2c.

Occupancy ALIVE on 101/202/303. System PASS on every seed (driven
block AUC > Brownian and > silent). MOMENTS `0.8333` and START_U
`0.8333` are declared leaks; the paper must say so. Field-only
`0.9583` is ≥ driven on seeds 202 and 303: the shapes are already
in the plume. Waveform1 Overall remains FAIL.

## Result hashes (after scoring)

| File | SHA-256 |
|---|---|
| results/u_only_baselines.md | `9917b7e5c095cf1921cab917000293f1c1e9714e6b255042d6c19efde945e873` |
| results/u_only_baselines.csv | `1826737a7c71613e5f49192e322b990b12768f9beff8ff238bdf39b3b6790889` |
| results/WAVEFORM_S_SCOUT.md | `023edfd28deb71397a907b8744032916fd021851c97b80635502ac2f0e40e8d1` |
| results/waveforms_scout.csv | `9d5f932a70afbf046ca313fd5bb9fdf44c75d805884fd4b813f55b92fa5f63a6` |
| results/waveforms_driven_seed101/run_status.txt | `2268ac23ed002a5002b267b2c720209294ec6c1cfa8a028543d915b9ed9e634e` |
| results/waveforms_driven_seed202/run_status.txt | `0fbd4f30048916af222896cfc0de06b4bc839816ae6f789d083547d336059a95` |
| results/waveforms_driven_seed303/run_status.txt | `5394d0596c55b664d0338896160b541c995032ce0eb1e5d5fbefac3b029088dc` |

Checker SHA-256 after the post-score scout writer (mean ± s.e. and
all-seed System lines only; templates and ridge unchanged):
`db3041d5afe1108cba4e340d38b1d128082c82d15c2567224a88de2ec53896bc`.
Pre-ridge checker hash remains in the table above.

## Checker

```
python examples/BSimReservoirPlanWaveformS/generate_waveforms_inputs.py
python examples/BSimReservoirPlanWaveformS/check_waveforms.py --u-only
python examples/BSimReservoirPlanWaveformS/check_waveforms.py
```
