# Plan Stage 3B — spatial fast AHL receiver

Stage 3B replaces only the per-cell receiver. The direct AHL source, field,
300 s windows, 75 s pulse, population mechanisms, and spatial grids are copied
forward into `BSimReservoirPlanStage3B`. The archived Stage 3 package remains
unchanged, and Stage 4 stays blocked.

## Benchmark provenance

`benchmark_receiver.py` replays the three recorded Stage 3 `voxels.csv` traces.
It does not run BSim or modify the AHL field. Each of the 200 recorded AHL voxel
traces drives:

`dR/dt = (C^n/(K^n + C^n) - R)/tau`

with `n=2`, explicit `tau`, exact first-order updates, and AHL converted from
molecules/um3 to uM using `602.2`. Outputs are weighted by the recorded voxel
cell densities. The initial receiver state is the Hill equilibrium of the first
recorded sample so receiver startup is not mistaken for input response.

The Hill form is a phenomenological LuxR/AHL receiver, not a retuned Danino
model. LuxR requires 3OC6-HSL for transcriptional activation and forms the
signal-dependent DNA-binding complex; purified LuxR/AHL binding is reversible
(Urbanowski et al., J Bacteriol 2004, DOI `10.1128/JB.186.3.631-637.2004`).
Wild-type LuxR 3OC6-HSL half-activation has been reported near 10–24 nM
(Collins et al., Mol Microbiol 2005, DOI
`10.1111/j.1365-2958.2004.04437.x`; Hawkins et al., AEM 2007,
PMCID `PMC2074898`). Those sensitivities are useful provenance, but the
recorded Stage 3 field is much higher; Stage 3B therefore benchmarks explicit
K values against the actual trace range rather than claiming a literature K
fits this operating point.

## Predeclared gate

All conditions are required:

1. `r(input, mean response) > 0.5` in every recorded permutation.
2. Density-weighted mid-input occupancy is 0.2–0.8.
3. Mean response is ordered low < mid < high.
4. Analytic `t95 = -ln(0.05)*tau <= 300 s`.

The deterministic benchmark selected and passed the predeclared receiver
`K=1.6 uM`, `n=2`, `tau=15 s`, analytic `t95=44.936 s`. The spatial Java integration
uses exactly those values without retuning. `K` is trace-matched and is not
claimed as a literature EC50. See `MANIFEST.md` for the production protocol,
output semantics, and evidence gate.

Closure used three new balanced holdout permutations with fresh seeds
`44/55/66` plus a seed-77 source-off control. All holdout gates passed without
parameter changes. Claims remain limited to a phenomenological, trace-matched
receiver; `t95` is calculated from the exact equation, not measured.
