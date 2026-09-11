# Manuscript — STATIC_HILL of the plume (not for submission)

IEEEtran methods note. Standing numbers only. Does not rewrite any
gate. **Not for journal submission.** Headline is the static-Hill
control, not a reservoir. Do not call this object a reservoir
anywhere in the tex.

## Build figures

```
python examples/PocketDish/manuscript_ieee/plot_manuscript_figures.py
```

## Build PDF (TeX Live / MiKTeX)

```
cd examples/PocketDish/manuscript_ieee
pdflatex two_living_layers
bibtex two_living_layers
pdflatex two_living_layers
pdflatex two_living_layers
```

`IEEEtran.cls` and `IEEEtran.bst` are required (standard TeX Live
`texlive-publishers` / MiKTeX `ieeetran`).

## Referee sentence (this draft)

On this occupied millimetre AHL dish, a linear map of the plume
is at chance on successive-bit XOR (field AUC 0.4677); a static
Hill of the same voxels, \(R\), and \(R\|L\) are all 1.0000; \(L\)
is 0.9927. The living layer's contribution over the linear field
is that Hill. Delay capacity already failed the \(L\)-memory
hypothesis. A 10-tap of the NARMA drive still beats \(R\|L\).

**Falsifier of the headline.** `STATIC_HILL` of the FIELD voxels
fails to match \(R\) or \(R\|L\), or the living output at a given
\(A\) is path-dependent (hysteresis). The second property cannot
hold on this identity; showing it requires a new named object.

**Falsifier of the memory half.** Hill \(R\) reconstructs
contemporaneous or delayed \(u\) better than the AHL map, or
\(R\|L\) test NRMSE on the frozen NARMA-10 \(u\) is at or below
the 10-tap of that drive.

AXIS_HOLDS (sign 5/5, mean NRMSE margin 0.1147) does not change
that sentence: living-versus-field signs hold; STATIC_HILL remains
the headline; the 10-tap remains the NARMA ceiling.

## Honesty locks (do not “fix” in the tex)

- STATIC_HILL is the headline, not a footnote. Do not reopen with
  “living computer,” “living supplies nonlinearity” as a mechanism,
  or “reservoir” as a description of this dish.
- Record the decomposition, not a tie: FIELD 0.4677 → one frozen
  Hill → 1.0000. Cells' marginal over Hill(A) is zero. \(L\) is
  unused (0.9927). `retune_Hill: false`, `second_Hill_K: false`,
  `static_hill_is_java_arm: false`.
- The \(L\)-memory hypothesis is closed on the Jaeger curve
  (MC_fading living 2.07 vs field 2.87; \(\tau_L\) did not extend
  past AHL 303 s). Do not list it as untested.
- Hysteresis at matched \(A\) is the unique remaining discriminator
  of STATIC_HILL. It is not a result in this note. Do not start
  Weber-on-dish Java from this draft. Job 5 remains PARTIAL
  (upper saddle-node unexplained) and is a different object.
- Lorenz and waveform stay living-versus-field reports. They are
  not a licence to drop the 10-tap and not evidence of a reservoir.
- Lane A extras are seed 101. No CHARC / IPC. Jaeger delay curve is
  a named extra on the frozen NARMA $u$; the $10$-tap MC $=9$ is a
  ceiling, not a CHARC number.
- Do not paste HybridDish 0.928 / 0.260 / 0.826 / 0.835 onto Lane A
  tables. Predecessor numbers stay in the HybridDish subsection.
- Do not paste Lane D $\gamma$ into this draft. `write/lane-d-note`
  was not started.
- Waveform: lead with field 0.820, not saturated AUC 1.000.
- Lane B occupation FAILs stay as a short negative. No Pe figure.
  No 90 s extra. CheY STEP_DOWN stays FAIL.
- Still **not for submission**. Framing moved; venue decision is
  separate. Paper 1 stays the peer's object.
