# Occupied living-layer reservoir — briefing for PI and collaborators

**Status:** ready to write from gated HybridDish / `bsim_clean` runs.  
**Audience:** Murat, Maryam, and anyone who has not lived in the gate files.  
**This manuscript:** Ceylin lead. Conference is acceptable; the evidence is aimed at ACS Synthetic Biology / Cell Systems / PLoS Computational Biology class venues.  
**Not this manuscript:** Maryam’s 6 s platform paper (`reservoir_new`). Scores are not mixed.

---

## 1. One-paragraph pitch

We have an in silico millimetre *E. coli* monolayer in which **artificial cells write a time-varying AHL command**, the bacteria **occupy a graded quorum-sensing receiver**, and a **slow luminescent reporter** plus a linear readout perform reservoir computing on **biological clocks** (5-minute windows, not 6 seconds). Matched controls show that the living readout uses the input rather than cell density or an empty dish. On NARMA-10 the intracellular state also beats a map of the AHL plume itself — weakly — and still loses to a delay line of the raw command. Smooth tasks (Mackey–Glass, waveforms) are already solved by the chemical field. That combination is the paper: **an occupied living layer that can compute, with the plume and the delay line on the same figure**, not a claim that bacteria outperform chemistry in general.

---

## 2. How this sits on the first paper

Maryam’s manuscript is the **hybrid architecture** paper: five programmable artificial cells, a flowing millimetre dish, attractant/repellent fields, voxel readout, layout and composition as design axes. It answers: *can this object be built and scored as a physical reservoir?*

This manuscript is the **living-layer** paper on a rebuilt operating point of the same idea. It answers: *when the receiver is actually on, and the window matches a reporter, do the cells add anything that density and the plume do not already contain?*

| | Platform paper (Maryam) | This paper (occupied dish) |
|---|---|---|
| Dish | 1000 × 500 × 10 µm | Same millimetre class |
| Clock | 6 s windows | **300 s windows, 75 s AHL pulse** |
| QS | Sampled channel; intracellular switch not the computer | **Occupied Hill receiver**, mean \(R \approx 0.2\) |
| Artificial cells | Five vesicle-style secretors (att/rep mix) | **One AHL write** + held acid bias |
| What is scored | Hybrid state (fields + bacteria) | **\(R\), \(L\), input-driven deaths** vs Brownian / silent / **field** / delay line of \(u\) |
| Contribution | Name the architecture and measure it as built | Isolate the living layer with controls the field cannot fake |

Nothing here cancels her paper. It is the experiment her Future Work already points to: occupied receiver, minutes not seconds, living vs carrier vs delay line.

---

## 3. What is on the dish

**Domain.** 1000 × 500 × 10 µm monolayer, no through-flow on the claim dish, conservative chemical transport (N0 kernel). Population spawned near a clamp (~1800, cap 2000) so warmup is for the reporter, not an exponential climb. Kinetics were **frozen before scores** (contrast with parameter-scanning active-matter RC; see Gaimann & Klopotek 2025).

### Artificial cells — what they do

They are **immobilized chemical writers**, not motile computers.

- **Claim-dish source (A0):** a point flux \(J = J_{\max}\,u(t)\) of AHL at the dish centre. This is an ideal pipette: command in, molecules out. It is **not** a transcription–translation vesicle.
- **Device-shaped envelope (A1):** a saturating Hill gate \(g(u)\) then \(J = J_{\max} g(u)\), payload-matched to A0. This is the literature-shaped *device* model we are willing to put in Methods: leak, saturation, a static input–output curve. Swapping A1 for A0 does **not** change the NARMA story (\(\lvert\Delta\rvert \approx 0.03\)). It is labelled a **hypothetical design envelope**, not a digital twin of Lentini vesicles.
- **Acid source:** held at a constant command so pH-dependent death is a bias, not a second message.

Lentini et al. (2014, 2017) remain the experimental citation for artificial–natural chemical communication (theophylline → IPTG; two-way AHL). Their vesicles run on **hours-scale TX–TL**. We do not simulate that chassis, and we do not have a laboratory command-to-flux curve. Until that curve exists, the honest AC sentence is: **one-way chemical write, A0 on the claim dish, A1 as the device envelope.**

### Bacteria — what they do (the focus is \(R\) and \(L\))

Simulated agents are motile *E. coli* with a **Lux-type AHL receiver** and a **slow reporter**. Native SdiA / AI-2 / AI-3 eavesdropping is not modelled.

**\(R\) — receptor occupancy (not light).**  
Local AHL concentration \(C\) (µM) drives a graded Hill occupancy

\[
\frac{\mathrm{d}R}{\mathrm{d}t} = \frac{1}{\tau_R}\left(\frac{C^{n}}{K^{n}+C^{n}} - R\right),
\quad K = 1.6~\mu\mathrm{M},\; n = 2,\; \tau_R = 15~\mathrm{s}.
\]

\(R \in [0,1]\) is “how on the receptor is.” \(K\) was chosen so the plume **sits on the Hill** (occupied). It is a LuxR-shaped phenomenological receiver, not Danino’s four-ODE genetic clock (Danino et al. 2010) and not a measured nM LuxR EC50. If mean \(R\) is below ~0.05, we do not report a living-layer score.

**\(L\) — the reporter a camera would see.**  
In the code this state is named luminescence:

\[
\frac{\mathrm{d}L}{\mathrm{d}t} = \frac{R - L}{\tau_L},\quad \tau_L = 1500~\mathrm{s}
\]

(\(t_{95} \approx 75\) min). \(L\) is a **slow optical stand-in** for a transcriptional reporter downstream of LuxR–AHL, so that \(L\) is not a duplicate of \(R\). \(\tau_L\) is a reporter clock, **not** luciferase biochemistry.

Growth (saturating media, ~30 min doubling), a population clamp, and acid-dependent death are on the dish so mass and mortality exist. Run-and-tumble chemotaxis is implemented (Berg) but wired to a **silent** attractant field: it is in Methods as present, **not** the NARMA encoder.

### What machine would see this — and is that practice?

You do **not** image \(R\). You image **\(L\)**: GFP under a fluorescence microscope or plate reader, or luxCDABE on a CCD / PMT. That is standard QS-sensor practice. The closest experimental cousins are Prindle et al. (2012) *biopixels* — *E. coli* carrying luxI / aiiA / luxCDABE in a microfluidic array, bioluminescence as the readout, chemical input, spatial coupling — and Danino et al. (2010) synchronized genetic clocks in a trapped monolayer. Our 300 s windows and slow \(L\) match that world. They do **not** match a 6 s window or a growth-OD reservoir (Ahavi et al. 2026).

The computational readout is a **linear ridge** on a 408-dimensional voxel state: spatial means of \(R\) and \(L\) (20 × 10) plus a few input-driven death bins. The reservoir is untrained; only the readout is trained (Jaeger 2001; Tanaka et al. 2019).

---

## 4. Benchmarks — which, why, and the honest sentence

Nulls on every task: **Brownian** (inert tracers, density only), **silent** (no AHL command), **field-only** (voxel AHL, no \(R\)/\(L\)), and a **delay line of \(u\)** where it is the right ceiling.

| Task | Why it is there | Result to say out loud |
|---|---|---|
| **NARMA-10** | Needs fading memory and nonlinearity; the classic RC stress test (Atiya & Parlos 2000). | Driven NRMSE **0.93** vs Brownian/silent **1.16** vs field **1.03** vs 10-tap of \(u\) **0.68**. Occupied. Weak absolute skill (\(R^2 \approx 0.15\)). **This is the living-layer exception.** Holds on average across 11 independent drives (9/11; two occupied losses kept). |
| **Mackey–Glass** | Standard one-step chaotic prediction (Mackey & Glass 1977). | Beats density; **field already solves** one-step (~0.02 vs driven ~0.26). The plume is the computer. |
| **Waveforms** | Can the dish tell temporal shapes (order; sine / square / triangle)? | Beats Brownian/silent; **field ≥ driven**. Moments of \(u\) already classify some shapes. |
| **Lorenz-63** | Chaotic multivariate drive (Lorenz 1963). | Beats nulls; does **not** beat field or a linear AR of \(x\). The clock lives in the input. |

We do not lead with “cells beat chemistry.” We lead with: **cells beat density; NARMA is the one task where they also beat the plume; everything smooth is in the carrier.** That is the Inubushi–Yoshimura point (2017): useful RC is a balance of memory and nonlinearity, not a trophy NRMSE.

Physical-reservoir neighbours we cite, not twin: photonic/mechanical PRC (Tanaka et al. 2019); a self-organizing chemical network (Baltussen et al. 2024); synthetic active particles (Wang & Cichos 2024); engineered bacterial consortia as perceptrons (Li et al. 2021); unmodified *E. coli* growth as a living reservoir (Ahavi et al. 2026). **What we add** is an **artificial-cell write + occupied intracellular receiver + field and delay-line controls** on a millimetre monolayer. Ahavi scores bulk growth; we score \(R\) and \(L\). Lentini demonstrates AC–bacteria chemistry; we score computation with the living layer isolated.

The simulator is BSim (Gorochowski et al. 2012), with a conservative field kernel and a frozen protocol.

---

## 5. What we are ready to write *now*

A self-contained methods-and-results manuscript (or a strong conference paper with the same spine):

1. Architecture figure: AC write → AHL field → \(R\) → \(L\) → linear readout; null arms on the same dish.  
2. Occupancy: mean \(R\), \(r(R,u)\), ALIVE on every scored seed.  
3. NARMA vs Brownian / silent / field / \(u\)-taps, including the 11-drive hold.  
4. Mackey–Glass and waveforms as the **carrier diagnostic** (field wins).  
5. A0 vs A1: device envelope in Methods; no new claim dish.  
6. Limitations in one paragraph: \(K\) and \(\tau_L\) are phenomenological; A1 is not a lab \(u\to J\) curve; chemotaxis is present and not the encoder; NARMA loses to a tap of \(u\); this is in silico.

Kinetics will not be retuned after seeing NRMSE. Closed FAILs (e.g. waveform vs field as a general shape-classifier; Danino 4-ODE as this dish’s 300 s receiver) stay closed.

---

## 6. Authorship and venue (plain)

- **Maryam** remains first on the 6 s hybrid-architecture manuscript.  
- **This occupied-dish manuscript is Ceylin first.** It uses a different clock, a different receiver, different controls, and different numbers. It can go out as a conference paper without waiting for a journal cycle; the same spine is what we would send to ACS Synth. Biol. / Cell Systems / PLoS Comput. Biol.  
- We do not paste HybridDish NRMSE into Maryam’s tables, and we do not paste her CHARC/6 s scores into this Results section.

---

## References (for the briefing)

Atiya, A. F. & Parlos, A. G. (2000). *IEEE Trans. Neural Netw.* **11**, 697–709.  
Ahavi, P. et al. (2026). Living bacterial reservoir computers. *Cell Syst.* **17**, 101654.  
Baltussen, M. G. et al. (2024). Chemical reservoir computation in a self-organizing reaction network. *Nature* **631**, 549–555.  
Danino, T. et al. (2010). A synchronized quorum of genetic clocks. *Nature* **463**, 326–330.  
Gaimann, M. U. & Klopotek, M. (2025). How to define, find and characterize robust strategies in the parameter space of active matter reservoir computing. Preprint / paper as cited in the group.  
Gorochowski, T. E. et al. (2012). BSim: an agent-based tool for modeling bacterial populations. *PLoS ONE* **7**, e42790.  
Inubushi, M. & Yoshimura, K. (2017). Reservoir computing beyond memory–nonlinearity trade-off. *Sci. Rep.* **7**, 10199.  
Jaeger, H. (2001). The “echo state” approach. GMD Report 148.  
Lentini, R. et al. (2014). Integrating artificial with natural cells. *Nat. Commun.* **5**, 4012.  
Lentini, R. et al. (2017). Two-way chemical communication between artificial and natural cells. *ACS Cent. Sci.* **3**, 117–123.  
Li, X. et al. (2021). Synthetic neural-like computing in microbial consortia. *Nat. Commun.* **12**, 3139.  
Lorenz, E. N. (1963). Deterministic nonperiodic flow. *J. Atmos. Sci.* **20**, 130–141.  
Mackey, M. C. & Glass, L. (1977). Oscillation and chaos in physiological control systems. *Science* **197**, 287–289.  
Prindle, A. et al. (2012). A sensing array of radically coupled genetic ‘biopixels’. *Nature* **481**, 39–44.  
Tanaka, G. et al. (2019). Recent advances in physical reservoir computing. *Neural Netw.* **115**, 100–123.  
Wang, X. & Cichos, F. (2024). Harnessing synthetic active particles for physical reservoir computing. *Nat. Commun.* **15**, 774.

Internal sources for numbers (not for the bibliography): `examples/HybridDish/WHAT_IS_ESTABLISHED.md`, `PROTOCOL.md`, Narma10b / C1 / BenchA / WaveformS scout reports.
