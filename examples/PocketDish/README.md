# PocketDish (Stage99+) — architecture freeze

**Date locked:** 2026-08-19.  
**Horizon:** two months.  
**This package does not execute BSim.**

Paper 1 is `reservoir_new`. Paper 2 evidence is the HybridDish claim
dish (already scored). PocketDish is a **new named device**, not a
retune of either.

| Object | Status in this chat |
|---|---|
| HybridDish claim dish, every `GATE_EVIDENCE.md` | **untouched** |
| `reservoir_new` / paper 1 Results | **not rewritten** |
| Stage99 evaluation (`external_analysis/Stage99_Trapped_Dish_Evaluation.md`) | **untrusted** (wrong AHL clock; “exact Danino blueprint”) |
| PocketDish | Groisman+AC ([`CLAIM_GROISMAN_AC.md`](CLAIM_GROISMAN_AC.md)); T1/T1b/T1c QS **DEAD**; D50 closed **ALIVE**; Flow CSTR **DEAD**; SI bulk **NO_PERIOD**; Mem **CLAMPED** |

## North star

A lab-plausible AC–living-cell reservoir whose **primary readout is
spatiotemporal reporter fluorescence**, in a **quiescent square pocket
opening onto a perfused bus** (Danino/Prindle *layout class*), with ACs
as cited transducers. Two-way only where Lentini-class chemistry exists.
Mechanical pucks are a **separate, cited** readout family — never a
silent replacement for \(L\).

If a mechanism has no paper or no units, it is
`HYPOTHETICAL_DESIGN_ENVELOPE` or it is cut.

## Frozen decisions (one line each)

1. **Geometry = PocketDish-A** (Danino square pocket + horizontal bus).
   Image 1 (1000×500 + vertical right-bus) is a labelled BSim variant,
   not the first living job. [`GEOMETRY_FREEZE.md`](GEOMETRY_FREEZE.md).
2. **Receiver = PocketHill** (HybridDish Hill \(K=1.6\), \(n=2\),
   \(\tau_R=15\) s, \(\tau_L=1500\) s). **Not** Stage 3 Danino 4-ODE.
   Oscillations at hours in a filled trap are a different named object
   (`PocketOsc`), not this backbone.
3. **AHL clocks** = HybridDish \(D=159\) µm²/s,
   \(k=0.0033\) s⁻¹ (\(\tau\approx 303\) s). **Not** Stage 8 ~4.2 h.
4. **Through-flow 8 µm/s through the pocket** stays the **negative
   control** (E0.3 / SweepS5 occupancy DEAD). Danino bus flow is
   *alongside* the open edge.
5. **Primary RC state** = optical fluorescence field. 408-D deaths are
   an optional labelled extra. Grober \(\omega\) is not gene-expression
   computing.

## Files

| File | What it freezes |
|---|---|
| [`CITATION_DOSSIER.md`](CITATION_DOSSIER.md) | papers, DOI, taken vs not, unverified numbers |
| [`GEOMETRY_FREEZE.md`](GEOMETRY_FREEZE.md) | PocketDish-A vs image 1 vs Danino 2010 |
| [`READOUT_STACK.md`](READOUT_STACK.md) | fluorescence / occupancy / optional 408 / pucks |
| [`CHIP_PLAN.md`](CHIP_PLAN.md) | Whole chip investigation; T4/T5 population **EMPTY**; Mem **CLAMPED** |
| [`POCKETNECK_T4_FROZEN_BUILDER_PROMPT.md`](POCKETNECK_T4_FROZEN_BUILDER_PROMPT.md) | T4 \(N(t)\) growth-on (already run) |
| [`FILL_DESIGN.md`](FILL_DESIGN.md) | How to keep a colony: monolayer **EMPTY** → membrane **CLAMPED** |
| [`POCKETMONOLAYER_T5_FROZEN_BUILDER_PROMPT.md`](POCKETMONOLAYER_T5_FROZEN_BUILDER_PROMPT.md) | T5 Danino height+repulsion (already run EMPTY) |
| [`T4_STANDING.md`](T4_STANDING.md) | T4 growth-on occupancy ALIVE; population EMPTY |
| [`T5_STANDING.md`](T5_STANDING.md) | T5 both arms ALIVE perfume, EMPTY garage; weir failed |
| [`MEM_STANDING.md`](MEM_STANDING.md) | Membrane occupancy ALIVE; population CLAMPED; fill exists |
| [`POCKETMEMBRANE_FROZEN_BUILDER_PROMPT.md`](POCKETMEMBRANE_FROZEN_BUILDER_PROMPT.md) | Membrane fill job (already run CLAMPED) |
| [`T3_STANDING.md`](T3_STANDING.md) | BSim field port = T2 class; living field ALIVE; \(N\to 0\) is T4 |
| [`GANTT_TWO_MONTH.md`](GANTT_TWO_MONTH.md) | week-by-week cut and kill criteria |
| [`T0_STANDING.md`](T0_STANDING.md) | T0 occupancy DEAD on the lab-match flush door |
| [`T1_STANDING.md`](T1_STANDING.md) | PocketFill-T1: volumetric LuxI DEAD; plug envelope not a claim |
| [`T2_STANDING.md`](T2_STANDING.md) | PocketNeck-T2: flush DEAD; \(W=50,20,10\) µm ALIVE |
| [`POCKETFILL_T1_FROZEN_BUILDER_PROMPT.md`](POCKETFILL_T1_FROZEN_BUILDER_PROMPT.md) | T1 occupancy builder prompt (already run DEAD) |
| [`HANDOFF_QS_AND_CONSORTIUM.md`](HANDOFF_QS_AND_CONSORTIUM.md) | Program for two new chats: QS occupancy vs consortium |
| [`CHAT_QS_BOOTSTRAP.md`](CHAT_QS_BOOTSTRAP.md) | Paste as the first message of the QS chat |
| [`CHAT_CONSORTIUM_BOOTSTRAP.md`](CHAT_CONSORTIUM_BOOTSTRAP.md) | Paste as the first message of the consortium chat (later) |
| [`QS_PACK_PROVENANCE.md`](QS_PACK_PROVENANCE.md) | Job 1 done: T1 pack ENGINEERING; cited \(d=0.5\) at SI height |
| [`POCKETOSC_D50_FROZEN_BUILDER_PROMPT.md`](POCKETOSC_D50_FROZEN_BUILDER_PROMPT.md) | Job 2 builder (already run); claim ALIVE |
| [`D50_STANDING.md`](D50_STANDING.md) | Closed \(d=0.5\), 1.65 µm optical LA **ALIVE**; T1c replay DEAD |
| [`POCKETOSC_FLOW_FROZEN_BUILDER_PROMPT.md`](POCKETOSC_FLOW_FROZEN_BUILDER_PROMPT.md) | Job 3 builder (already run); open flow DEAD |
| [`FLOW_STANDING.md`](FLOW_STANDING.md) | Closed ALIVE; 180 and 296 µm/min **DEAD**; no identity |
| [`DANINO_SI_MODEL.md`](DANINO_SI_MODEL.md) | TAKEN SI delay-DDE table; \(\mu\) free; not D1g |
| [`POCKETOSC_SI_FROZEN_BUILDER_PROMPT.md`](POCKETOSC_SI_FROZEN_BUILDER_PROMPT.md) | Bulk Fig. 4b twin (already run); all COVER NO_PERIOD |
| [`CHAT_QS_SI_BOOTSTRAP.md`](CHAT_QS_SI_BOOTSTRAP.md) | SI-DDE bootstrap (complete DEAD) |
| [`SI_STANDING.md`](SI_STANDING.md) | SI bulk delay-DDE **NO_PERIOD**; identity not scored |
| [`T1C_STANDING.md`](T1C_STANDING.md) | Closed LuxI DEAD at ENGINEERING pack; D50 is a different garage |
| [`T1B_STANDING.md`](T1B_STANDING.md) | T1b pulse DEAD; bath diagnostic ALIVE only |
| [`POCKETFILL_T1B_FROZEN_BUILDER_PROMPT.md`](POCKETFILL_T1B_FROZEN_BUILDER_PROMPT.md) | T1b QS seed occupancy (already run DEAD) |
| [`DESIGN_TROUBLESHOOT.md`](DESIGN_TROUBLESHOOT.md) | Why T0/T1 died; PocketNeck (small door), not hex-as-magic |
| [`GEOMETRY_NECK.md`](GEOMETRY_NECK.md) | PocketNeck door freeze (T2) |
| [`POCKETNECK_T2_FROZEN_BUILDER_PROMPT.md`](POCKETNECK_T2_FROZEN_BUILDER_PROMPT.md) | T2 neck-width occupancy (already run) |
| [`PROTOCOL.md`](PROTOCOL.md) | integrity, occupancy gates, no NRMSE hunting |
| [`POCKETDISH_TRANSPORT_FROZEN_BUILDER_PROMPT.md`](POCKETDISH_TRANSPORT_FROZEN_BUILDER_PROMPT.md) | **only** the first transport/occupancy job |
| [`refs/README.md`](refs/README.md) | the three user images (copy originals here) |

Visual briefing: open the Cursor canvases
`pocketdish-architecture.canvas.tsx`, `pocketneck-design.canvas.tsx`,
`pocketneck-chip-plan.canvas.tsx`, and `pocketneck-fill-design.canvas.tsx`
beside the chat.
