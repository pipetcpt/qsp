# MOG Antibody-Associated Disease (MOGAD) — QSP Model

> A quantitative systems pharmacology model of MOGAD built around one awkward
> observation: **rituximab lowers the MOG-IgG titre more than IVIG does, and
> prevents far fewer relapses.** Any model in which relapse risk simply tracks
> antibody concentration gets that backwards. This one gets it right, and the
> two structural choices that make it work — a short-lived CD20-negative
> antibody source, and a *threshold* hazard in functionally active antibody —
> then also predict why the oligodendrocyte survives, why serum GFAP stays low,
> and why steroid tapers buy time rather than remission.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`mogad_qsp_model.dot`](mogad_qsp_model.dot) |
| 🖼️ Map (SVG) | [`mogad_qsp_model.svg`](mogad_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi) | [`mogad_qsp_model.png`](mogad_qsp_model.png) |
| ⚙️ mrgsolve ODE model | [`mogad_mrgsolve_model.R`](mogad_mrgsolve_model.R) |
| 🐍 Executable reference implementation | [`mogad_reference_impl.py`](mogad_reference_impl.py) |
| 📊 Shiny dashboard | [`mogad_shiny_app.R`](mogad_shiny_app.R) |
| 📚 References (55, all PMID-verified) | [`mogad_references.md`](mogad_references.md) |

**Every simulated number on this page is reproduced by
`python3 mogad_reference_impl.py`** — a dependency-free RK4 twin of the same 34
ODEs with the same parameters. Nothing here is asserted from intuition.

---

## 1. The disease in one paragraph

MOGAD is an antibody-mediated demyelinating disease of the central nervous
system. Class-switched **IgG1 against myelin oligodendrocyte glycoprotein**
reaches the CNS through a permeabilised blood–brain barrier and binds MOG, which
sits on the **outermost lamella of the myelin sheath** — so the epitope is
directly accessible without the antibody ever entering a cell. Bound IgG1 fixes
complement (C1q → C4 → C3 → C5b-9) and opsonises myelin for phagocytosis by
microglia and monocyte-derived macrophages. Attacks present as optic neuritis
(often bilateral, with disc oedema and severe acuity loss), longitudinally
extensive myelitis with conus predilection, ADEM in children, or unilateral
cortical encephalitis with seizures (FLAMES). Recovery is usually good, relapses
occur in roughly half of patients, and — unlike multiple sclerosis — **there is
no secondary progressive phase**: disability accrues only in steps, from attacks.

## 2. The three structural commitments

Everything in the model follows from three choices. Each is stated here so it can
be attacked.

### Commitment 1 — the antibody source is short-lived and largely CD20-negative

The pathogenic IgG1 comes from a **plasmablast pool with a 7-day half-life**
(`KOUT_PB = 0.14/day`) sustained by an IL-6-dependent survival niche. Only part
of its generation is fed by CD20+ memory B cells; the parameter `FRAC_ESC = 0.60`
is the share that rituximab structurally cannot reach.

Consequences that fall out rather than being imposed: titres drop within days of
steroids or plasma exchange; they rebound over about a month once suppression
stops; and rituximab has a hard ceiling on how far it can lower them no matter
how completely it depletes CD20+ cells.

### Commitment 2 — myelin is destroyed, the oligodendrocyte survives

Myelin (`MYEL_*`) and axons (`AXON_*`) are separate states, and the
remyelination term is **multiplied by the surviving oligodendrocyte fraction
`OL`**. In MOGAD the oligodendrocyte death coefficient `KOLD` is small, matching
the neuropathology (Höftberger 2020, Takai 2020: demyelination with preserved
oligodendrocytes and intact AQP4).

Setting `AQP4_MODE = 1` multiplies oligodendrocyte death and astrocytic GFAP
release, converting the *same* model into an AQP4-IgG NMOSD comparator. This is
the model's built-in control experiment: identical antibody kinetics, identical
barrier, identical effectors, identical steroids — one switch.

| Identical attack, identical IVMP | MOGAD | AQP4-IgG comparator |
|---|---|---|
| Minimum surviving oligodendrocytes | **0.88** | **0.07** |
| Visual acuity at 1 year (logMAR) | 0.37 (~20/47) | 0.79 (~20/122) |
| RNFL at 1 year | 75 µm | 63 µm |
| Peak serum GFAP | 150 pg/mL | 3661 pg/mL |
| Peak serum NfL | 160 pg/mL | 162 pg/mL |

The recovery gap and the GFAP discrimination — both well documented clinically
(Stiebel-Kalish 2017; Marignier 2025) — come entirely from *whether the target
cell survives the attack*. Note the NfL columns are nearly identical: the model
says these two diseases differ in **glial** injury, not axonal injury, which is
exactly the reported dissociation.

Commitment 2 also forces a subtler feature. To fit a severe nadir *and* good
final acuity *and* marked RNFL thinning simultaneously, visual deficit has to
be split into reversible conduction block (proportional to demyelination) plus
irreversible axonal loss **beyond a functional reserve** (`AXRES = 0.15`).
Without that reserve term no parameterisation fits all three.

### Commitment 3 — relapse hazard is a threshold function of *active* antibody

Hazard is `LAM0 × Hill(AB_EFF, AB50 = 1.8, h = 3) × amplification`, where
`AB_EFF` is the **functionally active** antibody, not the measured titre. IVIG
lowers `AB_EFF` far more than it lowers the titre, because it also blocks Fcγ
receptors, competes for FcRn and scavenges complement.

Because the untreated titre sits on the Hill plateau, partial suppression of the
source buys very little. This is where the model earns its keep:

| Maintenance therapy | Titre vs none | **ARR vs none** |
|---|---|---|
| Maintenance IVIG 1 g/kg q4wk | 0.76 | **0.28** |
| Rituximab | **0.70** | **0.75** |
| IL-6R blockade | 0.59 | 0.35 |
| Mycophenolate / azathioprine | 0.79 | 0.77 |
| FcRn inhibitor | 0.47 | 0.39 |
| C5 inhibitor | **1.00** | **0.75** |
| IVIG + IL-6R blockade | 0.44 | 0.04 |

Rituximab lowers the titre **more** than IVIG and prevents **far fewer**
relapses. The C5 inhibitor lowers the titre not at all and still helps. Titre
reduction and relapse protection are not the same axis, and which node a drug
acts on decides the outcome.

**Falsifiable consequences.** Serial MOG-IgG titres should be a poor surrogate
for treatment benefit. IVIG's benefit should persist in patients whose titres do
not fall. Both are testable in cohorts that already exist.

## 3. What the model reproduces

Annualised relapse rate, 365-day horizon, index optic neuritis, IV
methylprednisolone from day 4 and a guideline-like 90-day oral wean:

| Therapy | Model ARR | Observed ARR (source) |
|---|---|---|
| None (steroid wean only) | **0.64** | 0.64, 95% CI 0.58–0.70 — Vilaseca 2026 |
| Maintenance IVIG 1 g/kg q4wk | **0.18** | 0.22, 95% CI 0.15–0.32 — Vilaseca 2026 |
| Maintenance IVIG 0.4 g/kg q8wk | **0.37** | dose-dependence, HR 3.31 — Chen 2022 |
| IL-6R blockade | **0.23** | 0.09, 95% CI 0.06–0.14 — Vilaseca 2026 |
| Rituximab | **0.48** | 0.59–0.63 — Chen 2020; Thakolwiboon 2021 |
| Mycophenolate | **0.49** | 0.67–0.84 — Chen 2020; Thakolwiboon 2021 |

`LAM0` was set from the 0.64 figure; the rest of the column is a consequence of
the structure rather than of per-drug fitting.

Attack course under standard treatment: visual acuity nadir logMAR 1.60,
recovering to 0.37 (~20/47) at one year; RNFL 98 → 75 µm; serum NfL 8 → peak
160 pg/mL; CSF peak 72 cells/µL; surviving oligodendrocytes 0.88. Total IgG at
steady state runs 15.0 g/L on maintenance IVIG, 5.7 g/L on rituximab
(hypogammaglobulinaemia) and 4.7 g/L on the FcRn inhibitor.

## 4. A negative result, reported as such

The model was built expecting taper duration to matter. It does not. At a matched
cumulative dose of 840 mg prednisone-equivalent:

| Taper length | Starting dose | P(relapse) 6 mo | P(relapse) 1 y | BMD at 1 y |
|---|---|---|---|---|
| 14 days | 120 mg/d | 0.308 | 0.500 | 0.996 |
| 28 days | 60 mg/d | 0.301 | 0.495 | 0.995 |
| 90 days | 18.7 mg/d | 0.287 | 0.484 | 0.993 |
| 180 days | 9.3 mg/d | 0.282 | 0.479 | 0.992 |

Stretching the wean from two weeks to six months moves the one-year relapse
probability from 0.500 to 0.479. Raising the *total* dose helps somewhat more
(2700 mg over 90 days → 0.46) at a measurable bone cost. The mechanistic reading:
**oral steroids buy time roughly in proportion to exposure and cannot deliver
remission** — which is an argument for starting a steroid-sparing agent early,
not for tapering slowly.

Steroid dependency itself *is* emergent. Nothing in the model says "relapse when
steroids stop". The hazard depends on antibody; antibody depends on a plasmablast
pool that glucocorticoids suppress reversibly; when the wean ends, plasmablasts
recover in a fortnight and the titre follows over a month, so the hazard climbs
back on its own.

## 5. Model contents

**Mechanistic map** — 181 nodes, 263 edges, 17 subgraph clusters: host
susceptibility and triggers · peripheral B lineage and the plasmablast source ·
antibody biology and FcRn handling · T-cell and cytokine amplification ·
blood–brain barrier transit · target antigen and oligodendrocyte biology ·
complement and phagocytic effectors · tissue outcome and repair · the AQP4
contrast arm · optic neuritis · myelitis · brain phenotypes (ADEM, FLAMES,
brainstem) · acute therapy · maintenance therapy · biomarkers · clinical
endpoints · treatment toxicity.

**mrgsolve model** — 34 ODE compartments, 119 parameters, 33 captured outputs,
17 prebuilt scenarios. Drug arms: IV methylprednisolone, oral prednisone with a
continuously tunable linear taper, plasma exchange, acute and maintenance IVIG,
rituximab, IL-6R blockade, antimetabolites, FcRn inhibition, C5 inhibition, plus
the AQP4 comparator switch and an intercurrent-infection trigger.

**Shiny dashboard** — 12 tabs, including a taper explorer that re-runs the
matched-dose sweep live and a side-by-side MOGAD/NMOSD comparator.

```r
library(mrgsolve)
mod <- mread_cache("mogad", project = ".", file = "mogad_mrgsolve_model.R")
# MOGAD_simulate_scenarios() and MOGAD_taper_sweep() are supplied as commented
# R code at the foot of the model file (so that mread() sees a clean spec) —
# copy that block into your session, then:
df <- MOGAD_simulate_scenarios(365)
shiny::runApp("mogad_shiny_app.R")
```

```bash
python3 mogad_reference_impl.py        # regenerates every number quoted above
dot -Tsvg mogad_qsp_model.dot -o mogad_qsp_model.svg
dot -Tpng -Gdpi=150 mogad_qsp_model.dot -o mogad_qsp_model.png
```

## 6. Known misses

Stated here rather than buried, because they are where the model is most likely
to be wrong.

1. **IL-6R blockade is under-predicted** — model 0.23 against a reported 0.09.
   Two readings, and the model cannot distinguish them: either the observed
   incidence-rate ratio of 0.08 is a pre-versus-on-treatment comparison inflated
   by regression to the mean (the trustworthy contemporaneous comparison, against
   IVIG ≥1 g/kg, was *not* significant — in which case the model is right), or
   the IL-6 axis carries more weight than `EMAX_IL6` and `KAMP` allow. Randomised
   data would settle it.
2. **`FRAC_ESC = 0.60` is fitted, not measured.** It is the single parameter that
   makes rituximab weaker than IVIG, so the whole source-versus-clearance
   conclusion rests on it. The competing explanation — that relapses on
   rituximab cluster in B-cell repopulation windows rather than arising from
   CD20-negative escape — corresponds to a low `FRAC_ESC` plus under-dosed
   redosing, and predicts that strict CD19-guided redosing rescues rituximab.
   That is the experiment which separates the two.
3. **Azathioprine and mycophenolate share one arm**, so the reported ordering
   between them (azathioprine apparently better) cannot be reproduced. That
   ordering is mechanistically unexplained and may be channelling in
   retrospective series; the model takes no position.
4. **Attack nadir is under-predicted and recovery is pessimistic** — nadir
   logMAR 1.60 against a reported average of count-fingers, final acuity 20/47
   against a reported average 20/30. Both would improve with an explicit
   conduction-block state separate from myelin content, which is not implemented.
5. **Relapses are a hazard, not events.** No attack actually fires from the
   hazard, so step-wise disability accrual across multiple attacks is not
   simulated. Drawing attack times from the hazard and injecting `TRIG` is the
   natural extension.
6. **No paediatric ADEM physiology.** Encephalopathy, seizures and cognitive
   outcome appear in the mechanistic map but have no ODE counterpart, and the
   much higher rate of transient seropositivity in children is not modelled.
7. **Titre-to-assay mapping is linear.** Real live cell-based assay dilutions are
   ordinal and step in twofold increments, so comparisons of simulated titres to
   reported dilutions are order-of-magnitude only.

## 7. Disclaimer

Educational and research model only. Parameters are approximations fitted to
published aggregate cohort data, not to individual patients; the model has not
been independently validated, and it must not be used for clinical decisions,
prescribing or regulatory submissions.
