# Systems-Based Quantitative Modeling for Clinical Safety Assessment
### Five QSP/PBPK studies for the KSOT 2026 talk — *Mechanistic Insights and Quantitative Approaches for Drug Development and Safety Assessment*

Sungpil Han, MD PhD · Dept. of Pharmacology, College of Medicine, The Catholic University of Korea

This folder holds the five studies built for the talk. Each study is a
self-contained, re-runnable R script under `R/` that implements a QSP/PBPK
model, calibrates it to published clinical data, draws its figures with
ggplot2 into `figures/`, and writes its headline numbers to `data/`. All
figures share one visual system (`R/theme_tox.R`: S-Core Dream typeface to
match the slide body font, navy/teal/gold palette, English-only figure text),
and `build_deck.py` reads the numbers from `data/*.json` so the slide text and
the figures cannot drift apart.

Deck language rule: cover, titles, technical terms, tables and figures in
English; body text in Korean-English mix. Every paragraph is written with
`latinLnBrk="0"` (the template's "allow Latin words to break" is switched off)
so English words never split mid-word.

The organizing thesis is the abstract's: **QSP turns an efficacy–safety
trade-off into a tunable optimization problem**, and does so at the scale of a
living, LLM-built model library. The five studies span autoimmune/neurologic,
oncologic, metabolic, and geriatric disease, plus the library itself.

| # | Study | Disease area | The optimization | Key figures |
|---|-------|--------------|------------------|-------------|
| 1 | QSP library safety atlas | all (416 models) | coverage of a drug-safety axis across the therapeutic map | `s1_areas`, `s1_quality`, `s1_atlas` |
| 2 | Multiple sclerosis DMT benefit–risk | autoimmune / neurology | natalizumab dosing interval vs PML risk (EID) | `s2_saturation`, `s2_pml_table`, `s2_benefit_risk`, `s2_dmt_map` |
| 3 | Checkpoint-inhibitor efficacy vs irAE | oncology | nivolumab/ipilimumab dose ratio vs grade 3–4 irAE | `s3_pareto`, `s3_timecourse`, `s3_organ_timing` |
| 4 | Incretin titration | obesity / metabolic + aging | target dose & escalation vs GI AE & lean-mass loss | `s4_titration`, `s4_frontier`, `s4_composition`, `s4_aging` |
| 5 | DOAC PBPK in the elderly | aging / renal–hepatic impairment | age→renal function→exposure→bleeding; dose individualization | `s5_age_crcl`, `s5_validation`, `s5_exposure_response`, `s5_individualization` |
| ★ | Merigolix flagship | gynecologic | GnRH-antagonist dose vs estradiol-threshold window | `s6_merigolix` |

---

## Study 1 — The QSP library as a safety atlas

**Question.** If QSP is the translational bridge for clinical safety, how much
of the disease map does a modern, LLM-built QSP library already cover, and how
much of it already encodes a drug-toxicity / adverse-event axis?

**Method.** `R/s1_library.R` walks every disease directory in the
repository, counts ODE compartments (`dxdt_` assignments) in each mrgsolve
model, nodes/clusters in each Graphviz map, and PubMed references in each
reference file; classifies each model into a therapeutic area from the README
gallery table; and scores each model for safety/toxicology relevance by a
~50-term toxicology lexicon (toxic, adverse, hepatotox, nephrotox, QT,
cytokine release, irAE, bleeding, PML, DILI, …). Terms are word-bounded and
acronyms are matched case-sensitively (an earlier substring version counted
e.g. "making" as AKI and inflated the share to ~74 %). A model is counted as
carrying a "safety axis" when its model code, mechanistic map, directory
README and references together contain ≥3 distinct lexicon terms.

**Results (mined live from the repository).**

| Metric | Value |
|---|---|
| Disease QSP models | **416** |
| Therapeutic areas | 22 |
| Models carrying an explicit safety/toxicity axis | **71 %** (293 / 415 mined) |
| Total ODE compartments across the library | ~13,400 (median 25/model, range 14–82) |
| Total mechanistic-map nodes | median 165/model |
| Total PubMed references | ~31,800 (median 59/model; 96 % clear the ≥30 gate) |

**Message.** The library is not a collection of efficacy toys — about 70 %
of 416 models already carry a toxicity/adverse-event axis, and every model
clears hard size and grounding gates. The safety infrastructure for
NAM-translation already exists, at scale, and grows by one model per day
(The talk deliberately shows no count-over-time chart.)

---

## Study 2 — Multiple sclerosis: natalizumab dosing interval vs PML

**Question.** Natalizumab (anti-α4-integrin) is among the most effective
relapse suppressors in MS but lifts CNS immune surveillance, risking
progressive multifocal leukoencephalopathy (PML) from JC-virus reactivation.
Can QSP convert this into a dosing optimization?

**Method.** `R/s2_ms.R` couples (i) a receptor-saturation model of
natalizumab over the dosing interval (monoexponential drug decline, t½ ≈ 11 d,
Emax occupancy), (ii) relapse-rate suppression as a function of time-averaged
α4 saturation, and (iii) a PML-hazard model built from the Ho (2017)
risk-stratification table (anti-JCV index, prior immunosuppressant, treatment
duration), with extended-interval dosing (EID, q6w) applying a hazard ratio of
0.24 (Ryerson 2019).

**Calibration & results.**

- Efficacy calibrated to AFFIRM: q4w relapse-rate reduction **68.5 %**
  (target 68 %); q6w EID retains **61.7 %** — a cost of only **0.30 relapse
  over 6 years**.
- PML risk, highest-risk stratum (JCV+ index >1.5, 6 years): **14 → 3.4 per
  1000** with EID; prior-immunosuppressant stratum **13.2 → 3.2** — a ~76 %
  reduction.
- Cross-DMT benefit–risk map places natalizumab-q6w near ocrelizumab on
  efficacy at a fraction of the serious-safety-event index.

**Message.** The dosing interval is a single controllable handle that slides
each patient left along the risk axis at almost no efficacy cost — a textbook
"efficacy–safety trade-off turned into optimization."

Anchors: Polman NEJM 2006 (AFFIRM); Foley Lancet Neurol 2022 (NOVA); Ho Lancet
Neurol 2017; Bloomgren NEJM 2012; Ryerson Neurology 2019; Kappos NEJM 2010
(FREEDOMS); Hauser NEJM 2017 (OPERA).

---

## Study 3 — Checkpoint inhibitors: efficacy vs immune-related toxicity

**Question.** Nivolumab + ipilimumab is the most effective and most toxic
melanoma regimen. The two antibodies are different handles — can a model show
why lowering ipilimumab (CheckMate 511) preserves efficacy while cutting
toxicity?

**Method.** `R/s3_ici.R` maps receptor occupancy (PD-1, CTLA-4) to two
effector pools — anti-tumor and self-reactive — with Bliss-type synergy, then
to ORR (Emax/Hill) and grade 3–4 irAE rate. The central biology: ipilimumab's
contribution to anti-tumor **priming saturates at a low dose**, while its
contribution to **toxicity keeps climbing super-linearly** — so efficacy and
toxicity read CTLA-4 occupancy through two different curves. Also an ODE
time-course (tumor + two effector pools) and an organ-by-organ irAE-onset
panel (Weber 2017).

**Calibration & results.** Monotherapy arms from CheckMate 067, combination
arms from the head-to-head CheckMate 511. The model reproduces the headline:
shifting NIVO1/IPI3 → NIVO3/IPI1 drops grade 3–4 irAE by **~12 points** at
essentially unchanged ORR. The dose-ratio sweep yields a continuous
efficacy–safety frontier on which the four real regimens sit; the irAE
time-course shows response and toxicity on two separate clocks.

Anchors: Larkin NEJM 2015/2019 (CM067); Lebbé JCO 2019 (CM511); Weber JCO 2017.

---

## Study 4 — Incretin titration in obesity (with an aging subgroup)

**Question.** For semaglutide and tirzepatide, the target dose and escalation
schedule trade weight loss against GI tolerability and against lean-mass loss —
the last of which is a safety concern in older adults.

**Method.** `R/s4_obesity.R` is an energy-balance ODE: exposure
accumulates toward the steady state for the current dose; appetite suppression
= Hill(exposure) → caloric intake; body weight changes are partitioned into fat
and fat-free mass by a Forbes-type rule that is age-sensitive (anabolic
resistance); and GI adverse events are a nausea state driven by the rate and
level of occupancy, blunted by a tolerance state.

**Calibration & results.**

- Weight loss at 72 wk: semaglutide 2.4 mg **15.5 %** (STEP 1: 14.9 %);
  tirzepatide 5/10/15 mg **14.3 / 19.3 / 21.7 %** (SURMOUNT-1: 15 / 19.5 /
  20.9 %).
- Faster (2-weekly) escalation adds large early nausea peaks for no efficacy
  gain — titration is the handle.
- Body-composition cost rises with age: lean fraction of weight lost
  **36 % at age 45 → 47 % at age 70**; a muscle-preserving arm reclaims much of
  it. Lean-mass loss climbs monotonically with age across both drugs — a
  dose-reduction signal for the elderly / sarcopenic.

Anchors: Wilding NEJM 2021 (STEP 1); Jastreboff NEJM 2022 (SURMOUNT-1); Hall
energy-balance / Forbes composition.

---

## Study 5 — DOAC PBPK/PD in the elderly: dosing the excluded patient

**Question.** Phase 1 characterizes DOAC exposure in young, normal-renal
volunteers; the patients who take DOACs are old and renally impaired. Can a
PBPK/PD workflow predict exposure and individualize dosing for patients never
studied?

**Method.** `R/s5_doac.R` links age → creatinine clearance (population
decline) → DOAC exposure (Rowland renal-impairment scaling from each drug's
renal-clearance fraction) → bleeding and stroke risk (exposure–response), and
compares fixed vs label-based dose reduction in a 4000-patient virtual elderly
cohort.

**Calibration & results.**

- Exposure prediction validated against renal-impairment PK from drug labels
  (dabigatran, apixaban, rivaroxaban, edoxaban × mild/moderate/severe):
  **GMFE 1.18×** with renal scaling alone.
- **Mechanistic result:** renal scaling under-predicts dabigatran (moderate RI
  observed 3.2× vs 1.9× predicted) — its RI exposure rises more than GFR alone,
  requiring a P-gp/absorption term; adding that mechanism-based correction
  brings GMFE to **1.14×** (the same "selective, mechanism-based correction"
  discipline as the gut–liver–kidney MPS case).
- Population mean CrCl crosses the dose-reduction threshold (50 mL/min) at
  **age 79**.
- In the virtual elderly dabigatran cohort, fixed full dosing over-exposes
  **68 %** of patients above the exposure ceiling; the label dose-reduction
  rule cuts this to **24 %** and lowers modelled major bleeding **6.3 → 4.1 per
  100 patient-years** — dosing an excluded population without a trial in it.

Anchors: US product labels (dabigatran, apixaban, rivaroxaban, edoxaban RI PK);
Cockcroft-Gault; Rowland–Tozer renal-impairment adjustment.

---

## Flagship — Merigolix (the abstract's named example)

`R/s6_merigolix.R`: GnRH-receptor antagonist for endometriosis. The figure
uses the outputs of the merigolix QSP model (24 states: 2-compartment PK +
HPG-axis indirect-response PD + menstrual cycle, lesion and BMD; virtual
population at week 12; github.com/pipetcpt/merigolix-qsp), built on the Phase 1
SAD/MAD PK/PD model (n = 51). Week-12 E2 by dose, the share of patients inside
the 20–40 pg/mL estrogen-threshold window versus below 20 pg/mL (full
suppression → BMD / vasomotor risk), and the composite therapeutic index, which
peaks at 160 mg QD (78.3).

---

## The NAMs bridge

The talk closes by connecting the library to new approach methodologies: the
gut–liver–kidney microphysiological system coupled to whole-body PBPK that
predicted human AUC for six drugs within 4-fold (GMFE 1.85) **without clinical
PK input**, extended to disease-state chips (DSS colitis, cisplatin injury) that
distinguish gut- from hepatic-driven exposure change. Study 5's DOAC validation
is deliberately built in the same honest, mechanism-based-correction style, so
the five new studies and the MPS/NAM work read as one program.

`R/s7_diagrams.R` draws the per-drug AUC GMFE figure from the validation
summary table (6 drugs × 4 observations; overall 1.85×; uniform vs selective
scaling rule 9.7× vs 1.85×) and the daily-build architecture diagram.

---

### Reproduce

```bash
cd presentation/2026-tox
Rscript R/run_all.R          # all studies: data/*.json + figures/*.png (ggplot2)
python3 build_deck.py        # assembles ../2026-tox.pptx on the uploaded template
```

R packages: ggplot2, dplyr, tidyr, purrr, readr, jsonlite, scales, patchwork,
ggrepel, deSolve, systemfonts, ragg. Deck assembly: python-pptx, Pillow, lxml.
