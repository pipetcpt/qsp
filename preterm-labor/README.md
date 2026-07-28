# Preterm Labor / Spontaneous Preterm Birth (sPTB) — QSP Model

> An integrated Quantitative Systems Pharmacology model of spontaneous preterm
> birth, built around one structural claim: parturition is a **common terminal
> pathway with three effector limbs** — myometrial contraction, cervical
> ripening, membrane weakening — and **every tocolytic in clinical use acts on
> exactly one of them.** The delivery hazard is written as the *product* of the
> three, so the benefit of tocolysis has a hard arithmetic bound; the model
> prints that bound as a number. The interventions that change outcome
> (antenatal corticosteroids, magnesium, in-utero transfer) do not treat the
> uterus at all — they treat the fetus.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`ptl_qsp_model.dot`](ptl_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`ptl_qsp_model.svg`](ptl_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`ptl_qsp_model.png`](ptl_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`ptl_mrgsolve_model.R`](ptl_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`ptl_shiny_app.R`](ptl_shiny_app.R) |
| 📚 References (100)       | [`ptl_references.md`](ptl_references.md) |

---

## 1. The disease in one paragraph

Preterm birth — delivery before 37+0 weeks — complicates roughly one in ten
pregnancies worldwide and is the leading cause of death in children under five.
It is not a disease. It is a **syndrome**: several mechanistically distinct
processes, each capable of starting the same terminal cascade early. Four of
them dominate. (1) The **placental CRH clock** runs exponentially from
mid-gestation with a doubling time near 3.3 weeks, feeding fetal adrenal
DHEA-S, estriol and a positive cortisol feedback loop that is unique to the
placenta; women destined to deliver preterm are already on a steeper trajectory
by 20 weeks. (2) **Functional progesterone withdrawal** — in humans there is no
fall in plasma progesterone before labour; withdrawal is a shift in the
myometrial PR-A/PR-B ratio that de-represses the contraction-associated protein
genes and, critically, breaks the reciprocal trans-repression between the
progesterone receptor and NF-κB. (3) **Intra-amniotic infection and sterile
inflammation** — ascending *Ureaplasma* and dysbiotic organisms, but also
DAMP-driven sterile inflammation which is about twice as common as
culture-proven microbial invasion; either drives NF-κB, IL-1β/IL-6/IL-8,
COX-2 and MMPs, and can provoke a fetal inflammatory response syndrome that
independently predicts cerebral palsy. (4) **Decidual haemorrhage and uterine
overdistension** — thrombin acting through PAR-1 is directly uterotonic, and
mechanical stretch induces the same CAP programme through AP-1.

All four converge. The myometrium moves from quiescence through activation to
stimulation; the cervix softens, ripens and dilates; the membranes weaken and
rupture. Which is why the pharmacology fails the way it does.

## 2. What this model is built to show

**The hazard is multiplicative.**

```
lambda(t) = LAMBDA0 × f_contr × f_cerv × f_memb × f_term
```

Atosiban blocks OXTR. Nifedipine blocks Ca_v1.2. Indomethacin blocks COX.
Magnesium antagonises calcium at the myofilament. Four different molecular
targets — **all inside `f_contr`.** A drug that acts on one factor of a product
can only reduce that product by the amount its factor contributes, and
`PTL_limb_decomposition()` prints exactly that:

```
  Hazard limb decomposition at 24 h (GA 29 wk, CL 14 mm)
  arm            f_contr  f_cerv  f_memb   lambda  lambda_rel   MVU   P_48h
  untreated        2.285   7.691   1.024  0.04329       1.000  85.4   0.081
  atosiban         0.972   7.535   1.024  0.01803       0.416  32.9   0.034
  nifedipine       0.913   7.520   1.024  0.01691       0.391  30.5   0.028
  indomethacin     1.880   7.674   1.024  0.03553       0.821  69.2   0.068
  MgSO4 2 g/h      1.722   7.656   1.024  0.03248       0.750  62.9   0.062
  nifed + Mg       0.712   7.489   1.024  0.01314       0.303  22.5   0.022
  terbutaline      1.705   7.649   1.024  0.03213       0.742  62.2   0.062

  Arithmetic floor on lambda_rel with f_contr at its minimum: 0.066
```

Note the middle two columns. `f_cerv` is 7.5–7.7 in **every arm** and `f_memb`
is 1.02 in **every arm**. The drugs move one column. That is the whole argument.

And the consequence, from `PTL_ledger()`: atosiban cuts the probability of
delivering within 48 hours from 0.081 to 0.034 — it halves the 48-hour event,
exactly as the Cochrane reviews report — but adds only **0.58 days** to the
expected gestational age at delivery. Halving a hazard for two days, while the
cervical limb runs untouched for the following three weeks, is worth about half
a day. Forty years of tocolytic trials found the same thing.

**So what is tocolysis for?** For the 48 hours. `PTL_acs_window()`:

| Delay to delivery | Lung maturity (none → ACS) | RDS (none → ACS) | ARR | NNT |
|---|---|---|---|---|
| 0.5 d | 0.05 → 0.20 | 0.621 → 0.549 | 0.071 | 14.0 |
| 1 d   | 0.05 → 0.31 | 0.612 → 0.486 | 0.125 | 8.0 |
| 2 d   | 0.05 → 0.48 | 0.593 → 0.385 | 0.208 | 4.8 |
| **3 d** | **0.05 → 0.51** | **0.574 → 0.350** | **0.224** | **4.5** |
| 7 d   | 0.06 → 0.35 | 0.495 → 0.357 | 0.138 | 7.2 |
| 14 d  | 0.08 → 0.21 | 0.358 → 0.301 | 0.056 | 17.8 |

RDS relative risk at the peak is 0.61, against 0.66 in Roberts 2020 Cochrane.
The benefit rises over the first two days and has more than halved by day ten —
the empirical 24 h to 7 d window, emerging from the PK/PD rather than being
imposed. Tocolysis is worth exactly the width of that window, which is also why
maintenance tocolysis has never shown benefit.

**And when the trigger is infection, the days cost something.**
`PTL_infection_paradox()`, MIAC at 29 weeks:

| Arm | Days vs none | Fetal inflammatory index | Cumulative IL-6 exposure | RDS | CP | Composite |
|---|---|---|---|---|---|---|
| no tocolysis | — | 0.398 | 31.1 | 0.340 | 0.0604 | 0.594 |
| atosiban 48 h | +0.54 | 0.419 | 33.0 | 0.332 | 0.0612 | 0.592 |
| prolonged tocolysis | +2.09 | 0.470 | 38.7 | 0.323 | 0.0625 | 0.591 |

Prolonging by two days buys a 5% relative reduction in RDS and pays for it with
a 25% rise in cumulative fetal inflammatory exposure and a rise in cerebral
palsy and early-onset sepsis. The composite comes out **flat**. That is not the
dramatic reversal one might hope for from a model, and it is reported here
unmassaged — but a wash is already a damning result for an intervention whose
entire justification is that the days it buys improve outcome.

## 3. Mechanistic map — 22 clusters, 264 nodes, 408 edges

1. Risk architecture & aetiologic heterogeneity (prior sPTB, short cervix, multifetal, Müllerian anomaly, LEEP, short interpregnancy interval, smoking, periodontal disease, social determinants, EBF1/EEFSEC/AGTR2/WNT4 GWAS loci)
2. **Placental CRH clock** (syncytiotrophoblast CRH, CRH-BP, CRHR1/2, fetal ACTH, fetal-zone DHEA-S, 16α-OH-DHEAS, aromatase, estriol surge, fetal cortisol → CRH *positive* feedback)
3. **Progesterone quiescence & functional withdrawal** (PR-A/PR-B ratio, PR-C, 20α-HSD, miR-200/ZEB1-2, STAT5B, PGRMC1/mPRα, PR↔NF-κB mutual trans-repression)
4. Oestrogen activation arm (E1/E2/E3, ERα, GPER, CAP gene transcription)
5. **MIAC & intra-amniotic inflammation** (Lactobacillus→CST IV dysbiosis, ascending *Ureaplasma*/*Mycoplasma*, *Fusobacterium*, TLR2/4, NLRP3/caspase-1, IL-1β/IL-6/IL-8/TNF-α, neutrophil influx, chorioamnionitis, FIRS, fetal microglia; **sterile** DAMP-driven inflammation as a separate node)
6. Decidual activation, senescence & haemorrhage (p38-MAPK, SASP, exosomes, tissue factor, thrombin, PAR-1, TAT complexes)
7. Uterine overdistension (Laplace wall tension, Piezo1/TRPV4, integrin–FAK, AP-1, stretch-induced CAPs)
8. Prostaglandin biosynthesis & the **15-PGDH disposal barrier** (cPLA2, COX-1/2, mPGES-1, AKR1C3, EP1-4, FP, PGI2; barrier breach as an explicit node)
9. Myometrial activation — the CAP programme (phases 0–3, OXTR ↑100–200×, connexin-43, Ca_v1.2, BK-β1 loss, pacemaker currents)
10. Excitation–contraction coupling (V_m, Ca²⁺ influx, IP3R/RyR, STIM1/Orai1, CaM–MLCK, MLC20, MLCP, RhoA/ROCK sensitisation, Montevideo units)
11. Relaxation pathways (β2–Gs–AC–cAMP–PKA with GRK/β-arrestin tachyphylaxis, NO–sGC–cGMP–PKG, PDE, relaxin, Mg²⁺ antagonism)
12. Oxytocin system (maternal/fetal/decidual OT, oxytocinase, Gq–PLC–IP3, PKC, Ferguson reflex, V1a cross-talk)
13. **Cervical remodelling** (softening → ripening → dilation → repair; collagen cross-linking, decorin/lumican, HAS2/hyaluronan, MMP-1/-8/-2/-9, TIMPs, neutrophil elastase, sonographic length, funnelling, Bishop score)
14. **Fetal membranes & PPROM** (amnion/chorion, ZAM weak zone, GDF-15, apoptosis, MMP-9/TIMP-1, tensile strength, latency, oligohydramnios)
15. Biomarkers & bedside prediction (quantitative fFN and its ~99% NPV, PAMG-1, phIGFBP-1, amniotic IL-6/MMP-8, serial CL, cell-free RNA clock, QUiPP-style composite)
16. Progestogen pharmacology (vaginal micronized P4 with the **first-uterine-pass** effect, 17-OHPC flip-flop IM depot, PR-B/NF-κB/PGDH mechanism, PREGNANT–OPPTIMUM–EPPPIC vs MEIS–PROLONG)
17. Acute tocolysis (atosiban, nifedipine, indomethacin, MgSO4, β2-agonists, NTG; the 48-hour goal; no maintenance tocolysis)
18. Tocolytic & steroid adverse effects (pulmonary oedema, tachycardia, hypotension, the Mg toxicity ladder, ductal constriction, fetal oliguria, NEC signal, maternal and neonatal hypoglycaemia, repeat-course growth penalty)
19. **Antenatal corticosteroids** (betamethasone phosphate + acetate, transplacental transfer past 11β-HSD2, fetal GR, ATII maturation, SP-A/B/C + ABCA3, ENaC lung-liquid clearance, the 24 h–7 d window, ALPS, WHO ACTION)
20. Fetal neuroprotection & adjuncts (MgSO4 <32 wk, NMDA blockade, erythromycin/ORACLE-I, the co-amoxiclav NEC signal, GBS prophylaxis, cerclage, pessary, in-utero transfer)
21. Delivery timing & neonatal outcomes (hazard → GA at delivery → RDS, BPD, severe IVH, PVL, NEC, sepsis, ROP, death, cerebral palsy, NICU stay)
22. Model observables & decision layer (`PTL_ledger()`, limb decomposition, NNT comparison)

```bash
dot -Tsvg ptl_qsp_model.dot -o ptl_qsp_model.svg
dot -Tpng -Gdpi=150 ptl_qsp_model.dot -o ptl_qsp_model.png
```

## 4. mrgsolve model — 45 ODE compartments, 186 parameters

| Block | Compartments |
|---|---|
| Progestogen PK | `P4V_DEP`, `P4_C`, `OHPC_DEP`, `OHPC_C` |
| Acute tocolytic PK | `ATO_C`, `ATO_P`, `NIF_GUT`, `NIF_C`, `IND_GUT`, `IND_C`, `MG_C`, `MG_P` |
| Fetal-directed PK | `BET_DP`, `BET_DA`, `BET_C`, `BET_F`, `TERB_DEP`, `TERB_C`, `DESENS` |
| Endocrine ignition | `CRH`, `DHEAS`, `CORTF`, `E3` |
| Progesterone withdrawal | `PRW` |
| Infection / inflammation | `BACT`, `NFKB`, `IL6`, `FIRS` |
| Uterotonin machinery | `COX2`, `PG`, `OXTR`, `CX43`, `CAMP` |
| Contraction | `CAI`, `CONTR` |
| Cervix | `COLL`, `HA`, `MMP`, `CLEN` |
| Membranes | `MEMB`, `FFN` |
| Fetal maturation | `SURF`, `MGBRAIN` |
| Outcome accounting | `CUMHAZ`, `CUMINFL` |

Three modelling decisions are worth calling out because they carry most of the
behaviour:

- **Withdrawal is a state, not a concentration.** `PRW` (the PR-A/PR-B index)
  is an ODE compartment because in humans plasma progesterone does not fall
  before labour. A model that derived withdrawal from plasma P4 would be
  modelling a species we are not.
- **Atosiban is competitive, not Emax.** It enters as `1/(1 + C/Ki)` on the
  oxytocin arm. Raise `OT_TONE` in the app and watch the blockade erode — an
  Emax form cannot express that, and it is what clinicians observe.
- **Receptor density is not signal.** OXTR rises 100-fold at term but the Gq
  response saturates, so the oxytocin arm enters as a fractional occupancy
  `OXTR/(KM_OXTR + OXTR)`. Getting this wrong (using raw fold-change) drives
  the contraction index to its ceiling in every scenario and destroys all drug
  effect — it was the first thing numerical validation caught.

Magnesium clearance is scaled **linearly by creatinine clearance**, because
magnesium is cleared renally and by no other route. Set `CRCL = 40` and a
standard 4 g + 2 g/h regimen peaks at **12.3 mg/dL** instead of 4.2 — through
the respiratory-depression rung of the toxicity ladder. That is the single most
clinically important PK fact about the drug and the model is built to make you
trip over it.

### Twelve prebuilt scenarios

| # | Scenario |
|---|---|
| 1 | Low-risk singleton, natural history to term |
| 2 | Asymptomatic short cervix (18 mm at 20 wk), untreated |
| 3 | Short cervix + vaginal progesterone 200 mg nightly |
| 4 | Prior sPTB + 17-OHPC 250 mg IM weekly |
| 4b | Prior sPTB, no progestogen (comparator for 4) |
| 5 | Short cervix: cerclage + vaginal progesterone |
| 6 | Threatened preterm labour at 29 wk, untreated |
| 7 | Threatened PTL 29 wk: atosiban 48 h + betamethasone |
| 8 | Threatened PTL 29 wk: nifedipine + betamethasone + MgSO4 |
| 9 | Threatened PTL 27 wk: indomethacin + betamethasone + MgSO4 |
| 10 | MIAC with intra-amniotic inflammation at 29 wk |
| 11 | PPROM at 30 wk: erythromycin + betamethasone, expectant |

## 5. Verification

There is no R toolchain in the build environment, so the ODE system was
verified by **independent numerical integration** of the same equations and the
same parsed parameter block (LSODA, rtol 1e-7). Every number quoted in this
README is an output of that run, not an aspiration.

| Check | Model | Literature |
|---|---|---|
| Atosiban plateau, 300 µg/min phase | 408 ng/mL | ~430 ng/mL (Goodwin 1995) |
| Nifedipine C_max, 20 mg q6h | 90 ng/mL | 50–100 ng/mL |
| Indomethacin C_max | 1606 ng/mL | 500–3000 ng/mL |
| Serum Mg, 4 g + 1 g/h | 4.2 mg/dL | 4–7 mg/dL therapeutic |
| Serum Mg, same regimen at CrCl 40 | 12.3 mg/dL | respiratory-depression rung |
| Betamethasone maternal C_max | 105 ng/mL | 60–120 ng/mL |
| Vaginal P4 plasma / local | 12 / 345 ng/mL-eq | plasma 10–20; local ≫ plasma |
| 17-OHPC steady state | 15 ng/mL | 10–20 ng/mL (Caritis 2011) |
| Low-risk P(birth < 37 wk) | 0.073 | ~0.06–0.08 population sPTB |
| Tocolysis effect on P(delivery ≤ 48 h) | 0.081 → 0.034 | roughly halved (Cochrane) |
| ACS effect on RDS at peak | RR 0.61 | RR 0.66 (Roberts 2020) |
| MgSO4 effect on cerebral palsy | RR 0.735 | RR ~0.68 (BEAM 2008) |
| Erythromycin effect on PPROM latency | P(≤7 d) 0.322 → 0.271 | prolonged (ORACLE-I) |

**Where it misses, stated plainly.** Scenario 2 gives P(birth < 34 wk) = 0.575
for an untreated 18 mm cervix, against ~0.34 in the Fonseca 2007 placebo arm
and ~0.16 in PREGNANT. The model treats a short mid-trimester cervix as more
deterministic than it is — there is no random-effects structure and no
non-progressor subgroup, so every short cervix in the model progresses. The
knock-on is that vaginal progesterone comes out at RR 0.44 for birth < 34 wk,
near Fonseca (0.56) but stronger than the EPPPIC individual-participant
estimate (~0.78). **Treat the prevention arm as qualitatively, not
quantitatively, calibrated, and do not quote its absolute rates.** The acute
arm, the corticosteroid window and the tocolytic head-to-head are the parts
that reproduce their literature anchors.

## 6. Shiny dashboard — 10 tabs

Patient & pregnancy · Drug PK · Uterine PD · Cervix & membranes · **Hazard
decomposition** · Clinical endpoints · Scenario comparison · Biomarkers ·
Safety ledger · Model notes.

```r
library(shiny); library(mrgsolve)
shiny::runApp("ptl_shiny_app.R")
```

Three things to try: set `CRCL` to 40 and watch the magnesium panel cross the
toxicity rungs; raise oxytocin tone above 2.5 and watch atosiban's competitive
blockade erode; move gestational age past 32 weeks with indomethacin selected
and watch the ductal-constriction index climb.

## 7. Using the model directly

```r
library(mrgsolve)
mod <- mread("ptl_mrgsolve_model.R")
list2env(as.list(mrgsolve::env_get(mod)), .GlobalEnv)   # load the helpers

PTL_ledger(mod)                  # days gained vs OUTCOME gained
PTL_limb_decomposition(mod)      # the arithmetic bound on tocolysis
PTL_acs_window(mod)              # the 24 h - 7 d steroid window
PTL_infection_paradox(mod)       # days up, outcome flat-to-worse
PTL_tocolytic_head_to_head(mod)  # efficacy and harm in one table
```

## 8. The assumption to argue with first

The multiplicative hazard is a **modelling choice, not a measured fact**. It
was chosen because it is the simplest structure in which a single-limb
intervention has a hard arithmetic bound, and because it reproduces the
observed dissociation between "tocolysis delays birth by 48 hours" (true) and
"tocolysis improves neonatal outcome" (false, unless steroids, magnesium or
transfer are delivered in the interval). An **additive** hazard would let a
sufficiently strong tocolytic drive lambda to zero — precisely the prediction
forty years of trials have falsified. Anyone re-using this model should treat
that structural choice, rather than any parameter value, as its principal
assumption.

Other limitations: contractions are a smooth index, not discrete events, with
no electromyography layer; fast states (Ca²⁺, cAMP) use large rate constants
rather than their true second-to-minute timescale; normalised indices are
relative-scale QSP states, not assayable concentrations; the neonatal outcome
layer is a set of gestational-age logistics with additive shifts, not a
mechanistic neonatal model; and the MEIS/PROLONG disagreement over 17-OHPC is
*encoded* (through a lower effective local exposure with no first-uterine pass),
not resolved — do not read scenario 4 as evidence either way.

---

## ⚠️ Disclaimer

Research, education and hypothesis generation only. This is **not** a tocolysis
protocol, **not** a dosing calculator, and has **not** been validated against
patient data. Threatened preterm labour is an obstetric emergency requiring
clinical judgment. Parameters are interpretations of published data for
simulation purposes.

*Part of the [QSP Disease Model Library](../README.md).*
