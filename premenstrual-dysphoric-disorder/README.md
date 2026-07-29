# Premenstrual Dysphoric Disorder (PMDD) — QSP Model

> A 60-compartment quantitative systems pharmacology model of PMDD posed **not**
> as "the luteal phase floods the brain with progesterone, so lower it" but as a
> **GAIN** problem: the neurosteroid-to-affect transduction is **non-monotonic**
> (an inverted U whose peak sits inside the physiological luteal range) and its
> gain is set by a **slow, symmetric change detector** in GABA-A receptor
> subunit composition. Hormone concentrations are *identical* in the model's
> PMDD and control phenotypes — only three parameters differ.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT, 222 nodes / 18 clusters / 302 edges) | [`pmdd_qsp_model.dot`](pmdd_qsp_model.dot) |
| 🖼️ Map (SVG) | [`pmdd_qsp_model.svg`](pmdd_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi) | [`pmdd_qsp_model.png`](pmdd_qsp_model.png) |
| ⚙️ mrgsolve ODE model (60 compartments, 30 scenarios) | [`pmdd_mrgsolve_model.R`](pmdd_mrgsolve_model.R) |
| 🐍 Dependency-free Python twin (38 self-checks, all passing) | [`pmdd_python_twin.py`](pmdd_python_twin.py) |
| 📊 Shiny dashboard (12 tabs) | [`pmdd_shiny_app.R`](pmdd_shiny_app.R) |
| 📚 References (59 PubMed-verified) | [`pmdd_references.md`](pmdd_references.md) |

---

## 1. The problem the model is built to solve

PMDD is a cyclical affective disorder confined to the luteal phase of the
ovulatory cycle, remitting within days of menses onset, affecting roughly
1.6–5.5% of reproductive-age women. Five findings about it are awkward
together, and any model of PMDD has to survive all five at once:

1. **Luteal E2, P4 and allopregnanolone (ALLO) concentrations in PMDD are
   indistinguishable from those of asymptomatic controls.** A level-based model
   therefore has no disease variable at all.
2. **Yet the ovarian steroids are causal.** Ovarian suppression abolishes the
   symptoms, and blinded add-back of estradiol *or* progesterone reinstates them
   in PMDD while doing nothing in controls (Schmidt 1998, [PMID 9435325](https://pubmed.ncbi.nlm.nih.gov/9435325/)).
3. **Lowering ALLO helps — but only a lot of lowering.** Dutasteride 2.5 mg
   blocks the luteal ALLO rise and improves symptoms; 0.5 mg does neither
   (Martinez 2016, [PMID 26272051](https://pubmed.ncbi.nlm.nih.gov/26272051/)).
4. **Blocking ALLO at its GABA-A site also helps.** Sepranolone
   (isoallopregnanolone), a neurosteroid-site antagonist, reduces symptoms
   (Bixo 2017, [PMID 28319848](https://pubmed.ncbi.nlm.nih.gov/28319848/)).
5. **And SSRIs, which RAISE brain ALLO**, are the best-evidenced treatment —
   working in 1–2 days rather than the 2–4 weeks they need in depression. They
   shift 3α-HSD toward the reductive direction (Griffin & Mellon 1999,
   [PMID 10441766](https://pubmed.ncbi.nlm.nih.gov/10441766/)).

A monotone "more ALLO is worse" model explains 3 and 4 and gets 5 backwards. A
monotone "more ALLO is better" model gets 3 and 4 backwards. Neither can
explain 1 and 2 at all, because in a monotone model the only way a normal
hormone level can be pathogenic is if the level is not normal.

## 2. The two structural commitments

**(1) The transduction is non-monotonic.** The affective drive produced by
effective neurosteroid load `L` at GABA-A is unimodal:

```
NS_DRIVE = SENS · (L/KP) · exp(1 − L/KP)          ← peaks at L = KP
```

Low-to-moderate ALLO is anxiogenic and dysphoric through extrasynaptic α4βδ
receptors; high load is sedative and anxiolytic, exactly as pregnanolone and
progesterone behave in humans (the bimodal mood association, Andréen 2006,
[PMID 16368261](https://pubmed.ncbi.nlm.nih.gov/16368261/)). **The disease
parameter is `KP` — where the peak sits relative to the physiological luteal
range — not the range itself.** In the PMDD phenotype `KP` = 5 nM, inside the
luteal ALLO range; in the control phenotype `KP` = 13 nM, above it.

**(2) The gain is a slow, symmetric change detector.** GABA-A subunit
composition moves on a timescale of **days** while ALLO turns over in
**minutes**:

* `DELTA` (δ-subunit abundance) follows the ALLO exposure *history* and is
  **tolerance**-like — its gain `GAMMA_D` is **negative**. Sustained exposure
  makes the system *less* reactive, as chronic-progesterone pseudopregnancy
  models show.
* `ALPHA4` (α4-subunit abundance) follows the **rate of change** of
  neurosteroid and total progestogen tone and is **sensitising** — its gain
  `GAMMA_A4` is positive. This is the withdrawal-induced α4 upregulation of
  Smith 1998 ([PMID 9572737](https://pubmed.ncbi.nlm.nih.gov/9572737/)).
* The detector is **symmetric**: the onset of exposure drives it as well as the
  withdrawal, which is why the add-back experiment behaves the way it does.

The whole difference between a PMDD phenotype and an asymptomatic control in
this model is three numbers: `KP` (5 vs 13 nM), `GAMMA_A4` (1.50 vs 0.25) and
`GAMMA_D` (−0.35 vs −0.05). Their hormones are byte-identical.

## 3. What the structure buys — five predictions and their falsification

Every number below is printed by `python3 pmdd_python_twin.py`. The
**falsified** column is the same model with `MONOTONE = 1`, `RATE_OFF = 1` and
both plasticity gains set to zero — a pure level detector, scaled so that it
reaches the *same* peak drive on an untreated cycle, so the comparison cannot be
won on scale.

| | Prediction | This model | Level detector |
|---|---|---|---|
| **[P1]** | the symptom peak **lags** the ALLO peak | ALLO peaks day 22.4, DRSP peaks day **28.0** → lag **+5.5 d** | lag **+1.0 d** — the peaks collapse together |
| **[P2]** | dutasteride works at 2.5 mg, not at 0.5 mg | **+45.4%** vs **+10.2%** | **+65.0%** vs **+36.1%** — an ordinary graded dose–response |
| **[P3]** | one neurosteroid analogue, two doses, **opposite signs** | sedative dose **+46.2%**, sub-sedative dose **−5.5%** | **−56.4%** and **−45.1%** — both harmful, monotonically |
| **[P4]** | a 4-day hormone-free interval beats a 7-day one | 24/4 **+59.2%** vs 21/7 **+50.8%** → **8.4-point** gap | gap **−1.0 point** — the regimens become equivalent |
| **[P5]** | add-back reinstates symptoms **transiently** | DRSP-11 20.1 → **30.8** in the first 10 days → **20.6** by day 112 of *identical* hormone input | 20.1 → 38.4 → **43.4** — it never fades |

Two of these deserve spelling out.

**[P2] is the flat top of the inverted U.** At `L = KP` the derivative
`d(NS_DRIVE)/dL` is exactly zero. Dutasteride 0.5 mg (31% SRD5A1 inhibition at
the simulated steady state of ~51 ng/mL) drops the luteal load from 5.9 to
4.1 nM and the position term barely moves — 0.986 to 0.983. Dutasteride 2.5 mg
(69% inhibition, ~255 ng/mL) drops it to 2.0 nM and falls off the peak onto the
low limb. **The dose threshold is not a potency threshold; it is a geometry
threshold.** A level model has no flat region and so predicts a graded response,
which is what the falsification run does.

**[P3] is the reason the SSRI story is not a contradiction.** In a
non-monotonic system the *sign* of ∂symptoms/∂ALLO depends on where you are on
the curve. Raising ALLO helps if you are above the peak and hurts if you are
below it. The model therefore predicts that the same neurosteroid analogue
should improve PMDD at a sedative dose and worsen it at a sub-sedative one —
which is also the cleanest available explanation of why brexanolone and
zuranolone treat a neurosteroid-withdrawal depression by giving *more*
neurosteroid ([PMID 29910040](https://pubmed.ncbi.nlm.nih.gov/29910040/),
[PMID 34190962](https://pubmed.ncbi.nlm.nih.gov/34190962/)). In this model the
SSRIs' own ALLO-raising arm (brain ALLO peak 5.15 → 6.58 nM) is close to
sign-neutral at the late-luteal symptom peak, and their benefit is carried by
the serotonergic arm; switching the transduction to a level detector makes the
ALLO rise unambiguously harmful and cuts the simulated sertraline effect from
**+56.5%** to **+32.9%**.

**[P5] is Schmidt 2017 in equations.** Symptoms are triggered by the *change*
in ovarian steroid levels and not by continuous stable levels
([PMID 28427285](https://pubmed.ncbi.nlm.nih.gov/28427285/)). A symmetric
change detector reproduces that directly: constant add-back on top of ovarian
suppression produces a flare while the signal is *moving* and then adapts away,
even though the hormone input never changes again. The same add-back leaves the
control phenotype asymptomatic (core luteal increase −1%).

## 4. The untreated cycle

| Quantity | Model | Target |
|---|---|---|
| E2 peak | 225.6 pg/mL | 180–320 |
| P4 peak | 10.9 ng/mL on cycle day 21.7 | 8–16, day 20–23 |
| plasma ALLO luteal peak | 5.57 nM | 4–8 |
| brain ALLO, follicular → luteal peak | 1.31 → 5.15 nM (peak day 22.4) | 1–2 → 4–8 |
| α4 / δ abundance at peak | 1.99× / 1.68× baseline | — |
| DRSP-11, follicular → luteal peak (PMDD) | 19.6 → 43.6 | — |
| DRSP-11, follicular → luteal peak (control) | 18.6 → 21.7 | — |
| core affective luteal increase, PMDD | **+134%** | DSM-5 requires ≥30% |
| core affective luteal increase, control | **+15%** | must stay <30% |

DRSP-11 is the sum of the eleven modelled DSM-5 domains, each scored 1–6, so 11
means no symptoms and 66 is the maximum. The DSM-5 ≥30% luteal-versus-follicular
criterion is evaluated on the five **core affective** domains (irritability,
depressed mood, anxiety, affective lability, anhedonia), which is where the
diagnostic weight actually sits.

## 5. Treatment scenarios

Effect = % reduction of the luteal DRSP-11 burden above the floor of the scale,
the model analogue of the primary endpoint of PMDD trials, on the fifth cycle.

| Regimen | Effect | Note |
|---|---|---|
| Sertraline 50 mg, luteal-only (cd 15–28) | **+56.5%** | **+49.9% in the very first treated cycle** — the fast onset that distinguishes PMDD from depression |
| Sertraline 50 mg, continuous | +62.6% | better on the depressive/anhedonic domains (core luteal 5.0 vs 5.9) via the slow BDNF branch |
| Sertraline 100 mg, continuous | +65.1% | flat dose–response above 50 mg |
| Sertraline 50 mg, symptom-onset (cd 22–) | +26.1% | starts too late to blunt the mid-luteal build-up |
| Fluoxetine 20 mg, continuous | +62.2% | long-t½ parent + norfluoxetine |
| Drospirenone/EE **24/4** | +59.2% | removes the corpus luteum and the endogenous ALLO surge |
| Drospirenone/EE **21/7** | +50.8% | the 7-day withdrawal window costs 8.4 points — **[P4]** |
| Drospirenone/EE continuous | +62.2% | no hormone-free interval at all |
| Leuprolide 3.75 mg q28d | +64.6% DRSP-11, **+80.0% core** | the residual is hypoestrogenic: hot flushes 2.13 → 5.54/day |
| Leuprolide + E2/P4 add-back | flare then remission | **[P5]**; protects the skeleton (+0.118 %/cycle vs −0.604 on leuprolide alone) |
| Sepranolone 10 mg SC q48h | +41.4% | competitive antagonism at the ALLO site (`KI_ISO` is fitted to this) |
| Dutasteride 2.5 mg/d | **+45.4%** | **[P2]** |
| Dutasteride 0.5 mg/d | **+10.2%** | **[P2]** — the flat top |
| Neurosteroid analogue, sedative dose | **+46.2%** | **[P3]** — descending limb |
| Neurosteroid analogue, sub-sedative dose | **−5.5%** | **[P3]** — ascending limb, symptoms worsen |
| Alprazolam 0.25 mg TID, luteal | +12.2% | limited because α4βδ receptors have no benzodiazepine site |
| Ulipristal acetate 5 mg/d | +68.6% | progesterone-receptor blockade sparing estradiol |
| Calcium carbonate 1200 mg/d | +28.8% | the adjunct with the best trial evidence |
| Spironolactone 100 mg/d, luteal | +3.8% | somatic domain only, as intended |
| CBT / mindfulness | +8.5% | prefrontal top-down training |
| Sertraline + drospirenone 24/4 | **+83.1%** | more than either alone — the two arms act on different nodes |
| Bilateral oophorectomy + transdermal E2 | +69.9% DRSP-11, **+83.4% core** | the cyclic affective burden is abolished; a constant estrogen-related somatic burden remains |

Placebo response is deliberately **not** modelled, so these numbers should be
compared with the *active-minus-placebo* differences of trials, not with raw
responder rates.

## 6. Mechanistic map — 18 clusters, 222 nodes

1. Predisposition — cellular sensitivity, not hormone level (ESC/E(Z), ESR1×COMT, SRD5A/AKR1C variants, GABRA4/GABRD, heritability)
2. HPO-axis cycle engine (KNDy → GnRH → FSH/LH → follicle → ovulation → corpus luteum → luteolysis → menses)
3. Ovarian steroids, binding proteins, and the two withdrawal windows
4. Neurosteroidogenesis: P4 → 5α-DHP → allopregnanolone, plus isoallopregnanolone, pregnanolone, THDOC, PregS/DHEAS
5. GABA-A receptor pharmacology — synaptic α1/γ2 versus extrasynaptic α4βδ, the neurosteroid site, the missing benzodiazepine site, KCC2/chloride
6. **GABA-A subunit plasticity** — the exposure memory, the withdrawal signal, and the time-constant mismatch that creates the lag
7. **The non-monotonic transduction** — ascending limb, flat top, descending limb, and where `KP` sits in cases versus controls
8. Serotonergic system — the fast (irritability) and slow (depressive) branches
9. Corticolimbic gain — amygdala, dACC, dlPFC, OFC, insula, ventral striatum
10. HPA axis and the blunted stress response of PMDD
11. Autonomic, sleep and circadian
12. Somatic symptom generation — RAAS/aldosterone/bloating, mastalgia, prostaglandins, craving, menstrual migraine, calcium-PTH
13. The eleven DSM-5 symptom domains and DRSP scoring
14. SSRI pharmacology — the two mechanisms and the three dosing strategies
15. Ovulation suppression — 24/4 versus 21/7, GnRH analogues, add-back, surgery
16. Neurosteroid-directed pharmacology — sepranolone, dutasteride 0.5 vs 2.5 mg, exogenous analogues, alprazolam
17. Adjunctive and non-pharmacological treatment
18. Clinical endpoints, safety, and the falsification test itself

## 7. The mrgsolve model — 60 compartments

| Block | Compartments |
|---|---|
| HPO cycle engine | FSH, LH, follicle, corpus luteum, E2, P4 |
| Neurosteroidogenesis | 5α-DHP, plasma ALLO, brain ALLO, brain isoALLO, exogenous analogue |
| GABA-A plasticity + detectors | ALLO memory, progestogen-tone memory, δ, α4, E2 memory |
| Serotonergic | extracellular 5-HT, autoreceptor desensitisation, BDNF/TrkB |
| Corticolimbic | amygdala, prefrontal control, reward |
| HPA axis | CRH, ACTH, cortisol |
| Somatic | aldosterone, ECF volume, breast tenderness, prostaglandin |
| Symptoms | 11 DRSP domains |
| Organ / safety | hot flushes, lumbar BMD, endometrial state |
| Drug PK | sertraline (+desmethyl), fluoxetine (+norfluoxetine), drospirenone, leuprolide, sepranolone, dutasteride, alprazolam, ulipristal |
| Pituitary | GnRH-receptor downregulation |

PK is parameterised from published values (sertraline t½ 26 h, fluoxetine 4 d /
norfluoxetine 9 d, drospirenone 31 h, dutasteride ~35 d, alprazolam 12 h), and
each drug attaches to a specific mechanistic node rather than to a generic
"effect" term: dutasteride to SRD5A1, sepranolone to the neurosteroid site,
drospirenone to gonadotropin suppression *and* the mineralocorticoid receptor
*and* the progestogen-tone detector, leuprolide through a flare-then-
downregulation pituitary state, SSRIs to SERT *and* to 3α-HSD.

**Two parameters are fitted rather than taken from the literature**: `KI_ISO`
(to the sepranolone effect size) and `W_SHT_AMY` (to the sertraline effect
size). Everything else is either published PK or is pinned by the physiological
ranges asserted in the twin's check table.

## 8. Reproducing the numbers

```bash
python3 pmdd_python_twin.py          # 30 scenarios, 38 checks — no dependencies
python3 pmdd_python_twin.py --csv out   # per-scenario time courses as CSV
dot -Tsvg pmdd_qsp_model.dot -o pmdd_qsp_model.svg
dot -Tpng -Gdpi=150 -Gsize=26,17 pmdd_qsp_model.dot -o pmdd_qsp_model.png
```

```r
library(mrgsolve); library(dplyr); library(shiny)
mod <- mread("pmdd_mrgsolve_model.R")
out <- PMDD_simulate_scenarios(mod)   # helper defined at the bottom of the file
PMDD_endpoint_table(out)
runApp("pmdd_shiny_app.R")            # 12-tab dashboard
```

The Python twin implements the identical 60-state system with a fixed-step RK4
integrator and the Python standard library only — no numpy, no scipy, no R — and
asserts 38 checks covering the cycle physiology, all five predictions, every
treatment effect quoted above and the safety readouts. **All 38 pass**; if a
parameter is edited, the failures say exactly which claim broke.

## 9. Shiny dashboard — 12 tabs

Patient profile · cycle engine · neurosteroids · **the inverted U with the
patient's own luteal trajectory drawn on it** · receptor plasticity · drug PK ·
the eleven DRSP domains · clinical endpoints · scenario comparison · **the five
predictions with a live falsification toggle** · safety and organ systems ·
model notes. The two presets (PMDD / control) change nothing but `KP` and the
two gains, which is the fastest way to see what the model is claiming.

## 10. Limitations

* An **educational and research** model. Not validated against patient-level
  data and not for clinical use.
* The cycle engine is a compact ODE oscillator tuned to reproduce the E2/P4/ALLO
  landmarks, not a full HPO model in the Röblitz sense
  ([PMID 23206386](https://pubmed.ncbi.nlm.nih.gov/23206386/)).
* Placebo response, adherence, comorbid depression and the perimenopausal
  transition are out of scope.
* `KP` is a phenotype knob, not a measurable quantity. Its testable content is
  the *pattern* of predictions [P1]–[P5], not its numerical value.
* The neurosteroid load `L` lumps allopregnanolone, pregnanolone, THDOC and
  exogenous analogues into one pool; the model does not distinguish their
  receptor-subtype selectivities.

---

**Disclaimer.** This model is for education and research only. It has not been
independently validated or certified, and must not be used for clinical
decisions, prescribing, or regulatory submission. Parameters and assumptions are
illustrative approximations that would require fitting and qualification against
real patient data before any applied use.
