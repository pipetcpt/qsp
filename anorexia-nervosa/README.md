# Anorexia Nervosa (AN) — QSP Model

> Integrated Quantitative Systems Pharmacology model of anorexia nervosa,
> linking eating-disorder cognition and compulsive restriction to a
> Forbes-partitioned **energy-balance / body-composition** core, the
> **starvation neuroendocrine cascade** (leptin · ghrelin · low-T3 ·
> hypercortisolemia · GnRH/estradiol · GH resistance with low IGF-1),
> **uncoupled bone remodeling** (P1NP down / CTX up → BMD Z-score),
> **cardiac adaptation** (bradycardia, QTc) and the hazard that nutritional
> rehabilitation itself creates — the **refeeding syndrome** — with PK/PD for
> olanzapine, fluoxetine (nutritionally gated), transdermal estradiol,
> teriparatide, risedronate and psychotherapy.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`an_qsp_model.dot`](an_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`an_qsp_model.svg`](an_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`an_qsp_model.png`](an_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`an_mrgsolve_model.R`](an_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`an_shiny_app.R`](an_shiny_app.R) |
| 📚 References (71)        | [`an_references.md`](an_references.md) |

---

## 1. Disease in one paragraph

Anorexia nervosa is the psychiatric illness with the highest mortality
(standardized mortality ratio ≈ 5.9). Cognitively driven restriction — fear of
weight gain, body-image distortion, habit-entrenched food avoidance and
compulsive exercise — holds intake below requirement, and the resulting
**negative energy balance** triggers a coordinated, largely *adaptive*
starvation response: fat and lean mass fall, resting energy expenditure is
suppressed *beyond* what body composition predicts (adaptive thermogenesis),
leptin collapses, free T3 falls with a normal TSH, cortisol rises, hepatic GH
resistance drops IGF-1, and the GnRH pulse generator switches off, producing
functional hypothalamic amenorrhea. Two of these arms converge on the skeleton:
low IGF-1 and low leptin suppress **bone formation** while estrogen deficiency
and hypercortisolemia accelerate **resorption**, giving AN its characteristic
*low-formation/high-resorption* uncoupling — unlike post-menopausal
high-turnover loss — with BMD Z-score decline, failure to accrue peak bone mass
and fragility fracture. Low T3 and vagal predominance produce bradycardia, and
electrolyte loss (purging) lengthens QTc. Critically, **starvation itself
perpetuates the cognition** (the Minnesota Semi-Starvation Experiment
reproduced AN-like rigidity and food preoccupation in healthy volunteers), so
nutritional restoration is a prerequisite for psychological recovery rather
than a consequence of it. But refeeding carries its own hazard: an insulin
surge drives phosphate, potassium and magnesium intracellularly onto an occult
total-body deficit and consumes thiamine — the **refeeding syndrome**.

## 2. Mechanistic map — 21 clusters, 200 nodes, 279 edges

1. Predisposition & risk architecture (polygenic liability, metabo-psychiatric GWAS loci, female sex, pubertal estrogen surge, perfectionism, harm avoidance, thin-ideal exposure, dieting trigger, microbiome, stressors)
2. Cognitive-affective neurocircuitry (ventral/dorsal striatum, restriction-as-habit, reward prediction error, food-cue anticipatory anxiety, insula interoception, body-image distortion, amygdala, dlPFC, OFC, set-shifting, **ED cognition drive**, **starvation perpetuation**)
3. Neurotransmitter & neurotrophic systems (5-HT1A/5-HT2A, tryptophan availability, dopamine D2/D3, noradrenergic tone, GABA/glutamate, endogenous opioids, BDNF)
4. Hypothalamic energy-balance circuitry (ARC NPY/AgRP vs POMC/CART, α-MSH, MC4R, PVN, LHA orexin, AMPK/mTOR, LepRb-STAT3, GHSR-1a, kisspeptin, starvation-induced hyperactivity)
5. Gut-brain peptides & GI physiology (acyl-ghrelin, PYY 3-36, GLP-1, CCK, insulin, gastroparesis, early satiety, constipation, SMA syndrome, microbiome/SCFA)
6. Energy balance & body composition (intake, meal plan, restriction, purging, exercise, REE, adaptive thermogenesis, TEF, TEE, FM, FFM, BMI, %mBMI)
7. Adipose endocrine axis (adipocyte mass, leptin, adiponectin, FFA, ketones, brown fat, hypothermia, lanugo)
8. HPA axis (CRH, ACTH, cortisol, urinary free cortisol, GR feedback resistance, catabolic drive)
9. HPG axis & functional hypothalamic amenorrhea (GnRH pulsatility, LH/FSH, follicular arrest, estradiol, progesterone, amenorrhea, male testosterone, subfertility)
10. GH/IGF-1 resistance axis (GH secretion, hepatic GHR-STAT5, IGF-1, IGFBP-1/3, growth arrest, FGF21)
11. Thyroid — low-T3 (euthyroid sick) syndrome (TRH/TSH, T4, deiodinase D1/D3 shift, free T3, reverse T3)
12. Bone remodeling & skeletal outcomes (osteoblast/osteoclast, RANKL/OPG, sclerostin, MSC-to-adipocyte switch, marrow adiposity, P1NP, CTX, BMD Z-score, microarchitecture, peak bone mass, fracture)
13. Cardiovascular & autonomic adaptation (vagal predominance, sympathetic withdrawal, bradycardia, orthostasis, LV mass loss, MVP, pericardial effusion, QTc, arrhythmia, refeeding cardiac failure)
14. Hematologic · hepatic · immune · renal (gelatinous marrow, cytopenias, starvation hepatitis, refeeding steatosis, hypoglycemia, immune dysfunction, pseudo-Bartter)
15. Electrolytes, refeeding syndrome & micronutrients (occult phosphate depletion, carbohydrate load, insulin surge, intracellular shift, phosphate/potassium/magnesium, ATP & 2,3-DPG, thiamine, Wernicke risk, edema, repletion)
16. Purging-subtype complications (vomiting, laxative abuse, metabolic alkalosis, hypokalemia, dental erosion, esophageal injury, Russell's sign)
17. Drug PK/PD — olanzapine (absorption, CYP1A2/UGT1A4, H1/5-HT2C appetite, D2 occupancy, metabolic AEs)
18. Drug PK/PD — fluoxetine (norfluoxetine, SERT occupancy, **nutritional efficacy gate**, relapse prevention)
19. Drug PK/PD — bone-targeted therapy (transdermal E2 vs oral contraceptive first pass, teriparatide, risedronate, vitamin D/calcium)
20. Psychological & nutritional interventions (FBT, CBT-E, MANTRA/SSCM, engagement, admission, higher-calorie vs conservative refeeding, post-discharge relapse)
21. Clinical endpoints (BMI/%mBMI, EDE-Q, menses resumption, BMD Z-score, HR/medical stability, phosphate nadir, weight restoration, length of stay, relapse, SMR, quality of life)

```bash
dot -Tsvg an_qsp_model.dot -o an_qsp_model.svg
dot -Tpng -Gdpi=150 an_qsp_model.dot -o an_qsp_model.png
```

## 3. mrgsolve model — 31 ODE compartments, 147 parameters

| Block | Compartments |
|---|---|
| Drug / intervention PK (9) | `OLZ_GUT`, `OLZ_CENT`, `FLX_CENT`, `FLX_NOR`, `E2_PATCH`, `E2_EXO`, `TPTD`, `BIS`, `THERAPY` |
| Behavior & energetics (6) | `DRIVE`, `EI`, `EI_ADAPT`, `FM`, `FFM`, `ADAPT` |
| Neuroendocrine axes (7) | `LEP`, `GHRL`, `T3`, `CORT`, `GNRH`, `E2`, `IGF1` |
| Skeletal (3) | `P1NP`, `CTX`, `BMD` |
| Cardiac (1) | `HR` (QTc derived algebraically) |
| Refeeding safety (5) | `DEPL`, `PHOS`, `POT`, `MG`, `THIA` |

Three structural choices carry most of the model's clinical content:

* **Forbes partition closed against energy conservation.** The lean:fat *mass*
  ratio follows `dFFM/dFM = C/FM` (C = 10.4 kg) and is solved together with
  `EB = ρ_FM·dFM + ρ_FFM·dFFM`. A patient presenting with 3-4 kg of fat
  therefore puts most early refeeding weight into lean/fluid mass, and the
  split shifts toward fat as fat mass is restored.
* **One starvation index drives every axis.** `fdef` combines a slow
  body-store term (%mBMI deficit) with a fast energy-flow term (intake vs
  requirement). This is why T3, cortisol and IGF-1 normalize within days-weeks
  of refeeding while leptin, GnRH/estradiol and BMD follow over months —
  and it reproduces the 3-6 month lag between weight restoration and return
  of menses via a leptin gate at **1.85 ng/mL** feeding a slow pulse generator.
* **Refeeding risk is set by escalation rate, not calorie count.** The insulin
  surge is `(EI − EI_ADAPT)`, i.e. intake above what the tissues are adapted
  to, multiplied by the occult depletion index. This reproduces the modern
  evidence that higher *starting* calories are safe when escalation is
  controlled and phosphate is replaced, while steep escalation onto a large
  deficit produces the classic day 3-7 phosphate nadir.

### 10 prebuilt scenarios (simulated output, 365 d, presenting BMI 14.0, adolescent)

| # | Scenario | Weight gain (first 30 d) | Final BMI | PO₄ nadir | Anchor |
|---|---|---|---|---|---|
| 1 | Untreated natural history | −0.26 kg/wk | 11.7 | 3.74 | Keys 1950; Arcelus 2011 |
| 2 | Inpatient higher-calorie refeeding (2000 kcal start) | **+1.17 kg/wk** | 21.3 | 2.87 (d 9) | Golden 2013; Garber 2021 |
| 3 | Conservative slow refeeding (1000 kcal start) | +0.30 kg/wk | 21.3 | 3.49 (d 23) | O'Connor 2016 |
| 4 | Steep escalation, no PO₄ prophylaxis (BMI 12.5, 4 y illness) | +1.49 kg/wk | 21.3 | **2.22 (d 6.5)** | Whitelaw 2010; Friedli 2017 |
| 5 | Same escalation **+ PO₄/Mg/thiamine prophylaxis** | +1.49 kg/wk | 21.3 | 3.34 | ASPEN 2020 |
| 6 | Outpatient FBT (adolescent) | +0.40 kg/wk | 21.3 | 3.56 | Lock 2010 |
| 7 | Outpatient CBT-E (adult) | +0.37 kg/wk | 21.3 | 3.57 | Fairburn 2013 |
| 8 | Olanzapine 10 mg + outpatient therapy | +0.62 kg/wk | 23.9 | 3.41 | Attia 2019 |
| 9 | Fluoxetine 60 mg (nutritional gate) | +0.40 kg/wk | 21.3 | 3.56 | Attia 1998; Walsh 2006 |
| 10 | Transdermal E2 + discharge at day 90 | +1.30 kg/wk | 17.0 (relapsing) | 3.16 | Misra 2011 |

Three results worth reading off the simulation:

* **Scenario 4 vs 5** — identical calorie escalation; phosphate nadir 2.22 vs
  3.34 mg/dL and composite refeeding-risk index 0.63 vs 0.25. The hazard is
  the escalation, the mitigation is repletion, not caloric timidity.
* **Scenario 9** — the SSRI gate rises from **0.10 at presentation (68 %mBMI)
  to 0.76 once restored (104 %mBMI)** at unchanged plasma exposure: the same
  drug does nothing in the underweight patient and something modest afterwards.
* **Scenario 10** — discharge at day 90 (BMI 19.6) with therapy withdrawn:
  drive climbs back from 0.45 to 0.88, intake falls below requirement and BMI
  drifts to 15.7 by day 540 while BMD resumes its decline. The relapse is
  driven by the cognition term, not by the meal plan.

Bone-arm comparison (adolescent, 540 d, BMD Z-score from −1.60):

| Arm | Final BMD Z | Slope, last 90 d |
|---|---|---|
| Weight restoration only | −1.75 | +0.025 |
| + transdermal estradiol | −1.66 | +0.038 |
| + teriparatide | −1.66 | +0.042 |
| + risedronate | −1.28 | +0.097 |
| + **oral** contraceptive | −1.88 | — (hepatic first pass drops IGF-1 to 232 ng/mL) |
| Adult, no bone therapy | −2.03 | −0.032 |

## 4. Shiny dashboard — 9 tabs

1. **Patient profile** — presenting BMI, median BMI for age/height, severity, illness duration, subtype, adolescent/adult, repletion sliders; baseline state table and value boxes.
2. **Energy & body composition** — prescribed plan vs actual intake vs expenditure, net energy balance, BMI/%mBMI, FM/FFM, adaptive thermogenesis and compulsive exercise, ED drive.
3. **Neuroendocrine axes** — leptin & ghrelin, free T3 & cortisol, GnRH index & estradiol, IGF-1, menses-resumption probability.
4. **Drug PK/PD** — olanzapine and fluoxetine/norfluoxetine exposure, the SSRI nutritional gate, therapy-engagement effect.
5. **Refeeding safety** — phosphate (with the 2.5 mg/dL line), potassium & magnesium, insulin-surge index vs occult depletion, thiamine and the composite risk index, HR and QTc.
6. **Bone health** — P1NP/CTX, BMD trajectory, and a five-arm bone-therapy comparison run live on the current patient.
7. **Clinical endpoints** — day 0/30/90/180/end endpoint table, weight-restoration and remission flags, medical-instability indicator.
8. **Scenario comparison** — all 10 scenarios on the chosen patient profile with a summary table.
9. **References** — the bibliography and disclaimer.

## 5. How to run

```r
install.packages(c("mrgsolve", "dplyr", "tidyr", "ggplot2", "purrr",
                   "shiny", "shinydashboard", "DT"))

library(mrgsolve)
mod <- mread(model = "an", file = "an_mrgsolve_model.R", project = ".")

# inpatient higher-calorie refeeding with weekly therapy
out <- mod |>
  param(EI_START = 2000, EI_ESCAL = 200, EI_PLAN_MAX = 3400, INPATIENT = 1) |>
  mrgsim(events = ev(amt = 1.4, cmt = "THERAPY", ii = 7, addl = 51),
         end = 365, delta = 0.25)

plot(out, "BMI+Phosphate+BMD_Zscore+Menses_probability")
```

```r
# interactive dashboard (from this directory)
shiny::runApp("an_shiny_app.R")
```

The model has been compiled and executed with mrgsolve; the numbers in the
tables above are simulator output, not hand-written estimates.

## 6. Known simplifications

* Drug compartments are dose-proportional **exposure** surrogates mapped
  through Emax functions onto their mechanistic targets, not literal plasma
  PK models — the library's standard convention.
* `RHO_FFM_GAIN` (5200 kcal/kg) is an **effective** tissue cost, larger than
  the ~1800 kcal/kg composition-only figure, because refeeding weight gain
  also pays for synthesis work; the remainder of the gap is carried explicitly
  by the `K_HYPERMET` refeeding-hypermetabolism term.
* Binge/purge behavior is a single lumped term (`PURGE`, `PURGE_LOSS`,
  `K_PURGE_LOSS`) rather than separate vomiting/laxative/diuretic pathways.
* Mortality, hospitalization and length of stay appear as risk indicators, not
  as simulated time-to-event outcomes.
* The meal plan tapers to maintenance between `GOAL_LO` and `GOAL_HI` %mBMI,
  which is a modeling stand-in for clinician down-titration at goal weight.

## 7. Disclaimer

Educational and research use only. This model has not been validated for
clinical decision-making and **nothing in it is a refeeding protocol**.
Refeeding a severely malnourished patient is an inpatient medical procedure
requiring specialist eating-disorder care and mandatory electrolyte
monitoring. Parameters are illustrative approximations calibrated to published
group-level data, not to any individual patient.
