# Familial Chylomicronemia Syndrome (FCS) — QSP Model

> An integrated Quantitative Systems Pharmacology model of familial
> chylomicronemia syndrome, built around one structural claim: plasma
> triglyceride is not set by a rate, it is set by a **saturated capacity**.
> Triglyceride-rich lipoprotein clearance is the **sum** of two limbs — a
> first-order, high-capacity lipolytic limb and a small, saturable,
> LPL-independent hepatic limb — and the FCS genotype multiplies the first one
> by **zero**. Everything clinically strange about this disease is a
> consequence of that single arithmetic fact: why a 20 g/day fat threshold
> behaves like a cliff rather than a slope, why forty years of fibrates and
> fish oil did nothing, why an antisense drug lowers triglyceride by 80% in a
> patient with undetectable lipase, and why a moderate percent reduction in
> mean triglyceride buys a much larger reduction in pancreatitis.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`fcs_qsp_model.dot`](fcs_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`fcs_qsp_model.svg`](fcs_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`fcs_qsp_model.png`](fcs_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`fcs_mrgsolve_model.R`](fcs_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`fcs_shiny_app.R`](fcs_shiny_app.R) |
| 📚 References (48)        | [`fcs_references.md`](fcs_references.md) |

**Scale:** 21 clusters · 222 nodes · 290 edges · 43 ODE compartments ·
5 drug PK models · 10 scenarios · 8 analysis functions.

---

## 1. The disease in one paragraph

Familial chylomicronemia syndrome is a monogenic, autosomal-recessive failure
of intravascular lipolysis affecting roughly 1–10 people per million. Biallelic
loss-of-function variants in any one of five genes — *LPL* itself (60–80% of
cases), *APOC2* (the obligate activator), *APOA5* (the tethering factor),
*GPIHBP1* (the endothelial platform that carries lipase to the capillary
lumen) or *LMF1* (the ER chaperone that folds it) — produce the same
phenotype, because they are five steps in the assembly of one machine. Fasting
triglyceride runs at 1 000–10 000 mg/dL from infancy, the plasma is visibly
milky, and the defining complication is **recurrent acute pancreatitis**,
which is more severe than biliary pancreatitis and, over a lifetime of
episodes, converts into chronic pancreatitis, exocrine failure and type 3c
diabetes. Between attacks patients carry a burden that trials have historically
not measured: daily abdominal pain, eruptive xanthomas, hepatosplenomegaly,
cognitive "brain fog", and a fat-restricted diet of under 20 g per day that
dominates every meal, holiday and social event of their lives. The disease is
also chronically **misdiagnosed** — as alcoholic or biliary pancreatitis, with
a median diagnostic delay measured in years — partly because lipaemic serum
interferes with the assays used to look for it, including amylase.

## 2. Why this model is built as a sum and not as a rate

Almost every published account of chylomicronemia describes the problem as
"impaired lipolysis". That description is true and, quantitatively, almost
useless: it does not predict a single one of the disease's clinical oddities.
This model instead writes clearance as an explicit sum of two structurally
different terms:

```
dTG/dt =  INPUT(diet, hepatic VLDL)
        - CL_LPL · f_geno · C              limb 1: first-order, capacity ~140 dL/h
        - Vmax_ind · g(apoC-III) · C/(Km+C) limb 2: saturable, Vmax ~1300 mg/h
        - CL_res · C                        residual scavenging
```

In health, limb 1 carries more than 95% of the flux and limb 2 is invisible.
In FCS, `f_geno ≈ 0`, so the entire dietary fat load is forced through a
Michaelis–Menten route whose Vmax (~31 g of triglyceride per day) is *the same
order of magnitude as the daily fat load itself*. Four consequences follow, and
each has a function in the R model that prints it as a number.

### Claim 1 — the saturation cliff (`FCS_saturation_curve`)

Because limb 2 saturates, steady-state triglyceride is a **hyperbola** in
dietary fat, not a line. In the calibrated LPL-null model:

| Dietary long-chain fat | Steady-state TG |
|---|---|
| 10 g/day | ~600 mg/dL (below the pancreatitis threshold) |
| 20 g/day | ~1 450 mg/dL |
| 60 g/day | ~9 500 mg/dL |

Tripling the fat intake multiplies triglyceride by more than six. This is why
the dietary prescription behaves as an **edge** rather than a dose-response,
why adherence in FCS is effectively all-or-nothing, and why a single restaurant
meal can be clinically consequential in a way that would be absurd in ordinary
hypertriglyceridemia. `FCS_binge_cost()` prices exactly that meal.

### Claim 2 — every conventional drug is arithmetically dead (`FCS_limb_decomposition`)

Fibrates work by inducing *LPL* transcription through PPAR-α. Insulin induces
LPL. Evinacumab works by removing ANGPTL3's inhibition *of LPL*. All three
multiply `f_geno`. Zero times anything is zero — so the model predicts their
failure in LPL-null FCS before any trial data are consulted, and
`FCS_limb_decomposition()` prints `limb1_LPL = 0.0 mg/h` in every such arm no
matter what is administered. The residual 10–20% triglyceride reduction seen
clinically with fibrate plus omega-3 is reproduced by the model, and it comes
entirely from the *input* side (omega-3 suppressing VLDL secretion) and from a
small apoC-III transcriptional effect — never from lipolysis.

The sharpest form of this prediction concerns evinacumab, and it has been
tested: in severe hypertriglyceridemia, the ANGPTL3 antibody lowered
triglyceride in patients **without** biallelic LPL-pathway mutations and did
little in patients **with** them (reference 22). `FCS_scenario_evinacumab()`
runs both arms side by side.

### Claim 3 — apoC-III knockdown works without any lipase at all

The mechanism that makes apoC-III antisense therapy effective in a patient with
zero lipase is that apoC-III does **two** jobs. It inhibits LPL — irrelevant
here — and it blocks hepatic remnant uptake via LDLR, LRP1 and syndecan-1
HSPG, a route that never needed lipase. In the model that second job is a
single multiplier on the Vmax of limb 2:

```
g(apoC-III) = 1 + IMAX_C3 · (1 − apoC-III / apoC-III_baseline)
```

With `IMAX_C3 = 2.2`, an 80% knockdown nearly triples the capacity of the only
working clearance route, and the model reproduces APPROACH (−77%), the Balance
open-label 12-month value (−73.7%) and PALISADE (−80%) — in a simulated patient
whose post-heparin lipase activity is identically zero.

### Claim 4 — the pancreatitis endpoint is convex, and that is where the benefit is

Acute pancreatitis hazard is modelled as a Hill function of triglyceride with
exponent 1.7 and half-maximal effect at 10 000 mg/dL, giving essentially no
hazard below the 880 mg/dL (10 mmol/L) threshold and a steeply rising hazard
above it. Two things follow.

First, because the hazard is **convex**, Jensen's inequality applies:
`E[λ(TG)] > λ(E[TG])`. A patient with no true fasting state, whose
triglyceride swings by thousands of mg/dL after every meal, carries more risk
than their fasting number implies. `FCS_jensen_gap()` quantifies the
under-report.

Second, convexity means that deleting the *upper tail* of the triglyceride
distribution removes disproportionately more risk than the change in the mean
suggests. That is the model's explanation for the otherwise puzzling
observation that apoC-III therapies reduce pancreatitis events by 83–88% —
numbers considerably larger than their headline percent triglyceride
reductions. `FCS_trial_ledger()` puts the two columns next to each other, and
argues that the endpoint these drugs should be scored on is **days per year
above 880 mg/dL**, not percent change in a fasting value that does not exist.

## 3. What is in the mechanistic map

`fcs_qsp_model.dot` renders 222 nodes across 21 clusters:

1. Genetic architecture — five genes, one machine, plus acquired autoantibody
   phenocopies and the *APOC3*/*ANGPTL3* loss-of-function natural experiments
2. Secondary amplifiers — alcohol, oestrogen, pregnancy, diabetes, drugs
3. Dietary fat as the input term, MCT portal bypass, adherence decay
4. Enterocyte chylomicron assembly (MTP, apoB48, Sar1b) and lymphatic transit
5. The LPL machinery — LMF1 folding, GPIHBP1 transcytosis, apoC-II activation,
   ANGPTL3/4/8 inhibition — **limb 1**
6. LPL-independent clearance — LDLR, LRP1, syndecan-1, RES — **limb 2**
7. apoC-III biology: transcription, both inhibitory jobs, VLDL promotion
8. The hepatic VLDL axis (the input the diet does not control)
9. The plasma compartment, lactescence, viscosity, assay interference
10. Acute pancreatitis: local lipolysis → unbound FFA → acinar Ca²⁺ overload →
    trypsinogen activation, plus the trypsin-independent ischaemic limb
11. Non-pancreatic burden: xanthomas, lipaemia retinalis, brain fog, QoL
12. Conventional TG drugs — all limb-1, all predicted to fail
13. Volanesorsen (naked 2′-MOE PS-ASO)
14. Olezarsen (GalNAc₃-ASO)
15. Plozasiran (GalNAc-siRNA, RISC effect compartment)
16. ANGPTL3 axis, gene therapy and the remaining pipeline
17. Safety architecture — why chemistry, not target, decided the class
18. Endpoints, and why percent TG change is the wrong one
19. Diagnosis: separating FCS from the 100× commoner multifactorial form
20. QSP model structure
21. The five structural claims, stated as claims

## 4. The ODE model

`fcs_mrgsolve_model.R` implements 43 compartments: gut/lymph fat transit,
four plasma lipoprotein-TG species, hepatic and adipose lipid pools, the
APOC3 mRNA→protein cascade, ANGPTL3, a functional LPL pool with turnover, the
pancreatic injury cascade, cumulative hazard and time-above-threshold
accumulators, four non-pancreatic burden scores, platelets and ALT, and full
PK for five drugs:

| Drug | Modality | PK structure | Regimen simulated |
|---|---|---|---|
| Volanesorsen | naked 2′-MOE PS-ASO | SC depot → plasma → liver + systemic tissue | 300 mg weekly |
| Olezarsen | GalNAc₃-ASO | as above, ASGPR-weighted to liver | 50 / 80 mg monthly |
| Plozasiran | GalNAc-siRNA | + RISC effect compartment (t½ ≈ 45 d) | 25 mg every 3 months |
| Evinacumab | anti-ANGPTL3 mAb | 2-compartment IV | 15 mg/kg every 4 weeks |
| Fenofibrate | small molecule | 1-compartment oral | 145 mg daily |

Ten scenarios are provided: natural history, the saturation curve, the
genotype gradient, conventional therapy, each of the four targeted agents,
dietary non-adherence (a single 60 g binge), and pregnancy — the highest-risk
clinical situation in FCS, where oestrogen raises VLDL secretion while lowering
what little lipase there is, and the model correctly back-loads the risk into
the third trimester.

Calibration targets and their sources are listed in a block comment inside the
R file and cross-indexed in `fcs_references.md`.

### Quick start

```r
source("fcs_mrgsolve_model.R")

d <- FCS_scenario_plozasiran(days = 365)   # one arm
FCS_plot_overview(d)

FCS_saturation_curve()      # claim 1 — the cliff
FCS_limb_decomposition()    # claim 2 — the dead limb
FCS_trial_ledger()          # claims 3 and 4 — knockdown and convexity
FCS_jensen_gap()            # what a fasting TG cannot see

results <- FCS_run_all()    # everything, with a printed summary
```

### Dashboard

```r
shiny::runApp("fcs_shiny_app.R")
```

Ten tabs: patient profile, lipid time course, the saturation cliff, limb
decomposition, drug PK and target engagement, clinical endpoints, the Jensen
gap, safety (platelets and ALT), scenario comparison, and biomarkers with the
literature calibration table.

## 5. Honest limitations

- **The hazard function is a fit, not a mechanism.** The Hill exponent of 1.7
  is chosen so that the model brackets the observed placebo pancreatitis rates
  of 0.20–0.30 events per patient-year in the Balance and PALISADE control
  arms. The *shape* (convex, threshold-like) is well supported; the exact
  exponent is not, and every conclusion about the size of the convexity
  dividend inherits that uncertainty. The Shiny app exposes `HILL_AP` as a
  slider for exactly this reason.
- **The model is deterministic and single-subject.** Real FCS trials are
  dominated by between-subject variability, which is why Balance's 6-month
  placebo-adjusted median (−43.5%) is so much smaller than its 12-month
  open-label value (−73.7%). The model is calibrated to the latter and should
  not be read as predicting the former.
- **`IMAX_C3 = 2.2` is the load-bearing parameter.** It encodes how much
  LPL-independent clearance capacity apoC-III is holding back, and it is
  inferred from clinical response rather than measured directly.
- **No cardiovascular endpoint.** The map takes the position that
  chylomicrons are too large to be atherogenic and that FCS therefore differs
  from remnant-rich multifactorial chylomicronemia in this respect. That
  position is debated and deliberately not modelled.
- Adipose and hepatic lipid pools are reporters with lumped constants; they do
  not close a whole-body energy balance and should not be used for one.

---

*Research and education only. Not for clinical decision-making, prescribing,
or regulatory submission.*
