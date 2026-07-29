# Idiopathic Intracranial Hypertension (IIH) — QSP Model

> A Quantitative Systems Pharmacology model of idiopathic intracranial
> hypertension built around a single equation and the bound hiding inside it:
>
> **ICP = I_f · R_out + P_sss**
>
> Every drug ever licensed or trialled for IIH — acetazolamide, topiramate,
> furosemide, GLP-1 receptor agonists, 11β-HSD1 inhibitors — acts on **I_f**,
> the CSF formation rate. **P_sss**, the dural sinus pressure, is *additive*.
> So as any secretion-blocking drug approaches perfection, ICP approaches
> P_sss **and stops**. That floor is a term no drug in the disease touches,
> and in a substantial minority of patients it sits **above** the pressure
> that defines remission.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`iih_qsp_model.dot`](iih_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`iih_qsp_model.svg`](iih_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`iih_qsp_model.png`](iih_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`iih_mrgsolve_model.R`](iih_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`iih_shiny_app.R`](iih_shiny_app.R) |
| 📚 References (72, PubMed-verified) | [`iih_references.md`](iih_references.md) |
| 🔬 Numerical reference implementation | [`iih_reference_check.py`](iih_reference_check.py) |
| 📄 Computed results        | [`iih_model_report.txt`](iih_model_report.txt) |

Every number in this README is printed by `iih_reference_check.py` — a pure
standard-library Python implementation of the same equations the R model
integrates. Nothing here is quoted from memory:

```bash
python3 iih_reference_check.py          # full report (~5 min)
python3 iih_reference_check.py --brief  # headline results only
```

---

## 1. The disease in one paragraph

Idiopathic intracranial hypertension is raised CSF pressure (> 25 cmH₂O)
without a mass lesion, hydrocephalus or venous thrombosis, occurring
overwhelmingly in women of reproductive age with central obesity and recent
weight gain. Its incidence tracks obesity prevalence. Papilloedema threatens
permanent visual field loss — the outcome that matters — while headache
dominates quality of life and often outlives the pressure that caused it.
Three mechanisms are pathologically entangled: **obesity raises
intra-abdominal and thence central venous pressure**, impeding cerebral venous
drainage; **the transverse sinuses are stenosed in almost every patient**,
adding a trans-stenotic pressure gradient upstream of CSF absorption; and
**CSF secretion is modulated** by an IIH-specific androgen signature acting on
choroid-plexus transporters. Treatment reflects the confusion: a carbonic
anhydrase inhibitor titrated to the limit of tolerability, weight loss up to
bariatric surgery, venous sinus stenting, CSF diversion, and optic nerve
sheath fenestration — interventions whose relative merits are argued largely
from single-arm pressure series.

## 2. Why the equation is the model

In steady state the Davson relation holds by conservation of volume:

```
ICP = I_f · R_out + P_sss
```

What makes IIH more than plumbing is that **P_sss is not a constant**. The
transverse sinus is a collapsible tube inside the skull, so CSF pressure
squeezes it:

```
P_sss = P_cv + G(ICP),      γ = dG/dICP > 0
```

Together these are a **positive feedback loop**. Closed-loop sensitivity to
any input is the open-loop sensitivity times 1/(1−γ), and the bound becomes
the self-consistent **venous floor**:

```
ICP_floor = P_cv + G(ICP_floor)
```

This is not a modelling convenience. Lalou 2020 ([PMID 31832847][l]) measured
it: in pseudotumour cerebri, CSF and sagittal sinus pressure track each other
(R = 0.96 at baseline, R = 0.92 during infusion), and during CSF drainage they
keep tracking **until the sinus pressure bottoms out while CSF pressure
continues to fall** — γ > 0 while the sinus is collapsible, γ → 0 once it is
open. The floor equation is the arithmetic of that observation.

[l]: https://pubmed.ncbi.nlm.nih.gov/31832847/

## 3. Model structure — 26 ODE compartments

| Block | Compartments |
|---|---|
| Drug PK/PD (9) | acetazolamide gut/central/CA-effect, topiramate gut/central, exenatide SC/central/effect, furosemide |
| Hydrodynamics (4) | ICP (Marmarou dynamics), P_sss (collapsible sinus), intra-abdominal pressure, central venous pressure |
| Host (4) | body weight, androgen index, R_out, meningeal-lymphatic reserve |
| Devices (2) | venous stent effect with restenosis, shunt patency |
| Optic nerve (2) | papilloedema (axoplasmic stasis), surviving axon fraction |
| Symptoms (3) | central sensitisation, HIT-6, pulsatile tinnitus |
| Exposure (2) | serum bicarbonate, cumulative injurious translaminar exposure |

The ICP relaxation time constant is ~10 minutes, so a lumbar puncture and a
CSF infusion study are *inside* the model rather than bolted on — which is
what makes the identifiability analysis (§6) possible.

The mechanistic map colours every therapy node by **which term of the equation
it acts on**, because that — not potency — decides whether it can work:
blue = secretion (I_f), green = outflow resistance (no licensed drug), orange
= the venous floor (stenting), pink = P_cv via the abdomen (weight loss),
purple = bypasses the equation entirely (shunt, LP, fenestration).

## 4. Calibration — one parameter per published anchor

| # | Parameter | Value | Anchored to |
|---|---|---|---|
| 1 | `ROUT0` | 27.92 mmHg/(mL/min) | IIHTT baseline CSF pressure 357.2 mmH₂O |
| 2 | `EC50_ACZ` | 21.55 mg/L | IIHTT between-arm ICP difference −59.9 mmH₂O ⇒ ε(2.5 g/day) = 0.283 |
| 3 | `KW_ACZ` | 17.29 kg/unit | IIHTT between-arm weight difference −4.05 kg |
| 4 | `KEDON` | 0.914 µm/(mmHg·d) | IIHTT OCT substudy: acetazolamide-arm RNFL −175 µm |
| 5 | `ED50_FRISEN` | 289 µm | IIHTT baseline Frisén grade 2.76 |
| 6 | `PMD_ED` | 2.91 dB | IIHTT acetazolamide-arm PMD gain +1.43 dB |
| 7 | `KRG_ED` | 0.0206 µm/µm | IIHTT acetazolamide-arm RGCL −3.6 µm |
| 8 | `ED_CRIT` | 392 µm | IIHTT mild cohort: no true axon loss at 6 months |
| 9 | `KAX` | 0.00227 /d | fulminant IIH: ~50% axon loss in 60 days |
| 10 | `EC50_GLP` | 0.011 µg/L | exenatide RCT: ICP −5.7 cmH₂O at 2.5 hours |

Physiological constants, not fitted: `IF0` (0.35 mL/min), `ELAST` (0.11 /mL,
Marmarou), `IOP` (15 mmHg), `KTRANS`. Structural choices with no data to fix
them, varied in the sensitivity analysis: the collapse curve (`GMAX`,
`ICPCOL`, `WCOL`), papilloedema resolution time constant (42 d, from the
clinical course), and the entire slow R_out/lymphatic remodelling loop — which
has **no human quantitative anchor in IIH** and is included only because
chronicity and relapse need one.

## 5. The central result: the venous floor is a hard bound

Each phenotype below is **pinned to its published baseline pressure**; what
differs between rows is only the *composition* of that pressure — and the
composition is exactly what no drug trial in IIH has ever measured.

| Phenotype | ICP₀ cmH₂O | R_out | P_sss cmH₂O | γ | floor cmH₂O | max drug ΔICP | remission by drug? |
|---|---|---|---|---|---|---|---|
| normal physiology | 12.0 | 10.6 | 6.9 | 0 | 6.9 | 5.1 | — |
| IIH, purely resistive | 35.7 | 52.2 | 10.9 | 0 | 10.9 | **24.8** | yes |
| IIHTT-like (moderate stenosis) | 35.7 | 27.9 | 22.4 | 0.23 | 18.1 | 17.6 | yes |
| stenting-referred (severe) | 37.0 | 9.9 | 32.3 | 0.37 | 29.2 | **7.8** | **NO** |
| fulminant | 55.0 | 13.4 | 48.6 | 0.28 | 45.5 | 9.5 | **NO** |

Hold one patient at exactly 37 cmH₂O and sweep only the gradient:

| trans-stenotic gradient (mmHg) | 0 | 3.4 | 6.7 | 10.1 | 13.4 | 16.8 |
|---|---|---|---|---|---|---|
| R_out | 52.0 | 42.5 | 32.9 | 23.3 | 13.8 | 4.2 † |
| venous floor (cmH₂O) | 12.2 | 14.7 | 17.8 | 21.6 | 26.7 | 33.5 |
| **max ΔICP by a *perfect* drug** | **24.8** | 22.3 | 19.2 | 15.4 | **10.3** | 3.5 |
| remission on a perfect drug? | yes | yes | yes | yes | **NO** | **NO** |

† R_out below the normal 6 mmHg/(mL/min) range — non-physiological, and that
constraint is what bounds how venous a patient can be.

**The therapeutic room available to every drug in this disease varies more
than two-fold across patients who are clinically indistinguishable without
venography.** In the same stenosis-dominated patient:

| intervention | ICP (cmH₂O) | remission |
|---|---|---|
| untreated | 37.0 | NO |
| acetazolamide 2.5 g/day (ε = 0.283) | 34.8 | NO |
| acetazolamide 4 g/day (ε = 0.378) | 34.1 | NO |
| **perfect secretion blocker (ε = 1.0)** | **29.2** | **NO** |
| venous sinus stenting alone | 18.2 | yes |
| stent + acetazolamide 2.5 g/day | 16.8 | yes |

A drug that abolishes CSF secretion *completely* does not reach remission; a
stent does, with no drug at all. They are not interchangeable at any potency,
because they act on different terms: **the drug scales I_f·R_out, the stent
moves the floor the drug is bounded by.** In a virtual population of 300
admissible patients, the ordering reverses across the composition:

| trans-stenotic gradient | n | acetazolamide 4 g/day | *perfect* drug | stenting |
|---|---|---|---|---|
| 0–5 mmHg (resistive) | 97 | 30% | **100%** | 11% |
| 5–12 mmHg (moderate) | 134 | 12% | 96% | 37% |
| > 12 mmHg (severe) | 69 | **0%** | 36% | **67%** |

A real stented cohort reached < 25 cmH₂O in 40/50 (80%, [PMID 29871989][p]) —
those patients were selected for a demonstrated gradient, i.e. the bottom row,
not the population average.

[p]: https://pubmed.ncbi.nlm.nih.gov/29871989/

## 6. Four further results

**The loop gain is measurable from published data, and it is not zero.**
Stenting lowers sagittal sinus pressure by 8.1 mmHg ([PMID 29922401][b]) and
CSF pressure by 16.8 cmH₂O = 12.4 mmHg ([PMID 29871989][p]) — a ratio of
**1.53**. A passive, non-collapsible sinus forces that ratio to be *exactly
1*: lowering P_sss by x lowers ICP by x and nothing else changes. A ratio
above 1 is only possible if the sinus is collapsible, and 1/(1−γ) = 1.53 gives
**γ ≈ 0.35**. Both known biases across that interval push the other way (the
acetazolamide dose fell from 950 to 300 mg/day and weight *rose* 1.1 kg), so
the gain is if anything under-estimated.

[b]: https://pubmed.ncbi.nlm.nih.gov/29922401/

**ICP is a state; vision is an integral.** Forced to the *same* day-180
pressure, a step intervention (drug) delivers 1.26× the visual gain of a ramp
(weight loss: 1.23 vs 0.97 dB), because permanent loss tracks the *integral*
of translaminar gradient above threshold, which the endpoint pressure has
already forgotten. A six-month ICP endpoint therefore systematically
over-values weight loss — and the weight component of GLP-1 agonists —
against drugs, stents and shunts.

**Timing beats potency where the floor is high.** In the fulminant phenotype,
stenting delayed by 120 days still leaves less permanent damage (79% axon
loss, PMD −15.2 dB) than maximal acetazolamide given immediately (96% axon
loss, PMD −20.5 dB), because the drug cannot cross the floor while the stent
moves it. This is the clinical rule — fulminant IIH needs mechanical
decompression, not dose escalation — arriving as a structural consequence
rather than an assumption.

**The ICP endpoint cannot separate drug potency from loop gain.** It
constrains only the product ε·I_f·R_out/(1−γ); every (γ, ε) pair on that curve
fits the trial exactly. Worse, the classical CSF infusion test measures
**R_out/(1−γ), not R_out**, so in IIH it over-states resistance by exactly the
amplification factor, and any drug effect predicted from it is over-stated
with it. What breaks the tie is not more pressure data but *one* simultaneous
sinus pressure trace, which yields dP_sss/dICP directly. That is this model's
practical recommendation.

## 7. Two discrepancies the model reports rather than absorbs

**① The IIHTT placebo arm falls ~5× more than any weight mechanism allows.**
Its 3.45 kg of weight loss buys at most 0.97 cmH₂O at the randomised IIH:WT
slope (0.28 cmH₂O/kg); the observed fall was 5.24 cmH₂O — an implied slope of
1.52 cmH₂O/kg against 0.28 (IIH:WT 12 mo), 0.31 (24 mo) and 0.51 (very-low-
calorie diet cohort). Fitted to each arm's *absolute* fall, the model's
shortfalls are **−49.7 mmH₂O (acetazolamide) and −46.8 mmH₂O (placebo)**.
Those being nearly equal is the finding: an arm-independent offset is what a
non-treatment effect looks like — enrolment required a raised opening
pressure, and lumbar-puncture opening pressure has large within-subject
variance, which is the textbook setup for regression to the mean.

The methodological consequence is sharp. Reading the trial at face value
requires acetazolamide to suppress CSF formation by **0.53**, at or above the
ceiling for carbonic anhydrase inhibition; the randomisation-protected
between-arm reading requires **0.28**. Same trial, two defensible readings, a
1.9-fold difference in inferred potency — and **a QSP model fitted to absolute
arm trajectories, which is the common practice, silently picks the first and
inherits the artefact as drug potency.** This model uses the second.

**② Weight mediates ~20% of the pressure effect but a reported 4% of the
vision effect, and the model cannot close that gap.** IIHTT's own mediation
analysis put the weight-mediated part of the visual benefit at 0.03 dB of
0.71 dB. The model reproduces the pressure share (18.7% against ~20%) and
**gets the vision share wrong (~16% against 4%)**. The step-versus-ramp timing
asymmetry is real but supplies only 1.26× of the ~5× required. Either the 4%
is a decomposition of noise — IIHTT's *total* PMD effect was 0.71 dB with a
95% CI of 0 to 1.43 and p = 0.050, and ~16% sits well inside that interval —
or weight loss lowers pressure without protecting the nerve proportionally,
which would require the ocular chain to be driven by something other than the
pressure fed to it here. The model cannot decide, and the discrepancy is left
standing.

## 8. Reproducing the trials

Calibrated on the IIHTT **between-arm** differences, because that is the
randomisation-protected estimand:

| between-arm difference | observed | model |
|---|---|---|
| CSF pressure (mmH₂O) | −59.9 | −56.9 |
| weight (kg) | −4.05 | −4.05 |
| Frisén grade | −0.70 | −0.78 |
| PMD (dB) | +0.71 | **+1.31** ✗ |

Out-of-sample tests (never fitted): IIH:WT's ~24% weight loss for remission →
model **29.7%** (average slope 0.272 cmH₂O/kg, within 3% of the trial slope
it was built on); IIHTT's RGCL numbers reproduced with essentially **zero**
true axon loss, which is the correct reading of a mild-visual-loss cohort.

## 9. What this model gets wrong

- **It over-predicts the visual value of pressure control** (between-arm PMD
  1.31 vs 0.71 dB) and of weight loss (~16% vs 4% mediation). The excess is
  the same size as the placebo-arm visual improvement it cannot generate
  (0.12 predicted vs 0.71 observed). Until that is resolved the ocular chain
  here is too responsive to pressure.
- **The exenatide anchor is only just reachable.** The observed −5.7 cmH₂O at
  2.5 hours forces the GLP-1R maximal secretory effect to ~0.5, essentially
  equal to acetazolamide's ceiling, and forces 10 µg bid to sit near
  saturation. That is falsifiable: higher doses or weekly analogues should
  **not** deepen the acute ICP fall, only the weight-mediated ramp. The
  trial's own time course (−5.7 at 2.5 h, −6.4 at 24 h, −5.6 at 12 weeks) is
  also non-monotone in a way this model structurally cannot produce while
  weight is falling.
- **Baseline papilloedema is too high** (456 µm RNFL against typical reported
  IIH values of ~250–350 µm). Only RNFL *change* is compared with the trial,
  and PMD/Frisén are normalised to baseline, so the pressure conclusions do
  not depend on it — but the absolute OCT scale should not be quoted.
- **The remission-threshold prediction is 24% too high** (29.7% vs 24%),
  because the model's pressure-per-kilogram slope is not constant: as pressure
  falls the sinus re-opens, the loop gain drops, and each further kilogram
  buys slightly less.
- **The slow outflow-remodelling loop is unanchored.** It produces chronicity
  and relapse, which the disease has, but no human IIH data constrain it.
- **Headache is a hypothesis, not a prediction.** The large central-
  sensitisation term exists because headache is famously poorly coupled to
  pressure; it is structurally motivated and quantitatively unanchored.
- **Per-patient floors must never be quoted without that patient's own
  manometry.** The sensitivity analysis shows the qualitative bound survives
  across the whole plausible parameter range, and its exact numerical value
  does not.

## 10. Scenarios

Fourteen prebuilt scenarios in `IIH_simulate_scenarios()`: natural history ·
IIHTT placebo · IIHTT acetazolamide 2.5 g/day · acetazolamide at the 4 g/day
target · topiramate · very-low-calorie diet · bariatric surgery · exenatide ·
venous sinus stenting · stent + acetazolamide · **maximal triple drug therapy
against the floor** · CSF shunt · fulminant IIH (drug vs urgent stent) ·
acetazolamide withdrawal and relapse. Plus `IIH_infusion_study()` and
`IIH_lumbar_puncture()`, the two diagnostic manoeuvres that would identify the
model, and `IIH_drug_ceiling_table()`, which computes the bound directly.

---

**Disclaimer.** This is an educational and research QSP model, assembled from
public literature and not independently validated. It must not be used for
clinical decisions, prescribing, or regulatory submission. Its parameters are
illustrative approximations, and its central quantity — an individual
patient's venous floor — requires catheter manometry that none of the drug
trials it is calibrated against performed.
