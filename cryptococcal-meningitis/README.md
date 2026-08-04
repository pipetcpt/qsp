# HIV-associated Cryptococcal Meningitis — QSP Model
### 크립토콕쿠스 수막염 (HIV 연관) 정량적 시스템 약리학 모델

<a href="cm_qsp_model.svg"><img src="cm_qsp_model.png" width="760" alt="Cryptococcal meningitis QSP mechanistic map"></a>

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (145 nodes, 18 clusters) | [`cm_qsp_model.dot`](cm_qsp_model.dot) · [SVG](cm_qsp_model.svg) · [PNG](cm_qsp_model.png) |
| ⚙️ mrgsolve ODE model (53 compartments, 16 scenarios) | [`cm_mrgsolve_model.R`](cm_mrgsolve_model.R) |
| 🐍 Independent Python reference implementation | [`cm_reference_model.py`](cm_reference_model.py) |
| 📊 Shiny dashboard (9 tabs) | [`cm_shiny_app.R`](cm_shiny_app.R) |
| 📚 References (44 verified PubMed citations) | [`cm_references.md`](cm_references.md) |

---

## The one-sentence version

**This disease runs on two clocks, and every trial endpoint measures the fast
one.** Viable yeast clear from the CSF in one to two weeks; the capsular
polysaccharide already dissolved there — which sets outflow resistance, hence
intracranial pressure, and which is the antigen stock immune reconstitution
reacts to — clears with a half-life of about thirteen days and is touched by no
antifungal drug. Early fungicidal activity, the pharmacodynamic endpoint of
every randomised trial in cryptococcal meningitis, cannot see the second clock
at all.

---

## 1. Why build this model

Cryptococcal meningitis kills roughly 112,000 people a year, most of them in
sub-Saharan Africa, and 10-week mortality in the best modern trial arm is still
24.8%. The therapeutic era since 2004 has been organised around a single
pharmacodynamic surrogate — **early fungicidal activity (EFA)**, the slope of
log₁₀ CSF colony counts over the first two weeks — because Brouwer *et al.*
showed that CSF clearance is exponential and that the slope is measurable,
reproducible and prognostic.

That surrogate has been extraordinarily productive. It identified flucytosine
as the right partner drug, it justified 1 mg/kg over 0.7 mg/kg of amphotericin,
and it licensed the single 10 mg/kg liposomal dose that WHO now recommends. But
it measures one variable, and this model exists to ask what the other variables
are doing while it does.

The answer the model gives is that the fungal burden and the antigen pool are
governed by rate constants differing by roughly twenty-fold, so that the
patient who is culture-negative at day 13 still has 85% of their presenting
antigen, still has an opening pressure above 250 mmH₂O, and still carries the
antigen load that will drive IRIS if antiretroviral therapy starts now.

---

## 2. The four structural commitments

Everything the model says follows from four decisions about how to write the
equations, each of which was forced by a published observation.

### 2.1 Intracranial pressure is not a state variable

It is the residual of a volume budget (Davson):

```
ICP        = Pss + Pel(Vex + oedema)          Pel = Pel0 · exp(Eel · V)
absorption = Pel / Rout
dVex/dt    = formation − absorption − leak − drainage
```

Nothing sets ICP directly, so pressure can only be changed by altering
formation, resistance or volume. This matters because it makes "treating the
pressure" a structurally different act from "treating the infection", and the
model can then say which interventions do which.

### 2.2 Outflow resistance is driven by antigen, not by yeast or by inflammation

`Rout` relaxes towards a target set by CSF GXM, with smaller contributions from
CSF leucocytes and IRIS. This is not an arbitrary choice. Bicanic *et al.*
(AIDS 2009) found that higher baseline fungal burden gave higher opening
pressure, but that **high burden was necessary and not sufficient**, and that
baseline pressure was *not* associated with CD4 count or with CSF
pro-inflammatory cytokines. A model in which pressure is driven by the yeast
cannot produce "necessary but not sufficient"; a model in which it is driven by
a slowly-clearing product of the yeast can.

### 2.3 The antigen pool clears about twenty times more slowly than the yeast

Sterility and antigen-negativity are different events, weeks apart. This is why
serum cryptococcal antigen stays positive long after cultures convert, and it
is the mechanism behind everything in section 4.

### 2.4 Death is a hazard integral with separable terms

Burden, time-unsterile, pressure, perfusion, anaemia, neutropenia,
hypokalaemia, renal function, IRIS, injury and steroid harm each contribute an
additive term. Because they are separable, the model can be asked *which clock
killed a virtual patient*, and different interventions answer differently.

---

## 3. Two shared nodes carry opposite signs

These are where the model earns its keep, because in each case a single
physiological variable is pushed in opposite directions by two drugs given
together.

### 3.1 Ergosterol

Fluconazole inhibits ERG11 and depletes membrane ergosterol. That depletion is
how the azole works — growth rate scales as `ERG^gERG`. It is also how the
azole removes amphotericin's binding *substrate*, because amphotericin kills by
extracting ergosterol into an extramembranous sterol sponge (Anderson 2014).
So the model writes:

```
growth        ∝ ERG^gERG                     (the azole benefit)
EC50(AmB)     ∝ 1 / ERG^aERG                 (the azole price)
```

One node, two signs, and no separate "drug interaction" term anywhere.

### 3.2 Glomerular filtration rate

Amphotericin lowers it. Flucytosine is cleared by it. The nephrotoxicity of one
drug is therefore the myelotoxicity of the other, and the interaction is a
shared clearance organ rather than additive marrow suppression.

---

## 4. What the model says

Every number below was computed twice — by `cm_mrgsolve_model.R` (mrgsolve,
LSODA) and by `cm_reference_model.py` (pure-Python fixed-step RK4, no shared
code) — and the two agree to three or four significant figures. See §6.

### 4.1 The two clocks, on the AMBITION regimen

| | value |
|---|---|
| EFA, days 0–14 | **−0.402** log₁₀ CFU/mL/day (trial: −0.40) |
| first negative CSF culture | day **13.5** |
| CSF GXM on that day | **78 µg/mL** = **85%** of the presenting 92 µg/mL |
| day ICP falls below 250 mmH₂O | day **17.9** |
| CSF GXM at week 10 | **3.5 µg/mL** (3.8% of presenting) |
| GXM half-life | **12.6 days** |

The culture is the first thing to normalise, by a margin of four days over
pressure and of weeks over antigen. Every trial stopped measuring at the first
of those three.

### 4.2 The model refutes the hypothesis it was built to test

The model was constructed to test the idea that **killing the fungus releases
capsule in a bolus and therefore drives the post-treatment rise in
intracranial pressure**. It does not, and the arithmetic is simple enough to
state:

| | value |
|---|---|
| CSF GXM pool at presentation | 91.8 µg/mL |
| total capsule carried by every living yeast in the CSF | 45.0 µg/mL |
| → maximum possible rise from complete sterilisation | **+49%** |

Observed in the simulation: GXM peaks at 132 µg/mL on day 5 (a 44% rise) and
ICP peaks at 266 mmH₂O — a rise of 16 mmH₂O over the presenting 250. Real and
measurable, but nowhere near enough to explain the 400–600 mmH₂O pressures that
require repeated drainage.

**The correct explanation is the opposite one.** Pressure does not spike on
treatment; it *persists*, because it is set by a large standing pool decaying on
its own 13-day clock while the yeast that made it are gone in five days. The
antigen pool has been accumulating for the weeks of subacute illness before
presentation; the yeast standing at any moment are a small fraction of what
has been shed. This is why guidelines call for *repeated* lumbar punctures over
weeks rather than one decompression at diagnosis, and why the most fungicidal
regimen is not the one with the least pressure trouble.

### 4.3 The partner-drug question is an ergosterol question

Day *et al.* (NEJM 2013) gave every patient amphotericin 1 mg/kg and found:

| partner | published EFA | model EFA |
|---|---|---|
| none | −0.31 | **−0.313** |
| fluconazole 800 mg | −0.32 | **−0.318** |
| flucytosine 100 mg/kg | −0.42 (Bicanic: −0.56) | **−0.448** |

Fluconazole adds essentially nothing to amphotericin, yet fluconazole alone is
clearly active (−0.11 published, −0.125 modelled). **No model with independent
additive drug effects can produce that pattern.** The ergosterol node produces
it directly: at day 7 the effective amphotericin EC50 is 3.16 µg/mL without
fluconazole and 5.05 µg/mL with it — a 60% loss of potency that almost exactly
cancels the azole's own contribution.

Carried through to mortality, the model gives AmB+5FC 31.8% versus AmB+FLU
36.3% at 10 weeks, a hazard ratio of **0.85** against ACTA's **0.62 (0.45–0.84)**.
So ergosterol antagonism reproduces the *direction* of ACTA's partner-drug
result and about **half its magnitude**. The other half is not explained by this
model, and that shortfall is reported rather than tuned away.

### 4.4 The single liposomal dose is an exposure-matching trick, not a potency one

| | 7 × AmB-d 1 mg/kg | single L-AmB 10 mg/kg |
|---|---|---|
| CNS amphotericin AUC₀₋₁₄ | 3.51 µg·d/mL | **3.39 µg·d/mL** |
| CNS Cmax (day) | 0.37 µg/mL (d7.0) | 0.28 µg/mL (d4.5) |
| days above the effect threshold | 17.2 | **20.2** |
| **renal cortical AUC₀₋₁₄** | **189 µg·d/mL** | **95 µg·d/mL** |
| GFR nadir | 45 mL/min | 51 mL/min |
| potassium nadir | 2.92 mmol/L | 3.04 mmol/L |
| haemoglobin nadir | 8.80 g/dL | 9.18 g/dL |
| EFA | −0.356 | −0.402 |
| 10-week mortality | 30.7% | 26.1% |

Equal CNS exposure, half the renal exposure, and a *longer* time above the
effect threshold from one infusion — which is AMBITION's result (24.8% vs
28.7%, grade 3/4 adverse events 50.0% vs 62.3%) with a mechanism attached.
The one thing that gets *worse* is neutrophils (nadir 1.01 vs 1.83), because
the liposomal arm gives flucytosine for fourteen days rather than seven.

### 4.5 Amphotericin nephrotoxicity is flucytosine myelotoxicity

The same 100 mg/kg/day of flucytosine, beside three different companions:

| regimen | GFR nadir | peak plasma 5FC | CSF 5FC AUC₀₋₁₄ | neutrophil nadir |
|---|---|---|---|---|
| fluconazole + 5FC (oral, 14 d) | 92 | **33.5 µg/mL** | 337 µg·d/mL | 1.68 |
| AmB-d 7 d + 5FC 7 d | 45 | **42.7 µg/mL** | 199 µg·d/mL | 1.83 |
| AmB-d 14 d + 5FC 14 d | 41 | **53.3 µg/mL** | 432 µg·d/mL | **0.90** |

A 59% higher flucytosine peak from an identical dose, because amphotericin
removed the organ that clears it. This is why the all-oral ACTA arm was better
tolerated than the two-week amphotericin arms despite giving *more* total
flucytosine, and it predicts that flucytosine dose reduction should be
indexed to measured GFR during amphotericin rather than to body weight alone.

### 4.6 ART timing is an antigen question, not a CD4 question

IRIS activity is written as the **product** of the rate of immune recovery and
the antigen stock present when recovery begins. Nothing in the model knows
about "ART timing":

| ART start | 26-week mortality | peak IRIS activity | CSF GXM at ART start |
|---|---|---|---|
| day 7 | **61.7%** | 0.322 | 129 µg/mL |
| day 10 | 59.4% | 0.282 | 112 µg/mL |
| day 14 | 56.6% | 0.236 | 90 µg/mL |
| day 21 | 52.2% | 0.169 | 62 µg/mL |
| day 28 | 48.5% | 0.120 | 42 µg/mL |
| day 35 | **45.6%** | 0.084 | 28 µg/mL |
| day 56 | 40.6% | 0.028 | 9 µg/mL |

COAT found 45% with ART at 1–2 weeks against 30% at 5 weeks — a 15-point gap.
The model gives 16 points across the same interval. Deferring ART does not make
the immune response smaller; it makes the antigen that response meets smaller.

The same structure predicts COAT's subgroup finding. In the paucicellular
phenotype (CD4 12, CSF leucocytes 8.7/µL) the peak IRIS drive on early ART is
0.57 against 0.28 in the median patient — because absent CSF inflammation means
a failure of the clearance machinery, so antigen persists, so the stock is
larger when ART arrives. The model's all-cause hazard ratio for early versus
deferred ART is 1.53 in that subgroup against 1.48 in the median patient, which
is a much weaker amplification than COAT's 3.87 versus 1.73; the model's
baseline hazard in that phenotype is already near-saturating, which compresses
the ratio.

### 4.7 Sertraline had to fail

| | value |
|---|---|
| plasma sertraline at 400 mg/day | 0.33 µg/mL |
| free brain concentration, peak | **33.9 ng/mL** |
| in-vitro EC50 for *C. neoformans* | 6 µg/mL |
| free brain concentration as % of EC50 | **0.56%** |
| model change in EFA | **−0.4015 → −0.4021** (0.0006 log₁₀/day) |
| model change in 10-week mortality | 26.15% → 26.14% |

ASTRO-CM enrolled 460 patients to find no benefit. The calculation that
predicts that result is one line of pharmacology — a free concentration against
an MIC — and it did not require the trial. This is the clearest example in the
model of a question that a QSP analysis answers before a phase 3 does.

### 4.8 Dexamethasone: the model reproduces the harm and locates part of it

| | AMBITION | + dexamethasone |
|---|---|---|
| EFA | −0.402 | **−0.308** |
| CSF leucocytes, day 14 | 8.9/µL | 1.0/µL |
| day of first negative culture | 13.5 | **23.3** |
| 10-week mortality | 26.1% | **33.8%** |
| permanent disability index | 0.059 | **0.089** (×1.51) |

CryptoDex was stopped for harm: mortality 47% vs 41%, good outcome 13% vs 25%,
and *slower CSF fungal clearance*. The slower clearance is the part the model
produces mechanistically, from glucocorticoid suppression of the macrophage
killing term. The mortality and disability signals required explicit
steroid-attributable harm terms, which are stated as such rather than hidden —
the model does not claim to derive them.

### 4.9 The one therapy that touches the slow clock is a needle

| | AMBITION | + 7 therapeutic LPs |
|---|---|---|
| peak ICP | 266 mmH₂O | **252 mmH₂O** |
| pressure-time integral above 250 mmH₂O | 174 mmH₂O·day | **0.8 mmH₂O·day** |
| CSF volume drained | 0 | 138 mL |
| GXM removed | 0 | **11.1 mg** |
| peak CSF GXM | 132 µg/mL | 113 µg/mL |
| EFA | −0.402 | −0.433 |
| 10-week mortality | 26.1% | **19.4%** |

Lumbar puncture is the only intervention in the model that removes antigen
rather than waiting for it to decay, which is why it is the only one whose
benefit does not depend on antifungal potency — and why Rolfes *et al.* found
the survival benefit of therapeutic LP to be present **regardless of opening
pressure**. In the high-burden phenotype the effect is larger still: 10-week
mortality 70.0% → 49.4%.

The model's 11-day relative risk for at least one therapeutic LP is 0.75,
against Rolfes' adjusted RR of 0.31. The observational estimate is confounded
in the direction of understating benefit (patients selected for LP had higher
burden and worse mental status), so the model is conservative relative to it,
but the gap is large and is the model's least-constrained major claim.

### 4.10 Four virtual patients, each its own equilibrium

Each phenotype is produced by relaxing the slow variables at a different
clamped presenting burden and CD4 count — none is a hand-edited state vector.

| phenotype | log₁₀ CFU/mL | CSF GXM | opening pressure | CSF WBC | CD4 | EFA | 10-wk death | with LPs | ICP integral |
|---|---|---|---|---|---|---|---|---|---|
| median trial participant | 5.00 | 92 | 250 | 7.7 | 25 | −0.402 | 26.1% | 19.4% | 173 |
| upper-quartile burden | 5.90 | 958 | 327 | 68 | 25 | −0.451 | 70.0% | **49.4%** | 3138 |
| paucicellular (COAT high-risk) | 5.30 | 186 | 284 | 8.7 | 12 | −0.400 | 43.1% | 26.1% | 919 |
| partially preserved immunity | 4.30 | 17 | 162 | 3.0 | 150 | −0.399 | 13.3% | 12.4% | 0 |

Note that EFA is nearly identical across all four (−0.40 to −0.45) while
10-week mortality ranges from 13% to 70%. **The surrogate that the trials
optimise is the variable that discriminates least between these patients.**
That is the single most compact statement of what this model is for.

---

## 5. Calibration — what was fitted to what

### 5.1 Pharmacodynamics: three parameters, three targets, six predictions

`CLbrD`, `aERG` and `CLbrL` were solved by sequential one-dimensional
bisection against three EFA values that each isolate one of them. Everything
else is a prediction.

| regimen | model EFA | published | source | status |
|---|---|---|---|---|
| AmB-d 1 mg/kg alone | **−0.313** | −0.31 | Day 2013 group 1 | **solved** → `CLbrD` |
| AmB-d 1 mg/kg + FLU 800 | **−0.318** | −0.32 | Day 2013 group 3 | **solved** → `aERG` |
| L-AmB 10 mg/kg single | **−0.402** | −0.40 | AMBITION | **solved** → `CLbrL` |
| fluconazole 800 | −0.101 | −0.07 | Longley 2008 | predicted |
| fluconazole 1200 | −0.125 | −0.11 | Nussbaum 2010 | predicted |
| fluconazole 1200 + 5FC | −0.273 | −0.28 | Nussbaum 2010 | predicted |
| AmB-d 0.7 mg/kg + 5FC | −0.383 | −0.40 | Bicanic 2008 (ratio) | predicted |
| AmB-d 1 mg/kg + 5FC | −0.448 | −0.42 … −0.56 | Day 2013 / Bicanic 2008 | predicted |
| AmB-d 1 mg/kg + 5FC ×7 d | **−0.356** | **−0.42** | AMBITION control | **worst residual** |

The last row is the model's largest EFA error and is not fixed. The 1-week
amphotericin arm underperforms because CNS amphotericin washes out with a
4.3-day half-life while fluconazole 1200 mg starts on day 8 and antagonises
what remains.

The amphotericin kill term is **nearly first-order in concentration**
(EC50 ≫ achievable C). This is not an aesthetic choice: Bicanic's within-trial
comparison of 0.7 against 1.0 mg/kg (EFA −0.45 vs −0.56) is a 1.43-fold dose
step producing a 1.24-fold rate step, and no saturating model can reproduce a
ratio that close to linear. It is the most exposed pharmacodynamic assumption
in the model.

### 5.2 Mortality: two-stage, and deliberately under-fitted

Seven hazard coefficients that no induction-therapy trial can identify — none
of them randomised pressure management, potassium replacement or transfusion —
were **fixed** at mechanistically anchored values. Only `hB` (burden), `hSTER`
(days culture-positive) and `hAMS` (neurological injury) were fitted, by
non-negative least squares, to 17 mortality endpoints from five randomised
trials.

| regimen | 2-wk model | 2-wk pub | 10-wk model | 10-wk pub | source |
|---|---|---|---|---|---|
| untreated | 44.8% | — | 96.9% | ~100% | historical |
| fluconazole 800 | 34.7% | — | **57.7%** | 60% | Longley 2008 |
| fluconazole 1200 | 32.7% | 30% | **53.7%** | 56% | Nussbaum 2010 |
| ACTA oral arm | 21.2% | 18.2% | **34.6%** | 35.1% | ACTA 2018 |
| AmB-d alone | 25.0% | 24.8% | **40.7%** | 43.6% | Day 2013 g1 |
| AmB-d + FLU | 23.1% | 19.7% | **36.3%** | 32.6% | Day 2013 g3 |
| AmB-d + 5FC (2 wk) | 19.2% | 15.0% | **31.8%** | 30.0% | Day 2013 g2 |
| AmB-d + 5FC (1 wk) | 19.0% | 15.5% | **30.6%** | 24.2 / 28.7% | ACTA / AMBITION |
| AMBITION L-AmB | 16.3% | ~13% | **26.1%** | 24.8% | AMBITION 2022 |

**Mean absolute error at the primary 10-week endpoint: 2.1 percentage
points**, from two free parameters across nine regimens and five trials. For
scale, ACTA and AMBITION disagree by 4.5 points about the *same* regimen, and
Day and Bicanic disagree by 0.14 log₁₀/day about the *same* EFA.

**A stated defect:** the model over-predicts 2-week mortality by 3–4 points
across the board. The cause is identified — the burden hazard is instantaneous
rather than integrated over the preceding days, so it front-loads deaths — and
it is not fixed, because the fix would require a lag structure that no
published data constrain.

---

## 6. Verification — two implementations, and the three bugs the comparison found

The model is implemented twice with no shared code: `cm_mrgsolve_model.R`
(mrgsolve, C++/LSODA) and `cm_reference_model.py` (pure-Python fixed-step RK4,
standard library only). Both were executed and compared.

| quantity | mrgsolve | Python |
|---|---|---|
| presenting CSF GXM | 91.5 | 91.78 µg/mL |
| presenting opening pressure | 250 | 250.1 mmH₂O |
| EFA, AmB-d 1 mg/kg alone | −0.3133 | −0.313 |
| EFA, AMBITION | −0.4016 | −0.402 |
| EFA, AmB-d + 5FC | −0.4482 | −0.448 |
| 10-week mortality, AMBITION | 26.09% | 26.1% |
| day of first negative culture | 13.50 | 13.5 |
| peak plasma 5FC, 2-wk AmB+5FC | 53.28 | 53.3 µg/mL |
| free brain sertraline | 33.85 | 33.8 ng/mL |
| 26-week mortality, ART day 7 | 61.7% | 61.6% |
| high-burden phenotype GXM | 955.9 | 958 µg/mL |

The comparison was not decorative. It exposed three real defects:

1. **Two `for` loops in the ODE block declaring the same index name.**
   mrgsolve hoists loop declarations into one scope and the compiler rejected
   it. Python never saw this.
2. **LSODA taking flucytosine very slightly negative after the last dose**,
   so that `pow(Cfc, 1.8)` in the marrow-suppression term returned `NaN` and
   *silently destroyed the mortality output of every arm that sterilised*. All
   fractional powers are now guarded. The fixed-step Python integrator never
   triggered this, so the bug was only visible from the other side.
3. **The Python burn-in was 70 days and had not converged.** The R
   implementation, started from the Python state and run again, drifted to
   exactly Python's 140-day values. Rather than paper over it, the 70-day
   clamp is documented as the modelling choice it is — and the 6 mmH₂O of extra
   pressure that full convergence adds (opening pressure 250 → 256 mmH₂O,
   GXM 91.8 → 101.6 µg/mL) is reported as a result: it is the model's estimate
   of how much presenting pressure depends on how long the patient has been ill.

One intentional difference remains: the therapeutic lumbar puncture is a
proportional drain over a 30-minute window in R and an instantaneous
pressure-target reset in Python. This is the only place the two disagree
materially — 138.1 mL drained versus 140.0, and 19.7% versus 19.4% 10-week
mortality.

---

## 7. Model contents

**53 ODE compartments** · **173 parameters** · **16 scenarios**

| group | compartments |
|---|---|
| Amphotericin B PK | `Ad` `Ad2` `Al` `Al2` `Abr` `Akid` — deoxycholate and liposomal disposition, separate CNS effect site and renal cortex |
| Flucytosine PK | `FCg` `FCc` `FCcsf` — GFR-dependent clearance |
| Fluconazole PK | `FLg` `FLc` `FLcsf` |
| Adjuncts | `DXg` `DXc` `SRg` `SRc` `SRbr` `IFNsc` |
| Fungus (clock 1) | `Fe` `Fres` `Ft` `Fi` `Fp` — extracellular, 5FC-resistant, persister, intracellular, parenchymal |
| Antigen (clock 2) | `GXM` |
| Drug target | `ERG` — the oppositely-signed node |
| Host immunity | `CD4` `VL` `MAC` `TH1` `IFNG` `PROIN` `IL10` `WBC` `IMML` `IRISa` |
| CNS hydrodynamics | `Vex` `Rout` `EDEMA` `LEAK` |
| Safety | `GFR` `Kser` `Hb` `ANC` `ALT` |
| Outcome | `NEUR` `DIS` `HAZ` + six bookkeeping integrals |

**Scenarios:** untreated · fluconazole 800 · fluconazole 1200 · ACTA oral arm ·
AmB-d alone (Day g1) · AmB-d + fluconazole (Day g3) · AmB-d 0.7 + 5FC (Bicanic
g1) · AmB-d 1.0 + 5FC (Day g2) · AmB-d + 5FC one week (ACTA/AMBITION control) ·
AMBITION single L-AmB · AMBITION + therapeutic LPs · AMBITION + dexamethasone ·
AMBITION + sertraline · AMBITION + IFN-γ · COAT early ART · COAT deferred ART.

---

## 8. Running it

```bash
# Mechanistic map
dot -Tsvg cm_qsp_model.dot -o cm_qsp_model.svg
dot -Tpng -Gdpi=150 cm_qsp_model.dot -o cm_qsp_model.png

# Python reference: full calibration report (~2 min, no dependencies)
python3 cm_reference_model.py
python3 cm_reference_model.py --brief

# mrgsolve model: calibration table
Rscript cm_mrgsolve_model.R --calibrate

# Shiny dashboard
Rscript -e 'shiny::runApp("cm_shiny_app.R")'
```

---

## 9. What would falsify this model most cheaply

The model's weight rests on one unmeasured quantity and one inferred one, and
both are testable:

1. **Serial CSF GXM concentration during induction therapy, alongside
   quantitative culture.** The 13-day antigen half-life that drives the entire
   slow clock is inferred from the well-documented lag of CrAg titre behind
   culture conversion — it has never been measured directly against paired
   colony counts. If CSF GXM in fact falls with the yeast, the model's central
   claim collapses. This is a straightforward assay on samples that are already
   being collected.

2. **Amphotericin concentration in human CSF or brain.** The CNS delivery
   parameters were solved *backwards* from early fungicidal activity, so the
   model's brain concentrations are an inference from effect, not an
   observation. The implied 36-fold difference in per-unit-plasma CNS delivery
   between deoxycholate and liposomal amphotericin is a prediction awaiting a
   measurement.

3. **A randomised trial of scheduled versus symptom-driven therapeutic lumbar
   puncture.** The pressure hazard coefficients are fixed rather than fitted
   because no trial has randomised pressure management. Every quantitative
   claim in §4.9 is an extrapolation constrained by a single confounded
   observational estimate, and it is the model's largest single intervention
   effect.

---

## 10. Limitations

- The ART-timing curve has **no optimum** because the model contains no
  competing risk from leaving HIV untreated. It must not be read beyond about
  eight weeks.
- The partner-drug hazard ratio (0.85) recovers only about half of ACTA's
  (0.62). Ergosterol antagonism is a partial explanation.
- The paucicellular subgroup's amplification of ART-timing harm (HR 1.53 vs
  1.48) is far weaker than COAT's (3.87 vs 1.73), because the model's baseline
  hazard in that phenotype is already near-saturating.
- Dexamethasone's mortality and disability harms required explicit added terms;
  only the slower fungal clearance is derived.
- The model is single-organ for the CNS: pulmonary disease, cryptococcomas
  outside the brain, and the eye are represented only as a lumped parenchymal
  reservoir.
- The perfusion-pressure hazard term is present but **unidentified** — no arm in
  the calibration set ever reaches a cerebral perfusion pressure below
  60 mmHg, so its coefficient is a prior, not a fit.

---

## ⚠️ Disclaimer

This is an educational and research model built from public literature. It has
not been independently validated or certified and **must not be used for
clinical decisions, prescribing, or regulatory submission.** Parameters are
illustrative approximations; fitting and validation against real patient data
would be required for any applied use.
