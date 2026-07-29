# 항암화학요법 유발 호중구감소증 · 발열성 호중구감소증 (CIN / FN)
### Chemotherapy-Induced Neutropenia and Febrile Neutropenia — QSP model

<a href="cin_qsp_model.svg"><img src="cin_qsp_model.png" width="640" alt="CIN/FN mechanistic map"></a>

---

## The one idea this model is built around

**The absolute neutrophil count is not a state you can dose against, because by
the time it moves the event that set it has already happened.** The ANC measured
on day 10 was decided by the kill in the proliferating compartment on days 0–2,
and everything in between is a delay line:

```
PROL  ->  TR1  ->  TR2  ->  TR3 (marrow storage)  ->  CIRC        (Friberg 2002)

   the DELAY (mean transit time, 89.3 h)   sets WHEN the nadir arrives
   the EXPOSURE (slope x Cp)               empties the pool
   the FEEDBACK exponent gamma (0.170)     sets HOW FAST it comes back
```

Those three are separate parameters, and nearly every clinical rule about
myelosuppression turns out to be a consequence of that separation rather than of
potency. The model exists so that the rules can be **derived** rather than
asserted — and so that the places where the derivation fails can be pointed at.

Every number below was produced by `cin_reference_impl.py`, a dependency-free
pure-stdlib re-implementation of `cin_mrgsolve_model.R`. Its full 742-line
output is committed as [`cin_numerical_report.txt`](cin_numerical_report.txt).
Reproduce with `python3 cin_reference_impl.py`.

---

## Deliverables

| File | What it is |
|------|-----------|
| [`cin_qsp_model.dot`](cin_qsp_model.dot) · [`.svg`](cin_qsp_model.svg) · [`.png`](cin_qsp_model.png) | Mechanistic map — **227 nodes, 323 edges, 19 clusters** |
| [`cin_mrgsolve_model.R`](cin_mrgsolve_model.R) | **53-compartment** mrgsolve model, 8 agents, 13 regimens, `mread()`-ready |
| [`cin_shiny_app.R`](cin_shiny_app.R) | **9-tab** dashboard, arranged so the ANC is never shown alone |
| [`cin_references.md`](cin_references.md) | **114 references**, every PMID resolved through the PubMed API |
| [`cin_reference_impl.py`](cin_reference_impl.py) | Pure-stdlib reference implementation + calibration + 15 analyses |
| [`cin_numerical_report.txt`](cin_numerical_report.txt) | Committed output — the source of every number here |

---

## Read this before any of the results

**Two parameters are fitted.** Both to Vogel 2005 (PMID 15718314): febrile
neutropenia 17% on placebo and 1% on pegfilgrastim after docetaxel 100 mg/m².
Both are hit exactly, which is what one parameter per number always does and is
evidence of nothing.

**There are only two because of a trade-off the model cannot escape.** The
overall myelotoxic scale can be set *either* so that the absolute durations
reported in the G-CSF trials come out right, *or* so that the model tells
regimens apart. Not both. The duration below the 0.5 ×10⁹/L threshold saturates
in this topology: the recovery limb starts from a floor set by φ that does not
depend on exposure, and climbs at a rate set by the feedback exponent, so once
the proliferative pool is empty at all, extra exposure buys almost no extra
duration. Raising the scale until CAE reaches Crawford 1991's 6-day median puts
**every** combination regimen at a nadir of 0.065–0.11 and a DSN of 5.4–6.8 d —
it stops telling AC apart from TAC. The scale is set for **discrimination**, and
the shortfall is reported.

### The six misses

| Reported | Source | Model | Error |
|---|---|---|---|
| CAE, no G-CSF, median DSN 6.0 d | Crawford 1991 | 4.64 d | −23% |
| CAE + filgrastim d4–17, median DSN 1.0 d | Crawford 1991 | 1.90 d | +90% |
| AT + pegfilgrastim, cycle-1 DSN 1.8 d | Green 2003 | 0.00 d | −100% |
| AT + daily filgrastim, cycle-1 DSN 1.6 d | Green 2003 | 1.46 d | −9% |
| Intact marrow + filgrastim ×5 d, ANC 4–6× | Lord 1989 / Dale 2018 | 1.9× | −62% |
| AC q3w P(FN), guideline band 10–20% | NCCN | 0.198 | top of band |

Read together they say two specific things:

1. **The model exaggerates how much the timing of G-CSF matters.** Given on day
   2, pegfilgrastim removes severe neutropenia completely (0.00 d against 1.8
   observed). Given on day 4, filgrastim leaves 1.90 d where 1.0 was observed.
   One is −100% and the other +90%, in opposite directions, and the axis they
   differ on is the start day. So the timing result in A4 is right in
   **direction** and too steep in **magnitude**.
2. **It ranks the two products in the wrong order.** Green found pegfilgrastim
   marginally *worse* than daily filgrastim (1.8 vs 1.6 d); the model has it
   better. The mechanism responsible is the one the model is proudest of —
   target-mediated clearance holding the long-acting product up for exactly as
   long as the ANC stays low — so the same term that produces the self-titration
   result also produces this error.

The intact-marrow miss is **forced**: the steady-state ANC multiple in a Friberg
loop is `((1−φ)(1+Emax_amp))^(1/γ)` with 1/γ = 5.9, so a 5-fold rise needs
`Emax_amp = 0.320` — and at that value pegfilgrastim abolishes severe
neutropenia altogether, which it demonstrably does not. `Emax_amp = 0.10` sits
between two observations the model cannot hold at once.

**What survives all of this is everything comparative.** A4–A14 change one thing
at a time within a fixed regimen and patient, and the quantity being compared —
timing, product, dose, cell-cycle dependence, count versus function — does not
depend on the absolute scale being right. Absolute P(FN) values do, and they are
anchored to a randomised-trial population, which Truong 2016 shows understates
real-world rates.

---

## Results

### 1. Three quantities, three owners — and the depth belongs to none of them

A 3-fold dose range on CAE:

| dose × | nadir ANC | nadir day | DSN (d) | P(FN) |
|---|---|---|---|---|
| 0.50 | 0.077 | 9.71 | 5.97 | 0.191 |
| 1.00 | 0.065 | 9.83 | 6.13 | 0.208 |
| 1.50 | 0.059 | 9.92 | 6.23 | 0.218 |

I expected the dose to move the depth. It moves it 1.31-fold across a 3-fold
range, moves the duration 1.04-fold, and moves the nadir **day** by 0.21 d.
Transit time, swept over the published 62–125 h range at fixed dose, moves the
nadir day 7.88 → 12.63 and the depth barely at all. So the transit time owns the
**day**, and the depth is owned by neither dose nor transit — it is floored by φ
and the feedback. That is why a growth factor given on day 8 cannot help much:
by day 8 the deficit has been manufactured and is in transit.

### 2. The endpoint is a duration, and the nadir ranks patients badly

Varying **only** the feedback exponent — changing how fast the patient recovers
and nothing about how hard the drug hit:

| γ | nadir ANC | DSN (d) | P(FN) |
|---|---|---|---|
| 0.119 | 0.048 | 7.97 | 0.272 |
| 0.170 | 0.065 | 6.13 | 0.208 |
| 0.221 | 0.086 | 4.99 | 0.163 |

A 1.8-fold spread in nadir maps to a 1.7-fold spread in risk that is driven by
the **integral**, exactly as Bodey's 1966 observation implies. Reporting a nadir
without a duration discards the variable that matters.

### 3. Pegfilgrastim is cleared by the cells it creates

Same fixed 6 mg dose in every row, ordered by how bad the regimen is without
support:

| regimen | DSN alone | DSN + peg | t½ day 5–8 | t½ day 15–18 | Cmax |
|---|---|---|---|---|---|
| docetaxel 75 | 3.94 | 0.00 | **27 h** | 21 h | 174 |
| AT | 4.41 | 0.00 | 42 h | 21 h | 174 |
| EP | 4.47 | 0.54 | 48 h | 21 h | 167 |
| CAE | 4.70 | 1.07 | 56 h | 21 h | 174 |
| topotecan | 5.36 | 2.53 | **92 h** | 21 h | 174 |

The half-life measured **during the nadir** rises 3.4-fold from the mildest to
the most severe regimen while the half-life measured after recovery does not
move. Reported pegfilgrastim half-lives span roughly 15–80 h in exactly this
way. This is not between-patient variability and not a covariate: it is the same
feedback loop read from the drug's side. A patient whose marrow is in trouble
automatically receives a longer exposure, because the mechanism that would have
cleared the drug is the mechanism that is missing — which is why pegfilgrastim
is a flat 6 mg regardless of weight while filgrastim is dosed per kilogram.
Filgrastim's own exposure is identical (Cmax 7.93, AUC 80) in a mild and a
severe regimen, because its renal route does not care about the ANC.

### 4. The ANC on a growth factor over-reads marrow output

Five G-CSF actions switched on one at a time, in an **intact** marrow:

| actions active | ANC peak | fold | storage pool at peak |
|---|---|---|---|
| none | 5.00 | 1.00 | 9.79 |
| storage release only | 7.76 | 1.55 | 2.92 |
| survival only | 5.92 | 1.18 | 9.62 |
| transit only | 7.81 | 1.56 | 9.15 |
| amplification only | 6.70 | 1.34 | 13.14 |
| all five | 9.53 | 1.91 | 2.28 |

Storage release alone accounts for more than half the rise and adds no cells: it
empties the marrow reserve into the blood, which is visible in the last column.
Repeat the decomposition in a marrow CAE has just emptied and the ordering
inverts — **transit** does almost all the therapeutic work (DSN 4.60 → 2.52 on
its own), because release and demargination are bookkeeping operations on a pool
that is empty. The pipeline, not the factory.

### 5. Timing beats dose by a wide margin

Pegfilgrastim 6 mg on CAE, one dose, start time swept:

| start (d) | nadir | DSN (d) | P(FN) | G-CSF AUC (×10³) |
|---|---|---|---|---|
| 0.0 (same day) | 0.367 | 1.58 | 0.035 | 9.3 |
| 1.0 | 0.433 | 1.07 | 0.016 | 20.8 |
| 2.0 | 0.483 | 0.53 | 0.010 | 35.3 |
| 4.0 | 0.560 | 0.00 | 0.010 | 49.4 |
| 7.0 | 0.327 | 1.91 | 0.042 | 37.7 |
| 10.0 | 0.169 | 3.79 | 0.245 | 20.1 |

Timing moves DSN across 0.53–3.79 d. A **16-fold** dose range at the best timing
moves it 1.74 → 1.03 d. Note also the AUC column: the same dose delivers 4-fold
more exposure when the marrow needs it, which is item 3 again.

**Same-day dosing** is worse in the model without any parameter saying so:
G-CSF raises the proliferation rate, the cytotoxic kills in proportion to that
rate, and if the two overlap the drug simply gets more cells. Switching the
amplification action off shrinks the penalty but does not abolish it —
amplification is the larger share, not the whole mechanism.

### 6. The actual decision — and it is not about neutrophils

Six cycles of CAE, expected number of FN episodes rather than P(at least one),
because the latter saturates above 0.95 for every strategy:

| strategy | RDI | DSN/cycle | expected FN | tumour @18 wk | Hb end |
|---|---|---|---|---|---|
| no support, full dose | 100% | 4.65 | 1.67 | 2.68 | 9.6 |
| dose 80% | 80% | 4.58 | 1.60 | **7.31** | 11.7 |
| dose 60% | 60% | 4.50 | 1.50 | **20.17** | 13.0 |
| pegfilgrastim d2, full dose | 100% | 1.10 | **0.10** | 2.68 | 9.6 |
| levofloxacin only | 100% | 4.42 | 1.08 | 2.68 | 9.6 |

Dose reduction buys very little FN risk and costs a great deal of tumour
control. The growth factor buys almost all of it and costs none. The tumour
layer is a deliberately thin log-kill and should be read as a **price list**,
not a survival prediction — but the direction of the trade is the whole argument
for prophylaxis, and it is a health-economic argument rather than a
haematological one.

**A13 makes that explicit.** Sweeping the baseline risk and adding
pegfilgrastim at each level: the *risk ratio* is roughly flat (9–23×) while the
*absolute* reduction, and therefore the NNT, scales with baseline risk. **There
is no kink at 20%, or anywhere else.** The threshold is a statement about the
price of a hospital admission divided by the price of a syringe, which is
exactly why it moved from 40% to 20% without any change in the biology
(Calhoun 2005; Cosler 2005).

### 7. Trilaciclib — and a prediction that corrected my expectation

| regimen | cycle-dep | DSN alone | DSN + trilaciclib | DSN saved |
|---|---|---|---|---|
| GCb | 0.75 | 6.53 | 3.90 | **2.64** |
| EP | 0.75 | 4.48 | 3.97 | 0.50 |
| topotecan | **0.95** | 5.35 | 5.31 | **0.05** |
| AC | 0.53 | 4.29 | 4.29 | 0.00 |

I expected the saving to rise monotonically with cell-cycle dependence. It does
not. Topotecan has the highest cycle dependence of any regimen here and gains
almost nothing. The missing variable is **where the regimen sits relative to the
threshold**: trilaciclib removes a *fraction* of the kill, and on a regimen whose
nadir is far below 0.5 a fractional reduction still leaves it below. The benefit
is largest for regimens sitting *near* the threshold. That is a sharper claim
than the one I started with, and it puts the licensed indication
(etoposide/carboplatin in small-cell lung cancer) inside the favourable window
while relapsed topotecan, where the drug is also used, sits outside it.

Trilaciclib is also the **only** intervention here that improves the cumulative
column: six cycles of EP end with marrow reserve 0.9958 versus 0.9934, because
the arrest protects the cycling stem fraction. A growth factor cannot do that —
it accelerates the pipeline downstream of the damage.

### 8. One insult, three lineages, three clocks

| regimen | ANC nadir day | platelet nadir day | Hb nadir day |
|---|---|---|---|
| CAE | 8.8 | 17.3 | 28.0 |
| EP | 8.8 | 19.9 | 28.0 |
| GCb | 11.1 | 21.1 | 28.0 |

The ordering holds for every regimen and no parameter says so — it follows from
transit times of 89/200/150 h and lifespans of 7.9 h / 10 d / 120 d. Across six
cycles the ANC nadir is essentially **reproducible** (0.169 → 0.151), the
platelet nadir is flat after cycle 1, and the haemoglobin falls every cycle and
never recovers (12.3 → 9.1 g/dL). A cell that lives 120 days cannot replace a
cycle's losses inside 21 days; one that lives 8 hours has no memory of the
previous cycle at all. Anaemia is the cumulative toxicity of a course,
neutropenia the acute toxicity of a cycle, and the growth factor that fixes the
second does nothing for the first.

### 9. The count and the defence move in opposite directions

Dexamethasone premedication, at 24 h:

| | ANC | marginated | function | **effective defence** | P(FN) |
|---|---|---|---|---|---|
| no dexamethasone | 4.87 | 4.90 | 0.872 | **8.52** | 0.1700 |
| dexamethasone | 5.09 | 1.99 | 0.729 | **5.16** | 0.1649 |

The measured ANC goes **up** and the effective defence goes **down** — exactly
the kind of number that gets a patient through a day-1 count gate. And in the
infection module, **levofloxacin changes no haematological number at all** and
still moves P(FN) from 0.265 to 0.215, because it acts on the other side of the
race. Reporting only the ANC would make the antibiotic look inert and the
steroid look protective. Both readings would be wrong. This is why the model
carries neutrophil **function** as a state distinct from neutrophil **count**,
and why the Shiny app never plots one without the other.

### 10. A clinician in the loop

Six cycles of CAE with a controller implementing the real rules (hold for
ANC < 1.5 or platelets < 100, reduce after two holds or after the cycle's FN
risk crosses 20%):

| policy | RDI | expected FN | tumour | reserve | Hb |
|---|---|---|---|---|---|
| no prophylaxis | **61.6%** | 1.47 | **16.96** | 0.992 | 13.6 |
| secondary | 94.7% | 0.65 | 2.90 | 0.985 | 11.6 |
| primary | 94.7% | 0.10 | 2.90 | 0.985 | 11.6 |

Note **which rule fires**. Neither patient trips the day-1 count gate often — by
day 21 the ANC has overshot its own baseline, which is what the feedback loop
does — so almost all of the lost dose intensity comes from a clinician looking
at a nadir and reducing the next dose. The intensity is being spent on a number
the growth factor could have fixed without touching the dose.

And the cost that the growth factor does **not** fix: in the compromised
patient, haemoglobin and marrow reserve are *worse* in the primary-prophylaxis
arm precisely because it delivered more chemotherapy. Preserving dose intensity
is not free; it is a decision to spend one toxicity to buy another.

---

## Numerical verification

Nothing above is believed before A0 passes.

* **Baseline steady state** — with no drug and no growth factor, every state is
  stationary to `0.00e+00` over 42 simulated days. Any drift would mean the
  results are a mixture of pharmacology and numerical error.
* **Step-size convergence** — the cytotoxic PK is stiff (docetaxel's central
  compartment empties at ~9.7 /h) so it is pre-integrated on a 0.025 h grid and
  interpolated; the biology is stepped with a **time-local** bound of
  1.2/λ(t) recomputed every step. dt 0.25 → 0.05 moves the CAE nadir from
  0.06520 to 0.06509 and DSN from 6.1250 to 6.1271 d.
* **DSN as a state versus the output grid** — agreement to 8.9e-16 on CAE.

Four defects found by these checks and fixed rather than absorbed, all recorded
in the source at the line that caused them:

1. G-CSF's transit acceleration was applied to the outflow but not the
   proliferation term — so *more* G-CSF deepened the nadir.
2. The `reserve` term multiplied the self-renewal gain, so through the 1/γ
   exponent a 12% stem-cell loss **halved** the baseline ANC.
3. The thrombopoietin ceiling was applied to both branches of its own equation,
   making it a sign trap: above the ceiling both factors were negative and TPO
   grew faster the further past its limit it went.
4. The megakaryocyte pool was self-renewing with gain > 1 behind an 8-day
   feedback delay — a delayed positive loop that drove the platelet rebound to
   7118 ×10⁹/L before the count caught up.

A calibration search was also silently fooled once: `min()` over a list
containing NaN returns the first element, which reads as "no neutropenia at
all", so a diverging integration pushed the search to a docetaxel slope of 2713.
Non-finite states now raise.

---

## Model structure

**53 ODE compartments.** 4 cytotoxic PK slots (A 3-compartment, B/C/D
2-compartment) covering 8 parameterised agents across 13 regimens, with
threshold-driven PD for paclitaxel; G-CSF depot/central/peripheral/**effect**
with target-mediated clearance; trilaciclib; dexamethasone; levofloxacin; a
quiescent stem pool feeding three lineages each with its own transit chain;
marrow storage and marginated pools kept separate from circulating; mucosal
barrier → bacterial translocation → IL-6 → CRP → temperature; the DSN counter
and FN hazard; a tumour with sensitive and resistant clones; and neutrophil
function as a state distinct from count.

**Structural choices worth arguing with**

* **φ = 2⁻⁸** is the reciprocal of the mitotic amplification factor and sets an
  ANC **floor** that no cytotoxic can push through. It is the single most
  consequential number in the model: n_div = 5 puts the floor at 0.236 and no
  grade 4 nadir is reachable at any exposure; n_div = 9 gives DSN 5.33 d. The
  value 8 is the smallest that reaches grade 4 nadirs on the regimens that
  produce them clinically. It should be read as **chosen**, not measured.
* **(1 − E_drug) is not clamped at zero.** At the published Friberg slopes the
  term goes negative during the infusion and the proliferative pool is actively
  depleted rather than merely stopped. Clamping caps the achievable nadir far
  above what is observed.
* **Reserve acts on the influx, not the self-renewal gain**, so reduced marrow
  reserve barely moves the baseline count and markedly slows recovery — which is
  what a previously-treated patient actually looks like.
* **γ_P = 1.0 and γ_E = 0.6**, larger than the neutrophil's 0.170. Reusing so
  weak an exponent for the slower lineages makes their steady state absurdly
  sensitive, because it enters as 1/γ.
* **The inflammatory gains are set analytically**, from what florid
  Gram-negative sepsis looks like at the bacterial carrying capacity (IL-6 500
  pg/mL, CRP 250 mg/L, +3.0 °C). The first version used gains an order of
  magnitude too high and reported a core temperature of 83 °C.
* **`kkill_B` is set so that neutrophil killing exactly balances bacterial
  growth at ANC 0.5**, i.e. the clinical threshold is a *property* of the model
  rather than a number written into it.

---

## What the model is not for

* It cannot reproduce an **agent-specific nadir day** — the nadir day is a
  function of MTT alone, and carboplatin's is reported near day 21.
* It **cannot assign a regimen to a guideline risk band** (A1, miss 6). Use it
  to compare within a regimen, not across the band boundary.
* It **does not reproduce which lineage is dose-limiting** for the
  platelet-limited agents. Carboplatin and gemcitabine come out
  neutrophil-limited because a global scale right for a three-drug regimen is
  too large for a mild single agent. A version that got this right would need
  per-agent slopes fitted to per-agent data rather than to ratios.
* The **tumour layer is a price list**, not a survival model. The
  dose-intensity/survival association is observational and confounded.
* Absolute P(FN) values are **trial-population** numbers and read low against
  registry data (Truong 2016).

---

## Running it

```bash
# Full numerical report (~2 min, standard library only)
python3 cin_reference_impl.py

# One analysis, or refit the calibration from the anchors
python3 cin_reference_impl.py --only A7
python3 cin_reference_impl.py --recalibrate

# Re-render the map
dot -Tsvg cin_qsp_model.dot -o cin_qsp_model.svg
dot -Tpng -Gdpi=150 cin_qsp_model.dot -o cin_qsp_model.png
```

```r
library(mrgsolve)
mod <- mread("cin_mrgsolve_model.R")
# agent/regimen tables and 18 prebuilt scenarios are in the usage block
# at the foot of the model file
shiny::runApp("cin_shiny_app.R")
```

---

## Disclaimer

Research, education and hypothesis generation only. Not validated for clinical
decisions, prescribing, or regulatory submission. Parameters are illustrative
approximations; the six named misses above are the model's own account of where
it disagrees with the literature.
