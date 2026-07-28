# IgE-mediated food allergy and anaphylaxis — QSP model

**IgE 매개 식품 알레르기 · 아나필락시스 · 유발용량(eliciting dose)을 하나의 계산값으로 다루는 모델**

| File | What it is |
|---|---|
| [`fana_qsp_model.dot`](fana_qsp_model.dot) · [`.svg`](fana_qsp_model.svg) · [`.png`](fana_qsp_model.png) | Mechanistic map — 20 modules, 219 labelled nodes |
| [`fana_mrgsolve_model.R`](fana_mrgsolve_model.R) | 35-ODE mrgsolve model, 12 scenarios, 16 diagnostics, virtual population |
| [`fana_shiny_app.R`](fana_shiny_app.R) | 9-tab interactive dashboard |
| [`fana_references.md`](fana_references.md) | 112 references + 5 out-of-scope pointers, every entry with a PubMed link |

Built and run under **mrgsolve 2.0.1 / R 4.3.3**.
`source("fana_mrgsolve_model.R"); run_diagnostics()` reproduces every number below.

---

## Why this disease needed a different kind of model

Most disease models have to invent their endpoint. Food allergy does not.
It has one number that the whole field agrees on:

> **the eliciting dose (ED)** — the milligrams of allergen protein that
> provoke an objective reaction.

Labelling law is written from its population distribution (VITAL reference
doses). Trials are powered on a single point of it ("tolerates a single
600 mg dose of peanut protein"). Patients organise their lives inside it.

So this model computes one thing and computes it properly. `find_ED()`
bisects on single-dose challenges until it finds the smallest dose that
produces a graded clinical reaction — the same operation a double-blind
food challenge performs, on the same observable. Nothing in the file
reports an ED that was assumed; every ED in every table below was
*measured off a simulated challenge*.

That choice forces the interesting question to the front. If the ED is
the output, **what is the input?**

---

## The engine: the specificity enters twice

FcεRI binds the **Fc** of IgE, not its paratope. The mast cell surface
therefore samples the serum IgE pool *without regard to specificity* — it
carries allergen-specific IgE at the same fraction

$$f = \frac{\text{sIgE}}{\text{total IgE}}$$

as the serum does. Degranulation requires one allergen molecule to
**bridge two adjacent receptors**. The probability that a randomly chosen
adjacent pair is specific-*and*-specific is $f^2$, so the density of
bridgeable pairs goes as

$$X \;\propto\; (\rho \cdot L \cdot f)^2 \cdot \frac{[A]}{K_A + [A]}$$

with $\rho$ = receptor density and $L$ = fractional occupancy. **The
exponent 2 is not fitted. It is combinatorics.** And it is the most
consequential number in the file, because it splits every drug in the
field into two classes:

| acts on | effect on the ED |
|---|---|
| **the surface** (ρ, L, f) — anti-IgE, anti-KIT | **two logs per log** |
| **the allergen** ([A]) — IgG4 from immunotherapy | **one log per log** |

Diagnostic **D3** checks the exponent numerically by sweeping sIgE and
regressing log ED on log surface index:

```
surface index (rel) :  0.250   0.500   1.000   2.000
ED (mg)             :  944.7   145.3    33.3     8.1
slope over the whole sweep  : -2.27      (theory: -2.00)
slope near the reference    : -2.13
```

The residual steepening is the allergen term saturating at the high end,
where the required dose pushes $[A]$ past $K_A$ — a dose region beyond
any food. Across the clinical range the law holds.

---

## Five structural commitments

**S1 — the threshold is set by a RATIO, not a titre.** Because $f$ is a
ratio, a patient with a big atopic background dilutes their own specific
IgE. Two patients with *identical* sIgE of 40 kU/L:

| sIgE | total IgE | f | model ED |
|---|---|---|---|
| 40 | 150 | 0.34 | **9.8 mg** |
| 40 | 300 | 0.17 | **33.3 mg** |
| 40 | 2000 | 0.025 | **6,944 mg** |

The third patient is sensitised and, for practical purposes, not
reactive. This is not a curiosity of the model: it is why the
sIgE/total-IgE **ratio** outperforms sIgE alone at predicting challenge
outcome (Gupta 2014; Sindher 2018), an observation with no explanation
under a "more IgE, more allergy" reading.

**S2 — receptor density is a slow state, not a constant.** Free IgE
stabilises FcεRI; remove the ligand and the receptor is degraded. So
anti-IgE has a **fast arm** (occupancy, days) and a **slow arm** (density,
months) that multiply *inside the square*:

```
wk  0  freeIgE 719.9 ng/mL   totalIgE  300 IU/mL   rho 0.931   SURFrel 1.000
wk  1  freeIgE   6.1          totalIgE  859        rho 0.791   SURFrel 0.213
wk  4  freeIgE  28.2          totalIgE 1846        rho 0.567   SURFrel 0.373
wk 16  freeIgE  20.7          totalIgE 2616        rho 0.348   SURFrel 0.200

ED at week  4 :   294 mg  (  8.8x baseline)
ED at week 16 : 2,601 mg  ( 78.1x baseline)
```

Free IgE is at target within **days**. The threshold is still climbing at
week 16. That gap is the receptor arm, and it is why every successful
anti-IgE food-allergy trial places its challenge at week 16–20 rather
than week 4. Nothing pharmacokinetic is happening in that window.

**S3 — immunotherapy and anti-IgE act at different points in series.**
IgG4 intercepts allergen in the interstitium *before* the surface
(dividing $[A]$: linear). Anti-IgE changes the surface (quadratic). Being
in series, they should compound:

| arm | model ED | fold |
|---|---|---|
| untreated | 33.3 mg | 1× |
| OIT 300 mg/d, 6 months | 496 mg | 14.9× |
| omalizumab, week 16–26 | 14,066 mg | 422× |
| **omalizumab + OIT** | **54,274 mg** | **1,629×** |
| *sum of the two singles* | | *437×* |
| *product of the two singles* | | *6,287×* |

Strongly **super-additive** but well short of the naive product — and the
shortfall is itself mechanistic. OIT raises specific IgE modestly and
total IgE with it, so the omalizumab dose (computed from the *baseline*
IgE, as it is in practice) becomes relatively under-matched; and the
allergen term saturates at the doses involved. The model's prospective
claim is the ordering and the super-additivity, not the number.

**S4 — the circulation has a reserve, so severity is a cliff.** About a
third of plasma volume can move into the interstitium before mean
arterial pressure changes at all. MAP is the product of a *graded* term
(vascular tone) and a *cliff* term (reserve spent):

```
dose     30 mg : PV deficit  6.9%   reserve spent 0.000   tone 1.06   MAP 93.1   grade 1
dose    100 mg : PV deficit 21.6%   reserve spent 0.027   tone 0.94   MAP 82.8   grade 3
dose    300 mg : PV deficit 30.2%   reserve spent 0.227   tone 0.87   MAP 69.1   grade 3
dose   1000 mg : PV deficit 32.8%   reserve spent 0.342   tone 0.85   MAP 58.3   grade 4
dose   3000 mg : PV deficit 33.9%   reserve spent 0.398   tone 0.85   MAP 53.2   grade 5
```

Watch the reserve column. The plasma deficit climbs steadily from 7% to
30% while the reserve column stays near zero — and then moves. Severity
grading scales look ordinal because a reserve is being spent, not because
a dial is being turned.

**S5 — the measured threshold is a random variable.** Repeat challenges
in one patient scatter by ~0.4–0.5 log10, and the variance is
cofactor-driven — acting on allergen **delivery**, not on immunology. The
model therefore refuses to let cofactors touch sIgE:

```
none             ED 33.3 mg  ( 0.00 log10)
exercise         ED 10.9 mg  (-0.48)
NSAID            ED  9.5 mg  (-0.54)
alcohol          ED 16.3 mg  (-0.30)
infection        ED 18.2 mg  (-0.26)
PPI              ED 15.1 mg  (-0.34)
exercise + NSAID ED  3.2 mg  (-1.02)
```

A patient's antibodies do not change between Tuesday and Wednesday. Their
threshold does.

---

## The informative negative: why dupilumab monotherapy failed

This is the result the model is proudest of, because the model was not
built to produce it.

Dupilumab blocks IL-4Rα, so it suppresses IgE class switching. Over
24 weeks the model gives:

```
sIgE       40.0  ->  10.6 IU/mL      (73% suppression)
total IgE   300  ->    75 IU/mL      (75% suppression)
f          0.133 -> 0.142            <-- barely moves
ED          33.3 ->  46.7 mg         ( 1.40x )
```

**It removes the numerator and the denominator at the same time.** The
term that gets squared is $f$, and $f$ is a ratio, so a drug that halves
both halves of the ratio does almost nothing to the threshold — no matter
how impressive the IgE suppression looks on a lab report. Dupilumab
monotherapy moved the peanut threshold in roughly 2 of 24 patients. The
model reproduces that as a structural consequence of S1, not as a fitted
weakness.

---

## The other negative: who should *not* be given anti-IgE alone

`omalizumab_dose()` implements the label table as the rule that underlies
it — the drug must supply IgE **binding sites** in stoichiometric excess,
about 0.008 mg/kg per IU/mL per fortnight — and it caps at the licensed
600 mg q2wk. Diagnostic **D16** holds the specific fraction fixed at
f = 0.125 so that all four patients start with essentially the same
threshold, and varies only how much total IgE the drug has to mop up:

| total IgE | dose | free IgE at wk16 | | ED | fold |
|---|---|---|---|---|---|
| 300 | 375 mg q4wk | 20.8 ng/mL | at target | 38 → 3,175 mg | **83×** |
| 1000 | 600 mg q2wk | 18.1 ng/mL | at target | 32 → 4,059 mg | **128×** |
| 2000 | 600 mg q2wk | 45.1 ng/mL | **missed** | 32 → 308 mg | **9.7×** |
| 4000 | 600 mg q2wk | 177.9 ng/mL | **missed** | 32 → 75 mg | **2.4×** |

The benefit does not taper. It **falls off the edge of the dosing table**,
because once free IgE stops reaching target the receptor arm never
engages at all and the quadratic term loses both of its factors at once.
This is a structural non-response, not a pharmacodynamic one, and it is
the model's most directly testable clinical prediction.

---

## Adrenaline: the argument is an integral

A 0.3 mg IM dose has a plasma half-life of ~2.5 minutes. It is a
15-minute window, and its *efficacy* does not decay with how long you
waited. What grows with delay is the deficit it has to repair:

```
protocol: 0.3 mg IM at the stated delay, repeated at +8 and +16 min

delay  2 min : MAP nadir 92.8   min below 65 mmHg  0.0   peak PV deficit 26.2%
delay 10 min : MAP nadir 84.6   min below 65 mmHg  0.0   peak PV deficit 25.1%
delay 30 min : MAP nadir 84.4   min below 65 mmHg  0.0   peak PV deficit 23.2%
no adrenaline: MAP nadir 62.1   min below 65 mmHg 38.1   peak PV deficit 33.5%
```

One caveat the model insists on: run the same sweep with a **single**
dose (`n_epi = 1`) and early adrenaline looks *worse* than late, because
with a 2.5-minute half-life one shot given at 2 minutes has washed out
before the leak peaks. That is a true statement about one injection and a
false statement about the guideline. The model is written with the
repeat-dosing protocol as the default for exactly this reason.

The antihistamine result is the mirror image. Cetirizine premedication:

```
dose    66 mg (2x threshold)
  no premedication   peak skin 0.69   MAP nadir 86.5   FEV1 nadir 81.3   grade 2
  cetirizine 10 mg   peak skin 0.26   MAP nadir 86.5   FEV1 nadir 81.3   grade 2
```

The skin score falls by more than half. **MAP and FEV1 do not move at
all.** The warning sign is removed and the lesion is left.

---

## Reproducing a trial as a *rate*, not a median

Trials report binary endpoints. A model that only produces a median
cannot be compared with one. So `vpop()` builds a virtual population,
applies the entry criterion, and reads the endpoint the way the trial
does — with the within-subject threshold variance of S5 applied
independently to each challenge:

```
population n = 160; meeting entry ED <= 100 mg: 75
entry ED  median 18.3 mg   IQR 9.3 - 55.3 mg
  observed entry (OUtMATCH): reacted at <=100 mg, median maximum
  tolerated dose 10 mg (i.e. eliciting dose ~30 mg)

placebo arm, second challenge, tolerating 600 mg :  2.7 %     observed:  7%
omalizumab wk16, tolerating 600 mg               : 54.5 %     observed: 67%
median ED  entry 18.3 mg  ->  omalizumab 775 mg
  observed  10 mg tolerated -> >1000 mg tolerated
```

The **placebo rate is not fitted**. It is what falls out of giving the
same patient two challenges when the threshold carries a 0.45 log10
within-subject standard deviation. A model with a deterministic threshold
cannot produce a placebo responder at all — which is a good reason not to
build one.

---

## What the model gets wrong

Stated here rather than buried.

1. **The omalizumab median shift is over-predicted.** 422× against an
   observed shift that is censored at ">1000 mg" but is probably nearer
   100×. Interestingly the *binary* endpoint is reproduced far better
   (54.5% vs 67%) than the median, because the endpoint is dominated by
   the low tail of the entry distribution rather than by the centre.
2. **The virtual population is not sensitive enough at entry** (median
   ED 18 mg vs ~30 mg observed), and the placebo rate consequently comes
   out at 2.7% rather than 7%. Both would move together if the
   within-subject SD were raised, and there was no attempt to tune them
   into agreement.
3. **Tryptase breaches the consensus rise criterion for essentially every
   objective reaction** in the model, whereas a large minority of real
   food-triggered anaphylaxis episodes do not. The reason is structural:
   there is one well-mixed mast cell compartment here, and the real
   explanation is that gut and skin mast cells are far from an
   antecubital vein. Fixing it needs spatial structure the model does
   not have.
4. **Basophil and mast cell are collapsed into one surface.** The model
   therefore cannot explain BAT/challenge discordance, because it has no
   mechanism for it. `RHO_FLOOR` is deliberately set from *skin mast
   cell* rather than *basophil* FcεRI downregulation — using the basophil
   number (a ~10–30× fall rather than ~3.6×) would over-predict anti-IgE
   efficacy by an order of magnitude, which is worth knowing about any
   model that uses basophils as its readout.
5. **The combination arm is the most aggressive claim in the file** and
   the one most worth testing when OUtMATCH stage 2 reports.
6. **`FINTACT` is the one free scaling parameter** (the intact-protein
   fraction crossing the gut is not measurable in humans) and it is
   confounded with `XL50`. Everything the model reports is a *ratio*, and
   ratios are insensitive to it — which is the point of building the
   model on thresholds rather than on concentrations.

---

## Model contents

**35 ODE compartments** — allergen (luminal, interstitial); mediators
(histamine, tryptase, PAF, cysLT, late-phase cells); granule content;
mast cell anergy; FcεRI density; plasma volume, interstitial shift and
the harm integral; bronchoconstriction, skin and GI states; Th2, memory
B, IgE plasma cells (specific and non-specific), specific and total IgE,
IgG4 plasma cells, IgG4, Treg; and PK for omalizumab (SC/central/
peripheral with target-mediated disposition), dupilumab, adrenaline
(IM and intranasal) and cetirizine.

**12 scenarios** — `sc_dbpcfc`, `sc_accidental`, `sc_omalizumab`,
`sc_oit`, `sc_combo`, `sc_dupilumab`, `sc_cofactor`, `sc_epi_timing`,
`sc_antihistamine`, `sc_high_ige`, `sc_natural`, `sc_biphasic`.

**16 diagnostics** — D1 baseline threshold · D2 the ratio · D3 the square
law · D4 the two arms of anti-IgE · D5 the total-IgE artefact · D6 the
dupilumab negative · D7 series architecture · D8 the adrenaline integral
· D9 cofactors · D10 antihistamine masking · D11 the MAP cliff · D12
tryptase · D13 immunotherapy on and off · D14 the age gate on remission ·
D15 the trial endpoint as a rate · D16 the anti-IgE non-responder.

**The Shiny app** re-derives the threshold by bisection every time a
slider moves. Nothing in it is a stored phenotype.

---

## Explicitly out of scope

FPIES · eosinophilic oesophagitis · alpha-gal syndrome · histamine
intolerance and scombroid · coeliac disease · MRGPRX2-mediated and
idiopathic anaphylaxis · clonal mast cell disease as a primary driver
(it enters only as `MCBURDEN`). These are drawn as a greyed-out cluster
on the map so that their absence is visible rather than implied, and
[`fana_references.md`](fana_references.md) points to the right literature
for each.

**This is a research model. It is not a clinical decision tool, it has
not been validated for any patient-level use, and no part of it should be
used to estimate a real person's threshold or to plan real exposure.**
