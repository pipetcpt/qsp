# 승모판 역류증 (Mitral Regurgitation) — QSP 모델

**One regurgitant volume, five denominators.**

A regurgitant orifice produces exactly one number that echocardiography reports —
the regurgitant volume — and that number is meaningless until it has been divided
by something. Mitral regurgitation is therefore not one disease but a family of
diseases indexed by the denominator:

| | denominator | what it decides |
|---|---|---|
| **1** | operating LA compliance `C_op` | **congestion.** A virgin atrium sits high on a steep pressure–volume curve, so 60 mL it has never seen produces a giant v wave and floods the lung. A remodelled atrium has shifted that curve right and flattened it, so the identical 60 mL is nearly silent. |
| **2** | total stroke volume | **the grade.** In a low-output patient the same leak is a larger fraction, so reported severity climbs as the ventricle fails. |
| **3** | LV end-diastolic volume | **valve or ventricle** — and therefore whether closing the valve can help at all. This is the COAPT / MITRA-FR axis. |
| **4** | the afterload the leak *removes* | **hidden contractility.** EF is computed against a stroke volume containing the leak, which is ejected into a low-pressure sink. |
| **5** | `k_PISA` | **the measurement's own bias.** PISA assumes a hemisphere; the functional orifice is a crescent. |

Every controversy in this disease — why two trials of the same device reached
opposite conclusions, why ejection fraction *falls* after a successful operation,
why annuloplasty recurs, why a hugely dilated ventricle sometimes has no
regurgitation at all — is a disagreement about which denominator applies. The
model is built so that these come out as arithmetic rather than as assertions.

43 ODEs. The beat-level circulation is not integrated but solved algebraically at
every derivative evaluation by two nested monotone root-finds, so a closed-loop
elastance circulation is coupled to growth, remodelling, neurohormonal and device
dynamics on a time base of days. Every equation was implemented twice — once in
`mr_mrgsolve_model.R` and once, independently, in Python/scipy — and the
cross-check found and fixed **eight** defects, listed at the end.

---

## Deliverables

| file | what it is |
|---|---|
| [`mr_qsp_model.dot`](mr_qsp_model.dot) · [SVG](mr_qsp_model.svg) · [PNG](mr_qsp_model.png) | mechanistic map: 187 nodes, 18 clusters, 305 edges |
| [`mr_mrgsolve_model.R`](mr_mrgsolve_model.R) | 43-ODE mrgsolve model, 17 annotated scenarios |
| [`mr_shiny_app.R`](mr_shiny_app.R) | 12-tab dashboard, one tab per denominator |
| [`mr_references.md`](mr_references.md) | 116 references, every PMID resolved against PubMed programmatically |
| [`mr_python_reference.py`](mr_python_reference.py) | the independent implementation — the source of every number below |
| [`mr_reference_output.txt`](mr_reference_output.txt) | its verbatim output |

---

## What is fitted, and what is predicted

This matters more than usual here, because the model's headline claim is about
two randomised trials.

**Fitted to outcome data: two numbers.** The baseline hazard for heart-failure
hospitalisation and the baseline hazard for death, both taken from the **COAPT
control arm alone**. The hazard *slopes* are set a priori from published hazard
ratios (about 1.22 per 5 mmHg of filling pressure), not fitted. Because the
hazard accumulators feed back into nothing else, the baseline hazards enter
strictly linearly and can be recovered exactly rather than by search.

**Fitted to baseline imaging: three knobs per virtual patient.** Infarct size,
the acute dilation it causes, and the leaflet area the patient started with —
targeted at the reported LVEDV, LVEF and PISA-EROA. Each patient is then built by
*forward simulation* of three years of ischaemic cardiomyopathy; atrial size and
stiffness, annular area, fibrosis, pulmonary vascular tone, filling pressure and
cardiac output are consequences, not targets.

**Everything else is prediction**, including all four trial arms except the one
used for calibration.

Baseline physiology is not fitted at all. The healthy subject is an **exact fixed
point of all 43 states** (largest ten-year drift 3.1 × 10⁻⁷), and its
haemodynamics are resting values:

| | model | textbook |
|---|---|---|
| LV end-diastolic volume | 114.0 mL | ~114 |
| LVEF | 0.620 | 0.62 |
| cardiac output | 4.94 L/min | 5.0 |
| mean arterial pressure | 95.1 mmHg | 93 |
| mean LA pressure | 8.13 mmHg | 8 |
| mean PA pressure | 14.1 mmHg | 14 |
| CVP | 4.76 mmHg | 4.5 |
| LA volume | 53.0 mL | 55 |
| end-diastolic wall stress | 12.7 mmHg | ~12 |
| end-systolic wall stress | 81.8 mmHg | ~80 |
| wall thickness | 0.94 cm | 0.95 |
| PVR | 1.20 WU | 1.2 |

---

## Results

### 1. Denominator one: the same regurgitant volume is two different diseases

A controlled experiment. Ventricle, orifice, blood volume, contractility and
vasculature are identical; only the atrium differs, and the orifice is adjusted so
the **regurgitant volume is matched at 60 mL**. The numerator is held constant and
only the denominator moves.

| | acute (virgin atrium) | chronic (remodelled atrium) |
|---|---|---|
| regurgitant volume | **60.0 mL** | **60.0 mL** |
| LA volume | 62 mL | 277 mL |
| **operating LA compliance** | **0.80 mL/mmHg** | **2.74 mL/mmHg** |
| v wave | **35.5 mmHg** | **10.9 mmHg** |
| effective wedge | **34.7 mmHg** | **18.9 mmHg** |
| alveolar oedema threshold | 20.0 | 20.0 |
| verdict | **pulmonary oedema** | compensated |

A 3.4-fold difference in one compliance turns the identical leak into flash
pulmonary oedema in one patient and a wedge pressure below the oedema threshold in
the other.

The model also charges for the compensation, which the textbook account tends to
skip. The dilated atrium is a volume reservoir, so it lowers mean filling pressure
(17.0 → 13.4 mmHg) *and* sequesters preload the ventricle no longer has: cardiac
index is slightly **lower** in the chronic arm (1.56 against 1.78) at identical
regurgitation. Atrial remodelling does not abolish the cost of the leak; it moves
that cost from pressure to flow — which is why the chronic patient is breathless
on exertion rather than drowning at rest.

### 2. Denominator four: the operative threshold of EF 60% is *derived*

Take a compensated severe primary MR and abolish the orifice **instantaneously**,
with contractility, chamber size, stiffness and blood volume unchanged. Nothing
biological happens in that step — only the parallel low-impedance path is removed.

| | pre-op | post-op |
|---|---|---|
| EROA | 0.418 cm² | 0 |
| regurgitant fraction | 0.479 | 0 |
| total stroke volume | 102.8 mL | 71.1 mL |
| **forward** stroke volume | 53.5 mL | **71.1 mL** |
| **LVEF** | **0.619** | **0.490** |
| end-systolic pressure | 93.6 mmHg | 124.2 mmHg |
| cardiac index | 2.11 | **2.80** |
| E_es (unchanged) | 2.877 | 2.877 |

Ejection fraction falls 13 points with contractility untouched, because the
denominator of EF contained a stroke volume that was being ejected into an atrium
at 21 mmHg instead of an aorta at 94 mmHg. Sweeping contractility and asking where
post-operative EF crosses 0.50:

> **post-operative EF reaches 0.50 exactly when pre-operative EF is 0.627.**

The model was never shown a guideline. That number falls out of the loading
arithmetic, and it sits on top of the 0.60 the guidelines specify. The threshold is
not a convention: it is where the afterload the leak had been removing stops being
affordable.

And note the column the debate usually omits — forward cardiac index **rises** at
every level of contractility. Abolishing the leak makes the ejection fraction look
worse and the patient better, and those are the same fact seen through different
denominators.

### 3. Denominator five: do the reported trial numbers close?

Take each trial's reported LVEDV, LVEF and PISA-EROA at face value; compute the
regurgitant volume that orifice implies at the patient's own pressures; subtract it
from the total stroke volume implied by LVEDV × LVEF; read off the cardiac index
the patient must have been living at.

| k_PISA | COAPT: EROA_true → CI | MITRA-FR: EROA_true → CI |
|---|---|---|
| 1.00 (reported taken as true) | 0.410 cm² → **0.71** | 0.310 cm² → 1.76 |
| 1.20 | 0.342 → 1.00 | 0.258 → 1.99 |
| 1.45 | 0.283 → 1.24 | 0.214 → 2.20 |
| 1.80 | 0.228 → 1.47 | 0.172 → 2.38 |

Read as a constraint diagram, not an accusation. Taken at face value the COAPT
triplet implies a cardiac index of 0.71 L/min/m² — which no ambulatory outpatient
sustains. The arithmetic only begins to close if a substantial part of the reported
orifice area is measurement inflation, which is exactly what the crescentic
geometry of a functional orifice predicts. The asymmetry is the point: the same
0.4 cm² threshold is doing very different work in the two cohorts, so **the
proportionality debate is being conducted in a currency whose exchange rate differs
between the two morphologies being compared.**

### 4. Denominator three: one calibration, and what it does and does not predict

Both device arms use the **same** device effect — 68% of the orifice abolished, the
same added mitral resistance — so any difference between trials comes only from the
denominator the device acts against.

| | | HF hospitalisation | | death | |
|---|---|---|---|---|---|
| | arm | model | observed | model | observed |
| COAPT | control | 0.679 | 0.679 | 0.461 | 0.461 · *calibrated* |
| COAPT | device | 0.180 | 0.358 | 0.178 | 0.291 |
| MITRA-FR | control | 0.276 | 0.474 | 0.156 | 0.224 |
| MITRA-FR | device | 0.147 | 0.487 | 0.088 | 0.243 |

**What is right — the ordering and its mechanism.** The identical procedure buys
the COAPT patient 8.5% more forward output and the MITRA-FR patient 2.2%, a
**3.8-fold difference**, purely because the orifice being closed is larger relative
to the ventricle behind it:

| | Δ EROA | Δ RVol | Δ wedge | Δ CI | Δ CI (%) |
|---|---|---|---|---|---|
| COAPT | −0.190 cm² | −19.1 mL | −5.96 mmHg | +0.120 | **+8.5%** |
| MITRA-FR | −0.153 cm² | −18.6 mL | −3.68 mmHg | +0.039 | **+2.2%** |

Two quantities that were never fitted also come out right: the post-device mean
mitral gradient (**3.2 mmHg** against a reported 3–4) and the COAPT LV
end-diastolic volume index (**106** against a reported 101 mL/m²). The direction of
every treatment effect is correct.

**What is wrong, and reported as a miss rather than fitted away.** The model
predicts a real benefit in MITRA-FR (0.276 → 0.147) where the trial found none.
The haemodynamic separation the model can generate is 3.8-fold; the separation the
trials show is a change of *sign*. An additive prognostic term cannot repair this —
additive terms cancel in a hazard ratio — so the gap is structural. Two honest
candidates, which the model cannot distinguish:

- the MITRA-FR device arm did not achieve what is assumed here (that trial required
  no core-laboratory acute success and reported substantially more residual
  regurgitation); or
- a large part of that cohort's risk lived in a compartment this model has no
  variable for — infarct burden, arrhythmia, renal and skeletal muscle disease — in
  which case no valve model of any sophistication will predict their outcome from
  valve geometry.

Against candidate (i): even with **no orifice reduction at all**, the model sits
below the observed MITRA-FR event rate, which points at (ii).

The fitted knobs say the same thing mechanistically. The COAPT patient's leaflet
area came out at 1.03× native, the MITRA-FR patient's at **1.35×**: the COAPT valve
was diseased relative to its ventricle, the MITRA-FR valve was not.

### 5. The vortex is entirely valve-mediated

Dilation enters its own derivative twice — directly (a bigger radius means more
wall stress) and through the valve (dilation displaces the papillary muscles,
raises coaptation demand, opens the orifice, adds volume load). Splitting
`G = d(dV0d/dt)/dV0d` numerically:

At every operating point where the total gain is positive, it is positive
**because of the valve**. The direct pathway turns *negative* as the chamber grows,
because mass tracks cavity size and filling pressure falls as the chamber becomes
more compliant. On its own, dilation is self-limiting. This is why abolishing the
orifice removes the engine of further dilation even in a patient whose
regurgitation is not what is killing them.

The gain also crosses zero — at an end-diastolic volume of about 286 mL for the
MITRA-FR-like ventricle — not because the patient is better but because the growth
law has saturated.

Stated plainly because the opposite is the tempting conclusion: these numbers do
**not** separate the two trials. The valve supplies essentially all the positive
gain in both. The separation lives in the *magnitude* of the haemodynamic gain a
procedure buys, not in the sign of this derivative.

### 6. The speed of dilation matters, not only its size — but two clocks oppose

Two ventricles are forced to exactly the same final chamber size, one over about
two months and one over about four years. Two time-dependent adaptations are in
play, so the experiment runs twice.

**(a) annulus clamped** — coaptation demand at the target size is then identical
(11.995 cm² in both), and the only difference is how much leaflet the valve had
time to grow:

| dilation time | leaflet area | reserve | EROA |
|---|---|---|---|
| fast (~2 months) | 10.093 cm² | −1.901 | **0.567 cm²** |
| slow (~4 years) | 10.329 cm² | −1.666 | **0.498 cm²** |

**(b) annulus free** — the slow arm has had years of annular dilation, which more
than cancels its leaflet advantage and *reverses* the ordering (0.000 vs
0.250 cm²).

So the speed argument holds only while the annulus is not itself remodelling —
which is the situation after an acute infarct, and not the situation in
long-standing atrial-functional disease. The model states this as a prediction
rather than hiding the arm that disagrees with the framing.

### 7. Heart rate: the per-beat and per-minute numbers disagree

Clamping the rate and changing nothing else, in severe MR:

| HR | RVol per beat | regurgitant L/min | forward SV | forward CO | wedge |
|---|---|---|---|---|---|
| 45 | 45.1 mL | 2.03 | 66.2 mL | 2.98 | 33.1 |
| 70 | 49.0 mL | 3.43 | 55.1 mL | 3.86 | 34.3 |
| 110 | 48.2 mL | 5.30 | 45.1 mL | 4.96 | 32.6 |

The regurgitant volume **per beat** is almost flat and does **not** fall as the
rate rises: two effects cancel, because the regurgitant period shortens but
`Ea = c·R_sys/T` means a shorter cycle raises end-systolic pressure and with it the
transmitral gradient. So slowing the heart does not buy the per-beat reduction the
geometric argument promises. Meanwhile regurgitant flow per minute rises steeply
and forward stroke volume falls — yet forward *output* rises, and the wedge is
unchanged.

**Model limit, stated rather than buried:** forward output rises monotonically to
the top of the swept range, so the model has no interior optimum and would
recommend tachycardia without bound. It contains no ischaemic penalty, no
atrioventricular optimisation and no arrhythmic cost. The robust content is the
per-beat/per-minute divergence, not the location of an optimum.

### 8. The second barrier: the same operation at different times

A patient is left on medical therapy, then has an identical, perfectly successful
procedure and twelve months of observation. The procedure never changes.

| deferral | PVR before | irreversible PVR | wedge at 12 mo | CI at 12 mo |
|---|---|---|---|---|
| 0 y | 1.55 WU | 0.22 | 20.4 | **1.47** |
| 2 y | 4.00 | 1.64 | 28.9 | 1.11 |
| 4 y | 5.61 | 3.34 | 28.7 | 0.88 |
| 6 y | 6.53 | 4.71 | 26.7 | **0.78** |

The irreversible component is a ratchet: it only goes up. An identical procedure
therefore buys progressively less, and past a point the lung has become the disease
the valve used to be. Timing is not a scheduling detail; it is part of the
intervention. Note that contractility *improves* slightly across the deferral rows
(medical therapy is working) — the deterioration is pulmonary vascular, not
contractile.

### 9. Afterload reduction is anti-regurgitant, in two separable parts

The leak and the aorta are in **parallel**, so lowering systemic impedance
redistributes flow forward with no change to the valve, the ventricle or the volume
status. Nitroprusside, with the orifice *identical* at 0.280 cm²:

| | R_sys | EROA | RF | forward SV | CI |
|---|---|---|---|---|---|
| baseline | 1.598 | 0.280 | 0.451 | 35.4 mL | 1.42 |
| nitroprusside | 1.061 | **0.280** | **0.379** | **43.1 mL** | **1.73** |

No remodelling, no time, no biology — purely a change in the competing impedance.

Chronic therapy adds a slower structural effect. Separating the two requires care:
zeroing a dosing *rate* does not remove a drug, because its concentration is a
state. So the decomposition is done against an untreated control arm, and by
stripping the drug compartments out of the treated arm while leaving everything the
drug has *done*:

| | EROA | RF | CI | wedge |
|---|---|---|---|---|
| day 0 | 0.280 | 0.451 | 1.42 | 26.6 |
| day 180, **untreated** | 0.492 | 0.589 | 1.28 | 56.1 |
| day 180 treated, drug stripped out | 0.437 | 0.595 | 1.38 | 37.0 |
| day 180 treated, drug present | 0.437 | 0.574 | 1.31 | 33.8 |

Both arms progress — this is a severe secondary regurgitation on the way up — so
the therapeutic effect is the difference from the untreated arm:

- structural component: **+0.006** regurgitant fraction
- instantaneous impedance component: **−0.021**

Note what that says about the endpoint. Therapy did a great deal of structural
good — it held EROA at 0.437 instead of 0.492 and the wedge at 33.8 instead of
56.1 mmHg — and almost none of that appears in the regurgitant *fraction*, because
RF has total stroke volume in its denominator and both move together. A trial
powered on regurgitant fraction is measuring mostly the impedance component, which
appears with the first dose and vanishes on withdrawal.

### 10. Annuloplasty: how far must a ring be undersized to hold?

Five procedures, identical patient, three years. The ring cuts the **annular** arm
of the loop and leaves the **tethering** arm intact, so whether it holds is a race
between how small the annulus was made and how much the ventricle goes on to
dilate. That race has a threshold:

| ring | vs native | EROA at 3 y | RF at 3 y | wedge | CI |
|---|---|---|---|---|---|
| 7.64 cm² | 1.00× | **0.427** | 0.588 | 27.1 | 1.31 |
| 7.25 | 0.95× | 0.198 | 0.370 | 16.5 | 1.53 |
| 6.87 | 0.90× | **0.002** | 0.006 | 8.8 | 1.70 |
| 5.50 | 0.72× | 0.000 | 0.000 | 8.7 | 1.70 |

A ring that merely *stabilises* the annulus at its current size does not hold: the
tethering arm re-opens the valve underneath it and the orifice is back to
0.427 cm² — recurrent severe regurgitation — within three years. Undersize by about
10% and the same operation holds completely. The threshold sits between 0.90 and
0.95 of native annular area, and it comes out of the geometry rather than being
assumed. That is the quantitative content of the surgical argument about
undersizing.

**Where the model is optimistic, and it should be said:** real series report
recurrence in a large minority even after aggressive undersizing, whereas here a
0.72× ring holds perfectly. The mechanism — recurrence concentrated in ventricles
that go on dilating — is present; the calibration of how hard tethering can pull
against a fixed annulus is evidently too forgiving.

For comparison, over the same three years: medical therapy alone lets EROA reach
1.479 cm² with a wedge of 46.3 mmHg and cardiac index 0.74; edge-to-edge repair
holds it at 0.318 cm², wedge 25.5, CI 1.37, at the price of a 4.6 mmHg transmitral
gradient; replacement holds at 0.022 cm² with a 1.2 mmHg gradient, and if the
chordae are divided the model charges a further 10% of contractility for the lost
ventricular scaffold.

---

## Verification

The model was implemented twice, independently, and the numbers above come from the
Python implementation.

- The healthy subject is an **exact fixed point** of all 43 states — no non-zero
  derivative above 10⁻¹², and a largest ten-year drift of 3.1 × 10⁻⁷.
- **Three independent integrators** (LSODA, RK45, DOP853) agree on a three-year
  trajectory to six significant figures. LSODA is 14× faster, confirming the
  system is genuinely stiff and that the fast integrators are not silently
  smoothing it.
- Tolerance sweep from rtol 10⁻⁵ to 10⁻⁸ changes results in the fifth significant
  figure.
- Both algebraic roots solved inside every derivative evaluation are **strictly
  monotone** — the end-systolic pressure residual decreasing in `P_es`, the
  volume-conservation residual increasing in `P_LA` — so each has a unique root and
  bracketed Brent iteration cannot converge to the wrong one.
- Every one of the 116 references was resolved against PubMed programmatically;
  where a query's top hit was not the intended paper it was re-labelled rather than
  silently mis-cited, and two junk hits were dropped.

### Defects the cross-check found and fixed

1. **A fibrosis→stiffness→low-output→RAAS→fibrosis loop with gain above one.** Even
   a trivial orifice (EROA 0.10 cm²) collapsed the ventricle to a 76 mL cavity
   within a year. Fibrotic stiffening was unbounded and driven by *diastolic*
   stress; it is now bounded and driven by *systolic* stress — which also makes the
   model correctly reproduce chronic MR as a dilating, thin-walled, low-fibrosis
   phenotype.
2. **The hypertrophy law had the wrong sign for volume overload.** Mass was driven
   by end-systolic stress, which *falls* in MR, so the ventricle atrophied. Mass now
   tracks cavity size (eccentric) with a separate concentric term.
3. **A linear atrial compliance could not produce a giant v wave**, which is the
   central observable of acute MR. Replaced with an exponential atrial
   pressure–volume law.
4. **That exponential then produced v waves above 500 mmHg** once atrial fibrosis
   stiffened the curve, because the excursion was integrated along the exponential.
   The v wave now comes from the *incremental* compliance at the operating point.
5. **Pressure natriuresis was ~17× too weak**, so any low-output state drove blood
   volume to 80 L.
6. **Diuretic natriuresis had no floor**, so guideline-directed therapy made a
   virtual patient hypovolaemic (cardiac index 1.94, wedge 3.2) instead of
   decongested.
7. **Two remodelling laws let organs shrink below their native size** — leaflets
   atrophied below native area and the atrium below its own volume — silently
   destroying coaptation reserve so that even mild MR appeared to progress. This one
   was caught only because it broke the healthy fixed point.
8. **A tanh cap on the hazard rate saturated**, making treated and untreated arms
   numerically identical and destroying the linearity that lets the baseline hazards
   be recovered exactly. The exponent is now clamped instead of the rate.

Verification also **refuted two of the author's framings**. Result 6 was written
expecting rate-limited leaflet growth to dominate; with a free annulus the ordering
reverses, and the disagreeing arm is reported. Result 5's loop-gain decomposition
was expected to separate proportionate from disproportionate MR; it does not, and
that is now stated explicitly.

### Honest limitation on the R model

**No R toolchain was available in the environment where this was built, so
`mr_mrgsolve_model.R` has not been executed.** The Python implementation is the
executed reference; the R file mirrors it equation for equation, but numerical
agreement between the two has *not* been demonstrated. This is a weaker claim than
the other models in this library make, and it is stated here rather than implied
away.

---

## Reproducing

```bash
# regenerate the map
dot -Tsvg mr_qsp_model.dot -o mr_qsp_model.svg
dot -Tpng -Gdpi=150 mr_qsp_model.dot -o mr_qsp_model.png

# reproduce every number in this README
python3 mr_python_reference.py > mr_reference_output.txt
```

```r
library(mrgsolve)
mod <- mread("mr_mrgsolve_model.R")

# healthy baseline: an exact fixed point
mod %>% mrgsim(end = 3650, delta = 10) %>% plot(EF + CI + Ppcw + EDV ~ time)

# acute papillary muscle rupture onto a virgin atrium
mod %>% init(EROApri = 0.60) %>% mrgsim(end = 14, delta = 0.1) %>%
    plot(RVol + vwave + ClaOp + Ppcw + EF + EFfwd ~ time)

shiny::runApp("mr_shiny_app.R")    # 12 tabs, one per denominator
```

---

> **면책 조항.** 본 모델은 교육 및 연구 목적의 QSP 모델입니다. 공개 문헌을 바탕으로
> 구성되었으나 독립적으로 검증·인증되지 않았으며, 실제 임상 의사결정, 처방, 또는
> 규제 제출에 직접 사용해서는 안 됩니다. 특히 승모판 중재술의 적응증과 시점에 관한
> 본 모델의 예측은 임상 근거가 아니라 가설이며, MITRA-FR 결과를 재현하지 못한다는
> 점(위 "Results 4" 참조)을 반드시 함께 읽어야 합니다.
