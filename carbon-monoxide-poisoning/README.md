# 일산화탄소 중독 (Carbon Monoxide Poisoning) — QSP 모델

**Two occupancies, two clocks.**

Carbon monoxide creates two distinct occupancies in the body, and they run on
different clocks:

| | pool | time constant | measured? |
|---|---|---|---|
| **fast** | carboxyhaemoglobin | t½ 313 min on air, **72 min** on a mask | **yes — this is the only one anyone measures** |
| **slow** | CO on myoglobin and on the reduced a3 haem of cytochrome c oxidase | off-rate τ ≈ **455 min**, and it loads *more* readily as tissue PO₂ falls | no |

Every clinical paradox in this disease is a consequence of monitoring and
treating the fast pool while the injury is being written by the slow one: a
reassuring pulse oximeter, a normal carboxyhaemoglobin in a patient who will
deteriorate three weeks later, and two randomised trials of hyperbaric oxygen
that reached opposite conclusions.

The model is 45 ODEs coupling a Coburn–Forster–Kane pulmonary mass balance to a
bistable autoimmune demyelination switch. Every equation was implemented twice —
once in `co_mrgsolve_model.R` and once, independently, in Python/scipy — and the
cross-check found and fixed six defects, listed at the end.

---

## Deliverables

| file | what it is |
|---|---|
| [`co_qsp_model.dot`](co_qsp_model.dot) · [SVG](co_qsp_model.svg) · [PNG](co_qsp_model.png) | mechanistic map: 150 nodes, 17 clusters, 239 edges |
| [`co_mrgsolve_model.R`](co_mrgsolve_model.R) | 45-ODE mrgsolve model, 17 annotated scenarios |
| [`co_shiny_app.R`](co_shiny_app.R) | 12-tab interactive dashboard |
| [`co_references.md`](co_references.md) | 132 references, every PMID resolved against PubMed programmatically |
| [`co_reference_output.txt`](co_reference_output.txt) | verbatim output of the independent Python implementation — the source of every number below |
| [`co_python_reference.py`](co_python_reference.py) | that independent implementation |

---

## Results

### 1. The half-life is a property of ventilation, and hyperbaric pressure runs into a wall the equation itself builds

The CFK transfer resistance decomposes as

```
B = 1/D_L,CO + P_L/V̇_A = 0.04000 + 0.16976 = 0.20976 mmHg·min/mL
```

so **81% of the resistance to carbon monoxide elimination is ventilation** and
only 19% is membrane diffusion. With constants taken from resting physiology and
not fitted (V̇_A 4.2 L/min, D_Lco 25 mL/min/mmHg, Haldane M = 245, V_b 5.5 L),
the model reproduces the measured half-lives:

| delivery | model | observed |
|---|---|---|
| room air | **312.9 min** | 320 min |
| non-rebreather mask | **71.7 min** | 74 min |
| nasal cannula 6 L/min | 124.1 min | — |
| intubated, FiO₂ 1.0 | 55.2 min | — |

One quantity in the whole gas-exchange block is fitted (`shunt` = 0.72, folding
venous admixture and mask leak into the mean capillary PO₂), and it was set once
on the 74-minute figure.

**Then pressure stops working.** Raising chamber pressure raises the driving
force `P_c,O₂` *and* the resistance term `P_L/V̇_A` in exact proportion, so

```
as ATA → ∞ :   P_c,O₂/B → V̇_A     and     k → V̇_A / (M·[O₂Hb]·V_b)
```

The rate constant has a **ceiling**, and therefore the half-life has a **floor**:

```
t½_min = ln2 · M · [O₂Hb] · V_b / V̇_A = 31.3 min
```

which depends only on alveolar ventilation, haemoglobin and blood volume — not
on pressure at all. The model duly predicts 42.4 min at 1.5 ATA, 38.6 at 3.0,
and still 35.6 at 20 ATA. The floor is real inside the equation and it is where
CFK breaks: **the observed half-life at 2.5–3 ATA is about 20 min, which is 36%
below a floor the equation says cannot be crossed.** That is a falsification, not
a bad parameter, and it is reported rather than fitted away — either an
extrapulmonary elimination route exists or the Haldane ratio is itself
pressure-dependent.

The corollary is practical: inside a chamber the *only* remaining lever is
ventilation, because the floor scales with it.

| V̇_A (L/min) | 4.2 | 6 | 8 | 10 | 12 |
|---|---|---|---|---|---|
| floor t½ (min) | 31.3 | 21.9 | 16.4 | 13.1 | 11.0 |

### 2. Hyperbaric oxygen works before it has removed any carbon monoxide at all

Dissolved oxygen is `0.003 × PaO₂`, and it is unaffected by how much
haemoglobin the CO has taken:

| pressure | PaO₂ | dissolved O₂ | vs whole-body A–V difference | vs cerebral |
|---|---|---|---|---|
| 1.0 ATA | 653 | 1.96 mL/dL | 0.43× | 0.31× |
| 2.0 ATA | 1413 | 4.24 | 0.92× | 0.67× |
| 2.4 ATA | 1717 | 5.15 | 1.12× | 0.82× |
| 3.0 ATA | 2173 | **6.52** | **1.42×** | **1.03×** |

Dissolved oxygen alone covers the entire cerebral demand at **2.90 ATA**. Above
that pressure the brain is fully oxygenated while the haemoglobin is still
saturated with carbon monoxide: **the therapeutic effect precedes the
pharmacokinetics.** This is why the argument about how fast a chamber clears
COHb — the argument the floor in result 1 says it cannot win — was the wrong
argument to be having.

### 3. Carbon monoxide is worse than the anaemia it imitates, and the model says by how much

CO does two things with one ligand: it removes carrying capacity, and it
left-shifts the curve for whatever capacity is left (P50_eff = P50·(1 − 0.75·F)).
Matching arterial oxygen *content* between CO and anaemia isolates the second:

| COHb | CaO₂ | equivalent anaemia | tissue PO₂ with CO | tissue PO₂ of that anaemia | deficit |
|---|---|---|---|---|---|
| 20% | 15.88 | 12.00 g/dL | 29.6 | 34.8 | 5.2 |
| 30% | 13.93 | 10.50 g/dL | 25.1 | 32.5 | 7.3 |
| **40%** | **11.98** | **9.00 g/dL** | **20.8** | **29.7** | **8.9 mmHg** |
| 50% | 10.03 | 7.50 g/dL | 16.5 | 26.4 | 9.9 |

At COHb 40% the arterial oxygen content equals that of a haemoglobin of 9.0
g/dL — a level nobody transfuses urgently — yet the tissue PO₂ is **8.9 mmHg
lower**, because the residual haemoglobin will not release what it holds. The
same COHb percentage is a different disease at a different haemoglobin, which is
why anaemia is the model's most dangerous comorbidity.

### 4. The two clocks, quantified

For a severe exposure (2298 ppm × 60 min, calibrated to a peak COHb of 40%)
treated with a mask, the pools separate cleanly:

| pool | t(50% of peak) | t(10% of peak) |
|---|---|---|
| COHb — *the measured one* | 173 min | 389 min |
| brain tissue CO | 173 min | 389 min |
| cardiac myoglobin-CO | 259 min | 562 min |
| skeletal muscle CO | 648 min | 6093 min |
| **brain cytochrome c oxidase** | **778 min** | **7519 min** |

COHb reaches the "normal" 5% at **346 min (5.8 h)**. At that moment brain
cytochrome c oxidase is still **19.9% CO-inhibited**, and it does not reach 10%
of its own peak until **7519 min (5.2 days)** — 19× later than COHb does.
(It never reaches 5% of peak, because the endogenous CO occupancy of the enzyme
is itself about 1%, which is roughly 5% of the poisoned peak: the floor is
physiological, not a numerical artefact.) **The instrument that decides
disposition is measuring the fast pool.** That is the model's account of the
long-standing observation that carboxyhaemoglobin does not predict outcome — it
is not that COHb is a noisy marker of severity, it is that it is a marker of the
wrong compartment.

Two structural features generate this, and neither was imposed:

- CO binds **only the reduced a3 haem**, so the drive is `∝ 1/(1 + PtO₂/K_O)`.
  Oxygen is a *competitive protector*, which means the tissue pool loads more
  avidly exactly when tissue oxygen is low — a vicious circle with no analogue
  in the blood compartment.
- Skeletal muscle, with its large myoglobin capacity and τ = 320 min, is a
  reservoir that back-diffuses. It is why COHb plateaus rather than continuing to
  fall when oxygen is withdrawn.

### 5. Delayed neurological sequelae are a bistable switch with a threshold you can write down

Peroxidation charge-modifies myelin basic protein into an antigen. The
autoreactive clone expands if and only if its proliferation beats its
contraction:

```
T_prol · H(MBPad) > T_death
  ⟹  H > T_death/T_prol = 0.4750
  ⟹  MBPad_crit = θ·((H/(1−H))^(1/n)) = 0.450 · (0.4750/0.5250)^(1/6) = 0.44256
```

and because demyelination liberates further antigen (epitope spreading), above
the threshold the loop **latches**. Below it, the transient adduct decays and the
clone contracts to its naive frequency.

The dose-response within one patient is consequently a **step**, and the model
puts it between COHb 26% and 33% for a 1-hour exposure treated promptly with a
mask:

| ambient CO | peak COHb | peak MBPad | ÷ threshold | clone | demyelination | cognitive | DNS |
|---|---|---|---|---|---|---|---|
| 700 ppm | 12.7% | 0.024 | 0.05× | 0.004 | 0.000 | 1.000 | no |
| 1100 | 19.5% | 0.078 | 0.18× | 0.004 | 0.000 | 1.000 | no |
| 1500 | 26.3% | 0.368 | 0.83× | 0.006 | 0.000 | 1.000 | no |
| **1900** | **32.9%** | **0.833** | **1.88×** | **0.517** | **0.495** | **0.693** | **YES** |
| 2300 | 39.5% | 0.833 | 1.88× | 0.517 | 0.495 | 0.693 | YES |
| 3200 | 53.7% | 0.833 | 1.88× | 0.517 | 0.495 | 0.693 | YES |

Note that every latched arm lands on the *same* attractor — 0.495 demyelination,
0.693 cognitive score — regardless of how much CO was inhaled. Severity above the
threshold is set by the switch, not by the dose.

The two things clinicians find most distinctive about DNS are consequences of
this structure rather than assumptions written into it:

- **the lucid interval** is the time the switch takes to climb (demyelination
  0.195 at 21 days, 0.466 at 42 days, plateau by 90);
- **the all-or-none character** is bistability.

The treatment implication is sharp and appears in the scenario table: oxygen by
any route — cannula, mask, 24-hour mask, intubation, carbogen — leaves this
patient above the threshold and does **not** prevent DNS, while a single early
hyperbaric session pulls the adduct peak to 0.300 and prevents it entirely.
Delaying the same three sessions to 20 h loses the effect completely. Oxygen
treats the fast pool; only the early chamber gets inside the window that matters.

The threshold depends only on the *ratio* `T_death/T_prol`, so the clone's
absolute speed could be calibrated to the clinical 2–6 week latency without
moving the threshold at all — which is how the latency was fixed.

### 6. Within a patient it is a step; the cohort incidence curve is a distribution of thresholds

This is the model's sharpest testable claim. Because the switch is bistable,
**no individual has a graded dose-response** — a patient either latches or does
not. The graded incidence seen in cohorts therefore cannot be a within-patient
dose-response; it must be between-patient variation in *where the step sits*. A
virtual population varying haemoglobin, ventilation, adduct-formation gain, clone
proliferation, terminal-oxidase sensitivity, cerebral flow reserve and xanthine
oxidase conversion reproduces incidences in the reported 10–40% band, and the
scatter of individual outcomes against the computed threshold is bimodal, not
continuous (see RESULT C, and tab 11 of the Shiny app).

### 7. Why the hyperbaric trials disagree: the window is the adduct, not the carboxyhaemoglobin

The benefit of hyperbaric oxygen in the model decays on the timescale of adduct
formation — near-maximal within the first hour or two, roughly halved by 6–8 h,
gone by 24 h. **Carboxyhaemoglobin is already below 5% in every arm beyond about
6 h**, so it cannot be what the window is tracking.

That is a structural account of a famous discordance. Weaver 2002 delivered the
first session at a median of roughly 4 h and gave three within 24 h: positive.
Scheinkestel 1999 treated after intensive-care stabilisation, frequently beyond
12 h: null. Same therapy, opposite sides of a window neither trial knew it was
straddling. The model does not settle whether hyperbaric oxygen works — it
offers a reason why two competent trials would disagree, which is a hypothesis
about the trials rather than evidence about the therapy.

The model also separates hyperbaric oxygen's three actions, which is worth doing
because they have different windows:

1. faster CO clearance — small, early, and capped by the floor in result 1;
2. dissolved oxygen covering demand outright — instant, and the mechanism in result 2;
3. **direct blockade of β2-integrin neutrophil adhesion** — neither CO clearance
   nor oxygen delivery, persisting for hours after the patient leaves the
   chamber, and in the model the action that carries most of the benefit.

If (3) dominates, a low-pressure chamber should retain most of the effect —
which the model predicts and which is testable.

### 8. The monitor becomes more reassuring as the patient becomes more poisoned

Carboxyhaemoglobin absorbs at 660 nm almost exactly as oxyhaemoglobin does and
is nearly invisible at 940 nm. Inverting the ratio-of-ratios on the CO-free
assumption therefore attributes the CO-occupied fraction to *oxygenated*
haemoglobin — SpO₂ reports (O₂Hb + COHb):

| true COHb | true SaO₂ | **SpO₂ displayed** | saturation gap |
|---|---|---|---|
| 20% | 77.6% | **94.6%** | 17.0 |
| 30% | 67.9% | **92.9%** | 25.0 |
| 40% | 58.2% | **90.8%** | 32.6 |
| 50% | 48.5% | **87.9%** | 39.4 |

The PaO₂ on a blood gas is normal too, because dissolved oxygen is untouched. So
the two instruments most likely to be reached for both read reassuringly, and
the saturation gap *is* the COHb — which makes multi-wavelength co-oximetry not
a refinement of pulse oximetry but the only instrument in the room that can see
the poison.

### 9. In pregnancy the treatment duration is set by the fetal compartment, not the maternal number

Fetal haemoglobin binds CO about 1.8× more avidly and the fetal compartment
equilibrates slowly, so the fetus lags, peaks higher, and clears far more slowly:

| therapy | maternal COHb < 5% | **fetal** COHb < 5% | ratio |
|---|---|---|---|
| room air | 904 min | 1592 min | 1.76× |
| non-rebreather mask | 204 min | **755 min** | **3.69×** |
| intubated, FiO₂ 1.0 | 157 min | 710 min | 4.51× |

The bedside teaching to continue oxygen for roughly five times the interval
needed to normalise the mother comes out of the fetal time constant rather than
being an independent rule. Note the direction of the effect: **the more
effectively you treat the mother, the larger the discrepancy becomes**, because
maternal clearance accelerates and fetal clearance barely does. This is also the
one place where prolonging normobaric oxygen past COHb clearance is not futile.

### 10. Fire smoke: two toxins on one axis, two antidotes with nothing in common

Cyanide inhibits the same terminal oxidase, so the occupancies compose on a
single axis, `f_total = 1 − (1 − f_CO)(1 − f_CN)`, and the two poisonings are
close to indistinguishable at the bedside. Their treatments share nothing:
oxygen cannot displace cyanide and hydroxocobalamin cannot displace carbon
monoxide. A lactate that will not fall on 100% oxygen is the discriminator, and
in the model it is the only one.

### 11. Hyperventilation buys pharmacokinetics and pays in perfusion

Raising ventilation shrinks `P_L/V̇_A` and speeds CO elimination, but it drops
PaCO₂ and constricts the cerebral circulation. In the model the two effects run
in opposite directions over the whole plausible range, so there is no ventilator
setting at which hyperventilating a CO-poisoned brain improves its oxygen
delivery. The escape is to decouple them: carbogen raises V̇_A while clamping
PaCO₂ with inspired CO₂, which converts the trade-off into a gain (scenario 06
against scenario 05).

---

## Scenario table

14 arms, from `co_reference_output.txt`. The severe exposure is 2298 ppm × 60 min
(calibrated peak COHb 40%); the critical exposure is 4198 ppm × 45 min (55%).
Peak COHb in this table is read off the 90-day output grid and so slightly
under-reads the true peak; the calibrated values are the authoritative ones.

| scenario | CcO peak | ATP brain | ATP watershed | lactate | MBPad | demyel. | troponin | CK | cognitive | DNS |
|---|---|---|---|---|---|---|---|---|---|---|
| 01 severe, no treatment | 35.0% | 0.616 | 0.597 | 11.8 | 0.946 | 0.495 | 17.3 | 7450 | 0.693 | **YES** |
| 02 severe, nasal cannula | 19.2% | 0.816 | 0.790 | 5.8 | 0.832 | 0.495 | 0.0 | 100 | 0.693 | **YES** |
| 03 severe, O₂ mask 6 h | 22.0% | 0.768 | 0.719 | 6.7 | 0.833 | 0.495 | 0.0 | 100 | 0.693 | **YES** |
| 04 severe, O₂ mask 24 h | 22.0% | 0.768 | 0.719 | 6.7 | 0.833 | 0.495 | 0.0 | 100 | 0.693 | **YES** |
| 05 severe, intubated + hyperventilation | 22.6% | 0.707 | 0.666 | 8.2 | 0.833 | 0.495 | 1.2 | 100 | 0.693 | **YES** |
| 06 severe, carbogen 95/5 | 19.9% | 0.763 | 0.713 | 6.0 | 0.833 | 0.495 | 0.0 | 100 | 0.693 | **YES** |
| 07 severe, HBO ×1 at 2 h, 3.0 ATA | 17.2% | 0.768 | 0.719 | 5.6 | **0.300** | **0.000** | 0.0 | 100 | **1.000** | – |
| 08 severe, HBO ×3 at 2 h (Weaver) | 17.2% | 0.768 | 0.719 | 5.6 | **0.292** | **0.000** | 0.0 | 100 | **1.000** | – |
| 09 severe, HBO ×3 delayed to 20 h | 22.0% | 0.768 | 0.719 | 6.7 | 0.833 | 0.495 | 0.0 | 100 | 0.693 | **YES** |
| 10 severe, HBO ×3 at 2.0 ATA | 17.8% | 0.768 | 0.719 | 5.9 | **0.378** | **0.000** | 0.0 | 100 | **1.000** | – |
| 11 critical, O₂ mask | 26.7% | 0.731 | 0.706 | 8.4 | 0.833 | 0.495 | 2.7 | 760 | 0.693 | **YES** |
| 12 critical, HBO ×3 early | 20.8% | 0.738 | 0.706 | 6.7 | **0.423** | **0.000** | 0.4 | 264 | **1.000** | – |
| 13 severe **+ anaemia Hb 9** | **47.6%** | **0.352** | **0.290** | **14.9** | 0.981 | **0.623** | 13.2 | 916 | **0.592** | **YES** |
| 14 severe + N-acetylcysteine | 22.0% | 0.762 | 0.714 | 6.7 | 0.833 | 0.495 | 0.0 | 100 | 0.693 | **YES** |

Four things in this table are worth reading carefully.

- **Arms 03 and 04 are identical.** Extending normobaric oxygen from 6 to 24
  hours changes nothing, because COHb is long gone either way. The one exception
  is pregnancy (result 9), where the fetal compartment is still loaded.
- **Arm 10 keeps almost the whole effect at 2.0 ATA**, which is what result 7
  predicts if the anti-adhesion action rather than the pressure is doing the work.
- **Arm 13 is the worst in the table** despite an unremarkable exposure. Anaemia
  is the comorbidity the model punishes hardest, and it is the only arm that
  produces frank necrosis.
- **Arm 01 is the only arm with myocardial injury and rhabdomyolysis**, and arm 02
  abolishes both — oxygen is highly effective at the things COHb drives directly,
  and ineffective at the thing that determines the neurological outcome. That
  dissociation is the whole model in one row-pair.

---

## Model structure

45 ODEs in 16 blocks. Time is in minutes throughout.

| block | states |
|---|---|
| environment | room CO |
| pulmonary exchange | CFK mass balance on COHb |
| tissue CO | brain, myocardium, skeletal muscle, splanchnic |
| myoglobin | cardiac, skeletal |
| terminal oxidase | brain, myocardium |
| energetics | brain, **watershed (globus pallidus)**, myocardium, lactate |
| haemodynamics | cerebral blood flow, PaCO₂ |
| reoxygenation injury | xanthine oxidase, ROS, NO/peroxynitrite, glutathione |
| innate | adherent neutrophils, myeloperoxidase, lipid peroxidation |
| adaptive | MBP adduct, autoreactive clone, microglia, demyelination |
| protective | haem oxygenase-1 (a positive feedback — it makes more CO), HIF-1α |
| organ injury | necrosis, oedema, ICP, troponin, ejection fraction, CK, creatinine |
| outcome | cognitive composite, cumulative DNS hazard |
| cyanide | blood, tissue, hydroxocobalamin |
| adjuncts | N-acetylcysteine (2-cpt), oxypurinol |
| fetus | fetal COHb |

Two modelling choices are worth flagging because they are not the obvious ones:

- **The tissue O₂ constant is 12 mmHg, not the mitochondrial ~1 mmHg.** Cells at
  the far edge of the Krogh diffusion field fail long before the enzyme runs out
  of substrate, so the supply function is a tissue-level, not an enzyme-level,
  quantity. With the enzyme Km the model produces almost no hypoxic injury at
  any COHb, which is wrong.
- **The globus pallidus watershed is not handicapped at baseline; it is
  steeper.** It sits lower on its own oxygen supply curve, so the same
  *fractional* fall in tissue PO₂ costs it more. Giving it a fixed baseline
  handicap instead would have made every healthy subject permanently ischaemic
  there.

### Calibration anchors

- CFK constants: resting physiology, not fitted. Reproduce 320 → **313** min on
  air and 74 → **72** min on a mask.
- Baseline COHb **0.618%** (endogenous haem catabolism 0.007 mL/min plus 1.5 ppm
  ambient) against a measured normal of 0.4–0.8%.
- Exposure scale calibrated by root-finding on peak COHb: 538 ppm → 10%,
  1409 → 25%, 2298 → 40%, 4198 ppm × 45 min → 55%. For scale, the OSHA 8-hour
  limit is 50 ppm, NIOSH IDLH is 1200 ppm, and a petrol engine in a closed
  garage exceeds 30 000 ppm.
- DNS incidence anchored on Weaver 2002 (46% with normobaric oxygen, 25% with
  hyperbaric).
- **The healthy unexposed subject is an exact fixed point of all 45 states.**

---

## What the independent implementation caught

Every equation was written twice. The Python implementation
(`co_python_reference.py`) is not a port — it was written first, executed, and
then transcribed to mrgsolve. It exposed six defects that would each have
produced confident, wrong results:

1. **The healthy subject self-destructed.** With injury drives written as
   absolute rather than excess-over-baseline quantities, an unexposed
   subject fully demyelinated in ten hours (lipid peroxidation reached 1412,
   demyelination 0.9999, cognitive score 0.38). Every driver now derives from a
   quantity that is *exactly* zero in health, and the fixed point is verified
   numerically on every run.
2. **The untreated arm was being treated.** The oxygen-therapy start defaulted
   to t = 0 rather than infinity, so the control arm in every comparison was
   silently receiving 85% oxygen. This is the defect that would have been
   hardest to see from the outputs, because every arm still looked plausible.
3. **Hyperbaric oxygen caused rhabdomyolysis.** Tissue CO loading was driven by
   the Haldane free tension `COHb·P_c,O₂/(M·[O₂Hb])`, which is a *lung-exchange*
   construct and rises about 20-fold inside a chamber at fixed COHb. The model
   therefore pumped CO into muscle during treatment and produced CK 4803 with
   HBO against 100 without — the wrong sign on the flagship therapy. Tissue
   loading is now driven by blood CO content; the Haldane term stays in the lung.
4. **Delayed sequelae were structurally unreachable.** With adduct clearance on a
   2.5-day time constant the primed clone could never outrun antigen loss, so DNS
   occurred at *no dose whatsoever* while every upstream variable looked correct.
   Peroxidation-modified myelin clears on a myelin-turnover timescale (weeks),
   not in days.
5. **The watershed compartment sat at 82% of baseline in health**, eating the
   reserve it was supposed to lose during poisoning. It needed normalising to its
   own reference, which also changed its meaning for the better (see above).
6. **Myocardial injury never occurred** because the energy threshold sat below
   anything reachable, and **N-acetylcysteine's peripheral compartment had the
   wrong volume scaling** in the two-compartment transfer term.

The cross-check also refuted a prior expectation of the author's: hyperbaric
oxygen's benefit was assumed to come mainly from accelerated CO clearance. The
floor derivation in result 1 and the mechanism split in result 7 both say it
does not, and that the anti-adhesion action — the one that is neither clearance
nor oxygenation — carries most of it.

---

## Known limits

- **CFK is falsified at hyperbaric pressure** (result 1). The model is simply
  wrong in the chamber, by about a factor of two on the half-life, and no attempt
  was made to hide this by refitting per pressure.
- **The DNS mechanism is a strong commitment to an animal model.** Immune-mediated
  demyelination after CO rests principally on rodent work from one group. If it
  is wrong, results 5–7 go with it.
- **The tolerance threshold's functional form is invented.** A Hill function with
  n = 6 gives bistability; the true shape is unknown. The sensitivity analysis
  shows the 90-day cognitive score depends more on this arm than on the CO
  pharmacokinetics, which is uncomfortable and honest.
- **N-acetylcysteine, allopurinol and carbogen are hypotheses**, included because
  they are the mechanistically obvious targets, not because human outcome
  evidence supports them.
- **Cyanide is lumped** into one tissue compartment with a single rhodanese term
  and has not been calibrated against fire-victim data.
- **The watershed is one compartment.** Real CO lesions are patchy and the model
  cannot speak to their distribution.
- The DNS "incidence" figures depend on an arbitrary cut (demyelination > 0.05)
  used to declare a case; the bimodality is a model result, the cut is not.

---

## Usage

```r
library(mrgsolve); library(dplyr)
mod <- mread("co_mrgsolve_model.R")

# severe accidental exposure, non-rebreather mask from 90 min
mod %>%
  param(ppm_fix = 2298, texp = 60, to2_start = 90, to2_stop = 450) %>%
  mrgsim(end = 90*1440, delta = 5) %>%
  plot(COHb_pct + CcOb + MBPad + Cog ~ time)

# the Weaver / Scheinkestel contrast: same therapy, 2 h against 20 h
early <- mod %>% param(ppm_fix = 2298, texp = 60, to2_start = 90,
                       to2_stop = 2970, thbo_start = 120, ATA = 3, nhbo = 3)
late  <- mod %>% param(ppm_fix = 2298, texp = 60, to2_start = 90,
                       to2_stop = 2970, thbo_start = 1200, ATA = 3, nhbo = 3)
```

```bash
# regenerate the map
dot -Tsvg co_qsp_model.dot -o co_qsp_model.svg
dot -Tpng -Gdpi=150 co_qsp_model.dot -o co_qsp_model.png

# reproduce every number in this README
python3 co_python_reference.py     # writes co_reference_output.txt
```

```r
shiny::runApp("co_shiny_app.R")    # 12 tabs
```

---

> **면책 조항.** 본 모델은 교육 및 연구 목적의 QSP 모델입니다. 공개 문헌을 바탕으로
> 구성되었으나 독립적으로 검증·인증되지 않았으며, 실제 임상 의사결정, 처방, 또는
> 규제 제출에 직접 사용해서는 안 됩니다. 특히 고압산소 치료의 적응증과 시점에 관한
> 본 모델의 예측은 임상 근거가 아니라 가설입니다.
