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
cross-check found and fixed nine defects, listed at the end.

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
| COHb — *the measured one* | 153 min | 366 min |
| brain tissue CO | 168 min | 383 min |
| cardiac myoglobin-CO | 237 min | 533 min |
| skeletal muscle CO | 620 min | 6075 min |
| **brain cytochrome c oxidase** | **750 min** | **7580 min** |

COHb reaches the "normal" 5% at **333 min (5.5 h)**. At that moment brain
cytochrome c oxidase is still **20.3% CO-inhibited**, and it does not reach 10%
of its own peak until **7580 min (5.3 days)** — **20.7× later** than COHb does.
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
puts it between COHb 30% and 33% for a 1-hour exposure treated promptly with a
mask. Note the 1700 ppm row: the adduct peak is 1.11× the threshold and there is
still no DNS — so it is not the *peak* that decides but the **time spent above
the threshold**, i.e. whether the primed clone can outrun antigen clearance. The
separatrix lives in the (peak × duration) plane, and `MBPad_crit` is a necessary
but not sufficient condition:

| ambient CO | peak COHb | peak MBPad | ÷ threshold | clone | demyelination | cognitive | DNS |
|---|---|---|---|---|---|---|---|
| 700 ppm | 12.8% | 0.024 | 0.05× | 0.004 | 0.000 | 1.000 | no |
| 1100 | 19.7% | 0.078 | 0.18× | 0.004 | 0.000 | 1.000 | no |
| 1500 | 26.6% | 0.368 | 0.83× | 0.006 | 0.000 | 1.000 | no |
| 1700 | 30.0% | 0.492 | 1.11× | 0.025 | 0.000 | 1.000 | no |
| **1900** | **33.3%** | **0.833** | **1.88×** | **0.517** | **0.495** | **0.693** | **YES** |
| 2300 | 40.0% | 0.833 | 1.88× | 0.517 | 0.495 | 0.693 | YES |
| 3200 | 54.7% | 0.833 | 1.88× | 0.517 | 0.495 | 0.693 | YES |

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
dose-response; it must be between-patient variation in *where the step sits*.

A virtual population of 150 (varying haemoglobin, alveolar ventilation,
adduct-formation gain, clone proliferation, terminal-oxidase sensitivity,
cerebral flow reserve and xanthine oxidase conversion), all treated with a mask:

| exposure | DNS incidence | mean cognitive score at 90 d |
|---|---|---|
| mild (COHb 10%) | 8.7% | 0.967 |
| moderate (COHb 25%) | **28.7%** | 0.903 |
| severe (COHb 40%) | 75.3% | 0.754 |
| critical (COHb 55%) | 96.0% | 0.669 |

The same populations given early hyperbaric oxygen:

| exposure | DNS, mask only | DNS, early HBO ×3 | mean cognitive score |
|---|---|---|---|
| moderate (COHb 25%) | 28.7% | **5.3%** | 0.980 |
| severe (COHb 40%) | 75.3% | **23.3%** | 0.921 |

**This is a partial calibration success, and the pattern of the miss is
informative.** Against Weaver 2002 — 46% cognitive sequelae at six weeks with
normobaric oxygen, 25% with hyperbaric — the model's **treated** arm is nearly
exact (23.3% against 25%) while its **untreated** arm badly over-predicts (75.3%
against 46%). So the model is not simply over-sensitive everywhere: it gets the
floor right and the ceiling wrong. **The switch is too easy to cross at high dose
in the absence of treatment**, which inflates the apparent treatment effect (69%
relative reduction against an observed 46%) even though the post-treatment rate
is right.

The mild arm (8.7%) and moderate arm (28.7%) both sit inside the reported 10–40%
band; the critical arm's 96% is not credible. No attempt has been made to tune
this away by widening the population variability until the numbers matched.

### 7. Why the hyperbaric trials disagree: the window is the adduct, not the carboxyhaemoglobin

The benefit of hyperbaric oxygen in the model decays on the timescale of adduct
formation, and because the outcome is bistable the decay presents as a **cliff**:

| first session at | COHb then | peak MBPad | demyelination | cognitive |
|---|---|---|---|---|
| none | – | 0.833 | 0.495 | 0.693 |
| 65 min | 48.3% | 0.282 | 0.000 | **1.000** |
| 120 min | 26.4% | 0.493 | 0.000 | **1.000** |
| 180 min | 15.6% | 0.563 | 0.002 | **0.999** |
| **210 min** | **12.2%** | **0.833** | **0.495** | **0.693** |
| 360 min | 4.1% | 0.833 | 0.495 | 0.693 |
| 1440 min | 0.5% | 0.833 | 0.495 | 0.693 |

Full benefit up to about 3 hours, none at all from 3.5 hours onward. Two things
about that table matter more than the cliff's exact position. First, **the COHb
column runs the wrong way**: it is *highest* (48.3%) in the arm that works best
and already below 5% in every arm that fails, so carboxyhaemoglobin cannot be
what the window is tracking. Second, the model's window is *narrower* than the
6-hour figure in practice guidelines — bistable systems give sharp edges, and the
real edge is presumably blurred by the between-patient variation in result 6.

That is a structural account of a famous discordance. Weaver 2002 delivered the
first session at a median of roughly 4 h and gave three within 24 h: positive.
Scheinkestel 1999 treated after intensive-care stabilisation, frequently beyond
12 h: null. Same therapy, opposite sides of a window neither trial knew it was
straddling. The model does not settle whether hyperbaric oxygen works — it
offers a reason why two competent trials would disagree, which is a hypothesis
about the trials rather than evidence about the therapy.

The model also separates hyperbaric oxygen's three actions:

1. faster CO clearance — small, early, and capped by the floor in result 1;
2. dissolved oxygen covering demand outright — instant, and the mechanism in result 2;
3. direct blockade of β2-integrin neutrophil adhesion — neither CO clearance nor
   oxygen delivery, persisting for hours after the patient leaves the chamber.

Knocking these out one at a time gives a result that **contradicted what the
author expected**, and the expectation was wrong rather than the model:
**the actions are redundant, not additive.** At 3.0 ATA, deleting the
anti-adhesion action entirely still prevents delayed sequelae (adduct peak rises
from 0.423 to 0.513, outcome unchanged). At 1.35 ATA — barely any extra dissolved
oxygen — keeping the anti-adhesion action still prevents them (0.542). Only
removing *both* fails (0.833, full latch). Either mechanism alone is sufficient,
which is why 2.0 ATA works as well as 3.0 and why the effect is robust to which
of the two you believe in.

The same sweep says something awkward about protocol: **session number makes no
difference at all** — one session, three and six give identical outcomes and
adduct peaks agreeing to three decimals. The Weaver protocol's three sessions
have no rationale inside this model. Either the model is missing whatever the
repeat sessions do, or they do nothing; the model cannot tell which, but it does
make the question sharp.

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
single axis, `f_total = 1 − (1 − f_CO)(1 − f_CN)`:

| arm | CcO from CO | CcO total | ATP brain | lactate | MBPad | cognitive |
|---|---|---|---|---|---|---|
| CO alone (COHb 30%) | 16.3% | 16.3% | 0.835 | 4.6 | 0.467 | **1.000** |
| CN alone | 3.1% | 66.8% | 0.351 | 18.1 | 0.833 | 0.684 |
| **CO + CN (fire smoke)** | 16.3% | **71.3%** | 0.318 | 20.0 | 0.833 | 0.669 |
| + hydroxocobalamin at 5 min | 16.3% | **24.0%** | 0.819 | 4.8 | 0.495 | **1.000** |
| + hydroxocobalamin at 25 min | 16.3% | 65.6% | 0.439 | 13.8 | 0.833 | 0.693 |
| + hydroxocobalamin at 50 min | 16.3% | 71.3% | 0.320 | 18.5 | 0.833 | 0.677 |
| **+ early HBO ×3, no cobalamin** | 12.8% | **70.1%** | 0.318 | 20.0 | 0.833 | **0.669** |

The last row is the one that matters. **Early hyperbaric oxygen — the
intervention that prevents delayed sequelae in every pure-CO arm in this
model — does nothing at all here** (0.669, identical to no treatment), because it
addresses the 16% of the enzyme that CO is holding and cannot touch the 67% that
cyanide is. Conversely cobalamin, which rescues completely at 5 minutes, cannot
displace carbon monoxide. The two antidotes are not interchangeable and not
additive; each is useless against the other's ligand, while the bedside picture
is nearly the same. A lactate of 20 that will not fall on 100% oxygen is the
discriminator, and in the model it is the only one.

Note also how fast the cobalamin window closes: complete rescue at 5 minutes,
none by 50. Like everything else in this disease, it is decided early.

### 11. Hyperventilation buys pharmacokinetics and pays in perfusion

Raising ventilation shrinks `P_L/V̇_A` and speeds CO elimination, but it drops
PaCO₂ and constricts the cerebral circulation. In the model the two effects run
in opposite directions over the whole plausible range, so there is no ventilator
setting at which hyperventilating a CO-poisoned brain improves its oxygen
delivery. The escape is to decouple them: carbogen raises V̇_A while clamping
PaCO₂ with inspired CO₂, which converts the trade-off into a gain (scenario 06
against scenario 05).

### 12. The long-term outcome is almost completely insensitive to the CO pharmacokinetics

Perturbing each parameter by ±20/25% and reading the 90-day cognitive score:

| parameter | \|ΔCog\|/Cog | −20% | +25% | what it is |
|---|---|---|---|---|
| `ATPthr` | **0.574** | 1.000 | 0.602 | energy reserve before the cascade ignites |
| `theta` | **0.446** | 0.690 | 0.999 | tolerance threshold for clonal expansion |
| `Hb` | **0.442** | 0.693 | 0.999 | haemoglobin |
| `kRepair` | 0.101 | 0.658 | 0.728 | remyelination |
| `Tprol` / `Tdeath` | 0.099 | — | — | clone kinetics |
| `kDemy` | 0.077 | 0.719 | 0.666 | demyelination rate |
| `kAdClr`, `kSpread` | 0.002 | — | — | adduct clearance, epitope spreading |
| `kon_cco`, `Kc`, `Ko`, `kappa`, `MHald`, `koff_cco` | **0.0002** | 0.693 | 0.693 | **the entire CO pharmacokinetic block** |
| `f_wshed`, `hbo_adh` | **0** | — | — | no influence at all |

Every parameter governing how carbon monoxide is taken up, distributed, bound and
eliminated — including the Haldane ratio itself — moves the 90-day cognitive
score by 0.02%. The parameters that matter are the energy-failure threshold, the
immunological tolerance threshold, and haemoglobin.

This needs one honest caveat: because the outcome is bistable, these derivatives
mostly measure *whether a parameter flips the switch*, and the reference arm is
already latched. So the correct reading is not "CO kinetics are unimportant in
CO poisoning" — they set the troponin, the lactate and the rhabdomyolysis, as the
scenario table shows. It is that **once the switch has been crossed, the
pharmacokinetics no longer have any purchase on the neurological outcome**, which
is the same conclusion as results 7 and 11 arrived at from different directions.
It also means the model's most consequential parameters are the two least
constrained by data — stated plainly in the limits below.

---

## Scenario table

14 arms, from `co_reference_output.txt`. The severe exposure is 2298 ppm × 60 min
and the critical one 4198 ppm × 45 min; both reproduce their calibration targets
(peak COHb 40.0% and 55.0%) exactly.

| scenario | COHb | CcO peak | ATP brain | ATP watershed | lactate | MBPad | demyel. | troponin | CK | cognitive | DNS |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 01 severe, no treatment | 40.0% | 35.0% | 0.616 | 0.597 | 11.8 | 0.946 | 0.495 | 17.4 | 7458 | 0.693 | **YES** |
| 02 severe, nasal cannula | 40.0% | 27.9% | 0.708 | 0.695 | 9.1 | 0.896 | 0.495 | 5.0 | 2353 | 0.693 | **YES** |
| 03 severe, O₂ mask 6 h | 40.0% | 22.0% | 0.762 | 0.714 | 6.8 | 0.833 | 0.495 | 0.0 | 100 | 0.693 | **YES** |
| 04 severe, O₂ mask 24 h | 40.0% | 22.0% | 0.762 | 0.714 | 6.8 | 0.833 | 0.495 | 0.0 | 100 | 0.693 | **YES** |
| 05 severe, intubated + hyperventilation | 40.0% | 22.7% | 0.705 | 0.656 | 8.3 | 0.833 | 0.495 | 1.3 | 100 | 0.693 | **YES** |
| 06 severe, carbogen 95/5 | 40.0% | 19.9% | 0.757 | 0.708 | 6.1 | 0.833 | 0.495 | 0.0 | 100 | 0.693 | **YES** |
| 07 severe, HBO ×1 at 2 h, 3.0 ATA | 40.0% | 18.3% | 0.762 | 0.714 | 5.9 | **0.300** | **0.000** | 0.0 | 100 | **1.000** | – |
| 08 severe, HBO ×3 at 2 h (Weaver) | 40.0% | 18.3% | 0.762 | 0.714 | 5.9 | **0.292** | **0.000** | 0.0 | 100 | **1.000** | – |
| 09 severe, HBO ×3 delayed to 20 h | 40.0% | 22.0% | 0.762 | 0.714 | 6.8 | 0.833 | 0.495 | 0.0 | 100 | 0.693 | **YES** |
| 10 severe, HBO ×3 at 2.0 ATA | 40.0% | 18.3% | 0.762 | 0.714 | 5.9 | **0.378** | **0.000** | 0.0 | 100 | **1.000** | – |
| 11 critical, O₂ mask | 55.0% | 26.8% | 0.691 | 0.621 | 8.4 | 0.833 | 0.495 | 2.8 | 761 | 0.693 | **YES** |
| 12 critical, HBO ×3 early | 55.0% | 21.4% | 0.691 | 0.621 | 6.9 | **0.423** | **0.000** | 0.5 | 267 | **1.000** | – |
| 13 severe **+ anaemia Hb 9** | **60.8%** | **47.7%** | **0.350** | **0.286** | **14.9** | 0.981 | **0.623** | 13.3 | 925 | **0.592** | **YES** |
| 14 severe + N-acetylcysteine | 40.0% | 22.0% | 0.762 | 0.714 | 6.8 | 0.833 | 0.495 | 0.0 | 100 | 0.693 | **YES** |

Five things in this table are worth reading carefully.

- **Oxygen is graded on the acute markers and useless on the neurological one.**
  Going from nothing (01) to cannula (02) to mask (03) walks CcO inhibition down
  35.0 → 27.9 → 22.0%, abolishes the troponin rise (17.4 → 5.0 → 0.0) and the
  rhabdomyolysis (7458 → 2353 → 100), and moves the neurological outcome not at
  all. Every one of those arms latches to the same 0.693. **That dissociation is
  the whole model in three rows.**
- **Arms 03 and 04 are identical.** Extending normobaric oxygen from 6 to 24
  hours changes nothing, because COHb is long gone either way. The one exception
  is pregnancy (result 9), where the fetal compartment is still loaded.
- **Arm 10 keeps the whole effect at 2.0 ATA**, consistent with the redundancy
  finding in result 7.
- **Arm 13 is the worst in the table.** The same exposure at Hb 9 reaches COHb
  60.8% rather than 40.0% — because the same *amount* of CO occupies a larger
  *fraction* of a smaller haemoglobin pool — and it is the only arm that produces
  frank necrosis. Anaemia is the comorbidity the model punishes hardest.
- **Arm 13 is also the only arm where anything moves the demyelination number
  upward** (0.623 against the universal 0.495). Above the switch, dose stops
  mattering; only added necrosis makes it worse.

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

Three more were found late, in the *analysis pipeline* rather than the model, and
they are the most instructive of the set because each produced a plausible number
rather than an obvious failure:

7. **The fire-smoke scenario had been running with no cyanide in it.** The driver
   passed `CN_rate` while the Python model read `CN_dose`; the dictionary update
   silently accepted the unknown key. Every CO+CN arm was numerically identical to
   CO alone and nobody would have noticed from the output. The fix was to make the
   parameter constructor **reject unknown keys**, and that check immediately found
   the next one.
8. **The "nasal cannula" arm was simulating a different exposure, not a different
   therapy.** The Python model had no `FiO2_trt` (the treatment oxygen fraction
   was hard-coded at 0.85) while the R file exposed it as a parameter, so passing
   `FiO2 = 0.44` changed what the patient breathed *during the poisoning*. The arm
   showed a lower peak COHb, which looked like a plausible result for a weaker
   oxygen mask and was in fact a different experiment. `shunt` had drifted the
   same way — a function default in one implementation, a parameter in the other.
9. **A uniform output grid was missing every acute peak.** Over a 90-day horizon a
   uniform grid samples every ~108 minutes, so all "peak" statistics were read off
   points that straddled the real maximum. This did not merely add noise — it
   inverted a conclusion, making a well-timed dose of hydroxocobalamin appear to
   *raise* peak tissue cyanide (11.22 with the drug against 5.57 without), because
   adding a dose event added grid points that caught more of the comparator's
   missed peak. The grid is now dense over the first 12 hours. After the fix, the
   exposure calibration reproduces its targets exactly (10.0 / 25.0 / 40.0 / 55.0%)
   where it had been reading 37.1% for a 40% target.

The cross-check also **refuted a prior expectation of the author's, twice.**
Hyperbaric oxygen's benefit was assumed to come mainly from accelerated CO
clearance; the floor derivation in result 1 says it cannot. It was then assumed to
come mainly from the anti-adhesion action instead; the mechanism split in result 7
says that is also wrong, because the two actions are *redundant* — either alone
suffices and only removing both loses the effect. Both hypotheses were stated
before the simulations were run, and both were wrong.

---

## Known limits

- **CFK is falsified at hyperbaric pressure** (result 1). The model is simply
  wrong in the chamber, by about a factor of two on the half-life, and no attempt
  was made to hide this by refitting per pressure.
- **The DNS mechanism is a strong commitment to an animal model.** Immune-mediated
  demyelination after CO rests principally on rodent work from one group. If it
  is wrong, results 5–7 go with it.
- **The model's two most influential parameters are its two least constrained.**
  The sensitivity analysis puts `ATPthr` (the energy reserve before the oxidative
  cascade ignites) and `theta` (the immunological tolerance threshold) at the top
  of the list, ahead of haemoglobin and far ahead of every CO pharmacokinetic
  parameter. Neither has direct experimental support; the Hill exponent of 6 was
  chosen because it gives bistability, not because the shape is known. Results
  5–7 and 11 should be read as conditional on these two numbers.
- **DNS incidence is over-predicted at high dose.** 75.3% at COHb 40% against
  Weaver's observed 46%, and 96% at COHb 55%. The switch is too easy to cross.
  The moderate-exposure arm (28.7%) is well calibrated; the severe one is not.
- **The HBO benefit is over-stated, but only through the control arm.** The
  treated rate at severe exposure (23.3%) nearly matches Weaver's 25%; the
  untreated rate (75.3%) over-predicts his 46%, which inflates the relative
  effect to 69% against an observed 46%.
- **The therapeutic window is too sharp.** Bistable systems give cliff edges, and
  the model's is at ~3.5 h against the ~6 h in practice guidelines. The real edge
  is presumably blurred by between-patient variation the single-subject runs do
  not show.
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
