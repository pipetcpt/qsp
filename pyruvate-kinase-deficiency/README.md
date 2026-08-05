# Pyruvate Kinase Deficiency (PKLR) — QSP Model

**Congenital non-spherocytic haemolytic anaemia · the commonest glycolytic enzymopathy**

<a href="pkd_qsp_model.svg"><img src="pkd_qsp_model.png" width="620" alt="PKD QSP mechanistic map"></a>

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (190 nodes, 283 edges, 20 clusters) | [`pkd_qsp_model.dot`](pkd_qsp_model.dot) · [SVG](pkd_qsp_model.svg) · [PNG](pkd_qsp_model.png) |
| ⚙️ mrgsolve model (80 ODEs, 24 scenarios) | [`pkd_mrgsolve_model.R`](pkd_mrgsolve_model.R) |
| 🐍 Executable reference implementation | [`pkd_reference_model.py`](pkd_reference_model.py) |
| 📊 Verified run output | [`pkd_reference_output.txt`](pkd_reference_output.txt) · [`pkd_population_results.json`](pkd_population_results.json) |
| 📱 Shiny dashboard (10 tabs) | [`pkd_shiny_app.R`](pkd_shiny_app.R) |
| 📚 References (267 PMIDs) | [`pkd_references.md`](pkd_references.md) |

---

## The one structural choice everything follows from

Pyruvate kinase is the **last ATP-generating step of glycolysis**, and it sits
**downstream** of the 1,3-bisphosphoglycerate branch point that feeds the
Rapoport–Luebering shunt. So a single enzyme lesion moves two quantities in
**opposite physiological directions**:

```
PK activity ↓  →  ATP ↓       →  cation pumps fail, the cell dies sooner
PK activity ↓  →  2,3-BPG ↑   →  the O₂ dissociation curve shifts RIGHT, so every
                                 surviving gram of haemoglobin unloads MORE oxygen
```

This model refuses to collapse those into one "severity" scalar. They are two
trunks out of the same node, and they reach the clinical endpoints from opposite
sides. The consequence is not rhetorical, and the model computes it:

> **A pyruvate kinase activator raises haemoglobin *by lowering 2,3-BPG*.**
> It adds oxygen carrier and subtracts oxygen unloading in the same act.
> Haemoglobin response is the registrational endpoint of every trial in this
> disease — ACTIVATE, ACTIVATE-T, DRIVE-PK — and it is **not a sufficient
> statistic** for the physiological benefit.

The second choice is that **cell age is an explicit state axis** (14 cohorts).
Mutant PK-R is thermolabile and an anucleate cell cannot replace it, so the
lesion *deepens* as a cell ages. Splenectomy, a PK activator, and gene therapy
act on three different parts of that axis, which is why they have qualitatively
different kinetics rather than merely different potencies.

---

## What the model gets, and what it costs

**Five fitted parameters.** Marrow output (set so a wild-type subject sits at
Hb 15.0 g/dL), the mutant PKR decay constant, the two haemolysis gains, and the
fraction of the ATP-dependent hazard an activator can reverse acutely.
Everything else is literature or back-calculated from the normal erythrocyte
operating point. Nineteen quantities of healthy physiology are then
**predictions** — and they are checked before any patient is simulated:

| Predicted (not fitted) | Model | Physiological |
|---|---|---|
| red cell count | 5.00 ×10¹²/L | 4.5–5.5 |
| reticulocytes | 0.96 % | 0.5–1.5 |
| mean red cell lifespan | 103.6 d | 100–120 |
| P50 | 26.8 mmHg | 26–27 |
| mixed venous PO₂ | 37.9 mmHg | 38–42 |
| haemoglobin catabolised | 7.2 g/d | 6–7 |
| iron recycled | 25.0 mg/d | 20–25 |
| total bilirubin | 0.88 mg/dL | 0.3–1.0 |
| erythroid amplification | 9.6 × | 8–12 |

The iron flux and the bilirubin production are the *same* destruction flux seen
from two organs, so they cannot be tuned independently — which is what forces
the haemoglobin, bilirubin, iron and reticulocyte predictions in the diseased
state to move together or not at all.

---

## Results

### 1 · The compensation knee is a closed-form quantity

Maximum erythroid amplification is 64×, amplification at basal erythropoietin is
9.6×, so production can rise at most **6.70-fold** and the critical red cell
lifespan is

```
L* = 120 d / 6.70 = 17.9 d
```

Above `L*` the marrow replaces what is destroyed and **haemoglobin is flat in the
genotype**; below it, production is capped and Hb falls in proportion to
lifespan. The simulated spectrum shows exactly that — lifespan slides smoothly
from 34 d to 6 d while Hb barely moves until it crosses 18 d:

| phenotype | α | Hb | retic | lifespan | bilirubin | P50 |
|---|---|---|---|---|---|---|
| very mild / compensated | 0.30 | 15.6 | 2.4 % | 34 d | 2.8 | 29.7 |
| mild | 0.22 | 13.8 | 7.2 % | 19 d | 4.6 | 30.9 |
| moderate | 0.16 | 11.5 | 10.6 % | 12 d | 6.5 | 32.2 |
| severe, not transfused | 0.12 | 8.1 | 18.0 % | 8 d | 6.6 | 33.6 |
| transfusion dependent | 0.09 | 5.5 | 25.6 % | 6 d | 6.3 | 34.9 |

Two things follow that are not obvious at the bedside. Haemoglobin is a
**saturating** readout of the lesion above `L*` and a **hypersensitive** one just
below it — and patients cluster near `L*`, which is why a modest intervention can
look dramatic. And reticulocytes and bilirubin keep rising monotonically
*through* the knee, so in the compensated range they, not Hb, carry the severity
information.

### 2 · Half the anaemia is already paid for before treatment

A normal subject at Hb 15.0 (P50 26.8) extracts 0.249 of the oxygen carried per
litre of blood down to a tissue PO₂ of 38 mmHg. The severe patient at Hb 8.11
(P50 33.6) extracts **0.360 — 45 % more per gram**. His *equivalent haemoglobin*
— the Hb a normal-P50 subject would need to unload the same oxygen — is
**11.75 g/dL, not 8.11**.

So of a 6.89 g/dL haemoglobin deficit, the shunt has already repaid **3.64 g/dL,
53 % of it**, with no treatment at all. That is the quantitative form of a
well-known bedside observation: PK-deficient patients tolerate haemoglobins that
would incapacitate other anaemic patients. It is not stoicism, it is a
right-shifted curve.

### 3 · The break-even 2,3-BPG reduction, in closed form

A drug delivering ΔHb leaves oxygen transport unchanged when
`Hb₁·extract(P50₁) = Hb₀·extract(P50₀)`. Solving for P50 and inverting
`P50 = P50ref·(DPG/DPG₀)^n` gives the **largest 2,3-BPG fall compatible with a
net gain**:

| ΔHb | Hb reached | break-even P50 | break-even 2,3-BPG | max 2,3-BPG fall |
|---|---|---|---|---|
| +1.0 | 9.11 | 31.06 | 7.61 mM | **21.5 %** |
| +1.5 | 9.61 | 30.05 | 6.86 mM | **29.2 %** |
| +2.0 | 10.11 | 29.15 | 6.24 mM | **35.6 %** |
| +3.0 | 11.11 | 27.63 | 5.28 mM | **45.5 %** |

The +1.5 g/dL row *is* the ACTIVATE primary endpoint. Simulating the actual
titration, mitapivat lowers 2,3-BPG by 11 % / 19 % / 23 % at 5 / 20 / 50 mg BID —
**inside** the band — so oxygen transport does improve. But it improves *less
than haemoglobin says*, and the gap widens with dose:

| dose | ΔHb | Δ equivalent Hb | 2,3-BPG |
|---|---|---|---|
| 5 mg BID | +1.59 | +1.54 | 8.59 mM |
| 20 mg BID | +2.49 | +2.12 | 7.82 mM |
| 50 mg BID | +2.79 | +2.20 | 7.49 mM |

**The dose that maximises the trial endpoint is not the dose that maximises
oxygen transport.** Nothing in a trial scoring ΔHb can see that, and the
break-even table above is given precisely so the conclusion can be re-scored
against a different 2,3-BPG estimate.

### 4 · Why the haemoglobin response is fast — and what that predicts

DRIVE-PK reported a **median 10 days** to the first >1.0 g/dL rise, with a range
out to 187 days. In an age-structured model this is arithmetic, not
pharmacokinetics: a step change in hazard relaxes the red cell pool with a time
constant equal to the **new mean lifespan**, and this patient's baseline lifespan
is 8 days, not 120.

> **Prediction:** time-to-response should be *inversely* related to baseline red
> cell lifespan. The sicker the patient, the **faster** the haemoglobin moves,
> and the 187-day outliers should be the *mildest* patients. Testable in the
> existing DRIVE-PK / ACTIVATE datasets without a single new sample.

The two drug limbs separate cleanly in time. Allosteric activation is at steady
state within a day and **cannot** rescue a cohort whose protein has already
decayed; thermostabilisation lengthens τ_PK from 50 to ~104 d over months and is
the only limb that reaches old cells. So the trajectory is biphasic, with a
second phase invisible in a 24-week trial — which is where the ACTIVATE extension
found responses sustained and deepening to week 96.

The **null-allele genotype gets nothing**, with no parameter added to make that
happen: an allosteric activator multiplies residual activity, and zero times
anything is zero. DRIVE-PK found responses *only* in patients with ≥1 missense
variant, correlated with baseline PK-R protein level. That is the same statement.

### 5 · A negative result the registry data can refute

PK deficiency is widely described as showing a *paradoxical rise* in reticulocyte
count after splenectomy, attributed to splenic reticulocyte sequestration. The
model says that cannot be the mechanism, and the argument is two lines. With
reversible pooling the circulating reticulocyte pool is

```
R = influx / ( k_exit + k_in·k_kill/(k_out+k_kill) )
```

Reversible *holding* cancels out entirely — cells that are released come back.
Only the **killing** term suppresses R. Writing `q` for the fraction of marrow
output killed in the spleen, that term equals `k_exit·q/(1−q)`, so splenectomy
raises reticulocytes by `q/(1−q)` **and raises haemoglobin by the same `q/(1−q)`**,
because the same `q` is the production that was being wasted.

> **The two fractional changes must be equal.** A large reticulocyte rise with a
> small haemoglobin rise is not available from this mechanism at any parameter
> value. Sixteen combinations of the two hazard gains were searched before this
> was recognised as structural rather than a calibration failure.

The model's median splenectomy gain is **+2.09 g/dL** (observed median +1.6) with a
median reticulocyte change of **+2.3 %** — small, and of the same order as the
haemoglobin change, exactly as the identity requires. So either the reported
reticulocyte rise is of that modest size, or its cause lies elsewhere — the
remaining candidates being maturation-factor inflation of the count and persisting
marrow drive, which make *opposite* predictions for the absolute reticulocyte
count.

### 6 · Gene therapy's time constant is set by the age axis

With 45 % of marrow output corrected from day 0, haemoglobin takes **months**.
The contrast with mitapivat's 10 days is not potency: an activator lowers the
hazard on cells that already exist, whereas gene therapy can only change cells
**not yet born**, so its time constant is the lifespan of the new *long-lived*
cohort (~120 d), not of the old short-lived one. A gene therapy trial powered on
a 24-week haemoglobin endpoint is reading its own transient.

### 7 · The benefit the haemoglobin endpoint cannot see

PK deficiency loads iron **without transfusion**: marrow expansion raises
erythroferrone, which suppresses hepcidin, which lifts the brake on duodenal
absorption. A PK activator reduces the destruction the marrow was compensating
for, so erythroferrone falls and hepcidin rises. The model reproduces every sign
van Beers 2024 measured — ERFE ↓, hepcidin ↑, EPO ↓, sTfR ↓, liver iron ↓.

The structural point: this benefit is driven by *marrow output*, so it is
**largest in patients whose haemoglobin barely moves** — the compensated ones
above the knee, whom a haemoglobin-response endpoint classifies as
non-responders.

### 8 · Why the diagnostic assay can read normal

Assayed red cell PK activity is a cohort-weighted mean, and reticulocytes carry
freshly made enzyme, so **the sicker the patient the more the assay over-reads**
(ratio 0.73 → 0.83 across the spectrum as reticulocytes rise 2 % → 26 %). This is
the mechanism behind the standard advice to interpret PK activity relative to
another age-dependent enzyme and to distrust a normal absolute value in a
reticulocytosis.

---

## Fourteen defects found by running the equations

Full log in [`pkd_reference_output.txt`](pkd_reference_output.txt) §7, each marked
at its fix site in both the Python and the R file. The instructive ones were not
coding slips but places where a physiologically plausible equation gave a
quantitatively **impossible** answer:

- **#1** A bare phosphoglycerate-kinase-equilibrium form makes 2,3-BPG ∝
  `[3-PG]·(ATP/ADP)`; a PK lesion raises the first and lowers the second and they
  cancel — predicting **+9 %** against an observed 2-fold rise. Fixed by adding
  3-PG inhibition of the 2,3-BPG **phosphatase**, so the block shuts the shunt's
  exit as well as feeding it.
- **#2** Solving for ATP and 2,3-BPG by simultaneous damped iteration **diverges**,
  because the FBP → PKR limb is positive feedback. Fixed by nesting.
- **#7** Michaelis–Menten bilirubin conjugation calibrated at the normal load
  saturates at ~3× normal haemolysis and assigned a moderately affected patient a
  bilirubin of **23 000 mg/dL**.
- **#12** Reticulocyte capture summed over two cohorts and debited to one: once
  the debited cohort hit its floor, the model **created** red cells and reported
  haemoglobins of **19–54 g/dL**. Mass balance in an age-structured model is
  per-cohort, not aggregate.
- **#13** Explicit RK4 at the burn-in step violates the stability limit (0.093 d)
  for the bilirubin state, which oscillated into its floor and reported
  **0.0 mg/dL in a haemolysing patient** — while the haemoglobin and reticulocyte
  columns of the step-convergence table still agreed to 1e-5. *A convergence
  check on slow outputs alone does not detect this.*
- **#14** With the whole hazard treated as acutely ATP-reversible, 5 mg BID raised
  haemoglobin by **+6.0 g/dL** and overshot to 22.7, against an observed mean of
  ~1.7. A cell that has already lost membrane and exported its adenylate pool is
  not rescued by refilling its ATP.

---

## Limitations — stated so they are not mistaken for results

- **pH is not a state.** The Bohr effect and the strong pH dependence of both
  BPGM activities are folded into fixed constants, so acid–base disturbance
  cannot be simulated.
- **The single most load-bearing parameter is the 3-PG inhibition constant of the
  2,3-BPG phosphatase (30 µM), and it is the least well pinned by data.** The
  oxygen-transport conclusion scales with the size of the 2,3-BPG excursion, so
  it inherits that uncertainty. The break-even table exists to let the conclusion
  be re-scored.
- Free vs Mg-bound ADP is not resolved; the PGK equilibrium uses total ADP, which
  overstates how far ATP/ADP falls.
- Splenic destruction and sequestration share one compartment with a mean-field
  transit; there is no distribution of transit times.
- The oxygen module is whole-body and single-tissue, so it cannot represent the
  regional differences in extraction where a right-shifted curve actually helps
  or hurts most.
- Cardiac output responds to haemoglobin rather than to tissue PO₂, so the model
  cannot arbitrate between the three closures of §2–3 — it reports all of them.
- The iron submodel reaches transferrin saturation in severe untreated disease;
  it captures the direction and the ranking of interventions, but the absolute
  magnitude of iron loading is bounded by that ceiling.
- **The model simulates one patient, not a distribution**, and several endpoints
  in this disease are defined on a *fraction of responders*. So it reproduces
  effect sizes but not response rates: on the ACTIVATE-T transfusion endpoint it
  predicts a 93 % reduction in units for the simulated patient, where the trial
  reported ≥33 % reduction in 37 % of participants. Turning the former into the
  latter requires a population layer (a distribution over residual activity, τ_PK
  and activatable fraction) that is not implemented here.

## Running it

```r
# mrgsolve model + 24 scenarios
source("pkd_mrgsolve_model.R")

# interactive dashboard
shiny::runApp("pkd_shiny_app.R")
```

```bash
# the executable reference: regenerates pkd_reference_output.txt from scratch
python3 pkd_reference_model.py      # no dependencies; ~6 min
```

`$ODE` solves the fast glycolytic subsystem **exactly**, for 29 cell populations,
at every derivative evaluation. That is affordable in compiled C++ but not free —
expect seconds to tens of seconds per simulated year. The Python reference
tabulates the same subsystem instead; §0.3 of the output file bounds the
difference between the two implementations.

---

> ⚠️ Educational and research model. Not validated against individual patient
> data; not for clinical decisions, prescribing, or regulatory submission.
