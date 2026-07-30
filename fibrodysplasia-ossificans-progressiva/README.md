# Fibrodysplasia Ossificans Progressiva (FOP, 진행성 골화성 섬유이형성증) — QSP Disease Model

> **Fibrodysplasia ossificans progressiva** is an ultra-rare autosomal-dominant
> disorder in which a single recurrent mutation in **ACVR1/ALK2** (c.617G>A,
> p.R206H, >95% of classic cases) turns **activin A** — normally a
> *non-signalling competitive inhibitor* of wild-type ALK2 — into a full
> agonist. Resident **fibro/adipogenic progenitors** in skeletal muscle
> connective tissue respond by running the complete **endochondral ossification**
> programme in the wrong place, converting muscle, tendon, ligament and fascia
> into mature, marrow-containing bone. The bone is never removed. Patients
> accumulate an irreversible skeleton-outside-the-skeleton, lose one joint after
> another, and die at a median age of **56 years**, most often of
> **cardiorespiratory failure from thoracic insufficiency syndrome**.

---

## What makes this model different from a generic disease model

Almost every mechanistic map of FOP ends at "heterotopic bone". That is where
the interesting problem *starts*, because FOP has an unusual property: **the
disease is a stock, and every drug and every trial endpoint is a rate.** The
whole model is organised around four consequences of that, and each one is
**derived from published numbers rather than asserted**.

### Axis 1 — the anti-flare ceiling: φ_FL = 0.30

The 36-month natural history study (NHS, PMID 36152026) reports four numbers in
one cohort. Multiplying three of them against the fourth settles a question the
paper never asks:

| quantity | value |
|---|---|
| flare-up rate | 229 flares / 114 patients / 26.8 mo = **0.90 per patient-year** |
| imaged flares with new HO at that site by day 84 | 14/52 = **26.9%** |
| mean new HO volume at those sites | **28.8 mL** |
| product = flare-attributable new HO | **6.97 mL/yr** |
| total mean annualised new HO | **23.6 mL/yr** |
| **φ_FL = flare-attributable fraction** | **6.97 / 23.6 = 0.295** |

**A drug that abolished 100% of flare-ups would remove ~30% of new heterotopic
bone and no more.** Two independent corroborations: 47% of patients report
progression *without* obvious flare-ups (PMID 27025942), and in LUMINA-1 placebo
only 12/29 (41%) of new lesions were co-located with a flare (PMID 39216107).
Most heterotopic bone in FOP forms quietly, without a flare to announce it.

The model therefore runs **two parallel routes** — episodic and smouldering —
and keeps them *exactly* separable by holding the receptor step linear in
activin A (justified: FOP serum activin A is ~434 pg/mL, PMID 37165433, far
below receptor saturation). Glucocorticoids, NSAIDs, mast-cell stabilisers and
IL-1 blockade enter only the episodic route, so they cannot beat φ_FL — **by
construction, which is the point**. Simulated: perfect abolition of every flare
gives 28.8%; real-world glucocorticoid bursts give 3.4%.

### Axis 2 — the ratchet: the rate falls 60%, the burden falls 11%

Heterotopic bone remodels internally but is not resorbed (`KRESORB = 0`), and
surgical excision provokes explosive recurrence. So the trial endpoint
(annualised new HO volume) is a **slope** while disability is the **area under
it**. One simulation of the MOVE regimen produces both:

| from the same 18-month run, age 15 | untreated | palovarotene |
|---|---|---|
| annualised new HO (mL/yr) — *the trial endpoint* | 37.3 | 14.8 (**−60.3%**) |
| **total** HO stock at 18 months (mL) — *the disability* | 306 | 272 (**−11.0%**) |

Both numbers are correct descriptions of the same trial. A trial that reports
only the first is not wrong; it is answering a different question from the one
the patient is asking.

### Axis 3 — age beats potency, and the reason is *competence*, not territory

The NHS velocity profile is non-monotonic: **21.9 mL/yr** at 8–<15, **41.5** at
15–<25, then a collapse to **4.6** at 25–65. The model reproduces all four band
values (see below) from a somatic growth curve (anchored to human growth, *not*
fitted to HO data) multiplied by a chondrogenic-competence factor that falls at
skeletal maturity.

Two candidate mechanisms for the adult collapse were tested and one was
**rejected by the fit**:

- *Territory exhaustion* — adults have run out of ossifiable tissue. A
  territory-limited model **cannot** reproduce the adolescent peak while keeping
  adults progressing; forcing it costs 18% on the 15–<25 velocity.
- *Loss of chondrogenic competence at skeletal maturity* — fits all four bands.
  It also explains, with the same parameter, why the physeal hazard of a retinoid
  disappears once the plates fuse. **One mechanism, two independent
  observations.**

This distinction is clinically decisive: under competence decline adults are
*not* out of substrate, so an adult who is operated on or injured can still
ossify catastrophically — which is what the clinical record shows, and what a
territory-exhaustion model would wrongly call safe.

The consequence for treatment, simulated across five starting ages:

| start age | HO already formed | % of age-56 burden still ahead | % of the age-40 burden averted | physis |
|---|---|---|---|---|
| 5 | 60 mL | 92.4% | **47.9%** | OPEN — PPC risk |
| 8 | 91 mL | 88.5% | **46.0%** | OPEN — PPC risk |
| 15 | 250 mL | 68.2% | **35.2%** | fused |
| 25 | 645 mL | 18.1% | **5.8%** | fused |
| 35 | 694 mL | 11.8% | **1.8%** | fused |

Note also what the model **refuses** to report: under any saturating model,
"lifetime HO averted" is the wrong endpoint, because slowing formation changes
*when* a burden is reached, not the final total. The honest currencies are
burden at a fixed age and **years of delay**.

### Axis 4 — efficacy and toxicity share one exposure axis

Palovarotene works by blocking chondrogenic commitment. **The growth plate *is*
chondrogenesis.** MOVE delivered 54–60% reduction in new HO *together with*
premature physeal closure in **21/57 (36.8%)** of children under 14, reduced
vertebral BMD and increased vertebral fracture risk (PMID 36583535); daily
dosing ablates growth plates in mice (PMID 30226468).

The model does not assume the two effects are separable — it **solves for the
implied selectivity** from those two published numbers and gets
**R_SEL = EC50(PPC)/EC50(HO) ≈ 3.7-fold** (conditional on the assumed maximal
physeal hazard, which is stated as a parameter). Then it sweeps the dose:

| palovarotene | 0.5 mg | 1 mg | 2 mg | 3 mg | **5 mg** | 7.5 mg | 10 mg |
|---|---|---|---|---|---|---|---|
| new HO suppressed | 15.4% | 25.8% | 39.1% | 47.4% | **57.2%** | 63.9% | 67.9% |
| physeal closure risk | 7.0% | 12.7% | 21.5% | 27.9% | **36.8%** | 43.7% | 48.2% |

At ~3.7-fold selectivity there is **no dose** that reaches 50% suppression at an
acceptable physeal hazard while the plates are open. The only separation axis is
skeletal maturity — i.e. *time* — and by then Axis 3 has cut the remaining
benefit to about a tenth. The `R_SEL` slider in the Shiny app answers the design
question directly: **a next-generation retinoid needs roughly an order of
magnitude more physeal selectivity to be usable in the age group that needs it.**

---

## Two trials explained, and one trial design questioned

**LUMINA-1 (garetosmab) reported a failed primary endpoint and a striking
secondary one.** Total lesion activity by PET-CT missed (p = 0.0741) while new
HO lesion development was suppressed (0% vs 40.9%, p = 0.0027). The model
produces both from one mechanism: **total** NaF activity sums the new
mineralising front *and* ongoing remodelling of the mature stock the drug cannot
touch. In the model's LUMINA-1 cohort only 39% of the NaF signal comes from new
lesions, so a 91% reduction in new lesions shows up as a **38% reduction in
total lesion activity** — a real effect on a diluted endpoint, in a trial of 20
versus 24 patients.

**The LUMINA-1 placebo arm formed bone ~12× faster than an NHS-average adult**
(1.21 new lesions per patient in 28 weeks versus ~0.1). The model cannot
reproduce the trial's event rate from population parameters without an explicit
12-fold activity enrichment — a quantitative measure of how selected that cohort
was, and a warning about the design MOVE used.

**MOVE was single-arm against untreated NHS controls.** Under this model the
apparent effect is (drug effect) × (intrinsic activity ratio between the
cohorts), and the two are not separable from the data:

| trial cohort activity relative to the external control | 1.00 | 0.80 | 0.60 | **0.40** | 0.25 |
|---|---|---|---|---|---|
| apparent reduction versus the external control | 57.0% | 65.6% | 74.1% | 82.7% | 89.2% |
| true drug effect | 57.0% | 57.1% | 57.1% | 57.2% | 57.2% |
| **apparent reduction with NO drug at all** | 0.0% | 19.8% | 39.7% | **59.7%** | 74.8% |

A cohort with 40% of the external control's intrinsic activity would show a
**~60% "reduction" with no drug whatsoever** — the entire MOVE effect size,
manufactured by selection. This does not show that palovarotene is inactive; it
shows the design cannot separate the two, which is exactly why the MOVE analysis
was so sensitive to a square-root transformation (99.4% versus 65.4% posterior
probability). The model reproduces that fragility too: with between-subject
variability the same simulated trial gives a **50.5% reduction in the mean but
30.1% after a square-root transform**, at a coefficient of variation of 0.92.

---

## Deliverables

| File | Contents |
|------|----------|
| `fop_qsp_model.dot` / `.svg` / `.png` | Graphviz mechanistic map — **21 clusters, 130 nodes, 203 edges** |
| `fop_mrgsolve_model.R` | **42-ODE** mrgsolve model (15 drug PK + 27 disease/clinical), **15 scenarios**, calibration report |
| `fop_shiny_app.R` | **11-tab** Shiny dashboard (all 21 outputs smoke-tested under `shiny::testServer`) |
| `fop_references.md` | **151 PubMed references**, every PMID resolved and title-verified through the E-utilities API |

**Reproduce everything:**

```bash
dot -Tsvg fop_qsp_model.dot -o fop_qsp_model.svg
dot -Tpng -Gdpi=150 fop_qsp_model.dot -o fop_qsp_model.png
Rscript fop_mrgsolve_model.R          # runs all 15 scenarios + calibration table
Rscript -e 'shiny::runApp("fop_shiny_app.R")'
```

---

## Model structure — 42 ODE compartments

**Drug PK (15)** — palovarotene (oral, 2-cpt, CYP3A4), garetosmab (IV, 2-cpt +
target-mediated disposition), ALK2 kinase inhibitor (fidrisertib/zilurgisertib
class), saracatinib, sirolimus, prednisone, anakinra.

**Signalling (12)** — trigger, inflammatory amplifier (mast cells, macrophages,
IL-1), HIF-1α, mTORC1, free activin A **split into smouldering and episodic
pools**, garetosmab–activin complex, pSMAD1/5/8 **split the same way**, activated
FAPs, fibroproliferative lesion, chondrogenic anlagen.

**Lesion maturation and bone (7)** — hypertrophic/vascular-invasion stage,
osteoblastic osteoid, mineralised thoracic HO, mineralised appendicular HO,
mature-stock NaF remodelling pool, cumulative new HO, cumulative lesion count.

**Clinical and safety (8)** — CAJIS, FVC % predicted, cumulative mortality
hazard, cumulative physeal-closure hazard, vertebral BMD Z, height Z,
mucocutaneous retinoid AE burden, cumulative flare count.

### Where each mechanism class caps out — computed from parameters, not from a dose

| mechanism class | analytic ceiling | set by |
|---|---|---|
| anti-inflammatory / anti-flare | **29.5%** | φ_FL, the flare-attributable fraction |
| RARγ commitment blockade | **80.6%** | the lesion-fate competition `KRESOLVE/KCOMMIT` |
| anti-activin A (ligand) | **93.0%** | `F_LEAK`, the ligand-independent receptor leak |
| ALK2 kinase inhibition | **95.0%** | `EMAX_AK` — the only class that reaches the leak |

The classes differ in **where they cap out**, not only in how potent they are.
Because R206H loses the FKBP12 clamp, some signal is ligand-independent, so no
anti-ligand antibody can reach it however completely it neutralises activin A.

### The lesion-fate checkpoint is a competition, not a switch

`KCOMMIT/(KCOMMIT + KRESOLVE) = 0.269` **is** the published day-84 ossification
probability of an imaged flare (14/52). A commitment blocker does not abolish
ossification — it shifts lesions from the ossifying branch into the resolving
branch. That is why 85% maximal blockade yields an 80.6% ceiling, and why
flare-triggered-only dosing achieves just **14%** against chronic dosing's
**57%**: by the time a flare is symptomatic, part of the lesion has already
crossed the checkpoint.

---

## Calibration — model versus published

| anchor | published | model |
|---|---|---|
| annualised new HO, 8–<15 yr (mL/yr) | 21.9 | 22.8 |
| annualised new HO, 15–<25 yr (mL/yr) | 41.5 | 39.4 |
| annualised new HO, 25–65 yr (mL/yr) | 4.6 | 4.6 |
| flare-attributable fraction of new HO | 0.295 | 0.295 |
| HO stock at 2–<8 yr (mL) | 68.8 | 68.6 |
| per-flare ossification probability | 0.269 | 0.269 |
| flare-up rate (per patient-year) | 0.90 | 0.90 |
| palovarotene reduction in annualised new HO | 54–60% | 60.3% |
| premature physeal closure, 18 mo, age <14 | 36.8% | 36.8% |
| retinoid-associated AE burden | 0.970 | 0.980 |
| FVC at age 30 (% predicted) | 44 ± 14 | 43.5 |
| median survival, untreated (yr) | 56 | 56.0 |
| garetosmab C<sub>trough</sub>, 10 mg/kg q4w (mg/L) | 105 ± 31 | 116.1 |
| **HO stock at 25–65 yr (mL)** | **575.2** | **739.4 — see below** |

Fitted parameters: `KF_HO`, `CMIN_C`, `A50C_C`, `AWC_C` (to the four velocity
bands); `H0_MORT` (to median survival); `EC50_PAL` (to the MOVE reduction);
`R_SEL` (to the PPC rate); `CL_GAR` (to C<sub>trough</sub>). Everything else is
either taken from a publication or fixed a priori.

### The FVC curve was not fitted, and it reproduces an odd observation

`FVC = 35 + 65 · K/(K + HO_thoracic)` with a small `K` saturates **early**. That
was chosen because Botman et al. (PMID 33748352) found restrictive physiology
already present in childhood and, in 4 of 5 patients, **thoracic HO progressing
into adulthood with no further FVC decline**. The model reproduces this
(FVC 71% at age 8 → 47% at 20 → 43.5% at 30 → 42.8% at 45) and it carries a
sobering implication that follows from the shape alone: the pulmonary damage is
front-loaded, so preventing *later* thoracic HO buys very little lung function.

---

## Failures reported rather than repaired

1. **The NHS dataset is internally inconsistent, and the model shows by how
   much.** Integrating the study's own longitudinal velocities from its own
   childhood baseline (68.8 mL at 2–<8, then 21.9 and 41.5 mL/yr) puts a patient
   at **~650 mL by age 25** — already above the same study's cross-sectional
   25–65 band mean of **575.2 mL**. The model was fitted to the *prospectively
   measured longitudinal* velocities (the study's primary outcome and the
   quantity every trial uses as its endpoint) and consequently overshoots the
   cross-sectional adult stock by **+29%**. This is not repairable by tuning: you
   cannot have 41.5 mL/yr running through the twenties and only 575 mL by
   mid-life. Plausible causes are on the data side — cross-sectional band means
   over a 40-year age range dominated by younger participants, survivor effects,
   and unequal follow-up in the velocity estimates.
2. **Territory is not identifiable from these data, and a free fit silently
   deleted it.** Allowed to float, the optimiser drove `HOMAX` to 10^16 and
   `GAMMA_T` to 10^5 — voiding the mechanism while fitting every velocity
   perfectly. `HOMAX` is therefore **fixed a priori** at 2000 mL (about 39%
   consumed by age 56) rather than fitted, and the model's late-life plateau is
   an assumption, not a result.
3. **The competence transition sits on its lower bound.** `AWC_C` fits to 0.50 yr
   — the width limit allowed. Two adjacent NHS band means of 41.5 then 4.6 mL/yr
   demand a near-*step* loss of chondrogenic competence in the mid-twenties,
   possibly sharper than this model permits. That is a strong biological claim
   and it is the data's, not the model's.
4. **CAJIS is not numerically calibrated.** It is implemented on its true 0–30
   scale, monotone in appendicular stock, saturating in late-stage disease, and
   validated in shape only. No published cohort mean was used to pin it. Read
   CAJIS output as ordinal; the model does not claim a CAJIS of 15.9 is
   measurably distinct from 15.5.
5. **Survival and life-years are extrapolation, not evidence.** No FOP trial has
   ever measured a functional or survival endpoint. The model's projection that a
   57%-effective commitment blocker started at age 8 raises median survival from
   56.0 to 68.5 years is the model reasoning out loud about what the volume
   endpoints imply — it is a hypothesis generator, not a result.

### Defects found by actually running the model

Nine, all fixed: multi-line `:` continuations that broke `@annotated` parsing; a
signalling steady state written with the leak term on the wrong side, so the
initial conditions sat off equilibrium; an `ACTAFL_0` expression with a term
multiplied by zero; the episodic/smouldering ratio initialised as **0.74 instead
of 0.295** — inverted, which would have made the anti-flare ceiling the opposite
of the finding; `$OMEGA` random effects drawn by default, making every
"deterministic" scenario silently stochastic; interpolation at a time point
outside the output grid returning `NA` and propagating into every reported
percentage; a garetosmab TMDD V<sub>max</sub> so large it eliminated 80% of each
dose, putting C<sub>trough</sub> 75-fold below the published value; scenarios
started at age 15 or 33 while carrying a 5-year-old's 60 mL of heterotopic bone,
which corrupted every percent-of-total-volume statement (Axis 2); and an adult
LUMINA-1 cohort rendered completely inert by territory exhaustion, producing
zeroes in every cell of that scenario.

---

## ⚠️ Disclaimer

This is an **educational and research-grade QSP model**, built from public
literature and not independently validated. It must not be used for clinical
decisions, prescribing, or regulatory submission. Parameters are illustrative
approximations; fitting and qualification against patient-level data would be
required for any applied use. Where the model extrapolates beyond published
endpoints — survival, life-years, CAJIS trajectories — it says so.
