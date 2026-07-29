# Tobacco Use Disorder (Nicotine Dependence) — QSP Model

> **담배 사용 장애 / 니코틴 의존** — a quantitative systems pharmacology model in
> which dependence is not "nicotine leaves and the smoker feels bad" but a
> **mismatch between two processes on very different time constants**: α4β2*
> nicotinic receptor **desensitization** (t½ ≈ 7 min on / 2 h off) and the
> **allostatic hedonic set-point** (t½ ≈ 14 days). Dopaminergic *supply*
> collapses within hours of the last cigarette; the *set-point* takes weeks to
> follow. Withdrawal is the transient between two matched states, not a state
> of its own.

| File | Description |
|------|-------------|
| `tud_qsp_model.dot` | Mechanistic map — 165 nodes, 13 functional clusters |
| `tud_qsp_model.svg` | Vector render of the map |
| `tud_qsp_model.png` | 150-dpi raster preview |
| `tud_mrgsolve_model.R` | mrgsolve QSP model — **31 ODE compartments, 10 therapy scenarios, 6 focused experiments** |
| `tud_shiny_app.R` | Shiny dashboard — 9 tabs |
| `tud_references.md` | 72 PubMed-indexed references, each PMID resolved against the live PubMed API |

## Mechanistic snapshot

![Tobacco Use Disorder QSP mechanistic map](tud_qsp_model.png)

**Thirteen functional clusters:**

1. **Product & delivery kinetics** — cigarette, VLNC, ENDS, patch, gum/lozenge; arterial spike; rate-of-rise reinforcement.
2. **Nicotine disposition & CYP2A6 pharmacogenetics** — 2-compartment PK, brain biophase, cotinine, 3′-hydroxycotinine, UGT2B10/FMO3, the NMR biomarker.
3. **nAChR pharmacology** — α4β2*, α6β2β3, α5, α7, α3β4; occupancy, desensitization, resensitization, upregulation, β2* PET.
4. **Mesolimbic dopamine reward circuit** — VTA DA/GABA/glutamate, phasic bursting, NAc shell release, DAT, MAO inhibition by smoke β-carbolines.
5. **Habenula–IPN aversion & intake titration** — MHb/IPN, CHRNA5 rs16969968, aversion set-point, inverted-U self-administration ceiling.
6. **Allostasis, withdrawal & negative affect** — CRF-CRF1, dynorphin-KOR, LC noradrenaline, HPA, orexin, MNWS.
7. **Craving, cue reactivity & habit** — dorsal striatum, insula, dlPFC, mGluR5, incentive salience, QSU-brief.
8. **Smoking behaviour & lapse/relapse dynamics** — CPD, FTND, lapse hazard, abstinence survival.
9. **Pharmacotherapy PK** — varenicline, cytisinicline, bupropion + hydroxybupropion, NRT patch/PRN, nortriptyline, clonidine, nicotine vaccine.
10. **Pharmacodynamics** — partial agonism, competitive blockade, DAT/NET inhibition, non-competitive nAChR block, 5-HT3 nausea.
11. **Non-nicotine tobacco toxicology** — CO/COHb, PAH → **CYP1A2 induction** (the clozapine/olanzapine/theophylline dose trap on quitting), acrolein, nitrosamines.
12. **Organ-system consequences & safety** — HR/BP, atherosclerosis, FEV1, appetite/weight, sleep, EAGLES neuropsychiatric safety, bupropion seizure risk.
13. **Biomarkers, covariates & clinical endpoints** — cotinine, NMR, exhaled CO, NNAL, β2* PET; CAR wk 9–12, 7-day PP, ΔQSU, ΔMNWS, Δweight.

## The 31-compartment mrgsolve model

| Group | Compartments |
|-------|--------------|
| Nicotine PK | `NICC`, `NICP`, `NICB` (biophase), `ASKIN` (patch depot), `AMOUTH` (buccal depot) |
| Metabolites | `COT`, `HCOT` — CYP2A6 scales **both** oxidation steps, which is what sets the NMR |
| Varenicline | `AGV`, `VARC`, `VARB` |
| Cytisinicline | `AGC`, `CYTC` |
| Bupropion | `AGB`, `BUPC`, `OHB` (hydroxybupropion, the main active moiety) |
| Receptor | `DES` (desensitized fraction), `RUP` (α4β2* pool, relative to never-smoker) |
| Reward / allostasis | `DA`, **`SETP`** (hedonic set-point), `ALLO` |
| Symptoms | `WD` (MNWS 0–4), `QSU` (craving 1–7), `HABIT` |
| Behaviour | `PABST` (abstinence survival) |
| Downstream | `COHB`, `WT`, `FEV`, `NAUS`, `TOL`, `SLP`, `PADH` |

The structural core is three lines:

```
DA      ← fast   : tracks activated, non-desensitized receptors (minutes)
SETP    ← slow   : dSETP/dt = KSET*(DA - SETP),  t½ = 14 d
DEFICIT = max(0, SETP - DA)          ← the withdrawal driver
```

In a smoker at steady state `DEFICIT ≈ 0` **by construction** — the set-point has
adapted to whatever tone the smoking pattern delivers. That is why nothing has to
be assumed about how bad quitting feels.

**Ten therapy scenarios**: continued smoking · unaided quit · patch 21 mg ·
combination NRT (patch + PRN gum) · varenicline standard · varenicline
preloading · bupropion SR · cytisinicline · varenicline + patch · patch in a
slow CYP2A6 metabolizer.

## What the model produces

All numbers below are **computed** from the ODE system, not asserted. They come
from an independent LSODA integration of the identical equations and are recorded
in section 7 of the R file so drift is detectable.

### Baseline smoker steady state (20 cig/day, normal CYP2A6, 90-day run-in)

| Quantity | Model | Literature |
|---|---|---|
| Plasma nicotine (24-h mean) | 12.7 ng/mL | 10–17 ng/mL |
| Cotinine | 194 ng/mL | 150–300 ng/mL |
| NMR (3HC/cotinine) | 0.39 | 0.20–0.45 |
| α4β2* occupancy (daytime) | 0.94 | >0.88 after one cigarette (Brody 2006) |
| β2* upregulation | **+70%** | +30 to +100% by PET |
| Exhaled CO | 25 ppm | 20–30 ppm |

### 1. Occupancy calibration is exact, and it is the reason the rest works

`KDNIC` = 5.4 nM was set so that **0.87 ng/mL gives 50% occupancy** — Brody's
reported EC50. The consequence is that occupancy is **saturated at every ordinary
exposure**: 2 ng/mL → 70%, 10 ng/mL → 92%, 20 ng/mL → 96%. Occupancy therefore
*cannot* distinguish a patch from a cigarette, which is why the model needs an
explicit route-resolved **phasic** term. This is a constraint the calibration
imposed, not a modelling preference.

### 2. Withdrawal timescales emerge from two rate constants

Nothing in the model is told when withdrawal peaks or how long it lasts.

| Readout | Model | Clinical |
|---|---|---|
| MNWS peak day | 5.2 | 1–3 |
| MNWS above ½ peak | 35 d | ~14–28 d |
| β2* pool: t50 → t90 | 15.0 d → 47.5 d | normalizes 3–4 wk (Cosgrove 2009) |
| Set-point 1.86 → wk 4 → wk 12 | 1.22 → 1.01 | — |

**The model runs roughly 1.5–2× slow.** The shape and the mechanism are right;
`KSET` and `KOUTA` are not fully calibrated. Reported rather than tuned away.

### 3. Fixed-dose NRT + nicotine titration reproduces the NMR result — with no interaction term

One scalar `F2A6` sets clearance. Intake titrates against clearance, so fast
metabolizers smoke *more* yet run *lower* plasma nicotine — and a **fixed** 21 mg
patch replaces a shrinking fraction of what they were used to:

| CYP2A6 | cig/day | NMR | pre-quit nicotine | patch wk 4 | **replacement ratio** | patch CAR | varenicline CAR |
|---|---|---|---|---|---|---|---|
| slow (0.35) | 14.3 | 0.14 | 17.6 ng/mL | 14.9 | **0.85** | **31.0%** | 37.7% |
| normal (1.0) | 20.0 | 0.39 | 12.7 | 7.3 | **0.57** | **27.0%** | 35.8% |
| very fast (2.0) | 26.5 | 0.78 | 9.6 | 4.0 | **0.42** | **23.3%** | 34.2% |

Patch efficacy falls **7.7 points** across the CYP2A6 range; varenicline —
renally cleared, CYP2A6-independent — moves only 3.5 points, and that residual is
the *shared dependence gradient*, not a drug effect. This is the Lerman 2015
treatment-matching interaction, produced by pharmacokinetics alone.

### 4. Partial agonism is **not** an efficacy optimum — and the model says why

Sweeping varenicline intrinsic activity 0 → 1:

| EMAXV | DA wk4 | MNWS peak | lapse reward | CAR wk 9–12 |
|---|---|---|---|---|
| 0.00 (antagonist) | 1.000 | 1.56 | 0.227 | 32.8% |
| 0.45 (varenicline) | 1.160 | 1.43 | 0.227 | 35.8% |
| 1.00 (full agonist) | 1.328 | 1.27 | **0.227** | **39.4%** |

CAR rises **monotonically** and the blockade term is **pinned at 0.227**, because
blockade depends on *occupancy* and occupancy is independent of intrinsic
activity. So the model **refutes** the common claim that partial agonism is
optimal for efficacy, and relocates the argument: partial activity plus
non-pulsatile oral delivery is what limits the drug's own reinforcing value —
something this model does not score. Consistent with the same arithmetic, adding
a patch (full-agonist activity) to varenicline lowers peak MNWS 1.43 → 1.31.

### 5. Reduced-nicotine cigarettes do nothing until nicotine falls very far

Because occupancy saturates, cutting nicotine per cigarette has almost no effect
until it drops by more than an order of magnitude:

| mg nicotine/cig | 1.10 | 0.55 | 0.25 | 0.10 | 0.050 | 0.030 | 0.015 |
|---|---|---|---|---|---|---|---|
| occupancy | 0.94 | 0.89 | 0.79 | 0.60 | 0.43 | 0.31 | 0.18 |
| receptor pool `RUP` | 1.70 | 1.70 | 1.69 | 1.66 | 1.63 | 1.59 | 1.52 |
| MNWS peak | 1.58 | 1.54 | 1.48 | 1.33 | 1.15 | 0.97 | 0.71 |
| CAR wk 9–12 | 12.5% | 13.1% | 14.3% | 17.2% | 20.6% | **23.8%** | 28.2% |

Halving nicotine buys 0.6 percentage points. The threshold sits **below ~0.10
mg/cigarette** — which is roughly where the US reduced-nicotine product standard
(0.4 mg/g ≈ 0.03 mg/cig, Donny 2015) actually lands. *Caveat: compensatory
smoking is not modelled, so these are upper bounds.*

### 6. A falsifiable prediction: varenicline + patch is not additive

| Arm | total occupancy | **varenicline's own** | nicotine's | MNWS peak | CAR |
|---|---|---|---|---|---|
| varenicline alone | 0.897 | **0.897** | 0.000 | 1.43 | 35.8% |
| varenicline + patch | 0.944 | **0.501** | 0.443 | 1.31 | 33.0% |

Adding a patch raises total occupancy by 5% but **nearly halves varenicline's own
occupancy** — and it is that term which carries the partial-agonist-specific
blockade. The model predicts the combination is roughly equal to, not better
than, varenicline alone. This **disagrees with Koegelenberg 2014** and agrees
with the later neutral trials, and it names the measurement that would settle it:
*measure varenicline β2* occupancy by PET with and without a concurrent patch. If
occupancy is not displaced, the model is wrong here.*

### Validation against trial abstinence rates

Behavioural hazard weights (`HAZ0`, `B1QSU`, `B2WD`, `B3HAB`, `GBLOCK`,
`GBLOCKV`, `GREINF`) were fitted to these seven arms and are marked `[FITTED]` in
`$PARAM`. **Everything upstream of the hazard** — PK, metabolism, occupancy,
desensitization, upregulation, dopamine, set-point — is fixed from mechanism and
never fitted to an abstinence rate.

| Arm | model | trial | Δ |
|---|---|---|---|
| Unaided quit (placebo) | 12.5% | 12.5% | fit target |
| Nicotine patch 21 mg | 27.0% | 23.4% | +3.6 |
| Combination NRT | 25.1% | 28.0% | −2.9 |
| Varenicline 1 mg BID | 35.8% | 33.5% | +2.3 |
| Bupropion SR | 22.2% | 22.6% | −0.4 |
| Cytisinicline 3 mg TID | 33.8% | 32.0% | +1.8 |
| Varenicline + patch | 33.0% | 38.0% | −5.0 |

**RMSE 2.80 percentage points.** Anchors: EAGLES (Anthenelli 2016), ORCA-2
(Rigotti 2023), Cochrane NRT reviews, Koegelenberg 2014.

## Known discrepancies

Stated, not tuned away:

1. **Withdrawal resolves ~1.5–2× slower** than the clinical 2–4 weeks.
2. **Combination NRT lands just below patch alone** (25.1 vs 27.0%) whereas Cochrane finds it clearly better. Cause: the fitted phasic gain `GREINF` makes PRN oral nicotine mildly hazard-increasing, which here outweighs its craving relief.
3. **Varenicline + patch predicted non-additive** — see above; a genuine disagreement with one trial.
4. **Bupropion relieves withdrawal more than the literature supports** (peak MNWS 1.06 vs 1.58 unaided) because `GDAT` enters the dopamine drive multiplicatively. Its CAR is nonetheless on target.
5. **`WT` is a population mean** over abstainers and relapsers, so the 52-week value (+2.0 kg) sits well below the ~5 kg seen in sustained abstainers; a sustained abstainer here reaches `WTSS` = 5 kg.
6. **Rate of rise is not resolved kinetically.** Smoking is a diurnal continuous input, so the arterial spike is represented as a per-product attribute (`FPHCIG` / `FPHORAL` / `FPHPATCH`) scaled by reinforcement-episode frequency — a modelling choice, not a derivation.
7. **Trajectories are conditional on remaining abstinent.** `PABST` does not feed back into nicotine input, so relief in relapsers is under-estimated.

## Shiny dashboard (9 tabs)

1. **Patient profile** — phenotype, CYP2A6 genotype, occupancy-calibration curve
2. **Nicotine PK** — nicotine / cotinine / 3HC / drug concentrations; diurnal zoom showing the overnight trough that produces time-to-first-cigarette craving
3. **Receptor pharmacology** — occupancy by ligand, desensitization, receptor pool, activation; live EMAXV sweep
4. **Reward & withdrawal** — DA vs set-point vs deficit, MNWS, QSU, habit
5. **Clinical endpoints** — abstinence survival curve + validation against trial CARs
6. **Scenario comparison** — arms side by side + the displacement analysis
7. **Biomarkers** — NMR treatment matching, cotinine, exhaled CO, β2* pool
8. **Safety & tolerability** — nausea/tolerance, sleep, weight, FEV1, persistence
9. **About** — structure, calibration anchors, disclaimer

## Running it

```r
# ODE model, all scenarios and experiments
Rscript tud_mrgsolve_model.R

# interactive dashboard
R -e 'shiny::runApp("tud_shiny_app.R")'

# re-render the mechanistic map
dot -Tsvg tud_qsp_model.dot -o tud_qsp_model.svg
dot -Tpng -Gdpi=150 tud_qsp_model.dot -o tud_qsp_model.png
```

Requires `mrgsolve`, `dplyr`, `tidyr`, `ggplot2`, `shiny`, and Graphviz.

## ⚠️ Disclaimer

Educational and research QSP model. Built from public literature and published
trial data, but **not independently validated or certified**. It must **not** be
used for clinical decision-making, prescribing, or regulatory submission.
Parameters are illustrative approximations; the behavioural hazard weights are
fitted to aggregate trial results and carry no individual predictive validity.
See [`tud_references.md`](tud_references.md) for the evidence behind every
parameter.
