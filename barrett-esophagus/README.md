# Barrett's Esophagus (BE) — QSP Model

> Integrated Quantitative Systems Pharmacology model of Barrett's esophagus and
> its progression to esophageal adenocarcinoma: acid **and bile** reflux →
> epithelial injury and oxidative DNA damage → NF-κB / IL-6-STAT3 inflammation →
> bile-acid-driven **CDX2 / SOX9 / HNF4A** reprogramming to intestinal
> metaplasia → a stepwise clonal cascade (**CDKN2A/p16 loss** early →
> **TP53 mutation** as the LGD→HGD gatekeeper → whole-genome doubling /
> aneuploidy) → LGD → HGD → invasive EAC — with the full chemoprevention and
> eradication stack (esomeprazole, vonoprazan, aspirin, UDCA, RFA, EMR/ESD,
> antireflux surgery) and the endpoints the trials actually measured.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`be_qsp_model.dot`](be_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`be_qsp_model.svg`](be_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`be_qsp_model.png`](be_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`be_mrgsolve_model.R`](be_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`be_shiny_app.R`](be_shiny_app.R) |
| 📚 References (74)        | [`be_references.md`](be_references.md) |

---

## 1. The question this model is built to answer

Barrett's esophagus is the only known precursor of esophageal adenocarcinoma,
and the standard of care for the last three decades — a proton pump inhibitor,
plus periodic endoscopy — treats it as an **acid** problem. If that framing were
right, profound acid suppression should nearly abolish progression. It does not.
The **AspECT** trial (2557 patients, median 8.9 years, the largest
chemoprevention RCT ever run in this disease) found a hazard ratio of **0.73**
for high-dose esomeprazole against low-dose, not 0.1. Aspirin alone did nothing
detectable (HR 0.93). But high-dose PPI **plus** aspirin reached **HR 0.59**.

So the model is built to answer a mechanistic question: *why is acid suppression
only partly effective, and why does a COX inhibitor add something a PPI cannot?*

## 2. The central tension the model encodes

**Acid suppression is not reflux suppression.** Three explicit terms carry this:

1. **`K_PH_BILE = 0.35`, saturating in pH.** A PPI removes the H⁺ arm of the
   refluxate and inactivates pepsin above pH 4, but it does not stop the
   *volume* reflux that delivers bile acids to the esophagus. And deoxycholate
   is **more** membrane-permeant and a **stronger** ROS/NF-κB→CDX2 inducer at
   neutral pH than in acid (Jenkins 2007 *Carcinogenesis* 28:136; Nehra 1999
   *Gut* 44:598). In the model `BILE` does not fall on a PPI at all, and the
   *effective* bile toxicity **rises**. The amplification is a Hill function of
   acid blockade rather than linear, because bile-acid ionization is sigmoid in
   luminal pH: the penalty is essentially fully incurred once the patient is
   above pH 4 for much of the day, and barely grows when the dose is doubled.
2. **`K_GAS_KI = 4.0`.** Profound acid suppression raises gastrin 2-4×, and
   gastrin is trophic to Barrett epithelium through CCK2R (Haigh 2003
   *Gastroenterology* 124:615) — set so a doubling of gastrin raises Ki-67 by
   20%. Acid control therefore buys a proliferative penalty.
3. **`ACID_POW = 1.99`.** Acid injury is *convex* in acid exposure time: damage
   is done by long, deep exposures, and pepsin is only proteolytic below pH 4.
   The consequence is a **saturating benefit curve** — most of the achievable
   acid benefit is already realized by low-dose PPI, so doubling the dose adds
   little. This is the single parameter that carries AspECT's HR of 0.73 rather
   than 0.1, and it was fitted to that number.
4. **Proliferation is split into a COX-*independent* arm and a COX-*dependent*
   one** (`K_INF_KI` = IL-6/STAT3, 35% of baseline Ki-67; `K_PGE2_KI` = PGE2,
   only 2%). Aspirin can only reach the second. That split is *why* the model
   says aspirin's effect must be modest — and it is consistent with the
   negative celecoxib trial (Heath 2007) as well as AspECT's non-significant
   aspirin main effect.

Set `K_PH_BILE` and `K_GAS_KI` to zero and PPI monotherapy drops to a ratio of
0.186 — better than any arm of the real trial (§5.5 shows the full dissection).
Those two parameters are each set from their own mechanistic literature, not
tuned to hit AspECT, and together they land the model on 0.750 against a
measured 0.73. They are the first thing a skeptical reader should attack.

## 3. Mechanistic map — 153 nodes, 250 edges, 15 clusters

1. Host predisposition & risk factors (male sex, age, central obesity, smoking, GWAS loci MHC/CRTC1/BARX1/FOXP1, *H. pylori* absence as a protective factor, hiatal hernia, GERD duration)
2. Antireflux barrier & reflux physiology (LES tone, TLESRs, crural diaphragm, parietal-cell acid output, pepsin, bile-acid pool → deoxycholate, DGER, peristaltic clearance, salivary bicarbonate, nocturnal acid breakthrough)
3. Luminal insult & barrier failure (AET, DeMeester score, Bilitec bile exposure, dilated intercellular spaces, claudin loss, permeability → a vicious cycle back onto exposure, erosive esophagitis)
4. Oxidative stress & genome damage (ROS/NOX5, iNOS→peroxynitrite, 8-OHdG, double-strand breaks, telomere shortening, APOBEC signature, ER stress, bile-acid mitochondrial injury)
5. Mucosal inflammatory signaling (NF-κB p65, IL-6→STAT3, IL-8, TNF-α, IL-1β/NLRP3, COX-2→PGE2→EP2/EP4)
6. Metaplastic reprogramming (BMP4→pSMAD1/5/8, SOX9, CDX1/CDX2, HNF4A, KLF4/5, SOX2/TP63 squamous-program loss, MUC2/TFF3 goblet differentiation, cell-of-origin candidates, Notch suppression)
7. Barrett segment phenotype (specialized intestinal metaplasia, Prague C&M, SSBE vs LSBE, metaplastic surface area, buried glands, islands and tongues)
8. Clonal evolution & genomic instability (CDKN2A/p16 loss and clonal fields, TP53 mutation/17p LOH, aberrant p53 IHC, whole-genome doubling, aneuploidy, SMAD4 as EAC-restricted, ERBB2/EGFR/KRAS/MYC amplification, TERT)
9. Histologic progression cascade (NDBE → indefinite → LGD ⇄ regression → HGD → intramucosal → invasive EAC → nodal metastasis) with the published per-year rates annotated on the nodes
10. Adiposity & growth-factor axis (visceral adipose, leptin/LepR, adiponectin loss, insulin/IGF-1, PI3K→AKT→mTOR, adipose IL-6, intra-abdominal pressure as a *separate* mechanical route to reflux)
11. Proliferation ⇄ apoptosis balance (Ki-67, cyclin D1/CDK4-6, caspase-3, Bcl-xL/survivin, crypt fission, mutation-fixation rate)
12. Acid-suppression pharmacology (enteric depot, CYP2C19 phenotype, covalent sulfenamide binding, pump resynthesis, pH>4 holding time, hypergastrinemia → CCK2R, long-term PPI safety signals, vonoprazan P-CAB)
13. Aspirin/NSAID chemoprevention (irreversible COX acetylation, mucosal COX recovery via turnover, PGE2 fall, apoptosis restoration, GI bleeding, statin and metformin adjuncts)
14. Endoscopic eradication & surveillance (Seattle protocol and its sampling error, WATS-3D, Cytosponge-TFF3, p53 IHC, TissueCypher, EMR/ESD, RFA, cryoablation, neosquamous re-epithelialization, CE-D/CE-IM, recurrence, stricture, fundoplication, baclofen, UDCA, alginate)
15. Clinical endpoints (heartburn, regurgitation, dysphagia, EAC incidence, AspECT composite, annual progression hazard, stage at diagnosis, QoL including cancer worry, surveillance cost-effectiveness)

## 4. mrgsolve model — 32 ODE compartments

| Group | Compartments |
|---|---|
| Esomeprazole PK/PD | `PPI_DEPOT`, `PPI_CENT`, `PUMP` (covalent inactivation vs 50-h resynthesis), `GASTRIN` |
| Aspirin PK/PD | `ASA_DEPOT`, `ASA_CENT`, `COX` (irreversible acetylation vs turnover recovery) |
| Other agents | `VPZ_DEPOT`, `VPZ_CENT` (P-CAB), `UDCA_CENT` |
| Procedures | `ABL_PULSE` (one RFA session), `EMR_PULSE` (one EMR/ESD) |
| Reflux/injury axis | `ACID` (AET %), `BILE`, `INJURY` (8-OHdG surrogate), `INFLAM`, `PGE2` |
| Metaplasia | `CDX2`, `BE_LEN` (Prague M, cm), `KI67`, `NEOSQ`, `STRICT` |
| Clonal cascade | `FP16`, `FLGD`, `FP53`\*, `FHGD`, `ANEUPL`, `EAC` |
| Host | `ADIPO` (visceral adiposity, weight-loss-responsive) |
| Endpoints | `SYMPT`, `CUM_EAC`, `CUM_COMP` (AspECT composite integrator) |

\* `FP53` is a fraction of the **dysplastic** clone pool, not of the whole
segment — see §10.

Two design choices worth flagging:

- **The mutational drive is a product, not a sum.**
  `MUT = (INJURY/INJURY_ref) × (KI67/KI67_ref)`. Mutations are fixed only when a
  damaged cell divides, so an intervention that halves damage *or* halves
  proliferation halves the drive, and one that does both multiplies the two
  reductions. Note what this does **not** buy: because the two arms enter
  multiplicatively, the model is exactly multiplicative in the drug effects and
  therefore **cannot** generate super-additivity. It predicts
  0.729 × 0.929 = 0.677 for the AspECT combination arm and delivers 0.675,
  against a trial estimate of 0.59. The gap is the model's one honest miss and
  a real statement: either the trial's combination estimate sits at the
  favourable end of its interval, or there is an acid × COX interaction this
  model does not contain.
- **`MUT` is normalized to a fixed reference patient, not to each patient's own
  baseline.** This matters more than it looks. Normalizing per-patient would
  make every covariate — obesity, smoking, baseline acid and bile exposure —
  cancel out of the risk calculation, so a lean non-smoker and an obese smoker
  with the same segment length would progress identically. (That bug was in the
  first draft of this model and the obesity anchor caught it: obesity came out
  *protective*.)
- **Post-ablation recurrence is not a fitted parameter.** There is no
  "recurrence rate" in the model. After RFA the segment's *target* length is
  still set by the patient's own CDX2 drive, so `SEED_BURIED` (residual/buried
  glands) regrows the segment whenever acid **and** bile exposure remain
  uncontrolled. Durability is therefore a reflux-control result, not an
  ablation-technique result — a prediction, not an assumption.

## 5. Verified simulation results

Every number below is a model **output**. The full 32-state system was
integrated for 32 arms over a 10-year horizon (RK4, dt = 0.01 d, 731 output
points per arm, 331 s), and the calibration itself was done numerically against
six independent published anchors — see the calibration record at the bottom of
[`be_mrgsolve_model.R`](be_mrgsolve_model.R).

**Numerical health check:** 32 arms × 32 states — 0 NaNs, 0 negative states,
0 fraction-type states exceeding 1.0. **PASS.**

### 5.1 Acid-suppression pharmacodynamics (steady state, year 5)

| Regimen | Active pump | Acid blockade | pH>4 holding | AET (% pH<4) | Gastrin | Symptom score |
|---|---|---|---|---|---|---|
| No therapy | 1.000 | 0% | 0% | 12.52% | 1.00× | 6.78 |
| Esomeprazole 20 mg OD | 0.491 | 50.9% | 50.9% | 5.66% | 2.02× | 6.05 |
| Esomeprazole 40 mg BID | 0.181 | 81.9% | 81.9% | 2.12% | 2.89× | 5.54 |
| Vonoprazan 20 mg OD | 1.000 | 40.6%\* | 40.6%\* | 5.15% | 2.18× | 5.95 |

Off-therapy AET of 12.5% falls to 5.7% on 20 mg OD and 2.1% on 40 mg BID — the
dose ladder reproduces the pH-metry literature (~51% vs ~82% of 24 h above
pH 4), and only the high-dose arm brings AET under the 4% normal threshold.
Gastrin rises 2.0-2.9×, matching the long-term-PPI literature.

\* Vonoprazan's `PUMP` stays at 1.000 because it is a reversible
K⁺-competitive blocker, not a covalent one; the blockade figure is the
instantaneous trough value, while the daily-averaged effect is what drives
`AET` down to 5.15%.

### 5.2 The AspECT 2×2 factorial (NDBE archetype, 10 years)

| Arm | Cumulative HGD/EAC | Ratio vs low-dose PPI | AspECT HR |
|---|---|---|---|
| No therapy | 21.94% | 3.29 | — (not randomized) |
| Esomeprazole 20 mg OD | 6.66% | 1.000 (ref) | 1.00 (ref) |
| Esomeprazole 40 mg BID | 5.00% | **0.750** | **0.73** |
| Eso 20 mg OD + aspirin | 6.19% | **0.929** | **0.93 (NS)** |
| Eso 40 mg BID + aspirin | 4.63% | **0.695** | **0.59** |
| Vonoprazan 20 mg OD | 6.02% | 0.904 | — |
| Eso 40 mg BID + UDCA | 3.18% | 0.477 | — |
| Antireflux surgery at y1 | 2.33% | 0.349 | — |

Absolute EAC incidence in the low-dose-PPI arm is **0.118 %/yr** against
Hvid-Jensen's nationwide 1.2 per 1000 person-years (0.12 %/yr). Note that this
anchor is matched to the *PPI-treated* arm, not the untreated one, because the
cohort it comes from was largely PPI-treated.

Two of the three AspECT ratios are hit to within 0.02. **The combination arm is
the model's one honest miss**: 0.695 against the trial's 0.59. That gap is
structural, not a tuning failure — the drug effects enter `MUT`
multiplicatively, so the model can only ever predict 0.750 × 0.929 = 0.697, and
it delivers 0.695. The trial found slightly *more* than multiplicative. Either
its combination estimate sits at the favourable end of its interval, or there is
a real acid × COX interaction this model does not contain.

The **UDCA** and **antireflux-surgery** arms are predictions with no trial to
check them against, and they are the two arms that act on the *bile* axis —
which is exactly why they outperform everything the acid-suppression ladder can
do. Treat their size with suspicion; that is where the model is most exposed.

### 5.3 SURF / AIM — endoscopic eradication in confirmed LGD

| Arm | HGD/EAC at 3 y | at 10 y | Prague M at 3 y | at 10 y | Peak stricture |
|---|---|---|---|---|---|
| Surveillance on high-dose PPI | **27.86%** | 35.80% | 4.95 cm | 4.88 cm | — |
| RFA ×3 promptly, on high-dose PPI | **3.77%** | 4.35% | 1.51 cm | 3.64 cm | 0.119 |
| RFA ×3 but **delayed to year 1** | 19.81% | 20.15% | 1.03 cm | 3.47 cm | 0.119 |
| RFA ×3 with **no acid control** | 5.54% | **14.85%** | 1.61 cm | 3.87 cm | 0.119 |
| RFA ×3 + antireflux surgery | 4.80% | **4.80%** | 0.62 cm | **1.48 cm** | 0.119 |
| HGD: EMR + RFA ×4 on high-dose PPI | 6.90% | 7.43% | 1.65 cm | 4.29 cm | 0.130 |
| HGD: surveillance on high-dose PPI | 36.69% | 48.17% | 5.93 cm | 5.83 cm | — |

The surveillance arm lands almost exactly on SURF (27.86% vs the trial's 26.5%
at 3 years). The ablation arm reaches 3.77% against the trial's 1.5% — the right
magnitude and direction (a 7.4-fold reduction where the trial saw 18-fold), but
conservative. The three-session course drives the LGD fraction from 0.50 to
0.0001 within the first year (CE-D), and peak stricture burden is 0.119, inside
the reported 6-12% band.

**Two results worth separating out:**

- **Delay is expensive.** The same three-session course started promptly gives
  3.77% at 3 years; started a year later, 19.81%. Almost the entire difference
  is composite accrued *before* the intervention, in a patient whose baseline
  hazard is ~9 %/yr. (This row exists because the first version of this
  scenario ablated at year 1 and therefore appeared to show that RFA barely
  works.)
- **Durability is a reflux-control result, not an ablation result.** All three
  RFA arms get the same ablation. On high-dose PPI, 10-year risk is 4.35% and
  the segment regrows from 1.51 cm to 3.64 cm. With no acid control, 14.85%.
  With antireflux surgery — the only intervention here that reduces the *bile*
  arm — 4.80% and the segment stays at 1.48 cm. There is no recurrence
  parameter in the model; this falls out of the segment's target length still
  being set by the patient's own CDX2 drive.

### 5.4 CYP2C19 pharmacogenomics (esomeprazole 20 mg OD)

| Phenotype | CL multiplier | Acid blockade | AET | Symptom | Gastrin | 10-y HGD/EAC |
|---|---|---|---|---|---|---|
| UM (\*17/\*17) | 1.80 | 36.8% | 7.49% | 6.28 | 1.64× | 8.77% |
| RM (\*1/\*17) | 1.40 | 42.7% | 6.71% | 6.19 | 1.79× | 7.76% |
| NM (\*1/\*1) | 1.00 | 50.9% | 5.66% | 6.05 | 2.02× | 6.66% |
| IM (\*1/\*2) | 0.55 | 65.3% | 3.90% | 5.81 | 2.43× | 5.47% |
| PM (\*2/\*2) | 0.25 | 80.9% | 2.13% | 5.54 | 2.89× | 5.00% |

A poor metabolizer on **standard-dose** esomeprazole achieves the acid control
(80.9%) that an ultrarapid metabolizer cannot reach even on high-dose BID. The
10-year risk spread between UM and PM on *identical prescriptions* is 8.77% vs
5.00% — a ratio of 0.57, i.e. **larger than the entire high-dose-vs-low-dose
effect the trial measured**. Vonoprazan, being largely CYP2C19-independent,
removes this spread. If the model is right about this, genotype matters more
than dose.

### 5.5 Mechanism dissection — removing the tension terms

| Configuration | `MUT` | Effective bile | Ki-67 | 10-y HGD/EAC | Ratio |
|---|---|---|---|---|---|
| As built | 0.458 | 48.2 | 23.42 | 5.00% | **0.750** |
| No pH-bile penalty (`K_PH_BILE`=0) | 0.368 | 36.5 | 22.94 | 3.09% | 0.463 |
| No gastrin trophic term (`K_GAS_KI`=0) | 0.310 | 48.2 | 15.85 | 2.11% | 0.316 |
| Both off | 0.246 | 36.5 | 15.37 | 1.24% | 0.186 |

With both terms removed, high-dose PPI monotherapy reaches 0.186 — a *quarter*
of the observed effect and better than the real trial's best combination arm.
Neither term was tuned to hit AspECT; each is set from its own mechanistic
literature, and together they land the model on 0.750 against a measured 0.73.
This table is the model's core claim in one place.

### 5.6 Covariates

| Arm | Ki-67 | `MUT` | Prague M at 10 y | 10-y EAC | 10-y HGD/EAC |
|---|---|---|---|---|---|
| Obesity 0.85, maintained | 22.40 | 0.587 | 3.02 cm | 1.642% | 8.34% |
| Obesity 0.85, 50% fat loss at y2 | 20.36 | 0.513 | 2.99 cm | 1.233% | 6.89% |
| Lean (0.15) | 19.04 | 0.467 | 3.01 cm | 0.822% | 5.20% |
| Non-smoker | 20.17 | 0.431 | 3.00 cm | 0.638% | 4.38% |
| SSBE 1 cm | 20.72 | 0.525 | 1.08 cm | 0.486% | 6.00% |
| LSBE 3 cm (reference) | 20.72 | 0.525 | 3.02 cm | 1.181% | 6.66% |
| LSBE 5 cm | 20.72 | 0.525 | 4.95 cm | 1.725% | 7.18% |
| LSBE 8 cm | 20.72 | 0.525 | 7.85 cm | 2.343% | 7.77% |

Segment-length scaling comes out at 1.46× from 3→5 cm and 1.98× from 3→8 cm,
against the ~1.7× per 2 cm reported by Anaparthy. Obese-vs-lean EAC risk is
2.00× (the fitted anchor). Weight loss at year 2 removes **24.9%** of the obese
arm's remaining 10-year EAC risk, acting through both the mechanical reflux
route and the adipokine/PI3K proliferative route — which is why its effect
survives full acid suppression.

### 5.7 The safety counterweight (year 5)

| Regimen | COX inhibition | UGI bleeding | Gastrin |
|---|---|---|---|
| No therapy | 0% | 0.350 %/yr | 1.00× |
| Esomeprazole 20 mg OD | 0% | 0.270 %/yr | 2.02× |
| Esomeprazole 40 mg BID | 0% | 0.221 %/yr | 2.89× |
| Eso 20 mg OD + aspirin | 73.3% | 0.705 %/yr | 2.02× |
| Eso 40 mg BID + aspirin | 73.3% | 0.577 %/yr | 2.89× |

Aspirin roughly doubles upper-GI bleeding risk; co-prescribed high-dose acid
suppression claws back about 18% of that excess. The model therefore reproduces
the actual trade-off of the AspECT combination arm rather than treating aspirin
as free — its 0.695 risk ratio is bought with +0.36 %/yr of bleeding.

## 6. Ten prebuilt scenarios

1. NDBE natural history (no therapy)
2. AspECT low-dose PPI — esomeprazole 20 mg OD
3. AspECT high-dose PPI — esomeprazole 40 mg BID
4. AspECT low-dose PPI + aspirin 300 mg
5. AspECT high-dose PPI + aspirin 300 mg (the best arm)
6. Vonoprazan 20 mg OD (P-CAB, CYP2C19-independent)
7. Confirmed LGD — surveillance on PPI (SURF control)
8. Confirmed LGD — RFA ×3 at year 1 (SURF ablation arm)
9. HGD — EMR of a visible lesion + RFA ×4
10. Antireflux surgery at year 1

Plus a CYP2C19 phenotype sweep (`BE_cyp2c19_sweep()`) and a weight-loss
comparison (`BE_weightloss()`).

## 7. Shiny dashboard — 9 tabs

Patient profile (with value boxes for AET, Prague M, current hazard and
cumulative risk) · Drug PK & acid PD (including a 30-day zoom showing pump
occupancy accumulating over ~5 days) · Reflux & injury axis (the tab that shows
bile *not* falling under acid suppression) · Metaplasia & growth · Clonal
evolution · Clinical endpoints · Scenario comparison (with a ratio-vs-low-dose-PPI
column, i.e. a simulated hazard ratio table) · Ablation & recurrence (interactive
session count, interval, and post-ablation reflux strategy) · PGx & safety.

```r
# from this directory
shiny::runApp("be_shiny_app.R")
```

## 8. Running the model

```r
library(mrgsolve); library(dplyr); library(ggplot2); library(purrr)
mod <- mread_cache("be_mrgsolve_model.R")

# AspECT high-dose arm, 10 years
ppi_high <- ev(amt = 2.0, cmt = "PPI_DEPOT", ii = 0.5, addl = 2*3650 - 1)
asa      <- ev(amt = 1.0, cmt = "ASA_DEPOT", ii = 1,   addl = 3650 - 1)
out <- mod %>% mrgsim(events = ppi_high + asa, end = 3650, delta = 5)
plot(out, CumInc_comp_pct + AET_pct + Mutational_drive + Prague_M_cm ~ time)

# all ten scenarios (helper defined in the model file's documentation block)
df <- BE_simulate_scenarios(3650)
```

Regenerate the map:

```bash
dot -Tsvg be_qsp_model.dot -o be_qsp_model.svg
dot -Tpng -Gdpi=150 be_qsp_model.dot -o be_qsp_model.png
```

## 9. What would falsify this model

| Claim | Falsified by |
|---|---|
| Effective bile toxicity rises under acid suppression | Showing bile-driven CDX2/ROS signaling *falls* at neutral pH |
| Hypergastrinemia is trophic to Barrett epithelium | CCK2R blockade not changing Barrett proliferation on PPI |
| High-dose PPI gives ≈0.74, not ≈0.1 | An RCT where PPI monotherapy abolishes progression |
| Aspirin's benefit is additive because it hits a different arm | A null aspirin effect *on top of* high-dose PPI |
| Post-ablation recurrence is a reflux-control phenomenon | Recurrence rates independent of post-ablation acid/bile exposure — the model predicts 4.35% vs 14.85% 10-year risk for the *same* ablation with vs without acid control |
| Genotype matters more than dose (UM 8.77% vs PM 5.00% on identical prescriptions) | A trial showing CYP2C19-guided dosing does not change progression |
| Delaying ablation in confirmed LGD is costly (3.77% → 19.81% at 3 y for a one-year delay) | A trial showing deferred ablation matches prompt ablation |
| Risk scales ~linearly with segment length | A length-independent progression rate |

## 10. Known limitations

- The clone "fractions" are deterministic field prevalences, not stochastic
  clones; the model cannot represent the variance that makes individual
  progression unpredictable, only the mean field. Curtius 2016 is the right
  companion for stochastic dynamics.
- **`FP53` is a fraction of the *dysplastic* clone pool, not of the whole
  segment.** It settles around 0.52-0.73 in the simulations, which is right for
  TP53 prevalence *within dysplastic clones* (~70% of HGD, >90% of EAC) but must
  **not** be read as aberrant-p53 IHC prevalence in non-dysplastic BE, which is
  5-15%. The hazard it feeds (`K_HGD × MUT × FP53 × FLGD`) only depends on the
  product, so the calibration is unaffected — but the state is easy to
  misinterpret on its own.
- The AspECT combination arm is under-predicted (0.695 vs 0.59) because the drug
  effects are strictly multiplicative in `MUT`. If the trial's super-additivity
  is real, the model is missing an acid × COX interaction.
- The untreated arm (21.9% composite over 10 years) is an extrapolation with no
  randomized comparator — every AspECT patient received a PPI. It implies a PPI
  hazard ratio of ~0.30 vs nothing, in the same direction as the confounded
  observational estimate of ~0.4 (Kastelein 2013) but stronger.
- Sampling error is drawn in the map but not in the ODEs: the model reports
  *true* dysplasia state, whereas a clinician sees a biopsy result. A detection
  layer (Seattle protocol miss rate, WATS-3D uplift) would be the natural next
  extension and would change the apparent LGD/HGD trajectories substantially.
- Drug exposure is in dose-proportional arbitrary units, not mg/L; the PK is
  structurally right (absorption, elimination, CYP2C19 scaling, irreversible
  target binding) but not a validated popPK model.
- EAC treatment (chemoradiation, esophagectomy, adjuvant nivolumab per
  CheckMate 577) is deliberately out of scope — the model stops at the
  transition to invasive cancer, which is where a *prevention* model should
  stop.
- No competing mortality. Over a 10-year horizon in a cohort with a mean age of
  58 this matters: the real AspECT composite included all-cause death, which the
  model's composite does not.

---

*Part of the QSP Disease Model Library. Educational and research use only — not
validated for clinical decision-making or regulatory submission.*
