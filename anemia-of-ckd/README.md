# Anemia of Chronic Kidney Disease (ACKD) — QSP Model

> Integrated Quantitative Systems Pharmacology model of renal anemia, linking
> **nephron loss → REP-cell attrition → blunted HIF-2α/EPO output**, the
> **hepcidin–ferroportin–erythroferrone iron-restriction axis** (including the
> loss of *renal* hepcidin clearance that is unique to CKD), and
> **inflammatory attenuation of EPOR–JAK2–STAT5 signalling**, to the full modern
> pharmacology stack — ESAs, HIF prolyl-hydroxylase inhibitors, oral and IV iron,
> anti-IL-6 adjuncts and transfusion — and to the endpoints the trials actually
> measured.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`ackd_qsp_model.dot`](ackd_qsp_model.dot) |
| 🖼️ Map (SVG, full detail) | [`ackd_qsp_model.svg`](ackd_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)     | [`ackd_qsp_model.png`](ackd_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`ackd_mrgsolve_model.R`](ackd_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`ackd_shiny_app.R`](ackd_shiny_app.R) |
| 📚 References             | [`ackd_references.md`](ackd_references.md) |

**Scale:** 219 nodes · 15 clusters · 333 edges in the map · 30 ODE compartments ·
15 prebuilt scenarios · 106 references.

---

## 1. Disease in one paragraph

Anemia is nearly universal in advanced CKD, and its central lesion is **not**
simple erythropoietin absence. As the kidney scars, the peritubular
interstitial fibroblast-like cells that make EPO — **REP cells** — are lost and
transdifferentiate into myofibroblasts under TGF-β/PDGF, so the *capacity* to
make EPO falls. On top of that, the shrinking nephron mass does less tubular
reabsorptive work and therefore consumes less oxygen, so cortical tissue pO₂ is
paradoxically **near-normal** and the hypoxic drive that should scream at the
remaining REP cells is **blunted**. The result is the signature laboratory
finding of renal anemia: an EPO level that is numerically normal (~10–20
mIU/mL) but grossly *inappropriate* for a haemoglobin of 9 g/dL, where an
intact kidney would produce hundreds. Layered on this are three amplifiers.
First, **functional iron deficiency**: hepcidin is driven up both by IL-6/STAT3
inflammation *and* by the loss of its own renal clearance, so ferroportin is
degraded, the ~11 mg/day of iron recycled from senescent red cells is trapped in
macrophages, and duodenal absorption is capped — stores can look full while the
marrow starves. Second, **shortened red-cell lifespan** (~115 days falling to
70–90 in ESKD) from uremic eryptosis plus dialysis-circuit and GI blood loss,
which raises the obligate production demand. Third, **direct cytokine
suppression** of BFU-E/CFU-E and of EPOR–JAK2–STAT5 signal transduction, which
is what "ESA hyporesponsiveness" actually means at the molecular level.
Therapeutically this is a field with a hard-won lesson: pushing haemoglobin
toward normal with escalating ESA doses produced **no benefit and excess
thrombosis** in the Normal-Hematocrit, CHOIR, CREATE and TREAT trials, so
guidelines now target a deliberately sub-normal 10–11.5 g/dL. The HIF-PHIs
(roxadustat, daprodustat, vadadustat) are interesting precisely because they act
on *all three* amplifiers at once — raising EPO only modestly, suppressing
hepcidin, and raising transferrin — instead of flooding the marrow with
supra-physiological EPO.

## 2. Mechanistic clusters (15 in the DOT map, 219 nodes)

1. **CKD substrate** — primary insult, nephron loss, hyperfiltration,
   tubulointerstitial fibrosis, peritubular capillary rarefaction, REP cells and
   their myofibroblast transdifferentiation, the cortical **O₂ paradox**, uremic
   milieu, dialysis modality, residual renal function.
2. **Oxygen sensing (PHD–VHL–HIF)** — O₂/2-OG/Fe(II)/ascorbate co-substrates,
   succinate–fumarate competition, PHD1/2/3, FIH-1, HIF-1α and HIF-2α,
   Pro402/Pro564 hydroxylation, pVHL–E3 ligase, proteasome, ARNT, p300/CBP, HRE,
   the HIF transcriptional program. *This cluster is the HIF-PHI drug target.*
3. **Erythropoietin** — EPO gene 3′ HRE enhancer, renal (~90%) vs hepatic
   (HIF-2α-inducible) synthesis, GATA-2/NF-κB repression, sialylation and
   half-life, EPOR homodimer, JAK2, STAT5, PI3K/AKT, MAPK/ERK, BCL-xL,
   GATA-1/KLF1, SOCS3 feedback, receptor-mediated clearance.
4. **Erythropoiesis cascade** — HSC → CMP/MEP → BFU-E → CFU-E →
   proerythroblast → basophilic/polychromatic/orthochromatic erythroblast →
   enucleation → marrow and blood reticulocytes → RBC mass; ineffective
   erythropoiesis, marrow fibrosis, haem/globin assembly, RBC lifespan,
   eryptosis, blood loss, splenic clearance.
5. **Systemic iron traffic** — dietary iron, DCYTB, DMT1, enterocyte labile
   pool, ferroportin, hephaestin/ceruloplasmin, transferrin/TIBC/TSAT, TfR1,
   STEAP3, mitochondrial iron, RE macrophage haem recycling (HMOX1), hepatocyte
   stores, ferritin, **functional vs absolute** iron deficiency, NTBI, overload.
6. **Hepcidin–ferroportin hub** — BMP6/BMP2, hemojuvelin, ALK2/ALK3–BMPR2,
   SMAD1/5/8, SMAD4, HAMP BMP-RE, HFE/TfR2 iron sensing, TMPRSS6, IL-6→STAT3,
   **renal hepcidin clearance**, ferroportin ubiquitination, erythroferrone and
   its BMP trap, macrophage/enterocyte iron sequestration.
7. **Inflammation & ESA hyporesponsiveness** — uremic toxins, membrane
   bioincompatibility, gut dysbiosis, access infection, IL-6/TNF-α/IFN-γ/IL-1β,
   hs-CRP, NF-κB, oxidative stress, EPOR-signal attenuation, progenitor
   suppression, secondary hyperparathyroidism, B12/folate/carnitine deficits.
8. **ESA pharmacology** — epoetin, darbepoetin, methoxy-PEG-epoetin beta; SC
   depot and flip-flop absorption, bioavailability, target-mediated disposition,
   half-life spread (8 h → 130 h), EPOR occupancy, protocolised escalation,
   ERFE induction, and the **iron-demand surge** that unmasks functional ID.
9. **HIF-PHI pharmacology** — roxadustat/daprodustat/vadadustat/molidustat/
   enarodustat; phosphate-binder chelation of absorption, albumin binding,
   CYP2C8/UGT1A9/OATP1B1, 2-OG-site PHD inhibition, HIF stabilisation, the
   **coordinated program** (EPO↑ + DMT1/DCYTB↑ + transferrin↑ + hepcidin↓),
   hepatic EPO re-activation, near-physiological EPO peaks, LDL-C lowering,
   off-target VEGF/glycolysis/EMT.
10. **Iron repletion therapy** — oral ferrous salts, ferric citrate,
    sucrosomial/ferric maltol; IV ferric carboxymaltose, iron sucrose,
    ferumoxytol, ferric derisomaltose, dialysate ferric pyrophosphate citrate;
    the hepcidin-gated **oral absorption ceiling**, GI intolerance, RE colloid
    processing, and the PIVOTAL proactive high-dose strategy.
11. **Emerging / adjunct** — anti-hepcidin, ferroportin inhibitors, anti-BMP6,
    anti-IL-6(R), IV ascorbate, transfusion (and its alloimmunisation and iron
    cost), nutrient repletion, calcimimetics.
12. **Haemoglobin, O₂ delivery & cardiac adaptation** — Hb, Hct/MCV/MCH, plasma
    volume dilution, O₂-carrying capacity, compensatory cardiac output, LVH,
    diastolic dysfunction, myocardial supply–demand mismatch, sympathetic/RAAS
    activation, the **cardio-renal-anaemia triad**, Hb cycling.
13. **Clinical endpoints** — Hb in the KDIGO band, Hb response, transfusion
    avoidance, FACIT-Fatigue, SF-36/KDQOL, exercise capacity, hospitalisation,
    MACE/MACE+, eGFR slope/ESKD, mortality, monthly ESA dose and **ERI**.
14. **Safety pathways** — Hb overshoot, viscosity, platelet activation,
    endothelial dysfunction, ESA hypertension, thrombosis (access/VTE/stroke),
    tumour progression, anti-EPO PRCA; HIF-PHI thrombosis/pulmonary-hypertension/
    retinal/cyst-growth concerns; FCM-induced FGF23 → hypophosphatemia, IV iron
    hypersensitivity, labile-iron oxidative injury and infection risk.
15. **Monitoring & titration policy** — Hb q2–4 weeks, TSAT/ferritin q1–3
    months, reticulocytes/RET-He, hepcidin, soluble TfR/ZPP/%HYPO, hs-CRP,
    phosphate, blood pressure, the KDIGO ±25% titration algorithm, and the
    formal hyporesponse workup.

## 3. The mrgsolve model (30 ODEs)

| Block | Compartments |
|---|---|
| Drug PK (6) | `ESA_SC`, `ESA_CENT`, `PHI_GUT`, `PHI_CENT`, `FE_GUT`, `FE_IV` |
| Kidney → EPO (3) | `REP`, `HIF`, `EPO` |
| Erythropoiesis (5) | `PROG`, `EB1`, `EB2`, `RET`, `RBC` |
| Iron (3) | `FE_TF`, `FE_STORE`, `TIBC` |
| Iron regulation (2) | `HEPC`, `ERFE` |
| Inflammation (2) | `IL6`, `CRP` |
| Mineral / lipid / haemodynamics (4) | `FGF23`, `PHOS`, `LDL`, `MAP` |
| Safety & outcomes (5) | `NTBI`, `THR`, `FATIGUE`, `LVMI`, `EGFR` |

### Two design decisions worth knowing about

**The baseline is self-calibrating.** Rather than hard-coding 30 initial
conditions, `[MAIN]` solves the stationary state analytically from the patient
descriptors (`EGFR0`, `BASE_HB`, `FIB`, `BLUNT`, `IL6_DRIVE`, `URE`, `TSAT0`,
`FERRITIN0`, `HEPC_BASE`, …) and then back-calculates the production constants
(`KIN_PROG`, `KIN_HEPC`, `KIN_IL6`, `KIN_TIBC`, `KDEP`, `LOSS_FE`, …) that hold
it. Build any virtual patient and the model still starts in balance — verified
numerically: every derivative is 0 at *t* = 0 except `EGFR` (which is *supposed*
to decline) and `THR` (a cumulative hazard integral).

**Iron balance is closed exactly, on purpose.** At steady state reticulocyte
delivery equals RBC destruction, so recycling equals erythron utilisation
(~11.2 mg/day here); plasma iron therefore balances *iff* obligate loss equals
baseline absorption, and the model sets `LOSS_FE = abs0` to enforce that.
Consequence: any net iron deficit you observe in a scenario comes from a
mechanism the model actually represents — dialysis/GI loss, hepcidin trapping,
or a drug-driven demand surge — and never from bookkeeping drift.

### Emergent behaviours (not fitted, they fall out of the structure)

- **Functional iron deficiency on ESA.** ERFE rises, hepcidin falls ~35–45%,
  and ferritin *drains* from 147 → 75–106 ng/mL while TSAT stays roughly flat.
  Stores pay for the new red cells; that store-drain-at-constant-TSAT pattern is
  the model's central iron-kinetic signature.
- **Oral iron futility.** Absorption is hepcidin-gated, and hepcidin rises with
  the very stores oral iron creates, so ferritin climbs to ~429 ng/mL while Hb
  moves +0.3 g/dL.
- **The HIF-PHI fingerprint.** EPO peaks only ~3–7× baseline (vs ~45–150× for a
  Q2W/Q4W ESA bolus), hepcidin drops hardest of any arm, TIBC rises, and LDL-C
  falls ~8–14%.
- **The overshoot penalty.** Cumulative thrombotic hazard scales ~4.3× from the
  maintenance arm to the Hb-normalisation arm, driven by viscosity, MAP and ESA
  exposure terms.
- **Hb cycling.** Stopping an ESA lets Hb fall on the RBC-lifespan time constant
  (weeks, not days) while hepcidin *rebounds above baseline* as ERFE collapses.

## 4. Validation

The model was checked two independent ways. First, the entire ODE system was
re-implemented from scratch in pure Python with an RK4 integrator and the two
implementations were compared — the epoetin-maintenance arm lands on **Hb 11.70
g/dL at day 365 in both**. Second, the spec was compiled and run under real
`mrgsolve` 1.x (R 4.3.3), which is where the table below comes from. Three
defects were found and fixed during this process, and they are worth recording
because each was a genuine modelling error rather than a typo:

| Defect found | Fix |
|---|---|
| TSAT reached 119% — physically impossible | Added a TSAT-dependent overflow term (`KOVER_DEP`) that shunts surplus iron into stores, bounding TSAT below 100% (max observed now 65%) |
| IV iron *worsened* anemia: ferritin 890 ng/mL with TSAT 6% | Hepcidin's store drive was a runaway power law and mobilisation capacity was too small, so high hepcidin became a one-way ratchet burying iron in stores. Added hepcidin-independent floors (`FPN_FLOOR`, `MOB_FLOOR`), capped the store drive, and raised `MOB_MAX` so loaded stores can actually rescue the marrow |
| Inflammation didn't blunt the drug response at all | Because `KIN_PROG` is back-calculated per patient, inflammation was absorbed into the baseline. Added `KIL6_EPOR`, an explicit cytokine attenuation of EPOR→JAK2→STAT5 transduction that scales the **drug** signal |
| `capture LVMI = LVMI` collided with the compartment name | Renamed captures to `LVMI_g_m2` / `eGFR_ml_min` (caught only by real mrgsolve) |
| `init()` was silently overridden by `[MAIN]` | Added `INIT_FROM_PARAM`; with it set to 0 a resumed run now reproduces a continuous run **exactly** — which is what makes the closed-loop titrator valid |

### All 15 scenarios at day 365 (real mrgsolve output)

`band` = days with Hb in 10–11.5 g/dL · `>13` = days above 13 g/dL ·
`Haz` = cumulative thrombotic hazard (a **ranking index**, not an event rate) ·
`EPOmax` = peak total erythropoietic signal (mIU/mL equivalent).

| Scenario | Hb | ΔHb | band | >13 | TSAT | Ferr | Hepc | EPOmax | LDL | MAP | P min | Haz | ΔFACIT | ΔLVMI |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1. Untreated natural history (CKD G4) | 9.41 | -0.09 | 0 | 0 | 23.6 | 148 | 61.7 | 9 | 100 | 95.0 | 4.50 | 1.00 | -0.3 | +0.5 |
| 2. Epoetin alfa 75 U/kg SC TIW (maintenance) | 11.70 | +2.20 | 170 | 0 | 25.0 | 106 | 50.3 | 49 | 104 | 98.9 | 4.50 | 1.28 | +6.1 | -14.2 |
| 3. Epoetin alfa 150 U/kg SC TIW (correction) | 13.18 | +3.68 | 60 | 108 | 26.1 | 91 | 37.3 | 89 | 104 | 102.3 | 4.50 | 1.70 | +9.3 | -13.0 |
| 4. Darbepoetin alfa Q2W | 13.94 | +4.44 | 47 | 213 | 26.3 | 83 | 32.2 | 404 | 104 | 103.9 | 4.50 | 2.43 | +10.6 | -11.8 |
| 5. Methoxy-PEG-epoetin beta Q4W | 14.78 | +5.28 | 37 | 252 | 26.4 | 75 | 26.8 | 1403 | 104 | 107.1 | 4.50 | 3.68 | +11.8 | -10.0 |
| 6. Epoetin + proactive IV iron (HD, PIVOTAL) | 13.44 | +3.94 | 34 | 225 | 52.3 | 801 | 97.3 | 85 | 104 | 102.4 | 4.50 | 2.11 | +9.5 | -12.9 |
| 7. Roxadustat TIW (dialysis) | 10.95 | +1.45 | 339 | 0 | 24.6 | 145 | 26.2 | 65 | 90 | 95.0 | 4.50 | 1.20 | +4.2 | -13.3 |
| 8. Daprodustat daily (non-dialysis) | 10.42 | +0.92 | 332 | 0 | 25.2 | 173 | 48.6 | 33 | 96 | 95.0 | 4.50 | 1.14 | +2.8 | -8.6 |
| 9. Oral iron alone (hepcidin-capped) | 9.79 | +0.29 | 88 | 0 | 24.4 | 429 | 121.5 | 9 | 101 | 95.0 | 4.50 | 1.00 | +1.1 | -4.4 |
| 10. Inflammatory hyporesponder, matched ESA dose | 12.74 | +3.24 | 71 | 0 | 25.7 | 96 | 40.1 | 85 | 104 | 101.6 | 4.50 | 1.51 | +8.5 | -13.2 |
| 11. Inflammatory hyporesponder, ESA escalated 3× | 14.53 | +5.03 | 38 | 245 | 26.1 | 78 | 27.4 | 246 | 104 | 105.9 | 4.50 | 3.25 | +11.5 | -10.9 |
| 12. Inflammatory hyporesponder, roxadustat | 10.35 | +0.85 | 323 | 0 | 26.8 | 202 | 38.4 | 40 | 92 | 95.0 | 4.50 | 1.17 | +2.5 | -7.7 |
| 13. Hb normalisation overshoot (safety) | 15.66 | +6.16 | 32 | 278 | 53.6 | 654 | 78.8 | 409 | 104 | 107.9 | 3.09 | 5.54 | +12.8 | -9.6 |
| 14. High-dose FCM → FGF23 → hypophosphatemia | 9.20 | -0.30 | 0 | 0 | 21.7 | 368 | 92.4 | 11 | 99 | 95.0 | **1.49** | 1.02 | -0.9 | +2.0 |
| 15. ESA interruption at day 120 (Hb cycling) | 9.49 | -0.01 | 145 | 0 | 23.7 | 140 | 61.1 | 89 | 100 | 95.0 | 4.50 | 1.14 | +0.1 | -6.3 |

**How to read this table.** The three arms that keep a patient inside the KDIGO
band for most of the year are the two HIF-PHI arms (339 and 332 days) and the
ESA *maintenance* dose (170 days). Every arm that pushes harder buys ΔHb and
ΔFACIT at the cost of days above 13 g/dL and a rising hazard index — scenario 13
is the model's rendering of the Normal-Hematocrit/CHOIR/TREAT result, and
scenario 11 shows the same trade-off arriving through dose escalation in a
hyporesponder. Scenario 14 is the only arm with a phosphate signal (nadir 1.49
mg/dL), and setting `FCM_FLAG = 0` abolishes it.

### Closed-loop titration

With the `INIT_FROM_PARAM` fix the Shiny app's KDIGO-style titrator (±25% steps
every 28 days, hold above the ceiling) behaves the way a clinic does: Hb climbs
to 11.77 g/dL by block 6, the dose is stepped down 0.625 → 0.198 au, and Hb
settles at **10.40 g/dL — inside the band, with zero days above 13**. The
open-loop scenarios overshoot precisely because they lack this feedback.

## 5. Using the model

```r
library(mrgsolve); library(dplyr)

# The spec is named *_mrgsolve_model.R, so load it by filename:
mod <- mread_file("ackd_mrgsolve_model.R", project = ".")

# Epoetin alfa 75 U/kg SC three times weekly for a year
tiw <- 7 / 3
out <- mod |>
  param(ESA_WK_UKG = 225) |>
  mrgsim(events = ev(amt = 0.5, cmt = "ESA_SC", ii = tiw,
                     addl = floor(365 / tiw)), end = 365) |>
  as_tibble()

# Switch molecules / classes by parameter, not by editing the model:
#   darbepoetin  : KA_ESA=0.25, F_ESA=0.45, KEL_ESA=0.65
#   roxadustat   : KA_PHI=3.0,  F_PHI=0.85, KEL_PHI=0.55  (dose into PHI_GUT)
#   iron sucrose : FCM_FLAG=0                              (no FGF23 signal)
```

Dose units: `ESA_SC` 1.0 au ≈ epoetin alfa 150 U/kg · `PHI_GUT` 1.0 au ≈
roxadustat 100 mg · `FE_GUT`/`FE_IV` in real mg elemental iron · `RBC` 0.33 ≈
one unit of packed cells.

Then launch the dashboard (8 tabs: patient profile · drug PK · erythropoiesis PD
· iron & hepcidin · clinical endpoints · scenario comparison · safety ·
titration and docs):

```r
shiny::runApp("ackd_shiny_app.R")
```

## 6. Limitations — read before drawing conclusions

- **The baseline is pinned.** `[MAIN]` back-calculates `KIN_PROG` so *every*
  virtual patient starts at `BASE_HB`, including the inflamed one. Raising
  `IL6_DRIVE` therefore does **not** lower the presenting Hb; inflammation shows
  up as a smaller Hb gain per unit ESA and a tighter iron gate. Always compare
  hyporesponders at **matched dose** (scenario 10 vs 3), never by their
  proportional reserve.
- Drug exposures are in dose-proportional `au`, first-order, with no
  target-mediated (EPOR) disposition — so the model cannot speak to
  concentration-based therapeutic monitoring.
- One lumped `HIF` activity stands in for HIF-1α and HIF-2α across kidney, liver
  and gut, so isoform-selective agents are out of scope.
- Hepcidin, ERFE and IL-6 are single well-mixed pools: no hepatic/marrow
  compartmentalisation, no diurnal rhythm.
- `FMAX_ABS = 0.28` comes from iron-deficient healthy volunteers; real CKD
  patients absorb less, so scenario 9 probably **overstates** oral-iron store
  loading.
- The thrombotic hazard **ranks** arms. It is not calibrated to events/100
  patient-years and must not be read as one.
- Long-horizon absolute Hb in the open-loop arms should not be taken literally:
  no clinician would hold a correction dose for 12 months. Use the titrator for
  realistic trajectories.

---

**Disclaimer.** Educational and research QSP model, built from public literature
and trial data. It has not been independently validated or certified, and must
not be used for clinical decisions, prescribing, or regulatory submissions.
Parameters are illustrative literature-anchored approximations; the references
support the *structure* of each mechanism rather than certifying the numerical
value attached to it. Fitting and qualification against real patient data would
be required for any applied use.
