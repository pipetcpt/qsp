# SIADH (Syndrome of Inappropriate Antidiuretic Hormone Secretion) — QSP Model

> Integrated Quantitative Systems Pharmacology model of SIADH, linking
> non-osmotic/autonomous arginine-vasopressin (AVP) drive to renal V2-receptor
> and aquaporin-2 (AQP2) functional expression, free-water retention and
> euvolemic dilutional hyponatremia, cerebral osmotic adaptation, and
> correction-rate-dependent osmotic-demyelination-syndrome (ODS) risk — with
> the modern pharmacology stack (tolvaptan, conivaptan, demeclocycline, oral
> urea, hypertonic saline, fluid restriction) and the clinical endpoints used
> to track and safely manage the disease (serum sodium, serum/urine
> osmolality, urine sodium, neuro-symptom severity, ODS risk).

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`siadh_qsp_model.dot`](siadh_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`siadh_qsp_model.svg`](siadh_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`siadh_qsp_model.png`](siadh_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`siadh_mrgsolve_model.R`](siadh_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`siadh_shiny_app.R`](siadh_shiny_app.R) |
| 📚 References             | [`siadh_references.md`](siadh_references.md) |

---

## 1. Disease in one paragraph

SIADH is the most common cause of hyponatremia in hospitalized patients and is
defined by **non-osmotic or autonomous release of arginine vasopressin (AVP)**
despite low plasma osmolality — a failure of the normal negative-feedback
suppression of AVP secretion. Causes include **ectopic tumoral production**
(classically small-cell lung cancer), **CNS and pulmonary disease**, and
numerous **drugs** (SSRIs, carbamazepine, cyclophosphamide/vincristine, MDMA),
as well as non-osmotic stimuli like pain, nausea, and stress. Persistently
elevated AVP binds renal **V2 receptors** on collecting-duct principal cells,
driving cAMP/PKA-mediated **aquaporin-2 (AQP2)** phosphorylation, apical
trafficking, and gene upregulation — increasing water permeability and
**free-water reabsorption** despite clinical **euvolemia** (RAAS suppressed,
mild natriuresis). The resulting **dilutional hyponatremia** causes brain-cell
swelling; over 24-48+ hours the brain adaptively **extrudes organic osmolytes**
to blunt further swelling, which paradoxically makes **chronic**, adapted
hyponatremia the setting most vulnerable to **osmotic demyelination syndrome
(ODS)** if serum sodium is corrected **too fast** — the reason guideline
correction-rate limits (≤8 mEq/L/24h, ≤18 mEq/L/48h) exist. Therapy targets
specific nodes: **fluid restriction** (first-line), **vaptans** (tolvaptan,
conivaptan — competitive V2 antagonism), **demeclocycline** (post-receptor
induced nephrogenic diabetes insipidus), **oral urea** (osmotic diuresis), and
**hypertonic saline** for acute/severe symptomatic disease, always balanced
against correction-rate safety.

## 2. Mechanistic clusters (10 in the DOT map, 109 nodes)

1. Etiologies & risk factors (SCLC/ectopic tumors, CNS/pulmonary disease, SSRIs, carbamazepine, cyclophosphamide/vincristine, MDMA, exogenous desmopressin, NSAIDs, postoperative stimuli, exercise-associated hyponatremia, low dietary solute)
2. Hypothalamic-neurohypophyseal AVP axis (osmoreceptors, baroreceptor reflex, non-osmotic AVP stimuli, ectopic tumoral production, reset-osmostat variant, Robertson subtype classification)
3. Renal V2-receptor / aquaporin-2 signaling (Gs-adenylate cyclase-cAMP-PKA, AQP2 phosphorylation/trafficking/membrane insertion/gene expression, NSIAD gain-of-function variant)
4. Water balance & osmoregulation (free-water reabsorption, urine/serum osmolality, electrolyte-free water clearance, total body water, dilutional hyponatremia, low-solute limit)
5. Volume status & renal sodium handling (euvolemia, RAAS suppression, ANP/BNP, natriuresis, urine sodium, cerebral-salt-wasting differential)
6. Neurological consequences of hyponatremia (cerebral edema, astrocyte swelling, brain organic-osmolyte efflux, chronic adaptation, seizures, coma, herniation, ODS risk)
7. Clinical manifestations & endpoints (asymptomatic/mild through severe symptoms, acute-vs-chronic classification, falls/fracture risk, composite neuro-symptom score, mortality association)
8. Drug PK/PD — treatment (tolvaptan, conivaptan, demeclocycline, hypertonic saline, fluid restriction, oral urea, loop diuretic + salt tablets, SGLT2 inhibitor, DDAVP overcorrection-prevention clamp)
9. Na-correction-rate safety & ODS prevention (24h/48h correction trackers, guideline limits, overcorrection flag, relowering protocol, high-risk-patient flag)
10. Biomarkers & diagnostic monitoring (copeptin, serum/urine osmolality, urine sodium, FENa/FEUrea, Bartter-Schwartz/Verbalis algorithm, hypouricemia, TSH/cortisol exclusion)

## 3. mrgsolve model (19 ODE compartments)

* **Drug PK / exposure (7 compartments)** — tolvaptan (oral depot + effective
  exposure), conivaptan (IV effective exposure), demeclocycline (oral depot +
  effective exposure), oral urea (effective exposure), and a hypertonic-saline
  3% infusion driver. As with other topical/local surrogates in this library,
  each is a day-resolution "sustained exposure" compartment mapped through an
  Emax function onto its mechanistic target, not a literal multi-compartment
  plasma PK model.
* **Renal V2R/AQP2 and water/sodium balance (8 compartments)** — AQP2
  functional expression, free-water clearance, total body water, serum
  sodium, serum osmolality, urine osmolality, urine sodium, brain
  organic-osmolyte adaptation index.
* **Correction-rate safety & outcome (4 compartments)** — ~24h- and
  ~48h-delayed sodium trackers (first-order lag compartments approximating a
  rolling correction-rate window), composite neuro-symptom severity score,
  and cumulative osmotic-demyelination-syndrome (ODS) risk score.

### 7 prebuilt scenarios

| # | Scenario | Calibration anchor |
|---|---|---|
| 1 | Untreated natural history (chronic, no restriction) | Bartter & Schwartz 1967 pathophysiology |
| 2 | Fluid restriction alone (<1 L/day, first-line) | Verbalis 2013; Spasovski 2014 guidelines |
| 3 | Tolvaptan 15mg PO QD (up-titrate to 60mg) | SALT-1/SALT-2 (Schrier 2006) |
| 4 | Conivaptan IV load + infusion (inpatient) | Ghali 2006; Verbalis 2008 |
| 5 | Demeclocycline 900mg/day (chronic, malignancy-associated) | Forrest 1978; De Troyer 1975 |
| 6 | Hypertonic saline 3% (acute severe, guideline-limited rate) | Sterns 2015; Spasovski 2014 |
| 7 | Overly rapid correction (unmonitored, ODS teaching example) | Sterns 1986; Sterns 1994 (ODS risk) |

## 4. Shiny dashboard (7 tabs)

1. **Patient profile** — severity slider, ectopic/malignancy flag, acute-vs-chronic flag, starting serum Na, adjunct therapy checkboxes, DDAVP-clamp flag.
2. **Drug PK** — fractional drug-effect trajectories for all four agents.
3. **Renal V2R/AQP2 PD** — AQP2 functional expression, free-water clearance, total body water.
4. **Sodium & osmolality** — serum Na, serum/urine osmolality, urine Na, and 24h/48h correction-rate vs guideline-limit lines.
5. **Neuro symptoms/ODS** — brain organic-osmolyte adaptation index, neuro-symptom score, cumulative ODS risk score.
6. **Scenario comparison** — runs all 7 built-in scenarios with the chosen patient profile.
7. **Biomarkers** — endpoint summary table (Day-0, Day-2, Day-end).

## 5. How to run

```bash
# 1) Render the mechanistic map
dot -Tsvg siadh_qsp_model.dot -o siadh_qsp_model.svg
dot -Tpng -Gdpi=150 siadh_qsp_model.dot -o siadh_qsp_model.png
```

```r
# 2) Simulate scenarios in R
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny",
                   "shinydashboard","DT"))
library(mrgsolve)
mod <- mread_cache("siadh_mrgsolve_model.R")
out <- mrgsim(mod, ev(amt=1.0, cmt="TOLVA_GUT", ii=1, addl=13), end=14)
plot(out, c("Serum_Na","Na_correction_rate_24h","ODS_risk_score","Neuro_symptom_score"))

# 3) Launch the dashboard
shiny::runApp("siadh_shiny_app.R")
```

## 6. Key clinical anchors used during calibration

| Endpoint | Comparator | Expected effect |
|---|---|---|
| Serum Na rise | Tolvaptan vs placebo, day 4 (SALT-1/2) | +3-4 mEq/L, sustained to day 30 |
| Serum Na rise | Conivaptan IV vs placebo | Similar/faster onset than oral vaptan (loading dose + infusion) |
| Serum Na rise | Demeclocycline, chronic SIADH | Slower onset (3-5 days) than vaptans; effective when vaptans unsuitable |
| Correction-rate limit | Guideline (Verbalis 2013; Spasovski 2014) | ≤8 mEq/L/24h (≤4-6 if high-risk), ≤18 mEq/L/48h |
| ODS risk | Chronic + rapid correction vs acute + rapid correction | Markedly higher in chronic/adapted brain (Sterns 1986, 1994) |
| Falls/attention deficits | Mild chronic hyponatremia vs normonatremia | Significantly increased even at Na 128-134 (Renneboog 2006) |

## 7. Caveats

* Designed for **research, education, and hypothesis generation** — not
  clinical decision support.
* Drug "exposure" compartments are pragmatic day-resolution surrogates for
  sustained mechanistic exposure, not literal multi-compartment plasma PK.
* The 24h/48h sodium-correction-rate trackers are first-order lag
  approximations of a rolling window, not an exact retrospective maximum —
  adequate for illustrating guideline-threshold dynamics, not for real
  bedside correction-rate calculations.
* The cerebral-adaptation → ODS-risk link is a simplified, qualitative
  representation of a well-established but incompletely quantified clinical
  phenomenon; absolute ODS_risk_score values are illustrative, not validated
  probabilities.
* Cerebral salt wasting (the key **hypovolemic** differential diagnosis,
  requiring salt/volume replacement rather than restriction) is represented
  only as a differential-diagnosis node, not simulated as its own phenotype.

## 8. License

Inherits the repository [LICENSE](../LICENSE).
