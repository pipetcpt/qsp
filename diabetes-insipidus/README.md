# Diabetes Insipidus — Arginine Vasopressin Disorders (AVP-D / AVP-R)

> Quantitative Systems Pharmacology model of central (AVP-deficiency, AVP-D)
> and nephrogenic (AVP-resistance, AVP-R) diabetes insipidus.
> Includes hypothalamic AVP secretion, V2-receptor / cAMP / AQP2 trafficking,
> renal water handling, thirst-driven drinking, and PK/PD for desmopressin
> (SC · intranasal · oral · sublingual lyophilisate), thiazide, amiloride,
> indomethacin, tolvaptan (comparator) and lithium (NDI inducer).

| Artifact | File |
|----------|------|
| 🗺️ Mechanistic map (Graphviz)   | [`di_qsp_model.dot`](di_qsp_model.dot) → [SVG](di_qsp_model.svg) · [PNG](di_qsp_model.png) |
| ⚙️ mrgsolve ODE model            | [`di_mrgsolve_model.R`](di_mrgsolve_model.R) |
| 📊 Shiny dashboard               | [`di_shiny_app.R`](di_shiny_app.R) |
| 📚 References                    | [`di_references.md`](di_references.md) |

## Disease overview

The 2022 international working group renamed *central diabetes insipidus* to
**Arginine Vasopressin Deficiency (AVP-D)** and *nephrogenic diabetes
insipidus* to **Arginine Vasopressin Resistance (AVP-R)** to prevent confusion
with diabetes mellitus. Both share the phenotype of hypotonic polyuria
(>50 mL/kg/day, urine osmolality <300 mOsm/kg) and compensatory polydipsia.

* **AVP-D**: loss of magnocellular AVP neurons in SON/PVN (autoimmune,
  craniopharyngioma/germinoma, transsphenoidal surgery, traumatic brain
  injury, infiltrative disease, *AVP* gene mutations) → low or absent
  copeptin.
* **AVP-R**: collecting-duct unresponsiveness — *AVPR2* (X-linked, V2R) or
  *AQP2* mutations, plus acquired forms (lithium >10 y use, hypercalcaemia,
  hypokalaemia, post-obstructive).
* **Gestational DI**: placental vasopressinase degrades native AVP but
  spares desmopressin.
* **Primary polydipsia**: behavioural excess water intake without AVP
  defect; differentiated by hypertonic-saline or arginine-stimulated
  copeptin.

## Mechanistic map (140+ nodes, 12 clusters)

Cluster layout: ① Etiology / genetics, ② Hypothalamic-neurohypophyseal axis,
③ Osmolality & volume sensing, ④ V2R / cAMP / PKA signalling, ⑤ AQP2
trafficking and long-term regulation, ⑥ Renal tubular water handling,
⑦ Thirst behaviour & volume status, ⑧ Drug PK/PD (desmopressin, thiazide,
amiloride, indomethacin, tolvaptan, lithium), ⑨ Co-existing pituitary
dysfunction, ⑩ Diagnostic biomarkers (copeptin, water deprivation), ⑪
Complications, ⑫ Special populations.

Render with:

```bash
dot -Tsvg di_qsp_model.dot -o di_qsp_model.svg
dot -Tpng -Gdpi=150 di_qsp_model.dot -o di_qsp_model.png
```

## mrgsolve model (20 ODE compartments)

| Compartment | Meaning |
|-------------|---------|
| `DEPOT_SC/IN/PO/SL` | Desmopressin depots (4 routes) |
| `CENT_DDAVP`, `PERI_DDAVP` | 2-cpt desmopressin disposition |
| `HCTZ_DEPOT/CENT`, `AMI_DEPOT/CENT`, `IND_DEPOT/CENT`, `TOL_DEPOT/CENT`, `LI_DEPOT/CENT` | Concomitant agents |
| `AVP_E` | Endogenous plasma AVP (pmol/L) — driven by osmolality, scaled by AVP-D severity, accelerated degradation when gestational vasopressinase is on |
| `AQP2_M`, `AQP2_A` | AQP2 mRNA pool and apical-membrane fraction (short-term shuttling + long-term transcription) |
| `TBW`, `NA_BODY` | Total body water (L) and exchangeable sodium (mmol) — yield plasma Na+, osmolality |
| `NDI_LI` | Lithium-induced NDI severity, slow accrual ↔ slow reversal |
| `CUM_URINE`, `CUM_HAZ` | Cumulative urine output, hyponatraemia hazard integral |

Key calibrations (see `di_references.md`):

* AVP-osmolality slope ≈ 0.4 pmol/L per mOsm/kg above 280 (Robertson 1976).
* AVP plasma t½ ≈ 10-30 min (Robertson 1973).
* DDAVP V2 potency ≈ 12× native AVP (Vavra 1968).
* DDAVP bioavailability: SC 100 % · IN ~4 % · PO ~0.16 % · sublingual MELT ~0.25 %
  (Lottermoser 1997, Steiner 2007, Fjellestad-Paulsen 1993).
* Maximum urine concentrating osmolality 1200 mOsm/kg; minimum 50 mOsm/kg
  (Knepper 2015).
* Lithium NDI develops slowly when plasma Li+ > 0.6 mmol/L; amiloride blocks
  ENaC-mediated Li+ entry (Batlle 1985 NEJM).
* Thiazide paradoxical antidiuresis in AVP-R: ~45 % reduction of free water
  loss (Earley 1962; Magaldi 2000).

### Pre-built scenarios (`run_scenario()`)

| Key | Description |
|-----|-------------|
| `untreated_CDI` | Complete AVP-D, no therapy |
| `DDAVP_SC_2ug_BID` | Desmopressin 2 µg SC q12h |
| `DDAVP_IN_10ug_BID` | Desmopressin 10 µg intranasal q12h |
| `DDAVP_PO_200ug_TID` | Desmopressin 200 µg oral q8h |
| `DDAVP_SL_120ug_TID` | Sublingual MELT 120 µg q8h |
| `NDI_lithium_HCTZ` | Severe AVP-R on hydrochlorothiazide 25 mg/d |
| `NDI_lithium_amiloride` | Lithium loading + amiloride 10 mg/d rescue |
| `NDI_indomethacin` | AVP-R + indomethacin 50 mg q8h |
| `tolvaptan_SIADH_comparator` | Tolvaptan 15 mg/d (V2 antagonist control) |
| `primary_polydipsia` | Behavioural +6 L/d intake, intact AVP axis |
| `gestational_DDAVP` | Pregnancy with vasopressinase + DDAVP IN |
| `pediatric_DDAVP_SC` | 20 kg child, DDAVP 0.3 µg SC q12h |

## Shiny app (8 tabs)

1. **Patient profile** — demographics & disease phenotype assignment
2. **Disease severity** — AVP-osmolality response curve and predicted basal labs
3. **Desmopressin PK** — plasma trajectory across the chosen route
4. **AQP2 / V2R PD** — V2 occupancy, apical AQP2, urinary osmolality
5. **Electrolytes & osmolality** — plasma Na+, osmolality, TBW
6. **Urine & thirst endpoints** — urine flow, cumulative output, thirst score
7. **Scenario comparison** — current vs untreated CDI reference
8. **Biomarkers & safety** — copeptin (Christ-Crain cut-off 4.9 pmol/L) and
   hyponatraemia hazard integral

Launch with:

```r
shiny::runApp("diabetes-insipidus/di_shiny_app.R")
```

## Highlighted endpoints

* **Urine output (L/d)** — primary efficacy endpoint; target ≤2.5 L/d on therapy.
* **Plasma Na+ (135-145 mmol/L)** — safety; hyponatraemia hazard accumulates
  below 130 mmol/L.
* **Urinary osmolality (mOsm/kg)** — concentration target ≥300 indicates
  V2/AQP2 response.
* **Stimulated copeptin (pmol/L)** — diagnostic cut-off 4.9 separates AVP-D
  from primary polydipsia (Christ-Crain 2019, Refardt 2023 NEJM).
* **V2 occupancy** — derived from total V2 agonist (AVP + DDAVP-equivalent).

## Caveats

* The model is *educational* — not validated for individual dosing decisions.
* AVP plasma kinetics use a single-compartment lumped representation;
  finer 2-compartment fitting requires individualized data.
* Long-term lithium NDI severity is modelled phenomenologically; recovery
  may take months in real patients.
* Thiazide and indomethacin effects are simplified empirical multipliers on
  the free-water excretion rate.
