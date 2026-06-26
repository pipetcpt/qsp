# Primary Aldosteronism (Conn's Syndrome) — QSP Model

[![Mechanistic Map](pa_qsp_model.png)](pa_qsp_model.svg)

---

## Disease Overview

**Primary Aldosteronism (PA)**, also known as Conn's Syndrome, is the most common form of endocrine hypertension, affecting 5–10% of hypertensive patients worldwide. It is characterised by autonomous overproduction of aldosterone from the adrenal cortex, independent of its normal physiological regulator (angiotensin II), leading to:

- Suppressed plasma renin activity (PRA) and markedly elevated aldosterone-to-renin ratio (ARR)
- Sodium retention and volume expansion → resistant hypertension
- Urinary potassium wasting → hypokalemia (in ~37% of cases)
- Metabolic alkalosis (bicarbonate retention)
- Disproportionate cardiovascular organ damage (cardiac fibrosis, LVH, stroke) beyond that explained by blood pressure alone

### Subtypes
| Subtype | Frequency | Treatment |
|---------|-----------|-----------|
| Aldosterone-Producing Adenoma (APA) | ~35% | Laparoscopic adrenalectomy |
| Bilateral Adrenal Hyperplasia (BAH) | ~60% | MR antagonist (lifelong) |
| Unilateral adrenal hyperplasia | ~3% | Adrenalectomy |
| PA with familial hyperaldosteronism | ~2% | Glucocorticoid-remediable (FH-I) |

---

## Mechanistic Map

The `.dot` / `.svg` / `.png` mechanistic diagram covers 10 subgraph clusters:

1. **RAAS Cascade** — Renin → AngI → AngII → Aldosterone (ACE, ACEi nodes)
2. **Adrenal Cortex** — KCNJ5/CACNA1D/ATP1A1 somatic mutations → Ca²⁺ influx → CYP11B2 upregulation → autonomous aldosterone synthesis
3. **Renal Collecting Duct** — MR activation → SGK1 → Nedd4-2 phosphorylation → ENaC surface expression → Na⁺ reabsorption / K⁺ and H⁺ secretion (ROMK, H⁺-ATPase)
4. **Cardiovascular System** — Volume expansion → MAP, TPR; aldosterone direct effect on vascular MR and fibroblast MR
5. **Ion Homeostasis** — Na⁺, K⁺, HCO₃⁻, Cl⁻ compartments; mineralocorticoid escape (ANP)
6. **Target Organ Damage** — Cardiac fibrosis (collagen I/III), LV mass index (LVH), renal fibrosis
7. **Drug PK** — Spironolactone (→ canrenone active metabolite), eplerenone, finerenone, ACEi (ramipril), CCB (amlodipine), adrenalectomy effect
8. **Drug PD** — MR competitive antagonism; ACEi suppression of AngII; CCB on TPR; finerenone extra anti-fibrotic
9. **Diagnostics** — ARR calculation, PAC (ng/dL), PRA (ng/mL/h), AVS lateralization, CT adrenal
10. **Clinical Endpoints** — BP normalisation, hypokalemia resolution, ARR normalisation, LVMi regression, CV event risk

---

## mrgsolve ODE Model

**File:** `pa_mrgsolve_model.R`

### Compartments (23 ODEs)

| Group | Compartments |
|-------|-------------|
| Drug PK | C_spiro, C_canrenone, C_eple, C_fine, C_acei, C_ccb |
| RAAS | Renin_c, AngII_c, Aldo_c |
| Renal/Ions | ENaC_act, Na_c, K_c, HCO3_c, Vol_c |
| Cardiovascular | MAP_c, TPR_c |
| Renal function | GFR_c |
| Target organ damage | CardFib, LVMi_c |
| Adrenal | APA_act, CYP11B2_c |
| Biomarkers | ARR_c, HOMA_proxy |

### Treatment Scenarios

| Scenario | Description |
|----------|-------------|
| 1. Untreated APA | APA_severity=1.0, no treatment, 2-year progression |
| 2. Adrenalectomy | Surgical cure (efficacy 90%), RAAS normalisation |
| 3. Spironolactone 100 mg | BAH phenotype, steroidal MRA |
| 4. Eplerenone 100 mg | Selective MRA, less anti-androgenic |
| 5. Finerenone 20 mg | Non-steroidal MRA, superior cardiac anti-fibrotic |
| 6. Spiro + Amlodipine | Combination for severe hypertension |
| 7. Normal Control | Baseline healthy RAAS |
| 8. ACEi (Ramipril 10 mg) | Limited efficacy in PA (AngII-independent aldosterone) |

### Key Model Parameters

- **Adrenalectomy ARR normalisation** at ~3 months (renin rebound simulated)
- **Spironolactone IC₅₀ for MR**: 1.2 μg/L; canrenone IC₅₀: 0.8 μg/L
- **Finerenone IC₅₀**: 0.65 μg/L (highest MR affinity)
- **Cardiac fibrosis rate**: proportional to Aldo × MAP product; reduced by non-steroidal MRA

---

## Shiny Dashboard

**File:** `pa_shiny_app.R`

### Tabs

| Tab | Content |
|-----|---------|
| 1. Patient Profile | Disease subtype, severity sliders, drug doses, baseline labs |
| 2. RAAS / Drug PK | Renin-AngII-Aldosterone trajectories, plasma drug levels, MR occupancy |
| 3. Aldosterone Panel | PAC, ARR over time, CYP11B2 activity, ENaC activation |
| 4. Ion Homeostasis | Na⁺/K⁺/HCO₃⁻ dynamics, volume, GFR, summary table |
| 5. Cardiovascular & TOD | MAP, cardiac fibrosis, LVMi, GFR trajectories |
| 6. Scenario Comparison | 7 predefined scenarios, endpoint selector, summary table |
| 7. Biomarker Explorer | ARR panel, dose-response curves (ARR & K⁺ vs spironolactone dose) |

---

## References

**File:** `pa_references.md`  
48 PubMed-linked references covering:
- Disease epidemiology (Conn 1955 to Brown 2020)
- RAAS physiology and pathophysiology
- APA somatic mutations (KCNJ5, CACNA1D, ATP1A1, ATP2B3)
- CYP11B2 aldosterone synthase biology
- Renal ENaC/ROMK/SGK1 mechanisms
- ARR screening, AVS subtyping
- Cardiovascular complications
- Adrenalectomy outcomes
- MR antagonists (spironolactone, eplerenone, finerenone FIDELIO/FIGARO)
- QSP/mathematical modeling (Guyton, Hallow PKPD)

---

## Files

| File | Description |
|------|-------------|
| `pa_qsp_model.dot` | Graphviz source (120+ nodes, 10 clusters) |
| `pa_qsp_model.svg` | Vector diagram (~170 KB) |
| `pa_qsp_model.png` | Raster diagram (150 dpi) |
| `pa_mrgsolve_model.R` | 23-ODE mrgsolve model, 8 treatment scenarios |
| `pa_shiny_app.R` | 7-tab interactive Shiny dashboard |
| `pa_references.md` | 48 curated PubMed references |

---

## Quick Start

```r
# Run the Shiny app
shiny::runApp("pa_shiny_app.R")

# Run the mrgsolve simulation
source("pa_mrgsolve_model.R")
```

---

*Model compiled 2026-06-24. Part of the QSP Disease Model Library.*
