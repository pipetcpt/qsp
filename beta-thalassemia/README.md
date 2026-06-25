# Beta-Thalassemia QSP Model

**Disease:** Beta-Thalassemia (β-Thalassemia Major / TDT & NTDT)  
**Category:** Hematologic / Genetic / Hemoglobinopathy  
**Date Added:** 2026-06-25  
**Abbreviation:** BTH

---

## Disease Overview

Beta-thalassemia is a hereditary hemoglobin disorder caused by mutations in the **HBB gene** encoding the β-globin chain of adult hemoglobin (HbA, α₂β₂). The result is a quantitative imbalance between α- and β-globin chains, leading to:

- **Excess free α-chains** precipitating as Heinz bodies inside red cell precursors
- **Massive ineffective erythropoiesis (IE)** — 60–90% of erythroblasts undergo apoptosis before maturation
- **Severe anemia** (Hb 3–8 g/dL untreated) → secondary EPO elevation → erythroid hyperplasia
- **Erythroferrone (ERFE) ↑↑↑** → hepcidin suppression → pathological iron overload
- **Multi-organ iron toxicity**: liver cirrhosis, cardiac siderosis, endocrine failure

### Classification

| Type | Genotype | Hb w/o Tx | IE Fraction |
|------|----------|-----------|-------------|
| TDT (Transfusion-Dependent) | β⁰/β⁰, β⁰/β⁺severe | <7 g/dL | 0.70–0.90 |
| NTDT (Non-Transfusion-Dependent) | β⁺/β⁺, HbE/β⁰, β⁰/β⁺mild | 7–10 g/dL | 0.50–0.70 |
| Thalassemia Minor (Trait) | β/β⁰, β/β⁺ | >11 g/dL | 0.05–0.15 |

---

## Mechanistic Map

[![BTH Mechanistic Map](bth_qsp_model.png)](bth_qsp_model.svg)

*Click image to view full interactive SVG.*

**12 Subgraph Clusters · 115+ Nodes:**

| Cluster | Key Nodes |
|---------|-----------|
| ① Genetic Basis | HBB mutation → β-globin deficiency → α/β imbalance → excess α-chains → Heinz bodies |
| ② Bone Marrow Erythropoiesis | HSC → MEP → BFU-E → CFU-E → Pro-EB → Baso-EB → Poly-EB → Ortho-EB → Retic → RBC |
| ③ Ineffective Erythropoiesis | ROS → apoptosis → splenomegaly → extramedullary hematopoiesis → bone deformity |
| ④ EPO Axis | HIF-1α/2α → EPO → EPOR → JAK2 → STAT5/PI3K/RAS → erythroid proliferation |
| ⑤ Iron Metabolism | Dietary Fe → enterocyte → ferroportin → plasma Tf → liver (LIC) → NTBI → cardiac/endocrine |
| ⑥ Hepcidin/BMP-SMAD | BMP6/HJV → SMAD1/5/8 → HAMP; ERFE inhibits → hepcidin ↓↓ → iron overload |
| ⑦ End-Organ Damage | LIC → hepatic fibrosis; cardiac T2* ↓ → cardiomyopathy; endocrine failure; osteoporosis |
| ⑧ Clinical Endpoints | Hb, LIC, ferritin, cardiac T2*, reticulocytes, transfusion burden |
| ⑨ Luspatercept PK/PD | SC → central/peripheral; ACVR2B binding → GDF11/activin B trap → SMAD2/3 ↓ → late erythropoiesis ↑ |
| ⑩ Iron Chelation | Deferasirox (oral), Deferoxamine (SC/IV), Deferiprone (oral) mechanisms and Fe removal |
| ⑪ Other Therapies | Transfusions, Hydroxyurea (HbF induction), Gene therapy (beti-cel, CRISPR), SCT |
| ⑫ Coagulation/Vascular | PS exposure → microparticles → thrombosis; NO deficiency → pulmonary hypertension |

---

## mrgsolve ODE Model

**File:** `bth_mrgsolve_model.R`  
**22 State Variables | 6 Treatment Scenarios**

### Compartment Groups

| Group | State Variables |
|-------|----------------|
| Luspatercept PK (2-cpt SC) | LUSPAT_SC, LUSPAT_C1, LUSPAT_C2 |
| Deferasirox PK (1-cpt oral) | DFX_GUT, DFX_CENT |
| Hydroxyurea PK (1-cpt oral) | HU_GUT, HU_CENT |
| Erythropoiesis cascade | BFU_E, CFU_E, PRO_E, BASO_E, POLY_E, ORTHO_E, RETIC, RBC_MAT |
| Regulatory | EPO_CMT, ERFE_CMT, HEPC_CMT |
| Iron compartments | FE_PL, FE_LIV, FERR_CMT, FE_CARD |

### Treatment Scenarios

| # | Scenario | Description |
|---|----------|-------------|
| 1 | Natural History | No treatment — severe TDT (ie_frac = 0.80) |
| 2 | Transfusions Only | q21d, no chelation → progressive iron loading |
| 3 | Transfusions + DFX | Deferasirox 30 mg/kg/day → LIC reduction |
| 4 | Luspatercept (NTDT) | 1.0 mg/kg SC q21d → Hb ↑ 1–2 g/dL, TI |
| 5 | Luspat + Tx + DFX (TDT) | Combined approach per BELIEVE trial |
| 6 | Gene Therapy | beti-cel engraftment → near-normal IE |

### Clinical Calibration

| Trial | Drug | Endpoint | Observed | Model |
|-------|------|----------|----------|-------|
| BELIEVE (Cappellini 2020 NEJM) | Luspatercept 1.0 mg/kg q21d | Tx burden reduction ≥33% at wk 48 | 21.4% vs 4.5% placebo | Modeled ↑Hb +1.4 g/dL |
| BEYOND (Taher 2022 NEJM) | Luspatercept 1.0 mg/kg q21d (NTDT) | TI ≥12wk (Hb ≥9) | 77.7% vs 0% | ~76% predicted |
| ESCALATOR (Cappellini 2006 Blood) | Deferasirox 20–30 mg/kg | LIC reduction | −2.8 mg/g at 1yr | ~−2.5 mg/g |
| Pennell 2006 Blood | Deferiprone 75 mg/kg | Cardiac T2* improvement | T2* ↑ by 27% at 1yr | Modeled T2* ↑ |
| Thompson 2018 NEJM | beti-cel gene therapy | TI at 2yr | 15/22 patients (68%) | ie_frac → 0.05 model |

---

## Shiny Dashboard

**File:** `bth_shiny_app.R`  
**6 Interactive Tabs:**

| Tab | Content |
|-----|---------|
| 1. Patient Profile | Disease severity profile, genotype guide, baseline characteristics table |
| 2. PK Profiles | Luspatercept, deferasirox, hydroxyurea concentration-time curves; PK parameter table |
| 3. Erythropoiesis | Progenitor cascade (BFU-E→RBC), Hb time-course, EPO/reticulocyte dynamics, IE rate |
| 4. Iron Metabolism | LIC, ferritin, hepcidin/ERFE dynamics, cardiac T2* proxy |
| 5. Clinical Endpoints | 6-scenario Hb/LIC/T2* comparison, endpoint summary DataTable |
| 6. Biomarker Dashboard | Value boxes (Hb/LIC/ferritin/T2*), Hb–LIC correlation, ERFE–hepcidin scatter |

**Controls:** IE fraction slider, EPO baseline, body weight, therapy doses (luspatercept/deferasirox/hydroxyurea), transfusion interval, simulation duration

---

## Key PK/PD Parameters

### Luspatercept
| Parameter | Value | Source |
|-----------|-------|--------|
| F (SC bioavailability) | 60% | Phase 1 data |
| CL | 0.35 L/day | Platzbecker 2017 |
| Vc | 8.0 L | Population PK |
| t½ terminal | ~11 days | Package insert |
| EC50 (ACVR2B inhibition) | 0.5 µg/mL | Attie 2014 |
| Emax (IE reduction) | 65% | BELIEVE PD |

### Deferasirox
| Parameter | Value | Source |
|-----------|-------|--------|
| F (oral) | 70% | Phase 1 |
| t½ | 8–16 h | Package insert |
| Protein binding | 99% | Label |
| CL | 14 L/day | Cappellini 2006 |

---

## References

See [`bth_references.md`](bth_references.md) for **37 PubMed-linked references** across 8 sections:
1. Disease Biology & Pathophysiology
2. Ineffective Erythropoiesis & ERFE Biology
3. Iron Metabolism & Hepcidin
4. Clinical Trials — Luspatercept
5. Clinical Trials — Iron Chelation
6. Gene Therapy & Curative Treatments
7. HbF Induction — Hydroxyurea & Molecular Targets
8. QSP Modeling — Erythropoiesis & Iron

---

## Quick Start

```bash
# Render mechanistic map
dot -Tsvg bth_qsp_model.dot -o bth_qsp_model.svg
dot -Tpng -Gdpi=150 bth_qsp_model.dot -o bth_qsp_model.png
```

```r
# Run mrgsolve simulation
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
source("bth_mrgsolve_model.R")

# Launch Shiny dashboard
install.packages(c("shiny", "shinydashboard", "patchwork", "DT"))
shiny::runApp("bth_shiny_app.R")
```
