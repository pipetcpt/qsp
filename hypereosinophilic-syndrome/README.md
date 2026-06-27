# Hypereosinophilic Syndrome (HES) QSP Model

[![Disease](https://img.shields.io/badge/Disease-HES-red)]()
[![Framework](https://img.shields.io/badge/Framework-mrgsolve%20%C2%B7%20Shiny%20%C2%B7%20Graphviz-blue)]()
[![ODEs](https://img.shields.io/badge/ODEs-20_compartments-green)]()
[![Drugs](https://img.shields.io/badge/Drugs-4_(Mepolizumab%20%C2%B7%20Benralizumab%20%C2%B7%20Imatinib%20%C2%B7%20Prednisolone)-purple)]()

---

## Overview

**Hypereosinophilic Syndrome (HES)** is defined by persistent peripheral blood eosinophilia
≥ 1500 cells/µL for ≥ 1 month with evidence of eosinophil-mediated end-organ damage,
in the absence of a secondary cause. Despite being a rare disease, it carries significant
mortality risk from cardiac (Löffler endocarditis → endomyocardial fibrosis → restrictive
cardiomyopathy) and other end-organ complications.

This QSP model integrates:

- **IL-5 / eosinophil kinetics** from bone marrow progenitor to tissue infiltration
- **FIP1L1-PDGFRA clonal pathway** (myeloid HES subtype)
- **Mechanistic PK/PD** for mepolizumab (TMDD), benralizumab (ADCC), imatinib (TKI), and prednisolone (GR-mediated)
- **Progressive cardiac fibrosis** and **pulmonary infiltration** scoring

---

## Disease Subtypes Modeled

| Subtype | Key Driver | Modeled Feature |
|---------|-----------|----------------|
| Lymphocytic HES (L-HES) | Aberrant Th2 T-cell clones → IL-5 overproduction | `kprod_IL5` ↑, responsive to anti-IL-5 |
| Myeloid/Clonal HES (M-HES) | FIP1L1-PDGFRA fusion → constitutive TK activation | `CLONAL_HES=1`, `CLONAL_FOLD` parameter |
| Reactive HES | Parasites, atopy, malignancy | Baseline IL-5 elevation |
| Idiopathic HES (IHES) | Unknown | Default model parameters |

---

## Model Files

| File | Description |
|------|-------------|
| `hes_qsp_model.dot` | Graphviz mechanistic map (115+ nodes, 14 clusters) |
| `hes_qsp_model.svg` | Rendered SVG of mechanistic map |
| `hes_qsp_model.png` | Rendered PNG (150 dpi) |
| `hes_mrgsolve_model.R` | mrgsolve ODE model — 20 compartments, 5 treatment scenarios |
| `hes_shiny_app.R` | Shiny dashboard — 6 tabs |
| `hes_references.md` | 50 annotated PubMed references |

---

## Mechanistic Map Highlights

The `.dot` file (14 clusters, 115+ nodes) covers:

1. **Upstream Immune Activation** — TSLP, IL-33, IL-25 alarmins → ILC2, mDC, mast cells
2. **Th2 Cytokine Axis** — IL-5 (central), IL-4, IL-13, GM-CSF, eotaxins (CCL11/24/26)
3. **Myeloid Clonal HES** — FIP1L1-PDGFRA fusion, JAK2 V617F, RAS/MAPK, PI3K/AKT, STAT5
4. **BM Eosinophilopoiesis** — CD34+ → EoP → EoImm → EoMat-BM → Blood
5. **Peripheral Blood Eo Kinetics** — AEC dynamics, CCR3 trafficking, SIGLEC-8 apoptosis, ADCC
6. **Tissue Infiltration & Granule Toxins** — MBP, ECP, EPO, EDN, LTC4, PAF, ROS
7. **Cardiac Pathophysiology** — Löffler endocarditis → EMF → RCM, troponin, BNP, LVEF
8. **Pulmonary Involvement** — Eosinophilic pneumonia, DLCO, FEV1, pulmonary hypertension
9. **Drug PK: Anti-IL-5 Biologics** — Mepolizumab (2-CMT + TMDD), Benralizumab (2-CMT + ADCC), Reslizumab
10. **Drug PK: Small Molecules** — Imatinib (1-CMT PO), Prednisolone (1-CMT PO)
11. **Glucocorticoid PD** — GR transactivation/transrepression, IL-5 suppression, eosinophil apoptosis
12. **GI & Skin Manifestations** — Eosinophilic esophagitis, gastroenteritis, urticaria
13. **Neurological Involvement** — Peripheral neuropathy, encephalopathy, CNS thromboembolism
14. **Clinical Endpoints** — AEC response, HES remission, organ damage score, QoL

---

## ODE Model Compartments (20 total)

### Drug PK (8 CMTs)
| Compartment | Description |
|-------------|-------------|
| `MEPO_DEPOT` | Mepolizumab SC depot (µg) |
| `MEPO_C1` | Mepolizumab central (µg) |
| `MEPO_C2` | Mepolizumab peripheral (µg) |
| `TMDD` | Mepolizumab–IL-5 TMDD complex (nM) |
| `BENRA_DEPOT` | Benralizumab SC depot (µg) |
| `BENRA_C1` | Benralizumab central (µg) |
| `BENRA_C2` | Benralizumab peripheral (µg) |
| `IMAT_GUT` | Imatinib gut absorption (mg) |
| `IMAT_C` | Imatinib central (mg/L) |
| `PRED_GUT` | Prednisolone gut absorption (mg) |
| `PRED_C` | Prednisolone central (mg/L) |

### Disease PD (9 CMTs)
| Compartment | Description |
|-------------|-------------|
| `IL5` | Serum IL-5 (pg/mL) |
| `EoP` | BM eosinophil progenitors |
| `EoI` | Immature BM eosinophils |
| `EoM_BM` | Mature BM eosinophils |
| `EO_BLOOD` | Peripheral blood AEC (cells/µL) |
| `FIBROSIS` | Cardiac fibrosis score (0–1) |
| `PULM_SCORE` | Pulmonary infiltration score (0–1) |

---

## Treatment Scenarios (5)

| Scenario | Drug | Dose / Regimen | HES Type |
|----------|------|----------------|----------|
| 1 | Untreated | — | Reactive (AEC 3000) |
| 2 | Prednisolone | 70 mg/day PO → taper | Reactive |
| 3 | Mepolizumab | 300 mg SC q4w | Reactive |
| 4 | Benralizumab | 30 mg SC q4w × 3 → q8w | Reactive |
| 5 | Imatinib | 100 mg/day PO | Clonal (FIP1L1-PDGFRA+, AEC 5000) |

---

## Key PK Parameters

| Drug | Route | t½ | SC F | Vc | CL | MOA |
|------|-------|-----|------|----|----|-----|
| Mepolizumab | SC | ~22 d | 80% | 3.5 L | 7.2 mL/h | Anti-IL-5 (TMDD) |
| Benralizumab | SC | ~15 d | 59% | 3.0 L | 4.0 mL/h | Anti-IL-5Rα (ADCC) |
| Imatinib | PO | ~18 h | 98% | 110 L | 12 L/h | PDGFRA/ABL TKI |
| Prednisolone | PO | ~20 h | 82% | 38 L | 3.5 L/h | GR agonist |

---

## Shiny Dashboard Tabs (6)

| Tab | Content |
|-----|---------|
| 1. Patient Profile | Disease overview, subtype, MoA comparison panels |
| 2. Pharmacokinetics | Drug concentration profiles, PK table, AUC/Cmax/Tmax metrics |
| 3. Eosinophil Kinetics | AEC time-course, BM compartments, IL-5 dynamics, % change |
| 4. Organ Damage | Cardiac fibrosis progression, pulmonary infiltration, severity interpretation |
| 5. Treatment Comparison | All 5 scenarios overlaid, week-24 summary table |
| 6. Biomarker Correlations | AEC vs fibrosis, IL-5 vs AEC, phase-space scatter plots |

---

## Usage

```r
# Install dependencies
install.packages(c("mrgsolve", "dplyr", "ggplot2", "patchwork", "tidyr"))

# Run the mrgsolve simulation
source("hes_mrgsolve_model.R")

# Launch interactive Shiny dashboard
install.packages(c("shiny", "shinydashboard", "plotly", "DT"))
shiny::runApp("hes_shiny_app.R")
```

```bash
# Render mechanistic map
dot -Tsvg hes_qsp_model.dot -o hes_qsp_model.svg
dot -Tpng -Gdpi=150 hes_qsp_model.dot -o hes_qsp_model.png
```

---

## Clinical Context

- **Diagnosis**: Peripheral AEC ≥ 1500/µL × 1 month + end-organ damage + exclusion of reactive causes
- **HES remission target**: AEC < 300 cells/µL (deep response); < 1500/µL (disease control)
- **Cardiac screening**: Echocardiography + troponin + NT-proBNP at baseline and q6 months
- **FIP1L1-PDGFRA testing**: FISH or RT-PCR in all HES patients; if positive → imatinib first-line

---

## References

See [`hes_references.md`](./hes_references.md) for 50 annotated citations covering disease pathogenesis,
clinical trials (MEPO, BENRA, imatinib), and PK/PD modeling.

Key papers:
- Rothenberg ME et al. *N Engl J Med* 2008 (mepolizumab)
- Cools J et al. *N Engl J Med* 2003 (FIP1L1-PDGFRA imatinib)
- Kuang FL et al. *N Engl J Med* 2019 (benralizumab)
- Gleich GJ et al. *Lancet* 2002 (imatinib for HES)
- Roufosse FA et al. *J Allergy Clin Immunol* 2020 (mepolizumab phase III)

---

*Generated by Claude Code Routine (CCR) — 2026-06-27*
