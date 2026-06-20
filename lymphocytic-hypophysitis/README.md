# Lymphocytic Hypophysitis (림프구성 뇌하수체염) — QSP Model

> **Disease**: Lymphocytic Hypophysitis (LyH)  
> **Category**: Autoimmune (Pituitary)  
> **ICD-10**: E23.6 (Other conditions of pituitary gland)  
> **Model version**: 1.0 · Date: 2026-06-20

---

## Overview

Lymphocytic Hypophysitis is a rare organ-specific autoimmune disease characterized by
lymphocytic infiltration of the pituitary gland leading to progressive destruction of
pituitary cells and **multi-axis hypopituitarism**. It is the most common form of primary
hypophysitis and predominantly affects women in the peripartum period. Immune checkpoint
inhibitors (anti-CTLA4, anti-PD1) represent a major iatrogenic cause.

### Key Pathophysiology

```
Triggers (pregnancy / ICI / infection)
  ↓ molecular mimicry + HLA-DR3/4 susceptibility
APCs → MHC-II antigen presentation (PIT-1, RABPHILIN-3L)
  ↓
CD4+ T cell polarization (Th1/Th17) + CD8+ cytotoxic T cells
  ↓ B cell activation → plasma cells → Anti-Pituitary Antibodies (APA)
     + ADCC + complement deposition
Lymphocytic infiltration → Pituitary mass → Hormone deficiencies:
  • ACTH ↓  → secondary adrenal insufficiency (~70%)
  • TSH ↓   → central hypothyroidism (~40%)
  • FSH/LH ↓ → hypogonadism (~40%)
  • GH ↓    → adult GH deficiency (~30%)
  • PRL ↑   → stalk compression hyperprolactinemia (~50%)
  • ADH ↓   → diabetes insipidus (posterior, ~25%)
```

---

## Mechanistic Map

[![QSP Mechanistic Map](lhyp_qsp_model.png)](lhyp_qsp_model.svg)

> Click the image to open the full-resolution SVG

**Map Statistics:**
- **154+ nodes** across 14 subgraph clusters
- **14 biological subgraph clusters** covering:
  1. Disease Triggers & Genetics
  2. Innate Immune System
  3. Adaptive Immune System (T & B cells)
  4. Cytokine & Chemokine Network
  5. Pituitary Anatomy & Pathology
  6. Hypothalamic Control
  7. HPA Axis
  8. HPT Axis
  9. HPG Axis
  10. GH/IGF-1 Axis
  11. Prolactin & ADH
  12. Drug PK
  13. Drug PD
  14. Clinical Endpoints

---

## mrgsolve ODE Model

**File**: `lhyp_mrgsolve_model.R`

### Model Structure: 22 ODE Compartments

| Subsystem | Compartments |
|-----------|-------------|
| **Drug PK** | Pred_gut, Pred_central, Pred_periph, AZA_gut, AZA_plasma, RTX_plasma |
| **Immune** | Naive T (Tn), Effector T (Te), Regulatory T (Tr), Naive B (Bn), Plasma cells (Bp), APA |
| **Pituitary** | PitInf (inflammatory mass), PitFunc (functional cell mass) |
| **HPA Axis** | ACTH, Cortisol |
| **HPT Axis** | TSH, fT4 |
| **GH Axis** | GH, IGF1 |
| **HPG Axis** | FSH, LH, E2 |
| **PRL/ADH** | Prolactin, ADH |

### Key Equations

**Cortisol Negative Feedback on ACTH (Hill kinetics):**
```
ACTH_prod = kprod_ACTH × [1 - Emax × Cortisol^n / (IC50^n + Cortisol^n)] × PitFunc
```

**Pituitary Functional Cell Mass:**
```
dPitFunc/dt = -kpit_inflam × (Te + APA×0.01) × (1 - ImmunoSupp) + kpit_repair × PitFunc
```

**Combined Immunosuppression:**
```
ImmunoSupp = 1 - (1 - E_pred) × (1 - E_aza × 0.5)
```

### 5 Treatment Scenarios

| # | Scenario | Description |
|---|----------|-------------|
| S1 | No Treatment | Natural history over 2 years |
| S2 | High-Dose Prednisolone | 60→40→30→20→10→5 mg/day taper |
| S3 | Pred + Azathioprine | Steroid-sparing combination |
| S4 | Rituximab | ICI-associated LyH, anti-CD20 therapy |
| S5 | Low-Dose Pred | Suboptimal 10 mg/day maintenance |

### Calibration Data

| Parameter | Source | Value |
|-----------|--------|-------|
| Prednisolone PK | Czock et al. Clin Pharmacokinet 2005 | ka=1.2/h, Vc=45L, CL=8.5L/h |
| ACTH deficiency prevalence | Honegger et al. JCEM 2015 | ~70% |
| TSH deficiency prevalence | Honegger et al. JCEM 2015 | ~40% |
| ICI-LyH steroid response | Faje et al. JCEM 2014 | ~50% remission |
| Rituximab Emax (B-depletion) | Maloney et al. Blood 1997 | 0.90 |

---

## Shiny App

**File**: `lhyp_shiny_app.R`

### 7 Interactive Tabs

| Tab | Content |
|-----|---------|
| **1. Patient Profile** | Disease severity settings, treatment selection, LyH subtype |
| **2. Drug PK** | Prednisolone, azathioprine, rituximab plasma concentrations |
| **3. Pituitary Function** | Function score, mass index, hormone deficiency heatmap |
| **4. Hormone Axes** | HPA (ACTH/Cortisol), HPT (TSH/fT4), GH/IGF-1, HPG (FSH/LH/E2) |
| **5. Immune Dynamics** | T cells, B cells, APA, immune-pituitary phase plot |
| **6. Clinical Endpoints** | Radar chart, adrenal crisis risk, summary table |
| **7. Biomarker Trajectories** | Any variable with reference ranges, PRL/ADH |

```r
# Launch the app
shiny::runApp("lhyp_shiny_app.R")
```

---

## References

**File**: `lhyp_references.md` — **60 references** in 9 sections

| Section | Count |
|---------|-------|
| Disease Overview & Epidemiology | 5 |
| Pathophysiology & Immunology | 7 |
| Clinical Presentation & Subtypes | 4 |
| Hormone Axis Dynamics & QSP | 11 |
| Treatment & Pharmacology | 8 |
| Checkpoint Inhibitor-Induced LyH | 4 |
| Biomarkers & Imaging | 5 |
| Mathematical & Systems Biology | 7 |
| Quality of Life & Outcomes | 5 |

---

## Files

```
lymphocytic-hypophysitis/
├── README.md                   ← This file
├── lhyp_qsp_model.dot          ← Graphviz mechanistic map (154+ nodes, 14 clusters)
├── lhyp_qsp_model.svg          ← SVG vector (scalable)
├── lhyp_qsp_model.png          ← PNG raster (150 dpi)
├── lhyp_mrgsolve_model.R       ← ODE model (22 compartments, 5 scenarios)
├── lhyp_shiny_app.R            ← Interactive Shiny dashboard (7 tabs)
└── lhyp_references.md          ← 60 references in 9 categories
```

---

## Key Clinical Insights from Model

1. **ACTH deficiency is the most common and earliest deficit** (~70% of LyH patients).
   Cortisol nadir often precedes clinical symptoms — screening is essential.

2. **Hyperprolactinemia from stalk compression** may precede overt hormone deficiency,
   providing an early biomarker of pituitary infiltration.

3. **High-dose prednisolone taper achieves optimal outcomes** by simultaneously
   suppressing lymphocytic infiltration AND reducing pituitary mass.

4. **Azathioprine (steroid-sparing)** reduces cumulative corticosteroid exposure with
   comparable immunosuppressive efficacy over long-term follow-up.

5. **ICI-associated LyH responds poorly to corticosteroids** alone — rituximab may be
   superior by depleting B cells producing APA.

6. **Posterior pituitary involvement** (ADH deficiency → DI) is a poor prognostic marker
   indicating more severe infiltration and less reversible damage.

---

## Dependencies

```r
# Required R packages
install.packages(c(
  "mrgsolve",
  "dplyr",
  "ggplot2",
  "tidyr",
  "shiny",
  "shinydashboard",
  "plotly",
  "DT"
))
```

---

*Part of the QSP Disease Model Library — automated daily model generation via Claude Code Routine*  
*Model: Lymphocytic Hypophysitis v1.0 | Generated: 2026-06-20*
