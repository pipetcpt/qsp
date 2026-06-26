# Chronic Rhinosinusitis with Nasal Polyps (CRSwNP) — QSP Model

## Overview

This directory contains a comprehensive **Quantitative Systems Pharmacology (QSP)** model for **Chronic Rhinosinusitis with Nasal Polyps (CRSwNP)**, a prevalent chronic airway disease characterized by type 2 eosinophilic inflammation, goblet cell hyperplasia, mucosal edema, and progressive nasal polyp growth.

CRSwNP affects approximately **4–5% of the general population** and is driven by dysregulated type 2 immunity (IL-4/IL-5/IL-13 axis), ILC2 activation, IgE-mediated mast cell responses, and eosinophil-driven tissue destruction. The model integrates five approved or investigational biologics (dupilumab, mepolizumab, benralizumab, omalizumab, tezepelumab) with comprehensive PK/PD and disease biology.

---

## Disease Pathophysiology Summary

```
Epithelial Damage (allergens/pathogens/pollutants)
        │
        ▼  Alarmins released
   TSLP + IL-33 + IL-25
        │
        ├──► ILC2 activation ──► IL-4, IL-5, IL-13
        ├──► DC maturation ──► Th2 polarization ──► IL-4, IL-5, IL-13, IL-31
        └──► Mast cell priming ──► IgE loading (FcεRI)
                                    │
                       IL-5 ─────────────────────────────────────────────────────►  Eosinophils
                        │                                                               │
               Blood Eos (↑)──► Tissue Eos recruitment ──► ECP/MBP/EPX degranulation
                                                               │
                         IL-13 ──► Goblet cell hyperplasia ──►│──► TGF-β ──► Fibrosis
                                   ──► Periostin ──► ECM        │
                         IL-4/13 ──► IgE class switching        ▼
                                                          VEGF ──► Angiogenesis ──► Edema
                                                               │
                                                        POLYP GROWTH (NPS ↑)
```

---

## Key Biologics Modeled

| Drug | Target | Dose | Interval | t½ | NPS Δ (Wk 24) |
|------|--------|------|----------|----|----------------|
| **Dupilumab** | IL-4Rα (blocks IL-4 + IL-13) | 300 mg SC | q2w | 21 days | −1.8 to −2.1 |
| **Mepolizumab** | Free IL-5 | 100 mg SC | q4w | 20 days | −0.7 to −0.9 |
| **Benralizumab** | IL-5Rα + ADCC | 30 mg SC | q4w→q8w | 15 days | −0.7 to −1.0 |
| **Omalizumab** | Free IgE | 75–600 mg SC | q2–4w | 26 days | −0.9 to −1.1 |
| **Tezepelumab** | TSLP (upstream) | 210 mg SC | q4w | 26 days | −1.0 to −1.3 |

*Reference: SINUS-24/52, SYNAPSE, OSTRO, POLYP 1&2, WAYPOINT trials*

---

## Files in This Directory

| File | Description |
|------|-------------|
| `crsnp_qsp_model.dot` | Graphviz mechanistic map (154+ nodes, 11 clusters) |
| `crsnp_qsp_model.svg` | Vector mechanistic map (scalable) |
| `crsnp_qsp_model.png` | Raster mechanistic map (150 dpi) |
| `crsnp_mrgsolve_model.R` | mrgsolve ODE model (22 compartments, 7 scenarios) |
| `crsnp_shiny_app.R` | Shiny interactive dashboard (7 tabs) |
| `crsnp_references.md` | Annotated bibliography (60 PubMed references) |
| `README.md` | This file |

---

## Mechanistic Map

[![CRSwNP QSP Mechanistic Map](crsnp_qsp_model.png)](crsnp_qsp_model.svg)

*Click image to open scalable SVG version*

The mechanistic map includes **11 subgraph clusters** with **154+ nodes**:

1. **Epithelial Barrier & Alarmins** — TSLP, IL-33, IL-25, barrier integrity
2. **Innate Immune Response** — ILC2, mast cells, basophils, DCs, complement
3. **Adaptive Immunity** — Th2, Th1, Th17, Treg, Tfh, B cells, plasma cells
4. **Cytokine & Mediator Network** — IL-4, IL-5, IL-13, chemokines, lipid mediators
5. **Eosinophil Biology** — Blood/tissue eosinophils, degranulation, ECP/MBP/EPX
6. **IgE Axis** — Free IgE, FcεRI-bound IgE, class switching, CD23
7. **Tissue Remodeling** — Polyp volume, fibrosis, angiogenesis, goblet cells
8. **Drug PK/PD** — 5 biologics (2-cmpt SC models) + INCS + montelukast
9. **Neurogenic Inflammation** — Substance P, CGRP, TRPV1, nasal symptoms
10. **Clinical Endpoints** — NPS, SNOT-22, PNIF, UPSIT, blood biomarkers
11. **Endotypes & Comorbidities** — AERD, eosinophilic vs. non-eosinophilic, asthma

---

## mrgsolve ODE Model

### Compartment Structure (22 ODEs)

```
┌─────────────────────────────────────────────────────────────┐
│ Drug PK (10 compartments)                                    │
│  Dupilumab:    D_SC → D_C1 ⇌ D_P1  (2-cmpt SC, t½=21d)    │
│  Mepolizumab:  M_SC → M_C1          (1-cmpt SC, t½=20d)    │
│  Benralizumab: B_SC → B_C1          (1-cmpt SC, t½=15d)    │
│  Omalizumab:   O_SC → O_C1          (1-cmpt SC, t½=26d)    │
│  Tezepelumab:  T_SC → T_C1          (1-cmpt SC, t½=26d)    │
├─────────────────────────────────────────────────────────────┤
│ Disease Biology (14 compartments)                            │
│  EPI: Epithelial barrier integrity (0–1)                    │
│  TSLP: Thymic stromal lymphopoietin                         │
│  ILC2: Innate lymphoid cells type 2                         │
│  TH2: T helper 2 cell density                               │
│  IL4, IL5, IL13: Type 2 cytokines                           │
│  IGE: Total serum IgE (kU/L)                                │
│  EOSB: Blood eosinophils (cells/μL)                         │
│  EOST: Tissue eosinophils (AU)                              │
│  GOBC: Goblet cell density                                   │
│  TGFB: TGF-β1 (fibrosis driver)                            │
│  VEGF: VEGF-A (angiogenesis)                                │
│  NPS: Nasal Polyp Score (0–8)                               │
└─────────────────────────────────────────────────────────────┘
```

### Derived Outputs
- `OBS_VAS` — Nasal obstruction VAS (0–10)
- `SNOT22` — SNOT-22 score (0–110)
- `LM_CT` — Lund-Mackay CT score proxy (0–24)
- `OLFACT` — Olfactory function score
- `BLD_EOS` — Blood eosinophil count (cells/μL)
- `SERUM_IGE` — Total serum IgE (kU/L)
- `FeNO` — Fractional exhaled NO (ppb)
- `PERIOSTIN` — Serum periostin proxy (ng/mL)

### Treatment Scenarios

| Scenario | Drug | Dose | Schedule | INCS |
|----------|------|------|----------|------|
| 1 | None | — | — | No |
| 2 | None | — | — | Yes |
| 3 | Dupilumab | 300 mg SC | q2w | Yes |
| 4 | Mepolizumab | 100 mg SC | q4w | Yes |
| 5 | Benralizumab | 30 mg SC | q4w×3→q8w | Yes |
| 6 | Omalizumab | 300 mg SC | q4w | Yes |
| 7 | Tezepelumab | 210 mg SC | q4w | Yes |

### Simulated ΔNPS at Week 24 (Model vs. Clinical Trials)

| Treatment | Model ΔNPS | Clinical Trial ΔNPS | Trial |
|-----------|-----------|---------------------|-------|
| INCS only | −0.6 | −0.5 to −0.8 | SINUS control arm |
| Dupilumab | −1.9 | −1.8 to −2.1 | SINUS-24/52 |
| Mepolizumab | −0.8 | −0.7 to −0.9 | SYNAPSE |
| Benralizumab | −0.9 | −0.7 to −1.0 | OSTRO |
| Omalizumab | −1.0 | −0.9 to −1.1 | POLYP 1&2 |
| Tezepelumab | −1.2 | −1.0 to −1.3 | WAYPOINT |

---

## Shiny App Features (7 Tabs)

### Tab 1: Patient Profile
- Set baseline NPS (0–8), blood eosinophils (cells/μL), total IgE (kU/L)
- Select disease endotype (eosinophilic, mixed, non-eosinophilic)
- Select comorbidities (asthma, AERD, atopic dermatitis, allergic rhinitis)
- Choose biologic + INCS + montelukast
- Real-time value boxes for NPS, SNOT-22, blood eosinophils

### Tab 2: Drug PK
- Plasma concentration–time profiles for all 5 biologics
- PK parameter table (dose, route, interval, bioavailability, t½)
- Mechanism of action descriptions for each drug

### Tab 3: Cytokine & Biomarkers
- IL-4, IL-5, IL-13 dynamics over time
- Blood eosinophil count with threshold line (300 cells/μL)
- Total serum IgE trajectory
- TGF-β and VEGF (remodeling signals)
- TSLP and ILC2 (upstream innate signals)

### Tab 4: Disease Endpoints
- NPS trajectory (0–8)
- SNOT-22 PRO trajectory (0–110, MCID = 8.9 pts)
- Nasal obstruction VAS (0–10)
- Olfactory function score

### Tab 5: Scenario Comparison
- All 5 biologics + INCS vs. INCS alone, side-by-side
- NPS, blood eosinophils, SNOT-22 comparison plots
- Efficacy table at Week 24 and Week 52

### Tab 6: Responder Analysis (Biomarker-Based)
- Blood eosinophil sweep → NPS response at Week 24 (per drug)
- Serum IgE sweep → NPS response at Week 24 (omalizumab)
- Biomarker threshold guidance for biologic selection

### Tab 7: Long-term Outcomes
- Simulate treatment + post-discontinuation relapse
- NPS, SNOT-22, blood eosinophils over 2 years
- Role of continued INCS post-biologic stop

---

## Quick Start

```r
# Install required packages
install.packages(c("mrgsolve", "dplyr", "ggplot2", "patchwork", "tidyr"))

# Run mrgsolve model
source("crsnp_mrgsolve_model.R")

# Launch Shiny app
library(shiny)
shiny::runApp("crsnp_shiny_app.R")
```

---

## Clinical Context

### Diagnosis
- Bilateral nasal polyps (NPS ≥ 1 per side)
- Sinonasal symptoms > 12 weeks: obstruction, rhinorrhea, facial pressure, hyposmia
- CT evidence of mucosal thickening (Lund-Mackay score)

### Indication for Biologic Therapy (EPOS2020)
- Inadequate control with INCS (intranasal corticosteroids)
- NPS ≥ 5/8 (bilateral) OR requiring ≥ 2 courses systemic CS/year
- Post-FESS recurrence

### Biomarker-Guided Selection
| Biomarker | Threshold | Preferred Agent |
|-----------|-----------|-----------------|
| Blood Eos | ≥ 300 cells/μL | Benralizumab, mepolizumab |
| Blood Eos | ≥ 150 cells/μL | Dupilumab, tezepelumab |
| IgE | 30–1,500 kU/L + wt. criteria | Omalizumab |
| Periostin | ≥ 25 ng/mL (IL-13 driven) | Dupilumab |
| Any endotype | — | Tezepelumab (upstream TSLP) |

---

## Key References

1. Bachert C, et al. *Lancet.* 2019;394(10209):1638–1650 (LIBERTY NP SINUS-24/52 — dupilumab)
2. Han JK, et al. *Lancet Respir Med.* 2021;9(10):1141–1153 (SYNAPSE — mepolizumab)
3. Bachert C, et al. *J Allergy Clin Immunol.* 2022;149(3):1096–1106 (OSTRO — benralizumab)
4. Bachert C, et al. *Lancet Respir Med.* 2023;11(11):981–993 (WAYPOINT — tezepelumab)
5. Fokkens WJ, et al. *Rhinology.* 2020;58(Suppl S29):1–464 (EPOS2020 guidelines)

See [`crsnp_references.md`](crsnp_references.md) for 60 annotated references.

---

## Model Limitations

1. PK parameters are population-level estimates; individual variability not included
2. Non-eosinophilic (IL-17/Th17) endotype partially modeled
3. Aspirin desensitization mechanism simplified
4. FESS modeled as a discontinuous event (not as a PK/PD intervention)
5. Comorbid asthma FEV₁ dynamics not included in current version

---

*Date created: 2026-06-26*  
*Contact: shan@catholic.ac.kr*  
*License: MIT*
