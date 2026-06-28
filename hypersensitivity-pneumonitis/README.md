# Hypersensitivity Pneumonitis (HP) — QSP Model

[![Model](https://img.shields.io/badge/QSP-mrgsolve%20%7C%20Shiny%20%7C%20Graphviz-blue)]()
[![Disease](https://img.shields.io/badge/Category-Chronic%20Lung%20Disease-green)]()
[![Status](https://img.shields.io/badge/Status-Complete-success)]()

## Disease Overview

**Hypersensitivity Pneumonitis (HP)** — also called Extrinsic Allergic Alveolitis — is an immune-mediated interstitial lung disease triggered by repeated inhalation of specific organic antigens (avian proteins, thermophilic bacteria, mold spores). The disease spans a spectrum from:

- **Acute HP**: flu-like symptoms 4–8h after exposure; fully reversible with antigen avoidance
- **Subacute HP**: insidious onset; BAL lymphocytosis + non-caseating granulomas
- **Chronic/Fibrotic HP**: progressive lung fibrosis resembling UIP/NSIP pattern; poor prognosis

### Epidemiology
- Prevalence: 3–30 per 100,000 (varies by antigen exposure type)
- Common causes: Bird-fancier's lung, Farmer's lung, Hot-tub lung, Humidifier lung
- 5-year mortality for fibrotic HP: ~20–40% (UIP pattern = worst prognosis)
- Nintedanib approved (2022) for progressive fibrosing ILD including fibrotic HP

---

## QSP Model Architecture

### Key Biological Pathways Modeled (10 Clusters, 100+ nodes)

| Cluster | Components |
|---------|-----------|
| ① Antigen Exposure | Inhaled antigen → alveolar deposition → soluble/particulate forms, mucociliary clearance |
| ② Innate Immunity | TLR2/4/9 → NF-κB → M1 macrophage, NLRP3 inflammasome, neutrophil recruitment |
| ③ Cytokine Network | TNF-α, IL-1β, IL-6, IL-12, IL-18, IFN-γ, TGF-β1, IL-17A, IL-10, CCL18, CXCL8, MMP-7 |
| ④ Adaptive Immunity | mDC → Lymph node → Th1/Th17/Treg differentiation; B cell → precipitin IgG, IC formation |
| ⑤ Granuloma | Non-caseating granuloma formation, lymphocyte cuffing, IL-10-driven resolution vs. fibrosis |
| ⑥ Fibrosis Cascade | AEC2 injury, EMT, fibroblast→myofibroblast, collagen deposition, ECM cross-linking, ROS |
| ⑦ Lung Pathophysiology | Alveolar wall thickening, UIP/NSIP patterns, FVC/DLCO decline, V/Q mismatch, pulm HTN |
| ⑧ Drug PK | Prednisolone (oral, 2-compt.), MMF/MPA, Azathioprine, Nintedanib, NAC |
| ⑨ Drug PD | GR transrepression (↓NF-κB), IMPDH inhibition, PDGFR/FGFR/VEGFR inhibition, ROS scavenging |
| ⑩ Clinical Endpoints | FVC%, DLCO%, KL-6, SP-D, BAL lymphocytosis, 6MWT, SGRQ, HRCT score, mortality |

---

## Model Files

| File | Description |
|------|-------------|
| `hp_qsp_model.dot` | Graphviz DOT source (100+ nodes, 10 clusters) |
| `hp_qsp_model.svg` | Vector mechanistic map |
| `hp_qsp_model.png` | Raster image (150 dpi) |
| `hp_mrgsolve_model.R` | mrgsolve ODE model + 7 treatment scenarios |
| `hp_shiny_app.R` | Interactive Shiny dashboard (8 tabs) |
| `hp_references.md` | 42 PubMed-linked references |

---

## ODE Compartments (22 States)

```
Antigen:   AG_lung
Innate:    M_M1, M_M2, Neutrophil
Cytokines: C_TNF, C_IL6, C_IL12, C_IFNg, C_TGFb, C_IL17, C_IL10
Adaptive:  T_Th1, T_Th17, T_Treg, Granuloma
Fibrosis:  Fibroblast, Myofib, Collagen, ROS
PK:        PDN_gut/cent, MPA_gut/cent, Nint_gut/cent
Function:  FVC, DLCO, KL6_serum
```

---

## Treatment Scenarios

| # | Scenario | Key Mechanism |
|---|----------|--------------|
| 1 | Untreated HP | Baseline disease progression |
| 2 | Antigen Avoidance (90%) | ↓ Antigen input → ↓ Th1/granuloma |
| 3 | Prednisolone 40mg/d | GR transrepression → ↓NF-κB/cytokines |
| 4 | MMF 1500mg BID | IMPDH inhibition → ↓T/B cell proliferation |
| 5 | Nintedanib 150mg BID | PDGFR/FGFR/VEGFR inhibition → ↓myofibroblast |
| 6 | PDN + MMF + partial avoidance | Combination immunosuppression + ↓antigen |
| 7 | Nintedanib + complete avoidance | Antifibrotic + primary prevention |

---

## Key PK Parameters

| Drug | F (%) | t½ | EC50 | Target |
|------|--------|-----|------|--------|
| Prednisolone | 82% | ~3.4 h | 150 ng/mL | GR / NF-κB |
| MMF (MPA) | 94% | ~17 h | 0.5 μg/mL | IMPDH |
| Nintedanib | 4.7% | ~10 h | 200 ng/mL | PDGFR/FGFR/VEGFR |

---

## Running the Model

```r
# Install dependencies
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "patchwork"))

# Run mrgsolve simulation (7 scenarios, 2-year)
source("hp_mrgsolve_model.R")

# Run interactive Shiny dashboard
install.packages(c("shiny", "shinydashboard", "plotly", "DT"))
shiny::runApp("hp_shiny_app.R")
```

```bash
# Render mechanistic map
neato -Tsvg hp_qsp_model.dot -o hp_qsp_model.svg
neato -Tpng -Gdpi=150 hp_qsp_model.dot -o hp_qsp_model.png
```

---

## Shiny Dashboard Tabs

1. **Patient Profile** — Disease overview, value boxes (FVC, DLCO, KL-6, Inflammation)
2. **PK Profiles** — Drug concentration-time curves for all agents
3. **Immune Response** — Th1/Th17/Treg dynamics, M1/M2 polarization, cytokine profiles
4. **Fibrosis & QSP** — Granuloma burden, fibroblast/myofibroblast, collagen deposition, ROS
5. **Clinical Endpoints** — FVC/DLCO trajectories, milestone table (Day 90/180/365/730)
6. **Scenario Comparison** — All 7 arms: FVC, collagen, inflammation side-by-side
7. **Biomarkers** — KL-6, inflammation index, biomarker summary table
8. **References** — Key publications

---

## Selected References

- Morisset J, et al. (2020) Mycophenolate vs Azathioprine in Fibrotic HP. *Lancet Respir Med* — [PMID 32169168](https://pubmed.ncbi.nlm.nih.gov/32169168/)
- Flaherty KR, et al. (2019) Nintedanib in Progressive Fibrosing ILD. *NEJM* — [PMID 31566307](https://pubmed.ncbi.nlm.nih.gov/31566307/)
- Raghu G, et al. (2020) ATS/JRS/ALAT HP Clinical Practice Guideline. *Am J Respir Crit Care Med* — [PMID 32706311](https://pubmed.ncbi.nlm.nih.gov/32706311/)
- Selman M, et al. (2012) HP Pathobiology. *Am J Respir Crit Care Med* — [PMID 22679007](https://pubmed.ncbi.nlm.nih.gov/22679007/)

---

*Model built by Claude Code Routine (CCR) | Date: 2026-06-28 | Category: Chronic Lung Disease*
