# Chronic Lymphocytic Leukemia (CLL) — QSP Model

[![Disease](https://img.shields.io/badge/Disease-CLL-red)]()
[![Drugs](https://img.shields.io/badge/Drugs-Ibrutinib%20%7C%20Venetoclax%20%7C%20Obinutuzumab-blue)]()
[![ODE Compartments](https://img.shields.io/badge/ODE%20Compartments-18-green)]()

## Overview

Chronic lymphocytic leukemia (CLL) is the most common adult leukemia in Western countries, characterized by accumulation of monoclonal CD19⁺CD5⁺CD23⁺ B cells in blood, bone marrow, and lymph nodes. This QSP model integrates:

- **BCR signaling** (BTK, PI3Kδ, PLCγ2, NF-κB)
- **Apoptotic machinery** (BCL-2 family: BCL-2, MCL-1, BIM, PUMA, BAX/BAK)
- **Tumor microenvironment** (CXCL12/CXCR4 BM retention, CXCL13/CXCR5 LN homing, CD40, BAFF)
- **Genetic risk stratification** (del(17p)/TP53, del(11q)/ATM, IGHV status, CLL-IPI)
- **Drug PK/PD**: BTK inhibitors (ibrutinib/acalabrutinib/zanubrutinib), venetoclax, obinutuzumab

## Mechanistic Map

[![CLL QSP Model](cll_qsp_model.png)](cll_qsp_model.svg)

*Click the image to open the full interactive SVG. 146 nodes, 10 subgraph clusters.*

| Cluster | Content |
|---------|---------|
| Normal B-cell development | HSC → CLP → Pro-B → Pre-B → Naive B → GC |
| CLL biology & staging | Rai 0–IV, Binet A–C, ALC, LDT, CD38/ZAP-70 |
| BCR signaling | LYN→SYK→BTK→PLCγ2→IP₃/DAG→NF-κB/PI3K/MAPK |
| Apoptosis (BCL-2 family) | BCL-2/MCL-1/BCL-XL vs BAX/BAK/BIM/PUMA/NOXA |
| Tumor microenvironment | BM/LN niche, CXCR4/CXCR5, NK, TAM, T-reg |
| Genetic risk | del(13q/11q/17p), Tri12, IGHV, TP53, SF3B1, CLL-IPI |
| BTK inhibitor PK/PD | Ibrutinib/acalabrutinib/zanubrutinib, BTK occupancy, lymphocytosis |
| BCL-2 inhibitor PK/PD | Venetoclax ramp-up, BCL-2 occupancy, TLS, MCL-1 resistance |
| Anti-CD20 PK/PD | Obinutuzumab TMDD, ADCC/CDC, NK activation |
| Clinical endpoints | IWCLL CR/PR/PD, MRD, PFS, OS, Richter transformation |

## mrgsolve ODE Model

**File:** `cll_mrgsolve_model.R`

| Feature | Detail |
|---------|--------|
| ODE compartments | 18 (≥15 required) |
| Drug PK | Ibrutinib (1-comp), Venetoclax (2-comp), Obinutuzumab (2-comp TMDD) |
| BTK occupancy | Covalent irreversible model, de novo synthesis (t½ ~2.9 d) |
| BCL-2 occupancy | Quasi-SS binding, Ki = 0.01 nM |
| Disease model | ALC (PB), BM infiltration (%), LN burden (%), logistic growth |
| Resistance | MCL-1 adaptive upregulation, NK activation |
| Scenarios | 6 treatment regimens (mono → triplet) |
| Calibration | RESONATE-2, CLL14, MURANO, SEQUOIA |

### 6 Treatment Scenarios

| # | Regimen | Key Reference |
|---|---------|---------------|
| 1 | Ibrutinib 420 mg QD | RESONATE-2 (Burger 2015 NEJM); ORR 86%, 2yr PFS 74% |
| 2 | Venetoclax 400 mg QD (ramp-up) | Roberts 2016 NEJM; ORR 79% (R/R) |
| 3 | Obinutuzumab ×6 cycles | CLL11 (Goede 2014 NEJM) |
| 4 | Venetoclax + Obinutuzumab | CLL14 (Fischer 2019 NEJM); uMRD blood 76%, 2yr PFS 88% |
| 5 | Ibrutinib + Venetoclax | CLARITY (Hillmen 2019 JCO); MRD-neg 53% |
| 6 | Triplet IB + VEN + OBI | Exploratory combination |

## Shiny App

**File:** `cll_shiny_app.R`

| Tab | Content |
|-----|---------|
| **Patient Profile** | CLL-IPI calculator, Rai/Binet staging, cytogenetics |
| **PK Profiles** | Ibrutinib, venetoclax, obinutuzumab concentration-time curves |
| **PD Biomarkers** | BTK occupancy, BCL-2 occupancy, CD20 occupancy, MCL-1 & NK |
| **Clinical Endpoints** | ALC over time, IWCLL response flags, tumor burden by compartment |
| **Scenario Comparison** | Head-to-head comparison of 6 regimens |
| **Genetic Risk** | del(17p)/TP53 and IGHV simulation, genomics table |

### Launch

```r
# Install required packages
install.packages(c("shiny","shinydashboard","mrgsolve","dplyr",
                   "ggplot2","plotly","tidyr","DT"))

# Run the app
shiny::runApp("cll_shiny_app.R")
```

## References

See [`cll_references.md`](cll_references.md) — 44 PubMed-linked references covering:
- CLL pathophysiology & staging
- BCR signaling & BTK biology
- Apoptosis & BCL-2 family
- Genetic risk factors (CLL-IPI, del17p, IGHV)
- BTK inhibitor trials (RESONATE-2, ELEVATE-TN, SEQUOIA, ALPINE)
- Venetoclax trials (MURANO, CLL14, CLARITY)
- Anti-CD20 mAb PK/PD
- MRD & response assessment
- QSP/PKPD modeling methods

## Key Biology Highlights

### Why CLL Cells Survive
CLL cells exploit multiple pro-survival signals:
1. **Autonomous BCR signaling** (IGHV-unmutated CLL) → NF-κB, PI3K
2. **BCL-2 overexpression** (del13q → miR-15a/16-1 loss) → apoptosis resistance
3. **Microenvironment support** (CXCL12/CXCR4 BM retention, BAFF, CD40L, IL-4)

### Drug Mechanisms
| Drug | Target | Key PD Effect |
|------|--------|---------------|
| Ibrutinib | BTK Cys481 (covalent) | Blocks NF-κB/PI3K; redistributes CLL from BM/LN → PB (lymphocytosis) |
| Venetoclax | BCL-2 BH3 groove | Frees BIM/BAX/BAK → MOMP → caspase cascade; MCL-1 compensatory rise |
| Obinutuzumab | CD20 (Type II) | ADCC (NK), direct cell death (PCD), less CDC vs rituximab |

### Resistance Mechanisms
- **BTK C481S mutation**: prevents covalent binding; switch to pirtobrutinib (non-covalent)
- **PLCγ2 gain-of-function**: bypasses BTK blockade downstream
- **BCL-2 G101V**: reduces venetoclax binding affinity
- **MCL-1 upregulation**: compensates for BCL-2 inhibition under venetoclax

## Model Parameters (Key Values)

| Parameter | Value | Source |
|-----------|-------|--------|
| Ibrutinib t½ | ~7 h (CL=980 L/h) | Pharmacyclics population PK |
| BTK protein t½ | ~2.9 d (kdeg=0.010 h⁻¹) | Honigberg 2010 |
| Venetoclax t½ | ~26 h (V1=250 L, CL=65 L/h) | AbbVie popPK |
| BCL-2 Ki venetoclax | 0.01 nM | Souers 2013 Nat Med |
| Obinutuzumab t½ | ~28 d | Mössner 2010 |
| CLL doubling time | ~9.7 mo (kprol=0.003 h⁻¹) | Typical newly-diagnosed |
