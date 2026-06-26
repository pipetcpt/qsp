# Hepatocellular Carcinoma (HCC) — QSP Model

[![Disease](https://img.shields.io/badge/Disease-Hepatocellular%20Carcinoma-red)](.)
[![Category](https://img.shields.io/badge/Category-Oncology%20%7C%20Liver-orange)](.)
[![Drugs](https://img.shields.io/badge/Drugs-Sorafenib%20%7C%20Lenvatinib%20%7C%20Atezo%2BBeva%20%7C%20Regorafenib-blue)](.)
[![ODEs](https://img.shields.io/badge/ODE%20Compartments-20-green)](hcc_mrgsolve_model.R)
[![Refs](https://img.shields.io/badge/References-43-purple)](hcc_references.md)

---

## Disease Overview

**Hepatocellular carcinoma (HCC)** is the most common primary liver malignancy, accounting for ~90% of primary liver cancers. It arises almost exclusively in the setting of chronic liver disease — hepatitis B/C infection, alcoholic cirrhosis, or NAFLD/NASH — with cirrhosis being the common soil from which HCC emerges. HCC is the **4th leading cause of cancer death worldwide** (~800,000 deaths/year).

The disease is staged using the **Barcelona Clinic Liver Cancer (BCLC)** system, which integrates tumor extent, vascular invasion, liver functional reserve (Child-Pugh/MELD), and performance status:

| BCLC Stage | Extent | Treatment | Median OS |
|------------|--------|-----------|-----------|
| 0 / A | Very early / Early | Resection, transplant, ablation | 5+ years |
| B | Intermediate | TACE | ~2.5 years |
| C | Advanced | Systemic therapy | 13–19 months |
| D | Terminal | Best supportive care | <3 months |

---

## QSP Model Overview

This model integrates **9 biological clusters** and covers the complete treatment landscape from molecular oncogenesis to clinical endpoints.

### Mechanistic Map Components

| Cluster | Key Nodes | Biological Role |
|---------|-----------|-----------------|
| ① Risk Factors & Etiology | HBV, HCV, NAFLD, alcohol, cirrhosis, ROS | Disease initiation |
| ② Oncogenic Events | TP53, CTNNB1, TERT, HBx, HCV NS5A | Transformation drive |
| ③ Intracellular Signaling | RAS/RAF/MEK/ERK, PI3K/AKT/mTOR, Wnt/β-catenin, STAT3, NF-κB | Proliferation/survival |
| ④ Angiogenesis & Hypoxia | VEGF-A/C, VEGFR1-3, PDGF-BB, HIF-1α/2α, TIE2 | Neovascularization |
| ⑤ Tumor Biology & EMT | Proliferation, apoptosis, CSC, SNAIL, EMT, MMP-2/9 | Invasion/metastasis |
| ⑥ Tumor Immune Microenvironment | CD8+ CTL, Tregs, TAMs, MDSCs, PD-1/PD-L1, CTLA-4 | Immune evasion |
| ⑦ Drug PK | Sorafenib, Lenvatinib, Atezolizumab, Bevacizumab, Nivolumab, Regorafenib, Cabozantinib, Ramucirumab | PK compartments |
| ⑧ Drug PD Targets | VEGFR2 inhibition, RAF inhibition, PD-L1 blockade, VEGF neutralization | Drug mechanisms |
| ⑨ Clinical Endpoints | BCLC, Child-Pugh, AFP, mRECIST, OS, PFS, TTP, ORR | Outcomes |

### Mechanistic Map

[![HCC QSP Mechanistic Map](hcc_qsp_model.png)](hcc_qsp_model.svg)

*Click image to open interactive SVG. 120+ nodes, 9 clusters, drug PK/PD wiring.*

---

## mrgsolve ODE Model

**File:** [`hcc_mrgsolve_model.R`](hcc_mrgsolve_model.R)

### Model Architecture (20 ODE compartments)

| Compartment Block | Compartments | Description |
|-------------------|-------------|-------------|
| Sorafenib PK | GUT_S, CENTRAL_S, PERIPH_S, MET_S | 2-cmt oral + active metabolite |
| Lenvatinib PK | GUT_L, CENTRAL_L, PERIPH_L | 2-cmt oral |
| Atezolizumab PK | CENTRAL_A, PERIPH_A | 2-cmt IV mAb |
| Bevacizumab PK | CENTRAL_B, PERIPH_B | 2-cmt IV mAb |
| Regorafenib PK | GUT_R, CENTRAL_R | 2-cmt oral (2nd-line) |
| Tumor Dynamics | TUMOR | Gompertz growth + drug kill + immune kill |
| Immune Cells | CD8T, TREG | CTL & regulatory T cell dynamics |
| Microenvironment | ANGIO, VEGF_FREE | Angiogenic drive + free VEGF-A |
| Biomarkers | AFP, LF | AFP ng/mL + liver function reserve |

### Key PD Model Features

- **Anti-angiogenic TGI**: Emax model for VEGFR2 inhibition (sorafenib/lenvatinib) + VEGF-A neutralization (bevacizumab); combined using Bliss independence
- **MAPK TGI**: Emax model for RAF kinase inhibition (sorafenib/regorafenib)
- **Immune checkpoint blockade**: CD8 T cell reinvigoration proportional to atezolizumab PD-L1 occupancy
- **Immune-tumor interactions**: CD8-mediated tumor lysis, T cell exhaustion by tumor, Treg-mediated suppression
- **Liver function decay**: Proportional to tumor burden, reversible

### 5 Treatment Scenarios

| # | Regimen | Reference Trial | Median OS |
|---|---------|----------------|-----------|
| 1 | **Sorafenib 400mg BID** | SHARP (Llovet 2008) | 10.7 mo |
| 2 | **Lenvatinib 8/12mg QD** | REFLECT (Kudo 2018) | 13.6 mo |
| 3 | **Atezo 1200mg + Beva 15mg/kg q3w** | IMbrave150 (Finn 2020) | 19.2 mo |
| 4 | **Regorafenib 160mg QD (3wk/1wk)** | RESORCE (Bruix 2017) | 10.6 mo |
| 5 | **Best Supportive Care (BSC)** | Historical | ~7 mo |

---

## Shiny Dashboard

**File:** [`hcc_shiny_app.R`](hcc_shiny_app.R)

### 6 Interactive Tabs

| Tab | Content |
|-----|---------|
| **① Patient Profile** | BCLC stage, Child-Pugh, AFP baseline, BCLC algorithm guide |
| **② Drug PK** | Time-concentration profiles for all active drugs |
| **③ Tumor Dynamics** | RECIST waterfall, TGI%, best response classification |
| **④ PD & Biomarkers** | AFP trajectory, liver function reserve, VEGF-A, angiogenesis state |
| **⑤ Scenario Comparison** | All 5 regimens overlaid + 6-month efficacy summary table |
| **⑥ Immune & Safety** | CD8/Treg dynamics, ICB effect, adverse event profile by drug |

**Run:** `shiny::runApp("hcc_shiny_app.R")`

---

## Key References

| Drug/Topic | Trial | Reference |
|------------|-------|-----------|
| Sorafenib | SHARP | [Llovet 2008 NEJM](https://pubmed.ncbi.nlm.nih.gov/18650514/) |
| Lenvatinib | REFLECT | [Kudo 2018 Lancet](https://pubmed.ncbi.nlm.nih.gov/29433850/) |
| Atezo+Beva | IMbrave150 | [Finn 2020 NEJM](https://pubmed.ncbi.nlm.nih.gov/32402160/) |
| Regorafenib | RESORCE | [Bruix 2017 Lancet](https://pubmed.ncbi.nlm.nih.gov/27932229/) |
| Cabozantinib | CELESTIAL | [Abou-Alfa 2018 NEJM](https://pubmed.ncbi.nlm.nih.gov/29972759/) |
| Ramucirumab | REACH-2 | [Zhu 2019 Lancet Oncol](https://pubmed.ncbi.nlm.nih.gov/30665869/) |
| mRECIST | — | [Lencioni 2010 Semin Liver Dis](https://pubmed.ncbi.nlm.nih.gov/20175033/) |
| BCLC Staging | — | [Llovet 1999 Semin Liver Dis](https://pubmed.ncbi.nlm.nih.gov/10518312/) |
| Guidelines | AASLD 2018 | [Marrero 2018 Hepatology](https://pubmed.ncbi.nlm.nih.gov/29624699/) |

Full reference list (43 sources): [`hcc_references.md`](hcc_references.md)

---

## Files in This Directory

| File | Description |
|------|-------------|
| `hcc_qsp_model.dot` | Graphviz source — full mechanistic map |
| `hcc_qsp_model.svg` | Interactive SVG (click to open) |
| `hcc_qsp_model.png` | PNG thumbnail (150 dpi) |
| `hcc_mrgsolve_model.R` | ODE model (20 compartments, 5 scenarios) |
| `hcc_shiny_app.R` | Interactive Shiny dashboard (6 tabs) |
| `hcc_references.md` | 43 PubMed references |
| `README.md` | This file |

---

*Generated by Claude Code Routine (CCR) · 2026-06-23 · QSP Disease Model Library*
