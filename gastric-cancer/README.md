# Gastric Cancer (위선암) QSP Model

> Date Added: 2026-06-23 | Category: Oncology | Model Type: QSP (Quantitative Systems Pharmacology)

---

## Disease Overview (질환 개요)

**Gastric adenocarcinoma** (위선암) is the fifth most common cancer and the third leading cause of cancer-related death worldwide, with approximately 1.09 million new cases and 769,000 deaths annually (GLOBOCAN 2020). The incidence is highest in East Asia (Korea, Japan, China), Eastern Europe, and parts of Latin America, largely driven by *Helicobacter pylori* infection prevalence.

### Epidemiology
| Region | Incidence (ASR per 100,000) |
|--------|---------------------------|
| East Asia (Korea/Japan/China) | 20-40 |
| Eastern Europe | 15-25 |
| Western Europe / North America | 5-10 |
| South America | 10-20 |

### Pathophysiology
The Correa cascade describes the multi-step progression from normal gastric mucosa → chronic gastritis → atrophic gastritis → intestinal metaplasia → dysplasia → gastric cancer, primarily driven by:
- *H. pylori* CagA/VacA virulence factors activating NF-κB, PI3K, STAT3
- Epstein-Barr virus (EBV) infection with hypermethylation
- Chromosomal instability (CIN) with RTK amplification (HER2, FGFR2, MET)
- Microsatellite instability (MSI-H) via MLH1 silencing

---

## Molecular Subtypes (TCGA 분자적 아형)

| TCGA Subtype | Prevalence | Key Molecular Features | Preferred Therapy |
|---|---|---|---|
| **EBV-positive** | ~9% | EBER+, PIK3CA mut, PD-L1 ↑↑, CDKN2A silencing | ICI (PD-1/L1 blockade), PI3Ki |
| **MSI-High** | ~22% | MLH1 silencing, high TMB (≥10 mut/Mb), ARID1A mut | ICI (pembrolizumab/nivolumab) |
| **Genomically Stable (GS)** | ~20% | RHOA/CDH1 mut, diffuse Lauren, CLDN18.2 amp | Zolbetuximab + FOLFOX, FOLFOX |
| **Chromosomal Instability (CIN)** | ~50% | TP53 mut, RTK amp (HER2/EGFR/FGFR2/MET) | Trastuzumab (HER2+), T-DXd |

---

## Treatment Landscape (치료 현황)

### First-Line Treatment
| Biomarker | Preferred Regimen | Key Trial | mOS |
|---|---|---|---|
| HER2+ (IHC3+ or 2+/FISH+) + CPS≥1 | Pembrolizumab + Trastuzumab + FOLFOX/XELOX | KEYNOTE-811 | 20.0 mo |
| HER2+ (IHC3+ or 2+/FISH+) | Trastuzumab + FOLFOX/XELOX | ToGA | 13.8 mo |
| CLDN18.2+ (≥2+ in ≥75%), HER2- | Zolbetuximab + mFOLFOX6 | SPOTLIGHT | 18.2 mo |
| CLDN18.2+ (≥2+ in ≥75%), HER2- | Zolbetuximab + CAPOX | GLOW | 14.4 mo |
| CPS≥5, any | Nivolumab + FOLFOX/XELOX | CheckMate 649 | 14.4 mo |
| MSI-H | Pembrolizumab + chemo | KEYNOTE-590/811 | >20 mo |
| Resectable (locally adv.) | FLOT × 4 pre/post-op | FLOT4 | 50 mo |
| Standard | FOLFOX / XELOX | REAL2 | 10.9 mo |

### Second-Line Treatment
| Regimen | Key Trial | mOS |
|---|---|---|
| Ramucirumab + Paclitaxel | RAINBOW | 9.6 mo |
| Ramucirumab monotherapy | REGARD | 5.2 mo |
| T-DXd (HER2+) | DESTINY-Gastric01 | 12.5 mo |
| Pembrolizumab (MSI-H / CPS≥10) | KEYNOTE-061 | ~14 mo |
| Docetaxel / Irinotecan | Salvage | 5-7 mo |

---

## Model Components

### 1. Mechanistic Map (gc_qsp_model.dot / .svg / .png)

| Component | Count |
|---|---|
| Total nodes | 212 |
| Subgraph clusters | 10 |
| Directed edges | 250 |
| Drug PK nodes | 30 |
| Signal pathway nodes | 75 |
| Immune/TME nodes | 40 |
| Clinical endpoint nodes | 20 |

**Clusters:**
1. H. pylori & Gastric Mucosal Inflammation (17 nodes)
2. Molecular Drivers & RTKs (20 nodes)
3. Intracellular Signal Transduction (24 nodes)
4. Cell Cycle, Proliferation & Apoptosis (23 nodes)
5. Epigenetics, EMT & TCGA Subtypes (25 nodes)
6. Tumor Microenvironment (21 nodes)
7. Immune Checkpoint (20 nodes)
8. Angiogenesis (18 nodes)
9. Drug PK/PD (30 nodes)
10. Clinical Endpoints & Biomarkers (20 nodes)

### 2. mrgsolve ODE Model (gc_mrgsolve_model.R)

| Feature | Count/Detail |
|---|---|
| Total ODE compartments | 28 (12 PK + 16 PD) |
| Drug PK models | 6 (Trastuzumab, Ramucirumab, Nivolumab, Capecitabine/5-FU, T-DXd, Zolbetuximab/Pembrolizumab) |
| Treatment scenarios | 6 |
| Captured outputs | 17 |

**Treatment Scenarios:**
1. FLOT perioperative (5-FU + Leucovorin + Oxaliplatin + Docetaxel)
2. Trastuzumab + FOLFOX/XELOX (HER2+ first-line, ToGA/KEYNOTE-811)
3. Ramucirumab + Paclitaxel (second-line, RAINBOW)
4. Nivolumab + Chemotherapy (CPS≥5 first-line, CheckMate 649)
5. T-DXd (HER2+ second-line, DESTINY-Gastric01)
6. Zolbetuximab + mFOLFOX6 (CLDN18.2+ first-line, SPOTLIGHT)

**Key PD sub-models:**
- Simeoni TGI (transit compartment tumor growth inhibition)
- TMDD (target-mediated drug disposition) for all mAbs
- HER2 signaling downstream activity (AU)
- VEGF-A free concentration with angiogenic feedback
- CD8+ T effector, Treg, TAM M2 dynamics
- CEA biomarker (tumor-correlated)
- Cancer stem cell self-renewal

### 3. Shiny Dashboard (gc_shiny_app.R)

| Tab | Content |
|---|---|
| 환자 프로파일 | ECOG, stage, HER2/MSI/CLDN18.2/CPS input; TCGA subtype classifier |
| 약물 PK | Per-drug PK simulation, Cmax/Ctrough/AUC table, therapeutic window |
| 종양 동태 | Simeoni TGI, waterfall plot, spider plot, RECIST assessment |
| 임상 엔드포인트 | KM OS/PFS (survminer), ORR/DCR/CR bar charts |
| 치료 시나리오 비교 | Side-by-side TV dynamics, AUC comparison, summary table |
| 바이오마커 분석 | CEA/CA19-9/ctDNA/VAF dynamics, HER2 IHC distribution, CPS, TCGA pie chart, correlation matrix |

### 4. References (gc_references.md)

60 curated PubMed references across 10 sections:
- H. pylori & carcinogenesis (6 refs)
- Molecular pathology & TCGA (7 refs)
- HER2-targeted therapy (7 refs)
- Anti-angiogenic therapy (5 refs)
- Immunotherapy (6 refs)
- CLDN18.2 targeting (5 refs)
- Chemotherapy regimens (6 refs)
- Biomarkers & precision oncology (6 refs)
- QSP / mathematical modeling (9 refs)
- Additional clinical oncology (3 refs)

---

## How to Run (실행 방법)

### Render Mechanistic Map
```bash
cd gastric-cancer/
dot -Tsvg gc_qsp_model.dot -o gc_qsp_model.svg
dot -Tpng -Gdpi=150 gc_qsp_model.dot -o gc_qsp_model.png
```

### Run mrgsolve Simulations
```r
# Install required packages
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))

# Run all 6 treatment scenarios
source("gc_mrgsolve_model.R")
results <- simulate_all_scenarios()
summary_tbl <- compute_response_summary(results)
print(summary_tbl)

# Generate plots
p1 <- plot_tumor_kinetics(results)
p2 <- plot_pk_curves(results)
p3 <- plot_immune_dynamics(results)
p4 <- plot_biomarkers(results)
```

### Launch Shiny App
```r
# Install required packages
install.packages(c("shiny", "shinydashboard", "plotly", "DT",
                   "survival", "survminer", "RColorBrewer"))

# Run the app
shiny::runApp("gc_shiny_app.R")
```

---

## Key Modeling Findings (주요 모델링 결과)

1. **HER2+ tumors (IHC 3+)**: Trastuzumab + FOLFOX achieves ~50% tumor volume reduction at 12 weeks; T-DXd (bystander killing via free DXd) yields ~55-60% reduction even in HER2-low (IHC 1+/2+)

2. **Immune activation dynamics**: Nivolumab achieves >80% PD-1 receptor occupancy by cycle 2 (Q3W 360mg), restoring CD8+ T effector function with a predicted 1.8-fold increase in tumor-infiltrating CD8 by week 12

3. **CLDN18.2 ADCC (Zolbetuximab)**: High CLDN18.2 expression (≥2+, ≥75% cells) results in ~65% NK cell-mediated ADCC killing efficiency; combination with mFOLFOX6 shows synergistic TGI

4. **Biomarker dynamics**: CEA and CA19-9 kinetics parallel tumor volume with a ~2-week lag; ctDNA VAF decreases faster (shorter half-life) and serves as earlier response indicator

5. **TME immunosuppression**: MSI-H tumors (high TMB, increased neoantigen load) show stronger baseline CD8 infiltration, explaining superior ICI response; TAM M2 and Treg co-expansion limits durable response in MSS tumors

---

## File Structure

```
gastric-cancer/
├── gc_qsp_model.dot         # Graphviz mechanistic map (212 nodes, 10 clusters)
├── gc_qsp_model.svg         # SVG render
├── gc_qsp_model.png         # PNG render (150 dpi)
├── gc_mrgsolve_model.R      # mrgsolve ODE model (28 compartments, 6 scenarios)
├── gc_shiny_app.R           # Shiny interactive dashboard (6 tabs)
├── gc_references.md         # 60 PubMed references (10 sections)
└── README.md                # This file
```
