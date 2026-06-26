# Breast Cancer QSP Model
# 유방암 정량적 시스템 약리학 모델

---

## 질환 개요 (Disease Overview)

유방암은 전 세계적으로 여성에서 가장 흔한 악성 종양으로, 2020년 기준 연간 약 230만 명의 신규 환자가 발생하며 전체 암 발생의 11.7%를 차지합니다. 분자 아형에 따라 크게 네 가지로 분류되며, 에스트로겐 수용체 양성(ER+)/HER2 음성이 전체의 약 65~70%, HER2 양성(HER2+)이 15~20%, 삼중음성 유방암(TNBC)이 약 15%를 차지합니다. 각 아형은 상이한 병태생리학적 기전, 치료 반응성 및 예후를 보이며, 이에 따른 맞춤형 치료 전략이 필수적입니다. 한국에서는 서구권과 유사하게 유방암 발생률이 꾸준히 증가하고 있으며 40~50대 여성에서의 유병률이 높습니다.

치료 패러다임은 분자 아형별로 크게 다릅니다. ER+ 유방암의 표준 치료는 아로마타제 억제제(AI) 또는 타목시펜을 기반으로 하는 내분비 요법이며, CDK4/6 억제제(팔보시클립, 리보시클립, 아베마시클립)의 병용으로 무진행 생존 기간이 유의하게 연장되었습니다. HER2 양성 유방암에서는 트라스투주맙, 퍼투주맙을 포함한 이중 HER2 차단 요법 및 항체-약물 접합체(T-DM1, T-DXd)가 생존 이득을 크게 향상시켰습니다. TNBC는 아직 표적 치료제가 제한적이나, 면역 관문 억제제(펨브롤리주맙), PARP 억제제(올라파립) 및 항체-약물 접합체의 등장으로 치료 옵션이 확대되고 있습니다.

Breast cancer is the most common malignancy in women worldwide, with an estimated 2.3 million new cases diagnosed in 2020, accounting for approximately 11.7% of all cancer diagnoses globally. It is classified into distinct molecular subtypes based on receptor expression: hormone receptor-positive (HR+)/HER2-negative (~65–70%), HER2-positive (~15–20%), and triple-negative breast cancer (TNBC, ~15%). Each subtype harbors unique oncogenic drivers, therapeutic sensitivities, and prognostic profiles. The HR+ subtype is driven by estrogen receptor (ER) signaling, HER2+ by overexpression/amplification of the ERBB2 receptor tyrosine kinase, and TNBC lacks these established targets and is associated with the poorest prognosis.

The treatment landscape for breast cancer has been transformed over the past two decades. In HR+/HER2- disease, CDK4/6 inhibitors combined with aromatase inhibitors or fulvestrant have become first-line standard of care in metastatic settings, significantly extending progression-free survival (mPFS 25–28 months vs. ~14–16 months for endocrine monotherapy). HER2+ disease benefits from highly effective anti-HER2 therapies including dual blockade with trastuzumab and pertuzumab, and next-generation antibody-drug conjugates (ADCs) such as T-DM1 and T-DXd. Despite these advances, resistance remains a major clinical challenge: ESR1 mutations develop under aromatase inhibitor pressure, CDK4/6 inhibitor resistance emerges via RB loss and cyclin E amplification, and HER2+ resistance is mediated by PI3K/AKT pathway activation. The QSP model developed here provides a mechanistic framework to simulate these dynamics, quantify treatment effects, and explore combination strategies.

---

## 모델 구조 (Model Structure)

| Component | Details |
|-----------|---------|
| Mechanistic Map Nodes | 110+ nodes, 9 subgraph clusters |
| ODE Compartments | 22 compartments |
| Drug Coverage | 14 drugs (6 drug classes) |
| Treatment Scenarios | 6 clinical regimens |
| Simulation Duration | Up to 1 year (8,760 hours) |
| Key Pathways Modeled | ER/PI3K/AKT/mTOR, HER2/MAPK, CDK4/6-RB, DNA damage/PARP, PD-1/PD-L1 |
| Biomarkers Simulated | Ki-67, CA15-3, RB phosphorylation, CD8+ TILs, PD-L1, ctDNA (ESR1) |
| Resistance Modules | ESR1 mutation, RB loss, PI3K activation, CDK4/6i bypass |

---

## 포함된 약물 (Drug Coverage)

| Drug | Class | Target | Key Trial |
|------|-------|--------|-----------|
| Palbociclib | CDK4/6 inhibitor | CDK4/CDK6 | PALOMA-2 |
| Ribociclib | CDK4/6 inhibitor | CDK4/CDK6 | MONALEESA-2 |
| Abemaciclib | CDK4/6 inhibitor | CDK4/CDK6 | MONARCH-3 |
| Letrozole | Aromatase inhibitor | CYP19A1 | PALOMA-2 |
| Anastrozole | Aromatase inhibitor | CYP19A1 | ATAC |
| Tamoxifen | SERM | ERα | EBCTCG meta-analysis |
| Trastuzumab | Anti-HER2 antibody | HER2/ERBB2 | CLEOPATRA |
| Pertuzumab | Anti-HER2 antibody | HER2 dimerization domain | CLEOPATRA |
| T-DM1 | ADC (HER2-targeted) | HER2 + tubulin | EMILIA |
| T-DXd | ADC (HER2-targeted) | HER2 + Topo I | DESTINY-Breast03 |
| Pembrolizumab | Anti-PD-1 mAb | PD-1 | KEYNOTE-522 |
| Olaparib | PARP inhibitor | PARP1/2 | OlympiAD |
| Alpelisib | PI3K inhibitor | PIK3Cα (PI3Kα) | SOLAR-1 |
| Everolimus | mTOR inhibitor | mTORC1 | BOLERO-2 |

---

## 주요 임상시험 시뮬레이션 (Key Clinical Trials Simulated)

| Trial | Regimen | Subtype | Key Result | Model Calibration Target |
|-------|---------|---------|------------|--------------------------|
| PALOMA-2 | Palbociclib + Letrozole | ER+/HER2- | mPFS 27.6 vs 14.5 mo | CDK4/6 inhibition depth, Ki-67 ≥50% suppression |
| MONALEESA-2 | Ribociclib + Letrozole | ER+/HER2- | mPFS 25.3 vs 16.0 mo | CDK4/6 + ER signaling synergy |
| MONARCH-3 | Abemaciclib + AI | ER+/HER2- | mPFS 28.2 vs 14.8 mo | CDK4/6 sustained inhibition (continuous dosing) |
| CLEOPATRA | Trastuzumab + Pertuzumab + Docetaxel | HER2+ | mPFS 18.7 mo, HR=0.62 | HER2 dimerization blockade, MAPK suppression |
| KEYNOTE-522 | Pembrolizumab + Chemo (neoadjuvant) | TNBC | pCR 64.8% vs 51.2% | CD8 expansion, PD-1 checkpoint release |
| OlympiAD | Olaparib monotherapy | BRCAm, HER2- | mPFS 7.0 vs 4.2 mo | PARP inhibition, synthetic lethality in BRCA-null |

---

## 모델 실행 방법 (How to Run)

### Prerequisites

```r
install.packages(c("mrgsolve", "ggplot2", "dplyr", "tidyr", "shiny",
                   "shinydashboard", "plotly"))
```

### mrgsolve ODE Model

```r
source("bc_mrgsolve_model.R")
# Simulates all 6 treatment scenarios
# Outputs: tumor volume, Ki-67, PK profiles, biomarker dynamics
```

### Shiny Interactive App

```r
shiny::runApp("bc_shiny_app.R")
# Opens interactive dashboard with 6 tabs:
# Tab 1: Patient Profile | Tab 2: PK | Tab 3: PD | Tab 4: Clinical Endpoints
# Tab 5: Scenario Comparison | Tab 6: Biomarkers
```

### Generate Mechanistic Map

```bash
# SVG (vector, for web/report embedding)
dot -Tsvg bc_qsp_model.dot -o bc_qsp_model.svg

# PNG (150 dpi raster, for README thumbnail)
dot -Tpng -Gdpi=150 bc_qsp_model.dot -o bc_qsp_model.png
```

---

## 파일 설명 (File Descriptions)

| File | Description |
|------|-------------|
| `bc_qsp_model.dot` | Graphviz DOT mechanistic map (110+ nodes, 9 clusters) |
| `bc_qsp_model.svg` | Rendered SVG of mechanistic map (vector format) |
| `bc_qsp_model.png` | Rendered PNG of mechanistic map (150 dpi) |
| `bc_mrgsolve_model.R` | mrgsolve ODE model with 22 compartments and 6 treatment scenarios |
| `bc_shiny_app.R` | Interactive Shiny dashboard with 6 tabs |
| `bc_references.md` | 50 PubMed-linked references organized by 9 sections |
| `README.md` | This file — model overview and usage guide |

---

## 주요 파라미터 (Key Parameters)

| Parameter | Symbol | Value | Unit | Source |
|-----------|--------|-------|------|--------|
| Tumor proliferation rate | kprol | 0.0008 | 1/hr | Calibrated to PALOMA-2 control arm |
| Tumor death rate | kdeath | 0.0002 | 1/hr | Literature (Simeoni model) |
| Tumor carrying capacity | Kmax | 1000 | relative units | Estimated |
| E2 EC50 (proliferation stimulation) | EC50_E2 | 50 | pmol/L | ER biology literature |
| Letrozole Emax (E2 suppression) | Emax_AI | 0.98 | fraction | ATAC/BIG 1-98 PD data |
| CDK4/6 inhibitor Emax | Emax_CDK | 0.85 | fraction | PALOMA-2 Ki-67 PD data |
| CDK4/6 inhibitor EC50 | EC50_CDK | 100 | ng/mL | Population PK/PD modeling |
| Palbociclib clearance | CL_palbo | 63 | L/hr | Population PK (Friberg) |
| Palbociclib volume of distribution | Vd_palbo | 2583 | L | Population PK (Friberg) |
| Trastuzumab clearance | CL_tras | 0.225 | L/day | Population PK (Bruno) |
| Trastuzumab Vd (central) | Vc_tras | 3.1 | L | Population PK (Bruno) |
| HER2 signaling Emax | Emax_HER2 | 0.70 | fraction | HER2 CLEOPATRA PD data |
| HER2 EC50 (trastuzumab) | EC50_HER2 | 50 | mg/L | PK/PD estimate |
| CD8+ T cell kill rate | kCD8_kill | 0.002 | 1/hr | Immunology model (Jiang) |
| PD-1 inhibition Emax | Emax_PD1 | 0.80 | fraction | KEYNOTE-522 pCR calibration |
| Olaparib clearance | CL_olap | 8.6 | L/hr | Population PK (Karlsson) |
| Olaparib Emax (PARP inhibition) | Emax_PARP | 0.95 | fraction | OlympiAD biomarker data |
| ESR1 mutation rate | mu_ESR1 | 1e-6 | 1/cell/hr | Estimated from ctDNA kinetics |
| RB loss rate (CDK4/6i resistance) | mu_RBlos | 5e-7 | 1/cell/hr | Turner 2019 (PALOMA-3 resistance) |
| PI3K activation rate (resistance) | kPI3K | 0.0003 | 1/hr | Juric 2015 (SOLAR-1) |

---

## 예상 시뮬레이션 출력 (Expected Simulation Outputs)

1. **Tumor Volume Trajectories** — All 6 regimens plotted as % change from baseline over 52 weeks, with uncertainty bands calibrated to trial median PFS times.
2. **Ki-67 Proliferation Index** — Dynamic change from baseline (%) at C1D14 and C2D1, calibrated to PALOMA-2 and MONARCH-3 biopsy data.
3. **CDK4/6 Inhibition Depth** — RB phosphorylation inhibition (%) over the dosing cycle; sustained inhibition shown for abemaciclib (continuous) vs. palbociclib/ribociclib (intermittent).
4. **CA15-3 Tumor Marker** — Serum biomarker response kinetics correlating with tumor burden; waterfall plots at Week 12.
5. **PK Profiles** — Drug concentration-time curves showing Cmax, Cmin, and AUC for each agent; steady-state achieved within 5–8 days for oral CDK4/6i.
6. **Immune Cell Dynamics** — CD8+ effector T cell expansion, regulatory T cell (Treg) suppression, and PD-L1 upregulation trajectories for TNBC/pembrolizumab simulation.
7. **PFS Curves** — Kaplan-Meier-like time-to-event curves for all scenarios, reproduced from model output; hazard ratios vs. control computed from tumor growth inhibition model.

---

## 한계 및 가정 (Limitations & Assumptions)

- **Spatial homogeneity**: Tumor spatial heterogeneity is not modeled; the well-mixed compartment assumption is used throughout (no intratumoral gradients in drug penetration or oxygen).
- **Binary resistance states**: Acquired resistance is modeled as deterministic transitions between sensitive and resistant subpopulations, not as continuous stochastic evolutionary processes.
- **Simplified immune microenvironment**: The immune module uses a 3-cell-type model (CD8+ effectors, Tregs, tumor cells) without dendritic cell, macrophage, or NK cell contributions.
- **DDI not modeled**: Drug-drug pharmacokinetic interactions (e.g., CYP3A4 induction/inhibition affecting CDK4/6i exposure) are not explicitly represented for combination regimens.
- **Cardiac toxicity**: LVEF decline (trastuzumab cardiotoxicity, ribociclib QTc prolongation) is modeled empirically as a dose-dependent function rather than mechanistically via cardiomyocyte biology.
- **Single tumor site**: Metastatic spread to multiple organ sites (bone, lung, liver, brain) is represented by a single virtual tumor compartment; site-specific PK differences are approximated by scaling parameters.
- **Population homogeneity**: The base model uses a single representative patient; population variability (BSV) must be introduced via IIV parameter distributions for population-level predictions.

---

## 개발 정보 (Development Information)

| Item | Details |
|------|---------|
| 생성일 (Created) | 2026-06-21 |
| 버전 (Version) | 1.0.0 |
| 작성자 (Author) | Claude Code Routine (CCR) — Automated QSP Library |
| 소프트웨어 (Software) | R ≥ 4.0, mrgsolve ≥ 1.0, Graphviz ≥ 2.40, Shiny ≥ 1.7 |
| 모델 프레임워크 (Framework) | ODE-based, deterministic; Simeoni-type tumor growth backbone |
| 검증 방법 (Validation) | Calibrated against PALOMA-2, MONALEESA-2, MONARCH-3, CLEOPATRA, KEYNOTE-522, OlympiAD primary endpoints |
| 참고문헌 (References) | See `bc_references.md` (50 PubMed-linked citations, 9 sections) |

---

*This model is part of the QSP Disease Model Library (CCR). See the main [README](../README.md) for the full library index.*
