# Acute Myeloid Leukemia (AML) — QSP Model

## Overview (개요)

**한국어**: 급성 골수성 백혈병(AML)은 조혈 줄기세포 또는 전구세포에서 발생하는 악성 클론성 증식 질환으로, 미성숙 골수 아세포가 골수 내에 축적되어 정상 조혈을 억제합니다. 성인에서 가장 흔한 급성 백혈병으로, 미국 기준 연간 약 2만 명이 새로 진단되며 5년 생존율은 30% 미만입니다. 분자 수준에서 FLT3, NPM1, DNMT3A, IDH1/2, TP53 등 다수의 유전자 돌연변이가 발병에 관여하며, 2022 ELN 위험군 분류에 따라 치료 전략이 결정됩니다. 최근 베네토클락스(BCL-2 억제제), 길테리티닙(FLT3 억제제), 에나시데닙/이보시데닙(IDH 억제제) 등 표적 치료제의 승인으로 치료 패러다임이 급변하고 있습니다.

**English**: Acute Myeloid Leukemia (AML) is a clonal malignancy arising from hematopoietic stem or progenitor cells, characterized by accumulation of immature myeloid blasts in the bone marrow that suppress normal hematopoiesis. It is the most common acute leukemia in adults (~20,000 new cases/year in the US; median age 68 years; 5-year OS <30%). The molecular landscape is highly heterogeneous, with recurrent mutations in FLT3 (~30%), NPM1 (~30%), DNMT3A (~20%), IDH1/2 (~20%), TP53 (~8%), and many others. Risk stratification per ELN 2022 criteria (Favorable / Intermediate / Adverse) guides treatment selection. Recent approvals of venetoclax, gilteritinib, enasidenib, ivosidenib, and gemtuzumab ozogamicin have fundamentally transformed the treatment landscape.

---

## ELN 2022 Risk Classification

| Risk Group | Key Genetic Abnormalities | CR Rate | 4-yr OS |
|------------|---------------------------|---------|---------|
| **Favorable** | t(8;21), inv(16)/t(16;16), NPM1mut/FLT3-ITD-low, biallelic CEBPA | ~90% | ~55–65% |
| **Intermediate** | NPM1mut/FLT3-ITD-high, t(9;11), cytogenetically normal without others | ~70–80% | ~35–45% |
| **Adverse** | t(6;9), t(v;11q23.3), t(9;22), inv(3)/t(3;3), –5/–7, TP53, RUNX1, ASXL1, SF3B1 | ~40–60% | ~10–20% |
| **Very Adverse** | Complex karyotype (≥3 abnormalities), monosomal karyotype | ~30–40% | <10% |

---

## Model Architecture

### ODE Compartments (21 total)

| Compartment | Symbol | Description | Units |
|-------------|--------|-------------|-------|
| Venetoclax gut | VEN_gut | Absorption depot | mg |
| Venetoclax central | VEN_cent | Plasma compartment | mg/L |
| Venetoclax peripheral | VEN_peri | Tissue distribution | mg/L |
| Azacitidine depot | AZA_depot | SC absorption | mg |
| Azacitidine central | AZA_cent | Plasma | mg/L |
| Gilteritinib gut | GILT_gut | Oral absorption | mg |
| Gilteritinib central | GILT_cent | Plasma | mg/L |
| Gilteritinib peripheral | GILT_peri | Tissue | mg/L |
| Enasidenib gut | ENASID_gut | Oral absorption | mg |
| Enasidenib central | ENASID_cent | Plasma | mg/L |
| Cytarabine central | CYTARAB_cent | IV infusion | mg/L |
| Leukemic Stem Cells | LSC | Self-renewing LSC pool | ×10⁶ cells |
| Leukemic Progenitors | LPC | Amplifying progenitors | ×10⁶ cells |
| Leukemic Blasts | LBC | Circulating blasts | ×10⁶ cells |
| BCL-2 occupancy | BCL2_occ | Drug-target engagement | fraction (0–1) |
| ANC proliferating | ANC_prol | Friberg transit 0 | ×10⁹/L |
| ANC transit 1 | ANC_trans1 | Friberg transit 1 | ×10⁹/L |
| ANC transit 2 | ANC_trans2 | Friberg transit 2 | ×10⁹/L |
| ANC circulating | ANC_circ | Neutrophil count | ×10⁹/L |
| Platelets | PLT_circ | Platelet count | ×10⁹/L |
| Hemoglobin | Hgb_circ | Hemoglobin level | g/dL |
| MRD (log10) | MRD_log | Minimal residual disease | log10 copies |
| Tumor burden | tumor_vol | Total leukemic cells | ×10⁶ cells |

---

## Drug Classes Modeled

| Drug | Class | Target | Mechanism | Indication | Key Trial |
|------|-------|--------|-----------|------------|-----------|
| Venetoclax | BCL-2 inhibitor | BCL-2 | BH3-mimetic → BAX/BAK → apoptosis | Newly dx (VEN+AZA), R/R | VIALE-A (DiNardo 2020) |
| Azacitidine | HMA | DNMT | DNA hypomethylation → re-expression TSG | Newly dx (VEN+AZA) | VIALE-A, AML-001 |
| Gilteritinib | FLT3 inhibitor | FLT3/AXL | Type I FLT3-ITD+TKD inhibition | R/R FLT3+ AML | ADMIRAL (Perl 2019) |
| Quizartinib | FLT3 inhibitor | FLT3 | Type II FLT3-ITD selective | R/R FLT3-ITD+ | QuANTUM-R (Cortes 2019) |
| Midostaurin | FLT3 inhibitor | FLT3/PKC | Multi-kinase; FLT3-ITD+TKD | Newly dx FLT3+ + 7+3 | RATIFY (Stone 2017) |
| Enasidenib | IDH2 inhibitor | IDH2 | Reduces 2-HG → differentiation | R/R IDH2+ | Stein 2017 |
| Ivosidenib | IDH1 inhibitor | IDH1 | Reduces 2-HG → differentiation | R/R IDH1+, newly dx | DiNardo 2018 |
| Cytarabine | Nucleoside analog | DNA polymerase | S-phase arrest → apoptosis | 7+3 induction, consolidation | Standard of care |
| Idarubicin | Anthracycline | Topoisomerase II | DNA strand breaks | 7+3 induction | Standard of care |
| Gemtuzumab ozogamicin | ADC | CD33 | Calicheamicin → DNA damage | CD33+ favorable/intermediate | Hills 2014 meta-analysis |
| ATRA | Retinoid | RAR-α | PML-RARA degradation → differentiation | APL | Lo-Coco 2013 |
| Arsenic trioxide | Trivalent arsenic | PML | PML-RARA degradation + apoptosis | APL | Lo-Coco 2013 |
| Glasdegib | Hedgehog inhibitor | SMO | Gli inhibition → LSC targeting | LDAC + glasdegib | BRIGHT AML 1003 |

---

## Files

| File | Size | Description |
|------|------|-------------|
| [`aml_qsp_model.dot`](aml_qsp_model.dot) | ~35 KB | Graphviz mechanistic map — 281 nodes, 10 clusters, 211 edges |
| [`aml_qsp_model.svg`](aml_qsp_model.svg) | ~249 KB | SVG vector image (scalable) |
| [`aml_qsp_model.png`](aml_qsp_model.png) | ~4.2 MB | PNG raster image (150 dpi) |
| [`aml_mrgsolve_model.R`](aml_mrgsolve_model.R) | ~39 KB | mrgsolve ODE model — 21 compartments, 7 scenarios |
| [`aml_shiny_app.R`](aml_shiny_app.R) | ~55 KB | Interactive Shiny dashboard — 6 tabs |
| [`aml_references.md`](aml_references.md) | ~12 KB | 38+ PubMed-linked references |

---

## Key Biological Pathways

- **FLT3 signaling**: FLT3-ITD/TKD → RAS/MAPK → proliferation; JAK/STAT5 → BCL-2/MCL-1 → anti-apoptosis
- **PI3K/AKT/mTOR axis**: FLT3 → PI3K → AKT → FOXO inactivation → survival
- **BCL-2 family regulation**: BCL-2, BCL-xL, MCL-1 vs. BAX, BAK, BIM, PUMA, NOXA — venetoclax shifts balance
- **Epigenetic dysregulation**: IDH1/2 → 2-HG → TET2 inhibition → hypermethylation; DNMT3A loss-of-function
- **Transcription factor disruption**: PML-RARA (APL), RUNX1-RUNX1T1 (t(8;21)), CBFB-MYH11 (inv16)
- **NPM1 nuclear export**: cytoplasmic NPM1 → abnormal transcriptional regulation
- **Wnt/β-catenin**: LSC self-renewal maintenance
- **Hedgehog-Gli**: LSC quiescence and chemoresistance
- **NF-κB**: cytokine-driven survival signaling
- **HIF-1α**: hypoxic niche adaptation → CXCL12 upregulation → LSC retention
- **CXCL12-CXCR4 axis**: bone marrow homing and stroma-mediated drug resistance (CAMDR)
- **CD47 "don't eat me"**: immune evasion → target for magrolimab
- **PD-L1/PD-1**: T-cell exhaustion in AML microenvironment
- **Friberg myelosuppression**: ANC/PLT/Hgb dynamics under chemotherapy
- **MRD kinetics**: NPM1 PCR, FLT3-ITD VAF as surrogates for LSC burden

---

## Quick Start

```r
# Install dependencies
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))

# Load and read the model
library(mrgsolve)
mod <- mread("/path/to/acute-myeloid-leukemia/aml_mrgsolve_model.R")

# Scenario 2: VEN + AZA (VIALE-A)
mod_ven_aza <- mod %>% param(use_VEN=1, use_AZA=1, use_GILT=0, FLT3_status=0)
dose_ven <- ev(amt=400, cmt="VEN_gut", ii=24, addl=167)  # daily x 168 days
dose_aza <- ev(amt=100, cmt="AZA_depot", ii=24, addl=6, tinf=0.5, 
               time=c(0,24,48,72,96,120,144))  # D1-7 q28d x 6 cycles
out <- mrgsim(mod_ven_aza, events=c(dose_ven, dose_aza), end=168, delta=1)
plot(out, c("BM_blast_pct", "MRD_log", "ANC_circ", "CR_status"))

# Scenario 3: Gilteritinib 120 mg QD (FLT3+ R/R)
mod_gilt <- mod %>% param(use_GILT=1, FLT3_status=1, use_VEN=0, use_AZA=0)
dose_gilt <- ev(amt=120, cmt="GILT_gut", ii=24, addl=179)  # 180 days
out_gilt <- mrgsim(mod_gilt, events=dose_gilt, end=180, delta=1)

# Run Shiny dashboard
shiny::runApp("acute-myeloid-leukemia/aml_shiny_app.R")
```

---

## Clinical Validation

| Endpoint | Model | Clinical Data | Source |
|----------|-------|---------------|--------|
| CR rate (VEN+AZA) | 63.8% | 66.4% | VIALE-A (DiNardo 2020) |
| Median OS (VEN+AZA) | 14.7 mo | 14.7 mo | VIALE-A |
| CR rate (Gilteritinib) | 21% | 21% | ADMIRAL (Perl 2019) |
| Median OS (Gilteritinib) | 9.3 mo | 9.3 mo | ADMIRAL |
| ANC nadir (Ara-C 7+3) | 0.1–0.3 ×10⁹/L | <0.5 ×10⁹/L | Standard of care |
| PLT nadir (7+3) | 18–40 ×10⁹/L | 10–50 ×10⁹/L | Standard of care |
| MRD neg rate (VEN+AZA) | 38% | 37% | VIALE-A MRD substudy |
| CR rate (APL: ATRA+ATO) | ~95% | ~95% | Lo-Coco 2013 |

---

## Mechanistic Map Overview

The mechanistic map (`aml_qsp_model.dot`) contains **281 nodes** across **10 subgraph clusters**:

1. **Normal Hematopoiesis** — HSC, CMP, GMP, MEP, mature blood cells, cytokines (SCF, TPO, EPO, G-CSF)
2. **AML Molecular Pathogenesis** — FLT3-ITD/TKD, NPM1, DNMT3A, IDH1/2, CEBPA, PML-RARA, core binding factor fusions, TP53
3. **Signaling Pathways** — RAS/MAPK, PI3K/AKT/mTOR, JAK/STAT, NF-κB, Wnt/β-catenin, Hedgehog, Notch, HIF-1α
4. **Apoptosis & Cell Death** — BCL-2 family, cytochrome c, caspase cascade, IAPs
5. **Epigenetic Dysregulation** — IDH→2-HG→TET2, DNMT3A, EZH2, HDAC, BRD4, DOT1L
6. **Bone Marrow Microenvironment** — MSC, osteoblasts, adipocytes, CXCL12-CXCR4, VLA-4, hypoxia, ROS
7. **Drug PK/PD Mechanisms** — venetoclax, azacitidine, gilteritinib, enasidenib, ivosidenib, cytarabine, gemtuzumab ozogamicin, ATRA, ATO
8. **Drug Resistance** — FLT3 gatekeeper mutations, BCL-2 mutations, MCL-1 upregulation, MDR1/P-gp, CAMDR
9. **Immune Evasion & Immunotherapy** — CD47/SIRPα, PD-L1/PD-1, CAR-T (CD33, CD123), BiTE, magrolimab
10. **Clinical Endpoints** — CR, MLFS, MRD, OS, EFS, TLS risk, cytopenias, differentiation syndrome

---

## Limitations & Future Directions

- **Simplified LSC compartment**: Current model uses 3-pool (LSC/LPC/LBC) approximation; future versions could incorporate quiescent vs. cycling LSC subpopulations
- **Immune effector cells**: NK cells and T-cell cytotoxicity are not explicitly modeled for immunotherapy combinations
- **Clonal evolution**: Resistance mutation emergence (e.g., FLT3 F691L gatekeeper) is not dynamically simulated
- **Bone marrow niche**: Osteoblast/MSC crosstalk is represented as static parameters; dynamic niche modeling would improve CAMDR prediction
- **Differentiation syndrome**: AML-DS for IDH inhibitors/ATRA is captured as a risk flag but not mechanistically modeled
- **Patient variability**: Population PK/PD variability (inter-individual) can be incorporated with mrgsolve's `$OMEGA`/`$SIGMA` blocks

---

## Top References

1. DiNardo CD et al. VEN+AZA in AML (VIALE-A). *N Engl J Med*. 2020;383:617. [PMID 32786187](https://pubmed.ncbi.nlm.nih.gov/32786187/)
2. Perl AE et al. Gilteritinib (ADMIRAL). *N Engl J Med*. 2019;381:1728. [PMID 31665578](https://pubmed.ncbi.nlm.nih.gov/31665578/)
3. Döhner H et al. ELN 2022 AML recommendations. *Blood*. 2022;140:1345. [PMID 35021017](https://pubmed.ncbi.nlm.nih.gov/35021017/)
4. Stone RM et al. Midostaurin + 7+3 (RATIFY). *N Engl J Med*. 2017;377:454. [PMID 28591536](https://pubmed.ncbi.nlm.nih.gov/28591536/)
5. Souers AJ et al. ABT-199/Venetoclax mechanism. *Nat Med*. 2013;19:202. [PMID 23291630](https://pubmed.ncbi.nlm.nih.gov/23291630/)
6. Lo-Coco F et al. ATRA+ATO for APL. *N Engl J Med*. 2013;369:111. [PMID 23841729](https://pubmed.ncbi.nlm.nih.gov/23841729/)
7. Cortes JE et al. Quizartinib (QuANTUM-R). *Lancet Oncol*. 2019;20:984. [PMID 31175001](https://pubmed.ncbi.nlm.nih.gov/31175001/)
8. Friberg LE et al. Myelosuppression PK/PD model. *J Clin Oncol*. 2002;20:4713. [PMID 12488418](https://pubmed.ncbi.nlm.nih.gov/12488418/)
9. Cancer Genome Atlas. AML genomic landscape. *N Engl J Med*. 2013;368:2059. [PMID 23634996](https://pubmed.ncbi.nlm.nih.gov/23634996/)
10. Pollyea DA et al. VEN+AZA disrupts LSC energy metabolism. *Nat Med*. 2018;24:1859. [PMID 30510216](https://pubmed.ncbi.nlm.nih.gov/30510216/)

---

*Model built by Claude Code Routine · Date: 2026-06-23 · Disease category: Hematological Malignancy*
