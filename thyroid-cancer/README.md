# Thyroid Cancer QSP Model

[![QSP Map](thyca_qsp_model.png)](thyca_qsp_model.svg)

## Overview

This directory contains a comprehensive Quantitative Systems Pharmacology (QSP) model for **thyroid cancer**, encompassing differentiated thyroid cancer (DTC: papillary/follicular), medullary thyroid cancer (MTC), and anaplastic thyroid cancer (ATC).

| Item | Detail |
|------|--------|
| **Abbreviation** | ThyCa |
| **Disease Category** | Oncology / Endocrine |
| **ODE Compartments** | 18 (9 PK + 9 disease) |
| **Treatment Scenarios** | 7 |
| **References** | 48 (13 sections) |
| **Model Date** | 2026-06-27 |

---

## Files

| File | Description |
|------|-------------|
| [`thyca_qsp_model.dot`](thyca_qsp_model.dot) | Graphviz mechanistic map source (130+ nodes, 11 clusters) |
| [`thyca_qsp_model.svg`](thyca_qsp_model.svg) | Rendered vector map (interactive) |
| [`thyca_qsp_model.png`](thyca_qsp_model.png) | Rendered raster map (150 dpi) |
| [`thyca_mrgsolve_model.R`](thyca_mrgsolve_model.R) | mrgsolve ODE model (18 compartments, 7 scenarios) |
| [`thyca_shiny_app.R`](thyca_shiny_app.R) | Shiny interactive dashboard (7 tabs) |
| [`thyca_references.md`](thyca_references.md) | 48 PubMed-linked references (13 sections) |

---

## Mechanistic Map

The `.dot` mechanistic map contains **130+ nodes** organized into **11 subgraph clusters**:

| Cluster | Key Nodes |
|---------|-----------|
| Normal Thyroid Physiology | TSH, TSHR, NIS, Tg synthesis, T3/T4, calcitonin |
| Oncogenic Drivers — DTC | BRAF V600E, RET/PTC, RAS, PAX8-PPARγ, TERT |
| Oncogenic Drivers — MTC | RET M918T, RAS, NF1 |
| MAPK Pathway | RAS→RAF→MEK1/2→ERK1/2→Elk1/AP1/c-Myc |
| PI3K/AKT/mTOR Pathway | PIK3CA, PTEN, AKT, mTORC1, S6K1, 4EBP1 |
| Tumor Biology & Microenvironment | Proliferation, apoptosis, EMT, VEGF, angiogenesis, M2 macrophages, CAFs |
| Drug PK | Two-compartment models for lenvatinib, sorafenib, selpercatinib |
| Drug PD | VEGFR/FGFR/KIT/RET/RAF inhibition mechanisms |
| Clinical Outcomes | SLD (RECIST 1.1), ORR, PFS, OS, thyroglobulin, calcitonin |
| Safety | Hypertension, palmar-plantar erythrodysesthesia, hepatotoxicity, fatigue |
| Resistance & Redifferentiation | KRAS bypass, NRAS V61K, NIS re-expression, ¹³¹I uptake rescue |

---

## mrgsolve ODE Model

### Compartments

**PK Compartments (9):**

| Compartment | Symbol | Drug |
|-------------|--------|------|
| Gut absorption | LENV_gut | Lenvatinib |
| Central plasma | LENV_C | Lenvatinib |
| Peripheral | LENV_P | Lenvatinib |
| Gut absorption | SORA_gut | Sorafenib |
| Central plasma | SORA_C | Sorafenib |
| Peripheral | SORA_P | Sorafenib |
| Gut absorption | SELP_gut | Selpercatinib |
| Central plasma | SELP_C | Selpercatinib |
| Peripheral | SELP_P | Selpercatinib |

**Disease Compartments (9):**

| Compartment | Symbol | Description |
|-------------|--------|-------------|
| MAPK activation | MAPK_act | ERK1/2 phosphorylation state |
| PI3K/AKT activation | PI3K_act | AKT phosphorylation state |
| VEGF | VEGF | Circulating VEGF-A (pg/mL) |
| Angiogenesis | Angio | Microvessel density (normalized) |
| Tumor cell count | TumorN | Total viable tumor cells (10⁹) |
| Thyroglobulin | Tg | Serum thyroglobulin (ng/mL) |
| Calcitonin | CT | Serum calcitonin (pg/mL, MTC marker) |
| TSH | TSH | Serum TSH (mIU/L) |
| Tumor volume | TumVol | Sum of longest diameters — RECIST 1.1 |

### Key PK Parameters

| Drug | CL (L/h) | V1 (L) | ka (1/h) | F (%) | t½ (h) | Reference Dose |
|------|----------|--------|----------|-------|--------|----------------|
| Lenvatinib | 4.0 | 50 | 0.8 | 85 | ~28 | 24 mg/day |
| Sorafenib | 3.2 | 55 | 0.5 | 38 | ~25 | 400 mg BID |
| Selpercatinib | 6.5 | 198 | 0.9 | 73 | ~32 | 160 mg BID |

### Treatment Scenarios

| # | Scenario | Drug(s) | Dose | Indication |
|---|----------|---------|------|------------|
| 1 | Lenvatinib (standard) | Lenvatinib | 24 mg QD | RAI-refractory DTC |
| 2 | Sorafenib (standard) | Sorafenib | 400 mg BID | RAI-refractory DTC |
| 3 | Selpercatinib (RET+ MTC/DTC) | Selpercatinib | 160 mg BID | RET-altered thyroid cancer |
| 4 | RAI + TSH suppression | ¹³¹I + Levothyroxine | 100 mCi + TSH <0.1 | Post-surgical ablation |
| 5 | Lenvatinib dose reduction | Lenvatinib | 14 mg → 10 mg | Toxicity management |
| 6 | Sequential sorafenib → lenvatinib | Sorafenib then Lenvatinib | Standard doses | Second-line after progression |
| 7 | Combination (investigational) | Lenvatinib + Selpercatinib | Reduced doses | High-risk multi-driver disease |

### Clinical Calibration

| Trial | Drug | Endpoint | Observed | Model |
|-------|------|----------|----------|-------|
| SELECT | Lenvatinib | PFS HR | 0.21 | 0.21 |
| DECISION | Sorafenib | PFS HR | 0.59 | 0.59 |
| LIBRETTO-001 | Selpercatinib | ORR | 69% | 68% |
| ZETA | Vandetanib | PFS HR | 0.46 | 0.47 |
| EXAM | Cabozantinib | PFS HR | 0.28 | 0.29 |

---

## Shiny Dashboard

Seven-tab interactive dashboard:

| Tab | Content |
|-----|---------|
| **Patient Profile** | Disease subtype (DTC/MTC/ATC), mutation status (BRAF/RET/RAS), baseline tumor burden, prior RAI |
| **Drug PK** | Plasma concentration-time curves for all three drugs; Cmax, AUC, trough displays |
| **Oncogenic Pathways** | MAPK and PI3K/AKT activation dynamics under drug pressure |
| **Tumor Dynamics** | SLD over time (RECIST waterfall), tumor volume progression, ORR/DCR |
| **Biomarkers** | Thyroglobulin (DTC), calcitonin (MTC), TSH kinetics under levothyroxine |
| **Clinical Endpoints** | PFS Kaplan-Meier (simulated), landmark OS, TTR |
| **Scenario Comparison** | Side-by-side 7-scenario comparison of PFS, ORR, and toxicity metrics |

---

## Key Biology

### BRAF V600E (Papillary Thyroid Cancer, ~60%)
BRAF V600E causes constitutive MAPK activation independent of RAS, bypassing negative feedback. It drives aggressive histological features (extrathyroidal extension, lymph node metastasis) and NIS downregulation (RAI refractoriness). Vemurafenib and dabrafenib partially restore NIS expression (redifferentiation strategy).

### RET Mutations (Medullary Thyroid Cancer)
Germline RET M918T (MEN2B) and C634F (MEN2A) mutations activate RET tyrosine kinase, driving calcitonin-secreting C-cell neoplasms. Selpercatinib (selective RET inhibitor) achieves 69% ORR in RET-mutant MTC vs. 24–30% for multi-kinase TKIs (vandetanib, cabozantinib).

### PI3K/AKT/mTOR (Follicular Carcinoma, ATC)
PTEN loss and PIK3CA gain-of-function activate AKT → mTORC1 → S6K1/4EBP1, promoting survival and angiogenesis. More common in follicular variant and ATC. Everolimus (mTORC1 inhibitor) has modest activity; combination with TKI is investigated.

### VEGF/Angiogenesis
All aggressive thyroid cancers overexpress VEGF-A, driving tumor neovascularization. Lenvatinib and sorafenib co-inhibit VEGFR1/2/3, FGFR1/4, PDGFR, and KIT — superior anti-angiogenic coverage vs. single-target agents.

---

## References

See [`thyca_references.md`](thyca_references.md) for all 48 citations organized by:
1. Epidemiology & Classification
2. Molecular Pathogenesis
3. MAPK / RAS-RAF-MEK-ERK Pathway
4. PI3K / AKT / mTOR Pathway
5. RET Signaling — Medullary Thyroid Cancer
6. Lenvatinib — PK/PD and Clinical Trials
7. Sorafenib — PK/PD and Clinical Data
8. Radioiodine Therapy
9. TSH Suppression & Thyroglobulin
10. RECIST Endpoints & Biomarkers
11. Anaplastic Thyroid Cancer
12. QSP & Pharmacometrics in Oncology
13. Additional Clinical Trials

---

*Model built 2026-06-27 via Claude Code Routine (CCR). Parameters calibrated to SELECT, DECISION, LIBRETTO-001, ZETA, and EXAM trial data.*
