# Diffuse Large B-Cell Lymphoma (DLBCL) — QSP Model

[![Mechanistic Map](dlbcl_qsp_model.png)](dlbcl_qsp_model.svg)

---

## Disease Overview

### Epidemiology
Diffuse Large B-Cell Lymphoma (DLBCL) is the most common aggressive non-Hodgkin lymphoma, accounting for approximately 25–30% of all NHL cases. In the United States, roughly 27,000 new cases are diagnosed annually. Median age at diagnosis is 64–70 years, with a slight male predominance (M:F ≈ 1.5:1). Globally, the annual incidence is estimated at 150,000–200,000 cases. With modern R-CHOP immunochemotherapy, approximately 60–65% of patients are cured; the remaining 30–40% relapse or are refractory (R/R DLBCL), with historically poor outcomes (median OS ~6 months in the pre-CAR-T era).

### Molecular Subtypes (Cell-of-Origin)
| Subtype | Frequency | Key Mutations | 5-yr OS | Pathway Dependence |
|---------|-----------|---------------|---------|-------------------|
| **GCB** (Germinal Center B-cell) | ~55% | EZH2, BCL2 translocation, KMT2D | ~65–75% | PI3K-AKT, BCL2 |
| **ABC** (Activated B-cell) | ~35% | MYD88 L265P, CD79B, NOTCH1 | ~45–55% | BCR-NF-κB, JAK-STAT |
| **Unclassifiable** | ~10% | Mixed | ~50–60% | Variable |

### Genetic Subtypes (Schmitz/Chapuy 2018 Classification)
| Genetic Subtype | Key Features | COO Enrichment |
|-----------------|--------------|----------------|
| **MCD** (MYD88/CD79B) | MYD88 L265P + CD79B mutations | ABC |
| **BN2** (BCL6/NOTCH2) | BCL6 fusions + NOTCH2 mutations | GCB/ABC border |
| **N1** (NOTCH1) | NOTCH1 mutations | ABC |
| **EZB** (EZH2/BCL2) | EZH2 mutations + BCL2 translocation | GCB |
| **ST2** (SGK1/TET2) | SGK1 + TET2 mutations | GCB |

### Double-Hit / Double-Expressor Lymphoma
- **Double-Hit Lymphoma (DHL)**: MYC + BCL2 (or BCL6) translocation; ~5–10% of DLBCL; poor prognosis with R-CHOP (2-yr PFS ~25–35%); treated with DA-EPOCH-R
- **Double-Expressor Lymphoma (DEL)**: MYC protein >40% + BCL2 protein >50% by IHC; ~30% of DLBCL; intermediate-poor prognosis

### IPI Risk Stratification
| IPI Score | Risk Group | 5-yr OS (R-CHOP era) |
|-----------|------------|---------------------|
| 0–1 | Low | ~80% |
| 2 | Low-Intermediate | ~70% |
| 3 | High-Intermediate | ~55% |
| 4–5 | High | ~35–45% |

IPI factors: Age >60, LDH > ULN, ECOG PS ≥2, Stage III/IV, Extranodal sites ≥2

---

## Key Pathway Summary

| Pathway | Role in DLBCL | Targetable Nodes | Drugs |
|---------|--------------|-----------------|-------|
| **BCR signaling → NF-κB** | Survival in ABC subtype; constitutive NF-κB | BTK, SYK, PLCγ2, PKCβ, IKK | Ibrutinib, zanubrutinib (BTK); R406 (SYK) |
| **PI3K-AKT-mTOR** | Proliferation, survival; activated by BCR & CD40 | PI3Kδ, AKT, mTORC1, FOXO1 | Idelalisib, copanlisib (PI3K); everolimus (mTOR) |
| **BCL2/apoptosis** | Anti-apoptotic barrier; BCL2 overexpression in GCB | BCL2, MCL1, BCL-XL, BAX/BAK | Venetoclax (BCL2); navitoclax (BCL2/BCL-XL) |
| **CD20 (surface antigen)** | B-cell surface marker; target for ADCC/CDC | CD20 (MS4A1) | Rituximab, obinutuzumab, ofatumumab |
| **GCB vs ABC determination** | BCL6 represses IRF4/BLIMP1 in GCB; IRF4 promotes plasma-cell in ABC | BCL6, IRF4, MYC, BLIMP1 | EZH2 inhibitors (tazemetostat) |
| **Epigenetics** | EZH2 gain-of-function (GCB); CREBBP/KMT2D loss-of-function | EZH2, CREBBP, EP300 | Tazemetostat (EZH2 inhibitor) |
| **DNA damage / cell cycle** | TP53 mutation in ~20%; CDK4/6 in proliferation | TP53, ATM, CDK4/6, RB | CDK4/6i (palbociclib); MDM2 inhibitors |
| **Tumor microenvironment** | PD-L1 upregulation, immune evasion, TAM M2 polarization | PD-1/PD-L1, LAG-3, TIM-3 | Pembrolizumab, nivolumab (anti-PD-1) |
| **ADC targets** | CD79b expressed in >90% DLBCL | CD79b | Polatuzumab vedotin (pola-vedotin) |

---

## Model Files

| File | Description | Size |
|------|-------------|------|
| `dlbcl_qsp_model.dot` | Graphviz mechanistic map (183 nodes, 10 clusters) | 39 KB |
| `dlbcl_qsp_model.svg` | Rendered SVG (scalable, interactive) | 220 KB |
| `dlbcl_qsp_model.png` | Rendered PNG (150 dpi) | 3.3 MB |
| `dlbcl_mrgsolve_model.R` | mrgsolve ODE PK/PD model + 6 scenarios | 32 KB |
| `dlbcl_shiny_app.R` | Interactive Shiny dashboard (6 tabs) | 52 KB |
| `dlbcl_references.md` | 43 curated PubMed references | 16 KB |

---

## Model Specifications

### ODE Compartments (18 total)
| # | Compartment | Description |
|---|-------------|-------------|
| 1 | `CRIT` | Rituximab central (mg/L) |
| 2 | `CPER` | Rituximab peripheral |
| 3 | `CPOLA` | Polatuzumab vedotin ADC |
| 4 | `CPOLA_PAY` | MMAE payload (released) |
| 5 | `CVEN_C` | Venetoclax central (mg/L) |
| 6 | `CVEN_P` | Venetoclax peripheral |
| 7 | `CIBRUT` | Ibrutinib plasma |
| 8 | `CRCHOP` | R-CHOP cytotoxic combined |
| 9 | `TUMOR_GCB` | GCB tumor burden |
| 10 | `TUMOR_ABC` | ABC tumor burden |
| 11 | `BCR_ACT` | BCR signaling activity (0–1) |
| 12 | `NFKB` | NF-κB transcription activity |
| 13 | `BCL2_EXP` | BCL2 expression level |
| 14 | `APOP_SIG` | Apoptotic signal accumulation |
| 15 | `NK_CELLS` | NK cell immune effector level |
| 16 | `CTL_CELLS` | CD8+ CTL level (includes CAR-T) |
| 17 | `PDL1_EXP` | PD-L1 expression |
| 18 | `LDH` | LDH tumor burden surrogate |

### Treatment Scenarios (6)
| # | Scenario | Regimen | Clinical Basis |
|---|----------|---------|----------------|
| 1 | No Treatment | Tumor progression only | Natural history control |
| 2 | R-CHOP | Rituximab 375 mg/m² d1 + CHOP d1–5, Q21d × 6 | GOYA, Coiffier 2002 |
| 3 | Pola-R-CHP | Pola 1.8 mg/kg + R-CHP Q21d × 6 | POLARIX (Tilly 2022) |
| 4 | R-CHOP + Venetoclax | R-CHOP + venetoclax 800 mg/day continuous | Morschhauser 2021 |
| 5 | Ibrutinib + R-CHOP | Ibrutinib 560 mg/day + R-CHOP (ABC subtype) | PHOENIX trial |
| 6 | CAR-T | Simplified effector:target expansion model | ZUMA-1, TRANSCEND |

### Shiny App Tabs (6)
| Tab | Contents |
|-----|---------|
| Patient Profile | Age, gender, IPI, COO subtype, MYC/BCL2 status, ECOG, treatment selection |
| Drug PK | Concentration-time profiles, Cmax/AUC, target occupancy for each drug |
| PD Key Indicators | Tumor burden, BCR activity, BCL2 occupancy, NK/CTL, PD-L1, LDH |
| Clinical Endpoints | Waterfall plot, spider plot, PFS/OS Kaplan-Meier curves |
| Scenario Comparison | Side-by-side tumor burden, response rates, 2yr PFS/OS bar charts |
| Biomarker Explorer | BCL2 vs venetoclax, BTK occupancy vs ibrutinib dose, PD-L1 scatter |

---

## How to Run the Model

### Prerequisites
```r
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "patchwork",
                   "shiny", "bslib", "plotly", "DT"))
```

### Run the mrgsolve ODE model
```r
# Navigate to the model directory
setwd("/path/to/diffuse-large-b-cell-lymphoma/")

# Source and run
source("dlbcl_mrgsolve_model.R")

# This will:
# 1. Define the 18-compartment ODE model
# 2. Simulate all 6 treatment scenarios
# 3. Generate PK, tumor burden, BCR, BCL2, immune, LDH, and comparison plots
# 4. Print clinical summary statistics
```

### Launch the Shiny dashboard
```r
setwd("/path/to/diffuse-large-b-cell-lymphoma/")
shiny::runApp("dlbcl_shiny_app.R")
```

### Render the mechanistic map
```bash
# SVG output
dot -Tsvg dlbcl_qsp_model.dot -o dlbcl_qsp_model.svg

# PNG at 150 dpi
dot -Tpng -Gdpi=150 dlbcl_qsp_model.dot -o dlbcl_qsp_model.png

# Alternative layout engine (for large graphs)
neato -Tsvg dlbcl_qsp_model.dot -o dlbcl_qsp_model_neato.svg
```

---

## Key Clinical Trial Results Referenced

| Trial | Treatment | N | Setting | PFS/EFS | OS | Key Finding |
|-------|-----------|---|---------|---------|-----|-------------|
| **Coiffier 2002** (GELA) | R-CHOP vs CHOP | 399 | 1L elderly | 2yr EFS 57% vs 38% | 2yr OS 70% vs 57% | Rituximab added to CHOP becomes SoC |
| **GOYA** (Vitolo 2017) | G-CHOP vs R-CHOP | 1418 | 1L | 3yr PFS 67% vs 67% | HR 0.92 | Obinutuzumab not superior to rituximab |
| **REMoDL-B** (Davies 2017) | R-CHOP±bortezomib | 806 | 1L DH/DE | 2yr PFS 68% vs 67% | Neutral | No benefit in double-expressor |
| **PHOENIX** (Younes 2019) | Ibr+R-CHOP vs R-CHOP | 838 | 1L | Neg overall; HR 0.75 in ≤60yr ABC | — | ABC younger pts may benefit |
| **POLARIX** (Tilly 2022) | Pola-R-CHP vs R-CHOP | 879 | 1L IPI≥2 | 2yr PFS 76.7% vs 70.2%, HR 0.73 | Neutral | Pola-R-CHP improves PFS in IPI≥2 |
| **ZUMA-1** (Neelapu 2017) | Axi-cel | 101 | R/R ≥2L | 18mo PFS 34% | 18mo OS 52% | First CAR-T approval; durable remissions |
| **JULIET** (Schuster 2019) | Tisa-cel | 93 | R/R ≥2L | 12mo PFS 34% | 12mo OS 49% | Confirmed CAR-T benefit |
| **TRANSCEND** (Abramson 2020) | Liso-cel | 256 | R/R ≥2L | 12mo PFS 44% | 12mo OS 58% | Lowest grade ≥3 neurotox; broad eligibility |
| **L-MIND** (Salles 2020) | Tafa+Len | 80 | R/R 2L | 12mo PFS 37% | 12mo OS 63% | Outpatient option for transplant-ineligible |

---

## Calibrated Model Parameters Summary

| Parameter | Value | Source |
|-----------|-------|--------|
| Rituximab CL | 0.23 L/day | Mould et al. 2007, CPT |
| Rituximab V1 | 3.1 L | Mould et al. 2007 |
| Venetoclax CL | 12 L/h | Freise et al. 2019, CPT:PSP |
| Ibrutinib CL/F | 60 L/h | De Zwart et al. 2016 |
| Pola ADC CL | 0.8 L/day | Deng et al. 2019 |
| GCB tumor kg | 0.035/day | Calibrated to GOYA SoC arm |
| ABC tumor kg | 0.048/day | Calibrated to PHOENIX control arm |
| BCL2 EC50 (venetoclax) | 0.5 µg/mL | Souers et al. 2013, Nature Med |
| BTK EC50 (ibrutinib) | 0.3 µg/mL | Herman et al. 2011, Blood |

---

## Mechanistic Map Clusters

The mechanistic map (`dlbcl_qsp_model.svg`) contains **183 nodes** across **10 subgraph clusters**:

1. **BCR Signaling & NF-κB** (rose, 23 nodes) — BCR → LYN → SYK → BTK → PLCγ2 → PKCβ → CARMA1/BCL10/MALT1 → IKK → NF-κB → IRF4/BLIMP1
2. **PI3K-AKT-mTOR** (light blue, 18 nodes) — PI3Kδ/PTEN/PIP2/PIP3 → AKT → TSC1/2 → RHEB → mTORC1/2 → S6K/4E-BP1/FOXO1
3. **Cell Survival & Apoptosis** (light purple, 19 nodes) — BCL2/BCL-XL/MCL1 ↔ BAX/BAK → cytochrome-c → APAF1 → caspase-9/3/7
4. **Germinal Center Biology** (light green, 17 nodes) — BCL6/IRF4 toggle, MYC, AID/AICDA, SHM/CSR, GCB/LZB/DZB zones
5. **Tumor Microenvironment** (light amber, 21 nodes) — TAM M1/M2, Treg, CD8+ CTL, NK, PD-1/PD-L1/LAG-3/TIM-3/CTLA-4, cytokines
6. **Drug PK Compartments** (cyan boxes, 12 nodes) — rituximab, pola-vedotin, venetoclax, ibrutinib, R-CHOP, obinutuzumab, lenalidomide
7. **DNA Damage & Cell Cycle** (light pink, 20 nodes) — TP53/ATM/ATR/CHK1-2 → p21 → CDK4/6 → RB/E2F1
8. **Clinical Endpoints** (lime diamonds, 15 nodes) — IPI, LDH, Ki-67, PET-CT, CR/PR/PD, OS/PFS, ECOG
9. **Cytokine Network** (light deep-purple, 18 nodes) — IL-6/IL-10/IL-21/IFN-γ/TNF-α/BAFF/APRIL/CXCL12-13
10. **Epigenetics & Mutations** (light deep-orange, 18 nodes) — EZH2/CREBBP/EP300/KMT2D, H3K marks, key driver mutations

---

## Citation

If using this model for research or education, please cite:

> QSP Disease Model Library, DLBCL QSP Model v1.0 (2026). Claude Code Routine (CCR) auto-generated.
> Key references: Alizadeh et al. 2000 (PMID: 10676951), Tilly et al. 2022 (PMID: 34904796), Abramson et al. 2020 (PMID: 33068767).

---

*For full reference list, see [dlbcl_references.md](dlbcl_references.md)*
