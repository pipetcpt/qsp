# Non-Hodgkin Lymphoma (DLBCL) — QSP Model

> **Disease:** Diffuse Large B-Cell Lymphoma (DLBCL) | **Category:** Hematologic Oncology  
> **Directory:** `non-hodgkin-lymphoma/` | **Abbreviation:** `nhl`  
> **Date Added:** 2026-06-25

[![DLBCL QSP Mechanistic Map](nhl_qsp_model.png)](nhl_qsp_model.svg)

---

## 1. Disease Overview

Diffuse Large B-Cell Lymphoma (DLBCL) is the **most common aggressive B-cell lymphoma**, accounting for 25–30% of all non-Hodgkin lymphoma cases worldwide (~150,000 new cases/year globally). It arises from germinal center (GC) B-cells or post-GC activated B-cells and can be cured with R-CHOP in ~60–65% of patients. The remaining ~35–40% who relapse/progress have poor outcomes (median OS ~6 months).

### Molecular Subtypes (Cell of Origin — COO)

| Subtype | Frequency | Key Drivers | R-CHOP Outcome |
|---------|-----------|-------------|----------------|
| **GCB-DLBCL** | ~50% | BCL-6+, t(14;18)/BCL-2, EZH2 Y641 | Better (~60% 5yr OS) |
| **ABC-DLBCL** | ~35% | MYD88 L265P, CD79A/B, NF-κB ↑ | Worse (~40% 5yr OS) |
| **Double-Hit** | ~8% | MYC + BCL-2/BCL-6 rearrangement | Very poor (~25% 5yr OS) |
| **PMBL** | ~3% | JAK-STAT, CD30+, 9p24.1 ampl | Intermediate |

---

## 2. Mechanistic Map

**File:** `nhl_qsp_model.dot` → `nhl_qsp_model.svg` / `nhl_qsp_model.png`

The mechanistic map contains **14 clusters** and **120+ nodes** covering:

| Cluster | Content |
|---------|---------|
| 1. B-Cell Development | Pro-B → Nadir B → GC B → Memory/Plasma; PAX5, BCL6, BLIMP1, IRF4, AID |
| 2. GC Reaction & DLBCL Origin | Dark/Light zone, FDC, Tfh, SHM, CSR; GCB/ABC/PMBL origin |
| 3. BCR Signaling (PI3K/AKT/mTOR) | Lyn→Syk→BTK→PLCγ2→Ca²⁺/NFAT; PI3K→AKT→mTORC1; RAS→MAPK |
| 4. NF-κB (Canonical/Non-canonical) | CBM complex→IKK→NF-κB p65/p50; NIK/BAFF-R non-canonical; MYD88 L265P |
| 5. MYC & Cell Cycle | c-MYC/t(8;14), Cyclin D/CDK4/6, RB1, E2F, p53/MDM2, DHL |
| 6. BCL-2 Family & Apoptosis | BCL-2/XL/MCL-1 vs BAX/BAK/BIM; Cytochrome C→Apaf1→Caspase cascade |
| 7. Epigenetics | EZH2 Y641, CREBBP/EP300, KMT2D, SWI/SNF, HDAC |
| 8. Tumor Microenvironment | TAM M1/M2, Tfh, NK, CD8+ CTL, Treg, CAF, VEGF, HIF-1α |
| 9. Immune Evasion | PD-1/PD-L1 (9p24.1), CD47/SIRPα, MHC-I loss (β2M), CD19/CD20 loss |
| 10. Rituximab PK/PD | 2-cmt + TMDD; CD20 binding; ADCC/CDC/direct apoptosis |
| 11. CHOP Components | Cyclophosphamide (4-OH-CPP), Doxorubicin, Vincristine, Prednisone |
| 12. Novel Agents | Polatuzumab vedotin (ADC), Ibrutinib (BTK), Venetoclax (BCL-2), Tafasitamab, Glofitamab |
| 13. CAR-T Therapy | Axi-cel/Tisa-cel/Liso-cel; CAR expansion; CRS; ICANS; CD19 loss resistance |
| 14. Clinical Endpoints | Lugano CR/PR/SD/PD; IPI/R-IPI; ctDNA MRD; PFS, OS, EFS24 |

---

## 3. mrgsolve ODE Model

**File:** `nhl_mrgsolve_model.R`

### Model Structure (22 Compartments)

**PK — 10 Compartments:**

| Compartment | Description |
|-------------|-------------|
| RuxCent | Rituximab central (mg) |
| RuxPeriph | Rituximab peripheral (mg) |
| CD20_free | Free CD20 antigen on tumor surface |
| CD20_RTX | CD20-Rituximab complex (TMDD) |
| CPP_active | 4-OH-Cyclophosphamide active metabolite (mg) |
| DoxCent | Doxorubicin central (mg) |
| DoxPeriph | Doxorubicin peripheral (mg) |
| VEN_gut | Venetoclax gut absorption (mg) |
| VEN_cent | Venetoclax central (mg) |
| IBR_cent | Ibrutinib central (mg) |

**PD — 12 Compartments:**

| Compartment | Description |
|-------------|-------------|
| Tumor | Tumor burden (logistic growth, T0=100 au) |
| BCR_signal | BCR/NF-κB signal activation (0-1) |
| BCL2_occ | BCL-2 occupancy by Venetoclax (0-1) |
| NK_cells | NK cell pool (ADCC effectors) |
| CD8_cells | CD8+ T-cell pool |
| ANC | Absolute Neutrophil Count (×10⁹/L) |
| Resistance | Drug resistance fraction (0-1) |
| Cum_RTX_dose | Cumulative rituximab exposure |
| Tumor_resp | Tumor response tracker |
| CRS_risk | CRS risk index |

### Key Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| Rituximab CL | 0.008 L/h | Berinstein 1998 |
| Rituximab Vc | 3.1 L | Maloney 1994 |
| RTX T½ (terminal) | ~21 days | Avivi 2014 |
| Tumor growth rate (kg) | 0.012/day (T2 ~58d) | Calibrated |
| R-CHOP CR rate (Day 126) | ~65% | Coiffier 2002 |
| kon(CD20-RTX) | 0.27/nM/h | Gibiansky 2014 |

### Treatment Scenarios (6)

| # | Scenario | Treatment | Clinical Reference |
|---|----------|-----------|-------------------|
| 1 | Untreated DLBCL | No treatment | Natural history |
| 2 | R-CHOP ×6 (Standard) | Rituximab + CHOP ×6 cycles | Coiffier 2002 NEJM |
| 3 | Pola-R-CHP ×6 (POLARIX) | Polatuzumab vedotin + R-CHP | Tilly 2022 NEJM |
| 4 | R-CHOP + Ibrutinib (ABC) | Ibrutinib + R-CHOP (ABC-DLBCL) | Younes 2019 NEJM |
| 5 | R-CHOP + Venetoclax | Venetoclax + R-CHOP (BCL-2 high) | Morschhauser 2021 JCO |
| 6 | R-CHOP (Double-Hit) | R-CHOP in DHL (high-risk subgroup) | Dunleavy 2013 JCO |

---

## 4. Shiny Dashboard

**File:** `nhl_shiny_app.R`

### 6-Tab Interactive Dashboard

| Tab | Content |
|-----|---------|
| **1. Patient Profile** | Disease subtype (GCB/ABC/DHL), IPI score, genomic drivers, BCL-2 status, prognosis by subtype |
| **2. Drug PK** | Rituximab 2-cmt + TMDD plot, CHOP concentration-time, CD20 receptor occupancy, add-on agent PK (venetoclax/ibrutinib), PK parameter table |
| **3. Tumor Dynamics** | SPD time course, response classification (CR/PR/SD/PD), NK/CD8+ immune effector dynamics, BCR/NF-κB signal suppression |
| **4. Clinical Endpoints** | Tumor reduction kinetics, endpoint summary table, waterfall plot (N=50), Lugano criteria reference, KM-like PFS curve |
| **5. Scenario Comparison** | 6-arm tumor burden comparison, scenario parameters, Day-126 ORR stacked bar, ANC myelosuppression by arm, resistance development |
| **6. Biomarkers & Toxicity** | ANC time profile, BCL-2 occupancy (venetoclax), CRS risk index, cumulative doxorubicin cardiotoxicity monitor, biomarker table |

### Key Interactive Controls

- **DLBCL Subtype**: GCB / ABC / Double-Hit (changes BCR signal, growth rate, prognosis)
- **BSA**: adjusts absolute doses for RTX, CPP, Dox
- **R-CHOP cycles**: 4–8 cycles
- **Dose sliders**: RTX, cyclophosphamide, doxorubicin per m²
- **Add-on agents**: toggle venetoclax (400–1200 mg) or ibrutinib (280–840 mg)
- **Simulation horizon**: 180–730 days

---

## 5. Model Calibration

| Clinical Trial | Treatment | Primary Endpoint | Observed | Model Target |
|----------------|-----------|-----------------|----------|--------------|
| Coiffier 2002 NEJM | R-CHOP vs CHOP | 3-yr EFS | 59% vs 38% | CR Day 126 ~65% |
| Tilly 2022 NEJM | Pola-R-CHP vs R-CHOP | 2-yr PFS | 76.7% vs 70.2% | ~6.5% PFS advantage |
| Younes 2019 NEJM (PHOENIX) | Ibr+R-CHOP | EFS (ABC) | HR 0.934 (ns) | Modest ABG benefit |
| Morschhauser 2021 JCO (CAVALLI) | VEN+R-CHOP | ORR (BCL-2+) | 88% vs 67% | ~21% ORR gain |
| Neelapu 2017 NEJM (ZUMA-1) | Axi-cel | ORR (refractory) | 82% | CAR-T module |

---

## 6. Running the Model

```r
# Requirements
install.packages(c("mrgsolve","dplyr","ggplot2","tidyr","shiny",
                   "shinydashboard","DT","scales"))

# Run simulation (deterministic)
source("nhl_mrgsolve_model.R")

# Launch Shiny dashboard
shiny::runApp("nhl_shiny_app.R")

# Render mechanistic map (Graphviz required)
# dot -Tsvg nhl_qsp_model.dot -o nhl_qsp_model.svg
# dot -Tpng -Gdpi=150 nhl_qsp_model.dot -o nhl_qsp_model.png
```

---

## 7. References Summary

50 PubMed-linked references organized by:
- Pivotal Clinical Trials (R-CHOP, POLARIX, PHOENIX, CAVALLI, ZUMA-1, JULIET)
- Disease Biology & Molecular Subtypes
- BCR/NF-κB Signaling
- BCL-2 & Apoptosis
- MYC & Double-Hit Lymphoma
- Novel Agents & CAR-T
- Rituximab/CHOP/Venetoclax/Ibrutinib PK-PD
- Response Assessment (Lugano criteria, ctDNA MRD)
- Epigenetics & Resistance

See [`nhl_references.md`](nhl_references.md) for full list.
