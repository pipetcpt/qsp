# Polycythemia Vera (진성 다혈증) — QSP Model

## Disease Overview

**Polycythemia Vera (PV)** is a BCR-ABL1-negative myeloproliferative neoplasm (MPN) characterized by clonal expansion of a hematopoietic stem cell harboring the **JAK2 V617F gain-of-function mutation** (>95% of patients). This leads to constitutive JAK-STAT5 signaling, resulting in:

- Uncontrolled erythropoiesis (elevated red cell mass)
- Thrombocytosis and leukocytosis
- Splenomegaly (extramedullary hematopoiesis)
- High thrombotic risk (DVT, stroke, MI, portal thrombosis)
- Risk of transformation to myelofibrosis (~10–20% at 20y) or AML (~2–7%)

**Prevalence:** ~2–3 per 100,000; median age at diagnosis ~60 years; M:F ratio ~1.2:1  
**Diagnosis (WHO 2016/2022):** Hct >49% (male) or >48% (female) + JAK2 mutation + BM hypercellularity

---

## Model Files

| File | Description |
|------|-------------|
| [`pv_qsp_model.dot`](pv_qsp_model.dot) | Graphviz mechanistic map (100+ nodes, 10 subgraph clusters) |
| [`pv_qsp_model.svg`](pv_qsp_model.svg) | Mechanistic map (SVG, scalable) |
| [`pv_qsp_model.png`](pv_qsp_model.png) | Mechanistic map (PNG, 150 dpi) |
| [`pv_mrgsolve_model.R`](pv_mrgsolve_model.R) | mrgsolve ODE model (16 compartments, 6 scenarios) |
| [`pv_shiny_app.R`](pv_shiny_app.R) | Interactive Shiny dashboard (7 tabs) |
| [`pv_references.md`](pv_references.md) | 58 PubMed references, categorized |

---

## Mechanistic Map

[![PV QSP Mechanistic Map](pv_qsp_model.png)](pv_qsp_model.svg)

### Subgraph Clusters (10 clusters, 100+ nodes)

| # | Cluster | Key Components |
|---|---------|---------------|
| 1 | **Genetic & Molecular Basis** | JAK2 V617F, TET2, ASXL1, SF3B1, IDH1/2, TP53, allele burden |
| 2 | **JAK-STAT Signaling** | JAK1, JAK2, STAT5A/B, STAT3, PI3K/AKT/mTOR, ERK1/2, SOCS1/3, BCL-XL |
| 3 | **Cytokine Receptors** | EPOR, MPL, G-CSFR, IL3R, IL6R/gp130, SCF/c-KIT, HIF-1α/2α, VHL |
| 4 | **BM Hematopoiesis** | HSC → MPP → CMP → MEP/GMP, BFU-E, CFU-E, erythroblasts, megakaryocytes |
| 5 | **Peripheral Blood** | RBC mass, Hgb, Hct, platelet pool, WBC, spleen sequestration |
| 6 | **Thrombosis** | Blood viscosity, platelet activation, NETs, coagulation, DVT/PE/stroke |
| 7 | **Ruxolitinib PK** | Absorption, central/peripheral distribution, CYP3A4 metabolism, metabolites |
| 8 | **HYU & IFN PK** | Hydroxyurea PK, PEG-IFN-α2a SC depot/central, phlebotomy |
| 9 | **Drug PD** | JAK inhibition (Emax), pSTAT5 reduction, ribonucleotide reductase inhibition |
| 10 | **Clinical Endpoints** | CHR, SVR35, MPN-SAF TSS, OS, MF/AML transformation, thrombosis risk |

---

## mrgsolve ODE Model

### Compartments (16 ODEs)

**PK Compartments:**
- `DEPOT_RUX` → `CENT_RUX` ⇌ `PERI_RUX` (Ruxolitinib, 2-compartment)
- `CENT_HYU` (Hydroxyurea, 1-compartment)
- `SC_IFN` → `CENT_IFN` (PEG-IFN-α2a, SC depot)

**Hematopoiesis ODEs:**
- `BFUE` → `CFUE` → `RETIC_BM` → `RETIC_C` → `RBC`
- `PLT` (platelet pool)
- `WBC` (white blood cell pool)
- `SPL` (spleen volume)
- `FIBRO` (BM fibrosis score, 0–3)
- `ALLELE` (JAK2 V617F allele burden, %)

### Drug Effect Models

```
Ruxolitinib JAK2 inhibition:
  INH_RUX = Imax × Cp^Hill / (IC50^Hill + Cp^Hill)
  IC50 = 0.86 ng/mL (JAK2), IC50_JAK1 = 1.01 ng/mL

Hydroxyurea ribonucleotide reductase inhibition:
  INH_HYU = Imax × CpμM / (IC50_μM + CpμM)
  IC50 = 150 μM, Imax = 0.85

PEG-IFN-α2a clonal suppression (allele burden):
  EFF_IFN = Emax × Cp / (EC50 + Cp)
  EC50 = 50 pg/mL, Emax = 0.70
```

### Treatment Scenarios

| # | Scenario | Key Result |
|---|----------|-----------|
| 1 | Untreated (natural history) | Hct rises progressively, allele burden expands |
| 2 | Phlebotomy + Aspirin | Hct controlled, allele burden unchanged |
| 3 | Hydroxyurea 500 mg/d | CHR achievable, allele burden minimally reduced |
| 4 | Ruxolitinib 10 mg BID (RESPONSE trial) | SVR35, pSTAT5 reduction, symptom improvement |
| 5 | PEG-IFN-α2a 45 μg/wk (PROUD-PV) | Allele burden reduction, fibrosis improvement |
| 6 | Ruxolitinib dose-response (5–20 mg BID) | Dose-dependent Hct reduction |

---

## Shiny Dashboard (7 Tabs)

1. **Overview** — Disease summary, pathophysiology diagram, model statistics
2. **Patient Profile** — Age/sex/allele burden inputs, risk stratification (ELN), treatment recommendation
3. **Pharmacokinetics** — Cp-time profiles, Cmax/Tmax/AUC, dose-response, JAK2 inhibition curve
4. **PD & Hematology** — Hct, PLT, WBC, erythroid progenitors, EPO, pSTAT5, CHR timeline, MPN-SAF TSS
5. **Clinical Endpoints** — Spleen volume (SVR35), thrombosis risk, BM fibrosis, MF/AML transformation
6. **Scenario Comparison** — Head-to-head: Untreated vs HYU vs Ruxolitinib vs PEG-IFN-α2a
7. **Biomarkers** — JAK2 allele burden, BM fibrosis grade, reticulocyte count, pSTAT5, annual thrombosis risk

---

## Key Calibration References

| Drug | Trial | Key Endpoint | Result |
|------|-------|-------------|--------|
| Ruxolitinib 10 mg BID | RESPONSE | SVR35 at week 32 | 38% vs 1% (BAT) |
| Ruxolitinib 10 mg BID | RESPONSE | Hct control at week 32 | 60% vs 19% (BAT) |
| PEG-IFN-α2a 45 μg/wk | PROUD-PV | CHR at week 52 | 43% vs 46% (HYU) |
| PEG-IFN-α2a | CONTINUUM | Molecular response (MR) | ~21% MR at 36 months |
| Hydroxyurea | ECLAP | Thrombotic events | Reduced vs control |
| Aspirin 81 mg/d | ECLAP | Thrombosis prevention | 60% risk reduction |

---

## Ruxolitinib PK Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| Oral bioavailability (F) | 95% | Shi et al. 2011 |
| ka (absorption) | 2.4 /h | Pop PK model |
| Vc (central volume) | 72 L | Pop PK model |
| Vp (peripheral volume) | 28 L | Pop PK model |
| CL (clearance) | 22 L/h | Pop PK model |
| IC50 JAK2 | 2.8 nM = 0.86 ng/mL | Biochemical assay |
| IC50 JAK1 | 3.3 nM = 1.01 ng/mL | Biochemical assay |
| t½ | ~3 hours | Multiple sources |

---

## Quick Start

```r
# Install dependencies
install.packages(c("mrgsolve", "dplyr", "ggplot2", "patchwork"))

# Run mrgsolve model
source("pv_mrgsolve_model.R")

# Launch Shiny app
shiny::runApp("pv_shiny_app.R")
```

---

## References

See [`pv_references.md`](pv_references.md) for 58 curated PubMed references covering:
- Disease epidemiology & molecular pathophysiology (refs 1–6)
- JAK2 V617F mutation biology (refs 7–13)
- JAK-STAT signaling & hematopoiesis (refs 14–18)
- Ruxolitinib clinical trials (refs 19–24)
- Hydroxyurea & interferon trials (refs 25–30)
- Thrombosis complications (refs 31–35)
- Disease progression & transformation (refs 36–38)
- PK/PD pharmacology (refs 39–44)
- QSP/PK-PD modeling methodology (refs 45–50)
- Clinical guidelines (refs 51–53)
- Natural history & biomarkers (refs 54–58)

---

*QSP Disease Model Library | Date: 2026-06-25 | Disease Category: Chronic / Myeloproliferative Neoplasm*
