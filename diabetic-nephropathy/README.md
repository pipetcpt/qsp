# Diabetic Nephropathy (DN) — QSP Model

> **당뇨병성 신병증 (Diabetic Nephropathy)** | ICD-10: N08, E11.21  
> 제2형 당뇨병의 가장 흔한 미세혈관 합병증으로, 전 세계 말기 신부전(ESKD)의 주요 원인

---

## Overview

Diabetic Nephropathy (DN) affects ~40% of patients with type 2 diabetes and is the leading cause of end-stage kidney disease (ESKD) worldwide. The pathophysiology involves multiple converging pathways:

| Pathway | Key Mediators | Outcome |
|---------|--------------|---------|
| **Hemodynamic** | AngII, RAAS, intraglomerular pressure ↑ | GBM damage, hyperfiltration |
| **Metabolic** | AGE, PKC-β, polyol, hexosamine | Oxidative stress, endothelial dysfunction |
| **Fibrotic** | TGF-β1, Smad2/3, CTGF, ECM | Glomerulosclerosis, tubulointerstitial fibrosis |
| **Inflammatory** | NF-κB, TNF-α, IL-1β, MCP-1, NLRP3 | Macrophage infiltration, cytokine amplification |
| **Podocyte** | Nephrin, Podocin, slit diaphragm disruption | Albuminuria, proteinuria |
| **Tubular** | SGLT2, tubular hypoxia, EMT | GFR decline, interstitial fibrosis |

---

## Mechanistic Map

[![DN QSP Mechanistic Map](dn_qsp_model.png)](dn_qsp_model.svg)

*Click image to open full SVG. 100+ nodes, 9 subclusters, complete drug PK/PD pathways.*

**Subclusters in the map:**
1. Hyperglycemia & Glucose Metabolism
2. AGE/RAGE Axis & Oxidative Stress
3. Renal Hemodynamics & RAAS
4. TGF-β Signaling & Glomerulosclerosis
5. Renal Inflammation
6. Podocyte Pathobiology & Glomerular Barrier
7. Tubular Compartment & CKD Progression
8. Drug PK/PD Mechanisms (4 drug classes)
9. Clinical Endpoints & Biomarkers

---

## Files

| File | Description |
|------|-------------|
| [`dn_qsp_model.dot`](dn_qsp_model.dot) | Graphviz mechanistic map source |
| [`dn_qsp_model.svg`](dn_qsp_model.svg) | Vector graphic (interactive) |
| [`dn_qsp_model.png`](dn_qsp_model.png) | Rasterized thumbnail (150 dpi) |
| [`dn_mrgsolve_model.R`](dn_mrgsolve_model.R) | mrgsolve ODE model + simulation |
| [`dn_shiny_app.R`](dn_shiny_app.R) | Interactive Shiny dashboard |
| [`dn_references.md`](dn_references.md) | 45 PubMed-linked references |

---

## mrgsolve Model Specifications

### Compartments (19 ODEs)

| # | Compartment | Units | Description |
|---|-------------|-------|-------------|
| 1–2 | GI_acei, CENT_acei | mg, mg | ACEi (enalapril) PK |
| 3–4 | GI_arb, CENT_arb | mg, mg | ARB (losartan) PK |
| 5–6 | GI_sglt2, CENT_sglt2 | mg, mg | SGLT2i (empagliflozin) PK |
| 7–8 | GI_fine, CENT_fine | mg, mg | Finerenone PK |
| 9 | BG | mmol/L | Blood glucose |
| 10 | AGE_cmpt | AU | AGE accumulation |
| 11 | AngII_cmpt | AU | Angiotensin II (RAAS activity) |
| 12 | TGF_cmpt | AU | TGF-β (fibrosis driver) |
| 13 | ROS_cmpt | AU | Oxidative stress |
| 14 | ECM_cmpt | AU | Extracellular matrix / glomerulosclerosis |
| 15 | Pod_cmpt | AU | Podocyte integrity (1=intact) |
| 16 | UACR_cmpt | mg/g | Urinary albumin-creatinine ratio |
| 17 | Tub_cmpt | AU | Tubular integrity |
| 18 | Fib_cmpt | AU | Interstitial fibrosis |
| 19 | GFR_cmpt | mL/min | Estimated GFR |

### Drug PK (Emax PD)

| Drug | Model | Emax | EC50 | Target |
|------|-------|------|------|--------|
| ACEi (Enalapril 10mg BID) | 1st-order oral, 1-cpt | 90% | 2.5 ng/mL | ACE → ↓AngII |
| ARB (Losartan 100mg QD) | 1st-order oral, 1-cpt | 85% | 80 ng/mL | AT1R blockade |
| SGLT2i (Empa 25mg QD) | 1st-order oral, 1-cpt | 85% | 15 ng/mL | SGLT2 → glucosuria, ↓tubular O₂ |
| Finerenone (20mg QD) | 1st-order oral, 1-cpt | 88% | 120 ng/mL | MR → ↓fibrosis, ↓inflammation |

### Treatment Scenarios

| ID | Regimen | Inspired by |
|----|---------|-------------|
| S0 | No treatment (natural history) | — |
| S1 | ACEi monotherapy | Lewis 1993 NEJM |
| S2 | ARB monotherapy | RENAAL / IDNT 2001 |
| S3 | SGLT2i monotherapy | CREDENCE / DAPA-CKD |
| S4 | ACEi + SGLT2i | EMPA-REG + CREDENCE combined analysis |
| S5 | SGLT2i + Finerenone | CONFIDENCE trial 2023 |
| S6 | ACEi + SGLT2i + Finerenone | Emerging triple combination |

---

## Key Calibration References

| Endpoint | Target | Reference |
|----------|--------|-----------|
| Natural GFR decline | −3 to −4 mL/min/yr | Perkins 2003, NEJM |
| ACEi UACR reduction | ~30–35% | Lewis 1993, NEJM |
| ARB UACR reduction | ~25–30% | RENAAL/IDNT 2001 |
| SGLT2i GFR slope | +2 mL/min/yr vs placebo | CREDENCE 2019 |
| Finerenone UACR | −31% | FIDELIO-DKD 2020 |
| Combo SGLT2i+Fine | Greater UACR reduction | CONFIDENCE 2023 |

---

## Shiny App Features (6 Tabs)

| Tab | Content |
|-----|---------|
| **1. Patient Profile** | Age, sex, DM duration, baseline eGFR/UACR/HbA1c, KDIGO CKD staging, drug dose inputs |
| **2. Drug PK** | Plasma concentration-time profiles, steady-state PK table, % inhibition by drug |
| **3. PD / Biomarkers** | TGF-β, ROS, ECM, Podocyte integrity, Interstitial fibrosis, AngII over time |
| **4. Clinical Endpoints** | eGFR trajectory, UACR, SBP, HbA1c, CKD stage progression |
| **5. Scenario Comparison** | Forest plots for eGFR and UACR% change, downloadable summary table |
| **6. GFR Slope & ESKD** | Annualized GFR slope bar chart, ESKD threshold timeline |

### Running the Shiny App

```r
# Install dependencies
install.packages(c("shiny","shinydashboard","mrgsolve","dplyr",
                   "ggplot2","tidyr","DT","plotly","scales"))

# Launch
shiny::runApp("dn_shiny_app.R")
```

---

## Running the mrgsolve Simulation

```r
# Install mrgsolve
install.packages(c("mrgsolve", "dplyr", "ggplot2", "patchwork"))

# Source and run
source("dn_mrgsolve_model.R")
# Produces:
# - 2-year multi-scenario simulation (S0-S6)
# - Combined 6-panel ggplot
# - GFR slope table
# - Summary outcome table at Year 2
```

---

## Disease Highlights

- **Prevalence**: ~40% of T2DM patients develop DN; ~100,000 new ESKD cases/year in the US attributable to DN
- **CKD Classification**: KDIGO G1-G5 × A1-A3 (eGFR × UACR matrix)
- **Breakthrough**: SGLT2 inhibitors reduce ESKD risk by ~30-40% independent of glycemic control
- **Novel MRA**: Finerenone (non-steroidal) provides renal and CV protection with lower hyperkalemia risk than spironolactone
- **Biomarker milestones**: UACR >300 mg/g = high risk; eGFR <30 = CKD G4; eGFR <15 = ESKD

---

*Model built by Claude Code Routine (CCR) | 2026-06-24*
