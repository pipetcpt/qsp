# NAFLD/MASLD QSP Model

**Non-Alcoholic Fatty Liver Disease / Metabolic-Associated Steatotic Liver Disease**  
Quantitative Systems Pharmacology Model — Version 1.0

---

## Overview

This directory contains a comprehensive QSP model for NAFLD/MASLD, covering the full spectrum from simple hepatic steatosis through MASH (metabolic-associated steatohepatitis) to advanced fibrosis and cirrhosis.

### Disease Summary

| Feature | Description |
|---------|-------------|
| **Global Prevalence** | ~25–30% of adults worldwide (~2 billion people) |
| **MASH Prevalence** | ~3–5% of general population; ~20% of NAFLD |
| **Fibrosis F3–F4** | ~5–10% of MASH patients; highest mortality risk |
| **Key Risk Factors** | T2DM, obesity, dyslipidemia, metabolic syndrome |
| **Primary Mechanism** | Insulin resistance → hepatic lipid accumulation → oxidative stress → inflammation → fibrosis |
| **FDA-Approved Tx** | Resmetirom (Rezdiffra™) — March 2024 |

---

## Mechanistic Map

[![NAFLD/MASLD QSP Mechanistic Map](nafld_qsp_model.png)](nafld_qsp_model.svg)

*Click image to view full-resolution SVG. The map covers 10 biological subsystems with 100+ nodes.*

### Subsystems Modeled

| # | Subsystem | Key Components |
|---|-----------|----------------|
| 1 | Adipose Tissue & Insulin Resistance | HSL/ATGL lipolysis, adipokines (adiponectin, leptin, resistin), ceramide |
| 2 | Pancreas & Glucose Homeostasis | β-cell function, GLP-1, glucagon, HOMA-IR |
| 3 | Hepatic Lipid Metabolism | DNL (SREBP-1c/ChREBP/ACC/FAS), FFA uptake, β-oxidation, VLDL secretion |
| 4 | Mitochondrial & Oxidative Stress | ETC, ROS, Nrf2/Keap1, lipid peroxidation (4-HNE, MDA), ER stress (UPR) |
| 5 | Hepatic Inflammation (MASH) | Kupffer cells, NF-κB, NLRP3, TNF-α, IL-1β, IL-6, neutrophil/monocyte recruitment |
| 6 | Hepatocyte Death | Apoptosis (caspase-3/8), necroptosis (RIPK3/MLKL), ballooning, CK-18 |
| 7 | Stellate Cell & Fibrogenesis | TGF-β1/SMAD2-3, PDGF, collagen I/III/IV, TIMP/MMP, LOXL2, YAP/TAZ |
| 8 | Gut-Liver Axis | Dysbiosis, LPS/TLR4, FXR/FGF-19, bile acid cycling, SCFA, TMAO |
| 9 | Drug PK/PD | FXR agonists, GLP-1 RAs, THRβ agonists, PPARα/δ, ACC inhibitors |
| 10 | Clinical Endpoints | NAS score, fibrosis stage, liver stiffness, ALT/AST, FIB-4, ELF score |

---

## Files

| File | Description |
|------|-------------|
| [`nafld_qsp_model.dot`](nafld_qsp_model.dot) | Graphviz mechanistic map source (100+ nodes, 10 subgraphs) |
| [`nafld_qsp_model.svg`](nafld_qsp_model.svg) | Vector mechanistic map (scalable) |
| [`nafld_qsp_model.png`](nafld_qsp_model.png) | Raster map at 150 dpi |
| [`nafld_mrgsolve_model.R`](nafld_mrgsolve_model.R) | mrgsolve ODE model (22 compartments, 5 scenarios) |
| [`nafld_shiny_app.R`](nafld_shiny_app.R) | Shiny dashboard (6 tabs) |
| [`nafld_references.md`](nafld_references.md) | 59 PubMed-linked references |

---

## mrgsolve Model

### Compartments (22 total)

**Drug PK (5):** `GUT`, `CENT`, `PERI` (FXR agonist 2-cmt), `GUT_GLP1`, `CENT_GLP1` (GLP-1 RA 1-cmt)

**Hepatic Lipid (3):** `LFFA` (free fatty acids), `LTAG` (triglycerides), `LDAG` (diacylglycerol)

**Oxidative/ER Stress (3):** `ROS_LVR`, `NRF2_ACT`, `ER_STRESS`

**Inflammation (6):** `KUP_ACT`, `TNF`, `IL6`, `IL1B`, `MCP1`, `NEUTRO`

**Cell Death (1):** `HEPATO_APOP`

**Fibrosis (3):** `HSC_ACTIV`, `TGF_B1`, `COLLAGEN`

### Treatment Scenarios

| Scenario | Drug | Dose | Route | Calibration |
|----------|------|------|-------|-------------|
| 1 | No treatment | — | — | Natural history baseline |
| 2 | Obeticholic acid (OCA) | 25 mg/day | PO | REGENERATE trial: F↓≥1 in 23% vs 12% |
| 3 | Semaglutide | 2.4 mg/wk | SC | NATIVE trial: NASH resolution 59% vs 17% |
| 4 | OCA + Semaglutide | Combination | PO + SC | SYNERGY projected |
| 5 | Resmetirom | 80 mg/day | PO | MAESTRO-NASH: NAS↓≥2 in 25.9%; F↓≥1 in 24.2% |

### Key Equations

```
dLTAG/dt = 0.15·LFFA           (esterification)
           - kVLDL_sec·LTAG     (VLDL export)
           - kLTAG_deg·LTAG     (lipolysis)

dHSC_ACTIV/dt = kHSC_act·TGF_B1·(1 − HSC)     (activation)
              + 0.002·HEPATO_APOP·(1 − HSC)   (apoptotic body signal)
              − kHSC_res·(1 + adiponectin)·HSC (reversion)

dCOLLAGEN/dt = kCol_prod·HSC·TGF_B1            (synthesis)
              − kCol_deg·(1 − 0.6·HSC)·COL     (MMP-mediated degradation)
```

---

## Shiny App

### Tabs

| Tab | Content |
|-----|---------|
| **Patient Profile** | BMI, HbA1c, adipokines, baseline NAS, mechanistic map image |
| **PK Profiles** | Plasma concentration–time curves for FXR agonist & GLP-1 RA |
| **PD Biomarkers** | Multi-marker time courses (TG, ROS, Kupffer, cytokines, collagen) |
| **Clinical Endpoints** | NAS score, fibrosis stage, ALT/AST, FIB-4, ELF score |
| **Scenario Comparison** | Side-by-side comparison of all 5 treatment scenarios |
| **Biomarker Panel** | Real-time value boxes, heatmap, waterfall plot, trial comparison table |

---

## Clinical Trial Calibration

| Trial | Intervention | Duration | Primary Endpoint | Response Rate |
|-------|-------------|----------|-----------------|---------------|
| MAESTRO-NASH | Resmetirom 80mg | 52 wk | NAS↓≥2 + F stable | **25.9%** vs 14.2% |
| MAESTRO-NASH | Resmetirom 100mg | 52 wk | NAS↓≥2 + F stable | **29.9%** vs 14.2% |
| REGENERATE | OCA 25mg | 18 mo | F↓≥1 no NAS worsen | **23%** vs 12% |
| CENTAUR | Cenicriviroc | 52 wk | F↓≥1 no NAS worsen | **20%** vs 10% |
| LEAN | Liraglutide 1.8mg | 48 wk | NAS↓≥2 + F stable | **39%** vs 9% |
| NATIVE | Semaglutide 2.4mg | 72 wk | NASH resolution | **59%** vs 17% |

---

## Running the Model

### Prerequisites

```r
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr",
                   "shiny", "shinydashboard", "shinyWidgets",
                   "plotly", "DT"))
```

### Quick Start

```r
# Run mrgsolve ODE simulation
source("nafld_mrgsolve_model.R")

# Launch Shiny dashboard
shiny::runApp("nafld_shiny_app.R")
```

---

## References

See [`nafld_references.md`](nafld_references.md) for 59 PubMed-linked citations covering:
- Disease pathophysiology & nomenclature
- Hepatic lipid metabolism & DNL
- Insulin resistance & adipokines
- Oxidative & ER stress
- Hepatic inflammation & fibrosis
- All major clinical trials
- QSP modeling methods
- Epidemiology & natural history

---

*Model created: 2026-06-24 | Category: Chronic Disease — Hepatology*
