# Sickle Cell Disease (SCD) — Quantitative Systems Pharmacology Model

[![Mechanistic Map](scd_qsp_model.png)](scd_qsp_model.svg)

> **Disease:** Sickle Cell Disease (HbSS, SCD) | **Category:** Hereditary Hemoglobinopathy  
> **Model version:** v1.0 | **Date:** 2026-06-21

---

## Disease Overview

Sickle Cell Disease is caused by a single nucleotide substitution (GAG→GTG, Glu→Val) in codon 6 of the **HBB gene**, producing **HbS** (α₂βS₂). Upon deoxygenation, HbS forms rigid 14-strand polymer fibers that distort red blood cells (RBCs) into a sickle shape. Repeated sickling leads to:

- **Hemolysis** (RBC lifespan 10–20 d vs. 120 d): cell-free Hb scavenges NO, elevating LDH, bilirubin
- **Vaso-occlusion**: P-selectin-mediated adhesion of sickle RBCs, neutrophils, and platelets causes microvascular stasis
- **End-organ damage**: pulmonary HTN, stroke, sickle nephropathy, AVN, retinopathy, priapism

Globally ~300,000 neonates are born with SCD annually; ~80% in sub-Saharan Africa.

---

## Model Files

| File | Description |
|------|-------------|
| [`scd_qsp_model.dot`](scd_qsp_model.dot) | Graphviz mechanistic map source |
| [`scd_qsp_model.svg`](scd_qsp_model.svg) | Vector mechanistic map |
| [`scd_qsp_model.png`](scd_qsp_model.png) | Rasterized map (150 dpi) |
| [`scd_mrgsolve_model.R`](scd_mrgsolve_model.R) | mrgsolve ODE model (24 compartments, 7 scenarios) |
| [`scd_shiny_app.R`](scd_shiny_app.R) | Shiny interactive dashboard (6 tabs) |
| [`scd_references.md`](scd_references.md) | 50 PubMed references |

---

## Mechanistic Map Summary

**Clusters (12 subgraphs):**

| # | Cluster | Key Nodes |
|---|---------|-----------|
| 1 | HbS Genetics & Molecular Pathophysiology | HBB gene → E6V mutation → HbS tetramer → deoxy-HbS → HbS polymer (nucleation, growth) |
| 2 | Erythropoiesis & RBC Biology | EPO/EPOR → BFU-E → CFU-E → RET → young RBC → dense RBC → sickle RBC |
| 3 | Hemolysis & Anemia | Intravascular/extravascular lysis → cell-free Hb → haptoglobin → free heme → HO-1 → bilirubin, LDH |
| 4 | Vaso-occlusion & Adhesion | P-selectin, E-selectin, VCAM-1, PSGL-1, α4β1 integrin, Lu/BCAM → WBC-RBC bridge → VOC rate |
| 5 | Inflammation & Oxidative Stress | TLR4 → NF-κB → TNF-α, IL-1β, IL-6, IL-8, NLRP3 → ROS cascade |
| 6 | NO Biology & Vasomotor Tone | eNOS → NO → sGC → cGMP → vasodilation; cell-free Hb scavenges NO |
| 7 | Coagulation & Thrombosis | TF → Xa → thrombin → fibrin → microthrombi; PS+ sickle RBC procoagulant |
| 8 | End-Organ Damage | ACS, stroke, sickle nephropathy, priapism, AVN, pulmonary HTN, retinopathy |
| 9 | HU PK | Oral absorption → plasma (μM) → RNR inhibition |
| 10 | Voxelotor PK | Oral → plasma → RBC binding (high affinity) |
| 11 | Crizanlizumab & L-Glu PK | IV 2-cmpt mAb PK; L-Glu amino acid pool |
| 12 | Drug PD & Clinical Endpoints | HbF↑, Hgb↑, VOC↓, annual crisis rate, QoL |

**Total nodes: ~150 | Total edges: ~170**

---

## mrgsolve ODE Model (24 Compartments)

### Compartments

| Module | Compartments |
|--------|-------------|
| Hydroxyurea PK | `HU_gut`, `HU_plasma` |
| Voxelotor PK | `VOX_plasma`, `VOX_RBC` |
| Crizanlizumab PK | `CRIZ_C`, `CRIZ_P` |
| L-Glutamine PK | `LG` |
| Erythropoiesis | `CFU_E`, `RET` |
| RBC dynamics | `RBC_S`, `RBC_N`, `HbF_frac`, `Hgb`, `free_Hb` |
| Hemolysis markers | `Haptoglobin`, `LDH`, `Bilirubin` |
| Vascular | `NO`, `P_selectin`, `VOC` |
| Oxidative/Metabolic | `NADH`, `Iron` |
| End-organ | `TRV`, `eGFR` |

### Treatment Scenarios

| # | Scenario | Key Mechanism |
|---|----------|---------------|
| 1 | No Treatment (baseline) | Untreated SCD dynamics |
| 2 | Hydroxyurea 20 mg/kg/d | RNR inhibition → HbF↑ |
| 3 | Voxelotor 1500 mg/d | HbS O₂ affinity↑ → sickling↓ |
| 4 | Crizanlizumab 5 mg/kg q4w | P-selectin blockade → VOC↓ |
| 5 | L-Glutamine 5g BID | NAD⁺ redox support → oxidative stress↓ |
| 6 | HU + Voxelotor | Additive: HbF↑ + sickling↓ |
| 7 | HU + VOX + CRIZ (Triple) | Maximum efficacy combination |

### Key PD Effects Modeled

| Drug | PD Model | Effect |
|------|----------|--------|
| Hydroxyurea | Emax Hill (EC₅₀=15 μM, Emax=18%) | HbF↑ (stress erythropoiesis) |
| Voxelotor | RBC occupancy → polymerization inhibition | Deoxy-HbS↓, Hgb↑ |
| Crizanlizumab | Imax model (EC₅₀=0.5 μg/mL) | P-selectin inhibition → VOC↓ 45% |
| L-Glutamine | Pool saturation | NADH↑, oxidative stress↓ |

---

## Shiny App (6 Tabs)

| Tab | Content |
|-----|---------|
| **1. Patient Profile** | Disease overview, QSP model structure, treatment MoA table, ValueBoxes |
| **2. Drug PK Profiles** | HU/VOX/CRIZ concentration-time plots; PK parameter table |
| **3. Hematological Response** | Hgb, HbF%, reticulocyte%, LDH trends; hematology summary table |
| **4. Vaso-occlusion & Biomarkers** | VOC rate, P-selectin, NO bioavailability, bilirubin |
| **5. Scenario Comparison** | Multi-scenario overlay with any endpoint; efficacy comparison table |
| **6. End-Organ Dashboard** | TRV (pulmonary HTN), eGFR (renal), ferritin, haptoglobin; KPI ValueBoxes |

---

## Key Clinical Parameters

| Parameter | Baseline (untreated HbSS) | Target (treated) |
|-----------|--------------------------|------------------|
| Hemoglobin | 6–9 g/dL | ≥9–11 g/dL |
| HbF% | 2–10% | ≥20% (with HU) |
| LDH | 400–700 U/L | <300 U/L |
| Annual VOC rate | 2–6 crises/year | <2 crises/year |
| Reticulocytes | 10–25% | 5–15% |
| TRV | >2.5 m/s (many patients) | <2.5 m/s |

---

## References

See [`scd_references.md`](scd_references.md) for 50 curated PubMed citations organized by:
- Pivotal clinical trials (MSH, HOPE, SUSTAIN, L-Glutamine phase III)
- HbS polymerization & RBC biology
- Hemolysis & vascular biology
- Vaso-occlusion & adhesion
- HbF & hydroxyurea mechanisms
- QSP/PK/PD modeling
- End-organ complications
- Gene therapy & emerging treatments
- Global burden & epidemiology

---

## Dependencies

```r
# R packages required
library(mrgsolve)   # ODE solver
library(dplyr)      # Data manipulation
library(ggplot2)    # Static visualization
library(tidyr)      # Data reshaping
library(shiny)      # Interactive dashboard
library(plotly)     # Interactive plots
library(DT)         # Interactive tables
library(shinydashboard) # Dashboard layout
library(scales)     # Color/axis scales
```

---

*Part of the QSP Disease Model Library — generated by Claude Code Routine (CCR)*
