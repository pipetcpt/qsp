# Alcoholic Liver Disease (ALD) — QSP Model

[![QSP Map](ald_qsp_model.png)](ald_qsp_model.svg)

## Overview

This module provides a comprehensive **Quantitative Systems Pharmacology (QSP)** model for **Alcoholic Liver Disease (ALD)**, a leading cause of liver-related mortality worldwide, encompassing a spectrum from simple hepatic steatosis through alcoholic hepatitis (AH) to cirrhosis and hepatocellular carcinoma.

The model integrates:
- **Ethanol Metabolism**: ADH, CYP2E1 (MEOS), ALDH2, acetaldehyde kinetics
- **Oxidative Stress**: ROS/GSH balance, Nrf2/KEAP1/ARE pathway, ferroptosis, lipid peroxidation
- **Gut–Liver Axis**: Dysbiosis, intestinal permeability, LPS/TLR4, NLRP3 inflammasome
- **Kupffer Cell Activation**: NF-κB, AP-1, IRF3, NLRP3/IL-1β/gasdermin D/pyroptosis
- **Neutrophil Infiltration**: CXCL8/CXCR2, NETs, MPO, elastase, tissue necrosis
- **Hepatocyte Death & Regeneration**: Apoptosis, necroptosis, ferroptosis, ER stress, HGF/c-Met
- **Liver Fibrosis**: HSC activation, TGF-β1/SMAD, PDGF-BB, ECM remodeling, MMP/TIMP balance
- **Drug PK/PD**: Prednisolone (2-cmt), N-acetylcysteine (IV), G-CSF (SC), pentoxifylline, anakinra
- **Clinical Endpoints**: MELD, Maddrey's DF, ABIC score, Lille score, 90-day mortality

---

## Disease Background

### Epidemiology
- Global alcohol-attributable deaths: ~3 million/year (5.3% of all deaths)
- ALD accounts for ~50% of liver disease-related mortality in high-income countries
- Severe alcoholic hepatitis (AH) carries 30–50% 90-day mortality without treatment (Thursz 2015)
- ALD now leads as indication for liver transplantation in the USA, surpassing HCV (2018)

### Disease Spectrum
| Stage | Description | Typical MELD | Reversibility |
|-------|-------------|--------------|---------------|
| Steatosis | Simple fat accumulation | 6–8 | Yes (abstinence) |
| Steatohepatitis | Inflammation + ballooning | 10–18 | Partial |
| Alcoholic Hepatitis | Acute inflammatory injury | 18–35 | Partial |
| Fibrosis (F1–F3) | Progressive scarring | 12–22 | Partial |
| Cirrhosis (F4) | Irreversible scarring | 18–30 | No |
| ACLF | Organ failure on cirrhosis | 25–40 | Low |

### Pathophysiology
1. **Ethanol metabolism**: ADH + CYP2E1 → acetaldehyde (toxic) + ROS; ALDH2 → acetate
2. **Oxidative stress**: CYP2E1 generates superoxide; depletes GSH; activates Nrf2/ARE defense
3. **Gut dysbiosis**: Ethanol disrupts tight junctions (ZO-1, occludin) → LPS translocation
4. **Kupffer cell priming**: TLR4/LPS → NF-κB; NLRP3 inflammasome → IL-1β/pyroptosis
5. **Neutrophil storm**: CXCL8 gradient → hepatic PMN infiltration → ROS, MPO, NETs, necrosis
6. **Hepatocyte death**: TNF-TNFR1 → apoptosis/necroptosis; ferroptosis (GPX4↓); ER stress
7. **Stellate cell activation**: TGF-β1/SMAD → collagen I/III deposition → fibrosis

---

## Model Files

| File | Description |
|------|-------------|
| `ald_qsp_model.dot` | Graphviz source — 200+ nodes, 10 subgraph clusters |
| `ald_qsp_model.svg` | High-resolution vector map |
| `ald_qsp_model.png` | Raster image (150 DPI) |
| `ald_mrgsolve_model.R` | ODE-based PK/PD model (22 compartments, 7 scenarios) |
| `ald_shiny_app.R` | Interactive Shiny dashboard (7 tabs) |
| `ald_references.md` | 61 curated PubMed references |

---

## Mechanistic Map Summary

10 subgraph clusters covering 200+ nodes:

| # | Cluster | Key Components |
|---|---------|----------------|
| 1 | Ethanol Metabolism | ADH1B/1C/4, CYP2E1, catalase, ALDH2, acetaldehyde, NADH/NAD+ redox shift |
| 2 | Oxidative Stress | ROS, GSH/GSSG, GPX1/4, SOD1/2, Nrf2/KEAP1/ARE, 4-HNE, MDA, ferroptosis |
| 3 | Gut–Liver Axis | Dysbiosis, LPS/LTA/flagellin, TLR4/2/5, MD-2/CD14/LBP, MyD88/TRIF, SCFA↓ |
| 4 | Kupffer Cell | NF-κB, IKK, NLRP3/ASC/caspase-1/GSDMD, IL-1β, TNF-α, MCP-1, CXCL1/8 |
| 5 | Neutrophil | CXCR2, NETs, MPO/HOCl, elastase, platelet aggregation, microthrombus |
| 6 | Hepatocyte | Apoptosis (casp-3/8/9), necroptosis (RIPK3/MLKL), ER stress, HGF/c-Met, STAT3 |
| 7 | Fibrosis/HSC | TGF-β1/SMAD2-3/SMAD7, PDGF-BB, collagen I/III, MMP2/9, TIMP1/2 |
| 8 | Drug PK | Prednisolone 2-cmt oral, NAC IV, G-CSF SC, pentoxifylline, anakinra SC |
| 9 | Drug PD | GRα/GRE/IκBα, NAC→GSH, G-CSF→granulopoiesis, IL-1Ra, PTX→cAMP |
| 10 | Endpoints | MELD, ALT/AST, bilirubin, INR, Maddrey's DF, Lille score, 90-day mortality |

---

## mrgsolve Model — State Variables (22 ODEs)

| Variable | Description | Initial (Severe AH) |
|----------|-------------|---------------------|
| `ETOH` | Blood ethanol (mg/dL) | 60 (active) / 0 |
| `AA` | Acetaldehyde (μM) | 2.0 |
| `ROS` | Reactive oxygen species (rel.) | 3.5 |
| `GSH` | Glutathione (mM) | 2.0 (depleted) |
| `LPS` | Portal LPS (rel.) | 4.0 |
| `KC` | Kupffer cell activation (rel.) | 3.5 |
| `TNF` | TNF-α (rel.) | 4.0 |
| `IL1B` | IL-1β (rel.) | 5.0 |
| `NEUT` | Liver neutrophils (rel.) | 3.0 |
| `H` | Hepatocyte viability (fraction) | 0.55 |
| `ALT` | Serum ALT (IU/L) | 180 |
| `BILI` | Bilirubin (mg/dL) | 10.0 |
| `INR` | INR | 1.8 |
| `F` | Fibrosis (Laennec 0–4) | 1.5 |
| `PRED_gut` | Prednisolone gut compartment | 0 |
| `PRED_C` | Prednisolone central (ng/mL) | 0 |
| `PRED_P` | Prednisolone peripheral | 0 |
| `NAC_C` | NAC plasma (μg/mL) | 0 |
| `GCSF_C` | G-CSF plasma (ng/mL) | 0 |
| `PTX_C` | Pentoxifylline plasma (ng/mL) | 0 |
| `ANK_C` | Anakinra plasma (ng/mL) | 0 |

---

## Treatment Scenarios

| # | Scenario | Regimen | Key Effect |
|---|----------|---------|------------|
| S1 | Active Drinking | No treatment, continued ethanol | Progressive liver injury |
| S2 | Abstinence | Ethanol cessation only | Partial spontaneous recovery |
| S3 | Prednisolone | 40 mg/day × 28 days | NF-κB inhibition; improves 28-day survival (STOPAH) |
| S4 | NAC IV | 150→50→100 mg/kg GET protocol | GSH replenishment; ROS reduction |
| S5 | Pred + NAC | Prednisolone + NAC combined | Best short-term survival in some cohorts |
| S6 | G-CSF | 5 μg/kg SC × 5 days | Granulopoiesis; hepatic regeneration; HGF↑ |
| S7 | Pred + Anakinra | Prednisolone + IL-1Ra 100 mg SC/d | Combined NF-κB + NLRP3 blockade (investigational) |

---

## Shiny App Features

| Tab | Content |
|-----|---------|
| 1 | Patient Profile: demographics, ALT, bilirubin, INR, fibrosis; MELD, Maddrey's DF, ABIC scores |
| 2 | Drug PK: prednisolone 2-cmt, NAC IV GET protocol, G-CSF SC kinetics |
| 3 | Disease Dynamics: ROS/GSH, LPS/Kupffer, TNF/IL-1β, neutrophil infiltration |
| 4 | Clinical Endpoints: MELD, ALT, bilirubin, INR, hepatocyte viability, fibrosis |
| 5 | Scenario Comparison: all 7 regimens head-to-head with summary table |
| 6 | Biomarkers & Risk: Lille score calculator, MELD/DF calculator, GSH sensitivity analysis, steroid responder stratification |

---

## Key Clinical Insights from Model

### Why Prednisolone Works — and Fails
Prednisolone suppresses NF-κB → reduces TNF-α, IL-1β, CXCL8 → decreases neutrophil recruitment → reduces hepatocyte necrosis. However, ~45% of patients are Lille non-responders (Lille score ≥0.45 at Day 7), predicted by persistent bilirubin elevation reflecting ongoing hepatocyte injury despite cytokine suppression.

### NAC's Complementary Role
GSH depletion is central to ALD pathogenesis (CYP2E1 generates ROS, acetaldehyde consumes GSH). NAC replenishes cysteine → GCLC-mediated GSH synthesis → ROS scavenging. The STOPAH trial showed no 90-day benefit for prednisolone alone vs. placebo, while the GET protocol combining prednisolone + NAC showed improved 6-month survival in the Nguyen-Khac study.

### G-CSF's Regenerative Mechanism
Beyond neutrophil mobilization, G-CSF stimulates CD34+ progenitor cells, which differentiate into hepatocytes (via HGF/IL-6 signaling) — explaining its hepatoregenerative benefits in ACLF independent of its anti-infective role.

### NLRP3 as a Therapeutic Target
Acetaldehyde and ROS provide both signal 1 (NF-κB priming via TLR4) and signal 2 (NLRP3 activation) for IL-1β maturation. Blocking IL-1R1 with anakinra interrupts IL-1β-mediated HSC activation, neutrophil recruitment, and hepatocyte apoptosis — rationale for the PILOT trial.

---

## Parameter Sources

| Parameter | Value | Source |
|-----------|-------|--------|
| CYP2E1 Km (ethanol) | 50 mg/dL | Lieber 2005 |
| Acetaldehyde half-life | ~20 min | Tuma 2003 |
| GSH hepatic baseline | 5 mM | Dey 2006 |
| Prednisolone CL | 8.5 L/h | Bergrem 1983 |
| NAC IV CL | 12 L/h | Borgström 1986 |
| G-CSF Vc | 4.5 L | Kuwabara 1994 |
| Fibrosis progression rate | 0.0003/h | Tsuchida 2017 |
| 90-day mortality MELD model | logistic | Kim 2003 |

---

## References

See [`ald_references.md`](ald_references.md) for 61 curated PubMed references:
- §1–2: Epidemiology + ethanol metabolism
- §3–4: Oxidative stress + gut–liver axis
- §5–7: Innate immunity + neutrophils + hepatocyte death
- §8–9: Fibrosis biology + prednisolone/NAC trials (STOPAH, GET)
- §10–11: G-CSF, anakinra, emerging therapies
- §12–14: Scoring systems + transplantation + drug PK

---

## How to Run

### mrgsolve Model
```r
install.packages(c("mrgsolve","dplyr","ggplot2","tidyr","patchwork"))
source("ald_mrgsolve_model.R")
```

### Shiny Dashboard
```r
install.packages(c("shiny","shinydashboard","mrgsolve","dplyr","ggplot2","tidyr","DT","plotly"))
shiny::runApp("ald_shiny_app.R")
```

### Regenerate Mechanistic Map
```bash
dot -Tsvg ald_qsp_model.dot -o ald_qsp_model.svg
dot -Tpng -Gdpi=150 ald_qsp_model.dot -o ald_qsp_model.png
```

---

*Model developed for the QSP Disease Model Library | Date: 2026-06-26*  
*Category: Hepatology / Gastroenterology / Alcohol-related liver disease*
