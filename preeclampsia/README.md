# Preeclampsia QSP Model

[![Disease](https://img.shields.io/badge/Disease-Preeclampsia-purple)]()
[![Compartments](https://img.shields.io/badge/ODE_Compartments-20-blue)]()
[![Scenarios](https://img.shields.io/badge/Scenarios-6-green)]()
[![References](https://img.shields.io/badge/References-60-orange)]()

## Overview

**Preeclampsia (PE)** is a pregnancy complication affecting 2–8% of pregnancies worldwide, defined by new-onset hypertension (SBP ≥140 or DBP ≥90 mmHg) after 20 weeks gestation combined with proteinuria (≥300 mg/24h) and/or end-organ damage. It is a leading cause of maternal and perinatal morbidity and mortality globally.

This QSP model integrates:
- Placental dysfunction and inadequate spiral artery remodeling
- Angiogenic imbalance: sFlt-1↑, PlGF↓, sEng↑
- Endothelial dysfunction: NO↓, ET-1↑, ROS↑
- Multiorgan involvement: cardiovascular, renal, coagulation, neurological, hepatic
- Drug PK/PD: Aspirin (prophylaxis), Labetalol, Nifedipine, MgSO4

---

## Mechanistic Map

[![Preeclampsia QSP Map](pe_qsp_model.png)](pe_qsp_model.svg)

> **Figure:** Mechanistic map of preeclampsia pathophysiology and pharmacology.
> Click image to open interactive SVG. Generated with Graphviz (`neato` layout).
>
> **15 subgraph clusters** covering:
> Maternal Risk Factors · Placental Dysfunction · Angiogenic Imbalance ·
> Endothelial Dysfunction · Cardiovascular/BP · Renal · Coagulation/HELLP ·
> Neurological · Hepatic · Complement · Fetal · Aspirin PK/PD · Labetalol PK/PD ·
> Nifedipine PK/PD · Magnesium Sulfate PK/PD · Clinical Endpoints

---

## Model Structure

### ODE Compartments (20 total)

| # | Compartment | Description |
|---|-------------|-------------|
| 1 | `DEPOT_ASP` | Aspirin gut absorption depot |
| 2 | `ASPIRIN` | Aspirin central plasma (mg) |
| 3 | `SALICYLATE` | Salicylate (hydrolysis product) |
| 4 | `COX1_INH` | COX-1 inhibition state (0–1, irreversible) |
| 5 | `DEPOT_LAB` | Labetalol gut depot |
| 6 | `LABETALOL` | Labetalol central plasma |
| 7 | `DEPOT_NIF` | Nifedipine gut depot |
| 8 | `NIFEDIPINE` | Nifedipine central plasma |
| 9 | `MG_PLASMA` | Magnesium plasma pool (mmol) |
| 10 | `SFLT1` | sFlt-1 plasma (pg/mL) |
| 11 | `PLGF` | PlGF plasma (pg/mL) |
| 12 | `SENG` | Soluble endoglin (ng/mL) |
| 13 | `NO_EA` | Nitric oxide bioavailability (relative) |
| 14 | `ET1` | Endothelin-1 (pg/mL) |
| 15 | `ROS` | Reactive oxygen species index |
| 16 | `SBP` | Systolic blood pressure (mmHg) |
| 17 | `DBP` | Diastolic blood pressure (mmHg) |
| 18 | `GFR_C` | Glomerular filtration rate (mL/min/1.73m²) |
| 19 | `PROTEINURIA` | Proteinuria (mg/24h) |
| 20 | `PLATELET` | Platelet count (×10³/µL) |
| +2 | `LDH_MK`, `SEIZURE_RISK` | HELLP/neurological markers |

### Treatment Scenarios (6)

| Scenario | Treatment | Rationale |
|----------|-----------|-----------|
| 1 | No treatment | Natural PE progression |
| 2 | Aspirin 75 mg/day (12 wk→) | ASPRE trial: 62% reduction in early PE |
| 3 | Labetalol 200 mg BID (24 wk→) | CHIPS trial: α1+β blocker for BP control |
| 4 | Nifedipine 30 mg MR/day (24 wk→) | L-type Ca²⁺ blocker, CHIPS trial |
| 5 | MgSO4 4 g load + 1 g/h (30 wk→) | Magpie trial: 58% ↓ eclampsia |
| 6 | Combination (Aspirin + Labetalol + MgSO4) | Optimal multidrug approach |

---

## Key Clinical Thresholds

| Biomarker | Normal | PE Threshold | Severe PE |
|-----------|--------|-------------|-----------|
| SBP (mmHg) | <140 | ≥140 | ≥160 |
| DBP (mmHg) | <90 | ≥90 | ≥110 |
| sFlt-1/PlGF ratio | <38 | 38–85 | >85 |
| Proteinuria (mg/24h) | <300 | ≥300 | ≥5000 |
| Platelets (×10³/µL) | >150 | 100–150 | <100 (HELLP) |
| LDH (IU/L) | <600 | 600–800 | >800 (HELLP) |
| Mg²⁺ plasma (mmol/L) | 0.7–1.0 | 1.7–3.5 (therapeutic) | >3.5 (toxic) |

---

## Files

| File | Description |
|------|-------------|
| `pe_qsp_model.dot` | Graphviz source (15 clusters, 150+ nodes) |
| `pe_qsp_model.svg` | Scalable vector mechanistic map |
| `pe_qsp_model.png` | PNG thumbnail (150 dpi) |
| `pe_mrgsolve_model.R` | mrgsolve ODE model + 6 scenario simulations |
| `pe_shiny_app.R` | 8-tab Shiny interactive dashboard |
| `pe_references.md` | 60 curated references (PubMed links) |
| `README.md` | This file |

---

## Quick Start

```bash
# Render mechanistic map
neato -Tsvg pe_qsp_model.dot -o pe_qsp_model.svg
neato -Tpng -Gdpi=150 pe_qsp_model.dot -o pe_qsp_model.png
```

```r
# Run mrgsolve simulation
install.packages(c("mrgsolve","dplyr","ggplot2","tidyr","patchwork"))
source("pe_mrgsolve_model.R")

# Launch Shiny dashboard
install.packages(c("shiny","shinydashboard","DT","gridExtra"))
shiny::runApp("pe_shiny_app.R")
```

---

## Shiny Dashboard Tabs

1. **Patient Profile** — Parameters, clinical status at GA 36wk, threshold table
2. **Drug PK** — Aspirin, Labetalol, Nifedipine, MgSO4 plasma concentrations over gestation
3. **Angiogenic Balance** — sFlt-1, PlGF, sEng, sFlt-1/PlGF ratio trajectories
4. **Cardio-Renal** — BP (SBP/DBP), endothelial markers, GFR, proteinuria
5. **HELLP & Neuro** — Platelet count, LDH, seizure risk, Mg safety window
6. **Scenario Comparison** — All 6 treatment arms side-by-side
7. **Biomarker Panel** — Heatmap + data table across gestation
8. **About** — Model description, calibration sources, thresholds

---

## Key Calibration Sources

| Trial | Finding | Model Parameter |
|-------|---------|----------------|
| ASPRE (Rolnik 2017, *Lancet*) | Aspirin 62% ↓ early PE | `F_Asp`, `kBP_Asp` |
| Maynard 2003 (*J Clin Invest*) | sFlt-1 overexpression → PE phenotype | `kprod_sFlt1`, `k_sFlt1_rise` |
| Verlohren 2010 (*AJOG*) | sFlt-1/PlGF ratio cutoff 38 | Threshold encoding |
| CHIPS 2015 (*NEJM*) | Labetalol/Nifedipine BP control | `Emax_Lab`, `Emax_Nif` |
| Magpie 2002 (*Lancet*) | MgSO4 58% ↓ eclampsia | `Emax_NMDA`, seizure ODE |

---

## Pathophysiology Summary

```
Placental Hypoxia / Insufficient Trophoblast Invasion
       ↓
HIF-1α → ↑sFlt-1 secretion ──→ Sequestration of VEGF & PlGF
       ↓                              ↓
   ↑sEng ──────────────────→  ↓ VEGFR2 signaling
                                      ↓
              Endothelial Dysfunction (↓NO, ↑ET-1, ↑ROS)
                          ↓
         ┌─────────────────┴──────────────────┐
         ↓                                     ↓
    ↑SVR → Hypertension              Glomerular endotheliosis
    ↑TXA2 → Platelet aggregation     → Podocyte injury
         ↓                                     ↓
    HELLP Syndrome                        Proteinuria
    Eclampsia (if MgSO4 ↓)           ↓GFR → ↑Creatinine
```

---

*Generated by Claude Code Routine (CCR) — QSP Disease Model Library | 2026-06-25*
