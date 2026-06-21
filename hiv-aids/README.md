# HIV/AIDS QSP Model

**Disease**: HIV/AIDS (Human Immunodeficiency Virus Infection & Acquired Immunodeficiency Syndrome)  
**Category**: Infectious Disease — Retroviral  
**Date Added**: 2026-06-21  
**Model Abbreviation**: `hiv`

---

## Mechanistic Map

[![HIV QSP Mechanistic Map](hiv_qsp_model.png)](hiv_qsp_model.svg)

> 11 subgraph clusters · 128 nodes · cross-cluster regulation edges

---

## Overview

HIV/AIDS is caused by infection with HIV-1 or HIV-2, retroviruses that primarily infect
CD4⁺ T lymphocytes, leading to progressive immune deficiency. Without ART,
median time from infection to AIDS is ~10 years. With modern integrase-based ART
(TDF/FTC/DTG or TAF/FTC/BIC), viral suppression rates exceed 90% at 48 weeks and
life expectancy approaches that of the general population.

### Key Pathophysiology

| Process | Mechanism | Clinical Consequence |
|---------|-----------|----------------------|
| Viral Entry | gp120–CD4 binding; CCR5/CXCR4 coreceptor | Tropism; target for Maraviroc/Ibalizumab |
| Reverse Transcription | RT (p51/p66) converts RNA→dsDNA | Error rate 3×10⁻⁵/bp; quasispecies diversity |
| Integration | Integrase strand transfer into host chromosome | Proviral DNA (target for INSTIs) |
| Viral Production | ~10⁹ virions/day in untreated HIV | Immune exhaustion; CD4 decline |
| Latent Reservoir | Resting CD4 memory cells; t½ ~44 years | Major barrier to cure |
| CD4 Depletion | Apoptosis, CTL killing, direct viral cytopathic effect | CD4 nadir <200 → AIDS |
| Immune Activation | Microbial translocation; chronic IL-6/TNF-α | CV risk, neurocognition, cancer |

---

## Antiretroviral Drug Classes

| Class | Mechanism | Key Drugs | IC50 (nM) | Resistance Mutations |
|-------|-----------|-----------|-----------|----------------------|
| NRTI | RT chain termination (intracell. TP/DP) | TDF, TAF, FTC, 3TC | TFV-DP ~100, FTC-TP ~1000 | K65R, M184V |
| NNRTI | RT allosteric inhibition | EFV, NVP, RPV | EFV ~3000 | K103N, E138K |
| PI | Protease cleavage inhibition | DRV, ATV, LPV | DRV ~0.5 | V82A, I50V |
| INSTI | Integrase strand transfer inhibition | DTG, BIC, RAL | DTG ~2, BIC ~1.7 | Q148H, R263K |
| CCR5 Antag. | Coreceptor blockade | Maraviroc | — | V3 loop (R5→X4) |
| Fusion | gp41 HR2 block | Enfuvirtide | — | gp41 HR1 mutations |

---

## Recommended First-Line Regimens (WHO/DHHS 2022)

| Regimen | Brand | Advantages | Limitations |
|---------|-------|-----------|-------------|
| TAF/FTC/BIC | Biktarvy | High barrier, once daily, no food req. | Cost, weight gain |
| TDF/FTC/DTG | Triumeq equiv. | High barrier, generic available | Neural tube defects (DTG) |
| TDF/FTC/EFV | Atripla equiv. | Cheapest, WHO LMIC | Neuropsychiatric, K103N |
| CAB-LA + RPV-LA | Cabenuva | Monthly/bimonthly injectable | NNRTI resistance, cost |

---

## mrgsolve Model Structure

### Compartments (18)

| # | Compartment | Description | Units |
|---|-------------|-------------|-------|
| 1 | TDF_GUT | TDF gut depot | mg |
| 2 | TDF_PLASMA | TFV plasma | µg/mL |
| 3 | TFV_DP | TFV-DP intracellular | nM |
| 4 | FTC_PLASMA | FTC plasma | µg/mL |
| 5 | FTC_TP | FTC-TP intracellular | nM |
| 6 | DTG_PLASMA | DTG plasma | µg/mL |
| 7 | BIC_PLASMA | BIC plasma | µg/mL |
| 8 | EFV_PLASMA | EFV plasma | µg/mL |
| 9 | DRV_PLASMA | DRV plasma | µg/mL |
| 10 | T_CELL | Uninfected CD4⁺ T cells | cells/µL |
| 11 | I_CELL | Productively infected CD4⁺ | cells/µL |
| 12 | V_FREE | Free plasma virus | copies/mL ×10⁻³ |
| 13 | L_CELL | Latently infected cells | cells/µL |
| 14 | E_CELL | CD8⁺ CTL effectors | cells/µL |
| 15 | INFLAM | Systemic IL-6 (inflammation) | pg/mL |
| 16 | RESIST | Drug resistance score | 0–1 |
| 17 | VL_LOG | Smoothed log₁₀ VL | — |
| 18 | CD4_SMOOTH | Smoothed CD4 count | cells/µL |

### Core ODE System (Perelson 1996)

```
dT/dt  = s_T - d_T·T - β_eff·V·T·(1-k_lat) + r_lat·L
dI/dt  = β_eff·V·T·(1-k_lat) - δ_I·I - k_kill·E·I
dV/dt  = p_V·I·(1-η_PI) - c_V·V + p_V·r_lat·L·0.1
dL/dt  = β_eff·V·T·k_lat - r_lat·L - k_kill·E·L·0.01
dE/dt  = s_E + p_E·I/(0.5+I) - d_E·E

β_eff = β·(1-η_total)·(1 + 2R)

η_total = 1 - (1-η_RT)·(1-η_INSTI)·(1-η_PI)
η_RT    = 1 - (1-η_NRTI)·(1-η_EFV)
η_NRTI  = 1 - (1-η_TFV)·(1-η_FTC)
η_drug  = C^n / (IC50^n + C^n)   [Emax Hill model]
```

---

## Treatment Scenarios (8)

| # | Scenario | Regimen | ART Start | Special |
|---|----------|---------|-----------|---------|
| 1 | No ART | — | Never | Natural progression |
| 2 | TDF/FTC/DTG | Triumeq-type | Day 0 | WHO 2022 preferred |
| 3 | TAF/FTC/BIC | Biktarvy | Day 0 | Contemporary 1st-line |
| 4 | TDF/FTC/EFV | WHO LMIC | Day 0 | NNRTI-based |
| 5 | DRV/r+DTG | Salvage | Day 0 | PI+INSTI 2nd-line |
| 6 | Delayed ART | TDF/FTC/DTG | Day 180 | Late treatment start |
| 7 | STI | TDF/FTC/DTG | Day 0→365 off→730 | Treatment interruption |
| 8 | PrEP | TDF/FTC | Day 0 | Low VL exposure prevention |

---

## Shiny App Tabs

| Tab | Content |
|-----|---------|
| ① Patient Profile | CDC/WHO HIV staging, baseline parameters, regimen pharmacology |
| ② Viral Kinetics | Plasma VL dynamics (log₁₀), biphasic decline, suppression probability |
| ③ CD4 & Immunity | CD4 count trajectory, CTL dynamics, IL-6 inflammation |
| ④ Drug PK | Plasma PK (µg/mL), intracellular TFV-DP/FTC-TP, overall efficacy η |
| ⑤ Scenario Compare | Multi-regimen VL/CD4 comparison, week-48 summary table |
| ⑥ Reservoir & Resistance | Latent reservoir decay, resistance score, virologic failure risk |

---

## Key Clinical Trial Evidence

| Trial | Regimen | n | VL <50 Week 48 | CD4 Rise |
|-------|---------|---|----------------|----------|
| SPRING-2 (2013) | DTG vs RAL (+TDF/FTC) | 822 | 88% vs 85% | +188 vs +192 |
| SINGLE (2013) | DTG/ABC/3TC vs EFV/TDF/FTC | 833 | 88% vs 81% (p=0.003) | +267 vs +208 |
| GEMINI-1/2 (2019) | DTG+3TC vs DTG+TDF/FTC | 1433 | 93% vs 93% | +259 vs +251 |
| ATLAS (2020) | CAB-LA+RPV-LA q4w vs oral | 308 | 92% vs 95% | Maintained |
| ATLAS-2M (2020) | CAB-LA+RPV-LA q8w vs q4w | 1045 | 94% vs 93% | Maintained |
| iPrEx (2010) | TDF/FTC PrEP | 2499 | 44% HIV reduction | N/A (prevention) |
| START (2015) | Early vs deferred ART | 4685 | AIDS/serious illness ↓57% | +215 vs +57 |

---

## Files

| File | Description |
|------|-------------|
| `hiv_qsp_model.dot` | Graphviz DOT mechanistic map source (11 clusters, 128 nodes) |
| `hiv_qsp_model.svg` | SVG rendered mechanistic map |
| `hiv_qsp_model.png` | PNG rendered mechanistic map (150 dpi) |
| `hiv_mrgsolve_model.R` | mrgsolve ODE model (18 CMT, 8 scenarios) |
| `hiv_shiny_app.R` | Shiny dashboard (6 tabs, interactive PK/PD simulation) |
| `hiv_references.md` | 62 PubMed references (16 sections) |
| `README.md` | This file |
