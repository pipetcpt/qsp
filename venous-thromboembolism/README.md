# Venous Thromboembolism (VTE) QSP Model

> **Deep Vein Thrombosis (DVT) · Pulmonary Embolism (PE)** — Comprehensive Quantitative Systems Pharmacology Model

[![Model](https://img.shields.io/badge/QSP-VTE%20%7C%20DVT%20%7C%20PE-blue)]()
[![ODEs](https://img.shields.io/badge/ODEs-19%20compartments-green)]()
[![Drugs](https://img.shields.io/badge/Drugs-5%20anticoagulants-orange)]()
[![Shiny](https://img.shields.io/badge/Shiny-6%20tabs-purple)]()

---

## Overview

Venous thromboembolism (VTE) — encompassing deep vein thrombosis (DVT) and pulmonary embolism (PE) — is the **third most common cardiovascular disease** after myocardial infarction and stroke, affecting 1-2 per 1,000 adults annually. PE carries a 30-day mortality of 1–3% (up to 30% for massive PE) and 2–4% of survivors develop chronic thromboembolic pulmonary hypertension (CTEPH).

This QSP model integrates:
- **Virchow's Triad** (stasis / endothelial injury / hypercoagulability)
- **Full coagulation cascade** (extrinsic + intrinsic + common pathways)
- **Platelet activation** (GPIb, GPVI, PAR1/4, GPIIb/IIIa)
- **Natural anticoagulants** (AT-III, Protein C/S, TFPI)
- **Fibrinolysis** (tPA/uPA-plasmin-D-dimer axis, PAI-1, TAFI)
- **Multi-drug PK/PD**: Rivaroxaban · Apixaban · Dabigatran · Warfarin · Enoxaparin

---

## Files

| File | Description |
|------|-------------|
| `vte_qsp_model.dot` | Graphviz mechanistic map (12 clusters, 140+ nodes) |
| `vte_qsp_model.svg` | SVG (interactive/vector) |
| `vte_qsp_model.png` | PNG 150 dpi |
| `vte_mrgsolve_model.R` | mrgsolve ODE model (19 compartments, 6 scenarios) |
| `vte_shiny_app.R` | Shiny dashboard (6 tabs) |
| `vte_references.md` | 57 PubMed references |

---

## Mechanistic Map Structure

The `.dot` file contains **12 subgraph clusters**:

| # | Cluster | Key Nodes |
|---|---------|-----------|
| ① | Virchow's Triad & Risk Factors | Stasis, Endothelial Injury, Hypercoagulability; FV Leiden, PT G20210A, APS, OCP |
| ② | Vascular Endothelium | TF expression, vWF release, PGI2/NO, TM, EPCR, tPA, PAI-1 |
| ③ | Extrinsic Pathway (TF) | TF exposed, FVII/FVIIa, TF-FVIIa tenase, TFPI inhibition |
| ④ | Intrinsic Pathway (Contact) | FXII → kallikrein → FXIa → FIXa, FVIII/FVIIIa, Intrinsic Tenase |
| ⑤ | Common Pathway / Thrombin | FXa, Prothrombinase, Thrombin (FIIa), Fibrinogen → Fibrin cross-linked |
| ⑥ | Platelet Activation | GPIb/GPVI, PAR1/4, ADP/P2Y12, TXA2, GPIIb/IIIa, PS exposure |
| ⑦ | Natural Anticoagulants | AT-III/Heparan sulfate, APC-Protein S, TFPI-FXa |
| ⑧ | Fibrinolysis | Plasminogen, tPA/uPA, Plasmin, PAI-1, α2-antiplasmin, TAFI, D-dimer |
| ⑨ | Drug PK | Rivaroxaban (2-cmt), Apixaban (1-cmt), Dabigatran, Warfarin, Enoxaparin, Alteplase |
| ⑩ | Drug PD | FXa inhibition (DOAC), DTI (Dabigatran), VKA (Warfarin VK cycle), Heparin-AT |
| ⑪ | Lab Monitoring | PT/INR, aPTT, Anti-Xa, TT/ECT, D-dimer, Thrombin generation |
| ⑫ | Clinical Endpoints | Proximal/Distal DVT, Low-risk PE, Submassive PE, Massive PE, CTEPH, PTS |

---

## mrgsolve Model: 19 ODE Compartments

### PK Compartments (7)
| Compartment | Drug | Key Parameters |
|-------------|------|---------------|
| `RIV_GUT`, `RIV_CENT`, `RIV_PERIPH` | Rivaroxaban | Ka=1.2/h, CL=4.8L/h, V1=33L, t½=5-9h |
| `APIX_CENT` | Apixaban | Ka=0.78/h, CL=3.3L/h, V1=23L, t½=12h |
| `DABI_CENT` | Dabigatran | Ka=0.35/h, F=6.5%, CL=8.5L/h (renal), t½=12-17h |
| `WARF_CENT` | Warfarin | Ka=0.9/h, CL=0.2L/h (CYP2C9), t½=36h |
| `ENOX_CENT` | Enoxaparin | Ka=0.23/h, CL=0.82L/h, F=92% SC |

### PD Compartments (12)
| Compartment | Represents |
|-------------|-----------|
| `FXa_ACT` | Active Factor Xa (inhibited by DOAC/LMWH-AT) |
| `FIIa_ACT` | Active Thrombin (inhibited by dabigatran/FXa cascade) |
| `FIBRIN_FORM` | Fibrin formation/lysis balance |
| `CLOT_SIZE` | Thrombus burden (0-100%) |
| `PLASMIN_ACT` | Fibrinolytic plasmin activity |
| `DDIMER_CONC` | D-dimer (ng/mL) — fibrinolysis marker |
| `VK_OX`, `VK_RED` | Vitamin K cycle (Warfarin indirect PD) |
| `FVII_POOL`, `FX_POOL`, `FII_POOL` | VK-dependent factor pools (%, t½ 5-60h) |

### PK-PD Linkage
- **FXa Inhibitors**: Hill Emax model; Rivaroxaban EC50=12 ng/mL; Apixaban EC50=5 ng/mL
- **Dabigatran**: Direct thrombin inhibition; EC50=35 ng/mL
- **Warfarin**: Indirect response via VK reduction → ↓ factor synthesis (γ-carboxylation)
- **LMWH**: ATIII activation → FXa inhibition; EC50=0.35 IU/mL Anti-Xa

### Renal Function
- Dabigatran CL adjusted: `CL_eff = CL_DABI × (eGFR/90)^0.85`
- Enoxaparin CL adjusted: `CL_eff = CL_ENOX × (eGFR/90)^0.65`

---

## Treatment Scenarios (6)

| # | Scenario | Regimen | Clinical Setting |
|---|----------|---------|-----------------|
| 1 | **DVT: Rivaroxaban** | 15 mg BID × 21d → 20 mg QD | Acute DVT (confirmed) |
| 2 | **PE: Apixaban** | 10 mg BID × 7d → 5 mg BID | Acute PE (stable) |
| 3 | **Warfarin + Bridge** | Warfarin 5 mg QD + Enoxaparin 1 mg/kg BID × 10d | Classic bridging therapy |
| 4 | **Prophylaxis: Enoxaparin** | 40 mg QD SC × 14d | Post-surgical VTE prevention |
| 5 | **Extended: Rivaroxaban** | 10 mg QD × 180d | Secondary VTE prevention |
| 6 | **CKD: Dabigatran** | 110 mg BID (GFR 30 vs. 90) | Renal impairment sensitivity |

---

## Shiny App: 6 Tabs

| Tab | Content |
|-----|---------|
| **① Patient Profile** | Demographics, Virchow's Triad risk factors, Wells score (DVT/PE), treatment duration recommendation, drug selection based on contraindications |
| **② Drug PK** | Rivaroxaban/Apixaban/Enoxaparin/Dabigatran/Warfarin plasma profiles; Cmax/trough info boxes; renal impairment comparison |
| **③ Coagulation PD** | FXa/FIIa inhibition, thrombin generation (ETP), Warfarin factor depletion/INR, Emax concentration-effect curves |
| **④ Thrombus Dynamics** | Clot resolution over time, fibrin-plasmin balance, D-dimer normalization; info boxes for Day 7/30/90 clot % |
| **⑤ Scenario Comparison** | Side-by-side 6-scenario comparison on any endpoint; efficacy/safety summary table |
| **⑥ Biomarker Dashboard** | INR, Anti-Xa, D-dimer, coagulation cascade activity; biomarker interpretation guide |

---

## Running the Model

```r
# 1. Install dependencies
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "shiny",
                   "shinydashboard", "plotly", "DT", "purrr"))

# 2. Run mrgsolve simulation (6 scenarios)
source("vte_mrgsolve_model.R")

# 3. Launch Shiny dashboard
shiny::runApp("vte_shiny_app.R")
```

```bash
# 4. Render mechanistic map
dot -Tsvg vte_qsp_model.dot -o vte_qsp_model.svg
dot -Tpng -Gdpi=150 vte_qsp_model.dot -o vte_qsp_model.png
```

---

## Key Clinical Parameters Calibrated to Trials

| Parameter | Value | Clinical Source |
|-----------|-------|-----------------|
| Rivaroxaban Cmax (15 mg fed) | ~270 ng/mL | Mueck et al., 2011 |
| Apixaban EC50 (FXa inhib.) | 5 ng/mL | Frost et al., 2013 |
| Dabigatran t½ (GFR 90) | 12-17h | Stangier et al., 2008 |
| Warfarin FVII t½ | ~5.7h (earliest INR rise) | Holford, 1986 |
| Enoxaparin Tmax Anti-Xa | 3-5h (SC) | Hulot et al., 2004 |
| DVT 3-month recurrence (unprovoked) | ~10% | Kyrle et al., 2004 |
| PE CTEPH rate | 2-4% | Pengo et al., 2004 |

---

## License

CC BY 4.0 — QSP Disease Model Library | [pipetcpt/qsp](https://github.com/pipetcpt/qsp)

*Model generated 2026-06-23 by Claude Code Routine (CCR)*
