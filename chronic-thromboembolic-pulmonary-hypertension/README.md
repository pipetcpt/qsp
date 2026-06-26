# CTEPH QSP Model

## Chronic Thromboembolic Pulmonary Hypertension (CTEPH)

[![Category](https://img.shields.io/badge/Category-Cardiovascular%2FPulmonary-blue)]()
[![Compartments](https://img.shields.io/badge/ODE%20Compartments-19-green)]()
[![Scenarios](https://img.shields.io/badge/Scenarios-7-orange)]()
[![References](https://img.shields.io/badge/References-60-red)]()

---

## Overview

**Chronic Thromboembolic Pulmonary Hypertension (CTEPH)** is a life-threatening condition arising from incomplete resolution of pulmonary emboli, leading to permanent mechanical obstruction of the pulmonary vasculature by organized fibrous thrombi. Superimposed on this fixed obstruction, progressive small-vessel disease driven by endothelin-1 overproduction, nitric oxide deficiency, and abnormal shear stress creates an additional variable component of elevated pulmonary vascular resistance (PVR). The combination causes progressive right ventricular pressure overload, remodeling, and ultimately failure.

CTEPH is classified as Group 4 pulmonary hypertension (ESC/ERS 2022). Diagnosis requires:
- mPAP > 20 mmHg at rest
- Pulmonary capillary wedge pressure ≤ 15 mmHg (pre-capillary)
- Mismatched perfusion defects on V/Q scan after ≥ 3 months of therapeutic anticoagulation
- Confirmation by CT pulmonary angiography or conventional angiography

---

## Disease Mechanism

```
Acute PE → Incomplete fibrinolysis → Organized thrombus → Fixed PVR obstruction
                 ↓                          ↓
    Thrombophilic state            Secondary vascular remodeling
         + recurrent PE              (ET-1, NO↓, PGI₂↓, SMC proliferation)
                                          ↓
                                   Variable PVR↑
                                          ↓
                                   Total PVR↑ → mPAP↑
                                          ↓
                               RV pressure overload
                                          ↓
                    Adaptive hypertrophy → Dilatation → Failure
                                          ↓
                               CO↓, SaO₂↓, BNP↑, 6MWD↓
```

### Key Pathological Processes
| Component | Mechanism | Clinical Consequence |
|-----------|-----------|---------------------|
| Fixed PVR | Organized thrombus, neointima, fibrosis | Mechanical obstruction, V/Q mismatch |
| Variable PVR | ET-1↑, NO↓, PGI₂↓, SMC proliferation | Vasomotor constriction, remodeling |
| RV remodeling | Pressure overload, ischemia, fibrosis | TAPSE↓, RV EF↓, RV dilation |
| Gas exchange | V/Q mismatch, dead space↑ | Hypoxemia, exertional dyspnea |
| Neurohumoral | BNP↑, RAAS activation, sympathetic | Sodium retention, arrhythmia |

---

## Model Files

| File | Description |
|------|-------------|
| `cteph_qsp_model.dot` | Graphviz mechanistic map source (12 clusters, 134+ nodes) |
| `cteph_qsp_model.svg` | High-resolution vector mechanistic map |
| `cteph_qsp_model.png` | 150 dpi raster mechanistic map |
| `cteph_mrgsolve_model.R` | mrgsolve ODE model (19 compartments, 7 scenarios) |
| `cteph_shiny_app.R` | Interactive Shiny dashboard (7 tabs) |
| `cteph_references.md` | 60 PubMed-linked references |

---

## Mechanistic Map

[![CTEPH QSP Mechanistic Map](cteph_qsp_model.png)](cteph_qsp_model.svg)

**12 Subgraph Clusters:**
1. **Thrombotic Origin & Acute PE** — Thrombophilia, platelet activation, recurrent PE
2. **Fibrinolysis Dysfunction** — tPA deficiency, PAI-1↑, plasmin↓, fibrin cross-linking
3. **Pulmonary Vascular Remodeling** — ET-1, NO↓, PDGF, TGF-β, SMC proliferation, plexiform lesions
4. **Hemodynamics** — mPAP, PVR (fixed + variable), CO, TPG, PCWP
5. **Right Ventricular Remodeling** — Hypertrophy, dilatation, fibrosis, TAPSE, RVEF
6. **Gas Exchange** — V/Q mismatch, dead space, PaO₂, SaO₂, DLCO
7. **Neurohumoral Biomarkers** — BNP, RAAS, ET-1, troponin, D-dimer
8. **Riociguat PK/PD** — sGC stimulation, cGMP↑, vasodilation, anti-proliferative
9. **ERA PK/PD** — Macitentan/ACT-132577, ETA/ETB blockade, anti-fibrotic
10. **Prostacyclin PK/PD** — Treprostinil, IP receptor, cAMP↑, platelet inhibition
11. **Surgical/Procedural** — PEA (≥70% PVR reduction), BPA (staged), anticoagulation
12. **Clinical Endpoints** — 6MWD, WHO FC, dyspnea, BNP, survival, TTCW

---

## mrgsolve ODE Model

### Compartments (19 total)

**Drug PK (7):**
| # | State Variable | Description |
|---|----------------|-------------|
| 1 | `C1_RIO` | Riociguat central (ng/mL) |
| 2 | `C2_RIO` | Riociguat peripheral |
| 3 | `MET_RIO` | M1 active metabolite |
| 4 | `C1_MAC` | Macitentan central (ng/mL) |
| 5 | `C2_MAC` | Macitentan peripheral |
| 6 | `MET_MAC` | ACT-132577 active metabolite |
| 7 | `C1_TREP` | Treprostinil central (ng/mL) |

**Disease PD (12):**
| # | State Variable | Description |
|---|----------------|-------------|
| 8 | `TB` | Thrombotic burden (0–1, normalized) |
| 9 | `PVR_fixed` | Fixed PVR — thrombus (dyn·s/cm⁵) |
| 10 | `PVR_var` | Variable PVR — vasomotor/remodeling |
| 11 | `ET1` | Plasma endothelin-1 (pg/mL) |
| 12 | `cGMP` | cGMP second messenger (pmol/mL) |
| 13 | `cAMP` | cAMP second messenger (pmol/mL) |
| 14 | `RV_work` | RV stroke work index (g·m/m²) |
| 15 | `mPAP` | Mean pulmonary artery pressure (mmHg) |
| 16 | `CO` | Cardiac output (L/min) |
| 17 | `SaO2` | Arterial O₂ saturation (%) |
| 18 | `BNP` | BNP plasma level (pg/mL) |
| 19 | `sixMWD` | 6-minute walk distance (m) |

### Treatment Scenarios

| Scenario | Intervention | Expected 6MWD Change | PVR Change |
|----------|-------------|----------------------|------------|
| 1 | No treatment | Decline ~−20 m/yr | Progressive ↑ |
| 2 | Anticoagulation only | Stable | Stable (no regression) |
| 3 | Riociguat 2.5 mg TID | +46 m (wk16, CHEST-1) | −31% |
| 4 | Macitentan 10 mg QD | +35–40 m | −25% |
| 5 | Riociguat + Macitentan | +60–75 m | −45% |
| 6 | BPA (5 sessions) + Riociguat | +70–90 m | −50% |
| 7 | Post-PEA + Riociguat + Macitentan | +120–150 m | −75% |

### Key PK Parameters

| Drug | CL (L/h) | Vc (L) | t½ (h) | Tmax (h) | F (%) |
|------|----------|--------|--------|----------|-------|
| Riociguat | 2.4 | 30 | 8.7 | 1.5 | 94 |
| Macitentan | 1.1 | 50 | 31 | 8 | 75 |
| ACT-132577 | 0.55 | 60 | 75 | — | — |
| Treprostinil | 4.0 | 14 | 2.4 | — | 79 (SC) |

---

## Shiny Dashboard (7 Tabs)

| Tab | Content |
|-----|---------|
| **1. Patient Profile** | Baseline hemodynamics, biomarkers, procedural history, disease severity table |
| **2. Pharmacokinetics** | PK concentration-time profiles for all 3 drugs + metabolites, PK parameter tables |
| **3. PD Signals** | cGMP/cAMP dynamics, ET-1, PVR decomposition (fixed vs variable), thrombotic burden |
| **4. Hemodynamics** | mPAP, total PVR, CO, SaO₂, RV stroke work, hemodynamic summary table |
| **5. Clinical Endpoints** | 6MWD, BNP, WHO FC, waterfall change plot, efficacy summary at Wk16/Wk52 |
| **6. Scenario Comparison** | Side-by-side comparison of 7 predefined clinical scenarios |
| **7. Biomarkers** | BNP/ET-1/TB trajectories, ESC/ERS risk stratification, current risk profile |

### Running the App
```r
install.packages(c("shiny","shinydashboard","mrgsolve","dplyr","tidyr",
                   "ggplot2","plotly","DT"))
library(shiny)
shiny::runApp("cteph_shiny_app.R")
```

---

## Clinical Background

### Epidemiology
- Incidence: ~3–5% of acute PE survivors develop CTEPH
- Prevalence: ~3–30 per million (likely underdiagnosed)
- Male:Female ratio ≈ 1:1
- Median survival without treatment: ~2–3 years from diagnosis

### Prognosis by Hemodynamics
| mPAP | 5-Year Survival |
|------|-----------------|
| 25–35 mmHg | ~70–80% |
| 35–55 mmHg | ~40–60% |
| >55 mmHg | <20% |

### Treatment Algorithm
```
Confirmed CTEPH diagnosis
         ↓
Lifelong anticoagulation (essential for all)
         ↓
CTEPH Expert Center multidisciplinary team assessment
         ↓
Is the patient operable? (proximal disease, good cardiopulmonary reserve)
    /           \
  YES            NO
   ↓              ↓
  PEA         Residual/inoperable CTEPH
(curative)         ↓
                Riociguat ± Macitentan
                   ± BPA (if accessible)
                   ± Prostacyclin (severe cases)
```

---

## References
60 references covering:
- ESC/ERS clinical guidelines (2022)
- CHEST-1/CHEST-2 (riociguat), SERAPHIN (macitentan), MERIT-1 (combination)
- PEA outcomes (Madani et al., Jamieson et al.)
- BPA technique and results (Mizoguchi, Wiedenroth, Ogawa)
- QSP/PK-PD modeling methodology

See [`cteph_references.md`](cteph_references.md) for the complete reference list.

---

*Built by Claude Code Routine | QSP Disease Model Library | 2026-06-26*
