# Myocarditis QSP Model

[![Model](myo_qsp_model.png)](myo_qsp_model.svg)

## Overview

**Myocarditis** (심근염) is an inflammatory disease of the myocardium caused primarily by viral infection (coxsackievirus B3, SARS-CoV-2, adenovirus, HHV-6), autoimmune activation, or toxic exposure. It is a leading cause of sudden cardiac death in young adults and a major pathway to dilated cardiomyopathy (DCM).

This Quantitative Systems Pharmacology (QSP) model comprehensively captures the three-phase pathophysiology of myocarditis:
1. **Viral phase**: Cardiomyocyte infection via CAR/ACE2 receptors, viral replication, and innate immune sensing (TLR3/7/9, RIG-I, MDA5, cGAS-STING)
2. **Inflammatory phase**: NK cell and macrophage-mediated injury, cytokine storm (TNF-α, IL-6, IL-1β, IFN-γ), CD4+/CD8+ T cell activation, molecular mimicry, and anti-cardiac antibody production
3. **Remodeling phase**: TGF-β–driven fibroblast activation, collagen deposition, ventricular dilation, and chronic heart failure progression

---

## Files

| File | Description |
|------|-------------|
| `myo_qsp_model.dot` | Graphviz mechanistic map (170+ nodes, 10 clusters) |
| `myo_qsp_model.svg` | Vector mechanistic map (interactive) |
| `myo_qsp_model.png` | Raster image (150 dpi) |
| `myo_mrgsolve_model.R` | mrgsolve ODE model (35 compartments, 5 treatment scenarios) |
| `myo_shiny_app.R` | Interactive Shiny dashboard (7 tabs) |
| `myo_references.md` | 60 PubMed-cited references (14 sections) |

---

## Mechanistic Map Structure

The `.dot` mechanistic map contains **10 subgraph clusters** and **170+ nodes**:

| Cluster | Nodes | Key Components |
|---------|-------|----------------|
| Viral Entry & Replication | 16 | CVB3, SARS-CoV-2, CAR receptor, ACE2, viral proteases 2A/3C |
| Pattern Recognition & Innate Signaling | 22 | TLR3/7/9, RIG-I, MDA5, cGAS-STING, NF-κB, IRF3/7, NLRP3 |
| Innate Immune Effectors & Cytokines | 24 | NK, M1/M2 macrophages, mDC, pDC, IFN-α/β/γ, TNF-α, IL-1β, IL-6, IL-12 |
| Adaptive Immune Response | 23 | CD4+ Th1/Th2/Th17/Treg/Tfh, CD8+ CTL, B cells, plasma cells, germinal center |
| Autoimmunity & Molecular Mimicry | 14 | Anti-myosin, anti-β1AR, anti-ANT, anti-TnI antibodies; ADCC, complement fixation |
| Cardiomyocyte Injury & Death | 18 | Healthy/infected/injured CMC, apoptosis, pyroptosis, troponin/BNP leak |
| Cardiac Remodeling & Fibrosis | 17 | TGF-β → myofibroblast → collagen, RAAS, ECM remodeling, LV dilation |
| Clinical Endpoints & Biomarkers | 16 | Troponin, BNP, LVEF, CMR-LGE, ECG, echocardiography, arrhythmia risk |
| Pharmacological Interventions | 19 | IVIG, prednisone, azathioprine, cyclosporine, colchicine, rituximab, JAK inhibitors |
| Giant Cell Myocarditis | 10 | Multinucleated giant cells, T cell pathology, AV block, combined IS therapy |

---

## mrgsolve ODE Model

### Compartments (35 ODEs)

| Category | Compartments |
|----------|-------------|
| Cardiomyocytes (3) | Healthy (H), Infected (I), Dead (D) |
| Viral load (1) | V (copies/mL) |
| Innate cells (3) | NK cells, M1 macrophages, M2 macrophages |
| Cytokines (7) | IFN-β, IFN-γ, TNF-α, IL-6, IL-1β, TGF-β, IL-10 |
| Adaptive immunity (9) | Naive CD4, Th1, Th17, Treg, Naive CD8, CTL, Naive B, Plasma cells, Antibodies |
| Cardiac remodeling (3) | Cardiac fibroblasts, Myofibroblasts, Collagen |
| Biomarkers (3) | Troponin I, BNP, LVEF |
| Drug PK (6) | IVIG, Prednisone, Azathioprine→6-MP, Cyclosporine, Colchicine |

### Key Mechanistic Equations

```
dV/dt = p_V·I − c_V·V·(1 + kIFN·IFNβ/(IFNβ50+IFNβ)) − NK_kill·NK·V
dH/dt = rH·H·(1−(H+I)/Hmax) − d_H·H − βV·V·H − CTL_bys·CTL·H·(1−E_IS)
dTroponin/dt = kTnLeak·(δI·I + kNec·IS_signal·H) − dTn·Troponin
dEF/dt = krec·(EF_target − EF)
```

### Drug PK Parameters

| Drug | Route | Bioavailability | Vd (L/kg) | t½ |
|------|-------|----------------|------------|-----|
| IVIG | IV infusion | 100% | 3.5 | 21 days |
| Prednisone | Oral | 80% | 0.97 | 2–3 h |
| Azathioprine → 6-MP | Oral | 47% | 0.8 | ~2 h |
| Cyclosporine | Oral | 35% | 4.0 | ~8–12 h |
| Colchicine | Oral | 45% | 250 | ~26–31 h |

### 5 Treatment Scenarios

| # | Scenario | Rationale |
|---|----------|-----------|
| 1 | No Treatment (Natural History) | Baseline disease course for comparison |
| 2 | IVIG monotherapy (2 g/kg IV) | Anti-antibody, Fc-R blockade; used in fulminant myocarditis |
| 3 | Prednisone + Azathioprine (TIMIC protocol) | Frustaci 2009 TIMIC trial; virus-negative inflammatory CM |
| 4 | Triple IS: IVIG + Pred + Aza + Cyclosporine | Giant cell myocarditis protocol (Cooper 2007) |
| 5 | IVIG + Colchicine | Myopericarditis; NLRP3/inflammasome suppression |

---

## Shiny Dashboard (7 Tabs)

| Tab | Content |
|-----|---------|
| **Overview** | Disease summary, value boxes (peak Tn, BNP, nadir EF, outcome), quick plots |
| **PK Profiles** | Drug dosing setup, concentration-time plots for all 5 drugs |
| **Viral & Innate** | Viral load, cardiomyocyte dynamics, innate cells, cytokine kinetics |
| **PD Biomarkers** | Troponin I, BNP kinetics, drug effect profiles, anti-cardiac antibodies |
| **Clinical Endpoints** | LVEF trajectory, outcome prediction, adaptive immune response, biomarker table |
| **Scenario Comparison** | 3-arm comparison of any combination of treatment scenarios |
| **Fibrosis & Remodeling** | Myofibroblast/collagen dynamics, fibrosis cytokines, EF–fibrosis correlation |

---

## Key Clinical Thresholds

| Biomarker | Threshold | Interpretation |
|-----------|-----------|----------------|
| Troponin I | > 0.04 ng/mL | Myocardial injury |
| Troponin I | > 10 ng/mL | Severe myocarditis |
| BNP | > 100 pg/mL | Heart failure |
| BNP | > 400 pg/mL | Decompensated HF |
| LVEF | ≥ 50% | Normal |
| LVEF | 35–50% | Mildly reduced |
| LVEF | < 35% | DCM threshold |
| CMR-LGE | Present | Fibrosis marker; predicts arrhythmia |

---

## Disease Epidemiology

- **Incidence**: 1–10/100,000 person-years (likely underdiagnosed)
- **Age**: Peaks in young adults (20–40 years), second peak > 60 years
- **Sex**: Male predominance (2:1 for viral myocarditis)
- **Etiology**: Viral (CVB3, SARS-CoV-2, adenovirus, parvovirus B19, HHV-6) most common; giant cell myocarditis rare (<1%)
- **Prognosis**:
  - Complete recovery: 40–50% of acute myocarditis cases
  - Progression to DCM: 20–30%
  - Giant cell myocarditis: 5-year survival < 50% without heart transplant

---

## Running the Model

### Prerequisites

```r
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "shiny",
                   "shinydashboard", "plotly", "DT"))
```

### mrgsolve ODE Simulation

```r
source("myo_mrgsolve_model.R")
# Compiles and runs all 5 scenarios automatically
# Generates plots for troponin, BNP, LVEF, immune cells, fibrosis
```

### Shiny Interactive Dashboard

```r
shiny::runApp("myo_shiny_app.R")
```

### Mechanistic Map Rendering

```bash
fdp -Tsvg myo_qsp_model.dot -o myo_qsp_model.svg
fdp -Tpng -Gdpi=150 myo_qsp_model.dot -o myo_qsp_model.png
```

---

## Key References

- Caforio ALP et al. ESC Position Statement on Myocarditis. *Eur Heart J.* 2013;34:2636–2648.
- Frustaci A et al. TIMIC Study (Pred+Aza). *Eur Heart J.* 2009;30:1995–2002.
- McNamara DM et al. IMAC-2 (IVIG RCT). *Circulation.* 2001;103:2254–2259.
- Cooper LT et al. Giant Cell Myocarditis. *J Card Fail.* 2007;13:620–625.
- Bozkurt B et al. Myocarditis Review. *J Am Coll Cardiol.* 2021;78:2208–2245.
- Fung G et al. Myocarditis Mechanisms. *Circ Res.* 2016;118:496–514.

> **Full reference list**: See [`myo_references.md`](myo_references.md) for 60 PubMed-cited references across 14 sections.
