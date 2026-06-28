# Calcific Aortic Valve Stenosis (CAVD / AS) — QSP Model

> **Disease category:** Chronic Cardiovascular Disease  
> **Abbreviation:** AS / CAVD  
> **Directory:** `aortic-stenosis/`  
> **Model date:** 2026-06-28

---

## Overview

Calcific Aortic Valve Disease (CAVD) leading to **Aortic Stenosis (AS)** is the most common valvular heart disease in adults over 65, affecting ~2.5% of the population and rising sharply with age. The disease progresses through a stereotyped sequence:

1. **Valve endothelial injury** → lipid infiltration → inflammation initiation
2. **Valve Interstitial Cell (VIC) osteogenic transformation** via BMP/Wnt/RANKL signaling
3. **Microcalcification → macrocalcification** → leaflet stiffening
4. **AVA reduction** → transvalvular gradient ↑ → LV pressure overload
5. **Concentric LV hypertrophy** → diastolic dysfunction → eventual systolic failure
6. **Symptom onset** (angina, syncope, dyspnea) → clinical decompensation → death

Without aortic valve replacement (TAVR/SAVR), median survival after symptom onset is 2–3 years. Current pharmacotherapy has not demonstrated efficacy in halting calcification progression (SEAS/ASTRONOMER trials for statins), but novel targets (RANKL/OPG, Lp(a), MGP) are under active investigation.

---

## Mechanistic Map

| File | Description |
|------|-------------|
| [`as_qsp_model.dot`](as_qsp_model.dot) | Full Graphviz mechanistic map (~110+ nodes, 12 subgraphs) |
| [`as_qsp_model.svg`](as_qsp_model.svg) | Vector SVG (scalable) |
| [`as_qsp_model.png`](as_qsp_model.png) | Raster PNG (150 dpi) |

### Map Structure (12 Subgraphs)

| # | Subgraph | Key Nodes |
|---|---------|-----------|
| 1 | Risk Factors & Initiating Stimuli | Age, LPA, LDL, BAV, NOTCH1, HTN, DM, CKD |
| 2 | Valve Endothelium & Initiating Injury | VEC injury, shear stress, oxLDL, eNOS, ROS |
| 3 | Inflammation & Immune Pathways | Macrophage, IL-1β, TNF-α, IL-6, MMP, complement |
| 4 | VIC Osteogenic Transformation | Runx2, Osterix, ALP, RANKL, OPG, Wnt/LRP5, MGP |
| 5 | Valve Calcification Progression | Ca score, AVA, leaflet mobility, macrocalcification |
| 6 | Hemodynamic Consequences | Mean gradient, Vmax, cardiac output, afterload |
| 7 | LV Remodeling | LVH, fibrosis, LVEF, diastolic dysfunction, BNP |
| 8 | Biomarkers & Clinical Endpoints | NYHA, 6MWD, MACE, survival, BNP, Troponin |
| 9 | Drug PK — Statin & PCSK9i | 2-compartment statin, 1-comp SC evolocumab |
| 10 | Drug PD Effects | LDL↓, Lp(a)↓, RANKL↓, MGP activation, anti-fibrosis |
| 11 | Interventional (TAVR/SAVR) | AVA restoration, reverse remodeling, complications |
| 12 | Systemic Consequences | HF, AF, pulmonary HTN, cardiorenal syndrome |

---

## mrgsolve ODE Model

**File:** [`as_mrgsolve_model.R`](as_mrgsolve_model.R)

### Compartments (20 state variables)

| Category | Compartments |
|----------|-------------|
| Statin PK | `STATIN_GUT`, `STATIN_CENTRAL`, `STATIN_PERIPH` |
| PCSK9i PK | `PCSK9I_DEPOT`, `PCSK9I_CENTRAL` |
| Denosumab PK | `DENO_DEPOT`, `DENO_CENTRAL` |
| Vitamin K2 PK | `VK2_DEPOT`, `VK2_CENTRAL` |
| ACEi PK | `ACEI_DEPOT`, `ACEI_CENTRAL` |
| Valve calcification | `CS` (calcium score), `RANKL`, `MGP_carbox` |
| Lipid biomarkers | `LDL_C`, `LPA` |
| Inflammation | `IL6` |
| Neurohormonal | `AngII` |
| LV remodeling | `LVMI`, `COLLAGEN`, `LVEF` |
| Clinical | `NTproBNP` |

### Key Equations

**Valve Calcium Score:**
```
dCS/dt = k_calc × calc_driver × (1 - calc_inhibit) × CS × (1 - CS/CS_max)
calc_driver = 0.3×(LDL/LDL₀) + 0.4×(RANKL/RANKL₀) + 0.3×(Lp(a)/LPA₀)
```

**AVA from Calcium Score (sigmoidal):**
```
AVA = AVA_min + (AVA₀ - AVA_min) / (1 + (CS/CS₅₀)^hill)
```

**Mean Transvalvular Gradient (Bernoulli):**
```
Vmax = CO / (AVA × 60 × 0.785)
ΔP_mean = 2.4 × Vmax²
```

**LV Fibrosis:**
```
dCOLLAGEN/dt = k_form × (AngII/AngII₀) × wall_stress × (1 - E_ACEi) - k_deg × COLLAGEN
```

**LVEF Decline:**
```
dLVEF/dt = -k_loss × (Collagen/Collagen_max) × max(afterload - 130, 0)/50 × LVEF
```

### Treatment Scenarios Simulated

| # | Scenario | Key Drugs |
|---|---------|-----------|
| 1 | No Treatment | — |
| 2 | Statin | Rosuvastatin 20mg QD |
| 3 | Statin + PCSK9i | + Evolocumab 140mg Q2W |
| 4 | Statin + Denosumab* | + Denosumab 60mg Q6M |
| 5 | Statin + Vitamin K2 | + MK-7 180μg QD |
| 6 | Max Medical Therapy | Statin + PCSK9i + VK2 + ACEi |

\* Denosumab for CAVD is investigational (hypothesis-generating)

### Key Simulation Findings

- **Statin alone**: Minimal effect on calcification progression (consistent with negative SEAS/ASTRONOMER/SALTIRE trials)
- **PCSK9i addition**: Meaningful LDL-C and modest Lp(a) reduction; clinically meaningful delay in calcification
- **Denosumab**: Largest hypothetical calcification reduction via RANKL blockade in VIC; awaiting clinical trial validation
- **Vitamin K2 (MK-7)**: Moderate benefit via restoration of carboxylated MGP, natural calcification inhibitor
- **Max medical therapy**: ~15-20% delay in severe AS onset estimated; AVR still required

---

## Shiny App

**File:** [`as_shiny_app.R`](as_shiny_app.R)

### 6 Interactive Tabs

| Tab | Content |
|-----|---------|
| **Patient Profile** | AS severity gauge (AVA, gradient, LVEF, NYHA), risk factors, prognosis |
| **Drug PK** | Time-concentration plots for statin, PCSK9i, denosumab; PK summary table |
| **PD Key Metrics** | LDL-C response, Lp(a), RANKL, IL-6, MGP carboxylation trajectories |
| **Clinical Endpoints** | AVA, mean gradient, LVEF, LVMI, NT-proBNP over time; time-to-severe prediction |
| **Scenario Comparison** | 6-scenario side-by-side comparison with summary table at year 10 |
| **Biomarker Panel** | Calcium score trajectory, collagen fraction, biomarker heatmap, AngII-LV axis |

### Interactive Controls

- Patient sliders: Age, LDL-C, Lp(a), LVEF, baseline AVA, baseline calcium score, simulation duration
- Drug toggles: Statin (dose), PCSK9i (dose, frequency), Denosumab (dose, frequency), Vitamin K2 (dose), ACEi/ARB (dose)
- Auto-running simulation with "Run Simulation" button

---

## Clinical Context

### AS Severity Classification (AHA/ACC 2021)

| Grade | AVA (cm²) | Mean PG (mmHg) | Vmax (m/s) |
|-------|-----------|----------------|------------|
| Mild  | >1.5      | <25            | <3.0       |
| Moderate | 1.0–1.5 | 25–40         | 3.0–4.0    |
| Severe | <1.0     | >40            | >4.0       |
| Very Severe | <0.6 | >60          | >5.0       |

### Typical Natural History (Untreated)

- Asymptomatic severe AS: ~50% develop symptoms within 5 years
- After symptom onset (triad: angina, syncope, HF):
  - Mean survival with angina: 5 years
  - Mean survival with syncope: 3 years
  - Mean survival with HF symptoms: 1-2 years

### Key Molecular Targets

| Target | Rationale | Drug | Trial Status |
|--------|-----------|------|-------------|
| HMG-CoA | LDL ↓ | Statin | Tested/Failed (SEAS) |
| PCSK9 | LDL+Lp(a) ↓ | Evolocumab | Phase III pending |
| RANK/RANKL | VIC osteogenesis ↓ | Denosumab | Phase II (SALTIRE-II) |
| MGP | Calcification inhibitor | Vitamin K2 | Observational |
| Lp(a) | oxPL-mediated calcification | Pelacarsen, Olpasiran | Phase III |
| ACE/AT1R | LV fibrosis ↓ | Ramipril | Supportive |
| IL-6 | Inflammation ↓ | Ziltivekimab | Hypothesis |

---

## References

See [`as_references.md`](as_references.md) for 50 curated PubMed-linked references organized by:
- Epidemiology & Natural History
- Pathophysiology
- RANKL/OPG/VIC transformation
- BMP/Notch/Wnt signaling
- Lp(a) and lipid-mediated calcification
- Matrix Gla Protein & Vitamin K
- Statin trials (SEAS, ASTRONOMER, SALTIRE)
- Novel targets (PCSK9i, Denosumab)
- TAVR/SAVR evidence (PARTNER, Evolut)
- LV remodeling & neurohormonal
- Hemodynamics & echocardiography
- CT calcium scoring
- QSP modeling

---

## Files Summary

| File | Description | Size |
|------|-------------|------|
| `as_qsp_model.dot` | Graphviz mechanistic map (12 subgraphs, 110+ nodes) | ~35KB |
| `as_qsp_model.svg` | Rendered vector image | ~150KB |
| `as_qsp_model.png` | Rendered raster image (150 DPI) | ~2.5MB |
| `as_mrgsolve_model.R` | ODE model (20 compartments, 6 scenarios) | ~15KB |
| `as_shiny_app.R` | Interactive Shiny dashboard (6 tabs) | ~20KB |
| `as_references.md` | 50 curated PubMed references | ~8KB |
| `README.md` | This file | ~8KB |

---

## Acknowledgments

Model parameters calibrated against:
- SEAS trial data (Rossebø et al., NEJM 2008)
- PARTNER / Evolut trial outcomes
- CT calcium score progression studies (Marechaux et al., Heart 2010)
- Echocardiographic valve area natural history (Pellikka et al., Circulation 2005)
- ESC/AHA AS management guidelines (2017/2021)

*This model is intended for research/educational purposes only and does not constitute clinical guidance.*
