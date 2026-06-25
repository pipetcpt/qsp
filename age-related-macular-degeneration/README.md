# Age-related Macular Degeneration (AMD) QSP Model

[![QSP](https://img.shields.io/badge/QSP-mrgsolve%20%C2%B7%20Shiny%20%C2%B7%20Graphviz-success)]()
[![Compartments](https://img.shields.io/badge/ODE%20compartments-20-blue)]()
[![Scenarios](https://img.shields.io/badge/clinical%20scenarios-6-orange)]()
[![References](https://img.shields.io/badge/references-55-green)]()

---

## Overview

Age-related macular degeneration (AMD) is the **leading cause of irreversible central vision loss** in adults over 50 in high-income countries, affecting ~200 million people worldwide. The disease is driven by the convergence of:

- **Drusen formation** and Bruch's membrane thickening (lipid deposition, RPE failure)
- **Complement system dysregulation** (CFH Y402H, ARMS2 A69S genetic variants)
- **VEGF-driven choroidal neovascularization (CNV)** — wet/neovascular AMD
- **Geographic atrophy (GA)** — progressive RPE and photoreceptor death (dry/atrophic AMD)

This QSP model integrates all four major biological axes into a **20-ODE compartmental system** capable of simulating anti-VEGF pharmacokinetics, VEGF/Ang-2 pathway suppression, complement activation, RPE cell dynamics, and clinically observable endpoints (BCVA, OCT-CST, CNV area, GA area).

---

## Mechanistic Map

[![AMD QSP Mechanistic Map](amd_qsp_model.png)](amd_qsp_model.svg)

*Click to open the full-resolution SVG. The map includes 10 subgraph clusters and ~130 nodes.*

### Key Clusters in the Map

| Cluster | Contents |
|---------|----------|
| Drug PK — Systemic | IVT injection → vitreous → retina/RPE → systemic drainage; 5 drugs |
| VEGF / Angiogenesis | HIF-1α/2α → VEGF-A165 → VEGFR-2 → PI3K/Akt/mTOR/ERK → EC proliferation/migration |
| Ang-2 / Tie2 | ANG1/ANG2 competition → Tie2 → vessel stabilization/destabilization (faricimab target) |
| Complement System | C3/C5 convertase (alternative) → MAC (C5b-9); CFH regulation; pegcetacoplan/avacopan |
| RPE / Bruch's | Lipofuscin A2E → oxidative stress → mtDNA damage → NLRP3/IL-1β → RPE apoptosis |
| CNV Formation | Choroidal vessel invasion → Type 1/2/3 CNV → SRF/IRF/PED → CST |
| Neuroinflammation | Microglia, macrophages, IL-6/TNFα/MCP-1 → VEGF amplification |
| Genetics | CFH Y402H, ARMS2 A69S, C3 R102G, HTRA1; smoking, UV, aging |
| Clinical Endpoints | BCVA (ETDRS letters), CST (μm), CNV area (mm²), GA area (mm²) |
| Pharmacological Rx | q4w/q8w/q12w dosing, T&E, PRN; AREDS2, PDT, gene therapy |

---

## mrgsolve Model (`amd_mrgsolve_model.R`)

### Compartments (20 ODEs)

| # | Compartment | Unit | Description |
|---|-------------|------|-------------|
| 1 | DRUG_VIT | nM | Drug in vitreous |
| 2 | DRUG_RET | nM | Drug at retina/RPE (effect site) |
| 3 | DRUG_SYS | mg | Drug in systemic circulation |
| 4 | VEGF_FREE | nM | Free VEGF-A in retina |
| 5 | VEGF_BOUND | nM | Drug:VEGF complex |
| 6 | VEGFR2_ACT | 0–1 | VEGFR-2 activation state |
| 7 | ANG2_FREE | nM | Free Ang-2 |
| 8 | ANG2_BOUND | nM | Drug:Ang-2 (faricimab) |
| 9 | C3_LOCAL | AU | Local retinal C3 complement |
| 10 | C5_LOCAL | AU | Local retinal C5 complement |
| 11 | MAC_LOCAL | AU | MAC (C5b-9) level |
| 12 | RPE_NORM | 0–1 | Normal RPE cell fraction |
| 13 | RPE_DAM | 0–1 | Damaged RPE fraction |
| 14 | LIPOFUSCIN | AU | Lipofuscin / A2E accumulation |
| 15 | DRUSEN | mm² | Drusen area |
| 16 | CNV_AREA | mm² | CNV lesion area |
| 17 | FLUID_EX | μm | Excess fluid above CST baseline |
| 18 | GA_AREA | mm² | Geographic atrophy area |
| 19 | BCVA_SCORE | letters | ETDRS best-corrected VA |
| 20 | PR_FRAC | 0–1 | Photoreceptor survival |

### Drug PK Parameters

| Drug | Dose | MW | Kd (VEGF) | t½ vitreous | Mechanism |
|------|------|----|-----------|-------------|-----------|
| Ranibizumab | 0.5 mg | 48 kDa | 0.04 nM | ~7.2 d | Anti-VEGF-A Fab |
| Aflibercept | 2 mg | 115 kDa | 0.0005 nM | ~7.0 d | VEGF trap (VEGF-A/B, PlGF) |
| Bevacizumab | 1.25 mg | 149 kDa | 0.2 nM | ~9.0 d | Anti-VEGF-A IgG |
| Faricimab | 6 mg | 150 kDa | 0.0003 nM | ~7.0 d | Dual VEGF-A + Ang-2 |
| Brolucizumab | 6 mg | 26 kDa | 0.06 nM | ~4.0 d | Anti-VEGF-A scFv |

### Clinical Scenarios

| # | Scenario | Regimen | Key Parameters |
|---|----------|---------|---------------|
| 1 | Ranibizumab standard | q4w ×3 → q8w ×10 | Kd=0.04 nM, 2-yr simulation |
| 2 | Aflibercept VIEW | q4w ×3 → q8w ×10 | Kd=0.0005 nM, VIEW1/2 calibrated |
| 3 | Faricimab TENAYA | q4w ×4 → q16w T&E | Dual VEGF/Ang-2; TENAYA/LUCERNE |
| 4 | Brolucizumab HAWK | q6w ×3 → q12w ×8 | Small scFv, fast clearance |
| 5 | Natural history | No treatment | Disease progression only |
| 6 | Dry AMD + AREDS2 | Supplements only | GA progression, 4-yr |

---

## Shiny Application (`amd_shiny_app.R`)

### Tabs

| Tab | Content |
|-----|---------|
| **1. Patient Profile** | AMD staging map, baseline parameter table, risk factor summary |
| **2. Drug PK** | Vitreous & retinal PK curves, free vs bound VEGF, PK summary table |
| **3. PD Key Markers** | VEGFR-2 suppression, free VEGF, Ang-2, PD–BCVA correlation plot |
| **4. Clinical Endpoints** | BCVA, CST, CNV area, GA area — treated vs untreated |
| **5. Scenario Comparison** | All 5 drugs side-by-side BCVA/CNV; injection burden vs efficacy plot |
| **6. Biomarker Explorer** | Complement C3/C5/MAC, RPE/lipofuscin, photoreceptor, drusen burden |

### Interactive Controls

- Drug selection (Ranibizumab / Aflibercept / Bevacizumab / Faricimab / Brolucizumab)
- AMD type (Wet/Dry)
- Simulation duration (180–1460 days)
- Loading doses, loading interval, maintenance interval
- Baseline BCVA, initial CNV area, initial excess fluid

### Launch

```r
install.packages(c("shiny","shinydashboard","mrgsolve","dplyr","ggplot2","tidyr","ggrepel"))
shiny::runApp("amd_shiny_app.R")
```

---

## Key Biological Insights Modeled

1. **Anti-VEGF binding affinity drives sustained suppression**: Aflibercept (Kd ~0.5 pM) achieves stronger and more durable VEGF capture than ranibizumab (Kd ~40 pM) per molecule, but the larger dose of brolucizumab/faricimab compensates with high molar drug concentration.

2. **Dual VEGF/Ang-2 blockade (faricimab)**: Ang-2 destabilizes vessels and amplifies VEGF-driven leakage; simultaneous neutralization synergistically reduces fluid and extends treatment intervals to q16w in ~45% of patients (TENAYA/LUCERNE).

3. **Complement amplification of GA**: Local C3/MAC deposition accelerates RPE atrophy, creating a vicious cycle (RPE death → impaired complement regulation → more MAC → more RPE death). Pegcetacoplan (C3 inhibitor) demonstrated ~22% slowing of GA in OAKS/DERBY trials.

4. **Lipofuscin A2E as drusen driver**: A2E accumulation in RPE lysosomes impairs phagocytosis → drusen growth → complement deposition. This feed-forward loop is the central driver of disease progression from intermediate to late AMD.

5. **BCVA as composite endpoint**: The model captures how visual acuity reflects contributions of CNV area (hemorrhage/fluid), GA area (photoreceptor loss), CST elevation (fluid distortion), and photoreceptor survival — providing mechanistic interpretation of clinical trial endpoints.

---

## Clinical Trial Calibration Notes

| Trial | Drug | 1-yr BCVA gain | Model prediction |
|-------|------|----------------|-----------------|
| MARINA/ANCHOR | Ranibizumab 0.5mg q4w | +7.2 letters | ~+8.2 letters |
| VIEW 1/2 | Aflibercept 2mg q8w | +8.4 letters | ~+9.1 letters |
| TENAYA/LUCERNE | Faricimab 6mg q4w→q16w | +10.7 letters | ~+11.0 letters |
| HAWK | Brolucizumab 6mg q12w | +6.6 letters | ~+7.8 letters |
| CATT (no Tx) | Natural history | −14.9 letters | ~−14.2 letters |

---

## File Structure

```
age-related-macular-degeneration/
├── README.md                    ← This file
├── amd_qsp_model.dot            ← Graphviz mechanistic map source
├── amd_qsp_model.svg            ← Vector graphic (scalable)
├── amd_qsp_model.png            ← Raster graphic (150 dpi)
├── amd_mrgsolve_model.R         ← 20-ODE mrgsolve model + 6 scenarios
├── amd_shiny_app.R              ← Interactive 6-tab Shiny dashboard
└── amd_references.md            ← 55 PubMed-linked references
```

---

## References

See [`amd_references.md`](amd_references.md) for 55 curated references covering epidemiology, pathophysiology, PK/PD, clinical trials, and QSP modeling methodology.

**Key papers**: Edwards et al. *Science* 2005 (CFH risk variant); Brown et al. *NEJM* 2006 (ANCHOR trial); CATT Research Group *NEJM* 2011; Khanani et al. *Ophthalmology* 2022 (TENAYA/LUCERNE); Hutton-Smith et al. *IOVS* 2018 (mechanistic AMD model).
