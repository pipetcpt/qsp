# Pancreatic Ductal Adenocarcinoma (PDAC) QSP Model

## Disease Overview

Pancreatic ductal adenocarcinoma (PDAC) is the most lethal solid malignancy, with a 5-year survival rate of ~12% (all stages combined). It accounts for ~3% of all cancer diagnoses but ~7% of cancer deaths in the United States. Approximately 80% of patients present with locally advanced or metastatic disease, making curative resection impossible.

**Key molecular characteristics:**
- **KRAS mutations**: ~95% (G12D ~44%, G12V ~26%, G12R ~14%, G12C ~2%)
- **TP53 loss/mutation**: ~75%
- **CDKN2A (p16) loss**: ~90%
- **SMAD4 loss**: ~55%
- **HRD (BRCA1/2 germline + somatic, PALB2, ATM)**: ~15-20%
- **MSI-H**: ~1%
- Dense desmoplastic stroma: 80–90% of tumor volume
- Profound immunosuppression (PD-L1 upregulation, MDSC, Treg)

## Mechanistic Map Preview

[![PDAC QSP Mechanistic Map](pdac_qsp_model.png)](pdac_qsp_model.svg)

*Click image to view full interactive SVG. Map contains 179 nodes across 10 biological clusters.*

## Deliverables

| File | Description | Format |
|------|-------------|--------|
| `pdac_qsp_model.dot` | Graphviz mechanistic map (179 nodes, 10 clusters) | DOT |
| `pdac_qsp_model.svg` | Vector mechanistic map | SVG |
| `pdac_qsp_model.png` | Raster mechanistic map (150 dpi) | PNG |
| `pdac_mrgsolve_model.R` | ODE PK/PD model (26 compartments, 7 scenarios) | R/mrgsolve |
| `pdac_shiny_app.R` | Interactive 6-tab dashboard | R/Shiny |
| `pdac_references.md` | Annotated bibliography (35+ references) | Markdown |

## Key Biological Pathways Modeled

| Pathway | Key Nodes | Drug Targets |
|---------|-----------|--------------|
| KRAS/MAPK | KRAS→BRAF→MEK1/2→ERK1/2→c-MYC | MRTX1133 (G12D), Adagrasib (G12C), SOS1-i |
| PI3K/AKT/mTOR | PI3K→PIP3→AKT(pT308)→mTORC1→S6K1/4EBP1 | Copanlisib, Everolimus |
| TGF-β/SMAD/EMT | TGFβ→SMAD2/3→SMAD4→SNAIL/ZEB→E-Cad loss | Galunisertib (TGFβR1-i) |
| Desmoplastic Stroma | PSC→myoCAF/iCAF→Collagen/HA→IFP↑→Drug penetration↓ | PEGPH20 (hyaluronidase), Pirfenidone |
| Immune TME | PD-L1↑, Treg↑, PMN-MDSC↑, TAM-M2↑, CD8↓ | Pembrolizumab, Ipilimumab, Relatlimab |
| Angiogenesis | HIF-1α→VEGF-A→VEGFR-2→Neovascularization | Bevacizumab, Ramucirumab |
| DNA Damage/HRD | BRCA1/2 loss→HRD→PARP dependency→Synthetic lethality | Olaparib, Rucaparib |
| Cell Cycle | CDK4/6-CyclinD1→pRb→E2F; CDKN2A loss→CDK4↑ | Abemaciclib (CDK4/6-i) |
| p53 Pathway | MDM2→p53 degradation; TP53 mut (GoF) | APR-246 (p53 reactivation) |
| Wnt/Notch/Hh | Wnt/β-Catenin, Hedgehog, NOTCH → stemness | Vismodegib (Hh-i) |

## Treatment Scenarios & Clinical Trial Reference Data

| # | Regimen | Doses | mPFS (mo) | mOS (mo) | ORR (%) | Key Trial |
|---|---------|-------|-----------|----------|---------|-----------|
| 1 | Untreated Control | — | 1.5 | 3.0 | 0 | — |
| 2 | Gemcitabine mono | 1000 mg/m² IV days 1,8,15 q28d | 3.7 | 6.7 | 7 | MPACT control arm |
| 3 | Gem + nab-Paclitaxel | Gem 1000 + nab-Pac 125 mg/m² days 1,8,15 q28d | 5.5 | 8.5 | 23 | MPACT (Von Hoff 2013, PMID: 24131140) |
| 4 | FOLFIRINOX | OHP 85 + CPT11 180 + 5-FU 400 bolus + 5-FU 2400 CI q14d | 6.4 | 11.1 | 31.6 | PRODIGE4 (Conroy 2011, PMID: 21561347) |
| 5 | mFOLFIRINOX | OHP 65 + CPT11 150 + 5-FU 2400 CI q14d (no bolus) | 6.0 | 10.5 | 28 | Clinical adaptation |
| 6 | MRTX1133 (KRAS G12D) | 100 mg BID oral | ~4.0 | ~8.0 | ~40 | Phase I/II ongoing |
| 7 | Olaparib (BRCA+) | 300 mg BID oral | 7.4 | NS | — | POLO (Golan 2019, PMID: 30945893) |

*mPFS = median PFS; mOS = median OS; ORR = objective response rate; NS = not significant; CI = continuous infusion*

## Model Architecture

### PK Components (mrgsolve ODE Model)

**Drug PK (12 compartments):**
- **Gemcitabine**: 2-compartment IV model (GEM_C1, GEM_C2) → intracellular dFdCTP (active)
- **nab-Paclitaxel**: 2-compartment IV model (NPAC_C1, NPAC_C2)
- **Oxaliplatin**: 1-compartment free plasma (OHP_FREE) → DNA adduct compartment (OHP_DNA)
- **Irinotecan/SN-38**: CPT11 plasma → SN38 (active) → SN38G (glucuronide)
- **5-FU**: 1-compartment plasma (FU5) → FdUMP (active metabolite)
- **MRTX1133**: oral 1-compartment (KRAS G12D inhibitor, F=45%)
- **Olaparib**: oral 1-compartment (PARP inhibitor, F=73%)

**Disease PD (14 compartments):**
- **Simeoni TGI model**: Proliferating cells (x0) → Transit compartments x1→x2→x3 → Tumor volume
- **KRAS signaling**: Normalized KRAS activity (0–1), inhibited by MRTX1133
- **Stroma resistance**: Dynamic stromal factor modulating effective drug penetration
- **CA19-9 biomarker**: Proportional to total tumor burden (x0 + TUMOR)
- **Friberg myelosuppression**: Prol → Tr1 → Tr2 → Tr3 → Circ (ANC, driven by Gem + SN-38)

### PD Model Features
- Bliss independence for drug combination effects
- Stroma-mediated drug penetration reduction (up to 60% at maximum stroma)
- KRAS G12D-specific enhancement by MRTX1133 on tumor growth rate
- Grade 3/4 neutropenia prediction via Friberg model nadir

## How to Run

### Prerequisites

```r
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "shiny", "bslib", "plotly"))
```

### Run mrgsolve Simulations

```r
source("pdac_mrgsolve_model.R")
# Outputs: tumor_volume_plot, ca199_plot, anc_plot, kras_signaling_plot
```

### Launch Shiny Dashboard

```r
shiny::runApp("pdac_shiny_app.R")
```

### Render Mechanistic Map

```bash
# From within pancreatic-cancer/ directory:
dot -Tsvg pdac_qsp_model.dot -o pdac_qsp_model.svg
dot -Tpng -Gdpi=150 pdac_qsp_model.dot -o pdac_qsp_model.png
```

## Shiny Dashboard Tabs

| Tab | Content |
|-----|---------|
| **Patient Profile** | Disease stage, molecular subtype (Classical/Basal-like/QM/Exocrine), KRAS mutation, HRD/MSI status, ECOG, BSA, risk radar chart |
| **Drug PK** | Plasma concentration-time profiles for all 7 drugs; log/linear scale; PK parameter table |
| **Tumor Dynamics** | Tumor volume over time (Simeoni TGI), CA19-9 trajectory, waterfall plot (best % change) |
| **Biomarkers** | CA19-9 kinetics, ctDNA VAF, ANC myelosuppression; toxicity risk table |
| **Clinical Endpoints** | OS/PFS parametric curves, ORR bar chart, G3/4 toxicity summary table |
| **Scenario Comparison** | Forest plot (HR for OS), comprehensive comparison table, sensitivity tornado plot |

## Clinical Context & Unmet Need

PDAC is primarily diagnosed at advanced/metastatic stage (~80% at diagnosis), with median OS of 3–4 months without treatment. Key challenges and emerging solutions:

1. **Late diagnosis**: No validated screening biomarker; CA19-9 is Lewis antigen–dependent (absent in 5–10% of patients)
2. **Chemoresistance**: Dense stroma (IFP > 40 mmHg) impedes drug delivery; KRAS drives intrinsic resistance via PI3K/AKT/mTOR survival signals
3. **Immune desert phenotype**: Low mutational burden (TMB ~1 mut/Mb); MDSC/Treg infiltration; CXCL1/5-driven exclusion of CD8+ T cells
4. **Metabolic reprogramming**: KRAS-driven macropinocytosis recycles amino acids; autophagy sustains survival under nutrient stress
5. **Emerging targeted therapies**:
   - **KRAS G12D inhibitors**: MRTX1133 (Phase I/II, Mirati), RMC-9805 (Revolution Medicines)
   - **SOS1 inhibitors**: BI-3406 (synergistic with MEK-i)
   - **PARP inhibitors + DDR**: Olaparib maintenance in BRCA+ (POLO trial)
   - **mRNA vaccines**: Personalized neoantigen vaccines (Moderna mRNA-4157 + pembrolizumab)
   - **Bispecific antibodies**: Anti-EGFR × anti-CD3 (KN026 + KN046)
   - **CAR-T cell therapy**: Anti-mesothelin, anti-MSLN

## References

See [`pdac_references.md`](pdac_references.md) for complete bibliography with 35+ annotated references organized by topic.

---

*Model created: 2026-06-23 | QSP Disease Model Library*
