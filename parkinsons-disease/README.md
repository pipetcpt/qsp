# Parkinson's Disease — QSP Model

> **Disease:** Parkinson's Disease (PD)  
> **Category:** Neurodegenerative Disease  
> **Model version:** 1.0 · Generated 2026-06-20

[![Mechanistic Map](pd_qsp_model.png)](pd_qsp_model.svg)

---

## Overview

Parkinson's Disease is the second most common neurodegenerative disorder, affecting ~10 million people worldwide. It is characterized by the **progressive loss of dopaminergic neurons in the substantia nigra pars compacta (SNpc)** and the accumulation of **Lewy bodies** (α-synuclein aggregates) in surviving neurons.

This QSP model integrates:
- α-Synuclein aggregation kinetics (monomer → oligomer → fibril → Lewy body)
- Mitochondrial dysfunction and oxidative stress cascades
- Neuroinflammation (microglia/astrocyte activation, cytokine networks)
- Dopamine synthesis, release, and metabolism
- Basal ganglia motor circuit (direct/indirect pathway imbalance)
- Population PK/PD for 4 drug classes

---

## Mechanistic Map

| File | Description |
|------|-------------|
| [`pd_qsp_model.dot`](pd_qsp_model.dot) | GraphViz source (16 subgraphs, 140+ nodes) |
| [`pd_qsp_model.svg`](pd_qsp_model.svg) | Scalable vector graphic (interactive) |
| [`pd_qsp_model.png`](pd_qsp_model.png) | Raster image (150 dpi) |

### Subgraph Clusters

| # | Cluster | Key Nodes |
|---|---------|-----------|
| 1 | Genetic Risk Factors | SNCA, LRRK2, PINK1, PRKN, DJ-1, GBA |
| 2 | α-Synuclein Aggregation | Monomer, Oligomers (toxic), Fibrils, Lewy Bodies |
| 3 | Mitochondrial Dysfunction | Complex I, ROS, mPTP, PINK1/Parkin mitophagy |
| 4 | Protein Quality Control | UPS, 26S Proteasome, CMA, Lysosome, mTOR/AMPK |
| 5 | Neuroinflammation | Microglia M1/M2, Astrocytes, TNF-α, IL-1β, NLRP3 |
| 6 | DA Neuron Pathophysiology | SNpc neuron pool, TH, AADC, MAO-B, COMT, DAT |
| 7 | Dopamine Metabolism | DA synthesis, D1/D2/D3/D4 receptors, DOPAC, HVA |
| 8 | Basal Ganglia Circuit | Striatum D1/D2 MSNs, GPi, GPe, STN, Thalamus |
| 9 | Clinical Manifestations | Tremor, Bradykinesia, Rigidity, LID, UPDRS, H&Y |
| 10 | Levodopa/Carbidopa PK | Absorption, BBB transport, wearing-off, ON-OFF |
| 11 | Dopamine Agonists | Pramipexole, Ropinirole, Rotigotine, Apomorphine |
| 12 | MAO-B/COMT Inhibitors | Rasagiline, Selegiline, Entacapone, Opicapone |
| 13 | Neuroprotection | GDNF, LRRK2 inhibitors, α-Syn immunotherapy, DBS |
| 14 | Gut-Brain Axis | Gut microbiome, Vagus nerve, ENS, Pesticides |
| 15 | Biomarkers | DaTscan, 18F-DOPA PET, CSF α-Syn, NfL, MRI SN |
| 16 | Population PK | Ka, Vd, CL, BBB transport, food-drug interaction |

---

## mrgsolve ODE Model

**File:** [`pd_mrgsolve_model.R`](pd_mrgsolve_model.R)

### Compartments (22 ODEs)

| # | Compartment | Description |
|---|-------------|-------------|
| 1 | `ASyn_M` | α-Syn monomer (nM) |
| 2 | `ASyn_O` | α-Syn oligomers — toxic species (nM) |
| 3 | `ASyn_F` | α-Syn fibrils/Lewy inclusions (nM) |
| 4 | `ROS` | Reactive oxygen species index (a.u.) |
| 5 | `NEUROINF` | Neuroinflammation index (a.u.) |
| 6 | `SNpc` | SNpc DA neuron pool (fraction 0–1) |
| 7 | `DA_syn` | Synaptic dopamine (nM) |
| 8 | `DA_brain` | Brain dopamine pool (nM) |
| 9 | `LD_gut` | Levodopa gut compartment (mg) |
| 10 | `LD_C` | Levodopa central plasma (mg/L) |
| 11 | `LD_P` | Levodopa peripheral compartment (mg/L) |
| 12 | `LD_brain` | Levodopa brain compartment (mg/L) |
| 13 | `PRAM_gut` | Pramipexole gut (mg) |
| 14 | `PRAM_C` | Pramipexole plasma (µg/L) |
| 15 | `PRAM_brain` | Pramipexole brain (µg/L) |
| 16 | `RAS_gut` | Rasagiline gut (mg) |
| 17 | `RAS_C` | Rasagiline plasma (µg/L) |
| 18 | `MAOB_act` | Active MAO-B enzyme fraction (0–1) |
| 19 | `ENT_gut` | Entacapone gut (mg) |
| 20 | `ENT_C` | Entacapone plasma (µg/L) |
| 21 | `UPDRS_III` | UPDRS-III motor score (continuous) |
| 22 | `LID_risk` | Cumulative LID risk index (a.u.) |

### Treatment Scenarios

| # | Scenario | Regimen | Key Finding |
|---|----------|---------|-------------|
| 1 | Untreated | None | 60% SNpc loss by year 10; UPDRS ↑ unchecked |
| 2 | Levodopa TID | 250 mg LD × 3/day from diagnosis | Best acute motor control; highest LID risk |
| 3 | Pramipexole TID | 0.75 mg × 3/day | Continuous D2/D3 stim → lower LID risk |
| 4 | Rasagiline QD | 1 mg QD | MAO-B 98% inhibition; possible neuroprotection (ADAGIO) |
| 5 | LD + Entacapone | 250 mg LD + 200 mg ENT × 3/day | Extended L-DOPA t½; reduced wearing-off |
| 6 | Triple Therapy | LD + Pramipexole + Rasagiline | Best long-term motor control (complex regimen) |
| 7 | Continuous Delivery | LD CR (q4h) + Pramipexole | Reduced pulsatile stim → lowest LID risk |

### Key Parameters Calibrated to Clinical Data

| Parameter | Value | Reference |
|-----------|-------|-----------|
| `kSN_death` | 0.00019 day⁻¹ | ~15 years from preclinical to 60% loss |
| `kTH` | 0.8 nmol mg⁻¹ h⁻¹ | Hornykiewicz 1998 |
| `Ka_LD` | 1.2 h⁻¹ | Homma et al. 2020 |
| `CL_LD` | 1.4 L/h/kg | Population PK studies |
| `kBBB_LD` | 0.3 h⁻¹ | LAT1 transport kinetics |
| `Imax_MAOB` | 1.0 | Rasagiline irreversible inhibition |
| `IC50_MAOB` | 0.0003 µM | Covalent MAO-B inhibitor |

---

## Shiny Dashboard

**File:** [`pd_shiny_app.R`](pd_shiny_app.R)

### Tab Structure

| Tab | Content |
|-----|---------|
| **Patient Profile** | Demographics, Hoehn & Yahr stage, genetic risk, Braak staging map |
| **Pharmacokinetics** | L-DOPA PK profiles, food-drug interaction, formulation comparison |
| **PD Biomarkers** | α-Syn oligomers, neuroinflammation, SNpc neuron survival, synaptic DA |
| **Motor Endpoints** | UPDRS-III trajectory, LID risk, GPi output, D2R stimulation |
| **Scenario Comparison** | Side-by-side multi-drug regimen comparison with outcome table |
| **Neuroprotection** | Disease-modifying interventions: GDNF, LRRK2 inhibitors, α-Syn immunotherapy, DBS |

**Launch:**
```r
shiny::runApp("pd_shiny_app.R")
```

---

## References

**File:** [`pd_references.md`](pd_references.md)

50 references organized into 14 sections:
1. Epidemiology · 2. α-Synuclein Aggregation · 3. Genetics  
4. Mitochondria & Oxidative Stress · 5. Neuroinflammation  
6. Dopamine Neurotransmission & Basal Ganglia · 7. Levodopa/Carbidopa  
8. Dopamine Agonists · 9. MAO-B & COMT Inhibitors  
10. QSP Modeling · 11. Neuroprotection · 12. Biomarkers & Clinical Trials  
13. Gut-Brain Axis · 14. Deep Brain Stimulation

---

## Model Limitations

- Simplified Euler integration used in Shiny app (use mrgsolve version for research)
- Spatial heterogeneity within SNpc not modeled
- Individual variability (CYP2D6 polymorphisms) requires PopPK extension
- Lewy body spreading to cortex modeled as scalar, not anatomical connectivity
- Non-motor symptoms (RBD, anosmia, depression) represented as surrogates

---

## Quick Start

```r
# 1. Install dependencies
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "shiny",
                   "shinydashboard", "plotly", "DT"))

# 2. Run mrgsolve model
source("pd_mrgsolve_model.R")

# 3. Launch Shiny app
shiny::runApp("pd_shiny_app.R")
```

---

*Generated by Claude Code Routine (CCR) | QSP Disease Model Library | 2026-06-20*
