# Gastroparesis (위마비) — QSP Model

[![Mechanistic Map](gp_qsp_model.png)](gp_qsp_model.svg)

---

## Disease Overview

**Gastroparesis** (위마비, delayed gastric emptying) is a chronic gastric motility disorder defined by objectively delayed gastric emptying in the absence of mechanical obstruction. The hallmark is reduced gastric emptying rate resulting in symptoms of nausea, vomiting, early satiety, bloating, and postprandial fullness.

**Prevalence**: ~5 million US adults; diabetic (~29% of T1DM, ~1% of T2DM), idiopathic (~36%), post-surgical (~13%)

**Diagnosis**: Gastric emptying scintigraphy (GES) — >10% retention at 4 hours is the gold standard

---

## Pathophysiology Summary

```
Hyperglycemia / DM
    ↓
Oxidative stress (ROS ↑, nNOS ↓)
    ↓
ICC loss (c-Kit+ pacemaker cells ↓) + nNOS neuron depletion
    ↓
Impaired gastric slow waves → ↓ antral contractions
    ↓
Pyloric hypertonicity (unchecked NO-mediated relaxation)
    ↓
Delayed gastric emptying (GER ↓, 4h retention ↑)
    ↓
GCSI symptoms: nausea, vomiting, early satiety, bloating
```

**Three Core Defects**:
1. **ICC Loss** — Interstitial cells of Cajal (pacemaker cells) depleted by 40–60% in DM-gastroparesis; drives slow-wave dysrhythmia
2. **nNOS Neuron Depletion** — Neuronal nitric oxide synthase-expressing neurons absent in ~85% of gastroparesis biopsies; loss of pyloric relaxation
3. **Vagal Neuropathy** — Impairs central regulation of antral contractility and gastric accommodation

---

## Drug Targets & Treatment Classes

| Drug | Class | Main Target | Prokinetic Mechanism |
|------|-------|-------------|----------------------|
| **Metoclopramide** | D2/5HT4 | D2 receptor (GI+CNS) | D2 block → ↑ACh; partial 5HT4 |
| **Domperidone** | D2 antagonist | Peripheral D2 | No BBB penetration (safer EPS) |
| **Erythromycin** | Motilin agonist | Motilin receptor | Phase III MMC induction |
| **Prucalopride** | 5-HT4 agonist | 5-HT4 receptor | Full agonist → ↑ACh, antral drive |
| **Relamorelin** | Ghrelin analogue | GHSR | Antral contractility ↑ |
| **Ondansetron** | 5-HT3 antagonist | CTZ 5-HT3 | Antiemetic (not prokinetic) |

---

## Model Components

### 1. Mechanistic Map (`gp_qsp_model.dot/.svg/.png`)

- **100+ nodes** across 10 pathophysiological clusters:
  - CNS & Vomiting Control (Area Postrema, NTS, DMV, CTZ)
  - Autonomic Nervous System (vagal efferents/afferents, sympathetic)
  - Enteric Nervous System & ICC (myenteric plexus, c-Kit/SCF, slow waves)
  - Neurotransmitter Signaling (D2, 5HT3, 5HT4, motilin, ghrelin, ACh)
  - Gastric Motor Function (antrum, pylorus, fundus, GER)
  - Core Pathophysiology (oxidative stress, nNOS, iNOS, macrophage M1/M2)
  - Diabetic Mechanisms (AGE-RAGE, microvascular, autonomic neuropathy)
  - Drug PK (7 drugs with distinct PK parameter sets)
  - Drug PD (receptor occupancy and downstream antral/pyloric effects)
  - Clinical Outcomes (GCSI, GES T½, 4h retention, nausea VAS)

### 2. mrgsolve ODE Model (`gp_mrgsolve_model.R`)

**16 ODE Compartments**:

| # | Compartment | Description |
|---|-------------|-------------|
| 1 | `MCP_GI` | Metoclopramide GI absorption depot |
| 2 | `MCP_plasma` | Metoclopramide central plasma (ng/mL) |
| 3 | `MCP_CNS` | Metoclopramide CNS/effect-site |
| 4 | `DOM_plasma` | Domperidone plasma |
| 5 | `ERY_plasma` | Erythromycin plasma |
| 6 | `PRU_plasma` | Prucalopride plasma |
| 7 | `REL_plasma` | Relamorelin plasma |
| 8 | `D2_effect` | D2 receptor occupancy (indirect response) |
| 9 | `HT4_effect` | 5-HT4 receptor agonism |
| 10 | `nNOS_act` | nNOS enzyme activity (normalised 0–1) |
| 11 | `ICC_dens` | ICC density (normalised 0–1) |
| 12 | `Antral_c` | Antral contractility (normalised 0–1) |
| 13 | `Pyloric_t` | Pyloric tone (0=open, 1=tight) |
| 14 | `GasVol` | Gastric content volume (mL) |
| 15 | `GER_cum` | Cumulative gastric emptying (mL) |
| 16 | `GCSI_dyn` | GCSI composite score (0–5) |

**7 Treatment Scenarios**:
- **S0**: Untreated gastroparesis (baseline)
- **S1**: Metoclopramide 10 mg QID
- **S2**: Domperidone 10 mg TID
- **S3**: Erythromycin 250 mg TID (short-term motilin agonism)
- **S4**: Prucalopride 2 mg QD
- **S5**: Relamorelin 100 mcg SC BID (investigational)
- **S6**: Prucalopride 2 mg QD + Ondansetron 8 mg TID

**Clinical Calibration**:
- APPROVE trial (Camilleri 2014): Metoclopramide GCSI ↓0.5 units vs PBO
- PRED trial (McCallum 2021): Prucalopride GES T½ improves ~20 min
- Relamorelin Phase 2b (Camilleri 2017, JAMA IM): Vomiting frequency ↓83% vs PBO

### 3. Shiny App (`gp_shiny_app.R`)

**7 Interactive Tabs**:

| Tab | Content |
|-----|---------|
| **Patient Profile** | Disease type, ICC/nNOS severity sliders, baseline disease state |
| **Pharmacokinetics** | Concentration-time profiles for all 6 drugs |
| **Pharmacodynamics** | D2/5HT4/motilin/GHSR receptor occupancy/activation |
| **GI Motility** | Gastric volume, GER, antral contractility, pyloric tone |
| **Clinical Endpoints** | GCSI score, 4h gastric retention, nausea VAS, value boxes |
| **Scenario Comparison** | All 7 scenarios side-by-side bar chart and table |
| **Biomarkers** | ICC density, nNOS activity, laboratory interpretation |

### 4. References (`gp_references.md`)

**53 PubMed-indexed references** across:
- Epidemiology & clinical features (5 refs)
- ICC loss & ENS degeneration (6 refs)
- nNOS & oxidative stress (4 refs)
- Autonomic neuropathy & diabetes (4 refs)
- Diagnosis & clinical assessment (4 refs)
- Treatment — prokinetics (8 refs)
- Emerging therapies & gastric electrical stimulation (5 refs)
- Motilin & ghrelin signaling (3 refs)
- QSP modeling approaches (4 refs)
- Biomarkers & pathology (5 refs)
- Clinical trials & outcomes (5 refs)

---

## Key Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| ICC density (DM-GP) | 35–45% of normal | Grover 2011 |
| nNOS neurons absent | ~85% of GP biopsies | Grover 2011 |
| Metoclopramide D2 IC50 | 5 ng/mL | Tonini 2004 |
| Prucalopride 5-HT4 EC50 | 3 ng/mL | Tack 2012 |
| Relamorelin GHSR EC50 | 15 ng/mL | Camilleri 2017 |
| Normal GER | ~0.025/min (~1.5%/min) | Siegel 1988 |
| Normal 4h retention | <10% | Abell 2008 |

---

## Files

| File | Description |
|------|-------------|
| `gp_qsp_model.dot` | Graphviz mechanistic map source |
| `gp_qsp_model.svg` | Vector graphic (interactive) |
| `gp_qsp_model.png` | Raster image (150 dpi) |
| `gp_mrgsolve_model.R` | ODE model with 16 compartments, 7 scenarios |
| `gp_shiny_app.R` | Interactive 7-tab Shiny dashboard |
| `gp_references.md` | 53 annotated PubMed references |

---

## Run Instructions

```r
# mrgsolve model
library(mrgsolve)
source("gp_mrgsolve_model.R")  # runs all 7 scenarios & generates PDF plots

# Shiny app
library(shiny)
shiny::runApp("gp_shiny_app.R")
```

---

*Model created: 2026-06-26 | QSP Disease Model Library*
