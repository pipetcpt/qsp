# Long COVID (PASC) — QSP Model

> **Post-Acute Sequelae of SARS-CoV-2 Infection (PASC)** — Quantitative Systems Pharmacology Model  
> Generated: 2026-06-26 | Abbreviation: `pasc`

---

## Disease Overview

**Long COVID**, formally termed Post-Acute Sequelae of SARS-CoV-2 (PASC), is a multi-system condition affecting **10–30% of COVID-19 survivors**. Symptoms persist beyond 12 weeks post-infection and span fatigue, cognitive impairment ("brain fog"), autonomic dysfunction (POTS), dyspnea, and musculoskeletal pain. With hundreds of millions of COVID-19 infections globally, PASC represents one of the largest emerging chronic disease burdens of the 21st century.

### Key Pathological Mechanisms Modeled

| Domain | Key Mechanism | Model Components |
|--------|--------------|-----------------|
| **Viral Persistence** | SARS-CoV-2 RNA/antigen persists in GI, lymph nodes, CNS | V_PLASMA, V_RES, V_AG |
| **Immune Dysregulation** | T/B cell exhaustion, IFN-I dysregulation, autoantibodies | IFN, CD8_exh, Auto_Ab |
| **Cytokine Storm Residue** | Persistent IL-6, TNF-α, IL-1β elevation | IL6, TNF |
| **Endothelial/Vascular** | Fibrin microthrombi, platelet hyperactivation, D-dimer↑ | Fibrin, Ddimer |
| **Neuroinflammation** | BBB disruption, microglial activation, serotonin depletion | BBB, Microglia, Serotonin |
| **Autonomic Dysfunction** | POTS, dysautonomia, anti-adrenergic autoantibodies | AutNom |
| **Mitochondrial Dysfunction** | ROS↑, ATP deficit, anaerobic shift, PEM | ROS, MitoDmg, Lactate |
| **Gut Dysbiosis** | Leaky gut, LPS translocation, SCFA depletion | (upstream inputs) |

---

## Model Architecture

### Mechanistic Map

[![PASC QSP Model](pasc_qsp_model.png)](pasc_qsp_model.svg)

- **Nodes**: 100+ components across 10 subgraph clusters
- **Clusters**: Viral Persistence · Immune Dysregulation · Endothelial/Vascular · Neuroinflammation · Autonomic · Mitochondrial/Energy · Gut-Immune Axis · Hormonal · Drug PK/PD · Clinical Endpoints

### mrgsolve ODE Model (`pasc_mrgsolve_model.R`)

- **28 compartments**: viral kinetics (4) + immune (5) + vascular (2) + neurological (4) + autonomic (1) + mitochondrial (3) + PK (8)
- **5+ treatment scenarios** with clinical trial benchmarks:

| Scenario | Treatment | Key Mechanism | Clinical Benchmark |
|----------|-----------|--------------|-------------------|
| S1 | No treatment | Natural history | — |
| S2 | Nirmatrelvir 300mg BID ×15d | Viral replication inhibition | RECOVER-VITAL (NCT05595369) |
| S3 | Metformin 500mg BID | AMPK↑, IL-6↓, mito protection | COVID-OUT: 41% LCC reduction (Bramante 2023 *Lancet Infect Dis*) |
| S4 | Low-dose naltrexone 4.5mg QD | Microglial suppression, TLR4↓ | LDN PASC pilot trials |
| S5 | Sertraline 50mg QD | Serotonin↑, platelet aggregation↓ | σ1R anti-inflammatory pathway |
| S6 | Nirmatrelvir + Metformin | Viral+metabolic dual target | Combination strategy |
| S7 | Full combination | All mechanisms | Theoretical optimal |

### Shiny Dashboard (`pasc_shiny_app.R`)

8 interactive tabs:

| Tab | Content |
|-----|---------|
| 1. Patient Profile | Disease domain radar chart, value boxes (FSS, VO2max, MoCA, POTS) |
| 2. Pharmacokinetics | Multi-drug PK curves, PK parameter table, MoA summaries |
| 3. Viral & Immune | Viral kinetics, immune biomarkers, cytokines, coagulation |
| 4. Neuro & Autonomic | Neuroinflammation cascade, POTS trajectory, mitochondrial, serotonin |
| 5. Clinical Endpoints | FSS, VO2max, MoCA, SF-36 PCS, week 52 summary table |
| 6. Scenario Comparison | All 7 scenarios overlaid, bar chart comparison, endpoint table |
| 7. Virtual Population | n=20–500 patient simulation, response rate analysis |
| 8. Biomarker Panel | NfL, CRP, D-dimer, biomarker reference ranges |

---

## Key Parameters & Calibration

| Parameter | Value | Source |
|-----------|-------|--------|
| Viral reservoir reactivation (kActiv) | 0.02/day | Proal & VanElzakker 2021 |
| IL-6 production rate (kIL6) | 0.50/day | Phetsouphanh et al. 2022 (Nat Immunol) |
| BBB disruption rate (kBBB) | 0.20/day | Fernández-Castañeda et al. 2022 (Cell) |
| Nirmatrelvir IC50 | 0.003 µg/mL | Hammond et al. 2022 (NEJM) |
| Metformin IL-6 EC50 | 500 ng/mL | Saisho 2015 |
| POTS orthostatic HR threshold | ≥30 bpm | Dani et al. 2021 |
| LDN EC50 (microglial suppression) | 2.0 ng/mL | Younger & Mackey 2009 |
| Sertraline SERT EC50 | ~50 ng/mL | Standard PK/PD |

---

## Clinical Endpoints

| Endpoint | Baseline (PASC) | Normal Range | Treatment Target |
|----------|----------------|--------------|-----------------|
| FSS (fatigue) | 5.8–6.5 | 1–3.5 | ≤4 (responder) |
| VO2max (% predicted) | 40–60% | 85–100% | ≥70% |
| MoCA (cognitive) | 20–25 | ≥26 | ≥26 |
| POTS ΔHR (bpm) | 35–50 | <30 | <30 |
| SF-36 PCS | 25–40 | 50+ | ≥45 |
| NfL (pg/mL) | 15–50 | <10 | <15 |

---

## Files

| File | Description |
|------|-------------|
| [`pasc_qsp_model.dot`](pasc_qsp_model.dot) | Graphviz mechanistic map (100+ nodes, 10 clusters) |
| [`pasc_qsp_model.svg`](pasc_qsp_model.svg) | Vector format (high resolution) |
| [`pasc_qsp_model.png`](pasc_qsp_model.png) | Raster format (150 dpi) |
| [`pasc_mrgsolve_model.R`](pasc_mrgsolve_model.R) | mrgsolve ODE model (28 compartments, 7 scenarios, VP analysis) |
| [`pasc_shiny_app.R`](pasc_shiny_app.R) | Shiny interactive dashboard (8 tabs) |
| [`pasc_references.md`](pasc_references.md) | 62 references across 15 sections |

---

## Usage

```bash
# Render mechanistic map
dot -Tsvg pasc_qsp_model.dot -o pasc_qsp_model.svg
dot -Tpng -Gdpi=150 pasc_qsp_model.dot -o pasc_qsp_model.png
```

```r
# Run mrgsolve model
install.packages(c("mrgsolve", "dplyr", "ggplot2", "purrr"))
source("pasc_mrgsolve_model.R")

# Launch Shiny dashboard
install.packages(c("shiny", "shinydashboard", "plotly", "DT"))
shiny::runApp("pasc_shiny_app.R")
```

---

## Key References

- **Davis et al. 2023** — *Nat Rev Microbiol*: Major mechanisms of long COVID
- **Klein et al. 2023** — *Nature*: Immune profiling distinguishing features of PASC
- **Bramante et al. 2023** — *Lancet Infect Dis*: Metformin 41% reduction in long COVID (COVID-OUT RCT)
- **Pretorius et al. 2021** — *Cardiovasc Diabetol*: Fibrin microclots in PASC
- **Fernández-Castañeda et al. 2022** — *Cell*: Mild COVID causes myelin and neural dysregulation
- **Wong et al. 2023** — *Cell*: Serotonin depletion in post-viral infection
- **Su et al. 2022** — *Cell*: Early predictors of PASC (autoantibodies, EBV reactivation, cort cortisol)
