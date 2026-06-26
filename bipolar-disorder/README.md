# Bipolar Disorder – QSP Model

> **Category:** Neuropsychiatric Disorders  
> **Abbreviation:** BD  
> **Date Added:** 2026-06-25  

---

## Overview

Bipolar disorder (BD) is a severe, recurrent psychiatric illness characterised by episodes of mania/hypomania alternating with periods of depression, affecting ~2.4 % of the global population. It is among the top causes of disability in working-age adults, with a lifetime risk of suicide of approximately 15–20 times the general population rate.

The disease involves dysregulation across multiple neurobiological systems simultaneously:

- **Monoamine neurotransmission** (dopamine hyperactivity in mania; serotonin/NE deficits in depression)
- **Signal transduction** (GSK-3β, PKC, cAMP/PKA, PI3K/AKT pathways)
- **Neuroplasticity** (BDNF ↓, reduced hippocampal volume, synaptic remodelling)
- **Neuroinflammation** (elevated IL-6, TNF-α, CRP across both phases)
- **HPA axis dysregulation** (hypercortisolaemia, blunted feedback)
- **Circadian rhythm disruption** (CLOCK/BMAL1, PER/CRY mutations)
- **Ion channel abnormalities** (CACNA1C risk allele → L-type Ca²⁺↑)

---

## Mechanistic Map

[![BD QSP Mechanistic Map](bd_qsp_model.png)](bd_qsp_model.svg)

*Click the image to open the interactive SVG.*

### Map Specifications

| Attribute | Value |
|-----------|-------|
| Total nodes | 120+ |
| Clusters | 12 |
| Drug targets shown | Lithium, Valproate, Quetiapine, Lamotrigine, Aripiprazole |
| PK compartments | 10 (drug-specific 1- and 2-cmt oral models) |
| PD pathways | Neurotransmitter, signal transduction, neuroplasticity, HPA, circadian, neuroinflammation, ion channels, gut-brain axis |

### Cluster Summary

| # | Cluster | Key Nodes |
|---|---------|-----------|
| 1 | Neurotransmitter Systems | DA, 5-HT, NE, GABA, Glu synaptic pools & receptors |
| 2 | Signal Transduction | GSK-3β, PKA/PKC, MAPK/ERK, AKT, mTOR, IP3/DAG |
| 3 | Ion Channels | Nav1.x, Cav1.2 (CACNA1C), Kv, HCN |
| 4 | HPA & Circadian | CRH→ACTH→Cortisol, CLOCK/BMAL1, PER/CRY, SCN |
| 5 | Neuroplasticity | BDNF/TrkB, neurogenesis, dendritic spine density, hippocampal volume |
| 6 | Neuroinflammation | Microglia, NF-κB, IL-6, TNF-α, NLRP3, IDO/kynurenine, ROS |
| 7 | Drug MOA | Inhibitor/agonist/antagonist annotations per drug |
| 8 | Drug PK | Gut → central (→ peripheral for Li) → concentration compartments |
| 9 | Clinical Endpoints | YMRS, MADRS, HAM-D, CGI-BP, GAF, response/remission thresholds |
| 10 | Genetic/Epigenetic | CACNA1C, ANK3, CLOCK, BDNF Val66Met, SLC6A4, COMT, miRNA-134/132 |
| 11 | Brain Circuits | PFC, amygdala, hippocampus, striatum/NAc, VTA, raphe, LC, hypothalamus |
| 12 | Gut-Brain & Metabolic | Gut microbiome, LPS, SCFA, insulin resistance, weight gain |

---

## mrgsolve ODE Model

**File:** `bd_mrgsolve_model.R`

### Compartments (22 ODEs)

| # | Compartment | Units |
|---|-------------|-------|
| 1–3 | Lithium gut, central, peripheral | mmol |
| 4–5 | Valproate gut, central | mg |
| 6–8 | Quetiapine gut, central + norquetiapine central | mg |
| 9–10 | Lamotrigine gut, central | mg |
| 11 | Dopamine neurotransmission index | dimensionless |
| 12 | Serotonin neurotransmission index | dimensionless |
| 13 | GSK-3β activity index | dimensionless |
| 14 | BDNF level index | dimensionless |
| 15 | IL-6 / neuroinflammation index | dimensionless |
| 16 | Cortisol / HPA index | dimensionless |
| 17 | YMRS score | points (0–60) |
| 18 | MADRS score | points (0–60) |
| 19 | Body weight change | kg |
| 20–21 | Circadian oscillator state & derivative | — |
| 22 | GAF score | points (0–100) |

### Key Parameters (selected)

| Parameter | Value | Source |
|-----------|-------|--------|
| Li CL (renal) | 1.80 L/h | Finley et al. 1995 |
| Li Vc | 30 L (0.4 L/kg) | Sproule 2002 |
| QTP F (bioavailability) | 9 % | DeVane & Nemeroff 2001 |
| VPA fu₀ (low Cc) | 10 % | Perucca 2002 |
| GSK-3β IC₅₀ Li | 0.70 mEq/L | Ryves & Harwood 2001 |
| BDNF EC₅₀ Li | 0.40 mEq/L | Frey et al. 2007 |
| YMRS Emax Li | 18 points | Bowden et al. 1994 / Cipriani 2013 |
| MADRS Emax QTP | 15 points | Young et al. 2010 EMBOLDEN |

### Treatment Scenarios

| # | Scenario | Primary Endpoint |
|---|----------|-----------------|
| 1 | Lithium 900 mg/d monotherapy (21 d) | Acute mania – YMRS |
| 2 | Valproate 1000 mg/d monotherapy (21 d) | Acute mania – YMRS |
| 3 | Quetiapine 300 mg QD (56 d) | BD depression – MADRS |
| 4 | Lithium + Quetiapine 300 mg QD (56 d) | BD depression – MADRS remission |
| 5 | Lithium maintenance (1 year) | BDNF, GSK-3β, long-term mood stability |
| 6 | Lamotrigine titration 25→50→100→200 mg/d (112 d) | BD-II depression – MADRS |

---

## Shiny Dashboard

**File:** `bd_shiny_app.R`

Run with:

```r
library(shiny)
shiny::runApp("bd_shiny_app.R")
```

### Tabs

| Tab | Content |
|-----|---------|
| 1. Patient Profile | Demographics, BD subtype, episode history, pharmacogenomics (CACNA1C / BDNF Val66Met / COMT / CYP3A4), CANMAT recommendations |
| 2. Pharmacokinetics | PK curves for all drugs, Css summary table, lithium therapeutic window |
| 3. PD Biomarkers | GSK-3β, BDNF, IL-6, cortisol dynamics; neuroplasticity & neuroinflammation panels |
| 4. Clinical Endpoints | YMRS + MADRS trajectories, GAF score, response/remission KPI value boxes |
| 5. Scenario Comparison | Side-by-side comparison of up to 8 predefined CANMAT-aligned scenarios |
| 6. Safety Monitor | Lithium toxicity alert, VPA monitoring, estimated QTc, weight change, safety checklist |

---

## References

See [`bd_references.md`](bd_references.md) — 46 citations organised in 14 sections:

- Disease overview & epidemiology (1–4)
- Neurobiology & pathophysiology (5–9)
- Signal transduction (10–13)
- Neuroplasticity / BDNF (14–17)
- Neuroinflammation & HPA (18–21)
- Ion channels & circadian (22–25)
- Lithium PK/PD (26–28)
- Valproate pharmacology (29–30)
- Quetiapine & atypicals (31–33)
- Lamotrigine (34–35)
- Combination therapies (36–38)
- Biomarkers & precision medicine (39–41)
- QSP/computational modelling (42–44)
- Gut-brain axis & metabolic (45–46)

---

## Files

| File | Description |
|------|-------------|
| `bd_qsp_model.dot` | Graphviz source (120+ nodes, 12 clusters) |
| `bd_qsp_model.svg` | Scalable vector image (interactive zoom) |
| `bd_qsp_model.png` | Raster image (150 dpi) |
| `bd_mrgsolve_model.R` | mrgsolve ODE model + 6 scenarios + plots |
| `bd_shiny_app.R` | Shiny dashboard (6 tabs, no external mrgsolve dep) |
| `bd_references.md` | 46 PubMed-curated references |
| `README.md` | This file |
