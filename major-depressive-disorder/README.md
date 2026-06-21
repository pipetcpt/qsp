# Major Depressive Disorder (MDD) — QSP Model

**Date Added:** 2026-06-20  
**Category:** Neuropsychiatric Disorders  
**Model Version:** 1.0  
**Abbreviation:** `mdd`

---

## Disease Overview

Major Depressive Disorder (MDD) is a highly prevalent, recurrent psychiatric disorder affecting approximately 280 million people worldwide (WHO, 2023). It is characterized by persistent depressed mood, anhedonia, cognitive impairment, and somatic symptoms lasting ≥ 2 weeks, causing significant functional impairment. MDD is the leading cause of disability globally and contributes substantially to the global burden of disease.

Despite decades of research, only ~33% of patients achieve remission with first-line antidepressant therapy (STAR*D trial). Treatment-resistant depression (TRD), defined as failure of ≥ 2 adequate antidepressant trials, affects approximately 30% of patients.

---

## Mechanistic Map Preview

[![MDD QSP Mechanistic Map](mdd_qsp_model.png)](mdd_qsp_model.svg)

*Click to view full-resolution SVG. 10 subsystem clusters, 100+ nodes.*

---

## Pathophysiology Summary

| System | MDD Pathology | Key Mediators | Consequence |
|--------|--------------|---------------|-------------|
| **Serotonergic** | ↓ synaptic 5-HT, ↑ SERT expression | 5-HT, SERT, 5-HT1A autoreceptor | Depressed mood, anhedonia |
| **Noradrenergic** | ↓ NE (locus coeruleus), NET upregulation | NE, NET, α2 autoreceptor | Fatigue, cognitive impairment |
| **Dopaminergic** | ↓ mesolimbic DA, reward circuit hypoactivity | DA, DAT, D1/D2 receptors | Anhedonia, motivational deficit |
| **HPA Axis** | Hypercortisolaemia, GR resistance, blunted DST | CRH, ACTH, cortisol, GR/MR | Hippocampal atrophy, neurotoxicity |
| **BDNF/Neuroplasticity** | ↓ BDNF, ↓ TrkB signalling, ↓ neurogenesis | BDNF, TrkB, CREB, mTORC1 | Synaptic loss, hippocampal shrinkage |
| **Neuroinflammation** | Microglia activation, ↑ pro-inflammatory cytokines | IL-1β, IL-6, TNF-α, CRP, IDO1 | Kynurenine pathway shunt, 5-HT depletion |
| **Kynurenine Pathway** | ↑ IDO1/KMO, ↑ QUIN (NMDA agonist), ↓ KA | QUIN, KA, kynurenine, TRP depletion | Excitotoxicity, neuroinflammation |
| **Glutamate/NMDA** | Excess synaptic glutamate, ↓ AMPA plasticity | NMDA-R, AMPA-R, mGluR2/3/5 | Excitotoxicity, LTP impairment |
| **Circadian Rhythm** | CLOCK/BMAL1 disruption, ↓ melatonin signalling | PER1/2, CRY1/2, SCN, MT1/MT2 | Sleep architecture disruption, diurnal variation |

---

## Drug Mechanism Table

| Drug Class | Example | Primary Target | Key Mechanism |
|-----------|---------|----------------|---------------|
| SSRI | Escitalopram 10-20 mg/day | SERT (Ki=1.1 nM) | Block 5-HT reuptake → ↑ synaptic 5-HT |
| SNRI | Venlafaxine 75-225 mg/day | SERT + NET | Block 5-HT & NE reuptake → ↑ both monoamines |
| TCA | Amitriptyline 75-150 mg/day | SERT + NET + H1/M1/α1 | Broad monoamine + receptor block |
| MAOI | Phenelzine 15-45 mg/day | MAO-A/B (irreversible) | ↓ 5-HT/NE/DA catabolism |
| NDRI | Bupropion 150-300 mg/day | NET + DAT | ↑ NE & DA (no SERT effect) |
| NMDA antagonist | Ketamine IV 0.5 mg/kg | NMDA-R (Ki≈3 µM) | Block NMDA → disinhibit AMPA → ↑ BDNF → mTOR activation → rapid synaptogenesis |
| Esketamine | Spravato 56-84 mg IN | NMDA-R | S-enantiomer; FDA approved for TRD (2019) |
| Lithium | Lithium carbonate 600-1200 mg | GSK-3β, inositol | ↓ GSK-3β → ↑ β-catenin → neuroprotection; ↑ BDNF |
| Atypical antipsychotic | Aripiprazole 2-15 mg | D2 partial agonist, 5-HT2A ant. | Augmentation; ↑ PFC DA/NE |
| Melatonergic | Agomelatine 25-50 mg | MT1/MT2 agonist, 5-HT2C ant. | Restore circadian rhythm; ↑ DA/NE in PFC |

---

## Model Files

| File | Description |
|------|-------------|
| `mdd_qsp_model.dot` | Graphviz DOT: 10 clusters, 100+ nodes, full MDD mechanistic map |
| `mdd_qsp_model.svg` | Rendered SVG (interactive zoom) |
| `mdd_qsp_model.png` | Rendered PNG (150 dpi) |
| `mdd_mrgsolve_model.R` | mrgsolve ODE model: 18 compartments, 6 treatment scenarios |
| `mdd_shiny_app.R` | 6-tab Shiny dashboard for interactive QSP simulation |
| `mdd_references.md` | 60 PubMed references organized by section |

---

## Model Parameters

### Drug PK Parameters

| Drug | Ka (h⁻¹) | CL (L/h) | Vd (L) | F (%) | t½ (h) | Ki_SERT (nM) | Ki_NET (nM) |
|------|---------|---------|-------|-------|--------|------------|-----------|
| Escitalopram | 0.46 | 37 | 1090 | 80 | 27 | 1.1 | — |
| Venlafaxine | 1.20 | 90 | 500 | 45 | 5 | 7.5 | 2.7 |
| Bupropion | 0.35 | 200 | 3000 | 85 | 21 | — | 52 |
| Ketamine IV | — | 122 | 210 | (IV) | 0.4 | — | — |

### Neurotransmitter Dynamics

| Parameter | Value | Units | Description |
|-----------|-------|-------|-------------|
| KSYN_5HT | 0.12 | nM/h | 5-HT synthesis rate |
| KDEG_5HT | 0.08 | h⁻¹ | 5-HT MAO-A degradation |
| KREUP_5HT | 0.15 | h⁻¹ | SERT-mediated reuptake |
| SS_5HT | 1.5 | nM | Healthy synaptic 5-HT |
| KSYN_NE | 0.08 | nM/h | NE synthesis rate |
| KREUP_NE | 0.12 | h⁻¹ | NET-mediated reuptake |
| SS_NE | 1.0 | nM | Healthy synaptic NE |
| KSYN_DA | 0.06 | nM/h | DA synthesis rate |
| SS_DA | 0.8 | nM | Healthy synaptic DA |

### HPA Axis Parameters

| Parameter | Value | Units | Description |
|-----------|-------|-------|-------------|
| SS_CRH | 0.67 | pg/mL | Baseline CRH |
| SS_ACTH | 1.67 | pg/mL | Baseline ACTH |
| SS_CORT | 15.0 | nmol/L | Baseline morning cortisol |
| CORT_FB | 0.04 | (nmol/L)⁻¹ | Cortisol neg. feedback strength |
| MDD_CORT_EX | 1.40 | — | MDD cortisol excess factor |

### MDD Disease Parameters (Baseline Perturbations)

| Parameter | Value | Description |
|-----------|-------|-------------|
| MDD_5HT_DEF | 0.60 | 5-HT deficit: 60% of normal |
| MDD_NE_DEF | 0.70 | NE deficit: 70% of normal |
| MDD_DA_DEF | 0.65 | DA deficit: 65% of normal |
| MDD_BDNF_DEF | 0.65 | BDNF deficit: 65% of normal |
| MDD_CORT_EX | 1.40 | Cortisol 140% of healthy |
| MDD_IL6_EX | 1.80 | IL-6 180% of healthy |
| HDRS_BASE | 22.0 | Baseline HDRS-17 (moderate-severe) |

---

## Clinical Scenarios

| # | Scenario | Drugs | HDRS at Week 8 (predicted) | Notes |
|---|---------|-------|--------------------------|-------|
| 1 | No Treatment | — | ~22 | Natural disease progression |
| 2 | Escitalopram 10mg/day | SSRI | ~16 | SERT occupancy ~75% |
| 3 | Venlafaxine 150mg/day | SNRI | ~15 | SERT + NET dual block |
| 4 | Ketamine IV 0.5mg/kg | NMDA antagonist | ~14 | Rapid onset (24-72h), mTOR-BDNF |
| 5 | ESC + Aripiprazole | SSRI + D2 partial | ~13 | Augmentation strategy |
| 6 | Bupropion 300mg/day | NDRI | ~17 | NE/DA focus, no SERT |
| 7 | TRD High-Stress | ESC 20mg | ~20 | High inflammation, cortisol dysreg |

**Response criterion:** ≥ 50% reduction in HDRS-17 from baseline  
**Remission criterion:** PHQ-9 < 5 (HDRS-17 ≤ 7)

---

## Shiny App Tabs

1. **Patient Profile** — Demographics, PHQ-9/HDRS/MADRS inputs, comorbidities, biomarker radar chart
2. **Drug PK** — Plasma concentration–time curves, SERT/NET occupancy dynamics, dose-response curve
3. **Neurotransmitter Dynamics** — 5-HT, NE, DA synaptic concentrations over time; monoamine bar chart at week 8
4. **HPA & Neuroinflammation** — Cortisol, IL-6, BDNF, mTOR/neurogenesis index, kynurenine:tryptophan ratio
5. **Clinical Endpoints** — HDRS-17/PHQ-9/MADRS trajectories, response/remission tiles, weekly summary table
6. **Scenario Comparison** — Multi-arm HDRS comparison, forest plot of effect sizes, time-to-response bar chart, BDNF comparison

---

## How to Run

### Prerequisites

```r
install.packages(c("mrgsolve","shiny","shinydashboard","ggplot2",
                   "dplyr","tidyr","plotly","DT"))
```

### mrgsolve Simulation

```r
source("mdd_mrgsolve_model.R")
# Runs all 6 scenarios and produces 6 ggplot2 figures
```

### Shiny App

```r
shiny::runApp("mdd_shiny_app.R")
```

### Render DOT → Images

```bash
dot -Tsvg mdd_qsp_model.dot -o mdd_qsp_model.svg
dot -Tpng -Gdpi=150 mdd_qsp_model.dot -o mdd_qsp_model.png
```

---

## Key References

- Rush AJ et al. STAR*D (2006). *Am J Psychiatry* 163:1905. PMID: 17074942
- Zarate CA et al. Ketamine RCT (2006). *Arch Gen Psychiatry* 63:856. PMID: 16894061
- Cipriani A et al. Network meta-analysis (2018). *Lancet* 391:1357. PMID: 29477251
- Duman RS et al. Neurotrophic model (2006). *Biol Psychiatry* 59:1116. PMID: 16631126
- Dowlati Y et al. Cytokines meta-analysis (2010). *Biol Psychiatry* 67:446. PMID: 20015486
