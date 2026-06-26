# ADHD QSP Model
## Attention Deficit Hyperactivity Disorder (ADHD)
### Quantitative Systems Pharmacology Model

> **Directory:** `adhd/` | **Abbreviation:** ADHD | **Date:** 2026-06-26  
> **Category:** Neuropsychiatric | **ICD-10:** F90 | **DSM-5:** 314.0x  
> **Prevalence:** ~5–7% children · ~2.5% adults worldwide

[![ADHD QSP Mechanistic Map](adhd_qsp_model.png)](adhd_qsp_model.svg)

---

## Disease Overview

Attention Deficit Hyperactivity Disorder (ADHD) is a neurodevelopmental disorder characterized by persistent patterns of inattention, hyperactivity, and impulsivity that interfere with functioning and development. It affects approximately 5–7% of school-age children and 2.5% of adults globally, making it one of the most prevalent psychiatric conditions across the lifespan.

The core pathophysiology involves dysregulation of catecholaminergic neurotransmission — primarily dopamine (DA) and norepinephrine (NE) — in prefrontal cortical circuits governing executive function, working memory, and inhibitory control. Genetic studies implicate multiple risk variants (DAT1, DRD4, DRD5, SNAP25, COMT) that collectively reduce dopaminergic tone in frontostriatal pathways, leading to the clinical syndrome.

---

## Core Pathophysiology

```
Genetic Risk Factors (DAT1, DRD4, COMT, SNAP25, BDNF Val66Met)
  │
  ├─► Reduced DA synthesis/release in VTA→PFC (mesocortical pathway)
  ├─► Reduced NE tone from LC→PFC projections
  └─► Delayed PFC cortical maturation (~3 year lag vs. neurotypical controls)
        │
        ▼
PFC Dopamine↓ + NE↓ → Suboptimal catecholamine tone
  │
  ├─► D1R understimulation → Weakened AMPAR synapses → Working Memory↓
  ├─► α2A-AR understimulation → HCN channel open → Ih current↑ → Network instability
  └─► Prefrontal inhibitory control ↓, Default Mode Network overactive
        │
        ▼
Clinical ADHD Syndrome:
  ├─► Inattention (ADHD-RS Inattention subscale ≥6/9 symptoms)
  ├─► Hyperactivity / Impulsivity (Hyperactivity subscale ≥6/9)
  ├─► Working memory deficits (BRIEF-2, digit span)
  └─► Executive function impairment (WCST, Stop-Signal Task)

Drug Treatment:
  Stimulants: MPH/AMP → DAT/NET inhibition + DA efflux → Synaptic DA/NE ↑
  Non-stimulants: ATX → NET selective inhibition → NE PFC ↑
                  GFN → α2A-AR agonist → HCN closure → PFC network ↑
                  VLX → NET inhibition + 5-HT modulation
```

---

## QSP Model Architecture

### Mechanistic Map
- **~140 nodes** · **9 subgraph clusters**
- Clusters: ① Genetic & Neurodevelopmental · ② Catecholamine Biosynthesis · ③ Dopamine System · ④ Norepinephrine System · ⑤ PFC Circuits & Cognition · ⑥ Pharmacokinetics · ⑦ Pharmacodynamics · ⑧ Clinical Endpoints · ⑨ Biomarkers

### mrgsolve ODE Model — 25 Compartments

| Module | Compartments | Description |
|--------|-------------|-------------|
| MPH PK | GUT1, CENT1, PER1 | 2-compartment oral PK |
| AMP PK | GUT2, CENT2, PER2 | 2-compartment oral PK |
| ATX PK | GUT3, CENT3, PER3 | CYP2D6-dependent clearance |
| GFN PK | GUT4, CENT4, PER4 | Extensive distribution (Vd=16 L/kg) |
| VLX PK | GUT5, CENT5 | 1-compartment |
| Neurotransmission | DA_syn, NE_syn | Synaptic catecholamines |
| Transporter Occupancy | DAT_occ, NET_occ | Competitive inhibition |
| PFC Tone | PFC_DA, PFC_NE | Inverted-U function |
| Cognitive | WM_idx, ExecFun | Working memory & executive function |
| Clinical | ADHD_RS, CGI_S, QoL_idx | Rating scale outcomes |

### Treatment Scenarios (7)

| # | Scenario | Drug | Dose | Mechanism |
|---|---------|------|------|-----------|
| 1 | Untreated | — | — | Natural history |
| 2 | MPH IR TID | Methylphenidate IR | 10 mg TID | DAT/NET inhibitor |
| 3 | MPH ER QD | Methylphenidate ER | 36 mg QD | DAT/NET (slow-release) |
| 4 | AMP XR QD | Amphetamine XR | 20 mg QD | DAT reverse transport |
| 5 | ATX QD | Atomoxetine | 80 mg QD | Selective NET inhibitor |
| 6 | GFN ER QD | Guanfacine ER | 4 mg QD | α2A-AR agonist |
| 7 | VLX ER QD | Viloxazine ER | 400 mg QD | NET + 5-HT modulator |

### Shiny Dashboard — 7 Tabs

| Tab | Content |
|-----|---------|
| ① Patient Profile | Demographics, ADHD subtype, comorbidities, genetic risk |
| ② Drug PK | Concentration–time profile, DAT/NET occupancy, PK parameters |
| ③ DA/NE Dynamics | Synaptic neurotransmitter levels, Inverted-U visualization |
| ④ PFC & Cognition | PFC tone indices, working memory gauge, executive function |
| ⑤ Clinical Endpoints | ADHD-RS-5 trajectory, CGI-S, response/remission rates |
| ⑥ Scenario Comparison | All 7 treatment arms, summary table (DataTable) |
| ⑦ Biomarker Panel | Multi-biomarker facet plot, DAT/NET vs. DA/NE correlations |

---

## Key Pharmacological Parameters

### PK Parameters (Reference Adults, 70 kg)

| Drug | t½ | Tmax | F% | CL (L/h) | Vc (L) |
|------|----|------|-----|----------|--------|
| MPH IR | 2.5 h | 1.5 h | 22% | 31.5 | 448 |
| MPH ER | ~8 h | 6–8 h | 22% | 31.5 | 448 |
| AMP XR | 9–14 h | 7 h | 75% | 39.2 | 245 |
| ATX | 5 h (EM) | 1–2 h | 63% | 24.5 | 59.5 |
| GFN ER | 17 h | 5 h | 80% | 3.5 | 196 |
| VLX ER | 7 h | 5 h | 88% | 19.6 | 105 |

### PD Parameters (Transporter Affinity)

| Drug | Ki(DAT) | Ki(NET) | Primary Target |
|------|---------|---------|---------------|
| MPH | 34 nM | 340 nM | DAT (primary) |
| AMP | ~100 nM | 40 nM | DAT reverse transport |
| ATX | — | 2 nM | NET (selective) |
| GFN | — | — | α2A-AR (1 nM) |
| VLX | — | 42 nM | NET + 5-HT2B |

---

## Clinical Trial Calibration

| Clinical Trial | Drug | Key Result | Reference |
|---------------|------|-----------|-----------|
| MTA Study (n=579) | MPH | ADHD-RS reduction ~10 pts vs. placebo | Arch Gen Psychiatry 1999 |
| Biederman 2002 (n=584) | AMP XR 20mg | ADHD-RS reduction ~12 pts | JAACAP 2002 |
| Michelson 2001 (n=297) | ATX 80mg | ADHD-RS reduction ~8 pts (EM) | Pediatrics 2001 |
| Sallee 2009 (n=345) | GFN ER 4mg | ADHD-RS reduction ~7 pts | JAACAP 2009 |
| Nasser 2021 (n=460) | VLX ER 400mg | ADHD-RS reduction ~8 pts | Neuropsychiatr 2021 |
| Volkow 1998 (PET) | MPH | DAT occupancy 72% at therapeutic dose | Am J Psychiatry 1998 |

---

## Running the Model

### Prerequisites
```r
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr",
                   "shiny", "shinydashboard", "plotly", "DT"))
```

### Run mrgsolve simulation
```r
source("adhd_mrgsolve_model.R")
```

### Launch Shiny app
```r
shiny::runApp("adhd_shiny_app.R")
```

### Render mechanistic map
```bash
dot -Tsvg adhd_qsp_model.dot -o adhd_qsp_model.svg
dot -Tpng -Gdpi=150 adhd_qsp_model.dot -o adhd_qsp_model.png
```

---

## File List

| File | Description |
|------|-------------|
| [`adhd_qsp_model.dot`](adhd_qsp_model.dot) | Graphviz mechanistic map source (~140 nodes, 9 clusters) |
| [`adhd_qsp_model.svg`](adhd_qsp_model.svg) | Vector mechanistic map (high resolution) |
| [`adhd_qsp_model.png`](adhd_qsp_model.png) | Raster map (150 dpi) |
| [`adhd_mrgsolve_model.R`](adhd_mrgsolve_model.R) | mrgsolve ODE QSP model (25 compartments, 7 scenarios) |
| [`adhd_shiny_app.R`](adhd_shiny_app.R) | Shiny interactive dashboard (7 tabs) |
| [`adhd_references.md`](adhd_references.md) | References: 53 PubMed citations (12 sections) |

---

## Key Concepts Modeled

1. **Inverted-U function**: PFC function is optimal at intermediate DA/NE tone — both too little (ADHD) and too much (excess stimulant) impair cognition
2. **CYP2D6 pharmacogenomics**: ATX clearance varies 4× between extensive (EM) and poor metabolizers (PM), affecting efficacy and side-effect profiles  
3. **Allometric scaling**: PK parameters scaled by body weight (WT^0.75 for clearance), enabling pediatric vs. adult comparisons
4. **DAT/NET occupancy**: Competitive inhibition modeled explicitly; therapeutic DAT occupancy for MPH is 50–80% (Volkow 1998)
5. **Comorbidity weighting**: Coexisting anxiety, depression, ODD are tracked as modulators of functional outcomes

---

*Model built by Claude Code Routine — 2026-06-26*
