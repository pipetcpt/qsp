# Acute Intermittent Porphyria (AIP) QSP Model

[![Disease](https://img.shields.io/badge/Disease-Acute%20Intermittent%20Porphyria-8B0000)](https://www.ncbi.nlm.nih.gov/medgen/10523)
[![Category](https://img.shields.io/badge/Category-Rare%20Metabolic%20Disease-orange)](.)
[![Gene](https://img.shields.io/badge/Gene-HMBS%20(PBGD)-purple)](https://www.omim.org/entry/176000)
[![Drug](https://img.shields.io/badge/Drug-Givosiran%20%7C%20Hemin%20IV-blue)](.)
[![References](https://img.shields.io/badge/References-57-brightgreen)](aip_references.md)

## Overview

Acute Intermittent Porphyria (AIP; OMIM #176000) is an autosomal dominant metabolic disorder caused by **≥50% deficiency of hydroxymethylbilane synthase (HMBS/PBGD)**, the third enzyme in the heme biosynthesis pathway. The hallmark is episodic accumulation of the neurotoxic heme precursor **δ-aminolevulinic acid (ALA)** and **porphobilinogen (PBG)** in the liver, triggering acute neurovisceral attacks.

This QSP model integrates:
- **Mechanistic heme biosynthesis pathway** (ALAS1 → ALA → PBG bottleneck → heme)
- **Givosiran siRNA PK/PD** (GalNAc-ASGPR hepatic delivery, RISC-mediated ALAS1 knockdown)
- **Hemin IV PK/PD** (hemopexin/albumin transport, HO-1 catabolism, ALAS1 feedback)
- **Neurotoxicity cascade** (GABA-A inhibition, oxidative stress, axonal degeneration)
- **Hormonal triggers** (menstrual cycle progesterone upregulation of ALAS1)
- **Clinical endpoints** (AAR, urinary ALA/PBG normalization, annual hospitalization rate)

---

## Mechanistic Map

[![AIP QSP Map](aip_qsp_model.png)](aip_qsp_model.svg)

*Click image to view full SVG. 11 subgraph clusters, 130+ nodes spanning heme biosynthesis, regulation, triggers, givosiran PK/PD, hemin PK/PD, neurotoxicity, and clinical endpoints.*

---

## Disease Pathophysiology

```
Trigger (drugs/hormones/fasting)
        │
        ↓
ALAS1 mRNA ↑↑↑  (rate-limiting enzyme)
        │
        ↓
ALA accumulation
        │
        ↓  (PBGD/HMBS: 50% deficient → bottleneck)
PBG accumulation ──────────────────────────────────────────┐
        │                                                   │
        ↓                                                   ↓
ALA neurotoxicity (GABA-A inhibition, ROS)          Dark urine (porphyrins)
        │
        ↓
Visceral pain · Tachycardia · Seizures · Paralysis
(Annual Attack Rate ↑ in untreated: ~10-15/year)
```

### Key Enzyme Deficiency: PBGD (50% residual activity)

| State | PBGD Activity | ALA (× normal) | PBG (× normal) |
|-------|--------------|----------------|----------------|
| Normal (wildtype) | 100% | 1× | 1× |
| AIP carrier (remission) | 50% | 2-5× | 3-10× |
| AIP acute attack | 50% | >20× | >50-200× |
| Post-givosiran | 50% + ALAS1↓ | <2× | <3× |

---

## Model Deliverables

| File | Description |
|------|-------------|
| [`aip_qsp_model.dot`](aip_qsp_model.dot) | Graphviz mechanistic map (130+ nodes, 11 clusters) |
| [`aip_qsp_model.svg`](aip_qsp_model.svg) | Scalable vector map |
| [`aip_qsp_model.png`](aip_qsp_model.png) | PNG thumbnail (150 dpi) |
| [`aip_mrgsolve_model.R`](aip_mrgsolve_model.R) | mrgsolve ODE model (17 compartments, 6 scenarios, VPop N=100) |
| [`aip_shiny_app.R`](aip_shiny_app.R) | Shiny dashboard (6 tabs) |
| [`aip_references.md`](aip_references.md) | 57 PubMed-cited references |

---

## mrgsolve Model Summary

### Compartments (17)

| # | Compartment | Description |
|---|-------------|-------------|
| 1 | `GIV_SC` | Givosiran SC depot (mg) |
| 2 | `GIV_C` | Givosiran plasma central (µg/L) |
| 3 | `GIV_P` | Givosiran peripheral (µg/L) |
| 4 | `GIV_LIV` | Givosiran liver (ng/g) |
| 5 | `ALAS1_mRNA` | ALAS1 mRNA (relative expression) |
| 6 | `ALAS1_PROT` | ALAS1 protein/activity (relative) |
| 7 | `ALA_LIV` | Hepatic ALA (nmol/g liver) |
| 8 | `ALA_PLAS` | Plasma ALA (µmol/L) |
| 9 | `PBG_LIV` | Hepatic PBG (nmol/g) |
| 10 | `PBG_PLAS` | Plasma PBG (µmol/L) |
| 11 | `HEME_LIV` | Hepatic free heme pool (nmol/g) |
| 12 | `HEM_C` | Hemin IV plasma (µg/L) |
| 13 | `HEM_LIV` | Hemin liver (nmol/g) |
| 14 | `NEUROTOX` | Neurotoxicity index |
| 15 | `ATK_DAY` | Attack-risk integral |
| 16 | `AUC_ALA` | Cumulative plasma ALA AUC |
| 17 | `AUC_PBG` | Cumulative plasma PBG AUC |

### Treatment Scenarios

| Scenario | Description |
|----------|-------------|
| 1 | Placebo (natural disease history, hormonal cycling) |
| 2 | Givosiran 2.5 mg/kg SC Q28d (standard ENVISION dose) |
| 3 | Hemin IV 3 mg/kg/d × 4 days (acute attack treatment) |
| 4 | Givosiran prophylaxis + breakthrough hemin at D90 |
| 5 | Givosiran 5.0 mg/kg Q28d (exploratory high dose) |
| 6 | Gene therapy – PBGD activity restoration to 95% |

### Clinical Trial Calibration

| Trial | Calibration Target |
|-------|-------------------|
| ENVISION (Balwani 2020 NEJM) | Givosiran 2.5 mg/kg: 74% AAR reduction vs placebo |
| ENVISION Month 6 | Urine ALA normalization: 73% of patients |
| ENVISION Month 6 | Urine PBG normalization: 63% of patients |
| Sardh 2019 NEJM | ALAS1 mRNA knockdown: ~87% at Month 3 trough |

---

## Shiny App Tabs

| Tab | Content |
|-----|---------|
| 1. Patient Profile | Demographics, HMBS genotype (PBGD activity), trigger factors, IV glucose |
| 2. Drug PK | Givosiran plasma/liver PK, Hemin IV plasma, siRNA knockdown efficiency |
| 3. PD Markers | ALAS1 mRNA/protein, ALA/PBG fold-change, heme pool, hormonal trigger |
| 4. Clinical Endpoints | AAR, urine ALA/PBG normalization, neurotoxicity index, attack-days |
| 5. Scenario Comparison | 6 treatment arms: time-course, AAR bar chart, spider plot, summary table |
| 6. Biomarkers & VPop | Virtual population simulation, normalization rates, E-R curve, eGFR trend |

---

## Key Pharmacology: Givosiran (GalNAc-siRNA)

```
SC injection (2.5 mg/kg Q28d)
    │
    ↓ absorption (t½ SC ~10h)
Plasma (2-compartment PK)
    │
    ↓ GalNAc–ASGPR receptor-mediated endocytosis
Hepatocytes (liver Cmax ~10,000 ng/g tissue)
    │
    ↓ RISC complex loading (Ago2)
ALAS1 mRNA catalytic cleavage (~87% knockdown at trough)
    │
    ↓
ALAS1 protein ↓ → ALA ↓ → PBG ↓ → AAR ↓74%
```

**FDA approval:** November 2019 | **Indication:** Adults with AIP

---

## Running the Model

### Prerequisites

```r
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "patchwork"))
# Shiny app also requires:
install.packages(c("shiny", "shinydashboard", "plotly", "DT"))
```

### Quick Start

```r
# Run all 6 scenarios + VPop
source("aip_mrgsolve_model.R")

# Launch interactive dashboard
shiny::runApp("aip_shiny_app.R")
```

### Render mechanistic map

```bash
dot -Tsvg aip_qsp_model.dot -o aip_qsp_model.svg
dot -Tpng -Gdpi=150 aip_qsp_model.dot -o aip_qsp_model.png
```

---

## References

See [`aip_references.md`](aip_references.md) for 57 PubMed citations organized in 12 sections:

1. Landmark Clinical Trials (ENVISION, Phase 1/2)
2. Heme Biosynthesis Pathway
3. ALAS1 Regulation
4. Givosiran Pharmacology
5. Hemin PK/PD
6. Neurotoxicity Mechanisms (ALA/PBG)
7. Epidemiology & Genetics
8. QSP/PK-PD Modeling
9. Long-term Complications
10. Hormonal Triggers
11. Gene Therapy
12. Drug Safety Databases

---

*Model generated: 2026-06-25 | Claude Code Routine (CCR)*
