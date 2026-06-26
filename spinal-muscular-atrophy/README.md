# Spinal Muscular Atrophy (SMA) — QSP Model

## Overview

Spinal Muscular Atrophy (SMA) is a severe autosomal recessive neuromuscular disease caused by homozygous deletion or mutation of the *SMN1* gene on chromosome 5q13.2. Loss of full-length Survival Motor Neuron (SMN) protein leads to progressive alpha motor neuron degeneration, denervation of skeletal muscle, and in severe cases respiratory failure and death in infancy.

This QSP model captures the complete mechanistic cascade from SMN2 pre-mRNA alternative splicing through motor neuron degeneration and skeletal muscle atrophy, integrated with three FDA-approved disease-modifying therapies.

**ICD-10:** G12.0 (infantile SMA), G12.1 (other inherited SMA), G12.9 (SMA unspecified)

---

## Disease Mechanism

```
SMN1 deletion (5q13.2)
        │
        ▼
SMN2 gene (backup)                 SMN2 copies: 1–4 (higher = milder)
        │
        ├─ Pre-mRNA alternative splicing ──────────────────────────────────┐
        │   ├─ C840T transition disrupts exon 7 ESE                        │
        │   ├─ ISS-N1 (intronic silencer) represses exon 7 inclusion       │
        │   └─ Baseline 10% exon 7 inclusion → 10% FL-SMN mRNA            │
        │                                                                   │
        ├─ FL-SMN mRNA → FL-SMN protein (full function)          Drug target: ISS-N1 & ESE
        └─ SMN-Δ7 mRNA → SMN-Δ7 protein (rapidly degraded)
                │
                ▼
         SMN protein pool < threshold
                │
                ▼
        snRNP assembly defect → mRNA splicing dysregulation
        Axonal transport defect → β-actin mRNA mis-targeting
                │
                ▼
         Alpha motor neuron degeneration (apoptosis)
                │
                ▼
         Neuromuscular junction immaturity + denervation
                │
                ▼
         Skeletal muscle atrophy (neurogenic)
                │
                ▼
         Loss of motor milestones → Respiratory failure → Death
```

---

## QSP Model Structure

### Biological Compartments

| Module | Components | Key Variables |
|--------|-----------|---------------|
| SMN2 Splicing | Exon 7 ESE/ESS, ISS-N1, hnRNP A1, SRSF1 | E7I (exon-7 inclusion fraction) |
| SMN mRNA | FL-SMN mRNA, SMN-Δ7 mRNA | mRNA concentrations |
| SMN Protein | FL-SMN pool, snRNP complex | SMN protein (normalized) |
| Motor Neuron | Alpha-MN pool, apoptotic cascade | MN fraction (0–1) |
| NMJ | Presynaptic terminal, AChR clustering | NMJ maturity score (0–1) |
| Muscle | Fiber mass, atrophy/reinnervation | Muscle mass (normalized) |
| Clinical | CMAP, HFMSE, CHOP-INTEND, FVC | Clinical endpoint scores |

### Drug Modules (PK/PD)

| Drug | Route | Mechanism | PK Model |
|------|-------|-----------|----------|
| **Nusinersen** (Spinraza) | Intrathecal | ASO blocks ISS-N1, ↑ exon 7 inclusion | 2-compartment CSF (lumbar→cervical→CNS tissue) |
| **Risdiplam** (Evrysdi) | Oral daily | Small molecule enhances SRSF1/Tra2β binding | 1-compartment + CNS penetration |
| **Onasemnogene** (Zolgensma) | IV single dose | AAV9 delivers SMN1 transgene to MNs | IV → plasma → BBB crossing → MN transduction |

---

## Files

| File | Description |
|------|-------------|
| `sma_qsp_model.dot` | Graphviz mechanistic map source (130+ nodes, 12 clusters) |
| `sma_qsp_model.svg` | Vector image of mechanistic map |
| `sma_qsp_model.png` | Raster image at 150 dpi |
| `sma_mrgsolve_model.R` | mrgsolve ODE model with 20 compartments, 6 scenarios |
| `sma_shiny_app.R` | 8-tab interactive Shiny dashboard |
| `sma_references.md` | 50 PubMed references with clinical trial summary |

---

## Mechanistic Map

[![SMA QSP Mechanistic Map](sma_qsp_model.png)](sma_qsp_model.svg)

The mechanistic map contains **130+ nodes** organized in **12 subgraph clusters**:
1. Genetic & Chromosomal Basis (5q13.2)
2. SMN2 Pre-mRNA Alternative Splicing
3. SMN Protein & snRNP Assembly Complex
4. Motor Neuron Biology & Degeneration
5. Neuromuscular Junction (NMJ) Biology
6. Skeletal Muscle Pathology
7. Clinical Endpoints & Biomarkers
8. Nusinersen (Spinraza) PK/PD
9. Risdiplam (Evrysdi) PK/PD
10. Onasemnogene (Zolgensma) PK/PD
11. SMA Disease Types & Natural History
12. Systemic Effects & Comorbidities

---

## mrgsolve ODE Model

### Compartments (20 ODE states)

**Nusinersen PK (3):** CSF lumbar, CSF cervical, CNS tissue  
**Risdiplam PK (3):** GI tract, plasma, CNS  
**Zolgensma PK (3):** Plasma vector, transduced MN vg load, transgene mRNA  
**Disease Biology (6):** FL-SMN mRNA, SMN-Δ7 mRNA, SMN protein, MN pool, NMJ score, Muscle mass  
**Cumulative (2):** AUC-SMN, MN lost

### Treatment Scenarios (6)

| Scenario | Drug | SMA Type | Clinical Trial Analogue |
|----------|------|----------|------------------------|
| 1 | None | Type I | Natural history |
| 2 | Nusinersen | Type I | ENDEAR trial |
| 3 | Risdiplam | Type II | SUNFISH trial |
| 4 | Zolgensma | Presymptomatic | SPR1NT trial |
| 5 | Nusinersen late | Type II | Late-start treatment |
| 6 | Risdiplam pediatric | Type II | Weight-based dosing |

### Running the Model

```r
library(mrgsolve)
library(dplyr)
library(ggplot2)

source("sma_mrgsolve_model.R")

# Build model
mod <- sma_model()

# Run single scenario: nusinersen in SMA Type I
mod_t1 <- param(mod, list(SMN2_copies = 2, k_MN_death = 0.004))
out <- mrgsim(mod_t1, ev_nusinersen(), delta = 1, end = 730)
plot(out, SMN_protein + MN_fraction + CMAP + CHOP_INTEND ~ time)

# Run all 6 scenarios
df_all <- run_scenarios()

# Virtual population (n=100)
vpc <- virtual_population(n = 100)
```

---

## Shiny App (8 Tabs)

| Tab | Content |
|-----|---------|
| 1. Patient Profile | Disease type selector, SMA classification table, mechanistic map |
| 2. Pharmacokinetics | Drug-specific PK parameters, concentration-time plots |
| 3. SMN Biology (PD) | Exon-7 inclusion dynamics, SMN protein, dose-response curves |
| 4. Motor Neuron & NMJ | MN pool dynamics, NMJ maturation, muscle mass |
| 5. Clinical Endpoints | CMAP, HFMSE, CHOP-INTEND, FVC, RULM value boxes + plots |
| 6. Scenario Comparison | Multi-treatment comparison, endpoint selection, summary table |
| 7. Biomarker Dashboard | NF-L proxy, SMN protein, correlation matrix |
| 8. Population Variability | Virtual population simulation, VPC plots, responder analysis |

### Launch App

```r
library(shiny)
runApp("sma_shiny_app.R")
```

---

## Key Parameters

### SMN2 Splicing

| Parameter | Value | Description |
|-----------|-------|-------------|
| E7I_base | 0.10 | Baseline exon-7 inclusion (10%) |
| Emax_NUS | 0.60 | Nusinersen max ΔE7I |
| EC50_NUS | 5.0 ng/g | Nusinersen EC50 (CNS tissue) |
| Emax_RIS | 0.50 | Risdiplam max ΔE7I |
| EC50_RIS | 80 ng/mL | Risdiplam EC50 (plasma) |

### Motor Neuron Dynamics

| Parameter | Value | Description |
|-----------|-------|-------------|
| k_MN_death | 0.002–0.004/day | MN death rate (type-dependent) |
| SMN_thresh | 0.30 | SMN protein threshold |
| MN_min | 0.05 | Irreducible MN fraction |

### Drug PK (Nusinersen)

| Parameter | Value | Description |
|-----------|-------|-------------|
| Dose | 12 mg IT | Standard adult/pediatric |
| Loading | Day 0,14,28,63 | 4 loading injections |
| Maintenance | Every 4 months | Ongoing |
| t½ (CNS) | 135–177 days | Long tissue half-life |

---

## Clinical Trial Calibration

Model parameters calibrated to reproduce key outcomes from:

- **ENDEAR** (NCT02193074): 51% motor milestone response vs 0% sham at Day 183
- **CHERISH** (NCT02292537): +4.0 HFMSE (nusinersen) vs −1.9 (sham) at 15 months
- **FIREFISH** (NCT02913482): 61% sitting ≥5 s at 12 months (risdiplam)
- **STR1VE** (NCT03461289): 91% alive without permanent ventilation at 14 months (Zolgensma)

---

## Disease Classification Reference

| SMA Type | Onset | SMN2 Copies | Maximum Ability | Untreated Survival |
|----------|-------|-------------|-----------------|-------------------|
| 0 | Prenatal/neonatal | 1 | None | <1 month |
| I (Werdnig-Hoffmann) | <6 months | 1–2 | Never sit | <2 years |
| II (Intermediate) | 6–18 months | 3 | Sit, never stand | Adult (reduced) |
| III (Kugelberg-Welander) | >18 months | 3–4 | Walk (may lose) | Normal |
| IV (Adult) | >21 years | ≥4 | Walk throughout | Normal |

---

## References

See [`sma_references.md`](sma_references.md) for 50 curated PubMed references covering:
- Disease pathophysiology & genetics (1–8)
- SMN protein biology & splicing (9–14)
- Motor neuron biology (15–18)
- NMJ pathology (19–21)
- Nusinersen clinical trials & PK (22–27)
- Risdiplam clinical trials & PK (28–32)
- Zolgensma gene therapy (33–37)
- Biomarkers (38–41)
- QSP & pharmacometric models (42–45)
- Supportive care & outcomes (46–50)

---

*Model built: 2026-06-23 | Disease: Spinal Muscular Atrophy | Category: Neuromuscular*
