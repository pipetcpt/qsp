# Duchenne Muscular Dystrophy (DMD) — QSP Model

**뒤시엔 근이영양증 정량적 시스템 약리학 모델**

---

## Disease Overview

| Feature | Details |
|---------|---------|
| **Full name** | Duchenne Muscular Dystrophy (DMD) |
| **OMIM** | #310200 |
| **Gene** | DMD (Xp21.2), largest human gene (2.4 Mb, 79 exons) |
| **Protein** | Dystrophin (427 kDa, rod-shaped cytoskeletal anchor) |
| **Inheritance** | X-linked recessive |
| **Prevalence** | ~1:3,500–5,000 live male births |
| **Onset** | Early childhood (2–5 years) |
| **LoA** | Median ~12yr untreated; ~15yr with corticosteroids |
| **Life expectancy** | Historically ~20s; now ≥30–40yr with modern care |

---

## Pathophysiology Summary

```
DMD gene mutation (frameshift/deletion)
        ↓
Dystrophin absent (<0.1% of normal)
        ↓
DAPC complex destabilized (sarcoglycans, dystroglycan, nNOS lost)
        ↓
Sarcolemmal fragility → microtears during contraction
        ↓
Ca²⁺ overload → calpain activation → necrosis
        ↓
ROS ↑ → DAMP release → TLR4 → NF-κB → M1 inflammation
        ↓
M2 macrophage TGF-β → myofibroblast → fibrosis ↑
        ↓
FAP fat infiltration + satellite cell exhaustion
        ↓
Progressive muscle function loss → LoA → respiratory failure
```

---

## Mechanistic Map

[![DMD QSP Mechanistic Map](dmd_qsp_model.png)](dmd_qsp_model.svg)

*12 subgraph clusters · 186+ nodes · Full PK/PD integration*

### Cluster Summary

| # | Cluster | Key Components |
|---|---------|---------------|
| 1 | Genetic Basis | DMD gene, 79 exons, deletion/nonsense/frameshift mutations, reading frame rule |
| 2 | DAPC Complex | Dystroglycan, sarcoglycans, syntrophins, nNOS, laminin-211 |
| 3 | Membrane Pathology | Sarcolemmal fragility, TRPC channels, Ca²⁺ influx, calpain activation |
| 4 | Oxidative Stress | ROS, mitochondrial dysfunction, Nrf2 antioxidant pathway |
| 5 | Inflammation | NF-κB, M1/M2 macrophages, TNF-α, IL-1β, IL-6, complement |
| 6 | Fibrosis | TGF-β1/Smad, CTGF, myofibroblasts, FAPs, collagen I/III |
| 7 | Regeneration | Satellite cells, Pax7/MyoD/myogenin, myostatin, IGF-1 |
| 8 | Corticosteroids | Deflazacort, prednisone, vamorolone, GR signaling, ADR profile |
| 9 | Exon-Skipping & Gene Therapy | Eteplirsen, golodirsen, casimersen, Elevidys (AAVrh74) |
| 10 | Novel Therapies | Givinostat, pamrevlumab, CRISPR, base editing, myostatin inhibitors |
| 11 | Clinical Endpoints | 6MWD, NSAA, FVC%, LVEF, CK, LoA, quality of life |
| 12 | Organ Systems | Cardiac, respiratory, GI, CNS (Dp71), bone |

---

## mrgsolve ODE Model

**File:** `dmd_mrgsolve_model.R`

### Model Architecture

| Category | Compartments |
|----------|-------------|
| **Corticosteroid PK** | DEPOT_CS, CENT_CS, PERIPH_CS |
| **ASO PK** | CENT_ASO, MUS_ASO, IC_ASO (intracellular) |
| **Gene Therapy** | AAV_CIRC, AAV_MUS |
| **Disease PD** | DYS (dystrophin), MEMI (membrane integrity), CAI (Ca²⁺), ROS |
| **Inflammation** | NFkB, M1 macrophages, M2 macrophages, TGFb |
| **Fibrosis/Regen** | FIB (fibrosis score), SC (satellite cells) |
| **Outcomes** | MF (muscle function), SWD (6MWD), FVC_pct, LVEF_pct |
| **Total** | **19 ODE compartments** |

### Treatment Scenarios

| # | Scenario | Drug/Dose |
|---|----------|-----------|
| 1 | **Natural History** | No treatment |
| 2 | **Deflazacort** | 0.9 mg/kg/day oral (FDA 2017) |
| 3 | **Prednisone** | 0.75 mg/kg/day oral |
| 4 | **Eteplirsen** | 30 mg/kg/wk IV (exon 51 skip) |
| 5 | **Gene Therapy** | Delandistrogene moxeparvovec 1.33×10¹⁴ vg/kg (single dose) |
| 6 | **Combination** | Deflazacort + Eteplirsen |

### Key Parameters (Calibrated From)

| Parameter | Value | Source |
|-----------|-------|--------|
| Deflazacort CL | 25 L/h (70 kg) | McDonald 2013 CINRG |
| ASO muscle t½ | ~4–6 wk | Mendell 2016 Ann Neurol |
| Gene therapy peak dystrophin | ~20–30% | Mendell 2020 JAMA Neurol |
| 6MWD natural decline | ~-6 m/yr (age 7–12) | Bello 2015 Neurology |
| FVC annual decline | ~4%/yr untreated | CINRG DNHS |

---

## Shiny App

**File:** `dmd_shiny_app.R`

### 6 Tabs

| Tab | Content |
|-----|---------|
| **1. Patient Profile** | Age, weight, genotype input; drug selection; disease overview; 6 value boxes |
| **2. PK Panel** | CS plasma conc., ASO plasma/intracellular, AAV muscle, PK summary table |
| **3. Dystrophin & Membrane** | Dystrophin%, membrane integrity, Ca²⁺, ROS |
| **4. Inflammation & Fibrosis** | NF-κB, M1/M2 macrophages, TGF-β1, fibrosis score, satellite cells |
| **5. Clinical Endpoints** | 6MWD, NSAA, FVC%, LVEF, serum CK, summary table |
| **6. Scenario Comparison** | All 6 scenarios overlay, bar chart at user-defined year, outcome table |

```r
# Launch the Shiny app
source("dmd_shiny_app.R")
```

---

## References

**File:** `dmd_references.md` — 57 PubMed-linked references

| Category | Count |
|----------|-------|
| Disease genetics | 6 |
| DAPC & membrane | 4 |
| Ca²⁺ dysregulation & ROS | 5 |
| Inflammation | 5 |
| Fibrosis | 4 |
| Satellite cells | 3 |
| Corticosteroids | 7 |
| Exon-skipping ASO | 5 |
| Gene therapy | 3 |
| Novel therapies | 5 |
| Clinical/Natural history | 5 |
| Cardiac | 2 |
| QSP/Modeling | 4 |
| **Total** | **57** |

---

## Approved & Late-Stage Drugs (2024)

| Drug | Class | Target | Approval | Patients |
|------|-------|--------|----------|---------|
| Deflazacort (Emflaza®) | Corticosteroid | GR, NF-κB | FDA 2017 | All DMD ≥2yr |
| Vamorolone (Agamree®) | Dissociative steroid | GR (selective) | FDA 2023 | ≥2yr |
| Eteplirsen (Exondys 51®) | PMO ASO | Exon 51 splicing | FDA 2016 | ~13% (exon 51 skip) |
| Golodirsen (Vyondys 53®) | PMO ASO | Exon 53 splicing | FDA 2019 | ~8% (exon 53 skip) |
| Viltolarsen (Viltepso®) | PMO ASO | Exon 53 splicing | FDA 2020 | ~8% (exon 53 skip) |
| Casimersen (Amondys 45®) | PMO ASO | Exon 45 splicing | FDA 2021 | ~8% (exon 45 skip) |
| Delandistrogene moxeparvovec (Elevidys®) | AAV gene therapy | Micro-dystrophin | FDA 2023 | 4–5yr (SRP-9001) |
| Givinostat (Duvyzat®) | HDAC inhibitor | HDAC I/II | FDA 2024 | ≥6yr |
| Ataluren (Translarna®) | Read-through | Stop codons (PTC) | EMA 2014 | ~10–15% nonsense |

---

## File Listing

```
duchenne-muscular-dystrophy/
├── README.md                   ← This file
├── dmd_qsp_model.dot           ← Graphviz mechanistic map source
├── dmd_qsp_model.svg           ← Interactive SVG (186+ nodes)
├── dmd_qsp_model.png           ← 150 dpi PNG thumbnail
├── dmd_mrgsolve_model.R        ← 19-compartment ODE model + 6 scenarios
├── dmd_shiny_app.R             ← 6-tab interactive Shiny dashboard
└── dmd_references.md           ← 57 PubMed references
```

---

*Model version: v1.0 | Created: 2026-06-25 | Part of the QSP Disease Model Library*
