# Gaucher Disease — QSP Model

[![Mechanistic Map](gcd_qsp_model.png)](gcd_qsp_model.svg)

## Overview

**Gaucher disease** is the most common lysosomal storage disorder, caused by biallelic pathogenic variants in **GBA1** (glucocerebrosidase gene, 1q22). Deficient lysosomal β-glucocerebrosidase (GBA) leads to accumulation of glucocerebroside (GC) and its highly toxic deacylated derivative **lyso-glucocerebroside (lyso-GL1 / GL-1)** primarily within tissue macrophages of the spleen, liver, and bone marrow — forming characteristic **Gaucher cells** (crumpled tissue-paper macrophages).

### Disease Classification

| Type | Involvement | GBA Residual Activity | Key Features |
|------|-------------|----------------------|--------------|
| **Type 1** (non-neuronopathic) | Visceral · Haematological · Skeletal | 1–15% | Splenomegaly, anaemia, thrombocytopenia, bone disease; most common (95%) |
| **Type 2** (acute neuronopathic) | + Acute CNS | <1% | Rapid progressive neurodegeneration; fatal by age 2–3 |
| **Type 3** (chronic neuronopathic) | + Chronic CNS | 1–10% | Oculomotor apraxia, seizures, cognitive decline; slower course |

### Treatment Approaches

| Strategy | Mechanism | Agents |
|----------|-----------|--------|
| **ERT** | M6P receptor → lysosomal GBA delivery | Imiglucerase, Velaglucerase α, Taliglucerase α |
| **SRT** | GCS inhibition → substrate synthesis ↓ | Eliglustat (selective), Miglustat (broad), Venglustat (CNS+) |
| **Pharmacological Chaperone** | Stabilise misfolded GBA → rescue from ERAD | Ambroxol (off-label, N370S-responsive) |
| **Gene Therapy** | Lentiviral GBA1 (phase 1–2) | In development |

---

## Model Components

### 1. Mechanistic Map (`gcd_qsp_model.dot / .svg / .png`)

| Cluster | Key Components |
|---------|---------------|
| Drug PK | ERT (2-compartment IV), Eliglustat / Miglustat / Venglustat (oral 1-comp) |
| Enzyme Biology | GBA gene → mRNA → ER glycosylation → M6P → lysosomal GBA; mutant ERAD; ERT uptake |
| Substrate | GCS synthesis → GC_macrophage → spleen / liver / bone marrow / CNS; lyso-GL1 biomarker |
| Gaucher Cell | Naive monocyte → tissue macrophage → GC engulfment → NF-κB activation |
| Cytokines | IL-1β, IL-6, TNF-α, MIP-1α, RANKL, OPG, chitotriosidase, CCL18, ferritin |
| Organ Manifestations | Splenomegaly, hepatomegaly, anaemia, thrombocytopenia, hepatic fibrosis |
| Bone | OB/OC balance, BMD, osteonecrosis |
| Neurological | CNS GBA deficiency, microglial activation, neurodegeneration (Type 2/3); GBA–Parkinson link |
| Hematopoiesis | HSC → erythropoiesis / megakaryopoiesis, bone marrow competition |
| PD Endpoints | GL-1 reduction, lyso-GL1, chitotriosidase, spleen/liver volume, Hb, platelets, BMD |

> **Map statistics:** ~115 nodes, 10 clusters, covers ERT + SRT + chaperone + gene therapy PD

---

### 2. mrgsolve ODE Model (`gcd_mrgsolve_model.R`)

#### Compartments (26 ODEs)

| Category | Compartments |
|----------|-------------|
| Drug PK | ERT_C, ERT_T (2-comp IV); ELIS_GUT/C; MIGS_GUT/C; VENG_GUT/C |
| Enzyme / Substrate | GBA, GC_MAC, GC_SP, GC_LV, GC_BM |
| Biomarkers | GL1, LYSOGL1, CHITR, FERRIT |
| Organ Volumes | SV (spleen), LV (liver) |
| Haematology | HGB (Hb), PLT (platelets) |
| Bone | BMD, OC (osteoclast), OB (osteoblast) |
| Inflammation | IL6 (composite cytokine), NFKB |

#### Treatment Scenarios

| Scenario | Treatment | Reference |
|----------|-----------|-----------|
| **S1** | Natural history | Charrow 2000 Arch Intern Med |
| **S2** | Imiglucerase 60 U/kg IV Q2W | Grabowski 2009 Genet Med |
| **S3** | Velaglucerase α 60 U/kg IV Q2W | Zimran 2010 Blood |
| **S4** | Eliglustat 84 mg BID (CYP2D6 EM) | Mistry 2015 JAMA |
| **S5** | Eliglustat 84 mg QD (CYP2D6 PM) | Balwani 2021 Am J Hematol |
| **S6** | Low-dose ERT 30 U/kg + Eliglustat 84 mg BID | Combination strategy |

#### Key PK Parameters

| Drug | CL | Vd | t½ | Source |
|------|----|----|----|--------|
| Imiglucerase | 1.4 L/h/kg | 0.73 L/kg | 15–45 min | Aerts 2003 |
| Velaglucerase α | ~1.3 L/h/kg | similar | 11–68 min | Zimran 2010 |
| Eliglustat (EM) | 38 L/h | 106 L | ~2h | Lukina 2014 |
| Eliglustat (PM) | 8 L/h | 106 L | ~6–8h | Balwani 2021 |
| Miglustat | 4.5 L/h | 28 L | ~6h | Cox 2000 |

---

### 3. Shiny Dashboard (`gcd_shiny_app.R`)

Nine interactive tabs:

| Tab | Content |
|-----|---------|
| **Overview** | Disease background + mechanistic map thumbnail |
| **Patient Profile** | Body weight, baseline labs, disease type, treatment selection |
| **Drug PK** | ERT plasma curve (first 4 weeks), SRT daily profiles, PK table |
| **Enzyme & Substrate** | GBA activity %, GC burden by compartment, GCS inhibition, NF-κB |
| **Organ / Haematology** | Spleen/liver volume, Hb, platelets with therapeutic goal lines |
| **Bone** | BMD T-score trajectory, osteoblast/osteoclast balance |
| **Scenario Comparison** | All 6 scenarios plotted together for selected endpoints |
| **Biomarkers** | GL-1, lyso-GL1, chitotriosidase, ferritin/IL-6 time courses |
| **Virtual Population** | Monte Carlo (n = 50–500), response rate summary |

---

### 4. References (`gcd_references.md`)

62 PubMed citations across 14 sections:
Pathophysiology · Enzyme Biology · Biomarkers · ERT (imiglucerase) · ERT (velaglucerase / taliglucerase) · SRT (eliglustat) · SRT (miglustat / venglustat) · Pharmacological Chaperone · Bone · Haematology · Neurology & GBA-PD · Clinical Monitoring · Epidemiology · QSP Modelling

---

## Running the Model

```bash
# Render the mechanistic map
dot -Tsvg gcd_qsp_model.dot -o gcd_qsp_model.svg
dot -Tpng -Gdpi=150 gcd_qsp_model.dot -o gcd_qsp_model.png
```

```r
# Install dependencies
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))

# Run mrgsolve model
source("gcd_mrgsolve_model.R")

# Launch Shiny dashboard
install.packages(c("shiny", "shinydashboard", "DT", "plotly"))
shiny::runApp("gcd_shiny_app.R")
```

---

## Key Clinical Thresholds

| Endpoint | Therapeutic Goal (Type 1) | Source |
|----------|--------------------------|--------|
| Hemoglobin | ≥11 g/dL (women), ≥12 g/dL (men) | Pastores 2004 |
| Platelets | >100 ×10⁹/L | Pastores 2004 |
| Spleen volume | <5× normal (≤1500 mL) | Weinreb 2004 |
| Liver volume | ≤1.25× normal | Pastores 2004 |
| GL-1 (plasma) | <1 μg/L | Murugesan 2016 |
| Lyso-GL1 | <1.5 ng/mL | Dekker 2011 |
| Chitotriosidase | <100 nmol/h/mL | Smid 2015 |

---

## GBA–Parkinson Disease Connection

A unique feature of this model is the mechanistic link between GBA deficiency and Parkinson's disease risk:
- **Heterozygous GBA variants** increase PD risk ~5× (Sidransky 2009 NEJM)
- Reduced lysosomal GBA → α-synuclein accumulation → Lewy body formation (Mazzulli 2011 Cell)
- **Bidirectional loop**: α-syn further inhibits GBA activity
- ERT/SRT/chaperones are being explored as disease-modifying agents in GBA-PD

---

*Generated by Claude Code Routine — 2026-06-24*
