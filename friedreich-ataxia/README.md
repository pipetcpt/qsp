# Friedreich Ataxia (FRDA) — QSP Model

Quantitative Systems Pharmacology model for **Friedreich Ataxia (FRDA)** — an autosomal-recessive neurodegenerative disorder caused by GAA-repeat expansion in intron 1 of the *FXN* gene, leading to deficient frataxin, disrupted mitochondrial Fe-S cluster biogenesis, iron accumulation, oxidative stress, and progressive ataxia, cardiomyopathy, and diabetes.

[![FRDA mechanistic map](frda_qsp_model.png)](frda_qsp_model.svg)

## Files

| File | Description |
|------|-------------|
| `frda_qsp_model.dot` | Graphviz source: 12 clusters covering genetics, Fe-S biogenesis, mitochondrial Fe handling, ETC/ATP, NRF2-KEAP1 axis, CNS (DRG / spinocerebellar / dentate), cardiomyopathy, β-cell, skeletal/systemic, omaveloxolone PK/PD, investigational therapies, clinical endpoints (≥160 nodes). |
| `frda_qsp_model.svg` | Vector mechanistic map. |
| `frda_qsp_model.png` | 150-dpi PNG. |
| `frda_mrgsolve_model.R` | mrgsolve model: 15 ODE compartments (FXN, Fe-S, mito Fe, ETC, ATP, ROS, AOX, DRG, cerebellar function, LVMI, β-cell, glucose, AAV-FXN, drug PK depots × 4) + 7 treatment scenarios. |
| `frda_shiny_app.R` | Shiny dashboard, 8 tabs (Patient · FXN/Genetics · Mito-Energetics · CNS · Cardiac · Pancreas · Therapy Comparison · Biomarkers/Safety). |
| `frda_references.md` | 70 PubMed references organized by topic (genetics, Fe-S biology, NRF2, cardiomyopathy, MOXIe trial, gene therapy). |

## Quick start

```bash
# Render the map
dot -Tsvg frda_qsp_model.dot -o frda_qsp_model.svg
dot -Tpng -Gdpi=150 frda_qsp_model.dot -o frda_qsp_model.png
```

```r
# Simulate omaveloxolone 150 mg QD for one year
library(mrgsolve); library(ggplot2)
source("frda_mrgsolve_model.R")
out <- mrgsim(mod, events = ev(amt = 150, cmt = "OMAV_GUT", ii = 24, addl = 364),
              end = 365, delta = 1)
plot(out, FXN + ETC + mFARS + LVMI ~ time)

# Launch the Shiny dashboard
shiny::runApp("frda_shiny_app.R")
```

## Disease biology — capsule

* **Mutation**: GAA·TTC trinucleotide repeat in intron 1 of *FXN* (9q21.11), 66–1700 repeats on both alleles in ~96% of patients; compound heterozygotes carry one point mutation.
* **Molecular**: triplex/R-loop DNA and H3K9me3 heterochromatin silence FXN, reducing frataxin to ~5–35% of wild type.
* **Biochemistry**: frataxin is an allosteric activator of NFS1-ISCU-ISD11-ACP, the mitochondrial Fe-S cluster scaffold; its loss cripples Complexes I/II/III and aconitase, dropping ATP and raising mitochondrial labile Fe and ROS.
* **NRF2 axis**: chronic ROS oxidizes KEAP1 Cys, but FRDA cells show paradoxically impaired NRF2 nuclear translocation — the rationale for omaveloxolone.
* **Clinical**: childhood-onset gait/limb ataxia, dysarthria, areflexia, loss of vibration & proprioception (DRG/dorsal-column path), 60–80% hypertrophic cardiomyopathy (HFpEF → HFrEF), ~30% diabetes, scoliosis, optic/auditory neuropathy.
* **Progression**: loss of ambulation typically 10–15 years from onset; leading cause of mortality is cardiac.

## Therapeutic landscape modelled

| Class | Agent | Mechanism captured |
|-------|-------|--------------------|
| **NRF2 activator** | Omaveloxolone (Skyclarys, 150 mg QD) | Cys151-KEAP1 modification, NRF2 stabilization → ARE-genes |
| **Ubiquinone analog** | Idebenone (450 mg TID) | Bypass Complex I deficiency |
| **Fe chelator** | Deferiprone (25 mg/kg/d) | Reduce mitochondrial labile Fe |
| **Protein replacement** | Nomlabofusp / CTI-1601 (SC) | TAT-FXN restoration of frataxin |
| **Gene therapy** | AAVrh10-FXN / LX2006 (single IV) | Durable FXN restoration |
| **Cardiac adjunct** | ACE-inhibitor | LV mass regression |

## Calibration anchors

* **MOXIe Part 2 (NCT02255435)**: mFARS Δ = -2.41 (placebo-corrected) at 48 w; delayed-start extension supports disease-modification (Lynch 2021/2022).
* **Boddaert 2007**: deferiprone 20 mg/kg/d × 6 mo → ↓cardiac Fe (T2*).
* **Pousset 2015**: cardiac progression to LVH within ~10 years of onset.
* **Cnop 2012**: β-cell oxidative-stress-driven apoptosis underlies FRDA diabetes.

## Clinical endpoints in the simulator

* **mFARS** (modified Friedreich Ataxia Rating Scale, 0–93)
* **T25FW** (timed 25-foot walk, seconds)
* **LVMI** (left-ventricular mass index, g/m²)
* **HbA1c**, plasma glucose
* **Plasma frataxin** (proxy biomarker)
* **ALT** (omaveloxolone class effect)
