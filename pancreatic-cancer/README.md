# Pancreatic Ductal Adenocarcinoma (PDAC) QSP Model

## Disease Overview

Pancreatic ductal adenocarcinoma (PDAC) is the most lethal solid malignancy, with a 5-year survival rate of ~12% (all stages combined). It accounts for ~3% of all cancer diagnoses but ~7% of cancer deaths. Key characteristics:
- **KRAS mutations**: ~95% (G12D 44%, G12V 26%, G12R 14%)
- **TP53 loss**: ~75%
- **CDKN2A loss**: ~90%
- **SMAD4 loss**: ~55%
- **HRD (BRCA1/2/PALB2)**: ~15–20% germline + somatic
- Dense desmoplastic stroma: 80–90% of tumor volume
- Profound immunosuppression (TME)

---

## Mechanistic Map Preview

[![PDAC QSP Mechanistic Map](pdac_qsp_model.png)](pdac_qsp_model.svg)

*Click image to view full SVG*

---

## Deliverables

| File | Description | Format |
|------|-------------|--------|
| `pdac_qsp_model.dot` | Graphviz mechanistic map | DOT |
| `pdac_qsp_model.svg` | Vector mechanistic map | SVG |
| `pdac_qsp_model.png` | Raster mechanistic map (150 dpi) | PNG |
| `pdac_mrgsolve_model.R` | ODE PK/PD model | R/mrgsolve |
| `pdac_shiny_app.R` | Interactive dashboard | R/Shiny |
| `pdac_references.md` | Annotated bibliography (40 refs) | Markdown |

---

## Key Biological Pathways Modeled

| Pathway | Key Nodes | Drug Targets |
|---------|-----------|--------------|
| KRAS/MAPK | KRAS→BRAF→MEK→ERK→MYC | MRTX1133 (G12D), Adagrasib (G12C), SOS1-i |
| PI3K/AKT/mTOR | PI3K→PIP3→AKT→mTORC1→S6K1 | Copanlisib, Everolimus |
| TGF-β/SMAD/EMT | TGFβ→SMAD4→SNAIL→E-Cad loss | Galunisertib (TGFβR1-i) |
| Desmoplastic Stroma | PSC→CAF→Collagen/HA→IFP↑→Drug↓ | PEGPH20, Pirfenidone |
| Immune TME | PD-L1↑, Treg↑, MDSC↑, CD8↓ | Pembrolizumab, Ipilimumab |
| Angiogenesis | HIF1α→VEGF→VEGFR2→Neovascularization | Bevacizumab, Ramucirumab |
| DNA Damage/HRD | BRCA1/2 loss→HRD→PARP dependency | Olaparib, Rucaparib |

---

## Treatment Scenarios & Clinical Trial Reference Data

| Scenario | Regimen | mPFS (mo) | mOS (mo) | ORR | Trial |
|----------|---------|-----------|----------|-----|-------|
| 1 | Untreated Control | 1.5 | 3.0 | 0% | — |
| 2 | Gemcitabine mono | 3.7 | 6.7 | 7% | MPACT control arm |
| 3 | Gem + nab-Paclitaxel | 5.5 | 8.5 | 23% | MPACT (Von Hoff 2013) |
| 4 | FOLFIRINOX | 6.4 | 11.1 | 31.6% | PRODIGE4 (Conroy 2011) |
| 5 | mFOLFIRINOX | 6.0 | 10.5 | 28% | Adapted from PRODIGE24 |
| 6 | MRTX1133 (KRAS G12D) | 4.0 | ~8.0 | ~40% | Phase I/II ongoing |
| 7 | Olaparib (BRCA+) | 7.4 | NS | — | POLO (Golan 2019) |

---

## Model Architecture

### PK Components (mrgsolve)

| Drug | Model | Route | Active Metabolite |
|------|-------|-------|-------------------|
| Gemcitabine | 2-compartment IV | IV infusion | dFdCTP (intracellular) |
| nab-Paclitaxel | 2-compartment IV | IV infusion | Paclitaxel (free) |
| Oxaliplatin | 1-compartment | IV infusion | Platinum-DNA adducts |
| Irinotecan/SN-38 | Sequential 1-cmt | IV infusion | SN-38 → SN-38G |
| 5-FU | 1-compartment | IV bolus/infusion | FdUMP |
| MRTX1133 | 1-compartment oral | PO | Parent compound |
| Olaparib | 1-compartment oral | PO | Parent compound |

### PD Components

- **Simeoni TGI model**: x0→x1→x2→x3→Tumor (transit damage compartments for cytotoxic effect)
- **KRAS signaling**: normalized ERK/RSK activity with Hill-type inhibition by MRTX1133
- **Stroma resistance factor**: modulates effective intratumoral drug penetration (k_pen × [1 − σ × Stroma])
- **CA19-9 biomarker**: linked to tumor burden via first-order kinetics
- **Friberg myelosuppression**: Prol→Tr1→Tr2→Tr3→Circ for neutrophils (ANC)

### ODE System Summary

| Compartment Block | # ODEs | Description |
|-------------------|--------|-------------|
| Drug PK | 14 | Plasma/tissue for 7 drugs |
| Active metabolites | 5 | dFdCTP, Pt-adduct, SN-38, FdUMP, free paclitaxel |
| TGI transit | 4 | x1–x4 damage compartments |
| Tumor volume | 1 | Net growth minus kill |
| KRAS/MAPK signaling | 3 | RAS-GTP, pERK, pRSK |
| Stroma (PSC/CAF) | 3 | Quiescent PSC, activated PSC, CAF pool |
| Immune TME | 4 | CD8+ T cells, Treg, MDSC, NK cells |
| Myelosuppression | 5 | Prol, Tr1, Tr2, Tr3, Circ (ANC) |
| Biomarkers | 2 | CA19-9, CEA |
| **Total** | **41** | |

---

## How to Run

### Prerequisites

```r
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr",
                   "shiny", "bslib", "plotly", "DT"))
```

### Run mrgsolve Simulations

```r
source("pdac_mrgsolve_model.R")
# Produces simulation outputs for all 7 treatment scenarios
```

### Launch Shiny App

```r
shiny::runApp("pdac_shiny_app.R")
```

### Render Mechanistic Map

```bash
dot -Tsvg pdac_qsp_model.dot -o pdac_qsp_model.svg
dot -Tpng -Gdpi=150 pdac_qsp_model.dot -o pdac_qsp_model.png
```

---

## Shiny App Tab Overview

| Tab | Content |
|-----|---------|
| 1. Patient Profile | ECOG status, mutation panel (KRAS/TP53/SMAD4/HRD), stage |
| 2. PK Profiles | Plasma concentration-time curves per drug, Css/AUC summary |
| 3. PD Key Metrics | KRAS inhibition %, stroma permeability, immune score |
| 4. Clinical Endpoints | Tumor volume trajectory, TTP, OS KM curves |
| 5. Scenario Comparison | Side-by-side overlay of all treatment arms |
| 6. Biomarkers | CA19-9 and ANC (neutrophils) time courses |

---

## Clinical Context

PDAC is primarily diagnosed at advanced/metastatic stage (~80% at diagnosis). Key challenges:

1. **Late diagnosis**: CA19-9 not specific for early-stage disease; no approved population screening
2. **Chemoresistance**: Dense desmoplastic stroma (80–90% of tumor volume) reduces drug penetration; KRAS drives intrinsic and acquired resistance
3. **Immune desert phenotype**: Low mutational burden (MSI-H ~1%), high immunosuppressive signaling (PD-L1, TGF-β, IL-10), abundant Tregs and MDSCs
4. **Metabolic reprogramming**: Macropinocytosis and autophagy support KRAS-driven metabolic needs under nutrient-poor conditions
5. **Emerging therapies**: KRAS G12D inhibitors (MRTX1133 Phase I/II), bispecific antibodies (anti-CEA × anti-CD3), mRNA neoantigen vaccines, CAR-T targeting mesothelin/MSLN

---

## References

See [`pdac_references.md`](pdac_references.md) for complete annotated bibliography (40 references).

Key citations:
- Conroy 2011 PRODIGE4 NEJM (PMID: 21561347) — FOLFIRINOX pivotal trial
- Von Hoff 2013 MPACT NEJM (PMID: 24131140) — nab-paclitaxel + gemcitabine
- Golan 2019 POLO NEJM (PMID: 30945893) — olaparib maintenance in BRCA+ PDAC
- Bailey 2016 Nature (PMID: 26909576) — molecular subtypes of PDAC
- Hallin 2020 Cancer Discov (PMID: 31998128) — MRTX849 KRAS G12C inhibitor
- Wang 2022 J Med Chem (PMID: 35505711) — MRTX1133 KRAS G12D inhibitor
