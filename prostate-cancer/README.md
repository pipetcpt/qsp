# Prostate Cancer — QSP Model

> **Disease**: Prostate Cancer (Localised → mHSPC → CRPC → mCRPC)
> **Category**: Oncology / Genitourinary Cancer
> **Date Added**: 2026-06-23
> **Abbreviation**: PC

---

## Disease Overview

Prostate cancer (PCa) is the most common non-skin malignancy and the second leading cause of cancer-related death in men in the United States (~300,000 new cases/year globally). The androgen signaling axis—through the hypothalamic-pituitary-gonadal (HPG) axis—drives the majority of PCa growth.

**Key pathological sequence:**
1. **Localised**: Organ-confined, driven by androgen receptor (AR) signaling
2. **Locally advanced / metastatic hormone-sensitive (mHSPC)**: Spread beyond prostate; still testosterone-sensitive
3. **Castration-Resistant (CRPC)**: PSA rises despite castrate testosterone (<50 ng/dL); AR reactivated by amplification, mutation, or splice variants
4. **Metastatic CRPC (mCRPC)**: Bone and soft tissue metastases; AR-independent pathways activated

---

## Mechanistic Map

[![Prostate Cancer QSP Mechanistic Map](pc_qsp_model.png)](pc_qsp_model.svg)

*Click image to view full SVG. Rendered by `fdp` (Graphviz).*

### Clusters (10 subgraphs, 130+ nodes)

| Cluster | Key Components |
|---------|---------------|
| **HPG Axis** | KiSS1/GPR54 → GnRH → LH/FSH → Testosterone production |
| **Androgen Biosynthesis** | Cholesterol → CYP11A1 → Pregnenolone → CYP17A1 → DHT |
| **AR Signaling** | DHT-AR binding → Nuclear translocation → PSA, TMPRSS2-ERG |
| **Tumor Cell Biology** | CDK4/6, Cyclin D1, pRb/E2F, BCL-2/BAX, p53/MDM2, PARP |
| **PI3K/AKT/mTOR** | PTEN loss → PIP3 → AKT → mTORC1/2 → S6K/4E-BP1 |
| **RAS/MAPK** | EGFR/HER2 → GRB2/SOS → RAS → RAF → MEK → ERK |
| **Bone Metastasis** | RANKL/RANK/OPG, DKK1/Wnt, PTHrP, Endothelin-1, BMD |
| **Drug PK/PD** | GnRH agents, ARPIs, Abiraterone, Docetaxel, PARP-i, Denosumab, Lu-PSMA |
| **Immune TME** | CTL/Treg, PD-1/PD-L1, CTLA-4, TAM, MDSC |
| **Clinical Endpoints** | PSA, Testosterone, Bone Scan Index, rPFS, OS, ARv7 |

---

## mrgsolve ODE Model

**File**: `pc_mrgsolve_model.R`

### Compartments (33 ODEs)

| Group | Compartments |
|-------|-------------|
| HPG Axis | LH, Testosterone, DHT |
| AR Signaling | AR_free, AR_DHT, AR_nuc, PSA |
| Tumor Kinetics | TC_proliferating, TC_quiescent, CRPC_fraction, ARv7_fraction |
| PI3K/AKT | AKT_activity |
| Bone | Osteoclasts, Osteoblasts, BMD, BoneMetastasis |
| Leuprolide PK | Depot, Plasma, Flare effect |
| Degarelix PK | SC depot, Plasma |
| Enzalutamide PK | Gut, Plasma (μM) |
| Abiraterone PK | Gut, Plasma (μM) |
| Docetaxel PK | Central, Peripheral (2-compartment) |
| Olaparib PK | Gut, Plasma |
| Denosumab PK | SC depot, Plasma |

### Treatment Scenarios (7)

| # | Scenario | Drugs |
|---|----------|-------|
| 1 | Untreated (natural history) | — |
| 2 | ADT Alone | Leuprolide 7.5 mg IM q28d |
| 3 | ADT + Enzalutamide | Leuprolide + Enzalutamide 160 mg QD |
| 4 | ADT + Abiraterone | Leuprolide + Abiraterone 1000 mg QD |
| 5 | ADT + Docetaxel ×6 | Leuprolide + Docetaxel 75 mg/m² q3w |
| 6 | Olaparib (HRR-def mCRPC) | Leuprolide + Olaparib 300 mg BID |
| 7 | Sequential ADT→ARPI→Docetaxel | Staged regimen |

### Key Clinical Calibration

| Target | Value | Reference |
|--------|-------|-----------|
| Castrate testosterone | <50 ng/dL within 4 weeks | Standard ADT criteria |
| PSA ≥50% decline (ADT) | ~90% response rate | STAMPEDE, CHAARTED |
| CRPC emergence (ADT alone) | ~12-24 months | Multiple RCTs |
| Enzalutamide OS benefit (mCRPC) | ~4 months | AFFIRM trial |
| Abiraterone PSA50 response | ~29% post-chemo | COU-AA-301 |
| Docetaxel OS benefit | +3.0 months | TAX 327 |
| Olaparib ORR (BRCA2-mut) | 33% | PROfound trial |
| Radium-223 OS benefit | +3.6 months | ALSYMPCA trial |

---

## Shiny Dashboard

**File**: `pc_shiny_app.R`

### Tabs (8)

| Tab | Content |
|-----|---------|
| **Patient Profile** | Demographics, baseline PSA/TC, disease stage, PTEN/HRR status |
| **Drug PK** | PK time-profiles for all 7 drug classes; hormone cascade |
| **HPG Axis & Hormones** | LH, Testosterone (ng/dL), DHT; castrate threshold |
| **AR Signaling & PSA** | Nuclear AR activity, PSA (linear + log), PSA50 waterfall |
| **Tumor Kinetics & Resistance** | TC burden, CRPC fraction, ARv7 emergence, AKT activity |
| **Bone Metastasis** | BMD, bone metastasis burden, OC/OB balance, bone markers |
| **Scenario Comparison** | Multi-regimen comparison with DT endpoints table |
| **Sensitivity Analysis** | 2D sensitivity + tornado plot (PTEN loss, HRR, ARv7) |

### Launch

```r
library(shiny)
library(mrgsolve)
shiny::runApp("prostate-cancer/pc_shiny_app.R")
```

---

## References

**File**: `pc_references.md`

63 curated references organized into 18 sections:

1. Epidemiology
2. AR Biology
3. HPG Axis / Testosterone
4. PI3K/AKT/mTOR
5. RAS/MAPK
6. Bone Metastasis
7. GnRH Agonists/Antagonists
8. ARPIs (Enzalutamide, Abiraterone, Apalutamide, Darolutamide)
9. Chemotherapy (Docetaxel, Cabazitaxel)
10. PARP Inhibitors (Olaparib, Rucaparib)
11. Bone-Targeted Agents
12. PSMA Imaging/Therapy
13. Immunotherapy
14. Resistance Mechanisms
15. QSP Modeling
16. Clinical Staging
17. Pharmacokinetics
18. Immune Microenvironment

---

## Key Biology — Disease Highlights

### Androgen Signaling Axis
The HPG axis maintains testosterone at ~432 ng/dL (15 nmol/L) in healthy men. In the prostate, testosterone is converted to the more potent DHT by 5α-reductase type 2 (SRD5A2). DHT binds AR with ~10× higher affinity (Kd ~1 nM) than testosterone. The AR-DHT complex translocates to the nucleus, dimerizes, and activates androgen response elements (AREs) driving PSA, TMPRSS2, NKX3-1, and pro-proliferative genes.

### Castration-Resistant Mechanisms
- **AR amplification**: Overexpression allows signaling at castrate T levels
- **AR mutation** (T878A, F876L): Converts antagonists to agonists
- **AR splice variants** (ARv7): Lack ligand-binding domain → ligand-independent activity; predicts enzalutamide/abiraterone resistance
- **PI3K/AKT activation** (PTEN loss ~70%): Bypasses AR, provides survival signals
- **Intratumoral androgen synthesis**: CYP17A1 activity in tumor cells

### Bone Metastasis Vicious Cycle
Prostate cancer preferentially metastasizes to bone (>90% of mCRPC patients). The vicious cycle involves:
1. Tumor cells secrete PTHrP → upregulates RANKL on osteoblasts → activates osteoclasts
2. Osteoclasts resorb bone → release stored TGF-β, IGF-1 → feeds back to tumor
3. Predominantly sclerotic lesions (~80%) driven by ET-1 and Wnt-mediated osteoblast stimulation

---

## File Structure

```
prostate-cancer/
├── README.md                  # This file
├── pc_qsp_model.dot           # Graphviz mechanistic map (130+ nodes, 10 clusters)
├── pc_qsp_model.svg           # Rendered SVG
├── pc_qsp_model.png           # Rendered PNG (150 dpi)
├── pc_mrgsolve_model.R        # ODE model (33 compartments, 7 scenarios)
├── pc_shiny_app.R             # Shiny dashboard (8 tabs)
└── pc_references.md           # 63 PubMed references (18 sections)
```

---

## Running the Model

```r
# 1. Render mechanistic map
system("fdp -Tsvg prostate-cancer/pc_qsp_model.dot -o pc_qsp_model.svg")
system("fdp -Tpng -Gdpi=150 prostate-cancer/pc_qsp_model.dot -o pc_qsp_model.png")

# 2. Run mrgsolve model
source("prostate-cancer/pc_mrgsolve_model.R")

# 3. Launch Shiny dashboard
shiny::runApp("prostate-cancer/pc_shiny_app.R")
```

**Requirements**: R ≥ 4.2, mrgsolve, shiny, dplyr, ggplot2, tidyr, DT, shinythemes, Graphviz
