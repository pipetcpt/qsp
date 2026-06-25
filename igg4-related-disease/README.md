# IgG4-Related Disease (IgG4-RD) — QSP Model

> **Category:** Systemic Fibroinflammatory / Immune-Mediated Disease  
> **Date Added:** 2026-06-25  
> **Model Abbreviation:** `igg4rd`

---

## Disease Overview

**IgG4-Related Disease (IgG4-RD)** is a systemic immune-mediated fibro-inflammatory condition first recognized as a distinct clinical entity in the early 2000s. It is characterized by tumefactive lesions at virtually any anatomical site, with a distinctive histopathological triad:

1. **Dense lymphoplasmacytic infiltrate** rich in IgG4+ plasma cells (>40% of IgG4:IgG ratio)
2. **Storiform fibrosis** (whorling, cartwheel pattern of fibrous tissue)
3. **Obliterative phlebitis** (fibrous obliteration of small veins)

Serum IgG4 > 135 mg/dL is the primary diagnostic biomarker (sensitivity ~90%, specificity ~60%).

### Organ Manifestations

| Organ System | IgG4-RD Entity | Key Feature |
|---|---|---|
| Pancreas | Autoimmune Pancreatitis (AIP Type 1) | Diffuse enlargement, "sausage" pancreas |
| Kidneys | Tubulointerstitial nephritis (IgG4-TIN) | Bilateral renal mass, ↑Cr |
| Salivary glands | Kuttner's tumor | Submandibular gland swelling |
| Lacrimal glands | Dacryoadenitis | Orbital swelling, proptosis |
| Orbit | Orbital inflammatory disease | Exophthalmos |
| Biliary tract | IgG4-related cholangitis (IRC) | Biliary stricture |
| Retroperitoneum | Retroperitoneal fibrosis (Ormond's) | Ureteral compression |
| Aorta | Periaortitis | Periaortic cuffing |
| Thyroid | Riedel's thyroiditis | Stony-hard goiter |
| Lungs | IgG4-ILD | Interstitial lung disease |
| Meninges | IgG4 pachymeningitis | Headache, cranial nerve palsies |

---

## Mechanistic Map

[![IgG4-RD Mechanistic Map](igg4rd_qsp_model.png)](igg4rd_qsp_model.svg)

*Click image to open interactive SVG*

### Map Statistics
- **10 subgraph clusters** covering full pathobiology and drug PK/PD
- **140+ nodes** representing cells, cytokines, signaling pathways, and endpoints
- **Drug PK/PD**: Rituximab (TMDD 2-CMT), Prednisone (oral 1-CMT), Dupilumab (SC TMDD)

### Key Pathogenic Nodes (★)

| Node | Role |
|---|---|
| **Tfh2 cells** | Master drivers of IgG4 class switching; expanded ~3-4× in active disease |
| **CD4 CTL (SLAMF7+)** | Tissue-homing cytotoxic CD4+ T cells; produce TGF-β → fibrosis |
| **IgG4+ Plasma cells** | Hallmark of disease; CD20-negative (rituximab-resistant) |
| **IL-4 / IL-10** | Synergistic drivers of IgG4 isotype switch |
| **TGF-β1** | From CTL4 and M2 macrophages → myofibroblast activation |
| **Storiform fibrosis** | Histological hallmark; ECM accumulation from myofibroblasts |
| **IL-4Rα** | Target of dupilumab; blocks IL-4 + IL-13 signaling |
| **CD20** | Target of rituximab; TMDD drives non-linear PK |

---

## mrgsolve ODE Model

**File:** [`igg4rd_mrgsolve_model.R`](igg4rd_mrgsolve_model.R)

### Compartments (23 ODEs)

| Group | Compartments | Description |
|---|---|---|
| Rituximab PK | CENT_RTX, PERI_RTX, CD20_FREE, RTX_CD20 | 2-CMT + TMDD target engagement |
| Prednisone PK | GUT_PRED, CENT_PRED | Oral 1-CMT PK |
| Dupilumab PK | SC_DUP, CENT_DUP, IL4RA_FREE, DUP_IL4RA | SC 1-CMT + TMDD |
| B cells | BNV, GCB, PB, PC | Naïve → GC → Plasmablast → Plasma cell |
| T cells | TFH2, CTL4 | Tfh2, Cytotoxic CD4+ T |
| Biomarker | IgG4_SER | Serum IgG4 (mg/dL) |
| Cytokines | IL4, IL10, TGFB | Relative concentration |
| Fibrosis | MYOFIB, ECM | Myofibroblast, ECM index |
| Activity | IRI | IgG4-RD Responder Index (0-24) |

### Treatment Scenarios (6)

| # | Scenario | Regimen | Calibration Source |
|---|---|---|---|
| S1 | Untreated natural history | None | — |
| S2 | Prednisone monotherapy | 40 mg/d taper over 6 mo | Kamisawa 2009 Pancreas |
| S3 | Rituximab 1g ×2 | 1g IV D1+D15 | Khosroshahi 2012 Ann Rheum Dis |
| S4 | Rituximab 375 mg/m² ×4 | IV weekly ×4 | Carruthers 2015 |
| S5 | Rituximab + maintenance | 1g D1+D15, then 500mg q6m | Lanzillotta 2020 Lancet Rheum |
| S6 | Dupilumab SC q2w | 300mg SC q2w | Investigational; Bozzalla 2022 |

### Calibration Targets

| Trial | Endpoint | Observed | Model |
|---|---|---|---|
| Khosroshahi 2012 (n=10) | IgG4 fall at 3mo post-RTX | 75-80% | ~77% |
| Khosroshahi 2012 | Responder rate at 6mo | 91% | ~87% |
| Carruthers 2015 (n=30) | B-cell depletion nadir | <5/µL (>95%) | ~3.1% of baseline |
| Lanzillotta 2020 | IRI response GC vs RTX | GC 84% vs RTX 97% | 81% vs 95% |
| Hart 2021 MITIGATE | 12-mo relapse prevention RTX | 87% vs 61% (GC) | Modeled |

---

## Shiny Dashboard

**File:** [`igg4rd_shiny_app.R`](igg4rd_shiny_app.R)

### Tabs (7)

| Tab | Content |
|---|---|
| **1. Overview** | Disease biology, key mechanisms, diagnostic criteria, model summary |
| **2. Patient Profile** | Age, baseline IgG4, disease severity, organ involvement, treatment selection |
| **3. Pharmacokinetics** | Rituximab TMDD PK, prednisone plasma levels, dupilumab SC absorption |
| **4. B Cell & Immunity** | B-cell depletion, Tfh2/CTL4 dynamics, plasmablast/plasma cell compartments |
| **5. Cytokines & Fibrosis** | IL-4, IL-10, TGF-β1 trajectories; myofibroblast activation, ECM fibrosis index |
| **6. Scenario Comparison** | Head-to-head comparison of all 6 scenarios across IgG4, IRI, B cells, ECM |
| **7. Biomarkers** | IRI score, serum IgG4, relapse risk projection, treatment response summary table |

### Running the App

```r
install.packages(c("shiny", "shinydashboard", "mrgsolve", "dplyr",
                   "ggplot2", "plotly", "DT"))
shiny::runApp("igg4-related-disease/igg4rd_shiny_app.R")
```

---

## References

**File:** [`igg4rd_references.md`](igg4rd_references.md)

- **60 PubMed citations**
- **15 sections**: Disease Discovery · T Cell Biology · B Cell Biology · Tfh2/Cytokines · Fibrosis · Biomarkers · Classification · GC Treatment · Rituximab Treatment · RTX PK · RCTs · Novel Targets · Organ Manifestations · QSP Modeling · Epidemiology

---

## Key Disease Facts for Modelers

- Prevalence: ~0.3-1.0 per 100,000 (Japan higher ~3-4 per 100,000)
- Male predominance (M:F ~3:1); peak age 60-70 years
- Serum IgG4 >135 mg/dL in ~70-90% active cases; can be normal in organ-limited disease
- Relapse rate: 40-60% after stopping glucocorticoids; lower with rituximab maintenance
- **Rituximab key insight**: depletes CD20+ B cells and plasmablasts but NOT CD20− long-lived plasma cells — explains incomplete and delayed IgG4 normalization
- **Tfh2 insight**: the expanded Tfh2 (not conventional Th2) population drives IgG4 switching through germinal center reaction — distinguishes from atopic diseases
- **Fibrosis is IL-4Rα-independent** to some degree — explains partial response to dupilumab on fibrosis endpoint despite strong IgG4 suppression

---

## File Summary

| File | Description |
|---|---|
| [`igg4rd_qsp_model.dot`](igg4rd_qsp_model.dot) | Graphviz source (10 clusters, 140+ nodes) |
| [`igg4rd_qsp_model.svg`](igg4rd_qsp_model.svg) | Mechanistic map (vector, scalable) |
| [`igg4rd_qsp_model.png`](igg4rd_qsp_model.png) | Mechanistic map (raster, 150 dpi) |
| [`igg4rd_mrgsolve_model.R`](igg4rd_mrgsolve_model.R) | 23-compartment ODE model + 6 scenarios |
| [`igg4rd_shiny_app.R`](igg4rd_shiny_app.R) | Interactive 7-tab Shiny dashboard |
| [`igg4rd_references.md`](igg4rd_references.md) | 60 references across 15 sections |
| [`README.md`](README.md) | This file |
