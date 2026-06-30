# Hodgkin Lymphoma (HL) — QSP Model

A quantitative systems pharmacology (QSP) model of classical and nodular
lymphocyte-predominant Hodgkin lymphoma, covering Reed-Sternberg (HRS) cell
biology, the tumor microenvironment, the PD-1/PD-L1 / 9p24.1 checkpoint axis,
modern frontline regimens (ABVD, AVD, escBEACOPP, A+AVD, N-AVD) and salvage
strategies (brentuximab vedotin and PD-1 blockade monotherapy or combinations).

## Files

| File | Description |
|------|-------------|
| `hl_qsp_model.dot` | Mechanistic map (160+ nodes across 12 subgraph clusters) |
| `hl_qsp_model.svg` | Vector render of the mechanistic map |
| `hl_qsp_model.png` | Raster render (150 dpi) of the mechanistic map |
| `hl_mrgsolve_model.R` | 24-compartment ODE model (`mrgsolve`) with 7 regimens encoded as event tables |
| `hl_shiny_app.R` | 8-tab Shiny dashboard (patient · PK · tumor/TARC · immune · hematology · toxicity · endpoints · biomarkers) |
| `hl_references.md` | 67 PubMed-linked references across 15 thematic sections |
| `README.md` | This document |

## Biological scope

The mechanistic map is organized into twelve subgraph clusters:

1. **Etiology / predisposition** — EBV (LMP1, LMP2A), HLA class II alleles,
   bimodal age, family history (e.g., KLHDC8B), immunosuppression, tobacco.
2. **Cell of origin → HRS** — germinal-center B-cell origin, crippled BCR,
   apoptosis escape (BCL2, c-FLIP), loss of B-cell program (PAX5 dim, OCT2↓,
   BOB1↓), CD30+ CD15+ HRS phenotype, NLPHL popcorn cell.
3. **Oncogenic signaling** — canonical & non-canonical NF-κB, TNFAIP3 (A20)
   loss, JAK2 amplification, STAT3/5/6, PI3K/AKT/mTOR, NOTCH1, AP-1 / JUNB,
   SOCS1 / PTPN1 loss, CIITA rearrangements.
4. **9p24.1 / checkpoint** — PDCD1LG1/2 + JAK2 amplification, PD-L1/PD-L2
   over-expression, PD-1+ TIM-3+ LAG-3+ exhausted T cells, MHC I/II
   modulation, FOXP3+ Treg pool.
5. **Tumor microenvironment** — Th2 / Th17 / TFH cells, eosinophils, plasma
   cells, mast cells (CD30L+), TAM (CD68/CD163), fibroblasts driving NS
   bands, neutrophils, NK cells, effector CD8 T cells.
6. **Cytokine / chemokine network** — IL-6, IL-10, IL-13 (autocrine STAT6),
   TGF-β, TNF-α, IFN-γ, CCL17 (TARC), CCL22 (MDC), CCL5, CXCL10, galectin-1,
   soluble IL-2R (sCD25), soluble CD30, ferritin.
7. **Clinical phenotype / staging** — cervical / mediastinal / axillary
   nodes, bone marrow / splenic involvement, B-symptoms, pruritus,
   alcohol-induced pain, Lugano stage, IPS, GHSG risk groups, histology
   (NS / MC / LR / LD / NLPHL).
8. **Imaging / biomarkers** — FDG-PET SUVmax, Deauville 1-5, iPET2, MTV/TLG,
   ctDNA, serum TARC, ESR/CRP, albumin, Hb, LDH, WBC/lymphopenia.
9. **Drug PK / PD** — PK compartments for Doxo, Bleo, Vinblastine,
   Dacarbazine, Etoposide, Cyclophosphamide, Vincristine, Prednisone,
   Procarbazine, Brentuximab vedotin + free MMAE payload, Nivolumab,
   Pembrolizumab, Bendamustine, Gemcitabine, G-CSF.
10. **Acute / late toxicity** — febrile neutropenia, mucositis, nausea /
    emesis, alopecia, peripheral neuropathy (MMAE/VBL/VCR), anthracycline
    cardiotoxicity, bleomycin pulmonary toxicity, gonadotoxicity / infertility,
    secondary malignancies, immune-related AEs, radiation late effects.
11. **Endpoints / outcomes** — end-of-treatment PET-CR, iPET2 negativity,
    PFS, OS, FFTF, DOR, ORR, QoL / FACT-Lym, HSCT eligibility.
12. **Therapy regimens** — ABVD, AVD, escBEACOPP, A+AVD, N-AVD, K-AVD,
    ISRT-RT, HDT+ASCT, BV mono, nivolumab mono, pembrolizumab mono,
    BV + nivolumab, allogeneic HSCT.

## mrgsolve ODE model

The R model (`hl_mrgsolve_model.R`) implements 24 ODE compartments:

* **PK** — 2-compartment doxorubicin; 1-compartment bleomycin,
  vinblastine, dacarbazine; brentuximab vedotin with explicit
  protease-mediated MMAE release into a free-payload compartment;
  nivolumab linear PK.
* **Tumor** — HRS mass with Simeoni-like growth (`lambda0` / `lambda1`),
  Emax-style kill from each cytotoxic drug, ADC-specific kill from BV /
  MMAE, and immune-mediated kill scaled by the effector CD8 T-cell pool.
* **Immune** — Effector vs exhausted CD8 T-cell compartments coupled to
  PD-L1 expression and the fractional PD-1 occupancy
  `Cnivo / (Cnivo + EC50)`; reinvigoration rate `k_reinvig` activates
  only when PD-1 is occupied. A Treg pool is driven by HRS antigen load.
* **Biomarkers** — Serum TARC, IL-6, IL-13 with production scaled to HRS
  mass; the IL-6/IL-13 surrogate also drives PD-L1 expression (an IFN-γ
  reflex pathway).
* **Hematology** — 5-compartment Friberg-Karlsson myelosuppression chain
  (PROL → TR1 → TR2 → TR3 → ANC) with drug-specific slope terms and an
  optional G-CSF factor.
* **Late toxicity** — Accumulators for doxorubicin cumulative dose,
  bleomycin cumulative exposure, MMAE / vinca neurotoxicity, and
  nivolumab AUC-driven irAE risk, each with a threshold or rate
  parameter calibrated against published toxicity literature.
* **Endpoints** — Continuous Deauville surrogate (function of HRS mass),
  pCR flag, PFS hazard with `h0_PFS + β_tum·HRS + β_TARC·TARC + β_imm·T_eff`.

Pre-built regimen event tables: `ev_ABVD()`, `ev_AVD()`, `ev_BV_AVD()`,
`ev_N_AVD()`, `ev_escBEACOPP()`, `ev_BV_mono()`, `ev_NIVO_mono()`,
`ev_BV_NIVO()`. Use `run_scenario("BV_AVD", tend_days = 365)` for a quick
1-year simulation.

## Shiny dashboard

`hl_shiny_app.R` exposes the model through eight tabs:

1. **Patient profile** — histology, Lugano stage, IPS, baseline MTV, EBV
   positivity, B-symptoms (radar-like risk summary).
2. **Drug PK** — concentration-time curves for the chosen regimen.
3. **Tumor & TARC** — HRS mass / MTV surrogate plus serum TARC.
4. **Immune dynamics** — effector vs exhausted CD8, Treg pool, PD-L1.
5. **Hematology** — Friberg ANC kinetics with grade-3/4 thresholds.
6. **Toxicity** — cumulative cardio, pulmonary, neuropathy, irAE indices.
7. **Endpoints** — surrogate Deauville, PFS hazard, 1-year endpoint table.
8. **Biomarkers** — TARC, IL-6, IL-13 panels.

## Calibration anchors

Parameter values were tuned to reproduce published clinical anchors:

| Anchor | Source | Model output |
|--------|--------|--------------|
| ABVD ×6 → 5-yr PFS ≈ 83% (advanced) | RATHL, NEJM 2016 | `PFS_surv(365 d)` ≈ 0.93 in baseline run |
| escBEACOPP ×6 → 5-yr PFS ≈ 88% | HD15, Lancet 2012 | deeper HRS depletion + higher cumulative toxicity |
| A+AVD vs ABVD: 6-yr PFS 82.3% vs 74.5% | ECHELON-1, NEJM 2018/2022 | BV/MMAE-driven additional HRS kill |
| N-AVD: 2-yr PFS 92% vs A+AVD 83% | SWOG S1826, NEJM 2024 | PD-1 occupancy + reinvigoration adds immune kill |
| Nivolumab R/R ORR ≈ 69%, CR ≈ 16% | CheckMate 205 | reinvigoration alone cannot fully clear bulky disease |
| Doxorubicin cardiotox threshold ≈ 400 mg/m² | Swain JCO 2003 | piecewise `cardio_drive` activates beyond threshold |

## How to run

```r
# Render the mechanistic map (Graphviz required)
system("dot -Tsvg hl_qsp_model.dot -o hl_qsp_model.svg")
system("dot -Tpng -Gdpi=150 hl_qsp_model.dot -o hl_qsp_model.png")

# Simulate a regimen
source("hl_mrgsolve_model.R")
sim <- run_scenario("BV_AVD", tend_days = 365)
head(sim)

# Launch the dashboard
shiny::runApp("hl_shiny_app.R")
```

## References

See `hl_references.md` (67 PubMed-linked entries across 15 thematic sections:
epidemiology, HRS biology, signaling, PD-1 axis, TME, cytokines, staging,
ABVD/AVD, escBEACOPP, BV, PD-1 inhibitors, salvage/HSCT, PK models, late
toxicity, QSP modeling).
