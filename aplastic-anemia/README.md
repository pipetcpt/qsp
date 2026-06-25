# Aplastic Anemia (재생불량성 빈혈) — QSP Model

[![Mechanistic Map](aa_qsp_model.png)](aa_qsp_model.svg)

> **Disease type:** Bone Marrow Failure Syndrome | Immune-Mediated HSC Destruction  
> **Date added:** 2026-06-25  
> **Abbreviation:** AA  

---

## Disease Overview

**Aplastic Anemia (AA)** is a life-threatening bone marrow failure syndrome caused by **T-cell–mediated autoimmune destruction of hematopoietic stem cells (HSCs)**. In ~70% of cases the trigger is idiopathic; recognized causes include viral infections (EBV, CMV, hepatitis, parvovirus B19), chemical exposures (benzene, pesticides), and drugs (chloramphenicol, carbamazepine). Genetic susceptibility via the **HLA-DR15** allele is well established.

### Key Pathophysiological Events
1. Environmental trigger → oligoclonal TCR expansion in autoreactive CD8+ T cells
2. CD8+ CTL activation (Th1 polarization, T-bet up-regulation, IFN-γ/TNF-α excess)
3. IFN-γ → JAK1-STAT1-IRF1 cascade → Fas up-regulation and HSC apoptosis
4. Regulatory T cells (Treg) are simultaneously depleted → loss of immune control
5. BM cellularity falls; adipocyte replacement of marrow (hostile niche)
6. Pancytopenia → infection (ANC-driven), bleeding (PLT-driven), anemia

### Disease Severity Classification (Camitta 1976)

| Category | ANC | Platelets | Reticulocytes |
|---|---|---|---|
| Non-Severe (nsAA) | >0.5×10⁹/L | — | — |
| Severe (sAA) | <0.5×10⁹/L | <20×10⁹/L | <20×10⁹/L |
| Very Severe (vsAA) | **<0.2×10⁹/L** | <20×10⁹/L | <20×10⁹/L |

### Response Criteria

| Response | ANC | Platelets | Hgb |
|---|---|---|---|
| Complete Response (CR) | >1.0×10⁹/L | >100×10⁹/L | >10 g/dL |
| Partial Response (PR) | Transfusion-independent; no longer meets sAA | | |

---

## Mechanistic Map

**File:** `aa_qsp_model.dot` → rendered as `aa_qsp_model.svg` / `aa_qsp_model.png`

### Clusters (12) and Coverage

| # | Cluster | Key Nodes |
|---|---|---|
| 1 | Environmental Triggers & Susceptibility | HLA-DR15, EBV/CMV, benzene, drugs, TCR expansion |
| 2 | Innate Immune Activation | pDC, mDC, NK cells, NKT, TLR3/7/9, NF-κB, IL-12, IL-15, IL-18 |
| 3 | Autoreactive T-Cell Biology | CD4+/CD8+ T cells, Treg (depleted), T-bet, NFATc1, perforin/granzyme, FasL |
| 4 | Cytokine Network | IFN-γ, TNF-α, JAK1-STAT1-IRF1, Caspase-8/3, Fas, EPO, TPO, G-CSF, SCF |
| 5 | BM Microenvironment | MSC, CXCL12/CXCR4, adipocytes, BM cellularity, Notch/Wnt signaling |
| 6 | HSC Compartment | LT-HSC → ST-HSC → MPP → CMP/CLP → GMP/MEP/BFU-E/CFU-E/CFU-Mk |
| 7 | Mature Blood Cells | RBC, Hgb, PLT, ANC, severity classification, CR/PR/NR, PNH clone |
| 8 | Telomere Biology | TERT/TERC/DKC1, telomerase, shelterin, short TL, DDR, HSC senescence |
| 9 | ATG PK | hATG/rATG 2-cpt IV PK, TMDD T-cell depletion, CDC, ADCC, serum sickness |
| 10 | CsA PK/PD | 2-cpt oral PK, CYP3A4, P-gp, CsA–cyclophilin–calcineurin–NFATc1 axis |
| 11 | Eltrombopag & Danazol | EPAG 1-cpt PK, MPL/JAK2-STAT5, HSC & Mk expansion; Danazol–AR–TERT |
| 12 | Supportive Care | G-CSF, transfusions, HSCT, iron chelation, sirolimus, romiplostim |

> **Total: 130+ nodes · 120+ edges**

---

## mrgsolve ODE Model

**File:** `aa_mrgsolve_model.R`

### Compartments (19 ODEs)

| Group | Compartments |
|---|---|
| Immune (3) | CD8 (autoreactive CD8+ T cells), TREG (regulatory T cells), CYTO (cytokine index) |
| HSC & Progenitors (3) | HSC (pool, AU), PROG (myeloid progenitors), CFU_E (erythroid progenitors) |
| Erythroid (2) | RETIC (reticulocytes), RBC (Hgb equivalent, g/dL) |
| Megakaryocyte/PLT (2) | MKP (megakaryocyte progenitors), PLT (platelets, ×10⁹/L) |
| Neutrophil (2) | NEUP (precursors), NEU (mature, ×10⁹/L) |
| Biology (1) | TELO (telomere length, relative) |
| ATG PK (2) | CATG (central), PATG (peripheral) |
| CsA PK (2) | CGUT (gut/oral), CCSA (central) |
| EPAG PK (2) | EGUT (gut/oral), CEPG (central) |

### Treatment Scenarios (6)

| # | Scenario | Basis / Calibration |
|---|---|---|
| 1 | Untreated sAA (natural history) | — |
| 2 | **hATG + CsA** (standard IST) | Young/Scheinberg 2011 NEJM; 6-mo CR ~68% |
| 3 | **rATG + CsA** (inferior IST) | Young/Scheinberg 2011 NEJM; 6-mo CR ~37% |
| 4 | **hATG + CsA + Eltrombopag** | Townsley 2017 NEJM; 6-mo CR ~58% |
| 5 | **Eltrombopag monotherapy** | Desmond 2013 JCI; response ~44% refractory |
| 6 | **CsA monotherapy** | nsAA/elderly; response ~30-40% |

### Key Calibration Targets

| Trial | Drug | Endpoint | Observed | Model Target |
|---|---|---|---|---|
| Young 2011 NEJM (n=120) | hATG+CsA | 6-mo CR | 68% | ~65-70% |
| Young 2011 NEJM (n=120) | rATG+CsA | 6-mo CR | 37% | ~35-40% |
| Townsley 2017 NEJM (n=92) | hATG+CsA+EPAG | 6-mo CR | 58% | ~55-60% |
| Olnes 2012 NEJM (n=25) | EPAG mono | 6-mo resp | 44% | ~40-45% |

---

## Shiny Dashboard

**File:** `aa_shiny_app.R`

### Tabs (6)

| Tab | Content |
|---|---|
| **Patient Profile** | Value boxes (ANC, Hgb, PLT, BM%); diagnosis criteria; 6-mo response gauge |
| **Pharmacokinetics** | ATG, CsA (Cmin), EPAG concentration–time; PK summary table |
| **BM & Immunity** | HSC pool, cytokine index, CD8+ T cells, Treg dynamics |
| **Clinical Endpoints** | ANC, Hgb, PLT over time; response status; clinical data table |
| **Scenario Comparison** | All 6 scenarios side-by-side; 6-mo response summary table |
| **Biomarkers** | Telomere length; CD8:Treg ratio; cytokine–HSC phase plot; risk indices |

### Interactive Controls
- Disease severity selector (sAA / vsAA / nsAA / Normal)
- Treatment checkboxes: hATG, rATG, CsA, EPAG, G-CSF
- Dose sliders: ATG (mg/admin), CsA (mg/day), EPAG (mg/day)
- Simulation duration: 90–730 days

---

## References

**File:** `aa_references.md` — **60 PubMed citations** organized in 13 sections:

1. Disease Overview & Pathophysiology (6)
2. Immunopathogenesis: T-Cell & Cytokine Biology (8)
3. Hematopoiesis & HSC Biology in AA (5)
4. ATG Pharmacology (5)
5. Cyclosporine A Pharmacology (3)
6. Eltrombopag Clinical Trials (5)
7. Eltrombopag Pharmacology (3)
8. Telomere Biology & Danazol (4)
9. PNH & Clonal Evolution (3)
10. Clinical Guidelines & Management (5)
11. Supportive Care & Complications (3)
12. Mathematical & QSP Modeling (6)
13. Recent Clinical Evidence 2020–2025 (4)

---

## Running the Model

```bash
# Render mechanistic map
dot -Tsvg aa_qsp_model.dot -o aa_qsp_model.svg
dot -Tpng -Gdpi=150 aa_qsp_model.dot -o aa_qsp_model.png
```

```r
# Run mrgsolve simulations
install.packages(c("mrgsolve","dplyr","ggplot2","tidyr"))
source("aa_mrgsolve_model.R")
# → outputs 6-scenario results & 4 publication-quality figures

# Launch Shiny dashboard
install.packages(c("shiny","shinydashboard","plotly","DT"))
shiny::runApp("aa_shiny_app.R")
```

---

## Key Biological Insights Captured

- **hATG superiority**: hATG (ATGAM) achieves higher 6-mo CR than rATG (Thymoglobulin) in treatment-naive adults despite comparable T-cell depletion — possibly due to different antigen specificity profiles.
- **EPAG's dual role**: Eltrombopag not only stimulates platelet production via MPL but also directly promotes HSC self-renewal through c-MPL on LT-HSCs, explaining its ability to restore trilineage hematopoiesis.
- **Telomere biology**: Patients with constitutional short telomeres (e.g., TERT/TERC mutations, dyskeratosis congenita) have inferior IST responses and may benefit from Danazol → androgen receptor → TERT transcription upregulation.
- **Treg depletion paradox**: ATG transiently depletes both pathogenic CD8+ T cells and Tregs; subsequent selective Treg recovery drives durable remission — CsA maintains the Treg homeostasis window during ATG-induced lymphopenia.
- **PNH clone evolution**: GPI-AP– cells (PIG-A mutant) that escape immune destruction expand as a "selective advantage" clone and are detectable in ~50% of AA patients; risk of PNH, MDS, AML evolution requires long-term monitoring.
