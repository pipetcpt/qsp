# qsp

## mrgsolve

- <https://vantage-research.net/qsp-in-r/>
- gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve
    - <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/>
    - <https://github.com/Genentech/gPKPDviz/>
    



## iqrtools

- <https://www.intiquan.com/acop2019_qsp/>

---

## QSP Disease Models

Each model lives in its own subdirectory (lowercase with dashes). Every directory contains:
- `*.dot` — Mechanistic map source (Graphviz DOT format)
- `*.svg` / `*.png` — Rendered mechanistic map
- `references.md` — Key literature (PubMed links)
- `*_model.R` — mrgsolve ODE model (compilable R code)
- `shiny_app/app.R` — Interactive Shiny PK/PD simulator

| Date | Disease | Category | Mechanism Summary | Map | Model | App |
|------|---------|----------|-------------------|-----|-------|-----|
| 2026-06-16 | [Rheumatoid Arthritis](#rheumatoid-arthritis) | 자가면역질환 | T/B cell–driven synovitis; TNF-α / IL-6 / JAK-STAT; bone erosion (RANKL/OPG); cDMARDs + biologics (TNFi, IL-6Ri, JAKi) | [DOT](rheumatoid-arthritis/ra_qsp.dot) · [SVG](rheumatoid-arthritis/ra_qsp.svg) | [R](rheumatoid-arthritis/ra_model.R) | [Shiny](rheumatoid-arthritis/shiny_app/app.R) |

---

### Rheumatoid Arthritis

> Directory: [`rheumatoid-arthritis/`](rheumatoid-arthritis/)

**Mechanistic Map** (130+ nodes, 13 pathway clusters):

| Cluster | Coverage |
|---------|----------|
| Immune Cell Activation | DC, Macrophage M0/M1, Th1/Th17/Treg, Tfh, B cells, Plasma cells, NK, Neutrophils, Mast cells, NETosis |
| Autoantibody Formation | PAD4, citrullination, ACPA/RF, immune complexes, complement (C1q/C3/C5a), FcγRIII |
| Cytokine Network | TNF-α, IL-1β, IL-6, IL-8, IL-10, IL-12, IL-17A/F, IL-21, IL-23, IFN-γ, GM-CSF, TGF-β, VEGF, OSM, CXCL13, CCL2/5 |
| Intracellular Signaling | JAK1/2/3/TYK2, STAT1/3/4/5/6, NF-κB/IKK, p38 MAPK, ERK, JNK, PI3K/Akt/mTOR, AP-1, HIF-1α, NFATc1, Wnt/β-catenin |
| Synovial Pathology | FLS (quiescent/activated), synovial macrophage, pannus, COX-2/PGE2/LTB4, MMPs 1/2/3/9/13, ADAMTS-4/5, TIMPs |
| Synovial Vasculature | Endothelium, ICAM-1/VCAM-1, E-selectin, VEGF-driven neo-angiogenesis, HIF-1α, eNOS |
| Bone & Cartilage | RANKL/RANK/OPG, osteoclast/osteoblast, bone erosion, cartilage ECM (type II collagen, aggrecan), Wnt/DKK-1/sclerostin, chondrocyte apoptosis |
| Pain Signaling | Nociceptor (C/Aδ), NGF/TrkA, TRPV1, Nav1.7/1.8, Substance P/CGRP, bradykinin, peripheral & central sensitization |
| Drug PK | Oral/SC/IV compartments, FcRn recycling, TMDD, ADA, protein binding |
| Drug Mechanisms | MTX (DHFR/adenosine), LEF (DHODH), HCQ (lysosomal), SSZ (NF-κB); 5×TNFi; 4×IL-6Ri; Abatacept; Rituximab; Denosumab; 4×JAKi |
| Biomarkers | CRP, ESR, IL-6, TNF-α, MMP-3, RF, Anti-CCP, fibrinogen, SAA, hepcidin, Hb, RANKL/OPG serum |
| Clinical Endpoints | DAS28-CRP/ESR, ACR20/50/70, CDAI, SDAI, HAQ-DI, Sharp/vdH score, RAMRIS, EULAR remission, MBDA |
| Systemic Complications | CV risk, atherosclerosis, BMD loss, RA-ILD, lymphoma, infections, VTE (JAKi), anemia, depression |

**Mechanistic Map Preview:**

[![RA QSP Map](rheumatoid-arthritis/ra_qsp.png)](rheumatoid-arthritis/ra_qsp.svg)

**mrgsolve Model Summary:**

| Compartment | States | Description |
|-------------|--------|-------------|
| TCZ PK | DEPOT_TCZ, C1_TCZ, C2_TCZ | 2-compartment SC/IV with first-order absorption |
| TMDD | R_FREE, RC | Free sIL-6Rα and TCZ-receptor complex |
| MTX PK | GI_MTX, C1_MTX | 1-compartment oral; active metabolite proxy |
| PD | TNFa, IL6, CRP, RANKL_pd | Cytokine and biomarker ODE dynamics |
| Scores | SJC28_ode, TJC28_ode | Inflammation-driven joint count adaptation |

**Key ODE relationships:**
- `dCRP/dt = ksyn_CRP × IL6_eff/(EC50+IL6_eff) − kout_CRP × CRP` (IL-6R blockade suppresses IL6_eff via TMDD)
- `DAS28-CRP = 0.56√TJC + 0.28√SJC + 0.36·ln(CRP+1) + 0.014·PtGA + 0.96`

**Files:**

| File | Description |
|------|-------------|
| [`ra_qsp.dot`](rheumatoid-arthritis/ra_qsp.dot) | Graphviz source (fdp layout, 130+ nodes) |
| [`ra_qsp.svg`](rheumatoid-arthritis/ra_qsp.svg) | Vector mechanistic map (511 KB) |
| [`ra_qsp.png`](rheumatoid-arthritis/ra_qsp.png) | Raster mechanistic map (3.3 MB, 150 dpi) |
| [`references.md`](rheumatoid-arthritis/references.md) | 66 annotated references (PubMed links) |
| [`ra_model.R`](rheumatoid-arthritis/ra_model.R) | mrgsolve ODE model + 4 dosing scenarios |
| [`shiny_app/app.R`](rheumatoid-arthritis/shiny_app/app.R) | Interactive Shiny simulator |
