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
- `*.dot` — Mechanistic map source (Graphviz DOT format, >100 nodes)
- `*.svg` / `*.png` — Rendered mechanistic map
- `*references.md` — Curated PubMed literature
- `*_model.R` / `*_mrgsolve_model.R` — mrgsolve ODE model + simulation scenarios
- `shiny_app/app.R` / `*_shiny_app.R` — Interactive Shiny PK/PD dashboard

| Date | Disease | Category | Mechanism Summary | Map | Model | Refs | App |
|------|---------|----------|-------------------|-----|-------|------|-----|
| 2026-06-16 | [**Rheumatoid Arthritis**](#rheumatoid-arthritis) | 자가면역질환 | T/B cell–driven synovitis; TNF-α / IL-6 / JAK-STAT; bone erosion (RANKL/OPG); cDMARDs + biologics (TNFi, IL-6Ri, JAKi) | [![RA](rheumatoid-arthritis/ra_qsp.png)](rheumatoid-arthritis/ra_qsp.svg) | [R](rheumatoid-arthritis/ra_model.R) | [refs](rheumatoid-arthritis/references.md) | [Shiny](rheumatoid-arthritis/shiny_app/app.R) |
| 2026-06-16 | [**Pulmonary Arterial Hypertension**](#pulmonary-arterial-hypertension-pah) | 만성질환 / 폐혈관 | EC dysfunction → ET-1↑/NO↓/PGI₂↓ → vasoconstriction + PASMC remodelling; BMPR2 loss; RV failure; ERA + PDE5i + PGI₂ | [![PAH](pulmonary-arterial-hypertension/pah_qsp_model.png)](pulmonary-arterial-hypertension/pah_qsp_model.svg) | [R](pulmonary-arterial-hypertension/pah_mrgsolve_model.R) | [refs](pulmonary-arterial-hypertension/pah_references.md) | [Shiny](pulmonary-arterial-hypertension/pah_shiny_app.R) |
| 2026-06-16 | [**Systemic Lupus Erythematosus**](#systemic-lupus-erythematosus-sle) | 자가면역질환 | NETs/cGAS-STING/TLR7/9 → pDC → IFN-α → BAFF↑ → B cell hyperactivation → anti-dsDNA → IC → complement consumption + lupus nephritis; HCQ + belimumab + anifrolumab + MMF + voclosporin | [![SLE](systemic-lupus-erythematosus/sle_qsp.png)](systemic-lupus-erythematosus/sle_qsp.svg) | [R](systemic-lupus-erythematosus/sle_model.R) | [refs](systemic-lupus-erythematosus/sle_references.md) | [Shiny](systemic-lupus-erythematosus/shiny_app/app.R) |

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

---

### Pulmonary Arterial Hypertension (PAH)

> Directory: [`pulmonary-arterial-hypertension/`](pulmonary-arterial-hypertension/)

> 폐동맥 내피세포 기능 부전 및 혈관 평활근 비대로 폐혈관 저항이 상승하여 우심실 부전을 유발하는 희귀 진행성 질환.

**Mechanistic Map** (130+ nodes, 10 pathway clusters):

| Cluster | Coverage |
|---------|----------|
| ET-1 Vasoconstriction | ET-1 → ETA/ETB → Gq/G12 → PLC/RhoA → IP₃/Ca²⁺/ROCK → MLCK/MLCP → vasoconstriction |
| NO-cGMP Vasodilation | eNOS/BH4/L-Arg → NO → sGC → cGMP → PKG → Kv-channel opening → vasodilation; PDE5 degradation |
| PGI₂-cAMP Vasodilation | AA → COX → PGIS → PGI₂ → IP-receptor → Gs → AC → cAMP → PKA → vasodilation |
| Growth Factor Signalling | BMPR2/ALK1/BMP9 → SMAD1/5/8 → ID1; TGF-β → SMAD2/3; PDGFR/VEGFR/FGFR → RAS-ERK/PI3K-AKT-mTOR/STAT3 |
| Vascular Remodelling | PASMC proliferation, EC dysfunction, EndMT, collagen synthesis, MMP-2/9, ECM remodelling, plexiform lesion |
| Inflammation & Immune | Th1/Th2/Th17/Treg, M1/M2 macrophage, mast cells, B cells, NK/CD8+, IL-6/IL-1β/TNF-α → NF-κB |
| Hypoxia & Metabolism | HIF-1α/2α, VHL/PHD axis, Warburg glycolysis, mitochondrial fission (DRP1/OPA1), Kv1.5↓, ROS/NOX4/SOD2 |
| RV-PA Haemodynamics | PVR → mPAP (Ohm's law); Ees/Ea coupling; Frank-Starling; RV hypertrophy → failure; TR; IVS bowing; RAAS/SNS |
| Biomarkers | BNP/NT-proBNP, hs-TnI, TAPSE, RVEF, RAP, ET-1, DLCO |
| Drug PK/PD | ERA (bosentan/ambrisentan/macitentan, IC₅₀ Emax Hill); PDE5i (sildenafil/tadalafil); sGC stimulator (riociguat); PGI₂ analogues (epoprostenol/treprostinil/selexipag) |

**Mechanistic Map Preview:**

[![PAH QSP Map](pulmonary-arterial-hypertension/pah_qsp_model.png)](pulmonary-arterial-hypertension/pah_qsp_model.svg)

**mrgsolve Model Summary:**

| Compartment | States | Description |
|-------------|--------|-------------|
| ERA PK | ERA_gut, ERA_central, ERA_periph, ERA_effect | 2-compartment oral + effect-site equilibration |
| PDE5i PK | PDE5_gut, PDE5_central, PDE5_effect | 1-compartment oral + effect-site |
| PGI₂ PK | PGI2_central, PGI2_effect | IV infusion (t½ ~3 min) |
| PD — Mediators | ET1, cGMP, cAMP | Turnover ODE for each second messenger |
| PD — Structure | VRI | Vascular Remodelling Index (logistic growth, drug reversal) |
| RV Function | Ees_RV | Adaptive hypertrophy → maladaptive decompensation |
| Biomarker | BNP_conc | Wall-stress–driven BNP production |

**Key ODE relationships:**
- `dVRI/dt = k_growth × VRI × (1 − VRI/VRI_max) − (k_ERA×ERA_Inh + k_PDE×ΔcGMP + k_PGI×PGI₂_Act) × VRI`
- `PVR = PVR_normal + (PVR_PAH0 − PVR_normal) × [tone_frac × ET1/ET1_PAH × (1−ERA_Inh) × cGMP₀/cGMP × cAMP₀/cAMP + remod_frac × VRI/VRI₀]`
- `mPAP = CO × PVR/80 + PAWP`
- Simulates 6 scenarios (no treatment → triple ERA+PDE5i+PGI₂) over 12 weeks

**Files:**

| File | Description |
|------|-------------|
| [`pah_qsp_model.dot`](pulmonary-arterial-hypertension/pah_qsp_model.dot) | Graphviz DOT source (736 lines, 130+ nodes, 8 subgraphs) |
| [`pah_qsp_model.svg`](pulmonary-arterial-hypertension/pah_qsp_model.svg) | Vector mechanistic map (256 KB) |
| [`pah_qsp_model.png`](pulmonary-arterial-hypertension/pah_qsp_model.png) | Raster mechanistic map (7.3 MB, 150 dpi) |
| [`pah_references.md`](pulmonary-arterial-hypertension/pah_references.md) | 40 annotated references with PubMed links |
| [`pah_mrgsolve_model.R`](pulmonary-arterial-hypertension/pah_mrgsolve_model.R) | mrgsolve ODE model + 6 treatment scenarios + dose-response |
| [`pah_shiny_app.R`](pulmonary-arterial-hypertension/pah_shiny_app.R) | 8-tab Shiny dashboard (patient profile, PK, DR, ESC/ERS risk) |

---

### Systemic Lupus Erythematosus (SLE)

> Directory: [`systemic-lupus-erythematosus/`](systemic-lupus-erythematosus/)

> NETs/cGAS-STING/TLR7/9 → pDC → IFN-α 과다 → BAFF↑ → B세포 과활성화 → 자가항체(anti-dsDNA) → 면역복합체 → 보체 소모 + 루푸스 신염을 유발하는 복잡한 다장기 자가면역질환.

**Mechanistic Map** (161 nodes, 9 pathway clusters):

| Cluster | Coverage |
|---------|----------|
| Cell Death & Antigen Release | UV → apoptosis, NETosis/PAD4 citrullination, secondary necrosis, HMGB1, mitochondrial DNA, defective phagocytosis, TREX1 deficiency |
| Innate Immunity & Type I IFN | pDC/mDC/monocytes, TLR7/TLR9/TLR4, cGAS-STING → cGAMP, IRF3/IRF5/IRF7/MyD88/TRIF, NF-κB, IFN-α/β, IFNAR1/2 → JAK1/TYK2 → ISGF3(STAT1/STAT2/IRF9) → ISGs (MX1/IFIT1/OAS/ISG15), DNase I, IL-1β/NLRP3/TNF-α/IL-12/IL-18 |
| B Cell Activation & Autoantibodies | Naïve B → GC → LLPC; BAFF/APRIL → BAFF-R/TACI/BCMA; BCR/SYK/BTK/PI3Kδ; CD40-CD40L; IL-21; anti-dsDNA/anti-Sm/anti-Ro/anti-La/anti-C1q/anti-nucleosome/ANA; IgG immune complexes |
| T Cell Dysregulation | Th1/Th2/Th17/Tfh/Treg/DN T cells/CD8+CTL; T-bet/RORγt/FoxP3; IFN-γ/IL-17A/IL-4/IL-6/IL-2/TGF-β; calcineurin/NFAT; mTOR; JAK3/STAT3/STAT5; CTLA-4/PD-1/PDL-1/ICOS |
| Complement System | Classical (C1q/C1r-C1s/C4/C2/C3), Alternative (Factor B/H/Properdin/C3 tick-over), Lectin (MBL); C3a/C5a anaphylatoxins; MAC (C5b-9); FcγRIIA/III; CR1; C4d deposition; C1q/C2/C4 gene deficiency |
| Lupus Nephritis | IC deposition → mesangial/endothelial/podocyte injury → GBM damage → proteinuria; eGFR decline → ESRD; TGF-β fibrosis; ISN/RPS Class I–VI; VEGF; IL-6; neutrophil infiltration; hematuria; hypertension |
| Multi-Organ Involvement | Malar rash/discoid lupus/SCLE/photosensitivity; BBB disruption → anti-NMDAR → NPSLE; accelerated atherosclerosis; endothelial dysfunction; aPL antibodies → APS/thrombosis; serositis; arthritis; alopecia; cytopenias |
| Drug PK/PD | HCQ (lysosomotropic TLR7/9 block; Vd ~5000 L, t½ ~50 days); Belimumab (anti-BAFF TMDD; t½ ~19 days; FcRn); Anifrolumab (anti-IFNAR1; 300 mg q4w); MMF/MPA (IMPDH inhibition; TMDD not applicable); Voclosporin/Tacrolimus (calcineurin-NFAT); Cyclophosphamide (DNA alkylation); Azathioprine (6-MP/PRPP); Corticosteroids (GR transrepression); Obinutuzumab (anti-CD20 B depletion) |
| Biomarkers & Clinical Endpoints | SLEDAI-2K; BILAG; SLICC damage index; anti-dsDNA titer; C3/C4 serum; ESR/CRP; lymphopenia; thrombocytopenia; UPCR; IFN score; BAFF level; SRI-4; BICLA; complete/partial renal response; LLDAS; DORIS remission; flare rate |

**Mechanistic Map Preview:**

[![SLE QSP Map](systemic-lupus-erythematosus/sle_qsp.png)](systemic-lupus-erythematosus/sle_qsp.svg)

**mrgsolve Model Summary:**

| Compartment | States | Description |
|-------------|--------|-------------|
| HCQ PK | HCQ_gut, HCQ_cent | 1-compartment oral; Vd ~5200 L, t½ ~50 days |
| Belimumab PK+TMDD | BELI_cent, BELI_periph, BAFF_free, BELI_cmplx | 2-compartment IV + TMDD for free BAFF |
| Anifrolumab PK | ANIF_cent | 1-compartment IV (300 mg q4w); IFNAR occupancy |
| MMF/MPA PK | MMF_gut, MPA_cent | 1-compartment oral; active metabolite MPA |
| Voclosporin PK | VOC_gut, VOC_cent | 1-compartment oral; IC50 ~0.8 ng/mL calcineurin |
| IFN pathway PD | TLR_act, IFNa_conc, IFNscore | TLR→pDC→IFN-α→ISG cascade; HCQ/anifrolumab effects |
| B cell axis | Bcell, PBlast, LLPC | BAFF-driven survival; MMF/CS/Cycloph suppression |
| Autoantibodies | AntiDsDNA, IC_burden | LLPC→anti-dsDNA→IC turnover |
| Complement | C3_serum, C4_serum | IC-driven consumption; liver synthesis |
| Kidney | Proteinuria, eGFR_ode | IC→glomerular injury→UPCR; voclosporin podocyte protection |

**Key ODE relationships:**
- `dIFNa/dt = k_pDC × TLR_act × (1 − E_ANIF) − k_IFNA × IFNa` (anifrolumab blocks downstream signaling via IFNAR1 occupancy)
- `dBAFF/dt = Rin_BAFF − kout_BAFF × BAFF − kon × (BELI/V1) × BAFF + koff × complex + k_IFN_BAFF × IFNa` (TMDD + IFN→BAFF upregulation feedback)
- `dProteinuria/dt = k_IC_prot × IC_burden − kout_prot × Prot × (1 + 1.5×E_VOC + 0.5×E_CS)` (voclosporin podocyte stabilization)
- Simulates 7 scenarios (no Tx → HCQ → +MMF → +Belimumab → +Anifrolumab → LN regimen → triple LN) over 52 weeks

**Files:**

| File | Description |
|------|-------------|
| [`sle_qsp.dot`](systemic-lupus-erythematosus/sle_qsp.dot) | Graphviz DOT source (fdp layout, 161 nodes, 9 subgraphs) |
| [`sle_qsp.svg`](systemic-lupus-erythematosus/sle_qsp.svg) | Vector mechanistic map (218 KB) |
| [`sle_qsp.png`](systemic-lupus-erythematosus/sle_qsp.png) | Raster mechanistic map (150 dpi) |
| [`sle_references.md`](systemic-lupus-erythematosus/sle_references.md) | 66 annotated references with PubMed links |
| [`sle_model.R`](systemic-lupus-erythematosus/sle_model.R) | mrgsolve ODE model + 7 treatment scenarios + population simulation (N=200) |
| [`shiny_app/app.R`](systemic-lupus-erythematosus/shiny_app/app.R) | 9-tab Shiny dashboard (patient profile, PK, IFN axis, anti-dsDNA, LN monitor, SLEDAI, population variability, map viewer) |
