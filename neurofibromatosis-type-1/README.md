# Neurofibromatosis Type 1 (NF1, 신경섬유종증 1형) — QSP Disease Model

> **Neurofibromatosis type 1 (NF1)** is an autosomal-dominant RASopathy
> (birth incidence ~1/2,700, ~50% de novo) caused by **loss-of-function
> mutations in NF1** (17q11.2), which encodes **neurofibromin**, a RAS-GAP
> that normally hydrolyzes RAS-GTP to RAS-GDP. Neurofibromin loss produces
> **constitutive RAS-GTP** and hyperactivation of the downstream
> **RAF-MEK-ERK (MAPK)** cascade in neural-crest-derived Schwann cells,
> driving **plexiform neurofibroma (PN)** growth, near-universal
> **cutaneous neurofibromas**, **optic pathway glioma (OPG)**, and a
> lifetime **8-13% risk of malignant peripheral nerve sheath tumor
> (MPNST)**, alongside skeletal (tibial pseudarthrosis, scoliosis),
> vascular (stenosis, hypertension), and neurocognitive manifestations.

The model captures: **NF1 biallelic loss** (germline + somatic 2nd hit)
→ **neurofibromin RAS-GAP deficiency** → **RAS-GTP accumulation** →
**RAF-MEK-ERK hyperactivation** → NF1-null Schwann-cell proliferation →
**plexiform/cutaneous neurofibroma growth** (puberty/pregnancy-accelerated)
and **optic pathway glioma**, with **atypical-neurofibroma → MPNST**
malignant progression risk. **Selumetinib** and **mirdametinib** are oral
**MEK1/2 inhibitors** that block the shared downstream node, suppressing
pERK and driving tumor-growth-inhibition of PN/OPG volume (REiNS ≥20%
volumetric response), subject to **adaptive RTK-feedback resistance** on
chronic dosing and **regrowth on discontinuation**.

---

## Deliverables

| File | Purpose |
|------|---------|
| `nf1_qsp_model.dot` / `.svg` / `.png` | Graphviz mechanistic map — **16 clusters, 108 nodes** |
| `nf1_mrgsolve_model.R` | **21-ODE** mrgsolve model (6 PK + 15 disease/PD/clinical) with 10 scenarios |
| `nf1_shiny_app.R` | Interactive **8-tab** Shiny dashboard |
| `nf1_references.md` | **54** PubMed/PMC citations grouped by section |

---

## Mechanistic Map — Cluster Index

1. **Genetics & etiology** — NF1 17q11.2 LOF, de novo (~50%), autosomal dominant, somatic 2nd-hit, microdeletion (SUZ12 co-deletion → severe phenotype), mosaicism, Legius syndrome (SPRED1) differential
2. **Neurofibromin — RAS-GAP function** — GRD domain, RAS-GTP↔GDP cycling, cAMP-PKA co-regulation
3. **RAS-RAF-MEK-ERK (MAPK) cascade** — RTK-GRB2-SOS, RAF dimerization, MEK1/2, ERK1/2, PI3K-AKT-mTOR parallel effector, ERK→RTK negative feedback & MEKi-induced adaptive resistance
4. **Schwann-cell / tumor microenvironment biology** — NF1-/- Schwann cell, mast cell/KIT recruitment, fibroblast/ECM remodeling, TAM
5. **Plexiform neurofibroma (PN) pathology** — congenital initiation, puberty/pregnancy-accelerated growth, internal/deep vs. superficial PN, nerve/airway/vascular compression, inoperability
6. **Cutaneous & subcutaneous neurofibroma** — puberty onset, hormone-receptor-driven accumulation, disfigurement
7. **Malignant transformation (MPNST)** — atypical neurofibroma (ANNUBP, CDKN2A), TP53/PRC2-SUZ12/EED loss, metastasis, FDG-PET surveillance
8. **Optic pathway & CNS glioma** — OPG onset (~20% of children), visual-pathway compression, other low-grade glioma, T2-hyperintensities (UBOs)
9. **Skeletal / orthopedic** — sphenoid wing dysplasia, dystrophic scoliosis, tibial dysplasia/pseudarthrosis, short stature, osteopenia
10. **Vasculopathy / cardiovascular** — arterial stenosis, aneurysm/moyamoya-like disease, hypertension, pheochromocytoma, congenital heart disease
11. **CNS / neurocognitive** — RAS-MAPK-dependent synaptic plasticity deficit, learning disability (30-70%), ADHD-like inattention, visuospatial deficits, headache
12. **Dermatologic & ophthalmic hallmarks** — café-au-lait macules, axillary/inguinal freckling, Lisch nodules, choroidal abnormalities
13. **Drug PK — MEK1/2 inhibitors** — selumetinib (25 mg/m2 BID), mirdametinib (2 mg/m2 BID, 3wk-on/1wk-off), trametinib/binimetinib (off-label), dose interruption/reduction
14. **Drug PD — tumor response** — MEK inhibition → pERK suppression → Simeoni-style tumor-growth-inhibition → PN/OPG regression, modest cNF response, adaptive resistance/rebound
15. **Safety / adverse events** — acneiform rash, paronychia, asymptomatic CPK elevation, GI AEs, reversible LVEF decline, retinal pigment epitheliopathy, theoretical growth-plate impact
16. **Clinical endpoints & outcomes** — tumor pain (NRS-11), pain interference, HRQoL, motor function, visual acuity, objective response rate (ORR, REiNS), MPNST-related mortality

---

## mrgsolve Model

### ODE Compartments (21)
**PK (6):** SEL_GUT, SEL_CENT (selumetinib); MIR_GUT, MIR_CENT (mirdametinib) — plus 2 derived plasma-concentration variables (SEL_CP, MIR_CP)

**Disease / PD / clinical (15):** PERK, RESIST, PN_PROLIF, PN_T1, PN_T2, PN_T3
(→ PN_TOTAL), OPG_VOL, CNF_BURDEN, PAIN, QOL, VISION, LVEF, DERM_AE, CPK_AE,
GROWTHZ

### Treatment Scenarios (10)
1. **Untreated (natural history, pediatric)** — logistic PN growth reference, no drug
2. **Selumetinib 25 mg/m2 BID (SPRINT, pediatric)** — approved regimen, target ≥20% REiNS response
3. **Selumetinib 20 mg/m2 BID (dose-reduced, AE)** — dermatologic-AE-driven dose reduction
4. **Mirdametinib 2 mg/m2 BID, 3wk-on/1wk-off (ReNeu, pediatric)** — approved pediatric regimen
5. **Mirdametinib 2 mg/m2 BID, 3wk-on/1wk-off (ReNeu, adult)** — approved adult regimen
6. **Selumetinib, drug holiday then rechallenge** — discontinuation-induced regrowth & re-response
7. **Selumetinib, OPG subgroup (pediatric)** — concurrent optic pathway glioma activated
8. **Selumetinib, poor adherence (60%)** — real-world adherence sensitivity
9. **Trametinib off-label approximation** — re-parameterized MEKi PK/PD block (illustrative)
10. **Mirdametinib adult, long-term (5 yr)** — chronic exposure, cutaneous-NF & safety accrual

### Calibration Anchors
| Trial / source | Endpoint | Modeled target |
|---|---|---|
| Dombi 2016 NEJM (PMID 28029918, phase 1) | 20-30 mg/m2 BID selumetinib, median PN volume Δ | -31% median best response |
| Gross 2020 NEJM (PMID 32187457, SPRINT stratum 1) | Selumetinib 25 mg/m2 BID pediatric ORR | 68% confirmed partial response |
| Moertel 2025 JCO (PMID 39514826, ReNeu) | Mirdametinib BICR-confirmed ORR | 52% pediatric / 41% adult |
| Dombi 2013 Neurology (PMID 24249804, REiNS) | Volumetric response threshold | ≥20% decrease = response |
| Dagalakis 2013 J Pediatr (PMID 24321536) | Puberty-accelerated PN growth | growth-rate multiplier ~2x |
| Cannon 2018 Orphanet J Rare Dis (PMID 29415745) | Cutaneous NF natural-history growth | monthly volumetric growth rate |
| Patel 2017 CPT:PSP (PMID 28326681) | Selumetinib population PK | t½ ~9 h central-compartment PK |

---

## Shiny App — 8 Tabs

1. **Patient & Overview** — covariate sidebar (BSA, baseline PN volume, puberty/OPG toggles) + mechanistic-map schematic
2. **Drug PK** — selumetinib and mirdametinib concentration-time
3. **Pathway PD** — MEK-inhibition fraction, pERK suppression, adaptive resistance
4. **Clinical endpoints** — PN/OPG/cNF volume, pain (NRS-11), HRQoL, visual acuity
5. **Scenario comparison** — all regimens overlaid + endpoint table
6. **Biomarkers** — adaptive resistance and pediatric growth Z-score trajectories
7. **Safety** — LVEF, dermatologic AE composite, CPK elevation
8. **References** — key trial citations

---

## Usage

```bash
# Render the mechanistic map
dot -Tsvg nf1_qsp_model.dot -o nf1_qsp_model.svg
dot -Tpng -Gdpi=150 nf1_qsp_model.dot -o nf1_qsp_model.png
```

```r
# Run a scenario in R
library(mrgsolve); library(dplyr); library(ggplot2)
source("nf1_mrgsolve_model.R")           # builds `nf1_mod` + `scenarios`
res <- run_scenario("2_Selumetinib_25mgm2_BID_Ped_SPRINT",
                     scenarios[["2_Selumetinib_25mgm2_BID_Ped_SPRINT"]])
plot(res$time/24/7, res$PN_RESPONSE_PCT, type = "l")

# Launch the dashboard
shiny::runApp("nf1_shiny_app.R")
```

---

## Disease-specific quick-reference (Korean / English)

| 한국어 | English |
|---|---|
| 신경섬유종증 1형 | Neurofibromatosis type 1 (NF1) |
| 신경섬유종 (뉴로피브로민) | Neurofibromin |
| 총상신경섬유종 | Plexiform neurofibroma (PN) |
| 피부신경섬유종 | Cutaneous neurofibroma (cNF) |
| 시신경로 신경교종 | Optic pathway glioma (OPG) |
| 악성 말초신경초종양 | Malignant peripheral nerve sheath tumor (MPNST) |
| 카페오레 반점 | Café-au-lait macules |
| 홍채 리쉬 결절 | Lisch nodules |
| 경골 형성이상/거짓관절증 | Tibial dysplasia / pseudarthrosis |
| MEK 억제제 | MEK1/2 inhibitor (selumetinib, mirdametinib, trametinib, binimetinib) |

---

*Built by Claude Code Routine on 2026-07-01 as part of the QSP Disease Model
Library. See root [README.md](../README.md) for the full model gallery.*
