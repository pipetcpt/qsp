# Malignant Pleural Mesothelioma (MPM, 악성 흉막 중피종) — QSP Model

> Integrated Quantitative Systems Pharmacology model of asbestos-driven
> malignant pleural mesothelioma, linking fiber-triggered frustrated
> phagocytosis and HMGB1-NLRP3 chronic inflammation to the BAP1/CDKN2A-
> p16-ARF/NF2-Hippo-YAP tumor-suppressor-loss landscape and a two-clone
> (epithelioid/sarcomatoid) tumor-burden network — together with the full
> modern therapeutic stack (cisplatin/pemetrexed ± bevacizumab, first- and
> second-line nivolumab ± ipilimumab, anti-mesothelin antibody-drug
> conjugate anetumab ravtansine, BAP1-synthetic-lethal PARP inhibitor
> rucaparib, Tumor-Treating Fields, and surgery/pleurodesis) and downstream
> clinical/biomarker consequences (pleural effusion, mesothelin/SMRP,
> myelosuppression, nephrotoxicity, ECOG, survival).

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`mpm_qsp_model.dot`](mpm_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`mpm_qsp_model.svg`](mpm_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`mpm_qsp_model.png`](mpm_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`mpm_mrgsolve_model.R`](mpm_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`mpm_shiny_app.R`](mpm_shiny_app.R) |
| 📚 References             | [`mpm_references.md`](mpm_references.md) |

---

## 1. Disease in one paragraph

Malignant pleural mesothelioma is an aggressive, long-latency (20-50 year)
cancer of the pleural mesothelial surface, overwhelmingly caused by
inhaled amphibole (crocidolite/amosite) or chrysotile asbestos fibers.
Macrophages that attempt but fail to fully phagocytose long fibers
("frustrated phagocytosis") generate ROS/RNS and DNA damage; necrotic
mesothelial cells release **HMGB1**, driving **NLRP3-inflammasome** and
**TNF-α/NF-κB** signaling that paradoxically rescues genomically damaged
cells from apoptosis — a chronic-inflammation "mutagenic sequence" that
seeds chromosomal instability. The resulting tumors are defined by loss of
**BAP1** (~60%, including a germline tumor-predisposition syndrome),
**CDKN2A/p16-ARF** homozygous deletion (~70-90%), and **NF2/Merlin**
loss (~40-50%) that dysregulates the **Hippo pathway**, driving nuclear
**YAP/TAZ-TEAD** oncogenic transcription; TP53 mutation is comparatively
rare. Three histologic subtypes — epithelioid (best prognosis, mesothelin-
high), sarcomatoid (worst, EMT-driven, drug-resistant), and biphasic — track
a shared molecular landscape. Standard care combines platinum/antifolate
chemotherapy (± bevacizumab), first-line dual checkpoint blockade
(nivolumab+ipilimumab, CheckMate 743), and multimodality
surgery/radiotherapy, with mesothelin-targeted and BAP1-synthetic-lethal
(PARP inhibitor) strategies in development.

## 2. Mechanistic clusters (14 in the DOT map, 115 nodes)

1. Fiber exposure & etiology (occupational/environmental asbestos, latency)
2. Frustrated phagocytosis, ROS/RNS & HMGB1-NLRP3 inflammation
3. Genomic landscape — BAP1/CDKN2A/NF2/TP53 tumor-suppressor loss
4. BAP1 loss — chromatin (PRC1-H2AK119Ub), HR-repair deficiency, Warburg shift
5. NF2/Merlin loss — Hippo pathway & YAP/TAZ-TEAD oncogenic transcription
6. Cell-cycle checkpoint (p16-Rb, p14ARF-p53) & apoptosis dysregulation
7. Histologic subtypes, tumor microenvironment & immune evasion
8. Angiogenic axis (HIF-1α/VEGF-A/VEGFR2)
9. Clinical phenotype, local invasion & IMIG/UICC TNM staging
10. Drug PK/PD — platinum/antifolate/anti-VEGF chemotherapy
11. Drug PK/PD — immune checkpoint inhibitors (nivolumab/ipilimumab/avelumab)
12. Drug PK/PD — mesothelin ADC (anetumab ravtansine) & BAP1-synthetic-lethal PARP inhibition (rucaparib)
13. Surgery, radiotherapy, Tumor-Treating Fields & pleurodesis
14. Clinical endpoints & biomarkers (SMRP, fibulin-3, ORR/PFS/OS, QoL)

## 3. mrgsolve model (25 ODE compartments)

* **Drug PK (7 agents, 15 compartments)** — cisplatin, pemetrexed,
  bevacizumab, nivolumab, ipilimumab, anetumab ravtansine (ADC), and oral
  rucaparib (literature PK anchors in `mpm_references.md`).
* **Two-clone tumor network (2 compartments)** — epithelioid and
  sarcomatoid tumor burden, each with distinct intrinsic growth rate and
  drug-sensitivity (chemo/IO/ADC/PARPi) gated by tumor-biology flags
  (`HIST_SARC_FRAC`, `BAP1_DEFICIENT`, `MSLN_HIGH`, `PDL1_HIGH`).
* **TME/organ readouts (5 compartments)** — pleural effusion volume,
  mesothelin (SMRP) biomarker, creatinine clearance (cisplatin
  nephrotoxicity), and a simplified 3-compartment Friberg-like
  myelosuppression chain (proliferating → transit → circulating ANC).
* **Clinical readouts (2 compartments)** — ECOG performance-status
  surrogate and cumulative mortality hazard (→ modeled survival
  probability).

### 10 prebuilt scenarios

| # | Scenario | Calibration anchor |
|---|---|---|
| 1 | Untreated / best supportive care        | Natural history |
| 2 | Cisplatin+Pemetrexed q3w x6              | Vogelzang 2003 EMPHACIS |
| 3 | Cisplatin+Pemetrexed+Bevacizumab q3w     | Zalcman 2016 MAPS |
| 4 | Nivolumab+Ipilimumab 1L                  | Baas 2021 CheckMate 743 |
| 5 | Nivolumab monotherapy 2L                 | Fennell 2021 CONFIRM |
| 6 | Nivolumab+Ipilimumab 2L                  | Scherpereel 2019 MAPS2 |
| 7 | TTFields + Cisplatin/Pemetrexed          | Ceresoli 2019 STELLAR |
| 8 | Anetumab ravtansine (mesothelin-high)    | Hassan 2020 Ph1 |
| 9 | Rucaparib (BAP1-deficient)               | Fennell 2021 MiST1 |
| 10 | Surgery (P/D) + chemo + talc pleurodesis | Treasure 2011 MARS / Lim 2024 MARS2 |

## 4. Shiny dashboard (8 tabs)

1. **Patient profile** — tumor-biology flags (sarcomatoid fraction, BAP1,
   mesothelin, PD-L1) & baseline ECOG.
2. **Drug PK** — log-scale plasma concentrations for the seven tracked drugs.
3. **Molecular pathway** — epithelioid vs. sarcomatoid clone trajectories
   and total tumor burden vs. carrying capacity.
4. **Immune/TME & effusion** — pleural effusion volume and ECOG trend.
5. **Clinical endpoints** — tumor volume, modeled survival probability,
   endpoint summary table.
6. **Scenario comparison** — runs all 10 scenarios with the chosen profile.
7. **Safety** — circulating ANC (myelosuppression) and creatinine clearance.
8. **Biomarkers** — serum mesothelin (SMRP) trajectory.

## 5. How to run

```bash
# 1) Render the mechanistic map
dot -Tsvg mpm_qsp_model.dot -o mpm_qsp_model.svg
dot -Tpng -Gdpi=150 mpm_qsp_model.dot -o mpm_qsp_model.png
```

```r
# 2) Simulate scenarios in R
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny",
                   "shinydashboard","DT"))
library(mrgsolve)
mod <- mread_cache("mpm_mrgsolve_model.R")
out <- mod %>% mrgsim(events = seq(ev(amt=240, cmt="NIVO_C1", ii=14, addl=25),
                                    ev(amt=63,  cmt="IPI_C1",  ii=42, addl=8)),
                       end = 252)
plot(out, c("Tumor_total","Survival_prob","Effusion_vol","SMRP_biomarker"))

# 3) Launch the dashboard
shiny::runApp("mpm_shiny_app.R")
```

## 6. Key clinical anchors used during calibration

| Endpoint | Comparator | Expected effect |
|---|---|---|
| OS | Cisplatin+Pemetrexed vs cisplatin alone | 12.1 vs 9.3 mo (EMPHACIS) |
| OS | +Bevacizumab vs chemo alone | 18.8 vs 16.1 mo (MAPS) |
| OS | Nivolumab+Ipilimumab 1L vs chemo | 18.1 vs 14.1 mo (CheckMate 743) |
| OS / PFS | Nivolumab 2L vs placebo | 9.2 vs 6.6 mo / PFS HR 0.67 (CONFIRM) |
| Disease control | Nivolumab±Ipilimumab 2L | ~44-50% at wk 12 (MAPS2) |
| OS | TTFields+chemo vs historical control | ~18.2 vs ~12 mo (STELLAR) |
| Disease control | Rucaparib, BAP1-deficient | ~58% at wk 12 (MiST1) |
| OS | Surgery(P/D)+chemo vs chemo alone | 19.3 vs 24.8 mo — surgery arm **worse** (MARS2) |

## 7. Caveats

* Designed for **research, education, and hypothesis generation** — not
  clinical decision support.
* PK parameters for all seven drugs are literature-typical population
  surrogates reconstructed from trial methods sections and public
  regulatory assessment reports, not a validated population-PK fit.
* The epithelioid/sarcomatoid two-clone structure is a simplification of a
  continuous histologic and molecular spectrum; biphasic tumors are
  represented as a mixed fraction (`HIST_SARC_FRAC`).
* Surgery and talc pleurodesis are modeled as instantaneous fractional
  reductions applied in the Shiny scenario logic rather than as ODE
  compartments, consistent with their nature as discrete procedural
  interventions.
* Immune-related adverse events, radiotherapy dosimetry, and detailed
  irAE-grade kinetics are referenced in the bibliography but intentionally
  **not** modeled as separate ODE compartments.

## 8. License

Inherits the repository [LICENSE](../LICENSE).
