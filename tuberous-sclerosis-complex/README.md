# Tuberous Sclerosis Complex (TSC) — QSP Model

> Multi-system, mTORC1-driven hamartomatous disease (TSC1/TSC2 loss-of-function) — SEGA · epilepsy · TAND · renal AML · pulmonary LAM · cutaneous angiofibromas · cardiac rhabdomyomas. Pharmacology centred on **everolimus**, **sirolimus**, **topical sirolimus**, **vigabatrin**, **cannabidiol**, with **EXIST-1/2/3 · MILES · GWPCARE6** trial anchors.

## Files
| File | Content |
|---|---|
| `tsc_qsp_model.dot` / `.svg` / `.png` | Mechanistic map · 16 clusters · 130+ nodes |
| `tsc_mrgsolve_model.R`                | mrgsolve model · 24 ODE compartments · 8 scenarios |
| `tsc_shiny_app.R`                     | 8-tab Shiny dashboard |
| `tsc_references.md`                   | 50 curated references with PubMed/DOI |

## Mechanistic map clusters
1. **Genetics** — TSC1 (9q34) / TSC2 (16p13.3); germline vs somatic two-hit; TSC2/PKD1 contiguous gene syndrome; mosaicism in "NMI" cases.
2. **mTORC1 core** — Hamartin · Tuberin · TBC1D7 → Rheb-GAP → mTORC1 → S6K1 / 4E-BP1 / ULK1 / TFEB; IRS-1 / Akt-S473 negative feedback.
3. **Upstream regulators** — Insulin/IGF-1 · PI3K-Akt; AMPK; amino-acid Ragulator; REDD1; ERK/RSK; GSK3β; Wnt; TNF-α/IKKβ.
4. **mTORC2** — RICTOR · mSIN1 · feedback Akt-S473 / SGK1 / FOXO; rebound under chronic mTORi.
5. **Neurological** — cortical tubers, SEN/SEGA → hydrocephalus; infantile spasms · focal seizures · LGS-like; TAND, autism (40–50%), intellectual disability (~50%); GABA / glutamate imbalance; eIF4E-FMRP synaptic plasticity.
6. **Renal** — AML (70–80%), aneurysm-driven Wunderlich haemorrhage, cysts, PKD1 contiguous ADPKD-like, chromophobe RCC.
7. **Pulmonary LAM** — TSC-LAM in 30–40% of adult women; HMB-45⁺/ERα⁺ LAM cells; MMP-driven cysts; VEGF-D > 800 pg/mL biomarker; estrogen amplification.
8. **Cutaneous** — facial angiofibromas, shagreen patch, ash-leaf macules, ungual fibromas, confetti lesions, fibrous cephalic plaque.
9. **Cardiac** — rhabdomyomas (>50% antenatal), WPW pre-excitation, LVOT/RVOT obstruction, postnatal regression.
10. **Ophthalmologic** — retinal astrocytic hamartoma (30–50%), achromic iris patch.
11. **Everolimus PK/PD** — F = 0.30, V1 = 110 L, CL = 8.8 L/h, t½ ~28 h; FKBP12 · everolimus complex → mTORC1 allosteric inhibition (IC50 ~2 nM).
12. **Sirolimus PK/PD** — F ~15%, t½ ~62 h; topical 0.1–1% for facial angiofibroma.
13. **Antiseizure / TAND therapy** — Vigabatrin GABA-T inhibitor (first-line infantile spasms), ACTH/prednisolone, **Cannabidiol** (Epidiolex; GWPCARE6), clobazam; CBD-CYP3A4 DDI raises EVE AUC.
14. **mTORi adverse effects** — stomatitis, dyslipidaemia, hyperglycaemia, infections, non-infectious pneumonitis, proteinuria, cytopenias.
15. **Clinical endpoints** — 2012 / 2021 consensus criteria, SEGA volume MRI, AML longest diameter, FASI/EFASI, seizure freq/28 d, QOLIE, TAND-L, FEV1.
16. **Pivotal trials** — EXIST-1/2/3, MILES, GWPCARE6, TOSCA registry.

## ODE compartments (mrgsolve)
| # | Compartment | Notes |
|---|---|---|
| 1 | EVE_GUT       | Everolimus depot |
| 2 | EVE_C         | Everolimus central (mg) |
| 3 | EVE_P         | Everolimus peripheral |
| 4 | SIR_GUT       | Sirolimus depot |
| 5 | SIR_C         | Sirolimus central |
| 6 | SIR_P         | Sirolimus peripheral |
| 7 | VGB_GUT       | Vigabatrin depot |
| 8 | VGB_C         | Vigabatrin central |
| 9 | CBD_GUT       | Cannabidiol depot |
| 10 | CBD_C        | CBD central |
| 11 | CBD_P        | CBD peripheral |
| 12 | TSIR_SKIN    | Topical sirolimus skin mass |
| 13 | MTOR_ACT     | mTORC1 activity (normalised) |
| 14 | SEGA         | SEGA volume (cm³) |
| 15 | AML          | Renal AML longest diameter (cm) |
| 16 | SKIN         | Facial angiofibroma FASI (0–100) |
| 17 | FEV1         | LAM FEV1 (%pred) |
| 18 | GABA         | GABA tone (normalised) |
| 19 | SZ           | Seizures / 28 d |
| 20–24 | H_STOMAT / H_LIPID / H_PNEU / H_VFD / H_HEPAT | Cumulative AE hazards |

## Treatment scenarios
1. **Untreated** natural history
2. **Everolimus 4.5 mg/d PO** (target trough 5–15 ng/mL, EXIST-1/2/3)
3. **Sirolimus 2 mg/d PO** (trough 6–14 ng/mL, MILES)
4. **Vigabatrin 1000 mg BID** (infantile spasms)
5. **Cannabidiol 25 mg/kg/d** (GWPCARE6)
6. **Everolimus + Vigabatrin + CBD** (refractory paediatric)
7. **Topical sirolimus 1% QD** (facial angiofibromas)
8. **Everolimus + topical sirolimus** (systemic + local)

## Calibration anchors
- **EXIST-1** — everolimus reduces SEGA volume ≥ 50% in 35% of patients at 6 mo (Franz NEJM 2013).
- **EXIST-2** — everolimus reduces AML diameter ≥ 50% in 42% at 12 mo (Bissler Lancet 2013).
- **EXIST-3** — everolimus 9 ng/mL trough → seizure-freq −40%, 15 ng/mL → −39% (French Lancet 2016).
- **MILES** — sirolimus stabilises FEV1 (vs −12 mL/yr placebo) at 12 mo (McCormack NEJM 2011).
- **GWPCARE6** — CBD 25 / 50 mg/kg/d → seizure-freq −48% / −47% at 16 wk (Thiele JAMA Neurol 2021).
- **Topical sirolimus** — facial angiofibroma severity reduces 40–60% at 6 mo (Wataya-Kaneda 2017, Koenig 2018).

## How to run

```r
# install.packages(c("mrgsolve","shiny","dplyr","tidyr","ggplot2"))
library(shiny)
shiny::runApp("tsc_shiny_app.R")
```

Render the map (Graphviz):

```bash
dot -Tsvg tsc_qsp_model.dot -o tsc_qsp_model.svg
dot -Tpng -Gdpi=150 tsc_qsp_model.dot -o tsc_qsp_model.png
```
