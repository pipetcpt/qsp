# NAFLD/NASH QSP Model

A quantitative systems pharmacology (QSP) model of **non-alcoholic fatty liver disease /
steatohepatitis (NAFLD/NASH)** in [`mrgsolve`](https://mrgsolve.org). It simulates hepatic fat,
inflammation, fibrosis and insulin resistance, and the effect of four therapies —
**resmetirom, obeticholic acid (OCA), semaglutide, empagliflozin** — with drug effect sizes
anchored to published clinical trials.

![Mechanism map](nafld_qsp_model.png)

## Files
| File | Purpose |
|---|---|
| `nafld_mrgsolve_model.R` | the QSP model (11 disease ODEs + 4 drug PK/PD modules, 6 scenarios) |
| `nafld_shiny_app.R` | interactive dashboard (embeds the same model) |
| `nafld_model_design_brief.md` | design rationale + recalibration & adversarial-review notes |
| `nafld_references.md` | references (DOI/PMID) |
| `nafld_qsp_model.dot/.svg/.png` | mechanism map (source + rendered) |
| `validation/validate.R` | reproducibility checks (exit 0 = all pass) |
| `validation/figures/` | validation evidence (overlay, before/after recalibration) |

## Validated effect sizes (model vs published, placebo-corrected)
| Drug | Model | Trial | Source |
|---|---|---|---|
| Resmetirom 100 mg — liver fat | **−34.0%** | −33.9% (wk52) | MAESTRO-NAFLD-1 (Harrison, *Nat Med* 2023) |
| Empagliflozin 10 mg — liver fat | **−24.6%** | −24.7% (between-group) | E-LIFT (Kuchay, *Diabetes Care* 2018) |
| Semaglutide 2.4 mg — body weight | **−11.3%** | −10.5% / −14.9% | ESSENCE 2025 / STEP-1 2021 |
| OCA 25 mg — fibrosis | ↓ (direction) | 23% vs 12% improvers | REGENERATE (Younossi, *Lancet* 2019) |

Structural anchors: liver-fat source split (NEFA 59% / DNL+diet 41%, **Donnelly 2005**);
NAS steatosis cutoffs (5 / 33 / 66 %, **Kleiner 2005**).

## Run
```r
install.packages(c("mrgsolve","dplyr","ggplot2","tidyr","patchwork",   # model
                   "shiny","bslib","plotly","DT"))                       # dashboard
```
```bash
Rscript nafld_mrgsolve_model.R          # run model + generate scenario figures
Rscript validation/validate.R           # reproducibility checks (placebo flat + effect sizes)
Rscript -e "shiny::runApp('nafld_shiny_app.R')"   # dashboard
```
(R 4.5/4.6 + Rtools on Windows; `mrgsolve` compiles C++ at model build.)

## What changed in this contribution
- **Resmetirom** recalibrated to MAESTRO-NAFLD-1 wk52 *placebo-corrected* value
  (−40% → **−34%**; DNL inhibition 0.40→0.30, efflux gain 0.60→0.48).
- **Empagliflozin**: added a direct hepatic-fat efflux term (`WEMP_LF`, calibrated to E-LIFT
  *placebo-corrected* −24.7%) and fixed a latent near-inert PK scaling (`EC50_EMP` 0.15→0.015).
  Liver-fat effect: **−1% → −24.6%**.
- Calibration anchored to **placebo-corrected (between-group)** trial effects throughout (the
  QSP placebo arm is flat, so within-arm changes would over-credit the drug).
- Citation/provenance corrections in the model header; design brief documents limitations
  surfaced by an adversarial multi-lens QSP review.
- Added a `validation/` harness so the effect sizes are reproducible (`Rscript validation/validate.R`).

## Limitations (see design brief)
- Steady-state model matches the trial **plateau (wk52)**, not the within-study time-course.
- NAS inflammation/ballooning sub-scores are crude proxies → use NAS **deltas**, not absolutes.
- OCA safety (pruritus, LDL) and the NASH-**resolution** endpoint are not modeled; empagliflozin
  weight/IR coefficients were not re-tuned after the EC50 fix.
- A continuous-output → **responder-rate (% NASH resolution)** bridge is future work.
