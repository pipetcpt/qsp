# Postural Orthostatic Tachycardia Syndrome (POTS, 체위성 기립성 빈맥 증후군) — QSP Model

> Integrated Quantitative Systems Pharmacology model of POTS, linking
> orthostatic venous pooling and baroreflex-arc control to excessive
> standing tachycardia across four overlapping pathophysiological subtypes
> — neuropathic (partial sympathetic denervation), hyperadrenergic (excess
> norepinephrine/NET dysfunction), hypovolemic (reduced plasma volume/RAAS
> paradox), and autoimmune/post-viral (incl. post-COVID, adrenergic/
> muscarinic autoantibodies) — together with the modern pharmacology stack
> (midodrine, beta-blockers, fludrocortisone, ivabradine, pyridostigmine,
> droxidopa) and non-pharmacologic therapy (salt/fluid loading, compression,
> exercise training).

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`pots_qsp_model.dot`](pots_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`pots_qsp_model.svg`](pots_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`pots_qsp_model.png`](pots_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`pots_mrgsolve_model.R`](pots_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`pots_shiny_app.R`](pots_shiny_app.R) |
| 📚 References             | [`pots_references.md`](pots_references.md) |

---

## 1. Disease in one paragraph

POTS is a form of chronic orthostatic intolerance defined by an excessive
heart-rate increase on standing (≥30 bpm in adults, ≥40 bpm in adolescents,
within 10 minutes of active standing or head-up tilt) **without** the
sustained blood-pressure drop that defines orthostatic hypotension
(Sheldon 2015 Heart Rhythm consensus). Upright posture displaces
~500-1000 mL of blood into leg and splanchnic venous capacitance beds;
normally the baroreflex arc (carotid/aortic baroreceptors → NTS → RVLM/
vagal nuclei) restores venous tone and cardiac filling within seconds.
In POTS this compensation is insufficient or the tachycardic response is
disproportionate, through four overlapping mechanisms: **neuropathic**
(length-dependent small-fiber sympathetic denervation impairs distal
venoconstriction), **hyperadrenergic** (NET dysfunction/excess synaptic
norepinephrine drives standing plasma NE >600 pg/mL), **hypovolemic**
(~10-15% reduced blood/plasma volume with a paradoxically low renin-
aldosterone response), and **autoimmune/post-viral** (agonist β1/β2 and
antagonist M2/M4 autoantibodies, or ganglionic α3-nAChR antibodies,
frequently triggered by viral illness including SARS-CoV-2/long-COVID).
Chronic activity avoidance drives cardiovascular deconditioning (small,
under-filled LV, reduced skeletal-muscle pump), which worsens tachycardia
and closes a vicious cycle. Treatment combines volume expansion (salt/
fluid, fludrocortisone), vasoconstriction (midodrine, compression),
selective heart-rate control (beta-blockers, ivabradine), ganglionic
augmentation (pyridostigmine), norepinephrine repletion (droxidopa, for
severe neuropathic cases), and — with the strongest evidence for durable
benefit — supervised exercise-training programs that reverse deconditioning.

## 2. Mechanistic clusters (16 in the DOT map, 114 nodes)

1. Orthostatic stress & venous pooling (gravity, venous capacitance, muscle pump)
2. Baroreflex arc & central autonomic control (baroreceptors → NTS/RVLM → vagal/sympathetic outflow)
3. Neuropathic subtype (distal sympathetic denervation, IENFD, QSART, splanchnic vasoconstrictor failure)
4. Hyperadrenergic subtype (NET deficiency, standing NE >600 pg/mL, α1/β1 excess)
5. Hypovolemic subtype (RAAS paradox, blunted AVP, renal salt wasting)
6. Autoimmune/post-viral subtype & overlap syndromes (adrenergic/muscarinic/ganglionic autoantibodies, MCAS, hEDS)
7. Cardiovascular hemodynamics (supine/standing HR, cardiac output, TPR, BP)
8. Cerebral blood flow & neuro-symptoms (autoregulation, presyncope, brain fog)
9. Deconditioning & skeletal muscle (cardiac atrophy, exercise intolerance)
10. GI/bladder/comorbid dysautonomia (gastroparesis, IBS, bladder, ME/CFS, migraine)
11. Drug PK/PD — midodrine (α1-agonist prodrug/desglymidodrine)
12. Drug PK/PD — beta-blocker (propranolol/bisoprolol)
13. Drug PK/PD — fludrocortisone (mineralocorticoid, plasma-volume expansion)
14. Drug PK/PD — ivabradine (If-channel blocker, BP-sparing)
15. Drug PK/PD — pyridostigmine & droxidopa (ganglionic AChE inhibition / NE precursor)
16. Non-pharmacologic therapy, clinical endpoints & biomarkers (stand test, tilt table, COMPASS-31, QoL, plasma NE/aldosterone/volume)

## 3. mrgsolve model (23 ODE compartments)

* **Drug PK (6 drug classes, 13 compartments)** — midodrine + active
  metabolite desglymidodrine, propranolol, fludrocortisone, ivabradine,
  pyridostigmine, droxidopa (literature-informed surrogate PK anchors in
  `pots_references.md`).
* **Autonomic/volume disease network (5 compartments)** — relative plasma
  volume, total peripheral resistance, standing norepinephrine, baroreflex
  sensitivity, deconditioning index.
* **Hemodynamic & clinical readouts (5 compartments)** — supine heart
  rate, standing heart rate (10-min stand-test surrogate), cerebral
  blood-flow index, COMPASS-31 composite score, quality-of-life score.

### 10 prebuilt scenarios

| # | Scenario | Calibration anchor |
|---|---|---|
| 1 | Untreated natural history (mixed subtype)                 | Sheldon 2015 HRS consensus |
| 2 | Salt/fluid loading + compression garments                 | Fu 2011 Heart Rhythm non-pharm review |
| 3 | Propranolol 20 mg PO TID (hyperadrenergic subtype)        | Raj 2009 Circulation |
| 4 | Midodrine 10 mg PO TID (neuropathic subtype)               | Sutton/Grubb open-label series |
| 5 | Fludrocortisone 0.1 mg PO QD + salt (hypovolemic subtype) | Raj 2005 Circulation blood-volume study |
| 6 | Ivabradine 5 mg PO BID (refractory tachycardia)           | Ruzieh 2017 / Moon 2018 case series |
| 7 | Pyridostigmine 60 mg PO TID (mild-moderate augmentation)  | Kanjwal 2011 Cardiol J |
| 8 | Droxidopa 100 mg PO TID (severe neuropathic subtype)      | Neurogenic-OI program, extrapolated |
| 9 | Combination (BB + midodrine + compression + exercise)     | Multimodal standard-of-care |
| 10 | Recumbent/semi-recumbent exercise-training program alone | Fu 2010/2011 "Levine protocol" |

## 4. Shiny dashboard (8 tabs)

1. **Patient profile** — adjustable age/weight/subtype & baseline severity.
2. **Drug PK** — log-scale plasma concentrations for the seven tracked drug moieties.
3. **Autonomic / Volume PD** — relative plasma volume, TPR, standing NE, baroreflex sensitivity.
4. **Hemodynamics / Stand test** — supine vs. standing HR and the ΔHR diagnostic threshold.
5. **Deconditioning / CBF** — cerebral blood-flow index and deconditioning trajectory.
6. **Symptom & QoL** — COMPASS-31 composite score and quality-of-life score.
7. **Scenario comparison** — runs all 10 scenarios with the chosen profile.
8. **References** — key citations and link to the full bibliography.

## 5. How to run

```bash
# 1) Render the mechanistic map
dot -Tsvg pots_qsp_model.dot -o pots_qsp_model.svg
dot -Tpng -Gdpi=150 pots_qsp_model.dot -o pots_qsp_model.png
```

```r
# 2) Simulate scenarios in R
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny",
                   "shinydashboard","DT"))
library(mrgsolve)
mod <- mread_cache("pots_mrgsolve_model.R")
out <- mrgsim(mod, ev(amt=10, cmt="GUT_MID", ii=8/24, addl=3*14), end=14*24)
plot(out, c("DeltaHR_bpm","HR_standing_bpm","PlasmaVolume_idx","COMPASS31_score"))

# 3) Launch the dashboard
shiny::runApp("pots_shiny_app.R")
```

## 6. Key clinical anchors used during calibration

| Endpoint | Comparator | Expected effect |
|---|---|---|
| ΔHR (10-min stand) | Diagnostic threshold | ≥30 bpm adult / ≥40 bpm adolescent (Sheldon 2015) |
| Standing HR | Propranolol 20 mg vs 80 mg vs placebo | low-dose (20 mg) most effective (Raj 2009) |
| Standing NE | Hyperadrenergic-subtype criterion | >600 pg/mL (Grubb/Low/Streeten) |
| Blood/plasma volume | POTS vs healthy controls | ~10-15% lower, RAAS paradox (Raj 2005) |
| Standing HR & QoL | Exercise training, 3 months | comparable to pharmacotherapy in many patients (Fu 2010/2011) |
| Standing HR | Ivabradine case series | selective HR reduction without BP effect |
| Standing HR | Pyridostigmine open-label | modest reduction; GI cramping dose-limiting |

## 7. Caveats

* Designed for **research, education, and hypothesis generation** — not
  clinical decision support.
* All PK/PD parameters are population means or literature-informed
  surrogates; inter-individual variability needs `omega()` blocks for
  population simulations.
* The four subtypes (neuropathic/hyperadrenergic/hypovolemic/autoimmune)
  overlap extensively in real patients; the `SUBTYPE` covariate selects a
  dominant-mechanism parameterization rather than modeling continuous
  co-occurrence.
* Supine-hypertension and other drug adverse effects are annotated in the
  mechanistic map but not tracked as explicit safety-endpoint compartments.

## 8. License

Inherits the repository [LICENSE](../LICENSE).
