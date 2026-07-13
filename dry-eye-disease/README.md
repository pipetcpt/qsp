# Dry Eye Disease (DED, 건성안) — QSP Model

> Integrated Quantitative Systems Pharmacology model of dry eye disease,
> linking tear-film hyperosmolarity, meibomian-gland (evaporative) and
> lacrimal-gland (aqueous-deficient/Sjögren) dysfunction, the innate→adaptive
> ocular-surface inflammatory cascade (dendritic cell → Th17/IL-17A → MMP-9),
> epithelial barrier breakdown, and corneal neurosensory dysfunction — with
> the modern topical pharmacology stack (cyclosporine, lifitegrast,
> loteprednol, perfluorohexyloctane, varenicline nasal spray, diquafosol) and
> the clinical endpoints used to track disease (OSDI, TBUT, corneal staining,
> Schirmer test, tear osmolarity, MMP-9).

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`ded_qsp_model.dot`](ded_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`ded_qsp_model.svg`](ded_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`ded_qsp_model.png`](ded_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`ded_mrgsolve_model.R`](ded_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`ded_shiny_app.R`](ded_shiny_app.R) |
| 📚 References             | [`ded_references.md`](ded_references.md) |

---

## 1. Disease in one paragraph

Dry eye disease (DED) affects 5-50% of adults worldwide (TFOS DEWS II) and is
defined as a multifactorial disease of the ocular surface characterized by a
**loss of tear-film homeostasis**, accompanied by ocular symptoms in which
**tear-film instability and hyperosmolarity, ocular surface inflammation and
damage, and neurosensory abnormalities play etiological roles**. Two
overlapping phenotypes converge on the same vicious cycle: **evaporative DED**
(meibomian gland dysfunction → lipid-layer deficiency → increased evaporation)
and **aqueous-deficient DED** (lacrimal-gland dysfunction, often Sjögren-
associated → reduced tear production). Either route raises tear-film
**osmolarity**, which activates MAPK/NF-κB/NLRP3 signaling in ocular-surface
epithelium, triggering an **innate inflammatory cascade** (IL-1β, TNF-α,
MMP-9) that matures dendritic cells and primes an **adaptive Th17/Th1
response**; effector cytokines (IL-17A, IFN-γ) home back to the ocular
surface, causing **goblet-cell loss, epithelial barrier breakdown, and
further tear-film destabilization** — closing a self-perpetuating loop.
Chronic disease also damages **corneal sub-basal nerves**, producing
peripheral/central sensitization and **neuropathic ocular pain** that can
persist even as signs improve. Modern therapy targets specific nodes of this
cycle: **cyclosporine A** and **lifitegrast** dampen the DC/Th17 axis;
**loteprednol** provides rapid but time-limited broad anti-inflammatory
pulse therapy; **perfluorohexyloctane** is a non-drug evaporation barrier;
**varenicline nasal spray** stimulates reflex tearing via the trigeminal
parasympathetic pathway; and **diquafosol** (P2Y2 agonist, approved in
Japan/Korea) stimulates mucin and aqueous secretion directly.

## 2. Mechanistic clusters (11 in the DOT map, 100+ nodes)

1. Epidemiologic & environmental drivers (age, sex/androgen, screen time, CL wear, refractive surgery, autoimmune disease, BAK preservative toxicity)
2. Trilaminar tear-film homeostasis (lipid/aqueous/mucin layers, TBUT, evaporation)
3. Meibomian gland dysfunction (evaporative DED — keratinization, obstruction, gland dropout)
4. Lacrimal gland dysfunction (aqueous-deficient DED — Sjögren infiltrate, acinar apoptosis, reflex-arc impairment)
5. Hyperosmolar stress signaling (MAPK/NF-κB/NLRP3 inflammasome, ROS)
6. Innate inflammatory cascade (IL-1β/TNF-α/IL-6/IL-8/MMP-9, neutrophil recruitment, DC activation)
7. Adaptive immunity — self-perpetuating vicious cycle (DC maturation → Th17/Th1 → IL-17A/IFN-γ → ocular homing → chronic inflammation feed-forward)
8. Corneal neurosensory dysfunction (TRPM8 cold receptors, nerve loss, peripheral/central sensitization, neuropathic pain, photophobia)
9. Ocular-surface epithelial damage (goblet-cell loss, squamous metaplasia, apoptosis, barrier dysfunction, filamentary keratitis)
10. Topical & systemic drug PK/PD (cyclosporine, lifitegrast, loteprednol, perfluorohexyloctane, varenicline, diquafosol, punctal plugs, autologous serum, omega-3, thermal pulsation/IPL)
11. Clinical endpoints & biomarkers (OSDI, eye-dryness VAS, TBUT, Schirmer, corneal/conjunctival staining, tear osmolarity, MMP-9 POC, meibography, NEI-VFQ-25)

## 3. mrgsolve model (22 ODE compartments)

* **Local/topical drug exposure (8 compartments)** — cyclosporine A
  (0.05%/0.09%), lifitegrast 5%, loteprednol 0.25% (+ cumulative-exposure
  tracker for IOP-risk), perfluorohexyloctane, varenicline nasal spray (nasal
  depot + trigeminal-reflex effect compartment), diquafosol 3%. Topical
  ophthalmic PK is dominated by tear-turnover on a minutes scale, so each
  drug is represented as a day-resolution surrogate "sustained local
  exposure" compartment rather than a literal plasma PK model.
* **Tear-film & ocular-surface disease network (14 compartments)** — tear
  osmolarity, lipid-layer quality, aqueous production drive, mucin/goblet
  density, dendritic-cell activation, IL-17A/Th17 tone, MMP-9, epithelial
  barrier integrity, corneal nerve density, neuropathic pain, TBUT, Schirmer,
  corneal staining, OSDI — the last four double as clinical-endpoint readouts.

### 9 prebuilt scenarios

| # | Scenario | Calibration anchor |
|---|---|---|
| 1 | Untreated natural history | TFOS DEWS II pathophysiology framework |
| 2 | Cyclosporine 0.05% BID (Restasis) | Sall 2000; slow (mo 3-6) onset |
| 3 | Cyclosporine 0.09% BID (Cequa) | Tauber 2018 OTX-101 |
| 4 | Lifitegrast 5% BID (Xiidra) | OPUS-2/3 (Tauber 2015, Holland 2017) |
| 5 | Loteprednol 0.25% QID, 2-wk induction | rapid onset, short-course pulse |
| 6 | Perfluorohexyloctane QID (Miebo) | GOBI/MOJAVE (Tauber 2023, Sheppard 2023) |
| 7 | Varenicline nasal spray BID (Tyrvaya) | ONSET-1/2 (Wirta 2022) |
| 8 | Diquafosol 3%, 6x/day | Matsumoto 2012 |
| 9 | Combo: loteprednol induction + cyclosporine maintenance | common real-world "soak and soothe" regimen |

## 4. Shiny dashboard (8 tabs)

1. **Patient profile** — MGD/ADDE severity sliders, Sjögren flag, refractive-surgery flag, baseline OSDI.
2. **Drug exposure** — fractional drug-effect trajectories for all six agents.
3. **Tear-film / PD** — tear osmolarity, lipid-layer quality, aqueous production, mucin/goblet density.
4. **Inflammatory cascade** — dendritic-cell activation, IL-17A tone, MMP-9, epithelial barrier integrity.
5. **Clinical endpoints** — TBUT, Schirmer test, corneal staining, OSDI symptom score.
6. **Neurosensory** — corneal nerve density and neuropathic ocular pain score.
7. **Scenario comparison** — runs all 9 built-in scenarios with the chosen patient profile.
8. **Biomarkers / safety** — loteprednol cumulative-exposure IOP-risk score, evaporation-rate trajectory, endpoint summary table.

## 5. How to run

```bash
# 1) Render the mechanistic map
dot -Tsvg ded_qsp_model.dot -o ded_qsp_model.svg
dot -Tpng -Gdpi=150 ded_qsp_model.dot -o ded_qsp_model.png
```

```r
# 2) Simulate scenarios in R
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny",
                   "shinydashboard","DT"))
library(mrgsolve)
mod <- mread_cache("ded_mrgsolve_model.R")
out <- mrgsim(mod, ev(amt=1.0, cmt="CSA_TF", ii=0.5, addl=168), end=84)
plot(out, c("OSDI_score","TBUT_s","Corneal_staining","Schirmer_mm"))

# 3) Launch the dashboard
shiny::runApp("ded_shiny_app.R")
```

## 6. Key clinical anchors used during calibration

| Endpoint | Comparator | Expected effect |
|---|---|---|
| Schirmer responder (≥10mm gain) | Cyclosporine 0.05% vs vehicle, 6 mo | ~15% vs ~5% (Sall 2000) |
| Schirmer responder | Cyclosporine 0.09% (Cequa), wk 12 | ~17% |
| Corneal fluorescein staining | Lifitegrast vs vehicle, wk 6 | −0.9 to −1.1 (OPUS-2/3) |
| Total corneal staining / eye-dryness VAS | Perfluorohexyloctane vs saline, wk 8 | tCFS & VAS significantly improved (GOBI/MOJAVE) |
| Schirmer score | Varenicline nasal spray vs vehicle, wk 4 | ≥10mm gain more frequent (ONSET-2) |
| TBUT / fluorescein staining | Diquafosol 3% vs placebo, wk 4-6 | significant improvement (Matsumoto 2012) |
| IOP / cataract risk | Loteprednol (any topical steroid) | rises with cumulative exposure — short-course use only |

## 7. Caveats

* Designed for **research, education, and hypothesis generation** — not
  clinical decision support.
* Topical ocular drug "exposure" compartments are pragmatic day-resolution
  surrogates for sustained local tissue concentration, not literal tear-film
  or aqueous-humor pharmacokinetics (which occur on a minutes-to-hours scale
  dominated by tear turnover).
* Perfluorohexyloctane and diquafosol are modeled with **no direct
  anti-inflammatory arm**, consistent with their non-immunomodulatory
  mechanisms (physical evaporation barrier and P2Y2 secretagogue,
  respectively).
* Corneal nerve density and neuropathic pain evolve on a much slower
  (months-scale) timescale than the inflammatory/tear-film compartments; the
  9 built-in scenarios (12-week default horizon) will show only partial
  recovery of these compartments.
* Loteprednol's IOP-risk score is a simplified linear function of cumulative
  exposure intended to flag the real-world need to limit corticosteroid
  courses to 2-4 weeks — it is not a validated IOP-elevation model.

## 8. License

Inherits the repository [LICENSE](../LICENSE).
