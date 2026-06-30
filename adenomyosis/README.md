# Adenomyosis (자궁선근증) — QSP Model

> Integrated Quantitative Systems Pharmacology model of uterine adenomyosis
> covering junctional-zone biology, local estrogen excess / progesterone
> resistance, neuroangiogenic pain, heavy menstrual bleeding (HMB), and the
> full modern pharmacology stack (GnRH antagonists with add-back, GnRH agonists,
> dienogest, LNG-IUS, aromatase inhibitors, SPRMs, NSAIDs, tranexamic acid,
> and surgical/interventional options).

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`adeno_qsp_model.dot`](adeno_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`adeno_qsp_model.svg`](adeno_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`adeno_qsp_model.png`](adeno_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`adeno_mrgsolve_model.R`](adeno_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`adeno_shiny_app.R`](adeno_shiny_app.R) |
| 📚 References             | [`adeno_references.md`](adeno_references.md) |

---

## 1. Disease in one paragraph

Adenomyosis is the presence of **endometrial glands and stroma within the
myometrium**, surrounded by hyperplastic/hypertrophic smooth-muscle.
Clinically it presents as **heavy menstrual bleeding (HMB)**, severe
**dysmenorrhea**, chronic pelvic pain, dyspareunia, and **subfertility / poor
obstetric outcome**. It is now best understood as a **junctional-zone disease**:
repeated peri-menstrual micro-trauma of the inner-myometrial *archi-myometrium*
plus disordered stem-cell biology drive basalis invagination, EMT, local
estrogen excess (aromatase induction by PGE2), progesterone resistance,
neuroangiogenesis, and periglandular fibrosis. Modern oral GnRH antagonists
with hormonal add-back, dienogest, and the LNG-IUS now provide
**uterus-sparing** alternatives to hysterectomy.

## 2. Mechanistic clusters (≥15 in the DOT map)

1. Etiogenesis (TIAR, basalis invagination, stem-cell, Müllerian remnant, EMT)
2. Junctional-zone hyperperistalsis & MUSA/MRI features
3. Local estrogen excess (StAR · CYP11A1 · CYP19A1 / aromatase · STS · 17βHSDs · ERα/β)
4. Progesterone resistance (PR-A↓ / PR-B imbalance · NF-κB co-repression · HOXA10)
5. Inflammation & prostaglandins (macrophages · IL-1β/6/8/TNF · NF-κB · COX-2 · PGE2/F2α)
6. Oxidative stress · free iron · Fenton · lipid peroxidation · NRF2
7. Fibrosis & EMT (TGF-β · SMAD2/3 · CTGF · collagen I/III · MMP/TIMP · myofibroblast)
8. Neuroangiogenesis & pain (HIF-1α · VEGF · NGF · TrkA · PIEZO2 · DRG/spinal · CNS sensitization)
9. HMB · fibrinolysis (PA/plasmin · PAI-1) · iron-deficiency anemia · QoL
10. Smooth-muscle hyperplasia/hypertrophy (OXTR · Ca²⁺ · MLCK · uterine enlargement)
11. Subfertility / implantation / obstetric outcome (HOXA10/LIF)
12. HPO axis (kisspeptin → GnRH → FSH/LH → ovary → systemic E2/P4)
13. Drug pharmacology (GnRH antagonists/agonists, progestins, AIs, SPRMs, NSAIDs, TXA, surgery)
14. Hypoestrogenic safety (BMD · hot flushes · lipids · LFTs · add-back)
15. Clinical endpoints & biomarkers (PBAC · VAS · NPRS · DPP-Q · PGIC · EHP-30 · JZ-MRI · CA-125 · Hb · ferritin)
16. Patient/population modifiers (age, parity, prior surgery, BMI, GWAS loci)

## 3. mrgsolve model (23 ODE compartments)

* **Drug PK (8 drugs, 16 compartments)** — relugolix-CT, elagolix, leuprolide
  depot, dienogest, LNG-IUS reservoir + systemic, letrozole, ibuprofen,
  tranexamic acid (literature PK anchors in `adeno_references.md`).
* **HPO axis (4 compartments)** — FSH, LH, systemic E2, P4 with kinetic
  feedback that captures GnRH antagonist suppression, GnRH-agonist flare→
  desensitization, and progestin-induced ovulation block.
* **Local lesion biology (7 compartments)** — aromatase, COX-2/PGE2,
  intralesional E2, TGF-β tone, fibrosis index, VEGF, NGF — with the canonical
  **PGE2 → cAMP → aromatase → local E2** positive feedback loop.
* **Clinical readouts (6 compartments)** — JZ thickness, PBAC HMB score,
  dysmenorrhea VAS, EHP-30 QoL, hemoglobin, lumbar BMD.

### 10 prebuilt scenarios

| # | Scenario | Calibration anchor |
|---|---|---|
| 1 | Untreated natural history | Vannuccini 2017 |
| 2 | NSAID + TXA symptomatic   | Lukes 2010 / NICE NG88 |
| 3 | Continuous COC            | Vannuccini 2018 |
| 4 | Dienogest 2 mg PO QD      | Osuga 2017 KIM-ADENO |
| 5 | LNG-IUS 52 mg             | Sheng 2009 / Cho 2008 |
| 6 | Leuprolide depot (no add-back) | Hornstein 1998 |
| 7 | Leuprolide + add-back     | Pierce 1999 / Surrey 2002 |
| 8 | Relugolix-CT 40 mg QD     | Giudice 2022 SPIRIT-1/2 |
| 9 | Elagolix 200 mg PO BID    | Taylor 2017 ELARIS-EM |
| 10 | Letrozole 2.5 mg QD       | Kim 2013 / Badawy 2012 |

## 4. Shiny dashboard (8 tabs)

1. **Patient profile** — adjustable age/BMI/parity/severity & baseline biomarkers.
2. **Drug PK** — log-scale plasma concentrations for the eight tracked drugs.
3. **HPO axis PD** — FSH/LH and systemic E2/P4 trajectories.
4. **Lesion biology** — aromatase, COX-2/PGE2, TGF-β, local E2, VEGF, NGF, fibrosis.
5. **Pain endpoint** — VAS and EHP-30 QoL.
6. **HMB/bleeding** — PBAC score and hemoglobin trajectory.
7. **Scenario comparison** — runs all 10 scenarios with the chosen profile.
8. **Biomarkers / safety** — BMD, JZ thickness, endpoint summary table.

## 5. How to run

```bash
# 1) Render the mechanistic map
dot -Tsvg adeno_qsp_model.dot -o adeno_qsp_model.svg
dot -Tpng -Gdpi=150 adeno_qsp_model.dot -o adeno_qsp_model.png
```

```r
# 2) Simulate scenarios in R
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny",
                   "shinydashboard","DT"))
library(mrgsolve)
mod <- mread_cache("adeno_mrgsolve_model.R")
out <- mrgsim(mod, ev(amt=2, cmt="GUT_DNG", ii=1, addl=365), end=365)
plot(out, c("VAS_pain","PBAC_score","JZ_mm","E2_pgmL"))

# 3) Launch the dashboard
shiny::runApp("adeno_shiny_app.R")
```

## 6. Key clinical anchors used during calibration

| Endpoint | Comparator | Expected effect |
|---|---|---|
| Dysmenorrhea VAS | Dienogest 2 mg vs placebo | −5.5 → −1.5 by 24 wk (Osuga 2017) |
| Dysmenorrhea VAS | Relugolix-CT vs placebo  | ≥ 30 % responder by 24 wk (SPIRIT) |
| Dysmenorrhea VAS | Elagolix 200 BID vs placebo | ≥ 50 % responder ~75 % (ELARIS) |
| PBAC             | LNG-IUS at 6 mo          | 270 → 55 (Sheng 2009) |
| Hb               | LNG-IUS at 12 mo         | +1.8 g/dL |
| JZ thickness     | Dienogest at 24 wk       | 11 → 8 mm |
| E2 systemic      | Leuprolide @ wk 4        | < 20 pg/mL |
| BMD loss         | Leuprolide 6 mo vs +add-back | −3.5 % → −0.5 to −1 % |
| Uterine vol      | Letrozole 12 wk          | −60 mL (Kim 2013) |

## 7. Caveats

* Designed for **research, education, and hypothesis generation** — not
  clinical decision support.
* All PK/PD parameters are population means; inter-individual variability
  needs to be added with `omega()` blocks for population simulations.
* The bidirectional PGE2 ↔ aromatase loop is intentionally tunable —
  adenomyosis biology is still actively being characterized.

## 8. License

Inherits the repository [LICENSE](../LICENSE).
