# Vasomotor Symptoms of Menopause (VMS, 폐경기 혈관운동증상) — QSP Model

> Integrated Quantitative Systems Pharmacology model of menopausal hot
> flushes and night sweats, linking ovarian follicle depletion and loss of
> estrogen negative feedback to arcuate KNDy-neuron (kisspeptin/neurokinin
> B/dynorphin) hypertrophy, NKB→NK3R-driven narrowing of the hypothalamic
> thermoneutral zone, and the resulting hot-flush frequency/severity —
> together with the modern pharmacology stack (NK3R antagonist
> fezolinetant, dual NK1/NK3 antagonist elinzanetant, hormone therapy,
> SSRIs/SNRIs, gabapentin, clonidine, oxybutynin) and downstream organ-
> system consequences (sleep, mood, bone, lipids, genitourinary syndrome).

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`vms_qsp_model.dot`](vms_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`vms_qsp_model.svg`](vms_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`vms_qsp_model.png`](vms_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`vms_mrgsolve_model.R`](vms_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`vms_shiny_app.R`](vms_shiny_app.R) |
| 📚 References             | [`vms_references.md`](vms_references.md) |

---

## 1. Disease in one paragraph

Vasomotor symptoms (VMS) — hot flushes and night sweats — affect up to
80% of women during the menopause transition and persist a median of
**7.4 years** (SWAN). Declining ovarian follicle numbers reduce estradiol
(E2) and inhibin B, removing negative feedback on the hypothalamus. Arcuate
(infundibular) **KNDy neurons**, which co-express kisspeptin, neurokinin B
(NKB), and dynorphin, undergo estrogen-withdrawal-driven **hypertrophy**
and increased NKB output. NKB acts on **NK3R** both autosynaptically (a
feed-forward loop) and on neurons of the **median preoptic nucleus (MnPO)**
thermoregulatory center, narrowing the **thermoneutral zone** so that small
elevations in core temperature trigger inappropriate heat-dissipation
responses (cutaneous vasodilation, sweating) — the hot flush. Downstream
consequences include sleep fragmentation, mood disturbance, accelerated
bone loss, adverse lipid changes, and genitourinary syndrome of menopause.
The 2023 approval of **fezolinetant**, a selective NK3R antagonist, and the
late-stage dual NK1/NK3 antagonist **elinzanetant** directly target this
mechanism, complementing hormone therapy and older nonhormonal options
(SSRIs/SNRIs, gabapentin, clonidine, oxybutynin).

## 2. Mechanistic clusters (14 in the DOT map, 100+ nodes)

1. Ovarian aging & follicle depletion (STRAW+10 staging)
2. HPG-axis feedback loss (GnRH pulse generator → FSH/LH → E2/E1)
3. Arcuate KNDy-neuron biology (Kiss1 · NKB/TAC3 · dynorphin · ERα brake)
4. NKB→NK3R signaling & median preoptic (MnPO) thermoregulatory center
5. Sympathetic effector arm (cutaneous vasodilation, sweating, flushing)
6. Clinical VMS phenotype (hot flush trigger → frequency → severity → night sweats)
7. Sleep architecture & CNS/mood consequences
8. Skeletal consequences (RANKL/OPG · osteoclast · CTX/P1NP · BMD · fracture)
9. Cardiometabolic & lipid consequences (LDL/HDL, endothelial function, adiposity)
10. Genitourinary syndrome of menopause (vaginal atrophy, pH, dyspareunia)
11. Drug PK/PD — NK3R antagonist (fezolinetant) & dual NK1/3 (elinzanetant)
12. Drug PK/PD — hormone therapy (oral/transdermal estradiol, CEE, progestogen)
13. Drug PK/PD — nonhormonal agents (SSRI/SNRI, gabapentinoid, α2-agonist)
14. Clinical endpoints, biomarkers & PROs (HF frequency/severity, MENQOL, PSQI, DXA)

## 3. mrgsolve model (24 ODE compartments)

* **Drug PK (5 drug classes, 17 compartments)** — fezolinetant, elinzanetant,
  oral estradiol, transdermal E2 patch, paroxetine, venlafaxine + active
  metabolite ODV, gabapentin (saturable absorption), clonidine (literature
  PK anchors in `vms_references.md`).
* **Ovarian/HPG axis (3 compartments)** — follicle pool, plasma E2, FSH —
  driving the estrogen-withdrawal signal.
* **KNDy/thermoregulatory network (2 compartments)** — KNDy/NKB tone and
  thermoneutral-zone (TNZ) half-width, with NK3R-antagonist and
  nonhormonal-drug reversal terms feeding directly into TNZ narrowing.
* **Clinical readouts (6 compartments)** — hot-flush frequency, hot-flush
  severity, PSQI sleep score, mood composite, bone resorption (CTX) and
  BMD, LDL-C.

### 10 prebuilt scenarios

| # | Scenario | Calibration anchor |
|---|---|---|
| 1 | Untreated natural history        | Avis 2015 SWAN |
| 2 | Fezolinetant 45 mg PO QD         | Lederman 2023 SKYLIGHT 1/2 |
| 3 | Fezolinetant 30 mg PO QD         | Johnson 2023 Ph2b |
| 4 | Elinzanetant 120 mg PO QD        | Pinkerton 2024 OASIS 1/2 |
| 5 | Oral estradiol 1 mg QD           | Kuhl 2005 / PEPI 1995 |
| 6 | Transdermal E2 patch 0.05 mg/d   | Notelovitz 2000 |
| 7 | Paroxetine 7.5 mg QD             | Simon 2013 / Stearns 2003 |
| 8 | Venlafaxine ER 75 mg QD          | Loprinzi 2000 / Joffe 2014 MsFLASH |
| 9 | Gabapentin 900 mg/d (300 mg TID) | Guttuso 2003 |
| 10 | Clonidine 0.1 mg QD              | Boekhout 2011 |

## 4. Shiny dashboard (8 tabs)

1. **Patient profile** — adjustable age/BMI/menopausal stage/severity & baselines.
2. **Drug PK** — log-scale plasma concentrations for the nine tracked drug moieties.
3. **HPG / KNDy axis** — follicle pool, E2, FSH, and KNDy/NKB tone trajectories.
4. **Thermoregulation** — thermoneutral-zone (TNZ) half-width vs. premenopausal reference.
5. **Hot-flush endpoint** — frequency (events/day) and severity (0-3).
6. **Sleep / Mood** — PSQI sleep-quality index and mood composite.
7. **Scenario comparison** — runs all 10 scenarios with the chosen profile.
8. **Biomarkers / safety** — bone (CTX, BMD), LDL-C, endpoint summary table.

## 5. How to run

```bash
# 1) Render the mechanistic map
dot -Tsvg vms_qsp_model.dot -o vms_qsp_model.svg
dot -Tpng -Gdpi=150 vms_qsp_model.dot -o vms_qsp_model.png
```

```r
# 2) Simulate scenarios in R
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny",
                   "shinydashboard","DT"))
library(mrgsolve)
mod <- mread_cache("vms_mrgsolve_model.R")
out <- mrgsim(mod, ev(amt=45, cmt="GUT_FEZ", ii=1, addl=168), end=168)
plot(out, c("HF_freq_perday","TNZ_degC","E2_pgmL","PSQI_sleep"))

# 3) Launch the dashboard
shiny::runApp("vms_shiny_app.R")
```

## 6. Key clinical anchors used during calibration

| Endpoint | Comparator | Expected effect |
|---|---|---|
| HF frequency | Fezolinetant 45 mg vs placebo | −63 to −64% at wk 12 (SKYLIGHT 1/2) |
| HF frequency | Fezolinetant 30 mg vs placebo | −60% at wk 12 |
| HF frequency | Elinzanetant 120 mg vs placebo | −74 to −84% at wk 12 (OASIS 1/2) |
| HF frequency | Oral/transdermal E2 vs placebo | −75 to −90% (gold-standard efficacy) |
| HF frequency | Paroxetine 7.5 mg vs placebo | −55 to −65% |
| HF frequency | Venlafaxine 75 mg vs placebo | −55% vs −15% placebo |
| HF frequency | Gabapentin 900 mg/d vs placebo | −45% vs −29% placebo; improves PSQI |
| HF frequency | Clonidine 0.1 mg vs placebo | −25 to −35%, modest |
| BMD          | HRT vs no HRT, 3 y | prevents ~2%/yr peri-menopausal loss (PEPI) |
| LDL-C        | Oral CEE 0.625 mg, 1 y | −10 to −15 mg/dL (PEPI) |

## 7. Caveats

* Designed for **research, education, and hypothesis generation** — not
  clinical decision support.
* All PK/PD parameters are population means; inter-individual variability
  needs to be added with `omega()` blocks for population simulations.
* The hot-flush "frequency" compartment is a continuous surrogate for what
  is physiologically a discrete, Poisson-like event process; it captures
  trends in mean daily rate, not individual event timing.
* Hormone-therapy safety outcomes (VTE, breast cancer, stroke risk from
  WHI) are referenced in the bibliography but intentionally **not**
  modeled as ODE compartments — this model focuses on VMS efficacy and the
  bone/lipid/GSM mechanistic axis.

## 8. License

Inherits the repository [LICENSE](../LICENSE).
