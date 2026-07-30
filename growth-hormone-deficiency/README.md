# 소아 성장호르몬 결핍증 (Paediatric Growth Hormone Deficiency, GHD) — QSP Model

> Quantitative Systems Pharmacology model of paediatric growth hormone
> deficiency. Linear growth is not a state to be prescribed — it is the
> **product** of a chain of independently breakable gains, and every lesion in
> the disease (and every drug) acts on exactly one link of that chain.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`ghd_qsp_model.dot`](ghd_qsp_model.dot) |
| 🖼️ Map (SVG, zoomable)   | [`ghd_qsp_model.svg`](ghd_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`ghd_qsp_model.png`](ghd_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`ghd_mrgsolve_model.R`](ghd_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`ghd_shiny_app.R`](ghd_shiny_app.R) |
| 📚 References (79, PubMed-verified) | [`ghd_references.md`](ghd_references.md) |

**Scale:** 197-node / 295-edge / 20-cluster mechanistic map · 42-ODE mrgsolve
model with 47 captured outputs · 19 prebuilt scenarios · 8-tab Shiny dashboard.

---

## 1. The organising idea

```
HV  =  [hypothalamic GHRH/SST/ghrelin drive]
       x [somatotroph secretory mass]
       x [GH pharmacokinetics and free fraction (GHBP)]
       x [GHR - JAK2 - STAT5b transduction, under a SOCS2 brake]
       x [hepatic IGF-1 output x its BINDING state (IGFBP-3/ALS ternary complex)]
       x [growth-plate responsiveness]
       x [REMAINING resting-zone proliferative reserve]
```

Writing growth as a product rather than a rate is what makes the model
diagnostic instead of merely descriptive. Three clinically important behaviours
then **fall out of the structure** rather than being coded:

**(1) Catch-up growth is a consumable, not a setpoint.** The resting-zone
reserve `RZ` multiplies the proliferation rate *and* is spent in proportion to
the growth actually realised. An untreated GHD child grows slowly and therefore
**saves** reserve; starting rhGH spends it. That single law reproduces, with no
fitted decay term: year-1 height velocity of ~10.5 cm/yr falling to ~8.0 cm/yr
by year 4 on an unchanged mg/kg dose, and the fact that **age at initiation is
the dominant modifiable predictor of near-adult height** (starting at age 10
instead of 5 ends 2.4 height-SDS lower in the model, on identical dosing).

**(2) Oestrogen enters twice, with opposite sign.** E2 amplifies GH pulse
amplitude — the pubertal spurt, cm/yr up — *and* accelerates reserve
consumption and epiphyseal fusion — growing years down. Whether a GnRH analogue
or an aromatase inhibitor helps is therefore an **emergent trade-off** in this
model, not a coded height gain. It also predicts the aromatase inhibitor's
bone-mineral penalty as part of the same trade-off (BMD Z −0.52 vs −0.06).

**(3) GHD is a double hit on IGF-1.** Hepatic IGF-1 synthesis falls *and* the
IGFBP-3/ALS carrier pool falls, so the free fraction rises and clearance of what
little is made accelerates by ~40%. The model therefore reproduces the awkward
clinical fact that **free IGF-1 is much better preserved than total IGF-1** in
GHD — which is precisely why growth in this model is driven mainly by the
**local, GH-dependent** growth-plate IGF-1 arm rather than by circulating free
IGF-1. Without that structure, a Laron-syndrome patient (high GH, no receptor)
would grow normally.

**Deficiency vs insensitivity.** The direct (dual-effector) GH action on the
resting zone is gated by receptor density **and** by STAT5b capacity. That one
gate is what separates GH *deficiency* (low GH, intact receptor — treat with
rhGH) from GH *insensitivity* (high GH, broken receptor — rhGH buys almost
nothing, and mecasermin is the only lever). Both phenotypes come out of the same
42 equations with two parameters changed.

---

## 2. Mechanistic map — 20 clusters, 197 nodes, 295 edges

Open [`ghd_qsp_model.svg`](ghd_qsp_model.svg) to zoom.

1. **유전·병인** — GH1, GHRHR, POU1F1, PROP1, HESX1, LHX3/4, SOX2/3, OTX2/GLI2/ARNT2, pituitary stalk interruption, craniopharyngioma, cranial irradiation, TBI/hypophysitis/haemochromatosis, idiopathic; and the GH-insensitivity genes (GHR, STAT5B, IGFALS, IGF1, PAPPA2)
2. **시상하부 오실레이터** — arcuate GHRH, periventricular somatostatin, reciprocal phase coupling, ghrelin/GHS-R1a, slow-wave sleep, exercise, hypoglycaemia, FFA inhibition, adiposity, leptin/NPY
3. **뇌하수체 somatotroph** — GHRH-R/Gs/cAMP/PKA/CREB/Pit-1, SSTR2/5, Ca²⁺ influx, granule pool, functional mass, pulse amplitude vs frequency, 22/20-kDa isoforms
4. **GH 순환·소실** — GHBP (shed GHR ectodomain), free GH, peripheral distribution, receptor-mediated vs linear clearance, flip-flop kinetics
5. **GHR–JAK2–STAT5b** — dimerisation, JAK2, STAT5b pY699, STAT1/3, SHC-ERK, IRS-1-PI3K-AKT, GAS elements, SOCS2/CIS brake, SHP-1/PTP1B, ligand-induced down-regulation, d3-GHR polymorphism
6. **IGF 시스템** — hepatic IGF-1, IGFBP-3, ALS, the 150-kDa ternary complex, 50-kDa binary complexes, IGFBP-1 (insulin-suppressed), PAPP-A2, free/bioactive IGF-1, IGF-1R/IRS-1, IGF-2R clearance, local growth-plate IGF-1
7. **성장판 엔진** — resting/proliferative/prehypertrophic/hypertrophic zones, IHH–PTHrP loop, FGF18–FGFR3 brake, CNP–NPR2, SOX9, RUNX2/MEF2C, COL10A1/MMP13/VEGF, vascular invasion, senescence, fusion, catch-up capacity
8. **뼈·미네랄** — osteoblast/osteoclast, RANKL/OPG, P1NP/osteocalcin/ALP, remodelling transient, BMD Z, cortical thickness, peak bone mass, SCFE, scoliosis
9. **사춘기 축** — KNDy kisspeptin/NKB/dynorphin, MKRN3 brake, GnRH pulse generator, LH/FSH, gonad, testosterone, aromatase, oestradiol, ERα, the spurt, fusion, hypogonadism in CPHD
10. **체성분·중간대사** — HSL/ATGL lipolysis, FFA flux, fat and lean mass, hepatic gluconeogenesis, insulin resistance, fasting insulin/glucose, infantile hypoglycaemia, lipids, REE, Na⁺/water retention
11. **다른 축과의 상호작용** — type-2 deiodinase and the free-T4 fall, 11β-HSD1 and unmasked adrenal insufficiency, glucocorticoid as a potent growth suppressor, prolactin/ADH
12. **일일 소마트로핀 PK** — SC depot, ka, CL/F, flip-flop, the single supraphysiological daily peak, adherence decay, devices
13. **주 1회 지속형 GH** — lonapegsomatropin (TransCon carrier), somatrogon (CTP fusion), somapacitan (albumin binder), release kinetics, non-physiological peak:trough, the within-week IGF-1 swing and why the sampling day matters
14. **메카세르민 (rhIGF-1)** — dosing, hypoglycaemia, lymphoid overgrowth, why it under-performs rhGH in GHD
15. **보조 약물** — GnRH analogue, aromatase inhibitor, oxandrolone, macimorelin (diagnostic), sex-steroid priming, CNP analogue, glucocorticoid minimisation
16. **진단 캐스케이드** — height velocity, height SDS, IGF-1/IGFBP-3 SDS, provocation tests, the assay-dependent cut-off, MRI, bone age, gene panel, IGF-1 generation test, retesting
17. **임상 엔드포인트** — annualised height velocity, height SDS, ΔHtSDS year 1, near-adult height, mid-parental target, predicted-adult-height gain, bone-age/chronological-age ratio, HRQoL
18. **안전성** — the IGF-1 SDS > +2 ceiling, benign intracranial hypertension, insulin resistance/T2DM, anti-drug antibodies, injection-site reactions, OSA in Prader-Willi, neoplasia surveillance, SAGhE/KIGS cohorts
19. **반응 예측 인자** — age at start, dose, provoked GH peak, height deficit vs target, birth weight/SGA, weight SDS, d3-GHR, first-year response, KIGS/Ranke models
20. **전환기·성인 GHD** — stopping criteria, retesting after washout, adult dosing, adult phenotype, cardio-metabolic legacy, structured transition

---

## 3. mrgsolve model — 42 ODE compartments

Time unit **days**; 4-year default horizon for growth endpoints, hours-scale
resolution for PK and provocation testing.

| Block | Compartments |
|---|---|
| Hypothalamic oscillator | `GHRH` `SST` `GHREL` |
| Somatotroph | `PITM` (functional mass) `GHPOOL` (releasable granules) |
| GH disposition | `GHC` `GHP` `GHBP` |
| Drug kinetics | `DEPD` (daily SC) · `DEPW` → `PROD` (weekly carrier → bioactive GH) · `DEPM` (mecasermin) · `LEUP` (GnRH analogue) · `AIC` (aromatase inhibitor) |
| Transduction | `GHRc` `PS5` (STAT5b signal) `SOCS` |
| IGF system | `IGF1` `BP3` `ALSC` `BP1` `IGFL` (local) |
| Growth plate | `RZ` (reserve) `PZ` `HZ` `GPS` (fusion program) |
| Skeleton | `HT` `BA` `BMDZ` |
| Puberty | `GT` `TST` `E2` |
| Body composition & metabolism | `FM` `LBM` `IRES` `INS` `GLU` `FFA` `LDL` |
| Safety | `ADA` `FT4` `ICPX` |

47 captured outputs including `HV_cm_yr`, `Height_SDS`, `IGF1_SDS`,
`Reserve_RZ`, `Fusion_prog`, `BA_CA_ratio`, `Ternary_frac`, `HOMA_IR`,
`Thyroid_tone`, `ICP_index`.

### 19 prebuilt scenarios

Natural history · healthy reference child · daily somatropin at 0.025 / 0.034 /
0.050 mg/kg/day · lonapegsomatropin / somatrogon / somapacitan at label dose ·
late start at age 10 · 60% adherence · GH insensitivity untreated and on
mecasermin · neutralising anti-drug antibodies · pubertal start alone, with a
GnRH analogue, and with an aromatase inhibitor · unreplaced central
hypothyroidism · high-dose glucocorticoid · IGFALS deficiency · plus a
closed-loop `GHD_titrate()` that reviews IGF-1 SDS annually and moves the dose
toward a target band.

---

## 4. Calibration — what the model actually produces

Read off this model as committed (mrgsolve 2.0.1 / R 4.3.3), default patient a
5-year-old boy with severe congenital isolated GHD, height SDS −2.3, 15.5 kg,
dose re-titrated +10%/yr for growth.

| Anchor | Target | Model |
|---|---|---|
| Healthy child: GH secretion rate | ~300 µg/day | 397 µg/day |
| Healthy child: IGF-1 / IGFBP-3 (age 5) | ~90 ng/mL / 2–4 mg/L | 97 / 3.3 |
| Healthy child: HV age 5→6→7→8 | 6.0 5.7 5.6 5.3 | 6.6 6.0 5.5 5.0 |
| Healthy child: height SDS drift over 4 yr | ~0 | +0.19 |
| Provoked peak GH: severe / partial / healthy | <5 / 5–9 / >10 µg/L | 0.4 / 2.2 / 10.1 |
| Untreated severe GHD: HV, IGF-1 SDS | 3–4 cm/yr, −2 to −3 | 3.3 cm/yr, −1.9 |
| Untreated severe GHD: height SDS lost/yr | 0.3–0.4 | 0.42 |
| Daily somatropin single dose: Cmax / t½ | 15–40 µg/L / 2–4 h | 17.3 µg/L / 4.4 h |
| Year-1 HV: 0.025 / 0.034 / 0.050 mg/kg/day | 9.3 / 10.3–10.8 / 11.3 | 9.6 / 10.5 / 11.5 |
| IGF-1 SDS on 0.034 / 0.050 mg/kg/day | ~+0.5 / ~+1.5 | +0.4 / +1.4 |
| ΔHeight SDS year 1 (0.034) | +0.8 to +1.1 | +1.06 |
| HV year 1→2→3→4 on unchanged mg/kg | 10.5 8.5 7.5 — | 10.5 9.5 8.7 8.0 |
| Lonapegsomatropin 0.24 mg/kg/wk (heiGHt) | 11.2 cm/yr | 10.3 |
| Somatrogon 0.66 mg/kg/wk | 10.1 cm/yr | 9.8 |
| Somapacitan 0.16 mg/kg/wk (REAL4) | 11.2 cm/yr | 10.2 |
| Weekly IGF-1 SDS peak day / trough day | 2–4 / 7 | 2.6–3.1 / 7 |
| Within-week IGF-1 SDS swing | 1.5–2.5 SDS | 1.7–2.6 by product |
| GH insensitivity untreated / + mecasermin | 3–4 / 8–9 cm/yr | 3.5 / 7.6 |
| IGFALS deficiency: IGF-1 SDS / HV | very low / mild deficit | −2.5 / 4.3 cm/yr |
| Unreplaced central hypothyroidism | blunted | 2.4 cm/yr, free T4 0.52, BA/CA 0.63 |
| High-dose glucocorticoid (12 mg/m²/day) | blunted | 4.0 cm/yr |
| 60% adherence | costs 1–2 cm/yr | 10.5 → 8.3 |
| Neutralising ADA | growth attenuation | 10.5 → 8.2 |
| Late start (age 10, same dose) | worse outcome | final height SDS −1.80 vs +0.59 |
| GH + GnRH analogue from age 12 | +4–5 cm | +1.14 height SDS |
| GH + anastrozole from age 12 | smaller gain, BMD cost | +0.84 SDS, BMD Z −0.52 vs −0.06 |

Building and validating this model changed it eight times. The defects the
validation runs exposed, and how each was fixed, are recorded in the model file:
a bone-mineral equation that integrated without a homeostatic term and drifted
to Z = +60; an IGF-1 distribution volume seven times too small, which made
mecasermin produce IGF-1 of 3600 ng/mL; a free-T4 reference that made every
euthyroid child mildly thyrotoxic; neutralising antibodies written into the free
fraction, where — because receptor-mediated clearance also acts on free GH —
they **cancelled exactly** and had no effect at steady state; a direct
growth-plate GH term that was not gated by the receptor, so Laron syndrome grew
normally; a severity index that read only the pituitary, so Laron and IGFALS
patients started at normal height; `$MAIN` re-deriving initial conditions and
silently defeating state carry-forward in the annual titration loop (now guarded
by `INITMODE`); and a senescence law keyed to the growth *drive* rather than to
realised growth, which could not separate a healthy child from a child in
catch-up because their drives differ by only ~7% while their velocities differ
by 60%.

---

## 5. Known limitations — read before using any output

1. **Ultradian GH pulses are not resolved.** Only a circadian (nocturnal-surge)
   modulation is carried, so the simulated endogenous profile is smoother than
   physiology: peaks underestimated (~11 vs 15–25 µg/L), troughs overestimated
   (~1.9 vs <0.2 µg/L), and the 24-h mean therefore sits above the reported
   2–3 µg/L. Any question that turns on pulse *shape* rather than 24-h exposure
   is out of scope.
2. The GH → STAT5b step uses an **apparent, system-level EC50 of 60 µg/L**,
   above the molecular GHR Kd, so that the IGF-1 dose–response stays close to
   proportional across the therapeutic range as clinical data require. It is a
   lumped transduction gain, not a binding constant.
3. Height and IGF-1 SDS use **smooth internal reference tables**, not a
   validated national growth chart or an assay-specific IGF-1 reference. SDS
   values are internally consistent, not transferable.
4. Provoked GH peaks are compressed and saturate in the secretagogue dose, so
   the model separates severe / partial / normal secretion but must not be used
   to explore diagnostic cut-offs or to compare stimulation agents.
5. Body composition, insulin resistance, lipids, BMD Z, the intracranial-pressure
   index and the ADA titre are **semi-mechanistic target-tracking**
   representations: directions and time courses are calibrated, absolute values
   are illustrative. `ICP_index` is a relative risk index, not a pressure.
6. Adherence is deterministic (every n-th dose taken), not a stochastic
   missed-dose process.
7. The initial reserve `RZ_0` integrates the same consumption law over the
   child's life so far but using a **lifetime-average** growth rate, so starting
   a run at age 8 gives ~15–20% less reserve than integrating the same child
   from age 5 to 8. Compare arms at the same starting age.
8. `POT_W` (bioactive GH equivalents per mg of construct) is **fitted** per
   product so that each label dose reproduces its own phase-3 height velocity.
   It is not an independently measured potency.

---

## 6. How to run

```r
install.packages(c("mrgsolve", "dplyr", "tidyr", "ggplot2", "shiny", "DT"))

library(mrgsolve); library(dplyr)
mod <- mread_cache("ghd_mrgsolve_model.R")

# untreated natural history vs daily somatropin 0.034 mg/kg/day, 4 years
none  <- ev(amt = 0, cmt = "DEPD")
daily <- ev(time = 0:1459, amt = 0.034 * 15.5, cmt = "DEPD")
mrgsim(mod, daily, end = 1460, delta = 1, hmax = 0.25) %>% as_tibble()

# GH provocation test (bolus into the ghrelin/secretagogue compartment)
mrgsim(mod, ev(amt = 25, cmt = "GHREL"), end = 0.25, delta = 0.001)

# the interactive dashboard
shiny::runApp("ghd_shiny_app.R")
```

Use `hmax = 0.25` on multi-year runs so the 24-h circadian term is not stepped
over; set `AMPC = 0` to switch circadian modulation off for faster screening.
The commented helper block at the end of the model file contains all 19
scenarios and the annual IGF-1-titration loop.

To re-render the map:

```bash
dot -Tsvg ghd_qsp_model.dot -o ghd_qsp_model.svg
dot -Tpng -Gdpi=150 -Gsize=60,27 ghd_qsp_model.dot -o ghd_qsp_model.png
```

---

## ⚠️ 면책 조항 (Disclaimer)

교육·연구 목적의 QSP 모델입니다. 공개 문헌을 바탕으로 구성했으나 독립적으로
검증·인증되지 않았으며 **실제 임상 의사결정, 처방, 규제 제출에 사용해서는 안
됩니다.** 파라미터는 설명을 위한 근사치이며 개별 환자 데이터에 적합되지
않았습니다. 특히 이 모델이 출력하는 신장 SDS·IGF-1 SDS 값은 내부 참조표에
기반한 것으로, 검증된 성장도표나 검사법별 IGF-1 참조범위를 대체하지 않습니다.

This model is for education, research and hypothesis generation only. It has not
been independently validated or certified and must not be used for clinical
decisions, prescribing, or regulatory submission.
