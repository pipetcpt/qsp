# Diabetic Foot Ulcer (DFU) — QSP Model

> An integrated Quantitative Systems Pharmacology model of the diabetic foot
> ulcer, linking the three axes that make a foot ulcerate — **peripheral
> neuropathy with loss of protective sensation**, **peripheral arterial disease
> and microvascular dysfunction**, and **repetitive plantar pressure on an
> insensate, deformed foot** — to the wound biology that then refuses to close:
> a protease-dominant, biofilm-colonised, senescent-fibroblast bed in which the
> M1→M2 macrophage switch never happens and keratinocytes proliferate but do not
> migrate. The model carries the full modern therapeutic stack (offloading,
> sharp debridement, becaplermin, TLC-NOSF sucrose octasulfate, esmolol gel,
> cellular tissue products, NPWT, hyperbaric/topical oxygen, systemic
> antibiotics with explicit bone penetration, revascularisation, glucose-lowering
> therapy) and — unusually — it keeps running **after** the wound closes, because
> in this disease closure is a remission and the recurrence is the illness.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`dfu_qsp_model.dot`](dfu_qsp_model.dot) |
| 🖼️ Map (SVG) | [`dfu_qsp_model.svg`](dfu_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi) | [`dfu_qsp_model.png`](dfu_qsp_model.png) |
| ⚙️ mrgsolve ODE model (44 compartments) | [`dfu_mrgsolve_model.R`](dfu_mrgsolve_model.R) |
| 🐍 Dependency-free reference implementation | [`dfu_reference_implementation.py`](dfu_reference_implementation.py) |
| 📊 Shiny dashboard (12 tabs) | [`dfu_shiny_app.R`](dfu_shiny_app.R) |
| 📚 References (95 entries) | [`dfu_references.md`](dfu_references.md) |

**Every number in this README is computed, not asserted.** Reproduce them with:

```bash
python3 dfu_reference_implementation.py            # the 23-scenario table + analysis
python3 dfu_reference_implementation.py --anchors  # model vs published anchors
```

That script integrates the identical 44-state ODE system with fixed-step RK4
using nothing but the Python standard library, so the claims below can be
checked without R, without mrgsolve, and without installing anything.

---

## 1. The disease in one paragraph

A diabetic foot ulcer is not primarily a wound-healing disease; it is a
**sensory-motor-vascular-mechanical** disease that happens to present as a
wound. Chronic hyperglycaemia drives four parallel biochemical pathways (polyol,
hexosamine, PKC-β, advanced glycation) that destroy small and large peripheral
nerve fibres. Small-fibre loss removes pain — the warning signal — while motor
neuropathy produces claw-toe deformity and autonomic neuropathy produces dry,
fissured skin and arteriovenous shunting. Glycated collagen stiffens the joints
and shortens the Achilles tendon, driving peak plantar pressure up at exactly
the metatarsal heads that the atrophied fat pad no longer protects. Callus forms
over the pressure point, a subkeratotic haematoma develops, and the skin breaks
— painlessly, so the patient keeps walking on it. Once open, the wound enters a
chronic inflammatory state it cannot leave: neutrophils persist and undergo
NETosis (primed by hyperglycaemia), the M1→M2 macrophage switch fails, IL-1β
sustains its own production, MMP-8 and MMP-9 overwhelm TIMP-1 and digest both
the matrix and any growth factor applied to it, fibroblasts become senescent,
and keratinocytes at the wound edge proliferate vigorously while failing to
migrate. Layered on top is infection: >60% of chronic wounds carry a biofilm
that returns to full antibiotic tolerance within 2–3 days of debridement, and
~20% of infected ulcers reach bone. Roughly 19–34% of people with diabetes will
develop a foot ulcer; ~20% of those with an infected ulcer will lose part of the
limb; and among those who heal, **~40% re-ulcerate within a year and ~60% within
three** — because closing the wound changes nothing about the neuropathy, the
deformity or the pressure that produced it.

---

## 2. Mechanistic map — 290 nodes, 21 clusters, 421 edges

[![DFU mechanistic map](dfu_qsp_model.png)](dfu_qsp_model.svg)

| # | Cluster | What it contains |
|---|---|---|
| 1 | Systemic diabetic milieu | HbA1c, polyol/hexosamine/PKC-β/AGE-RAGE, dyslipidaemia, CKD, nutrition, anaemia, smoking |
| 2 | Peripheral neuropathy | small/large-fibre loss, IENFD, monofilament, VPT, LOPS, motor & autonomic neuropathy, loss of neurogenic inflammation |
| 3 | Biomechanics | limited joint mobility, equinus, claw toe, MT-head prominence, Charcot, peak plantar pressure, callus, repetitive stress |
| 4 | PAD & macrovascular perfusion | infrapopliteal disease, Mönckeberg calcification, ABI/toe pressure/TcPO₂, WIfI, revascularisation, restenosis |
| 5 | Microvascular dysfunction | capillary BM thickening, eNOS uncoupling, ET-1, pericyte loss, functional ischaemia with a normal ABI |
| 6 | Oxygenation & HIF | wound pO₂, prolyl hydroxylase, HBOT/topical O₂, and the methylglyoxal block of HIF-1α transactivation |
| 7 | Ulcer initiation & classification | area, perimeter, depth, probe-to-bone, Wagner/UT/SINBAD/IWGDF |
| 8 | Haemostasis & acute onset | platelet PDGF/TGF-β release, DAMPs, TLR4, complement, early recruitment |
| 9 | Chronic stalled inflammation | NETosis, M1 persistence, the failed M1→M2 switch, NLRP3/IL-1β autoloop, FOXO1, miR dysregulation |
| 10 | Protease imbalance | MMP-1/2/8/9, TIMP-1/2, the MMP-9/TIMP-1 predictor, proteolysis of applied growth factor |
| 11 | Growth-factor axis | PDGF-BB, VEGF, FGF-2, EGF, KGF, TGF-β1/3, SDF-1α and its DPP-4 cleavage |
| 12 | Angiogenesis | tip/stalk cells, VEGFR2, Ang/Tie2, EPC mobilopathy, granulation vascular density |
| 13 | Fibroblasts & ECM | proliferation, senescence + SASP, myofibroblast contraction, collagen I/III, LOX cross-linking |
| 14 | Keratinocytes | the migration paradox (c-Myc/β-catenin), epibole, edge callus, epithelial advance rate |
| 15 | Infection & biofilm | bioburden, EPS matrix, quorum sensing, persisters, tolerance, osteomyelitis, bone penetration |
| 16 | Impaired innate defence | chemotaxis/phagocytosis/killing, glycated complement, defensins — the inflamed-but-defenceless paradox |
| 17 | Offloading | TCC, instant TCC, removable walker, half-shoe, **adherence**, felted foam, Achilles lengthening |
| 18 | Pharmacology & procedures | debridement, becaplermin, TLC-NOSF, esmolol, NPWT, CTPs, LeucoPatch, antibiotics, rifampicin, GLP-1/SGLT2/DPP-4i, amputation |
| 19 | Clinical endpoints | PAR₄, closure, time-to-heal, recurrence, amputation-free survival, hospitalisation, QoL, 5-year mortality |
| 20 | Remission substrate | scar tensile strength, unchanged neuropathy and deformity, temperature monitoring, surveillance |
| 21 | Safety | becaplermin boxed warning, nephrotoxicity, *C. difficile*, resistance, hypoglycaemia, TCC-induced lesions, HBOT barotrauma |

---

## 3. ODE model — 44 compartments

`dfu_mrgsolve_model.R` (mrgsolve ≥ 1.5) and `dfu_reference_implementation.py`
implement the same system, in days:

| Block | States |
|---|---|
| Therapy PK (11) | becaplermin depot & wound tissue, TLC-NOSF, esmolol depot & wound, antibiotic central/wound/bone, cellular tissue product, oxygen, debridement action |
| Slow systemic axes (3) | HbA1c, neuropathy, effective perfusion |
| Infection (4) | planktonic bioburden (log₁₀ CFU/g), biofilm biomass, osteomyelitis, exposed-bone nidus |
| Inflammation (5) | neutrophils, M1, M2, IL-1β, ROS |
| Proteases (2) | MMP-9, TIMP-1 |
| Growth factors / angiogenesis (5) | PDGF-BB, VEGF-A, SDF-1α, EPC drive, granulation vascular density |
| Repair tissue (5) | fibroblasts, senescent fibroblasts, collagen, granulation fill, keratinocyte migratory competence |
| Wound (3) | **AREA**, depth, scar maturation |
| Outcomes (6) | re-ulceration hazard, amputation hazard, ulcer-free days, antibiotic-days, treatment burden, cumulative edge advance |

Twenty-three prebuilt scenarios cover offloading devices and adherence,
debridement frequency, four topical/device therapies alone and combined,
ischaemia and revascularisation timing, infection with and without debridement,
osteomyelitis managed medically vs surgically, HBOT, NPWT, cellular tissue
products, glycaemic intensification, and an optimal bundle with and without
post-closure remission care.

---

## 4. What the model actually says

### 4.1 Closure is a perimeter process — and that is why PAR₄ works

The wound closes by epithelium advancing inwards from the margin, so

```
dA/dt = −k(t) · P ,      P = 2√(πA)      ⟺      dr/dt = −k(t)
```

The **radius** falls linearly; all biology enters through `k` (cm/day) and
nothing enters through `A`. Two consequences fall out.

First, `k` is a property of the wound bed, not of the wound size. Across a
16-fold range of presenting area under identical therapy the model's mid-course
`k` varies by only 60% — i.e. across a 4-fold range of *radius* — while
time-to-heal varies 2.6-fold, almost entirely as geometry `t = r₀/k`:

| A₀ (cm²) | r₀ (cm) | PAR₄ | k mid-course (cm/day) | t_close observed | predicted from k | predicted from PAR₄ |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.399 | 57.6% | 0.00832 | 41.3 d | 58.1 d | 80.3 d |
| 1.0 | 0.564 | 47.0% | 0.01042 | 52.9 d | 66.3 d | 103.1 d |
| 2.0 | 0.798 | 39.3% | 0.01246 | 67.5 d | 77.7 d | 126.8 d |
| 4.0 | 1.128 | 34.3% | 0.01409 | 85.6 d | 93.5 d | 147.8 d |
| 8.0 | 1.596 | 31.7% | 0.01555 | 107.7 d | 113.2 d | 161.3 d |

Second, given a 4-week percent-area-reduction the perimeter law has an exact,
**size-free** closed form:

```
t_heal = 28 / (1 − √(1 − PAR₄))
```

PAR₄ itself is size-dependent (a big wound cannot shed 50% of its area in four
weeks at the same `k`), but the *mapping from PAR₄ to time-to-heal is not* —
which is precisely why Sheehan's 50% threshold works as a universal bedside
rule rather than needing a nomogram. PAR₄ = 50% forecasts closure at 95.6 days,
just past the 12-week trial window; that is the arithmetic behind the rule.

In the model this closed form systematically **over-predicts** (right-hand
column above vs the observed times), because
`k` is still rising during weeks 1–4 as the wound bed converts. PAR₄ is
therefore a *conservative* screen: wounds that pass it beat its own forecast,
and the ones that fail it are the ones whose `k` never rose.

### 4.2 Offloading is device efficacy × adherence — and adherence is the whole story

| Arm | Device efficacy | Adherence | Effective | t_close |
|---|---:|---:|---:|---:|
| No offloading | 0.00 | — | 0.00 | **> 540 d (never closes)** |
| Half-shoe | 0.40 | 0.55 | 0.22 | 72.7 d |
| Removable cast walker | **0.87** | **0.28** | 0.24 | 67.5 d |
| Total contact cast | 0.85 | 1.00 | 0.85 | 41.5 d |
| **The same removable walker at 100% wear** | **0.87** | **1.00** | **0.87** | **41.3 d** |

The walker-vs-TCC gap is 26.0 days. Hold adherence at 1.0 and it becomes
**−0.1 days**: 99% of the gap is behavioural. The removable walker is
mechanically the *better* device — 87% vs 85% peak-pressure reduction — and it
loses because it comes off during 72% of daily activity (Armstrong 2003). The
model reproduces the entire published TCC advantage without ever giving the TCC
a better pressure profile, which reframes "use a total contact cast" as
"make the device irremovable," and explains why an instant-TCC rendered
irremovable performs like a TCC.

### 4.3 Topical growth factor is a race against protease

Becaplermin is cleared from the wound at `KOUT_BEC + KDEG_BEC_PROT · (MMP-9/TIMP-1)`.

| Arm | t_close | MMP-9/TIMP-1 | wound becaplermin (d28) | PDGF drive |
|---|---:|---:|---:|---:|
| Standard care | 67.5 d | 1.43 | — | 0.203 |
| + becaplermin | 62.9 d | 1.46 | 0.134 | 0.345 |
| + TLC-NOSF | 55.7 d | 0.82 | — | 0.254 |
| + both | 51.7 d | 0.85 | **0.226** | **0.484** |

**At the receptor the interaction is genuinely superadditive.** Protease
modulation raises wound becaplermin 1.69-fold at an unchanged dose, and the
PDGF-drive gains are +0.142 (drug alone), +0.051 (dressing alone) and +0.281
together — **1.45× the sum of the singles**.

**At the whole wound it washes out.** Days gained are 4.5, 11.8 and 15.8 against
a sum-of-singles of 16.3, and in rate space the arms are almost exactly
multiplicative: 1.072 × 1.212 = 1.299 predicted vs 1.306 observed. The reason is
mechanistic and worth stating plainly: MMP-9 excess damages the wound through
*several parallel channels* — matrix degradation, granulation loss, inhibition
of keratinocyte migration — and the dressing repairs all of them whether or not
a growth factor is present. So a combination that is unambiguously working
better at the target looks disappointing at the endpoint. This is a concrete
mechanistic account of why combination wound-care trials so often fail to beat
the sum of their parts.

### 4.4 Oxygen gates everything, so revascularisation is not an adjunct

One Hill gate `OG = pO₂³/(pO₂³ + 26³)` multiplies collagen synthesis (prolyl
hydroxylase is O₂-dependent), angiogenesis, the M1→M2 switch, the neutrophil
oxidative burst and the epithelial edge. Below it, nothing anabolic happens.

| Arm | area at 12 weeks | t_close |
|---|---:|---:|
| Ischaemic wound (perfusion 0.32) + becaplermin | 3.140 cm² (**grew**) | never |
| Revascularise on day 14, then the identical drug | 0.767 cm² | 119.2 d |

Same drug, same dose, same dressing schedule; 4.1-fold difference in 12-week
wound area, and one arm never closes at all. On top of this the model carries
the **HIF paradox**: the diabetic wound is hypoxic, so the hypoxia signal is
high, but methylglyoxal adducts on HIF-1α/p300 block transactivation, so the
hypoxia response the wound is entitled to never arrives (Thangarajah 2009). The
wound is starved of oxygen *and* deaf to the starvation.

### 4.5 The biofilm makes the antibiotic a partial drug

| Arm | mean biofilm (0–56 d) | realised kill | bioburden nadir | t_close |
|---|---:|---:|---:|---:|
| 14-day antibiotic alone | 0.954 | **12.7%** of potential | 7.26 log | 77.6 d |
| + weekly sharp debridement | 0.828 | **24.3%** of potential | 6.27 log | 68.8 d |

**87% of the antibiotic course is spent on an organism it cannot reach** — 76%
even with weekly debridement, because the biofilm returns to tolerance within
2–3 days and weekly debridement therefore controls it for well under half the
time. The corollary the model makes explicit is that **debridement's benefit is
not durable bioburden control**: it is the removal of the senescent,
epibolised, hyperkeratotic edge (which the model tracks separately, and which is
where most of the 8.8-day gain actually comes from) plus a short pharmacological
window. If you want the antibiotic to work, dose it into that window.

### 4.6 Osteomyelitis: the nidus beats the exposure

| Arm | osteomyelitis burden at 12 wk | amputation risk at 1 y | antibiotic-days |
|---|---:|---:|---:|
| 6-week antibiotic only | 0.731 | 9.4% | 42 |
| Bone resection d7 + 3-week antibiotic | **0.245** | **3.2%** | **21** |

Half the antibiotic exposure, a third of the amputation risk, and a bone burden
three times lower — because the model's bone-infection term regrows from an
*exposed-bone nidus* state, and resection removes the substrate the regrowth
feeds on. Doubling the duration of a drug that reaches 22% of its potential
effect in bone cannot substitute for removing what it is failing to sterilise.

### 4.7 HbA1c is a slow-axis drug: it treats the *next* ulcer

| | t_close | recurrence @ 12 mo | neuropathy at 1 y | at 5 y | recurrence in year 4 alone |
|---|---:|---:|---:|---:|---:|
| Standard care | 67.5 d | 41.8% | 0.781 | 0.781 | 40.5% |
| + HbA1c 9.0 → 7.0 | 62.2 d | 40.8% | 0.718 | 0.542 | 35.2% |

Intensive glycaemic control buys **5 days** on the ulcer in front of you and
1 percentage point of 12-month recurrence, because neuropathy relaxes with a
~4.6-year time constant. Its payoff arrives in year 4 and later, as the
neuropathy state finally moves from 0.78 to 0.54. This is exactly the pattern in
the literature — a Cochrane review finds the effect of glycaemic control on
index-ulcer healing uncertain, while HbA1c is a robust predictor of ulcer
*incidence* — and the model reproduces both facts from one mechanism rather than
two assumptions.

### 4.8 Closure is remission, not cure

| Arm | t_close | recurrence @ 12 mo |
|---|---:|---:|
| Optimal bundle (TCC + weekly debridement + TLC-NOSF + becaplermin) | 33.6 d | 41.9% |
| The same bundle + remission care (footwear, temperature monitoring, education, surveillance, HbA1c 7.0) | 32.6 d | **16.2%** |

The two arms close the ulcer within **1.0 day** of each other and then diverge by
**26 percentage points**. Essentially the entire benefit of the second arm is
earned after the wound is closed — at the point where most wound-care trials
stop measuring. The model's post-closure hazard is driven by the substrate that
closure did not touch: neuropathy (unchanged), deformity and plantar pressure
(unchanged), and a scar that plateaus at ~80% of normal tensile strength. Under
this structure "healed" is a state the model *keeps simulating*, and the
question "which dressing closed it faster" turns out to be worth about a day
against a question nobody randomised.

---

## 5. Calibration

`python3 dfu_reference_implementation.py --anchors`

| Anchor | Published | Model |
|---|---|---|
| Armstrong 2001 — TCC, mean days to healing (A₀ = 1.4 cm²) | 33.5 d | 36.3 d |
| Armstrong 2001 — removable cast walker | 50.4 d | 59.6 d |
| Armstrong 2001 — half-shoe | 61.0 d | 64.4 d |
| Armstrong 2003 — RCW wear fraction | 28% of daily activity | encoded as `ADHERENCE = 0.28` |
| Armstrong 2017 — recurrence at 1 year | ~40% | 41.8% |
| Wieman 1998 — becaplermin | 50% vs 35% closure at 20 wk | 7% faster healing rate |
| Edmonds 2018 (Explorer) — TLC-NOSF | 48% vs 30% closure at 20 wk | 17% faster healing rate |
| Wolcott 2010 — biofilm reformation after debridement | 24–72 h | 2–3 d (`KG_BIOF = 1.45 /day`) |
| Ceri 1999 — biofilm antibiotic tolerance | 100–1000× planktonic MIC | `TOL_BIOF = 0.92` |

**Known limitations.** (i) Trajectories are deterministic and single-patient, so
population endpoints like "89.5% healed at 12 weeks" appear as time-to-closure
rather than as a proportion. (ii) The post-closure re-ulceration hazard is
constant in time-since-closure: this matches the 1-year figure (41.8% vs ~40%)
but over-predicts 3-year recurrence (79.4% vs ~60%), because in reality the
hazard declines as the highest-risk patients recur early and the scar matures.
(iii) Drug amounts are dose-proportional surrogate units; the PK structure
(depot → wound tissue → bone) is mechanistically correct but not fitted to
measured tissue concentrations. (iv) Deep anatomical structure — tendon sheath
spread, compartment anatomy, angiosome-specific revascularisation — is absent.

---

## 6. Running it

```r
library(mrgsolve); library(dplyr); library(ggplot2)
mod <- mread("dfu_mrgsolve_model.R")

# standard care: removable cast walker at real-world adherence, debride q2wk
debride <- ev(time = seq(0, 140, by = 14), amt = 1, cmt = 11)   # DEBR compartment
out <- mod %>%
  param(DEV_EFF = 0.87, ADHERENCE = 0.28, AREA0 = 2.0) %>%
  mrgsim(events = debride, end = 540, delta = 0.25)

plot(out, AREA + RADIUS + TCPO2 + PROT_o ~ time)

# now change one number: put the same device on the foot all day
mod %>% param(DEV_EFF = 0.87, ADHERENCE = 1.00, AREA0 = 2.0) %>%
  mrgsim(events = debride, end = 540) %>% plot(AREA ~ time)
```

```bash
# 12-tab interactive dashboard
Rscript -e 'shiny::runApp("dfu_shiny_app.R")'

# no R? the reference implementation reproduces every number above
python3 dfu_reference_implementation.py
```

---

## ⚠️ Disclaimer

This is an educational and research QSP model. It is built from public
literature and clinical-trial summary data, has not been independently validated
or fitted to patient-level data, and **must not be used for clinical
decision-making, prescribing, or regulatory submission.** Parameters are
illustrative approximations.
