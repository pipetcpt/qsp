# Chronic Insomnia Disorder (CID) — QSP Model
# 만성 불면장애 QSP 모델

> Quantitative Systems Pharmacology model of chronic insomnia built around the
> claim that **insomnia is not a deficiency of sleep drive**. Process S and
> Process C are close to normal in most patients. What is abnormal is a *third*
> process — **24-hour hyperarousal** — that raises the wake side of the
> VLPO/ascending-arousal flip-flop switch, and a **learned behavioural loop**
> that keeps arousal there. Hypnotics push the switch from the other side.
> CBT-I enters the loop. That difference is why two treatments with the same
> effect size at week 4 behave completely differently at week 12, and
> reproducing it is what this model is for.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`ins_qsp_model.dot`](ins_qsp_model.dot) |
| 🖼️ Map (SVG, full detail) | [`ins_qsp_model.svg`](ins_qsp_model.svg) |
| 🖼️ Map (PNG)              | [`ins_qsp_model.png`](ins_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`ins_mrgsolve_model.R`](ins_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`ins_shiny_app.R`](ins_shiny_app.R) |
| 📚 References             | [`ins_references.md`](ins_references.md) |

**Scale:** 217 nodes · 19 clusters · 329 edges in the map · 63 ODE compartments ·
19 prebuilt scenarios (17 interventions + 2 matched untreated controls) ·
9 dashboard tabs · 105 PubMed-verified references.

---

## 1. Disease in one paragraph

Chronic insomnia disorder is difficulty initiating or maintaining sleep, with
daytime consequences, on at least three nights a week for at least three
months, **despite adequate opportunity to sleep**. That last clause is the
whole problem: these patients are not sleep-deprived, they are sleep-*unable*.
Polysomnography in insomnia is often nearly normal, and the objective deficit
(20–60 min of lost sleep) is far smaller than the subjective one. What is
reliably abnormal is arousal, measured every way it has been measured:
high-frequency beta and gamma EEG intruding on NREM, whole-brain metabolic
rate that fails to fall from wake to NREM on FDG-PET, a flattened evening
cortisol nadir with nyctohemeral HPA activation, raised muscle sympathetic
nerve activity and reduced heart-rate variability. Insomnia is therefore best
modelled as a **disorder of the wake side of the sleep switch**, not of sleep
drive.

The second half of the disease is learned. Spielman's 3-P model says a
predisposed person (trait anxiety, ruminative style, short-sleep phenotype)
is tipped into insomnia by a precipitating stressor, and then kept there by
what they do about it: they go to bed earlier and lie there longer, which
spreads the same homeostatic sleep pressure across more hours and *lowers*
sleep efficiency; the extra wakefulness in bed conditions the bed itself as a
wake cue; and trying harder to sleep is itself arousing. The precipitant can
resolve entirely and the insomnia persists, because the loop is now
self-sustaining.

한국어 요약: 만성 불면장애는 잠잘 기회가 충분한데도 입면 또는 수면 유지가 어려운
상태가 주 3회 이상 3개월 넘게 지속되는 질환입니다. 이 모델의 출발점은 "불면은
수면 구동이 부족해서 생기는 것이 아니다"라는 관점입니다. 실제로 항상성 수면압
(Process S)과 일주기 리듬(Process C)은 대체로 정상이며, 이상이 있는 것은
**24시간 지속되는 과각성**(피질 β/γ 침입, 전뇌 대사율 상승, 야간 코르티솔 nadir
소실, 교감신경 항진)입니다. 여기에 스필만의 3-P 모형이 말하는 **학습된 유지
고리**(침상시간 연장 → 수면압 희석 → 침대의 각성 조건화 → 수면 노력 증가 →
과각성 강화)가 더해지면서 촉발 요인이 사라진 뒤에도 불면이 스스로 지속됩니다.
수면제는 스위치를 반대쪽에서 밀 뿐 이 고리에 개입하지 않고, CBT-I는 고리 자체를
표적합니다. 이 비대칭을 재현하는 것이 본 모델의 핵심 목적입니다.

---

## 2. Mechanistic map

[![CID QSP map](ins_qsp_model.png)](ins_qsp_model.svg)

19 clusters, 217 nodes, 329 edges. Click the image for the full-detail SVG.

| # | Cluster | What it contains |
|---|---|---|
| ① | Process C — SCN pacemaker | light → ipRGC/melanopsin → RHT → SCN; CLOCK:BMAL1 / PER-CRY TTFL; CK1δ/ε; SCN→PVN→IML→SCG→pineal; AANAT/ASMT; MT1/MT2; phase-response curve; wake-maintenance zone |
| ② | Process S — homeostatic drive | waking metabolism → ATP → CD73 → adenosine; A1/A2A; PGD2/L-PGDS; SWA (EEG δ); synaptic homeostasis |
| ③ | Ascending arousal system | LC-NE, DR-5HT, TMN-histamine, LDT/PPT-ACh, vPAG-DA, basal forebrain, parabrachial glutamate, thalamic relay |
| ④ | Sleep-active nuclei & flip-flop | VLPO GABA/galanin, MnPO, parafacial zone, cortical nNOS; GABA-A subunits; the mutual-inhibition bistable switch; spindles/K-complexes |
| ⑤ | Orexin | LH orexin neurons, OX1R/OX2R, switch stabilisation, MCH antagonism, metabolic inputs |
| ⑥ | **24-h hyperarousal** | cognitive/somatic/emotional arousal, β-γ EEG intrusion, raised cerebral metabolic rate, amygdala and salience network, local sleep, HR/HRV, MSNA |
| ⑦ | HPA axis | CRH → ACTH → cortisol, impaired GR/MR feedback, the arousal↔CRH feed-forward loop, IL-6/TNF |
| ⑧ | Sleep architecture | W/N1/N2/N3/REM transitions, SLD–vlPAG REM flip-flop, atonia, ~90-min cycling, fragmentation index |
| ⑨ | **3-P perpetuating loop** | predisposing/precipitating/perpetuating; TIB extension, naps, sleep effort, conditioned arousal, DBAS, monitoring, misperception |
| ⑩ | CBT-I | sleep restriction, stimulus control, cognitive restructuring, relaxation, mindfulness/ACT, TIB re-titration |
| ⑪ | BzRA PK | zolpidem (IR/ER/SL), eszopiclone, triazolam/temazepam, CYP3A4 |
| ⑫ | DORA PK | suvorexant, lemborexant, daridorexant, receptor occupancy, food effect |
| ⑬ | Melatonergic / H1 / 5-HT2A | ramelteon + M-II (and the fluvoxamine interaction), exogenous melatonin, low-dose doxepin, trazodone, mirtazapine, OTC antihistamines, gabapentinoids |
| ⑭ | Exogenous modifiers | caffeine, alcohol (biphasic), shift work, evening blue light, morning bright light, exercise timing, bedroom temperature, nicotine |
| ⑮ | Adverse effects | next-morning residual sedation, psychomotor/driving impairment, falls, complex sleep behaviours, amnesia, tolerance, rebound, dependence, respiratory depression, DORA-specific effects |
| ⑯ | Comorbidity | depression, anxiety, chronic pain, COMISA, RLS, vasomotor symptoms, ageing, ADHD, substance use — all bidirectional |
| ⑰ | Long-term consequences | hypertension/non-dipping, cardiovascular risk, insulin resistance, glymphatic clearance, cognition, accidents, productivity, suicidality |
| ⑱ | Clinical endpoints | SOL/LPS, WASO, TST, SE, awakenings, ISI, PSQI, subjective sleep, N3/REM, daytime function, quality of life, remission |

---

## 3. The mrgsolve model

63 ODE compartments. `time` is in hours and clock time is `TIME mod 24`, so a
schedule (bedtime, rise time, light, dosing) is a first-class part of the model
rather than a covariate.

### 3.1 The sleep switch

Everything reported by the model — SOL, WASO, TST, sleep efficiency, N3 and REM
minutes, awakening count — is read out of **one** continuous variable:

```
drive = wS·S + wA·adenosine_eff + E_BzRA + E_melatonergic + E_5HT2A
        − wT·ΔcoreTemp − C(t) − AROUSAL − W_ascending − ultradian − θ

SLP   = opportunity × SLPmax / (1 + exp(−drive / k))
```

`W_ascending` is itself driven by `1 − SLP`, so the noradrenergic, histaminergic
and orexinergic tones collapse once sleep begins and rise once it ends. That
feedback makes the switch **bistable**, which is the point: sleep onset and
awakening are transitions, not gradual slides, and the model shows genuine
hysteresis.

### 3.2 Where WASO comes from

Awakenings are not stochastic. An ultradian oscillator whose phase advances
**only during sleep** (so it effectively resets at each sleep onset) produces a
~90-minute trough in sleep drive that deepens as Process S dissipates. A good
sleeper's margin over threshold absorbs it; an insomniac's does not. This
reproduces the clinical pattern that fragmentation concentrates in the last
third of the night — when S has dissipated, cortisol is rising and the
ultradian troughs are deepest.

### 3.3 Compartment map

| Block | Compartments |
|---|---|
| PK (23) | zolpidem, eszopiclone, suvorexant, lemborexant, daridorexant, ramelteon + M-II, exogenous melatonin, doxepin, trazodone, caffeine, ethanol (Michaelis–Menten) — each with an absorption compartment |
| Circadian (3) | pineal melatonin `MELP` (pg/mL), phase shift `PHI`, core-temperature deviation `CBTD` |
| Two-process (2) | `S`, extracellular adenosine `ADO` |
| Ultradian (2) | `ULTA`, `ULTB` — a sleep-gated harmonic oscillator |
| Flip-flop (3) | `NE`, `HA`, `OXA` |
| Arousal & HPA (5) | `AROU`, precipitating stressor `STR`, `CRH`, `ACTH`, `CORT` |
| 3-P loop (4) | conditioned arousal `COND`, sleep effort `SEFF`, beliefs `DBAS`, time in bed `TIBS` |
| Drug adaptation (3) | BzRA tolerance `TOLB`, withdrawal/rebound `WDR`, ethanol rebound `ETRB` |
| Sleep measurement (8) | onset flag `ONS`, running `SOLBAR`/`WASOBR`/`TSTBAR`/`SEBAR`/`N3BAR`/`REMBAR`/`RESBAR` |
| Consequences (5) | misperception `MISP`, `ISI`, depressive load `DEP`, `IL6`, glymphatic debt `GLYM` |
| Counters (5) | cumulative SOL, WASO, TST, awakenings, residual-sedation exposure |

### 3.4 Nine patient phenotypes

`good_sleeper`, `chronic_insomnia`, `severe_insomnia`, `sleep_onset`,
`sleep_maintenance`, `elderly`, `comorbid_depression`, `delayed_phase`,
`shift_worker`. A phenotype is **not** a set of target sleep numbers — it is a
starting point on the 3-P loop, defined mostly by constitutive arousal `A0`
plus, where relevant, circadian period, pineal output and SWS capacity. The
sleep metrics are outputs.

### 3.5 Usage

```r
source("ins_mrgsolve_model.R")

# one arm: daridorexant 50 mg for 12 weeks in a chronic-insomnia patient
r  <- ins_run_scenario("7. Daridorexant 50 mg nightly", weeks = 12)
head(r$nightly)          # night-by-night SOL / WASO / TST / SE / ISI ...

# a custom arm
sim <- ins_simulate(
  phenotype      = "elderly",
  weeks          = 8,
  param          = list(PAIN = 0.4),                       # patient covariates
  treat          = list(BLTON = 1, BLTLUX = 10000),        # starts at randomisation
  events         = ev_nightly("DARD", 25, days = 56),
  cbti_start_day = 15
)
ins_nightly(sim)

# every prebuilt arm, summarised at weeks 1, 4 and 12
ins_run_all()
```

Every simulation runs a **150-day untreated lead-in first**, so the patient is
at their chronic equilibrium at randomisation rather than mid-transient. Only
the post-lead-in window is returned, with time re-zeroed.

---

## 4. Scenarios and what they are for

| # | Arm | The question it answers |
|---|---|---|
| 1 | Untreated chronic insomnia | Where does the 3-P loop settle, and does it stay there? |
| 2 | Zolpidem 10 mg × 4 weeks, then stop | Tolerance build-up and rebound insomnia after discontinuation |
| 3 | Zolpidem 3 nights/week | Does intermittent dosing stay under the tolerance driver? |
| 4 | Eszopiclone 3 mg × 12 weeks | A longer half-life buys maintenance and costs residual sedation |
| 5–8 | Suvorexant 20, lemborexant 10, daridorexant 50 and 25 mg | Three DORAs whose half-lives differ ~2-fold by design |
| 9 | Ramelteon 8 mg | A latency drug with no maintenance effect |
| 10 | Melatonin 2 mg, elderly | Replacement when pineal output is 35% of young (control: arm 18) |
| 11 | Low-dose doxepin 6 mg | H1 blockade only: last-third sleep, no onset effect (control: arm 19) |
| 12 | Trazodone 50 mg | Off-label 5-HT2A/H1/α1, raises N3, blunts REM |
| 13 | CBT-I alone | The only arm that changes the perpetuating loop |
| 14 | CBT-I + zolpidem, drug stopped at week 6 | Does the drug add anything the behaviour therapy does not? |
| 15 | Caffeine 200 mg at 16:00 | Adenosine antagonism seven hours before bed |
| 16 | Alcohol 40 g at 22:00 | Shorter latency, second-half rebound, REM suppression |
| 17 | Morning bright light, delayed phase | Treating the clock rather than the switch |
| 18–19 | Untreated elderly / untreated maintenance-type | Matched controls for arms 10 and 11 |

---

## 5. Validation

Full numbers, including week-1/4/12 tables and the change from each arm's
matched control, are in the calibration block at the end of
[`ins_mrgsolve_model.R`](ins_mrgsolve_model.R). The headline comparison,
week 4, change from the matched untreated control:

| Intervention | ΔSOL | ΔWASO | ΔTST | ΔISI | Literature (placebo-subtracted) |
|---|---:|---:|---:|---:|---|
| Zolpidem 10 mg | −21 | −51 | +70 | −6.6 | LPS −20…−30, WASO −20, TST +30…+50 |
| Eszopiclone 3 mg | −20 | −71 | +91 | −8.6 | LPS −15…−20, WASO −25, TST +45…+60 |
| Suvorexant 20 mg | −23 | −47 | +70 | −6.3 | LPS −10…−20, WASO −20…−25 |
| Lemborexant 10 mg | −25 | −44 | +70 | −6.0 | LPS −20, WASO −25, TST +40…+60 |
| Daridorexant 50 mg | −24 | −47 | +70 | −6.5 | LPS −11, WASO −18, sTST +22 |
| Ramelteon 8 mg | −18 | −53 | +71 | −6.0 | LPS −9…−13, no WASO/TST effect |
| Melatonin 2 mg (elderly) | −14 | **+5** | +15 | −1.5 | sSOL −20, **no** objective WASO effect |
| Doxepin 6 mg | −13 | −19 | +40 | −3.0 | WASO −20…−25, TST +25…+30, no SOL effect |
| CBT-I | −27 | −76 | +46 | −12.5 | SOL −20, WASO −25, SE +10, ISI −8…−10 |
| Caffeine 200 mg 16:00 | **+44** | −19 | −24 | +3.2 | 400 mg 6 h pre-bed costs ~1 h of TST |
| Alcohol 40 g | −16 | **+42** | −24 | +3.4 | shorter latency, second-half fragmentation |

The untreated baseline the model settles on for the chronic-insomnia phenotype
— SOL 47 min, WASO 100 min, TST 399 min, SE 74.6%, ISI 18.8, subjective TST
292 min — matches the pooled entry criteria and baselines of the recent Phase 3
hypnotic programmes. Sleep architecture comes out in the right direction too:
DORAs raise REM (130–137 min vs 88 untreated), ethanol suppresses it (56),
trazodone raises N3 (97 vs 86), BzRA suppress N3 slightly, and the elderly
phenotype retains about 40% of young N3 (37 min).

### Where the model is honest about missing

- **There is no placebo response.** Real placebo arms improve LPS by 10–20 min
  and WASO by 15–25 min, so every Δ above is larger than the corresponding
  placebo-subtracted trial delta by roughly that amount. Comparisons *between
  active arms* are the ones this model is built to make.
- **There is a TST ceiling.** Total sleep time cannot exceed `SLPMAX × TIB`, so
  once a hypnotic consolidates sleep to near continuity, every effective drug
  converges on the same TST at a fixed time in bed. The model therefore
  **cannot** rank suvorexant, lemborexant, daridorexant 50 mg and zolpidem on
  total sleep time — only on latency, on the shape of the night, on residual
  sedation, and on what happens when the drug stops.
- **CBT-I is a responder, not an average.** Arm 13 reaches ISI 6.4 at week 4
  and 4.9 at week 12. Real CBT-I trials report ~60–70% response and ~40–50%
  remission; a single deterministic patient cannot represent that spread.
- **Awakenings repeat.** The ultradian mechanism gives the same patient the
  same awakening pattern every night.
- **The loop is genuinely bistable**, which is deliberate (it *is* Spielman's
  model taken literally) but means outcomes are steep in `A0` near the
  threshold — small covariate changes there produce large predicted differences.

---

## 6. The claim the model exists to test

Both CBT-I and a hypnotic drop ISI by a similar amount at week 4. At week 12,
with nothing being taken, CBT-I has held and zolpidem has not — and arm 2 ends
*worse* than the untreated arm (ISI 22.0 vs 19.7) because the rebound nights
re-conditioned the bed.

No parameter was set to produce that. It follows from **where the two act**:
CBT-I removes terms from the arousal equation (`COND`, `SEFF`, `DBAS`) and
compresses `TIBS`, all of which are state variables that persist; the drug adds
a term proportional to a concentration, which is gone when the concentration
is. If that asymmetry did not fall out of the structure, the model would be
wrong about what insomnia is.

---

## 7. Shiny dashboard

```r
shiny::runApp("ins_shiny_app.R")
```

Nine tabs: patient profile and arousal decomposition · two-process core and the
switch (with the hysteresis loop) · circadian phase, light and melatonin · drug
PK and target occupancy · nightly sleep metrics against a matched control ·
sleep architecture · the 3-P perpetuating loop · clinical endpoints and safety ·
side-by-side scenario comparison. Every run simulates the matched untreated
control alongside the active arm.

---

## 8. Files

```
chronic-insomnia-disorder/
├── README.md                  이 문서
├── ins_qsp_model.dot          기계론적 지도 소스 (Graphviz)
├── ins_qsp_model.svg          지도 (벡터, 전체 해상도)
├── ins_qsp_model.png          지도 (150 dpi)
├── ins_mrgsolve_model.R       63-ODE QSP 모델 + 19개 시나리오 + 검증 블록
├── ins_shiny_app.R            9-탭 인터랙티브 대시보드
└── ins_references.md          105편 (PubMed 검증 완료)
```

Rebuild the map with:

```bash
dot -Tsvg ins_qsp_model.dot -o ins_qsp_model.svg
dot -Tpng -Gdpi=150 ins_qsp_model.dot -o ins_qsp_model.png
```

---

## ⚠️ Disclaimer

교육 및 연구 목적의 정성적·반정량적 QSP 모델입니다. 공개 문헌과 임상시험 데이터를
바탕으로 구성되었으나 독립적으로 검증·인증되지 않았으며, **실제 임상 의사결정,
처방, 또는 규제 제출에 직접 사용해서는 안 됩니다.** 특히 수면제 선택·용량·중단은
반드시 의료진과 상의하십시오.
