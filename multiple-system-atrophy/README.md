# 다계통 위축 (Multiple System Atrophy, MSA) — QSP Model

> Quantitative Systems Pharmacology model of multiple system atrophy.
> MSA's dominant clinical problem is not a rate — it is a **broken control
> loop**. Upright blood pressure is written as the **product** of seven
> independently breakable gains, and the whole therapeutic logic of the
> disease follows from *which factor* is broken rather than from how severe
> the hypotension is.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`msa_qsp_model.dot`](msa_qsp_model.dot) |
| 🖼️ Map (SVG, zoomable)   | [`msa_qsp_model.svg`](msa_qsp_model.svg) |
| 🖼️ Map (PNG)             | [`msa_qsp_model.png`](msa_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`msa_mrgsolve_model.R`](msa_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`msa_shiny_app.R`](msa_shiny_app.R) |
| 📚 References (139, PubMed-verified) | [`msa_references.md`](msa_references.md) |

**Scale:** 208-node / 329-edge / 22-cluster mechanistic map · **61-ODE**
mrgsolve model with 36 captured outputs · 19 prebuilt scenarios · **52
self-test assertions, all passing** under mrgsolve 2.0.1 · **10-tab** Shiny
dashboard · 139 references whose PMIDs were programmatically verified against
NCBI E-utilities.

---

## 1. The organising idea

```
ΔMAP(upright) = [baroreceptor afferent signal]        ← INTACT in MSA
              × [central integration NTS→CVLM→RVLM]   ← G_CENT: dying
              × [preganglionic IML output]            ← G_CENT: dying
              × [POSTGANGLIONIC releasable NE pool]   ← POSTG: SPARED
              × [vascular α1-adrenoceptor density]    ← A1R: UP-regulated
              × [effective circulating volume]        ← ECF: eroded at night
              × [venous capacitance / splanchnic pooling]
```

MSA drives the **central** terms toward zero and leaves the
**postganglionic** term comparatively intact. Pure autonomic failure (PAF)
and Parkinson disease (PD) do the opposite. In this model all four
phenotypes are produced by changing **only the `VULN_*` vulnerability
parameters** — never by changing an equation. Four clinically important
behaviours then *fall out of the structure* rather than being coded.

### (1) Drug selectivity is a map of *where* the product is broken

Atomoxetine enters the model at exactly one place: as a multiplier on NET
clearance of synaptic norepinephrine. Its pressor effect is therefore
proportional to what is already being released,
`NEREL = SNA × NEVES × POSTG`. Midodrine's active metabolite instead enters
the **agonist term directly**, downstream of every lesion.

Nothing in the code says "atomoxetine is MSA-selective". Simulating the same
18 mg dose in the two phenotypes gives it anyway — and at **matched
orthostatic severity** (ΔSBP −24 vs −21 mmHg), so the split cannot be an
artefact of one phenotype simply being sicker:

| Phenotype | `G_CENT` | `POSTG` | supine NE | midodrine 10 mg | droxidopa 300 mg | atomoxetine 18 mg |
|---|---|---|---|---|---|---|
| **MSA-P**, yr 7 | 0.016 | **0.91** | 138 pg/mL | +36.4 mmHg | +29.5 mmHg | **+21.2 mmHg** |
| **PAF**, yr 5 (postganglionic mirror image) | **0.63** | 0.075 | 31 pg/mL | +21.1 mmHg | +6.7 mmHg | **+3.4 mmHg** |

Midodrine works in both because it is downstream of every lesion. Droxidopa
and atomoxetine both need the terminal — droxidopa for its AADC, atomoxetine
for the transporter and the release it amplifies — so both collapse in PAF.
The low supine norepinephrine in PAF against a normal value in MSA falls out
of the same structure, because plasma NE in this model is spillover from
release rather than a synaptic concentration.

Two terms carry the whole dissociation:

```c
NEUP = KUP*POSTG*(1.0 - NETI)*NES;        // NET lives ON the terminal
NEPL = NEPL0*(NEREL/KNEREF)*(1.0 + KSPILL*NETI);   // plasma NE is spillover
```

Gating reuptake capacity on `POSTG` is what makes NET inhibition
phenotype-selective: where the terminal is gone there is no transporter left
to block. Remove that gate and a PAF patient gets nearly the same atomoxetine
response as an MSA patient — the opposite of what the crossover trials show.

### (2) Carbidopa antagonises droxidopa — one shared enzyme, written once

Droxidopa becomes norepinephrine only through

```c
DROXNE = KAADC * CDRX * POSTG * (1.0 - CARBI);
```

`POSTG` is the residual terminal supplying peripheral AADC; `CARBI` is
carbidopa's occupancy of **that same enzyme**. Carbidopa is given for the
parkinsonism, not for the blood pressure — and it costs **48% of the
droxidopa response** (+29.5 → +15.5 mmHg at year 7) without either drug
"knowing" about the other.

### (3) The pressor ceiling is set at night, not by daytime potency

Both sodium and water excretion rise exponentially with mean pressure, so
the hours spent recumbent spend the very volume the patient needs at 07:00:

```
supine hypertension → pressure natriuresis → overnight ECF loss
   → smaller central volume on standing → WORSE morning orthostatic
   hypotension → more pressor → …
```

Because the loop is closed in the equations, the timing rules come out as
model *outputs*. Day 7 of a 7-day run, MSA-P at disease year 7:

| arm | night supine peak SBP | overnight Na⁺ loss | overnight urine | dawn ECF | morning standing SBP |
|---|---|---|---|---|---|
| untreated | 98.6 | 51.4 mmol | 0.600 L | 15.24 L | 79.0 |
| midodrine t.i.d. (**correct** timing, last dose 11:00) | 98.3 | 48.7 mmol | 0.581 L | 15.14 L | 78.7 |
| **midodrine at 22:00 (wrong timing)** | **125.8** | **72.4 mmol** | 0.700 L | **14.11 L** | 86.6 |
| head-up tilt sleeping | **93.5** | 47.0 mmol | 0.575 L | **15.41 L** | 79.3 |
| bedtime desmopressin | 98.6 | 51.7 mmol | **0.475 L** | 15.27 L | 79.1 |
| head-up tilt + desmopressin | **93.5** | 47.2 mmol | **0.457 L** | **15.44 L** | 79.4 |

Three things worth reading off this table:

1. **Correctly timed midodrine causes no supine hypertension at all** (98.3 vs
   98.6 untreated) while the same total daily dose given at 22:00 produces
   frank supine hypertension of 126 mmHg. The >4 h-before-recumbency rule is
   an output here, not an assumption.
2. **The harm of bedtime dosing is mechanistically visible**: overnight sodium
   loss rises 51 → 72 mmol and the patient arrives at dawn 1.1 L down. This
   does not settle: the pressor is on board every night, so natriuresis never
   returns to intake and the deficit keeps widening.
3. **Head-up tilt and desmopressin act only on the night and still help the
   day** — they are the only arms that leave dawn volume *higher* than
   untreated, and they achieve it with no daytime pressor whatsoever.

One honest qualification: in the model, 22:00 midodrine still leaves enough
drug on board at 07:30 to raise morning standing pressure (86.6 vs 79.0), so
the model does **not** reproduce "bedtime dosing makes the next morning
worse" as a net blood-pressure statement. What it does reproduce is the
mechanism of the harm — supine hypertension and progressive volume depletion
— which is the reason the rule exists.

### (4) Levodopa failure in MSA-P is postsynaptic

Striatal output is `STRIAT = G_POST × DA/(EC50DA + DA)` with
`G_POST = NMSN^1.15`. Levodopa can raise only `DA`. In the PD phenotype the
medium spiny neuron survives (`NMSN ≈ 1`) and levodopa nearly normalises
output; in MSA-P `NMSN` falls, so the **same brain exposure buys a shrinking
benefit**:

| disease year | MSA-P `NMSN` | PD `NMSN` | levodopa benefit, MSA-P | levodopa benefit, PD |
|---|---|---|---|---|
| 2 | 0.96 | 1.00 | 85 % | 86 % |
| 4 | 0.80 | 0.99 | 74 % | 77 % |
| 6 | 0.54 | 0.96 | 37 % | 73 % |
| 8 | 0.29 | 0.93 | **17 %** | **72 %** |

An early partial response that wanes over 1–2 years, and no capacity for
dyskinesia because there is no supersensitive target to over-drive — both
red flags, one equation.

---

## 2. Why MSA is modelled as a **glial** disease

α-synuclein aggregates in **oligodendrocytes**, not neurons. The cascade in
this model is self-amplifying and its trigger is myelin failure, not
neuronal stress:

```
MYE ↓  →  p25α/TPPP leaves the sheath for the soma  →  templates nucleation
       →  oligomer  →  GCI  →  releases seed  →  taken up by the NEXT
       oligodendrocyte  →  MYE ↓ …
```

Neuronal death is **downstream**, driven largely by the
`WTROPH × (1 − MYE)` term — the trophic support the oligodendrocyte stopped
providing. One shared `NDRIVE` is then applied to nine neuronal populations
with different vulnerabilities, which is what makes the phenotypes a
parameter choice rather than nine separate models. It also states the
therapeutic implication plainly: a neuron-directed agent is aimed at the
wrong compartment.

---

## 3. Two timescales, one model

Neurodegeneration has time constants of **years**; the baroreflex has time
constants of **seconds** (`TAU_MAP = 0.002 h ≈ 7 s`). Rather than fight the
stiffness, the driver uses a two-stage architecture:

| function | what it does |
|---|---|
| `msa_history(years, phenotype)` | circadian posture, 14-year natural history |
| `msa_state_at(year, phenotype)` | freezes the full 61-state vector at a chosen disease year |
| `msa_tilt(state, …)` | restarts from that state, 40-min tilt at 30-s resolution → ΔSBP₃ₘᵢₙ, ΔHR/ΔSBP |
| `msa_day(state, …)` | restarts from that state, 48-h circadian run with real dosing |

**The tilt test therefore becomes a readout of disease stage**, computed from
the same 61 equations that generate the 14-year trajectory:

| disease year | supine SBP | 3-min SBP | ΔSBP | ΔHR | ΔHR/ΔSBP | upright NE rise |
|---|---|---|---|---|---|---|
| 1 | 116.5 | 109.6 | −6.9 | +17.3 | 2.51 | +22 % |
| 3 | 114.8 | 107.5 | −7.4 | +13.3 | 1.81 | +28 % |
| 5 | 108.1 | 95.5 | −12.6 | +7.0 | **0.55** | +37 % |
| 7 | 96.4 | 72.0 | **−24.4** | +1.4 | **0.06** | +16 % |
| 9 | 94.7 | 65.7 | **−28.9** | +0.1 | **0.004** | +1 % |

The ΔHR/ΔSBP ratio crossing below 0.5 — the number that makes orthostatic
hypotension *neurogenic* — is produced, not prescribed, and it crosses
because vagal withdrawal on standing is gated on medullary survival in the
same way the sympathetic arm is. A healthy control run through the identical
code holds at ΔSBP −6.9 mmHg with a ratio of 2.51.

Note the non-monotone NE column: the *percentage* increment peaks in
mid-course, because a partially surviving reflex acts on an already-low
baseline. The model therefore does **not** support "blunted percentage NE
rise" as an early marker — only the late collapse (+1% at year 9) and the
normal *supine* value are robust.

### Two engineering notes that mattered

**`TSCALE`.** A 61-state stiff system with a numerical Jacobian spends
essentially all of its time resolving a 7-second transient that no annual
endpoint depends on. `TSCALE` stretches *only* the fast haemodynamic time
constants, leaving every quasi-steady value untouched; history runs use
`TSCALE = 40`, tilt and 48-h runs use `TSCALE = 1`. The self-test asserts
the two agree — the disease layer and every clinical endpoint match to
**better than 10⁻¹⁰ relative difference**, and an 8-year run drops from
minutes to ~30 s.

**`INITMAIN`.** mrgsolve applies `$MAIN` *after* `init()`, so unguarded
`X_0 =` assignments silently reset a restarted patient back to healthy
physiology. Before this was gated, every tilt returned identical numbers for
a healthy control, a year-8 MSA patient and a PAF patient — a bug that looks
exactly like a working model until you compare phenotypes.

---

## 4. Model contents

**61 ODEs.** Glial proteinopathy (10: p25α, monomer, oligomer, GCI, seed,
microglia, myeloperoxidase, oxidative stress, CoQ10, myelin) · neuronal
populations (9: IML, medullary cardiovascular, cardiovagal, nigral, striatal
MSN, olivopontocerebellar, Onuf, respiratory, **postganglionic**) ·
neurofilament light · fast cardiovascular loop (10: SNA, synaptic NE,
vesicular NE, α1 density, TPR, splanchnic pool, central volume, HR, vagal
tone, MAP) · renal/volume/endocrine (6: sodium, ECF, AQP2, AVP, renin,
aldosterone) · striatal dopamine · drug PK (25 compartments across
midodrine/desglymidodrine, droxidopa, fludrocortisone, atomoxetine,
pyridostigmine, desmopressin, levodopa/carbidopa/brain levodopa, ubiquinol,
verdiperstat, anti-α-synuclein antibody, yohimbine, octreotide, water bolus)
· clinical endpoints (5: OHSA, UMSARS-I, UMSARS-II, SCOPA-AUT, cumulative
hazard).

**19 scenarios**, including the four phenotype comparisons, the
droxidopa–carbidopa interaction, the atomoxetine MSA-vs-PAF dissociation,
six arms of the nocturnal loop, pyridostigmine's reflex-only amplification,
levodopa in MSA-P vs PD across disease years, ubiquinol in COQ2-wild-type vs
COQ2-deficient patients, a myeloperoxidase inhibitor that must come out
**null**, and an anti-α-synuclein antibody started early vs late.

**Natural history calibration.** UMSARS-II progresses **6.1 points/year**
(cohorts report 5–8), median survival **8.8 years** from pathology onset
(cohorts report 6–10 from motor onset), post-void residual reaches 135 mL by
year 5, stridor emerges by year 8, and plasma NfL rises 12 → 31 pg/mL and
then declines — consistent with NfL indexing the *rate* of loss rather than
cumulative damage. MSA-C comes out ataxia-dominant (18 vs 11 UMSARS-II
points) and MSA-P parkinsonism-dominant (22 vs 10) from the vulnerability
parameters alone.

---

## 5. Running it

```bash
# 61-ODE model: builds, runs 52 self-test assertions, then 19 scenarios
Rscript msa_mrgsolve_model.R

# Shiny dashboard (10 tabs)
MSA_SHINY_RUN=1 Rscript -e 'shiny::runApp("msa_shiny_app.R", port = 8080)'

# Head-less test of all 10 dashboard panels — needs mrgsolve, not a browser
Rscript -e 'source("msa_shiny_app.R"); msa_shiny_selftest()'

# Re-render the mechanistic map
dot -Tsvg msa_qsp_model.dot -o msa_qsp_model.svg
dot -Tpng -Gdpi=150 -Gsize=72,36 msa_qsp_model.dot -o msa_qsp_model.png
```

Verified under **R 4.3.3 · mrgsolve 2.0.1 · shiny 1.14.0 · Graphviz 2.43**:
**52/52** model assertions and **13/13** dashboard panel checks pass.

---

## 6. Dashboard

| Tab | Content |
|---|---|
| 1 환자 프로파일 | the **seven multiplicands** shown separately, with the drug implication stated |
| 2 기립경사 검사 | tilt BP/HR/NE/cerebral-autoregulation margin + the diagnostic metric table |
| 3 24시간 혈압·야간 루프 | 48-h SBP with night shading, natriuresis, urine output, ECF |
| 4 약물 PK | plasma concentrations plus α1 / NET / AADC target occupancy |
| 5 자율신경 축 | SNA, NE handling, RAAS, AVP–AQP2, and α1 up-regulation over years |
| 6 신경병리 진행 | glial cascade, myelin/mitochondria/inflammation, regional neuronal loss |
| 7 운동증상·레보도파 | phenotype composition and the waning levodopa response vs PD |
| 8 임상 엔드포인트 | UMSARS I/II/IV, SCOPA-AUT, PVR, survival |
| 9 시나리오 비교 | 13 regimens ranked by morning standing SBP against supine hypertension |
| 10 바이오마커·감별 | NfL, seed, MIBG correlate (`POSTG`), four-phenotype discrimination panel |

---

## 7. Known limitations

- **Untreated supine hypertension is under-represented.** Treatment-emergent
  supine hypertension is reproduced sharply (126 mmHg on a 22:00 pressor), but
  in the default parameterisation an *untreated* patient keeps a preserved
  rather than hypertensive supine pressure. The untreated supine-hypertensive
  phenotype — roughly half of MSA patients — is reached by raising `FTON`, the
  lesion-independent tonic sympathetic drive. It is exposed as a slider rather
  than baked in, because the model has nothing to say about which patients
  have it.
- **The blunted upright NE increment is not an early marker in this model.**
  The percentage rise is non-monotone (it peaks mid-course, because a partly
  surviving reflex acts on an already-low baseline). Only the late collapse
  and the normal *supine* value are robust; absolute NE values are calibrated
  to a single reference release rate and should not be over-interpreted.
- **Late-stage pressor responses are large** (+36 mmHg from 10 mg midodrine at
  year 7). Severe nOH genuinely does respond steeply, but the upper end of
  the model's dose-response is less well constrained than the mid-range.
- The **verdiperstat phase-3 (M-STAR) null result could not be verified as a
  PubMed-indexed paper** and is therefore cited only by registry number. The
  scenario should be read as a structural prediction, not as a reproduction
  of a published trial — see §13 of the references.
- Cognitive impairment, dysphagia and cerebellar oculomotor findings are
  present as map nodes and hazard terms but are not separately calibrated.
- Parameters are illustrative approximations fitted to *published summary
  statistics*, not to individual patient data.

---

## ⚠️ 면책 조항 (Disclaimer)

본 모델은 **교육 및 연구 목적의 정량적 시스템 약리학(QSP) 모델**입니다.
공개 문헌과 임상시험 요약값을 바탕으로 보정했으나 독립적으로 검증·인증되지
않았으며, **실제 임상 의사결정, 처방, 또는 규제 제출에 직접 사용해서는 안
됩니다.**
