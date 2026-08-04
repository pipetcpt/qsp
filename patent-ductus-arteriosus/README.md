# Patent Ductus Arteriosus of Prematurity — QSP Model
### 미숙아 동맥관 개존증 정량적 시스템 약리학 모델

| 산출물 | 파일 |
|--------|------|
| 기계론적 지도 | [`pda_qsp_model.dot`](pda_qsp_model.dot) · [SVG](pda_qsp_model.svg) · [PNG](pda_qsp_model.png) |
| mrgsolve 모델 | [`pda_mrgsolve_model.R`](pda_mrgsolve_model.R) |
| Python 참조 구현 | [`pda_reference_model.py`](pda_reference_model.py) |
| Shiny 대시보드 | [`pda_shiny_app.R`](pda_shiny_app.R) |
| 참고문헌 | [`pda_references.md`](pda_references.md) |

---

## The one-sentence version

The ductus arteriosus is **not a lesion** — it is a normal fetal vessel held open
by its own prostaglandin E₂ — and once you write it that way, the three things
that have made this field intractable become three separate variables:
**constriction** (which drugs cause), **closure** (which needs ductal wall
hypoxia, and which the preterm wall is too thin to supply), and **burden**
(the time-integral of shunt, which is what the outcomes actually see).

동맥관은 병변이 아니라 자기 자신의 프로스타글란딘 E₂가 열어 두는 정상 태아
혈관입니다. 그렇게 쓰는 순간 이 분야를 40년간 교착시킨 세 가지가 서로 다른
변수로 분리됩니다 — **수축**(약이 만드는 것), **폐쇄**(관벽 저산소증이
필요하고, 미숙아의 얇은 벽은 그것을 만들지 못함), **부담**(단락의 시간 적분,
결과가 실제로 보는 것).

---

## 1. Why this disease needs a QSP model

Every controversy here is quantitative, not qualitative.

COX inhibitors close the duct. They do not improve the outcomes the duct is
blamed for. **Baby-OSCAR** (2024, n=653) closed ducts with early ibuprofen and
got death-or-moderate/severe-BPD of **69.2% vs 63.5%** — numerically *worse*,
with death 13.6% vs 10.3%. **BeNeDuctus** (2023, n=273) found expectant
management non-inferior, 46.3% vs 63.5% on its composite. **TIPP** (2001,
n=1202) halved PDA (24% vs 50%) and cut severe IVH (9% vs 13%) and changed
18-month death-or-disability not at all (47% vs 46%).

A model that cannot reproduce *"the drug works and the patient does not
benefit"* is not modelling this disease.

And the reason it happens is measurable. **Semberova 2017** followed 280 infants
≤1500 g under a genuinely non-interventional policy: **85% closed
spontaneously before discharge**, with medians of **71 d below 26+0 weeks,
13 d at 26+0–27+6, 8 d at 28+0–29+6 and 6 d at ≥30 weeks**. If nearly every
duct closes eventually, then a drug's effect on any fixed-time closure endpoint
is *bounded above by the untreated non-closure fraction at that time*. That
ceiling — not drug potency — is what the negative trials measured.

---

## 2. The four structural commitments

### 2.1 Two enzyme sites, not one drug class

Indomethacin and ibuprofen compete with arachidonate in the cyclooxygenase
**channel**. Acetaminophen reduces the ferryl-protoporphyrin radical at the
physically separate **peroxidase site**, so it competes with **peroxide**, not
with substrate. Written correctly — channel occupancy additive, peroxidase term
multiplicative with an IC₅₀ that scales with peroxide tone —

```
I_chan  = Σ(Cu/Ki) / (1 + Σ(Cu/Ki))                    ibuprofen + indomethacin
I_perox = Cu / (Cu + IC50·(1 + peroxide/K))             acetaminophen only
I_COX   = 1 − (1 − I_chan)(1 − I_perox)                 one catalytic cycle
```

the model **predicts** that acetaminophen loses potency under
chorioamnionitis while ibuprofen does not. That is a falsifiable consequence of
getting the enzymology right, not a fitted effect.

### 2.2 Constriction is not closure — and wall O₂ is the switch

Permanent closure requires ductal **wall** hypoxia: the constricted lumen stops
supplying O₂ by diffusion, the media goes hypoxic, and HIF-1α/VEGF/TGF-β1 build
a neointimal cushion. Following **Kajino 2002**, the media has two supplies with
*opposite* thickness dependence:

- **diffusion** from lumen and adventitia — adequate only across a **thin** wall,
  falling with the square of diffusion distance;
- **vasa vasorum** — present only in a **thick** wall, and abolished by the
  constriction itself.

So a thin preterm wall stays oxygenated while constricted → no cushion → the
duct **reopens** when drug clears. A thick term wall loses its vasa and goes
profoundly hypoxic → permanent closure. Because the wall also thickens
postnatally, the same duct constricted at day 20 remodels where it would not
have at day 2 — so the model predicts that **early treatment reopens more often
than late treatment**.

### 2.3 Tone is a sigmoid of net drive, not a product of gains

This one was a bug before it was a commitment, and it is worth stating plainly.
The first version wrote

```
tone = Tmax(GA) × O2gain(GA) × (1 − relax_PGE2) × (1 − relax_NO)     # WRONG
```

which makes achievable tone the *arithmetic product of every immaturity factor*.
A 26-week duct then could not close at any drug exposure — peak tone was **0.02
at 92% ductal COX inhibition**. That is flatly wrong: roughly half of treated
26-week infants do close. The gestational limit belongs somewhere else entirely
— in `TMAXGA`, the maximal occlusion a thin immature media can generate:

| GA | unconstricted d_max | occlusion ceiling | residual lumen at full tone | closable by drug? |
|----|--------------------|-------------------|------------------------------|-------------------|
| 24 wk | 1.89 mm | 0.798 | **0.435 mm** | **no** — above the 0.30 mm threshold |
| 25 wk | 2.03 mm | 0.873 | **0.320 mm** | knife-edge |
| 26 wk | 2.16 mm | 0.922 | **0.238 mm** | yes |
| 28 wk | 2.43 mm | 0.971 | 0.083 mm | yes |
| 38 wk | 3.78 mm | 1.000 | 0.005 mm | yes |

The steep clinical GA gradient in treatment success therefore falls **out** of
the model rather than being fitted **into** it. And it locates
Semberova's near-discontinuity — 71 d below 26 weeks versus 13 d just above —
at the same place the model puts its structural threshold, ~25.5 weeks.

### 2.4 Burden, not closure, is what the outcomes see

```
d(PDA burden)/dt = max(0, Q_shunt − Q_significant) / 24        (mL/min/kg)·day
```

Shunt itself is hydraulic: ductal resistance ∝ **1/d⁴**, which is why
"moderate-to-large ≥1.5 mm" is a real threshold and not a convention — a 2 mm
duct and a 1 mm duct differ 16-fold in resistance. And the shunt is *unmasked*
by the postnatal fall in PVR over 2–5 days, which is why the "day-3 large PDA"
exists and why prophylaxis within 24 h treats ducts that do not yet matter.

---

## 3. The model

**51 ODEs**, time in hours of postnatal age, drug amounts per kg.

| block | compartments |
|-------|-------------|
| Drug PK | ibuprofen (2-cmt + ductal effect site), indomethacin (2-cmt + slowly-reversible effect site), acetaminophen (gut + 2-cmt) |
| Hepatic safety | glutathione, ALT (NAPQI route) |
| Prostanoid | ductal PGE₂, circulating PGE₂, pulmonary 15-PGDH, cAMP, NO tone, EP4 density, ductal peroxide tone, SMC Ca²⁺ |
| Ductus | tone, neointimal cushion, wall thickness, VEGF, TGF-β1 |
| Haemodynamics | PVR, SVR, LA pressure, LV dilation, lung water, compliance, PDA burden, mesenteric and cerebral deficits |
| Kidney | GFR, creatinine, urine output, fluid balance, renal prostanoid tone |
| Gut / platelet / bilirubin | mucosal PGE₂, thromboxane, bleeding time, total and free bilirubin |
| Development / outcome | alveolarisation, ventilator exposure, cumulative COX exposure, and 5 hazard integrals (BPD, NEC, IVH, SIP, death) |

COX appears **four times** — ductus, kidney, gut mucosa, platelet — because
efficacy and toxicity are the same molecular event in different organs. The only
things that separate the two NSAIDs are organ-specific prostanoid dependence and
the **non-COX vasoconstriction** that indomethacin has and ibuprofen does not.

---

## 4. What the model says

### 4.1 Drug exposure, and a half-life that is shorter than it should be

At 26 weeks and 0.80 kg:

| | ibuprofen 10-5-5 q24h | indomethacin 0.2-0.1-0.1 q24h | acetaminophen 15 q6h |
|---|---|---|---|
| C_max (total) | 71.4 mg/L | 0.800 mg/L | 35.4 mg/L |
| unbound C_max | 4.15 µM (f_u 1.2%) | 0.0224 µM (f_u 1%) | 199 µM (f_u 85%) |
| peak ductal COX inhibition | 91.8% | 92.0% | 94.3% |
| AUC₀-∞ | 113 mg·h/L | — | — |

The three drugs reach comparable ductal COX inhibition from unbound
concentrations spanning **four orders of magnitude**, which is the entire point
of writing the enzymology rather than the dose.

One diagnostic is worth reporting because it looks like an error and is not.
Ibuprofen's instantaneous half-life at birth is **36.3 h**, inside the reported
20–43 h envelope — but the *apparent* terminal half-life measured from a
log-linear slope across days 4–10 is **16.7 h**. Clearance matures ~20%/day, so
the fitted slope is steepened by maturation occurring inside the sampling
window. Published preterm ibuprofen half-lives estimated without a maturation
term are biased the same way, which is one reason the reported ranges are so
wide and mutually inconsistent.

### 4.2 The same regimen is a smaller dose at day 7 than at day 2

| regimen | AUC over the 6 d after dose 1 | peak ductal COX inhibition |
|---------|-------------------------------|----------------------------|
| 10-5-5 from day 2 | 112 mg·h/L | 91.8% |
| 10-5-5 from day 7 | **74 mg·h/L (0.66×)** | 91.2% |
| 20-10-10 from day 7 | 148 mg·h/L (1.32×) | — |

Clearance maturation alone removes a third of the exposure by the second week.
This is the model's account of why late rescue needs a higher dose — and note
that peak COX inhibition barely moves, so the *reason* a fixed regimen fails
later is not loss of peak effect.

### 4.3 Spontaneous closure, and the discontinuity at the structural threshold

| GA | modelled closure day | observed (Semberova 2017) | duct at day 3 | shunt at day 3 | Q_p:Q_s |
|----|---------------------|---------------------------|---------------|----------------|---------|
| 24 wk | 57 d | 71 d (<26+0 group) | 1.89 mm | 124 mL/min/kg | 1.69 |
| 25 wk | 27 d | 71 d (<26+0 group) | 2.02 mm | 139 | 1.77 |
| 26 wk | **17 d** | **13 d** (26+0–27+6) | 2.15 mm | 151 | 1.84 |
| 27 wk | 12 d | 13 d (26+0–27+6) | 2.25 mm | 158 | 1.88 |
| 28 wk | 8 d | 8 d (28+0–29+6) | 2.19 mm | 150 | 1.83 |
| 30 wk | 4 d | 6 d (≥30 wk) | 0.87 mm | 10 | 1.05 |
| 32 wk | 2 d | — | 0.19 mm | 0 | 1.00 |

The 26–30 week range is reproduced closely. The **<26-week medians are
under-predicted** (57 d and 27 d against an observed group median of 71 d) — see
§9; the cause is identified and left unfixed.

What the model does get is the *shape*: a sharp break between 25 and 26 weeks,
sitting exactly where §2.3's occlusion ceiling crosses the closure threshold.
Above it, tone alone closes the duct in days. Below it, closure requires a
neointimal cushion the immature duct is barely competent to build, and takes
weeks.

### 4.4 Sixteen scenarios

Closure is counted only if it occurs within 8 days of the first dose —
otherwise the table would credit each drug with the spontaneous closure that was
going to happen at day 17 anyway. All deltas are against a matched no-drug
control.

| scenario (26 wk unless stated) | closure | spont. | reopens | Δburden | ΔSCr | Δ(death/BPD) |
|---|---|---|---|---|---|---|
| S1 expectant management | — | d17 | — | 0 | 0.00 | 0.0 |
| S2 early ibuprofen d2 (Baby-OSCAR) | d4.4 | d17 | **yes** | −227 | +0.05 | −2.5 |
| S3 late ibuprofen d7, standard | d7.6 | d17 | — | **−9** | +0.03 | **+15.3** |
| S4 late ibuprofen d7, high dose | d7.6 | d17 | — | **−9** | +0.06 | **+21.3** |
| S5 continuous ibuprofen infusion d2 | d3.3 | d17 | yes | −222 | +0.09 | +6.2 |
| S6 indomethacin d2 | d4.3 | d17 | yes | −226 | +0.05 | −8.7 |
| S7 prophylactic indomethacin from 8 h | — | d17 | — | −231 | +0.05 | −6.6 |
| S8 acetaminophen d2 | d4.9 | d17 | yes | −164 | **+0.00** | −10.5 |
| S9 acetaminophen d2 + chorioamnionitis | **fails** | d17 | — | −107 | +0.00 | −9.1 |
| S10 ibuprofen d2 + chorioamnionitis | d4.4 | d17 | yes | −227 | +0.05 | −2.5 |
| S11 ibuprofen + acetaminophen d2 | d2.7 | d17 | yes | −228 | +0.06 | +0.2 |
| S12 ibuprofen d2 + second course d7 | d4.4 | d17 | yes | −227 | +0.05 | +8.0 |
| S13 indomethacin d2 + hydrocortisone | d4.3 | d17 | yes | −226 | +0.05 | −8.7 |
| S14 ibuprofen d2, **24 wk** | **fails** | d57 | — | −56 | +0.06 | +4.2 |
| S15 ibuprofen d2, **29 wk** | d2.5 | **d6** | — | −23 | +0.07 | **+15.0** |
| S16 targeted: high-dose at d10 | d10.5 | d17 | — | **−0** | +0.04 | **+22.7** |

### 4.5 The finding that contradicts the hypothesis this model was built to test

The model was built expecting to show that **targeted late treatment** — treat
only the ducts echocardiography identifies as high-burden — captures most of the
benefit at a fraction of the drug exposure. It says the opposite, for a reason
visible in one column.

Look at Δburden for S3, S4 and S16: **−9, −9 and −0**. Treating at day 7 or
day 10 reduces PDA burden essentially not at all. The reason is that burden
accrues almost entirely in the **first week**: the untreated duct is already
narrowing by day 5, and by day 7 its shunt has fallen below the significance
threshold, so there is no burden left to prevent. What late treatment still
delivers in full is the COX exposure — and so Δ(death/BPD) is **+15.3, +21.3 and
+22.7**: pure harm.

This inverts the clinical proposal. If burden matters, it must be reduced in the
first week, and by the time a duct can be *confirmed* high-burden on serial
echoes at day 10, the damage is done and only the toxicity remains. Targeting
therefore needs early **prediction**, not late confirmation — a much harder
task than the one the strategy literature usually poses. It is also a direct
model-level explanation of PDA-TOLERATE's null result at 6–14 days.

S15 makes the same point from the other side: a 29-week duct that would have
closed on its own in **6 days** is treated, closes 3.5 days sooner, saves 23
units of burden, and ends up **+15.0 points worse**. That is the competing-risk
ceiling of §1 expressed as a number.

### 4.6 Why the model needs a drug harm, and how big it has to be

Calibrated to Baby-OSCAR's population (23–28 wk, n=90 virtual subjects, mean GA
25.9 wk), the two arms match the trial exactly:

| | modelled | observed |
|---|---|---|
| expectant, death or moderate/severe BPD | **63.5%** | 63.5% |
| early ibuprofen, death or moderate/severe BPD | **69.2%** | 69.2% |
| PDA burden, expectant → ibuprofen | 241 → 117 (**−51%**) | — |
| cumulative systemic COX exposure | 0 → 2.09 fraction·days | — |

But the interesting number is the one from *before* the harm term was fitted.
With burden as the only mechanism, early ibuprofen was predicted to **improve**
the composite by **10.4 points** (63.5% → 53.1%). The trial found it **worsened**
by 5.7. So the model has to find **16.1 points** somewhere, and the honest
reading is that this quantity — not the fit — is the result:

> Either the burden→BPD link is far weaker than the value assumed here, or early
> COX inhibition carries a cost of roughly 16 percentage points of death-or-BPD
> that the burden mechanism does not see. **The trial data cannot distinguish
> these**, because `B_BUR` and `B_BPD_COX` enter the composite in opposite
> directions and only their difference is observed.

This model resolves the ambiguity by fixing `B_BUR` a priori and attributing the
entire gap to a COX-exposure harm on BPD hazard, whose mechanism
(prostaglandin-dependent alveolar septation, plus renal fluid retention into an
already wet lung) is plausible but whose magnitude is *inferred from the
discrepancy rather than measured*. A reader who believes burden matters less
should read the same 16 points as evidence for a smaller `B_BUR`. Both readings
are consistent with every trial cited here, and saying so is more useful than
picking one silently.

One decomposition does **not** match and is reported rather than tuned away: the
model routes the whole harm through BPD, giving death 15.7% (ibuprofen) vs 17.0%
(expectant) — the *wrong direction*, since Baby-OSCAR observed 13.6% vs 10.3%.
The composite is right and its split between death and BPD is not.

### 4.7 Acetaminophen is peroxide-sensitive; ibuprofen is not

| | peak ductal COX inhibition | peak tone | attributable closure |
|---|---|---|---|
| acetaminophen, no sepsis | 94.3% | 0.89 | **yes, d4.9** |
| acetaminophen + chorioamnionitis | 88.8% | 0.80 | **NO** |
| ibuprofen, no sepsis | 91.8% | 0.89 | yes, d4.4 |
| ibuprofen + chorioamnionitis | **91.8%** | **0.89** | yes, d4.4 |

Raising ductal peroxide tone costs acetaminophen 5.5 points of COX inhibition
and its entire closure effect, while ibuprofen does not move by a thousandth.
This is the model's cleanest falsifiable prediction and it comes purely from
where the two drug classes bind — nothing here was fitted to any
inflammation-stratified dataset.

The honest caveat: closure near the threshold is knife-edge. At `IC50_APAP` = 8 µM
acetaminophen fails to close even without sepsis; at 6 µM it closes. The
*direction* of the peroxide effect is structural, but the *width* of the window
in which it flips closure is not well identified.

### 4.8 Reopening — a prediction, with an unexpected corollary

Nothing below was fitted to reopening data.

| GA | wall index | wall PO₂ at peak tone | neointima at d20 | reopens |
|----|-----------|----------------------|------------------|---------|
| 24 wk | 0.60 | **11.1 mmHg** | 0.00 | no (never closed) |
| 26 wk | 0.75 | 4.6 | **0.73** | **yes** |
| 28 wk | 0.90 | 4.1 | 0.98 | no |
| 30 wk | 1.05 | 6.2 | 0.99 | no |
| 34 wk | 1.35 | 5.2 | 0.99 | no |
| 38 wk | 1.65 | 4.5 | 1.00 | no |

The reopening phenotype appears exactly where clinicians see it — the treated
duct that constricts, then comes back — and disappears with advancing
gestation. The corollary the mechanism forces is testable and, as far as we can
tell, untested: because the ductal wall thickens postnatally, **the same duct
constricted late remodels where it would not have constricted early**, so early
treatment should reopen more often than late treatment at matched gestational
age.

### 4.9 The two NSAIDs diverge in the organs, not at the duct

| | cerebral flow nadir | mesenteric nadir | GFR nadir | urine nadir | bleeding-time index |
|---|---|---|---|---|---|
| ibuprofen | 78% | 70% | 87% | 2.2 mL/kg/h | 1.61 |
| indomethacin | **63%** | **57%** | 84% | 2.2 | 1.54 |
| acetaminophen | 79% | 71% | **100%** | 2.6 | **1.00** |

At essentially identical ductal COX inhibition (91.8 vs 92.0%), indomethacin
takes 15 more points off cerebral flow and 13 off mesenteric flow. That
separation is one parameter — the non-COX vasoconstriction of §2 — and it
reproduces the direction of the head-to-head Doppler trials. Acetaminophen
spares the kidney entirely (GFR nadir 100%, ΔSCr +0.00) and leaves platelet
COX-1 untouched, which is why it is the agent of choice when renal function or
bleeding is the constraint.

What the model does **not** support is a ranking of the two NSAIDs on the
composite endpoint. Their relative BPD hazard here depends on the ratio of each
drug's tissue Kᵢ values, and those ratios are not independently identified.

---

## 5. Calibration — what was fitted to what

Nine parameters were fitted, to eleven targets, in three separable stages. Every
other number in the model is fixed a priori from the literature.

| stage | fitted | to | result |
|-------|--------|----|--------|
| 1. ductal physiology | `NET50`, `NETW`, `TAUSYN0`, `KTAUSYN`, `KINVGA`, `KEP4`, `KREMOD` | 8 gestational spontaneous-closure times **and** 3 treated closure times, jointly | §4.3 |
| 2. drug potency | `KI_IBU`, `KI_IND`, `IC50_APAP` | comparable day-7 closure rates across the three drugs | §4.1, §4.7 |
| 3. outcome hazard | `H_BPD0` (level), `B_BPD_COX` (arm separation) | Baby-OSCAR, both arms | §4.6 |

Two points about stage 1 matter more than the numbers. First, untreated and
treated closure had to be fitted **jointly**: fitting the threshold to
spontaneous closure alone put it beyond the reach of any drug (§6, defect 4).
Second, `B_BUR` was deliberately **not** fitted, which is what makes §4.6's
16-point gap an interpretable quantity instead of an artefact of a free
parameter.

Nothing below was fitted, and nothing was adjusted after these were computed:

- reopening and its gestational dependence, plus the early-versus-late corollary (§4.8)
- the loss of acetaminophen effect under high peroxide tone (§4.7)
- the failure of late and "targeted" treatment to reduce burden at all (§4.5)
- BeNeDuctus, PDA-TOLERATE, TIPP
- the renal / cerebral / mesenteric separation of the two NSAIDs (§4.9)
- the entire 24-week structural failure (§4.4, S14)
---

## 6. Verification — two implementations, and the bugs the comparison found

Every equation is written twice: once in `pda_mrgsolve_model.R` (C++ via
mrgsolve) and once in `pda_reference_model.py` (SciPy/LSODA). A script diffs the
two parameter sets entry by entry; **176 shared parameters agree exactly**, and
the only two Python-only entries are the continuous-infusion rates, which
mrgsolve expresses through the event object's `rate` column instead.

Within the R model, all shared algebra lives in **one C++ macro**, `PDA_ALG(TT)`,
invoked from both `$ODE` and `$TABLE`. Writing that algebra twice — once for the
derivatives and once for the reported outputs — is the classic way a model this
size silently drifts, so that what gets reported is no longer what was
integrated. A macro cannot drift.

The comparison and the calibration together found **six real defects**, all of
which are documented in the source at the point where they were fixed:

1. **Ductal resistance gain 56× too large.** Qp:Qs was capped at 1.44 no matter
   how large the duct, so every shunt-driven output silently vanished — burden
   was identically zero in all 16 scenarios.
2. **A lung-water → PVR feedback strong enough to self-limit the shunt.** The
   duct throttled its own flow through pulmonary oedema.
3. **The product-form tone equation** of §2.3, which made preterm closure
   arithmetically impossible at any dose.
4. **A threshold fitted to untreated data alone.** `NET50` was first fitted to
   spontaneous-closure times only. Because most of the untreated rise in net
   drive came from contractile maturation and loss of NO — two terms *no drug
   touches* — the fit placed the threshold above anything PGE₂ suppression could
   reach, and every drug in the model became inert. The fix was to fit untreated
   **and** treated closure jointly, and to recognise that the dominant driver of
   postnatal constriction is the decline of the PGE₂/EP4 drive, which is exactly
   the axis the drugs act on.
5. **A wall-O₂ formulation with the right sign and the wrong magnitude.** A
   single supply divided by thickness gave only a 1.8-fold preterm-to-term
   hypoxia ratio and produced no reopening at any gestation, making §2.2 a claim
   the equations did not actually support.
6. **ΔSCr reported as a difference of window maxima.** Neonatal creatinine falls
   monotonically from the maternal value, so both arms' maxima are just the birth
   value and the difference is identically zero however nephrotoxic the drug is.
   The drug's real signature is a *blunted fall*, which needs an elementwise
   difference against a matched control.

Two errors of recollection were also caught by fetching the source abstracts
rather than trusting memory: Baby-OSCAR's primary outcome is **69.2%**, not
69.4%; and — far more consequential — Semberova's 26-week median is **13 days**,
not the ~48 days a smooth curve through the gestational ages would suggest. The
second correction turned out to *support* the architecture, because the 5.5-fold
jump it revealed sits exactly at the model's structural occlusion threshold.

---

## 7. Running it

```bash
# Mechanistic map (192 nodes, 21 clusters, 252 edges)
dot -Tsvg pda_qsp_model.dot -o pda_qsp_model.svg
dot -Tpng -Gdpi=150 pda_qsp_model.dot -o pda_qsp_model.png

# Python reference: full calibration report
python3 pda_reference_model.py            # all scenarios + virtual trials
python3 pda_reference_model.py --quick     # skip the population runs

# mrgsolve model: calibration tables
Rscript pda_mrgsolve_model.R

# Shiny dashboard (9 tabs)
Rscript -e 'shiny::runApp("pda_shiny_app.R")'
```

---

## 8. What would falsify this model most cheaply

The claims are ordered by how little work it would take to kill them.

1. **Acetaminophen's peroxide sensitivity.** Stratify any acetaminophen PDA
   cohort by chorioamnionitis or CRP and compare closure rates, with an
   ibuprofen arm as the control for severity confounding. The model says
   acetaminophen's ductal COX inhibition falls substantially while ibuprofen's
   does not move at all. If acetaminophen holds up in inflamed infants, §2.1 is
   wrong.
2. **Early-versus-late reopening.** The wall-thickening mechanism predicts
   reopening is *more* common after early treatment than after late treatment at
   the same gestational age. This is extractable from existing cohorts.
3. **The structural ceiling.** The model says no exposure closes a 24-week duct
   because the achievable occlusion leaves 0.44 mm. A dose-escalation study
   below 25 weeks that achieves closure in a substantial fraction refutes it —
   and note the 25-week row above sits only 0.02 mm above threshold, so the
   model is deliberately fragile exactly where clinical practice is uncertain.
4. **Burden versus closure as the outcome driver.** If PDA burden is computed
   from serial echoes and does *not* outperform binary closure status in
   predicting BPD, §2.4 is wrong.

---

## 9. Limitations

- **The <26-week spontaneous-closure median is under-predicted.** The model gives
  a pooled ~40 days against Semberova's 71. The cause is identified and left
  unfixed: the 25-week duct sits only 0.02 mm above the closure threshold, so it
  is hypersensitive to a small amount of remodelling and closes far sooner than
  a 24-week duct. The model reproduces the *existence and location* of the
  gestational discontinuity but not its full depth.
- **The drug-harm term is inferred, not measured.** §2.4's burden mechanism
  alone makes early ibuprofen beneficial; Baby-OSCAR found it slightly harmful.
  The model resolves this with an explicit COX-exposure harm on BPD hazard whose
  *mechanism* (prostaglandin-dependent alveolar septation, plus renal
  fluid retention) is plausible but whose *magnitude* is fitted to the arm
  separation. The useful output is therefore the trade-off size, not the fit.
- Hazards are deterministic integrals over a virtual population, not a
  time-to-event likelihood on patient-level data.
- Ligation and transcatheter closure appear on the mechanistic map but are not
  implemented as ODE interventions.
- Ibuprofen PK in preterm infants is genuinely under-determined in the
  literature: the commonly cited CL, Vd and t½ triplets are mutually
  inconsistent. A self-consistent set inside the reported envelope was adopted,
  and because only unbound drug reaches COX, the absolute total-concentration
  scale is absorbed into the fitted Kᵢ.
- No enantioselectivity (ibuprofen is given as a racemate with in-vivo R→S
  inversion), no genotype covariates, no CYP2C9 variation.

---

## ⚠️ Disclaimer

교육 및 연구 목적의 QSP 모델입니다. 공개 문헌과 임상시험 데이터로 구성했으나
독립적으로 검증·인증되지 않았으며, **실제 임상 의사결정, 처방, 규제 제출에
사용해서는 안 됩니다.** 특히 미숙아 동맥관 개존증은 "치료해야 하는가" 자체가
미해결 문제이며, 이 모델은 그 불확실성을 해소하는 것이 아니라 정량적으로
서술하기 위한 도구입니다.

This is an educational and research model. It is built from public literature
and trial data but has not been independently validated or certified, and **must
not be used for clinical decisions, prescribing, or regulatory submission.**
Whether preterm PDA should be treated at all is an open question; this model
exists to describe that uncertainty quantitatively, not to resolve it.
