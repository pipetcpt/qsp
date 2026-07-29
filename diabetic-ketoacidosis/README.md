# 당뇨병성 케톤산증 · 고혈당 고삼투압 상태 QSP 모델
# Diabetic Ketoacidosis & Hyperglycaemic Hyperosmolar State — QSP Model

> **한 줄 요약.** 케톤산증을 "인슐린이 없어 혈당이 오르고 케톤이 생긴다"로 놓지
> 않고, **배출 경로의 수가 다른 두 개의 결핍**으로 놓았다. 포도당에는 출구가
> 둘(인슐린 매개 처리 + 신장 배출)이고 그 신장 출구의 컨덕턴스는 **세포외액량**
> 이다. 케토산에는 출구가 하나이고 그것은 **포화**한다. 우리가 재는 혈당은
> 상당 부분 용적(volume) 신호이고, 우리가 잘 재지 않는 케톤이 곧 질병이다.
>
> **One-line thesis.** DKA is posed here not as "insulin is missing, so glucose
> rises and ketones appear", but as **two deficits with different numbers of
> exits**. Glucose has two (insulin-mediated disposal and a renal escape valve
> whose conductance *is* the extracellular volume); ketoacid has one, and it
> saturates. The glucose we measure is largely a volume signal; the ketone we
> usually do not measure is the disease.

---

## 파일 (Files)

| 파일 | 내용 |
|------|------|
| [`dka_qsp_model.dot`](dka_qsp_model.dot) | 기계론적 지도 소스 — 218 노드, 297 에지, 19 클러스터 |
| [`dka_qsp_model.svg`](dka_qsp_model.svg) · [`dka_qsp_model.png`](dka_qsp_model.png) | 렌더링된 지도 (SVG 벡터 / PNG 150 dpi) |
| [`dka_mrgsolve_model.R`](dka_mrgsolve_model.R) | mrgsolve 모델 — 42 ODE, 폐루프 프로토콜 실행기, 12 시나리오 |
| [`dka_shiny_app.R`](dka_shiny_app.R) | Shiny 대시보드 — 11 탭 |
| [`dka_references.md`](dka_references.md) | 문헌 132편, 섹션별 분류, 각 항목이 지지하는 모델 부분 표기 |
| [`dka_reference_check.py`](dka_reference_check.py) | 동일한 42개 방정식의 **독립 재구현** (순수 표준 라이브러리 Python, RK4) |
| [`dka_reference_output.txt`](dka_reference_output.txt) | 위 스크립트의 전체 출력 — **아래 모든 수치의 출처** |

아래에 인용된 숫자는 전부 `dka_reference_check.py`가 계산한 것이며, 그 전체
출력이 `dka_reference_output.txt`에 커밋되어 있습니다. 주장하지 않고 계산했습니다.
Every number below was computed by `dka_reference_check.py`, whose complete
output is committed. Nothing here is asserted.

---

## 구조적으로 특이한 세 가지 선택 (Three structural choices)

### 1. 중탄산염은 상태변수가 아니다 — 나머지(residual)다

Plasma [HCO₃⁻] is never integrated. At every derivative evaluation the model
solves the physicochemical (Stewart) electroneutrality condition

```
[HCO3-](pH, PCO2) + [A-](pH) + [ketoanion] + [lactate] + [SIG]  =  SID_app
SID_app = Na + K + (Ca,Mg) - Cl - (unmetabolised infused organic anion)
```

for pH by bisection, then reads off the bicarbonate. This costs work per step and
buys the following **for free, with no book-keeping term anywhere in the code**:

- making a ketoacid from a neutral triglyceride adds a strong anion → HCO₃⁻ falls
  1:1 with the ketoanion;
- **oxidising** a ketoanion removes it → HCO₃⁻ is regenerated 1:1, so retained
  ketoanions are literally *potential bicarbonate* (437 mmol at presentation);
- excreting a ketoanion **with Na⁺/K⁺** removes the anion *and* a strong cation →
  the anion gap closes but the bicarbonate does not recover: an organic acidosis
  is silently converted into a hyperchloraemic one;
- excreting it **with NH₄⁺** removes only the anion → base is preserved, and
  NH₄Cl excretion (which also removes Cl⁻) is *how the kidney repairs* a
  hyperchloraemic acidosis;
- 0.9% saline (SID 0) diluting a plasma SID of ~49 is acidifying, automatically;
- acetate/gluconate/lactate in balanced crystalloids count as strong anions until
  metabolised, then vanish — transient dip, then base gain.

### 2. 내인성 인슐린은 간에서 특혜를 받는다 (portal privilege)

Secreted insulin passes the liver first (~50% first-pass extraction), so portal
concentration is several-fold peripheral. The CPT-1 gate therefore sees
`PORTF·20·SEC` more insulin than the adipocyte does; **exogenous insulin has no
such privilege**. This single asymmetry is what makes HHS a different disease from
DKA in the same equations.

### 3. 내원 시점의 검사값은 입력이 아니라 출력이다

No laboratory value is typed in anywhere. Every scenario integrates a prodrome
from a healthy steady state with insulin withdrawn, and treatment starts from
wherever that lands. The 24-hour prodrome of a 70 kg adult with type 1 diabetes
and a moderate precipitating illness comes out as:

| | computed | typical reported |
|---|---|---|
| glucose | **470 mg/dL** (26.1 mmol/L) | 350–600 |
| pH | **7.132** | 7.00–7.24 (moderate–severe) |
| bicarbonate | **7.5 mmol/L** | 5–12 |
| anion gap | **28.7 mEq/L** | 25–35 |
| β-hydroxybutyrate | **14.8 mmol/L** | 10–20 |
| BHB : acetoacetate | **8.6** | 6–10 |
| Na measured / corrected | **137 / 146** | 130–138 / 142–150 |
| potassium | **5.10 mmol/L** | 4.5–5.8 |
| chloride | **101 mmol/L** | 95–105 |
| PCO₂ | **23 mmHg** | 18–28 |
| effective osmolality | **301 mOsm/kg** | 295–320 |
| BUN | **31 mg/dL** | 25–40 |
| creatinine (true / Jaffe assay) | **1.65 / 1.80 mg/dL** | the artefact is real |
| GFR | **52 mL/min** | prerenal |
| water deficit | **5.5 L = 7.9% of body weight** | ~100 mL/kg |
| **potassium deficit** | **301 mmol = 4.3 mmol/kg** | consensus 3–5 mmol/kg |
| cumulative urine | **7.5 L**, containing 192 g of glucose | — |

The classical deficit table is therefore *reproduced*, not entered.

---

## 계산된 결과 (Results — computed, not asserted)

### R1. 우리가 재는 지표들은 정보량의 역순으로 정상화된다

Standard protocol (15 mL/kg then 250 mL/h of 0.9% saline, insulin 0.1 U/kg/h,
KCl 40 mmol/L, dextrose titrated by a closed loop on the simulated glucose):

| resolution criterion | crosses at |
|---|---|
| glucose < 250 mg/dL (ADA target) | **2.9 h** |
| pH > 7.30 (ADA) | 3.7 h |
| bicarbonate ≥ 15 (ADA) | 4.8 h |
| bicarbonate ≥ 18 | 7.1 h |
| anion gap ≤ 12 (ADA) | **10.1 h** |
| **BHB < 0.6 mmol/L (JBDS — the actual disease)** | **10.7 h** |

Glucose is at target **7.8 hours before the ketosis is gone** — a factor of 3.7 in
time. First-hour rates: glucose −75 mg/dL/h, BHB −2.76 mmol/L/h.

### R2. 첫 한 시간의 혈당 강하는 어디서 오는가 — 두 개의 레버는 공유되지 않는다

| arm | Δglucose 1 h | Δglucose 2 h | ΔBHB 1 h | ΔBHB 2 h |
|---|---|---|---|---|
| fluid only, no insulin | **−37** | −61 | **−0.38** | −0.44 |
| insulin only, no fluid | −36 | −101 | **−2.30** | −5.32 |
| both (standard) | −75 | −156 | −2.76 | −5.80 |

Fluid alone lowers glucose by 37 mg/dL in an hour while barely touching the
ketone, because saline **reopens the renal escape valve**: GFR 52 → 64 → 83
mL/min over three hours, renal glucose clearance 1.50 → 1.74 L/h, urinary loss
39 mmol/h at a time when hepatic output is 105 mmol/h and total disposal 63. A
large part of the early glucose fall is renal, not metabolic — which is exactly
why the glucose is a poor guide to the disease.

### R3. 저용량 인슐린이 통하는 이유 — 두 IC50의 순서만으로 도출됨

`IC50` for suppression of lipolysis is 15 µU/mL; `EC50` for glucose disposal is
60. Nothing else in the model distinguishes the two arms.

| insulin | steady state | lipolysis suppression | disposal stimulation | t BHB<0.6 | t glucose<250 | K nadir | glucose nadir |
|---|---|---|---|---|---|---|---|
| 0.025 U/kg/h | 35 µU/mL | 0.345 | 0.130 | *not by 40 h* | *never* | 4.08 | 250 |
| 0.05 | 70 | 0.501 | 0.230 | *not by 40 h* | 4.1 h | 3.75 | 243 |
| **0.10** | 140 | 0.649 | 0.375 | **10.7 h** | 2.9 h | 3.28 | 207 |
| 0.14 | 196 | 0.709 | 0.456 | 9.5 h | 2.5 h | 3.24 | 204 |
| 0.20 | 280 | 0.761 | 0.545 | 9.4 h | 2.2 h | 3.20 | 201 |
| 0.40 | 560 | 0.833 | 0.705 | **8.9 h** | 1.9 h | 3.00 | 192 |

Across the 4-fold range 0.10 → 0.40 U/kg/h the **antiketogenic arm moves only
0.649 → 0.833** (nearly saturated), so resolution is **17% faster for four times
the insulin**; the **disposal arm moves 0.375 → 0.705** and is not saturated, so
the glucose fall, the glucose nadir and the potassium nadir all scale with dose.
Below the recommended range the lipolytic arm genuinely lapses (0.025–0.05
U/kg/h never clear the ketosis) — the other half of the same result. This is
Kitabchi's 1976 finding **derived from two half-maximal constants in the ratio
15 : 60**, rather than assumed.

### R4. 고혈당은 용적(volume) 질환이다 — 신장 탈출 밸브

Identical insulin deficiency; only the capacity to drink varies. 36 h prodrome:

| water access | glucose | eff. osmolality | GFR | urinary glucose | BHB | GCS |
|---|---|---|---|---|---|---|
| 1.00 (intact) | 540 mg/dL | 303 | 43 mL/min | 37 mmol/h | 16.3 | 11.6 |
| 0.50 | 637 | 332 | 31 | 30 | 17.0 | 4.4 |
| 0.25 | 788 | 366 | 20 | 21 | 17.6 | 3.0 |
| 0.05 | 882 | 385 | 16 | 17 | 17.8 | 3.0 |
| 0.00 (none) | 902 | 389 | 16 | 16 | 17.8 | 3.0 |

Nothing in the model knows the word "HHS". Glucose is set by how much of the
filtered load the kidney is still permitted to throw away, and the ketone barely
changes across the whole range — the two axes are independent.

### R5. DKA와 HHS는 같은 방정식, 두 개의 노브

| phenotype | knobs changed | glucose | pH | HCO₃ | gap | BHB | osm | Na corr | GCS |
|---|---|---|---|---|---|---|---|---|---|
| classic DKA (type 1) | β 0, water 1.0 | 540 | 7.072 | 6.9 | 30.8 | 16.3 | 303 | 147 | 11.6 |
| DKA, poor intake | β 0, water 0.25 | 788 | 7.068 | 9.7 | 34.4 | 17.6 | 366 | 178 | 3.0 |
| mixed DKA/HHS | β 0.10, water 0.25 | 499 | 7.340 | 18.0 | 19.6 | 4.7 | 317 | 154 | 13.1 |
| **HHS, early** | β 0.26, water 0.08 | 549 | 7.375 | 23.0 | 18.8 | **2.3** | 351 | 171 | 3.1 |
| **HHS, advanced** | β 0.22, water 0.06 | **694** | **7.376** | 23.2 | 20.4 | **3.0** | **372** | 181 | 3.0 |
| euglycaemic (SGLT2i) | SGLT2 on | **289** | 7.346 | 17.7 | 18.6 | 4.1 | 296 | 144 | 15.0 |

Residual β-cell function restrains **lipolysis** (IC50 15) and, through portal
privilege, the **CPT-1 gate**, long before it restrains glucose (EC50 60). That
ordering is the whole difference between a ketoacidotic and a hyperosmolar
presentation.

### R6. 정상혈당 케톤산증 — 밸브를 약으로 열어둔 상태

| | glucose | BHB | anion gap | pH | urinary glucose |
|---|---|---|---|---|---|
| no SGLT2 inhibitor | 422 mg/dL | 4.3 | 17.9 | 7.343 | 31 mmol/h |
| SGLT2 inhibitor on board | **289 mg/dL** | **4.1** | 18.6 | 7.346 | 36 mmol/h |

The **same ketone burden at a glucose 133 mg/dL lower**. A glucose-triggered
diagnostic rule (>250 mg/dL) fires at 5.0 h without the drug and at **22.4 h**
with it, while the anion gap has already crossed 16 at **7.3 h**. The diagnostic
delay is 15 hours, and it is a property of the trigger, not of the disease.

### R7. 중탄산염이 음이온차보다 늦게 회복되는 이유 — 원인이 두 개이고 분리 가능하다

| fluid | t gap<12 | HCO₃ at gap closure | Cl at 12 h | HCO₃ at 12 h | HCO₃ at 24 h | pH 24 h |
|---|---|---|---|---|---|---|
| 0.9% saline | 10.1 h | 19.8 | **113** | 20.1 | **20.3** | 7.368 |
| Plasma-Lyte 148 | 11.2 h | 22.9 | 107 | 22.9 | **23.8** | 7.392 |
| Ringer's lactate | 12.1 h | 21.7 | 108 | 21.7 | 23.5 | 7.390 |
| 0.9% saline + AKI | 4.8 h | 15.1 | 112 | 20.1 | 21.3 | 7.377 |
| Plasma-Lyte + AKI | 5.9 h | 21.0 | 108 | 24.0 | 24.5 | 7.396 |

The ledger: **102 mmol** of ketoanion had already been excreted with Na⁺/K⁺ before
arrival — over a 35 L bicarbonate space that is 2.9 mmol/L of base that cannot be
recovered by oxidising anything, because the carbon carrying it has left the body.
Against that, **625 mmol of chloride** is infused during treatment. Removing the
chloride load alone is worth **+3.6 mmol/L** of bicarbonate at 24 h.

And a result that inverts the intuition: **impose AKI and the anion gap closes
5 hours EARLIER**, because a ketoanion that cannot be filtered away has to be
*oxidised* — and oxidising it regenerates the bicarbonate that excreting it
loses. A failing kidney makes the gap look better and the base ledger worse.

### R8. 중탄산염 자체에 대한 정직한 경고 — 조건수(conditioning)가 나쁘다

At presentation the electroneutrality ledger reads:

```
strong ion difference             +48.79 mEq/L
weak acid (albumin, phosphate)    -13.45
ketoanion                         -16.55
lactate                            -4.43
other unmeasured strong anion      -6.82
--------------------------------  ------
                                   +7.53   = bicarbonate
```

The bicarbonate is **15% of the largest term in its own definition**. Perturb only
the fractional excretion of chloride by ±15%:

| FE(Cl) | Cl | HCO₃ | pH | **anion gap** | **BHB** |
|---|---|---|---|---|---|
| 0.0119 | 102.1 | 6.9 | 7.083 | 28.5 | 14.9 |
| 0.0140 | 101.0 | 7.5 | 7.132 | 28.7 | 14.8 |
| 0.0161 | 100.1 | 8.1 | 7.169 | 28.9 | 14.8 |

A 30% swing in one renal parameter moves the bicarbonate by 1.2 mmol/L and the pH
by 0.086, while the gap and the ketone — which involve no cancellation — move by
1.4% and 0.7%. **This is a property of the chemistry, not of the model.** It is
why the anion gap and the point-of-care ketone are the robust bedside variables,
why "bicarbonate recovery" is a poor endpoint once chloride is being infused, and
why *this model's bicarbonate and pH predictions deserve less trust than its gap
and ketone predictions*. Stated as a warning label rather than buried.

### R9. 칼륨 — 내원 시 혈청값은 결핍량에 대해 거의 정보가 없다

| presentation | serum K | pH | eff. osm | **K deficit** | mmol/kg |
|---|---|---|---|---|---|
| standard, 24 h | 5.10 | 7.132 | 301 | 301 mmol | 4.3 |
| milder illness | 5.06 | 7.162 | 300 | 286 | 4.1 |
| severe illness | 5.14 | 7.105 | 302 | 306 | 4.4 |
| hyperosmolar | 6.08 | 7.043 | 345 | 348 | 5.0 |
| **long prodrome, 40 h** | **5.29** | 7.072 | 303 | **509** | **7.3** |
| HHS-like | 5.82 | 7.376 | 372 | 303 | 4.3 |

Serum K spans 20% while the deficit spans 78%. Two patients here differ by
**223 mmol of body potassium at serum values 0.19 mmol/L apart**.

What the ADA "hold insulin if K < 3.3" rule buys, quantified:

| KCl in fluid | K nadir *with* rule | K nadir *without* | hours insulin held | t BHB<0.6 with rule | without |
|---|---|---|---|---|---|
| 0 mmol/L | 3.30 | **2.75** | 0.2 | 11.2 h | 8.4 h |
| 20 | 3.30 | 2.75 | 0.0 | 10.7 | 8.4 |
| 40 | 3.29 | 2.75 | 3.0 | 10.7 | 8.4 |

The rule trades **~2.3 hours of extra ketoacidosis for 0.55 mmol/L of potassium**.
That is the trade, and it is why DKA cannot be treated with insulin alone even in
a patient whose potassium looks high.

### R10. 뇌부종 — 모델이 음성 시험(PECARN FLUID 2018)을 재현한다

| fluid-rate arm | Vbr at presentation | peak | swelling | peak ICP | min GCS |
|---|---|---|---|---|---|
| slow (10 mL/kg + 125 mL/h) | 1.0055 | 1.0107 | 0.52% | 20 mmHg | 13.0 |
| standard (15 mL/kg + 250) | 1.0055 | 1.0106 | 0.51% | 20 | 12.9 |
| fast (20 mL/kg + 500) | 1.0055 | 1.0116 | 0.61% | 20 | 12.7 |
| very fast (30 mL/kg + 750) | 1.0055 | 1.0126 | 0.71% | 21 | 12.4 |

A **6-fold range of fluid rate moves peak brain volume by 0.2 percentage points**,
which is what the randomised trial found. Decomposing the swelling in a severe
presentation: **both mechanisms 0.47%, osmotic alone 0.00%, injury alone 0.94%**.
The osmotic rate-mismatch is the *smaller* contribution; the larger is an
ischaemic insult that accrues **before treatment starts**, driven by hypocapnia
and hypoperfusion — which is why the reported risk factors are a low PCO₂ and a
high urea *at presentation* rather than anything done afterwards.

### R11. 니트로프루시드 역설 — 환자가 좋아지는데 소변 케톤은 나빠진다

| h | BHB | AcAc | total | BHB:AcAc | hepatic redox | acetone |
|---|---|---|---|---|---|---|
| 0 | 14.83 | 1.72 | 16.55 | **8.60** | 2.87 | 0.34 |
| 2 | 9.03 | **2.01** | 11.04 | 4.50 | 1.49 | 0.36 |
| 4 | 5.00 | 1.39 | 6.39 | 3.60 | 1.20 | 0.38 |
| 12 | 0.47 | 0.16 | 0.63 | 3.03 | 1.01 | 0.35 |
| 24 | 0.32 | 0.11 | 0.43 | 3.00 | 1.00 | 0.26 |

Total ketone falls monotonically, but the acetoacetate **concentration rises** for
the first two hours and its fraction goes from 10.4% to 21.7% as hepatic redox
normalises. At 4 h the total ketone is down **61%** while the
nitroprusside-visible pool (AcAc + acetone) is down **14%**. Acetone, with a 23 h
half-life, is still 0.26 mmol/L at 24 h — days of ketotic breath and positive
strips after the disease has gone.

### R12. 중탄산염 투여 — 무엇을 사는가

Severe presentation (pH 7.069, HCO₃ 6.8, gap 30.5, BHB 16.1):

| arm | pH 1 h | pH 4 h | pH 12 h | pH 24 h | K nadir | PCO₂ 1 h |
|---|---|---|---|---|---|---|
| none | 7.073 | 7.293 | 7.366 | 7.368 | 3.30 | 25 |
| 100 mmol / 2 h | 7.206 | 7.345 | 7.376 | 7.372 | 3.30 | 25 |
| 200 mmol / 2 h | **7.306** | 7.393 | 7.394 | 7.379 | 3.29 | 27 |

Bicarbonate buys **+0.233 pH units at 1 h**, +0.029 at 12 h and +0.011 at 24 h.
The acid was never the reservoir — the ketoanion was, and insulin sets the rate at
which it is oxidised back into bicarbonate. The costs are computed too: a higher
PCO₂ (the respiratory drive is switched off by the very pH it corrected) and
200 mmol of extra sodium in a patient whose brain is already adapting to a falling
osmolality.

### R13. 정맥주입 중단 — 겹침(overlap)은 형식이 아니다

| transition | gap at 2 h | gap at 6 h | BHB at 6 h | glucose at 6 h |
|---|---|---|---|---|
| insulin stopped at gap closure, no s.c. | 11.7 | **15.4** | **5.52** | 425 |
| stopped, 10 U s.c. lispro given | 11.5 | 13.2 | 3.71 | 365 |
| stopped 2 h after gap closure, 10 U | 11.4 | 13.0 | 3.59 | 362 |

Plasma insulin has a **5-minute** half-life; the CPT-1 gate reopens within an hour
of stopping, and ketogenesis restarts from a fat store that treatment never
emptied. The gap *reopens* to 15.4 within six hours.

### R14. 피하 인슐린 — 흥미로운 숫자는 정점이 아니라 저점, 그리고 두 개의 저점이다

Mild-to-moderate presentation, 0.3 U/kg s.c. load then 0.2 U/kg every 2 h:

| route | t BHB<0.6 | t glucose<250 | K nadir | peak plasma | **trough plasma** | **trough effect site** |
|---|---|---|---|---|---|---|
| IV infusion | 9.1 h | 2.6 h | 3.28 | 70 µU/mL | 23 | 48 |
| s.c. q2h | 9.1 h | 2.2 h | 2.81 | 163 | **13** | **26** |

The **plasma** trough falls to 13 µU/mL — *below* the 15 µU/mL lipolysis IC50,
because plasma insulin has a 5-minute half-life and no bolus regimen can avoid a
trough. But the **effect site** only falls to 26 µU/mL (1.73× the IC50), because
it is a 21-minute low-pass filter on that profile. **That filter is the whole
reason an intermittent subcutaneous protocol can suppress ketogenesis at all** —
and it says exactly where it breaks: lengthen the interval, or slow absorption
(shock, oedema, regular insulin at kₐ 0.5 instead of 1.5 /h), and it is the
effect-site trough, not the peak, that crosses the constant.

### R15. 표현형 지도 — 같은 모델, 여러 질병

| phenotype | glucose | pH | HCO₃ | gap | BHB | BHB:AcAc | osm | GCS |
|---|---|---|---|---|---|---|---|---|
| adult T1D, moderate illness | 506 | 7.088 | 7.0 | 29.9 | 15.7 | 8.6 | 302 | 12.0 |
| adult T1D, severe sepsis | 510 | 7.069 | 6.8 | 30.5 | 16.1 | 8.6 | 303 | 11.6 |
| child 30 kg (allometric) | 424 | 7.207 | 9.0 | 27.5 | 13.2 | 4.1 | 299 | 15.0 |
| DKA in pregnancy | 469 | **6.985** | 5.8 | 31.0 | **20.6** | 10.2 | 304 | 9.7 |
| **alcoholic ketoacidosis** | **184** | 7.351 | 18.3 | 19.0 | 7.9 | **18.4** | 295 | 15.0 |
| euglycaemic (SGLT2i, T2D) | 289 | 7.345 | 17.6 | 18.7 | 4.1 | 7.3 | 296 | 15.0 |
| HHS, elderly type 2 | 681 | 7.376 | 23.1 | 20.2 | 3.0 | 6.9 | **370** | 3.0 |
| DKA with established AKI | **732** | 7.209 | 9.2 | 24.8 | 14.1 | 8.6 | 299 | 15.0 |

Alcoholic ketoacidosis comes out with **BHB:AcAc 18.4** (against 8.6 in DKA) at a
glucose of 184 mg/dL: the *redox* term, not the insulin term, is doing the work —
which is why a nitroprusside strip, blind to BHB, is at its most misleading in
exactly this patient.

### R16. 민감도 (elasticities of the presenting phenotype, ±25%)

| parameter | Δglucose | ΔBHB | ΔHCO₃ |
|---|---|---|---|
| LIPMAX (lipolytic capacity) | +58% | **+232%** | **−271%** |
| KGSCALE (ketoacid per NEFA) | +58% | +220% | −269% |
| VMAX_KOX (ketone oxidation) | −34% | −143% | +167% |
| HGP0 (hepatic output) | **+133%** | +8% | −49% |
| FENA (renal Na handling) | +39% | +9% | −88% |
| TMGLU (renal glucose Tm) | +35% | −2% | −1% |
| GFRMAX | −37% | +0.4% | −27% |
| KMAL (CPT-1 gate constant) | +6% | +44% | −52% |
| FECL (renal Cl handling) | −6% | −2% | **+44%** |
| IC50_LIP, EC50_UP, VMI_UP, PORTF | ~0 | ~0 | ~0 |

The last row is not a bug and is worth reading: in an **untreated** type 1 patient
plasma insulin is essentially zero, so the insulin IC50s and the portal ratio have
no leverage over the *presentation*. They have all their leverage over the
*response to treatment* (R3) and over the *HHS phenotype* (R5). The parameters
that set what walks through the door and the parameters that set what happens next
are almost disjoint.

### R17. 보존 검사 (conservation — these must hold or the model is wrong)

- electroneutrality residual, max over 24 h: **1.9 × 10⁻¹² mEq/L**
- anion gap reconstructed from its parts at 6 h: **14.573 vs 14.573**
- sodium: ECF content 1551 → 1819 mmol = infused 491 − urinary 223 (residual 0)
- potassium: total body **+6.7 mmol** vs (infused 133.6 − urinary 126.9) = +6.7;
  closure residual **3.8 × 10⁻¹² mmol** (0.000% of throughput)
- ketone: retained 437 → 12 mmol, 30 mmol excreted, the remaining **395 mmol
  oxidised — regenerating bicarbonate 1:1**, which is where the recovered base
  came from.

---

## 모델 구조 (Model structure)

**42 ODEs.** Volume and electrolytes (VECF, VICF, NAE, CLE, KE, KI, PHOSE, ORGA)
· substrates and ketones (GLU, GLYCO, FFA, LAC, ACAC, BHB, ACET) · nitrogen
(UREA, CREA) · respiratory (PCO₂) · insulin (INSSC, INSP, INSEF, INSES) ·
counter-regulatory (GCG, CORT, EPI, IR, REDOX) · renal (GFRR, ALDO, NH4C) ·
β-cell (BETAF) · illness (ILL) · brain (OSMB, VBR, INJ) · CNS (MENT) · six
cumulative ledgers (UKET, UKETN, UGLU, UKCUM, UVOL, CLIN) plus two mass-balance
audit trackers.

**The load-bearing parameters** are four half-maximal insulin concentrations taken
from dose-response literature rather than fitted — lipolysis 15, potassium shift
25, hepatic output 30, peripheral disposal 60 µU/mL — and the *ordering* of those
four numbers generates R1, R3, R5 and R14.

**Closed-loop protocols.** The dextrose switch, the insulin taper and the ADA
potassium rule are feedback on the simulated patient, not events on a clock, so
the protocol is an output of the model as much as an input to it. That is what
makes R9's "hours insulin held" a computed quantity.

**Allometry.** Flux, clearance and volume parameters scale with body weight; rate
constants, concentrations, IC50s and stoichiometries do not. The 30 kg child in
R15 is a 30 kg child, not a 70 kg adult with a relabelled weight.

---

## 알려진 한계와 편향 (Known limitations, recorded rather than tuned away)

1. **중탄산염의 조건수.** Documented quantitatively in R8. The gap and the ketone
   are the trustworthy outputs; the bicarbonate and pH are not, and the model says
   so with numbers rather than a disclaimer.
2. **뇌부종의 분해는 추론이다.** The 0.00%/0.94% split between osmotic and
   ischaemic mechanisms in R10 reproduces the reported risk factors and the
   negative fluid-rate trial, but the partition itself cannot be measured in
   patients. It should be read as "the osmotic story is insufficient", not as a
   calibrated ratio.
3. **문맥:말초 인슐린 비(PORTF)가 DKA/HHS 구분 전체를 지고 있다.** It is taken
   from fasting first-pass extraction studies; whether the same ratio holds during
   a hyperglycaemic crisis with hepatic congestion is unknown. This is the model's
   single largest structural bet.
4. **케톤 산화 용량의 포화/비포화 분할** (VMAX_KOX vs KLIN_KOX) is constrained by
   steady-state concentrations, not measured. The non-saturable arm is what keeps
   the system bounded, so this is not a cosmetic choice.
5. **0.05 U/kg/h에서 케톤산증이 40시간 내에 해소되지 않는다.** Reported as-is in
   R3. Whether real patients on half the recommended dose fail this badly is not
   established; the model is extrapolating below the range its constants were
   measured in.
6. **선행 질환은 단일 스칼라**(ILL, decaying with τ 30 h). There is no infection
   dynamics, no source control, no antibiotic effect. A patient whose sepsis is
   worsening cannot be represented.
7. **응고·심전기생리·약물유전학 없음.** Venous thromboembolism and arrhythmia
   appear on the mechanistic map as outcomes but are not dynamical states; the
   potassium panel reports a concentration, not a rhythm.
8. **첫 시간 혈당 강하가 다소 빠르다** (−75 mg/dL/h against the ADA's stated
   expectation of 50–75). The model sits at the upper edge of the guideline band.
   The cause is identifiable — renal clearance at a presenting GFR of 52 mL/min —
   and was left alone rather than tuned, because R2 depends on that flux being
   honest.

---

## 재현 (Reproducing this)

```bash
# the mechanistic map
dot -Tsvg dka_qsp_model.dot -o dka_qsp_model.svg
dot -Tpng -Gdpi=150 dka_qsp_model.dot -o dka_qsp_model.png

# every number in this README, from an independent implementation
python3 dka_reference_check.py            # ~6 min, no dependencies at all
python3 dka_reference_check.py --quick    # just the presentation table
```

```r
# the mrgsolve model
source("dka_mrgsolve_model.R")
pr <- present();  print(t(presentation_row(pr$sim)))
print(resolution_table(sc1_standard()))
print(lapply(sc10_phenotypes(), presentation_row))

# the dashboard
shiny::runApp("dka_shiny_app.R")
```

---

## ⚠️ 면책 (Disclaimer)

교육 및 연구 목적의 QSP 모델입니다. 개별 환자 데이터에 적합·검증되지 않았으며
실제 임상 의사결정, 처방, 규제 제출에 사용해서는 안 됩니다. 파라미터는 공개
문헌에 근거한 근사치이고, 위에 기록한 한계는 축소하지 않고 그대로 남겨 두었습니다.

Educational and research use only. This model has not been fitted or validated
against individual patient data and must not be used to guide treatment. The
limitations above are recorded at full strength rather than minimised.
