# 알라질 증후군 (Alagille Syndrome, ALGS) — QSP 모델

**JAG1/NOTCH2 반수부전 → 담관 형성 부전 → 만성 담즙 정체 → 소양증·황색종·성장 부전·담즙성 간경변**

| 산출물 | 파일 | 규모 |
|---|---|---|
| 기계론적 지도 | [`algs_qsp_model.dot`](algs_qsp_model.dot) · [SVG](algs_qsp_model.svg) · [PNG](algs_qsp_model.png) | 184 노드 · 251 엣지 · 22 클러스터 |
| mrgsolve ODE 모델 | [`algs_mrgsolve_model.R`](algs_mrgsolve_model.R) | 44 구획 · 16 시나리오 |
| Shiny 대시보드 | [`algs_shiny_app.R`](algs_shiny_app.R) | 12 탭 · 23 출력 (`testServer` 전수 통과) |
| 참고문헌 | [`algs_references.md`](algs_references.md) | 154 PMID (전부 PubMed API로 실제 확인) |

---

## The one structural fact everything follows from

ALGS cholestasis is not one quantity. It is the collision of two fluxes that
are set **independently** and that every clinical measurement mixes together:

- **`J_DUCT`** — how much bile the liver can physically push into the gut.
  Set developmentally by Notch dose, slowly repaired postnatally by ductular
  reaction, and moved by **no drug in this model**.
- **`R`** — how much of that bile comes straight back through ASBT. The
  **only** thing an IBAT inhibitor touches.

Write those two separately and seven consequences fall out of published
numbers rather than being asserted. Every figure below is printed by
`run_all_algs()`.

---

## Axis 1 — the drug's ceiling is a duct property, not a bile-acid property

Define **Φ_EHC = R / (S + R)**: the share of the hepatocyte's bile-acid input
that arrived back from the gut rather than being freshly synthesised. An ideal
IBAT inhibitor removes Φ_EHC of that input and not one percent more.

In an unobstructed liver Φ_EHC = **0.98**, because each molecule is recycled
~50 times per molecule synthesised. Duct paucity caps biliary output, which
caps `R` with it. Calibrated **only** to the ASSERT bile-acid ratio, the model
puts the average trial patient at

> **Φ_EHC = 0.42** (duct capacity `JCAP` = J_DUCT/S = 1.99)

Duct paucity has destroyed more than half the drug's target before the first
dose. Sweeping duct capacity moves the achievable 24-week response
monotonically from **−92%** (JCAP 6.1) to **−9%** (JCAP 0.38). The drug does
not fail because bile acids are too high. It fails because there are no ducts.

## Axis 2 — a hypothesis the model rejected, and the collision that replaced it

The design intent was that below some duct capacity an IBAT inhibitor would
**invert**: blocking ASBT removes ileal FGF19, derepresses CYP7A1, and the
newly synthesised bile acid would have no duct to leave by, so serum bile acid
would *rise* on treatment.

**The model says that does not happen** — for a reason worth more than the
hypothesis was. In ALGS the FGF19 signal is *already floored before treatment*
(baseline FGF19 = 0.04 of normal; synthesis already derepressed 2.14-fold), so
the synthesis reserve that would drive an inversion has already been spent.
Blockade in severe paucity is **futile, not harmful**.

What survives is a graded threshold, and it lands somewhere unexpected:

| model boundary | duct capacity | total bilirubin | GALA's independent cut-point |
|---|---|---|---|
| response falls below 30% | JCAP 1.25 | **5.42 mg/dL** | **5.0 mg/dL** |
| response falls below 15% (futile) | JCAP 0.60 | **10.70 mg/dL** | **10.0 mg/dL** |

GALA is a natural-history cohort keyed on native liver survival with no drug
in it. This model never saw those thresholds. The **ratio** of the two
boundaries is **1.97** here against GALA's **2.00** — and unlike the
boundaries themselves, the ratio is completely independent of the single
assumed calibration anchor (see `FAILURE 2`). Two datasets with nothing in
common appear to be measuring the same duct-capacity threshold from opposite
directions.

## Axis 3 — itch is not a bile-acid measurement, and ASSERT's placebo arm proves it arithmetically

In ASSERT the **placebo** arm's itch *fell* by 0.8 of the 1.7 points the drug
arm fell — 47% of the on-drug improvement is not drug — while the same placebo
arm's bile acids **rose** by 22 µmol/L. One axis has a large favourable
placebo response; the other has a negative one. They cannot be the same
variable. The consequence is quantitative:

| slope (itch points per µmol/L) | value | ×controlled |
|---|---|---|
| placebo-controlled (ASSERT) | 0.9/113 = 0.0080 | 1.00 |
| single-arm (ASSERT) | 1.7/90 = 0.0189 | 2.37 |
| single-arm (ICONIC) | 1.6/96 = 0.0167 | 2.09 |

The two **single-arm** slopes agree with each other to 13% and are both ~2.1×
the placebo-controlled truth. Any model calibrated on single-arm cholestasis
data over-attributes itch to bile acids about two-fold. This model carries the
components separately and reproduces the controlled slope (0.0083 vs 0.0080)
and the drug-attributable itch difference (**−0.90 vs −0.90**) exactly.

A psychophysical note that changed the pharmacology: fitting ASSERT **rejected
a logarithmic (Weber–Fechner) itch law.** Forced through a saturating 0–4 map
with a baseline at 2.8, the trial demands a Stevens exponent of **1.66** —
supralinear, which no itch psychophysics supports. Treating the 0–4 scale as a
ceiling on the *instrument* rather than on the *sensation* puts the operating
point in the near-linear range and the solved exponent falls to **0.95**. The
reported scale is not the latent variable, and pretending otherwise distorts
the drug effect.

## Axis 4 — the survival benefit is larger than its own bile-acid effect can explain

GALA supplies its own exposure–hazard gradient: relative to TB < 5.0 mg/dL the
transplant hazard is 4.8× at 5–10 and 15.6× at ≥10, i.e. a power law of
exponent **n = 1.5–1.85** (n = 1.85 reproduces the severe stratum best, 17.3
predicted vs 15.6 observed).

Running six years of maralixibat against this model's own natural history:

| exponent n | model HR | published HR | GALA TB≥10 back-predicted | observed |
|---|---|---|---|---|
| 1.60 | 0.445 | 0.305 | 11.8 | 15.6 |
| **1.85** | **0.414** | **0.305** | **17.3** | **15.6** |
| 2.40 | 0.356 | 0.305 | 40.3 | 15.6 |
| 2.90 | 0.314 | 0.305 | **87.1** | 15.6 |

**No single exponent fits both.** At the exponent that fits GALA's own strata,
the model predicts HR 0.414 where 0.305 was published; at the exponent that
reproduces 0.305, GALA's severe stratum is over-predicted **5.6-fold**. A
residual hazard ratio of 0.305/0.414 = **0.74 is not explained by bile acids.**

It gets worse, not better: gamma frailty (Axis 7) *attenuates* population
hazard ratios toward 1 over time, so the individual-level effect implied by an
observed population HR of 0.305 is larger still. The candidates — six-year
drug-persistence selection against a comparator matched on baseline labs only,
a bile-acid-independent benefit, or curvature the cross-sectional strata cannot
see — are **not separable by any published data**. The model ships with n = 1.6
so that it *under*-predicts the trial, and scenario 12 prints both errors side
by side.

## Axis 5 — where in the gut a drug acts decides whether it starves the child

ASBT sits in the **terminal ileum**, downstream of the duodenal micellar window
where fat and vitamins A/D/E/K are absorbed. At steady state a duct-limited
liver delivers `J_DUCT` to the duodenum whatever ASBT is doing. Two years, one
patient, severe paucity:

| | serum bile acid | fat absorption | vitamin D | INR | height z |
|---|---|---|---|---|---|
| untreated | 236 | 0.662 | 19.1 | 1.47 | −1.51 |
| **odevixibat** | −46% | **0.652** | **18.8** | 1.48 | **−1.42** (improves) |
| **cholestyramine** | −37% | **0.294** | **8.5** | **2.12** | **−2.38** |

Both drugs remove bile acids; only the proximal-acting one causes
steatorrhoea, and the contrast is derived from anatomy rather than asserted.
The model also finds **its own boundary case**: in *mild* paucity the duct is
not saturated, duodenal delivery does fall with the drug (6.2 → 4.3 mM), fat
absorption drops 0.844 → 0.752, and the contrast narrows. The safety argument
for IBAT inhibitors depends on the liver being duct-limited.

## Axis 6 — the population ceiling

Cardiac disease dominates ALGS mortality in year one and
vasculopathy/intracranial haemorrhage later; both are Notch-dose diseases no
cholestasis drug touches. To age 18 the liver carries **77.9%** of the fatal
hazard, cardiac 10.4%, vascular 11.7% — so a therapy achieving a
liver-specific HR of **0.386** moves all-cause hazard only to **0.522**.

## Axis 7 — the GALA curve is heterogeneity, not biology

No single trajectory fits GALA native liver survival: its hazard is strongly
front-loaded and every mechanistic trajectory here is not (best
single-trajectory fit SSE 3.3×10⁻², visibly wrong at age 5). Admitting a
**gamma frailty of variance 2.94** fits all three time points **65-fold
better** (SSE 5.0×10⁻⁴):

| age | model (population) | GALA |
|---|---|---|
| 5 yr | 67.7% | 66.8% |
| 10 yr | 52.6% | 54.4% |
| 18 yr | 41.3% | 40.3% |

A variance of 2.94 means the **standard deviation of individual hazard is 1.7×
its mean** — the notorious ALGS phenotype spread (same variant → infant
transplant or asymptomatic parent), measured rather than described.
Consequence: population curves must use the `_POP` outputs. Comparing the
individual-level `exp(-H)` outputs to a cohort study is a category error that
this file made once (reporting 1.5% 18-year native liver survival against a
true 40.3%) before it was caught.

---

## Calibration summary

| target | published | model |
|---|---|---|
| ASSERT sBA, drug ÷ placebo at wk 21–24 | 0.550 | **0.543** |
| ASSERT itch difference (drug − placebo) | −0.90 | **−0.896** |
| ASSERT sBA difference | −113 µmol/L | −108 µmol/L |
| ASSERT placebo-controlled itch slope | 0.0080 | **0.0083** |
| ICONIC week-48 sBA change | −96 µmol/L | −107 µmol/L |
| ICONIC week-48 ItchRO(Obs) reported | −1.6 | −1.49 |
| GALA native liver survival 5/10/18 yr | 66.8/54.4/40.3% | **67.7/52.6/41.3%** |
| GALA event-free 10/18 yr | 48.5/34.0% | 46.6/36.8% |
| daily-bolus vs continuous-input dosing | — | 4.0% on sBA, 5.7% on blockade |

---

## 시나리오 (16)

| # | 시나리오 | 목적 |
|---|---|---|
| 0 | 평형 검증 | 24주 시험 재현이 평형화 과도현상이 아님을 확인 |
| 1 | 자연 경과 (1→18세) | GALA 대비 보정 |
| 2 | ICONIC 재현 (maralixibat 380 µg/kg/d) | 단일군 기울기 |
| 2b | 일일 볼루스 vs 연속 투여 등가성 | 연속 입력 근사 검증 |
| 3 | ASSERT 재현 — **양쪽 군** (odevixibat 120 µg/kg/d) | Axis 3의 산술적 증명 |
| 4 | 담즙산당 소양증 기울기 | 대조군 유무에 따른 2.1배 차이 |
| 5 | 담관 용량 스윕 | Axis 2 — GALA 절단점과의 충돌 |
| 6 | Φ_EHC 식별 가능성 | FAILURE 1의 정량화 |
| 6b | 빌리루빈 앵커 민감도 | FAILURE 2의 정량화 |
| 7 | PEBD (외과적 담즙 전환) vs IBAT 억제제 | 완전 차단의 상한 |
| 8 | 약제 패널 (24주) | 9개 요법 동일 환자 비교 |
| 9 | 미셀 창 / 작용 부위 | Axis 5 및 그 경계 사례 |
| 10 | 조기 vs 지연 치료 시작 | 적분 논증 (노출 회피 39.7% → 17.4%) |
| 11 | 무사건 생존 (6년) | Axis 4의 핵심 검정 |
| 12 | 지수 모순 스윕 | 어떤 n도 양쪽을 만족하지 않음 |
| 13 | 경쟁 위험 상한 | Axis 6 |
| 14 | 비반응자 표현형 | DPR 0.32 vs 0.15 전체 궤적 |
| 15 | 출생 후 담관 재생 | FAILURE 4의 분해 |

```r
source("algs_mrgsolve_model.R")
run_all_algs()          # 전 시나리오 실행
shiny::runApp("algs_shiny_app.R")
```

렌더링:
```bash
dot -Tsvg algs_qsp_model.dot -o algs_qsp_model.svg
dot -Tpng -Gdpi=150 algs_qsp_model.dot -o algs_qsp_model.png
```

---

## 이 모델이 하지 못하는 것 (수리하지 않고 보고함)

1. **Φ_EHC는 단일 시험으로 개별 식별되지 않는다.** 담관 용량, 달성된 ASBT
   차단률, 합성 예비능이 정상상태 담즙산 비율에 하나의 조합으로만 들어간다.
   가정된 최대 차단률을 0.50→0.92로 바꾸며 ASSERT를 계속 정확히 재현하도록
   담관 용량을 재적합하면 JCAP은 3.58→1.95 (1.8배), Φ_EHC는 0.57→0.42
   (1.36배) 이동한다. Φ_EHC가 더 견고하지만 어느 쪽도 고정되지 않는다.
2. **빌리루빈 경계는 가정된 앵커 하나에 의존한다.** ASSERT 집단의 평균 기저
   총빌리루빈을 원논문 초록에서 찾지 못해 3.5 mg/dL로 두었다. 앵커
   3.0/3.5/4.0은 30% 경계를 4.58/5.42/6.26, 15% 경계를 8.93/10.70/12.48로
   옮긴다. GALA 절단점을 괄호로 감싸지만 점추정치는 증거가 아니다. **두
   경계의 비(1.97)만이 앵커 독립적**이며, 그래서 Axis 2는 비를 앞세운다.
3. **담즙산 비의존 소양증 분율은 풀린 것이 아니라 가정된 것이다.** ASSERT
   위약군은 *비약물* 성분(47%)을 정량할 뿐, 담즙산이 0이 되었을 때 남는
   소양증에 대해서는 아무 말도 하지 않는다. 완화 요인 하나: 시도한 모든
   Stevens 지수(0.75–1.05)와 가정 분율(8–18%) 조합에서 잔여치가 12–17%로
   나와, 결론은 파라미터보다 덜 민감하다. 그래도 가정이다.
4. **자발적 호전은 절반만 재현된다.** 섬유화 항을 끄면 담관 재생만으로
   빌리루빈이 1→18세에 3.50→2.20으로 떨어진다. 전체 모델에서는 섬유화가
   담관 회복을 앞질러 3.50→4.48로 오른다. 실제 아동이 어느 쪽을 따르는지가
   임상적 핵심 질문이며 이 모델은 결정하지 못한다 — 재생 속도를 적합할
   연속 생검 DPR 시계열이 존재하지 않기 때문이다. 혈청 담즙산은 어느
   쪽이든 거의 움직이지 않으며(236.9→233.4), 이 자체가 예측이다: 부분적
   담관 회복은 **담즙산보다 빌리루빈을 훨씬 먼저 개선**시킨다.
5. **위험도는 코호트로 보정하고 개인에게 적용한다.** GALA 적합에 필요한
   frailty 분산이 2.94로 매우 커서, 코호트 곡선은 상당 부분 환자 간 분산에
   대한 진술이다. 시나리오 1은 보정 점검이지 예후가 아니다.
6. **간경변 약리학과 위약군 진행이 없다.** FIB 3 이상에서도 동일한
   담관/합성 구조를 적용하지만 비대상성 ALGS 간의 간세포 부전 생리는 담고
   있지 않다. 또 ASSERT 위약군 담즙산은 24주에 22 µmol/L **상승**했는데
   이 모델의 무치료군은 평탄(−0.9)하므로, 그 상승을 만드는 무언가를 담고
   있지 않고 절대 변화량이 그만큼 과대평가된다(−109 vs −90).

### 실행해서 발견한 결함 9건 (결론을 바꿨을 것들)

1. `$MAIN`이 매 실행마다 초기조건을 재설정해 외부 `init()`을 조용히
   덮어썼고, 그 결과 **모든 담관 용량 시나리오가 동일한 기저 빌리루빈
   3.5를 보고**해 Axis 2의 임계값 자체가 측정 불가능했다.
2. 시나리오 1이 개인 수준 `exp(-H)` 생존을 코호트 곡선과 비교해 18년 자가
   간 생존을 실제 40.3%에 대해 **1.5%**로 보고했다.
3. Φ_EHC가 처음에 간세포→혈장→간세포 무익회로를 포함해, 모든 환자에서
   ~0.98로 고정되며 그 양의 의미를 파괴했다.
4. `EPS`와 `THETA`는 mrgsolve 예약어라 컴파일에 실패했다.
5. 주석 속 아포스트로피가 작은따옴표 R 모델 문자열을 종료시켰다.
6. `$TABLE`이 출력 시점에 재계산하지 않고 마지막 `$ODE` 도함수 호출이 남긴
   지역변수를 캡처했다.
7. 첫 섬유화 속도는 400일 만에 FIB를 2로 몰았고, 수정본은 17년간 1.0에
   얼려버렸다.
8. 볼루스 대 연속 투여 검증이 정확히 투여 주기인 정수 일에 샘플링해
   **39%의 순수 에일리어싱 불일치**를 보고했다 (세밀 격자에서 5.7%).
9. 콜레스티라민 군이 미셀 형성에 대해서만 담즙산을 결합하고 ASBT 흡수에
   대해서는 결합하지 않아, 담즙산 효과가 전혀 없는 순수 흡수장애 약처럼
   보였다.

---

## 면책

연구 및 교육 목적의 모델입니다. 임상적으로 검증되지 않았으며 진료 의사결정
도구가 아닙니다. 어떤 파라미터도 환자 관리에 사용해서는 안 됩니다.
