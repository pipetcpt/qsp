# 남성 성선기능저하증 QSP 모델 (Male Hypogonadism)

**Male Hypogonadism · MHG** — 시상하부-뇌하수체-생식샘(HPG) 축, Vermeulen 결합
평형, 고환내 테스토스테론, 적혈구 생성, 골, 체성분을 하나의 기계론적 모델로 묶은
정량적 시스템 약리학(QSP) 모델.

> ⚠️ **교육 및 연구 목적 전용.** 임상 의사결정, 처방, 규제 제출에 사용하지 마십시오.

---

## 이 모델이 하려는 말 (The claim)

**측정하는 값과 작용하는 값이 다르고, 그 사이의 지도는 비선형이다.**

임상에서 재는 것은 총 테스토스테론(TT)이지만 실제로 수용체에 닿는 것은 유리
테스토스테론(FT)이다. 둘 사이는 SHBG가 지배하는 질량작용 평형이고, SHBG 결합은
**포화한다**. 그래서 FT는 TT의 **볼록함수**다 (SHBG 35 nmol/L, 알부민 4.3 g/dL):

| 총 T (ng/dL) | 100 | 300 | 600 | 1000 | 1500 | 2000 |
|---|---|---|---|---|---|---|
| 유리 T (pg/mL) | 17.5 | 56.0 | 122.9 | 228.2 | 381.1 | 550.8 |
| 유리 분율 (%) | 1.75 | 1.87 | 2.05 | 2.28 | 2.54 | **2.75** |

유리 분율이 총 T와 함께 **올라간다**. 이 한 줄에서 아래 세 가지가 따라 나오고,
모델의 분석 함수들은 각각을 산문이 아니라 **숫자로** 출력한다.

### ① 진단 역치는 프레임 의존적이다

단일 총 T 역치(300 ng/dL)는 단일 유리 T 역치(65 pg/mL)가 **아니다**.
`MHG_diagnostic_frame()` 출력:

| SHBG (nmol/L) | TT=300에서의 FT (pg/mL) | FT=65가 되는 TT (ng/dL) | 판정 |
|---|---|---|---|
| 15 | 83.5 | **237** | TT 역치가 과잉진단 |
| 25 | 67.3 | 290 | 거의 일치 |
| 35 | 56.0 | 344 | TT 역치가 과소진단 |
| 55 | 41.5 | 450 | TT 역치가 과소진단 |
| 90 | 28.2 | **635** | TT 역치가 과소진단 |

같은 유리 T 역치가 총 T 기준으로는 **237 → 635 ng/dL, 2.7배** 움직인다.
어느 쪽이 성선기능저하인지를 결정하는 것이 생식샘이 아니라 **결합 단백질**인
구간이 존재한다는 뜻이다. 방향도 양쪽이다 — 비만/인슐린 저항성 남성(낮은 SHBG)은
과잉진단되고, 고령 남성(높은 SHBG)은 과소진단된다. `MHG_shbg_paradox()`가 두
환자를 나란히 시뮬레이션한다.

### ② 파형은 예상보다 훨씬 덜 중요하다 — 모델이 자기 가설을 반증한다

Jensen 부등식에 따라 E[FT(TT)] ≥ FT(E[TT]) 이므로, 평균 총 T가 같아도 첨두-저점
진폭이 큰 요법이 시간평균 유리 T를 더 많이 전달한다. 이 모델은 원래 그 볼록성이
**볼록한 적혈구 반응과 합성되어** 근주 대 경피 적혈구증가증 격차(약 3배)를
설명하리라는 가설로 만들어졌다. **모델을 돌려보니 그렇지 않았다.**

`MHG_convexity_decomposition()`이 실제로 출력하는 산술:

| 요소 | 모델 출력 | 해석 |
|---|---|---|
| (a) 볼록 **결합** | 겔 +0.3% → 근주 q2주 **+4.4%** | 실재하며 진폭에 비례해 커진다. 수학이 요구하는 그대로 |
| (b) 볼록 **적혈구 반응**과의 합성 | **≈ 0 또는 음수** | 첨두가 EC50를 넘어 Hill의 **오목한** 상단으로 진입. 아래에서 볼록한 반응은 위에서 오목하다 |
| (c) 용량을 맞춘 **파형 효과** | Hct **~0.1점** | 주간 용량 동일·진폭 2배(100 mg 주1회 vs 200 mg 격주)에서의 실제 차이 |
| (d) 고정 역치(Hct>54%)를 **집단 분포**에 적용 | 발생률 **~2.5배** | 평균 ~1점 이동이 발생률 수 배 차이로 증폭 |

즉 문헌의 근주-경피 격차는 **대부분 용량 차이이고, 나머지는 판정 규칙이 만든다.**
(b)는 이 모델이 표현하려던 가설을 스스로 반증한 지점이며, 튜닝으로 지우지 않고
그대로 출력한다.

> **노출된 가정**: `EC50_HEPC = 300 pg/mL`(정상 유리 T 50-210보다 높게)는
> **보정값이며 측정값이 아니다.** `MHG_hepcidin_sensitivity()`가 이 값을 150에서
> 600까지 움직였을 때 위 결론이 얼마나 살아남는지 표로 출력한다.

### ③ 혈청 T와 고환내 T는 서로 다른 변수이고, 하나만 측정된다

고환내 테스토스테론(ITT)은 혈청의 약 35-100배(약 700 nmol/L)이며, 정자형성은
ITT가 정상의 약 20-30% 아래로 떨어지면 정지한다. 외인성 테스토스테론은 혈청 T를
정상화하면서 LH를 끄고 ITT를 몇 퍼센트로 붕괴시킨다.

`MHG_ITT_collapse()`가 두 값을 같은 표에 넣어 해리를 보여준다. **혈청 T가 정상인데
ITT가 한 자릿수인 행들이 곧 "호르몬을 대체하면서 환자를 불임으로 만드는" 요법이다.**
LHCGR 신호를 남겨두는 요법(hCG 병용, 클로미펜)만이 ITT를 역치 위로 유지한다.
대체(replacement)와 회복(restoration)은 다른 개입이다.

### ④ 보너스: "테스토스테론의 효과" 중 상당수는 에스트라디올의 효과다

`MHG_finkelstein()`이 NEJM 2013 설계(GnRH 작용제 배경 + 단계적 T 겔 ± 아나스트로졸)를
재현한다. 제지방량과 근력은 T를 따라가지만 **체지방은 E2를 따라간다.** 그래서
여성유방증을 막으려고 아로마타제를 억제하는 것은 공짜가 아니다 —
`MHG_aromatase_cost()`가 그 대가를 요추 해면골 vBMD로 환산해서 청구한다.

---

## 파일 구성 (Files)

| 파일 | 내용 |
|------|------|
| [`mhg_qsp_model.dot`](mhg_qsp_model.dot) | 기계론적 지도 소스 — **224 노드 · 268 엣지 · 22 클러스터** |
| [`mhg_qsp_model.svg`](mhg_qsp_model.svg) | 벡터 지도 (확대 가능) |
| [`mhg_qsp_model.png`](mhg_qsp_model.png) | 래스터 지도 (150 dpi) |
| [`mhg_mrgsolve_model.R`](mhg_mrgsolve_model.R) | **49 ODE 구획** · 8개 T 제형 + 5개 비안드로겐 약물 · 14 시나리오 · 11 분석 함수 |
| [`mhg_shiny_app.R`](mhg_shiny_app.R) | 10개 탭 인터랙티브 대시보드 |
| [`mhg_references.md`](mhg_references.md) | **78개 참고문헌** (PubMed 링크) + "모델이 재현하지 못하는 것" 부록 |

---

## 기계론적 지도 (Mechanistic Map)

[![MHG QSP map](mhg_qsp_model.png)](mhg_qsp_model.svg)

22개 클러스터:

| # | 클러스터 | # | 클러스터 |
|---|---|---|---|
| ① | 시상하부 — GnRH 박동 생성기 (KNDy) | ⑫ | 적혈구 생성 (헵시딘·EPO·철) |
| ② | 뇌하수체 전엽 — 성선자극세포 | ⑬ | 골 (E2 우세 경로) |
| ③ | Leydig 세포 스테로이드생성 | ⑭ | 근육·체성분 |
| ④ | Sertoli 세포·정자형성 | ⑮ | 중추신경 — 성욕·발기·기분 |
| ⑤ | **순환 결합 평형 (수학적 중심)** | ⑯ | 전립선·LUTS (포화 모델) |
| ⑥ | SHBG 조절 | ⑰ | 심혈관·안전성 |
| ⑦ | 말초 대사 (5α-환원·방향화) | ⑱ | 테스토스테론 제형 PK |
| ⑧ | AR / ER 수용체 신호전달 | ⑲ | 비안드로겐 전략 |
| ⑨ | 원발성(고환성) 원인 | ⑳ | 외인성 안드로겐에 의한 축 차단 |
| ⑩ | 이차성(중추성) 원인 | ㉑ | 진단·모니터링 |
| ⑪ | 기능성 — 비만·염증 악순환 | ㉒ | 임상 종점·주요 시험 |

지도에서 분홍색 노트(⚑)는 모델의 정량적 핵심 주장이며, 각각 R 모델의 분석 함수
하나에 대응한다.

---

## mrgsolve 모델 (49 ODE · 132 파라미터)

```r
source("mhg_mrgsolve_model.R")
MHG_run_all()                      # 11개 분석 함수 전부 실행

d <- MHG_scenario_im_q2wk()        # 근주 시피온산 200 mg q2wk
MHG_plot_overview(d)
MHG_plot_waveforms()               # 제형별 파형 비교
```

### 구획 구성

| 모듈 | 구획 |
|------|------|
| 테스토스테론 PK | `DEP_IM` `DEP_TU` `DEP_SC` `DEP_GEL` `DEP_ORAL` `DEP_PEL` `DEP_NAS` `CENT` `PERIPH` |
| 비안드로겐 약물 PK | `HCG_D/C` `FSHD_D/C` `CLO_D/C` `ANA_D/C` |
| HPG 축 | `GNRHD` `LH` `FSH` `LEYCAP` `ITT` `INHB` |
| 호르몬·결합 | `SHBG` `E2` `DHT` |
| 정자형성 (74+14일 전달 사슬) | `SG1` `SG2` `SG3` `SG4` `EPID` |
| 적혈구 생성 | `HEPC` `EPO` `PROG_E` `RETIC` `RBC` |
| 골 | `SCLERO` `OB` `OC` `BMD_TR` `BMD_CO` `RSP`(리모델링 공간) |
| 체성분·기타 | `LEAN` `FAT` `PSA` `LIBIDO` `VITAL` |
| 노출 적분 | `CUMFT` `CUMDRIVE` |

### 환자 원형 (patient archetypes) — 무치료 기저값

모델을 실제로 돌려 얻은 값이며, 각 원인별로 임상에서 기대되는 검사 패턴을 재현한다.

| 원형 | TT (ng/dL) | FT (pg/mL) | SHBG | LH | FSH | ITT (%) | 패턴 |
|---|---|---|---|---|---|---|---|
| `PT_ORGANIC` 원발성 고환부전 | 232 | 38.8 | 40 | **14.8** | 13.6 | 29 | 고성선자극 |
| `PT_FUNCTIONAL` 비만 기능성 | 300 | 58.9 | **32** | 1.5 | 1.9 | 47 | 저성선자극 · SHBG 낮음 |
| `PT_ELDERLY` 노화 관련 | 272 | 40.4 | **49** | 6.1 | 6.2 | 31 | 혼합형 · SHBG 높음 |
| `PT_SECONDARY` 뇌하수체성 | 276 | 51.1 | 35 | 1.3 | 1.7 | 40 | 저성선자극 |
| `PT_KLINEFELTER` 47,XXY | 172 | 30.7 | 35 | 15.7 | **31.4** | 22 | FSH ≫ LH |
| `PT_OPIOID` 아편유사제 | 261 | 44.6 | 40 | 1.1 | 1.5 | 34 | 저성선자극 · 가역적 |

### 치료 시나리오 (14)

| # | 시나리오 | 함수 |
|---|---|---|
| 1 | 무치료 기능성 성선저하 자연경과 (3년) | `MHG_scenario_natural()` |
| 2 | 근주 시피온산 200 mg q2wk | `MHG_scenario_im_q2wk()` |
| 3 | 근주 시피온산 100 mg 주 1회 (동일 주간 용량) | `MHG_scenario_im_weekly()` |
| 4 | 피하 자동주사기 75 mg 주 1회 | `MHG_scenario_sc()` |
| 5 | 경피 겔 1.62% 81 mg/day | `MHG_scenario_gel()` |
| 6 | 경구 운데칸산 237 mg BID | `MHG_scenario_oral()` |
| 7 | 근주 운데칸산 1000 mg q12wk | `MHG_scenario_tu_im()` |
| 8 | 피하 펠릿 750 mg q4개월 | `MHG_scenario_pellet()` |
| 9 | hCG 단독 1500 IU 주 3회 (생식능 보존) | `MHG_scenario_hcg_mono()` |
| 10 | 겔 + hCG 500 IU EOD (ITT 구제) | `MHG_scenario_t_plus_hcg()` |
| 11 | 클로미펜 25 mg 매일 | `MHG_scenario_clomiphene()` |
| 12 | 체중 감량 단독 (약물 없음) | `MHG_scenario_weight_loss()` |
| 13 | 아편유사제 유발 → 180일에 중단 | `MHG_scenario_opioid()` |
| 14 | 3년 TRT 후 중단 — 축·정자 회복 | `MHG_scenario_cessation()` |

### 분석 함수 (11)

| 함수 | 출력 |
|------|------|
| `MHG_free_T_nomogram()` | TT × SHBG 격자의 유리 T, 유리 분율, 볼록성 검증 |
| `MHG_diagnostic_frame()` | 총 T 역치가 SHBG에 따라 2.7배 움직이는 표 |
| `MHG_shbg_paradox()` | 비만 vs 고령 두 환자 — TT와 FT 판정이 엇갈림 |
| `MHG_convexity_decomposition()` | **핵심 원장**: (a) 결합 · (b) 반응 · (c) 용량맞춤 파형 · (d) 역치 4요소 분해 |
| `MHG_hepcidin_sensitivity()` | `EC50_HEPC`에 대한 결론의 민감도 |
| `MHG_ITT_collapse()` | 혈청 T와 ITT를 같은 표에 — 대체 vs 회복 |
| `MHG_recovery_curve()` | 중단 후 정자 회복 (Liu 2006 대조) |
| `MHG_finkelstein()` | T 대 E2 기여 해리 (NEJM 2013 설계 재현) |
| `MHG_aromatase_cost()` | 아로마타제 억제의 대가를 vBMD로 환산 |
| `MHG_formulation_ledger()` | 7개 제형 × 12개 종점 1년 대조표 |
| `MHG_trial_ledger()` | 모델 vs 발표 종점 — **재현 실패 항목 포함** |

---

## Shiny 대시보드 (10 탭)

```r
setwd("male-hypogonadism"); shiny::runApp("mhg_shiny_app.R")
```

① 환자 프로파일 · ② 결합 평형(SHBG 노모그램) · ③ 진단 역치의 프레임 의존성 ·
④ PK 파형 · ⑤ 적혈구증가증 볼록성 원장 · ⑥ 고환내 T·생식능 ·
⑦ 골·체성분 · ⑧ 바이오마커 · ⑨ 시나리오 비교 · ⑩ 임상시험 대조표

사이드바에서 6개 환자 유형, 12개 요법(최대 3개 동시 비교), 나이·SHBG·인슐린
저항성·지방량, 그리고 **`EC50_HEPC` 가정 자체**를 직접 움직여 볼 수 있다.

---

## 모델이 재현하지 **못하는** 것

정직한 QSP 모델은 자신이 틀리는 지점을 명시해야 한다. 가장 중요한 항목:

> **TRAVERSE 골절 하위연구 (NEJM 2024): 임상 골절 3.50% vs 2.46%, HR 1.43.**
> 테스토스테론을 준 쪽에서 **골절이 더 많았다.** 모델은 T-Trials의 요추 해면골
> vBMD 상승(발표 +7.5%, 모델 +5.7~+7.7%)을 재현하도록 보정되어 있으므로, 구조상
> 골절 증가를 예측할 수 **없다.** 이것은 파라미터 조정으로 고칠 문제가 아니라
> **BMD를 골절의 대리지표로 삼는 가정 자체의 실패**이며, `MHG_trial_ledger()`가
> 이를 숨기지 않고 출력한다.
>
> 덧붙여, 모델이 그 BMD 상승을 만들어내는 기전 자체가 이 경고를 강화한다. 상승분의
> 상당 부분은 새 뼈가 아니라 **리모델링 공간(remodelling space)의 회복**이다 —
> 성선기능저하 상태의 높은 골교체가 만들어 놓은, 아직 메워지지 않은 흡수와(cavity)가
> 안드로겐 회복으로 골교체가 느려지면서 채워지는 것이다(항흡수제 초기 BMD 상승과
> 같은 기전). 밀도계가 기록하는 이득이 곧 강도의 이득이라는 보장은 없다.

그 밖에: TRAVERSE MACE(결정론적 모델에 사건 과정 없음), LH 박동성(평균 농도로
단순화), SHBG 알로스테릭 결합(Vermeulen 단순 모델 채택), 개체간 변이(AR CAG 반복
등 미구현). 전체 목록은 [`mhg_references.md` 부록](mhg_references.md#부록-모델이-재현하지-못하는-것-what-this-model-does-not-reproduce) 참조.

---

## 보정 표적 요약 (Calibration Targets)

| 종점 | 출처 | 모델에서의 역할 |
|------|------|-----------------|
| 생산율 6 mg/day, MCR ~1000 L/day | Southren 1965 | `KSPILL × ITT0` 항등식의 앵커 |
| K_SHBG 1×10⁹, K_Alb 3.6×10⁴ M⁻¹ | Vermeulen 1999 | 결합 평형 그 자체 |
| ITT ≈ 혈청의 35-100배 | Jarow 2001 · Roth 2010 | `ITT0 = 700 nmol/L` |
| T 200 mg 주간 투여 → ITT −94% | Coviello 2005 | `HCG_POT` · `ITT50_S` |
| 정자형성 주기 74일 | Heller & Clermont 1964 | `TAU_SPG` |
| 회복: 6개월 67% · 12개월 90% | Liu 2006 (Lancet) | `MHG_recovery_curve()` |
| 요추 해면골 vBMD +7.5% (1년) | T-Trials Bone (JAMA IM 2017) | `KFORM`/`KRES` |
| 체지방 증가는 E2 결핍 탓 | Finkelstein (NEJM 2013) | `SFAT_E` |
| 근주 vs 경피 적혈구증가증 ~3배 | Ohlander 2018 | `EC50_HEPC` (**보정값**) |
| MACE HR 0.96 | TRAVERSE (NEJM 2023) | 재현 시도 안 함 (명시) |
| 골절 HR 1.43 | TRAVERSE 골절 (NEJM 2024) | **재현 실패 (명시)** |

---

## 재현 (Reproduce)

```bash
# 지도 렌더링
dot -Tsvg mhg_qsp_model.dot -o mhg_qsp_model.svg
dot -Tpng -Gdpi=150 mhg_qsp_model.dot -o mhg_qsp_model.png

# 모델 실행
Rscript -e 'source("mhg_mrgsolve_model.R"); MHG_run_all()'
```

필요 패키지: `mrgsolve`, `shiny`. 지도 렌더링에는 Graphviz.

---

> **면책 조항** — 본 모델은 공개 문헌과 임상시험 데이터를 바탕으로 구성된
> 교육·연구용 QSP 모델입니다. 독립적으로 검증·인증되지 않았으며 실제 임상
> 의사결정에 사용해서는 안 됩니다.
