# 악성 흉막 중피종 QSP 모델 (Malignant Pleural Mesothelioma)

[![Disease](https://img.shields.io/badge/Disease-Malignant%20Pleural%20Mesothelioma-blue)]()
[![Category](https://img.shields.io/badge/Category-Thoracic%20Oncology-orange)]()
[![Compartments](https://img.shields.io/badge/ODE%20Compartments-53-green)]()
[![Scenarios](https://img.shields.io/badge/Treatment%20Scenarios-20-purple)]()
[![References](https://img.shields.io/badge/References-213-yellow)]()
[![Map](https://img.shields.io/badge/Map-252%20nodes%20%C2%B7%2019%20clusters-lightgrey)]()

## 개요 (Overview)

**악성 흉막 중피종(MPM)** 은 석면 노출 후 30–50년의 잠복기를 거쳐 발생하는
흉막의 원발성 악성종양입니다. 그러나 이 모델이 붙잡는 것은 예후의 나쁨이 아니라
**종양의 모양**입니다.

**Malignant pleural mesothelioma is not a mass. It is a rind** — a sheet of
tumour spread over the parietal and visceral pleural surfaces whose only
geometric degree of freedom is its **thickness**. Almost every quantitative
intuition imported from other solid tumours — what a 30% response means, how
deep a drug gets, what shrinking the tumour does for the patient — breaks on
that one fact. This model is built to make the breakage computable.

---

## 이 모델의 하나의 구조적 약속 (The one structural commitment)

**종양 부담은 상태변수가 아닙니다.** 상태변수는 세 흉막면(벽측·장측·엽간열)
위의 **질량**이고, 두께는 거기서 계산됩니다:

```
h_i = ( T_i + M_i + N_i ) / ( rho x A_i )
        생존세포  기질   괴사물          면적
```

그리고 임상적으로 의미 있는 모든 양이 **같은 껍질의 서로 다른 범함수**입니다:

```
부피        V = A h                          -> h에 선형
mRECIST     = 4 h_par + 2 h_vis              -> h에 선형
세포에 닿는 약  f = (lam/h)(1 - exp(-h/lam))  -> h에 열등선형 (포화)
폐 포획      E_add = k x phi x h_vis          -> 장측 잎에만
흉수        벽측 껍질이 림프 stomata를 덮는 정도의 함수
```

구(sphere)에서는 `dV/dr = 4(pi)r²`, 껍질에서는 **`dV/dh = A`로 상수**입니다.

| 측정된 30% 선형 감소 | 함축되는 부피 살상 |
|---|---|
| 구형 종양 (RECIST가 설계된 대상) | **65.7 %** |
| 껍질 (중피종, mRECIST) | **30.0 %** |
| 비율 | **2.19 배** |

같은 "부분 관해(partial response)"라는 단어가 두 질환에서 **2.2배 다른 세포
살상**을 뜻합니다. 이것은 가정이 아니라 산수이고, 모델은 이 변환을 출력합니다
(벽측 잎 1 mm = 45 mL, 세 잎 합계 1 mm = 96 mL).

---

## 하나의 콜라겐 상태, 세 가지 결과 (One collagen state, three consequences)

껍질의 콜라겐 분율 `phi = M / (T + M + N)` 은 모델에서 **정확히 세 곳**에만
들어갑니다.

1. **두께 안에** — 세포가 죽어도 기질은 사라지지 않는다 → 측정 편향
2. **침투 길이 안에** — `lambda(phi) = lambda_0 (1-phi)^1.5` → 전달 저하
3. **장측 잎의 부가 탄성 안에** — `E_add = k phi h_vis` → 폐 포획

보정된 모델에 페메트렉시드+시스플라틴 6주기를 투여하면 셋이 동시에 나타납니다:

| | 기저 | 반응 최저점 (117일) |
|---|---|---|
| 생존 종양세포 | — | **−60.3 %** |
| mRECIST 합 | 32.0 mm | **−30.1 %** |
| 콜라겐 분율 phi | 0.265 | **0.499** |
| 소분자 깊이평균 노출 f | 0.257 | **0.215** |
| 부가 탄성 E_add | 8.15 cmH₂O/L | **10.01** |
| FVC | 2.27 L | **2.11 L** |

세 가지가 한꺼번에 읽힙니다.

- **측정이 살상의 절반을 숨긴다.** 생존세포는 60% 줄었는데 두께는 30%만
  줄었습니다. 차이 27.7 %p는 콜라겐과 아직 청소되지 않은 괴사물입니다. 이 간극은
  요법에 거의 무관하게 일정합니다 — pem/cis 27.7 · +bev 28.2 ·
  nivo+ipi 28.4 · gem/cis 24.4 %p.
- **반응하는 과정에서 약이 덜 들어간다.** 껍질이 얇아졌는데도 침투 분율이
  0.257 → 0.215로 **떨어집니다.** 얇아진 이득보다 섬유화된 손해가 큽니다.
- **반응해도 숨이 편해지지 않는다.** 30% 영상 반응을 거치는 동안 부가 탄성은
  8.15 → 10.01로 **올라가고** FVC는 2.27 → 2.11 L로 오히려 떨어집니다.
  큰 영상 반응이 작은 생존 이득으로만 이어지는 이 질환의 특징이,
  가정이 아니라 결과로 나옵니다.

부수적으로 모델은 `d(부피)%` 와 `d(mRECIST)%` 가 소수점 첫째 자리까지
일치함을 보여 (−32.1 vs −32.1) 껍질 기하학의 선형성을 스스로 확인합니다.

---

## 두 방향의 전달, 둘 다 수 밀리미터 (Two routes, opposite gradients)

전신 투여는 껍질의 **혈관이 있는 바닥**에서, 흉강내 투여는 **자유 표면**에서
들어옵니다. 방향은 반대지만 식은 같습니다.

| phi = 0.25 | 침투 길이 | 6 mm 껍질에서의 f |
|---|---|---|
| 소분자 (전신) | 1.62 mm | 0.257 |
| IgG 항체 (전신) | 0.39 mm | **0.065** |
| T세포 침윤 | 0.31 mm | 0.052 |
| 흉강내 (자유 표면) | 2.3 mm | 0.66 |

6 mm 껍질이 보는 항체 간질 농도는 혈장 평형값의 **6.5%**, 1 mm 껍질이 보는 값의
18%입니다. `f`는 `lambda/h`로 위가 막혀 있으므로 **용량 증량으로 고칠 수 없는
전달 문제**입니다. 흥미롭게도 이것이 제한하는 것은 수용체 점유가 아닙니다 —
니볼루맙은 혈장 농도만으로 PD-1을 포화시킵니다. 제한되는 것은 **T세포 침윤**
입니다.

흉강내 경로는 잔존 두께 12 mm에서 1 mm 잔존 대비 노출의 **23%** 밖에 받지
못합니다. 모든 양성 흉강내 치료 결과가 **육안적 완전 절제 이후**에만 보고되고
절제 대신으로는 보고되지 않는 이유가 여기서 나옵니다.

---

## 보정 (Calibration): 5개 상수, 5개 종점

각 상수는 **일차원 이분법으로 단 하나의 종점**에만 맞춰졌고, 앞 단계를
흐트러뜨리지 않는 순서로 배열됩니다 (`mpm_calibration.py`).

```
HZ0      = 0.00045625    <- 최선지지요법 중앙 생존 7.0개월
EMAX_CIS = 0.148438      <- 시스플라틴 단독 9.3개월 (EMPHACIS 대조군)
EMAX_PEM = 0.0283594     <- 페메트렉시드+시스플라틴 12.1개월 (EMPHACIS)
BEV_LAM  = 0.295625      <- 베바시주맙 증분 +2.7개월 (MAPS)
EMAX_IO  = 0.624424      <- 니볼루맙+이필리무맙/화학요법 중앙비 1.28 (CM-743)
```

| 요법 | 관측 OS | 모델 OS | 역할 |
|---|---|---|---|
| 최선지지요법 (역사적) | 7.0 mo | 7.0 mo | 보정 |
| 시스플라틴 단독 | 9.3 mo | 9.3 mo | 보정 |
| 페메트렉시드 + 시스플라틴 | 12.1 mo | 12.1 mo | 보정 |
| + 베바시주맙 (MAPS 증분) | 14.8 mo | 14.8 mo | 보정 |
| 니볼루맙 + 이필리무맙 (CM-743 비) | 15.5 mo | 15.5 mo | 보정 |
| 젬시타빈 + 시스플라틴 | 11.2 mo | 9.7 mo | **보류** |
| 페메트렉시드 + 카보플라틴 | 12.7 mo | 14.9 mo | **보류** |
| 니볼루맙 단독 2차 (CONFIRM) | 10.2 mo | 7.7 mo | **보류** |

8개 군 중앙 생존의 상대 RMSE **11.6%**. 절대 생존기간은 EMPHACIS에 고정되어
있습니다 — MAPS와 CheckMate 743은 더 좋은 환자군을 등록했으므로 증분과 비로만
맞췄습니다.

---

## 예측이 틀린 곳 (Where the model is wrong)

이 저장소의 다른 모델들과 달리, 이 모델의 **사전 등록된 예측은 실패했습니다.**
숨기지 않고 그대로 보고합니다.

### 1. CheckMate 743 조직형 하위군 — 방향이 반대다

보정에 쓴 것은 CM-743 전체 집단의 비 **한 점**뿐이고, 상피양/비상피양 분리는
예측으로 남겨두었습니다. 그 예측이 **부호가 아니라 방향에서** 틀렸습니다.

| | 관측 (HR) | 모델 (중앙비) |
|---|---|---|
| 상피양 | 0.86 | 0.70 |
| 비상피양 | **0.46** | **0.87** |

관측은 육종양에서 면역관문 억제의 이점이 **커지는데**, 모델은 **작아집니다.**

`mpm_emt_sensitivity.py`는 이것이 기울기 하나를 잘못 놓아서가 아님을 보입니다.
화학요법 저항을 `EMT_CHEMO = 0.97`(육종양에서 화학요법이 살상력의 3%만 유지)까지
밀어도 비상피양 비는 0.79에 그쳐 상피양 0.68과 여전히 **순서가 반대**이고, PD-L1을 분율이
허용하는 한계까지 올려도 두 하위군이 같이 움직일 뿐 갈라지지 않습니다.

원인은 구조적입니다. 모델은 조직형 벌점을 **양쪽 군 모두에** 매깁니다 —
`HZ_EMT`는 기저 위험을, `EMT_KCOL`은 콜라겐을 통해 백금과 **똑같이 T세포의**
침투 길이를, `EMT_KG`는 증식을 올립니다. 그래서 x=0.85에서 두 군이 4.5개월과
5.7개월로 함께 무너지고, 작은 두 수의 비는 1에서 멀어질 수 없습니다.
CheckMate 743은 정반대를 말합니다 — 비상피양 환자가 면역요법에서 **18.1개월**로
상피양(18.7)과 사실상 같았고, 무너진 것은 화학요법(8.8개월)뿐이었습니다.
따라서 필요한 것은 **면역 효과기가 그 벌점을 면제받는 기전**이며, 그것이
무엇인지는 이 코드가 아니라 중피종에 대한 질문입니다.

### 2. MARS2 — 수술의 부호가 반대다

MARS2(2024)는 확대 흉막절제/박피술 + 화학요법이 화학요법 단독보다 **나빴다**고
보고했습니다 (OS HR 1.28). 모델은 반대로 큰 이득을 예측합니다 (26.5 vs
12.1개월). 이 모델에서 수술은 **두께만 깨끗이 되돌리는 조작**이고, 시험이
말하는 해악의 원천이 전부 빠져 있기 때문입니다 — 수술 사망(MARS2 30일 3.4%),
개흉술이 남기는 수개월의 수행능력·폐용적 부채, 그리고 결정적으로 **"육안적 완전
절제"가 남기는 미세 잔존 병변은 45% 피복률의 얇은 껍질이 아니라 흉막면 거의
전체에 걸친 세포**라는 사실입니다. 두께가 0.1 mm로 갈 때 피복률이 1로 가야
하는데 이 모델에서는 그렇지 않습니다.

### 3. 그 외 알려진 불일치

- **집단 관해율이 높다.** 가상 집단(n=24)에서 pem/cis ORR 71% (관측 41%),
  중앙 OS 16.0개월 (관측 12.1). 개별 환자의 중앙 생존값을 모아 그 중앙값을
  취하는 방식이 오른쪽으로 치우칩니다. 결정론적 대표 환자는 종점에 정확히
  맞으므로, 집단 층은 예시용으로만 읽어야 합니다.
- **면역요법의 PFS가 과도하게 길다** (35.3개월 vs 관측 6.8). 대표 환자의
  반응이 너무 깊어 최저점에서 20% 반등에 도달하지 않습니다.
- **화학요법 2년 생존이 낮다** (5% vs 관측 27%). 면역요법 쪽은 39% vs 41%로
  잘 맞습니다.

---

## 그 외 검증에서 나온 것 (Other results)

**호흡곤란의 두 원인은 분리 가능하다.** 흉수는 제거 가능하고 폐 포획은 그렇지
않습니다. 흉막 압력계 문헌의 포획 역치(탄성 14.5 cmH₂O/L)는 이 모델에서
`phi = 0.29`일 때 **장측 껍질 7.0 mm**에서 교차합니다. 흉막유착술은 폐가 펴져야
성립하므로 모델에서는 유착 형성 속도에 apposition 항으로 들어갑니다 —
90일째 탈크군은 유착 0.86 · 흉수 81 mL · FVC 2.04 L, 유치도관군은 흉수 337 mL ·
FVC 1.96 L입니다.

**두 약이 신장을 통해 서로를 증폭한다.** 시스플라틴 신독성으로 eGFR이
90 → 74 mL/min으로 떨어지면 페메트렉시드 청소율이 기저의 82%로 떨어져 노출이
올라가고 후반 주기의 골수 최저점이 깊어집니다. 1주기 ANC 최저 1.46 ×10⁹/L
(11일째, 관측 8–12일), 엽산·B12 보충을 빼면 **0.52로 64% 더 깊어집니다** —
EMPHACIS가 시험 도중 보충을 추가하고 3–4등급 독성이 떨어진 것과 같은 방향입니다.

**SMRP는 두께가 아니라 살아있는 세포를 따라간다.** 90일째 생존세포 −52.6%,
mRECIST −22.6%, SMRP −45.2%. 다만 SMRP는 신장으로 배설되므로 위의 eGFR 하락이
값을 끌어올리는 교란이 있고, 모델은 그 교란도 함께 보여줍니다.

**수치 위생.** dt = 0.08일은 골수 통과 사슬(ktr = 0.87/d)의 안정 영역 밖이며
그 1년 종점은 무의미합니다. dt = 0.04 이하에서는 수렴하여 스텝을 반으로 줄여도
1년 종점이 유효숫자 4자리에서만 움직입니다. 365일 동안 음수로 가는 상태변수는
없습니다.

---

## 파일 (Files)

| 파일 | 내용 |
|---|---|
| [`mpm_qsp_model.dot`](mpm_qsp_model.dot) · [`.svg`](mpm_qsp_model.svg) · [`.png`](mpm_qsp_model.png) | 기계론적 지도 — 252 노드 · 335 엣지 · 19 클러스터 |
| [`mpm_mrgsolve_model.R`](mpm_mrgsolve_model.R) | mrgsolve 모델 — 53 ODE 구획, 20 시나리오, 종점 추출·교차곡선 함수 |
| [`mpm_shiny_app.R`](mpm_shiny_app.R) | Shiny 대시보드 — 12개 탭 |
| [`mpm_reference_model.py`](mpm_reference_model.py) | 의존성 없는 Python RK4 참조 구현 + 12개 구조 검사 |
| [`mpm_calibration.py`](mpm_calibration.py) | 단계별 이분법 보정 + 가상 집단 |
| [`mpm_calibration_output.txt`](mpm_calibration_output.txt) | 위 스크립트의 실제 실행 출력 (README 수치의 출처) |
| [`mpm_emt_sensitivity.py`](mpm_emt_sensitivity.py) · [`_output.txt`](mpm_emt_sensitivity_output.txt) | 실패한 하위군 예측에 대한 사후 민감도 분석 |
| [`mpm_fetch_references.py`](mpm_fetch_references.py) | PubMed E-utilities 조회 스크립트 |
| [`mpm_references.md`](mpm_references.md) | 문헌 213편 (모두 실제 조회 결과) |

### 실행 (Running)

```bash
dot -Tsvg mpm_qsp_model.dot -o mpm_qsp_model.svg
dot -Tpng -Gdpi=150 mpm_qsp_model.dot -o mpm_qsp_model.png

python3 mpm_reference_model.py        # 참조 구현 + 12개 구조 검사 (R 불필요)
python3 mpm_calibration.py 24         # 보정 + 가상 집단
python3 mpm_emt_sensitivity.py        # 실패한 예측의 사후 분석
python3 mpm_fetch_references.py       # 문헌 재수집
```

```r
source("mpm_mrgsolve_model.R")
endpoint_table(); crossover_curve()
shiny::runApp("mpm_shiny_app.R")
```

---

## 모델 구조 (Model structure)

### 53개 ODE 구획

| 그룹 | 구획 |
|---|---|
| 세포독성 PK (8) | `CIS1` `CIS2` `CIST` `PEM1` `PEM2` `PEMT` `GEM1` `GEMT` |
| 생물학적 제제 PK (7) | `BEV1` `BEV2` `NIV1` `NIV2` `NIVT` `IPI1` `IPI2` |
| 흉강내 투여 (2) | `IPPL` `IPT` |
| 껍질 (9) | `TPAR` `TVIS` `TFIS` · `MPAR` `MVIS` `MFIS` · `NPAR` `NVIS` `NFIS` |
| 면역 (6) | `TEFF` `TREG` `PRIME` `IFNG` `PDL1` `TAM` |
| 신호전달 (3) | `VEGF` `TGFB` `IL6` |
| 흉막강 (2) | `PLV` `SYMPH` |
| 호흡역학 (2) | `VEXP` `DYSP` |
| 독성 (9) | `PROL` `TR1` `TR2` `TR3` `ANC` · `PTK` `GFR` `IRAE` `NEURO` |
| 숙주·바이오마커·생존 (5) | `CACHEX` `SMRP` `ARG` `CH` + 예약 슬롯 `ECOGs` |

### 20개 치료 시나리오

`bsc` · `cis` · `pemcis` · `pemcarbo` · `gemcis` · `pemcisbev` ·
`nivoipi` · `nivo` · `pembrochemo` · `chemo_then_io` ·
`adi_pemcis` (ASS1 결핍) · `ttf_chemo` ·
`pd_chemo` (MARS2) · `epp_chemo` (MARS) · `surg_hithoc` · `ip_only` ·
`talc` · `ipc` · `talc_chemo` · `pemcis_nofolate`

### 지도의 19개 클러스터

1 석면 노출·섬유 생체잔류성 · 2 만성 염증(HMGB1·NLRP3·IL-1β) ·
3 유전체(BAP1·NF2·CDKN2A) · 4 Hippo-YAP과 세포주기 · 5 조직형 EMT 축 ·
6 섬유화 기질 · **7 껍질 기하학** · **8 두 방향의 약물 전달** ·
9 백금/항엽산 PK-PD · 10 면역 미세환경 · 11 면역관문 억제 PK-PD ·
12 항혈관·표적 치료 · 13 흉막강(흉수·유착술) · 14 호흡역학과 증상 ·
15 수술·방사선 · 16 독성 · 17 바이오마커 · 18 임상 종점 ·
19 세 가지 구조적 주장

---

## 한계 (Limitations)

- **표면 피복률이 고정 파라미터입니다.** 진단 시 피복률(벽측 0.45 · 장측 0.35 ·
  엽간열 0.50)을 상수로 두므로, 종양이 흉막면을 따라 퍼져나가는 단계와
  절제 후 미세 잔존이 전 표면에 남는 상황을 기술하지 못합니다. 위의 MARS2
  불일치가 정확히 여기서 나옵니다.
- **깊이 평균 근사.** `f = (lam/h)(1-exp(-h/lam))`는 껍질을 한 구획으로 두고
  노출의 깊이 평균을 취한 평균장 근사이며, 바깥층이 먼저 죽고 안쪽이 남는
  공간 구조를 명시적으로 갖지 않습니다.
- **골수·신장 노출을 종양 구획의 잔류 화학종으로 대리합니다.** 혈장 반감기가
  2시간 미만이어서 혈장 농도로는 독성의 시간척도가 나오지 않기 때문이며,
  골수가 실제로 보는 폴리글루타메이트·백금 부가체와 같은 종이라는 근사입니다.
- **생존은 단일 위험함수**이고 중앙 생존은 누적 위험이 ln2에 도달하는 시점으로
  정의됩니다. 개별 환자 생존시간은 `mpm_calibration.py`의 가상 집단을 쓰십시오.
- 복막·심막·고환초막 중피종은 다루지 않습니다. 기하학이 다릅니다.

---

## 면책 (Disclaimer)

본 모델은 **교육 및 연구 목적의 정성적·반정량적 QSP 모델**입니다. 공개 문헌과
임상시험 요약값을 바탕으로 구성되었으나 독립적으로 검증·인증되지 않았으며,
**실제 임상 의사결정, 처방, 또는 규제 제출에 직접 사용해서는 안 됩니다.**
