# 트랜스티레틴 아밀로이드증 (ATTR) QSP 모델

> **디렉토리:** `transthyretin-amyloidosis/` | **약어:** ATTR | **날짜:** 2026-06-24  
> **분류:** 희귀질환 / 단백질 접힘 이상 / 아밀로이드증

[![ATTR QSP 기계론적 지도](attr_qsp_model.png)](attr_qsp_model.svg)

---

## 질환 개요

**트랜스티레틴 아밀로이드증(Transthyretin Amyloidosis, ATTR)** 은 간에서 합성되는 수송 단백질 **트랜스티레틴(TTR)** 이 정상 사량체(tetramer) 구조에서 해리되어 잘못 접힌 단량체 → 독성 올리고머 → 성숙 아밀로이드 섬유로 집합하는 단백질 접힘 이상 질환입니다.

| 특성 | ATTRwt (야생형) | ATTRv (유전성) |
|------|----------------|----------------|
| 원인 | 노화 관련 TTR 불안정화 | *TTR* 유전자 미스센스 변이 |
| 발병 연령 | ≥60세, 남성 우세 | 30–70세, 변이마다 다름 |
| 주요 표현형 | 심근병증(HFpEF→HFrEF) | 다발신경병증 / 심장 혼합 |
| 대표 변이 | — | V30M(포르투갈/일본/스웨덴), V122I(아프리카계) |
| 유병률 | ~10–13% 고령 HF 부검 | ~50,000명/전 세계 |

---

## 핵심 발병기전 — 4단계

| 단계 | 기전 | 핵심 매개자 |
|------|------|------------|
| **1. TTR 합성** | 간 헤파토사이트에서 TTR 단량체 합성 → β-시트 접힘 → 동형사량체 조립 → 혈장 분비 | TTR mRNA, pre-TTR → 성숙 TTR(14 kDa), 사량체(55 kDa) |
| **2. 사량체 해리** | 열·pH 스트레스 또는 ATTRv 변이에 의한 사량체 해리 → 잘못 접힌 단량체 생성 → 핵 형성(nucleation) → 올리고머 | 해리속도(kdis), 올리고머(독성 중간체), 프로토피브릴 |
| **3. 아밀로이드 침착** | 성숙 섬유가 심장(ATTRwt≫), 말초신경(ATTRv≫), 신장·비장·GI에 선택적 침착 | FIB_HRT, FIB_NRV, FIB_SYS |
| **4. 장기 손상** | 침착 → 세포 독성(NLRP3·IL-1β·TNF-α) → 심근세포 아포토시스 → 심실 비후·이완 장애 / 슈반세포 압박 → 축삭 퇴화 → 다발신경병증 | LVEF↓, NT-proBNP↑, NIS↑, mBMI↓ |

---

## 치료 기전 개요

| 약물 | 기전 | 표적 | 주요 임상시험 |
|------|------|------|-------------|
| **타파미디스** (61mg PO QD) | T4 결합부위 점유 → 사량체 운동적 안정화 (Emax ~80%) | TTR 사량체 해리 ↓ | ATTR-ACT (Maurer 2018 NEJM) |
| **파티시란** (0.3mg/kg IV Q3W) | LNP-siRNA → ApoE-LDLR 간 흡수 → RISC/Ago2 → TTR mRNA 절단 | TTR mRNA ↓80% | APOLLO (Adams 2018 NEJM) |
| **부트리시란** (25mg SC Q3M) | GalNAc-siRNA → ASGR1 → RISC → mRNA 분해 | TTR mRNA ↓83% | HELIOS-A (Gillmore 2021 NEJM) |
| **이노테르센** (300mg SC QW) | 2'-MOE ASO → RNase H1 → TTR mRNA 절단 | TTR mRNA ↓72% | NEURO-TTR (Benson 2018 Lancet) |
| **아코라미디스** (800mg PO BID) | 고선택적 T4부위 결합 → 사량체 안정화 | TTR 사량체 해리 ↓ | ATTRiBUTE-CM (Elliott 2023 NEJM) |
| **디플루니살** (250mg PO BID) | NSAID + T4부위 약한 결합 | TTR 안정화(비선택적) | Berk 2013 JAMA |

---

## 모델 파일 구성

| 파일 | 사양 | 설명 |
|------|------|------|
| [`attr_qsp_model.dot`](attr_qsp_model.dot) | **116 노드, 10 클러스터** | Graphviz 기계론적 지도 (DOT 소스) |
| [`attr_qsp_model.svg`](attr_qsp_model.svg) | 벡터 그래픽 | 확대 가능 SVG |
| [`attr_qsp_model.png`](attr_qsp_model.png) | 150 dpi PNG | 미리보기용 래스터 이미지 |
| [`attr_mrgsolve_model.R`](attr_mrgsolve_model.R) | **25구획 ODE, 7치료 시나리오** | mrgsolve PK/PD 모델 |
| [`attr_shiny_app.R`](attr_shiny_app.R) | **8탭 대시보드** | Shiny 인터랙티브 앱 |
| [`attr_references.md`](attr_references.md) | **60개 PubMed 인용 (10섹션)** | 근거 문헌 목록 |

---

## mrgsolve 모델 상세 — 25개 ODE 구획

### 약물 PK 구획 (9개)

| 구획 | 약물 | 경로 | 핵심 파라미터 |
|------|------|------|-------------|
| A_TAF_GUT, A_TAF_C, A_TAF_P | 타파미디스 | 경구 2구획 | ka=0.42/h, CL=0.96L/h, t½~55h |
| A_VUT_SC, A_VUT_C | 부트리시란 | 피하 1구획 | ka=0.08/h, F=82%, t½~5d |
| A_INO_SC, A_INO_C | 이노테르센 | 피하 1구획 | ka=0.10/h, F=70%, CL=0.04L/h, t½~30d |
| A_PAT_C, A_PAT_P | 파티시란 | 정맥 2구획 | CL=0.18L/h, Vc=3.3L, t½~3d |

### 질환 PD 구획 (16개)

| 구획 | 생물학적 의미 | 주요 동역학 |
|------|-------------|-----------|
| TTR_MRNA | 간 TTR mRNA (정규화) | kin(1-ERNA) - kout·mRNA; E_RNA = E_VUT + E_PAT + E_INO |
| TTR_TET | 혈장 TTR 사량체 | ksyn·mRNA - kdis(1-Estab)·TET - kout·TET |
| TTR_MONO | 잘못 접힌 단량체 | 2·kdis·TET - kagg·MONO - kdeg·MONO |
| TTR_OLIGO | 독성 올리고머 | kagg·MONO - kfib·OLIGO - kdeg·OLIGO |
| FIB_HRT | 심장 아밀로이드 부하 | kfib·OLIGO·f_heart - kdeg_FIB·FIB_HRT |
| FIB_NRV | 말초신경 아밀로이드 | kfib·OLIGO·f_nerve - kdeg_FIB·FIB_NRV |
| FIB_SYS | 전신 아밀로이드 | kfib·OLIGO·f_sys - kdeg_FIB·FIB_SYS |
| INFLAM | 심장 염증 지수 | kin·FIB_HRT/(FIB50+FIB_HRT) - kout·INFLAM |
| LVEF | 좌심실 구혈률 (%) | krec·(EF_base-EF) - kdet·FIB·INFLAM·EF |
| NT_proBNP | NT-proBNP (pg/mL) | kin·(INFLAM + 1/EF) - kout·NT_proBNP |
| NIS | 신경병증 손상 점수 | kin_NIS·FIB_NRV - kout_NIS·NIS |
| mBMI | 수정 체질량지수 | -kdet·NIS/(NIS+NIS50)·mBMI |
| eGFR | 사구체여과율 | -kdet_eGFR·FIB_SYS·eGFR |
| SYMP_CARD | 심장 증상 누적 | NYHA 유사 지표 통합 |
| SYMP_NEURO | 신경 증상 누적 | FAP stage 유사 지표 통합 |

---

## 7가지 치료 시나리오

| # | 시나리오 | 약물·용량 | 대상 | 임상 근거 |
|---|----------|----------|------|----------|
| S1 | 자연경과 ATTRwt | 없음 | ATTRwt-CM | Maurer 2018 (위약군) |
| S2 | 자연경과 ATTRv | 없음 | ATTRv-NP | Adams 2018 (위약군) |
| S3 | 타파미디스 | 61mg PO QD | ATTRwt-CM | ATTR-ACT: CV사망+입원 30% ↓ |
| S4 | 파티시란 | 0.3mg/kg IV Q3W | ATTRv-NP | APOLLO: mNIS+7 34점 차이 |
| S5 | 부트리시란 | 25mg SC Q3M | ATTRv-NP | HELIOS-A: NIS 17점 개선 |
| S6 | 이노테르센 | 300mg SC QW | ATTRv-NP | NEURO-TTR: mNIS+7 19점 차이 |
| S7 | 타파미디스+부트리시란 | 병용 | 가상 | 이중 기전 (탐색용) |

---

## Shiny 앱 탭 구성 (8개)

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | ATTR 표현형·치료 설정·바이오마커 빠른 요약 |
| 2. 약물 PK | 4종 약물 혈중 농도·약물 효과(Estab, ERNA) 시각화 |
| 3. TTR 접힘 이상 | mRNA 감소·사량체·올리고머·섬유 축적 동역학 |
| 4. 심장 결과 | LVEF·NT-proBNP·심장 섬유 부하·NYHA 추정 |
| 5. 신경 결과 | NIS·mBMI·신경 섬유 부하·eGFR |
| 6. 시나리오 비교 | 5개 치료 옵션 직접 비교·18개월 요약 표 |
| 7. 바이오마커 대시보드 | 심장/신경 바이오마커 통합 패널·시점별 요약 |
| 8. 모델 정보 | ODE 구조·파라미터 근거·핵심 참고문헌 |

---

## 실행 방법

```r
# 1) mrgsolve ODE 모델 실행
install.packages(c("mrgsolve","dplyr","ggplot2","tidyr","patchwork"))
source("transthyretin-amyloidosis/attr_mrgsolve_model.R")

# 2) Shiny 대시보드 실행
install.packages(c("shiny","shinydashboard","plotly","DT"))
shiny::runApp("transthyretin-amyloidosis/attr_shiny_app.R")

# 3) Graphviz 지도 렌더링
dot -Tsvg attr_qsp_model.dot -o attr_qsp_model.svg
dot -Tpng -Gdpi=150 attr_qsp_model.dot -o attr_qsp_model.png
```

---

## 주요 임상시험 결과 요약

| 임상시험 | 약물 | N | 1차 결과 | 효과크기 |
|---------|------|---|---------|---------|
| **ATTR-ACT** (2018) | 타파미디스 61mg | 441 | CV사망+HF입원 | HR 0.70 (95%CI 0.51-0.96) |
| **APOLLO** (2018) | 파티시란 0.3mg/kg | 225 | mNIS+7 개선 | -34점 vs 위약 (p<0.001) |
| **HELIOS-A** (2021) | 부트리시란 25mg SC | 164 | mNIS+7 개선 | -17점 vs 위약 적응적비교 |
| **NEURO-TTR** (2018) | 이노테르센 300mg | 172 | mNIS+7 개선 | -19점 vs 위약 (p<0.001) |
| **ATTRiBUTE-CM** (2023) | 아코라미디스 800mg | 632 | 6MWT+KCCQ | p=0.0066 (계층적 복합) |

---

*QSP Disease Model Library — 매일 Claude Code Routine이 새 모델을 추가합니다.*
