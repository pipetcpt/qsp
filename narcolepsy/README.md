# Narcolepsy Type 1 — QSP Model

> **질환 분류**: 신경계 / 수면장애 (Neurological / Sleep Disorder)  
> **디렉토리**: `narcolepsy/`  
> **약어**: `narc`

---

## 병태생리 개요 (Pathophysiology)

나르콜렙시 1형(NT1)은 **시상하부 외측의 오렉신(히포크레틴) 신경세포가 자가면역 기전으로 선택적으로 파괴**되어 발생하는 만성 신경계 질환입니다. 정상인에서 약 70,000개 존재하는 오렉신 신경세포가 NT1 환자에서는 **85–95% 소실**되며, 뇌척수액 히포크레틴-1(OXA) 농도가 110 pg/mL 미만으로 떨어집니다.

### 핵심 병리 경로

| 경로 | 세부 기전 |
|------|----------|
| **자가면역 병인** | HLA-DQB1\*06:02 감수성 유전자 + H1N1 인플루엔자/Pandemrix 백신 분자 모방 → CD4+/CD8+ T세포·NK세포가 오렉신 신경세포 파괴 |
| **오렉신 수용체 시스템** | OXA(OX1R≫OX2R), OXB(OX2R≥OX1R) → Gq/Gi 이중 신호 → 각성 아민계 활성화 |
| **각성 시스템 붕괴** | LC(NE)·TMN(히스타민)·VTA(DA)·DRN(5-HT)·기저전뇌(ACh) 활성 불안정 |
| **플립-플롭 스위치 불안정** | 오렉신 결핍 → VLPO–각성계 상호억제 균형 파괴 → 각성/수면 상태 전환 빈번 |
| **카탈렉시** | 편도체 감정 자극 → REM 근긴장소실 회로(SubC 글루타메이트계) 불수의 활성화 |
| **대사 합병증** | 오렉신 결핍 → 렙틴 저항성·에너지 균형 이상 → 비만·T2DM 위험 증가 |

---

## 모델 파일 구성 (Model Files)

| 파일 | 설명 |
|------|------|
| [`narc_qsp_model.dot`](narc_qsp_model.dot) | Graphviz 기계론적 지도 (소스) |
| [`narc_qsp_model.svg`](narc_qsp_model.svg) | 벡터 그래픽 (확대 가능) |
| [`narc_qsp_model.png`](narc_qsp_model.png) | 래스터 이미지 (150 dpi) |
| [`narc_mrgsolve_model.R`](narc_mrgsolve_model.R) | mrgsolve ODE QSP 모델 |
| [`narc_shiny_app.R`](narc_shiny_app.R) | Shiny 인터랙티브 대시보드 |
| [`narc_references.md`](narc_references.md) | 참고문헌 (35+ PubMed 인용) |

---

## 기계론적 지도 (Mechanistic Map)

[![Narcolepsy QSP Map](narc_qsp_model.png)](narc_qsp_model.svg)

### 클러스터 구성 (12개 서브그래프)

1. **자가면역 병인** — HLA-DQB1\*06:02, T세포, 분자 모방, 자가항체
2. **시상하부 오렉신 신경세포** — HCRT 유전자, OXA/OXB, 신경세포 소실
3. **오렉신 수용체 시스템** — OX1R, OX2R, Gq/Gi 신호전달
4. **각성 모노아민 시스템** — LC, TMN, VTA, DRN, 기저전뇌
5. **수면 촉진 시스템** — VLPO, 아데노신, 일주기 리듬(SCN)
6. **REM 수면 조절 및 카탈렉시** — PPT/LDT, SubC, 편도체
7. **임상 증상 4징** — EDS, 카탈렉시, 입면 환각, 수면 마비
8. **약물 PK 구획** — 소디움 옥시베이트, 모다피닐, 피톨리산트, 솔리암페톨
9. **약리 표적 및 효과** — GABA-B, DAT/NET, H3R, SERT, OX2R
10. **바이오마커 및 엔드포인트** — CSF OXA, MSLT, PSG, ESS
11. **플립-플롭 스위치** — 각성/수면 이중안정 상태, NT1 불안정성
12. **대사 합병증** — 비만, T2DM, 우울/불안, ADHD 유사 증상

---

## mrgsolve ODE 모델 (Pharmacological Model)

### 구획 (23개)

| 구획 | 변수명 | 설명 |
|------|--------|------|
| **소디움 옥시베이트 PK** | GUT_OXY, CENT_OXY, PERI_OXY | 3구획 모델, t½ ≈ 30–60분 |
| **모다피닐 PK** | GUT_MOD, CENT_MOD | 2구획, t½ ≈ 15시간 |
| **피톨리산트 PK** | GUT_PIT, CENT_PIT | 2구획, t½ ≈ 10–12시간 |
| **솔리암페톨 PK** | GUT_SOL, CENT_SOL | 2구획, t½ ≈ 7.1시간 |
| **벤라팍신 PK** | GUT_VEN, CENT_VEN | 2구획, 항카탈렉시 |
| **각성 신경계** | WAKE_LC, WAKE_TMN, WAKE_VTA, WAKE_DRN | 모노아민 발화율 (0–1) |
| **수면 시스템** | SLEEP_P, VLPO_ACT, ADENOSINE | 항상성 수면압·VLPO·아데노신 |
| **상태 변수** | WAKE_STATE, REM_STATE, NREM_STATE | 수면/각성/REM 플립-플롭 |
| **임상 축적** | EDS_ACC, CATAPLEXY_ACC | 누적 졸림·카탈렉시 |

### 치료 시나리오 (7개)

| 시나리오 | 약물 | 용량 |
|----------|------|------|
| 1 | 무치료 NT1 기저 | — |
| 2 | 소디움 옥시베이트 | 4.5 g + 2.75 g 분할 (취침 시) |
| 3 | 모다피닐 | 200 mg 1일 1회 |
| 4 | 피톨리산트 | 18 mg 1일 1회 |
| 5 | 솔리암페톨 | 150 mg 1일 1회 |
| 6 | 소디움 옥시베이트 + 피톨리산트 | 병합 |
| 7 | 벤라팍신 | 75 mg (카탈렉시 중점) |

### 주요 임상시험 보정

| 참조 | 치료 | 결과 |
|------|------|------|
| Black 2010 (Sleep Med) | 소디움 옥시베이트 | 카탈렉시 -69–75% |
| HARMONY I (Szakacs 2017, Lancet Neurol) | 피톨리산트 | ESS -5.8 vs 위약 -3.4 |
| TONES 3 (Schweitzer 2019, Sleep) | 솔리암페톨 | ESS -7.7점 |
| US MMSG 2000 | 모다피닐 | ESS -4.3점 |
| Mignot 2002 (Lancet) | CSF 오렉신-1 | 진단 기준 <110 pg/mL |

---

## Shiny 대시보드 (Interactive Dashboard)

```r
shiny::runApp("narcolepsy/narc_shiny_app.R")
```

### 탭 구성 (8개)

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 오렉신 신경세포 생존율, ESS, 카탈렉시 빈도, HLA 상태 |
| 2. Drug PK | 약물 농도-시간 곡선, Cmax/Tmax/AUC 테이블 |
| 3. 수면-각성 조절 | 모노아민 시스템 활성, 플립-플롭 스위치, 24시간 일주기 |
| 4. 임상 엔드포인트 | ESS·카탈렉시·MSLT·수면 구조 추적 |
| 5. 치료 시나리오 비교 | 임상시험 벤치마크 폭포수 차트·산림 그림 |
| 6. 바이오마커 | CSF 히포크레틴-1, MSLT, HLA 진단 시각화 |
| 7. 자가면역 기전 | HLA 연관성, 신경세포 파괴 타임라인, T세포 |
| 8. 참고문헌 | 주요 임상시험 요약 테이블 |

---

## 실행 방법 (Usage)

```bash
# 기계론적 지도 렌더링 (Graphviz 필요)
dot -Tsvg narcolepsy/narc_qsp_model.dot -o narc_qsp_model.svg
dot -Tpng -Gdpi=150 narcolepsy/narc_qsp_model.dot -o narc_qsp_model.png
```

```r
# mrgsolve 모델 실행 (R ≥ 4.0)
install.packages(c("mrgsolve", "dplyr", "ggplot2"))
source("narcolepsy/narc_mrgsolve_model.R")

# Shiny 대시보드 실행
install.packages(c("shiny", "shinydashboard", "plotly"))
shiny::runApp("narcolepsy/narc_shiny_app.R")
```

---

## 핵심 파라미터 (Key Parameters)

| 파라미터 | 값 | 출처 |
|---------|-----|------|
| 정상 CSF 오렉신-1 | 200–300 pg/mL | Mignot 2002 |
| NT1 진단 기준 | <110 pg/mL | ICSD-3 |
| 오렉신 신경세포 소실 | 85–95% | Thannickal 2000 |
| HLA-DQB1\*06:02 유병률 | ~95% (NT1) vs 25% (일반) | Mignot 1997 |
| 유병률 | 25–50/100,000 | Longstreth 2007 |
| 평균 진단 지연 | 10–15년 | Thorpy 2015 |
| Emax (피톨리산트 → ESS) | -5.8점 | HARMONY I |
| Emax (솔리암페톨 → ESS) | -7.7점 | TONES 3 |

---

## 면책 조항

본 모델은 **교육 및 연구 목적**으로 제작된 정성적·반정량적 QSP 모델입니다. 공개 문헌 기반으로 구성되었으며 임상 의사결정에 직접 사용할 수 없습니다.
