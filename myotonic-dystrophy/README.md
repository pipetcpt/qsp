# 근긴장성 이영양증 제1형 (Myotonic Dystrophy Type 1) — QSP Model

> **디렉토리:** `myotonic-dystrophy/` | **약어:** DM1 | **날짜:** 2026-06-26  
> **분류:** 신경근육 (Neuromuscular) | **OMIM:** #160900

[![DM1 QSP Mechanistic Map](dm1_qsp_model.png)](dm1_qsp_model.svg)

---

## 병태생리 핵심 경로

```
DMPK Gene (Chr 19q13.32)
  └─► CTG Repeat Expansion (>50 → up to 10,000+ repeats)
        └─► Mutant DMPK mRNA (CUGn) retained in nucleus
              └─► CUGn RNA Hairpin → Nuclear RNA Foci
                    ├─► MBNL1 Sequestration (loss-of-function)
                    └─► PKCβII Activation → CUGBP1/CELF1 Hyperphosphorylation (gain-of-function)
                          │
                          ▼ MBNL1↓ + CUGBP1↑ = Fetal Splicing Program Reversion
                          │
                          ├─► CLCN1 fetal isoform → ClC-1 ↓ → Myotonia
                          ├─► INSR fetal isoform  → Insulin resistance
                          ├─► SERCA1 fetal isoform → Ca²⁺ dysregulation
                          ├─► TNNT2 fetal isoform → Cardiac contractile defect
                          ├─► MAPT 4R/3R imbalance → Tau CNS pathology
                          └─► CAMK2D/KCNQ1 splicing → Arrhythmia risk
```

---

## QSP 모델 구성

### 기계론적 지도 (Mechanistic Map)
- **노드 수:** ~158개 (10 서브그래프 클러스터)
- **클러스터:** ①유전·분자 ②RNA결합단백질 ③선택적스플라이싱 ④골격근 ⑤심장 ⑥CNS·인지 ⑦내분비·대사 ⑧약동학 ⑨약력학 ⑩임상엔드포인트

### mrgsolve ODE 모델 (22 구획)

| 구획 | 변수명 | 설명 |
|------|--------|------|
| PK-1 | MEX_GUT | 멕실레틴 장관 흡수 |
| PK-2 | MEX_CENT | 멕실레틴 중앙 혈장 |
| PK-3 | MEX_PERI | 멕실레틴 주변 조직 |
| PK-4 | ASO_PLASMA | ASO 혈장 농도 |
| PK-5 | ASO_MUSCLE | ASO 근육 조직 |
| PK-6 | ASO_NUCL | ASO 핵내 농도 (활성) |
| BIO-1 | CUG_FOCI | CUG RNA foci 부담 (0-1) |
| BIO-2 | MBNL1_FREE | 유리 MBNL1 분율 |
| BIO-3 | CUGBP1_ACT | CUGBP1 활성화 수준 |
| BIO-4 | CLCN1_FETAL | CLCN1 태아형 스플라이싱 분율 |
| BIO-5 | SERCA_FETAL | SERCA1 태아형 스플라이싱 분율 |
| BIO-6 | INSR_FETAL | INSR 태아형 스플라이싱 분율 |
| END-1 | MYOTONIA | 근긴장 VAS 점수 (0-10) |
| END-2 | GRIP_STR | 악력 (kg) |
| END-3 | MUSCLE_MASS | 골격근 질량 (kg) |
| END-4 | PR_INT | PR 간격 (ms) |
| END-5 | QTc_INT | QTc 간격 (ms) |
| END-6 | HOMA_IR | HOMA-IR (인슐린 저항성) |
| END-7 | FVC_PCT | FVC % 예측치 |

### 치료 시나리오 (7가지)

| # | 시나리오 | 투여 | 주요 표적 |
|---|----------|------|-----------|
| 1 | 자연 경과 | 없음 | — |
| 2 | 멕실레틴 200 mg TID | 경구 8시간마다 | Nav1.4 차단 |
| 3 | 멕실레틴 300 mg TID (MELT 용량) | 경구 8시간마다 | Nav1.4 차단 |
| 4 | ASO 4주 간격 (DYNE-101 regimen) | 피하/정맥 28일마다 | CUG RNA 분해 |
| 5 | 멕실레틴 300 mg + ASO 복합 | 복합 | Nav1.4 + CUG RNA |
| 6 | 유전자 치료 (AAV-MBNL1, 실험) | 근육내 단회 | MBNL1 복원 |
| 7 | 중증 DM1 (CTG=1200, 무치료) | 없음 | 중증 자연 경과 |

### Shiny 대시보드 (7탭)

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | CTG 반복 크기 · 중증도 분류 · 초기 상태 |
| 2. 약물 PK | 멕실레틴 혈장 농도 · ASO 조직 농도 |
| 3. 근육·근긴장 | VAS 점수 · 악력 · Nav1.4/ClC-1 활성 |
| 4. 심장 안전성 | QTc · PR 간격 · 심장 위험 분류 |
| 5. 시나리오 비교 | 다중 치료 비교 · 1년 결과 표 |
| 6. 바이오마커 패널 | CLCN1/INSR/SERCA1 스플라이싱 지수 |
| 7. CNS·대사 | HOMA-IR · FVC% · 전신 합병증 |

---

## 주요 임상시험 보정 (Clinical Trial Calibration)

| 임상시험 | 치료 | 결과 | 참고문헌 |
|---------|------|------|---------|
| Logigian 2010 RCT | 멕실레틴 150/200 mg TID | VAS 근긴장 −2.6점 vs 위약 | Neurology 2010 |
| MELT 2018 | 멕실레틴 300 mg TID | 안전성 확인, 근긴장 개선 | Muscle Nerve 2018 |
| Ionis-DMPK-2.5Rx 2015 | ISIS 598769 SC | DMPK mRNA 50–80% 감소 | Cunningham 2015 |
| DYNE-101 Ph2 | AOC 1001 IV | CLCN1 스플라이싱 지수 ~40pp 개선 | 2023 |
| Groh 2008 NEJM | 관찰 연구 | HV > 70 ms → SCD 위험 5배 | NEJM 2008 |

---

## 파일 목록

| 파일 | 설명 |
|------|------|
| [`dm1_qsp_model.dot`](dm1_qsp_model.dot) | Graphviz 기계론적 지도 (~158 노드, 10 클러스터) |
| [`dm1_qsp_model.svg`](dm1_qsp_model.svg) | 벡터 형식 지도 (고해상도) |
| [`dm1_qsp_model.png`](dm1_qsp_model.png) | 래스터 형식 지도 (150 dpi) |
| [`dm1_mrgsolve_model.R`](dm1_mrgsolve_model.R) | mrgsolve ODE QSP 모델 (22 구획, 7 시나리오) |
| [`dm1_shiny_app.R`](dm1_shiny_app.R) | Shiny 인터랙티브 대시보드 (7탭) |
| [`dm1_references.md`](dm1_references.md) | 참고문헌 40편 (10개 섹션) |

---

## 실행 방법

```bash
# 기계론적 지도 렌더링
dot -Tsvg myotonic-dystrophy/dm1_qsp_model.dot -o myotonic-dystrophy/dm1_qsp_model.svg
dot -Tpng -Gdpi=150 myotonic-dystrophy/dm1_qsp_model.dot -o myotonic-dystrophy/dm1_qsp_model.png
```

```r
# mrgsolve 모델
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
source("myotonic-dystrophy/dm1_mrgsolve_model.R")

# Shiny 대시보드
install.packages(c("shiny", "shinydashboard", "plotly", "DT", "scales"))
shiny::runApp("myotonic-dystrophy/dm1_shiny_app.R")
```

---

## 핵심 약어

| 약어 | 설명 |
|------|------|
| DMPK | Dystrophia Myotonica Protein Kinase |
| CTG | Cytosine-Thymine-Guanine (trinucleotide repeat) |
| MBNL1 | Muscleblind-Like Protein 1 |
| CUGBP1/CELF1 | CUG-Binding Protein / CELF Family |
| ClC-1 | Chloride Channel 1 (skeletal muscle) |
| Nav1.4 | Voltage-gated Sodium Channel, skeletal muscle |
| INSR | Insulin Receptor |
| SERCA1 | Sarco/Endoplasmic Reticulum Ca²⁺-ATPase 1 |
| ASO | Antisense Oligonucleotide |
| MELT | Mexiletine Evaluation of Toxicity trial |
| EDS | Excessive Daytime Sleepiness |
| HOMA-IR | Homeostatic Model Assessment for Insulin Resistance |
| VAS | Visual Analogue Scale |
| MIRS | Muscular Impairment Rating Scale |
