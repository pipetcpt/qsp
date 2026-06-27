# 신세포암 (ccRCC) QSP 모델 — Renal Cell Carcinoma

[![Nodes](https://img.shields.io/badge/Nodes-125%2B-blue)](rcc_qsp_model.dot)
[![ODE Compartments](https://img.shields.io/badge/ODE%20Compartments-18-green)](rcc_mrgsolve_model.R)
[![Regimens](https://img.shields.io/badge/Regimens-7-orange)](rcc_mrgsolve_model.R)
[![References](https://img.shields.io/badge/References-60-red)](rcc_references.md)

## 개요 (Overview)

투명세포 신세포암(ccRCC)은 신장암의 약 75%를 차지하며, **VHL 종양억제유전자 돌연변이**(60–90%)로 인한 pVHL 소실 → HIF-1α/2α 축적 → VEGF·CA9·GLUT1 과발현이 핵심 발병기전입니다. 이 모델은 VHL/HIF/VEGF 산소 감지 경로, PI3K/AKT/mTOR, MAPK/RAS, 종양 면역 미세환경(TME), 약물 PK/PD를 통합하는 정량적 시스템 약리학(QSP) 프레임워크입니다.

## 발병기전 요약 (Pathophysiology)

| 경로 | 핵심 이벤트 |
|------|------------|
| VHL/HIF 산소 감지 | VHL 돌연변이 → pVHL 소실 → HIF-2α 유비퀴틴화 실패 → HRE 전사 활성화 |
| 혈관신생 | HIF-2α → VEGF-A↑ → VEGFR2 인산화 → PI3K/MAPK → 내피세포 증식·이동 |
| mTOR 피드백 | PI3K/AKT → mTORC1 활성화 → HIF 번역 증가 (양성 피드백) |
| 면역 회피 | VEGF → MDSC 모집, Treg 유도, PD-L1 발현 → CD8+ T세포 기능 억제 |
| 대사 재프로그래밍 | HIF-1α → GLUT1·LDHA↑ → Warburg 효과 |

## 모델 파일 (Model Files)

| 파일 | 내용 |
|------|------|
| [`rcc_qsp_model.dot`](rcc_qsp_model.dot) | Graphviz 기계론적 지도 (125+ 노드, 11 클러스터) |
| [`rcc_qsp_model.svg`](rcc_qsp_model.svg) | SVG 벡터 지도 |
| [`rcc_qsp_model.png`](rcc_qsp_model.png) | PNG 썸네일 (150 dpi) |
| [`rcc_mrgsolve_model.R`](rcc_mrgsolve_model.R) | mrgsolve ODE 모델 (18 구획, 7 치료 시나리오) |
| [`rcc_shiny_app.R`](rcc_shiny_app.R) | 7탭 Shiny 인터랙티브 대시보드 |
| [`rcc_references.md`](rcc_references.md) | 60개 PubMed 참고문헌 (14 섹션) |

## 기계론적 지도 (Mechanistic Map)

[![ccRCC QSP Model](rcc_qsp_model.png)](rcc_qsp_model.svg)

*클릭하면 확대 SVG 보기*

## ODE 구획 (ODE Compartments — 18)

| 그룹 | 구획 | 설명 |
|------|------|------|
| Sunitinib PK | DEPOT_SUN, CENT_SUN, PERI_SUN, MET_SUN | 2-구획 PK + SU12662 활성 대사체 |
| Nivolumab TMDD | CENT_NIV, PERI_NIV, PD1_FREE, PD1_BOUND | 표적 매개 약동학 (PD-1 결합) |
| Belzutifan | DEPOT_BEZ, CENT_BEZ | HIF-2α 억제제 PK |
| VHL/HIF/VEGF | pVHL, HIF2A, VEGF, VEGFR2_ACT | 산소 감지 경로 |
| mTOR | mTOR_ACT | mTORC1 활성도 (AU) |
| Tumor (Simeoni) | TUM_W1, TUM_W2, TUM_W3, TUM_VOL | 3-전달 TGI 모델 |
| Immune TME | CD8_T, TREG, MDSC | 종양 면역 미세환경 |

## 치료 시나리오 (Treatment Scenarios)

| # | 요법 | 임상 근거 | 핵심 파라미터 |
|---|------|----------|--------------|
| 1 | Untreated | 대조군 | — |
| 2 | Sunitinib 50 mg (4/2 주기) | NEJM 2007 (Motzer) | CL=51.8 L/h, SU12662 활성 대사체 |
| 3 | Nivolumab + Ipilimumab | CheckMate 214 (Motzer 2018) | TMDD: kon=0.32 nM⁻¹h⁻¹ |
| 4 | Pembrolizumab + Axitinib | KEYNOTE-426 (Rini 2019) | 이중 VEGFR+PD-1 차단 |
| 5 | Cabozantinib + Nivolumab | CheckMate 9ER (Choueiri 2021) | MET/VEGFR2 + PD-1 |
| 6 | Cabozantinib 60 mg | METEOR (Choueiri 2015) | IC50=0.006 µM |
| 7 | Everolimus 10 mg | RECORD-1 (Motzer 2008) | mTOR IC50=0.15 nM |
| 8 | Belzutifan 120 mg | LITESPARK-005 (Choueiri 2023) | HIF-2α IC50=0.018 µM |

## Shiny 앱 탭 구성 (Shiny App Tabs)

| 탭 | 내용 |
|----|------|
| 1. Patient Profile | IMDC 위험도, 종양 초기 부피, 시뮬레이션 기간 설정 |
| 2. PK | Sunitinib/SU12662/Nivolumab/Belzutifan 농도-시간 곡선 |
| 3. VHL/HIF/VEGF Pathway | HIF-2α·VEGF·VEGFR2·mTOR 동역학, HIF-2α SS vs pVHL 곡선 |
| 4. Tumor Dynamics | Simeoni TGI 종양 부피 곡선, 폭포 차트(BOR), TGI 통계 |
| 5. Immune TME | CD8·Treg·MDSC 동역학, PD-1 점유율, MDSC vs VEGF 산점도 |
| 6. Scenario Comparison | 8개 요법 동시 비교 (종양 부피 + 엔드포인트 테이블) |
| 7. Biomarker Dashboard | 바이오마커 Z-score 히트맵 (주 12 스냅샷) |

## 실행 방법 (How to Run)

### mrgsolve 모델

```r
library(mrgsolve)
source("rcc_mrgsolve_model.R")
# 위 스크립트가 자동으로 컴파일 → 시뮬레이션 → 플롯 생성
```

### Shiny 앱

```r
library(shiny)
shiny::runApp("rcc_shiny_app.R")
```

### Graphviz 렌더링

```bash
dot -Tsvg rcc_qsp_model.dot -o rcc_qsp_model.svg
dot -Tpng -Gdpi=150 rcc_qsp_model.dot -o rcc_qsp_model.png
```

## 임상 데이터 보정 요약 (Clinical Calibration Summary)

| 임상시험 | 요법 | mPFS (실제) | 모델 최적화 |
|----------|------|------------|------------|
| CheckMate 214 | Nivo+Ipi | 11.6 mo | PD-1 occ, CD8 kill |
| KEYNOTE-426 | Pembro+Axitinib | 15.1 mo | Axitinib IC50, 이중 차단 |
| CheckMate 9ER | Cabo+Nivo | 16.6 mo | Cabo MET/VEGFR IC50 |
| CLEAR | Len+Pembro | 23.9 mo | (lenvatinib 별도 미모델링) |
| METEOR | Cabozantinib | 7.4 mo | Cabo 단독 용량반응 |
| RECORD-1 | Everolimus | 4.0 mo | mTOR IC50 보정 |
| LITESPARK-005 | Belzutifan | ORR 22% | HIF-2α 억제 IC50 |
