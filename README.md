# qsp

## mrgsolve

- <https://vantage-research.net/qsp-in-r/>
- gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve
    - <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/>
    - <https://github.com/Genentech/gPKPDviz/>
    



## iqrtools

- <https://www.intiquan.com/acop2019_qsp/>

---

## QSP 질환 모델 라이브러리 — Alpha-1 Antitrypsin Deficiency (AATD)

> 알파-1 항트립신 결핍증 QSP 모델 (2026-06-24 추가)

[![Disease](https://img.shields.io/badge/Disease-Alpha--1%20Antitrypsin%20Deficiency-red)]()
[![Nodes](https://img.shields.io/badge/Map-130%2B%20nodes%20%7C%2010%20clusters-blue)]()
[![ODE](https://img.shields.io/badge/ODE-20%20compartments%20%7C%206%20scenarios-green)]()

### 질환 개요

**SERPINA1 Glu342Lys(Z 대립유전자) → Z-AAT 소포체 내 루프-시트 중합체 축적 → 간(기능획득 독성) + 폐(기능소실) 이중 병리**

- **간 경로**: Z-AAT 중합체 → ER 스트레스/UPR → HSC 활성화(TGF-β1) → 섬유화 F0→F4 → 간경변 → HCC(10-40배↑)
- **폐 경로**: 혈청 AAT < 11 µM → 폐포 NE 무제한 활성 + MMP-12 → 탄력소 분해 → 범소엽성 폐기종 → COPD (FEV1 50-200 mL/yr 감소)

### 모델 파일

| 파일 | 설명 |
|------|------|
| [`aatd_qsp_model.dot`](alpha1-antitrypsin-deficiency/aatd_qsp_model.dot) · [`.svg`](alpha1-antitrypsin-deficiency/aatd_qsp_model.svg) · [`.png`](alpha1-antitrypsin-deficiency/aatd_qsp_model.png) | 기계론적 지도 (130+ 노드, 10 클러스터) |
| [`aatd_mrgsolve_model.R`](alpha1-antitrypsin-deficiency/aatd_mrgsolve_model.R) | mrgsolve ODE 모델 (20구획, 6치료 시나리오) |
| [`aatd_shiny_app.R`](alpha1-antitrypsin-deficiency/aatd_shiny_app.R) | Shiny 대시보드 (6탭) |
| [`aatd_references.md`](alpha1-antitrypsin-deficiency/aatd_references.md) | 참고문헌 (54편, 13섹션) |
| [`README.md`](alpha1-antitrypsin-deficiency/README.md) | 상세 모델 설명 |

### 기계론적 지도

[![AATD QSP Map](alpha1-antitrypsin-deficiency/aatd_qsp_model.png)](alpha1-antitrypsin-deficiency/aatd_qsp_model.svg)

### 6가지 치료 시나리오 (5년 시뮬레이션)

| 시나리오 | 치료 | 임상 근거 |
|---------|------|---------|
| S1 | 무치료 (자연 경과) | Tanash 2010, Janus 1985 |
| S2 | AAT 보충요법 — Prolastin-C 60 mg/kg/wk IV | RAPID Trial (Chapman 2015 *Lancet*) |
| S3 | Fazirsiran siRNA — 200 mg SQ q12wk | ARO-AAT Phase 2 (Strnad 2022 *NEJM*) |
| S4 | NE 억제제 — Alvelestat 60 mg BID | Phase 2 (McElvaney 2020 *AJRCCM*) |
| S5 | 유전자치료 — rAAV-SERPINA1 단회 | Flotte 2004, Mueller 2008 |
| S6 | 보충요법 + NE 억제제 조합 | 시뮬레이션 기반 예측 |
