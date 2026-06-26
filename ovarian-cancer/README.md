# 난소암 (Ovarian Cancer) QSP 모델

> **분류**: 부인종양학(Gynecologic Oncology) | **약어**: OC (HGSOC) | **날짜**: 2026-06-24

[![QSP 기계론적 지도](oc_qsp_model.png)](oc_qsp_model.svg)

---

## 질환 개요

난소암은 여성 암 사망의 약 5위를 차지하는 부인과 암으로, 연간 전 세계 약 32만 명이 진단된다. 가장 흔하고 치명적인 아형은 **고등급 장액성 난소암(HGSOC, High-Grade Serous Ovarian Carcinoma)**으로 전체 난소암의 약 70%를 차지한다. HGSOC는 다음 특징으로 정의된다:

- **TP53 변이**: 96% 이상에서 확인되는 거의 필수적인 변이
- **HRD(Homologous Recombination Deficiency)**: 약 50%에서 확인 (BRCA1/2 변이 15%+7%, 기타 HRR 유전자 변이 포함)
- **진단 지연**: 70-80%가 FIGO III/IV 병기에서 진단
- **복막 전이**: 종양 세포가 복강 내로 탈락하여 복막·대망(omentum)에 착상
- **백금 감수성**: 초기 치료에 양호하게 반응하나 내성 발생이 흔함

### 치료 패러다임 전환

2018-2019년 PARP 억제제 임상시험 결과(SOLO-1, PRIMA, PAOLA-1)로 유지요법이 표준화되어, HRD 상태에 따른 층화 치료가 핵심이 되었다.

---

## 10대 기계론적 서브시스템

| # | 클러스터 | 핵심 구성요소 |
|---|---------|------------|
| **1** | **DDR/HRR** — BRCA-HRR 경로 | BRCA1/2, RAD51, PARP1/2, ATM/ATR, CHK1/2, HRD score, NER, NHEJ |
| **2** | **PI3K/AKT/mTOR** — 증식·생존 | PIK3CA, PTEN, AKT, mTORC1/2, S6K1, ERK, RAS/RAF/MEK, CDK4/6, RB1, E2F |
| **3** | **VEGF/혈관신생** | HIF-1α, VEGF-A/B/C, VEGFR1/2, 내피세포, 주피세포, DLL4/Notch, Bevacizumab |
| **4** | **종양 미세환경(TME)** | CAF, TAM(M1/M2), MDSC, NK, CD8+ T, Treg, IL-6, TGF-β, IL-10, MMP-2/9, LPA, STAT3 |
| **5** | **면역회피** | PD-L1/PD-1, CTLA-4, IDO1, LAG-3, TIM-3, TIGIT, FoxP3, TLS, Pembrolizumab |
| **6** | **복막 전이** | 원발 종양, 탈락, 구상체, 복막세포, 대망, CA-125/MUC16, HE4, EMT, FIGO 병기 |
| **7** | **백금계 항암제 PK/PD** | 카보플라틴(Calvert AUC), 파클리탁셀(3구획), Pt-DNA 부가물, G2/M 정지, MDR1, GST-π |
| **8** | **PARP 억제제 PK/PD** | 오라파립(300mg BID), 니라파립(300mg QD), PARP 트래핑, 합성 치사, BRCA 역변이 내성 |
| **9** | **종양 성장/세포 생물학** | Gompertz 성장, CSC(ALDH1+), BCL-2/BAX, 카스파제 캐스케이드, Wnt/Notch, c-Myc |
| **10** | **임상 엔드포인트** | CA-125, HE4, ROMA, PFS, OS, RECIST 1.1, PFI, ctDNA, HRD 검사, BRCA 검사 |

---

## 18구획 ODE 모델

| # | 구획 (기호) | 설명 |
|---|-----------|------|
| 1 | `CAR_C1` | 카보플라틴 중심 구획 |
| 2 | `CAR_C2` | 카보플라틴 말초 구획 |
| 3 | `PAC_C1` | 파클리탁셀 중심 구획 |
| 4 | `PAC_C2` | 파클리탁셀 말초 구획 |
| 5 | `PAC_C3` | 파클리탁셀 깊은 말초 구획 |
| 6 | `OLA_gut` | 오라파립 위장관 흡수 구획 |
| 7 | `OLA_C1` | 오라파립 중심 구획 (Cmax≈5µM) |
| 8 | `OLA_C2` | 오라파립 말초 구획 |
| 9 | `NIRA_C1` | 니라파립 중심 구획 (t½≈36h) |
| 10 | `NIRA_C2` | 니라파립 말초 구획 |
| 11 | `BEV_C1` | 베바시주맙 중심 구획 |
| 12 | `BEV_C2` | 베바시주맙 말초 구획 |
| 13 | `VEGF` | 유리 VEGF-A 농도 (ng/mL) |
| 14 | `TV` | 종양 부피 (cm³, Gompertz 모델) |
| 15 | `CA125` | CA-125 혈청 (U/mL) |
| 16 | `Pt_DNA` | 백금-DNA 부가물 (상대값) |
| 17 | `CD8T` | CD8+ T세포 (상대값) |
| 18 | `HRD` | PARP 억제제 HRD 손상 축적 (0-1) |

---

## 6가지 치료 시나리오 — 2년 시뮬레이션

| # | 시나리오 | 약물 | 근거 임상시험 | 적응증 |
|---|---------|------|------------|--------|
| **S1** | 무치료 | — | 자연 경과 | — |
| **S2** | 카보플라틴+파클리탁셀 ×6사이클 | Carbo AUC6 + Pacli 175mg/m² | ICON3 (Parmar 2003 Lancet) | 표준 1차 |
| **S3** | Carbo+Pacli+베바시주맙 → 베바시주맙 유지 | +Bev 15mg/kg q3w | ICON7/GOG218 | 고위험 1차 |
| **S4** | Carbo+Pacli → 오라파립 유지 | Ola 300mg BID 2년 | **SOLO-1** (mPFS NR) | **BRCA 변이** |
| **S5** | Carbo+Pacli → 니라파립 유지 | Nira 200-300mg QD | **PRIMA** (mPFS 13.8mo HRD+) | **HRD 양성** |
| **S6** | Carbo+Pacli+Bev → 오라파립+Bev 유지 | Ola+Bev maint | **PAOLA-1** (mPFS 22.1mo HRD+) | **HRD+, Bev 적합** |

---

## 주요 파라미터 보정

| 파라미터 | 값 | 출처 |
|---------|-----|------|
| 카보플라틴 CL | GFR×0.134+0.00571×BW (L/h) | Chatelut 1995 JNCI |
| 파클리탁셀 CL | 13.2 L/h (비선형 PK) | Gianni 1995 JCO |
| 오라파립 t½ | 11.9h (300mg BID) | Doherty 2014 Clin Pharmacokinet |
| 니라파립 t½ | 36h (QD 투여) | Sandhu 2013 JCO |
| 베바시주맙 t½ | ~20일 (IgG1) | Lu 2008 Cancer Chemother Pharmacol |
| CA-125 t½ | ~23일 (혈청 반감기) | Rustin 1996 JCO |
| 종양 배가시간 | ~60일 (무치료, Gompertz) | Oza 2015 Lancet Oncol |
| SOLO-1 mPFS | NR vs 13.8mo (HR 0.30, BRCA+) | Moore 2018 NEJM |
| PRIMA mPFS | 13.8mo vs 8.2mo (HR 0.43, HRD+) | Gonzalez-Martin 2019 NEJM |
| PAOLA-1 mPFS | 22.1mo vs 16.6mo (HR 0.33, HRD+) | Ray-Coquard 2019 NEJM |

---

## Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| **① 환자 프로파일** | BRCA 상태, HRD 점수, 초기 CA-125, GFR, FIGO 병기, 치료 적합성 매트릭스 |
| **② 약물 PK** | 카보플라틴·파클리탁셀·오라파립·니라파립·베바시주맙 시간-농도 곡선 |
| **③ PD 바이오마커** | CA-125 동역학, 백금-DNA 부가물, HRD 손상 축적, CD8+ T세포 침윤 |
| **④ 종양 반응** | 종양 부피 Gompertz 곡선, RECIST 분류, 최선 반응%, 추정 PFS |
| **⑤ 시나리오 비교** | 6가지 치료 시나리오 종양·CA-125 비교, 요약 테이블 |
| **⑥ 바이오마커 패널** | 종합 바이오마커 6개 패널, BRCA/HRD 치료 결정 트리, 임상시험 참조 수치 |

---

## QSP 모델 파일

| 구성요소 | 파일 | 사양 |
|---------|------|------|
| 🗺️ 기계론적 지도 (DOT) | [`oc_qsp_model.dot`](oc_qsp_model.dot) | **180+ 노드, 10클러스터** (fdp 레이아웃) |
| 🖼️ 기계론적 지도 (SVG) | [`oc_qsp_model.svg`](oc_qsp_model.svg) | 벡터 이미지 |
| 🖼️ 기계론적 지도 (PNG) | [`oc_qsp_model.png`](oc_qsp_model.png) | 150 dpi 래스터 이미지 |
| ⚙️ mrgsolve ODE 모델 | [`oc_mrgsolve_model.R`](oc_mrgsolve_model.R) | **18구획 ODE**, **6치료 시나리오**, SOLO-1/PRIMA/PAOLA-1 보정 |
| 📊 Shiny 대시보드 | [`oc_shiny_app.R`](oc_shiny_app.R) | **6탭** 인터랙티브 대시보드 |
| 📚 참고문헌 | [`oc_references.md`](oc_references.md) | **55개 PubMed 인용** (14개 섹션) |

---

## Shiny 앱 실행 방법

```r
# 필수 패키지 설치
install.packages(c("shiny", "shinydashboard", "mrgsolve", "dplyr",
                   "ggplot2", "tidyr", "DT", "plotly", "patchwork"))

# 앱 실행
shiny::runApp("ovarian-cancer/oc_shiny_app.R")
```

---

*난소암 QSP 모델 | Ovarian Cancer (HGSOC) QSP | 2026-06-24 | Claude Code Routine*
