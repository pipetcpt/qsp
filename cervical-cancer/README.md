# 자궁경부암 (Cervical Cancer) QSP 모델

> **분류**: 부인종양학(Gynecologic Oncology) | **약어**: CC | **날짜**: 2026-07-01

[![QSP 기계론적 지도](cc_qsp_model.png)](cc_qsp_model.svg)

---

## 질환 개요

자궁경부암은 전 세계 여성암 4위를 차지하며, 연간 약 66만 명이 진단되고 35만 명이 사망한다(GLOBOCAN 2020). 사실상 모든 증례가 **고위험군 인유두종바이러스(High-risk HPV, 특히 HPV16/18)**의 지속 감염에서 기인하며, 다음의 병태생리적 특징을 갖는다:

- **HPV E6/E7 종양유전자**: E6은 E6AP를 매개로 p53을 유비퀴틴-프로테아좀 경로로 분해시켜 G1 체크포인트를 소실시키고, E7은 pRb에 결합해 E2F를 방출함으로써 세포주기를 통제 불능 상태로 만든다.
- **바이러스 유전체 통합(Integration)**: 에피솜 상태에서 숙주 유전체로 통합되며 E2 억제인자가 소실되어 E6/E7이 과발현된다.
- **CIN 진행 연속체**: CIN1(경도)→CIN2(중등도)→CIN3/CIS(고도/상피내암)를 거쳐 기저막을 침범하는 침윤암으로 진행하며, p16INK4a 과발현이 대리 표지자로 활용된다.
- **선별검사로 예방 가능**: Pap도말세포검사·HPV 병행검사와 HPV 백신 접종으로 발생률을 획기적으로 낮출 수 있는 대표적 예방 가능 암종이다.
- **진단 시 병기 분포**: 선별검사 미비 지역에서는 국소진행성(FIGO IB3-IVA) 상태로 진단되는 경우가 많아 동시화학방사선(CCRT)이 표준치료의 근간이 된다.

### 치료 패러다임 전환

2014년 GOG-240(Tewari, NEJM)이 재발/전이성 질환에서 베바시주맙의 생존 이득을 입증한 이후, 2021년 KEYNOTE-826(Colombo, NEJM)이 PD-L1 CPS≥1 1차 치료에 펨브롤리주맙을 편입시켰고, 2024년 KEYNOTE-A18(Lorusso, Lancet)이 국소진행성 고위험군의 동시화학방사선에 펨브롤리주맙을 추가하여 표준치료 패러다임을 근본적으로 바꾸었다. 2차 이상 재발/전이성 질환에서는 조직인자(Tissue Factor) 표적 항체-약물 접합체 **티소투맙 베도틴**(innovaTV 204/301)이 새로운 치료 축으로 자리잡았다.

---

## 11대 기계론적 서브시스템

| # | 클러스터 | 핵심 구성요소 |
|---|---------|------------|
| **1** | **HPV 종양유전 경로** | HPV16/18 E6/E7, E6AP, p53 분해, pRb-E2F 축, 바이러스 통합, CIN1/2/3, p16INK4a, 텔로머라제 |
| **2** | **세포주기/증식/아폽토시스 회피** | Cyclin E/D1, CDK2/4/6, p21, p27, PIK3CA, AKT/mTOR, BCL-2/BAX, Survivin, 캐스파제 |
| **3** | **종양 미세환경/면역회피** | PD-L1/PD-1, CD8+T, Treg, MDSC, TAM, IDO1, MHC-I 하향조절(E7 매개), T세포 소진 |
| **4** | **혈관신생** | HIF-1α, VEGF-A, VEGFR1/2, 종양 저산소증, 베바시주맙 |
| **5** | **침습/전이/병기** | EMT, MMP-2/9, 림프혈관침범(LVSI), 골반·대동맥주위 림프절, FIGO 2018 병기 |
| **6** | **백금 기반 동시화학방사선 PK/PD** | 시스플라틴(주간 40mg/m²), 백금-DNA 부가물, NER, 방사선 증감 |
| **7** | **방사선치료** | EBRT+근접치료, LQ 모델(α/β=10Gy), TCP, 재산소화, 가속 재증식 |
| **8** | **면역관문억제제 PK/PD** | 펨브롤리주맙, PD-1 수용체 점유, T세포 재활성화, CPS 점수 |
| **9** | **항체-약물 접합체 PK/PD** | 티소투맙 베도틴, 조직인자(TF) 표적, MMAE 방출, 미세소관 붕괴, 방관자 살해 |
| **10** | **임상 엔드포인트/바이오마커** | SCC-Ag, HPV 바이러스 부하, RECIST 1.1, PFS, OS, 병리학적 완전관해 |
| **11** | **선별검사/예방(맥락)** | HPV 백신, Pap/HPV 병행검사, 콜포스코피, LEEP/원추절제술 |

---

## 19구획 ODE 모델

| # | 구획 (기호) | 설명 |
|---|-----------|------|
| 1 | `CIS_C1` | 시스플라틴 중심 구획 |
| 2 | `CIS_C2` | 시스플라틴 말초 구획 |
| 3 | `PAC_C1` | 파클리탁셀 중심 구획 (재발/전이성 간헐 화학요법) |
| 4 | `PAC_C2` | 파클리탁셀 말초 구획 |
| 5 | `BEV_C1` | 베바시주맙 중심 구획 |
| 6 | `BEV_C2` | 베바시주맙 말초 구획 |
| 7 | `PEMBRO_C1` | 펨브롤리주맙 중심 구획 |
| 8 | `PEMBRO_C2` | 펨브롤리주맙 말초 구획 |
| 9 | `TV_ADC_C1` | 티소투맙 베도틴(ADC) 중심 구획 |
| 10 | `TV_ADC_C2` | 티소투맙 베도틴(ADC) 말초 구획 |
| 11 | `MMAE_free` | 유리 종양내 MMAE 페이로드 (상대값) |
| 12 | `VEGF` | 유리 VEGF-A 농도 (ng/mL) |
| 13 | `Pt_DNA` | 백금-DNA 부가물 (상대값, 0-1) |
| 14 | `RT_SF` | 누적 방사선 손상 (LQ 모델, 상대값) |
| 15 | `TV` | 종양 부피 (cm³, Gompertz 모델) |
| 16 | `SCCAg` | SCC-Ag 혈청 (ng/mL) |
| 17 | `HPVload` | HPV 바이러스 부하 (상대 log10 copies) |
| 18 | `CD8T` | CD8+ T세포 (상대값) |
| 19 | `PDL1_exp` | 종양 PD-L1 발현 (상대값, CPS 유사) |

---

## 6가지 치료 시나리오 — 2년 시뮬레이션

| # | 시나리오 | 약물/방사선 | 근거 임상시험 | 적응증 |
|---|---------|------|------------|--------|
| **S1** | 무치료 | — | 자연 경과 | — |
| **S2** | 시스플라틴 동시화학방사선(CCRT) | Cisplatin 40mg/m² 주간×6 + EBRT/근접치료 | **RTOG-90-01** (Morris 1999/Eifel 2004 NEJM/JCO) | 국소진행성 표준 |
| **S3** | CCRT + 펨브롤리주맙 (동시+유지) | +Pembro 200mg q3w | **KEYNOTE-A18** (Lorusso 2024 Lancet) | 국소진행성 고위험 |
| **S4** | 백금+파클리탁셀+베바시주맙 | Cis/Carbo+Pacli+Bev 15mg/kg q3w | **GOG-240** (Tewari 2014 NEJM) | 재발/전이성 1차 |
| **S5** | 티소투맙 베도틴 단독 | TV-ADC 2.0mg/kg q3w | **innovaTV 301** (Vergote 2024 NEJM) | 재발/전이성 2차 이상 |
| **S6** | 백금+파클리탁셀+베바시주맙+펨브롤리주맙 | 4제 병용 | **KEYNOTE-826** (Colombo 2021 NEJM) | 재발/전이성 1차, CPS≥1 |

---

## 주요 파라미터 보정

| 파라미터 | 값 | 출처 |
|---------|-----|------|
| 시스플라틴 표준 용량 | 40mg/m² 주간 IV ×5-6 (CCRT 동시) | Rose 1999 NEJM (GOG-120) |
| 시스플라틴 CL | ~30 L/h (총 백금 기준) | Reece 1987 J Clin Oncol |
| 방사선 LQ α/β | 10 Gy (자궁경부 편평상피암) | Fowler 1989 Br J Radiol |
| EBRT + 근접치료 총선량 | EQD2 ≥85Gy (HR-CTV) | Pötter 2018/2021 (EMBRACE-I/II) |
| 베바시주맙 t½ | ~20일 (IgG1) | Lu 2008 Cancer Chemother Pharmacol |
| 펨브롤리주맙 t½ | ~22일, 선형 PK | Ahamadi 2017 CPT:PSP |
| 티소투맙 베도틴 용량 | 2.0mg/kg q3w IV | Coleman 2021 Lancet Oncol (innovaTV 204) |
| SCC-Ag t½ | ~2.8일 | Gaarenstroom 2000 Int J Gynecol Cancer |
| RTOG-90-01 결과 | 동시화학방사선 OS 우위 | Eifel 2004 J Clin Oncol |
| GOG-240 mOS | 16.8 vs 13.3개월 (HR 0.71) | Tewari 2014 NEJM |
| KEYNOTE-826 결과 | mOS 개선, CPS≥1군 최대 이득 | Colombo 2021 NEJM |
| KEYNOTE-A18 mPFS | HR 0.70 (고위험 국소진행성) | Lorusso 2024 Lancet |
| innovaTV 301 mOS | 11.5 vs 9.5개월 (HR 0.70) | Vergote 2024 NEJM |

---

## Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| **① 환자 프로파일** | FIGO 병기, PD-L1 CPS 상태, 초기 SCC-Ag, CCRT 여부, 치료 적합성 매트릭스 |
| **② 약물 PK** | 시스플라틴·파클리탁셀·펨브롤리주맙·티소투맙 베도틴·베바시주맙 시간-농도 곡선 |
| **③ PD 주요지표** | 백금-DNA 부가물, 누적 방사선 손상(LQ), CD8+ T세포, 종양 PD-L1 발현 |
| **④ 임상 엔드포인트** | 종양 부피 Gompertz 곡선, RECIST 분류, 최선 반응%, 추정 PFS |
| **⑤ 시나리오 비교** | 6가지 치료 시나리오 종양·SCC-Ag 비교, 요약 테이블 |
| **⑥ 바이오마커** | 종합 바이오마커 6개 패널, 병기별 치료 알고리즘, 임상시험 참조 수치 |

---

## QSP 모델 파일

| 구성요소 | 파일 | 사양 |
|---------|------|------|
| 🗺️ 기계론적 지도 (DOT) | [`cc_qsp_model.dot`](cc_qsp_model.dot) | **150+ 노드, 11클러스터** (fdp 레이아웃) |
| 🖼️ 기계론적 지도 (SVG) | [`cc_qsp_model.svg`](cc_qsp_model.svg) | 벡터 이미지 |
| 🖼️ 기계론적 지도 (PNG) | [`cc_qsp_model.png`](cc_qsp_model.png) | 150 dpi 래스터 이미지 |
| ⚙️ mrgsolve ODE 모델 | [`cc_mrgsolve_model.R`](cc_mrgsolve_model.R) | **19구획 ODE**, **6치료 시나리오**, RTOG-90-01/GOG-240/KEYNOTE-A18/826/innovaTV 301 보정 |
| 📊 Shiny 대시보드 | [`cc_shiny_app.R`](cc_shiny_app.R) | **6탭** 인터랙티브 대시보드 |
| 📚 참고문헌 | [`cc_references.md`](cc_references.md) | **63개 PubMed 인용** (15개 섹션) |

---

## Shiny 앱 실행 방법

```r
# 필수 패키지 설치
install.packages(c("shiny", "shinydashboard", "mrgsolve", "dplyr",
                   "ggplot2", "tidyr", "DT", "plotly", "patchwork"))

# 앱 실행
shiny::runApp("cervical-cancer/cc_shiny_app.R")
```

---

*자궁경부암 QSP 모델 | Cervical Cancer (HPV-driven) QSP | 2026-07-01 | Claude Code Routine*
