# Uterine Leiomyoma (Fibroid) QSP Model
# 자궁근종 정량적 시스템 약리학 모델

[![Mechanistic Map](ufl_qsp_model.png)](ufl_qsp_model.svg)

---

## 개요 (Overview)

**질환명:** 자궁근종 (Uterine Leiomyoma / Uterine Fibroids)  
**유병률:** 가임기 여성의 70–80% (50세 이전 누적 발생률)  
**주요 증상:** 과다 월경출혈(HMB), 골반통, 골반 압박감, 불임  
**약물 표적:** GnRH 수용체, 에스트로겐/프로게스테론 수용체, ECM 리모델링

---

## 모델 구성 (Model Components)

### 기계론적 지도 (Mechanistic Map)
15개 클러스터, 120+ 노드의 상호작용 경로 포함:

| 클러스터 | 내용 |
|---------|------|
| HPG 축 | GnRH → LH/FSH 분비, KNDy 신경세포, 뇌하수체 탈감작 |
| 난소 스테로이드 생성 | 콜레스테롤 → E2/P4 생합성 (CYP11A1, CYP17A1, CYP19A1) |
| 호르몬 피드백 | E2/P4 → GnRH 음성피드백, 인히빈 → FSH 억제 |
| 자궁 생물학 | 자궁근층/내막 ERα/PR 신호, 자궁혈관 VEGF |
| 근종 발병기전 | MED12 돌연변이, HMGA2, 아로마타제 과발현, ECM 축적 |
| 세포 내 신호전달 | MAPK/ERK, PI3K/AKT/mTOR, Wnt/β-catenin, NF-κB |
| ECM 리모델링 | 콜라겐 I/III, MMP/TIMP 균형, LOX 교차결합 |
| 염증/면역 | M2 대식세포, 비만세포, PGE2, COX-2 |
| GnRH 작용제 PK | 류프로라이드 데포 흡수-분포-제거 |
| GnRH 길항제 PK | 엘라고릭스/렐루고릭스 경구 PK |
| SPRM PK | 울리프리스탈 (UPA) PK/PD |
| 골 건강 | BMD, RANKL/OPG, 저에스트로겐 골흡수 |
| 임상 엔드포인트 | MBL, PBAC 점수, 헤모글로빈, UFS-QoL |
| 기타 치료 | LNG-IUD, OCP, 트라넥사믹산, 수술, UAE, MRgFUS |
| 위험인자 | 인종, 초경 연령, 미산부, 비만, 가족력 |

### mrgsolve ODE 모델
**18개 ODE 구획:**

| 구획 | 설명 |
|-----|------|
| GnRH_C | GnRH 펄스 농도 |
| LH_C, FSH_C | 황체화/난포자극호르몬 |
| E2_C, P4_C | 에스트라디올, 프로게스테론 |
| V_fib, ECM_fib | 근종 용적 (세포성 + ECM) |
| MBL_cum | 주기당 월경혈량 |
| Hgb_C | 헤모글로빈 |
| BMD_C | 골밀도 (정규화) |
| Leu_depot/plasma | 류프로라이드 데포 PK |
| Ela_gut/plasma | 엘라고릭스 경구 PK |
| Rel_gut/plasma | 렐루고릭스 경구 PK |
| UPA_gut/plasma | UPA 경구 PK |

### 치료 시나리오 (6 Scenarios)

| 시나리오 | 치료 | 임상시험 근거 |
|---------|------|-------------|
| S1 | 무치료 (자연 경과) | — |
| S2 | 류프로라이드 3.75mg 데포 (q4w, 24주) | Friedman 1989 Fertil Steril |
| S3 | 엘라고릭스 150mg BID (24주) | ELARIS UF-I, Simon 2020 NEJM |
| S4 | 엘라고릭스 200mg BID + 호르몬 보충 | ELARIS UF-I/II, Simon/Schlaff NEJM 2020 |
| S5 | 렐루고릭스 복합제 (40mg+E2/NET) QD | LIBERTY 1, Lukes 2021 NEJM |
| S6 | 울리프리스탈 5mg QD (13주×2 코스) | PEARL I/II, Donnez 2012 NEJM |

---

## 주요 임상시험 결과 (Key Clinical Trial Results)

| 임상시험 | 약물 | 1차 엔드포인트 달성률 | 참고문헌 |
|---------|------|-------------------|---------|
| ELARIS UF-I | Elagolix 200mg BID+AB | **68.5%** (Week 24) | Simon JA, NEJM 2020;382:328 |
| ELARIS UF-II | Elagolix 200mg BID+AB | **76.5%** (Week 24) | Schlaff WD, NEJM 2020;382:317 |
| LIBERTY 1 | Relugolix combination | **71.2%** (Week 24) | Lukes AS, NEJM 2021;384:630 |
| LIBERTY 2 | Relugolix combination | **70.6%** (Week 24) | Al-Hendy A, NEJM 2021;384:630 |
| PRIMROSE 1 | Linzagolix 200mg+AB | **93.9%** (Week 24) | Murji A, NEJM 2022;387:1767 |
| PEARL I | UPA 5mg×13wk | **91%** controlled bleeding | Donnez J, NEJM 2012;366:409 |

*1차 엔드포인트: 월경혈량 < 80 mL/cycle (HMB 정의) AND ≥ 50% 감소*

---

## 핵심 약동학 파라미터 (Key PK Parameters)

| 약물 | 용량 | t½ | 생체이용률 | Tmax |
|-----|------|----|---------|------|
| Leuprolide depot | 3.75 mg Q4W | 3–4주 (데포) | ~95% | 3–4주 (지속) |
| Elagolix | 150/200 mg BID | 4–6 h | 56% | ~1 h |
| Relugolix | 40 mg QD | ~60 h | 12% | ~2 h |
| UPA | 5 mg QD | 32–38 h | 87% | ~1 h |

---

## Shiny 앱 구성 (Shiny App Tabs)

| 탭 | 내용 |
|----|------|
| ① 환자 프로파일 | 초기 파라미터, 위험인자 매트릭스, 치료 결정 가이드 |
| ② 약물 PK | 혈중 약물 농도-시간 그래프, PK 파라미터 테이블 |
| ③ PD 주요 지표 | E2, P4, LH/FSH, 근종 용적 동역학 |
| ④ 임상 엔드포인트 | MBL, 헤모글로빈, BMD, 안면홍조 점수 |
| ⑤ 시나리오 비교 | 6가지 치료 시나리오 전체 비교 그래프 & 요약 테이블 |
| ⑥ 바이오마커 패널 | 임상시험 결과, PBAC 점수, 치료 반응 지표 |

---

## 파일 목록 (Files)

| 파일 | 설명 |
|------|------|
| `ufl_qsp_model.dot` | Graphviz 기계론적 지도 소스 |
| `ufl_qsp_model.svg` | 기계론적 지도 SVG (벡터, 고해상도) |
| `ufl_qsp_model.png` | 기계론적 지도 PNG (150 dpi) |
| `ufl_mrgsolve_model.R` | mrgsolve ODE 모델 + 6가지 시나리오 + 시각화 |
| `ufl_shiny_app.R` | Shiny 인터랙티브 대시보드 (6탭) |
| `ufl_references.md` | 참고문헌 60개 (PubMed 링크 포함) |
| `README.md` | 이 파일 |

---

## 실행 방법 (How to Run)

```r
# mrgsolve 모델 실행
source("ufl_mrgsolve_model.R")

# Shiny 앱 실행
shiny::runApp("ufl_shiny_app.R")
```

### 필요 R 패키지

```r
install.packages(c("mrgsolve", "ggplot2", "dplyr", "tidyr",
                   "shiny", "shinydashboard", "plotly", "DT"))
```

---

## 질환 개요 (Disease Overview)

### 발병기전 핵심

```
HPG 축 활성
    ↓
E2 (에스트라디올) / P4 (프로게스테론) 분비
    ↓
자궁근층 내 ERα/PR 과발현 (MED12/HMGA2 돌연변이)
    ↓
세포 증식 ↑ (MAPK/ERK, PI3K/AKT, Wnt/β-catenin)
+ ECM 축적 ↑ (TGF-β → 콜라겐 I/III)
+ 아로마타제 과발현 → 국소 E2 생산 (양성 피드백)
    ↓
자궁근종 형성 및 성장
    ↓
AUB (과다 월경출혈) + 골반통 + 골반 압박감
    ↓
철 결핍 빈혈 + 불임 + 삶의 질 저하
```

### GnRH 길항제 vs 작용제 비교

| 특성 | GnRH 작용제 (류프로라이드) | GnRH 길항제 (엘라고릭스/렐루고릭스) |
|-----|------------------------|--------------------------------|
| 작용 기전 | 뇌하수체 GnRHR 하향조절 (지속 노출) | 경쟁적 GnRHR 차단 (즉각적) |
| 초기 반응 | 자극 효과 (flare) 1–2주 | 즉각적 억제 (flare 없음) |
| LH/FSH 억제 시작 | 2–4주 후 억제 | 수일 내 억제 |
| 투여 경로 | 주사 (데포) | 경구 |
| 가역성 | 완전 회복 (중단 후 3–6개월) | 빠른 회복 (단기 t½) |
| 호르몬 보충 필요성 | 6개월 이상 시 필요 | 장기 치료 시 필요 |

---

*모델 생성일: 2026-06-25 | QSP Library CCR Auto-generated*  
*참고문헌: 60개 | 구획 수: 18개 ODE | Shiny 탭: 6개*
