# Glioblastoma Multiforme (GBM) — QSP Model

[![Disease](https://img.shields.io/badge/Disease-Glioblastoma-red)](.) [![Category](https://img.shields.io/badge/Category-CNS%20Oncology-blueviolet)](.) [![Drugs](https://img.shields.io/badge/Drugs-TMZ%20%C2%B7%20Bevacizumab%20%C2%B7%20Pembrolizumab%20%C2%B7%20TTF-blue)](.) [![Status](https://img.shields.io/badge/Status-Complete-brightgreen)](.)

## 개요 (Overview)

**교모세포종(Glioblastoma Multiforme, GBM)**은 WHO Grade IV 악성 신경교종으로, 성인에서 가장 흔하고 예후가 불량한 원발성 뇌종양이다. 중앙 생존기간은 약 14–22개월(MGMT 메틸화 여부에 따라)이며, 5년 생존율은 5% 미만이다.

이 모델은 **Stupp 프로토콜** (TMZ + 방사선치료 → 보조 TMZ), **베바시주맙**, **펨브롤리주맙**, **TTF** 등 주요 치료 전략을 망라하는 포괄적 QSP 모델을 제공한다.

---

## 기계론적 지도 (Mechanistic Map)

[![GBM Mechanistic Map](gbm_qsp_model.png)](gbm_qsp_model.svg)

*전체 해상도 SVG를 보려면 이미지를 클릭하세요.*

| 구성 | 내용 |
|------|------|
| **총 노드** | 220+ |
| **서브그래프 클러스터** | 13 |
| **주요 경로** | EGFR→PI3K/AKT/mTOR, RAS-MAPK, p53-RB, TMZ-MGMT-MMR, VEGF-혈관신생, 면역관문, GSC 줄기세포 유지 |

### 클러스터 목록

| # | 클러스터 | 주요 구성 요소 |
|---|---------|--------------|
| 1 | 유전체 변이 | IDH1/2, MGMT, EGFR amp, PTEN 소실, TP53, TERT, CDKN2A 결실 |
| 2 | RTK/RAS-MAPK | EGFR, PDGFR, MET, GRB2-SOS, RAS-RAF-MEK-ERK |
| 3 | PI3K/AKT/mTOR | PI3K, PIP3, PTEN, AKT, TSC1/2, mTORC1/2, S6K1, 4E-BP1 |
| 4 | 세포주기 & 아포토시스 | CDK4/6-CycD, RB-E2F, p16/p21, BCL-2 family, Caspase cascade |
| 5 | DNA 손상 & 복구 | O6-MeG, MGMT, MMR, DSB, ATM/ATR-CHK1/2, NHEJ/HR |
| 6 | 종양 미세환경 | TAM-M1/M2, MDSCs, CAFs, TGF-β, IL-6/STAT3, IDO1, ARG1 |
| 7 | 혈관신생 | HIF-1α, VEGF-A/B/C, VEGFR1/2, ANG1/2, FGF2, 신생혈관 |
| 8 | 면역반응 & 관문 | CD8+ CTL, Treg, PD-1/PD-L1, CTLA-4, TIM-3, LAG-3, TIGIT |
| 9 | 교모세포종 줄기세포 | CD133+, SOX2, Notch, Wnt/β-cat, SHH/GLI, BMI1, EZH2 |
| 10 | 약물 PK | TMZ 2-cmt+BBB, 베바시주맙 2-cmt+VEGF 결합, 펨브롤리주맙 |
| 11 | 혈뇌장벽 | BBB tight junction, P-gp, BCRP, Kp,brain, 혈액-종양 장벽(BTB) |
| 12 | 방사선치료 | LQ 모델(α/β=10 Gy), OER, TMZ-RT 시너지, SRS/SBRT |
| 13 | 임상 엔드포인트 | 종양 부피, RANO 기준, PFS/OS, KPS, MGMT/IDH 예후인자 |

---

## mrgsolve ODE 모델

### 구획 설계 (18구획)

| 카테고리 | 구획 | 설명 |
|----------|------|------|
| **TMZ PK** | Gut, Cp_tmz, Cp2_tmz, Cbrain, O6MeG | GI 흡수→2-구획 PK→BBB 투과→O6-MeG 병변 |
| **BEV PK** | BEV_Cp, BEV_Cp2, VEGF_free, BEV_VEGF | 2-구획 PK + VEGF 결합/해리 동력학 |
| **Anti-PD1** | APD1_Cp, PD1_occ | 1-구획 PK + PD-1 수용체 점유 |
| **종양 세포** | Ts, Tr, GSC | 민감/내성/줄기세포 아집단 |
| **면역** | CD8_eff, Treg_c, TAM_M2 | CD8+ T세포, 조절 T세포, M2-TAM |
| **혈관신생** | NV | VEGF 구동 신생혈관 지수 |

### 핵심 방정식

**종양 성장 (Gompertz 모델)**
```
dTs/dt = Ts × [−kg × ln(N/K) − kill_TMZ − kill_RT − kill_CD8 − k_resist]
```

**O6-MeG 병변 (TMZ 기전)**
```
dO6MeG/dt = k_O6 × C_brain − kMGMT × O6MeG − k_O6deg × O6MeG
```
*MGMT 메틸화: kMGMT = 0.05 h⁻¹ (low repair)*
*MGMT 비메틸화: kMGMT = 0.40 h⁻¹ (high repair)*

**RT 효과 (Linear-Quadratic 모델)**
```
SF = exp(−αD − βD²)  [α=0.30 Gy⁻¹, β=0.030 Gy⁻², α/β=10 Gy]
```

### 파라미터 보정 근거

| 파라미터 | 값 | 출처 |
|----------|---|------|
| TMZ CL | 11.4 L/h | Ostermann 2004 Clin Cancer Res |
| TMZ V₁ | 22.5 L | Baker 2003 Clin Cancer Res |
| TMZ Kp,brain | 0.28 | Ostermann 2004 |
| BEV CL | 0.207 L/day | Lu 2008 Cancer Chemother Pharmacol |
| GBM α/β | 10 Gy | Joiner & van der Kogel 2009 |
| kg (Gompertz) | 0.003 day⁻¹ | Calibrated to Stupp 2005 |

### 치료 시나리오 (7가지)

| # | 시나리오 | 기반 임상시험 | 주요 기전 |
|---|---------|-------------|---------|
| S1 | 무치료 대조군 | — | Gompertz 자연 경과 |
| S2 | Stupp (MGMT 메틸화) | Stupp 2005 NEJM, Hegi 2005 NEJM | TMZ+RT→보조 TMZ; O6-MeG 축적 |
| S3 | Stupp (MGMT 비메틸화) | Hegi 2005 NEJM | MGMT 고발현→O6-MeG 신속 복구→내성 |
| S4 | Stupp + 베바시주맙 | AVAGLIO 2014, RTOG0825 | VEGF 중화→혈관신생 억제 |
| S5 | Stupp + TTF | EF-14 2017 JAMA | 전기장→분열기 세포 방해 |
| S6 | 펨브롤리주맙 + TMZ | Keynote-028, Reardon 2020 | PD-1 차단→CD8 재활성화 |
| S7 | 베바시주맙 단독 (구제) | Friedman 2009 JCO | VEGF 중화 (재발 GBM) |

---

## Shiny 대시보드

### 탭 구성 (7탭)

| 탭 | 내용 |
|---|------|
| **1. Patient Profile** | 연령·KPS·MGMT/IDH 상태·절제 범위·치료 선택 |
| **2. Drug PK** | TMZ 혈장/뇌 PK, 베바시주맙 PK, 펨브롤리주맙 PK |
| **3. DNA Damage & MGMT** | O6-MeG 병변 동태, MGMT 메틸화 비교 |
| **4. Tumor Dynamics** | Ts/Tr/GSC 세포 집단, 종양 부피 추세 |
| **5. Scenario Comparison** | 7가지 치료 전략 동시 비교 |
| **6. TME & Biomarkers** | CD8/Treg/TAM 면역 미세환경, VEGF, PD-1 점유율 |
| **7. Clinical Endpoints** | 종양 직경(RANO), 치료 kill rate 구성요소 |

### 실행 방법

```r
# R 패키지 설치
install.packages(c("shiny", "mrgsolve", "dplyr", "ggplot2",
                   "tidyr", "shinydashboard", "DT", "gridExtra"))

# 앱 실행
shiny::runApp("gbm_shiny_app.R")
```

---

## 핵심 임상 예후인자

| 인자 | 좋은 예후 | 나쁜 예후 | 중앙 OS 차이 |
|------|----------|---------|------------|
| MGMT 메틸화 | 메틸화 (+) | 비메틸화 (−) | 21.7 vs 12.6개월 (Hegi 2005) |
| IDH 돌연변이 | IDH-mutant | IDH-WT | ~31 vs ~15개월 |
| 절제 범위 | GTR | 조직검사 | ~3–4개월 차이 |
| KPS | ≥70 | <70 | 독립 예후인자 |
| 연령 | <50세 | >70세 | 연령 독립 예후인자 |

---

## 파일 목록

| 파일 | 설명 |
|------|------|
| [`gbm_qsp_model.dot`](gbm_qsp_model.dot) | Graphviz 기계론적 지도 원본 (220+ 노드, 13 클러스터) |
| [`gbm_qsp_model.svg`](gbm_qsp_model.svg) | 벡터 형식 지도 (브라우저에서 확대 가능) |
| [`gbm_qsp_model.png`](gbm_qsp_model.png) | 래스터 형식 지도 (150 dpi) |
| [`gbm_mrgsolve_model.R`](gbm_mrgsolve_model.R) | mrgsolve ODE QSP 모델 (18구획, 7 시나리오) |
| [`gbm_shiny_app.R`](gbm_shiny_app.R) | Shiny 인터랙티브 대시보드 (7탭) |
| [`gbm_references.md`](gbm_references.md) | 참고문헌 60편+ (13개 섹션) |

---

## 병태생리 요약

```
EGFR amp / PTEN loss / IDH1 WT
        ↓
PI3K/AKT/mTOR 과활성화 → 세포 증식
RAS/MAPK 과활성화 → 침윤
        ↓
종양 성장 + 혈관신생 (VEGF-A↑ → HIF-1α)
        ↓
면역억제 (M2-TAM + Treg + PD-L1↑)
        ↓
TMZ 저항 (MGMT 비메틸화) + GSC 재발
```

**치료 전략**:
- **1차**: TMZ+RT(방사선) (Stupp protocol) → 보조 TMZ × 6 사이클
- **MGMT 메틸화**: TMZ 반응 우수 (HR 0.45, Hegi 2005)
- **재발 시**: 베바시주맙 ± 이리노테칸, 펨브롤리주맙, 로무스틴
- **신규**: TTF (Optune®) → EF-14 연구에서 중앙 OS 20.9개월

---

*모델 생성: Claude Code Routine (CCR), 2026-06-26*
*참고 가이드라인: EANO 2021, NCCN CNS Tumors v2.2024*
