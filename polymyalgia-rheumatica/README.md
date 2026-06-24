# Polymyalgia Rheumatica (PMR) — QSP Model

> **류마티카 다발성 근통** | 자가면역/염증 질환 | IL-6 경로 · 코르티코스테로이드 PK/PD · 토실리주맙

[![Graphviz](https://img.shields.io/badge/Map-130+%20nodes%2C%2012%20clusters-blue)](pmr_qsp_model.svg)
[![mrgsolve](https://img.shields.io/badge/ODE-22%20compartments%2C%207%20scenarios-green)](pmr_mrgsolve_model.R)
[![Shiny](https://img.shields.io/badge/Shiny-6%20tabs%20dashboard-orange)](pmr_shiny_app.R)
[![References](https://img.shields.io/badge/Refs-55%20PubMed%20citations-red)](pmr_references.md)

---

## 질환 개요 (Disease Overview)

류마티카 다발성 근통(PMR)은 **50세 이상**에서 발생하는 흔한 염증성 류마티스 질환으로, 어깨와 골반 거들의 **양측성 근통과 조조강직**을 특징으로 합니다.

| 특성 | 내용 |
|------|------|
| **발병률** | 50–100/100,000/년 (50세 이상) |
| **성별** | 여:남 ≈ 2–3:1 |
| **발병 연령** | 주로 70–80대 |
| **지역 편향** | 북유럽 계통에서 높음 |
| **GCA 중복** | PMR의 약 15–20%에서 거대세포 동맥염 동반 |
| **치료 반응** | 코르티코스테로이드(prednisolone)에 극적 반응 — **진단적 특징** |

### 핵심 증상
- 양측 어깨·골반 거들 통증 및 무력감
- 조조강직 ≥45분
- CRP/ESR 상승 (급성기 반응물)
- 체중 감소, 발열, 피로감

---

## 기계론적 지도 (Mechanistic Map)

[![PMR QSP Map](pmr_qsp_model.png)](pmr_qsp_model.svg)

> *이미지 클릭 시 SVG 전체 화면으로 확인 가능*

### 12개 클러스터 구성

| 클러스터 | 핵심 내용 | 주요 노드 |
|---------|----------|----------|
| **역학/위험인자** | HLA-DRB1*04, PTPN22, 환경 유발 | HLA_DRB1_04, PTPN22_variant |
| **선천 면역** | DC, 대식세포 M1, 중성구, NLRP3 | NLRP3, NF_kB_innate, MacM1 |
| **적응 면역** | Th1/Th17/Treg 균형 이상 | Th17_cells, Treg_cells, RORgt |
| **사이토카인** | IL-6 중심, TNF-α, IL-17A | IL6, sIL6R, IL17A, TNF_alpha |
| **JAK-STAT 신호** | JAK1/2-STAT3, SOCS3 피드백 | STAT3, JAK1, SOCS3 |
| **표적 조직** | 어깨·고관절 활막/점액낭 | Subacromial_bursa, FLS_synov |
| **Prednisolone PK/PD** | GR 결합, 전위억제, GILZ | GR_Pred_complex, Transrepression |
| **Tocilizumab PK/PD** | mIL-6R/sIL-6R 차단, TMDD | TCZ_mIL6R_cpx, IL6_signal_blk |
| **HPA 축** | 코르티솔 억제, 부신 위축 | ACTH, Cortisol_endog, HPA_feedback |
| **골 효과** | RANK/RANKL/OPG, BMD 감소 | Osteoclast, RANKL, BMD_lumbar |
| **혈관/GCA** | 측두동맥, 대동맥 염증 | Temporal_artery, PMR_GCA_overlap |
| **임상 엔드포인트** | PMR-AS, CRP, ESR | PMR_AS, Remission, Relapse_event |

---

## 발병기전 요약 (Key Pathophysiology)

```
환경 유발(감염 등) + 유전 소인(HLA-DRB1*04)
        ↓
  선천 면역 활성화 (mDC, NLRP3, TLR4)
        ↓
  IL-23 → Th17 분화 / IFN-γ → Th1 편향
     ↓                    ↓
  IL-17A/F             IFN-γ, TNF-α
     ↓                    ↓
      ──→ IL-6 폭풍 (근육·관절 조직)
              ↓
    JAK1/2 → STAT3 → 급성기 단백(CRP, 피브리노겐)
              ↓
      통증·조조강직 (PGE2, 조직 부종)
              ↓
          PMR-AS ↑

치료:
  Prednisolone → GR 활성화
    → NF-κB·AP-1 전위억제 → IL-6·TNF-α↓
    → GILZ·Annexin A1↑ → 항염증
  Tocilizumab (IL-6R 차단) → STAT3 억제 → CRP 정상화
```

---

## mrgsolve 모델 사양 (ODE Model Specifications)

**파일**: [`pmr_mrgsolve_model.R`](pmr_mrgsolve_model.R)

### 22개 구획 (Compartments)

| 구획 그룹 | 구획 수 | 내용 |
|----------|--------|------|
| Prednisolone PK | 3 | Depot → Central ↔ Peripheral |
| Tocilizumab PK | 3 | Depot → Central ↔ Peripheral (TMDD) |
| HPA 축 | 1 | 코르티솔 (억제 반응) |
| IL-6 경로 | 2 | IL-6, 가용성 IL-6Rα |
| 급성기 반응물 | 2 | CRP, ESR |
| 골 효과 | 1 | BMD (정규화) |
| 질환 활성도 | 2 | PMR-AS, 재발 위험 |

### 7개 치료 시나리오

| # | 시나리오 | 약물·용량 | 근거 임상시험 |
|---|---------|----------|------------|
| **S1** | 무치료 (자연경과) | — | 자연경과 코호트 |
| **S2** | Pred 15mg → 테이퍼 2.5mg/mo | ACR 표준 | Dejaco 2015 ACR/EULAR |
| **S3** | Pred 22.5mg → 급속 테이퍼 4mg/mo | 중증례 | BSR 가이드라인 |
| **S4** | Pred 15mg → 완만 테이퍼 1mg/mo | 재발 방지 | 관찰 코호트 |
| **S5** | TCZ 162mg SC QW + Pred 12.5mg | GC 절약 | GiACTA (Stone 2017 NEJM) |
| **S6** | TCZ 162mg SC Q2W + Pred 12.5mg | GC 절약 | SEMAPHORE/SAPHYR 참고 |
| **S7** | TCZ QW만 (스테로이드 무병 유도) | 탐색적 | PMR-SPARE Phase 2 |

### 주요 파라미터 (임상시험 보정)

| 파라미터 | 값 | 출처 |
|---------|-----|------|
| Pred CL | 14.0 L/h | Bergmann 2012 |
| Pred V1 | 30.0 L | Bergmann 2012 |
| Pred Fu | 0.28 | Buttgereit 2005 |
| TCZ CL (선형) | 0.29 L/h | Nishimoto 2008 |
| IL-6 기저 (PMR) | 15 pg/mL | Roche 1993 |
| CRP 기저 (PMR) | 35 mg/L | 임상 관찰 |
| EC50 (Pred→PMR-AS) | 180 ng/mL free | 모델 보정 |
| EC50 (TCZ→PMR-AS) | 50 nM | 모델 보정 |

---

## Shiny 대시보드 (Interactive Dashboard)

**파일**: [`pmr_shiny_app.R`](pmr_shiny_app.R)

### 6개 탭 구성

| 탭 | 내용 |
|----|------|
| **① 환자 프로파일** | 진단 기준, 기저 파라미터, PMR-AS/CRP/ESR/IL-6 value box |
| **② 약물 PK** | Prednisolone 총/유리 농도, Tocilizumab nM 추이, PK 요약 표 |
| **③ 염증 마커** | IL-6, sIL-6Rα, CRP, ESR 시계열 플롯 (정상 기준선 포함) |
| **④ 질환 활성도** | PMR-AS 추이, HPA 축 코르티솔 억제, 재발 위험 점수 |
| **⑤ 시나리오 비교** | 7개 치료군 CRP·PMR-AS·IL-6·BMD 비교 (색상 구분) |
| **⑥ 바이오마커 탐색기** | BMD 추이, GC 절약 분석, IL-6↔CRP 상관, 참조값 표 |

### 실행 방법
```r
install.packages(c("shiny","shinydashboard","mrgsolve","dplyr","ggplot2","plotly","DT"))
shiny::runApp("pmr_shiny_app.R")
```

---

## 참고문헌 (References)

**파일**: [`pmr_references.md`](pmr_references.md)

| 섹션 | 인용 수 |
|------|--------|
| 역학 및 임상 특징 | 5 |
| 진단 기준 | 4 |
| 치료 가이드라인 | 4 |
| 발병기전 및 면역학 | 7 |
| IL-6 경로 및 사이토카인 | 5 |
| 코르티코스테로이드 PK/PD | 5 |
| 토실리주맙 PK/PD 및 임상시험 | 6 |
| HPA 축 억제 | 3 |
| 골다공증 및 골 효과 | 4 |
| 질환 활성도 점수 | 4 |
| 혈관 침범 및 GCA 중복 | 4 |
| QSP/PK-PD 모델링 방법론 | 4 |
| **합계** | **55** |

---

## 산출물 요약 (Deliverables)

| 구성요소 | 파일 | 사양 |
|---------|------|------|
| 🗺️ 기계론적 지도 | [`pmr_qsp_model.dot`](pmr_qsp_model.dot) · [`.svg`](pmr_qsp_model.svg) · [`.png`](pmr_qsp_model.png) | **130+ 노드, 12 클러스터** |
| ⚙️ mrgsolve ODE | [`pmr_mrgsolve_model.R`](pmr_mrgsolve_model.R) | **22구획 ODE, 7치료 시나리오, VPop 200명** |
| 📊 Shiny 앱 | [`pmr_shiny_app.R`](pmr_shiny_app.R) | **6탭** (환자프로파일·PK·염증마커·질환활성도·시나리오비교·바이오마커탐색기) |
| 📚 참고문헌 | [`pmr_references.md`](pmr_references.md) | **55개 PubMed 인용** (12개 섹션) |

---

## 임상적 시사점 (Clinical Implications)

| 주제 | 시사점 |
|------|--------|
| **조기 진단** | CRP >5mg/L + ESR >50mm/hr + 어깨 양측성 통증 + 조조강직 ≥45min → PMR 강력 의심 |
| **GC 시작** | Pred 12.5–25mg/d → 빠른 반응(24–72h) 없으면 진단 재고 |
| **테이퍼링** | ACR 권고: 4주 후 매월 2.5mg 감량 (재발 시 이전 용량 복귀) |
| **TCZ 추가** | 재발성/GC 의존성 PMR에서 TCZ 162mg QW → 관해 유지 + GC 절약 |
| **골다공증 예방** | GC ≥3개월 사용 시 bisphosphonate + Vit D/Ca 필수 |
| **GCA 감시** | 두통·측두동맥 압통·시력 이상 → 즉각 GCA 평가 (고용량 GC 긴급 투여) |

---

*생성일: 2026-06-24 | Claude Code Routine (CCR) 자동 생성*
