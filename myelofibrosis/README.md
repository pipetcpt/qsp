# 골수섬유증 (Myelofibrosis · MF) QSP 모델

> **디렉토리**: `myelofibrosis/` | **날짜**: 2026-06-23 | **카테고리**: 종양·혈액

[![MF QSP Map](mf_qsp_model.png)](mf_qsp_model.svg)

---

## 질환 개요 (Disease Overview)

**골수섬유증(Myelofibrosis, MF)** 은 조혈 줄기세포(HSC)에서 유래하는 만성 골수증식성 종양(MPN)으로, 골수 내 콜라겐·레티큘린 섬유화, 무효조혈(ineffective hematopoiesis), 비정상적 거핵구 증식, 그리고 비장·간을 중심으로 하는 **수외 조혈(extramedullary hematopoiesis, EMH)** 을 특징으로 합니다. 일차성(Primary MF, PMF)과 다혈증성·본태성 혈소판증가증 후 이행성(Post-PV/ET MF)으로 구분됩니다.

- **발생률**: 미국 기준 연 약 3,000~4,000명, 중앙 발병 연령 67세
- **5년 생존율**: 중간위험-2/고위험군 ~35%, 저위험군 ~70% (DIPSS Plus 기준)
- **AML 전환율**: 연 5~10%, 고위험 분자 이상(ASXL1, EZH2, SRSF2, IDH1/2, TP53) 보유 시 현저히 증가

### 분자 병인 (Molecular Pathogenesis)

| 돌연변이 | 빈도 | 기전 |
|---------|-----|------|
| **JAK2 V617F** | ~55–60% | JAK2 구성적 활성화 → STAT3/5 과인산화 → 무제한 증식 |
| **CALR exon 9 (Type 1/2)** | ~25–30% | MPL 수용체 이상 활성화 → JAK/STAT 신호 항진 |
| **MPL W515L/K** | ~5–10% | TPO 수용체 돌연변이 → JAK/STAT 신호 항진 |
| **Triple-negative** | ~10% | 기타 미확인 기전 |
| **HMR 돌연변이** | 다양 | ASXL1(30–40%), EZH2(5%), SRSF2(10%), IDH1/2(5%) — 고위험 |

---

## QSP 모델 구성 (Model Architecture)

### 4종 산출물 (Four Deliverables)

| 산출물 | 파일 | 세부 사양 |
|--------|------|----------|
| 🗺️ **기계론적 지도** | [`mf_qsp_model.dot`](mf_qsp_model.dot) / [`.svg`](mf_qsp_model.svg) / [`.png`](mf_qsp_model.png) | **204 노드, 12 클러스터, 224 엣지** |
| ⚙️ **mrgsolve ODE 모델** | [`mf_mrgsolve_model.R`](mf_mrgsolve_model.R) | **23 구획 ODE**, 6 치료 시나리오 |
| 📊 **Shiny 대시보드** | [`mf_shiny_app.R`](mf_shiny_app.R) | **6 탭** 인터랙티브 앱 |
| 📚 **참고문헌** | [`mf_references.md`](mf_references.md) | **36개** PubMed 인용 |

---

### 기계론적 지도 — 12개 클러스터 (Mechanistic Map)

| # | 클러스터 | 주요 구성요소 |
|---|---------|-------------|
| 1 | **유전·돌연변이 드라이버** | JAK2V617F, CALR T1/T2, MPL W515L/K, ASXL1, EZH2, SRSF2, IDH1/2, TP53, TET2, DNMT3A, DIPSS/MIPSS70 위험도 점수 |
| 2 | **JAK/STAT 신호** | JAK1/2/TYK2, STAT1/3/5a/5b, pSTAT3/pSTAT5, SOCS1/3, SHP2, BCL2/MCL1/CCND1/MYC/PIM1 |
| 3 | **HSC 니치** | LT-HSC → ST-HSC → MPP → CMP/GMP/MEP → MkP/BFU-E; CXCL12/CXCR4, SCF/c-KIT, TPO, EPO |
| 4 | **골수 미세환경·섬유화** | MSC, 조골/파골세포, 비정상 거핵구, TGF-β1/PDGF/bFGF/CTGF, 콜라겐 I/III, MF grade 0–3, 골경화 |
| 5 | **사이토카인 폭풍** | IL-1β/6/8/10/12/13, TNF-α, IFN-γ, CXCL10, CCL2, NF-κB, mTOR, PI3K, RAS/ERK, NLRP3 |
| 6 | **수외 조혈** | 순환 CD34+, 비장/간 EMH, 비장·간 부피, 골수 부전 |
| 7 | **혈액학적 결과** | Hgb 역학, RBC 생산/파괴, 혈소판, WBC/중성구, 순환 아세포, 수혈 의존성 |
| 8 | **혈전·혈관** | 혈소판 활성화, TXA2, 트롬빈 생성, TF/FXa 연쇄, PAI-1, 내피 기능 이상, DVT/PE/BCS/PVT |
| 9 | **약물 PK** | Ruxolitinib 2구획(ka/CL/V1/V2/Q/F1), Fedratinib, Pacritinib, Momelotinib PK |
| 10 | **약물 PD** | JAK1/2 억제(Emax/IC50), pSTAT3/5 억제, SVR35/TSS50 출력, 빈혈/혈소판감소증 부작용 |
| 11 | **임상 엔드포인트** | SVR35, TSS50, 골수 조직 반응, CHR, OS, PFS, AML 전환, IWG-MRT/ELN 2023 기준 |
| L | **범례** | 노드 형태·색상 안내 |

---

### mrgsolve ODE 모델 — 23개 구획

```
Ruxolitinib PK  : DEPOT_RUX, CENT_RUX, PERI_RUX             (3구획)
Fedratinib PK   : DEPOT_FED, CENT_FED                        (2구획)
Pacritinib PK   : DEPOT_PAC, CENT_PAC                        (2구획)
BET 억제제 PK   : CENT_BET  (pelabresib)                     (1구획)
JAK/STAT PD     : pSTAT3, pSTAT5                             (2구획)
클론 역학        : NHSC (종양 HSC), NHSC_N (정상 HSC)         (2구획)
적혈계           : PROG_E, RET, RBC                           (3구획)
거핵계           : MEG_P, PLT                                 (2구획)
거시적 지표      : SPLEEN, FIBROSIS                           (2구획)
사이토카인       : IL6, TNF                                   (2구획)
증상 점수        : TSS                                        (1구획)
```

**임상 보정 참고 임상시험**

| 시나리오 | 임상시험 | 주요 결과 |
|---------|---------|---------|
| Ruxolitinib 20 mg BID | COMFORT-I (Verstovsek 2012) | SVR35 41.9%, TSS50 45.9% |
| Ruxolitinib 15 mg BID | COMFORT-I 저혈소판 코호트 | SVR35 ~28% |
| Fedratinib 400 mg QD | JAKARTA (Pardanani 2015) | SVR35 36%, TSS50 36% |
| Pacritinib 200 mg BID | PERSIST-2 (Mesa 2017) | SVR35 18%, 혈소판 <50×10⁹/L |
| Ruxolitinib + Pelabresib | MANIFEST-2 (Pemmaraju 2024) | SVR35 66% vs 35% |
| 무치료 | 자연 경과 | 비장 연 ~10% 증가 |

---

### Shiny 대시보드 — 6개 탭

| 탭 | 내용 |
|----|------|
| **① 환자 프로파일** | 나이/성별/진단, DIPSS Plus 위험도, 돌연변이(JAK2/CALR/MPL), 기저 비장 부피·Hgb·혈소판, 증상 점수 |
| **② 약동학 (PK)** | 약물 선택, 용량·투여 간격 입력, 혈중 농도-시간 곡선(Cp vs time), Cmax/Cmin/AUC/t½ 테이블 |
| **③ PD 바이오마커** | pSTAT3/pSTAT5 억제 프로파일, JAK2 V617F VAF 추이, IL-6/TNF-α 사이토카인 동태 |
| **④ 임상 엔드포인트** | 비장 부피 (SVR35 기준선), Hgb 변화, 혈소판 수, TSS 점수, 골수 섬유화 등급 |
| **⑤ 치료 비교** | 다중 약물 체크박스, SVR35/TSS50/Hgb 나비 플롯, waterfall 비장 반응 |
| **⑥ 바이오마커 역학** | JAK2 VAF–비장 반응 상관관계, 사이토카인 히트맵, AML 전환 위험 KM 곡선 |

---

## 실행 방법 (How to Run)

### 1. 기계론적 지도 렌더링

```bash
# Graphviz 설치 필요
dot -Tsvg myelofibrosis/mf_qsp_model.dot -o myelofibrosis/mf_qsp_model.svg
dot -Tpng -Gdpi=150 myelofibrosis/mf_qsp_model.dot -o myelofibrosis/mf_qsp_model.png
```

### 2. mrgsolve 모델 실행

```r
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
library(mrgsolve)
source("myelofibrosis/mf_mrgsolve_model.R")
# 시뮬레이션 및 플롯이 자동 실행됩니다
```

### 3. Shiny 앱 실행

```r
install.packages(c("shiny", "shinydashboard", "ggplot2", "plotly", "DT"))
shiny::runApp("myelofibrosis/mf_shiny_app.R")
```

---

## 주요 임상시험 요약 (Key Clinical Trials)

| 시험명 | 약물 | 대상 | SVR35 | TSS50 | 중앙 OS | PMID |
|-------|------|------|-------|-------|---------|------|
| **COMFORT-I** | Ruxolitinib 20mg BID vs 위약 | Int-2/High PMF | 41.9% vs 0.7% | 45.9% vs 5.3% | — | [22375971](https://pubmed.ncbi.nlm.nih.gov/22375971/) |
| **COMFORT-II** | Ruxolitinib vs BAT | Int-2/High PMF | 28% vs 0% | — | NR vs 27.9mo | [22375970](https://pubmed.ncbi.nlm.nih.gov/22375970/) |
| **JAKARTA** | Fedratinib 400mg QD | Int-2/High, Rux-naive | 36% vs 1% | 36% vs 6% | — | [26003172](https://pubmed.ncbi.nlm.nih.gov/26003172/) |
| **PERSIST-2** | Pacritinib 200mg BID | PLT <100×10⁹/L | 18% vs 3% | — | — | [29049469](https://pubmed.ncbi.nlm.nih.gov/29049469/) |
| **SIMPLIFY-1** | Momelotinib 200mg QD | Rux-naive | 26.5% vs 29% | 28.4% vs 42.2% | — | [28930484](https://pubmed.ncbi.nlm.nih.gov/28930484/) |
| **MANIFEST-2** | Pelabresib + Ruxolitinib | Rux-naive | 65.9% vs 35.2% | 52.3% vs 37.5% | NR | [39504566](https://pubmed.ncbi.nlm.nih.gov/39504566/) |

---

## 참고문헌 (References)

전체 참고문헌 목록은 [`mf_references.md`](mf_references.md)를 참조하세요 (36개 PubMed 인용).

주요 인용:
- Verstovsek S et al. *N Engl J Med* 2012;366:799–807 (COMFORT-I) · PMID [22375971](https://pubmed.ncbi.nlm.nih.gov/22375971/)
- Pardanani A et al. *J Clin Oncol* 2015;33:2771–2779 (JAKARTA) · PMID [26003172](https://pubmed.ncbi.nlm.nih.gov/26003172/)
- Klampfl T et al. *N Engl J Med* 2013;369:2379–2390 (CALR mutations) · PMID [24325356](https://pubmed.ncbi.nlm.nih.gov/24325356/)
- James C et al. *Nature* 2005;434:1144–1148 (JAK2 V617F discovery) · PMID [15793561](https://pubmed.ncbi.nlm.nih.gov/15793561/)

---

*이 모델은 Claude Code Routine(CCR)에 의해 2026-06-23에 자동 생성되었습니다.*  
*교육·연구 목적으로 제작되었으며, 임상 의사결정에 직접 사용해서는 안 됩니다.*
