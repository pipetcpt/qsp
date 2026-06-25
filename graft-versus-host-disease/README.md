# Graft-versus-Host Disease (GvHD) QSP Model

[![Disease](https://img.shields.io/badge/Disease-GvHD-red)](.) [![Category](https://img.shields.io/badge/Category-Transplant%20Immunology-blue)](.) [![Drugs](https://img.shields.io/badge/Drugs-CsA%20%7C%20TAC%20%7C%20RUX%20%7C%20BELU%20%7C%20MMF-orange)](.)

## 질환 개요 (Disease Overview)

**이식편대숙주병(GvHD, Graft-versus-Host Disease)** 은 동종 조혈모세포이식(Allo-HSCT) 후 공여자 면역세포가 수혜자의 조직을 공격하는 중증 합병증입니다. 급성 GvHD(aGvHD)는 이식 후 100일 이내, 만성 GvHD(cGvHD)는 그 이후에 주로 발생합니다.

**Graft-versus-Host Disease (GvHD)** is a severe complication following allogeneic hematopoietic stem cell transplantation (Allo-HSCT), where donor immune cells attack recipient tissues. Acute GvHD typically occurs within 100 days post-HSCT, while chronic GvHD develops later and can persist for years.

---

## 기계론적 지도 (Mechanistic Map)

[![GvHD Mechanistic Map](gvhd_qsp_model.png)](gvhd_qsp_model.svg)

*Click the image to open the full SVG with interactive zoom*

### 주요 클러스터 (Key Clusters)

| # | Cluster | 핵심 내용 |
|---|---------|---------|
| 1 | **HSCT Context & Conditioning** | MAC/RIC 전처치, TBI, 조직손상, DAMP/PAMP 방출, 장벽 손상 |
| 2 | **Antigen Presentation & Priming** | Host DC, MHC mismatch, 직접/간접 동종반응, CD28-B7 공자극 |
| 3 | **Donor T Cell Differentiation** | Th1/Th17/Treg/CD8 분화, NFAT·NF-κB·JAK-STAT·ROCK2 신호전달 |
| 4 | **Cytokine Network** | TNF-α, IFN-γ, IL-6, IL-17A, IL-10, TGF-β, BAFF |
| 5 | **Target Organ: Skin** | 모양세포 소양증, lichenoid/sclerotic 병변, mLSS 점수 |
| 6 | **Target Organ: Gut** | 장상피세포 아포토시스, 크립트 손상, ST2/REG3α 바이오마커 |
| 7 | **Target Organ: Liver** | 담관 손상, 담즙정체, Glucksberg 등급 |
| 8 | **Target Organ: Lung** | 폐쇄성 세기관지염(BOS), FEV1 감소, CLAD 점수 |
| 9 | **B Cell Pathology** | Tfh-B GC 반응, 자가항체, BTK 경로 |
| 10 | **Fibrosis (TGF-β/ROCK2)** | SMAD2/3, EMT, 근섬유아세포 활성화, ROCK2/IRF4 |
| 11 | **Drug PK: CNI** | CsA 2구획, TAC 2구획, CYP3A4/5 대사 |
| 12 | **Drug PK/PD: Ruxolitinib** | JAK1/2 저해, STAT3/5 차단, Treg 확장 |
| 13 | **Other Drugs** | 스테로이드(NF-κB), Belumosudil(ROCK2), MMF(IMPDH), Ibrutinib(BTK) |
| 14 | **Clinical Endpoints & Biomarkers** | Glucksberg/NIH 점수, ORR, FFS, OS, NRM |

**총 노드 수: ~130+ | 서브그래프 클러스터: 14**

---

## mrgsolve ODE 모델 (Pharmacokinetic/Pharmacodynamic Model)

**파일**: `gvhd_mrgsolve_model.R`

### 구획 구성 (Compartments: 32 total)

**약동학 (PK) — 16 compartments:**
| 약물 | 구획 | 특징 |
|------|------|------|
| Cyclosporine A (CsA) | 3 (Gut/Central/Peripheral) | F=30%, CYP3A4, C₀ target 100-300 ng/mL |
| Tacrolimus (TAC) | 3 (Gut/Central/Peripheral) | F=25%, CYP3A5, C₀ target 5-15 ng/mL |
| Prednisone (PRED) | 2 (Gut/Central) | F=99%, GRα 결합 |
| Ruxolitinib (RUX) | 3 (Gut/Central/Peripheral) | F=95%, T½~3h, V=72L |
| Belumosudil (BELU) | 2 (Gut/Central) | F=80%, T½~20h, ROCK2 선택적 저해 |
| MMF/MPA | 2 (Gut/Central) | F=94%, IMPDH 저해 |

**면역/생물학 (PD) — 16 compartments:**
- T 세포: Th1, Th17, Treg, CD8 effector
- B 세포 풀
- 사이토카인 6종: TNF-α, IFN-γ, IL-17A, IL-10, TGF-β, IL-6
- 장기 손상 5종: Skin, Gut, Liver, Lung, Fibrosis

### 약물 효과 모델 (PD Effects)

| 약물 | 표적 | 모델 |
|------|------|------|
| CsA/TAC | Calcineurin → NFAT → IL-2 | Emax (Hill n=1.5) |
| Ruxolitinib | JAK1/2 → STAT3/5 | Emax, Treg 확장 효과 포함 |
| Belumosudil | ROCK2 → IRF4/STAT3 → Th17↓/Treg↑ | Emax, 섬유화 억제 |
| MMF/MPA | IMPDH → 림프구 증식 억제 | Emax |
| Prednisone | GRα → NF-κB 억제 → broad cytokine↓ | Emax |

### 치료 시나리오 (6 Treatment Scenarios)

1. **No Prophylaxis** — 기저 상태 (historical control)
2. **CsA Monoprophylaxis** — 단독 CsA 예방요법
3. **CsA + MMF** — 표준 조합 예방요법
4. **TAC + MMF** — 현재 가장 흔한 표준요법 (NMDP/EBMT 권고)
5. **CsA → Ruxolitinib** — 스테로이드 불응성 cGvHD (REACH3 기반)
6. **CsA → Belumosudil** — 2선 이상 cGvHD (ROCKstar 기반)

---

## Shiny 대시보드 (Interactive Dashboard)

**파일**: `gvhd_shiny_app.R`

### 탭 구성 (8 Tabs)

| 탭 | 내용 |
|----|------|
| 1. **Patient & HSCT Profile** | 환자 프로파일, GvHD 위험도 레이더 차트, 약물 표적 설명 |
| 2. **Drug PK Dashboard** | CsA/TAC/RUX/BELU 농도-시간 곡선, PK 요약 테이블, PD 효과 |
| 3. **Immune Cell Dynamics** | Th1/Th17/Treg/CD8 동태, Th17/Treg 비율, B세포 |
| 4. **Cytokine Network** | 전염증성(TNF-α, IFN-γ, IL-17A, IL-6) vs 항염증성(IL-10, TGF-β) |
| 5. **Organ Damage & Endpoints** | 피부/장/간/폐 손상 점수, aGvHD/cGvHD 등급, FFS |
| 6. **Scenario Comparison** | 6개 치료 시나리오 병렬 비교 |
| 7. **Biomarkers** | ST2, REG3α, sTNFR1 바이오마커 시계열 |
| 8. **Mechanistic Map** | 전체 기전 지도 PNG/SVG |

---

## 임상적 맥락 (Clinical Context)

### 역학 (Epidemiology)
- 동종 HSCT 후 **aGvHD 발생률**: 형제 공여 30-50%, 비혈연 공여 50-70%
- **cGvHD**: aGvHD 생존자의 40-70%에서 발생
- GvHD는 비재발 사망(NRM)의 가장 주요한 원인

### 발병 기전 3단계 (Three-Phase Pathophysiology)
1. **Phase 1 (Afferent)**: 전처치에 의한 조직 손상 → DAMP/PAMP 방출 → 숙주 DC 활성화
2. **Phase 2 (Efferent)**: 공여자 T세포의 동종항원 인식 → Th1/Th17 극화, CD8 활성화
3. **Phase 3 (Effector)**: TNF-α/IFN-γ/IL-17A 매개 표적 장기 손상

### FDA 승인 약물
| 약물 | 적응증 | 근거 임상시험 |
|------|--------|------------|
| **Ruxolitinib** (Jakafi) | 스테로이드 불응성 aGvHD, cGvHD | REACH1, REACH2, REACH3 |
| **Belumosudil** (Rezurock) | 2선 이상 cGvHD | ROCKstar |
| **Ibrutinib** (Imbruvica) | 1선 이후 cGvHD | 단일군 연구 |

---

## 주요 바이오마커 (Key Biomarkers)

| 바이오마커 | 정상치 | GvHD 예측 | 표적 장기 |
|-----------|--------|-----------|---------|
| ST2 (sST2) | <33 ng/mL | >33 ng/mL → grade 3-4 GI GvHD | 장 |
| REG3α | <10 ng/mL | >23 ng/mL → 불량 예후 | 장 상피 |
| sTNFR1 | <2 ng/mL | 상승 → aGvHD 중증도 | 전신/피부 |
| CXCL9 | <100 pg/mL | IFN-γ 유도; T세포 귀소 | 전신 |
| Elafin | <10 ng/mL | 피부 GvHD 마커 | 피부 |

---

## 실행 방법 (How to Run)

```bash
# 1. Graphviz 렌더링
dot -Tsvg gvhd_qsp_model.dot -o gvhd_qsp_model.svg
dot -Tpng -Gdpi=150 gvhd_qsp_model.dot -o gvhd_qsp_model.png
```

```r
# 2. mrgsolve 모델 실행 (R)
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
source("gvhd_mrgsolve_model.R")

# 3. Shiny 대시보드
install.packages(c("shiny", "shinydashboard", "plotly", "DT"))
shiny::runApp("gvhd_shiny_app.R")
```

---

## 파일 목록 (File List)

| 파일 | 설명 |
|------|------|
| `gvhd_qsp_model.dot` | Graphviz 기계론적 지도 소스 (130+ 노드, 14 클러스터) |
| `gvhd_qsp_model.svg` | 벡터 형식 지도 (브라우저에서 대화형 확대 가능) |
| `gvhd_qsp_model.png` | 래스터 형식 지도 (150 dpi) |
| `gvhd_mrgsolve_model.R` | mrgsolve ODE 모델 + 6개 치료 시나리오 + 시각화 |
| `gvhd_shiny_app.R` | Shiny 대화형 시뮬레이터 (8 탭) |
| `gvhd_references.md` | 참고문헌 60편 (섹션별 분류) |
| `README.md` | 이 파일 |

---

## 참고문헌 주요 인용 (Key References)

- Zeiser R et al. *NEJM* 2020;382:1800 — Ruxolitinib for SR aGvHD (REACH2)
- Zeiser R et al. *NEJM* 2021;385:228 — Ruxolitinib for SR cGvHD (REACH3)
- Cutler C et al. *Blood* 2021;138:2278 — Belumosudil ROCKstar trial
- Ferrara JL et al. *Lancet* 2009;373:1550 — GvHD pathophysiology review
- Vander Lugt MT et al. *NEJM* 2013;369:529 — ST2 biomarker
- Peled JU et al. *NEJM* 2020;382:822 — Microbiome & GvHD outcome

전체 참고문헌 60편 → [`gvhd_references.md`](gvhd_references.md)

---

*Generated by Claude Code Routine (CCR) — 2026-06-25*  
*Disease category: Transplant Immunology / HSCT Complication*
