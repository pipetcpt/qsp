# 패혈증 (Sepsis & Septic Shock) — QSP Model

> **Directory:** `sepsis/` | **Abbreviation:** SEP | **Date:** 2026-06-24  
> **Category:** Acute Critical Care / Infectious Disease  
> **Global Burden:** ~49 million cases/year; ~11 million deaths (22% of all global deaths)

[![SEP QSP Mechanistic Map](sep_qsp_model.png)](sep_qsp_model.svg)

---

## 질환 개요 (Disease Overview)

**패혈증(Sepsis)**은 감염에 대한 숙주의 조절 장애 반응(dysregulated host response)으로 생명을 위협하는 장기부전(organ dysfunction)이 발생하는 질환입니다 (Sepsis-3, Singer et al. 2016). 단순한 전신염증반응(SIRS)을 넘어, 면역계의 동시적 과활성화(pro-inflammatory storm)와 면역억제(immunosuppression/CARS)가 병존하는 복잡한 병태생리가 특징입니다.

**패혈성 쇼크(Septic Shock)**는 패혈증의 하위 집단으로, 혈압저하(MAP <65 mmHg)와 함께 혈중 젖산>2 mmol/L가 지속되어 28일 사망률이 40% 이상에 이르는 중증 상태입니다.

---

## 병태생리 요약 (Pathophysiology Summary)

| 단계 | 핵심 기전 | 주요 노드 |
|------|----------|----------|
| **1. 병원체 인식** | LPS/PGN/DAMP → TLR4/TLR2/NLRP3 → NF-κB 활성화 | PAMPs, DAMPs, PRRs, MyD88, NFkB |
| **2. 사이토카인 폭풍** | TNFα·IL-1β·IL-6·IL-8 과분비 → 이차 세포 활성화 | TNF, IL6, IL1B, IL10, IFNγ |
| **3. 선천면역 과활성** | 호중구 조직 침윤·NET·ROS·MMP-9 | Neut_T, NET, ROS, Elastase |
| **4. 보체 활성화** | C3→C5a(아나필라톡신)→백혈구·혈관 반응 | C5a, C5aR, MAC |
| **5. 응고/DIC** | TF 발현↑→트롬빈↑→피브린↑+PAI-1↑→DIC | Thrombin, Fibrin, PAI1 |
| **6. 내피세포 장애** | VE-cadherin 분리·혈관 투과성↑·NO↑→혈관확장쇼크 | ENDOT, VascPerm, MAP_vasc |
| **7. 다장기부전(MODS)** | 폐(ARDS)·신장(AKI)·간·뇌·순환부전 | SOFA = 0–24 |
| **8. 후기 면역억제** | T세포 아포프토시스·PD-1↑·MDSC↑→CARS | Immunosuppression, CARS |

---

## 치료 시나리오 (Treatment Scenarios)

| 시나리오 | 치료 | 기전 | 주요 임상시험 |
|---------|------|------|-------------|
| **S1** | 무치료 (자연경과) | — | 관찰 코호트 |
| **S2** | 항생제 단독 (Meropenem 1g q8h) | fT>MIC 달성 → 균 사멸 | Craig 1998 CID |
| **S3** | 항생제 + 노르에피네프린 | α1-adrenoceptor→MAP 회복 | SOAP II NEJM 2010 |
| **S4** | 번들 (항생제+NE+수액 30 mL/kg) | RAAS 보정·전부하 회복 | EGDT NEJM 2001 |
| **S5** | 번들 + 하이드로코티손 200 mg/day | GR→사이토카인 억제·혈관수축제 민감도↑ | ADRENAL/APROCCHSS NEJM 2018 |
| **S6** | 번들 + HC + 토실리주맙 8 mg/kg | IL-6R 차단→STAT3 억제→CRP↓·PCT↓ | REMAP-CAP NEJM 2021 |
| **S7** | 면역저하 환자 (항생제+NE, 고균량) | 비정상적 면역반응, 높은 이환·사망 | 임상 코호트 |

---

## QSP 모델 사양 (Model Specifications)

| 구성요소 | 파일 | 사양 |
|---------|------|-----|
| 🗺️ 기계론적 지도 | [`sep_qsp_model.dot`](sep_qsp_model.dot) / [`.svg`](sep_qsp_model.svg) / [`.png`](sep_qsp_model.png) | **130+ 노드, 11 클러스터** |
| ⚙️ mrgsolve ODE 모델 | [`sep_mrgsolve_model.R`](sep_mrgsolve_model.R) | **24구획 ODE**, **7치료 시나리오** |
| 📊 Shiny 대시보드 | [`sep_shiny_app.R`](sep_shiny_app.R) | **8탭** (환자프로파일·항생제PK·사이토카인/면역·혈역학/SOFA·장기기능·시나리오비교·바이오마커·About) |
| 📚 참고문헌 | [`sep_references.md`](sep_references.md) | **55개 PubMed 인용** (14개 섹션) |

---

## 기계론적 지도 클러스터 (Mechanistic Map Clusters)

| # | 클러스터 | 노드 수 | 핵심 내용 |
|---|---------|---------|---------|
| 1 | 병원체 인식 (Pathogen Recognition) | 20 | LPS, PGN, TLR4/2/5/9, NLRP3, NF-κB |
| 2 | 선천면역 활성화 (Innate Immunity) | 15 | Neutrophils, NETs, Macrophages, ROS |
| 3 | 사이토카인 네트워크 (Cytokine Storm) | 18 | TNFα, IL-1β, IL-6, IL-8, IL-10, HMGB1 |
| 4 | 보체 계통 (Complement) | 12 | C3, C5a, MAC, C5aR |
| 5 | 응고/DIC (Coagulation) | 15 | TF, Thrombin, Fibrin, PAI-1, APC |
| 6 | 내피세포 기능부전 (Endothelial) | 16 | VE-cadherin, ICAM-1, NO, Edema |
| 7 | 다장기부전 (Organ Failure/MODS) | 22 | ARDS, AKI, SOFA 6-domain |
| 8 | 약물 PK/PD (Drug PK/PD) | 15 | Abx, NE, Hydrocortisone, Tocilizumab |
| 9 | 임상 바이오마커 (Biomarkers) | 13 | PCT, CRP, Lactate, SOFA |
| 10 | 대사/미토콘드리아 (Metabolic) | 8 | Lactate, ROS, Mitochondria |
| 11 | 적응면역/면역억제 (Adaptive/CARS) | 11 | Treg, MDSC, PD-1, Apoptosis |

---

## mrgsolve 모델 구획 (24 ODE Compartments)

| 그룹 | 구획 | 설명 |
|------|------|------|
| **병원체** | BACT | 균 부담 (CFU/mL, 로지스틱 성장 + 항생제 사멸) |
| **항생제 PK** | ABX1, ABX2 | 2구획 메로페넴 PK (중심·조직) |
| **사이토카인** | TNF, IL6, IL10, IL1B | 4대 사이토카인 (Emax-억제 피드백) |
| **면역세포** | NEUT_B, NEUT_T, MACS | 혈중/조직 호중구, 활성화 대식세포 |
| **보체** | C5A | C5a 효과기 |
| **응고** | THROMBIN, FIBRIN, PAI1 | 트롬빈·피브린·PAI-1 |
| **내피** | ENDOT | 내피 손상 지수 (0–1) |
| **장기** | PF_RATIO, CREATININE, BILIRUBIN, LACTATE, MAP_val, PLT_COUNT | SOFA 6개 도메인 기반 |
| **약물 PK** | NE_C, HC_C, TOCI_C | 노르에피네프린·하이드로코티손·토실리주맙 |

---

## SOFA 계산 (Sepsis-3 기반)

| 도메인 | 변수 | 0점 | 1점 | 2점 | 3점 | 4점 |
|--------|------|-----|-----|-----|-----|-----|
| 폐 | PaO₂/FiO₂ | ≥400 | <400 | <300 | <200 | <100 |
| 신장 | Creatinine (mg/dL) | <1.2 | 1.2–2.0 | 2.0–3.5 | 3.5–5.0 | >5.0 |
| 간 | Bilirubin (mg/dL) | <1.2 | 1.2–2.0 | 2.0–6.0 | 6.0–12.0 | >12.0 |
| 심혈관 | MAP + vasopressor | ≥70 | <70 | MAP<65 | NE low dose | NE high dose |
| 응고 | Platelet (×10⁹/L) | ≥150 | <150 | <100 | <50 | <20 |
| CNS | 의식/Enc. proxy | 정상 | 경도 | 중등도 | 중증 | 혼수 |
| **합계** | **SOFA 0–24** | — | — | **≥2 = Sepsis** | — | **≥11 ~40% 사망** |

---

## 실행 방법 (How to Run)

```r
# 1) mrgsolve 모델 실행
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
source("sep_mrgsolve_model.R")   # 7개 시나리오 자동 시뮬레이션
plot(p_sofa)                     # SOFA score 비교
plot(p_mort)                     # 28-day mortality
print(summary_72h)               # 72h 요약 표

# 2) Shiny 대시보드 실행
install.packages("shiny", "shinydashboard")
shiny::runApp("sep_shiny_app.R")

# 3) 기계론적 지도 렌더링 (Graphviz)
# dot -Tsvg sep_qsp_model.dot -o sep_qsp_model.svg
# dot -Tpng -Gdpi=150 sep_qsp_model.dot -o sep_qsp_model.png
```

---

## 주요 임상 파라미터 보정 근거 (Parameter Calibration Notes)

| 파라미터 | 값 | 근거 |
|---------|-----|------|
| Bacterial growth rate | kgrow = 1.2/h | E. coli doubling ~35 min in vivo |
| Antibiotic MIC | 0.5 mg/L (meropenem) | EUCAST breakpoint |
| TNFα t1/2 | ~70 min | Beutler 1985 (kdeg=0.6/h) |
| IL-6 t1/2 | ~4.6 h | Taniguchi 1999 (kdeg=0.15/h) |
| Meropenem CL | 10 L/h | Roberts 2009 critical care PK |
| Hydrocortisone Emax | 0.65 (65% cytokine suppression) | Annane 2002 in vitro |
| Tocilizumab EC50 | 1.5 mcg/mL | IL-6R Kd ~1-2 μg/mL |
| SOFA mortality | logit = -6.5 + 0.45×SOFA | Seymour 2017 JAMA external validation |
| Antibiotic delay effect | ~7%/h survival decrease | Kumar 2006 Crit Care Med |

---

*Generated by Claude Code Routine (CCR) — 2026-06-24*  
*For educational and research purposes. Not for clinical use.*
