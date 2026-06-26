# 당뇨병성 망막병증 (Diabetic Retinopathy) QSP 모델

> **디렉토리:** `diabetic-retinopathy/` | **약어:** DR | **날짜:** 2026-06-24  
> **분류:** 만성질환 / 안과 / 당뇨합병증

---

[![DR QSP 기계론적 지도](dr_qsp_model.png)](dr_qsp_model.svg)

---

## 질환 개요

**당뇨병성 망막병증(Diabetic Retinopathy, DR)**은 당뇨병의 가장 흔한 미세혈관 합병증으로, 전 세계 성인 실명의 주요 원인입니다.

| 항목 | 내용 |
|------|------|
| **유병률** | 전체 당뇨 환자의 약 34.6% (~1억 4,600만 명) |
| **DME** | 성인 노동 연령층 실명의 주원인; 당뇨 환자의 ~7% |
| **PDR 위험** | 중증 NPDR에서 1년 내 ~50%가 PDR로 진행 |
| **전 세계 부담** | 2045년까지 ~2억 2,400만 명 DR 환자 예상 (Teo 2021) |
| **임상 결과 지표** | ETDRS 시력 (글자수), OCT CRT (µm), DR 중증도 단계 |

---

## 주요 병태생리 경로

### 1단계: 과혈당 → 4가지 생화학 경로
| 경로 | 핵심 효소/분자 | 결과 |
|------|--------------|------|
| **폴리올 경로** | 알도스환원효소(AR) | 소르비톨↑, NADPH 고갈, 유사저산소증 |
| **헥소사민 경로** | GFAT, O-GlcNAc | Sp1 당화 → TGF-β, PAI-1 ↑ |
| **PKC 경로** | PKCβ1/β2/δ | VEGF↑, NF-κB↑, eNOS↓, ET-1↑ |
| **AGE-RAGE 경로** | 메틸글리옥살 | 단백질 교차결합, 기저막 비후, NF-κB↑ |

### 2단계: 산화-니트로화 스트레스
- 미토콘드리아 ETC 과부하 → O₂•⁻ 생성 → PARP 활성화 → GAPDH 억제 → 경로 1–4 증폭 피드백
- BH4 고갈 → eNOS 탈결합 → NO 감소 + ONOO⁻ 증가 → 세포자사

### 3단계: VEGF/신생혈관 형성
- HIF-1α(저산소) + NF-κB(염증) → VEGF-A165 과발현
- VEGFR2 신호 → PI3K/AKT(투과성↑) + ERK(증식↑) + PLCγ(PKC)
- Ang2 증가 / Tie-2 탈안정화 → 주피세포 지지 소실 → 혈관 취약

### 4단계: 혈관 구조 병변
| 병변 | 단계 | 기전 |
|------|------|------|
| 미세동맥류 | 초기 NPDR | 주피세포 소실 → 혈관 팽창 |
| 경성 삼출물 | NPDR | 지질 누출 |
| 면화반 | 중증 NPDR | 신경섬유 허혈 |
| IRMA | 중증 NPDR | 모세혈관 폐쇄 → 측부 순환 |
| 신생혈관(NVE/NVD) | PDR | VEGF 구동 혈관신생 |
| 당뇨황반부종(DME) | 모든 단계 | BRB 파괴 → 황반 내 액체 |

---

## QSP 모델 파일

| 파일 | 설명 |
|------|------|
| [`dr_qsp_model.dot`](dr_qsp_model.dot) | Graphviz 기계론적 지도 (소스) |
| [`dr_qsp_model.svg`](dr_qsp_model.svg) | SVG 벡터 이미지 |
| [`dr_qsp_model.png`](dr_qsp_model.png) | PNG 이미지 (150 dpi) |
| [`dr_mrgsolve_model.R`](dr_mrgsolve_model.R) | mrgsolve ODE 모델 |
| [`dr_shiny_app.R`](dr_shiny_app.R) | Shiny 대시보드 |
| [`dr_references.md`](dr_references.md) | 참고문헌 57개 |

---

## 모델 사양

### 기계론적 지도 (DOT)
- **노드 수:** 210+ (9개 클러스터)
- **클러스터:**
  1. 전신 위험인자 (Systemic Risk Factors)
  2. 과혈당 유발 생화학 경로 (Biochemical Pathways)
  3. 산화-니트로화 스트레스 (Oxidative Stress)
  4. VEGF/혈관신생 신호 (VEGF Signaling)
  5. 신경염증 경로 (Neuroinflammation)
  6. 망막 혈관 구조 병변 (Retinal Vascular Pathology)
  7. 망막 신경퇴행 (Retinal Neurodegeneration)
  8. 약물 PK/PD (Drug PK/PD)
  9. 임상 엔드포인트 (Clinical Endpoints)

### mrgsolve ODE 모델
- **구획:** 18개 ODE
  - 약물 PK: `DRUG_VIT`, `DRUG_CENT`, `DRUG_PERIPH`, `CORT_VIT`
  - 혈당: `BG`, `HBA1C`
  - VEGF: `VEGF_FREE`, `VEGF_BOUND`, `VEGF_PLANT` (PlGF)
  - 산화 스트레스: `ROS`, `AGE`
  - 염증: `CYT`, `ICAM`
  - 세포: `PERICYTE`, `EC_COUNT`
  - 구조: `PERM`, `NV`, `CRT`
  - 시력: `VA`
- **치료 시나리오:** 6개
  | 시나리오 | 설명 | 근거 임상시험 |
  |---------|------|-------------|
  | S0 | 무치료 (불량 혈당 조절) | DCCT 대조군 |
  | S1 | 혈당 조절만 (HbA1c → 7%) | DCCT/EDIC |
  | S2 | 아플리버셉트 2mg IVT q4w×5→q8w | PROTOCOL T, PANORAMA |
  | S3 | 라니비주맙 0.5mg IVT q4w | RISE/RIDE |
  | S4 | 파리시맙 6mg IVT q4w×4→q16w | TENAYA/LUCERNE |
  | S5 | 아플리버셉트 + 혈당 조절 병용 | CLARITY + DCCT 기반 |

### Shiny 대시보드 (8탭)
| 탭 | 내용 |
|----|------|
| ① Patient Profile | 환자 파라미터 설정, DR 단계 자동 분류, 위험도 추정 |
| ② Drug PK | 유리체 내 약물 농도, PK 파라미터, AUC/Cmax |
| ③ VEGF / Angiogenesis | 자유 VEGF, 약물-VEGF 복합체, NV 지수, 주피세포 |
| ④ Oxidative/Inflammation | ROS, AGE, 사이토카인, ICAM-1 |
| ⑤ Retinal Structure | CRT (OCT), 혈관 투과성, 내피세포 수 |
| ⑥ Visual Outcomes | BCVA (ETDRS), VA 변화량, 결과 요약 |
| ⑦ Scenario Comparison | 6개 시나리오 동시 비교 플롯 + 엔드포인트 표 |
| ⑧ Biomarkers & About | 바이오마커 프로파일, 모델 개요, 보정 데이터 |

---

## 실행 방법

```bash
# 1) 기계론적 지도 재렌더링 (Graphviz 필요)
dot -Tsvg dr_qsp_model.dot -o dr_qsp_model.svg
dot -Tpng -Gdpi=150 dr_qsp_model.dot -o dr_qsp_model.png
```

```r
# 2) mrgsolve 모델 실행 (R 패키지 필요)
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
source("dr_mrgsolve_model.R")

# 3) Shiny 대시보드 실행
install.packages(c("shiny", "shinydashboard", "plotly", "DT"))
shiny::runApp("dr_shiny_app.R")
```

---

## 주요 임상시험 요약

| 시험명 | 약물 | N | 기간 | 1차 결과 |
|--------|------|---|------|---------|
| **PROTOCOL T** | AFL vs RBZ vs Bev | 660 | 1년 | VA: +13.3/+11.2/+9.7 글자 |
| **RISE/RIDE** | Ranibizumab | 382+382 | 2년 | VA: +10.9 글자 |
| **CLARITY** | Aflibercept vs PRP | 232 | 1년 | VA: +3.3 글자 (AFL 우세) |
| **PANORAMA** | Aflibercept | 402 | 2년 | 2-step 개선: 65% vs 15% |
| **TENAYA** | Faricimab | 331+327 | 1년 | VA: +5.8 글자; CRT −189µm |
| **LUCERNE** | Faricimab | 330+338 | 1년 | VA: +6.6 글자; CRT −194µm |
| **DCCT** | 인슐린 집중 치료 | 1,441 | 6.5년 | 신규 DR 76% 감소 |

---

## 참고문헌
총 **57개** PubMed 인용 → [`dr_references.md`](dr_references.md)

- 역학/분류 (4편)
- 과혈당 경로 (6편)
- VEGF/혈관신생 (6편)
- 항-VEGF 임상시험 (8편)
- 혈당 조절 (4편)
- 산화 스트레스/AGE (4편)
- 신경퇴행 (3편)
- PK/QSP 모델링 (5편)
- 인플라마솜/염증 (3편)
- Ang/Tie2 신호 (4편)
- OCT 바이오마커 (3편)
- 스테로이드 치료 (2편)
- 신경보호/신규 치료제 (3편)
- GLP-1 작용제 (2편)
