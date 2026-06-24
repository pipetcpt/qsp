# Kawasaki Disease (KD) — QSP Model

[![Model](kd_qsp_model.png)](kd_qsp_model.svg)

## 개요 (Overview)

**가와사키병(Kawasaki Disease, KD)**은 주로 5세 미만 소아에서 발생하는 원인 미상의 급성 전신 혈관염(vasculitis)으로, 선진국에서 **소아 후천성 심장질환의 가장 흔한 원인**입니다. 적절히 치료받지 못할 경우 관상동맥류(coronary artery aneurysm, CAA) 형성으로 이어져 평생 허혈성 심질환 위험을 초래합니다.

**병태생리 핵심:**  
미확인 감염 트리거 → 선천 면역 활성화(대식세포, 중성구) → NLRP3 인플라마좀 → 사이토카인 폭풍(IL-1β, IL-6, TNF-α) → 혈관 내피 손상 → 관상동맥 중막 파괴 → 동맥류 형성

---

## 질환 특성 (Disease Features)

| 항목 | 내용 |
|------|------|
| **분류** | 소아 전신 혈관염 (Pediatric Systemic Vasculitis) |
| **유병률** | 일본 ~ 300/100,000명, 한국 ~ 200/100,000명, 미국 ~ 20/100,000명 (5세 미만) |
| **호발 연령** | 6개월–5세 (median ~ 2세) |
| **진단 기준** | AHA 2017: 5일 이상 발열 + 5가지 임상 증상 중 4개 이상 |
| **CAA 발생률** | 치료 없이 ~ 25%, IVIG 후 ~ 3–5% |
| **IVIG 저항성** | 약 10–20% |
| **재발률** | 약 3% |

### 핵심 경로 (Key Pathways)

```
Unknown Trigger
  ↓
TLR/NLR Activation → Macrophage Activation
  ↓
NLRP3 Inflammasome → Caspase-1 → IL-1β (mature)
  ↓
Cytokine Storm: IL-1β + IL-6 + TNF-α
  ↓
Endothelial Activation (VCAM-1, ICAM-1, Tissue Factor ↑)
  ↓
Coronary Artery Inflammation → Medial Destruction → Aneurysm
  ↓
Thrombocytosis (peak Week 2–3) → Thrombosis Risk
  ↓
CAA: small (z 2.5–5) → medium (z 5–10) → giant (z ≥ 10)
```

---

## 치료 알고리즘 (Treatment Algorithm)

```
Diagnosis confirmed
├── Standard-risk (Kobayashi score < 4)
│   └── IVIG 2 g/kg × 1 dose + Aspirin 80–100 mg/kg/day
│       ├── Responder (80–90%) → Switch to low-dose aspirin (3–5 mg/kg/day)
│       └── Non-responder → IVIG-resistant protocol
│
├── High-risk (Kobayashi score ≥ 4 or Egami score ≥ 3)
│   └── IVIG + Aspirin + Prednisolone 2 mg/kg/day × 5–15 days
│
└── IVIG-Resistant (~10–20%)
    ├── Option A: 2nd IVIG 2 g/kg
    ├── Option B: Infliximab 5 mg/kg IV (TNF-α blockade)
    ├── Option C: IV Methylprednisolone pulse
    └── Option D: Anakinra 4 mg/kg/day SC (IL-1 blockade)
```

---

## QSP 모델 구조 (Model Structure)

### 기계론적 지도 (Mechanistic Map)

| 구성 요소 | 내용 |
|-----------|------|
| **총 노드 수** | 134개 |
| **클러스터 수** | 14개 |
| **클러스터 목록** | 감염 트리거, 선천 면역, NLRP3 인플라마좀, 사이토카인 네트워크, 적응 면역, 내피세포 활성화, 관상동맥 병리, 혈소판 생물학, 발열·급성기 반응, IVIG PK, 아스피린 PK, 코르티코스테로이드 PK, 생물학적 제제 PK, 임상 엔드포인트 |

### mrgsolve ODE 모델

| 구획 | 설명 |
|------|------|
| A_IVIG_c/p | IVIG 중심/말초 구획 (2-구획 + FcRn 재순환) |
| A_ASA_gut/c, A_SA_c | 아스피린 흡수/혈장/살리실산 대사체 |
| A_MP_c/p | 메틸프레드니솔론 2-구획 |
| A_IFX_c/p | 인플릭시맙 2-구획 |
| A_ANK_gut/c | 아나킨라 SC 흡수 + 혈장 |
| IL1b, IL6, TNFa | 사이토카인 동역학 |
| Mac_act, EC_act | 대식세포·내피세포 활성 |
| Fever | 체온 (발열 역학) |
| CRP | C-반응성 단백 |
| PLT_c | 혈소판 수 (혈소판 증가증) |
| CAL_Z | 관상동맥 Z-점수 |

**총 21개 구획, 5가지 치료 시나리오**

### 5가지 치료 시나리오

| 시나리오 | 설명 |
|----------|------|
| **S1** | IVIG 2 g/kg + 고용량 아스피린 → 저용량 아스피린 (표준 치료) |
| **S2** | IVIG + 아스피린 + 메틸프레드니솔론 (Kobayashi 고위험군) |
| **S3** | IVIG 저항성 → 2차 IVIG 투여 |
| **S4** | IVIG 저항성 → 인플릭시맙 구제 요법 |
| **S5** | IVIG 저항성 → 아나킨라 구제 요법 (IL-1 차단) |

### Shiny 대시보드 (6개 탭)

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 진단 기준, 치료 알고리즘, Value Box |
| 2. 약물 동태 (PK) | IVIG / 아스피린 / 스테로이드 / 생물학적 제제 PK 곡선 |
| 3. 사이토카인/염증 | IL-1β, IL-6, TNF-α, 대식세포·내피세포 활성화 |
| 4. 임상 엔드포인트 | 발열, CRP, 관상동맥 Z-점수, 혈소판 |
| 5. 시나리오 비교 | 5가지 치료 시나리오 비교 |
| 6. 바이오마커/위험도 | Kobayashi 점수 계산기, CAA 위험 확률 |

---

## 파일 목록 (Files)

| 파일 | 설명 |
|------|------|
| `kd_qsp_model.dot` | Graphviz 기계론적 지도 소스 (134 노드, 14 클러스터) |
| `kd_qsp_model.svg` | 벡터 이미지 |
| `kd_qsp_model.png` | 래스터 이미지 (150 dpi) |
| `kd_mrgsolve_model.R` | mrgsolve ODE QSP 모델 (5 시나리오) |
| `kd_shiny_app.R` | Shiny 인터랙티브 대시보드 (6 탭) |
| `kd_references.md` | 참고문헌 (60개 PubMed 링크, 14개 섹션) |

---

## 실행 방법 (Usage)

```bash
# 기계론적 지도 렌더링
dot -Tsvg kd_qsp_model.dot -o kd_qsp_model.svg
dot -Tpng -Gdpi=150 kd_qsp_model.dot -o kd_qsp_model.png
```

```r
# mrgsolve 모델 실행
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
source("kd_mrgsolve_model.R")
results <- run_all_scenarios(wt = 15)

# Shiny 앱 실행
install.packages(c("shiny", "shinydashboard", "plotly", "DT"))
shiny::runApp("kd_shiny_app.R")
```

---

## 핵심 파라미터 (Key Parameters)

| 파라미터 | 값 | 출처 |
|----------|----|------|
| IVIG 용량 | 2 g/kg IV | AHA 2017 |
| IVIG t½ | 21–28일 | Tremoulet 2015 |
| Infliximab EC50 | 2.5 μg/mL | Tremoulet 2020 |
| Anakinra EC50 | 1.0 μg/mL | Ouldali 2019 |
| IL-6 기저치 | 2.5 ng/mL | Matsubara 2013 |
| 관상동맥 Z-점수 임계값 | ≥ 2.5 (CAL), ≥ 10 (거대 동맥류) | AHA 2017 |

---

## 참고문헌 요약 (References Summary)

60개 참고문헌 — 14개 섹션: 역학(5) · 병태생리(5) · 사이토카인(6) · 관상동맥 병리(5) · IVIG 약리(6) · 위험 점수(4) · 코르티코스테로이드(3) · 인플릭시맙(3) · 아나킨라/생물학적 제제(3) · 혈소판/아스피린(3) · 심초음파/Z-점수(3) · QSP 모델링(4) · COVID-19/MIS-C(3) · 장기 예후(5)

자세한 내용은 [`kd_references.md`](kd_references.md) 참조.
