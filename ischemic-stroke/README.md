# Ischemic Stroke — QSP Model
# 허혈성 뇌졸중 정량적 시스템 약리학(QSP) 모델

[![Disease](https://img.shields.io/badge/Disease-Ischemic%20Stroke-red)]()
[![Category](https://img.shields.io/badge/Category-Neurovascular-blue)]()
[![Compartments](https://img.shields.io/badge/ODE%20Compartments-18-green)]()
[![Scenarios](https://img.shields.io/badge/Treatment%20Scenarios-5-orange)]()

---

## Disease Overview (질환 개요)

**허혈성 뇌졸중(Ischemic Stroke)**은 뇌혈관이 혈전 또는 색전으로 막혀 해당 뇌 영역에 혈류가 차단되어 신경세포가 손상·사멸하는 응급 신경과 질환이다. 전체 뇌졸중의 약 80%를 차지하며, 전 세계 사망 원인 2위, 장기 장애 원인 1위이다.

### Pathophysiology (병태생리)

```
혈전/색전 형성
      ↓
뇌혈관 폐색 → 뇌혈류(CBF) 감소
      ↓
허혈 핵심부 (Core, CBF < 10 mL/100g/min)  →  불가역적 괴사
허혈 주변부 (Penumbra, 10–20 mL/100g/min) →  시간-의존적 회복 가능
      ↓
에너지 실패 → Na⁺/K⁺-ATPase 기능 상실
      ↓
탈분극 → 흥분독성 글루타메이트 방출 → NMDA 수용체 → Ca²⁺ 과부하
      ↓
ROS 생성 · 미토콘드리아 기능부전 · Caspase 활성화
      ↓
신경 염증 (소교세포, 성상세포, 호중구 침윤) → MMP-9 → BBB 파괴
      ↓
혈관성 부종 · 출혈성 전환 · 뇌압 상승
```

### Key Drug Targets (핵심 치료 표적)

| Target | Drug | Indication |
|--------|------|-----------|
| Fibrin/Plasminogen | IV tPA (Alteplase) | Acute thrombolysis (≤4.5h) |
| Mechanical clot | EVT (thrombectomy) | Large vessel occlusion (≤24h, imaging-selected) |
| COX-1 / TXA₂ | Aspirin | Secondary prevention |
| P2Y12 (ADP receptor) | Clopidogrel | Secondary prevention (+ aspirin) |
| Factor Xa | Apixaban / Rivaroxaban | AF-related stroke prevention |
| Thrombin | Dabigatran | AF-related stroke prevention |
| HMG-CoA Reductase | Statins | Plaque stabilization + pleiotropic |
| RAAS / Sympathetic | ACEi / ARB / CCB | BP control, endothelial protection |

---

## Mechanistic Map (기계론적 지도)

[![Mechanistic Map Preview](is_qsp_model.png)](is_qsp_model.svg)

> 클릭하면 고해상도 SVG가 열립니다 / Click to view interactive SVG

**12 Subgraph Clusters · 141 Nodes · 120+ Edges**

| Cluster | Contents |
|---------|----------|
| Risk Factors & Comorbidities | HTN, DM, AF, Dyslipidemia, Smoking, Obesity |
| Vascular Pathology & Thrombosis | Atherosclerosis → Plaque → Platelet → Coagulation |
| Acute Treatment PK/PD | IV tPA (2-compartment), EVT, PAI-1, Fibrinolysis |
| Ischemic Core & Penumbra | CBF dynamics, ATP, PIDs, DWI/PWI mismatch |
| Excitotoxicity & Ion Dysregulation | Glutamate, NMDA/AMPA, Ca²⁺, VGCC, Apoptosis |
| Oxidative Stress & NO | ROS, XO, NADPH oxidase, eNOS/iNOS, Peroxynitrite |
| Neuroinflammation & BBB | Microglia, Astrocytes, Cytokines, MMP-9, BBB |
| Reperfusion Injury | ROS burst, mPTP, No-reflow, Hemorrhagic transformation |
| Secondary Prevention PK/PD | Aspirin, Clopidogrel, Apixaban, Statin, Antihypertensive |
| Neuroprotective & Emerging | NMDA antagonists, Edaravone, BDNF, Stem cell |
| Clinical Outcomes | NIHSS, mRS, Barthel, Mortality, Recurrence |
| Biomarkers | GFAP, UCH-L1, NSE, S100β, IL-6, D-Dimer |

---

## mrgsolve ODE Model (수학적 ODE 모델)

**File:** `is_mrgsolve_model.R`

### Compartments (18 ODEs)

| # | Compartment | Description | Units |
|---|-------------|-------------|-------|
| 1 | THROMBUS | Thrombus burden | 0–1 (normalized) |
| 2 | CBF_CORE | CBF in ischemic core | mL/100g/min |
| 3 | CBF_PEN | CBF in penumbra | mL/100g/min |
| 4 | TPA_CENT | tPA central compartment | mg |
| 5 | TPA_PERI | tPA peripheral compartment | mg |
| 6 | ASP_GUT | Aspirin gut depot | mg |
| 7 | ASP_CENT | Aspirin central | mg |
| 8 | NOAC_GUT | Apixaban gut depot | mg |
| 9 | NOAC_CENT | Apixaban central | mg |
| 10 | NOAC_PERI | Apixaban peripheral | mg |
| 11 | ATP_PEN | Penumbral ATP (normalized) | 0–1 |
| 12 | GLUT | Extracellular glutamate | mmol/L |
| 13 | CA2 | Intracellular Ca²⁺ | mmol/L |
| 14 | ROS | Reactive oxygen species | a.u. |
| 15 | IL6 | Serum IL-6 | pg/mL |
| 16 | BBB | BBB integrity | 0–1 |
| 17 | INFARCT | Infarct core volume | mL |
| 18 | NIHSS | NIHSS score (continuous proxy) | 0–42 |

### 5 Treatment Scenarios

| Scenario | Intervention | Clinical Reference |
|----------|-------------|-------------------|
| 1 | IV tPA at 2h + Aspirin (standard care) | NINDS Trial 1995 |
| 2 | Late IV tPA at 4.5h + Aspirin | ECASS-3 2008 |
| 3 | Antiplatelet only (no thrombolysis) | IST 1997 |
| 4 | IV tPA + Apixaban (AF patient) | ARISTOTLE 2011 |
| 5 | EVT (mechanical thrombectomy) simulation | DEFUSE-3 2018 |

### Key PK Parameters

| Drug | CL | V1 | Half-life | Source |
|------|----|----|-----------|--------|
| tPA (alteplase) | 550 mL/min | 3,500 mL | ~5 min | Tanswell 2002 |
| Aspirin (oral) | 10 L/h | 12 L | ~0.8 h | Levy 1985 |
| Apixaban (oral) | 3.3 L/h | 21 L | ~12 h | Frost 2008 |

---

## Shiny Dashboard (대시보드)

**File:** `is_shiny_app.R`

### 8 Tabs

| Tab | Content |
|-----|---------|
| 1. Patient Profile | Demographics, risk factors, severity indicators |
| 2. Acute Treatment PK | tPA plasma concentration, recanalization kinetics |
| 3. Ischemic Cascade | CBF dynamics, ATP depletion, glutamate/Ca²⁺ |
| 4. Neuroinflammation & BBB | IL-6 kinetics, BBB integrity, edema risk |
| 5. Clinical Endpoints | NIHSS/mRS trajectory, infarct volume, 90d outcome |
| 6. Secondary Prevention | Aspirin/apixaban PK, COX-1/Xa inhibition curves |
| 7. Scenario Comparison | Side-by-side 5-scenario comparison |
| 8. Biomarker Dynamics | GFAP, UCH-L1, NSE, S100β, D-Dimer proxies |

---

## Usage (실행 방법)

### Render Mechanistic Map
```bash
# SVG (vector)
dot -Tsvg is_qsp_model.dot -o is_qsp_model.svg

# PNG (150 dpi)
dot -Tpng -Gdpi=150 is_qsp_model.dot -o is_qsp_model.png
```

### Run mrgsolve ODE Model
```r
install.packages(c("mrgsolve", "dplyr", "ggplot2", "patchwork"))
source("ischemic-stroke/is_mrgsolve_model.R")
# → 5-scenario summary table + 6-panel plot
```

### Run Shiny Dashboard
```r
install.packages(c("shiny", "shinydashboard", "plotly", "DT", "mrgsolve", "dplyr", "ggplot2"))
shiny::runApp("ischemic-stroke/is_shiny_app.R")
```

---

## Key Model Insights (모델의 주요 시사점)

1. **Time is brain:** tPA at 2h vs. 4.5h — earlier treatment saves more penumbra. Each 1h delay ≈ 1.9 million neurons lost.
2. **EVT superiority:** Mechanical thrombectomy achieves near-complete recanalization in hours, especially for large vessel occlusion (DAWN/DEFUSE-3 window up to 24h).
3. **BBB disruption peaks at 24–48h** — this window coincides with the highest hemorrhagic transformation risk; MMP-9–driven model captures this dynamics.
4. **Secondary prevention critical:** Aspirin (IST) + statin (SPARCL) reduce recurrence by ~25–30%; apixaban for AF prevents ~50% of cardioembolic recurrences (ARISTOTLE).
5. **Penumbra salvage window is narrow:** ATP in penumbra begins recovering within 30 min of recanalization; untreated, it converts to infarct at ~0.15/h.

---

## References (참고문헌)

Full reference list with 50 PubMed citations: [`is_references.md`](is_references.md)

Key references: NINDS 1995 · ECASS-3 2008 · DEFUSE-3 2018 · ARISTOTLE 2011 · SPARCL 2006 · Tanswell 2002

---

*Generated by Claude Code Routine (CCR) | 2026-06-25*
