# Primary Open-Angle Glaucoma (POAG) — QSP Model

**Abbreviation:** POAG · **Category:** 안과 (Ophthalmology / Chronic Disease)
**Date Created:** 2026-06-27 · **Directory:** `glaucoma/`

---

## 개요 (Overview)

원발개방각녹내장(POAG)은 전 세계 비가역적 실명의 주요 원인으로,
방수(aqueous humor) 유출 장애로 인한 안압(IOP) 상승이 시신경 손상을 유발하는 만성 진행성 안질환입니다.
이 QSP 모델은 방수역학, 섬유주망(TM) 생물학, 시신경 유두(ONH) 손상, 망막신경절세포(RGC) 사멸 경로를
통합하여 5가지 약물군 및 레이저·수술 치료의 IOP 및 시야 보존 효과를 정량적으로 예측합니다.

Primary open-angle glaucoma (POAG) is the leading cause of irreversible blindness worldwide.
This QSP model integrates aqueous humor dynamics (Goldmann equation), trabecular meshwork (TM)
biology (ECM accumulation, ROCK signaling), optic nerve head (ONH) biomechanics (lamina cribrosa
deformation, BDNF retrograde transport), and retinal ganglion cell (RGC) apoptosis to quantitatively
predict the long-term effects of five drug classes and surgical interventions on IOP and visual field
(VF-MD) preservation.

---

## 기계론적 지도 (Mechanistic Map)

[![POAG QSP Model](poag_qsp_model.png)](poag_qsp_model.svg)

> Click the image to open the interactive SVG map.
> Rendered with `sfdp` layout engine (Graphviz 2.42).

### 10개 서브그래프 클러스터:
| # | 클러스터 | 주요 구성 요소 |
|---|---------|--------------|
| 1 | 유전·위험인자 | MYOC, OPTN, CDKN2B-AS1, CAV1/2, LOXL1 (XFS) |
| 2 | 모양체·방수 생성 | CA-II/IV, Na⁺/K⁺-ATPase, β₂-AR, α₂-AR, AQP1 |
| 3 | 섬유주·전통 유출로 | RhoA→ROCK→actin, TGF-β2→Smad2/3→ECM, TM senescence |
| 4 | 포도막공막 유출로 | FP-R→PKC→MMP-1,3→ciliary ECM remodeling |
| 5 | 안압 역학 | Goldmann equation, diurnal variation, OHT, NTG |
| 6 | 약물 PK/PD | PGA, BB, CAI, A2A, ROCK-I, SLT, MIGS, Trabeculectomy |
| 7 | 시신경 유두 생물학 | LC biomechanics, TPG, axonal transport, BDNF, TrkB/p75NTR |
| 8 | RGC 사멸 경로 | NMDA-Ca²⁺-nNOS, ROS, Bax/Bcl-2→Casp9→Casp3, DLK→JNK |
| 9 | 임상 종점·모니터링 | RNFL, GCL, VF-MD, VFI, CDR, GPA |
| 10 | 전신·혈관 인자 | OPP, autoregulation, XFS, PDG, neuroprotection |

---

## mrgsolve ODE 모델 (mrgsolve Model)

**파일:** [`poag_mrgsolve_model.R`](poag_mrgsolve_model.R)

### 구획 (16 ODE Compartments):
| # | 상태변수 | 설명 |
|---|---------|------|
| 1–5 | C_PGA, C_BB, C_CAI, C_A2A, C_ROCK | Drug PK in aqueous (ng/mL) |
| 6 | F_aq | Aqueous production rate (μL/min) |
| 7 | F_uv | Uveoscleral outflow rate (μL/min) |
| 8 | C_tm | TM outflow facility (μL/min/mmHg) |
| 9 | IOP | Intraocular pressure (mmHg) — Goldmann dynamic |
| 10 | ECM_TM | ECM accumulation index in TM |
| 11 | BDNF | BDNF in ONH (pg/mL) |
| 12 | Casp3 | Caspase-3 apoptotic index (0–1) |
| 13 | RGC | Retinal ganglion cell count (millions) |
| 14 | RNFL | RNFL thickness (μm) |
| 15 | VF_MD | Visual field mean deviation (dB) |

### 치료 시나리오 (8 Treatment Scenarios):
1. **Untreated POAG** — natural history (IOP ≈ 24 mmHg)
2. **Latanoprost QD** — PGA, ↑F_uv +100%
3. **Timolol BID** — β₂-blocker, ↓F_prod 30%
4. **Dorzolamide TID** — CAI, ↓F_prod 25%
5. **Brimonidine BID** — A2A + neuroprotection (↑BDNF)
6. **Netarsudil QD** — ROCK-I, ↑C_tm 35% + ↓P_ep 3 mmHg
7. **Latanoprost + Timolol** — Fixed-dose combination
8. **Triple therapy** — PGA + BB + CAI

### 핵심 수식 (Key Equations):
```
Goldmann:  IOP_eq = (F_aq - F_uv) / C_tm + P_ep
TM ECM:    dECM/dt = k_prod + k_IOP·max(IOP−18,0) − k_clear·ECM − ROCK_effect
BDNF:      dBDNF/dt = k_prod·(1+A2A_boost) − (k_deg + k_IOP·max(IOP−18,0))·BDNF
RGC loss:  dRGC/dt = −(k_base + k_Casp3·Casp3_eff) · RGC
VF-MD:     VF_MD_eq = −25·(1 − RGC/RGC₀)^2.5  dB
```

---

## Shiny 앱 (Shiny App)

**파일:** [`poag_shiny_app.R`](poag_shiny_app.R)

### 7개 탭 구성:
| 탭 | 내용 |
|----|------|
| 1. Patient Profile | Demographics, diagnosis stage, baseline IOP/RNFL/VF-MD |
| 2. Drug PK | Aqueous concentration profiles for 5 drug classes |
| 3. IOP Dynamics | Real-time IOP, F_aq, F_uv, C_tm with treatment |
| 4. Clinical Endpoints | VF-MD, RNFL, RGC count over 10 years |
| 5. Scenario Comparison | 8-scenario comparison with IOP, VF, RNFL tables |
| 6. Biomarkers & Neuroprotection | BDNF, Caspase-3, ECM_TM dynamics |
| 7. Sensitivity Analysis | Tornado plot, IOP vs VF scatter, parameter sensitivity |

---

## 참고문헌 (References)

**파일:** [`poag_references.md`](poag_references.md) — **66개** 문헌

### 주요 임상시험:
- **AGIS** (2000) — IOP control prevents VF deterioration
- **OHTS** (2002) — Ocular Hypertension Treatment Study
- **EMGT** (2003) — Early Manifest Glaucoma Trial
- **CIGTS** (2001–2011) — Collaborative Initial Glaucoma Treatment Study
- **LiGHT** (2019) — SLT first-line vs. eye drops (Lancet)
- **UKGTS** (2015) — Latanoprost RCT (Lancet)
- **TVT** (2012) — Tube vs. Trabeculectomy

---

## 파라미터 출처 (Key Parameter Sources)

| 파라미터 | 값 | 출처 |
|---------|---|------|
| F_prod baseline | 2.5 μL/min | Brubaker 1991 IOVS |
| F_uv baseline | 0.4 μL/min | Brubaker 1991 IOVS |
| C_tm normal | 0.30 μL/min/mmHg | Goldmann 1951 |
| P_ep | 8.0 mmHg | Goldmann 1951 |
| PGA ↑F_uv | +80–120% | Stjernschantz 2001 IOVS |
| BB ↓F_prod | 25–35% | Liu 2003 IOVS |
| ROCK-I ↑C_tm | +30–40% | Tian 2005 Br J Ophthalmol |
| RGC normal loss | 0.5%/yr | Harwerth 2004 IOVS |

---

## 실행 방법 (Usage)

```r
# mrgsolve 시뮬레이션
source("poag_mrgsolve_model.R")

# Shiny 앱 실행
library(shiny)
runApp("poag_shiny_app.R")
```

### 필요 패키지:
```r
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr",
                   "shiny", "shinydashboard", "plotly", "DT"))
```

---

## 파일 목록 (Files)

| 파일 | 설명 |
|------|------|
| `poag_qsp_model.dot` | Graphviz DOT — 150+ 노드, 10 서브그래프 |
| `poag_qsp_model.svg` | 벡터 이미지 (인터랙티브) |
| `poag_qsp_model.png` | 래스터 이미지 (150 dpi) |
| `poag_mrgsolve_model.R` | ODE 모델 + 8개 시나리오 시뮬레이션 |
| `poag_shiny_app.R` | 7탭 Shiny 대시보드 |
| `poag_references.md` | 66개 문헌 목록 |
| `README.md` | 이 파일 |
