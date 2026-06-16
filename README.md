# qsp

## mrgsolve

- <https://vantage-research.net/qsp-in-r/>
- gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve
    - <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/>
    - <https://github.com/Genentech/gPKPDviz/>
    



## iqrtools

- <https://www.intiquan.com/acop2019_qsp/>

---

## QSP Disease Models

Each model directory contains:
- **`.dot`** — Graphviz mechanistic map source (>100 nodes)
- **`.svg` / `.png`** — Rendered pathway diagram
- **`_references.md`** — Curated PubMed literature
- **`_mrgsolve_model.R`** — mrgsolve ODE model + simulation scenarios
- **`_shiny_app.R`** — Interactive Shiny dashboard skeleton

---

| Date | 질환 (Disease) | 분류 | 메커니즘 맵 | 모델 파일 | 레퍼런스 | Shiny 앱 |
|------|---------------|------|------------|----------|----------|----------|
| 2026-06-16 | **[폐동맥 고혈압 (PAH)](pulmonary-arterial-hypertension/)** | 만성질환 / 폐혈관 | [![PAH map](pulmonary-arterial-hypertension/pah_qsp_model.png)](pulmonary-arterial-hypertension/pah_qsp_model.svg) | [pah_mrgsolve_model.R](pulmonary-arterial-hypertension/pah_mrgsolve_model.R) | [References](pulmonary-arterial-hypertension/pah_references.md) | [pah_shiny_app.R](pulmonary-arterial-hypertension/pah_shiny_app.R) |

### 폐동맥 고혈압 (Pulmonary Arterial Hypertension, PAH)

> 폐동맥 내피세포 기능 부전 및 혈관 평활근 비대로 폐혈관 저항이 상승하여 우심실 부전을 유발하는 희귀 진행성 질환.

**모델 구성 요소 (>100 nodes):**

| 축 | 주요 성분 |
|----|----------|
| 내피소 축 (ET-1) | ET-1 → ETA/ETB → Gq → PLC → IP₃/DAG → Ca²⁺ → MLCK → 혈관수축 |
| NO-cGMP 축 | eNOS → NO → sGC → cGMP → PKG → MLCP 활성화 → 혈관이완 |
| PGI₂-cAMP 축 | 프로스타시클린 → IP수용체 → Gs → AC → cAMP → PKA → 혈관이완 |
| 성장인자 축 | BMPR2↓ / TGF-β / PDGF / VEGF / FGF → ERK/AKT/mTOR/STAT3 → PASMC 증식·생존 |
| 면역/염증 | Th1/Th2/Th17/Treg, M1/M2 대식세포, IL-6/IL-1β/TNF-α → NF-κB |
| 대사/저산소 | HIF-1α/HIF-2α → Warburg 효과, 미토콘드리아 이상, Kv1.5↓, ROS↑ |
| 우심실 모델 | RV-PA coupling (Ees/Ea), Frank-Starling, BNP, TAPSE |
| 임상 지표 | PVR, mPAP, CO, 6MWD, WHO-FC, BNP, 생존율 |
| 약물 PK/PD | ERA (bosentan/ambrisentan/macitentan), PDE5i (sildenafil/tadalafil), sGC (riociguat), PGI₂ (epoprostenol/treprostinil/selexipag) |

**파일 목록:**

| 파일 | 설명 |
|------|------|
| [`pah_qsp_model.dot`](pulmonary-arterial-hypertension/pah_qsp_model.dot) | Graphviz DOT 소스 (130+ 노드, 8개 서브그래프) |
| [`pah_qsp_model.svg`](pulmonary-arterial-hypertension/pah_qsp_model.svg) | 벡터 경로 다이어그램 |
| [`pah_qsp_model.png`](pulmonary-arterial-hypertension/pah_qsp_model.png) | 래스터 이미지 (150 dpi) |
| [`pah_references.md`](pulmonary-arterial-hypertension/pah_references.md) | 40개 핵심 논문 (PubMed 링크 포함) |
| [`pah_mrgsolve_model.R`](pulmonary-arterial-hypertension/pah_mrgsolve_model.R) | mrgsolve ODE 모델 + 6가지 치료 시나리오 시뮬레이션 |
| [`pah_shiny_app.R`](pulmonary-arterial-hypertension/pah_shiny_app.R) | Shiny 대시보드 (PK/PD, 용량-반응, 위험도 평가) |
