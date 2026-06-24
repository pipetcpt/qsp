# Multiple Myeloma (MM) QSP Model

[![MM QSP Mechanistic Map](mm_qsp_model.png)](mm_qsp_model.svg)

## 개요 (Overview)

다발성 골수종(Multiple Myeloma, MM)은 골수 내 형질세포(plasma cell)의 악성 클론이 증식하는 혈액 종양으로, 미국에서 연간 약 35,000명이 새로 진단되는 두 번째로 흔한 혈액암입니다. 오스테오리틱 골 병변, M-단백(M-protein) 과분비, 신부전, 빈혈, 고칼슘혈증(CRAB 기준)이 특징적입니다.

본 QSP 모델은 **골수 미세환경(BMME) · 종양세포 역학 · 골 리모델링 · 약물 PK/PD**를 통합하는 13-클러스터 기계론적 지도와 mrgsolve ODE 시뮬레이션을 포함합니다.

---

## 주요 병태생리 경로 (Key Pathophysiology)

| 경로 | 핵심 메커니즘 | 임상 결과 |
|------|-------------|---------|
| Proteasome (UPS) 과부하 | 형질세포의 대량 면역글로불린 합성 → ER stress → UPS 과활성 | 약물 표적 (BTZ/CFZ) |
| IL-6 / JAK-STAT3 | BMSC가 IL-6 분비 → MM cell 증식 · 항아포토시스 | M-protein↑, Hgb↓ |
| RANKL/OPG 불균형 | MM cell이 RANKL↑·OPG↓ 유도 → OC 과활성 | 골용해 병변, 병리적 골절 |
| DKK1 / Wnt 억제 | MM cell이 DKK1·Sclerostin 분비 → OB 기능 억제 | 골 형성 억제 → 순수 골 소실 |
| IKZF1/IKZF3 (Ikaros/Aiolos) | IMiD가 Cereblon(CRBN) E3 ligase 탈선 → IKZF1/3 분해 → IRF4↓ | LEN/POM 효능 기전 |
| CD38 / 면역 | 형질세포 특이 CD38 과발현 → DARA ADCC/CDC 활성 | 다라투무맙 타겟 |
| BCL-2 / BH3 mimetics | t(11;14) 환자 BCL-2 의존성↑ → VEN이 BIM 방출 | venetoclax 감수성 |
| PD-1/PD-L1 축 | MM cell PD-L1 발현 → CTL 소진 | 면역 회피 |

---

## 파일 구성 (File List)

| 파일 | 설명 |
|------|------|
| [`mm_qsp_model.dot`](mm_qsp_model.dot) | Graphviz 기계론적 지도 (13 클러스터, 130+ 노드) |
| [`mm_qsp_model.svg`](mm_qsp_model.svg) | SVG 벡터 이미지 (클릭하면 확대 가능) |
| [`mm_qsp_model.png`](mm_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [`mm_mrgsolve_model.R`](mm_mrgsolve_model.R) | mrgsolve ODE 모델 (30 구획, 6 치료 시나리오) |
| [`mm_shiny_app.R`](mm_shiny_app.R) | Shiny 대시보드 (6탭 인터랙티브) |
| [`mm_references.md`](mm_references.md) | 참고문헌 50편 (PubMed 링크 포함) |

---

## 기계론적 지도 클러스터 (Mechanistic Map Clusters)

1. **Plasma Cell Differentiation** — HSC → CLP → B cell → GC → Plasma cell
2. **Myeloma Cell Biology** — MGUS → SMM → MM cell; 종양 유전체 이상
3. **Bone Marrow Microenvironment (BMME)** — BMSC, BMEC, Osteoblast/Osteoclast, MSC, 접착분자
4. **Cytokine & Growth Factor Signaling** — IL-6, APRIL/BAFF, RANKL/OPG, DKK1, MIP-1α
5. **Intracellular Signaling (MM Cell)** — JAK-STAT3, PI3K/AKT/mTOR, NF-κB, Proteasome, BCL-2 family
6. **Immune Evasion** — CTL, NK cell, Treg, MDSC, PD-1/PD-L1, CD38, BCMA, SLAMF7
7. **Myeloma Bone Disease** — Osteolytic lesion, NTX/PINP biomarkers, Hypercalcemia, SRE
8. **M-Protein & Clinical Markers** — M-protein, sFLC, β₂M, LDH, Albumin, BMPC%
9. **Drug PK** — Bortezomib (3-cmt), Carfilzomib (2-cmt), Lenalidomide, Daratumumab (TMDD), Dexamethasone, Venetoclax, Zoledronic acid
10. **Drug PD** — Proteasome inhibition, CRBN/IKZF1/3 degradation, ADCC/CDC, BCL-2 BH3 mimetics, CAR-T/BiTE
11. **Clinical Endpoints** — ISS/R-ISS staging, CR/VGPR/PR/SD/PD, MRD negativity, PFS, OS
12. **Renal Complications** — Cast nephropathy, AL amyloidosis, eGFR
13. **Hematologic Complications** — Anemia, Neutropenia, Thrombocytopenia

---

## mrgsolve 모델 구획 (Model Compartments)

| 카테고리 | 구획 | 수 |
|---------|------|---|
| 종양세포 | MM_S (감수성), MM_R (내성) | 2 |
| 질환 마커 | MP (M-protein), FLC (sFLC), IL6, VEGF | 4 |
| 골 리모델링 | OB, OC, BV, NTX, PINP | 5 |
| 전신 바이오마커 | Hgb, B2M | 2 |
| Bortezomib PK | BTZ1, BTZ2, BTZ3 | 3 |
| Lenalidomide PK | LEN_gut, LEN1 | 2 |
| Daratumumab PK | DARA1, DARA2, DARA_CD38 | 3 |
| Dexamethasone PK | DEX_gut, DEX1 | 2 |
| Venetoclax PK | VEN_gut, VEN1 | 2 |
| Zoledronic Acid PK | ZOL1, ZOL2, ZOL_bone | 3 |
| **합계** | | **28** |

---

## 치료 시나리오 (Treatment Scenarios)

| 시나리오 | 레지멘 | 임상시험 근거 |
|---------|--------|------------|
| 1 | 무치료 (자연경과) | — |
| 2 | **VRd** (Bortezomib + Lenalidomide + Dex) | SWOG S0777 (Durie 2017) |
| 3 | **DRd** (Daratumumab + Lenalidomide + Dex) | MAIA (Facon 2019) |
| 4 | **KRd** (Carfilzomib + Lenalidomide + Dex) | ASPIRE (Stewart 2015) |
| 5 | **VenDex** (Venetoclax + Dex, t(11;14)) | BELLINI (Kumar 2020) |
| 6 | **DVRd** (Dara + VRd, 고위험군) | PERSEUS (Sonneveld 2024) |

---

## 약물 PK 파라미터 요약 (Drug PK Parameters)

| 약물 | 투여 경로 | CL (L/hr) | Vd (L) | 반감기 | 특이사항 |
|------|---------|-----------|--------|--------|--------|
| Bortezomib (BTZ) | IV/SC | 9.0 | 4.7 (central) | ~76h (deep) | 3-구획, SC F=83% |
| Carfilzomib (CFZ) | IV | 245 | 8.0 | ~1h | 매우 빠른 소실, 비가역적 |
| Lenalidomide (LEN) | PO | 3.2 | 65 | ~3h | F=90%, 신배설 |
| Daratumumab (DARA) | IV/SC | 0.2 (linear) | 3.0 (70kg) | ~14d | TMDD 모델 |
| Dexamethasone (DEX) | PO | 20 | 130 | ~4h | F=78% |
| Venetoclax (VEN) | PO | 14 | 256 | ~12h | F=55% (고지방식사) |
| Zoledronic acid (ZOL) | IV | 3.8 | 4.0 | ~105h (bone) | 골 결합 후 완속 방출 |

---

## 반응 기준 (Response Criteria — IMWG 2016)

| 반응 | M-protein 기준 | 기타 기준 |
|------|-------------|--------|
| **sCR** | 혈청/소변 M-protein 음성 | 정상 FLC ratio + BM <5% PC + NGS MRD 음성 |
| **CR** | 혈청/소변 M-protein 음성 | BM <5% PC |
| **VGPR** | ≥90% M-protein 감소 | 소변 <100 mg/24hr |
| **PR** | ≥50% M-protein 감소 | 소변 ≥90% 감소 (or <200 mg/24hr) |
| **SD** | PR 기준 미충족, PD 기준 미도달 | — |
| **PD** | ≥25% M-protein 증가 (최저치 대비) | 새 병변 또는 기존 병변 25% 증가 |

---

## Shiny 앱 탭 구성 (Shiny Dashboard Tabs)

| 탭 | 내용 |
|----|------|
| **환자 프로파일** | 체중·나이·ISS 병기·세포유전학적 위험·기저치 바이오마커 입력 |
| **PK 프로파일** | 레지멘별 약물 농도-시간 곡선, EC50 참조선, 약물 효과 시각화 |
| **PD / 종양 부담** | MM 세포 동역학 (감수성/내성 클론), M-protein 반응 |
| **임상 엔드포인트** | 반응 분류, 골 리모델링, 헤모글로빈 추이, 반응 테이블 |
| **시나리오 비교** | 여러 레지멘 동시 비교, 최적 반응 waterfall 도표 |
| **바이오마커** | 골 전환 마커, 사이토카인, β₂M vs M-protein 상관관계 |

---

## 주요 임상 시험 데이터 기반 파라미터 보정 (Parameter Calibration)

| 파라미터 | 값 | 근거 |
|---------|------|------|
| MM 세포 배가 시간 | ~35일 (kgrow=0.02/day) | Turesson 2004 |
| M-protein 반감기 | ~9일 (kMP_elim=0.08) | 임상 관찰 |
| VRd ORR (PR이상) | ~82% | SWOG S0777 |
| DRd PFS 개선 | 중앙값 61.9개월 vs 34.4개월 | MAIA trial |
| DVRd MRD 음성률 | ~60% (10⁻⁵) | PERSEUS trial |
| ZOL 골 재흡수 억제 | ~60-70% | 임상 NTX 감소 |

---

*생성일: 2026-06-21 | Claude Code Routine (CCR) 자동 생성*
