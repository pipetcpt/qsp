# 대장암 (Colorectal Cancer) — QSP Model

> **분류**: 종양·소화기 | **약어**: CRC | **버전**: 1.0 | **날짜**: 2026-06-23

---

## 질환 개요

대장암(Colorectal Cancer, CRC)은 결장 또는 직장에서 발생하는 악성 종양으로, 전 세계에서 세 번째로 흔하고 두 번째로 사망률이 높은 암입니다 (GLOBOCAN 2020: 연간 약 193만 신규 환자, 93만 명 사망). 대부분의 CRC는 정상 상피 → 선종(adenoma) → 침윤성 암으로 진행되는 **Vogelstein 염기서열**(APC→KRAS→SMAD4→TP53)을 따르며, 약 15%는 MSI-H 경로를 통해 빠르게 진행됩니다.

---

## 병태생리 주요 경로

| 경로 | 핵심 분자 | 임상 의미 |
|------|-----------|-----------|
| **Wnt/β-catenin** | APC(mut), β-catenin, TCF/LEF | APC 돌연변이 ~80%, c-Myc/Cyclin D1 과발현 |
| **RAS/MAPK** | KRAS(40%), NRAS(5%), BRAF V600E(10%) | 항-EGFR 치료 반응 예측 |
| **PI3K/AKT** | PIK3CA(20%), PTEN 소실 | mTOR 활성화, 내성 기전 |
| **TP53** | R175H/R248W/R273H hotspot | 세포주기 정지, 아포토시스 장애 |
| **TGF-β/EMT** | SMAD4(35%), ZEB1/TWIST | 전이, 불량 예후 |
| **MSI/MMR** | MLH1/MSH2/MSH6/PMS2 결손 | 높은 종양변이부담→면역항암제 감수성 |

---

## 분자 아형 (Consensus Molecular Subtypes, CMS)

| CMS | 주요 특성 | 치료 전략 |
|-----|-----------|-----------|
| **CMS1** (MSI Immune) | MSI-H, BRAF V600E, 높은 면역활성 | 항-PD-1/PD-L1 (pembrolizumab) |
| **CMS2** (Canonical) | CIN 높음, WNT/MYC 활성화 | FOLFOX + anti-EGFR (RAS-WT) |
| **CMS3** (Metabolic) | KRAS 돌연변이, 지질대사 이상 | FOLFIRI/FOLFOX ± BEV |
| **CMS4** (Mesenchymal) | TGF-β 활성, EMT, 간질 침윤 | 항혈관신생 ± 항-EGFR |

---

## 모델 구성 (Model Architecture)

### 기계론적 지도 (Mechanistic Map)

[![CRC QSP Map](crc_qsp_model.png)](crc_qsp_model.svg)

- **노드 수**: 130+ (12개 서브그래프 클러스터)
- **클러스터**: Wnt · RAS/MAPK · PI3K/AKT · TP53/Apoptosis · Cell Cycle · TGF-β/EMT · TME · Immune Checkpoints · Angiogenesis · MSI/MMR · Drug PK/PD · Clinical Endpoints

### mrgsolve ODE 모델

**20 ODE 구획:**

| 구획 그룹 | 변수 | 모델링 내용 |
|-----------|------|-------------|
| 5-FU PK | FU1, FU2, FU_ic | 2-구획 + 세포내 활성화, DPD 대사 |
| Oxaliplatin | OX1, OX_DNA | 백금-DNA 부가물 형성/복구 |
| Irinotecan | IRI1, SN38, SN38G | UGT1A1*28 다형성 반영 |
| Bevacizumab | BEV1, BEV_VEGF | TMDD 모델, VEGF 중화 |
| Cetuximab | CTX1 | EGFR 점유율 |
| Pembrolizumab | PEM1 | PD-1 점유율 |
| Tumor | Ts, Tr | 민감세포/내성세포 로지스틱 성장 |
| Biomarkers | CEA, ctDNA | 종양 부담 연동 |
| Immune | CD8eff | T세포 활성도 (PD-1 차단) |
| VEGF, EGFR_occ, PD1_occ | 3 | 표적 분자 동역학 |

**치료 시나리오 (7개):**

| # | 시나리오 | 근거 임상시험 |
|---|---------|--------------|
| 1 | 무치료 (자연경과) | — |
| 2 | FOLFOX6 | MOSAIC (André 2004) |
| 3 | FOLFIRI | GERCOR (Tournigand 2004) |
| 4 | FOLFOX + Bevacizumab | NO16966 (Saltz 2008) |
| 5 | FOLFIRI + Cetuximab (RAS-WT) | CRYSTAL (Van Cutsem 2009) |
| 6 | FOLFIRI + Bevacizumab | TRIBE (Falcone 2013) |
| 7 | Pembrolizumab (MSI-H) | KEYNOTE-177 (André 2020) |

### Shiny 앱 구조 (6 탭)

| 탭 | 내용 |
|----|------|
| 환자 프로파일 | 나이·BSA·RAS/BRAF/MSI 상태·CMS 아형·치료 선택 |
| PK 프로파일 | 5-FU/옥살리플라틴/이리노테칸/베바시주맙/세툭시맙/펨브롤리주맙 농도-시간 |
| PD 주요지표 | 종양 직경 RECIST 변화, 민감/내성 세포, CD8 T세포 활성도 |
| 임상 엔드포인트 | 최우수 반응, PFS 추정, RECIST 분류 |
| 시나리오 비교 | 다중 요법 동시 비교 그래프 및 요약 표 |
| 바이오마커 | CEA, ctDNA, 유리 VEGF, EGFR/PD-1 점유율 |

---

## 약물 요약

| 약물 | 기전 | 적응증 |
|------|------|--------|
| 5-FU/LV | TS 억제 → dTTP 고갈 | 모든 CRC 1차 |
| 옥살리플라틴 | 백금-DNA 부가물 | FOLFOX 백본 |
| 이리노테칸 | Top1 독 (SN-38) | FOLFIRI 백본 |
| 베바시주맙 | 항-VEGF-A (TMDD) | 전이성 CRC 1차 |
| 세툭시맙 | 항-EGFR (RAS-WT) | CRYSTAL/OPUS |
| 펨브롤리주맙 | 항-PD-1 | MSI-H/dMMR CRC 1차 |
| 레고라페닙 | 다중-키나제 TKI | 3차 이상 |

---

## 참고 파일

| 파일 | 설명 |
|------|------|
| [`crc_qsp_model.dot`](crc_qsp_model.dot) | Graphviz 기계론적 지도 원본 |
| [`crc_qsp_model.svg`](crc_qsp_model.svg) | 벡터 그래픽 (확대 가능) |
| [`crc_qsp_model.png`](crc_qsp_model.png) | 150 dpi 래스터 이미지 |
| [`crc_mrgsolve_model.R`](crc_mrgsolve_model.R) | mrgsolve ODE 모델 (20구획, 7시나리오) |
| [`crc_shiny_app.R`](crc_shiny_app.R) | Shiny 대화형 대시보드 (6탭) |
| [`crc_references.md`](crc_references.md) | 참고문헌 40개 (섹션별 분류) |
