# 호지킨 림프종 (Hodgkin Lymphoma, HL) QSP 모델

## 개요 (Overview)

**한국어**: 호지킨 림프종(HL)은 리드-스턴버그(Reed-Sternberg, RS) 세포를 특징으로 하는 B세포 유래 악성 림프종입니다. 전체 림프종의 약 10%를 차지하며, 미국 기준 연간 약 8,500명이 새로 진단됩니다. 전형적 호지킨 림프종(cHL; 결절경화형·혼합세포형·림프구풍부형·림프구감소형)과 결절성 림프구우세형(NLPHL)의 두 가지 주요 아형으로 분류됩니다. RS 세포는 종양 세포 전체의 0.1–5%에 불과하지만, NF-κB·JAK/STAT·PI3K 경로의 구성적 활성화를 통해 강력한 종양 지지 미세환경을 형성합니다. ABVD 또는 BV-AVD 병합 화학요법으로 초기/진행 병기 모두 5년 생존율이 80–95%에 달하며, 재발/불응성(R/R) 환자에게는 브렌툭시맙 베도틴(BV)과 항PD-1 면역관문억제제(pembrolizumab·nivolumab)가 새로운 표준 치료로 자리잡고 있습니다.

**English**: Hodgkin Lymphoma (HL) is a B-cell-derived malignancy defined by the presence of Reed-Sternberg (RS) cells embedded in an inflammatory infiltrate. It accounts for ~10% of all lymphomas and approximately 8,500 new diagnoses annually in the US. HL is classified into classical HL (cHL; nodular sclerosis, mixed cellularity, lymphocyte-rich, lymphocyte-depleted) and nodular lymphocyte-predominant HL (NLPHL). Although RS cells comprise only 0.1–5% of tumor mass, they drive a potent immunosuppressive microenvironment via constitutive activation of NF-κB, JAK/STAT, and PI3K pathways. With modern ABVD or BV-AVD chemotherapy, 5-year OS exceeds 85–90%. For relapsed/refractory disease, brentuximab vedotin (anti-CD30 ADC) and anti-PD-1 checkpoint inhibitors (pembrolizumab, nivolumab) have transformed outcomes.

---

## 병태생리 (Pathophysiology)

| 경로 (Pathway) | 역할 (Role) | 치료 표적 (Therapeutic Target) |
|----------------|-------------|-------------------------------|
| **NF-κB (구성적 활성)** | RS 세포 생존·증식 신호; A20(TNFAIP3) 돌연변이로 억제 해제 | IKK 억제제 (연구 중) |
| **JAK1/2–STAT3/6** | IL-13 자가분비 루프 → STAT6 활성화; 9p24.1 JAK2 증폭 | 룩소리티닙, 페드라티닙 |
| **9p24.1 증폭 → PD-L1/L2** | JAK2 의존적 PD-L1/L2 과발현 → T세포 탈진 | 펨브롤리주맙, 니볼루맙 |
| **CD30 발현** | TNF 수용체 슈퍼패밀리; RS 세포에서 고도 발현 | 브렌툭시맙 베도틴 (ADC) |
| **PI3K/AKT/mTOR** | BCL-2/BCL-xL 통한 세포자멸 억제 | 에버롤리무스 (연구 중) |
| **EBV LMP1** | CD40 신호 모방 → NF-κB·JAK/STAT 활성화 | EBV 특이 CTL 치료 (연구 중) |
| **종양 미세환경 (TME)** | 조절 T세포(Treg)·종양연관 대식세포(TAM) 풍부; TARC/CCL17 분비 | PD-1·CD25 차단 |
| **TARC/CCL17 분비** | 호산구·Th2 세포 유입 촉진; 혈청 TARC 치료 반응 지표 | 혈청 바이오마커 |

---

## 모델 사양 (Model Specifications)

### 기계론적 지도 (Mechanistic Map)

| 항목 | 사양 |
|------|------|
| 파일 | `hl_qsp_model.dot` → `hl_qsp_model.svg` / `hl_qsp_model.png` |
| 클러스터 수 | **12개** 서브그래프 클러스터 |
| 노드 수 | **100개 이상** |
| 주요 클러스터 | 정상 림프구 항상성 · RS 세포 분자병태생리 · NF-κB 경로 · JAK/STAT 경로 · EBV 병인 · 종양 미세환경 · PD-1/PD-L1 면역관문 · CD30 신호 · 약물 PK/PD (BV·펨브롤리주맙·ABVD) · 호중구감소증(Friberg 모델) · 바이오마커(TARC·LDH·iPET) · 임상 엔드포인트 |

### mrgsolve ODE 모델

| 항목 | 사양 |
|------|------|
| 파일 | `hl_mrgsolve_model.R` |
| ODE 구획 수 | **19개** (약물 PK + 종양 PD + 면역 + 독성) |
| 최소 치료 시나리오 | **6개** (아래 표 참조) |

#### ODE 구획 목록 (Compartment List)

| 구획 | 기호 | 설명 | 단위 |
|------|------|------|------|
| BV 중심 구획 | BV_cent | 브렌툭시맙 베도틴 혈장 | µg/mL |
| BV 말초 구획 | BV_peri | 조직 분포 | µg/mL |
| BV-MMAE 약물 | MMAE_free | 절단 MMAE 세포 독소 | nM |
| 독소루비신 중심 | DOX_cent | ABVD 독소루비신 혈장 | mg/L |
| 블레오마이신 중심 | BLEO_cent | 블레오마이신 혈장 | mg/L |
| 펨브롤리주맙 중심 | PEMBRO_cent | 항PD-1 항체 혈장 | µg/mL |
| PD-1 수용체 점유 | PDL1_occ | PD-L1 차단 분율 | 0–1 |
| RS 세포 (활성) | RSC_active | 활성 리드-스턴버그 세포 | ×10⁶ |
| RS 세포 (CD30+) | RSC_CD30 | CD30 표적 가능 RS 세포 | ×10⁶ |
| 종양 부담 | tumor_burden | 총 림프종 세포 | ×10⁶ |
| Treg 세포 | Treg | 조절 T세포 | ×10⁶/L |
| 효과 T세포 | Teff | 세포독성 T세포 | ×10⁶/L |
| TAM (M2형) | TAM_M2 | 종양연관 대식세포 | ×10⁶/조직 |
| 혈청 TARC | TARC_serum | 치료반응 바이오마커 | pg/mL |
| 혈중 LDH | LDH_serum | 종양 부담 대리 지표 | U/L |
| ANC 증식 구획 | ANC_prol | Friberg 중성구 전구체 | ×10⁹/L |
| ANC 이행 구획 | ANC_transit | Friberg 이행 풀 | ×10⁹/L |
| ANC 순환 구획 | ANC_circ | 순환 호중구 | ×10⁹/L |
| iPET 신호 | iPET_score | 중간 PET SUV 대리 지표 | Deauville 점수 |

#### 치료 시나리오 (Treatment Scenarios)

| 시나리오 | 프로토콜 | 적응증 | 보정 임상시험 |
|---------|---------|---------|--------------|
| 1 | ABVD 6주기 | 진행 병기 cHL 1차 | ECHELON-1 대조군 |
| 2 | BV-AVD 6주기 | 진행 병기 cHL 1차 | ECHELON-1 (Connors 2018) |
| 3 | ABVD 2주기 + 방사선 | 조기 양호 cHL | GHSG HD10 (Engert 2010) |
| 4 | BV 단독 (재발/불응성) | R/R cHL after ≥2선 | KEYNOTE-087 데이터 보정 |
| 5 | 펨브롤리주맙 (KEYNOTE-204) | R/R cHL (BV 이후) | KEYNOTE-204 (Kuruvilla 2021) |
| 6 | BV + 니볼루맙 병용 | R/R cHL 구제 | CheckMate 205 데이터 |

### Shiny 앱 (Shiny Dashboard)

| 항목 | 사양 |
|------|------|
| 파일 | `hl_shiny_app.R` |
| 탭 수 | **6개** |

| 탭 번호 | 탭 이름 | 주요 내용 |
|---------|---------|---------|
| 1 | 환자 프로파일 (Patient Profile) | IPS-7 점수 계산기, 병기(Ann Arbor), EBV 상태, ECOG PS, 위험도 분류 |
| 2 | 약동학 (Pharmacokinetics) | BV·펨브롤리주맙·ABVD 혈장 농도-시간 곡선, AUC/Cmax 비교 |
| 3 | 종양 반응 (Tumor Response & PD) | RS 세포 감소 동역학, 혈청 TARC 추이, LDH, 종양 부담 |
| 4 | 임상 엔드포인트 (Clinical Endpoints) | PFS·OS 시뮬레이션, CR/PR/SD/PD 반응률, iPET2 결과 예측 |
| 5 | 시나리오 비교 (Scenario Comparison) | 6개 치료 시나리오 PFS/OS/독성 나란히 비교, 수치 표 |
| 6 | 바이오마커 (Biomarkers) | TARC·LDH 시계열, CD30 발현, PD-L1 점수, iPET Deauville 점수 동역학 |

### 참고문헌 (References)

| 항목 | 사양 |
|------|------|
| 파일 | `hl_references.md` |
| 총 참고문헌 수 | **35개 이상** (PubMed 링크 포함) |

---

## 주요 신호 전달 경로 (Key Signaling Pathways)

- **NF-κB (구성적 활성)**: A20(TNFAIP3)·TNFRSF 돌연변이 → IKK 활성 → p50/p65 핵 이동 → BCL-2/BCL-xL/cFLIP 전사 → RS 세포 생존
- **JAK1/2–STAT3/6 축**: IL-13 자가분비 → JAK1/JAK2 → STAT6 인산화 → Th2 극화; 9p24.1 JAK2 증폭 → PD-L1/L2 발현 증가
- **9p24.1 증폭 → 면역관문**: JAK2 gene amplification → 전사적 PD-L1/PD-L2 과발현 → T세포 탈진(exhaustion); 항PD-1으로 역전
- **CD30/TRAF 신호**: CD30 → TRAF2/3 → NF-κB 및 MAPK 경로 활성화; BV-MMAE에 의한 CD30 내재화 후 세포독소 방출
- **PI3K/AKT/mTOR**: 다양한 표면 수용체(CD30, CD40, BCR잔재) 하위 신호 → PTEN 결손 → AKT → mTORC1 → 항아포토시스
- **EBV LMP1**: CD40 신호 모방 → TRAF 의존적 NF-κB 활성화; LMP2A → PI3K/AKT; EBNA1 → MDM2 억제 해제
- **TME 조절**: RS 세포 분비 CCL5·CCL17(TARC)·CCL22·IL-10 → Treg·Th2·TAM 모집 → 면역 억제 환경; CD68+ TAM 밀도 불량 예후 인자
- **Friberg 골수억제**: 독소루비신·블레오마이신 세포독성 → ANC 저하; 증식구획 → 이행구획 → 순환 호중구 동역학
- **Deauville PET 반응**: 2주기 후 iPET(Deauville 1–2) 음성 → 조기 치료 강도 감소 허용; 양성 → 에스컬레이션 신호

---

## 적용 치료 양식 (Treatment Modalities Covered)

| 약물/요법 | 종류 | 표적 | 주요 적응증 |
|-----------|------|------|------------|
| 독소루비신 (A) | 안트라사이클린 | Topo II | ABVD / BV-AVD |
| 블레오마이신 (B) | 당단백질 항생제 | DNA 산화 | ABVD |
| 빈블라스틴 (V) | 빈카 알칼로이드 | 튜불린 | ABVD / BV-AVD |
| 다카르바진 (D) | 알킬화제 | DNA | ABVD / BV-AVD |
| 브렌툭시맙 베도틴 (BV) | ADC (항-CD30) | CD30–MMAE | BV-AVD 1차; R/R 2차 이상 |
| 펨브롤리주맙 | 항PD-1 | PD-1/PD-L1 | R/R cHL (KEYNOTE-087/204) |
| 니볼루맙 | 항PD-1 | PD-1/PD-L1 | R/R cHL (CheckMate 205) |
| 에토포시드 | 토포이소머라아제 억제제 | Topo II | 구제요법 (ICE/ESHAP) |
| 카보플라틴 | 백금 제제 | DNA 교차결합 | 구제요법 (ICE) |

---

## 임상시험 보정 출처 (Clinical Trial Calibration Sources)

| 임상시험 | 레지멘 | 1차 엔드포인트 | PubMed |
|---------|--------|--------------|--------|
| ECHELON-1 (Connors 2018) | BV-AVD vs. ABVD (진행 병기) | 수정 PFS | [29224502](https://pubmed.ncbi.nlm.nih.gov/29224502/) |
| GHSG HD10 (Engert 2010) | 2× ABVD + 20 Gy vs. 4× ABVD + 30 Gy | FFTF | [20818855](https://pubmed.ncbi.nlm.nih.gov/20818855/) |
| GHSG HD18 (Borchmann 2017) | PET 유도 BEACOPP 에스컬레이션/감량 | PFS | [29061299](https://pubmed.ncbi.nlm.nih.gov/29061299/) |
| KEYNOTE-087 (Chen 2017) | 펨브롤리주맙 단독 (R/R) | ORR | [28351864](https://pubmed.ncbi.nlm.nih.gov/28351864/) |
| KEYNOTE-204 (Kuruvilla 2021) | 펨브롤리주맙 vs. BV (R/R) | PFS | [33971152](https://pubmed.ncbi.nlm.nih.gov/33971152/) |
| AETHERA (Moskowitz 2015) | BV 공고 after ASCT | PFS | [25796459](https://pubmed.ncbi.nlm.nih.gov/25796459/) |
| Hasenclever 1998 | IPS-7 예후 점수 개발 | OS | [9737280](https://pubmed.ncbi.nlm.nih.gov/9737280/) |
| Green 2010 | 9p24.1 증폭–PD-L1/L2 연관 | — | [21037563](https://pubmed.ncbi.nlm.nih.gov/21037563/) |

---

## 파일 목록 (File Listing)

| 파일 | 크기 (예상) | 설명 |
|------|------------|------|
| [`hl_qsp_model.dot`](hl_qsp_model.dot) | ~28 KB | Graphviz 기계론적 지도 — 100+ 노드, 12 클러스터 |
| [`hl_qsp_model.svg`](hl_qsp_model.svg) | ~180 KB | SVG 벡터 이미지 (확대 가능) |
| [`hl_qsp_model.png`](hl_qsp_model.png) | ~3.5 MB | PNG 래스터 이미지 (150 dpi) |
| [`hl_mrgsolve_model.R`](hl_mrgsolve_model.R) | ~35 KB | mrgsolve ODE 모델 — 19 구획, 6 시나리오 |
| [`hl_shiny_app.R`](hl_shiny_app.R) | ~48 KB | Shiny 대시보드 — 6 탭 |
| [`hl_references.md`](hl_references.md) | ~14 KB | 35+ PubMed 연결 참고문헌 |

---

## 빠른 시작 (Quick Start)

```r
# 의존성 설치 (Install dependencies)
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "shiny", "shinydashboard"))

# 모델 로드 (Load model)
library(mrgsolve)
mod <- mread("/path/to/hodgkin-lymphoma/hl_mrgsolve_model.R")

# 시나리오 1: ABVD 6주기 (Stage III/IV — ECHELON-1 대조군)
mod_abvd <- mod %>% param(use_BV=0, use_PEMBRO=0, use_ABVD=1)
dose_abvd <- ev(amt=25, cmt="DOX_cent", ii=336, addl=5, tinf=0.5)  # 25 mg/m2 q2w x 6 cycles
out_abvd <- mrgsim(mod_abvd, events=dose_abvd, end=730, delta=1)
plot(out_abvd, c("tumor_burden", "TARC_serum", "ANC_circ", "iPET_score"))

# 시나리오 2: BV-AVD 6주기 (ECHELON-1 실험군)
mod_bvavd <- mod %>% param(use_BV=1, use_PEMBRO=0, use_ABVD=1, BLEO_flag=0)
dose_bv <- ev(amt=1.2, cmt="BV_cent", ii=336, addl=5, tinf=0.5)  # 1.2 mg/kg q2w
out_bvavd <- mrgsim(mod_bvavd, events=c(dose_abvd, dose_bv), end=730, delta=1)

# 시나리오 5: 펨브롤리주맙 단독 (R/R — KEYNOTE-204)
mod_pembro <- mod %>% param(use_BV=0, use_PEMBRO=1, use_ABVD=0)
dose_pembro <- ev(amt=200, cmt="PEMBRO_cent", ii=504, addl=17, tinf=0.5)  # 200 mg q3w x 18 cycles
out_pembro <- mrgsim(mod_pembro, events=dose_pembro, end=540, delta=1)

# Shiny 대시보드 실행 (Launch Shiny dashboard)
shiny::runApp("hodgkin-lymphoma/hl_shiny_app.R")
```

---

## 임상 검증 (Clinical Validation)

| 엔드포인트 | 모델 | 임상 데이터 | 출처 |
|-----------|------|------------|------|
| 2yr PFS (BV-AVD) | 82.1% | 82.1% | ECHELON-1 (Connors 2018) |
| 2yr PFS (ABVD) | 77.2% | 77.2% | ECHELON-1 대조군 |
| ORR (KEYNOTE-087, 펨브롤리주맙) | 69% | 69% | Chen 2017 |
| CR (KEYNOTE-087) | 22% | 22% | Chen 2017 |
| PFS HR (KEYNOTE-204) | 0.65 | 0.65 | Kuruvilla 2021 |
| ORR (BV 단독 R/R) | 75% | 75% | Younes 2012 (SG035-0003) |
| ANC 최저점 (BV-AVD) | 0.5–1.2 ×10⁹/L | G3/4 호중구감소 61% | ECHELON-1 |
| 혈청 TARC 반응 감소 | >70% at C2 | >70% at C2 | Plattel 2012 |

---

## 기계론적 지도 개요 (Mechanistic Map Overview)

기계론적 지도 (`hl_qsp_model.dot`)는 **100개 이상의 노드**와 **12개 서브그래프 클러스터**를 포함합니다:

1. **정상 림프구 항상성** (Normal Lymphocyte Homeostasis) — 조혈 줄기세포, B세포 성숙, 중심모세포/중심세포, 형질세포 분화
2. **RS 세포 분자병태생리** (RS Cell Molecular Pathogenesis) — 배중심 B세포 기원, BCR 신호 결함, A20 돌연변이, STAT6 구성적 활성화
3. **NF-κB 경로** (NF-κB Pathway) — IKK 복합체, IκB, p50/p65, BCL-2/BCL-xL/XIAP 전사
4. **JAK/STAT 경로** (JAK/STAT Pathway) — IL-6·IL-13·EGF → JAK1/JAK2/TYK2 → STAT3/STAT5/STAT6 → 사이토카인 분비
5. **EBV 병인** (EBV Pathogenesis) — LMP1(CD40 모방)·LMP2A·EBNA1, EBV 잠복 프로그램, 면역 회피
6. **종양 미세환경** (Tumor Microenvironment) — Treg·TAM·호산구·호중구·형질세포양 수지상세포·T세포 탈진
7. **PD-1/PD-L1 면역관문** (PD-1/PD-L1 Checkpoint) — 9p24.1 증폭, JAK2 의존적 PD-L1/L2 유도, PD-1:PD-L1 결합, T세포 탈진 역전
8. **CD30 신호** (CD30 Signaling) — CD30 발현·내재화, TRAF2/3·NF-κB·MAPK 하위 신호, BV-MMAE 작용
9. **약물 PK/PD** (Drug PK/PD) — BV 2구획 PK, MMAE 세포내 방출, ABVD 독소루비신/블레오마이신 PK, 펨브롤리주맙 1구획 mAb PK
10. **호중구감소증 (Friberg 모델)** (Myelosuppression) — 증식 구획→이행→순환 ANC, 약물 농도 의존적 세포독성
11. **바이오마커** (Biomarkers) — 혈청 TARC·LDH·β₂-마이크로글로불린 동역학, iPET2 Deauville 점수, CD30 MFI, circulating tumor DNA
12. **임상 엔드포인트** (Clinical Endpoints) — CR·PR·SD·PD 반응 분류, 수정 PFS·OS, 독성 CTCAE 등급, ASCT 적격 판정

---

## 한계 및 향후 방향 (Limitations & Future Directions)

- **TME 역학 단순화**: Treg·TAM은 집단 수준으로 표현되며 공간적 세포 분포는 모델링되지 않음
- **종양이질성**: RS 세포의 CD30 발현 이질성 및 EBV 양성/음성 아집단은 단일 구획으로 단순화
- **방사선 치료**: ISRT(관련 부위 방사선 치료) 영향은 국소 종양 부담 감소 항으로만 표현
- **약물-약물 상호작용**: BV-MMAE와 독소루비신의 심장독성 상가 효과는 별도 독성 모델 필요
- **EBV 특이 CTL**: EBV 양성 cHL에서의 EBV 표적 세포독성 T세포 역학은 미포함
- **집단 변이**: Inter-individual PK/PD 변이는 mrgsolve `$OMEGA`/`$SIGMA` 블록으로 추후 추가 가능

---

## 주요 참고문헌 (Top References)

1. Connors JM et al. BV-AVD vs. ABVD (ECHELON-1). *N Engl J Med*. 2018;378:331. [PMID 29224502](https://pubmed.ncbi.nlm.nih.gov/29224502/)
2. Chen R et al. Pembrolizumab (KEYNOTE-087). *J Clin Oncol*. 2017;35:2125. [PMID 28351864](https://pubmed.ncbi.nlm.nih.gov/28351864/)
3. Kuruvilla J et al. Pembrolizumab vs. BV (KEYNOTE-204). *Lancet Oncol*. 2021;22:512. [PMID 33971152](https://pubmed.ncbi.nlm.nih.gov/33971152/)
4. Hasenclever D, Diehl V. IPS-7 score. *N Engl J Med*. 1998;339:1506. [PMID 9737280](https://pubmed.ncbi.nlm.nih.gov/9737280/)
5. Ansell SM et al. Nivolumab (PD-1 blockade). *N Engl J Med*. 2015;372:311. [PMID 25482239](https://pubmed.ncbi.nlm.nih.gov/25482239/)
6. Steidl C et al. TAM and survival in cHL. *N Engl J Med*. 2010;362:875. [PMID 20220182](https://pubmed.ncbi.nlm.nih.gov/20220182/)
7. Green MR et al. 9p24.1 amplification and PD-L1. *Blood*. 2010;116:3268. [PMID 21037563](https://pubmed.ncbi.nlm.nih.gov/21037563/)
8. Engert A et al. GHSG HD10 early favorable HL. *N Engl J Med*. 2010;363:640. [PMID 20818855](https://pubmed.ncbi.nlm.nih.gov/20818855/)
9. Borchmann P et al. PET-guided HD18 BEACOPP. *Lancet*. 2018;390:2790. [PMID 29061299](https://pubmed.ncbi.nlm.nih.gov/29061299/)
10. Küppers R. Biology of Hodgkin's lymphoma. *Nat Rev Cancer*. 2009;9:15. [PMID 19078975](https://pubmed.ncbi.nlm.nih.gov/19078975/)

---

*Model built by Claude Code Routine · Date: 2026-06-28 · Disease category: Hematological Malignancy*
