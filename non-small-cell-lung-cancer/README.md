# Non-Small Cell Lung Cancer (NSCLC) — QSP Model

> **비소세포 폐암(NSCLC)** 정량적 시스템 약리학(QSP) 모델  
> Driver mutations · Targeted therapies · Immune checkpoint inhibitors · Chemotherapy · Resistance mechanisms

[![NSCLC QSP Model](nsclc_qsp_model.png)](nsclc_qsp_model.svg)

---

## 질환 개요 (Disease Overview)

비소세포 폐암(Non-Small Cell Lung Cancer, NSCLC)은 폐암의 약 85%를 차지하며 전 세계적으로 암 사망 원인 1위입니다.

**핵심 병태생리 (Key Pathophysiology):**
- **Driver mutations:** EGFR(exon19 del/L858R/T790M), KRAS G12C, ALK/ROS1 fusion, BRAF V600E, MET exon14, RET, NTRK
- **Oncogenic signaling:** RAS-RAF-MEK-ERK, PI3K-AKT-mTOR, JAK-STAT3 → 무제한 증식
- **Immune evasion:** PD-L1 과발현 → CD8+ T세포 기능 억제; Treg·MDSC 축적; IDO1·TGF-β
- **Tumor microenvironment (TME):** HIF1α → VEGF → 종양혈관신생; EMT → 전이
- **Drug resistance:** T790M/C797S (EGFR), bypass signaling (MET amp, HER2 amp), SCLC transformation

---

## 모델 구성 (Model Structure)

### 기계론적 지도 (Mechanistic Map)
| 항목 | 값 |
|------|-----|
| 노드 수 | 231 |
| 클러스터 수 | 10 |
| 엣지 수 | 301 |
| 렌더 엔진 | sfdp |

**10개 클러스터:**
1. **Oncogene drivers** — EGFR/KRAS/ALK/ROS1/BRAF/MET/HER2/RET/NTRK/NRG1/PIK3CA/STK11/KEAP1/TP53
2. **Oncogenic signaling** — RAS/RAF/MEK/ERK, PI3K/AKT/mTOR, JAK/STAT3, SRC/FAK
3. **Tumor biology** — HIF1α-VEGF-angiogenesis, EMT, CDK4/6-RB1 cell cycle, BCL2/BAX apoptosis, TERT
4. **Immune checkpoint** — PD-L1/PD-1, CTLA4/CD28/B7, CD8+/NK/Treg/TAM, IDO1, LAG3/TIM3/TIGIT
5. **Drug PK/PD** — Osimertinib, Alectinib, Sotorasib, Pembrolizumab, Atezolizumab, Cisplatin, Pemetrexed
6. **Resistance mechanisms** — T790M, C797S, MET/HER2 amp, AXL/FGFR bypass, SCLC transformation
7. **Biomarkers** — PD-L1 TPS, TMB, ctDNA, ALK FISH, CEA, CYFRA 21-1, OS/PFS/ORR
8. **Tumor dynamics** — Sensitive/resistant cell pools, Gompertz growth, RECIST response categories
9. **Toxicity** — irAE (pneumonitis/colitis/nephritis), TKI rash/ILD, chemo myelosuppression
10. **Clinical endpoints** — OS, PFS, ORR, DCR, RECIST CR/PR/SD/PD, Stage IIIA-IVB, CNS/bone mets

### mrgsolve ODE 모델
| 항목 | 값 |
|------|-----|
| ODE 구획 수 | 19 |
| 약물 PK 구획 | 9 (Osimertinib 2 + Alectinib 2 + Sotorasib 2 + Pembrolizumab 1 + Cisplatin 1 + Pemetrexed 1) |
| 질환 PD 구획 | 10 (tumor sensitive/resistant/total, CD8+ Teff, PD1 occupancy, Treg, CEA, ctDNA, ANC, PD-L1) |
| 치료 시나리오 | 7 |
| 시뮬레이션 기간 | 24개월 |

**7가지 치료 시나리오:**
| 시나리오 | 적응증 | 근거 임상시험 |
|----------|--------|-------------|
| S1 | 무치료 (자연 진행) | — |
| S2 | Osimertinib 80mg/일 | FLAURA (PFS 18.9 mo) |
| S3 | Alectinib 600mg BID | ALEX (PFS 34.8 mo) |
| S4 | Carboplatin + Pemetrexed Q3W | ECOG 1594 |
| S5 | Pembrolizumab 200mg Q3W (PD-L1≥50%) | KEYNOTE-024 (PFS 10.3 mo) |
| S6 | Pembrolizumab + Carboplatin + Pemetrexed | KEYNOTE-189 (PFS 9.0 mo) |
| S7 | Sotorasib 960mg/일 (KRAS G12C+) | CodeBreaK100 (ORR 37%) |

### Shiny 대시보드 (6개 탭)
| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 나이/ECOG/Stage/조직형/분자 프로파일/PD-L1/TMB; 치료 적응증 매트릭스 |
| 2. 약동학(PK) | 약물별 다중 cycle 농도-시간 곡선; Cmax/AUC/t½ 표; 신·간기능 용량 조정 |
| 3. 종양 반응 | 감수성/저항성 세포 궤적; RECIST waterfall; spider plot |
| 4. 바이오마커 | CEA, ctDNA 대립유전자 분율, ANC nadir, PFS KM 곡선; 내성 예측 |
| 5. 시나리오 비교 | 7가지 치료 다중 비교; ORR/DCR/PFS 표; forest plot + 임상시험 benchmark |
| 6. 독성·안전성 | 혈액독성 시뮬레이션; irAE 확률; 장기별 신호등; 용량 수정 권고 |

---

## 파일 목록 (Files)

| 파일 | 설명 |
|------|------|
| [`nsclc_qsp_model.dot`](nsclc_qsp_model.dot) | Graphviz DOT 기계론적 지도 소스 |
| [`nsclc_qsp_model.svg`](nsclc_qsp_model.svg) | 벡터 그래픽 |
| [`nsclc_qsp_model.png`](nsclc_qsp_model.png) | 래스터 이미지 |
| [`nsclc_mrgsolve_model.R`](nsclc_mrgsolve_model.R) | mrgsolve ODE 모델 + 시뮬레이션 코드 |
| [`nsclc_shiny_app.R`](nsclc_shiny_app.R) | Shiny 인터랙티브 대시보드 |
| [`nsclc_references.md`](nsclc_references.md) | 참고문헌 78편 (섹션별 분류) |

---

## 실행 방법 (How to Run)

```bash
# 기계론적 지도 렌더링
sfdp -Tsvg nsclc_qsp_model.dot -o nsclc_qsp_model.svg
sfdp -Tpng -Gdpi=96 nsclc_qsp_model.dot -o nsclc_qsp_model.png
```

```r
# mrgsolve 시뮬레이션
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "survminer"))
library(mrgsolve)
source("nsclc_mrgsolve_model.R")

# Shiny 대시보드
install.packages(c("shiny", "shinydashboard", "ggplot2", "dplyr", "plotly", "DT"))
shiny::runApp("nsclc_shiny_app.R")
```

---

## 주요 파라미터 (Key Parameters)

| 파라미터 | 값 | 단위 | 근거 |
|----------|-----|------|------|
| 기저 종양 성장률 (kg) | 0.002 | /day | NSCLC TGI 문헌 |
| Osimertinib t½ | 48 | h | FLAURA PK |
| Alectinib t½ | 33 | h | ALEX PK |
| Pembrolizumab t½ | 27 | days | KEYNOTE biologic PK |
| 내성 출현률 (kr) | 0.0001 | /day | 임상 내성 발생 모델 |
| Osimertinib Emax | 0.92 | — | FLAURA 중앙 PFS 보정 |
| PD-L1 TPS 임계값 | 50% | % | KEYNOTE-024 선정 기준 |

---

## 참고문헌 핵심 (Key References)

- Mok TS et al. (2009) IPASS: gefitinib in EGFR+ NSCLC. *NEJM*. PMID: 19692680
- Soria JC et al. (2018) FLAURA: osimertinib vs 1st-gen TKI. *NEJM*. PMID: 29151359
- Peters S et al. (2017) ALEX: alectinib vs crizotinib. *NEJM*. PMID: 28586279
- Reck M et al. (2016) KEYNOTE-024: pembrolizumab vs chemo (PD-L1≥50%). *NEJM*. PMID: 27718847
- Gandhi L et al. (2018) KEYNOTE-189: pem+chemo. *NEJM*. PMID: 29658856
- Skoulidis F et al. (2021) CodeBreaK100: sotorasib. *NEJM*. PMID: 34096690

전체 78편 참고문헌 → [`nsclc_references.md`](nsclc_references.md) 참조

---

*Generated by Claude Code Routine (CCR) — 2026-06-23*
