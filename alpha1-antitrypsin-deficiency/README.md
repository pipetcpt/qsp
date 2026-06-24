# Alpha-1 Antitrypsin Deficiency (AATD) — QSP Model

[![Disease](https://img.shields.io/badge/Disease-AATD%20%2F%20Alpha--1%20Antitrypsin%20Deficiency-red)]()
[![Category](https://img.shields.io/badge/Category-Rare%20Genetic%20%2F%20Pulmonary%20%2F%20Hepatic-orange)]()
[![Nodes](https://img.shields.io/badge/Mechanistic%20Map-130%2B%20nodes%20%7C%2010%20clusters-blue)]()
[![ODE](https://img.shields.io/badge/ODE%20Model-20%20compartments%20%7C%206%20scenarios-green)]()
[![Shiny](https://img.shields.io/badge/Shiny%20App-6%20tabs-purple)]()
[![Refs](https://img.shields.io/badge/References-54%20PubMed-brightgreen)]()

---

## 질환 개요 (Disease Overview)

**알파-1 항트립신 결핍증(Alpha-1 Antitrypsin Deficiency, AATD)**은 *SERPINA1* 유전자 돌연변이로 인한 유전성 대사 질환으로, 간과 폐를 동시에 침범하는 독특한 이중 병태생리를 가집니다. 가장 중증 표현형인 **ZZ 유전자형(PIZZZ)**은 전 세계 약 1/3,000–5,000명에서 나타납니다.

### 핵심 병태생리 이중 축 (Dual Pathology)

| 경로 | 기전 | 결과 |
|------|------|------|
| **간 (기능획득 독성)** | Glu342Lys 돌연변이 → Z-AAT 단백질 소포체 내 루프-시트 중합체 축적 → ER 스트레스 → 간세포 손상/섬유화 | 간경변, 간세포암(HCC 위험 10-40배↑) |
| **폐 (기능소실)** | 혈청 AAT 부족(<11 µM 보호 역치) → 폐포 중성구 엘라스타제(NE) 무제한 활성 → 탄력소 분해 → 범소엽성 폐기종 | FEV1 급격 감소(정상의 2-3배), 중증 COPD |

---

## 모델 구성 (Model Components)

| 산출물 | 파일 | 사양 |
|--------|------|------|
| 🗺️ 기계론적 지도 | [`aatd_qsp_model.dot`](aatd_qsp_model.dot) / [`.svg`](aatd_qsp_model.svg) / [`.png`](aatd_qsp_model.png) | **130+ 노드, 10 클러스터** |
| ⚙️ mrgsolve ODE 모델 | [`aatd_mrgsolve_model.R`](aatd_mrgsolve_model.R) | **20 구획, 6 치료 시나리오** |
| 📊 Shiny 대시보드 | [`aatd_shiny_app.R`](aatd_shiny_app.R) | **6탭 인터랙티브 앱** |
| 📚 참고문헌 | [`aatd_references.md`](aatd_references.md) | **54개 PubMed 인용, 13섹션** |

---

## 기계론적 지도 (Mechanistic Map)

[![AATD QSP Mechanistic Map](aatd_qsp_model.png)](aatd_qsp_model.svg)

### 10개 서브그래프 클러스터

| # | 클러스터 | 핵심 요소 |
|---|---------|----------|
| 1 | 유전학 & 유전자형 | SERPINA1, M/Z/S/Null 대립유전자, MM/MZ/SZ/ZZ/SS 유전자형 |
| 2 | ER 단백질 품질 관리 | 소포체 접힘 (칼넥신, BiP/GRP78), Z-AAT 중합체, UPR (IRE1α/PERK/ATF6), ERAD, 자가포식 |
| 3 | 간 병리 | 간세포 손상, ER 스트레스, 간성상세포(HSC) 활성화, TGF-β1, 콜라겐 침착, Metavir F0-F4 섬유화, 간경변, HCC |
| 4 | AAT 생물학 & 분포 | M-AAT/Z-AAT 분비, 혈청/ELF 분포, NE/PR3/카텝신G 억제, 항염증 기능 |
| 5 | 프로테아제-항프로테아제 균형 | 폐 PMN 모집, NE 방출/활성화, MPO/ROS, MMP-12, MMP-9, SLPI, 엘라핀, TIMP-3, AAT 산화 |
| 6 | 폐 병리 | 탄력소/피브로넥틴/콜라겐/라미닌 분해, ECM 파괴, 범소엽성 폐기종, FEV1, FVC, DLCO, RV/TLC, 기관지확장증 |
| 7 | 폐 염증 cascade | IL-8, LTB4, TNF-α, IL-1β, IL-17, IFN-γ, IL-10, NF-κB, AP-1, 보체 C5a, CD8/Th17/NK 세포 |
| 8 | 약물 PK/PD | 보충요법(Prolastin-C/Zemaira/Aralast/Glassia), siRNA(Fazirsiran/Belcesiran), 교정제(VX-864/GSK3117391), NE 억제제(Alvelestat/Lonodelestat/Brensocatib), 유전자요법(rAAV/CRISPR), LABA/LAMA/ICS |
| 9 | 치료 PD 효과 | 혈청 AAT↑, ELF AAT↑, 폐 NE↓, Z-AAT 중합체↓, FEV1 감소 억제, CT 폐밀도 보존, 악화 감소 |
| 10 | 임상 엔드포인트 | SGRQ, CAT, mMRC, 6MWD, 혈청 AAT, ELF AAT, CT 폐기종 지수, 간 바이오마커, 사망률 |

---

## mrgsolve ODE 모델

### 구획 구성 (20개)

| 그룹 | 구획 | 설명 |
|------|------|------|
| **간 Z-AAT** | `ZAAT_ER`, `ZAAT_Poly`, `HSC_act`, `Liver_coll`, `Liver_fib` | ER내 Z-AAT, 중합체, 간성상세포, 콜라겐, 섬유화 |
| **혈청 AAT** | `AAT_C1`, `AAT_C2` | 혈청 AAT 2-구획 PK (V1=3.76 L/kg, t½=4.5일) |
| **폐 경로** | `PMN_lung`, `IL8_lung`, `NE_free`, `MMP12_lung`, `Elastin`, `FEV1_pct` | PMN, IL-8, 자유 NE, MMP-12, 탄력소, FEV1 |
| **약물 PK** | `AUG_C1`, `AUG_C2`, `NEi_A`, `NEi_C`, `siRNA_Eff`, `Gene_Eff` | 보충요법 2-구획, NE억제제, siRNA 효과 구획, 유전자요법 |

### 6개 치료 시나리오

| 시나리오 | 치료 | 근거 임상시험 |
|---------|------|-------------|
| **S1** 무치료 | — | 자연 경과 (Tanash 2010, Janus 1985) |
| **S2** AAT 보충요법 | Prolastin-C 60 mg/kg IV 주1회 | RAPID Trial (Chapman 2015, *Lancet*) |
| **S3** Fazirsiran (siRNA) | GalNAc-siRNA 200 mg SQ q12주 | ARO-AAT Phase 2 (Strnad 2022, *NEJM*) |
| **S4** NE 억제제 | Alvelestat 60 mg BID PO | Phase 2 (McElvaney 2020, *AJRCCM*) |
| **S5** 유전자요법 | rAAV-SERPINA1 단회 투여 | Flotte 2004, Mueller 2008 |
| **S6** 보충요법 + NE 억제제 | 조합 치료 | 시뮬레이션 기반 예측 |

### 주요 파라미터 보정 근거

| 파라미터 | 값 | 출처 |
|---------|-----|------|
| ZZ 혈청 AAT | 6-8 mg/dL (~2-7 µM) | Crystal 1990 |
| 보호 역치 | >11 µM (57 mg/dL) | RAPID Trial 2015 |
| ELF/혈청 비율 | 0.10 | Hubbard 1991 |
| AAT 반감기 | 4.5일 | Prolastin-C PK label |
| FEV1 연간 감소 (ZZ) | 50-200 mL/yr | Dirksen 1999/2009 |
| siRNA mRNA 억제 | ~88% | Strnad 2022 |
| Alvelestat NE 억제 | ~75% | McElvaney 2020 |

---

## Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 유전자형, 흡연 상태, FEV1, CT 폐기종 지수, 위험도 요약, 약물 표적 |
| 2. 약물 PK / AAT 수준 | 혈청 AAT 시간 경과, ELF AAT, ELF NE 억제율(%), PK 파라미터 테이블 |
| 3. 폐 PD | 자유 NE 활성, 탄력소 함량(%), MMP-12, 폐 PMN 부하 |
| 4. 임상 엔드포인트 | FEV1 (% 예측), CT 폐기종 지수, SGRQ 삶의 질, 연간 악화 위험 |
| 5. 시나리오 비교 | 6개 치료군 FEV1/혈청 AAT/Z-중합체 비교 그래프, 5년 결과 요약 표 |
| 6. 바이오마커 | 간 Z-AAT 중합체, Metavir 섬유화, 간성상세포 활성화, 데스모신(탄력소 분해지표) |

---

## 실행 방법 (Quick Start)

```r
# 1. 기계론적 지도 렌더링
# dot -Tsvg aatd_qsp_model.dot -o aatd_qsp_model.svg
# dot -Tpng -Gdpi=150 aatd_qsp_model.dot -o aatd_qsp_model.png

# 2. mrgsolve 시뮬레이션
install.packages(c("mrgsolve", "dplyr", "ggplot2", "patchwork"))
source("aatd_mrgsolve_model.R")   # 자동으로 6개 시나리오 실행 및 시각화

# 3. Shiny 대시보드
install.packages(c("shiny", "shinydashboard", "plotly", "DT", "shinycssloaders"))
shiny::runApp("aatd_shiny_app.R")
```

---

## 주요 임상 특징 요약

- **유병률**: ZZ 유전자형 약 1/3,000-5,000 (유럽계 > 아시아계)
- **진단 지연**: 평균 5-7년 (희귀질환 특성)
- **폐 표현형**: ZZ 성인의 80%, 하엽 우세 범소엽성 폐기종
- **간 표현형**: 신생아 담즙 정체, 성인 간경변/HCC
- **흡연 효과**: 흡연 시 FEV1 감소 2-3배 가속
- **ERS/ATS 지침**: FEV1 <80% + ZZ 확인 시 AAT 보충요법 권고

---

## 참고 문헌

54개 PubMed 인용 — [aatd_references.md](aatd_references.md) 참조  
주요 임상시험: RAPID Trial (2015, *Lancet*), ARO-AAT Phase 2 (2022, *NEJM*), Alvelestat Phase 2 (2020, *AJRCCM*)

---

*모델 생성: Claude Code Routine (CCR) | 2026-06-24*
