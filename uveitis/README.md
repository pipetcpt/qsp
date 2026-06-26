# Uveitis QSP Model

[![Disease](https://img.shields.io/badge/Disease-Uveitis-blue)]()
[![Category](https://img.shields.io/badge/Category-Autoimmune%20%7C%20Ocular-orange)]()
[![Compartments](https://img.shields.io/badge/ODE%20Compartments-20-green)]()
[![Scenarios](https://img.shields.io/badge/Treatment%20Scenarios-7-purple)]()
[![References](https://img.shields.io/badge/References-60-yellow)]()

## 개요 (Overview)

**포도막염(Uveitis)**은 안구의 포도막(홍채·모양체·맥락막)에 발생하는 염증 질환으로, 전방/중간/후방/범포도막염으로 분류됩니다. 비감염성 포도막염은 자가면역 기전에 의해 발생하며, 치료하지 않을 경우 낭포황반부종(CME), 이차 녹내장, 백내장, 망막 박리를 통해 영구적 시력 손실로 진행합니다.

**Uveitis** is inflammation of the uveal tract (iris, ciliary body, choroid) and is classified anatomically as anterior, intermediate, posterior, or panuveitis. Non-infectious uveitis is driven by autoimmune mechanisms involving Th1/Th17 cells, TNF-α, IL-6, and VEGF-mediated blood-ocular barrier (BOB) disruption. The primary site-specific PK challenge is the blood-aqueous barrier (BAB) and blood-retinal barrier (BRB).

---

## 기계론적 지도 (Mechanistic Map)

[![Mechanistic Map](uvt_qsp_model.png)](uvt_qsp_model.svg)

*클릭하면 고해상도 SVG 파일로 이동합니다.*

### 주요 클러스터 (12 Clusters)

| 클러스터 | 내용 |
|---------|------|
| Ocular Anatomy & Immune Privilege | 전방·후방 안구 구조, 면역특권 소실 기전 |
| Triggering Mechanisms | HLA-B27, IRBP/S-Ag 자가항원, 감염성 유발, VKH/Birdshot |
| Innate Immune Activation | DC, 대식세포, 중성구, NLRP3 인플라마솜, NF-κB |
| T Cell Adaptive Immunity | Th1/Th17, Treg, CD8+, TCR/CD28 공동자극, FoxP3 |
| B Cell & Humoral Immunity | 형질세포, 자가항체(anti-IRBP, anti-Recoverin), 면역복합체 |
| Cytokine & Chemokine Network | TNF-α, IL-6/STAT3, IL-17A, VEGF-A, COX-2/PGE2 |
| Blood-Ocular Barrier Disruption | BAB/BRB 파괴, 전방세포/플레어, CME, 신생혈관 |
| Clinical Manifestations | 시력(logMAR), IOP, SUN 등급, 백내장, 이차녹내장 |
| Drug PK | 점안제, 안와주위/유리체내 주사, 전신 경구/IV/SC |
| Drug PD Mechanisms | GR 점유, NF-κB 억제, TNF 중화, VEGF 차단, 칼시뉴린 억제 |
| Complications & Disease Course | 급성↔만성 전환, 재발률, 구조적 손상 |
| Biomarkers & Monitoring | SUN 등급, OCT-CST, FFA, 혈청 TNF/IL-6, 약물농도/ADA |

---

## ODE 모델 구성 (Model Structure)

### 구획 (20 Compartments)

| 모듈 | 구획 | 설명 |
|------|------|------|
| **Drug PK** | Cgut, Cp, Cperiph | 경구/전신 PK (1구획 흡수 + 2구획) |
| | C_ant, C_post | 전방방·후방 안구 약물 농도 |
| | C_depot | 유리체내 서방형 임플란트 (Ozurdex) |
| **Immune** | T_eff, T_reg | Th1/Th17 효과 T세포, Treg |
| | APC_act, Macro | 활성화 APC/DC, M1 대식세포 |
| **Cytokines** | TNF, IL6, IL17, VEGF | 4가지 핵심 사이토카인 |
| **Barrier** | BAB_int, BRB_int | 혈-안방장벽, 혈-망막장벽 무결성 (0-1) |
| **Clinical PD** | Cells_AH, CME | 전방세포, 낭포황반부종 |
| | VA_def, IOP_e | 시력 결손 (logMAR), IOP 과잉 |
| | GR_occ | 글루코코르티코이드 수용체 점유율 |

### 치료 시나리오 (7 Scenarios)

| ID | 치료법 | 용량/용법 | 주요 기전 | 근거 임상시험 |
|----|--------|----------|----------|--------------|
| S1 | 무치료 | — | 자연 경과 | — |
| S2 | 점안 프레드니솔론 1% | 매 6시간 | 국소 NF-κB 억제 | Cunningham 1990 |
| S3 | 안와주위 트리암시놀론 40mg | 단회 | 지속성 코르티코스테로이드 | Sallam 2011 |
| S4 | 유리체내 덱사메타손 임플란트 (Ozurdex) | 0.7mg 서방형 | 후방 BRB 회복, CME ↓ | HURON Trial (Lowder 2011) |
| S5 | 전신 프레드니손 1mg/kg/일 (점감) | 매일 → 유지 | 전신 GR 활성화, T세포 억제 | MUST Trial (Kempen 2011) |
| S6 | 아달리무맙 40mg SC q2w | 격주 피하주사 | TNF-α 중화, leukostasis ↓ | VISUAL I (Jaffe 2016 NEJM) |
| S7 | 병합 (프레드니손 + 아달리무맙) | S5 + S6 | 이중 기전, 빠른 장벽 회복 | VISUAL I/II |

---

## 주요 수식 (Key Equations)

### BAB 무결성
```
dBAB_int/dt = k_BABrep × (1 + 2·GR_occ) × (1 - BAB_int)
            - k_BABdeg × (TNF + 0.5·IL6) × BAB_int
```

### BRB 무결성
```
dBRB_int/dt = k_BRBrep × (1 + 1.5·GR_occ + 2·aVEGF_eff) × (1 - BRB_int)
            - k_BRBdeg × (VEGF + 0.3·TNF) × BRB_int
```

### 낭포황반부종 (CME)
```
dCME/dt = k_CMEform × (1-BRB_int) × VEGF
        - k_CMEres × (1 + 2·aVEGF_eff + GR_occ) × CME
```

### GR 점유 (Hill 함수)
```
GR_occ = Emax_cs × [Drug]^H / (EC50_cs^H + [Drug]^H)
```

---

## 임상 파라미터 보정 (Clinical Calibration)

| 파라미터 | 모델 값 | 임상 데이터 | 출처 |
|---------|--------|-----------|------|
| 아달리무맙 효과 onset | ~2-4주 | 4주 내 유의한 염증 감소 | VISUAL I |
| Ozurdex CME 해소율 | 60-70% | Day 60 완전 해소 65% | HURON Trial |
| 스테로이드 반응자 IOP 상승 | 30-40% | ~38% ≥5mmHg 상승 | Jaffe 2006 |
| Treg/Teff 비율 | 0.67 at flare | 건강 대조군의 약 1/3 | Caspi 2010 |
| BAB 파괴 → 전방세포 | ~48-72h | Flare onset 1-3일 | Nussenblatt 1985 |

---

## Shiny 앱 구조 (8 Tabs)

| 탭 | 내용 |
|----|------|
| 1. Patient Profile | 환자 설정, 질환 중증도, 치료 선택, 실행 |
| 2. Pharmacokinetics | 혈장/안구 약물 농도, GR 점유율 |
| 3. Immune & Cytokines | T세포 집단, APC/대식세포, TNF/IL-6/IL-17/VEGF |
| 4. Barrier Integrity | BAB/BRB 무결성, 전방세포(SUN 등급) |
| 5. Clinical Endpoints | 시력(logMAR), OCT-CST, IOP, 임상 마일스톤 |
| 6. Scenario Comparison | 7개 치료법 비교 (시각/CME/TNF/요약 표) |
| 7. Virtual Patients | VP 집단 시뮬레이션 (n=최대 500), 아형별 반응률 |
| 8. Biomarker Monitor | 사이토카인/OCT/약물 농도 추적, 모니터링 스케줄 |

---

## 파일 목록 (Files)

| 파일 | 설명 |
|------|------|
| [`uvt_qsp_model.dot`](uvt_qsp_model.dot) | Graphviz 기계론적 지도 (130+ 노드, 12 클러스터) |
| [`uvt_qsp_model.svg`](uvt_qsp_model.svg) | 벡터 형식 지도 (고해상도) |
| [`uvt_qsp_model.png`](uvt_qsp_model.png) | 래스터 형식 지도 (150 dpi) |
| [`uvt_mrgsolve_model.R`](uvt_mrgsolve_model.R) | mrgsolve ODE QSP 모델 (20구획, 7 시나리오, VP n=200) |
| [`uvt_shiny_app.R`](uvt_shiny_app.R) | Shiny 인터랙티브 대시보드 (8탭) |
| [`uvt_references.md`](uvt_references.md) | 참고문헌 60편 (15개 섹션) |
| [`README.md`](README.md) | 이 파일 |

---

## 실행 방법 (Usage)

```r
# 1. mrgsolve 모델 실행
install.packages(c("mrgsolve","dplyr","ggplot2","tidyr"))
source("uvt_mrgsolve_model.R")

# 2. Shiny 앱 실행
install.packages(c("shiny","shinydashboard","DT"))
shiny::runApp("uvt_shiny_app.R")
```

```bash
# Graphviz 지도 재렌더링
dot -Tsvg uvt_qsp_model.dot -o uvt_qsp_model.svg
dot -Tpng -Gdpi=150 uvt_qsp_model.dot -o uvt_qsp_model.png
```

---

## 주요 결론 (Key Findings)

1. **아달리무맙(S6)** 은 TNF 중화를 통해 BAB 회복 및 전방 염증 억제에 가장 효과적이며, VISUAL I/II 시험 결과를 재현합니다.
2. **유리체내 임플란트(S4)** 는 CME 해소에 탁월하나, 장기 스테로이드 IOP 상승 위험이 있습니다.
3. **병합 요법(S7)** 이 가장 빠른 장벽 회복을 보이며, 구조적 손상 예방에 우수합니다.
4. **가상 환자 분석** 에서 Th17 우세형은 항-TNF 단독 요법에 상대적으로 느린 반응을 보여, IL-17 차단제 추가 고려가 필요합니다.
5. **스테로이드 점안제(S2)** 는 전방 포도막염에 효과적이나 후방 병변에는 안구 내 농도가 불충분합니다.

---

*생성일: 2026-06-26 | QSP Disease Model Library | Claude Code Routine (CCR)*
