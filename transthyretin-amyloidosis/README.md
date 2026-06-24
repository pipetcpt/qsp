# 트랜스티레틴 아밀로이드증 (Transthyretin Amyloidosis, ATTR) — QSP 모델

[![DOT nodes](https://img.shields.io/badge/기계론적%20지도-120+%20노드-blue)](attr_qsp_model.svg)
[![ODE compartments](https://img.shields.io/badge/ODE%20구획-20개-green)](attr_mrgsolve_model.R)
[![Shiny tabs](https://img.shields.io/badge/Shiny%20탭-8개-orange)](attr_shiny_app.R)
[![References](https://img.shields.io/badge/참고문헌-45편-purple)](attr_references.md)

---

## 질환 개요

**트랜스티레틴(TTR) 아밀로이드증**은 간에서 합성되는 TTR 단백질의 사량체가 불안정해져 해리되고, 잘못 폴딩된 단량체가 올리고머 → 피브릴 → 아밀로이드로 응집하여 심장·말초신경 등 장기에 침착하는 진행성 전신 단백질 미스폴딩 질환입니다.

| 특성 | 내용 |
|------|------|
| **유병률** | ATTRwt: 75세 이상 HFpEF 환자의 ~13% / ATTRv: 유전성, V122I는 아프리카계 3.4% |
| **발병 기전** | TTR 사량체 해리 → 단량체 미스폴딩 → 올리고머(독성) → 아밀로이드 피브릴 → 장기 침착 |
| **주요 표현형** | ① ATTRwt 심근병증 (심부전, 전도 이상) ② ATTRv 다발신경병증 (V30M 등) |
| **진단 방법** | Tc-99m DPD/PYP 신티그래피 ± 유전자검사 (비생검 진단 가능) |
| **치료제** | 타파미디스(사량체 안정화), 파티시란·뷔트리시란·이노테르센(siRNA/ASO) |
| **예후** | ATTRwt 심장형 중앙 생존기간 2.5-5년 (무치료), 타파미디스로 30% 사망률↓ |

---

## 기계론적 지도 미리보기

[![ATTR QSP 기계론적 지도](attr_qsp_model.png)](attr_qsp_model.svg)

*(클릭하면 확대 가능한 SVG 파일로 이동)*

---

## 모델 구성: 10개 서브그래프 클러스터

| 클러스터 | 설명 | 주요 노드 |
|---------|------|---------|
| 1. TTR 유전학·간 생산 | 야생형/변이형 유전자 → mRNA → 사량체 분비 | TTR_gene, TTR_mRNA, Choroid_plexus, siRNA 표적 |
| 2. TTR 응집 연쇄 | 사량체 해리 → 단량체 → 올리고머 → 피브릴 | Tet_dissociation, TTR_misfolded, TTR_oligomers, TTR_amyloid |
| 3. 심장 ATTR 병리 | 아밀로이드 침착 → LV비후 → HFpEF → BNP/TnT↑ | Cardiac_dep, Diastolic_dysfxn, BNP_biomarker, AFib_ATTR |
| 4. 신경 ATTR 병리 | 말초신경 침착 → 소/대섬유 신경병증 → NIS↑ | PNS_dep, Axonal_degen, NIS_score_n, mNIS7_score |
| 5. 다장기 침범 | 신장·위장관·간·연조직 | Renal_dep, GI_submucosal, Soft_tissue |
| 6. 타파미디스 PK/PD | 경구 흡수 → T4 결합부위 점유 → 사량체 안정화 | TAF_plasma_Cp, TAF_TTR_complex, Dissoc_inhibit |
| 7. siRNA/ASO 치료 | LNP/GalNAc 전달 → RISC 로딩 → TTR mRNA 분해 | PAT_RISC_load, VUT_TTR_87, siRNA_net_TTR |
| 8. 임상 병기·결과 | NAC 병기분류, 생존율, 입원, QoL | NAC_staging, Stage1-3, Overall_surv |
| 9. 세포 스트레스 | 올리고머 독성 → Ca²⁺↑ → NLRP3 → 염증 | Oligomer_tox_s, NLRP3_inflamm, ER_stress_s |
| 10. 진단 워크업 | DPD 신티그래피·심초음파·CMR·생검 | DPD_PYP_scan, Perugini_grade, ECV_elevated |

---

## mrgsolve ODE 모델 (20개 구획)

### 구획 목록

| 그룹 | 구획 | 설명 |
|------|------|------|
| TTR PK | TTR_C, TTR_P | 혈장/말초 TTR 농도 (2구획 모델) |
| 응집 | TTR_MF, TTR_OL | 잘못폴딩 단량체, 가용성 올리고머 |
| 아밀로이드 | AMY_H, AMY_N, AMY_GI | 심장·신경·위장관 아밀로이드 침착 |
| 심장 | LV_THICK, BNP_C, TROP_C, CARD_FUNC | LV벽두께·NT-proBNP·TnT·심기능지수 |
| 신경 | NIS_TOT, AUTO_NP | NIS 총점, 자율신경 장애 지수 |
| 약물 PK | TAF_GUT, TAF_C, PAT_EFF, VUT_EFF, INO_EFF | 타파미디스 PK, siRNA 효과 구획 |
| 복합 결과 | SixMWT, KCCQ_IDX | 6분보행거리, KCCQ-OS |

### 핵심 방정식

| 모듈 | 방정식 |
|------|--------|
| TTR 생산-소실 | `dTTR_C/dt = kprod × (1 − Inh_siRNA) − kel × TTR_C − kdiss_eff × TTR_C` |
| 사량체 해리 | `kdiss_eff = kdiss_base × mut_factor × (1 − Emax_TAF × [TAF] / (EC50_TAF + [TAF]))` |
| 응집 연쇄 | `dTTR_MF/dt = kconf × kdiss_eff × TTR_C − kolig × TTR_MF` |
| 아밀로이드 침착 | `dAMY_H/dt = frac_heart × kfib × TTR_OL − kclear × AMY_H` |
| BNP 동역학 | `dBNP/dt = 0.1 × (BNP_base + BNP_scale × AMY_H + 50 × TTR_OL − BNP)` |

### 치료 시나리오 (6개)

| # | 시나리오 | 설계 | 임상시험 참조 |
|---|---------|------|--------------|
| 1 | 자연 경과 (무치료) | ATTRwt 심근병증, 무치료 5년 | Ruberg 2019 JACC |
| 2 | 타파미디스 61mg QD | 1일 1회 경구, Emax PD | ATTR-ACT (Maurer 2018 NEJM) |
| 3 | 파티시란 0.3mg/kg Q3W IV | LNP siRNA, TTR↓80% | APOLLO (Adams 2018 NEJM) |
| 4 | 뷔트리시란 25mg Q3M SC | GalNAc-siRNA, TTR↓87% | HELIOS-A (Solomon 2022 NEJM) |
| 5 | 이노테르센 284mg QW SC | 2'-MOE ASO, TTR↓75% | NEURO-TTR (Benson 2018 NEJM) |
| 6 | 타파미디스 + HF 적극 관리 | 병용 전략 | 임상 실무 모델 |

---

## Shiny 대시보드 (8개 탭)

| 탭 | 기능 |
|----|------|
| 1. 환자 프로파일 | 나이·성별·진단유형·NAC 병기 추정, ATTR 역학 요약 |
| 2. TTR PK/응집 | TTR 농도·응집 경로·타파미디스 억제율 실시간 시뮬레이션 |
| 3. 심장 ATTR | LV벽두께·NT-proBNP·트로포닌·심장기능 시뮬레이션, ATTR-ACT 참조 |
| 4. 신경 ATTR | NIS 점수·자율신경 지표, FAP 임상시험 비교 테이블 |
| 5. 시나리오 비교 | 6개 치료 시나리오 동시 비교 (지표 선택 가능), 말기 요약 표 |
| 6. 바이오마커 추적 | 6개 바이오마커 실시간 패널, NAC 병기 변화 |
| 7. 민감도 분석 | 타파미디스 EC50/Emax 등 파라미터 민감도, 토네이도 차트 |
| 8. 가상 환자 집단 | N=10-200명 몬테카를로, 스파게티·VPC·NAC 분포 |

---

## 실행 방법

```r
# 1. 기계론적 지도 SVG/PNG 렌더링
# dot -Tsvg attr_qsp_model.dot -o attr_qsp_model.svg
# dot -Tpng -Gdpi=150 attr_qsp_model.dot -o attr_qsp_model.png

# 2. mrgsolve 모델 실행 (R 4.0+)
install.packages(c("mrgsolve", "dplyr", "ggplot2", "patchwork", "tidyr"))
source("attr_mrgsolve_model.R")
results <- run_all_scenarios(mod_attr)
plot_attr_results(results)

# 3. Shiny 대시보드 실행
install.packages(c("shiny", "shinydashboard", "plotly", "DT", "shinyWidgets"))
shiny::runApp("attr_shiny_app.R")
```

---

## 파일 목록

| 파일 | 설명 |
|------|------|
| [`attr_qsp_model.dot`](attr_qsp_model.dot) | Graphviz 기계론적 지도 소스 (120+ 노드, 10 클러스터) |
| [`attr_qsp_model.svg`](attr_qsp_model.svg) | 렌더링된 SVG 벡터 이미지 |
| [`attr_qsp_model.png`](attr_qsp_model.png) | 렌더링된 PNG 이미지 (150 dpi) |
| [`attr_mrgsolve_model.R`](attr_mrgsolve_model.R) | mrgsolve ODE 모델 (20구획, 6시나리오) |
| [`attr_shiny_app.R`](attr_shiny_app.R) | Shiny 대시보드 (8탭) |
| [`attr_references.md`](attr_references.md) | 참고문헌 (45편 PubMed, 9섹션) |
| [`README.md`](README.md) | 본 문서 |

---

> *생성일: 2026-06-24 · Claude Code Routine (CCR)*
