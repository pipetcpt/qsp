# ALS (Amyotrophic Lateral Sclerosis) — QSP Model

[![ALS Mechanistic Map](als_qsp_model.png)](als_qsp_model.svg)

---

## 개요 (Overview)

**근위축성 측삭경화증 (ALS)**은 상위 및 하위 운동신경원 모두에 영향을 미치는 치명적인 신경퇴행성 질환입니다. 증상 발현 후 중앙값 생존기간은 2–5년이며, 호흡부전이 주요 사망 원인입니다. 전 세계 유병률은 약 4–6/100,000명이며, 90%는 산발성, 10%는 가족성(SOD1, C9orf72, TDP-43/FUS 돌연변이)입니다.

**Amyotrophic lateral sclerosis (ALS)** is a fatal progressive neurodegenerative disease affecting both upper and lower motor neurons. Median survival from symptom onset is 2–5 years, with respiratory failure as the leading cause of death (~65%). Global prevalence is ~4–6 per 100,000. ~90% are sporadic; ~10% familial (SOD1 2%, C9orf72 5–10%, FUS/TDP-43 others).

---

## 모델 구성 (Model Components)

### 기계론적 지도 (Mechanistic Map)

**196개 이상 노드 · 10개 서브그래프 클러스터**

| 클러스터 | 구성 요소 수 | 내용 |
|---------|------------|------|
| 1. Genetic & Molecular | 18 | SOD1, TDP-43, FUS, C9orf72, UBQLN2, VCP, NEK1 등 |
| 2. Protein Quality Control | 18 | 26S Proteasome, Autophagy, ER Stress (PERK/IRE1/ATF6) |
| 3. Oxidative Stress & Mitochondria | 20 | ROS, GSH, Nrf2, 미토콘드리아 기능, 아팝토시스 |
| 4. Glutamate Excitotoxicity | 18 | EAAT2, AMPA/NMDA 수용체, Ca²⁺, nNOS, ONOO⁻ |
| 5. Motor Neuron Survival | 23 | BDNF/TrkB, IGF-1, GDNF, PI3K/Akt/mTOR, ERK1/2 |
| 6. Neuroinflammation | 25 | 미세아교세포 M1/M2, TNF-α, IL-1β, NLRP3, Treg |
| 7. Axonal Transport | 18 | NF-L/M/H, Kinesin, Dynein, 미토콘드리아 수송 |
| 8. NMJ & Motor Unit | 18 | ACh, nAChR, Agrin/MuSK, EMG, MUNIX |
| 9. Drug PK/PD | 21 | Riluzole, Edaravone, Tofersen, AMX0035, Masitinib |
| 10. Clinical Endpoints | 17 | ALSFRS-R, FVC%, NfL, pNF-H, 생존 |

### mrgsolve ODE 모델 (ALS ODE Model)

**26 ODEs · 60+ 파라미터 · 7 치료 시나리오**

```r
# 주요 구획 목록
# PK : DEPOT_RIL, C1_RIL, C2_RIL, IV_EDA, DEPOT_TOF, C1_TOF, C2_TOF,
#      DEPOT_PB, C_PB  (9개)
# Disease: MN_upper, MN_lower, SOD1_wt, SOD1_mis, TDP43_nuc, TDP43_cyto,
#          Glut_syn, Ca_i, ROS, GSH, Mito, Mic_act, TNFa, BDNF  (14개)
# Clinical: NfL_CSF, ALSFRS, FVC  (3개)
```

#### 핵심 발병기전 수식

```
# 운동신경원 사망률 (다인자 모델)
MN_death_rate = k_MN_death × (1 + 3·SOD1_burden + 2·TDP43_cyto_frac
                              + 1.5·ROS_norm + 1.2·Glu_excess + 0.8·TNFa
                              - BDNF_prot_wt × (BDNF/BDNF₀))
```

### Shiny 대시보드 탭 (6 Tabs)

| 탭 | 내용 |
|---|------|
| Patient Profile | 환자 파라미터, 질환 개요, PK 요약표, 임상 milestone |
| Drug PK | Riluzole/Edaravone/Tofersen 농도-시간 곡선 및 약물 효과 |
| Biomarkers | NfL, TNF-α, ROS/GSH, 미토콘드리아 기능, BDNF |
| Clinical Endpoints | ALSFRS-R, FVC%, 운동신경원 생존%, 글루탐산 독성 |
| Scenario Comparison | 7가지 치료군 동시 비교 (ALSFRS-R · MN 생존 · NfL) |
| Mechanistic Pathways | TDP-43 핵/세포질 동역학, SOD1 단백질 misfolding |

---

## 약물 PK/PD 파라미터 (Drug PK/PD Parameters)

| 약물 | 용법 | F(%) | t½ | CL | Vd | 작용기전 | IC50/EC50 |
|-----|-----|------|----|----|----|---------|----|
| **Riluzole** | 50 mg PO BID | 60 | 12 h | 28 L/h | 245 L | ↓ Glu release via Na⁺ channel | IC50=0.5 μg/mL |
| **Edaravone** | 60 mg IV/day (cycle) | ~100 | 4.5 h | 18 L/h | 120 L | Free-radical scavenger | IC50=1.2 μg/mL |
| **Tofersen** | 100 mg SC q4w | — | 7 days | 0.5 L/h | 20 L | SOD1 mRNA knockdown (ASO) | EC50=0.1 μg/mL CSF |
| **AMX0035** | PB 3g + TUDCA 1g PO BID | 85 | 3 h | 12 L/h | 50 L | ↓ ER stress / mitochondrial protection | EC50=50 μmol/L |
| **Masitinib** | 4.5 mg/kg/day PO | 58 | 40 h | 25 L/h | 190 L | c-Kit/PDGFR-β → microglia ↓ | Phase 3 ongoing |

---

## 7가지 치료 시나리오 (Treatment Scenarios)

| 시나리오 | 치료 | 예상 ALSFRS-R 슬로다운 |
|---------|------|-----------------|
| 1. Untreated | 없음 | 기준선 |
| 2. Riluzole | 50 mg BID | ~12% 개선 |
| 3. Edaravone | 60 mg/day (cycle) | ~10% 개선 |
| 4. Riluzole + Edaravone | 병용 | ~20% 개선 |
| 5. Tofersen | SC q4w (SOD1-ALS) | ~30-40% 개선 (SOD1 subtype) |
| 6. AMX0035 | PB+TUDCA BID | ~25% 개선 |
| 7. All Drugs | 모든 병용 | ~45-50% 개선 |

---

## 병태생리 요약 (Pathophysiology Summary)

```
유전적 소인 (SOD1/TDP-43/FUS/C9orf72)
    ↓
단백질 항상성 붕괴 → 독성 집적체 형성
    ↓
① 글루탐산 독성: EAAT2↓ → 시냅스 Glu↑ → AMPA/NMDA과활성 → Ca²⁺↑
② 산화 스트레스: SOD1 기능↓ → ROS↑ → GSH↓ → 미토콘드리아 손상
③ 신경염증: 미세아교세포 M1 편향 → TNF-α/IL-1β↑ → 신경독성
④ 영양인자 결핍: BDNF↓ / GDNF↓ → 생존신호 감소
⑤ 축삭수송 장애: NF 축적 → 미토콘드리아 수송 차단
    ↓
운동신경원 소실 → 탈신경위축 → ALSFRS-R 저하 → 호흡부전 → 사망
```

---

## 파일 목록 (File List)

| 파일 | 설명 |
|------|------|
| [als_qsp_model.dot](als_qsp_model.dot) | Graphviz 기계론적 지도 소스 (196+ 노드, 10 클러스터) |
| [als_qsp_model.svg](als_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [als_qsp_model.png](als_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [als_mrgsolve_model.R](als_mrgsolve_model.R) | mrgsolve ODE 모델 (26 구획, 60+ 파라미터, 7 시나리오) |
| [als_shiny_app.R](als_shiny_app.R) | Shiny 대시보드 (6탭, plotly, DT) |
| [als_references.md](als_references.md) | 참고문헌 40편 (PubMed 링크 포함) |

---

## 주요 임상시험 참고 (Key Clinical Trials)

| 시험 | 약물 | 결과 |
|-----|------|------|
| Bensimon 1994 (N Engl J Med) | Riluzole | 생존기간 +3개월 |
| Lacomblez 1996 (Lancet) | Riluzole | 확인 + 용량최적화 |
| Abe 2017 (Lancet Neurol) | Edaravone | ALSFRS-R 저하 +33% 개선 |
| ATLAS/Miller 2022 (NEJM) | Tofersen | SOD1 단백질 ↓50%, NfL ↓60% |
| CENTAUR/Paganoni 2020 (NEJM) | AMX0035 | ALSFRS-R +2.3점/yr (45 weeks) |
| AB Science Phase 3 | Masitinib | Phase 3 ALSFRS-R primary endpoint 진행 중 |

---

*모델 작성일: 2026-06-21 | ALS QSP Library v1.0*
