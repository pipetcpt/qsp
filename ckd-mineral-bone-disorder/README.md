# CKD-Mineral Bone Disorder (CKD-MBD) QSP Model

## 개요 (Overview)

만성신장병-미네랄골질환(CKD-MBD)은 사구체여과율(GFR) 감소로 인한 FGF23-Klotho-PTH-비타민D 축의 조절 이상, 골이영양증, 혈관 석회화를 포함하는 복잡한 전신 증후군입니다. 이 모델은 해당 병태생리의 핵심 기전을 17개 ODE 구획으로 표현하며, 7가지 치료 시나리오를 포함합니다.

Chronic Kidney Disease–Mineral Bone Disorder (CKD-MBD) is a systemic syndrome characterized by dysregulation of the FGF23–Klotho–PTH–vitamin D axis, renal osteodystrophy, and vascular calcification as GFR declines. This QSP model captures these mechanisms in 17 ODE compartments with 7 treatment scenarios.

---

## 파일 목록 (Files)

| 파일 | 설명 |
|------|------|
| `ckdmbd_qsp_model.dot` | Graphviz 기계론적 지도 소스 (115+ nodes, 10 clusters) |
| `ckdmbd_qsp_model.svg` | SVG 기계론적 지도 (벡터) |
| `ckdmbd_qsp_model.png` | PNG 기계론적 지도 (150 dpi) |
| `ckdmbd_mrgsolve_model.R` | mrgsolve ODE 모델 (17구획, 7시나리오) |
| `ckdmbd_shiny_app.R` | Shiny 인터랙티브 대시보드 (7탭) |
| `ckdmbd_references.md` | 참고문헌 42편 (9섹션) |

---

## 기계론적 지도 (Mechanistic Map)

[![CKD-MBD QSP Mechanistic Map](ckdmbd_qsp_model.png)](ckdmbd_qsp_model.svg)

---

## 모델 구조 (Model Architecture)

### 10개 서브그래프 클러스터

| 클러스터 | 주요 노드 |
|---------|-----------|
| ① Kidney Function | GFR, CKD_Stage, NaPi-IIa/IIc, CYP27B1, CYP24A1, Klotho |
| ② Phosphate & FGF23 | Pi_plasma, FGF23_plasma, FGF23_bone, PHEX, DMP1, GALNT3 |
| ③ Vitamin D Axis | VitD3_skin, 25-OH-D, 1,25-OH₂D, VDR, RXR, VDRE, CaBP9k |
| ④ Calcium Homeostasis | Ca_plasma, Ca_ionized, CaSR_kidney, TRPV5, TRPV6, Calcitonin |
| ⑤ Parathyroid Gland | PTH_plasma, CaSR_ptg, PTG_mass, SecHPT, TerHPT, Nodular |
| ⑥ Bone Remodeling | Osteoblast, Osteoclast, RANKL, OPG, BMD, Sclerostin, DKK1 |
| ⑦ Vascular Calcification | VSMC, MGP, ucMGP, Fetuin-A, Hydroxyapatite, CarotidIMT |
| ⑧ Drug PK | Sevelamer, Cinacalcet, Paricalcitol, Etelcalcetide, Denosumab |
| ⑨ Drug PD | CaSR_act, VDR_act, Pi_bind_GI, RANKL_inh, Emax models |
| ⑩ Clinical Outcomes | iPTH_lab, Pi_lab, Ca_lab, FGF23_lab, DXA_BMD, CV_event |

---

## ODE 구획 (17 ODEs)

| # | 구획 | 설명 |
|---|------|------|
| 1 | Pi | 혈청 인산염 (mg/dL) |
| 2 | FGF23 | 혈장 FGF23 (pg/mL) |
| 3 | Klotho | 가용성 Klotho (상대 단위) |
| 4 | PTH | 완전형 PTH (pg/mL) |
| 5 | VitD25 | 25-OH-비타민D (nmol/L) |
| 6 | VitD_act | 1,25-OH₂D 칼시트리올 (pg/mL) |
| 7 | Ca | 혈청 칼슘 (mg/dL) |
| 8 | OB | 골아세포 활성 |
| 9 | OC | 파골세포 활성 |
| 10 | BMD | 골밀도 (정상 대비) |
| 11 | VC | 혈관 석회화 점수 |
| 12 | CIN_GUT | 시나칼셋 장관 (mg) |
| 13 | CIN_PLASMA | 시나칼셋 혈장 (ng/mL) |
| 14 | PAR_PLASMA | 파리칼시톨 혈장 (ng/mL) |
| 15 | ETEL_PLASMA | 에텔칼세타이드 혈장 (ng/mL) |
| 16 | DEN_DEPOT | 데노수맙 SC 저장소 (mg) |
| 17 | DEN_PLASMA | 데노수맙 혈장 (mg/L) |

---

## 치료 시나리오 (Treatment Scenarios)

| 시나리오 | 치료 | 대상 |
|---------|------|------|
| S1 | 무치료 (CKD G5 자연 경과) | 기준선 |
| S2 | Sevelamer 2400 mg/day | 인산염 결합제 |
| S3 | Cinacalcet 60 mg/day | 칼시미메틱 |
| S4 | Paricalcitol 4 mcg 3x/주 IV | 비타민D 수용체 작용제 |
| S5 | Sevelamer + Cinacalcet | 병합 |
| S6 | Etelcalcetide 5 mg 3x/주 IV | 차세대 칼시미메틱 |
| S7 | Sevelamer + Etelcalcetide + Denosumab | 삼중 치료 |

---

## KDIGO 2017 치료 목표 (Treatment Targets)

| 지표 | 목표 범위 |
|------|----------|
| 완전형 PTH (iPTH) | 150–600 pg/mL (CKD G5D) |
| 혈청 인산염 | < 5.5 mg/dL |
| 혈청 칼슘 | 8.4–10.2 mg/dL |
| Ca×Pi 곱 | < 55 mg²/dL² |
| 25-OH-비타민D | > 75 nmol/L |

---

## Shiny 앱 탭 구성 (Shiny App Tabs)

1. **환자 프로파일** — CKD 병기, 기저 검사 수치, 시뮬레이션 기간
2. **약물 PK** — 인산염 결합제 · 칼시미메틱 · 비타민D · 데노수맙 PK 프로파일
3. **PTH & 미네랄** — iPTH, Pi, Ca, Ca×Pi, FGF23, 비타민D, KDIGO 달성 여부
4. **골 질환** — 골밀도, 골아/파골세포, P1NP/CTX, 골절 위험
5. **심혈관** — 혈관 석회화, GFR 추이, Klotho, CV 위험 평가
6. **시나리오 비교** — 5개 표준 치료 전략 비교 (iPTH, Pi, BMD, VC)
7. **바이오마커** — 정규화 대시보드, 스파이더 차트, FGF23-Klotho 축

---

## 주요 참고문헌 (Key References)

- KDIGO CKD-MBD Guidelines 2017 · Gutierrez 2008 (FGF23-mortality, NEJM) · Hu 2011 (Klotho-VC, JASN)
- Block 2004 (mineral mortality, JASN) · Tentori 2008 (DOPPS) · Chertow 2002 (cinacalcet ACHIEVE)
- Luo 1997 (MGP knockout calcification, Nature) · Ketteler 2003 (fetuin-A, Lancet)
- Peterson & Riggs 2010 (calcium-bone QSP model, Bone)

전체 참고문헌: [`ckdmbd_references.md`](ckdmbd_references.md)
