# Epilepsy (뇌전증) QSP Model

**분류**: 신경학 | 만성 신경계 질환  
**날짜**: 2026-06-21  
**버전**: 1.0

---

## 개요 (Overview)

뇌전증(간질)은 전 세계 약 5천만 명에게 영향을 미치는 가장 흔한 만성 신경계 질환 중 하나로,
두 번 이상의 비유발성 발작 또는 한 번의 발작 이후 10년 재발 위험이 높은 경우로 정의됩니다 (ILAE 2014).
성인의 뇌전증 발생률은 약 50/100,000/년이며, 전 세계 유병률은 0.5–1%입니다.

이 QSP 모델은 뇌전증의 분자적·세포적·네트워크 수준 병태생리로부터 항뇌전증약물(AED)의
PK/PD 모델, 약물 내성 기전, 그리고 임상 결과(발작 빈도, SUDEP 위험, 삶의 질)까지를
정량적으로 연결하는 통합 시스템 모델입니다.

---

## 기계론적 지도 (Mechanistic Map)

[![Epilepsy QSP Model](epi_qsp_model.png)](epi_qsp_model.svg)

**10개 서브그래프 클러스터 (136개 노드):**

| 클러스터 | 내용 |
|---------|------|
| 1. 유전적 원인 | SCN1A/SCN2A/SCN8A, KCNQ2/3, GABRG2/A1/B3, DEPDC5, TSC1/2, CDKL5, STXBP1, PCDH19 |
| 2. 후천적 원인 | 내측두엽경화증(MTS), 외상성 뇌손상, 뇌졸중, 피질이형성증(FCD), 자가면역뇌염 |
| 3. 전압개폐 이온채널 | Nav1.1/1.2/1.6, Kv7.2/7.3, HCN1/2, T/N형 Ca2+ 채널, BK/SK 채널 |
| 4. GABA성 억제 시스템 | GAD65/67, GABA-A (αβγ 서브유닛), GABA-B, PV+/SST+ 사이신경세포 |
| 5. 글루타메이트성 흥분 시스템 | NMDA/AMPA/Kainate 수용체, mGluR1/2, EAAT1/2 재흡수 수송체 |
| 6. 발작 네트워크 역학 | 해마 CA3→CA1→DG, 시상-피질 루프, 발작역치(Θ), 뇌전증지속증(SE) |
| 7. AED 약동학 | VPA(2구획), LEV(1구획, 신장), CBZ(자가유도+CBZ-E 대사체), LTG(DDI 민감) |
| 8. AED 약력학 | GABA-T 억제, SV2A 결합, Nav 차단, AMPA 길항, Cav α2δ-1, mTORC1 억제 |
| 9. 약물 내성 기전 | P-gp(ABCB1), MRP2, BCRP, 혈뇌장벽 배출, Nav/GABA-A 수용체 변이 |
| 10. 임상 결과 | 발작 빈도, ILAE 분류, 반응률, 발작소실, SUDEP 위험, QoL(QOLIE-89) |

---

## 주요 병태생리 경로 (Key Pathophysiological Pathways)

### 1. SCN1A 유전자 돌연변이 → Dravet 증후군
```
SCN1A LOF variant
    → Nav1.1 단백질 기능 소실
    → PV+ 억제성 사이신경세포의 활동전위 발화 장애
    → E/I 불균형 (억제 감소 > 흥분 증가)
    → 발작역치 감소 → 발열 유발 발작 → Dravet 증후군
```

### 2. 내측두엽경화증 (MTS) → 측두엽 뇌전증 (TLE)
```
해마 CA1/CA3 신경세포 소실
    → 이끼섬유 발아(mossy fiber sprouting)
    → CA3 재귀 회로 형성 (이상 회로)
    → 동기화 돌발파 (burst firing)
    → 측두엽 발작 발생
```

### 3. GABA-A 수용체 내재화 (SE 기전)
```
5분 이상 지속 발작 (확립된 SE)
    → GABA-A 수용체 클라트린 매개 내재화
    → 시냅스 후막의 GABA-A 밀도 감소
    → BZD 반응성 감소 (수용체 수 감소)
    → 불응성 SE 전환
    → NMDA 수용체 과활성화 → 흥분독성 신경세포 사멸
```

### 4. P-gp 과발현 → 약물 내성 뇌전증 (DRE)
```
반복적 발작 → NF-κB 활성화
    → ABCB1(MDR1/P-gp) 전사 증가
    → 혈뇌장벽 AED 배출 증가
    → 뇌 내 AED 노출 감소
    → 치료 실패 → DRE 정의 충족 (≥2 AED 실패)
```

---

## 약물 PK/PD 파라미터 요약

| 약물 | 작용 기전 | Vc (L) | CL (L/h) | t½ | 치료범위 |
|------|---------|--------|----------|-----|---------|
| VPA | GABA-T 억제 + Nav | 9.1 | 0.47 | 9–17h | 50–100 mcg/mL |
| LEV | SV2A 결합 | 42 | 3.8 (신장) | 7h | 12–46 mcg/mL |
| CBZ | Nav 서서히 불활성화 | 51 | 3→6.5 (자가유도) | 8–12→5–8h | 4–12 mcg/mL |
| LTG | Nav 상태의존 차단 | 77 | 1.5 (UGT1A4) | 25–36h | 3–14 mcg/mL |
| PHT | Nav 불활성화 안정화 | 50 | 비선형(Km~4) | 22h | 10–20 mcg/mL |
| GBP | Cav2.2 α2δ-1 | 58 | 6.6 (신장) | 5–7h | 2–20 mcg/mL |

### 주요 약물상호작용 (DDI)

| 조합 | 기전 | 임상 영향 |
|------|------|---------|
| VPA + LTG | VPA → UGT1A4 억제 | LTG t½ 2배 증가 (24→48h), 독성 위험↑ |
| CBZ + LTG | CBZ → CYP3A4 유도 | LTG 청소율 2배 증가, 용량 2배 필요 |
| VPA + PHT | VPA → CYP2C9 억제 | PHT 독성 위험 (포화 동력학 주의) |
| CBZ + 경구피임약 | CBZ → CYP3A4 유도 | 에스트로겐 청소율 증가 → 실패 위험 |

---

## mrgsolve ODE 모델 구조

### 구획 (Compartments, 16개)
```
VPA:  AGUT → ACENT ⇌ APER → 대사
LEV:  BGUT → BCENT → 신장 배설
CBZ:  CGUT → CCENT → CMETA(CBZ-E) → 배설
LTG:  DGUT → DCENT → UGT1A4 배설
PD:   GABA(간접반응) · SYNAP(글루타메이트) · SV2A_OCC · NAV_BLOCK
동적:  STHRES(발작역치) · PGP(P-gp 발현)
```

### 간접반응 모델 (GABA, VPA)
```
dGABA/dt = kin × (1 + I_max_VPA) - kout × GABA
여기서 I_max_VPA = Imax × (CNS_VPA^n) / (IC50^n + CNS_VPA^n)
```

### 발작 빈도 모델
```
SeizFreq(t) = SeizBasal × exp[-k_seiz × (STHRES(t) - STHRES0)]
```

---

## 치료 시나리오 (10개)

| 시나리오 | 약물 | 용량 | 핵심 내용 |
|---------|------|------|---------|
| 1 | 없음 | — | 기저 발작 빈도 (8회/월) |
| 2 | VPA | 1,000 mg/day BID | GABA-T 억제 → GABA↑ → 역치↑ |
| 3 | LEV | 3,000 mg/day BID | SV2A 결합 → 시냅스 소포 방출↓ |
| 4 | CBZ | 600 mg/day BID | Nav 자가유도 포함, 2-4주 후 CL↑ |
| 5 | LTG | 200 mg/day BID | 서서히 증량 필요 |
| 6 | VPA+LTG | 500/100 mg/day | DDI: LTG t½ 2배 → 낮은 LTG 용량 |
| 7 | CBZ+LTG | 600/400 mg/day | DDI: LTG CL 2배 → 높은 LTG 용량 필요 |
| 8 | VPA (DRE) | 1,000 mg/day | P-gp 3배 → CNS 노출↓ → 효과↓ |
| 9 | IV BZD (SE) | 응급 | 초기 SE에서 BZD로 발작역치 즉각 회복 |
| 10 | VPA+에버롤리무스 | TSC | mTORC1 억제 → 피질 이형성 경로 차단 |

---

## Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 연령·체중·발작 유형·유전자·MRI·병력 입력 |
| 2. PK 프로파일 | VPA/LEV/CBZ/LTG 혈장 농도-시간 곡선 (정상상태 포함) |
| 3. PD 바이오마커 | GABA 수준·글루타메이트·SV2A 점유율·Nav 차단율 |
| 4. 임상 결과 | 발작 빈도·반응률·발작소실 확률·발작역치 동태 |
| 5. 시나리오 비교 | 단독요법 vs 병용요법 vs DRE 비교 대시보드 |
| 6. 내성·위험도 | P-gp 발현 추이·혈뇌장벽 노출비·SUDEP 위험 계산·mTOR |

---

## 파일 목록 (File Index)

| 파일 | 설명 |
|------|------|
| [epi_qsp_model.dot](epi_qsp_model.dot) | Graphviz 기계론적 지도 소스 (136 노드, 10 클러스터) |
| [epi_qsp_model.svg](epi_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [epi_qsp_model.png](epi_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [epi_mrgsolve_model.R](epi_mrgsolve_model.R) | mrgsolve ODE 모델 (16구획, 10 시나리오) |
| [epi_shiny_app.R](epi_shiny_app.R) | Shiny 대시보드 (6탭 인터랙티브 앱) |
| [epi_references.md](epi_references.md) | 참고문헌 62편 (PubMed 링크 포함) |

---

## 주요 임상시험 근거

| 약물 | 임상시험 | 결과 |
|------|---------|------|
| LEV | KEEPER (Cereghino 2000) | 발작 빈도 26% 감소, 반응률 33% |
| LTG | Matsuo 1993 | 위약 대비 발작 빈도 25% 감소 |
| CBZ | Mattson 1985 | 단순부분발작·GTCS에서 우수 |
| VPA | Chadwick 1999 | 전신 발작에서 CBZ와 동등 효과 |
| 에버롤리무스 (TSC) | EXIST-3 (Curatolo 2016) | 발작 빈도 41% 감소 |
| 펜플루라민 (Dravet) | PHENOMENON 2019 | 발작 빈도 71% 감소 |

---

## 기술 노트

- 전체 모델은 ODE 기반 (mrgsolve) + 발작역치 간접반응 모델
- CBZ 자가유도는 Michaelis-Menten 포화 함수로 근사
- P-gp 발현은 발작 부하에 의해 구동되는 1차 역학
- DDI 시뮬레이션에서 VPA-LTG와 CBZ-LTG 상호작용 구현
- SUDEP 위험은 반정량적 점수 시스템 (야간 GTCS 유무, 발작 조절 상태)

---

*Generated by QSP Disease Model Library | Claude Code Routine | 2026-06-21*
