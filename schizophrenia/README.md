# Schizophrenia QSP Model

**분류**: 신경정신질환 (Neuropsychiatric Disorder)  
**발병기전**: 도파민 과활성(메솔림빅) / 저활성(메소코티컬) + NMDA 기능저하 + GABAergic PV 인터뉴런 결핍 + 세로토닌 이상 + 신경염증  
**날짜**: 2026-06-20

---

## 1. 기계론적 지도 (Mechanistic Map)

[![Schizophrenia QSP Model](sch_qsp_model.png)](sch_qsp_model.svg)

> 클릭하면 전체 해상도 SVG를 볼 수 있습니다.

### 주요 클러스터 (10개)

| # | 클러스터 | 핵심 내용 |
|---|----------|-----------|
| ① | 신경발달 위험인자 | DISC1, NRG1, DTNBP1, COMT Val158Met, C4A; 환경 스트레스, 대마초 |
| ② | 도파민 경로 | 메솔림빅(↑), 메소코티컬(↓), 니그로선조체(EPS), 관상하구뇌하수체(PRL); D1/D2/D3 수용체 |
| ③ | 글루타메이트/NMDA | NMDA 기능저하 → PV 인터뉴런 탈억제 → 피질 탈억제 → 피질하 DA ↑ |
| ④ | 세로토닌 | DRN, 5-HT2A/2C/1A; SGA의 5-HT2A 차단이 메소코티컬 DA 회복에 핵심 |
| ⑤ | GABAergic 인터뉴런 | PV+ 세포 결핍 (GAD67 ↓ 25-50%), 감마 진동 소실 → 인지 결핍 |
| ⑥ | 신경염증/산화스트레스 | IL-6, IL-1β, TNF-α ↑; C4A 보체 → 시냅스 과잉 제거; BDNF ↓ |
| ⑦ | 항정신병약물 PK | HAL, RIS/PALI, CLZ, ARI — 2구획 모델, 대사물 포함 |
| ⑧ | 약력학/수용체 점유 | D2 점유 65-80% = 치료 윈도우; 5-HT2A 점유 >80% = SGA 음성증상 개선 |
| ⑨ | 임상 엔드포인트 | PANSS 양성/음성/일반, RBANS 인지, 재발율, 기능 회복 |
| ⑩ | 부작용 | EPS, 지연성 운동장애, 대사증후군, 고프로락틴혈증, QTc 연장, 무과립구증 |

---

## 2. mrgsolve ODE 모델 (`sch_mrgsolve_model.R`)

### 구획 (22개)

| 구획 유형 | 구획 이름 |
|-----------|-----------|
| **HAL PK** | GUT_HAL, CENT_HAL, PERI_HAL |
| **RIS/PALI PK** | GUT_RIS, CENT_RIS, PERI_RIS, CENT_PALI |
| **CLZ PK** | CENT_CLZ |
| **ARI/dARI PK** | GUT_ARI, CENT_ARI, CENT_dARI |
| **도파민 PD** | DA_MESOLIM, DA_MESOCORT, DA_NIGROSTR |
| **기타 PD** | PRL_CMPT, PV_ACT |
| **임상** | PANSS_POS, PANSS_NEG, PANSS_GEN |
| **바이오마커** | BDNF_CMPT, IL6_CMPT, EPS_RISK |

### 약물 PK 파라미터 (주요)

| 약물 | F (%) | CL (L/h) | Vc (L) | t½ | Brain Kp |
|------|--------|-----------|--------|-----|----------|
| Haloperidol (FGA) | 65 | 15 | 20 | 18-24h | 12 |
| Risperidone (SGA) | 74 | 25 | 30 | 3h(→21h PALI) | 7 |
| Clozapine (TRS) | 55 | 30 | 50 | 12h | 6 |
| Aripiprazole (부분 D2) | 87 | 3.6 | 245 | 75h | ~15 |

### 치료 시나리오 (7개)

| # | 시나리오 | 대상 |
|---|----------|------|
| 1 | 무치료 (자연경과) | 비교 기준 |
| 2 | Haloperidol 10 mg/d (FGA 표준) | 1세대 |
| 3 | Haloperidol 5 mg/d (저용량) | 저용량 FGA |
| 4 | Risperidone 4 mg/d (SGA) | 2세대 표준 |
| 5 | Clozapine 300 mg/d (TRS) | 치료 저항성 |
| 6 | Aripiprazole 15 mg/d (부분 D2) | 부분 효능제 |
| 7 | Risperidone 2 mg/d (저용량 SGA) | 저용량 SGA |

### 주요 임상 보정 출처
- **CATIE 2005** (Lieberman JA, NEJM): 항정신병약물 비교 효능
- **EUFEST 2008** (Kahn RS, Lancet): 초발 조현병 약물 비교
- **Kapur 2000** (AJP): D2 점유 치료 윈도우 PET
- **Nordstrom 1995** (AJP): Clozapine 저 D2 점유 PET

---

## 3. Shiny 대시보드 (`sch_shiny_app.R`)

**6개 탭 구성**:

| 탭 | 내용 |
|----|------|
| ① 환자 프로파일 | 인구학적 정보, 기저 PANSS 설정, 치료 선택, 발병기전 요약 |
| ② PK 프로파일 | 혈장 농도-시간 곡선, 뇌 농도, 정상상태 접근 |
| ③ D2/5-HT2A 점유 | 수용체 점유율 vs 시간, 치료 윈도우 시각화, 약물 fingerprint |
| ④ PANSS 임상 엔드포인트 | PANSS 총점/하위척도, 도파민 경로 역학, 반응률 |
| ⑤ 시나리오 비교 | 7개 치료법 동시 비교 (PANSS, D2 점유, EPS) |
| ⑥ 바이오마커 | 프로락틴, EPS 위험지수, BDNF, IL-6, PV 인터뉴런 활성 |

**실행 방법**:
```r
library(shiny); runApp("sch_shiny_app.R")
```

---

## 4. 참고문헌 (`sch_references.md`)

45개 PubMed 참고문헌 (11개 섹션):
- 역학 및 임상 개요
- 도파민 가설
- 글루타메이트/NMDA 수용체 가설
- GABAergic 인터뉴런 및 회로 기능 이상
- 세로토닌계
- 신경염증 및 산화스트레스
- 유전학 및 신경발달
- 항정신병약물 약동학
- D2/5-HT2A 점유 (PET 연구)
- 임상시험 및 치료 효능
- QSP/전산 모델링

---

## 5. 파일 목록

```
schizophrenia/
├── sch_qsp_model.dot        # Graphviz 기계론적 지도 (10 clusters, 160+ nodes)
├── sch_qsp_model.svg        # SVG 벡터 이미지
├── sch_qsp_model.png        # PNG 이미지 (150 dpi)
├── sch_mrgsolve_model.R     # mrgsolve ODE 모델 (22구획, 7 시나리오)
├── sch_shiny_app.R          # Shiny 대시보드 (6탭)
├── sch_references.md        # 참고문헌 45개
└── README.md                # 이 파일
```

---

## 6. 조현병 병태생리 핵심 요약

### 도파민 이중 가설 (Dual Dopamine Hypothesis)

```
메솔림빅 경로:  VTA → NAc    (↑ SCZ) → 양성증상 (환각, 망상)
메소코티컬 경로: VTA → DLPFC  (↓ SCZ) → 음성증상 + 인지결핍
니그로선조체:    SNc → 선조체  (D2 >80% 차단) → EPS
관상하구-하수체: 시상하부→뇌하수체 (D2차단) → 고프로락틴혈증
```

### SGA vs FGA의 차이

| | FGA (예: Haloperidol) | SGA (예: Risperidone) |
|-|------------------------|------------------------|
| D2 점유 | 높음 (~78%) | 높음 (~75%) |
| 5-HT2A 점유 | 낮음 (~30%) | 높음 (~96%) |
| 양성증상 | 우수 | 우수 |
| 음성증상 | 제한적 | 개선 (5-HT2A → 메소코티컬 DA ↑) |
| EPS | 높음 | 낮음 (5-HT2A → 니그로선조체 DA 회복) |
| 대사 | 낮음 | 중간-높음 (5-HT2C, H1, M1 차단) |

### Aripiprazole의 특이성
- **부분 D2 효능제**: DA 과활성 시 → 길항 효과; DA 저활성 시 → 효능제 효과
- **도파민 안정화**: 메솔림빅(↓) + 메소코티컬(↑) 균형
- **5-HT1A 부분 효능제**: PV 인터뉴런 회복 + 메소코티컬 DA ↑
