# qsp

## mrgsolve

- <https://vantage-research.net/qsp-in-r/>
- gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve
    - <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/>
    - <https://github.com/Genentech/gPKPDviz/>
    



## iqrtools

- <https://www.intiquan.com/acop2019_qsp/>

---

## QSP Disease Model Library

매일 Claude Code Routine이 추가하는 질환별 QSP 모델 라이브러리입니다.  
각 모델은 기계론적 지도(.dot/.svg/.png) · mrgsolve ODE 모델(.R) · Shiny 앱(.R) · 참고문헌(.md)으로 구성됩니다.

| 날짜 | 분류 | 질환명 | 디렉토리 | 모델 미리보기 |
|------|------|--------|----------|--------------|
| 2026-06-19 | 만성질환 | 요로결석 (만성 재발성) | [urolithiasis/](urolithiasis/) | [![URI QSP](urolithiasis/uri_qsp_model.png)](urolithiasis/uri_qsp_model.svg) |
| 2026-06-19 | 자가면역/혈관염 | 호산구 육아종증 다발혈관염 (EGPA) | [egpa/](egpa/) | [![EGPA QSP](egpa/egpa_qsp_model.png)](egpa/egpa_qsp_model.svg) |
| 2026-06-19 | 만성질환 | 원발성 부갑상선 기능 항진증 (PHPT) | [primary-hyperparathyroidism/](primary-hyperparathyroidism/) | [![PHPT QSP](primary-hyperparathyroidism/phpt_qsp_model.png)](primary-hyperparathyroidism/phpt_qsp_model.svg) |

---

## 요로결석 (만성 재발성, Chronic Recurrent Urolithiasis)

### 개요

요로결석은 미국 성인 인구의 약 8.8%(남성 10.6%, 여성 7.1%)에서 발생하며, 5년 재발률은 50%에 달하는 만성 재발성 질환입니다. 발병기전은 소변 내 결석 형성 물질(칼슘, 수산염, 요산)의 과포화(supersaturation)와 결석 형성 억제 물질(구연산, 마그네슘, THP)의 부족이 복합적으로 작용합니다.

### 주요 병태생리 경로

| 경로 | 핵심 메커니즘 | 임상 이상 |
|------|------------|---------|
| 칼슘 항상성 | PTH↑ → 1,25(OH)₂D↑ → 장관 Ca²⁺ 흡수↑ | 고칼슘뇨증 (>300mg/day) |
| 수산염 대사 | AGXT 결핍(PH1) / 지방 흡수장애 → 장관 OX↑ | 고수산뇨증 (>45mg/day) |
| 요산 대사 | XO 과활성 / URAT1 변이 → UA↑ | 고요산뇨증 (>800mg/day) |
| 구연산 처리 | 대사산증 / RTA → 신세뇨관 구연산 재흡수↑ | 저구연산뇨증 (<320mg/day) |
| 소변 과포화 | Ca × OX / Ksp > 1 → 핵화(nucleation) | CaOx SS > 1.0 |
| Randall's 플라크 | 상피하 인회석 침착 → 결석 핵 형성 nidus | 결석 성장 |

### 약물 PK/PD 파라미터

| 약물 | 작용기전 | 주요 PK | 임상 효과 |
|------|---------|---------|---------|
| HCTZ 25mg/day | NCC 억제 → 원위세뇨관 Ca²⁺ 재흡수↑ | F=0.65, t½=6-15h, CL=18L/h | 요중 Ca 30-45% 감소 |
| Allopurinol 300mg/day | XO 기전불활성화(mechanism-based) | F=0.90, Oxypurinol t½=18-30h | UA 생성 40-60% 감소 |
| K-Citrate 60mEq/day | 요중 구연산↑ + 요 pH 알칼리화 | F=0.95, CL=15L/h | CaOx SS 50% 감소 |
| Tamsulosin 0.4mg/day | α₁A/D-수용체 차단 → 요관 이완 | F=0.90, t½=9-16h | 자연 배출율 28% 증가 |
| Lumasiran (siRNA) | GalNAc-간 표적 → HAGO1 억제 | 월 1회 피하주사 | 요중 OX 53% 감소 (PH1) |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [uri_qsp_model.dot](urolithiasis/uri_qsp_model.dot) | Graphviz 기계론적 지도 소스 (100+ 노드, 10 클러스터) |
| [uri_qsp_model.svg](urolithiasis/uri_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [uri_qsp_model.png](urolithiasis/uri_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [uri_mrgsolve_model.R](urolithiasis/uri_mrgsolve_model.R) | mrgsolve ODE 모델 (17 구획, 7 시나리오) |
| [uri_shiny_app.R](urolithiasis/uri_shiny_app.R) | Shiny 대시보드 (6탭: 환자/PK/요화학/결석위험/시나리오비교/바이오마커) |
| [uri_references.md](urolithiasis/uri_references.md) | 참고문헌 38편 (PubMed 링크 포함) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **미치료 CaOx 결석 형성자**: 5년간 결석 성장 및 GFR 저하 추적
2. **HCTZ 25mg/day**: 고칼슘뇨증 환자에서 요중 Ca 30-45% 감소
3. **K-Citrate 60mEq/day**: 저구연산뇨증/UA 결석에서 CaOx SS 50% 감소
4. **Allopurinol 300mg/day**: 대사증후군/고요산뇨증에서 UA 생성 억제
5. **Lumasiran (PH1)**: 원발성 고수산뇨증에서 요중 OX 53% 감소
6. **생활습관 + HCTZ + K-Citrate 병용**: 종합 치료에서 결석 성장 억제
7. **MetSyn + UA 결석 + Allopurinol + K-Citrate**: 복합 요산 결석 관리

---

## 호산구 육아종증 다발혈관염 (EGPA, Eosinophilic Granulomatosis with Polyangiitis)

### 개요

EGPA(구 Churg-Strauss 증후군)는 유병률 약 1-3/백만명/년의 희귀 ANCA 연관 혈관염으로, 3단계 특징적 경과를 보입니다:
- **전구기 (Prodromal)**: 성인 발병 천식 + 알레르기 비염/비용종 (수년 지속)
- **호산구기 (Eosinophilic)**: 혈중 호산구 >10% (>1.5×10⁹/L) + 폐/심장/위장관 침윤
- **혈관염기 (Vasculitic)**: 소-중형 혈관 괴사성 혈관염 (단발신경염, 심장, 신장, 피부)

ANCA(항MPO)는 약 40%에서 양성이며 신장/신경 침범 위험이 높고, ANCA 음성에서는 심장·폐 호산구 침윤이 두드러집니다.

### 주요 병태생리 경로

| 경로 | 핵심 분자 | 임상 결과 |
|------|----------|---------|
| Th2 분극화 | TSLP/IL-33 → ILC2 → IL-5/IL-13/IL-4 | 호산구증가, IgE 상승, 기도 과반응 |
| 호산구 생산 | IL-5 → 골수 생산 → 혈중 이동 (CCR3/Eotaxin) | 혈중 호산구 3,000-15,000/µL |
| 호산구 과립단백 | MBP/ECP/EPX/EDN 방출 | 심근염, 신경독성, 혈관 손상 |
| ANCA 경로 | 항MPO → FcγRIIa → 호중구 활성화 → NETs/ROS | 괴사성 혈관염, RPGN |
| 육아종 형성 | Th1/M1 대식세포 → 상피양 세포 → 거대세포 | 폐 침윤, 조직 섬유화 |
| 혈관 손상 | 내피 세포 손상 → 투과성↑ → 미세혈전 → 허혈 | 단발신경염, GI 경색, 신장 손상 |

### Five-Factor Score (FFS) — 예후

| 항목 | 점수 |
|------|------|
| Creatinine > 150 µmol/L | +1 |
| 단백뇨 > 1g/일 | +1 |
| GI 침범 (출혈/천공) | +1 |
| 심근병증 | +1 |
| CNS 침범 | +1 |

FFS 0: 5년 사망률 ~11% | FFS ≥2: ~26%

### 약물 PK/PD 파라미터

| 약물 | 용량/경로 | 반감기 | 표적 | 주요 효과 |
|------|---------|--------|------|---------|
| Prednisolone | 50mg/일 경구 → 테이퍼 | ~2.5h | GR (광범위) | 호산구 60-80%↓, 혈관염 억제 |
| Mepolizumab | 300mg SC q4w | ~16-22일 | IL-5 (TMDD) | 호산구 90%↓, 재발률 50%↓ |
| Benralizumab | 30mg SC q4w→q8w | ~15일 | IL-5Rα (ADCC) | 호산구 근제로 감소 |
| Cyclophosphamide | 15mg/kg IV ×6 | ~7h (활성대사체) | DNA 알킬화 | ANCA 80%↓, B세포 소멸 |
| Rituximab | 1000mg IV ×2 | ~22일 | CD20 (ADCC/CDC) | B세포/형질세포 제거, ANCA↓ |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [egpa_qsp_model.dot](egpa/egpa_qsp_model.dot) | Graphviz 기계론적 지도 (150+ 노드, 10 클러스터) |
| [egpa_qsp_model.svg](egpa/egpa_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [egpa_qsp_model.png](egpa/egpa_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [egpa_mrgsolve_model.R](egpa/egpa_mrgsolve_model.R) | mrgsolve ODE 모델 (22 구획, 6 시나리오) |
| [egpa_shiny_app.R](egpa/egpa_shiny_app.R) | Shiny 대시보드 (7탭) |
| [egpa_references.md](egpa/egpa_references.md) | 참고문헌 60편 (PubMed 링크 포함) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **무치료 (자연 경과)**: 혈관염 악화 및 장기 손상 누적
2. **Prednisolone 단독**: 초기 50mg/일 → 26주 테이퍼 (표준 1차)
3. **Mepolizumab + Prednisolone**: MIRRA 임상시험 기반, 항IL-5 + 스테로이드 절감
4. **Benralizumab + Prednisolone**: 항IL-5Rα ADCC 매개 호산구 근제로 감소
5. **Cyclophosphamide + Prednisolone (중증 ANCA+)**: FFS ≥2 환자의 유도 치료
6. **Rituximab + Prednisolone (불응성 ANCA+)**: CD20 표적 B세포 소멸로 ANCA 억제

### Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | EGPA 개요, ACR/EULAR 2022 분류기준, FFS, 치료 알고리즘 |
| 2. 약물 PK | 혈청 약물 농도 추적 (Mepo/Benra/Pred/Cyclo/Ritu) |
| 3. 호산구/IL-5 | 혈중·조직 호산구, 유리 IL-5, 총IgE 동태 |
| 4. 혈관염/장기 | 혈관염 활성, ANCA, 심장 손상, 신경 손상 점수 |
| 5. 임상 엔드포인트 | BVAS, FEV1%, LVEF%, eGFR, Week 52 요약 |
| 6. 시나리오 비교 | 4개 표준 치료 시나리오 동시 비교 |
| 7. 바이오마커 | 호산구 억제율, 관해 상태, 바이오마커 히트맵 |

---

## 원발성 부갑상선 기능 항진증 (Primary Hyperparathyroidism, PHPT)

### 개요

PHPT는 부갑상선 선종(85%), 다발성 선종(3%), 증식(12%), 암종(1%)에 의한 PTH 자율과다분비로 유발되는 질환입니다. 미국 성인 인구의 약 0.3-1%에서 발생하며, 폐경 후 여성에서 유병률이 높습니다 (약 0.5%). 현재 대부분(약 75-80%)은 무증상으로 발견되며, "Bones, Groans, Moans, Stones"로 요약되는 고전적 임상 증후군은 드물어졌습니다.

### 주요 병태생리 경로

| 경로 | 핵심 분자 | 임상 결과 |
|------|----------|---------|
| CaSR-PTH 피드백 상실 | 선종 → CaSR 발현↓ → PTH 자율분비 | iPTH > 65 pg/mL (정상: 15-65) |
| 뼈 칼슘 동원 | PTH → RANKL↑/OPG↓ → OC 활성↑ | 혈청 Ca > 2.6 mM, CTX↑, P1NP↑ |
| 신장 Ca 재흡수 | PTH → TRPV5↑ → DCT Ca 재흡수↑ | 요중 Ca > 300 mg/day (역설적 고칼슘뇨) |
| 1,25D 합성 | PTH → CYP27B1↑ → Calcitriol↑ | 장관 Ca 흡수↑ → 고칼슘혈증 악화 |
| 인산뇨 | PTH → NaPi2a 내재화 → FEPi↑ | 저인산혈증 (PO4 < 0.8 mM) |
| FGF23-Klotho | PTH↑ → FGF23↑ (음성 피드백 손상) | 1,25D 자기조절 부전 |
| BMD 감소 | OC > OB 불균형 → 피질골 우선 손실 | 전완골(요골) T-score ≤ -2.5 조기 |

### 약물 PK/PD 파라미터

| 약물 | 용량 | 기전 | 주요 PK | 임상 효과 |
|------|------|------|---------|---------|
| Cinacalcet | 30-180 mg/day po | CaSR allosteric potentiator | F=22%, Vc=55L, CL=125L/d, t½=6-8h | Ca↓ ~0.5 mg/dL, PTH↓ 30-50%; BMD 변화 없음 |
| Denosumab | 60 mg SC q6mo | 항RANKL mAb (TMDD) | F=62%, Vd=3.1L, t½=28일 | CTX↓ 70-80%, BMD LS↑ 5%/yr; PTH/Ca 영향 없음 |
| Alendronate | 70 mg/wk po | Farnesyl-PP synthase 억제 | F=0.6%, 골 반감기 ~10년 | CTX↓ 50-70%, BMD 안정화; Ca 변화 없음 |
| Zoledronate | 5 mg IV/year | 동상 | F=100% (IV), t½=146h | 골절 위험 감소; PHPT 데이터 제한 |
| Calcitriol | 0.25-1 µg/day | VDR 직접 활성화 | F~100%, t½=5-8h | CKD-PHPT에서 Ca/PO4 보정 |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [phpt_qsp_model.dot](primary-hyperparathyroidism/phpt_qsp_model.dot) | Graphviz 기계론적 지도 (110+ 노드, 12 클러스터) |
| [phpt_qsp_model.svg](primary-hyperparathyroidism/phpt_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [phpt_qsp_model.png](primary-hyperparathyroidism/phpt_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [phpt_mrgsolve_model.R](primary-hyperparathyroidism/phpt_mrgsolve_model.R) | mrgsolve ODE 모델 (20 구획, 7+1 시나리오) |
| [phpt_shiny_app.R](primary-hyperparathyroidism/phpt_shiny_app.R) | Shiny 대시보드 (6탭: 프로파일/PK/Ca·PTH·VitD/골리모델링/시나리오비교/바이오마커) |
| [phpt_references.md](primary-hyperparathyroidism/phpt_references.md) | 참고문헌 65편 (PubMed 링크 포함) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **정상 (healthy baseline)**: Ca 1.20 mM, iPTH 52 pg/mL, BMD LS 0.960 g/cm²
2. **미치료 PHPT (경증, 5년 추적)**: Ca ~2.8 mM, PTH ~250 pg/mL, BMD 감소 1-2%/yr
3. **미치료 PHPT (중증, 5년 추적)**: Ca ~3.2 mM, PTH ~2800 pg/mL, Hungry-bone 위험
4. **Cinacalcet 60 mg/day**: 1년 후 Ca 정상화, PTH 40% 감소; BMD는 보전 효과 없음
5. **Denosumab 60 mg q6mo**: BMD LS +5%/yr; Ca/PTH 직접 효과 없음 (간접 OC 억제)
6. **부갑상선절제술 (day 90)**: 수술 후 72h 내 Ca 정상화, 장기 BMD 회복
7. **Cinacalcet + Denosumab 병용**: Ca 정상화 + BMD 보전 최적 조합 (수술 거부 환자)
8. **CKD-PHPT (eGFR 30) + Cinacalcet 90 mg**: 이차성 HPT와 동반 시 용량 증가 필요

### Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | PHPT 개요, NIH 2022 수술 기준, CaSR-PTH 시그모이드 곡선, 약물 요약 |
| 2. 약물 PK | Cinacalcet 혈중 농도, Denosumab TMDD PK, EC50 dose-response, PK 파라미터 표 |
| 3. Ca/PTH/VitD 동태 | iPTH·혈청 Ca·PO4·Calcitriol·요중 Ca·eGFR 시계열 |
| 4. 골 리모델링 & BMD | 요추/대퇴경부 BMD, T-score, OB/OC 역학, CTX/P1NP 마커 |
| 5. 시나리오 비교 | 7개 표준 치료 시나리오 선택 비교, Year-1 요약 테이블 |
| 6. 바이오마커 대시보드 | 실시간 valueBox (PTH/Ca/BMD/GFR), 전체 마커 테이블, Z-score 히트맵 |
