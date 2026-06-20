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
| 2026-06-19 | 만성질환 | 만성 정맥 부전 (CVI) | [chronic-venous-insufficiency/](chronic-venous-insufficiency/) | [![CVI QSP](chronic-venous-insufficiency/cvi_qsp_model.png)](chronic-venous-insufficiency/cvi_qsp_model.svg) |
| 2026-06-19 | 자가면역질환/혈관염 | IgA 혈관염 (헤녹-쇤라인 자반증, IgAV/HSP) | [iga-vasculitis/](iga-vasculitis/) | [![IgAV QSP](iga-vasculitis/igav_qsp_model.png)](iga-vasculitis/igav_qsp_model.svg) |
| 2026-06-19 | 자가면역질환/혈액 | 자가면역 용혈성 빈혈 (AIHA) | [autoimmune-hemolytic-anemia/](autoimmune-hemolytic-anemia/) | [![AIHA QSP](autoimmune-hemolytic-anemia/aiha_qsp_model.png)](autoimmune-hemolytic-anemia/aiha_qsp_model.svg) |
| 2026-06-19 | 자가면역질환/신장 | 굿파스처 증후군 (GPS) | [goodpasture-syndrome/](goodpasture-syndrome/) | [![GPS QSP](goodpasture-syndrome/gps_qsp_model.png)](goodpasture-syndrome/gps_qsp_model.svg) |
| 2026-06-19 | 만성질환 | 파젯병 (Paget's Disease of Bone, PBD) | [pagets-disease/](pagets-disease/) | [![PBD QSP](pagets-disease/pbd_qsp_model.png)](pagets-disease/pbd_qsp_model.svg) |
| 2026-06-20 | 자가면역/혈관염 | 현미경적 다발혈관염 (MPA) | [microscopic-polyangiitis/](microscopic-polyangiitis/) | [![MPA QSP](microscopic-polyangiitis/mpa_qsp_model.png)](microscopic-polyangiitis/mpa_qsp_model.svg) |
| 2026-06-20 | 자가면역질환/연골 | 재발성 다발연골염 (RP) | [relapsing-polychondritis/](relapsing-polychondritis/) | [![RPC QSP](relapsing-polychondritis/rpc_qsp_model.png)](relapsing-polychondritis/rpc_qsp_model.svg) |

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

---

## 만성 정맥 부전 (Chronic Venous Insufficiency, CVI)

### 개요

만성 정맥 부전은 전 세계 성인 인구의 약 35%가 경험하는 매우 흔한 만성 혈관 질환으로, 여성(55%)이 남성(45%)보다 유병률이 높습니다. 정맥 판막의 기능 부전으로 인한 역류(reflux)와 보행 시 정맥압 상승(elevated ambulatory venous pressure, AVP)이 핵심 발병기전이며, 장기 방치 시 정맥성 하지 궤양(venous leg ulcer, VLU)으로 진행할 수 있습니다. CVI는 CEAP 분류(C0-C6)에 따라 임상 중증도를 구분하며, C6는 활동성 궤양을 의미합니다.

### 주요 병태생리 경로

| 경로 | 핵심 메커니즘 | 임상 결과 |
|------|------------|---------|
| 정맥 판막 기능 부전 | 판막 무능증 → 역류 → AVP↑ (정상 <45 mmHg) | 정맥류, 하지 무거움 |
| 내피세포 기능 이상 | 비정상 전단 응력 → eNOS↓, ROS↑, ET-1↑ | ICAM-1/VCAM-1 발현↑ |
| 백혈구 트래핑 | ICAM-1 → PMN 포착 → 탈과립 → MBP·elastase | 모세혈관 손상, 조직 괴사 |
| 피리빈 섬유소 커프 | 피리빈 삼출 → 모세혈관 주위 섬유소 침착 | O₂/영양소 확산 장벽 |
| 부종 & 피부 변화 | 모세혈관 고혈압 → 혈장 삼출 → 림프 과부하 | 발목 부종, 지방피부경화증 |
| 궤양 형성 | 지방피부경화증 → 피부 장벽 파괴 → VLU | C6 분류, 삶의 질 저하 |

### CEAP 분류 체계

| CEAP 등급 | 임상 소견 | 유병률(성인) |
|----------|---------|-----------|
| C0 | 증상은 있으나 가시적 이상 없음 | 14% |
| C1 | 모세혈관 확장증/망상 정맥류 | 30% |
| C2 | 정맥류 | 23% |
| C3 | 부종 (피부 변화 없음) | 7% |
| C4 | 피부 변화 (색소침착, 습진, 지방피부경화증) | 4% |
| C5 | 치유된 궤양 + 피부 변화 | 1% |
| C6 | 활동성 정맥성 궤양 | 0.5% |

### 약물 PK/PD 파라미터

| 약물 | 작용기전 | 주요 PK | 임상 효과 |
|------|---------|---------|---------|
| MPFF 1000mg/day (Daflon) | 정맥 긴장도↑, 백혈구 부착↓, 혈관 투과성↓ | ka=1.5/h, CL=20L/h, V=80L, F=90% | VCSS 3.7점↓, 부종↓ (RELIEF meta-analysis) |
| Pentoxifylline 400mg TID | 혈액 점도↓, 피브리노겐↓, PAI-1↓, TNF-α↓ | ka=2.0/h, CL=45L/h, V=55L, F=25% | 궤양 치유율 1.7배↑ (Cochrane NNT=5) |
| Enoxaparin 40mg SC qd | 항Xa → 트롬빈↓, 피브린↓, DVT 예방 | ka=0.5/h, CL=1.2L/h, V=6L, F=90% | DVT 재발 50%↓, PTS 위험↓ |
| Rutosides 1000mg/day | 모세혈관 투과성↓, 산화 스트레스↓ | ka=1.0/h, 생체이용률 ~30% | 부종 23%↓, 다리 무거움 개선 |
| 압박 요법 20-40 mmHg | AVP 40-50%↓, VRT↑, 림프 배액↑ | 비약리적 기계적 효과 | 궤양 치유율 2배↑ (ESCHAR trial) |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [cvi_qsp_model.dot](chronic-venous-insufficiency/cvi_qsp_model.dot) | Graphviz 기계론적 지도 (146+ 노드, 11 클러스터) |
| [cvi_qsp_model.svg](chronic-venous-insufficiency/cvi_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [cvi_qsp_model.png](chronic-venous-insufficiency/cvi_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [cvi_mrgsolve_model.R](chronic-venous-insufficiency/cvi_mrgsolve_model.R) | mrgsolve ODE 모델 (16 구획, 7 시나리오) |
| [cvi_shiny_app.R](chronic-venous-insufficiency/cvi_shiny_app.R) | Shiny 대시보드 (6탭: 환자프로파일/PK/정맥혈역학/염증바이오마커/임상엔드포인트/시나리오비교) |
| [cvi_references.md](chronic-venous-insufficiency/cvi_references.md) | 참고문헌 54편 (PubMed 링크 포함) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **미치료 CVI (자연 경과, 5년)**: 판막 기능 저하 → AVP 점진적 상승 → 궤양 발생 위험 증가
2. **MPFF 1000mg/day 단독**: VCSS 3.7점↓, 부종↓, 삶의 질↑ (RELIEF meta-analysis 기준)
3. **Pentoxifylline + 압박 요법**: 활동성 궤양 환자에서 치유율 1.7배↑ (Cochrane 2012)
4. **LMWH + 압박**: DVT 후 PTS 예방, 피브린 커프↓
5. **삼중 병용 (MPFF + PTX + 압박)**: 중증 CVI/활동성 궤양 최적 조합
6. **압박 요법 단독**: ESCHAR trial 표준 치료
7. **MPFF + 압박 (중등도 CVI)**: 증상 개선 + 피부 변화 억제

### Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | CEAP 등급 입력, 위험인자 설정, 기저 질환 특성 요약 |
| 2. 약물 PK | MPFF/PTX/LMWH/루토사이드 혈중 농도 시계열, Cmax/AUC/t½ 표 |
| 3. 정맥 혈역학 | AVP·부종·VRT 시계열 (압박 유무 비교) |
| 4. 염증 바이오마커 | 백혈구 활성화·내피 기능·혈관 투과성·피브린 커프 추적 |
| 5. 임상 엔드포인트 | VCSS·CIVIQ-20 QoL·궤양 면적 시계열, 치료 시나리오별 개선율 |
| 6. 시나리오 비교 | 8개 치료군 레이더 차트·막대 그래프·요약 테이블 |

---

## IgA 혈관염 (IgA Vasculitis / Henoch-Schönlein Purpura, IgAV)

### 개요

IgA 혈관염(IgAV, 구칭: 헤녹-쇤라인 자반증/HSP)은 IgA 면역복합체가 소혈관(피부·신장·위장관·관절)에 침착되어 발생하는 전신 소혈관 혈관염입니다. 소아에서 가장 흔한 전신 혈관염이며(연간 발생률 10-20/100,000), 성인에서는 더 심한 신장 침범(IgAV 신염)을 보입니다. 병태생리는 "다중 타격(multi-hit)" 모델로 설명됩니다:

1. **Gd-IgA1 과잉 생산**: IgA1 힌지 영역의 비정상 O-당화 → 갈락토스 결핍 IgA1(Gd-IgA1) 과잉
2. **자가항체 형성**: 항-Gd-IgA1 IgG 자가항체 생성 → Gd-IgA1 면역복합체(IC) 형성
3. **IC 침착**: 피부 진피 혈관, 신장 메산지움, 위장관 점막, 활막에 침착
4. **보체 활성화**: 주로 렉틴 경로(MBL/MASP-2) → C3/C5/MAC 활성화
5. **혈관 염증**: 호중구·단핵구 동원, 내피 활성화, 백혈구 파괴성 혈관염(LCV)

### 주요 병태생리 경로

| 경로 | 핵심 분자 | 임상 결과 |
|------|----------|---------|
| Gd-IgA1 생산 | C1GalT1 결핍/ST6GalNAc II 과발현 → 힌지 갈락토스 결핍 | 혈청 Gd-IgA1 >2.0 mg/L |
| 면역복합체 형성 | Gd-IgA1 + 항-Gd-IgA1 IgG → 대형/소형 IC | 혈청 IC 상승 |
| 렉틴 보체 경로 | MBL-MASP-2 → C4/C2 → C3b → C5a/MAC | sC5b-9 상승, 사구체 손상 |
| 메산지움 침착 | FcαRI → NF-κB → PDGF-B/TGF-β → 메산지움 증식 | 단백뇨/혈뇨 |
| 족세포 손상 | TGF-β → nephrin/podocin↓ → 여과 장벽 파괴 | 신증후군 범위 단백뇨 |
| 피부 혈관 침착 | 진피 IgA 침착 → LCV → 자반 | 촉지성 자반 (하지) |
| 위장관 침착 | GI 점막 IgA 침착 → 복통/장중첩 위험 | 복통, GI 출혈 |

### 약물 PK/PD 파라미터

| 약물 | 용량 | 기전 | 주요 PK | 임상 효과 |
|------|------|------|---------|---------|
| Prednisolone | 1 mg/kg/day → 테이퍼 | GR 작용제 → NF-κB↓ | F=80%, t½=2.5h, CL=15L/h | 자반·GI 증상 호전, 신염 예방 불확실 |
| Mycophenolate mofetil (MMF) | 1000mg BID | IMPDH 억제 → 림프구 증식↓ | F=94%, MPA t½=18h | 단백뇨 40-55%↓ (Ren et al. 2012) |
| Rituximab | 375 mg/m² ×4 | 항CD20 → B세포 고갈 → Gd-IgA1↓ | t½=22일, TMDD 동역학 | B세포 소멸 → Gd-IgA1↓, 관해 유도 |
| Sparsentan | 400 mg QD | 이중 AT1R/ETA-R 길항 → TGF-β↓ | F=70%, t½=12h | 단백뇨 65%↓ (PROTECT 임상시험) |
| Dapagliflozin | 10 mg QD | SGLT2 억제 → TGF-β↓, GFR 보호 | F=78%, t½=12.9h | eGFR 저하율 50%↓ (DAPA-CKD) |
| ACEi/ARB | 에날라프릴 10mg/일 | 안지오텐신 II↓ → 수입세동맥 이완 | F=55%, t½=11h | 단백뇨 15-25%↓, TGF-β↓ |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [igav_qsp_model.dot](iga-vasculitis/igav_qsp_model.dot) | Graphviz 기계론적 지도 (200 노드, 10 클러스터) |
| [igav_qsp_model.svg](iga-vasculitis/igav_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [igav_qsp_model.png](iga-vasculitis/igav_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [igav_mrgsolve_model.R](iga-vasculitis/igav_mrgsolve_model.R) | mrgsolve ODE 모델 (25 구획, 8 시나리오) |
| [igav_shiny_app.R](iga-vasculitis/igav_shiny_app.R) | Shiny 대시보드 (7탭: 환자프로파일/PK/PD지표/임상엔드포인트/시나리오비교/바이오마커/About) |
| [igav_references.md](iga-vasculitis/igav_references.md) | 참고문헌 70편 (PubMed 링크 포함) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **미치료 (자연 경과)**: 52주간 단백뇨 진행, eGFR 저하, CKD 위험 누적
2. **Prednisolone 단독 (1 mg/kg/day × 4주 테이퍼)**: 자반·GI·관절 호전, 신염 예방 불확실
3. **Prednisolone + MMF 1000mg BID**: 단백뇨 55%↓, 면역복합체 억제
4. **Rituximab 375mg/m² × 4회**: B세포 고갈, Gd-IgA1↓, 관해 유도
5. **Prednisolone + Rituximab 병용**: 강력 면역억제 (중증 IgAV 신염)
6. **Sparsentan 400mg QD (PROTECT 기반)**: 단백뇨 65%↓, eGFR 보호
7. **ACEi/ARB (신보호 단독)**: 단백뇨 20%↓, 진행 억제
8. **Sparsentan + SGLT2 억제제 + 스테로이드 (병용 최적)**: 복합 신보호 + 면역억제

### Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 인구통계, 기저 IgAV 중증도 설정, 치료 선택, 질환 활성도 요약 |
| 2. 약물 PK | Prednisolone/MMF/Rituximab 혈중 농도 시계열 및 PK 파라미터 표 |
| 3. PD 핵심 지표 | 면역복합체·보체 동태, BAFF/IL-6/TNF-α 사이토카인 추적 |
| 4. 임상 엔드포인트 | 단백뇨·eGFR·혈뇨·자반점수·GI 점수, CKD 단계 위험 |
| 5. 시나리오 비교 | 8개 치료군 단백뇨/eGFR/IC 동시 비교, 최종 아웃컴 요약 테이블 |
| 6. 바이오마커 | Gd-IgA1·anti-Gd-IgA1·sC5b-9 기준범위, 바이오마커 추적 차트 |
| 7. About | IgAV 질환 개요, 병태생리 설명, QSP 모델 요약, 핵심 참고문헌 |

---

## 자가면역 용혈성 빈혈 (Autoimmune Hemolytic Anemia, AIHA)

### 개요

자가면역 용혈성 빈혈(AIHA)은 자가항체가 적혈구(RBC) 표면 항원에 결합하여 조기 파괴를 유발하는 후천성 혈액 질환입니다. 연간 발생률 약 1-3/100,000명이며, 여성에서 약간 많습니다. 크게 두 가지 주요 아형으로 나뉩니다:

1. **온형 AIHA (Warm AIHA, wAIHA, ~70%)**: IgG 자가항체가 37°C에서 Rh·Band3·CD47 항원에 결합 → FcγR 매개 비장 대식세포의 식세포 작용(체외 용혈)
2. **한랭 응집소 질환 (Cold Agglutinin Disease, CAD, ~15–20%)**: IgM 단클론 항체가 <15°C에서 I/i 항원에 결합 → 고전 보체 경로 활성화 (C1q→C3b→C5b-9 MAC) → 혈관내 용혈 + C3b 옵소닌화

병태생리 핵심: 면역관용 소실 → B세포/형질세포 자가항체 생산 → RBC 표면 항원 결합 → 식세포/보체 매개 용혈 → 보상적 적혈구생성(EPO↑, 망상적혈구↑)

### 주요 병태생리 경로

| 경로 | 핵심 분자 | 임상 결과 |
|------|----------|---------|
| IgG 생산 (warm) | GC반응 → SHM → 형질세포 → IgG1/IgG3 | DAT IgG+ |
| FcγR 매개 식세포 작용 | FcγRI/IIA/IIIA → Syk → 탐식 | 구상 적혈구, 비장 비대 |
| 보체 고전 경로 (cold) | IgM → C1q/C1r/C1s → C3b → MAC | DAT C3d+, 혈관내 용혈 |
| 보상 적혈구생성 | 저산소 → HIF-1α → EPO↑ → CFU-E↑ | 망상적혈구증가증 |
| 헤모글로빈 분해 (체외) | HO-1 → 비리베르딘 → 비결합 빌리루빈 | 황달, 빌리루빈↑ |
| 자유 헤모글로빈 (혈관내) | Hb → 합토글로빈 고갈 → NO 소거 | 혈색소뇨, 혈관수축 |

### 약물 PK/PD 파라미터

| 약물 | 용량 | 기전 | 주요 PK | 임상 효과 |
|------|------|------|---------|---------|
| Prednisolone | 1~1.5 mg/kg/day → 테이퍼 | GR → FcγR↓, Ab↓ | F=80%, t½=2.5h, CL=15L/h | CR ~55-70% at 3 weeks |
| Rituximab | 375 mg/m² IV ×4 주간 | 항CD20 → B세포 고갈 | t½=22일, TMDD 동역학 | CR ~65% at 1 year (Barcellini 2018) |
| Sutimlimab (Enjaymo) | 6.5-7.5g IV q2w | 항C1s → 고전 보체 경로 차단 | t½=20일, Cmax~1300μg/mL | Hb +1.5-2.3 g/dL (CADENZA trial) |
| Fostamatinib (Tavalia) | 150mg PO BID | Syk 억제 → FcγR 신호↓ | R406 t½=14h, IC50=41nM | 38% OR (Phase 3 trial 2023) |
| Dexamethasone | 40mg/day ×4d 사이클 | GR (Kd=1nM) → 강력 면역억제 | F=78%, t½=36-72h | 단기 신속 반응 (CR 60%) |
| MMF | 1000mg PO BID | IMPDH 억제 → B세포 증식↓ | MPA t½=18h, F=94% | 재발 방지 유지요법 |
| IVIG | 1 g/kg IV ×2일 | FcRn 포화→ Ab 분해↑ + FcγR 차단 | t½=21일, Vd=3.5L | 급성 구제 (단기 효과) |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [aiha_qsp_model.dot](autoimmune-hemolytic-anemia/aiha_qsp_model.dot) | Graphviz 기계론적 지도 (15 클러스터, 100+ 노드) |
| [aiha_qsp_model.svg](autoimmune-hemolytic-anemia/aiha_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [aiha_qsp_model.png](autoimmune-hemolytic-anemia/aiha_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [aiha_mrgsolve_model.R](autoimmune-hemolytic-anemia/aiha_mrgsolve_model.R) | mrgsolve ODE 모델 (26 구획, 8 시나리오) |
| [aiha_shiny_app.R](autoimmune-hemolytic-anemia/aiha_shiny_app.R) | Shiny 대시보드 (7탭: 환자프로파일/PK/PD지표/임상엔드포인트/시나리오비교/바이오마커/About) |
| [aiha_references.md](autoimmune-hemolytic-anemia/aiha_references.md) | 참고문헌 76편 (PubMed 링크 포함) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **미치료 (자연 경과)**: 1년간 hemoglobin 지속 저하, LDH↑, 빌리루빈↑, 망상적혈구증가증
2. **Prednisolone 1 mg/kg/day → 테이퍼**: 3주내 빠른 반응, 4주차 테이퍼, CR ~55-70%
3. **Prednisolone + Rituximab**: 복합 면역억제, 1년 CR ~68%, 재발률 감소 (Birgens 2013)
4. **Dexamethasone 펄스 + Rituximab**: 강력 초기 반응 + 지속 B세포 고갈
5. **Fostamatinib 150mg BID (불응성 wAIHA)**: Syk 경로 억제, OR ~38% (Phase 3)
6. **Prednisolone + MMF (유지요법)**: 스테로이드 감량 후 MMF 유지, 재발 방지
7. **IVIG 1 g/kg ×2d (급성 구제)**: 즉각적 FcγR 차단, 단기 효과, 수혈 회피
8. **Sutimlimab 6.5g q2w (CAD)**: 고전 보체 경로 완전 차단, Hb +2.3 g/dL (CADENZA)

### Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 인구통계, 기저 AIHA 아형(온형/한랭) 설정, 치료 선택, 기저치 요약 표 |
| 2. 약물 PK | Prednisolone/Rituximab/Sutimlimab/Fostamatinib/MMF/IVIG 혈중 농도 시계열 |
| 3. PD 핵심 지표 | B세포·형질세포·자가항체 동태, C3b 보체, EPO·망상적혈구 반응 |
| 4. 임상 엔드포인트 | 헤모글로빈 궤적, 완전반응(CR≥10g/dL)/부분반응/수혈 필요 시점 |
| 5. 시나리오 비교 | A/B/C/D 4개 치료군 동시 비교 (Hb, LDH, 자가항체 그래프 + 결과 테이블) |
| 6. 바이오마커 | LDH·합토글로빈·빌리루빈·DAT 점수 추적 및 기준 범위 참고표 |
| 7. About | AIHA 질환 개요, 모델 구조 설명, 참고문헌, 모델 한계 |

---

## 굿파스처 증후군 (Goodpasture Syndrome, GPS)

### 개요

굿파스처 증후군은 제4형 콜라겐 α3 사슬의 NC1 도메인(Goodpasture 항원)에 대한 자가항체(anti-GBM IgG)가 사구체 기저막(GBM)과 폐포 기저막을 공격하여 급속진행성 사구체신염(RPGN) 및 폐포 출혈(DAH)을 유발하는 희귀 자가면역 질환입니다. 연간 발생률은 100만 명당 0.5–1명으로 매우 드물며, 치료받지 않으면 수 주 내 신부전·호흡부전으로 사망에 이릅니다.

### 주요 병태생리 경로

| 경로 | 핵심 메커니즘 | 임상 이상 |
|------|------------|---------|
| 자가항원 노출 | 흡연·감염·용매 → GBM α3(IV)NC1 노출 | 면역 반응 개시 |
| B세포 활성화 | 항원 제시 → B세포 → 형질세포 분화 | anti-GBM IgG 생산 |
| 보체 활성화 | IgG → C1q → 고전경로 → C5a 생성 | 호중구 유입·조직 파괴 |
| 신장 손상 | 반월형 사구체신염 → GFR 급감 | RPGN (혈청 Cr ↑↑) |
| 폐 손상 | 폐포 기저막 IgG → 폐포 출혈 | DAH (DLCO ↓, 객혈) |
| 염증 증폭 | 호중구 → ROS·MMP → 기저막 파괴 | 단백뇨·혈뇨·폐섬유화 |

### 치료 시나리오 (mrgsolve 시뮬레이션)

1. **미치료 (자연 경과)**: anti-GBM 항체 지속, GFR 급속 저하, 폐출혈 지속
2. **혈장교환(Plex) + 사이클로포스파미드(CY) + 프레드니솔론(표준요법)**: 혈중 항체 신속 제거, B세포 억제, 3개월 내 anti-GBM 음전율 ~85%
3. **Plex + 리툭시맙(RTX) + 프레드니솔론**: CD20+ B세포 고갈, 재발성/불응성 환자에 적용
4. **Plex + 아바코판(Avacopan) + 프레드니솔론**: C5aR 차단으로 호중구 매개 손상 억제
5. **CY + 프레드니솔론 (혈장교환 없음)**: 혈장교환 불가 시 대안
6. **프레드니솔론 단독**: 증등증 폐 침범만 있는 경우 또는 대안

### 모델 구조 (mrgsolve ODE)

| 구획 | 상태변수 | 생물학적 의미 |
|------|---------|------------|
| Drug PK | CY_C, OHCY_C | 사이클로포스파미드 + 활성대사체 4-OH-CY (간 활성화) |
| Drug PK | PRED_C | 프레드니솔론 (1구획 경구 PK) |
| Drug PK | RTX_C, RTX_P | 리툭시맙 2구획 PK (중심/말초) |
| Drug PK | AVA_C | 아바코판 (C5aR1 저해제) |
| 면역 | AntiGBM | 혈중 anti-GBM IgG 항체 농도 |
| 면역 | B_cells, Plasma_cells | B세포·형질세포 동태 |
| 보체 | C5a | C5a 아나필라톡신 농도 |
| 염증 | Neutrophil_kidney | 사구체 호중구 침윤 |
| 신장 | GBM_damage, GFR_c | GBM 손상 지수, GFR (정규화) |
| 신장 | Proteinuria_c, Hematuria_c | 단백뇨·혈뇨 |
| 폐 | Lung_damage, DLCO_c | 폐포 손상, DLCO (정규화) |
| 염증 | CRP_c, T_regs | CRP, 조절T세포 |

### 기계론적 지도 파일

| 파일 | 설명 |
|------|------|
| [gps_qsp_model.dot](goodpasture-syndrome/gps_qsp_model.dot) | Graphviz 기계론적 지도 (9 클러스터, 100+ 노드) |
| [gps_qsp_model.svg](goodpasture-syndrome/gps_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [gps_qsp_model.png](goodpasture-syndrome/gps_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [gps_mrgsolve_model.R](goodpasture-syndrome/gps_mrgsolve_model.R) | mrgsolve ODE 모델 (19 구획, 6 시나리오) |
| [gps_shiny_app.R](goodpasture-syndrome/gps_shiny_app.R) | Shiny 대시보드 (7탭: 환자프로파일/PK/항체·보체/신장/폐/시나리오비교/참고문헌) |
| [gps_references.md](goodpasture-syndrome/gps_references.md) | 참고문헌 60편 (13 섹션, PubMed 링크 포함) |

### Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 질환 개요, 진단 기준, 치료 알고리즘, 예후 지표 설정 |
| 2. 약물 PK | 혈장교환 일정, CY/PRED/RTX/아바코판 혈중 농도 시계열 |
| 3. 항체 & 보체 | anti-GBM IgG 역가, B세포·형질세포 동태, C5a, 호중구 침윤 |
| 4. 신장 엔드포인트 | GFR 궤적, GBM 손상 지수, 단백뇨, 투석 위험도 |
| 5. 폐 엔드포인트 | DLCO 궤적, 폐포 손상 지수, 폐 완해율 |
| 6. 시나리오 비교 | 6개 치료군 동시 비교 (anti-GBM, GFR, DLCO, CRP 그래프 + 결과 테이블) |
| 7. 참고문헌 | 역학/항원/유전학/자가항체/보체/신장/폐/혈장교환/CY/RTX/아바코판/PK/PD/예후 별 분류 |

---

## 파젯병 (Paget's Disease of Bone, PBD)

### 개요

파젯병(Paget's Disease of Bone)은 국소적으로 비정상적으로 가속화된 골 개조(bone remodeling)를 특징으로 하는 만성 대사성 골질환입니다. 파젯 파골세포는 비정상적으로 크고(핵 최대 100개), 과활성 상태로 기계적 강도가 낮은 woven/mosaic bone을 생성합니다. SQSTM1/p62 돌연변이(가족성 사례의 30–50%), RANK 경로 조절장애, paramyxovirus 가설이 병인에 관여합니다.

### 병태생리 핵심 경로

| 경로 | 핵심 분자 | 임상 결과 |
|------|----------|---------|
| RANKL-OPG 불균형 | RANKL↑ / OPG↓ → NFATc1 과활성 | 파골세포 과형성, 골흡수↑↑ |
| SQSTM1/p62 돌연변이 | UBA 도메인 손상 → NF-κB 과활성 → 파골세포 과민화 | 가족성 PBD (30-50%) |
| 파골세포 형태 이상 | DC-STAMP↑, 핵 100개 이상, TRAP↑↑ | 국소 골용해 → 병적 골절 |
| 골형성 결합 과잉 | TGF-β1/IGF-1 방출 → OBpre 과자극 | Woven/mosaic bone 형성 |
| 혈관 합병증 | VEGF↑ → AV shunting → 심박출량↑ | 고박출 심부전 |

### 약물 PK/PD 파라미터

| 약물 | 용량/경로 | 반감기 | 표적 | 주요 효과 |
|------|---------|--------|------|---------|
| Zoledronic acid | 5mg IV 1회 | 골 결합 t½ ~10년 | FPPS (mevalonate 경로) | bsALP 정상화 89% (6개월) |
| Alendronate | 40mg/일 경구 × 6개월 | 골 결합 t½ ~10년 | FPPS | bsALP 정상화 ~70% |
| Calcitonin | 100IU SC 매일 | ~1h | CALCR (Gs/cAMP/PKA) | OC 억제 ~25-30% |
| Denosumab | 60mg SC Q6M | ~25-30일 | RANKL 중화 (IgG2) | NTX ~65-70%↓ |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [pbd_qsp_model.dot](pagets-disease/pbd_qsp_model.dot) | Graphviz 기계론적 지도 (154 노드, 12 클러스터, 201 엣지) |
| [pbd_qsp_model.svg](pagets-disease/pbd_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [pbd_qsp_model.png](pagets-disease/pbd_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [pbd_mrgsolve_model.R](pagets-disease/pbd_mrgsolve_model.R) | mrgsolve ODE 모델 (22 구획, 7 시나리오) |
| [pbd_shiny_app.R](pagets-disease/pbd_shiny_app.R) | Shiny 대시보드 (6탭: 환자프로파일/PK/PD지표/임상엔드포인트/시나리오비교/바이오마커) |
| [pbd_references.md](pagets-disease/pbd_references.md) | 참고문헌 53편 (PubMed 링크, 9섹션) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **미치료 (자연 경과)**: 2년간 bsALP·NTX·CTX 상승 및 BMD 변화 추적
2. **Zoledronic acid 5mg IV 1회**: HORIZON 임상시험 기반, 6개월 bsALP 정상화 89%
3. **Alendronate 40mg/일 × 6개월**: 경구 비스포스포네이트 표준 치료
4. **Calcitonin 100IU SC 매일 × 6개월**: CALCR 매개 파골세포 억제
5. **Denosumab 60mg SC Q6M**: RANKL 중화, 강력한 NTX 억제
6. **ZA 5mg IV + 통증 관리**: 증상 중심 병용 치료
7. **순차 치료**: Alendronate 6개월 → ZA 전환

### 주요 보정 기준 (Calibration Target)

- HORIZON 임상시험: ZA 투여 6개월 후 bsALP 정상화 89% (PMID: 21434807)
- Siris et al. 역학: 55세 이상 0.5-3% 유병률, 남성 다발

---

## 현미경적 다발혈관염 (Microscopic Polyangiitis, MPA)

### 개요

현미경적 다발혈관염(MPA)은 소혈관(arterioles, venules, capillaries)을 주로 침범하는 ANCA 연관 혈관염(AAV)입니다. 주요 ANCA는 항MPO IgG(pANCA, ~75%)이며, 면역복합체 침착이 거의 없는 **pauci-immune** 괴사성 사구체신염과 **미만성 폐포출혈(DAH)**이 대표적 합병증입니다. 육아종은 형성되지 않아 GPA와 구별됩니다. 연간 발생률은 유럽에서 100만 명당 2-10명이며, 일본에서 더 높습니다.

### 병태생리 핵심 경로

| 경로 | 핵심 분자 | 임상 결과 |
|------|----------|---------|
| ANCA 생성 | B세포 내성 소실 → 형질세포 → anti-MPO IgG (pANCA) | ANCA 역가↑ (>100 EU) |
| 호중구 프라이밍 | TNF-α/GM-CSF → MPO 세포표면 발현 | ANCA 결합 가능 상태 |
| ANCA-호중구 활성화 | FcγRIIa 교차결합 → NADPH oxidase → ROS/NETosis | 혈관벽 파괴 |
| 보체 대체경로 | C5a 생성 → C5aR1 → 호중구 추가 프라이밍 | 양성 피드백 증폭 |
| 사구체 모세혈관염 | pauci-immune 괴사 → RPGN → 초승달(crescent) 형성 | GFR↓, ESRD 위험 |
| 폐 모세혈관염 | 폐포 모세혈관 파열 → DAH | 객혈, DLCO↓ |

### 약물 PK/PD 핵심 파라미터

| 약물 | 용량/경로 | 반감기 | 표적 | 주요 효과 |
|------|---------|--------|------|---------|
| Cyclophosphamide | 2mg/kg/d oral 또는 15mg/kg IV Q3W | 3-10시간 (4-OH-CY t½~4h) | DNA 알킬화 | B/T세포 독성, 완전관해 ~53% (RAVE) |
| Rituximab | 375mg/m²×4 or 1000mg×2 IV | ~350시간 (3주 반감기) | CD20 (ADCC/CDC) | B세포 고갈, 완전관해 ~64% (RAVE) |
| Prednisolone | 1mg/kg/day → taper | 2-4시간 | GR (transrepression) | 염증 억제, 관해 유도 보조 |
| Azathioprine | 2mg/kg/day oral | ~1시간 (6-TGN t½ 3-13일) | HGPRT/퓨린 합성 | B/T세포 증식 억제 (유지) |
| Mycophenolate | 3g/day oral | ~12-17시간 (MPA) | IMPDH (IMP→GMP) | 형질세포 억제 (유지, 2선) |
| Avacopan | 30mg BID oral | ~16시간 | C5aR1 역작용제 | BVAS 관해: 비열등 vs PRED (ADVOCATE) |
| Plasma Exchange | 7회/2주 | 즉각 | ANCA/complement 제거 | ANCA 50-70%↓/1회, ESRD 영향 미미 (PEXIVAS) |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [![MPA QSP](microscopic-polyangiitis/mpa_qsp_model.png)](microscopic-polyangiitis/mpa_qsp_model.svg) | 기계론적 지도 미리보기 (클릭 시 SVG 확대) |
| [mpa_qsp_model.dot](microscopic-polyangiitis/mpa_qsp_model.dot) | Graphviz 기계론적 지도 (160+ 노드, 11 클러스터) |
| [mpa_qsp_model.svg](microscopic-polyangiitis/mpa_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [mpa_qsp_model.png](microscopic-polyangiitis/mpa_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [mpa_mrgsolve_model.R](microscopic-polyangiitis/mpa_mrgsolve_model.R) | mrgsolve ODE 모델 (21 구획, 7 시나리오) |
| [mpa_shiny_app.R](microscopic-polyangiitis/mpa_shiny_app.R) | Shiny 대시보드 (7탭) |
| [mpa_references.md](microscopic-polyangiitis/mpa_references.md) | 참고문헌 64편 (14섹션, PubMed 링크) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **미치료 (자연 경과)**: ANCA 지속, GFR 급속 저하, DAH 지속, ESRD 진행
2. **CY 경구 + Prednisolone (CYCLOPS 프로토콜)**: 경구 CY 2mg/kg/d + PRED 1mg/kg/d → 완전관해 ~53% (6개월)
3. **Rituximab + Prednisolone (RAVE 프로토콜)**: RTX 375mg/m²×4 + PRED → 완전관해 ~64% (6개월), 특히 재발 MPA에 우월
4. **CY + Prednisolone + 혈장교환 (PEXIVAS)**: PLEX 7회 추가 → ANCA 신속 제거, 그러나 ESRD/사망률 개선 미확인
5. **Rituximab + Avacopan (GC-free, ADVOCATE)**: 스테로이드 없이 C5aR1 차단 + B세포 고갈 → BVAS 비열등, 52주 지속 관해 우월
6. **AZA 유지 (post-induction, IMPROVE)**: 유도 후 AZA 2mg/kg/d → 재발률 ~35%(5년) — RTX 유지보다 열등
7. **RTX 유지 500mg Q6M (MAINRITSAN)**: 고정 간격 RTX 유지 → 재발률 대폭 감소 (AZA 대비 HR 0.36)

### Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 질환 개요, BVAS 구성, 바이오마커 기준치, 치료 알고리즘, 예후 |
| 2. 약물 PK | CY/4-OH-CY, RTX, Prednisolone, Avacopan 혈중농도 시계열 |
| 3. 면역학·ANCA | anti-MPO ANCA 역가, B세포·형질세포, C5a, 호중구 활성화, 혈관내피 손상 |
| 4. 신장 엔드포인트 | GFR 궤적, 혈청 크레아티닌, 사구체 염증 지수, 신장 섬유화 |
| 5. 폐 엔드포인트 | DAH 지수, DLCO, CRP (전신 염증) |
| 6. 시나리오 비교 | 7개 치료군 동시 비교 (BVAS, GFR, ANCA, DLCO 그래프 + 결과표) |
| 7. 참고문헌 | 핵심 임상시험(RAVE/RITUXVAS/PEXIVAS/ADVOCATE/MAINRITSAN) + 기전 문헌 |

### 주요 보정 기준 (Calibration Targets)

- RAVE 임상시험: RTX 완전관해 64% vs CY 53% (6개월) — PMID 20647198
- PEXIVAS: 혈장교환이 ESRD/사망률에 추가 이득 없음 — PMID 32053298
- ADVOCATE: 아바코판이 BVAS 관해에서 Pred 비열등, 52주 지속 관해 우월 — PMID 33596356
- MAINRITSAN: RTX 유지 시 재발률 HR 0.36 vs AZA — PMID 25372085

---

## 재발성 다발연골염 (Relapsing Polychondritis, RP)

### 개요

재발성 다발연골염(RP)은 귀·코·기관기관지 연골, 관절, 눈, 심혈관계, 내이를 반복적으로 침범하는 희귀 자가면역질환입니다. 유병률은 약 100만 명당 3.5명이며, 주요 자가항원은 **제2형 콜라겐(Type II collagen)**, Matrilin-1, 제9형·11형 콜라겐, COMP입니다. CD4+ Th1/Th17 세포 및 항-CII IgG가 연골에 침착되어 보체를 활성화시키고, MMPs/ADAMTS가 연골 ECM을 파괴합니다. 기도 허탈(tracheomalacia)과 심혈관 침범이 주요 사망 원인입니다.

### 병태생리 핵심 경로

| 경로 | 핵심 분자 | 임상 결과 |
|------|----------|---------|
| 자가항원 노출 | Type II 콜라겐/Matrilin-1 → DC 항원제시 | 자가반응 T세포·B세포 활성화 |
| Th1/Th17 분화 | IL-12/IL-23 → Th1(IFN-γ) / IL-6+TGF-β → Th17(IL-17A/F) | 연골 침윤 |
| 항-CII IgG 생성 | 형질세포 → Anti-CII IgG → IC 형성 | 연골 내 면역복합체 침착 |
| 보체 활성화 | IC → C1q → C3a/C5a anaphylatoxin → MAC | 연골세포 용해 |
| MMP/ADAMTS 활성화 | TNF-α, IL-17A → MMP-1/3/9/13, ADAMTS-4/5↑ | 콜라겐/프로테오글리칸 분해 |
| 연골세포 소실 | ROS, 퍼포린/그란자임 B, MAC → apoptosis | 섬유연골 대체 → 석회화 |
| 기도 합병증 | 기관기관지 연골 소실 → tracheomalacia | 기도 허탈·협착·호흡부전 |

### 약물 PK/PD 핵심 파라미터

| 약물 | 용량/경로 | 반감기 | 표적 | 주요 효과 |
|------|---------|--------|------|---------|
| Prednisolone | 40-60mg/d PO (유도) → 10mg/d (유지) | 2.5-4h (MW 360.4) | GR-α (NF-κB transrepression) | TNF/IL-6/IL-1β↓, 관해 유도 1선 |
| Dapsone | 100mg/d PO | 20-30h | MPO 억제 / ROS↓ | 경증 RP 1선, 호중구 기능 억제 |
| Methotrexate | 7.5-25mg/wk PO/SC | 3-10h (MTX-PG 수일) | DHFR/ATIC 억제 → adenosine↑ | Th1/Th17↓, B세포↓, steroid sparer |
| Tocilizumab | 8mg/kg q4w IV | ~240h (t½ ~10d) | IL-6R 차단 → JAK-STAT3↓ | CRP 정상화, 재발 빈도↓, 기도 안정 |
| Abatacept | 750mg q4w IV | ~400h (t½ ~17d) | CTLA4-Ig → CD28-B7 차단 | T세포 활성화↓, 연골 파괴 속도↓ |
| TNF 억제제 | 표준 용량 | 주 특성에 따라 | TNF-α 중화 | 재발 빈도↓, 일부 불응 |
| Rituximab | 375mg/m²×4 IV | ~350h | CD20 → B세포 고갈 | Anti-CII Ab↓, 불응 RP에 유용 |
| Colchicine | 0.5-1.5mg/d | 9-30h | 튜불린 중합/NLRP3↓ | 급성 발작 예방, 보조 |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [![RP QSP](relapsing-polychondritis/rpc_qsp_model.png)](relapsing-polychondritis/rpc_qsp_model.svg) | 기계론적 지도 미리보기 (클릭 시 SVG 확대) |
| [rpc_qsp_model.dot](relapsing-polychondritis/rpc_qsp_model.dot) | Graphviz 기계론적 지도 (154 노드, 10 클러스터) |
| [rpc_qsp_model.svg](relapsing-polychondritis/rpc_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [rpc_qsp_model.png](relapsing-polychondritis/rpc_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [rpc_mrgsolve_model.R](relapsing-polychondritis/rpc_mrgsolve_model.R) | mrgsolve ODE 모델 (20 구획, 7 치료 시나리오) |
| [rpc_shiny_app.R](relapsing-polychondritis/rpc_shiny_app.R) | Shiny 대시보드 (7탭) |
| [rpc_references.md](relapsing-polychondritis/rpc_references.md) | 참고문헌 64편 (12섹션, PubMed 링크) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **미치료 (자연경과)**: Anti-CII IgG 지속 축적, 보체 활성화, MMP 과발현, 연골 무결성 점진적 감소
2. **Prednisone 유도 → 유지 (60→10mg/d)**: GR 점유율 상승, NF-κB 억제, TNF/IL-6↓, 연골 안정화
3. **Prednisone + Methotrexate**: Th17 및 B세포 추가 억제, IL-17A↓, 스테로이드 절감 효과
4. **Tocilizumab 8mg/kg q4w**: IL-6R 완전 차단, JAK-STAT3 억제, CRP 정상화, 기도 합병증↓
5. **Abatacept**: CD28-B7 공동자극 차단, Tact/Th17↓, 연골 파괴 속도 감소
6. **Dapsone 100mg/d**: 경증-중등도 RP, MPO/ROS 억제, 호중구 매개 조직 손상↓
7. **Prednisone + Tocilizumab 병용**: 가장 강력한 항염 효과, RPDAI 최저

### Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 질환 개요, RPDAI 구성, 역학, 기관 침범 패턴, 예후인자 |
| 2. 약물 PK | Prednisolone 2CMT, GR 점유율, Tocilizumab 2CMT 혈중농도 시계열 |
| 3. 면역학·사이토카인 | T세포(Th1/Th17/Treg), B세포·Anti-CII Ab, 사이토카인(TNF/IL-6/IL-17/IL-1β), 보체 |
| 4. 연골 & RPDAI | 연골 무결성 (0-100%), RPDAI 대리지표, MMP 활성도, 기관별 위험 점수 |
| 5. 시나리오 비교 | 7개 치료군 동시 비교 (CartPct, RPDAI, CRP, Anti-CII Ab 등 선택 가능) |
| 6. 바이오마커 | CRP, Anti-CII Ab, IL-6/TNF/IL-17/IL-1β, IC, 보체, MMP 시계열 + 주요 시점 표 |
| 7. 민감도 분석 | 10개 파라미터에 대한 단변량 민감도 (연골무결성/RPDAI/Anti-CII Ab 평가) |

### 주요 보정 기준 (Calibration Targets)

- Mathian et al. 2019: Tocilizumab 치료 9예, 관해 유지 중앙값 24개월 (PMID: 30642828)
- Arnaud et al. 2012: RPDAI 타당성 검증, 5점 이상 활성 질환 — PMID 22072555
- Shimizu et al. 2019: TCZ 투여 후 CRP 정상화 및 재발 빈도 감소 (PMID: 30620291)
- Dion et al. 2007: 예후인자: 기도 침범, 심혈관 침범이 사망률과 연관

