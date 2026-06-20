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
| 2026-06-20 | 만성질환/결정관절염 | 가성통풍 (Pseudogout / CPPD Crystal Deposition Disease) | [pseudogout/](pseudogout/) | [![CPPD QSP](pseudogout/cppd_qsp_model.png)](pseudogout/cppd_qsp_model.svg) |
| 2026-06-20 | 만성질환/내분비 | 만성 갑상선 기능 저하증 (Chronic Hypothyroidism) | [chronic-hypothyroidism/](chronic-hypothyroidism/) | [![HYPO QSP](chronic-hypothyroidism/hypo_qsp_model.png)](chronic-hypothyroidism/hypo_qsp_model.svg) |
| 2026-06-20 | 만성질환/소화기 | 만성 위염 (Chronic Gastritis – H. pylori & Correa Cascade) | [chronic-gastritis/](chronic-gastritis/) | [![CGAST QSP](chronic-gastritis/cgast_qsp_model.png)](chronic-gastritis/cgast_qsp_model.svg) |
| 2026-06-20 | 자가면역질환/혈관염 | 결절성 다발동맥염 (Polyarteritis Nodosa, PAN) | [polyarteritis-nodosa/](polyarteritis-nodosa/) | [![PAN QSP](polyarteritis-nodosa/pan_qsp_model.png)](polyarteritis-nodosa/pan_qsp_model.svg) |
| 2026-06-20 | 자가면역질환/혈액 | 악성 빈혈 (Pernicious Anemia, PA) | [pernicious-anemia/](pernicious-anemia/) | [![PA QSP](pernicious-anemia/pna_qsp_model.png)](pernicious-anemia/pna_qsp_model.svg) |
| 2026-06-20 | 자가면역질환/반응성관절염 | 반응성 관절염 (Reactive Arthritis, ReA) | [reactive-arthritis/](reactive-arthritis/) | [![ReA QSP](reactive-arthritis/rea_qsp_model.png)](reactive-arthritis/rea_qsp_model.svg) |
| 2026-06-20 | 만성질환/직업성폐질환 | 진폐증 (Pneumoconiosis – Silicosis / CWP / Asbestosis) | [pneumoconiosis/](pneumoconiosis/) | [![PNM QSP](pneumoconiosis/pnm_qsp_model.png)](pneumoconiosis/pnm_qsp_model.svg) | 
| 2026-06-20 | 자가면역질환 / 혈액학 | 에반스 증후군 (Evans Syndrome – AIHA+ITP) | [evans-syndrome/](evans-syndrome/) | [![ES QSP](evans-syndrome/es_qsp_model.png)](evans-syndrome/es_qsp_model.svg) |
| 2026-06-20 | 만성질환/소화기 | 게실병 (Diverticular Disease — Diverticulosis / Diverticulitis) | [diverticular-disease/](diverticular-disease/) | [![DIV QSP](diverticular-disease/div_qsp_model.png)](diverticular-disease/div_qsp_model.svg) |
| 2026-06-20 | 자가면역질환/다발내분비 | 자가면역 다발성 내분비병증 (APS Type 1 / APECED — AIRE Mutation) | [autoimmune-polyendocrinopathy/](autoimmune-polyendocrinopathy/) | [![APS QSP](autoimmune-polyendocrinopathy/aps_qsp_model.png)](autoimmune-polyendocrinopathy/aps_qsp_model.svg) |
| 2026-06-20 | 자가면역질환/뇌하수체 | 림프구성 뇌하수체염 (Lymphocytic Hypophysitis, LyH) | [lymphocytic-hypophysitis/](lymphocytic-hypophysitis/) | [![LyH QSP](lymphocytic-hypophysitis/lhyp_qsp_model.png)](lymphocytic-hypophysitis/lhyp_qsp_model.svg) |
| 2026-06-20 | 만성질환/신장-요로계 | 만성 신우신염 (Chronic Pyelonephritis, CPN) | [chronic-pyelonephritis/](chronic-pyelonephritis/) | [![CPN QSP](chronic-pyelonephritis/cpn_qsp_model.png)](chronic-pyelonephritis/cpn_qsp_model.svg) |
| 2026-06-20 | 만성질환/대사 | 무증상 고요산혈증 (Asymptomatic Hyperuricemia, AHY) | [asymptomatic-hyperuricemia/](asymptomatic-hyperuricemia/) | [![AHY QSP](asymptomatic-hyperuricemia/ahy_qsp_model.png)](asymptomatic-hyperuricemia/ahy_qsp_model.svg) |
| 2026-06-20 | 만성질환/부인과 | 자궁내막증 (Endometriosis) | [endometriosis/](endometriosis/) | [![ENDO QSP](endometriosis/endo_qsp_model.png)](endometriosis/endo_qsp_model.svg) |
| 2026-06-20 | 신경퇴행성질환 | 알츠하이머병 (Alzheimer's Disease, AD) | [alzheimers-disease/](alzheimers-disease/) | [![AD QSP](alzheimers-disease/ad_qsp_model.png)](alzheimers-disease/ad_qsp_model.svg) |
| 2026-06-20 | 만성질환/피부/알레르기 | 아토피 피부염 (Atopic Dermatitis, AD) | [atopic-dermatitis/](atopic-dermatitis/) | [![ATAD QSP](atopic-dermatitis/ad_qsp_model.png)](atopic-dermatitis/ad_qsp_model.svg) |
| 2026-06-20 | 신경퇴행성질환 | 파킨슨병 (Parkinson's Disease, PD) | [parkinsons-disease/](parkinsons-disease/) | [![PD QSP](parkinsons-disease/pd_qsp_model.png)](parkinsons-disease/pd_qsp_model.svg) |
| 2026-06-20 | 신경정신질환 | 주요우울장애 (Major Depressive Disorder, MDD) | [major-depressive-disorder/](major-depressive-disorder/) | [![MDD QSP](major-depressive-disorder/mdd_qsp_model.png)](major-depressive-disorder/mdd_qsp_model.svg) |
| 2026-06-20 | 종양혈액학 | 만성 골수성 백혈병 (Chronic Myeloid Leukemia, CML) — 첫 번째 종양 모델 | [chronic-myeloid-leukemia/](chronic-myeloid-leukemia/) | [![CML QSP](chronic-myeloid-leukemia/cml_qsp_model.png)](chronic-myeloid-leukemia/cml_qsp_model.svg) |
| 2026-06-20 | 신경학 / 두통 | 편두통 (Migraine — Trigeminovascular, CGRP, CSD, Triptans/mAbs/Gepants) | [migraine/](migraine/) | [![MGR QSP](migraine/mgr_qsp_model.png)](migraine/mgr_qsp_model.svg) |
| 2026-06-20 | 만성 자가면역 피부질환 | 건선 (Psoriasis — IL-23/IL-17 Axis · TNF · Biologics/JAKi/PDE4i) | [psoriasis/](psoriasis/) | [![PSO QSP](psoriasis/pso_qsp_model.png)](psoriasis/pso_qsp_model.svg) |

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


---

## 가성통풍 (Pseudogout / CPPD Crystal Deposition Disease)

### 개요

가성통풍(CPPD 결정 침착 질환)은 칼슘 피로인산염(Ca₂P₂O₇, CPPD) 결정이 관절 연골·섬유연골·활막에 침착되어 반복성 급성 관절염과 만성 관절병증을 유발하는 결정성 관절염입니다. 유병률은 60세 이상 성인의 4-7%, 80세 이상 사망 후 부검에서 최대 35%에서 발견됩니다. 핵심 병태생리는 **세포외 무기 피로인산염(PPi) 과다 축적** → CPPD 결정 핵형성 → 결정 탈락(shedding) → **NLRP3 인플라마솜 활성화** → IL-1β 과다 분비 → 급성 관절염입니다.

### 병태생리 핵심 경로

| 경로 | 핵심 분자 | 임상 결과 |
|------|----------|---------|
| PPi 생성 | ENPP1/NPP1 (ATP→AMP+PPi) | 세포외 PPi ↑ (>6 µM) |
| PPi 유출 | ANKH 채널 (세포내→세포외) | 관절 내 PPi 축적 |
| PPi 가수분해 부전 | TNAP/ALPL↓ (PPi→2Pi 차단) | Fe²⁺/Cu²⁺/ALPL 돌연변이 |
| 결정 핵형성 | PPi + Ca²⁺ + Mg²⁺↓ → CPPD | 섬유연골 석회화 |
| 결정 탈락 | 외상/수술 → SF 내 결정 방출 | 급성 발작 유발 |
| NLRP3 활성화 | 카텝신 B 유출 → NLRP3-ASC-Casp-1 | IL-1β 성숙·분비 |
| 호중구 동원 | IL-8/CXCL8, C5a, LTB4 → SF 내 WBC↑ | SF WBC >50,000/µL |
| MMP/ADAMTS | IL-1β → MMP-1/3/13, ADAMTS-4/5 | 연골 기질 파괴 |

### 약물 PK/PD 핵심 파라미터

| 약물 | 용량/경로 | 반감기 | 표적 | 주요 효과 |
|------|---------|--------|------|---------|
| Colchicine | 0.5mg BID PO (급성/예방) | 26-31h (F=45%) | 튜불린 중합·NLRP3-ASC | 호중구 이동↓, IL-1β 분비↓ |
| Indomethacin | 50mg TID PO × 7일 | 4.5h (F=98%) | COX-1/2 비선택적 억제 | PGE2↓, 통증·부종 개선 |
| Naproxen | 500mg BID PO | 12-17h | COX-1/2 | GI 부작용 적음, NSAID 대안 |
| Prednisolone | 30-40mg/d PO → 테이퍼 | 2.5-4h (F=82%) | GR-α (NF-κB transrepression) | 다중 사이토카인↓, 신기능 저하 환자 선택 |
| IA Triamcinolone | 40mg 관절강 내 주사 | 국소 효과 지속 2-4주 | 국소 GR | 1개 관절 침범 시 1차 선택 |
| Anakinra | 100mg/d SC × 3-5일 | 4-6h (SC F=95%) | IL-1R 경쟁적 길항 | 불응성 또는 NSAIDs 금기 시 |
| Canakinumab | 150mg SC q3m | 26일 (항IL-1β mAb) | IL-1β 중화 | 반복 발작 예방, 고비용 |
| Hydroxychloroquine | 200-400mg/d PO (만성 예방) | 40-50일 | 리소솜 pH↑ → NLRP3↓ | 제한적 증거 |
| Methotrexate | 7.5-15mg/wk PO/SC (만성) | 3-10h (MTX-PG 수일) | ATIC/DHFR → adenosine↑ | 만성 CPPD 관절병증 예방 |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [![CPPD QSP](pseudogout/cppd_qsp_model.png)](pseudogout/cppd_qsp_model.svg) | 기계론적 지도 미리보기 (클릭 시 SVG 확대) |
| [cppd_qsp_model.dot](pseudogout/cppd_qsp_model.dot) | Graphviz 기계론적 지도 (12 클러스터, 140+ 노드) |
| [cppd_qsp_model.svg](pseudogout/cppd_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [cppd_qsp_model.png](pseudogout/cppd_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [cppd_mrgsolve_model.R](pseudogout/cppd_mrgsolve_model.R) | mrgsolve ODE 모델 (20 구획, 7 치료 시나리오) |
| [cppd_shiny_app.R](pseudogout/cppd_shiny_app.R) | Shiny 대시보드 (7탭) |
| [cppd_references.md](pseudogout/cppd_references.md) | 참고문헌 62편 (12섹션, PubMed 링크) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **미치료 급성 발작**: CPPD 결정 탈락 → NLRP3 최대 활성 → IL-1β 급등 → SF WBC 50,000/µL → 자연 소실 7-10일
2. **Colchicine 0.5mg BID**: NLRP3-ASC speck 억제 + 호중구 이동 차단 → 3일 내 통증 VAS 50% 감소
3. **Indomethacin 50mg TID**: COX-2 90% 억제 → PGE2↓ → 통증·발열 신속 개선
4. **Prednisolone 30mg/d (테이퍼)**: GR 점유율 88%+ → NF-κB transrepression → 모든 사이토카인 동반 감소
5. **Anakinra 100mg/d SC**: IL-1R 차단(KD≈0.5µg/mL) → IL-1β 증폭 루프 차단 → 불응성 CPPD에 유효
6. **Colchicine + Indomethacin 병용**: NLRP3·호중구(Colch) + PGE2·통증(Indo) 이중 표적
7. **Colchicine 0.5mg/day 장기 예방 (180일)**: 연간 발작 빈도↓, 연골 무결성 보전, 결정 부담 안정화

### Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 역학, 위험인자, 임상 양상, EULAR 진단기준, PPi 대사 요약표 |
| 2. 약물 PK | 콜히친/인도메타신/프레드니솔론/아나킨라 혈중농도 · GR/NLRP3/COX/IL-1R 점유율 |
| 3. 염증 패널 | IL-1β, NLRP3 활성, IL-6, PGE2, SF 호중구, Lipoxin A4 시계열 |
| 4. 결정 동역학 | PPi 세포외 농도, 연골 내 CPPD 부담, SF 내 CPPD 결정, 핵형성 해설 |
| 5. 임상 엔드포인트 | 통증 VAS, CRP, 연골 무결성(%), 주요 시점 요약 표 |
| 6. 시나리오 비교 | 6개 치료군 동시 비교 (IL-1β/통증/CRP/NLRP3 선택) |
| 7. 바이오마커 & 민감도 | 시점별 바이오마커 표, 10개 파라미터 단변량 민감도 분석 (Day 7 IL-1β) |

### 주요 보정 기준 (Calibration Targets)

- Martinon et al. 2006 (*Nature*): CPPD 결정이 NLRP3 인플라마솜 활성화 → IL-1β 분비 — PMID 16407889
- Announ et al. 2009: 말기 신기능 저하 환자에서 아나킨라로 가성통풍 발작 관해 — PMID 19345128
- Nuki 2008: 콜히친 0.5mg BID로 3일 내 통증 50% 이상 감소 — PMID 18638431
- Zhang et al. 2011 (EULAR Part I/II): CPPD 진단 및 치료 권고 — PMID 21216756/21216757
- Pascart & Richette 2019: CPPD 관리 업데이트, HCQ/MTX 만성 예방 근거 — PMID 31621572

---

## 만성 갑상선 기능 저하증 (Chronic Hypothyroidism)

### 개요

만성 갑상선 기능 저하증은 전 세계 성인의 약 2–4%(여성 5–7%, 남성 1–2%)에서 발생하는 흔한 내분비 질환입니다. 갑상선에서 티록신(T4) 및 트리요오드티로닌(T3) 합성·분비가 지속적으로 부족하여 다장기 대사 저하, 심혈관 효과, 신경정신계 증상이 나타납니다. 주요 원인은 하시모토 갑상선염(자가면역), 갑상선 절제술 후, 방사성 요오드(I-131) 치료 후, 약물 유발(아미오다론·리튬·IFN-α) 등입니다.

### 주요 병태생리 경로

| 경로 | 핵심 메커니즘 | 임상 이상 |
|------|------------|---------|
| HPT 축 피드백 | fT4/fT3 음성 피드백 소실 → TSH ↑↑ | 혈청 TSH >4.5 mIU/L (오현성 갑저: >10) |
| T4 합성 장애 | NIS/TPO 결핍, Tg 파괴 → T4 분비 ↓ | 총 T4 <58 nmol/L, fT4 <12 pmol/L |
| 말초 대사 이상 | D1·D2 활성↓ → T4→T3 변환↓ | fT3 <3.5 pmol/L, rT3 상대적 ↑ |
| 심혈관 효과 | T3 ↓ → β-AR 발현↓, MHC-β 전환 | 서맥, 심박출량↓, LDL↑, 이완기 고혈압 |
| 대사 저하 | UCP 발현↓ → BMR ↓15-30% | 체중 증가 2-10 kg, 추위 불내성 |
| 신경정신계 | 세로토닌·도파민 대사 ↓ | 피로, 우울증, 인지 저하, 건반사 이완 지연 |
| 점액수종 | GAG(글리코사미노글리칸) 피하 축적 | 비함요부종, 안주위 부종, 거대 혀 |
| 생식 이상 | TRH↑ → 프로락틴↑ → 무배란 | 월경 이상, 불임, 습관성 유산 |

### 약물 PK/PD 핵심 파라미터

| 약물 | 용량/경로 | PK 파라미터 | 임상 효과 |
|------|---------|------------|---------|
| 레보티록신 (LT4) | 25–200 μg/일 경구 (표준: 1.6 μg/kg/day) | F=0.65–0.80(공복), Ka=0.35 h⁻¹, Vss=700L, t½=7일 | TSH 정상화 6-8주, T4 풀 보충 |
| 리오티로닌 (LT3) | 5–25 μg/일 (2회 분할) | F=0.95, Ka=0.80 h⁻¹, V1=15L, t½=24h | 즉각적 T3 상승, 증상 잔존 환자 보조 |
| 탄산칼슘 (병용) | — | LT4 흡수 20–40% 감소 | 복용 간격 ≥4시간 필요 |
| 황산철 (병용) | — | LT4 흡수 30–40% 감소 | 복용 간격 ≥4시간 필요 |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [![HYPO QSP](chronic-hypothyroidism/hypo_qsp_model.png)](chronic-hypothyroidism/hypo_qsp_model.svg) | 기계론적 지도 미리보기 (클릭 시 SVG 확대) |
| [hypo_qsp_model.dot](chronic-hypothyroidism/hypo_qsp_model.dot) | Graphviz 기계론적 지도 (13 클러스터, 148+ 노드) |
| [hypo_qsp_model.svg](chronic-hypothyroidism/hypo_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [hypo_qsp_model.png](chronic-hypothyroidism/hypo_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [hypo_mrgsolve_model.R](chronic-hypothyroidism/hypo_mrgsolve_model.R) | mrgsolve ODE 모델 (15 구획, 7 치료 시나리오) |
| [hypo_shiny_app.R](chronic-hypothyroidism/hypo_shiny_app.R) | Shiny 대시보드 (6탭) |
| [hypo_references.md](chronic-hypothyroidism/hypo_references.md) | 참고문헌 40편 (10섹션, PubMed 링크) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **정상 갑상선 (Euthyroid 기준)**: HPT 축 정상 상태 — TSH 2.0 mIU/L, fT4 20 pmol/L, fT3 5.4 pmol/L
2. **미치료 갑저 (갑상선 전절제)**: TSH 급등 >30 mIU/L, T4/T3 소진, 모든 PD 효과 최악 (LDL↑, BMR↓, 증상↑)
3. **LT4 100 μg/day 표준 치료**: 6-8주 후 TSH 정상화, fT4 목표 범위 달성
4. **LT4 50 μg/day (준임상적 갑저)**: 갑상선 잔존 30% + 저용량 LT4 → TSH 정상 상단, FT4 낮은 정상
5. **LT4 175 μg/day 과치료 (TSH 억제)**: TSH <0.1 mIU/L, fT4 과잉 → 심방세동·골밀도 감소 위험
6. **T4+T3 병용 (LT4 100 μg + LT3 10 μg/day)**: LT3의 즉각적 T3 상승, 증상 잔존 환자 추가 이점
7. **갑상선 전절제 후 LT4 125 μg/day**: 체중 기반(1.8 μg/kg, 70 kg) 완전 대체 시뮬레이션

### Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 질환 유형(6종), 초기 검사 결과 입력, 병태생리 요약, 초기 TSH/FT4/FT3 신호등 |
| 2. HPT 축 | TRH·TSH·FT4·FT3·rT3 시계열, TSH 목표 달성 시점 예측, HPT 정상화 표 |
| 3. 약물 PK | LT4·LT3 PK 프로파일, 72시간 초기 흡수, 약물 상호작용 테이블 |
| 4. 장기 효과 | 심박수, LDL, BMR, 증상 점수, 체중 변화, 골밀도 변화 추적 |
| 5. 시나리오 비교 | 7가지 치료 시나리오 동시 비교 (TSH/FT4) + 정상상태 요약 표 |
| 6. 바이오마커·용량 최적화 | 용량-TSH 반응 곡선, 목표 TSH 기반 최적 용량 계산기, 모니터링 일정 |

### 주요 보정 기준 (Calibration Targets)

- Berberich et al. 2017 (*Thyroid*): 수학 모델 — Hill 피드백 계수 및 HPT 축 SS 조건 보정 — PMID 29095128
- Jonklaas et al. 2014 (*Thyroid*, ATA 가이드라인): LT4 표준 용량 1.6 μg/kg/day, SS 달성 6-8주 — PMID 25266247
- Bianco & Kim 2006 (*J Clin Invest*): D1/D2/D3 탈요오드화효소 상대 활성 비율 — PMID 17016550
- Klein & Ojamaa 2001 (*NEJM*): 갑저 심혈관 효과 — 심박수·심박출량·LDL 보정 기준 — PMID 11172193
- Pearce et al. 2013 (ETA 가이드라인): 준임상적 갑저 치료 결정 및 TSH 목표 범위 — PMID 24783053

---

## 만성 위염 (Chronic Gastritis – H. pylori Pathogenesis & Correa Cascade)

### 개요

만성 위염은 *Helicobacter pylori* (Hp) 감염을 주원인으로 하며, 전 세계 성인의 약 44%가 감염되어 있습니다. 감염자의 10–15%는 소화성 궤양, 1–3%는 위암으로 진행하며, Correa cascade 개념에 따라 **정상 점막 → 표재성 위염 → 위축성 위염 → 장상피화생 → 이형성 → 위선암**의 단계적 경로를 밟습니다. 본 QSP 모델은 Hp 세균 역학, 독성인자(CagA/VacA/OipA), 면역 반응 cascade, 위산 분비 생리, Correa 진행 스코어, 제균요법 PK/PD를 통합 시뮬레이션합니다.

### 주요 병태생리 경로

| 경로 | 핵심 메커니즘 | 임상 이상 |
|------|------------|---------|
| Hp 점막 집락 | 요소분해효소(암모니아) → pH 완충 + 편모 운동 → 점액층 침투 → BabA/HopQ 접착 | 점막 밀도 10³–10⁸ CFU/생검 |
| CagA 독성 | T4SS 주입 → EPIYA 인산화 → SHP2-CagA 복합체 → MAPK/ERK → 세포 극성 소실 | 위선암 위험 2-3배 ↑ |
| VacA 독성 | s1/m1 균주 → 공포화 + 미토콘드리아 손상 → ROS 생성 + NLRP3 인플라마솜 → IL-1β | 상피세포 아폽토시스 |
| 선천면역 활성화 | TLR4/5/9, NOD1 → NF-κB → IL-8(호중구 모집), IL-1β, TNF-α, IL-6 | 급성 위염 → 만성화 |
| 적응면역 | Th1(IFN-γ) 우세 → 위산 분비 억제 (IL-1β 가장 강력); Treg ↓ 시 Th17 증가 | 위축성 진행 가속 |
| 위산 분비 이상 | G세포(가스트린)↑, D세포(소마토스타틴)↓ → H⁺/K⁺-ATPase 과활성 → 위산↑(초기) → 위축 시↓ | 십이지장궤양(초기) vs. 위암(말기) |
| Correa cascade | IL-1β/TNF-α/IFN-γ → OLGA 위축 + E-cadherin 소실 + TP53 변이 → 선암 | OLGA III-IV: 1.2-2.5%/년 |

### 약물 PK/PD 핵심 파라미터

| 약물 | 용량/경로 | PK 파라미터 | 임상 효과 |
|------|---------|------------|---------|
| 오메프라졸 (PPI) | 20 mg BID 경구 | F=0.65, Ka=0.50/h, t½=1.5h | 3-5일 후 최대 산분비 억제, 위내 pH↑ |
| 아목시실린 | 1000 mg BID 경구 | F=0.90, Ka=1.0/h, t½=1.2h | MIC 0.016-0.5 mg/L vs Hp, pH 의존성 낮음 |
| 클래리스로마이신 | 500 mg BID 경구 | F=0.55, Ka=0.6/h, t½=7h | MIC 0.016-64 mg/L, 23S rRNA 결합, 내성률 상승 중 |
| 메트로니다졸 | 500 mg TID 경구 | F=0.99, Ka=0.8/h, t½=6h | rdxA nitroreduction → DNA 이중가닥 절단 |
| 비스무스 (BSS) | 525 mg QID 경구 | 국소 작용 | Hp 접착 차단 + 세포막 손상 |
| 보노프라잔 (VPZ) | 20 mg QD 경구 | F=0.80, Ka=1.2/h, t½=7h | K⁺-경쟁적 산분비 차단, 첫날부터 완전 억제 |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [![CGAST QSP](chronic-gastritis/cgast_qsp_model.png)](chronic-gastritis/cgast_qsp_model.svg) | 기계론적 지도 미리보기 (클릭 시 SVG 확대) |
| [cgast_qsp_model.dot](chronic-gastritis/cgast_qsp_model.dot) | Graphviz 기계론적 지도 (10 클러스터, 153+ 노드) |
| [cgast_qsp_model.svg](chronic-gastritis/cgast_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [cgast_qsp_model.png](chronic-gastritis/cgast_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [cgast_mrgsolve_model.R](chronic-gastritis/cgast_mrgsolve_model.R) | mrgsolve ODE 모델 (22 구획, 7 치료 시나리오) |
| [cgast_shiny_app.R](chronic-gastritis/cgast_shiny_app.R) | Shiny 대시보드 (7탭) |
| [cgast_references.md](chronic-gastritis/cgast_references.md) | 참고문헌 60편 (12섹션, PubMed 링크) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **무치료 (자연경과)**: Hp 지속 집락 → NF-κB/IL-8 만성 상승 → 위축 점진 진행
2. **PPI 단독 (오메프라졸 20mg BID × 8주)**: 산분비 억제 → 증상 완화 but Hp 지속
3. **표준 삼제 요법 (PPI+AMX+CLR × 14일)**: ITT 제균율 80-85%, CLR 내성 시 급격히 감소
4. **비스무스 사제 요법 (PPI+AMX+CLR+BSS × 14일)**: ITT 85-92%, CLR 내성 지역 우선
5. **메트로니다졸 사제 요법 (PPI+AMX+MTZ+BSS × 14일)**: CLR 내성 지역 대안
6. **보노프라잔 삼제 (VPZ+AMX+CLR × 14일)**: 첫날부터 완전 산억제 → AMX 안정성↑, ITT 88-94%
7. **제균 후 5년 추적 (Correa cascade 역행)**: 위축 점수 서서히 감소, IM은 부분 역행, 위암 위험 감소

### Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 질환 개요, Correa 단계 표, 독성인자 설명, 치료 추천, 위암 위험도 평가 |
| 2. H. pylori 역학 | Hp 점막 밀도 시계열, 제균 임계치 표시, 제균율 표, 약물 농도 그래프 |
| 3. 염증 반응 | IL-8·IL-1β·TNF-α·IFN-γ·IL-10 사이토카인 동태, 면역세포(호중구·Th1·Treg) |
| 4. 위산 생리 | 위산 출력(mmol/h), 혈청 가스트린(G17), 점액층 두께(μm), PPI PD |
| 5. Correa Cascade | 위축·장상피화생 점수 시계열, OLGA 단계, PGI/PGII 비율 추적 |
| 6. 치료 비교 | 5가지 요법 동시 비교 (Hp 제균, 증상 점수), 요약 표 |
| 7. 바이오마커 & 임상 지표 | GastroPanel®(PGI+PGII+G17+IgG), 제균 검사(¹³C-UBT/HpSA), 기준치 표 |

### 주요 보정 기준 (Calibration Targets)

- Malfertheiner et al. 2022 (*Gut*, Maastricht VI): 표준 삼제 ITT 80-85%, 비스무스 사제 85-92% — PMID 35469816
- Graham & Dore 2020 (*Expert Rev Gastroenterol Hepatol*): VPZ 삼제 ITT 88-94% — PMID 32356979
- El-Omar et al. 2000 (*Nature*): IL-1β 다형성 → 위산 억제 + 위암 위험 → IL-1β 용량반응 보정 — PMID 10746728
- Correa 1992 (*Cancer Res*): Correa cascade 단계별 진행 속도 — PMID 1612357
- Rugge et al. 2010 (*Aliment Pharmacol Ther*): OLGA III-IV: 1.2-2.5%/년 위암 위험 — PMID 20199491

---

## 결절성 다발동맥염 (Polyarteritis Nodosa, PAN)

### 개요

결절성 다발동맥염(PAN)은 중간 크기 동맥의 전층에 걸친 괴사성 염증으로 특징지어지는 전신 혈관염입니다. 연간 발생률은 100만 명당 약 2–9명이며, 30–40%는 B형 간염 바이러스(HBV)와 연관됩니다. 핵심 병태생리는 면역복합체(HBV 항원–항체 복합체 등)의 혈관벽 침착 → 보체 활성화 → 호중구 동원 → 혈관 전층 괴사성 염증 → 미세동맥류 형성 → 장기 허혈의 연쇄 반응입니다. 주요 침범 장기는 신장(신동맥 미세동맥류), 말초신경(단발신경염), 위장관(장간막 허혈), 피부(그물 모양 청피증, 결절), 고환(고환통) 등입니다. ANCA는 일반적으로 음성이며, 이것이 ANCA 연관 혈관염(GPA, MPA, EGPA)과의 핵심 감별 포인트입니다.

### 주요 병태생리 경로

| 경로 | 핵심 메커니즘 | 임상 결과 |
|------|------------|---------|
| 면역복합체 형성 | HBV 표면항원 + IgG/IgM → 면역복합체 → 혈관벽 침착 | 보체 활성화, 혈관 염증 유발 |
| 보체 고전경로 | C1q/C4/C3 활성화 → C5a(강력한 화학유인물질) + MAC | 호중구 대량 동원 |
| 호중구 침윤 | C5a → 호중구 혈관벽 이동 → MPO/elastase/ROS 방출 → NET 형성 | 혈관내피 괴사(fibrinoid necrosis) |
| 전층 동맥 염증 | Th1/CD4+ T세포 + 대식세포 → TNF-α/IL-6/IL-1β | 혈관 전층(transmural) 파괴 |
| 미세동맥류 | 혈관벽 약화 → 국소 확장 → 파열 위험 | 내장 출혈, 신장 경색 |
| 혈전 형성 | 내피세포 손상 → vWF↑/TF↑/TXA₂↑ → 혈소판 활성화 | 허혈성 장기 손상 |
| 섬유화 | TGF-β1 + 섬유아세포 → 콜라겐 침착 → 동맥 협착 | 만성 허혈, 신동맥 협착성 고혈압 |
| HBV-PAN 특이 | HBV 표면항원(HBsAg) 지속 → 면역복합체 지속 생성 | 항바이러스 치료로 관해 유도 가능 |

### 약물 PK/PD 핵심 파라미터

| 약물 | 용량/경로 | PK 파라미터 | 임상 효과 |
|------|---------|------------|---------|
| 프레드니솔론 (Prednisolone) | 1 mg/kg/day PO 4주 후 감량 | F=0.82, t½=3–4h, Vd=45L, CL=10.5L/h | NF-κB 억제 → 사이토카인↓ (EMAX=85%) |
| 사이클로포스파미드 정맥 펄스 (CYC IV) | 750 mg/m² q3w × 6회 (CYCLOPS 기반) | Vd=38L, CL=5.8L/h; 활성화: CYP2B6 → 4-OH-CYC | B세포/T세포 고갈, 면역억제 |
| 아자티오프린 (AZA) | 2 mg/kg/day PO (유지요법) | F=~0.88, Ka=0.9/h, CL=8.2L/h | 6-MP → 퓨린 합성 억제 → T세포↓ |
| 엔테카비어/테노포비어 (HBV-PAN) | 표준 항바이러스 용량 | 간 표적, 경구 | HBV DNA↓ → 면역복합체 생성↓ → 관해 |
| 혈장교환 (PE) | 3–4L × 6회 (HBV-PAN 초기) | 면역복합체 직접 제거 | 급성 관해 유도, 항바이러스와 병용 |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [![PAN QSP](polyarteritis-nodosa/pan_qsp_model.png)](polyarteritis-nodosa/pan_qsp_model.svg) | 기계론적 지도 미리보기 (클릭 시 SVG 확대) |
| [pan_qsp_model.dot](polyarteritis-nodosa/pan_qsp_model.dot) | Graphviz 기계론적 지도 (13 클러스터, 234+ 노드, 329 엣지) |
| [pan_qsp_model.svg](polyarteritis-nodosa/pan_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [pan_qsp_model.png](polyarteritis-nodosa/pan_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [pan_mrgsolve_model.R](polyarteritis-nodosa/pan_mrgsolve_model.R) | mrgsolve ODE 모델 (19 구획, 6 치료 시나리오) |
| [pan_shiny_app.R](polyarteritis-nodosa/pan_shiny_app.R) | Shiny 대시보드 (6탭: 환자프로파일·약물PK·PD바이오마커·임상엔드포인트·시나리오비교·바이오마커상관) |
| [pan_references.md](polyarteritis-nodosa/pan_references.md) | 참고문헌 44편 (9섹션, PubMed 링크) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **무치료 (자연경과)**: 면역복합체 지속 → 보체 활성화 → 혈관 염증 진행 → 장기 손상 축적
2. **고용량 프레드니솔론 단독 (1 mg/kg/day × 4주 후 감량)**: NF-κB 억제 → 사이토카인↓, 혈관 염증 완화; 2년 생존율 약 55% (1950~1960년대 데이터)
3. **프레드니솔론 + CYC IV 펄스 (CYCLOPS 기반)**: 림프구 고갈 → BVAS 급속 저하, 18개월 관해율 ~90%; CYCLOPS 시험 (de Groot 2009)
4. **프레드니솔론 + 아자티오프린 유지 (Pagnoux 2008 프로토콜)**: 관해 유도 후 AZA로 전환, 재발율 감소
5. **HBV-PAN: 항바이러스 + 저용량 프레드니솔론 + 혈장교환**: Guillevin 1995/2005 프로토콜 — HBsAg 소실율 50–67%, 면역억제제 최소화
6. **관해 유도 → 유지 순차 치료 (EULAR 2016 지침 기반)**: 고용량 Pred+CYC 유도 → 관해 후 AZA/MTX 유지, 24개월 추적

### Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 나이/성별/체중/HBV 상태/기저 BVAS/CRP/eGFR 입력, 장기 침범 체크박스, Five Factor Score 자동 계산, 위험 계층화 (저/중/고, 5년 사망률 추정) |
| 2. 약물 PK | 치료 선택(4가지 레짐), 용량/빈도/기간 슬라이더, 프레드니솔론·CYC·AZA 혈중 농도 시계열, 치료/독성 범위 표시, Cmax/Tmax/AUC 요약 |
| 3. PD 바이오마커 | EC50/Emax 민감도 슬라이더, 면역복합체·C3 보체·사이토카인(IL-6/CRP/ESR)·호중구/B세포/T세포 동태 시계열 |
| 4. 임상 엔드포인트 | BVAS(0–63, 색상 구역 표시), 혈관 염증 지수, eGFR(CKD 단계 구역), 신경 손상 점수(단발신경염), 미세동맥류 부담 지수 |
| 5. 시나리오 비교 | 5가지 치료군 동시 비교(무치료/Pred단독/Pred+CYC/Pred+AZA/HBV+Pred), BVAS·CRP·eGFR 시계열, 관해율 포레스트 플롯, 결과 요약 DT 표 |
| 6. 바이오마커 상관 | BVAS vs CRP 산점도, eGFR vs 면역복합체 산점도, 상관행렬 히트맵, 기저 대비 변화율 폭포 차트, 아형별 BVAS 궤적 |

### 주요 보정 기준 (Calibration Targets)

- Guillevin et al. 1995 (*Medicine*): HBV-PAN 항바이러스+혈장교환 전략 — HBsAg 소실율 50–67%, 사망률 14% → 6% — PMID 7776987
- de Groot et al. 2009 (*Ann Intern Med*, CYCLOPS): CYC IV 펄스 vs. 경구 — 관해율 동등(88% vs 82%), 백혈구감소증 부작용 감소 — PMID 19246992
- Guillevin et al. 2011 (*Arthritis Rheum*): 개정 Five Factor Score — 심각한 위장관·심근·중추신경 침범 시 FFS≥2 예측 인자 — PMID 21200183
- Samson et al. 2014 (*J Rheumatol*): PAN 장기 추적(10년) — 재발율 31%, 장기 손상 57%, 표준화 사망비 2.7 — PMID 24429165
- Pagnoux et al. 2008 (*NEJM*): 관해 후 MTX vs AZA 유지 — 18개월 재발율 동등(36% vs 33%) — PMID 19052126

---

## 악성 빈혈 (Pernicious Anemia, PA)

### 개요

악성 빈혈(PA)은 위 벽세포(parietal cell)와 내인성 인자(intrinsic factor, IF)에 대한 자가면역 공격으로 인해 코발라민(비타민 B12) 흡수가 차단되어 발생하는 자가면역 질환입니다. 진행된 위체부 위축(corpus gastric atrophy)으로 인해 거대아세포 빈혈(megaloblastic anemia)과 아급성 복합 변성(subacute combined degeneration, SCD) 등 신경학적 합병증이 유발됩니다.

| 항목 | 내용 |
|------|------|
| **병인** | 항위벽세포항체(Anti-PC Ab, ~90%) 및 항내인성인자항체(Anti-IF Ab Type I 70%, Type II 35%) |
| **유병률** | 일반인구의 ~0.1%; 60세 이상 여성에서 높음 |
| **주요 발현** | 거대아세포 빈혈 · 무위산증 · 신경병증 · SCD |
| **진단 바이오마커** | 혈청 B12 <200 pg/mL, HoloTC <35 pmol/L, MMA >0.4 µmol/L, Hcy >15 µmol/L |
| **치료** | 시아노코발라민/하이드록소코발라민 IM 주사 또는 고용량 경구 보충 |

### 기계론적 지도 (Mechanistic Map)

[![PA QSP Map](pernicious-anemia/pna_qsp_model.png)](pernicious-anemia/pna_qsp_model.svg)

- **10개 서브그래프 클러스터**: 유전/환경 위험인자 · 위자가면역 · 식이 B12 처리 · 회장 흡수 · 혈장 수송 · 세포내 대사 · 조혈 및 거대아세포 빈혈 · 신경계 합병증 · 약물 PK/PD · 임상 지표
- **105+ 노드**: 항체(Anti-PC, Anti-IF T1/T2), CD4/CD8 T세포, 큐빌린/암니오니스 수용체, 트랜스코발라민 I/II/III, 메틸코발라민, 아데노실코발라민, 메치오닌 합성효소, MMA, Hcy 등
- **시아노코발라민/하이드록소코발라민 PK**: IM depot → 혈장 분포 → 간 저장 → 골수/신경계 전달

### mrgsolve 모델 구성

| 구획 (Compartment) | 역할 |
|--------------------|------|
| DEPOT | IM 주사 depot (µg) |
| ORAL_GI | 경구 투여 B12 (GI 내) |
| IF_POOL | 기능적 내인성 인자 (0-1) |
| CBA_PORT | 문맥혈 코발라민 (µg) |
| PLASMA | 혈장 총 코발라민 (pg/mL) |
| HOLOTC | 홀로트랜스코발라민 II (pmol/L) |
| LIVER | 간 코발라민 저장 (µg) |
| BONE_MRW | 골수 코발라민 풀 |
| NERVE | 신경계 코발라민 풀 |
| AUTO_AB | 자가항체 부담 (0-1) |
| PARIETAL | 위 벽세포 기능 지수 (0-1) |
| HGB | 헤모글로빈 (g/dL) |
| MCV | 평균 적혈구 용적 (fL) |
| RETIC | 망상적혈구 (×10⁹/L) |
| NEURO | 신경학적 장애 점수 (0-10) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **무치료 (S1)**: 자가면역 진행 → IF 지속 감소 → B12 흡수 차단 → Hb 저하, MCV 상승, 신경 손상 누적
2. **IM 표준 (S2 — 시아노코발라민 1000 µg/day × 7일 후 월 1회 유지)**: 4–10일 내 망상적혈구 위기 → 6–8주 내 Hb 회복 → MCV 2–3개월 내 정상화
3. **IM 유지만 (S3 — 월 1회 1000 µg, 로딩 없음)**: Hb 회복 지연 (~10–12주); 로딩 없이도 1년 내 완전 관해 가능
4. **고용량 경구 B12 (S4 — 1000 µg/day, 수동흡수 1%)**: Kuzminski 1998 임상시험 보정; 90일 내 IM 동등 효과; 경구 순응도 우수
5. **적극적 IM (S5 — 2000 µg/day × 14일 → 격주 유지)**: 빠른 신경계 회복 목표; 심각한 SCD 환자에 권장

### Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 기저 Hb·MCV·신경 점수·벽세포 기능(PA 중증도) 입력, 위 흡수 효율 다이어그램 |
| 2. 약물 PK | 혈장 B12 시계열, HoloTC, 간 저장량; IM/경구 경로 선택, 용량/빈도 슬라이더 |
| 3. 혈액 PD | Hb, MCV 시계열; 망상적혈구 위기 (첫 30일 확대); 혈액학적 이정표 테이블 |
| 4. 신경계 PD | 신경학적 장애 점수 (0-10), SCD 단계 표시; 신경 B12 풀 플롯; 회복 평가 테이블 |
| 5. 바이오마커 | MMA·Hcy 동태; 벽세포 기능 및 자가항체 부담; 바이오마커 요약 DT 표 |
| 6. 시나리오 비교 | 5가지 치료군 동시 비교; 엔드포인트 선택; 6개월·12개월 결과 요약 |

### 주요 보정 기준 (Calibration Targets)

- Kuzminski et al. 1998 (*Blood*): 고용량 경구(2000 µg/day) vs. IM 1000 µg — 90일 MMA·Hcy·Hb 회복 동등 — PMID 9694707
- Andrès et al. 2004 (*Am J Med*): 망상적혈구 위기 Day 4–7(범위 3–10일), Hb 주당 ~1–2 g/dL 상승 — PMID 15235333
- Wolffenbuttel et al. 2019 (*Mayo Clin Proc*): HoloTC <35 pmol/L = 가장 조기 결핍 마커; 총 B12 <200 pg/mL = 후기 고갈 — PMID 31193809
- Healton et al. 1991 (*Medicine*): SCD 임상 특성; 신경학적 악화/회복 시간경과 — PMID 2072028
- Stabler 2013 (*NEJM*): MMA >0.4 µmol/L 및 tHcy >15 µmol/L를 기능적 결핍의 진단 기준으로 사용 — PMID 23301732

### 파일 목록

| 파일 | 설명 |
|------|------|
| [pna_qsp_model.dot](pernicious-anemia/pna_qsp_model.dot) | Graphviz 기계론적 지도 (105+ 노드, 10개 클러스터) |
| [pna_qsp_model.svg](pernicious-anemia/pna_qsp_model.svg) | SVG 벡터 이미지 |
| [pna_qsp_model.png](pernicious-anemia/pna_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [pna_mrgsolve_model.R](pernicious-anemia/pna_mrgsolve_model.R) | mrgsolve ODE 모델 (15 구획, 5 시나리오) |
| [pna_shiny_app.R](pernicious-anemia/pna_shiny_app.R) | Shiny 대시보드 (6탭: 환자프로파일·약물PK·혈액PD·신경PD·바이오마커·시나리오비교) |
| [pna_references.md](pernicious-anemia/pna_references.md) | 참고문헌 35편 (9섹션, PubMed 링크) |

---

## 반응성 관절염 (Reactive Arthritis, ReA)

### 개요

반응성 관절염(ReA)은 비뇨생식기(Chlamydia trachomatis) 또는 위장관(Salmonella, Campylobacter, Yersinia, Shigella) 감염 후 발생하는 무균성 염증성 관절염입니다. 감염 병원체가 관절 자체에 존재하지 않음에도 불구하고 분자 모방(molecular mimicry), HLA-B27 연관 T세포 활성화, Th17 축 과활성화에 의해 주로 하지 대관절을 침범하는 비대칭 올리고관절염이 발생합니다.

| 항목 | 내용 |
|------|------|
| **병인** | 감염 후 HLA-B27 연관 면역 반응; Chlamydia 활막 지속성; Th17/IL-17 축 |
| **삼징후(고전 트라이어드)** | 요도염 · 관절염 · 결막염 (Reiter 증후군) |
| **HLA-B27 양성률** | 환자의 ~75% (일반인 8%) |
| **유병률** | GI 감염 후 약 1–4%, GU Chlamydia 감염 후 약 1–3% |
| **만성화** | 약 15–30%에서 6개월 이상 지속 |
| **주요 치료** | NSAIDs → 항생제(Chlamydia) → DMARDs (SSZ, MTX) → 생물학적 제제 (ETN, SECU) |

### 기계론적 지도 (Mechanistic Map)

[![ReA QSP Map](reactive-arthritis/rea_qsp_model.png)](reactive-arthritis/rea_qsp_model.svg)

- **10개 서브그래프 클러스터**: 유발 감염 · 선천면역 · HLA-B27/분자모방 · 적응면역(T·B세포) · 활막관절 병리 · 관절외 증상 · 골/연골 개조 · NSAID PK/PD · DMARD/생물학적제제 PK/PD · 임상지표
- **145+ 노드**: LPS·TLR4/2/9·NLRP3·NF-κB·Th17·IL-17A·TNF-α·RANKL·OPG·MMP-3/9/13·VEGF·PGE2·Etanercept·Secukinumab 등
- **치료 경로**: NSAIDs(COX-2 억제), 설파살라진(NF-κB·IL-6 억제), MTX(AICAR·DHFR 경로), ETN(TNF-α 중화), Secukinumab(IL-17A 중화), 독시사이클린(Chlamydia 박멸)

### mrgsolve 모델 구성

| 구획 (Compartment) | 역할 |
|--------------------|------|
| PATH | 급성 병원균 부하 (정규화, 0–1) |
| CHRON_PATH | Chlamydia 활막 지속성 (만성화 인자) |
| INNATE | 선천면역 활성화 점수 |
| NEUT | 호중구 활성화 |
| MACRO | 대식세포 활성화 (M1) |
| TH1 | Th1 세포 (IFN-γ, TNF-α 생산) |
| TH17 | Th17 세포 (IL-17A, IL-22 생산) |
| TREG | 조절 T세포 (IL-10, TGF-β) |
| TNF_c | TNF-α (pg/mL) |
| IL17_c | IL-17A (pg/mL) |
| IL6_c | IL-6 (pg/mL) |
| IL10_c | IL-10 (pg/mL) |
| IFNg_c | IFN-γ (pg/mL) |
| SYNOV | 활막 염증 점수 (0–10) |
| CARTDMG | 누적 연골 손상 (0–1, 비가역적) |
| PAIN | VAS 통증 (0–100) |
| CRP_c | CRP (mg/L) |
| JCOUNT | 종창 관절 수 (0–28) |
| NSAID_DEPOT / NSAID_C | 나프록센 PK (Vd 0.12 L/kg, t½ 2.1h) |
| SSZ_DEPOT / SSZ_C | 설파살라진 PK (F 10%, t½ 21h) |
| MTX_DEPOT / MTX_C | 메토트렉세이트 PK (SC F 78%, t½ 10h) |
| TNFi_DEPOT / TNFi_C | 에타너셉트 PK (SC F 76%, t½ 125h) |
| IL17i_DEPOT / IL17i_C | 세쿠키누맙 PK (SC F 73%, t½ 27d) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **무치료 (S1)**: HLA-B27(+) Chlamydia ReA 자연 경과 → 12주 만성화 위험 예측; CRP·활막염 지속, 연골 손상 누적
2. **NSAIDs 단독 (S2 — 나프록센 500mg BID × 12주)**: 통증·활막 염증 부분 억제, CRP 부분 감소; 근본 면역 기전 미해결
3. **항생제 + NSAIDs (S3 — 독시사이클린 12주 + 나프록센)**: Carter 2010 임상시험 보정; Chlamydia 지속성 억제 → 관해율 향상
4. **NSAIDs + 설파살라진 (S4 — SSZ 2–3 g/day)**: Clegg 1996 RCT 보정; 만성 관절 증상에서 관절지수 26% 개선
5. **NSAIDs + SSZ → 에타너셉트 (S5 — 12주 후 ETN 50mg/wk)**: Braun 2014 근거; BASDAI·VAS 50% 이상 감소, 12주 ETN 시작 후 빠른 반응

### Shiny 앱 탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 나이/성별/체중/HLA-B27 상태/유발 감염 종류/기저 CRP·VAS·관절 수 입력, 만성화 위험 점수 자동 계산, 예후 예측 그래프 |
| 2. 약물 PK | NSAID·SSZ·MTX·에타너셉트·세쿠키누맙 용량/기간 슬라이더, 혈장 농도 시계열, IC50 기준선 표시, PK 지표 테이블 |
| 3. 질환 역학 | 면역세포(Th1/Th17/Treg/호중구) 동태, 사이토카인(TNF·IL-17·IL-6·IL-10·IFN-γ) 시계열, 병원균 제거 및 활막 염증 |
| 4. 임상 엔드포인트 | VAS 통증·CRP·종창 관절 수 시계열, 정상화 시점 자동 계산, 연골 손상 누적 그래프, 임상 결과 요약 테이블 |
| 5. 시나리오 비교 | 5가지 치료군 동시 비교, 엔드포인트 선택, 지정 주차 결과 막대 그래프, 주요 시점(4·12·24·52주) 결과 테이블 |
| 6. 바이오마커 패널 | HLA-B27 배율 민감도 분석, 사이토카인 상관관계 히트맵, IC50 민감도 분석(NSAID·SSZ·ETN), 만성화 위험 지표 |

### 주요 보정 기준 (Calibration Targets)

- Carter JD 2010 (*Ann Intern Med*): 항생제 병합(Doxy+Azithro) 6개월 후 관해율 ~40% vs. 위약 ~11% — PMID 20957185
- Clegg DO 1996 (*Arthritis Rheum*): SSZ 6개월 관절지수 개선 26% vs. 위약 9% (n=134) — PMID 8961909
- Braun J 2014 (*J Rheumatol*): ETN 12주 후 BASDAI <4 달성 ~60%, VAS 50% 감소 — PMID 25128038
- Baeten D 2015 (*NEJM*): Secukinumab AS ASAS20: 61% vs. 위약 28% (wk 16) — PMID 26699169

### 파일 목록

| 파일 | 설명 |
|------|------|
| [![ReA QSP](reactive-arthritis/rea_qsp_model.png)](reactive-arthritis/rea_qsp_model.svg) | 기계론적 지도 미리보기 (클릭 시 SVG 확대) |
| [rea_qsp_model.dot](reactive-arthritis/rea_qsp_model.dot) | Graphviz 기계론적 지도 (10 클러스터, 145+ 노드) |
| [rea_qsp_model.svg](reactive-arthritis/rea_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [rea_qsp_model.png](reactive-arthritis/rea_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [rea_mrgsolve_model.R](reactive-arthritis/rea_mrgsolve_model.R) | mrgsolve ODE 모델 (26 구획, 5 치료 시나리오) |
| [rea_shiny_app.R](reactive-arthritis/rea_shiny_app.R) | Shiny 대시보드 (6탭: 환자프로파일·약물PK·질환역학·임상엔드포인트·시나리오비교·바이오마커패널) |
| [rea_references.md](reactive-arthritis/rea_references.md) | 참고문헌 51편 (11섹션, PubMed 링크) |
---

## 진폐증 (Pneumoconiosis) — 상세 설명

### 질환 개요

진폐증(Pneumoconiosis)은 규소(silica), 석탄 분진, 석면(asbestos) 등 무기 광물성 분진의 장기 흡입으로 인해 폐 대식세포가 분진을 탐식하는 과정에서 만성 염증이 유발되고, 이것이 TGF-β1 매개 섬유화 경로를 활성화하여 폐 조직에 결절성·미만성 섬유화가 발생하는 직업성 만성 폐질환입니다.

주요 아형:
- **규폐증(Silicosis)**: 결정형 이산화규소(SiO₂) 흡입 → NLRP3 인플라마솜 활성화 → 규폐성 결절 → 진행성 대규모 섬유증(PMF)
- **탄광부 진폐증(CWP, Coal Workers' Pneumoconiosis)**: 석탄 분진 흡입 → 단순 탄광부 진폐증 → PMF 진행 가능
- **석면증(Asbestosis)**: 석면 섬유 → 좌절된 식세포 작용 → 광범위한 간질 섬유화 → 악성중피종 위험 증가

### 발병 기전 핵심 경로

```
분진 흡입 → 폐포 침착
    ↓
폐포 대식세포 탐식 → 용해소체 파열 → Cathepsin B 방출
    ↓
NLRP3 인플라마솜 활성화 (K⁺ 유출, ROS 자극)
    ↓
Caspase-1 → IL-1β, IL-18 성숙 → 파이롭토시스
    ↓
TNF-α, TGF-β1, IL-6 방출 (NF-κB 경로)
    ↓                        ↓
산화 스트레스 (ROS ↑)    TGF-β1 → Smad2/3 → 섬유아세포 증식
    ↓                        ↓
GSH 고갈               근섬유아세포 분화 (α-SMA+) → 콜라겐 과다 침착
    ↓                        ↓
DNA 손상, 세포사       규폐성 결절 → PMF (진행성 대규모 섬유증)
                             ↓
                    FVC↓, DLCO↓, TLC↓ → 폐동맥 고혈압 → 폐심장증 → 사망
```

### QSP 모델 특징

| 항목 | 내용 |
|------|------|
| **기계론적 지도** | 11개 서브그래프 클러스터, 120+ 노드 |
| **ODE 구획 수** | 20개 (약물 PK 4종 + 질환 PD 16개) |
| **치료 시나리오** | 6개 (자연경과, NAC, 피르페니돈, 닌테다닙, 분진 제거+NAC, 병합요법) |
| **바이오마커** | KL-6, SP-D, FVC, DLCO, mPAP, PVR |
| **참고문헌** | 52편 (PubMed 링크 포함) |

### 치료 목표 및 약물 PK/PD

| 약물 | 기전 | 용량 | ka (/h) | CL (L/h) | Vd (L) | Emax |
|------|------|------|---------|---------|--------|------|
| N-Acetylcysteine (NAC) | 항산화 (GSH 보충) | 600 mg TID | 1.5 | 15.0 | 25 | 0.80 |
| Pirfenidone (피르페니돈) | TGF-β1↓, 콜라겐↓ | 2403 mg/d | 2.0 | 3.5 | 70 | 0.75 |
| Nintedanib (닌테다닙) | PDGFR/FGFR/VEGFR 억제 | 300 mg/d | 0.8 | 1.3 | 640 | 0.70 |
| Tetrandrine (테트란드린) | NF-κB↓, NLRP3↓ | 120 mg/d | 1.2 | 8.0 | 120 | 0.65 |

### Shiny 앱 탭 구성

1. **환자 프로파일** — 나이, 직업 노출력, 흡연력, 기저 FVC, 질환 아형 설정
2. **약물 PK** — NAC·피르페니돈·닌테다닙·테트란드린 혈장 농도 프로파일
3. **PD 바이오마커** — 사이토카인(IL-1β, TNF-α, TGF-β), ROS/GSH, 대식세포 동태
4. **임상 엔드포인트** — FVC, DLCO, mPAP, KL-6, 호흡곤란 점수 시계열 추이
5. **시나리오 비교** — 6개 치료 전략 동시 비교
6. **가상환자 집단** — 몬테카를로 시뮬레이션 (n=20~200), VP 분포 및 반응자 분석

---

## Evans Syndrome (에반스 증후군) QSP Model

| 항목 | 내용 |
|------|------|
| **날짜** | 2026-06-20 |
| **분류** | 자가면역질환 / 혈액학 |
| **디렉토리** | [evans-syndrome/](evans-syndrome/) |
| **기계론적 지도** | [![ES QSP](evans-syndrome/es_qsp_model.png)](evans-syndrome/es_qsp_model.svg) |

### 질환 개요

에반스 증후군(Evans Syndrome, ES)은 **자가면역 용혈성 빈혈(AIHA)과 면역성 혈소판 감소증(ITP)이 동시에 또는 순차적으로 발생**하는 희귀 자가면역 혈액 질환입니다. Robert Evans가 1951년 최초로 기술하였으며, 다클론 IgG 자가항체가 적혈구 표면(Rh, Glycophorin A, Band 3)과 혈소판 표면(GPIbα, GPIIb/IIIa)을 동시에 표적으로 합니다.

### 핵심 병태생리 경로

| 경로 | 핵심 분자 | 임상 결과 |
|------|----------|---------|
| 자가항체 생성 (Anti-RBC) | Treg 결핍 → GC 과활성화 → IgG class-switch | DAT 양성, 온난자가항체, 혈관외 용혈 |
| 자가항체 생성 (Anti-Plt) | FcγRIIB 결핍 → B세포 억제 감소 → anti-GPIbα/GPIIbIIIa IgG | PAIgG 상승, 비장 혈소판 탐식 |
| 보체 활성화 | C1q → Classical Pathway → C3b opsonization → C3d (DAT C3d+) | 혈관내 용혈(MAC), C3d+ AIHA |
| Treg/Tfh 불균형 | mTOR↑ → Treg↓; IL-21↑ → GC Tfh 과활성화 | 면역관용 붕괴, 자가항체 지속 |
| 골수 보상 | EPO↑ → BFU-E 증가 → 망상적혈구↑; TPO↑ → 거핵구 생성↑ | 망상적혈구증가, 골수 고증식 |
| 비장 기능 | FcγRI/III 큰포식세포 → 적혈구·혈소판 탐식 | 비장비대, 세망내피계 파괴 증가 |

### 치료 시나리오 (6개)

| 시나리오 | 약물 | 주요 기전 | 예상 CR률 |
|---------|------|---------|---------|
| 무치료 (Natural Course) | — | — | — |
| 1차: Prednisone 단독 | 1–2 mg/kg/day | GR→B세포 사멸, FcγR↓ | ~40–60% (초기) |
| 2차: Pred + Rituximab | 0.5 mg/kg + 375mg/m² ×4 | B세포 CD20 표적 탈진 (6–12개월) | ~60–70% |
| 3차: Pred + MMF | 0.25mg/kg + 2g/day | IMPDH 억제 → 림프구 증식↓ | ~50–65% |
| 불응성: Sirolimus + Eltrombopag | 4mg/day + 50mg/day | mTOR↓→Treg↑; TPO-R 자극→Plt↑ | ~55–70% (Plt 우선) |
| 외과: 비장절제 | 수술 | 비장 FcγR 대식세포 제거 | ~60–80% (Plt) |

### 약물 PK 요약

| 약물 | 용량/경로 | Bioavailability | t½ | 표적 | 주요 PD 효과 |
|------|---------|---------|------|------|---------|
| Prednisolone | 1–2 mg/kg/day PO | 82% | 2.5 h | GR (NR3C1) | B세포 사멸 70%, FcγR 60%↓ |
| IVIg | 1–2 g/kg IV ×2일 | 100% | 3–4주 | FcRn / FcγR | FcγR blockade 75%, IgG 이화↑ |
| Rituximab | 375 mg/m² IV ×4주 | 100% | 22일 | CD20 | B세포 탈진 95% (6–12개월) |
| MMF/MPA | 1–3 g/day PO | 94% | 17 h | IMPDH | 림프구 항증식 65% |
| Eltrombopag | 50 mg/day PO | 52% | 26 h | c-Mpl (TPO-R) | 거핵구 2.5× 자극 → Plt↑ |
| Sirolimus | 2–6 mg/day PO | 14% | 62 h | mTORC1 | Treg 2×↑, T효과세포↓ |

### 모델 파일

| 파일 | 설명 |
|------|------|
| [es_qsp_model.dot](evans-syndrome/es_qsp_model.dot) | Graphviz 기계론적 지도 소스 (10개 클러스터, 138+ 노드) |
| [es_qsp_model.svg](evans-syndrome/es_qsp_model.svg) | SVG 벡터 이미지 |
| [es_qsp_model.png](evans-syndrome/es_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [es_mrgsolve_model.R](evans-syndrome/es_mrgsolve_model.R) | mrgsolve ODE 모델 (20 구획, 6 치료 시나리오) |
| [es_shiny_app.R](evans-syndrome/es_shiny_app.R) | Shiny 대시보드 (6탭) |
| [es_references.md](evans-syndrome/es_references.md) | 참고문헌 57편 (PubMed 링크 포함) |

### Shiny 앱 탭 구성

1. **환자 프로파일** — 인구통계, 기저 혈액 수치, 면역 지표, 이전 치료력
2. **약물 PK** — Prednisolone·Rituximab·MPA·Eltrombopag·Sirolimus·IVIg 혈장 농도 프로파일
3. **혈액학적 반응** — Hgb, 혈소판, 망상적혈구, 거핵구 풀 시계열 + CR 기준 시각화
4. **면역학적 지표** — Anti-RBC IgG, Anti-Plt IgG, Treg, 자가반응 B세포, C3b 보체 동태
5. **시나리오 비교** — 6개 치료 전략 동시 비교; Day 180/365 반응 요약 표
6. **바이오마커 대시보드** — DAT, PAIgG, LDH, Haptoglobin, Treg%, B cell% 실시간 평가 + 치료 권고


---

## 게실병 (Diverticular Disease — Diverticulosis / Diverticulitis)

### 개요

게실병은 서양 성인 인구의 약 35–50%에서 발생하는 만성 소화기 질환으로, 대장 점막이 근육층 약화 부위를 통해 탈출하여 게실(diverticulum)을 형성하는 질환입니다. 대부분은 무증상(게실증)으로 지내지만, 약 4–15%에서 급성 게실염(acute diverticulitis)이 발생하며, 복잡성 합병증(천공, 농양, 누공, 출혈)으로 진행할 수 있습니다.

### 주요 병태생리 경로

| 경로 | 핵심 메커니즘 | 임상 이상 |
|------|------------|---------|
| 대장 생역학 | 저섬유질 식이 → 소변량 감소 → 분절 수축↑ → 장내강 압력↑ | 게실 형성 (La Place 법칙) |
| ECM 리모델링 | 노화/TNF-α → MMP-1/9↑, TIMP-1↓ → 콜라겐 분해 → 근육층 약화 | 대장벽 취약화 |
| 장내 미생물 교란 | 저섬유질 → Faecalibacterium↓, E.coli↑ → LPS↑ | 만성 저강도 염증 |
| NF-κB 경로 | LPS → TLR4 → IKKβ → NF-κB → TNF-α, IL-6, IL-1β, COX-2 | 급성 게실염 |
| 급성 게실염 | 미세천공 → 주변 지방 염증 → 국소 농양 | Hinchey I–IV 분류 |
| 만성 염증/섬유화 | M1→M2 전환, TGF-β → 근섬유아세포 → 콜라겐 침착 | 협착, 누공 |

### 약물 PK/PD 파라미터

| 약물 | 작용기전 | 주요 PK | 임상 효과 |
|------|---------|---------|---------|
| Rifaximin 400mg TID × 7d/월 | RNA 중합효소 억제 (장관 내 작용, F<0.4%) | GI t½≈1h, 최소 전신 흡수 | 미생물 다양성↑, LPS↓, 재발↓ 40% |
| Mesalamine (5-ASA) 1.6g/day | NF-κB, COX-2 억제 | pH 방출형, t½≈0.5–2h | CRP↓, 재발↓ 32% |
| Ciprofloxacin 500mg BID × 10d | DNA 회전효소 억제 (그람음성) | F=70%, t½=4–6h | 급성 게실염 해소율 93% |
| Metronidazole 500mg TID × 10d | 혐기성균 DNA 손상 | F=99%, t½=6–8h | 혐기성균 제거 |
| Psyllium 15g/day | 부피형성 완하제, 단쇄지방산 전구체 | 발효 → 부티레이트 | 장내강 압력↓, 게실 형성↓ |

### Hinchey 분류 및 치료 원칙

| Hinchey | 소견 | 1차 치료 |
|---------|------|---------|
| Stage 0 | SUDD (증상성 합병증 없는 게실병) | 식이 조절 + Rifaximin ± Mesalamine |
| Stage I | 국소 농양 (pericolic) | 항생제 ± 경피적 배액 |
| Stage II | 골반 농양 | 항생제 + 경피적 배액 |
| Stage III | 화농성 복막염 | 응급수술 (Hartmann's procedure) |
| Stage IV | 대변성 복막염 | 응급수술 |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [div_qsp_model.dot](diverticular-disease/div_qsp_model.dot) | Graphviz 기계론적 지도 소스 (140+ 노드, 10 클러스터) |
| [div_qsp_model.svg](diverticular-disease/div_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [div_qsp_model.png](diverticular-disease/div_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [div_mrgsolve_model.R](diverticular-disease/div_mrgsolve_model.R) | mrgsolve ODE 모델 (20 구획, 6 시나리오) |
| [div_shiny_app.R](diverticular-disease/div_shiny_app.R) | Shiny 대시보드 (6 탭) |
| [div_references.md](diverticular-disease/div_references.md) | 참고문헌 (66개 PubMed 링크, 섹션별 분류) |
| [div_references.md](diverticular-disease/div_references.md) | 참고문헌 (66개 PubMed 링크, 섹션별 분류) |

---

## 자가면역 다발성 내분비병증 (APS Type 1 / APECED)

### 개요

자가면역 다발성 내분비병증(APS, Autoimmune Polyendocrinopathy Syndrome)은 두 개 이상의 내분비선이 자가면역 기전에 의해 손상되는 증후군입니다. APS Type 1(APECED)은 **AIRE 유전자 돌연변이**로 인한 중추 면역관용 파괴가 원인이며, 핀란드 1/25,000의 유병률을 보이는 희귀 단일 유전자 질환입니다. 3대 주요 구성 요소는 만성 피부칸디다증(CMC), 부갑상선기능저하증, 부신기능저하증(애디슨병)입니다.

### 주요 병태생리 경로

| 경로 | 핵심 기전 | 임상 결과 |
|------|---------|---------|
| AIRE 소실 → 흉선 | mTEC에서 TSA 발현 불능 → 자가반응 T세포 음성선택 실패 | 자가반응 T세포 말초 방출 |
| 자가항체 생성 | 항21-OH(부신), 항NALP5(부갑상선), 항GAD65(췌장), 항TPO(갑상선) | 다장기 자가면역 |
| 항IFN-α 중화항체 | APS1 병리표지 항체 (>95% 양성) | CMC 위험 증가 |
| 부신 파괴 | CTL + 자가항체 → 부신피질 위축 | 코르티솔/알도스테론 결핍 |
| 부갑상선 파괴 | 항NALP5 Ab → PTG 기능 소실 | 저칼슘혈증, 경련, QTc 연장 |
| 췌장 베타세포 파괴 | 항GAD65/IA-2 Ab + CD8+ CTL | T1DM (APS1 환자 20%) |
| 갑상선 파괴 | 항TPO Ab + Th1 매개 염증 | 하시모토 갑상선염, 기능저하증 |

### 약물 PK/PD 파라미터

| 약물 | 기전 | 주요 PK | 임상 효과 |
|------|------|---------|---------|
| 히드로코르티손 20mg/day | GR 작용제 (부신 대체) | F=0.96, t½=1.5h, CL=1.1L/min | AM 코르티솔 정상화 (8-20 µg/dL) |
| 플루드로코르티손 100µg/day | MR 작용제 (광물질코르티코이드) | F=0.99, t½=3.5h | Na+/K+ 균형, 기립성 저혈압 방지 |
| 칼시트리올 0.5µg/day | VDR 작용제 | t½=5-8h | 장관 Ca²⁺ 흡수 증가 → 저칼슘 교정 |
| 레보티록신 75-100µg/day | 갑상선호르몬 대체 | F=0.75-0.80, t½=9일 | TSH 정상화 (0.4-4.0 mIU/L) |
| 사이클로스포린 A 3.5mg/kg | 칼시뉴린 억제 → NFAT 차단 | F=0.30(가변), t½=8-24h, IC50=150ng/mL | AutoT 억제 ~50% |
| 아바타셉트 125mg SC/주 | CD28/B7 공자극 차단 | F=0.79, t½=13일, IC50=2µg/mL | Treg 기능 강화, AutoT 억제 |
| 리툭시맙 375mg/m² q6mo | 항CD20 → B세포 고갈 | IV, t½=21일, IC50=5µg/mL | 자가항체 ~85% 감소 |
| 토파시티닙 10mg/day | JAK1/3 억제 → STAT1/3/5 차단 | F=0.74, t½=3h, IC50=8ng/mL | IFN-γ/IL-17 신호 70% 차단 |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [aps_qsp_model.dot](autoimmune-polyendocrinopathy/aps_qsp_model.dot) | Graphviz 기계론적 지도 소스 (140+ 노드, 11 클러스터) |
| [aps_qsp_model.svg](autoimmune-polyendocrinopathy/aps_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [aps_qsp_model.png](autoimmune-polyendocrinopathy/aps_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [aps_mrgsolve_model.R](autoimmune-polyendocrinopathy/aps_mrgsolve_model.R) | mrgsolve ODE 모델 (20 구획, 7 시나리오) |
| [aps_shiny_app.R](autoimmune-polyendocrinopathy/aps_shiny_app.R) | Shiny 대시보드 (7탭: 환자/면역/PK/내분비/끝점/시나리오/바이오마커) |
| [aps_references.md](autoimmune-polyendocrinopathy/aps_references.md) | 참고문헌 60편 (PubMed 링크, 14 섹션 분류) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **자연경과 (중증 AIRE 돌연변이)**: 치료 없이 5년간 4개 장기 기능 상실 추적
2. **HRT 단독**: 히드로코르티손 20mg/day → 코르티솔 보충, 질병 진행은 지속
3. **HRT + 사이클로스포린 A**: T세포 억제로 자가반응 T세포 감소
4. **HRT + 아바타셉트**: CD28 공자극 차단 → Treg 기능 회복 + 장기 파괴 감소 (최선)
5. **HRT + 리툭시맙**: B세포 고갈로 자가항체 85% 감소 (6개월마다 투여)
6. **HRT + 토파시티닙**: JAK 억제로 IFN-γ/IL-17 신호 차단 → 염증 감소
7. **조기 개입 (경증 돌연변이)**: AIRE 30% 기능 소실 + 조기 HRT → 5년 후 장기 기능 유지



---

## 림프구성 뇌하수체염 (Lymphocytic Hypophysitis, LyH)

### 개요

림프구성 뇌하수체염(LyH)은 뇌하수체에 림프구가 침윤하여 뇌하수체 세포를 파괴하는 희귀 자가면역 질환으로, **다축 뇌하수체 기능 저하증**을 초래합니다. 주로 임신 주산기 여성에서 발생하며, 면역관문억제제(anti-CTLA4, anti-PD1) 투여 시에도 발생할 수 있습니다. 이 QSP 모델은 T세포·B세포·항뇌하수체 항체 매개 면역 반응과 다축 내분비 피드백(HPA·HPT·HPG·GH·PRL·ADH)을 포함한 복합 시스템을 구현합니다.

### 주요 병태생리 경로

| 경로 | 핵심 기전 | 임상 결과 |
|------|---------|---------|
| 자가항원 제시 | PIT-1, RABPHILIN-3L 항원 → MHC-II → CD4+ T세포 활성화 | 뇌하수체 특이적 자가면역 |
| Th1/Th17 극성화 | IL-12(→Th1), IL-23(→Th17) → IFN-γ, IL-17 분비 | CD8+ CTL 매개 뇌하수체세포 파괴 |
| B세포 활성화 | Tfh 도움 → 형질세포 → 항뇌하수체 항체(APA) | ADCC·보체 매개 추가 파괴 |
| HPA 축 손상 | 부신피질자극세포(Corticotrophs) 파괴 | ACTH 결핍 → 이차성 부신부전 (~70%) |
| HPT 축 손상 | 갑상선자극세포(Thyrotrophs) 파괴 | 중추성 갑상선기능저하증 (~40%) |
| HPG 축 손상 | 생식선자극세포(Gonadotrophs) 파괴 | 중추성 성선기능저하증 (~40%) |
| GH 축 손상 | 성장호르몬세포(Somatotrophs) 파괴 | 성인형 GH 결핍 (~30%) |
| 줄기 압박 효과 | 뇌하수체 줄기(인두돌기) 도파민 차단 | 고프로락틴혈증 (~50%) |
| 뇌하수체 후엽 침범 | AVP 신경세포 파괴 | 요붕증 (DI, ~25%) |
| 면역관문억제제 LyH | CTLA-4 차단 → 뇌하수체 CTLA-4 발현 표적 | ICI 치료 후 LyH (~10%) |

### 약물 PK/PD 파라미터

| 약물 | 기전 | 주요 PK | 임상 효과 |
|------|------|---------|---------|
| 프레드니솔론 60mg/day → taper | GR 활성화 → NF-κB 억제, 림프구 아폽토시스 | F=0.82, Vc=45L, CL=8.5L/h, t½=3h | 뇌하수체 용적 감소, 면역 억제 ~85% |
| 아자티오프린 150mg/day | 6-TGN → 퓨린 합성 억제 → T/B 증식 억제 | F=0.50, ka=0.80/h, CL=12L/h | T/B 증식 70% 억제 (steroid-sparing) |
| 리툭시맙 1000mg IV ×2 | 항CD20 → ADCC/CDC → B세포 고갈 | Vc=3.1L, CL=0.015L/h, t½=21일 | APA 생성 90% 억제 |
| 히드로코르티손 (대체) | 부신부전 호르몬 대체 | AM 20mg / PM 10mg | 코르티솔 정상화 (200-700nmol/L) |
| 데스모프레신 (DDAVP) | V2R 작용제 → AQP2 → 수분 재흡수 | IN 0.1mg, t½=75min | 요붕증 교정, 소변삼투압↑ |
| 레보티록신 (T4 대체) | T4 보충 → 갑상선 자극 호르몬 수용체 | F=0.75-0.80, t½=7일 | fT4 정상화 (12-22 pmol/L) |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [lhyp_qsp_model.dot](lymphocytic-hypophysitis/lhyp_qsp_model.dot) | Graphviz 기계론적 지도 소스 (154+ 노드, 14 클러스터) |
| [lhyp_qsp_model.svg](lymphocytic-hypophysitis/lhyp_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [lhyp_qsp_model.png](lymphocytic-hypophysitis/lhyp_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [lhyp_mrgsolve_model.R](lymphocytic-hypophysitis/lhyp_mrgsolve_model.R) | mrgsolve ODE 모델 (22 구획, 5 시나리오) |
| [lhyp_shiny_app.R](lymphocytic-hypophysitis/lhyp_shiny_app.R) | Shiny 대시보드 (7탭) |
| [lhyp_references.md](lymphocytic-hypophysitis/lhyp_references.md) | 참고문헌 60편 (PubMed 링크, 9 섹션 분류) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션, 22 ODE 구획)

1. **자연경과**: 치료 없이 2년간 다축 기능 저하 진행 추적
2. **고용량 프레드니솔론 taper**: 60→40→30→20→10→5 mg/day 점감 — 최적 치료
3. **프레드니솔론 + 아자티오프린**: 스테로이드 절약 병합 요법 — 장기 부작용↓
4. **리툭시맙**: ICI 연관 LyH에서 항CD20 B세포 고갈 — APA 90% 감소
5. **저용량 프레드니솔론 유지**: 10mg/day 유지 요법 — 준최적 치료

### Shiny 대시보드 7탭 구성

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 질환 중증도, 치료 선택, LyH 아형 설정 |
| 2. 약물 PK | 프레드니솔론·아자티오프린·리툭시맙 혈중 농도 |
| 3. 뇌하수체 기능 | 기능점수, 염증 용적, 호르몬 결핍 히트맵 |
| 4. 호르몬 축 | HPA(ACTH/코르티솔), HPT(TSH/fT4), GH/IGF-1, HPG(FSH/LH/E2) |
| 5. 면역 역학 | T세포·B세포·APA 동태, 면역-뇌하수체 위상도 |
| 6. 임상 끝점 | 방사형 차트, 부신 위기 위험, 핵심 시점 요약표 |
| 7. 바이오마커 | 선택 바이오마커 추적, 정상 참조 범위 표시 |

---

## 만성 신우신염 (Chronic Pyelonephritis, CPN)

### 개요

만성 신우신염은 반복적인 상부 요로 세균 감염 또는 방광요관역류(VUR)로 인해 신장 실질에 영구적인 반흔(scar)이 형성되는 만성 신장 질환입니다. 미국에서 매년 약 100만 건의 입원 관련 신우신염이 발생하며, 소아에서 CKD 원인의 10–20%를 차지합니다.

### 주요 병태생리 경로

| 경로 | 핵심 메커니즘 | 임상 이상 |
|------|------------|---------|
| UPEC 병원성 | FimH/PapG 부착 → HlyA 독소 → 조직 침입 | 신우 감염 → 실질 손상 |
| VUR 매개 역류 | VUR Grade I–V → 신내 역류 → 복합 유두 취약 | 반복 감염 → 신장 반흔 |
| TLR4/NFκB | LPS → TLR4 → MyD88 → NFκB | IL-6/IL-8/TNF-α 과다 생성 |
| NLRP3 인플라마좀 | HlyA → NLRP3 → Caspase-1 → IL-1β | 세뇨관 세포 사멸 (pyroptosis) |
| TGF-β1 섬유화 | TGF-β1 → EMT → 근섬유모세포 → Collagen I/III | DMSA 흡수 결손, GFR 저하 |
| RAAS 활성화 | Scar → RAAS → Ang II → TGF-β1 | 이차성 고혈압, 진행 CKD |

### 항생제 PK/PD 파라미터

| 약물 | 작용기전 | 주요 PK | PK/PD 목표 |
|------|---------|---------|----------|
| Ciprofloxacin 500 mg BID | DNA Gyrase/TopoisIV 억제 (AUC-dep.) | F=70%, t½=4h, CL=25 L/h | fAUC/MIC > 125 |
| TMP-SMX 160/800 mg BID | DHFR 억제 (시간-의존형) | F=95%/85%, t½=10h | T>MIC > 40% |
| Nitrofurantoin 100 mg QD | DNA/단백 손상 (요중 농도 의존) | F=75%, t½<1h | Cmax/MIC(요중) > 4 |
| Fosfomycin 3g 단회 | 세포벽 합성 억제 (MurA) | 요중 30–40% 회수 | Cmax > 4× MIC (요중) |
| Gentamicin 5 mg/kg IV QD | 30S 리보솜 억제 (농도-의존형) | Vd=0.3 L/kg, t½=2h | Cmax/MIC > 10 |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [cpn_qsp_model.dot](chronic-pyelonephritis/cpn_qsp_model.dot) | Graphviz 기계론적 지도 소스 (10 클러스터, 140+ 노드) |
| [cpn_qsp_model.svg](chronic-pyelonephritis/cpn_qsp_model.svg) | SVG 벡터 이미지 |
| [cpn_qsp_model.png](chronic-pyelonephritis/cpn_qsp_model.png) | PNG 이미지 (150 dpi) |
| [cpn_mrgsolve_model.R](chronic-pyelonephritis/cpn_mrgsolve_model.R) | mrgsolve ODE 모델 (18 구획, 7 시나리오) |
| [cpn_shiny_app.R](chronic-pyelonephritis/cpn_shiny_app.R) | Shiny 대시보드 (7탭) |
| [cpn_references.md](chronic-pyelonephritis/cpn_references.md) | 참고문헌 38편 (PubMed 링크 포함) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션, 18 ODE 구획)

1. **미치료**: 세균 부하 지속 → 빠른 GFR 저하 및 신장 반흔 형성
2. **Ciprofloxacin 14일**: fAUC/MIC > 125 달성 → 7–10일 내 세균 청소, GFR 보존
3. **TMP-SMX 14일**: 감수성 균주에서 동등한 효과 (내성률 증가 시 주의)
4. **Ciprofloxacin 7일**: 14일 대비 약간 열등하나 내성 선택압 감소
5. **Nitrofurantoin 6개월 예방**: 재발률 60% 감소, 장기 GFR 보존 효과
6. **Cipro 14일 + Nitrofurantoin 6개월**: 최우수 전략 — 급성 치료 후 재발 예방
7. **TMP-SMX + 내성균 (MIC × 4)**: 치료 실패 시뮬레이션 → 빠른 섬유화 진행

### Shiny 대시보드 7탭 구성

| 탭 | 내용 |
|----|------|
| ① 환자 프로파일 | GFR, VUR 등급, 동반질환, 항생제 선택 |
| ② PK 모니터링 | Cipro/TMP/NIT 혈중 농도, fAUC/MIC |
| ③ 세균 동태 | log₁₀ CFU, 생물막, Kill Rate, 면역세포 |
| ④ 신기능 (GFR) | GFR 추이, 크레아티닌, 반흔, CKD 단계 전환 |
| ⑤ 시나리오 비교 | 5가지 항생제 전략 동시 비교 |
| ⑥ 바이오마커 | IL-6, TGF-β1, 콜라겐, 요로패혈증 위험 |
| ⑦ 참고문헌 | 섹션별 분류 PubMed 링크 |

---

## 무증상 고요산혈증 (Asymptomatic Hyperuricemia, AHY)

### 개요

무증상 고요산혈증(AHY)은 혈청 요산(SUA) ≥7.0mg/dL (남성) / ≥6.0mg/dL (여성)이면서 통풍 발작·토파이·요산 신석증이 없는 상태입니다. 전 세계 성인의 15-25%에서 발생하며, 아시아인에서 더 높은 빈도를 보입니다. 장기적으로 통풍, 만성 신부전(CKD), 고혈압, 심혈관 질환의 중요한 위험인자입니다.

### 주요 병태생리 경로

| 경로 | 핵심 메커니즘 | 임상 이상 |
|------|-------------|---------|
| 퓨린 이화 | XO/XDH: Hypoxanthine → Xanthine → UA (+ H₂O₂/ROS) | 과생성형 고요산혈증 (~10-15%) |
| 신장 재흡수 | URAT1(SLC22A12), GLUT9(SLC2A9) 과다 재흡수 | 배설 감소형 (~85%) |
| 장관 분비 | ABCG2 Q141K 변이: 장관+신장 분비 50%↓ | 유전성 고요산혈증 |
| MSU 결정 형성 | SUA > 6.8mg/dL → 핵화·성장·NLRP3 활성 | 통풍 전단계 |
| 심혈관 영향 | UA → eNOS↓ → NO↓ → 혈압↑ → 죽상경화↑ | 고혈압, 심혈관 이벤트 +27% |
| 신장 영향 | UA 직접독성 + ROS → 세뇨관 손상 → GFR↓ | CKD 가속 진행 |
| 대사 영향 | UA → 인슐린 저항성↑ → HOMA-IR↑ | 대사증후군, NAFLD |

### 약물 PK/PD 파라미터

| 약물 | 기전 | 주요 PK | 임상 효과 |
|------|------|---------|---------|
| 알로퓨리놀 300mg/day | XO 경쟁적/비가역 억제 (옥시퓨리놀) | F=90%, t½(oxypurinol)=18-30h | SUA↓ 30-40% |
| 페북소스타트 80mg/day | XO 비퓨린 선택적 억제 (IC50=0.6nM) | F=84%, t½=5-8h | SUA↓ 40-53% |
| 레시누라드 200mg/day | URAT1+OAT4 이중 억제 | F≈100%, t½=5h | SUA↓ +15-30% (병용) |
| 도티누라드 4mg/day | URAT1 선택적 억제 | F≈100%, t½=14-17h | SUA↓ 15-20% |
| 페글로티카제 8mg IV q2주 | UA→알란토인 (uricase, 용해도 10배↑) | t½=6-14일 (PEG화) | SUA↓ >80% (24h 내) |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [ahy_qsp_model.dot](asymptomatic-hyperuricemia/ahy_qsp_model.dot) | Graphviz 기계론적 지도 소스 (15 클러스터, 100+ 노드) |
| [ahy_qsp_model.svg](asymptomatic-hyperuricemia/ahy_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [ahy_qsp_model.png](asymptomatic-hyperuricemia/ahy_qsp_model.png) | PNG 이미지 (150 dpi) |
| [ahy_mrgsolve_model.R](asymptomatic-hyperuricemia/ahy_mrgsolve_model.R) | mrgsolve ODE 모델 (19 구획, 7 시나리오, 용량-반응 분석) |
| [ahy_shiny_app.R](asymptomatic-hyperuricemia/ahy_shiny_app.R) | Shiny 대시보드 (6탭 인터랙티브 시뮬레이터) |
| [ahy_references.md](asymptomatic-hyperuricemia/ahy_references.md) | 참고문헌 45편 (PubMed 링크, 11개 섹션) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션, 19 ODE 구획)

1. **미치료 AHY (SUA=7.5)**: 2년간 SUA 점진적 상승, GFR 감소, MSU 결정 누적, 심혈관 위험 악화
2. **알로퓨리놀 300mg/day**: XO 억제 → SUA 5.2mg/dL (치료 목표 달성), GFR 안정
3. **페북소스타트 80mg/day**: 고선택적 XO 억제 → SUA 4.5mg/dL (더 강력한 효과)
4. **알로퓨리놀 + 요산배설촉진제 병용**: 이중 기전 → SUA 4.0mg/dL (적극적 treat-to-target)
5. **고과당식이 + 알코올**: SUA 9.5mg/dL 급격 상승, 통풍 발작 위험 10배 증가
6. **ABCG2 Q141K 변이 + 페북소스타트 120mg**: 유전적 고위험군 최대용량 치료
7. **적극적 치료목표 SUA<6mg/dL**: 알로퓨리놀 600mg + 생활습관 → 심혈관/신장 보호 극대화

### Shiny 대시보드 6탭 구성

| 탭 | 내용 |
|----|------|
| ① 환자 프로파일 | 질환 개요, 위험 평가, SUA 시계열, 초기 value boxes |
| ② 약물 PK | 혈중 농도-시간 곡선 (알로·페북소·요산배설촉진제), XO 억제율 |
| ③ 요산 역학 | SUA, 요중 UA, MSU 결정 침착량, 통풍 발작 위험 지수 |
| ④ 심혈관·신장 영향 | eGFR, MAP, 누적 CV 위험, 내피기능/NO 수준 |
| ⑤ 시나리오 비교 | 7개 치료 시나리오 SUA 동시 비교, 2년 결과표 (DT) |
| ⑥ 바이오마커 | IL-1β, hs-CRP, HOMA-IR, SUA-CRP 산점도, 최종 바이오마커 요약 |

---

## 자궁내막증 (Endometriosis)

### 개요

자궁내막증은 전 세계 가임기 여성의 약 10%(약 1억 9천만 명)에서 발생하는 만성 에스트로겐 의존성 질환입니다. 자궁내막 조직이 자궁 밖(복막, 난소, 더글라스와)에 이소성으로 착상하여 월경 주기에 반응하고, 만성 골반통 · 월경통 · 성교통 · 불임을 유발합니다. 평균 진단까지 7–10년이 소요되는 진단 지연이 큰 문제입니다.

Endometriosis affects ~10% of reproductive-age women (190 million globally). Ectopic endometrial tissue implants outside the uterus, driven by estrogen, inflammation, and immune dysregulation, causing chronic pelvic pain and infertility.

### 주요 병태생리 경로

| 경로 | 핵심 메커니즘 | 임상 이상 |
|------|------------|---------|
| HPO 축 이상 | GnRH 펄스 → FSH/LH → E2 과잉 자극 이소성 내막 | 고에스트로겐 환경 유지 |
| 역행성 월경 | 자궁내막 세포 복강 내 착상 → 이소성 병소 형성 | rAFS Stage I–IV |
| 아로마타제 과발현 | PGE2 → SF-1 → CYP19A1↑ → 국소 E2 자가증폭 | 병소 내 E2 300–3000 pmol/L |
| 프로게스테론 저항성 | PR-B↓, PR-A 우위 → HOXA10 발현 억제 | 내막 분화 불능, 불임 |
| PGE2 양성 피드백 | COX-2↑ → PGE2↑ → cAMP → 아로마타제↑ → E2↑ | 통증·병소 동시 악화 |
| NK세포 기능 저하 | E2·복강액 → NK 세포독성↓ → 이소성 세포 생존 | 병소 제거 실패 |
| NGF/통증 경로 | 병소 → NGF↑ → TrkA → TRPV1 → 중추감작 | 중증 월경통 · 만성골반통 |
| 산화 스트레스 | 역행 혈액 → 철분 과부하 → Fenton → ROS↑ | DNA 손상, 병소 성장 촉진 |

### 약물 PK/PD 파라미터

| 약물 | 기전 | 주요 PK | 임상 효과 |
|------|------|---------|---------|
| Leuprolide Depot 3.75mg/4wk | GnRH 수용체 탈감작(desensitization) | ka_dep=0.005/h, t½_depot≈14d | E2 <20 pg/mL, 병소 ↓40–60% |
| Elagolix 150mg/d | GnRH 경쟁적 길항 (저용량 부분 억제) | ka=1.2/h, t½=4–6h, F=57% | E2 ~70 pg/mL, 월경통 ↓45% |
| Elagolix 200mg BID | GnRH 경쟁적 길항 (고용량 완전 억제) | AUC 2× 저용량, ka=1.2/h | E2 <12 pg/mL, 월경통 ↓76% |
| Dienogest 2mg/d | PR 작용 (항증식, 항염증, 항혈관신생) | ka=1.5/h, t½=11h, F=91% | 병소 부피 ↓50–70%, IL-8↓ |
| Letrozole 2.5mg/d + 노르에티스테론 add-back | 아로마타제 비가역 억제 + 에스트로겐 보충 | ka_let=1.8/h, t½=34h | 국소 E2↓90%, BMD 보호 |
| Combined OCP (EE+progestin) | 혈중 E2 억제 + 탈락막화 | EE 30mcg/d 주기적/연속적 | 월경통↓ 40–50% |
| NSAID (Ibuprofen 400–800mg TID) | COX-1/2 억제 → PGE2↓ | ka=1.5/h, t½=2h | 통증 NRS↓ 2–3점 (일시적) |

### 치료 시나리오 (mrgsolve 시뮬레이션)

| 시나리오 | 약물 | 목표 | 12개월 예측 효과 |
|---------|------|------|----------------|
| 0 | 미치료 | 자연 경과 | 병소 2→6cm³, 통증 NRS 4→7 |
| 1 | Leuprolide Depot 3.75mg q28d | E2<20 pg/mL, 무배란 | 병소↓55%, 통증↓65%, BMD↓2.5% |
| 2 | Elagolix 150mg/d | 부분 억제, 월경 유지 | 병소↓30%, 통증↓45%, BMD↓1% |
| 3 | Elagolix 200mg BID | 완전 억제, 무월경 | 병소↓50%, 통증↓75%, BMD↓2% |
| 4 | Dienogest 2mg/d | 항증식·항염증 | 병소↓40%, 통증↓55%, BMD 변화 미미 |
| 5 | Letrozole 2.5mg/d + 노르에티스테론 5mg/d | 국소 E2↓, 전신 E2 보호 | 병소↓60%, 통증↓70%, BMD 보호 |
| 6 | Leuprolide Depot + 노르에티스테론 Add-back | E2 억제 + BMD 보호 | 병소↓55%, 통증↓65%, BMD↓0.5% |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [endo_qsp_model.dot](endometriosis/endo_qsp_model.dot) | Graphviz DOT 소스 (176노드, 9클러스터, 267엣지) |
| [endo_qsp_model.svg](endometriosis/endo_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [endo_qsp_model.png](endometriosis/endo_qsp_model.png) | PNG 이미지 (150 dpi) |
| [endo_mrgsolve_model.R](endometriosis/endo_mrgsolve_model.R) | mrgsolve ODE 모델 (22구획, 7시나리오) |
| [endo_shiny_app.R](endometriosis/endo_shiny_app.R) | Shiny 대시보드 (6탭: 환자프로파일·PK·HPO호르몬·통증임상·병소역학·시나리오비교) |
| [endo_references.md](endometriosis/endo_references.md) | 참고문헌 44편 (PubMed 링크, 10개 섹션) |
| [README.md](endometriosis/README.md) | 디렉토리 상세 문서 |

### 핵심 참고문헌

1. Giudice LC & Kao LC. "Endometriosis." *Lancet* 2004;364:1789–99. [PMID 15488215](https://pubmed.ncbi.nlm.nih.gov/15488215/)
2. Zondervan KT et al. "Endometriosis." *N Engl J Med* 2020;382:1244–56. [PMID 31777792](https://pubmed.ncbi.nlm.nih.gov/31777792/)
3. Taylor HS et al. "Treatment of Endometriosis-Associated Pain with Elagolix." *N Engl J Med* 2017;377:28–40. [PMID 28985706](https://pubmed.ncbi.nlm.nih.gov/28985706/)
4. Bulun SE et al. "Aromatase in endometriosis." *Semin Reprod Med* 2004;22:45–50. [PMID 15083381](https://pubmed.ncbi.nlm.nih.gov/15083381/)
5. Strowitzki T et al. "Dienogest is as effective as leuprolide acetate in treating the painful symptoms of endometriosis." *Hum Reprod* 2010;25:633–41. [PMID 20027944](https://pubmed.ncbi.nlm.nih.gov/20027944/)

---

## 알츠하이머병 (Alzheimer's Disease, AD)

### 개요 / Overview

알츠하이머병은 전 세계 치매 원인의 약 60–70%를 차지하는 가장 흔한 신경퇴행성 질환입니다.  
전 세계 약 5,500만 명의 치매 환자 중 다수가 AD 환자이며, 2050년까지 3배 증가 예상.  
아밀로이드 폭포 가설(Amyloid Cascade Hypothesis), 타우 병리(NFT), 신경염증, 콜린성 결핍이 주요 병인입니다.

Alzheimer's disease accounts for 60–70% of all dementia cases globally. Core pathological mechanisms: amyloid cascade (APP→BACE1→Aβ42 overproduction/aggregation), tau hyperphosphorylation/NFT formation, microglial neuroinflammation (M1/DAM activation, NLRP3 inflammasome), and loss of basal forebrain cholinergic neurons (BFCN).

### 병태생리 핵심 경로

| 경로 | 핵심 분자 | 임상 결과 |
|------|----------|---------|
| APP 가공/아밀로이드 생성 | BACE1, γ-세크레타제(PSEN1/2), Aβ42 | 노인반(Senile Plaque) 형성 |
| 아밀로이드 응집/독성 | Aβ 올리고머 → 프로토피브릴 → 섬유 → 플라크 | 시냅스 손상, LTP 억제, AMPAR 소실 |
| 아밀로이드 소거 | IDE, NEP, LRP1(efflux), RAGE(influx), ApoE4(↓clearance) | CSF Aβ42 감소, Amyloid PET SUVr↑ |
| 타우 과인산화 | CDK5/p25, GSK-3β, PP2A↓, pTau181/217/231 | NFT 형성, Braak stage 진행 |
| 타우 확산(prion-like) | tau oligomer 세포간 전달, BIN1 | 피질 전반 확산 |
| 신경염증 | 소교세포 M1/DAM, TREM2, NLRP3, TNF-α, IL-1β | 시냅스 제거 (C1q/C3), NeuroCellDeath |
| 콜린성 기능 저하 | ChAT↓, ACh↓, BFCN 손실, nAChR/mAChR↓ | 인지·기억 기능 저하 |
| 시냅스 소실 | Dendritic spine↓, LTP↓, CaMKII, BDNF/TrkB↓ | MMSE↓, ADAS-Cog↑ |
| 미토콘드리아 기능 이상 | Complex I 손상, ROS, mPTP, Caspase-3 | ATP 결핍, 신경세포 사멸 |
| 혈뇌장벽 이상 | BBB 투과성↑, LRP1↓, CAA, 글림프계 이상 | ISF Aβ 축적, ARIA 위험 |

### 약물 PK/PD 파라미터

| 약물 | 작용기전 | 주요 PK | 임상 효과 |
|------|---------|---------|---------|
| Donepezil 10 mg QD | AChE 가역 억제 (Ki≈6.7 nM) | F≈100%, t½≈70h, Vd=594L, Kp_brain=0.18 | ADAS-Cog −3.1점, MMSE +2–3점 (24주) |
| Rivastigmine 9.5 mg/24h 패치 | AChE + BuChE 유사비가역 억제 | 패치 F=36%, t½≈1h | ADAS-Cog −2.6점 |
| Galantamine 24 mg | AChE 경쟁적 억제 + α7-nAChR 알로스테릭 조절 | F≈100%, t½≈7h | ADAS-Cog −3.0점 |
| Memantine 20 mg QD | NMDA 수용체 개방채널 차단 (IC50≈0.8 mg/L) | t½≈60–80h, Vd=520L | SIB +5.7점, 중등-중증 AD |
| Lecanemab 10 mg/kg Q2W IV | Aβ 프로토피브릴 고친화 결합 (Kd≈1 nM) | t½≈5일, CNS/plasma≈0.1% | CDR-SB −27%, Amyloid PET −59.1 CL (CLARITY-AD, 2023) |
| Donanemab Q4W IV | 피로글루타메이트-Aβ(pE3-Aβ) 플라크 표적 | t½≈7일 (mAb) | CDR-SB −35% (TRAILBLAZER-ALZ2, 2023) |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [`ad_qsp_model.dot`](alzheimers-disease/ad_qsp_model.dot) | 기계론적 지도: 12 클러스터, 140+ 노드, fdp 레이아웃 |
| [`ad_qsp_model.svg`](alzheimers-disease/ad_qsp_model.svg) | 기계론적 지도 (SVG, 벡터) |
| [`ad_qsp_model.png`](alzheimers-disease/ad_qsp_model.png) | 기계론적 지도 (PNG, 150 dpi) |
| [`ad_mrgsolve_model.R`](alzheimers-disease/ad_mrgsolve_model.R) | mrgsolve ODE 모델: 21 구획, 7 시나리오 (대조군/donepezil/memantine/병용/lecanemab/APOE4/FAD) |
| [`ad_shiny_app.R`](alzheimers-disease/ad_shiny_app.R) | Shiny 대시보드 6탭: 환자 프로파일 · 약물 PK · 아밀로이드/타우 바이오마커 · 신경염증/시냅스 · 인지 엔드포인트 · 시나리오 비교 |
| [`ad_references.md`](alzheimers-disease/ad_references.md) | 참고문헌 41편 (PubMed PMID): 아밀로이드 폭포 · 타우 병리 · 신경염증 · 콜린성계 · 약물 PK/PD · 임상시험 · 유전위험인자 · QSP 모델링 · 바이오마커 |
| [`README.md`](alzheimers-disease/README.md) | 질환 및 모델 상세 설명 |

### 핵심 참고문헌

1. Hardy J & Selkoe DJ. "The amyloid hypothesis of Alzheimer's disease." *Science* 2002;297:353–6. [PMID 12130773](https://pubmed.ncbi.nlm.nih.gov/12130773/)
2. Heneka MT et al. "Neuroinflammation in Alzheimer's disease." *Lancet Neurol* 2015;14:388–405. [PMID 25792098](https://pubmed.ncbi.nlm.nih.gov/25792098/)
3. van Dyck CH et al. "Lecanemab in Early Alzheimer's Disease." *N Engl J Med* 2023;388:9–21. [PMID 36449413](https://pubmed.ncbi.nlm.nih.gov/36449413/)
4. Sims JR et al. "Donanemab in Early Symptomatic Alzheimer's Disease." *N Engl J Med* 2023;389:1919–31. [PMID 37459141](https://pubmed.ncbi.nlm.nih.gov/37459141/)
5. Iqbal K et al. "Tau and neurodegenerative disease: the story so far." *Nat Rev Neurol* 2016;12:15–27. [PMID 26635538](https://pubmed.ncbi.nlm.nih.gov/26635538/)

---

## 아토피 피부염 (Atopic Dermatitis, AD)

### 개요 / Overview

아토피 피부염(AD)은 전 세계 소아의 15–20%, 성인의 5–10%에서 발생하는 가장 흔한 만성 염증성 피부질환입니다.
주요 병태생리: FLG 유전자 변이에 의한 표피 장벽 손상 → TSLP/IL-33/IL-25 알라민 분비 → ILC2 및 Th2 세포 활성화 → IL-4/IL-13/IL-31 과다 분비 → STAT6 신호전달 → FLG 추가 억제(악순환), IgE 상승, 가려움증(IL-31→DRG 뉴런→TRPV1/TRPA1).

Atopic dermatitis is driven by a Th2/ILC2-polarized immune response with IL-4 and IL-13 as the core pathogenic cytokines acting through the shared IL-4Rα chain. JAK1/TYK2 downstream STAT6 phosphorylation suppresses filaggrin expression creating an itch-barrier vicious cycle. Key clinical targets: dupilumab (IL-4Rα), tralokinumab/lebrikizumab (IL-13), nemolizumab (IL-31Ra), upadacitinib/abrocitinib (JAK1), baricitinib (JAK1/2).

### 병태생리 핵심 경로

| 경로 | 핵심 분자 | 임상 결과 |
|------|----------|---------|
| 표피 장벽 손상 | FLG↓, 세라마이드↓, CLDN1↓, TEWL↑, pH↑, KLK5↑ | 알레르겐 침투 촉진, S. aureus 집락화 |
| 알라민 분비 (KC) | TSLP, IL-33/ST2, IL-25 | ILC2, DC, 비만세포 활성화 |
| Th2 분화 | GATA-3, OX40L, IL-4, IL-13 | STAT6 인산화, IgE 생산 |
| JAK-STAT 신호 | JAK1/TYK2→pSTAT6 | FLG 발현↓, TARC↑, IgE↑ |
| 가려움 신호 | IL-31→IL-31Ra(DRG)→TRPV1/TRPA1 | NRS 가려움↑, 긁기-장벽 악순환 |
| 호산구 활성화 | IL-5, Eotaxin-3(CCL26), ECP, MBP | 조직 손상, 혈중 호산구↑ |
| IgE/비만세포 | IgE→FcεRI→비만세포 탈과립, PGD2, LTC4 | 즉시형 알레르기 반응 |
| 만성화 | IL-22(Th22)→표피 증식, TGF-β→태선화 | EASI 악화, 삶의 질↓ |

### 약물 PK/PD 파라미터

| 약물 | 작용기전 | 주요 PK | 임상 효과 (wk16) |
|------|---------|---------|-----------------|
| Dupilumab 300mg Q2W SC | IL-4Rα 차단 (IL-4+IL-13 동시 억제) | t½≈21일, F_SC=64%, CL=0.21 L/d | EASI-75: 44%, IGA0/1: 38% (SOLO-1/2) |
| Tralokinumab 300mg Q2W SC | 항-IL-13 단일 차단 | t½≈22일, F_SC=72% | EASI-75: 38%, IGA0/1: 25% (ECZTRA-1/2) |
| Lebrikizumab 250mg Q2W SC | 항-IL-13 (느린 해리, KD<0.1 nM) | t½≈25일 | EASI-75: 43%, IGA0/1: 31% (ADvocate 1/2) |
| Nemolizumab 60mg Q4W SC | 항-IL-31Ra (가려움 선택적) | t½≈22일, RO>80% | NRS -49%, EASI-75: 26% (AD-301) |
| Upadacitinib 30mg QD oral | JAK1 선택적 억제 (IC50=0.045 ng/mL) | t½≈8시간, Vd=166 L | EASI-75: 62%, NRS -55% (Rising Up) |
| Abrocitinib 200mg QD oral | JAK1 선택적 억제 (IC50~0.07 ng/mL) | t½≈3–5시간 | EASI-75: 63%, IGA0/1: 44% (JADE MONO-1/2) |
| Baricitinib 4mg QD oral | JAK1/JAK2 억제 | t½≈12시간 | EASI-75: 48%, IGA0/1: 31% (BREEZE-AD) |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [`ad_qsp_model.dot`](atopic-dermatitis/ad_qsp_model.dot) | 기계론적 지도: 13 클러스터, 130+ 노드, ortho 레이아웃 |
| [`ad_qsp_model.svg`](atopic-dermatitis/ad_qsp_model.svg) | 기계론적 지도 (SVG, 벡터) |
| [`ad_qsp_model.png`](atopic-dermatitis/ad_qsp_model.png) | 기계론적 지도 (PNG, 150 dpi) |
| [`ad_mrgsolve_model.R`](atopic-dermatitis/ad_mrgsolve_model.R) | mrgsolve ODE 모델: 20 구획 (dupilumab 2-cmt, upadacitinib 1-cmt, nemolizumab), 6 시나리오 |
| [`ad_shiny_app.R`](atopic-dermatitis/ad_shiny_app.R) | Shiny 대시보드 6탭: 환자프로파일 · 약물PK/RO · PD바이오마커 · 임상엔드포인트 · 시나리오비교 · 장벽&가려움 |
| [`ad_references.md`](atopic-dermatitis/ad_references.md) | 참고문헌 52편 (PubMed 링크, 15개 섹션) |
| [`README.md`](atopic-dermatitis/README.md) | 질환 및 모델 상세 설명 |

### 핵심 참고문헌

1. Langan SM, Irvine AD, Weidinger S. "Atopic dermatitis." *Lancet* 2020;396:345–60. [PMID 32738956](https://pubmed.ncbi.nlm.nih.gov/32738956/)
2. Simpson EL et al. (SOLO 1&2). "Two phase 3 trials of dupilumab versus placebo in atopic dermatitis." *N Engl J Med* 2016;375:2335–48. [PMID 27690741](https://pubmed.ncbi.nlm.nih.gov/27690741/)
3. Guttman-Yassky E et al. (Rising Up). "Upadacitinib in adults with moderate-to-severe atopic dermatitis." *Lancet* 2020;395:1754–68. [PMID 32416769](https://pubmed.ncbi.nlm.nih.gov/32416769/)
4. Palmer CN et al. "Common loss-of-function variants of filaggrin are a major predisposing factor for atopic dermatitis." *Nat Genet* 2006;38:441–6. [PMID 16550169](https://pubmed.ncbi.nlm.nih.gov/16550169/)
5. Cevikbas F et al. "A sensory neuron-expressed IL-31 receptor mediates T helper cell-dependent itch." *J Allergy Clin Immunol* 2014;133:448–60. [PMID 24373349](https://pubmed.ncbi.nlm.nih.gov/24373349/)


---

## 파킨슨병 (Parkinson's Disease) — 2026-06-20

파킨슨병(PD)은 흑질(substantia nigra pars compacta, SNpc)의 도파민 신경세포가 진행성으로 소실되고, 남은 신경세포에 루이소체(Lewy body, α-시누클레인 응집체)가 형성되는 신경퇴행성 질환입니다.

Parkinson's disease is a progressive neurodegenerative disorder characterized by dopaminergic neuron loss in the SNpc (~60% at motor symptom onset) and Lewy body pathology driven by α-synuclein aggregation. Key drug targets: dopamine replacement (levodopa), dopamine agonists (D2/D3: pramipexole, ropinirole), MAO-B inhibitors (rasagiline, safinamide), COMT inhibitors (entacapone, opicapone), and DBS (STN/GPi).

### 병태생리 핵심 경로

| 경로 | 핵심 분자 | 임상 결과 |
|------|----------|---------|
| α-시누클레인 응집 | SNCA, LRRK2, pS129-aSyn, 루이소체 | 뉴런 독성, Braak 병기 진행 |
| 미토콘드리아 기능장애 | Complex I 억제, ROS, PINK1/Parkin, mPTP | 산화 스트레스, 세포자멸사 |
| 단백질 품질 관리 붕괴 | UPS 장애, CMA 차단(LAMP-2A), 자가포식 | aSyn 축적 가속화 |
| 신경염증 | 미세아교세포 M1 활성화, NLRP3, TNF-α, IL-1β | 도파민 뉴런 추가 손상 |
| 도파민 합성 감소 | TH 활성↓, AADC↓, DAT↓, VMAT2↓ | 선조체 DA 결핍 |
| 기저핵 회로 불균형 | GPi 과활성, STN 과흥분, 직접/간접 경로 불균형 | 떨림, 경직, 서동증 |
| 루이소체 전파 | 분비된 aSyn, 프리온 유사 전파, 미주신경 | Braak 단계 1→6 |
| 장-뇌 축 | 장내 미생물 불균형, ENS, 미주신경 역행 전파 | 비운동 증상, 자율신경 장애 |

### 약물 PK/PD 파라미터

| 약물 | 작용기전 | 주요 PK | 임상 효과 |
|------|---------|---------|---------|
| 레보도파/카르비도파 250mg TID | AADC → DA 전환, D1/D2 자극 | Ka=1.2 h⁻¹, t½≈1.5h, BBB: LAT1 | UPDRS-III −20pts; wearing-off/LID 위험 |
| 프라미펙솔 0.75mg TID | D3>D2 수용체 완전효현제 | Ka=0.8 h⁻¹, Vd=490L, t½=8–12h | UPDRS-III −11pts; LID 발생 지연 |
| 로피니롤 2mg TID | D3>D2 수용체 효현제 | t½=6h, CYP1A2 대사 | 초기 단독요법 또는 병용 |
| 라사질린 1mg QD | MAO-B 비가역적 억제 (~98%) | t½=1.3h, 효소회복 t½≈7d | UPDRS-III −3.6pts (ADAGIO) |
| 엔타카폰 200mg (각 LD와 함께) | COMT 가역적 억제 (~65%) | t½=0.4h, 단기 Imax 작용 | LD t½ 연장, wearing-off 감소 |
| 오피카폰 50mg QD | COMT 강력/지속 억제 (~70%) | t½=0.7h, 24h COMT 억제 유지 | 1일 1회 복용, LD 안정화 |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [`pd_qsp_model.dot`](parkinsons-disease/pd_qsp_model.dot) | 기계론적 지도: 16 클러스터, 140+ 노드 |
| [`pd_qsp_model.svg`](parkinsons-disease/pd_qsp_model.svg) | 기계론적 지도 (SVG, 벡터) |
| [`pd_qsp_model.png`](parkinsons-disease/pd_qsp_model.png) | 기계론적 지도 (PNG, 150 dpi) |
| [`pd_mrgsolve_model.R`](parkinsons-disease/pd_mrgsolve_model.R) | mrgsolve ODE 모델: 22 구획, 7 치료 시나리오 |
| [`pd_shiny_app.R`](parkinsons-disease/pd_shiny_app.R) | Shiny 대시보드 6탭: 환자프로파일 · PK · PD바이오마커 · 운동엔드포인트 · 시나리오비교 · 신경보호 |
| [`pd_references.md`](parkinsons-disease/pd_references.md) | 참고문헌 50편 (PubMed 링크, 14개 섹션) |
| [`README.md`](parkinsons-disease/README.md) | 질환 및 모델 상세 설명 |

### 핵심 참고문헌

1. Braak H et al. "Staging of brain pathology related to sporadic Parkinson's disease." *Neurobiol Aging* 2003;24:197–211. [PMID 12498954](https://pubmed.ncbi.nlm.nih.gov/12498954/)
2. Spillantini MG et al. "Alpha-synuclein in Lewy bodies." *Nature* 1997;388:839–40. [PMID 9278044](https://pubmed.ncbi.nlm.nih.gov/9278044/)
3. Luk KC et al. "Pathological α-synuclein transmission initiates Parkinson-like neurodegeneration." *Science* 2012;338:949–53. [PMID 23161999](https://pubmed.ncbi.nlm.nih.gov/23161999/)
4. Rascol O et al. "A five-year study of the incidence of dyskinesia: ropinirole or levodopa." *N Engl J Med* 2000;342:1484–91. [PMID 10816186](https://pubmed.ncbi.nlm.nih.gov/10816186/)
5. Athauda D et al. "Exenatide once weekly versus placebo in Parkinson's disease." *Lancet* 2017;390:1664–75. [PMID 28781108](https://pubmed.ncbi.nlm.nih.gov/28781108/)

---

## 주요우울장애 (Major Depressive Disorder, MDD) — 2026-06-20

[![MDD QSP Mechanistic Map](major-depressive-disorder/mdd_qsp_model.png)](major-depressive-disorder/mdd_qsp_model.svg)

Major Depressive Disorder (MDD) is the world's leading cause of disability, affecting ~280 million people. The QSP model integrates 10 mechanistic subsystems: serotonergic, noradrenergic, dopaminergic, HPA axis, BDNF/neuroplasticity, neuroinflammation/kynurenine, glutamate/NMDA, circadian rhythm, mTOR/neurogenesis, and clinical endpoints. Drug classes modelled include SSRIs (escitalopram), SNRIs (venlafaxine), NRDIs (bupropion), NMDA antagonists (ketamine/esketamine), lithium, TCAs, MAOIs, and atypical antipsychotic augmentation.

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [`mdd_qsp_model.dot`](major-depressive-disorder/mdd_qsp_model.dot) | 기계론적 지도: 10 클러스터, 100+ 노드 |
| [`mdd_qsp_model.svg`](major-depressive-disorder/mdd_qsp_model.svg) | 기계론적 지도 (SVG, 벡터) |
| [`mdd_qsp_model.png`](major-depressive-disorder/mdd_qsp_model.png) | 기계론적 지도 (PNG, 150 dpi) |
| [`mdd_mrgsolve_model.R`](major-depressive-disorder/mdd_mrgsolve_model.R) | mrgsolve ODE 모델: 18 구획, 6 치료 시나리오 |
| [`mdd_shiny_app.R`](major-depressive-disorder/mdd_shiny_app.R) | Shiny 대시보드 6탭: 환자프로파일 · Drug PK · 신경전달물질 · HPA/신경염증 · 임상엔드포인트 · 시나리오비교 |
| [`mdd_references.md`](major-depressive-disorder/mdd_references.md) | 참고문헌 60편 (PubMed 링크, 10개 섹션) |
| [`README.md`](major-depressive-disorder/README.md) | 질환 및 모델 상세 설명 |

### 핵심 참고문헌

1. Rush AJ et al. "STAR*D outcomes in depressed outpatients." *Am J Psychiatry* 2006;163:1905. [PMID 17074942](https://pubmed.ncbi.nlm.nih.gov/17074942/)
2. Zarate CA et al. "A randomized trial of ketamine in treatment-resistant depression." *Arch Gen Psychiatry* 2006;63:856. [PMID 16894061](https://pubmed.ncbi.nlm.nih.gov/16894061/)
3. Cipriani A et al. "Comparative efficacy of 21 antidepressants." *Lancet* 2018;391:1357. [PMID 29477251](https://pubmed.ncbi.nlm.nih.gov/29477251/)
4. Duman RS & Monteggia LM. "A neurotrophic model for stress-related mood disorders." *Biol Psychiatry* 2006;59:1116. [PMID 16631126](https://pubmed.ncbi.nlm.nih.gov/16631126/)
5. Dowlati Y et al. "A meta-analysis of cytokines in major depression." *Biol Psychiatry* 2010;67:446. [PMID 20015486](https://pubmed.ncbi.nlm.nih.gov/20015486/)

---

## 만성 골수성 백혈병 (Chronic Myeloid Leukemia, CML) — 첫 번째 종양 모델

[![CML QSP Model](chronic-myeloid-leukemia/cml_qsp_model.png)](chronic-myeloid-leukemia/cml_qsp_model.svg)

### 개요

만성 골수성 백혈병(CML)은 t(9;22)(q34;q11) 전좌에 의해 형성된 **필라델피아(Ph) 염색체**와 그 산물인 **BCR-ABL1 융합 단백질(p210)**에 의해 구동되는 클론성 골수증식 종양입니다. BCR-ABL 티로신 키나아제의 구성적 활성화가 RAS-MAPK, PI3K-AKT-mTOR, JAK-STAT 경로를 통해 골수 세포의 무제한 증식과 아폽토시스 억제를 일으킵니다. TKI(티로신 키나아제 억제제) 시대에 10년 생존율은 약 85%로 개선되었으며, 지속적 깊은 분자반응(DMR) 달성 시 치료 없는 관해(TFR) 시도가 가능합니다.

### 병기 및 주요 병태생리 경로

| 병기 | 아세포 % | 특징 |
|------|---------|------|
| 만성기 (CP) | <10% | BCR-ABL IS% 높음, 치료 반응 우수 |
| 가속기 (AP) | 10-19% | 클론 진화, 호염기구 ≥20% |
| 급성전환기 (BC) | ≥20% | 골수성/림프구성, 예후 불량 |

### TKI 약물 파라미터

| TKI | 세대 | 용량 | IC50 | T315I | 주요 임상시험 |
|-----|------|------|------|-------|------------|
| 이마티닙 | 1세대 | 400mg/day | 0.25 µM | 내성 | IRIS (O'Brien 2003) |
| 다사티닙 | 2세대 | 100mg/day | 0.025 µM | 내성 | DASISION (Kantarjian 2010) |
| 닐로티닙 | 2세대 | 300mg BID | 0.020 µM | 내성 | ENESTnd (Saglio 2010) |
| 포나티닙 | 3세대 | 15-45mg/day | 0.37 nM | 감수성 | OPTIC (Cortes 2022) |
| 아시미닙 | STAMP | 40mg BID | 0.5 nM | 200mg BID | ASCEMBL (Réa 2021) |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [`cml_qsp_model.dot`](chronic-myeloid-leukemia/cml_qsp_model.dot) | 기계론적 지도: 12 클러스터, 174 노드 |
| [`cml_qsp_model.svg`](chronic-myeloid-leukemia/cml_qsp_model.svg) | 기계론적 지도 (SVG, 벡터) |
| [`cml_qsp_model.png`](chronic-myeloid-leukemia/cml_qsp_model.png) | 기계론적 지도 (PNG, 150 dpi) |
| [`cml_mrgsolve_model.R`](chronic-myeloid-leukemia/cml_mrgsolve_model.R) | mrgsolve ODE 모델: 22 구획, 7 치료 시나리오 (Michor 2005 기반) |
| [`cml_shiny_app.R`](chronic-myeloid-leukemia/cml_shiny_app.R) | Shiny 대시보드 6탭: 환자프로파일 · PK · 분자반응 · 세포역학 · 시나리오비교 · 내성/TFR |
| [`cml_references.md`](chronic-myeloid-leukemia/cml_references.md) | 참고문헌 60편 (PubMed 링크, 9개 섹션) |
| [`README.md`](chronic-myeloid-leukemia/README.md) | 질환 및 모델 상세 설명 |

### 핵심 참고문헌

1. Michor F et al. "Dynamics of chronic myeloid leukaemia." *Nature* 2005;435:1267. [PMID 15988530](https://pubmed.ncbi.nlm.nih.gov/15988530/)
2. O'Brien SG et al. (IRIS) "Imatinib vs interferon for CML." *NEJM* 2003;348:994. [PMID 12637609](https://pubmed.ncbi.nlm.nih.gov/12637609/)
3. Saglio G et al. (ENESTnd) "Nilotinib vs imatinib for CML." *NEJM* 2010;362:2251. [PMID 20525993](https://pubmed.ncbi.nlm.nih.gov/20525993/)
4. Kantarjian H et al. (DASISION) "Dasatinib vs imatinib for CML." *NEJM* 2010;362:2260. [PMID 20525994](https://pubmed.ncbi.nlm.nih.gov/20525994/)
5. Réa D et al. (ASCEMBL) "Asciminib vs bosutinib in CML." *Blood* 2021;138:2031. [PMID 34407542](https://pubmed.ncbi.nlm.nih.gov/34407542/)
6. Wylie AA et al. "The allosteric inhibitor ABL001 (asciminib)." *Nature* 2017;543:733. [PMID 28329763](https://pubmed.ncbi.nlm.nih.gov/28329763/)
7. Mahon FX et al. (STIM) "TKI discontinuation in CML." *Lancet Oncol* 2010;11:1029. [PMID 20965785](https://pubmed.ncbi.nlm.nih.gov/20965785/)
8. Hochhaus A et al. (ELN 2020) "ELN 2020 recommendations for CML." *Leukemia* 2020;34:966. [PMID 32127639](https://pubmed.ncbi.nlm.nih.gov/32127639/)

---

## 편두통 (Migraine — Trigeminovascular, CGRP, CSD)

[![Migraine QSP Model](migraine/mgr_qsp_model.png)](migraine/mgr_qsp_model.svg)

편두통은 전 세계 인구의 약 15%(약 10억 명)에 영향을 미치는 신경학적 질환으로, GBD 2016에서 두 번째로 높은 YLD 원인입니다. 삼차신경혈관 경로(TGV)와 CGRP(Calcitonin Gene-Related Peptide) 신호 경로가 핵심 발병기전이며, 피질 확산 억제(CSD)가 아우라를 설명합니다.

### 주요 병태생리 경로

| 경로 | 핵심 메커니즘 | 임상 표현 |
|------|------------|---------|
| CSD (피질 확산 억제) | K⁺ 유출 → NMDA 활성화 → 신경세포 탈분극 파동 (3–5 mm/min) | 시각 아우라(편두통 전조) |
| 삼차신경혈관 | TG 뉴런 → CGRP/SP/PACAP 방출 → 경막 혈관 확장 | 박동성 편두통 통증 |
| CGRP 신호 | CLR/RAMP1 수용체 → cAMP↑ → KATP → 혈관 확장 + 비만세포 탈과립 | 두개내압↑, 통증 강화 |
| 중추감작 | TNC wind-up → PKCε → MAPK/ERK → BDNF | 피부 이질통, 광·음과민증 |
| PGE2/COX-2 | PLA2 → AA → COX-2 → PGE2 → EP1/TRPV1 | 말초 감작, 열통증 |
| 산화질소 | NOS → NO → sGC → cGMP → KATP | 혈관확장, TG 감작 |

### 약물 PK/PD 파라미터

| 약물 | 작용기전 | 주요 PK | 임상 효과 |
|------|---------|---------|---------|
| 수마트립탄 SC 6 mg | 5-HT1B/1D 작용제 | F=97%, t½=2h, CL=72 L/h | 2h 통증소실 ~35% |
| 라스미디탄 200 mg | 5-HT1F 작용제 (혈관수축 無) | F=38%, t½=5h | 2h 통증소실 ~32% |
| 리메게판트 75 mg | CGRP-R 길항제 (게판트) | F=64%, t½=11h, Ki=0.027 nM | 2h 통증소실 ~21%; 예방 MMD -1.75일 |
| 에레누맙 140 mg SC | 항-CGRP 수용체 mAb | Tmax≈6일, t½=28일, Vd=3.86 L | MMD -3.7일, 50% 반응률 47% |
| 프레마네주맙 225 mg | 항-CGRP 리간드 mAb | t½=31일 | MMD -3.7일 |
| 갈카네주맙 120 mg | 항-CGRP 리간드 mAb | t½=27일 | MMD -4.7일 (50% 반응률 52%) |
| 토피라메이트 100 mg/d | AMPA 차단 + CA 억제 + GABA↑ | F=80%, t½=21h | MMD -2.1일, CSD 억제 60% |
| 프로프라놀롤 160 mg/d | β1/β2 차단 + 5-HT2A 차단 | F=26%, t½=4h | MMD -1.8일 |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [`mgr_qsp_model.dot`](migraine/mgr_qsp_model.dot) | Graphviz 기계론적 지도 소스 (100+ 노드, 12 클러스터) |
| [`mgr_qsp_model.svg`](migraine/mgr_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [`mgr_qsp_model.png`](migraine/mgr_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [`mgr_mrgsolve_model.R`](migraine/mgr_mrgsolve_model.R) | mrgsolve ODE 모델 (18 구획, 7 치료 시나리오) |
| [`mgr_shiny_app.R`](migraine/mgr_shiny_app.R) | Shiny 대시보드 6탭: 환자프로파일 · PK · CGRP/PD · 임상엔드포인트 · 시나리오비교 · 바이오마커 |
| [`mgr_references.md`](migraine/mgr_references.md) | 참고문헌 50편 (PubMed 링크, 5개 섹션) |
| [`README.md`](migraine/README.md) | 질환 및 모델 상세 설명 |

### 핵심 참고문헌

1. Edvinsson L et al. "CGRP as the target of new migraine therapies." *Nat Rev Neurol* 2018;14:338. [PMID 29691490](https://pubmed.ncbi.nlm.nih.gov/29691490/)
2. Goadsby PJ et al. "Pathophysiology of migraine: a disorder of sensory processing." *Physiol Rev* 2017;97:553. [PMID 28179394](https://pubmed.ncbi.nlm.nih.gov/28179394/)
3. Goadsby PJ et al. (STRIVE) "Erenumab for episodic migraine." *NEJM* 2017;377:2123. [PMID 29171821](https://pubmed.ncbi.nlm.nih.gov/29171821/)
4. Kuca B et al. (SAMURAI) "Lasmiditan for migraine." *Neurology* 2018;91:e2222. [PMID 30446598](https://pubmed.ncbi.nlm.nih.gov/30446598/)
5. Croop R et al. (ARTISAN-EM) "Rimegepant for migraine prevention." *Lancet* 2021;397:51. [PMID 33340483](https://pubmed.ncbi.nlm.nih.gov/33340483/)
6. Ashina M. "Migraine." *NEJM* 2020;383:1866. [PMID 33211930](https://pubmed.ncbi.nlm.nih.gov/33211930/)

---

## 건선 (Psoriasis — IL-23/IL-17 Axis · TNF · Biologic & Small-Molecule Therapies)

### 개요

건선은 전 세계 인구의 약 2–3%, 한국 약 1.5%에서 발생하는 만성 면역 매개 피부 질환으로,
수지상세포-Th17-IL-17A 축이 각질형성세포를 과다증식시켜 특징적인 인설성 홍반 판(plaque)을 형성합니다.
중등도–중증 기준(PASI≥10 또는 BSA≥10%)은 전체 환자의 약 20–30%를 차지하며,
생물학적 제제(IL-17A/IL-23 억제제)로 PASI90 달성률이 75% 이상에 이릅니다.

[![Psoriasis QSP Model](psoriasis/pso_qsp_model.png)](psoriasis/pso_qsp_model.svg)

### 주요 병태생리 경로

| 경로 | 핵심 메커니즘 | 임상 결과 |
|------|------------|---------|
| IL-23/Th17 축 | mDC→IL-23→STAT3/RORγt→Th17→IL-17A | 각질세포 NF-κB 활성화, 과증식 |
| IL-17A/NF-κB | IL-17RA/RC→NF-κB→KC 과증식, CXCL8, β-defensin | acanthosis, parakeratosis, PASI↑ |
| TNF-α 경로 | 대식세포/Th1→TNF-α→내피세포 활성화, 면역세포 유출 | 혈관 확장, 피부 홍반 |
| IL-22 경로 | Th17/Th22→IL-22→STAT3→KC Bcl-2↑, 세포사멸↓ | 과증식 유지 |
| 피부 장벽 결함 | Loricrin/Filaggrin↓, TJ 결함 | 항원 진입, 재발 촉진 |
| 건선관절염 | IL-17A/TNF-α → 건·부착부위 염증 | 30%에서 PsA 동반 |

### 약물 PK/PD 파라미터

| 약물 | 작용기전 | 생체이용률 | t½ | PASI75 Wk16 | 임상시험 |
|------|---------|---------|-----|------------|--------|
| 아달리무맙 40mg q2w SC | anti-TNF-α mAb | F=0.64 | 14일 | 71% | CHAMPION |
| 세쿠키누맙 300mg q4w SC | anti-IL-17A mAb | F=0.73 | 27일 | 77-80% | FIXTURE |
| 리산키주맙 150mg q12w SC | anti-IL-23p19 mAb | F=0.89 | 28일 | 88-91% | UltIMMa |
| 우스테키누맙 45mg q12w SC | anti-IL-12/23 p40 | F=0.57 | 21일 | 67-71% | PHOENIX-1 |
| 아프레밀라스트 30mg BID PO | PDE4 억제제 | F=0.73 | 9h | 33-40% | ESTEEM-1 |
| 토파시티닙 10mg BID PO | JAK1/3 억제제 | F=0.74 | 3h | 39-46% | OPT Pivotal |
| 메토트렉세이트 20mg qw PO | DHFR 억제제 | F=0.70 | 3-10h | 26-36% | Heydendael 2003 |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [`pso_qsp_model.dot`](psoriasis/pso_qsp_model.dot) | Graphviz 기계론적 지도 소스 (100+ 노드, 11 클러스터) |
| [`pso_qsp_model.svg`](psoriasis/pso_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [`pso_qsp_model.png`](psoriasis/pso_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [`pso_mrgsolve_model.R`](psoriasis/pso_mrgsolve_model.R) | mrgsolve ODE 모델 (25 구획, 7 치료 시나리오) |
| [`pso_shiny_app.R`](psoriasis/pso_shiny_app.R) | Shiny 대시보드 6탭 |
| [`pso_references.md`](psoriasis/pso_references.md) | 참고문헌 42편 (PubMed 링크) |

### 핵심 참고문헌

1. Griffiths CEM et al. "Psoriasis." *Lancet* 2021;397:1301. [PMID 33812489](https://pubmed.ncbi.nlm.nih.gov/33812489/)
2. Langley RG et al. "Secukinumab in plaque psoriasis." *N Engl J Med* 2014;371:326. [PMID 25007392](https://pubmed.ncbi.nlm.nih.gov/25007392/)
3. Gordon KB et al. "Risankizumab in moderate-to-severe plaque psoriasis (UltIMMa)." *Lancet* 2018;392:650. [PMID 30097359](https://pubmed.ncbi.nlm.nih.gov/30097359/)
4. Papp K et al. "Apremilast in moderate-to-severe psoriasis (ESTEEM-1)." *J Am Acad Dermatol* 2015;73:37. [PMID 26089047](https://pubmed.ncbi.nlm.nih.gov/26089047/)
5. Armstrong AW, Read C. "Psoriasis: a review." *JAMA* 2020;323:1945. [PMID 32427307](https://pubmed.ncbi.nlm.nih.gov/32427307/)
