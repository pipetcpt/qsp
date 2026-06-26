# 쿠싱 증후군 (Cushing's Syndrome) — QSP 모델

> **디렉토리:** `cushings-syndrome/` | **약어:** CS | **날짜:** 2026-06-24  
> **분류:** 내분비 질환 (Endocrine Disease)

[![CS QSP 기계론적 지도](cs_qsp_model.png)](cs_qsp_model.svg)

---

## 질환 개요

**쿠싱 증후군(Cushing's Syndrome)**은 내인성 또는 외인성 글루코코르티코이드 과잉에 의해 발생하는 내분비 질환입니다.

| 항목 | 내용 |
|------|------|
| 유병률 | 100만 명당 10-15명/년 (내인성) |
| 가장 흔한 원인 | 쿠싱병(뇌하수체 ACTH 선종) — 70% |
| 기타 원인 | 이소성 ACTH(10%), 부신 선종(15%), 부신암(5%) |
| 표준 치료 | 경접형동 수술 (뇌하수체 선종) / 부신절제술 (부신 원인) |
| 2차 치료 | 파시레오티드, 케토코나졸, 메티라폰, 오실로드로스탯, 미페프리스톤 |
| 핵심 생화학 지표 | UFC >150 μg/24h, LNSC >10 nmol/L, 1mg DST >1.8 μg/dL |

---

## 핵심 기전 (13개 클러스터, 140+ 노드)

| 클러스터 | 핵심 기전 |
|---------|-----------|
| 1. 시상하부 | CLOCK/BMAL1 일주기 리듬, CRH/AVP 분비, 소마토스타틴·도파민 억제 |
| 2. 뇌하수체 전엽 | CRHR1 → Gs → PKA → CREB → POMC → ACTH; USP8 돌연변이(~50% 쿠싱병); CDK4/6-Rb 증식 |
| 3. 부신피질 생합성 | MC2R → cAMP → PKA → StAR → CYP11A1 → ... → CYP11B1 → 코르티솔 |
| 4. GR 신호 | GR-α/β, HSP90, FKBP51/52, 핵이동, GRE/nGRE, AP-1/NF-κB 접촉억제, GILZ, MKP1, SGK1 |
| 5. 대사 합병증 | PEPCK/G6Pase↑(간), GLUT4↓(말초), 인슐린저항성, 내장지방↑, 근육위축, 골다공증 |
| 6. 심혈관·신장 | RAAS 과활성(Ang II→AT1R→알도스테론), 나트륨저류, 고혈압, VTE, 이상지질혈증 |
| 7. 면역 억제 | NF-κB/AP-1 억제→사이토카인↓, 림프구감소증, NK세포↓, 감염위험↑ |
| 8. CNS 효과 | 해마위축, BDNF↓, NMDA흥분독성, 세로토닌↓, 우울증·인지기능저하 |
| 9. 파시레오티드 PK/PD | 2구획 SC; SSTR5>SSTR2; Gi→cAMP↓→ACTH억제(최대 65%) |
| 10. 스테로이드 합성 억제제 | 케토코나졸(CYP17A1+11B1), 메티라폰(CYP11B1), 오실로드로스탯(CYP11B1/B2) |
| 11. GR 길항제·카버골린 | 미페프리스톤(GR 경쟁적 길항), 카버골린(D2R→ACTH억제), 미토탄(부신독성) |
| 12. 임상 진단·평가 | UFC, LNSC, 1mg DST, 고용량 DST, CRH 자극 검사, IPSS |
| 13. 병인 분류 | 쿠싱병, 이소성 ACTH, 부신 선종/암, PBMAH, 맥쿤-올브라이트, 주기성 쿠싱 |

---

## mrgsolve ODE 모델 (21구획)

| 모듈 | 구획 | 핵심 동역학 |
|------|------|------------|
| HPA 축 | CRH, ACTH_PIT, ACTH_PL | 일주기 CRH 합성 + 종양 ACTH 분비 + GR 음성 피드백 (Hill n=2) |
| 부신 스테로이드 | F_ADR, F_PL | Michaelis-Menten (ACTH→코르티솔); 약물 Bliss 결합 억제 |
| GR 동역학 | GR_FREE, GR_BOUND, GR_NUC | 2-단계(세포질→핵) 결합/해리 ODE |
| 대사 | GLUCOSE, INSULIN, VAT, MUSCLE, BMD, BP | GR 구동 대사 변화; 인슐린 피드백 |
| 임상 출력 | UFC_ACC, LNSC | 비례 분비 모델 |
| 파시레오티드 PK | A_PAS_C, A_PAS_P | 2구획; SC 흡수; Emax=0.62 |
| 케토코나졸 PK | A_KETO | 1구획; CYP17A1+11B1 이중 억제; Emax=0.72 |
| 메티라폰 PK | A_METY | 1구획; CYP11B1 선택 억제; Emax=0.82 |
| 오실로드로스탯 PK | A_OSILO | 1구획; 고 효능(EC50=0.15 μg/mL); Emax=0.85 |
| 미페프리스톤 PK | A_MIFE | 1구획; 고 Vd (115 L); GR 점유 Emax=0.82 |

---

## 치료 시나리오 임상 보정 데이터

| 시나리오 | 약물 | 임상시험 | 주요 결과 |
|---------|------|----------|----------|
| 자연경과 (무치료) | — | 역사적 코호트 | UFC >500 μg/24h, 합병증 누적 |
| 파시레오티드 0.6mg BID | 파시레오티드 | PASPORT-CUSHINGS (Colao 2012 NEJM) | UFC 정상화 22-24% (6개월) |
| 케토코나졸 400mg BID | 케토코나졸 | Castinetti 2014 Eur J Endocrinol | UFC 정상화 49%, 간독성 주의 |
| 오실로드로스탯 5mg BID | 오실로드로스탯 | LINC 3/4 (Pivonello 2020 Lancet DE) | UFC 정상화 86% (유지기) |
| 미페프리스톤 600mg QD | 미페프리스톤 | SEISMIC (Fleseriu 2012 JCEM) | 혈당/BP 임상반응 87% (GR 길항) |
| 수술 후 관해 | — | 메타분석 | 재발률 ~20% (5년 내) |

---

## QSP 모델 구성 파일

| 구성요소 | 파일 | 사양 |
|---------|------|-----|
| 🗺️ 기계론적 지도 (DOT) | [`cs_qsp_model.dot`](cs_qsp_model.dot) | **140+ 노드, 13 클러스터** |
| 🖼️ SVG 벡터 이미지 | [`cs_qsp_model.svg`](cs_qsp_model.svg) | 확대 가능 고해상도 |
| 🖼️ PNG 래스터 이미지 | [`cs_qsp_model.png`](cs_qsp_model.png) | 150 dpi |
| ⚙️ mrgsolve ODE 모델 | [`cs_mrgsolve_model.R`](cs_mrgsolve_model.R) | **21구획 ODE, 6치료 시나리오** |
| 📊 Shiny 앱 | [`cs_shiny_app.R`](cs_shiny_app.R) | **8탭 대시보드** |
| 📚 참고문헌 | [`cs_references.md`](cs_references.md) | **55개 PubMed 인용 (10개 섹션)** |

---

## Shiny 앱 탭 구성 (8탭)

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 기저 코르티솔/ACTH/혈당 설정, 병인 선택, 진단 기준 |
| 2. HPA축/PK 동역학 | 일주기 리듬, ACTH 동역학, 코르티솔 시계열, 약물 PK |
| 3. 스테로이드 생합성 | 효소 억제 동역학, 경로 개요, 코르티솔 합성 억제 |
| 4. 임상 지표 | UFC, LNSC, 관해 기준 판정, 덱사메타손 억제 검사 |
| 5. 시나리오 비교 | 5가지 치료 병렬 비교 (코르티솔/ACTH/혈당/요약 표) |
| 6. 바이오마커 패널 | GR 핵 활성화, HPA 피드백, 덱사메타손 검사 시뮬레이션 |
| 7. 대사 합병증 | 혈당/인슐린, 체성분(VAT+근육), BMD, 혈압 |
| 8. 가상 집단 분석 | N=10-500 가상 환자, 반응률, UFC 분포 |

---

## 진단 바이오마커 참고값

| 바이오마커 | 정상 | 쿠싱 증후군 | 단위 |
|-----------|------|------------|------|
| UFC 24h | < 50 | > 150 (흔히 >500) | μg/24h |
| LNSC (자정 타액) | < 4 | > 10 | nmol/L |
| 1mg DST 코르티솔 | < 1.8 | > 1.8 (미억제) | μg/dL |
| 혈장 ACTH | 10–46 | ↑(뇌하수체/이소성) / ↓(부신) | pg/mL |
| 공복혈당 | < 5.6 | 5.6–11.1+ | mmol/L |
| BMD T-score | > -1.0 | < -1.0 (흔히 < -2.5) | |
| 수축기 혈압 | < 140 | > 140 (흔히 150-180) | mmHg |

---

## 참조 문헌 (핵심 5편)

1. Colao et al. (2012) *N Engl J Med* — 파시레오티드 Phase III
2. Pivonello et al. (2020) *Lancet Diabetes Endocrinol* — 오실로드로스탯 LINC 3
3. Nieman et al. (2018) *J Clin Endocrinol Metab* — 치료 가이드라인
4. Fleseriu et al. (2012) *J Clin Endocrinol Metab* — 미페프리스톤 SEISMIC
5. Miller & Auchus (2011) *Endocr Rev* — 스테로이드 생합성 기전

---

*모델 생성: 2026-06-24 | Claude Code Routine | QSP Library [pipetcpt/qsp](https://github.com/pipetcpt/qsp)*
