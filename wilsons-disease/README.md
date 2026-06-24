# Wilson's Disease (윌슨병) QSP Model

> **디렉토리**: `wilsons-disease/` | **약어**: WD | **날짜**: 2026-06-24

[![WD QSP 기계론적 지도](wd_qsp_model.png)](wd_qsp_model.svg)

---

## 질환 개요

**윌슨병(Wilson's Disease, WD)**은 *ATP7B* 유전자의 상염색체 열성 돌연변이에 의해 발생하는 구리 대사 장애입니다. ATP7B는 간세포에서 구리를 담즙으로 배출하고 세룰로플라스민에 구리를 적재하는 P형 Cu-ATPase로, 기능 소실 시 구리가 간·뇌·각막·신장에 축적되어 진행성 간 및 신경 손상을 유발합니다.

| 특성 | 값 |
|------|-----|
| 유병률 | 1/30,000 (유전자 보인자 1/90) |
| 발병 연령 | 5–35세 (5–45세까지 보고) |
| 유전 양식 | 상염색체 열성 (AR) |
| 유전자 위치 | 13q14.3 (ATP7B, 21개 엑손) |
| 흔한 돌연변이 | p.His1069Gln (유럽, 35%), p.Arg778Leu (아시아, 20%) |
| Leipzig 진단 점수 | ≥4점 = 확진 |

---

## 병태생리 요약

```
구리 섭취 → GI 흡수 (CTR1/ATP7A) → 문맥 → 간세포 흡수
                                              ↓
                                    ATP7B 기능 ×
                                              ↓
                              ┌── 담즙 배출 ↓ ──┐
                              │                  │
                    간 Cu 축적 (>250 μg/g)    세룰로플라스민 ↓
                              │
                    Cu 과잉 → Fenton 반응 → ROS → 산화 스트레스
                              │
                    MT 포화 → NCBC 증가 → 전신 분포
                              │
                  ┌───────────┼────────────┐
                  ↓           ↓            ↓
               뇌 Cu 축적  신장 독성     각막 Cu
             (기저핵 우선)  Fanconi    KF Ring
                  ↓
          신경/정신증상 (UWDRS↑)
```

---

## 치료 약물 비교

| 약물 | 기전 | 주요 적응증 | NCBC 감소 | 주요 부작용 |
|------|------|-----------|---------|-----------|
| **D-Penicillamine** | 구리 킬레이션 → 요중 배설 | 간형 WD, 1차 치료 | ~60–70% | 신독성, SLE-유사, 신경악화 (50%) |
| **Trientine** | 구리 킬레이션 (DPA 2nd-line) | DPA 부작용 시 | ~50–60% | 위장관 불편 (경미) |
| **Zinc Acetate** | 장관 MT 유도 → Cu 흡수 차단 | 유지요법, 임신부, 소아 | ~30–40% | 위장관 자극, Fe 상호작용 |
| **ALXN1840 (TTM)** | TTM-Cu-Albumin 삼중복합체 형성 | 신경형 WD (1차), 간형 | **~98%** | 경미 (ATLAS trial) |

---

## QSP 모델 파일

| 구성요소 | 파일 | 사양 |
|---------|------|-----|
| 🗺️ 기계론적 지도 | [`wd_qsp_model.dot`](wd_qsp_model.dot) | **119 노드, 11 클러스터** |
| 🖼️ SVG | [`wd_qsp_model.svg`](wd_qsp_model.svg) | 벡터 포맷 |
| 🖼️ PNG | [`wd_qsp_model.png`](wd_qsp_model.png) | 150 DPI |
| ⚙️ mrgsolve ODE | [`wd_mrgsolve_model.R`](wd_mrgsolve_model.R) | **24구획 ODE**, **8치료 시나리오** |
| 📊 Shiny 앱 | [`wd_shiny_app.R`](wd_shiny_app.R) | **8탭** 인터랙티브 대시보드 |
| 📚 참고문헌 | [`wd_references.md`](wd_references.md) | **60개 문헌** (13섹션) |

---

## 기계론적 지도 클러스터 (11개)

| # | 클러스터 | 주요 노드 |
|---|---------|---------|
| 1 | GI 구리 흡수 | CTR1, DMT1, MT_gut, ATP7A, 문맥 |
| 2 | 간세포 구리 샤페론 | ATOX1, CCS, COX17, SOD1, MT_hepatic |
| 3 | ATP7B & 담즙 배출 | ATP7B_WT/Mutant, TGN, 담즙 소포, 세룰로플라스민 |
| 4 | 전신 구리 & 바이오마커 | Cp_serum, NCBC, 24h 요중 Cu, Leipzig 점수 |
| 5 | 간 병태 & 섬유화 | Cu-ROS, Kupffer, HSC, TGF-β, 섬유화, ALT/AST |
| 6 | 뇌 구리 & 신경정신과 | 기저핵, 도파민신경, 진전/근긴장이상증, UWDRS |
| 7 | 다장기 독성 | KF Ring, Fanconi 증후군, 용혈성 빈혈 |
| 8 | DPA PK/PD | 흡수→킬레이션→요중 배설, 신경악화 역설 |
| 9 | Zinc PK/PD | 장관 MT 유도, Cu 흡수 차단 |
| 10 | Trientine & ALXN1840 PK/PD | 삼중복합체, 분변 배설 |
| 11 | 임상 결과 | Leipzig 점수, 치료 반응, 간이식 적응증 |

---

## mrgsolve ODE 모델 구획 (24개)

```
[약물 PK — 8구획]
GUT_DPA → CENT_DPA           D-Penicillamine 2구획 PK
GUT_ZN  → CENT_ZN            Zinc 2구획 PK
GUT_TRI → CENT_TRI           Trientine 2구획 PK
GUT_TTM → CENT_TTM           ALXN1840 2구획 PK

[구리 동역학 — 7구획]
CU_GI                        GI 구리 (흡수 풀)
CU_HEP                       간 구리 (핵심 축적 부위)
MT_HEP                       간 메탈로티오닌-결합 구리
CU_NCBC                      NCBC (자유 구리, 독성 분획)
CP_SERUM                     혈청 세룰로플라스민
CU_URINE                     요중 구리 배설률
CU_BRAIN / CU_KIDNEY / CU_CORNEA  장기 구리

[간 병태 — 3구획]
ROS_HEP                      간 ROS 지수
ALT_SERUM                    혈청 ALT
FIBROSIS                     섬유화 점수 (Metavir F0–F4)

[신경 — 1구획]
NEURODEGENERATION            신경퇴행 지수
```

---

## 치료 시나리오 (8개)

| 시나리오 | 약물 | 임상 근거 | 핵심 결과 |
|---------|------|---------|---------|
| S1 | 무치료 WD | 자연경과 | 간 Cu ↑↑, 섬유화 진행, 신경퇴행 |
| S2 | DPA 500mg TID | Walshe 1956 | NCBC ↓60%, 초기 신경악화 가능 |
| S3 | Zinc 50mg TID | Brewer 1998 | 흡수 차단, NCBC ↓30–40% |
| S4 | Trientine 500mg TID | Weiss 2013 | NCBC ↓50%, DPA보다 부작용 적음 |
| S5 | ALXN1840 15mg QD | ATLAS 2022 | **NCBC ↓98%** |
| S6 | DPA→Zinc 전환 (1년 후) | AASLD 가이드라인 | 초기 킬레이션 후 유지 전략 |
| S7 | ALXN1840 + Trientine 병용 | 가상 탐색 | NCBC 극대 감소 |
| S8 | 정상 WT 대조 | 참조 | 정상 Cp, NCBC<10, ALT 정상 |

---

## Shiny 앱 탭 구성 (8탭)

1. **환자 프로파일**: Leipzig 점수 계산기, ATP7B 돌연변이 분포, 유전형-표현형
2. **약물 PK**: DPA/Zinc/Trientine/ALXN1840 혈중 농도 시뮬레이션
3. **구리 동역학**: 혈청 구리·NCBC·세룰로플라스민·장기 분포 시계열
4. **간 결과**: ALT·섬유화·ROS 경과, 간이식 위험 지수
5. **신경/안과 결과**: 뇌 구리·UWDRS·KF Ring 경과, DPA 역설적 악화
6. **시나리오 비교**: 8개 시나리오 병렬 비교, ATLAS 임상시험 대비
7. **바이오마커 탐색기**: 구리 바이오마커 상관관계, 진단 패널
8. **모델 정보**: 파라미터 출처, 임상시험 요약

---

## 참고문헌 주요 인용

- Członkowska A, et al. *Nat Rev Dis Primers* 2018 (**포괄적 질환 리뷰**)
- Bandmann O, et al. *Lancet Neurol* 2015 (**신경병 기전**)
- Schilsky ML, et al. *NEJM Evid* 2022 (**ATLAS Trial — ALXN1840**)
- EASL. *J Hepatol* 2012 (**진단·치료 가이드라인**)
- Brewer GJ, et al. *J Lab Clin Med* 1998 (**Zinc 15년 추적**)
