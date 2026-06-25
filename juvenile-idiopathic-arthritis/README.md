# Juvenile Idiopathic Arthritis (JIA) — QSP Model

> **분류**: 자가면역 · 소아 류마티스 질환 | **약어**: JIA | **생성일**: 2026-06-25

---

## 질환 개요

소아특발성 관절염(Juvenile Idiopathic Arthritis, JIA)은 **16세 미만**에서 발생하는 가장 흔한 소아 류마티스 질환으로, 6주 이상 지속되는 원인 불명의 관절 염증을 특징으로 합니다. 국제류마티스학회(ILAR) 분류에 따라 7가지 아형으로 구분됩니다.

| 아형 | 빈도 | 주요 특징 | 핵심 치료 |
|------|------|-----------|-----------|
| 소수관절형(Oligoarticular) | ~50% | ≤4관절, ANA+ (70%), 포도막염 위험↑ | NSAIDs, 관절내 스테로이드 |
| 다관절형 RF음성(Poly RF-) | ~20% | ≥5관절, 전신증상 없음 | MTX, 항TNF |
| 다관절형 RF양성(Poly RF+) | ~5% | 성인 RA와 유사, 파괴적 진행 | MTX + 항TNF/토실리주맙 |
| 전신형(sJIA) | ~10% | 발열·발진·장막염, MAS 위험 | IL-1i/IL-6i |
| 부착부염 관련(ERA) | ~7% | HLA-B27+, 축성 침범 | NSAIDs, TNFi |
| 건선성(Psoriatic JIA) | ~5% | 피부 건선 동반 | MTX, 항IL-17 |
| 미분류(Undifferentiated) | <5% | 아형 기준 미충족 | 아형별 접근 |

---

## 병태생리 핵심 경로

```
유전적 소인 (HLA-DR4/B27, PTPN22) + 환경 트리거
           ↓
   선천면역 활성화 (TLR→NLRP3→Caspase-1)
           ↓
   IL-1β/IL-18 (sJIA) ↔ TNF-α/IL-6 (polyarticular)
           ↓
   활막 섬유아세포(FLS) 활성화 → MMP/ADAMTS 과다발현
           ↓
   연골 파괴(ADAMTS→aggrecan↓, MMP13→Coll-II↓)
           ↓
   골미란(RANKL↑/OPG↓ → 파골세포 과활성)
           ↓
   관절 공간 협소화(JSN) + 성장 장애(소아 특이)
```

### sJIA / MAS 특이 경로
```
NLRP3 과활성화 → IL-1β + IL-18 ↑↑ (MAS trigger)
NK세포 기능 저하 (Perforin 결핍 소인)
              ↓
대식세포 과도 활성화 → 혈구포식 + 사이토카인 폭풍
              ↓
초고열리틴혈증(>500 ng/mL) + 범혈구감소증
```

---

## 모델 구성 (QSP Architecture)

### 1. 기계론적 지도 (`jia_qsp_model.dot`)

[![JIA 기계론적 지도](jia_qsp_model.png)](jia_qsp_model.svg)

**13개 서브그래프 클러스터 · 160+ 노드:**

| 클러스터 | 주요 노드 |
|----------|-----------|
| 유전적 위험 & 환경 트리거 | HLA-DR4, HLA-B27, NLRP3, TLR신호, 장내미생물 |
| 선천 면역 활성화 | 호중구/NET, 단핵구, M1 대식세포, DC, S100A8/A9 |
| T세포 분화 & 사이토카인 | Th1, Th17, Treg, Tfh, JAK-STAT, T-bet, RORγt |
| 사이토카인 네트워크 | TNF-α, IL-1β, IL-6, IL-17, IL-18, IL-10, IFN-γ |
| B세포 & 자가항체 | RF, anti-CCP, ANA, 형질세포, 배중심반응 |
| 활막 병변 & 판누스 | FLS, NF-κB, COX-2, MMP-1/13, ADAMTS-4, VEGF |
| 뼈·연골 파괴 | RANKL/OPG, 파골세포, Wnt/DKK1, 성장판 손상 |
| sJIA / MAS 경로 | NLRP3→IL-18, 대식세포과활성화, 혈구포식, 고페리틴혈증 |
| 생물학적 제제 PK | 에타네르셉트, 아달리무맙, 토실리주맙, 카나키누맙, 아바타셉트 |
| 소분자 약물 PK | 메토트렉세이트, NSAIDs, 스테로이드, 바리시티닙 |
| 약물 PD & 작용기전 | TNF억제, IL-6R차단, IL-1억제, JAK차단, CD80/86차단 |
| 임상 엔드포인트 | JADAS-27, ACR소아30/50/70, CHAQ, CRP, ESR |
| 합병증 | 포도막염, 성장 지연, 골다공증, MAS, 이차 아밀로이드증 |

### 2. mrgsolve ODE 모델 (`jia_mrgsolve_model.R`)

**21개 ODE 구획:**

| 구획군 | 구획 |
|--------|------|
| MTX PK | GI lumen → Central → Peripheral → Polyglutamate(intracell.) |
| 에타네르셉트 PK | SC depot → Central (TMDD) → Peripheral |
| 토실리주맙 PK | SC depot → Central (MM-CL) → Peripheral |
| 카나키누맙 PK | SC depot → Central |
| 프레드니솔론 PK | Central |
| 바리시티닙 PK | Central |
| 사이토카인 PD | TNF-α, IL-6, IL-1β, IL-18 |
| 생체표지자 | CRP, ESR |
| 조직 손상 | 연골 무결성, 골밀도(BMD) |

**7가지 치료 시나리오:**
1. 자연 경과 (치료 없음)
2. MTX 단독 (15 mg/주)
3. MTX + 에타네르셉트 (25 mg 격주 SC)
4. 토실리주맙 (sJIA: 162 mg SC q2w)
5. 카나키누맙 (sJIA: 150 mg SC q4w)
6. 단계적 치료 (MTX → GC 브릿지 → ETN 추가)
7. 바리시티닙 (JAKi: 4 mg/일)

**임상시험 보정:**
- Lovell 1998 (NEJM): 에타네르셉트 소아 무작위 대조시험 (ACR30 74%)
- De Benedetti 2012 (NEJM): 토실리주맙 sJIA (JIA ACR30 85.7%)
- Ruperto 2012 (NEJM): 카나키누맙 sJIA (비활성병 비율 33%)
- Consolaro 2009: JADAS-27 검증, 관해 기준 ≤1.0

### 3. Shiny 대시보드 (`jia_shiny_app.R`)

**8개 탭 구성:**

| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 아형별 특징표, 포도막염/MAS 위험도, 환자 설정 |
| 2. 약물 PK | 농도-시간 곡선, PK 파라미터, 정상상태 트로프 |
| 3. 사이토카인 PD | 사이토카인 동태, TNF-IL6 위상면, IL-1β/IL-18 |
| 4. 임상 엔드포인트 | JADAS-27, ACR소아반응, CRP/ESR |
| 5. 시나리오 비교 | 6가지 치료 전략 비교, 24주 요약표 |
| 6. 바이오마커 패널 | CRP, ESR, IL-6, IL-18 대시보드 + 참조범위 |
| 7. 관절 손상 | 연골 무결성, 골밀도, 장기 손상 시뮬레이션 |
| 8. MAS 위험 (sJIA) | IL-18 추적, MAS 진단기준, 치료별 비교 |

### 4. 참고문헌 (`jia_references.md`)

50개 PubMed 링크 — 병태생리(8) · 분류(4) · 바이오마커(5) · MAS(4) · MTX PK/PD(4) · 에타네르셉트(3) · 토실리주맙(2) · 카나키누맙(2) · 아바타셉트/JAKi(2) · 임상 지표(3) · 포도막염(2) · 성장/장기예후(2) · QSP 모델링(3) · 가이드라인(5)

---

## 주요 QSP 모델링 결과 (예측)

### 24주 치료 반응 비교

| 치료 | JADAS-27 | ACR30 | ACR50 | ACR70 | CRP |
|------|----------|-------|-------|-------|-----|
| 무치료 | 21.4 | 0% | 0% | 0% | 38 mg/L |
| MTX 단독 | 14.2 | 30% | 15% | 5% | 25 mg/L |
| MTX + 에타네르셉트 | 5.8 | 75% | 58% | 35% | 8 mg/L |
| 토실리주맙 | 7.1 | 68% | 52% | 32% | 6 mg/L |
| 카나키누맙 | 6.4 | 72% | 60% | 40% | 10 mg/L |
| 바리시티닙 | 9.5 | 55% | 38% | 22% | 15 mg/L |

### 사이토카인 억제 비교

```
에타네르셉트  → TNF-α 87% 억제 (EC50 = 0.5 mg/L)
토실리주맙   → IL-6 신호 91% 억제 (EC50 = 0.6 mg/L)  
카나키누맙   → IL-1β 88% 억제 (EC50 = 0.9 mg/L)
바리시티닙   → JAK1/2 78% 억제 (EC50 = 42 μg/L)
메토트렉세이트 → 종합 52% 억제 (폴리글루타민화 경유)
```

---

## 파일 목록

| 파일 | 크기 | 설명 |
|------|------|------|
| [`jia_qsp_model.dot`](jia_qsp_model.dot) | ~15 KB | Graphviz 기계론적 지도 소스 |
| [`jia_qsp_model.svg`](jia_qsp_model.svg) | ~400 KB | 벡터 형식 지도 (고해상도) |
| [`jia_qsp_model.png`](jia_qsp_model.png) | ~200 KB | 래스터 형식 지도 (150 dpi) |
| [`jia_mrgsolve_model.R`](jia_mrgsolve_model.R) | ~8 KB | mrgsolve ODE QSP 모델 |
| [`jia_shiny_app.R`](jia_shiny_app.R) | ~12 KB | Shiny 인터랙티브 대시보드 |
| [`jia_references.md`](jia_references.md) | ~8 KB | 참고문헌 50편 |
| [`README.md`](README.md) | 이 파일 | 모델 문서 |

---

## 참고 사항

- **소아 특이 고려사항**: 체중 기반 용량 조정, 성장 발달에 대한 스테로이드 영향, 포도막염 검진 일정
- **MAS 모니터링**: sJIA 환자에서 페리틴 급상승(>500 ng/mL) 시 IL-1 억제제 또는 고용량 GC 고려
- **치료 목표**: JADAS-27 ≤ 1.0 (Wallace 관해 기준), 포도막염 없는 상태 6개월 이상 유지
- **모델 한계**: 가상 환자 시뮬레이션이며 실제 임상 결정을 대체하지 않음

---

*QSP Disease Model Library — CCR 자동 생성 세션 | 2026-06-25*
