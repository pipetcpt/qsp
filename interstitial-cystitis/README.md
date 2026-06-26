# 간질성 방광염 / 방광 통증 증후군 (Interstitial Cystitis / Bladder Pain Syndrome, IC/BPS)

> **디렉토리:** `interstitial-cystitis/` | **약어:** IC/BPS | **날짜:** 2026-06-26  
> **유병률:** 여성의 2–7%, 남성의 0.5%에서 발생; 미국 300–800만 명 추정 이환  
> **핵심 키워드:** GAG layer · mast cells · neurogenic inflammation · central sensitization · TRPV1 · bladder pain

---

## 질환 개요

**간질성 방광염 / 방광 통증 증후군(IC/BPS)** 은 요로 감염이나 다른 명백한 원인 없이 방광 통증, 절박뇨, 빈뇨를 특징으로 하는 만성 비뇨기 질환입니다. 삶의 질에 심각한 영향을 미치며, 다양한 기전이 복합적으로 작용하는 이질적(heterogeneous) 질환입니다.

### 주요 아형 (Subtypes)

| 특징 | **Hunner형 (염증형)** | **비Hunner형 (다수)** |
|------|----------------------|----------------------|
| 비율 | 5–15% | 85–95% |
| 방광경 소견 | Hunner 궤양 (mucosal lesion) | 점상출혈, 사구체화 |
| 주요 병리 | 림프구/형질세포 침윤, IgG4+ 세포 | GAG층 결핍, 비만세포 침윤 |
| 염증지표 | IL-6, IFN-γ 매우 높음 | 중등도 |
| 권장 치료 | CsA, 전기소작술 | PPS, BoNTA, 방광 내 주입치료 |

---

## 기계론적 경로 요약

```
트리거 (스트레스 · 감염 · 자가면역 · 유전 · 알레르겐)
  │
  ▼
GAG층 결핍 (헤파란황산 · 콘드로이틴황산 · 히알루론산 결핍)
  │  → 요로상피 투과성 ↑ → 요중 K+ 누출 → 신경 탈분극
  │
  ├─ 비만세포 활성화 → 히스타민 · 트립타아제 · PGE2 · TNF-α
  │    └─ PAR2 활성화 → 요로상피 장벽↓ (정피드백)
  │
  ├─ C-섬유 (TRPV1+, P2X3+) 감작
  │    └─ SP · CGRP 역행성 분비 → 신경성 염증 → 비만세포 추가 활성화
  │
  ├─ NGF ↑ → TrkA → C-섬유 증식 + TRPV1↑
  │
  ├─ 척수 감작 (Wind-up · NMDA · BDNF)
  │    └─ 중추 감작 (ACC · 도서엽 · 전전두엽) → 통증 증폭
  │
  └─ 방광벽 리모델링
       ├─ TGF-β1 → 콜라겐 침착 → 섬유화
       └─ 방광 용량↓ → 빈뇨 · OLS 점수↑
```

---

## 기계론적 지도 (Mechanistic Map)

[![IC/BPS QSP Mechanistic Map](ic_bps_qsp_model.png)](ic_bps_qsp_model.svg)

*클릭하면 고해상도 SVG 지도로 이동합니다 (140+ 노드, 9 클러스터)*

---

## mrgsolve ODE 모델 (`ic_bps_mrgsolve_model.R`)

### 구획 구성 (22 ODEs)

| 범주 | 구획 | 설명 |
|------|------|------|
| **약물 PK (8개)** | PPS_GUT, PPS_CENT | 펜토산 폴리설페이트 (F=6%) |
| | HYD_GUT, HYD_CENT | 히드록시진 (F=80%) |
| | CSA_GUT, CSA_CENT | 시클로스포린 A (F=35%) |
| | AMI_GUT, AMI_CENT | 아미트립틸린 (F=50%) |
| **질환 PD (14개)** | GAG | GAG층 완전성 (0–1) |
| | PERM | 요로상피 투과도 (0–1) |
| | MC | 비만세포 활성 지수 |
| | HIST | 히스타민 수준 |
| | SP | 서브스턴스 P |
| | NGF | 신경성장인자 |
| | IL6 | 인터루킨-6 |
| | TNF | TNF-α |
| | C_FIBER | C-섬유 감작 지수 |
| | SPINAL | 척수 감작 지수 |
| | CENTRAL | 중추 감작 지수 |
| | CAP | 기능적 방광 용량 (mL) |
| | PAIN | VAS 통증 점수 (0–10) |
| | OLS | O'Leary-Sant 점수 (0–20) |

### 치료 시나리오 7가지

| # | 요법 | 용량 | 근거 임상시험 |
|---|------|------|--------------|
| S1 | 무치료 (자연경과) | — | 자연 경과 관찰 |
| S2 | **PPS (Elmiron) 구강** | 100mg TID | Nickel 2005 *Urology*; Hanno 2003 *Urology* |
| S3 | **히드록시진** | 25mg QD | Sant 2003 *J Urol*; Theoharides 1991 |
| S4 | **방광 내 DMSO** | 50% 50mL, q2주×6회 | Sant 1987 *Urology*; Fowler 1981 |
| S5 | **시클로스포린 A** | 3mg/kg/day (Hunner형) | Sairanen 2005 *J Urol*; Forrest 2012 |
| S6 | **보툴리눔독소 A** | 100U 방광 내 주입 | Gottsch 2011 *J Urol*; Kuo 2010 |
| S7 | **삼중 병용요법** | PPS + 히드록시진 + 아미트립틸린 | Foster 2010 *J Urol* |

### 주요 파라미터 근거

- **GAG 합성/분해 상수**: Parsons 2007 *Urology* — GAG 결핍 기전
- **비만세포 활성화 속도**: Peeker 2000 *J Urol* — mast cell density IC/BPS
- **C-섬유 감작**: Nazif 2007 *Urology* — neural upregulation in IC/BPS
- **PPS 생체이용률 6%**: Nickel 2005 — PPS PK 데이터
- **CsA 반감기 24h**: Forrest 2012 — calcineurin inhibitor IC/BPS

---

## Shiny 대시보드 (`ic_bps_shiny_app.R`)

### 탭 구성 (8탭)

| 탭 | 내용 |
|----|------|
| 1. **Overview** | 질환 소개, QSP 구조, 아형 비교표, value boxes |
| 2. **Patient Profile** | 기저 바이오마커 레이더 차트, UPOINT 도메인 점수, 아형 분류 |
| 3. **PK — Drug Levels** | 약물 농도-시간 프로파일, PK 파라미터 표 |
| 4. **PD — Biomarkers** | 개별 PD 바이오마커 동태, 히트맵(기저 vs 치료 후) |
| 5. **Clinical Endpoints** | VAS 통증, OLS 점수, 방광 용량, 배뇨 빈도 |
| 6. **Scenario Comparison** | 7개 시나리오 동시 비교, 요약표 |
| 7. **Subtype Explorer** | Hunner vs 비Hunner 면역 프로파일, CsA 반응 비교 |
| 8. **Sensitivity Analysis** | 1-way 민감도 분석, 토네이도 다이어그램 |

### 실행 방법

```r
# R 패키지 설치
install.packages(c("shiny", "shinydashboard", "mrgsolve", "dplyr",
                   "ggplot2", "tidyr", "DT", "plotly", "patchwork"))

# Shiny 앱 실행
shiny::runApp("ic_bps_shiny_app.R")
```

---

## 참고문헌 요약 (`ic_bps_references.md`)

총 **62편** PubMed 인용 — 17개 섹션 분류:

1. 임상지침 및 질환 정의 (4편)
2. 역학 및 유병률 (3편)
3. GAG층 결핍 및 요로상피 장벽 (5편)
4. 비만세포 병태생리 (4편)
5. 신경성 염증 및 감각 경로 (6편)
6. 면역병리 및 아형 (5편)
7. 중추 감작 및 통증 기전 (4편)
8. 진단 및 결과 지표 (4편)
9. 치료: PPS (3편)
10. 치료: 방광 내 DMSO (3편)
11. 치료: 히드록시진, 아미트립틸린, 다중모달 (4편)
12. 치료: 시클로스포린 A (3편)
13. 치료: 보툴리눔독소 A (3편)
14. 치료: 신경조절술 (2편)
15. 신규/연구 중 치료 (4편)
16. 동반질환 및 전신 연관성 (3편)
17. QSP 모델링 방법론 (2편)

---

## 파일 목록

| 파일 | 설명 |
|------|------|
| [`ic_bps_qsp_model.dot`](ic_bps_qsp_model.dot) | Graphviz 기계론적 지도 (140+ 노드, 9 클러스터) |
| [`ic_bps_qsp_model.svg`](ic_bps_qsp_model.svg) | 벡터 형식 지도 (고해상도) |
| [`ic_bps_qsp_model.png`](ic_bps_qsp_model.png) | 래스터 형식 지도 (150 dpi) |
| [`ic_bps_mrgsolve_model.R`](ic_bps_mrgsolve_model.R) | mrgsolve ODE QSP 모델 (22구획, 7 시나리오, VP n=200) |
| [`ic_bps_shiny_app.R`](ic_bps_shiny_app.R) | Shiny 인터랙티브 대시보드 (8탭) |
| [`ic_bps_references.md`](ic_bps_references.md) | 참고문헌 62편 (17개 섹션) |
| [`README.md`](README.md) | 이 문서 |

---

## 임상 적용 노트

### 약물 선택 가이드 (Evidence-Based)

```
IC/BPS 진단 확정
  │
  ├─ 방광경 → Hunner 병변 확인?
  │    ├─ YES (Hunner형)
  │    │    → 전기소작/레이저 소작술 (1차)
  │    │    → CsA 3mg/kg/day (재발 방지)
  │    │    → 스테로이드 + 면역억제
  │    │
  │    └─ NO (비Hunner형)
  │         → 1단계: 행동치료 + 식이관리 + 물리치료
  │         → 2단계: PPS 100mg TID + 히드록시진 25–50mg QD
  │         → 3단계: 방광 내 주입 (DMSO/헤파린/리도카인)
  │         → 4단계: 아미트립틸린 25–75mg QD (중추 감작)
  │         → 5단계: BoNTA 100U 방광 내 주입 (q6개월)
  │         → 6단계: 신경조절술 (천골신경조절술)
  └─
```

### 치료 반응 예측 인자 (QSP 모델)

| 인자 | 반응 좋은 군 | 반응 나쁜 군 |
|------|------------|------------|
| GAG층 지수 | > 50% | < 30% |
| 비만세포 활성 | 중등도 | 매우 높음 |
| 중추 감작 지수 | 낮음 | 높음 (전신 통증 동반) |
| 아형 | 비Hunner (PPS) | Hunner (CsA 필요) |
| 방광 용량 | > 200 mL | < 150 mL (BoNTA 우선) |

---

*모델 버전 1.0 · Claude Code Routine · 2026-06-26*
