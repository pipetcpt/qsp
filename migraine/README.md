# Migraine QSP Model

## 개요 (Overview)

편두통(Migraine)은 전 세계 인구의 약 15%(약 10억 명)에게 영향을 미치는 가장 흔한 신경학적 질환 중 하나입니다.
세계 장애 부담(GBD) 분석에서 두 번째로 높은 장애보정수명손실(YLD) 원인으로 꼽히며, 특히 20–50세 여성에서 최고 유병률을 보입니다.

| 항목 | 수치 |
|------|------|
| 전 세계 유병률 | ~15% (여성 20%, 남성 10%) |
| 만성 편두통(≥15 MHD) | 유병률의 약 8% |
| 평균 발작 지속시간 | 4–72시간 |
| 연간 직접 비용 (미국) | > $36억 달러 |
| DALY (global) | 45.1 million DALY (GBD 2016) |

---

## 주요 병태생리 경로 (Key Pathophysiological Pathways)

### 1. 삼차신경혈관 경로 (Trigeminovascular Pathway)

| 구성요소 | 역할 |
|----------|------|
| 삼차신경절 (TG) | CGRP/SP/PACAP 합성 및 방출 |
| 삼차신경핵 꼬리 (TNC) | 2차 뉴런, 중추감작 개시 |
| 삼차신경경부복합체 (TCC) | C1–C2 수렴, 목·뒷머리 통증 |
| 시상 (Thalamus) | 3차 뉴런, 통증 인식 |
| 대뇌 피질 | 광과민증·음과민증 처리 |

### 2. CGRP 신호 경로 (CGRP Signaling)

편두통 발작 시 삼차신경 말단에서 CGRP가 과다 방출되어:
- **CLR/RAMP1 수용체** 활성화 → cAMP↑ → KATP 채널 개방 → **경막 혈관 확장**
- 경막 비만세포 탈과립 → 히스타민·브래디키닌 방출 → 삼차신경 말초 감작
- CGRP 항체(에레누맙·프레마네주맙·갈카네주맙) 및 게판트(리메게판트·유브로게판트)의 치료 표적

### 3. 피질 확산 억제 (CSD — Cortical Spreading Depression)

CSD는 아우라(시각 섬광암점)의 신경생리학적 기전이며, 삼차신경 말초 감작을 통해 두통 발작을 유발합니다:
- K⁺ 유출(30–80 mM) → NMDA 수용체 활성화 → Ca²⁺ 유입 → 신경세포 탈분극
- 전파속도: 3–5 mm/분 (후두엽 → 전두엽)
- **CSD → Panx1 채널 활성화 → 경막 삼차신경 활성화 → CGRP 방출**

### 4. 중추 감작 (Central Sensitization)

발작 개시 10분 이상 경과 시 TNC 2차 뉴런의 과활성화:
- 피부 이질통(Cutaneous allodynia) — 접촉·온도에 대한 과민
- PKCε → MAPK/ERK → BDNF 상향조절 → 만성 편두통으로 진행
- **트립탄: 중추감작 전(10분 내)에 복용 시 가장 효과적**

---

## 약물 PK/PD 파라미터 (Drug PK/PD Parameters)

### 급성기 치료 (Acute Treatments)

| 약물 | 작용기전 | 주요 PK | 2h 통증소실률 |
|------|---------|---------|--------------|
| 수마트립탄 SC 6 mg | 5-HT1B/1D 작용제 | F=97%, t½=2h, CL=72L/h | ~35–40% |
| 수마트립탄 경구 100 mg | 5-HT1B/1D 작용제 | F=14%, t½=2h | ~25–30% |
| 라스미디탄 200 mg | 5-HT1F 작용제 (혈관수축 없음) | F=38%, t½=5h | ~32–39% |
| 리메게판트 75 mg | CGRP-R 길항제 (게판트) | F=64%, t½=11h | ~21% |
| 유브로게판트 100 mg | CGRP-R 길항제 (게판트) | F=44%, t½=7h | ~19% |

### 예방 치료 (Preventive Treatments)

| 약물 | 작용기전 | 투여 | MMD 감소 | 50% 반응률 |
|------|---------|------|----------|-----------|
| 에레누맙 140 mg | 항-CGRP 수용체 mAb | SC 월 1회 | 3.7일 | 47% |
| 프레마네주맙 225 mg | 항-CGRP 리간드 mAb | SC 월 1회 | 3.7일 | 43% |
| 갈카네주맙 120 mg | 항-CGRP 리간드 mAb | SC 월 1회 | 4.7일 | 52% |
| 리메게판트 75 mg QOD | CGRP-R 길항제 | 경구 격일 | 1.75일 | 28% |
| 토피라메이트 100 mg | AMPA/CA 억제 | 경구 분할 | 2.1일 | 37% |
| 프로프라놀롤 160 mg | β1/β2-차단제 | 경구 분할 | 1.8일 | 35% |
| 아미트립틸린 50 mg | SERT/NET 억제 (TCA) | 경구 야간 | 1.6일 | 30% |

---

## 모델 구조 (Model Structure)

### mrgsolve ODE 구획 (18개 상태변수)

```
약물 PK 구획 (10개):
  DEPOT_SUM, CENT_SUM, PERI_SUM     — 수마트립탄 SC 2구획
  DEPOT_ERE, CENT_ERE, PERI_ERE     — 에레누맙 SC 2구획
  CENT_RIM                           — 리메게판트 1구획
  DEPOT_TOP, CENT_TOP                — 토피라메이트 2구획
  CGRPR_FREE                         — 유리 CGRP 수용체

질환 PD 구획 (8개):
  CGRP_TG     — 삼차신경절/혈장 CGRP (pmol/L)
  CSD_ACT     — CSD 활성도 (0–1)
  TG_ACT      — 삼차신경 활성화 (0–1)
  CS_STATE    — 중추감작 상태 (0–1)
  PGE2_COMP   — 조직 PGE2 (pg/mL)
  NO_COMP     — 산화질소 수준 (pmol/L)
  SEROTONIN   — 혈소판 5-HT (ng/mL)
  PAIN_SCORE  — VAS 통증 점수 (0–10)
```

### 치료 시나리오 (7개)

1. **무치료 급성 발작** — CSD 유발 후 24h 통증 추적
2. **수마트립탄 SC 6 mg** — 5-HT1B/1D 작용 → CGRP 방출 억제 + 혈관 수축
3. **라스미디탄 200 mg** — 선택적 5-HT1F 작용 (CNS 침투 우수)
4. **리메게판트 75 mg** — 급성 + 예방 이중 효과 (CGRP-R 길항)
5. **에레누맙 140 mg 월 1회 × 3개월** — CGRP-R 점유율 및 MMD 감소
6. **토피라메이트 100 mg/day SS** — CSD 억제 + AMPA 차단
7. **만성 편두통 vs 에레누맙 1년** — 장기 진행 모델링

---

## 파일 목록 (File List)

| 파일 | 설명 |
|------|------|
| [mgr_qsp_model.dot](mgr_qsp_model.dot) | Graphviz 기계론적 지도 소스 (100+ 노드, 12 서브그래프 클러스터) |
| [mgr_qsp_model.svg](mgr_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [mgr_qsp_model.png](mgr_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [mgr_mrgsolve_model.R](mgr_mrgsolve_model.R) | mrgsolve ODE 모델 (18 구획, 7 시나리오) |
| [mgr_shiny_app.R](mgr_shiny_app.R) | Shiny 대시보드 (6탭) |
| [mgr_references.md](mgr_references.md) | 참고문헌 50편 (PubMed 링크 포함) |

---

## 기계론적 지도 미리보기 (Mechanistic Map Preview)

[![Migraine QSP Model](mgr_qsp_model.png)](mgr_qsp_model.svg)

*클릭하면 확대 가능한 SVG 이미지로 연결됩니다.*

---

## 주요 임상시험 벤치마크 (Key Clinical Trial Benchmarks)

| 임상시험 | 약물 | 1차 평가변수 | 결과 |
|----------|------|-------------|------|
| STRIVE | 에레누맙 140mg | MMD 감소 | -3.7일 (vs -1.8 위약) |
| EVOLVE-1 | 갈카네주맙 120mg | MMD 감소 | -4.7일 (vs -2.8 위약) |
| HALO-EM | 프레마네주맙 225mg | MMD 감소 | -3.7일 (vs -2.5 위약) |
| SAMURAI | 라스미디탄 200mg | 2h 통증소실 | 32.2% (vs 15.3% 위약) |
| ARTISAN-EM | 리메게판트 75mg | 2h 통증소실 | 21.2% (vs 10.9% 위약) |
| ACHIEVE-I | 유브로게판트 100mg | 2h 통증소실 | 19.2% (vs 11.8% 위약) |

---

## 생성일 (Generated)

2026-06-20 · QSP Disease Model Library (CCR)
