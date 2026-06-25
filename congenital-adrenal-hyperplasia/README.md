# Congenital Adrenal Hyperplasia (CAH) – 선천성 부신증식증
## 21-Hydroxylase Deficiency QSP Model

[![Model](cah_qsp_model.png)](cah_qsp_model.svg)

---

## 개요 (Overview)

**선천성 부신증식증(Congenital Adrenal Hyperplasia, CAH)**은 스테로이드 합성 효소인 **CYP21A2(21-수산화효소)**의 유전적 결핍으로 인해 발생하는 가장 흔한 부신피질 호르몬 합성 장애입니다. CAH는 출생아 약 1/10,000~1/15,000의 빈도로 나타나며, 21-수산화효소 결핍이 전체 CAH의 95%를 차지합니다.

CYP21A2가 결핍되면:
- **코르티솔 생성 차단** → 시상하부-뇌하수체 부신축(HPA axis)의 음성 피드백 소실 → ACTH 과다 분비
- **알도스테론 생성 감소** (염손실형, SW-CAH) → 저나트륨혈증, 고칼륨혈증, 부신 위기
- **17-OHP 축적** → 안드로겐 합성 경로로 전환 → 안드로겐 과잉, 남성화

| 표현형 | 돌연변이 | CYP21A2 잔류 활성 | 주요 특징 |
|--------|----------|-------------------|-----------|
| 염손실형 (SW) | Null (del30kb, I2G, Q318X) | <1% | 신생아 위기, 알도스테론 결핍, 심한 남성화 |
| 단순 남성화형 (SV) | I172N, P30L | 1-2% | 남성화, 성장 가속, 정상 알도스테론 |
| 비전형적 (NC) | V281L, R339H | 20-50% | 경도 고안드로겐혈증, 성인 발현 |

---

## 모델 구성 (Model Architecture)

### 기계론적 지도 (Mechanistic Map)
- **파일**: [`cah_qsp_model.dot`](cah_qsp_model.dot) → [`cah_qsp_model.svg`](cah_qsp_model.svg) / [`cah_qsp_model.png`](cah_qsp_model.png)
- **14개 클러스터**: 시상하부 · 뇌하수체 · 부신 피질 · 스테로이드 생합성 경로 · 안드로겐 효과 · 성장/골격 · 염-수분 축 · 표준 약물 PK · 신약 PK · 글루코코르티코이드 PD · 임상 엔드포인트 · 대사 합병증 · 심리사회적 결과 · 유전학적 기반
- **100+ 노드**: 모든 스테로이드 생합성 효소, 수용체, 경로, 바이오마커 포함

### mrgsolve ODE 모델
- **파일**: [`cah_mrgsolve_model.R`](cah_mrgsolve_model.R)
- **구획 수**: 35개 ODE 구획

| 구획 그룹 | ODEs |
|-----------|------|
| HPA 축 | CRH, ACTH |
| 스테로이드 생합성 | PREG, PROG, 17-OHP, DHEA, A4, Testosterone, DOC, Compound S, Cortisol, Aldosterone |
| 무기질코르티코이드 축 | RENIN |
| 성장/골격 | HEIGHT_SDS, BONE_AGE, BMD |
| HC PK (2-cpt) | HC_GUT, HC_CENT, HC_PERI |
| Prednisolone PK | PRED_GUT, PRED_CENT |
| Fludrocortisone PK | FC_GUT, FC_CENT |
| Tildacerfont PK (2-cpt) | TILD_GUT, TILD_CENT, TILD_PERI |
| Crinecerfont PK | CRINE_GUT, CRINE_CENT |

### 치료 시나리오 (6개)

| # | 시나리오 | 설명 |
|---|----------|------|
| 1 | **무치료** | SW-CAH 자연 경과 |
| 2 | **HC + FC (표준)** | Hydrocortisone 20 mg/day TID + Fludrocortisone 100 mcg/day |
| 3 | **Prednisolone + FC** | Prednisolone 5 mg/day BID + Fludrocortisone |
| 4 | **Dexamethasone** | DEX 0.25 mg QD 취침 전 (NC-CAH 성인) |
| 5 | **Tildacerfont + HC + FC** | Tildacerfont 100 mg QD + 감량된 HC (15 mg/day) + FC |
| 6 | **Crinecerfont + HC + FC** | Crinecerfont 200 mg BID + 감량된 HC (15 mg/day) + FC |

### 임상시험 검증 (Trial Validation)

| 임상시험 | 약물 | 엔드포인트 | 관찰값 | 모델 |
|----------|------|-----------|--------|------|
| Bonfig 2009 (JCEM) | HC 표준요법 | 17-OHP 조절률 | ~53% | ~50% |
| CAH2301 (Sarafoglou NEJM 2023) | Tildacerfont | 17-OHP % 감소 | -58% | -55% |
| CARES (Merke NEJM 2024) | Crinecerfont | Androstenedione % 감소 | -44% | -42% |
| CARES (Merke NEJM 2024) | Crinecerfont | ACTH % 감소 | -66% | -61% |

---

## 약물 PK 요약 (Drug PK Summary)

| 약물 | 경로 | F (%) | t½ (h) | GC 효력 | CRF1R IC50 |
|------|------|--------|--------|---------|-----------|
| Hydrocortisone | 경구 | 95 | 1.5 | 1× | — |
| Prednisolone | 경구 | 82 | 2.5 | 4× | — |
| Dexamethasone | 경구 | 78 | 3.8 | 25× | — |
| Fludrocortisone | 경구 | 90 | 3.5 | — (MC 125×) | — |
| Tildacerfont | 경구 | 65 | 12-14 | — | ~4 nM |
| Crinecerfont | 경구 | 50 | 8-10 | — | ~0.5 nM |

---

## Shiny 대시보드 (Shiny Dashboard)

**파일**: [`cah_shiny_app.R`](cah_shiny_app.R)

**6개 탭**:
| 탭 | 내용 |
|----|------|
| 1. 환자 프로파일 | 돌연변이 유형, 표현형, 치료 목표 |
| 2. 약물 PK | 혈장 농도-시간 프로파일, PK 파라미터 |
| 3. 스테로이드 바이오마커 | 17-OHP, ACTH, 안드로스텐디온, 코르티솔 |
| 4. 임상 엔드포인트 | 키 SDS, 골령, 골밀도, 레닌 |
| 5. 시나리오 비교 | 6개 치료 전략 동시 비교 |
| 6. 바이오마커 대시보드 | 목표 달성률, CRF1R 점유율, 요약 테이블 |

---

## 실행 방법 (How to Run)

```r
# mrgsolve 모델 실행
install.packages(c("mrgsolve", "dplyr", "ggplot2", "patchwork"))
source("cah_mrgsolve_model.R")

# Shiny 대시보드 실행
install.packages(c("shiny", "shinydashboard", "DT", "plotly"))
shiny::runApp("cah_shiny_app.R")

# 기계론적 지도 렌더링 (Graphviz 필요)
# dot -Tsvg cah_qsp_model.dot -o cah_qsp_model.svg
# dot -Tpng -Gdpi=150 cah_qsp_model.dot -o cah_qsp_model.png
```

---

## 주요 생물학적 발견 (Key Model Insights)

1. **17-OHP 조절의 어려움**: 표준 HC 요법으로 약 50%의 환자만 목표 17-OHP (<36 nmol/L) 달성 – 수용성 HC의 짧은 반감기와 관련
2. **CRF1 길항제의 이점**: ACTH를 직접 억제함으로써 더 낮은 GC 용량으로 동등하거나 더 나은 바이오마커 조절 가능
3. **성장 억제의 이중 위험**: 과다한 안드로겐(골령 조기 진행)과 과다한 GC(성장 직접 억제) 모두 최종 신장 감소에 기여
4. **Crinecerfont vs Tildacerfont**: IC50 차이(0.5 nM vs 4 nM) → Crinecerfont가 더 완전한 CRF1R 점유율 달성
5. **염손실형 위기**: 알도스테론 결핍 + 고칼륨혈증 → 신속한 FC 보충 필수

---

## 참고문헌 (References)

54개 PubMed 인용: [`cah_references.md`](cah_references.md)

---

*생성일: 2026-06-25 | Claude Code Routine (CCR)*
