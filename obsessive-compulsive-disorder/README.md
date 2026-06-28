# OCD — Obsessive-Compulsive Disorder QSP Model

> **강박장애 (강박증)** 에 대한 정량적 시스템 약리학(QSP) 모델: CSTC 회로 역학, SERT 점유 기반 SSRI/클로미프라민 PK/PD, ERP 치료 효과, Y-BOCS 임상 엔드포인트를 통합한 시뮬레이션

---

## 기계론적 지도 (Mechanistic Map)

[![OCD QSP Mechanistic Map](ocd_qsp_model.png)](ocd_qsp_model.svg)

*클릭하면 확대 가능한 SVG 버전을 볼 수 있습니다.*

---

## 개요 (Overview)

**강박장애(OCD)**는 전 세계 인구의 약 2–3%에서 발생하는 만성 신경정신 질환으로, 반복적이고 원치 않는 침습적 생각(강박 사고)과 이를 완화하기 위한 반복 행동(강박 행동)이 특징입니다.

### 핵심 병태생리

| 구성 요소 | 역할 |
|-----------|------|
| **OFC (안와전두피질) 과활성** | 침습적 사고, 오류 신호 생성 |
| **미상핵 과대사** | 강박 행동 충동 유지 |
| **직접 경로(Direct/Go) 우세** | 시상 탈억제 → 회로 과부하 |
| **간접 경로(Indirect/No-Go) 약화** | 행동 억제 실패 |
| **SERT 기능 저하** | 시냅스 5-HT 감소 |
| **OFC → 미상핵 글루타메이트 과활성** | 피질-선조체 구동력 증가 |

---

## 파일 구성 (File Structure)

| 파일 | 내용 |
|------|------|
| `ocd_qsp_model.dot` | Graphviz DOT 기계론적 지도 (130+ 노드, 13개 클러스터) |
| `ocd_qsp_model.svg` | 벡터 그래픽 (확대 가능) |
| `ocd_qsp_model.png` | 래스터 이미지 (150 dpi) |
| `ocd_mrgsolve_model.R` | mrgsolve ODE PK/PD 모델 (6개 치료 시나리오, 가상 환자 100명) |
| `ocd_shiny_app.R` | 인터랙티브 Shiny 대시보드 (8개 탭) |
| `ocd_references.md` | 참고문헌 50편 (PubMed 링크 포함) |
| `README.md` | 이 파일 |

---

## 모델 구성 (Model Architecture)

### 기계론적 지도 — 13개 클러스터

```
1. SSRI PK           ← sertraline 2구획 + BBB 수송 + CYP2D6/3A4 대사
2. Clomipramine PK   ← 초회 통과 효과 + 활성 대사체(DCMI)
3. 증강 요법 PK      ← risperidone, aripiprazole, memantine, D-cycloserine
4. 세로토닌 시스템   ← TPH2/AADC 합성 → SERT → MAO-A → 5-HT1A/2A/2C 수용체
5. 도파민 시스템     ← VTA/SNc → 선조체/OFC → D1R(직접) / D2R(간접)
6. 글루타메이트/GABA ← OFC-선조체 글루타메이트 시냅스 → NMDA/AMPA/mGluR
7. CSTC 회로         ← OFC → 미상핵 → GPi/GPe/STN → 시상 → OFC 루프
8. HPA 축/스트레스   ← CRH → ACTH → 코르티솔 → 해마 위축 → 불안 루프
9. 신경염증          ← 미세아교세포 → TNF-α/IL-6 → IDO → 키누레닌/QUIN
10. 신경가소성        ← BDNF → TrkB → ERK/AKT → 수상돌기 리모델링 → LTP/LTD
11. CBT/ERP          ← 공포소거 → OFC 정상화 → 습관 회로 재형성
12. 임상 엔드포인트  ← Y-BOCS 총점/소척도 → 반응/관해 기준
13. 유전/바이오마커  ← SLC6A4/HTR2A/COMT/SLC1A1 변이 → 영상/PET 바이오마커
```

### mrgsolve ODE 모델 — 23개 구획

```
약동학(PK):
  AG_SSRI → A1_SSRI → A2_SSRI / A_CNS     (sertraline 2구획 + CNS)
  AG_CMI  → A1_CMI  → A_DCMI               (clomipramine + 활성 대사체)
  AG_RISP → A1_RISP                         (risperidone)

약력학(PD):
  SERT_OCC   ← SSRI/CMI CNS 농도 Emax 모델
  HT5_SYN    ← 시냅스 5-HT (방출-재흡수-분해 역학)
  DES_5HT1   ← 5-HT1A 자가수용체 탈감작 (지연 반응 기전)
  OFC_ACT    ← OFC 활동도 (CSTC 루프 + 5-HT 조절)
  CAUD_ACT   ← 미상핵 활동도
  THAL_ACT   ← 시상 활동도
  DIR_PATH   ← 직접 경로 활성화
  IND_PATH   ← 간접 경로 활성화
  D2R_OCC    ← 항정신병약물 D2R 점유율
  BDNF_LV    ← BDNF 수준 (신경가소성 마커)
  ERP_EFF    ← ERP 누적 효과 (0–1)
  YBOCS      ← Y-BOCS 점수 (0–40, 임상 엔드포인트)
  ANXIETY    ← 불안 상태 (정규화)
```

---

## 치료 시나리오 (6개)

| 시나리오 | 치료법 | 특징 |
|----------|--------|------|
| 1 | 무처치 (기저선) | Y-BOCS = 28 지속 |
| 2 | **Sertraline 200 mg/일** | SERT ≥80% 목표, 6–12주 지연 효과 |
| 3 | **Clomipramine 250 mg/일** | 가장 강력한 SERT 억제, 부작용 ↑ |
| 4 | Sertraline + **Risperidone 1.5 mg** | 12주 SSRI 불충분 반응 후 D2 증강 |
| 5 | Sertraline + **ERP** (병용) | 최고 효능, OFC 정상화 기전 2중 |
| 6 | **ERP 단독** | 심리치료 단독, 관해율 ↓ |

---

## 핵심 파라미터 (Key Parameters)

| 파라미터 | 값 | 근거 |
|----------|----|----|
| SERT EC50 (sertraline) | 1.2 ng/mL | Zitterl et al. 2008 |
| SERT EC50 (clomipramine) | 0.4 ng/mL | 추정값 (CMI 더 강력) |
| SERT 점유율 필요 임계값 | ≥80% | OCD vs 우울증(60%) |
| 5-HT1A 탈감작 t½ | ~4주 | SSRI 지연 반응 설명 |
| Y-BOCS 반응 기준 | ≥35% 감소 | Goodman et al. 1989 |
| 관해 기준 | Y-BOCS ≤12 | APA 가이드라인 |
| SSRI 반응률 | 40–60% | Soomro et al. 2008 |
| ERP 최대 OFC 억제 | 45% | Foa et al. 2005 기반 |

---

## 실행 방법 (How to Run)

### mrgsolve 모델

```r
library(mrgsolve)
library(dplyr)
library(ggplot2)

# 모델 실행 (6개 치료 시나리오)
source("ocd_mrgsolve_model.R")
```

### Shiny 대시보드

```r
library(shiny)
library(shinydashboard)
library(mrgsolve)
library(plotly)
library(DT)

shiny::runApp("ocd_shiny_app.R")
```

### Graphviz 렌더링

```bash
dot -Tsvg ocd_qsp_model.dot -o ocd_qsp_model.svg
dot -Tpng -Gdpi=150 ocd_qsp_model.dot -o ocd_qsp_model.png
```

---

## Shiny 대시보드 탭 (8개)

| 탭 | 내용 |
|----|------|
| 🧑 **Patient Profile** | 환자 특성, 약물유전체(CYP2D6/SLC6A4/COMT), 중증도 |
| 💊 **PK — Drug Levels** | SSRI/클로미프라민 혈중·뇌 내 농도, SERT 점유율 곡선 |
| 🧠 **PD — Neurotransmitters** | 시냅스 5-HT 역학, 5-HT1A 탈감작, SERT-5HT 상관 |
| 🔄 **CSTC Circuit** | OFC/미상핵/시상 활동도, 직접/간접 경로 균형, BDNF |
| 📊 **Clinical Endpoints** | Y-BOCS 궤적, 반응/관해 valueBox, 불안 상태 |
| ⚖️ **Scenario Comparison** | 6개 시나리오 병렬 비교, 결과 표 |
| 🔬 **Biomarkers** | 영상 바이오마커, 약물유전체 영향, DBS 시뮬레이션 |
| ℹ️ **About** | 모델 설명, 한계, 핵심 참고문헌 |

---

## 임상적 의의 (Clinical Significance)

- **SSRI 고용량 필요성**: 우울증(SERT 60%)과 달리 OCD는 SERT ≥80% 점유 필요 → 최대 허용 용량까지 증량
- **지연 반응 기전**: 5-HT1A 자가수용체 탈감작이 2–4주 걸려 임상 효과 지연 → 조기 중단 주의
- **ERP + SSRI 병용**: 단독 치료보다 우월, 각각 다른 기전(약물: 5-HT ↑ / ERP: OFC 하향 조절)
- **증강 요법**: SSRI 단독 불충분 시 리스페리돈(D2 차단)으로 간접 경로 복원

---

## 참고문헌 (Key References)

- Soomro GM et al. (2008) Cochrane Review: SSRIs for OCD
- Foa EB et al. (2005) JAMA: ERP vs Clomipramine RCT
- Zitterl W et al. (2008) Neuropsychopharmacology: SERT occupancy
- Goodman WK et al. (1989) Arch Gen Psychiatry: Y-BOCS development
- Saxena S & Rauch SL (2004) Annu Rev Neurosci: CSTC circuit
- Bloch MH et al. (2006) Mol Psychiatry: augmentation meta-analysis

→ 전체 50편: [`ocd_references.md`](ocd_references.md)

---

*Generated by Claude Code Routine (CCR) — 2026-06-28*
