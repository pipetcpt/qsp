# 만성 폐쇄성 폐질환 (Chronic Obstructive Pulmonary Disease, COPD) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 호흡기

[![COPD QSP Model](copd_qsp.png)](copd_qsp.svg)

## 개요 (Overview)

COPD는 전 세계적으로 약 3억 9,000만 명이 이환된 주요 만성 호흡기 질환으로, 흡연이 가장 흔한 원인이다. 흡연에 의한 산화 스트레스와 단백분해효소(neutrophil elastase·MMP) 과활성이 기도 상피 손상과 폐포 파괴(기종)를 유발하며, IL-8·TNF-α를 매개로 한 지속적 기도 염증이 기도 협착을 진행시킨다. GOLD 가이드라인에 따라 LAMA·LABA 기관지확장제가 기본 치료이고, 호산구형 또는 잦은 악화형에는 ICS 추가, 중증 폐기종형에는 PDE4 억제제(로플루밀라스트)가 병용된다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 흡연·산화 스트레스 | ROS → NF-κB → IL-8·TNF-α 분비 | 호중구·대식세포 기도 침윤 |
| 단백분해효소/항단백분해효소 불균형 | NE·MMP-9 ↑ / α1-AT·TIMP ↓ | 폐포 격벽 파괴 → 기종 |
| 기도 점액 과분비 | 배상세포 증식, MUC5AC 과발현 | 만성 기침·객담 |
| 기도 리모델링 | 섬유아세포 활성화, 기도벽 비후 | 기류 제한(FEV₁ 감소) |
| 호산구형 기도 염증 | IL-5·IL-13 → 호산구 집적 | ICS 반응 예측, 악화 위험 |
| 폐동맥 혈관 리모델링 | 저산소증 → 폐혈관 수축 → 폐동맥 고혈압 | 우심부전(폐성심) |
| 급성 악화(AECOPD) | 바이러스/세균 감염 + 기저 염증 | 폐기능 급속 악화, 입원 |

## 주요 약물 표적 (Drug Targets)

- **LAMA(장시간 지속 무스카린 길항제)**: M3 수용체 차단 → 기관지 확장 → FEV₁ 개선 (티오트로피움·글리코피로니움)
- **LABA(장시간 지속 베타₂ 작용제)**: β₂ 수용체 자극 → 평활근 이완 (살메테롤·올로다테롤)
- **ICS(흡입 코르티코스테로이드)**: IL-8·호산구 염증 억제, 악화 예방 (플루티카손·부데소니드)
- **PDE4 억제제(로플루밀라스트)**: cAMP 분해 억제 → 항염증, FEV₁ 개선, 급성 악화 감소
- **항생제(아지트로마이신)**: 잦은 악화형에 장기 예방적 투여, 항균 + 항염증

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [copd_qsp.dot](copd_qsp.dot) | Graphviz 기계론적 지도 소스 (약 430 노드 / 12 클러스터) |
| [copd_qsp.svg](copd_qsp.svg) | SVG 벡터 이미지 (확대 가능) |
| [copd_qsp.png](copd_qsp.png) | PNG 이미지 (150 dpi) |
| [copd_mrgsolve_model.R](copd_mrgsolve_model.R) | mrgsolve ODE 모델 (약 26 구획 / 6개 치료 시나리오) |
| [copd_shiny_app.R](copd_shiny_app.R) | Shiny 대시보드 |
| [copd_references.md](copd_references.md) | 참고문헌 (약 52편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: LAMA·LABA·ICS (폐·중심·말초 2구획) + PDE4i 혈장 PK 구획 + IL-8, 호중구 탄력효소(NE), CRP, 호산구, FEV₁, 폐기종 지수(Emph), 폐혈관 저항(PVR), 급성 악화 누적·연율, CAT 점수 PD 구획
- **주요 치료 시나리오**: (1) 무치료, (2) LAMA 단독, (3) LAMA + LABA 병용, (4) LAMA + LABA + ICS 3제, (5) LAMA + LABA + ICS + PDE4i 4제, (6) 흡연 중단 단독
- **보정/근거**: FEV₁ 변화는 UPLIFT(티오트로피움), TRILOGY/TRINITY(ICS-LABA-LAMA) 임상시험 데이터 기반

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) 환자 프로파일 — GOLD 단계·흡연력·호산구 수 설정; (2) PK 탭 — LAMA/LABA/ICS/PDE4i 폐·혈장 농도; (3) PD 주요지표 — IL-8, 호중구 탄력효소, CRP, 호산구; (4) 임상 엔드포인트 — FEV₁ 추이, 급성 악화 누적, CAT 점수; (5) 시나리오 비교 — 6가지 치료 전략 1년 결과; (6) 바이오마커 — 혈청 피브리노겐, 혈중 호산구, 폐기종 진행 지수

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("copd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("copd_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg copd_qsp.dot -o copd_qsp.svg
```

## 참고문헌 (References)

자세한 인용은 [copd_references.md](copd_references.md) 참조 (약 52편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
