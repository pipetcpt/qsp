# 만성 췌장염 (Chronic Pancreatitis, CP) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도

[![CP QSP Model](cp_qsp_model.png)](cp_qsp_model.svg)

## 개요 (Overview)

만성 췌장염은 반복적인 염증 손상에 의해 췌장 실질이 비가역적으로 섬유화되는 진행성 질환으로, 서구권에서는 알코올, 아시아에서는 특발성이 주요 원인이다. 전 세계 유병률은 인구 10만 명당 약 50명으로 추정되며, 통증·외분비 부전(지방변증)·내분비 부전(당뇨병)이 3대 임상 결과이다. 췌성상세포(PSC) 활성화가 섬유화의 핵심 구동 인자이며, 효소 보충요법(PERT)으로 소화 불량을 교정하고 통증 관리에는 opioid·gabapentin을 사용한다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 세포 내 트립신 조기 활성화 | 리소좀 불안정 → 자가 소화 | 반복 급성 악화 |
| 췌성상세포(PSC) 활성화 | TGF-β, IL-6 → PSC 활성 → 콜라겐 과잉 | 췌장 실질 섬유화 |
| 신경 감작 | TNF-α, IL-6 → 말초/중추 신경 과민 | 만성 내장 통증 |
| 외분비 부전 | 선세포 소실, 주췌관 협착 | 지방변증, 영양 불량 |
| 내분비 부전 | 베타세포 섬유화 침윤 | 췌장인성 당뇨병(T3cDM) |
| 산화 스트레스 | ROS 축적 → 선세포 손상 | 세포사, 염증 지속 |
| 췌관 고압·결석 | 점액 마개·췌석 → 담관 폐쇄 | 통증 악화, 감염 위험 |

## 주요 약물 표적 (Drug Targets)

- **췌장 효소 보충(PERT — 판크레아제/리파제)**: 외분비 부전 교정 → 지방 흡수율 개선, 영양 상태 회복
- **아편유사제(Opioid)**: 중추·말초 mu 수용체 → 내장 통증 감소, 중독 위험 동반
- **Gabapentin/Pregabalin**: 전압 의존성 Ca²⁺ 채널(α2δ) 억제 → 중추 신경 감작 완화
- **항섬유화제(Pirfenidone/로사르탄)**: TGF-β 신호 억제 → PSC 활성·콜라겐 합성 감소 (연구용)
- **항산화 보충제**: 셀레늄·비타민 C/E·메티오닌 → ROS 감소 (EUROPAC-2 등)

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [cp_qsp_model.dot](cp_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 400 노드 / 12 클러스터) |
| [cp_qsp_model.svg](cp_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [cp_qsp_model.png](cp_qsp_model.png) | PNG 이미지 (150 dpi) |
| [cp_mrgsolve_model.R](cp_mrgsolve_model.R) | mrgsolve ODE 모델 (약 20 구획 / 5개 치료 시나리오) |
| [cp_shiny_app.R](cp_shiny_app.R) | Shiny 대시보드 |
| [cp_references.md](cp_references.md) | 참고문헌 (약 55편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: Opioid (gut·중심·말초 2구획) + PERT 십이지장 구획 PK + TNF-α, IL-6, TGF-β, ROS, PSC, 섬유화 지수(FIB), 외분비 기능(EXO), 베타세포(BETA), 혈당(GLUC), 말초 신경 감작(pSENS), 중추 신경 감작(cSENS), 통증 점수 PD 구획
- **주요 치료 시나리오**: (1) 무치료 자연경과, (2) PERT 단독, (3) Opioid + Gabapentin 병용, (4) 항섬유화 치료(Pirfenidone·로사르탄), (5) PERT + 진통 + 항섬유화 복합 치료
- **보정/근거**: PERT 효과는 Lévy P et al. (Gut 2012), 통증 관리 시나리오는 USPS 패널 권고(Pancreas 2020) 기반

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) 환자 프로파일 — 원인·중증도·췌장 잔여 기능 설정; (2) PK 탭 — Opioid 혈장 농도·PERT 십이지장 농도; (3) PD 주요지표 — 섬유화 지수, TNF-α, TGF-β, ROS; (4) 임상 엔드포인트 — 통증 점수, 지방 흡수율, 혈당; (5) 시나리오 비교 — 5가지 치료 전략 2년 결과; (6) 바이오마커 — 외분비 기능(분변 탄력효소), 베타세포 기능, 혈청 아밀라제/리파제

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("cp_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("cp_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg cp_qsp_model.dot -o cp_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [cp_references.md](cp_references.md) 참조 (약 55편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
