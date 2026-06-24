# 악성 빈혈 (Pernicious Anemia, PNA) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 혈액

[![PNA QSP Model](pna_qsp_model.png)](pna_qsp_model.svg)

## 개요 (Overview)
악성 빈혈은 자가면역 위염에 의한 내인자(intrinsic factor, IF) 결핍으로 비타민 B12(코발라민) 흡수가 차단되어 발생하는 거대적아구성 빈혈이다. 일반 인구 약 0.1%에서 유병하며, 60세 이상 노인과 여성에서 더 흔하고 HLA-DR3/DR4 소인을 가진다. 항내인자 항체(IF 차단항체 70%)와 항위벽세포 항체(APC-Ab 90%)가 핵심 발병기전으로, 위체부 벽세포 파괴 → IF 결핍 → B12 흡수 불가 → 조혈 장애 및 신경계 합병증(아급성 연합 척수 변성)으로 이어진다. 비타민 B12 근육주사 또는 고용량 경구 보충이 치료의 근간이다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 자가면역 위염 경로 | APC-Ab·IF-Ab, CD4 Th1 자가반응성 T세포 | 벽세포 파괴·위산 무산증 |
| 내인자 결핍 경로 | IF 분비 소실 → cobalamin-IF 복합체 형성 불가 | 회장 말단 B12 흡수 차단 |
| B12 전신 분포 경로 | 트랜스코발라민 II(TC-II) 매개 운반, 장간 순환 | 혈청 B12·holotranscobalamin 감소 |
| 거대적아구성 빈혈 경로 | DNA 합성 장애(thymidylate 경로), 세포 분열 이상 | MCV 증가·Hb 감소·망상적혈구 반응 |
| 신경독성 경로 | 메틸말론산·호모시스테인 축적, 미엘린 손상 | 아급성 연합 척수 변성·말초신경병증 |
| 경구 수동 흡수 경로 | IF 독립적 수동 확산 (~1%/용량) | 고용량 경구 요법 효과 |

## 주요 약물 표적 (Drug Targets)
- **비타민 B12 근육주사 (IM)**: 시아노코발라민/하이드록소코발라민 — 급속 보충, 신경 회복 우선 시 선호
- **고용량 경구 비타민 B12**: 1,000 µg/일 — 수동 흡수(1%) 이용, 순응도 우수
- **설하·비강 제제**: 대안 투여 경로 — IF 결핍 우회
- **면역억제 (실험적)**: 코르티코스테로이드 단기 투여 — 항체 억제, 임상 사용 제한적
- **엽산 보충**: 동반 결핍 교정 시 사용 (단독 투여 시 신경합병증 가릴 수 있어 주의)

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [pna_qsp_model.dot](pna_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 154 노드 / 11 클러스터) |
| [pna_qsp_model.svg](pna_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [pna_qsp_model.png](pna_qsp_model.png) | PNG 이미지 (150 dpi) |
| [pna_mrgsolve_model.R](pna_mrgsolve_model.R) | mrgsolve ODE 모델 (약 18 구획 / 5 치료 시나리오) |
| [pna_shiny_app.R](pna_shiny_app.R) | Shiny 대시보드 |
| [pna_references.md](pna_references.md) | 참고문헌 (약 35편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: IM 데포·경구 위장관 흡수·IF pool·장간문맥 B12·혈장·holotranscobalamin·간·골수·신경 조직 구획, 자가항체·위벽세포·Hb·MCV·망상적혈구·신경 합병증 지수 PD 구획
- **주요 치료 시나리오**: ① 무치료, ② IM 표준(1,000 µg/일 1주 후 월 1회), ③ IM 유지요법, ④ 경구 고용량 요법, ⑤ 집중 공격적 보충 요법
- **보정/근거**: 혈청 B12·MCV 정상화 시간 및 신경합병증 회복 경과를 문헌(Stabler 2013, Carmel 2008)에서 참고 보정

## Shiny 대시보드 (Dashboard)
환자 프로파일(항체 양성, 위염 단계, 신경증상 유무) · B12 약동학 PK 프로파일 · 혈액학 지표(Hb·MCV·망상적혈구) PD · 신경합병증 진행도 · 투여 경로 시나리오 비교 · 바이오마커(holotranscobalamin, 메틸말론산, 호모시스테인) 탭으로 구성

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("pna_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("pna_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg pna_qsp_model.dot -o pna_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [pna_references.md](pna_references.md) 참조 (약 35편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
