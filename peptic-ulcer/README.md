# 소화성 궤양 (Peptic Ulcer Disease, PUD) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도

[![PUD QSP Model](pud_qsp_model.png)](pud_qsp_model.svg)

## 개요 (Overview)
소화성 궤양은 위·십이지장 점막에 발생하는 만성 궤양으로, 전 세계 유병률은 약 10%에 달하며 매년 상당한 입원 및 출혈 합병증을 유발한다. 핵심 발병기전은 *Helicobacter pylori* 감염 및 NSAIDs 사용에 의한 점막 방어-공격 불균형으로, H. pylori의 CagA·VacA 독소가 점막 손상을 유발하고 NSAIDs는 COX-1 억제를 통해 프로스타글란딘 생성을 감소시킨다. 주요 치료 표적은 위산 분비 펌프(H⁺/K⁺-ATPase), H. pylori 균주, 및 점막 보호 기전이다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| H. pylori 독소 경로 | CagA T4SS, VacA, 산화 스트레스 | 점막 상피 손상·염증 |
| 위산 분비 경로 | H⁺/K⁺-ATPase (프로톤 펌프), 가스트린·히스타민 | 위내 산도 상승 |
| 점막 방어 경로 | 뮤신(mucus), 프로스타글란딘 E2, 중탄산염 | 점막 보호막 유지 |
| NSAIDs 경로 | COX-1 억제 → PGE2 감소 | 점막 혈류 감소·궤양 유발 |
| 염증 증폭 경로 | IL-1β, TNF-α, NF-κB 활성화 | 호중구 침윤·조직 손상 |
| 궤양 치유 경로 | EGF, TGF-β, 세포 재생 | 궤양 면적 감소 |

## 주요 약물 표적 (Drug Targets)
- **프로톤 펌프 억제제 (PPI)**: 오메프라졸·에소메프라졸 — H⁺/K⁺-ATPase 비가역적 억제, 위산 분비 차단
- **H2 수용체 길항제 (H2RA)**: 라니티딘·파모티딘 — 히스타민 H2 수용체 차단, 야간 산 분비 억제
- **H. pylori 제균 항생제**: 아목시실린(AMX)·클래리스로마이신(CLR) — 삼제·사제 요법
- **사이토프로텍티브**: 미소프로스톨 — PGE1 유사체로 점막 보호
- **항균 병용**: 메트로니다졸·비스무스 — 내성 균주 대응 사제 요법

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [pud_qsp_model.dot](pud_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 174 노드 / 9 클러스터) |
| [pud_qsp_model.svg](pud_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [pud_qsp_model.png](pud_qsp_model.png) | PNG 이미지 (150 dpi) |
| [pud_mrgsolve_model.R](pud_mrgsolve_model.R) | mrgsolve ODE 모델 (약 18 구획 / 5 치료 시나리오) |
| [pud_shiny_app.R](pud_shiny_app.R) | Shiny 대시보드 |
| [pud_references.md](pud_references.md) | 참고문헌 (약 55편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: PPI/H2RA/AMX/CLR/NSAID 각 약물의 위장관·중심·말초 PK 구획, H. pylori 부하(log10), 프로톤 펌프 활성 분율, 위내 pH, 점막층(뮤신), 프로스타글란딘, 염증 점수, 궤양 면적
- **주요 치료 시나리오**: ① PPI BID 단독, ② H2RA BID 단독, ③ 삼제요법(PPI+AMX+CLR), ④ NSAID 단독 투여, ⑤ NSAID + PPI 병용
- **보정/근거**: MAASTRICHT V 가이드라인 및 주요 H. pylori 제균 임상시험(MLST 데이터)에서 제균율·궤양 치유율을 참고하여 파라미터 보정

## Shiny 대시보드 (Dashboard)
환자 프로파일(H. pylori 상태, NSAID 사용, 위험인자) · 약물 PK 시뮬레이션 · 위내 pH 및 산 억제 정도 · 궤양 면적 변화 및 치유율 · 치료 시나리오 비교(제균율·재발률) · 주요 바이오마커(CRP, PG, HP 부하) 탭으로 구성된 인터랙티브 대시보드

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("pud_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("pud_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg pud_qsp_model.dot -o pud_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [pud_references.md](pud_references.md) 참조 (약 55편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
