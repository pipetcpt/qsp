# 만성 위염 (Chronic Gastritis, CGAST) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도

[![CGAST QSP Model](cgast_qsp_model.png)](cgast_qsp_model.svg)

## 개요 (Overview)

만성 위염은 *Helicobacter pylori* 감염이 가장 흔한 원인으로, 전 세계 인구의 약 44%가 감염된 것으로 추정된다. Correa 연쇄(Cascade)에 따라 만성 위염 → 위축성 위염 → 장상피화생 → 이형성 → 위암 순으로 진행하며, NF-κB 경로를 통한 IL-8·IL-1β·TNF-α 등의 사이토카인 폭포가 점막 손상을 주도한다. PPI 기반 제균 3제·4제요법은 H. pylori를 박멸하여 점막 염증을 해소하고 Correa 연쇄를 차단하는 핵심 치료 전략이다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| H. pylori 감염·NF-κB 활성화 | CagA/VacA 독소 → NF-κB → IL-8 분비 | 호중구 침윤, 점막 손상 |
| Th1/Treg 불균형 | IFN-γ ↑, IL-10 저하 | 만성 염증 지속, 위선 손실 |
| 위산 조절 이상 | 가스트린 ↑ (체부 위염), 위산 분비 ↑ | 소화성 궤양 위험 |
| 산화 스트레스 | ROS 생성, 점막 방어막(뮤신) 손상 | 세포 손상, 위축 촉진 |
| Correa 연쇄 — 위축 | 벽세포·주세포 소실, PG I/II 비 저하 | 위산 분비 감소, 세균 과증식 |
| 장상피화생 진행 | CDX2 발현, 장형 상피 치환 | 악성 전환 위험 |
| 제균 후 회복 | 염증 해소, 위축 부분 역전(조기 단계) | 암 위험 감소 |

## 주요 약물 표적 (Drug Targets)

- **PPI (프로톤 펌프 억제제)**: 오메프라졸·에소메프라졸 — H⁺/K⁺-ATPase 억제, 위산 분비 차단 → 항생제 효능 보조
- **Vonoprazan(P-CAB)**: 칼륨 경쟁적 위산 차단제, PPI보다 빠르고 강력한 산 억제
- **아목시실린**: 세포벽 합성 억제, H. pylori 1차 살균
- **클래리트로마이신**: 리보솜 50S 소단위 억제, 제균 3제요법의 핵심
- **메트로니다졸**: 항혐기성 항균, 클래리트로마이신 내성 대체제
- **비스무트(BSS)**: 점막 보호, 항 H. pylori 작용 → 4제요법 구성

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [cgast_qsp_model.dot](cgast_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 439 노드 / 10 클러스터) |
| [cgast_qsp_model.svg](cgast_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [cgast_qsp_model.png](cgast_qsp_model.png) | PNG 이미지 (150 dpi) |
| [cgast_mrgsolve_model.R](cgast_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 7개 치료 시나리오) |
| [cgast_shiny_app.R](cgast_shiny_app.R) | Shiny 대시보드 |
| [cgast_references.md](cgast_references.md) | 참고문헌 (약 61편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: PPI·아목시실린·클래리트로마이신·메트로니다졸 PK 구획 + H. pylori 균량, NF-κB, IL-8, IL-1β, TNF-α, IFN-γ, IL-10, 호중구, Th1, Treg, 가스트린, 위산, 점막, 위축 점수, 장상피화생 점수, 증상 PD 구획
- **주요 치료 시나리오**: (1) 무치료 자연경과, (2) PPI 단독, (3) 표준 3제요법(PPI+AMX+CLR), (4) 비스무트 4제요법, (5) 메트로니다졸 4제요법, (6) Vonoprazan 3제요법, (7) 제균 후 5년 추적
- **보정/근거**: 각 시나리오의 H. pylori 제균율은 Malfertheiner et al. (Gut 2017) 등 주요 가이드라인 데이터에 기반

## Shiny 대시보드 (Dashboard)

6개 이상의 탭으로 구성: (1) 환자 프로파일 — 감염 단계·항생제 내성 설정; (2) PK 탭 — PPI·항생제 혈장 농도; (3) PD 주요지표 — H. pylori 균량, 사이토카인(IL-8·IL-1β·TNF-α), 가스트린·위산; (4) 임상 엔드포인트 — 위축/장상피화생 점수, 증상 점수; (5) 시나리오 비교 — 7가지 제균 요법 결과; (6) 바이오마커 — PGI/PGII 비율(위축 지표), CagA 상태

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("cgast_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("cgast_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg cgast_qsp_model.dot -o cgast_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [cgast_references.md](cgast_references.md) 참조 (약 61편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
