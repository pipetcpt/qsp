# 국소분절사구체경화증 (Focal Segmental Glomerulosclerosis, FSGS) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 신장·비뇨

[![FSGS QSP Model](fsgs_qsp_model.png)](fsgs_qsp_model.svg)

## 개요 (Overview)

국소분절사구체경화증(FSGS)은 족세포(podocyte) 손상을 특징으로 하는 신증후군의 중요한 원인 질환으로, 성인 신증후군의 약 35~50%를 차지합니다. 원발성 FSGS에서는 순환 투과인자(permeability factor, 가설적으로 suPAR·CLCF1 등)가 족세포 표면 수용체를 자극하여 발족 돌기 소실(foot process effacement)과 단백뇨를 유발합니다. 치료하지 않으면 수년 내 말기신부전으로 진행하며, RAAS 억제·스테로이드·칼시뉴린억제제·리툭시맙이 주요 치료 옵션입니다. 최근 스파르센탄(DUPLEX 시험)이 새로운 표적치료제로 승인되었습니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 순환 투과인자 | CLCF1·suPAR → 족세포 β3 인테그린 활성화 | 발족 돌기 소실, 단백뇨 |
| 족세포 손상 | nephrin·podocin·synaptopodin 소실 → 사구체 여과장벽 파괴 | 단백뇨 3.5 g/일 이상 |
| 사구체 경화 | TGF-β → 메산지움 세포 활성화 → 콜라겐 침착 | 사구체 분절 경화, GFR 저하 |
| 보체 경로 | C3b·MAC 형성 → 족세포 세포사 | 사구체 손상 가속 |
| RAAS 과활성화 | Ang II → 사구체내압 상승 → 단백뇨 가중 | 신기능 저하 촉진 |
| 만성 신장 염증 | 단백뇨 → 세관주위 염증 → 간질 섬유화 | eGFR 진행성 감소 |
| 이차 FSGS | 비만·단일신·고혈압 → 사구체 비대·경화 | 단백뇨, 만성신부전 |

## 주요 약물 표적 (Drug Targets)

- **RAAS 억제제** (ACEi/ARB): 사구체내압 감소, 단백뇨 억제 — 모든 FSGS의 기본 치료
- **프레드니솔론(고용량 스테로이드)**: 족세포 보호, 단백뇨 감소 — 원발성 FSGS 1차
- **칼시뉴린억제제** (타크로리무스, 사이클로스포린): 족세포 직접 보호 + 면역억제 — 스테로이드 저항성
- **리툭시맙**: B세포 고갈, 족세포 직접 안정화(sphingomyelin phosphodiesterase 경로)
- **스파르센탄(Sparsentan)**: 이중 길항제(AT1R + 엔도텔린 ETA 수용체) — DUPLEX 시험

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [fsgs_qsp_model.dot](fsgs_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 143 노드 / 12 클러스터) |
| [fsgs_qsp_model.svg](fsgs_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [fsgs_qsp_model.png](fsgs_qsp_model.png) | PNG 이미지 (150 dpi) |
| [fsgs_mrgsolve_model.R](fsgs_mrgsolve_model.R) | mrgsolve ODE 모델 (약 23 구획 / 6 치료 시나리오) |
| [fsgs_shiny_app.R](fsgs_shiny_app.R) | Shiny 대시보드 |
| [fsgs_references.md](fsgs_references.md) | 참고문헌 (약 67편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 질환 상태 구획(순환 투과인자, 족세포 분율, 발족 돌기 소실 지수, 단백뇨, eGFR, 경화 분율, TGF-β, 보체 활성, 사구체 염증) + 약물 PK 구획(프레드니솔론·타크로리무스·리툭시맙·스파르센탄 PK)
- **주요 치료 시나리오**: ① 자연 경과(무치료), ② 프레드니솔론 단독, ③ 프레드니솔론 + 타크로리무스(스테로이드 저항성 프로토콜), ④ 프레드니솔론 + 타크로리무스 + 리툭시맙(불응성 FSGS), ⑤ 스파르센탄 단독(DUPLEX 시험 기반), ⑥ 완전 복합 요법
- **보정/근거**: DUPLEX 시험(스파르센탄), Cattran FSGS 국제 등록자료, D'Agati 분류 및 FSGS 자연 경과 문헌 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① 환자 프로파일(기저 단백뇨·eGFR·suPAR·FSGS 아형 설정), ② PK 탭(각 약물 혈중 농도·타크로리무스 트로프 추적), ③ 족세포·신장 PD 탭(족세포 분율·단백뇨·eGFR 추이), ④ 임상 엔드포인트(완전/부분 관해·eGFR 50% 감소 도달 시간), ⑤ 시나리오 비교(6개 치료 전략 동시 비교), ⑥ 바이오마커(TGF-β·보체·경화 분율·세관 손상 지수 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("fsgs_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("fsgs_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg fsgs_qsp_model.dot -o fsgs_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [fsgs_references.md](fsgs_references.md) 참조 (약 67편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
