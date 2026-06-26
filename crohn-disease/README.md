# 크론병 (Crohn's Disease, CD) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도

[![CD QSP Model](cd_qsp_model.png)](cd_qsp_model.svg)

## 개요 (Overview)

크론병은 소장·대장을 포함한 소화관 어느 부위에도 발생할 수 있는 전층성(transmural) 만성 염증성 장 질환으로, 전 세계적으로 약 300만 명 이상이 이환되어 있으며 젊은 성인에서 호발한다. Th1/Th17 주도의 과도한 면역 반응이 TNF-α·IL-12/23·IL-17의 사이토카인 네트워크를 통해 장벽 투과성 증가와 전층성 육아종 형성을 유발한다. 항TNF 생물학제제(인플릭시맙·아달리무맙)가 생물학적 치료의 기반이며, 항IL-12/23(우스테키누맙), 항인테그린(베돌리주맙), JAK 억제제(우파다시티닙)가 추가 치료 옵션으로 확립되어 있다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| Th1 과활성화 | IL-12 → Th1 분화 → IFN-γ, TNF-α 분비 | 전층성 육아종 염증 |
| Th17 경로 | IL-23 → Th17 → IL-17, IL-22 분비 | 장벽 손상, 점막 면역 이상 |
| TNF-α 중심 염증 | TNF-α → NF-κB → 사이토카인 캐스케이드 | 궤양, 협착, 누공 형성 |
| 장벽 기능 저하 | 밀착연접 단백질(ZO-1·클로딘) 감소 | 세균 전위, 내독소혈증 |
| 골수 이동·호중구 침윤 | IL-8, MIP-1α → 호중구 모집 | 점막 손상, 분변 칼프로텍틴 ↑ |
| 장간막 지방조직 염증 | 크리핑 지방 → 아디포카인·사이토카인 분비 | 질환 위치 결정, 협착 촉진 |
| 골 대사 이상 | 스테로이드 + 만성 염증 → BMD 감소 | 골다공증 위험 |

## 주요 약물 표적 (Drug Targets)

- **항TNF 항체(인플릭시맙·아달리무맙·세르톨리주맙)**: TNF-α 중화 → NF-κB 억제, 점막 치유 촉진
- **우스테키누맙(항IL-12/23)**: p40 소단위 차단 → Th1/Th17 분화 억제 (UNIFI/CERTIFI)
- **베돌리주맙(항α4β7 인테그린)**: 장 선택적 림프구 귀소 차단 (GEMINI 시험)
- **우파다시티닙(JAK1 억제제)**: IL-6·IL-12·IL-23 신호 억제 → 범사이토카인 억제 (U-EXCEL)
- **아자티오프린/6-MP(티오퓨린)**: 퓨린 합성 억제 → 증식 림프구 억제, 생물학제제와 병용으로 항체 형성 억제

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [cd_qsp_model.dot](cd_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 392 노드 / 10 클러스터) |
| [cd_qsp_model.svg](cd_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [cd_qsp_model.png](cd_qsp_model.png) | PNG 이미지 (150 dpi) |
| [cd_mrgsolve_model.R](cd_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 8개 치료 시나리오) |
| [cd_shiny_app.R](cd_shiny_app.R) | Shiny 대시보드 |
| [cd_references.md](cd_references.md) | 참고문헌 (약 61편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 인플릭시맙·아달리무맙·우스테키누맙·베돌리주맙·티오구아닌·프레드니솔론·우파다시티닙 PK 구획(10구획) + TNF-α, IL-12/23, IL-17, Th17, Th1, Treg, 호중구, 점막 염증 지수(MI), CRP, 분변 칼프로텍틴(FC), 골밀도, 헤모글로빈 PD 구획(12구획)
- **주요 치료 시나리오**: (1) 무치료, (2) 스테로이드 단독, (3) 인플릭시맙 유도-유지, (4) 아달리무맙 SC, (5) 우스테키누맙 IV→SC, (6) 베돌리주맙 IV, (7) 아자티오프린 + 인플릭시맙, (8) 우파다시티닙 경구
- **보정/근거**: 인플릭시맙 PK는 Ng CM et al. (Clin Pharmacokinet 2010), 임상 반응률은 ACCENT I/II, CHARM, UNIFI 시험 데이터 기반

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) 환자 프로파일 — 병변 위치·질환 중증도·면역원성 위험 설정; (2) PK 탭 — 생물학제제·소분자 혈장 농도 시계열; (3) PD 주요지표 — TNF-α, IL-12/23, Th1/Th17/Treg, CRP; (4) 임상 엔드포인트 — CDAI, 분변 칼프로텍틴, 점막 치유; (5) 시나리오 비교 — 8가지 치료 전략 1년 결과; (6) 바이오마커 — CRP, 분변 칼프로텍틴, 헤모글로빈, 골밀도, 약물 농도 TDM

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("cd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("cd_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg cd_qsp_model.dot -o cd_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [cd_references.md](cd_references.md) 참조 (약 61편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
