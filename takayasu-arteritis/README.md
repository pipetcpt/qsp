# 다카야스 동맥염 (Takayasu Arteritis, TA) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 혈관염

[![TA QSP Model](ta_qsp_model.png)](ta_qsp_model.svg)

## 개요 (Overview)

다카야스 동맥염은 대동맥 및 주요 분지를 침범하는 만성 육아종성 대혈관염으로, 전 세계 연간 발생률은 인구 100만 명당 약 1~3명이며 40세 미만 젊은 여성에서 주로 발생한다. 혈관벽 전층에 T세포와 대식세포가 침윤하여 육아종을 형성하고, 내막 증식과 중막 섬유화로 혈관 협착·폐쇄·동맥류를 초래한다. IL-6/JAK-STAT 축과 TNF-α가 핵심 염증 매개체이며, 토실리주맙과 스테로이드가 현재 주요 치료 전략이다. 치료 없이 방치하면 뇌졸중, 심부전, 신혈관성 고혈압 등 심각한 혈관 합병증이 발생한다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 선천 면역 활성화 | 수지상세포·대식세포 혈관벽 침윤, Toll-like receptor 활성화 | 초기 혈관벽 염증 |
| Th1/Th17 분화 | IL-6, IL-17, IL-12, IFN-γ 과분비 | 육아종 형성, 혈관 손상 지속 |
| IL-6/STAT3 축 | sIL-6R 매개 trans-signaling, CRP 유도 | 전신 급성기 반응, 혈관염 활성도 지표 |
| TNF-α 경로 | TNF-α → NF-κB → 내막 증식, VCAM-1 발현 | 혈관 협착 진행, 조직 허혈 |
| 혈관 리모델링 | TGF-β → 중막 섬유화, neointima 형성 | 내강 협착, 측부 순환 발달 |
| 혈류역학 변화 | 협착 → PET/MRA 활성 병변, 혈압 비대칭 | NIH 점수·혈관 활성도 지표 상승 |

## 주요 약물 표적 (Drug Targets)

- **IL-6R 차단제** — 토실리주맙(Tocilizumab) 162 mg SC qw: IL-6 trans-signaling 차단, 재발률 감소 (TAKT trial)
- **코르티코스테로이드** — 프레드니솔론 1 mg/kg/day → 점감: 급성 염증 억제, 육아종 형성 감소
- **메토트렉세이트(MTX)** — 면역억제 보조요법: 스테로이드 감량 가능케 하는 steroid-sparing effect
- **인플릭시맙(Infliximab)** — TNF-α 차단, 불응성 환자에서 사용
- **JAK 억제제** — 바리시티닙·토파시티닙(연구 단계): JAK1/2 억제로 IL-6·IFN-γ 신호 차단

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [ta_qsp_model.dot](ta_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 373 노드 / 11 클러스터) |
| [ta_qsp_model.svg](ta_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [ta_qsp_model.png](ta_qsp_model.png) | PNG 이미지 (150 dpi) |
| [ta_mrgsolve_model.R](ta_mrgsolve_model.R) | mrgsolve ODE 모델 (약 23 구획 / 5 치료 시나리오) |
| [ta_shiny_app.R](ta_shiny_app.R) | Shiny 대시보드 |
| [ta_references.md](ta_references.md) | 참고문헌 (약 58편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK 구획(프레드니솔론 3구획, 토실리주맙 SC 3구획, MTX 경구·혈장·폴리글루타메이트 3구획, 인플릭시맙 2구획) + 질환 PD 구획(IL-6, sIL-6R, IL-6복합체, TNF-α, Th1, Th17, Treg, 혈관벽 침윤, ST, CRP, PET활성도, 혈관벽두께)으로 총 약 23개 구획
- **주요 치료 시나리오**: (1) 스테로이드 단독(자연경과), (2) 프레드니솔론 + MTX, (3) 프레드니솔론 + 토실리주맙(TAKT trial 설계), (4) 프레드니솔론 + 인플릭시맙, (5) 치료 없음(자연경과)
- **보정/근거**: TAKT 임상시험(토실리주맙) 재발률·PET 활성도 데이터, ACR/EULAR 가이드라인 기준 반응 정의를 근거로 파라미터 보정

## Shiny 대시보드 (Dashboard)

환자 프로파일(연령·성별·질환 아형·기저 NIH 점수) · PK 시각화(각 약물 혈중 농도-시간 곡선) · 염증 PD 지표(CRP·IL-6·PET 활성도 시계열) · 혈관 구조 지표(혈관벽 두께·협착 정도) · 치료 시나리오 비교(단독 vs 병용요법 재발률) · 바이오마커 패널(Th1/Th17/Treg 균형) 등 6개 이상의 탭으로 구성된다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("ta_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("ta_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg ta_qsp_model.dot -o ta_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [ta_references.md](ta_references.md) 참조 (약 58편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
