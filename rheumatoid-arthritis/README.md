# 류마티스 관절염 (Rheumatoid Arthritis, RA) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 자가면역·류마티스

[![RA QSP Model](ra_qsp_model.png)](ra_qsp_model.svg)

## 개요 (Overview)

류마티스 관절염(RA)은 전 세계 성인 인구의 약 0.5~1%에서 발생하는 만성 전신 자가면역 관절염으로, 여성에서 2~3배 많습니다. 활막 내 판누스(pannus) 형성, TNF-α·IL-6·IL-1β에 의한 지속적 관절 염증, RANKL 유도 파골세포 활성화로 인한 골·연골 미란이 핵심 발병기전입니다. 혈청 RF·anti-CCP 항체가 진단 마커이자 병인 기전의 일부를 구성합니다. 항TNF·항IL-6·CTLA-4-Ig·CD20·JAK 억제제 등 다양한 생물학제제·표적치료제가 표준 치료로 확립되어 있습니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| TNF-α/NF-κB 경로 | 대식세포·FLS의 TNF → NF-κB → IL-6·IL-1β·MMP 연쇄 | 활막 비후, 판누스, 골미란 |
| IL-6/JAK-STAT3 경로 | IL-6 → IL-6R/gp130 → JAK1-STAT3 → 급성기 반응 | CRP 상승, 전신 증상, 빈혈 |
| T세포 공자극 (CD28/CTLA-4) | CD4+ Th1·Th17 → 시노비오사이트·B세포 활성화 | 자가항체 생성, 관절 침윤 |
| B세포/자가항체 | RF·anti-CCP 면역복합체 → 보체 활성화 → FcγR 신호 | 관절 내 보체 손상 |
| RANKL/OPG 불균형 | IL-17·TNF → FLS의 RANKL 과발현 → 파골세포 | 방사선 골 미란, 골다공증 |
| FLS 침습성 증식 | Akt/PI3K·Wnt 신호 → 활막 섬유아세포(FLS) 판누스 형성 | 연골 침식, 관절 변형 |
| JAK-STAT 통합 | IL-2·IFN-γ·GM-CSF → JAK1/2/3 → 다중 사이토카인 생성 | 만성 활막 염증 유지 |

## 주요 약물 표적 (Drug Targets)

- **항TNF 제제 (메토트렉세이트 병용: 에타너셉트, 인플릭시맙, 아달리무맙, 세르톨리주맙, 골리무맙)**: TNF 직접 중화; 방사선 진행 억제 근거 강력
- **항IL-6R (토실리주맙, 사릴루맙)**: IL-6 신호 차단 → CRP·급성기 반응 정상화; MTX 병용 또는 단독
- **CTLA-4-Ig (아바타셉트)**: CD80/86 차단 → T세포 공자극 억제; RF·anti-CCP 양성에서 효과 우수
- **항CD20 (리툭시맙)**: B세포 고갈 → 자가항체 감소; 항TNF 실패 후 주요 대안
- **JAK 억제제 (토파시티닙, 바리시티닙, 우파다시티닙, 필고티닙)**: JAK1/2/3 선택 억제; 경구 투여 편의성
- **csDMARD (메토트렉세이트, 레플루노미드, 하이드록시클로로퀸, 설파살라진)**: 기반 치료, 조기 RA 1차

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [ra_qsp_model.dot](ra_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 100+ 노드 / 11 클러스터) |
| [ra_qsp_model.svg](ra_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [ra_qsp_model.png](ra_qsp_model.png) | PNG 이미지 (150 dpi) |
| [ra_mrgsolve_model.R](ra_mrgsolve_model.R) | mrgsolve ODE 모델 (약 20 구획 / 약 29개 시나리오) |
| [ra_shiny_app.R](ra_shiny_app.R) | Shiny 대시보드 |
| [ra_references.md](ra_references.md) | 참고문헌 (약 42편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 항TNF/항IL-6/CTLA-4-Ig/JAK억제제/MTX의 PK 구획 + TNF·IL-6·IL-17 농도 동태, CD4+ T세포·FLS 활성화 모듈, 관절 손상 지수(DAS28·HAQ), 방사선 진행 mTSS 구획 포함
- **주요 치료 시나리오**: 무치료 자연 경과, MTX 단독, 에타너셉트+MTX, 아달리무맙+MTX, 토실리주맙 단독/병용, 아바타셉트, 리툭시맙, 토파시티닙, 바리시티닙, 우파다시티닙 등
- **보정/근거**: TEMPO(에타너셉트+MTX), OPTION(토실리주맙), AIM(아바타셉트), REFLEX(리툭시맙), RA-BEAM(바리시티닙), ORAL Scan(토파시티닙) 임상 데이터 기반

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) **환자 프로파일** — 혈청 양성·DAS28 기저값·RF/anti-CCP 설정; (2) **PK 프로파일** — 생물학제제/소분자 혈중 농도 시간 경과; (3) **PD 주요지표** — TNF·IL-6 억제율, CRP 동태; (4) **임상 엔드포인트** — DAS28·ACR20/50/70·HAQ·mTSS 경시 변화; (5) **시나리오 비교** — 치료 전략 간 관해율·방사선 진행 비교; (6) **바이오마커** — RF, anti-CCP, ESR, CRP, 관절 초음파 스코어.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("ra_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("ra_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg ra_qsp_model.dot -o ra_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [ra_references.md](ra_references.md) 참조 (약 42편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
