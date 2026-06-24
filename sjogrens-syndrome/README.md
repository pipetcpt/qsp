# 쇼그렌 증후군 (Sjögren's Syndrome, SS) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 자가면역·류마티스

[![SS QSP Model](ss_qsp_model.png)](ss_qsp_model.svg)

## 개요 (Overview)

쇼그렌 증후군은 주로 침샘·눈물샘 등 외분비선을 CD4+ T림프구와 B세포가 침윤하여 만성 건조증(입마름·눈마름)을 일으키는 전신 자가면역 질환으로, 유병률은 약 0.1~0.6%이며 여성에서 9:1로 압도적으로 많습니다. I형 인터페론(IFN-α/β) 신호 과활성 및 B세포 과활성화(RF·항SSA/Ro·항SSB/La 자가항체 생성)가 핵심 병인입니다. 림프종(주로 MALT형) 합병 위험이 건강인 대비 약 15~20배 증가하는 것이 주요 장기 합병증입니다. 현재 승인된 질환 수정 치료는 없으며, 하이드록시클로로퀸·피로카르핀·세비멜린 대증 치료, 전신 침범에는 면역억제제·리툭시맙이 오프라벨로 사용됩니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| I형 IFN 경로 | pDC의 TLR7/9 → IFN-α/β 과생성 → ISG 발현 → 선 손상 | 내피·외분비선 기능 저하 |
| B세포 과활성화 | BAFF 과발현 → B세포 생존·증식 → RF·항SSA/SSB 자가항체 | 면역복합체 침착, 보체 활성화 |
| CD4+ T세포(Th1/Tfh) 침윤 | IFN-γ·IL-21 → 침샘 내 림프구 병소(focus) 형성 | 외분비선 파괴 진행 |
| BAFF/APRIL 경로 | BAFF 과발현 → B세포 분화·자가항체 지속 | B세포 림프종 위험 증가 |
| 머스카린 수용체 억제 | 항M3R 자가항체 → 침샘·눈물샘 분비 억제 | 구강 건조증, 안구 건조증 |
| 신장·신경 침범 | 면역복합체 → 세뇨관 산증, 신경초 침범 | 저칼륨혈증, 말초신경병증 |
| 림프종 전환 경로 | 만성 B세포 활성화 → MALT 림프종 클론 증식 | 치명적 합병증 |

## 주요 약물 표적 (Drug Targets)

- **하이드록시클로로퀸 (Hydroxychloroquine)**: TLR9·IFN-α 억제; 피로·관절통·경증 전신 증상 완화
- **피로카르핀/세비멜린**: M3R 작용제 → 침샘·눈물샘 분비 자극; 구강·안구 건조증 대증 치료
- **리툭시맙 (Rituximab)**: 항CD20 → B세포 고갈 → 자가항체 감소; 혈관염·신증·신경 침범 적응증
- **벨리무맙 (Belimumab)**: 항BAFF → B세포 생존 억제; pSS 적응증 임상시험 진행 중
- **JAK 억제제 (바리시티닙, 우파다시티닙)**: JAK1/2 → I형 IFN·IL-6 신호 억제; 임상시험 단계
- **아바타셉트**: CD80/86 차단 → Tfh·형질세포 활성화 억제; pSS 2상 임상 결과 일부 긍정적

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [ss_qsp_model.dot](ss_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 100+ 노드 / 12 클러스터) |
| [ss_qsp_model.svg](ss_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [ss_qsp_model.png](ss_qsp_model.png) | PNG 이미지 (150 dpi) |
| [ss_mrgsolve_model.R](ss_mrgsolve_model.R) | mrgsolve ODE 모델 (약 19 구획 / 약 29개 시나리오) |
| [ss_shiny_app.R](ss_shiny_app.R) | Shiny 대시보드 |
| [ss_references.md](ss_references.md) | 참고문헌 (약 47편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 리툭시맙/벨리무맙/하이드록시클로로퀸 PK 구획 + I형 IFN 생성 동태, BAFF·B세포 수 모듈, 자가항체(항SSA) 역가 구획, 침샘 기능(USFR, unstimulated salivary flow rate) 예측 모듈, ESSDAI 점수 추정 구획 포함
- **주요 치료 시나리오**: 무치료 자연 경과, 하이드록시클로로퀸, 피로카르핀, 메토트렉세이트, 아자티오프린, 리툭시맙 1회, 리툭시맙 반복, 벨리무맙, 바리시티닙, 리툭시맙+벨리무맙 병용 등
- **보정/근거**: TRACTISS(리툭시맙 RCT), BELISS(벨리무맙), EULAR pSS 권고안, Vitali 분류기준 임상 코호트 데이터 기반

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) **환자 프로파일** — 항SSA/SSB 상태·ESSDAI 기저값·림프구 침윤 초점 점수 설정; (2) **PK 프로파일** — 리툭시맙/벨리무맙/소분자 혈중 농도; (3) **PD 주요지표** — IFN 점수·BAFF·B세포 수 억제 동태; (4) **임상 엔드포인트** — ESSDAI·ESSPRI·침샘 유량·안구 건조 점수 경시 변화; (5) **시나리오 비교** — 치료 전략별 선 기능·전신 활성도 비교; (6) **바이오마커** — 항SSA, IgG, 보체 C3/C4, IFN 신호 점수.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("ss_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("ss_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg ss_qsp_model.dot -o ss_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [ss_references.md](ss_references.md) 참조 (약 47편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
