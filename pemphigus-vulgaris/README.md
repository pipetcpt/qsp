# 심상성 천포창 (Pemphigus Vulgaris, PV) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 피부

[![PV QSP Model](pv_qsp_model.png)](pv_qsp_model.svg)

## 개요 (Overview)

심상성 천포창(PV)은 표피 세포 간 접착 단백질인 데스모글레인-3(Dsg3), 경우에 따라 데스모글레인-1(Dsg1)에 대한 자가항체(주로 IgG4)가 케라티노사이트 세포-세포 접착을 파괴하여 표피내 수포 및 점막 미란을 일으키는 희귀 자가면역 수포성 질환이다. 전 세계 연간 발생률은 약 100만 명당 1–5명으로 추정되며, 아슈케나지 유대인 및 지중해계에서 유전적으로 더 흔하다. 치료받지 않으면 광범위한 피부 박탈과 패혈증으로 사망률이 높으므로, 스테로이드와 리툭시맙이 핵심 치료제이다. 최근 에프가르티기모드(항FcRn)의 적응증 확대 연구가 진행 중이다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| TFH-GC B세포 축 | TFH(IL-21)→ GC B세포 → 기억 B세포·형질세포 → 항Dsg3 IgG | 항체 지속 생산, 재발 |
| 항Dsg3 IgG4 결합 | IgG4 항체 → Dsg3 표적 파괴, p38-MAPK 활성화 → 데스모좀 분해 | 표피내 수포(acantholysis) |
| 보체 활성화 | Dsg3-IgG1 복합체 → C1q → C3b → 염증 증폭 | 표피 손상 가속 |
| Treg 기능 저하 | CD4+Foxp3+ Treg↓ → TFH·GC B세포 억제 실패 | 자가면역 관용 파괴 |
| FcRn-IgG 재순환 | FcRn이 IgG 분해 방지 → 항Dsg3 항체 반감기 연장 | 고역가 항체 유지 |
| 각질세포 신호 | Dsg3 소실 → EGFR·p38 → 세포 수축·아폽토시스 | 수포 형성, 미란 확대 |
| 장기 형질세포(LLPC) | 골수 LLPC → 스테로이드·리툭시맙 내성 항체 생산 | 재발, 치료 저항성 |

## 주요 약물 표적 (Drug Targets)

- **고용량 프레드니솔론 (Prednisolone)**: 전반적 면역억제; TFH·GC B세포 억제, Treg 확장 → 1차 치료
- **리툭시맙 (Rituximab)**: 항CD20 B세포 고갈; RITUX3 임상시험에서 스테로이드 단독 대비 CR율 우월
- **마이코페놀레이트 모페틸 (MMF)**: 퓨린 합성 억제 → B세포·T세포 증식 차단; 유지 면역억제
- **에프가르티기모드 (Efgartigimod)**: 항FcRn → IgG(항Dsg3 항체 포함) 혈중 농도 신속 감소; BALLAD 연구
- **IV 메틸프레드니솔론 펄스**: 중증 급성기 신속 관해 유도; 고용량 스테로이드 이후 감량

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [pv_qsp_model.dot](pv_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 159 노드 / 12 클러스터) |
| [pv_qsp_model.svg](pv_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [pv_qsp_model.png](pv_qsp_model.png) | PNG 이미지 (150 dpi) |
| [pv_mrgsolve_model.R](pv_mrgsolve_model.R) | mrgsolve ODE 모델 (약 23 구획 / 6개 치료 시나리오) |
| [pv_shiny_app.R](pv_shiny_app.R) | Shiny 대시보드 |
| [pv_references.md](pv_references.md) | 참고문헌 (약 40편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 프레드니솔론 경구 2구획+말초, 리툭시맙 IV 2구획, MMF 경구 1구획, IVIg 1구획, 에프가르티기모드 1구획+FcRn; 미성숙 B세포·GC B세포·기억 B세포·SLPC·LLPC·TFH·Treg, 항Dsg3 IgG, Dsg3 발현, 수포/PDAI 점수, 보체 활성화, 코르티코스테로이드 누적 골손실
- **주요 치료 시나리오**: ① 고용량 CS 단독(기존 표준) ② 리툭시맙+저용량 CS(RITUX3) ③ MMF+중등도 CS ④ 리툭시맙+MMF 병용 ⑤ IV 메틸프레드니솔론 펄스+리툭시맙(중증) ⑥ 에프가르티기모드+저용량 CS
- **보정/근거**: RITUX3(리툭시맙+저용량 CS vs 고용량 CS), BALLAD(에프가르티기모드) 임상시험 데이터 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① **환자 프로파일** (기저 항Dsg3 역가·PDAI 점수·점막 침범 여부 설정), ② **PK** (혈장 약물 농도 및 FcRn 점유율, 리툭시맙 CD20 점유율), ③ **PD 주요지표** (항Dsg3 IgG·GC B세포·LLPC 추이), ④ **임상 엔드포인트** (PDAI·CR/PR율·재발 시점), ⑤ **시나리오 비교** (6개 치료 전략 직접 비교), ⑥ **바이오마커** (항Dsg3 IgG4·IgG1·누적 스테로이드 독성 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("pv_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("pv_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg pv_qsp_model.dot -o pv_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [pv_references.md](pv_references.md) 참조 (약 40편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
