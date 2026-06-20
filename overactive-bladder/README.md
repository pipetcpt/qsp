# 과민성 방광 (Overactive Bladder, OAB) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 신장·비뇨

[![OAB QSP Model](oab_qsp_model.png)](oab_qsp_model.svg)

## 개요 (Overview)

과민성 방광(OAB)은 절박뇨를 주증상으로 빈뇨·야간뇨·절박성 요실금을 동반하는 하부요로증상 증후군으로, 40세 이상 성인의 약 11–16%에서 발생하고 고령일수록 유병률이 증가한다. 핵심 발병기전은 방광 배뇨근의 과활동(involuntary detrusor contraction)으로, 무스카린 수용체(주로 M2/M3)의 과잉 자극과 β3-아드레날린 수용체 신호 저하, 방광 요로상피에서의 ATP-퓨린수용체 및 NGF-TrkA 신호 이상이 기여한다. 일차 치료는 행동 요법 및 골반저 근육 운동이며, 약물로는 항무스카린제와 미라베그론(β3 작용제)이 사용된다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| M2/M3 과활성 | ACh→M3 수용체 → IP3/DAG → 배뇨근 수축 | 절박뇨, 절박성 요실금 |
| β3-아드레날린 신호 저하 | NE→β3-AR → cAMP↓ → 배뇨근 이완 감소 | 방광 용량 감소 |
| 요로상피-방광 감각 | 요로상피 ATP→P2X3 → 구심성 C/Aδ섬유 활성화 | 방광 과민, 절박감 |
| NGF 상승 | 배뇨근 과활동 → NGF↑ → TrkA → 신경 과성장 | 감각 신경 민감화 |
| 방광벽 섬유화 | 만성 허혈·염증 → TGF-β → 콜라겐↑ | 방광 순응도 감소 |
| 척수-뇌 배뇨 반사 | 척수 PAG-PMC 회로 과활성 | 배뇨 억제 실패, 절박감 |
| 비교적 방광 용량 감소 | 배뇨근 두께 증가, 방광 용량↓ | 빈뇨, OABq 점수 증가 |

## 주요 약물 표적 (Drug Targets)

- **옥시부티닌 (Oxybutynin)**: M1/M2/M3 비선택적 항무스카린제; 절박뇨·빈뇨 감소, 구강건조·인지 부작용
- **톨테로딘 (Tolterodine)**: M2/M3 선택적 항무스카린제; 방광 선택성 개선
- **솔리페나신 (Solifenacin)**: M3 선택적 항무스카린제; 빈뇨 및 절박성 요실금 감소
- **미라베그론 (Mirabegron)**: β3-AR 선택적 작용제 → 배뇨근 이완·방광 용량 증가; 항무스카린 부작용 없음
- **미라베그론+솔리페나신 복합**: 서로 다른 기전 병용으로 상가 또는 상승 효과 (SYNERGY 연구)

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [oab_qsp_model.dot](oab_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 146 노드 / 13 클러스터) |
| [oab_qsp_model.svg](oab_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [oab_qsp_model.png](oab_qsp_model.png) | PNG 이미지 (150 dpi) |
| [oab_mrgsolve_model.R](oab_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 6개 치료 시나리오) |
| [oab_shiny_app.R](oab_shiny_app.R) | Shiny 대시보드 |
| [oab_references.md](oab_references.md) | 참고문헌 (약 42편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 옥시부티닌·톨테로딘·솔리페나신·미라베그론·솔리페나신2 경구 각 2구획; M3 수용체 점유율(RO_M3)·β3 수용체 점유율(RO_B3)·배뇨근 활성(DetAct)·방광 용량(BladCap)·배뇨 빈도·절박감·절박성 요실금(UUI)·NGF·ATP·야간뇨·OABq 점수
- **주요 치료 시나리오**: ① 무치료 ② 옥시부티닌 IR 5mg TID ③ 톨테로딘 ER 4mg QD ④ 솔리페나신 10mg QD ⑤ 미라베그론 50mg QD ⑥ 미라베그론+솔리페나신 병용
- **보정/근거**: ARIES(솔리페나신), BEYOND(미라베그론), SYNERGY(복합요법) 임상시험 데이터 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① **환자 프로파일** (기저 배뇨 빈도·절박성 요실금 횟수·성별·연령 설정), ② **PK** (혈장 약물 농도-시간 곡선 및 M3/β3 수용체 점유율), ③ **PD 주요지표** (배뇨근 활성·방광 용량 추이), ④ **임상 엔드포인트** (배뇨 빈도·UUI 횟수·야간뇨·OABq 점수 변화), ⑤ **시나리오 비교** (6개 치료 전략 직접 비교), ⑥ **바이오마커** (NGF·ATP 방광 점막 신호, 방광 순응도 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("oab_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("oab_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg oab_qsp_model.dot -o oab_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [oab_references.md](oab_references.md) 참조 (약 42편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
