# 폐쇄성 수면 무호흡 (Obstructive Sleep Apnea, OSA) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 호흡기

[![OSA QSP Model](osa_qsp_model.png)](osa_qsp_model.svg)

## 개요 (Overview)

폐쇄성 수면 무호흡(OSA)은 수면 중 반복적인 상기도 허탈로 인해 무호흡(≥10초) 및 저호흡이 발생하는 질환으로, 성인 유병률은 약 10–30%(중증 이상 약 4–7%)에 달한다. 비만, 상기도 구조적 이상(작은 하악·편도 비대), 상기도 확장근 긴장도 저하가 핵심 위험 인자이다. 간헐적 저산소증(IH)은 교감신경 활성화, HIF-1α 유도, 전신 염증, 산화 스트레스를 통해 고혈압·심방세동·대사 증후군·인지 기능 저하의 위험을 높인다. 양압기(CPAP)가 1차 치료이며, 비CPAP 표현형(중추성·혼합성·낮은 각성 역치)에는 약물 보조 요법이 점차 활용되고 있다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 상기도 허탈 | Pcrit 상승 (critical closing pressure) → 흡기 시 음압에 의한 허탈 | AHI 증가, 무호흡·저호흡 |
| 간헐적 저산소 | SpO2 반복 저하 → HIF-1α↑ → VEGF, EPO | 교감신경 활성, 야간 저산소 |
| 교감신경 과활성 | NE 분비↑ → 혈관 수축, HR·BP 상승 | 고혈압, 심방세동 위험 |
| 전신 염증 | IL-6, CRP, TNF-α 상승 → 내피 기능 장애 | 동맥경화, 심혈관 질환 |
| 인슐린 저항성 | HIF-1α → 글루코코르티코이드 증가, 지방분해↑ | HOMA-IR 증가, 대사 증후군 |
| 수면 분절 | 각성 반응(arousal) → 서파·REM 수면 감소 | 주간 졸림(ESS), 인지 저하 |
| 비만-OSA 악순환 | 체지방↑ → 상기도 지방 침착 → Pcrit↑ | 비만과 OSA 상호 악화 |

## 주요 약물 표적 (Drug Targets)

- **CPAP (양압기)**: 지속적 양압으로 상기도 개방 유지; AHI 95% 이상 감소, 혈압 저하
- **모다피닐/솔리암페톨 (각성 촉진제)**: 도파민 재흡수 차단·NE 증가 → 잔류 주간 졸림(EDS) 개선
- **티르제파타이드 (GLP-1/GIPR 이중 작용제)**: 비만 OSA에서 체중 감량 → AHI 유의미 감소 (SURMOUNT-OSA)
- **에스조피클론 (비벤조디아제핀계 수면제)**: 각성 역치 증가 → 일부 OSA 표현형에서 AHI 개선
- **아세타졸아미드 (탄산탈수효소 억제제)**: 중추성 무호흡·혼합형 OSA 표현형에서 고리통기(loop gain) 감소

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [osa_qsp_model.dot](osa_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 178 노드 / 10 클러스터) |
| [osa_qsp_model.svg](osa_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [osa_qsp_model.png](osa_qsp_model.png) | PNG 이미지 (150 dpi) |
| [osa_mrgsolve_model.R](osa_mrgsolve_model.R) | mrgsolve ODE 모델 (약 24 구획 / 6개 치료 시나리오) |
| [osa_shiny_app.R](osa_shiny_app.R) | Shiny 대시보드 |
| [osa_references.md](osa_references.md) | 참고문헌 (약 65편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 모다피닐·솔리암페톨 경구 2구획씩, 티르제파타이드 SC 3구획, 에스조피클론·아세타졸아미드 경구 각 2구획; Pcrit(상기도 임계압)·상기도 이득(LG)·각성 역치(AT)·AHI·SpO2·HIF-1α·교감신경 활성(SNA)·수축기혈압·심박수·CRP·HOMA-IR·ESS·체중
- **주요 치료 시나리오**: ① 무치료(AHI=35) ② CPAP 10 cmH2O ③ CPAP+모다피닐 200mg ④ CPAP+솔리암페톨 150mg ⑤ 티르제파타이드 10mg QW(비CPAP) ⑥ 아세타졸아미드+에스조피클론(비CPAP 표현형)
- **보정/근거**: SURMOUNT-OSA(티르제파타이드), TONES(솔리암페톨), CPAP AHI 반응 문헌 데이터 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① **환자 프로파일** (기저 AHI·BMI·고혈압 유무·OSA 표현형 설정), ② **PK** (혈장 약물 농도-시간 곡선), ③ **PD 주요지표** (AHI·SpO2·Pcrit·교감신경 활성 추이), ④ **임상 엔드포인트** (ESS·혈압·HOMA-IR·CRP 변화), ⑤ **시나리오 비교** (6개 치료 전략 직접 비교), ⑥ **바이오마커** (HIF-1α·체중·심박수·수면 분절 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("osa_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("osa_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg osa_qsp_model.dot -o osa_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [osa_references.md](osa_references.md) 참조 (약 65편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
