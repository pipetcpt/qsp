# 성인형 스틸병 (Adult-Onset Still's Disease, AOSD) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 자가면역·류마티스

[![AOSD QSP Model](aosd_qsp_model.png)](aosd_qsp_model.svg)

## 개요 (Overview)

성인형 스틸병(AOSD)은 원인 불명의 전신 자가염증 질환으로, 인구 100만 명당 약 1.5–2.2명의 발생률을 보이는 희귀 질환이다. 핵심 발병기전은 NLRP3 인플라마솜 활성화에 의한 IL-1β/IL-18 과잉 생산, 이로 인한 IL-6·TNF-α·IFN-γ 사이토카인 폭풍이며, 전신 거대포식세포 활성화가 특징이다. 임상적으로 고열(39°C 이상, 일중 1–2회), 일시적 연어색 피부발진, 관절염, 고 페리틴혈증(종종 10,000 ng/mL 이상)이 삼징으로 나타난다. 심각한 합병증으로 대식세포 활성화 증후군(MAS)이 약 10–15%에서 발생하며 생명을 위협한다. IL-1 차단제(아나킨라, 카나키누맙) 및 IL-6 차단제(토실리주맙)가 표준 생물학적 치료이다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| NLRP3 인플라마솜 | 병원체/위험신호 → Caspase-1 → pro-IL-1β 절단 | IL-1β 과잉 생산, 발열 |
| IL-18 축 | 대식세포 IL-18 분비 → NK세포·CD8⁺ T세포 활성 | MAS 유도, IFN-γ 증가 |
| IL-6 경로 | STAT3 활성화 → 급성기 반응 단백 합성 | CRP·페리틴 급상승, 빈혈 |
| 대식세포 활성화 | M1 분극화, 거식작용 항진 | 간비종대, 혈구감소증(MAS) |
| TNF-α 신호 | NF-κB → 다수 염증 유전자 발현 | 전신 염증, 관절 손상 |
| 페리틴 과잉 | 철 저장 단백 과잉 → 직접 면역 조절 이상 | 진단 바이오마커 역할 |
| NK세포 기능 장애 | IL-18에 의한 NK세포 탈진 | 세포독성 면역감시 소실 |

## 주요 약물 표적 (Drug Targets)

- **아나킨라 (Anakinra, IL-1Ra)**: IL-1 수용체 길항제; 급성기 조기 반응 우수, 반감기 짧아 매일 SC 투여
- **카나키누맙 (Canakinumab, 항-IL-1β 항체)**: 4주 1회 SC 투여; AOSD 및 MAS 예방에 효과적
- **토실리주맙 (Tocilizumab, 항-IL-6R 항체)**: IL-6 신호 차단; CRP·페리틴 정상화, 관절염 개선
- **코르티코스테로이드**: 1차 치료; 급성기 억제에 효과적이나 장기 사용 부작용 및 MAS 예방 효과 제한
- **토파시티닙 (Tofacitinib, JAK 억제제)**: 불응성 AOSD 및 MAS에서 IFN-γ 경로 차단 근거 증가

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [aosd_qsp_model.dot](aosd_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 161 노드 / 13 클러스터) |
| [aosd_qsp_model.svg](aosd_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [aosd_qsp_model.png](aosd_qsp_model.png) | PNG 이미지 (150 dpi) |
| [aosd_mrgsolve_model.R](aosd_mrgsolve_model.R) | mrgsolve ODE 모델 (약 21 구획 / 다수 치료 시나리오) |
| [aosd_shiny_app.R](aosd_shiny_app.R) | Shiny 대시보드 |
| [aosd_references.md](aosd_references.md) | 참고문헌 (약 67편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 아나킨라(SC 2구획 PK) + 카나키누맙(SC 2구획) + 토실리주맙(IV 2구획) + 코르티코스테로이드 + 토파시티닙 PK 구획 + IL-1β·IL-6·IL-18·IFN-γ·TNF-α·페리틴·CRP·활성 대식세포·NK세포·AOSD 활성도·MAS 위험 PD 구획
- **주요 치료 시나리오**: ① 무치료, ② 코르티코스테로이드 단독, ③ 아나킨라 + 스테로이드, ④ 카나키누맙, ⑤ 토실리주맙, ⑥ 토파시티닙(불응성), ⑦ 병용(스테로이드 + 생물학적 제제)
- **보정/근거**: SOBI ANAKIN 임상연구(아나킨라), Yamaguchi 분류기준, Fautrel 기준 파라미터 참조

## Shiny 대시보드 (Dashboard)

환자 프로파일(발열 패턴·피부발진·페리틴 초기치) 탭, 약물 PK 프로파일 탭, IL-1/IL-6/IL-18/IFN-γ 사이토카인 동태 탭, 임상 엔드포인트(발열·관절염 점수·페리틴) 탭, 치료 시나리오 비교 탭, MAS 위험 바이오마커 탭으로 구성된다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("aosd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("aosd_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg aosd_qsp_model.dot -o aosd_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [aosd_references.md](aosd_references.md) 참조 (약 67편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
