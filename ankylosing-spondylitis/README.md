# 강직성 척추염 (Ankylosing Spondylitis, AS) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 자가면역·류마티스

[![AS QSP Model](as_qsp_model.png)](as_qsp_model.svg)

## 개요 (Overview)

강직성 척추염(AS)은 척추와 천장관절을 주로 침범하는 만성 염증성 관절염으로, HLA-B27 양성률이 약 85–95%에 달하며 인구 0.1–0.5%의 유병률을 보인다. 발병기전의 핵심은 HLA-B27 관련 미접힘 단백질 반응(UPR)과 IL-23/IL-17 축의 부착부염(enthesitis)이다. 장내 미생물 불균형과 장 점막 염증이 IL-23 과잉 생산을 촉진하고, 이는 부착부 조직의 Th17세포를 활성화하여 IL-17A를 분비시킨다. 만성 염증은 이중적 결과를 초래하는데, 한편으로는 뼈 침식(손상), 다른 한편으로는 TNF에 의한 Wnt/BMP 신호 활성화로 인한 이소성 골 형성(신생골, syndesmophyte)이 진행된다. TNF 억제제와 IL-17A 억제제가 BASDAI/ASDAS 개선에 유효하며, JAK 억제제도 활용된다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| HLA-B27/UPR | 변이 HLA-B27 이중합체화 → ER 스트레스 → IL-23 | Th17 축 활성화 |
| IL-23/IL-17 축 | IL-23 → RORγt⁺ T세포 → IL-17A/F 분비 | 부착부 및 척추 염증 |
| TNF-α 신호 | TNF → NF-κB → 파골세포·연골세포 활성화 | 골 미란, 디스크 손상 |
| 신생골 형성 | Wnt/DKK1 균형 → BMP → 섬유아세포 골화 | 신데스모파이트, 척추 강직 |
| 장-척추 축 | 장내 미생물 이상 → 점막 IL-23 증가 | 척추관절염 유발 |
| RANKL/OPG | 염증 → RANKL↑ OPG↓ → 파골세포 증식 | 골미란, 골다공증 |
| CRP/질환 활성도 | 전신 염증 → ASDAS/BASDAI 상승 | 삶의 질 저하, 기능 손실 |

## 주요 약물 표적 (Drug Targets)

- **TNF 억제제 (아달리무맙, 에타너셉트, 세르톨리주맙 등)**: BASDAI 50% 개선(BASDAI50) 약 40–50%; 영상 염증 억제 우수, 신생골 억제 불완전
- **IL-17A 억제제 (세쿠키누맙, 익세키주맙)**: AS에서 TNF 억제제 동등 이상 효과; 장 크론병 동반 시 주의
- **JAK 억제제 (토파시티닙, 우파다시티닙)**: 경구 투여, TNF/IL-17 차단에 불응성 환자에서 유효
- **IL-23 억제제 (리산키주맙, 구셀쿠맙)**: 말초 관절염 우세형에서 효과; 축성 척추염 효과 연구 진행 중
- **NSAID**: 1차 증상 치료; 장기 복용 시 신생골 형성 억제 가능성

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [as_qsp_model.dot](as_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 191 노드 / 10 클러스터) |
| [as_qsp_model.svg](as_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [as_qsp_model.png](as_qsp_model.png) | PNG 이미지 (150 dpi) |
| [as_mrgsolve_model.R](as_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 다수 치료 시나리오) |
| [as_shiny_app.R](as_shiny_app.R) | Shiny 대시보드 |
| [as_references.md](as_references.md) | 참고문헌 (약 56편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 아달리무맙·에타너셉트(SC PK) + 세쿠키누맙(SC 2구획) + 토파시티닙·우파다시티닙(경구) + NSAID PK 구획 + TNF·IL-17A·IL-23·IL-6·CRP·RANKL·OPG·파골세포·골미란·신생골/mSASSS·BASDAI 복합 질환 활성도 PD 구획
- **주요 치료 시나리오**: ① 무치료, ② NSAID 단독, ③ 아달리무맙, ④ 에타너셉트, ⑤ 세쿠키누맙, ⑥ 토파시티닙, ⑦ TNF 실패 후 IL-17 전환
- **보정/근거**: MEASURE 임상시험(세쿠키누맙), ATLAS/ABILITY 아달리무맙/세르톨리주맙 데이터, mSASSS 진행 모델 파라미터 참조

## Shiny 대시보드 (Dashboard)

환자 프로파일(HLA-B27·BASDAI 초기치·영상 소견) 탭, 약물 PK 및 TNF/IL-17 억제율 탭, 질환 활성도(BASDAI/ASDAS) 변화 탭, 구조적 손상(mSASSS·신생골) 탭, 치료 시나리오 비교 탭, 바이오마커(CRP·RANKL·mSASSS) 탭으로 구성된다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("as_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("as_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg as_qsp_model.dot -o as_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [as_references.md](as_references.md) 참조 (약 56편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
