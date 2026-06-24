# 게실병 (Diverticular Disease, DIV) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도

[![DIV QSP Model](div_qsp_model.png)](div_qsp_model.svg)

## 개요 (Overview)

게실병(Diverticular Disease)은 대장 벽에 소낭(게실)이 형성되는 질환으로, 서구화된 식이 패턴과 밀접하게 연관됩니다. 60세 이상 서양인의 약 60%에서 게실이 발견되며, 이 중 약 25%가 게실염으로 진행합니다. 핵심 발병기전은 저섬유식에 의한 장내압 상승·장벽 구조 약화, 장내 미생물 불균형(dysbiosis), 점막 방어력 감소로, 게실 내 세균 과증식이 NF-κB 경로를 통해 염증을 유발합니다. 치료는 고섬유식·리팍시민(cyclic)·메살라민(5-ASA)·급성기 항생제로 이루어집니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 장내압 상승 | 저섬유 → 장 내강 용적 감소 → 분절성 수축 과항진 | 게실 형성(점막 탈출) |
| 장벽 구조 약화 | 콜라겐 이상, 기저 근육층 결손 | 게실 수 증가 |
| 장내 미생물 불균형 | 보호균 감소, 병원성균 과증식 → LPS 생성 | 점막 방어력 저하 |
| NF-κB·염증 경로 | LPS → TLR4 → NF-κB → TNF-α, IL-6, IL-1β | 게실주위 염증, CRP 상승 |
| 점막 손상 | 부티레이트 감소 → 점막 상피 재생 저하 | 만성 염증 지속 |
| 내장 과민성 | 만성 염증 → 내장 감각 신경 과민 | 복통, IBS 유사 증상 |
| 급성 합병증 | 게실 천공·농양·누공·폐쇄 | 입원·수술 필요 |

## 주요 약물 표적 (Drug Targets)

- **리팍시민(Rifaximin)**: 장관 선택적 항생제 → 병원성균 억제, LPS 생성 감소
- **메살라민(5-ASA)**: NF-κB 억제 → 점막 염증 완화, 재발 예방
- **고섬유식/차전자피(Psyllium)**: 장내압 저하, 보호균 증식 촉진
- **시프로플록사신 + 메트로니다졸**: 급성 게실염 항균 치료
- **프로바이오틱스**: 장내 미생물 생태계 회복

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [div_qsp_model.dot](div_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 152 노드 / 10 클러스터) |
| [div_qsp_model.svg](div_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [div_qsp_model.png](div_qsp_model.png) | PNG 이미지 (150 dpi) |
| [div_mrgsolve_model.R](div_mrgsolve_model.R) | mrgsolve ODE 모델 (약 20 구획 / 6 치료 시나리오) |
| [div_shiny_app.R](div_shiny_app.R) | Shiny 대시보드 |
| [div_references.md](div_references.md) | 참고문헌 (약 61편, DOI 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 장내 섬유 농도, 장내압, 점막 방어력, 보호균·병원성균 밀도, LPS, NF-κB, TNF-α/IL-6/IL-1β, 호중구, CRP, 게실 수, 만성 염증, 약물(리팍시민·메살라민·시프로·메트로) PK 구획, 내장 과민성, 콜라겐
- **주요 치료 시나리오**: ① 자연 경과(서양식 식이), ② 고섬유 식이(차전자피 +15 g/일), ③ 리팍시민 cyclic(400 mg TID × 7일/월), ④ 메살라민 유지요법(1.6 g/일), ⑤ 급성 게실염 — 시프로플록사신 + 메트로니다졸, ⑥ 복합 요법(고섬유 + 리팍시민 + 메살라민)
- **보정/근거**: PREVENT 리팍시민 연구, DIVER 메살라민 시험, Strate 등의 역학 데이터 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① 환자 프로파일(식이 섬유 섭취·BMI·기저 게실 수 설정), ② PK 탭(항생제·메살라민 혈중 농도), ③ 미생물·염증 PD 탭(보호균/병원성균·NF-κB·사이토카인 추이), ④ 임상 엔드포인트(게실 수·CRP·게실염 발생 위험), ⑤ 시나리오 비교(6개 치료 전략 동시 비교), ⑥ 바이오마커(LPS·내장 과민성·점막 방어력 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("div_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("div_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg div_qsp_model.dot -o div_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [div_references.md](div_references.md) 참조 (약 61편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
