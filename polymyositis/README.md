# 다발성 근염 (Polymyositis, PM) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 자가면역·류마티스

[![PM QSP Model](pm_qsp_model.png)](pm_qsp_model.svg)

## 개요 (Overview)
다발성 근염(PM)은 CD8+ 세포독성 T세포가 근섬유를 직접 공격하는 자가면역 근염으로, 인구 10만 명당 약 1~5명에서 발생하며 여성에서 더 흔하다. MHC-I 과발현 및 CD8+ T세포의 perforin·granzyme B 매개 근섬유 침습이 주된 조직 손상 기전이며, 항Jo-1 항체(항합성효소 항체) 등 근염 특이 자가항체가 진단 및 예후 인자로 활용된다. 근력 저하(proximal myopathy), 혈청 CK 상승, 근전도 이상, 근생검에서 CD8 T세포 침윤이 진단 기준을 이룬다. 고용량 스테로이드가 1차 치료이며 스테로이드 절약 면역억제제(MTX, AZA), 의례의 경우 리툭시맙이나 JAK 억제제가 추가된다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| MHC-I 발현 경로 | IFN-β/α 자극 → STAT1 활성화 → MHC-I 과발현 | CD8 T세포 근섬유 인식 가능 |
| CD8 세포독성 경로 | Perforin·Granzyme B, TRAIL·FasL | 근섬유 괴사·CK 방출 |
| Th1 염증 경로 | IFN-γ, TNF-α, IL-6 | 근육 내 염증 증폭 |
| B세포·자가항체 경로 | 항Jo-1·항Mi-2 등 MSA, 형질세포 | 항합성효소 증후군, ILD 위험 |
| 근육 재생 장애 경로 | 만성 염증 → 위성세포 고갈 | 비가역적 근력 저하 |
| ILD 동반 경로 | IFN, TGF-β, IL-6 → 폐 섬유화 | 항Jo-1 양성에서 간질성 폐질환 |

## 주요 약물 표적 (Drug Targets)
- **프레드니솔론**: 광범위 항염·면역억제 — NF-κB·AP-1 억제, 1차 치료
- **메토트렉세이트 (MTX)**: 엽산 대사 길항, 다핵형 PG 축적 — 스테로이드 절약제
- **아자티오프린 (AZA)/6-TGN**: 퓨린 합성 억제, 관해 유지
- **리툭시맙 (RTX)**: 항CD20 단클론항체 → B세포·형질세포 고갈, 내화 MSA(항Jo-1 등) 양성 PM
- **JAK 억제제**: JAK1/2 억제 → IFN 신호 차단 → MHC-I 발현 억제

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [pm_qsp_model.dot](pm_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 193 노드 / 12 클러스터) |
| [pm_qsp_model.svg](pm_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [pm_qsp_model.png](pm_qsp_model.png) | PNG 이미지 (150 dpi) |
| [pm_mrgsolve_model.R](pm_mrgsolve_model.R) | mrgsolve ODE 모델 (약 28 구획 / 6 치료 시나리오) |
| [pm_shiny_app.R](pm_shiny_app.R) | Shiny 대시보드 |
| [pm_references.md](pm_references.md) | 참고문헌 (약 53편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 프레드니솔론·MTX·6-TGN·RTX(2구획+CD20결합)·IVIG·JAK 억제제 PK 구획, CD8N/CD8E·CD4Th1·B세포·형질세포·자가항체·IFN-γ·TNF-α·IL-6·MHC-I·근육 염증·CK·MMT8(근력) PD 구획
- **주요 치료 시나리오**: ① 무치료, ② 프레드니솔론 단독, ③ 프레드니솔론+MTX, ④ 프레드니솔론+AZA, ⑤ 리툭시맙+프레드니솔론, ⑥ JAK억제제+프레드니솔론
- **보정/근거**: RIM 시험(리툭시맙 PM/DM) 및 IMACS 국제 네트워크 임상 데이터를 기반으로 MMT8 근력 점수 및 CK 정상화 시간 보정

## Shiny 대시보드 (Dashboard)
환자 프로파일(MSA 상태, 근력 점수, ILD 동반) · 약물 PK(RTX CD20 포화 포함) · 면역 PD(CD8·B세포·사이토카인) · 근육 임상 엔드포인트(MMT8·CK) · 치료 시나리오 비교(관해율·스테로이드 용량) · 바이오마커(IFN-γ·자가항체) 탭으로 구성

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("pm_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("pm_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg pm_qsp_model.dot -o pm_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [pm_references.md](pm_references.md) 참조 (약 53편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
