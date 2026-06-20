# 결절성 다발동맥염 (PAN) (Polyarteritis Nodosa, PAN) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 혈관염

[![PAN QSP Model](pan_qsp_model.png)](pan_qsp_model.svg)

## 개요 (Overview)
결절성 다발동맥염(PAN)은 중형 동맥의 괴사성 혈관염으로, ANCA 음성이며 신사구체·폐 모세혈관을 침범하지 않는 점이 특징적이다. 유병률은 인구 100만 명당 약 4~9명 수준이며, B형 간염 바이러스(HBV) 연관 PAN과 특발성 PAN으로 구분된다. 핵심 발병기전은 면역복합체(특히 HBV surface Ag)의 혈관벽 침착 → 보체 활성화 → 중성구·대식세포 침윤 → 혈관 괴사 → 동맥류·경색으로 이어지는 과정이다. 고용량 스테로이드 단독(비HBV) 또는 스테로이드+항바이러스+혈장교환(HBV 연관)이 표준치료이며, 재발성·난치성에는 시클로포스파미드가 추가된다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 면역복합체 형성 경로 | HBsAg-항체 복합체, 보체 C3·C4 소모 | 혈관벽 침착·보체 활성화 |
| 보체·중성구 경로 | C5a 주화성, 중성구 활성화·산화 폭발 | 혈관 내피 손상·괴사 |
| 대식세포 활성화 경로 | TNF-α·IL-6·IL-1β, 산화질소 | 염증 증폭·조직 파괴 |
| T세포 활성화 경로 | CD4 Th1, 세포성 면역 | 자가면역 혈관염 유지 |
| 혈관 괴사·섬유화 경로 | 혈관 중막 괴사, 섬유소 침착 | 동맥류 형성·혈관 폐쇄 |
| 장기 손상 경로 | 신장(신부전), 말초신경(신경병증), 피부(결절) | 다장기 기능 저하 |

## 주요 약물 표적 (Drug Targets)
- **고용량 프레드니솔론**: 항염·면역억제 — NF-κB 억제, 혈관 염증 급속 억제
- **시클로포스파미드 (CYC)**: 알킬화 면역억제 — T·B세포 증식 억제, 중증·재발 PAN
- **아자티오프린 (AZA)**: 퓨린 대사 억제 — 관해 유지
- **항바이러스제 (HBV 연관)**: 엔테카비르·테노포비르 — HBsAg 제거로 근본 원인 치료
- **혈장교환**: 순환 면역복합체 제거 — HBV 연관 PAN 급성기 보조

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [pan_qsp_model.dot](pan_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 237 노드 / 13 클러스터) |
| [pan_qsp_model.svg](pan_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [pan_qsp_model.png](pan_qsp_model.png) | PNG 이미지 (150 dpi) |
| [pan_mrgsolve_model.R](pan_mrgsolve_model.R) | mrgsolve ODE 모델 (약 19 구획 / 6 치료 시나리오) |
| [pan_shiny_app.R](pan_shiny_app.R) | Shiny 대시보드 |
| [pan_references.md](pan_references.md) | 참고문헌 (약 44편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 프레드니솔론·시클로포스파미드·활성대사체·아자티오프린 PK 구획, 면역복합체·보체·중성구·대식세포·B세포·T세포·사이토카인(TNF-α)·혈관 염증·섬유화·동맥류·장기 손상(신장·신경·범발성) PD 구획
- **주요 치료 시나리오**: ① 무치료(자연경과), ② 스테로이드 단독(유도+유지), ③ 스테로이드+시클로포스파미드, ④ 스테로이드+아자티오프린, ⑤ 스테로이드 단독(감량), ⑥ HBV 연관 PAN(스테로이드+CYC+AZA)
- **보정/근거**: French Vasculitis Study Group(FVSG) 코호트 데이터 및 Five Factor Score(FFS) 예후 지표를 참고하여 관해율·재발률 파라미터 보정

## Shiny 대시보드 (Dashboard)
환자 프로파일(HBV 상태, FFS 점수, 장기 침범) · 약물 PK 시뮬레이션 · 면역 PD(면역복합체·사이토카인) · 혈관 손상·장기 기능 임상 엔드포인트(BVAS) · 치료 시나리오 비교(관해율·재발률) · 바이오마커(CRP, 보체, HBsAg) 탭으로 구성

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("pan_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("pan_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg pan_qsp_model.dot -o pan_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [pan_references.md](pan_references.md) 참조 (약 44편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
