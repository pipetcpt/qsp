# 자가면역 췌장염 (Autoimmune Pancreatitis, AIP) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도

[![AIP QSP Model](aip_qsp_model.png)](aip_qsp_model.svg)

## 개요 (Overview)
자가면역 췌장염(AIP)은 IgG4 관련 질환(IgG4-RD)의 췌장 표현형으로, 림프형질세포 침윤과 섬유화를 특징으로 합니다. 1형 AIP(IgG4-관련)는 혈청 IgG4 상승과 다발장기 침범을 보이며, 2형 AIP는 과립구 상피 병변을 특징으로 하고 IgG4와는 무관합니다. 발병기전의 핵심은 Th2/Treg 분화 편향, IgG4+ 형질세포 과증식, 췌장 성상세포 활성화에 의한 섬유화입니다. 스테로이드가 관해 유도의 표준이며 재발 시 리툭시맙이 효과적입니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| Th2/Treg 편향 | IL-4, IL-10, TGF-β 과생성 | IgG4 class switching 촉진 |
| IgG4+ 형질세포 침윤 | Tfh2-B세포 상호작용, 생발 중심 반응 | 혈청 IgG4 상승 (>135 mg/dL 진단 기준) |
| 췌장 성상세포 활성화 | TGF-β, PDGF → α-SMA 발현 | 췌장 섬유화·위축 |
| 담관 침범 | 림프형질세포성 경화성 담관염 | 폐쇄성 황달, ALP/GGT 상승 |
| 외분비 기능 부전 | 선포세포 손상, 섬유화 | 소화효소 분비 감소, 지방변 |
| 내분비 기능 부전 | β세포 손실, 인슐린 분비 저하 | 췌장원성 당뇨병 |
| IgG4 면역복합체 | 활성화 Fcγ수용체 저해(블로킹 항체 가설) | 염증 지속/조직 손상 |

## 주요 약물 표적 (Drug Targets)
- **프레드니솔론 (0.6 mg/kg/일 유도 → 감량)**: GR 매개 Th1/Th2 균형 조절, 림프형질세포 감소 — 1차 관해 유도 표준
- **리툭시맙 (항-CD20)**: B세포·IgG4+ 형질아세포 고갈 → IgG4 급격히 감소 — 재발성/스테로이드 의존성 AIP 2차 치료
- **아자티오프린 → 6-TGN**: 림프구 증식 억제 → 스테로이드 절약 유지 치료
- **MMF/MPA**: IMPDH 억제 → Th2 B세포 증식 억제 — 아자티오프린 대체

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [aip_qsp_model.dot](aip_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 167 노드 / 11 클러스터) |
| [aip_qsp_model.svg](aip_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [aip_qsp_model.png](aip_qsp_model.png) | PNG 이미지 (150 dpi) |
| [aip_mrgsolve_model.R](aip_mrgsolve_model.R) | mrgsolve ODE 모델 (약 24 구획 / 7개 치료 시나리오) |
| [aip_shiny_app.R](aip_shiny_app.R) | Shiny 대시보드 |
| [aip_references.md](aip_references.md) | 참고문헌 (약 44편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 약물 PK(프레드니솔론 3구획, 리툭시맙 2구획, 아자티오프린→6-MP→6-TGN), 면역 세포(naïve CD4, Th2, Treg, naïve B, GC B, IgG4+ 형질아세포, 혈청 IgG4), 섬유화 단계(정지·활성화 성상세포, 콜라겐), 장기 기능(외분비·β세포), 사이토카인(IL-4, IL-10, TGF-β, TNF-α)
- **주요 치료 시나리오**: ① 프레드니솔론 유도+감량, ② 프레드니솔론 유지(5 mg/일×3년), ③ 리툭시맙 1000 mg×2, ④ 프레드니솔론+아자티오프린, ⑤ RTX 유도+재발 시 재치료, ⑥ 프레드니솔론+아자티오프린 병합, ⑦ 비치료(자연경과)
- **보정/근거**: Kamisawa, Hart, Khosroshahi 등 IgG4-RD 스테로이드 반응 연구 및 리툭시맙 임상 데이터 참고

## Shiny 대시보드 (Dashboard)
환자 프로파일(AIP 아형·IgG4 기저치·당뇨 동반 여부), 약동학(스테로이드·RTX 혈장 농도), PD 면역 지표(IgG4·형질아세포·Th2/Treg), 장기 기능 탭(외분비·내분비·섬유화 지수), 치료 시나리오 비교, 재발 예측 바이오마커 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("aip_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("aip_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg aip_qsp_model.dot -o aip_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [aip_references.md](aip_references.md) 참조 (약 44편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
