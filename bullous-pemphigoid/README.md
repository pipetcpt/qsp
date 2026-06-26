# 수포성 유천포창 (Bullous Pemphigoid, BP) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 피부

[![BP QSP Model](bp_qsp_model.png)](bp_qsp_model.svg)

## 개요 (Overview)
수포성 유천포창(BP)은 표피-진피 접합부(DEJ)의 BP180(collagen XVII)과 BP230에 대한 자가항체(IgG, IgE)가 표피하 수포를 유발하는 가장 흔한 자가면역 수포 질환입니다. 주로 70세 이상 고령에서 발생하며(유병률 인구 10만 명당 약 20~40명), Th2 편향 면역과 호산구 침윤, 비만세포 활성화, 보체 C5a 매개 호중구 모집이 조직 손상의 핵심입니다. 표준 치료는 전신 코르티코스테로이드이며, 노령 부작용 감소를 위해 독시사이클린+국소 스테로이드, 두필루맙, 오말리주맙, 리툭시맙이 사용됩니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| Th2 편향 및 B세포 활성화 | IL-4/IL-13 → class switching → IgE/IgG 생성 | 항-BP180/BP230 자가항체 생성 |
| 항-BP180 IgG 보체 활성화 | C1q 활성화 → C3b → C5a 생성 | 호중구·호산구 DEJ 모집 |
| 항-BP180 IgE 비만세포 활성화 | FcεRI 교차결합 → 탈과립 → 히스타민 | 소양증, 두드러기양 전구 병변 |
| 호산구 침윤 | IL-5, eotaxin → 호산구 MBP/ECP 방출 | DEJ 추가 손상, 수포 확대 |
| 단백분해효소 매개 DEJ 손상 | MMP-9, NE → BP180 절단 | 표피하 수포 형성 |
| 장기지속 형질세포 | 골수 LLPC → 지속적 자가항체 생성 | 재발 및 불완전 관해 |

## 주요 약물 표적 (Drug Targets)
- **프레드니솔론 (고용량 taper)**: 광범위 항염 → 급속 수포 억제 — 표준 1차 치료 (BLISTER 시험)
- **저용량 프레드니솔론+독시사이클린**: MMP 억제, 항염 — 고령 환자 코르티코스테로이드 절약
- **두필루맙 (항-IL-4Rα)**: IL-4/IL-13 차단 → Th2·IgE 감소 — 생물학제제 옵션 (여러 사례 보고)
- **오말리주맙 (항-IgE)**: 유리 IgE 중화 → 비만세포·호산구 활성화 억제 — IgE 고가 BP
- **리툭시맙 (항-CD20)**: B세포·형질세포 고갈 → 자가항체 감소 — 난치성/재발성 BP

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [bp_qsp_model.dot](bp_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 162 노드 / 13 클러스터) |
| [bp_qsp_model.svg](bp_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [bp_qsp_model.png](bp_qsp_model.png) | PNG 이미지 (150 dpi) |
| [bp_mrgsolve_model.R](bp_mrgsolve_model.R) | mrgsolve ODE 모델 (약 23 구획 / 6개 치료 시나리오) |
| [bp_shiny_app.R](bp_shiny_app.R) | Shiny 대시보드 |
| [bp_references.md](bp_references.md) | 참고문헌 (약 65편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 약물 PK(프레드니솔론 3구획, 두필루맙 SC 2구획, 오말리주맙 SC 2구획, 리툭시맙, 독시사이클린 2구획), 면역 PD(B세포 naïve/활성/기억/SLPC/LLPC, Th2, 항-BP180 IgG/IgE, 혈액 호산구, 피부 호산구, 비만세포, 보체 C5a, DEJ 완전성)
- **주요 치료 시나리오**: ① 비치료, ② 고용량 프레드니솔론 감량, ③ 저용량 프레드니솔론+독시사이클린, ④ 두필루맙+저용량 프레드니솔론, ⑤ 오말리주맙+저용량 프레드니솔론, ⑥ 리툭시맙 1g×2+단기 스테로이드
- **보정/근거**: BLISTER 시험(독시사이클린), MVP(베르탈리주맙) 등 주요 BP 임상시험 및 항체 역가·BPDAI 점수 데이터 참고

## Shiny 대시보드 (Dashboard)
환자 프로파일(연령·기저 IgE·호산구·BPDAI), 약동학(각 약물 혈장 농도), 면역 PD 탭(B세포 분획·IgG/IgE 항체 역가·호산구), 피부 손상 탭(DEJ 완전성·수포 형성·BPDAI), 치료 시나리오 비교(6개), 재발 예측 바이오마커 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("bp_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("bp_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg bp_qsp_model.dot -o bp_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [bp_references.md](bp_references.md) 참조 (약 65편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
