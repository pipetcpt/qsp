# 자가면역 용혈성 빈혈 (Autoimmune Hemolytic Anemia, AIHA) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 혈액

[![AIHA QSP Model](aiha_qsp_model.png)](aiha_qsp_model.svg)

## 개요 (Overview)
자가면역 용혈성 빈혈(AIHA)은 자가항체가 적혈구 표면 항원에 결합하여 조기 파괴를 유발하는 희귀 혈액질환으로, 인구 10만 명당 약 1~3명의 유병률을 보입니다. 온형(warm) AIHA는 IgG 자가항체(주로 Rh 항원 대상)가 비장에서 식세포 제거를 촉진하며, 한랭(cold) AIHA는 IgM 항체가 보체 C1q를 활성화하여 혈관 내 용혈을 유발합니다. 핵심 치료 표적은 자가항체 생성 B세포·형질세포, C1s/C3b 보체 경로, 비장 대식세포의 Fc수용체입니다. 최근 승인된 수티림맙(항C1s)은 한랭응집소병(CAD)의 혈관 내 용혈을 직접 차단하여 큰 관심을 받고 있습니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| B세포 활성화 및 자가항체 생성 | BCR 신호, TLR7/9, T-B 상호작용 | 항-RBC IgG/IgM 자가항체 생성 |
| 비장 식세포 제거 (세포외 용혈) | FcγRI/III, C3b 옵소닌화 | 비장비대, 빈혈 악화 |
| 보체 활성화 (혈관 내 용혈) | C1q/C1s 활성화 → C3b/MAC 형성 | 헤모글로빈뇨, 급성 용혈 위기 |
| EPO-망상적혈구 보상 | EPO 수용체, BFU-E/CFU-E 증식 | 망상적혈구증가증, 골수 과형성 |
| 조직 철 재활용 장애 | 헤모글로빈 · 합토글로빈 고갈 | 혈청 LDH 상승, 빌리루빈 증가 |
| Syk/BTK 신호 (FcR 하류) | 비장 대식세포 활성화 | 포스타마티닙(Syk 억제제) 표적 |

## 주요 약물 표적 (Drug Targets)
- **코르티코스테로이드 (프레드니솔론, 덱사메타손)**: 스테로이드 수용체 → T세포 억제, 대식세포 Fc수용체 발현 감소 — 1차 표준 치료
- **리툭시맙 (항-CD20)**: B세포 고갈 → 자가항체 생성 억제 — 2차 또는 1차 병합 치료
- **수티림맙 (항-C1s)**: 보체 C1s 차단 → 혈관 내 용혈 억제 — CAD 승인 치료제 (CADENZA trial)
- **포스타마티닙 (Syk 억제제)**: FcR 하류 신호 차단 → 비장 식세포 활성 감소 — 난치성 온형 AIHA
- **MMF (미코페놀레이트 모페틸)**: 림프구 증식 억제 → 자가항체 감소 — 유지 치료
- **IVIG**: FcRn 포화 및 Fc수용체 차단 → 급성 구조 치료

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [aiha_qsp_model.dot](aiha_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 219 노드 / 15 클러스터) |
| [aiha_qsp_model.svg](aiha_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [aiha_qsp_model.png](aiha_qsp_model.png) | PNG 이미지 (150 dpi) |
| [aiha_mrgsolve_model.R](aiha_mrgsolve_model.R) | mrgsolve ODE 모델 (약 31 구획 / 8개 치료 시나리오) |
| [aiha_shiny_app.R](aiha_shiny_app.R) | Shiny 대시보드 |
| [aiha_references.md](aiha_references.md) | 참고문헌 (약 76편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 면역 구획(B세포·형질세포·IgG/IgM 자가항체·보체 C1s/C3b/MAC), RBC 구획(RBC·망상적혈구·EPO), 바이오마커(LDH·합토글로빈·빌리루빈), 약물 PK 구획(리툭시맙 TMDD, 수티림맙, 프레드니솔론, 포스타마티닙/R406, 덱사메타손, MMF/MPA, IVIG)
- **주요 치료 시나리오**: ① 비치료(자연경과), ② 프레드니솔론 1 mg/kg/일→감량, ③ 프레드니솔론+리툭시맙, ④ 덱사메타손 펄스+리툭시맙, ⑤ 포스타마티닙 150 mg BID(난치성), ⑥ 프레드니솔론+MMF, ⑦ IVIG 급성 구조, ⑧ 수티림맙 6.5 g q2w(CAD/CADENZA)
- **보정/근거**: CADENZA 시험(수티림맙), FACILITATE 시험(포스타마티닙), LAZAR 등 리툭시맙 메타분석 데이터를 참고하여 파라미터 보정

## Shiny 대시보드 (Dashboard)
환자 프로파일(온형/한랭 아형 선택, 체중, 기저 Hb), 약동학(각 약물 혈장 농도-시간 곡선), 혈액학 지표(Hb·RBC·망상적혈구 변화), 보체·면역 바이오마커(C3b·LDH·빌리루빈·합토글로빈), 치료 시나리오 비교(8개 시나리오 overlay), 임상 엔드포인트(빈혈 관해율·용혈 지표) 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("aiha_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("aiha_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg aiha_qsp_model.dot -o aiha_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [aiha_references.md](aiha_references.md) 참조 (약 76편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
