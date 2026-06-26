# 항인지질항체 증후군 (Antiphospholipid Syndrome, APS) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 자가면역·류마티스

[![APS QSP Model](aps_qsp_model.png)](aps_qsp_model.svg)

## 개요 (Overview)

항인지질항체 증후군(APS)은 항인지질항체(aPL: 루푸스 항응고인자·항카르디올리핀 항체·항β₂-당단백 I 항체)에 의해 동맥·정맥 혈전증 및 임신 합병증이 반복적으로 발생하는 자가면역 혈전증 질환이다. 전신 홍반 루푸스(SLE) 환자의 약 30%에서 동반되며(이차성 APS), 독립적으로도 발생한다(원발성 APS). aPL은 β₂-GPI에 결합하여 내피세포·혈소판·단핵구를 활성화하고, 보체 경로(C5a)와 조직인자(TF) 발현을 촉진하여 혈전 형성을 유도한다. 삼중 양성(triple positive: 세 가지 aPL 항체 모두 양성) 환자에서 혈전 재발 위험이 가장 높다. 항응고 치료(와파린/LMWH)와 히드록시클로로퀸이 핵심 치료이며, 난치성에는 리툭시맙·벨리무맙이 시도된다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| aPL-β₂GPI 복합체 | aPL이 β₂-GPI에 결합 → 내피세포·혈소판 활성화 | 혈전 선행 조건 형성 |
| 내피세포 활성화 | NF-κB → TF, ICAM-1, VCAM-1 발현 | 동맥·정맥 혈전증 |
| 보체 활성화 | C5a → 호중구 활성화, NETs 형성 | 혈전·임신 합병증 |
| 혈소판 활성화 | aPL → GpIbα 수용체 → 혈소판 응집 | 동맥 혈전증, 혈소판 감소증 |
| 응고계 과잉 | 트롬빈 생성 증가, 프로테인 C/S 억제 | 혈전증 위험 증폭 |
| mTOR/신장 | 내피세포 mTOR 활성화 → 신 미세혈전증 | APS 신증(CAPS 위험) |
| 태반 보체 | 태반 영양막 세포 보체 침착 → 세포 손상 | 임신 소실, 전자간증 |

## 주요 약물 표적 (Drug Targets)

- **와파린 (Warfarin, 비타민 K 길항제)**: 정맥 혈전 APS의 표준 항응고; INR 목표 2.0–3.0(고위험 시 3.0–4.0)
- **저분자량 헤파린 (LMWH)**: 임신 중 APS의 핵심 치료; 태반 보호 효과 포함
- **히드록시클로로퀸 (HCQ)**: TLR 억제 → aPL 생산 감소, 혈소판 응집 억제; SLE-APS에서 혈전 예방
- **리바록사반 (Rivaroxaban, 직접 Xa 인자 억제제)**: 삼중 양성 환자에서 와파린보다 열등; 단순 정맥 혈전에서만 대안
- **리툭시맙 (Rituximab, 항-CD20 항체)**: 난치성 APS, CAPS(재앙성 APS)에서 B세포 고갈로 aPL 감소
- **저용량 아스피린**: 동맥 혈전 고위험 및 임신 APS 예방 병용

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [aps_qsp_model.dot](aps_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 186 노드 / 13 클러스터) |
| [aps_qsp_model.svg](aps_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [aps_qsp_model.png](aps_qsp_model.png) | PNG 이미지 (150 dpi) |
| [aps_mrgsolve_model.R](aps_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 다수 치료 시나리오) |
| [aps_shiny_app.R](aps_shiny_app.R) | Shiny 대시보드 |
| [aps_references.md](aps_references.md) | 참고문헌 (약 58편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 와파린(경구 2구획 + 효과 구획) + LMWH + HCQ(2구획) + 리바록사반(2구획) + 아스피린 PK 구획 + 리툭시맙(2구획) + aPL IgG·B세포·보체 C5a·내피세포 TF·혈소판 활성·트롬빈 생성·심부정맥혈전 위험·임신 생존율·mTOR 신장·INR PD 구획
- **주요 치료 시나리오**: ① 무치료, ② 와파린 단독, ③ LMWH 단독(임신), ④ 와파린 + HCQ, ⑤ 리바록사반(저위험), ⑥ 아스피린 + HCQ, ⑦ 리툭시맙 + 항응고(불응성 APS)
- **보정/근거**: RAPS 임상시험(리바록사반 vs. 와파린), Crowther et al. NEJM 와파린 INR 목표치, Beppu et al. CAPS 레지스트리 파라미터 참조

## Shiny 대시보드 (Dashboard)

환자 프로파일(aPL 프로파일·혈전 기왕력·임신 여부) 탭, 항응고제 PK 및 INR/항Xa 수준 탭, 혈전 형성 위험 동태 탭, 임상 엔드포인트(혈전 재발·임신 결과) 탭, 치료 시나리오 비교 탭, 바이오마커(aPL 역가·보체·TF) 탭으로 구성된다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("aps_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("aps_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg aps_qsp_model.dot -o aps_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [aps_references.md](aps_references.md) 참조 (약 58편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
