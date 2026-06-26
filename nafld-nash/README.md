# 비알코올 지방간/지방간염 (NAFLD/NASH) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도

[![NAFLD QSP Model](nafld_qsp_model.png)](nafld_qsp_model.svg)

## 개요 (Overview)

비알코올 지방간질환(NAFLD)은 알코올 섭취 없이 간세포 내 지방이 축적되는 질환으로, 전 세계 유병률은 약 25%에 달하고 제2형 당뇨병·비만과 밀접하게 연관된다. 단순 지방간(NAFL)에서 지방간염(NASH), 간섬유화, 간경변 및 간세포암(HCC)으로 진행할 수 있으며, NASH는 전 세계 간이식의 주요 원인 중 하나로 부상하고 있다. 발병기전은 '다중 히트(multiple-hit)' 모델로, 인슐린 저항성에 의한 간 지방 축적, 지방독성(lipotoxicity), Kupffer 세포 활성화 및 TNF-α/IL-6 매개 염증, TGF-β/HSC 활성화에 의한 섬유화 진행으로 이어진다. 2024년 FDA 승인을 받은 레스메티롬(THR-β 선택적 작용제)과 GLP-1 수용체 작용제가 현재 주요 치료제로 부상하고 있다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 인슐린 저항성 | IRS-1 인산화 장애 → 지방분해 억제 실패, 간 지방 축적 | 지방간(>5% 간세포 지방화) |
| 지방독성 | 유리지방산·DAG·세라미드 → ER 스트레스, 미토콘드리아 기능장애 | 간세포 손상, ALT 상승 |
| Kupffer 세포 활성화 | LPS/TLR4 → NF-κB → TNF-α, IL-6 분비 | 간 염증(NASH 등급) |
| 성상세포(HSC) 활성화 | TGF-β1 → SMAD2/3 → 콜라겐 I/III 합성 | 간섬유화(F0→F4) |
| 지질 대사 | THR-β → 간 지방산 산화↑, LDL 콜레스테롤↓ | 간 지방 감소 |
| 장-간 축 | 장내 미생물군 이상 → 단쇄지방산·LPS 유입 | 전신 염증 악화 |
| 지방세포 내분비 | 아디포넥틴↓, 렙틴↑, FGF21 변화 | 인슐린 저항성 악화 |

## 주요 약물 표적 (Drug Targets)

- **레스메티롬 (Resmetirom, THR-β 작용제)**: 간 선택적 갑상선 호르몬 수용체 β 활성화 → 지방산 산화 촉진, 간 지방 감소 (MAESTRO-NASH 3상)
- **세마글루타이드 (Semaglutide, GLP-1RA)**: 음식 섭취 감소, 체중 감량, 인슐린 감수성 개선 → 간 지방 및 NAS 감소
- **엠파글리플로진 (SGLT-2 억제제)**: 요당 배설로 에너지 음성 균형 → 간지방·체중 감소
- **오베티콜산 (Obeticholic acid, FXR 작용제)**: 담즙산 수용체 활성화 → 지질대사·간섬유화 개선
- **오를리스탯 (Orlistat, 지방분해효소 억제)**: 식이 지방 흡수 억제 → 체중 감량 보조

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [nafld_qsp_model.dot](nafld_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 140 노드 / 10 클러스터) |
| [nafld_qsp_model.svg](nafld_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [nafld_qsp_model.png](nafld_qsp_model.png) | PNG 이미지 (150 dpi) |
| [nafld_mrgsolve_model.R](nafld_mrgsolve_model.R) | mrgsolve ODE 모델 (약 21 구획 / 5개 치료 시나리오) |
| [nafld_shiny_app.R](nafld_shiny_app.R) | Shiny 대시보드 |
| [nafld_references.md](nafld_references.md) | 참고문헌 (약 47편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 레스메티롬 경구 3구획 PK(gut→중심→간), 오베티콜산 2구획, 세마글루타이드 SC 2구획, 엠파글리플로진 경구 2구획; 체중·체지방·아디포넥틴·인슐린 저항성·간 지방(LIVER_FAT)·Kupffer 세포·TNF-α·IL-6·TGF-β·HSC·콜라겐·ALT 구획
- **주요 치료 시나리오**: ① 위약(placebo) ② 레스메티롬 100mg QD ③ 오베티콜산 25mg QD ④ 세마글루타이드 2.4mg SC QW ⑤ 엠파글리플로진 10mg QD
- **보정/근거**: MAESTRO-NASH(레스메티롬), REGENERATE(오베티콜산), STEP-1(세마글루타이드) 임상시험 파라미터 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① **환자 프로파일** (BMI·당뇨 유무·기저 섬유화 단계 설정), ② **PK** (혈장/간 약물 농도-시간 곡선), ③ **PD 주요지표** (간 지방 분율, ALT, 인슐린 저항성 HOMA-IR), ④ **임상 엔드포인트** (NAS 점수, 섬유화 단계, 콜라겐 정량화), ⑤ **시나리오 비교** (5개 치료 전략 직접 비교), ⑥ **바이오마커** (TNF-α, TGF-β, 아디포넥틴, 체중 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("nafld_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("nafld_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg nafld_qsp_model.dot -o nafld_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [nafld_references.md](nafld_references.md) 참조 (약 47편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
