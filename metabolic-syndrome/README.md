# 대사 증후군 (Metabolic Syndrome, MS) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![MS QSP Model](ms_qsp_model.png)](ms_qsp_model.svg)

## 개요 (Overview)
대사 증후군은 인슐린 저항성을 공통 기반으로 복부 비만·고혈당·이상지질혈증·고혈압이 군집하는 대사 이상 복합체입니다. 전 세계 성인 유병률은 약 25%(일부 국가 30% 이상)이며, 제2형 당뇨병 발생 위험을 5배, 심혈관 질환 위험을 2~3배 높입니다. 핵심 발병기전은 내장 지방 과잉 축적 → 유리지방산·아디포카인 불균형(렙틴 증가/아디포넥틴 감소) → 간·골격근·지방조직 인슐린 저항성 → 보상적 고인슐린혈증 → 다장기 대사 이상의 악순환입니다. 치료는 생활습관 교정(식이·운동)이 기본이며, 메트포르민·GLP-1 수용체 작용제·SGLT-2 억제제·스타틴·ARB를 조합하여 각 성분을 표적합니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 내장지방 축적 | 지방세포 비대·과형성, 지방독성(FFA 과잉) | 인슐린 저항성 유도 |
| 인슐린 저항성 | IRS-1 세린 인산화·PI3K/Akt 신호 감소, JNK 활성화 | 고혈당·고인슐린혈증 |
| 지방간 및 이상지질혈증 | VLDL 과생산, LDL 소입자화, HDL 감소, TG 증가 | 죽상경화 촉진 |
| 저등급 만성 염증 | 지방조직 대식세포 침윤, TNF-α·IL-6·CRP 상승 | 인슐린 저항성 악화 |
| RAAS/교감신경 활성화 | 안지오텐신 II 증가·교감신경 과활성 | 혈압 상승·신장 손상 |
| AMPK 경로 억제 | 에너지 센서 AMPK 활성 저하 | 지방산 산화 감소·포도당 신합성 증가 |
| 아디포카인 불균형 | 렙틴 저항성, 아디포넥틴 감소 | 식욕·에너지 조절 이상 |

## 주요 약물 표적 (Drug Targets)
- **메트포르민**: 간 포도당 신합성 억제(AMPK 활성화), 인슐린 감수성 개선 (1차)
- **GLP-1 수용체 작용제** (세마글루타이드·리라글루타이드): 인슐린 분비 촉진·글루카곤 억제·체중 감소·심혈관 보호
- **SGLT-2 억제제** (엠파글리플로진·다파글리플로진): 신장 포도당 재흡수 억제, 체중·혈압 감소, 심신 보호
- **스타틴** (아토르바스타틴·로수바스타틴): HMG-CoA 환원효소 억제 → LDL-C 감소, 심혈관 위험 감소
- **ARB** (발사르탄·이르베사르탄): 안지오텐신 II 수용체 차단, 혈압·인슐린 저항성 개선
- **피오글리타존** (TZD/PPARγ 작용제): 말초 인슐린 감수성 개선, 아디포넥틴 증가

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [ms_qsp_model.dot](ms_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 174 노드 / 16 클러스터) |
| [ms_qsp_model.svg](ms_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [ms_qsp_model.png](ms_qsp_model.png) | PNG 이미지 (150 dpi) |
| [ms_mrgsolve_model.R](ms_mrgsolve_model.R) | mrgsolve ODE 모델 (약 32 구획 / 6개 치료 시나리오) |
| [ms_shiny_app.R](ms_shiny_app.R) | Shiny 대시보드 |
| [ms_references.md](ms_references.md) | 참고문헌 (약 50편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 약물 PK 구획(메트포르민 3구획, GLP-1 작용제 2구획, SGLT-2 억제제 2구획, 스타틴 2구획, ARB 2구획) + 질환 PD 구획(혈당·인슐린·GLP-1·내장지방·피하지방·렙틴·아디포넥틴·VLDL-C·LDL-C·HDL-C·TG·TNF-α·IL-6·IL-1β·CRP·안지오텐신 II·MAP·AMPK)
- **주요 치료 시나리오**: ① 생활습관(무약물), ② 메트포르민 단독, ③ GLP-1 수용체 작용제 단독, ④ 메트포르민 + GLP-1 병용, ⑤ SGLT-2 억제제 + 메트포르민, ⑥ 다중 대사 표적 병용(메트포르민 + GLP-1 + 스타틴 + ARB)
- **보정/근거**: UKPDS(메트포르민), LEADER(리라글루타이드), EMPA-REG OUTCOME, SUSTAIN-6(세마글루타이드) 주요 심혈관 결과 시험 데이터 기반

## Shiny 대시보드 (Dashboard)
환자 프로파일 입력(체중·허리둘레·혈압·HbA1c·지질·기저 심혈관 위험), 약물 PK 농도 추이, 주요 PD 바이오마커(혈당·인슐린·지질·체중·CRP), 임상 엔드포인트(HbA1c·혈압·LDL-C·체중 목표 달성 예측), 치료 시나리오 비교, 심혈관 위험 감소 시각화 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("ms_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("ms_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg ms_qsp_model.dot -o ms_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [ms_references.md](ms_references.md) 참조 (약 50편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
