# 과민성 장증후군 (IBS) (Irritable Bowel Syndrome, IBS) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도

[![IBS QSP Model](ibs_qsp_model.png)](ibs_qsp_model.svg)

## 개요 (Overview)
과민성 장증후군(IBS)은 기질적 이상 없이 복통과 배변 습관 변화를 특징으로 하는 만성 기능성 장 질환으로, 전 세계 유병률은 약 11%에 달합니다. 핵심 발병기전은 뇌-장 축(brain-gut axis) 이상, 내장 과민성, 장 운동 이상, 저등급 점막 염증, 장내 미생물 불균형(dysbiosis), 장 투과성 증가가 복잡하게 상호작용하는 다인성 장애입니다. 변비 우세형(IBS-C), 설사 우세형(IBS-D), 혼합형(IBS-M), 분류불가형(IBS-U)의 4가지 아형으로 나뉩니다. 치료는 식이 조절(저 FODMAP), 신경조절제(삼환계 항우울제·SSRI), 세로토닌 조절제(알로세트론·테가세로드), 장 특이 약물(리나클로타이드·루비프로스톤)을 포함합니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 뇌-장 축 이상 | CRF/스트레스 → HPA 축 과활성화, 장 신경계 과민 | 복통·복부 팽만 악화 |
| 5-HT(세로토닌) 신호 이상 | 장 EC 세포의 5-HT 분비 이상, SERT 기능 저하 | 운동·분비·내장 감각 조절 장애 |
| 내장 과민성 | TRPV1·TRPA1 발현 증가, 척수 감작 | 통증 역치 저하 |
| 장 투과성 증가 | Tight junction 단백(occludin·claudin) 이상, 비만세포 활성화 | 세균 성분(LPS) 유입·저등급 염증 |
| 장내 미생물 불균형 | Firmicutes/Bacteroidetes 비율 변화, SCFA 생성 감소 | 운동·면역·신경 조절 이상 |
| 장 운동 이상 | 5-HT4/5-HT3 수용체·클로라이드 채널 조절 이상 | IBS-D/IBS-C 아형 분류 |

## 주요 약물 표적 (Drug Targets)
- **리나클로타이드 / 플레카나타이드** (GC-C 작용제): cGMP 상승 → 장액 분비 촉진, 통증 억제 (IBS-C)
- **루비프로스톤** (클로라이드 채널 활성화제): 장액 분비 증가 (IBS-C)
- **알로세트론** (5-HT3 길항제): 결장 운동 억제, 내장 감각 억제 (중증 IBS-D)
- **삼환계 항우울제(TCA)** (아미트립틸린): 내장 통증 조절, 5-HT/NE 재흡수 억제
- **리팍시민** (비흡수성 항생제): 장내 미생물 불균형 교정 (IBS-D, 비변비형)
- **메베베린 / 히오신부틸브로마이드** (항경련제): 장 평활근 이완, 복통 완화

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [ibs_qsp_model.dot](ibs_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 136 노드 / 9 클러스터) |
| [ibs_qsp_model.svg](ibs_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [ibs_qsp_model.png](ibs_qsp_model.png) | PNG 이미지 (150 dpi) |
| [ibs_mrgsolve_model.R](ibs_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 5개 치료 시나리오) |
| [ibs_shiny_app.R](ibs_shiny_app.R) | Shiny 대시보드 |
| [ibs_references.md](ibs_references.md) | 참고문헌 (약 55편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 뇌-장 축 모듈(스트레스·CRF·코르티솔), 세로토닌 모듈(장 5-HT·SERT 점유), 비만세포 활성화, 염증, 장 투과성, 미생물(dysbiosis), SCFA, 내장 과민성, 장 운동, 통증·팽만·배변 불편감 PD 출력 구획 + 약물 PK 구획(GC-C 작용제 2구획, 리팍시민 1구획, 알로세트론 2구획)
- **주요 치료 시나리오**: ① 자연경과(무치료), ② 리나클로타이드(IBS-C), ③ 알로세트론(IBS-D), ④ 리팍시민 단기 치료, ⑤ TCA + 심리치료 병용
- **보정/근거**: LINACLOTIDE-301/302 시험, TARGET 1/2 시험(알로세트론), TARGET 3 시험(리팍시민) 등 주요 3상 임상 데이터 기반

## Shiny 대시보드 (Dashboard)
환자 프로파일 입력(아형·기저 증상 중증도·스트레스 수준), 약물 PK 농도 추이, 뇌-장 축 PD 바이오마커(5-HT·CRF·코르티솔), 임상 엔드포인트(IBS-SSS 점수·복통·배변 빈도 변화), 치료 시나리오 비교, 아형별 반응 예측 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("ibs_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("ibs_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg ibs_qsp_model.dot -o ibs_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [ibs_references.md](ibs_references.md) 참조 (약 55편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
