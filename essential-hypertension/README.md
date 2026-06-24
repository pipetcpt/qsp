# 본태성 고혈압 (Essential Hypertension, EH) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 심혈관

[![EH QSP Model](eh_qsp_model.png)](eh_qsp_model.svg)

## 개요 (Overview)

본태성 고혈압(Essential Hypertension)은 이차적 원인 없이 발생하는 혈압 상승(수축기 ≥ 130 mmHg 또는 이완기 ≥ 80 mmHg)으로, 전 세계 성인의 약 30~45%가 이환되어 있습니다. 심근경색·뇌졸중·만성 신부전의 주요 원인으로, 전 세계 사망 원인 1위 위험인자입니다. 핵심 병태생리는 RAAS 과활성화, 교감신경 긴장도 상승, 신장 나트륨·체액 항상성 이상이 심장 출력(CO)과 말초혈관저항(TPR)을 높이는 복합 기전입니다. ACEi/ARB·CCB·이뇨제·베타차단제가 1차 치료제입니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| RAAS 축 | 레닌 → 안지오텐시노겐 → Ang I → Ang II → AT1R | 혈관수축, 알도스테론 분비, Na 저류 |
| 교감신경 활성화 | 노르에피네프린 → α1R(혈관수축)·β1R(심박수·수축력) | 심박출량·TPR 상승 |
| 신장 나트륨 항상성 | Ang II·알도스테론 → 집합관 Na 재흡수 | 체액 과잉, 혈압 상승 |
| 내피 기능장애 | Ang II → ROS → NO 감소 | 혈관 이완능 저하 |
| 압력-나트륨뇨 관계 | 혈압 상승 → 신압 증가 → 나트륨뇨 지연 | 혈압 재설정(set-point 이상) |
| 혈관 리모델링 | 만성 고혈압 → 중막 비후 → TPR 증가 | 고혈압 고착화 |
| 장기 손상 | 심실 비대, 단백뇨, 망막 병변 | 심혈관·신장 합병증 |

## 주요 약물 표적 (Drug Targets)

- **ACE억제제** (라미프릴, 에날라프릴): ACE 억제 → Ang II 감소, 브라디키닌 축적 → 혈관확장
- **ARB** (로사르탄, 발사르탄): AT1R 직접 차단 → Ang II 효과 차단
- **CCB** (암로디핀): L형 칼슘채널 차단 → 혈관평활근 이완, TPR 감소
- **베타차단제** (비소프롤롤): β1AR 차단 → 심박수·수축력 감소, 레닌 분비 억제
- **티아지드 이뇨제** (HCTZ): 원위세뇨관 Na-Cl 공동수송체 억제 → 나트륨뇨 증가

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [eh_qsp_model.dot](eh_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 173 노드 / 10 클러스터) |
| [eh_qsp_model.svg](eh_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [eh_qsp_model.png](eh_qsp_model.png) | PNG 이미지 (150 dpi) |
| [eh_mrgsolve_model.R](eh_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 6 치료 시나리오) |
| [eh_shiny_app.R](eh_shiny_app.R) | Shiny 대시보드 |
| [eh_references.md](eh_references.md) | 참고문헌 (약 41편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK 구획(ACEi·ARB·CCB·베타차단제·HCTZ 각 중심·말초) + PD 구획(레닌, 안지오텐신 II, 알도스테론, 체액량, 심박출량, 말초혈관저항, 수축기/이완기혈압, 나트륨뇨, 심박수)
- **주요 치료 시나리오**: ① 무치료(미치료 고혈압), ② ACE억제제 단독(라미프릴 10 mg QD), ③ ARB 단독(로사르탄 100 mg QD), ④ CCB 단독(암로디핀 10 mg QD), ⑤ 베타차단제 단독(비소프롤롤 10 mg QD), ⑥ 3제 병용(ACEi + CCB + 티아지드) — 표준 1차 복합요법
- **보정/근거**: ALLHAT(암로디핀/리시노프릴/클로르탈리돈), HOT 시험, HOPE 시험(라미프릴) 데이터 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① 환자 프로파일(기저 혈압·신장 기능·염분 섭취·BMI 설정), ② PK 탭(5개 약물 혈중 농도 및 활성 대사체), ③ RAAS/혈역학 PD 탭(Ang II·알도스테론·TPR·CO 추이), ④ 임상 엔드포인트(수축기·이완기 혈압·심박수·혈압 목표 달성), ⑤ 시나리오 비교(6개 치료 전략 동시 비교), ⑥ 바이오마커(레닌·BNP·신기능·나트륨뇨 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("eh_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("eh_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg eh_qsp_model.dot -o eh_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [eh_references.md](eh_references.md) 참조 (약 41편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
