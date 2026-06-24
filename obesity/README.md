# 비만 (Obesity, OB) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![OB QSP Model](ob_qsp_model.png)](ob_qsp_model.svg)

## 개요 (Overview)

비만은 BMI ≥ 30 kg/m²으로 정의되는 만성 대사질환으로, 전 세계 성인 유병률은 약 16%(약 8억 명)에 달하며 계속 증가하고 있다. 핵심 발병기전은 에너지 섭취 및 소비 항상성의 조절이상으로, 시상하부 멜라노코르틴-4 수용체(MC4R) 경로의 렙틴 저항성, GLP-1/GIP 분비 감소, 위 배출 가속화, 인슐린 저항성 등이 상호작용한다. GLP-1 및 GIP/GLP-1 이중 수용체 작용제(세마글루타이드·티르제파타이드)가 임상에서 10–22%의 체중 감량을 달성하며 치료 패러다임을 바꾸고 있다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 렙틴-멜라노코르틴 축 | 렙틴 저항성 → MC4R 활성 저하 → 식욕 억제 실패 | 과식, 체중 증가 |
| GLP-1/GIP 분비 감소 | L·K세포 기능 저하 → 인크레틴 효과 감소 | 인슐린 분비 불충분, 위 배출 가속 |
| 인슐린 저항성 | IRS-1/PI3K/Akt 신호 장애 → 간·근육·지방 인슐린 반응 저하 | 공복혈당 상승, HOMA-IR 증가 |
| 지방조직 만성 염증 | M1 대식세포 침윤 → TNF-α, IL-6 → 지방세포 인슐린 저항성 | 전신 만성 염증 |
| 에너지 소비 저하 | 적응성 열발생(adaptive thermogenesis) 감소 | 체중 감량 후 반동 체중 증가 |
| 장내 미생물군 이상 | Firmicutes/Bacteroidetes 비율 증가 → 에너지 추출 증가 | 비만 악화, 지방간 |
| 그렐린·NPY/AgRP 상승 | 식욕 촉진 펩타이드 상승 | 식욕 증가, 체중 증가 지속 |

## 주요 약물 표적 (Drug Targets)

- **세마글루타이드 (Semaglutide, GLP-1RA)**: GLP-1 수용체 활성화 → 식욕 억제·위 배출 지연·인슐린 분비 촉진 (STEP-1: 약 15% 체중 감량)
- **티르제파타이드 (Tirzepatide, GLP-1/GIPR 이중 작용제)**: GLP-1+GIP 수용체 동시 활성화 → 약 21% 체중 감량 (SURMOUNT-1)
- **오를리스탯 (Orlistat, 지방분해효소 억제제)**: 위장관 지방 흡수 약 30% 차단
- **펜터민/토피라메이트 (Phentermine/Topiramate)**: 중추 교감신경 활성화 + GABA 수용체 조절 → 식욕 억제
- **날트렉손/부프로피온 (Naltrexone/Bupropion)**: 멜라노코르틴 경로 강화, 보상 경로 억제

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [ob_qsp_model.dot](ob_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 199 노드 / 11 클러스터) |
| [ob_qsp_model.svg](ob_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [ob_qsp_model.png](ob_qsp_model.png) | PNG 이미지 (150 dpi) |
| [ob_mrgsolve_model.R](ob_mrgsolve_model.R) | mrgsolve ODE 모델 (약 20 구획 / 5개 치료 시나리오) |
| [ob_shiny_app.R](ob_shiny_app.R) | Shiny 대시보드 |
| [ob_references.md](ob_references.md) | 참고문헌 (약 45편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 세마글루타이드 SC 2구획, 티르제파타이드 SC 2구획, 오를리스탯 경구 2구획, CNS(시상하부) 1구획, GLP1R/GIPR 점유율 구획; 음식 섭취·위 배출·그렐린·인슐린·혈당·체지방·체중·렙틴·중성지방·HbA1c·전신 염증·HOMA-IR
- **주요 치료 시나리오**: ① 위약 ② 세마글루타이드 2.4mg SC QW ③ 티르제파타이드 15mg SC QW ④ 오를리스탯 120mg TID ⑤ 펜터민/토피라메이트 경구 QD
- **보정/근거**: STEP-1(세마글루타이드), SURMOUNT-1(티르제파타이드), XENDOS(오를리스탯) 임상시험 체중 감량 곡선 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① **환자 프로파일** (기저 BMI·체중·당뇨 유무·렙틴 저항성 설정), ② **PK** (혈장 약물 농도 및 GLP-1R/GIPR 수용체 점유율), ③ **PD 주요지표** (식욕·혈당·인슐린·HOMA-IR 추이), ④ **임상 엔드포인트** (체중 감량률·BMI·허리둘레 변화), ⑤ **시나리오 비교** (5개 치료 전략 직접 비교), ⑥ **바이오마커** (렙틴·중성지방·HbA1c·전신 염증 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("ob_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("ob_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg ob_qsp_model.dot -o ob_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [ob_references.md](ob_references.md) 참조 (약 45편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
