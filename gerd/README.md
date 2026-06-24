# 위식도 역류질환 (Gastroesophageal Reflux Disease, GERD) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도

[![GERD QSP Model](gerd_qsp_model.png)](gerd_qsp_model.svg)

## 개요 (Overview)

위식도 역류질환(GERD)은 위산과 위 내용물이 식도로 역류하여 점막 손상 및 증상(속쓰림·역류)을 일으키는 질환으로, 서구 성인의 약 10~20%에서 유병률을 보입니다. 핵심 병태생리는 하부식도괄약근(LES)의 일과성 이완(TLESRs) 또는 기능 저하, 식도 산 제거 능력 감소, 식도 점막 방어 기전 이상입니다. 미란성 역류식도염(ERD)에서 바렛 식도까지 스펙트럼이 넓으며, 위산 분비를 강력히 억제하는 PPI(양성자펌프억제제)와 P-CAB(칼륨 경쟁적 산 분비 억제제)이 주요 치료제입니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| LES 기능 이상 | TLESRs 빈도 증가·기저 LES압 저하 | 산 역류 횟수 증가 |
| 위산 분비 | H+/K+-ATPase → 위산 pH 1~2 유지 | 역류 시 식도 손상 |
| 식도 산 노출 | 역류 산 → 식도 점막 상피 손상 | DeMeester 점수 상승 |
| 점막 방어 기전 | 점액·중탄산염·상피 증식 감소 | 미란·궤양 형성 |
| 염증 캐스케이드 | 산·펩신 → IL-8, TNF-α → 호중구 침윤 | 점막 염증, 미란 |
| 바렛 식도 진행 | 만성 산 노출 → 장형 화생 | 식도선암 위험 증가 |
| 식도 운동 이상 | 식도 연동 운동 약화 → 산 제거 지연 | 야간 역류 증가 |

## 주요 약물 표적 (Drug Targets)

- **PPI** (오메프라졸, 에소메프라졸, 판토프라졸): 비가역적 H+/K+-ATPase 억제 → 위산 pH ≥ 4 유지, 표준 치료
- **P-CAB** (보노프라잔, 테고프라잔): 칼륨 경쟁적 가역적 H+/K+-ATPase 억제 → 더 빠른 효과·야간 위산 억제 우수
- **H2 수용체 길항제** (파모티딘, 라니티딘): 히스타민 H2R 차단 → 보조적 위산 억제
- **프로키네틱제** (메토클로프라미드, 도파민 길항제): 위 배출 촉진, LES 압력 상승
- **알긴산염/항산제**: 물리적 산 중화 및 산 역류 방벽

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [gerd_qsp_model.dot](gerd_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 176 노드 / 11 클러스터) |
| [gerd_qsp_model.svg](gerd_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [gerd_qsp_model.png](gerd_qsp_model.png) | PNG 이미지 (150 dpi) |
| [gerd_mrgsolve_model.R](gerd_mrgsolve_model.R) | mrgsolve ODE 모델 (약 23 구획 / 6 치료 시나리오) |
| [gerd_shiny_app.R](gerd_shiny_app.R) | Shiny 대시보드 |
| [gerd_references.md](gerd_references.md) | 참고문헌 (약 41편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK 구획(PPI·H2RA·P-CAB 각 GI/중심 compartment) + PD 구획(H+/K+-ATPase 활성, 위내 pH, 위산 분비속도, LES 압력/TLESRs 빈도, 식도 산 노출 시간, 점막 손상 지수, 바렛 진행 위험, 염증 지수)
- **주요 치료 시나리오**: ① 무치료, ② PPI 표준 용량(오메프라졸 20 mg QD), ③ PPI 고용량(에소메프라졸 40 mg QD), ④ P-CAB(보노프라잔 20 mg QD), ⑤ H2RA(파모티딘 20 mg BID), ⑥ PPI + H2RA 병용(야간 증상 돌파 억제)
- **보정/근거**: Metz 등 PPI PD 모델, 보노프라잔 PHALCON-EE 시험, 식도 pH-임피던스 모니터링 정상값 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① 환자 프로파일(기저 위산 분비·LES 기능·미란 등급·H. pylori 상태 설정), ② PK 탭(PPI/P-CAB/H2RA 혈중 농도·CYP2C19 대사 유형), ③ 위산/식도 PD 탭(위내 pH·산 노출 시간·TLESRs 빈도 추이), ④ 임상 엔드포인트(증상 완해율·미란 치유율·바렛 위험), ⑤ 시나리오 비교(6개 치료 전략 동시 비교), ⑥ 바이오마커(H+/K+-ATPase 활성·점막 손상 지수·염증 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("gerd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("gerd_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg gerd_qsp_model.dot -o gerd_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [gerd_references.md](gerd_references.md) 참조 (약 41편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
