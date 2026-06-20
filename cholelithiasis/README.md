# 담석증 (Cholelithiasis, CHOL) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도

[![CHOL QSP Model](chol_qsp_model.png)](chol_qsp_model.svg)

## 개요 (Overview)

담석증은 전 세계 성인의 약 10~20%에서 발생하는 흔한 소화기 질환으로, 특히 서구 국가와 고령 여성에서 유병률이 높다. 대부분(~80%)은 콜레스테롤 과포화 담즙에서 결정핵이 형성되는 콜레스테롤 담석이며, 담낭 운동 장애(정체)가 결석 성장의 핵심 촉진 인자이다. 우르소데옥시콜산(UDCA)은 담즙 내 콜레스테롤 포화지수를 낮추고 결석을 용해하는 유일한 경구 약물 치료로, 스타틴·에제티미브와의 병용을 통해 콜레스테롤 합성 및 흡수를 추가 억제한다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 담즙 콜레스테롤 과포화 | ABCG5/G8 과발현, HMGCR 활성 증가 | 콜레스테롤 포화지수(CSI) 상승 |
| 핵형성 촉진 | 뮤신 글리코단백질·프로뉴클레이팅 인자 분비 | 콜레스테롤 결정 형성 |
| 담낭 정체·운동 장애 | CCK 분비 감소, 평활근 수축 저하 | 결정 체류 시간 연장 → 결석 성장 |
| 장간 순환 교란 | 담즙산 재흡수 감소, 콜레스테롤/인지질 비율 증가 | 담즙산 풀(pool) 축소 |
| 만성 염증·점막 반응 | IL-6, CRP 상승, 담낭벽 비후 | 담낭염 유발, 통증 |
| 스타틴 효과 | HMGCR 억제 → 간내 콜레스테롤 합성 감소 | 담즙 콜레스테롤 농도 저하 |
| UDCA 기전 | 콜레스테롤 결정화 억제, 담즙 친수성 증가 | 소형 결석 용해 촉진 |

## 주요 약물 표적 (Drug Targets)

- **우르소데옥시콜산(UDCA)**: 담즙산 조성 정상화, CSI 감소, 소형 방사선 투과성 결석 용해 (Actigall 등)
- **스타틴(HMG-CoA 환원효소 억제제)**: 로수바스타틴·심바스타틴 — 간내 콜레스테롤 합성 억제
- **에제티미브**: NPC1L1 차단 → 장내 콜레스테롤 흡수 억제
- **담즙산 격리제(콜레스티라민)**: 장간 순환 차단, 담즙산 손실 유도
- **CCK 수용체 작용제(Sincalide)**: 담낭 수축 촉진, 담즙 배출 개선

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [chol_qsp_model.dot](chol_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 431 노드 / 10 클러스터) |
| [chol_qsp_model.svg](chol_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [chol_qsp_model.png](chol_qsp_model.png) | PNG 이미지 (150 dpi) |
| [chol_mrgsolve_model.R](chol_mrgsolve_model.R) | mrgsolve ODE 모델 (약 25 구획 / 5개 치료 시나리오) |
| [chol_shiny_app.R](chol_shiny_app.R) | Shiny 대시보드 |
| [chol_references.md](chol_references.md) | 참고문헌 (약 46편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: UDCA/스타틴/에제티미브 각각의 흡수·혈장·간장·담즙·담낭 PK 구획 + 담즙산 풀, 간내 콜레스테롤, 담즙 콜레스테롤/인지질, 담낭 용적, 결정 질량, 결석 용적, IL-6, CRP PD 구획
- **주요 치료 시나리오**: (1) 무치료 자연경과, (2) UDCA 단독, (3) 스타틴 단독, (4) UDCA + 스타틴 병용, (5) UDCA + 에제티미브 병용
- **보정/근거**: UDCA 용해 효과는 Portincasa et al. (Lancet 2006), 스타틴의 CSI 감소 효과는 관련 RCT 데이터를 기반으로 파라미터 보정

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) 환자 프로파일 — 위험인자(비만·연령·여성·식이) 설정; (2) PK 탭 — UDCA/스타틴/에제티미브 혈장 농도 시계열; (3) PD 주요지표 — CSI, 담즙산 풀, 담낭 운동; (4) 임상 엔드포인트 — 결석 용적 변화, 용해율; (5) 시나리오 비교 — 치료 전략별 1년 결과 비교; (6) 바이오마커 — IL-6, CRP, 담즙 인지질/콜레스테롤 비율

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("chol_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("chol_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg chol_qsp_model.dot -o chol_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [chol_references.md](chol_references.md) 참조 (약 46편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
