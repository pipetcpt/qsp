# 만성 정맥 부전 (Chronic Venous Insufficiency, CVI) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 심혈관

[![CVI QSP Model](cvi_qsp_model.png)](cvi_qsp_model.svg)

## 개요 (Overview)

만성 정맥 부전은 성인 인구의 약 25~40%가 하지 정맥류를 경험하고, 이 중 약 1~2%가 정맥성 궤양으로 진행하는 흔한 혈관 질환이다. 핵심 병태생리는 판막 기능 부전으로 인한 역류성 정맥 고혈압이며, 내피세포 기능 장애 → 백혈구 트래핑 → 세포외기질 리모델링(섬유성 커프 형성)의 연쇄가 부종과 궤양을 유발한다. MPFF(마이크로나이즈드 정제 플라보노이드 분획)와 압박 치료가 핵심 보존 치료이며, 심부정맥 혈전 후 증후군(PTS)에는 LMWH가 추가된다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 판막 부전·역류 | 정맥 판막 기능 저하 → 직립 정맥압(AVP) 상승 | 하지 부종, 피부 변색 |
| 백혈구 트래핑·활성화 | 내피 P-셀렉틴 ↑ → 호중구·단핵구 부착 | 염증 매개인자 분비, 조직 손상 |
| 내피세포 기능 장애 | 산화 스트레스, ICAM-1 ↑ → 혈관 투과성 증가 | 단백질·적혈구 누출 → 피부 침착 |
| 섬유성 커프(Fibrin Cuff) 형성 | 피브리노겐 누출 → 피부 주위 섬유화 | 산소·영양 확산 장애 → 궤양 취약성 |
| 만성 피부 염증·리모델링 | MMP/TIMP 불균형 → 피부 섬유화 | 지방 피부 경화증, 색소 침착 |
| 정맥성 궤양 형성 | 국소 허혈 + 염증 → 피부 파괴 | VCSS 상승, 삶의 질 저하 |
| PTS(혈전 후 증후군) | 심부정맥 혈전 → 재소통 불완전·판막 손상 | 중증 CVI, 궤양 재발 |

## 주요 약물 표적 (Drug Targets)

- **MPFF(Daflon/Detralex)**: 정맥 긴장도 증가, P-셀렉틴 억제, 백혈구 트래핑 감소 → 부종 감소, 궤양 치유 촉진
- **Pentoxifylline(Trental)**: 혈액 점도 개선, TNF-α 억제, 궤양 치유 보조
- **LMWH(저분자 헤파린)**: 항응고 + 항염증 → PTS 예방 및 치료
- **압박 스타킹/붕대(압박 치료)**: 정맥 내강 축소 → AVP 감소, 부종 해소 (비약물 주요 치료)
- **Sulodexide**: 내피 보호, 헤파린 황산 보충 → 단백뇨 감소, 미세혈관 보호

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [cvi_qsp_model.dot](cvi_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 501 노드 / 11 클러스터) |
| [cvi_qsp_model.svg](cvi_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [cvi_qsp_model.png](cvi_qsp_model.png) | PNG 이미지 (150 dpi) |
| [cvi_mrgsolve_model.R](cvi_mrgsolve_model.R) | mrgsolve ODE 모델 (약 20 구획 / 7개 치료 시나리오) |
| [cvi_shiny_app.R](cvi_shiny_app.R) | Shiny 대시보드 |
| [cvi_references.md](cvi_references.md) | 참고문헌 (약 54편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: MPFF·Pentoxifylline·LMWH PK (gut/SC 흡수·혈장) 구획 + 직립 정맥압(PRESS), 백혈구 트래핑(LEU), 내피세포 기능 장애(EC_DYS), 혈관 투과성(PERM), 섬유성 커프(FIBCUFF), 염증(INFLAM), 부종(EDEMA), 피부 섬유화(FIBROS), 정맥성 궤양(ULCER), VCSS 점수, 삶의 질(QOL) PD 구획
- **주요 치료 시나리오**: (1) 무치료(압박 없음), (2) MPFF 단독, (3) Pentoxifylline + 압박, (4) LMWH + 압박(PTS), (5) MPFF + Pentoxifylline + 압박, (6) 경증(압박 단독), (7) MPFF + 압박
- **보정/근거**: MPFF 효과는 MPFF ESCHAR 후속 메타분석(Eur J Vasc Endovasc Surg), Pentoxifylline은 Cochrane 리뷰(2012) 기반

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) 환자 프로파일 — CEAP 분류·기저 정맥압·궤양 크기 설정; (2) PK 탭 — MPFF/Pentoxifylline/LMWH 혈장 농도; (3) PD 주요지표 — 백혈구 트래핑, 내피 기능, 섬유성 커프; (4) 임상 엔드포인트 — VCSS 점수, 궤양 치유율, 부종 감소; (5) 시나리오 비교 — 7가지 치료 전략 1년 결과; (6) 바이오마커 — CRP, D-다이머, 삶의 질 점수

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("cvi_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("cvi_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg cvi_qsp_model.dot -o cvi_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [cvi_references.md](cvi_references.md) 참조 (약 54편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
