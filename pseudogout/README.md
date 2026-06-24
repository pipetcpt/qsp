# 가성통풍 (CPPD) (Pseudogout (CPPD), CPPD) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![CPPD QSP Model](cppd_qsp_model.png)](cppd_qsp_model.svg)

## 개요 (Overview)

가성통풍(Pseudogout)은 칼슘피로인산(CPPD, calcium pyrophosphate dihydrate) 결정이 관절 연골 및 활막에 침착되어 발생하는 결정 유발 관절염입니다. 주로 무릎, 손목 등 대관절에 갑작스러운 염증 발작을 일으키며, 유병률은 60세 이상에서 현저히 증가합니다. NLRP3 인플라마솜이 CPPD 결정을 인식하여 IL-1β·IL-18을 대량 방출하는 것이 핵심 발병기전이며, 콜히친·NSAID·스테로이드가 주된 치료 표적입니다. 이차성 CPPD는 부갑상선 기능 항진증, 저마그네슘혈증, 혈색소증 등 대사 이상과 연관됩니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| CPPD 결정 형성 | ANKH·ENPP1에 의한 무기 피로인산(PPi) 과생성, 트랜스글루타미나제 2 활성화 | 연골 석회화, 방사선 연골석회증 |
| NLRP3 인플라마솜 활성화 | CPPD 결정 → K+ 유출 → NLRP3/ASC/Caspase-1 복합체 조립 | IL-1β·IL-18 성숙·분비 |
| 중성구 침윤 | IL-1β 유발 CXCL8·C5a → 관절 내 중성구 대량 동원 | 급성 관절 염증, 통증 |
| COX-2/프로스타글란딘 경로 | 대식세포·시노비오사이트의 PGE2 생성 증가 | 혈관 확장, 통증 증폭 |
| 기질 금속단백분해효소 (MMP) | IL-1β·TNF → MMP-3·MMP-13 발현 | 연골 분해, 관절 손상 |
| 대사 이상 연관 경로 | PTH 과잉 → 세포외 Ca²⁺/PPi 상승 | 이차 CPPD 촉진 |
| 만성 저등급 염증 | M1 대식세포·TGF-β → 섬유화, 만성 관절병증 | 퇴행성 관절 손상 |

## 주요 약물 표적 (Drug Targets)

- **콜히친 (Colchicine)**: 튜불린 중합 억제 → 중성구 화학주성·NLRP3 활성화 차단; 발작 예방 및 급성기 치료
- **NSAID (나프록센, 인도메타신)**: COX-1/2 억제 → PGE2 감소; 급성 발작 1차 치료
- **코르티코스테로이드 (프레드니솔론, 트리암시놀론 관절강내 주사)**: 광범위 항염 → NF-κB 억제; NSAID 금기 시 대안
- **IL-1 차단제 (아나킨라, 카나키누맙)**: IL-1β 직접 중화; 재발성·난치성 발작
- **메토트렉세이트/하이드록시클로로퀸**: 만성 재발 억제 (오프라벨)

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [cppd_qsp_model.dot](cppd_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 100+ 노드 / 12 클러스터) |
| [cppd_qsp_model.svg](cppd_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [cppd_qsp_model.png](cppd_qsp_model.png) | PNG 이미지 (150 dpi) |
| [cppd_mrgsolve_model.R](cppd_mrgsolve_model.R) | mrgsolve ODE 모델 (약 21 구획 / 약 27개 시나리오) |
| [cppd_shiny_app.R](cppd_shiny_app.R) | Shiny 대시보드 |
| [cppd_references.md](cppd_references.md) | 참고문헌 (약 62편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 콜히친/NSAID/스테로이드/IL-1Ra의 2구획 PK + NLRP3 활성화 → IL-1β 생성 → 중성구 동원 → 관절 염증 지수(JII)로 이어지는 PD 모듈, 연골 PPi 축적 구획 포함
- **주요 치료 시나리오**: 무치료 자연 경과, 콜히친 단독, NSAID 단독, 스테로이드 관절강내 주사, 아나킨라 피하, 카나키누맙 단회, 콜히친+NSAID 병용, 예방적 저용량 콜히친 장기투여 등
- **보정/근거**: Dalbeth 등 *Lancet* 2021(CPPD 역학), Ea 등 *Arthritis Rheum* NLRP3 기전 연구, CRESCENT 임상 데이터(카나키누맙), 콜히친 PK 문헌(Terkeltaub 등) 기반으로 파라미터 설정

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) **환자 프로파일** — 나이·동반 질환·유발 인자 설정; (2) **PK 프로파일** — 약물별 혈중/조직 농도 시간 경과; (3) **염증 지표** — IL-1β·중성구·CRP 동태; (4) **임상 엔드포인트** — 관절 통증 점수·발작 빈도 시뮬레이션; (5) **시나리오 비교** — 단일/병용 치료 간 비교; (6) **바이오마커** — 혈청 PPi, uric acid, 관절액 결정 지수.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("cppd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("cppd_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg cppd_qsp_model.dot -o cppd_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [cppd_references.md](cppd_references.md) 참조 (약 62편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
