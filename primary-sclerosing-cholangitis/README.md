# 원발성 경화성 담관염 (PSC) (Primary Sclerosing Cholangitis, PSC) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도

[![PSC QSP Model](psc_qsp_model.png)](psc_qsp_model.svg)

## 개요 (Overview)
원발성 경화성 담관염(PSC)은 간내·외 담관의 다발성 협착과 섬유화가 진행되는 만성 담즙정체성 간질환으로, 인구 10만 명당 약 5~10명에서 유병한다. 70~80%에서 염증성 장질환(IBD, 주로 궤양성 대장염)이 동반되며, 진단 후 중앙 생존기간은 12~15년이다. 장내 세균총 이상에 의한 LPS 장간 이동(gut-liver axis) → 담관 상피세포 손상 · TH17-면역 불균형 → LOXL2 매개 간성상세포 활성화 → 담관 주위 섬유화·협착의 기전이 핵심이다. 담관암(CCA) 발생 위험이 연간 1.5~2%로 크게 높아 정기 감시가 필수적이다. 현재 효과적인 약물 치료가 없어 증상 완화와 최종적으로 간이식이 유일한 근치법이다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 장-간 축(Gut-Liver Axis) | 장벽 기능 저하 → LPS 문맥 이동 → TLR4 활성화 | 담관 상피 염증 유발 |
| 담관 상피 손상 경로 | TH17/IL-17A, TNF-α, 담관 상피 노화(senescence) | 담관 상피 세포 탈락 |
| FXR-담즙산 조절 경로 | FXR 활성 저하 → 담즙산 합성 과다 → 독성 축적 | 담관 독성·담즙정체 |
| 간성상세포 활성화 경로 | TGF-β·LOXL2 → HSC 활성화·콜라겐 가교 | 담관 주위 섬유화·협착 |
| 담관 섬유화 진행 경로 | Col1a1 축적, Fibroscan 경도 증가 | 문맥 고혈압·간경변 |
| 담관암 위험 경로 | 만성 담즙정체·DNA 손상·염증 | CCA 위험 연 1.5~2% |

## 주요 약물 표적 (Drug Targets)
- **UDCA**: 담즙산 친수성 개선·세포 보호 — PSC에서 생존 개선 미입증, 고용량에서 해로울 수 있어 가이드라인별 의견 분분
- **NorUDCA (24-노르우르소데옥시콜산)**: 담관 친수성 담즙산 대체, PSC 특이 Phase 2/3 시험
- **오베티콜산 (FXR 작용제)**: FXR 활성 → 담즙산 합성 억제·항섬유화 — PSC Phase 3 진행 중
- **베자피브레이트 (PPARα)**: 담즙산 독성 감소·항염 — PSC 임상 탐색 중
- **항생제 (메트로니다졸·반코마이신)**: 장내 세균총 조절 — PSC에서 소규모 연구 진행
- **간이식**: 말기 PSC·담관암 고위험 — 근치적 치료

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [psc_qsp_model.dot](psc_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 213 노드 / 11 클러스터) |
| [psc_qsp_model.svg](psc_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [psc_qsp_model.png](psc_qsp_model.png) | PNG 이미지 (150 dpi) |
| [psc_mrgsolve_model.R](psc_mrgsolve_model.R) | mrgsolve ODE 모델 (약 27 구획 / 5 치료 시나리오) |
| [psc_shiny_app.R](psc_shiny_app.R) | Shiny 대시보드 |
| [psc_references.md](psc_references.md) | 참고문헌 (약 43편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: UDCA·OCA·NorUDCA·베자피브레이트 PK 구획(담즙·혈장 분리), LPS·장벽·FXR·담즙산 풀·담즙 친수성지수·IL-17A·TNF-α·IL-6·Treg/IL-10·담관 상피 건강·담관 노화·간성상세포·Col1a1·LOXL2·ALP·빌리루빈·Fibroscan·문맥압·담관암 위험 PD 구획
- **주요 치료 시나리오**: ① 무치료 자연경과, ② UDCA 15 mg/kg/일, ③ OCA 10 mg/일, ④ UDCA+OCA 병용, ⑤ 베자피브레이트 400 mg/일
- **보정/근거**: PSC 자연경과 코호트(Boonstra 2013) 및 PRIMROSE(NorUDCA), AESOP(OCA for PSC) 임상시험 데이터를 기반으로 ALP·섬유화 진행 속도 파라미터 보정

## Shiny 대시보드 (Dashboard)
환자 프로파일(IBD 동반 여부, 담관 협착 위치, Amsterdam PSC 위험 점수) · 담즙산/FXR PK/PD · 담관 염증·섬유화 PD · 간 기능 임상 엔드포인트(ALP·빌리루빈·Fibroscan) · 치료 시나리오 비교(ALP 반응·섬유화 억제) · 담관암 위험·바이오마커(LOXL2·IL-17A·Col1a1) 탭으로 구성

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("psc_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("psc_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg psc_qsp_model.dot -o psc_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [psc_references.md](psc_references.md) 참조 (약 43편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
