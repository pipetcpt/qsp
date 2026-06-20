# 반응성 관절염 (Reactive Arthritis, ReA) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 자가면역·류마티스

[![ReA QSP Model](rea_qsp_model.png)](rea_qsp_model.svg)

## 개요 (Overview)

반응성 관절염(ReA)은 요도·장·인두 등 원격 감염 후 수 주 내에 발생하는 무균성 관절염으로, HLA-B27 양성자에서 빈발합니다. 주요 원인균은 *Chlamydia trachomatis*, *Salmonella*, *Shigella*, *Yersinia*, *Campylobacter* 등이며, 이들의 항원 성분이 관절 내에서 T세포를 지속 자극하여 비대칭 하지 관절염·부착부염·요도염·결막염의 삼주징(舊 라이터 증후군)을 초래합니다. 급성기는 NSAID와 항생제가 기반이 되며, 만성·재발 병례에서 설파살라진, 면역억제제, 생물학제제(항TNF)가 사용됩니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| HLA-B27 연관 면역 | HLA-B27 이중 사슬 오접힘(UPR) → IL-23 과발현 | Th17 경로 활성화, 부착부염 소인 |
| 감염 항원 지속 | 세균 LPS·펩타이드 → 관절 대식세포 TLR 활성화 | IL-1β·TNF·IL-6 분비 지속 |
| Th1/Th17 활성화 | CD4+ Th1(IFN-γ)·Th17(IL-17) → 시노비오사이트·RANKL 유도 | 활막 염증, 골미란 |
| CD8+ T세포 세포독성 | HLA-B27 → 세균 펩타이드 제시 → CTL 활성화 → 관절 조직 손상 | 만성 활막 손상 |
| RANKL/OPG 불균형 | 염증성 사이토카인 → RANKL 상승 → 파골세포 활성화 | 관절 주위 골 손실 |
| 비대칭 부착부염 | IL-17 + 역학적 스트레스 → 부착부 섬유아세포 활성화 | 종골·무릎 부착부 통증 |
| 항균 면역-자가면역 교차 | 분자 모방(molecular mimicry) → 자가항원 반응 | 재발·만성화 경향 |

## 주요 약물 표적 (Drug Targets)

- **NSAID (나프록센, 인도메타신, 디클로페낙)**: COX 억제 → PGE2 감소; 급성 관절 통증·종창의 1차 치료
- **항생제 (독시사이클린, 아지트로마이신)**: 지속 세균 항원 제거; *Chlamydia* 관련 ReA에서 효과적
- **코르티코스테로이드**: 관절강내 주사 또는 전신; NSAID 불충분 시 단기 보조
- **설파살라진**: 말초 관절 만성화 시 DMARD; 임상적 근거 중등도
- **메토트렉세이트/아자티오프린**: 만성 재발성 병례 면역억제
- **항TNF 제제 (에타너셉트, 아달리무맙)**: HLA-B27 양성 척추 침범·난치성 병례 오프라벨 사용

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [rea_qsp_model.dot](rea_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 100+ 노드 / 10 클러스터) |
| [rea_qsp_model.svg](rea_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [rea_qsp_model.png](rea_qsp_model.png) | PNG 이미지 (150 dpi) |
| [rea_mrgsolve_model.R](rea_mrgsolve_model.R) | mrgsolve ODE 모델 (약 29 구획 / 약 47개 시나리오) |
| [rea_shiny_app.R](rea_shiny_app.R) | Shiny 대시보드 |
| [rea_references.md](rea_references.md) | 참고문헌 (약 55편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: NSAID/항생제/설파살라진/항TNF의 PK 구획 + 세균 항원 부하 감쇠 모듈, Th1·Th17 활성화 동태, IL-17·TNF·IL-1β 생성 구획, 관절 염증 지수(JII), 부착부 점수 구획 포함
- **주요 치료 시나리오**: 무치료 자연 경과, NSAID 단독, 독시사이클린 단독, NSAID+항생제, 스테로이드 관절강내 주사, 설파살라진, 메토트렉세이트, 에타너셉트, 에타너셉트+메토트렉세이트 병용 등 다수
- **보정/근거**: Hannu 등 Best Pract Res Clin Rheumatol 2011, Yli-Kerttula 등 *Ann Rheum Dis* 항생제 RCT, ASAS 척추관절염 분류기준 관련 임상 데이터 기반

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) **환자 프로파일** — 선행 감염 유형·HLA-B27 상태·동반 증상; (2) **PK 프로파일** — NSAID/항생제/생물학제제 혈중 농도; (3) **PD 주요지표** — 세균 항원 부하·Th17·IL-17 억제 동태; (4) **임상 엔드포인트** — 관절 통증(VAS)·CRP·부착부 스코어; (5) **시나리오 비교** — 치료 전략별 염증 해소 비교; (6) **바이오마커** — ESR, CRP, HLA-B27, 활막액 세포수.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("rea_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("rea_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg rea_qsp_model.dot -o rea_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [rea_references.md](rea_references.md) 참조 (약 55편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
