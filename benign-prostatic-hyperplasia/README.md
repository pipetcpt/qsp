# 양성 전립선 비대증 (Benign Prostatic Hyperplasia, BPH) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 신장·비뇨

[![BPH QSP Model](bph_qsp_model.png)](bph_qsp_model.svg)

## 개요 (Overview)
양성 전립선 비대증(BPH)은 50세 이상 남성의 약 50%, 80세 이상에서는 80% 이상이 경험하는 가장 흔한 비뇨기과 질환입니다. DHT(디하이드로테스토스테론)가 안드로겐 수용체(AR)를 활성화하여 전립선 기질·상피세포 증식을 유도하고, α1 아드레날린 수용체 매개 평활근 긴장이 더해져 하부요로증상(LUTS)을 초래합니다. 핵심 치료 표적은 5α-환원효소(SRD5A1/2, DHT 생성 차단)와 α1-아드레날린 수용체(평활근 이완)이며, PDE5 억제제가 방광경부 및 전립선 평활근에 추가 효과를 나타냅니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| DHT-AR 매개 증식 | 테스토스테론 → 5α-환원효소(SRD5A1/2) → DHT → AR 활성화 | 전립선 부피 증가, PSA 상승 |
| α1-AR 매개 긴장 | α1a/α1d 수용체, 평활근 수축 | 기능성 요도 폐쇄, 배뇨 지연 |
| 기질-상피 상호작용 | IGF-1, FGF, EGF, TGF-β 신호 | 전립선 기질 과증식 |
| 전립선 염증 | COX-2, IL-8, NF-κB 경로 | 만성 염증성 BPH 가속 |
| 방광 기능 변화 | 배뇨근 비대, 불안정 수축 | 야간뇨, 빈뇨, 절박성 |
| cGMP-PDE5 경로 | NO → sGC → cGMP 생성; PDE5 분해 | 방광경부·요도 이완 (타달라필 표적) |

## 주요 약물 표적 (Drug Targets)
- **탐수로신 (α1-차단제)**: α1a/α1d 수용체 선택적 길항 → 전립선·방광경부 평활근 이완 — 즉각 LUTS 개선
- **피나스테리드 (SRD5A2 억제제)**: DHT 합성 40~70% 감소 → 전립선 위축 — PLESS 시험 기반
- **두타스테리드 (이중 SRD5A1/2 억제제)**: DHT 90% 이상 감소 → 더 완전한 AR 차단 — COMBAT 시험 기반
- **타달라필 (PDE5 억제제)**: cGMP 증가 → 평활근 이완 — BPH-LUTS에 FDA 승인 (5 mg QD)
- **두타스테리드+탐수로신 병합**: 구조 및 기능 동시 개선 — COMBAT, CombAT 시험

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [bph_qsp_model.dot](bph_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 167 노드 / 11 클러스터) |
| [bph_qsp_model.svg](bph_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [bph_qsp_model.png](bph_qsp_model.png) | PNG 이미지 (150 dpi) |
| [bph_mrgsolve_model.R](bph_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 6개 치료 시나리오) |
| [bph_shiny_app.R](bph_shiny_app.R) | Shiny 대시보드 |
| [bph_references.md](bph_references.md) | 참고문헌 (약 46편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 약물 PK(탐수로신 3구획, 피나스테리드 3구획, 두타스테리드 3구획, 타달라필 3구획), 호르몬(테스토스테론·DHT 혈장·DHT 전립선 조직), AR 활성화 복합체, 구조/기능 엔드포인트(전립선 부피, cGMP, α1-AR 점유율, IPSS, 최대요속 Qmax, 잔뇨량 PVR, PSA, 전립선 염증 지수)
- **주요 치료 시나리오**: ① 경과 관찰(무치료), ② 탐수로신 0.4 mg QD, ③ 피나스테리드 5 mg QD, ④ 두타스테리드 0.5 mg QD, ⑤ 두타스테리드+탐수로신 병합, ⑥ 타달라필 5 mg QD
- **보정/근거**: PLESS(피나스테리드 4년), CombAT(병합 4년), MTOPS 시험 IPSS/Qmax/전립선 부피 변화 데이터 참고

## Shiny 대시보드 (Dashboard)
환자 프로파일(연령·전립선 부피·PSA·기저 IPSS), 약동학(각 약물 혈장 농도-시간), 호르몬/AR 탭(테스토스테론·DHT·AR 점유율), 임상 엔드포인트(IPSS·Qmax·PVR·전립선 부피), 치료 시나리오 비교(6개 오버레이), 바이오마커(PSA·DHT 전립선·염증 지수) 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("bph_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("bph_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg bph_qsp_model.dot -o bph_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [bph_references.md](bph_references.md) 참조 (약 46편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
