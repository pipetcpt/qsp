# 거대세포 동맥염 (Giant Cell Arteritis, GCA) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 혈관염

[![GCA QSP Model](gca_qsp_model.png)](gca_qsp_model.svg)

## 개요 (Overview)

거대세포 동맥염(GCA)은 50세 이상에서 발생하는 대혈관(측두동맥, 대동맥 가지)의 육아종성 혈관염으로, 50세 이상 인구 10만 명당 약 10~30명의 발생률을 보입니다. 측두동맥 압통·두통·턱 파행·실명(전방허혈성 시신경병)이 주요 증상이며, 즉각적인 스테로이드 치료가 실명 예방에 필수적입니다. 병태생리는 혈관벽 수지상세포가 항원을 제시하여 Th1/Th17 세포를 활성화하고, 대식세포·거대세포가 혈관벽에 침윤하여 IL-6 과발현과 VEGF 매개 혈관 신생이 일어납니다. 토실리주맙(항IL-6R)이 스테로이드 감량에 효과적입니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 수지상세포 항원 제시 | 혈관벽 DC → CD4+ T세포 활성화 | Th1/Th17 분화 촉진 |
| Th17/IL-17 경로 | IL-17A → 호중구 유입·혈관 염증 | ESR/CRP 상승, 혈관 손상 |
| IL-6 과발현 | 대식세포 유래 IL-6 → 급성기 반응 | ESR·CRP 상승, 전신 증상 |
| 대식세포·거대세포 육아종 | 혈관 외막·중막 침윤 → 내막 증식 | 혈관 폐쇄, 허혈 |
| VEGF 매개 혈관 신생 | Th17 → VEGF → 신혈관 형성 | 내막 증식 가속 |
| 코르티코스테로이드 반응 | GR → NF-κB 억제 → 사이토카인 감소 | 증상 빠른 개선 |
| 골밀도 손실 | 장기 스테로이드 → BMD 감소 | 골다공증·골절 위험 |

## 주요 약물 표적 (Drug Targets)

- **고용량 스테로이드** (프레드니솔론 40~60 mg/일): 즉각적인 염증 억제, 실명 예방 — 1차 치료
- **토실리주맙(Tocilizumab)**: 항IL-6R → IL-6 신호 차단 → 스테로이드 감량 가능(GiACTA 시험)
- **메토트렉세이트**: 스테로이드 감량 보조, 재발 예방(소규모 RCT)
- **아바타셉트**: CTLA4-Ig → T세포 공동자극 차단 → 유지 치료 연구 중
- **자누브루티닙(BTK 억제제)**: 새로운 치료 표적 연구 중

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [gca_qsp_model.dot](gca_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 188 노드 / 10 클러스터) |
| [gca_qsp_model.svg](gca_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [gca_qsp_model.png](gca_qsp_model.png) | PNG 이미지 (150 dpi) |
| [gca_mrgsolve_model.R](gca_mrgsolve_model.R) | mrgsolve ODE 모델 (약 23 구획 / 5 치료 시나리오) |
| [gca_shiny_app.R](gca_shiny_app.R) | Shiny 대시보드 |
| [gca_references.md](gca_references.md) | 참고문헌 (약 42편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK 구획(프레드니솔론 GI/혈장/조직, 토실리주맙 SC depot/중심/말초, sIL-6R·IL-6 free/bound TMDD 구획) + PD 구획(IL-6 유리형/결합형, CRP, ESR, 수지상세포 활성화, 대식세포 활성화, Th17, VEGF, 골밀도, 누적 스테로이드 용량)
- **주요 치료 시나리오**: ① GC 단독 완속 감량(GiACTA 대조군), ② GC 단독 급속 감량(GiACTA 위약군), ③ 토실리주맙 IV 8 mg/kg q4w + GC 완속 감량(GiACTA arm 1), ④ 토실리주맙 SC 162 mg qw + GC 완속 감량(GiACTA arm 2), ⑤ 토실리주맙 SC 162 mg q2w + GC 완속 감량(GiACTA arm 3)
- **보정/근거**: GiACTA 시험(Stone et al. NEJM 2017), 토실리주맙 IL-6R PK/PD 문헌, ESR/CRP 급성기 반응 모델 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① 환자 프로파일(기저 ESR·CRP·VDI 손상 지수·골밀도 설정), ② PK 탭(프레드니솔론·토실리주맙 혈중 농도 및 TMDD), ③ 면역 PD 탭(IL-6·sIL-6R·CRP·ESR·대식세포 활성 추이), ④ 임상 엔드포인트(관해 유지율·재발 횟수·스테로이드 누적 용량), ⑤ 시나리오 비교(5개 GiACTA 치료 전략 동시 비교), ⑥ 바이오마커(Th17·VEGF·골밀도·누적 GC 용량 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("gca_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("gca_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg gca_qsp_model.dot -o gca_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [gca_references.md](gca_references.md) 참조 (약 42편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
