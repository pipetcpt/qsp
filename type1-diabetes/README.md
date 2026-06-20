# 제1형 당뇨병 (Type 1 Diabetes, T1DM) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![T1DM QSP Model](t1dm_qsp_model.png)](t1dm_qsp_model.svg)

## 개요 (Overview)

제1형 당뇨병(T1DM)은 자가반응성 CD8+ T세포(CTL)가 췌장 베타세포를 선택적으로 파괴하여 절대적 인슐린 결핍을 초래하는 자가면역 질환이다. 전 세계 유병률은 약 850만 명(2022년)이며 연간 약 18만 명의 소아·청소년에서 새로 진단된다. 베타세포 소실은 진단 전 수년에 걸쳐 서서히 진행되며, 발현 시점에는 이미 베타세포의 70~90%가 소실된 상태이다. 인슐린 집중 요법이 혈당 조절의 근간이며, CD3 항체 테플리주맙이 최근 발병 지연 효능으로 최초 승인된 면역조절 치료제이다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 자가면역 베타세포 공격 | GAD·IA-2·ZnT8 자가항체, CD8+ CTL의 MHC-I 매개 베타세포 용해 | 베타세포 질량 감소, C-펩타이드 소실 |
| Treg 기능 부전 | Foxp3+ Treg 활성 저하, 면역 관용 붕괴 | CTL 억제 실패, 염증 지속 |
| 인슐린 결핍 | 베타세포 소실 → 기저·식후 인슐린 분비 중단 | 혈당 조절 불능, 케톤산증 위험 |
| 포도당-인슐린 동역학 | 간 포도당 생산(EGP) 조절 불능, 말초 포도당 이용 감소 | 공복·식후 고혈당 |
| 글루카곤 역조절 이상 | 인슐린 억제 부재 → 글루카곤 과분비, EGP 증가 | 저혈당 후 회복 불량 |
| HbA1c 동역학 | 평균 혈당 → 적혈구 당화 축적 | 만성 합병증(망막·신·신경·혈관병증) 위험 |

## 주요 약물 표적 (Drug Targets)

- **기저 인슐린** — 인슐린 데글루덱·글라진: 공복 혈당 조절, 야간 저혈당 예방
- **식사 인슐린(볼루스)** — 초속효성 인슐린 아스파트·리스프로: 식후 고혈당 억제
- **인슐린 펌프(CSII) + 하이브리드 폐쇄 루프(HCL)** — 연속혈당측정(CGM) 연동 자동 인슐린 주입
- **CD3 항체 — 테플리주맙(Teplizumab)**: 자가반응 T세포 Teff 재프로그래밍, 발병 지연 (TrialNet TN10)
- **GAD 백신·항-CD20(리툭시맙)** — 연구 단계: 베타세포 자가면역 경로 차단

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [t1dm_qsp_model.dot](t1dm_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 236 노드 / 12 클러스터) |
| [t1dm_qsp_model.svg](t1dm_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [t1dm_qsp_model.png](t1dm_qsp_model.png) | PNG 이미지 (150 dpi) |
| [t1dm_mrgsolve_model.R](t1dm_mrgsolve_model.R) | mrgsolve ODE 모델 (약 20 구획 / 6 치료 시나리오) |
| [t1dm_shiny_app.R](t1dm_shiny_app.R) | Shiny 대시보드 |
| [t1dm_references.md](t1dm_references.md) | 참고문헌 (약 52편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 베타세포·CTL·Treg(면역) 3구획 + C-펩타이드 1구획 + 포도당 동역학(Gp·Gt) 2구획 + 인슐린 PK(SC1·SC2·중심·말초·원격작용) 5구획 + 글루카곤 1구획 + HbA1c 1구획 + 테플리주맙 PK 1구획 + CGM·APC 2구획 + 식사 포도당(Qsto1·Qsto2·Qgut) 3구획 + 평균 혈당 1구획으로 총 약 20개 구획
- **주요 치료 시나리오**: (1) 미치료 T1DM, (2) MDI(기저+볼루스 인슐린), (3) CSII(인슐린 펌프), (4) 하이브리드 폐쇄 루프(HCL), (5) 테플리주맙(2단계 발병 예방), (6) 테플리주맙 + MDI 병용
- **보정/근거**: Bergman minimal model 기반 포도당-인슐린 동역학; TrialNet TN10 임상시험(테플리주맙 2년 지연 효과); Dalla Man 식사 흡수 모델; DCCT 연구 HbA1c-합병증 관계 참조

## Shiny 대시보드 (Dashboard)

환자 프로파일(체중·연령·베타세포 잔존 기능·HbA1c 기저치 설정) · PK 시각화(인슐린·테플리주맙 혈중 농도-시간 곡선) · 혈당 PD 지표(공복·식후 혈당·CGM 시뮬레이션) · 임상 엔드포인트(HbA1c·저혈당 빈도·C-펩타이드 보존) · 치료 시나리오 비교(베타세포 질량·HbA1c 궤적) · 바이오마커 패널(CTL·Treg·C-펩타이드·글루카곤) 등 6개 이상의 탭으로 구성된다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("t1dm_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("t1dm_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg t1dm_qsp_model.dot -o t1dm_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [t1dm_references.md](t1dm_references.md) 참조 (약 52편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
