# 원형 탈모증 (Alopecia Areata, AA) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 피부

[![AA QSP Model](aa_qsp_model.png)](aa_qsp_model.svg)

## 개요 (Overview)

원형 탈모증(AA)은 모낭을 표적으로 하는 T세포 매개 자가면역 질환으로, 전 세계 인구의 약 2%에서 평생 1회 이상 발생하는 흔한 탈모증이다. 핵심 발병기전은 모낭이 정상적으로 유지하는 '면역 특권(immune privilege)' 붕괴로, MHC-I 발현 억제 기전이 소실되면 CD8⁺ NKG2D⁺ T세포와 NK세포가 모낭 세포를 공격한다. IFN-γ·IL-15가 JAK1/2-STAT1/STAT5 신호를 통해 이 면역 공격을 증폭시키며, CXCL10을 통한 추가 T세포 모집이 지속된다. 두피 원형 탈모에서 전두 탈모(전체 두피), 전신 탈모(전신 체모)까지 중증도가 다양하다. JAK 억제제(바리시티닙, 리틀레시티닙, 데우크라바시티닙)가 최근 승인되어 임상에 도입됐다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 면역 특권 붕괴 | MHC-I 발현 억제 해제 → 모낭 항원 노출 | CD8⁺ T세포 공격 개시 |
| IL-15/NKG2D 축 | IL-15 → NK세포·CD8⁺ T세포 NKG2D 발현 증가 | 모낭 세포독성 |
| JAK-STAT 신호 | IFN-γ·IL-15 → JAK1/2 → STAT1/3/5 | 염증 유전자 발현, 면역 공격 지속 |
| IFN-γ/CXCL10 피드백 | CXCL10 분비 → CXCR3⁺ T세포 추가 모집 | 염증 증폭 순환 |
| Treg 기능 장애 | Treg 감소·기능 손상 → 면역 억제 실패 | 자기 관용 소실 |
| 모낭 주기 | 면역 공격 → 성장기(anagen) 조기 종료 | 모발 탈락·빈모 |
| JAK3/TYK2 의존 경로 | IL-4·IL-13(아토피 동반) → Th2 성분 | 아토피성 AA에서 중요 |

## 주요 약물 표적 (Drug Targets)

- **바리시티닙 (Baricitinib, JAK1/2 억제제)**: BRAVE-AA1/2 임상시험에서 중등도-중증 AA의 SALT 점수 유의 개선; 2022년 FDA 승인
- **리틀레시티닙 (Ritlecitinib, JAK3/TEC 억제제)**: ALLEGRO 임상시험; 12세 이상 적응증
- **루소리티닙 (Ruxolitinib, JAK1/2 억제제)**: 외용제 및 경구제 형태; IFN-γ/CXCL10 하향 조절
- **두필루맙 (Dupilumab, 항-IL-4Rα)**: 아토피 동반 AA에서 일부 효과; Th2 경로 차단
- **코르티코스테로이드**: 국소·병변내·전신; 일시적 효과, 장기 치료 제한

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [aa_qsp_model.dot](aa_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 164 노드 / 11 클러스터) |
| [aa_qsp_model.svg](aa_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [aa_qsp_model.png](aa_qsp_model.png) | PNG 이미지 (150 dpi) |
| [aa_mrgsolve_model.R](aa_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 다수 치료 시나리오) |
| [aa_shiny_app.R](aa_shiny_app.R) | Shiny 대시보드 |
| [aa_references.md](aa_references.md) | 참고문헌 (약 55편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 바리시티닙·루소리티닙 경구 PK(2구획) + 두필루맙 SC PK + JAK3 결합 구획 + NKG2DL·NK세포·CD8 나이브·CD8 이펙터·Treg·IFN-γ·IL-15·CXCL10·pSTAT1·pSTAT5·면역 특권 지수·성장기 모발·모발 밀도·전신 염증 PD 구획
- **주요 치료 시나리오**: ① 무치료(자연 경과), ② 바리시티닙 2 mg/day, ③ 바리시티닙 4 mg/day, ④ 루소리티닙, ⑤ 두필루맙(아토피 동반), ⑥ 스테로이드 + JAK 억제제 병용
- **보정/근거**: BRAVE-AA1/2(바리시티닙), Xing et al. Nat Med 2014 기초 면역 파라미터 참조

## Shiny 대시보드 (Dashboard)

환자 프로파일(SALT 점수·탈모 범위·아토피 동반 여부) 탭, JAK 억제제 PK 및 표적 점유율 탭, 면역세포 동태(CD8·Treg·NK) 탭, 모발 밀도 및 SALT 점수 변화 탭, 치료 시나리오 비교 탭, 사이토카인 바이오마커(IFN-γ·CXCL10) 탭으로 구성된다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("aa_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("aa_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg aa_qsp_model.dot -o aa_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [aa_references.md](aa_references.md) 참조 (약 55편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
