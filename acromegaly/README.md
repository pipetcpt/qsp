# 말단비대증 (Acromegaly, ACRO) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![ACRO QSP Model](acro_qsp_model.png)](acro_qsp_model.svg)

## 개요 (Overview)

말단비대증은 성장호르몬(GH)을 과다 분비하는 뇌하수체 선종(GH-종분비 선종)에 의해 발생하며, 인구 100만 명당 약 60명의 유병률을 보이는 희귀 내분비 질환이다. GH 과잉 분비는 간에서 IGF-1 생성을 지속적으로 촉진하여 말단 연조직·골·내장 과성장, 심혈관 합병증, 당대사 이상, 관절병증 등을 유발한다. 진단 지연이 평균 8–10년으로 길고, 치료받지 않으면 표준화 사망비가 2–3배 증가한다. 수술적 절제가 1차 치료이며, 잔존 종양에는 소마토스타틴 유사체(SSA)·페그비소만트·도파민 작용제가 약물 표적이다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| GH 과잉 분비 | 뇌하수체 선종 자율 GH 분비, GHRH 과자극 | GH 기저치 및 포도당부하 최저치 상승 |
| IGF-1 축 | JAK2-STAT5b 활성화 → 간 IGF-1 합성 증가 | 연조직·장기 비대, 뼈 과성장 |
| SSTR 신호 | SSTR2/5에 SSA 결합 → cAMP 억제, GH 분비 감소 | 혈중 GH/IGF-1 정상화 |
| GHR 차단 | 페그비소만트(GHR 길항제) → IGF-1 생성 차단 | GH 상승에도 IGF-1 정상화 |
| 심혈관 리모델링 | IGF-1 수용체 → 좌심실 비대, 부정맥 | 심부전, 심혈관 사망 위험 |
| 대사 이상 | GH의 인슐린 길항 효과 | 당뇨, 고혈압, 고지혈증 |
| 종양 성장 | MAPK/PI3K/AKT → 세포 증식 | 종양 크기·침습성 결정 |

## 주요 약물 표적 (Drug Targets)

- **소마토스타틴 유사체 (Octreotide LAR, Lanreotide)**: SSTR2/5 작용 → GH·IGF-1 감소, 종양 크기 소폭 축소
- **페그비소만트 (Pegvisomant)**: GHR 길항제 → IGF-1 정상화율 약 90%; GH 상승에 의한 종양 크기 모니터링 필요
- **도파민 작용제 (Cabergoline)**: D2 수용체 작용, 단독 효과 제한적이나 SSA와 병용 시 IGF-1 추가 억제
- **수술 (경접형동 수술)**: 미세선종 치유율 약 80%; 거대선종 완치율 낮음
- **방사선 치료**: 서서히 GH 저하; 뇌하수체 기능 저하증 위험

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [acro_qsp_model.dot](acro_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 147 노드 / 10 클러스터) |
| [acro_qsp_model.svg](acro_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [acro_qsp_model.png](acro_qsp_model.png) | PNG 이미지 (150 dpi) |
| [acro_mrgsolve_model.R](acro_mrgsolve_model.R) | mrgsolve ODE 모델 (약 17 구획 / 다수 치료 시나리오) |
| [acro_shiny_app.R](acro_shiny_app.R) | Shiny 대시보드 |
| [acro_references.md](acro_references.md) | 참고문헌 (약 39편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: SSA(Octreotide) 2구획 PK + 페그비소만트 2구획 PK + SSTR 결합·GHR 차단 PD 구획 + GH 선종 분비·GH 혈장·STAT5b 활성·IGF-1 간·IGF-1 혈장·종양 부피·좌심실비대·혈당·관절 점수 등 PD 구획
- **주요 치료 시나리오**: ① 무치료, ② Octreotide LAR 단독, ③ 페그비소만트 단독, ④ SSA + 페그비소만트 병용, ⑤ 수술 후 잔존 선종 SSA, ⑥ 도파민 작용제 추가 병용
- **보정/근거**: ACROSTUDY 레지스트리(페그비소만트), PRIMARYS 임상시험(Lanreotide), Colao et al. 메타분석 파라미터 참조

## Shiny 대시보드 (Dashboard)

환자 프로파일(성별·나이·종양 크기) 탭, GH/IGF-1 PK-PD 동태 탭, 선종 부피 변화 탭, 심혈관·대사 합병증 지표 탭, 치료 시나리오 비교 탭, 바이오마커(SSTR 점유율·IGF-1 정상화율) 탭으로 구성된다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("acro_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("acro_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg acro_qsp_model.dot -o acro_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [acro_references.md](acro_references.md) 참조 (약 39편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
