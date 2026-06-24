# 호산구 육아종증 다발혈관염 (Eosinophilic GPA, EGPA) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 혈관염

[![EGPA QSP Model](egpa_qsp_model.png)](egpa_qsp_model.svg)

## 개요 (Overview)

호산구 육아종증 다발혈관염(EGPA, 구 Churg-Strauss 증후군)은 천식·호산구증가증·소혈관 육아종성 염증을 3대 특징으로 하는 희귀 ANCA 연관 혈관염으로, 연간 발생률은 100만 명당 약 0.5~3명입니다. 병인은 IL-5 매개 호산구 과증식과 Th2 편향 면역 반응, 일부 환자에서 ANCA(MPO) 양성 소혈관염이 복합됩니다. 활동성 장기 손상(심장·신경·신장)이 없는 경증은 스테로이드 단독으로, 중등증 이상 또는 스테로이드 의존성은 메폴리주맙(항IL-5)이 표준 치료입니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| IL-5 경로 | IL-5 → 호산구 골수 성숙·혈중 유입·조직 침윤 | 혈중 호산구 수 급증 |
| Th2 편향 면역 | IL-4, IL-13 과발현 → IgE 생성, 점막 염증 | 천식 악화, 비용종 |
| 호산구 조직 손상 | ECP, MBP, EPO 방출 → 혈관·심근·신경 독성 | 심근염, 말초신경병 |
| ANCA(MPO) 경로 | anti-MPO ANCA → 호중구 활성화 → 소혈관 염증 | 사구체신염, 폐포출혈 |
| 육아종 형성 | 대식세포·거대세포 → 혈관 주위 육아종 | 장기 기능 손상 |
| 천식 기도 염증 | 호산구·비만세포·IgE → 기관지 과반응성 | FEV1 저하, 천식 발작 |
| 심장 침범 | 호산구성 심근염·심내막염(Löffler) | 심부전, 부정맥 위험 |

## 주요 약물 표적 (Drug Targets)

- **프레드니솔론(전신 코르티코스테로이드)**: 광범위 면역억제, 호산구 세포사 유도
- **메폴리주맙(Mepolizumab)**: 항IL-5 단일클론항체 → 호산구 감소, MIRRA 시험에서 관해율 향상
- **벤랄리주맙(Benralizumab)**: 항IL-5Rα → 호산구·호염기구 직접 고갈 (NK세포 ADCC 기전)
- **리툭시맙**: ANCA 양성 중증 EGPA에서 B세포 고갈
- **아자티오프린/메토트렉세이트**: 스테로이드 감량 유지요법

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [egpa_qsp_model.dot](egpa_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 166 노드 / 10 클러스터) |
| [egpa_qsp_model.svg](egpa_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [egpa_qsp_model.png](egpa_qsp_model.png) | PNG 이미지 (150 dpi) |
| [egpa_mrgsolve_model.R](egpa_mrgsolve_model.R) | mrgsolve ODE 모델 (약 32 구획 / 6 치료 시나리오) |
| [egpa_shiny_app.R](egpa_shiny_app.R) | Shiny 대시보드 |
| [egpa_references.md](egpa_references.md) | 참고문헌 (약 61편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK 구획(프레드니솔론·메폴리주맙·벤랄리주맙 depot/중심/말초) + PD 구획(IL-5 free/bound, 호산구 전구체·혈중·조직, Th2 사이토카인, ANCA, 혈관 손상 지수, 천식 FEV1, 심장 침범 지수, IgE, 육아종 부담)
- **주요 치료 시나리오**: ① 무치료(자연 경과), ② 프레드니솔론 단독, ③ 메폴리주맙 + 프레드니솔론(MIRRA 시험 기반), ④ 벤랄리주맙 + 프레드니솔론, ⑤ 리툭시맙 + 프레드니솔론(ANCA 양성 중증), ⑥ 스테로이드 감량 유지요법(아자티오프린)
- **보정/근거**: MIRRA 시험(메폴리주맙), MANDARA 시험(벤랄리주맙), Guillevin 등 EGPA 코호트 데이터 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① 환자 프로파일(ANCA 상태·기저 호산구 수·FEV1·장기 침범 설정), ② PK 탭(생물학적 제제 및 스테로이드 혈중 농도), ③ 면역 PD 탭(IL-5·호산구·Th2 사이토카인 추이), ④ 임상 엔드포인트(관해율·재발 횟수·FEV1·심장 지표), ⑤ 시나리오 비교(6개 치료 전략 동시 비교), ⑥ 바이오마커(IgE·ANCA 역가·혈관 손상 지수 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("egpa_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("egpa_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg egpa_qsp_model.dot -o egpa_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [egpa_references.md](egpa_references.md) 참조 (약 61편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
