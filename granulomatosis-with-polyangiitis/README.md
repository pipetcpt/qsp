# 육아종증 다발혈관염 (GPA) (Granulomatosis with Polyangiitis, GPA) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 혈관염

[![GPA QSP Model](gpa_qsp_model.png)](gpa_qsp_model.svg)

## 개요 (Overview)

육아종증 다발혈관염(구 베게너 육아종증)은 PR3-ANCA(프로테이나제 3에 대한 항호중구세포질항체)에 의해 호중구가 활성화되고 소혈관에 괴사성 육아종성 염증이 발생하는 전신 혈관염입니다. 연간 발생률은 100만 명당 약 10~12명이며, 상기도(비부비동), 폐, 신장이 주된 침범 장기입니다. PR3-ANCA가 호중구 막 표면 PR3에 결합해 중성구를 프라이밍하고, 보체 C5a가 이를 더욱 증폭시켜 NETs(호중구 세포외 덫)와 혈관 내피 손상을 일으킵니다. 리툭시맙과 사이클로포스파마이드가 관해 유도에 사용되며, 신규 C5aR 차단제 아바코판이 스테로이드 대체 치료로 승인되었습니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| PR3-ANCA 형성 | B세포→형질세포 → PR3-ANCA IgG 분비; LAMP-2 교차반응 | ANCA 역가 상승 |
| 호중구 프라이밍·활성화 | TNF-α/LPS 프라이밍 → ANCA 결합 → Fc수용체·FcγRIIIb 활성화 → 산화적 폭발 | 혈관 내피 손상 |
| NETs 형성 | 활성 호중구 → 크로마틴·엘라스타제·MPO 방출 → 혈전 형성 | 신사구체·폐 손상 |
| 보체-C5a 축 | ANCA-결합 호중구 → 대체경로 C3b/C5 분열 → C5a 방출 | C5aR1 활성·호중구 증폭 |
| 육아종 형성 | Th1 CD4 T세포·대식세포 → IL-12·IFN-γ → 다핵 거대세포·상피양 세포 | 비강·폐 육아종, 연골 파괴 |
| 신장 침범 | 반월체형 사구체신염(ANCA-GN) → GFR 급감 | 급성 신부전, ESRD 진행 |
| 리툭시맙 B세포 고갈 | 항-CD20 → B세포·기억 B세포 제거 → ANCA 감소 | 관해 유지, 재발 감소 |

## 주요 약물 표적 (Drug Targets)

- **리툭시맙 (Rituximab, RTX)**: 항-CD20 B세포 고갈; RAVE 시험에서 CYC 대비 비열등(관해 유도), MAINRITSAN 시험에서 유지 요법 우월성 입증
- **사이클로포스파마이드 (Cyclophosphamide, CYC)**: 알킬화제; IV 펄스(CYCLOPS 시험) 및 경구 요법으로 관해 유도
- **아바코판 (Avacopan)**: 경구 C5aR1 선택적 차단제; ADVOCATE 시험에서 고용량 스테로이드 대체 달성
- **아자티오프린/메토트렉세이트 (Azathioprine/MTX)**: 스테로이드 감량 후 유지 요법
- **고용량 글루코코르티코이드**: 관해 유도 보조; 아바코판으로 용량 감소 추세
- **마이코페놀레이트 모페틸 (MMF)**: 유지 요법 대안 (신독성 우려 환자)

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [gpa_qsp_model.dot](gpa_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 110+ 노드 / 13 클러스터) |
| [gpa_qsp_model.svg](gpa_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [gpa_qsp_model.png](gpa_qsp_model.png) | PNG 이미지 (150 dpi) |
| [gpa_mrgsolve_model.R](gpa_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 5개 치료 시나리오) |
| [gpa_shiny_app.R](gpa_shiny_app.R) | Shiny 대시보드 |
| [gpa_references.md](gpa_references.md) | 참고문헌 (약 52편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK(RTX 2구획+TMDD 결합, CYC 장·혈장·활성대사체, GC 장·혈장, 아바코판 장·혈장) + 질환 PD(순진 B세포, 기억 B세포, 장기생존 형질세포, ANCA, C5a, 휴지·활성 호중구, NETs, 내피 손상, 육아종 지수, GFR, BVAS 질환 활성도 점수) 구획 포함
- **주요 치료 시나리오**: ① 무치료 기저선, ② RTX+GC 관해 유도(RAVE 프로토콜), ③ CYC+GC 관해 유도, ④ 아바코판+RTX(ADVOCATE 설계), ⑤ 재발 후 RTX 재유도
- **보정/근거**: RAVE(Stone 2010), RITUXVAS(Jones 2010), ADVOCATE(Jayne 2021), MAINRITSAN3(Charles 2023) 임상시험 데이터를 기반으로 BVAS 관해율·재발률 곡선을 정성적으로 보정

## Shiny 대시보드 (Dashboard)

환자 프로파일(기저 BVAS, ANCA 역가, GFR, 장기 침범 선택) 탭, 약물 PK 시뮬레이션(RTX·CYC 혈중 농도), 면역 PD 지표(B세포, ANCA, C5a, NETs), 신장 임상 엔드포인트(GFR, 단백뇨), 폐·비강 지표, 5개 치료 시나리오 비교(BVAS 궤적), 바이오마커(PR3-ANCA·MPO-ANCA·CRP) 탭으로 구성됩니다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("gpa_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("gpa_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg gpa_qsp_model.dot -o gpa_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [gpa_references.md](gpa_references.md) 참조 (약 52편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
