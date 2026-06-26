# 굿파스처 증후군 (Goodpasture Syndrome, GPS) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 신장·비뇨

[![GPS QSP Model](gps_qsp_model.png)](gps_qsp_model.svg)

## 개요 (Overview)

굿파스처 증후군(항-GBM 질환)은 제IV형 콜라겐 α3 사슬 NC1 도메인에 대한 자가항체가 사구체기저막(GBM)과 폐포기저막에 결합하여 보체 활성 및 호중구 침윤을 유발하는 희귀 자가면역 질환입니다. 연간 발생률은 100만 명당 약 0.5~1명으로 매우 드물며, 급속진행성 사구체신염과 폐포출혈이 동시에 나타나는 폐신증후군이 특징입니다. HLA-DR15 대립유전자가 주요 유전적 위험인자이며, 흡연과 탄화수소계 용매 노출이 폐 출혈의 촉발 인자로 알려져 있습니다. 혈장교환술로 순환 항체를 제거하고 사이클로포스파마이드(CY) 및 고용량 코르티코스테로이드로 면역을 억제하는 것이 표준 치료이며, 난치성 또는 재발 사례에서 리툭시맙이 사용됩니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 항-GBM 항체 형성 | HLA-DR15 제한 항원 제시 → B세포 활성 → 형질세포 → IgG 항α3(IV) 항체 | 사구체 선형 IgG 침착 |
| 보체 활성 | 항체-GBM 결합 → 고전 경로 C3 활성 → C5a·MAC 생성 | 호중구 모집, 내피 손상 |
| 호중구 침윤 | C5a 수용체 자극 → 산화적 폭발 → 프로테아제 분비 | 사구체 초승달 형성, GFR 저하 |
| Th1 면역 반응 | CD4 T세포 → IFN-γ · TNF-α → 대식세포 활성 | 염증 지속 및 조직 파괴 |
| 폐포 손상 | 흡연/흡입 노출 → 폐포기저막 항원 노출 → Ab 결합 | 폐포출혈·객혈 |
| 조절 T세포 결핍 | Treg 기능 저하 → 자가내성 파탄 | 자가면역 지속 |

## 주요 약물 표적 (Drug Targets)

- **혈장교환술 (Plasmapheresis)**: 순환 항-GBM 항체 직접 제거; 신경한 초기 치료
- **사이클로포스파마이드 (Cyclophosphamide, CY)**: B세포·형질세포 알킬화 억제; 항체 생성 감소
- **고용량 코르티코스테로이드 (Prednisolone)**: 염증 억제, 호중구 이동 차단, 보체 하류 효과 완화
- **리툭시맙 (Rituximab, RTX)**: 항-CD20 → B세포 고갈; 재발·난치성 항-GBM 질환에 적용
- **아바코판 (Avacopan)**: C5aR1 차단제; 병존 ANCA 혈관염 시 신독성 감소 기대
- **면역흡착 (Immunoadsorption)**: Protein A·G 컬럼으로 항-GBM IgG 선택적 제거

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [gps_qsp_model.dot](gps_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 100+ 노드 / 10 클러스터) |
| [gps_qsp_model.svg](gps_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [gps_qsp_model.png](gps_qsp_model.png) | PNG 이미지 (150 dpi) |
| [gps_mrgsolve_model.R](gps_mrgsolve_model.R) | mrgsolve ODE 모델 (약 20 구획 / 6개 치료 시나리오) |
| [gps_shiny_app.R](gps_shiny_app.R) | Shiny 대시보드 |
| [gps_references.md](gps_references.md) | 참고문헌 (약 59편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK(CY·활성대사체, 프레드니솔론, RTX 2구획, 아바코판) + 질환 PD(항-GBM 항체, B세포, 형질세포, C5a, 신장 호중구, GBM 손상, GFR, 폐 손상, DLCO, 단백뇨, 혈뇨, CRP, Treg) 구획 포함
- **주요 치료 시나리오**: ① 무치료 기저선, ② 혈장교환+CY+스테로이드(표준), ③ RTX 단독, ④ 혈장교환+RTX+스테로이드, ⑤ 아바코판 병합, ⑥ 조기 vs. 지연 치료 비교
- **보정/근거**: Levy et al. (Ann Intern Med 2001), McAdoo & Pusey (Clin J Am Soc Nephrol 2017) 등 항-GBM 질환 코호트 데이터를 기반으로 GFR 회복 곡선과 항체 역가 감소 동역학을 정성적으로 보정

## Shiny 대시보드 (Dashboard)

환자 프로파일(초기 GFR, 항-GBM 역가, 폐 침범 여부) 설정 탭, 약물 PK 시뮬레이션(RTX·CY 혈중 농도), 주요 PD 지표(항체 역가, 보체 C5a, 호중구), 신장 임상 엔드포인트(GFR 추이, 단백뇨, 혈뇨), 폐 임상 엔드포인트(DLCO, 폐출혈 지표), 6개 치료 시나리오 비교, 바이오마커 탭으로 구성됩니다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("gps_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("gps_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg gps_qsp_model.dot -o gps_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [gps_references.md](gps_references.md) 참조 (약 59편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
