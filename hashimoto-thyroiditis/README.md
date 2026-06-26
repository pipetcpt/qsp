# 하시모토 갑상선염 (Hashimoto's Thyroiditis, HT) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![HT QSP Model](ht_qsp_model.png)](ht_qsp_model.svg)

## 개요 (Overview)

하시모토 갑상선염(만성 자가면역 갑상선염)은 항-TPO(갑상선과산화효소) 및 항-Tg(티로글로불린) 자가항체와 자기반응 T세포에 의한 갑상선 파괴로 점진적 기능저하가 발생하는 가장 흔한 자가면역 갑상선 질환입니다. 전 세계 유병률은 약 1~2%이며, 여성에서 7~10배 높은 빈도로 발생합니다. HLA-DR3/DR5, CTLA-4, PTPN22 등의 유전인자와 요오드 과잉·흡연·스트레스 등 환경 인자가 복합적으로 작용합니다. 자가면역 파괴가 진행되면 T4 생성이 감소하고 TSH가 상승하는 현성 갑상선기능저하증으로 이행하며, 레보티록신(LT4) 보충이 표준 치료입니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| T세포 자가면역 | Th1 우세 → IFN-γ·TNF-α → 갑상선 세포 아포토시스(Fas/FasL) | 갑상선 세포 파괴 |
| B세포·자가항체 | 항-TPO IgG → 보체 활성·ADCC → 추가 조직 손상 | 항-TPO/Tg 상승 |
| HPT 축 탈조절 | 갑상선 파괴 → T4↓ → 뇌하수체 TSH 분비 증가 → 갑상선종 | TSH 상승, 갑상선종 |
| 항산화 결핍 | 셀레늄 결핍 → GPx 활성 감소 → ROS 축적 → 갑상선 세포 손상 | 염증 지속·항-TPO 지속 |
| Treg 기능 이상 | FOXP3+ Treg 감소 → Th17/Th1 균형 파괴 | 자가면역 지속 |
| 갑상선 기능저하 | 갑상선 실질 소실 → T4·T3 생성 불충분 | 피로·서맥·부종·고지혈증 |

## 주요 약물 표적 (Drug Targets)

- **레보티록신 (Levothyroxine, LT4)**: T4 보충 → 말초 T3 전환; TSH 정상화가 치료 목표
- **셀레늄 (Selenium)**: 셀레노단백질(GPx·탈요오드효소) 보조 → 항-TPO 감소 효과 (SELENOIT 시험 근거)
- **리오티로닌 (Liothyronine, LiT3)**: T3 직접 보충; LT4 단독 치료 시 증상 지속 환자의 병합 요법
- **고용량 LT4 TSH 억제 요법**: 갑상선종 축소 목적으로 TSH를 낮게 유지
- **글루코코르티코이드**: 급성 통증성 갑상선염(아급성 하시모토) 일시적 사용
- **메티마졸 (탄화분기)**: 하시모토 독성증(Hashitoxicosis) 초기 갑상선기능항진기에 단기 사용

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [ht_qsp_model.dot](ht_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 100+ 노드 / 10 클러스터) |
| [ht_qsp_model.svg](ht_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [ht_qsp_model.png](ht_qsp_model.png) | PNG 이미지 (150 dpi) |
| [ht_mrgsolve_model.R](ht_mrgsolve_model.R) | mrgsolve ODE 모델 (약 19 구획 / 7개 치료 시나리오) |
| [ht_shiny_app.R](ht_shiny_app.R) | Shiny 대시보드 |
| [ht_references.md](ht_references.md) | 참고문헌 (약 44편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: HPT 축(TRH·TSH), 갑상선 내 T4·T3 합성, 혈장 T4·T3·fT4·fT3·rT3, 조직 T4·T3, 면역(Th1·Treg·B세포·항-TPO 항체·항-Tg 항체·갑상선 손상 지수), 약물 PK(LT4 장·중심·말초, LiT3, 셀레늄) 구획 포함
- **주요 치료 시나리오**: ① 무치료 기저선(진행성 기능저하), ② LT4 100 mcg/day, ③ 셀레늄 200 mcg/day, ④ LT4+셀레늄 병합, ⑤ LT4+LiT3 병합 요법, ⑥ 고용량 LT4(TSH 억제), ⑦ LT4 75 mcg + 셀레늄 200 mcg
- **보정/근거**: SELENOIT 시험(Ventura 2017), Gartner et al. (JCEM 2002) 셀레늄 데이터, Celi et al. (JCEM 2011) T4+T3 병합 요법 데이터를 기반으로 TSH·fT4 정상화 곡선을 정성적으로 보정

## Shiny 대시보드 (Dashboard)

환자 프로파일(기저 TSH, 항-TPO 역가, 갑상선 손상 정도, 셀레늄 결핍 여부) 탭, 갑상선 호르몬 PK·PD 동역학, 면역 지표(항-TPO·항-Tg 항체, T세포 균형), 임상 엔드포인트(TSH 정상화, fT4 수준), 7개 치료 시나리오 비교, 바이오마커(TSH·항-TPO·셀레늄) 탭으로 구성됩니다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("ht_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("ht_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg ht_qsp_model.dot -o ht_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [ht_references.md](ht_references.md) 참조 (약 44편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
