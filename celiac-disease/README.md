# 셀리악병 (Celiac Disease, CD) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도

[![CD QSP Model](cd_qsp_model.png)](cd_qsp_model.svg)

## 개요 (Overview)
셀리악병은 HLA-DQ2/DQ8 양성 유전적 소인을 가진 개인에서 글루텐(밀·보리·호밀의 글리아딘/글루테닌)이 소장 점막의 자가면역 반응을 유발하여 융모 위축과 흡수 불량을 초래하는 전신 자가면역 질환입니다. 서양에서 유병률은 인구의 약 1%(진단율은 훨씬 낮음)에 달하며, 조직 트랜스글루타미나제(tTG)에 의한 탈아미드 글루텐 펩타이드(DGP)가 HLA-DQ에 제시되어 CD4 Th1/Th17 반응을 촉발합니다. 완전한 글루텐 제거식(GFD)이 유일한 승인 치료이며, 장벽 투과성 차단제(라라조타이드), tTG2 억제제(ZED1227), 항-IL-15 항체(AMG714) 등 신약이 개발 중입니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 장 투과성 증가 | 조누린 상향 조절, 밀착연접 손상 → 글리아딘 통과 | 고유판 글루텐 유입 |
| tTG2 탈아미드화 | tTG2 → 글리아딘 Gln→Glu → 항원성 극대화 | HLA-DQ2/8 결합 효율 증가 |
| HLA-DQ2/8 제시 | DC → CD4+ Th1/Th17 활성화 | IFN-γ·IL-17 기반 점막 손상 |
| IL-15 매개 IEL 활성화 | IL-15 → NK-T 유사 IEL → 상피세포 직접 세포독성 | 융모 상피 파괴 |
| B세포·자가항체 생성 | tTG2 IgA 자가항체, 항DGP IgA/IgG | 진단 바이오마커, 추가 조직 손상 |
| 영양소 흡수 불량 | 융모 위축 → 흡수 면적 감소 | 철·엽산·칼슘·비타민D 결핍 |

## 주요 약물 표적 (Drug Targets)
- **글루텐 제거식 (GFD)**: 항원 제거 → 면역 반응 소실 — 유일한 표준 치료
- **라라조타이드 아세테이트 (장벽 보호제)**: 밀착연접 강화 → 장 투과성 감소 — 글루텐 노출 완충
- **ZED1227 (tTG2 억제제)**: tTG2 경쟁적 억제 → DGP 생성 차단 — 2상 시험에서 융모 위축 개선
- **AMG714 (항-IL-15)**: IL-15 중화 → IEL 활성화 억제 — 불응성 셀리악병(RCD) 대상

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [cd_qsp_model.dot](cd_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 133 노드 / 10 클러스터) |
| [cd_qsp_model.svg](cd_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [cd_qsp_model.png](cd_qsp_model.png) | PNG 이미지 (150 dpi) |
| [cd_mrgsolve_model.R](cd_mrgsolve_model.R) | mrgsolve ODE 모델 (약 20 구획 / 6개 치료 시나리오) |
| [cd_shiny_app.R](cd_shiny_app.R) | Shiny 대시보드 |
| [cd_references.md](cd_references.md) | 참고문헌 (약 50편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 항원 경로(장관 글루텐·고유판 DGP·tTG2 활성화·장 투과성), 면역(IEL·CD4 T세포·IFN-γ·IL-17·IL-21·B세포·항-tTG IgA), 조직 병리(융모 높이 VH·융모/음와 비율 CrD), 흡수 기능(흡수 면적·철 저장·BMD), 약물 PK(DrugGut→DrugPlasma)
- **주요 치료 시나리오**: ① 비치료(정상 식이), ② 완전 GFD, ③ 부분 GFD(5% 오염), ④ GFD+라라조타이드 2 mg TID, ⑤ GFD+ZED1227(tTG2 억제제), ⑥ GFD+AMG714(항-IL-15, RCD)
- **보정/근거**: ZED1227 2상 시험(Tye-Din 등), Kelly 융모 위축 회복 데이터, Fernández-Bañares 항체 역가 데이터 참고

## Shiny 대시보드 (Dashboard)
환자 프로파일(HLA 아형·글루텐 섭취량·진단 시 항-tTG IgA), 약동학(신약 후보 혈장 농도), 면역 반응 탭(IEL·CD4T·항체 역가), 조직 병리 탭(융모 높이·V:C 비율), 치료 시나리오 비교(6개), 영양·장기 예후 바이오마커(철·BMD·흡수 면적) 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("cd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("cd_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg cd_qsp_model.dot -o cd_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [cd_references.md](cd_references.md) 참조 (약 50편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
