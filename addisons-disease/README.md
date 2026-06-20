# 애디슨병 (원발성 부신부전) (Addison's Disease, ADD) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![ADD QSP Model](add_qsp_model.png)](add_qsp_model.svg)

## 개요 (Overview)

애디슨병(원발성 부신부전)은 부신피질의 자가면역 파괴로 코르티솔과 알도스테론 분비가 만성적으로 결핍되는 질환으로, 서구에서 인구 10만 명당 약 100–140명의 유병률을 보인다. 자가면역 기전에서는 T세포 매개 부신피질 세포 파괴와 21-히드록실라제(CYP21A2)에 대한 자가항체가 핵심이다. 코르티솔 결핍은 HPA 축의 음성 피드백을 소실시켜 ACTH·CRH가 현저히 상승하고, 알도스테론 결핍은 저나트륨혈증·고칼륨혈증·저혈압을 유발한다. 부신위기는 생명을 위협하는 응급으로, 스트레스 상황에서 스테로이드 용량을 증량해야 한다. 히드로코르티손과 플루드로코르티손 보충 요법의 PK/PD가 모델의 핵심이다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| HPA 축 조절 이상 | CRH → ACTH 상승 (음성 피드백 소실) | 피부 색소 침착(ACTH-MSH 과다) |
| 코르티솔 결핍 | GR 점유 감소 → 면역 억제 해제 | 피로, 저혈압, 저혈당 |
| 알도스테론 결핍 | MR 점유 감소 → 신장 Na⁺ 재흡수 감소 | 저나트륨, 고칼륨, 탈수 |
| 자가면역 파괴 | CD4⁺/CD8⁺ T세포, 21-OHase 항체 | 부신피질 섬유화·소실 |
| 스트레스 반응 결핍 | 코르티솔 급증 불가 → 부신위기 | 쇼크, 사망 위험 |
| 코르티솔 일주기 | HPA축 일주기 리듬 소실 | 수면 장애, 대사 이상 |
| DHEA 결핍 | 부신 안드로겐 부족 | 피로감, 성기능 저하(특히 여성) |

## 주요 약물 표적 (Drug Targets)

- **히드로코르티손 (Hydrocortisone)**: 코르티솔 보충; 1일 15–25 mg 분할 투여, 일주기 리듬 모방이 목표
- **플루드로코르티손 (Fludrocortisone)**: 미네랄코르티코이드 보충; 혈압·전해질 조절
- **DHEA (Dehydroepiandrosterone)**: 여성 환자에서 삶의 질 개선 목적 보충
- **스트레스 용량 (Sick day rules)**: 발열·수술 시 용량 2–3배 증량 프로토콜
- **히드로코르티손 연속 주입 펌프**: 생리적 일주기 패턴 재현 연구 중

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [add_qsp_model.dot](add_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 248 노드 / 10 클러스터) |
| [add_qsp_model.svg](add_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [add_qsp_model.png](add_qsp_model.png) | PNG 이미지 (150 dpi) |
| [add_mrgsolve_model.R](add_mrgsolve_model.R) | mrgsolve ODE 모델 (약 23 구획 / 다수 치료 시나리오) |
| [add_shiny_app.R](add_shiny_app.R) | Shiny 대시보드 |
| [add_references.md](add_references.md) | 참고문헌 (약 49편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 히드로코르티손(경구 2구획 PK) + 플루드로코르티손(2구획) + DHEA 구획 + CRH·ACTH·코르티솔(내인성)·GR-free·GR-bound·GR-mRNA·나트륨·칼륨·혈압·혈당·BMD·ACTH-MSH 신호·부신 잔존 조직·부신위기 위험 PD 구획
- **주요 치료 시나리오**: ① 무치료(진단 전), ② 표준 HC 3회 분복, ③ 아침 편중 HC, ④ 일주기 펌프 주입, ⑤ 스트레스 용량 대응, ⑥ HC + Fludrocortisone + DHEA 삼중 보충
- **보정/근거**: Bornstein et al. 내분비학회 가이드라인, Johannsson et al. 일주기 HC 임상시험 파라미터 참조

## Shiny 대시보드 (Dashboard)

환자 프로파일(체중·성별·잔존 부신 기능) 탭, HC/Fludrocortisone PK 및 코르티솔 일주기 프로파일 탭, 전해질·혈압 동태 탭, HPA 축 바이오마커(ACTH·CRH) 탭, 치료 시나리오 비교(표준 vs. 일주기 vs. 스트레스 용량) 탭, 장기 결과(BMD·대사) 탭으로 구성된다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("add_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("add_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg add_qsp_model.dot -o add_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [add_references.md](add_references.md) 참조 (약 49편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
