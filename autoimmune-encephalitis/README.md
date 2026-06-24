# 자가면역 뇌염 (Autoimmune Encephalitis, AIE) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 신경

[![AIE QSP Model](aie_qsp_model.png)](aie_qsp_model.svg)

## 개요 (Overview)

자가면역 뇌염(AIE)은 신경 세포표면 항원(주로 NMDA 수용체, LGI1, CASPR2, GABA-B 수용체 등)에 대한 항체가 뇌 기능을 직접 손상시키는 자가면역 질환군이다. 항NMDAR 뇌염이 가장 흔한 형태로, 젊은 여성에서 주로 발생하며 난소 기형종과 동반될 수 있다. 발병기전은 B세포가 생산한 자가항체가 혈액뇌장벽(BBB)을 통과하여 해마·대뇌피질의 NMDA 수용체 GluN1 아단위에 결합하고, 수용체를 내재화·감소시켜 신경전달 이상을 유발한다. 임상적으로 신경정신과적 증상(정신증·인지장애)→발작→운동이상→의식 저하의 전형적 진행을 보인다. 1차 면역치료(스테로이드·IVIG·혈장교환)에 불충분한 반응 시 2차 치료(리툭시맙·사이클로포스파마이드)를 시행한다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 항체 생산 | 배중심 B세포 → 형질세포 → anti-NMDAR IgG | CSF 항체 상승 |
| BBB 투과 | 염증 → 타이트 정션 손상 → IgG CSF 유입 | CNS 항체 농도 증가 |
| NMDAR 내재화 | 항체-GluN1 결합 → 클라트린 매개 내재화 | 시냅스 NMDAR 밀도↓ |
| 미세아교세포 활성 | NMDAR 손실 → MG 활성화 → IL-6·TNF-α 분비 | 신경염증 지속 |
| 도파민 균형 이상 | NMDAR↓ → DA 탈억제 → 과다 도파민 | 정신증 증상 |
| 흥분독성 | 과도한 글루타메이트 → Ca²⁺ 유입 | 인지기능 손상 |
| 발작 역치 | NMDAR↓ → 피질 흥분/억제 균형 이상 | 간질발작 |

## 주요 약물 표적 (Drug Targets)

- **메틸프레드니솔론 (Methylprednisolone)**: 1차 면역치료; BBB 안정화, 염증 억제; 고용량 펄스 후 경구 테이퍼링
- **IVIG (정맥 면역글로불린)**: 항체 중화, Fc 수용체 포화; 3–5일 주입; 스테로이드와 병용 가능
- **혈장교환술 (Plasmapheresis)**: aPL 항체 신속 제거; 치명적 증례에서 우선 고려
- **리툭시맙 (Rituximab, 항-CD20)**: B세포 고갈 → 항체 생산 억제; 재발 예방 및 2차 치료 표준
- **토실리주맙 (Tocilizumab, 항-IL-6R)**: 난치성·재발성 AIE에서 IL-6 매개 BBB 손상 억제 근거 증가
- **종양 제거**: 기형종 동반 시 종양 절제가 면역 자극원 제거로 핵심 치료

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [aie_qsp_model.dot](aie_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 183 노드 / 12 클러스터) |
| [aie_qsp_model.svg](aie_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [aie_qsp_model.png](aie_qsp_model.png) | PNG 이미지 (150 dpi) |
| [aie_mrgsolve_model.R](aie_mrgsolve_model.R) | mrgsolve ODE 모델 (약 25 구획 / 다수 치료 시나리오) |
| [aie_shiny_app.R](aie_shiny_app.R) | Shiny 대시보드 |
| [aie_references.md](aie_references.md) | 참고문헌 (약 65편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: B세포(배중심·형질세포·장수명 형질세포·기억세포) + 혈청/CSF 항-NMDAR 항체 + BBB 투과도 + 미세아교세포·IL-6·GFAP + NMDAR 발현·글루타메이트·인지 기능·정신증·발작 PD 구획 + IVIG(2구획) + 메틸프레드니솔론(2구획) + 리툭시맙(2구획) + 토실리주맙(2구획) PK 구획
- **주요 치료 시나리오**: ① 무치료(자연 경과), ② 스테로이드 단독, ③ IVIG 단독, ④ 스테로이드 + IVIG 병용(1차 표준), ⑤ 혈장교환 + 스테로이드, ⑥ 2차 치료(리툭시맙), ⑦ 토실리주맙(난치성)
- **보정/근거**: Titulaer et al. Lancet Neurol 2013 예후 코호트, Dalmau et al. 초기 기술 논문, 리툭시맙 케이스 시리즈 파라미터 참조

## Shiny 대시보드 (Dashboard)

환자 프로파일(항체 역가·종양 동반 여부·초기 중증도) 탭, 면역치료 PK 및 B세포 동태 탭, NMDAR 발현·항체 CSF 농도 탭, 임상 엔드포인트(mRS·인지·발작) 탭, 치료 전략 비교(1차 vs. 2차 면역치료) 탭, 신경염증 바이오마커(GFAP·IL-6·항체 역가) 탭으로 구성된다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("aie_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("aie_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg aie_qsp_model.dot -o aie_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [aie_references.md](aie_references.md) 참조 (약 65편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
