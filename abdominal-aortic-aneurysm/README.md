# 복부 대동맥류 (Abdominal Aortic Aneurysm, AAA) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 심혈관

[![AAA QSP Model](aaa_qsp_model.png)](aaa_qsp_model.svg)

## 개요 (Overview)

복부 대동맥류(AAA)는 복부 대동맥의 최대 직경이 3 cm 이상으로 비정상적으로 확장된 상태로, 70세 이상 남성에서 유병률이 약 4–8%에 달한다. 핵심 발병기전은 대동맥 중막의 평활근세포(VSMC) 소실, 엘라스틴·콜라겐 등 세포외기질(ECM)의 MMP(기질금속단백분해효소) 매개 분해, 그리고 대식세포·T세포 주도의 만성 벽내 염증이다. 대동맥 직경이 5.5 cm를 초과하거나 급속 확장(>1 cm/년)이 있으면 파열 위험이 현저히 증가하며, 파열 시 사망률은 80% 이상이다. 독시사이클린(MMP 억제), 스타틴(항염증·ROS 저감), 베타차단제(혈역학적 부하 완화)가 주요 약물 개입 표적이다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| ECM 분해 | MMP-9, MMP-2 활성화 → 엘라스틴·콜라겐 분해 | 대동맥 벽 강도 저하, 직경 확장 |
| 만성 염증 | 대식세포 활성, TNF-α, ROS 생성 | VSMC 세포사멸, 벽내 염증 지속 |
| VSMC 소실 | 세포사멸 가속 · 증식 저하 | 중막 약화 |
| 혈역학적 스트레스 | 혈류 난류, 복강내압 전달 | 벽 응력 증가, ILT 형성 |
| 산화 스트레스 | ROS → NF-κB 활성 | MMP 전사 촉진, VSMC 손상 |
| 혈전내강 형성 | 혈소판 활성, 피브린 축적 | ILT 성장, 벽 영양 결핍 |
| 스타틴 효과 | NF-κB 억제, ROS 저감, MMP 하향 | 진행 억제 |

## 주요 약물 표적 (Drug Targets)

- **독시사이클린 (Doxycycline)**: MMP-2/9 직접 억제 → ECM 보전, 동물모델에서 AAA 진행 억제 확인
- **스타틴 (Statins)**: 플레이오트로픽 항염증 효과, NF-κB 억제 및 ROS 감소, VSMC 보호
- **베타차단제 (Beta-blockers, 예: Propranolol)**: 수축기 혈압 및 ΔP/Δt 감소 → 벽 응력 저감
- **ACE 억제제/ARB**: 안지오텐신 II 매개 MMP 유도 억제, 혈압 관리
- **항혈소판제/항응고제**: ILT 확장 억제

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [aaa_qsp_model.dot](aaa_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 137 노드 / 13 클러스터) |
| [aaa_qsp_model.svg](aaa_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [aaa_qsp_model.png](aaa_qsp_model.png) | PNG 이미지 (150 dpi) |
| [aaa_mrgsolve_model.R](aaa_mrgsolve_model.R) | mrgsolve ODE 모델 (약 25 구획 / 다수 치료 시나리오) |
| [aaa_shiny_app.R](aaa_shiny_app.R) | Shiny 대시보드 |
| [aaa_references.md](aaa_references.md) | 참고문헌 (약 48편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 독시사이클린·스타틴·프로프라놀롤 각각 2–3구획 PK + 대식세포·TNF-α·ROS·MMP-9·MMP-2·엘라스틴·콜라겐·VSMC·ILT·대동맥직경 10개 PD 구획
- **주요 치료 시나리오**: ① 무치료 자연경과, ② 독시사이클린 단독, ③ 스타틴 단독, ④ 독시사이클린 + 스타틴 병용, ⑤ 베타차단제 + 스타틴, ⑥ 삼중 병용 요법
- **보정/근거**: ADAM 연구(베타차단제), 스타틴 관찰 코호트, Longo et al. MMP-9 동물모델 파라미터 참조

## Shiny 대시보드 (Dashboard)

환자 프로파일(연령·초기 직경·위험인자) 설정 탭, 약물 PK 농도-시간 곡선 탭, 대동맥 직경 변화 탭, MMP/ECM 동태 탭, 치료 시나리오 비교 탭, 파열 위험도 및 바이오마커 탭으로 구성되어 있으며, 슬라이더로 용량·투여 간격을 실시간 조정할 수 있다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("aaa_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("aaa_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg aaa_qsp_model.dot -o aaa_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [aaa_references.md](aaa_references.md) 참조 (약 48편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
