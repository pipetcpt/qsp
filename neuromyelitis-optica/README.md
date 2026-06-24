# 시신경척수염 (NMOSD, Neuromyelitis Optica) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 신경

[![NMO QSP Model](nmo_qsp_model.png)](nmo_qsp_model.svg)

## 개요 (Overview)

시신경척수염 스펙트럼장애(NMOSD)는 항AQP4(아쿠아포린-4) IgG 항체가 성상세포의 AQP4를 공격하여 보체-매개 세포독성(CDC) 및 항체-의존 세포독성(ADCC)을 유발하는 희귀 자가면역 CNS 질환이다. 전 세계 유병률은 약 100만 명당 1–10명으로 추정되며 여성(~9:1)과 아시아계에서 더 흔하다. 임상적으로 심한 시신경염(실명 위험)과 척수염(사지 마비·자율신경 장애)이 재발성으로 나타나며, 치료 없이는 누적 장애가 빠르게 악화된다. IL-6 수용체 신호는 B세포 생존 및 항체 생산을 촉진하므로 사트랄리주맙, 에쿨리주맙(C5), 이네빌리주맙(항CD19) 등이 표준 치료로 자리잡았다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 항AQP4 IgG 생산 | BAFF/APRIL 의존 B세포 활성화 → 형질세포 → IgG 분비 | 혈청 AQP4-Ab 역가 상승 |
| 보체 활성화 | C1q → C3 → C5 → MAC(C5b-9) 형성 | 성상세포(Astrocyte) 용해 |
| ADCC | NK세포·호중구 FcγR → 성상세포 파괴 | 병변 확대, 백질 손상 |
| IL-6 신호 | IL-6R/JAK-STAT3 → 형질모세포 생존, Th17 분화 | 재발 빈도 증가 |
| 성상세포 손상 | AQP4 소실 → 물·이온 항상성 파괴 | 조직 부종, 희소돌기아교세포 이차 손상 |
| 탈수초·축삭 손상 | 이차 희소돌기아교세포 손상 | EDSS 장애 축적 |
| CD19+ B세포 축 | CD19 발현 세포(형질모세포 포함) | 항체 생산 지속 |

## 주요 약물 표적 (Drug Targets)

- **에쿨리주맙 (Eculizumab)**: C5 보체 억제 → MAC 형성 차단; NMOSD 재발을 위약 대비 약 94% 감소 (PREVENT 시험)
- **이네빌리주맙 (Inebilizumab)**: 항CD19 → B세포·형질모세포 고갈 (N-MOmentum 시험)
- **사트랄리주맙 (Satralizumab)**: 항IL-6R(재활용형) → IL-6 신호 차단 (SAkuraStar/SAkuraSky 시험)
- **리툭시맙 (Rituximab)**: 항CD20 B세포 고갈; 장기 유지 면역억제 (오프라벨)
- **프레드니솔론**: 급성 재발 시 고용량 IV 스테로이드 펄스; 면역억제 유지
- **미코페놀레이트 모페틸 (MMF)**: 푸린 합성 억제; 리툭시맙과 병용 유지요법

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [nmo_qsp_model.dot](nmo_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 176 노드 / 14 클러스터) |
| [nmo_qsp_model.svg](nmo_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [nmo_qsp_model.png](nmo_qsp_model.png) | PNG 이미지 (150 dpi) |
| [nmo_mrgsolve_model.R](nmo_mrgsolve_model.R) | mrgsolve ODE 모델 (약 29 구획 / 6개 치료 시나리오) |
| [shiny_app/](shiny_app/) | Shiny 대시보드 (`shiny_app/app.R`) |
| [nmo_references.md](nmo_references.md) | 참고문헌 (약 44편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 에쿨리주맙 2구획 IV PK+C5 결합, 이네빌리주맙 2구획 IV PK, 사트랄리주맙 SC 3구획+IL-6R 복합체, 리툭시맙 2구획 IV PK, 프레드니솔론 경구 2구획, MPA 1구획; B세포 성숙 단계(Bnaive→Bact→PB→PC), 항AQP4 IgG, C5·MAC 활성화, 성상세포·희소돌기아교세포 생존율, EDSS, NfL, IL-6·TNF-α
- **주요 치료 시나리오**: ① 무치료 ② 에쿨리주맙 900mg IV Q2W ③ 이네빌리주맙 300mg IV Q6M ④ 사트랄리주맙 120mg SC Q8W ⑤ 리툭시맙+MMF ⑥ 프레드니솔론 펄스(급성기)
- **보정/근거**: PREVENT(에쿨리주맙), N-MOmentum(이네빌리주맙), SAkuraStar(사트랄리주맙) 임상시험 데이터 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① **환자 프로파일** (항체 역가·기저 EDSS·재발력 설정), ② **PK** (혈장 약물 농도-시간 곡선 및 보체/IL-6R 점유율), ③ **PD 주요지표** (B세포·형질세포·항AQP4 IgG 추이), ④ **임상 엔드포인트** (재발률, EDSS, MRI 병변), ⑤ **시나리오 비교** (6개 치료 전략 직접 비교), ⑥ **바이오마커** (NfL, GFAP, IL-6 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("nmo_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("shiny_app/")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg nmo_qsp_model.dot -o nmo_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [nmo_references.md](nmo_references.md) 참조 (약 44편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
