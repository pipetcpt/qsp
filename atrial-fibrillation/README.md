# 심방세동 (Atrial Fibrillation, AF) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 심혈관

[![AF QSP Model](af_qsp_model.png)](af_qsp_model.svg)

## 개요 (Overview)

심방세동(AF)은 가장 흔한 지속성 부정맥으로, 전 세계 성인의 약 1–2% 유병률을 보이며 75세 이상에서는 10%를 초과한다. 발병기전은 전기적 리모델링(유효불응기 단축, 이온채널 발현 변화)과 구조적 리모델링(심방 섬유화, 확장)의 악순환('AF가 AF를 낳는다')이다. 폐정맥 기시부의 이소성 전기활동이 AF를 촉발하고, 심방 내 재진입 회로가 이를 유지시킨다. 전기적 리모델링은 시간에 따라 비가역적으로 진행하여 발작성 → 지속성 → 영구성 AF로 이행한다. 심방 혈전 형성으로 뇌졸중 위험이 5배 증가한다. 율동 조절(cardioversion, 항부정맥제), 심박수 조절(베타차단제, 칼슘채널차단제), 항응고(NOAC)가 세 축의 치료 전략이다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 전기적 리모델링 | ICa,L↓, IK1↑, IKur↓ → ERP 단축 | AF 유지·재진입 회로 형성 |
| 구조적 리모델링 | AngII·TGF-β·ROS → 심방 섬유화 | 전도 지연, AF 영속화 |
| Ca²⁺ 과부하 | SR Ca²⁺ 누출 → 지연후탈분극(DAD) | 이소성 박동, AF 유발 |
| 자율신경계 | 교감·부교감 균형 이상 | 박동수 변동, AF 발작 촉발 |
| NLRP3 염증 | 심방 심근세포 인플라마솜 활성 | 섬유화 촉진, AF 지속 |
| 응고 활성화 | AF 혈류 정체 → 혈소판 활성, 트롬빈 생성 | 심방이 혈전, 뇌졸중 |
| 레닌-안지오텐신 | AngII → SMAD2/3 → 콜라겐 합성 | 심방 섬유화 악순환 |

## 주요 약물 표적 (Drug Targets)

- **아미오다론 (Amiodarone, Vaughan Williams III/I)**: 광범위 이온채널 차단; 율동 유지 효과 우수하나 장기 독성 (갑상선·폐·간)
- **플레카이나이드/프로파페논 (Class IC)**: 발작성 AF 율동 전환; 구조적 심질환 금기
- **메토프롤롤/비소프롤롤 (베타차단제)**: 심박수 조절; 수술후 AF 예방
- **아픽사반/리바록사반/다비가트란 (NOAC)**: 비판막성 AF 뇌졸중 예방; 와파린 대비 두개내 출혈 위험 감소
- **베라파밀/딜티아젬 (칼슘채널차단제)**: 심박수 조절; HFrEF에서 사용 주의

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [af_qsp_model.dot](af_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 184 노드 / 10 클러스터) |
| [af_qsp_model.svg](af_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [af_qsp_model.png](af_qsp_model.png) | PNG 이미지 (150 dpi) |
| [af_mrgsolve_model.R](af_mrgsolve_model.R) | mrgsolve ODE 모델 (약 24 구획 / 다수 치료 시나리오) |
| [af_shiny_app.R](af_shiny_app.R) | Shiny 대시보드 |
| [af_references.md](af_references.md) | 참고문헌 (약 56편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 아미오다론(경구 2구획 PK + 심실 분포) + 아픽사반(2구획) + 메토프롤롤(경구) PK 구획 + 유효불응기(ERP)·AF 부담·심방 섬유화·QTc·심박수·AngII·ROS·SMAD2/3·IL-6·FXa·트롬빈·뇌졸중 위험·노르에피네프린 PD 구획
- **주요 치료 시나리오**: ① 무치료(AF 자연 진행), ② 베타차단제 심박수 조절, ③ 아미오다론 율동 조절, ④ NOAC 항응고 단독, ⑤ 율동 조절 + NOAC, ⑥ 심박수 조절 + NOAC(영구성 AF)
- **보정/근거**: AFFIRM 임상시험(율동 vs. 심박수 조절), ARISTOTLE(아픽사반), ENGAGE AF-TIMI 48 파라미터 참조

## Shiny 대시보드 (Dashboard)

환자 프로파일(CHA₂DS₂-VASc 점수·HAS-BLED·AF 유형) 탭, 항부정맥제 PK 및 이온채널 점유율 탭, AF 부담·ERP 변화 탭, 뇌졸중/출혈 위험 동태 탭, 치료 전략 비교(율동 vs. 심박수 조절) 탭, 바이오마커(섬유화 지표·CRP) 탭으로 구성된다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("af_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("af_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg af_qsp_model.dot -o af_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [af_references.md](af_references.md) 참조 (약 56편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
