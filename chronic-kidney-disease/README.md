# 만성 신부전 (Chronic Kidney Disease, CKD) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 신장·비뇨

[![CKD QSP Model](ckd_qsp_model.png)](ckd_qsp_model.svg)

## 개요 (Overview)

만성 신부전(CKD)은 전 세계 성인의 약 10~13%에서 발생하며, 당뇨병성 신병증과 고혈압이 가장 흔한 원인이다. 사구체 과여과(초기) → 네프론 손실 → 보상성 비대 → 단백뇨·섬유화의 악순환이 eGFR을 수년에 걸쳐 진행성으로 감소시킨다. RAAS 억제제(ACEi/ARB), SGLT2 억제제(다파글리플로진, DAPA-CKD), 비스테로이드성 MRA(피네레논, FIDELIO-DKD)의 3제 병용이 현재 가장 강력한 신보호 전략으로 확립되어 있다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 사구체 과여과·단백뇨 | 안지오텐신 II → 수출세동맥 수축 → 사구체 내압 상승 | UACR 증가, 발세포 손상 |
| RAAS 과활성 | Ang II → 알도스테론 → TGF-β 분비 | 신섬유화, 혈압 상승, 염증 |
| 산화 스트레스·NF-κB | ROS, TNF-α, IL-6 증가 | 튜불로인터스티셜 염증 |
| 간질 섬유화 | TGF-β → 근섬유아세포 → 콜라겐 과다 침착 | 네프론 대체 → eGFR 감소 |
| CKD-MBD | FGF-23 ↑, 클로토 ↓, VitD 저하 → PTH ↑ | 혈관 석회화, 이차성 부갑상선 기능 항진 |
| 신성 빈혈 | EPO 생산 저하, 헵시딘 ↑ → 헤모글로빈 저하 | 피로, 심혈관 부하 증가 |
| 좌심실 비대(LVH) | 혈압·빈혈·RAAS 과활성 → 심근 비대 | 심혈관 사망 위험 증가 |

## 주요 약물 표적 (Drug Targets)

- **ACE 억제제/ARB**: 라미프릴·로사르탄 — Ang II 억제 → 사구체 내압 감소, 단백뇨 감소, TGF-β 억제
- **SGLT2 억제제(다파글리플로진)**: 근위세뇨관 소디움-포도당 재흡수 차단 → 사구체 내압 감소, 신보호(DAPA-CKD)
- **피네레논(비스테로이드성 MRA)**: 미네랄로코르티코이드 수용체 차단 → 섬유화·염증 억제(FIDELIO/FIGARO-DKD)
- **ESA(적혈구 자극제)/다르베포에틴**: EPO 수용체 자극 → 빈혈 교정
- **인산 결합제/활성 VitD**: 인·PTH 조절 → CKD-MBD 관리

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [ckd_qsp_model.dot](ckd_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 484 노드 / 10 클러스터) |
| [ckd_qsp_model.svg](ckd_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [ckd_qsp_model.png](ckd_qsp_model.png) | PNG 이미지 (150 dpi) |
| [ckd_mrgsolve_model.R](ckd_mrgsolve_model.R) | mrgsolve ODE 모델 (약 30 구획 / 5개 치료 시나리오) |
| [ckd_shiny_app.R](ckd_shiny_app.R) | Shiny 대시보드 |
| [ckd_references.md](ckd_references.md) | 참고문헌 (약 37편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: ACEi/ARB/피네레논/SGLT2i/ESA/인산 결합제 PK 구획 + 네프론 수, eGFR, UACR, Ang II, 알도스테론, 혈압, 대식세포, IL-6, TNF-α, TGF-β, 콜라겐, 인산, 클로토, FGF-23, VitD, PTH, EPO, 헵시딘, 헤모글로빈, 좌심실 비대 지수, 혈관 석회화 PD 구획
- **주요 치료 시나리오**: (1) 자연경과(무치료), (2) ACEi 단독(라미프릴 10 mg), (3) ACEi + 피네레논 20 mg, (4) ACEi + 다파글리플로진 10 mg (DAPA-CKD 방식), (5) 3제 병용(ACEi + Dapa + 피네레논)
- **보정/근거**: eGFR 감소 속도 및 UACR 반응은 DAPA-CKD (NEJM 2020), FIDELIO-DKD (NEJM 2020) 임상시험 데이터 기반

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) 환자 프로파일 — CKD 단계·당뇨 유무·기저 단백뇨 설정; (2) PK 탭 — ACEi/Dapagliflozin/피네레논 혈장 농도; (3) PD 주요지표 — eGFR 추이, UACR, Ang II; (4) 임상 엔드포인트 — 신부전 진행, 헤모글로빈, CKD-MBD 지표; (5) 시나리오 비교 — 5가지 치료 전략 5년 결과; (6) 바이오마커 — FGF-23, 클로토, PTH, 혈관 석회화 점수

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("ckd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("ckd_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg ckd_qsp_model.dot -o ckd_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [ckd_references.md](ckd_references.md) 참조 (약 37편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
