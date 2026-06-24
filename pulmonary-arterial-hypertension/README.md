# 폐동맥 고혈압 (PAH) (Pulmonary Arterial Hypertension, PAH) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 심혈관

[![PAH QSP Model](pah_qsp_model.png)](pah_qsp_model.svg)

## 개요 (Overview)

폐동맥 고혈압(PAH)은 평균 폐동맥압(mPAP) ≥20 mmHg, 폐혈관저항(PVR) ≥3 Wood 단위를 특징으로 하는 진행성 폐혈관 질환입니다. 전 세계 유병률은 약 15~60명/100만 명이며, 진단 후 치료 없이는 중앙 생존 기간이 2.8년에 불과합니다. 핵심 발병기전은 엔도텔린-1 과발현, 프로스타사이클린 결핍, 산화질소(NO) 생성 저하에 의한 폐혈관 수축·리모델링으로, 폐동맥 평활근세포(PASMC) 증식과 내피기능 장애가 RV 후부하를 점진적으로 증가시킵니다. 엔도텔린 수용체 길항제(ERA), PDE5 억제제, 프로스타사이클린 유사체, sGC 자극제가 현재 승인 치료 표적입니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 엔도텔린-1 경로 | ET-1 과발현 → ETA/ETB 수용체 → PASMC 수축·증식 | PVR 상승, 폐혈관 리모델링 |
| NO-sGC-cGMP 경로 | eNOS 감소 → NO 결핍 → sGC 활성 저하 → 혈관 수축 | 내피 기능 장애, 혈소판 활성화 |
| 프로스타사이클린 경로 | PGI2 합성효소 감소 → cAMP 감소 → 혈관 수축·혈소판 응집 | 폐혈관 긴장도 증가 |
| BMPR2/TGF-β 신호 | BMPR2 기능 소실 돌연변이 → TGF-β 우세 → PASMC 과증식 억제 실패 | 혈관 내막·중막 비후 |
| 염증 및 면역 경로 | 대식세포 침윤·IL-6·CCL2 → 혈관주위 염증 | 병변 진행 가속, 예후 악화 |
| 우심실 적응 및 부전 | 후부하 증가 → RV 비대 → RV-PA 결합 이상 → 우심부전 | 6MWD 감소, 사망 |
| 저산소증 유도 경로 | HIF-1α → VEGF·ET-1 상향 조절, Warburg 대사 | PASMC 생존 촉진, 관강 폐색 |

## 주요 약물 표적 (Drug Targets)

- **엔도텔린 수용체 길항제 ERA (보센탄, 암브리센탄, 마시텐탄)**: ETA(±ETB) 차단 → 혈관 수축·증식 억제; SERAPHIN 임상 근거
- **PDE5 억제제 (실데나필, 타달라필)**: cGMP 분해 억제 → NO 경로 강화; SUPER-1·PHIRST 임상 근거
- **sGC 자극제 (리오시구아트)**: NO 비의존적 sGC 직접 활성화; PATENT 임상 근거
- **프로스타사이클린 유사체 (에포프로스테놀, 트레프로스티닐, 일로프로스트)**: IP 수용체 → cAMP 상승; 정맥·피하·흡입 경로
- **IP 수용체 작용제 (셀렉시팍)**: 선택적 IP 수용체 활성화; GRIPHON 임상 근거
- **병용 치료 (ERA+PDE5i, ERA+프로스타사이클린)**: AMBITION·TRITON 임상 근거

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [pah_qsp_model.dot](pah_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 100+ 노드 / 10 클러스터) |
| [pah_qsp_model.svg](pah_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [pah_qsp_model.png](pah_qsp_model.png) | PNG 이미지 (150 dpi) |
| [pah_mrgsolve_model.R](pah_mrgsolve_model.R) | mrgsolve ODE 모델 (약 17 구획 / 약 19개 시나리오) |
| [pah_shiny_app.R](pah_shiny_app.R) | Shiny 대시보드 |
| [pah_references.md](pah_references.md) | 참고문헌 (약 48편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: ERA/PDE5i/프로스타사이클린/sGC 자극제의 1~2구획 PK + ET-1·NO·PGI2 농도 동태 모듈, PVR·mPAP 출력 구획, RV 기능(EF, CI) 예측 모듈 포함
- **주요 치료 시나리오**: 무치료 자연 경과, 보센탄, 암브리센탄, 실데나필, 타달라필, 리오시구아트, 에포프로스테놀, 셀렉시팍, 암브리센탄+타달라필(AMBITION), ERA+프로스타사이클린 삼중 병용 등
- **보정/근거**: SERAPHIN(마시텐탄), PATENT-1(리오시구아트), GRIPHON(셀렉시팍), AMBITION(암브리센탄+타달라필), TRITON(마시텐탄+타달라필) 임상 데이터 기반 파라미터 보정

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) **환자 프로파일** — WHO 기능 등급·혈역학·동반 질환 설정; (2) **PK 프로파일** — 약물 혈중 농도 시간 경과; (3) **PD 주요지표** — ET-1·NO·PGI2 억제율, PVR 변화; (4) **임상 엔드포인트** — 6분 보행 거리(6MWD), mPAP, NT-proBNP 경시 변화; (5) **시나리오 비교** — 단일/병용 요법 비교; (6) **바이오마커** — 우심실 기능, TAPSE, 심박출량.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("pah_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("pah_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg pah_qsp_model.dot -o pah_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [pah_references.md](pah_references.md) 참조 (약 48편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
