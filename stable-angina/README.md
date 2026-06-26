# 안정형 협심증 (Stable Angina, SA) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 심혈관

[![SA QSP Model](sa_qsp_model.png)](sa_qsp_model.svg)

## 개요 (Overview)

안정형 협심증(만성 관상동맥 증후군, CCS)은 죽상경화성 관상동맥 협착으로 인해 심근의 산소 수요가 공급을 초과할 때 발생하는 예측 가능한 흉통으로, 전 세계적으로 약 1억 1천만 명이 이환되어 있습니다. 핵심 발병기전은 지질 침착·염증·플라크 형성에 의한 관상동맥 협착이며, 운동·정서적 스트레스 시 산소 수급 불균형으로 증상이 재현됩니다. 항허혈 치료(베타차단제·칼슘길항제·질산염·이바브라딘·라놀라진), 항혈소판(아스피린·클로피도그렐), 스타틴·ACEI/ARB가 표준 치료 근간을 이룹니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 죽상경화 플라크 형성 | LDL 산화 → 거품세포 → 지질 코어·섬유성 피막 형성 | 관상동맥 협착, 혈류 예비능 감소 |
| 심근 산소 수급 불균형 | 심박수·수축력·후부하 증가 → MVO₂ 상승 vs. 협착 제한 공급 | 허혈성 흉통, ST 변화 |
| 내피 기능 장애 | eNOS 감소 → NO 결핍 → 혈관 수축·혈소판 활성화 | 관상동맥 긴장도 증가 |
| 교감신경 활성화 | 카테콜라민 → β1 수용체 → 심박수·수축력·산소 소비 증가 | 협심 증상 유발 역치 저하 |
| 후기 나트륨 전류 (INaL) | 허혈 → INaL 상승 → 세포내 Ca²⁺ 과부하 → 이완기 긴장 | 심근 산소 소비 증가, 부정맥 |
| 혈소판 활성화·응집 | ADP·TXA2 → GP IIb/IIIa → 혈전 형성 위험 | 급성 관상동맥 증후군 전환 위험 |
| 염증·지질 경로 | hsCRP·IL-6·ox-LDL → 플라크 불안정화 장기 기전 | 심혈관 사건 위험 증가 |

## 주요 약물 표적 (Drug Targets)

- **베타차단제 (메토프롤롤, 비소프롤롤, 아테놀롤)**: β1 차단 → 심박수·수축력 감소 → MVO₂ 감소; 증상·예후 개선
- **칼슘길항제 (암로디핀, 딜티아젬, 베라파밀)**: L형 Ca²⁺ 채널 차단 → 혈관 확장·심박수 감소
- **질산염 (이소소르비드 모노니트레이트·다이니트레이트, 설하 NTG)**: eNOS 독립적 NO 공급 → 정맥·관상동맥 확장
- **이바브라딘 (Ivabradine)**: If 채널 차단 → 순수 심박수 감소; 동율동 유지 환자 β차단제 대안
- **라놀라진 (Ranolazine)**: 후기 INaL 억제 → Ca²⁺ 과부하 감소; 병용 항협심증
- **아스피린/클로피도그렐**: COX-1/P2Y12 억제 → 혈소판 응집 억제; 심혈관 사건 이차 예방
- **스타틴 (아토르바스타틴, 로수바스타틴)**: HMG-CoA 억제 → LDL 감소·pleiotropic 항염; 플라크 안정화

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [sa_qsp_model.dot](sa_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 100+ 노드 / 15 클러스터) |
| [sa_qsp_model.svg](sa_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [sa_qsp_model.png](sa_qsp_model.png) | PNG 이미지 (150 dpi) |
| [sa_mrgsolve_model.R](sa_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 약 11개 시나리오) |
| [sa_shiny_app.R](sa_shiny_app.R) | Shiny 대시보드 |
| [sa_references.md](sa_references.md) | 참고문헌 (약 49편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 베타차단제/CCB/질산염/이바브라딘/라놀라진의 1~2구획 PK + 심박수·혈압·MVO₂ 동태 모듈, INaL 억제 효과 구획, 협심 발작 빈도 예측, LDL·플라크 면적 장기 경과 모듈 포함
- **주요 치료 시나리오**: 무치료, 베타차단제 단독, CCB 단독, 질산염, 베타차단제+CCB, 이바브라딘, 라놀라진, 아스피린+스타틴, 최적 내과 치료(OMT) 병용 등
- **보정/근거**: COURAGE(OMT vs. PCI), BEAUTIFUL(이바브라딘), MERLIN-TIMI 36(라놀라진), TNT(아토르바스타틴 집중 치료) 임상 데이터 기반

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) **환자 프로파일** — 관상동맥 협착도·기저 심박수·혈압·위험 인자 설정; (2) **PK 프로파일** — 항협심증 약물 혈중 농도 경시 변화; (3) **PD 주요지표** — 심박수·혈압·MVO₂ 감소 동태; (4) **임상 엔드포인트** — 주간 협심 발작 횟수·운동 부하 지속 시간 변화; (5) **시나리오 비교** — 치료 전략별 증상 조절·LDL 감소 비교; (6) **바이오마커** — LDL-C, hsCRP, 혈당, NT-proBNP, 허혈 역치.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("sa_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("sa_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg sa_qsp_model.dot -o sa_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [sa_references.md](sa_references.md) 참조 (약 49편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
