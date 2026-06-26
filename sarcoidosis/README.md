# 사르코이드증 (Sarcoidosis, SARC) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 호흡기

[![SARC QSP Model](sarc_qsp_model.png)](sarc_qsp_model.svg)

## 개요 (Overview)

사르코이드증은 원인 불명의 비건락성(non-caseating) 육아종이 폐·림프절·피부·눈·심장 등 다장기에 형성되는 전신 염증성 질환으로, 유병률은 10~40/100,000명이며 흑인과 스칸디나비아인에서 높습니다. Th1/IFN-γ 편향 면역반응과 활성화 대식세포가 상호 자극하여 육아종을 형성·유지하는 것이 핵심 발병기전입니다. 약 2/3은 자연 관해되지만, 나머지는 폐섬유화·폐고혈압·심장 사르코이드증으로 진행하며 생명을 위협합니다. 코르티코스테로이드가 1차 치료이며, 메토트렉세이트·아자티오프린, 난치성에서 항TNF 제제가 사용됩니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 항원 제시 및 Th1 편향 | 불명 항원(미코박테리아 유사?) → APC → CD4+ Th1 → IFN-γ·IL-2 과생성 | 육아종 개시, 폐침윤 |
| 대식세포 M1 활성화 | IFN-γ → 대식세포 M1 분극 → TNF·IL-12·IL-18 → 상피세포양 거대세포 | 비건락성 육아종 형성 |
| 육아종 형성·유지 | TNF·IFN-γ → CXCL10·CCL2 → CD4+ T세포·단핵구 지속 동원 | 육아종 확대·영속화 |
| TGF-β/섬유화 경로 | 만성 육아종 → TGF-β1 분비 → 근섬유아세포 활성화 | 진행성 폐섬유화 |
| 칼슘 대사 이상 | 육아종 대식세포의 1α-수산화효소 과활성 → 1,25(OH)₂D₃ 상승 | 고칼슘혈증, 신석증 |
| 심장·전도 사르코이드증 | 심근 육아종·섬유화 → 전도 장애·심실 부정맥 | 심장 돌연사 위험 |
| NF-κB·mTOR 경로 | TNF·IL-1 → NF-κB; 영양 신호 → mTOR → 대식세포 증식 | 만성 육아종 생존 신호 |

## 주요 약물 표적 (Drug Targets)

- **코르티코스테로이드 (프레드니솔론)**: 광범위 항염·Th1 억제; 증상성 폐 침범·신경·심장 사르코이드증 1차 치료
- **메토트렉세이트**: 엽산 길항 → T세포·대식세포 기능 억제; 스테로이드 절감제
- **아자티오프린/마이코페놀레이트**: 퓨린·피리미딘 합성 억제; 만성 유지 치료
- **하이드록시클로로퀸**: TLR9 억제·리소솜 산성화 차단; 피부·고칼슘혈증 사르코이드증
- **항TNF 제제 (인플릭시맙, 아달리무맙)**: TNF 중화 → 육아종 유지 신호 차단; 난치성 폐·신경·심장 사르코이드증
- **테트라사이클린계/클로파자민**: 소규모 증거; 신경 사르코이드 일부에서 보조

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [sarc_qsp_model.dot](sarc_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 100+ 노드 / 12 클러스터) |
| [sarc_qsp_model.svg](sarc_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [sarc_qsp_model.png](sarc_qsp_model.png) | PNG 이미지 (150 dpi) |
| [sarc_mrgsolve_model.R](sarc_mrgsolve_model.R) | mrgsolve ODE 모델 (약 21 구획 / 약 17개 시나리오) |
| [sarc_shiny_app.R](sarc_shiny_app.R) | Shiny 대시보드 |
| [sarc_references.md](sarc_references.md) | 참고문헌 (약 50편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 스테로이드/메토트렉세이트/인플릭시맙 PK 구획 + Th1 활성화·IFN-γ 동태, 육아종 크기 지수, TGF-β·섬유화 모듈, 혈청 ACE·칼슘 예측 구획, 폐기능(FVC·DLCO) 모사 모듈 포함
- **주요 치료 시나리오**: 무치료 자연 관해/진행 이분 경과, 프레드니솔론 단독, 메토트렉세이트, 아자티오프린, 프레드니솔론+메토트렉세이트, 인플릭시맙, 인플릭시맙+메토트렉세이트, 하이드록시클로로퀸 등
- **보정/근거**: Baughman 등 *Am J Respir Crit Care Med* 스테로이드 RCT, GRAPPA 인플릭시맙 연구, 유럽 사르코이드증 가이드라인 임상 데이터 기반

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) **환자 프로파일** — 병기(Scadding)·침범 장기·혈청 ACE·칼슘 설정; (2) **PK 프로파일** — 스테로이드/면역억제제/생물학제제 혈중 농도; (3) **PD 주요지표** — IFN-γ·TNF 억제율, 육아종 크기 동태; (4) **임상 엔드포인트** — FVC·DLCO·6MWD·혈청 ACE 경시 변화; (5) **시나리오 비교** — 치료 전략별 관해율·섬유화 진행 비교; (6) **바이오마커** — 혈청 ACE, 칼슘, sIL-2R, 흉부 X선 병기.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("sarc_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("sarc_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg sarc_qsp_model.dot -o sarc_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [sarc_references.md](sarc_references.md) 참조 (약 50편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
