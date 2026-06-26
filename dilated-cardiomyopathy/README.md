# 확장성 심근병증 (Dilated Cardiomyopathy, DCM) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 심혈관

[![DCM QSP Model](dcm_qsp_model.png)](dcm_qsp_model.svg)

## 개요 (Overview)

확장성 심근병증(Dilated Cardiomyopathy, DCM)은 심실 확장과 수축기능 저하(EF < 40%)를 특징으로 하는 심근 질환으로, 심부전의 가장 흔한 원인 중 하나입니다. 유병률은 인구 10만 명당 약 36명이며, 유전성(티틴 변이 등), 바이러스성, 독성, 자가면역성 등 다양한 원인이 존재합니다. 핵심 병태생리는 심근 손상 후 신경호르몬(RAAS·교감신경) 활성화가 심실 역리모델링을 유발하는 악순환으로, GDMT(ARNI·베타차단제·MRA·SGLT2i)가 신경호르몬 억제를 통해 역리모델링을 유도합니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| RAAS 활성화 | 안지오텐신 II → AT1R → 심근 비대·섬유화 | LV 확장, LVEDP 상승 |
| 교감신경 활성화 | 노르에피네프린 → β1AR → 심박수·수축력 증가 | 심근 산소소비 증가, 세포사 |
| 알도스테론 과잉 | 나트륨 저류·심근 섬유화 촉진 | 심실 경직, 부정맥 위험 |
| 나트리우레틱 펩타이드 | ANP/BNP → 혈관확장·나트륨 배설 | BNP 상승(심부전 진행 지표) |
| 심근 섬유화 | TGF-β → 콜라겐 합성 | 심실 경직, 박출 기능 저하 |
| 미토콘드리아 기능장애 | 에너지 대사 이상, ROS 증가 | 심근세포 자멸사 |
| 전신 염증 | TNF-α, IL-6 → 심근 억압 효과 | 심부전 진행 가속 |

## 주요 약물 표적 (Drug Targets)

- **ARNI** (사쿠비트릴/발사르탄): 네프릴리신 억제(BNP 증가) + ARB(AT1R 차단) → 이중 신경호르몬 억제
- **베타차단제** (카르베딜올, 메토프롤롤): β1AR 차단 → 심박수 감소, 심근 보호
- **MRA** (스피로노락톤, 에플레레논): 알도스테론 수용체 차단 → 섬유화 억제
- **SGLT2억제제** (다파글리플로진, 엠파글리플로진): 오스몰 삼투압·에너지 대사 개선
- **ACE억제제** (에날라프릴): 안지오텐신 II 생성 억제, 브라디키닌 축적

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [dcm_qsp_model.dot](dcm_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 148 노드 / 11 클러스터) |
| [dcm_qsp_model.svg](dcm_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [dcm_qsp_model.png](dcm_qsp_model.png) | PNG 이미지 (150 dpi) |
| [dcm_mrgsolve_model.R](dcm_mrgsolve_model.R) | mrgsolve ODE 모델 (약 24 구획 / 5 치료 시나리오) |
| [dcm_shiny_app.R](dcm_shiny_app.R) | Shiny 대시보드 |
| [dcm_references.md](dcm_references.md) | 참고문헌 (약 52편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK 구획(에날라프릴·카르베딜올·스피로노락톤·사쿠비트릴·다파글리플로진 각 depot+중심) + PD 구획(안지오텐신 II, 알도스테론, BNP, LVEF, LV 용적, 심박수, 섬유화 지수)
- **주요 치료 시나리오**: ① 무치료(자연 경과), ② ACE억제제 단독, ③ ACE억제제 + 베타차단제, ④ ARNI + 베타차단제 + MRA(3제), ⑤ 완전 GDMT 4제(ARNI + 베타차단제 + MRA + SGLT2i)
- **보정/근거**: PARADIGM-HF(LCZ696), DAPA-HF, EMPEROR-Reduced 임상시험 데이터 및 Konstam et al. 역리모델링 문헌 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① 환자 프로파일(기저 LVEF·BNP·NYHA 등급 설정), ② PK 탭(5개 약물 혈중 농도 추이), ③ 심장 PD 탭(LVEF·LV 용적·BNP 변화), ④ 임상 엔드포인트(6분 보행 거리·NYHA 개선·입원율 추정), ⑤ 시나리오 비교(5개 치료 전략 동시 비교), ⑥ 바이오마커(안지오텐신 II·알도스테론·섬유화 지수 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("dcm_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("dcm_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg dcm_qsp_model.dot -o dcm_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [dcm_references.md](dcm_references.md) 참조 (약 52편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
