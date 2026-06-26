# 전신 홍반 루푸스 (SLE) (Systemic Lupus Erythematosus, SLE) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 자가면역·류마티스

[![SLE QSP Model](sle_qsp.png)](sle_qsp.svg)

## 개요 (Overview)

전신 홍반 루푸스(SLE)는 I형 인터페론(IFN-α/β) 과활성화와 항이중나선DNA(항dsDNA) 항체·면역복합체에 의한 다장기 손상을 특징으로 하는 만성 자가면역 질환으로, 전 세계 유병률은 약 20~150/100,000명이며 가임기 여성에서 10배 이상 많습니다. 신장(루푸스 신염), 피부, 신경계, 혈액계, 심혈관계 등 전신을 침범하며 재발-관해를 반복합니다. 항말라리아제(하이드록시클로로퀸)가 모든 SLE의 기반 치료이며, 벨리무맙(항BAFF)·아니프롤루맙(항IFNAR)이 표준 면역억제 치료 위에 추가되는 표적 치료제로 승인되어 있습니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| I형 IFN 과활성화 | pDC의 TLR7/9 → IFN-α/β → ISG 상향조절 → 면역 활성화 증폭 | 전신 염증 지속, 장기 손상 가속 |
| 항dsDNA·면역복합체 형성 | 세포 사멸 크로마틴 → 자가항체(항dsDNA·항Sm) → IC 형성 | 신장·피부·혈관 IC 침착 |
| 보체 활성화 | IC → C1q·C3a·C5a → 사구체 손상·호중구 동원 | 루푸스 신염, 혈관염 |
| B세포 과활성화 (BAFF 경로) | BAFF 과발현 → B세포 생존·형질세포 분화 → 자가항체 대량 생성 | 항dsDNA 역가 상승, 보체 소비 |
| Tfh–형질세포 축 | IL-21·IL-6 → Tfh → GC 반응 → 고친화성 자가항체 성숙 | 자가면역 지속·악화 |
| NETosis 경로 | 호중구 저등급 과립구(LDG) → NETosis → dsDNA 방출 → pDC 활성화 | IFN-α 분비 양성 피드백 |
| T조절세포(Treg) 결핍 | Treg 수·기능 저하 → 자가반응 T세포 억제 실패 | 면역 관용 붕괴, 재발 |

## 주요 약물 표적 (Drug Targets)

- **하이드록시클로로퀸 (HCQ)**: TLR7/9·IFN-α 억제; 모든 SLE 기반 치료, 플레어 예방·장기 예후 개선
- **코르티코스테로이드**: 광범위 항염·면역억제; 급성 플레어·장기 침범 단기 조절
- **마이코페놀레이트/아자티오프린/사이클로포스파마이드**: 세포 증식 억제; 루푸스 신염 유도·유지
- **벨리무맙 (Belimumab)**: 항BAFF → B세포 생존 억제; 활성 SLE(혈청 양성)에서 승인
- **아니프롤루맙 (Anifrolumab)**: 항IFNAR1 → I형 IFN 신호 차단; 중등도~중증 활성 SLE 승인
- **보코로소스(Voclosporin)/타크롤리무스**: 칼시뉴린 억제 → Th 세포 활성화 억제; 루푸스 신염

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [sle_qsp.dot](sle_qsp.dot) | Graphviz 기계론적 지도 소스 (약 100+ 노드 / 9 클러스터) |
| [sle_qsp.svg](sle_qsp.svg) | SVG 벡터 이미지 (확대 가능) |
| [sle_qsp.png](sle_qsp.png) | PNG 이미지 (150 dpi) |
| [sle_model.R](sle_model.R) | mrgsolve ODE 모델 (약 26 구획 / 약 38개 시나리오) |
| [shiny_app/](shiny_app/) | Shiny 대시보드 (`shiny_app/app.R`) |
| [sle_references.md](sle_references.md) | 참고문헌 (약 66편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: HCQ/벨리무맙/아니프롤루맙/MMF/사이클로포스파마이드 PK 구획 + I형 IFN 생성 동태, pDC 활성화 구획, BAFF·B세포·형질세포 모듈, 항dsDNA 역가 구획, 보체 C3/C4 동태, 루푸스 신염 GFR 모듈, SLEDAI 점수 추정 구획 포함
- **주요 치료 시나리오**: 무치료 자연 경과, HCQ 단독, 프레드니솔론, MMF, 아자티오프린, 사이클로포스파마이드, 벨리무맙, 아니프롤루맙, HCQ+벨리무맙, HCQ+아니프롤루맙, MMF+벨리무맙+HCQ 삼중 병용 등 다수
- **보정/근거**: BLISS-52/-76(벨리무맙), TULIP-1/-2(아니프롤루맙), AURA-LV·AURORA(보코로소스), ALMS(MMF vs. IV-CYC 루푸스 신염) 임상 데이터 기반

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) **환자 프로파일** — 항dsDNA·보체 기저치·침범 장기·SLEDAI 설정; (2) **PK 프로파일** — 생물학제제/소분자 혈중 농도 시간 경과; (3) **PD 주요지표** — IFN 점수·BAFF·항dsDNA 억제 동태; (4) **임상 엔드포인트** — SLEDAI·BILAG·GFR·단백뇨 경시 변화; (5) **시나리오 비교** — 치료 전략별 관해율·장기 손상 비교; (6) **바이오마커** — 항dsDNA, 보체 C3/C4, IFN 신호 점수, 혈구수, 단백뇨.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("sle_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("shiny_app/")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg sle_qsp.dot -o sle_qsp.svg
```

## 참고문헌 (References)

자세한 인용은 [sle_references.md](sle_references.md) 참조 (약 66편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
