# 만성 B형 간염 (Chronic Hepatitis B, CHB) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도

[![CHB QSP Model](chb_qsp_model.png)](chb_qsp_model.svg)

## 개요 (Overview)

만성 B형 간염(CHB)은 전 세계 약 2억 9,600만 명이 감염된 공중보건 문제로, 간경변 및 간세포암(HCC)의 주요 원인이다. HBV는 감염 간세포 핵 내에 cccDNA(공유결합 폐환형 DNA)를 형성하여 항바이러스제 치료 후에도 지속되며, 면역 관용기·면역 활성기·비활성기 등 자연경과가 복잡하다. 뉴클레오시드/뉴클레오타이드 유사체(엔테카비르·테노포비르)는 바이러스 복제를 억제하여 간섬유화 진행을 늦추고, 페그인터페론-α2a는 면역 조절과 직접 항바이러스 작용을 병행한다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| cccDNA 지속·전사 | 핵내 cccDNA → pgRNA → HBV DNA 복제 사이클 | 바이러스 지속, 면역 회피 |
| 면역 관용·T세포 소진 | CTL 소진, 면역 관문(PD-1/Tim-3) 발현 증가 | 만성화, 면역 활성기 전환 지연 |
| 선천 면역 억제 | HBV의 IFN 신호 차단(cGAS-STING) | 초기 바이러스 제거 실패 |
| 간세포 손상 | CTL 매개 세포용해, ALT 상승 | 간염 급성 악화 위험 |
| 간성상세포(HSC) 활성화 | TGF-β, TNF-α → HSC → 콜라겐 합성 | 간섬유화(F0→F4) 진행 |
| HCC 발생 위험 | 고바이러스혈증, 간경변, HBsAg 지속 | 연간 HCC 위험 0.5~3% |
| HBsAg 생산·분비 | cccDNA 전사 + 통합 HBV DNA | 혈청 HBsAg 양성 지속 |

## 주요 약물 표적 (Drug Targets)

- **엔테카비르(ETV)**: dNTP 경쟁 → HBV 역전사효소 억제, 내성 장벽 높음
- **테노포비르 디소프록실(TDF) / TAF**: 뉴클레오타이드 유사체 → HBV DNA 중합효소 억제
- **페그인터페론-α2a(Peg-IFN)**: 선천·적응 면역 활성화 + 직접 항바이러스, HBsAg 소실 가능성
- **siRNA/ASO (연구 중)**: HBsAg mRNA 분해 → HBsAg 혈청 농도 감소 → 면역 회복 유도
- **캡시드 조립 억제제(CAM)**: 새로운 기전, cccDNA 보충 차단

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [chb_qsp_model.dot](chb_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 438 노드 / 10 클러스터) |
| [chb_qsp_model.svg](chb_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [chb_qsp_model.png](chb_qsp_model.png) | PNG 이미지 (150 dpi) |
| [chb_mrgsolve_model.R](chb_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 6개 치료 시나리오) |
| [chb_shiny_app.R](chb_shiny_app.R) | Shiny 대시보드 |
| [chb_references.md](chb_references.md) | 참고문헌 (약 54편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: ETV/TDF PK (gut·plasma·인트라셀룰러 삼인산) + Peg-IFN SC·plasma + 표적세포(T), 감염세포(I), HBV DNA(V), cccDNA, HBsAg, CTL 반응, T세포 소진, 선천IFN, ALT, HSC 활성, 간섬유화, HCC 누적 위험 PD 구획
- **주요 치료 시나리오**: (1) 무치료, (2) 엔테카비르 0.5 mg QD, (3) TDF 300 mg QD, (4) Peg-IFN-α2a 48주, (5) ETV + Peg-IFN 병용, (6) ETV + siRNA 병용
- **보정/근거**: 엔테카비르 및 테노포비르 바이러스 억제 동태는 Perelson AS 등의 바이러스동태 모델(Hepatology 2012), HCC 위험 보정은 PAGE-B 점수 데이터 기반

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) 환자 프로파일 — 면역 단계·기저 바이러스량·간섬유화 단계 설정; (2) PK 탭 — ETV/TDF/Peg-IFN 농도 시계열; (3) PD 주요지표 — HBV DNA(log10), cccDNA, HBsAg, CTL; (4) 임상 엔드포인트 — ALT 정상화, 간섬유화, HCC 위험 누적; (5) 시나리오 비교 — 6가지 치료 전략의 5년 결과; (6) 바이오마커 — HBsAg 정량, HBeAg, 혈소판(간경변 대리지표)

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("chb_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("chb_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg chb_qsp_model.dot -o chb_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [chb_references.md](chb_references.md) 참조 (약 54편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
