# 다낭성 난소 증후군 (PCOS) (Polycystic Ovary Syndrome, PCOS) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 내분비·대사

[![PCOS QSP Model](pcos_qsp_model.png)](pcos_qsp_model.svg)

## 개요 (Overview)
다낭성 난소 증후군(PCOS)은 가임기 여성의 5~15%에서 유병하는 가장 흔한 내분비 질환으로, 고안드로겐혈증·만성 무배란·다낭성 난소 소견을 특징으로 한다(Rotterdam 2003 기준). 핵심 발병기전은 ① 시상하부 GnRH 박동 빈도 증가(키스펩틴-뉴로키닌-다이노르핀[KNDy] 뉴런 이상) → LH 과다·FSH 상대적 감소, ② 난소 난포막세포의 CYP17A1 과활성 → 안드로겐 과생성, ③ 인슐린저항성 → 고인슐린혈증 → SHBG 감소 → 유리 안드로겐 증가의 상호 증폭 악순환이다. 주요 치료 표적은 HPO축 조절, 인슐린감수성 개선, 안드로겐 과잉 억제이다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| GnRH 박동 경로 | KNDy 뉴런 이상, GnRH 박동 빈도 증가 | LH/FSH 비율 증가(>2) |
| 난소 스테로이드 경로 | CYP17A1 과활성, LH 과자극 | 안드로스텐디온·테스토스테론 과잉 |
| 인슐린저항성 경로 | INSR 후 신호 이상(PI3K-Akt), IRS-1 인산화 | 고인슐린혈증·보상성 과인슐린 |
| SHBG-유리 안드로겐 경로 | 간 SHBG 합성 감소(인슐린 억제) | 생체이용 테스토스테론 증가 |
| 난포 성숙 장애 경로 | AMH 과다, FSH 부족, 난포 정체 | 다낭성 난소·무배란·불임 |
| 대사 합병증 경로 | 내장지방 축적, 이상지질혈증, 염증 | 제2형 당뇨·심혈관 위험 증가 |

## 주요 약물 표적 (Drug Targets)
- **메트포르민**: AMPK 활성화 → 인슐린감수성 개선, 간 포도당 신생 억제, 난소 안드로겐 생성 간접 억제
- **클로미펜/레트로졸**: 에스트로겐 수용체 길항 / 아로마타제 억제 → FSH 증가 → 배란 유도
- **경구피임약 (OCP)**: 에티닐에스트라디올+프로게스틴 → LH 억제·SHBG 증가·안드로겐 감소
- **스피로노락톤/시프로테론**: 안드로겐 수용체 차단 → 다모증·여드름 개선
- **GnRH 유사체**: HPO축 하향 조절 — 체외수정 시 과자극 방지

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [pcos_qsp_model.dot](pcos_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 196 노드 / 14 클러스터) |
| [pcos_qsp_model.svg](pcos_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [pcos_qsp_model.png](pcos_qsp_model.png) | PNG 이미지 (150 dpi) |
| [pcos_mrgsolve_model.R](pcos_mrgsolve_model.R) | mrgsolve ODE 모델 (약 24 구획 / 6 치료 시나리오) |
| [pcos_shiny_app.R](pcos_shiny_app.R) | Shiny 대시보드 |
| [pcos_references.md](pcos_references.md) | 참고문헌 (약 43편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: GnRH 구동·LH·FSH·총테스토스테론·에스트라디올·프로게스테론·AMH·AFC상태·우성난포·인슐린·혈당·SHBG·유리테스토스테론·IGFBP1·CRP·BMI·다모증 점수 PD 구획, 메트포르민·레트로졸·에티닐에스트라디올·스피로노락톤·클로미펜 PK 구획
- **주요 치료 시나리오**: ① S1 무치료 PCOS, ② S2 메트포르민, ③ S3 레트로졸, ④ S4 경구피임약, ⑤ S5 메트포르민+레트로졸, ⑥ S6 스피로노락톤
- **보정/근거**: PCOS 표준 임상시험(NICHD PPCOSII) 및 메트포르민·레트로졸 무작위 비교 연구에서 배란율·호르몬 변화 데이터를 참고 보정

## Shiny 대시보드 (Dashboard)
환자 프로파일(Rotterdam 표현형, BMI, 인슐린저항성 지수) · HPO축 호르몬 PK/PD · 난포 성숙·배란 PD · 대사 지표(HOMA-IR, 혈당, 지질) 임상 엔드포인트 · 치료 시나리오 비교(배란율·안드로겐 변화) · 바이오마커(AMH, SHBG, 유리 테스토스테론) 탭으로 구성

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("pcos_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("pcos_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg pcos_qsp_model.dot -o pcos_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [pcos_references.md](pcos_references.md) 참조 (약 43편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
