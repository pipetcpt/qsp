# 말초동맥질환 (PAD) (Peripheral Arterial Disease, PAD) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 심혈관

[![PAD QSP Model](pad_qsp_model.png)](pad_qsp_model.svg)

## 개요 (Overview)
말초동맥질환(PAD)은 하지 동맥의 죽상경화로 인한 협착·폐쇄로 발생하며, 전 세계 약 2억 명의 성인에서 유병이 확인된다. 발목상완지수(ABI) 0.9 미만이 진단 기준이며, 간헐성 파행에서 중증 사지 허혈(CLI)까지 다양한 임상 스펙트럼을 보인다. 주요 발병기전은 LDL 산화·내피세포 활성화·단핵구 침윤을 통한 죽상반(plaque) 형성 및 혈소판 과활성화로 인한 혈전성 협착이다. 주요 치료 표적은 혈소판 활성화 경로(P2Y12, COX-1), LDL 수치 조절, 및 혈관 평활근 이완(PDE3 억제)이다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 죽상경화 경로 | ox-LDL, 거품세포, MMP, 내피기능 장애 | 죽상반 형성·ABI 감소 |
| 혈소판 활성화 경로 | P2Y12·TXA2, GPIIb/IIIa | 혈전 형성·급성 혈관 폐쇄 |
| 응고 경로 | 트롬빈·피브린 생성, FXa | 혈전 증폭 |
| 허혈-재관류 경로 | ROS, 근육 허혈, 신경 손상 | 파행·조직 괴사 |
| 측부 혈관 형성 | VEGF, FGF, HIF-1α | 허혈 보상 회복 |
| 전신 염증 경로 | hsCRP, IL-6, 내피 접착 분자 | 심혈관 위험 증폭 |
| 지질 경로 | LDL-C, HDL-C, apoB100 | 죽상반 성장·불안정화 |

## 주요 약물 표적 (Drug Targets)
- **항혈소판제 — 클로피도그렐**: P2Y12 수용체 비가역적 차단 (ADP 경로 억제)
- **항혈소판제 — 아스피린**: COX-1 억제 → TXA2 생성 감소
- **항혈소판제 — 티카그렐러**: P2Y12 가역적 차단, DAPT 구성
- **항응고제 — 리바록사반 저용량**: FXa 억제 (COMPASS 시험 근거)
- **혈관이완제 — 실로스타졸**: PDE3 억제 → cAMP 증가, 파행 거리 개선
- **스타틴 — 아토르바스타틴**: HMG-CoA 환원효소 억제, LDL-C 감소·죽상반 안정화

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [pad_qsp_model.dot](pad_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 193 노드 / 12 클러스터) |
| [pad_qsp_model.svg](pad_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [pad_qsp_model.png](pad_qsp_model.png) | PNG 이미지 (150 dpi) |
| [pad_mrgsolve_model.R](pad_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 7 치료 시나리오) |
| [pad_shiny_app.R](pad_shiny_app.R) | Shiny 대시보드 |
| [pad_references.md](pad_references.md) | 참고문헌 (약 54편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 클로피도그렐/활성대사체·아스피린·티카그렐러·리바록사반·실로스타졸·아토르바스타틴 PK 구획, 혈소판 응집·트롬빈 활성·LDL-C·죽상반 부피·ABI·측부 혈관·보행 거리·내피 기능·hsCRP PD 구획
- **주요 치료 시나리오**: ① 무치료, ② 아스피린 100 mg, ③ 클로피도그렐 75 mg, ④ DAPT(클로피도그렐+아스피린), ⑤ COMPASS(리바록사반 2.5 mg BID + 아스피린), ⑥ 실로스타졸+아스피린, ⑦ 최적 복합요법(클로피도그렐+아스피린+리바록사반+스타틴)
- **보정/근거**: CAPRIE, COMPASS, EUCLID, CASPAR 임상시험 데이터를 기반으로 ABI 변화율 및 MACE 위험 감소 파라미터 보정

## Shiny 대시보드 (Dashboard)
환자 프로파일(Fontaine 분류, 동반 질환, 흡연력) · 약물 PK 프로파일 · 혈소판 응집·응고 억제 PD · ABI·보행거리 임상 엔드포인트 · 치료 시나리오 비교(MACE 위험, 파행 개선) · 혈관 바이오마커(hsCRP, LDL-C, 죽상반 부피) 탭으로 구성

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("pad_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("pad_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg pad_qsp_model.dot -o pad_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [pad_references.md](pad_references.md) 참조 (약 54편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
