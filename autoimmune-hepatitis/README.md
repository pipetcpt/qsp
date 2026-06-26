# 자가면역 간염 (Autoimmune Hepatitis, AIH) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 소화기·간담도

[![AIH QSP Model](aih_qsp_model.png)](aih_qsp_model.svg)

## 개요 (Overview)
자가면역 간염(AIH)은 자가반응 T세포가 간세포를 공격하고 자가항체(ANA, ASMA, LKM-1)와 고감마글로불린혈증을 동반하는 만성 자가면역 간질환입니다. 전 세계 유병률은 인구 10만 명당 약 10~17명으로 여성에서 더 흔하며, 치료하지 않으면 간경변·간부전으로 진행합니다. 핵심 발병기전은 자가반응 Th1/Th17 세포에 의한 간세포 직접 손상과 조절 T세포(Treg)의 기능 저하입니다. 1차 치료는 프레드니솔론 단독 또는 아자티오프린 병합이며, 난치성 환자에는 MMF 또는 리툭시맙이 사용됩니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 자가반응 T세포 활성화 | HLA-DR3/DR4, 분자 모방, CTLA-4 다형성 | 간세포 직접 세포독성 |
| Th1 주도 염증 | IFN-γ, TNF-α 과생성 | 간세포 괴사, ALT 급등 |
| Treg 기능 부전 | FoxP3 발현 감소, TGF-β 결핍 | 자가내성 소실, 만성 활성화 |
| B세포·자가항체 | ANA/ASMA/LKM-1, IgG 증가 | 간세포 추가 손상, 진단 마커 |
| 섬유화 진행 | TGF-β, stellate cell 활성화 | 간경변 위험 증가 |
| GR 매개 면역 억제 | 글루코코르티코이드 수용체, NF-κB 억제 | 스테로이드 유도 관해 |

## 주요 약물 표적 (Drug Targets)
- **프레드니솔론**: 글루코코르티코이드 수용체 → NF-κB·AP-1 억제, Treg 증진 — 1차 유도 치료
- **아자티오프린 → 6-TGN**: HGPRT 경로 → DNA 합성 억제, 림프구 증식 감소 — 유지 치료 표준
- **부데소니드**: 1차 통과 대사 고도 → 간 선택 작용, 전신 부작용 감소 — 비간경변 AIH
- **MMF/MPA**: IMPDH 억제 → T·B세포 증식 차단 — 아자티오프린 불내성 2차 치료
- **리툭시맙 (항-CD20)**: B세포 고갈 → 자가항체 감소 — 난치성·재발성 AIH

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [aih_qsp_model.dot](aih_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 212 노드 / 12 클러스터) |
| [aih_qsp_model.svg](aih_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [aih_qsp_model.png](aih_qsp_model.png) | PNG 이미지 (150 dpi) |
| [aih_mrgsolve_model.R](aih_mrgsolve_model.R) | mrgsolve ODE 모델 (약 22 구획 / 6개 치료 시나리오) |
| [aih_shiny_app.R](aih_shiny_app.R) | Shiny 대시보드 |
| [aih_references.md](aih_references.md) | 참고문헌 (약 45편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 약물 PK 구획(프레드니솔론 2구획, 아자티오프린→6-MP→6-TGN, MMF→MPA, 리툭시맙 2구획 + TMDD CD20), 면역 구획(GR 점유율, Th1, Treg, B세포, 자가항체), 사이토카인(IFN-γ, TGF-β, IL-6), 간 손상 지수, 혈청 ALT
- **주요 치료 시나리오**: ① 비치료(자연경과), ② 프레드니솔론 단독, ③ 프레드니솔론+아자티오프린(표준), ④ 부데소니드+아자티오프린, ⑤ 프레드니솔론+MMF, ⑥ 리툭시맙+프레드니솔론(난치성)
- **보정/근거**: IAIHG(국제 자가면역 간염 그룹) 진단 기준 및 주요 임상시험(Czaja, Manns, Zachou 등) ALT·IgG 시계열 데이터 참고

## Shiny 대시보드 (Dashboard)
환자 프로파일(아형·중증도·체중), 약동학 탭(각 약물 혈장 농도), 간 기능 지표(ALT·IgG·빌리루빈 변화), 면역 바이오마커(Th1/Treg 비율·자가항체 역가), 치료 시나리오 비교(6개 오버레이), 임상 엔드포인트(완전 생화학 관해율·재발 예측) 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("aih_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("aih_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg aih_qsp_model.dot -o aih_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [aih_references.md](aih_references.md) 참조 (약 45편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
