# 비후성 심근병증 (Hypertrophic Cardiomyopathy, HCM) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 심혈관

[![HCM QSP Model](hcm_qsp_model.png)](hcm_qsp_model.svg)

## 개요 (Overview)

비후성 심근병증(HCM)은 유병률 500명 중 1명으로 가장 흔한 유전성 심근병증이며, MYH7·MYBPC3 등 근절(sarcomere) 단백 유전자 변이에 의한 과수축(hypercontractility)이 핵심 발병기전입니다. 변이 근절은 수축-이완 주기에서 저속 이완(SRX) 상태에서 고속 활성(DRX) 상태로의 전환을 증가시켜 ATP 소비 증가·열 발생·심근 비대·섬유화를 유발합니다. 좌심실유출로 폐쇄(LVOTO)가 동반된 폐쇄형 HCM은 전체의 약 70%를 차지하며, 운동 불내성·실신·심방세동·돌연사의 원인이 됩니다. 마바캄텐(mavacamten)은 최초의 심근 마이오신 억제제로 SRX 상태를 증가시켜 과수축을 억제합니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 근절 과수축 | SRX→DRX 전환 증가 → 과도한 가교 형성 → 과수축 | 좌심실 비대, 이완 불량 |
| 좌심실유출로 폐쇄 | 비대 중격 + 수축기 전방 운동(SAM) → 동적 LVOTO | 압력차 증가, 호흡 곤란 |
| 칼슘 신호 이상 | 세포 내 Ca²⁺ 과부하 → 칼시뉴린-NFAT → 비대 유전자 활성 | 심근 비대 유전자 발현 |
| 심근 섬유화 | TGF-β → 심근 섬유아세포 활성 → 콜라겐 침착 | 간질 섬유화, 이완 장애 |
| ERK/MAPK 신호 | 기계적 스트레스 → ERK1/2 → 세포 성장 경로 활성 | 비대 유지 |
| 심방세동·돌연사 | LVEDP 상승·섬유화 → 심방 확대·심실 빈맥 | 뇌졸중, 급사 위험 |

## 주요 약물 표적 (Drug Targets)

- **마바캄텐 (Mavacamten)**: β-심근 마이오신 ATPase 억제 → SRX 상태 복원 → 과수축 감소·LVOTO 완화; EXPLORER-HCM, VALOR-HCM 시험 근거
- **아피캄텐 (Aficamten)**: 차세대 마이오신 억제제; SEQUOIA-HCM 3상에서 EXPLORER-HCM과 유사한 LVOT 압력차 감소
- **베타차단제 (Metoprolol·Atenolol)**: 심박수·수축력 감소 → LVOTO 및 증상 완화; 1차 약물 치료
- **비라파밀/딜티아젬**: 음성 변시성 작용 → 이완 시간 연장; 베타차단제 불내성 환자
- **디소피라미드**: 음성 수축성 + 항부정맥; LVOTO 지속 시 베타차단제 병합
- **중격 감소 요법**: 외과적 심근절제술(septal myectomy) 또는 알코올 중격 절제술(ASA); 약물 불응성 LVOTO

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [hcm_qsp_model.dot](hcm_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 110+ 노드 / 12 클러스터) |
| [hcm_qsp_model.svg](hcm_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [hcm_qsp_model.png](hcm_qsp_model.png) | PNG 이미지 (150 dpi) |
| [hcm_mrgsolve_model.R](hcm_mrgsolve_model.R) | mrgsolve ODE 모델 (약 19 구획 / 5+개 치료 시나리오) |
| [hcm_shiny_app.R](hcm_shiny_app.R) | Shiny 대시보드 |
| [hcm_references.md](hcm_references.md) | 참고문헌 (약 56편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK(마바캄텐 장·중심·말초, 베타차단제 단일 구획) + 세포 신호(세포 내 Ca²⁺·SR Ca²⁺, 칼시뉴린, NFAT 핵 내 이동, ERK 활성화) + 심장 구조(IVS 두께·LV 질량·TGF-β1·콜라겐·LVOT 압력차·LVEDP·HR·NT-proBNP·트로포닌 I·심방세동 위험) 구획 포함
- **주요 치료 시나리오**: ① 무치료 기저선, ② 마바캄텐 단독, ③ 마바캄텐+베타차단제 병합, ④ 베타차단제 단독, ⑤ 아피캄텐 유사 시나리오, ⑥+ 추가 용량-반응 분석
- **보정/근거**: EXPLORER-HCM(Olivotto 2020, Lancet) — LVOT 압력차 감소 및 KCCQ 개선, VALOR-HCM(Desai 2023, JAMA Cardiol) — 중격 감소 적응증 감소, SEQUOIA-HCM(Nagueh 2024, Lancet) 데이터를 기반으로 정성적 보정

## Shiny 대시보드 (Dashboard)

환자 프로파일(기저 LVOT 압력차, IVS 두께, NYHA 등급, 유전자 변이 선택) 탭, 약물 PK 동역학(마바캄텐 혈중 농도), 근절·세포 신호 PD 지표(SRX 비율·Ca²⁺·칼시뉴린), 심장 구조 임상 엔드포인트(LVOT 압력차·LV 비대·섬유화), NT-proBNP·트로포닌 바이오마커, 5개+ 치료 시나리오 비교, 심방세동·돌연사 위험 지표 탭으로 구성됩니다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("hcm_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("hcm_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg hcm_qsp_model.dot -o hcm_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [hcm_references.md](hcm_references.md) 참조 (약 56편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
