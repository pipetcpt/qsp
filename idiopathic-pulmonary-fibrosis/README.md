# 특발성 폐섬유화증 (IPF) (Idiopathic Pulmonary Fibrosis, IPF) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 호흡기

[![IPF QSP Model](ipf_qsp_model.png)](ipf_qsp_model.svg)

## 개요 (Overview)

특발성 폐섬유화증(IPF)은 폐포 상피 손상 후 이상 회복 반응으로 섬유아세포·근섬유아세포가 활성화되고 과도한 ECM(세포외기질)이 침착되어 폐 실질이 비가역적으로 파괴되는 진행성 섬유성 간질성 폐질환입니다. 중앙 생존기간이 진단 후 약 3~5년에 불과하며, 연간 발생률은 10만 명당 약 10~20명입니다. 상피-중간엽 전환(EMT)과 TGF-β1이 핵심 구동 인자로, AEC2(제2형 폐포 상피세포) 손상 → TGF-β 신호 → 섬유아세포 활성 → 콜라겐 침착의 악순환이 FVC 및 DLCO의 진행성 감소를 초래합니다. 닌테다닙과 피르페니돈이 현재 진행 속도를 늦출 수 있는 유일한 항섬유화제입니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 폐포 상피 손상 | 노화·흡연·유전(MUC5B·TERT)·환경 → AEC2 손상·ER 스트레스 | 반복 미세 손상 |
| TGF-β1 신호 활성 | AEC2 손상 → TGF-β1 분비 → Smad2/3 → 섬유아세포 활성화·EMT | 섬유아세포 증식 |
| M2 대식세포 분극 | TGF-β·IL-4·IL-13 → M2 대식세포 → 항섬유화 억제·TGF-β 추가 분비 | 섬유화 촉진 피드백 |
| ROS 산화 스트레스 | 미토콘드리아 기능 저하·NADPH산화효소 → ROS 과잉 → AEC2 추가 손상 | 섬유화 가속 |
| 근섬유아세포 분화 | TGF-β → α-SMA+ 근섬유아세포 → 수축·콜라겐 대량 합성 | ECM 침착, 폐 경직 |
| MMP/TIMP 불균형 | TIMP 과잉 → MMP 억제 → ECM 분해 감소 → 콜라겐 축적 | 폐 구조 파괴, 확산 능력 저하 |

## 주요 약물 표적 (Drug Targets)

- **피르페니돈 (Pirfenidone)**: TGF-β 신호 억제·항산화·항염증 다중 기전; ASCEND·CAPACITY 3상에서 FVC 감소 억제 입증
- **닌테다닙 (Nintedanib)**: VEGFR·FGFR·PDGFR 동시 억제(3중 타이로신 키나제 억제) → 섬유아세포 증식·이동 차단; INPULSIS 3상 근거
- **피르페니돈+닌테다닙 병합**: 상호보완 기전; 단독 대비 추가 FVC 보존 효과 탐색 중
- **TGF-β1 직접 차단 (bintrafusp alfa 등)**: 개발 중; Smad 경로 상류 차단
- **LPA1 수용체 길항제 (BMS-986020)**: 섬유아세포 이동 억제; 임상 2상 데이터 존재
- **항섬유화 조합 전략**: IL-13·IL-4Ra 차단(lebrikizumab), galunisertib(TGFβ-RI 억제제) 연구 중

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [ipf_qsp_model.dot](ipf_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 110+ 노드 / 12 클러스터) |
| [ipf_qsp_model.svg](ipf_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [ipf_qsp_model.png](ipf_qsp_model.png) | PNG 이미지 (150 dpi) |
| [ipf_mrgsolve_model.R](ipf_mrgsolve_model.R) | mrgsolve ODE 모델 (약 17 구획 / 5개 치료 시나리오) |
| [ipf_shiny_app.R](ipf_shiny_app.R) | Shiny 대시보드 |
| [ipf_references.md](ipf_references.md) | 참고문헌 (약 50편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK(피르페니돈 장·중심·말초, 닌테다닙 장·중심·말초) + 질환 PD(AEC2, TGF-β1, M2 대식세포, ROS, 활성 섬유아세포, 근섬유아세포, 콜라겐, MMP, TIMP) + 임상 지표(FVC%, DLCO%) 구획 포함
- **주요 치료 시나리오**: ① 위약(자연 경과), ② 피르페니돈 801 mg TID, ③ 닌테다닙 150 mg BID, ④ 피르페니돈+닌테다닙 병합, ⑤ 피르페니돈 저용량 267 mg TID
- **보정/근거**: ASCEND(King 2014, NEJM) — 연간 FVC 감소 −235 mL → −125 mL, INPULSIS(Richeldi 2014, NEJM) — 연간 FVC 감소 −223 mL → −115 mL 데이터를 기반으로 FVC·DLCO 감소 곡선을 정성적으로 보정

## Shiny 대시보드 (Dashboard)

환자 프로파일(기저 FVC%, DLCO%, 흡연력, GAP 점수) 탭, 약물 PK 동역학(피르페니돈·닌테다닙 혈중 농도), TGF-β·섬유화 PD 지표(콜라겐·근섬유아세포), 폐기능 임상 엔드포인트(FVC%·DLCO% 추이), 5개 치료 시나리오 비교(연간 FVC 감소 비교), 바이오마커(TGF-β1·MMP-7·SP-D) 탭으로 구성됩니다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("ipf_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("ipf_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg ipf_qsp_model.dot -o ipf_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [ipf_references.md](ipf_references.md) 참조 (약 50편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
