# 피부근염 (Dermatomyositis, DM) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 자가면역·류마티스

[![DM QSP Model](dm_qsp_model.png)](dm_qsp_model.svg)

## 개요 (Overview)

피부근염(Dermatomyositis, DM)은 I형 인터페론(IFN) 경로와 보체계가 주도하는 근육 및 피부의 염증성 자가면역 질환으로, 성인에서 연간 발생률은 인구 100만 명당 약 5~10명입니다. 근위부 근육 약화와 특징적 피부 발진(헬리오트로프 발진, Gottron 구진)이 동반되며, 근육 미세혈관의 보체 매개 손상이 핵심 병태생리입니다. 주요 치료 표적은 IFN-α/β 신호, B세포 매개 자가항체 생성, CD8+ T세포 활성화 및 보체 경로이며, 스테로이드·IVIG·리툭시맙·JAK억제제가 대표 치료제입니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| I형 IFN 신호 | pDC 유래 IFN-α/β → MxA, IFIT 과발현 | IFN 점수 상승, 근육 손상 지속 |
| 보체 매개 혈관병증 | C3b/MAC 형성 → 근육 미세혈관 내피 손상 | 모세혈관 소실, 근섬유 허혈 |
| B세포·자가항체 | 항Jo-1, 항Mi-2, 항MDA5 자가항체 생성 | 근육·폐 섬유화 유발 |
| CD8+ T세포 침윤 | 근육 내 세포독성 T세포 → 근섬유 직접 손상 | CK 상승, 근력 저하 |
| Th17/Treg 불균형 | IL-17 과발현, Treg 감소 | 만성 염증 지속 |
| JAK-STAT 경로 | JAK1/2 → STAT1 활성화 | IFN 반응 유전자 발현 증폭 |
| TGF-β 섬유화 | 만성 염증 → TGF-β → 근육·폐 섬유화 | FVC 저하, 기능 손실 |

## 주요 약물 표적 (Drug Targets)

- **코르티코스테로이드** (프레드니솔론): 광범위 면역억제, NF-κB 경로 억제
- **IVIG** (정맥 면역글로불린): 보체 억제, Fc수용체 봉쇄, 자가항체 중화
- **메토트렉세이트(MTX)**: 엽산 길항, T세포·B세포 증식 억제
- **리툭시맙(RTX)**: 항CD20 → B세포 고갈, 자가항체 생성 억제
- **JAK억제제** (룩소리티닙, 바리시티닙): JAK1/2 억제 → IFN 신호 차단
- **하이드록시클로로퀸(HCQ)**: toll-like receptor 신호 억제, 피부 병변 개선

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [dm_qsp_model.dot](dm_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 175 노드 / 13 클러스터) |
| [dm_qsp_model.svg](dm_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [dm_qsp_model.png](dm_qsp_model.png) | PNG 이미지 (150 dpi) |
| [dm_mrgsolve_model.R](dm_mrgsolve_model.R) | mrgsolve ODE 모델 (약 24 구획 / 5 치료 시나리오) |
| [dm_shiny_app.R](dm_shiny_app.R) | Shiny 대시보드 |
| [dm_references.md](dm_references.md) | 참고문헌 (약 57편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK 구획(프레드니솔론·IVIG·MTX·리툭시맙·JAK억제제 각 1~3구획) + PD 구획(IFN 점수, 보체, B세포, 자가항체, 근육 손상, CK, MMT8, CDASI, FVC, Treg, Th17, CD8 활성화, 모세혈관)
- **주요 치료 시나리오**: ① 무치료(자연 경과), ② 프레드니솔론 단독, ③ 프레드니솔론 + IVIG, ④ 리툭시맙 + 스테로이드, ⑤ JAK억제제 + 스테로이드
- **보정/근거**: MERI 및 IMACS 코어 세트 임상 데이터, Oddis et al. ArthritisRheum 리툭시맙 RCT, 바리시티닙 2상 시험(MyoPath) 등 참고

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① 환자 프로파일(인구통계·자가항체·기저 상태), ② PK 탭(각 약물 혈중 농도 시뮬레이션), ③ 면역 PD 탭(IFN 점수·보체·B세포·자가항체 추이), ④ 임상 엔드포인트(MMT8 근력·CK·CDASI 피부·FVC), ⑤ 시나리오 비교(5개 치료군 동시 비교), ⑥ 바이오마커 패널(Treg/Th17 비율·CD8 활성화 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("dm_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("dm_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg dm_qsp_model.dot -o dm_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [dm_references.md](dm_references.md) 참조 (약 57편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
