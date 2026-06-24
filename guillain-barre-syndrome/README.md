# 길랭-바레 증후군 (Guillain-Barré Syndrome, GBS) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 신경

[![GBS QSP Model](gbs_qsp_model.png)](gbs_qsp_model.svg)

## 개요 (Overview)

길랭-바레 증후군은 선행 감염(주로 캄필로박터 제주니, CMV, EBV)에 대한 면역반응이 분자모방(molecular mimicry)을 통해 말초신경 항강글리오시드 항체(항-GM1, 항-GD1a, 항-GQ1b)를 유발하고, 보체 매개 축삭 또는 수초 손상으로 급성 이완성 마비가 발생하는 자가면역 신경병증입니다. 연간 발생률은 10만 명당 약 1~2명으로 선행 감염 후 2~4주 내에 발병합니다. 두 가지 주요 아형인 AIDP(급성 염증성 탈수초 다발신경병증, 서구)과 AMAN(급성 운동 축삭 신경병증, 아시아)은 발병기전과 예후가 다릅니다. IVIG와 혈장교환이 유일하게 효과가 입증된 치료이며, 심각한 경우 기계환기가 필요합니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 분자모방·항강글리오시드 항체 | 캄필로박터 LOS-GM1 교차반응 → B세포·형질세포 → IgG/IgM 항강글리오시드 Ab | 항체 역가 상승 |
| 보체 활성·MAC 형성 | Ab-강글리오시드 결합 → 고전 경로 C3b → C5a + MAC(C5b-9) | 수초 또는 축삭막 손상 |
| AIDP 수초 손상 | Ab+보체 → Schwann 세포막 공격 → 수초 박리 → 전도 차단 | 사지 마비·감각 소실 |
| AMAN 축삭 손상 | 항-GM1/GD1a + MAC → Ranvier 결절 손상 → 축삭 변성 | 운동 마비·회복 지연 |
| T세포 염증 반응 | Th1·Th17 → IFN-γ·IL-17 → 대식세포 신경 침윤 | 탈수초·축삭 파괴 촉진 |
| 호흡근 침범 | 상행성 마비 → 횡격막·호흡보조근 약화 | FVC 감소, 기계환기 필요 |

## 주요 약물 표적 (Drug Targets)

- **정맥 면역글로불린 (IVIG)**: Fc수용체 포화·보체 억제·항이디오타입 작용으로 항체 매개 신경 손상 차단; 표준 치료(0.4 g/kg×5일)
- **혈장교환 (Plasmapheresis, PE)**: 순환 항강글리오시드 항체 직접 제거; IVIG 대안 치료
- **에쿨리주맙 (Eculizumab)**: 항-C5 단클론항체 → MAC 형성 차단; AMAN에서 실험적 적용, 전임상 효과 확인
- **지지 치료**: 기계환기(FVC < 20 mL/kg), 심박 모니터링, 재활 치료

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [gbs_qsp_model.dot](gbs_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 110+ 노드 / 12 클러스터) |
| [gbs_qsp_model.svg](gbs_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [gbs_qsp_model.png](gbs_qsp_model.png) | PNG 이미지 (150 dpi) |
| [gbs_mrgsolve_model.R](gbs_mrgsolve_model.R) | mrgsolve ODE 모델 (약 26 구획 / 6개 치료 시나리오) |
| [gbs_shiny_app.R](gbs_shiny_app.R) | Shiny 대시보드 |
| [gbs_references.md](gbs_references.md) | 참고문헌 (약 55편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 병원체·DC·대식세포 활성, B세포·형질세포, 항-GM1/GD1a/GQ1b 항체, T세포(Th1·Th17·Treg), 사이토카인(IL-6·TNF-α·IL-10), 보체(C3b·MAC·C5a), 수초 손상·축삭 손상·신경기능, GBS 점수·Hughes 등급·FVC%, IVIG PK(중심·말초), 혈장교환 효과, 에쿨리주맙 PK 구획 포함
- **주요 치료 시나리오**: ① 무치료 AIDP 기저선, ② IVIG 조기 투여(7일), ③ IVIG 지연 투여(14일), ④ 혈장교환, ⑤ IVIG+에쿨리주맙, ⑥ AMAN 무치료
- **보정/근거**: Dutch GBS Study Group, IGOS 코호트, Halstead et al. (eculizumab 전임상, Brain 2008), Dutch 혈장교환 시험 데이터를 기반으로 Hughes 등급 회복 곡선을 정성적으로 보정

## Shiny 대시보드 (Dashboard)

환자 프로파일(선행 감염, 아형 AIDP/AMAN, 기저 Hughes 등급, 항체 아형) 탭, IVIG·혈장교환 PK 시뮬레이션, 항체·보체·T세포 PD 지표, 신경기능 및 GBS 임상 엔드포인트(Hughes 등급, GBS 점수, FVC%), 6개 치료 시나리오 비교, 바이오마커(NfL, 항강글리오시드 역가) 탭으로 구성됩니다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("gbs_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("gbs_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg gbs_qsp_model.dot -o gbs_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [gbs_references.md](gbs_references.md) 참조 (약 55편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
