# 다발성 경화증 (Multiple Sclerosis, MS) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 신경

[![MS QSP Model](ms_qsp.png)](ms_qsp.svg)

## 개요 (Overview)

다발성 경화증(MS)은 중추신경계(CNS)의 자가면역성 만성 염증·탈수초·축삭 손상 질환으로, 전 세계 약 280만 명이 이환되어 있으며 여성에서 약 3배 더 흔하다. 자가반응 Th1/Th17 세포 및 B세포가 혈액뇌관문(BBB)을 통과하여 중추신경계 백질로 침투하고, 미엘린을 분해하여 재발-완화 또는 진행성 신경장애를 유발한다. 핵심 발병기전으로는 Th1(IFN-γ)/Th17(IL-17) 축에 의한 신경염증, B세포 및 형질세포의 국소 항체 생성, 미세아교세포 활성화, 희소돌기아교세포(OPC) 손상, 축삭 뉴로필라멘트 방출이 포함된다. 치료 표적으로는 림프구 재순환(S1P 수용체), CD20+ B세포(항CD20), 나탈리주맙(VLA-4/α4-인테그린), 인터페론-β(면역조절) 등이 주요하다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| Th1/Th17 분화 | HLA-DRB1*15:01, IL-12, IL-23, IFN-γ, IL-17 | 백질 병변, 탈수초 |
| BBB 투과성 증가 | α4-인테그린/VCAM-1 결합, MMP 활성화 | 림프구 CNS 침윤 |
| B세포·항체 기전 | CD20+ B세포, 형질세포, IgG 과두 밴드 | 뇌척수액 올리고클론 밴드 |
| 미세아교세포 활성화 | NF-κB, TNF-α, NO, ROS | 수초·축삭 손상 진행 |
| 탈수초·재수초화 실패 | 희소돌기아교세포 전구체(OPC) 분화 장애 | T-score 저하, 장애 축적 |
| 축삭 손상 | NfL 혈청 상승, GFAP 방출 | 비가역적 장애(EDSS 증가) |
| 신경보호 결핍 | BDNF, IGF-1 감소 | 진행성 뇌 위축 |

## 주요 약물 표적 (Drug Targets)

- **인터페론-β (IFN-β1a/1b)**: Th1 억제, BBB 안정화; 재발률 약 30% 감소
- **나탈리주맙 (Natalizumab)**: α4-인테그린 차단 → BBB 통과 차단; 재발률 약 68% 감소 (AFFIRM)
- **오크렐리주맙 (Ocrelizumab)**: 항CD20 B세포 고갈; RRMS/PPMS 재발·진행 억제 (OPERA I/II, ORATORIO)
- **시포니모드/피네고리모드 (S1P 수용체 조절제)**: 림프구 이차 림프기관 격리; 말초 림프구 60–80% 감소
- **디메틸 푸마레이트 (DMF)**: Nrf2 경로 활성화, Th1→Th2 면역 전환
- **클라드리빈 (Cladribine)**: 선택적 림프구 감소; 장기 면역 재구성

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [ms_qsp.dot](ms_qsp.dot) | Graphviz 기계론적 지도 소스 (약 195 노드 / 12 클러스터) |
| [ms_qsp.svg](ms_qsp.svg) | SVG 벡터 이미지 (확대 가능) |
| [ms_qsp.png](ms_qsp.png) | PNG 이미지 (150 dpi) |
| [ms_mrgsolve_model.R](ms_mrgsolve_model.R) | mrgsolve ODE 모델 (약 27 구획 / 7개 치료 시나리오) |
| [ms_shiny_app.R](ms_shiny_app.R) | Shiny 대시보드 |
| [ms_references.md](ms_references.md) | 참고문헌 (약 40편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 인터페론-β 2구획 PK, 나탈리주맙 2구획+수용체 결합, 오크렐리주맙 2구획, S1P 조절제 2구획, DMF·클라드리빈 단일구획; 말초 Th1·Th17·Treg·B세포, BBB 투과성, CNS 내 cTh1·cTh17·미세아교세포·희소돌기아교세포·OPC·미엘린·축삭, NfL·GFAP 바이오마커
- **주요 치료 시나리오**: ① 무치료(placebo) ② IFN-β-1a 30μg IM QW ③ 나탈리주맙 300mg IV Q4W ④ 오크렐리주맙 600mg IV Q6M ⑤ 시포니모드 경구 QD ⑥ DMF 240mg BID ⑦ 클라드리빈 경구(누적 3.5mg/kg)
- **보정/근거**: AFFIRM(나탈리주맙), OPERA I/II(오크렐리주맙), TRANSFORMS(피네고리모드), DEFINE(DMF), CLARITY(클라드리빈) 임상시험 데이터 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① **환자 프로파일** (성별·나이·체중·기저 EDSS 설정), ② **PK** (혈장 약물 농도-시간 곡선, 수용체 점유율), ③ **PD 주요지표** (말초 Th1/Th17/Treg/B세포 변화), ④ **임상 엔드포인트** (ARR, EDSS, T2 병변 수), ⑤ **시나리오 비교** (7개 약물 또는 병용요법 직접 비교), ⑥ **바이오마커** (혈청 NfL, GFAP 추이).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("ms_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("ms_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg ms_qsp.dot -o ms_qsp.svg
```

## 참고문헌 (References)

자세한 인용은 [ms_references.md](ms_references.md) 참조 (약 40편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
