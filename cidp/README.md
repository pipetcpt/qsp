# 만성 염증성 탈수초성 다발신경병증 (Chronic Inflammatory Demyelinating Polyneuropathy, CIDP) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 신경

[![CIDP QSP Model](cidp_qsp_model.png)](cidp_qsp_model.svg)

## 개요 (Overview)

CIDP는 말초신경 말이집(수초)을 자가면역 기전으로 파괴하는 진행성/재발성 신경병증으로, 인구 10만 명당 약 1~8명에서 발생한다. 병원성 자가항체(항-노도파라노딘, 항-컨탁틴 등)와 T세포 매개 염증이 협력하여 슈반세포와 수초를 공격하며, 축삭 손상이 누적되면 비가역적 장애로 이어진다. IVIG(정맥 내 면역글로불린)가 1차 표준 치료이며, 스테로이드·혈장교환·리툭시맙이 대안이고, FcRn 억제제(에프가르티지모드)가 최신 승인 치료제이다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 자가항체 생성·보체 활성화 | 병원성 IgG → 보체 경로 → 막 공격 복합체 | 수초 손상, 신경 전도 차단 |
| T세포(Th1/Th17) 침윤 | IFN-γ, IL-17 → 내피 투과성 ↑ → 신경내막 침윤 | 탈수초, 부종 |
| FcRn 매개 IgG 항상성 | FcRn ↑ → 병원성 IgG 반감기 연장 | 자가항체 지속 |
| 수초 손상·탈락 | 슈반세포 손상 → 수초 분절 탈락 | 신경전도속도(NCV) 감소 |
| 재수초화 억제 | 염증 지속 → 슈반세포 분화 억제 | 비가역적 장애 진행 |
| 축삭 변성 | 만성 탈수초 → 축삭 허혈·에너지 고갈 | 근력 저하, 위축 |
| 신경섬유 손상 마커 | NfL(혈청 신경필라멘트 경쇄) 상승 | 질환 활성도·치료 반응 모니터링 |

## 주요 약물 표적 (Drug Targets)

- **IVIG(정맥 내 면역글로불린)**: FcγR 차단, Fc-FcRn 경쟁 → 자가항체 분해 촉진, Treg 증가, 보체 억제
- **코르티코스테로이드(프레드니솔론)**: T세포·B세포 억제, 혈액신경장벽 투과성 감소
- **혈장교환(PLEX)**: 직접 자가항체·보체 제거
- **리툭시맙(RTX)**: CD20 표적 → B세포 고갈 → 항체 생산 감소 (난치성 CIDP)
- **에프가르티지모드(Efgartigimod)**: FcRn 차단 → 전체 IgG(자가항체 포함) 반감기 단축 (ADHERE 시험)

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [cidp_qsp_model.dot](cidp_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 358 노드 / 12 클러스터) |
| [cidp_qsp_model.svg](cidp_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [cidp_qsp_model.png](cidp_qsp_model.png) | PNG 이미지 (150 dpi) |
| [cidp_mrgsolve_model.R](cidp_mrgsolve_model.R) | mrgsolve ODE 모델 (약 25 구획 / 6개 치료 시나리오) |
| [cidp_shiny_app.R](cidp_shiny_app.R) | Shiny 대시보드 |
| [cidp_references.md](cidp_references.md) | 참고문헌 (약 55편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: IVIG (중심·말초 2구획) + 코르티코스테로이드 PK + 리툭시맙(RTX 혈장·CD20 결합) + 에프가르티지모드 PK + PLEX 효과 구획 + Th1, Th17, Treg, B세포, 형질세포, 대식세포, 보체, 병원성 자가항체, 탈수초화, 재수초화, 축삭 밀도, NfL, NCV, INCAT 장애 점수 PD 구획
- **주요 치료 시나리오**: (1) 무치료, (2) IVIG 2 g/kg q4주, (3) 프레드니솔론 테이퍼, (4) 혈장교환 + IVIG 유지, (5) 리툭시맙 2회 투여, (6) 에프가르티지모드 2사이클
- **보정/근거**: IVIG 반응은 PATH 시험(NEJM 2017), 에프가르티지모드는 ADHERE 시험(Lancet Neurol 2023) 기반

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: (1) 환자 프로파일 — 발병 형태(재발-완화형/진행형)·기저 INCAT 점수 설정; (2) PK 탭 — IVIG/RTX/에프가르티지모드 혈장 농도; (3) PD 주요지표 — 자가항체, 보체, Th1/Th17/Treg; (4) 임상 엔드포인트 — INCAT 장애 점수, NCV, NfL; (5) 시나리오 비교 — 6가지 치료 전략 1년 결과; (6) 바이오마커 — 혈청 NfL, 탈수초화 지수, 재수초화 진행

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("cidp_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("cidp_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg cidp_qsp_model.dot -o cidp_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [cidp_references.md](cidp_references.md) 참조 (약 55편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
