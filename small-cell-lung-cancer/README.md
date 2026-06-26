# 소세포 폐암 (Small Cell Lung Cancer, SCLC) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 종양학

[![SCLC QSP Model](sclc_qsp_model.png)](sclc_qsp_model.svg)

## 개요 (Overview)

소세포 폐암(SCLC)은 전체 폐암의 약 15%를 차지하는 고도로 공격적인 신경내분비 종양으로, 진단 시 약 70%가 광범위기(Extensive Stage, ES)에 해당합니다. 무치료 시 중앙 생존 기간은 2~4개월에 불과하며, 1차 치료(카보플라틴/에토포사이드 + 아테졸리주맙 또는 더발루맙)로도 중앙 전체 생존 기간이 약 12~14개월에 불과합니다. SCLC의 핵심 발병기전은 TP53·RB1의 동시 소실, ASCL1/NEUROD1/YAP1/POU2F3 전사인자에 의한 신경내분비(NE) 아형 결정, DLL3의 비정상적 표면 발현(Notch 억제 배위자), 그리고 고도로 면역 억제적인 종양 미세환경(TME)입니다. DLL3 양성율이 ~80%로 높아 신규 BiTE 항체인 tarlatamab의 표적이 됩니다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| TP53/RB1 소실 | 세포주기 정지 불능 → 급속 증식·분열 | 고도 악성, 빠른 전이 |
| NE 아형 전사 프로그램 | ASCL1/NEUROD1 (NE-high), YAP1/POU2F3 (NE-low) | 치료 반응성 및 내성 결정 |
| DLL3/Notch 경로 | DLL3 과발현 → Notch 억제 → NE 표현형 유지 | Tarlatamab BiTE 표적 |
| DNA 손상 응답 | Topo-II 독소(에토포사이드), 백금계(카보플라틴) → DSB | 1차 치료 근거 |
| 면역 회피 | PD-L1/TIGIT/CD47, MHC-I 소실, MDSC | IO 치료 표적 |
| 세포주기 (CDK4/6) | Rb 비활성화 → E2F → 급속 분열 | 트릴라시클립 골수 보호 근거 |
| 우심실 적응 | 폐 침윤 → 저산소증 → CO 감소 | 임상 악화 지표 |

## 주요 약물 표적 (Drug Targets)

- **백금 + 에토포사이드 (CE)**: DNA 이중나선 절단 → 종양세포 사멸; 1차 표준치료
- **아테졸리주맙 (Atezo)**: PD-L1 차단; IMpower133 OS HR 0.70 (13.9 vs 12.3개월)
- **더발루맙 (Durva)**: PD-L1 차단; CASPIAN OS HR 0.73 (13.0 vs 10.3개월)
- **Tarlatamab**: DLL3 × CD3 이중특이항체(BiTE) → DLL3+ 종양세포로 T세포 유도; DeLLphi-301 ORR 40%, mDOR 9.7개월
- **Lurbinectedin**: RNA 중합효소 II 분해, DNA 알킬화; 2차 ORR 35.2%
- **Topotecan**: TOP1 독소 → SSB; 표준 2차 치료
- **트릴라시클립 (Trilaciclib)**: CDK4/6 억제 → G1 정지 → 조혈모세포 보호 (골수억제 ↓); 화학요법 전 투여

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [sclc_qsp_model.dot](sclc_qsp_model.dot) | Graphviz 기계론적 지도 (100+ 노드, 12 클러스터) |
| [sclc_qsp_model.svg](sclc_qsp_model.svg) | SVG 벡터 이미지 |
| [sclc_qsp_model.png](sclc_qsp_model.png) | PNG 이미지 (150 dpi) |
| [sclc_mrgsolve_model.R](sclc_mrgsolve_model.R) | mrgsolve ODE 모델 (21 구획, 8 시나리오) |
| [sclc_shiny_app.R](sclc_shiny_app.R) | Shiny 대시보드 (6 탭) |
| [sclc_references.md](sclc_references.md) | 참고문헌 (48편+, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조 (21 ODEs)**:
  - PK (13 구획): 아테졸리주맙(2-cmt), 더발루맙(2-cmt), 카보플라틴(1-cmt), 에토포사이드(2-cmt), lurbinectedin(1-cmt), topotecan(1-cmt), tarlatamab(3-cmt), trilaciclib(1-cmt)
  - PD (8 구획): DNA 손상 지표, 민감 종양세포(Ts), G1 정지세포(Ta), 내성세포(Tr), CD8+ T세포, PD-L1 수용체 점유율, NK 세포, IFN-γ
- **주요 치료 시나리오 (8가지)**:
  1. 무치료 자연 경과
  2. CE 화학요법 단독 (카보플라틴 + 에토포사이드)
  3. CE + 아테졸리주맙 (IMpower133)
  4. CE + 더발루맙 (CASPIAN)
  5. Tarlatamab 단독 (2차 치료, DeLLphi-301)
  6. Lurbinectedin 단독 (2차)
  7. Topotecan 단독 (2차)
  8. CE + Trilaciclib (골수보호)
- **보정 근거**: IMpower133(아테졸리주맙), CASPIAN(더발루맙), DeLLphi-301(tarlatamab), ATLANTIS(lurbinectedin), Trigo 2020

## Shiny 대시보드 (Dashboard)

6개 탭:
1. **환자 프로파일** — 병기·ECOG·NE 아형·DLL3 발현·TMB·동반질환 설정
2. **Drug PK** — 면역항암제·화학요법·trilaciclib 혈중농도 경시 변화
3. **분자 경로** — DNA 손상 지표·BiTE 효과·PD-L1 RO·IFN-γ·NK
4. **임상 엔드포인트** — SLD 종양 부담, 반응 분류(PR/SD/PD), 세포 구획별 동태
5. **시나리오 비교** — 8가지 치료 전략 동시 비교 + 효능 요약표
6. **바이오마커 & TME** — CD8T·NK·IFN-γ, NE 마커(ProGRP·NSE), DLL3 발현 수준별 tarlatamab 효능

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("sclc_mrgsolve_model.R")
out <- mrgsim(mod, end = 365 * 24, delta = 24)
plot(out)
# Shiny 대시보드:
shiny::runApp("sclc_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg sclc_qsp_model.dot -o sclc_qsp_model.svg
dot -Tpng -Gdpi=150 sclc_qsp_model.dot -o sclc_qsp_model.png
```

## 참고문헌 (References)

자세한 인용은 [sclc_references.md](sclc_references.md) 참조 (48편+).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
