# 백반증 (Vitiligo, VIT) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 피부

[![VIT QSP Model](vit_qsp_model.png)](vit_qsp_model.svg)

## 개요 (Overview)

백반증은 CD8+ T세포 매개 멜라닌세포 파괴로 피부에 탈색반이 형성되는 자가면역 피부 질환으로, 전 세계 유병률은 약 0.5~2%이다. 자가반응성 CD8+ T세포가 멜라닌세포 특이 항원(Melan-A, PMEL17 등)을 인식하여 세포독성 파괴를 유발하며, IFN-γ → JAK1/2-STAT1 → CXCL10 분비 → CD8+ T세포 피부 귀소의 자기 강화 피드백 루프가 병변 확장의 핵심이다. JAK1/2 억제제 룩솔리티닙(외용) 크림이 최초로 FDA 승인을 받아 색소 재침착에 효능을 보였으며, 아파멜라노타이드(MC1R 작용제)와 NB-UVB 광선치료가 병용된다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 산화 스트레스·멜라닌세포 취약성 | H₂O₂ 축적, NKG2D 리간드 발현 유도 → NK/NKT 세포 활성화 | 초기 멜라닌세포 손상 및 항원 노출 |
| IFN-γ — JAK-STAT1 축 | CD8+ T세포·NK세포 IFN-γ → JAK1/2 → pSTAT1 → CXCL9/CXCL10 분비 | CXCR3+ CD8+ T세포 피부 귀소 증폭 |
| CD8+ CTL 멜라닌세포 살해 | CXCL10 구배 → CD8+ T세포 병변 귀소 → 퍼포린·그랜자임 B 세포독성 | 멜라닌세포 밀도 감소, 탈색반 확장 |
| Treg 기능 부전 | Foxp3+ Treg 피부 침윤 감소 → IFN-γ/CXCL10 축 억제 실패 | 자가면역 반응 지속·확장 |
| MITF 신호 감소 | IFN-γ → MITF 억제 → 멜라닌 합성 감소, 멜라닌세포 생존 저하 | 멜라닌 생산 중단, 피부 탈색 |
| 모낭 멜라닌세포 저장소 | NB-UVB → 모낭 저장 멜라닌세포 재활성화·이동 | 모낭 주위부터 색소 재침착 시작 |

## 주요 약물 표적 (Drug Targets)

- **JAK1/2 억제제 (외용) — 룩솔리티닙 1.5% 크림**: pSTAT1 억제 → CXCL10 감소 → CD8+ T세포 귀소 차단, 색소 재침착 (TRuE-V1/2 trials)
- **JAK1/2 억제제 (경구) — 룩솔리티닙 10 mg BID**: 전신 IFN-γ/CXCL10 축 억제, 광범위 백반증에 고려
- **MC1R 작용제 — 아파멜라노타이드(Afamelanotide)**: α-MSH 유사체, 모낭 멜라닌세포 증식·이동 촉진, NB-UVB와 병용 시 색소 재침착 가속
- **NB-UVB 광선치료**: 면역 조절(Treg 증가) + 모낭 저장 멜라닌세포 동원
- **국소 칼시뉴린 억제제 — 타크로리무스**: T세포 활성화 차단, 안면·굴측부 백반증 2차 치료

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [vit_qsp_model.dot](vit_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 262 노드 / 11 클러스터) |
| [vit_qsp_model.svg](vit_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [vit_qsp_model.png](vit_qsp_model.png) | PNG 이미지 (150 dpi) |
| [vit_mrgsolve_model.R](vit_mrgsolve_model.R) | mrgsolve ODE 모델 (약 21 구획 / 5 치료 시나리오) |
| [vit_shiny_app.R](vit_shiny_app.R) | Shiny 대시보드 |
| [vit_references.md](vit_references.md) | 참고문헌 (약 32편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 약물 PK(룩솔리티닙 경구 흡수·혈장·피부 3구획, 아파멜라노타이드 SC·혈장 2구획) + 면역 구획(NKG2D 리간드·NKG2D 활성·CD8+ CTL·Treg·IFN-γ·CXCL10·pSTAT1) + 멜라닌세포·MITF·멜라닌 합성·모낭 저장소 + 복합 지표(염증 지수·VASI·누적 색소 재침착)로 총 약 21개 구획
- **주요 치료 시나리오**: (1) 위약(자연 진행), (2) 룩솔리티닙 크림 BID, (3) 룩솔리티닙 크림 QD, (4) 룩솔리티닙 경구 10 mg BID, (5) 아파멜라노타이드 + NB-UVB 병용
- **보정/근거**: TRuE-V1/2(룩솔리티닙 외용 24주 F-VASI 반응), RECAL(아파멜라노타이드 + NB-UVB), VASI 점수 역학 데이터 참조

## Shiny 대시보드 (Dashboard)

환자 프로파일(체중·VASI 기저치·병변 분포·Fitzpatrick 피부 타입 설정) · PK 시각화(룩솔리티닙 피부·혈장 농도-시간 곡선) · 면역 PD 지표(IFN-γ·CXCL10·pSTAT1·CD8+ T세포 시계열) · 임상 엔드포인트(VASI 변화·색소 재침착 면적 비율) · 치료 시나리오 비교(5개 요법 장기 멜라닌세포 회복) · 바이오마커 패널(멜라닌 함량·모낭 저장소·Treg 수준) 등 6개 이상의 탭으로 구성된다.

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("vit_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("vit_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg vit_qsp_model.dot -o vit_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [vit_references.md](vit_references.md) 참조 (약 32편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
