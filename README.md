# qsp

## mrgsolve

- <https://vantage-research.net/qsp-in-r/>
- gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve
    - <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/>
    - <https://github.com/Genentech/gPKPDviz/>
    



## iqrtools

- <https://www.intiquan.com/acop2019_qsp/>

---

## QSP Disease Model Library

매일 Claude Code Routine이 추가하는 질환별 QSP 모델 라이브러리입니다.  
각 모델은 기계론적 지도(.dot/.svg/.png) · mrgsolve ODE 모델(.R) · Shiny 앱(.R) · 참고문헌(.md)으로 구성됩니다.

| 날짜 | 분류 | 질환명 | 디렉토리 | 모델 미리보기 |
|------|------|--------|----------|--------------|
| 2026-06-19 | 만성질환 | 요로결석 (만성 재발성) | [urolithiasis/](urolithiasis/) | [![URI QSP](urolithiasis/uri_qsp_model.png)](urolithiasis/uri_qsp_model.svg) |

---

## 요로결석 (만성 재발성, Chronic Recurrent Urolithiasis)

### 개요

요로결석은 미국 성인 인구의 약 8.8%(남성 10.6%, 여성 7.1%)에서 발생하며, 5년 재발률은 50%에 달하는 만성 재발성 질환입니다. 발병기전은 소변 내 결석 형성 물질(칼슘, 수산염, 요산)의 과포화(supersaturation)와 결석 형성 억제 물질(구연산, 마그네슘, THP)의 부족이 복합적으로 작용합니다.

### 주요 병태생리 경로

| 경로 | 핵심 메커니즘 | 임상 이상 |
|------|------------|---------|
| 칼슘 항상성 | PTH↑ → 1,25(OH)₂D↑ → 장관 Ca²⁺ 흡수↑ | 고칼슘뇨증 (>300mg/day) |
| 수산염 대사 | AGXT 결핍(PH1) / 지방 흡수장애 → 장관 OX↑ | 고수산뇨증 (>45mg/day) |
| 요산 대사 | XO 과활성 / URAT1 변이 → UA↑ | 고요산뇨증 (>800mg/day) |
| 구연산 처리 | 대사산증 / RTA → 신세뇨관 구연산 재흡수↑ | 저구연산뇨증 (<320mg/day) |
| 소변 과포화 | Ca × OX / Ksp > 1 → 핵화(nucleation) | CaOx SS > 1.0 |
| Randall's 플라크 | 상피하 인회석 침착 → 결석 핵 형성 nidus | 결석 성장 |

### 약물 PK/PD 파라미터

| 약물 | 작용기전 | 주요 PK | 임상 효과 |
|------|---------|---------|---------|
| HCTZ 25mg/day | NCC 억제 → 원위세뇨관 Ca²⁺ 재흡수↑ | F=0.65, t½=6-15h, CL=18L/h | 요중 Ca 30-45% 감소 |
| Allopurinol 300mg/day | XO 기전불활성화(mechanism-based) | F=0.90, Oxypurinol t½=18-30h | UA 생성 40-60% 감소 |
| K-Citrate 60mEq/day | 요중 구연산↑ + 요 pH 알칼리화 | F=0.95, CL=15L/h | CaOx SS 50% 감소 |
| Tamsulosin 0.4mg/day | α₁A/D-수용체 차단 → 요관 이완 | F=0.90, t½=9-16h | 자연 배출율 28% 증가 |
| Lumasiran (siRNA) | GalNAc-간 표적 → HAGO1 억제 | 월 1회 피하주사 | 요중 OX 53% 감소 (PH1) |

### 모델 파일 목록

| 파일 | 설명 |
|------|------|
| [uri_qsp_model.dot](urolithiasis/uri_qsp_model.dot) | Graphviz 기계론적 지도 소스 (100+ 노드, 10 클러스터) |
| [uri_qsp_model.svg](urolithiasis/uri_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [uri_qsp_model.png](urolithiasis/uri_qsp_model.png) | PNG 래스터 이미지 (150 dpi) |
| [uri_mrgsolve_model.R](urolithiasis/uri_mrgsolve_model.R) | mrgsolve ODE 모델 (17 구획, 7 시나리오) |
| [uri_shiny_app.R](urolithiasis/uri_shiny_app.R) | Shiny 대시보드 (6탭: 환자/PK/요화학/결석위험/시나리오비교/바이오마커) |
| [uri_references.md](urolithiasis/uri_references.md) | 참고문헌 38편 (PubMed 링크 포함) |

### 주요 치료 시나리오 (mrgsolve 시뮬레이션)

1. **미치료 CaOx 결석 형성자**: 5년간 결석 성장 및 GFR 저하 추적
2. **HCTZ 25mg/day**: 고칼슘뇨증 환자에서 요중 Ca 30-45% 감소
3. **K-Citrate 60mEq/day**: 저구연산뇨증/UA 결석에서 CaOx SS 50% 감소
4. **Allopurinol 300mg/day**: 대사증후군/고요산뇨증에서 UA 생성 억제
5. **Lumasiran (PH1)**: 원발성 고수산뇨증에서 요중 OX 53% 감소
6. **생활습관 + HCTZ + K-Citrate 병용**: 종합 치료에서 결석 성장 억제
7. **MetSyn + UA 결석 + Allopurinol + K-Citrate**: 복합 요산 결석 관리
