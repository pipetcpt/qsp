# qsp

## mrgsolve

- <https://vantage-research.net/qsp-in-r/>
- gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve
    - <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/>
    - <https://github.com/Genentech/gPKPDviz/>
    



## iqrtools

- <https://www.intiquan.com/acop2019_qsp/>

---

## QSP Disease Models

| # | 질환명 | 디렉토리 | 기전 요약 | Mechanistic Map | 날짜 |
|---|--------|----------|-----------|-----------------|------|
| 1 | **폐동맥 고혈압 (PAH)**<br>Pulmonary Arterial Hypertension | [`폐동맥-고혈압/`](폐동맥-고혈압/) | ET-1(내피세포 수축인자) / NO-sGC-cGMP / PGI2-cAMP 세 핵심 경로의 불균형으로 폐동맥 혈관 수축 및 리모델링 유발; BMPR2 돌연변이로 항증식 신호 소실; 우심실 후부하 증가 → RV 비후/확장 → 우심부전 | [![PAH Map](폐동맥-고혈압/pah_qsp_map.png)](폐동맥-고혈압/pah_qsp_map.svg) | 2026-06-16 |

### 폐동맥 고혈압 (PAH) 상세

| 항목 | 내용 |
|------|------|
| **질환 분류** | 만성질환 / 심혈관-폐 |
| **ICD-10** | I27.0 (원발성 폐동맥 고혈압) |
| **주요 병태생리** | 내피세포 기능 부전 → ET-1↑ / NO↓ / PGI2↓; PASMC 증식·항세포사멸; 혈관 리모델링(내막비후·중막비대·외막섬유화·총상병변); 혈관 내 혈전 형성; 우심실 비후 → 탈보상 |
| **핵심 경로 수** | 14개 서브경로 (ET-1, NO/sGC/cGMP, PGI2/cAMP, 세로토닌, 성장인자, TGF-β/BMP, 염증, 산화스트레스, RhoA/ROCK, RAAS, 혈관 리모델링, 혈역학, RV 역학, 임상 엔드포인트) |
| **노드 수 (DOT)** | 130+ 노드, 160+ 엣지 |
| **약물 PK** | Sildenafil (2-cmt oral) · Bosentan (2-cmt oral) · Treprostinil (1-cmt SC) |
| **약물 PD** | ET-1 경로 → ERA (보센탄/앰브리센탄/마시텐탄) ; cGMP 경로 → PDE5i (실데나필/타달라필) + sGC 자극제(리오시구아트) ; cAMP 경로 → PGI2 유사체(에포프로스테놀/트레프로스티닐/일로프로스트) + IP 수용체 작용제(셀렉시파그/MRE-269) |
| **임상 엔드포인트** | PVR, mPAP, CO, Ees/Ea (RV-PA coupling), NT-proBNP, 6MWD, WHO FC |
| **파일 목록** | [DOT](폐동맥-고혈압/pah_qsp_map.dot) · [SVG](폐동맥-고혈압/pah_qsp_map.svg) · [PNG](폐동맥-고혈압/pah_qsp_map.png) · [R 모델](폐동맥-고혈압/pah_mrgsolve.R) · [Shiny 앱](폐동맥-고혈압/shiny_app.R) · [참고문헌](폐동맥-고혈압/references.md) |
| **참고문헌 수** | 43편 (PubMed 링크 포함) |
| **주요 임상시험** | SUPER-1 (sildenafil), PATENT-1 (riociguat), SERAPHIN (macitentan), ARIES-1/2 (ambrisentan), GRIPHON (selexipag), AMBITION (combination) |
