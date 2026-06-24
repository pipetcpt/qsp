# QSP 질환 모델 라이브러리 (QSP Disease Model Library)

> 매일 **Claude Code Routine(CCR)** 이 질환 하나를 선택해 **정량적 시스템 약리학(Quantitative Systems Pharmacology, QSP)** 모델을 처음부터 끝까지 구축하고 `main`에 직접 커밋하는, **살아 있는(living) 오픈 모델 라이브러리**입니다.

![models](https://img.shields.io/badge/models-139-blue) ![framework](https://img.shields.io/badge/QSP-mrgsolve%20%C2%B7%20Shiny%20%C2%B7%20Graphviz-success) ![automation](https://img.shields.io/badge/built%20by-Claude%20Code%20Routine-orange)

현재 **139개 질환**에 대한 완성된 QSP 모델이 수록되어 있으며, 각 모델은 ①기계론적 지도, ②mrgsolve ODE 모델, ③Shiny 대시보드, ④참고문헌의 네 가지 산출물로 구성됩니다. 아래 [모델 갤러리](#-모델-갤러리-model-gallery)에서 전체 목록을 확인할 수 있습니다.
---

## 1. 프로젝트 소개 (Overview)

이 저장소는 사람이 한 번에 설계한 정적인 모델 모음이 아니라, **자동화된 AI 에이전트가 매일 한 편씩 새로운 질환 모델을 추가하며 성장하는 라이브러리**입니다. 각 세션에서 Claude Code는 다음을 수행합니다.

1. 아직 다루지 않은 질환을 선택한다.
2. 최신 문헌과 임상시험 데이터를 바탕으로 질환의 **기계론적 병태생리 지도**를 그린다.
3. 약물 PK/PD와 질환 진행을 연결하는 **mrgsolve 기반 ODE 모델**을 작성한다.
4. 파라미터를 탐색·비교할 수 있는 **Shiny 대시보드**를 만든다.
5. 모든 가정과 파라미터의 근거가 되는 **참고문헌**을 정리한다.
6. 위 산출물을 커밋·푸시한다.

목표는 다양한 치료 영역에 걸쳐 **재현 가능하고, 투명하며, 교육적으로 활용 가능한 QSP 모델의 참조 컬렉션**을 구축하는 것입니다.

## 2. QSP란 무엇인가 (What is QSP?)

**정량적 시스템 약리학(QSP)** 은 시스템 생물학과 약동/약력학(PK/PD)을 결합하여, **약물–표적–경로–질환–환자**로 이어지는 인과 사슬을 수학적(주로 상미분방정식, ODE)으로 표현하는 모델링 분야입니다. 단순한 통계적 용량–반응 곡선을 넘어, **"왜 그리고 어떻게" 약효와 독성이 나타나는지**를 기전 수준에서 기술합니다.

QSP 모델은 다음과 같은 질문에 답하려 합니다.

- 이 표적을 조절하면 하류 경로와 임상 바이오마커가 어떻게 움직이는가?
- 작용기전이 다른 약물을 병용하면 어떤 시너지/상쇄가 나타나는가?
- 어떤 환자 아형(유전형·중증도·동반질환)에서 반응이 달라지는가?
- 어떤 용량·투여 간격이 효능과 안전성의 균형을 최적화하는가?

## 3. 신약개발에서의 중요성 (Why QSP Matters in Drug Development)

QSP는 규제기관이 권장하는 **모델 기반 신약개발(Model-Informed Drug Development, MIDD)** 패러다임의 핵심 도구로 자리 잡았습니다. 신약개발 전 주기에서 다음과 같은 가치를 제공합니다.

| 단계 | QSP의 기여 |
|------|------------|
| **표적 발굴·검증** | 경로 수준 시뮬레이션으로 표적 조절이 질환 표현형에 미치는 영향을 사전 평가하고, 우회 경로·내성 기전을 예측 |
| **작용기전(MoA) 규명** | in vitro/in vivo 데이터를 통합해 약효의 인과 구조를 정량화하고 바이오마커 전략을 수립 |
| **용량·용법 선택** | First-in-human 용량, 적정 용법, 치료 범위를 기전 기반으로 예측하여 임상 1상 설계를 합리화 |
| **임상시험 설계** | 환자 아집단별 반응을 시뮬레이션해 환자 선정·계층화·엔드포인트·샘플 크기 최적화 |
| **중개연구(Translation)** | 동물–인간, 성인–소아, 건강인–환자 간 외삽으로 종간/집단간 차이를 정량화 |
| **병용요법 전략** | 서로 다른 표적의 조합 효과(시너지/상쇄)를 가상 환자군에서 탐색 |
| **규제 소통** | FDA·EMA의 MIDD 프레임워크에서 근거 패키지로 활용되어 임상시험 면제·라벨 확장을 뒷받침 |

신약개발 실패의 상당수는 **효능 부족**과 **예상치 못한 안전성 문제**에서 비롯됩니다. QSP는 후보물질이 임상에 진입하기 전에 *in silico* 로 가설을 검증하여 **개발 후기 실패(late-stage attrition)를 줄이고**, 실패할 프로그램은 더 일찍·더 싸게 중단(fail fast)하도록 돕습니다. 본 라이브러리는 이러한 QSP 워크플로를 다양한 질환에 적용한 **재현 가능한 교육·연구용 예제**를 제공합니다.

## 4. 각 모델의 구성 (Four Deliverables per Model)

모든 질환 디렉토리는 동일한 4종 산출물 체계를 따릅니다.

| 산출물 | 파일 | 설명 |
|--------|------|------|
| 🗺️ **기계론적 지도** | `*_qsp_model.dot` / `.svg` / `.png` | Graphviz로 그린 병태생리·약물 작용 네트워크 (100+ 노드, 8+ 클러스터) |
| ⚙️ **mrgsolve 모델** | `*_mrgsolve_model.R` | 약물 PK + 질환 PD를 연결하는 ODE 모델 (15+ 구획, 5+ 치료 시나리오, 임상시험 보정 메모) |
| 📊 **Shiny 대시보드** | `*_shiny_app.R` | 환자 프로파일·PK·PD·임상 엔드포인트·시나리오 비교·바이오마커를 탐색하는 인터랙티브 앱 (6+ 탭) |
| 📚 **참고문헌** | `*_references.md` | 섹션별로 분류된 30+ PubMed 인용 |

각 디렉토리에는 위 산출물을 요약한 **개별 `README.md`** 가 포함되어, 해당 질환의 개요·핵심 경로·모델 사양·실행 방법을 한눈에 볼 수 있습니다.

## 5. 작업 방식 (How It Is Built)

- **자동 생성**: Claude Code Routine이 매일 1개 질환을 선택하여 4종 산출물과 디렉토리 README를 작성한 뒤 커밋·푸시합니다.
- **중복 방지**: 이미 존재하는 질환/디렉토리는 건너뛰고, 카테고리를 번갈아 선택합니다.
- **품질 기준**: 기계론적 지도(100+ 노드/8+ 클러스터), mrgsolve(15+ 구획/5+ 시나리오), Shiny(6+ 탭), 참고문헌(30+편)을 최소 기준으로 합니다.
- **명명 규칙**: 디렉토리는 소문자-하이픈 영문명, 파일은 `<약어>_qsp_model.*`, `<약어>_mrgsolve_model.R`, `<약어>_shiny_app.R`, `<약어>_references.md` 형식을 따릅니다.

## 6. 기술 스택 (Technology Stack)

| 도구 | 용도 |
|------|------|
| **Graphviz** (`dot`) | 기계론적 지도 렌더링 (`.dot` → `.svg`/`.png`) |
| **mrgsolve** (R) | ODE 기반 PK/PD/QSP 시뮬레이션 |
| **Shiny** (R) | 인터랙티브 시뮬레이션 대시보드 |
| **Claude Code Routine** | 매일 자동 모델 생성·문서화·커밋 |

## 7. 사용 방법 (Usage)

```bash
# 1) 기계론적 지도 렌더링 (Graphviz 필요)
dot -Tsvg <disease>/<abbr>_qsp_model.dot -o <abbr>_qsp_model.svg
dot -Tpng -Gdpi=150 <disease>/<abbr>_qsp_model.dot -o <abbr>_qsp_model.png
```

```r
# 2) mrgsolve 모델 실행 (R)
install.packages(c("mrgsolve", "dplyr", "ggplot2"))
library(mrgsolve)
mod <- mread("<disease>/<abbr>_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)

# 3) Shiny 대시보드 실행
install.packages("shiny")
shiny::runApp("<disease>/<abbr>_shiny_app.R")
```

## 8. 디렉토리 구조 (Repository Layout)

```
qsp/
├── README.md                     # 본 문서 (전체 모델 갤러리)
├── CLAUDE.md                     # 라이브러리 운영·생성 지침
├── <disease>/                    # 질환별 디렉토리 (총 118개)
│   ├── README.md                 # 질환별 요약 문서
│   ├── <abbr>_qsp_model.dot      # 기계론적 지도 소스
│   ├── <abbr>_qsp_model.svg/.png # 렌더링 이미지
│   ├── <abbr>_mrgsolve_model.R   # mrgsolve ODE 모델
│   ├── <abbr>_shiny_app.R        # Shiny 대시보드
│   └── <abbr>_references.md      # 참고문헌
└── ...
```

---

## 📚 모델 갤러리 (Model Gallery)

전체 **134개** QSP 모델입니다. 모델명을 클릭하면 해당 디렉토리로, 그림을 클릭하면 확대 가능한 SVG 지도로 이동합니다. 각 행의 링크에서 기계론적 지도(🗺️), mrgsolve 모델(⚙️), 참고문헌(📚), 상세 README(📄)에 바로 접근할 수 있습니다.

**분류별 모델 수**: 내분비·대사 21 · 소화기·간담도 17 · 자가면역·류마티스 13 · 심혈관 11 · 신장·비뇨 11 · 신경 11 · 호흡기 8 · 혈관염 8 · 소아 혈관염 1 · 피부 6 · 혈액 5 · 종양 3 · 정신·신경 2 · 감염 2 · 신경근육 1

| # | 분류 | 모델 | 미리보기 | 요약 및 링크 |
|---|------|------|----------|--------------|
| 1 | 심혈관 | [**복부 대동맥류**<br><sub>Abdominal Aortic Aneurysm · AAA</sub>](abdominal-aortic-aneurysm/) | <a href="abdominal-aortic-aneurysm/aaa_qsp_model.svg"><img src="abdominal-aortic-aneurysm/aaa_qsp_model.png" width="190" alt="AAA"></a> | 대동맥 중막 ECM 분해(MMP)·만성 염증으로 인한 진행성 확장과 파열 위험. 베타차단·독시사이클린·항염증 표적을 모델링.<br>[🗺️ 지도](abdominal-aortic-aneurysm/aaa_qsp_model.svg) · [⚙️ mrgsolve](abdominal-aortic-aneurysm/aaa_mrgsolve_model.R) · [📚 문헌](abdominal-aortic-aneurysm/aaa_references.md) · [📄 README](abdominal-aortic-aneurysm/README.md) |
| 2 | 내분비·대사 | [**말단비대증**<br><sub>Acromegaly · ACRO</sub>](acromegaly/) | <a href="acromegaly/acro_qsp_model.svg"><img src="acromegaly/acro_qsp_model.png" width="190" alt="ACRO"></a> | 뇌하수체 선종의 GH 과다분비 → IGF-1 상승. 소마토스타틴 유사체·페그비소만트·도파민작용제 PK/PD.<br>[🗺️ 지도](acromegaly/acro_qsp_model.svg) · [⚙️ mrgsolve](acromegaly/acro_mrgsolve_model.R) · [📚 문헌](acromegaly/acro_references.md) · [📄 README](acromegaly/README.md) |
| 3 | 내분비·대사 | [**애디슨병 (원발성 부신부전)**<br><sub>Addison's Disease · ADD</sub>](addisons-disease/) | <a href="addisons-disease/add_qsp_model.svg"><img src="addisons-disease/add_qsp_model.png" width="190" alt="ADD"></a> | 부신피질 자가면역 파괴로 코르티솔·알도스테론 결핍. HPA축·히드로코르티손/플루드로코르티손 보충 시뮬레이션.<br>[🗺️ 지도](addisons-disease/add_qsp_model.svg) · [⚙️ mrgsolve](addisons-disease/add_mrgsolve_model.R) · [📚 문헌](addisons-disease/add_references.md) · [📄 README](addisons-disease/README.md) |
| 4 | 신장·비뇨 | [**상염색체 우성 다낭신 (ADPKD)**<br><sub>Autosomal Dominant PKD · ADPKD</sub>](adpkd/) | <a href="adpkd/adpkd_qsp_model.svg"><img src="adpkd/adpkd_qsp_model.png" width="190" alt="ADPKD"></a> | PKD1/2 변이 → cAMP 증가로 낭종 성장·신기능 저하. 톨밥탄(V2 수용체 길항)의 TKV/eGFR 효과.<br>[🗺️ 지도](adpkd/adpkd_qsp_model.svg) · [⚙️ mrgsolve](adpkd/adpkd_mrgsolve_model.R) · [📚 문헌](adpkd/adpkd_references.md) · [📄 README](adpkd/README.md) |
| 5 | 자가면역·류마티스 | [**성인형 스틸병**<br><sub>Adult-Onset Still's Disease · AOSD</sub>](adult-onset-stills-disease/) | <a href="adult-onset-stills-disease/aosd_qsp_model.svg"><img src="adult-onset-stills-disease/aosd_qsp_model.png" width="190" alt="AOSD"></a> | IL-1/IL-6/IL-18 매개 자가염증, 발열·관절염·페리틴 급상승. 아나킨라·카나키누맙·토실리주맙 표적.<br>[🗺️ 지도](adult-onset-stills-disease/aosd_qsp_model.svg) · [⚙️ mrgsolve](adult-onset-stills-disease/aosd_mrgsolve_model.R) · [📚 문헌](adult-onset-stills-disease/aosd_references.md) · [📄 README](adult-onset-stills-disease/README.md) |
| 6 | 피부 | [**원형 탈모증**<br><sub>Alopecia Areata · AA</sub>](alopecia-areata/) | <a href="alopecia-areata/aa_qsp_model.svg"><img src="alopecia-areata/aa_qsp_model.png" width="190" alt="AA"></a> | 모낭 면역특권 붕괴, IFN-γ/IL-15 JAK-STAT 신호로 모발 손실. 바리시티닙 등 JAK 억제제 반응.<br>[🗺️ 지도](alopecia-areata/aa_qsp_model.svg) · [⚙️ mrgsolve](alopecia-areata/aa_mrgsolve_model.R) · [📚 문헌](alopecia-areata/aa_references.md) · [📄 README](alopecia-areata/README.md) |
| 7 | 신경 | [**알츠하이머병**<br><sub>Alzheimer's Disease · AD</sub>](alzheimers-disease/) | <a href="alzheimers-disease/ad_qsp_model.svg"><img src="alzheimers-disease/ad_qsp_model.png" width="190" alt="AD"></a> | 아밀로이드-β 침착·타우 신경섬유엉킴·신경염증에 의한 진행성 인지저하. 항아밀로이드 항체(레카네맙/도나네맙)·콜린에스터분해효소 억제제.<br>[🗺️ 지도](alzheimers-disease/ad_qsp_model.svg) · [⚙️ mrgsolve](alzheimers-disease/ad_mrgsolve_model.R) · [📚 문헌](alzheimers-disease/ad_references.md) · [📄 README](alzheimers-disease/README.md) |
| 8 | 신경 | [**근위축성 측삭경화증 (ALS)**<br><sub>Amyotrophic Lateral Sclerosis · ALS</sub>](amyotrophic-lateral-sclerosis/) | <a href="amyotrophic-lateral-sclerosis/als_qsp_model.svg"><img src="amyotrophic-lateral-sclerosis/als_qsp_model.png" width="190" alt="ALS"></a> | 상·하위 운동신경세포의 진행성 변성(SOD1/TDP-43/C9orf72). 릴루졸·에다라본·토퍼센(SOD1).<br>[🗺️ 지도](amyotrophic-lateral-sclerosis/als_qsp_model.svg) · [⚙️ mrgsolve](amyotrophic-lateral-sclerosis/als_mrgsolve_model.R) · [📚 문헌](amyotrophic-lateral-sclerosis/als_references.md) · [📄 README](amyotrophic-lateral-sclerosis/README.md) |
| 9 | 자가면역·류마티스 | [**강직성 척추염**<br><sub>Ankylosing Spondylitis · AS</sub>](ankylosing-spondylitis/) | <a href="ankylosing-spondylitis/as_qsp_model.svg"><img src="ankylosing-spondylitis/as_qsp_model.png" width="190" alt="AS"></a> | HLA-B27·IL-23/IL-17·TNF 축의 부착부염·골증식. TNF/IL-17 억제제의 BASDAI 효과.<br>[🗺️ 지도](ankylosing-spondylitis/as_qsp_model.svg) · [⚙️ mrgsolve](ankylosing-spondylitis/as_mrgsolve_model.R) · [📚 문헌](ankylosing-spondylitis/as_references.md) · [📄 README](ankylosing-spondylitis/README.md) |
| 10 | 자가면역·류마티스 | [**항인지질항체 증후군**<br><sub>Antiphospholipid Syndrome · APS</sub>](antiphospholipid-syndrome/) | <a href="antiphospholipid-syndrome/aps_qsp_model.svg"><img src="antiphospholipid-syndrome/aps_qsp_model.png" width="190" alt="APS"></a> | 항인지질항체에 의한 동·정맥 혈전 및 임신 합병증. 항응고(와파린/헤파린)·히드록시클로로퀸.<br>[🗺️ 지도](antiphospholipid-syndrome/aps_qsp_model.svg) · [⚙️ mrgsolve](antiphospholipid-syndrome/aps_mrgsolve_model.R) · [📚 문헌](antiphospholipid-syndrome/aps_references.md) · [📄 README](antiphospholipid-syndrome/README.md) |
| 11 | 내분비·대사 | [**무증상 고요산혈증**<br><sub>Asymptomatic Hyperuricemia · AHU</sub>](asymptomatic-hyperuricemia/) | <a href="asymptomatic-hyperuricemia/ahy_qsp_model.svg"><img src="asymptomatic-hyperuricemia/ahy_qsp_model.png" width="190" alt="AHU"></a> | 요산 과포화(무증상)와 통풍·신질환 진행 위험. 생활습관 및 선택적 요산저하(잔틴산화효소 억제).<br>[🗺️ 지도](asymptomatic-hyperuricemia/ahy_qsp_model.svg) · [⚙️ mrgsolve](asymptomatic-hyperuricemia/ahy_mrgsolve_model.R) · [📚 문헌](asymptomatic-hyperuricemia/ahy_references.md) · [📄 README](asymptomatic-hyperuricemia/README.md) |
| 12 | 피부 | [**아토피 피부염**<br><sub>Atopic Dermatitis · AD</sub>](atopic-dermatitis/) | <a href="atopic-dermatitis/ad_qsp_model.svg"><img src="atopic-dermatitis/ad_qsp_model.png" width="190" alt="AD"></a> | 피부장벽 기능장애·Th2(IL-4/IL-13/IL-31) 염증. 두필루맙·JAK 억제제·국소 칼시뉴린억제제.<br>[🗺️ 지도](atopic-dermatitis/ad_qsp_model.svg) · [⚙️ mrgsolve](atopic-dermatitis/ad_mrgsolve_model.R) · [📚 문헌](atopic-dermatitis/ad_references.md) · [📄 README](atopic-dermatitis/README.md) |
| 13 | 심혈관 | [**심방세동**<br><sub>Atrial Fibrillation · AF</sub>](atrial-fibrillation/) | <a href="atrial-fibrillation/af_qsp_model.svg"><img src="atrial-fibrillation/af_qsp_model.png" width="190" alt="AF"></a> | 심방 전기적·구조적 리모델링과 재진입 회로. 율동/심박수 조절 및 NOAC 항응고 전략.<br>[🗺️ 지도](atrial-fibrillation/af_qsp_model.svg) · [⚙️ mrgsolve](atrial-fibrillation/af_mrgsolve_model.R) · [📚 문헌](atrial-fibrillation/af_references.md) · [📄 README](atrial-fibrillation/README.md) |
| 14 | 신경 | [**자가면역 뇌염**<br><sub>Autoimmune Encephalitis · AIE</sub>](autoimmune-encephalitis/) | <a href="autoimmune-encephalitis/aie_qsp_model.svg"><img src="autoimmune-encephalitis/aie_qsp_model.png" width="190" alt="AIE"></a> | 항NMDAR 등 신경표면 항체에 의한 뇌염. 스테로이드·IVIG·혈장교환·리툭시맙 면역치료.<br>[🗺️ 지도](autoimmune-encephalitis/aie_qsp_model.svg) · [⚙️ mrgsolve](autoimmune-encephalitis/aie_mrgsolve_model.R) · [📚 문헌](autoimmune-encephalitis/aie_references.md) · [📄 README](autoimmune-encephalitis/README.md) |
| 15 | 혈액 | [**자가면역 용혈성 빈혈**<br><sub>Autoimmune Hemolytic Anemia · AIHA</sub>](autoimmune-hemolytic-anemia/) | <a href="autoimmune-hemolytic-anemia/aiha_qsp_model.svg"><img src="autoimmune-hemolytic-anemia/aiha_qsp_model.png" width="190" alt="AIHA"></a> | 온형 IgG·한랭 IgM 자가항체에 의한 적혈구 파괴. 스테로이드·리툭시맙·수티림맙(항C1s).<br>[🗺️ 지도](autoimmune-hemolytic-anemia/aiha_qsp_model.svg) · [⚙️ mrgsolve](autoimmune-hemolytic-anemia/aiha_mrgsolve_model.R) · [📚 문헌](autoimmune-hemolytic-anemia/aiha_references.md) · [📄 README](autoimmune-hemolytic-anemia/README.md) |
| 16 | 소화기·간담도 | [**자가면역 간염**<br><sub>Autoimmune Hepatitis · AIH</sub>](autoimmune-hepatitis/) | <a href="autoimmune-hepatitis/aih_qsp_model.svg"><img src="autoimmune-hepatitis/aih_qsp_model.png" width="190" alt="AIH"></a> | T세포 매개 간세포 손상과 자가항체·고감마글로불린혈증. 스테로이드·아자티오프린 유도/유지.<br>[🗺️ 지도](autoimmune-hepatitis/aih_qsp_model.svg) · [⚙️ mrgsolve](autoimmune-hepatitis/aih_mrgsolve_model.R) · [📚 문헌](autoimmune-hepatitis/aih_references.md) · [📄 README](autoimmune-hepatitis/README.md) |
| 17 | 소화기·간담도 | [**자가면역 췌장염**<br><sub>Autoimmune Pancreatitis · AIP</sub>](autoimmune-pancreatitis/) | <a href="autoimmune-pancreatitis/aip_qsp_model.svg"><img src="autoimmune-pancreatitis/aip_qsp_model.png" width="190" alt="AIP"></a> | IgG4 관련 림프형질세포 침윤·섬유염증. 스테로이드 관해유도·리툭시맙 재발관리.<br>[🗺️ 지도](autoimmune-pancreatitis/aip_qsp_model.svg) · [⚙️ mrgsolve](autoimmune-pancreatitis/aip_mrgsolve_model.R) · [📚 문헌](autoimmune-pancreatitis/aip_references.md) · [📄 README](autoimmune-pancreatitis/README.md) |
| 18 | 내분비·대사 | [**자가면역 다발내분비병증 (APECED)**<br><sub>Autoimmune Polyendocrinopathy (APS-1) · APS-1</sub>](autoimmune-polyendocrinopathy/) | <a href="autoimmune-polyendocrinopathy/aps_qsp_model.svg"><img src="autoimmune-polyendocrinopathy/aps_qsp_model.png" width="190" alt="APS-1"></a> | AIRE 변이로 중추관용 소실 → 다장기 내분비 자가면역. 호르몬 보충 및 면역조절 모델.<br>[🗺️ 지도](autoimmune-polyendocrinopathy/aps_qsp_model.svg) · [⚙️ mrgsolve](autoimmune-polyendocrinopathy/aps_mrgsolve_model.R) · [📚 문헌](autoimmune-polyendocrinopathy/aps_references.md) · [📄 README](autoimmune-polyendocrinopathy/README.md) |
| 19 | 혈관염 | [**베체트병**<br><sub>Behçet's Disease · BD</sub>](behcet-disease/) | <a href="behcet-disease/bd_qsp_model.svg"><img src="behcet-disease/bd_qsp_model.png" width="190" alt="BD"></a> | 구강/생식기 궤양·포도막염, 호중구·IL-1 매개 가변혈관 염증. 콜히친·TNF 억제·아프레밀라스트.<br>[🗺️ 지도](behcet-disease/bd_qsp_model.svg) · [⚙️ mrgsolve](behcet-disease/bd_mrgsolve_model.R) · [📚 문헌](behcet-disease/bd_references.md) · [📄 README](behcet-disease/README.md) |
| 20 | 신장·비뇨 | [**양성 전립선 비대증**<br><sub>Benign Prostatic Hyperplasia · BPH</sub>](benign-prostatic-hyperplasia/) | <a href="benign-prostatic-hyperplasia/bph_qsp_model.svg"><img src="benign-prostatic-hyperplasia/bph_qsp_model.png" width="190" alt="BPH"></a> | DHT 매개 전립선 증식과 평활근 긴장에 의한 하부요로증상. α차단제·5α환원효소 억제제.<br>[🗺️ 지도](benign-prostatic-hyperplasia/bph_qsp_model.svg) · [⚙️ mrgsolve](benign-prostatic-hyperplasia/bph_mrgsolve_model.R) · [📚 문헌](benign-prostatic-hyperplasia/bph_references.md) · [📄 README](benign-prostatic-hyperplasia/README.md) |
| 21 | 종양 | [**유방암**<br><sub>Breast Cancer · BC</sub>](breast-cancer/) | <a href="breast-cancer/bc_qsp_model.svg"><img src="breast-cancer/bc_qsp_model.png" width="190" alt="BC"></a> | ER+/HER2+/TNBC 아형별 증식 신호. 내분비요법·CDK4/6 억제제·항HER2·면역항암제.<br>[🗺️ 지도](breast-cancer/bc_qsp_model.svg) · [⚙️ mrgsolve](breast-cancer/bc_mrgsolve_model.R) · [📚 문헌](breast-cancer/bc_references.md) · [📄 README](breast-cancer/README.md) |
| 22 | 호흡기 | [**기관지 천식**<br><sub>Bronchial Asthma · BA</sub>](bronchial-asthma/) | <a href="bronchial-asthma/ba_qsp_model.svg"><img src="bronchial-asthma/ba_qsp_model.png" width="190" alt="BA"></a> | Th2/호산구 기도 염증과 기관지 과민성. ICS·LABA 및 항IL-5/IL-4Rα 생물학제제.<br>[🗺️ 지도](bronchial-asthma/ba_qsp_model.svg) · [⚙️ mrgsolve](bronchial-asthma/ba_mrgsolve_model.R) · [📚 문헌](bronchial-asthma/ba_references.md) · [📄 README](bronchial-asthma/README.md) |
| 23 | 호흡기 | [**기관지 확장증**<br><sub>Bronchiectasis · BEX</sub>](bronchiectasis/) | <a href="bronchiectasis/bex_qsp_model.svg"><img src="bronchiectasis/bex_qsp_model.png" width="190" alt="BEX"></a> | 감염-염증 악순환(Cole vicious cycle)과 호중구 엘라스타제 기도 파괴. 기도청결·항생제·항염.<br>[🗺️ 지도](bronchiectasis/bex_qsp_model.svg) · [⚙️ mrgsolve](bronchiectasis/bex_mrgsolve_model.R) · [📚 문헌](bronchiectasis/bex_references.md) · [📄 README](bronchiectasis/README.md) |
| 24 | 피부 | [**수포성 유천포창**<br><sub>Bullous Pemphigoid · BP</sub>](bullous-pemphigoid/) | <a href="bullous-pemphigoid/bp_qsp_model.svg"><img src="bullous-pemphigoid/bp_qsp_model.png" width="190" alt="BP"></a> | 항BP180/BP230 항체에 의한 표피하 수포 형성. 스테로이드·리툭시맙·오말리주맙.<br>[🗺️ 지도](bullous-pemphigoid/bp_qsp_model.svg) · [⚙️ mrgsolve](bullous-pemphigoid/bp_mrgsolve_model.R) · [📚 문헌](bullous-pemphigoid/bp_references.md) · [📄 README](bullous-pemphigoid/README.md) |
| 25 | 소화기·간담도 | [**셀리악병**<br><sub>Celiac Disease · CD</sub>](celiac-disease/) | <a href="celiac-disease/cd_qsp_model.svg"><img src="celiac-disease/cd_qsp_model.png" width="190" alt="CD"></a> | 글루텐-tTG 면역반응(HLA-DQ2/8)으로 융모 위축. 글루텐 제거식 및 글루텐분해효소 신약.<br>[🗺️ 지도](celiac-disease/cd_qsp_model.svg) · [⚙️ mrgsolve](celiac-disease/cd_mrgsolve_model.R) · [📚 문헌](celiac-disease/cd_references.md) · [📄 README](celiac-disease/README.md) |
| 26 | 소화기·간담도 | [**담석증**<br><sub>Cholelithiasis · CHOL</sub>](cholelithiasis/) | <a href="cholelithiasis/chol_qsp_model.svg"><img src="cholelithiasis/chol_qsp_model.png" width="190" alt="CHOL"></a> | 담즙 콜레스테롤 과포화·담낭 정체에 의한 결석 형성. 우르소데옥시콜산 용해 요법.<br>[🗺️ 지도](cholelithiasis/chol_qsp_model.svg) · [⚙️ mrgsolve](cholelithiasis/chol_mrgsolve_model.R) · [📚 문헌](cholelithiasis/chol_references.md) · [📄 README](cholelithiasis/README.md) |
| 27 | 소화기·간담도 | [**만성 위염**<br><sub>Chronic Gastritis · CGAST</sub>](chronic-gastritis/) | <a href="chronic-gastritis/cgast_qsp_model.svg"><img src="chronic-gastritis/cgast_qsp_model.png" width="190" alt="CGAST"></a> | H. pylori-Correa 연쇄(위축→장상피화생→이형성). 제균요법·PPI의 점막 회복 효과.<br>[🗺️ 지도](chronic-gastritis/cgast_qsp_model.svg) · [⚙️ mrgsolve](chronic-gastritis/cgast_mrgsolve_model.R) · [📚 문헌](chronic-gastritis/cgast_references.md) · [📄 README](chronic-gastritis/README.md) |
| 28 | 소화기·간담도 | [**만성 B형 간염**<br><sub>Chronic Hepatitis B · CHB</sub>](chronic-hepatitis-b/) | <a href="chronic-hepatitis-b/chb_qsp_model.svg"><img src="chronic-hepatitis-b/chb_qsp_model.png" width="190" alt="CHB"></a> | HBV cccDNA 지속과 면역 매개 간손상. 뉴클레오시드유사체·페그IFN의 바이러스 억제.<br>[🗺️ 지도](chronic-hepatitis-b/chb_qsp_model.svg) · [⚙️ mrgsolve](chronic-hepatitis-b/chb_mrgsolve_model.R) · [📚 문헌](chronic-hepatitis-b/chb_references.md) · [📄 README](chronic-hepatitis-b/README.md) |
| 29 | 내분비·대사 | [**만성 갑상선 기능 저하증**<br><sub>Chronic Hypothyroidism · HYPO</sub>](chronic-hypothyroidism/) | <a href="chronic-hypothyroidism/hypo_qsp_model.svg"><img src="chronic-hypothyroidism/hypo_qsp_model.png" width="190" alt="HYPO"></a> | 시상하부-뇌하수체-갑상선 축 저하(TSH↑/FT4↓). 레보티록신 보충의 용량-반응.<br>[🗺️ 지도](chronic-hypothyroidism/hypo_qsp_model.svg) · [⚙️ mrgsolve](chronic-hypothyroidism/hypo_mrgsolve_model.R) · [📚 문헌](chronic-hypothyroidism/hypo_references.md) · [📄 README](chronic-hypothyroidism/README.md) |
| 30 | 신장·비뇨 | [**만성 신부전**<br><sub>Chronic Kidney Disease · CKD</sub>](chronic-kidney-disease/) | <a href="chronic-kidney-disease/ckd_qsp_model.svg"><img src="chronic-kidney-disease/ckd_qsp_model.png" width="190" alt="CKD"></a> | 사구체 과여과·섬유화로 eGFR 진행성 저하. RAAS 억제·SGLT2i·피네레논 신보호.<br>[🗺️ 지도](chronic-kidney-disease/ckd_qsp_model.svg) · [⚙️ mrgsolve](chronic-kidney-disease/ckd_mrgsolve_model.R) · [📚 문헌](chronic-kidney-disease/ckd_references.md) · [📄 README](chronic-kidney-disease/README.md) |
| 31 | 종양 | [**만성 골수성 백혈병 (CML)**<br><sub>Chronic Myeloid Leukemia · CML</sub>](chronic-myeloid-leukemia/) | <a href="chronic-myeloid-leukemia/cml_qsp_model.svg"><img src="chronic-myeloid-leukemia/cml_qsp_model.png" width="190" alt="CML"></a> | BCR-ABL 융합 티로신키나아제 구동 백혈병. 이매티닙 등 TKI·T315I 내성·무치료 관해(TFR).<br>[🗺️ 지도](chronic-myeloid-leukemia/cml_qsp_model.svg) · [⚙️ mrgsolve](chronic-myeloid-leukemia/cml_mrgsolve_model.R) · [📚 문헌](chronic-myeloid-leukemia/cml_references.md) · [📄 README](chronic-myeloid-leukemia/README.md) |
| 32 | 소화기·간담도 | [**만성 췌장염**<br><sub>Chronic Pancreatitis · CP</sub>](chronic-pancreatitis/) | <a href="chronic-pancreatitis/cp_qsp_model.svg"><img src="chronic-pancreatitis/cp_qsp_model.png" width="190" alt="CP"></a> | 반복 손상 → 췌성상세포 활성·섬유화·외분비 부전. 효소보충·통증관리 전략.<br>[🗺️ 지도](chronic-pancreatitis/cp_qsp_model.svg) · [⚙️ mrgsolve](chronic-pancreatitis/cp_mrgsolve_model.R) · [📚 문헌](chronic-pancreatitis/cp_references.md) · [📄 README](chronic-pancreatitis/README.md) |
| 33 | 신장·비뇨 | [**만성 신우신염**<br><sub>Chronic Pyelonephritis · CPN</sub>](chronic-pyelonephritis/) | <a href="chronic-pyelonephritis/cpn_qsp_model.svg"><img src="chronic-pyelonephritis/cpn_qsp_model.png" width="190" alt="CPN"></a> | 반복 신장 감염·역류로 인한 신실질 흉터·신기능 저하. 항생제·역류 교정.<br>[🗺️ 지도](chronic-pyelonephritis/cpn_qsp_model.svg) · [⚙️ mrgsolve](chronic-pyelonephritis/cpn_mrgsolve_model.R) · [📚 문헌](chronic-pyelonephritis/cpn_references.md) · [📄 README](chronic-pyelonephritis/README.md) |
| 34 | 심혈관 | [**만성 정맥 부전**<br><sub>Chronic Venous Insufficiency · CVI</sub>](chronic-venous-insufficiency/) | <a href="chronic-venous-insufficiency/cvi_qsp_model.svg"><img src="chronic-venous-insufficiency/cvi_qsp_model.png" width="190" alt="CVI"></a> | 판막 부전·정맥압 상승·백혈구 트래핑에 의한 부종·궤양. 정맥활성약(MPFF)·압박요법.<br>[🗺️ 지도](chronic-venous-insufficiency/cvi_qsp_model.svg) · [⚙️ mrgsolve](chronic-venous-insufficiency/cvi_mrgsolve_model.R) · [📚 문헌](chronic-venous-insufficiency/cvi_references.md) · [📄 README](chronic-venous-insufficiency/README.md) |
| 35 | 신경 | [**만성 염증성 탈수초성 다발신경병증 (CIDP)**<br><sub>Chronic Inflammatory Demyelinating Polyneuropathy · CIDP</sub>](cidp/) | <a href="cidp/cidp_qsp_model.svg"><img src="cidp/cidp_qsp_model.png" width="190" alt="CIDP"></a> | 자가면역 말이집 파괴에 의한 진행/재발성 신경병증. IVIG·스테로이드·혈장교환.<br>[🗺️ 지도](cidp/cidp_qsp_model.svg) · [⚙️ mrgsolve](cidp/cidp_mrgsolve_model.R) · [📚 문헌](cidp/cidp_references.md) · [📄 README](cidp/README.md) |
| 36 | 호흡기 | [**만성 폐쇄성 폐질환 (COPD)**<br><sub>Chronic Obstructive Pulmonary Disease · COPD</sub>](copd/) | <a href="copd/copd_qsp.svg"><img src="copd/copd_qsp.png" width="190" alt="COPD"></a> | 흡연-산화스트레스·단백분해효소 불균형에 의한 기도/폐포 파괴. LABA/LAMA·ICS 기관지확장.<br>[🗺️ 지도](copd/copd_qsp.svg) · [⚙️ mrgsolve](copd/copd_mrgsolve_model.R) · [📚 문헌](copd/copd_references.md) · [📄 README](copd/README.md) |
| 37 | 소화기·간담도 | [**크론병**<br><sub>Crohn's Disease · CD</sub>](crohn-disease/) | <a href="crohn-disease/cd_qsp_model.svg"><img src="crohn-disease/cd_qsp_model.png" width="190" alt="CD"></a> | Th1/Th17 경유 전층성 장 염증과 장벽 투과성. 항TNF·항IL-12/23·항인테그린 생물학제제.<br>[🗺️ 지도](crohn-disease/cd_qsp_model.svg) · [⚙️ mrgsolve](crohn-disease/cd_mrgsolve_model.R) · [📚 문헌](crohn-disease/cd_references.md) · [📄 README](crohn-disease/README.md) |
| 38 | 호흡기 | [**낭성 섬유증**<br><sub>Cystic Fibrosis · CF</sub>](cystic-fibrosis/) | <a href="cystic-fibrosis/cf_qsp_model.svg"><img src="cystic-fibrosis/cf_qsp_model.png" width="190" alt="CF"></a> | CFTR 변이(ΔF508)에 의한 점액 점성↑·만성 기도 감염. CFTR 조절제(엘렉사/테자/이바카프토르, Trikafta).<br>[🗺️ 지도](cystic-fibrosis/cf_qsp_model.svg) · [⚙️ mrgsolve](cystic-fibrosis/cf_mrgsolve_model.R) · [📚 문헌](cystic-fibrosis/cf_references.md) · [📄 README](cystic-fibrosis/README.md) |
| 39 | 자가면역·류마티스 | [**피부근염**<br><sub>Dermatomyositis · DM</sub>](dermatomyositis/) | <a href="dermatomyositis/dm_qsp_model.svg"><img src="dermatomyositis/dm_qsp_model.png" width="190" alt="DM"></a> | I형 IFN·보체 매개 근육·피부 미세혈관병증. 스테로이드·IVIG·리툭시맙.<br>[🗺️ 지도](dermatomyositis/dm_qsp_model.svg) · [⚙️ mrgsolve](dermatomyositis/dm_mrgsolve_model.R) · [📚 문헌](dermatomyositis/dm_references.md) · [📄 README](dermatomyositis/README.md) |
| 40 | 심혈관 | [**확장성 심근병증**<br><sub>Dilated Cardiomyopathy · DCM</sub>](dilated-cardiomyopathy/) | <a href="dilated-cardiomyopathy/dcm_qsp_model.svg"><img src="dilated-cardiomyopathy/dcm_qsp_model.png" width="190" alt="DCM"></a> | 심실 확장·수축기능 저하·신경호르몬 활성화. GDMT(ARNI/BB/MRA/SGLT2i) 역리모델링.<br>[🗺️ 지도](dilated-cardiomyopathy/dcm_qsp_model.svg) · [⚙️ mrgsolve](dilated-cardiomyopathy/dcm_mrgsolve_model.R) · [📚 문헌](dilated-cardiomyopathy/dcm_references.md) · [📄 README](dilated-cardiomyopathy/README.md) |
| 41 | 소화기·간담도 | [**게실병**<br><sub>Diverticular Disease · DIV</sub>](diverticular-disease/) | <a href="diverticular-disease/div_qsp_model.svg"><img src="diverticular-disease/div_qsp_model.png" width="190" alt="DIV"></a> | 장벽 구조·미생물·저섬유식에 의한 게실/게실염. 섬유·항생제·항염 치료.<br>[🗺️ 지도](diverticular-disease/div_qsp_model.svg) · [⚙️ mrgsolve](diverticular-disease/div_mrgsolve_model.R) · [📚 문헌](diverticular-disease/div_references.md) · [📄 README](diverticular-disease/README.md) |
| 42 | 내분비·대사 | [**이상지질혈증**<br><sub>Dyslipidemia · DYSLIP</sub>](dyslipidemia/) | <a href="dyslipidemia/dyslip_qsp_model.svg"><img src="dyslipidemia/dyslip_qsp_model.png" width="190" alt="DYSLIP"></a> | 간 콜레스테롤 합성·LDL 수용체·역수송 항상성. 스타틴·에제티미브·PCSK9 억제제.<br>[🗺️ 지도](dyslipidemia/dyslip_qsp_model.svg) · [⚙️ mrgsolve](dyslipidemia/dyslip_mrgsolve_model.R) · [📚 문헌](dyslipidemia/dyslip_references.md) · [📄 README](dyslipidemia/README.md) |
| 43 | 혈관염 | [**호산구 육아종증 다발혈관염 (EGPA)**<br><sub>Eosinophilic GPA · EGPA</sub>](egpa/) | <a href="egpa/egpa_qsp_model.svg"><img src="egpa/egpa_qsp_model.png" width="190" alt="EGPA"></a> | IL-5 매개 호산구증가와 ANCA 연관 소혈관염. 스테로이드·메폴리주맙·벤랄리주맙.<br>[🗺️ 지도](egpa/egpa_qsp_model.svg) · [⚙️ mrgsolve](egpa/egpa_mrgsolve_model.R) · [📚 문헌](egpa/egpa_references.md) · [📄 README](egpa/README.md) |
| 44 | 내분비·대사 | [**자궁내막증**<br><sub>Endometriosis · ENDO</sub>](endometriosis/) | <a href="endometriosis/endo_qsp_model.svg"><img src="endometriosis/endo_qsp_model.png" width="190" alt="ENDO"></a> | 자궁내막 조직의 자궁외 이식·에스트로겐 의존 염증·통증. GnRH 길항제·프로게스틴.<br>[🗺️ 지도](endometriosis/endo_qsp_model.svg) · [⚙️ mrgsolve](endometriosis/endo_mrgsolve_model.R) · [📚 문헌](endometriosis/endo_references.md) · [📄 README](endometriosis/README.md) |
| 45 | 신경 | [**뇌전증**<br><sub>Epilepsy · EPI</sub>](epilepsy/) | <a href="epilepsy/epi_qsp_model.svg"><img src="epilepsy/epi_qsp_model.png" width="190" alt="EPI"></a> | 이온채널·흥분/억제 불균형에 의한 반복 발작. 항발작제(VPA/LEV/CBZ/LTG).<br>[🗺️ 지도](epilepsy/epi_qsp_model.svg) · [⚙️ mrgsolve](epilepsy/epi_mrgsolve_model.R) · [📚 문헌](epilepsy/epi_references.md) · [📄 README](epilepsy/README.md) |
| 46 | 심혈관 | [**본태성 고혈압**<br><sub>Essential Hypertension · EH</sub>](essential-hypertension/) | <a href="essential-hypertension/eh_qsp_model.svg"><img src="essential-hypertension/eh_qsp_model.png" width="190" alt="EH"></a> | RAAS·교감신경·나트륨/체액 항상성에 의한 혈압 조절. ACEi/ARB·CCB·이뇨제.<br>[🗺️ 지도](essential-hypertension/eh_qsp_model.svg) · [⚙️ mrgsolve](essential-hypertension/eh_mrgsolve_model.R) · [📚 문헌](essential-hypertension/eh_references.md) · [📄 README](essential-hypertension/README.md) |
| 47 | 혈액 | [**에반스 증후군**<br><sub>Evans Syndrome · ES</sub>](evans-syndrome/) | <a href="evans-syndrome/es_qsp_model.svg"><img src="evans-syndrome/es_qsp_model.png" width="190" alt="ES"></a> | AIHA와 ITP가 동반된 다계열 자가면역 혈구감소. 스테로이드·리툭시맙·TPO-RA.<br>[🗺️ 지도](evans-syndrome/es_qsp_model.svg) · [⚙️ mrgsolve](evans-syndrome/es_mrgsolve_model.R) · [📚 문헌](evans-syndrome/es_references.md) · [📄 README](evans-syndrome/README.md) |
| 48 | 신장·비뇨 | [**국소분절사구체경화증 (FSGS)**<br><sub>Focal Segmental Glomerulosclerosis · FSGS</sub>](fsgs/) | <a href="fsgs/fsgs_qsp_model.svg"><img src="fsgs/fsgs_qsp_model.png" width="190" alt="FSGS"></a> | 순환 투과인자·족세포 손상에 의한 신증후군 단백뇨. RAAS 억제·스테로이드·칼시뉴린억제제.<br>[🗺️ 지도](fsgs/fsgs_qsp_model.svg) · [⚙️ mrgsolve](fsgs/fsgs_mrgsolve_model.R) · [📚 문헌](fsgs/fsgs_references.md) · [📄 README](fsgs/README.md) |
| 49 | 소화기·간담도 | [**위식도 역류질환 (GERD)**<br><sub>Gastroesophageal Reflux Disease · GERD</sub>](gerd/) | <a href="gerd/gerd_qsp_model.svg"><img src="gerd/gerd_qsp_model.png" width="190" alt="GERD"></a> | 하부식도괄약근 이완·산 역류에 의한 식도 점막 손상. PPI·P-CAB 위산 억제.<br>[🗺️ 지도](gerd/gerd_qsp_model.svg) · [⚙️ mrgsolve](gerd/gerd_mrgsolve_model.R) · [📚 문헌](gerd/gerd_references.md) · [📄 README](gerd/README.md) |
| 50 | 혈관염 | [**거대세포 동맥염**<br><sub>Giant Cell Arteritis · GCA</sub>](giant-cell-arteritis/) | <a href="giant-cell-arteritis/gca_qsp_model.svg"><img src="giant-cell-arteritis/gca_qsp_model.png" width="190" alt="GCA"></a> | 대혈관의 IL-6/Th17 육아종성 동맥염과 실명 위험. 스테로이드·토실리주맙(항IL-6R).<br>[🗺️ 지도](giant-cell-arteritis/gca_qsp_model.svg) · [⚙️ mrgsolve](giant-cell-arteritis/gca_mrgsolve_model.R) · [📚 문헌](giant-cell-arteritis/gca_references.md) · [📄 README](giant-cell-arteritis/README.md) |
| 51 | 신장·비뇨 | [**굿파스처 증후군**<br><sub>Goodpasture Syndrome · GPS</sub>](goodpasture-syndrome/) | <a href="goodpasture-syndrome/gps_qsp_model.svg"><img src="goodpasture-syndrome/gps_qsp_model.png" width="190" alt="GPS"></a> | 항GBM 항체에 의한 급속진행성 사구체신염·폐포출혈. 혈장교환·CY·리툭시맙.<br>[🗺️ 지도](goodpasture-syndrome/gps_qsp_model.svg) · [⚙️ mrgsolve](goodpasture-syndrome/gps_mrgsolve_model.R) · [📚 문헌](goodpasture-syndrome/gps_references.md) · [📄 README](goodpasture-syndrome/README.md) |
| 52 | 내분비·대사 | [**통풍**<br><sub>Gout · GOUT</sub>](gout/) | <a href="gout/gout_qsp_model.svg"><img src="gout/gout_qsp_model.png" width="190" alt="GOUT"></a> | 요산 과포화·MSU 결정·NLRP3-IL-1β 염증성 관절염. 요산저하제(알로푸리놀/페북소스타트)·콜히친.<br>[🗺️ 지도](gout/gout_qsp_model.svg) · [⚙️ mrgsolve](gout/gout_mrgsolve_model.R) · [📚 문헌](gout/gout_references.md) · [📄 README](gout/README.md) |
| 53 | 혈관염 | [**육아종증 다발혈관염 (GPA)**<br><sub>Granulomatosis with Polyangiitis · GPA</sub>](granulomatosis-with-polyangiitis/) | <a href="granulomatosis-with-polyangiitis/gpa_qsp_model.svg"><img src="granulomatosis-with-polyangiitis/gpa_qsp_model.png" width="190" alt="GPA"></a> | PR3-ANCA 호중구 활성·육아종성 소혈관염. 리툭시맙·아바코판(C5aR 차단).<br>[🗺️ 지도](granulomatosis-with-polyangiitis/gpa_qsp_model.svg) · [⚙️ mrgsolve](granulomatosis-with-polyangiitis/gpa_mrgsolve_model.R) · [📚 문헌](granulomatosis-with-polyangiitis/gpa_references.md) · [📄 README](granulomatosis-with-polyangiitis/README.md) |
| 54 | 내분비·대사 | [**그레이브스병**<br><sub>Graves' Disease · GD</sub>](graves-disease/) | <a href="graves-disease/gd_qsp_model.svg"><img src="graves-disease/gd_qsp_model.png" width="190" alt="GD"></a> | TSH 수용체 자극항체에 의한 갑상선기능항진. 항갑상선제·방사성요오드·수술.<br>[🗺️ 지도](graves-disease/gd_qsp_model.svg) · [⚙️ mrgsolve](graves-disease/gd_mrgsolve_model.R) · [📚 문헌](graves-disease/gd_references.md) · [📄 README](graves-disease/README.md) |
| 55 | 신경 | [**길랭-바레 증후군**<br><sub>Guillain-Barré Syndrome · GBS</sub>](guillain-barre-syndrome/) | <a href="guillain-barre-syndrome/gbs_qsp_model.svg"><img src="guillain-barre-syndrome/gbs_qsp_model.png" width="190" alt="GBS"></a> | 분자모방 항강글리오시드 항체·보체 매개 급성 신경병증. IVIG·혈장교환.<br>[🗺️ 지도](guillain-barre-syndrome/gbs_qsp_model.svg) · [⚙️ mrgsolve](guillain-barre-syndrome/gbs_mrgsolve_model.R) · [📚 문헌](guillain-barre-syndrome/gbs_references.md) · [📄 README](guillain-barre-syndrome/README.md) |
| 56 | 내분비·대사 | [**하시모토 갑상선염**<br><sub>Hashimoto's Thyroiditis · HT</sub>](hashimoto-thyroiditis/) | <a href="hashimoto-thyroiditis/ht_qsp_model.svg"><img src="hashimoto-thyroiditis/ht_qsp_model.png" width="190" alt="HT"></a> | 항TPO/Tg 항체·T세포 매개 갑상선 파괴와 기능저하. 레보티록신 보충.<br>[🗺️ 지도](hashimoto-thyroiditis/ht_qsp_model.svg) · [⚙️ mrgsolve](hashimoto-thyroiditis/ht_mrgsolve_model.R) · [📚 문헌](hashimoto-thyroiditis/ht_references.md) · [📄 README](hashimoto-thyroiditis/README.md) |
| 57 | 심혈관 | [**심부전 (보존 박출률, HFpEF)**<br><sub>Heart Failure with Preserved EF · HFpEF</sub>](heart-failure-hfpef/) | <a href="heart-failure-hfpef/hfpef_qsp_model.svg"><img src="heart-failure-hfpef/hfpef_qsp_model.png" width="190" alt="HFpEF"></a> | 심근 경직·전신 염증·미세혈관 기능부전에 의한 확장기 부전. SGLT2i·MRA·이뇨제.<br>[🗺️ 지도](heart-failure-hfpef/hfpef_qsp_model.svg) · [⚙️ mrgsolve](heart-failure-hfpef/hfpef_mrgsolve_model.R) · [📚 문헌](heart-failure-hfpef/hfpef_references.md) · [📄 README](heart-failure-hfpef/README.md) |
| 58 | 심혈관 | [**심부전 (감소 박출률, HFrEF)**<br><sub>Heart Failure with Reduced EF · HFrEF</sub>](heart-failure-hfref/) | <a href="heart-failure-hfref/hfref_qsp_model.svg"><img src="heart-failure-hfref/hfref_qsp_model.png" width="190" alt="HFrEF"></a> | 신경호르몬(RAAS/SNS) 과활성과 심실 리모델링. ARNI·BB·MRA·SGLT2i 4대 요법.<br>[🗺️ 지도](heart-failure-hfref/hfref_qsp_model.svg) · [⚙️ mrgsolve](heart-failure-hfref/hfref_mrgsolve_model.R) · [📚 문헌](heart-failure-hfref/hfref_references.md) · [📄 README](heart-failure-hfref/README.md) |
| 59 | 감염 | [**HIV/AIDS**<br><sub>HIV/AIDS · HIV</sub>](hiv-aids/) | <a href="hiv-aids/hiv_qsp_model.svg"><img src="hiv-aids/hiv_qsp_model.png" width="190" alt="HIV"></a> | HIV의 CD4 T세포 감염·고갈(Perelson 바이러스 동역학). 항레트로바이러스 병합요법(ART).<br>[🗺️ 지도](hiv-aids/hiv_qsp_model.svg) · [⚙️ mrgsolve](hiv-aids/hiv_mrgsolve_model.R) · [📚 문헌](hiv-aids/hiv_references.md) · [📄 README](hiv-aids/README.md) |
| 60 | 심혈관 | [**비후성 심근병증**<br><sub>Hypertrophic Cardiomyopathy · HCM</sub>](hypertrophic-cardiomyopathy/) | <a href="hypertrophic-cardiomyopathy/hcm_qsp_model.svg"><img src="hypertrophic-cardiomyopathy/hcm_qsp_model.png" width="190" alt="HCM"></a> | 근절 변이에 의한 과수축·좌심실유출로 폐쇄. 마바캄텐(마이오신 억제)·베타차단제.<br>[🗺️ 지도](hypertrophic-cardiomyopathy/hcm_qsp_model.svg) · [⚙️ mrgsolve](hypertrophic-cardiomyopathy/hcm_mrgsolve_model.R) · [📚 문헌](hypertrophic-cardiomyopathy/hcm_references.md) · [📄 README](hypertrophic-cardiomyopathy/README.md) |
| 61 | 호흡기 | [**특발성 폐섬유화증 (IPF)**<br><sub>Idiopathic Pulmonary Fibrosis · IPF</sub>](idiopathic-pulmonary-fibrosis/) | <a href="idiopathic-pulmonary-fibrosis/ipf_qsp_model.svg"><img src="idiopathic-pulmonary-fibrosis/ipf_qsp_model.png" width="190" alt="IPF"></a> | 상피손상-섬유아세포 활성·ECM 침착에 의한 진행성 섬유화. 닌테다닙·피르페니돈.<br>[🗺️ 지도](idiopathic-pulmonary-fibrosis/ipf_qsp_model.svg) · [⚙️ mrgsolve](idiopathic-pulmonary-fibrosis/ipf_mrgsolve_model.R) · [📚 문헌](idiopathic-pulmonary-fibrosis/ipf_references.md) · [📄 README](idiopathic-pulmonary-fibrosis/README.md) |
| 62 | 신장·비뇨 | [**IgA 신병증**<br><sub>IgA Nephropathy · IgAN</sub>](iga-nephropathy/) | <a href="iga-nephropathy/igan_qsp_model.svg"><img src="iga-nephropathy/igan_qsp_model.png" width="190" alt="IgAN"></a> | Gd-IgA1 면역복합체 메산지움 침착·보체 활성. 스테로이드·스파르센탄·보체억제제(이프타코판).<br>[🗺️ 지도](iga-nephropathy/igan_qsp_model.svg) · [⚙️ mrgsolve](iga-nephropathy/igan_mrgsolve_model.R) · [📚 문헌](iga-nephropathy/igan_references.md) · [📄 README](iga-nephropathy/README.md) |
| 63 | 혈관염 | [**IgA 혈관염 (HSP)**<br><sub>IgA Vasculitis · IgAV</sub>](iga-vasculitis/) | <a href="iga-vasculitis/igav_qsp_model.svg"><img src="iga-vasculitis/igav_qsp_model.png" width="190" alt="IgAV"></a> | Gd-IgA1 면역복합체 소혈관 침착에 의한 촉지성 자반·신염. 스테로이드·면역억제.<br>[🗺️ 지도](iga-vasculitis/igav_qsp_model.svg) · [⚙️ mrgsolve](iga-vasculitis/igav_mrgsolve_model.R) · [📚 문헌](iga-vasculitis/igav_references.md) · [📄 README](iga-vasculitis/README.md) |
| 64 | 혈액 | [**면역혈소판감소자반증 (ITP)**<br><sub>Immune Thrombocytopenic Purpura · ITP</sub>](immune-thrombocytopenic-purpura/) | <a href="immune-thrombocytopenic-purpura/itp_qsp_model.svg"><img src="immune-thrombocytopenic-purpura/itp_qsp_model.png" width="190" alt="ITP"></a> | 항혈소판 항체 매개 파괴와 생성 저하. 스테로이드·TPO 수용체작용제·리툭시맙.<br>[🗺️ 지도](immune-thrombocytopenic-purpura/itp_qsp_model.svg) · [⚙️ mrgsolve](immune-thrombocytopenic-purpura/itp_mrgsolve_model.R) · [📚 문헌](immune-thrombocytopenic-purpura/itp_references.md) · [📄 README](immune-thrombocytopenic-purpura/README.md) |
| 65 | 소화기·간담도 | [**과민성 장증후군 (IBS)**<br><sub>Irritable Bowel Syndrome · IBS</sub>](irritable-bowel-syndrome/) | <a href="irritable-bowel-syndrome/ibs_qsp_model.svg"><img src="irritable-bowel-syndrome/ibs_qsp_model.png" width="190" alt="IBS"></a> | 뇌-장축 이상·내장과민·미생물 변화. 식이·신경조절·장특이 약물(리나클로타이드 등).<br>[🗺️ 지도](irritable-bowel-syndrome/ibs_qsp_model.svg) · [⚙️ mrgsolve](irritable-bowel-syndrome/ibs_mrgsolve_model.R) · [📚 문헌](irritable-bowel-syndrome/ibs_references.md) · [📄 README](irritable-bowel-syndrome/README.md) |
| 66 | 소화기·간담도 | [**간경변증**<br><sub>Liver Cirrhosis · LC</sub>](liver-cirrhosis/) | <a href="liver-cirrhosis/lc_qsp_model.svg"><img src="liver-cirrhosis/lc_qsp_model.png" width="190" alt="LC"></a> | 만성 손상 → 성상세포 섬유화·문맥압항진·합병증. 원인치료 및 합병증 관리.<br>[🗺️ 지도](liver-cirrhosis/lc_qsp_model.svg) · [⚙️ mrgsolve](liver-cirrhosis/lc_mrgsolve_model.R) · [📚 문헌](liver-cirrhosis/lc_references.md) · [📄 README](liver-cirrhosis/README.md) |
| 67 | 내분비·대사 | [**림프구성 뇌하수체염**<br><sub>Lymphocytic Hypophysitis · LHY</sub>](lymphocytic-hypophysitis/) | <a href="lymphocytic-hypophysitis/lhyp_qsp_model.svg"><img src="lymphocytic-hypophysitis/lhyp_qsp_model.png" width="190" alt="LHY"></a> | 뇌하수체 자가면역 침윤·뇌하수체 기능저하(면역관문억제제 연관 포함). 스테로이드·호르몬 보충.<br>[🗺️ 지도](lymphocytic-hypophysitis/lhyp_qsp_model.svg) · [⚙️ mrgsolve](lymphocytic-hypophysitis/lhyp_mrgsolve_model.R) · [📚 문헌](lymphocytic-hypophysitis/lhyp_references.md) · [📄 README](lymphocytic-hypophysitis/README.md) |
| 68 | 정신·신경 | [**주요우울장애 (MDD)**<br><sub>Major Depressive Disorder · MDD</sub>](major-depressive-disorder/) | <a href="major-depressive-disorder/mdd_qsp_model.svg"><img src="major-depressive-disorder/mdd_qsp_model.png" width="190" alt="MDD"></a> | 모노아민·HPA축·신경가소성 이상. SSRI/SNRI·케타민·신경자극.<br>[🗺️ 지도](major-depressive-disorder/mdd_qsp_model.svg) · [⚙️ mrgsolve](major-depressive-disorder/mdd_mrgsolve_model.R) · [📚 문헌](major-depressive-disorder/mdd_references.md) · [📄 README](major-depressive-disorder/README.md) |
| 69 | 신장·비뇨 | [**막성 신병증**<br><sub>Membranous Nephropathy · MN</sub>](membranous-nephropathy/) | <a href="membranous-nephropathy/mn_qsp_model.svg"><img src="membranous-nephropathy/mn_qsp_model.png" width="190" alt="MN"></a> | 항PLA2R 항체에 의한 상피하 면역침착·신증후군. 리툭시맙 B세포 고갈.<br>[🗺️ 지도](membranous-nephropathy/mn_qsp_model.svg) · [⚙️ mrgsolve](membranous-nephropathy/mn_mrgsolve_model.R) · [📚 문헌](membranous-nephropathy/mn_references.md) · [📄 README](membranous-nephropathy/README.md) |
| 70 | 내분비·대사 | [**대사 증후군**<br><sub>Metabolic Syndrome · MS</sub>](metabolic-syndrome/) | <a href="metabolic-syndrome/ms_qsp_model.svg"><img src="metabolic-syndrome/ms_qsp_model.png" width="190" alt="MS"></a> | 인슐린 저항성·내장지방·이상지질·고혈압 군집. 생활습관 및 대사 표적 약물.<br>[🗺️ 지도](metabolic-syndrome/ms_qsp_model.svg) · [⚙️ mrgsolve](metabolic-syndrome/ms_mrgsolve_model.R) · [📚 문헌](metabolic-syndrome/ms_references.md) · [📄 README](metabolic-syndrome/README.md) |
| 71 | 혈관염 | [**현미경적 다발혈관염 (MPA)**<br><sub>Microscopic Polyangiitis · MPA</sub>](microscopic-polyangiitis/) | <a href="microscopic-polyangiitis/mpa_qsp_model.svg"><img src="microscopic-polyangiitis/mpa_qsp_model.png" width="190" alt="MPA"></a> | MPO-ANCA 매개 소혈관염·급속진행성 사구체신염. 리툭시맙·아바코판.<br>[🗺️ 지도](microscopic-polyangiitis/mpa_qsp_model.svg) · [⚙️ mrgsolve](microscopic-polyangiitis/mpa_mrgsolve_model.R) · [📚 문헌](microscopic-polyangiitis/mpa_references.md) · [📄 README](microscopic-polyangiitis/README.md) |
| 72 | 신경 | [**편두통**<br><sub>Migraine · MGR</sub>](migraine/) | <a href="migraine/mgr_qsp_model.svg"><img src="migraine/mgr_qsp_model.png" width="190" alt="MGR"></a> | 삼차혈관계·CGRP·피질확산성억제(CSD). 트립탄·항CGRP 항체·게판트.<br>[🗺️ 지도](migraine/mgr_qsp_model.svg) · [⚙️ mrgsolve](migraine/mgr_mrgsolve_model.R) · [📚 문헌](migraine/mgr_references.md) · [📄 README](migraine/README.md) |
| 73 | 신장·비뇨 | [**미세변화 신증후군**<br><sub>Minimal Change Disease · MCD</sub>](minimal-change-disease/) | <a href="minimal-change-disease/mcd_qsp_model.svg"><img src="minimal-change-disease/mcd_qsp_model.png" width="190" alt="MCD"></a> | T세포 매개 순환인자에 의한 족세포 손상·대량 단백뇨. 스테로이드·칼시뉴린억제제.<br>[🗺️ 지도](minimal-change-disease/mcd_qsp_model.svg) · [⚙️ mrgsolve](minimal-change-disease/mcd_mrgsolve_model.R) · [📚 문헌](minimal-change-disease/mcd_references.md) · [📄 README](minimal-change-disease/README.md) |
| 74 | 자가면역·류마티스 | [**혼합결합조직병 (MCTD)**<br><sub>Mixed Connective Tissue Disease · MCTD</sub>](mixed-connective-tissue-disease/) | <a href="mixed-connective-tissue-disease/mctd_qsp_model.svg"><img src="mixed-connective-tissue-disease/mctd_qsp_model.png" width="190" alt="MCTD"></a> | 항U1-RNP 항체 양성, SLE/SSc/PM 중복 양상. 스테로이드·면역억제.<br>[🗺️ 지도](mixed-connective-tissue-disease/mctd_qsp_model.svg) · [⚙️ mrgsolve](mixed-connective-tissue-disease/mctd_mrgsolve_model.R) · [📚 문헌](mixed-connective-tissue-disease/mctd_references.md) · [📄 README](mixed-connective-tissue-disease/README.md) |
| 75 | 종양 | [**다발골수종**<br><sub>Multiple Myeloma · MM</sub>](multiple-myeloma/) | <a href="multiple-myeloma/mm_qsp_model.svg"><img src="multiple-myeloma/mm_qsp_model.png" width="190" alt="MM"></a> | 골수 형질세포 클론 증식·골 파괴. 프로테아좀 억제제·IMiD·항CD38.<br>[🗺️ 지도](multiple-myeloma/mm_qsp_model.svg) · [⚙️ mrgsolve](multiple-myeloma/mm_mrgsolve_model.R) · [📚 문헌](multiple-myeloma/mm_references.md) · [📄 README](multiple-myeloma/README.md) |
| 76 | 신경 | [**다발성 경화증**<br><sub>Multiple Sclerosis · MS</sub>](multiple-sclerosis/) | <a href="multiple-sclerosis/ms_qsp.svg"><img src="multiple-sclerosis/ms_qsp.png" width="190" alt="MS"></a> | 자가반응 T/B세포에 의한 중추신경 탈수초·축삭 손상. 질병조절치료(항CD20·S1P 조절제).<br>[🗺️ 지도](multiple-sclerosis/ms_qsp.svg) · [⚙️ mrgsolve](multiple-sclerosis/ms_mrgsolve_model.R) · [📚 문헌](multiple-sclerosis/ms_references.md) · [📄 README](multiple-sclerosis/README.md) |
| 77 | 신경 | [**중증 근무력증**<br><sub>Myasthenia Gravis · MG</sub>](myasthenia-gravis/) | <a href="myasthenia-gravis/mg_qsp_model.svg"><img src="myasthenia-gravis/mg_qsp_model.png" width="190" alt="MG"></a> | 항AChR 항체·보체 매개 신경근접합 차단. 콜린에스터분해효소 억제·FcRn/보체 억제제.<br>[🗺️ 지도](myasthenia-gravis/mg_qsp_model.svg) · [⚙️ mrgsolve](myasthenia-gravis/mg_mrgsolve_model.R) · [📚 문헌](myasthenia-gravis/mg_references.md) · [📄 README](myasthenia-gravis/README.md) |
| 78 | 소화기·간담도 | [**비알코올 지방간/지방간염 (NAFLD/NASH)**<br><sub>NAFLD/NASH · NAFLD</sub>](nafld-nash/) | <a href="nafld-nash/nafld_qsp_model.svg"><img src="nafld-nash/nafld_qsp_model.png" width="190" alt="NAFLD"></a> | 지방독성·염증·섬유화 진행(지방간→지방간염). 레스메티롬(THR-β)·GLP-1 작용제.<br>[🗺️ 지도](nafld-nash/nafld_qsp_model.svg) · [⚙️ mrgsolve](nafld-nash/nafld_mrgsolve_model.R) · [📚 문헌](nafld-nash/nafld_references.md) · [📄 README](nafld-nash/README.md) |
| 79 | 신경 | [**시신경척수염 (NMOSD)**<br><sub>Neuromyelitis Optica · NMO</sub>](neuromyelitis-optica/) | <a href="neuromyelitis-optica/nmo_qsp_model.svg"><img src="neuromyelitis-optica/nmo_qsp_model.png" width="190" alt="NMO"></a> | 항AQP4 항체·보체 매개 성상세포 손상. 에쿨리주맙·항IL-6R(사트랄리주맙)·항CD19(이네빌리주맙).<br>[🗺️ 지도](neuromyelitis-optica/nmo_qsp_model.svg) · [⚙️ mrgsolve](neuromyelitis-optica/nmo_mrgsolve_model.R) · [📚 문헌](neuromyelitis-optica/nmo_references.md) · [📄 README](neuromyelitis-optica/README.md) |
| 80 | 내분비·대사 | [**비만**<br><sub>Obesity · OB</sub>](obesity/) | <a href="obesity/ob_qsp_model.svg"><img src="obesity/ob_qsp_model.png" width="190" alt="OB"></a> | 에너지 항상성·식욕조절(렙틴/GLP-1/멜라노코르틴) 조절이상. GLP-1/GIP 작용제 체중감량.<br>[🗺️ 지도](obesity/ob_qsp_model.svg) · [⚙️ mrgsolve](obesity/ob_mrgsolve_model.R) · [📚 문헌](obesity/ob_references.md) · [📄 README](obesity/README.md) |
| 81 | 호흡기 | [**폐쇄성 수면 무호흡 (OSA)**<br><sub>Obstructive Sleep Apnea · OSA</sub>](obstructive-sleep-apnea/) | <a href="obstructive-sleep-apnea/osa_qsp_model.svg"><img src="obstructive-sleep-apnea/osa_qsp_model.png" width="190" alt="OSA"></a> | 상기도 허탈·간헐적 저산소·교감신경 활성. CPAP·체중감량·약물 보조.<br>[🗺️ 지도](obstructive-sleep-apnea/osa_qsp_model.svg) · [⚙️ mrgsolve](obstructive-sleep-apnea/osa_mrgsolve_model.R) · [📚 문헌](obstructive-sleep-apnea/osa_references.md) · [📄 README](obstructive-sleep-apnea/README.md) |
| 82 | 내분비·대사 | [**골다공증**<br><sub>Osteoporosis · OP</sub>](osteoporosis/) | <a href="osteoporosis/op_qsp_model.svg"><img src="osteoporosis/op_qsp_model.png" width="190" alt="OP"></a> | RANKL-OPG 불균형으로 골흡수↑/형성↓. 비스포스포네이트·데노수맙·로모소주맙.<br>[🗺️ 지도](osteoporosis/op_qsp_model.svg) · [⚙️ mrgsolve](osteoporosis/op_mrgsolve_model.R) · [📚 문헌](osteoporosis/op_references.md) · [📄 README](osteoporosis/README.md) |
| 83 | 신장·비뇨 | [**과민성 방광 (OAB)**<br><sub>Overactive Bladder · OAB</sub>](overactive-bladder/) | <a href="overactive-bladder/oab_qsp_model.svg"><img src="overactive-bladder/oab_qsp_model.png" width="190" alt="OAB"></a> | 배뇨근 과활동(무스카린/β3 수용체 불균형). 항무스카린제·미라베그론(β3 작용제).<br>[🗺️ 지도](overactive-bladder/oab_qsp_model.svg) · [⚙️ mrgsolve](overactive-bladder/oab_mrgsolve_model.R) · [📚 문헌](overactive-bladder/oab_references.md) · [📄 README](overactive-bladder/README.md) |
| 84 | 내분비·대사 | [**파젯병 (골)**<br><sub>Paget's Disease of Bone · PBD</sub>](pagets-disease/) | <a href="pagets-disease/pbd_qsp_model.svg"><img src="pagets-disease/pbd_qsp_model.png" width="190" alt="PBD"></a> | SQSTM1/RANKL 과활성에 의한 비정상 골개조. 졸레드론산 등 비스포스포네이트.<br>[🗺️ 지도](pagets-disease/pbd_qsp_model.svg) · [⚙️ mrgsolve](pagets-disease/pbd_mrgsolve_model.R) · [📚 문헌](pagets-disease/pbd_references.md) · [📄 README](pagets-disease/README.md) |
| 85 | 신경 | [**파킨슨병**<br><sub>Parkinson's Disease · PD</sub>](parkinsons-disease/) | <a href="parkinsons-disease/pd_qsp_model.svg"><img src="parkinsons-disease/pd_qsp_model.png" width="190" alt="PD"></a> | 흑질 도파민 신경세포 변성·α-시누클레인 응집. 레보도파·도파민작용제·MAO-B 억제제.<br>[🗺️ 지도](parkinsons-disease/pd_qsp_model.svg) · [⚙️ mrgsolve](parkinsons-disease/pd_mrgsolve_model.R) · [📚 문헌](parkinsons-disease/pd_references.md) · [📄 README](parkinsons-disease/README.md) |
| 86 | 피부 | [**심상성 천포창**<br><sub>Pemphigus Vulgaris · PV</sub>](pemphigus-vulgaris/) | <a href="pemphigus-vulgaris/pv_qsp_model.svg"><img src="pemphigus-vulgaris/pv_qsp_model.png" width="190" alt="PV"></a> | 항데스모글레인 항체에 의한 표피내 수포(천포창). 리툭시맙·스테로이드.<br>[🗺️ 지도](pemphigus-vulgaris/pv_qsp_model.svg) · [⚙️ mrgsolve](pemphigus-vulgaris/pv_mrgsolve_model.R) · [📚 문헌](pemphigus-vulgaris/pv_references.md) · [📄 README](pemphigus-vulgaris/README.md) |
| 87 | 소화기·간담도 | [**소화성 궤양**<br><sub>Peptic Ulcer Disease · PUD</sub>](peptic-ulcer/) | <a href="peptic-ulcer/pud_qsp_model.svg"><img src="peptic-ulcer/pud_qsp_model.png" width="190" alt="PUD"></a> | H. pylori·NSAID에 의한 점막 방어-공격 불균형. PPI·제균요법 점막 치유.<br>[🗺️ 지도](peptic-ulcer/pud_qsp_model.svg) · [⚙️ mrgsolve](peptic-ulcer/pud_mrgsolve_model.R) · [📚 문헌](peptic-ulcer/pud_references.md) · [📄 README](peptic-ulcer/README.md) |
| 88 | 심혈관 | [**말초동맥질환 (PAD)**<br><sub>Peripheral Arterial Disease · PAD</sub>](peripheral-arterial-disease/) | <a href="peripheral-arterial-disease/pad_qsp_model.svg"><img src="peripheral-arterial-disease/pad_qsp_model.png" width="190" alt="PAD"></a> | 죽상경화에 의한 하지 허혈·파행. 항혈소판·스타틴·실로스타졸.<br>[🗺️ 지도](peripheral-arterial-disease/pad_qsp_model.svg) · [⚙️ mrgsolve](peripheral-arterial-disease/pad_mrgsolve_model.R) · [📚 문헌](peripheral-arterial-disease/pad_references.md) · [📄 README](peripheral-arterial-disease/README.md) |
| 89 | 혈액 | [**악성 빈혈**<br><sub>Pernicious Anemia · PNA</sub>](pernicious-anemia/) | <a href="pernicious-anemia/pna_qsp_model.svg"><img src="pernicious-anemia/pna_qsp_model.png" width="190" alt="PNA"></a> | 항내인자/벽세포 항체에 의한 비타민 B12 흡수장애. B12 보충요법.<br>[🗺️ 지도](pernicious-anemia/pna_qsp_model.svg) · [⚙️ mrgsolve](pernicious-anemia/pna_mrgsolve_model.R) · [📚 문헌](pernicious-anemia/pna_references.md) · [📄 README](pernicious-anemia/README.md) |
| 90 | 호흡기 | [**진폐증**<br><sub>Pneumoconiosis · PNM</sub>](pneumoconiosis/) | <a href="pneumoconiosis/pnm_qsp_model.svg"><img src="pneumoconiosis/pnm_qsp_model.png" width="190" alt="PNM"></a> | 분진(실리카/석탄/석면) 대식세포 활성·진행성 폐섬유화. 노출차단·대증치료.<br>[🗺️ 지도](pneumoconiosis/pnm_qsp_model.svg) · [⚙️ mrgsolve](pneumoconiosis/pnm_mrgsolve_model.R) · [📚 문헌](pneumoconiosis/pnm_references.md) · [📄 README](pneumoconiosis/README.md) |
| 91 | 혈관염 | [**결절성 다발동맥염 (PAN)**<br><sub>Polyarteritis Nodosa · PAN</sub>](polyarteritis-nodosa/) | <a href="polyarteritis-nodosa/pan_qsp_model.svg"><img src="polyarteritis-nodosa/pan_qsp_model.png" width="190" alt="PAN"></a> | 중형 동맥의 괴사성 염증(HBV 연관 포함). 스테로이드·CY·항바이러스.<br>[🗺️ 지도](polyarteritis-nodosa/pan_qsp_model.svg) · [⚙️ mrgsolve](polyarteritis-nodosa/pan_mrgsolve_model.R) · [📚 문헌](polyarteritis-nodosa/pan_references.md) · [📄 README](polyarteritis-nodosa/README.md) |
| 92 | 내분비·대사 | [**다낭성 난소 증후군 (PCOS)**<br><sub>Polycystic Ovary Syndrome · PCOS</sub>](polycystic-ovary-syndrome/) | <a href="polycystic-ovary-syndrome/pcos_qsp_model.svg"><img src="polycystic-ovary-syndrome/pcos_qsp_model.png" width="190" alt="PCOS"></a> | 고안드로겐·인슐린저항·LH/FSH 이상. 메트포르민·항안드로겐·배란유도.<br>[🗺️ 지도](polycystic-ovary-syndrome/pcos_qsp_model.svg) · [⚙️ mrgsolve](polycystic-ovary-syndrome/pcos_mrgsolve_model.R) · [📚 문헌](polycystic-ovary-syndrome/pcos_references.md) · [📄 README](polycystic-ovary-syndrome/README.md) |
| 93 | 자가면역·류마티스 | [**다발성 근염**<br><sub>Polymyositis · PM</sub>](polymyositis/) | <a href="polymyositis/pm_qsp_model.svg"><img src="polymyositis/pm_qsp_model.png" width="190" alt="PM"></a> | CD8 T세포 매개 근섬유 침습·근력저하. 스테로이드·면역억제제.<br>[🗺️ 지도](polymyositis/pm_qsp_model.svg) · [⚙️ mrgsolve](polymyositis/pm_mrgsolve_model.R) · [📚 문헌](polymyositis/pm_references.md) · [📄 README](polymyositis/README.md) |
| 94 | 소화기·간담도 | [**원발성 담즙성 담관염 (PBC)**<br><sub>Primary Biliary Cholangitis · PBC</sub>](primary-biliary-cholangitis/) | <a href="primary-biliary-cholangitis/pbc_qsp_model.svg"><img src="primary-biliary-cholangitis/pbc_qsp_model.png" width="190" alt="PBC"></a> | 항미토콘드리아항체에 의한 담관 파괴·담즙정체. UDCA·오베티콜산.<br>[🗺️ 지도](primary-biliary-cholangitis/pbc_qsp_model.svg) · [⚙️ mrgsolve](primary-biliary-cholangitis/pbc_mrgsolve_model.R) · [📚 문헌](primary-biliary-cholangitis/pbc_references.md) · [📄 README](primary-biliary-cholangitis/README.md) |
| 95 | 내분비·대사 | [**원발성 부갑상선 기능 항진증 (PHPT)**<br><sub>Primary Hyperparathyroidism · PHPT</sub>](primary-hyperparathyroidism/) | <a href="primary-hyperparathyroidism/phpt_qsp_model.svg"><img src="primary-hyperparathyroidism/phpt_qsp_model.png" width="190" alt="PHPT"></a> | 부갑상선 선종의 PTH 자율과다·고칼슘혈증·골소실. 시나칼셋·부갑상선절제술.<br>[🗺️ 지도](primary-hyperparathyroidism/phpt_qsp_model.svg) · [⚙️ mrgsolve](primary-hyperparathyroidism/phpt_mrgsolve_model.R) · [📚 문헌](primary-hyperparathyroidism/phpt_references.md) · [📄 README](primary-hyperparathyroidism/README.md) |
| 96 | 소화기·간담도 | [**원발성 경화성 담관염 (PSC)**<br><sub>Primary Sclerosing Cholangitis · PSC</sub>](primary-sclerosing-cholangitis/) | <a href="primary-sclerosing-cholangitis/psc_qsp_model.svg"><img src="primary-sclerosing-cholangitis/psc_qsp_model.png" width="190" alt="PSC"></a> | 담관 섬유화·다발 협착(IBD 연관)·담관암 위험. 대증치료·간이식.<br>[🗺️ 지도](primary-sclerosing-cholangitis/psc_qsp_model.svg) · [⚙️ mrgsolve](primary-sclerosing-cholangitis/psc_mrgsolve_model.R) · [📚 문헌](primary-sclerosing-cholangitis/psc_references.md) · [📄 README](primary-sclerosing-cholangitis/README.md) |
| 97 | 내분비·대사 | [**가성통풍 (CPPD)**<br><sub>Pseudogout (CPPD) · CPPD</sub>](pseudogout/) | <a href="pseudogout/cppd_qsp_model.svg"><img src="pseudogout/cppd_qsp_model.png" width="190" alt="CPPD"></a> | 칼슘피로인산 결정 침착·NLRP3 염증성 관절염. 콜히친·NSAID·스테로이드.<br>[🗺️ 지도](pseudogout/cppd_qsp_model.svg) · [⚙️ mrgsolve](pseudogout/cppd_mrgsolve_model.R) · [📚 문헌](pseudogout/cppd_references.md) · [📄 README](pseudogout/README.md) |
| 98 | 피부 | [**건선**<br><sub>Psoriasis · PSO</sub>](psoriasis/) | <a href="psoriasis/pso_qsp_model.svg"><img src="psoriasis/pso_qsp_model.png" width="190" alt="PSO"></a> | IL-23/IL-17 축에 의한 각질세포 과증식. 항IL-17/IL-23 생물학제제·항TNF.<br>[🗺️ 지도](psoriasis/pso_qsp_model.svg) · [⚙️ mrgsolve](psoriasis/pso_mrgsolve_model.R) · [📚 문헌](psoriasis/pso_references.md) · [📄 README](psoriasis/README.md) |
| 99 | 자가면역·류마티스 | [**건선성 관절염**<br><sub>Psoriatic Arthritis · PsA</sub>](psoriatic-arthritis/) | <a href="psoriatic-arthritis/psa_qsp_model.svg"><img src="psoriatic-arthritis/psa_qsp_model.png" width="190" alt="PsA"></a> | TNF·IL-23/IL-17 매개 부착부염·관절염·피부 건선. 생물학제제·JAK 억제제.<br>[🗺️ 지도](psoriatic-arthritis/psa_qsp_model.svg) · [⚙️ mrgsolve](psoriatic-arthritis/psa_mrgsolve_model.R) · [📚 문헌](psoriatic-arthritis/psa_references.md) · [📄 README](psoriatic-arthritis/README.md) |
| 100 | 심혈관 | [**폐동맥 고혈압 (PAH)**<br><sub>Pulmonary Arterial Hypertension · PAH</sub>](pulmonary-arterial-hypertension/) | <a href="pulmonary-arterial-hypertension/pah_qsp_model.svg"><img src="pulmonary-arterial-hypertension/pah_qsp_model.png" width="190" alt="PAH"></a> | 폐혈관 리모델링(엔도텔린·NO·프로스타사이클린 경로). ERA·PDE5i·프로스타사이클린.<br>[🗺️ 지도](pulmonary-arterial-hypertension/pah_qsp_model.svg) · [⚙️ mrgsolve](pulmonary-arterial-hypertension/pah_mrgsolve_model.R) · [📚 문헌](pulmonary-arterial-hypertension/pah_references.md) · [📄 README](pulmonary-arterial-hypertension/README.md) |
| 101 | 자가면역·류마티스 | [**반응성 관절염**<br><sub>Reactive Arthritis · ReA</sub>](reactive-arthritis/) | <a href="reactive-arthritis/rea_qsp_model.svg"><img src="reactive-arthritis/rea_qsp_model.png" width="190" alt="ReA"></a> | 감염 후 HLA-B27 연관 무균성 관절염. NSAID·설파살라진·생물학제제.<br>[🗺️ 지도](reactive-arthritis/rea_qsp_model.svg) · [⚙️ mrgsolve](reactive-arthritis/rea_mrgsolve_model.R) · [📚 문헌](reactive-arthritis/rea_references.md) · [📄 README](reactive-arthritis/README.md) |
| 102 | 자가면역·류마티스 | [**재발성 다발연골염**<br><sub>Relapsing Polychondritis · RPC</sub>](relapsing-polychondritis/) | <a href="relapsing-polychondritis/rpc_qsp_model.svg"><img src="relapsing-polychondritis/rpc_qsp_model.png" width="190" alt="RPC"></a> | 연골(II형 콜라겐) 자가면역 염증(귀·코·기도). 스테로이드·면역억제·생물학제제.<br>[🗺️ 지도](relapsing-polychondritis/rpc_qsp_model.svg) · [⚙️ mrgsolve](relapsing-polychondritis/rpc_mrgsolve_model.R) · [📚 문헌](relapsing-polychondritis/rpc_references.md) · [📄 README](relapsing-polychondritis/README.md) |
| 103 | 자가면역·류마티스 | [**류마티스 관절염**<br><sub>Rheumatoid Arthritis · RA</sub>](rheumatoid-arthritis/) | <a href="rheumatoid-arthritis/ra_qsp_model.svg"><img src="rheumatoid-arthritis/ra_qsp_model.png" width="190" alt="RA"></a> | TNF/IL-6 매개 활막 판누스·골미란. 항TNF·항IL-6·JAK 억제제·CTLA-4-Ig.<br>[🗺️ 지도](rheumatoid-arthritis/ra_qsp_model.svg) · [⚙️ mrgsolve](rheumatoid-arthritis/ra_mrgsolve_model.R) · [📚 문헌](rheumatoid-arthritis/ra_references.md) · [📄 README](rheumatoid-arthritis/README.md) |
| 104 | 호흡기 | [**사르코이드증**<br><sub>Sarcoidosis · SARC</sub>](sarcoidosis/) | <a href="sarcoidosis/sarc_qsp_model.svg"><img src="sarcoidosis/sarc_qsp_model.png" width="190" alt="SARC"></a> | Th1/IFN-γ 매개 비건락성 육아종(폐 우세 다장기). 스테로이드·메토트렉세이트·항TNF.<br>[🗺️ 지도](sarcoidosis/sarc_qsp_model.svg) · [⚙️ mrgsolve](sarcoidosis/sarc_mrgsolve_model.R) · [📚 문헌](sarcoidosis/sarc_references.md) · [📄 README](sarcoidosis/README.md) |
| 105 | 정신·신경 | [**조현병**<br><sub>Schizophrenia · SCH</sub>](schizophrenia/) | <a href="schizophrenia/sch_qsp_model.svg"><img src="schizophrenia/sch_qsp_model.png" width="190" alt="SCH"></a> | 도파민/글루탐산/GABA/세로토닌 신경전달 이상. 항정신병약(D2 길항/부분작용).<br>[🗺️ 지도](schizophrenia/sch_qsp_model.svg) · [⚙️ mrgsolve](schizophrenia/sch_mrgsolve_model.R) · [📚 문헌](schizophrenia/sch_references.md) · [📄 README](schizophrenia/README.md) |
| 106 | 혈액 | [**겸상적혈구병**<br><sub>Sickle Cell Disease · SCD</sub>](sickle-cell-disease/) | <a href="sickle-cell-disease/scd_qsp_model.svg"><img src="sickle-cell-disease/scd_qsp_model.png" width="190" alt="SCD"></a> | HbS 중합·적혈구 겸상화·혈관폐색. 하이드록시유레아·복셀로토르·크리잔리주맙·L-글루타민.<br>[🗺️ 지도](sickle-cell-disease/scd_qsp_model.svg) · [⚙️ mrgsolve](sickle-cell-disease/scd_mrgsolve_model.R) · [📚 문헌](sickle-cell-disease/scd_references.md) · [📄 README](sickle-cell-disease/README.md) |
| 107 | 자가면역·류마티스 | [**쇼그렌 증후군**<br><sub>Sjögren's Syndrome · SS</sub>](sjogrens-syndrome/) | <a href="sjogrens-syndrome/ss_qsp_model.svg"><img src="sjogrens-syndrome/ss_qsp_model.png" width="190" alt="SS"></a> | 외분비선 림프구 침윤·I형 IFN에 의한 건조증. 대증치료·전신 면역조절.<br>[🗺️ 지도](sjogrens-syndrome/ss_qsp_model.svg) · [⚙️ mrgsolve](sjogrens-syndrome/ss_mrgsolve_model.R) · [📚 문헌](sjogrens-syndrome/ss_references.md) · [📄 README](sjogrens-syndrome/README.md) |
| 108 | 심혈관 | [**안정형 협심증**<br><sub>Stable Angina · SA</sub>](stable-angina/) | <a href="stable-angina/sa_qsp_model.svg"><img src="stable-angina/sa_qsp_model.png" width="190" alt="SA"></a> | 죽상경화 관상동맥의 산소 수급 불균형. 항허혈제·항혈소판·스타틴.<br>[🗺️ 지도](stable-angina/sa_qsp_model.svg) · [⚙️ mrgsolve](stable-angina/sa_mrgsolve_model.R) · [📚 문헌](stable-angina/sa_references.md) · [📄 README](stable-angina/README.md) |
| 109 | 자가면역·류마티스 | [**전신 홍반 루푸스 (SLE)**<br><sub>Systemic Lupus Erythematosus · SLE</sub>](systemic-lupus-erythematosus/) | <a href="systemic-lupus-erythematosus/sle_qsp.svg"><img src="systemic-lupus-erythematosus/sle_qsp.png" width="190" alt="SLE"></a> | I형 IFN·항dsDNA 면역복합체에 의한 다장기 손상. 항말라리아·벨리무맙·아니프롤루맙.<br>[🗺️ 지도](systemic-lupus-erythematosus/sle_qsp.svg) · [⚙️ mrgsolve](systemic-lupus-erythematosus/sle_model.R) · [📚 문헌](systemic-lupus-erythematosus/sle_references.md) · [📄 README](systemic-lupus-erythematosus/README.md) |
| 110 | 자가면역·류마티스 | [**전신경화증**<br><sub>Systemic Sclerosis · SSc</sub>](systemic-sclerosis/) | <a href="systemic-sclerosis/ssc_qsp_model.svg"><img src="systemic-sclerosis/ssc_qsp_model.png" width="190" alt="SSc"></a> | 혈관병증·자가면역·섬유화의 삼중 병태. 면역억제·항섬유화·혈관확장제.<br>[🗺️ 지도](systemic-sclerosis/ssc_qsp_model.svg) · [⚙️ mrgsolve](systemic-sclerosis/ssc_mrgsolve_model.R) · [📚 문헌](systemic-sclerosis/ssc_references.md) · [📄 README](systemic-sclerosis/README.md) |
| 111 | 혈관염 | [**다카야스 동맥염**<br><sub>Takayasu Arteritis · TA</sub>](takayasu-arteritis/) | <a href="takayasu-arteritis/ta_qsp_model.svg"><img src="takayasu-arteritis/ta_qsp_model.png" width="190" alt="TA"></a> | 대동맥/주요 분지의 육아종성 대혈관염·협착. 스테로이드·토실리주맙.<br>[🗺️ 지도](takayasu-arteritis/ta_qsp_model.svg) · [⚙️ mrgsolve](takayasu-arteritis/ta_mrgsolve_model.R) · [📚 문헌](takayasu-arteritis/ta_references.md) · [📄 README](takayasu-arteritis/README.md) |
| 112 | 내분비·대사 | [**제1형 당뇨병**<br><sub>Type 1 Diabetes · T1DM</sub>](type1-diabetes/) | <a href="type1-diabetes/t1dm_qsp_model.svg"><img src="type1-diabetes/t1dm_qsp_model.png" width="190" alt="T1DM"></a> | 자가면역 베타세포 파괴에 의한 인슐린 결핍. 인슐린 요법·테플리주맙(발병 지연).<br>[🗺️ 지도](type1-diabetes/t1dm_qsp_model.svg) · [⚙️ mrgsolve](type1-diabetes/t1dm_mrgsolve_model.R) · [📚 문헌](type1-diabetes/t1dm_references.md) · [📄 README](type1-diabetes/README.md) |
| 113 | 내분비·대사 | [**제2형 당뇨병**<br><sub>Type 2 Diabetes · T2DM</sub>](type2-diabetes/) | <a href="type2-diabetes/t2dm_qsp_model.svg"><img src="type2-diabetes/t2dm_qsp_model.png" width="190" alt="T2DM"></a> | 인슐린 저항성·베타세포 기능부전. 메트포르민·GLP-1·SGLT2i.<br>[🗺️ 지도](type2-diabetes/t2dm_qsp_model.svg) · [⚙️ mrgsolve](type2-diabetes/t2dm_mrgsolve_model.R) · [📚 문헌](type2-diabetes/t2dm_references.md) · [📄 README](type2-diabetes/README.md) |
| 114 | 소화기·간담도 | [**궤양성 대장염**<br><sub>Ulcerative Colitis · UC</sub>](ulcerative-colitis/) | <a href="ulcerative-colitis/uc_qsp_model.svg"><img src="ulcerative-colitis/uc_qsp_model.png" width="190" alt="UC"></a> | 대장 점막의 Th2/장벽 염증. 5-ASA·항TNF·항인테그린·JAK 억제제.<br>[🗺️ 지도](ulcerative-colitis/uc_qsp_model.svg) · [⚙️ mrgsolve](ulcerative-colitis/uc_mrgsolve_model.R) · [📚 문헌](ulcerative-colitis/uc_references.md) · [📄 README](ulcerative-colitis/README.md) |
| 115 | 신장·비뇨 | [**요로결석 (만성 재발성)**<br><sub>Urolithiasis · URI</sub>](urolithiasis/) | <a href="urolithiasis/uri_qsp_model.svg"><img src="urolithiasis/uri_qsp_model.png" width="190" alt="URI"></a> | 소변 과포화·결정화·Randall 플라크에 의한 결석. 티아지드·구연산칼륨·알로푸리놀.<br>[🗺️ 지도](urolithiasis/uri_qsp_model.svg) · [⚙️ mrgsolve](urolithiasis/uri_mrgsolve_model.R) · [📚 문헌](urolithiasis/uri_references.md) · [📄 README](urolithiasis/README.md) |
| 116 | 피부 | [**백반증**<br><sub>Vitiligo · VIT</sub>](vitiligo/) | <a href="vitiligo/vit_qsp_model.svg"><img src="vitiligo/vit_qsp_model.png" width="190" alt="VIT"></a> | CD8 T세포·IFN-γ-CXCL10 축에 의한 멜라닌세포 파괴. JAK 억제제(룩솔리티닙) 색소 재침착.<br>[🗺️ 지도](vitiligo/vit_qsp_model.svg) · [⚙️ mrgsolve](vitiligo/vit_mrgsolve_model.R) · [📚 문헌](vitiligo/vit_references.md) · [📄 README](vitiligo/README.md) |
| 117 | 종양·호흡기 | [**비소세포 폐암 (NSCLC)**<br><sub>Non-Small Cell Lung Cancer · NSCLC</sub>](non-small-cell-lung-cancer/) | <a href="non-small-cell-lung-cancer/nsclc_qsp_model.svg"><img src="non-small-cell-lung-cancer/nsclc_qsp_model.png" width="190" alt="NSCLC"></a> | EGFR/KRAS/ALK 드라이버 돌연변이 → oncogenic signaling. TKI(오시머티닙·알렉티닙·소토라십)·면역관문억제제(펨브롤리주맙)·항암화학요법 PK/PD. FLAURA·ALEX·KEYNOTE-189 보정.<br>[🗺️ 지도](non-small-cell-lung-cancer/nsclc_qsp_model.svg) · [⚙️ mrgsolve](non-small-cell-lung-cancer/nsclc_mrgsolve_model.R) · [📚 문헌](non-small-cell-lung-cancer/nsclc_references.md) · [📄 README](non-small-cell-lung-cancer/README.md) |
| 118 | 감염/간담도 | [**만성 C형 간염**<br><sub>Chronic Hepatitis C · CHC/HCV</sub>](chronic-hepatitis-c/) | <a href="chronic-hepatitis-c/HCV_qsp_model.svg"><img src="chronic-hepatitis-c/HCV_qsp_model.png" width="190" alt="HCV"></a> | Perelson 표적세포 제한 바이러스 동역학(T/I/V ODE)에 DAA PK/PD 통합. SOF/LED·SOF/VEL·GLE/PIB·PEG-IFN/RBV 7개 시나리오. NS5B·NS5A·NS3 억제 효능(εp/εi), CTL 소진, 간섬유화(Metavir F-score), HCC 위험 모델링. ION·ASTRAL·ENDURANCE 임상시험 보정.<br>[🗺️ 지도](chronic-hepatitis-c/HCV_qsp_model.svg) · [⚙️ mrgsolve](chronic-hepatitis-c/HCV_mrgsolve_model.R) · [📚 문헌](chronic-hepatitis-c/HCV_references.md) · [📄 README](chronic-hepatitis-c/README.md) |
| 119 | 종양·간담도 | [**간세포암종 (HCC)**<br><sub>Hepatocellular Carcinoma · HCC</sub>](hepatocellular-carcinoma/) | <a href="hepatocellular-carcinoma/hcc_qsp_model.svg"><img src="hepatocellular-carcinoma/hcc_qsp_model.png" width="190" alt="HCC"></a> | HBV/HCV·NAFLD 기반 간암. RAS/RAF/MEK/ERK·PI3K/AKT/mTOR·Wnt/β-catenin·VEGF/혈관신생·종양면역 미세환경(PD-L1/CD8/Treg/TAM) 통합. 소라페닙·렌바티닙·아테조리주맙+베바시주맙(IMbrave150)·레고라페닙 PK/PD. 20구획 ODE, 5치료시나리오, AFP/간기능 바이오마커.<br>[🗺️ 지도](hepatocellular-carcinoma/hcc_qsp_model.svg) · [⚙️ mrgsolve](hepatocellular-carcinoma/hcc_mrgsolve_model.R) · [📊 Shiny](hepatocellular-carcinoma/hcc_shiny_app.R) · [📚 문헌](hepatocellular-carcinoma/hcc_references.md) · [📄 README](hepatocellular-carcinoma/README.md) |

| 120 | 신경계·유전 | [**헌팅턴병 (Huntington's Disease)**<br><sub>Huntington's Disease · HD</sub>](huntingtons-disease/) | <a href="huntingtons-disease/hd_qsp_model.svg"><img src="huntingtons-disease/hd_qsp_model.png" width="190" alt="HD"></a> | CAG 반복 확장(≥36) → mHTT 생성·집적 → 선조체 중간가시신경세포(MSN) 변성. mHTT 집적 폭포(단량체→올리고머→섬유)·BDNF-TrkB 결핍(REST/NRSF·HAP1 축)·흥분독성(eNMDAR/Ca²⁺/calpain)·미토콘드리아 기능부전(PGC-1α·Complex I/II/III)·신경염증(NLRP3/IL-1β)·아포토시스(Casp-3/6). VMAT2 억제제(TBZ·DTBZ·VBZ)·ASO(토미너센, CSF mHTT ↓74%)·스플라이싱 조절제(브라나플람·PTC518) PK/PD. 20구획 ODE, 7치료시나리오(TETRA-HD·FIRST-HD·KINECT-HD·GENERATION-HD1 임상 보정), UHDRS-TMS·TFC·cUHDRS·CAP score·CSF NfL/mHTT 바이오마커.<br>[🗺️ 지도](huntingtons-disease/hd_qsp_model.svg) · [⚙️ mrgsolve](huntingtons-disease/hd_mrgsolve_model.R) · [📊 Shiny](huntingtons-disease/hd_shiny_app.R) · [📚 문헌](huntingtons-disease/hd_references.md) · [📄 README](huntingtons-disease/README.md) |
| 121 | 종양·소화기 | [**대장암**<br><sub>Colorectal Cancer · CRC</sub>](colorectal-cancer/) | <a href="colorectal-cancer/crc_qsp_model.svg"><img src="colorectal-cancer/crc_qsp_model.png" width="190" alt="CRC"></a> | APC/Wnt·KRAS/MAPK·PI3K/TP53 경로 통합, MSI-H/CMS 아형, 종양 미세환경(CD8/Treg/MDSC/TAM), 혈관신생(VEGF/VEGFR). 5-FU/옥살리플라틴/이리노테칸(UGT1A1)·베바시주맙(TMDD)·세툭시맙·펨브롤리주맙(MSI-H) PK/PD. 20구획 ODE, 7치료시나리오(MOSAIC·CRYSTAL·NO16966·TRIBE·KEYNOTE-177 임상 보정), CEA·ctDNA·VEGF·EGFR/PD-1 점유율 바이오마커.<br>[🗺️ 지도](colorectal-cancer/crc_qsp_model.svg) · [⚙️ mrgsolve](colorectal-cancer/crc_mrgsolve_model.R) · [📊 Shiny](colorectal-cancer/crc_shiny_app.R) · [📚 문헌](colorectal-cancer/crc_references.md) · [📄 README](colorectal-cancer/README.md) |
| 122 | 종양·비뇨기 | [**전립선암**<br><sub>Prostate Cancer · PC</sub>](prostate-cancer/) | <a href="prostate-cancer/pc_qsp_model.svg"><img src="prostate-cancer/pc_qsp_model.png" width="190" alt="PC"></a> | HPG 축·AR 신호·PI3K/AKT·골전이(RANKL 악순환)와 CRPC/ARv7 내성. GnRH 제제·ARPI·도세탁셀·PARP 억제제·Lu-PSMA.<br>[🗺️ 지도](prostate-cancer/pc_qsp_model.svg) · [⚙️ mrgsolve](prostate-cancer/pc_mrgsolve_model.R) · [📚 문헌](prostate-cancer/pc_references.md) · [📄 README](prostate-cancer/README.md) |
| 123 | 종양·혈액 | [**급성 골수성 백혈병**<br><sub>Acute Myeloid Leukemia · AML</sub>](acute-myeloid-leukemia/) | <a href="acute-myeloid-leukemia/aml_qsp_model.svg"><img src="acute-myeloid-leukemia/aml_qsp_model.png" width="190" alt="AML"></a> | FLT3/NPM1/IDH/DNMT3A 돌연변이·BCL-2 패밀리·백혈병 줄기세포·골수 미세환경·후성유전학 이상. 베네토클락스(BCL-2)·아자시티딘·길테리티닙(FLT3)·에나시데닙(IDH2)·시타라빈 PK/PD. 21구획 ODE, Friberg 골수억제, 7치료 시나리오(VIALE-A·ADMIRAL·QuANTUM-R·RATIFY 임상 보정).<br>[🗺️ 지도](acute-myeloid-leukemia/aml_qsp_model.svg) · [⚙️ mrgsolve](acute-myeloid-leukemia/aml_mrgsolve_model.R) · [📊 Shiny](acute-myeloid-leukemia/aml_shiny_app.R) · [📚 문헌](acute-myeloid-leukemia/aml_references.md) · [📄 README](acute-myeloid-leukemia/README.md) |
| 124 | 종양·소화기 | [**위선암**<br><sub>Gastric Cancer · GC</sub>](gastric-cancer/) | <a href="gastric-cancer/gc_qsp_model.svg"><img src="gastric-cancer/gc_qsp_model.png" width="190" alt="GC"></a> | H. pylori CagA/VacA·HER2/FGFR2/MET·VEGF 혈관신생·CLDN18.2·PD-L1/PD-1 면역관문·TCGA 아형(EBV/MSI-H/GS/CIN)·TME. 트라스투주맙(TMDD)·라무시루맙(VEGFR2)·니볼루맙(PD-1)·T-DXd(HER2 ADC)·졸베툭시맙(CLDN18.2 ADCC)·FLOT/FOLFOX PK/PD. 28구획 ODE, Simeoni TGI, 6치료 시나리오(ToGA·RAINBOW·CheckMate649·SPOTLIGHT·DESTINY-Gastric01·FLOT4 임상 보정), CEA/CA19-9/ctDNA 바이오마커.<br>[🗺️ 지도](gastric-cancer/gc_qsp_model.svg) · [⚙️ mrgsolve](gastric-cancer/gc_mrgsolve_model.R) · [📊 Shiny](gastric-cancer/gc_shiny_app.R) · [📚 문헌](gastric-cancer/gc_references.md) · [📄 README](gastric-cancer/README.md) |
| 125 | 종양·혈액 | [**골수섬유증**<br><sub>Myelofibrosis · MF</sub>](myelofibrosis/) | <a href="myelofibrosis/mf_qsp_model.svg"><img src="myelofibrosis/mf_qsp_model.png" width="190" alt="MF"></a> | JAK2V617F/CALR/MPL 돌연변이·JAK/STAT3/5 신호·HSC 니치·골수 섬유화(TGF-β/PDGF/bFGF)·수외조혈(비장종대)·사이토카인 폭풍(IL-6/TNF-α)·혈전/출혈 위험. Ruxolitinib(2구획 PK)·Fedratinib·Pacritinib·Momelotinib·Pelabresib 병용 PK/PD. 23구획 ODE, 6치료 시나리오(COMFORT-I·JAKARTA·PERSIST-2·MANIFEST-2 임상 보정), SVR35/TSS50/VAF/pSTAT3·5 바이오마커.<br>[🗺️ 지도](myelofibrosis/mf_qsp_model.svg) · [⚙️ mrgsolve](myelofibrosis/mf_mrgsolve_model.R) · [📊 Shiny](myelofibrosis/mf_shiny_app.R) · [📚 문헌](myelofibrosis/mf_references.md) · [📄 README](myelofibrosis/README.md) |
| 126 | 종양·소화기 | [**췌장 선암**<br><sub>Pancreatic Ductal Adenocarcinoma · PDAC</sub>](pancreatic-cancer/) | <a href="pancreatic-cancer/pdac_qsp_model.svg"><img src="pancreatic-cancer/pdac_qsp_model.png" width="190" alt="PDAC"></a> | KRAS G12D/V/R(~95%)·TP53(~75%)·CDKN2A(~90%)·SMAD4(~55%) 드라이버 돌연변이, RAS/MAPK/PI3K/AKT/mTOR 신호, TGF-β/SMAD/EMT 축, 치밀 섬유화 기질(PSC·CAF·HA·IFP↑→약물침투↓), 고면역억제 TME(PD-L1·Treg·MDSC·TAM-M2), HIF-1α/VEGF 혈관신생, BRCA1/2/PALB2 HRD. 젬시타빈/nab-파클리탁셀·FOLFIRINOX/mFOLFIRINOX·KRAS G12D 억제제(MRTX1133)·올라파립 PK/PD. 264노드 10클러스터, 7치료 시나리오(MPACT·PRODIGE4·POLO 임상 보정), CA19-9·ctDNA·Friberg 골수억제 바이오마커.<br>[🗺️ 지도](pancreatic-cancer/pdac_qsp_model.svg) · [⚙️ mrgsolve](pancreatic-cancer/pdac_mrgsolve_model.R) · [📊 Shiny](pancreatic-cancer/pdac_shiny_app.R) · [📚 문헌](pancreatic-cancer/pdac_references.md) · [📄 README](pancreatic-cancer/README.md) |
| 127 | 종양·혈액 | [**골수이형성 증후군**<br><sub>Myelodysplastic Syndrome · MDS</sub>](myelodysplastic-syndrome/) | <a href="myelodysplastic-syndrome/mds_qsp_model.svg"><img src="myelodysplastic-syndrome/mds_qsp_model.png" width="190" alt="MDS"></a> | CHIP→MDS 클론 진화; SF3B1/SRSF2/U2AF1/ZRSR2 스플라이싱 돌연변이(고리 철아세포 형성·ABCB7 소실); TET2/DNMT3A/IDH1·2/ASXL1/EZH2 후성유전학 이상; del(5q)→RPS14 반수부족·miR-145/146a 소실·TIRAP↑; TP53 복잡핵형; 골수 미세환경 (CXCL12·SCF·TPO·EPO·GCSF 니치); GDF11/TGF-β1 상승→비효율 적혈구생성; 헵시딘↑/철과부하. 아자시티딘(2구획 SC/IV PK)·데시타빈(IV/oral ASTX727)·레날리도마이드(CRBN→CK1α 분해→RPS14 복원)·루스파터셉트(ActRIIB-Fc TGF-β 포획→Smad2/3↓)·다르베포에틴·베네토클락스(BCL-2 BH3 모방체) PK/PD. 324노드 10클러스터, 18구획 ODE, 7치료 시나리오(COMMANDS·MEDALIST·MDS-003/004·ASTX727·VIALE-A 임상 보정), Hgb/PLT/ANC/수혈 의존성·GDF11·Hepcidin·IronStore·IPSS-R/IPSS-M 바이오마커.<br>[🗺️ 지도](myelodysplastic-syndrome/mds_qsp_model.svg) · [⚙️ mrgsolve](myelodysplastic-syndrome/mds_mrgsolve_model.R) · [📊 Shiny](myelodysplastic-syndrome/mds_shiny_app.R) · [📚 문헌](myelodysplastic-syndrome/mds_references.md) · [📄 README](myelodysplastic-syndrome/README.md) |
| 128 | 근골격·신경 | [**섬유근통**<br><sub>Fibromyalgia · FM</sub>](fibromyalgia/) | <a href="fibromyalgia/fm_qsp_model.svg"><img src="fibromyalgia/fm_qsp_model.png" width="190" alt="FM"></a> | 중추감작(척수 WDR뉴런 LTP·wind-up·NMDA 수용체 탈억제); 하행성 통증조절계 결함(PAG-RVM-LC/Raphe 축, DPMS↓, DNIC 소실); 신경전달물질 불균형(CSF substance P↑·NE↓·5-HT↓); 신경염증(척수 미세아교세포 활성화·IL-1β/TNF-α·NLRP3 인플라마좀·KCC2↓→GABA 탈억제); HPA 축 이상(저코르티솔혈증·GH/IGF-1 결핍); 자율신경계 불균형(교감항진·HRV↓); 비회복성 수면(α파-δ파 침범·SWS↓); 섬유근통 뇌회로(ACC·섬엽·PFC 과활성화·기본모드네트워크). 둘록세틴(2구획 PK·SERT/NET IC50)·프레가발린(α2δ-1 Ca채널 차단·Emax 모델)·밀나시프란(SERT/NET)·아미트립틸린(H1 차단→수면 개선) PK/PD. 100+ 노드 10클러스터, **30구획 ODE**, 6치료 시나리오, Pain NRS·FIQR·CSF SP·미세아교세포·SWS·코르티솔·SNS tone 바이오마커.<br>[🗺️ 지도](fibromyalgia/fm_qsp_model.svg) · [⚙️ mrgsolve](fibromyalgia/fm_mrgsolve_model.R) · [📊 Shiny](fibromyalgia/fm_shiny_app.R) · [📚 문헌](fibromyalgia/fm_references.md) · [📄 README](fibromyalgia/README.md) |
| 129 | 혈액·응고 | [**혈우병 A**<br><sub>Hemophilia A · HA</sub>](hemophilia-a/) | <a href="hemophilia-a/ha_qsp_model.svg"><img src="hemophilia-a/ha_qsp_model.png" width="190" alt="HA"></a> | X염색체 연관 FVIII 결핍(F8 유전자 돌연변이) → 내인성 Xase 복합체(FIXa·FVIIIa) 형성 불능 → 트롬빈 생성 급감 → 불안정 피브린 클롯. 중증(<1 IU/dL) ABR ~30/년; 억제항체(30% 중증) 발생으로 대체치료 실패. 에미시주맙(HAVEN 1/3/4: FIXa–FX 이중특이항체·FVIII 모방; ABR 1.5/년)·피투시란(ATLAS-INH: siRNA·항트롬빈 mRNA 녹다운·AT 감소→트롬빈↑; ABR ~0)·마스타시맙(항TFPI·외인성경로 증폭)·유전자치료(AAV5-FVIII·valoctocogene roxaparvovec). SHL/EHL FVIII PK(2구획)·에미시주맙 SC PK(3구획·t½ 4–5주)·피투시란 PK/AT 간접반응 모델(mRNA/단백질 2단계 녹다운). **167 노드 10클러스터**, **16구획 ODE**, **7 치료 시나리오**(무예방·SHL-FVIII 3×/wk·EHL Q3-4d·에미시주맙 Q1W/Q4W·피투시란 Q1M·FVIII+에미시주맙 병용), ETP·ABR·Pettersson 관절점수·QoL·억제항체 역가 바이오마커. 55개 PubMed 인용(Manco-Johnson 2007·HAVEN·ATLAS·A-LONG·HOPE-B·WFH Guidelines 포함).<br>[🗺️ 지도](hemophilia-a/ha_qsp_model.svg) · [⚙️ mrgsolve](hemophilia-a/ha_mrgsolve_model.R) · [📊 Shiny](hemophilia-a/ha_shiny_app.R) · [📚 문헌](hemophilia-a/ha_references.md) · [📄 README](hemophilia-a/README.md) |
| 130 | 혈액·종양 | [**만성 림프구성 백혈병**<br><sub>Chronic Lymphocytic Leukemia · CLL</sub>](chronic-lymphocytic-leukemia/) | <a href="chronic-lymphocytic-leukemia/cll_qsp_model.svg"><img src="chronic-lymphocytic-leukemia/cll_qsp_model.png" width="190" alt="CLL"></a> | CD19⁺CD5⁺CD23⁺ 단클론 B세포 축적(서구권 성인 가장 흔한 백혈병). BCR 자율 신호전달(IGHV 비변이형·항원 자극)→LYN/SYK/BTK/PLCγ2→NF-κB/PI3Kδ/MAPK 경로 활성화; BCL-2 과발현(del13q→miR-15a/16-1 소실); 미세환경 의존(CXCL12/CXCR4 골수 억류·CXCL13/CXCR5 림프절 귀소·CD40L·BAFF·IL-4·NK 억제). CLL-IPI 예후(del(17p)/TP53·del(11q)/ATM·IGHV·β2M·임상병기). 이브루티닙(BTK Cys481 공유결합·420mg QD·RESONATE-2: ORR 86%·2yr PFS 74%)·아칼라브루티닙·자누브루티닙(선택성↑); 베네토클락스(BCL-2 BH3 모방·Ki~0.01nM·CLL14: uMRD 76%·2yr PFS 88%)·오비누투주맙(Type II 항CD20·ADCC↑·CDC↓·PCD). **146 노드 10클러스터**, **18구획 ODE**(이브루티닙 1구획·BTK 공유결합 모델·베네토클락스 2구획·BCL-2 준정상상태·오비누투주맙 TMDD·ALC/BM/LN 질환 구획·MCL-1 내성·NK 활성), **6치료 시나리오**(이브루티닙·베네토클락스·오비누투주맙·VEN+OBI CLL14·IB+VEN·삼중 병용). RESONATE-2·CLL14·MURANO·SEQUOIA·ALPINE 보정. 44개 PubMed 인용.<br>[🗺️ 지도](chronic-lymphocytic-leukemia/cll_qsp_model.svg) · [⚙️ mrgsolve](chronic-lymphocytic-leukemia/cll_mrgsolve_model.R) · [📊 Shiny](chronic-lymphocytic-leukemia/cll_shiny_app.R) · [📚 문헌](chronic-lymphocytic-leukemia/cll_references.md) · [📄 README](chronic-lymphocytic-leukemia/README.md) |
| 132 | 신경근육 | [**척수성 근위축증 (SMA)**<br><sub>Spinal Muscular Atrophy · SMA</sub>](spinal-muscular-atrophy/) | <a href="spinal-muscular-atrophy/sma_qsp_model.svg"><img src="spinal-muscular-atrophy/sma_qsp_model.png" width="190" alt="SMA"></a> | *SMN1* 5q13.2 결손→SMN 단백질 소실→알파 운동신경세포 진행성 사멸·NMJ 미성숙·신경원성 근위축. SMN2 대체 스플라이싱(엑손7 포함율 10%→90%)·SMN 단백질 역치·MN pool·NMJ 성숙도·근육량 ODE. 누시너센(IT ASO·ISS-N1 차단)·리스디플람(경구 스플라이싱 조절제)·오나셈노진(AAV9 유전자치료) 3종 PK/PD. **130+ 노드 12클러스터**, **20구획 ODE**, **6치료 시나리오**(ENDEAR·CHERISH·FIREFISH·SPR1NT 임상 보정). CHOP-INTEND·HFMSE·CMAP·FVC·혈청 NF-L 바이오마커. 50개 PubMed 인용.<br>[🗺️ 지도](spinal-muscular-atrophy/sma_qsp_model.svg) · [⚙️ mrgsolve](spinal-muscular-atrophy/sma_mrgsolve_model.R) · [📊 Shiny](spinal-muscular-atrophy/sma_shiny_app.R) · [📚 문헌](spinal-muscular-atrophy/sma_references.md) · [📄 README](spinal-muscular-atrophy/README.md) |
| 131 | 혈액·응고 | [**정맥 혈전색전증 (DVT/PE)**<br><sub>Venous Thromboembolism · VTE</sub>](venous-thromboembolism/) | <a href="venous-thromboembolism/vte_qsp_model.svg"><img src="venous-thromboembolism/vte_qsp_model.png" width="190" alt="VTE"></a> | 심부정맥 혈전증(DVT)과 폐색전증(PE)을 통합한 QSP 모델. Virchow's Triad(혈류정체·내피 손상·과응고), 외인성 경로(TF-FVIIa-TFPI), 내인성 경로(접촉활성화-FXIIa-FIXa-FVIIIa), 공통 경로(Prothrombinase-트롬빈-피브린 가교), 혈소판 활성화(GPIb/GPVI/PAR1/4·GPIIb/IIIa), 자연 항응고(AT-III·단백C/S·TFPI), 섬유용해(tPA/uPA·플라스민·PAI-1·TAFI·D-이량체). 리바록사반(2구획 PK·FXa EC50=12 ng/mL)·아픽사반(EC50=5 ng/mL)·다비가트란(직접트롬빈억제·EC50=35 ng/mL)·와파린(VK 사이클 간접반응·FVII/FX/FII 풀 반감기)·에녹사파린(AT-III 활성화·항Xa) 5종 약물 PK/PD. **140+ 노드 12클러스터**, **19구획 ODE**(PK 7·FXa/FIIa·피브린·혈전크기·플라스민·D-이량체·VK산화/환원·FVII/FX/FII풀), **6치료 시나리오**(DVT:리바록사반 15→20mg·PE:아픽사반 10→5mg BID·와파린+LMWH 브리지·수술예방:에녹사파린 40mg QD·확장예방:리바록사반 10mg QD·신부전:다비가트란 110mg BID GFR30 vs 90). EINSTEIN/AMPLIFY/RE-COVER/ROCKET-AF/Mueck 2011/Frost 2015 임상 파라미터 보정. INR·Anti-Xa·aPTT·D-이량체·혈전잔여% 바이오마커. Wells 점수(DVT/PE) 사전확률 계산기, 금기사항별 약물 추천 로직. 57개 PubMed 인용.<br>[🗺️ 지도](venous-thromboembolism/vte_qsp_model.svg) · [⚙️ mrgsolve](venous-thromboembolism/vte_mrgsolve_model.R) · [📊 Shiny](venous-thromboembolism/vte_shiny_app.R) · [📚 문헌](venous-thromboembolism/vte_references.md) · [📄 README](venous-thromboembolism/README.md) |
| 133 | 소아 혈관염 | [**가와사키병**<br><sub>Kawasaki Disease · KD</sub>](kawasaki-disease/) | <a href="kawasaki-disease/kd_qsp_model.svg"><img src="kawasaki-disease/kd_qsp_model.png" width="190" alt="KD"></a> | 원인 불명 트리거 → TLR/NLR 선천 면역 활성화 → NLRP3 인플라마좀(Caspase-1·IL-1β 성숙) → 사이토카인 폭풍(IL-1β·IL-6·TNF-α) → 혈관 내피 활성화(VCAM-1·ICAM-1·TF↑) → 관상동맥 중막 파괴·동맥류(AHA Z-점수 분류: small z≥2.5, medium z≥5, giant z≥10) → 혈소판 증가증(2주 피크) → 혈전위험. IVIG 2 g/kg(2구획+FcRn 재순환·Emax NF-κB 억제)·고용량 아스피린→저용량(COX-1/2 비가역)·메틸프레드니솔론(NF-κB 억제·GR 경로)·인플릭시맙 5 mg/kg(TNF-α 중화·n=1.8)·아나킨라 4 mg/kg/day(IL-1R 경쟁차단) 5종 PK/PD. **134 노드 14클러스터**, **21구획 ODE**(PK 11·IL1β·IL6·TNFα·대식세포·내피세포·발열·CRP·혈소판·관상동맥 Z-점수), **5치료 시나리오**(S1 표준IVIG·S2 고위험+스테로이드·S3 IVIG저항-2차IVIG·S4 인플릭시맙구제·S5 아나킨라구제). Kobayashi/Egami 위험점수 계산기·관상동맥 Z-점수 추적·IVIG 저항성 확률 모델. McCrindle/Kobayashi/KIDCARE Trial 보정. 60개 PubMed 인용.<br>[🗺️ 지도](kawasaki-disease/kd_qsp_model.svg) · [⚙️ mrgsolve](kawasaki-disease/kd_mrgsolve_model.R) · [📊 Shiny](kawasaki-disease/kd_shiny_app.R) · [📚 문헌](kawasaki-disease/kd_references.md) · [📄 README](kawasaki-disease/README.md) |
| 134 | 내분비·대사 | [**쿠싱 증후군**<br><sub>Cushing's Syndrome · CS</sub>](cushings-syndrome/) | <a href="cushings-syndrome/cs_qsp_model.svg"><img src="cushings-syndrome/cs_qsp_model.png" width="190" alt="CS"></a> | 뇌하수체 ACTH 선종(쿠싱병 70%)·이소성 ACTH(10%)·부신 선종(15%) 등에 의한 만성 고코르티솔혈증. 시상하부 일주기 CRH 리듬·CRHR1-PKA-CREB-POMC-ACTH 경로; USP8 탈유비퀴틴화(~50% 쿠싱병)·CDK4/6 세포증식; CYP11A1→CYP17A1→CYP21A2→CYP11B1 스테로이드 생합성; GR-α/HSP90/FKBP51·핵이동·GRE/nGRE·AP1/NF-κB 접촉억제·GILZ·SGK1; PEPCK/G6Pase↑·인슐린저항성·내장지방·근육위축·골다공증(RANKL↑/OPG↓)·RAAS 고혈압. **140+ 노드 13클러스터**, **21구획 ODE**(HPA 3+부신코르티솔 2+GR 3+대사 6+임상출력 1+약물PK 6), **6치료 시나리오**(자연경과·파시레오티드 0.6mg BID·케토코나졸 400mg BID·오실로드로스탯 5mg BID·미페프리스톤 600mg QD·수술 후 관해). PASPORT-CUSHINGS(Colao 2012 NEJM)·LINC 3/4(Feelders 2019/Pivonello 2020)·SEISMIC(Fleseriu 2012 JCEM) 임상 보정. UFC·LNSC·1mg DST·덱사메타손억제검사·BMD·HDRS-17. 8탭 Shiny 대시보드(환자프로파일·HPA/PK·스테로이드생합성·임상지표·시나리오비교·바이오마커·대사합병증·가상집단). 55개 PubMed 인용.<br>[🗺️ 지도](cushings-syndrome/cs_qsp_model.svg) · [⚙️ mrgsolve](cushings-syndrome/cs_mrgsolve_model.R) · [📊 Shiny](cushings-syndrome/cs_shiny_app.R) · [📚 문헌](cushings-syndrome/cs_references.md) · [📄 README](cushings-syndrome/README.md) |
| 135 | 희귀·유전질환 | [**유전성 혈관부종**<br><sub>Hereditary Angioedema · HAE</sub>](hereditary-angioedema/) | <a href="hereditary-angioedema/hae_qsp_model.svg"><img src="hereditary-angioedema/hae_qsp_model.png" width="190" alt="HAE"></a> | *SERPING1* 돌연변이 → C1-INH 결핍/기능이상(Type I/II) 또는 *F12* Thr328Lys 이득기능돌연변이(Type III) → 칼리크레인-키닌계(KKS) 무제한 활성화 → 브라디키닌(BK) 과잉 → B2R-Gq/IP3/Ca²⁺/eNOS/NO 경로 → VE-cadherin 소실·혈장삼출 → 피하/후두 혈관부종. FXII→FXIIa 접촉활성화, FXIIa→Prekal→Kal→HMWK 절단, BK→B2R/B1R 이중수용체 신호, C1-INH:FXIIa/Kal 복합체, 보체C4 소진, IL-1β·B1R 상향조절 염증 증폭. 이카티반트(B2R Ki=0.47nM)·C1-INH IV(Berinert 20 IU/kg)·에칼란티드·재조합 C1-INH 급성치료; 베로트랄스탓(IC50=3.7nM 칼리크레인 경구)·라나델루맙(KD<100pM prekallikrein SC)·C1-INH SC(Haegarda 60IU/kg 2×/wk) 예방치료. **120+ 노드 12클러스터**, **20구획 ODE**(PK 11+생물학적 9), **6치료 시나리오**(무치료·이카티반트·C1-INH IV·베로트랄스탓·라나델루맙·C1-INH SC). FAST-1/3(Cicardi 2010·Lumry 2011)·HELP OLE(Banerji 2020 87% 감소)·BELO(Farkas 2020 44% 감소)·CONFIDENT(Craig 2017 95% 감소) 보정. BK 농도·B2R 점유율·VP 지수·부종점수·C4·C1-INH% 바이오마커. 58개 PubMed 인용.<br>[🗺️ 지도](hereditary-angioedema/hae_qsp_model.svg) · [⚙️ mrgsolve](hereditary-angioedema/hae_mrgsolve_model.R) · [📊 Shiny](hereditary-angioedema/hae_shiny_app.R) · [📚 문헌](hereditary-angioedema/hae_references.md) · [📄 README](hereditary-angioedema/README.md) |
| 136 | 희귀·유전질환 | [**파브리병**<br><sub>Fabry Disease · FBR</sub>](fabry-disease/) | <a href="fabry-disease/fbr_qsp_model.svg"><img src="fabry-disease/fbr_qsp_model.png" width="190" alt="FBR"></a> | X-연관 리소소말 저장 질환. *GLA* 변이 → α-갈락토시다제 A(α-Gal A) 결핍 → 글로보트리아오실세라미드(Gb3)·lyso-Gb3 조직축적 → 신세뇨관·심근·신경·혈관내피 손상. Gb3 합성(UGCG·B4GALT5)·M6PR/리소솜 가수분해·GCS 억제(SRT). 신장(eGFR↓·UPCR↑·FSGS)·심장(LVH·HCM→DCM·부정맥)·신경(신경병성 통증·뇌졸중)·피부(혈관각화종) 다장기 병증. 아갈시다제 베타(1 mg/kg Q2W·ERT 표준)·아갈시다제 알파(0.2 mg/kg Q2W)·페구니알시다제 알파(1 mg/kg Q4W)·미갈라스타트(150 mg QOD·적합 변이 전용 경구 샤페론)·루세라스탓(1000 mg TID·GCS억제제 SRT) 5종 PK/PD. **138 노드 14클러스터**, **22구획 ODE**, **6치료 시나리오**(자연경과·아갈시다제 베타·알파·미갈라스타트·페구니알시다제 알파·ERT+루세라스탓). FABRY-001(Eng 2001 NEJM)·ATTRACT(Germain 2016 NEJM)·BRIGHT(Schiffmann 2021 JAMA)·MODIFY(Lenders 2022) 보정. eGFR·UPCR·LVMi·lyso-Gb3·통증 바이오마커. 60개 PubMed 인용.<br>[🗺️ 지도](fabry-disease/fbr_qsp_model.svg) · [⚙️ mrgsolve](fabry-disease/fbr_mrgsolve_model.R) · [📊 Shiny](fabry-disease/fbr_shiny_app.R) · [📚 문헌](fabry-disease/fbr_references.md) · [📄 README](fabry-disease/README.md) |
| 137 | 심혈관 | [**심근염**<br><sub>Myocarditis · MYO</sub>](myocarditis/) | <a href="myocarditis/myo_qsp_model.svg"><img src="myocarditis/myo_qsp_model.png" width="190" alt="MYO"></a> | 바이러스(CVB3·SARS-CoV-2·아데노바이러스·HHV-6) 심근세포 감염 → CAR/ACE2 수용체 진입 → 바이러스 복제·TLR3/7/9·RIG-I·MDA5·cGAS-STING 선천면역 활성화 → NF-κB·IRF3/7·NLRP3 인플라마좀 → IFN-α/β/γ·TNF-α·IL-1β·IL-6 사이토카인 폭풍 → NK세포·M1대식세포 심근 손상. CD4+ Th1/Th17 및 CD8+ CTL 적응 면역; 분자 유사성(molecular mimicry) → 항심근미오신·항β1-AR·항ANT·항TnI 자가항체 → ADCC·보체 활성화. TGF-β→근섬유아세포→콜라겐 침착→LV 확장·DCM. IVIG 2g/kg(IV 주입·2구획 PK·t½=21일·Fc-R 차단)·프레드니솔론(경구 F=80%·t½=2-3h·GR 억제)·아자티오프린→6-MP(전구약물 F=47%·HPRT/TPMT경로)·사이클로스포린(F=35%·칼시뉴린억제·t½=8-12h)·콜히친(F=45%·NLRP3/튜불린억제·Vd=250 L/kg) 5종 PK/PD. **170+ 노드 10클러스터**, **35구획 ODE**(심근세포 3+바이러스 1+선천면역 3+사이토카인 7+적응면역 9+섬유화 3+바이오마커 3+약물PK 6), **5치료 시나리오**(자연경과·IVIG 단독·프레드니솔론+아자티오프린 TIMIC·IVIG+Pred+Aza+CsA 거대세포·IVIG+콜히친). TIMIC(Frustaci 2009)·IMAC-2(McNamara 2001)·Cooper 2007 거대세포 프로토콜 보정. 트로포닌 I·BNP·LVEF·CMR-LGE·부정맥 위험 바이오마커. 7탭 Shiny 대시보드. 60개 PubMed 인용.<br>[🗺️ 지도](myocarditis/myo_qsp_model.svg) · [⚙️ mrgsolve](myocarditis/myo_mrgsolve_model.R) · [📊 Shiny](myocarditis/myo_shiny_app.R) · [📚 문헌](myocarditis/myo_references.md) · [📄 README](myocarditis/README.md) |
| 138 | 희귀·유전질환 | [**고셔병**<br><sub>Gaucher Disease · GCD</sub>](gaucher-disease/) | <a href="gaucher-disease/gcd_qsp_model.svg"><img src="gaucher-disease/gcd_qsp_model.png" width="190" alt="GCD"></a> | *GBA1* 이중대립변이 → 리소솜 β-글루코세레브로시다제(GBA) 결핍 → 글루코세레브로사이드(GC) 및 고독성 탈아실화 유도체 lyso-GL1 조직대식세포 축적 → 고셔세포 형성. GBA합성(ER→Golgi·M6P수용체·ERAD)·GCS 기질합성·M6P수용체 매개 ERT 리소솜 전달; GC→비장/간/골수/CNS 구획 이동; NF-κB→IL-1β·IL-6·TNF-α·MIP-1α·RANKL; 비장비대·간비대·빈혈·혈소판감소·골밀도감소. GBA-파킨슨 연계(α-시누클레인 축적). 이미글루세라제(60U/kg Q2W·2구획·CL=1.4L/h/kg·EMAX=85%)·벨라글루세라제α; 엘리글루스타트(84mg BID·GCS IC50=10nM·CYP2D6·EM/PM)·미글루스타트(IC50=50μM); 암브록솔(샤페론). **115+ 노드 10클러스터**, **26구획 ODE**(약물PK 8+효소/기질 5+바이오마커 4+장기용적 2+혈액 2+골 3+염증 2), **6치료 시나리오**(자연경과·이미글루세라제·벨라글루세라제α·엘리글루스타트EM·엘리글루스타트PM·저용량ERT+엘리글루스타트). Barton 1991 NEJM·Mistry 2015 JAMA·Zimran 2010 Blood·Balwani 2021 AJH 보정. GL-1·lyso-GL1·키토트리오시다제·페리틴·SV·LV·Hb·PLT·BMD 바이오마커. 9탭 Shiny 대시보드. 62개 PubMed 인용.<br>[🗺️ 지도](gaucher-disease/gcd_qsp_model.svg) · [⚙️ mrgsolve](gaucher-disease/gcd_mrgsolve_model.R) · [📊 Shiny](gaucher-disease/gcd_shiny_app.R) · [📚 문헌](gaucher-disease/gcd_references.md) · [📄 README](gaucher-disease/README.md) |
| 139 | 내분비·대사 | [**원발성 알도스테론증**<br><sub>Primary Aldosteronism · PA</sub>](primary-aldosteronism/) | <a href="primary-aldosteronism/pa_qsp_model.svg"><img src="primary-aldosteronism/pa_qsp_model.png" width="190" alt="PA"></a> | **Conn 증후군** — 부신 피질에서의 자율적 알도스테론 과분비(APA·BAH) → 레닌 억제·ARR 상승·Na⁺ 저류·K⁺ 소실·대사성 알칼리증. KCNJ5/CACNA1D/ATP1A1/ATP2B3 체성 돌연변이 → Ca²⁺ 내유 → CYP11B2(알도스테론 합성효소) 과발현 → 자율적 알도스테론 생성; RAAS 캐스케이드(레닌→AngI→AngII→알도스테론) + APA 자율분비. MR→SGK1→Nedd4-2 인산화 → ENaC 세포표면 발현↑ → Na⁺ 재흡수·ROMK K⁺ 분비·H⁺ 분비 → 저칼륨혈증·대사성 알칼리증. 부피팽창→MAP↑, 심근/혈관 MR 직접 활성→심근섬유화·LVH. 진단: ARR(≥30)·PAC(>15 ng/dL)·부신정맥 채혈(AVS). 복강경 부신절제술(APA 단측); 스피로노락톤(IC50=1.2 μg/L·활성대사체 카렌오논 t½≈20h)·에플레레논(IC50=2.5 μg/L)·파이네레논(IC50=0.65 μg/L·비스테로이드성·심장보호 우월) MR 길항; ACEi·CCB 병용. **120+ 노드 10클러스터**, **23구획 ODE**(약물PK 6+RAAS 3+신장/이온 5+심혈관 2+신기능 1+장기손상 2+부신 2+바이오마커 2), **8치료 시나리오**(무치료 APA 2년 진행·부신절제술·스피로노락톤 100mg·에플레레논 100mg·파이네레논 20mg·스피로+암로디핀·정상 대조·ACEi). Choi 2011 Science(KCNJ5)·Rossi 2006 JACC·Milliez 2005 JACC·Pitt 1999 NEJM(RALES)·Bakris 2020 NEJM(FIDELIO) 보정. ARR·PAC·PRA·K⁺·HCO₃⁻·MAP·LVMi·심근섬유화·GFR·HOMA proxy 바이오마커. 7탭 Shiny 대시보드(환자프로파일·RAAS/PK·알도스테론 패널·이온 항상성·심혈관/장기손상·시나리오 비교·바이오마커 탐색기). 48개 PubMed 인용.<br>[🗺️ 지도](primary-aldosteronism/pa_qsp_model.svg) · [⚙️ mrgsolve](primary-aldosteronism/pa_mrgsolve_model.R) · [📊 Shiny](primary-aldosteronism/pa_shiny_app.R) · [📚 문헌](primary-aldosteronism/pa_references.md) · [📄 README](primary-aldosteronism/README.md) |

---

## 🧬 파브리병 (Fabry Disease) — 최신 모델 상세 (2026-06-24)

> **디렉토리:** [`fabry-disease/`](fabry-disease/) | **약어:** FBR | **날짜:** 2026-06-24

[![FBR QSP 기계론적 지도](fabry-disease/fbr_qsp_model.png)](fabry-disease/fbr_qsp_model.svg)

**질환**: 파브리병(Fabry Disease, FBR) | **분류**: X-연관 리소소말 저장 질환 | **유병률**: 1:40,000 (고전형) | **효소**: α-갈락토시다제 A (α-Gal A)

### 핵심 기전 (14개 클러스터)

| 클러스터 | 핵심 기전 |
|---------|-----------|
| 1. 유전적 기반 | GLA 유전자(Xq22.1) 변이(미스센스 ~60%); X-연관 유전; 반접합 남성 고전형 vs 이형접합 여성 보인자(Lyon X-불활성화); 적합 변이 ~40%(미갈라스타트 적응증) |
| 2. α-Gal A 효소 생물학 | ER 합성 → 골지체 M6P 인산화 → CI-MPR/CD-MPR 수용체 → 리소솜 전달 → pH 4.5–5.0 최적 Gb3 가수분해 |
| 3. 당스핑고지질 대사 | 세라미드 → GCS/UGCG → GlcCer → LacCer → A4GalT → Gb3 합성; Gb3 탈아실화 → lyso-Gb3(독성 신호분자, 민감 바이오마커) |
| 4. 신장 병리 | 족세포 Gb3 → 족돌기 소실 → UPCR↑ → TGF-β 섬유화 → FSGS → eGFR 감소(-3~-12/yr 무치료) → ESRD |
| 5. 심장 병리 | 심근세포 Gb3 → LVMi 증가 → 이완기 기능장애 → 심근 섬유화(LGE) → 부정맥(SVT·VF·완전방실차단) → 급사 |
| 6. 신경계 병리 | DRG Gb3 → 소섬유 신경병증 → 신경병성 통증(BPI-SF) + CNS 혈관내피 → 백질 병변 → TIA/뇌경색 |
| 7. 기타 장기 | 혈관각화종(피부), 무한증, 각막 소용돌이(cornea verticillata), 위장관 운동장애, 감각신경성 난청 |
| 8. 염증 폭포 | Lyso-Gb3 → TLR4 → NF-κB → IL-6/TNF-α → NLRP3 인플라마좀 → eNOS↓ → 내피세포 활성화 |
| 9. ERT PK/PD | 아갈시다제 베타(1 mg/kg Q2W, t½ ~45min, M6P→리소솜, Emax ~80%), 아갈시다제 알파(0.2 mg/kg Q2W), 페구니알시다제 알파(1 mg/kg Q4W, PEG t½ ~80h) |
| 10. 미갈라스타트 (샤페론) | 150 mg QOD 경구; 잘못 접힌 α-Gal A 안정화; EC50 ~0.25 μg/mL; 적합 변이 전용; ATTRACT 임상(ERT 비열등) |
| 11. 기질감소요법 (SRT) | 루세라스탓(GCS IC50 ~0.18 μg/mL, Emax 42%), 벵글루스탓(CNS 투과); Gb3 상류 기질 감소 |
| 12. 바이오마커 | 혈장 lyso-Gb3(μg/L, 가장 민감), 소변 Gb3(nmol/mg Cr), 백혈구 α-Gal A 활성, DBS 신생아 선별, eGFR, LVMi, UPCR |
| 13. 임상 엔드포인트 | eGFR 기울기(mL/min/1.73m²/yr), UPCR, LVMi, BPI-SF 통증(0–10), EQ-5D QoL, MSSI 중증도(0–138) |
| 14. 자연 경과 | 고전형 남성(소아기 발현) vs 후기 발현형(심장·신장형) vs 여성 보인자; Fabry Registry/FOS 데이터; 진단 지연 평균 10–20년 |

### mrgsolve ODE 모델 (22구획)

| 모듈 | 구획 | 핵심 동역학 |
|------|------|------------|
| 아갈시다제 베타 PK | A_AGAB_C, A_AGAB_P, A_AGAB_LYS | 2구획+리소솜 전달; CL=0.42 L/h; M6PR-매개 k_lys=0.35/h |
| 아갈시다제 알파 PK | A_AGAA_C, A_AGAA_P, A_AGAA_LYS | 2구획; CL=0.55 L/h |
| 미갈라스타트 PK | A_MIG_GUT, A_MIG_C | 경구 ka=0.82/h, F=75%, t½~3.5h |
| 루세라스탓 PK | A_LUC_GUT, A_LUC_C | 경구 IC50_GCS=0.18 μg/mL, Emax=42% |
| α-Gal A 효소 | E_GalA | ERT(Emax_ERT=70) + 미갈라스타트(Emax_MIG=6, 적합 변이 전용) + 기저 잔여 |
| 당스핑고지질 | GB3_PLM, GB3_KID, GB3_HRT, LGB3_PLM | 합성-효소분해 ODE; SRT→GCS 억제→상류 감소 |
| 염증 | INFLAM | lyso-Gb3 구동 k_in – k_out ODE |
| 장기 기능 | eGFR, UPCR, LVMi, PAIN | Gb3 축적 의존 감소; ERT/샤페론 보호 효과 |

### 6가지 치료 시나리오 임상 근거

| 시나리오 | 약물 | 임상시험 | 주요 결과 |
|---------|------|---------|---------|
| S1: 자연경과 | 없음 | Mehta 2009; Fabry Registry | eGFR -3~-12/yr, lyso-Gb3 30–80 μg/L, ESRD·심장사 위험 |
| S2: 아갈시다제 베타 | 1 mg/kg IV Q2W | FABRY-001 (Eng 2001 NEJM); Banikazemi 2007 AIM | 복합 신장·심장·뇌 사건 61% 감소 |
| S3: 아갈시다제 알파 | 0.2 mg/kg IV Q2W | Schiffmann 2001 Ann Intern Med | 신경병성 통증 개선, 신기능 안정, 소변 Gb3 감소 |
| S4: 미갈라스타트 | 150 mg PO QOD | ATTRACT (Germain 2016 NEJM); Hughes 2017 Lancet | ERT 비열등 (적합 변이); eGFR 기울기 -0.3 vs -1.0 mL/min/yr |
| S5: 페구니알시다제 알파 | 1 mg/kg IV Q4W | BRIGHT (Schiffmann 2021 JAMA) | eGFR 안정, lyso-Gb3 -50%, 4주 1회 투여 편의성 |
| S6: ERT + 루세라스탓 | AgaB + 1000 mg TID | MODIFY (Lenders 2022 Lancet DE) | Gb3 추가 감소, BPI-SF 신경병성 통증 -1.5점 |

### QSP 모델 파일

| 구성요소 | 파일 | 사양 |
|---------|------|-----|
| 🗺️ 기계론적 지도 | [`fbr_qsp_model.dot`](fabry-disease/fbr_qsp_model.dot) | **138 노드, 14클러스터** |
| ⚙️ mrgsolve ODE | [`fbr_mrgsolve_model.R`](fabry-disease/fbr_mrgsolve_model.R) | **22구획 ODE**, **6치료 시나리오** |
| 📊 Shiny 앱 | [`fbr_shiny_app.R`](fabry-disease/fbr_shiny_app.R) | **8탭** (환자 프로파일·PK/효소·Gb3 동역학·신장·심장·시나리오 비교·바이오마커·가상 집단) |
| 📚 참고문헌 | [`fbr_references.md`](fabry-disease/fbr_references.md) | **60개 PubMed 인용** (14개 섹션) |
---

## 🫀 심근염 (Myocarditis) — 최신 모델 상세 (2026-06-24)

> **디렉토리:** [`myocarditis/`](myocarditis/) | **약어:** MYO | **날짜:** 2026-06-24

[![MYO QSP 기계론적 지도](myocarditis/myo_qsp_model.png)](myocarditis/myo_qsp_model.svg)

**질환**: 심근염(Myocarditis, MYO) | **분류**: 심근의 염증성 질환 | **발병률**: 1–10/100,000인년 | **주요 원인**: 바이러스(CVB3·SARS-CoV-2·아데노바이러스·HHV-6·파보바이러스 B19)

---

### 기전 개요 — 3단계 병태생리

| 단계 | 기전 | 핵심 매개자 |
|------|------|------------|
| **1. 바이러스 단계** (1–7일) | CAR/ACE2 수용체 진입 → 복제 → 심근세포 직접 손상 | CVB3 프로테아제 2A/3C, dsRNA, 트로포닌 유출 |
| **2. 면역 단계** (1–4주) | PRR 활성화 → 선천면역·적응면역 폭주 → 분자 유사성 자가항체 | IFN-β/γ, TNF-α, IL-1β, IL-6, NK세포, M1, CTL |
| **3. 리모델링 단계** (4주+) | TGF-β → 근섬유아세포 활성화 → 콜라겐 침착 → LV 확장 | TGF-β, MMP/TIMP, RAAS, 피브로넥틴, DCM |

---

### ODE 구획 구조 (35개)

| 범주 | 구획 | 핵심 동역학 |
|------|------|------------|
| 심근세포 (3) | H(정상)·I(감염)·D(사망) | 바이러스 감염률·CTL 용해·허용 세포 재생 |
| 바이러스 (1) | V (copies/mL) | 복제(p_V·I) – 제거(IFN-β·NK 살상) |
| 선천면역 (3) | NK·M1·M2 | V/IFN-γ 구동 팽창; M2→TGF-β |
| 사이토카인 (7) | IFN-β·IFN-γ·TNF-α·IL-6·IL-1β·TGF-β·IL-10 | 생성(면역세포)–분해(1차 동역학) |
| 적응면역 (9) | 나이브 CD4·Th1·Th17·Treg·나이브 CD8·CTL·나이브 B·형질세포·항체 | Th1→IFN-γ; CTL→직접 세포 용해; B→자가항체 |
| 섬유화 (3) | 심장 섬유아세포·근섬유아세포·콜라겐 | TGF-β→MF 전환 →콜라겐 합성 |
| 바이오마커 (3) | 트로포닌 I·BNP·LVEF | 세포 손상 유출; 압력 과부하; EF 회복 |
| 약물 PK (6) | IVIG·PRED_A/C·AZA_A/C·MP6_C·CSA_A/C·COLC_A/C | 2구획/경구 1구획; 전구약물 변환 |

### 핵심 기계론적 방정식

```
dV/dt  = p_V·I − c_V·V·(1 + kIFN·IFNβ/(IFNβ50+IFNβ)) − NK_kill·NK·V
dH/dt  = r_H·H·(1−(H+I)/Hmax) − d_H·H − βV·V·H − CTL_bys·CTL·H·(1−E_IS)
dTnI/dt = kLeak·(δI·I + kNec·IS·H) − dTn·TnI
dEF/dt = krec·(EF_target − EF)
```

---

### 5가지 치료 시나리오 임상 근거

| 시나리오 | 약물 | 임상시험 | 주요 결과 |
|---------|------|---------|---------|
| S1: 자연경과 | 없음 | Bozkurt 2021 JACC | 완전 회복 40–50%, DCM 진행 20–30% |
| S2: IVIG 단독 | 2 g/kg IV 1회 | IMAC-2 (McNamara 2001) | 바이러스성·급성 심부전; Fc-R 차단 |
| S3: 프레드니솔론+아자티오프린 | PRED 1 mg/kg/d + AZA 2 mg/kg/d | TIMIC (Frustaci 2009) | 바이러스 음성 염증성 CM; EF +8.5% |
| S4: 삼중 면역억제 | IVIG+PRED+AZA+CsA | Cooper 2007 (거대세포) | 거대세포 심근염; 생존 연장 |
| S5: IVIG+콜히친 | IVIG 2g/kg + 콜히친 0.5mg BID | 심근심낭염 프로토콜 | NLRP3·인플라마좀 억제; 재발 감소 |

---

### QSP 모델 파일

| 구성요소 | 파일 | 사양 |
|---------|------|-----|
| 🗺️ 기계론적 지도 | [`myo_qsp_model.dot`](myocarditis/myo_qsp_model.dot) | **170+ 노드, 10클러스터** |
| ⚙️ mrgsolve ODE | [`myo_mrgsolve_model.R`](myocarditis/myo_mrgsolve_model.R) | **35구획 ODE**, **5치료 시나리오** |
| 📊 Shiny 앱 | [`myo_shiny_app.R`](myocarditis/myo_shiny_app.R) | **7탭** (개요·PK 프로파일·바이러스/선천면역·PD 바이오마커·임상 엔드포인트·시나리오 비교·섬유화/리모델링) |
| 📚 참고문헌 | [`myo_references.md`](myocarditis/myo_references.md) | **60개 PubMed 인용** (14개 섹션) |

---

## 🧬 고셔병 (Gaucher Disease) — 최신 모델 상세 (2026-06-24)

> **디렉토리:** [`gaucher-disease/`](gaucher-disease/) | **약어:** GCD | **날짜:** 2026-06-24

[![GCD QSP 기계론적 지도](gaucher-disease/gcd_qsp_model.png)](gaucher-disease/gcd_qsp_model.svg)

**질환**: 고셔병(Gaucher Disease) | **분류**: 희귀·유전질환(리소솜 축적 질환) | **유병률**: 1/40,000 (일반); 1/800 (아슈케나지 유대인) | **원인유전자**: *GBA1* (1q22)

---

### 핵심 발병기전

| 단계 | 기전 | 핵심 매개자 |
|------|------|------------|
| **1. 효소 결핍** | GBA1 변이 → 미스폴딩 → ERAD 분해 or M6P 리소솜 전달 감소 | UGGT, 칼넥신, ERAD, M6PR |
| **2. 기질 축적** | GCS 합성 > GBA 분해 → GC·lyso-GL1 리소솜 축적 | GC_MAC, GC_SP, GC_LV, GC_BM |
| **3. 고셔세포 형성** | 조직 대식세포 GC 불완전 소화 → NF-κB 활성화 | TGF-β, M2 편극 |
| **4. 전신 염증** | IL-1β·IL-6·TNF-α·MIP-1α·RANKL 분비 | 키토트리오시다제, 페리틴 |
| **5. 장기 손상** | 비장/간 비대, 골수 조혈 억제, RANKL→OC→골밀도↓ | SV, LV, Hb, PLT, BMD |
| **6. GBA-PD 연계** | GBA↓ → α-syn 클리어런스↓ → 루이체 | α-Synuclein, Lewy Body |

---

### 치료 시나리오

| 시나리오 | 약물 | 임상시험 | 주요 결과 |
|---------|------|---------|---------|
| **S1** | 자연경과 | Charrow 2000 | 비장비대·빈혈·골파괴 진행 |
| **S2** | 이미글루세라제 60U/kg Q2W | Barton 1991 NEJM | GL-1 -70%, SV -30–50%, Hb +2 g/dL |
| **S3** | 벨라글루세라제α 60U/kg Q2W | Zimran 2010 Blood | 유사 효능, 천연형 만노스 |
| **S4** | 엘리글루스타트 84mg BID (CYP2D6 EM) | Mistry 2015 JAMA | SV -28%, Hb +1.2, PLT +32% |
| **S5** | 엘리글루스타트 84mg QD (CYP2D6 PM) | Balwani 2021 AJH | AUC 4–5×↑, 효능 유사 |
| **S6** | 저용량 ERT 30U/kg + 엘리글루스타트 BID | 병용 전략 | 이중 억제, 주사 빈도 절감 |

---

### QSP 모델 파일

| 구성요소 | 파일 | 사양 |
|---------|------|-----|
| 🗺️ 기계론적 지도 | [`gcd_qsp_model.dot`](gaucher-disease/gcd_qsp_model.dot) | **115+ 노드, 10클러스터** |
| ⚙️ mrgsolve ODE | [`gcd_mrgsolve_model.R`](gaucher-disease/gcd_mrgsolve_model.R) | **26구획 ODE**, **6치료 시나리오** |
| 📊 Shiny 앱 | [`gcd_shiny_app.R`](gaucher-disease/gcd_shiny_app.R) | **9탭** (개요·환자프로파일·PK·효소/기질·장기/혈액·골·시나리오비교·바이오마커·가상집단) |
| 📚 참고문헌 | [`gcd_references.md`](gaucher-disease/gcd_references.md) | **62개 PubMed 인용** (14개 섹션) |

---

## 🫀 원발성 알도스테론증 (Primary Aldosteronism) — 최신 모델 상세 (2026-06-24)

> **디렉토리:** [`primary-aldosteronism/`](primary-aldosteronism/) | **약어:** PA | **날짜:** 2026-06-24

[![PA QSP 기계론적 지도](primary-aldosteronism/pa_qsp_model.png)](primary-aldosteronism/pa_qsp_model.svg)

**질환**: 원발성 알도스테론증(Primary Aldosteronism, PA) / Conn 증후군 | **분류**: 내분비·대사(부신) | **유병률**: 고혈압 환자의 5–10% | **원인**: APA(선종 ~35%)·BAH(양측 과형성 ~60%)

---

### 병태생리 요약

| 단계 | 핵심 기전 | 주요 노드 |
|------|----------|----------|
| **1. 체성 돌연변이** | KCNJ5/CACNA1D/ATP1A1/ATP2B3 이온 채널·펌프 변이 → Ca²⁺ 내유 → CYP11B2 과발현 | KCNJ5, CACNA1D, CYP11B2 |
| **2. 자율 알도스테론 과분비** | APA·BAH 에서 AngII 비의존적 알도스테론 생성 → 레닌 억제 | Aldo_c, Renin_c, ARR_c |
| **3. ENaC/ROMK 신장 효과** | MR→SGK1→Nedd4-2 인산화 → ENaC 세포표면 발현↑ → Na⁺ 재흡수·K⁺ 분비 | ENaC_act, K_c, HCO3_c |
| **4. 심혈관 표적 장기 손상** | 알도스테론 직접 심근섬유아세포 MR 활성 → 콜라겐 침착·LVH | CardFib, LVMi_c |
| **5. 대사 합병증** | Na⁺ 과잉·IR 악화, 심혈관 위험 4배 상승 | HOMA_proxy, MAP_c |

---

### 치료 시나리오

| 시나리오 | 약물·중재 | 임상시험 | 주요 결과 |
|---------|---------|---------|---------|
| **S1** | 무치료 APA | 자연경과 | MAP↑, K⁺↓, ARR>>30, LVMi↑ |
| **S2** | 복강경 부신절제술 | Rossi 2013 Hypertension | MAP -15 mmHg, ARR 정상화, LVMi -20% |
| **S3** | 스피로노락톤 100 mg/d | Monticone 2015 JH | K⁺ 회복, ARR 정상화, BP -10–15 mmHg |
| **S4** | 에플레레논 100 mg/d | Pitt 2003 NEJM | 선택적 MRA, 항안드로겐 부작용 없음 |
| **S5** | 파이네레논 20 mg/d | Bakris 2020 NEJM(FIDELIO) | 비스테로이드성·심장섬유화 우월 억제 |
| **S6** | 스피로노락톤 + 암로디핀 | 병용 전략 | 추가 BP 강하 (-5–8 mmHg) |
| **S7** | 정상 대조 | 참조 | Renin=1.0, Aldo=8 ng/dL, K=4.0 |
| **S8** | 라미프릴(ACEi) 10 mg/d | 제한적 효과 | AngII↓ but APA 알도스테론 불변 |

---

### QSP 모델 파일

| 구성요소 | 파일 | 사양 |
|---------|------|-----|
| 🗺️ 기계론적 지도 | [`pa_qsp_model.dot`](primary-aldosteronism/pa_qsp_model.dot) | **120+ 노드, 10클러스터** |
| ⚙️ mrgsolve ODE | [`pa_mrgsolve_model.R`](primary-aldosteronism/pa_mrgsolve_model.R) | **23구획 ODE**, **8치료 시나리오** |
| 📊 Shiny 앱 | [`pa_shiny_app.R`](primary-aldosteronism/pa_shiny_app.R) | **7탭** (환자프로파일·RAAS/PK·알도스테론패널·이온항상성·심혈관/TOD·시나리오비교·바이오마커탐색기) |
| 📚 참고문헌 | [`pa_references.md`](primary-aldosteronism/pa_references.md) | **48개 PubMed 인용** (12개 섹션) |
