# QSP 질환 모델 라이브러리 (QSP Disease Model Library)

> 매일 **Claude Code Routine(CCR)** 이 질환 하나를 선택해 **정량적 시스템 약리학(Quantitative Systems Pharmacology, QSP)** 모델을 처음부터 끝까지 구축하고 `main`에 직접 커밋하는, **살아 있는(living) 오픈 모델 라이브러리**입니다.

![models](https://img.shields.io/badge/models-127-blue) ![framework](https://img.shields.io/badge/QSP-mrgsolve%20%C2%B7%20Shiny%20%C2%B7%20Graphviz-success) ![automation](https://img.shields.io/badge/built%20by-Claude%20Code%20Routine-orange)

현재 **127개 질환**에 대한 완성된 QSP 모델이 수록되어 있으며, 각 모델은 ①기계론적 지도, ②mrgsolve ODE 모델, ③Shiny 대시보드, ④참고문헌의 네 가지 산출물로 구성됩니다. 아래 [모델 갤러리](#-모델-갤러리-model-gallery)에서 전체 목록을 확인할 수 있습니다.

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

전체 **122개** QSP 모델입니다. 모델명을 클릭하면 해당 디렉토리로, 그림을 클릭하면 확대 가능한 SVG 지도로 이동합니다. 각 행의 링크에서 기계론적 지도(🗺️), mrgsolve 모델(⚙️), 참고문헌(📚), 상세 README(📄)에 바로 접근할 수 있습니다.

**분류별 모델 수**: 내분비·대사 20 · 소화기·간담도 17 · 자가면역·류마티스 13 · 심혈관 11 · 신장·비뇨 11 · 신경 11 · 호흡기 8 · 혈관염 8 · 피부 6 · 혈액 5 · 종양 3 · 정신·신경 2 · 감염 2

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


---

## 🧬 헌팅턴병 (Huntington's Disease) — 최신 모델 상세 (2026-06-23)

[![헌팅턴병 QSP 기계론적 지도](huntingtons-disease/hd_qsp_model.png)](huntingtons-disease/hd_qsp_model.svg)

**질환**: 헌팅턴병(HD) | **유전자**: *HTT* (CAG 반복 ≥36) | **OMIM**: [143100](https://omim.org/entry/143100)

### 주요 경로 클러스터 (14개)

| 클러스터 | 핵심 기전 |
|---------|-----------|
| 1. 유전적 기반 | CAG 길이 의존적 mHTT mRNA 생성, polyQ 번역 |
| 2. mHTT 집적 폭포 | 단량체 → 올리고머(독성) → 원섬유 → 핵내 포함체 |
| 3. 단백질 품질관리 | UPS(포화), 자가포식(LC3/Beclin-1/p62), 리소솜/TFEB |
| 4. 전사 조절이상 | REST/NRSF, CBP/HDAC, PGC-1α, BDNF 유전자 억제 |
| 5. BDNF-TrkB 신경영양 신호 | HAP1 매개 축삭 수송, TrkB/PI3K/Akt 생존 경로 |
| 6. 미토콘드리아 기능부전 | Complex I/II/III 손상, ROS, ΔΨm 저하, 사이토크롬 c 방출 |
| 7. 흥분독성 | 시냅스외 NMDAR 과민화, Ca²⁺/calpain/nNOS 연쇄 |
| 8. 선조체 회로 | D1/D2-MSN 도파민·GABA 회로 파괴 |
| 9. 신경염증 | 미세아교세포 M1 TLR4/NF-κB → IL-1β/TNF-α → NLRP3 인플라마솜 |
| 10. 아포토시스 | BAX/Casp-3/6, p53/PUMA, AIF → 선조체 용적 소실 |
| 11. 약물 PK | TBZ/DTBZ/VBZ(VMAT2), 토미너센(IT-ASO), 브라나플람, 리루졸 |
| 12. 증상 완화 PD | VMAT2 억제 → DA 고갈 → 무도증 감소 |
| 13. 질환 수정 치료 | HTT 발현 저하: ASO(mRNA ↓74%), 스플라이싱(mRNA ↓50%) |
| 14. 임상 지표 | UHDRS-TMS, TFC, cUHDRS, CAP score, CSF NfL/mHTT, MRI 위축 |

### 치료 시나리오 (7종)

| # | 시나리오 | 약물·용법 | 기전 | 임상 근거 |
|---|---------|----------|------|----------|
| 1 | 자연경과 | — | — | ENROLL-HD · TRACK-HD |
| 2 | TBZ 25 mg/일 | 테트라베나진 TID | VMAT2 억제 → 도파민 고갈 | TETRA-HD (NEJM 2008) |
| 3 | DTBZ 30 mg/일 | 도이테트라베나진 BID | d-KIE 중간대사체 → VMAT2 억제 | FIRST-HD (JAMA 2016) |
| 4 | VBZ 80 mg QD | 발베나진 | VMAT2 억제 (1일 1회) | KINECT-HD (NEJM 2023) |
| 5 | 토미너센 120 mg Q8W | 척수강내 ASO | RNase H1 → mHTT mRNA ↓74% | GENERATION-HD1 (NEJM 2022) |
| 6 | 브라나플람 50 mg QW | 경구 스플라이싱 조절제 | Exon 49 스킵 → NMD → mHTT ↓50% | NCT04000594 (진행 중) |
| 7 | 병용: DTBZ + 토미너센 | 두 약물 병용 | 증상 완화 + 질환 수정 | 모델 시뮬레이션 |

### 산출물

| 파일 | 내용 |
|------|------|
| [`hd_qsp_model.dot`](huntingtons-disease/hd_qsp_model.dot) | Graphviz 소스 (14 클러스터, 110+ 노드) |
| [`hd_qsp_model.svg`](huntingtons-disease/hd_qsp_model.svg) | 벡터 기계론적 지도 |
| [`hd_qsp_model.png`](huntingtons-disease/hd_qsp_model.png) | 래스터 지도 (150 dpi) |
| [`hd_mrgsolve_model.R`](huntingtons-disease/hd_mrgsolve_model.R) | 20구획 ODE (7 치료 시나리오) |
| [`hd_shiny_app.R`](huntingtons-disease/hd_shiny_app.R) | 6탭 Shiny 대시보드 |
| [`hd_references.md`](huntingtons-disease/hd_references.md) | 53개 PubMed 참고문헌 |

---

## ⚠️ 면책 조항 (Disclaimer)

본 라이브러리의 모든 모델은 **교육 및 연구 목적의 정성적·반정량적 QSP 모델**입니다. 공개 문헌과 임상시험 데이터를 바탕으로 구성되었으나 독립적으로 검증·인증되지 않았으며, **실제 임상 의사결정, 처방, 또는 규제 제출에 직접 사용해서는 안 됩니다.** 파라미터와 가정은 설명을 위한 근사치이며, 실제 환자 데이터에 대한 적합·검증이 별도로 필요합니다.

## 📖 참고 자료 (References & Tools)

- mrgsolve를 이용한 R 기반 QSP: <https://vantage-research.net/qsp-in-r/>
- gPKPDviz — mrgsolve 기반 PK/PD 시뮬레이션 Shiny 도구
  - 논문: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/>
  - 코드: <https://github.com/Genentech/gPKPDviz/>

## 📄 라이선스 (License)

본 저장소의 라이선스는 [LICENSE](LICENSE) 파일을 참조하세요.

---

## 🧬 119. 간세포암종 (Hepatocellular Carcinoma, HCC)

**간세포암종(HCC)** 은 전 세계 원발성 간암의 약 90%를 차지하며, 만성 간질환(HBV·HCV·알코올성 간경변·NAFLD/NASH)이 있는 조직에서 발생한다. 연간 약 80만 명이 사망하는 세계 4위 암 사망 원인이다.

### 핵심 병태생리 경로

| 경로 | 주요 구성요소 | 치료 표적 |
|------|-------------|----------|
| **MAPK 신호** | RAS→RAF→MEK→ERK | 소라페닙, 레고라페닙 (RAF 억제) |
| **PI3K/AKT/mTOR** | EGFR/MET/IGF-1R 하류 | 내성 우회 경로 (병용 필요) |
| **Wnt/β-catenin** | CTNNB1 돌연변이, APC 파괴복합체 | ~30% 간암에서 활성 |
| **VEGF/혈관신생** | HIF-1α→VEGF-A/VEGFR2, PDGF-BB/PDGFRβ | 렌바티닙, 베바시주맙 |
| **종양면역 미세환경** | PD-1/PD-L1, CD8+ CTL 소진, TAM-M2, MDSC, Treg | 아테조리주맙 (IMbrave150) |

### QSP 모델 구성 (4종 산출물)

| 산출물 | 파일 | 내용 |
|--------|------|------|
| 🗺️ 기계론적 지도 | [`hcc_qsp_model.dot/.svg/.png`](hepatocellular-carcinoma/) | 120+ 노드, 9개 클러스터 (위험인자·발암·신호전달·혈관신생·종양·면역·PK·PD·임상지표) |
| ⚙️ mrgsolve ODE | [`hcc_mrgsolve_model.R`](hepatocellular-carcinoma/hcc_mrgsolve_model.R) | 20구획 (5약물 PK + 종양/면역/VEGF/AFP/간기능), 5치료 시나리오 |
| 📊 Shiny 앱 | [`hcc_shiny_app.R`](hepatocellular-carcinoma/hcc_shiny_app.R) | 6탭 (환자 프로파일·PK·종양·PD/바이오마커·시나리오 비교·면역/안전성) |
| 📚 참고문헌 | [`hcc_references.md`](hepatocellular-carcinoma/hcc_references.md) | 43개 PubMed 인용 (SHARP·REFLECT·IMbrave150·RESORCE·CELESTIAL·REACH-2 포함) |

### 치료 시나리오 임상 데이터

| 요법 | 임상시험 | 중앙 OS | ORR | 비고 |
|------|----------|---------|-----|------|
| 소라페닙 400mg BID | SHARP 2008 (Llovet NEJM) | 10.7개월 | 2% | 최초 승인 전신요법 |
| 렌바티닙 8/12mg QD | REFLECT 2018 (Kudo Lancet) | 13.6개월 | 24% | REFLECT: non-inferior to sora |
| 아테조+베바 q3w | IMbrave150 2020 (Finn NEJM) | 19.2개월 | 30% | 현 1차 표준요법 |
| 레고라페닙 160mg QD | RESORCE 2017 (Bruix Lancet) | 10.6개월 | 11% | 소라페닙 진행 후 2차 |
| 최적 지지요법(BSC) | 역사적 대조군 | ~7개월 | — | 무치료 자연경과 |
---

## 121. 대장암 (Colorectal Cancer · CRC)

> **디렉토리**: [`colorectal-cancer/`](colorectal-cancer/) | **날짜**: 2026-06-23

[![CRC QSP Map](colorectal-cancer/crc_qsp_model.png)](colorectal-cancer/crc_qsp_model.svg)

### 병태생리 요약

대장암은 **Vogelstein 순차 돌연변이 모델**(APC→KRAS→SMAD4→TP53)에 따라 정상 상피에서 선종을 거쳐 침윤암으로 진행됩니다. 약 40%는 KRAS 돌연변이로 항-EGFR 치료에 내성이 있으며, 약 15%는 MMR 결손(dMMR/MSI-H)으로 높은 신항원 부담을 보여 면역관문억제제에 우수한 반응을 나타냅니다. **CMS(Consensus Molecular Subtypes) 1~4**로 분류되며 각 아형은 상이한 치료 전략을 요구합니다.

### 핵심 경로

| 클러스터 | 주요 분자 |
|----------|-----------|
| Wnt/β-catenin | APC, β-catenin, TCF/LEF, RSPO/LGRS5 |
| RAS/MAPK | KRAS, BRAF V600E, MEK1/2, ERK1/2 |
| PI3K/AKT/mTOR | PIK3CA, PTEN, AKT, mTORC1/2, S6K |
| TP53/Apoptosis | TP53 hotspot, MDM2, BAX/BCL2, Casp-3/7 |
| Cell Cycle | CDK4/6, Cyclin D/E, RB1, E2F, p21/p27 |
| TGF-β/EMT | SMAD4, SNAIL/ZEB1/TWIST, E-cadherin↓, MMP-2/9 |
| TME | CD8+ T, Treg, MDSC, M1/M2-TAM, CAF, NK |
| Immune Checkpoints | PD-1/PD-L1, CTLA-4, TIM-3/LAG-3/TIGIT |
| Angiogenesis | HIF-1α, VEGF-A, VEGFR1/2, Ang/Tie2 |
| MSI/MMR | MLH1/MSH2/MSH6/PMS2, MSI-H/TMB-H |

### QSP 모델 구성 (4종 산출물)

| 산출물 | 파일 | 내용 |
|--------|------|------|
| 🗺️ 기계론적 지도 | [`crc_qsp_model.dot/.svg/.png`](colorectal-cancer/) | 130+ 노드, 12개 클러스터 |
| ⚙️ mrgsolve ODE | [`crc_mrgsolve_model.R`](colorectal-cancer/crc_mrgsolve_model.R) | 20구획 ODE, 7치료 시나리오 |
| 📊 Shiny 앱 | [`crc_shiny_app.R`](colorectal-cancer/crc_shiny_app.R) | 6탭 대화형 대시보드 |
| 📚 참고문헌 | [`crc_references.md`](colorectal-cancer/crc_references.md) | 40개 PubMed 인용 |

### 치료 시나리오 임상 데이터

| 요법 | 임상시험 | 주요 결과 | 적응 대상 |
|------|----------|-----------|-----------|
| FOLFOX6 | MOSAIC (André 2004, NEJM) | mOS 17.9개월, adjuvant DFS ↑ | 3기 보조·전이성 |
| FOLFIRI | GERCOR (Tournigand 2004, JCO) | mPFS 8.5개월 | FOLFOX 대안 1차 |
| FOLFOX + Bevacizumab | NO16966 (Saltz 2008, JCO) | mPFS 9.4개월 | 전이성 1차 |
| FOLFIRI + Cetuximab | CRYSTAL (Van Cutsem 2009, NEJM) | RAS-WT ORR 57%, mPFS 9.9개월 | RAS-WT 전이성 |
| FOLFIRI + Bevacizumab | TRIBE (Falcone 2013, JCO) | mPFS 9.7개월 | 전이성 1차 |
| Pembrolizumab | KEYNOTE-177 (André 2020, NEJM) | mPFS 16.5개월, PFS HR 0.60 | MSI-H/dMMR 1차 |
| Regorafenib | CORRECT (Grothey 2013, Lancet) | mOS 6.4개월 vs 5.0개월 | 3차 이상 |

---

## 122. 전립선암 (Prostate Cancer · PC)

> **디렉토리**: [`prostate-cancer/`](prostate-cancer/) | **날짜**: 2026-06-23

[![Prostate Cancer QSP Map](prostate-cancer/pc_qsp_model.png)](prostate-cancer/pc_qsp_model.svg)

### 병태생리 요약

전립선암은 남성에서 가장 흔한 악성종양(전 세계 2위 암 사망원인)으로, **안드로겐 수용체(AR) 신호**가 종양 증식의 핵심 구동력입니다. 시상하부-뇌하수체-고환 축(HPG축)으로 생성된 테스토스테론은 5α-환원효소에 의해 더 강력한 DHT로 변환되어 AR을 활성화하고 PSA 발현·세포 증식을 촉진합니다. 안드로겐 박탈요법(ADT)이 표준 치료이나, 대부분 12~24개월 내에 거세저항성(CRPC)으로 진행하며, 이 시점에서 AR 증폭·돌연변이·스플라이스 변이체(ARv7), PI3K/AKT 활성화(PTEN 소실 ~70%), 골전이 악순환(RANKL-RANK-OPG)이 주요 내성 기전으로 작동합니다.

### 핵심 병태생리 경로

| 클러스터 | 핵심 분자 | 치료 표적 |
|----------|----------|----------|
| **HPG 축** | KiSS1→GnRH→LH→테스토스테론 | GnRH 작용제(루프롤라이드)/길항제(데가렐릭스·렐루골릭스) |
| **안드로겐 생합성** | CYP17A1(17α-hydroxylase/C17,20-lyase), DHT, SRD5A1/2 | 아비라테론(CYP17A1 억제) |
| **AR 신호** | AR-DHT 복합체→핵 전좌→ARE→PSA | 엔잘루타마이드·아팔루타마이드·다로루타마이드 (2세대 ARPI) |
| **AR 내성** | AR 증폭, AR-V7 스플라이스 변이체, F876L/T878A 돌연변이 | 예후 바이오마커(cfDNA·ctRNA) |
| **PI3K/AKT/mTOR** | PTEN 소실(~70%) → PIP3 → AKT → mTORC1/2 | AKT 억제제·mTOR 억제제(연구 중) |
| **RAS/MAPK** | EGFR/HER2 → GRB2/SOS → ERK → ETS 전사인자 | EGFR/MEK 억제제(후기 연구) |
| **세포 주기·아포토시스** | CDK4/6·Cyclin D1·pRb/E2F, BCL-2/BAX, p53, PARP/BRCA1/2 | PARP 억제제(올라파립·루카파립; BRCA2 변이체) |
| **골전이 악순환** | RANKL/RANK/OPG, PTHrP, ET-1(경화성), DKK1, CXCL12/CXCR4 | 데노수맙(항RANKL), 졸레드론산, 라듐-223(α방사선) |
| **종양 면역 미세환경** | PD-1/PD-L1, CTLA-4, TAM-M2, MDSC, Treg | 펨브롤리주맙, 이필리무맙, 시풀루셀-T |
| **임상 지표** | PSA, 테스토스테론, 골스캔(BSI), PSMA PET, rPFS, OS | 진단·모니터링·치료 반응 평가 |

### QSP 모델 구성 (4종 산출물)

| 산출물 | 파일 | 내용 |
|--------|------|------|
| 🗺️ 기계론적 지도 | [`pc_qsp_model.dot/.svg/.png`](prostate-cancer/) | 130+ 노드, 10개 클러스터 (HPG축·생합성·AR·종양·PI3K·MAPK·골전이·약물·면역·임상지표) |
| ⚙️ mrgsolve ODE | [`pc_mrgsolve_model.R`](prostate-cancer/pc_mrgsolve_model.R) | 33구획 ODE — HPG/AR/종양운동/PI3K-AKT/골전이 + 8가지 약물군 PK; 7치료 시나리오 |
| 📊 Shiny 앱 | [`pc_shiny_app.R`](prostate-cancer/pc_shiny_app.R) | 8탭 (환자 프로파일·약물PK·HPG축·AR/PSA·종양·골전이·시나리오 비교·민감도 분석) |
| 📚 참고문헌 | [`pc_references.md`](prostate-cancer/pc_references.md) | 63개 PubMed 인용 (AFFIRM·PREVAIL·COU-AA·CHAARTED·TAX327·PROfound·ALSYMPCA·VISION 포함) |

### 치료 시나리오 임상 데이터

| 요법 | 임상시험 | OS/PFS 혜택 | PSA50 반응 | 적응 단계 |
|------|----------|------------|-----------|---------|
| 루프롤라이드 (ADT) | — | 표준 기준 | ~90% (초기) | mHSPC 1차 |
| 엔잘루타마이드 160 mg QD | AFFIRM (Scher 2012, NEJM) | OS +4개월 | 54% | mCRPC (화학요법 후) |
| 엔잘루타마이드 160 mg QD | PREVAIL (Beer 2014, NEJM) | OS +2.2개월 | 78% | mCRPC (화학요법 전) |
| 아비라테론+프레드니손 | COU-AA-301 (de Bono 2011) | OS +3.9개월 | 29% | mCRPC (화학요법 후) |
| 도세탁셀 75 mg/m² q3w×6 | TAX 327 (Tannock 2004) | OS +3.0개월 | 45% | mCRPC 1차 |
| 올라파립 300 mg BID | PROfound (de Bono 2020) | rPFS HR=0.34 | 33%(BRCA1/2) | mCRPC HRR변이 |
| 라듐-223 55 kBq/kg q4w×6 | ALSYMPCA (Parker 2013) | OS +3.6개월 | — | 증상성 골전이 mCRPC |
| Lu-PSMA-617 7.4 GBq q6w | VISION (Sartor 2021) | OS +4.0개월 | 46% | PSMA+ mCRPC ≥2L |

---

## 123. 급성 골수성 백혈병 (Acute Myeloid Leukemia · AML)

> **디렉토리**: [`acute-myeloid-leukemia/`](acute-myeloid-leukemia/) | **날짜**: 2026-06-23

[![AML QSP Map](acute-myeloid-leukemia/aml_qsp_model.png)](acute-myeloid-leukemia/aml_qsp_model.svg)

### 병태생리 요약

급성 골수성 백혈병(AML)은 조혈 줄기세포 또는 전구세포에서 발생하는 클론성 악성 증식 질환으로, 미성숙 골수 아세포가 골수 내에 축적되어 정상 조혈을 억제합니다. 성인 급성 백혈병 중 가장 흔하며(미국 기준 연 ~2만 명), 중앙 연령은 68세, 5년 생존율은 30% 미만입니다. **FLT3-ITD/TKD**(~30%), **NPM1**(~30%), **DNMT3A**(~20%), **IDH1/2**(~20%), **TP53**(~8%) 등 다수의 재발성 돌연변이가 병태생리를 구동하며, ELN 2022 위험군 분류에 따라 치료 전략이 결정됩니다. 베네토클락스·FLT3 억제제·IDH 억제제 등 표적 치료제의 연이은 승인으로 치료 패러다임이 급격히 변화하고 있습니다.

### 핵심 병태생리 경로

| 클러스터 | 핵심 분자 | 치료 표적 |
|----------|----------|----------|
| **정상 조혈** | HSC→CMP/GMP/MEP→성숙 혈구, SCF/TPO/EPO/G-CSF | 골수억제 모니터링 (Friberg 모델) |
| **AML 분자 병인** | FLT3-ITD→RAS/RAF→MAPK 증식; FLT3-TKD→JAK/STAT5→BCL-2 항아포프토시스 | 길테리티닙(FLT3/AXL), 퀴자르티닙(FLT3-ITD), 미도스타우린 |
| **PI3K/AKT/mTOR** | FLT3/KIT→PI3K→AKT→FOXO 비활성화→생존 | 연구 중인 AKT·mTOR 억제제 |
| **BCL-2 패밀리** | BCL-2·BCL-xL·MCL-1 대 BAX·BAK·BIM·PUMA·NOXA | 베네토클락스(BH3 모방체·BCL-2 선택적) |
| **후성유전학** | IDH1/2→2-HG→TET2 억제→과메틸화; DNMT3A 소실 | 에나시데닙(IDH2), 이보시데닙(IDH1), 아자시티딘(DNMT 억제) |
| **골수 미세환경** | CXCL12-CXCR4(LSC 유지), VLA-4-피브로넥틴(CAMDR), HIF-1α(저산소 적응) | 플레릭사포(CXCR4), 키잘리(Hedgehog) |
| **면역 회피** | CD47("don't eat me")→SIRPα; PD-L1→PD-1 T세포 고갈 | 마그롤리맙(CD47), 펨브롤리주맙 |
| **임상 지표** | 골수 아세포 %, MRD(NPM1 PCR·FLT3 VAF), ANC/PLT, TLS 위험 | CR<5%, MRD음성, 골수 억제 모니터링 |

### QSP 모델 구성 (4종 산출물)

| 산출물 | 파일 | 내용 |
|--------|------|------|
| 🗺️ 기계론적 지도 | [`aml_qsp_model.dot/.svg/.png`](acute-myeloid-leukemia/) | 281 노드, 10개 클러스터, 211 엣지 |
| ⚙️ mrgsolve ODE | [`aml_mrgsolve_model.R`](acute-myeloid-leukemia/aml_mrgsolve_model.R) | 21구획 ODE — 9 약물 PK + 12 질환 PD, Friberg 골수억제, 7치료 시나리오 |
| 📊 Shiny 앱 | [`aml_shiny_app.R`](acute-myeloid-leukemia/aml_shiny_app.R) | 6탭 (환자 프로파일·약물PK·백혈병 역학·임상 엔드포인트·시나리오 비교·바이오마커) |
| 📚 참고문헌 | [`aml_references.md`](acute-myeloid-leukemia/aml_references.md) | 38개 PubMed 인용 (VIALE-A·ADMIRAL·QuANTUM-R·RATIFY·APL ATO 포함) |

### 치료 시나리오 임상 데이터

| 요법 | 임상시험 | CR 율 | 중앙 OS | 적응 단계 |
|------|----------|-------|---------|---------|
| 7+3 (Ara-C+Ida) | 표준 기준 | ~65–75% | ~1–2년 | 신규 집중 요법 가능 |
| 베네토클락스+아자시티딘 | VIALE-A (DiNardo 2020) | 66.4% | 14.7개월 | 신규 집중 요법 불가 |
| 길테리티닙 120 mg/일 | ADMIRAL (Perl 2019) | 21% | 9.3개월 | 재발/불응 FLT3+ |
| 퀴자르티닙 | QuANTUM-R (Cortes 2019) | 4.3% CRc | 6.2개월 | 재발/불응 FLT3-ITD+ |
| 에나시데닙 100 mg/일 | Stein 2017 | 19.3% | 9.3개월 | 재발/불응 IDH2+ |
| ATRA + 삼산화비소 | APL0406 (Lo-Coco 2013) | ~95% | >90% 2년 | APL (PML-RARA) |

---

## 124. 위선암 (Gastric Cancer · GC)

> **디렉토리**: [`gastric-cancer/`](gastric-cancer/) | **날짜**: 2026-06-23

[![GC QSP Map](gastric-cancer/gc_qsp_model.png)](gastric-cancer/gc_qsp_model.svg)

### 병태생리 요약

위선암(Gastric Adenocarcinoma)은 전 세계 5번째로 흔한 암이자 암 사망 원인 3위로, 연간 약 109만 명의 신규 환자와 77만 명의 사망이 발생합니다(GLOBOCAN 2020). 한국·일본·중국 등 동아시아에서 발생률이 가장 높으며, **헬리코박터 파일로리(H. pylori)** 감염이 가장 중요한 위험 인자입니다. Correa 폭포(정상 위 점막 → 만성 위염 → 위축성 위염 → 장상피화생 → 이형성증 → 위선암)가 장형 위암의 발생기전을 설명합니다. **TCGA 분자 아형**은 EBV 양성(9%), MSI-H(22%), 유전체 안정(GS, 20%), 염색체 불안정(CIN, 50%)으로 분류되며, 각 아형별 치료 전략이 상이합니다.

### 핵심 병태생리 경로

| 클러스터 | 핵심 분자 | 치료 표적 |
|----------|----------|----------|
| **H. pylori 염증** | CagA→NF-κB/STAT3/PI3K; VacA→ROS; IL-8/IL-6/TNF-α | H. pylori 제균(아목시실린+클래리스로마이신+PPI) |
| **수용체 티로신키나제** | HER2(IHC3+, ~15%)·FGFR2 amp·MET amp·VEGFR2 | 트라스투주맙·라무시루맙·아파티닙·졸베툭시맙 |
| **신호전달** | PI3K/AKT/mTOR·RAS/RAF/MEK/ERK·JAK1/2/STAT3 | mTOR 억제제(연구 중), MEK 억제제 |
| **세포주기/아포프토시스** | TP53 mut·CDK4/6·BCL-2·BAX·Caspase-3/9 | CDK4/6i(연구 중), 베네토클락스(연구 중) |
| **후성유전학·EMT** | EBV 과메틸화·MLH1 소실·EZH2·TGF-β/SMAD·SNAIL/ZEB1 | EZH2i(연구 중), HDAC 억제제 |
| **종양 미세환경(TME)** | CAF·TAM M2·MDSC·Treg·IL-10·TGF-β | TME 재형성 병용요법 |
| **면역관문** | PD-L1/PD-1(CPS ≥5)·CTLA-4·LAG-3·TIM-3·CD47 | 니볼루맙·펨브롤리주맙(CPS ≥5/MSI-H) |
| **혈관신생** | VEGF-A→VEGFR2·HIF-1α·Ang-2/Tie2 | 라무시루맙(VEGFR2 mAb) |

### QSP 모델 구성 (4종 산출물)

| 산출물 | 파일 | 내용 |
|--------|------|------|
| 🗺️ 기계론적 지도 | [`gc_qsp_model.dot/.svg/.png`](gastric-cancer/) | 212 노드, 10개 클러스터, 250 엣지 |
| ⚙️ mrgsolve ODE | [`gc_mrgsolve_model.R`](gastric-cancer/gc_mrgsolve_model.R) | 28구획 ODE (12 약물 PK + 16 질환 PD), Simeoni TGI, TMDD, 6치료 시나리오 |
| 📊 Shiny 앱 | [`gc_shiny_app.R`](gastric-cancer/gc_shiny_app.R) | 6탭 (환자 프로파일·약물PK·종양 동태·임상 엔드포인트·시나리오 비교·바이오마커) |
| 📚 참고문헌 | [`gc_references.md`](gastric-cancer/gc_references.md) | 60개 PubMed 인용 (ToGA·RAINBOW·CheckMate649·SPOTLIGHT·DESTINY-Gastric01·FLOT4·KEYNOTE-811 포함) |

### 치료 시나리오 임상 데이터

| 요법 | 임상시험 | ORR | 중앙 OS | 적응 단계 |
|------|----------|-----|---------|---------|
| FLOT × 4 pre/post-op | FLOT4 (Al-Batran 2019) | — | 50개월 | 절제 가능 국소 진행성 |
| 트라스투주맙 + FOLFOX/XELOX | ToGA (Bang 2010) | 47% | 13.8개월 | HER2+ 1차 |
| 펨브롤리주맙 + 트라스투주맙 + 화학요법 | KEYNOTE-811 (2024) | 74.4% | 20.0개월 | HER2+ CPS≥1 1차 |
| 라무시루맙 + 파클리탁셀 | RAINBOW (Wilke 2014) | 28% | 9.6개월 | 2차 |
| 니볼루맙 + FOLFOX/XELOX | CheckMate 649 (CPS≥5) | 60% | 14.4개월 | CPS≥5 1차 |
| T-DXd 6.4 mg/kg Q3W | DESTINY-Gastric01 | 51.3% | 12.5개월 | HER2+ 2차 |
| 졸베툭시맙 + mFOLFOX6 | SPOTLIGHT (Shitara 2023) | — | 18.2개월 | CLDN18.2+ HER2- 1차 |

---

## 125. 골수섬유증 (Myelofibrosis · MF)

> **디렉토리**: [`myelofibrosis/`](myelofibrosis/) | **날짜**: 2026-06-23

[![MF QSP Map](myelofibrosis/mf_qsp_model.png)](myelofibrosis/mf_qsp_model.svg)

### 병태생리 요약

골수섬유증(Myelofibrosis, MF)은 조혈 줄기세포(HSC)에서 발생하는 만성 골수증식성 종양(MPN)으로, 골수 내 레티큘린·콜라겐 섬유화, 비정상적 거핵구 증식, 비장·간을 중심으로 하는 수외 조혈(EMH)을 특징으로 합니다. 일차성(PMF)과 다혈증성·본태성 혈소판증가증 후 이행성(Post-PV/ET MF)으로 구분되며, 미국 기준 연 약 3,000~4,000명이 발생하고 중앙 연령은 67세입니다. **JAK2 V617F**(~55–60%), **CALR exon 9 mutations**(~25–30%), **MPL W515L/K**(~5–10%)이 주요 드라이버 돌연변이이며, 모두 JAK/STAT 신호 경로를 구성적으로 활성화합니다. **ASXL1, EZH2, SRSF2, IDH1/2, TP53** 등 고위험 분자 이상(HMR)이 동반되면 AML 전환 위험이 현저히 증가합니다. JAK 억제제(루소리티닙, 페드라티닙, 파크리티닙, 모멜로티닙)와 BET 억제제 병용(pelabresib + ruxolitinib, MANIFEST-2)이 현재 치료 패러다임을 형성하고 있습니다.

### 핵심 병태생리 경로

| 클러스터 | 핵심 분자 | 치료 표적 |
|----------|----------|----------|
| **분자 드라이버** | JAK2V617F·CALR T1/T2·MPL W515L/K·ASXL1·EZH2·SRSF2·IDH1/2 | DIPSS-Plus/MIPSS70 위험도 분류 |
| **JAK/STAT 신호** | JAK1/2→STAT3/5 인산화→BCL-2/MCL-1/MYC/PIM1 표적 유전자 | 루소리티닙(JAK1/2), 페드라티닙(JAK2/FLT3), 파크리티닙(JAK2/FLT3/IRAK1/CSF1R) |
| **HSC 니치** | LT-HSC→ST-HSC→CMP/GMP/MEP; CXCL12/CXCR4, SCF/c-KIT, TPO | CXCR4 축, c-KIT 신호 |
| **골수 섬유화** | 비정상 거핵구→TGF-β1/PDGF/bFGF/CTGF 분비→MSC 활성화→콜라겐 과침착 | TGF-β 억제제(연구 중), 갈루니서팁 |
| **사이토카인 폭풍** | IL-1β/6/8/10, TNF-α, IFN-γ, CXCL10, NF-κB, mTOR | BET 억제제(펠라브레십): BRD4→c-MYC/IL-6/TNF-α 억제 |
| **수외 조혈** | 순환 CD34+→비장/간 EMH; 비장 종대(비장 부피 >450 mL) | 비장용적 축소(SVR35) 치료 목표 |
| **혈액학적 결과** | 빈혈(Hgb <10 g/dL), 혈소판감소증, 순환 아세포, 수혈 의존성 | 모멜로티닙(ACVR1 억제→빈혈 개선), EPO 제제 |
| **혈전·출혈** | 혈소판 활성화, 트롬빈 생성, PAI-1, 내피 기능 이상 | 혈전 예방(저용량 아스피린), JAK 억제제 |
| **임상 엔드포인트** | SVR35(비장 35% 축소), TSS50(증상 50% 개선), 골수 조직 반응, AML 전환, OS | IWG-MRT 2013 / ELN 2023 반응 기준 |

### QSP 모델 구성 (4종 산출물)

| 산출물 | 파일 | 내용 |
|--------|------|------|
| 🗺️ 기계론적 지도 | [`mf_qsp_model.dot/.svg/.png`](myelofibrosis/) | **204 노드, 12 클러스터, 224 엣지** |
| ⚙️ mrgsolve ODE | [`mf_mrgsolve_model.R`](myelofibrosis/mf_mrgsolve_model.R) | **23구획 ODE** (Ruxolitinib 2구획 PK + Fedratinib/Pacritinib/Pelabresib PK + JAK/STAT PD + 적혈·거핵계 + 비장·섬유화 + 사이토카인 + TSS), **6 치료 시나리오** (무치료·Rux 20mg·Rux 15mg·Fedratinib·Rux+Pelabresib·Pacritinib), COMFORT-I·JAKARTA·PERSIST-2·MANIFEST-2 임상 보정 |
| 📊 Shiny 앱 | [`mf_shiny_app.R`](myelofibrosis/mf_shiny_app.R) | **6탭** (환자 프로파일·DIPSS Plus·PK 프로파일·PD 바이오마커·임상 엔드포인트·치료 비교·바이오마커 역학), plotly 인터랙티브, bslib 다크 테마 |
| 📚 참고문헌 | [`mf_references.md`](myelofibrosis/mf_references.md) | **36개 PubMed 인용** (COMFORT-I/II·JAKARTA·PERSIST-2·SIMPLIFY-1·MANIFEST-2·JAK2V617F 발견·CALR 돌연변이·QSP 모델링 포함) |

### 치료 시나리오 임상 데이터

| 요법 | 임상시험 | SVR35 | TSS50 | 중앙 OS | PMID |
|------|----------|-------|-------|---------|------|
| 루소리티닙 20mg BID | COMFORT-I (Verstovsek 2012) | 41.9% vs 0.7% | 45.9% vs 5.3% | NR vs 27.4개월 | [22375971](https://pubmed.ncbi.nlm.nih.gov/22375971/) |
| 루소리티닙 vs BAT | COMFORT-II (Harrison 2012) | 28.5% vs 0% | — | — | [22375970](https://pubmed.ncbi.nlm.nih.gov/22375970/) |
| 페드라티닙 400mg QD | JAKARTA (Pardanani 2015) | 36% vs 1% | 36% vs 6% | — | [26003172](https://pubmed.ncbi.nlm.nih.gov/26003172/) |
| 파크리티닙 200mg BID | PERSIST-2 (Mesa 2017) | 18% vs 3% | — | — | [29049469](https://pubmed.ncbi.nlm.nih.gov/29049469/) |
| 모멜로티닙 200mg QD | SIMPLIFY-1 (Mesa 2017) | 26.5% vs 29% | 28% vs 42% | — | [28930484](https://pubmed.ncbi.nlm.nih.gov/28930484/) |
| Pelabresib + Rux | MANIFEST-2 (Pemmaraju 2024) | **65.9% vs 35.2%** | 52.3% vs 37.5% | NR | [39504566](https://pubmed.ncbi.nlm.nih.gov/39504566/) |

---

## 126. 췌장 선암 (Pancreatic Ductal Adenocarcinoma · PDAC)

> **디렉토리**: [`pancreatic-cancer/`](pancreatic-cancer/) | **날짜**: 2026-06-23

[![PDAC QSP Map](pancreatic-cancer/pdac_qsp_model.png)](pancreatic-cancer/pdac_qsp_model.svg)

### 병태생리 요약

췌장 선암(PDAC)은 고형 악성종양 중 가장 치명적인 암으로, 전 병기 합산 5년 생존율이 약 12%에 불과합니다. 전 세계적으로 연간 약 50만 명이 진단받으며, 암 사망 원인 7위를 차지합니다. 대부분(~80%)이 절제 불가능한 국소 진행성 또는 전이성 상태로 진단됩니다. **KRAS 돌연변이**(~95%, G12D 44%·G12V 26%·G12R 14%)가 핵심 발암 구동인자이며, **TP53 소실**(~75%), **CDKN2A/p16 소실**(~90%), **SMAD4 소실**(~55%)이 동반됩니다. 종양 부피의 80~90%를 차지하는 **치밀 섬유화 기질(desmoplastic stroma)**이 약물 침투를 차단하고 면역억제 미세환경을 형성하는 독특한 병태생리를 보입니다.

### 핵심 병태생리 경로

| 클러스터 | 핵심 분자 | 치료 표적 |
|----------|----------|----------|
| **KRAS/MAPK 신호** | KRAS G12D→BRAF→MEK→ERK→MYC·ETS | MRTX1133(KRAS G12D), Adagrasib(G12C), SOS1 억제제(BI-3406) |
| **PI3K/AKT/mTOR** | PI3K→PIP3→AKT→mTORC1/2→S6K·4EBP1 | 에베롤리무스, 코판리십(연구 중) |
| **TGF-β/SMAD/EMT** | TGFβ→SMAD2/3→SMAD4→SNAIL/ZEB1→E-Cad소실·MMP2/9→침윤·전이 | 갈루니서팁(TGFβRI 억제), SMAD4 biomarker |
| **섬유화 기질** | PSC→myoCAF/iCAF→콜라겐·히알루론산→IFP↑→약물침투↓ | PEGPH20(히알루로니다제), 피르페니돈, FAP CAR-T |
| **면역억제 TME** | PD-L1·CTLA-4·LAG-3·TIM-3, Treg·MDSC·TAM-M2, IL-10·TGF-β | 펨브롤리주맙(MSI-H), 이필리무맙, 렐라틀리맙 |
| **혈관신생** | HIF-1α→VEGF-A→VEGFR2·Ang-2→신생혈관(미성숙) | 베바시주맙(제한적 효과), 라무시루맙 |
| **DNA 손상·HRD** | BRCA1/2·PALB2·ATM 돌연변이→HRD→PARP 의존성 | 올라파립(POLO 승인), 루카파립, 니라파립 |
| **임상 지표** | CA 19-9·CEA·ctDNA(KRAS VAF)·CA 19-9 응답, RECIST, PFS/OS | 치료 반응 조기 예측·내성 감시 |

### QSP 모델 구성 (4종 산출물)

| 산출물 | 파일 | 내용 |
|--------|------|------|
| 🗺️ 기계론적 지도 | [`pdac_qsp_model.dot/.svg/.png`](pancreatic-cancer/) | **264노드, 10개 클러스터** (정상 췌장·드라이버 돌연변이·KRAS 신호·TGF-β/EMT·섬유화 기질·면역 TME·혈관신생·약물 PK·종양 생물학·임상 지표) |
| ⚙️ mrgsolve ODE | [`pdac_mrgsolve_model.R`](pancreatic-cancer/pdac_mrgsolve_model.R) | 젬시타빈(2구획)·nab-파클리탁셀(2구획)·옥살리플라틴·이리노테칸/SN-38·5-FU·MRTX1133·올라파립 PK + Simeoni TGI(4-전이구획) + CA19-9 + Friberg 골수억제 + 기질저항; **7치료 시나리오**(무치료·Gem mono·Gem+nab-Pac·FOLFIRINOX·mFOLFIRINOX·MRTX1133·올라파립), MPACT·PRODIGE4·POLO 임상 보정 |
| 📊 Shiny 앱 | [`pdac_shiny_app.R`](pancreatic-cancer/pdac_shiny_app.R) | **6탭** (환자 프로파일·약물 PK·종양 역학·바이오마커·임상 엔드포인트·시나리오 비교), plotly 인터랙티브, bslib darkly 테마, 1,357줄 |
| 📚 참고문헌 | [`pdac_references.md`](pancreatic-cancer/pdac_references.md) | **51개 PubMed 인용** (MPACT·PRODIGE4/ACCORD11·POLO·CodeBreak·MRTX1133 Phase I·PEGPH20·ICI in PDAC 포함) |

### 치료 시나리오 임상 데이터

| 요법 | 임상시험 | ORR | 중앙 PFS | 중앙 OS | 적응 |
|------|----------|-----|---------|---------|------|
| 젬시타빈 mono | Burris 1997 (NEJM) | 7% | 3.7개월 | 6.7개월 | 1차 표준 기준 |
| Gem + nab-Paclitaxel | MPACT (Von Hoff 2013, NEJM) | 23% | 5.5개월 | **8.5개월** | 1차 전이성 |
| FOLFIRINOX | PRODIGE4/ACCORD11 (Conroy 2011) | 31.6% | 6.4개월 | **11.1개월** | 1차 PS 0–1 |
| mFOLFIRINOX | PRODIGE24 (Conroy 2018, adjuvant) | 28% | 6.0개월 | 10.5개월 | 1차 또는 수술 후 보조 |
| MRTX1133 (KRAS G12D) | Phase I/II (NCT05737706, 진행 중) | ~40% (예측) | ~4개월 | ~8개월 | G12D+ 2차 이상 |
| 올라파립 유지 | POLO (Golan 2019, NEJM) | — | **7.4 vs 3.8개월** (HR 0.53) | NS (OS) | gBRCA+ 1차 후 유지 |

---

## 127. 골수이형성 증후군 (Myelodysplastic Syndrome · MDS)

> **디렉토리**: [`myelodysplastic-syndrome/`](myelodysplastic-syndrome/) | **날짜**: 2026-06-23

[![MDS QSP Map](myelodysplastic-syndrome/mds_qsp_model.png)](myelodysplastic-syndrome/mds_qsp_model.svg)

### 병태생리 요약

골수이형성 증후군(MDS)은 조혈 줄기세포(HSC)의 클론성 이상 증식으로 인해 골수에서 하나 이상의 혈구 계열에 이형성(dysplasia)이 나타나고 비효율적 조혈(ineffective hematopoiesis)이 발생하는 혈액 종양입니다. 서구 기준 연 발생률은 10만 명당 약 4–5명이며, 65세 이상에서 발생률이 급격히 높아져 중앙 연령은 약 70세입니다. 전체 MDS의 약 30%는 급성 골수성 백혈병(AML)으로 진행합니다.

**핵심 병태생리:** (1) CHIP(Clonal Hematopoiesis of Indeterminate Potential) 상태의 초기 돌연변이(TET2, DNMT3A, ASXL1)가 클론성 이점을 제공 → (2) 스플라이싱 인자 돌연변이(SF3B1 고리 철아세포·SRSF2 CMML·U2AF1·ZRSR2), 후성유전 이상(IDH1/2 → 2-HG → TET2 억제, EZH2 H3K27me3↑), 전사인자 돌연변이(RUNX1, TP53), del(5q)·단염색체 7번 등 세포유전학 이상이 중첩되어 조혈 분화가 왜곡 → (3) GDF11/TGF-β1 상승으로 적혈구 성숙이 차단(Smad2/3 과활성화)되어 비효율 적혈구생성 → (4) 헵시딘↑/페로포틴↓에 의한 철 이용 장애 + 수혈 의존성에 의한 이차 철과부하 → (5) 사이토카인 과활성(TNF-α, IFN-γ, IL-6)에 의한 정상 HSC 억제.

### 핵심 병태생리 경로

| 클러스터 | 핵심 분자 | 치료 표적 |
|----------|----------|----------|
| **스플라이싱 인자** | SF3B1(~30% MDS-RS)·SRSF2·U2AF1·ZRSR2 → 비정상 인트론 스플라이싱·mRNA 이형체·NMD | SF3B1 mut → ABCB7 소실 → 고리 철아세포; H3B-8800(스플라이싱 조절제, 연구 중) |
| **후성유전학** | TET2(5mC→5hmC)·DNMT3A·IDH1/2(2-HG)·ASXL1/EZH2(H3K27me3↑) | 아자시티딘·데시타빈(DNMT1 포획→탈메틸화); 이보시데닙(IDH1), 에나시데닙(IDH2) |
| **클론 진화** | CHIP→MDS clone: TP53(복잡핵형)·RUNX1·NPM1/FLT3(AML 전환) | IPSS-R/IPSS-M 위험 분층; allo-SCT |
| **GDF11/TGF-β** | MDS 골수 GDF11/ActB/GDF8 상승 → Smad2/3↑ → 적혈구 성숙 차단 | **루스파터셉트**(ActRIIB-Fc TGF-β 리간드 포획 → Smad2/3↓ → 후기 적혈구 성숙 회복) |
| **del(5q) 기전** | RPS14 반수부족→리보솜 스트레스→p53 활성화; miR-145/146a 소실→TIRAP·TRAF6↑→NF-κB·IL-6 | **레날리도마이드**(CRBN→CK1α 분해→MDM2 억제→p53 재활성화 + RPS14 복원) |
| **BM 미세환경** | CXCL12/CXCR4(HSC 유착)·TNF-α/IFN-γ(정상 HSC 억제)·MDSC·Treg·NK세포 기능 저하 | CXCR4 길항제(플레릭사포): 클론 동원 |
| **철 항상성** | 헵시딘↑(IL-6→STAT3→BMP/SMAD)→페로포틴↓→철 이용 감소; 수혈의존성→이차 철과부하 | 데페록사민·데페라시록스(철 킬레이트); 다르베포에틴(헵시딘 억제) |
| **임상 지표** | Hgb <10 g/dL·PLT <50K·ANC <0.5K·BM 아세포%·IPSS-R(Very Low~Very High)·IPSS-M(분자) | 수혈 독립성(TI), HI-E/HI-P/HI-N, CR/mCR, AML 전환 |

### QSP 모델 구성 (4종 산출물)

| 산출물 | 파일 | 내용 |
|--------|------|------|
| 🗺️ 기계론적 지도 | [`mds_qsp_model.dot/.svg/.png`](myelodysplastic-syndrome/) | **324 노드, 10개 클러스터** (정상 조혈·클론 진화·후성유전학·스플라이싱·BM 미세환경·사이토카인·이형성·임상 지표·약물 PK/PD·AlloSCT) |
| ⚙️ mrgsolve ODE | [`mds_mrgsolve_model.R`](myelodysplastic-syndrome/mds_mrgsolve_model.R) | **18구획 ODE** (6 약물 PK: AZA 2구획·DEC·LEN·LUSP·Darbe + 12 질환 PD: Blast·VAF·DNAmeth·IneffErythro·EryProg·Hgb·PLT·ANC·GDF11·Hepcidin·IronStore·TGFb_signal), **7치료 시나리오** (BSC·AZA SC·DEC IV·Oral-DEC/CED·Lenalidomide·Luspatercept·VEN+AZA), COMMANDS·MEDALIST·MDS-003·ASTX727·VIALE-A 임상 보정 |
| 📊 Shiny 앱 | [`mds_shiny_app.R`](myelodysplastic-syndrome/mds_shiny_app.R) | **6탭** (환자 프로파일·IPSS-R 계산기·Drug PK·Disease Dynamics·Hematologic Endpoints·Scenario Comparison·Biomarker Tracker), plotly 인터랙티브, bslib darkly 테마, 1,074줄 |
| 📚 참고문헌 | [`mds_references.md`](myelodysplastic-syndrome/mds_references.md) | **30+ PubMed 인용** (COMMANDS·MEDALIST·MDS-003·MDS-004·ASTX727·QUAZAR·VIALE-A·SF3B1/TET2 발견·IPSS-R·IPSS-M 포함) |

### 치료 시나리오 임상 데이터

| 요법 | 임상시험 | 주요 결과 | 적응 단계 |
|------|----------|----------|---------|
| 아자시티딘 75 mg/m² SC d1-7 q28d | CALGB 9221 (Silverman 2002) | CR 7%·PR 16%·HI 37%; OS 혜택(고위험 MDS) | 중등도~고위험 MDS 1차 |
| 아자시티딘 vs CCR | AZA-001 (Fenaux 2009, Lancet Oncol) | OS **24.5 vs 15.0개월** (HR 0.58) | IPSS High/Int-2 |
| 데시타빈 20 mg/m² IV d1-5 q28d | D-0007 (Kantarjian 2006) | ORR 17%; MRD CR 9% | 중등도~고위험 MDS |
| Oral-DEC/cedazuridine (ASTX727) | ASCERTAIN (Garcia-Manero 2020, JCO) | 5d AUC 노출 IV 동등; ORR 25% | 중등도~고위험 MDS |
| 레날리도마이드 10 mg QD 21d q28d | MDS-003 (List 2006, NEJM) | TI 67%·세포유전학 CR 45% | del(5q) 저위험 MDS |
| 루스파터셉트 1.0 mg/kg q3w | COMMANDS (Platzbecker 2023, NEJM) | **TI ≥12주 58.5% vs 31.1%** (vs ESA); Hgb ↑1.5 g/dL 이상 | 저~중등도위험 MDS-RS (1차) |
| 루스파터셉트 1.0 mg/kg q3w | MEDALIST (Fenaux 2020, NEJM) | TI ≥8주 **38% vs 13%** (vs 위약); HR-QoL 개선 | ESA 실패 MDS-RS 수혈의존성 |
| 베네토클락스 + 아자시티딘 | VIALE-A 하위군/고위험 MDS 외삽 | CR ~15-20% (MDS 추정); MDS TP53 반응 제한적 | 고위험/집중 요법 불가 MDS |
