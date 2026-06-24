# QSP 질환 모델 라이브러리 (QSP Disease Model Library)

> 매일 **Claude Code Routine(CCR)** 이 질환 하나를 선택해 **정량적 시스템 약리학(Quantitative Systems Pharmacology, QSP)** 모델을 처음부터 끝까지 구축하고 `main`에 직접 커밋하는, **살아 있는(living) 오픈 모델 라이브러리**입니다.

![models](https://img.shields.io/badge/models-136-blue) ![framework](https://img.shields.io/badge/QSP-mrgsolve%20%C2%B7%20Shiny%20%C2%B7%20Graphviz-success) ![automation](https://img.shields.io/badge/built%20by-Claude%20Code%20Routine-orange)

현재 **136개 질환**에 대한 완성된 QSP 모델이 수록되어 있으며, 각 모델은 ①기계론적 지도, ②mrgsolve ODE 모델, ③Shiny 대시보드, ④참고문헌의 네 가지 산출물로 구성됩니다. 아래 [모델 갤러리](#-모델-갤러리-model-gallery)에서 전체 목록을 확인할 수 있습니다.

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
| 136 | 중환자·패혈증 | [**패혈증 / 전신 염증 반응 증후군**<br><sub>Sepsis / SIRS · SEP</sub>](sepsis/) | <a href="sepsis/sep_qsp_model.svg"><img src="sepsis/sep_qsp_model.png" width="190" alt="Sepsis"></a> | LPS/LTA→TLR4/TLR2→MyD88→IRAK1/4→TRAF6→IKK→NF-κB 활성화 → TNF-α·IL-1β·IL-6·IL-8 사이토카인 폭풍 → 보체 활성화(C3a·C5a) → TF 발현·응고 항진(Thrombin→Fibrin·DIC) → iNOS↑·NO↑·혈관 확장→MAP 감소(혈관 마비) → 내피세포 기능장애·글리코칼릭스 탈락→모세혈관 누출 → HPA 축 활성화(ACTH→Cortisol·CIRCI) → 다장기부전(심기능 저하·AKI·ARDS·간기능 장애·장 장벽 손상). 피페라실린/타조박탐(2구획 IV PK; %fT>MIC 살균)·노르에피네프린(Emax=30mmHg EC50=0.15)·하이드로코르티손(200mg/day; Imax=65% 사이토카인 억제)·바소프레신(V1a 수용체, 0.03 units/min) 4종 PK/PD. **124+ 노드 11클러스터**, **20구획 ODE**(세균·호중구·대식세포·TNF/IL-6/IL-10/IL-1β·트롬빈·피브린·혈소판·NO·장기손상·젖산·MAP·크레아티닌·항생제·승압제·스테로이드), **6치료 시나리오**(S1 무치료·S2 조기항생제·S3 항생제+NE·S4 풀번들·S5 지연항생제·S6 불응성패혈성쇼크). Rivers 2001(EGDT)·Kumar 2006(항생제 타이밍)·VASST 2008·ADRENAL 2018 임상시험 보정. SOFA 점수·젖산·프로칼시토닌·CRP·혈소판·INR·크레아티닌·MAP 바이오마커. 55개 PubMed 인용.<br>[🗺️ 지도](sepsis/sep_qsp_model.svg) · [⚙️ mrgsolve](sepsis/sep_mrgsolve_model.R) · [📊 Shiny](sepsis/sep_shiny_app.R) · [📚 문헌](sepsis/sep_references.md) · [📄 README](sepsis/README.md) |


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

---

## 128. 섬유근통 (Fibromyalgia · FM)

> **디렉토리**: [`fibromyalgia/`](fibromyalgia/) | **날짜**: 2026-06-23

[![FM QSP Map](fibromyalgia/fm_qsp_model.png)](fibromyalgia/fm_qsp_model.svg)

### 병태생리 요약

섬유근통(FM)은 전신 만성 통증과 피로, 수면 장애, 인지 기능 저하를 특징으로 하는 복합적 중추신경계 통증 장애입니다. 유병률은 약 2–5%로 여성에서 더 흔하며(여:남 = 7:1), 2016년 ACR 기준으로 진단합니다(WPI ≥ 7 + SS ≥ 5, 또는 WPI 4–6 + SS ≥ 9). 섬유근통의 핵심 기전은 **중추감작(central sensitization)** — 척수 후각에서 정상 자극에 대해 비정상적으로 증폭된 통증 신호 처리입니다.

**핵심 병태생리:** (1) 말초 조직의 미세손상 및 NGF/PGE2 상승으로 C-섬유/Aδ 구심 활동 증가 → (2) 척수 후각 WDR 뉴런에서 substance P/CGRP 분비 증가 → NMDA 수용체 활성화 및 Mg²⁺ 차단 해제 → wind-up(시간적 합산) → 척수 장기강화(LTP)로 중추감작 확립 → (3) 하행성 통증조절계(PAG-RVM-LC) 기능 저하(NE↓·5-HT↓)로 내인성 통증억제 소실(DNIC 결함) → (4) 척수 미세아교세포 활성화(IL-1β/TNF-α/BDNF 분비) → NMDA 수용체 인산화 강화·KCC2↓→GABA 탈억제 → (5) HPA 축 이상(저코르티솔혈증 → 항염증 제동 소실) + GH/IGF-1 결핍(조직 회복 저하) → (6) 교감신경 과항진·HRV↓(부교감 저하) → (7) α파-δ파 수면 침범으로 서파수면(SWS) 결손 → 수면-통증 악순환.

### 핵심 병태생리 경로

| 클러스터 | 핵심 분자 | 치료 표적 |
|----------|----------|----------|
| **말초 감작** | NGF↑·PGE2·TRPV1·NaV1.7/1.8·DRG C/Aδ 섬유 | 말초 TRPV1 길항제(연구 중) |
| **척수 중추감작** | SP(CSF↑3배)·Glu·NMDA·AMPA·NK1R·PKCε·ERK1/2·wind-up·LTP | 프레가발린(α2δ-1 차단→Ca²+↓→Glu/SP 방출↓) |
| **하행성 통증조절** | PAG-RVM axis·LC(NE)·Raphe(5-HT)·DPMS·DNIC 소실 | 둘록세틴/밀나시프란(SERT+NET 억제→NE/5-HT↑→하행억제↑) |
| **신경염증** | 척수 미세아교세포·BDNF/TrkB·IL-1β·TNF-α·CX3CR1/fractalkine·NLRP3·KCC2↓ | 미세아교세포 조절 연구(저용량 날트렉손 등) |
| **HPA 축** | CRH-ACTH-코르티솔 음성피드백 약화·저코르티솔혈증·GH/IGF-1 결핍 | 운동 요법(HPA 정상화 효과) |
| **자율신경** | 교감신경 과항진·HRV↓·기립성 이상·발한 기능 이상 | 자율신경 훈련; HRV 바이오피드백 |
| **수면** | 아데노신·서파수면(SWS)·α-δ 침범·비회복성 수면·수면GH 분비↓ | 아미트립틸린 25 mg QHS(H1 차단→진정→SWS↑) |
| **임상 지표** | Pain NRS·FIQR(0–100)·Fatigue VAS·PHQ-9·GAD-7·PPT·PGIC·SF-36 | 30% 통증 감소 = 임상적 의미 있는 반응 |

### QSP 모델 구성 (4종 산출물)

| 산출물 | 파일 | 내용 |
|--------|------|------|
| 🗺️ 기계론적 지도 | [`fm_qsp_model.dot/.svg/.png`](fibromyalgia/) | **100+ 노드, 10개 클러스터** (말초감작·척수후각·뇌통증회로·하행조절·HPA 축·자율신경계·수면·신경염증·약물 PK/PD·임상 지표) |
| ⚙️ mrgsolve ODE | [`fm_mrgsolve_model.R`](fibromyalgia/fm_mrgsolve_model.R) | **30구획 ODE** (4 약물 PK: 둘록세틴 2구획·프레가발린·밀나시프란·아미트립틸린 + 26 질환 PD), **6치료 시나리오** (미치료·DUL·PRE·MIL·DUL+PRE 병용·TCA) |
| 📊 Shiny 앱 | [`fm_shiny_app.R`](fibromyalgia/fm_shiny_app.R) | **6탭** (환자 프로파일·Drug PK·PD 핵심지표·임상 엔드포인트·시나리오 비교·바이오마커), shinydashboard purple 테마 |
| 📚 참고문헌 | [`fm_references.md`](fibromyalgia/fm_references.md) | **60개 PubMed 인용** (중추감작·하행조절·HPA·자율신경·수면·신경염증·4대 약물 임상시험) |

### 치료 시나리오 임상 데이터

| 요법 | 임상시험 | 주요 결과 | 적응 |
|------|----------|----------|------|
| 둘록세틴 60 mg QD | Arnold 2004 (Arthritis Rheum) | Pain ≥30% 개선 55% vs 33% (위약); FIQR ↓26% | FDA 승인 FM |
| 프레가발린 300–450 mg/d | Crofford 2005 (Arthritis Rheum) | Pain ≥50% 개선 29% vs 13% (위약); 수면 개선 | FDA 승인 FM |
| 밀나시프란 100–200 mg/d | Clauw 2008 (Arch Intern Med) | Pain ≥50% 24% vs 10% (위약); NRS ↓33% | FDA 승인 FM |
| 아미트립틸린 25 mg QHS | Nishishinya 2008 (Rheumatology) | 통증·수면·FIQ 개선(위약 대비); 단기 효과 강함 | 오프라벨 |
| 둘록세틴 + 프레가발린 병용 | 모델 시뮬레이션 | 상보적 기전(SERT/NET + α2δ-1)으로 추가 효과 예측 | 임상 가이드라인 미확립 |
| 아미트립틸린 수면 표적 | Moldofsky 1975 이후 연구 | SWS 회복 → 통증·피로 개선 선행 | 수면 장애 동반 시 |

---

## 모델 #129 상세: 혈우병 A (Hemophilia A)

> **디렉토리:** [`hemophilia-a/`](hemophilia-a/) | **약어:** HA | **날짜:** 2026-06-23

혈우병 A는 F8 유전자 돌연변이로 인한 응고인자 VIII 결핍이 주원인인 X염색체 연관 출혈 질환이다. FVIII는 내인성 응고경로에서 FIXa와 결합해 Xase 복합체를 형성하고 FX→FXa 전환을 증폭시키는 핵심 보조인자다. FVIII 결핍 시 트롬빈 생성이 심각히 감소하여 중증 환자에서 연간 약 30회 자연 출혈이 발생하고, 반복적인 혈관절증이 불가역적 관절 손상으로 이어진다.

### 핵심 병태생리 경로

| 클러스터 | 핵심 분자 | 치료 표적 |
|---------|---------|---------|
| **혈관 손상·혈소판** | TF 노출·vWF·GPIb·GPIIbIIIa·TXA2·ADP | 혈소판 초기 지혈 플러그 |
| **외인성 경로(TF/FVIIa)** | TF·FVIIa·TFPI 억제·FIX/FX 활성화 | 마스타시맙(항TFPI·외인성 경로 증폭) |
| **FVIII 생물학·내인성 경로** | FVIIIa·FIXa·Xase 복합체·VWF 보호 | SHL/EHL FVIII 대체·에미시주맙 모방 |
| **공통 경로·트롬빈 생성** | FXa·Prothrombinase·Thrombin burst·FXIIIa | Thrombin ETP 임상 지표 |
| **자연 항응고 기전** | AT·Protein C/S·EPCR·tPA/Plasmin | 피투시란(siRNA→AT 감소) |
| **억제항체 면역학** | TH세포·B세포·항FVIII IgG·BU 역가 | ITI·리투시맙·에미시주맙 우회 |
| **약물 PK/PD** | SHL/EHL FVIII(2구획)·에미시주맙(SC, t½ 4–5주)·피투시란(siRNA 간접반응) | HAVEN 1/3/4·ATLAS-INH 임상 보정 |
| **출혈 표현형** | 혈관절증·근육 혈종·두개내출혈·GI 출혈·ABR | 예방요법 목표 ABR <3 |
| **혈우병성 관절병증** | 활막 철 침착·ROS·연골파괴·Pettersson 점수 | 관절 보호 예방요법·물리치료 |
| **임상 엔드포인트** | FVIII trough ≥1%·ETP·ABR·QoL(EQ-5D) | Zero-bleed 표현형 달성 |

### QSP 모델 구성 (4종 산출물)

| 산출물 | 파일 | 내용 |
|--------|------|------|
| 🗺️ 기계론적 지도 | [`ha_qsp_model.dot/.svg/.png`](hemophilia-a/) | **167 노드, 10개 클러스터** (혈관손상·외인성경로·FVIII/내인성경로·공통경로·항응고기전·억제항체면역학·약물PK/PD·출혈표현형·관절병증·임상엔드포인트) |
| ⚙️ mrgsolve ODE | [`ha_mrgsolve_model.R`](hemophilia-a/ha_mrgsolve_model.R) | **16구획 ODE** (FVIII 2구획·에미시주맙 3구획·피투시란 2구획·AT mRNA/단백질·억제항체·ETP·CumBleeds·관절점수·QoL·Synovitis·FVIII_eff), **7 치료 시나리오** (무예방·SHL-FVIII·EHL-FVIII·에미시주맙 Q1W/Q4W·피투시란·병용) |
| 📊 Shiny 앱 | [`ha_shiny_app.R`](hemophilia-a/ha_shiny_app.R) | **6탭** (환자 프로파일·FVIII PK·PD 핵심지표·출혈 위험/ABR·시나리오 비교·바이오마커), bslib darkly 테마, plotly 인터랙티브, 내장 ODE 시뮬레이터 |
| 📚 참고문헌 | [`ha_references.md`](hemophilia-a/ha_references.md) | **55개 PubMed 인용** (HAVEN 1/3/4·ATLAS-INH·A-LONG·HOPE-B·Manco-Johnson 2007 NEJM·WFH Guidelines 2020) |

### 치료 시나리오 임상 데이터

| 요법 | 임상시험 | 주요 결과 | 비고 |
|------|----------|----------|------|
| SHL-FVIII 25 IU/kg 3×/주 | Manco-Johnson 2007 (NEJM) | ABR ~3.0 vs 32.4 (on-demand) | 소아 예방요법 표준 |
| EHL-FVIII Fc 50 IU/kg Q3-4d | A-LONG (Mahlangu 2014 Blood) | ABR 1.6 (개인화 예방) | t½ ~18–19 h |
| 에미시주맙 1.5 mg/kg Q1W | HAVEN 3 (Mahlangu 2018 NEJM) | ABR 1.5 vs 38.2 (미치료); 96% 감소 | 비억제항체 HA |
| 에미시주맙 6 mg/kg Q4W | HAVEN 4 (Pipe 2019 Blood Adv) | ABR 2.4; Q4W 투여 편의성 ↑ | 월 1회 피하주사 |
| 에미시주맙 3 mg/kg Q1W | HAVEN 1 (Oldenburg 2017 NEJM) | ABR 2.9 vs 23.3; 87% 감소 | 억제항체 HA |
| 피투시란 80 mg SC Q1M | ATLAS-INH (Young 2023 NEJM) | ABR 0.0 vs 17.8; 99% 감소 | siRNA; AT ↓~75% |
| valoctocogene roxaparvovec | HOPE-B (Ozelo 2022 NEJM) | FVIII 수준 23 IU/dL 달성; ABR 근접 zero | AAV5 유전자 치료 |

---

## 모델 #130 상세: 만성 림프구성 백혈병 (Chronic Lymphocytic Leukemia, CLL)

> **디렉토리:** [`chronic-lymphocytic-leukemia/`](chronic-lymphocytic-leukemia/) | **약어:** CLL | **날짜:** 2026-06-23

만성 림프구성 백혈병(CLL)은 CD19⁺CD5⁺CD23⁺ 단클론 B세포가 혈액·골수·림프절에 축적되는, 서구권 성인에서 가장 흔한 백혈병이다. BCR 자율 신호전달, BCL-2 과발현, 종양 미세환경 의존성이 핵심 병리이며, BTK 억제제(이브루티닙/아칼라브루티닙/자누브루티닙), BCL-2 억제제(베네토클락스), 항CD20 항체(오비누투주맙/리툭시맙) 3종 계열의 표적치료로 패러다임이 완전히 바뀐 질환이다.

[![CLL QSP Map](chronic-lymphocytic-leukemia/cll_qsp_model.png)](chronic-lymphocytic-leukemia/cll_qsp_model.svg)

### 핵심 병태생리 경로

| 클러스터 | 핵심 분자 | 치료 표적 |
|---------|---------|---------|
| **정상 B세포 발생** | HSC→CLP→Pro-B→Pre-B→미성숙→나이브 B세포→GC | 악성 전환 이전 단계 |
| **CLL 세포 생물학·병기** | CD19+CD5+CD23+; ALC; LDT; Rai 0-IV; Binet A-C; CD38; ZAP-70; CD49d | 병기·예후 지표 |
| **BCR 신호전달 경로** | LYN→SYK→BLNK→BTK→PLCγ2→IP₃/DAG→Ca²⁺/NFAT·PKCβ/NF-κB; PI3Kδ→PIP₃→AKT→mTOR; RAS/RAF→MEK/ERK | 이브루티닙(BTK)·이델라리십(PI3Kδ) |
| **세포사멸 조절 (BCL-2 패밀리)** | BCL-2·BCL-XL·MCL-1·A1 대 BAX/BAK; BIM·PUMA·NOXA·BAD·HRK BH3-only; MOMP→사이토크롬C→카스파제 3/7 | 베네토클락스(BCL-2) |
| **종양 미세환경** | CXCL12/CXCR4 골수 억류; CXCL13/CXCR5 림프절 귀소; CD40L·BAFF·IL-4·IL-21 생존 신호; NK 세포 ADCC; T조절세포 면역억제 | 신호전달 차단 |
| **유전·분자 예후** | IGHV 변이/비변이; del(13q14)·del(11q22)·del(17p13)·Trisomy 12; TP53/ATM/NOTCH1/SF3B1/BIRC3 돌연변이; miR-15a/16-1; CLL-IPI | 치료 선택 지도 |
| **BTK 억제제 PK/PD** | 이브루티닙 420mg QD(1구획 PK; Cys481 공유결합·t½~7h); BTK 단백질 재합성(t½~2.9일); 림프구 재분포(골수/림프절→말초혈액); off-target(ITK/TEC/EGFR)→AF·출혈 위험; C481S/PLCγ2 획득내성 | BTK 95% 이상 점유 필요 |
| **베네토클락스 PK/PD** | 400mg QD 서서히 증량(20→50→100→200→400mg); 2구획 PK(t½~26h); BCL-2 Ki~0.01nM; BH3 결합홈 차단→BIM/BAX/BAK 유리; TLS 위험(ALC·LN·신기능); MCL-1 보상 상향 내성; BCL-2 G101V 획득내성 | BCL-2 점유·TLS 예방 |
| **항CD20 mAb PK/PD** | 오비누투주맙(Type II·당화공학·ADCC↑·PCD↑·CDC↓·낮은 내재화); TMDD 2구획 PK; NK 세포 FcγRIIIA 활성화; C1q→C3d 보체 경로; CD20 항원 shedding(내성) | ADCC·PCD 세포 사멸 |
| **임상 엔드포인트** | IWCLL 2018 CR/PR/SD/PD; MRD 미검출(<10⁻⁴); PFS·OS; Richter 전환(DLBCL/HL ~5%/5년); TTNT; QoL(FACT-G) | 반응 지표 |

### QSP 모델 구성 (4종 산출물)

| 산출물 | 파일 | 내용 |
|--------|------|------|
| 🗺️ 기계론적 지도 | [`cll_qsp_model.dot/.svg/.png`](chronic-lymphocytic-leukemia/) | **146 노드, 10개 클러스터** |
| ⚙️ mrgsolve ODE | [`cll_mrgsolve_model.R`](chronic-lymphocytic-leukemia/cll_mrgsolve_model.R) | **18구획 ODE**, 6치료 시나리오, 4대 임상시험 보정 |
| 📊 Shiny 앱 | [`cll_shiny_app.R`](chronic-lymphocytic-leukemia/cll_shiny_app.R) | **6탭** (환자 프로파일·PK·PD 바이오마커·임상 엔드포인트·시나리오 비교·유전 위험) |
| 📚 참고문헌 | [`cll_references.md`](chronic-lymphocytic-leukemia/cll_references.md) | **44개 PubMed 인용** (9개 섹션) |

### 치료 시나리오 임상 데이터

| 요법 | 임상시험 | 주요 결과 | 비고 |
|------|----------|----------|------|
| 이브루티닙 420mg QD | RESONATE-2 (Burger 2015 NEJM) | ORR 86%, 2yr PFS 74% vs 클로람부실 43% | 1차 치료 기준 |
| 베네토클락스 단독 | Roberts 2016 NEJM | R/R CLL ORR 79%, MRD-neg 5% | 서서히 증량 프로토콜 |
| 오비누투주맙 6사이클 | CLL11 (Goede 2014 NEJM) | ORR 75%, PFS HR 0.16 vs chlorambucil mono | 동반질환 노인 |
| **VEN+OBI (CLL14)** | Fischer 2019 NEJM | 2yr PFS 88.2% vs 64.1%; uMRD 혈액 76%, 골수 57% | 1차, 동반질환 노인 표준 |
| 이브루티닙+베네토클락스 | CLARITY (Hillmen 2019 JCO) | R/R CLL MRD-neg 53%(골수); 2yr PFS 89% | MRD 가이드 고정기간 |
| 자누브루티닙 vs 이브루티닙 | ALPINE (Brown 2023 NEJM) | PFS NI 달성; AF 2.5% vs 10.1% | 선택성 개선 |

---

## 모델 #131 상세: 정맥 혈전색전증 (Venous Thromboembolism, VTE)

> **디렉토리:** [`venous-thromboembolism/`](venous-thromboembolism/) | **약어:** VTE | **날짜:** 2026-06-23

정맥 혈전색전증(VTE)은 심부정맥 혈전증(DVT)과 폐색전증(PE)을 통칭하는 질환으로, 연간 1-2명/1,000명에서 발생하는 제3위 심혈관 질환이다. Virchow's Triad(혈류정체·내피손상·과응고 상태)가 병리기전의 핵심이며, FV Leiden·PT G20210A·단백C/S 결핍 등 유전 위험인자와 수술·악성종양·경구피임약 등 후천적 위험인자가 복합된다. 응고 연쇄 반응(Coagulation Cascade), 혈소판 활성화, 자연 항응고 시스템(AT-III/단백C-S/TFPI), 섬유용해(플라스민-D-이량체 축)가 혈전 생성과 분해의 균형을 결정한다.

[![VTE QSP Map](venous-thromboembolism/vte_qsp_model.png)](venous-thromboembolism/vte_qsp_model.svg)

### 핵심 병태생리 경로

| 클러스터 | 핵심 분자/세포 | 치료 표적 |
|---------|--------------|---------|
| Virchow's Triad | 혈류정체·내피손상·과응고 | 위험 인자 제거/예방 |
| 내피세포 | TF 발현·vWF 방출·PGI2/NO·TM/EPCR·tPA/PAI-1 | 항혈전 내피 기능 보존 |
| 외인성 경로 | TF-FVIIa tenase → FXa/FIXa | TFPI |
| 내인성 경로 | FXIIa→FXIa→FIXa-FVIIIa Intrinsic Tenase | 접촉인자 억제 |
| 공통 경로 | FXa-FVa Prothrombinase → FIIa(트롬빈) → 피브린 | FXa 억제제·DTI |
| 혈소판 | GPIb/GPVI/PAR1·ADP/P2Y12·GPIIb/IIIa·PS 노출 | P2Y12 억제제(PE 방지) |
| 자연 항응고 | AT-III·APC-단백S·TFPI | 결핍 시 보충/치료 |
| 섬유용해 | tPA/uPA → 플라스민 → FDP/D-이량체 | 혈전용해제(알테플라제) |
| 약물 PK | 리바록사반·아픽사반·다비가트란·와파린·에녹사파린 | 5종 항응고제 PK/PD |
| 임상 결과 | 근위부 DVT → PE → 만성 혈전색전성 폐고혈압(CTEPH) | 기간별 항응고 전략 |

### QSP 모델 구성

| 구성요소 | 파일 | 사양 |
|---------|------|-----|
| 🗺️ 기계론적 지도 | [`vte_qsp_model.dot`](venous-thromboembolism/vte_qsp_model.dot) | **140+ 노드, 12클러스터** |
| ⚙️ mrgsolve ODE | [`vte_mrgsolve_model.R`](venous-thromboembolism/vte_mrgsolve_model.R) | **19구획 ODE**, **6치료 시나리오** |
| 📊 Shiny 앱 | [`vte_shiny_app.R`](venous-thromboembolism/vte_shiny_app.R) | **6탭** (환자 프로파일·PK·응고 PD·혈전 동태·시나리오 비교·바이오마커) |
| 📚 참고문헌 | [`vte_references.md`](venous-thromboembolism/vte_references.md) | **57개 PubMed 인용** (15개 섹션) |

### 치료 시나리오 임상 데이터

| 요법 | 임상시험 | 주요 결과 | 비고 |
|------|----------|----------|------|
| 리바록사반 15mg BID×21d → 20mg QD | EINSTEIN-DVT/PE (NEJM 2010/2012) | DVT재발/PE 방지 비열등성 NNT ~30 | 브리지 불필요 |
| 아픽사반 10mg BID×7d → 5mg BID | AMPLIFY (Agnelli 2013 NEJM) | 재발 위험 RR 0.84; 주요 출혈 ↓ 69% | 신장 친화적 |
| 다비가트란 150mg BID (LMWH 이후) | RE-COVER (Schulman 2009 NEJM) | 재발 VTE 비열등성 HR 1.10; 출혈 유사 | 신기능 의존 |
| 와파린 INR 2-3 + 에녹사파린 브리지 | 역사적 표준 | DVT 재발 3개월 2-3% | INR 모니터링 필수 |
| 에녹사파린 40mg QD (수술예방) | MEDENOX (Samama 1999 NEJM) | 증상 VTE↓ 63% vs 위약 | 입원 내과 환자 |
| 리바록사반 10mg QD (연장예방) | EINSTEIN-EXT (Bauersachs 2010 NEJM) | 재발 82% 감소 vs 위약 | 출혈 위험 주의 |

---

## 🧬 척수성 근위축증 (Spinal Muscular Atrophy) — 최신 모델 상세 (2026-06-23)

[![SMA QSP 기계론적 지도](spinal-muscular-atrophy/sma_qsp_model.png)](spinal-muscular-atrophy/sma_qsp_model.svg)

**질환**: 척수성 근위축증(SMA) | **유전자**: *SMN1* (5q13.2 동형접합 결손) | **OMIM**: [253300](https://omim.org/entry/253300) | **유병률**: 1/6,000–1/10,000 출생

### 주요 경로 클러스터 (12개)

| 클러스터 | 핵심 기전 |
|---------|-----------|
| 1. 유전적 기반 | SMN1 결손, SMN2 복사수(1–4), C840T 전환변이, 5q13.2 위치 |
| 2. SMN2 스플라이싱 | 엑손7 ESE/ESS/ISS-N1, hnRNP A1/A2, SRSF1/Tra2β, 10% 기저 FL-SMN |
| 3. SMN 단백질·복합체 | FL-SMN 올리고머, Gemin2–8, Sm 단백질, snRNP 조립, 스플라이소좀 |
| 4. 운동신경세포 생물학 | α-MN pool, BCL-2/BAX/Casp-3 아포토시스, BDNF/VEGF 생존 신호 |
| 5. 신경근육접합부 | AChR 클러스터링(MuSK·Agrin·Lrp4), SNAP-25, 종판 전위 |
| 6. 골격근 | Type I/II 섬유, 신경원성 위축, UPS/자가포식, IGF-1/mTORC1 |
| 7. 임상 지표 | CMAP, CHOP-INTEND, HFMSE, RULM, MFM-32, FVC, 운동 마일스톤 |
| 8. 누시너센 PK/PD | IT 주사 → CSF 허리/경추 → CNS 조직, ISS-N1 차단, 엑손7 포함율 ↑ |
| 9. 리스디플람 PK/PD | 경구 흡수 → 혈장/CNS, SRSF1/Tra2β 향상, 전신적 FL-SMN ↑ |
| 10. 오나셈노진 PK/PD | IV → AAV9 BBB 통과 → MN 형질도입 → SMN1 트랜스진 지속 발현 |
| 11. 질환 유형·자연경과 | Type 0–IV (SMN2 복사수 역비례), 신생아 선별검사, 호흡/척추측만 |
| 12. 전신 효과 | 심장 기형, 대사 이상, 골다공증, 삼킴 장애, 삶의 질 |

### mrgsolve ODE 모델 (20구획)

| 모듈 | 구획 | 핵심 방정식 |
|------|------|------------|
| 누시너센 PK | A_CSF_L, A_CSF_C, A_CNS_NUS | CSF bulk flow(619 mL/day) + 비선형 CNS 흡수(Michaelis-Menten) |
| 리스디플람 PK | A_gut_RIS, A_plasma_RIS, A_CNS_RIS | 1차 흡수(ka=1.5/h), Kp_brain=0.5 |
| 오나셈노진 PK | A_plasma_ZOL, A_MN_ZOL, A_tg_mRNA | 빠른 소거(t½~7h) + MN 형질도입 + 트랜스진 mRNA 동역학 |
| SMN 생물학 | FL_SMN_mRNA, dSMN_mRNA, SMN_prot | dSMN/dt = k_syn×E7I − k_deg×SMN |
| 운동신경·NMJ·근육 | MN_pool, NMJ_score, Muscle_mass | MN 사멸(SMN 역치 이하 Emax), NMJ 성숙도, 신경원성 위축 |
| 누적 지표 | AUC_SMN, MN_lost | 시간 적분 |

### 치료 시나리오 임상 데이터

| 시나리오 | 약물 | 임상시험 | 주요 결과 |
|---------|------|----------|----------|
| SMA Type I 자연경과 | 없음 | 역사적 대조 | 2세 이전 사망/영구환기 ~100% |
| SMA Type I — 누시너센 | 누시너센 IT 12mg | ENDEAR (2017 NEJM) | 운동 마일스톤 반응 51% vs 0% |
| SMA Type II — 리스디플람 | 리스디플람 5mg/day | SUNFISH (2022 Lancet Neurol) | MFM32 +1.36 vs −0.19 |
| 무증상 SMA — 오나셈노진 | Zolgensma 1.1×10¹⁴ vg/kg | SPR1NT (2022 Nat Med) | 독립 보행 달성률 77%(3복사수) |
| Type II 늦은 시작(1년 후) | 누시너센 | CHERISH 후속 | HFMSE 개선 지연·약화 |
| 소아 체중 기반 | 리스디플람 0.2mg/kg | FIREFISH (2021 NEJM) | 좌위 유지 ≥5초: 61%(12개월) |

### QSP 모델 구성

| 구성요소 | 파일 | 사양 |
|---------|------|-----|
| 🗺️ 기계론적 지도 | [`sma_qsp_model.dot`](spinal-muscular-atrophy/sma_qsp_model.dot) | **130+ 노드, 12클러스터** |
| ⚙️ mrgsolve ODE | [`sma_mrgsolve_model.R`](spinal-muscular-atrophy/sma_mrgsolve_model.R) | **20구획 ODE**, **6치료 시나리오** |
| 📊 Shiny 앱 | [`sma_shiny_app.R`](spinal-muscular-atrophy/sma_shiny_app.R) | **8탭** (환자 프로파일·PK·SMN 생물학·운동신경/NMJ·임상 지표·시나리오 비교·바이오마커·가상 집단) |
| 📚 참고문헌 | [`sma_references.md`](spinal-muscular-atrophy/sma_references.md) | **50개 PubMed 인용** (10개 섹션) |

---

## 🫀 가와사키병 (Kawasaki Disease) — 최신 모델 상세 (2026-06-24)

> **디렉토리:** [`kawasaki-disease/`](kawasaki-disease/) | **약어:** KD | **날짜:** 2026-06-24

[![KD QSP 기계론적 지도](kawasaki-disease/kd_qsp_model.png)](kawasaki-disease/kd_qsp_model.svg)

**질환**: 가와사키병(Kawasaki Disease, KD) | **분류**: 소아 전신 혈관염 | **호발 연령**: 6개월–5세 | **CAA 위험**: 미치료 시 ~25%

### 핵심 기전 (14개 클러스터)

| 클러스터 | 핵심 기전 |
|---------|-----------|
| 1. 감염·환경 트리거 | RNA/DNA 바이러스, 세균 초항원, 환경 인자, ITPKC/CASP3/HLA-B15 유전 소인 |
| 2. 선천 면역 활성화 | TLR2/4·NLR → 대식세포·수지상세포·중성구·NK세포·보체(C3/C5/MAC) |
| 3. NLRP3 인플라마좀 | NLRP3→ASC→Caspase-1→Pro-IL-1β 절단→성숙 IL-1β; Gasdermin D 공극형성; 파이롭토시스 |
| 4. 사이토카인 네트워크 | IL-1β·IL-1R→NF-κB; IL-6→JAK1/2-STAT3; TNF-α→TNFR1→NF-κB; IL-8·IL-18·IFN-γ·MCP-1 |
| 5. 적응 면역 | Th1/Th17/Treg 분화; 형질아세포→IgG·IgA; FcγR; MHC II 항원 제시 |
| 6. 내피세포 활성화 | 정지 EC → 활성 EC; VCAM-1·ICAM-1·E-셀렉틴·조직인자 발현; eNOS↓; 혈관 투과성↑ |
| 7. 관상동맥 병리 | 관상동맥 염증→중막 파괴→탄성막 파열→동맥류(small/medium/giant); 혈전→심근경색 |
| 8. 혈소판 생물학 | 혈소판 활성화·TXA2·GPIIb/IIIa·vWF; 혈소판 증가증(2주 피크) |
| 9. 발열·급성기 반응 | COX-2→PGE2→시상하부 발열 중추; CRP·ESR·페리틴·프로칼시토닌·알부민↓ |
| 10. IVIG PK | 2구획·FcRn 재순환; t½ 21–28일; EC50=8 g/L; Emax=80% 사이토카인 억제 |
| 11. 아스피린 PK | 고용량(80–100 mg/kg/day) → 저용량(3–5 mg/kg/day); 살리실산 대사체; COX-1/2 |
| 12. 스테로이드 PK | 메틸프레드니솔론 2구획; GR-α → GRE 전사활성화/NF-κB 억제 |
| 13. 생물학적 제제 PK | 인플릭시맙 2구획(TNF-α 중화); 아나킨라 SC(IL-1R 경쟁적 차단); 사이클로스포린 |
| 14. 임상 엔드포인트 | 발열 기간; Kobayashi/Egami 위험점수; CAA Z-점수(AHA 분류); IVIG 저항성 ~15%; 재발률 3% |

### mrgsolve ODE 모델 (21 구획)

| 모듈 | 구획 | 핵심 동역학 |
|------|------|------------|
| IVIG PK | A_IVIG_c, A_IVIG_p | 2구획·FcRn 재순환(F=60%)·CL=0.0033 L/h/kg |
| 아스피린 PK | A_ASA_gut, A_ASA_c, A_SA_c | ka=0.80/h·COX-1 비가역·살리실산 대사(CL_SA=0.01) |
| 메틸프레드니솔론 | A_MP_c, A_MP_p | 2구획·CL=0.48 L/h/kg·t½~2h |
| 인플릭시맙 | A_IFX_c, A_IFX_p | 2구획·EC50=2.5 μg/mL·Hill n=1.8 |
| 아나킨라 | A_ANK_gut, A_ANK_c | SC 흡수(ka=0.30/h)·EC50=1.0 μg/mL |
| 사이토카인 | IL1b, IL6, TNFa | 대식세포 구동 생산·약물 Emax 억제 ODE |
| 활성화 상태 | Mac_act, EC_act | 로지스틱 성장(Mac)·사이토카인 구동(EC) |
| 임상 지표 | Fever, CRP, PLT_c, CAL_Z | PGE2/IL-6/IL-1β 구동; IL-6 혈소판 생성; 관상동맥 Z-점수 |

### 5가지 치료 시나리오 임상 근거

| 시나리오 | 치료법 | 임상시험/근거 | 주요 결과 |
|---------|--------|-------------|----------|
| S1: 표준 | IVIG 2 g/kg + 아스피린 | Newburger 1991 (NEJM) | CAA 발생 3–5%; 발열 소실률 85% |
| S2: 고위험 + 스테로이드 | + 메틸프레드니솔론 | RAISE Trial (Kobayashi 2012 Lancet) | CAA 위험 0% vs 3% (스테로이드) |
| S3: IVIG 저항 → 2차 IVIG | IVIG 2 g/kg × 2회 | Burns 1998 (PIDJ) | 2차 반응률 ~50% |
| S4: IVIG 저항 → 인플릭시맙 | 인플릭시맙 5 mg/kg | KIDCARE (Tremoulet 2019 Lancet) | 2차IVIG vs IFX: IFX 비열등; CAA 유사 |
| S5: IVIG 저항 → 아나킨라 | 아나킨라 4 mg/kg/day | Ouldali 2019 (J Pediatr) | 발열 소실 96%; 관상동맥 안정화 |

### QSP 모델 파일

| 구성요소 | 파일 | 사양 |
|---------|------|-----|
| 🗺️ 기계론적 지도 | [`kd_qsp_model.dot`](kawasaki-disease/kd_qsp_model.dot) | **134 노드, 14클러스터** |
| ⚙️ mrgsolve ODE | [`kd_mrgsolve_model.R`](kawasaki-disease/kd_mrgsolve_model.R) | **21구획 ODE**, **5치료 시나리오** |
| 📊 Shiny 앱 | [`kd_shiny_app.R`](kawasaki-disease/kd_shiny_app.R) | **6탭** (환자 프로파일·PK·사이토카인/염증·임상 엔드포인트·시나리오 비교·바이오마커/위험도) |
| 📚 참고문헌 | [`kd_references.md`](kawasaki-disease/kd_references.md) | **60개 PubMed 인용** (14개 섹션) |

---

## 🫀 쿠싱 증후군 (Cushing's Syndrome) — 최신 모델 상세 (2026-06-24)

> **디렉토리:** [`cushings-syndrome/`](cushings-syndrome/) | **약어:** CS | **날짜:** 2026-06-24

[![CS QSP 기계론적 지도](cushings-syndrome/cs_qsp_model.png)](cushings-syndrome/cs_qsp_model.svg)

**질환**: 쿠싱 증후군(Cushing's Syndrome, CS) | **분류**: 내분비·부신 질환 | **유병률**: 100만 명당 10–15명/년 | **주요 합병증**: 고혈압, 당뇨, 골다공증, 우울증

### 핵심 기전 (13개 클러스터)

| 클러스터 | 핵심 기전 |
|---------|-----------|
| 1. 시상하부 | CLOCK/BMAL1 일주기 리듬, CRH/AVP 분비, 소마토스타틴·도파민 억제, 스트레스 입력(NE) |
| 2. 뇌하수체 전엽 | CRHR1 → Gs → PKA → CREB → POMC → PC1/2 → ACTH; USP8 돌연변이(~50% 쿠싱병); CDK4/6-Rb 증식; 뇌하수체 선종 |
| 3. 부신피질 생합성 | MC2R → cAMP → PKA → StAR; 콜레스테롤 → CYP11A1 → 프레그네놀론 → CYP17A1 → 17-OHP → CYP21A2 → CYP11B1 → 코르티솔 |
| 4. GR 신호 | GR-α/β, HSP90, FKBP51(음성)/FKBP52(양성), 핵이동 이량체, GRE/nGRE, AP-1·NF-κB 접촉억제, GILZ, MKP1, SGK1 |
| 5. 대사 합병증 | PEPCK/G6Pase↑(간당신생), GLUT4↓, 인슐린저항성, 내장지방↑, Atrogin-1/MuRF1 근육위축, RANKL↑/OPG↓ 골다공증 |
| 6. 심혈관·신장 | RAAS 과활성(Ang II→AT1R→알도스테론→Na저류), 내피기능장애(NO↓/ET-1↑), 고혈압, VTE, 이상지질혈증 |
| 7. 면역 억제 | NF-κB/AP-1 접촉억제→IL-6/TNF-α↓, 림프구감소증, NK세포↓, 호중구 탈변연화, 감염 위험↑ |
| 8. CNS 효과 | 해마위축·BDNF↓·NMDA 흥분독성, 세로토닌계↓, 우울증(50-70%)·인지장애·불면 |
| 9. 파시레오티드 PK/PD | 2구획 SC; SSTR5>SSTR2; Gi→cAMP↓→ACTH억제(최대 65%); 혈당상승 부작용(SSTR5) |
| 10. 스테로이드 합성 억제제 | 케토코나졸(CYP17A1+11B1, Emax 72%), 메티라폰(CYP11B1, Emax 82%), 오실로드로스탯(CYP11B1/B2, Emax 85%, EC50=0.15 μg/mL) |
| 11. GR 길항제·카버골린 | 미페프리스톤(GR 경쟁적 길항, Emax 82%), 카버골린(D2R→ACTH 억제), 미토탄(부신피질 세포독성) |
| 12. 임상 진단·평가 | UFC 24h, LNSC, 1mg/8mg DST, CRH 자극 검사, IPSS, ACTH, BMD-DXA, HDRS-17 |
| 13. 병인 분류 | 쿠싱병(USP8 변이, 뇌하수체 선종), 이소성 ACTH(소세포폐암/카르시노이드), 부신 선종/암, PBMAH, 주기성 쿠싱 |

### mrgsolve ODE 모델 (21구획)

| 모듈 | 구획 | 핵심 방정식 |
|------|------|------------|
| HPA 축 | CRH, ACTH_PIT, ACTH_PL | 일주기 CRH 합성(cos 파형) + 종양 ACTH 추가 + GR 음성 피드백(Hill n=2) |
| 부신 스테로이드 | F_ADR, F_PL | Michaelis-Menten(ACTH→코르티솔); 약물 Bliss 결합 억제율 |
| GR 동역학 | GR_FREE, GR_BOUND, GR_NUC | 2단계(세포질→핵) ODE; 미페프리스톤 길항 감쇠 |
| 대사 | GLUCOSE, INSULIN, VAT, MUSCLE, BMD, BP | GR-eff 구동 대사 변화; 인슐린 피드백 루프 |
| 임상 출력 | UFC_ACC | UFC_coef × 유리코르티솔 누적 |
| 약물 PK | A_PAS_C/P, A_KETO, A_METY, A_OSILO, A_MIFE | 1–2구획 모델 (5종 약물) |

### 치료 시나리오 임상 데이터

| 시나리오 | 약물 | 임상시험 | 주요 결과 |
|---------|------|----------|----------|
| 자연경과 (쿠싱병 무치료) | — | 역사적 코호트 | UFC >500 μg/24h, 합병증 누적, 심혈관 사망 위험 4배 |
| 파시레오티드 0.6mg BID | 파시레오티드 SC | PASPORT-CUSHINGS (Colao 2012 NEJM) | UFC 정상화 22–24% (6개월), 혈당상승 부작용 |
| 케토코나졸 400mg BID | 케토코나졸 | Castinetti 2014 Eur J Endocrinol | UFC 정상화 49%, 신속 효과, 간독성 주의 |
| 오실로드로스탯 5mg BID | 오실로드로스탯 | LINC 3 (Pivonello 2020 Lancet DE) | UFC 정상화 86%(유지기), 가장 높은 효능 |
| 미페프리스톤 600mg QD | 미페프리스톤 | SEISMIC (Fleseriu 2012 JCEM) | 혈당/BP 임상반응 87%; 코르티솔 오히려 상승 (GR 길항) |
| 수술 후 관해 | — | 메타분석 | 재발률 ~20% (5년 내); HPA 축 회복 수개월 소요 |

### QSP 모델 구성

| 구성요소 | 파일 | 사양 |
|---------|------|-----|
| 🗺️ 기계론적 지도 | [`cs_qsp_model.dot`](cushings-syndrome/cs_qsp_model.dot) | **140+ 노드, 13클러스터** |
| ⚙️ mrgsolve ODE | [`cs_mrgsolve_model.R`](cushings-syndrome/cs_mrgsolve_model.R) | **21구획 ODE**, **6치료 시나리오** |
| 📊 Shiny 앱 | [`cs_shiny_app.R`](cushings-syndrome/cs_shiny_app.R) | **8탭** (환자 프로파일·HPA/PK·스테로이드생합성·임상지표·시나리오비교·바이오마커·대사합병증·가상집단) |
| 📚 참고문헌 | [`cs_references.md`](cushings-syndrome/cs_references.md) | **55개 PubMed 인용** (10개 섹션) |

---

## 🩸 유전성 혈관부종 (Hereditary Angioedema) — 최신 모델 상세 (2026-06-24)

> **디렉토리:** [`hereditary-angioedema/`](hereditary-angioedema/) | **약어:** HAE | **날짜:** 2026-06-24

[![HAE QSP 기계론적 지도](hereditary-angioedema/hae_qsp_model.png)](hereditary-angioedema/hae_qsp_model.svg)

**질환**: 유전성 혈관부종(Hereditary Angioedema, HAE) | **분류**: 희귀·유전질환 | **유병률**: 1:50,000 | **주요 매개체**: 브라디키닌(B2R 의존성)

### 핵심 기전 (12개 클러스터)

| 클러스터 | 핵심 기전 |
|---------|-----------|
| 1. 유전적 기반 | *SERPING1* 유전자 돌연변이(>700종) → C1-INH 단백질 결핍(Type I)/기능이상(Type II); *F12* Thr328Lys 이득기능돌연변이(Type III·에스트로겐 감수성) |
| 2. 접촉 활성화 | FXII → FXIIa (음전하 표면 자동활성화); FXIIa ↔ Prekallikrein·HMWK 증폭 루프; 아연·열쇼크 단백질 보조인자 |
| 3. 칼리크레인-키닌계(KKS) | 혈장 칼리크레인, HMWK 절단, 브라디키닌(BK) 생성, ACE·NEP·APP 키니나아제 분해 |
| 4. C1-INH 생물학 | SERPIN 기전, C1-INH:FXIIa·C1-INH:Kal 복합체, 보체 C1r/C1s/MASP1/2 억제, C4 소진(진단 마커) |
| 5. BK 수용체 신호 | B2R(구성적·고친화성)·B1R(유도성·염증 시 상향조절); Gq/IP3/Ca²⁺/eNOS/NO/PGI2 경로 |
| 6. 혈관 효과 | EC 타이트정션 파괴, VE-카데린 내재화, 혈장 삼출 → 피하/복부/후두 부종 |
| 7. 발작 병인 | 트리거(감정적 스트레스·수술·ACEi·에스트로겐), 발작 타임라인(전구증상→피크→소실 48–72h), 후두 사망 위험 |
| 8. 급성 치료 | 이카티반트(B2R Ki=0.47nM, 30mg SC), C1-INH IV(Berinert/Ruconest), 에칼란티드(칼리크레인) |
| 9. 예방 치료 | 베로트랄스탓(IC50=3.7nM 칼리크레인 경구 QD), 라나델루맙(KD<100pM prekallikrein SC Q2W), C1-INH SC(Haegarda), 다나졸 |
| 10. 염증 증폭 | IL-1β→B1R 상향조절, 응고 크로스토크(FXIa/트롬빈), vWF/P-셀렉틴, HMWK 절단 지표 |
| 11. 임상 엔드포인트 | 발작 빈도, ACE 검사(Angioedema Control Test), AE-QoL, 후두 위험도, 바이오마커(C4/C1-INH%) |
| 12. 진단 알고리즘 | C4 → C1-INH 항원 → C1-INH 기능 → C1q → FXII 유전자 검사 |

### mrgsolve ODE 모델 (20 구획)

| 모듈 | 구획 | 핵심 동역학 |
|------|------|------------|
| 이카티반트 PK | A_ICA_depot, A_ICA_C, A_ICA_P | SC 데포 → 2구획; ka=0.74/h, CL=15.5 L/h, Vc=29 L, t½=1.3h |
| C1-INH IV PK | A_C1INH_IV | 1구획; CL=0.051 L/h, Vd=3.3 L, t½=45h |
| C1-INH SC PK | A_C1INH_SC, A_C1INH_SC_C | SC 데포; F=43%, ka=0.025/h |
| 베로트랄스탓 PK | A_BER_gut, A_BER_C | 경구 흡수 F=57%; Vd=268 L, t½=93h |
| 라나델루맙 PK | A_LAN_depot, A_LAN_C, A_LAN_P | SC → 2구획; F=61%, t½=17d |
| C1-INH 생물학 | C1INH_free | 간 합성 – 키니나아제 소비 – 분해 |
| 접촉 활성화 | FXII_act | FXIIa 동역학; C1-INH 억제; 트리거 강제 |
| KKS | Kallikrein_act | FXIIa 구동; C1-INH·베로트랄스탓·라나델루맙 억제 |
| 브라디키닌 | BK_plasma | 칼리크레인 구동 합성; ACE 분해 |
| B2R 동역학 | B2R_free, B2R_bound | BK 결합(kon/koff); 이카티반트 경쟁적 차단; 내재화 |
| 혈관 투과성 | VP | B2R_bound → NO/PGI2 구동; Emax 모델(Hill=1.8) |
| 부종 점수 | SW_score | VP 역치 구동 부종 형성; 소실 동역학 |

### 약물 PK/PD 파라미터

| 약물 | 투여 경로 | 용량 | t½ | 주요 PD |
|------|---------|------|-----|--------|
| 이카티반트 | SC | 30 mg | 1.3 h | Ki(B2R) = 0.47 nM |
| C1-INH IV (Berinert) | IV | 20 IU/kg | 45 h | C1-INH 보충 → 억제 복원 |
| C1-INH SC (Haegarda) | SC | 60 IU/kg 2×/wk | 45 h | 정상상태 C1-INH >40% 유지 |
| 베로트랄스탓 | PO | 150 mg QD | 93 h | IC50(Kal) = 3.7 nM, Emax = 92% |
| 라나델루맙 | SC | 300 mg Q2W | 17 d | KD(prekallikrein) < 100 pM, Emax = 93% |

### 6가지 치료 시나리오 임상 근거

| 시나리오 | 약물 | 임상시험 | 주요 결과 |
|---------|------|----------|----------|
| S1: 무치료 HAE | — | 역사 코호트 | BK ↑10배, 최대 부종, 후두 사망 위험 |
| S2: 이카티반트 SC | 이카티반트 30 mg | FAST-1/3 (NEJM 2010·2011) | B2R ~98% 차단; 2–4h 내 소실 |
| S3: C1-INH IV | Berinert 20 IU/kg | Cicardi 2012 NEJM | C1-INH 복원; BK ↓↓ |
| S4: 베로트랄스탓 QD | 150 mg 경구 | BELO 2020 | 발작 44% 감소; Kal IC50 정상상태 |
| S5: 라나델루맙 Q2W | 300 mg SC | HELP OLE 2020 | 발작 87% 감소; 32% 무발작 |
| S6: C1-INH SC | Haegarda 60 IU/kg | CONFIDENT 2017 | 발작 95% 감소; C4 정상화 |

### QSP 모델 구성

| 구성요소 | 파일 | 사양 |
|---------|------|-----|
| 🗺️ 기계론적 지도 | [`hae_qsp_model.dot`](hereditary-angioedema/hae_qsp_model.dot) | **120+ 노드, 12클러스터**, fdp 레이아웃 |
| ⚙️ mrgsolve ODE | [`hae_mrgsolve_model.R`](hereditary-angioedema/hae_mrgsolve_model.R) | **20구획 ODE**, **6치료 시나리오** |
| 📊 Shiny 앱 | [`hae_shiny_app.R`](hereditary-angioedema/hae_shiny_app.R) | **6탭** (환자 프로파일·PK·KKS 생물학·임상 엔드포인트·시나리오 비교·바이오마커) |
| 📚 참고문헌 | [`hae_references.md`](hereditary-angioedema/hae_references.md) | **58개 PubMed 인용** (12개 섹션) |

---

## #136 — 패혈증 / Sepsis (SEP)

> **디렉토리:** [`sepsis/`](sepsis/) | **약어:** SEP | **날짜:** 2026-06-24

[![Sepsis QSP 기계론적 지도](sepsis/sep_qsp_model.png)](sepsis/sep_qsp_model.svg)

**패혈증(Sepsis)**은 감염에 대한 숙주 반응의 조절 장애로 인해 생명을 위협하는 장기 기능 장애가 발생하는 증후군입니다(Sepsis-3, Singer et al. JAMA 2016). 전 세계적으로 연간 약 4,900만 명이 이환되고 1,100만 명이 사망하는 중환자 의학의 핵심 과제입니다.

### 병태생리 클러스터 (10 clusters)

| 클러스터 | 핵심 기전 |
|---------|-----------|
| 1. 감염 & PAMP | 세균(Gram⁻ LPS·Gram⁺ LTA)·진균(β-glucan)·바이러스(ssRNA/dsRNA) → PAMP/DAMP 생성 |
| 2. 선천면역 인식 | TLR4/TLR2/TLR9 → MD2 → MyD88 → IRAK1/4 → TRAF6 → IKK → IκB 인산화 → NF-κB 핵 이동 |
| 3. 사이토카인 네트워크 | NF-κB → TNF-α·IL-1β·IL-6·IL-8(G-CSF); 음성 피드백: IL-10·TGF-β1; HMGB1 후기 매개체 |
| 4. 응고 & DIC | TF↑ → FVIIa → FXa → Prothrombin → Thrombin → Fibrin → D-dimer; PAI-1↑·tPA↓·Plasmin↓ → 혈소판 감소·DIC |
| 5. 보체 활성화 | C1q/MBL·MASP → C3 분열 → C3a/C3b → C5 분열 → C5a(호중구 동원)·MAC(용균) |
| 6. 혈관 기능 장애 | iNOS↑ → NO↑ → sGC/cGMP → 혈관 이완 → SVR↓ → MAP↓; 글리코칼릭스 탈락 → 모세혈관 누출 → 부종 |
| 7. HPA 축 | CRH→ACTH→코르티솔 합성 → GR 활성화 → 항염증 유전자; CIRCI(30~40%에서 상대적 부신 기능 저하) |
| 8. 다장기부전 | 심기능 저하·AKI·ARDS·간기능 장애·뇌병증·장 장벽 손상·DIC·미토콘드리아 기능 장애·MODS |
| 9. 약물 PK | 피페라실린/타조박탐(2구획)·반코마이신·노르에피네프린·하이드로코르티손·바소프레신·인슐린 |
| 10. 임상 엔드포인트 | SOFA 점수·qSOFA·젖산·프로칼시토닌·CRP·WBC·혈소판·INR·크레아티닌·P/F 비·패혈성 쇼크 |

### mrgsolve ODE 모델 (20 구획)

| 구획 | 변수 | 핵심 동역학 |
|------|------|------------|
| 감염 | B | 세균 부하(CFU/mL); 로지스틱 성장 kb=0.9/h; 호중구·항생제 제거 |
| 선천면역 | N, M | 호중구(cells/µL)·활성화 대식세포; TNF 동원·IL-10 억제 |
| 사이토카인 | TNF, IL6, IL10, IL1b | pg/mL; t½=2h~6h; 코르티솔 억제(NF-κB 경로) |
| 응고 | Th, F, Plt | 트롬빈(nM)·피브린(µg/mL)·혈소판(×10³/µL); 트롬빈-매개 소비 |
| 혈관 | NO | NO(µM); iNOS(대식세포·사이토카인 구동); MAP에 직접 영향 |
| 장기 | D_tissue, Lac, MAP, Cr | 장기손상 지표(0–1)·젖산·MAP(mmHg)·크레아티닌; GFR 손상 연계 |
| 항생제 PK | AB_C, AB_P | 피페라실린 2구획; CL=15 L/h, Vc=10 L, t½=1h |
| 승압제 | NE_C | 노르에피네프린; CL=150 L/h, t½≈2min |
| 스테로이드 | HC_C, Cort | 하이드로코르티손 + 내인성 코르티솔; 총 코르티솔 풀 |

### 약물 PK/PD

| 약물 | 경로 | 용량 | t½ | 주요 PD |
|------|------|------|----|--------|
| 피페라실린/타조박탐 | IV | 4.5 g q6h | 1.0 h | %fT>MIC; Emax=0.95, MIC=16 µg/mL |
| 노르에피네프린 | IV 지속 | 0.1–0.5 µg/kg/min | ~2 min | α1-작용 → SVR↑ → MAP↑ (Emax=30 mmHg) |
| 하이드로코르티손 | IV 지속 | 200 mg/day | 1.5 h | GR → 사이토카인 억제 (Imax=65%, IC50=5 µg/dL) |
| 바소프레신 | IV 지속 | 0.03 units/min | ~15 min | V1a → 혈관 수축 → MAP↑ (Emax=15 mmHg) |

### 6가지 치료 시나리오 임상 근거

| 시나리오 | 치료 | 임상시험 | 주요 결과 |
|---------|------|----------|----------|
| S1: 무치료 | — | 역사 코호트 | 진행성 MAP↓, 다장기부전, 높은 사망률 |
| S2: 조기 항생제(1h) | PipTazo 4.5g q6h | Kumar 2006 Crit Care Med | 1h 지연마다 사망률 7% 증가 |
| S3: 항생제 + NE | PipTazo + NE 0.2 µg/kg/min | De Backer NEJM 2010 | MAP ≥65 복원; 젖산 청소 |
| S4: 완전 번들 | AB + NE + HC 200mg/day | ADRENAL 2018 NEJM | 쇼크 역전 단축; 스테로이드 의존도 ↓ |
| S5: 지연 항생제(6h) | PipTazo 6h 후 시작 | Kumar 2006 (지연 코호트) | 6h 지연 → 사망률 +14% 추정 |
| S6: 불응성 패혈성 쇼크 | NE 0.5 + VP 0.03 + HC | VASST 2008 NEJM | 바소프레신이 NE 절약; 스테로이드로 승압제 이탈 가속 |

### QSP 모델 구성

| 구성요소 | 파일 | 사양 |
|---------|------|-----|
| 🗺️ 기계론적 지도 | [`sep_qsp_model.dot`](sepsis/sep_qsp_model.dot) | **124+ 노드, 11클러스터**, 211 엣지 |
| ⚙️ mrgsolve ODE | [`sep_mrgsolve_model.R`](sepsis/sep_mrgsolve_model.R) | **20구획 ODE**, **6치료 시나리오** |
| 📊 Shiny 앱 | [`sep_shiny_app.R`](sepsis/sep_shiny_app.R) | **6탭** (환자 프로파일·PK·사이토카인·혈역학·시나리오 비교·바이오마커) |
| 📚 참고문헌 | [`sep_references.md`](sepsis/sep_references.md) | **55개 PubMed 인용** (10개 섹션) |
