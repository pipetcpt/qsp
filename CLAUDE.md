# CLAUDE.md — QSP Disease Model Library

## 프로젝트 목적 (Project Purpose)

이 저장소는 **Claude Code Routine (CCR)**을 통해 매일 자동으로 질환별 정량적 시스템 약리학(QSP) 모델이 추가되는 살아있는 라이브러리입니다.
각 세션에서 Claude가 질환 하나를 선택하고, 기계론적 지도 / mrgsolve 모델 / Shiny 앱 / 참고문헌을 완전히 작성한 뒤, main 브랜치에 직접 커밋 및 푸시합니다.

This repository is a living library where a new disease-specific Quantitative Systems Pharmacology (QSP) model is added daily via Claude Code Routine (CCR). Each session selects one disease, builds all four deliverables (mechanistic map · mrgsolve ODE model · Shiny app · references), then commits and pushes directly to main.

---

## 세션 종료 시 반드시 수행할 사항 (Session Stop Requirements)

The stop hook (`~/.claude/stop-hook-git-check.sh`) enforces that every session ends with:
1. All changes committed
2. All commits pushed to the remote branch

이 조건이 충족되지 않으면 세션이 종료되지 않습니다.

---

## 새 모델 추가 지침 (Guidelines for Adding a New Model)

### 질환 선택 (Disease Selection)
- 아래 질환 목록(`분류/질환명/발병기전`)에서 **아직 만들지 않은** 것을 선택
- 오늘 날짜와 이전에 만든 것을 참고해 **매일 다른 카테고리**에서 선택
- 이미 추가된 질환은 README.md 표를 참조

### 디렉토리 규칙 (Directory Convention)
- 소문자 + 하이픈: `질환명-영어/` (예: `iga-nephropathy/`, `crohn-disease/`)
- 중복 디렉토리 불가 — `ls` 로 확인 후 생성

### 파일 명명 (File Naming)
| 파일 | 규칙 |
|------|------|
| DOT  | `<abbr>_qsp_model.dot` or `<abbr>_qsp.dot` |
| SVG  | `<abbr>_qsp_model.svg` |
| PNG  | `<abbr>_qsp_model.png` (150 dpi, `dot -Tpng -Gdpi=150`) |
| R model | `<abbr>_mrgsolve_model.R` |
| Shiny | `<abbr>_shiny_app.R` (or `shiny_app/app.R`) |
| References | `<abbr>_references.md` |

### 모델 품질 기준 (Model Quality Standards)
1. **기계론적 지도 (.dot)**: 100개 이상 노드, 최소 8개 서브그래프 클러스터, 약물 PK/PD 포함
2. **mrgsolve 모델 (.R)**: 
   - 최소 15개 ODE 구획 (약물 PK + 질환 PD)
   - 최소 5개 치료 시나리오
   - 주요 임상시험 데이터로 파라미터 보정 메모 포함
3. **Shiny 앱**: 최소 6개 탭 (환자 프로파일 · PK · PD 주요지표 · 임상 엔드포인트 · 시나리오 비교 · 바이오마커)
4. **참고문헌**: 최소 30개 PubMed 링크, 각 섹션별 분류

### README 업데이트 (README Update)

⚠️ **2026-07-02 사고 재발 방지 — 반드시 읽을 것.** 예전에는 표 아래에 질환별
상세 섹션을 추가했지만, 이것이 누적되어 파일이 수천 줄로 불어났고 결국 표의
일부 행에 남은 인용문헌·볼드 서식이 `**` 짝이 안 맞는 상태로 굳어지면서
**표 렌더링 자체가 깨지는 사고**가 발생했습니다. 아래 규칙은 그 재발을
막기 위한 것이며 예외 없이 지켜야 합니다.

- **표 아래에 질환 상세 섹션을 절대 추가하지 않는다.** README는 "①소개
  섹션 → ②모델 갤러리 표 → ③면책/참고자료/라이선스"로 끝나야 합니다.
  질환에 대한 상세 설명은 오직 `<disease-dir>/README.md` (디렉토리별
  README)에만 작성합니다.
- 표에 새 행을 추가할 때 (마지막 행 뒤, 번호는 항상 연속):
  - **요약 셀은 한 문장**으로만 작성합니다. 대략 190자 이내.
  - **볼드(`**`)·이탤릭(`*`)을 요약 셀 안에 쓰지 않습니다.** (모델명 셀의
    `[**한글명**...]`은 예외 — 그건 항상 정확히 한 쌍입니다.) 셀 하나의
    `**` 개수가 홀수가 되면 그 행 이후의 표 전체 렌더링이 깨집니다.
  - **인용문헌·저자명·연도·저널명·PMID를 요약 셀에 넣지 않습니다.**
    (`(Author 2020 Journal[PMID 12345])` 같은 것 전부 금지 — 근거는
    `<abbr>_references.md`에 있습니다.)
  - 영문 부제(`<sub>...</sub>`)는 `English Name · ABBR` 형태만 유지하고,
    기전·약물 키워드를 줄줄이 덧붙이지 않습니다.
  - 표 이미지: `<a href="path/to/svg"><img src="path/to/png" width="190" alt="ABBR"></a>` 형식 (기존 행과 동일한 패턴을 그대로 따라 작성).
- **행을 추가한 뒤, 커밋 전에 반드시 정리 스크립트를 실행합니다**:
  ```bash
  python3 scripts/fix_readme_table.py
  ```
  이 스크립트가 방금 추가한 행의 서식(볼드 제거·인용 제거·길이 제한·부제
  정리·카테고리 정규화)을 자동으로 정리하고, 표 번호를 다시 매기고,
  분류별 개수와 "N개 모델" 문구를 전부 갱신하고, 꼬리 섹션(면책/참고자료/
  라이선스)이 없으면 복원합니다. 이 스크립트는 표 형식만 정규화할 뿐 다른
  행의 내용(질환 설명)은 건드리지 않습니다.
- 위 스크립트를 실행한 뒤 `python3 scripts/fix_readme_table.py --check`
  가 `PASS`를 출력하는지 확인하고 커밋합니다. `FAIL`이 나오면 (표 번호가
  끊김·열 개수 불일치·`**` 짝 불일치·깨진 링크·꼬리 섹션 누락 중 하나)
  원인을 고친 뒤 다시 확인합니다.
- 기존 다른 행의 내용은 수정 금지 (정리 스크립트에 의한 서식 정규화는 예외).

### 커밋 & 푸시 (Commit & Push)
```bash
git add -A
git commit -m "Add <Disease> QSP model: mechanistic map, mrgsolve ODE, Shiny app, references"
git push -u origin HEAD
```
PR은 생성하지 말고 main에 직접 병합하거나, 사용자 지시에 따라 진행.

---

## 기술 스택 (Technology Stack)

| 도구 | 용도 |
|------|------|
| **Graphviz** (`dot`) | 기계론적 지도 렌더링 (.dot → .svg/.png) |
| **mrgsolve** (R) | ODE 기반 PK/PD 모델 |
| **Shiny** (R) | 인터랙티브 대시보드 |
| **Claude Code Routine** | 매일 자동 모델 생성 및 커밋 |

---

## 질환 목록 (Disease List)

세션마다 아래 목록에서 선택. 이미 추가된 것은 README 표 참조.

### 자가면역질환
류마티스 관절염 · 전신 홍반 루푸스 · 쇼그렌 증후군 · 전신경화증 · 다발성 근염 · 피부근염 · 베체트병 · 강직성 척추염 · 건선성 관절염 · 반응성 관절염 · 혼합결합조직병(MCTD) · 항인지질항체 증후군 · 육아종증 다발혈관염(GPA) · 호산구 육아종증 다발혈관염(EGPA) · 현미경적 다발혈관염(MPA) · 결절성 다발동맥염(PAN) · 다카야스 동맥염 · 거대세포 동맥염 · IgA 혈관염 · 재발성 다발연골염 · 크론병 · 궤양성 대장염 · 자가면역 간염 · 원발성 담즙성 담관염(PBC) · 원발성 경화성 담관염(PSC) · 셀리악병 · 자가면역 췌장염 · 악성 빈혈 · 제1형 당뇨병 · 하시모토 갑상선염 · 그레이브스병 · 애디슨병 · 다발성 경화증 · 중증 근무력증 · 시신경척수염(NMO) · 길랭-바레 증후군 · CIDP · 자가면역 뇌염 · ITP · AIHA · 에반스 증후군 · 심상성 천포창 · 수포성 유천포창 · 백반증 · 원형 탈모증 · 굿파스처 증후군 · **IgA 신병증** ✓ · 성인형 스틸병

### 만성질환
제2형 당뇨병 · 이상지질혈증 · 대사 증후군 · 비만 · 골다공증 · 파젯병 · 통풍 · 가성통풍 · 다낭성 난소 증후군(PCOS) · 말단비대증 · 원발성 부갑상선 기능 항진증 · 만성 갑상선 기능 저하증 · **본태성 고혈압** ✓ · 안정형 협심증 · 심부전(HFrEF) · 심부전(HFpEF) · 심방세동 · 말초동맥질환(PAD) · 만성 정맥 부전 · 복부 대동맥류 · 비후성 심근병증(HCM) · 확장성 심근병증(DCM) · COPD · 기관지 천식 · 기관지 확장증 · **폐동맥 고혈압(PAH)** ✓ · 특발성 폐섬유화증(IPF) · 진폐증 · 폐쇄성 수면 무호흡증(OSA) · 사르코이드증 · GERD · 만성 위염 · 소화성 궤양 · NAFLD · 간경변증 · 만성 B형 간염 · 만성 췌장염 · 과민성 장증후군(IBS) · 게실병 · 담석증 · 만성 신부전(CKD) · 미세변화 신증후군 · FSGS · 막성 신병증 · ADPKD · 양성 전립선 비대증(BPH) · 과민성 방광 · 요로결석

(✓ = 완성됨 / completed)
