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
- 표의 **마지막 행**에 새 모델 추가 (날짜순 유지)
- 표 이미지: `[![Alt](path/to/png)](path/to/svg)` 형식으로 썸네일 링크
- README 맨 아래에 해당 질환 상세 섹션 추가
- 기존 내용 수정 금지

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
제2형 당뇨병 · 이상지질혈증 · 대사 증후군 · 비만 · 골다공증 · 파젯병 · 통풍 · 가성통풍 · 다낭성 난소 증후군(PCOS) · 말단비대증 · 원발성 부갑상선 기능 항진증 · 만성 갑상선 기능 저하증 · 본태성 고혈압 · 안정형 협심증 · 심부전(HFrEF) · 심부전(HFpEF) · 심방세동 · 말초동맥질환(PAD) · 만성 정맥 부전 · 복부 대동맥류 · 비후성 심근병증(HCM) · 확장성 심근병증(DCM) · COPD · 기관지 천식 · 기관지 확장증 · **폐동맥 고혈압(PAH)** ✓ · 특발성 폐섬유화증(IPF) · 진폐증 · 폐쇄성 수면 무호흡증(OSA) · 사르코이드증 · GERD · 만성 위염 · 소화성 궤양 · NAFLD · 간경변증 · 만성 B형 간염 · 만성 췌장염 · 과민성 장증후군(IBS) · 게실병 · 담석증 · 만성 신부전(CKD) · 미세변화 신증후군 · FSGS · 막성 신병증 · ADPKD · 양성 전립선 비대증(BPH) · 과민성 방광 · 요로결석

(✓ = 완성됨 / completed)
