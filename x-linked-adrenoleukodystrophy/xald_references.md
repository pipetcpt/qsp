# X-연관 부신백질형성장애 (X-ALD) — 참고문헌
### X-linked Adrenoleukodystrophy · References for the QSP model

**PMID 검증 방법.** 아래 모든 항목의 PMID는 NCBI E-utilities(`esearch` + `esummary`)로
조회하여 **저자·저널·연도·제목이 실제로 일치하는지 한 건씩 확인**했습니다. 기억에
의존해 적은 PMID는 한 건도 없습니다(초안 작성 중 확인 없이 적어본 두 개의 PMID는
전혀 다른 논문 — 하나는 *Cryptococcus* 페로몬 유전자, 하나는 TIPS 간관류 CT — 이었고,
그 시점부터 전수 조회로 전환했습니다). 검증되지 않은 인용은 이 문서에 넣지 않았습니다.

각 항목의 **[모델]** 표시는 그 문헌이 `xald_mrgsolve_model.R`의 어느 파라미터·구조·
앵커를 뒷받침하는지를 가리킵니다.

---

## 1. 유전자 · ALDP 단백 · 변이 스펙트럼

1. Mosser J, et al. Putative X-linked adrenoleukodystrophy gene shares unexpected homology with ABC transporters. **Nature** 1993.
   <https://pubmed.ncbi.nlm.nih.gov/8441467/>
   — *ABCD1* 원인 유전자 동정. [모델] `ABCD1`→`ALDP` 구조, `MUTRES`(잔존 수송 기능).

2. Kemp S, et al. ABCD1 mutations and the X-linked adrenoleukodystrophy mutation database: role in diagnosis and clinical correlations. **Hum Mutat** 2001.
   <https://pubmed.ncbi.nlm.nih.gov/11748843/>
   — 변이 스펙트럼과 **유전형–표현형 상관의 부재**. [모델] `MUTRES`를 표현형 결정에서
   분리하고 `SUSC`에 결정권을 준 근거.

3. Fu H, et al. A Novel Missense Variant of the ABCD1 Gene in X-Linked Adrenoleukodystrophy in a Chinese Family. **Mol Genet Genomic Med** 2025.
   <https://pubmed.ncbi.nlm.nih.gov/41041882/>
   — missense 변이의 잔존 기능 가변성. [모델] `MUTRES` > 0 설정 가능성.

4. Pujol A, et al. Functional overlap between ABCD1 (ALD) and ABCD2 (ALDR) transporters: a therapeutic target for X-adrenoleukodystrophy. **Hum Mol Genet** 2004.
   <https://pubmed.ncbi.nlm.nih.gov/15489218/>
   — ABCD2의 부분적 기능 중복. [모델] `ABCD2F`, `FABCD2`, 그리고 ABCD2가 **기저에서는
   보상하지 않고 유도되었을 때만** 기여하도록 한 구조(`rect(ABCD2F-1)`).

---

## 2. VLCFA 대사 — 신장(elongation) vs 퍼옥시좀 분해

5. Schackmann MJ, et al. Enzymatic characterization of ELOVL1, a key enzyme in very long-chain fatty acid synthesis. **Biochim Biophys Acta** 2015.
   <https://pubmed.ncbi.nlm.nih.gov/25499606/>
   — ELOVL1이 C26:0 합성의 율속 효소. [모델] `KELONG`, `ELOVL1` 경쟁적 저해 구조.

6. Kemp S, et al. The complex lipidome as a driver of tissue-specific pathology in adrenoleukodystrophy. **Mol Neurodegener** 2026.
   <https://pubmed.ncbi.nlm.nih.gov/42469918/>
   — 조직별 지질종 차이가 조직별 병리를 만든다. [모델] 부신/뇌/척수/고환을 **별개
   구획**으로 두고 각 조직에 고유 유입·처리 상수를 준 근거.

7. Hama K, et al. Phosphatidylcholine with a C26:0 moiety, a precursor of a diagnostic marker for X-ALD, is synthesized by LPLAT10/LPEAT2. **J Lipid Res** 2026.
   <https://pubmed.ncbi.nlm.nih.gov/41478358/>
   — C26:0-lysoPC가 C26:0-PC 리모델링에서 생성됨. [모델] `KLPC`와 **초선형 지수
   `HLPC` = 1.7** — lysoPC가 C26:0 자체보다 민감한 표지자인 이유.

8. McGuinness MC, et al. Evaluation of pharmacological induction of fatty acid beta-oxidation in X-linked adrenoleukodystrophy. **Mol Genet Metab** 2001.
   <https://pubmed.ncbi.nlm.nih.gov/11592822/>
   — β-산화 약리적 유도의 한계. [모델] `EIND2` = 0.35(작게) 설정.

9. Come JH, et al. Discovery and Optimization of Pyrazole Amides as Inhibitors of ELOVL1. **J Med Chem** 2021.
   <https://pubmed.ncbi.nlm.nih.gov/34748351/>
   — 직접 ELOVL1 저해제. [모델] `ELOVLID`.

10. Holley S, et al. A comprehensive discovery platform for ELOVL1 small-molecule inhibitors targeting very long-chain fatty acid synthesis in adrenoleukodystrophy. **J Biol Chem** 2026.
    <https://pubmed.ncbi.nlm.nih.gov/42103229/>
    — 기질 감소 요법(substrate reduction) 전략. [모델] `ELOVLID` 시나리오.

---

## 3. 신생아 선별 · 진단 표지자

11. Tang C, et al. A pilot study of newborn screening for X-linked adrenoleukodystrophy based on liquid chromatography-tandem mass spectrometry. **Clin Chim Acta** 2024.
    <https://pubmed.ncbi.nlm.nih.gov/37977233/>
    — DBS C26:0-lysoPC 선별. [모델] 앵커 A5·A6(출생 시 lysoPC 상승 배수).

12. Chen HA, et al. High incidence of null variants identified from newborn screening of X-linked adrenoleukodystrophy in Taiwan. **Mol Genet Metab Rep** 2022.
    <https://pubmed.ncbi.nlm.nih.gov/36046390/>
    — 선별에서 확인되는 null 변이 비율. [모델] `MUTRES` 기본값 0.02.

13. Hricovec M, et al. A Qualitative Study on Parental Experiences with Genetic Counseling After a Positive Newborn Screen for Recently Added Conditions on the RUSP. **Int J Neonatal Screen** 2025.
    <https://pubmed.ncbi.nlm.nih.gov/41283363/>
    — RUSP 등재 이후 실제 운용. [모델] 선별→감시 경로(지도의 cluster 14).

14. Pierpont EI, et al. Neurofilament light chain as a prognostic marker in cerebral adrenoleukodystrophy. **Brain** 2026.
    <https://pubmed.ncbi.nlm.nih.gov/40988568/>
    — 축삭 손상 혈액 표지자. [모델] `AXB`(축삭 완전성)를 상태변수로 둔 근거.

---

## 4. 표현형 · 자연사 — "같은 변이, 다른 병"

15. Engelen M, et al. X-linked adrenoleukodystrophy (X-ALD): clinical presentation and guidelines for diagnosis, follow-up and management. **Orphanet J Rare Dis** 2012.
    <https://pubmed.ncbi.nlm.nih.gov/22889154/>
    — 표현형 분포·경과·관리 표준. [모델] 앵커 A8·A11, 시나리오 02–06.

16. Lier J, et al. Natural History of Clinical Phenotypes and Their Biochemical Correlates in Adult X-Linked Adrenoleukodystrophy. **J Inherit Metab Dis** 2026.
    <https://pubmed.ncbi.nlm.nih.gov/41853938/>
    — 생화학 지표와 임상 표현형의 관계. [모델] **혈장 C26:0가 표현형을 예측하지 못한다**는
    이 모델의 중심 주장(시나리오 20 스윕)의 직접 근거.

17. Mallack EJ, et al. Clinical and radiographic course of arrested cerebral adrenoleukodystrophy. **Neurology** 2020.
    <https://pubmed.ncbi.nlm.nih.gov/32482842/>
    — 저절로 정지하는 뇌형이 존재한다. [모델] `MGP`의 **이중안정성**(자연적 소멸 가능성).

18. Carlson AM, et al. Five men with arresting and relapsing cerebral adrenoleukodystrophy. **J Neurol** 2021.
    <https://pubmed.ncbi.nlm.nih.gov/32995952/>
    — 정지 후 재발. [모델] 스위치가 임계 근방에서 왕복할 수 있는 구조.

19. Gupta AO, et al. Differential outcomes for frontal versus posterior demyelination in childhood cerebral adrenoleukodystrophy. **J Inherit Metab Dis** 2021.
    <https://pubmed.ncbi.nlm.nih.gov/34499753/>
    — 전두엽형 vs 후두-두정형. [모델] 지도의 `SPLENIUM`/`FRONTAL` 노드.

20. Jiang W, et al. Initial frontal lobe involvement in adult cerebral X-linked adrenoleukodystrophy. **Acta Neurol Belg** 2023.
    <https://pubmed.ncbi.nlm.nih.gov/37247117/>
    — 성인 뇌형의 전두엽 우세. [모델] 성인 전환 아형.

21. Terrone G, et al. Cerebral X-Linked Adrenoleukodystrophy Associated with Hemophilia A: a case report. **Horm Res Paediatr** 2026.
    <https://pubmed.ncbi.nlm.nih.gov/42328992/>
    — Xq28 인접 유전자 동시 침범. [모델] 유전좌위 맥락(지도 cluster 1).

---

## 5. 미세아교세포 — 스위치의 방아쇠 (모델의 핵심)

22. Eichler FS, et al. Is microglial apoptosis an early pathogenic change in cerebral X-linked adrenoleukodystrophy? **Ann Neurol** 2008.
    <https://pubmed.ncbi.nlm.nih.gov/18571777/>
    — **미세아교세포 세포자멸사가 탈수초보다 먼저 일어난다.** [모델] 점화 게이트를
    "잔해"가 아니라 **미세아교세포 결손 `MGDEF`** 로 잡은 직접 근거. (잔해로 게이트를
    걸었던 첫 판은 스위치가 자기 시동을 걸 수 없었습니다 — 모델 파일에 기록.)

23. Bergner CG, et al. Microglia damage precedes major myelin breakdown in X-linked adrenoleukodystrophy and metachromatic leukodystrophy. **Glia** 2019.
    <https://pubmed.ncbi.nlm.nih.gov/30980503/>
    — 부검 기반 확증. 병변 주변부 미세아교세포 고갈. [모델] `KMGINF`(병변 내 추가 소실),
    `MGLDEP` 노드.

24. Lombroso SI, et al. Microglia replacement by peripheral delivery of CSF1R inhibitor-resistant hematopoietic cells. **Mol Ther** 2026.
    <https://pubmed.ncbi.nlm.nih.gov/41232525/>
    — 말초 조혈세포에 의한 미세아교세포 치환의 동역학. [모델] `KMGREP` = 0.0020/day
    (**치료 지체 창의 원천**), 앵커 A18.

25. Berghoff SA, et al. Microglia facilitate repair of demyelinated lesions via post-squalene sterol synthesis. **Nat Neurosci** 2021.
    <https://pubmed.ncbi.nlm.nih.gov/33349711/>
    — 미세아교세포의 **수복** 기능. [모델] `PHEFFN`/`PHEFFC`(교정 세포의 우월한 처리
    능력)와 `KMYER`(재수초화).

26. Parasar P, et al. iPSC-Derived Astrocytes to Model Neuroinflammatory and Metabolic Responses in X-linked Adrenoleukodystrophy. **J Biotechnol Biomed** 2023.
    <https://pubmed.ncbi.nlm.nih.gov/38077449/>
    — 성상교세포의 자율적 염증 반응. [모델] `AST`, `FAST`(CCL2 증폭 기여).

---

## 6. 신경염증 · 혈액뇌장벽

27. McGuinness MC, et al. Cerebral inflammation in X-linked adrenoleukodystrophy. **Arch Immunol Ther Exp** 1999.
    <https://pubmed.ncbi.nlm.nih.gov/10604233/>
    — TNF-α 중심의 병변 내 염증. [모델] `TNF`, `KOLGDT`, `KBBBD`.

28. Selmaj K, et al. Identification of lymphotoxin and tumor necrosis factor in multiple sclerosis lesions. **J Clin Invest** 1991.
    <https://pubmed.ncbi.nlm.nih.gov/1999503/>
    — TNF의 희소교세포 독성(탈수초 병변 일반). [모델] `KOLGDT` 크기의 참조.

29. Rodríguez-Pascau L, et al. The brain penetrant PPARγ agonist leriglitazone restores multiple altered pathways in models of X-linked adrenoleukodystrophy. **Sci Transl Med** 2021.
    <https://pubmed.ncbi.nlm.nih.gov/34078742/>
    — PPARγ 활성화가 염증·미토콘드리아·장벽을 동시에 건드린다. [모델] `KSPPAR`,
    `KPPBBB`, `PPARA` 경로.

---

## 7. 산화 스트레스 · 미토콘드리아 · Nrf2

30. Fourcade S, et al. Oxidative stress, mitochondrial and proteostasis malfunction in adrenoleukodystrophy: A paradigm for axonal degeneration. **Free Radic Biol Med** 2015.
    <https://pubmed.ncbi.nlm.nih.gov/26073123/>
    — 산화 스트레스 → 축삭 변성 축. [모델] `ROSB`, `ATPD`, `KAXDSC`(AMN 축삭병증).

31. Petrillo S, et al. Antioxidant Response in Human X-Linked Adrenoleukodystrophy Fibroblasts. **Antioxidants (Basel)** 2022.
    <https://pubmed.ncbi.nlm.nih.gov/36358497/>
    — X-ALD 세포의 항산화 반응 감쇠. [모델] `NRF2` 저하 구조.

32. Seminotti B, et al. Nuclear Factor Erythroid-2-Related Factor 2 Signaling in the Neuropathophysiology of Inherited Metabolic Disorders. **Front Cell Neurosci** 2021.
    <https://pubmed.ncbi.nlm.nih.gov/34955754/>
    — Nrf2 경로. [모델] `NRF`, `KNRF`, `DMFD`.

33. Kartha RV, et al. Mechanisms of Antioxidant Induction with High-Dose N-Acetylcysteine in Childhood Cerebral Adrenoleukodystrophy. **CNS Drugs** 2015.
    <https://pubmed.ncbi.nlm.nih.gov/26670322/>
    — 고용량 NAC. [모델] `NACD`, `ENAC`, 시나리오 21(**실패로 기록**).

34. Dann J, et al. Monomethyl Fumarate Modulates Iron Metabolism and Mitochondrial Function in Microglia. **Cell Mol Neurobiol** 2026.
    <https://pubmed.ncbi.nlm.nih.gov/42474538/>
    — 푸마레이트의 미세아교세포 작용. [모델] `DMFD`.

---

## 8. 부신 · 고환 축

35. Sheng Y, et al. Lessons from the gonadotropin-regulated long chain acyl-CoA synthetase (GR-LACS) null mouse model: a role in steroidogenesis. **J Steroid Biochem Mol Biol** 2009.
    <https://pubmed.ncbi.nlm.nih.gov/19167491/>
    — 지방산 활성화 효소와 스테로이드 생성의 연결. [모델] `STAR`/`STEROID` 축,
    `KIMC2R`(막 미세도메인 → MC2R 커플링 저하).

*(부신부족 유병률·발병 연령의 1차 근거는 위 15번 Engelen 2012 가이드라인과 16번
Lier 2026 자연사 연구입니다 — 앵커 A8·A9·A10.)*

---

## 9. 척수형 (AMN) · 여성 보유자

36. van de Stadt SIW, et al. Spinal cord atrophy as a measure of severity of myelopathy in adrenoleukodystrophy. **J Inherit Metab Dis** 2020.
    <https://pubmed.ncbi.nlm.nih.gov/32077106/>
    — 척수 위축의 정량화. [모델] `CSTI`/`DCI`, `LENCST`/`LENDC`(길이 의존 취약성),
    앵커 A15(EDSS).

37. Engelen M, et al. X-linked adrenoleukodystrophy in women: a cross-sectional cohort study. **Brain** 2014.
    <https://pubmed.ncbi.nlm.nih.gov/24480483/>
    — 여성 보유자의 척수형 유병률. [모델] `FMOS`, 시나리오 05·06.

38. Wang Z, et al. Familial skewed X chromosome inactivation in adrenoleukodystrophy manifesting heterozygotes from a Chinese pedigree. **PLoS One** 2013.
    <https://pubmed.ncbi.nlm.nih.gov/23469258/>
    — X 불활성화 편향이 발현을 결정. [모델] `FMOS`↔`XINACT`, `CRSC`(세포 간 상보).

39. Schäfer L, et al. Self-reported quality of life in symptomatic and asymptomatic women with X-linked adrenoleukodystrophy. **Brain Behav** 2023.
    <https://pubmed.ncbi.nlm.nih.gov/36748403/>
    — 여성 보유자의 임상 부담. [모델] 시나리오 05·06의 임상적 의미.

---

## 10. 로렌조 오일 — 표지자는 좋아지고 병은 진행한다

40. Aubourg P, et al. A two-year trial of oleic and erucic acids ("Lorenzo's oil") as treatment for adrenomyeloneuropathy. **N Engl J Med** 1993.
    <https://pubmed.ncbi.nlm.nih.gov/8350883/>
    — **혈장 VLCFA는 정상화되고 신경학적 진행은 막지 못했다.** [모델] 앵커 A16(−45%)와
    **A17(Loes 감소율 0.0%)** — 시나리오 07 vs 08.

41. Moser HW, et al. Follow-up of 89 asymptomatic patients with adrenoleukodystrophy treated with Lorenzo's oil. **Arch Neurol** 2005.
    <https://pubmed.ncbi.nlm.nih.gov/16009761/>
    — 무증상 단계 장기 추적. [모델] 시나리오 09(AMN + 식이제한), `DIETSC`.

42. Terre'Blanche G, et al. Treatment of an adrenomyeloneuropathy patient with Lorenzo's oil and supplementation with docosahexaenoic acid. **Lipids Health Dis** 2011.
    <https://pubmed.ncbi.nlm.nih.gov/21871076/>
    — 개별 환자의 혈장 지방산 반응. [모델] `KIERU`, 에루크산 정상상태 농도(약 150 μmol/L).

43. Hartley MD, et al. A Thyroid Hormone-Based Strategy for Correcting the Biochemical Abnormality in X-Linked Adrenoleukodystrophy. **Endocrinology** 2017.
    <https://pubmed.ncbi.nlm.nih.gov/28200172/>
    — CNS 침투 티로미메틱(소베티롬 계열)이 **뇌** VLCFA를 낮춘다. [모델] `SOBB`,
    `ESOBBR`, 시나리오 22 — 로렌조 오일과 대비되는 유일한 대사 경로.

---

## 11. 동종 HSCT · 자가 유전자치료 · 안전성

44. Peters C, et al. Cerebral X-linked adrenoleukodystrophy: the international hematopoietic cell transplantation experience from 1982 to 1999. **Blood** 2004.
    <https://pubmed.ncbi.nlm.nih.gov/15073029/>
    — 국제 코호트. **치료 창(Loes·NFS)의 중요성.** [모델] `INWINDOW`, 시나리오 10 vs 12.

45. Miller WP, et al. Outcomes after allogeneic hematopoietic cell transplantation for childhood cerebral adrenoleukodystrophy. **Blood** 2011.
    <https://pubmed.ncbi.nlm.nih.gov/21586746/>
    — 최대 단일기관 코호트. 이식 후에도 일정 기간 진행. [모델] **`LAG` 지체 창**, 앵커 A18.

46. Raymond GV, et al. Survival and Functional Outcomes in Boys with Cerebral Adrenoleukodystrophy with and without Hematopoietic Stem Cell Transplantation. **Biol Blood Marrow Transplant** 2019.
    <https://pubmed.ncbi.nlm.nih.gov/30292747/>
    — 이식군 vs 비이식군 직접 비교. [모델] 시나리오 10 vs **짝지은 대조군 11**.

47. Cartier N, et al. Hematopoietic stem cell gene therapy with a lentiviral vector in X-linked adrenoleukodystrophy. **Science** 2009.
    <https://pubmed.ncbi.nlm.nih.gov/19892975/>
    — 자가 유전자치료 최초 보고. [모델] `LENTID`, `VCNS`, `MONOC`.

48. Eichler F, et al. Lentiviral Gene Therapy for Cerebral Adrenoleukodystrophy. **N Engl J Med** 2024.
    <https://pubmed.ncbi.nlm.nih.gov/39383459/>
    — eli-cel(elivaldogene autotemcel) 장기 성적, 무MFD 생존. [모델] 앵커 A19·A22.

49. Duncan CN, et al. Hematologic Cancer after Gene Therapy for Cerebral Adrenoleukodystrophy. **N Engl J Med** 2024.
    <https://pubmed.ncbi.nlm.nih.gov/39383458/>
    — **삽입 돌연변이에 의한 MDS/AML** — 실측 안전성 신호. [모델] `CLN`, `HMDS`, `PMDS`,
    `KCLNB`(부설판 노출 증폭), 시나리오 14 vs 15.

50. Siddique N, et al. Comparison of Methods to Estimate Busulfan Exposure in Pediatric Hematopoietic Stem Cell Transplant Recipients. **Ther Drug Monit** 2026.
    <https://pubmed.ncbi.nlm.nih.gov/42297573/>
    — 부설판 노출(AUC) 목표 관리. [모델] `KEBUS`, `AUCB`, 앵커 A23(90 mg·h/L).

---

## 12. QSP 방법론

51. Elmokadem A, et al. Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial. **CPT Pharmacometrics Syst Pharmacol** 2019.
    <https://pubmed.ncbi.nlm.nih.gov/31652028/>
    — 이 저장소가 쓰는 도구. [모델] 파일 전체 구조.

52. Dogné JM, et al. From animal testing to model-informed drug development: building on ICH M15 and EMA initiatives. **Front Pharmacol** 2026.
    <https://pubmed.ncbi.nlm.nih.gov/42499501/>
    — MIDD 규제 맥락(ICH M15). [모델] 모델의 용도 규정.

---

## 13. 파라미터 출처 요약 (LIT / FIT / STRUCT)

| 구분 | 뜻 | 예 |
|------|-----|-----|
| **LIT** | 문헌에서 직접 가져온 값 | 정상/환자 혈장 C26:0, 부설판 AUC 목표, 코르티솔 반감기 |
| **FIT** | `xald_calibrate()`의 좌표별 반복으로 앵커에 맞춘 값 | `KELONG`, `KOMEG`, `KLPC`, `KADMD`, `KIGN`, `KMGD50`, `KIERU`, `KMGREP`, `KAXDSC` |
| **STRUCT** | 직접 측정값이 없는 구조적 선택 (그렇다고 명시) | `KAMP50`/`HAMP`(증폭 포화), `KSPREAD`(병변 전선 확산), `HSUSC`(감수성 Hill), `FMGL`, `KBRSYN`/`KINBR` 분배비 |

**모델이 명시적으로 경험적인 관계 하나.** Loes 점수와 탈수초 분율의 관계
(`LOESMAX * hup(DEMYF, KLOES50, HLOES)`)는 기전식이 아니라 **판독 척도의 보정
곡선**입니다. Loes는 부위별 가중 합산 점수이므로 부피 분율과의 관계가 선형이 아니며,
저값 구간에서 평평해야 순수 AMN 환자가 뇌 점수를 얻지 않습니다. 이 한 곳은 숨기지
않고 여기에 표시합니다. (지수 포화 곡선을 먼저 써봤을 때 탈수초 1%에서 Loes가 1점을
넘어 모든 AMN 환자가 뇌 병변 점수를 받았고, 그래서 Hill로 교체했습니다.)

---

## 14. 이 문서에 넣지 않은 것

- **미검증 PMID**: 한 건도 없습니다.
- **레리글리타존 ADVANCE 시험의 세부 사건 수**: 1차 평가변수 미충족은 확실하나(위 5번,
  Köhler 2023, *Lancet Neurol* — PMID 36681445), 뇌 병변 발생 건수의 정확한 분모/분자를
  본문 확인 없이 옮기지 않았습니다. 모델은 이 시험을 **"부분적 지연, 예방은 아님"**
  으로만 반영했고(시나리오 18·19), 수치 앵커로 쓰지 않았습니다.
- **국내 유병률 통계**: 신뢰할 수 있는 1차 출처를 확인하지 못해 넣지 않았습니다.
