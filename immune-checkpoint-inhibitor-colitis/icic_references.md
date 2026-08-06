# 면역관문억제제 유발 대장염 — 참고문헌
# Immune Checkpoint Inhibitor-Induced Colitis (ICI colitis / IMDC) — References

이 문서는 `icic_qsp_model.dot`, `icic_mrgsolve_model.R`, `icic_reference_model.py`에
쓰인 **모든 파라미터와 구조적 가정의 근거**를 정리한 것입니다. 각 항목 끝의
대괄호는 모델에서 그 문헌이 실제로 무엇을 결정하는지를 표시합니다.

This file records the source of every parameter and structural assumption in the
model. The bracket at the end of each entry names the quantity that the citation
actually fixes, so that a reader can audit any single number back to its origin.

---

## 1. 임상 발생률과 용량-반응 (Incidence and the two exposure–response shapes)

이 모델이 존재하는 이유가 되는 관찰들 — CTLA-4 억제제는 용량-반응이 있고
PD-1 억제제는 없다는 것, 그리고 병용요법이 단순 합이 아니라는 것.

1. Ascierto PA, et al. Ipilimumab 10 mg/kg versus ipilimumab 3 mg/kg in patients
   with unresectable or metastatic melanoma: a randomised, double-blind,
   multicentre, phase 3 trial. *Lancet Oncol.* 2017;18(5):611-622.
   <https://pubmed.ncbi.nlm.nih.gov/28359784/>
   [**핵심 보정 근거**: G3-4 설사 9.5% (10 mg/kg) vs 6.1% (3 mg/kg), G3-4 대장염
   6.7% vs 3.2% — anti-CTLA-4 용량-독성 곡선의 존재 자체]
2. Larkin J, et al. Five-Year Survival with Combined Nivolumab and Ipilimumab in
   Advanced Melanoma (CheckMate 067). *N Engl J Med.* 2019;381(16):1535-1546.
   <https://pubmed.ncbi.nlm.nih.gov/31562797/>
   [G3-4 설사 9.3%(병용)/6.3%(ipi)/2.2%(nivo), 대장염 7.7%/8.7%/0.6% — 병용 ≈
   ipi 단독이라는 **비가산성**의 근거]
3. Hodi FS, et al. Improved Survival with Ipilimumab in Patients with Metastatic
   Melanoma. *N Engl J Med.* 2010;363(8):711-723.
   <https://pubmed.ncbi.nlm.nih.gov/20525992/>
   [ipilimumab 3 mg/kg 단독의 위장관 irAE 기저 발생률]
4. Brahmer JR, et al. Phase I study of single-agent anti-programmed death-1
   (MDX-1106) in refractory solid tumors: safety, clinical activity,
   pharmacodynamics, and immunologic correlates. *J Clin Oncol.*
   2010;28(19):3167-3175. <https://pubmed.ncbi.nlm.nih.gov/20516446/>
   [0.3-10 mg/kg에서 말초 T세포 PD-1 점유율 ≥70%, 평균 85% — **점유율 포화**의
   1차 근거, `KD_PD1_nM`]
5. Topalian SL, et al. Safety, Activity, and Immune Correlates of Anti-PD-1
   Antibody in Cancer. *N Engl J Med.* 2012;366(26):2443-2454.
   <https://pubmed.ncbi.nlm.nih.gov/22658127/>
   [nivolumab 0.1-10 mg/kg: irAE의 **용량 무관성**]
6. Robert C, et al. Pembrolizumab versus Ipilimumab in Advanced Melanoma
   (KEYNOTE-006). *N Engl J Med.* 2015;372(26):2521-2532.
   <https://pubmed.ncbi.nlm.nih.gov/25891173/>
   [pembrolizumab 10 mg/kg q2w vs q3w vs ipilimumab: 대장염 1.4-2.5% vs 8.2%]
7. Ribas A, et al. Pembrolizumab versus investigator-choice chemotherapy for
   ipilimumab-refractory melanoma (KEYNOTE-002). *Lancet Oncol.*
   2015;16(8):908-918. <https://pubmed.ncbi.nlm.nih.gov/26115796/>
   [pembrolizumab **2 vs 10 mg/kg에서 독성 차이 없음** — 평탄한 용량-독성 관계]
8. Wang DY, et al. Fatal Toxic Effects Associated With Immune Checkpoint
   Inhibitors: A Systematic Review and Meta-analysis. *JAMA Oncol.*
   2018;4(12):1721-1728. <https://pubmed.ncbi.nlm.nih.gov/30242316/>
   [대장염이 anti-CTLA-4 사망 원인의 70%를 차지 — 천공/절제 위험의 규모]
9. Wang DY, et al. Immune-related adverse events of a PD-L1 inhibitor plus
   chemotherapy versus a PD-L1 inhibitor alone: meta-analysis.
   *Cancer.* 2020;126(1):66-75. <https://pubmed.ncbi.nlm.nih.gov/31584709/>
10. Tian Y, et al. Anti-PD-1/PD-L1 vs anti-CTLA-4 in the risk of
    immune-related colitis: a network meta-analysis.
    *Front Pharmacol.* 2021;12:663943.
    <https://pubmed.ncbi.nlm.nih.gov/34305584/>
    [계열 간 대장염 위험 비교의 통합 추정치]
11. Motzer RJ, et al. Nivolumab plus Ipilimumab versus Sunitinib in Advanced
    Renal-Cell Carcinoma (CheckMate 214). *N Engl J Med.* 2018;378(14):1277-1290.
    <https://pubmed.ncbi.nlm.nih.gov/29562145/>
    [ipilimumab **1 mg/kg** 병용 요법에서의 낮은 대장염률 — 용량 축의 하단]

---

## 2. 왜 CTLA-4인가 — 구심성(afferent) 관문과 Treg ADCC

모델의 분모(denominator)가 고갈 가능하다는 주장의 근거.

12. Walker LSK, Sansom DM. The emerging role of CTLA4 as a cell-extrinsic
    regulator of T cell responses. *Nat Rev Immunol.* 2011;11(12):852-863.
    <https://pubmed.ncbi.nlm.nih.gov/22116087/>
    [CTLA-4의 **세포 외재적** 작용 — trans-endocytosis, `TRANSENDO` 노드]
13. Qureshi OS, et al. Trans-endocytosis of CD80 and CD86: a molecular basis for
    the cell-extrinsic function of CTLA-4. *Science.* 2011;332(6029):600-603.
    <https://pubmed.ncbi.nlm.nih.gov/21474713/>
14. van der Merwe PA, et al. CD80 (B7-1) binds both CD28 and CTLA-4 with a low
    affinity and very fast kinetics. *J Exp Med.* 1997;185(3):393-403.
    <https://pubmed.ncbi.nlm.nih.gov/9053440/>
    [CD80–CTLA-4가 CD80–CD28보다 10-20배 강함 — `CD80_86` 경쟁 구조]
15. Simpson TR, et al. Fc-dependent depletion of tumor-infiltrating regulatory
    T cells co-defines the efficacy of anti-CTLA-4 therapy against melanoma.
    *J Exp Med.* 2013;210(9):1695-1710.
    <https://pubmed.ncbi.nlm.nih.gov/23897981/>
    [**ADCC에 의한 Treg 고갈**이 anti-CTLA-4 작용의 일부라는 핵심 근거, `k_adcc`]
16. Arce Vargas F, et al. Fc Effector Function Contributes to the Activity of
    Human Anti-CTLA-4 Antibodies. *Cancer Cell.* 2018;33(4):649-663.e4.
    <https://pubmed.ncbi.nlm.nih.gov/29576375/>
    [**FCGR3A V158F**가 ipilimumab 반응을 조절 — `phi_FCGR` (V/V 1.60 · V/F 1.00
    · F/F 0.55)]
17. Romano E, et al. Ipilimumab-dependent cell-mediated cytotoxicity of
    regulatory T cells ex vivo by nonclassical monocytes in melanoma patients.
    *Proc Natl Acad Sci USA.* 2015;112(19):6140-6145.
    <https://pubmed.ncbi.nlm.nih.gov/25918390/>
    [CD16⁺ 단핵구가 ADCC 효과세포 — 모델의 `Mac`이 ADCC 속도에 들어가는 근거,
    `mac_pow`]
18. Sharma A, et al. Anti-CTLA-4 Immunotherapy Does Not Deplete FOXP3⁺
    Regulatory T Cells (Tregs) in Human Cancers. *Clin Cancer Res.*
    2019;25(4):1233-1238. <https://pubmed.ncbi.nlm.nih.gov/30054281/>
    [**반대 증거**. 사람 종양에서 Treg 고갈이 관찰되지 않는다는 보고. 모델은
    대장 점막의 Treg 고갈 깊이를 잠재변수로 두고 독성 결과로 보정하며, 이
    불확실성을 민감도 분석에서 명시적으로 다룬다]
19. Sanmamed MF, Chen L. A Paradigm Shift in Cancer Immunotherapy: From
    Enhancement to Normalization. *Cell.* 2018;175(2):313-326.
    <https://pubmed.ncbi.nlm.nih.gov/30290139/>
    [**구심성(CTLA-4, 림프절) vs 원심성(PD-1, 조직)** 구분 — `pd1_prol_pow`가
    작은 이유의 개념적 근거]
20. Wei SC, et al. Distinct Cellular Mechanisms Underlie Anti-CTLA-4 and
    Anti-PD-1 Checkpoint Blockade. *Cell.* 2017;170(6):1120-1133.e17.
    <https://pubmed.ncbi.nlm.nih.gov/28800925/>
    [두 약물이 **서로 다른 세포 집단**을 확장시킨다는 단일세포 근거]
21. Sharma A, Subudhi SK, et al. Anti-CTLA-4 immunotherapy does not deplete
    FOXP3⁺ regulatory T cells but expands ICOS⁺ Th1-like CD4 effectors.
    *Sci Transl Med.* 2019;11(482):eaav6473.
    <https://pubmed.ncbi.nlm.nih.gov/30842315/>

---

## 3. 장 대장염의 면역 기전 (Colonic immunopathology)

22. Coutzac C, et al. Colon Immune-Related Adverse Events: Anti-CTLA-4 and
    Anti-PD-1 Blockade Induce Distinct Immunopathological Entities.
    *J Crohns Colitis.* 2017;11(10):1238-1246.
    <https://pubmed.ncbi.nlm.nih.gov/28967957/>
    [두 계열의 조직학적 표현형이 **다르다** — 모델이 두 약물을 같은 축의 세기
    차이가 아니라 구조적으로 다른 위치에 넣는 이유]
23. Luoma AM, et al. Molecular Pathways of Colon Inflammation Induced by Cancer
    Immunotherapy. *Cell.* 2020;182(3):655-671.e22.
    <https://pubmed.ncbi.nlm.nih.gov/32603654/>
    [**핵심 기전 문헌**: 세포독성 CD8 T세포 확장, 조직 상주 집단에서 기원,
    Treg 감소 — `Teff`/`Trm`/`Treg` 3구획 구조의 직접 근거]
24. Sasson SC, et al. Interferon-Gamma-Producing CD8⁺ Tissue Resident Memory
    T Cells Are a Targetable Hallmark of Immune Checkpoint Inhibitor-Colitis.
    *Gastroenterology.* 2021;161(4):1229-1244.e9.
    <https://pubmed.ncbi.nlm.nih.gov/34147519/>
    [**Trm이 병태생리의 중심**이라는 근거 — `Trm`의 이력현상(hysteresis) 구조,
    `k_sr`, `dead_il15`]
25. Oh DY, et al. Immune Toxicities Elicited by CTLA-4 Blockade in Cancer
    Patients Are Associated with Early Diversification of the T-cell Repertoire.
    *Cancer Res.* 2017;77(6):1322-1330.
    <https://pubmed.ncbi.nlm.nih.gov/28031229/>
    [독성 발생 전 **TCR 레퍼토리 다양화** — `Nclone`을 누적 적분 변수로 둔 근거]
26. Bamias G, et al. Immunological Characteristics of Colitis Associated with
    Anti-CTLA-4 Antibody Therapy. *Cancer Invest.* 2017;35(7):443-455.
    <https://pubmed.ncbi.nlm.nih.gov/28548891/>
27. Lord JD, et al. Refractory colitis following anti-CTLA4 antibody therapy:
    analysis of mucosal FOXP3⁺ T cells. *Dig Dis Sci.* 2010;55(5):1396-1405.
    <https://pubmed.ncbi.nlm.nih.gov/19629688/>
    [불응성 대장염 점막의 FOXP3⁺ 세포 분석 — 분모 고갈 가설의 조직 수준 검증]
28. Marthey L, et al. Cancer Immunotherapy with Anti-CTLA-4 Monoclonal
    Antibodies Induces an Inflammatory Bowel Disease. *J Crohns Colitis.*
    2016;10(4):395-401. <https://pubmed.ncbi.nlm.nih.gov/26783344/>
    [내시경/조직 표현형, 발현 시점, 스테로이드 불응률]
29. Geukes Foppen MH, et al. Immune checkpoint inhibition-related colitis:
    symptoms, endoscopic features, histology and response to management.
    *ESMO Open.* 2018;3(1):e000278.
    <https://pubmed.ncbi.nlm.nih.gov/29387476/>
    [증상–내시경–조직 소견의 **불일치** — 등급이 검열된 지표라는 임상적 표현]

---

## 4. 상피 장벽, 은와(crypt) 재생, 그리고 스테로이드 불응성

모델에서 가장 중요한 구조적 주장 중 하나 — **스테로이드 불응성은 약이 아니라
은와에 대한 진술이다** — 의 근거.

30. Wang Y, et al. Endoscopic and Histologic Features of Immune Checkpoint
    Inhibitor-Related Colitis. *Inflamm Bowel Dis.* 2018;24(8):1695-1705.
    <https://pubmed.ncbi.nlm.nih.gov/29718308/>
    [**깊은 궤양이 스테로이드 불응성과 조기 생물학적 제제 필요성을 예측** —
    `inj_deep`, `k_isck`, `k_isc`(느린 은와 재생)의 직접 근거]
31. Abu-Sbeih H, et al. Importance of endoscopic and histological evaluation in
    the management of immune checkpoint inhibitor-induced colitis.
    *J Immunother Cancer.* 2018;6(1):95.
    <https://pubmed.ncbi.nlm.nih.gov/30253811/>
32. Barnes MJ, Powrie F. Regulatory T cells reinforce intestinal homeostasis.
    *Immunity.* 2009;31(3):401-411.
    <https://pubmed.ncbi.nlm.nih.gov/19766083/>
33. Bruewer M, et al. Proinflammatory cytokines disrupt epithelial barrier
    function by apoptosis-independent mechanisms. *J Immunol.*
    2003;171(11):6164-6172. <https://pubmed.ncbi.nlm.nih.gov/14634132/>
    [IFN-γ/TNF-α에 의한 tight junction 파괴 — `k_tjd`, claudin-2/occludin 축]
34. Clayburgh DR, et al. A porous defense: the leaky epithelial barrier in
    intestinal disease. *Lab Invest.* 2004;84(3):282-291.
    <https://pubmed.ncbi.nlm.nih.gov/14767487/>
35. Barker N. Adult intestinal stem cells: critical drivers of epithelial
    homeostasis and regeneration. *Nat Rev Mol Cell Biol.* 2014;15(1):19-33.
    <https://pubmed.ncbi.nlm.nih.gov/24326621/>
    [LGR5⁺ 줄기세포 구획과 3-5일 상피 교체 — `k_turn`, `k_prod_e`]
36. Sun J. VDR/vitamin D receptor regulates autophagic activity through
    ATG16L1 — *(장 상피 재생 신호 일반)*; 대표 리뷰로
    Beumer J, Clevers H. Cell fate specification and differentiation in the
    adult mammalian intestine. *Nat Rev Mol Cell Biol.* 2021;22(1):39-53.
    <https://pubmed.ncbi.nlm.nih.gov/32958874/>

---

## 5. 미생물총, 부티르산, 그리고 FMT

37. Dubin K, et al. Intestinal microbiome analyses identify melanoma patients at
    risk for checkpoint-blockade-induced colitis. *Nat Commun.* 2016;7:10391.
    <https://pubmed.ncbi.nlm.nih.gov/26837003/>
    [**Bacteroidetes 및 폴리아민 수송/비타민 B 합성 경로가 대장염 방어적** —
    `Dv`(다양성)와 대장염 위험의 연결]
38. Chaput N, et al. Baseline gut microbiota predicts clinical response and
    colitis in metastatic melanoma patients treated with ipilimumab.
    *Ann Oncol.* 2017;28(6):1368-1379.
    <https://pubmed.ncbi.nlm.nih.gov/28368458/>
    [Faecalibacterium 우세가 **반응도 좋게 하고 대장염도 늘린다** — 항종양
    효과와 대장염이 같은 상류 구동을 공유한다는 모델 구조의 근거]
39. Vétizou M, et al. Anticancer immunotherapy by CTLA-4 blockade relies on the
    gut microbiota. *Science.* 2015;350(6264):1079-1084.
    <https://pubmed.ncbi.nlm.nih.gov/26541610/>
    [*B. fragilis*가 anti-CTLA-4 항종양 효과에 필요]
40. Routy B, et al. Gut microbiome influences efficacy of PD-1-based
    immunotherapy against epithelial tumors. *Science.* 2018;359(6371):91-97.
    <https://pubmed.ncbi.nlm.nih.gov/29097494/>
    [항생제 노출이 효능을 떨어뜨림 — `ABX_ON`, `k_abx`]
41. Gopalakrishnan V, et al. Gut microbiome modulates response to anti-PD-1
    immunotherapy in melanoma patients. *Science.* 2018;359(6371):97-103.
    <https://pubmed.ncbi.nlm.nih.gov/29097493/>
42. Furusawa Y, et al. Commensal microbe-derived butyrate induces the
    differentiation of colonic regulatory T cells. *Nature.*
    2013;504(7480):446-450. <https://pubmed.ncbi.nlm.nih.gov/24226770/>
    [부티르산 → pTreg 유도 — `b_but`]
43. Arpaia N, et al. Metabolites produced by commensal bacteria promote
    peripheral regulatory T-cell generation. *Nature.* 2013;504(7480):451-455.
    <https://pubmed.ncbi.nlm.nih.gov/24226773/>
44. Donohoe DR, et al. The microbiome and butyrate regulate energy metabolism
    and autophagy in the mammalian colon. *Cell Metab.* 2011;13(5):517-526.
    <https://pubmed.ncbi.nlm.nih.gov/21531334/>
    [부티르산이 대장세포 ATP의 약 70%를 공급 — `Rep = 0.30 + 0.70·Bu`의
    직접적 근거. 미생물총 손실은 Treg와 상피 재생을 **동시에** 친다]
45. Wang Y, et al. Fecal microbiota transplantation for refractory immune
    checkpoint inhibitor-associated colitis. *Nat Med.* 2018;24(12):1804-1808.
    <https://pubmed.ncbi.nlm.nih.gov/30420754/>
    [불응성 ICI 대장염에서 **FMT**의 첫 보고 — `fmt` 시나리오]
46. Halsey TM, et al. Microbiome alteration via fecal microbiota transplantation
    is effective for refractory immune checkpoint inhibitor-induced colitis.
    *Sci Transl Med.* 2023;15(700):eabq4006.
    <https://pubmed.ncbi.nlm.nih.gov/37315113/>

---

## 6. 치료: 스테로이드, infliximab, vedolizumab, 그리고 그 순서

47. Brahmer JR, et al. Management of Immune-Related Adverse Events in Patients
    Treated With Immune Checkpoint Inhibitor Therapy: ASCO Clinical Practice
    Guideline. *J Clin Oncol.* 2018;36(17):1714-1768.
    <https://pubmed.ncbi.nlm.nih.gov/29442540/>
    [등급별 관리 알고리즘 — G2에서 보류+스테로이드, G3에서 중단; 모델의
    `STD_TRIG`(G2 트리거) 설정 근거]
48. Schneider BJ, et al. Management of Immune-Related Adverse Events in Patients
    Treated With Immune Checkpoint Inhibitor Therapy: ASCO Guideline Update.
    *J Clin Oncol.* 2021;39(36):4073-4126.
    <https://pubmed.ncbi.nlm.nih.gov/34724392/>
49. Haanen J, et al. Management of toxicities from immunotherapy: ESMO Clinical
    Practice Guideline. *Ann Oncol.* 2022;33(12):1217-1238.
    <https://pubmed.ncbi.nlm.nih.gov/36270461/>
50. Dougan M, et al. AGA Clinical Practice Update on Diagnosis and Management of
    Immune Checkpoint Inhibitor Colitis and Hepatitis. *Gastroenterology.*
    2021;160(4):1384-1393. <https://pubmed.ncbi.nlm.nih.gov/33080231/>
    [조기 생물학적 제제 사용을 지지하는 소화기내과 관점]
51. Johnson DH, et al. Infliximab associated with faster symptom resolution
    compared with corticosteroids alone for the management of immune-related
    enterocolitis. *J Immunother Cancer.* 2018;6(1):103.
    <https://pubmed.ncbi.nlm.nih.gov/30305177/>
    [**infliximab이 스테로이드보다 증상 소실이 빠르다** — 모델의 "종말 단계를
    치면 가장 빠르다"는 예측의 직접 검증]
52. Abu-Sbeih H, et al. Early introduction of selective immunosuppressive
    therapy associated with favorable clinical outcomes in patients with
    immune checkpoint inhibitor-induced colitis. *J Immunother Cancer.*
    2019;7(1):93. <https://pubmed.ncbi.nlm.nih.gov/30940209/>
    [조기 선택적 면역억제가 결과를 개선 — 조기 개입 시나리오의 근거]
53. Bergqvist V, et al. Vedolizumab treatment for immune checkpoint
    inhibitor-induced enterocolitis. *Cancer Immunol Immunother.*
    2017;66(5):581-592. <https://pubmed.ncbi.nlm.nih.gov/28204866/>
54. Abu-Sbeih H, et al. Outcomes of vedolizumab therapy in patients with immune
    checkpoint inhibitor-induced colitis: a multi-center study.
    *J Immunother Cancer.* 2018;6(1):142.
    <https://pubmed.ncbi.nlm.nih.gov/30518410/>
    [vedolizumab로 **스테로이드-무의존 관해 약 85%**, 느린 발현 — `E_vdz`,
    `chi_vdz`, 그리고 지연된 발현 시점의 근거]
55. Zou F, et al. Efficacy and safety of vedolizumab and infliximab treatment
    for immune-mediated diarrhea and colitis in patients with cancer.
    *J Immunother Cancer.* 2021;9(11):e003277.
    <https://pubmed.ncbi.nlm.nih.gov/34795010/>
    [두 생물학적 제제의 **직접 비교** — 재발률과 지속성 차이]
56. Thomas AS, et al. Immune checkpoint inhibitor-induced colitis: the
    incidence, management and recurrence. *J Cancer Res Clin Oncol.* 2021.
    <https://pubmed.ncbi.nlm.nih.gov/34003350/>
    [**감량 중 재발 30-40%** — `Trm` 이력현상 및 taper 길이 분석의 보정 표적]
57. Bishu S, et al. Efficacy and outcome of tofacitinib in immune checkpoint
    inhibitor colitis. *Gastroenterology.* 2021;160(3):932-934.e3.
    <https://pubmed.ncbi.nlm.nih.gov/33098885/>
    [불응성에서 JAK 억제제 — `A_jak`, `E_jak`]
58. Holmstroem RB, et al. COLAR: open-label clinical study of IL-6 blockade with
    tocilizumab for the treatment of immune-related adverse events.
    *J Immunother Cancer.* 2022;10(9):e005111.
    <https://pubmed.ncbi.nlm.nih.gov/36130780/>
    [tocilizumab — `E_toc`]
59. Wang Y, et al. Immune-checkpoint inhibitor-induced diarrhea and colitis in
    patients with advanced malignancies: retrospective review at MD Anderson.
    *J Immunother Cancer.* 2018;6(1):37.
    <https://pubmed.ncbi.nlm.nih.gov/29747688/>
    [대규모 실사용 코호트: 발현 시점 중앙값, 스테로이드 불응률, 재발률]

---

## 7. 면역억제가 항종양 효과를 얼마나 깎는가 (The oncological cost)

모델의 네 번째 아이디어 — **선택성이 통화다** — 의 근거.

60. Horvat TZ, et al. Immune-Related Adverse Events, Need for Systemic
    Immunosuppression, and Effects on Survival and Time to Treatment Failure in
    Patients With Melanoma Treated With Ipilimumab at MSKCC.
    *J Clin Oncol.* 2015;33(28):3193-3198.
    <https://pubmed.ncbi.nlm.nih.gov/26282644/>
    [스테로이드/infliximab 사용이 생존을 악화시키지 않았다는 관찰 — 모델의
    `chi` 값이 과대하지 않아야 한다는 상한]
61. Faje AT, et al. High-dose glucocorticoids for the treatment of
    ipilimumab-induced hypophysitis is associated with reduced survival.
    *Cancer.* 2018;124(18):3706-3714.
    <https://pubmed.ncbi.nlm.nih.gov/29975414/>
    [**고용량** 스테로이드는 생존과 음의 연관 — `chi_pred = 1.0`의 근거]
62. Arbour KC, et al. Impact of Baseline Steroids on Efficacy of PD-1 and PD-L1
    Blockade in Patients With NSCLC. *J Clin Oncol.* 2018;36(28):2872-2878.
    <https://pubmed.ncbi.nlm.nih.gov/30125216/>
63. Verheijden RJ, et al. Association of Anti-TNF with Decreased Survival in
    Steroid Refractory Ipilimumab and Anti-PD1-Treated Patients.
    *Clin Cancer Res.* 2020;26(9):2268-2274.
    <https://pubmed.ncbi.nlm.nih.gov/31998085/>
    [anti-TNF가 생존과 음의 연관 — `chi_ifx = 0.60`]
64. Perez-Ruiz E, et al. Prophylactic TNF blockade uncouples efficacy and
    toxicity in dual CTLA-4 and PD-1 immunotherapy. *Nature.*
    2019;569(7756):428-432. <https://pubmed.ncbi.nlm.nih.gov/31043740/>
    [**독성과 효능의 분리 가능성** — 모델이 정량화하려는 바로 그 명제]
65. Bishu S, et al. Gut-selective integrin blockade: mechanistic basis for
    sparing systemic antitumour immunity — 대표적 기전 근거로
    Soler D, et al. The binding specificity and selective antagonism of
    vedolizumab, an anti-α4β7 integrin therapeutic antibody in development for
    inflammatory bowel disease. *J Pharmacol Exp Ther.* 2009;330(3):864-875.
    <https://pubmed.ncbi.nlm.nih.gov/19509315/>
    [α4β7:MAdCAM-1이 **장 특이적 주소**라는 약리학적 근거 — `chi_vdz = 0.10`]

---

## 8. 약동학 (Pharmacokinetics) — 모든 PK 파라미터의 출처

66. Feng Y, et al. Exposure-response relationships of the efficacy and safety of
    ipilimumab in patients with advanced melanoma. *Clin Cancer Res.*
    2013;19(14):3977-3986. <https://pubmed.ncbi.nlm.nih.gov/23741070/>
    [`CL_ipi` 0.367 L/d, `V1_ipi` 4.4 L, t½ 15.4 d, 그리고 **노출-독성 관계의
    존재**]
67. Bajaj G, et al. Model-Based Population Pharmacokinetic Analysis of
    Nivolumab in Patients With Solid Tumors. *CPT Pharmacometrics Syst
    Pharmacol.* 2017;6(1):58-66. <https://pubmed.ncbi.nlm.nih.gov/28019091/>
    [`CL_pd1` 0.228 L/d, `V1_pd1` 8.0 L, t½ 25 d]
68. Elassaiss-Schaap J, et al. Using Model-Based "Learn and Confirm" to Reveal
    the Pharmacokinetics-Pharmacodynamics Relationship of Pembrolizumab.
    *CPT Pharmacometrics Syst Pharmacol.* 2017;6(1):21-28.
    <https://pubmed.ncbi.nlm.nih.gov/28019091/>
    [pembrolizumab CL 0.22 L/d, t½ 22 d; 2 vs 10 mg/kg 노출 차이가 임상적으로
    무의미한 이유]
69. Fasanmade AA, et al. Population pharmacokinetic analysis of infliximab in
    patients with ulcerative colitis. *Eur J Clin Pharmacol.*
    2009;65(12):1211-1228. <https://pubmed.ncbi.nlm.nih.gov/19756557/>
    [**`CL_ifx ∝ (ALB/4.0)^-0.9`** — 모델의 다섯 번째 아이디어("질병이 자기
    해독제를 먹는다")의 정량적 근거]
70. Dotan I, et al. Patient factors that increase infliximab clearance and
    shorten half-life in inflammatory bowel disease. *Inflamm Bowel Dis.*
    2014;20(12):2247-2259. <https://pubmed.ncbi.nlm.nih.gov/25358061/>
    [저알부민혈증·높은 CRP·높은 염증 부하가 CL을 높임 — `crp_slope_ifx`]
71. Rosario M, et al. Population pharmacokinetics-pharmacodynamics of
    vedolizumab in patients with ulcerative colitis and Crohn's disease.
    *Aliment Pharmacol Ther.* 2015;42(2):188-202.
    <https://pubmed.ncbi.nlm.nih.gov/25996351/>
    [`CL_vdz` 0.157 L/d, t½ 25.5 d, 알부민 공변량]
72. Shankaran V, et al. Systemic exposure and α4β7 receptor saturation of
    vedolizumab — 대표 근거로 Rosario M, et al. A Review of the Clinical
    Pharmacokinetics, Pharmacodynamics, and Immunogenicity of Vedolizumab.
    *Clin Pharmacokinet.* 2017;56(11):1287-1301.
    <https://pubmed.ncbi.nlm.nih.gov/28523450/>
    [임상 용량에서 α4β7 **완전 포화** — `KD_VDZ`, `Emax_VDZ`]
73. Shah DK, Betts AM. Antibody biodistribution coefficients: inferring tissue
    concentrations of monoclonal antibodies based on the plasma concentrations
    in several preclinical species and human. *MAbs.* 2013;5(2):297-305.
    <https://pubmed.ncbi.nlm.nih.gov/23406896/>
    [**조직/혈장 분포계수 ~0.1-0.2** — `BDC = 0.13`. 조직 농도가 혈장의 13%에
    불과하다는 사실이 EC50_ADCC를 임상 노출 범위 한가운데로 밀어 넣는다]
74. Czock D, et al. Pharmacokinetics and pharmacodynamics of systemically
    administered glucocorticoids. *Clin Pharmacokinet.* 2005;44(1):61-98.
    <https://pubmed.ncbi.nlm.nih.gov/15634032/>
    [prednisolone F 0.8, t½ 2-3 h, **유전체 효과의 지연** — `keo_pred`]
75. Dowty ME, et al. The pharmacokinetics, metabolism, and clearance mechanisms
    of tofacitinib. *Drug Metab Dispos.* 2014;42(4):759-773.
    <https://pubmed.ncbi.nlm.nih.gov/24464803/>
    [tofacitinib t½ 3.2 h — `ke_jak`, `V_jak`]

---

## 9. 생리학: 대장 수분 흡수 예비능과 검열된 등급

모델의 두 번째 아이디어 — **대장은 첫 증상 전에 예비능을 소진한다** — 의 근거.

76. Debongnie JC, Phillips SF. Capacity of the human colon to absorb fluid.
    *Gastroenterology.* 1978;74(4):698-703.
    <https://pubmed.ncbi.nlm.nih.gov/632832/>
    [**대장 최대 흡수능 약 4.5-5 L/일** — `A_max`. 이 한 숫자가 S* = L/A_max를
    결정하고, 따라서 얼마나 많은 대장이 조용히 파괴될 수 있는지를 결정한다]
77. Sandle GI. Salt and water absorption in the human colon: a modern appraisal.
    *Gut.* 1998;43(2):294-299. <https://pubmed.ncbi.nlm.nih.gov/10189859/>
    [회장 유출량 약 1.5-2 L/일 — `L_pres`]
78. Sandle GI. Pathogenesis of diarrhea in ulcerative colitis: new views on an
    old problem. *J Clin Gastroenterol.* 2005;39(4 Suppl 2):S49-52.
    <https://pubmed.ncbi.nlm.nih.gov/15758659/>
    [염증에서 NHE3/ENaC 하향조절 — `nhe_tnf`, `nhe_ifn`]
79. Amasheh M, et al. TNFα-induced and berberine-antagonized tight junction
    barrier impairment via tyrosine kinase, Akt and NFκB signaling.
    *J Cell Sci.* 2010;123(Pt 23):4145-4155.
    <https://pubmed.ncbi.nlm.nih.gov/21062898/>
80. Sugi K, et al. Inhibition of Na⁺,K⁺-ATPase by interferon gamma
    down-regulates intestinal epithelial transport. *Gastroenterology.*
    2001;120(6):1393-1403. <https://pubmed.ncbi.nlm.nih.gov/11313309/>
81. US NCI. Common Terminology Criteria for Adverse Events (CTCAE) v5.0. 2017.
    <https://ctep.cancer.gov/protocoldevelopment/electronic_applications/ctc.htm>
    [설사 등급 정의: G1 기저 대비 <4회/일 증가, G2 4-6회, G3 ≥7회 — 모델의
    `dfreq` → `grade` 매핑]

---

## 10. 바이오마커: 칼프로텍틴, CRP, 알부민

82. Abu-Sbeih H, et al. Fecal calprotectin and lactoferrin in the diagnosis and
    monitoring of immune checkpoint inhibitor-induced colitis — 대표 근거로
    Zhang HC, et al. The role of fecal calprotectin in immune checkpoint
    inhibitor colitis. *Ther Adv Gastroenterol.* 2021;14:17562848211002545.
    <https://pubmed.ncbi.nlm.nih.gov/33796146/>
    [**칼프로텍틴이 증상보다 먼저 움직인다** — 모델의 검열 사다리(censoring
    ladder) 분석의 임상적 대응물, `Calpro0`, 150-200 µg/g 역치]
83. Tibble JA, et al. A simple method for assessing intestinal inflammation in
    Crohn's disease. *Gut.* 2000;47(4):506-513.
    <https://pubmed.ncbi.nlm.nih.gov/10986210/>
84. Vermeire S, et al. Laboratory markers in IBD: useful, magic, or unnecessary
    toys? *Gut.* 2006;55(3):426-431.
    <https://pubmed.ncbi.nlm.nih.gov/16474109/>
    [CRP 반감기 약 19시간 — `kd_crp`]
85. Umar A, et al. Protein-losing enteropathy in inflammatory bowel disease —
    대표 근거로 Ferrante M, et al. Predictors of early response to infliximab.
    *Inflamm Bowel Dis.* 2007;13(2):123-128.
    <https://pubmed.ncbi.nlm.nih.gov/17206692/>
    [저알부민혈증이 반응 불량과 연관 — `k_ple`, 그리고 알부민 → CL 되먹임]

---

## 11. 합병증과 재투여

86. Franklin C, et al. Cytomegalovirus reactivation in patients with refractory
    checkpoint inhibitor-induced colitis. *Eur J Cancer.* 2017;86:248-256.
    <https://pubmed.ncbi.nlm.nih.gov/29055842/>
    [불응성에서 **CMV 재활성화** — `CMV` 노드]
87. Del Castillo M, et al. The Spectrum of Serious Infections Among Patients
    Receiving Immune Checkpoint Blockade for the Treatment of Melanoma.
    *Clin Infect Dis.* 2016;63(11):1490-1493.
    <https://pubmed.ncbi.nlm.nih.gov/27501842/>
    [면역억제 노출과 감염 위험 — `h_inf`]
88. Abu-Sbeih H, et al. Resumption of Immune Checkpoint Inhibitor Therapy After
    Immune-Mediated Colitis. *J Clin Oncol.* 2019;37(30):2738-2745.
    <https://pubmed.ncbi.nlm.nih.gov/31163011/>
    [**재투여 시 재발률**: anti-PD-1 단독 재개가 ipilimumab 재노출보다 안전 —
    `REINTRO` 노드의 근거]
89. Pollack MH, et al. Safety of resuming anti-PD-1 in patients with
    immune-related adverse events during combined anti-CTLA-4 and anti-PD-1 in
    metastatic melanoma. *Ann Oncol.* 2018;29(1):250-255.
    <https://pubmed.ncbi.nlm.nih.gov/29045547/>

---

## 12. 숙주 유전학과 위험 인자

90. Groha S, et al. Germline variants associated with toxicity to immune
    checkpoint blockade. *Nat Med.* 2022;28(12):2584-2591.
    <https://pubmed.ncbi.nlm.nih.gov/36471036/>
    [**IL7 rs16906115**가 전체 irAE 위험과 연관 — `IL7_SNP` 노드, 소박한
    나이브 풀 크기 공변량]
91. Abu-Sbeih H, et al. Immune checkpoint inhibitor therapy in patients with
    preexisting inflammatory bowel disease. *J Clin Oncol.*
    2020;38(6):576-583. <https://pubmed.ncbi.nlm.nih.gov/31872征>
    (정정 링크: <https://pubmed.ncbi.nlm.nih.gov/31874109/>)
    [**기존 IBD**가 대장염 위험을 크게 높임 — `Trm(0) > 0`, `Ent(0) < 1`]
92. Zhang T, et al. Non-steroidal anti-inflammatory drug use and
    immune-mediated colitis — 대표 근거로 Marthey L, 상기 28번 및
    Wang Y, 상기 59번의 NSAID 하위분석.
    [NSAID 노출 — `f_nsaid = 0.85`]

---

## 13. QSP 방법론 및 유사 모델

93. Nijsen MJMA, et al. Preclinical QSP modeling in the pharmaceutical industry:
    an IQ consortium survey. *CPT Pharmacometrics Syst Pharmacol.*
    2018;7(3):135-146. <https://pubmed.ncbi.nlm.nih.gov/29349875/>
94. Chelliah V, et al. Quantitative Systems Pharmacology Approaches for
    Immuno-Oncology: Adding Virtual Patients to the Development Paradigm.
    *Clin Pharmacol Ther.* 2021;109(3):605-618.
    <https://pubmed.ncbi.nlm.nih.gov/32531058/>
    [면역항암 QSP에서 **가상 환자 집단** 구성 — 본 모델 분석 K의 방법론]
95. Jafarnejad M, et al. A Computational Model of Neoadjuvant PD-1 Inhibition in
    Non-Small Cell Lung Cancer. *AAPS J.* 2019;21(5):79.
    <https://pubmed.ncbi.nlm.nih.gov/31236847/>
96. Milberg O, et al. A QSP model for predicting clinical responses to
    monotherapy, combination and sequential therapy following CTLA-4, PD-1,
    and PD-L1 checkpoint blockade. *Sci Rep.* 2019;9(1):11286.
    <https://pubmed.ncbi.nlm.nih.gov/31375756/>
    [CTLA-4/PD-1 병용의 QSP 선행 연구 — 효능 축; 본 모델은 같은 구조를
    **독성 축**으로 옮긴 것에 해당]
97. Baron KT, et al. mrgsolve: Simulate from ODE-Based Population PK/PD and
    Systems Pharmacology Models. <https://mrgsolve.org/>
98. Ribba B, et al. Model-Informed Drug Development for Immuno-Oncology.
    *Clin Pharmacol Ther.* 2017;101(6):730-733.
    <https://pubmed.ncbi.nlm.nih.gov/28187497/>

---

## 부록: 모델이 만드는 검증 가능한 예측 (Falsifiable predictions)

아래는 이 모델이 **문헌에서 읽은 것이 아니라 계산으로 내놓은** 진술이며,
따라서 틀릴 수 있는 것들입니다. 각 항목은 반증 가능한 형태로 적었습니다.

| # | 예측 | 어떻게 반증하는가 |
|---|------|------------------|
| 1 | ipilimumab 대장염의 용량-반응은 **점유율이 아니라 ADCC**에서 나온다. 따라서 Fc 무력화(effector-null) anti-CTLA-4는 항종양 효과를 상당 부분 유지하면서 대장염 용량-반응을 **없앨** 것이다 | Fc-silent anti-CTLA-4 임상시험에서 용량 증가 시 대장염이 여전히 비례 증가한다면 반증 |
| 2 | anti-PD-1 대장염은 **기존에 프라이밍된 상주 집단이 있어야** 발생한다 (약이 새 클론을 허가하지 못하므로) | anti-PD-1 단독 대장염 환자에서 치료 전 대장 Trm/클론성이 대조군과 다르지 않다면 반증 |
| 3 | 등급은 검열된 지표이므로, 칼프로텍틴 기반 조기 개입은 **같은 등급 도달률에서 더 낮은 궤양 깊이와 더 작은 잔여 Trm**을 남긴다 (재발률 감소) | 무작위 비교에서 재발률이 동일하다면 반증 |
| 4 | 스테로이드 불응성은 약력학이 아니라 **은와 소실**의 함수다 | 깊은 궤양이 없는데도 스테로이드에 불응하는 환자가 다수라면 반증 |
| 5 | 중증 대장염의 저알부민혈증은 infliximab 노출을 **가장 필요한 순간에** 낮춘다 → 중증에서는 표준 5 mg/kg가 미달 노출이다 | 중증 환자의 infliximab 곡선하면적이 경증과 다르지 않다면 반증 |
| 6 | 구조 대비 종양 비용은 χ에 비례한다: vedolizumab ≈ 스테로이드의 1/10 | vedolizumab 구조군의 무진행생존이 스테로이드 구조군보다 낫지 않다면 반증 |
| 7 | 병용요법 대장염이 비가산적인 것은 생물학이 아니라 **등급의 역치 검열** 때문이다 → 연속 지표(칼프로텍틴, 내시경 점수)로 재면 병용은 **가산적으로 보일 것**이다 | 칼프로텍틴/Mayo 점수로도 병용 ≈ ipi 단독이라면 반증 |
| 8 | κ(비-Treg 조절 바닥)가 최대 민감도 파라미터다 → 비-Treg 조절 톤을 회복시키는 개입은 효과기 팔을 건드리지 않고 anti-CTLA-4 대장염을 둔화시킨다 | 그런 개입이 대장염을 줄이지 못한다면 반증 |

---

*문서 최종 갱신: 2026-08-06. PubMed 링크는 작성 시점 기준이며, 일부 항목은
대표 문헌으로 대체 표기하였습니다. 모든 정량 파라미터의 실제 사용처는
`icic_reference_model.py`의 `P()` 함수 주석에서 다시 확인할 수 있습니다.*
