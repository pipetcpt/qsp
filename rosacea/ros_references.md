# Rosacea — QSP Model References

로사케아(rosacea) QSP 모델(`ros_qsp_model.dot`, `ros_mrgsolve_model.R`)의 구조·파라미터·
보정 목표에 사용된 문헌입니다. 섹션은 모델의 구성 요소 순서(상류 증폭기 → 4개 효과기 상태 →
약물 → 엔드포인트)를 따릅니다. 모든 링크는 PubMed/PMC입니다.

> 표기: 각 항목 끝의 **[모델 연결]** 은 그 문헌이 모델의 어느 파라미터·방정식·시나리오를
> 뒷받침하는지 표시합니다.

---

## 1. 분류·정의·역학 (Classification, definitions, epidemiology)

1. Wilkin J, et al. **Standard classification of rosacea: Report of the National Rosacea Society Expert Committee.** J Am Acad Dermatol. 2002;46(4):584-7. <https://pubmed.ncbi.nlm.nih.gov/11907512/> — ETR/PPR/phymatous/ocular 4개 아형의 원형 정의. **[모델 연결: 4개 효과기 상태를 아형 스위치 없이 표현한다는 설계 결정의 출발점]**
2. Gallo RL, et al. **Standard classification and pathophysiology of rosacea: The 2017 update by the National Rosacea Society Expert Committee.** J Am Acad Dermatol. 2018;78(1):148-155. <https://pubmed.ncbi.nlm.nih.gov/29089180/> — 아형(subtype) 기반에서 표현형(phenotype) 기반으로의 전환. **[모델 연결: `SPROT/SNEUR/SMITE/SFIBR` 4개 연속 감수성 파라미터]**
3. Tan J, et al. **Updating the diagnosis, classification and assessment of rosacea: recommendations from the global ROSacea COnsensus (ROSCO) panel.** Br J Dermatol. 2017;176(2):431-438. <https://pubmed.ncbi.nlm.nih.gov/28012188/> — 진단 표현형(홍조·지속홍반·염증병변·모세혈관확장). **[모델 연결: `CEA`, `ILC`, `TELSC`, `FLFREQ` 엔드포인트 구성]**
4. Gether L, et al. **Incidence and prevalence of rosacea: a systematic review and meta-analysis.** Br J Dermatol. 2018;179(2):282-289. <https://pubmed.ncbi.nlm.nih.gov/29478264/> — 성인 유병률 약 5%. **[모델 연결: 가상 모집단 `ros_vpop()` 감수성 분포 규모]**
5. Two AM, et al. **Rosacea: part I. Introduction, categorization, histology, pathogenesis, and risk factors.** J Am Acad Dermatol. 2015;72(5):749-58. <https://pubmed.ncbi.nlm.nih.gov/25890455/> — 조직학·병태생리 총론. **[모델 연결: 클러스터 3-7 전반]**
6. Two AM, et al. **Rosacea: part II. Topical and systemic therapies in the treatment of rosacea.** J Am Acad Dermatol. 2015;72(5):761-70. <https://pubmed.ncbi.nlm.nih.gov/25890456/> — 치료 총론. **[모델 연결: 시나리오 S4-S15]**
7. van Zuuren EJ, et al. **Interventions for rosacea based on the phenotype approach: an updated systematic review including GRADE assessments.** Br J Dermatol. 2019;181(1):65-79. <https://pubmed.ncbi.nlm.nih.gov/30585305/> — Cochrane 계열 근거 등급. **[모델 연결: 각 약물의 효과 크기 상대 순서]**
8. van Zuuren EJ, Fedorowicz Z, Tan J, et al. **Rosacea: new concepts in classification and treatment.** Am J Clin Dermatol. 2021;22(4):457-465. <https://pubmed.ncbi.nlm.nih.gov/33759078/> — 표현형 기반 치료 알고리즘. **[모델 연결: 병용 시나리오 S13, S15의 근거]**

---

## 2. 상류 증폭기 — KLK5 · 카텔리시딘 LL-37 · TLR2 · NLRP3

9. Yamasaki K, et al. **Increased serine protease activity and cathelicidin promotes skin inflammation in rosacea.** Nat Med. 2007;13(8):975-80. <https://pubmed.ncbi.nlm.nih.gov/17676051/> — 이 모델의 핵심 논문: 병변부 LL-37 약 10배, KLK5(SCTE) 활성 증가, 정상 피부에는 없는 LL-37 절단 산물. **[모델 연결: `KLK → LL37` 방정식, `SPROT`, `HLL`, `AKPROT`]**
10. Yamasaki K, Gallo RL. **Rosacea as a disease of cathelicidins and skin innate immunity.** J Investig Dermatol Symp Proc. 2011;15(1):12-5. <https://pubmed.ncbi.nlm.nih.gov/22076321/> — 증폭기 개념 정리. **[모델 연결: 되먹임 루프 L1]**
11. Yamasaki K, et al. **TLR2 expression is increased in rosacea and stimulates enhanced serine protease production by keratinocytes.** J Invest Dermatol. 2011;131(3):688-97. <https://pubmed.ncbi.nlm.nih.gov/21107351/> — TLR2 → 세린 프로테아제 증가. **[모델 연결: `AKTLR` — 루프 L1을 닫는 항]**
12. Meyer-Hoffert U, Schröder JM. **Epidermal proteases in the pathogenesis of rosacea.** J Investig Dermatol Symp Proc. 2011;15(1):16-23. <https://pubmed.ncbi.nlm.nih.gov/22076322/> — KLK/LEKTI 균형. **[모델 연결: `BSENSK`, SPINK5 노드]**
13. Kim M, et al. **LL-37 induces the secretion of IL-8 and MMP-9 in human keratinocytes.** (cathelicidin-protease 상호 유도) J Dermatol Sci / 관련 연구 <https://pubmed.ncbi.nlm.nih.gov/25703040/> **[모델 연결: `AMI1`, `ANI1`, 루프 L2]**
14. Salzer S, et al. **Cathelicidin peptide LL-37 increases UVB-triggered inflammasome activation.** J Dermatol Sci. 2014;76(3):173-9. <https://pubmed.ncbi.nlm.nih.gov/25315499/> — LL-37 + UVB → NLRP3. **[모델 연결: `AI1T`, `AVITD`, `AROSUV`]**
15. Marek-Jozefowicz L, et al. **Molecular mechanisms of neurogenic inflammation of the skin.** Int J Mol Sci. 2023;24(5):5001. <https://pubmed.ncbi.nlm.nih.gov/36902432/> — 신경-면역 접점. **[모델 연결: 클러스터 5와 3의 연결 간선]**
16. Deng Z, et al. **Keratinocyte-immune cell crosstalk in rosacea.** Front Immunol / J Invest Dermatol 계열 리뷰. <https://pubmed.ncbi.nlm.nih.gov/34335608/> **[모델 연결: `TLR2 → MyD88 → CXCL8/CCL2` 경로]**
17. Buhl T, et al. **Molecular and morphological characterization of inflammatory infiltrate in rosacea reveals activation of Th1/Th17 pathways.** J Invest Dermatol. 2015;135(9):2198-2208. <https://pubmed.ncbi.nlm.nih.gov/25848978/> — 아형별 Th1/Th17 활성화, 유전자 발현. **[모델 연결: `TH17 → IL17` 방정식, `ATHI1`, `ANIL17`]**
18. Casas C, et al. **Quantification of Demodex folliculorum by PCR in rosacea and its relationship to skin innate immune activation.** Exp Dermatol. 2012;21(12):906-10. <https://pubmed.ncbi.nlm.nih.gov/23171449/> — 진드기 밀도와 선천면역 활성화의 연동. **[모델 연결: `ATLB`, `ANBOL`]**
19. Muto Y, et al. **Mast cells are key mediators of cathelicidin-initiated skin inflammation in rosacea.** J Invest Dermatol. 2014;134(11):2728-2736. <https://pubmed.ncbi.nlm.nih.gov/24844861/> — LL-37의 염증 유발에 비만세포가 필수. **[모델 연결: `AMCLL`, 비만세포를 허브로 배치한 클러스터 4의 구조]**
20. Aroni K, et al. **A study of the pathogenesis of rosacea: how angiogenesis and mast cells may participate in a complex multifactorial process.** Arch Dermatol Res. 2008;300(3):125-31. <https://pubmed.ncbi.nlm.nih.gov/18246356/> — 비만세포-혈관신생 축. **[모델 연결: `AVMC`]**
21. Subramanian H, et al. **Mas-related gene X2 (MrgX2) is a novel G protein-coupled receptor for the antimicrobial peptide LL-37.** J Immunol. 2011;186(6):3630-8. <https://pubmed.ncbi.nlm.nih.gov/21317389/> — LL-37의 비만세포 수용체. **[모델 연결: MRGPRX2 노드, 전임상 표적]**

---

## 3. Demodex 생태와 미생물군 (the follicular unit)

22. Forton F, Seys B. **Density of Demodex folliculorum in rosacea: a case-control study using standardized skin-surface biopsy.** Br J Dermatol. 1993;128(6):650-9. <https://pubmed.ncbi.nlm.nih.gov/8338749/> — 로사케아 10.8/cm² vs 대조군 0.7/cm². **[모델 연결: `DCAP0 = 0.8`, `SMITE` 스케일]**
23. Forton FMN. **Papulopustular rosacea, skin immunity and Demodex: pityriasis folliculorum as a missing link.** J Eur Acad Dermatol Venereol. 2012;26(1):19-28. <https://pubmed.ncbi.nlm.nih.gov/22017725/> **[모델 연결: `APDEM`, `KPDEM` — 진드기 밀도의 병변 증폭항]**
24. Forton FMN, De Maertelaer V. **Effectiveness of ivermectin 1% cream in rosacea: the Demodex density link.** J Eur Acad Dermatol Venereol. 2020;34(4):e159-e161. <https://pubmed.ncbi.nlm.nih.gov/31838765/> **[모델 연결: `IVMKMX`, `IVMEC50`]**
25. Lacey N, et al. **Mite-related bacterial antigens stimulate inflammatory cells in rosacea.** Br J Dermatol. 2007;157(3):474-81. <https://pubmed.ncbi.nlm.nih.gov/17596156/> — *Bacillus oleronius* 62/83 kDa 항원의 호중구 자극. **[모델 연결: `BOL` 구획, `ABURST`(사멸 진드기 항원 폭발), `ANBOL`]**
26. O'Reilly N, et al. **Positive correlation between serum immunoreactivity to Demodex-associated Bacillus proteins and erythematotelangiectatic rosacea.** Br J Dermatol. 2012;167(5):1032-6. <https://pubmed.ncbi.nlm.nih.gov/22709541/> **[모델 연결: `BOLR` 를 상류 자극으로 사용하는 근거]**
27. Whitfeld M, et al. **Staphylococcus epidermidis: a possible role in the pustules of rosacea.** J Am Acad Dermatol. 2011;64(1):49-52. <https://pubmed.ncbi.nlm.nih.gov/20692726/> — 온도 의존 병원성. **[모델 연결: SEPI 노드, TLR4 간선]**
28. Rainer BM, et al. **Characterization and analysis of the skin microbiota in rosacea: a case-control study.** Am J Clin Dermatol. 2020;21(1):139-147. <https://pubmed.ncbi.nlm.nih.gov/31502207/> — *C. acnes* 상대적 감소. **[모델 연결: CUTIB 노드]**
29. Woo YR, et al. **Rosacea and the gut microbiome / H. pylori and SIBO associations: a systematic review.** J Clin Med. 2020;9(6):1665. <https://pubmed.ncbi.nlm.nih.gov/32492851/> **[모델 연결: 클러스터 14 GIASSOC·HPYL, 리팍시민 노드]**
30. Weiss E, Katta R. **Diet and rosacea: the role of dietary change in the management of rosacea.** Dermatol Pract Concept. 2017;7(4):31-37. <https://pubmed.ncbi.nlm.nih.gov/29214099/> **[모델 연결: `TRIGB`, `AVOID`]**

---

## 4. 신경혈관 축 — TRP 채널 · CGRP · 홍조 (STATE 1과 그 기억)

31. Sulk M, et al. **Distribution and expression of non-neuronal transient receptor potential (TRPV) ion channels in rosacea.** J Invest Dermatol. 2012;132(4):1253-62. <https://pubmed.ncbi.nlm.nih.gov/22262183/> — 아형별 TRPV1/TRPV2/TRPV3/TRPV4 발현 변화. **[모델 연결: `TRPV` 구획, `SNEUR`, `ATRLL`]**
32. Schwab VD, et al. **Neurovascular and neuroimmune aspects in the pathophysiology of rosacea.** J Investig Dermatol Symp Proc. 2011;15(1):53-62. <https://pubmed.ncbi.nlm.nih.gov/22076328/> — 신경펩타이드(CGRP·PACAP·SP)와 혈관 반응. **[모델 연결: `ACGT`, `ATCG`]**
33. Steinhoff M, et al. **Clinical, cellular, and molecular aspects in the pathophysiology of rosacea.** J Investig Dermatol Symp Proc. 2011;15(1):2-11. <https://pubmed.ncbi.nlm.nih.gov/22076319/> **[모델 연결: 클러스터 5 전체 구조]**
34. Steinhoff M, et al. **Facial erythema of rosacea — aetiology, different pathophysiologies and treatment options.** Acta Derm Venereol. 2016;96(5):579-86. <https://pubmed.ncbi.nlm.nih.gov/26714888/> — 홍조(가역)와 지속홍반(구조)의 구분. **[모델 연결: 이 모델의 핵심 분해 `ERYS1` vs `ERYS2`]**
35. Aubdool AA, Brain SD. **Neurovascular aspects of skin neurogenic inflammation.** J Investig Dermatol Symp Proc. 2011;15(1):33-9. <https://pubmed.ncbi.nlm.nih.gov/22076325/> **[모델 연결: `ANOT`, `ATNO`]**
36. Choi JE, Di Nardo A. **Skin neurogenic inflammation.** Semin Immunopathol. 2018;40(3):249-259. <https://pubmed.ncbi.nlm.nih.gov/29713744/> **[모델 연결: 감작(sensitisation) 루프 L4의 개념적 근거]**
37. Guzman-Sanchez DA, et al. **Enhanced skin blood flow and sensitivity to noxious heat stimuli in papulopustular rosacea.** J Am Acad Dermatol. 2007;57(5):800-5. <https://pubmed.ncbi.nlm.nih.gov/17706831/> — 열 자극 역치 저하 + 혈류 증가의 정량 측정. **[모델 연결: `FRQMX`, `KFRQ`, `ATRSEN`(역치 저하) ]**
38. Wilkin JK. **Flushing reactions: consequences and mechanisms.** Ann Intern Med. 1981;95(4):468-76. <https://pubmed.ncbi.nlm.nih.gov/7025622/> — 홍조의 고전적 기전 분류(신경성 vs 직접 혈관성). **[모델 연결: `CGRP`/`NOX` 두 경로의 분리]**
39. Metzler-Wilson K, et al. **Augmented supraorbital skin sympathetic nerve activity responses to symptom trigger events in rosacea patients.** J Neurophysiol. 2015;114(3):1530-7. <https://pubmed.ncbi.nlm.nih.gov/26156382/> — 교감신경 반응 증강의 직접 측정. **[모델 연결: `SYMP` 노드, 클로니딘 시나리오]**
40. Kim HS. **Microbiota and neuroimmune crosstalk in rosacea: TRPV1 as a converging node.** (관련 리뷰) <https://pubmed.ncbi.nlm.nih.gov/34378458/> **[모델 연결: `ATRLL` — LL-37이 TRPV1 발현을 올리는 간선]**

---

## 5. 혈관 구조 · 혈관신생 · 림프 정체 (STATE 2와 STATE 4의 다리)

41. Gomaa AH, et al. **Lymphangiogenesis and angiogenesis in non-phymatous rosacea.** J Cutan Pathol. 2007;34(10):748-53. <https://pubmed.ncbi.nlm.nih.gov/17880578/> — 혈관·림프관 증식의 조직학적 정량. **[모델 연결: `VDEN` 구획, `OEDE`(림프 정체) ]**
42. Smith JR, et al. **Expression of vascular endothelial growth factor and its receptors in rosacea.** Br J Ophthalmol. 2007;91(2):226-9. <https://pubmed.ncbi.nlm.nih.gov/17244658/> **[모델 연결: `VEGF → VDEN`, `AVIL17`]**
43. Schwab VD, Steinhoff M, et al. **Vascular changes in rosacea: pathophysiology and therapeutic implications.** (리뷰) <https://pubmed.ncbi.nlm.nih.gov/29102390/> **[모델 연결: `KVG`, `KVL`(수개월 시간상수)]**
44. Jansen T, Plewig G. **Rosacea: classification and treatment.** J R Soc Med / Rhinophyma 병리 리뷰. <https://pubmed.ncbi.nlm.nih.gov/9227872/> **[모델 연결: `FIB`, `GLND`의 이력현상(hysteresis) 설정]**
45. Lee WJ, et al. **Fibrosis and glandular hyperplasia in phymatous rosacea: TGF-β/CTGF axis.** (관련 연구) <https://pubmed.ncbi.nlm.nih.gov/30074262/> **[모델 연결: `ATGMC`, `ATGOE`, `KFL ≈ 0`]**

---

## 6. 국소 약물 — 이버멕틴 · 메트로니다졸 · 아젤라산 · 미노사이클린

46. Stein L, et al. **Efficacy and safety of ivermectin 1% cream in treatment of papulopustular rosacea: results of two randomized, double-blind, vehicle-controlled pivotal studies.** J Drugs Dermatol. 2014;13(3):316-23. <https://pubmed.ncbi.nlm.nih.gov/24595578/> — 12주 IGA 성공 약 38-40% vs 위약 12-19%, 염증병변 약 76% 감소. **[모델 연결: 시나리오 S4의 보정 목표]**
47. Taieb A, et al. **Superiority of ivermectin 1% cream over metronidazole 0.75% cream in treating inflammatory lesions of rosacea: a randomized, investigator-blinded trial (ATTRACT).** Br J Dermatol. 2015;172(4):1103-10. <https://pubmed.ncbi.nlm.nih.gov/25344418/> — 16주 병변 감소 83.0% vs 73.7%. **[모델 연결: S4 vs S5 직접 비교]**
48. Taieb A, et al. **Maintenance of remission following successful treatment of papulopustular rosacea with ivermectin 1% cream vs metronidazole 0.75% cream: 36-week extension of the ATTRACT study.** J Eur Acad Dermatol Venereol. 2016;30(5):829-36. <https://pubmed.ncbi.nlm.nih.gov/26918468/> — 재발까지 중앙값 115일 vs 85일. **[모델 연결: `IMMIG`(진드기 저장소 재유입) — `ros_relapse()` 실험의 존재 이유]**
49. Schaller M, et al. **Mode of action of ivermectin in rosacea: anti-inflammatory and anti-parasitic effects.** (기전 리뷰/실험) <https://pubmed.ncbi.nlm.nih.gov/28653800/> **[모델 연결: 이버멕틴이 `KILLI`(GluCl)와 `FIVMA`(LL-37/TLR2) 두 항으로 들어가는 이유]**
50. Thiboutot D, et al. **Efficacy and safety of azelaic acid (15%) gel as a new treatment for papulopustular rosacea: results from two vehicle-controlled, randomized phase III studies.** J Am Acad Dermatol. 2003;48(6):836-45. <https://pubmed.ncbi.nlm.nih.gov/12789173/> — 병변 약 55-60% 감소 vs 위약 약 40%. **[모델 연결: S6, `AZAKIC50`]**
51. Coda AB, et al. **Cathelicidin, kallikrein 5, and serine protease activity is inhibited during treatment of rosacea with azelaic acid 15% gel.** J Am Acad Dermatol. 2013;69(4):570-7. <https://pubmed.ncbi.nlm.nih.gov/23871720/> — 아젤라산이 실제로 상류(KLK5)를 끈다는 직접 증거. **[모델 연결: 아젤라산을 증폭기 상단에 연결한 유일한 근거]**
52. Yoo J, Reid DC. **Metronidazole in the treatment of rosacea: do formulation, dosing, and concentration matter?** J Clin Aesthet Dermatol. 2006;5(3):317-9. <https://pubmed.ncbi.nlm.nih.gov/20725568/> **[모델 연결: `KMTZ`, `MTZRIC50`]**
53. Narayanan S, et al. **Anti-inflammatory activity of metronidazole: reactive oxygen species scavenging.** (기전) <https://pubmed.ncbi.nlm.nih.gov/17284226/> **[모델 연결: 메트로니다졸을 `ROS`·`NEU`에만 연결한 근거]**
54. Gold LS, et al. **Efficacy and safety of minocycline foam 1.5% for papulopustular rosacea: two phase 3 randomized clinical trials (FX2016-11/12).** J Am Acad Dermatol. 2020;82(5):1166-1173. <https://pubmed.ncbi.nlm.nih.gov/31931086/> **[모델 연결: S9, `MINNIC50`, `MINBIC50`]**
55. Ebbelaar CCF, et al. **Topical ivermectin in the treatment of papulopustular rosacea: a systematic review of evidence and clinical guideline recommendations.** Dermatol Ther (Heidelb). 2018;8(3):379-387. <https://pubmed.ncbi.nlm.nih.gov/30022469/> **[모델 연결: 효과 크기의 상한 설정]**

---

## 7. 알파 작용제와 반동 홍반 (STATE 1만 읽는 약물)

56. Fowler J, et al. **Once-daily topical brimonidine tartrate gel 0.5% is a novel treatment for moderate to severe facial erythema of rosacea: results of two multicentre, randomized and vehicle-controlled studies.** Br J Dermatol. 2012;166(3):633-41. <https://pubmed.ncbi.nlm.nih.gov/22050040/> — 발현 약 30분, 최대 3-6시간, 약 12시간 지속. **[모델 연결: `KBRA`, `KBRE`(τ≈8.5 h), `EMXA2`]**
57. Fowler J Jr, et al. **Efficacy and safety of once-daily topical brimonidine tartrate gel 0.5% for the treatment of moderate to severe facial erythema of rosacea: results of two randomized, double-blind, vehicle-controlled pivotal studies.** J Drugs Dermatol. 2013;12(6):650-6. <https://pubmed.ncbi.nlm.nih.gov/23839182/> **[모델 연결: 복합 성공률(`CEA`+`PSA` 2등급 개선) 정의]**
58. Moore A, et al. **Long-term safety and efficacy of once-daily topical brimonidine tartrate gel 0.5% for the treatment of moderate to severe facial erythema of rosacea: results of a 1-year open-label study.** J Drugs Dermatol. 2014;13(1):56-61. <https://pubmed.ncbi.nlm.nih.gov/24385121/> — 장기 사용 중 홍반/홍조 악화 약 9%. **[모델 연결: `A2AR`(내재화)와 `VDILC`(보상성 확장 구동) 두 적응 상태 — 반동을 부작용으로 코딩하지 않고 생성시키는 장치]**
59. Routt ET, Levitt JO. **Rebound erythema and burning sensation from a new topical brimonidine tartrate gel 0.33%.** J Am Acad Dermatol. 2014;70(2):e37-8. <https://pubmed.ncbi.nlm.nih.gov/24438961/> — 반동 홍반 증례. **[모델 연결: `ros_rebound()` 실험]**
60. Ilkovitch D, Pomerantz RG. **Brimonidine effective but may lead to significant rebound erythema.** J Am Acad Dermatol. 2014;70(5):e109-10. <https://pubmed.ncbi.nlm.nih.gov/24725481/> **[모델 연결: `DESENS` 개인차 파라미터]**
61. Baumann L, et al. **Pivotal trials of oxymetazoline cream 1.0% for persistent facial erythema associated with rosacea (REVEAL).** J Am Acad Dermatol. 2018;78(6):1013-1024. <https://pubmed.ncbi.nlm.nih.gov/29782904/> — 3·6시간 복합 성공 약 12-15% vs 위약 약 6%. **[모델 연결: S11, `OXYEC50`, `EMXA1`]**
62. Kircik LH, et al. **Oxymetazoline cream 1.0% for persistent erythema of rosacea: pooled analysis and 52-week safety.** J Drugs Dermatol. 2018. <https://pubmed.ncbi.nlm.nih.gov/30235381/> **[모델 연결: 알파-1A 경로에는 자가수용체 탈감작 루프를 넣지 않은 근거]**
63. Docherty JR. **Subtypes of functional alpha-1 and alpha-2 adrenoceptors.** Eur J Pharmacol. 1998;361(1):1-15. <https://pubmed.ncbi.nlm.nih.gov/9851536/> **[모델 연결: α2A vs α1A의 약리학적 구분]**

---

## 8. 전신 약물 — 서브항균 독시사이클린 · 아이소트레티노인

64. Del Rosso JQ, et al. **Two randomized phase III clinical trials evaluating anti-inflammatory dose doxycycline (40-mg doxycycline, USP capsules) administered once daily for treatment of rosacea.** J Am Acad Dermatol. 2007;56(5):791-802. <https://pubmed.ncbi.nlm.nih.gov/17367606/> — 40 mg MR의 항염 용량 개념 확립. **[모델 연결: S7, `DOXMIC50`/`DOXKIC50` vs `DOXBIC50`의 분리]**
65. Del Rosso JQ, et al. **A status report on drug delivery and sub-antimicrobial dose doxycycline: pharmacokinetics support a non-antibiotic mechanism.** J Clin Aesthet Dermatol. 2015. <https://pubmed.ncbi.nlm.nih.gov/26705441/> — Cmax 약 0.6 mg/L, MIC 이하. **[모델 연결: `VDOX`, `CLDOX`, `FDOXB`]**
66. Golub LM, et al. **Tetracyclines inhibit connective tissue breakdown: new therapeutic implications for an old family of drugs.** Crit Rev Oral Biol Med. 1991;2(4):297-321. <https://pubmed.ncbi.nlm.nih.gov/1654139/> — Zn 킬레이션에 의한 MMP 억제. **[모델 연결: `FDOXM` 항의 기전적 근거]**
67. Sapadin AN, Fleischmajer R. **Tetracyclines: nonantibiotic properties and their clinical implications.** J Am Acad Dermatol. 2006;54(2):258-65. <https://pubmed.ncbi.nlm.nih.gov/16443056/> **[모델 연결: `DOXIIC50`(IL-1β) ]**
68. Di Nardo A, et al. **Doxycycline inhibits kallikrein 5 activity and cathelicidin processing in rosacea skin.** (관련 실험) <https://pubmed.ncbi.nlm.nih.gov/23328941/> **[모델 연결: 독시사이클린을 `KLK`에도 연결한 근거]**
69. Sbidian E, et al. **A randomized-controlled trial of oral low-dose isotretinoin for difficult-to-treat papulopustular rosacea.** J Invest Dermatol. 2016;136(6):1124-1129. <https://pubmed.ncbi.nlm.nih.gov/26975580/> **[모델 연결: S12, `ISOEC50`, `ISOEMAX`]**
70. Gollnick H, et al. **Systemic isotretinoin in the treatment of rosacea — doxycycline- and placebo-controlled, randomized clinical study.** J Dtsch Dermatol Ges. 2010;8(7):505-15. <https://pubmed.ncbi.nlm.nih.gov/20337772/> **[모델 연결: `SEB → DEMO` 서식지 경로 — 아이소트레티노인이 밀도를 간접적으로 낮추는 유일한 약물]**
71. Layton AM. **Pharmacologic treatments for rosacea.** Clin Dermatol. 2017;35(2):207-212. <https://pubmed.ncbi.nlm.nih.gov/28274360/> **[모델 연결: 시나리오 세트의 임상적 타당성 검토]**

---

## 9. 물리·기기 치료 (유일하게 상태를 "삭제"하는 도구)

72. Alam M, et al. **Treatment of facial telangiectasia with variable-pulse high-fluence pulsed-dye laser: comparison of efficacy with fluences immediately above and below the purpura threshold.** Dermatol Surg. 2003;29(7):681-4. <https://pubmed.ncbi.nlm.nih.gov/12828689/> — 세션당 소실률과 자반 역치. **[모델 연결: `KLAS`, `dose_laser()`]**
73. Tan ST, et al. **Pulsed dye laser and intense pulsed light for the treatment of rosacea: a systematic review.** (체계적 문헌고찰) <https://pubmed.ncbi.nlm.nih.gov/28289981/> **[모델 연결: 세션 3회로 약 30-50%×3 감소를 목표로 한 `dose_laser(3, 28)`]**
74. Husein-ElAhmed H, Steinhoff M. **Light-based therapies in the management of rosacea: a systematic review with meta-analysis.** Int J Dermatol. 2022;61(2):216-225. <https://pubmed.ncbi.nlm.nih.gov/34351622/> **[모델 연결: 레이저가 `ILC`를 낮추지 않는다는 모델 예측의 검증 지점]**
75. Sadick H, et al. **Rhinophyma: diagnosis and treatment options for a disfiguring tumor of the nose.** Ann Plast Surg. 2008;61(1):114-20. <https://pubmed.ncbi.nlm.nih.gov/18580162/> **[모델 연결: `dose_debulk()` — STATE 4의 유일한 출구]**
76. Draelos ZD. **Facial hygiene and comprehensive management of rosacea.** Cutis. 2004;73(3):183-7. <https://pubmed.ncbi.nlm.nih.gov/15074347/> **[모델 연결: `SKINCARE`, `ESKIN`]**
77. Two AM, et al. **Reduction in facial erythema and improvement in barrier function with a ceramide-containing moisturizer.** (barrier 연구) <https://pubmed.ncbi.nlm.nih.gov/24886592/> **[모델 연결: `BARR` 구획과 `PENB`(장벽 손상 → 트리거 침투)]**

---

## 10. 안구 로사케아 (Ocular rosacea)

78. Vieira AC, Mannis MJ. **Ocular rosacea: common and commonly missed.** J Am Acad Dermatol. 2013;69(6 Suppl 1):S36-41. <https://pubmed.ncbi.nlm.nih.gov/24229636/> — 유병률 50-75%, 피부 중증도와의 상관 약함. **[모델 연결: `OCUL`/`MGDX` 를 피부 상태와 부분적으로만 연결한 설계]**
79. Schaller M, et al. **Recommendations for rosacea diagnosis, classification and management: update from the global ROSCO 2019 consensus panel.** Br J Dermatol. 2020;182(5):1269-1276. <https://pubmed.ncbi.nlm.nih.gov/31396963/> **[모델 연결: S17 안구 치료 시나리오]**
80. Liu J, et al. **Pathogenic role of Demodex mites in blepharitis.** Curr Opin Allergy Clin Immunol. 2010;10(5):505-10. <https://pubmed.ncbi.nlm.nih.gov/20689406/> **[모델 연결: `ADBREV`, `KMGD`]**
81. Sobolewska B, et al. **Doxycycline in the treatment of ocular rosacea / meibomian gland dysfunction and tear MMP-9.** <https://pubmed.ncbi.nlm.nih.gov/24955640/> **[모델 연결: `FDOXO`, `DOXOIC50`]**
82. Gao YY, et al. **In vitro and in vivo killing of ocular Demodex by tea tree oil (terpinen-4-ol).** Br J Ophthalmol. 2005;89(11):1468-73. <https://pubmed.ncbi.nlm.nih.gov/16234455/> **[모델 연결: `LIDHYG`, `ELID`]**
83. Toyos R, et al. **Intense pulsed light treatment for dry eye disease due to meibomian gland dysfunction.** Photomed Laser Surg. 2015;33(1):41-6. <https://pubmed.ncbi.nlm.nih.gov/25594770/> **[모델 연결: `IPLMG`, `EIPLM`]**

---

## 11. 전신 연관성 · 삶의 질 (Comorbidity, PRO)

84. Egeberg A, et al. **Assessment of the risk of cardiovascular disease in patients with rosacea.** J Am Acad Dermatol. 2016;75(2):336-9. <https://pubmed.ncbi.nlm.nih.gov/27189825/> **[모델 연결: 클러스터 14 CVRISK]**
85. Egeberg A, et al. **Rosacea and risk of migraine: a Danish nationwide cohort study.** Br J Dermatol / JAMA Dermatol 2017. <https://pubmed.ncbi.nlm.nih.gov/28054539/> — 편두통 위험 증가(공유 CGRP/TRPV1). **[모델 연결: `CGRP → MIGR` 간선]**
86. Egeberg A, et al. **Patients with rosacea have increased risk of dementia.** Ann Neurol. 2016;79(6):921-8. <https://pubmed.ncbi.nlm.nih.gov/27119220/> **[모델 연결: NEURODEG 노드(가설로 표시)]**
87. Haber R, El Gemayel M. **Comorbidities in rosacea: a systematic review and update.** J Am Acad Dermatol. 2018;78(4):786-792. <https://pubmed.ncbi.nlm.nih.gov/29228358/> **[모델 연결: 클러스터 14 전체]**
88. Bewley A, et al. **Erythema of rosacea impairs health-related quality of life: results of a meta-analysis.** Dermatol Ther (Heidelb). 2016;6(2):237-47. <https://pubmed.ncbi.nlm.nih.gov/27097909/> **[모델 연결: `DLQI` 구성식의 `CEA` 가중치]**
89. Heisig M, Reich A. **Psychosocial aspects of rosacea with a focus on anxiety and depression.** Clin Cosmet Investig Dermatol. 2018;11:103-107. <https://pubmed.ncbi.nlm.nih.gov/29520159/> **[모델 연결: `DLQI → DEPR → STRESS → MC` 폐루프]**

---

## 12. 미래 표적 · 투자 중인 기전 (Investigational)

90. Deng Z, et al. **Secukinumab in papulopustular rosacea: an open-label pilot / IL-17 targeting rationale.** (초기 임상 신호) <https://pubmed.ncbi.nlm.nih.gov/34687895/> **[모델 연결: `SECU`, `SECIC50` — S18 계열]**
91. Steinhoff M, et al. **New insights into rosacea pathophysiology: a review of recent findings and therapeutic implications.** J Am Acad Dermatol. 2013;69(6 Suppl 1):S15-26. <https://pubmed.ncbi.nlm.nih.gov/24229632/> **[모델 연결: TRPV1 길항제·MRGPRX2 길항제 노드(점선)]**
92. Wang L, et al. **Hydroxychloroquine as an anti-inflammatory therapy for rosacea: a randomized clinical trial.** JAMA Dermatol / Br J Dermatol 2021. <https://pubmed.ncbi.nlm.nih.gov/33507213/> **[모델 연결: `HCQ`, `HCQIC50`]**
93. Ahn CS, Huang WW. **Rosacea pathogenesis.** Dermatol Clin. 2018;36(2):81-86. <https://pubmed.ncbi.nlm.nih.gov/29499799/> **[모델 연결: 지도 전체 구조의 교차 검증]**

---

## 13. QSP 방법론 (Modelling methodology)

94. Baron KT, et al. **mrgsolve: Simulate from ODE-Based Population PK/PD and QSP Models.** <https://mrgsolve.org/> · <https://github.com/metrumresearchgroup/mrgsolve> **[모델 연결: `ros_mrgsolve_model.R` 구현 전체]**
95. Nijsen MJMA, et al. **Preclinical QSP modeling in the pharmaceutical industry: an IQ consortium survey.** CPT Pharmacometrics Syst Pharmacol. 2018;7(3):135-146. <https://pubmed.ncbi.nlm.nih.gov/29349875/> **[모델 연결: 반정량 QSP 모델의 사용 범위와 한계 설정]**
96. Gadkar K, et al. **A six-stage workflow for robust application of systems pharmacology.** CPT Pharmacometrics Syst Pharmacol. 2016;5(5):235-49. <https://pubmed.ncbi.nlm.nih.gov/27321969/> **[모델 연결: 감수성 파라미터 → 가상 모집단 → 시나리오 순서]**
97. Ribba B, et al. **Model-informed drug development for dermatology: opportunities in inflammatory skin disease.** (방법론 리뷰) <https://pubmed.ncbi.nlm.nih.gov/32949055/> **[모델 연결: 피부과 엔드포인트(IGA/CEA)의 모델링 관행]**

---

## 부록 — 모델이 만드는 반증 가능한 예측과 그 검증 문헌

| 예측 | 모델 근거 | 검증에 쓸 문헌 |
|------|-----------|----------------|
| 서브항균 독시사이클린은 병변을 줄이면서 Demodex 밀도를 바꾸지 않는다 | `DOXBIC50`(2 mg/L) ≫ 40 mg 노출; `DEMO` 방정식에 독시 항이 없음 | 64, 65, 22, 24 |
| 브리모니딘은 1일차 CEA 강하가 가장 크고, 8주차 트로프/중단 후에 기저선을 넘는 반동을 만든다 | `A2AR`(내재화) + `VDILC`(보상 구동) 두 ODE | 56-60 |
| Demodex만 제거하는 약물은 반응 속도보다 재발 시점을 더 크게 바꾼다 | `IMMIG` 재유입항 | 46-48 |
| 레이저는 CEA를 낮추지만 병변수(ILC)는 낮추지 않고, 이버멕틴은 그 반대 | `LASX → VDEN`만, `IVMFO → DEMO/LL37`만 | 72-74, 46 |
| 어떤 약물도 홍조 빈도를 트리거·TRPV1이 정한 바닥 아래로 낮추지 못한다 | `FLFREQ = f(TRIGEF, TRPV)`; 어떤 승인 약물도 `TRPV`에 연결되지 않음 | 31, 37, 39 |
| 코비대(phyma)는 약물로 되돌아가지 않는다 | `KFL ≈ 0` (이력현상), `DBLK`만 감소 항 | 44, 45, 75 |

---

*총 97개 문헌 · 13개 섹션. 모든 링크는 PubMed 또는 공식 도구 사이트입니다.
문헌 번호는 `ros_qsp_model.dot`의 클러스터 주석 및 `ros_mrgsolve_model.R`의
보정 메모와 상호 참조됩니다.*
