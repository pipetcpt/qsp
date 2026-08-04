# 니만-피크병 C형 (NPC) — 참고문헌
# Niemann-Pick Disease Type C — Annotated References

이 모델(`npc_qsp_model.dot` · `npc_mrgsolve_model.R` · `npc_reference_model.py` ·
`npc_shiny_app.R`)의 모든 구조적 가정과 파라미터의 근거 문헌입니다.

**모든 PMID는 NCBI E-utilities로 개별 조회하여 제목·저자·연도·저널을 확인했습니다.**
링크 형식: `https://pubmed.ncbi.nlm.nih.gov/<PMID>/`

모델의 **정량적 보정 목표(calibration target)** 로 직접 사용된 문헌은 🎯 로,
보정에 쓰이지 않고 **모델 예측의 검증(validation)** 에만 쓰인 문헌은 🔍 로 표시했습니다.

---

## 목차

| 절 | 주제 | 편수 |
|----|------|------|
| [1](#1-유전자--분자-기반-genetics-and-molecular-basis) | 유전자 · 분자 기반 | 10 |
| [2](#2-리소좀-콜레스테롤-배출-기전-lysosomal-cholesterol-egress) | 리소좀 콜레스테롤 배출 기전 | 13 |
| [3](#3-2차-지질-축적--리소좀-칼슘-secondary-storage-and-lysosomal-calcium) | 2차 지질 축적 · 리소좀 칼슘 | 8 |
| [4](#4-오토파지--mtorc1--tfeb-autophagy-and-nutrient-signalling) | 오토파지 · mTORC1 · TFEB | 8 |
| [5](#5-신경병리--신경염증-neuropathology-and-neuroinflammation) | 신경병리 · 신경염증 | 13 |
| [6](#6-바이오마커-biomarkers) | 바이오마커 | 14 |
| [7](#7-자연사--중증도-척도-natural-history-and-severity-scales) | 자연사 · 중증도 척도 | 19 |
| [8](#8-미글루스타트-miglustat) | 미글루스타트 | 13 |
| [9](#9-아리모클로몰-arimoclomol) | 아리모클로몰 | 7 |
| [10](#10-레바세틸류신-levacetylleucine--n-acetyl-l-leucine) | 레바세틸류신 | 11 |
| [11](#11-사이클로덱스트린-2-hpβcd--adrabetadex) | 사이클로덱스트린 | 12 |
| [12](#12-개발-중--실험적-접근-emerging-approaches) | 개발 중 접근 | 8 |
| [13](#13-qsp-방법론-qsp-methodology) | QSP 방법론 | 8 |
| | **합계** | **144** |

---

## 1. 유전자 · 분자 기반 (Genetics and molecular basis)

모델 클러스터 2(`cluster_gene`)와 상태변수 `NPC1_ER`·`NPC1_L`, 유전형 표
`GENOTYPES`의 근거.

1. Carstea ED et al. **Niemann-Pick C1 disease gene: homology to mediators of cholesterol homeostasis.** *Science* 1997. — NPC1 최초 클로닝, 13-TM 막단백질 및 sterol-sensing domain 상동성. [PMID 9211849](https://pubmed.ncbi.nlm.nih.gov/9211849/)
2. Naureckiene S et al. **Identification of HE1 as the second gene of Niemann-Pick C disease.** *Science* 2000. — NPC2(=HE1) 동정. 모델에서 `f_npc2`로 반영. [PMID 11125141](https://pubmed.ncbi.nlm.nih.gov/11125141/)
3. Vanier MT. **Niemann-Pick disease type C.** *Orphanet J Rare Dis* 2010. — 표준 종설. 발병 연령형 분류(주산기·조기영아·후기영아·청소년·성인)의 근거. [PMID 20525256](https://pubmed.ncbi.nlm.nih.gov/20525256/)
4. Wassif CA et al. **High incidence of unrecognized visceral/neurological late-onset Niemann-Pick disease, type C1, predicted by analysis of massively parallel sequencing data sets.** *Genet Med* 2016. — 성인형 과소진단, 보인자 빈도. `mild/mild` 유전형의 근거. [PMID 25764212](https://pubmed.ncbi.nlm.nih.gov/25764212/)
5. Nakasone N et al. **Endoplasmic reticulum-associated degradation of Niemann-Pick C1: evidence for the role of heat shock proteins and identification of lysine residues that accept ubiquitin.** *J Biol Chem* 2014. — 🎯 변이 NPC1의 ERAD와 **HSP 의존적 접힘 구제**. 모델의 `kerad`·`Emax_fold`·`theta_eff` 구조의 직접 근거. [PMID 24891511](https://pubmed.ncbi.nlm.nih.gov/24891511/)
6. Burton BK et al. **Estimating the prevalence of Niemann-Pick disease type C (NPC) in the United States.** *Mol Genet Metab* 2021. [PMID 34304992](https://pubmed.ncbi.nlm.nih.gov/34304992/)
7. Labrecque M et al. **Estimated prevalence of Niemann-Pick type C disease in Quebec.** *Sci Rep* 2021. [PMID 34799641](https://pubmed.ncbi.nlm.nih.gov/34799641/)
8. Geberhiwot T et al. **Consensus clinical management guidelines for Niemann-Pick disease type C.** *Orphanet J Rare Dis* 2018. — 진단·치료 표준. [PMID 29625568](https://pubmed.ncbi.nlm.nih.gov/29625568/)
9. Yoon HJ et al. **The point mutation of the cholesterol trafficking membrane protein NPC1 may affect its proper function.** *Comput Biol Chem* 2022. [PMID 35850050](https://pubmed.ncbi.nlm.nih.gov/35850050/)
10. Elghobashi-Meinhardt N. **Cholesterol Transport in Wild-Type NPC1 and P691S: Molecular Dynamics Simulations Reveal Changes in Dynamical Behavior.** *Int J Mol Sci* 2020. — 변이별 잔여 기능 차이. [PMID 32331453](https://pubmed.ncbi.nlm.nih.gov/32331453/)

---

## 2. 리소좀 콜레스테롤 배출 기전 (Lysosomal cholesterol egress)

모델 클러스터 3(`cluster_egress`)과 상태변수 `CHOL_V`·`CHOL_C`,
Michaelis-Menten 배출항 `Vmax * f_eg * CHOL/(Km + CHOL)`의 근거.

11. Pentchev PG et al. **A defect in cholesterol esterification in Niemann-Pick disease (type C) patients.** *Proc Natl Acad Sci U S A* 1985. — 원발 결함의 최초 생화학적 기술. [PMID 3865225](https://pubmed.ncbi.nlm.nih.gov/3865225/)
12. Liscum L, Faust JR. **Low density lipoprotein (LDL)-mediated suppression of cholesterol synthesis and LDL uptake is defective in Niemann-Pick type C fibroblasts.** *J Biol Chem* 1987. — 🎯 **저장은 넘치는데 ER은 굶는다**는 SREBP2 역설(클러스터 4)의 직접 근거. [PMID 3680287](https://pubmed.ncbi.nlm.nih.gov/3680287/)
13. Slotte JP et al. **Intracellular transport of cholesterol in type C Niemann-Pick fibroblasts.** *Biochim Biophys Acta* 1989. [PMID 2804059](https://pubmed.ncbi.nlm.nih.gov/2804059/)
14. Infante RE et al. **NPC2 facilitates bidirectional transfer of cholesterol between NPC1 and lipid bilayers, a step in cholesterol egress from lysosomes.** *Proc Natl Acad Sci U S A* 2008. — NPC2→NPC1 직렬 전달. 모델에서 `f_eg = f_NPC1 × f_npc2`(곱셈, 직렬)로 반영. [PMID 18772377](https://pubmed.ncbi.nlm.nih.gov/18772377/)
15. Wang ML et al. **Identification of surface residues on Niemann-Pick C2 essential for hydrophobic handoff of cholesterol to NPC1 in lysosomes.** *Cell Metab* 2010. — 소수성 손넘김. [PMID 20674861](https://pubmed.ncbi.nlm.nih.gov/20674861/)
16. Li X et al. **Clues to the mechanism of cholesterol transfer from the structure of NPC1 middle lumenal domain bound to NPC2.** *Proc Natl Acad Sci U S A* 2016. [PMID 27551080](https://pubmed.ncbi.nlm.nih.gov/27551080/)
17. Qian H et al. **Structural Basis of Low-pH-Dependent Lysosomal Cholesterol Egress by NPC1 and NPC2.** *Cell* 2020. — 🎯 배출이 **낮은 pH에 의존**한다는 구조적 근거. 모델에서 리소좀 pH 상승(`Kph`)이 가수분해효소 활성뿐 아니라 배출 자체를 손상시키는 되먹임의 근거. [PMID 32544384](https://pubmed.ncbi.nlm.nih.gov/32544384/)
18. Pfeffer SR. **NPC intracellular cholesterol transporter 1 (NPC1)-mediated cholesterol export from lysosomes.** *J Biol Chem* 2019. [PMID 30710017](https://pubmed.ncbi.nlm.nih.gov/30710017/)
19. Sandhu J et al. **Aster Proteins Facilitate Nonvesicular Plasma Membrane to ER Cholesterol Transport in Mammalian Cells.** *Cell* 2018. — 클러스터 3의 `EGRESS_PM` 경로. [PMID 30220461](https://pubmed.ncbi.nlm.nih.gov/30220461/)
20. Naito T, Saheki Y. **GRAMD1-mediated accessible cholesterol sensing and transport.** *Biochim Biophys Acta Mol Cell Biol Lipids* 2021. [PMID 33932585](https://pubmed.ncbi.nlm.nih.gov/33932585/)
21. Ferrari A et al. **Aster Proteins Regulate the Accessible Cholesterol Pool in the Plasma Membrane.** *Mol Cell Biol* 2020. [PMID 32719109](https://pubmed.ncbi.nlm.nih.gov/32719109/)
22. Long T et al. **Structural basis for itraconazole-mediated NPC1 inhibition.** *Nat Commun* 2020. — 약물에 의한 NPC1 억제(역방향 검증). [PMID 31919352](https://pubmed.ncbi.nlm.nih.gov/31919352/)
23. Elghobashi-Meinhardt N. **Niemann-Pick type C disease: a QM/MM study of conformational changes in cholesterol in the NPC1(NTD) and NPC2 binding pockets.** *Biochemistry* 2014. [PMID 25251378](https://pubmed.ncbi.nlm.nih.gov/25251378/)

---

## 3. 2차 지질 축적 · 리소좀 칼슘 (Secondary storage and lysosomal calcium)

모델 클러스터 5·6과 상태변수 `SPH`·`CA_LY`·`GSL_V`·`GSL_C`의 근거.

24. Lloyd-Evans E et al. **Niemann-Pick disease type C1 is a sphingosine storage disease that causes deregulation of lysosomal calcium.** *Nat Med* 2008. — 🎯 **스핑고신이 최초 축적 지질이며 산성 Ca²⁺ 저장 재충전을 차단**한다는 핵심 논문. 모델의 `SPH → CA_LY` 억제항(`Ksph`)과 `SPH`가 `f_NPC1` 의존적 배출을 갖는 구조의 근거. [PMID 18953351](https://pubmed.ncbi.nlm.nih.gov/18953351/)
25. Lloyd-Evans E, Platt FM. **Lysosomal Ca²⁺ homeostasis: role in pathogenesis of lysosomal storage diseases.** *Cell Calcium* 2011. [PMID 21724254](https://pubmed.ncbi.nlm.nih.gov/21724254/)
26. Zervas M et al. **Critical role for glycosphingolipids in Niemann-Pick disease type C.** *Curr Biol* 2001. — GSL 합성 감소가 마우스 표현형을 개선 → 미글루스타트의 근거. [PMID 11525744](https://pubmed.ncbi.nlm.nih.gov/11525744/)
27. Vanier MT. **Lipid changes in Niemann-Pick disease type C brain: personal experience and review of the literature.** *Neurochem Res* 1999. — 🔍 NPC 뇌의 GM2/GM3 상승 배수. 모델의 CNS GSL 축적(모델 예측 ~2.6배)의 검증 근거. [PMID 10227680](https://pubmed.ncbi.nlm.nih.gov/10227680/)
28. Shen D et al. **Lipid storage disorders block lysosomal trafficking by inhibiting a TRP channel and lysosomal calcium release.** *Nat Commun* 2012. — 콜레스테롤에 의한 TRPML1 억제(`Ktr_chol`). [PMID 22415822](https://pubmed.ncbi.nlm.nih.gov/22415822/)
29. Pagano RE. **Endocytic trafficking of glycosphingolipids in sphingolipid storage diseases.** *Philos Trans R Soc Lond B Biol Sci* 2003. [PMID 12803922](https://pubmed.ncbi.nlm.nih.gov/12803922/)
30. Walkley SU et al. **Initiation and growth of ectopic neurites and meganeurites during postnatal cortical development in ganglioside storage disease.** *Brain Res Dev Brain Res* 1990. — GM2에 의한 메가뉴라이트(클러스터 11 `MEGANEURITE`). [PMID 2108821](https://pubmed.ncbi.nlm.nih.gov/2108821/)
31. Kuech EM et al. **Alterations in membrane trafficking and pathophysiological implications in lysosomal storage disorders.** *Biochimie* 2016. [PMID 27664461](https://pubmed.ncbi.nlm.nih.gov/27664461/)

---

## 4. 오토파지 · mTORC1 · TFEB (Autophagy and nutrient signalling)

모델 클러스터 7과 상태변수 `AUTOPH`·`HYD`의 근거.

32. Elrick MJ et al. **Impaired proteolysis underlies autophagic dysfunction in Niemann-Pick type C disease.** *Hum Mol Genet* 2012. — 🎯 오토파지 **유도는 정상이고 분해가 실패**한다는 점. 모델이 `AUTOPH`를 유입-상수·분해-감소 구조로 쓴 근거. [PMID 22872701](https://pubmed.ncbi.nlm.nih.gov/22872701/)
33. Settembre C et al. **TFEB links autophagy to lysosomal biogenesis.** *Science* 2011. [PMID 21617040](https://pubmed.ncbi.nlm.nih.gov/21617040/)
34. Castellano BM et al. **Lysosomal cholesterol activates mTORC1 via an SLC38A9-Niemann-Pick C1 signaling complex.** *Science* 2017. — 리소좀 콜레스테롤 → mTORC1. [PMID 28336668](https://pubmed.ncbi.nlm.nih.gov/28336668/)
35. Lim CY et al. **ER-lysosome contacts enable cholesterol sensing by mTORC1 and drive aberrant growth signalling in Niemann-Pick type C.** *Nat Cell Biol* 2019. [PMID 31548609](https://pubmed.ncbi.nlm.nih.gov/31548609/)
36. Davis OB et al. **NPC1-mTORC1 Signaling Couples Cholesterol Sensing to Organelle Homeostasis and Is a Targetable Pathway in Niemann-Pick Type C.** *Dev Cell* 2021. [PMID 33308480](https://pubmed.ncbi.nlm.nih.gov/33308480/)
37. Kataura T et al. **Targeting the autophagy-NAD axis protects against cell death in Niemann-Pick type C1 disease models.** *Cell Death Dis* 2024. [PMID 38821960](https://pubmed.ncbi.nlm.nih.gov/38821960/)
38. Dai S et al. **Methyl-β-cyclodextrin restores impaired autophagy flux in Niemann-Pick C1-deficient cells through activation of AMPK.** *Autophagy* 2017. [PMID 28613987](https://pubmed.ncbi.nlm.nih.gov/28613987/)
39. Lee H et al. **Pathological roles of the VEGF/SphK pathway in Niemann-Pick type C neurons.** *Nat Commun* 2014. [PMID 25417698](https://pubmed.ncbi.nlm.nih.gov/25417698/)

---

## 5. 신경병리 · 신경염증 (Neuropathology and neuroinflammation)

모델 클러스터 11·12·13과 상태변수 `PC`·`PC_S`·`PC_LOST`·`INFL`·`SYN`·`CBL`의 근거.

40. Higashi Y et al. **Cerebellar degeneration in the Niemann-Pick type C mouse.** *Acta Neuropathol* 1993. — 푸르킨예 세포 소실의 최초 기술. [PMID 8382896](https://pubmed.ncbi.nlm.nih.gov/8382896/)
41. Sarna JR et al. **Patterned Purkinje cell degeneration in mouse models of Niemann-Pick type C disease.** *J Comp Neurol* 2003. — 전엽→후엽 구배. [PMID 12528192](https://pubmed.ncbi.nlm.nih.gov/12528192/)
42. German DC et al. **Neurodegeneration in the Niemann-Pick C mouse: glial involvement.** *Neuroscience* 2002. [PMID 11823057](https://pubmed.ncbi.nlm.nih.gov/11823057/)
43. Elrick MJ et al. **Conditional Niemann-Pick C mice demonstrate cell autonomous Purkinje cell neurodegeneration.** *Hum Mol Genet* 2010. — 🎯 **뉴런 자율적 사멸**. 모델이 `CHOL_CNS → PC_S`를 세포자율 경로로, 신경염증을 *증폭기*로만 쓴 근거. [PMID 20007718](https://pubmed.ncbi.nlm.nih.gov/20007718/)
44. Dinkel L et al. **Myeloid cell-specific loss of NPC1 in mice recapitulates microgliosis and neurodegeneration in patients.** *Sci Transl Med* 2024. — 미세아교세포 단독 결손도 신경퇴행을 낳는다 → `INFL → PC_S` 증폭 고리의 근거. [PMID 39630885](https://pubmed.ncbi.nlm.nih.gov/39630885/)
45. Love S et al. **Neurofibrillary tangles in Niemann-Pick disease type C.** *Brain* 1995. — AD형 NFT. [PMID 7894998](https://pubmed.ncbi.nlm.nih.gov/7894998/)
46. Bu B et al. **Niemann-Pick disease type C yields possible clue for why cerebellar neurons do not form neurofibrillary tangles.** *Neurobiol Dis* 2002. [PMID 12505421](https://pubmed.ncbi.nlm.nih.gov/12505421/)
47. Malnar M et al. **Bidirectional links between Alzheimer's disease and Niemann-Pick type C disease.** *Neurobiol Dis* 2014. [PMID 24907492](https://pubmed.ncbi.nlm.nih.gov/24907492/)
48. Mattsson N et al. **Gamma-secretase-dependent amyloid-beta is increased in Niemann-Pick type C: a cross-sectional study.** *Neurology* 2011. [PMID 21205675](https://pubmed.ncbi.nlm.nih.gov/21205675/)
49. Woś M et al. **Mitochondrial dysfunction in fibroblasts derived from patients with Niemann-Pick type C disease.** *Arch Biochem Biophys* 2016. — 모델의 `MITO`·`ROS`. [PMID 26869201](https://pubmed.ncbi.nlm.nih.gov/26869201/)
50. Takikita S et al. **Perturbed myelination process of premyelinating oligodendrocyte in Niemann-Pick type C mouse.** *J Neuropathol Exp Neurol* 2004. — 🎯 발달 취약성(`v_dev`·`tau_dev`)의 근거: 같은 생화학적 손상이 활발한 수초화 시기에 더 큰 비용을 치른다. [PMID 15217094](https://pubmed.ncbi.nlm.nih.gov/15217094/)
51. Burbulla LF et al. **Modeling Brain Pathology of Niemann-Pick Disease Type C Using Patient-Derived Neurons.** *Mov Disord* 2021. [PMID 33438272](https://pubmed.ncbi.nlm.nih.gov/33438272/)
52. Wheeler S, Sillence DJ. **Niemann-Pick type C disease: cellular pathology and pharmacotherapy.** *J Neurochem* 2020. — 세포병리·약물 종설. [PMID 31608980](https://pubmed.ncbi.nlm.nih.gov/31608980/)

---

## 6. 바이오마커 (Biomarkers)

모델 클러스터 9와 상태변수 `TRIOL`·`PPCS`·`TCG`·`NFL`·`CALB`의 근거.
**모델의 1번 구조적 주장(구획 주장)이 여기에 걸려 있습니다.**

53. Porter FD et al. **Cholesterol oxidation products are sensitive and specific blood-based biomarkers for Niemann-Pick C1 disease.** *Sci Transl Med* 2010. — 7-KC·C-triol의 최초 확립. [PMID 21048217](https://pubmed.ncbi.nlm.nih.gov/21048217/)
54. Jiang X et al. **A sensitive and specific LC-MS/MS method for rapid diagnosis of Niemann-Pick C1 disease from human plasma.** *J Lipid Res* 2011. [PMID 21518695](https://pubmed.ncbi.nlm.nih.gov/21518695/)
55. Kuchar L et al. **Quantitation of plasmatic lysosphingomyelin and lysosphingomyelin-509 for differential screening of Niemann-Pick A/B and C diseases.** *Anal Biochem* 2017. [PMID 28259515](https://pubmed.ncbi.nlm.nih.gov/28259515/)
56. Sidhu R et al. **N-acyl-O-phosphocholineserines: structures of a novel class of lipids that are biomarkers for Niemann-Pick C1 disease.** *J Lipid Res* 2019. — lysoSM-509의 실제 구조는 PPCS. 모델 상태변수 `PPCS`. [PMID 31201291](https://pubmed.ncbi.nlm.nih.gov/31201291/)
57. Jiang X et al. **Development of a bile acid-based newborn screen for Niemann-Pick disease type C.** *Sci Transl Med* 2016. — 담즙산 B(TCG). 모델 상태변수 `TCG`. [PMID 27147587](https://pubmed.ncbi.nlm.nih.gov/27147587/)
58. Jiang X et al. **Diagnosis of Niemann-Pick C1 by measurement of bile acid biomarkers in archived newborn dried blood spots.** *Mol Genet Metab* 2019. — 신생아기부터 이미 상승 → 모델이 표지자를 출생 직후 상승하는 것으로 쓴 근거. [PMID 30172462](https://pubmed.ncbi.nlm.nih.gov/30172462/)
59. Mazzacuva F et al. **Identification of novel bile acids as biomarkers for the early diagnosis of Niemann-Pick C disease.** *FEBS Lett* 2016. [PMID 27139891](https://pubmed.ncbi.nlm.nih.gov/27139891/)
60. Bradbury A et al. **Cerebrospinal Fluid Calbindin D Concentration as a Biomarker of Cerebellar Disease Progression in Niemann-Pick Type C1 Disease.** *J Pharmacol Exp Ther* 2016. — 🎯 CSF 칼빈딘이 **푸르킨예 사멸 flux**의 판독값. 모델에서 `CALB`를 사멸 속도(`die`)에 비례시킨 근거 — 저장량이 아니라 *사멸 속도*라는 점이 결정적. [PMID 27307499](https://pubmed.ncbi.nlm.nih.gov/27307499/)
61. Alam MS et al. **Plasma signature of neurological disease in the monogenetic disorder Niemann-Pick Type C.** *J Biol Chem* 2014. [PMID 24488491](https://pubmed.ncbi.nlm.nih.gov/24488491/)
62. Sidhu R et al. **A validated LC-MS/MS assay for quantification of 24(S)-hydroxycholesterol in plasma and cerebrospinal fluid.** *J Lipid Res* 2015. — 뇌 유래 산화 스테롤(모델의 `OHC24` 노드). [PMID 25866316](https://pubmed.ncbi.nlm.nih.gov/25866316/)
63. Agrawal N et al. **Neurofilament light chain in cerebrospinal fluid as a novel biomarker in evaluating both clinical severity and therapeutic response in Niemann-Pick disease type C1.** *Genet Med* 2023. — 모델 상태변수 `NFL`. [PMID 36470574](https://pubmed.ncbi.nlm.nih.gov/36470574/)
64. Eratne D et al. **Plasma neurofilament light chain is increased in Niemann-Pick Type C but glial fibrillary acidic protein is not.** *Acta Neuropsychiatr* 2024. [PMID 38533577](https://pubmed.ncbi.nlm.nih.gov/38533577/)
65. Stern S et al. **Evaluation of the landscape of pharmacodynamic biomarkers in Niemann-Pick Disease Type C (NPC).** *Orphanet J Rare Dis* 2024. — 🔍 표지자–임상 연결의 현재 한계를 정리. 모델의 구획 주장과 직접 대응. [PMID 39061081](https://pubmed.ncbi.nlm.nih.gov/39061081/)
66. Pataj Z et al. **Quantification of oxysterols in human plasma and red blood cells by liquid chromatography high-resolution tandem mass spectrometry.** *J Chromatogr A* 2016. [PMID 26607314](https://pubmed.ncbi.nlm.nih.gov/26607314/)

---

## 7. 자연사 · 중증도 척도 (Natural history and severity scales)

모델의 임상 척도 매핑(`SARA`·`NPCCSS5`·`NPCCSS4`·`NPCCSS17`)과
**2번 구조적 주장(예비능 주장)** 의 근거.

67. Yanjanin NM et al. **Linear clinical progression, independent of age of onset, in Niemann-Pick disease, type C.** *Am J Med Genet B Neuropsychiatr Genet* 2010. — 🎯🔍 **모델 구조를 결정한 논문.** 진행이 (a)선형이고 (b)발병 연령과 무관하다는 두 사실 때문에 손상을 *적분*으로, 사멸 속도를 *포화하는 관문*으로 썼습니다. 선형성은 보정에, 발병연령-무관성은 검증에 사용. [PMID 19415691](https://pubmed.ncbi.nlm.nih.gov/19415691/)
68. Patterson MC et al. **Validation of the 5-domain Niemann-Pick type C Clinical Severity Scale.** *Orphanet J Rare Dis* 2021. — 5영역 척도. [PMID 33579322](https://pubmed.ncbi.nlm.nih.gov/33579322/)
69. Mengel E et al. **Clinical disease progression and biomarkers in Niemann-Pick disease type C: a prospective cohort study.** *Orphanet J Rare Dis* 2020. — 🎯 **모델의 4개 정량 목표의 출처**: 5영역 연 1.5점, 17영역 연 ~2.9점, 혈장 triol 환자 88.31 vs 대조 5.97 ng/mL, triol–5영역 Spearman ρ = 0.265. [PMID 33228797](https://pubmed.ncbi.nlm.nih.gov/33228797/)
70. Mengel E et al. **Correction to: Clinical disease progression and biomarkers in Niemann-Pick disease type C.** *Orphanet J Rare Dis* 2021. [PMID 34074315](https://pubmed.ncbi.nlm.nih.gov/34074315/)
71. Cortina-Borja M et al. **Annual severity increment score as a tool for stratifying patients with Niemann-Pick disease type C and for recruitment to clinical trials.** *Orphanet J Rare Dis* 2018. [PMID 30115089](https://pubmed.ncbi.nlm.nih.gov/30115089/)
72. Imrie J et al. **The natural history of Niemann-Pick disease type C in the UK.** *J Inherit Metab Dis* 2007. [PMID 17160617](https://pubmed.ncbi.nlm.nih.gov/17160617/)
73. Mengel E et al. **Niemann-Pick disease type C symptomatology: an expert-based clinical description.** *Orphanet J Rare Dis* 2013. [PMID 24135395](https://pubmed.ncbi.nlm.nih.gov/24135395/)
74. Stampfer M et al. **Niemann-Pick disease type C clinical database: cognitive and coordination deficits are early disease indicators.** *Orphanet J Rare Dis* 2013. [PMID 23433426](https://pubmed.ncbi.nlm.nih.gov/23433426/)
75. Walterfang M et al. **Dysphagia as a risk factor for mortality in Niemann-Pick disease type C: systematic literature review and evidence from studies with miglustat.** *Orphanet J Rare Dis* 2012. — 🎯 모델의 생존 위험함수를 **연하기능의 제곱**에 비례시킨 근거(`h_swal`). [PMID 23039766](https://pubmed.ncbi.nlm.nih.gov/23039766/)
76. Bianconi SE et al. **Evaluation of age of death in Niemann-Pick disease, type C: Utility of disease support group websites to understand natural history.** *Mol Genet Metab* 2019. — 🔍 사망 연령 분포. [PMID 30850267](https://pubmed.ncbi.nlm.nih.gov/30850267/)
77. Gardin A et al. **A Retrospective Multicentric Study of 34 Patients with Niemann-Pick Type C Disease and Early Liver Involvement.** *J Pediatr* 2023. — 모델의 `h_liver`·주산기형. [PMID 36265573](https://pubmed.ncbi.nlm.nih.gov/36265573/)
78. Schmitz-Hübsch T et al. **Scale for the assessment and rating of ataxia: development of a new clinical scale.** *Neurology* 2006. — SARA 원논문(0-40점). [PMID 16769946](https://pubmed.ncbi.nlm.nih.gov/16769946/)
79. Bremova-Ertl T et al. **A cross-sectional, prospective ocular motor study in 72 patients with Niemann-Pick disease type C.** *Eur J Neurol* 2021. — 사카드 속도 정량. [PMID 34096670](https://pubmed.ncbi.nlm.nih.gov/34096670/)
80. Hopf S et al. **Vertical saccadic palsy and foveal retinal thinning in Niemann-Pick disease type C.** *PLoS One* 2021. [PMID 34086834](https://pubmed.ncbi.nlm.nih.gov/34086834/)
81. Grillini A et al. **Measuring saccades in patients with Niemann-Pick type C: A comparison between video-oculography and a novel method.** *Clin Park Relat Disord* 2022. [PMID 36338825](https://pubmed.ncbi.nlm.nih.gov/36338825/)
82. King KA et al. **Auditory phenotype of Niemann-Pick disease, type C1.** *Ear Hear* 2014. — 🎯 질환 자체의 청력 저하 기저치. 사이클로덱스트린 이독성과 **구분**하기 위해 필요. [PMID 24225652](https://pubmed.ncbi.nlm.nih.gov/24225652/)
83. King KA et al. **Hearing loss is an early consequence of Npc1 gene deletion in the mouse model of Niemann-Pick disease, type C.** *J Assoc Res Otolaryngol* 2014. [PMID 24839095](https://pubmed.ncbi.nlm.nih.gov/24839095/)
84. Ong LT et al. **Psychosis symptoms associated with Niemann-Pick disease type C.** *Psychiatr Genet* 2021. — 성인형의 정신증상 초발. [PMID 34133410](https://pubmed.ncbi.nlm.nih.gov/34133410/)
85. Patterson MC et al. **Recommendations for the diagnosis and management of Niemann-Pick disease type C: an update.** *Mol Genet Metab* 2012. [PMID 22572546](https://pubmed.ncbi.nlm.nih.gov/22572546/)

---

## 8. 미글루스타트 (Miglustat)

모델 클러스터 16과 `mig_*` 파라미터의 근거.

86. Patterson MC et al. **Miglustat for treatment of Niemann-Pick C disease: a randomised controlled study.** *Lancet Neurol* 2007. — 🎯 1차 지표는 수평 사카드 속도(HSEM). 연하기능 개선·청력 안정·보행 악화 지연. 모델이 미글루스타트 효과를 **내장 > CNS 비대칭**으로 예측하는 것과 대응. [PMID 17689147](https://pubmed.ncbi.nlm.nih.gov/17689147/)
87. Wraith JE et al. **Miglustat in adult and juvenile patients with Niemann-Pick disease type C: long-term data from a clinical trial.** *Mol Genet Metab* 2010. [PMID 20045366](https://pubmed.ncbi.nlm.nih.gov/20045366/)
88. Patterson MC et al. **Long-term miglustat therapy in children with Niemann-Pick disease type C.** *J Child Neurol* 2010. [PMID 19822772](https://pubmed.ncbi.nlm.nih.gov/19822772/)
89. Patterson MC et al. **Stable or improved neurological manifestations during miglustat therapy in patients from the international disease registry for Niemann-Pick disease type C.** *Orphanet J Rare Dis* 2015. [PMID 26017010](https://pubmed.ncbi.nlm.nih.gov/26017010/)
90. Patterson MC et al. **Treatment outcomes following continuous miglustat therapy in patients with Niemann-Pick disease Type C.** *Orphanet J Rare Dis* 2020. [PMID 32334605](https://pubmed.ncbi.nlm.nih.gov/32334605/)
91. Pineda M et al. **Miglustat in Niemann-Pick disease type C patients: a review.** *Orphanet J Rare Dis* 2018. [PMID 30111334](https://pubmed.ncbi.nlm.nih.gov/30111334/)
92. Solomon BI et al. **Association of Miglustat With Swallowing Outcomes in Niemann-Pick Disease, Type C1.** *JAMA Neurol* 2020. — 🔍 연하 결과에 대한 가장 정량적인 근거. 모델의 `SWALLOW` 경로 검증. [PMID 32897301](https://pubmed.ncbi.nlm.nih.gov/32897301/)
93. Platt FM et al. **Prevention of lysosomal storage in Tay-Sachs mice treated with N-butyldeoxynojirimycin.** *Science* 1997. — 기질 감소 요법의 원리. [PMID 9103204](https://pubmed.ncbi.nlm.nih.gov/9103204/)
94. Butters TD et al. **Imino sugar inhibitors for treating the lysosomal glycosphingolipidoses.** *Glycobiology* 2005. — 🎯 GCS 억제 IC₅₀ 범위(`mig_IC50`). [PMID 15901676](https://pubmed.ncbi.nlm.nih.gov/15901676/)
95. Andersson U et al. **N-butyldeoxygalactonojirimycin: a more selective inhibitor of glycosphingolipid biosynthesis than N-butyldeoxynojirimycin, in vitro and in vivo.** *Biochem Pharmacol* 2000. — 장 이당분해효소 억제(설사)의 기전적 근거. [PMID 10718340](https://pubmed.ncbi.nlm.nih.gov/10718340/)
96. Shayman JA. **The development and use of small molecule inhibitors of glycosphingolipid metabolism for lysosomal storage diseases.** *J Lipid Res* 2014. [PMID 24534703](https://pubmed.ncbi.nlm.nih.gov/24534703/)
97. Maegawa GH et al. **Pharmacokinetics, safety and tolerability of miglustat in the treatment of pediatric patients with GM2 gangliosidosis.** *Mol Genet Metab* 2009. — 🎯 소아 PK. 모델의 `mig_V`·`mig_CL`·체표면적 기준 용량. [PMID 19447653](https://pubmed.ncbi.nlm.nih.gov/19447653/)
98. Belmatoug N et al. **Gastrointestinal disturbances and their management in miglustat-treated patients.** *J Inherit Metab Dis* 2011. — 모델의 `MIG_GI_AE` → 용량 감량 경로. [PMID 21779792](https://pubmed.ncbi.nlm.nih.gov/21779792/)

---

## 9. 아리모클로몰 (Arimoclomol)

모델 클러스터 17과 `ari_*`·`Emax_fold` 파라미터의 근거.

99. Mengel E et al. **Efficacy and safety of arimoclomol in Niemann-Pick disease type C: Results from a double-blind, randomised, placebo-controlled, multinational phase 2/3 trial of a novel treatment.** *J Inherit Metab Dis* 2021. — 🎯 12개월 5영역 NPCCSS 차이 −1.40 (95% CI −2.76, −0.03; p = 0.046); 위약군 진행 2.15; **미글루스타트 병용 하위군 −2.06 (p = 0.006)**. 모델의 `Emax_fold` 보정 목표이자, 곱셈적 상호작용 예측의 검증 대상. [PMID 34418116](https://pubmed.ncbi.nlm.nih.gov/34418116/)
100. Mengel E et al. **Efficacy results from a 12-month double-blind randomized trial of arimoclomol for treatment of Niemann-Pick disease type C (NPC): Presenting a rescored 4-domain NPC Clinical Severity Scale.** *Mol Genet Metab Rep* 2025. — 🔍 R4DNPCCSS 차이 −1.70 (95% CI −3.05, −0.34; p = 0.016). 인지영역 제외·연하영역 재채점 이유. 모델의 `n4_frac`·`n4_gain` 예측 검증. [PMID 40520915](https://pubmed.ncbi.nlm.nih.gov/40520915/)
101. Kirkegaard T et al. **Hsp70 stabilizes lysosomes and reverts Niemann-Pick disease-associated lysosomal pathology.** *Nature* 2010. — 🎯 HSP70에 의한 리소좀 안정화(`e_hyd_hsp`). [PMID 20111001](https://pubmed.ncbi.nlm.nih.gov/20111001/)
102. Petersen NH, Kirkegaard T. **Connecting Hsp70, sphingolipid metabolism and lysosomal stability.** *Cell Cycle* 2010. [PMID 20519957](https://pubmed.ncbi.nlm.nih.gov/20519957/)
103. Gray J et al. **Heat shock protein amplification improves cerebellar myelination in the Npc1^nih mouse model.** *EBioMedicine* 2022. [PMID 36455410](https://pubmed.ncbi.nlm.nih.gov/36455410/)
104. Keam SJ. **Arimoclomol: First Approval.** *Drugs* 2025. — 🎯 승인 용량(체중별 47/62/93/124 mg tid)과 PK. [PMID 39715913](https://pubmed.ncbi.nlm.nih.gov/39715913/)
105. Cudkowicz ME et al. **Arimoclomol at dosages up to 300 mg/day is well tolerated and safe in amyotrophic lateral sclerosis.** *Muscle Nerve* 2008. — 🎯 선형 PK, t½ ~4 h, **용량 의존적 CSF 농도 상승**. 모델의 `ari_Kp_csf`. [PMID 18551622](https://pubmed.ncbi.nlm.nih.gov/18551622/)

---

## 10. 레바세틸류신 (Levacetylleucine / N-acetyl-L-leucine)

모델 클러스터 18과 `nal_*` 파라미터의 근거.

106. Bremova-Ertl T et al. **Trial of N-Acetyl-l-Leucine in Niemann-Pick Disease Type C.** *N Engl J Med* 2024. — 🎯 IB1001-301. 12주 SARA 변화 −1.97 ± 2.43(약) vs −0.60 ± 2.39(위약), LS 평균차 **−1.28 (95% CI −1.91, −0.65; p < 0.001)**; 기저 SARA 15.91 ± 7.65. 모델의 `nal_Emax_sym` 보정 목표. [PMID 38294974](https://pubmed.ncbi.nlm.nih.gov/38294974/)
107. Bremova-Ertl T et al. **Efficacy and safety of N-acetyl-L-leucine in Niemann-Pick disease type C.** *J Neurol* 2022. [PMID 34387740](https://pubmed.ncbi.nlm.nih.gov/34387740/)
108. Patterson MC et al. **Disease-Modifying, Neuroprotective Effect of N-Acetyl-l-Leucine in Adult and Pediatric Patients With Niemann-Pick Disease Type C.** *Neurology* 2025. — 🔍 장기 연장. 기저치 대비 SARA 변화 12개월 −1.88 ± 2.89, 18개월 −1.64 ± 3.24. **모델이 재현하지 못하는 항목입니다**: 이 값에 `nal_Emax_dm`을 맞추려면 CNS 지질 유입을 247% 차단해야 하므로(물리적으로 불가능) 보정에 쓰지 않고, 마우스 데이터로 지지되는 8%를 쓰고 불일치를 그대로 보고했습니다. README 7절 ② 참조. [PMID 40513057](https://pubmed.ncbi.nlm.nih.gov/40513057/)
109. Fields T et al. **N-acetyl-L-leucine for Niemann-Pick type C: a multinational double-blind randomized placebo-controlled crossover study.** *Trials* 2023. — 🎯 **교차설계·12주 기간·세척기**. 모델의 3번 구조적 주장(설계 주장)의 직접 근거. [PMID 37248494](https://pubmed.ncbi.nlm.nih.gov/37248494/)
110. Fields T et al. **A master protocol to investigate a novel therapy acetyl-L-leucine for three ultra-rare neurodegenerative diseases.** *Trials* 2021. [PMID 33482890](https://pubmed.ncbi.nlm.nih.gov/33482890/)
111. Bremova T et al. **Acetyl-dl-leucine in Niemann-Pick type C: A case series.** *Neurology* 2015. — 🎯 **수일~수주 내 효과 발현, 중단 시 재악화**. 모델이 효과부위 `ke0`를 0.15/d(t½ ~4.6일)로 잡고 효과를 `D_rev`에 걸은 근거. [PMID 26400580](https://pubmed.ncbi.nlm.nih.gov/26400580/)
112. Kaya E et al. **Acetyl-leucine slows disease progression in lysosomal storage disorders.** *Brain Commun* 2021. — 🎯 마우스 진행 지연 → `nal_Emax_dm`(질병조절 성분)의 존재 근거. [PMID 33738443](https://pubmed.ncbi.nlm.nih.gov/33738443/)
113. Günther L et al. **N-acetyl-L-leucine accelerates vestibular compensation after unilateral labyrinthectomy by action in the cerebellum and thalamus.** *PLoS One* 2015. — 대증 기전(`NAL_VEST`). [PMID 25803613](https://pubmed.ncbi.nlm.nih.gov/25803613/)
114. Vibert N, Vidal PP. **In vitro effects of acetyl-DL-leucine (tanganil) on central vestibular neurons and vestibulo-ocular networks of the guinea-pig.** *Eur J Neurosci* 2001. [PMID 11207808](https://pubmed.ncbi.nlm.nih.gov/11207808/)
115. Schniepp R et al. **Acetyl-DL-leucine improves gait variability in patients with cerebellar ataxia - a case series.** *Cerebellum Ataxias* 2016. [PMID 27073690](https://pubmed.ncbi.nlm.nih.gov/27073690/)
116. Martakis K et al. **Safety and efficacy of levacetylleucine in ataxia-telangiectasia: a phase 3, randomised, double-blind trial.** *Lancet Neurol* 2026. — 🔍 다른 질환에서 같은 설계·같은 크기의 효과. 대증 기전 해석의 외부 검증. [PMID 42309084](https://pubmed.ncbi.nlm.nih.gov/42309084/)

---

## 11. 사이클로덱스트린 (2-HPβCD / adrabetadex)

모델 클러스터 19와 `cd_*` 파라미터의 근거.
**기전-분리불가 독성(efficacy와 ototoxicity가 같은 기전)의 근거가 여기 있습니다.**

117. Ory DS et al. **Intrathecal 2-hydroxypropyl-β-cyclodextrin decreases neurological disease progression in Niemann-Pick disease, type C1: a non-randomised, open-label, phase 1-2 trial.** *Lancet* 2017. — 🎯 공개라벨. 진행 21/21(역사적 대조) vs 7/14(치료). **"NSS minus hearing"으로 채점** — 약이 청력을 망가뜨리기 때문에 청력을 지표에서 빼야 했다는 사실이 모델의 독성 항의 직접 근거. [PMID 28803710](https://pubmed.ncbi.nlm.nih.gov/28803710/)
118. Liu B et al. **Reversal of defective lysosomal transport in NPC disease ameliorates liver dysfunction and neurodegeneration in the npc1-/- mouse.** *Proc Natl Acad Sci U S A* 2009. — 단회 투여로도 저장 감소. [PMID 19171898](https://pubmed.ncbi.nlm.nih.gov/19171898/)
119. Vite CH et al. **Intracisternal cyclodextrin prevents cerebellar dysfunction and Purkinje cell death in feline Niemann-Pick type C1 disease.** *Sci Transl Med* 2015. — 🎯 고양이 모델에서 푸르킨예 사멸 예방. 모델에서 사이클로덱스트린이 `CHOL_C`를 NPC1 무관하게 직접 감소시키는 항의 근거. [PMID 25717099](https://pubmed.ncbi.nlm.nih.gov/25717099/)
120. Vance JE, Peake KB. **Function of the Niemann-Pick type C proteins and their bypass by cyclodextrin.** *Curr Opin Lipidol* 2011. — 🎯 **NPC1 우회**. 모델의 유전형-무관 예측(null/null·NPC2에서도 효과 유지)의 근거. [PMID 21412152](https://pubmed.ncbi.nlm.nih.gov/21412152/)
121. Crumling MA et al. **Hearing loss and hair cell death in mice given the cholesterol-chelating agent hydroxypropyl-β-cyclodextrin.** *PLoS One* 2012. — 🎯 외유모세포 사멸. 모델의 `cd_koto`·`OHC`. [PMID 23285273](https://pubmed.ncbi.nlm.nih.gov/23285273/)
122. Crumling MA et al. **Cyclodextrins and Iatrogenic Hearing Loss: New Drugs with Significant Risk.** *Front Cell Neurosci* 2017. [PMID 29163061](https://pubmed.ncbi.nlm.nih.gov/29163061/)
123. Ward S et al. **2-hydroxypropyl-beta-cyclodextrin raises hearing threshold in normal cats and in cats with Niemann-Pick type C disease.** *Pediatr Res* 2010. — 🎯 정상 동물에서도 발생 → 질환과 무관한 약물 독성. [PMID 20357695](https://pubmed.ncbi.nlm.nih.gov/20357695/)
124. Cronin S et al. **Hearing Loss and Otopathology Following Systemic and Intracerebroventricular Delivery of 2-Hydroxypropyl-Beta-Cyclodextrin.** *J Assoc Res Otolaryngol* 2015. — 🎯 **투여 경로별 차이**. 모델이 이독성을 CSF 농도(와우수도관 경유)에 걸고 정맥 경로와 분리한 근거. [PMID 26055150](https://pubmed.ncbi.nlm.nih.gov/26055150/)
125. Takahashi S et al. **Susceptibility of outer hair cells to cholesterol chelator 2-hydroxypropyl-β-cyclodextrine is prestin-dependent.** *Sci Rep* 2016. — 🎯 prestin 의존성 → **효능(막 콜레스테롤 추출)과 독성이 같은 기전**. [PMID 26903308](https://pubmed.ncbi.nlm.nih.gov/26903308/)
126. Zhou Y et al. **The susceptibility of cochlear outer hair cells to cyclodextrin is not related to their electromotile activity.** *Acta Neuropathol Commun* 2018. — 🔍 위 논문과 상반되는 결과. 기전 세부는 미해결임을 모델에 명시. [PMID 30249300](https://pubmed.ncbi.nlm.nih.gov/30249300/)
127. Matsuo M et al. **Effects of intracerebroventricular administration of 2-hydroxypropyl-β-cyclodextrin in a patient with Niemann-Pick Type C disease.** *Mol Genet Metab Rep* 2014. [PMID 27896112](https://pubmed.ncbi.nlm.nih.gov/27896112/)
128. El-Darzi N et al. **2-Hydroxypropyl-β-cyclodextrin reduces retinal cholesterol in wild-type and Cyp27a1−/−Cyp46a1−/− mice.** *Br J Pharmacol* 2021. — 표적 외 조직에서도 콜레스테롤을 뽑는다. [PMID 32698250](https://pubmed.ncbi.nlm.nih.gov/32698250/)

---

## 12. 개발 중 · 실험적 접근 (Emerging approaches)

모델 클러스터 20의 근거.

129. Chandler RJ et al. **Systemic AAV9 gene therapy improves the lifespan of mice with Niemann-Pick disease, type C1.** *Hum Mol Genet* 2017. [PMID 27798114](https://pubmed.ncbi.nlm.nih.gov/27798114/)
130. Pipalia NH et al. **Histone deacetylase inhibitor treatment dramatically reduces cholesterol accumulation in Niemann-Pick type C1 mutant human fibroblasts.** *Proc Natl Acad Sci U S A* 2011. [PMID 21436030](https://pubmed.ncbi.nlm.nih.gov/21436030/)
131. Pipalia NH et al. **Histone deacetylase inhibitors correct the cholesterol storage defect in most Niemann-Pick C1 mutant cells.** *J Lipid Res* 2017. — 🔍 **변이별 반응 차이** → 모델의 유전형×기전 상호작용을 뒷받침. [PMID 28193631](https://pubmed.ncbi.nlm.nih.gov/28193631/)
132. Pipalia NH et al. **HSP90 inhibitors reduce cholesterol storage in Niemann-Pick type C1 mutant fibroblasts.** *J Lipid Res* 2021. — 샤페론 경로의 대체 접근. [PMID 34481829](https://pubmed.ncbi.nlm.nih.gov/34481829/)
133. Maceyka M, Spiegel S. **The potential of histone deacetylase inhibitors in Niemann-Pick type C disease.** *FEBS J* 2013. [PMID 23992240](https://pubmed.ncbi.nlm.nih.gov/23992240/)
134. Hovakimyan M et al. **Combined therapy with cyclodextrin/allopregnanolone and miglustat improves motor but not cognitive functions in Niemann-Pick Type C1 mice.** *Neuroscience* 2013. — 🔍 병용의 상가/시너지 여부. [PMID 23948640](https://pubmed.ncbi.nlm.nih.gov/23948640/)
135. Maass F et al. **Reduced cerebellar neurodegeneration after combined therapy with cyclodextrin/allopregnanolone and miglustat in NPC1 mutant mice.** *J Neurosci Res* 2015. [PMID 25400034](https://pubmed.ncbi.nlm.nih.gov/25400034/)
136. Ebner L et al. **Evaluation of Two Liver Treatment Strategies in a Mouse Model of Niemann-Pick-Disease Type C1.** *Int J Mol Sci* 2018. [PMID 29587349](https://pubmed.ncbi.nlm.nih.gov/29587349/)

---

## 13. QSP 방법론 (QSP methodology)

137. Elmokadem A et al. **Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.** *CPT Pharmacometrics Syst Pharmacol* 2019. — mrgsolve 표준 인용. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
138. Helmlinger G et al. **Quantitative Systems Pharmacology: An Exemplar Model-Building Workflow With Applications in Cardiovascular, Metabolic, and Oncology Drug Development.** *CPT Pharmacometrics Syst Pharmacol* 2019. [PMID 31087533](https://pubmed.ncbi.nlm.nih.gov/31087533/)
139. Allen RJ et al. **Efficient Generation and Selection of Virtual Populations in Quantitative Systems Pharmacology Models.** *CPT Pharmacometrics Syst Pharmacol* 2016. — Shiny 앱의 가상 환자군 생성 방식. [PMID 27069777](https://pubmed.ncbi.nlm.nih.gov/27069777/)
140. Rieger TR et al. **Improving the generation and selection of virtual populations in quantitative systems pharmacology models.** *Prog Biophys Mol Biol* 2018. [PMID 29902482](https://pubmed.ncbi.nlm.nih.gov/29902482/)
141. Bai JPF et al. **FDA-Industry Scientific Exchange on assessing quantitative systems pharmacology models in clinical drug development.** *AAPS J* 2021. [PMID 33931790](https://pubmed.ncbi.nlm.nih.gov/33931790/)
142. Bai JPF et al. **Translational Quantitative Systems Pharmacology in Drug Development: from Current Landscape to Good Practices.** *AAPS J* 2019. [PMID 31161268](https://pubmed.ncbi.nlm.nih.gov/31161268/)
143. Braniff N et al. **An integrated quantitative systems pharmacology virtual population approach for calibration with oncology efficacy data.** *CPT Pharmacometrics Syst Pharmacol* 2025. [PMID 39508122](https://pubmed.ncbi.nlm.nih.gov/39508122/)
144. Traynard P et al. **Logic Modeling in Quantitative Systems Pharmacology.** *CPT Pharmacometrics Syst Pharmacol* 2017. [PMID 28681552](https://pubmed.ncbi.nlm.nih.gov/28681552/)

---

## 정량적 보정 목표 요약 (Calibration targets, at a glance)

| # | 목표 | 출처 | 값 | 용도 |
|---|------|------|-----|------|
| T1 | 혈장 C-triol, 환자 | PMID 33228797 | 88.31 ng/mL | 🎯 보정 (`ktri`, `Ktri`) |
| T2 | 혈장 C-triol, 대조 | PMID 33228797 | 5.97 ng/mL | 🎯 보정 (`ktri`, `Ktri`) |
| T3 | 5영역 NPCCSS 진행 | PMID 33228797 | 1.5 점/년 | 🎯 보정 (`kdie`) |
| T4 | 17영역 NPCCSS 진행 | PMID 33228797 | 2.7–2.9 점/년 | 🔍 예측 검증 |
| T5 | 아리모클로몰 12개월 5영역 차 | PMID 34418116 | −1.40 점 | 🎯 보정 (`Emax_fold`) |
| T6 | 아리모클로몰 12개월 4영역 차 | PMID 40520915 | −1.70 점 | 🔍 예측 검증 |
| T7 | 아리모클로몰 위약군 진행 | PMID 34418116 | +2.11–2.15 점/년 | 🔍 예측 검증 |
| T8 | 레바세틸류신 12주 SARA 차 | PMID 38294974 | −1.28 점 | 🎯 보정 (`nal_Emax_sym`) |
| T9 | C-triol vs 5영역 Spearman ρ | PMID 33228797 | 0.265 | 🔍 예측 검증 |
| T10 | 아드라베타덱스 공개라벨 진행 | PMID 28803710 | 7/14 vs 21/21 | 🔍 예측 검증 |
| T11 | 기저 SARA (IB1001-301) | PMID 38294974 | 15.91 ± 7.65 | 🎯 보정 (`sara_k`) |
| T12 | 레바세틸류신 18개월 SARA 변화 | PMID 40513057 | −1.64 점 | 🔍 검증 — **불일치** |
| T13 | 레바세틸류신 12개월 SARA 변화 | PMID 40513057 | −1.88 점 | 🔍 검증 — **불일치** |
| T14 | 진행의 발병연령 무관 선형성 | PMID 19415691 | — | 🔍 구조적 검증 |
| T15 | 아리모클로몰 미글루스타트 하위군 | PMID 34418116 | −2.06 점 | 🔍 예측 검증 |

보정에 쓰인 목표는 7개(T1·T2·T3·T8·T11 및 이들과 결합된 자연사 앵커), 검증에만 쓰인 목표는 8개입니다. 검증 목표 중 T5·T6·T12·T13·T15는 **재현되지 않았고**, 조정하지 않고 그대로 보고했습니다. 결과는
[README.md](README.md)의 "모델–데이터 대조표"에 있으며, **맞지 않은 항목도
조정하지 않고 그대로 보고**했습니다.

---

## 면책 (Disclaimer)

본 참고문헌 목록과 모델은 교육·연구 목적입니다. 개별 환자의 진단·치료 결정에
사용해서는 안 됩니다. 인용된 수치는 원문에서 확인 가능한 값이며, 모델 파라미터는
그 수치를 재현하도록 조정된 근사치입니다.
