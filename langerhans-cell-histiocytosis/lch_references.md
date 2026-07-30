# 랑게르한스 세포 조직증식증 (Langerhans Cell Histiocytosis, LCH) — 참고문헌

이 모델(`lch_qsp_model.dot`, `lch_mrgsolve_model.R`, `lch_shiny_app.R`)의 각 구조·파라미터
근거가 되는 문헌 목록입니다.

**검증 방법 (verification method).** 모든 항목의 PMID는 NCBI E-utilities
(`esearch` → `esummary`)로 조회하여, 반환된 제목을 질의 제목과 토큰 중첩 비율로 대조한 뒤
확정했습니다. 대조에 실패한 후보는 목록에서 제외했습니다 — 즉 아래 PMID는 모두 실재하는
레코드이며 표시된 저자·저널·연도와 일치합니다.

---

## A. 클론 기원과 드라이버 변이 (Clonality & MAPK Driver Lesions)

모델 대응: `cluster_driver`, 파라미터 `GENO`

1. Willman CL, Busque L, Griffith BB, et al. **Langerhans'-cell histiocytosis (histiocytosis X) — a clonal proliferative disease.** N Engl J Med. 1994. [PMID 8008029](https://pubmed.ncbi.nlm.nih.gov/8008029/) — LCH가 반응성 증식이 아니라 클론성 질환임을 확립한 논문.
2. Badalian-Very G, Vergilio JA, Degar BA, et al. **Recurrent BRAF mutations in Langerhans cell histiocytosis.** Blood. 2010. [PMID 20519626](https://pubmed.ncbi.nlm.nih.gov/20519626/) — BRAF V600E ~57% 최초 보고. 모델의 `GENO=1` 경로 근거.
3. Chakraborty R, Hampton OA, Shen X, et al. **Mutually exclusive recurrent somatic mutations in MAP2K1 and BRAF support a central role for ERK activation in LCH pathogenesis.** Blood. 2014. [PMID 25202140](https://pubmed.ncbi.nlm.nih.gov/25202140/) — "모든 길은 pERK로 통한다"는 모델의 중심 노드 근거.
4. Brown NA, Furtado LV, Betz BL, et al. **High prevalence of somatic MAP2K1 mutations in BRAF V600E-negative Langerhans cell histiocytosis.** Blood. 2014. [PMID 24982505](https://pubmed.ncbi.nlm.nih.gov/24982505/) — `GENO=2`(BRAFi 무효, MEKi 유효) 시나리오 근거.
5. Nelson DS, Quispel W, Badalian-Very G, et al. **Somatic activating ARAF mutations in Langerhans cell histiocytosis.** Blood. 2014. [PMID 24652991](https://pubmed.ncbi.nlm.nih.gov/24652991/)
6. Nelson DS, van Halteren A, Quispel WT, et al. **MAP2K1 and MAP3K1 mutations in Langerhans cell histiocytosis.** Genes Chromosomes Cancer. 2015. [PMID 25899310](https://pubmed.ncbi.nlm.nih.gov/25899310/)
7. Liu R, et al. **Somatic ARAF mutations in pediatric Langerhans cell histiocytosis: clinicopathological features.** Clin Exp Med. 2023. [PMID 37572153](https://pubmed.ncbi.nlm.nih.gov/37572153/)
8. Diamond EL, Durham BH, Haroche J, et al. **Diverse and targetable kinase alterations drive histiocytic neoplasms.** Cancer Discov. 2016. [PMID 26566875](https://pubmed.ncbi.nlm.nih.gov/26566875/)
9. Durham BH, Lopez Rodrigo E, Picarsic J, et al. **Activating mutations in CSF1R and additional receptor tyrosine kinases in histiocytic neoplasms.** Nat Med. 2019. [PMID 31768065](https://pubmed.ncbi.nlm.nih.gov/31768065/) — `GF_RTK` 노드.
10. Mourah S, How-Kit A, Meignin V, et al. **Recurrent NRAS mutations in pulmonary Langerhans cell histiocytosis.** Eur Respir J. 2016. [PMID 27076591](https://pubmed.ncbi.nlm.nih.gov/27076591/) — 폐 LCH의 별개 유전형.
11. Jouenne F, Chevret S, Bugnet E, et al. **Genetic landscape of adult Langerhans cell histiocytosis with lung involvement.** Eur Respir J. 2020. [PMID 31806714](https://pubmed.ncbi.nlm.nih.gov/31806714/)
12. Emile JF, Abla O, Fraitag S, et al. **Revised classification of histiocytoses and neoplasms of the macrophage-dendritic cell lineages.** Blood. 2016. [PMID 26966089](https://pubmed.ncbi.nlm.nih.gov/26966089/)
13. Kemps PG, Picarsic J, Durham BH, et al. **ALK-positive histiocytosis: a new clinicopathologic spectrum highlighting neurologic involvement.** Blood. 2022. [PMID 34727172](https://pubmed.ncbi.nlm.nih.gov/34727172/)

## B. 기원세포 분화단계 — 모델의 첫 번째 구조적 선택 (Cell of Origin)

모델 대응: `cluster_origin`, 파라미터 `PRECM0`, `THB/THS/THR/THP/THC/THL`

14. Berres ML, Lim KP, Peters T, et al. **BRAF-V600E expression in precursor versus differentiated dendritic cells defines clinically distinct LCH risk groups.** J Exp Med. 2014. [PMID 25646268](https://pubmed.ncbi.nlm.nih.gov/25646268/) — **모델의 파티션(θ) 구조를 직접 뒷받침하는 핵심 논문.** 변이가 발생한 분화단계가 질환 범위를 결정.
15. Milne P, Bigley V, Bacon CM, et al. **Hematopoietic origin of Langerhans cell histiocytosis and Erdheim-Chester disease in adults.** Blood. 2017. [PMID 28512190](https://pubmed.ncbi.nlm.nih.gov/28512190/)
16. Xiao Y, van Halteren AGS, Lei X, et al. **Bone marrow-derived myeloid progenitors as driver mutation carriers in high- and low-risk Langerhans cell histiocytosis.** Blood. 2020. [PMID 32750121](https://pubmed.ncbi.nlm.nih.gov/32750121/) — 골수 전구세포 저장고(`PRECM`)의 존재 근거.
17. Badalian-Very G. **A common progenitor cell in LCH and ECD.** Blood. 2014. [PMID 25124781](https://pubmed.ncbi.nlm.nih.gov/25124781/)
18. Halbritter F, Farlik M, Schwentner R, et al. **Epigenomics and single-cell sequencing define a developmental hierarchy in Langerhans cell histiocytosis.** Cancer Discov. 2019. [PMID 31345789](https://pubmed.ncbi.nlm.nih.gov/31345789/)
19. Gruber TA. **Single-cell RNA sequencing reveals a developmental hierarchy in Langerhans cell histiocytosis.** Cancer Discov. 2019. [PMID 31575563](https://pubmed.ncbi.nlm.nih.gov/31575563/)
20. Shi H, He H, Cui L, et al. **Transcriptomic landscape of circulating mononuclear phagocytes in Langerhans cell histiocytosis at single-cell level.** Blood. 2021. [PMID 34132762](https://pubmed.ncbi.nlm.nih.gov/34132762/) — `CIRC` 구획(순환 변이 전구세포).
21. Allen CE, Li L, Peters TL, et al. **Cell-specific gene expression in Langerhans cell histiocytosis lesions reveals a distinct profile compared with epidermal Langerhans cells.** J Immunol. 2010. [PMID 20220088](https://pubmed.ncbi.nlm.nih.gov/20220088/) — 지도에서 `Epidermal_LC`를 점선(비-기원)으로 그린 근거.
22. Schwentner R, Jug G, Kauer MO, et al. **JAG2 signaling induces differentiation of CD14+ monocytes into Langerhans cell histiocytosis-like cells.** J Leukoc Biol. 2019. [PMID 30296338](https://pubmed.ncbi.nlm.nih.gov/30296338/)
23. Kvedaraite E, Milne P, Khalilnezhad A, et al. **Notch-dependent cooperativity between myeloid lineages promotes Langerhans cell histiocytosis pathology.** Sci Immunol. 2022. [PMID 36525505](https://pubmed.ncbi.nlm.nih.gov/36525505/)

## C. ERK 출력 · 세네센스 · 병변 내 포획 (MAPK Output, Senescence, DC Trapping)

모델 대응: `ERK`, `CCND`, `BCL`, `SASP`, `PFMAX`

24. Bollag G, Hirth P, Tsai J, et al. **Clinical efficacy of a RAF inhibitor needs broad target blockade in BRAF-mutant melanoma.** Nature. 2010. [PMID 20823850](https://pubmed.ncbi.nlm.nih.gov/20823850/) — **pERK를 80% 이상 억제해야 반응이 나타난다**는 관찰. 모델의 Hill 계수(`HC=2`, `EC50C=0.35`)가 만드는 문턱의 근거.
25. Bigenwald C, Le Berichel J, Wilk CM, et al. **BRAF V600E-induced senescence drives Langerhans cell histiocytosis pathophysiology.** Nat Med. 2021. [PMID 33958797](https://pubmed.ncbi.nlm.nih.gov/33958797/) — **`SASP` 구획의 직접 근거.** 낮은 증식률과 높은 사이토카인 분비가 공존하는 이유.
26. Hogstad B, Berres ML, Chakraborty R, et al. **RAF/MEK/extracellular signal-related kinase pathway suppresses dendritic cell migration and traps dendritic cells in Langerhans cell histiocytosis lesions.** J Exp Med. 2018. [PMID 29263218](https://pubmed.ncbi.nlm.nih.gov/29263218/) — `CCR6_CCR7` "떠날 수 없어 축적된다" 노드.
27. Brabencova E, Tazi A, Lorenzato M, et al. **Langerhans cells in Langerhans cell granulomatosis are not actively proliferating cells.** Am J Pathol. 1998. [PMID 9588881](https://pubmed.ncbi.nlm.nih.gov/9588881/) — **`PFMAX=0.09`(낮은 증식분획)의 근거이자, S기 의존 약물의 한계와 클라드리빈의 우위를 설명하는 축.**
28. Price JG, Idoyaga J, Salmon H, et al. **CDKN1A regulates Langerhans cell survival and promotes Treg cell generation upon exposure to ionizing irradiation.** Nat Immunol. 2015. [PMID 26343536](https://pubmed.ncbi.nlm.nih.gov/26343536/) — `p16`/`p21` 축.
29. Olsson Åkefeldt S, Ismail MB, Valentin H, et al. **Targeting BCL2 family in human myeloid dendritic cells: a challenge to cure diseases with chronic inflammations associated with bone loss.** Clin Dev Immunol. 2013. [PMID 23762095](https://pubmed.ncbi.nlm.nih.gov/23762095/) — `BCL2A1` 생존신호 노드.
30. Annels NE, da Costa CE, Prins FA, et al. **Aberrant chemokine receptor expression and chemokine production by Langerhans cells underlies the pathogenesis of Langerhans cell histiocytosis.** J Exp Med. 2003. [PMID 12743170](https://pubmed.ncbi.nlm.nih.gov/12743170/) — `CCL20`-CCR6 자가증폭 루프(`FREC`).

## D. 병변 면역 미세환경과 사이토카인 (Lesional Microenvironment & Secretome)

모델 대응: `cluster_micro`, `cluster_secretome`, `IL1B/TNFA/IL6/OSM/MMP9/RANKL`, `TREG`

31. Egeler RM, Favara BE, van Meurs M, et al. **Differential in situ cytokine profiles of Langerhans-like cells and T cells in Langerhans cell histiocytosis: abundant expression of cytokines relevant to disease and treatment.** Blood. 1999. [PMID 10590064](https://pubmed.ncbi.nlm.nih.gov/10590064/) — 분비체 목록의 1차 근거.
32. Arenzana-Seisdedos F, Barbey S, Virelizier JL, et al. **Histiocytosis X. Purified (T6+) cells from bone granuloma produce interleukin 1 and prostaglandin E2 in culture.** J Clin Invest. 1986. [PMID 2418061](https://pubmed.ncbi.nlm.nih.gov/2418061/) — IL-1이 골 병변에서 직접 생산됨을 보인 고전 논문.
33. Rosso DA, Ripoli MF, Roy A, et al. **Serum levels of interleukin-1 receptor antagonist and tumor necrosis factor-alpha are elevated in children with Langerhans cell histiocytosis.** J Pediatr Hematol Oncol. 2003. [PMID 12794527](https://pubmed.ncbi.nlm.nih.gov/12794527/)
34. Senechal B, Elain G, Jeziorski E, et al. **Expansion of regulatory T cells in patients with Langerhans cell histiocytosis.** PLoS Med. 2007. [PMID 17696642](https://pubmed.ncbi.nlm.nih.gov/17696642/) — `TREG`가 면역 청소(`KIMM`)를 억제하는 구조.
35. Sengal A, Velazquez J, Hahne M, et al. **Overcoming T-cell exhaustion in LCH: PD-1 blockade and targeted MAPK inhibition are synergistic in a mouse model of LCH.** Blood. 2021. [PMID 33075814](https://pubmed.ncbi.nlm.nih.gov/33075814/)
36. Coury F, Annels N, Rivollier A, et al. **Langerhans cell histiocytosis reveals a new IL-17A-dependent pathway of dendritic cell fusion.** Nat Med. 2008. [PMID 18157139](https://pubmed.ncbi.nlm.nih.gov/18157139/) — `IL17ON` 스위치(논쟁 중인 경로를 켜고 끌 수 있게 둔 이유).
37. da Costa CE, Annels NE, Faaij CM, et al. **Presence of osteoclast-like multinucleated giant cells in the bone and nonostotic lesions of Langerhans cell histiocytosis.** J Exp Med. 2005. [PMID 15753204](https://pubmed.ncbi.nlm.nih.gov/15753204/) — `MNGC` → `OCL`.
38. Rivollier A, Mazzorana M, Tebib J, et al. **Immature dendritic cell transdifferentiation into osteoclasts: a novel pathway sustained by the rheumatoid arthritis microenvironment.** Blood. 2004. [PMID 15308576](https://pubmed.ncbi.nlm.nih.gov/15308576/)
39. Grosjean F, Nasi S, Schneider P, et al. **Dendritic cells cause bone lesions in a new mouse model of histiocytosis.** PLoS One. 2015. [PMID 26247358](https://pubmed.ncbi.nlm.nih.gov/26247358/)

## E. 골 대사 · RANKL/OPG (Bone Remodelling)

모델 대응: `RANKL`, `OCL`, `BVOL`, `KRES`, `KHEAL`

40. Makras P, Tsoli M, Anastasilakis AD, et al. **Serum osteoprotegerin, RANKL, and Dkk-1 levels in adults with Langerhans cell histiocytosis.** J Clin Endocrinol Metab. 2012. [PMID 22278426](https://pubmed.ncbi.nlm.nih.gov/22278426/) — `RANKL0`, `KR_RANKL` 보정.
41. Makras P, Papadogias D, Kaltsas G, et al. **Rationale for the application of RANKL inhibition in the treatment of Langerhans cell histiocytosis.** J Clin Endocrinol Metab. 2015. [PMID 25375981](https://pubmed.ncbi.nlm.nih.gov/25375981/)
42. Makras P, Tsoli M, Kaltsas G. **Denosumab for the treatment of adult multisystem Langerhans cell histiocytosis.** Metabolism. 2017. [PMID 28285639](https://pubmed.ncbi.nlm.nih.gov/28285639/)
43. Anastasilakis AD, Tsoli M, Kaltsas G, Makras P. **Bone metabolism in Langerhans cell histiocytosis.** Endocr Connect. 2018. [PMID 29967185](https://pubmed.ncbi.nlm.nih.gov/29967185/)

## F. 임상 스펙트럼 · 분류 · 예후 (Clinical Spectrum, Classification, Prognosis)

모델 대응: `cluster_class`, `DAS`, `Wk6_response`

44. Allen CE, Merad M, McClain KL. **Langerhans-cell histiocytosis.** N Engl J Med. 2018. [PMID 30157397](https://pubmed.ncbi.nlm.nih.gov/30157397/)
45. Rodriguez-Galindo C, Allen CE. **Langerhans cell histiocytosis.** Blood. 2020. [PMID 32106306](https://pubmed.ncbi.nlm.nih.gov/32106306/)
46. McClain KL, Bigenwald C, Collin M, et al. **Histiocytic disorders.** Nat Rev Dis Primers. 2021. [PMID 34620874](https://pubmed.ncbi.nlm.nih.gov/34620874/)
47. Badalian-Very G, Vergilio JA, Fleming M, Rollins BJ. **Pathogenesis of Langerhans cell histiocytosis.** Annu Rev Pathol. 2013. [PMID 22906202](https://pubmed.ncbi.nlm.nih.gov/22906202/)
48. Haupt R, Minkov M, Astigarraga I, et al. **Langerhans cell histiocytosis (LCH): guidelines for diagnosis, clinical work-up, and treatment for patients till the age of 18 years.** Pediatr Blood Cancer. 2013. [PMID 23109216](https://pubmed.ncbi.nlm.nih.gov/23109216/) — RO+/RO−, special site 정의.
49. Goyal G, Tazi A, Go RS, et al. **International expert consensus recommendations for the diagnosis and treatment of Langerhans cell histiocytosis in adults.** Blood. 2022. [PMID 35271698](https://pubmed.ncbi.nlm.nih.gov/35271698/)
50. Donadieu J, Piguet C, Bernard F, et al. **A new clinical score for disease activity in Langerhans cell histiocytosis.** Pediatr Blood Cancer. 2004. [PMID 15390280](https://pubmed.ncbi.nlm.nih.gov/15390280/) — **`DAS` 출력의 근거 척도.**
51. Minkov M, Grois N, Heitger A, et al. **Response to initial treatment of multisystem Langerhans cell histiocytosis: an important prognostic indicator.** Med Pediatr Oncol. 2002. [PMID 12376981](https://pubmed.ncbi.nlm.nih.gov/12376981/) — 6주 반응이 최강 예후인자.
52. Simko SJ, Garmezy B, Abhyankar H, et al. **Differentiating skin-limited and multisystem Langerhans cell histiocytosis.** J Pediatr. 2014. [PMID 25441388](https://pubmed.ncbi.nlm.nih.gov/25441388/)
53. Rigaud C, Barkaoui MA, Thomas C, et al. **Langerhans cell histiocytosis: therapeutic strategy and outcome in a 30-year nationwide cohort of 1478 patients under 18 years of age.** Br J Haematol. 2016. [PMID 27273725](https://pubmed.ncbi.nlm.nih.gov/27273725/) — 재활성화율·후유증 빈도 보정.
54. Liu H, Stiller CA, Crooks CJ, et al. **Incidence, prevalence and survival in patients with Langerhans cell histiocytosis: a national registry study.** Br J Haematol. 2022. [PMID 36122574](https://pubmed.ncbi.nlm.nih.gov/36122574/)
55. Héritier S, Emile JF, Barkaoui MA, et al. **BRAF mutation correlates with high-risk Langerhans cell histiocytosis and increased resistance to first-line therapy.** J Clin Oncol. 2016. [PMID 27382093](https://pubmed.ncbi.nlm.nih.gov/27382093/)
56. Thalhammer J, et al. **Childhood Langerhans cell histiocytosis hematological involvement: severity associated with BRAF V600E.** Blood. 2025. [PMID 39486044](https://pubmed.ncbi.nlm.nih.gov/39486044/)
57. Aricò M, Astigarraga I, Braier J, et al. **Lack of bone lesions at diagnosis is associated with inferior outcome in multisystem Langerhans cell histiocytosis of childhood.** Br J Haematol. 2015. [PMID 25522229](https://pubmed.ncbi.nlm.nih.gov/25522229/)
58. Minkov M. **Multisystem Langerhans cell histiocytosis in children: current treatment and future directions.** Paediatr Drugs. 2011. [PMID 21351807](https://pubmed.ncbi.nlm.nih.gov/21351807/)

## G. 1차 화학요법 임상시험 (Front-Line Chemotherapy Trials)

모델 대응: 시나리오 S2/S3, `r_vblpred()`

59. Gadner H, Grois N, Arico M, et al. **A randomized trial of treatment for multisystem Langerhans' cell histiocytosis.** J Pediatr. 2001. [PMID 11343051](https://pubmed.ncbi.nlm.nih.gov/11343051/) — LCH-II.
60. Gadner H, Grois N, Pötschger U, et al. **Improved outcome in multisystem Langerhans cell histiocytosis is associated with therapy intensification.** Blood. 2008. [PMID 18089850](https://pubmed.ncbi.nlm.nih.gov/18089850/)
61. Gadner H, Minkov M, Grois N, et al. **Therapy prolongation improves outcome in multisystem Langerhans cell histiocytosis.** Blood. 2013. [PMID 23589673](https://pubmed.ncbi.nlm.nih.gov/23589673/) — LCH-III: 12개월 치료가 재활성화를 줄임 → 모델의 `MNTC` 유지요법 항.
62. Morimoto A, Shioda Y, Imamura T, et al. **Intensified and prolonged therapy comprising cytarabine, vincristine and prednisolone for Langerhans cell histiocytosis (JLSG-02).** Int J Hematol. 2016. [PMID 27040279](https://pubmed.ncbi.nlm.nih.gov/27040279/)

## H. 구제 요법 (Salvage Therapy)

모델 대응: 시나리오 S3/S4/S8, `r_cladarac()`

63. Donadieu J, Bernard F, van Noesel M, et al. **Cladribine and cytarabine in refractory multisystem Langerhans cell histiocytosis: results of an international phase 2 study.** Blood. 2015. [PMID 26194764](https://pubmed.ncbi.nlm.nih.gov/26194764/) — **2-CdA 9 mg/m²/d + Ara-C 500 mg/m²/d × 5일 q28일 용법과 반응률·혈액독성의 1차 근거.**
64. Bernard F, Thomas C, Bertrand Y, et al. **Multi-centre pilot study of 2-chlorodeoxyadenosine and cytosine arabinoside combined chemotherapy in refractory Langerhans cell histiocytosis with haematological dysfunction.** Eur J Cancer. 2005. [PMID 16291085](https://pubmed.ncbi.nlm.nih.gov/16291085/)
65. Simko SJ, Tran HD, Jones J, et al. **Clofarabine salvage therapy in refractory multifocal histiocytic disorders, including Langerhans cell histiocytosis, juvenile xanthogranuloma and Rosai-Dorfman disease.** Pediatr Blood Cancer. 2014. [PMID 24106153](https://pubmed.ncbi.nlm.nih.gov/24106153/)
66. Parekh D, et al. **Clofarabine monotherapy in aggressive, relapsed and refractory Langerhans cell histiocytosis.** Br J Haematol. 2024. [PMID 38501389](https://pubmed.ncbi.nlm.nih.gov/38501389/)
67. Veys PA, Nanduri V, Baker KS, et al. **Haematopoietic stem cell transplantation for refractory Langerhans cell histiocytosis: outcome by intensity of conditioning.** Br J Haematol. 2015. [PMID 25817915](https://pubmed.ncbi.nlm.nih.gov/25817915/)

## I. MAPK 표적치료와 중단 후 재발 — 모델의 두 번째 구조적 선택

모델 대응: 시나리오 S5/S6/S7/S8, `SL_MAPKI_KILL=0`, `Reservoir_persist` → `Rebound`

68. Haroche J, Cohen-Aubart F, Emile JF, et al. **Dramatic efficacy of vemurafenib in both multisystemic and refractory Erdheim-Chester disease and Langerhans cell histiocytosis harboring BRAF V600E mutation.** Blood. 2013. [PMID 23258922](https://pubmed.ncbi.nlm.nih.gov/23258922/)
69. Héritier S, Jehanne M, Leverger G, et al. **Vemurafenib use in an infant for high-risk Langerhans cell histiocytosis.** JAMA Oncol. 2015. [PMID 26180941](https://pubmed.ncbi.nlm.nih.gov/26180941/)
70. Donadieu J, Larabi IA, Tardieu M, et al. **Vemurafenib for refractory multisystem Langerhans cell histiocytosis in children: an international observational study.** J Clin Oncol. 2019. [PMID 31513482](https://pubmed.ncbi.nlm.nih.gov/31513482/) — **소아 20 mg/kg/d 용법, 며칠 내 임상반응, 중단 후 높은 재활성화율 — 모델의 세포정지성/재발 구조의 핵심 근거.**
71. Diamond EL, Subbiah V, Lockhart AC, et al. **Vemurafenib for BRAF V600-mutant Erdheim-Chester disease and Langerhans cell histiocytosis (VE-BASKET).** JAMA Oncol. 2018. [PMID 29188284](https://pubmed.ncbi.nlm.nih.gov/29188284/)
72. Subbiah V, Puzanov I, Blay JY, et al. **Pan-cancer efficacy of vemurafenib in BRAF V600-mutant non-melanoma cancers.** Cancer Discov. 2020. [PMID 32029534](https://pubmed.ncbi.nlm.nih.gov/32029534/)
73. Diamond EL, Durham BH, Ulaner GA, et al. **Efficacy of MEK inhibition in patients with histiocytic neoplasms.** Nature. 2019. [PMID 30867592](https://pubmed.ncbi.nlm.nih.gov/30867592/) — 유전형과 무관한 MEKi 효과 → `MEKi_genotype` 노드.
74. Eckstein OS, Visser J, Rodriguez-Galindo C, Allen CE. **Clinical responses and persistent BRAF V600E+ blood cells in children with LCH treated with MAPK pathway inhibition.** Blood. 2019. [PMID 30718231](https://pubmed.ncbi.nlm.nih.gov/30718231/) — **MAPKi 하에서 변이 세포가 사라지지 않는다는 직접 증거(`Reservoir_persist`).**
75. Evseev D, Kalinina I, Raykina E, et al. **Vemurafenib provides a rapid and robust clinical response in pediatric Langerhans cell histiocytosis with the BRAF V600E mutation but does not eliminate low-level minimal residual disease.** Int J Hematol. 2021. [PMID 34383272](https://pubmed.ncbi.nlm.nih.gov/34383272/) — cfDNA/MRD 플래토의 근거.
76. Cohen Aubart F, Emile JF, Carrat F, et al. **Targeted therapies in 54 patients with Erdheim-Chester disease, including follow-up after interruption (the LOVE study).** Blood. 2017. [PMID 28667012](https://pubmed.ncbi.nlm.nih.gov/28667012/) — **중단 후 재발률의 정량적 근거.**
77. Whitlock JA, Geoerger B, Dunkel IJ, et al. **Dabrafenib, alone or in combination with trametinib, in BRAF V600-mutated pediatric Langerhans cell histiocytosis.** Blood Adv. 2023. [PMID 36884302](https://pubmed.ncbi.nlm.nih.gov/36884302/) — 소아 dabrafenib/trametinib 용량.
78. Cournoyer E, Ferrell J, Sharp S, et al. **Dabrafenib and trametinib in Langerhans cell histiocytosis and other histiocytic disorders.** Haematologica. 2024. [PMID 37731389](https://pubmed.ncbi.nlm.nih.gov/37731389/)
79. Donadieu J, et al. **Long-term MAPK inhibition of childhood refractory Langerhans cell histiocytosis: an observational study.** Blood Adv. 2026. [PMID 41678955](https://pubmed.ncbi.nlm.nih.gov/41678955/)
80. Lorillon G, Jouenne F, Baroudjian B, et al. **Response to trametinib of a pulmonary Langerhans cell histiocytosis harboring a MAP2K1 deletion.** Am J Respir Crit Care Med. 2018. [PMID 29694792](https://pubmed.ncbi.nlm.nih.gov/29694792/)
81. Tardieu M, Néron A, Duvert-Lehembre S, et al. **Cutaneous adverse events in children treated with vemurafenib for refractory BRAF V600E mutated Langerhans cell histiocytosis.** Pediatr Blood Cancer. 2021. [PMID 34109735](https://pubmed.ncbi.nlm.nih.gov/34109735/) — `SKTOX`(역설적 피부독성) 보정.
82. Long GV, Stroyakovskiy D, Gogas H, et al. **Combined BRAF and MEK inhibition versus BRAF inhibition alone in melanoma.** N Engl J Med. 2014. [PMID 25265492](https://pubmed.ncbi.nlm.nih.gov/25265492/) — 병용이 역설 활성화와 피부독성을 줄인다는 근거.

## J. 바이오마커 · 반응평가 (Biomarkers & Response Assessment)

모델 대응: `CFDNA`, `SCD163`, `CRP`, `PET`

83. Héritier S, Hélias-Rodzewicz Z, Lapillonne H, et al. **Circulating cell-free BRAF V600E as a biomarker in children with Langerhans cell histiocytosis.** Br J Haematol. 2017. [PMID 28444728](https://pubmed.ncbi.nlm.nih.gov/28444728/) — **`CFDNA` 구획의 직접 근거.**
84. Shimizu S, et al. **Detection of BRAF V600E mutation in radiological Langerhans cell histiocytosis.** Int J Hematol. 2023. [PMID 37010809](https://pubmed.ncbi.nlm.nih.gov/37010809/)
85. Ji X, et al. **18F-FDG PET/CT in pediatric Langerhans cell histiocytosis: relation to BRAF V600E mutation and risk stratification.** Eur J Radiol. 2025. [PMID 40819625](https://pubmed.ncbi.nlm.nih.gov/40819625/)
86. Baratto L, Nyalakonda R, Hawk KE, et al. **Comparison of whole-body DW-MRI with 2-[18F]FDG PET for staging and treatment monitoring of children with Langerhans cell histiocytosis.** Eur J Nucl Med Mol Imaging. 2023. [PMID 36717409](https://pubmed.ncbi.nlm.nih.gov/36717409/)

## K. CNS · 내분비 영구 후유증 — 모델의 세 번째 구조적 선택

모델 대응: `AVPN`, `ANTPIT`, `NEUR`, `TTET`, 시나리오 S9/S9b

87. Grois N, Fahrner B, Arceci RJ, et al. **Central nervous system disease in Langerhans cell histiocytosis.** J Pediatr. 2010. [PMID 20434166](https://pubmed.ncbi.nlm.nih.gov/20434166/)
88. Grois N, Pötschger U, Prosch H, et al. **Risk factors for diabetes insipidus in Langerhans cell histiocytosis.** Pediatr Blood Cancer. 2006. [PMID 16047354](https://pubmed.ncbi.nlm.nih.gov/16047354/) — **두개안면 CNS-risk 병변과 다장기 침범이 CDI 위험을 3–4배 높인다 → 지도의 `Site_CNSrisk → Stalk_thick` 엣지와 `THP` 파티션.**
89. Prosch H, Grois N, Prayer D, et al. **Central diabetes insipidus as presenting symptom of Langerhans cell histiocytosis.** Pediatr Blood Cancer. 2004. [PMID 15382278](https://pubmed.ncbi.nlm.nih.gov/15382278/)
90. Abla O, Weitzman S, Minkov M, et al. **Diabetes insipidus in Langerhans cell histiocytosis: when is treatment indicated?** Pediatr Blood Cancer. 2009. [PMID 19142995](https://pubmed.ncbi.nlm.nih.gov/19142995/) — 진단 지연이 CDI 확정에 미치는 영향(`TTET`).
91. Donadieu J, Rolon MA, Thomas C, et al. **Endocrine involvement in pediatric-onset Langerhans' cell histiocytosis: a population-based study.** J Pediatr. 2004. [PMID 15001940](https://pubmed.ncbi.nlm.nih.gov/15001940/) — CDI 누적 발생률 및 전엽 결손 동반율(≈50%).
92. Vaiani E, Felizzia G, Lubieniecki F, et al. **Paediatric Langerhans cell histiocytosis disease: long-term sequelae in the hypothalamic endocrine system.** Horm Res Paediatr. 2021. [PMID 34167121](https://pubmed.ncbi.nlm.nih.gov/34167121/)
93. Kurtulmus N, Mert M, Tanakol R, Yarman S. **The pituitary gland in patients with Langerhans cell histiocytosis: a clinical and radiological evaluation.** Endocrine. 2015. [PMID 25209890](https://pubmed.ncbi.nlm.nih.gov/25209890/)
94. Haupt R, Nanduri V, Calevo MG, et al. **Permanent consequences in Langerhans cell histiocytosis patients: a pilot study from the Histiocyte Society–Late Effects Study Group.** Pediatr Blood Cancer. 2004. [PMID 15049016](https://pubmed.ncbi.nlm.nih.gov/15049016/) — **영구 후유증 빈도(≈30–40%)의 근거.**
95. Mittheisz E, Seidl R, Prayer D, et al. **Central nervous system-related permanent consequences in patients with Langerhans cell histiocytosis.** Pediatr Blood Cancer. 2007. [PMID 16470521](https://pubmed.ncbi.nlm.nih.gov/16470521/)
96. Wnorowski M, Prosch H, Prayer D, et al. **Pattern and course of neurodegeneration in Langerhans cell histiocytosis.** J Pediatr. 2008. [PMID 18571550](https://pubmed.ncbi.nlm.nih.gov/18571550/) — 치상핵·뇌교 T2 신호와 임상 ND의 시간차(`ND_radiologic`).
97. Héritier S, Barkaoui MA, Miron J, et al. **Incidence and risk factors for clinical neurodegenerative Langerhans cell histiocytosis: a longitudinal cohort study.** Br J Haematol. 2018. [PMID 30421536](https://pubmed.ncbi.nlm.nih.gov/30421536/) — CDI가 ND의 최강 위험인자(`W_CDI_ND`).
98. McClain KL, Picarsic J, Chakraborty R, et al. **CNS Langerhans cell histiocytosis: common hematopoietic origin for LCH-associated neurodegeneration and mass lesions.** Cancer. 2018. [PMID 29624648](https://pubmed.ncbi.nlm.nih.gov/29624648/) — `BRAF_mono_CNS` 노드.
99. Wilk CM, Cathomas F, Török O, et al. **Circulating senescent myeloid cells infiltrate the brain and cause neurodegeneration in histiocytic disorders.** Immunity. 2023. [PMID 38091952](https://pubmed.ncbi.nlm.nih.gov/38091952/) — **`LCNS` → `NEUR` 축의 기계론적 근거.**
100. Vicario R, Fragkogianni S, Pokrovskii M, et al. **Role of clonal inflammatory microglia in histiocytosis-associated neurodegeneration.** Neuron. 2025. [PMID 40081365](https://pubmed.ncbi.nlm.nih.gov/40081365/)
101. Yeh EA, Greenberg J, Abla O, et al. **Evaluation and treatment of Langerhans cell histiocytosis patients with central nervous system abnormalities: current views and new vistas.** Pediatr Blood Cancer. 2018. [PMID 28944988](https://pubmed.ncbi.nlm.nih.gov/28944988/)
102. Imashuku S, Ishida S, Koike K, et al. **Follow-up of pediatric patients treated by IVIG for Langerhans cell histiocytosis (LCH)-related neurodegenerative CNS disease.** Int J Hematol. 2015. [PMID 25491495](https://pubmed.ncbi.nlm.nih.gov/25491495/) — `IVIG` 노드(안정화, 회복은 드묾).
103. Imashuku S. **High dose immunoglobulin (IVIG) may reduce the incidence of Langerhans cell histiocytosis (LCH)-associated neurodegenerative disease.** CNS Neurol Disord Drug Targets. 2009. [PMID 19702569](https://pubmed.ncbi.nlm.nih.gov/19702569/)
104. Baek C, et al. **Quantitative brain MRI analysis in neurodegenerative Langerhans cell histiocytosis.** Eur J Neurol. 2025. [PMID 40522095](https://pubmed.ncbi.nlm.nih.gov/40522095/)

## L. 간 · 폐 침범 (Liver & Lung Involvement)

모델 대응: `LRO`, `BILF`, `LLUNG`, `LUNGC`, `SMOKE`, 시나리오 S10/S10b

105. Braier J, Ciocca M, Latella A, et al. **Cholestasis, sclerosing cholangitis, and liver transplantation in Langerhans cell histiocytosis.** Med Pediatr Oncol. 2002. [PMID 11836717](https://pubmed.ncbi.nlm.nih.gov/11836717/) — `BILF`(비가역 담관 섬유화)의 근거.
106. Carrere X, Barkaoui MA, Charlotte F, et al. **High prevalence of BRAF V600E in patients with cholestasis, sclerosing cholangitis or liver fibrosis in Langerhans cell histiocytosis.** Pediatr Blood Cancer. 2021. [PMID 33991404](https://pubmed.ncbi.nlm.nih.gov/33991404/)
107. Ziogas IA, Wu WK, Matsuoka LK, et al. **Liver transplantation for Langerhans cell histiocytosis: a US population-based analysis and systematic review.** Liver Transpl. 2021. [PMID 33484600](https://pubmed.ncbi.nlm.nih.gov/33484600/)
108. Vassallo R, Ryu JH, Schroeder DR, et al. **Clinical outcomes of pulmonary Langerhans'-cell histiocytosis in adults.** N Engl J Med. 2002. [PMID 11844849](https://pubmed.ncbi.nlm.nih.gov/11844849/)
109. Vassallo R, Ryu JH, Colby TV, et al. **Pulmonary Langerhans'-cell histiocytosis.** N Engl J Med. 2000. [PMID 10877650](https://pubmed.ncbi.nlm.nih.gov/10877650/)
110. Tazi A. **Adult pulmonary Langerhans' cell histiocytosis.** Eur Respir J. 2006. [PMID 16772390](https://pubmed.ncbi.nlm.nih.gov/16772390/) — 흡연이 >90%에서 동반(`SMOKE`, `ASMK`).
111. Tazi A, Marc K, Dominique S, et al. **Serial computed tomography and lung function testing in pulmonary Langerhans' cell histiocytosis.** Eur Respir J. 2012. [PMID 22441752](https://pubmed.ncbi.nlm.nih.gov/22441752/) — 결절은 소실되지만 낭성 변화는 비가역(`LUNGC`).
112. Benattia A, Bugnet E, Walter-Petrich A, et al. **Long-term outcomes of adult pulmonary Langerhans cell histiocytosis: a prospective cohort.** Eur Respir J. 2022. [PMID 34675043](https://pubmed.ncbi.nlm.nih.gov/34675043/)
113. Elia D, Torre O, Cassandro R, et al. **Pulmonary Langerhans cell histiocytosis: a comprehensive analysis of 40 patients and literature review.** Eur J Intern Med. 2015. [PMID 25899682](https://pubmed.ncbi.nlm.nih.gov/25899682/)
114. Ronceray L, Pötschger U, Janka G, et al. **Pulmonary involvement in pediatric-onset multisystem Langerhans cell histiocytosis: effect on course and outcome.** J Pediatr. 2012. [PMID 22284564](https://pubmed.ncbi.nlm.nih.gov/22284564/)
115. Braier J, Latella A, Balancini B, et al. **Outcome in children with pulmonary Langerhans cell histiocytosis.** Pediatr Blood Cancer. 2004. [PMID 15390304](https://pubmed.ncbi.nlm.nih.gov/15390304/)

## M. 약물 PK/PD 파라미터 출처 (Drug PK/PD Sources)

모델 대응: `$PARAM` PK 블록 전체

116. Zhang W, Heinzmann D, Grippo JF. **Clinical pharmacokinetics of vemurafenib.** Clin Pharmacokinet. 2017. [PMID 28255850](https://pubmed.ncbi.nlm.nih.gov/28255850/) — `CL_VEM`, `V2_VEM`, t½ ≈ 50–60 h, 단백결합 >99%(`FU_VEM`), 자가유도(`KENZ`, `EMAX_IND`).
117. Grippo JF, Zhang W, Heinzmann D, et al. **A phase I, randomized, open-label study of the multiple-dose pharmacokinetics of vemurafenib in patients with BRAF V600E mutation-positive metastatic melanoma.** Cancer Chemother Pharmacol. 2014. [PMID 24178368](https://pubmed.ncbi.nlm.nih.gov/24178368/)
118. Flaherty KT, Puzanov I, Kim KB, et al. **Inhibition of mutated, activated BRAF in metastatic melanoma.** N Engl J Med. 2010. [PMID 20818844](https://pubmed.ncbi.nlm.nih.gov/20818844/)
119. Falchook GS, Long GV, Kurzrock R, et al. **Dabrafenib in patients with melanoma, untreated brain metastases, and other solid tumours: a phase 1 dose-escalation trial.** Lancet. 2012. [PMID 22608338](https://pubmed.ncbi.nlm.nih.gov/22608338/) — `FCNS_DAB`(CNS 투과).
120. Ouellet D, Gibiansky E, Leonowens C, et al. **Population pharmacokinetics of dabrafenib, a BRAF inhibitor: effect of dose, time, covariates, and relationship with its metabolites.** J Clin Pharmacol. 2014. [PMID 24408395](https://pubmed.ncbi.nlm.nih.gov/24408395/) — `CL_DAB`, `FM_DAB`, 하이드록시 대사체 처리.
121. Balakirouchenane D, Guégan S, Csajka C, et al. **Population pharmacokinetics/pharmacodynamics of dabrafenib plus trametinib in patients with BRAF-mutated metastatic melanoma.** Cancers (Basel). 2020. [PMID 32283865](https://pubmed.ncbi.nlm.nih.gov/32283865/)
122. Ouellet D, Kassir N, Chiu J, et al. **Population pharmacokinetics and exposure-response of trametinib, a MEK inhibitor, in patients with BRAF V600 mutation-positive melanoma.** Cancer Chemother Pharmacol. 2016. [PMID 26940938](https://pubmed.ncbi.nlm.nih.gov/26940938/) — `CL_TRA`, 긴 t½(≈4–5일)와 정상상태 지연.
123. Infante JR, Fecher LA, Falchook GS, et al. **Safety, pharmacokinetic, pharmacodynamic, and efficacy data for the oral MEK inhibitor trametinib: a phase 1 dose-escalation trial.** Lancet Oncol. 2012. [PMID 22805291](https://pubmed.ncbi.nlm.nih.gov/22805291/)
124. Liliemark J. **The clinical pharmacokinetics of cladribine.** Clin Pharmacokinet. 1997. [PMID 9068927](https://pubmed.ncbi.nlm.nih.gov/9068927/) — `CL_CLAD`, `V_CLAD`, CSF 투과 ≈25%(`FCNS_CLAD`).
125. Albertioni F, Lindemalm S, Reichelova V, et al. **Pharmacokinetics of cladribine in plasma and its 5'-monophosphate and 5'-triphosphate in leukemic cells of patients with chronic lymphocytic leukemia.** Clin Cancer Res. 1998. [PMID 9533533](https://pubmed.ncbi.nlm.nih.gov/9533533/) — 세포내 Cd-ATP 축적/소실(`KIN_CLAD`, `KOUT_CLATP`).
126. Lindemalm S, Liliemark J, Juliusson G, et al. **Application of population pharmacokinetics to cladribine.** BMC Pharmacol. 2005. [PMID 15757511](https://pubmed.ncbi.nlm.nih.gov/15757511/)
127. Slevin ML, Piall EM, Aherne GW, et al. **Effect of dose and schedule on pharmacokinetics of high-dose cytosine arabinoside in plasma and cerebrospinal fluid.** J Clin Oncol. 1983. [PMID 6583325](https://pubmed.ncbi.nlm.nih.gov/6583325/) — `CL_ARAC`, `FCNS_ARAC`.
128. Heinemann V, Hertel LW, Grindey GB, Plunkett W. **Comparison of the cellular pharmacokinetics and toxicity of 2',2'-difluorodeoxycytidine and 1-beta-D-arabinofuranosylcytosine.** Cancer Res. 1988. [PMID 3383195](https://pubmed.ncbi.nlm.nih.gov/3383195/) — ara-CTP 세포내 반감기(`KOUT_ARATP`).
129. Lamba JK. **Genetic factors influencing cytarabine therapy.** Pharmacogenomics. 2009. [PMID 19842938](https://pubmed.ncbi.nlm.nih.gov/19842938/) — dCK 활성화 경로.
130. Balis FM, Holcenberg JS, Bleyer WA. **Clinical pharmacokinetics of commonly used anticancer drugs.** Clin Pharmacokinet. 1983. [PMID 6189661](https://pubmed.ncbi.nlm.nih.gov/6189661/) — 빈블라스틴 다구획 처리(`CL_VBL`, `V2_VBL`).
131. Petersen KB, Jusko WJ, Rasmussen M, Schmiegelow K. **Population pharmacokinetics of prednisolone in children with acute lymphoblastic leukemia.** Cancer Chemother Pharmacol. 2003. [PMID 12698270](https://pubmed.ncbi.nlm.nih.gov/12698270/) — `CL_PRE`, `V_PRE`.
132. Möhlmann JE, et al. **Population pharmacokinetics of total and protein-unbound prednisolone in children with immune-mediated disease.** Clin Pharmacokinet. 2026. [PMID 41379279](https://pubmed.ncbi.nlm.nih.gov/41379279/)
133. Friberg LE, Henningsson A, Maas H, et al. **Model of chemotherapy-induced myelosuppression with parameter consistency across drugs.** J Clin Oncol. 2002. [PMID 12488418](https://pubmed.ncbi.nlm.nih.gov/12488418/) — **`PROL`/`TR1`/`TR2`/`ANC` 구조와 `KTR_F`, `GAM`, `KCIRC`의 출처.**

## N. 중추성 요붕증 진단 (Diagnosis of Central Diabetes Insipidus)

모델 대응: `AVPN` → `CDI` 판정, `Water_dep`/copeptin 노드

134. Fenske W, Refardt J, Chifu I, et al. **A copeptin-based approach in the diagnosis of diabetes insipidus.** N Engl J Med. 2018. [PMID 30067922](https://pubmed.ncbi.nlm.nih.gov/30067922/)
135. Christ-Crain M, Bichet DG, Fenske WK, et al. **Diabetes insipidus.** Nat Rev Dis Primers. 2019. [PMID 31395885](https://pubmed.ncbi.nlm.nih.gov/31395885/)
136. Gippert S, et al. **Arginine-stimulated copeptin-based diagnosis of central diabetes insipidus in children and adolescents.** Horm Res Paediatr. 2024. [PMID 37607514](https://pubmed.ncbi.nlm.nih.gov/37607514/)

---

## 문헌이 모델의 어느 부분을 결정했는가 (Traceability Summary)

| 모델 요소 | 결정한 문헌 |
|---|---|
| 기원세포 파티션 θ (구조적 선택 ①) | 14, 15, 16, 18, 19, 20 |
| pERK 단일 수렴 노드 | 3, 4, 5, 8, 73 |
| >80% pERK 억제 문턱 (`EC50C`, `HC`) | 24 |
| 낮은 증식분획 `PFMAX=0.09` | 27, 25 |
| SASP·세네센스 분비체 | 25, 28, 31, 32 |
| 병변 내 포획(CCR6/CCR7) | 26, 30 |
| MAPKi = 세포정지성, 저장고 잔존 (구조적 선택 ②) | 70, 74, 75, 76, 79 |
| 중단 후 재발률 | 70, 76, 79 |
| 클라드리빈이 비분열 세포에도 작용 | 63, 64, 124, 125 |
| RANKL/OPG → 골 흡수 | 40, 41, 42, 43, 37, 38 |
| 후유증 = 활성 질환의 시간 적분 (구조적 선택 ③) | 88, 90, 91, 94, 95, 96, 97 |
| CDI → ND 위험 증폭 `W_CDI_ND` | 97, 99, 100 |
| 담관 섬유화 비가역성 | 105, 106, 107 |
| 폐 낭성 변화 비가역성 | 110, 111, 112 |
| 골수억제(Friberg) 파라미터 | 133, 63 |
| 표적치료 피부독성 | 81, 82 |
| DAS 척도 | 50, 51 |

**총 136편** (모두 PMID 실재 검증 완료).
