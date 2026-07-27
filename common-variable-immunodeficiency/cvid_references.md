# 보통 가변 면역결핍증 (CVID) — 참고문헌 (References)

**Common Variable Immunodeficiency · QSP Model Reference List**

이 목록의 **모든 PMID는 NCBI E-utilities로 조회하여 저자·연도·저널·제목을
확인한 것**입니다. 기억에 의존해 적은 PMID는 사용하지 않았습니다 (초안 작성
과정에서 기억으로 적은 5개 중 4개가 전혀 다른 논문을 가리켰고, 모두 조회 결과로
교체했습니다).

---

## 이 모델이 정량적으로 의존하는 핵심 앵커 (The Load-bearing Anchors)

모델의 파라미터가 **직접 고정되어 있는** 문헌은 다음 6개입니다. 나머지 문헌은
구조(어떤 노드가 어떤 노드에 연결되는지)를 정당화하며, 아래 6개는 숫자를
정합니다.

| 앵커 | 문헌 | 관찰값 | 모델값 |
|---|---|---|---|
| 폐렴 노출-반응 기울기 | Orange 2010 [PMID 20675197](https://pubmed.ncbi.nlm.nih.gov/20675197/) | IgG trough 100 mg/dL당 −27% | −29.0% |
| 폐렴 발생률 앵커점 | Orange 2010 | trough 500 mg/dL에서 0.113/년 | 0.113/년 |
| SC:IV 용량 보정계수 | EMA / Berger·Jolles [PMID 25384609](https://pubmed.ncbi.nlm.nih.gov/25384609/) | 1.37 (AUC 매칭) | 1.341 |
| 비감염 합병증 사망 RR | Resnick 2012 [PMID 22180439](https://pubmed.ncbi.nlm.nih.gov/22180439/) | 약 11배 | 10.5배 |
| smB− 전환기억 B세포 | Wehr 2008 EUROclass [PMID 17898316](https://pubmed.ncbi.nlm.nih.gov/17898316/) | B세포의 <2% | 1.67% |
| 기관지확장증 진행과 trough | Lucas 2010 [PMID 20471071](https://pubmed.ncbi.nlm.nih.gov/20471071/) | trough ≥600–800에서 진행 억제 | 0.049 → 0.002 Reiff점/년 |

**과소예측을 그대로 보고하는 항목:** 레니올리십의 12주 지표 병변 반응은
관찰값 −39% (Rao 2023 Blood [PMID 36399712](https://pubmed.ncbi.nlm.nih.gov/36399712/))
대비 모델 −21%입니다. 임의의 림프절 축소 항을 제거하고 PI3Kδ 억제만 남긴
결과이므로 맞추지 않고 그대로 둡니다.

---

## 1. 정의·진단기준·역학·레지스트리 (Definition, Diagnostic Criteria, Epidemiology & Registries)

1. **Cunningham-Rundles C, Bodian C** Common variable immunodeficiency: clinical and immunological features of 248 patients. *Clin Immunol 1999 92:34-48*. [PMID 10413651](https://pubmed.ncbi.nlm.nih.gov/10413651/)
2. **Quinti I et al.** Long-term follow-up and outcome of a large cohort of patients with common variable immunodeficiency. *J Clin Immunol 2007 27:308-16*. [PMID 17510807](https://pubmed.ncbi.nlm.nih.gov/17510807/)
3. **Chapel H et al.** Common variable immunodeficiency disorders: division into distinct clinical phenotypes. *Blood 2008 112:277-86*. [PMID 18319398](https://pubmed.ncbi.nlm.nih.gov/18319398/)
4. **Resnick ES et al.** Morbidity and mortality in common variable immune deficiency over 4 decades. *Blood 2012 119:1650-7*. [PMID 22180439](https://pubmed.ncbi.nlm.nih.gov/22180439/)
5. **Gathmann B et al.** Clinical picture and treatment of 2212 patients with common variable immunodeficiency. *J Allergy Clin Immunol 2014 134:116-26*. [PMID 24582312](https://pubmed.ncbi.nlm.nih.gov/24582312/)
6. **Bonilla FA et al.** International Consensus Document (ICON): Common Variable Immunodeficiency Disorders. *J Allergy Clin Immunol Pract 2016 4:38-59*. [PMID 26563668](https://pubmed.ncbi.nlm.nih.gov/26563668/)
7. **Odnoletkova I et al.** The burden of common variable immunodeficiency disorders: a retrospective analysis of the European Society for Immunodeficiency (ESID) registry data. *Orphanet J Rare Dis 2018 13:201*. [PMID 30419968](https://pubmed.ncbi.nlm.nih.gov/30419968/)
8. **Seidel MG et al.** The European Society for Immunodeficiencies (ESID) Registry Working Definitions for the Clinical Diagnosis of Inborn Errors of Immunity. *J Allergy Clin Immunol Pract 2019 7:1763-1770*. [PMID 30776527](https://pubmed.ncbi.nlm.nih.gov/30776527/)
9. **Ameratunga R, Woon ST** Perspective: Evolving Concepts in the Diagnosis and Understanding of Common Variable Immunodeficiency Disorders (CVID). *Clin Rev Allergy Immunol 2020 59:109-121*. [PMID 31720921](https://pubmed.ncbi.nlm.nih.gov/31720921/)
10. **Ho HE, Cunningham-Rundles C** Non-infectious Complications of Common Variable Immunodeficiency: Updated Clinical Spectrum, Sequelae, and Insights to Pathogenesis. *Front Immunol 2020 11:149*. [PMID 32117289](https://pubmed.ncbi.nlm.nih.gov/32117289/)
11. **Pulvirenti F et al.** Health-Related Quality of Life in Common Variable Immunodeficiency Italian Patients Switched to Remote Assistance During the COVID-19 Pandemic. *J Allergy Clin Immunol Pract 2020 8:1894-1899.e2*. [PMID 32278865](https://pubmed.ncbi.nlm.nih.gov/32278865/)
12. **Cunningham-Rundles C et al.** Genetics and clinical phenotypes in common variable immunodeficiency. *Front Genet 2023 14:1272912*. [PMID 38274105](https://pubmed.ncbi.nlm.nih.gov/38274105/)

## 2. B세포 발달·계열전환 차단·EUROclass 분류 (B-cell Development, the Class-switch Block & Classification)

13. **Warnatz K et al.** Expansion of CD19(hi)CD21(lo/neg) B cells in common variable immunodeficiency (CVID) patients with autoimmune cytopenia. *Immunobiology 2002 206:502-13*. [PMID 12607725](https://pubmed.ncbi.nlm.nih.gov/12607725/)
14. **Piqueras B et al.** Common variable immunodeficiency patient classification based on impaired B cell memory differentiation correlates with clinical aspects. *J Clin Immunol 2003 23:385-400*. [PMID 14601647](https://pubmed.ncbi.nlm.nih.gov/14601647/)
15. **Wehr C et al.** The EUROclass trial: defining subgroups in common variable immunodeficiency. *Blood 2008 111:77-85*. [PMID 17898316](https://pubmed.ncbi.nlm.nih.gov/17898316/)
16. **Rakhmanov M et al.** Circulating CD21low B cells in common variable immunodeficiency resemble tissue homing, innate-like B cells. *Proc Natl Acad Sci U S A 2009 106:13451-6*. [PMID 19666505](https://pubmed.ncbi.nlm.nih.gov/19666505/)
17. **Isnardi I et al.** Complement receptor 2/CD21- human naive B cells contain mostly autoreactive unresponsive clones. *Blood 2010 115:5026-36*. [PMID 20231422](https://pubmed.ncbi.nlm.nih.gov/20231422/)
18. **Romberg N et al.** Patients with common variable immunodeficiency with autoimmune cytopenias exhibit hyperplastic yet inefficient germinal center responses. *J Allergy Clin Immunol 2019 143:258-265*. [PMID 29935219](https://pubmed.ncbi.nlm.nih.gov/29935219/)
19. **Victora GD, Nussenzweig MC** Germinal Centers. *Annu Rev Immunol 2022 40:413-442*. [PMID 35113731](https://pubmed.ncbi.nlm.nih.gov/35113731/)

## 3. T세포 도움·조절 실패·CTLA-4 (T-cell Help, Regulatory Failure & CTLA-4)

20. **Giovannetti A et al.** Unravelling the complexity of T cell abnormalities in common variable immunodeficiency. *J Immunol 2007 178:3932-43*. [PMID 17339494](https://pubmed.ncbi.nlm.nih.gov/17339494/)
21. **Malphettes M et al.** Late-onset combined immune deficiency: a subset of common variable immunodeficiency with severe T cell defect. *Clin Infect Dis 2009 49:1329-38*. [PMID 19807277](https://pubmed.ncbi.nlm.nih.gov/19807277/)
22. **Qureshi OS et al.** Trans-endocytosis of CD80 and CD86: a molecular basis for the cell-extrinsic function of CTLA-4. *Science 2011 332:600-3*. [PMID 21474713](https://pubmed.ncbi.nlm.nih.gov/21474713/)
23. **Walker LS, Sansom DM** The emerging role of CTLA4 as a cell-extrinsic regulator of T cell responses. *Nat Rev Immunol 2011 11:852-63*. [PMID 22116087](https://pubmed.ncbi.nlm.nih.gov/22116087/)
24. **Unger S et al.** The T(H)1 phenotype of follicular helper T cells indicates an IFN-γ-associated immune dysregulation in patients with CD21low common variable immunodeficiency. *J Allergy Clin Immunol 2018 141:730-740*. [PMID 28554560](https://pubmed.ncbi.nlm.nih.gov/28554560/)
25. **Crotty S** T Follicular Helper Cell Biology: A Decade of Discovery and Diseases. *Immunity 2019 50:1132-1148*. [PMID 31117010](https://pubmed.ncbi.nlm.nih.gov/31117010/)
26. **Fernando SL et al.** The Immune Dysregulation of Common Variable Immunodeficiency Disorders. *Immunol Lett 2021 230:21-26*. [PMID 33333111](https://pubmed.ncbi.nlm.nih.gov/33333111/)

## 4. 유전적 기반 (Genetic Substrate — TACI, BAFF-R, CTLA4, LRBA, NFKB1/2, PIK3CD)

27. **Grimbacher B et al.** Homozygous loss of ICOS is associated with adult-onset common variable immunodeficiency. *Nat Immunol 2003 4:261-8*. [PMID 12577056](https://pubmed.ncbi.nlm.nih.gov/12577056/)
28. **Castigli E et al.** TACI is mutant in common variable immunodeficiency and IgA deficiency. *Nat Genet 2005 37:829-34*. [PMID 16007086](https://pubmed.ncbi.nlm.nih.gov/16007086/)
29. **Warnatz K et al.** B-cell activating factor receptor deficiency is associated with an adult-onset antibody deficiency syndrome in humans. *Proc Natl Acad Sci U S A 2009 106:13945-50*. [PMID 19666484](https://pubmed.ncbi.nlm.nih.gov/19666484/)
30. **Orange JS et al.** Genome-wide association identifies diverse causes of common variable immunodeficiency. *J Allergy Clin Immunol 2011 127:1360-7.e6*. [PMID 21497890](https://pubmed.ncbi.nlm.nih.gov/21497890/)
31. **Lopez-Herrera G et al.** Deleterious mutations in LRBA are associated with a syndrome of immune deficiency and autoimmunity. *Am J Hum Genet 2012 90:986-1001*. [PMID 22608502](https://pubmed.ncbi.nlm.nih.gov/22608502/)
32. **Angulo I et al.** Phosphoinositide 3-kinase δ gene mutation predisposes to respiratory infection and airway damage. *Science 2013 342:866-71*. [PMID 24136356](https://pubmed.ncbi.nlm.nih.gov/24136356/)
33. **Chen K et al.** Germline mutations in NFKB2 implicate the noncanonical NF-κB pathway in the pathogenesis of common variable immunodeficiency. *Am J Hum Genet 2013 93:812-24*. [PMID 24140114](https://pubmed.ncbi.nlm.nih.gov/24140114/)
34. **Kuehn HS et al.** Immune dysregulation in human subjects with heterozygous germline mutations in CTLA4. *Science 2014 345:1623-1627*. [PMID 25213377](https://pubmed.ncbi.nlm.nih.gov/25213377/)
35. **Lucas CL et al.** Dominant-activating germline mutations in the gene encoding the PI(3)K catalytic subunit p110δ result in T cell senescence and human immunodeficiency. *Nat Immunol 2014 15:88-97*. [PMID 24165795](https://pubmed.ncbi.nlm.nih.gov/24165795/)
36. **Schubert D et al.** Autosomal dominant immune dysregulation syndrome in humans with CTLA4 mutations. *Nat Med 2014 20:1410-1416*. [PMID 25329329](https://pubmed.ncbi.nlm.nih.gov/25329329/)
37. **van Zelm MC et al.** Human CD19 and CD40L deficiencies impair antibody selection and differentially affect somatic hypermutation. *J Allergy Clin Immunol 2014 134:135-44*. [PMID 24418477](https://pubmed.ncbi.nlm.nih.gov/24418477/)
38. **Fliegauf M et al.** Haploinsufficiency of the NF-κB1 Subunit p50 in Common Variable Immunodeficiency. *Am J Hum Genet 2015 97:389-403*. [PMID 26279205](https://pubmed.ncbi.nlm.nih.gov/26279205/)
39. **Li J et al.** Association of CLEC16A with human common variable immunodeficiency disorder and role in murine B cells. *Nat Commun 2015 6:6804*. [PMID 25891430](https://pubmed.ncbi.nlm.nih.gov/25891430/)
40. **Lo B et al.** AUTOIMMUNE DISEASE. Patients with LRBA deficiency show CTLA4 loss and immune dysregulation responsive to abatacept therapy. *Science 2015 349:436-40*. [PMID 26206937](https://pubmed.ncbi.nlm.nih.gov/26206937/)
41. **Bogaert DJ et al.** Genes associated with common variable immunodeficiency: one diagnosis to rule them all?. *J Med Genet 2016 53:575-90*. [PMID 27250108](https://pubmed.ncbi.nlm.nih.gov/27250108/)
42. **Coulter TI et al.** Clinical spectrum and features of activated phosphoinositide 3-kinase δ syndrome: A large patient cohort study. *J Allergy Clin Immunol 2017 139:597-606.e4*. [PMID 27555459](https://pubmed.ncbi.nlm.nih.gov/27555459/)
43. **Boutboul D et al.** Dominant-negative IKZF1 mutations cause a T, B, and myeloid cell combined immunodeficiency. *J Clin Invest 2018 128:3071-3087*. [PMID 29889099](https://pubmed.ncbi.nlm.nih.gov/29889099/)
44. **Maccari ME et al.** Disease Evolution and Response to Rapamycin in Activated Phosphoinositide 3-Kinase δ Syndrome: The European Society for Immunodeficiencies-Activated Phosphoinositide 3-Kinase δ Syndrome Registry. *Front Immunol 2018 9:543*. [PMID 29599784](https://pubmed.ncbi.nlm.nih.gov/29599784/)
45. **Schwab C et al.** Phenotype, penetrance, and treatment of 133 cytotoxic T-lymphocyte antigen 4-insufficient subjects. *J Allergy Clin Immunol 2018 142:1932-1946*. [PMID 29729943](https://pubmed.ncbi.nlm.nih.gov/29729943/)
46. **Salzer U, Grimbacher B** TACI deficiency - a complex system out of balance. *Curr Opin Immunol 2021 71:81-88*. [PMID 34247095](https://pubmed.ncbi.nlm.nih.gov/34247095/)
47. **Rao VK et al.** A randomised, placebo-controlled, phase III trial of leniolisib in activated phosphoinositide 3-kinase delta (PI3Kδ) syndrome (APDS): Adolescent and adult subgroup analysis. *Clin Immunol 2025 270:110400*. [PMID 39561927](https://pubmed.ncbi.nlm.nih.gov/39561927/)

## 5. BAFF / APRIL 축과 B세포 sink (The BAFF-APRIL Axis and the B-cell Sink)

48. **Vallerskog T et al.** Differential effects on BAFF and APRIL levels in rituximab-treated patients with systemic lupus erythematosus and rheumatoid arthritis. *Arthritis Res Ther 2006 8:R167*. [PMID 17092341](https://pubmed.ncbi.nlm.nih.gov/17092341/)
49. **Knight AK et al.** High serum levels of BAFF, APRIL, and TACI in common variable immunodeficiency. *Clin Immunol 2007 124:182-9*. [PMID 17556024](https://pubmed.ncbi.nlm.nih.gov/17556024/)
50. **Kreuzaler M et al.** Soluble BAFF levels inversely correlate with peripheral B cell numbers and the expression of BAFF receptors. *J Immunol 2012 188:497-503*. [PMID 22124120](https://pubmed.ncbi.nlm.nih.gov/22124120/)
51. **Cambridge G et al.** The effect of B-cell depletion therapy on serological evidence of B-cell and plasmablast activation in patients with rheumatoid arthritis over multiple cycles of rituximab treatment. *J Autoimmun 2014 50:67-76*. [PMID 24365380](https://pubmed.ncbi.nlm.nih.gov/24365380/)

## 6. 면역글로불린 보충 약동학 (Immunoglobulin Replacement Pharmacokinetics)

52. **Koleba T, Ensom MH** Pharmacokinetics of intravenous immunoglobulin: a systematic review. *Pharmacotherapy 2006 26:813-27*. [PMID 16716135](https://pubmed.ncbi.nlm.nih.gov/16716135/)
53. **Roopenian DC, Akilesh S** FcRn: the neonatal Fc receptor comes of age. *Nat Rev Immunol 2007 7:715-25*. [PMID 17703228](https://pubmed.ncbi.nlm.nih.gov/17703228/)
54. **Borte M et al.** Efficacy and safety of subcutaneous vivaglobin® replacement therapy in previously untreated patients with primary immunodeficiency: a prospective, multicenter study. *J Clin Immunol 2011 31:952-61*. [PMID 21932110](https://pubmed.ncbi.nlm.nih.gov/21932110/)
55. **Shapiro RS** Why I use subcutaneous immunoglobulin (SCIG). *J Clin Immunol 2013 33 Suppl 2:S95-8*. [PMID 23264027](https://pubmed.ncbi.nlm.nih.gov/23264027/)
56. **Shapiro RS** Subcutaneous immunoglobulin therapy given by subcutaneous rapid push vs infusion pump: a retrospective analysis. *Ann Allergy Asthma Immunol 2013 111:51-5*. [PMID 23806460](https://pubmed.ncbi.nlm.nih.gov/23806460/)
57. **Jolles S et al.** Current treatment options with immunoglobulin G for the individualization of care in patients with primary immunodeficiency disease. *Clin Exp Immunol 2015 179:146-60*. [PMID 25384609](https://pubmed.ncbi.nlm.nih.gov/25384609/)
58. **Wasserman RL** Recombinant human hyaluronidase-facilitated subcutaneous immunoglobulin infusion in primary immunodeficiency diseases. *Immunotherapy 2017 9:1035-1050*. [PMID 28871852](https://pubmed.ncbi.nlm.nih.gov/28871852/)
59. **Jolles S et al.** Long-Term Efficacy and Safety of Hizentra® in Patients with Primary Immunodeficiency in Japan, Europe, and the United States: a Review of 7 Phase 3 Trials. *J Clin Immunol 2018 38:864-875*. [PMID 30415311](https://pubmed.ncbi.nlm.nih.gov/30415311/)
60. **Ward ES, Ober RJ** Targeting FcRn to Generate Antibody-Based Therapeutics. *Trends Pharmacol Sci 2018 39:892-904*. [PMID 30143244](https://pubmed.ncbi.nlm.nih.gov/30143244/)

## 7. 노출-반응: trough IgG와 감염 (Exposure-Response — Trough IgG and Infection)

61. **Roifman CM et al.** High-dose versus low-dose intravenous immunoglobulin in hypogammaglobulinaemia and chronic lung disease. *Lancet 1987 1:1075-7*. [PMID 2883406](https://pubmed.ncbi.nlm.nih.gov/2883406/)
62. **Eijkhout HW et al.** The effect of two different dosages of intravenous immunoglobulin on the incidence of recurrent infections in patients with primary hypogammaglobulinemia. A randomized, double-blind, multicenter crossover trial. *Ann Intern Med 2001 135:165-74*. [PMID 11487483](https://pubmed.ncbi.nlm.nih.gov/11487483/)
63. **Busse PJ et al.** Efficacy of intravenous immunoglobulin in the prevention of pneumonia in patients with common variable immunodeficiency. *J Allergy Clin Immunol 2002 109:1001-4*. [PMID 12063531](https://pubmed.ncbi.nlm.nih.gov/12063531/)
64. **Lucas M et al.** Infection outcomes in patients with common variable immunodeficiency disorders: relationship to immunoglobulin therapy over 22 years. *J Allergy Clin Immunol 2010 125:1354-1360.e4*. [PMID 20471071](https://pubmed.ncbi.nlm.nih.gov/20471071/)
65. **Orange JS et al.** Impact of trough IgG on pneumonia incidence in primary immunodeficiency: A meta-analysis of clinical studies. *Clin Immunol 2010 137:21-30*. [PMID 20675197](https://pubmed.ncbi.nlm.nih.gov/20675197/)
66. **Bonagura VR** Illustrative cases on individualizing immunoglobulin therapy in primary immunodeficiency disease. *Ann Allergy Asthma Immunol 2013 111:S10-3*. [PMID 24267400](https://pubmed.ncbi.nlm.nih.gov/24267400/)
67. **Perez EE et al.** Update on the use of immunoglobulin in human disease: A review of evidence. *J Allergy Clin Immunol 2017 139:S1-S46*. [PMID 28041678](https://pubmed.ncbi.nlm.nih.gov/28041678/)

## 8. 백신 반응·특이항체·옵소닌 기능 (Vaccine Response, Specific Antibody & Opsonic Function)

68. **Goldacker S et al.** Active vaccination in patients with common variable immunodeficiency (CVID). *Clin Immunol 2007 124:294-303*. [PMID 17602874](https://pubmed.ncbi.nlm.nih.gov/17602874/)
69. **Orange JS et al.** Use and interpretation of diagnostic vaccination in primary immunodeficiency: a working group report of the Basic and Clinical Immunology Interest Section of the American Academy of Allergy, Asthma & Immunology. *J Allergy Clin Immunol 2012 130:S1-24*. [PMID 22935624](https://pubmed.ncbi.nlm.nih.gov/22935624/)

## 9. 점막 면역·IgA 공백·미생물 전위 (Mucosal Immunity, the IgA Gap & Microbial Translocation)

70. **Mantis NJ et al.** Secretory IgA's complex roles in immunity and mucosal homeostasis in the gut. *Mucosal Immunol 2011 4:603-11*. [PMID 21975936](https://pubmed.ncbi.nlm.nih.gov/21975936/)
71. **Shulzhenko N et al.** Crosstalk between B lymphocytes, microbiota and the intestinal epithelium governs immunity versus metabolism in the gut. *Nat Med 2011 17:1585-93*. [PMID 22101768](https://pubmed.ncbi.nlm.nih.gov/22101768/)
72. **Brandtzaeg P** Secretory IgA: Designed for Anti-Microbial Defense. *Front Immunol 2013 4:222*. [PMID 23964273](https://pubmed.ncbi.nlm.nih.gov/23964273/)
73. **Perreau M et al.** Exhaustion of bacteria-specific CD4 T cells and microbial translocation in common variable immunodeficiency disorders. *J Exp Med 2014 211:2033-45*. [PMID 25225461](https://pubmed.ncbi.nlm.nih.gov/25225461/)
74. **Jørgensen SF et al.** Altered gut microbiota profile in common variable immunodeficiency associates with levels of lipopolysaccharide and markers of systemic immune activation. *Mucosal Immunol 2016 9:1455-1465*. [PMID 26982597](https://pubmed.ncbi.nlm.nih.gov/26982597/)
75. **Berbers RM et al.** Gut microbial dysbiosis, IgA, and Enterococcus in common variable immunodeficiency with immune dysregulation. *Microbiome 2025 13:12*. [PMID 39819634](https://pubmed.ncbi.nlm.nih.gov/39819634/)

## 10. 폐질환: 기관지확장증과 구조적 손상 (Lung Disease — Bronchiectasis & Structural Damage)

76. **Bhalla M et al.** Cystic fibrosis: scoring system with thin-section CT. *Radiology 1991 179:783-8*. [PMID 2027992](https://pubmed.ncbi.nlm.nih.gov/2027992/)
77. **Reiff DB et al.** CT findings in bronchiectasis: limited value in distinguishing between idiopathic and specific types. *AJR Am J Roentgenol 1995 165:261-7*. [PMID 7618537](https://pubmed.ncbi.nlm.nih.gov/7618537/)
78. **Kainulainen L et al.** Pulmonary abnormalities in patients with primary hypogammaglobulinemia. *J Allergy Clin Immunol 1999 104:1031-6*. [PMID 10550749](https://pubmed.ncbi.nlm.nih.gov/10550749/)
79. **Thickett KM et al.** Common variable immune deficiency: respiratory manifestations, pulmonary function and high-resolution CT scan findings. *QJM 2002 95:655-62*. [PMID 12324637](https://pubmed.ncbi.nlm.nih.gov/12324637/)
80. **Hampson FA et al.** Respiratory disease in common variable immunodeficiency and other primary immunodeficiency disorders. *Clin Radiol 2012 67:587-95*. [PMID 22226567](https://pubmed.ncbi.nlm.nih.gov/22226567/)
81. **Altenburg J et al.** Effect of azithromycin maintenance treatment on infectious exacerbations among patients with non-cystic fibrosis bronchiectasis: the BAT randomized controlled trial. *JAMA 2013 309:1251-9*. [PMID 23532241](https://pubmed.ncbi.nlm.nih.gov/23532241/)
82. **Chalmers JD et al.** The bronchiectasis severity index. An international derivation and validation study. *Am J Respir Crit Care Med 2014 189:576-85*. [PMID 24328736](https://pubmed.ncbi.nlm.nih.gov/24328736/)

## 11. GLILD·육아종·간질성 폐질환 (GLILD, Granuloma & Interstitial Lung Disease)

83. **Bates CA et al.** Granulomatous-lymphocytic lung disease shortens survival in common variable immunodeficiency. *J Allergy Clin Immunol 2004 114:415-21*. [PMID 15316526](https://pubmed.ncbi.nlm.nih.gov/15316526/)
84. **Chase NM et al.** Use of combination chemotherapy for treatment of granulomatous and lymphocytic interstitial lung disease (GLILD) in patients with common variable immunodeficiency (CVID). *J Clin Immunol 2013 33:30-9*. [PMID 22930256](https://pubmed.ncbi.nlm.nih.gov/22930256/)
85. **Gregersen S et al.** Lung disease, T-cells and inflammation in common variable immunodeficiency disorders. *Scand J Clin Lab Invest 2013 73:514-22*. [PMID 23957371](https://pubmed.ncbi.nlm.nih.gov/23957371/)
86. **Maglione PJ et al.** Progression of Common Variable Immunodeficiency Interstitial Lung Disease Accompanies Distinct Pulmonary and Laboratory Findings. *J Allergy Clin Immunol Pract 2015 3:941-50*. [PMID 26372540](https://pubmed.ncbi.nlm.nih.gov/26372540/)
87. **Hurst JR et al.** British Lung Foundation/United Kingdom Primary Immunodeficiency Network Consensus Statement on the Definition, Diagnosis, and Management of Granulomatous-Lymphocytic Interstitial Lung Disease in Common Variable Immunodeficiency Disorders. *J Allergy Clin Immunol Pract 2017 5:938-945*. [PMID 28351785](https://pubmed.ncbi.nlm.nih.gov/28351785/)
88. **Verbsky JW et al.** Rituximab and antimetabolite treatment of granulomatous and lymphocytic interstitial lung disease in common variable immunodeficiency. *J Allergy Clin Immunol 2021 147:704-712.e17*. [PMID 32745555](https://pubmed.ncbi.nlm.nih.gov/32745555/)
89. **Wobma H et al.** Genetic diagnosis of immune dysregulation can lead to targeted therapy for interstitial lung disease: A case series and single center approach. *Pediatr Pulmonol 2022 57:1577-1587*. [PMID 35426264](https://pubmed.ncbi.nlm.nih.gov/35426264/)
90. **Bintalib HM et al.** Clinical Practice Guideline for the diagnosis of Granulomatous-Lymphocytic Interstitial Lung Disease (GLILD) in patients with Common Variable Immunodeficiency Disorders (CVID) - an ERS Clinical Research Collaboration. *Eur Respir J 2026*. [PMID 42350065](https://pubmed.ncbi.nlm.nih.gov/42350065/)

## 12. 자가면역 세포감소증 (Autoimmune Cytopenias — ITP, AIHA, Evans)

91. **Michel M et al.** Autoimmune thrombocytopenic purpura and common variable immunodeficiency: analysis of 21 cases and review of the literature. *Medicine (Baltimore) 2004 83:254-263*. [PMID 15232313](https://pubmed.ncbi.nlm.nih.gov/15232313/)
92. **Boileau J et al.** Autoimmunity in common variable immunodeficiency: correlation with lymphocyte phenotype in the French DEFI study. *J Autoimmun 2011 36:25-32*. [PMID 21075598](https://pubmed.ncbi.nlm.nih.gov/21075598/)
93. **Gobert D et al.** Efficacy and safety of rituximab in common variable immunodeficiency-associated immune cytopenias: a retrospective multicentre study on 33 patients. *Br J Haematol 2011 155:498-508*. [PMID 21981575](https://pubmed.ncbi.nlm.nih.gov/21981575/)
94. **Podjasek JC, Abraham RS** Autoimmune cytopenias in common variable immunodeficiency. *Front Immunol 2012 3:189*. [PMID 22837758](https://pubmed.ncbi.nlm.nih.gov/22837758/)
95. **Seidel MG** Autoimmune and other cytopenias in primary immunodeficiencies: pathomechanisms, novel differential diagnoses, and treatment. *Blood 2014 124:2337-44*. [PMID 25163701](https://pubmed.ncbi.nlm.nih.gov/25163701/)
96. **Azizi G et al.** Autoimmunity in common variable immunodeficiency: epidemiology, pathophysiology and management. *Expert Rev Clin Immunol 2017 13:101-115*. [PMID 27636680](https://pubmed.ncbi.nlm.nih.gov/27636680/)
97. **Feuille EJ et al.** Autoimmune Cytopenias and Associated Conditions in CVID: a Report From the USIDNET Registry. *J Clin Immunol 2018 38:28-34*. [PMID 29080979](https://pubmed.ncbi.nlm.nih.gov/29080979/)

## 13. 장병증·단백소실·간질환 (Enteropathy, Protein Loss & Liver Disease)

98. **Daniels JA et al.** Gastrointestinal tract pathology in patients with common variable immunodeficiency (CVID): a clinicopathologic study and review. *Am J Surg Pathol 2007 31:1800-12*. [PMID 18043034](https://pubmed.ncbi.nlm.nih.gov/18043034/)
99. **Ward C et al.** Abnormal liver function in common variable immunodeficiency disorders due to nodular regenerative hyperplasia. *Clin Exp Immunol 2008 153:331-7*. [PMID 18647320](https://pubmed.ncbi.nlm.nih.gov/18647320/)
100. **Malamut G et al.** The enteropathy associated with common variable immunodeficiency: the delineated frontiers with celiac disease. *Am J Gastroenterol 2010 105:2262-75*. [PMID 20551941](https://pubmed.ncbi.nlm.nih.gov/20551941/)
101. **Umar SB, DiBaise JK** Protein-losing enteropathy: case illustrations and clinical review. *Am J Gastroenterol 2010 105:43-9; quiz 50*. [PMID 19789526](https://pubmed.ncbi.nlm.nih.gov/19789526/)
102. **Woodward JM et al.** The role of chronic norovirus infection in the enteropathy associated with common variable immunodeficiency. *Am J Gastroenterol 2015 110:320-7*. [PMID 25623655](https://pubmed.ncbi.nlm.nih.gov/25623655/)
103. **Jørgensen SF et al.** A Cross-Sectional Study of the Prevalence of Gastrointestinal Symptoms and Pathology in Patients With Common Variable Immunodeficiency. *Am J Gastroenterol 2016 111:1467-1475*. [PMID 27527747](https://pubmed.ncbi.nlm.nih.gov/27527747/)
104. **Pecoraro A et al.** Heterogeneity of Liver Disease in Common Variable Immunodeficiency Disorders. *Front Immunol 2020 11:338*. [PMID 32184784](https://pubmed.ncbi.nlm.nih.gov/32184784/)
105. **Hercun J et al.** Development of hepatic fibrosis in common variable immunodeficiency-related porto-sinusoidal vascular disorder. *Aliment Pharmacol Ther 2024 60:888-896*. [PMID 39090843](https://pubmed.ncbi.nlm.nih.gov/39090843/)

## 14. 악성종양 (Malignancy)

106. **Dhalla F et al.** Review of gastric cancer risk factors in patients with common variable immunodeficiency disorders, resulting in a proposal for a surveillance programme. *Clin Exp Immunol 2011 165:1-7*. [PMID 21470209](https://pubmed.ncbi.nlm.nih.gov/21470209/)
107. **Tak Manesh A et al.** Epidemiology and pathophysiology of malignancy in common variable immunodeficiency?. *Allergol Immunopathol (Madr) 2017 45:602-615*. [PMID 28411962](https://pubmed.ncbi.nlm.nih.gov/28411962/)
108. **Kiaee F et al.** Malignancy in common variable immunodeficiency: a systematic review and meta-analysis. *Expert Rev Clin Immunol 2019 15:1105-1113*. [PMID 31452405](https://pubmed.ncbi.nlm.nih.gov/31452405/)

## 15. 표적치료·유전자형 지향 치료 (Targeted and Genotype-directed Therapy)

109. **Wehr C et al.** Multicenter experience in hematopoietic stem cell transplantation for serious complications of common variable immunodeficiency. *J Allergy Clin Immunol 2015 135:988-997.e6*. [PMID 25595268](https://pubmed.ncbi.nlm.nih.gov/25595268/)
110. **Rao VK et al.** Effective "activated PI3Kδ syndrome"-targeted therapy with the PI3Kδ inhibitor leniolisib. *Blood 2017 130:2307-2316*. [PMID 28972011](https://pubmed.ncbi.nlm.nih.gov/28972011/)
111. **Taraldsrud E et al.** Defective IL-4 signaling in T cells defines severe common variable immunodeficiency. *J Autoimmun 2017 81:110-119*. [PMID 28476239](https://pubmed.ncbi.nlm.nih.gov/28476239/)
112. **Kiykim A et al.** Abatacept as a Long-Term Targeted Therapy for LRBA Deficiency. *J Allergy Clin Immunol Pract 2019 7:2790-2800.e15*. [PMID 31238161](https://pubmed.ncbi.nlm.nih.gov/31238161/)
113. **Tesch VK et al.** Long-term outcome of LRBA deficiency in 76 patients after various treatment modalities as evaluated by the immune deficiency and dysregulation activity (IDDA) score. *J Allergy Clin Immunol 2020 145:1452-1463*. [PMID 31887391](https://pubmed.ncbi.nlm.nih.gov/31887391/)
114. **Egg D et al.** Therapeutic options for CTLA-4 insufficiency. *J Allergy Clin Immunol 2022 149:736-746*. [PMID 34111452](https://pubmed.ncbi.nlm.nih.gov/34111452/)
115. **Rao VK et al.** A randomized, placebo-controlled phase 3 trial of the PI3Kδ inhibitor leniolisib for activated PI3Kδ syndrome. *Blood 2023 141:971-983*. [PMID 36399712](https://pubmed.ncbi.nlm.nih.gov/36399712/)

## 16. 보충요법의 이상반응 (Adverse Effects of Immunoglobulin Replacement)

116. **Orbach H et al.** Intravenous immunoglobulin: adverse effects and safe administration. *Clin Rev Allergy Immunol 2005 29:173-84*. [PMID 16391392](https://pubmed.ncbi.nlm.nih.gov/16391392/)
117. **Kahwaji J et al.** Acute hemolysis after high-dose intravenous immunoglobulin therapy in highly HLA sensitized patients. *Clin J Am Soc Nephrol 2009 4:1993-7*. [PMID 19833910](https://pubmed.ncbi.nlm.nih.gov/19833910/)
118. **Stiehm ER** Adverse effects of human immunoglobulin therapy. *Transfus Med Rev 2013 27:171-8*. [PMID 23835249](https://pubmed.ncbi.nlm.nih.gov/23835249/)
119. **Ammann EM et al.** Intravenous immune globulin and thromboembolic adverse events: A systematic review and meta-analysis of RCTs. *Am J Hematol 2016 91:594-605*. [PMID 26973084](https://pubmed.ncbi.nlm.nih.gov/26973084/)

## 17. 기타 기전·리뷰·보조 문헌 (Additional Mechanism, Review & Supporting Literature)

120. **Gregersen H et al.** The impact of M-component type and immunoglobulin concentration on the risk of malignant transformation in patients with monoclonal gammopathy of undetermined significance. *Haematologica 2001 86:1172-9*. [PMID 11694403](https://pubmed.ncbi.nlm.nih.gov/11694403/)
121. **El-Shanawany T, Jolles S** Intravenous immunoglobulin and autoimmune disease. *Ann N Y Acad Sci 2007 1110:507-15*. [PMID 17911466](https://pubmed.ncbi.nlm.nih.gov/17911466/)
122. **Magri G et al.** Innate lymphoid cells integrate stromal and immunological signals to enhance antibody production by splenic marginal zone B cells. *Nat Immunol 2014 15:354-364*. [PMID 24562309](https://pubmed.ncbi.nlm.nih.gov/24562309/)
123. **Rubinstein A et al.** Long-Term Safety of Facilitated Subcutaneous Immunoglobulin 10% Treatment in US Clinical Practice in Patients with Primary Immunodeficiency Diseases: Results from a Post-Authorization Safety Study. *J Clin Immunol 2024 44:181*. [PMID 39158670](https://pubmed.ncbi.nlm.nih.gov/39158670/)
124. **Bintalib HM et al.** Investigating pulmonary and non-infectious complications in common variable immunodeficiency disorders: a UK national multi-centre study. *Front Immunol 2024 15:1451813*. [PMID 39318627](https://pubmed.ncbi.nlm.nih.gov/39318627/)

---

## 문헌 활용 요약 (How the Literature Maps onto the Model)

| 모델 요소 | 근거 섹션 |
|---|---|
| 계열전환 차단 (`FCSR`), EUROclass 표현형 | 2 |
| Tfh 확장·Treg/CTLA-4 실패 (ARM 2 뿌리) | 3, 4 |
| BAFF sink 소실 → 리툭시맙 후 BAFF 상승 | 5 |
| 2구획 IgG PK, FcRn 포화, SC 생체이용률 | 6 |
| 폐렴 노출-반응 Hill 함수, 개인별 생물학적 IgG 수준 | 7 |
| 점막 IgA 공백 → 부비동염의 용량-무관 바닥 | 8, 9 |
| 기관지확장증 래칫, Reiff/Bhalla 점수, 아지트로마이신 | 10 |
| GLILD 구조적 영(null), 리툭시맙+아자티오프린 | 11 |
| 자가항체·혈소판·비장절제 트레이드오프 | 12 |
| 단백소실 장병증 → IgG 청소율 되먹임 | 13 |
| 누적 악성종양 위험 | 14 |
| 아바타셉트(CTLA4/LRBA), 레니올리십(APDS), HSCT | 15 |
| 주입 반응·혈전·용혈 등 치료 대가 | 16 |

## 검색 전략 (Search Strategy)

PubMed E-utilities(`esearch` → `esummary`)로 제목·저자 필드 질의를 던져
PMID를 확정한 뒤, 반환된 메타데이터(저자·연도·저널·제목)를 그대로 인용에
사용했습니다. 반환된 논문이 주제에서 벗어난 경우 인용에서 제외했습니다.
따라서 이 목록의 인용 문자열은 기억이 아니라 조회 결과에서 생성된 것입니다.

## 면책 (Disclaimer)

본 QSP 모델은 **교육 및 연구 목적**입니다. 파라미터는 공개 문헌의 집계값에
맞추어 보정되었으나 개별 환자 데이터로 검증되지 않았으며, 임상 의사결정·처방·
규제 제출에 직접 사용해서는 안 됩니다.
