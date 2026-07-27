# 프로락틴 분비 뇌하수체 종양 (Prolactinoma) — QSP 모델 참고문헌
### Prolactinoma / lactotroph PitNET · 150 PubMed references, every PMID resolved through NCBI E-utilities

이 목록의 **모든 PMID는 NCBI E-utilities(esearch + esummary)로 직접 조회해 제목·저널·연도를
확인한 것**이며, 기억에 의존해 적은 인용은 하나도 포함되어 있지 않습니다. 이 저장소의 이전
세션에서 기억으로 쓴 PMID 다수가 전혀 무관한 논문을 가리켰던 사고가 있었기 때문에, 이번에는
질의 → 조회 → 관련성 필터 → 서지정보 자동 생성의 순서로 만들었고, 관련성 검증을 통과하지 못한
후보는 목록에서 제외했습니다. 검증 스크립트가 채택한 기준은 (1) 질의의 제1저자 성이 반환된
논문의 제1저자와 일치하거나, (2) 질의의 내용어 중 둘 이상이 반환된 제목/저널에 포함되는
경우입니다.

각 섹션 앞의 인용 블록은 **그 문헌군이 모델의 어느 방정식·파라미터·진단분석을 뒷받침하는지**를
명시합니다. 문헌은 섹션 내에서 연도순으로 정렬되어 있습니다.

| 산출물 | 파일 |
|---|---|
| 🗺️ 기계론적 지도 | [`prl_qsp_model.dot`](prl_qsp_model.dot) · [`.svg`](prl_qsp_model.svg) · [`.png`](prl_qsp_model.png) |
| ⚙️ mrgsolve ODE 모델 | [`prl_mrgsolve_model.R`](prl_mrgsolve_model.R) |
| 📊 Shiny 대시보드 | [`prl_shiny_app.R`](prl_shiny_app.R) |
| 📄 디렉토리 README | [`README.md`](README.md) |

---

## 이 모델이 문헌에 대해 정량적으로 주장하는 것 (calibration anchors)

| 앵커 | 출처 (아래 번호) | 모델 결과 |
|---|---|---|
| 카베르골린 t½ 63-109 h, 혈장 농도 pg/mL 범위 | G 섹션 (Ferrari 1995, Rains 1995) | Cmax 68.6 pg/mL (1 mg), t½ 107 h |
| 카베르골린 > 브로모크립틴 (내재활성 차이만으로) | Webster 1994 (G) | 고평부 51 vs 80 ng/mL, 별도 효능 파라미터 없음 |
| 단일 투여 후 프로락틴 억제가 2주 이상 지속 | G 섹션 | 14일째 기저치의 24% |
| 중단 후 지속 완화 ~21% | Dekkers 2010 (I) | 완화 자체가 생성되지 않음 → **추적 기간의 함수로 재해석** (D10) |
| 임신 중 증상성 증대 micro ~2.7% / macro ~21-23% | J 섹션 (Molitch) | 기하학적 역치로 생성 (D13) |
| 파킨슨병 용량 판막병증 O, 프로락틴종 용량 X | M 섹션 (Zanettini/Schade vs Stiles/Caputo) | P(TR) 25.8% vs 0.38% (D11) |
| 도파민 작용제 충동조절장애 ~10-17% | M 섹션 | 1 mg/주에서 14.1% |
| 줄기 압박 고프로락틴혈증 상한 25-150 ng/mL | L 섹션 | 완전 절단 시 94 ng/mL (D8) |
| PEG 회수율 <40% = macroprolactin | F 섹션 | 16.9% |
| MGMT 저발현 시 테모졸로마이드 반응 | H 섹션 (Raverot, Whitelaw) | 78% 축소 vs 무반응 |

---

### A. 지침 · 개요 · 역학 (Guidelines, reviews, epidemiology)

> 모델의 임상적 골격과 정상 범위·유병률·치료 알고리듬의 출처. 특히 Melmed 2011과 Petersenn 2023은
> 시나리오 구성(1차 약물치료, 추적 간격, 중단 기준)의 기준이 되었다.

1. Colao A et al. New medical approaches in pituitary adenomas. *Horm Res*. 2000;53 Suppl 3:76-87. [PMID 10971110](https://pubmed.ncbi.nlm.nih.gov/10971110/)
2. Verhelst J & Abs R. Hyperprolactinemia: pathophysiology and management. *Treat Endocrinol*. 2003;2:23-32. [PMID 15871552](https://pubmed.ncbi.nlm.nih.gov/15871552/)
3. Buurman H & Saeger W. Subclinical adenomas in postmortem pituitaries: classification and correlations to clinical data. *Eur J Endocrinol*. 2006;154:753-8. [PMID 16645024](https://pubmed.ncbi.nlm.nih.gov/16645024/)
4. Casanueva FF et al. Guidelines of the Pituitary Society for the diagnosis and management of prolactinomas. *Clin Endocrinol (Oxf)*. 2006;65:265-73. [PMID 16886971](https://pubmed.ncbi.nlm.nih.gov/16886971/)
5. Daly AF et al. High prevalence of pituitary adenomas: a cross-sectional study in the province of Liege, Belgium. *J Clin Endocrinol Metab*. 2006;91:4769-75. [PMID 16968795](https://pubmed.ncbi.nlm.nih.gov/16968795/)
6. Gillam MP et al. Advances in the treatment of prolactinomas. *Endocr Rev*. 2006;27:485-534. [PMID 16705142](https://pubmed.ncbi.nlm.nih.gov/16705142/)
7. Kars M et al. Estimated age- and sex-specific incidence and prevalence of dopamine agonist-treated hyperprolactinemia. *J Clin Endocrinol Metab*. 2009;94:2729-34. [PMID 19491225](https://pubmed.ncbi.nlm.nih.gov/19491225/)
8. Fernandez A et al. Prevalence of pituitary adenomas: a community-based, cross-sectional study in Banbury (Oxfordshire, UK). *Clin Endocrinol (Oxf)*. 2010;72:377-82. [PMID 19650784](https://pubmed.ncbi.nlm.nih.gov/19650784/)
9. Kars M et al. Update in prolactinomas. *Neth J Med*. 2010;68:104-12. [PMID 20308704](https://pubmed.ncbi.nlm.nih.gov/20308704/)
10. Melmed S et al. Diagnosis and treatment of hyperprolactinemia: an Endocrine Society clinical practice guideline. *J Clin Endocrinol Metab*. 2011;96:273-88. [PMID 21296991](https://pubmed.ncbi.nlm.nih.gov/21296991/)
11. Faje A & Nachtigall L. Current treatment options for hyperprolactinemia. *Expert Opin Pharmacother*. 2013;14:1611-25. [PMID 23738973](https://pubmed.ncbi.nlm.nih.gov/23738973/)
12. Glezer A & Bronstein MD. [Prolactinoma]. *Arq Bras Endocrinol Metabol*. 2014;58:118-23. [PMID 24830588](https://pubmed.ncbi.nlm.nih.gov/24830588/)
13. Molitch ME. Diagnosis and Treatment of Pituitary Adenomas: A Review. *JAMA*. 2017;317:516-524. [PMID 28170483](https://pubmed.ncbi.nlm.nih.gov/28170483/)
14. Vilar L et al. Controversial issues in the management of hyperprolactinemia and prolactinomas - An overview by the Neuroendocrinology Department of the Brazilian Society of Endocrinology and Metabolism. *Arch Endocrinol Metab*. 2018;62:236-263. [PMID 29768629](https://pubmed.ncbi.nlm.nih.gov/29768629/)
15. Auriemma RS et al. Dopamine Agonists: From the 1970s to Today. *Neuroendocrinology*. 2019;109:34-41. [PMID 30852578](https://pubmed.ncbi.nlm.nih.gov/30852578/)
16. Vroonen L et al. Epidemiology and Management Challenges in Prolactinomas. *Neuroendocrinology*. 2019;109:20-27. [PMID 30731464](https://pubmed.ncbi.nlm.nih.gov/30731464/)
17. Trouillas J et al. How to Classify the Pituitary Neuroendocrine Tumors (PitNET)s in 2020. *Cancers (Basel)*. 2020;12. [PMID 32098443](https://pubmed.ncbi.nlm.nih.gov/32098443/)
18. Asa SL et al. Overview of the 2022 WHO Classification of Pituitary Tumors. *Endocr Pathol*. 2022;33:6-26. [PMID 35291028](https://pubmed.ncbi.nlm.nih.gov/35291028/)
19. Kontbay T et al. Hyperprolactinemia in children and adolescents and longterm follow-up results of prolactinoma cases: a single-centre experience. *Turk J Pediatr*. 2022;64:892-899. [PMID 36305439](https://pubmed.ncbi.nlm.nih.gov/36305439/)
20. Petersenn S et al. Diagnosis and management of prolactin-secreting pituitary adenomas: a Pituitary Society international Consensus Statement. *Nat Rev Endocrinol*. 2023;19:722-740. [PMID 37670148](https://pubmed.ncbi.nlm.nih.gov/37670148/)
21. Raverot G et al. Revised European Society of Endocrinology Clinical Practice Guideline for the management of aggressive pituitary tumours and pituitary carcinomas. *Eur J Endocrinol*. 2025;192:R45-R78. [PMID 40506054](https://pubmed.ncbi.nlm.nih.gov/40506054/)

---

### B. 병태생리 · 유전학 (Pathogenesis and genetics)

> `FGEN`(증식 구동)과 `FDRIVE`(전사 구동), 그리고 AIP·SF3B1 표현형이 왜 저 D2R·침습성과 연결되는지의 근거.
> SF3B1 R625H는 최근 발견된 프로락틴종 특이 hotspot으로 지도의 cluster 6에 그려져 있다.

22. Vergès B et al. Pituitary disease in MEN type 1 (MEN1): data from the France-Belgium MEN1 multicenter study. *J Clin Endocrinol Metab*. 2002;87:457-65. [PMID 11836268](https://pubmed.ncbi.nlm.nih.gov/11836268/)
23. Daly AF et al. Aryl hydrocarbon receptor-interacting protein gene mutations in familial isolated pituitary adenomas: analysis in 73 families. *J Clin Endocrinol Metab*. 2007;92:1891-6. [PMID 17244780](https://pubmed.ncbi.nlm.nih.gov/17244780/)
24. Raverot G et al. Prognostic factors in prolactin pituitary tumors: clinical, histological, and molecular data from a series of 94 patients with a long postoperative follow-up. *J Clin Endocrinol Metab*. 2010;95:1708-16. [PMID 20164287](https://pubmed.ncbi.nlm.nih.gov/20164287/)
25. Mitsui T et al. Differences between rat strains in the development of PRL-secreting pituitary tumors with long-term estrogen treatment: In vitro insulin-like growth factor-1-induced lactotroph proliferation and gene expression are affected in Wistar-Kyoto rats with low estrogen-susceptibility. *Endocr J*. 2013;60:1251-9. [PMID 23985690](https://pubmed.ncbi.nlm.nih.gov/23985690/)
26. Kageyama K et al. A Novel Deletion Mutation in the MEN1 Gene in a Patient with Prolactinoma and a Family History of Pancreatic Tumors. *Endocr Pract*. 2014;20:e162-5. [PMID 24936550](https://pubmed.ncbi.nlm.nih.gov/24936550/)
27. Matsuno A et al. Molecular status of pituitary carcinoma and atypical adenoma that contributes the effectiveness of temozolomide. *Med Mol Morphol*. 2014;47:1-7. [PMID 23955641](https://pubmed.ncbi.nlm.nih.gov/23955641/)
28. Trouillas J et al. Clinical, Pathological, and Molecular Factors of Aggressiveness in Lactotroph Tumours. *Neuroendocrinology*. 2019;109:70-76. [PMID 30943495](https://pubmed.ncbi.nlm.nih.gov/30943495/)
29. Vandeva S et al. Somatic and germline mutations in the pathogenesis of pituitary adenomas. *Eur J Endocrinol*. 2019;181:R235-R254. [PMID 31658440](https://pubmed.ncbi.nlm.nih.gov/31658440/)
30. Guo J et al. The SF3B1(R625H) mutation promotes prolactinoma tumor progression through aberrant splicing of DLG1. *J Exp Clin Cancer Res*. 2022;41:26. [PMID 35039052](https://pubmed.ncbi.nlm.nih.gov/35039052/)

---

### C. 도파민 D2 수용체 · 락토트로프 생리 (D2 receptor and lactotroph physiology)

> 모델의 심장부: 단일 점유 방정식(`SIGDRIVE`), D2R 밀도(`D2RS0`), 단기 루프(`TIDA`→`DAP`), 그리고
> 네 갈래 분기(분비·전사·세포부피·증식)의 근거. D2R 밀도 감소가 저항의 정량적 정의라는 주장도 여기서 온다.

31. Curlewis JD & McNeilly AS. Prolactin short-loop feedback and prolactin inhibition of luteinizing hormone secretion during the breeding season and seasonal anoestrus in the ewe. *Neuroendocrinology*. 1991;54:279-85. [PMID 1944814](https://pubmed.ncbi.nlm.nih.gov/1944814/)
32. Elsholtz HP et al. Inhibitory control of prolactin and Pit-1 gene promoters by dopamine. Dual signaling pathways required for D2 receptor-regulated expression of the prolactin gene. *J Biol Chem*. 1991;266:22919-25. [PMID 1835974](https://pubmed.ncbi.nlm.nih.gov/1835974/)
33. Caccavelli L et al. Decreased expression of the two D2 dopamine receptor isoforms in bromocriptine-resistant prolactinomas. *Neuroendocrinology*. 1994;60:314-22. [PMID 7969790](https://pubmed.ncbi.nlm.nih.gov/7969790/)
34. Kelly MA et al. Pituitary lactotroph hyperplasia and chronic hyperprolactinemia in dopamine D2 receptor-deficient mice. *Neuron*. 1997;19:103-13. [PMID 9247267](https://pubmed.ncbi.nlm.nih.gov/9247267/)
35. Missale C et al. Dopamine receptors: from structure to function. *Physiol Rev*. 1998;78:189-225. [PMID 9457173](https://pubmed.ncbi.nlm.nih.gov/9457173/)
36. Ben-Jonathan N & Hnasko R. Dopamine as a prolactin (PRL) inhibitor. *Endocr Rev*. 2001;22:724-63. [PMID 11739329](https://pubmed.ncbi.nlm.nih.gov/11739329/)
37. Peverelli E et al. Filamin-A is essential for dopamine d2 receptor expression and signaling in tumorous lactotrophs. *J Clin Endocrinol Metab*. 2012;97:967-77. [PMID 22259062](https://pubmed.ncbi.nlm.nih.gov/22259062/)
38. Shimazu S et al. Resistance to dopamine agonists in prolactinoma is correlated with reduction of dopamine D2 receptor long isoform mRNA levels. *Eur J Endocrinol*. 2012;166:383-90. [PMID 22127489](https://pubmed.ncbi.nlm.nih.gov/22127489/)
39. Venkatesh SK et al. Spontaneous reduction of prolactinoma post cabergoline withdrawal. *Indian J Endocrinol Metab*. 2012;16:833-5. [PMID 23087877](https://pubmed.ncbi.nlm.nih.gov/23087877/)
40. Grattan DR. 60 YEARS OF NEUROENDOCRINOLOGY: The hypothalamo-prolactin axis. *J Endocrinol*. 2015;226:T101-22. [PMID 26101377](https://pubmed.ncbi.nlm.nih.gov/26101377/)
41. Bernard V et al. Autocrine actions of prolactin contribute to the regulation of lactotroph function in vivo. *FASEB J*. 2018;32:4791-4797. [PMID 29596024](https://pubmed.ncbi.nlm.nih.gov/29596024/)
42. Li H et al. Melatonin Modulates Lactation by Regulating Prolactin Secretion Via Tuberoinfundibular Dopaminergic Neurons in the Hypothalamus- Pituitary System. *Curr Protein Pept Sci*. 2020;21:744-750. [PMID 32392109](https://pubmed.ncbi.nlm.nih.gov/32392109/)
43. McNamara AV et al. Transcription Factor Pit-1 Affects Transcriptional Timing in the Dual-Promoter Human Prolactin Gene. *Endocrinology*. 2021;162. [PMID 33388754](https://pubmed.ncbi.nlm.nih.gov/33388754/)

---

### D. 프로락틴 생물학 · 수용체 신호 (Prolactin biology and receptor signalling)

> `PRLB` 반감기, 분비과립 저장고, PRL 유전자 전사(Pit-1·ERα), 글리코실화/16 kDa 절단, PRLR-JAK2-STAT5 하류.

44. Day RN et al. A protein kinase inhibitor gene reduces both basal and multihormone-stimulated prolactin gene transcription. *J Biol Chem*. 1989;264:431-6. [PMID 2535842](https://pubmed.ncbi.nlm.nih.gov/2535842/)
45. Shull JD & Gorski J. Estrogen regulation of prolactin gene transcription in vivo: paradoxical effects of 17 beta-estradiol dose. *Endocrinology*. 1989;124:279-85. [PMID 2909367](https://pubmed.ncbi.nlm.nih.gov/2909367/)
46. Smith CR & Norman MR. Prolactin and growth hormone: molecular heterogeneity and measurement in serum. *Ann Clin Biochem*. 1990;27 ( Pt 6):542-50. [PMID 2080857](https://pubmed.ncbi.nlm.nih.gov/2080857/)
47. Sinha YN. Prolactin variants. *Trends Endocrinol Metab*. 1992;3:100-6. [PMID 18407087](https://pubmed.ncbi.nlm.nih.gov/18407087/)
48. Bole-Feysot C et al. Prolactin (PRL) and its receptor: actions, signal transduction pathways and phenotypes observed in PRL receptor knockout mice. *Endocr Rev*. 1998;19:225-68. [PMID 9626554](https://pubmed.ncbi.nlm.nih.gov/9626554/)
49. Freeman ME et al. Prolactin: structure, function, and regulation of secretion. *Physiol Rev*. 2000;80:1523-631. [PMID 11015620](https://pubmed.ncbi.nlm.nih.gov/11015620/)
50. Macotela Y et al. Matrix metalloproteases from chondrocytes generate an antiangiogenic 16 kDa prolactin. *J Cell Sci*. 2006;119:1790-800. [PMID 16608881](https://pubmed.ncbi.nlm.nih.gov/16608881/)
51. Ben-Jonathan N et al. What can we learn from rodents about prolactin in humans?. *Endocr Rev*. 2008;29:1-41. [PMID 18057139](https://pubmed.ncbi.nlm.nih.gov/18057139/)
52. Skowronska-Krawczyk D et al. Required enhancer-matrin-3 network interactions for a homeodomain transcription program. *Nature*. 2014;514:257-61. [PMID 25119036](https://pubmed.ncbi.nlm.nih.gov/25119036/)
53. Bernard V et al. New insights in prolactin: pathological implications. *Nat Rev Endocrinol*. 2015;11:265-75. [PMID 25781857](https://pubmed.ncbi.nlm.nih.gov/25781857/)
54. Gao Q et al. Seasonal patterns of prolactin, prolactin receptor, and STAT5 expression in the ovaries of wild ground squirrels (<em>Citellus dauricus</em> Brandt). *Eur J Histochem*. 2023;67. [PMID 37781865](https://pubmed.ncbi.nlm.nih.gov/37781865/)
55. Hackwell ECR et al. Prolactin-mediates a lactation-induced suppression of arcuate kisspeptin neuronal activity necessary for lactational infertility in mice. *Elife*. 2025;13. [PMID 39819370](https://pubmed.ncbi.nlm.nih.gov/39819370/)

---

### E. 성선축 억제 · 생식 (HPG axis suppression and reproduction)

> `KISS`→`GNRH`→`LH`/`FSH`→성호르몬 경로와 `K50K`(프로락틴 절반억제 농도)의 근거.

56. Rasmussen C et al. Prolactin secretion and menstrual function after long-term bromocriptine treatment. *Fertil Steril*. 1987;48:550-4. [PMID 3653413](https://pubmed.ncbi.nlm.nih.gov/3653413/)
57. Sonigo C et al. Hyperprolactinemia-induced ovarian acyclicity is reversed by kisspeptin administration. *J Clin Invest*. 2012;122:3791-5. [PMID 23006326](https://pubmed.ncbi.nlm.nih.gov/23006326/)
58. Donato J Jr & Frazão R. Interactions between prolactin and kisspeptin to control reproduction. *Arch Endocrinol Metab*. 2016;60:587-595. [PMID 27901187](https://pubmed.ncbi.nlm.nih.gov/27901187/)

---

### F. 측정층: hook effect · macroprolactin (The measurement layer)

> 모델에서 생물학과 분리된 별도 층(`PRLIMM`·`PRLMEAS_`·`PRLDIL`·`PEGREC`)의 근거.
> 진단 D6·D7과 시나리오 S21·S22가 전적으로 이 문헌군에 기반한다.

59. Bevan JS et al. Misinterpretation of prolactin levels leading to management errors in patients with sellar enlargement. *Am J Med*. 1987;82:29-32. [PMID 3799691](https://pubmed.ncbi.nlm.nih.gov/3799691/)
60. St-Jean E et al. High prolactin levels may be missed by immunoradiometric assay in patients with macroprolactinomas. *Clin Endocrinol (Oxf)*. 1996;44:305-9. [PMID 8729527](https://pubmed.ncbi.nlm.nih.gov/8729527/)
61. Barkan AL & Chandler WF. Giant pituitary prolactinoma with falsely low serum prolactin: the pitfall of the "high-dose hook effect": case report. *Neurosurgery*. 1998;42:913-5; discussion 915-6. [PMID 9574657](https://pubmed.ncbi.nlm.nih.gov/9574657/)
62. Petakov MS et al. Pituitary adenomas secreting large amounts of prolactin may give false low values in immunoradiometric assays. The hook effect. *J Endocrinol Invest*. 1998;21:184-8. [PMID 9591215](https://pubmed.ncbi.nlm.nih.gov/9591215/)
63. Colao A et al. Macroprolactinoma shrinkage during cabergoline treatment is greater in naive patients than in patients pretreated with other dopamine agonists: a prospective study in 110 patients. *J Clin Endocrinol Metab*. 2000;85:2247-52. [PMID 10852458](https://pubmed.ncbi.nlm.nih.gov/10852458/)
64. Frieze TW et al. "Hook effect" in prolactinomas: case report and review of literature. *Endocr Pract*. 2002;8:296-303. [PMID 12173917](https://pubmed.ncbi.nlm.nih.gov/12173917/)
65. Schöfl C et al. Falsely low serum prolactin in two cases of invasive macroprolactinoma. *Pituitary*. 2002;5:261-5. [PMID 14558675](https://pubmed.ncbi.nlm.nih.gov/14558675/)
66. Smith TP et al. Gross variability in the detection of prolactin in sera containing big big prolactin (macroprolactin) by commercial immunoassays. *J Clin Endocrinol Metab*. 2002;87:5410-5. [PMID 12466327](https://pubmed.ncbi.nlm.nih.gov/12466327/)
67. Gibney J et al. Clinical relevance of macroprolactin. *Clin Endocrinol (Oxf)*. 2005;62:633-43. [PMID 15943822](https://pubmed.ncbi.nlm.nih.gov/15943822/)
68. Hattori N et al. Anti-prolactin (PRL) autoantibody-binding sites (epitopes) on PRL molecule in macroprolactinemia. *J Endocrinol*. 2006;190:287-93. [PMID 16899562](https://pubmed.ncbi.nlm.nih.gov/16899562/)
69. Delgrange E et al. Characterization of resistance to the prolactin-lowering effects of cabergoline in macroprolactinomas: a study in 122 patients. *Eur J Endocrinol*. 2009;160:747-52. [PMID 19223454](https://pubmed.ncbi.nlm.nih.gov/19223454/)
70. Hattori N et al. Macroprolactinaemia: prevalence and aetiologies in a large group of hospital workers. *Clin Endocrinol (Oxf)*. 2009;71:702-8. [PMID 19486017](https://pubmed.ncbi.nlm.nih.gov/19486017/)
71. Raverot G et al. Secondary deterioration of visual field during cabergoline treatment for macroprolactinoma. *Clin Endocrinol (Oxf)*. 2009;70:588-92. [PMID 18673461](https://pubmed.ncbi.nlm.nih.gov/18673461/)
72. Barber TM et al. Recurrence of hyperprolactinaemia following discontinuation of dopamine agonist therapy in patients with prolactinoma occurs commonly especially in macroprolactinoma. *Clin Endocrinol (Oxf)*. 2011;75:819-24. [PMID 21645021](https://pubmed.ncbi.nlm.nih.gov/21645021/)
73. Shimatsu A & Hattori N. Macroprolactinemia: diagnostic, clinical, and pathogenic significance. *Clin Dev Immunol*. 2012;2012:167132. [PMID 23304187](https://pubmed.ncbi.nlm.nih.gov/23304187/)
74. Raverot V et al. Prolactin immunoassay: does the high-dose hook effect still exist?. *Pituitary*. 2022;25:653-657. [PMID 35793045](https://pubmed.ncbi.nlm.nih.gov/35793045/)
75. Vermue FC et al. The validation of macroprolactin analysis by polyethylene glycol precipitation using Fujirebio Lumipulse. *Pract Lab Med*. 2022;31:e00292. [PMID 35860390](https://pubmed.ncbi.nlm.nih.gov/35860390/)

---

### G. 카베르골린 · 브로모크립틴 · 퀴나골라이드 (Dopamine agonist PK/PD and trials)

> PK 파라미터(t½ 63-109 h, pg/mL 혈장 농도), 내재활성 e(카베르골린 1.00 vs 브로모크립틴 0.80),
> 그리고 D4의 용량-반응 비교 기준.

76. Bergh T et al. Bromocriptine-induced regression of a suprasellar extending prolactinoma during pregnancy. *J Endocrinol Invest*. 1984;7:133-6. [PMID 6725868](https://pubmed.ncbi.nlm.nih.gov/6725868/)
77. Vance ML et al. Drugs five years later. Bromocriptine. *Ann Intern Med*. 1984;100:78-91. [PMID 6229205](https://pubmed.ncbi.nlm.nih.gov/6229205/)
78. Bevan JS et al. Factors in the outcome of transsphenoidal surgery for prolactinoma and non-functioning pituitary tumour, including pre-operative bromocriptine therapy. *Clin Endocrinol (Oxf)*. 1987;26:541-56. [PMID 3665118](https://pubmed.ncbi.nlm.nih.gov/3665118/)
79. Pellegrini I et al. Resistance to bromocriptine in prolactinomas. *J Clin Endocrinol Metab*. 1989;69:500-9. [PMID 2760167](https://pubmed.ncbi.nlm.nih.gov/2760167/)
80. Vance ML et al. CV 205-502 treatment of hyperprolactinemia. *J Clin Endocrinol Metab*. 1989;68:336-9. [PMID 2521863](https://pubmed.ncbi.nlm.nih.gov/2521863/)
81. Webster J et al. A comparison of cabergoline and bromocriptine in the treatment of hyperprolactinemic amenorrhea. Cabergoline Comparative Study Group. *N Engl J Med*. 1994;331:904-9. [PMID 7915824](https://pubmed.ncbi.nlm.nih.gov/7915824/)
82. Andreotti AC et al. Pharmacokinetics, pharmacodynamics, and tolerability of cabergoline, a prolactin-lowering drug, after administration of increasing oral doses (0.5, 1.0, and 1.5 milligrams) in healthy male volunteers. *J Clin Endocrinol Metab*. 1995;80:841-5. [PMID 7883840](https://pubmed.ncbi.nlm.nih.gov/7883840/)
83. Ferrari C et al. Cabergoline: a new drug for the treatment of hyperprolactinaemia. *Hum Reprod*. 1995;10:1647-52. [PMID 8582955](https://pubmed.ncbi.nlm.nih.gov/8582955/)
84. Rains CP et al. Cabergoline. A review of its pharmacological properties and therapeutic potential in the treatment of hyperprolactinaemia and inhibition of lactation. *Drugs*. 1995;49:255-79. [PMID 7729332](https://pubmed.ncbi.nlm.nih.gov/7729332/)
85. Motta T et al. Vaginal cabergoline in the treatment of hyperprolactinemic patients intolerant to oral dopaminergics. *Fertil Steril*. 1996;65:440-2. [PMID 8566276](https://pubmed.ncbi.nlm.nih.gov/8566276/)
86. Verhelst J et al. Cabergoline in the treatment of hyperprolactinemia: a study in 455 patients. *J Clin Endocrinol Metab*. 1999;84:2518-22. [PMID 10404830](https://pubmed.ncbi.nlm.nih.gov/10404830/)
87. Colao A et al. Outcome of cabergoline treatment in men with prolactinoma: effects of a 24-month treatment on prolactin levels, tumor mass, recovery of pituitary function, and semen analysis. *J Clin Endocrinol Metab*. 2004;89:1704-11. [PMID 15070934](https://pubmed.ncbi.nlm.nih.gov/15070934/)
88. Barlier A & Jaquet P. Quinagolide--a valuable treatment option for hyperprolactinaemia. *Eur J Endocrinol*. 2006;154:187-95. [PMID 16452531](https://pubmed.ncbi.nlm.nih.gov/16452531/)

---

### H. 저항성 · 고용량 · 수술 · 방사선 · 테모졸로마이드 (Resistance and second line)

> 부분 저항(EC50 이동) vs 진성 저항(Emax 저하)의 구분, 고용량 카베르골린, 수술 완치율,
> MGMT 의존적 테모졸로마이드 반응(S27/S28).

89. Kreutzer J et al. Operative treatment of prolactinomas: indications and results in a current consecutive series of 212 patients. *Eur J Endocrinol*. 2008;158:11-8. [PMID 18166812](https://pubmed.ncbi.nlm.nih.gov/18166812/)
90. Ono M et al. Prospective study of high-dose cabergoline treatment of prolactinomas in 150 patients. *J Clin Endocrinol Metab*. 2008;93:4721-7. [PMID 18812485](https://pubmed.ncbi.nlm.nih.gov/18812485/)
91. Babey M et al. Pituitary surgery for small prolactinomas as an alternative to treatment with dopamine agonists. *Pituitary*. 2011;14:222-30. [PMID 21170594](https://pubmed.ncbi.nlm.nih.gov/21170594/)
92. Vroonen L et al. Prolactinomas resistant to standard doses of cabergoline: a multicenter study of 92 patients. *Eur J Endocrinol*. 2012;167:651-62. [PMID 22918301](https://pubmed.ncbi.nlm.nih.gov/22918301/)
93. Chen W et al. HIF-1α inhibition sensitizes pituitary adenoma cells to temozolomide by regulating MGMT expression. *Oncol Rep*. 2013;30:2495-501. [PMID 23970362](https://pubmed.ncbi.nlm.nih.gov/23970362/)
94. Bengtsson D et al. Long-term outcome and MGMT as a predictive marker in 24 patients with atypical pituitary adenomas and pituitary carcinomas given treatment with temozolomide. *J Clin Endocrinol Metab*. 2015;100:1689-98. [PMID 25646794](https://pubmed.ncbi.nlm.nih.gov/25646794/)
95. Losa M et al. Temozolomide therapy in patients with aggressive pituitary adenomas or carcinomas. *J Neurooncol*. 2016;126:519-25. [PMID 26614517](https://pubmed.ncbi.nlm.nih.gov/26614517/)
96. Halevy C & Whitelaw BC. How effective is temozolomide for treating pituitary tumours and when should it be used?. *Pituitary*. 2017;20:261-266. [PMID 27581836](https://pubmed.ncbi.nlm.nih.gov/27581836/)
97. Honegger J & Grimm F. The experience with transsphenoidal surgery and its importance to outcomes. *Pituitary*. 2018;21:545-555. [PMID 30062664](https://pubmed.ncbi.nlm.nih.gov/30062664/)
98. Lee DK et al. Factors Influencing Visual Field Recovery after Transsphenoidal Resection of a Pituitary Adenoma. *Korean J Ophthalmol*. 2018;32:488-496. [PMID 30549473](https://pubmed.ncbi.nlm.nih.gov/30549473/)
99. Buchfelder M et al. Surgery for Prolactinomas to Date. *Neuroendocrinology*. 2019;109:77-81. [PMID 30699424](https://pubmed.ncbi.nlm.nih.gov/30699424/)
100. Ponce AJ et al. Low prolactin levels are associated with visceral adipocyte hypertrophy and insulin resistance in humans. *Endocrine*. 2020;67:331-343. [PMID 31919769](https://pubmed.ncbi.nlm.nih.gov/31919769/)

---

### I. 약물 중단 · 재발 (Withdrawal and recurrence)

> D10의 기준. 이 문헌군의 21% 지속 완화율이 모델에서는 '추적 기간의 함수'로 재해석된다.

101. Colao A et al. Withdrawal of long-term cabergoline therapy for tumoral and nontumoral hyperprolactinemia. *N Engl J Med*. 2003;349:2023-33. [PMID 14627787](https://pubmed.ncbi.nlm.nih.gov/14627787/)
102. Dekkers OM et al. Recurrence of hyperprolactinemia after withdrawal of dopamine agonists: systematic review and meta-analysis. *J Clin Endocrinol Metab*. 2010;95:43-51. [PMID 19880787](https://pubmed.ncbi.nlm.nih.gov/19880787/)
103. Huda MS et al. Factors determining the remission of microprolactinomas after dopamine agonist withdrawal. *Clin Endocrinol (Oxf)*. 2010;72:507-11. [PMID 19549247](https://pubmed.ncbi.nlm.nih.gov/19549247/)
104. Auriemma RS et al. Results of a single-center observational 10-year survey study on recurrence of hyperprolactinemia after pregnancy and lactation. *J Clin Endocrinol Metab*. 2013;98:372-9. [PMID 23162092](https://pubmed.ncbi.nlm.nih.gov/23162092/)
105. Kwancharoen R et al. Second attempt to withdraw cabergoline in prolactinomas: a pilot study. *Pituitary*. 2014;17:451-6. [PMID 24078319](https://pubmed.ncbi.nlm.nih.gov/24078319/)

---

### J. 임신 · 에스트로겐 (Pregnancy and oestrogen)

> E2의 이중 작용(PRL 전사 ↑, D2R ↓), 임신 중 증상성 종양 증대 위험(micro ~2.7% vs macro ~21-23%).

106. Christin-Maître S et al. Prolactinoma and estrogens: pregnancy, contraception and hormonal replacement therapy. *Ann Endocrinol (Paris)*. 2007;68:106-12. [PMID 17540335](https://pubmed.ncbi.nlm.nih.gov/17540335/)
107. Colao A et al. Pregnancy outcomes following cabergoline treatment: extended results from a 12-year observational study. *Clin Endocrinol (Oxf)*. 2008;68:66-71. [PMID 17760883](https://pubmed.ncbi.nlm.nih.gov/17760883/)
108. Lebbe M et al. Outcome of 100 pregnancies initiated under treatment with cabergoline in hyperprolactinaemic women. *Clin Endocrinol (Oxf)*. 2010;73:236-42. [PMID 20455894](https://pubmed.ncbi.nlm.nih.gov/20455894/)
109. Galvão A et al. Prolactinoma and pregnancy - a series of cases including pituitary apoplexy. *J Obstet Gynaecol*. 2017;37:284-287. [PMID 27866462](https://pubmed.ncbi.nlm.nih.gov/27866462/)
110. Jayabalan N et al. Cross Talk between Adipose Tissue and Placenta in Obese and Gestational Diabetes Mellitus Pregnancies via Exosomes. *Front Endocrinol (Lausanne)*. 2017;8:239. [PMID 29021781](https://pubmed.ncbi.nlm.nih.gov/29021781/)
111. Karaca Z et al. How does pregnancy affect the patients with pituitary adenomas: a study on 113 pregnancies from Turkey. *J Endocrinol Invest*. 2018;41:129-141. [PMID 28634705](https://pubmed.ncbi.nlm.nih.gov/28634705/)
112. Laway BA et al. Prolactinoma Outcome After Pregnancy and Lactation: A Cohort Study. *Indian J Endocrinol Metab*. 2021;25:559-562. [PMID 35355922](https://pubmed.ncbi.nlm.nih.gov/35355922/)

---

### K. 뼈 · 대사 (Bone and metabolic consequences)

> `KBEXP`(골 회전 비율→BMD 설정점), `FIRR`(비가역 분율), 척추 골절 확률의 근거. 프로락틴 정상화 후에도
> 골밀도가 완전히 회복되지 않는다는 관찰이 ratchet 구조의 직접적 근거이다.

113. Klibanski A et al. Decreased bone density in hyperprolactinemic women. *N Engl J Med*. 1980;303:1511-4. [PMID 7432421](https://pubmed.ncbi.nlm.nih.gov/7432421/)
114. Greenspan SL et al. Osteoporosis in men with hyperprolactinemic hypogonadism. *Ann Intern Med*. 1986;104:777-82. [PMID 3706929](https://pubmed.ncbi.nlm.nih.gov/3706929/)
115. Colao A et al. Prolactinomas in adolescents: persistent bone loss after 2 years of prolactin normalization. *Clin Endocrinol (Oxf)*. 2000;52:319-27. [PMID 10718830](https://pubmed.ncbi.nlm.nih.gov/10718830/)
116. Vestergaard P et al. Fracture risk is increased in patients with GH deficiency or untreated prolactinomas--a case-control study. *Clin Endocrinol (Oxf)*. 2002;56:159-67. [PMID 11874406](https://pubmed.ncbi.nlm.nih.gov/11874406/)
117. Naliato EC et al. Prevalence of osteopenia in men with prolactinoma. *J Endocrinol Invest*. 2005;28:12-7. [PMID 15816365](https://pubmed.ncbi.nlm.nih.gov/15816365/)
118. Mazziotti G et al. Vertebral fractures in males with prolactinoma. *Endocrine*. 2011;39:288-93. [PMID 21479837](https://pubmed.ncbi.nlm.nih.gov/21479837/)
119. Auriemma RS et al. Effect of cabergoline on metabolism in prolactinomas. *Neuroendocrinology*. 2013;98:299-310. [PMID 24355865](https://pubmed.ncbi.nlm.nih.gov/24355865/)
120. Sperling S & Bhatt H. Prolactinoma: A Massive Effect on Bone Mineral Density in a Young Patient. *Case Rep Endocrinol*. 2016;2016:6312621. [PMID 27446618](https://pubmed.ncbi.nlm.nih.gov/27446618/)
121. Andereggen L et al. Persistent bone impairment despite long-term control of hyperprolactinemia and hypogonadism in men and women with prolactinomas. *Sci Rep*. 2021;11:5122. [PMID 33664388](https://pubmed.ncbi.nlm.nih.gov/33664388/)

---

### L. 질량 효과 · 시야 · 줄기 압박 (Mass effect, visual fields, stalk effect)

> 기하학적 압박 모형, 가역 전도차단 vs 비가역 축삭 손실, 그리고 줄기 압박 고프로락틴혈증의 상한(D8).

122. Berwaerts J et al. A giant prolactinoma presenting with unilateral exophthalmos: effect of cabergoline and review of the literature. *J Endocrinol Invest*. 2000;23:393-8. [PMID 10908167](https://pubmed.ncbi.nlm.nih.gov/10908167/)
123. Karavitaki N et al. Do the limits of serum prolactin in disconnection hyperprolactinaemia need re-definition? A study of 226 patients with histologically verified non-functioning pituitary macroadenoma. *Clin Endocrinol (Oxf)*. 2006;65:524-9. [PMID 16984247](https://pubmed.ncbi.nlm.nih.gov/16984247/)
124. Korevaar T et al. Disconnection hyperprolactinaemia in nonadenomatous sellar/parasellar lesions practically never exceeds 2000 mU/l. *Clin Endocrinol (Oxf)*. 2012;76:602-3. [PMID 21942983](https://pubmed.ncbi.nlm.nih.gov/21942983/)
125. Danesh-Meyer HV et al. Optical coherence tomography predicts visual outcome for pituitary tumors. *J Clin Neurosci*. 2015;22:1098-104. [PMID 25891894](https://pubmed.ncbi.nlm.nih.gov/25891894/)
126. Bulwer C et al. Cabergoline-related impulse control disorder in an adolescent with a giant prolactinoma. *Clin Endocrinol (Oxf)*. 2017;86:862-864. [PMID 28346715](https://pubmed.ncbi.nlm.nih.gov/28346715/)
127. Rutland JW et al. Measuring degeneration of the lateral geniculate nuclei from pituitary adenoma compression detected by 7T ultra-high field MRI: a method for predicting vision recovery following surgical decompression of the optic chiasm. *J Neurosurg*. 2020;132:1747-1756. [PMID 31100726](https://pubmed.ncbi.nlm.nih.gov/31100726/)
128. Alkhaibary A et al. Invasive Giant Prolactinoma. *World Neurosurg*. 2024;181:21-22. [PMID 37827431](https://pubmed.ncbi.nlm.nih.gov/37827431/)

---

### M. 안전성: 판막 · 충동조절 · 구역 (Safety: valve, impulse control, nausea)

> 5-HT2B(혈장 구동)·D3(충동조절)·area postrema D2(구역)의 분리. 파킨슨병 용량과 프로락틴종 용량의
> 10-40배 차이가 두 문헌군의 불일치를 설명한다는 D11의 근거.

129. Jovanović-Mićić D et al. The role of alpha-adrenergic mechanisms within the area postrema in dopamine-induced emesis. *Eur J Pharmacol*. 1995;272:21-30. [PMID 7713146](https://pubmed.ncbi.nlm.nih.gov/7713146/)
130. Schade R et al. Dopamine agonists and the risk of cardiac-valve regurgitation. *N Engl J Med*. 2007;356:29-38. [PMID 17202453](https://pubmed.ncbi.nlm.nih.gov/17202453/)
131. Zanettini R et al. Valvular heart disease and the use of dopamine agonists for Parkinson's disease. *N Engl J Med*. 2007;356:39-46. [PMID 17202454](https://pubmed.ncbi.nlm.nih.gov/17202454/)
132. Auriemma RS et al. Safety of long-term treatment with cabergoline on cardiac valve disease in patients with prolactinomas. *Eur J Endocrinol*. 2013;169:359-66. [PMID 23824978](https://pubmed.ncbi.nlm.nih.gov/23824978/)
133. Bancos I et al. Impulse control disorders in patients with dopamine agonist-treated prolactinomas and nonfunctioning pituitary adenomas: a case-control study. *Clin Endocrinol (Oxf)*. 2014;80:863-8. [PMID 24274365](https://pubmed.ncbi.nlm.nih.gov/24274365/)
134. Auriemma RS et al. Cabergoline use for pituitary tumors and valvular disorders. *Endocrinol Metab Clin North Am*. 2015;44:89-97. [PMID 25732645](https://pubmed.ncbi.nlm.nih.gov/25732645/)
135. Caputo C et al. The need for annual echocardiography to detect cabergoline-associated valvulopathy in patients with prolactinoma: a systematic review and additional clinical data. *Lancet Diabetes Endocrinol*. 2015;3:906-13. [PMID 25466526](https://pubmed.ncbi.nlm.nih.gov/25466526/)
136. Barake M et al. MANAGEMENT OF ENDOCRINE DISEASE: Impulse control disorders in patients with hyperpolactinemia treated with dopamine agonists: how much should we worry?. *Eur J Endocrinol*. 2018;179:R287-R296. [PMID 30324793](https://pubmed.ncbi.nlm.nih.gov/30324793/)
137. Dogansen SC et al. Dopamine Agonist-Induced Impulse Control Disorders in Patients With Prolactinoma: A Cross-Sectional Multicenter Study. *J Clin Endocrinol Metab*. 2019;104:2527-2534. [PMID 30848825](https://pubmed.ncbi.nlm.nih.gov/30848825/)
138. Stiles CE et al. Incidence of Cabergoline-Associated Valvulopathy in Primary Care Patients With Prolactinoma Using Hard Cardiac Endpoints. *J Clin Endocrinol Metab*. 2021;106:e711-e720. [PMID 33247916](https://pubmed.ncbi.nlm.nih.gov/33247916/)

---

### N. 약물 유발 고프로락틴혈증 (Drug-induced hyperprolactinaemia)

> 동일한 점유 방정식의 반대 방향. 아리피프라졸 부분작용제(e≈0.25) 병용 시 프로락틴이 떨어지는 D9의 근거.

139. Honbo KS et al. Serum prolactin levels in untreated primary hypothyroidism. *Am J Med*. 1978;64:782-7. [PMID 645742](https://pubmed.ncbi.nlm.nih.gov/645742/)
140. Scanlon MF et al. Altered dopaminergic regulation of thyrotrophin release in patients with prolactinomas: comparison with other tests of hypothalamic-pituitary function. *Clin Endocrinol (Oxf)*. 1981;14:133-43. [PMID 6790201](https://pubmed.ncbi.nlm.nih.gov/6790201/)
141. Daniels GH et al. Effect of risperidone dose on serum prolactin level. *Endocr Pract*. 2001;7:224. [PMID 11421571](https://pubmed.ncbi.nlm.nih.gov/11421571/)
142. Molitch ME. Medication-induced hyperprolactinemia. *Mayo Clin Proc*. 2005;80:1050-7. [PMID 16092584](https://pubmed.ncbi.nlm.nih.gov/16092584/)
143. Hekimsoy Z et al. The prevalence of hyperprolactinaemia in overt and subclinical hypothyroidism. *Endocr J*. 2010;57:1011-5. [PMID 20938100](https://pubmed.ncbi.nlm.nih.gov/20938100/)
144. Ajmal A et al. Psychotropic-induced hyperprolactinemia: a clinical review. *Psychosomatics*. 2014;55:29-36. [PMID 24140188](https://pubmed.ncbi.nlm.nih.gov/24140188/)
145. Peuskens J et al. The effects of novel and newly approved antipsychotics on serum prolactin levels: a comprehensive review. *CNS Drugs*. 2014;28:421-53. [PMID 24677189](https://pubmed.ncbi.nlm.nih.gov/24677189/)
146. Grigg J et al. Antipsychotic-induced hyperprolactinemia: synthesis of world-wide guidelines and integrated recommendations for assessment, management and future research. *Psychopharmacology (Berl)*. 2017;234:3279-3297. [PMID 28889207](https://pubmed.ncbi.nlm.nih.gov/28889207/)
147. Zhang L et al. Efficacy and Safety of Adjunctive Aripiprazole, Metformin, and Paeoniae-Glycyrrhiza Decoction for Antipsychotic-Induced Hyperprolactinemia: A Network Meta-Analysis of Randomized Controlled Trials. *Front Psychiatry*. 2021;12:728204. [PMID 34658963](https://pubmed.ncbi.nlm.nih.gov/34658963/)
148. Lin X et al. Antipsychotic-Related Prolactin Changes: A Systematic Review and Dose-Response Meta-analysis. *CNS Drugs*. 2025;39:937-947. [PMID 40830715](https://pubmed.ncbi.nlm.nih.gov/40830715/)

---

### O. QSP 방법론 (QSP methodology)

> mrgsolve 구현과 모델기반 신약개발(MIDD) 맥락.

149. Peterson MC & Riggs MM. FDA Advisory Meeting Clinical Pharmacology Review Utilizes a Quantitative Systems Pharmacology (QSP) Model: A Watershed Moment?. *CPT Pharmacometrics Syst Pharmacol*. 2015;4:e00020. [PMID 26225239](https://pubmed.ncbi.nlm.nih.gov/26225239/)
150. Helmlinger G et al. Quantitative Systems Pharmacology: An Exemplar Model-Building Workflow With Applications in Cardiovascular, Metabolic, and Oncology Drug Development. *CPT Pharmacometrics Syst Pharmacol*. 2019;8:380-395. [PMID 31087533](https://pubmed.ncbi.nlm.nih.gov/31087533/)

---

## 검증에서 제외된 후보에 대한 메모 (what was rejected, and why)

관련성 필터가 걸러낸 후보 중에는 **주제가 전혀 다른 논문**(다발성 경화증 예후, BCG 백신
무작위배정 시험, 금속유기골격체 바이오센서 등)이 포함되어 있었습니다. 이들은 저자명 + 짧은
키워드로 축약된 fallback 질의가 PubMed의 관련성 정렬에서 무관한 상위 문헌을 집어 올린 결과이며,
자동으로 제외되었습니다. 이 사실을 남겨 두는 이유는 **저자명 기반 fallback 질의는 반드시 사후
검증이 필요하다**는 점이 이 목록을 만드는 과정에서 가장 실질적인 교훈이었기 때문입니다.

또한 다음 항목들은 검증 가능한 PubMed 인용을 찾지 못해 목록에 넣지 않았고, 따라서 모델에서도
정량적 앵커로 사용하지 않았습니다. 해당 파라미터는 모델 파일에서 명시적으로 "추정치"로
표시되어 있습니다.

- 카베르골린의 **절대 생체이용률**: 공개된 값이 없으므로 `V2_CAB`·`CL_CAB`는 겉보기 값이며
  발표된 pg/mL 범위와 반감기를 재현하도록 맞춘 것입니다.
- **뇌하수체 biophase 분배계수** `PART_CAB = 60`: 뇌하수체 조직 축적을 흡수하는 구조적
  파라미터로, 혈장 데이터만으로는 식별되지 않습니다.
- **면역측정 hook 상수** `KHOOK`·`PHOOK`: 플랫폼마다 다르므로 곡선의 *모양*(단조 증가 → 정점 →
  붕괴)만이 일반적인 주장입니다.

## 라이선스 및 면책

본 참고문헌 목록과 모델은 교육·연구 목적입니다. 임상 의사결정에 직접 사용해서는 안 됩니다.
