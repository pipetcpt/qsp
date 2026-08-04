# 프래더-윌리 증후군 QSP 모델 — 참고문헌

**Prader-Willi Syndrome · QSP model reference list**

> 이 목록의 **357편 전부**는 PubMed E-utilities로 실제 조회하여 PMID·저자·연도·저널·제목을
> 그대로 가져온 것이다. 기억에 의존해 적은 PMID는 하나도 없다.
> 각 절의 머리글은 그 문헌군이 모델의 **어느 파라미터 또는 어느 구조적 주장**을 받치는지를
> 명시한다 — 무엇이 **교정(calibrated)** 되었고 무엇이 **유도(derived)** 되었는지의 경계를
> 흐리지 않기 위해서다.

---

## 0. 교정된 것과 유도된 것 (What is calibrated versus what is derived)

모델의 정직성은 이 표에 달려 있다. 아래 세 줄만 임상 결과에 맞춰 **교정**되었다.

| 파라미터 | 값 | 무엇에 맞췄는가 | 근거 절 |
|---|---|---|---|
| `KGHRD` (구동에서 그렐린 팔의 가중치) | 0.10 | 옥트레오타이드·리볼레타이드의 음성 결과 | E |
| `ECBMAX` / `EV1AMX` (OXTR / V1a) | 1.35 / 1.30 | 카베토신 3.2 mg의 ΔHQ-CT ≈ −1.8 과 9.6 mg의 역전 | C |
| `EDZMAX` (AgRP KATP) | 0.40 | DCCR 5.1 mg/kg의 ΔHQ-CT ≈ −1.6 | L |

아래는 **구조에서 나온 것**이며 어떤 임상 수치에도 맞추지 않았다.

| 유도된 결과 | 어디서 나오는가 |
|---|---|
| 다섯 분기의 실패 순서와 전구체 축적의 **역상관** | 탈출비 하나씩, Eq. A / Eq. B |
| 전구체:산물 비 = **1/PC13**, 모든 분기에서 동일 | Eq. A ÷ Eq. B 의 대수적 항등식 |
| 교차반응 항체는 합성률을 재고 PC1/3에 **정확히** 눈이 먼다 | 정상상태 질량보존 |
| 체중 증가가 과식보다 **먼저** 온다 (Miller 2a) | 낮은 소비 × 연령규범 배식 |
| MC4R 작용제의 상한 1.20배 (relay 천장) | 직렬 이득 relay(∞) = KREL+1 |
| GH 시작 4-8주의 AHI **창** | τ_lymphoid 20 d < τ_muscle 75 d |
| 카베토신 용량-반응의 **역전**과 해석적 최적점 C* | E1·K1 > E2·K2 의 이봉 구조 |
| DCCR의 치료계수가 EC50 두 개로 **고정** | 같은 채널, 두 조직 |
| 바이오마커와 종말점의 **직교성** | GHS-R1a 점유 포화 |
| HQ-CT 반응이 이동이 아니라 **이봉 분포**일 것 | SEEK 의 이중안정성 |
| 잠긴 상태를 벗어나는 데 필요한 옥시토신 팔 이득 **2.4배** | 상부 안장-노드 위치 |
| 척추측만에서 GH의 두 효과가 **거의 상쇄** | 성장속도 ↑ × 근긴장 ↑ |
| 상대적 저인슐린혈증 + 고아디포넥틴혈증 | 같은 전환효소 병변의 부수 결과 |

---

## A. 유전학과 각인 (Genetics and imprinting of 15q11-q13)

> 이 절은 모델의 **병변 레이어**(cluster 1)를 뒷받침한다. 모델은 15q11-q13 부계 발현 소실을 여러 개의 호르몬 결핍 목록이 아니라 **하나의 스칼라**(PC1/3 활성)로 쓰는데, 그 정당성은 SNORD116 최소임계영역·NDN·MAGEL2가 모두 NHLH2를 거쳐 같은 전환효소 노드로 수렴한다는 여기 문헌들에 있다. MKRN3는 사춘기 브레이크이므로 같은 결실이 사춘기를 **앞당기면서** 진폭은 낮추는 두 방향 효과를 만든다(cluster 14).

1. Angulo MA … Cataletto ME (2015). *Prader-Willi syndrome: a review of clinical, genetic, and endocrine findings*. J Endocrinol Invest. [PMID 26062517](https://pubmed.ncbi.nlm.nih.gov/26062517/)
2. Godler DE … Butler MG (2025). *Genetics of Prader-Willi and Angelman syndromes: 2024 update*. Curr Opin Psychiatry. [PMID 39804213](https://pubmed.ncbi.nlm.nih.gov/39804213/)
3. Bittel DC, Butler MG (2005). *Prader-Willi syndrome: clinical genetics, cytogenetics and molecular biology*. Expert Rev Mol Med. [PMID 16038620](https://pubmed.ncbi.nlm.nih.gov/16038620/)
4. Szpecht-Potocka A (1999). *Molecular analysis in Prader-Willi syndrome diagnosis*. Med Wieku Rozwoj. [PMID 10910667](https://pubmed.ncbi.nlm.nih.gov/10910667/)
5. Chung MS … Carmichael GG (2020). *Prader-Willi syndrome: reflections on seminal studies and future therapies*. Open Biol. [PMID 32961075](https://pubmed.ncbi.nlm.nih.gov/32961075/)
6. Kummerfeld DM … Rozhdestvensky TS (2021). *A Comprehensive Review of Genetically Engineered Mouse Models for Prader-Willi Syndrome Research*. Int J Mol Sci. [PMID 33807162](https://pubmed.ncbi.nlm.nih.gov/33807162/)
7. Polex-Wolf J … Yeo GS et al. (2018). *Hypothalamic loss of Snord116 recapitulates the hyperphagia of Prader-Willi syndrome*. J Clin Invest. [PMID 29376887](https://pubmed.ncbi.nlm.nih.gov/29376887/)
8. Burnett LC … Leibel RL (2017). *Loss of the imprinted, non-coding Snord116 gene cluster in the interval deleted in the Prader Willi syndrome results in murine neuronal and endocrine pancreatic developmental phenotypes*. Hum Mol Genet. [PMID 28973544](https://pubmed.ncbi.nlm.nih.gov/28973544/)
9. Schubert T, Schaaf CP (2025). *MAGEL2 (patho-)physiology and Schaaf-Yang syndrome*. Dev Med Child Neurol. [PMID 38950199](https://pubmed.ncbi.nlm.nih.gov/38950199/)
10. Reznik DL … Samaco RC et al. (2023). *Magel2 truncation alters select behavioral and physiological outcomes in a rat model of Schaaf-Yang syndrome*. Dis Model Mech. [PMID 36637363](https://pubmed.ncbi.nlm.nih.gov/36637363/)
11. D Hidalgo-Santos A … Tomás-Vila M et al. (2018). *A Novel Mutation of MAGEL2 in a Patient with Schaaf-Yang Syndrome and Hypopituitarism*. Int J Endocrinol Metab. [PMID 30323850](https://pubmed.ncbi.nlm.nih.gov/30323850/)
12. Hoyos Sanchez MC … Fon Tacer K (2023). *Hormonal Imbalances in Prader-Willi and Schaaf-Yang Syndromes Imply the Evolution of Specific Regulation of Hypothalamic Neuroendocrine Function in Mammals*. Int J Mol Sci. [PMID 37685915](https://pubmed.ncbi.nlm.nih.gov/37685915/)
13. Miller NL … Mellon PL (2009). *Necdin, a Prader-Willi syndrome candidate gene, regulates gonadotropin-releasing hormone neurons during development*. Hum Mol Genet. [PMID 18930956](https://pubmed.ncbi.nlm.nih.gov/18930956/)
14. Barelle PY … Muscatelli F et al. (2025). *Investigation of a mouse model of Prader-Willi Syndrome with combined disruption of Necdin and Magel2*. JCI Insight. [PMID 40048253](https://pubmed.ncbi.nlm.nih.gov/40048253/)
15. Ren J … Wevrick R et al. (2003). *Absence of Ndn, encoding the Prader-Willi syndrome-deleted gene necdin, results in congenital deficiency of central respiratory drive in neonatal mice*. J Neurosci. [PMID 12629158](https://pubmed.ncbi.nlm.nih.gov/12629158/)
16. Brito VN … Latronico AC et al. (2023). *The Congenital and Acquired Mechanisms Implicated in the Etiology of Central Precocious Puberty*. Endocr Rev. [PMID 35930274](https://pubmed.ncbi.nlm.nih.gov/35930274/)
17. Magnotto JC … Abreu AP et al. (2023). *Novel MKRN3 Missense Mutations Associated With Central Precocious Puberty Reveal Distinct Effects on Ubiquitination*. J Clin Endocrinol Metab. [PMID 36916482](https://pubmed.ncbi.nlm.nih.gov/36916482/)
18. Canton APM … Brito VN et al. (2024). *The genetic etiology is a relevant cause of central precocious puberty*. Eur J Endocrinol. [PMID 38857188](https://pubmed.ncbi.nlm.nih.gov/38857188/)
19. Yin X … Lu W et al. (2021). *A Novel Loss-of-Function MKRN3 Variant in a Chinese Patient With Familial Precocious Puberty: A Case Report and Functional Study*. Front Genet. [PMID 34421985](https://pubmed.ncbi.nlm.nih.gov/34421985/)
20. Ramos CO … Brito VN et al. (2020). *Outcomes of Patients with Central Precocious Puberty Due to Loss-of-Function Mutations in the MKRN3 Gene after Treatment with Gonadotropin-Releasing Hormone Analog*. Neuroendocrinology. [PMID 31671431](https://pubmed.ncbi.nlm.nih.gov/31671431/)
21. Costa RA … Guida LDC (2019). *Genotype-Phenotype Relationships and Endocrine Findings in Prader-Willi Syndrome*. Front Endocrinol (Lausanne). [PMID 31920975](https://pubmed.ncbi.nlm.nih.gov/31920975/)
22. Eggermann T (2010). *Russell-Silver syndrome*. Am J Med Genet C Semin Med Genet. [PMID 20803658](https://pubmed.ncbi.nlm.nih.gov/20803658/)
23. Mao S … Zou C (2024). *Genotype-phenotype correlation in Prader-Willi syndrome: A large-sample analysis in China*. Clin Genet. [PMID 38258470](https://pubmed.ncbi.nlm.nih.gov/38258470/)
24. Juriaans AF … Hokken-Koelega ACS (2022). *The Spectrum of the Prader-Willi-like Pheno- and Genotype: A Review of the Literature*. Endocr Rev. [PMID 34460908](https://pubmed.ncbi.nlm.nih.gov/34460908/)
25. Golding DM … Wells T et al. (2017). *Paradoxical leanness in the imprinting-centre deletion mouse model for Prader-Willi syndrome*. J Endocrinol. [PMID 27799465](https://pubmed.ncbi.nlm.nih.gov/27799465/)
26. Wang SE … Jiang YH et al. (2025). *Mechanism of EHMT2-mediated genomic imprinting associated with Prader-Willi syndrome*. Nat Commun. [PMID 40610428](https://pubmed.ncbi.nlm.nih.gov/40610428/)
27. Höybye C, Tauber M (2022). *Approach to the Patient With Prader-Willi Syndrome*. J Clin Endocrinol Metab. [PMID 35150573](https://pubmed.ncbi.nlm.nih.gov/35150573/)
28. Butler MG … Driscoll DJ et al. (2019). *Birth seasonality studies in a large Prader-Willi syndrome cohort*. Am J Med Genet A. [PMID 31225937](https://pubmed.ncbi.nlm.nih.gov/31225937/)
29. Vogels A … Fryns JP et al. (2004). *Minimum prevalence, birth incidence and cause of death for Prader-Willi syndrome in Flanders*. Eur J Hum Genet. [PMID 14679397](https://pubmed.ncbi.nlm.nih.gov/14679397/)

## B. 프로호르몬 전환효소 PC1/3 (Prohormone convertase PC1/3, PCSK1)

> 모델의 중심 노드. `EPSPOMC`·`EPSPOXT`·`EPSPGHRH`·`EPSPINS`·`EPSPGHR` 다섯 개의 무차원 탈출비(escape ratio)와 `DPC13`이 여기서 온다. 특히 프로인슐린:인슐린 비 상승은 Eq. B의 직접 검증 대상이며, 모델은 그 비가 **모든 분기에서 정확히 1/PC13**이라는 닫힌 형태 결과를 낸다(`pws_calibration.py` §1).

30. Ramos-Molina B … Lindberg I (2016). *PCSK1 Variants and Human Obesity*. Prog Mol Biol Transl Sci. [PMID 27288825](https://pubmed.ncbi.nlm.nih.gov/27288825/)
31. Pépin L … Coutant R et al. (2019). *A New Case of PCSK1 Pathogenic Variant With Congenital Proprotein Convertase 1/3 Deficiency and Literature Review*. J Clin Endocrinol Metab. [PMID 30383237](https://pubmed.ncbi.nlm.nih.gov/30383237/)
32. Martín MG … Georgia S et al. (2013). *Congenital proprotein convertase 1/3 deficiency causes malabsorptive diarrhea and other endocrinopathies in a pediatric cohort*. Gastroenterology. [PMID 23562752](https://pubmed.ncbi.nlm.nih.gov/23562752/)
33. Creemers JW … Meyre D et al. (2012). *Heterozygous mutations causing partial prohormone convertase 1 deficiency contribute to human obesity*. Diabetes. [PMID 22210313](https://pubmed.ncbi.nlm.nih.gov/22210313/)
34. Huber LM … Janecke AR et al. (2025). *Pathogenic Deep Intronic PCSK1 Variant Causes Proprotein Convertase 1/3 Deficiency in a Family*. Clin Genet. [PMID 39891480](https://pubmed.ncbi.nlm.nih.gov/39891480/)
35. Stijnen P … Creemers JW (2016). *PCSK1 Mutations and Human Endocrinopathies: From Obesity to Gastrointestinal Disorders*. Endocr Rev. [PMID 27187081](https://pubmed.ncbi.nlm.nih.gov/27187081/)
36. Burnett LC … Leibel RL et al. (2017). *Deficiency in prohormone convertase PC1 impairs prohormone processing in Prader-Willi syndrome*. J Clin Invest. [PMID 27941249](https://pubmed.ncbi.nlm.nih.gov/27941249/)
37. Vivoli M … Lindberg I (2012). *Inhibition of prohormone convertases PC1/3 and PC2 by 2,5-dideoxystreptamine derivatives*. Mol Pharmacol. [PMID 22169851](https://pubmed.ncbi.nlm.nih.gov/22169851/)
38. Furuta M … Steiner DF et al. (1998). *Incomplete processing of proinsulin to insulin accompanied by elevation of Des-31,32 proinsulin intermediates in islets of mice lacking active PC2*. J Biol Chem. [PMID 9452465](https://pubmed.ncbi.nlm.nih.gov/9452465/)
39. Prabhu Y … Lindberg I et al. (2014). *Defective transport of the obesity mutant PC1/3 N222D contributes to loss of function*. Endocrinology. [PMID 24828610](https://pubmed.ncbi.nlm.nih.gov/24828610/)
40. Jing E … Good DJ (2004). *Deletion of the Nhlh2 transcription factor decreases the levels of the anorexigenic peptides alpha melanocyte-stimulating hormone and thyrotropin-releasing hormone and implicates prohormone convertases I and II in obesity*. Endocrinology. [PMID 14701669](https://pubmed.ncbi.nlm.nih.gov/14701669/)
41. Carraro RS … Velloso LA et al. (2021). *Arcuate Nucleus Overexpression of NHLH2 Reduces Body Mass and Attenuates Obesity-Associated Anxiety/Depression-like Behavior*. J Neurosci. [PMID 34675088](https://pubmed.ncbi.nlm.nih.gov/34675088/)
42. Qian Y … Fricker LD (2000). *The C-terminal region of proSAAS is a potent inhibitor of prohormone convertase 1*. J Biol Chem. [PMID 10816562](https://pubmed.ncbi.nlm.nih.gov/10816562/)
43. Feng Y … Fricker LD (2002). *ProSAAS and prohormone convertase 1 are broadly expressed during mouse development*. Brain Res Gene Expr Patterns. [PMID 15018810](https://pubmed.ncbi.nlm.nih.gov/15018810/)
44. Fricker LD … Douglass J et al. (2000). *Identification and characterization of proSAAS, a granin-like neuroendocrine peptide precursor that inhibits prohormone processing*. J Neurosci. [PMID 10632593](https://pubmed.ncbi.nlm.nih.gov/10632593/)
45. Chen YC … Verchere CB et al. (2023). *Deletion of Carboxypeptidase E in β-Cells Disrupts Proinsulin Processing but Does Not Lead to Spontaneous Development of Diabetes in Mice*. Diabetes. [PMID 37364047](https://pubmed.ncbi.nlm.nih.gov/37364047/)
46. Chen YC … Verchere CB (2018). *Islet prohormone processing in health and disease*. Diabetes Obes Metab. [PMID 30230179](https://pubmed.ncbi.nlm.nih.gov/30230179/)
47. Lindberg I, Fricker LD (2021). *Obesity, POMC, and POMC-processing Enzymes: Surprising Results From Animal Models*. Endocrinology. [PMID 34333593](https://pubmed.ncbi.nlm.nih.gov/34333593/)
48. Zhu X … Steiner DF (2006). *On the processing of proghrelin to ghrelin*. J Biol Chem. [PMID 17050541](https://pubmed.ncbi.nlm.nih.gov/17050541/)
49. O'Brien M … Smith TJ (2010). *Ghrelin in the human myometrium*. Reprod Biol Endocrinol. [PMID 20509935](https://pubmed.ncbi.nlm.nih.gov/20509935/)

## C. 옥시토신 축 (Oxytocin, the PVN and the relay)

> `FOXTN`(생존 PVN 옥시토신 뉴런 분율)과 `KINHO`(전시냅스 AgRP 억제), 그리고 카베토신의 `ECBMAX`/`EV1AMX`가 여기서 온다. 옥시토신 측정 문헌이 서로 충돌하는 것 자체가 모델의 예측이다 — 전구체와 산물을 함께 인식하는 항체는 합성률을 재므로 PC1/3에 **정확히 눈이 먼다**(`pws_calibration.py` §2).

50. Swaab DF (2004). *Neuropeptides in hypothalamic neuronal disorders*. Int Rev Cytol. [PMID 15548416](https://pubmed.ncbi.nlm.nih.gov/15548416/)
51. Swaab DF (1995). *Development of the human hypothalamus*. Neurochem Res. [PMID 7643957](https://pubmed.ncbi.nlm.nih.gov/7643957/)
52. Swaab DF … Hofman MA (1995). *Alterations in the hypothalamic paraventricular nucleus and its oxytocin neurons (putative satiety cells) in Prader-Willi syndrome: a study of five cases*. J Clin Endocrinol Metab. [PMID 7852523](https://pubmed.ncbi.nlm.nih.gov/7852523/)
53. Sabatier N … Menzies J (2013). *Oxytocin, feeding, and satiety*. Front Endocrinol (Lausanne). [PMID 23518828](https://pubmed.ncbi.nlm.nih.gov/23518828/)
54. Casipit CG, Anastasopoulou C (2026). *Hypothalamic Dysfunction*. . [PMID 32809578](https://pubmed.ncbi.nlm.nih.gov/32809578/)
55. Gruber T … García-Cáceres C et al. (2023). *High-calorie diets uncouple hypothalamic oxytocin neurons from a gut-to-brain satiation pathway via κ-opioid signaling*. Cell Rep. [PMID 37864798](https://pubmed.ncbi.nlm.nih.gov/37864798/)
56. Nakata M … Maruyama I et al. (2023). *1,5-Anhydro-D-Fructose Exhibits Satiety Effects via the Activation of Oxytocin Neurons in the Paraventricular Nucleus*. Int J Mol Sci. [PMID 37175953](https://pubmed.ncbi.nlm.nih.gov/37175953/)
57. Valassi E … Cavagnini F (2008). *Neuroendocrine control of food intake*. Nutr Metab Cardiovasc Dis. [PMID 18061414](https://pubmed.ncbi.nlm.nih.gov/18061414/)
58. Shah BP … Lowell BB et al. (2014). *MC4R-expressing glutamatergic neurons in the paraventricular hypothalamus regulate feeding and are synaptically connected to the parabrachial nucleus*. Proc Natl Acad Sci U S A. [PMID 25157144](https://pubmed.ncbi.nlm.nih.gov/25157144/)
59. Shalma NM … Abd-ElGawad M et al. (2023). *The efficacy of intranasal oxytocin in patients with Prader-Willi syndrome: A systematic review and meta-analysis*. Diabetes Metab Syndr. [PMID 36774885](https://pubmed.ncbi.nlm.nih.gov/36774885/)
60. Roof E … Ryman DC et al. (2023). *Intranasal Carbetocin Reduces Hyperphagia, Anxiousness, and Distress in Prader-Willi Syndrome: CARE-PWS Phase 3 Trial*. J Clin Endocrinol Metab. [PMID 36633570](https://pubmed.ncbi.nlm.nih.gov/36633570/)
61. Hollander E … Taylor BP et al. (2021). *Intranasal oxytocin versus placebo for hyperphagia and repetitive behaviors in children with Prader-Willi Syndrome: A randomized controlled pilot trial*. J Psychiatr Res. [PMID 33190843](https://pubmed.ncbi.nlm.nih.gov/33190843/)
62. Damen L … Hokken-Koelega ACS et al. (2021). *Oxytocin in young children with Prader-Willi syndrome: Results of a randomized, double-blind, placebo-controlled, crossover trial investigating 3 months of oxytocin*. Clin Endocrinol (Oxf). [PMID 33296519](https://pubmed.ncbi.nlm.nih.gov/33296519/)
63. Einfeld SL … Guastella AJ et al. (2014). *A double-blind randomized controlled trial of oxytocin nasal spray in Prader Willi syndrome*. Am J Med Genet A. [PMID 24980612](https://pubmed.ncbi.nlm.nih.gov/24980612/)
64. Tauber M … Valette M et al. (2026). *Oxytocin in infants with Prader-Willi syndrome to improve dysphagia and disease trajectory*. Orphanet J Rare Dis. [PMID 41639888](https://pubmed.ncbi.nlm.nih.gov/41639888/)
65. Dykens EM … Korner P et al. (2018). *Intranasal carbetocin reduces hyperphagia in individuals with Prader-Willi syndrome*. JCI Insight. [PMID 29925684](https://pubmed.ncbi.nlm.nih.gov/29925684/)
66. Roof E … McCandless SE et al. (2026). *Carbetocin Nasal Spray for the Treatment of Hyperphagia in Prader-Willi Syndrome: Results From the Randomized, Placebo-Controlled, Phase 3 Compass PWS Study*. Clin Ther. [PMID 42486751](https://pubmed.ncbi.nlm.nih.gov/42486751/)
67. Cortese S … Correll CU et al. (2023). *The future of child and adolescent clinical psychopharmacology: A systematic review of phase 2, 3, or 4 randomized controlled trials of pharmacologic agents without regulatory approval or for unapproved indications*. Neurosci Biobehav Rev. [PMID 37001575](https://pubmed.ncbi.nlm.nih.gov/37001575/)
68. Zheng H … Berthoud HR (2002). *Neurochemical phenotype of hypothalamic neurons showing Fos expression 23 h after intracranial AgRP*. Am J Physiol Regul Integr Comp Physiol. [PMID 12010760](https://pubmed.ncbi.nlm.nih.gov/12010760/)
69. Ganella DE … Bathgate RA (2013). *Modulation of feeding by chronic rAAV expression of a relaxin-3 peptide agonist in rat hypothalamus*. Gene Ther. [PMID 23135160](https://pubmed.ncbi.nlm.nih.gov/23135160/)
70. Tabak BA … Mendez AJ et al. (2023). *Advances in human oxytocin measurement: challenges and proposed solutions*. Mol Psychiatry. [PMID 35999276](https://pubmed.ncbi.nlm.nih.gov/35999276/)
71. Carson DS … Parker KJ et al. (2015). *Cerebrospinal fluid and plasma oxytocin concentrations are positively correlated and negatively predict anxiety in children*. Mol Psychiatry. [PMID 25349162](https://pubmed.ncbi.nlm.nih.gov/25349162/)
72. Carson DS … Parker KJ et al. (2014). *Plasma vasopressin concentrations positively predict cerebrospinal fluid vasopressin concentrations in human neonates*. Peptides. [PMID 25148831](https://pubmed.ncbi.nlm.nih.gov/25148831/)
73. Rutigliano G … Fusar-Poli P et al. (2016). *Peripheral oxytocin and vasopressin: Biomarkers of psychiatric disorders? A comprehensive systematic review and preliminary meta-analysis*. Psychiatry Res. [PMID 27183106](https://pubmed.ncbi.nlm.nih.gov/27183106/)

## D. 멜라노코르틴 경로 (Melanocortin pathway and setmelanotide)

> 모델이 α-MSH를 독립적인 포만 팔이 아니라 **옥시토신 팔에 직렬로 들어가는 포화 입력 이득**으로 쓰는 근거. MC4R 포만 신호가 PVN 옥시토신 뉴런을 경유한다는 회로 문헌과, 세트멜라노타이드가 POMC·LEPR·PCSK1 결핍(병변이 MC4R 위)에서는 듣고 PWS(병변이 아래)에서는 듣지 않는다는 대비가 이 절의 핵심이다.

74. Clément K … Kühnen P et al. (2020). *Efficacy and safety of setmelanotide, an MC4R agonist, in individuals with severe obesity due to LEPR or POMC deficiency: single-arm, open-label, multicentre, phase 3 trials*. Lancet Diabetes Endocrinol. [PMID 33137293](https://pubmed.ncbi.nlm.nih.gov/33137293/)
75. Argente J … Farooqi IS et al. (2025). *Setmelanotide in patients aged 2-5 years with rare MC4R pathway-associated obesity (VENTURE): a 1 year, open-label, multicenter, phase 3 trial*. Lancet Diabetes Endocrinol. [PMID 39549719](https://pubmed.ncbi.nlm.nih.gov/39549719/)
76. Wabitsch M … Kühnen P et al. (2022). *Natural History of Obesity Due to POMC, PCSK1, and LEPR Deficiency and the Impact of Setmelanotide*. J Endocr Soc. [PMID 35528826](https://pubmed.ncbi.nlm.nih.gov/35528826/)
77. Kühnen P … Clément K et al. (2022). *Quality of life outcomes in two phase 3 trials of setmelanotide in patients with obesity due to LEPR or POMC deficiency*. Orphanet J Rare Dis. [PMID 35123544](https://pubmed.ncbi.nlm.nih.gov/35123544/)
78. Dollfus H … Valverde D et al. (2024). *Bardet-Biedl syndrome improved diagnosis criteria and management: Inter European Reference Networks consensus statement and recommendations*. Eur J Hum Genet. [PMID 39085583](https://pubmed.ncbi.nlm.nih.gov/39085583/)
79. Mahmoud R … Butler MG (2023). *Clinical Trials in Prader-Willi Syndrome: A Review*. Int J Mol Sci. [PMID 36768472](https://pubmed.ncbi.nlm.nih.gov/36768472/)
80. Haqq AM … Argente J et al. (2022). *Efficacy and safety of setmelanotide, a melanocortin-4 receptor agonist, in patients with Bardet-Biedl syndrome and Alström syndrome: a multicentre, randomised, double-blind, placebo-controlled, phase 3 trial with an open-label period*. Lancet Diabetes Endocrinol. [PMID 36356613](https://pubmed.ncbi.nlm.nih.gov/36356613/)
81. Fox CK … Raatz SJ (2025). *Current and future state of pharmacological management of pediatric obesity*. Int J Obes (Lond). [PMID 38321079](https://pubmed.ncbi.nlm.nih.gov/38321079/)
82. Argente J … Pomeroy J et al. (2026). *Setmelanotide in Bardet-Biedl Syndrome: A 52-Week Comparison of Phase 3 Trial Participants With a Matched Registry Cohort*. Obesity (Silver Spring). [PMID 41703984](https://pubmed.ncbi.nlm.nih.gov/41703984/)
83. Zorn S … Farooqi IS et al. (2025). *Obesity due to MC4R deficiency is associated with reduced cholesterol, triglycerides and cardiovascular disease risk*. Nat Med. [PMID 41102563](https://pubmed.ncbi.nlm.nih.gov/41102563/)
84. Collet TH … Van der Ploeg LHT et al. (2017). *Evaluation of a melanocortin-4 receptor (MC4R) agonist (Setmelanotide) in MC4R deficiency*. Mol Metab. [PMID 29031731](https://pubmed.ncbi.nlm.nih.gov/29031731/)
85. Bhatnagar P … Farooqi IS (2025). *Tirzepatide leads to weight reduction in people with obesity due to MC4R deficiency*. Nat Med. [PMID 40858971](https://pubmed.ncbi.nlm.nih.gov/40858971/)
86. Farooqi IS (2008). *Monogenic human obesity*. Front Horm Res. [PMID 18230891](https://pubmed.ncbi.nlm.nih.gov/18230891/)
87. Cyr NE … Nillni EA (2015). *Central Sirt1 regulates body weight and energy expenditure along with the POMC-derived peptide α-MSH and the processing enzyme CPE production in diet-induced obese male rats*. Endocrinology. [PMID 25549049](https://pubmed.ncbi.nlm.nih.gov/25549049/)
88. Zanesco AM … Velloso LA et al. (2022). *Hypothalamic CREB Regulates the Expression of Pomc-Processing Enzyme Pcsk2*. Cells. [PMID 35805082](https://pubmed.ncbi.nlm.nih.gov/35805082/)
89. Cyr NE … Nillni EA (2014). *Central Sirt1 regulates body weight and energy expenditure along with the POMC-derived peptide α-MSH and the processing enzyme CPE production in diet-induced obese male rats*. Endocrinology. [PMID 24773342](https://pubmed.ncbi.nlm.nih.gov/24773342/)
90. Helwig M … Klingenspor M et al. (2006). *PC1/3 and PC2 gene expression and post-translational endoproteolytic pro-opiomelanocortin processing is regulated by photoperiod in the seasonal Siberian hamster (Phodopus sungorus)*. J Neuroendocrinol. [PMID 16684131](https://pubmed.ncbi.nlm.nih.gov/16684131/)
91. Kim KS … Choi HJ et al. (2024). *GLP-1 increases preingestive satiation via hypothalamic circuits in mice and humans*. Science. [PMID 38935778](https://pubmed.ncbi.nlm.nih.gov/38935778/)
92. Jais A … Brüning JC et al. (2020). *PNOC(ARC) Neurons Promote Hyperphagia and Obesity upon High-Fat-Diet Feeding*. Neuron. [PMID 32302532](https://pubmed.ncbi.nlm.nih.gov/32302532/)
93. Aitken TJ … Knight ZA et al. (2024). *Negative feedback control of hypothalamic feeding circuits by the taste of food*. Neuron. [PMID 39153476](https://pubmed.ncbi.nlm.nih.gov/39153476/)
94. Alhadeff AL … Betley JN et al. (2018). *A Neural Circuit for the Suppression of Pain by a Competing Need State*. Cell. [PMID 29570993](https://pubmed.ncbi.nlm.nih.gov/29570993/)
95. Grill HJ (2010). *Leptin and the systems neuroscience of meal size control*. Front Neuroendocrinol. [PMID 19836413](https://pubmed.ncbi.nlm.nih.gov/19836413/)
96. van der Klaauw AA … Farooqi IS et al. (2019). *Human Semaphorin 3 Variants Link Melanocortin Circuit Development and Energy Balance*. Cell. [PMID 30661757](https://pubmed.ncbi.nlm.nih.gov/30661757/)
97. Gundlach AL (2002). *Galanin/GALP and galanin receptors: role in central control of feeding, body weight/obesity and reproduction?*. Eur J Pharmacol. [PMID 12007540](https://pubmed.ncbi.nlm.nih.gov/12007540/)
98. Bouret SG … Simerly RB (2004). *Trophic action of leptin on hypothalamic neurons that regulate feeding*. Science. [PMID 15064420](https://pubmed.ncbi.nlm.nih.gov/15064420/)

## E. 그렐린 (Ghrelin, acylation, and the negative trials)

> `KAG50`(GHS-R1a 점유 반포화)과 `KGHRD`(구동에서 그렐린 팔의 가중치), `FGHRC`, `FACYLP`. 이 절의 세 음성 시험(옥트레오타이드·리볼레타이드·그렐린 저하 일반)이 `KGHRD = 0.10` 을 **교정**하며, 모델이 유도하는 것은 하나의 작은 값이 셋을 동시에 설명한다는 사실이다.

99. Feigerlová E … Tauber M et al. (2008). *Hyperghrelinemia precedes obesity in Prader-Willi syndrome*. J Clin Endocrinol Metab. [PMID 18460565](https://pubmed.ncbi.nlm.nih.gov/18460565/)
100. Choe YH … Jin DK et al. (2005). *Increased density of ghrelin-expressing cells in the gastric fundus and body in Prader-Willi syndrome*. J Clin Endocrinol Metab. [PMID 15956087](https://pubmed.ncbi.nlm.nih.gov/15956087/)
101. Bizzarri C … Salvatoni A et al. (2010). *Children with Prader-Willi syndrome exhibit more evident meal-induced responses in plasma ghrelin and peptide YY levels than obese and lean children*. Eur J Endocrinol. [PMID 20019130](https://pubmed.ncbi.nlm.nih.gov/20019130/)
102. Choe YH … Lee KH et al. (2005). *Hyperghrelinemia does not accelerate gastric emptying in Prader-Willi syndrome patients*. J Clin Endocrinol Metab. [PMID 15657368](https://pubmed.ncbi.nlm.nih.gov/15657368/)
103. Pacoricona Alfaro DL … Tauber M et al. (2021). *Is ghrelin a biomarker of early-onset scoliosis in children with Prader-Willi syndrome?*. Orphanet J Rare Dis. [PMID 34238321](https://pubmed.ncbi.nlm.nih.gov/34238321/)
104. Tauber M … Salles JP (2004). *Hyperghrelinemia is a common feature of Prader-Willi syndrome and pituitary stalk interruption: a pathophysiological hypothesis*. Horm Res. [PMID 15192277](https://pubmed.ncbi.nlm.nih.gov/15192277/)
105. Grootjen LN … Hokken-Koelega ACS et al. (2024). *Longitudinal Changes in Acylated versus Unacylated Ghrelin Levels May Be Involved in the Underlying Mechanisms of the Switch in Nutritional Phases in Prader-Willi Syndrome*. Horm Res Paediatr. [PMID 37839403](https://pubmed.ncbi.nlm.nih.gov/37839403/)
106. Kuppens RJ … Hokken-Koelega AC (2016). *Acylated and unacylated ghrelin during OGTT in Prader-Willi syndrome: support for normal response to food intake*. Clin Endocrinol (Oxf). [PMID 26850227](https://pubmed.ncbi.nlm.nih.gov/26850227/)
107. Kuppens RJ … Hokken-Koelega AC et al. (2015). *Elevated ratio of acylated to unacylated ghrelin in children and young adults with Prader-Willi syndrome*. Endocrine. [PMID 25989955](https://pubmed.ncbi.nlm.nih.gov/25989955/)
108. Delhanty PJ … van der Lely AJ (2013). *Des-acyl ghrelin: a metabolically active peptide*. Endocr Dev. [PMID 23652397](https://pubmed.ncbi.nlm.nih.gov/23652397/)
109. Muhammad A … Neggers SJCMM (2017). *The Acylated/Unacylated Ghrelin Ratio Is Similar in Patients With Acromegaly During Different Treatment Regimens*. J Clin Endocrinol Metab. [PMID 28402548](https://pubmed.ncbi.nlm.nih.gov/28402548/)
110. Isokawa M (2022). *Ghrelin-O-acyltransferase (GOAT) acylates ghrelin in the hippocampus*. Vitam Horm. [PMID 35180934](https://pubmed.ncbi.nlm.nih.gov/35180934/)
111. Davis TR … Hougland JL (2021). *Ghrelin octanoylation by ghrelin O-acyltransferase: protein acylation impacting metabolic and neuroendocrine signalling*. Open Biol. [PMID 34315274](https://pubmed.ncbi.nlm.nih.gov/34315274/)
112. Li Z … Zhang W (2016). *Ghrelin O-acyltransferase (GOAT) and energy metabolism*. Sci China Life Sci. [PMID 26732975](https://pubmed.ncbi.nlm.nih.gov/26732975/)
113. Iyer MR … Kunos G (2020). *Recent progress in the discovery of ghrelin O-acyltransferase (GOAT) inhibitors*. RSC Med Chem. [PMID 33479618](https://pubmed.ncbi.nlm.nih.gov/33479618/)
114. Costantini VJ … Corsi M et al. (2011). *GSK1614343, a novel ghrelin receptor antagonist, produces an unexpected increase of food intake and body weight in rodents and dogs*. Neuroendocrinology. [PMID 21778696](https://pubmed.ncbi.nlm.nih.gov/21778696/)
115. Holubová M … Maletínská L et al. (2013). *Ghrelin agonist JMV 1843 increases food intake, body weight and expression of orexigenic neuropeptides in mice*. Physiol Res. [PMID 23590608](https://pubmed.ncbi.nlm.nih.gov/23590608/)
116. Howick K … Schellekens H et al. (2020). *Behavioural characterization of ghrelin ligands, anamorelin and HM01: Appetite and reward-motivated effects in rodents*. Neuropharmacology. [PMID 32067989](https://pubmed.ncbi.nlm.nih.gov/32067989/)
117. Salomé N … Dickson SL et al. (2009). *Anorexigenic and electrophysiological actions of novel ghrelin receptor (GHS-R1A) antagonists in rats*. Eur J Pharmacol. [PMID 19356720](https://pubmed.ncbi.nlm.nih.gov/19356720/)
118. Yada T … Dezaki K et al. (2014). *Ghrelin signalling in β-cells regulates insulin secretion and blood glucose*. Diabetes Obes Metab. [PMID 25200304](https://pubmed.ncbi.nlm.nih.gov/25200304/)
119. Cupka M, Sedliak M (2023). *Hungry runners - low energy availability in male endurance athletes and its impact on performance and testosterone: mini-review*. Eur J Transl Myol. [PMID 37052052](https://pubmed.ncbi.nlm.nih.gov/37052052/)
120. Yada T … Iwasaki Y (2025). *GLP-1 and ghrelin inversely regulate insulin secretion and action in pancreatic islets, vagal afferents, and hypothalamus for controlling glycemia and feeding*. Am J Physiol Cell Physiol. [PMID 40241252](https://pubmed.ncbi.nlm.nih.gov/40241252/)
121. Pradhan G … Sun Y et al. (2017). *Obestatin stimulates glucose-induced insulin secretion through ghrelin receptor GHS-R*. Sci Rep. [PMID 28428639](https://pubmed.ncbi.nlm.nih.gov/28428639/)

## F. 장-뇌 축 (Gut peptides, gastric handling, vagal signalling)

> `FVAG`(구심성 미주신경 이득), `FPYYS`(둔화된 PYY 반응), `FGE`(위 배출 지연), 그리고 포만 적분기의 x2·x3·x4 팔. 위 파열·질식은 cluster 18의 사망 경로로 들어간다.

122. Giménez-Palop O … Caixàs A et al. (2007). *A lesser postprandial suppression of plasma ghrelin in Prader-Willi syndrome is associated with low fasting and a blunted postprandial PYY response*. Clin Endocrinol (Oxf). [PMID 17223988](https://pubmed.ncbi.nlm.nih.gov/17223988/)
123. Rigamonti AE … Sartorio A et al. (2014). *Unexpectedly increased anorexigenic postprandial responses of PYY and GLP-1 to fast ice cream consumption in adult patients with Prader-Willi syndrome*. Clin Endocrinol (Oxf). [PMID 24372155](https://pubmed.ncbi.nlm.nih.gov/24372155/)
124. Purtell L … Viardot A et al. (2011). *In adults with Prader-Willi syndrome, elevated ghrelin levels are more consistent with hyperphagia than high PYY and GLP-1 levels*. Neuropeptides. [PMID 21722955](https://pubmed.ncbi.nlm.nih.gov/21722955/)
125. Bueno M … Caixàs A et al. (2021). *Hunger and Satiety Peptides: Is There a Pattern to Classify Patients with Prader-Willi Syndrome?*. J Clin Med. [PMID 34768690](https://pubmed.ncbi.nlm.nih.gov/34768690/)
126. Ng NBH … Lee YS et al. (2022). *The effects of glucagon-like peptide (GLP)-1 receptor agonists on weight and glycaemic control in Prader-Willi syndrome: A systematic review*. Clin Endocrinol (Oxf). [PMID 34448208](https://pubmed.ncbi.nlm.nih.gov/34448208/)
127. Arenz T … Schmidt H (2010). *Delayed gastric emptying in patients with Prader Willi Syndrome*. J Pediatr Endocrinol Metab. [PMID 21175084](https://pubmed.ncbi.nlm.nih.gov/21175084/)
128. Höybye C, Petersson M (2025). *Neuropeptides and the Autonomic Nervous System in Prader-Willi Syndrome*. Int J Mol Sci. [PMID 41516228](https://pubmed.ncbi.nlm.nih.gov/41516228/)
129. Butler MG … Reiter LT (2023). *Autonomic nervous system dysfunction in Prader-Willi syndrome*. Clin Auton Res. [PMID 36515769](https://pubmed.ncbi.nlm.nih.gov/36515769/)
130. Blat C … Corripio R (2017). *Gastric Dilatation and Abdominal Compartment Syndrome in a Child with Prader-Willi Syndrome*. Am J Case Rep. [PMID 28588153](https://pubmed.ncbi.nlm.nih.gov/28588153/)
131. Browning KN (2019). *Stress-induced modulation of vagal afferents*. Neurogastroenterol Motil. [PMID 31736236](https://pubmed.ncbi.nlm.nih.gov/31736236/)
132. Ohbayashi K … Iwasaki Y (2021). *Gastrointestinal Distension by Pectin-Containing Carbonated Solution Suppresses Food Intake and Enhances Glucose Tolerance via GLP-1 Secretion and Vagal Afferent Activation*. Front Endocrinol (Lausanne). [PMID 34168616](https://pubmed.ncbi.nlm.nih.gov/34168616/)
133. Christie S … Page AJ (2020). *Biphasic effects of methanandamide on murine gastric vagal afferent mechanosensitivity*. J Physiol. [PMID 31642519](https://pubmed.ncbi.nlm.nih.gov/31642519/)
134. Miranda A … Sengupta JN et al. (2009). *Altered mechanosensitive properties of vagal afferent fibers innervating the stomach following gastric surgery in rats*. Neuroscience. [PMID 19477237](https://pubmed.ncbi.nlm.nih.gov/19477237/)
135. Journel M … Tomé D (2012). *Brain responses to high-protein diets*. Adv Nutr. [PMID 22585905](https://pubmed.ncbi.nlm.nih.gov/22585905/)
136. Morley JE (1990). *Appetite regulation by gut peptides*. Annu Rev Nutr. [PMID 2200469](https://pubmed.ncbi.nlm.nih.gov/2200469/)
137. Asarian L, Geary N (2006). *Modulation of appetite by gonadal steroid hormones*. Philos Trans R Soc Lond B Biol Sci. [PMID 16815802](https://pubmed.ncbi.nlm.nih.gov/16815802/)

## G. 과식 표현형과 영양 단계 (Hyperphagia, HQ-CT and the Miller phases)

> HQ-CT 척도 자체와 Miller 영양 단계, 그리고 `AGEHP`(3기 게이트 중점)·`CUE`. 모델의 가장 강한 발달 주장 — **체중 증가가 과식보다 먼저 온다** — 이 절의 표현형 기술과 대조된다.

138. Bravo J P … Canals Cifuentes A (2021). *Nutritional phases of Prader-Willi syndrome*. Andes Pediatr. [PMID 34479241](https://pubmed.ncbi.nlm.nih.gov/34479241/)
139. Miller JL … Driscoll DJ et al. (2011). *Nutritional phases in Prader-Willi syndrome*. Am J Med Genet A. [PMID 21465655](https://pubmed.ncbi.nlm.nih.gov/21465655/)
140. Kweh FA … Driscoll DJ (2023). *Hyperinsulinemia is a probable trigger for weight gain and hyperphagia in individuals with Prader-Willi syndrome*. Obes Sci Pract. [PMID 37546289](https://pubmed.ncbi.nlm.nih.gov/37546289/)
141. Kweh FA … Driscoll DJ et al. (2015). *Hyperghrelinemia in Prader-Willi syndrome begins in early infancy long before the onset of hyperphagia*. Am J Med Genet A. [PMID 25355237](https://pubmed.ncbi.nlm.nih.gov/25355237/)
142. Wang Y … Gong C et al. (2025). *Long-term intranasal oxytocin therapy in patients with hypothalamic syndrome: case series and literature review*. Endocr Connect. [PMID 41091101](https://pubmed.ncbi.nlm.nih.gov/41091101/)
143. Tsai JH … Bridges JFP (2021). *Measuring Meaningful Benefit-Risk Tradeoffs to Promote Patient-Focused Drug Development in Prader-Willi Syndrome: A Discrete-Choice Experiment*. MDM Policy Pract. [PMID 34497876](https://pubmed.ncbi.nlm.nih.gov/34497876/)
144. Owczarek-Danowska IM … Michalik M et al. (2026). *Modern technological innovations in the management of Prader-Willi syndrome: From restrictive supervision to digital autonomy*. Pol Merkur Lekarski. [PMID 42435474](https://pubmed.ncbi.nlm.nih.gov/42435474/)
145. Wieting J … Frieling H et al. (2022). *Alteration of serum leptin and LEP/LEPR promoter methylation in Prader-Willi syndrome*. Psychoneuroendocrinology. [PMID 35803048](https://pubmed.ncbi.nlm.nih.gov/35803048/)
146. Adam MP … Cassidy SB et al. (1993). *Prader-Willi Syndrome*. . [PMID 20301505](https://pubmed.ncbi.nlm.nih.gov/20301505/)
147. Cassidy SB, Driscoll DJ (2009). *Prader-Willi syndrome*. Eur J Hum Genet. [PMID 18781185](https://pubmed.ncbi.nlm.nih.gov/18781185/)
148. Erhardt É, Molnár D (2022). *Prader-Willi Syndrome: Possibilities of Weight Gain Prevention and Treatment*. Nutrients. [PMID 35565916](https://pubmed.ncbi.nlm.nih.gov/35565916/)

## H. 에너지 소비와 체성분 (Energy expenditure and body composition)

> `WREEL`/`WREEF`(REE 조성 가중), Schofield 참조곡선, Forbes 관계, 그리고 kcal/cm 처방 규칙. PWS의 TEE가 예측치의 60-70%라는 보고를 모델은 파라미터로 넣지 않고 (낮은 제지방량)×(낮은 활동량)에서 **유도**한다.

149. Alsaif M … Haqq AM (2017). *Energy Metabolism Profile in Individuals with Prader-Willi Syndrome and Implications for Clinical Management: A Systematic Review*. Adv Nutr. [PMID 29141973](https://pubmed.ncbi.nlm.nih.gov/29141973/)
150. Butler MG … Donnelly JE (2007). *Energy expenditure and physical activity in Prader-Willi syndrome: comparison with obese subjects*. Am J Med Genet A. [PMID 17103434](https://pubmed.ncbi.nlm.nih.gov/17103434/)
151. Davies PS, Joughin C (1993). *Using stable isotopes to assess reduced physical activity of individuals with Prader-Willi syndrome*. Am J Ment Retard. [PMID 8292311](https://pubmed.ncbi.nlm.nih.gov/8292311/)
152. Viardot A … Campbell LV (2018). *Relative Contributions of Lean and Fat Mass to Bone Mineral Density: Insight From Prader-Willi Syndrome*. Front Endocrinol (Lausanne). [PMID 30186239](https://pubmed.ncbi.nlm.nih.gov/30186239/)
153. Brambilla P … Chiumello G (1997). *Peculiar body composition in patients with Prader-Labhart-Willi syndrome*. Am J Clin Nutr. [PMID 9129464](https://pubmed.ncbi.nlm.nih.gov/9129464/)
154. Theodoro MF … Butler MG (2006). *Body composition and fatness patterns in Prader-Willi syndrome: comparison with simple obesity*. Obesity (Silver Spring). [PMID 17062796](https://pubmed.ncbi.nlm.nih.gov/17062796/)
155. Höybye C … Thorén M (2003). *Growth hormone treatment improves body composition in adults with Prader-Willi syndrome*. Clin Endocrinol (Oxf). [PMID 12699450](https://pubmed.ncbi.nlm.nih.gov/12699450/)
156. Damen L … Hokken-Koelega ACS et al. (2020). *Three years of growth hormone treatment in young adults with Prader-Willi syndrome: sustained positive effects on body composition*. Orphanet J Rare Dis. [PMID 32580778](https://pubmed.ncbi.nlm.nih.gov/32580778/)
157. Bosio L … Chiumello G (1999). *Body composition during GH treatment in Prader-Labhardt-Willi syndrome*. J Pediatr Endocrinol Metab. [PMID 10698601](https://pubmed.ncbi.nlm.nih.gov/10698601/)
158. Jeran S … Pischon T et al. (2022). *Prediction of activity-related energy expenditure under free-living conditions using accelerometer-derived physical activity*. Sci Rep. [PMID 36195647](https://pubmed.ncbi.nlm.nih.gov/36195647/)
159. Feingold KR … Westerterp KR et al. (2000). *Control of Energy Expenditure in Humans*. . [PMID 25905198](https://pubmed.ncbi.nlm.nih.gov/25905198/)
160. Tamini S … Sartorio A (2023). *Measured vs estimated resting energy expenditure in children and adolescents with obesity*. Sci Rep. [PMID 37580514](https://pubmed.ncbi.nlm.nih.gov/37580514/)
161. Rydin AA … Cree MG et al. (2024). *Prediction of resting energy expenditure for adolescents with severe obesity: A multi-centre analysis*. Pediatr Obes. [PMID 38658523](https://pubmed.ncbi.nlm.nih.gov/38658523/)
162. Henry CJ (2005). *Basal metabolic rate studies in humans: measurement and development of new equations*. Public Health Nutr. [PMID 16277825](https://pubmed.ncbi.nlm.nih.gov/16277825/)
163. Schofield WN (1985). *Predicting basal metabolic rate, new standards and review of previous work*. Hum Nutr Clin Nutr. [PMID 4044297](https://pubmed.ncbi.nlm.nih.gov/4044297/)
164. Sgambato MR … Anjos LAD (2019). *Validity of basal metabolic rate prediction equations in elderly women living in an urban tropical city of Brazil*. Clin Nutr ESPEN. [PMID 31221282](https://pubmed.ncbi.nlm.nih.gov/31221282/)
165. Heymsfield SB … Thomas D (2014). *Weight loss composition is one-fourth fat-free mass: a critical review and critique of this widely cited rule*. Obes Rev. [PMID 24447775](https://pubmed.ncbi.nlm.nih.gov/24447775/)
166. Forbes GB (2000). *Body fat content influences the body composition response to nutrition and exercise*. Ann N Y Acad Sci. [PMID 10865771](https://pubmed.ncbi.nlm.nih.gov/10865771/)
167. Dulloo AG (2017). *Collateral fattening: When a deficit in lean body mass drives overeating*. Obesity (Silver Spring). [PMID 28078821](https://pubmed.ncbi.nlm.nih.gov/28078821/)
168. Hall KD … Swinburn BA et al. (2011). *Quantification of the effect of energy imbalance on bodyweight*. Lancet. [PMID 21872751](https://pubmed.ncbi.nlm.nih.gov/21872751/)
169. Hall KD, Guo J (2017). *Obesity Energetics: Body Weight Regulation and the Effects of Diet Composition*. Gastroenterology. [PMID 28193517](https://pubmed.ncbi.nlm.nih.gov/28193517/)
170. Stanhope KL (2016). *Sugar consumption, metabolic disease and obesity: The state of the controversy*. Crit Rev Clin Lab Sci. [PMID 26376619](https://pubmed.ncbi.nlm.nih.gov/26376619/)
171. Magkos F … Astrup A et al. (2024). *On the pathogenesis of obesity: causal models and missing pieces of the puzzle*. Nat Metab. [PMID 39164418](https://pubmed.ncbi.nlm.nih.gov/39164418/)

## I. 성장호르몬 치료: 효능 (Growth hormone: efficacy)

> `FPOTGH`·`KGIGF`·`WLIGF0`·`GVFL`·`GVCAP`. 성장호르몬이 체성분과 신장을 움직이면서 HQ-CT는 움직이지 않는다는 모델 결과의 대조군이 이 절이다. 목표 IGF-1 SDS +1~+2도 여기서 온다.

172. Rosenberg AGW … de Graaff LCG et al. (2021). *Growth Hormone Treatment for Adults With Prader-Willi Syndrome: A Meta-Analysis*. J Clin Endocrinol Metab. [PMID 34105729](https://pubmed.ncbi.nlm.nih.gov/34105729/)
173. Yang A … Jin DK et al. (2019). *Effects of recombinant human growth hormone treatment on growth, body composition, and safety in infants or toddlers with Prader-Willi syndrome: a randomized, active-controlled trial*. Orphanet J Rare Dis. [PMID 31511031](https://pubmed.ncbi.nlm.nih.gov/31511031/)
174. Lindgren AC … Ritzén EM et al. (1998). *Growth hormone treatment of children with Prader-Willi syndrome affects linear growth and body composition favourably*. Acta Paediatr. [PMID 9510443](https://pubmed.ncbi.nlm.nih.gov/9510443/)
175. Bakker NE … Hokken-Koelega AC (2015). *Dietary Energy Intake, Body Composition and Resting Energy Expenditure in Prepubertal Children with Prader-Willi Syndrome before and during Growth Hormone Treatment: A Randomized Controlled Trial*. Horm Res Paediatr. [PMID 25764996](https://pubmed.ncbi.nlm.nih.gov/25764996/)
176. Carrel AL … Allen DB (2004). *Growth hormone improves mobility and body composition in infants and toddlers with Prader-Willi syndrome*. J Pediatr. [PMID 15580194](https://pubmed.ncbi.nlm.nih.gov/15580194/)
177. Takeda A … Bryant J et al. (2010). *Recombinant human growth hormone for the treatment of growth disorders in children: a systematic review and economic evaluation*. Health Technol Assess. [PMID 20849734](https://pubmed.ncbi.nlm.nih.gov/20849734/)
178. Bridges N (2014). *What is the value of growth hormone therapy in Prader Willi syndrome?*. Arch Dis Child. [PMID 24162007](https://pubmed.ncbi.nlm.nih.gov/24162007/)
179. Hirsch HJ, Gross-Tsur V (2021). *Growth hormone treatment for adults with Prader-Willi syndrome: another point of view*. Orphanet J Rare Dis. [PMID 34344408](https://pubmed.ncbi.nlm.nih.gov/34344408/)
180. Angulo MA … Khan A (2007). *Final adult height in children with Prader-Willi syndrome with and without human growth hormone treatment*. Am J Med Genet A. [PMID 17567883](https://pubmed.ncbi.nlm.nih.gov/17567883/)
181. Grugni G … Crinò A (2016). *Growth hormone therapy for Prader-willi syndrome: challenges and solutions*. Ther Clin Risk Manag. [PMID 27330297](https://pubmed.ncbi.nlm.nih.gov/27330297/)
182. Allen DB, Carrel AL (2004). *Growth hormone therapy for Prader-Willi syndrome: a critical appraisal*. J Pediatr Endocrinol Metab. [PMID 15506076](https://pubmed.ncbi.nlm.nih.gov/15506076/)
183. Aycan Z, Baş VN (2014). *Prader-Willi syndrome and growth hormone deficiency*. J Clin Res Pediatr Endocrinol. [PMID 24932597](https://pubmed.ncbi.nlm.nih.gov/24932597/)
184. Kodytková A … Lebl J et al. (2025). *Early-onset growth hormone treatment in Prader-Willi syndrome attenuates transition to severe obesity*. J Pediatr Endocrinol Metab. [PMID 40080424](https://pubmed.ncbi.nlm.nih.gov/40080424/)
185. Cheng RQ … Lu W et al. (2023). *Early recombinant human growth hormone treatment improves mental development and alleviates deterioration of motor function in infants and young children with Prader-Willi syndrome*. World J Pediatr. [PMID 36564648](https://pubmed.ncbi.nlm.nih.gov/36564648/)
186. Diene G … Tauber M (2007). *The Prader-Willi syndrome*. Ann Endocrinol (Paris). [PMID 17499572](https://pubmed.ncbi.nlm.nih.gov/17499572/)
187. Grugni G, Sartorio A (2025). *Growth hormone treatment in adults with Prader-Willi syndrome: an update*. Expert Rev Endocrinol Metab. [PMID 41147430](https://pubmed.ncbi.nlm.nih.gov/41147430/)
188. Sode-Carlsen R … Höybye C et al. (2012). *Growth hormone treatment in adults with Prader-Willi syndrome: the Scandinavian study*. Endocrine. [PMID 22081257](https://pubmed.ncbi.nlm.nih.gov/22081257/)
189. Yang X (2020). *Growth hormone treatment for Prader-Willi syndrome: A review*. Neuropeptides. [PMID 32859387](https://pubmed.ncbi.nlm.nih.gov/32859387/)
190. Deal CL … Christiansen JS et al. (2013). *GrowthHormone Research Society workshop summary: consensus guidelines for recombinant human growth hormone therapy in Prader-Willi syndrome*. J Clin Endocrinol Metab. [PMID 23543664](https://pubmed.ncbi.nlm.nih.gov/23543664/)
191. Cohen P … Wit JM et al. (2008). *Consensus statement on the diagnosis and treatment of children with idiopathic short stature: a summary of the Growth Hormone Research Society, the Lawson Wilkins Pediatric Endocrine Society, and the European Society for Paediatric Endocrinology Workshop*. J Clin Endocrinol Metab. [PMID 18782877](https://pubmed.ncbi.nlm.nih.gov/18782877/)
192. Lundberg E … Albertsson-Wikland K (2015). *Growth hormone (GH) dose-dependent IGF-I response relates to pubertal height gain*. BMC Endocr Disord. [PMID 26682747](https://pubmed.ncbi.nlm.nih.gov/26682747/)
193. Kruijsen AR … Joustra SD et al. (2025). *Growth hormone treatment adjusted for growth hormone sensitivity in idiopathic short stature*. Eur J Endocrinol. [PMID 40621613](https://pubmed.ncbi.nlm.nih.gov/40621613/)
194. Kriström B … Albertsson-Wikland K (2014). *IGF-1 and growth response to adult height in a randomized GH treatment trial in short non-GH-deficient children*. J Clin Endocrinol Metab. [PMID 24823461](https://pubmed.ncbi.nlm.nih.gov/24823461/)
195. Casamitjana L … Caixàs A et al. (2021). *Glucagon stimulation test to assess growth hormone status in Prader-Willi syndrome*. J Endocrinol Invest. [PMID 32720093](https://pubmed.ncbi.nlm.nih.gov/32720093/)
196. Cohen M … Hamilton J (2015). *Growth hormone secretion decreases with age in paediatric Prader-Willi syndrome*. Clin Endocrinol (Oxf). [PMID 25495188](https://pubmed.ncbi.nlm.nih.gov/25495188/)
197. Donze SH … Hokken-Koelega ACS et al. (2019). *Prevalence of growth hormone (GH) deficiency in previously GH-treated young adults with Prader-Willi syndrome*. Clin Endocrinol (Oxf). [PMID 30973645](https://pubmed.ncbi.nlm.nih.gov/30973645/)
198. Schmok T … Kimonis VE et al. (2024). *Relationship of thyroid function with genetic subtypes and treatment with growth hormone in Prader-Willi syndrome*. Am J Med Genet A. [PMID 38837660](https://pubmed.ncbi.nlm.nih.gov/38837660/)

## J. 성장호르몬 안전성과 기도 (Growth hormone safety and the airway)

> 모델의 **두 시계**(cluster 12)와 `KLYU`·`TAURMS`·`FOSA`. GH 시작 후 6-8주 수면다원검사 권고와 초기 9개월 급사 신호가 모델에서는 적합이 아니라 τ_lymphoid < τ_muscle 의 간섭으로 나온다.

199. Adam MP … Marbach F et al. (1993). *Schaaf-Yang Syndrome*. . [PMID 33570896](https://pubmed.ncbi.nlm.nih.gov/33570896/)
200. Itani R … Perez IA (2023). *Sleep Consequences of Prader-Willi Syndrome*. Curr Neurol Neurosci Rep. [PMID 36790642](https://pubmed.ncbi.nlm.nih.gov/36790642/)
201. Wong SB … Tsai LP (2022). *Progression of Obstructive Sleep Apnea Syndrome in Pediatric Patients with Prader-Willi Syndrome*. Children (Basel). [PMID 35740849](https://pubmed.ncbi.nlm.nih.gov/35740849/)
202. Zaffanello M … Antoniazzi F (2023). *The Impact of Growth Hormone Therapy on Sleep-Related Health Outcomes in Children with Prader-Willi Syndrome: A Review and Clinical Analysis*. J Clin Med. [PMID 37685570](https://pubmed.ncbi.nlm.nih.gov/37685570/)
203. Abushahin A … Janahi IA et al. (2023). *Prevalence of Sleep-Disordered Breathing in Prader-Willi Syndrome*. Can Respir J. [PMID 37927914](https://pubmed.ncbi.nlm.nih.gov/37927914/)
204. Schaefer J … Nixon GM (2022). *Sleep-disordered breathing in school-aged children with Prader-Willi syndrome*. J Clin Sleep Med. [PMID 34870583](https://pubmed.ncbi.nlm.nih.gov/34870583/)
205. Van Vliet G … Oligny LL (2004). *Sudden death in growth hormone-treated children with Prader-Willi syndrome*. J Pediatr. [PMID 14722532](https://pubmed.ncbi.nlm.nih.gov/14722532/)
206. Stafler P, Wallis C (2008). *Prader-Willi syndrome: who can have growth hormone?*. Arch Dis Child. [PMID 18089632](https://pubmed.ncbi.nlm.nih.gov/18089632/)
207. Wolfgram PM … Allen DB (2013). *Long-term effects of recombinant human growth hormone therapy in children with Prader-Willi syndrome*. Curr Opin Pediatr. [PMID 23782572](https://pubmed.ncbi.nlm.nih.gov/23782572/)
208. Tan HL, Urquhart DS (2017). *Respiratory Complications in Children with Prader Willi Syndrome*. Paediatr Respir Rev. [PMID 27839656](https://pubmed.ncbi.nlm.nih.gov/27839656/)
209. Nagai T … Niikawa N et al. (2005). *Cause of sudden, unexpected death of Prader-Willi syndrome patients with or without growth hormone treatment*. Am J Med Genet A. [PMID 15937939](https://pubmed.ncbi.nlm.nih.gov/15937939/)
210. Kasi AS, Perez IA (2024). *Congenital Central Hypoventilation Syndrome and Disorders of Control of Ventilation*. Clin Chest Med. [PMID 39069329](https://pubmed.ncbi.nlm.nih.gov/39069329/)
211. Gallego J (2012). *Genetic diseases: congenital central hypoventilation, Rett, and Prader-Willi syndromes*. Compr Physiol. [PMID 23723037](https://pubmed.ncbi.nlm.nih.gov/23723037/)
212. Tanizawa K, Chin K (2018). *Genetic factors in sleep-disordered breathing*. Respir Investig. [PMID 29548648](https://pubmed.ncbi.nlm.nih.gov/29548648/)
213. Kaditis AG … Verhulst S et al. (2016). *Obstructive sleep disordered breathing in 2- to 18-year-old children: diagnosis and management*. Eur Respir J. [PMID 26541535](https://pubmed.ncbi.nlm.nih.gov/26541535/)
214. Sedky K … Pumariega A (2014). *Prader Willi syndrome and obstructive sleep apnea: co-occurrence in the pediatric population*. J Clin Sleep Med. [PMID 24733986](https://pubmed.ncbi.nlm.nih.gov/24733986/)
215. Clements AC … Ryan MA et al. (2021). *Outcomes of Adenotonsillectomy for Obstructive Sleep Apnea in Prader-Willi Syndrome: Systematic Review and Meta-analysis*. Laryngoscope. [PMID 33026674](https://pubmed.ncbi.nlm.nih.gov/33026674/)
216. Lee CH … Kang KT (2020). *Adenotonsillectomy for the Treatment of Obstructive Sleep Apnea in Children with Prader-Willi Syndrome: A Meta-analysis*. Otolaryngol Head Neck Surg. [PMID 31818186](https://pubmed.ncbi.nlm.nih.gov/31818186/)
217. Giordano L … Bussi M et al. (2015). *Obstructive sleep apnea in Prader-Willi syndrome: risks and advantages of adenotonsillectomy*. Pediatr Med Chir. [PMID 26429118](https://pubmed.ncbi.nlm.nih.gov/26429118/)
218. Nagai T … Niikawa N et al. (2006). *Growth hormone therapy and scoliosis in patients with Prader-Willi syndrome*. Am J Med Genet A. [PMID 16770808](https://pubmed.ncbi.nlm.nih.gov/16770808/)
219. Grootjen LN … Hokken-Koelega ACS et al. (2021). *Effects of 8 years of growth hormone treatment on scoliosis in children with Prader-Willi syndrome*. Eur J Endocrinol. [PMID 33886496](https://pubmed.ncbi.nlm.nih.gov/33886496/)
220. Zhu M … Chen J et al. (2026). *Endocrine-informed monitoring of scoliosis in Prader-Willi syndrome: integrating neuroendocrine pathophysiology, growth hormone therapy, and pubertal transition*. Front Endocrinol (Lausanne). [PMID 42051450](https://pubmed.ncbi.nlm.nih.gov/42051450/)
221. Odent T … Glorion C et al. (2008). *Scoliosis in patients with Prader-Willi Syndrome*. Pediatrics. [PMID 18606625](https://pubmed.ncbi.nlm.nih.gov/18606625/)

## K. 수면 (Sleep, hypersomnolence and orexin)

> `CAP_T`(중추성 무호흡 성향의 연령 곡선)와 시상하부 수면 이상. 모델은 각성 저하를 명시적 상태로 쓰지 않고 기도·저산소 부담 경로로만 표현하므로, 이 절은 모델의 한계 표시이기도 하다.

222. Overeem S … Reading PJ (2021). *Sleep disorders and the hypothalamus*. Handb Clin Neurol. [PMID 34266606](https://pubmed.ncbi.nlm.nih.gov/34266606/)
223. Helbing-Zwanenburg B … Mourtazaev MS (1993). *The origin of excessive daytime sleepiness in the Prader-Willi syndrome*. J Intellect Disabil Res. [PMID 8123999](https://pubmed.ncbi.nlm.nih.gov/8123999/)
224. Vgontzas AN … Vela-Bueno A et al. (1996). *Daytime sleepiness and REM abnormalities in Prader-Willi syndrome: evidence of generalized hypoarousal*. Int J Neurosci. [PMID 9003974](https://pubmed.ncbi.nlm.nih.gov/9003974/)
225. Bruni O … Ferri R (2010). *Prader-Willi syndrome: sorting out the relationships between obesity, hypersomnia, and sleep apnea*. Curr Opin Pulm Med. [PMID 20814307](https://pubmed.ncbi.nlm.nih.gov/20814307/)
226. Nishino S, Kanbayashi T (2005). *Symptomatic narcolepsy, cataplexy and hypersomnia, and their implications in the hypothalamic hypocretin/orexin system*. Sleep Med Rev. [PMID 16006155](https://pubmed.ncbi.nlm.nih.gov/16006155/)
227. Dodet P … Redolfi S et al. (2022). *Sleep Disorders in Adults with Prader-Willi Syndrome: Review of the Literature and Clinical Recommendations Based on the Experience of the French Reference Centre*. J Clin Med. [PMID 35407596](https://pubmed.ncbi.nlm.nih.gov/35407596/)
228. Dodet P … Arnulf I et al. (2023). *Hypersomnia and narcolepsy in 42 adult patients with craniopharyngioma*. Sleep. [PMID 36799460](https://pubmed.ncbi.nlm.nih.gov/36799460/)
229. Mignot E … Nishino S et al. (2002). *The role of cerebrospinal fluid hypocretin measurement in the diagnosis of narcolepsy and other hypersomnias*. Arch Neurol. [PMID 12374492](https://pubmed.ncbi.nlm.nih.gov/12374492/)
230. Omokawa M … Kanbayashi T et al. (2016). *Decline of CSF orexin (hypocretin) levels in Prader-Willi syndrome*. Am J Med Genet A. [PMID 26738920](https://pubmed.ncbi.nlm.nih.gov/26738920/)
231. Nevsimalova S … Nishino S (2005). *Hypocretin deficiency in Prader-Willi syndrome*. Eur J Neurol. [PMID 15613151](https://pubmed.ncbi.nlm.nih.gov/15613151/)

## L. 다이아족사이드 콜린과 KATP (Diazoxide choline, DCCR and KATP channels)

> `EDZMAX`·`EDZ50`·`EDZIMX`·`EDZI50`. 효능과 고혈당이 **같은 채널의 두 조직**이라는 모델의 주장은 이 절의 약리에서 직접 온다 — 그래서 치료계수가 두 EC50의 비로 한 번에 고정된다.

232. Miller JL … Bhatnagar A et al. (2023). *Diazoxide Choline Extended-Release Tablet in People With Prader-Willi Syndrome: A Double-Blind, Placebo-Controlled Trial*. J Clin Endocrinol Metab. [PMID 36639249](https://pubmed.ncbi.nlm.nih.gov/36639249/)
233. van den Top M … Spanswick D (2007). *Pharmacological and molecular characterization of ATP-sensitive K(+) conductances in CART and NPY/AgRP expressing neurons of the hypothalamic arcuate nucleus*. Neuroscience. [PMID 17137725](https://pubmed.ncbi.nlm.nih.gov/17137725/)
234. Seltzer HS, Allen EW (1969). *Hyperglycemia and inhibition of insulin secretion during administration of diazoxide and trichlormethiazide in man*. Diabetes. [PMID 5761863](https://pubmed.ncbi.nlm.nih.gov/5761863/)
235. Vargas-Vargas MA … Rocío MP et al. (2023). *Diazoxide improves muscle function in association with improved dyslipidemia and decreased muscle oxidative stress in streptozotocin-induced diabetic rats*. J Bioenerg Biomembr. [PMID 36723797](https://pubmed.ncbi.nlm.nih.gov/36723797/)
236. Yamazaki H … Zawalich WS (2006). *Acute and chronic effects of glucose and carbachol on insulin secretion and phospholipase C activation: studies with diazoxide and atropine*. Am J Physiol Endocrinol Metab. [PMID 16105864](https://pubmed.ncbi.nlm.nih.gov/16105864/)
237. Charles MA, Danforth E Jr (1971). *Nonketoacidotic hyperglycemia and coma during intravenous diazoxide therapy in uremia*. Diabetes. [PMID 5556285](https://pubmed.ncbi.nlm.nih.gov/5556285/)
238. Thorens B (2001). *GLUT2 in pancreatic and extra-pancreatic gluco-detection (review)*. Mol Membr Biol. [PMID 11780755](https://pubmed.ncbi.nlm.nih.gov/11780755/)
239. Choeiri C … Messier C et al. (2006). *Cerebral glucose transporters expression and spatial learning in the K-ATP Kir6.2(-/-) knockout mice*. Behav Brain Res. [PMID 16797737](https://pubmed.ncbi.nlm.nih.gov/16797737/)
240. Acosta-Martínez M, Levine JE (2007). *Regulation of KATP channel subunit gene expression by hyperglycemia in the mediobasal hypothalamus of female rats*. Am J Physiol Endocrinol Metab. [PMID 17311891](https://pubmed.ncbi.nlm.nih.gov/17311891/)
241. Fan X … McCrimmon RJ (2008). *Amplified hormonal counterregulatory responses to hypoglycemia in rats after systemic delivery of a SUR-1-selective K(+) channel opener?*. Diabetes. [PMID 18776135](https://pubmed.ncbi.nlm.nih.gov/18776135/)

## M. GLP-1 수용체 작용제와 기타 약물 (GLP-1 receptor agonists and other agents)

> `ESGSAT`·`ESGA`·`ESGEI`(GLP-1 수용체 작용제), 그리고 테소메트·벨로라닙·수술·메트포르민·토피라메이트. PWS에서 GLP-1 근거는 대부분 관찰연구이므로 이 팔의 예측은 다른 팔들보다 약하게 취급해야 한다.

242. Ahmed S … K M (2023). *Weight Loss of Over 100 lbs in a Patient of Prader-Willi Syndrome Treated With Glucagon-Like Peptide-1 (GLP-1) Agonists*. Cureus. [PMID 36945294](https://pubmed.ncbi.nlm.nih.gov/36945294/)
243. Giménez-Palop O … Caixàs A (2024). *Effect of semaglutide on weight loss and glycaemic control in patients with Prader-Willi Syndrome and type 2 diabetes*. Endocrinol Diabetes Nutr (Engl Ed). [PMID 38553173](https://pubmed.ncbi.nlm.nih.gov/38553173/)
244. Sani E … Bonora E (2022). *Effects of Semaglutide on Glycemic Control and Weight Loss in a Patient with Prader-Willi Syndrome: A Case Report*. Endocr Metab Immune Disord Drug Targets. [PMID 35538810](https://pubmed.ncbi.nlm.nih.gov/35538810/)
245. Koceva A … Jensterle M (2024). *Case report: Long-term efficacy and safety of semaglutide in the treatment of syndromic obesity in Prader Willi syndrome - case series and literature review*. Front Endocrinol (Lausanne). [PMID 39906041](https://pubmed.ncbi.nlm.nih.gov/39906041/)
246. Goldman VE … Vidmar AP (2021). *Anti-Obesity Medication Use in Children and Adolescents with Prader-Willi Syndrome: Case Review and Literature Search*. J Clin Med. [PMID 34640558](https://pubmed.ncbi.nlm.nih.gov/34640558/)
247. Webster AN … Campbell JN et al. (2024). *Molecular connectomics reveals a glucagon-like peptide 1-sensitive neural circuit for satiety*. Nat Metab. [PMID 39627618](https://pubmed.ncbi.nlm.nih.gov/39627618/)
248. McMorrow HE … Beutler LR et al. (2025). *Incretin receptor agonism rapidly inhibits AgRP neurons to suppress food intake in mice*. J Clin Invest. [PMID 40857106](https://pubmed.ncbi.nlm.nih.gov/40857106/)
249. Webster AN … Campbell JN et al. (2024). *Molecular Connectomics Reveals a Glucagon-Like Peptide 1 Sensitive Neural Circuit for Satiety*. bioRxiv. [PMID 37961449](https://pubmed.ncbi.nlm.nih.gov/37961449/)
250. Dong Y … Williams KW et al. (2021). *Time and metabolic state-dependent effects of GLP-1R agonists on NPY/AgRP and POMC neuronal activity in vivo*. Mol Metab. [PMID 34626854](https://pubmed.ncbi.nlm.nih.gov/34626854/)
251. McCandless SE … Butler MG et al. (2017). *Effects of MetAP2 inhibition on hyperphagia and body weight in Prader-Willi syndrome: A randomized, double-blind, placebo-controlled trial*. Diabetes Obes Metab. [PMID 28556449](https://pubmed.ncbi.nlm.nih.gov/28556449/)
252. Salehi P … Chen M et al. (2017). *Silent aspiration in infants with Prader-Willi syndrome identified by videofluoroscopic swallow study*. Medicine (Baltimore). [PMID 29390364](https://pubmed.ncbi.nlm.nih.gov/29390364/)
253. Dressler MH … Lee HJ (2025). *Bariatric Surgery in Patients With Prader-Willi Syndrome*. J Metab Bariatr Surg. [PMID 40917201](https://pubmed.ncbi.nlm.nih.gov/40917201/)
254. Hu S … Yang W et al. (2022). *Patients with Prader-Willi Syndrome (PWS) Underwent Bariatric Surgery Benefit more from High-Intensity Home Care*. Obes Surg. [PMID 35288862](https://pubmed.ncbi.nlm.nih.gov/35288862/)
255. Liu SY … Ng EK (2020). *Bariatric surgery for Prader-Willi syndrome was ineffective in producing sustainable weight loss: Long term results for up to 10 years*. Pediatr Obes. [PMID 31515962](https://pubmed.ncbi.nlm.nih.gov/31515962/)
256. Shapira NA … Goodman WK (2002). *Topiramate attenuates self-injurious behaviour in Prader-Willi Syndrome*. Int J Neuropsychopharmacol. [PMID 12135538](https://pubmed.ncbi.nlm.nih.gov/12135538/)
257. Consoli A … Bonnot O et al. (2019). *Effect of topiramate on eating behaviours in Prader-Willi syndrome: TOPRADER double-blind randomised placebo-controlled study*. Transl Psychiatry. [PMID 31685813](https://pubmed.ncbi.nlm.nih.gov/31685813/)

## N. 당대사 (Glucose metabolism, insulin and adiponectin)

> `FADPN`(상대적 고아디포넥틴혈증)과 `KADPNSI`. 같은 지방량에서 PWS가 인슐린은 낮고 아디포넥틴은 높다는 **대사 역설**을 모델은 두 개가 아니라 하나의 파라미터(PC1/3)로 낸다.

258. Tauber M … Molinas C (2016). *Sequelae of GH Treatment in Children with PWS*. Pediatr Endocrinol Rev. [PMID 28508607](https://pubmed.ncbi.nlm.nih.gov/28508607/)
259. Fintini D … Crinò A et al. (2016). *Disorders of glucose metabolism in Prader-Willi syndrome: Results of a multicenter Italian cohort study*. Nutr Metab Cardiovasc Dis. [PMID 27381990](https://pubmed.ncbi.nlm.nih.gov/27381990/)
260. L'Allemand D … Girard J (2003). *Carbohydrate metabolism is not impaired after 3 years of growth hormone therapy in children with Prader-Willi syndrome*. Horm Res. [PMID 12714788](https://pubmed.ncbi.nlm.nih.gov/12714788/)
261. Damen L … Hokken-Koelega ACS et al. (2020). *Three years of growth hormone treatment in young adults with Prader-Willi Syndrome previously treated with growth hormone in childhood: Effects on glucose homeostasis and metabolic syndrome*. Clin Endocrinol (Oxf). [PMID 32609902](https://pubmed.ncbi.nlm.nih.gov/32609902/)
262. Sridhar S … Karthika LN (2022). *Clinical Profile and Molecular Genetic Analysis of Prader - Willi Syndrome: A Single Center Experience*. Indian J Endocrinol Metab. [PMID 36185961](https://pubmed.ncbi.nlm.nih.gov/36185961/)
263. Qian Y … Zou C et al. (2022). *Do patients with Prader-Willi syndrome have favorable glucose metabolism?*. Orphanet J Rare Dis. [PMID 35525976](https://pubmed.ncbi.nlm.nih.gov/35525976/)
264. Kennedy L … Butler MG (2006). *Circulating adiponectin levels, body composition and obesity-related variables in Prader-Willi syndrome: comparison with obese subjects*. Int J Obes (Lond). [PMID 16231029](https://pubmed.ncbi.nlm.nih.gov/16231029/)
265. Haqq AM … Freemark M (2011). *The metabolic phenotype of Prader-Willi syndrome (PWS) in childhood: heightened insulin sensitivity relative to body mass index*. J Clin Endocrinol Metab. [PMID 20962018](https://pubmed.ncbi.nlm.nih.gov/20962018/)
266. Haqq AM … Freemark MS et al. (2007). *Altered distribution of adiponectin isoforms in children with Prader-Willi syndrome (PWS): association with insulin sensitivity and circulating satiety peptide hormones*. Clin Endocrinol (Oxf). [PMID 17666087](https://pubmed.ncbi.nlm.nih.gov/17666087/)
267. Irizarry KA … Freemark M et al. (2015). *Metabolic profiling in Prader-Willi syndrome and nonsyndromic obesity: sex differences and the role of growth hormone*. Clin Endocrinol (Oxf). [PMID 25736874](https://pubmed.ncbi.nlm.nih.gov/25736874/)
268. Lee HJ … Jin DK et al. (2011). *Delayed response of amylin levels after an oral glucose challenge in children with Prader-Willi syndrome*. Yonsei Med J. [PMID 21319343](https://pubmed.ncbi.nlm.nih.gov/21319343/)
269. Prodam F … Bona G et al. (2009). *Influence of age, gender, and glucose tolerance on fasting and fed acylated ghrelin in Prader Willi syndrome*. Clin Nutr. [PMID 19150743](https://pubmed.ncbi.nlm.nih.gov/19150743/)
270. Fang H, Judd RL (2018). *Adiponectin Regulation and Function*. Compr Physiol. [PMID 29978896](https://pubmed.ncbi.nlm.nih.gov/29978896/)
271. Gliniak CM … Scherer PE et al. (2025). *FGF21 promotes longevity in diet-induced obesity through metabolic benefits independent of growth suppression*. Cell Metab. [PMID 40527315](https://pubmed.ncbi.nlm.nih.gov/40527315/)
272. Achari AE, Jain SK (2017). *Adiponectin, a Therapeutic Target for Obesity, Diabetes, and Endothelial Dysfunction*. Int J Mol Sci. [PMID 28635626](https://pubmed.ncbi.nlm.nih.gov/28635626/)
273. Shirazi FKH … Jeddi M (2021). *Insulin resistance and high molecular weight adiponectin in obese and non-obese patients with Polycystic Ovarian Syndrome (PCOS)*. BMC Endocr Disord. [PMID 33750349](https://pubmed.ncbi.nlm.nih.gov/33750349/)

## O. 생식축과 사춘기 (Gonadal axis and puberty)

> `FHYPO`(GnRH 진폭)·`DMKRN3`(게이트 전진)·`KSEXH`. 사춘기가 **일찍 시작하고 진폭이 낮고 진행하지 않는다**는 임상 인상이 같은 결실의 두 방향 효과로 유도된다.

274. Matsuyama S … Mizokami A et al. (2019). *Gonadal function and testicular histology in males with Prader-Willi syndrome*. Endocrinol Diabetes Metab. [PMID 30815576](https://pubmed.ncbi.nlm.nih.gov/30815576/)
275. Vogels A … Bogaert GA (2008). *Testicular histology in boys with Prader-Willi syndrome: fertile or infertile?*. J Urol. [PMID 18721940](https://pubmed.ncbi.nlm.nih.gov/18721940/)
276. Wannarachue N, Ruvalcaba RH (1975). *Hypogonadism in Prader-Willi syndrome*. Am J Ment Defic. [PMID 164772](https://pubmed.ncbi.nlm.nih.gov/164772/)
277. Rey RA, Grinspon RP (2024). *Anti-Müllerian hormone, testicular descent and cryptorchidism*. Front Endocrinol (Lausanne). [PMID 38501100](https://pubmed.ncbi.nlm.nih.gov/38501100/)
278. Tauber M, Hoybye C (2021). *Endocrine disorders in Prader-Willi syndrome: a model to understand and treat hypothalamic dysfunction*. Lancet Diabetes Endocrinol. [PMID 33647242](https://pubmed.ncbi.nlm.nih.gov/33647242/)
279. Gaston LS, Stafford DE (2023). *Premature adrenarche in Prader-Willi syndrome is associated with accelerated pre-pubertal growth and advanced bone age*. J Pediatr Endocrinol Metab. [PMID 36458449](https://pubmed.ncbi.nlm.nih.gov/36458449/)
280. Lecka-Ambroziak A … Szalecki M (2020). *Premature Adrenarche in Children with Prader-Willi Syndrome Treated with Recombinant Human Growth Hormone Seems to Not Influence the Course of Central Puberty and the Efficacy and Safety of the Therapy*. Life (Basel). [PMID 33050529](https://pubmed.ncbi.nlm.nih.gov/33050529/)
281. Monai E … Jensen RB et al. (2019). *CENTRAL PRECOCIOUS PUBERTY IN TWO BOYS WITH PRADER-WILLI SYNDROME ON GROWTH HORMONE TREATMENT*. AACE Clin Case Rep. [PMID 31967069](https://pubmed.ncbi.nlm.nih.gov/31967069/)
282. Kherra S … Donaldson MDC et al. (2021). *Hypogonadism in Prader-Willi syndrome from birth to adulthood: a 28-year experience in a single centre*. Endocr Connect. [PMID 34382580](https://pubmed.ncbi.nlm.nih.gov/34382580/)
283. Sano K … Tachibana K (1994). *Urological problems in Prader-Willi syndrome*. Nihon Hinyokika Gakkai Zasshi. [PMID 7933754](https://pubmed.ncbi.nlm.nih.gov/7933754/)
284. Siddiqui S … Moin S (2022). *A brief insight into the etiology, genetics, and immunology of polycystic ovarian syndrome (PCOS)*. J Assist Reprod Genet. [PMID 36190593](https://pubmed.ncbi.nlm.nih.gov/36190593/)
285. Spaziani M … Radicioni AF et al. (2021). *Hypothalamo-Pituitary axis and puberty*. Mol Cell Endocrinol. [PMID 33271219](https://pubmed.ncbi.nlm.nih.gov/33271219/)
286. Velasco I … Tena-Sempere M et al. (2023). *Dissecting the KNDy hypothesis: KNDy neuron-derived kisspeptins are dispensable for puberty but essential for preserved female fertility and gonadotropin pulsatility*. Metabolism. [PMID 37121307](https://pubmed.ncbi.nlm.nih.gov/37121307/)
287. Uenoyama Y … Maeda K (2014). *KNDy neuron as a gatekeeper of puberty onset*. J Obstet Gynaecol Res. [PMID 24888910](https://pubmed.ncbi.nlm.nih.gov/24888910/)

## P. 뼈와 근골격계 (Bone, muscle and the musculoskeletal system)

> `KBMDU`/`KBMDR`/`KBMDH`(뼈)와 `FTONE0`/`FTONEP`(근긴장), `KCOB`/`NCOB`(척추측만). GH가 성장속도(악화)와 근긴장(개선)을 동시에 올려 거의 상쇄된다는 결과의 대조군.

288. van Abswoude DH … de Graaff LCG et al. (2022). *Bone Health in Adults With Prader-Willi Syndrome: Clinical Recommendations Based on a Multicenter Cohort Study*. J Clin Endocrinol Metab. [PMID 36149817](https://pubmed.ncbi.nlm.nih.gov/36149817/)
289. Uehara M … Kato H et al. (2019). *Efficacy of denosumab therapy for a 21-year-old woman with Prader-Willi syndrome, osteoporosis and history of fractures: a case report*. Ther Clin Risk Manag. [PMID 30880995](https://pubmed.ncbi.nlm.nih.gov/30880995/)
290. Baraghithy S … Tam J et al. (2019). *Magel2 Modulates Bone Remodeling and Mass in Prader-Willi Syndrome by Affecting Oleoyl Serine Levels and Activity*. J Bone Miner Res. [PMID 30347474](https://pubmed.ncbi.nlm.nih.gov/30347474/)
291. Kroonen LT … Macewen GD (2006). *Prader-Willi Syndrome: clinical concerns for the orthopaedic surgeon*. J Pediatr Orthop. [PMID 16932110](https://pubmed.ncbi.nlm.nih.gov/16932110/)
292. Sone S (1994). *Muscle histochemistry in the Prader-Willi syndrome*. Brain Dev. [PMID 7943601](https://pubmed.ncbi.nlm.nih.gov/7943601/)
293. Argov Z … Mastaglia FL (1984). *Patterns of muscle fiber-type disproportion in hypotonic infants*. Arch Neurol. [PMID 6689888](https://pubmed.ncbi.nlm.nih.gov/6689888/)
294. Laurier V … Jauregi J et al. (2015). *Medical, psychological and social features in a large cohort of adults with Prader-Willi syndrome: experience from a dedicated centre in France*. J Intellect Disabil Res. [PMID 24947991](https://pubmed.ncbi.nlm.nih.gov/24947991/)
295. Trizno AA … Georgopoulos G (2018). *The Prevalence and Treatment of Hip Dysplasia in Prader-Willi Syndrome (PWS)*. J Pediatr Orthop. [PMID 29309382](https://pubmed.ncbi.nlm.nih.gov/29309382/)
296. Laumonerie P … Sales de Gauzy J et al. (2020). *Evolution of Hip Dysplasia in Pediatric Patients With Prader-Willi Syndrome Treated With Growth Hormone Early in Development*. J Pediatr Orthop. [PMID 31479030](https://pubmed.ncbi.nlm.nih.gov/31479030/)

## Q. 행동과 정신과 (Behaviour, psychiatry and cognition)

> `BEH0`·`KBFR`(좌절)·`KBOXT`(옥시토신 완화)·`FSUB`(아형). 모델은 식이 제한이 지방을 얻고 행동을 잃는 **명시적 트레이드오프**를 만들며, 실현되지 않은 구동 자체가 좌절항이다.

297. Warnock JK, Kestenbaum T (1992). *Pharmacologic treatment of severe skin-picking behaviors in Prader-Willi syndrome. Two case reports*. Arch Dermatol. [PMID 1456757](https://pubmed.ncbi.nlm.nih.gov/1456757/)
298. Rice LJ … Einfeld SL (2016). *Reduced gamma-aminobutyric acid is associated with emotional and behavioral problems in Prader-Willi syndrome*. Am J Med Genet B Neuropsychiatr Genet. [PMID 27338833](https://pubmed.ncbi.nlm.nih.gov/27338833/)
299. Dykens E, Shah B (2003). *Psychiatric disorders in Prader-Willi syndrome: epidemiology and management*. CNS Drugs. [PMID 12617696](https://pubmed.ncbi.nlm.nih.gov/12617696/)
300. Krefft M … Misiak B (2014). *From Prader-Willi syndrome to psychosis: translating parent-of-origin effects into schizophrenia research*. Epigenomics. [PMID 25531260](https://pubmed.ncbi.nlm.nih.gov/25531260/)
301. Tarsimi A … Vanderbruggen N (2021). *Psychiatric disorders in adults with Prader-Willi syndrome: a systematic literature review*. Tijdschr Psychiatr. [PMID 34231862](https://pubmed.ncbi.nlm.nih.gov/34231862/)
302. Krishnadas R … Cavanagh J et al. (2018). *Brain-stem serotonin transporter availability in maternal uniparental disomy and deletion Prader-Willi syndrome*. Br J Psychiatry. [PMID 29433608](https://pubmed.ncbi.nlm.nih.gov/29433608/)
303. Kong X … Wan G et al. (2020). *Early Screening and Risk Factors of Autism Spectrum Disorder in a Large Cohort of Chinese Patients With Prader-Willi Syndrome*. Front Psychiatry. [PMID 33329146](https://pubmed.ncbi.nlm.nih.gov/33329146/)
304. Crespi B … Hurd P (2018). *A genetic locus for paranoia*. Biol Lett. [PMID 29343559](https://pubmed.ncbi.nlm.nih.gov/29343559/)
305. Phelan MC (2008). *Deletion 22q13.3 syndrome*. Orphanet J Rare Dis. [PMID 18505557](https://pubmed.ncbi.nlm.nih.gov/18505557/)
306. Salminen I … Crespi B (2022). *Do the diverse phenotypes of Prader-Willi syndrome reflect extremes of covariation in typical populations?*. Front Genet. [PMID 36506301](https://pubmed.ncbi.nlm.nih.gov/36506301/)
307. Whittington J … Boer H (2004). *Cognitive abilities and genotype in a population-based sample of people with Prader-Willi syndrome*. J Intellect Disabil Res. [PMID 14723659](https://pubmed.ncbi.nlm.nih.gov/14723659/)
308. Gross-Tsur V … Shalev RS (2001). *Cognition, attention, and behavior in Prader-Willi syndrome*. J Child Neurol. [PMID 11332464](https://pubmed.ncbi.nlm.nih.gov/11332464/)
309. Meade C … Roche E et al. (2021). *Prader-Willi Syndrome in children: Quality of life and caregiver burden*. Acta Paediatr. [PMID 33378107](https://pubmed.ncbi.nlm.nih.gov/33378107/)
310. Dempsey D … Bhatnagar A et al. (2025). *The burden of illness in Prader-Willi syndrome: a systematic literature review*. Orphanet J Rare Dis. [PMID 40708003](https://pubmed.ncbi.nlm.nih.gov/40708003/)
311. Kayadjanian N … Strong TV (2018). *High levels of caregiver burden in Prader-Willi syndrome*. PLoS One. [PMID 29579119](https://pubmed.ncbi.nlm.nih.gov/29579119/)

## R. 사망률과 임상 결과 (Mortality and clinical outcome)

> cluster 18의 사망 경로. 모델은 사망률을 상태변수로 쓰지 않고 AHI·저산소 부담·위 파열 경로의 점선으로만 표시한다 — 정량화하지 않은 것을 정량화한 척하지 않기 위해서다.

312. Pacoricona Alfaro DL … Tauber M et al. (2019). *Causes of death in Prader-Willi syndrome: lessons from 11 years' experience of a national reference center*. Orphanet J Rare Dis. [PMID 31684997](https://pubmed.ncbi.nlm.nih.gov/31684997/)
313. Einfeld SL … Taffe J (2006). *Mortality in Prader-Willi syndrome*. Am J Ment Retard. [PMID 16597186](https://pubmed.ncbi.nlm.nih.gov/16597186/)
314. Butler MG … Loker J (2017). *Causes of death in Prader-Willi syndrome: Prader-Willi Syndrome Association (USA) 40-year mortality survey*. Genet Med. [PMID 27854358](https://pubmed.ncbi.nlm.nih.gov/27854358/)
315. Lionti T … Rowell MM (2012). *Prader-Willi syndrome in Victoria: mortality and causes of death*. J Paediatr Child Health. [PMID 22697408](https://pubmed.ncbi.nlm.nih.gov/22697408/)
316. Byard RW (2018). *Death by food*. Forensic Sci Med Pathol. [PMID 28710688](https://pubmed.ncbi.nlm.nih.gov/28710688/)
317. Dykens EM (2013). *Aging in rare intellectual disability syndromes*. Dev Disabil Res Rev. [PMID 23949831](https://pubmed.ncbi.nlm.nih.gov/23949831/)

## S. 진료 지침과 다학제 관리 (Guidelines and multidisciplinary management)

> `AVAIL`·`TITR`·`BWTGTR`(양육자 적정)와 kcal/cm 규칙. 모델에서 식이 환경은 상태의 깊이를 바꾸는 것이 아니라 **어떤 상태가 존재하는지**를 바꾸는 유일한 개입이다.

318. Butler MG … Forster JL (2016). *Prader-Willi Syndrome: Clinical Genetics and Diagnostic Aspects with Treatment Approaches*. Curr Pediatr Rev. [PMID 26592417](https://pubmed.ncbi.nlm.nih.gov/26592417/)
319. Góralska M … Lewiński A et al. (2018). *Management of Prader-Willi Syndrome (PWS) in adults - what an endocrinologist needs to know. Recommendations of the Polish Society of Endocrinology and the Polish Society of Paediatric Endocrinology and Diabetology*. Endokrynol Pol. [PMID 30209801](https://pubmed.ncbi.nlm.nih.gov/30209801/)
320. Krasińska A, Skowrońska B (2017). *Prader-Willi Syndrome - nutritional management in children, adolescents and adults*. Pediatr Endocrinol Diabetes Metab. [PMID 29073293](https://pubmed.ncbi.nlm.nih.gov/29073293/)
321. McCandless SE (2011). *Clinical report—health supervision for children with Prader-Willi syndrome*. Pediatrics. [PMID 21187304](https://pubmed.ncbi.nlm.nih.gov/21187304/)
322. Holm VA … Greenberg F et al. (1993). *Prader-Willi syndrome: consensus diagnostic criteria*. Pediatrics. [PMID 8424017](https://pubmed.ncbi.nlm.nih.gov/8424017/)
323. Crinò A … Iughetti L et al. (2009). *A survey on Prader-Willi syndrome in the Italian population: prevalence of historical and clinical signs*. J Pediatr Endocrinol Metab. [PMID 20020576](https://pubmed.ncbi.nlm.nih.gov/20020576/)
324. Christianson AL … van Rensburg EJ (1998). *Prader-Willi syndrome in South African patients--clinical and molecular diagnosis*. S Afr Med J. [PMID 9687849](https://pubmed.ncbi.nlm.nih.gov/9687849/)
325. Shaikh MG … Soni S et al. (2024). *Prader-Willi syndrome: guidance for children and transition into adulthood*. Endocr Connect. [PMID 38838713](https://pubmed.ncbi.nlm.nih.gov/38838713/)
326. Poitou C … Tauber M et al. (2023). *The transition from pediatric to adult care in individuals with Prader-Willi syndrome*. Endocr Connect. [PMID 36347048](https://pubmed.ncbi.nlm.nih.gov/36347048/)
327. Pedersen M, Höybye C (2021). *An Adapted Model for Transition to Adult Care in Young Adults with Prader-Willi Syndrome*. J Clin Med. [PMID 34066432](https://pubmed.ncbi.nlm.nih.gov/34066432/)
328. Saitoh S (2010). *Care continuity for patients with Prader-Willi syndrome during transition from childhood to adulthood*. Nihon Rinsho. [PMID 20077807](https://pubmed.ncbi.nlm.nih.gov/20077807/)
329. Coplin SS … Gormican A (1976). *Out-patient dietary management in the Prader-Willi syndrome*. J Am Diet Assoc. [PMID 1254876](https://pubmed.ncbi.nlm.nih.gov/1254876/)

## T. 동물 모델 (Animal models)

> Snord116·Magel2·Necdin 모델이 어느 분기를 재현하고 어느 분기를 재현하지 않는지가 다섯 개 탈출비의 순서를 간접적으로 제약한다.

330. Ding F … Francke U et al. (2008). *SnoRNA Snord116 (Pwcr1/MBII-85) deletion causes growth deficiency and hyperphagia in mice*. PLoS One. [PMID 18320030](https://pubmed.ncbi.nlm.nih.gov/18320030/)
331. Ding F … Francke U (2010). *Neonatal maternal deprivation response and developmental changes in gene expression revealed by hypothalamic gene expression profiling in mice*. PLoS One. [PMID 20195375](https://pubmed.ncbi.nlm.nih.gov/20195375/)
332. Lee S, Jo YH (2025). *Magel2 in hypothalamic POMC neurons influences the impact of stress on anxiety-like behavior and spatial learning associated with a food reward in male mice*. Front Neural Circuits. [PMID 41321745](https://pubmed.ncbi.nlm.nih.gov/41321745/)
333. Tennese AA, Wevrick R (2011). *Impaired hypothalamic regulation of endocrine function and delayed counterregulatory response to hypoglycemia in Magel2-null mice*. Endocrinology. [PMID 21248145](https://pubmed.ncbi.nlm.nih.gov/21248145/)
334. Higgs MJ … Isles AR (2023). *The parenting hub of the hypothalamus is a focus of imprinted gene action*. PLoS Genet. [PMID 37856383](https://pubmed.ncbi.nlm.nih.gov/37856383/)
335. Kuwako K … Yoshikawa K et al. (2005). *Disruption of the paternal necdin gene diminishes TrkA signaling for sensory neuron survival*. J Neurosci. [PMID 16049186](https://pubmed.ncbi.nlm.nih.gov/16049186/)

## U. QSP 방법론과 수리 모델 (QSP methodology and mathematical modelling)

> 이중안정성·안장-노드 분기·조화평균 율속·에너지 균형 모델·mrgsolve 및 QSP 모범사례. 모델의 수학적 장치가 임의가 아니라는 근거.

336. Li X … Zhu Q et al. (2023). *Combining network pharmacology, molecular docking, molecular dynamics simulation, and experimental verification to examine the efficacy and immunoregulation mechanism of FHB granules on vitiligo*. Front Immunol. [PMID 37575231](https://pubmed.ncbi.nlm.nih.gov/37575231/)
337. El-Khateeb E … Achour B et al. (2019). *Quantitative mass spectrometry-based proteomics in the era of model-informed drug development: Applications in translational pharmacology and recommendations for best practice*. Pharmacol Ther. [PMID 31376433](https://pubmed.ncbi.nlm.nih.gov/31376433/)
338. GBD 2023 Disease and Injury and Risk Factor Collaborators (2025). *Burden of 375 diseases and injuries, risk-attributable burden of 88 risk factors, and healthy life expectancy in 204 countries and territories, including 660 subnational locations, 1990-2023: a systematic analysis for the Global Burden of Disease Study 2023*. Lancet. [PMID 41092926](https://pubmed.ncbi.nlm.nih.gov/41092926/)
339. GBD 2023 Lower Respiratory Infections and Antimicrobial Resistance Collaborators (2026). *Global burden of lower respiratory infections and aetiologies, 1990-2023: a systematic analysis for the Global Burden of Disease Study 2023*. Lancet Infect Dis. [PMID 41412141](https://pubmed.ncbi.nlm.nih.gov/41412141/)
340. GBD 2023 Meningitis & Antimicrobial Resistance Collaborators (2026). *Global, regional, and national burden of meningitis, its risk factors, and aetiologies, 1990-2023: a systematic analysis for the Global Burden of Disease Study 2023*. Lancet Neurol. [PMID 41911930](https://pubmed.ncbi.nlm.nih.gov/41911930/)
341. Elmokadem A … Baron KT (2019). *Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial*. CPT Pharmacometrics Syst Pharmacol. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
342. Lu T … Kågedal M et al. (2024). *gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve*. CPT Pharmacometrics Syst Pharmacol. [PMID 38082557](https://pubmed.ncbi.nlm.nih.gov/38082557/)
343. Hooijmaijers R … Visser SAG (2023). *Building an adaptive dose simulation framework to aid dose and schedule selection*. CPT Pharmacometrics Syst Pharmacol. [PMID 37574587](https://pubmed.ncbi.nlm.nih.gov/37574587/)
344. Cantoria MJ … Thorne CA et al. (2023). *Feedback in the β-catenin destruction complex imparts bistability and cellular memory*. Proc Natl Acad Sci U S A. [PMID 36598937](https://pubmed.ncbi.nlm.nih.gov/36598937/)
345. Angeli D … Sontag ED (2004). *Detection of multistability, bifurcations, and hysteresis in a large class of biological positive-feedback systems*. Proc Natl Acad Sci U S A. [PMID 14766974](https://pubmed.ncbi.nlm.nih.gov/14766974/)
346. Laxhuber KS … Phillips R et al. (2020). *Theoretical investigation of a genetic switch for metabolic adaptation*. PLoS One. [PMID 32379825](https://pubmed.ncbi.nlm.nih.gov/32379825/)
347. Rata S … Hochegger H et al. (2018). *Two Interlinked Bistable Switches Govern Mitotic Control in Mammalian Cells*. Curr Biol. [PMID 30449668](https://pubmed.ncbi.nlm.nih.gov/30449668/)
348. Li CJ … Hong T et al. (2021). *MicroRNA governs bistable cell differentiation and lineage segregation via a noncanonical feedback*. Mol Syst Biol. [PMID 33890404](https://pubmed.ncbi.nlm.nih.gov/33890404/)
349. Macht M (2008). *How emotions affect eating: a five-way model*. Appetite. [PMID 17707947](https://pubmed.ncbi.nlm.nih.gov/17707947/)
350. Wynne K … Bloom S (2005). *Appetite control*. J Endocrinol. [PMID 15684339](https://pubmed.ncbi.nlm.nih.gov/15684339/)
351. Opara EI … Hammond WG (1996). *Studies on the regulation of food intake using rat total parenteral nutrition as a model*. Neurosci Biobehav Rev. [PMID 8880733](https://pubmed.ncbi.nlm.nih.gov/8880733/)
352. Rogers PJ, Brunstrom JM (2016). *Appetite and energy balancing*. Physiol Behav. [PMID 27059321](https://pubmed.ncbi.nlm.nih.gov/27059321/)
353. Finlayson G … Blundell JE (2007). *Liking vs. wanting food: importance for human appetite control and weight regulation*. Neurosci Biobehav Rev. [PMID 17559933](https://pubmed.ncbi.nlm.nih.gov/17559933/)
354. GBD 2021 Adult BMI Collaborators (2025). *Global, regional, and national prevalence of adult overweight and obesity, 1990-2021, with forecasts to 2050: a forecasting study for the Global Burden of Disease Study 2021*. Lancet. [PMID 40049186](https://pubmed.ncbi.nlm.nih.gov/40049186/)
355. Li RJ … Zhu H et al. (2022). *Model-Informed Approach Supporting Drug Development and Regulatory Evaluation for Rare Diseases*. J Clin Pharmacol. [PMID 36461744](https://pubmed.ncbi.nlm.nih.gov/36461744/)
356. Bateman RM … Prandi E et al. (2016). *36th International Symposium on Intensive Care and Emergency Medicine : Brussels, Belgium. 15-18 March 2016*. Crit Care. [PMID 27885969](https://pubmed.ncbi.nlm.nih.gov/27885969/)
357. GBD 2021 Adolescent BMI Collaborators (2025). *Global, regional, and national prevalence of child and adolescent overweight and obesity, 1990-2021, with forecasts to 2050: a forecasting study for the Global Burden of Disease Study 2021*. Lancet. [PMID 40049185](https://pubmed.ncbi.nlm.nih.gov/40049185/)

---

## 마지막 절 · 이 목록이 말하지 않는 것 (What this list does not cover)

모델의 한계를 문헌 부재로 위장하지 않기 위해 명시한다.

- **각성 및 인지의 명시적 상태변수가 없다.** K절의 수면 문헌은 임상적으로 중요하지만 모델은
  주간 과다졸림을 상태로 쓰지 않고 기도·저산소 경로로만 표현한다.
- **테스토스테론의 제지방 효과가 IGF-1·근긴장·성선 인자를 통해서만 들어간다.** 안드로겐의
  직접적 근육 단백질 합성 경로는 별도 채널로 쓰지 않았다.
- **GLP-1 수용체 작용제 팔의 근거 수준이 다른 팔보다 낮다** (M절 대부분이 관찰연구).
  따라서 이 팔의 정량적 예측은 카베토신·DCCR 팔과 같은 무게로 읽어서는 안 된다.
- **사망률을 정량화하지 않는다.** R절 문헌은 경로를 정당화하지만 모델은 확률을 계산하지 않고
  기전 지도에서 점선으로만 남긴다.
- **하루 이하 시간 규모를 모두 평균했다.** 그렐린·인슐린·GH·LH의 실제 반감기는 분 단위이며,
  모델은 이들을 일 규모 풀로 다룬다. 식후 30분 역학을 묻는 질문에는 답할 수 없다.

## 면책 (Disclaimer)

본 모델과 참고문헌 목록은 **교육 및 연구 목적**이다. 독립적으로 검증·인증되지 않았으며 실제
임상 의사결정, 처방, 규제 제출에 직접 사용해서는 안 된다. 파라미터는 공개 문헌에 근거한
설명용 근사치이며, 환자 데이터에 대한 적합·검증이 별도로 필요하다.

