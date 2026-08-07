# Neonatal opioid withdrawal syndrome (NOWS / NAS) — reference list

Every PMID below was resolved against the PubMed E-utilities API while this
model was being written, and the author / year / journal / title shown here is
the record PubMed returned — not a recollection. Where a paper is used to fix a
specific number in `nows_mrgsolve_model.R`, the parameter it constrains is
named in the annotation.

`https://pubmed.ncbi.nlm.nih.gov/<PMID>/`

---

## 1. Epidemiology, definitions and clinical overview

| # | Reference | PMID |
|---|---|---|
| 1 | Patrick SW et al. Neonatal abstinence syndrome and associated health care expenditures: United States, 2000–2009. *JAMA* 2012 | [22546608](https://pubmed.ncbi.nlm.nih.gov/22546608/) |
| 2 | Patrick SW et al. Increasing incidence and geographic distribution of neonatal abstinence syndrome: United States 2009 to 2012. *J Perinatol* 2015 | [26219703](https://pubmed.ncbi.nlm.nih.gov/26219703/) |
| 3 | Hudak ML, Tan RC. Neonatal drug withdrawal (AAP clinical report). *Pediatrics* 2012 | [22291123](https://pubmed.ncbi.nlm.nih.gov/22291123/) |
| 4 | Wachman EM, Schiff DM, Silverstein M. Neonatal abstinence syndrome: advances in diagnosis and treatment. *JAMA* 2018 | [29614184](https://pubmed.ncbi.nlm.nih.gov/29614184/) |
| 5 | Devlin LA et al. Neonatal opioid withdrawal syndrome: a review of the science and a look toward the future. 2021 | [34556799](https://pubmed.ncbi.nlm.nih.gov/34556799/) |
| 6 | Devlin LA, Davis JM. A practical approach to neonatal opiate withdrawal syndrome. 2018 | [29100261](https://pubmed.ncbi.nlm.nih.gov/29100261/) |
| 7 | Cheng FY et al. Neonatal opioid withdrawal syndrome. *Pediatr Clin North Am* 2025 | [40619192](https://pubmed.ncbi.nlm.nih.gov/40619192/) |
| 8 | Spence K et al. Non-pharmacologic and pharmacologic care of the neonate with opioid withdrawal syndrome. *Semin Perinatol* 2025 | [39706694](https://pubmed.ncbi.nlm.nih.gov/39706694/) |
| 9 | Plouffe R et al. Neonatal abstinence syndrome hospitalizations in Canada: a descriptive study. *Can J Public Health* 2023 | [36482143](https://pubmed.ncbi.nlm.nih.gov/36482143/) |
| 10 | Weikel BW et al. State initiatives to improve care for infants with prenatal substance exposure. *Neoreviews* 2025 | [40164206](https://pubmed.ncbi.nlm.nih.gov/40164206/) |
| 11 | Sarkar S, Donn SM. Management of neonatal abstinence syndrome in neonatal intensive care units: a national survey. *J Perinatol* 2006 | [16355103](https://pubmed.ncbi.nlm.nih.gov/16355103/) |
| 12 | Dickes L et al. Potential for Medicaid savings: a state and national comparison of an innovative neonatal abstinence syndrome intervention. *Popul Health Manag* 2017 | [28409699](https://pubmed.ncbi.nlm.nih.gov/28409699/) |

---

## 2. The measurement instruments — Finnegan, its successors, and Eat-Sleep-Console

The model treats scoring as a separate organ (`SCNS`, `SANS`, `SGI`, `KFENV`,
`KFS`, `THRSTART`, `ESCFAIL`). These papers are why.

| # | Reference | PMID |
|---|---|---|
| 13 | Finnegan LP et al. Neonatal abstinence syndrome: assessment and management. *Addict Dis* 1975 — the 21-item instrument itself | [1163358](https://pubmed.ncbi.nlm.nih.gov/1163358/) |
| 14 | Timpson W et al. A quality improvement initiative to increase scoring consistency and accuracy of the Finnegan tool. *Adv Neonatal Care* 2018 — inter-rater variance, constrains the threshold-as-random-variable assumption | [29045256](https://pubmed.ncbi.nlm.nih.gov/29045256/) |
| 15 | Kocherlakota P et al. A new scoring system for the assessment of neonatal abstinence syndrome. *Am J Perinatol* 2020 | [31777045](https://pubmed.ncbi.nlm.nih.gov/31777045/) |
| 16 | Casavant SG et al. Integrated review of the assessment of newborns with neonatal abstinence syndrome. *J Obstet Gynecol Neonatal Nurs* 2021 | [34116058](https://pubmed.ncbi.nlm.nih.gov/34116058/) |
| 17 | Grossman MR et al. An initiative to improve the quality of care of infants with neonatal abstinence syndrome. *Pediatrics* 2017 — the origin of the function-based approach | [28562267](https://pubmed.ncbi.nlm.nih.gov/28562267/) |
| 18 | Young LW et al. Eat, Sleep, Console approach or usual care for neonatal opioid withdrawal. *N Engl J Med* 2023 — **ESC-NOW**; 52.0% → 19.5% pharmacotherapy, 14.9 → 8.2 days to medically ready | [37125831](https://pubmed.ncbi.nlm.nih.gov/37125831/) |
| 19 | Hein S et al. Eat, Sleep, Console and adjunctive buprenorphine improved outcomes in neonatal opioid withdrawal syndrome. *Adv Neonatal Care* 2021 | [33278102](https://pubmed.ncbi.nlm.nih.gov/33278102/) |
| 20 | Ouyang L et al. Population-level changes in infant outcomes associated with a state Eat Sleep Console initiative. *Hosp Pediatr* 2025 | [40721219](https://pubmed.ncbi.nlm.nih.gov/40721219/) |
| 21 | Ghadiali M et al. Implementing the Eat, Sleep, Console model in a small safety-net hospital. *Cureus* 2026 | [42434680](https://pubmed.ncbi.nlm.nih.gov/42434680/) |
| 22 | Devlin LA et al. Temporal trends in outcomes for infants with neonatal opioid withdrawal managed with the Finnegan scoring system. *J Perinatol* 2026 | [41991673](https://pubmed.ncbi.nlm.nih.gov/41991673/) |

---

## 3. Maternal opioid use disorder treatment and the MOTHER comparison

| # | Reference | PMID |
|---|---|---|
| 23 | Jones HE et al. Neonatal abstinence syndrome after methadone or buprenorphine exposure. *N Engl J Med* 2010 — **MOTHER**; buprenorphine-exposed neonates needed less morphine and a shorter course at a similar treated fraction | [21142534](https://pubmed.ncbi.nlm.nih.gov/21142534/) |
| 24 | Jones HE et al. The relationship between maternal methadone dose at delivery and neonatal outcome: methodological and design considerations. *Neurotoxicol Teratol* 2013 | [24099621](https://pubmed.ncbi.nlm.nih.gov/24099621/) |
| 25 | Cleary BJ et al. Methadone dose and neonatal abstinence syndrome — systematic review and meta-analysis. *Addiction* 2010 — the null that the saturating adaptation map (`AGAIN`, `EA50A`) is built to explain | [20840198](https://pubmed.ncbi.nlm.nih.gov/20840198/) |
| 26 | Ordean A, Tubman-Broeren M. Safety and efficacy of buprenorphine–naloxone in pregnancy: a systematic review. *Pathophysiology* 2023 | [36810423](https://pubmed.ncbi.nlm.nih.gov/36810423/) |
| 27 | Akbarzadeh F et al. Buprenorphine-based treatment outcomes in pregnant opioid-dependent women: systematic review and meta-analysis. *Addiction* 2026 | [41345707](https://pubmed.ncbi.nlm.nih.gov/41345707/) |
| 28 | Wilson LA et al. Comparative effects of prenatal exposure to methadone versus buprenorphine. *BMC Pregnancy Childbirth* 2026 | [42204512](https://pubmed.ncbi.nlm.nih.gov/42204512/) |
| 29 | Anderer S. Extended-release buprenorphine found safe and effective during pregnancy. *JAMA* 2026 | [41931283](https://pubmed.ncbi.nlm.nih.gov/41931283/) |
| 30 | Extended-release vs sublingual buprenorphine in pregnancy through 12 months. 2026 | [41837971](https://pubmed.ncbi.nlm.nih.gov/41837971/) |
| 31 | Naltrexone compared with buprenorphine or methadone in pregnancy: a systematic review. 2024 | [38227945](https://pubmed.ncbi.nlm.nih.gov/38227945/) |
| 32 | Atoui Z et al. Evidence-informed approaches to medication management for opioid use disorder in special populations. *J Addict Dis* 2026 | [42411921](https://pubmed.ncbi.nlm.nih.gov/42411921/) |

---

## 4. Placental transfer, cord concentrations and pregnancy pharmacokinetics

Constrains `CORD_D`, `CORD_B`, `KMPD`, `KMPB`, `CNB0`.

| # | Reference | PMID |
|---|---|---|
| 33 | Coles LD et al. Distribution of saquinavir, methadone and buprenorphine in maternal brain, placenta and fetus. *J Pharm Sci* 2009 | [19116954](https://pubmed.ncbi.nlm.nih.gov/19116954/) |
| 34 | Samiee-Zafarghandy S et al. Pharmacometric evaluation of umbilical cord blood concentration-based early initiation of treatment in methadone-exposed neonates. *Children* 2021 | [33668712](https://pubmed.ncbi.nlm.nih.gov/33668712/) |
| 35 | Zaidi SS et al. Pharmacometric modeling of opioid disposition in pregnancy: a systematic review of PopPK and PBPK approaches. *Front Pharmacol* 2026 | [42088578](https://pubmed.ncbi.nlm.nih.gov/42088578/) |
| 36 | Dickmann LJ, Isoherranen N. Quantitative prediction of CYP2B6 induction by estradiol during pregnancy. *Drug Metab Dispos* 2013 — why maternal trough falls in the third trimester | [22815312](https://pubmed.ncbi.nlm.nih.gov/22815312/) |
| 37 | Esposito DB et al. Characteristics of prescription opioid analgesics in pregnancy and risk of neonatal opioid withdrawal syndrome. *JAMA Netw Open* 2022 | [36001312](https://pubmed.ncbi.nlm.nih.gov/36001312/) |

---

## 5. Neonatal pharmacokinetics — morphine, M6G, methadone, buprenorphine, clonidine, phenobarbital

Constrains `CLMREF`, `V1MREF`, `V2MREF`, `FORAL`, `FM6G`, `CLGREF`, `CLDREF`,
`VDREF`, `CLBREF`, `VBREF`, `CLNREF`, `CLCREF`, `CLPREF`.

| # | Reference | PMID |
|---|---|---|
| 38 | Bouwmeester NJ et al. Developmental pharmacokinetics of morphine and its metabolites in neonates, infants and young children. *Br J Anaesth* 2004 | [14722170](https://pubmed.ncbi.nlm.nih.gov/14722170/) |
| 39 | Knibbe CAJ et al. Morphine glucuronidation in preterm neonates, infants and children younger than 3 years. *Clin Pharmacokinet* 2009 — the PMA-driven glucuronidation model behind `PMA50U`/`HILLU` | [19650676](https://pubmed.ncbi.nlm.nih.gov/19650676/) |
| 40 | Anand KJS et al. Morphine pharmacokinetics and pharmacodynamics in preterm and term neonates (NEOPAIN secondary results). *Br J Anaesth* 2008 | [18723857](https://pubmed.ncbi.nlm.nih.gov/18723857/) |
| 41 | Krekels EHJ et al. Predictive performance of a recently developed population pharmacokinetic model for morphine and its metabolites. *Clin Pharmacokinet* 2011 | [21142267](https://pubmed.ncbi.nlm.nih.gov/21142267/) |
| 42 | Verscheijden LFM et al. PBPK/PD model for the prediction of morphine brain disposition and analgesia in adults and children. *PLoS Comput Biol* 2021 — brain penetration of morphine vs M6G, behind `RPG` | [33661919](https://pubmed.ncbi.nlm.nih.gov/33661919/) |
| 43 | van Donge T et al. Methadone dosing strategies in preterm neonates can be simplified. *Br J Clin Pharmacol* 2019 | [30805946](https://pubmed.ncbi.nlm.nih.gov/30805946/) |
| 44 | Ng CM et al. Population pharmacokinetic model of sublingual buprenorphine in neonatal abstinence syndrome. *Pharmacotherapy* 2015 | [26172282](https://pubmed.ncbi.nlm.nih.gov/26172282/) |
| 45 | Kraft WK et al. Revised dose schema of sublingual buprenorphine in the treatment of the neonatal opioid abstinence syndrome. *Addiction* 2011 | [20925688](https://pubmed.ncbi.nlm.nih.gov/20925688/) |
| 46 | Alsmadi MM et al. Salivary therapeutic monitoring of buprenorphine in neonates after maternal sublingual dosing guided by PBPK. *Ther Drug Monit* 2024 | [38366333](https://pubmed.ncbi.nlm.nih.gov/38366333/) |
| 47 | Xie HG et al. Clonidine clearance matures rapidly during the early postnatal period: a population pharmacokinetic analysis. *J Clin Pharmacol* 2011 — behind `CLCREF` and its GFR coupling | [20484620](https://pubmed.ncbi.nlm.nih.gov/20484620/) |
| 48 | van Hoogdalem MW et al. Pharmacotherapy of neonatal opioid withdrawal syndrome: a review of pharmacokinetics and pharmacodynamics. *Expert Opin Drug Metab Toxicol* 2021 | [33049155](https://pubmed.ncbi.nlm.nih.gov/33049155/) |
| 49 | McPhail BT et al. Opioid treatment for neonatal opioid withdrawal syndrome: current challenges and future approaches. *J Clin Pharmacol* 2021 | [33382111](https://pubmed.ncbi.nlm.nih.gov/33382111/) |
| 50 | Mahdy WYB et al. Evaluation of fentanyl-emerged adverse events and pharmacokinetics in neonates: a PBPK approach. *Clin Pharmacokinet* 2025 | [40999188](https://pubmed.ncbi.nlm.nih.gov/40999188/) |
| 51 | Olofsen E et al. Morphine and hydromorphone pharmacodynamics in human volunteers. *Br J Anaesth* 2026 — effect-site equilibration, behind `KE0M` | [41656122](https://pubmed.ncbi.nlm.nih.gov/41656122/) |
| 52 | Morse JD et al. Pharmacokinetic modeling and simulation to understand diamorphine dose–response in neonates and children. *Paediatr Anaesth* 2022 | [35212432](https://pubmed.ncbi.nlm.nih.gov/35212432/) |
| 53 | Ing Lorenzini K et al. Pharmacokinetic–pharmacodynamic modelling of opioids in healthy human volunteers: a minireview. *Basic Clin Pharmacol Toxicol* 2012 | [21995512](https://pubmed.ncbi.nlm.nih.gov/21995512/) |
| 54 | Favié LMA et al. Prediction of drug exposure in critically ill encephalopathic neonates. *Clin Pharmacol Ther* 2020 | [32463940](https://pubmed.ncbi.nlm.nih.gov/32463940/) |
| 55 | Yalcin N et al. Population pharmacokinetics in critically ill neonates and infants. *BMJ Paediatr Open* 2022 | [36437518](https://pubmed.ncbi.nlm.nih.gov/36437518/) |

---

## 6. Maturation and allometry — the arithmetic behind the hidden taper

| # | Reference | PMID |
|---|---|---|
| 56 | Anderson BJ, Holford NHG. Mechanism-based concepts of size and maturity in pharmacokinetics. *Annu Rev Pharmacol Toxicol* 2008 — the (WT/70)^0.75 × sigmoid-PMA structure used throughout | [17914927](https://pubmed.ncbi.nlm.nih.gov/17914927/) |
| 57 | Holford N et al. A pharmacokinetic standard for babies and adults. *J Pharm Sci* 2013 | [23650116](https://pubmed.ncbi.nlm.nih.gov/23650116/) |
| 58 | Germovsek E et al. Scaling clearance in paediatric pharmacokinetics: all models are wrong, which are useful? *Br J Clin Pharmacol* 2017 | [27767204](https://pubmed.ncbi.nlm.nih.gov/27767204/) |
| 59 | Rhodin MM et al. Human renal function maturation: a quantitative description using weight and postmenstrual age. *Pediatr Nephrol* 2009 — behind `PMA50R`/`HILLR` | [18846389](https://pubmed.ncbi.nlm.nih.gov/18846389/) |
| 60 | Wu Y et al. Prediction of glomerular filtration rate maturation across preterm and term neonates and young infants. *AAPS J* 2022 | [35212832](https://pubmed.ncbi.nlm.nih.gov/35212832/) |
| 61 | Badée J et al. Characterization of the ontogeny of hepatic UDP-glucuronosyltransferase enzymes. *J Clin Pharmacol* 2019 | [31502688](https://pubmed.ncbi.nlm.nih.gov/31502688/) |

---

## 7. Neurobiology — adenylyl cyclase superactivation, the locus coeruleus, GIRK, and receptor regulation

This is the cluster that justifies writing withdrawal as `GAP = A − ITONE`
rather than as a function of concentration.

| # | Reference | PMID |
|---|---|---|
| 62 | Nestler EJ, Aghajanian GK. Molecular and cellular basis of addiction. *Science* 1997 | [9311927](https://pubmed.ncbi.nlm.nih.gov/9311927/) |
| 63 | Nestler EJ. Reflections on: "A general role for adaptations in G-proteins and the cyclic AMP system in mediating the chronic actions of morphine and cocaine". *Brain Res* 2016 | [26740398](https://pubmed.ncbi.nlm.nih.gov/26740398/) |
| 64 | Kogan JH, Nestler EJ, Aghajanian GK. Elevated basal firing rates and enhanced responses to 8-Br-cAMP in locus coeruleus neurons in brain slices from opiate-dependent rats. *Eur J Pharmacol* 1992 — the direct measurement of the set-point A | [1618268](https://pubmed.ncbi.nlm.nih.gov/1618268/) |
| 65 | Kogan JH, Aghajanian GK. Long-term glutamate desensitization in locus coeruleus neurons and its role in opiate withdrawal. *Brain Res* 1995 | [8528694](https://pubmed.ncbi.nlm.nih.gov/8528694/) |
| 66 | Christie MJ. Cellular neuroadaptations to chronic opioids: tolerance, withdrawal and addiction. *Br J Pharmacol* 2008 | [18414400](https://pubmed.ncbi.nlm.nih.gov/18414400/) |
| 67 | Williams JT et al. Regulation of µ-opioid receptors: desensitization, phosphorylation, internalization and tolerance. *Pharmacol Rev* 2013 — behind `KRDOWN`/`KRREC` and `KSHIFT` | [23321159](https://pubmed.ncbi.nlm.nih.gov/23321159/) |
| 68 | Underwood O et al. Key phosphorylation sites for robust β-arrestin2 binding at the MOR revisited. *Commun Biol* 2024 | [39095612](https://pubmed.ncbi.nlm.nih.gov/39095612/) |
| 69 | Groom S et al. A novel G protein-biased agonist at the µ opioid receptor induces substantial receptor desensitisation. *Br J Pharmacol* 2023 | [33245558](https://pubmed.ncbi.nlm.nih.gov/33245558/) |
| 70 | Zhao H et al. EGFR-dependent subcellular communication was responsible for morphine-mediated AC superactivation. *Cell Signal* 2013 | [23142605](https://pubmed.ncbi.nlm.nih.gov/23142605/) |
| 71 | Jalali Mashayekhi F et al. Expression levels of the tyrosine hydroxylase gene and histone modifications around its promoter in the locus coeruleus. *Eur Addict Res* 2018 | [30517913](https://pubmed.ncbi.nlm.nih.gov/30517913/) |
| 72 | Llorca-Torralba M et al. Opioid activity in the locus coeruleus is modulated by chronic pain. *Mol Neurobiol* 2019 | [30284123](https://pubmed.ncbi.nlm.nih.gov/30284123/) |
| 73 | Kwok CHT et al. Pannexin-1 channel inhibition alleviates opioid withdrawal in rodents by modulating locus coeruleus to spinal cord circuitry. *Nat Commun* 2024 | [39048565](https://pubmed.ncbi.nlm.nih.gov/39048565/) |
| 74 | Foster SL et al. Cell-type specific expression and behavioral impact of galanin and GalR1 in the locus coeruleus during opioid withdrawal. *Addict Biol* 2021 | [33768673](https://pubmed.ncbi.nlm.nih.gov/33768673/) |
| 75 | Akbarian S et al. BDNF is essential for opiate-induced plasticity of noradrenergic neurons. *J Neurosci* 2002 | [12019333](https://pubmed.ncbi.nlm.nih.gov/12019333/) |
| 76 | The buprenorphine paradox: how buprenorphine triggers and resolves opioid withdrawal. 2026 | [41802339](https://pubmed.ncbi.nlm.nih.gov/41802339/) |
| 77 | Kozińska RB et al. Can buprenorphine be overdosed? The ceiling effect and its clinical implications. *Pharmaceuticals* 2026 — behind `EMAXB` | [42356521](https://pubmed.ncbi.nlm.nih.gov/42356521/) |
| 78 | Bahji A et al. Reframing buprenorphine as a pharmacologic modifier of opioid-induced respiratory depression. *Pharmaceuticals* 2026 | [42198473](https://pubmed.ncbi.nlm.nih.gov/42198473/) |

---

## 8. Pharmacologic treatment trials and comparative effectiveness

| # | Reference | PMID |
|---|---|---|
| 79 | Kraft WK et al. Buprenorphine for the treatment of the neonatal abstinence syndrome. *N Engl J Med* 2017 — median treatment 15 vs 28 days, length of stay 21 vs 33 days | [28468518](https://pubmed.ncbi.nlm.nih.gov/28468518/) |
| 80 | Kraft WK et al. Buprenorphine for the neonatal abstinence syndrome (correspondence). *N Engl J Med* 2017 | [28877016](https://pubmed.ncbi.nlm.nih.gov/28877016/) |
| 81 | Hall ES et al. Comparison of neonatal abstinence syndrome treatment with sublingual buprenorphine versus conventional opioids. *Am J Perinatol* 2018 | [29112997](https://pubmed.ncbi.nlm.nih.gov/29112997/) |
| 82 | Agthe AG et al. Clonidine as an adjunct therapy to opioids for neonatal abstinence syndrome: a randomized, controlled trial. *Pediatrics* 2009 — 11 vs 15 days, behind `EMAXA2`/`EC50A2` | [19398463](https://pubmed.ncbi.nlm.nih.gov/19398463/) |
| 83 | D'Abaco E et al. Does the addition of clonidine to opioid therapy improve outcomes in infants with neonatal abstinence syndrome? *J Paediatr Child Health* 2021 | [33493373](https://pubmed.ncbi.nlm.nih.gov/33493373/) |
| 84 | Brusseau C et al. Clonidine versus phenobarbital as adjunctive therapy for neonatal abstinence syndrome. *J Perinatol* 2020 | [32424335](https://pubmed.ncbi.nlm.nih.gov/32424335/) |
| 85 | Coyle MG et al. Diluted tincture of opium and phenobarbital versus DTO alone for neonatal opiate withdrawal. *J Pediatr* 2002 | [12032522](https://pubmed.ncbi.nlm.nih.gov/12032522/) |
| 86 | Sutter MB et al. Morphine versus methadone for neonatal opioid withdrawal syndrome: a randomized controlled pilot study. *BMC Pediatr* 2022 | [35705944](https://pubmed.ncbi.nlm.nih.gov/35705944/) |
| 87 | Zankl A et al. Opioid treatment for opioid withdrawal in newborn infants. *Cochrane Database Syst Rev* 2021 | [34231914](https://pubmed.ncbi.nlm.nih.gov/34231914/) |
| 88 | Devlin LA et al. Symptom-based dosing for neonatal opioid withdrawal: the OPTimize NOW randomized clinical trial. *JAMA* 2026 | [42033722](https://pubmed.ncbi.nlm.nih.gov/42033722/) |
| 89 | Optimizing pharmacologic treatment for neonatal opioid withdrawal syndrome (OPTimize NOW) — protocol. 2025 | [40866977](https://pubmed.ncbi.nlm.nih.gov/40866977/) |
| 90 | Pérez-Jiménez JM et al. Therapeutic update in neonatal opioid withdrawal syndrome: comparative effectiveness of pharmacological treatments. *Health Sci Rep* 2026 | [42005669](https://pubmed.ncbi.nlm.nih.gov/42005669/) |
| 91 | Ahmad M et al. Interventions to reduce pharmacologic opioid exposure in neonatal opioid withdrawal syndrome: a systematic review. *Front Pain Res* 2026 | [42137907](https://pubmed.ncbi.nlm.nih.gov/42137907/) |
| 92 | Urlesberger B et al. Acupuncture for neonatal abstinence syndrome in newborn infants. *Cochrane Database Syst Rev* 2025 | [39981752](https://pubmed.ncbi.nlm.nih.gov/39981752/) |

---

## 9. Weaning protocols and practice variation — where the 3-fold duration spread comes from

| # | Reference | PMID |
|---|---|---|
| 93 | Hall ES et al. Implementation of a neonatal abstinence syndrome weaning protocol: a multicenter cohort study. *Pediatrics* 2015 | [26371196](https://pubmed.ncbi.nlm.nih.gov/26371196/) |
| 94 | Wachman EM et al. Standard fixed-schedule methadone taper versus symptom-triggered methadone approach. *Hosp Pediatr* 2019 | [31270130](https://pubmed.ncbi.nlm.nih.gov/31270130/) |
| 95 | Isemann B et al. Maternal and neonatal factors impacting response to methadone therapy in infants treated for neonatal abstinence syndrome. *J Perinatol* 2011 | [20508596](https://pubmed.ncbi.nlm.nih.gov/20508596/) |
| 96 | Improving care for neonatal abstinence syndrome. 2016 | [27244809](https://pubmed.ncbi.nlm.nih.gov/27244809/) |
| 97 | Parlaman J et al. Improving care for infants with neonatal abstinence syndrome: a multicenter, community hospital-based study. *Hosp Pediatr* 2019 | [31308049](https://pubmed.ncbi.nlm.nih.gov/31308049/) |
| 98 | Overview of perinatal quality collaboratives and their activities. 2026 | [41932817](https://pubmed.ncbi.nlm.nih.gov/41932817/) |
| 99 | Hughes MV et al. Team-led empowerment: a toolkit enhancing NOWS care. *J Perinat Neonatal Nurs* 2025 | [39325948](https://pubmed.ncbi.nlm.nih.gov/39325948/) |

---

## 10. Non-pharmacologic care, rooming-in and breastfeeding

Constrains `CARE`, `KCG`, `KCARE`, `RIDD`, `RIDB`.

| # | Reference | PMID |
|---|---|---|
| 100 | Pahl A et al. Non-pharmacological care for opioid withdrawal in newborns. *Cochrane Database Syst Rev* 2020 | [33348423](https://pubmed.ncbi.nlm.nih.gov/33348423/) |
| 101 | Holmes AV et al. Rooming-in to treat neonatal abstinence syndrome: improved family-centered care at lower cost. *Pediatrics* 2016 | [27194629](https://pubmed.ncbi.nlm.nih.gov/27194629/) |
| 102 | McKnight S et al. Rooming-in for infants at risk of neonatal abstinence syndrome. *Am J Perinatol* 2016 | [26588259](https://pubmed.ncbi.nlm.nih.gov/26588259/) |
| 103 | Rooming-in care for infants of opioid-dependent mothers: implementation and evaluation. 2016 | [27035006](https://pubmed.ncbi.nlm.nih.gov/27035006/) |
| 104 | Rooming-in for infants at risk for neonatal abstinence syndrome: outcomes. 2021 | [33202425](https://pubmed.ncbi.nlm.nih.gov/33202425/) |
| 105 | Lembeck AL et al. Outcome differences in neonates exposed in utero to opioids managed in the NICU versus the pediatric floor. *J Addict Med* 2019 | [30252690](https://pubmed.ncbi.nlm.nih.gov/30252690/) |
| 106 | Jansson LM et al. Methadone maintenance and breastfeeding in the neonatal period. *Pediatrics* 2008 | [18166563](https://pubmed.ncbi.nlm.nih.gov/18166563/) |
| 107 | Jansson LM et al. Methadone maintenance and long-term lactation. 2008 | [18333767](https://pubmed.ncbi.nlm.nih.gov/18333767/) |
| 108 | Begg EJ et al. Distribution of R- and S-methadone into human milk during multiple, medium to high oral dosing. *Br J Clin Pharmacol* 2001 — the 2–9 ng/mL infant plasma concentrations behind `RIDD` | [11736879](https://pubmed.ncbi.nlm.nih.gov/11736879/) |
| 109 | Cantin C et al. Examining the effect of newborn feeding method on postnatal length of stay and healthcare utilisation. *BMJ Paediatr Open* 2026 | [42552081](https://pubmed.ncbi.nlm.nih.gov/42552081/) |
| 110 | Yonke N et al. Breastfeeding motivators and barriers in women receiving medications for opioid use disorder. *Breastfeed Med* 2020 | [31692370](https://pubmed.ncbi.nlm.nih.gov/31692370/) |
| 111 | Short VL et al. A pilot study to assess breastfeeding knowledge, attitudes and perceptions of perinatal staff. *Breastfeed Med* 2019 | [30888210](https://pubmed.ncbi.nlm.nih.gov/30888210/) |
| 112 | Rinaldi K et al. Verbal behavior of mothers with opioid use disorder while feeding infants with neonatal opioid withdrawal syndrome. *Adv Neonatal Care* 2023 | [37011182](https://pubmed.ncbi.nlm.nih.gov/37011182/) |

---

## 11. Gestational age, sex and the preterm phenotype

Constrains `GA50A`, `HGAA` and the preterm isolation experiment.

| # | Reference | PMID |
|---|---|---|
| 113 | Dysart K et al. Sequela of preterm versus term infants born to mothers on a methadone maintenance program: differential course of neonatal abstinence syndrome. *J Perinat Med* 2007 | [17511598](https://pubmed.ncbi.nlm.nih.gov/17511598/) |
| 114 | Anderson VA et al. Sex-related differences in the severity of neonatal opioid withdrawal syndrome. *Am J Perinatol* 2025 | [39074806](https://pubmed.ncbi.nlm.nih.gov/39074806/) |

---

## 12. Polysubstance exposure — the score components an opioid cannot reach

| # | Reference | PMID |
|---|---|---|
| 115 | Winklbaur B et al. Association between prenatal tobacco exposure and outcome of neonates born to opioid-maintained mothers. *Eur Addict Res* 2009 | [19420947](https://pubmed.ncbi.nlm.nih.gov/19420947/) |
| 116 | Wachman EM et al. Impact of psychiatric medication co-exposure on neonatal abstinence syndrome severity. *Drug Alcohol Depend* 2018 | [30205307](https://pubmed.ncbi.nlm.nih.gov/30205307/) |
| 117 | Okoye NC et al. Patterns of neonatal co-exposure to gabapentin and commonly abused drugs. *J Anal Toxicol* 2021 | [32860706](https://pubmed.ncbi.nlm.nih.gov/32860706/) |
| 118 | DeLisle A et al. Gabapentin use during pregnancy and lactation with and without concurrent opioid exposure. *J Addict Med* 2023 | [36069804](https://pubmed.ncbi.nlm.nih.gov/36069804/) |
| 119 | Umer A et al. Gabapentin and opioid co-exposure during pregnancy and adverse perinatal outcomes. *Am J Perinatol* 2026 | [41950956](https://pubmed.ncbi.nlm.nih.gov/41950956/) |
| 120 | Cooper B et al. Association of prenatal fentanyl exposure with neonatal opioid withdrawal syndrome severity. *Hosp Pediatr* 2026 | [42014095](https://pubmed.ncbi.nlm.nih.gov/42014095/) |
| 121 | Sullivan K et al. Outcomes associated with illicit fentanyl exposure during pregnancy. *Hosp Pediatr* 2026 | [42463150](https://pubmed.ncbi.nlm.nih.gov/42463150/) |
| 122 | Illicit fentanyl in the prenatal period: a significant emerging risk for NOWS. 2024 | [39471848](https://pubmed.ncbi.nlm.nih.gov/39471848/) |
| 123 | Hull I et al. A case of maternal and neonatal withdrawal after exposure to fentanyl adulterated with medetomidine. *J Addict Med* 2025 — an α2 agonist in the supply, and therefore in the gap | [41152206](https://pubmed.ncbi.nlm.nih.gov/41152206/) |

---

## 13. Genetics and epigenetics of NOWS severity

| # | Reference | PMID |
|---|---|---|
| 124 | Wachman EM et al. Association of OPRM1 and COMT single-nucleotide polymorphisms with hospital length of stay and treatment of neonatal abstinence syndrome. *JAMA* 2013 | [23632726](https://pubmed.ncbi.nlm.nih.gov/23632726/) |
| 125 | Wachman EM et al. Variations in opioid receptor genes in neonatal abstinence syndrome. *Drug Alcohol Depend* 2015 | [26233486](https://pubmed.ncbi.nlm.nih.gov/26233486/) |
| 126 | Wachman EM et al. Placental OPRM1 DNA methylation and associations with neonatal opioid withdrawal syndrome. 2020 | [33763662](https://pubmed.ncbi.nlm.nih.gov/33763662/) |
| 127 | Camerota M et al. Effects of pharmacologic treatment for neonatal abstinence syndrome on DNA methylation and neurobehavior. *J Pediatr* 2022 | [34971656](https://pubmed.ncbi.nlm.nih.gov/34971656/) |
| 128 | Baldo BA. Neonatal opioid toxicity: opioid withdrawal (abstinence) syndrome with emphasis on pharmacogenomics. *Arch Toxicol* 2023 | [37537419](https://pubmed.ncbi.nlm.nih.gov/37537419/) |

---

## 14. Growth, discharge and longer-term outcome

| # | Reference | PMID |
|---|---|---|
| 129 | Cheng FY et al. Early weight loss percentile curves and feeding practices in opioid-exposed infants. *Hosp Pediatr* 2022 | [36073203](https://pubmed.ncbi.nlm.nih.gov/36073203/) |
| 130 | Incidence of faltering growth in infants with neonatal opioid withdrawal syndrome. 2026 | [42091124](https://pubmed.ncbi.nlm.nih.gov/42091124/) |
| 131 | Early postnatal weight changes in opioid-exposed infants managed using Eat Sleep Console. 2025 | [39973526](https://pubmed.ncbi.nlm.nih.gov/39973526/) |
| 132 | Gaither JR et al. Hospital readmissions among infants with neonatal opioid withdrawal syndrome. *JAMA Netw Open* 2024 | [39316398](https://pubmed.ncbi.nlm.nih.gov/39316398/) |
| 133 | Balalian AA et al. Prenatal exposure to opioids and neurodevelopment in infancy and childhood: a systematic review. *Front Pediatr* 2023 | [36896405](https://pubmed.ncbi.nlm.nih.gov/36896405/) |
| 134 | Rajaprakash M et al. Neurodevelopmental outcomes of prenatal opioid exposure and neonatal opioid withdrawal syndrome: a systematic review. *J Perinatol* 2026 | [41398105](https://pubmed.ncbi.nlm.nih.gov/41398105/) |
| 135 | Prenatal exposure to methadone or buprenorphine: early childhood developmental outcomes. 2018 | [29413437](https://pubmed.ncbi.nlm.nih.gov/29413437/) |
| 136 | Wu Y et al. Antenatal opioid exposure and cerebral cortical maturation in newborns. *JAMA Netw Open* 2026 | [42172030](https://pubmed.ncbi.nlm.nih.gov/42172030/) |
| 137 | Vishnubhotla RV et al. Fetal brain volumes and brain gyrification index associated with opioid exposure. *Brain Commun* 2026 | [42164953](https://pubmed.ncbi.nlm.nih.gov/42164953/) |
| 138 | Mascarenhas M et al. Engagement with recommended developmental follow-up among infants with intrauterine opioid exposure. *J Dev Behav Pediatr* 2025 | [40601929](https://pubmed.ncbi.nlm.nih.gov/40601929/) |

---

## 15. Iatrogenic withdrawal elsewhere — the same physiology, a different patient

Included because the re-inducible adaptation pool `AT` is the same object that
produces iatrogenic withdrawal in the paediatric ICU, and that literature is
where its time constant is best observed.

| # | Reference | PMID |
|---|---|---|
| 139 | Kim CS et al. Evaluation of dexmedetomidine withdrawal and management after prolonged infusion. *Clin Ther* 2024 — α2 tachyphylaxis and rebound, behind `KTACH`/`KTOFF` | [39379223](https://pubmed.ncbi.nlm.nih.gov/39379223/) |
| 140 | Playfor S et al. Sedation in critically ill children. *J Clin Med* 2025 | [40944032](https://pubmed.ncbi.nlm.nih.gov/40944032/) |
| 141 | Zhan Y et al. Correlation between analgesic and sedative drug withdrawal and salivary cortisol levels in children in the PICU. *Pharmacotherapy* 2025 | [41251283](https://pubmed.ncbi.nlm.nih.gov/41251283/) |
| 142 | Men S, Wang H. Phenobarbital in nuclear receptor activation: an update. *Drug Metab Dispos* 2023 — the CYP/UGT induction arm of `PHENO` | [36351837](https://pubmed.ncbi.nlm.nih.gov/36351837/) |
| 143 | Tien YC et al. Phenobarbital treatment at a neonatal age results in decreased efficacy of omeprazole in adult mice. *Drug Metab Dispos* 2017 | [28062542](https://pubmed.ncbi.nlm.nih.gov/28062542/) |

---

## 16. Modelling tools

- **mrgsolve** — <https://mrgsolve.org/>
- **QSP in R with mrgsolve** — <https://vantage-research.net/qsp-in-r/>
- **gPKPDviz** (mrgsolve-based Shiny PK/PD simulation): paper
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/> · code
  <https://github.com/Genentech/gPKPDviz/>
- PubMed E-utilities (used to verify every PMID above) —
  <https://www.ncbi.nlm.nih.gov/books/NBK25501/>

---

## A note on what these references do and do not support

The pharmacokinetic, receptor-pharmacology and maturation parameters above are
taken from the cited measurements. The **neuroadaptation block is not
measured in human neonates**. Its structure (a durable pool set in utero plus a
re-inducible pool driven by current exposure) is inferred from the rodent locus
coeruleus literature (refs 62–66, 71–75) and from the clinical time course, and
its two time constants — 12 days for the durable pool, 60 hours for the
re-inducible one — are the model's most exposed assumptions. They are
identifiable in principle from serial scores under different weaning rules, and
`README.md` states what data would falsify them.
