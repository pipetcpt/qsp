# Takotsubo Syndrome (TTS) — Reference List

References for `tts_qsp_model.dot`, `tts_mrgsolve_model.R` and `tts_shiny_app.R`.

**How this list was built.** Every PMID below was resolved by querying the NCBI
E-utilities PubMed API during the session that produced this model, and the
title / journal / year printed here are the ones PubMed returned for that PMID.
Nothing was written from memory. Where a query returned no defensible match the
claim was left uncited in the model files rather than attached to a guessed
identifier — the reverse is how a bibliography silently fills with PMIDs
pointing at unrelated papers.

**125 references, 15 sections.**

A note on what these citations do and do not support. The model's central
structural claim — that apical ballooning is the apex crossing an
agonist-driven Gs-to-Gi threshold first because its beta2-adrenoceptor density
is highest — rests on section 2, and in particular on the experimental
demonstration that adrenaline (but not noradrenaline) recruits Gi through
beta2-AR and that this produces apical, not basal, cardiodepression. The
quantitative apex-to-base receptor density ratio used in the model (1.40, with
beta2 fraction 0.42 apically vs 0.24 basally) is a *modelling assumption
consistent with* the innervation and receptor literature in section 4, not a
directly measured human value; treat it as the model's main falsifiable input.

---

### 1. Definition, diagnostic criteria, registries and epidemiology

1. Takotsubo (Stress) Cardiomyopathy. *N Engl J Med* 2015. PMID [26716925](https://pubmed.ncbi.nlm.nih.gov/26716925/)
2. International Expert Consensus Document on Takotsubo Syndrome (Part I): Clinical Characteristics, Diagnostic Criteria, and Pathophysiology. *Eur Heart J* 2018. PMID [29850871](https://pubmed.ncbi.nlm.nih.gov/29850871/)
3. International Expert Consensus Document on Takotsubo Syndrome (Part II): Diagnostic Workup, Outcome, and Management. *Eur Heart J* 2018. PMID [29850820](https://pubmed.ncbi.nlm.nih.gov/29850820/)
4. Pathophysiology of Takotsubo syndrome - a joint scientific statement from the Heart Failure Association Takotsubo Syndrome Study Group and Myocardial Function Working Group of the European Society of Cardiology - Part 1: overview and the central role for catecholamines and sympathetic nervous system. *Eur J Heart Fail* 2022. PMID [34907620](https://pubmed.ncbi.nlm.nih.gov/34907620/)
5. Age-Related Variations in Takotsubo Syndrome. *J Am Coll Cardiol* 2020. PMID [32327096](https://pubmed.ncbi.nlm.nih.gov/32327096/)
6. Apical ballooning syndrome or takotsubo cardiomyopathy: a systematic review. *Eur Heart J* 2006. PMID [16720686](https://pubmed.ncbi.nlm.nih.gov/16720686/)
7. Happy heart syndrome: role of positive emotional stress in takotsubo syndrome. *Eur Heart J* 2016. PMID [26935270](https://pubmed.ncbi.nlm.nih.gov/26935270/)
8. Hemodynamic Vulnerability in Male Patients with Takotsubo Syndrome: Pathophysiological Roles of Tachycardia. *Mayo Clin Proc* 2026. PMID [42520964](https://pubmed.ncbi.nlm.nih.gov/42520964/)
9. Non-ACS Presentation in Takotsubo Syndrome: A High-Risk Phenotype with Distinct Treatment Patterns and Worse Outcomes. *JACC Asia* 2026. PMID [42496624](https://pubmed.ncbi.nlm.nih.gov/42496624/)
10. Stress cardiomyopathy: clinical and ventriculographic characteristics in 107 North American subjects. *Int J Cardiol* 2010. PMID [19155079](https://pubmed.ncbi.nlm.nih.gov/19155079/)
11. Long-Term Prognosis, Risk Assessment, and Management of Patients Diagnosed with Takotsubo Syndrome: A Narrative Review. *J Pers Med* 2025. PMID [41003128](https://pubmed.ncbi.nlm.nih.gov/41003128/)
12. Trigger-Negative Takotsubo Syndrome: A Distinct Phenotype Prone to Recurrence. *Circ Rep* 2026. PMID [42110870](https://pubmed.ncbi.nlm.nih.gov/42110870/)

### 2. Catecholamines and the beta2-adrenoceptor Gs-to-Gi switch (the model's core premise)

13. Stress cardiomyopathy: a syndrome of catecholamine-mediated myocardial stunning?. *Cell Mol Neurobiol* 2012. PMID [22297544](https://pubmed.ncbi.nlm.nih.gov/22297544/)
14. High levels of circulating epinephrine trigger apical cardiodepression in a β2-adrenergic receptor/Gi-dependent manner: a new model of Takotsubo cardiomyopathy. *Circulation* 2012. PMID [22732314](https://pubmed.ncbi.nlm.nih.gov/22732314/)
15. Epinephrine regulation of hemodynamics in catecholamine knockouts and the pithed mouse. *Ann N Y Acad Sci* 2008. PMID [19120125](https://pubmed.ncbi.nlm.nih.gov/19120125/)
16. Epinephrine activates both Gs and Gi pathways, but norepinephrine activates only the Gs pathway through human beta2-adrenoceptors overexpressed in mouse heart. *Mol Pharmacol* 2004. PMID [15102960](https://pubmed.ncbi.nlm.nih.gov/15102960/)
17. Switching of the coupling of the beta2-adrenergic receptor to different G proteins by protein kinase A. *Nature* 1997. PMID [9363896](https://pubmed.ncbi.nlm.nih.gov/9363896/)
18. Pertussis toxin suppresses carbachol-evoked cardiodepression but does not modify cardiostimulation mediated through beta1- and putative beta4-adrenoceptors in mouse left atria: no evidence for beta2- and beta3-adrenoreceptor function. *Naunyn Schmiedebergs Arch Pharmacol* 2000. PMID [10685868](https://pubmed.ncbi.nlm.nih.gov/10685868/)
19. Ligand-directed signalling at beta-adrenoceptors. *Br J Pharmacol* 2010. PMID [20132209](https://pubmed.ncbi.nlm.nih.gov/20132209/)
20. Roles of PKA, PI3K, and cPLA2 in the NO-mediated negative inotropic effect of beta2-adrenoceptor agonists in guinea pig right papillary muscles. *Am J Physiol Cell Physiol* 2008. PMID [17942637](https://pubmed.ncbi.nlm.nih.gov/17942637/)
21. Pathophysiological mechanisms of catecholamine and cocaine-mediated cardiotoxicity. *Heart Fail Rev* 2014. PMID [24398587](https://pubmed.ncbi.nlm.nih.gov/24398587/)
22. Pathophysiology of Takotsubo Syndrome. *Circulation* 2017. PMID [28606950](https://pubmed.ncbi.nlm.nih.gov/28606950/)
23. Catecholamine-induced acute myocardial stunning after accidental intra-operative noradrenaline bolus. *Anaesth Rep* 2022. PMID [36246420](https://pubmed.ncbi.nlm.nih.gov/36246420/)
24. Adrenaline, Takotsubo Cardiomyopathy, Anaphylaxis, and Kounis Syndrome (ATAK) Complex: Clinical Phenotypes, Differential Diagnosis, and Management. *JACC Case Rep* 2026. PMID [42132731](https://pubmed.ncbi.nlm.nih.gov/42132731/)
25. miR-181a post-transcriptionally targets GRK2 to limit maladaptive signaling in cardiomyocytes. *Front Cardiovasc Med* 2026. PMID [42099781](https://pubmed.ncbi.nlm.nih.gov/42099781/)
26. Non-canonical signaling initiated by hormone-responsive G protein-coupled receptors from subcellular compartments. *Pharmacol Ther* 2025. PMID [39722422](https://pubmed.ncbi.nlm.nih.gov/39722422/)
27. Beta1- and beta2-adrenoceptor polymorphisms and cardiovascular diseases. *Fundam Clin Pharmacol* 2008. PMID [18353108](https://pubmed.ncbi.nlm.nih.gov/18353108/)
28. Adenylyl cyclase isoforms 5 and 6 in the cardiovascular system: complex regulation and divergent roles. *Front Pharmacol* 2024. PMID [38633617](https://pubmed.ncbi.nlm.nih.gov/38633617/)

### 3. Brain-heart axis, sex and oestradiol

29. Altered limbic and autonomic processing supports brain-heart axis in Takotsubo syndrome. *Eur Heart J* 2019. PMID [30831580](https://pubmed.ncbi.nlm.nih.gov/30831580/)
30. Stress-induced takotsubo syndrome: dynamic changes in regional cerebral metabolism revealed by quantitative PET imaging. *Neuroimage* 2026. PMID [41679568](https://pubmed.ncbi.nlm.nih.gov/41679568/)
31. Estradiol mitigates stress-induced cardiac injury and inflammation by downregulating ADAM17 via the GPER-1/PI3K signaling pathway. *Cell Mol Life Sci* 2023. PMID [37572114](https://pubmed.ncbi.nlm.nih.gov/37572114/)
32. Cardiac biomarkers combined with neuroimaging localization predict long-term outcomes in acute stroke patients. *Sci Rep* 2026. PMID [42260095](https://pubmed.ncbi.nlm.nih.gov/42260095/)
33. Autonomic Imbalance in Cardiomyopathy and Heart Failure: From Neurobiology to Precision Neuromodulation. *Curr Cardiol Rep* 2025. PMID [41117884](https://pubmed.ncbi.nlm.nih.gov/41117884/)

### 4. Regional receptor field, innervation and the apex-to-base gradient

34. Assessment of takotsubo (ampulla) cardiomyopathy using iodine-123 metaiodobenzylguanidine scintigraphy. *Acta Radiol* 2006. PMID [17135004](https://pubmed.ncbi.nlm.nih.gov/17135004/)
35. Cardiac Sympathetic Positron Emission Tomography Imaging with Meta-[(18)F]Fluorobenzylguanidine is Sensitive to Uptake-1 in Rats. *ACS Chem Neurosci* 2021. PMID [34714061](https://pubmed.ncbi.nlm.nih.gov/34714061/)
36. Chronic nicotine exposure is associated with electrophysiological and sympathetic remodeling in the intact rabbit heart. *Am J Physiol Heart Circ Physiol* 2024. PMID [38551482](https://pubmed.ncbi.nlm.nih.gov/38551482/)

### 5. Imaging: oedema, strain, energetics and the perfusion-contraction mismatch

37. Clinical characteristics and cardiovascular magnetic resonance findings in stress (takotsubo) cardiomyopathy. *JAMA* 2011. PMID [21771988](https://pubmed.ncbi.nlm.nih.gov/21771988/)
38. CMR Findings Across Disease Phases in Takotsubo Syndrome: Insights From the Multicenter EVOLUTION Registry. *Am J Cardiol* 2026. PMID [42086089](https://pubmed.ncbi.nlm.nih.gov/42086089/)
39. Response to Letters Regarding Article, "Persistent Long-Term Structural, Functional, and Metabolic Changes After Stress-Induced (Takotsubo) Cardiomyopathy". *Circulation* 2018. PMID [30354449](https://pubmed.ncbi.nlm.nih.gov/30354449/)
40. Response by Scally and Dawson to Letters Regarding Article, "Myocardial and Systemic Inflammation in Acute Stress-Induced (Takotsubo) Cardiomyopathy". *Circulation* 2019. PMID [31545684](https://pubmed.ncbi.nlm.nih.gov/31545684/)
41. CMR Reveals the Influence of Trigger and Classification on the Myocardial Tissue Response in Takotsubo Syndrome. *Circ Cardiovasc Imaging* 2026. PMID [42381632](https://pubmed.ncbi.nlm.nih.gov/42381632/)
42. Fluorodeoxyglucose positron emission tomography/computed tomography in Takotsubo cardiomyopathy complicated by sinus arrest and third-degree atrioventricular block: a case report. *Eur Heart J Case Rep* 2026. PMID [41669316](https://pubmed.ncbi.nlm.nih.gov/41669316/)
43. Reversibility of myocardial oedema and regional wall motion abnormalities in a non-Takotsubo-like pattern in a lightning strike survivor detected by cardiac MRI: a case report. *Eur Heart J Case Rep* 2026. PMID [41978764](https://pubmed.ncbi.nlm.nih.gov/41978764/)
44. Right Ventricular Impairment Prevalence in Takotsubo Syndrome and Associated Clinical Characteristics and Outcomes: EVOLUTION Registry Results. *Radiol Cardiothorac Imaging* 2026. PMID [41989283](https://pubmed.ncbi.nlm.nih.gov/41989283/)
45. Persistent left atrial mechanical dysfunction and long-term risk in Takotsubo syndrome. *Int J Cardiol* 2026. PMID [42035838](https://pubmed.ncbi.nlm.nih.gov/42035838/)
46. Mid-Ventricular Takotsubo Cardiomyopathy. ** 2026. PMID [32491438](https://pubmed.ncbi.nlm.nih.gov/32491438/)
47. Atypical Pediatric Presentation of Takotsubo Syndrome. *JACC Case Rep* 2026. PMID [42484563](https://pubmed.ncbi.nlm.nih.gov/42484563/)

### 6. Biomarkers: the small-troponin / large-NT-proBNP discordance

48. The Role of Cardiac Biomarkers in Evaluating Takotsubo Cardiomyopathy: A Systematic Review. *Cureus* 2025. PMID [40677457](https://pubmed.ncbi.nlm.nih.gov/40677457/)
49. Cardiac Troponin release, myocardial function and inflammation in patients with takotsubo syndrome: a cardiac magnetic resonance study. *Int J Cardiovasc Imaging* 2026. PMID [42370982](https://pubmed.ncbi.nlm.nih.gov/42370982/)
50. Reply: Potential Mechanisms of Troponin Release in Stable Coronary Artery Disease. *J Am Coll Cardiol* 2023. PMID [37993208](https://pubmed.ncbi.nlm.nih.gov/37993208/)
51. Left ventricular function recovery in Takotsubo syndrome, clinical and pathophysiological insights: a state-of-the-art review. *Curr Opin Cardiol* 2025. PMID [40900430](https://pubmed.ncbi.nlm.nih.gov/40900430/)
52. Takotsubo Syndrome Mimicking Acute Myocardial Infarction: A Case Report. *Cureus* 2026. PMID [42022711](https://pubmed.ncbi.nlm.nih.gov/42022711/)
53. Correlation of Right Ventricular Wall Stress With Plasma B-Type Natriuretic Peptide Levels in Patients With Pulmonary Hypertension. *Circ J* 2019. PMID [30971626](https://pubmed.ncbi.nlm.nih.gov/30971626/)

### 7. Coronary microvascular function

54. Transient left ventricular apical ballooning without coronary artery stenosis: a novel heart syndrome mimicking acute myocardial infarction. Angina Pectoris-Myocardial Infarction Investigations in Japan. *J Am Coll Cardiol* 2001. PMID [11451258](https://pubmed.ncbi.nlm.nih.gov/11451258/)
55. Correlation Between Index of Microcirculatory Resistance and Angiography-Derived Microcirculatory Resistance in Takotsubo Syndrome. *Catheter Cardiovasc Interv* 2026. PMID [41332115](https://pubmed.ncbi.nlm.nih.gov/41332115/)
56. Vasospasm-Induced Takotsubo Cardiomyopathy: An Underrecognized Phenotype of Ischemia With Nonobstructive Coronary Arteries. *Cureus* 2026. PMID [42255814](https://pubmed.ncbi.nlm.nih.gov/42255814/)
57. Coronary Microvascular Dysfunction in Stress Cardiomyopathy: At the Heart of the Problem. *Life (Basel)* 2026. PMID [42514160](https://pubmed.ncbi.nlm.nih.gov/42514160/)
58. Coronary Flow Reserve in Adults: Pathophysiology, Assessment Modalities, Clinical Applications, and Prognostic Significance. *Medicina (Kaunas)* 2026. PMID [42356048](https://pubmed.ncbi.nlm.nih.gov/42356048/)
59. Coronary microvascular dysfunction in menopausal women. *Heart* 2026. PMID [42331611](https://pubmed.ncbi.nlm.nih.gov/42331611/)
60. Evaluation of non-invasive imaging parameters in coronary microvascular disease: a systematic review. *BMC Med Imaging* 2021. PMID [33407208](https://pubmed.ncbi.nlm.nih.gov/33407208/)
61. A microstructurally motivated framework to study autoregulation in the coronary circulation. *J Physiol* 2026. PMID [42303295](https://pubmed.ncbi.nlm.nih.gov/42303295/)

### 8. Dynamic LVOT obstruction, cardiogenic shock and mechanical support

62. Cardiogenic shock in Takotsubo syndrome: Insights into phenotype-tailored management. *Med Intensiva (Engl Ed)* 2026. PMID [42251021](https://pubmed.ncbi.nlm.nih.gov/42251021/)
63. Timing of cardiogenic shock and clinical outcomes in takotsubo syndrome: A multicenter cohort study. *Int J Cardiol* 2026. PMID [41765143](https://pubmed.ncbi.nlm.nih.gov/41765143/)
64. Left ventricular outflow tract obstruction in critically ill patients: from pathophysiology and diagnosis to the management with the "LVOTO" bundle. *Anaesth Crit Care Pain Med* 2026. PMID [41713684](https://pubmed.ncbi.nlm.nih.gov/41713684/)
65. Takotsubo syndrome complicated by dynamic left ventricular outflow tract obstruction after left bundle branch area pacing. *HeartRhythm Case Rep* 2026. PMID [42472203](https://pubmed.ncbi.nlm.nih.gov/42472203/)
66. Management of cardiogenic shock by circulatory support during reverse Tako-Tsubo following amphetamine exposure: A report of two cases. *Heart Lung* 2021. PMID [33243478](https://pubmed.ncbi.nlm.nih.gov/33243478/)
67. Successful Use of Impella Support to Treat Cardiogenic Shock Secondary to Takotsubo Cardiomyopathy Owing to Alcohol Withdrawal. *Cureus* 2025. PMID [40951064](https://pubmed.ncbi.nlm.nih.gov/40951064/)
68. Takotsubo syndrome treated with VA-ECMO in plastic surgery: echocardiography first-case series. *Eur Heart J Case Rep* 2026. PMID [42428664](https://pubmed.ncbi.nlm.nih.gov/42428664/)
69. Shock and Awe: The Tactical Trade-Offs of Impella(®) Versus Intra-Aortic Balloon Pump in Takotsubo Cardiomyopathy. *Reports (MDPI)* 2025. PMID [40710834](https://pubmed.ncbi.nlm.nih.gov/40710834/)
70. Could Esmolol/Landiolol Infusion During Pulsed Field Ablation of Atrial Fibrillation Have Averted Takotsubo Syndrome?. *JACC Case Rep* 2026. PMID [41954320](https://pubmed.ncbi.nlm.nih.gov/41954320/)

### 9. Electrophysiology, QT and arrhythmia

71. Takotsubo Syndrome and Sudden Cardiac Death. *Angiology* 2023. PMID [35668627](https://pubmed.ncbi.nlm.nih.gov/35668627/)
72. Ventricular arrhythmias in patients with Takotsubo syndrome. *J Arrhythm* 2018. PMID [30167007](https://pubmed.ncbi.nlm.nih.gov/30167007/)
73. Electrocardiographic Findings in Takotsubo Cardiomyopathy: ECG Evolution and Its Difference from the ECG of Acute Coronary Syndrome. *Clin Med Insights Cardiol* 2014. PMID [24653650](https://pubmed.ncbi.nlm.nih.gov/24653650/)
74. Malignant arrhythmia in apical ballooning syndrome: risk factors and outcomes. *Indian Pacing Electrophysiol J* 2008. PMID [18679529](https://pubmed.ncbi.nlm.nih.gov/18679529/)
75. QTc Prolongation in Stress-Induced Cardiomyopathy: A Case of Stabilization With Temporary High-Rate Pacing. *JACC Case Rep* 2026. PMID [42207064](https://pubmed.ncbi.nlm.nih.gov/42207064/)
76. An unexpected sequel to happiness: happy heart syndrome. *Oxf Med Case Reports* 2026. PMID [42422269](https://pubmed.ncbi.nlm.nih.gov/42422269/)

### 10. Left ventricular thrombus and embolism

77. Left Ventricular Thrombi in Takotsubo Syndrome: Incidence, Predictors, and Management: Results From the GEIST (German Italian Stress Cardiomyopathy) Registry. *J Am Heart Assoc* 2017. PMID [29203578](https://pubmed.ncbi.nlm.nih.gov/29203578/)
78. Intraventricular thrombus in Takotsubo syndrome: Incidence, predictors, management, and prognosis. Insights from the RETAKO registry. *Int J Cardiol* 2025. PMID [39826576](https://pubmed.ncbi.nlm.nih.gov/39826576/)
79. Comparison of treatment outcomes of direct oral anticoagulants and heparin for patients with Takotsubo cardiomyopathy: A nationwide cohort analysis. *PLoS One* 2025. PMID [41231880](https://pubmed.ncbi.nlm.nih.gov/41231880/)

### 11. Treatment, outcome and recurrence

80. Stress Cardiomyopathy Diagnosis and Treatment: JACC State-of-the-Art Review. *J Am Coll Cardiol* 2018. PMID [30309474](https://pubmed.ncbi.nlm.nih.gov/30309474/)
81. Survival benefit of secondary prevention medical therapy in takotsubo cardiomyopathy: a Bayesian network meta-analysis. *Eur Heart J Open* 2025. PMID [40357262](https://pubmed.ncbi.nlm.nih.gov/40357262/)
82. SGLT2 inhibitors are associated with improved long-term survival in Takotsubo syndrome: insights from large-scale real-world data. *Eur Heart J Cardiovasc Pharmacother* 2026. PMID [41549637](https://pubmed.ncbi.nlm.nih.gov/41549637/)
83. Prognostic relevance of GRACE risk score in Takotsubo syndrome. *Eur Heart J Acute Cardiovasc Care* 2020. PMID [31642689](https://pubmed.ncbi.nlm.nih.gov/31642689/)
84. Long-term prognosis in Takotsubo syndrome compared to heart failure: observations from a global federated research network. *ESC Heart Fail* 2026. PMID [41761829](https://pubmed.ncbi.nlm.nih.gov/41761829/)
85. Serotonin Norepinephrine Reuptake Inhibitor Is Associated With Lower Mortality Among Patients Presenting With Takotsubo Cardiomyopathy. *J Am Heart Assoc* 2025. PMID [40551298](https://pubmed.ncbi.nlm.nih.gov/40551298/)
86. Cancer-Associated Takotsubo Syndrome: An Emerging Intersection of Cardio-Oncology. *Cardiol Rev* 2026. PMID [41941451](https://pubmed.ncbi.nlm.nih.gov/41941451/)
87. Beyond the Apex: A Case Series of Mid-ventricular Takotsubo Cardiomyopathy. *Cureus* 2026. PMID [41737107](https://pubmed.ncbi.nlm.nih.gov/41737107/)
88. Recurrent Takotsubo Cardiomyopathy Triggered by COVID-19 Infection Complicated by Ventricular Tachycardia Arrest and Cardiogenic Shock: A Case Report. *Cureus* 2026. PMID [42273520](https://pubmed.ncbi.nlm.nih.gov/42273520/)
89. Retrospective Diagnosis of Recurrent Takotsubo Syndrome Episodes in a 40-Year-Old Woman. *Diagnostics (Basel)* 2026. PMID [42449864](https://pubmed.ncbi.nlm.nih.gov/42449864/)

### 12. Drug pharmacokinetics and pharmacodynamics used in the model

90. Population Pharmacokinetics of Levosimendan and its Metabolites OR-1855 and OR-1896 in Critically Ill Adults, Neonates and Infants on Veno-Arterial ECMO. *Clin Pharmacokinet* 2026. PMID [41288922](https://pubmed.ncbi.nlm.nih.gov/41288922/)
91. Drug sequestration and metabolite formation: key pharmacokinetic challenges for levosimendan use during ECMO support in cardiogenic shock. *J Intensive Care* 2026. PMID [42032778](https://pubmed.ncbi.nlm.nih.gov/42032778/)
92. Pharmacokinetics of levosimendan in critically Ill children on extracorporeal membrane oxygenation: a prospective observational study. *Front Pediatr* 2025. PMID [40809383](https://pubmed.ncbi.nlm.nih.gov/40809383/)
93. Use of Levosimendan in Intensive Care Unit Settings: An Opinion Paper. *J Cardiovasc Pharmacol* 2019. PMID [30489437](https://pubmed.ncbi.nlm.nih.gov/30489437/)
94. Targeting cardiac myosin in HFrEF: mechanism, clinical evidence, and the role of omecamtiv mecarbil. *J Cardiothorac Surg* 2026. PMID [42289718](https://pubmed.ncbi.nlm.nih.gov/42289718/)
95. Population pharmacokinetics and pharmacodynamics of dobutamine in neonates on the first days of life. *Br J Clin Pharmacol* 2020. PMID [31657867](https://pubmed.ncbi.nlm.nih.gov/31657867/)
96. Time-Varying Clearance in Milrinone Pharmacokinetics from Premature Neonates to Adolescents. *Clin Pharmacokinet* 2024. PMID [38613610](https://pubmed.ncbi.nlm.nih.gov/38613610/)
97. Novel Nebulized Milrinone Formulation for the Treatment of Acute Heart Failure Requiring Inotropic Therapy: A Phase 1 Study. *J Card Fail* 2024. PMID [37871843](https://pubmed.ncbi.nlm.nih.gov/37871843/)
98. Metoprolol Population Pharmacokinetics in Older Chinese Patients With CKM Syndrome: Joint Effects of rs1065852 and CKM(2)S(2)-BAG Score on Clearance. *Clin Transl Sci* 2026. PMID [42496093](https://pubmed.ncbi.nlm.nih.gov/42496093/)
99. Pharmacokinetics and Pharmacodynamics of a Novel Formulation of Furosemide Administered as a Single Subcutaneous Injection. *Clin Ther* 2025. PMID [41093674](https://pubmed.ncbi.nlm.nih.gov/41093674/)
100. Intravenous Bolus Phenylephrine and Intravenous Bolus Norepinephrine for Treatment of Maternal Hypotension in Spinal Anesthesia During Cesarean Section: A Prospective, Randomized, Comparative Study. *Cureus* 2026. PMID [41835774](https://pubmed.ncbi.nlm.nih.gov/41835774/)
101. An Overview of the Pharmacokinetics and Pharmacodynamics of Landiolol (an Ultra-Short Acting β1 Selective Antagonist) in Atrial Fibrillation. *Pharmaceutics* 2024. PMID [38675178](https://pubmed.ncbi.nlm.nih.gov/38675178/)
102. A multicentre observational study on landiolol use, efficacy, and safety in European patients with supraventricular arrhythmia (Landi-UP). *Eur Heart J Acute Cardiovasc Care* 2026. PMID [41052283](https://pubmed.ncbi.nlm.nih.gov/41052283/)

### 13. Signalling, calcium handling, myofilaments and stunning

103. Differential changes in cyclic adenosine 3'-5' monophosphate (cAMP) effectors and major Ca(2+) handling proteins during diabetic cardiomyopathy. *J Cell Mol Med* 2023. PMID [36967707](https://pubmed.ncbi.nlm.nih.gov/36967707/)
104. Phosphoproteomics distinguishes disease-specific mechanisms for human phospholamban cardiomyopathy reversible by RNA therapy. *Signal Transduct Target Ther* 2026. PMID [42204136](https://pubmed.ncbi.nlm.nih.gov/42204136/)
105. Epac1 increases myosin regulatory light-chain phosphorylation, energetic cost of contraction, and susceptibility to heart failure. *PLoS One* 2025. PMID [40526593](https://pubmed.ncbi.nlm.nih.gov/40526593/)
106. Proteolytic degradation of atrial sarcomere proteins underlies contractile defects in atrial fibrillation. *Am J Physiol Heart Circ Physiol* 2024. PMID [38940916](https://pubmed.ncbi.nlm.nih.gov/38940916/)
107. Myocardial stunning: mechanisms, molecular insights, and gaps in knowledge. *Biosci Rep* 2025. PMID [41400626](https://pubmed.ncbi.nlm.nih.gov/41400626/)
108. Delayed adaptation of the heart to stress: late preconditioning. *Stroke* 2004. PMID [15459441](https://pubmed.ncbi.nlm.nih.gov/15459441/)

### 14. Differential diagnosis and mimics

109. Diagnostic Challenges in Takotsubo Syndrome: Bridging Mimics, Mechanisms, and Management. *J Clin Med* 2026. PMID [42452549](https://pubmed.ncbi.nlm.nih.gov/42452549/)
110. The Role of Serum Biomarkers for the Differential Diagnosis and Prognostic Assessment of Myocardial Infarction with Non-Obstructive Coronary Arteries: A Narrative Review. *J Clin Med* 2026. PMID [41976894](https://pubmed.ncbi.nlm.nih.gov/41976894/)
111. When Myocarditis Masquerades as ST-Elevation Myocardial Infarction: A Case of Coxsackie B-induced Acute Heart Failure With Rapid Recovery. *Cureus* 2025. PMID [41552188](https://pubmed.ncbi.nlm.nih.gov/41552188/)
112. Pressor use and its impact on outcomes in aneurysmal subarachnoid hemorrhage patients with takotsubo cardiomyopathy: a quantitative analysis. *Acta Neurochir (Wien)* 2026. PMID [41701379](https://pubmed.ncbi.nlm.nih.gov/41701379/)
113. Successful Management of Subarachnoid Hemorrhage Complicated by Takotsubo Cardiomyopathy Using Distal Transradial Access (dTRA) Coiling and Integrated Pharmacotherapy Under Intra-aortic Balloon Pumping (IABP) Support: A Stroke-Heart Syndrome Case. *Cureus* 2025. PMID [40671977](https://pubmed.ncbi.nlm.nih.gov/40671977/)
114. Neurogenic shock to the heart: a rare case of meningitis-triggered reverse takotsubo cardiomyopathy. *Proc (Bayl Univ Med Cent)* 2026. PMID [42269061](https://pubmed.ncbi.nlm.nih.gov/42269061/)
115. A Case of Peripartum Cardiogenic Shock Resulting From Reverse Takotsubo Cardiomyopathy. *J Med Cases* 2026. PMID [42327864](https://pubmed.ncbi.nlm.nih.gov/42327864/)
116. Severe Biventricular Takotsubo Syndrome Triggered by Status Epilepticus With Pulmonary Embolism and Critical Illness Myopathy. *JACC Case Rep* 2026. PMID [42489612](https://pubmed.ncbi.nlm.nih.gov/42489612/)
117. Inverted Takotsubo Cardiomyopathy Triggered by Corticosteroid Administration in a Patient with Undiagnosed Pheochromocytoma. *J Investig Allergol Clin Immunol* 2026. PMID [42522684](https://pubmed.ncbi.nlm.nih.gov/42522684/)
118. Reverse Takotsubo Syndrome in Women During High Hormonal States Related to In Vitro Fertilization. *JACC Case Rep* 2026. PMID [42138668](https://pubmed.ncbi.nlm.nih.gov/42138668/)
119. Cardiomyocyte-enriched OTUD5 alleviates septic cardiomyopathy by promoting NLRP3 deubiquitination and inhibiting NLRP3 inflammasome activation. *Clin Transl Med* 2026. PMID [42437953](https://pubmed.ncbi.nlm.nih.gov/42437953/)
120. Oxidative Stress in Takotsubo Syndrome: Insights into Extracellular Vesicles and Their Potential Clinical Relevance. *Antioxidants (Basel)* 2026. PMID [41897449](https://pubmed.ncbi.nlm.nih.gov/41897449/)
121. Biventricular takotsubo syndrome complicated with cardiogenic shock and shark fin sign requiring ECPELLA: a case report. *Eur Heart J Case Rep* 2025. PMID [40800555](https://pubmed.ncbi.nlm.nih.gov/40800555/)

### 15. Quantitative systems pharmacology and modelling methodology

122. Development of a Quantitative Systems Pharmacology Model to Interrogate Mitochondrial Metabolism in Heart Failure. *bioRxiv* 2025. PMID [40672261](https://pubmed.ncbi.nlm.nih.gov/40672261/)
123. Development of a population pharmacokinetics model and an R-Shiny simulation platform for moxifloxacin pharmacokinetics in non-human primates. *J Pharmacokinet Pharmacodyn* 2026. PMID [41920374](https://pubmed.ncbi.nlm.nih.gov/41920374/)
124. OpenPMX Software for Nonlinear Mixed-Effect Models in Pharmacometrics: Precision Compared With NONMEM First-Order Conditional Estimation. *CPT Pharmacometrics Syst Pharmacol* 2026. PMID [42166222](https://pubmed.ncbi.nlm.nih.gov/42166222/)
125. Stress-strain phenotyping, reverse remodelling and therapeutic responsiveness in HFrEF: Predictive value of end-systolic wall stress and global longitudinal strain. *Int J Cardiol* 2026. PMID [42269879](https://pubmed.ncbi.nlm.nih.gov/42269879/)

---

## Not cited, and why

* No PMID is attached to the specific numerical values of `RHO_AP`, `FB2_AP`,
  `INN_AP` or `TH_AP`. They are the model's structural inputs and are chosen to
  be *consistent with* the innervation and receptor literature above, not read
  off a single paper. Section 4 supports the existence and direction of the
  gradients; it does not fix their magnitude.
* Reverse (basal) takotsubo is represented in the model only by inverting those
  same parameters. The papers in section 14 establish that the phenotype
  exists; none of them measures an inverted receptor gradient. The model
  therefore states the inverted gradient as a *prediction*, not a citation.
* The oedema time constant was set from the QTc / T-wave time course (peak on
  days 2-3) rather than from CMR T2 normalisation (weeks to months). Both are in
  sections 5 and 9; they disagree about the timescale, and the model follows the
  electrophysiological one. This is recorded as a limitation in the README, not
  resolved by citation.

## Licence

Educational and research use only. See the repository `LICENSE`. These models
are not validated for clinical or regulatory use.
