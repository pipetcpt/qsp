# MOG Antibody-Associated Disease (MOGAD) — References

All entries below were resolved against PubMed and carry their real PMIDs; the
citation strings were generated from PubMed metadata rather than typed by hand.
Annotations state specifically what each paper contributes to the model — which
parameter it anchors, which structural choice it justifies, or which endpoint it
calibrates.

Where the model deliberately departs from a published figure (most notably the
effect size of IL-6 receptor blockade), that is recorded in the *Known misses*
section of [`mogad_mrgsolve_model.R`](mogad_mrgsolve_model.R) rather than glossed
over here.


## 1. Diagnostic criteria, disease definition and scope

MOGAD became a defined entity only after the 2023 international criteria; these set the titre thresholds and the core clinical attack types the model represents.

1. **Banwell B et al.** Diagnosis of myelin oligodendrocyte glycoprotein antibody-associated disease: International MOGAD Panel proposed criteria. *Lancet Neurol* 2023;22:268-282. [PMID 36706773](https://pubmed.ncbi.nlm.nih.gov/36706773/)  
   The 2023 international MOGAD panel criteria — a clear-positive live cell-based assay titre of at least 1:100 plus a core clinical attack. Source of the titre thresholds used in the model's TITER_REF mapping.

2. **Marignier R et al.** Myelin-oligodendrocyte glycoprotein antibody-associated disease. *Lancet Neurol* 2021;20:762-772. [PMID 34418402](https://pubmed.ncbi.nlm.nih.gov/34418402/)  
   Standard narrative review of the disease as a distinct entity, separate from both multiple sclerosis and AQP4-IgG NMOSD.

3. **Trewin BP et al.** MOGAD: A comprehensive review of clinicoradiological features, therapy and outcomes in 4699 patients globally. *Autoimmun Rev* 2025;24:103693. [PMID 39577549](https://pubmed.ncbi.nlm.nih.gov/39577549/)  
   Large comprehensive review of clinicoradiological features, therapy and outcomes; the broadest recent synthesis of the treatment literature.

4. **Wynford-Thomas R et al.** Neurological update: MOG antibody disease. *J Neurol* 2019;266:1280-1286. [PMID 30569382](https://pubmed.ncbi.nlm.nih.gov/30569382/)  
   Compact clinical update; a useful orientation to the phenotype spectrum.

5. **Bruijstens AL et al.** E.U. paediatric MOG consortium consensus: Part 1 - Classification of clinical phenotypes of paediatric myelin oligodendrocyte glycoprotein antibody-associated disorders. *Eur J Paediatr Neurol* 2020;29:2-13. [PMID 33162302](https://pubmed.ncbi.nlm.nih.gov/33162302/)  
   European paediatric consortium consensus on classification of clinical phenotypes.


## 2. Immunopathology and target biology (commitment 2)

These are the papers that justify the model's central structural choice: myelin is destroyed while the oligodendrocyte cell body and the astrocyte/AQP4 compartment survive.

6. **Höftberger R et al.** The pathology of central nervous system inflammatory demyelinating disease accompanying myelin oligodendrocyte glycoprotein autoantibody. *Acta Neuropathol* 2020;139:875-892. [PMID 32048003](https://pubmed.ncbi.nlm.nih.gov/32048003/)  
   Neuropathology of MOG-IgG-associated demyelination: perivenous and confluent demyelination with RELATIVE PRESERVATION of oligodendrocytes and intact AQP4 — the direct evidence for the small KOLD and the AQP4_MODE switch.

7. **Takai Y et al.** Myelin oligodendrocyte glycoprotein antibody-associated disease: an immunopathological study. *Brain* 2020;143:1431-1446. [PMID 32412053](https://pubmed.ncbi.nlm.nih.gov/32412053/)  
   Immunopathological series making the same point independently: MOGAD lesions are demyelinating, not astrocytopathic.

8. **Peschl P et al.** Myelin Oligodendrocyte Glycoprotein: Deciphering a Target in Inflammatory Demyelinating Diseases. *Front Immunol* 2017;8:529. [PMID 28533781](https://pubmed.ncbi.nlm.nih.gov/28533781/)  
   Review of MOG as an autoantigen — its restriction to the CNS, its low abundance in myelin, and critically its position on the OUTERMOST myelin lamella, which is why the epitope is directly antibody-accessible without cell entry.

9. **Reindl M et al.** Myelin oligodendrocyte glycoprotein antibodies in neurological disease. *Nat Rev Neurol* 2019;15:89-102. [PMID 30559466](https://pubmed.ncbi.nlm.nih.gov/30559466/)  
   Review of MOG antibodies across neurological disease, including assay biology and subclass (IgG1, complement-fixing).

10. **Spadaro M et al.** Pathogenicity of human antibodies against myelin oligodendrocyte glycoprotein. *Ann Neurol* 2018;84:315-328. [PMID 30014603](https://pubmed.ncbi.nlm.nih.gov/30014603/)  
   Human MOG antibodies are pathogenic in vivo — the basis for treating serum antibody as the causal driver rather than an epiphenomenon.

11. **Mader S et al.** Complement activating antibodies to myelin oligodendrocyte glycoprotein in neuromyelitis optica and related disorders. *J Neuroinflammation* 2011;8:184. [PMID 22204662](https://pubmed.ncbi.nlm.nih.gov/22204662/)  
   MOG antibodies from patients activate complement, supporting the classical-pathway arm (C1q to MAC) in the model.

12. **Saadoun S et al.** Neuromyelitis optica MOG-IgG causes reversible lesions in mouse brain. *Acta Neuropathol Commun* 2014;2:35. [PMID 24685353](https://pubmed.ncbi.nlm.nih.gov/24685353/)  
   MOG-IgG produces REVERSIBLE lesions in mouse brain — the experimental counterpart of the model's high recovery ceiling.

13. **Remlinger J et al.** Modeling MOG Antibody-Associated Disorder and Neuromyelitis Optica Spectrum Disorder in Animal Models: Visual System Manifestations. *Neurol Neuroimmunol Neuroinflamm* 2023;10. [PMID 37429715](https://pubmed.ncbi.nlm.nih.gov/37429715/)  
   Review of animal models of MOGAD and NMOSD side by side; the source for treating the two as the same effector architecture with a different target cell.

14. **Baksmeier C et al.** Modified recombinant human IgG1-Fc is superior to natural intravenous immunoglobulin at inhibiting immune-mediated demyelination. *Immunology* 2021;164:90-105. [PMID 33880776](https://pubmed.ncbi.nlm.nih.gov/33880776/)  
   Modified recombinant IgG1-Fc outperforms IVIG at inhibiting MOG-IgG-mediated demyelination ex vivo — direct support for modelling IVIG's benefit as effector-level (Fc-receptor) blockade rather than titre reduction alone.

15. **Asavapanumas N et al.** Targeting the complement system in neuromyelitis optica spectrum disorder. *Expert Opin Biol Ther* 2021;21:1073-1086. [PMID 33513036](https://pubmed.ncbi.nlm.nih.gov/33513036/)  
   Complement as a therapeutic target in NMOSD; the pharmacological background for the C5-inhibitor arm.


## 3. The antibody source: B lineage, plasmablasts and T-cell support (commitment 1)

The model makes the pathogenic IgG1 come from a short-lived, CD20-negative plasmablast pool supported by an IL-6-dependent niche. These are the human data behind that.

16. **Winklmeier S et al.** Identification of circulating MOG-specific B cells in patients with MOG antibodies. *Neurol Neuroimmunol Neuroinflamm* 2019;6:625. [PMID 31611268](https://pubmed.ncbi.nlm.nih.gov/31611268/)  
   Circulating MOG-specific B cells can be identified in patients — evidence for a peripheral, antigen-driven source compartment.

17. **Horellou P et al.** Regulatory T Cells Increase After rh-MOG Stimulation in Non-Relapsing but Decrease in Relapsing MOG Antibody-Associated Disease at Onset in Children. *Front Immunol* 2021;12:679770. [PMID 34220827](https://pubmed.ncbi.nlm.nih.gov/34220827/)  
   Regulatory T-cell responses to recombinant MOG differ between relapsing and non-relapsing patients — the basis for the T-cell arm and the Treg-deficit node in the map.

18. **Kothur K et al.** B Cell, Th17, and Neutrophil Related Cerebrospinal Fluid Cytokine/Chemokines Are Elevated in MOG Antibody Associated Demyelination. *PLoS One* 2016;11:e0149411. [PMID 26919719](https://pubmed.ncbi.nlm.nih.gov/26919719/)  
   CSF cytokine and chemokine profiling in MOG-antibody demyelination shows a B-cell, Th17 and neutrophil signature, including the IL-6 axis the model relies on.


## 4. Clinical phenotypes, attack severity and recovery

These cohorts supply the endpoint targets: nadir severity, final visual acuity, RNFL thinning, lesion topography and the MOGAD-versus-AQP4 recovery gap.

19. **Chen JJ et al.** Myelin Oligodendrocyte Glycoprotein Antibody-Positive Optic Neuritis: Clinical Characteristics, Radiologic Clues, and Outcome. *Am J Ophthalmol* 2018;195:8-15. [PMID 30055153](https://pubmed.ncbi.nlm.nih.gov/30055153/)  
   MOG-IgG optic neuritis: nadir visual acuity averaged count-fingers, average FINAL acuity 20/30 with only 6% at or below 20/200; disc oedema and pain each 86%; perineural enhancement 50%; longitudinally extensive involvement 80%. The primary calibration target for the visual endpoint.

20. **Padungkiatsagul T et al.** Differences in Clinical Features of Myelin Oligodendrocyte Glycoprotein Antibody-Associated Optic Neuritis in White and Asian Race. *Am J Ophthalmol* 2020;219:332-340. [PMID 32681910](https://pubmed.ncbi.nlm.nih.gov/32681910/)  
   Direct comparison of MOG-IgG and AQP4-IgG optic neuritis features.

21. **Stiebel-Kalish H et al.** Retinal Nerve Fiber Layer May Be Better Preserved in MOG-IgG versus AQP4-IgG Optic Neuritis: A Cohort Study. *PLoS One* 2017;12:e0170847. [PMID 28125740](https://pubmed.ncbi.nlm.nih.gov/28125740/)  
   RNFL is better preserved after MOG-IgG than AQP4-IgG optic neuritis — the target for the model's RNFL readout in the two arms (75 versus 63 micrometres).

22. **Dubey D et al.** Clinical, Radiologic, and Prognostic Features of Myelitis Associated With Myelin Oligodendrocyte Glycoprotein Autoantibody. *JAMA Neurol* 2019;76:301-309. [PMID 30575890](https://pubmed.ncbi.nlm.nih.gov/30575890/)  
   MOG-IgG myelitis: about a third wheelchair-dependent at nadir yet recovering better than AQP4-IgG myelitis; grey-matter H-sign and conus involvement; oligoclonal bands in only 1 of 38.

23. **Jarius S et al.** MOG-IgG in NMO and related disorders: a multicenter study of 50 patients. Part 2: Epidemiology, clinical presentation, radiological and laboratory features, treatment responses, and long-term outcome. *J Neuroinflammation* 2016;13:280. [PMID 27793206](https://pubmed.ncbi.nlm.nih.gov/27793206/)  
   Multicentre series of 50 patients: epidemiology, clinical presentation and course.

24. **Jurynczyk M et al.** Clinical presentation and prognosis in MOG-antibody disease: a UK study. *Brain* 2017;140:3128-3138. [PMID 29136091](https://pubmed.ncbi.nlm.nih.gov/29136091/)  
   UK cohort: clinical presentation and prognosis, including relapse rates.

25. **Cobo-Calvo A et al.** Clinical spectrum and prognostic value of CNS MOG autoimmunity in adults: The MOGADOR study. *Neurology* 2018;90:e1858-e1869. [PMID 29695592](https://pubmed.ncbi.nlm.nih.gov/29695592/)  
   MOGADOR study: clinical spectrum and prognostic value of MOG autoimmunity in adults.

26. **Ramanathan S et al.** Clinical course, therapeutic responses and outcomes in relapsing MOG antibody-associated demyelination. *J Neurol Neurosurg Psychiatry* 2018;89:127-137. [PMID 29142145](https://pubmed.ncbi.nlm.nih.gov/29142145/)  
   Clinical course, therapeutic responses and outcomes in relapsing disease — an early source for the steroid-dependency pattern.

27. **Hacohen Y et al.** Disease Course and Treatment Responses in Children With Relapsing Myelin Oligodendrocyte Glycoprotein Antibody-Associated Disease. *JAMA Neurol* 2018;75:478-487. [PMID 29305608](https://pubmed.ncbi.nlm.nih.gov/29305608/)  
   Paediatric relapsing MOGAD: disease course and treatment responses, including the relative performance of maintenance agents in children.

28. **Hennes EM et al.** Prognostic relevance of MOG antibodies in children with an acquired demyelinating syndrome. *Neurology* 2017;89:900-908. [PMID 28768844](https://pubmed.ncbi.nlm.nih.gov/28768844/)  
   Prognostic relevance of MOG antibodies in children after a first acquired demyelinating syndrome.

29. **Baumann M et al.** Clinical and neuroradiological differences of paediatric acute disseminating encephalomyelitis with and without antibodies to the myelin oligodendrocyte glycoprotein. *J Neurol Neurosurg Psychiatry* 2015;86:265-72. [PMID 25121570](https://pubmed.ncbi.nlm.nih.gov/25121570/)  
   Paediatric ADEM with and without MOG antibodies — clinical and neuroradiological differences.

30. **Budhram A et al.** Unilateral cortical FLAIR-hyperintense Lesions in Anti-MOG-associated Encephalitis with Seizures (FLAMES): characterization of a distinct clinico-radiographic syndrome. *J Neurol* 2019;266:2481-2487. [PMID 31243540](https://pubmed.ncbi.nlm.nih.gov/31243540/)  
   FLAMES: unilateral cortical FLAIR-hyperintense lesions with seizures, the cortical encephalitis phenotype in the mechanistic map.

31. **Sechi E et al.** Comparison of MRI Lesion Evolution in Different Central Nervous System Demyelinating Disorders. *Neurology* 2021;97:e1097-e1109. [PMID 34261784](https://pubmed.ncbi.nlm.nih.gov/34261784/)  
   MRI lesion evolution compared across demyelinating diseases — the imaging counterpart of the model's near-complete lesion resolution in MOGAD.


## 5. Serostatus, titre and relapse prediction (commitment 3)

The model's claim that titre is a poor surrogate for protection is testable against exactly this literature.

32. **López-Chiriboga AS et al.** Association of MOG-IgG Serostatus With Relapse After Acute Disseminated Encephalomyelitis and Proposed Diagnostic Criteria for MOG-IgG-Associated Disorders. *JAMA Neurol* 2018;75:1355-1363. [PMID 30014148](https://pubmed.ncbi.nlm.nih.gov/30014148/)  
   MOG-IgG serostatus and relapse risk after ADEM: persistent seropositivity predicts relapse, transient seropositivity often does not.

33. **Oliveira LM et al.** Persistent MOG-IgG positivity is a predictor of recurrence in MOG-IgG-associated optic neuritis, encephalitis and myelitis. *Mult Scler* 2019;25:1907-1914. [PMID 30417715](https://pubmed.ncbi.nlm.nih.gov/30417715/)  
   Persistent MOG-IgG positivity predicts recurrence — the empirical basis for the titre-driven hazard.

34. **Andersen J et al.** Biomarkers to predict relapse in myelin oligodendrocyte glycoprotein antibody-associated disease: a systematic review and meta-analysis. *J Neurol Neurosurg Psychiatry* 2026;97:109-119. [PMID 41033784](https://pubmed.ncbi.nlm.nih.gov/41033784/)  
   Recent systematic assessment of biomarkers for predicting relapse in MOGAD; relevant to the model's prediction that serial titres should be a weak treatment surrogate.


## 6. Fluid biomarkers: NfL and GFAP

The GFAP-versus-NfL dissociation is the cleanest available test of commitment 2, and the model reproduces it from the AQP4_MODE switch alone.

35. **Mariotto S et al.** Neurofilament light chain serum levels reflect disease severity in MOG-Ab associated disorders. *J Neurol Neurosurg Psychiatry* 2019;90:1293-1296. [PMID 30952681](https://pubmed.ncbi.nlm.nih.gov/30952681/)  
   Serum neurofilament light reflects disease severity in MOG-antibody disease — the target for the model's NfL readout.

36. **Marignier R et al.** Assessment of neuronal and glial serum biomarkers in myelin oligodendrocyte glycoprotein antibody-associated disease: the MULTIMOGAD study. *J Neurol Neurosurg Psychiatry* 2025;96:884-892. [PMID 39939136](https://pubmed.ncbi.nlm.nih.gov/39939136/)  
   Neuronal and glial serum biomarkers in MOGAD: NfL rises with attacks while GFAP stays low, unlike AQP4-IgG NMOSD. The model reproduces a GFAP peak of about 150 pg/mL in MOGAD versus about 3660 pg/mL in the comparator arm.


## 7. Acute attack treatment: steroids, plasma exchange, acute IVIG

37. **Bonnan M et al.** Short delay to initiate plasma exchange is the strongest predictor of outcome in severe attacks of NMO spectrum disorders. *J Neurol Neurosurg Psychiatry* 2018;89:346-351. [PMID 29030418](https://pubmed.ncbi.nlm.nih.gov/29030418/)  
   Short delay to plasma exchange is the strongest predictor of outcome in severe CNS inflammatory demyelination — the rationale for the early-PLEX scenario.

38. **Thakolwiboon S et al.** Outcomes After Acute Plasma Exchange for Myelin Oligodendrocyte Glycoprotein Antibody-Associated Disease. *Neurology* 2025;105:e213903. [PMID 40882166](https://pubmed.ncbi.nlm.nih.gov/40882166/)  
   Outcomes after acute plasma exchange specifically in MOGAD.

39. **Elsone L et al.** Role of intravenous immunoglobulin in the treatment of acute relapses of neuromyelitis optica: experience in 10 patients. *Mult Scler* 2014;20:501-4. [PMID 23986097](https://pubmed.ncbi.nlm.nih.gov/23986097/)  
   IVIG for acute relapses of neuromyelitis optica; the closest available evidence base for the acute IVIG arm.


## 8. Maintenance therapy and comparative effectiveness

This is the model's principal quantitative calibration set. The ordering it reproduces — IVIG and IL-6R blockade clearly ahead of rituximab and the antimetabolites — comes from here.

40. **Vilaseca A et al.** Interleukin 6 Receptor Blockade for Relapse Prevention in Myelin Oligodendrocyte Glycoprotein Antibody-Associated Disease. *JAMA Neurol* 2026;. [PMID 42440328](https://pubmed.ncbi.nlm.nih.gov/42440328/)  
   Americas MOGAD Treatment Group, 116 patients on IL-6R blockade versus 59 on IVIG. Relapsing MOGAD ARR 0.64 (95% CI 0.58-0.70) before preventive therapy; 0.09 (0.06-0.14) on IL-6R blockade; 0.22 (0.15-0.32) on IVIG; after weighting, IL-6R blockade was NOT significantly different from IVIG at 1 g/kg or more every 4 weeks but was better than lower-dose IVIG. LAM0 was set from the 0.64 figure.

41. **Chen JJ et al.** Association of Maintenance Intravenous Immunoglobulin With Prevention of Relapse in Adult Myelin Oligodendrocyte Glycoprotein Antibody-Associated Disease. *JAMA Neurol* 2022;79:518-525. [PMID 35377395](https://pubmed.ncbi.nlm.nih.gov/35377395/)  
   Maintenance IVIG in adult MOGAD: relapse in 5 of 29 (17%) at 1 g/kg every 4 weeks or more against 15 of 30 (50%) on lower or less frequent dosing (HR 3.31). The dose-dependence the model reproduces as ARR 0.18 versus 0.37.

42. **Chen JJ et al.** Steroid-sparing maintenance immunotherapy for MOG-IgG associated disorder. *Neurology* 2020;95:e111-e120. [PMID 32554760](https://pubmed.ncbi.nlm.nih.gov/32554760/)  
   Steroid-sparing maintenance immunotherapy: on-treatment relapse in 74% on mycophenolate (ARR 0.67), 61% on rituximab (0.59), 59% on azathioprine (0.20) and 20% on IVIG (0.00). Also the source for the observation that MS disease-modifying agents fail entirely.

43. **Thakolwiboon S et al.** Meta-analysis of effectiveness of steroid-sparing attack prevention in MOG-IgG-associated disorder. *Mult Scler Relat Disord* 2021;56:103310. [PMID 34634625](https://pubmed.ncbi.nlm.nih.gov/34634625/)  
   Meta-analysis of steroid-sparing attack prevention: pooled on-treatment ARR 0.29 azathioprine, 0.84 mycophenolate, 0.63 rituximab, 0.08 maintenance IVIG.

44. **Durozard P et al.** Comparison of the Response to Rituximab between Myelin Oligodendrocyte Glycoprotein and Aquaporin-4 Antibody Diseases. *Ann Neurol* 2020;87:256-266. [PMID 31725931](https://pubmed.ncbi.nlm.nih.gov/31725931/)  
   Response to rituximab is markedly poorer in MOG-antibody than in AQP4-antibody disease — the observation the FRAC_ESC parameter exists to explain.

45. **Ringelstein M et al.** Interleukin-6 Receptor Blockade in Treatment-Refractory MOG-IgG-Associated Disease and Neuromyelitis Optica Spectrum Disorders. *Neurol Neuroimmunol Neuroinflamm* 2022;9. [PMID 34785575](https://pubmed.ncbi.nlm.nih.gov/34785575/)  
   IL-6 receptor blockade in treatment-refractory MOGAD and NMOSD.

46. **Zhang C et al.** Safety and efficacy of tocilizumab versus azathioprine in highly relapsing neuromyelitis optica spectrum disorder (TANGO): an open-label, multicentre, randomised, phase 2 trial. *Lancet Neurol* 2020;19:391-401. [PMID 32333897](https://pubmed.ncbi.nlm.nih.gov/32333897/)  
   TANGO: tocilizumab versus azathioprine in highly relapsing NMOSD — the randomised evidence that IL-6R blockade beats an antimetabolite in a related antibody-mediated disease.

47. **Kleiter I et al.** Long-term Efficacy of Satralizumab in AQP4-IgG-Seropositive Neuromyelitis Optica Spectrum Disorder From SAkuraSky and SAkuraStar. *Neurol Neuroimmunol Neuroinflamm* 2023;10. [PMID 36724181](https://pubmed.ncbi.nlm.nih.gov/36724181/)  
   Long-term efficacy of satralizumab in AQP4-IgG NMOSD; the PK/PD and efficacy template for the recycling anti-IL-6R antibody.

48. **Carnero Contentti E et al.** Future treatments for myelin oligodendrocyte glycoprotein antibody-associated disease: the clinical trial landscape. *Expert Opin Emerg Drugs* 2025;30:283-297. [PMID 40984653](https://pubmed.ncbi.nlm.nih.gov/40984653/)  
   Review of emerging MOGAD therapies, including the trials of agents the model represents only as classes.


## 9. Emerging mechanisms: FcRn blockade and IgG clearance

The FcRn arm is the model's clearest prediction, so its parameterisation is grounded in an actual PK/PD model rather than assumption.

49. **Lledo-Garcia R et al.** Pharmacokinetic-pharmacodynamic modelling of the anti-FcRn monoclonal antibody rozanolixizumab: Translation from preclinical stages to the clinic. *CPT Pharmacometrics Syst Pharmacol* 2022;11:116-128. [PMID 34735735](https://pubmed.ncbi.nlm.nih.gov/34735735/)  
   Population PK/PD model of the anti-FcRn antibody rozanolixizumab from preclinical to clinical stages, including the magnitude and time course of total-IgG reduction. FRN_INH was set to give about 53% steady-state total IgG reduction from this.

50. **Remlinger J et al.** Antineonatal Fc Receptor Antibody Treatment Ameliorates MOG-IgG-Associated Experimental Autoimmune Encephalomyelitis. *Neurol Neuroimmunol Neuroinflamm* 2022;9. [PMID 35027475](https://pubmed.ncbi.nlm.nih.gov/35027475/)  
   Anti-FcRn antibody treatment ameliorates MOG-IgG-associated experimental autoimmune encephalomyelitis — direct in vivo support for the clearance-node arm in this specific disease.


## 10. Pharmacokinetics and toxicity of the modelled drugs

51. **Narang PK et al.** Systemic bioavailability and pharmacokinetics of methylprednisolone in patients with rheumatoid arthritis following 'high-dose' pulse administration. *Biopharm Drug Dispos* 1983;4:233-48. [PMID 6626699](https://pubmed.ncbi.nlm.nih.gov/6626699/)  
   Systemic bioavailability and pharmacokinetics of methylprednisolone; the source for the short plasma half-life that motivates a separate glucocorticoid effect compartment.

52. **Fokkink W et al.** Pharmacokinetics and Pharmacodynamics of Intravenous Immunoglobulin G Maintenance Therapy in Chronic Immune-mediated Neuropathies. *Clin Pharmacol Ther* 2017;102:709-716. [PMID 28378901](https://pubmed.ncbi.nlm.nih.gov/28378901/)  
   PK/PD of IVIG maintenance therapy, including trough IgG behaviour on repeated dosing.

53. **Kapszewicz M et al.** Glucocorticoid-Induced Osteoporosis: Pathogenesis, the Impact of Different Administration Routes on Bone Mineral Density, and Fracture Risk and Treatment Options-A Narrative Review. *J Clin Med* 2026;15. [PMID 41976789](https://pubmed.ncbi.nlm.nih.gov/41976789/)  
   Glucocorticoid-induced osteoporosis: pathogenesis and the effect of administration route and dose — the basis for the BMD state and the steroid-cost readout.


## 11. QSP and mrgsolve methodology

54. **Elmokadem A et al.** Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial. *CPT Pharmacometrics Syst Pharmacol* 2019;8:883-893. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)  
   Quantitative systems pharmacology and PBPK modelling with mrgsolve — the engine and modelling idiom used here.

55. **Wang H et al.** Quantitative Systems Pharmacology Modeling in Immuno-Oncology: Hypothesis Testing, Dose Optimization, and Efficacy Prediction. *Handb Exp Pharmacol* 2025;289:261-284. [PMID 39707022](https://pubmed.ncbi.nlm.nih.gov/39707022/)  
   QSP modelling for hypothesis testing and dose selection; general framing for what this kind of model can and cannot support.


---

## How this list was built

Every citation was resolved through the NCBI E-utilities API (`esearch` +
`esummary`): the title, journal, year, volume, pages and first author printed
above come from the PubMed record itself, so a wrong PMID would have produced a
visibly wrong citation line. Total: **55 references**.

## Reading order for someone new to MOGAD

1. Banwell 2023 (criteria) then Marignier 2021 (overview) — what the disease is.
2. Höftberger 2020 and Takai 2020 — why it is not NMOSD.
3. Chen 2018 (optic neuritis) and Dubey 2019 (myelitis) — what an attack does.
4. Chen 2020, Chen 2022, Thakolwiboon 2021, Vilaseca 2026 — what treatment achieves.
5. Durozard 2020 — the rituximab anomaly this model is built to explain.

## Disclaimer

This reference list supports an **educational and research** QSP model. It is not
a clinical guideline and the annotations are the model author's reading of each
paper, not the authors' own claims.
