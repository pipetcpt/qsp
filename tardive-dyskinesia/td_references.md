# 지연성 운동이상증 (Tardive Dyskinesia) QSP 모델 — 참고문헌 (References)

지연성 운동이상증의 기계론적 지도(`td_qsp_model.dot`), mrgsolve ODE 모델
(`td_mrgsolve_model.R`), Shiny 대시보드(`td_shiny_app.R`) 구축에 사용한 문헌
목록입니다. 총 72편.

> **링크 형식에 관한 註**: 아래 링크는 모두 **PubMed 검색 질의(query) URL**입니다.
> PMID를 잘못 기재하면 전혀 다른 논문을 인용하게 되므로, 저자·연도·제목
> 키워드로 검색되는 링크를 사용했습니다. 모델 파라미터가 특정 문헌에
> 정박(anchor)된 부분은 `td_mrgsolve_model.R` 상단의 CALIBRATION ANCHORS
> 절에 대응 관계가 정리되어 있습니다.
>
> **중요**: README와 모델 파일에 인용된 수치 중 "모델:" 로 표시된 값은 모두
> 이 모델의 계산 결과이며 문헌 값이 아닙니다. 문헌은 파라미터 보정의
> 근거로만 사용되었습니다.

---

## 1. 진단 기준 · 평가 척도 (Diagnosis & rating scales)

1. Guy W. **ECDEU Assessment Manual for Psychopharmacology — Abnormal
   Involuntary Movement Scale (AIMS).** 1976.
   <https://pubmed.ncbi.nlm.nih.gov/?term=abnormal+involuntary+movement+scale+AIMS+Guy+ECDEU+assessment+manual>
2. Schooler NR, Kane JM. **Research diagnoses for tardive dyskinesia.**
   *Arch Gen Psychiatry.* 1982.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Schooler+Kane+research+diagnoses+tardive+dyskinesia+1982>
3. Stacy M, et al. **Assessment of interrater and intrarater reliability of
   the Abnormal Involuntary Movement Scale in tardive dyskinesia.**
   *Mov Disord.* 2007.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Stacy+interrater+intrarater+reliability+abnormal+involuntary+movement+scale+tardive+dyskinesia>
4. Kane JM, et al. **Tardive dyskinesia: prevalence, incidence, and risk
   factors.** *J Clin Psychopharmacol.* 1988.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Kane+tardive+dyskinesia+prevalence+incidence+risk+factors+prospective>
5. Carbon M, et al. **Tardive dyskinesia prevalence in the period of
   second-generation antipsychotic use: a meta-analysis.**
   *J Clin Psychiatry.* 2017.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Carbon+tardive+dyskinesia+prevalence+second-generation+antipsychotic+meta-analysis+2017>
6. Correll CU, Schenk EM. **Tardive dyskinesia and new antipsychotics.**
   *Curr Opin Psychiatry.* 2008.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Correll+Schenk+tardive+dyskinesia+new+antipsychotics+annualized+incidence>
7. Woods SW, et al. **Incidence of tardive dyskinesia with atypical versus
   conventional antipsychotic medications: a prospective cohort study.**
   *J Clin Psychiatry.* 2010.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Woods+incidence+tardive+dyskinesia+atypical+conventional+antipsychotic+prospective+cohort>

## 2. 병태생리 (1) — D2 수용체 초민감성 가설

8. Klawans HL, Rubovits R. **An experimental model of tardive dyskinesia.**
   *J Neural Transm.* 1972.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Klawans+Rubovits+experimental+model+tardive+dyskinesia+1972>
9. Burt DR, Creese I, Snyder SH. **Antischizophrenic drugs: chronic
   treatment elevates dopamine receptor binding in brain.** *Science.* 1977.
   <https://pubmed.ncbi.nlm.nih.gov/?term=Burt+Creese+Snyder+antischizophrenic+drugs+chronic+treatment+elevates+dopamine+receptor+binding>
10. Seeman P. **All roads to schizophrenia lead to dopamine
    supersensitivity and elevated dopamine D2High receptors.**
    *CNS Neurosci Ther.* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Seeman+all+roads+schizophrenia+dopamine+supersensitivity+D2High+receptors>
11. Samaha AN, et al. **"Breakthrough" dopamine supersensitivity during
    ongoing antipsychotic treatment leads to treatment failure over time.**
    *J Neurosci.* 2007.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Samaha+breakthrough+dopamine+supersensitivity+ongoing+antipsychotic+treatment+failure+over+time>
12. Silvestri S, et al. **Increased dopamine D2 receptor binding after
    long-term treatment with antipsychotics in humans: a clinical PET
    study.** *Psychopharmacology.* 2000.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Silvestri+increased+dopamine+D2+receptor+binding+long-term+antipsychotics+humans+PET>
13. Chouinard G, Jones BD. **Neuroleptic-induced supersensitivity
    psychosis: clinical and pharmacologic characteristics.**
    *Am J Psychiatry.* 1980.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Chouinard+Jones+neuroleptic-induced+supersensitivity+psychosis+clinical+pharmacologic+characteristics>
14. Iyo M, et al. **Optimal extent of dopamine D2 receptor occupancy by
    antipsychotics for treatment of dopamine supersensitivity psychosis.**
    *J Clin Psychopharmacol.* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Iyo+optimal+extent+dopamine+D2+receptor+occupancy+antipsychotics+dopamine+supersensitivity+psychosis>
15. Teo JT, Edwards MJ, Bhatia K. **Tardive dyskinesia is caused by
    maladaptive synaptic plasticity: a hypothesis.** *Mov Disord.* 2012.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Teo+Edwards+Bhatia+tardive+dyskinesia+maladaptive+synaptic+plasticity+hypothesis>

## 3. 병태생리 (2) — 산화 스트레스 · 구조적 손상 · 개재뉴런 소실

16. Cadet JL, Lohr JB. **Free radicals and the pathobiology of brain
    dopamine systems.** *Ann N Y Acad Sci.* 1987.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Cadet+Lohr+free+radicals+pathobiology+brain+dopamine+systems+tardive+dyskinesia>
17. Andreassen OA, Jørgensen HA. **Neurotoxicity associated with
    neuroleptic-induced oral dyskinesias in rats: implications for tardive
    dyskinesia?** *Prog Neurobiol.* 2000.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Andreassen+Jorgensen+neurotoxicity+neuroleptic-induced+oral+dyskinesias+rats+tardive+dyskinesia>
18. Miller R, Chouinard G. **Loss of striatal cholinergic neurons as a
    basis for tardive and L-dopa-induced dyskinesias, neuroleptic-induced
    supersensitivity psychosis and refractory schizophrenia.**
    *Biol Psychiatry.* 1993.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Miller+Chouinard+loss+striatal+cholinergic+neurons+tardive+dyskinesia+supersensitivity+psychosis>
19. Sachdev PS. **The current status of tardive dyskinesia.**
    *Aust N Z J Psychiatry.* 2000.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Sachdev+current+status+tardive+dyskinesia+pathophysiology+review>
20. Lohr JB, et al. **Oxidative mechanisms and tardive dyskinesia.**
    *CNS Drugs.* 2003.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Lohr+oxidative+mechanisms+tardive+dyskinesia+CNS+Drugs+2003>
21. Zhang XY, et al. **The novel oxidative stress marker and superoxide
    dismutase activity in patients with tardive dyskinesia.**
    *Schizophr Res.* 2007.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Zhang+superoxide+dismutase+activity+oxidative+stress+patients+tardive+dyskinesia>
22. Konradi C, Heckers S. **Antipsychotic drugs and neuroplasticity:
    insights into the treatment and neurobiology of schizophrenia.**
    *Biol Psychiatry.* 2001.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Konradi+Heckers+antipsychotic+drugs+neuroplasticity+treatment+neurobiology+schizophrenia>

## 4. 선조체 미세회로 · 기저핵 루프

23. Albin RL, Young AB, Penney JB. **The functional anatomy of basal
    ganglia disorders.** *Trends Neurosci.* 1989.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Albin+Young+Penney+functional+anatomy+basal+ganglia+disorders+1989>
24. DeLong MR. **Primate models of movement disorders of basal ganglia
    origin.** *Trends Neurosci.* 1990.
    <https://pubmed.ncbi.nlm.nih.gov/?term=DeLong+primate+models+movement+disorders+basal+ganglia+origin+1990>
25. Gerfen CR, Surmeier DJ. **Modulation of striatal projection systems by
    dopamine.** *Annu Rev Neurosci.* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Gerfen+Surmeier+modulation+striatal+projection+systems+by+dopamine>
26. Bordia T, Perez XA. **Cholinergic control of striatal neurons to
    modulate L-dopa-induced dyskinesias.** *Eur J Neurosci.* 2019.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Bordia+Perez+cholinergic+control+striatal+neurons+modulate+dyskinesias>
27. Aquino CC, Lang AE. **Tardive dyskinesia syndromes: current concepts.**
    *Parkinsonism Relat Disord.* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Aquino+Lang+tardive+dyskinesia+syndromes+current+concepts>

## 5. D2 점유율 · 항정신병제 PK/PD

28. Farde L, et al. **Positron emission tomographic analysis of central D1
    and D2 dopamine receptor occupancy in patients treated with classical
    neuroleptics and clozapine.** *Arch Gen Psychiatry.* 1992.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Farde+positron+emission+tomographic+analysis+central+D1+D2+dopamine+receptor+occupancy+clozapine>
29. Kapur S, et al. **Relationship between dopamine D2 occupancy, clinical
    response, and side effects: a double-blind PET study of
    first-episode schizophrenia.** *Am J Psychiatry.* 2000.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Kapur+relationship+between+dopamine+D2+occupancy+clinical+response+side+effects+first-episode+schizophrenia>
30. Nyberg S, et al. **Suggested minimal effective dose of risperidone
    based on PET-measured D2 and 5-HT2A receptor occupancy.**
    *Am J Psychiatry.* 1999.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Nyberg+suggested+minimal+effective+dose+risperidone+PET+D2+5-HT2A+receptor+occupancy>
31. Kapur S, Seeman P. **Does fast dissociation from the dopamine D2
    receptor explain the action of atypical antipsychotics?**
    *Am J Psychiatry.* 2001.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Kapur+Seeman+fast+dissociation+dopamine+D2+receptor+atypical+antipsychotics>
32. Sherwood M, et al. **Model-based meta-analysis of relationships between
    dopamine D2 receptor occupancy and antipsychotic efficacy.**
    *Br J Clin Pharmacol.* 2012.
    <https://pubmed.ncbi.nlm.nih.gov/?term=model-based+meta-analysis+dopamine+D2+receptor+occupancy+antipsychotic+efficacy>

## 6. VMAT2 억제제 — 발베나진 (valbenazine)

33. Hauser RA, et al. **KINECT 3: a phase 3 randomized, double-blind,
    placebo-controlled trial of valbenazine for tardive dyskinesia.**
    *Am J Psychiatry.* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Hauser+KINECT+3+phase+3+randomized+valbenazine+tardive+dyskinesia+2017>
34. Factor SA, et al. **The effects of valbenazine in participants with
    tardive dyskinesia: results of the 1-year KINECT 3 extension study.**
    *J Clin Psychiatry.* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Factor+valbenazine+KINECT+3+extension+study+1-year+tardive+dyskinesia>
35. Marder SR, et al. **Long-term safety and tolerability of valbenazine
    (KINECT 4) in adults with tardive dyskinesia.** *J Clin Psychiatry.*
    2019.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Marder+long-term+safety+tolerability+valbenazine+KINECT+4+tardive+dyskinesia>
36. Grigoriadis DE, et al. **Pharmacologic characterization of
    valbenazine (NBI-98854) and its metabolites.**
    *J Pharmacol Exp Ther.* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Grigoriadis+pharmacologic+characterization+valbenazine+NBI-98854+metabolites+VMAT2>
37. O'Brien CF, et al. **NBI-98854, a selective monoamine transport
    inhibitor for the treatment of tardive dyskinesia: a randomized,
    double-blind, placebo-controlled study.** *Mov Disord.* 2015.
    <https://pubmed.ncbi.nlm.nih.gov/?term=OBrien+NBI-98854+selective+monoamine+transport+inhibitor+tardive+dyskinesia+randomized>

## 7. VMAT2 억제제 — 듀테트라베나진 · 테트라베나진

38. Anderson KE, et al. **Deutetrabenazine for treatment of involuntary
    movements in patients with tardive dyskinesia (AIM-TD): a double-blind,
    randomised, placebo-controlled, phase 3 trial.** *Lancet Psychiatry.*
    2017.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Anderson+deutetrabenazine+AIM-TD+double-blind+randomised+placebo-controlled+phase+3+tardive+dyskinesia>
39. Fernandez HH, et al. **Randomized controlled trial of deutetrabenazine
    for tardive dyskinesia: the ARM-TD study.** *Neurology.* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Fernandez+randomized+controlled+trial+deutetrabenazine+tardive+dyskinesia+ARM-TD>
40. Schneider F, et al. **Deutetrabenazine for tardive dyskinesia: results
    from an open-label, long-term extension study.**
    *J Clin Psychiatry / Mov Disord Clin Pract.* 2020.
    <https://pubmed.ncbi.nlm.nih.gov/?term=deutetrabenazine+tardive+dyskinesia+open-label+long-term+extension+study+3-year>
41. Ondo WG, et al. **Tetrabenazine treatment for tardive dyskinesia:
    assessment by randomized videotape protocol.** *Am J Psychiatry.* 1999.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ondo+tetrabenazine+treatment+tardive+dyskinesia+randomized+videotape+protocol>
42. Stamler D, et al. **Pharmacokinetics of deutetrabenazine and
    tetrabenazine: dose proportionality and food effect.**
    *Clin Pharmacol Drug Dev.* 2018.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Stamler+pharmacokinetics+deutetrabenazine+tetrabenazine+dose+proportionality+food+effect>

## 8. CYP2D6 · 약물유전학 · 노출-반응

43. Thai HT, et al. **Population pharmacokinetics of valbenazine and its
    active metabolite in subjects with tardive dyskinesia.**
    *Clin Pharmacokinet / J Clin Pharmacol.* 2018.
    <https://pubmed.ncbi.nlm.nih.gov/?term=population+pharmacokinetics+valbenazine+active+metabolite+NBI-98782+tardive+dyskinesia+CYP2D6>
44. Zhu H, et al. **Effect of CYP2D6 phenotype on exposure to
    deutetrabenazine metabolites and dose recommendations.**
    *Clin Transl Sci / Br J Clin Pharmacol.* 2019.
    <https://pubmed.ncbi.nlm.nih.gov/?term=CYP2D6+poor+metabolizer+deutetrabenazine+alpha+beta+dihydrotetrabenazine+exposure+dose+recommendation>
45. Patsopoulos NA, et al. **CYP2D6 polymorphisms and the risk of tardive
    dyskinesia in schizophrenia: a meta-analysis.**
    *Pharmacogenet Genomics.* 2005.
    <https://pubmed.ncbi.nlm.nih.gov/?term=CYP2D6+polymorphisms+risk+tardive+dyskinesia+schizophrenia+meta-analysis>
46. Zai CC, et al. **Genetics of tardive dyskinesia: a systematic review.**
    *J Psychiatr Res / Mol Neuropsychiatry.* 2018.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Zai+genetics+tardive+dyskinesia+systematic+review+DRD2+SOD2+HTR2A>
47. Bakker PR, et al. **Antipsychotic-induced tardive dyskinesia and the
    Ser9Gly polymorphism in the DRD3 gene: a meta-analysis.**
    *Schizophr Res.* 2006.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Bakker+antipsychotic-induced+tardive+dyskinesia+DRD3+Ser9Gly+polymorphism+meta-analysis>

## 9. 원인 약물 전략 — 감량 · 중단 · 전환

48. Glazer WM, et al. **Tardive dyskinesia and the course of neuroleptic
    withdrawal: outcome after discontinuation.** *Br J Psychiatry.* 1990.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Glazer+tardive+dyskinesia+neuroleptic+withdrawal+outcome+after+discontinuation>
49. Zutshi D, Cloud LJ, Factor SA. **Tardive syndromes are rarely
    reversible after discontinuing dopamine receptor blocking agents.**
    *Tremor Other Hyperkinet Mov.* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Zutshi+Cloud+Factor+tardive+syndromes+rarely+reversible+after+discontinuing+dopamine+receptor+blocking+agents>
50. Bergman H, et al. **Antipsychotic reduction and/or cessation and
    antipsychotics as specific treatments for tardive dyskinesia.**
    *Cochrane Database Syst Rev.* 2018.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Bergman+antipsychotic+reduction+cessation+specific+treatments+tardive+dyskinesia+Cochrane>
51. Hazari N, et al. **Clozapine and tardive movement disorders: a
    review.** *Asian J Psychiatr.* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Hazari+clozapine+tardive+movement+disorders+review>
52. Lieberman JA, et al. **The effects of clozapine on tardive
    dyskinesia.** *Br J Psychiatry.* 1991.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Lieberman+effects+of+clozapine+on+tardive+dyskinesia+1991>
53. Emsley R, et al. **Antipsychotic discontinuation and relapse in
    schizophrenia: a systematic review of long-term outcomes.**
    *BMC Psychiatry / Schizophr Res.* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Emsley+antipsychotic+discontinuation+relapse+schizophrenia+systematic+review+long-term>

## 10. 보조 치료 · 비약물 치료

54. Angus S, et al. **A controlled trial of amantadine in drug-induced
    extrapyramidal disorders / tardive dyskinesia.**
    *J Clin Psychopharmacol.* 1997.
    <https://pubmed.ncbi.nlm.nih.gov/?term=amantadine+randomized+controlled+trial+tardive+dyskinesia+clinical+psychopharmacology>
55. Zhang WF, et al. **Extract of Ginkgo biloba treatment for tardive
    dyskinesia in schizophrenia: a randomized, double-blind,
    placebo-controlled trial.** *J Clin Psychiatry.* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Zhang+extract+Ginkgo+biloba+treatment+tardive+dyskinesia+schizophrenia+randomized+double-blind>
56. Soares-Weiser K, et al. **Vitamin E for antipsychotic-induced tardive
    dyskinesia.** *Cochrane Database Syst Rev.* 2018.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Soares-Weiser+vitamin+E+antipsychotic-induced+tardive+dyskinesia+Cochrane>
57. Lerner V, et al. **Vitamin B6 treatment for tardive dyskinesia: a
    randomized, double-blind, placebo-controlled, crossover study.**
    *J Clin Psychiatry.* 2007.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Lerner+vitamin+B6+treatment+tardive+dyskinesia+randomized+double-blind+crossover>
58. Rodrigues FB, et al. **Botulinum toxin for tardive syndromes and
    oromandibular dystonia.** *Cochrane / Mov Disord Clin Pract.* 2020.
    <https://pubmed.ncbi.nlm.nih.gov/?term=botulinum+toxin+tardive+dyskinesia+oromandibular+dystonia+review>
59. Macerollo A, Deuschl G. **Deep brain stimulation for tardive
    syndromes: systematic review and meta-analysis.** *J Neurol Sci.* 2018.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Macerollo+Deuschl+deep+brain+stimulation+tardive+syndromes+systematic+review+meta-analysis>
60. Bergman H, Soares-Weiser K. **Anticholinergic medication for
    antipsychotic-induced tardive dyskinesia.**
    *Cochrane Database Syst Rev.* 2018.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Bergman+Soares-Weiser+anticholinergic+medication+antipsychotic-induced+tardive+dyskinesia+Cochrane>

## 11. 위험인자 · 특수 인구

61. Woerner MG, et al. **Prospective study of tardive dyskinesia in the
    elderly: rates and risk factors.** *Am J Psychiatry.* 1998.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Woerner+prospective+study+tardive+dyskinesia+elderly+rates+risk+factors>
62. Rana AQ, et al. **Metoclopramide-induced tardive dyskinesia and other
    tardive syndromes from non-psychiatric dopamine blockers.**
    *Int J Neurosci / Clin Neuropharmacol.* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/?term=metoclopramide-induced+tardive+dyskinesia+non-psychiatric+dopamine+receptor+blocking+agents>
63. Ballesteros J, et al. **Tardive dyskinesia associated with higher
    mortality in psychiatric patients: results of a meta-analysis.**
    *J Clin Psychopharmacol.* 2000.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ballesteros+tardive+dyskinesia+higher+mortality+psychiatric+patients+meta-analysis>
64. McEvoy J, et al. **Effect of tardive dyskinesia on quality of life and
    functioning: a real-world analysis.**
    *Neuropsychiatr Dis Treat / J Clin Psychiatry.* 2019.
    <https://pubmed.ncbi.nlm.nih.gov/?term=McEvoy+tardive+dyskinesia+quality+of+life+functioning+real-world+burden>

## 12. 가이드라인 · 종합 리뷰

65. Bhidayasiri R, et al. **Evidence-based guideline: treatment of tardive
    syndromes — report of the Guideline Development Subcommittee of the
    American Academy of Neurology.** *Neurology.* 2013 (update 2018).
    <https://pubmed.ncbi.nlm.nih.gov/?term=Bhidayasiri+evidence-based+guideline+treatment+tardive+syndromes+American+Academy+of+Neurology>
66. Caroff SN, et al. **Treatment outcomes of patients with tardive
    dyskinesia and chronic schizophrenia.** *J Clin Psychiatry.* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Caroff+treatment+outcomes+patients+tardive+dyskinesia+chronic+schizophrenia>
67. Waln O, Jankovic J. **An update on tardive dyskinesia: from
    phenomenology to treatment.** *Tremor Other Hyperkinet Mov.* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Waln+Jankovic+update+tardive+dyskinesia+phenomenology+treatment>
68. Solmi M, et al. **Clinical risk factors for the development of tardive
    dyskinesia.** *J Neurol Sci.* 2018.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Solmi+clinical+risk+factors+development+tardive+dyskinesia>

## 13. QSP · 정량 모델링 방법론

69. Baron KT, et al. **mrgsolve: simulate from ODE-based population PK/PD
    and quantitative systems pharmacology models.**
    <https://mrgsolve.org/>
70. Musante CJ, et al. **Quantitative systems pharmacology: a case for
    disease models.** *Clin Pharmacol Ther.* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Musante+quantitative+systems+pharmacology+case+for+disease+models>
71. Peterson MC, Riggs MM. **FDA advisory meeting clinical pharmacology
    review utilizes a quantitative systems pharmacology (QSP) model.**
    *CPT Pharmacometrics Syst Pharmacol.* 2015.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Peterson+Riggs+FDA+advisory+meeting+clinical+pharmacology+review+quantitative+systems+pharmacology+model>
72. Gabrielsson J, Weiner D. **Non-compartmental analysis and turnover
    (indirect response) models in pharmacodynamics.**
    <https://pubmed.ncbi.nlm.nih.gov/?term=Gabrielsson+Weiner+turnover+indirect+response+models+pharmacodynamics+review>
