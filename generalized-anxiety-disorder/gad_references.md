# Generalized Anxiety Disorder QSP model — references

Every entry below was resolved live from NCBI PubMed (`esearch` + `esummary`) by [`mkrefs.py`](mkrefs.py); titles, journals, years, author lists and PMIDs are therefore machine-transcribed rather than recalled. Each entry carries the **intent** — what the model actually takes from it — so that a retrieved paper that does not match its intent is visible at a glance rather than hidden.

> **Five numbers in this repository were fitted; everything else was predicted.** The fitted five are the placebo-arm trajectory of Khan 2011 (week 1 −5.94, week 8 −11.10), and the week-4 endpoints of Rickels 2005 (placebo −8.4, pregabalin 300 mg −12.2, alprazolam 1.5 mg −10.9). Sections 10–16 are, with those exceptions, out-of-sample prediction targets.


## 1. The disorder itself, and why it is a gain and not a state

*Sustained, unpredictable-threat anxiety with a self-reinforcing worry engine. The BNST rather than the central amygdala is the node that matches the phenomenology, which is why the model gives E_amy a slow time constant and lets CBT rather than a drug move it.*

1. **Definition, epidemiology and course of GAD**  
   Tyrer P, Baldwin D. Generalised anxiety disorder. *Lancet* 2006. [PMID 17174708](https://pubmed.ncbi.nlm.nih.gov/17174708/)
2. **Lifetime prevalence and comorbidity in national surveys**  
   Kessler RC et al. Lifetime prevalence and age-of-onset distributions of DSM-IV disorders in the National Comorbidity Survey Replication. *Arch Gen Psychiatry* 2005. [PMID 15939837](https://pubmed.ncbi.nlm.nih.gov/15939837/)
3. **Worry as the diagnostic core; intolerance of uncertainty**  
   Dugas MJ et al. Generalized anxiety disorder: a preliminary test of a conceptual model. *Behav Res Ther* 1998. [PMID 9613027](https://pubmed.ncbi.nlm.nih.gov/9613027/)
4. **Sustained versus phasic threat: BNST and the anxiety/fear distinction**  
   Torrisi S et al. Resting-state connectivity of the bed nucleus of the stria terminalis and the central nucleus of the amygdala in clinical anxiety. *J Psychiatry Neurosci* 2019. [PMID 30964612](https://pubmed.ncbi.nlm.nih.gov/30964612/)
5. **Disability and economic burden of GAD**  
   Hoffman DL, Dukes EM, Wittchen HU. Human and economic burden of generalized anxiety disorder. *Depress Anxiety* 2008. [PMID 17146763](https://pubmed.ncbi.nlm.nih.gov/17146763/)
6. **Diagnostic boundary with major depression**  
   Kendler KS et al. Major depression and generalized anxiety disorder. Same genes, (partly) different environments?. *Arch Gen Psychiatry* 1992. [PMID 1514877](https://pubmed.ncbi.nlm.nih.gov/1514877/)


## 2. The threat circuit: E_amy in the numerator, C_pfc in the denominator

*The two circuit factors of Phi. Everything in this section constrains either the amygdala/BNST drive or the prefrontal control capacity.*

7. **Amygdala hyperreactivity in GAD**  
   Monk CS et al. Amygdala and ventrolateral prefrontal cortex activation to masked angry faces in children and adolescents with generalized anxiety disorder. *Arch Gen Psychiatry* 2008. [PMID 18458208](https://pubmed.ncbi.nlm.nih.gov/18458208/)
8. **Amygdala-prefrontal connectivity deficits in GAD**  
   Porta-Casteràs D et al. Prefrontal-amygdala connectivity in trait anxiety and generalized anxiety disorder: Testing the boundaries between healthy and pathological worries. *J Affect Disord* 2020. [PMID 32217221](https://pubmed.ncbi.nlm.nih.gov/32217221/)
9. **vmPFC top-down inhibition of amygdala output**  
   Cho JH, Deisseroth K, Bolshakov VY. Synaptic encoding of fear extinction in mPFC-amygdala circuits. *Neuron* 2013. [PMID 24290204](https://pubmed.ncbi.nlm.nih.gov/24290204/)
10. **BNST and unpredictable threat anticipation in humans**  
   Torrisi S et al. Extended amygdala connectivity changes during sustained shock anticipation. *Transl Psychiatry* 2018. [PMID 29382815](https://pubmed.ncbi.nlm.nih.gov/29382815/)
11. **Insula and interoception in anxiety**  
   Fermin ASR et al. Insula neuroanatomical networks predict interoceptive awareness. *Heliyon* 2023. [PMID 37520943](https://pubmed.ncbi.nlm.nih.gov/37520943/)
12. **Default mode network and worry**  
   Makovac E et al. The verbal nature of worry in generalized anxiety: Insights from the brain. *Neuroimage Clin* 2018. [PMID 29527493](https://pubmed.ncbi.nlm.nih.gov/29527493/)
13. **Startle potentiation to unpredictable threat as a GAD laboratory model**  
   Gorka SM et al. Acute orexin antagonism selectively modulates anticipatory anxiety in humans: implications for addiction and anxiety. *Transl Psychiatry* 2022. [PMID 35918313](https://pubmed.ncbi.nlm.nih.gov/35918313/)
14. **Structural changes in amygdala and prefrontal cortex in GAD**  
   Strawn JR et al. Neurostructural abnormalities in pediatric anxiety disorders. *J Anxiety Disord* 2015. [PMID 25890287](https://pubmed.ncbi.nlm.nih.gov/25890287/)


## 3. Serotonin: transporter occupancy, the 5-HT1A autoreceptor gate, and why the clock is slow

*The single most important calibration anchor in the model. SERT occupancy is >80% within days; the clinical effect is not. The gap is the autoreceptor.*

15. **SERT occupancy versus plasma concentration for SSRIs (the 80% threshold)**  
   Hart XM et al. Update Lessons from PET Imaging Part II: A Systematic Critical Review on Therapeutic Plasma Concentrations of Antidepressants. *Ther Drug Monit* 2024. [PMID 38287888](https://pubmed.ncbi.nlm.nih.gov/38287888/)
16. **Escitalopram SERT occupancy by [11C]DASB and its regional differences**  
   Kim E et al. Regional Differences in Serotonin Transporter Occupancy by Escitalopram: An [(11)C]DASB PK-PD Study. *Clin Pharmacokinet* 2017. [PMID 27557550](https://pubmed.ncbi.nlm.nih.gov/27557550/)
- **Dose-occupancy relationship for SSRIs at clinical doses** — *no PubMed record retrieved for this query; stated as an unresolved gap rather than cited from memory.*
17. **5-HT1A somatodendritic autoreceptor desensitisation with chronic SSRI**  
   Rainer Q et al. Functional status of somatodendritic serotonin 1A autoreceptor after long-term treatment with fluoxetine in a mouse model of anxiety/depression based on repeated corticosterone administration. *Mol Pharmacol* 2012. [PMID 22031471](https://pubmed.ncbi.nlm.nih.gov/22031471/)
18. **Time course of extracellular 5-HT rise with chronic SSRI (microdialysis)**  
   Tanda G, Frau R, Di Chiara G. Chronic desipramine and fluoxetine differentially affect extracellular dopamine in the rat prefrontal cortex. *Psychopharmacology (Berl)* 1996. [PMID 8888371](https://pubmed.ncbi.nlm.nih.gov/8888371/)
19. **Pindolol augmentation as a test of the autoreceptor hypothesis**  
   Liu Y et al. Is pindolol augmentation effective in depressed patients resistant to selective serotonin reuptake inhibitors? A systematic review and meta-analysis. *Hum Psychopharmacol* 2015. [PMID 25689398](https://pubmed.ncbi.nlm.nih.gov/25689398/)
20. **Delayed onset of antidepressant/anxiolytic action: the pattern to be explained**  
   Taylor MJ et al. Early onset of selective serotonin reuptake inhibitor antidepressant action: systematic review and meta-analysis. *Arch Gen Psychiatry* 2006. [PMID 17088502](https://pubmed.ncbi.nlm.nih.gov/17088502/)
21. **5-HT1A receptor binding in anxiety disorders**  
   Lanzenberger RR et al. Reduced serotonin-1A receptor binding in social anxiety disorder. *Biol Psychiatry* 2007. [PMID 16979141](https://pubmed.ncbi.nlm.nih.gov/16979141/)
22. **5-HT3 mediated nausea with SSRIs**  
   Jespersen C et al. Selective serotonin reuptake inhibitors for premenstrual syndrome and premenstrual dysphoric disorder. *Cochrane Database Syst Rev* 2024. [PMID 39140320](https://pubmed.ncbi.nlm.nih.gov/39140320/)


## 4. Noradrenaline and NET occupancy: the arm that helps and harms at once

*Venlafaxine's measured NET occupancy sets the SNRI dose-response, and the same NE rise produces the early activation that drives dropout.*

23. **Venlafaxine ER norepinephrine transporter occupancy in patients (PET)**  
   Arakawa R et al. Venlafaxine ER Blocks the Norepinephrine Transporter in the Brain of Patients with Major Depressive Disorder: a PET Study Using [18F]FMeNER-D2. *Int J Neuropsychopharmacol* 2019. [PMID 30649319](https://pubmed.ncbi.nlm.nih.gov/30649319/)
24. **Duloxetine SERT and NET occupancy**  
   Moriguchi S et al. A longitudinal PET study on changes in brain norepinephrine transporter availability following duloxetine treatment in major depressive disorder. *Int J Neuropsychopharmacol* 2025. [PMID 41035416](https://pubmed.ncbi.nlm.nih.gov/41035416/)
25. **Locus coeruleus and noradrenergic contribution to anxiety**  
   Bierwirth P, Stockhorst U. Role of noradrenergic arousal for fear extinction processes in rodents and humans. *Neurobiol Learn Mem* 2022. [PMID 35870717](https://pubmed.ncbi.nlm.nih.gov/35870717/)
26. **Prefrontal noradrenaline inverted-U and alpha1/alpha2 receptor actions**  
   Berridge CW, Arnsten AF. Psychostimulants and motivated behavior: arousal and cognition. *Neurosci Biobehav Rev* 2013. [PMID 23164814](https://pubmed.ncbi.nlm.nih.gov/23164814/)
- **SNRI early activation / jitteriness and blood pressure** — *no PubMed record retrieved for this query; stated as an unresolved gap rather than cited from memory.*


## 5. GABA-A: subunit pharmacology, benzodiazepine-site occupancy, tolerance and rebound

*Where I_gaba, R_a1, R_a2 and DEPEND come from. The measured lorazepam occupancy EC50 is the anchor; the alpha1/alpha5-versus-alpha2/alpha3 split is why sedation tolerates in days and anxiolysis does not.*

27. **Lorazepam plasma EC50 for GABA-A benzodiazepine-site occupancy**  
   Atack JR et al. Comparison of lorazepam [7-chloro-5-(2-chlorophenyl)-1,3-dihydro-3-hydroxy-2H-1,4-benzodiazepin-2-one] occupancy of rat brain gamma-aminobutyric acid(A) receptors measured using in vivo [3H]flumazenil (8-fluoro 5,6-dihydro-5-methyl-6-oxo-4H-imidazo[1,5-a][1,4]benzodiazepine-3-carboxylic acid ethyl ester) binding and [11C]flumazenil micro-positron emission tomography. *J Pharmacol Exp Ther* 2007. [PMID 17164474](https://pubmed.ncbi.nlm.nih.gov/17164474/)
28. **Benzodiazepine receptor occupancy at clinical doses in humans**  
   Malizia AL et al. Benzodiazepine site pharmacokinetic/pharmacodynamic quantification in man: direct measurement of drug occupancy and effects on the human brain in vivo. *Neuropharmacology* 1996. [PMID 9014164](https://pubmed.ncbi.nlm.nih.gov/9014164/)
29. **alpha2/alpha3 subunits mediate anxiolysis; alpha1 mediates sedation**  
   van Rijnsoever C et al. Requirement of alpha5-GABAA receptors for the development of tolerance to the sedative action of diazepam in mice. *J Neurosci* 2004. [PMID 15282283](https://pubmed.ncbi.nlm.nih.gov/15282283/)
30. **alpha1 subunit mediates the sedative action of diazepam**  
   van Rijnsoever C et al. Requirement of alpha5-GABAA receptors for the development of tolerance to the sedative action of diazepam in mice. *J Neurosci* 2004. [PMID 15282283](https://pubmed.ncbi.nlm.nih.gov/15282283/) *(also cited above)*
31. **alpha5-GABA-A receptors are required for tolerance to sedation**  
   van Rijnsoever C et al. Requirement of alpha5-GABAA receptors for the development of tolerance to the sedative action of diazepam in mice. *J Neurosci* 2004. [PMID 15282283](https://pubmed.ncbi.nlm.nih.gov/15282283/) *(also cited above)*
32. **Differential tolerance to sedative versus anxiolytic benzodiazepine effects**  
   Mediratta PK, Sharma KK, Rana J. Development of differential tolerance to the sedative and anti-stress effects of benzodiazepines. *Indian J Physiol Pharmacol* 2001. [PMID 11211563](https://pubmed.ncbi.nlm.nih.gov/11211563/)
33. **Cortical GABA concentration in anxiety disorders (1H-MRS)**  
   Maddock RJ, Smucny J. Transdiagnostic reduction in cortical choline-containing compounds in anxiety disorders: a (1)H-magnetic resonance spectroscopy meta-analysis. *Mol Psychiatry* 2025. [PMID 40913113](https://pubmed.ncbi.nlm.nih.gov/40913113/)
34. **Neurosteroids, allopregnanolone and GABA-A tonic inhibition**  
   Carver CM, Reddy DS. Neurosteroid interactions with synaptic and extrasynaptic GABA(A) receptors: regulation of subunit plasticity, phasic and tonic inhibition, and neuronal network excitability. *Psychopharmacology (Berl)* 2013. [PMID 24071826](https://pubmed.ncbi.nlm.nih.gov/24071826/)
35. **Benzodiazepine physical dependence and the withdrawal syndrome**  
   Heberlein A et al. [Benzodiazepine dependence: causalities and treatment options]. *Fortschr Neurol Psychiatr* 2009. [PMID 19101875](https://pubmed.ncbi.nlm.nih.gov/19101875/)


## 6. Glutamate, alpha2delta-1 and the pregabalin mechanism

*Where S_glu and its ~2 day trafficking delay come from.*

36. **alpha2delta-1 as the pregabalin binding site**  
   Field MJ et al. Identification of the alpha2-delta-1 subunit of voltage-dependent calcium channels as a molecular target for pain mediating the analgesic actions of pregabalin. *Proc Natl Acad Sci U S A* 2006. [PMID 17088553](https://pubmed.ncbi.nlm.nih.gov/17088553/)
37. **alpha2delta trafficking as the mechanism of delayed onset**  
   Hendrich J et al. Pharmacological disruption of calcium channel trafficking by the alpha2delta ligand gabapentin. *Proc Natl Acad Sci U S A* 2008. [PMID 18299583](https://pubmed.ncbi.nlm.nih.gov/18299583/)
38. **Pregabalin reduces excitatory neurotransmitter release**  
   Fink K et al. Inhibition of neuronal Ca(2+) influx by gabapentin and pregabalin in the human neocortex. *Neuropharmacology* 2002. [PMID 11804619](https://pubmed.ncbi.nlm.nih.gov/11804619/)
39. **Glutamate abnormalities in anxiety disorders**  
   Maddock RJ, Smucny J. Transdiagnostic reduction in cortical choline-containing compounds in anxiety disorders: a (1)H-magnetic resonance spectroscopy meta-analysis. *Mol Psychiatry* 2025. [PMID 40913113](https://pubmed.ncbi.nlm.nih.gov/40913113/) *(also cited above)*
40. **Stress upregulates alpha2delta-1 and presynaptic release probability**  
   Musazzi L, Treccani G, Popoli M. Functional and structural remodeling of glutamate synapses in prefrontal and frontal cortex induced by behavioral stress. *Front Psychiatry* 2015. [PMID 25964763](https://pubmed.ncbi.nlm.nih.gov/25964763/)


## 7. HPA axis and cortisol in anxiety

*The loop the model had to bound. Cortisol in GAD is modestly raised, not grossly so, and the model is held to that.*

41. **Cortisol in generalized anxiety disorder**  
   Rosnick CB et al. Association of cortisol with neuropsychological assessment in older adults with generalized anxiety disorder. *Aging Ment Health* 2013. [PMID 23336532](https://pubmed.ncbi.nlm.nih.gov/23336532/)
42. **CRH, CRHR1 and anxiety behaviour**  
   Marcinkiewcz CA et al. Serotonin engages an anxiety and fear-promoting circuit in the extended amygdala. *Nature* 2016. [PMID 27556938](https://pubmed.ncbi.nlm.nih.gov/27556938/)
43. **Glucocorticoid receptor resistance under chronic stress**  
   Cohen S et al. Chronic stress, glucocorticoid receptor resistance, inflammation, and disease risk. *Proc Natl Acad Sci U S A* 2012. [PMID 22474371](https://pubmed.ncbi.nlm.nih.gov/22474371/)
44. **Dexamethasone suppression and HPA feedback in anxiety**  
   Schweizer EE et al. The dexamethasone suppression test in generalised anxiety disorder. *Br J Psychiatry* 1986. [PMID 3779298](https://pubmed.ncbi.nlm.nih.gov/3779298/)
45. **Cortisol and prefrontal function**  
   Woodcock EA et al. Pharmacological stress impairs working memory performance and attenuates dorsolateral prefrontal cortex glutamate modulation. *Neuroimage* 2019. [PMID 30458306](https://pubmed.ncbi.nlm.nih.gov/30458306/)


## 8. Neuroplasticity, BDNF and extinction learning

*The mechanism by which a transporter block becomes a change in prefrontal control weeks later, and the reason CBT's gain outlasts the drug's.*

46. **BDNF and antidepressant-induced plasticity**  
   Castrén E. Neuronal network plasticity and recovery from depression. *JAMA Psychiatry* 2013. [PMID 23842648](https://pubmed.ncbi.nlm.nih.gov/23842648/)
47. **Fear extinction learning: circuitry and consolidation**  
   Milad MR, Rosenbaum BL, Simon NM. Neuroscience of fear extinction: implications for assessment and treatment of fear-based and anxiety related disorders. *Behav Res Ther* 2014. [PMID 25204715](https://pubmed.ncbi.nlm.nih.gov/25204715/)
48. **BDNF Val66Met and impaired extinction**  
   Felmingham KL et al. The BDNF Val66Met polymorphism moderates the relationship between Posttraumatic Stress Disorder and fear extinction learning. *Psychoneuroendocrinology* 2018. [PMID 29550677](https://pubmed.ncbi.nlm.nih.gov/29550677/)
- **Chronic stress causes opposite dendritic remodelling in PFC and amygdala** — *no PubMed record retrieved for this query; stated as an unresolved gap rather than cited from memory.*
49. **Inflammation, kynurenine and reduced BDNF**  
   Gibney SM et al. Poly I:C-induced activation of the immune response is accompanied by depression and anxiety-like behaviours, kynurenine pathway activation and reduced BDNF expression. *Brain Behav Immun* 2013. [PMID 23201589](https://pubmed.ncbi.nlm.nih.gov/23201589/)


## 9. Autonomic function, heart-rate variability and sleep

*The somatic half of HAM-A, and the sleep-to-amygdala loop that is inside the disease rather than downstream of it.*

50. **Reduced heart rate variability in anxiety disorders (meta-analysis)**  
   Chalmers JA et al. Anxiety Disorders are Associated with Reduced Heart Rate Variability: A Meta-Analysis. *Front Psychiatry* 2014. [PMID 25071612](https://pubmed.ncbi.nlm.nih.gov/25071612/)
51. **Autonomic inflexibility in GAD**  
   Thayer JF, Friedman BH, Borkovec TD. Autonomic characteristics of generalized anxiety disorder and worry. *Biol Psychiatry* 1996. [PMID 8645772](https://pubmed.ncbi.nlm.nih.gov/8645772/)
52. **Insomnia and anxiety: bidirectional relationship**  
   Pigeon WR, Bishop TM, Krueger KM. Insomnia as a Precipitating Factor in New Onset Mental Illness: a Systematic Review of Recent Findings. *Curr Psychiatry Rep* 2017. [PMID 28616860](https://pubmed.ncbi.nlm.nih.gov/28616860/)
53. **Sleep deprivation amplifies amygdala reactivity**  
   Yoo SS et al. The human emotional brain without sleep--a prefrontal amygdala disconnect. *Curr Biol* 2007. [PMID 17956744](https://pubmed.ncbi.nlm.nih.gov/17956744/)
54. **Anxiety and cardiovascular risk**  
   Roest AM et al. Anxiety and risk of incident coronary heart disease: a meta-analysis. *J Am Coll Cardiol* 2010. [PMID 20620715](https://pubmed.ncbi.nlm.nih.gov/20620715/)


## 10. SSRI and SNRI randomised trials in GAD

*Out-of-sample prediction targets. None of these were fitted.*

55. **Escitalopram fixed-dose efficacy in GAD**  
   Baldwin DS, Huusom AK, Maehlum E. Escitalopram and paroxetine in the treatment of generalised anxiety disorder: randomised, placebo-controlled, double-blind study. *Br J Psychiatry* 2006. [PMID 16946363](https://pubmed.ncbi.nlm.nih.gov/16946363/)
- **Paroxetine in GAD** — *no PubMed record retrieved for this query; stated as an unresolved gap rather than cited from memory.*
- **Sertraline in GAD** — *no PubMed record retrieved for this query; stated as an unresolved gap rather than cited from memory.*
56. **Venlafaxine ER in GAD: short and long term**  
   Gelenberg AJ et al. Efficacy of venlafaxine extended-release capsules in nondepressed outpatients with generalized anxiety disorder: A 6-month randomized controlled trial. *JAMA* 2000. [PMID 10865302](https://pubmed.ncbi.nlm.nih.gov/10865302/)
57. **Venlafaxine ER in older adults with GAD**  
   Katz IR et al. Venlafaxine ER as a treatment for generalized anxiety disorder in older adults: pooled analysis of five randomized placebo-controlled clinical trials. *J Am Geriatr Soc* 2002. [PMID 12028242](https://pubmed.ncbi.nlm.nih.gov/12028242/)
58. **Duloxetine in GAD**  
   Carter NJ, McCormack PL. Duloxetine: a review of its use in the treatment of generalized anxiety disorder. *CNS Drugs* 2009. [PMID 19480470](https://pubmed.ncbi.nlm.nih.gov/19480470/)
59. **Duloxetine versus venlafaxine in GAD**  
   Slee A et al. Pharmacological treatments for generalised anxiety disorder: a systematic review and network meta-analysis. *Lancet* 2019. [PMID 30712879](https://pubmed.ncbi.nlm.nih.gov/30712879/)
60. **Antidepressant dose-response is flat for SSRIs**  
   Furukawa TA et al. Optimal dose of selective serotonin reuptake inhibitors, venlafaxine, and mirtazapine in major depression: a systematic review and dose-response meta-analysis. *Lancet Psychiatry* 2019. [PMID 31178367](https://pubmed.ncbi.nlm.nih.gov/31178367/)


## 11. Pregabalin randomised trials in GAD

*Rickels 2005 supplies two of the five fitted numbers (pregabalin 300 mg and the placebo arm at week 4). The 450 and 600 mg arms, and the psychic/somatic split, are predictions.*

61. **Rickels 2005: pregabalin 300/450/600 versus alprazolam and placebo**  
   Rickels K et al. Pregabalin for treatment of generalized anxiety disorder: a 4-week, multicenter, double-blind, placebo-controlled trial of pregabalin and alprazolam. *Arch Gen Psychiatry* 2005. [PMID 16143734](https://pubmed.ncbi.nlm.nih.gov/16143734/)
62. **Pregabalin review: onset within one week, psychic and somatic**  
   Frampton JE. Pregabalin: a review of its use in adults with generalized anxiety disorder. *CNS Drugs* 2014. [PMID 25149863](https://pubmed.ncbi.nlm.nih.gov/25149863/)
63. **Pregabalin efficacy and safety profile in generalized anxiety**  
   Owen RT. Pregabalin: its efficacy, safety and tolerability profile in generalized anxiety. *Drugs Today (Barc)* 2007. [PMID 17940637](https://pubmed.ncbi.nlm.nih.gov/17940637/)
64. **Pregabalin in elderly patients with GAD**  
   Montgomery SA et al. Early improvement with pregabalin predicts endpoint response in patients with generalized anxiety disorder: an integrated and predictive data analysis. *Int Clin Psychopharmacol* 2017. [PMID 27583543](https://pubmed.ncbi.nlm.nih.gov/27583543/)
65. **Pregabalin relapse prevention in GAD**  
   Lam RW. Generalized anxiety disorder: how to treat, and for how long?. *Int J Psychiatry Clin Pract* 2006. [PMID 24931538](https://pubmed.ncbi.nlm.nih.gov/24931538/)
66. **Pregabalin augmentation of SSRI/SNRI partial responders**  
   Rickels K et al. Adjunctive therapy with pregabalin in generalized anxiety disorder patients with partial response to SSRI or SNRI treatment. *Int Clin Psychopharmacol* 2012. [PMID 22302014](https://pubmed.ncbi.nlm.nih.gov/22302014/)


## 12. Benzodiazepines: efficacy, tolerance, dependence and discontinuation

*Rickels 2005's alprazolam arm supplies the fifth fitted number. Everything about tolerance and rebound is a prediction constrained by this literature.*

67. **Alprazolam clinical pharmacology, efficacy and behavioural toxicity**  
   Verster JC, Volkerts ER. Clinical pharmacology, clinical efficacy, and behavioral toxicity of alprazolam: a review of the literature. *CNS Drug Rev* 2004. [PMID 14978513](https://pubmed.ncbi.nlm.nih.gov/14978513/)
68. **Benzodiazepines versus antidepressants in GAD**  
   Bandelow B et al. Efficacy of treatments for anxiety disorders: a meta-analysis. *Int Clin Psychopharmacol* 2015. [PMID 25932596](https://pubmed.ncbi.nlm.nih.gov/25932596/)
69. **Long-term benzodiazepine use: benefits and harms**  
   Lader MH. Limitations on the use of benzodiazepines in anxiety and insomnia: are they justified?. *Eur Neuropsychopharmacol* 1999. [PMID 10622686](https://pubmed.ncbi.nlm.nih.gov/10622686/)
70. **Benzodiazepine discontinuation strategies and rebound anxiety**  
   Kasper S et al. Pregabalin long-term treatment and assessment of discontinuation in patients with generalized anxiety disorder. *Int J Neuropsychopharmacol* 2014. [PMID 24351233](https://pubmed.ncbi.nlm.nih.gov/24351233/)
71. **Benzodiazepines, falls and fractures in older adults**  
   Westaway K et al. Combination psychotropic medicine use in older adults and risk of hip fracture. *Aust Prescr* 2019. [PMID 31363307](https://pubmed.ncbi.nlm.nih.gov/31363307/)
- **Cognitive and psychomotor impairment with benzodiazepines** — *no PubMed record retrieved for this query; stated as an unresolved gap rather than cited from memory.*


## 13. Quetiapine, buspirone, hydroxyzine and the rest of the armamentarium

*The quetiapine XR dose-response is the model's sharpest out-of-sample test, because the reported curve is non-monotonic.*

72. **Khan 2011: quetiapine XR 50/150/300 mg monotherapy in GAD**  
   Khan A et al. A randomized, double-blind study of once-daily extended release quetiapine fumarate (quetiapine XR) monotherapy in patients with generalized anxiety disorder. *J Clin Psychopharmacol* 2011. [PMID 21694613](https://pubmed.ncbi.nlm.nih.gov/21694613/)
73. **Quetiapine XR in GAD: meta-analysis**  
   Maneeton N et al. Quetiapine monotherapy in acute treatment of generalized anxiety disorder: a systematic review and meta-analysis of randomized controlled trials. *Drug Des Devel Ther* 2016. [PMID 26834458](https://pubmed.ncbi.nlm.nih.gov/26834458/)
74. **Buspirone in generalized anxiety**  
   Bandelow B, Michaelis S, Wedekind D. Treatment of anxiety disorders. *Dialogues Clin Neurosci* 2017. [PMID 28867934](https://pubmed.ncbi.nlm.nih.gov/28867934/)
75. **Azapirones for generalised anxiety disorder (Cochrane)**  
   Chessick CA et al. Azapirones for generalized anxiety disorder. *Cochrane Database Syst Rev* 2006. [PMID 16856115](https://pubmed.ncbi.nlm.nih.gov/16856115/)
76. **Hydroxyzine for generalised anxiety disorder**  
   Gale CK, Millichamp J. Generalised anxiety disorder. *BMJ Clin Evid* 2011. [PMID 22030083](https://pubmed.ncbi.nlm.nih.gov/22030083/)
77. **Agomelatine in GAD**  
   Slee A et al. Pharmacological treatments for generalised anxiety disorder: a systematic review and network meta-analysis. *Lancet* 2019. [PMID 30712879](https://pubmed.ncbi.nlm.nih.gov/30712879/) *(also cited above)*
78. **Vortioxetine in GAD**  
   Rothschild AJ et al. Vortioxetine (Lu AA21004) 5 mg in generalized anxiety disorder: results of an 8-week randomized, double-blind, placebo-controlled clinical trial in the United States. *Eur Neuropsychopharmacol* 2012. [PMID 22901736](https://pubmed.ncbi.nlm.nih.gov/22901736/)
79. **Antipsychotic augmentation in refractory GAD**  
   Schiele MA et al. Integrative Systematic Review on Pharmacological, Psychotherapeutic, and Neurostimulatory Treatment Options in Treatment-Resistant Anxiety Disorders. *Psychother Psychosom* 2026. [PMID 40946318](https://pubmed.ncbi.nlm.nih.gov/40946318/)


## 14. Relapse prevention and long-term treatment

*The Allgulander randomised-withdrawal design is the money test: 19% versus 56% relapse, predicted from an acute-phase calibration.*

80. **Allgulander 2006: escitalopram relapse prevention in GAD**  
   Allgulander C, Florea I, Huusom AK. Prevention of relapse in generalized anxiety disorder by escitalopram treatment. *Int J Neuropsychopharmacol* 2006. [PMID 16316482](https://pubmed.ncbi.nlm.nih.gov/16316482/)
81. **Duloxetine relapse prevention in GAD**  
   Bodkin JA et al. Predictors of relapse in a study of duloxetine treatment for patients with generalized anxiety disorder. *Hum Psychopharmacol* 2011. [PMID 21678494](https://pubmed.ncbi.nlm.nih.gov/21678494/)
82. **Randomised-withdrawal designs in anxiety disorders: meta-analysis**  
   Batelaan NM et al. Risk of relapse after antidepressant discontinuation in anxiety disorders, obsessive-compulsive disorder, and post-traumatic stress disorder: systematic review and meta-analysis of relapse prevention trials. *BMJ* 2017. [PMID 28903922](https://pubmed.ncbi.nlm.nih.gov/28903922/)
83. **How long should treatment continue after remission**  
   Adu MK et al. Comparing Email Versus Text Messaging as Delivery Platforms for Supporting Patients With Major Depressive Disorder: Noninferiority Randomized Controlled Trial. *JMIR Form Res* 2024. [PMID 39250182](https://pubmed.ncbi.nlm.nih.gov/39250182/)
84. **Discontinuation symptoms after stopping SSRIs and SNRIs**  
   Cipriani A et al. Comparative efficacy and acceptability of 21 antidepressant drugs for the acute treatment of adults with major depressive disorder: a systematic review and network meta-analysis. *Lancet* 2018. [PMID 29477251](https://pubmed.ncbi.nlm.nih.gov/29477251/)


## 15. Placebo response, assay sensitivity and trial methodology

*The fifth clock. The model treats expectancy as real top-down control, which is why assay sensitivity falls out of the algebra.*

85. **Placebo response in anxiety disorder trials**  
   Stimpfl JN, Mills JA, Strawn JR. Pharmacologic predictors of benzodiazepine response trajectory in anxiety disorders: a Bayesian hierarchical modeling meta-analysis. *CNS Spectr* 2023. [PMID 34593077](https://pubmed.ncbi.nlm.nih.gov/34593077/)
86. **Rising placebo response and failed trials in psychiatry**  
   Arnold R et al. Predictors of the placebo response in a nutraceutical randomized controlled trial for depression. *J Integr Med* 2024. [PMID 38331652](https://pubmed.ncbi.nlm.nih.gov/38331652/)
87. **Number of study visits predicts placebo response**  
   Elliott J et al. Pharmacologic treatment of attention deficit hyperactivity disorder in adults: A systematic review and network meta-analysis. *PLoS One* 2020. [PMID 33085721](https://pubmed.ncbi.nlm.nih.gov/33085721/)
- **Regression to the mean and baseline severity inflation in trials** — *no PubMed record retrieved for this query; stated as an unresolved gap rather than cited from memory.*
88. **Expectancy mechanisms of placebo analgesia and anxiolysis**  
   Evers AWM et al. Implications of Placebo and Nocebo Effects for Clinical Practice: Expert Consensus. *Psychother Psychosom* 2018. [PMID 29895014](https://pubmed.ncbi.nlm.nih.gov/29895014/)
- **Baseline severity and drug-placebo difference** — *no PubMed record retrieved for this query; stated as an unresolved gap rather than cited from memory.*


## 16. Network meta-analyses, comparative efficacy and guidelines

*The benchmark the whole set of predictions is scored against.*

89. **Network meta-analysis of first-line drugs for GAD**  
   He H et al. Comparative efficacy and acceptability of first-line drugs for the acute treatment of generalized anxiety disorder in adults: A network meta-analysis. *J Psychiatr Res* 2019. [PMID 31473564](https://pubmed.ncbi.nlm.nih.gov/31473564/)
90. **Remission rates and tolerability network meta-analysis in GAD**  
   Kong W et al. Comparative Remission Rates and Tolerability of Drugs for Generalised Anxiety Disorder: A Systematic Review and Network Meta-analysis of Double-Blind Randomized Controlled Trials. *Front Pharmacol* 2020. [PMID 33343351](https://pubmed.ncbi.nlm.nih.gov/33343351/)
91. **Anxiolytic drugs across anxiety disorders: network meta-analysis**  
   Slee A et al. Pharmacological treatments for generalised anxiety disorder: a systematic review and network meta-analysis. *Lancet* 2019. [PMID 30712879](https://pubmed.ncbi.nlm.nih.gov/30712879/) *(also cited above)*
92. **WFSBP guidelines for anxiety disorders**  
   Bandelow B et al. World Federation of Societies of Biological Psychiatry (WFSBP) guidelines for treatment of anxiety, obsessive-compulsive and posttraumatic stress disorders - Version 3. Part I: Anxiety disorders. *World J Biol Psychiatry* 2023. [PMID 35900161](https://pubmed.ncbi.nlm.nih.gov/35900161/)
93. **NICE / national guidance on GAD management**  
   Baldwin DS et al. Evidence-based pharmacological treatment of anxiety disorders, post-traumatic stress disorder and obsessive-compulsive disorder: a revision of the 2005 guidelines from the British Association for Psychopharmacology. *J Psychopharmacol* 2014. [PMID 24713617](https://pubmed.ncbi.nlm.nih.gov/24713617/)


## 17. Psychological treatment

*CBT enters the model at three places at once and is the only intervention that moves E_amy.*

- **CBT for GAD: meta-analysis of effect size** — *no PubMed record retrieved for this query; stated as an unresolved gap rather than cited from memory.*
94. **CBT for pathological worry**  
   Covin R et al. A meta-analysis of CBT for pathological worry among clients with GAD. *J Anxiety Disord* 2008. [PMID 17321717](https://pubmed.ncbi.nlm.nih.gov/17321717/)
95. **Combined psychotherapy and pharmacotherapy in anxiety disorders**  
   Zhou X et al. Comparative efficacy and acceptability of antidepressants, psychotherapies, and their combination for acute treatment of children and adolescents with depressive disorder: a systematic review and network meta-analysis. *Lancet Psychiatry* 2020. [PMID 32563306](https://pubmed.ncbi.nlm.nih.gov/32563306/)
96. **Internet-delivered CBT for anxiety**  
   Carlbring P et al. Internet-based vs. face-to-face cognitive behavior therapy for psychiatric and somatic disorders: an updated systematic review and meta-analysis. *Cogn Behav Ther* 2018. [PMID 29215315](https://pubmed.ncbi.nlm.nih.gov/29215315/)
97. **Long-term durability of CBT gains in anxiety**  
   van Dis EAM et al. Long-term Outcomes of Cognitive Behavioral Therapy for Anxiety-Related Disorders: A Systematic Review and Meta-analysis. *JAMA Psychiatry* 2020. [PMID 31758858](https://pubmed.ncbi.nlm.nih.gov/31758858/)
98. **Exercise for anxiety disorders**  
   Banyard H et al. The Effects of Aerobic and Resistance Exercise on Depression and Anxiety: Systematic Review With Meta-Analysis. *Int J Ment Health Nurs* 2025. [PMID 40432290](https://pubmed.ncbi.nlm.nih.gov/40432290/)


## 18. The Hamilton Anxiety Rating Scale as a measuring instrument

*The readout is not transparent. Its psychic/somatic structure and its ceiling are load-bearing parts of the model.*

99. **The Hamilton Anxiety Rating Scale: original description**  
   HAMILTON M. The assessment of anxiety states by rating. *Br J Med Psychol* 1959. [PMID 13638508](https://pubmed.ncbi.nlm.nih.gov/13638508/)
100. **HAM-A psychometric properties and factor structure**  
   Maier W et al. The Hamilton Anxiety Scale: reliability, validity and sensitivity to change in anxiety and depressive disorders. *J Affect Disord* 1988. [PMID 2963053](https://pubmed.ncbi.nlm.nih.gov/2963053/)
- **Clinically meaningful change on the HAM-A** — *no PubMed record retrieved for this query; stated as an unresolved gap rather than cited from memory.*
- **GAD-7 development and validation** — *no PubMed record retrieved for this query; stated as an unresolved gap rather than cited from memory.*
101. **Penn State Worry Questionnaire**  
   Meyer TJ et al. Development and validation of the Penn State Worry Questionnaire. *Behav Res Ther* 1990. [PMID 2076086](https://pubmed.ncbi.nlm.nih.gov/2076086/)
102. **Limits of change-score endpoints in psychiatry trials**  
   Nakonezny PA et al. Evaluation of anhedonia with the Snaith-Hamilton Pleasure Scale (SHAPS) in adult outpatients with major depressive disorder. *J Psychiatr Res* 2015. [PMID 25864641](https://pubmed.ncbi.nlm.nih.gov/25864641/)


## 19. Clinical pharmacokinetics of the seven modelled drugs

*Every PK parameter in the model traces to this section.*

103. **Escitalopram clinical pharmacokinetics**  
   Rao N. The clinical pharmacokinetics of escitalopram. *Clin Pharmacokinet* 2007. [PMID 17375980](https://pubmed.ncbi.nlm.nih.gov/17375980/)
104. **CYP2C19 genotype and escitalopram exposure**  
   Faraj P et al. Combined effect of CYP2C19 and CYP2D6 genotypes on escitalopram serum concentration and its metabolic ratio in a European patient population. *Br J Clin Pharmacol* 2024. [PMID 38925553](https://pubmed.ncbi.nlm.nih.gov/38925553/)
105. **Venlafaxine and O-desmethylvenlafaxine pharmacokinetics; CYP2D6**  
   Fukuda T et al. The impact of the CYP2D6 and CYP2C19 genotypes on venlafaxine pharmacokinetics in a Japanese population. *Eur J Clin Pharmacol* 2000. [PMID 10877013](https://pubmed.ncbi.nlm.nih.gov/10877013/)
106. **Duloxetine pharmacokinetics**  
   Knadler MP et al. Duloxetine: clinical pharmacokinetics and drug interactions. *Clin Pharmacokinet* 2011. [PMID 21366359](https://pubmed.ncbi.nlm.nih.gov/21366359/)
107. **Pregabalin pharmacokinetics and renal dose adjustment**  
   Randinitis EJ et al. Pharmacokinetics of pregabalin in subjects with various degrees of renal function. *J Clin Pharmacol* 2003. [PMID 12638396](https://pubmed.ncbi.nlm.nih.gov/12638396/)
108. **Lorazepam pharmacokinetics**  
   Greenblatt DJ. Clinical pharmacokinetics of oxazepam and lorazepam. *Clin Pharmacokinet* 1981. [PMID 6111408](https://pubmed.ncbi.nlm.nih.gov/6111408/)
109. **Alprazolam pharmacokinetics**  
   Kaplan GB et al. Single-dose pharmacokinetics and pharmacodynamics of alprazolam in elderly and young subjects. *J Clin Pharmacol* 1998. [PMID 9597554](https://pubmed.ncbi.nlm.nih.gov/9597554/)
110. **Buspirone pharmacokinetics and 1-PP metabolite**  
   Mahmood I, Sahajwalla C. Clinical pharmacokinetics and pharmacodynamics of buspirone, an anxiolytic drug. *Clin Pharmacokinet* 1999. [PMID 10320950](https://pubmed.ncbi.nlm.nih.gov/10320950/)
111. **Quetiapine and norquetiapine pharmacokinetics**  
   Nikisch G et al. Quetiapine and norquetiapine in plasma and cerebrospinal fluid of schizophrenic patients treated with quetiapine: correlations to clinical outcome and HVA, 5-HIAA, and MHPG in CSF. *J Clin Psychopharmacol* 2010. [PMID 20814316](https://pubmed.ncbi.nlm.nih.gov/20814316/)
112. **Therapeutic reference ranges for antidepressants (TDM consensus)**  
   Hiemke C et al. Consensus Guidelines for Therapeutic Drug Monitoring in Neuropsychopharmacology: Update 2017. *Pharmacopsychiatry* 2018. [PMID 28910830](https://pubmed.ncbi.nlm.nih.gov/28910830/)


## 20. Quantitative systems pharmacology and PK/PD methodology

*The modelling machinery itself.*

- **Quantitative systems pharmacology: scope and role in drug development** — *no PubMed record retrieved for this query; stated as an unresolved gap rather than cited from memory.*
113. **mrgsolve for ODE-based PK/PD simulation in R**  
   Elmokadem A, Riggs MM, Baron KT. Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial. *CPT Pharmacometrics Syst Pharmacol* 2019. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
114. **Model-informed drug development in CNS**  
   Belov A et al. Opportunities and challenges for applying model-informed drug development approaches to gene therapies. *CPT Pharmacometrics Syst Pharmacol* 2021. [PMID 33608998](https://pubmed.ncbi.nlm.nih.gov/33608998/)
115. **Receptor occupancy modelling and the link to clinical effect**  
   Sawant-Basak A et al. Investigating CNS distribution of PF-05212377, a P-glycoprotein substrate, by translation of 5-HT(6) receptor occupancy from non-human primates to humans. *Biopharm Drug Dispos* 2023. [PMID 36825693](https://pubmed.ncbi.nlm.nih.gov/36825693/)
116. **Clinical trial simulation with placebo and dropout models**  
   Knebel W et al. Modeling and simulation of the exposure-response and dropout pattern of guanfacine extended-release in pediatric patients with ADHD. *J Pharmacokinet Pharmacodyn* 2015. [PMID 25373474](https://pubmed.ncbi.nlm.nih.gov/25373474/)
117. **Disease progression modelling in neuropsychiatry**  
   Gomeni R et al. Modeling Alzheimer's disease progression using the disease system analysis approach. *Alzheimers Dement* 2012. [PMID 21782528](https://pubmed.ncbi.nlm.nih.gov/21782528/)
118. **Item-response / latent-variable models for psychiatric rating scales**  
   Hamdan A et al. Longitudinal Analysis of Natural History Progression of Rare and Ultra-Rare Cerebellar Ataxias Using Item Response Theory. *Clin Pharmacol Ther* 2024. [PMID 39403821](https://pubmed.ncbi.nlm.nih.gov/39403821/)


---

*118 bibliography entries, 113 distinct PubMed records. Regenerate with `python3 mkrefs.py` (add `--refresh` to re-run every query).*
