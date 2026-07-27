# Chronic Insomnia Disorder — 참고문헌 / References

만성 불면장애 QSP 모델(`ins_qsp_model.dot` · `ins_mrgsolve_model.R` · `ins_shiny_app.R`)의
근거 문헌입니다. 아래 **105편**은 모두 PubMed E-utilities로 조회하여 PMID·저널·연도를
확인한 것이며, 링크는 해당 PMID로 직접 연결됩니다.

문헌은 모델의 구성 순서를 따라 배열되어 있습니다:
질환 정의 → 두 과정 모형 → 항상성/일주기 기전 → flip-flop 스위치 → 과각성 →
3-P 유지 고리 → CBT-I → 약물별 PK/PD → 외인성 조절인자 → 결과 및 정량 모델링.

> ⚠️ 본 목록은 교육·연구 목적의 모델 문서화를 위한 것이며, 진료 지침을 대체하지 않습니다.

---

## A. 역학·정의·진단

1. **Morin CM et al.** Insomnia disorder. *Nat Rev Dis Primers* 2015. [PMID 27189779](https://pubmed.ncbi.nlm.nih.gov/27189779/)
2. **Ohayon MM et al.** Epidemiology of insomnia: what we know and what we still need to learn. *Sleep Med Rev* 2002. [PMID 12531146](https://pubmed.ncbi.nlm.nih.gov/12531146/)
3. **Morin CM et al.** The Insomnia Severity Index: psychometric indicators to detect insomnia cases and evaluate treatment response. *Sleep* 2011. [PMID 21532953](https://pubmed.ncbi.nlm.nih.gov/21532953/)
4. **Buysse DJ et al.** The Pittsburgh Sleep Quality Index: a new instrument for psychiatric practice and research. *Psychiatry Res* 1989. [PMID 2748771](https://pubmed.ncbi.nlm.nih.gov/2748771/)
5. **Riemann D et al.** The European Insomnia Guideline: An update on the diagnosis and treatment of insomnia 2023. *J Sleep Res* 2023. [PMID 38016484](https://pubmed.ncbi.nlm.nih.gov/38016484/)
6. **Edinger JD et al.** Behavioral and psychological treatments for chronic insomnia disorder in adults: an American Academy of Sleep Medicine clinical practice guideline. *J Clin Sleep Med* 2021. [PMID 33164742](https://pubmed.ncbi.nlm.nih.gov/33164742/)
7. **Sateia MJ et al.** Clinical Practice Guideline for the Pharmacologic Treatment of Chronic Insomnia in Adults: An American Academy of Sleep Medicine Clinical Practice Guideline. *J Clin Sleep Med* 2017. [PMID 27998379](https://pubmed.ncbi.nlm.nih.gov/27998379/)
8. **Qaseem A et al.** Management of Chronic Insomnia Disorder in Adults: A Clinical Practice Guideline From the American College of Physicians. *Ann Intern Med* 2016. [PMID 27136449](https://pubmed.ncbi.nlm.nih.gov/27136449/)
9. **Bastien CH et al.** Validation of the Insomnia Severity Index as an outcome measure for insomnia research. *Sleep Med* 2001. [PMID 11438246](https://pubmed.ncbi.nlm.nih.gov/11438246/)
10. **Buysse DJ et al.** Insomnia. *JAMA* 2013. [PMID 23423416](https://pubmed.ncbi.nlm.nih.gov/23423416/)

## B. 두 과정 모형 (Process S x C)

11. **Borbély AA et al.** The two-process model of sleep regulation: a reappraisal. *J Sleep Res* 2016. [PMID 26762182](https://pubmed.ncbi.nlm.nih.gov/26762182/)
12. **Borbély AA et al.** A two process model of sleep regulation. *Hum Neurobiol* 1982. [PMID 7185792](https://pubmed.ncbi.nlm.nih.gov/7185792/)
13. **Daan S et al.** Timing of human sleep: recovery process gated by a circadian pacemaker. *Am J Physiol* 1984. [PMID 6696142](https://pubmed.ncbi.nlm.nih.gov/6696142/)
14. **Dijk DJ et al.** Contribution of the circadian pacemaker and the sleep homeostat to sleep propensity, sleep structure, electroencephalographic slow waves, and sleep spindle activity in humans. *J Neurosci* 1995. [PMID 7751928](https://pubmed.ncbi.nlm.nih.gov/7751928/)
15. **Phillips AJ et al.** A quantitative model of sleep-wake dynamics based on the physiology of the brainstem ascending arousal system. *J Biol Rhythms* 2007. [PMID 17440218](https://pubmed.ncbi.nlm.nih.gov/17440218/)
16. **Booth V et al.** Physiologically-based modeling of sleep-wake regulatory networks. *Math Biosci* 2014. [PMID 24530893](https://pubmed.ncbi.nlm.nih.gov/24530893/)
17. **Achermann P et al.** Simulation of human sleep: ultradian dynamics of electroencephalographic slow-wave activity. *J Biol Rhythms* 1990. [PMID 2133124](https://pubmed.ncbi.nlm.nih.gov/2133124/)
18. **Skeldon AC et al.** Mathematical models for sleep-wake dynamics: comparison of the two-process model and a mutual inhibition neuronal model. *PLoS One* 2014. [PMID 25084361](https://pubmed.ncbi.nlm.nih.gov/25084361/)

## C. 항상성 수면압 · 아데노신

19. **Porkka-Heiskanen T et al.** Adenosine: a mediator of the sleep-inducing effects of prolonged wakefulness. *Science* 1997. [PMID 9157887](https://pubmed.ncbi.nlm.nih.gov/9157887/)
20. **Huang ZL et al.** Prostaglandins and adenosine in the regulation of sleep and wakefulness. *Curr Opin Pharmacol* 2007. [PMID 17129762](https://pubmed.ncbi.nlm.nih.gov/17129762/)
21. **Bjorness TE et al.** Adenosine and sleep. *Curr Neuropharmacol* 2009. [PMID 20190965](https://pubmed.ncbi.nlm.nih.gov/20190965/)
22. **Rétey JV et al.** A genetic variation in the adenosine A2A receptor gene (ADORA2A) contributes to individual sensitivity to caffeine effects on sleep. *Clin Pharmacol Ther* 2007. [PMID 17329997](https://pubmed.ncbi.nlm.nih.gov/17329997/)
23. **Tononi G et al.** Sleep function and synaptic homeostasis. *Sleep Med Rev* 2006. [PMID 16376591](https://pubmed.ncbi.nlm.nih.gov/16376591/)

## D. 일주기 시계 · 멜라토닌

24. **Reppert SM et al.** Coordination of circadian timing in mammals. *Nature* 2002. [PMID 12198538](https://pubmed.ncbi.nlm.nih.gov/12198538/)
25. **Takahashi JS et al.** Transcriptional architecture of the mammalian circadian clock. *Nat Rev Genet* 2017. [PMID 27990019](https://pubmed.ncbi.nlm.nih.gov/27990019/)
26. **Berson DM et al.** Phototransduction by retinal ganglion cells that set the circadian clock. *Science* 2002. [PMID 11834835](https://pubmed.ncbi.nlm.nih.gov/11834835/)
27. **Lewy AJ et al.** Exogenous melatonin's phase-shifting effects on the endogenous melatonin profile in sighted humans: a brief review and critique of the literature. *J Biol Rhythms* 1997. [PMID 9406034](https://pubmed.ncbi.nlm.nih.gov/9406034/)
28. **Czeisler CA et al.** Stability, precision, and near-24-hour period of the human circadian pacemaker. *Science* 1999. [PMID 10381883](https://pubmed.ncbi.nlm.nih.gov/10381883/)
29. **Zeitzer JM et al.** Sensitivity of the human circadian pacemaker to nocturnal light: melatonin phase resetting and suppression. *J Physiol* 2000. [PMID 10922269](https://pubmed.ncbi.nlm.nih.gov/10922269/)
30. **Liu J et al.** MT1 and MT2 Melatonin Receptors: A Therapeutic Perspective. *Annu Rev Pharmacol Toxicol* 2016. [PMID 26514204](https://pubmed.ncbi.nlm.nih.gov/26514204/)
31. **Cajochen C et al.** Role of melatonin in the regulation of human circadian rhythms and sleep. *J Neuroendocrinol* 2003. [PMID 12622846](https://pubmed.ncbi.nlm.nih.gov/12622846/)
32. **Khalsa SB et al.** A phase response curve to single bright light pulses in human subjects. *J Physiol* 2003. [PMID 12717008](https://pubmed.ncbi.nlm.nih.gov/12717008/)
33. **St Hilaire MA et al.** Human phase response curve to a 1 h pulse of bright white light. *J Physiol* 2012. [PMID 22547633](https://pubmed.ncbi.nlm.nih.gov/22547633/)

## E. 각성계 · Flip-flop 스위치 · 오렉신

34. **Saper CB et al.** The sleep switch: hypothalamic control of sleep and wakefulness. *Trends Neurosci* 2001. [PMID 11718878](https://pubmed.ncbi.nlm.nih.gov/11718878/)
35. **Scammell TE et al.** Neural Circuitry of Wakefulness and Sleep. *Neuron* 2017. [PMID 28231463](https://pubmed.ncbi.nlm.nih.gov/28231463/)
36. **Sherin JE et al.** Activation of ventrolateral preoptic neurons during sleep. *Science* 1996. [PMID 8539624](https://pubmed.ncbi.nlm.nih.gov/8539624/)
37. **Sakurai T et al.** The neural circuit of orexin (hypocretin): maintaining sleep and wakefulness. *Nat Rev Neurosci* 2007. [PMID 17299454](https://pubmed.ncbi.nlm.nih.gov/17299454/)
38. **de Lecea L et al.** The hypocretins: hypothalamus-specific peptides with neuroexcitatory activity. *Proc Natl Acad Sci U S A* 1998. [PMID 9419374](https://pubmed.ncbi.nlm.nih.gov/9419374/)
39. **Mieda M et al.** Differential roles of orexin receptor-1 and -2 in the regulation of non-REM and REM sleep. *J Neurosci* 2011. [PMID 21525292](https://pubmed.ncbi.nlm.nih.gov/21525292/)
40. **Rudolph U et al.** Beyond classical benzodiazepines: novel therapeutic potential of GABAA receptor subtypes. *Nat Rev Drug Discov* 2011. [PMID 21799515](https://pubmed.ncbi.nlm.nih.gov/21799515/)
41. **Nutt DJ et al.** Searching for perfect sleep: the continuing evolution of GABAA receptor modulators as hypnotics. *J Psychopharmacol* 2010. [PMID 19942638](https://pubmed.ncbi.nlm.nih.gov/19942638/)
42. **Saper CB et al.** Hypothalamic regulation of sleep and circadian rhythms. *Nature* 2005. [PMID 16251950](https://pubmed.ncbi.nlm.nih.gov/16251950/)

## F. 과각성 가설 (핵심 병태생리)

43. **Riemann D et al.** The hyperarousal model of insomnia: a review of the concept and its evidence. *Sleep Med Rev* 2010. [PMID 19481481](https://pubmed.ncbi.nlm.nih.gov/19481481/)
44. **Bonnet MH et al.** Hyperarousal and insomnia: state of the science. *Sleep Med Rev* 2010. [PMID 19640748](https://pubmed.ncbi.nlm.nih.gov/19640748/)
45. **Nofzinger EA et al.** Functional neuroimaging evidence for hyperarousal in insomnia. *Am J Psychiatry* 2004. [PMID 15514418](https://pubmed.ncbi.nlm.nih.gov/15514418/)
46. **Perlis ML et al.** Beta/Gamma EEG activity in patients with primary and secondary insomnia and good sleeper controls. *Sleep* 2001. [PMID 11204046](https://pubmed.ncbi.nlm.nih.gov/11204046/)
47. **Vgontzas AN et al.** Chronic insomnia is associated with nyctohemeral activation of the hypothalamic-pituitary-adrenal axis: clinical implications. *J Clin Endocrinol Metab* 2001. [PMID 11502812](https://pubmed.ncbi.nlm.nih.gov/11502812/)
48. **Vgontzas AN et al.** Insomnia with objective short sleep duration: the most biologically severe phenotype of the disorder. *Sleep Med Rev* 2013. [PMID 23419741](https://pubmed.ncbi.nlm.nih.gov/23419741/)
49. **Riemann D et al.** "Hyperarousal and insomnia: state of the science". *Sleep Med Rev* 2010. [PMID 19945890](https://pubmed.ncbi.nlm.nih.gov/19945890/)
50. **Fernandez-Mendoza J et al.** Insomnia Phenotypes, Cardiovascular Risk and Their Link to Brain Health. *Circ Res* 2025. [PMID 40811499](https://pubmed.ncbi.nlm.nih.gov/40811499/)
51. **Colombo MA et al.** Wake High-Density Electroencephalographic Spatiospectral Signatures of Insomnia. *Sleep* 2016. [PMID 26951395](https://pubmed.ncbi.nlm.nih.gov/26951395/)
52. **Perlis ML et al.** Psychophysiological insomnia: the behavioural model and a neurocognitive perspective. *J Sleep Res* 1997. [PMID 9358396](https://pubmed.ncbi.nlm.nih.gov/9358396/)

## G. 3-P 모형 · 인지행동 이론

53. **Spielman AJ et al.** A behavioral perspective on insomnia treatment. *Psychiatr Clin North Am* 1987. [PMID 3332317](https://pubmed.ncbi.nlm.nih.gov/3332317/)
54. **Spielman AJ et al.** Treatment of chronic insomnia by restriction of time in bed. *Sleep* 1987. [PMID 3563247](https://pubmed.ncbi.nlm.nih.gov/3563247/)
55. **Harvey AG et al.** A cognitive model of insomnia. *Behav Res Ther* 2002. [PMID 12186352](https://pubmed.ncbi.nlm.nih.gov/12186352/)
56. **Espie CA et al.** The attention-intention-effort pathway in the development of psychophysiologic insomnia: a theoretical review. *Sleep Med Rev* 2006. [PMID 16809056](https://pubmed.ncbi.nlm.nih.gov/16809056/)
57. **Morin CM et al.** Dysfunctional beliefs and attitudes about sleep (DBAS): validation of a brief version (DBAS-16). *Sleep* 2007. [PMID 18041487](https://pubmed.ncbi.nlm.nih.gov/18041487/)

## H. CBT-I 임상시험 · 메타분석

58. **Trauer JM et al.** Cognitive Behavioral Therapy for Chronic Insomnia: A Systematic Review and Meta-analysis. *Ann Intern Med* 2015. [PMID 26054060](https://pubmed.ncbi.nlm.nih.gov/26054060/)
59. **van Straten A et al.** Cognitive and behavioral therapies in the treatment of insomnia: A meta-analysis. *Sleep Med Rev* 2018. [PMID 28392168](https://pubmed.ncbi.nlm.nih.gov/28392168/)
60. **Morin CM et al.** Cognitive behavioral therapy, singly and combined with medication, for persistent insomnia: a randomized controlled trial. *JAMA* 2009. [PMID 19454639](https://pubmed.ncbi.nlm.nih.gov/19454639/)
61. **Espie CA et al.** Effect of Digital Cognitive Behavioral Therapy for Insomnia on Health, Psychological Well-being, and Sleep-Related Quality of Life: A Randomized Clinical Trial. *JAMA Psychiatry* 2019. [PMID 30264137](https://pubmed.ncbi.nlm.nih.gov/30264137/)
62. **Miller CB et al.** The evidence base of sleep restriction therapy for treating insomnia disorder. *Sleep Med Rev* 2014. [PMID 24629826](https://pubmed.ncbi.nlm.nih.gov/24629826/)
63. **Kyle SD et al.** Sleep restriction therapy for insomnia is associated with reduced objective total sleep time, increased daytime somnolence, and objectively impaired vigilance: implications for the clinical management of insomnia disorder. *Sleep* 2014. [PMID 24497651](https://pubmed.ncbi.nlm.nih.gov/24497651/)
64. **Riemann D et al.** The treatments of chronic insomnia: a review of benzodiazepine receptor agonists and psychological and behavioral therapies. *Sleep Med Rev* 2009. [PMID 19201632](https://pubmed.ncbi.nlm.nih.gov/19201632/)

## I. BzRA (Z-drug) PK/PD · 안전성

65. **Greenblatt DJ et al.** Gender differences in pharmacokinetics and pharmacodynamics of zolpidem following sublingual administration. *J Clin Pharmacol* 2014. [PMID 24203450](https://pubmed.ncbi.nlm.nih.gov/24203450/)
66. **Krystal AD et al.** Sustained efficacy of eszopiclone over 6 months of nightly treatment: results of a randomized, double-blind, placebo-controlled study in adults with chronic insomnia. *Sleep* 2003. [PMID 14655910](https://pubmed.ncbi.nlm.nih.gov/14655910/)
67. **Walsh JK et al.** Nightly treatment of primary insomnia with eszopiclone for six months: effect on sleep, quality of life, and work limitations. *Sleep* 2007. [PMID 17702264](https://pubmed.ncbi.nlm.nih.gov/17702264/)
68. **Buscemi N et al.** The efficacy and safety of drug treatments for chronic insomnia in adults: a meta-analysis of RCTs. *J Gen Intern Med* 2007. [PMID 17619935](https://pubmed.ncbi.nlm.nih.gov/17619935/)
69. **Glass J et al.** Sedative hypnotics in older people with insomnia: meta-analysis of risks and benefits. *BMJ* 2005. [PMID 16284208](https://pubmed.ncbi.nlm.nih.gov/16284208/)
70. **Verster JC et al.** Residual effects of sleep medication on driving ability. *Sleep Med Rev* 2004. [PMID 15233958](https://pubmed.ncbi.nlm.nih.gov/15233958/)
71. **Soyka M et al.** Treatment of Benzodiazepine Dependence. *N Engl J Med* 2017. [PMID 28328330](https://pubmed.ncbi.nlm.nih.gov/28328330/)
72. **Roth T et al.** Issues in the use of benzodiazepine therapy. *J Clin Psychiatry* 1992. [PMID 1613014](https://pubmed.ncbi.nlm.nih.gov/1613014/)
73. **Drover DR et al.** Comparative pharmacokinetics and pharmacodynamics of short-acting hypnosedatives: zaleplon, zolpidem and zopiclone. *Clin Pharmacokinet* 2004. [PMID 15005637](https://pubmed.ncbi.nlm.nih.gov/15005637/)

## J. DORA (오렉신 길항제)

74. **Herring WJ et al.** Suvorexant in Patients With Insomnia: Results From Two 3-Month Randomized Controlled Clinical Trials. *Biol Psychiatry* 2016. [PMID 25526970](https://pubmed.ncbi.nlm.nih.gov/25526970/)
75. **Michelson D et al.** Safety and efficacy of suvorexant during 1-year treatment of insomnia with subsequent abrupt treatment discontinuation: a phase 3 randomised, double-blind, placebo-controlled trial. *Lancet Neurol* 2014. [PMID 24680372](https://pubmed.ncbi.nlm.nih.gov/24680372/)
76. **Rosenberg R et al.** Comparison of Lemborexant With Placebo and Zolpidem Tartrate Extended Release for the Treatment of Older Adults With Insomnia Disorder: A Phase 3 Randomized Clinical Trial. *JAMA Netw Open* 2019. [PMID 31880796](https://pubmed.ncbi.nlm.nih.gov/31880796/)
77. **Kärppä M et al.** Long-term efficacy and tolerability of lemborexant compared with placebo in adults with insomnia disorder: results from the phase 3 randomized clinical trial SUNRISE 2. *Sleep* 2020. [PMID 32585700](https://pubmed.ncbi.nlm.nih.gov/32585700/)
78. **Mignot E et al.** Safety and efficacy of daridorexant in patients with insomnia disorder: results from two multicentre, randomised, double-blind, placebo-controlled, phase 3 trials. *Lancet Neurol* 2022. [PMID 35065036](https://pubmed.ncbi.nlm.nih.gov/35065036/)
79. **Muehlan C et al.** Accelerated Development of the Dual Orexin Receptor Antagonist ACT-541468: Integration of a Microtracer in a First-in-Human Study. *Clin Pharmacol Ther* 2018. [PMID 29446069](https://pubmed.ncbi.nlm.nih.gov/29446069/)
80. **Ufer M et al.** Abuse potential assessment of the new dual orexin receptor antagonist daridorexant in recreational sedative drug users as compared to suvorexant and zolpidem. *Sleep* 2022. [PMID 34480579](https://pubmed.ncbi.nlm.nih.gov/34480579/)
81. **Kishi T et al.** Comparative efficacy and safety of daridorexant, lemborexant, and suvorexant for insomnia: a systematic review and network meta-analysis. *Transl Psychiatry* 2025. [PMID 40555730](https://pubmed.ncbi.nlm.nih.gov/40555730/)
82. **Dauvilliers Y et al.** Daridorexant, a New Dual Orexin Receptor Antagonist to Treat Insomnia Disorder. *Ann Neurol* 2020. [PMID 31953863](https://pubmed.ncbi.nlm.nih.gov/31953863/)

## K. 멜라토닌계 · H1 · 5-HT2A 약물

83. **Kuriyama A et al.** Ramelteon for the treatment of insomnia in adults: a systematic review and meta-analysis. *Sleep Med* 2014. [PMID 24656909](https://pubmed.ncbi.nlm.nih.gov/24656909/)
84. **Karim A et al.** Disposition kinetics and tolerance of escalating single doses of ramelteon, a high-affinity MT1 and MT2 melatonin receptor agonist indicated for treatment of insomnia. *J Clin Pharmacol* 2006. [PMID 16432265](https://pubmed.ncbi.nlm.nih.gov/16432265/)
85. **Ferracioli-Oda E et al.** Meta-analysis: melatonin for the treatment of primary sleep disorders. *PLoS One* 2013. [PMID 23691095](https://pubmed.ncbi.nlm.nih.gov/23691095/)
86. **Wade AG et al.** Prolonged release melatonin in the treatment of primary insomnia: evaluation of the age cut-off for short- and long-term response. *Curr Med Res Opin* 2011. [PMID 21091391](https://pubmed.ncbi.nlm.nih.gov/21091391/)
87. **Stahl SM et al.** Mechanism of action of trazodone: a multifunctional drug. *CNS Spectr* 2009. [PMID 20095366](https://pubmed.ncbi.nlm.nih.gov/20095366/)

## L. 카페인 · 알코올 · 광 · 교대근무

88. **Drake C et al.** Caffeine effects on sleep taken 0, 3, or 6 hours before going to bed. *J Clin Sleep Med* 2013. [PMID 24235903](https://pubmed.ncbi.nlm.nih.gov/24235903/)
89. **Clark I et al.** Coffee, caffeine, and sleep: A systematic review of epidemiological studies and randomized controlled trials. *Sleep Med Rev* 2017. [PMID 26899133](https://pubmed.ncbi.nlm.nih.gov/26899133/)
90. **Roehrs T et al.** Sleep, sleepiness, and alcohol use. *Alcohol Res Health* 2001. [PMID 11584549](https://pubmed.ncbi.nlm.nih.gov/11584549/)
91. **van Maanen A et al.** The effects of light therapy on sleep problems: A systematic review and meta-analysis. *Sleep Med Rev* 2016. [PMID 26606319](https://pubmed.ncbi.nlm.nih.gov/26606319/)
92. **Chang AM et al.** Evening use of light-emitting eReaders negatively affects sleep, circadian timing, and next-morning alertness. *Proc Natl Acad Sci U S A* 2015. [PMID 25535358](https://pubmed.ncbi.nlm.nih.gov/25535358/)
93. **Wright KP Jr et al.** Shift work and the assessment and management of shift work disorder (SWD). *Sleep Med Rev* 2013. [PMID 22560640](https://pubmed.ncbi.nlm.nih.gov/22560640/)
94. **Ebrahim IO et al.** Alcohol and sleep I: effects on normal sleep. *Alcohol Clin Exp Res* 2013. [PMID 23347102](https://pubmed.ncbi.nlm.nih.gov/23347102/)
95. **Auger RR et al.** Clinical Practice Guideline for the Treatment of Intrinsic Circadian Rhythm Sleep-Wake Disorders: Advanced Sleep-Wake Phase Disorder (ASWPD), Delayed Sleep-Wake Phase Disorder (DSWPD), Non-24-Hour Sleep-Wake Rhythm Disorder (N24SWD), and Irregular Sleep-Wake Rhythm Disorder (ISWRD). An Update for 2015: An American Academy of Sleep Medicine Clinical Practice Guideline. *J Clin Sleep Med* 2015. [PMID 26414986](https://pubmed.ncbi.nlm.nih.gov/26414986/)

## M. 결과 · 동반질환 · 정량 모델링

96. **Baglioni C et al.** Insomnia as a predictor of depression: a meta-analytic evaluation of longitudinal epidemiological studies. *J Affect Disord* 2011. [PMID 21300408](https://pubmed.ncbi.nlm.nih.gov/21300408/)
97. **Sofi F et al.** Insomnia and risk of cardiovascular disease: a meta-analysis. *Eur J Prev Cardiol* 2014. [PMID 22942213](https://pubmed.ncbi.nlm.nih.gov/22942213/)
98. **Anothaisintawee T et al.** Sleep disturbances compared to traditional risk factors for diabetes development: Systematic review and meta-analysis. *Sleep Med Rev* 2016. [PMID 26687279](https://pubmed.ncbi.nlm.nih.gov/26687279/)
99. **Irwin MR et al.** Sleep Disturbance, Sleep Duration, and Inflammation: A Systematic Review and Meta-Analysis of Cohort Studies and Experimental Sleep Deprivation. *Biol Psychiatry* 2016. [PMID 26140821](https://pubmed.ncbi.nlm.nih.gov/26140821/)
100. **Sweetman A et al.** Co-Morbid Insomnia and Sleep Apnea (COMISA): Prevalence, Consequences, Methodological Considerations, and Recent Randomized Controlled Trials. *Brain Sci* 2019. [PMID 31842520](https://pubmed.ncbi.nlm.nih.gov/31842520/)
101. **Bizzotto R et al.** Multinomial logistic estimation of Markov-chain models for modeling sleep architecture in primary insomnia patients. *J Pharmacokinet Pharmacodyn* 2010. [PMID 20052524](https://pubmed.ncbi.nlm.nih.gov/20052524/)
102. **Ohayon MM et al.** Meta-analysis of quantitative sleep parameters from childhood to old age in healthy individuals: developing normative sleep values across the human lifespan. *Sleep* 2004. [PMID 15586779](https://pubmed.ncbi.nlm.nih.gov/15586779/)
103. **Bizzotto R et al.** Multinomial logistic functions in markov chain models of sleep architecture: internal and external validation and covariate analysis. *AAPS J* 2011. [PMID 21691915](https://pubmed.ncbi.nlm.nih.gov/21691915/)
104. **Postnova S et al.** Prediction of Cognitive Performance and Subjective Sleepiness Using a Model of Arousal Dynamics. *J Biol Rhythms* 2018. [PMID 29671707](https://pubmed.ncbi.nlm.nih.gov/29671707/)
105. **Phillips AJ et al.** Probing the mechanisms of chronotype using quantitative modeling. *J Biol Rhythms* 2010. [PMID 20484693](https://pubmed.ncbi.nlm.nih.gov/20484693/)
---

## 문헌과 모델 파라미터의 대응 (Mapping to model parameters)

| 모델 구성요소 | 파라미터 | 주요 근거 |
|---|---|---|
| Process S 상승·감쇠 시상수 | `TAUR` 18.2 h · `TAUD` 4.2 h | Daan/Beersma/Borbély 1984 · Borbély 2016 |
| 일주기 각성 신호 진폭·위상 | `AMPC` · `TPKC` | Dijk & Czeisler 1995 |
| 내인성 주기 | `TAUINT` 24.18 h | Czeisler 1999 |
| 광 위상반응곡선 | `KLIGHT` · `prcL` | Khalsa/St Hilaire PRC 연구 |
| 멜라토닌 위상반응곡선 | `KMELPH` · `prcM` | Lewy melatonin PRC |
| 야간 멜라토닌 광억제 | `LUX50` 90 lux | Zeitzer 2000 |
| 아데노신·카페인 경쟁 | `KICAF` 5 mg/L | Porkka-Heiskanen 1997 · Landolt |
| Flip-flop 상호억제 | `WNE` · `WHA` · `WOX` | Saper 2005 · Scammell 2017 |
| 과각성 상태변수 | `A0` · `AROU` | Riemann 2010 · Bonnet & Arand 2010 |
| HPA 야간 상승 | `GCORT` · `AMPCOR` | Vgontzas 2001 |
| 3-P 유지 고리 | `COND` · `SEFF` · `DBAS` · `TIBS` | Spielman 1987 · Harvey 2002 · Espie 2006 |
| 수면 제한 요법 | `KTIBDN` · `SETGT` | Spielman 1987 · Miller 2014 |
| 졸피뎀 PK | `KAZOL` · `VZOL` · `CLZOL` | zolpidem PK 문헌 · Greenblatt |
| DORA 반감기 대비 | `CLSUV` · `CLLEM` · `CLDAR` | Herring 2016 · Rosenberg 2019 · Mignot 2022 |
| 라멜테온 + M-II | `FMM2` · `CLRM2` | Karim ramelteon 용량증가 연구 |
| 저용량 독세핀 H1 | `EC50DOX` · `KADOX` | Krystal doxepin 35일 수면검사실 시험 |
| BzRA 내성·반동 | `KTOLIN` · `KTOLOUT` · `GREB` | Roehrs & Roth · Soyka |
| 글림프 청소 | `KGIN` · `N3REF` | Xie 2013 계열 문헌 |
| 불면–우울 양방향 | `GDEPISI` | Baglioni 2011 |
| 마르코프 수면구조 모델링 참조 | 전체 아키텍처 | Bizzotto · Karlsson 수면 마르코프 모델 |
