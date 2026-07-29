# Idiopathic Intracranial Hypertension — annotated references

Every entry below was retrieved from PubMed programmatically (NCBI
E-utilities) while this model was being built, so the author, year,
journal, volume, pages and PMID are the record's own and not transcribed
from memory. The script that produced this file is
`build_refs.py`-style: it takes a list of PMIDs plus the annotation
explaining what each one is used FOR in the model, and fills in the
citation from the database.

**72 references**, grouped by the role they play in the model.
Entries marked **[CALIBRATION]** supply a number that a parameter was fitted
to; **[TEST]** entries supply a number used to test the model without being
fitted to; **[STRUCTURE]** entries justify a piece of the equations.

---

## 1. Definition, diagnosis, epidemiology and natural history

1. Friedman DI et al. (2013). *Revised diagnostic criteria for the pseudotumor cerebri syndrome in adults and children.* Neurology 81:1159-65. [PMID 23966248](https://pubmed.ncbi.nlm.nih.gov/23966248/)
   The revised diagnostic criteria used throughout this model, including the >25 cmH2O (250 mmH2O) opening-pressure threshold that defines 'remission' in every scenario here.

2. Mollan SP et al. (2018). *Idiopathic intracranial hypertension: consensus guidelines on management.* J Neurol Neurosurg Psychiatry 89:1088-1100. [PMID 29903905](https://pubmed.ncbi.nlm.nih.gov/29903905/)
   Consensus management guidance; the source for the therapy ladder encoded in the scenario list.

3. Mollan SP et al. (2018). *Evaluation and management of adult idiopathic intracranial hypertension.* Pract Neurol 18:485-488. [PMID 30154235](https://pubmed.ncbi.nlm.nih.gov/30154235/)
   Practical companion to the consensus guidance.

4. Binder DK et al. (2004). *Idiopathic intracranial hypertension.* Neurosurgery 54:538-51; discussion 551-2. [PMID 15028127](https://pubmed.ncbi.nlm.nih.gov/15028127/)
   Older but still useful synthesis of mechanism and surgical management.

5. Thambisetty M et al. (2007). *Fulminant idiopathic intracranial hypertension.* Neurology 68:229-32. [PMID 17224579](https://pubmed.ncbi.nlm.nih.gov/17224579/)
   Fulminant IIH: the phenotype the model reproduces as a high venous floor with days-to-weeks to blindness, in which dose escalation cannot work.

6. Virdee J et al. (2020). *Reviewing the Recent Developments in Idiopathic Intracranial Hypertension.* Ophthalmol Ther 9:767-781. [PMID 32902722](https://pubmed.ncbi.nlm.nih.gov/32902722/)
   Review of developments up to 2020, including the shift towards venous and metabolic mechanisms.

7. Adderley NJ et al. (2019). *Association Between Idiopathic Intracranial Hypertension and Risk of Cardiovascular Diseases in Women in the United Kingdom.* JAMA Neurol 76:1088-1098. [PMID 31282950](https://pubmed.ncbi.nlm.nih.gov/31282950/)
   Long-term cardiovascular risk in women with IIH — the reason the weight axis in this model is not only a pressure lever.

8. Botrous M et al. (2026). *Incidence of idiopathic intracranial hypertension (IIH) in adults and children: a 14-year population-based study.* J Neurol 273. [PMID 42287435](https://pubmed.ncbi.nlm.nih.gov/42287435/)
   Recent population-based incidence data; IIH incidence tracks obesity prevalence.

9. Wall M (2026). *What Have We Learned From the Idiopathic Intracranial Hypertension Treatment Trial the William F. Hoyt Lecture.* J Neuroophthalmol 46:271-278. [PMID 42133960](https://pubmed.ncbi.nlm.nih.gov/42133960/)
   Retrospective on what the IIH Treatment Trial did and did not establish.


## 2. CSF secretion — the I_f term

10. Damkier HH et al. (2013). *Cerebrospinal fluid secretion by the choroid plexus.* Physiol Rev 93:1847-92. [PMID 24137023](https://pubmed.ncbi.nlm.nih.gov/24137023/) **[CALIBRATION]** **[STRUCTURE]**
   The definitive review of choroid-plexus CSF secretion: Na+/K+-ATPase, NKCC1, carbonic anhydrase, AQP1. Source of the 0.35 mL/min formation rate used as IF0.

11. Sinclair AJ et al. (2007). *Corticosteroids, 11beta-hydroxysteroid dehydrogenase isozymes and the rabbit choroid plexus.* J Neuroendocrinol 19:614-20. [PMID 17620103](https://pubmed.ncbi.nlm.nih.gov/17620103/)
   11beta-HSD isozymes in choroid plexus — the mechanistic basis for the glucocorticoid/androgen input to secretion.

12. Sinclair AJ et al. (2010). *Cerebrospinal fluid corticosteroid levels and cortisol metabolism in patients with idiopathic intracranial hypertension: a link between 11beta-HSD1 and intracranial pressure regulation?.* J Clin Endocrinol Metab 95:5348-56. [PMID 20826586](https://pubmed.ncbi.nlm.nih.gov/20826586/)
   CSF corticosteroid levels and cortisol metabolism in IIH; links 11beta-HSD1 to intracranial pressure regulation.

13. O'Reilly MW et al. (2019). *A unique androgen excess signature in idiopathic intracranial hypertension is linked to cerebrospinal fluid dynamics.* JCI Insight 4. [PMID 30753168](https://pubmed.ncbi.nlm.nih.gov/30753168/)
   The androgen-excess signature specific to IIH, linked to CSF dynamics — the anchor for this model's androgen->secretion channel.

14. Libien J et al. (2017). *Role of vitamin A metabolism in IIH: Results from the idiopathic intracranial hypertension treatment trial.* J Neurol Sci 372:78-84. [PMID 28017254](https://pubmed.ncbi.nlm.nih.gov/28017254/)
   Vitamin A/retinol metabolism in the IIHTT cohort; an alternative secretion/absorption input, included in the map but not quantified in the ODEs.


## 3. CSF outflow — the R_out term

15. Børgesen SE et al. (1992). *Computerized infusion test compared to steady pressure constant infusion test in measurement of resistance to CSF outflow.* Acta Neurochir (Wien) 119:12-6. [PMID 1481738](https://pubmed.ncbi.nlm.nih.gov/1481738/) **[STRUCTURE]**
   Infusion-test methodology for measuring resistance to CSF outflow — the measurement this model argues is biased by the venous loop.

16. Aspelund A et al. (2015). *A dural lymphatic vascular system that drains brain interstitial fluid and macromolecules.* J Exp Med 212:991-9. [PMID 26077718](https://pubmed.ncbi.nlm.nih.gov/26077718/) **[STRUCTURE]**
   Discovery of a dural lymphatic vascular system draining brain interstitial fluid.

17. Louveau A et al. (2015). *Structural and functional features of central nervous system lymphatic vessels.* Nature 523:337-41. [PMID 26030524](https://pubmed.ncbi.nlm.nih.gov/26030524/) **[STRUCTURE]**
   Structural and functional characterisation of CNS lymphatic vessels.

18. Iliff JJ et al. (2012). *A paravascular pathway facilitates CSF flow through the brain parenchyma and the clearance of interstitial solutes, including amyloid β.* Sci Transl Med 4:147ra111. [PMID 22896675](https://pubmed.ncbi.nlm.nih.gov/22896675/) **[STRUCTURE]**
   The glymphatic (paravascular) clearance pathway.

19. Da Mesquita S et al. (2018). *Functional aspects of meningeal lymphatics in ageing and Alzheimer's disease.* Nature 560:185-191. [PMID 30046111](https://pubmed.ncbi.nlm.nih.gov/30046111/) **[STRUCTURE]**
   Meningeal lymphatic function declines with ageing — the motivation for the slow outflow-reserve loop, which has no human IIH anchor and is flagged as such.


## 4. The collapsible sinus — the venous floor and the loop gain

20. Farb RI et al. (2003). *Idiopathic intracranial hypertension: the prevalence and morphology of sinovenous stenosis.* Neurology 60:1418-24. [PMID 12743224](https://pubmed.ncbi.nlm.nih.gov/12743224/) **[STRUCTURE]**
   Bilateral transverse sinus stenosis is present in almost all IIH and almost no controls — the anatomical basis of the amplifier.

21. King JO et al. (1995). *Cerebral venography and manometry in idiopathic intracranial hypertension.* Neurology 45:2224-8. [PMID 8848197](https://pubmed.ncbi.nlm.nih.gov/8848197/) **[STRUCTURE]**
   Cerebral venography with manometry in IIH: the original demonstration of elevated sinus pressures and trans-stenotic gradients.

22. King JO et al. (2002). *Manometry combined with cervical puncture in idiopathic intracranial hypertension.* Neurology 58:26-30. [PMID 11781401](https://pubmed.ncbi.nlm.nih.gov/11781401/) **[STRUCTURE]**
   Manometry combined with cervical puncture: sinus pressure falls when CSF pressure is lowered, i.e. the stenosis is a consequence as well as a cause.

23. Rohr A et al. (2007). *Reversibility of venous sinus obstruction in idiopathic intracranial hypertension.* AJNR Am J Neuroradiol 28:656-9. [PMID 17416816](https://pubmed.ncbi.nlm.nih.gov/17416816/) **[STRUCTURE]**
   Venous sinus obstruction in IIH is REVERSIBLE after CSF pressure is lowered — direct evidence the sinus behaves as a collapsible tube.

24. Lalou AD et al. (2020). *Coupling of CSF and sagittal sinus pressure in adult patients with pseudotumour cerebri.* Acta Neurochir (Wien) 162:1001-1009. [PMID 31832847](https://pubmed.ncbi.nlm.nih.gov/31832847/) **[STRUCTURE]**
   The key quantitative anchor for the loop gain: CSF and sagittal sinus pressure are coupled (R=0.96 baseline, R=0.92 during infusion), and during CSF drainage they track until sinus pressure bottoms out while CSF pressure keeps falling. The paper explicitly discusses the implications for calculating outflow resistance.

25. Bateman GA et al. (2009). *A mathematical model of idiopathic intracranial hypertension incorporating increased arterial inflow and variable venous outflow collapsibility.* J Neurosurg 110:446-56. [PMID 18847344](https://pubmed.ncbi.nlm.nih.gov/18847344/) **[STRUCTURE]**
   A mathematical model of IIH incorporating variable venous outflow collapsibility — the closest published antecedent of the structure used here.

26. Boddu SR et al. (2018). *Pressure variations in cerebral venous sinuses of idiopathic intracranial hypertension patients.* J Vasc Interv Neurol 10:25-30. [PMID 29922401](https://pubmed.ncbi.nlm.nih.gov/29922401/) **[CALIBRATION]**
   Sinus manometry before and after stenting: sagittal sinus pressure -8.1 mmHg, trans-stenotic gradient -15.7 mmHg. Numerator of this model's loop-gain calculation.


## 5. Obesity, intra-abdominal pressure and the venous base P_cv

27. Sugerman HJ et al. (1997). *Increased intra-abdominal pressure and cardiac filling pressures in obesity-associated pseudotumor cerebri.* Neurology 49:507-11. [PMID 9270586](https://pubmed.ncbi.nlm.nih.gov/9270586/) **[CALIBRATION]**
   Direct measurement in obesity-associated pseudotumour cerebri: bladder pressure 22+-3 cmH2O, central venous pressure 20+-6 mmHg, pleural pressure 15+-10 mmHg. The anchor for the abdomen->thorax->venous chain.

28. Sugerman HJ (2001). *Effects of increased intra-abdominal pressure in severe obesity.* Surg Clin North Am 81:1063-75, vi. [PMID 11589245](https://pubmed.ncbi.nlm.nih.gov/11589245/)
   Broader review of the haemodynamic consequences of raised intra-abdominal pressure in severe obesity.

29. Sugerman HJ et al. (2001). *Continuous negative abdominal pressure device to treat pseudotumor cerebri.* Int J Obes Relat Metab Disord 25:486-90. [PMID 11319651](https://pubmed.ncbi.nlm.nih.gov/11319651/) **[TEST]**
   A continuous negative abdominal pressure device lowered bladder pressure from 19.1 to 12.5 cmH2O and relieved headache and pulsatile tinnitus within FIVE MINUTES — a near-experimental confirmation of the abdominal-venous pathway and of its short time constant.

30. Berdahl JP et al. (2012). *Body mass index has a linear relationship with cerebrospinal fluid pressure.* Invest Ophthalmol Vis Sci 53:1422-7. [PMID 22323469](https://pubmed.ncbi.nlm.nih.gov/22323469/) **[TEST]**
   Body mass index has a linear relationship with CSF pressure — the population-level version of this model's weight->P_cv->ICP chain.


## 6. Weight-loss interventions — the pressure/weight slope

31. Sinclair AJ et al. (2010). *Low energy diet and intracranial pressure in women with idiopathic intracranial hypertension: prospective cohort study.* BMJ 341:c2701. [PMID 20610512](https://pubmed.ncbi.nlm.nih.gov/20610512/) **[TEST]**
   Very-low-calorie diet cohort acting as its own control: -15.7 kg, ICP -8.0 cmH2O, HIT-6 -7.6. Implies 0.51 cmH2O/kg.

32. Mollan SP et al. (2021). *Effectiveness of Bariatric Surgery vs Community Weight Management Intervention for the Treatment of Idiopathic Intracranial Hypertension: A Randomized Clinical Trial.* JAMA Neurol 78:678-686. [PMID 33900360](https://pubmed.ncbi.nlm.nih.gov/33900360/) **[CALIBRATION]**
   IIH:WT randomised trial of bariatric surgery vs community weight management: ICP -6.0 cmCSF at 12 months and -8.2 at 24 months for -21.4 and -26.6 kg. Implies 0.28-0.31 cmH2O/kg — the slope this model uses.

33. Mollan SP et al. (2022). *Association of Amount of Weight Lost After Bariatric Surgery With Intracranial Pressure in Women With Idiopathic Intracranial Hypertension.* Neurology 99:e1090-e1099. [PMID 35790425](https://pubmed.ncbi.nlm.nih.gov/35790425/) **[TEST]**
   IIH:WT substudy relating the AMOUNT of weight lost to ICP (R2 = 0.47), with ~24% weight loss associated with remission. Used as an out-of-sample test, not a calibration target.

34. Abbott S et al. (2023). *Weight Management Interventions for Adults With Idiopathic Intracranial Hypertension: A Systematic Review and Practice Recommendations.* Neurology 101:e2138-e2150. [PMID 37813577](https://pubmed.ncbi.nlm.nih.gov/37813577/)
   Systematic review of weight-management interventions and practice recommendations in IIH.

35. Merola J et al. (2020). *Cerebrospinal fluid diversion versus bariatric surgery in the management of idiopathic intracranial hypertension.* Br J Neurosurg 34:9-12. [PMID 31805794](https://pubmed.ncbi.nlm.nih.gov/31805794/)
   CSF diversion versus bariatric surgery — the comparison the venous-floor argument reframes.


## 7. Acetazolamide and other secretion-blocking drugs

36. NORDIC Idiopathic Intracranial Hypertension Study Group Writing Committee et al. (2014). *Effect of acetazolamide on visual function in patients with idiopathic intracranial hypertension and mild visual loss: the idiopathic intracranial hypertension treatment trial.* JAMA 311:1641-51. [PMID 24756514](https://pubmed.ncbi.nlm.nih.gov/24756514/) **[CALIBRATION]**
   The IIH Treatment Trial. Primary calibration source: opening pressure 357.2 -> 244.9 (acetazolamide) vs 304.8 mmH2O (placebo), between-arm difference -59.9 mmH2O; PMD +1.43 vs +0.71 dB; Frisen 2.76 -> 1.45 vs 2.15; weight -7.50 vs -3.45 kg; mean achieved dose 2.5 g/day against a 4 g/day target; and the mediation analysis putting the weight-mediated part of the visual benefit at 0.03 dB.

37. ten Hove MW et al. (2016). *Safety and Tolerability of Acetazolamide in the Idiopathic Intracranial Hypertension Treatment Trial.* J Neuroophthalmol 36:13-9. [PMID 26587993](https://pubmed.ncbi.nlm.nih.gov/26587993/) **[TEST]**
   Safety and tolerability of acetazolamide in IIHTT — the empirical basis for modelling adherence, not potency, as the binding constraint on dose.

38. Celebisoy N et al. (2007). *Treatment of idiopathic intracranial hypertension: topiramate vs acetazolamide, an open-label study.* Acta Neurol Scand 116:322-7. [PMID 17922725](https://pubmed.ncbi.nlm.nih.gov/17922725/) **[TEST]**
   Topiramate versus acetazolamide, open-label: the basis for topiramate's weak carbonic-anhydrase effect plus a weight-loss component.

39. Mitchell JL et al. (2023). *The effect of GLP-1RA exenatide on idiopathic intracranial hypertension: a randomized clinical trial.* Brain 146:1821-1830. [PMID 36907221](https://pubmed.ncbi.nlm.nih.gov/36907221/) **[CALIBRATION]**
   Exenatide randomised trial: ICP -5.7 cmH2O at 2.5 hours, -6.4 at 24 hours, -5.6 at 12 weeks, from a baseline of 30.6 cmCSF. The 2.5-hour value calibrates EC50_GLP and forces the GLP-1R secretory ceiling close to acetazolamide's.

40. Botfield HF et al. (2017). *A glucagon-like peptide-1 receptor agonist reduces intracranial pressure in a rat model of hydrocephalus.* Sci Transl Med 9. [PMID 28835515](https://pubmed.ncbi.nlm.nih.gov/28835515/)
   Mechanistic origin of the GLP-1 hypothesis: a GLP-1 receptor agonist lowers ICP in rats, with choroid-plexus Na+/K+-ATPase as the target.

41. Grech O et al. (2024). *Effect of glucagon like peptide-1 receptor agonist exenatide, used as an intracranial pressure lowering agent, on cognition in Idiopathic Intracranial Hypertension.* Eye (Lond) 38:1374-1379. [PMID 38212401](https://pubmed.ncbi.nlm.nih.gov/38212401/)
   Cognitive outcomes alongside the exenatide ICP result.

42. Markey K et al. (2020). *11β-Hydroxysteroid dehydrogenase type 1 inhibition in idiopathic intracranial hypertension: a double-blind randomized controlled trial.* Brain Commun 2:fcz050. [PMID 32954315](https://pubmed.ncbi.nlm.nih.gov/32954315/)
   AZD4017, an 11beta-HSD1 inhibitor, in a double-blind randomised trial — the negative result that keeps the androgen/glucocorticoid channel modest in this model.

43. Markey KA et al. (2017). *Assessing the Efficacy and Safety of an 11β-Hydroxysteroid Dehydrogenase Type 1 Inhibitor (AZD4017) in the Idiopathic Intracranial Hypertension Drug Trial, IIH:DT: Clinical Methods and Design for a Phase II Randomized Controlled Trial.* JMIR Res Protoc 6:e181. [PMID 28923789](https://pubmed.ncbi.nlm.nih.gov/28923789/)
   Design and methods of that trial (IIH:DT).

44. Hardy RS et al. (2021). *11βHSD1 Inhibition with AZD4017 Improves Lipid Profiles and Lean Muscle Mass in Idiopathic Intracranial Hypertension.* J Clin Endocrinol Metab 106:174-187. [PMID 33098644](https://pubmed.ncbi.nlm.nih.gov/33098644/)
   Metabolic effects of AZD4017 in the same cohort.

45. Stefanou MI et al. (2025). *Efficacy and Safety of GLP-1 and Dual GIP/GLP-1 Receptor Agonists in Idiopathic Intracranial Hypertension: A Systematic Review and Meta-Analysis.* Eur J Neurol 32:e70358. [PMID 40937960](https://pubmed.ncbi.nlm.nih.gov/40937960/)
   Meta-analysis of GLP-1 and dual GIP/GLP-1 receptor agonists in IIH.

46. Ognard J et al. (2025). *Use of glucagon-like peptide-1 receptor agonists in idiopathic intracranial hypertension : a systematic review.* J Headache Pain 26:202. [PMID 41057780](https://pubmed.ncbi.nlm.nih.gov/41057780/)
   Systematic review of GLP-1 receptor agonist use in IIH.

47. de Oliveira HM et al. (2026). *Effect of GLP-1 receptor agonists on idiopathic intracranial hypertension: A systematic review.* Headache 66:286-297. [PMID 41246926](https://pubmed.ncbi.nlm.nih.gov/41246926/)
   A second, independent systematic review of the same question.


## 8. Venous sinus stenting — moving the floor

48. Higgins JN et al. (2002). *Venous sinus stenting for refractory benign intracranial hypertension.* Lancet 359:228-30. [PMID 11812561](https://pubmed.ncbi.nlm.nih.gov/11812561/)
   The first report of venous sinus stenting for refractory raised CSF pressure.

49. Patsalides A et al. (2019). *Venous sinus stenting lowers the intracranial pressure in patients with idiopathic intracranial hypertension.* J Neurointerv Surg 11:175-178. [PMID 29871989](https://pubmed.ncbi.nlm.nih.gov/29871989/) **[CALIBRATION]**
   The quantitative anchor for stenting: CSF opening pressure 37.0 -> 20.2 cmH2O at 3 months in 50 patients, WITH acetazolamide dose falling from 950 to 300 mg/day and weight RISING 1.1 kg — so the pressure fall cannot be attributed to drug or weight. Denominator of the loop-gain calculation.

50. Dinkin MJ & Patsalides A (2017). *Venous Sinus Stenting in Idiopathic Intracranial Hypertension: Results of a Prospective Trial.* J Neuroophthalmol 37:113-121. [PMID 27556959](https://pubmed.ncbi.nlm.nih.gov/27556959/)
   Prospective stenting trial with visual and headache outcomes.

51. Nicholson P et al. (2019). *Venous sinus stenting for idiopathic intracranial hypertension: a systematic review and meta-analysis.* J Neurointerv Surg 11:380-385. [PMID 30166333](https://pubmed.ncbi.nlm.nih.gov/30166333/)
   Systematic review and meta-analysis of stenting outcomes.

52. Kahan J et al. (2021). *Predicting the need for retreatment in venous sinus stenting for idiopathic intracranial hypertension.* J Neurointerv Surg 13:574-579. [PMID 32895320](https://pubmed.ncbi.nlm.nih.gov/32895320/)
   Predictors of the need for retreatment — the basis for the restenosis term.

53. Boddu S et al. (2016). *Resolution of Pulsatile Tinnitus after Venous Sinus Stenting in Patients with Idiopathic Intracranial Hypertension.* PLoS One 11:e0164466. [PMID 27768690](https://pubmed.ncbi.nlm.nih.gov/27768690/) **[TEST]**
   Resolution of pulsatile tinnitus after stenting, which is why tinnitus in this model tracks the trans-stenotic gradient rather than ICP.

54. Gurney SP et al. (2020). *Exploring The Current Management Idiopathic Intracranial Hypertension, And Understanding The Role Of Dural Venous Sinus Stenting.* Eye Brain 12:1-13. [PMID 32021528](https://pubmed.ncbi.nlm.nih.gov/32021528/)
   Review of current management with specific attention to where stenting fits.


## 9. Optic nerve, papilloedema and the two OCT channels

55. Hayreh SS (2016). *Pathogenesis of optic disc edema in raised intracranial pressure.* Prog Retin Eye Res 50:108-44. [PMID 26453995](https://pubmed.ncbi.nlm.nih.gov/26453995/) **[STRUCTURE]**
   Comprehensive account of the pathogenesis of optic disc oedema in raised intracranial pressure: axoplasmic stasis at the lamina cribrosa driven by the translaminar pressure gradient.

56. Tso MO & Hayreh SS (1977). *Optic disc edema in raised intracranial pressure. IV. Axoplasmic transport in experimental papilledema.* Arch Ophthalmol 95:1458-62. [PMID 70201](https://pubmed.ncbi.nlm.nih.gov/70201/) **[STRUCTURE]**
   Experimental demonstration that papilloedema is a disorder of axoplasmic transport.

57. Jóhannesson G et al. (2018). *Intracranial and Intraocular Pressure at the Lamina Cribrosa: Gradient Effects.* Curr Neurol Neurosci Rep 18:25. [PMID 29651628](https://pubmed.ncbi.nlm.nih.gov/29651628/) **[STRUCTURE]**
   Intracranial and intraocular pressure at the lamina cribrosa: the gradient this model uses as the actual insult.

58. Liu KC et al. (2020). *Current concepts of cerebrospinal fluid dynamics and the translaminar cribrosa pressure gradient: a paradigm of optic disk disease.* Surv Ophthalmol 65:48-66. [PMID 31449832](https://pubmed.ncbi.nlm.nih.gov/31449832/) **[STRUCTURE]**
   CSF dynamics and the translaminar pressure gradient across optic disc disease.

59. Frisén L (2017). *Swelling of the Optic Nerve Head: A Backstage View of a Staging Scheme.* J Neuroophthalmol 37:3-6. [PMID 28187078](https://pubmed.ncbi.nlm.nih.gov/28187078/)
   Frisen's own account of his papilloedema staging scheme.

60. Sinclair AJ et al. (2012). *Rating papilloedema: an evaluation of the Frisén classification in idiopathic intracranial hypertension.* J Neurol 259:1406-12. [PMID 22237821](https://pubmed.ncbi.nlm.nih.gov/22237821/)
   Evaluation of the Frisen classification in IIH specifically, including its reliability.

61. Optical Coherence Tomography Substudy Committee & NORDIC Idiopathic Intracranial Hypertension Study Group (2015). *Papilledema Outcomes from the Optical Coherence Tomography Substudy of the Idiopathic Intracranial Hypertension Treatment Trial.* Ophthalmology 122:1939-45.e2. [PMID 26198807](https://pubmed.ncbi.nlm.nih.gov/26198807/) **[CALIBRATION]**
   IIHTT OCT substudy. Calibration source for the ocular module: RNFL -175 vs -89 um, ONH volume -4.9 vs -2.1 mm3, and RGCL thinning of 3.6 vs 2.1 um — greater in the BETTER-treated arm, which this model reproduces as oedema contamination rather than damage.

62. OCT Sub-Study Committee for NORDIC Idiopathic Intracranial Hypertension Study Group et al. (2014). *Baseline OCT measurements in the idiopathic intracranial hypertension treatment trial, part I: quality control, comparisons, and variability.* Invest Ophthalmol Vis Sci 55:8180-8. [PMID 25370510](https://pubmed.ncbi.nlm.nih.gov/25370510/)
   Baseline OCT methodology and variability in the same substudy.

63. Wang JK et al. (2017). *Peripapillary Retinal Pigment Epithelium Layer Shape Changes From Acetazolamide Treatment in the Idiopathic Intracranial Hypertension Treatment Trial.* Invest Ophthalmol Vis Sci 58:2554-2565. [PMID 28492874](https://pubmed.ncbi.nlm.nih.gov/28492874/)
   Peripapillary RPE shape changes with acetazolamide treatment.

64. Kupersmith MJ et al. (2017). *The Effect of Treatment of Idiopathic Intracranial Hypertension on Prevalence of Retinal and Choroidal Folds.* Am J Ophthalmol 176:77-86. [PMID 28040526](https://pubmed.ncbi.nlm.nih.gov/28040526/)
   Effect of treatment on retinal and choroidal folds.

65. Sibony PA et al. (2015). *Retinal and Choroidal Folds in Papilledema.* Invest Ophthalmol Vis Sci 56:5670-80. [PMID 26335066](https://pubmed.ncbi.nlm.nih.gov/26335066/)
   Retinal and choroidal folds in papilloedema.

66. Prokop K et al. (2024). *Effectiveness of optic nerve sheath fenestration in preserving vision in idiopathic intracranial hypertension: an updated meta-analysis and systematic review.* Acta Neurochir (Wien) 166:476. [PMID 39585430](https://pubmed.ncbi.nlm.nih.gov/39585430/) **[STRUCTURE]**
   Updated meta-analysis of optic nerve sheath fenestration — a procedure that protects the nerve without lowering ICP, which is why it appears in this model as a modifier of the swelling drive only.

67. Anzeljc AJ et al. (2018). *A 15-year review of secondary and tertiary optic nerve sheath fenestration for idiopathic intracranial hypertension.* Orbit 37:266-272. [PMID 29313398](https://pubmed.ncbi.nlm.nih.gov/29313398/)
   Fifteen-year single-centre experience of fenestration.


## 10. CSF diversion and surgical management

68. Tsermoulas G et al. (2022). *The Birmingham Standardized Idiopathic Intracranial Hypertension Shunt Protocol: Technical Note.* World Neurosurg 167:147-151. [PMID 36089279](https://pubmed.ncbi.nlm.nih.gov/36089279/)
   A standardised shunt protocol for IIH.

69. Greener DL et al. (2020). *Idiopathic Intracranial Hypertension: Shunt Failure and the Role of Obesity.* World Neurosurg 137:e83-e88. [PMID 31954904](https://pubmed.ncbi.nlm.nih.gov/31954904/)
   Shunt failure and the role of obesity.

70. Patel J et al. (2024). *Cerebrospinal Fluid Diversion from the Cisterna Magna in Patients with Idiopathic Intracranial Hypertension and Slit Ventricles: Long-Term Effectiveness, Revision Rates, and Clinical Outcomes.* World Neurosurg 186:e326-e334. [PMID 38548048](https://pubmed.ncbi.nlm.nih.gov/38548048/)
   Long-term effectiveness and revision rates of CSF diversion in IIH with slit ventricles.


## 11. Intracranial pressure dynamics and modelling

71. Marmarou A et al. (1978). *A nonlinear analysis of the cerebrospinal fluid system and intracranial pressure dynamics.* J Neurosurg 48:332-44. [PMID 632857](https://pubmed.ncbi.nlm.nih.gov/632857/) **[CALIBRATION]** **[STRUCTURE]**
   Marmarou's nonlinear analysis of the CSF system: the pressure-volume relationship and dP/dt = E*P*(I_f - I_out) used verbatim as this model's ICP equation, and the source of the elastance value.

72. Bruce BB et al. (2016). *Quality of life at 6 months in the Idiopathic Intracranial Hypertension Treatment Trial.* Neurology 87:1871-1877. [PMID 27694262](https://pubmed.ncbi.nlm.nih.gov/27694262/)
   Quality of life at 6 months in IIHTT — the endpoint the headache module targets.


---

## Notes on how these sources are used, and where they disagree

**The pressure-per-kilogram slope is not agreed in the literature.** The
IIH:WT randomised trial implies 0.28 cmH2O/kg at 12 months and 0.31 at 24
months; the very-low-calorie-diet cohort implies 0.51; the IIHTT placebo
arm, read at face value, implies 1.52. This model uses the randomised
value (0.28) and treats the IIHTT placebo arm's excess as non-mechanistic —
see `iih_model_report.txt`, section A4, where the model's shortfall against
each IIHTT arm turns out to be nearly identical (-49.7 and -46.8 mmH2O),
which is what an arm-independent artefact looks like and what a treatment
effect does not.

**The loop gain is derived from two papers that were not written to be
combined.** The trans-stenotic and sagittal sinus pressure falls come from
PMID 29922401 (n=45) and the CSF opening-pressure fall from PMID 29871989
(n=50), both from the same centre and overlapping in time but not the same
cohort, and with different measurement intervals (immediate vs 3 months).
The derived ratio (1.53, giving gamma ~ 0.35) should be read as an
order-of-magnitude argument that gamma is not zero, which is the claim that
matters, rather than as a precise estimate. PMID 31832847 supports the same
conclusion by direct simultaneous measurement in a smaller series.

**No source measures the composition of an individual patient's pressure.**
The central bound in this model — how much of a given opening pressure is
venous floor and how much is resistive — requires catheter manometry that
none of the drug trials performed. That is the model's principal
recommendation, and its principal limitation.

