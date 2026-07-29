# Premenstrual Dysphoric Disorder (PMDD) — Annotated References

Every entry was resolved against PubMed by title and the returned record was
checked against the intended title before being listed, so each PMID link
points at the paper named beside it. Annotated entries say what the paper
contributes to *this* model; the remainder are supporting literature.

The five predictions referred to as **[P1]**-**[P5]** are defined in
[`README.md`](README.md) and asserted numerically in
[`pmdd_python_twin.py`](pmdd_python_twin.py).

---

## 1. Definition, diagnosis, epidemiology, burden

1. **Kendler KS** Longitudinal population-based twin study of retrospectively reported premenstrual symptoms and lifetime major depression. *Am J Psychiatry* 1998. [PMID 9734548](https://pubmed.ncbi.nlm.nih.gov/9734548/)
   - Twin evidence that premenstrual symptom liability is substantially heritable and partly shared with major depression.
2. **Wittchen H -U** Prevalence, incidence and stability of premenstrual dysphoric disorder in the community. *Psychol Med* 2002. [PMID 11883723](https://pubmed.ncbi.nlm.nih.gov/11883723/)
   - Community prevalence, incidence and diagnostic stability of PMDD.
3. **Eisenlohr-Moul TA** Toward the Reliable Diagnosis of DSM-5 Premenstrual Dysphoric Disorder: The Carolina Premenstrual Assessment Scoring System (C-PASS). *Am J Psychiatry* 2017. [PMID 27523500](https://pubmed.ncbi.nlm.nih.gov/27523500/)
   - The prospective two-cycle DRSP scoring rules that the model's DSM-5 >=30% criterion check imitates.
4. **Eisenlohr-Moul TA** Perimenstrual exacerbation of symptoms in borderline personality disorder: evidence from multilevel models and the Carolina Premenstrual Assessment Scoring System. *Psychol Med* 2018. [PMID 29860953](https://pubmed.ncbi.nlm.nih.gov/29860953/)
   - Perimenstrual exacerbation outside PMDD — hormone sensitivity behaves like a dimension, not a category.

---

## 2. The central experiment — hormone LEVELS are normal, SENSITIVITY is not

5. **Schmidt PJ** Differential behavioral effects of gonadal steroids in women with and in those without premenstrual syndrome. *N Engl J Med* 1998. [PMID 9435325](https://pubmed.ncbi.nlm.nih.gov/9435325/)
   - **The anchor experiment.** Ovarian suppression abolishes the symptoms; blinded add-back of estradiol OR progesterone reinstates them in PMDD and does nothing in controls. This is why every disease parameter in this model is a gain and not a hormone level.
6. **Soares CN** Reproductive hormone sensitivity and risk for depression across the female life cycle: a continuum of vulnerability?. *J Psychiatry Neurosci* 2008. [PMID 18592034](https://pubmed.ncbi.nlm.nih.gov/18592034/)
   - The reproductive-hormone-sensitivity account across the female life cycle, which places PMDD alongside perinatal and perimenopausal depression.
7. **Schmidt PJ** Premenstrual Dysphoric Disorder Symptoms Following Ovarian Suppression: Triggered by Change in Ovarian Steroid Levels But Not Continuous Stable Levels. *Am J Psychiatry* 2017. [PMID 28427285](https://pubmed.ncbi.nlm.nih.gov/28427285/)
   - Symptoms are triggered by the CHANGE in ovarian steroid levels and not by continuous stable levels — the direct empirical basis for the symmetric change detector and for prediction [P5].
8. **Dubey N** The ESC/E(Z) complex, an effector of response to ovarian steroids, manifests an intrinsic difference in cells from women with premenstrual dysphoric disorder. *Mol Psychiatry* 2017. [PMID 28044059](https://pubmed.ncbi.nlm.nih.gov/28044059/)
   - A cell-autonomous difference in the RESPONSE to ovarian steroids, demonstrated in cells cultured away from the ovary — sensitivity is intrinsic, not a consequence of exposure.
9. **Yen JY** Estrogen levels, emotion regulation, and emotional symptoms of women with premenstrual dysphoric disorder: The moderating effect of estrogen receptor 1α polymorphism. *Prog Neuropsychopharmacol Biol Psychiatry* 2018. [PMID 29146473](https://pubmed.ncbi.nlm.nih.gov/29146473/)
10. **Peters JR** Dimensional Affective Sensitivity to Hormones across the Menstrual Cycle (DASH-MC): A transdiagnostic framework for ovarian steroid influences on psychopathology. *Mol Psychiatry* 2025. [PMID 39143323](https://pubmed.ncbi.nlm.nih.gov/39143323/)
   - The DASH-MC framework: affective sensitivity to hormone CHANGE as a transdiagnostic dimension. This model is one numerical implementation of that idea.

---

## 3. Neurosteroidogenesis and allopregnanolone in PMDD

11. **Uzunov DP** Fluoxetine-elicited changes in brain neurosteroid content measured by negative ion mass fragmentography. *Proc Natl Acad Sci U S A* 1996. [PMID 8901628](https://pubmed.ncbi.nlm.nih.gov/8901628/)
   - The companion finding that fluoxetine changes brain neurosteroid content at clinically relevant exposures.
12. **Griffin LD** Selective serotonin reuptake inhibitors directly alter activity of neurosteroidogenic enzymes. *Proc Natl Acad Sci U S A* 1999. [PMID 10557352](https://pubmed.ncbi.nlm.nih.gov/10557352/)
   - SSRIs shift 3alpha-HSD toward the reductive direction and RAISE brain allopregnanolone. The mechanism that makes the sign of the SSRI's neurosteroid arm a genuine question — see [P3].
13. **Andréen L** Allopregnanolone concentration and mood--a bimodal association in postmenopausal women treated with oral progesterone. *Psychopharmacology (Berl)* 2006. [PMID 16724185](https://pubmed.ncbi.nlm.nih.gov/16724185/)
   - Allopregnanolone and mood are related BIMODALLY, not monotonically. The single most important empirical support for the inverted-U transduction in this model.
14. **Nyberg S** Allopregnanolone decrease with symptom improvement during placebo and gonadotropin-releasing hormone agonist treatment in women with severe premenstrual syndrome. *Gynecol Endocrinol* 2007. [PMID 17558683](https://pubmed.ncbi.nlm.nih.gov/17558683/)
   - Symptom improvement tracks the FALL in allopregnanolone during GnRH-agonist treatment and during placebo response.
15. **Martinez PE** 5α-Reductase Inhibition Prevents the Luteal Phase Increase in Plasma Allopregnanolone Levels and Mitigates Symptoms in Women with Premenstrual Dysphoric Disorder. *Neuropsychopharmacology* 2016. [PMID 26272051](https://pubmed.ncbi.nlm.nih.gov/26272051/)
   - Dutasteride 2.5 mg blocks the luteal allopregnanolone rise and improves symptoms; 0.5 mg does neither. The dose threshold that motivates prediction [P2] — the flat top of the inverted U.
16. **Bixo M** Treatment of premenstrual dysphoric disorder with the GABA(A) receptor modulating steroid antagonist Sepranolone (UC1010)-A randomized controlled trial. *Psychoneuroendocrinology* 2017. [PMID 28319848](https://pubmed.ncbi.nlm.nih.gov/28319848/)
   - An allopregnanolone-site ANTAGONIST improves PMDD — the mirror image of the dutasteride result, and the calibration target for the model's KI_ISO.
17. **Hantsoo L** Allopregnanolone in premenstrual dysphoric disorder (PMDD): Evidence for dysregulated sensitivity to GABA-A receptor modulating neuroactive steroids across the menstrual cycle. *Neurobiol Stress* 2020. [PMID 32435664](https://pubmed.ncbi.nlm.nih.gov/32435664/)
   - Review of the evidence that PMDD is a disorder of SENSITIVITY to GABA-A-modulating neuroactive steroids rather than of their concentrations.
18. **Hamidovic A** Neuroactive steroid hormone trajectories across the menstrual cycle in premenstrual dysphoric disorder (PMDD): the PHASE study. *Mol Psychiatry* 2024. [PMID 38664491](https://pubmed.ncbi.nlm.nih.gov/38664491/)
   - The PHASE study: neuroactive steroid trajectories across the cycle in PMDD — trajectory and sensitivity rather than absolute level.

---

## 4. GABA-A receptor pharmacology and subunit plasticity

19. **Smith SS** GABA(A) receptor alpha4 subunit suppression prevents withdrawal properties of an endogenous steroid. *Nature* 1998. [PMID 9582073](https://pubmed.ncbi.nlm.nih.gov/9582073/)
   - Suppressing the GABA-A alpha4 subunit prevents the withdrawal properties of an endogenous steroid. The experiment the model's ALPHA4 state is built on.
20. **Smith SS** Withdrawal from 3alpha-OH-5alpha-pregnan-20-One using a pseudopregnancy model alters the kinetics of hippocampal GABAA-gated current and increases the GABAA receptor alpha4 subunit in association with increased anxiety. *J Neurosci* 1998. [PMID 9651210](https://pubmed.ncbi.nlm.nih.gov/9651210/)
   - Progesterone withdrawal changes GABA-A kinetics and raises alpha4 in association with increased anxiety — the withdrawal arm of the plasticity block.
21. **Epperson CN** Cortical gamma-aminobutyric acid levels across the menstrual cycle in healthy women and those with premenstrual dysphoric disorder: a proton magnetic resonance spectroscopy study. *Arch Gen Psychiatry* 2002. [PMID 12215085](https://pubmed.ncbi.nlm.nih.gov/12215085/)
   - Cortical GABA moves in OPPOSITE directions across the cycle in PMDD and in controls (1H-MRS) — human evidence of inverted regulation, not altered supply.
22. **Maguire JL** Ovarian cycle-linked changes in GABA(A) receptors mediating tonic inhibition alter seizure susceptibility and anxiety. *Nat Neurosci* 2005. [PMID 15895085](https://pubmed.ncbi.nlm.nih.gov/15895085/)
   - Cycle-linked GABA-A subunit changes alter seizure susceptibility and anxiety in mice — the delta/alpha4 plasticity is a real, measured ovarian-cycle phenomenon.
23. **Hosie AM** Endogenous neurosteroids regulate GABAA receptors through two discrete transmembrane sites. *Nature* 2006. [PMID 17108970](https://pubmed.ncbi.nlm.nih.gov/17108970/)
   - Endogenous neurosteroids act at two discrete transmembrane sites on GABA-A — the molecular target of the load variable L.
24. **Shen H** Reversal of neurosteroid effects at alpha4beta2delta GABAA receptors triggers anxiety at puberty. *Nat Neurosci* 2007. [PMID 17351635](https://pubmed.ncbi.nlm.nih.gov/17351635/)
   - At alpha4-beta2-delta receptors the sign of neurosteroid modulation REVERSES — the receptor-level basis for a non-monotonic transduction.
25. **Maguire J** Neurosteroid synthesis-mediated regulation of GABA(A) receptors: relevance to the ovarian cycle and stress. *J Neurosci* 2007. [PMID 17329412](https://pubmed.ncbi.nlm.nih.gov/17329412/)
   - Neurosteroid-synthesis-driven regulation of GABA-A receptors across the ovarian cycle and under stress.

---

## 5. Serotonergic pharmacology — SSRIs, tryptophan, dosing strategy

26. **Rapkin AJ** Whole-blood serotonin in premenstrual syndrome. *Obstet Gynecol* 1987. [PMID 3627623](https://pubmed.ncbi.nlm.nih.gov/3627623/)
   - The origin of the luteal serotonergic-hypofunction hypothesis.
27. **Menkes DB** Acute tryptophan depletion aggravates premenstrual syndrome. *J Affect Disord* 1994. [PMID 7798465](https://pubmed.ncbi.nlm.nih.gov/7798465/)
   - Depleting tryptophan makes premenstrual symptoms worse — the serotonergic arm is causal, not incidental.
28. **Steiner M** Fluoxetine in the treatment of premenstrual dysphoria. Canadian Fluoxetine/Premenstrual Dysphoria Collaborative Study Group. *N Engl J Med* 1995. [PMID 7739706](https://pubmed.ncbi.nlm.nih.gov/7739706/)
   - Fluoxetine 20 mg in PMDD (Canadian collaborative trial) — one of the two effect sizes the SSRI arm is calibrated against.
29. **Yonkers KA** Symptomatic improvement of premenstrual dysphoric disorder with sertraline treatment. A randomized controlled trial. Sertraline Premenstrual Dysphoric Collaborative Study Group. *JAMA* 1997. [PMID 9307345](https://pubmed.ncbi.nlm.nih.gov/9307345/)
   - Sertraline in PMDD, randomised and placebo-controlled — the other SSRI calibration anchor.
30. **Steinberg S** A placebo-controlled study of the effects of L-tryptophan in patients with premenstrual dysphoria. *Adv Exp Med Biol* 1999. [PMID 10721042](https://pubmed.ncbi.nlm.nih.gov/10721042/)
31. **Freeman EW** Differential response to antidepressants in women with premenstrual syndrome/premenstrual dysphoric disorder: a randomized controlled trial. *Arch Gen Psychiatry* 1999. [PMID 10530636](https://pubmed.ncbi.nlm.nih.gov/10530636/)
   - Serotonergic versus non-serotonergic antidepressants in PMS/PMDD — the pharmacological specificity of the serotonergic arm.
32. **Yonkers KA** Symptom-Onset Dosing of Sertraline for the Treatment of Premenstrual Dysphoric Disorder: A Randomized Clinical Trial. *JAMA Psychiatry* 2015. [PMID 26351969](https://pubmed.ncbi.nlm.nih.gov/26351969/)
   - Symptom-onset dosing — the trial behind the model's symptom-onset scenario, which the model predicts is weaker than luteal-phase dosing.
33. **Barone JC** Luteal phase sertraline treatment of premenstrual dysphoric disorder (PMDD): Effects on markers of hypothalamic pituitary adrenal (HPA) axis activation and inflammation. *Psychoneuroendocrinology* 2024. [PMID 39096755](https://pubmed.ncbi.nlm.nih.gov/39096755/)

---

## 6. Ovulation suppression — combined pills, GnRH analogues, add-back, surgery

34. **Wyatt KM** The effectiveness of GnRHa with and without 'add-back' therapy in treating premenstrual syndrome: a meta analysis. *BJOG* 2004. [PMID 15198787](https://pubmed.ncbi.nlm.nih.gov/15198787/)
   - Meta-analysis of GnRH agonists with and without add-back — the effect size and the bone/vasomotor cost the model tracks as BMD and hot flushes.
35. **Cronje WH** Hysterectomy and bilateral oophorectomy for severe premenstrual syndrome. *Hum Reprod* 2004. [PMID 15229203](https://pubmed.ncbi.nlm.nih.gov/15229203/)
   - Definitive surgical treatment for severe, refractory PMS — the model's oophorectomy scenario.
36. **Pearlstein TB** Treatment of premenstrual dysphoric disorder with a new drospirenone-containing oral contraceptive formulation. *Contraception* 2005. [PMID 16307962](https://pubmed.ncbi.nlm.nih.gov/16307962/)
   - The 24/4 regimen in PMDD: shortening the hormone-free interval from 7 to 4 days is the design change, which is exactly what prediction [P4] is about.
37. **Segebladh B** Evaluation of different add-back estradiol and progesterone treatments to gonadotropin-releasing hormone agonist treatment in patients with premenstrual dysphoric disorder. *Am J Obstet Gynecol* 2009. [PMID 19398092](https://pubmed.ncbi.nlm.nih.gov/19398092/)
   - Which add-back regimens reinstate symptoms during GnRH-agonist treatment — dose guidance for the add-back scenario.
38. **Lopez LM** Oral contraceptives containing drospirenone for premenstrual syndrome. *Cochrane Database Syst Rev* 2012. [PMID 22336820](https://pubmed.ncbi.nlm.nih.gov/22336820/)
   - Cochrane review of drospirenone-containing pills for premenstrual syndrome.
39. **Fu Y** [Efficacy and safety of a combined oral contraceptive containing drospirenone 3 mg and ethinylestradiol 20 µg in the treatment of premenstrual dysphoric disorder: a randomized, double blind placebo-controlled study]. *Zhonghua Fu Chan Ke Za Zhi* 2014. [PMID 25327732](https://pubmed.ncbi.nlm.nih.gov/25327732/)
40. **Comasco E** Ulipristal Acetate for Treatment of Premenstrual Dysphoric Disorder: A Proof-of-Concept Randomized Controlled Trial. *Am J Psychiatry* 2021. [PMID 33297719](https://pubmed.ncbi.nlm.nih.gov/33297719/)
   - A selective progesterone-receptor modulator improves PMDD while sparing estradiol — progesterone-receptor signalling is part of the causal path.
41. **Ciritel AA** Effect of hormonal contraception in individuals with anxiety and mood (affective) disorders: a rapid review. *BMJ Sex Reprod Health* 2026. [PMID 42259606](https://pubmed.ncbi.nlm.nih.gov/42259606/)

---

## 7. Corticolimbic imaging, HPA axis, sleep and autonomic function

42. **Girdler SS** Dysregulation of cardiovascular and neuroendocrine responses to stress in premenstrual dysphoric disorder. *Psychiatry Res* 1998. [PMID 9858034](https://pubmed.ncbi.nlm.nih.gov/9858034/)
   - The blunted stress response of PMDD, modelled here as loss of the allopregnanolone brake on the HPA axis.
43. **Parry BL** Sleep, rhythms and women's mood. Part I. Menstrual cycle, pregnancy and postpartum. *Sleep Med Rev* 2006. [PMID 16460973](https://pubmed.ncbi.nlm.nih.gov/16460973/)
44. **Protopopescu X** Toward a functional neuroanatomy of premenstrual dysphoric disorder. *J Affect Disord* 2008. [PMID 18031826](https://pubmed.ncbi.nlm.nih.gov/18031826/)
   - The corticolimbic pattern that the AMY / PFC / REW block of this model reproduces.
45. **Baker FC** Daytime sleepiness, psychomotor performance, waking EEG spectra and evoked potentials in women with severe premenstrual syndrome. *J Sleep Res* 2010. [PMID 19840240](https://pubmed.ncbi.nlm.nih.gov/19840240/)
46. **Gingnell M** Menstrual cycle effects on amygdala reactivity to emotional stimulation in premenstrual dysphoric disorder. *Horm Behav* 2012. [PMID 22814368](https://pubmed.ncbi.nlm.nih.gov/22814368/)
   - Cycle-phase-dependent amygdala reactivity in PMDD — the imaging counterpart of the model's AMY state.
47. **Baller EB** Abnormalities of dorsolateral prefrontal function in women with premenstrual dysphoric disorder: a multimodal neuroimaging study. *Am J Psychiatry* 2013. [PMID 23361612](https://pubmed.ncbi.nlm.nih.gov/23361612/)
   - Prefrontal dysfunction in PMDD across modalities — the PFC state.
48. **Comasco E** Neuroimaging the Menstrual Cycle and Premenstrual Dysphoric Disorder. *Curr Psychiatry Rep* 2015. [PMID 26272540](https://pubmed.ncbi.nlm.nih.gov/26272540/)

---

## 8. Adjunctive and non-pharmacological treatment

49. **Wang M** Treatment of premenstrual syndrome by spironolactone: a double-blind, placebo-controlled study. *Acta Obstet Gynecol Scand* 1995. [PMID 8533564](https://pubmed.ncbi.nlm.nih.gov/8533564/)
   - Spironolactone relieves the somatic/fluid domain, which is why the model's spironolactone scenario moves S_PHY and almost nothing else.
50. **Freeman EW** A double-blind trial of oral progesterone, alprazolam, and placebo in treatment of severe premenstrual syndrome. *JAMA* 1995. [PMID 7791258](https://pubmed.ncbi.nlm.nih.gov/7791258/)
   - The alprazolam effect size, together with the negative result for oral progesterone.
51. **Thys-Jacobs S** Calcium carbonate and the premenstrual syndrome: effects on premenstrual and menstrual symptoms. Premenstrual Syndrome Study Group. *Am J Obstet Gynecol* 1998. [PMID 9731851](https://pubmed.ncbi.nlm.nih.gov/9731851/)
   - Calcium 1200 mg/d — the adjunct with the best trial evidence, and the calibration target for CA_SUPP.
52. **Schellenberg R** Treatment for the premenstrual syndrome with agnus castus fruit extract: prospective, randomised, placebo controlled study. *BMJ* 2001. [PMID 11159568](https://pubmed.ncbi.nlm.nih.gov/11159568/)
   - Vitex agnus-castus in premenstrual syndrome.
53. **Hunter MS** A randomized comparison of psychological (cognitive behavior therapy), medical (fluoxetine) and combined treatment for women with premenstrual dysphoric disorder. *J Psychosom Obstet Gynaecol* 2002. [PMID 12436805](https://pubmed.ncbi.nlm.nih.gov/12436805/)
54. **Bharati M** Comparing the Effects of Yoga & Oral Calcium Administration in Alleviating Symptoms of Premenstrual Syndrome in Medical Undergraduates. *J Caring Sci* 2016. [PMID 27752483](https://pubmed.ncbi.nlm.nih.gov/27752483/)

---

## 9. The neurosteroid contrast case — postpartum depression

55. **Meltzer-Brody S** Brexanolone injection in post-partum depression: two multicentre, double-blind, randomised, placebo-controlled, phase 3 trials. *Lancet* 2018. [PMID 30177236](https://pubmed.ncbi.nlm.nih.gov/30177236/)
   - Intravenous allopregnanolone is THERAPEUTIC in postpartum depression — a clinical existence proof for the descending limb in prediction [P3].
56. **Deligiannidis KM** Effect of Zuranolone vs Placebo in Postpartum Depression: A Randomized Clinical Trial. *JAMA Psychiatry* 2021. [PMID 34190962](https://pubmed.ncbi.nlm.nih.gov/34190962/)
   - An oral neurosteroid at a sedative dose treats a neurosteroid-withdrawal depression — the same argument as brexanolone, at a different exposure.
57. **Schiller CE** Effects of gonadal steroids on reward circuitry function and anhedonia in women with a history of postpartum depression. *J Affect Disord* 2022. [PMID 35777494](https://pubmed.ncbi.nlm.nih.gov/35777494/)

---

## 10. Quantitative modelling of the menstrual cycle, and QSP practice

58. **Röblitz S** A mathematical model of the human menstrual cycle for the administration of GnRH analogues. *J Theor Biol* 2013. [PMID 23206386](https://pubmed.ncbi.nlm.nih.gov/23206386/)
   - The published precedent for the HPO cycle engine used here, including GnRH-analogue administration.
59. **Geerts H** Quantitative Systems Pharmacology Development and Application in Neuroscience. *Handb Exp Pharmacol* 2025. [PMID 40111539](https://pubmed.ncbi.nlm.nih.gov/40111539/)
   - General QSP practice in neuroscience: how far mechanistic models can be pushed and where they stop being identifiable.

---

## 11. Tools and software

- **mrgsolve** — ODE-based PK/PD and QSP simulation in R: <https://mrgsolve.org/>
- **Graphviz** — mechanistic-map rendering (`dot -Tsvg`): <https://graphviz.org/>
- **Shiny** — interactive dashboards in R: <https://shiny.posit.co/>
- **gPKPDviz** — an mrgsolve-based PK/PD simulation Shiny tool:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/> · code at
  <https://github.com/Genentech/gPKPDviz/>

## 12. Clinical guidance

- ACOG clinical guidance on premenstrual syndrome:
  <https://www.acog.org/clinical/clinical-guidance>
- International Society for Premenstrual Disorders (ISPMD) consensus on
  classification, quantification and management: <https://ispmd.org/>
- DSM-5-TR, Premenstrual Dysphoric Disorder, criteria A-G.

_59 PubMed-verified references._
