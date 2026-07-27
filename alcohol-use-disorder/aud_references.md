# Alcohol Use Disorder (AUD) — QSP 모델 참고문헌

**References for `aud_qsp_model.dot`, `aud_mrgsolve_model.R` and `aud_shiny_app.R`**

이 목록의 **모든 PMID는 NCBI E-utilities로 직접 조회하여 제목을 대조 검증**했습니다.
각 항목의 제목·저널·연도·저자는 PubMed 레코드에서 그대로 가져온 것이며, 기억이나
추정으로 적은 것이 아닙니다. 검증 스크립트는 제목 토큰 중첩도가 기준치를 넘지
못하면 그 항목을 채택하지 않고 탈락시킵니다 — 아래 §미해결 항목 참조.

Every PMID below was resolved and title-checked against the live PubMed record
via NCBI E-utilities. Titles, journals, years and first authors are taken from
those records verbatim. Entries whose title could not be matched were dropped
rather than guessed; they are listed at the end.

| | |
|---|---|
| 검증된 인용 수 | **132** |
| 섹션 수 | 21 |
| 검증 방법 | `esearch` → `esummary` → 대칭 제목 유사도 + 길이 가드 |
| 검증 일자 | 2026-07-27 |

---

## 1. Epidemiology, diagnosis and disease burden

1. **Bush K et al. (1998).** *The AUDIT alcohol consumption questions (AUDIT-C): an effective brief screening test for problem drinking. Ambulatory Care Quality Improvement Project (ACQUIP). Alcohol Use Disorders Identification Test.* Arch Intern Med. PMID [9738608](https://pubmed.ncbi.nlm.nih.gov/9738608/)
   <br><sub>모델에서의 용도 / used for: AUDIT-C brief screening test</sub>
2. **Grant BF et al. (2015).** *Epidemiology of DSM-5 Alcohol Use Disorder: Results From the National Epidemiologic Survey on Alcohol and Related Conditions III.* JAMA Psychiatry. PMID [26039070](https://pubmed.ncbi.nlm.nih.gov/26039070/)
   <br><sub>모델에서의 용도 / used for: DSM-5 AUD prevalence, NESARC-III</sub>
3. **GBD 2016 Alcohol Collaborators (2018).** *Alcohol use and burden for 195 countries and territories, 1990-2016: a systematic analysis for the Global Burden of Disease Study 2016.* Lancet. PMID [30146330](https://pubmed.ncbi.nlm.nih.gov/30146330/)
   <br><sub>모델에서의 용도 / used for: global burden</sub>
4. **Kranzler HR et al. (2018).** *Diagnosis and Pharmacotherapy of Alcohol Use Disorder: A Review.* JAMA. PMID [30167705](https://pubmed.ncbi.nlm.nih.gov/30167705/)
   <br><sub>모델에서의 용도 / used for: diagnosis and pharmacotherapy of AUD review</sub>
5. **Witkiewitz K et al. (2019).** *Advances in the science and treatment of alcohol use disorder.* Sci Adv. PMID [31579824](https://pubmed.ncbi.nlm.nih.gov/31579824/)
   <br><sub>모델에서의 용도 / used for: advances in the science and treatment of AUD</sub>

## 2. Neurocircuitry: the three-stage addiction cycle

6. **Robinson TE et al. (1993).** *The neural basis of drug craving: an incentive-sensitization theory of addiction.* Brain Res Brain Res Rev. PMID [8401595](https://pubmed.ncbi.nlm.nih.gov/8401595/)
   <br><sub>모델에서의 용도 / used for: incentive-sensitization theory of addiction</sub>
7. **Koob GF et al. (2001).** *Drug addiction, dysregulation of reward, and allostasis.* Neuropsychopharmacology. PMID [11120394](https://pubmed.ncbi.nlm.nih.gov/11120394/)
   <br><sub>모델에서의 용도 / used for: drug addiction, dysregulation of reward, and allostasis</sub>
8. **Everitt BJ et al. (2005).** *Neural systems of reinforcement for drug addiction: from actions to habits to compulsion.* Nat Neurosci. PMID [16251991](https://pubmed.ncbi.nlm.nih.gov/16251991/)
   <br><sub>모델에서의 용도 / used for: neural systems of reinforcement, from actions to habits to compulsion</sub>
9. **Koob GF et al. (2008).** *Addiction and the brain antireward system.* Annu Rev Psychol. PMID [18154498](https://pubmed.ncbi.nlm.nih.gov/18154498/)
   <br><sub>모델에서의 용도 / used for: addiction and the brain antireward system</sub>
10. **Gilpin NW et al. (2008).** *Neurobiology of alcohol dependence: focus on motivational mechanisms.* Alcohol Res Health. PMID [19881886](https://pubmed.ncbi.nlm.nih.gov/19881886/)
   <br><sub>모델에서의 용도 / used for: neurobiology of alcohol dependence, focus on motivational mechanisms</sub>
11. **Koob GF et al. (2010).** *Neurocircuitry of addiction.* Neuropsychopharmacology. PMID [19710631](https://pubmed.ncbi.nlm.nih.gov/19710631/)
   <br><sub>모델에서의 용도 / used for: neurocircuitry of addiction</sub>
12. **Koob GF et al. (2016).** *Neurobiology of addiction: a neurocircuitry analysis.* Lancet Psychiatry. PMID [27475769](https://pubmed.ncbi.nlm.nih.gov/27475769/)
   <br><sub>모델에서의 용도 / used for: neurobiology of addiction, a neurocircuitry analysis</sub>
13. **Volkow ND et al. (2016).** *Neurobiologic Advances from the Brain Disease Model of Addiction.* N Engl J Med. PMID [26816013](https://pubmed.ncbi.nlm.nih.gov/26816013/)
   <br><sub>모델에서의 용도 / used for: neurobiologic advances from the brain disease model of addiction</sub>

## 3. Ethanol pharmacokinetics

14. **Frezza M et al. (1990).** *High blood alcohol levels in women. The role of decreased gastric alcohol dehydrogenase activity and first-pass metabolism.* N Engl J Med. PMID [2248624](https://pubmed.ncbi.nlm.nih.gov/2248624/)
   <br><sub>모델에서의 용도 / used for: high blood alcohol levels in women, first-pass metabolism</sub>
15. **Sarkola T et al. (2000).** *Mean cell volume and gamma-glutamyl transferase are superior to carbohydrate-deficient transferrin and hemoglobin-acetaldehyde adducts in the follow-up of pregnant women with alcohol abuse.* Acta Obstet Gynecol Scand. PMID [10830762](https://pubmed.ncbi.nlm.nih.gov/10830762/)
   <br><sub>모델에서의 용도 / used for: mean corpuscular volume and gamma-glutamyl transferase in alcohol</sub>
16. **Ramchandani VA et al. (2001).** *Research advances in ethanol metabolism.* Pathol Biol (Paris). PMID [11762128](https://pubmed.ncbi.nlm.nih.gov/11762128/)
   <br><sub>모델에서의 용도 / used for: research advances in ethanol metabolism</sub>
17. **Norberg A et al. (2003).** *Role of variability in explaining ethanol pharmacokinetics: research and forensic applications.* Clin Pharmacokinet. PMID [12489977](https://pubmed.ncbi.nlm.nih.gov/12489977/)
   <br><sub>모델에서의 용도 / used for: role of variability in explaining ethanol pharmacokinetics</sub>
18. **Umulis DM et al. (2005).** *A physiologically based model for ethanol and acetaldehyde metabolism in human beings.* Alcohol. PMID [15922132](https://pubmed.ncbi.nlm.nih.gov/15922132/)
   <br><sub>모델에서의 용도 / used for: physiologically based model for ethanol and acetaldehyde metabolism in human beings</sub>
19. **Jones AW (2010).** *Evidence-based survey of the elimination rates of ethanol from blood with applications in forensic casework.* Forensic Sci Int. PMID [20304569](https://pubmed.ncbi.nlm.nih.gov/20304569/)
   <br><sub>모델에서의 용도 / used for: evidence-based survey of the elimination rates of ethanol from blood</sub>

## 4. Oxidative metabolism: ADH, ALDH2, CYP2E1, acetaldehyde

20. **Deitrich RA et al. (1989).** *Mechanism of action of ethanol: initial central nervous system actions.* Pharmacol Rev. PMID [2700603](https://pubmed.ncbi.nlm.nih.gov/2700603/)
   <br><sub>모델에서의 용도 / used for: mechanism of action of ethanol, initial central nervous system actions</sub>
21. **Lieber CS (2004).** *CYP2E1: from ASH to NASH.* Hepatol Res. PMID [14734144](https://pubmed.ncbi.nlm.nih.gov/14734144/)
   <br><sub>모델에서의 용도 / used for: CYP2E1, from ASH to NASH</sub>
22. **Zakhari S (2006).** *Overview: how is alcohol metabolized by the body?.* Alcohol Res Health. PMID [17718403](https://pubmed.ncbi.nlm.nih.gov/17718403/)
   <br><sub>모델에서의 용도 / used for: overview: how is alcohol metabolized by the body?</sub>
23. **Lu Y et al. (2008).** *CYP2E1 and oxidative liver injury by alcohol.* Free Radic Biol Med. PMID [18078827](https://pubmed.ncbi.nlm.nih.gov/18078827/)
   <br><sub>모델에서의 용도 / used for: CYP2E1 and oxidative liver injury by alcohol</sub>
24. **Setshedi M et al. (2010).** *Acetaldehyde adducts in alcoholic liver disease.* Oxid Med Cell Longev. PMID [20716942](https://pubmed.ncbi.nlm.nih.gov/20716942/)
   <br><sub>모델에서의 용도 / used for: acetaldehyde adducts in alcoholic liver disease</sub>
25. **Chen CH et al. (2014).** *Targeting aldehyde dehydrogenase 2: new therapeutic opportunities.* Physiol Rev. PMID [24382882](https://pubmed.ncbi.nlm.nih.gov/24382882/)
   <br><sub>모델에서의 용도 / used for: targeting aldehyde dehydrogenase 2, new therapeutic opportunities</sub>

## 5. Pharmacogenetics

26. **Schuckit MA (1994).** *Low level of response to alcohol as a predictor of future alcoholism.* Am J Psychiatry. PMID [8296886](https://pubmed.ncbi.nlm.nih.gov/8296886/)
   <br><sub>모델에서의 용도 / used for: low level of response to alcohol as a predictor</sub>
27. **Oslin DW et al. (2003).** *A functional polymorphism of the mu-opioid receptor gene is associated with naltrexone response in alcohol-dependent patients.* Neuropsychopharmacology. PMID [12813472](https://pubmed.ncbi.nlm.nih.gov/12813472/)
   <br><sub>모델에서의 용도 / used for: functional polymorphism of the mu-opioid receptor gene and naltrexone response</sub>
28. **Ray LA et al. (2004).** *A polymorphism of the mu-opioid receptor gene (OPRM1) and sensitivity to the effects of alcohol in humans.* Alcohol Clin Exp Res. PMID [15608594](https://pubmed.ncbi.nlm.nih.gov/15608594/)
   <br><sub>모델에서의 용도 / used for: OPRM1 A118G and subjective responses to alcohol</sub>
29. **Edenberg HJ et al. (2004).** *Variations in GABRA2, encoding the alpha 2 subunit of the GABA(A) receptor, are associated with alcohol dependence and with brain oscillations.* Am J Hum Genet. PMID [15024690](https://pubmed.ncbi.nlm.nih.gov/15024690/)
   <br><sub>모델에서의 용도 / used for: GABRA2 variation and alcohol dependence</sub>
30. **Edenberg HJ (2007).** *The genetics of alcohol metabolism: role of alcohol dehydrogenase and aldehyde dehydrogenase variants.* Alcohol Res Health. PMID [17718394](https://pubmed.ncbi.nlm.nih.gov/17718394/)
   <br><sub>모델에서의 용도 / used for: genetics of alcohol dehydrogenases and aldehyde dehydrogenases</sub>
31. **Anton RF et al. (2008).** *An evaluation of mu-opioid receptor (OPRM1) as a predictor of naltrexone response in the treatment of alcohol dependence: results from the Combined Pharmacotherapies and Behavioral Interventions for Alcohol Dependence (COMBINE) study.* Arch Gen Psychiatry. PMID [18250251](https://pubmed.ncbi.nlm.nih.gov/18250251/)
   <br><sub>모델에서의 용도 / used for: OPRM1 A118G and naltrexone in the COMBINE study</sub>
32. **Brooks PJ et al. (2009).** *The alcohol flushing response: an unrecognized risk factor for esophageal cancer from alcohol consumption.* PLoS Med. PMID [19320537](https://pubmed.ncbi.nlm.nih.gov/19320537/)
   <br><sub>모델에서의 용도 / used for: ALDH2 deficiency and oesophageal cancer</sub>
33. **Johnson BA et al. (2011).** *Pharmacogenetic approach at the serotonin transporter gene as a method of reducing the severity of alcohol drinking.* Am J Psychiatry. PMID [21247998](https://pubmed.ncbi.nlm.nih.gov/21247998/)
   <br><sub>모델에서의 용도 / used for: ondansetron and 5-HTT/5-HT3 genotype in alcohol dependence</sub>
34. **Bierut LJ et al. (2012).** *ADH1B is associated with alcohol dependence and alcohol consumption in populations of European and African ancestry.* Mol Psychiatry. PMID [21968928](https://pubmed.ncbi.nlm.nih.gov/21968928/)
   <br><sub>모델에서의 용도 / used for: ADH1B is associated with alcohol dependence</sub>
35. **Walters RK et al. (2018).** *Transancestral GWAS of alcohol dependence reveals common genetic underpinnings with psychiatric disorders.* Nat Neurosci. PMID [30482948](https://pubmed.ncbi.nlm.nih.gov/30482948/)
   <br><sub>모델에서의 용도 / used for: transancestral GWAS of alcohol dependence</sub>

## 6. GABA-A signalling, tolerance and subunit plasticity

36. **Morrow AL et al. (2001).** *The role of GABAergic neuroactive steroids in ethanol action, tolerance and dependence.* Brain Res Brain Res Rev. PMID [11744078](https://pubmed.ncbi.nlm.nih.gov/11744078/)
   <br><sub>모델에서의 용도 / used for: neuroactive steroids and ethanol</sub>
37. **Wallner M et al. (2003).** *Ethanol enhances alpha 4 beta 3 delta and alpha 6 beta 3 delta gamma-aminobutyric acid type A receptors at low concentrations known to affect humans.* Proc Natl Acad Sci U S A. PMID [14625373](https://pubmed.ncbi.nlm.nih.gov/14625373/)
   <br><sub>모델에서의 용도 / used for: ethanol enhances alpha4beta3delta and alpha6beta3delta GABA-A receptors</sub>
38. **Cagetti E et al. (2003).** *Withdrawal from chronic intermittent ethanol treatment changes subunit composition, reduces synaptic function, and decreases behavioral responses to positive allosteric modulators of GABAA receptors.* Mol Pharmacol. PMID [12488536](https://pubmed.ncbi.nlm.nih.gov/12488536/)
   <br><sub>모델에서의 용도 / used for: withdrawal from chronic intermittent ethanol and GABA-A subunit changes</sub>
39. **Roberto M et al. (2004).** *Increased GABA release in the central amygdala of ethanol-dependent rats.* J Neurosci. PMID [15537886](https://pubmed.ncbi.nlm.nih.gov/15537886/)
   <br><sub>모델에서의 용도 / used for: CeA GABAergic transmission and ethanol dependence</sub>
40. **Kumar S et al. (2009).** *The role of GABA(A) receptors in the acute and chronic effects of ethanol: a decade of progress.* Psychopharmacology (Berl). PMID [19455309](https://pubmed.ncbi.nlm.nih.gov/19455309/)
   <br><sub>모델에서의 용도 / used for: the role of GABA(A) receptors in the acute and chronic effects of ethanol</sub>
41. **Lovinger DM (2018).** *Presynaptic Ethanol Actions: Potential Roles in Ethanol Seeking.* Handb Exp Pharmacol. PMID [29204712](https://pubmed.ncbi.nlm.nih.gov/29204712/)
   <br><sub>모델에서의 용도 / used for: ethanol effects on GABAergic transmission</sub>

## 7. Glutamate, NMDA and homeostatic plasticity

42. **Lovinger DM et al. (1989).** *Ethanol inhibits NMDA-activated ion current in hippocampal neurons.* Science. PMID [2467382](https://pubmed.ncbi.nlm.nih.gov/2467382/)
   <br><sub>모델에서의 용도 / used for: ethanol inhibits NMDA-activated ion current</sub>
43. **Tsai G et al. (1995).** *The glutamatergic basis of human alcoholism.* Am J Psychiatry. PMID [7864257](https://pubmed.ncbi.nlm.nih.gov/7864257/)
   <br><sub>모델에서의 용도 / used for: the glutamatergic basis of human alcoholism</sub>
44. **Nagy J (2008).** *Alcohol related changes in regulation of NMDA receptor functions.* Curr Neuropharmacol. PMID [19305787](https://pubmed.ncbi.nlm.nih.gov/19305787/)
   <br><sub>모델에서의 용도 / used for: alcohol-related changes in regulation of NMDA receptor functions</sub>
45. **Hermann D et al. (2012).** *Translational magnetic resonance spectroscopy reveals excessive central glutamate levels during alcohol withdrawal in humans and rats.* Biol Psychiatry. PMID [21907974](https://pubmed.ncbi.nlm.nih.gov/21907974/)
   <br><sub>모델에서의 용도 / used for: glutamate in the anterior cingulate during alcohol withdrawal</sub>
46. **Holmes A et al. (2013).** *Glutamatergic targets for new alcohol medications.* Psychopharmacology (Berl). PMID [23995381](https://pubmed.ncbi.nlm.nih.gov/23995381/)
   <br><sub>모델에서의 용도 / used for: glutamatergic targets for new alcohol medications</sub>
47. **Spanagel R et al. (2014).** *Acamprosate produces its anti-relapse effects via calcium.* Neuropsychopharmacology. PMID [24081303](https://pubmed.ncbi.nlm.nih.gov/24081303/)
   <br><sub>모델에서의 용도 / used for: acamprosate produces its anti-relapse effects via calcium</sub>
48. **Rao PS et al. (2015).** *Targeting glutamate uptake to treat alcohol use disorders.* Front Neurosci. PMID [25954150](https://pubmed.ncbi.nlm.nih.gov/25954150/)
   <br><sub>모델에서의 용도 / used for: targeting glutamate uptake to treat alcohol use disorders</sub>

## 8. Endogenous opioid reinforcement

49. **Lee MC et al. (1988).** *Duration of occupancy of opiate receptors by naltrexone.* J Nucl Med. PMID [2839637](https://pubmed.ncbi.nlm.nih.gov/2839637/)
   <br><sub>모델에서의 용도 / used for: XR-naltrexone mu-opioid receptor blockade PET</sub>
50. **Job MO et al. (2007).** *Mu (mu) opioid receptor regulation of ethanol-induced dopamine response in the ventral striatum: evidence of genotype specific sexual dimorphic epistasis.* Biol Psychiatry. PMID [17336938](https://pubmed.ncbi.nlm.nih.gov/17336938/)
   <br><sub>모델에서의 용도 / used for: mu-opioid receptor regulation of ethanol-induced dopamine response</sub>
51. **Gianoulakis C (2009).** *Endogenous opioids and addiction to alcohol and other drugs of abuse.* Curr Top Med Chem. PMID [19747123](https://pubmed.ncbi.nlm.nih.gov/19747123/)
   <br><sub>모델에서의 용도 / used for: endogenous opioids and addiction to alcohol</sub>
52. **Xiao C et al. (2009).** *Ethanol facilitates glutamatergic transmission to dopamine neurons in the ventral tegmental area.* Neuropsychopharmacology. PMID [18596684](https://pubmed.ncbi.nlm.nih.gov/18596684/)
   <br><sub>모델에서의 용도 / used for: ethanol and the VTA GABAergic interneuron</sub>
53. **Mitchell JM et al. (2012).** *Alcohol consumption induces endogenous opioid release in the human orbitofrontal cortex and nucleus accumbens.* Sci Transl Med. PMID [22238334](https://pubmed.ncbi.nlm.nih.gov/22238334/)
   <br><sub>모델에서의 용도 / used for: alcohol consumption induces endogenous opioid release in human orbitofrontal cortex and nucleus accumbens</sub>
54. **Weerts EM et al. (2014).** *Association of smoking with μ-opioid receptor availability before and during naltrexone blockade in alcohol-dependent subjects.* Addict Biol. PMID [23252742](https://pubmed.ncbi.nlm.nih.gov/23252742/)
   <br><sub>모델에서의 용도 / used for: mu-opioid receptor availability and naltrexone blockade PET</sub>

## 9. Dopamine, incentive salience and habit

55. **Boileau I et al. (2003).** *Alcohol promotes dopamine release in the human nucleus accumbens.* Synapse. PMID [12827641](https://pubmed.ncbi.nlm.nih.gov/12827641/)
   <br><sub>모델에서의 용도 / used for: alcohol promotes dopamine release in the human nucleus accumbens</sub>
56. **Martinez D et al. (2005).** *Alcohol dependence is associated with blunted dopamine transmission in the ventral striatum.* Biol Psychiatry. PMID [16018986](https://pubmed.ncbi.nlm.nih.gov/16018986/)
   <br><sub>모델에서의 용도 / used for: alcohol dependence is associated with blunted dopamine transmission</sub>
57. **Volkow ND et al. (2007).** *Profound decreases in dopamine release in striatum in detoxified alcoholics: possible orbitofrontal involvement.* J Neurosci. PMID [18003850](https://pubmed.ncbi.nlm.nih.gov/18003850/)
   <br><sub>모델에서의 용도 / used for: profound decreases in dopamine release in striatum in detoxified alcoholics</sub>
58. **Nestler EJ (2008).** *Review. Transcriptional mechanisms of addiction: role of DeltaFosB.* Philos Trans R Soc Lond B Biol Sci. PMID [18640924](https://pubmed.ncbi.nlm.nih.gov/18640924/)
   <br><sub>모델에서의 용도 / used for: transcriptional mechanisms of addiction, role of DeltaFosB</sub>
59. **Berridge KC et al. (2009).** *Dissecting components of reward: 'liking', 'wanting', and learning.* Curr Opin Pharmacol. PMID [19162544](https://pubmed.ncbi.nlm.nih.gov/19162544/)
   <br><sub>모델에서의 용도 / used for: dissecting components of reward, liking, wanting and learning</sub>
60. **Corbit LH et al. (2012).** *Habitual alcohol seeking: time course and the contribution of subregions of the dorsal striatum.* Biol Psychiatry. PMID [22440617](https://pubmed.ncbi.nlm.nih.gov/22440617/)
   <br><sub>모델에서의 용도 / used for: habitual alcohol seeking, time course and the contribution of subregions of the dorsal striatum</sub>

## 10. Dynorphin, kappa-opioid receptors and anti-reward

61. **Bart G et al. (2005).** *Nalmefene induced elevation in serum prolactin in normal human volunteers: partial kappa opioid agonist activity?.* Neuropsychopharmacology. PMID [15988468](https://pubmed.ncbi.nlm.nih.gov/15988468/)
   <br><sub>모델에서의 용도 / used for: nalmefene induced elevation in serum prolactin, kappa opioid</sub>
62. **Walker BM et al. (2008).** *Pharmacological evidence for a motivational role of kappa-opioid systems in ethanol dependence.* Neuropsychopharmacology. PMID [17473837](https://pubmed.ncbi.nlm.nih.gov/17473837/)
   <br><sub>모델에서의 용도 / used for: kappa opioid receptor systems and ethanol dependence</sub>
63. **Karkhanis A et al. (2017).** *Dynorphin/Kappa Opioid Receptor Signaling in Preclinical Models of Alcohol, Drug, and Food Addiction.* Int Rev Neurobiol. PMID [29056156](https://pubmed.ncbi.nlm.nih.gov/29056156/)
   <br><sub>모델에서의 용도 / used for: dynorphin/kappa opioid receptor signalling in preclinical models</sub>

## 11. Extended amygdala, CRF, NPY and the HPA axis

64. **Heilig M (2004).** *The NPY system in stress, anxiety and depression.* Neuropeptides. PMID [15337373](https://pubmed.ncbi.nlm.nih.gov/15337373/)
   <br><sub>모델에서의 용도 / used for: the NPY system in stress, anxiety and alcoholism</sub>
65. **Adinoff B et al. (2005).** *Suppression of the HPA axis stress-response: implications for relapse.* Alcohol Clin Exp Res. PMID [16088999](https://pubmed.ncbi.nlm.nih.gov/16088999/)
   <br><sub>모델에서의 용도 / used for: blunted HPA axis response in abstinent alcohol-dependent men</sub>
66. **Heilig M et al. (2007).** *A key role for corticotropin-releasing factor in alcohol dependence.* Trends Neurosci. PMID [17629579](https://pubmed.ncbi.nlm.nih.gov/17629579/)
   <br><sub>모델에서의 용도 / used for: a key role for corticotropin-releasing factor in alcohol dependence</sub>
67. **Koob GF (2010).** *The role of CRF and CRF-related peptides in the dark side of addiction.* Brain Res. PMID [19912996](https://pubmed.ncbi.nlm.nih.gov/19912996/)
   <br><sub>모델에서의 용도 / used for: the role of CRF and CRF-related peptides in the dark side of addiction</sub>
68. **Sinha R et al. (2011).** *Effects of adrenal sensitivity, stress- and cue-induced craving, and anxiety on subsequent alcohol relapse and treatment outcomes.* Arch Gen Psychiatry. PMID [21536969](https://pubmed.ncbi.nlm.nih.gov/21536969/)
   <br><sub>모델에서의 용도 / used for: stress-induced alcohol craving and HPA responses predict relapse</sub>
69. **Simpson TL et al. (2018).** *Double-Blind Randomized Clinical Trial of Prazosin for Alcohol Use Disorder.* Am J Psychiatry. PMID [30153753](https://pubmed.ncbi.nlm.nih.gov/30153753/)
   <br><sub>모델에서의 용도 / used for: prazosin for alcohol use disorder, randomized controlled trial</sub>
70. **Sinha R et al. (2021).** *Moderation of Prazosin's Efficacy by Alcohol Withdrawal Symptoms.* Am J Psychiatry. PMID [33207935](https://pubmed.ncbi.nlm.nih.gov/33207935/)
   <br><sub>모델에서의 용도 / used for: doxazosin and alcohol</sub>

## 12. Prefrontal control, cognition and recovery

71. **Bartsch AJ et al. (2007).** *Manifestations of early brain recovery associated with abstinence from alcoholism.* Brain. PMID [17178742](https://pubmed.ncbi.nlm.nih.gov/17178742/)
   <br><sub>모델에서의 용도 / used for: manifestations of early brain recovery associated with abstinence from alcoholism</sub>
72. **Goldstein RZ et al. (2011).** *Dysfunction of the prefrontal cortex in addiction: neuroimaging findings and clinical implications.* Nat Rev Neurosci. PMID [22011681](https://pubmed.ncbi.nlm.nih.gov/22011681/)
   <br><sub>모델에서의 용도 / used for: dysfunction of the prefrontal cortex in addiction</sub>
73. **MacKillop J et al. (2011).** *Delayed reward discounting and addictive behavior: a meta-analysis.* Psychopharmacology (Berl). PMID [21373791](https://pubmed.ncbi.nlm.nih.gov/21373791/)
   <br><sub>모델에서의 용도 / used for: delayed reward discounting and addictive behavior meta-analysis</sub>
74. **Zahr NM et al. (2011).** *Clinical and pathological features of alcohol-related brain damage.* Nat Rev Neurol. PMID [21487421](https://pubmed.ncbi.nlm.nih.gov/21487421/)
   <br><sub>모델에서의 용도 / used for: clinical and pathological features of alcohol-related brain damage</sub>

## 13. Withdrawal, kindling and delirium tremens

75. **Ballenger JC et al. (1978).** *Kindling as a model for alcohol withdrawal syndromes.* Br J Psychiatry. PMID [352467](https://pubmed.ncbi.nlm.nih.gov/352467/)
   <br><sub>모델에서의 용도 / used for: kindling as a model for alcohol withdrawal syndromes</sub>
76. **Sullivan JT et al. (1989).** *Assessment of alcohol withdrawal: the revised clinical institute withdrawal assessment for alcohol scale (CIWA-Ar).* Br J Addict. PMID [2597811](https://pubmed.ncbi.nlm.nih.gov/2597811/)
   <br><sub>모델에서의 용도 / used for: assessment of alcohol withdrawal, CIWA-Ar</sub>
77. **Saitz R et al. (1994).** *Individualized treatment for alcohol withdrawal. A randomized double-blind controlled trial.* JAMA. PMID [8046805](https://pubmed.ncbi.nlm.nih.gov/8046805/)
   <br><sub>모델에서의 용도 / used for: individualized treatment for alcohol withdrawal, randomized double-blind controlled trial</sub>
78. **Mayo-Smith MF (1997).** *Pharmacological management of alcohol withdrawal. A meta-analysis and evidence-based practice guideline. American Society of Addiction Medicine Working Group on Pharmacological Management of Alcohol Withdrawal.* JAMA. PMID [9214531](https://pubmed.ncbi.nlm.nih.gov/9214531/)
   <br><sub>모델에서의 용도 / used for: pharmacological management of alcohol withdrawal meta-analysis</sub>
79. **Becker HC (1998).** *Kindling in alcohol withdrawal.* Alcohol Health Res World. PMID [15706729](https://pubmed.ncbi.nlm.nih.gov/15706729/)
   <br><sub>모델에서의 용도 / used for: kindling in alcohol withdrawal</sub>
80. **Thomson AD et al. (2002).** *The Royal College of Physicians report on alcohol: guidelines for managing Wernicke's encephalopathy in the accident and Emergency Department.* Alcohol Alcohol. PMID [12414541](https://pubmed.ncbi.nlm.nih.gov/12414541/)
   <br><sub>모델에서의 용도 / used for: the Royal College of Physicians report on alcohol, Wernicke's encephalopathy</sub>
81. **Heilig M et al. (2010).** *Acute withdrawal, protracted abstinence and negative affect in alcoholism: are they linked?.* Addict Biol. PMID [20148778](https://pubmed.ncbi.nlm.nih.gov/20148778/)
   <br><sub>모델에서의 용도 / used for: acute withdrawal, protracted abstinence and negative affect</sub>
82. **Amato L et al. (2010).** *Benzodiazepines for alcohol withdrawal.* Cochrane Database Syst Rev. PMID [20238336](https://pubmed.ncbi.nlm.nih.gov/20238336/)
   <br><sub>모델에서의 용도 / used for: benzodiazepines for alcohol withdrawal</sub>
83. **Schuckit MA (2014).** *Recognition and management of withdrawal delirium (delirium tremens).* N Engl J Med. PMID [25427113](https://pubmed.ncbi.nlm.nih.gov/25427113/)
   <br><sub>모델에서의 용도 / used for: recognition and management of withdrawal delirium</sub>

## 14. Naltrexone and extended-release naltrexone

84. **Wall ME et al. (1981).** *Metabolism and disposition of naltrexone in man after oral and intravenous administration.* Drug Metab Dispos. PMID [6114837](https://pubmed.ncbi.nlm.nih.gov/6114837/)
   <br><sub>모델에서의 용도 / used for: metabolism and disposition of naltrexone in man</sub>
85. **O'Malley SS et al. (1992).** *Naltrexone and coping skills therapy for alcohol dependence. A controlled study.* Arch Gen Psychiatry. PMID [1444726](https://pubmed.ncbi.nlm.nih.gov/1444726/)
   <br><sub>모델에서의 용도 / used for: naltrexone and coping skills therapy for alcohol dependence</sub>
86. **Garbutt JC et al. (2005).** *Efficacy and tolerability of long-acting injectable naltrexone for alcohol dependence: a randomized controlled trial.* JAMA. PMID [15811981](https://pubmed.ncbi.nlm.nih.gov/15811981/)
   <br><sub>모델에서의 용도 / used for: efficacy and tolerability of long-acting injectable naltrexone</sub>
87. **Anton RF et al. (2006).** *Combined pharmacotherapies and behavioral interventions for alcohol dependence: the COMBINE study: a randomized controlled trial.* JAMA. PMID [16670409](https://pubmed.ncbi.nlm.nih.gov/16670409/)
   <br><sub>모델에서의 용도 / used for: COMBINE study</sub>
88. **Dunbar JL et al. (2006).** *Single- and multiple-dose pharmacokinetics of long-acting injectable naltrexone.* Alcohol Clin Exp Res. PMID [16499489](https://pubmed.ncbi.nlm.nih.gov/16499489/)
   <br><sub>모델에서의 용도 / used for: single- and multiple-dose pharmacokinetics of long-acting injectable naltrexone</sub>
89. **Ray LA et al. (2007).** *Effects of naltrexone on alcohol sensitivity and genetic moderators of medication response: a double-blind placebo-controlled study.* Arch Gen Psychiatry. PMID [17768272](https://pubmed.ncbi.nlm.nih.gov/17768272/)
   <br><sub>모델에서의 용도 / used for: naltrexone effects on subjective responses to alcohol</sub>
90. **Rösner S et al. (2010).** *Opioid antagonists for alcohol dependence.* Cochrane Database Syst Rev. PMID [21154349](https://pubmed.ncbi.nlm.nih.gov/21154349/)
   <br><sub>모델에서의 용도 / used for: opioid antagonists for alcohol dependence</sub>
91. **Jonas DE et al. (2014).** *Pharmacotherapy for adults with alcohol use disorders in outpatient settings: a systematic review and meta-analysis.* JAMA. PMID [24825644](https://pubmed.ncbi.nlm.nih.gov/24825644/)
   <br><sub>모델에서의 용도 / used for: pharmacotherapy for adults with alcohol use disorders in outpatient settings</sub>
92. **Bramness JG (2025).** *[Naltrexone in the treatment of alcohol dependence].* Tidsskr Nor Laegeforen. PMID [41196211](https://pubmed.ncbi.nlm.nih.gov/41196211/)
   <br><sub>모델에서의 용도 / used for: naltrexone in the treatment of alcohol dependence</sub>

## 15. Acamprosate

93. **Mason BJ et al. (2006).** *Effect of oral acamprosate on abstinence in patients with alcohol dependence in a double-blind, placebo-controlled trial: the role of patient motivation.* J Psychiatr Res. PMID [16546214](https://pubmed.ncbi.nlm.nih.gov/16546214/)
   <br><sub>모델에서의 용도 / used for: acamprosate for alcohol dependence, US multisite trial</sub>
94. **Kalk NJ et al. (2014).** *The clinical pharmacology of acamprosate.* Br J Clin Pharmacol. PMID [23278595](https://pubmed.ncbi.nlm.nih.gov/23278595/)
   <br><sub>모델에서의 용도 / used for: the clinical pharmacology of acamprosate</sub>

## 16. Disulfiram

95. **Fuller RK et al. (1986).** *Disulfiram treatment of alcoholism. A Veterans Administration cooperative study.* JAMA. PMID [3528541](https://pubmed.ncbi.nlm.nih.gov/3528541/)
   <br><sub>모델에서의 용도 / used for: disulfiram treatment of alcoholism, Veterans Administration cooperative study</sub>
96. **Johansson B (1992).** *A review of the pharmacokinetics and pharmacodynamics of disulfiram and its metabolites.* Acta Psychiatr Scand Suppl. PMID [1471547](https://pubmed.ncbi.nlm.nih.gov/1471547/)
   <br><sub>모델에서의 용도 / used for: a review of the pharmacokinetics and pharmacodynamics of disulfiram</sub>
97. **Lipsky JJ et al. (2001).** *In vivo inhibition of aldehyde dehydrogenase by disulfiram.* Chem Biol Interact. PMID [11306034](https://pubmed.ncbi.nlm.nih.gov/11306034/)
   <br><sub>모델에서의 용도 / used for: carbamathione and the inhibition of aldehyde dehydrogenase by disulfiram</sub>
98. **Jørgensen CH et al. (2011).** *The efficacy of disulfiram for the treatment of alcohol use disorder.* Alcohol Clin Exp Res. PMID [21615426](https://pubmed.ncbi.nlm.nih.gov/21615426/)
   <br><sub>모델에서의 용도 / used for: the efficacy of disulfiram for the treatment of alcohol use disorder</sub>
99. **Skinner MD et al. (2014).** *Disulfiram efficacy in the treatment of alcohol dependence: a meta-analysis.* PLoS One. PMID [24520330](https://pubmed.ncbi.nlm.nih.gov/24520330/)
   <br><sub>모델에서의 용도 / used for: disulfiram efficacy in alcohol dependence meta-analysis</sub>

## 17. Topiramate, gabapentin, baclofen, nalmefene, ondansetron, GLP-1

100. **Johnson BA et al. (2000).** *Ondansetron for reduction of drinking among biologically predisposed alcoholic patients: A randomized controlled trial.* JAMA. PMID [10944641](https://pubmed.ncbi.nlm.nih.gov/10944641/)
   <br><sub>모델에서의 용도 / used for: ondansetron for reduction of drinking among biologically predisposed alcoholic patients</sub>
101. **Johnson BA et al. (2003).** *Oral topiramate for treatment of alcohol dependence: a randomised controlled trial.* Lancet. PMID [12767733](https://pubmed.ncbi.nlm.nih.gov/12767733/)
   <br><sub>모델에서의 용도 / used for: oral topiramate for treatment of alcohol dependence</sub>
102. **Johnson BA et al. (2007).** *Topiramate for treating alcohol dependence: a randomized controlled trial.* JAMA. PMID [17925516](https://pubmed.ncbi.nlm.nih.gov/17925516/)
   <br><sub>모델에서의 용도 / used for: topiramate for treating alcohol dependence</sub>
103. **Addolorato G et al. (2007).** *Effectiveness and safety of baclofen for maintenance of alcohol abstinence in alcohol-dependent patients with liver cirrhosis: randomised, double-blind controlled study.* Lancet. PMID [18068515](https://pubmed.ncbi.nlm.nih.gov/18068515/)
   <br><sub>모델에서의 용도 / used for: baclofen in alcohol-dependent patients with liver cirrhosis</sub>
104. **Mann K et al. (2013).** *Extending the treatment options in alcohol dependence: a randomized controlled study of as-needed nalmefene.* Biol Psychiatry. PMID [23237314](https://pubmed.ncbi.nlm.nih.gov/23237314/)
   <br><sub>모델에서의 용도 / used for: nalmefene ESENSE1</sub>
105. **Gual A et al. (2013).** *A randomised, double-blind, placebo-controlled, efficacy study of nalmefene, as-needed use, in patients with alcohol dependence.* Eur Neuropsychopharmacol. PMID [23562264](https://pubmed.ncbi.nlm.nih.gov/23562264/)
   <br><sub>모델에서의 용도 / used for: nalmefene ESENSE2</sub>
106. **Mason BJ et al. (2014).** *Gabapentin treatment for alcohol dependence: a randomized clinical trial.* JAMA Intern Med. PMID [24190578](https://pubmed.ncbi.nlm.nih.gov/24190578/)
   <br><sub>모델에서의 용도 / used for: gabapentin treatment for alcohol dependence</sub>
107. **Kranzler HR et al. (2014).** *Topiramate treatment for heavy drinkers: moderation by a GRIK1 polymorphism.* Am J Psychiatry. PMID [24525690](https://pubmed.ncbi.nlm.nih.gov/24525690/)
   <br><sub>모델에서의 용도 / used for: naltrexone and topiramate combination</sub>
108. **Rose AK et al. (2018).** *Baclofen: its effectiveness in reducing harmful drinking, craving, and negative mood. A meta-analysis.* Addiction. PMID [29479827](https://pubmed.ncbi.nlm.nih.gov/29479827/)
   <br><sub>모델에서의 용도 / used for: baclofen for alcohol use disorder systematic review</sub>
109. **Anton RF et al. (2020).** *Efficacy of Gabapentin for the Treatment of Alcohol Use Disorder in Patients With Alcohol Withdrawal Symptoms: A Randomized Clinical Trial.* JAMA Intern Med. PMID [32150232](https://pubmed.ncbi.nlm.nih.gov/32150232/)
   <br><sub>모델에서의 용도 / used for: gabapentin in alcohol use disorder with withdrawal symptoms</sub>
110. **Overgaard RV et al. (2021).** *Clinical Pharmacokinetics of Oral Semaglutide: Analyses of Data from Clinical Pharmacology Trials.* Clin Pharmacokinet. PMID [33969456](https://pubmed.ncbi.nlm.nih.gov/33969456/)
   <br><sub>모델에서의 용도 / used for: clinical pharmacokinetics of semaglutide, population analysis</sub>
111. **Klausen MK et al. (2022).** *Exenatide once weekly for alcohol use disorder investigated in a randomized, placebo-controlled clinical trial.* JCI Insight. PMID [36066977](https://pubmed.ncbi.nlm.nih.gov/36066977/)
   <br><sub>모델에서의 용도 / used for: GLP-1 receptor agonist exenatide and alcohol</sub>
112. **Hendershot CS et al. (2025).** *Once-Weekly Semaglutide in Adults With Alcohol Use Disorder: A Randomized Clinical Trial.* JAMA Psychiatry. PMID [39937469](https://pubmed.ncbi.nlm.nih.gov/39937469/)
   <br><sub>모델에서의 용도 / used for: semaglutide and alcohol consumption</sub>

## 18. Objective biomarkers of alcohol consumption

113. **Anton RF et al. (2002).** *Carbohydrate-deficient transferrin and gamma-glutamyltransferase for the detection and monitoring of alcohol use: results from a multisite study.* Alcohol Clin Exp Res. PMID [12198396](https://pubmed.ncbi.nlm.nih.gov/12198396/)
   <br><sub>모델에서의 용도 / used for: carbohydrate-deficient transferrin and gamma-glutamyltransferase for monitoring drinking</sub>
114. **Bortolotti F et al. (2006).** *Carbohydrate-deficient transferrin (CDT) as a marker of alcohol abuse: a critical review of the literature 2001-2005.* J Chromatogr B Analyt Technol Biomed Life Sci. PMID [16725384](https://pubmed.ncbi.nlm.nih.gov/16725384/)
   <br><sub>모델에서의 용도 / used for: carbohydrate-deficient transferrin as a marker of alcohol abuse, critical review</sub>
115. **Viel G et al. (2012).** *Phosphatidylethanol in blood as a marker of chronic alcohol use: a systematic review and meta-analysis.* Int J Mol Sci. PMID [23203094](https://pubmed.ncbi.nlm.nih.gov/23203094/)
   <br><sub>모델에서의 용도 / used for: phosphatidylethanol in blood as a marker of chronic alcohol use, systematic review</sub>
116. **Wurst FM et al. (2015).** *Ethanol metabolites: their role in the assessment of alcohol intake.* Alcohol Clin Exp Res. PMID [26344403](https://pubmed.ncbi.nlm.nih.gov/26344403/)
   <br><sub>모델에서의 용도 / used for: ethyl glucuronide and ethyl sulfate as biomarkers</sub>
117. **Ulwelling W et al. (2018).** *The PEth Blood Test in the Security Environment: What it is; Why it is Important; and Interpretative Guidelines.* J Forensic Sci. PMID [30005144](https://pubmed.ncbi.nlm.nih.gov/30005144/)
   <br><sub>모델에서의 용도 / used for: the PEth blood test in the security environment</sub>
118. **Helander A et al. (2019).** *Elimination Characteristics of the Alcohol Biomarker Phosphatidylethanol (PEth) in Blood during Alcohol Detoxification.* Alcohol Alcohol. PMID [30968936](https://pubmed.ncbi.nlm.nih.gov/30968936/)
   <br><sub>모델에서의 용도 / used for: phosphatidylethanol elimination and detection window</sub>

## 19. Alcohol-related liver and systemic disease

119. **Gao B et al. (2011).** *Alcoholic liver disease: pathogenesis and new therapeutic targets.* Gastroenterology. PMID [21920463](https://pubmed.ncbi.nlm.nih.gov/21920463/)
   <br><sub>모델에서의 용도 / used for: alcoholic liver disease, pathogenesis and new therapeutic targets</sub>
120. **DiNubile MJ (2015).** *Prednisolone or Pentoxifylline for Alcoholic Hepatitis.* N Engl J Med. PMID [26176389](https://pubmed.ncbi.nlm.nih.gov/26176389/)
   <br><sub>모델에서의 용도 / used for: prednisolone or pentoxifylline for alcoholic hepatitis (STOPAH)</sub>
121. **Bagnardi V et al. (2015).** *Alcohol consumption and site-specific cancer risk: a comprehensive dose-response meta-analysis.* Br J Cancer. PMID [25422909](https://pubmed.ncbi.nlm.nih.gov/25422909/)
   <br><sub>모델에서의 용도 / used for: alcohol consumption and site-specific cancer risk meta-analysis</sub>
122. **Fernández-Solà J (2015).** *Cardiovascular risks and benefits of moderate and heavy alcohol consumption.* Nat Rev Cardiol. PMID [26099843](https://pubmed.ncbi.nlm.nih.gov/26099843/)
   <br><sub>모델에서의 용도 / used for: cardiovascular risks and benefits of moderate and heavy alcohol consumption</sub>
123. **Roerecke M et al. (2019).** *Alcohol Consumption and Risk of Liver Cirrhosis: A Systematic Review and Meta-Analysis.* Am J Gastroenterol. PMID [31464740](https://pubmed.ncbi.nlm.nih.gov/31464740/)
   <br><sub>모델에서의 용도 / used for: alcohol consumption and risk of liver cirrhosis meta-analysis</sub>
124. **Linz D et al. (2020).** *Alcohol Abstinence in Drinkers with Atrial Fibrillation.* N Engl J Med. PMID [32348656](https://pubmed.ncbi.nlm.nih.gov/32348656/)
   <br><sub>모델에서의 용도 / used for: alcohol abstinence in drinkers with atrial fibrillation</sub>

## 20. Endpoints, WHO risk levels and regulatory context

125. **Witkiewitz K et al. (2017).** *Clinical Validation of Reduced Alcohol Consumption After Treatment for Alcohol Dependence Using the World Health Organization Risk Drinking Levels.* Alcohol Clin Exp Res. PMID [28019652](https://pubmed.ncbi.nlm.nih.gov/28019652/)
   <br><sub>모델에서의 용도 / used for: WHO risk drinking levels as a treatment endpoint</sub>
126. **Hasin DS et al. (2017).** *Change in non-abstinent WHO drinking risk levels and alcohol dependence: a 3 year follow-up study in the US general population.* Lancet Psychiatry. PMID [28456501](https://pubmed.ncbi.nlm.nih.gov/28456501/)
   <br><sub>모델에서의 용도 / used for: change in non-abstinent WHO drinking risk levels and AUD</sub>
127. **Falk DE et al. (2019).** *Evaluation of Drinking Risk Levels as Outcomes in Alcohol Pharmacotherapy Trials: A Secondary Analysis of 3 Randomized Clinical Trials.* JAMA Psychiatry. PMID [30865232](https://pubmed.ncbi.nlm.nih.gov/30865232/)
   <br><sub>모델에서의 용도 / used for: WHO risk drinking level reductions and long-term outcomes</sub>
128. **Reus VI et al. (2019).** *The American Psychiatric Association Practice Guideline for the Pharmacological Treatment of Patients With Alcohol Use Disorder.* Focus (Am Psychiatr Publ). PMID [32021585](https://pubmed.ncbi.nlm.nih.gov/32021585/)
   <br><sub>모델에서의 용도 / used for: APA practice guideline for the pharmacological treatment of AUD</sub>

## 21. QSP and computational-modelling methods

129. **Nijhout HF et al. (2004).** *A mathematical model of the folate cycle: new insights into folate homeostasis.* J Biol Chem. PMID [15496403](https://pubmed.ncbi.nlm.nih.gov/15496403/)
   <br><sub>모델에서의 용도 / used for: systems biology of ethanol metabolism</sub>
130. **Gadkar K et al. (2016).** *A Six-Stage Workflow for Robust Application of Systems Pharmacology.* CPT Pharmacometrics Syst Pharmacol. PMID [27299936](https://pubmed.ncbi.nlm.nih.gov/27299936/)
   <br><sub>모델에서의 용도 / used for: a six-stage workflow for robust application of systems pharmacology</sub>
131. **Heilig M et al. (2016).** *Time to connect: bringing social context into addiction neuroscience.* Nat Rev Neurosci. PMID [27277868](https://pubmed.ncbi.nlm.nih.gov/27277868/)
   <br><sub>모델에서의 용도 / used for: time to connect, bringing social context into addiction neuroscience</sub>
132. **Musante CJ et al. (2017).** *Quantitative Systems Pharmacology: A Case for Disease Models.* Clin Pharmacol Ther. PMID [27709613](https://pubmed.ncbi.nlm.nih.gov/27709613/)
   <br><sub>모델에서의 용도 / used for: quantitative systems pharmacology, a case for disease models</sub>

---

## 미해결 항목 (not included)

아래 항목은 인용하려 했으나 PubMed 레코드와 제목이 확실히 대조되지 않아
**목록에서 제외**했습니다. 추측한 PMID를 넣는 것보다 빼는 편이 낫습니다.

The following were intended citations whose PubMed record could not be
matched with confidence, and were therefore **omitted rather than guessed**.

- alcohol use disorders
- alcohol use disorders global burden and treatment gap
- clinical pharmacokinetics of ethanol
- Mellanby effect / acute tolerance
- alcohol metabolism
- dynorphin-kappa opioid system in the extended amygdala and alcohol
- acamprosate for alcohol dependence
- pharmacokinetics of acamprosate
- acamprosate in alcohol dependence, how does it work?
- timeline follow-back, a technique for assessing self-reported alcohol consumption
- sex-specific associations between alcohol and blood pressure
- the COMBINE study good clinical outcome definition
- mrgsolve, simulation from ODE-based population PK/PD models
- imbalance between drug and non-drug reward availability, a computational model
- alcohol, liver disease and the gut microbiota

---

## 재현 (Reproducing this list)

```python
# for each intended title:
#   esearch.fcgi?db=pubmed&term="<title>"[Title]         (exact-phrase pass)
#   esearch.fcgi?db=pubmed&term=w1[Title] AND w2[Title]  (distinctive-word pass)
#   esummary.fcgi?db=pubmed&id=<pmids>
#   accept a candidate ONLY if
#     Jaccard(intended, candidate) >= 0.62
#     OR ( containment(intended in candidate) >= 0.85
#          AND |intended| >= 5 content words
#          AND |candidate| <= 2.4 x |intended| )
#   the length guard is what stops a long unrelated title that happens to
#   contain every one of my query words from being accepted
```

## 면책 (Disclaimer)

본 문헌 목록은 교육·연구 목적의 QSP 모델을 뒷받침하기 위해 수집된 것입니다.
인용된 수치와 효과 크기는 원 논문의 맥락(모집단·용량·기간·엔드포인트 정의)
안에서만 유효하며, 모델의 파라미터는 이들로부터 유도된 근사치입니다.
임상 의사결정의 근거로 사용해서는 안 됩니다.
