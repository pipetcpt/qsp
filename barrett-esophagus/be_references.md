# Barrett's Esophagus (BE) — Annotated References

References supporting the QSP model in this directory. Each entry notes **what
the paper contributes to the model** — a parameter, a structural assumption, or
a calibration target. PubMed / PMC links throughout.

---

## 1. Epidemiology & natural history — the progression rates the model must hit

1. **Hvid-Jensen F, et al. Incidence of adenocarcinoma among patients with Barrett's esophagus. N Engl J Med. 2011;365(15):1375-83.**
   <https://pubmed.ncbi.nlm.nih.gov/21995385/>
   Danish nationwide cohort, 11,028 BE patients. EAC incidence **1.2 per 1000
   person-years (0.12%/yr)** — an order of magnitude below the 0.5%/yr figure
   assumed by earlier surveillance guidelines. This is the primary calibration
   target, and it is matched against the model's **low-dose-PPI** NDBE arm
   rather than its untreated arm, because the cohort this rate comes from was
   largely PPI-treated. The model reproduces it at 0.120 %/yr.

2. **Desai TK, et al. The incidence of oesophageal adenocarcinoma in non-dysplastic Barrett's oesophagus: a meta-analysis. Gut. 2012;61(7):970-6.**
   <https://pubmed.ncbi.nlm.nih.gov/21997553/>
   Meta-analytic NDBE → EAC rate **0.33%/yr**; brackets the upper end of the
   model's NDBE hazard.

3. **Duits LC, et al. Barrett's oesophagus patients with low-grade dysplasia can be accurately risk-stratified after histological review by an expert pathology panel. Gut. 2015;64(5):700-6.**
   <https://pubmed.ncbi.nlm.nih.gov/25034523/>
   **73% of community-diagnosed LGD is downstaged** on expert review; confirmed
   LGD progresses to HGD/EAC at **~9%/yr**. Sets `INIT_FLGD`/`INIT_FP53` for
   the confirmed-LGD archetype and motivates the `KREG_LGD` regression term.

4. **Rastogi A, et al. Incidence of esophageal adenocarcinoma in patients with Barrett's esophagus and high-grade dysplasia: a meta-analysis. Gastrointest Endosc. 2008;67(3):394-8.**
   <https://pubmed.ncbi.nlm.nih.gov/18045592/>
   HGD → EAC incidence **~6%/yr**; the calibration target for `K_EAC`.

5. **Anaparthy R, et al. Association between length of Barrett's esophagus and risk of high-grade dysplasia or adenocarcinoma in patients without dysplasia. Clin Gastroenterol Hepatol. 2013;11(11):1430-6.**
   <https://pubmed.ncbi.nlm.nih.gov/23707463/>
   Risk rises with segment length (**~1.7× per additional 2 cm**). Basis for
   `SEG_POW = 1.0` about `SEG_REF = 3 cm`.

6. **Pohl H, et al. Length of Barrett's oesophagus and cancer risk: implications from a large sample of patients with early oesophageal adenocarcinoma. Gut. 2016;65(2):196-201.**
   <https://pubmed.ncbi.nlm.nih.gov/25731871/>
   Independent confirmation of the length–risk gradient, including short-segment
   BE cancers.

7. **Sharma P, et al. The development and validation of an endoscopic grading system for Barrett's esophagus: the Prague C & M criteria. Gastroenterology. 2006;131(5):1392-9.**
   <https://pubmed.ncbi.nlm.nih.gov/17101315/>
   Defines the `Prague_M_cm` state variable and the SSBE/LSBE 3-cm threshold.

8. **Ronkainen J, et al. Prevalence of Barrett's esophagus in the general population: an endoscopic study. Gastroenterology. 2005;129(6):1825-31.**
   <https://pubmed.ncbi.nlm.nih.gov/16344051/>
   Kalixanda population endoscopy study: BE prevalence ~1.6%, most undiagnosed —
   context for the case-finding arm of the map (Cytosponge).

9. **Coleman HG, et al. The epidemiology of esophageal adenocarcinoma. Gastroenterology. 2018;154(2):390-405.**
   <https://pubmed.ncbi.nlm.nih.gov/28780073/>
   Risk-factor magnitudes used to set `SMOKE`, `OBESITY`, `SEX_M` effects.

---

## 2. Reflux physiology, acid and bile exposure

10. **Kahrilas PJ, et al. Approaches to the diagnosis and grading of hiatal hernia. Best Pract Res Clin Gastroenterol. 2008;22(4):601-16.**
    <https://pubmed.ncbi.nlm.nih.gov/18656819/>
    Hiatal hernia axial length as a determinant of the antireflux barrier —
    the `HH_SIZE` term.

11. **Roman S, et al. Ambulatory reflux monitoring for diagnosis of gastro-esophageal reflux disease: update of the Porto consensus. Neurogastroenterol Motil. 2017;29(10):1-15.**
    <https://pubmed.ncbi.nlm.nih.gov/28370768/>
    Acid exposure time (**AET, % of 24 h with pH<4**) as the canonical metric;
    normal <4%, pathologic >6%. Defines the `ACID` compartment's units.

12. **Nehra D, et al. Toxic bile acids in gastro-oesophageal reflux disease: influence of gastric acidity. Gut. 1999;44(5):598-602.**
    <https://pubmed.ncbi.nlm.nih.gov/10205192/>
    Aspirate study: bile-acid concentrations rise with disease severity, and
    **the toxicity of the refluxate is pH-dependent** — direct empirical basis
    for the model's `K_PH_BILE` term.

13. **Kauer WK, et al. Mixed reflux of gastric and duodenal juices is more harmful to the esophagus than gastric juice alone. Ann Surg. 1995;222(4):525-33.**
    <https://pubmed.ncbi.nlm.nih.gov/7574932/>
    Duodenogastroesophageal reflux is more injurious than acid alone — the
    justification for giving bile a large share of reference injury (22%, vs
    70% acid and 8% smoking) *and* for letting that share grow under acid
    suppression rather than shrink.

14. **Vaezi MF, Richter JE. Role of acid and duodenogastroesophageal reflux in gastroesophageal reflux disease. Gastroenterology. 1996;111(5):1192-9.**
    <https://pubmed.ncbi.nlm.nih.gov/8898632/>
    Bilitec bile-exposure measurements; the `BILE` compartment's units.

15. **Tack J. Review article: the role of bile and pepsin in the pathophysiology and treatment of gastro-oesophageal reflux disease. Aliment Pharmacol Ther. 2006;24 Suppl 2:10-6.**
    <https://pubmed.ncbi.nlm.nih.gov/16939429/>
    Pepsin is inactivated above pH 4 while bile-acid delivery is unchanged by
    acid suppression — the asymmetry the model encodes.

16. **Farré R, et al. Critical role of stress in increased oesophageal mucosa permeability and dilated intercellular spaces. Gut. 2007;56(9):1191-7.**
    <https://pubmed.ncbi.nlm.nih.gov/17272649/>
    Dilated intercellular spaces / permeability as the barrier-failure node
    feeding back onto symptom generation.

---

## 3. Bile acids, NF-κB and CDX2 — the metaplasia switch

17. **Kazumori H, et al. Bile acids directly augment caudal related homeobox gene Cdx2 expression in oesophageal keratinocytes in Barrett's epithelium. Gut. 2006;55(1):16-25.**
    <https://pubmed.ncbi.nlm.nih.gov/16118348/>
    Bile acids directly induce **CDX2** in esophageal keratinocytes. The
    structural basis of the `K_BILE_CDX2` edge.

18. **Huo X, et al. Acid and bile salt-induced CDX2 expression differs in esophageal squamous cells from patients with and without Barrett's esophagus. Gastroenterology. 2010;139(1):194-203.**
    <https://pubmed.ncbi.nlm.nih.gov/20303354/>
    CDX2 induction requires NF-κB activation and differs by host susceptibility
    — supports the `K_INF_CDX2` term alongside the bile term.

19. **Wang DH. The esophageal squamous epithelial cell — still a reasonable candidate for the Barrett's esophagus cell of origin? Cell Mol Gastroenterol Hepatol. 2017;4(1):157-160.**
    <https://pubmed.ncbi.nlm.nih.gov/28848795/>
    Cell-of-origin debate underlying the `Progenitor_origin` node.

20. **Que J, et al. Pathogenesis and cells of origin of Barrett's esophagus. Gastroenterology. 2019;157(2):349-364.**
    <https://pubmed.ncbi.nlm.nih.gov/31082367/>
    Comprehensive review of BMP4/SOX9/CDX2/HNF4A reprogramming — the cluster-6
    architecture of the map.

21. **Kong J, et al. Induction of intestinalization in human esophageal keratinocytes is a multistep process. Carcinogenesis. 2011;32(6):922-30.**
    <https://pubmed.ncbi.nlm.nih.gov/21444358/>
    CDX2 alone is insufficient; multistep reprogramming justifies the sigmoidal
    (Hill) `CDX2` induction rather than a linear one.

22. **Milano F, et al. Bone morphogenetic protein 4 expressed in esophagitis induces a columnar phenotype in esophageal squamous cells. Gastroenterology. 2007;133(4):1198-209.**
    <https://pubmed.ncbi.nlm.nih.gov/17919494/>
    BMP4 → columnar phenotype; the BMP4 → SOX9 → CDX2 chain in the map.

---

## 4. Oxidative stress, inflammation and proliferation

23. **Jenkins GJ, et al. Deoxycholic acid at neutral and acid pH is genotoxic to oesophageal cells through the induction of ROS: the potential role of anti-oxidants in Barrett's oesophagus. Carcinogenesis. 2007;28(1):136-42.**
    <https://pubmed.ncbi.nlm.nih.gov/16966446/>
    **Deoxycholate is genotoxic at neutral pH via ROS** — the single most
    important citation for the model's central tension (`K_PH_BILE`).

24. **Dvorak K, et al. Bile acids in combination with low pH induce oxidative stress and oxidative DNA damage: relevance to the pathogenesis of Barrett's oesophagus. Gut. 2007;56(6):763-71.**
    <https://pubmed.ncbi.nlm.nih.gov/17145738/>
    8-OHdG oxidative DNA adducts as the injury readout — the `INJURY`
    compartment's biological referent.

25. **Fitzgerald RC, et al. Inflammatory gradient in Barrett's oesophagus: implications for disease complications. Gut. 2002;51(3):316-22.**
    <https://pubmed.ncbi.nlm.nih.gov/12171950/>
    IL-8/IL-1β gradient along the segment; supports the inflammation cluster.

26. **Souza RF, et al. Gastroesophageal reflux might cause esophagitis through a cytokine-mediated mechanism rather than caustic acid injury. Gastroenterology. 2009;137(5):1776-84.**
    <https://pubmed.ncbi.nlm.nih.gov/19660463/>
    Reflux esophagitis is cytokine-mediated, not simply a chemical burn — why
    `INFLAM` is an explicit state between exposure and metaplasia.

27. **Konturek PC, et al. Prostaglandins as mediators of COX-2 derived carcinogenesis in gastrointestinal tract. J Physiol Pharmacol. 2005;56 Suppl 5:57-73.**
    <https://pubmed.ncbi.nlm.nih.gov/16247188/>
    COX-2 → PGE2 → EP2/EP4 → cyclin D1 proliferative axis, the target of the
    aspirin arm.

---

## 5. Clonal evolution & genomics — the ordering of the hits

28. **Weaver JMJ, et al. Ordering of mutations in preinvasive disease stages of esophageal carcinogenesis. Nat Genet. 2014;46(8):837-843.**
    <https://pubmed.ncbi.nlm.nih.gov/24952744/>
    **TP53 mutation marks the transition to HGD/EAC while SMAD4 is
    EAC-restricted** — the ordering that the cascade `FP16 → FLGD → FP53 →
    FHGD → EAC` implements.

29. **Stachler MD, et al. Paired exome analysis of Barrett's esophagus and adenocarcinoma. Nat Genet. 2015;47(9):1047-55.**
    <https://pubmed.ncbi.nlm.nih.gov/26192918/>
    Genome doubling after TP53 loss; punctuated rather than gradual evolution —
    the `ANEUPL` state and its `K_ANEU_EAC` amplification of the EAC transition.

30. **Ross-Innes CS, et al. Whole-genome sequencing provides new insights into the clonal architecture of Barrett's esophagus and esophageal adenocarcinoma. Nat Genet. 2015;47(9):1038-46.**
    <https://pubmed.ncbi.nlm.nih.gov/26192915/>
    Clonal diversity of the metaplastic field; why risk scales with segment
    area, not just grade.

31. **Wong DJ, et al. p16(INK4a) lesions are common, early abnormalities that undergo clonal expansion in Barrett's metaplastic epithelium. Cancer Res. 2001;61(22):8284-9.**
    <https://pubmed.ncbi.nlm.nih.gov/11719461/>
    **CDKN2A/p16 loss is early and common**, forming large clonal fields —
    justifies the high default `INIT_FP16 = 0.45`.

32. **Reid BJ, et al. Predictors of progression to cancer in Barrett's esophagus: baseline histology and flow cytometry identify low- and high-risk patient subsets. Am J Gastroenterol. 2000;95(7):1669-76.**
    <https://pubmed.ncbi.nlm.nih.gov/10925966/>
    Aneuploidy/tetraploidy by flow cytometry predicts progression — the
    `Aneuploidy` risk multiplier.

33. **Kastelein F, et al. Aberrant p53 protein expression is associated with an increased risk of neoplastic progression in patients with Barrett's oesophagus. Gut. 2013;62(12):1676-83.**
    <https://pubmed.ncbi.nlm.nih.gov/23256952/>
    Aberrant p53 IHC as a clinically usable risk stratifier — the `p53_IHC`
    node and the `TP53_mut_frac` output.

34. **Killcoyne S, et al. Genomic copy number predicts esophageal cancer years before transformation. Nat Med. 2020;26(11):1726-1732.**
    <https://pubmed.ncbi.nlm.nih.gov/33020648/>
    Copy-number burden predicts progression years in advance — supports
    modeling genomic instability as a slow accumulating state.

---

## 6. Acid-suppression pharmacology (PPI, P-CAB) and its limits

35. **Shin JM, Sachs G. Pharmacology of proton pump inhibitors. Curr Gastroenterol Rep. 2008;10(6):528-34.**
    <https://pubmed.ncbi.nlm.nih.gov/19006606/>
    Covalent sulfenamide binding and **pump resynthesis half-life ~50 h** —
    the `KSYN_PUMP = 0.333/day` and `KINH_PUMP` structure.

36. **Klotz U, et al. Clinical impact of CYP2C19 polymorphism on the action of proton pump inhibitors: a review of a special problem. Int J Clin Pharmacol Ther. 2006;44(7):297-302.**
    <https://pubmed.ncbi.nlm.nih.gov/16961157/>
    CYP2C19 phenotype changes PPI exposure several-fold — basis for `CL2C19`
    (UM 1.8 … PM 0.25).

37. **Furuta T, et al. Influence of CYP2C19 pharmacogenetic polymorphism on proton pump inhibitor-based therapies. Drug Metab Pharmacokinet. 2005;20(3):153-67.**
    <https://pubmed.ncbi.nlm.nih.gov/15988117/>
    Genotype-dependent acid control and treatment failure in extensive/ultrarapid
    metabolizers.

38. **Lima JJ, et al. Clinical Pharmacogenetics Implementation Consortium (CPIC) guideline for CYP2C19 and proton pump inhibitor dosing. Clin Pharmacol Ther. 2021;109(6):1417-1423.**
    <https://pubmed.ncbi.nlm.nih.gov/32770672/>
    Actionable phenotype definitions used for the pharmacogenomic scenario.

39. **Sachs G, et al. The pharmacology of the gastric acid pump: the H+,K+ ATPase. Annu Rev Pharmacol Toxicol. 1995;35:277-305.**
    <https://pubmed.ncbi.nlm.nih.gov/7598495/>
    Pump biology underlying both PPI (covalent) and P-CAB (K+-competitive)
    mechanisms.

40. **Oshima T, Miwa H. Potent potassium-competitive acid blockers: a new era for the treatment of acid-related diseases. J Neurogastroenterol Motil. 2018;24(3):334-344.**
    <https://pubmed.ncbi.nlm.nih.gov/29739175/>
    Vonoprazan pharmacology: t½ ~7 h, rapid and sustained acid inhibition,
    minimal CYP2C19 dependence — the `VPZ` arm.

41. **Lundell L, et al. Systematic review: the effects of long-term proton pump inhibitor use on serum gastrin levels and gastric histology. Aliment Pharmacol Ther. 2015;42(6):649-63.**
    <https://pubmed.ncbi.nlm.nih.gov/26177572/>
    Long-term PPI raises gastrin severalfold — the `K_GASTRIN` calibration.

42. **Haigh CR, et al. Gastrin induces proliferation in Barrett's metaplasia through activation of the CCK2 receptor. Gastroenterology. 2003;124(3):615-25.**
    <https://pubmed.ncbi.nlm.nih.gov/12612900/>
    **Gastrin is trophic to Barrett epithelium via CCK2R** — the empirical
    warrant for `K_GAS_KI`, the second half of the model's central tension.

43. **Freedberg DE, et al. The risks and benefits of long-term use of proton pump inhibitors: expert review and best practice advice from the AGA. Gastroenterology. 2017;152(4):706-715.**
    <https://pubmed.ncbi.nlm.nih.gov/28257716/>
    Long-term PPI safety signals summarized in the `PPI_safety_index` readout.

---

## 7. Chemoprevention trials — the model's hazard-ratio targets

44. **Jankowski JAZ, et al. Esomeprazole and aspirin in Barrett's oesophagus (AspECT): a randomised factorial trial. Lancet. 2018;392(10145):400-408.**
    <https://pubmed.ncbi.nlm.nih.gov/30060998/>
    **The single most important calibration target.** 2557 patients, 2×2
    factorial, median 8.9 y. Composite (all-cause mortality, EAC, HGD): high-dose
    esomeprazole **HR 0.73**, aspirin **HR 0.93 (NS alone)**, high-dose PPI +
    aspirin vs low-dose PPI alone **HR 0.59**.

45. **Kastelein F, et al. Proton pump inhibitors reduce the risk of neoplastic progression in patients with Barrett's esophagus. Clin Gastroenterol Hepatol. 2013;11(4):382-8.**
    <https://pubmed.ncbi.nlm.nih.gov/23200977/>
    Observational cohort: PPI use associated with reduced progression (HR ~0.4),
    an upper bound on the plausible PPI effect.

46. **Singh S, et al. Statins are associated with reduced risk of esophageal cancer, particularly in patients with Barrett's esophagus: a systematic review and meta-analysis. Clin Gastroenterol Hepatol. 2013;11(6):620-9.**
    <https://pubmed.ncbi.nlm.nih.gov/23357487/>
    Statin association (~30-40% reduction) — the adjunct node in the map.

47. **Corley DA, et al. Protective association of aspirin/NSAIDs and esophageal cancer: a systematic review and meta-analysis. Gastroenterology. 2003;124(1):47-56.**
    <https://pubmed.ncbi.nlm.nih.gov/12512029/>
    Pre-trial observational basis for the aspirin arm.

48. **Chak A, et al. Metformin does not reduce markers of cell proliferation in esophageal tissues of patients with Barrett's esophagus. Clin Gastroenterol Hepatol. 2015;13(4):665-72.**
    <https://pubmed.ncbi.nlm.nih.gov/25218668/>
    Negative phase 2 (no pS6K reduction) — why metformin appears in the map as
    a dashed, weak edge rather than an active treatment arm.

49. **Heath EI, et al. Secondary chemoprevention of Barrett's esophagus with celecoxib: results of a randomized trial. J Natl Cancer Inst. 2007;99(7):545-57.**
    <https://pubmed.ncbi.nlm.nih.gov/17405999/>
    Celecoxib did **not** reduce dysplasia progression — an important negative
    result constraining how strong the PGE2 arm alone can be.

---

## 8. Endoscopic eradication therapy

50. **Shaheen NJ, et al. Radiofrequency ablation in Barrett's esophagus with dysplasia (AIM Dysplasia trial). N Engl J Med. 2009;360(22):2277-88.**
    <https://pubmed.ncbi.nlm.nih.gov/19474425/>
    RFA: **CE-D 90.5%, CE-IM 77.4%**, progression 3.6% vs 16.3%, stricture 6%.
    Calibrates `K_ABL_DYS` and `K_STRICT`.

51. **Phoa KN, et al. Radiofrequency ablation vs endoscopic surveillance for patients with Barrett esophagus and low-grade dysplasia (SURF): a randomized clinical trial. JAMA. 2014;311(12):1209-17.**
    <https://pubmed.ncbi.nlm.nih.gov/24668102/>
    **Progression to HGD/EAC 1.5% (RFA) vs 26.5% (surveillance) at 3 years.**
    The target for scenarios 7 and 8.

52. **Pech O, et al. Long-term efficacy and safety of endoscopic resection for patients with mucosal adenocarcinoma of the esophagus. Gastroenterology. 2014;146(3):652-60.**
    <https://pubmed.ncbi.nlm.nih.gov/24269290/>
    EMR for T1a EAC: long-term complete remission ~93% — the `EMR_PULSE` arm.

53. **Cotton CC, et al. Development of evidence-based surveillance intervals after radiofrequency ablation of Barrett's esophagus. Gastroenterology. 2018;155(2):316-326.**
    <https://pubmed.ncbi.nlm.nih.gov/29674001/>
    Post-CE-IM recurrence kinetics (~8-10%/yr early) — the behavior the model
    generates from persistent reflux rather than from a fitted recurrence rate.

54. **Krishnamoorthi R, et al. Factors associated with progression of Barrett's esophagus: a systematic review and meta-analysis. Clin Gastroenterol Hepatol. 2018;16(7):1046-1055.**
    <https://pubmed.ncbi.nlm.nih.gov/29199144/>
    Pooled risk factors for progression, used to sanity-check relative effect
    sizes among model covariates.

55. **Haidry RJ, et al. Improvement over time in outcomes for patients undergoing endoscopic therapy for Barrett's oesophagus-related neoplasia: 6-year experience from the UK RFA registry. Gut. 2015;64(8):1192-9.**
    <https://pubmed.ncbi.nlm.nih.gov/25539672/>
    Real-world CE-D/CE-IM rates and the number of sessions needed (median 1
    EMR + 2-3 RFA) — the 3-4 session `rfa()` courses in the scenarios.

56. **Canto MI, et al. Multifocal cryoballoon ablation for eradication of Barrett's esophagus-related neoplasia. Gastrointest Endosc. 2020;92(6):1200-1212.**
    <https://pubmed.ncbi.nlm.nih.gov/32544553/>
    Cryoablation as an alternative modality (map node `Cryo`).

---

## 9. Antireflux surgery and non-drug reflux control

57. **Maret-Ouda J, et al. Antireflux surgery and risk of esophageal adenocarcinoma: a systematic review and meta-analysis. Ann Surg. 2016;263(2):251-7.**
    <https://pubmed.ncbi.nlm.nih.gov/26501714/>
    Fundoplication reduces but does not abolish EAC risk — the calibration
    ceiling for `EFF_FUNDO_A/B`.

58. **Zaninotto G, et al. Long-term results of Nissen fundoplication in reflux disease. Br J Surg. 2005 / Attwood SE, et al. Long-term safety of proton pump inhibitor therapy assessed under controlled, randomised clinical trial conditions: data from the SOPRAN and LOTUS studies. Aliment Pharmacol Ther. 2015;41(11):1162-74.**
    <https://pubmed.ncbi.nlm.nih.gov/25858591/>
    Head-to-head long-term PPI vs surgery data; supports surgery's distinct
    effect on **volume/bile** reflux, the mechanistic point of scenario 10.

59. **Zhang Q, et al. Effect of baclofen on transient lower esophageal sphincter relaxations. Gastroenterology / Aliment Pharmacol Ther reviews of GABA-B agonists.**
    <https://pubmed.ncbi.nlm.nih.gov/12105842/>
    Baclofen reduces TLESRs and reflux episodes — the `BACLOFEN` parameter.

60. **Leiman DA, et al. Alginate therapy is effective treatment for GERD symptoms: a systematic review and meta-analysis. Dis Esophagus. 2017;30(5):1-9.**
    <https://pubmed.ncbi.nlm.nih.gov/28375448/>
    Alginate raft postprandial effect (`ALGINATE`).

61. **Banerjee B, et al. Effect of ursodeoxycholic acid on oxidative DNA damage in Barrett's esophagus. / Peng S, et al. In Barrett's esophagus patients and Barrett's cell lines, ursodeoxycholic acid increases antioxidant expression and prevents DNA damage by bile acids. Am J Physiol Gastrointest Liver Physiol. 2014;307(2):G129-39.**
    <https://pubmed.ncbi.nlm.nih.gov/24852565/>
    UDCA shifts the bile pool and reduces bile-induced DNA damage — the
    `EMAX_UDCA` arm.

---

## 10. Obesity, metabolic drivers

62. **Kendall BJ, et al. Leptin and the risk of Barrett's oesophagus. Gut. 2008;57(4):448-54.**
    <https://pubmed.ncbi.nlm.nih.gov/18178609/>
    Leptin associated with BE risk independent of reflux — the `K_ADI_KI` term.

63. **Ryan AM, et al. Obesity, metabolic syndrome and esophageal adenocarcinoma: epidemiology, etiology and new targets. Cancer Epidemiol. 2011;35(4):309-19.**
    <https://pubmed.ncbi.nlm.nih.gov/21470937/>
    Central obesity acts through both mechanical reflux and adipokine/IGF
    signaling — the two distinct obesity edges in the model
    (`K_OBES_REFL`, `K_ADI_INF`/`K_ADI_KI`).

64. **Singh S, et al. Central adiposity is associated with increased risk of esophageal inflammation, metaplasia, and adenocarcinoma: a systematic review and meta-analysis. Clin Gastroenterol Hepatol. 2013;11(11):1399-1412.**
    <https://pubmed.ncbi.nlm.nih.gov/23707461/>
    Dose-response for visceral adiposity across the whole cascade.

---

## 11. Surveillance, biomarkers and guidelines

65. **Shaheen NJ, et al. Diagnosis and management of Barrett's esophagus: an updated ACG guideline. Am J Gastroenterol. 2022;117(4):559-587.**
    <https://pubmed.ncbi.nlm.nih.gov/35354777/>
    Current surveillance intervals, LGD management, endoscopic eradication
    indications — the clinical scaffolding of scenarios 7-9.

66. **Fitzgerald RC, et al. British Society of Gastroenterology guidelines on the diagnosis and management of Barrett's oesophagus. Gut. 2014;63(1):7-42.**
    <https://pubmed.ncbi.nlm.nih.gov/24165758/>
    Segment-length-stratified surveillance and p53 IHC recommendation.

67. **Fitzgerald RC, et al. Cytosponge-trefoil factor 3 versus usual care to identify Barrett's oesophagus in a primary care setting (BEST3): a multicentre, pragmatic, randomised controlled trial. Lancet. 2020;396(10247):333-344.**
    <https://pubmed.ncbi.nlm.nih.gov/32738958/>
    Non-endoscopic case finding — the `Cytosponge` node.

68. **Vennalaganti PR, et al. Increased detection of Barrett's esophagus-associated neoplasia using wide-area transepithelial sampling: a multicenter, prospective, randomized trial. Gastrointest Endosc. 2018;87(2):348-355.**
    <https://pubmed.ncbi.nlm.nih.gov/28757316/>
    WATS-3D reduces sampling error — the `Sampling_error` node.

69. **Critchley-Thorne RJ, et al. A tissue systems pathology assay for high-risk Barrett's esophagus. Cancer Epidemiol Biomarkers Prev. 2016;25(6):958-68.**
    <https://pubmed.ncbi.nlm.nih.gov/27197290/>
    Objective multi-marker risk prediction (`TissueCypher`).

70. **Sharma P, et al. AGA Clinical Practice Update on Endoscopic Treatment of Barrett's Esophagus with Dysplasia and/or Early Cancer: Expert Review. Gastroenterology. 2020;158(3):760-769.**
    <https://pubmed.ncbi.nlm.nih.gov/31730766/>
    Best-practice advice for EET, post-ablation surveillance and recurrence
    management.

---

## 12. Modeling methodology

71. **Baio FE, et al. / Kong CY, et al. Exploring the recent trend in esophageal adenocarcinoma incidence and mortality using comparative simulation modeling. Cancer Epidemiol Biomarkers Prev. 2014;23(6):997-1006.**
    <https://pubmed.ncbi.nlm.nih.gov/24692500/>
    CISNET comparative microsimulation of the BE → EAC pathway; the
    population-level counterpart to this mechanistic model.

72. **Curtius K, et al. A molecular clock infers heterogeneous tissue age among patients with Barrett's esophagus. PLoS Comput Biol. 2016;12(5):e1004919.**
    <https://pubmed.ncbi.nlm.nih.gov/27168458/>
    Quantitative clonal-dynamics modeling of the Barrett field — methodological
    basis for representing dysplasia as clone fractions rather than discrete
    states.

73. **Hazelton WD, et al. / Kroep S, et al. An accurate cancer incidence in Barrett's esophagus: a best estimate using published data and modeling. Gastroenterology. 2015;149(3):577-585.**
    <https://pubmed.ncbi.nlm.nih.gov/26010929/>
    Model-based reconciliation of divergent published EAC incidence estimates —
    the reason the model targets ~0.1-0.3%/yr rather than 0.5%/yr.

74. **Baker CD, et al. mrgsolve: Simulate from ODE-Based Population PK/PD and Systems Pharmacology Models.**
    <https://mrgsolve.org/>
    Simulation engine used by `be_mrgsolve_model.R`.

---

## How the model's falsifiable claims map onto this list

| Model element | Would be falsified by |
|---|---|
| `K_PH_BILE = 0.35` (effective bile toxicity rises with acid blockade, saturating in pH) | Showing bile-driven CDX2/ROS signaling is *reduced*, not increased, at neutral pH (contra refs 12, 23) |
| `K_GAS_KI = 4.0` (hypergastrinemia is trophic: a doubling of gastrin raises Ki-67 by 20%) | Showing CCK2R blockade does not change Barrett proliferation on PPI (contra ref 42) |
| High-dose PPI cumulative-composite ratio 0.73, not ~0.1 | An RCT showing PPI monotherapy abolishes progression (contra ref 44) |
| Aspirin reaches only the COX-dependent (PGE2) fraction of proliferation, so its effect must be small | Celecoxib's negative trial already limits this arm (ref 49); a null aspirin effect *on top of* high-dose PPI would falsify the additivity |
| Post-ablation recurrence is a reflux-control phenomenon | Showing recurrence rates are independent of post-ablation acid/bile exposure (contra refs 53, 57) |

---

*Compiled for the QSP Disease Model Library. Educational and research use only.*
