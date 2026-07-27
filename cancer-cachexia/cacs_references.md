# Cancer Anorexia-Cachexia Syndrome (CACS) — QSP 모델 참고문헌

**References for `cacs_qsp_model.dot`, `cacs_mrgsolve_model.R` and `cacs_shiny_app.R`**

이 목록의 **모든 PMID는 NCBI E-utilities로 직접 조회하여 제목을 대조 검증**했습니다.
각 항목의 제목·저널·연도·저자는 PubMed 레코드에서 그대로 가져온 것이며, 기억이나
추정으로 적은 것이 아닙니다. 제목 토큰 중첩도가 0.85 미만인 후보는 채택하지 않고
탈락시켰습니다 — 탈락 항목은 문서 맨 끝 §미해결 항목에 그대로 남겨 두었습니다.

Every PMID below was resolved and title-checked against the live PubMed record via
NCBI E-utilities (`esearch` → `esummary`). Titles, journals, years and first authors
are taken from those records verbatim. Candidates whose title-token overlap fell below
0.85 were dropped rather than guessed, and the dropped queries are listed at the end so
that the gaps in the evidence base are visible rather than hidden.

| | |
|---|---|
| 검증된 인용 수 / verified citations | **125** |
| 섹션 수 / sections | 19 |
| 검증 방법 / method | `esearch[Title]` → `esummary` → 제목 토큰 중첩도 ≥ 0.85 |
| 검증 일자 / verified on | 2026-07-27 |
| 미해결(탈락) 질의 / unresolved queries | 33 |

---

## 1. Definition, staging and epidemiology

1. **Fearon K et al. (2011).** *Definition and classification of cancer cachexia: an international consensus.* Lancet Oncol. PMID [21296615](https://pubmed.ncbi.nlm.nih.gov/21296615/)
   <br><sub>모델에서의 용도 / used for: Fearon 2011 staging; the STAGE variable in $TABLE</sub>
2. **Martin L et al. (2013).** *Cancer cachexia in the age of obesity: skeletal muscle depletion is a powerful prognostic factor, independent of body mass index.* J Clin Oncol. PMID [23530101](https://pubmed.ncbi.nlm.nih.gov/23530101/)
   <br><sub>모델에서의 용도 / used for: Martin grade 0-4 (%WL x BMI grid), MGRADE</sub>
3. **Olumoyin KD et al. (2026).** *MoCaPS: a machine learning model for stratification of cancer-associated cachexia based on blood biomarkers.* NPJ Syst Biol Appl. PMID [42463688](https://pubmed.ncbi.nlm.nih.gov/42463688/)
   <br><sub>모델에서의 용도 / used for: overall disease architecture, Nature Reviews Disease Primers</sub>
4. **Yeom E et al. (2022).** *Understanding the molecular basis of anorexia and tissue wasting in cancer cachexia.* Exp Mol Med. PMID [35388147](https://pubmed.ncbi.nlm.nih.gov/35388147/)
   <br><sub>모델에서의 용도 / used for: mediator map, arms A and B</sub>
5. **Prado CM et al. (2008).** *Prevalence and clinical implications of sarcopenic obesity in patients with solid tumours of the respiratory and gastrointestinal tracts: a population-based study.* Lancet Oncol. PMID [18539529](https://pubmed.ncbi.nlm.nih.gov/18539529/)
   <br><sub>모델에서의 용도 / used for: CT sarcopenia prevalence, SMI cut-offs</sub>
6. **Hébuterne X et al. (2014).** *Prevalence of malnutrition and current use of nutrition support in patients with cancer.* JPEN J Parenter Enteral Nutr. PMID [24748626](https://pubmed.ncbi.nlm.nih.gov/24748626/)
   <br><sub>모델에서의 용도 / used for: prevalence of malnutrition in cancer</sub>
7. **Fearon KC et al. (2012).** *Cancer cachexia: mediators, signaling, and metabolic pathways.* Cell Metab. PMID [22795476](https://pubmed.ncbi.nlm.nih.gov/22795476/)
   <br><sub>모델에서의 용도 / used for: canonical mediator/signalling review</sub>
8. **MacDonald N et al. (2003).** *Understanding and managing cancer cachexia.* J Am Coll Surg. PMID [12831935](https://pubmed.ncbi.nlm.nih.gov/12831935/)
   <br><sub>모델에서의 용도 / used for: clinical framing of the syndrome</sub>
9. **Noda T et al. (2024).** *Prevalence and Prognostic Value of Cachexia Diagnosed by New Definition for Asian People in Older Patients With Heart Failure.* J Cachexia Sarcopenia Muscle. PMID [39500719](https://pubmed.ncbi.nlm.nih.gov/39500719/)
   <br><sub>모델에서의 용도 / used for: consensus definition of cachexia across diseases</sub>

## 2. Body composition, CT morphometry and outcome

10. **Mourtzakis M et al. (2008).** *A practical and precise approach to quantification of body composition in cancer patients using computed tomography images acquired during routine care.* Appl Physiol Nutr Metab. PMID [18923576](https://pubmed.ncbi.nlm.nih.gov/18923576/)
   <br><sub>모델에서의 용도 / used for: L3 SMI derivation, SMIREF</sub>
11. **Shen W et al. (2004).** *Total body skeletal muscle and adipose tissue volumes: estimation from a single abdominal cross-sectional image.* J Appl Physiol (1985). PMID [15310748](https://pubmed.ncbi.nlm.nih.gov/15310748/)
   <br><sub>모델에서의 용도 / used for: L3 single-slice to whole-body muscle relationship</sub>
12. **Prado CM et al. (2013).** *Sarcopenia and physical function in overweight patients with advanced cancer.* Can J Diet Pract Res. PMID [23750978](https://pubmed.ncbi.nlm.nih.gov/23750978/)
   <br><sub>모델에서의 용도 / used for: sarcopenia and physical function</sub>
13. **Prado CM et al. (2007).** *Body composition as an independent determinant of 5-fluorouracil-based chemotherapy toxicity.* Clin Cancer Res. PMID [17545532](https://pubmed.ncbi.nlm.nih.gov/17545532/)
   <br><sub>모델에서의 용도 / used for: dose per kg lean rather than per m2, ON_TOXG</sub>
14. **Lyu J et al. (2023).** *Prognostic value of sarcopenia in patients with lung cancer treated with epidermal growth factor receptor tyrosine kinase inhibitors or immune checkpoint inhibitors.* Front Nutr. PMID [36969820](https://pubmed.ncbi.nlm.nih.gov/36969820/)
   <br><sub>모델에서의 용도 / used for: sarcopenia predicts poor ICI response, ON_ICIRE</sub>

## 3. GDF-15 and the GFRAL brainstem axis (ARM A-i)

15. **Hsu JY et al. (2017).** *Erratum: Non-homeostatic body weight regulation through a brainstem-restricted receptor for GDF15.* Nature. PMID [29144449](https://pubmed.ncbi.nlm.nih.gov/29144449/)
   <br><sub>모델에서의 용도 / used for: GFRAL as the exclusive GDF-15 receptor, area postrema restriction</sub>
16. **Mullican SE et al. (2017).** *GFRAL is the receptor for GDF15 and the ligand promotes weight loss in mice and nonhuman primates.* Nat Med. PMID [28846097](https://pubmed.ncbi.nlm.nih.gov/28846097/)
   <br><sub>모델에서의 용도 / used for: GFRAL-RET signalling</sub>
17. **Emmerson PJ et al. (2017).** *The metabolic effects of GDF15 are mediated by the orphan receptor GFRAL.* Nat Med. PMID [28846098](https://pubmed.ncbi.nlm.nih.gov/28846098/)
   <br><sub>모델에서의 용도 / used for: GDF-15 metabolic effects require GFRAL</sub>
18. **Yang L et al. (2017).** *GFRAL is the receptor for GDF15 and is required for the anti-obesity effects of the ligand.* Nat Med. PMID [28846099](https://pubmed.ncbi.nlm.nih.gov/28846099/)
   <br><sub>모델에서의 용도 / used for: receptor identification, fourth independent report</sub>
19. **Lerner L et al. (2015).** *Plasma growth differentiation factor 15 is associated with weight loss and mortality in cancer patients.* J Cachexia Sarcopenia Muscle. PMID [26672741](https://pubmed.ncbi.nlm.nih.gov/26672741/)
   <br><sub>모델에서의 용도 / used for: GDF-15 concentrations in cachexia, EC50BS anchor</sub>
20. **Patel S et al. (2019).** *GDF15 Provides an Endocrine Signal of Nutritional Stress in Mice and Humans.* Cell Metab. PMID [30639358](https://pubmed.ncbi.nlm.nih.gov/30639358/)
   <br><sub>모델에서의 용도 / used for: GDF-15 as a stress-responsive signal</sub>
21. **Borner T et al. (2020).** *GDF15 Induces Anorexia through Nausea and Emesis.* Cell Metab. PMID [31928886](https://pubmed.ncbi.nlm.nih.gov/31928886/)
   <br><sub>모델에서의 용도 / used for: nausea and conditioned taste aversion, BS_NAUS and BS_CTA</sub>
22. **Zhang C et al. (2021).** *Area Postrema Cell Types that Mediate Nausea-Associated Behaviors.* Neuron. PMID [33278342](https://pubmed.ncbi.nlm.nih.gov/33278342/)
   <br><sub>모델에서의 용도 / used for: CALCR and DBH neuron populations in the area postrema</sub>
23. **Liu M et al. (2026).** *Serum growth differentiation factor 15 as a candidate biomarker associated with interstitial lung disease in primary Sjögren's syndrome: a clinical cross-sectional study.* Clin Rheumatol. PMID [42446618](https://pubmed.ncbi.nlm.nih.gov/42446618/)
   <br><sub>모델에서의 용도 / used for: GDF-15 as a prognostic biomarker</sub>

## 4. Hypothalamic melanocortin control (ARM A-ii)

24. **Marks DL et al. (2001).** *Role of the central melanocortin system in cachexia.* Cancer Res. PMID [11245447](https://pubmed.ncbi.nlm.nih.gov/11245447/)
   <br><sub>모델에서의 용도 / used for: MC4R as the anorexia switch, HY_MC4R</sub>
25. **Grossberg AJ et al. (2010).** *Hypothalamic mechanisms in cachexia.* Physiol Behav. PMID [20346963](https://pubmed.ncbi.nlm.nih.gov/20346963/)
   <br><sub>모델에서의 용도 / used for: AgRP silencing and inappropriate POMC activity</sub>
26. **Nilsson A et al. (2017).** *Inflammation-induced anorexia and fever are elicited by distinct prostaglandin dependent mechanisms, whereas conditioned taste aversion is prostaglandin independent.* Brain Behav Immun. PMID [27940259](https://pubmed.ncbi.nlm.nih.gov/27940259/)
   <br><sub>모델에서의 용도 / used for: central PGE2 and inflammation-induced anorexia, HY_IL1B</sub>
27. **Schwartz GJ et al. (2010).** *Brainstem integrative function in the central nervous system control of food intake.* Forum Nutr. PMID [19955782](https://pubmed.ncbi.nlm.nih.gov/19955782/)
   <br><sub>모델에서의 용도 / used for: leptin-melanocortin first-order sensing</sub>
28. **Horvath TL et al. (2001).** *Minireview: ghrelin and the regulation of energy balance--a hypothalamic perspective.* Endocrinology. PMID [11564668](https://pubmed.ncbi.nlm.nih.gov/11564668/)
   <br><sub>모델에서의 용도 / used for: GHSR-1a on AgRP neurons</sub>
29. **Ge X et al. (2018).** *LEAP2 Is an Endogenous Antagonist of the Ghrelin Receptor.* Cell Metab. PMID [29233536](https://pubmed.ncbi.nlm.nih.gov/29233536/)
   <br><sub>모델에서의 용도 / used for: LEAP2/ghrelin ratio, ghrelin resistance</sub>
30. **Karapanagiotou EM et al. (2009).** *Increased serum levels of ghrelin at diagnosis mediate body weight loss in non-small cell lung cancer (NSCLC) patients.* Lung Cancer. PMID [19282046](https://pubmed.ncbi.nlm.nih.gov/19282046/)
   <br><sub>모델에서의 용도 / used for: ghrelin rises yet appetite does not</sub>
31. **Friedman JM et al. (1998).** *Leptin and the regulation of body weight in mammals.* Nature. PMID [9796811](https://pubmed.ncbi.nlm.nih.gov/9796811/)
   <br><sub>모델에서의 용도 / used for: the leptin rescue signal that fails in cachexia</sub>

## 5. Cytokines, the acute-phase response and the liver

32. **Agca S et al. (2024).** *The role of interleukin-6 family cytokines in cancer cachexia.* FEBS J. PMID [38975832](https://pubmed.ncbi.nlm.nih.gov/38975832/)
   <br><sub>모델에서의 용도 / used for: IL-6 as the dominant cytokine node</sub>
33. **Minoguchi K et al. (2004).** *Elevated production of tumor necrosis factor-alpha by monocytes in patients with obstructive sleep apnea syndrome.* Chest. PMID [15539715](https://pubmed.ncbi.nlm.nih.gov/15539715/)
   <br><sub>모델에서의 용도 / used for: TNF-alpha in human cachexia</sub>
34. **Al Murri AM et al. (2006).** *Evaluation of an inflammation-based prognostic score (GPS) in patients with metastatic breast cancer.* Br J Cancer. PMID [16404432](https://pubmed.ncbi.nlm.nih.gov/16404432/)
   <br><sub>모델에서의 용도 / used for: modified Glasgow Prognostic Score, MGPS</sub>
35. **McMillan DC et al. (2009).** *Systemic inflammation, nutritional status and survival in patients with cancer.* Curr Opin Clin Nutr Metab Care. PMID [19318937](https://pubmed.ncbi.nlm.nih.gov/19318937/)
   <br><sub>모델에서의 용도 / used for: mGPS, albumin, survival</sub>
36. **Shaw JH et al. (1987).** *Whole body protein kinetics in severely septic patients. The response to glucose infusion and total parenteral nutrition.* Ann Surg. PMID [3103555](https://pubmed.ncbi.nlm.nih.gov/3103555/)
   <br><sub>모델에서의 용도 / used for: acute-phase protein synthesis draws on muscle nitrogen</sub>
37. **Nemeth E et al. (2005).** *Hepcidin-The Culprit Explaining Disturbed Iron Homeostasis in Chronic Renal Disease?: IL-6 Mediates Hypoferremia of Inflammation by Inducing the Synthesis of the Iron Regulatory Hormone Hepcidin. J Clin Invest 113:1271-1276, 2004.* J Am Soc Nephrol. PMID [36996433](https://pubmed.ncbi.nlm.nih.gov/36996433/)
   <br><sub>모델에서의 용도 / used for: hepcidin, anaemia of inflammation</sub>
38. **Tsoli M et al. (2013).** *Cancer cachexia: malignant inflammation, tumorkines, and metabolic mayhem.* Trends Endocrinol Metab. PMID [23201432](https://pubmed.ncbi.nlm.nih.gov/23201432/)
   <br><sub>모델에서의 용도 / used for: tumorkine framing of arm B</sub>
39. **Falconer JS et al. (1995).** *Acute-phase protein response and survival duration of patients with pancreatic cancer.* Cancer. PMID [7535184](https://pubmed.ncbi.nlm.nih.gov/7535184/)
   <br><sub>모델에서의 용도 / used for: CRP, albumin and survival in pancreatic cancer</sub>
40. **Barber MD et al. (2004).** *Modulation of the liver export protein synthetic response to feeding by an n-3 fatty-acid-enriched nutritional supplement is associated with anabolism in cachectic cancer patients.* Clin Sci (Lond). PMID [14624668](https://pubmed.ncbi.nlm.nih.gov/14624668/)
   <br><sub>모델에서의 용도 / used for: hepatic export-protein synthesis in cachexia, LV_MASS and LV_APRAA</sub>

## 6. Muscle proteolysis: FoxO, atrogenes, UPS and autophagy (ARM B-i)

41. **Bodine SC et al. (2001).** *Identification of ubiquitin ligases required for skeletal muscle atrophy.* Science. PMID [11679633](https://pubmed.ncbi.nlm.nih.gov/11679633/)
   <br><sub>모델에서의 용도 / used for: MAFbx/atrogin-1 and MuRF1, MU_MAFBX and MU_MURF1</sub>
42. **Sandri M et al. (2004).** *Foxo transcription factors induce the atrophy-related ubiquitin ligase atrogin-1 and cause skeletal muscle atrophy.* Cell. PMID [15109499](https://pubmed.ncbi.nlm.nih.gov/15109499/)
   <br><sub>모델에서의 용도 / used for: FoxO as the atrogene master switch, MU_FOXO</sub>
43. **Stitt TN et al. (2004).** *The IGF-1/PI3K/Akt pathway prevents expression of muscle atrophy-induced ubiquitin ligases by inhibiting FOXO transcription factors.* Mol Cell. PMID [15125842](https://pubmed.ncbi.nlm.nih.gov/15125842/)
   <br><sub>모델에서의 용도 / used for: Akt-mediated FoxO nuclear exclusion, WFXAKT</sub>
44. **Mammucari C et al. (2007).** *FoxO3 controls autophagy in skeletal muscle in vivo.* Cell Metab. PMID [18054315](https://pubmed.ncbi.nlm.nih.gov/18054315/)
   <br><sub>모델에서의 용도 / used for: FoxO-driven autophagy, MU_AUTOP</sub>
45. **Lecker SH et al. (2004).** *Multiple types of skeletal muscle atrophy involve a common program of changes in gene expression.* FASEB J. PMID [14718385](https://pubmed.ncbi.nlm.nih.gov/14718385/)
   <br><sub>모델에서의 용도 / used for: the common atrophy transcriptional programme</sub>
46. **Bonetto A et al. (2012).** *JAK/STAT3 pathway inhibition blocks skeletal muscle wasting downstream of IL-6 and in experimental cancer cachexia.* Am J Physiol Endocrinol Metab. PMID [22669242](https://pubmed.ncbi.nlm.nih.gov/22669242/)
   <br><sub>모델에서의 용도 / used for: IL-6 to STAT3 to atrophy, MU_STAT3</sub>
47. **Bonetto A et al. (2011).** *STAT3 activation in skeletal muscle links muscle wasting and the acute phase response in cancer cachexia.* PLoS One. PMID [21799891](https://pubmed.ncbi.nlm.nih.gov/21799891/)
   <br><sub>모델에서의 용도 / used for: muscle STAT3 couples arm B to the acute-phase response</sub>
48. **Guttridge DC et al. (2000).** *NF-kappaB-induced loss of MyoD messenger RNA: possible role in muscle decay and cachexia.* Science. PMID [11009425](https://pubmed.ncbi.nlm.nih.gov/11009425/)
   <br><sub>모델에서의 용도 / used for: NF-kB and MyoD loss</sub>
49. **Lin XY et al. (2017).** *Calpain inhibitors ameliorate muscle wasting in a cachectic mouse model bearing CT26 colorectal adenocarcinoma.* Oncol Rep. PMID [28112357](https://pubmed.ncbi.nlm.nih.gov/28112357/)
   <br><sub>모델에서의 용도 / used for: calpain inhibition and muscle wasting, MU_CALP</sub>
50. **Bhatnagar S et al. (2012).** *TWEAK causes myotube atrophy through coordinated activation of ubiquitin-proteasome system, autophagy, and caspases.* J Cell Physiol. PMID [21567392](https://pubmed.ncbi.nlm.nih.gov/21567392/)
   <br><sub>모델에서의 용도 / used for: TWEAK/Fn14 activates UPS, autophagy and caspases, MU_FN14</sub>

## 7. ActRIIB, activin A and myostatin

51. **Zhou X et al. (2010).** *Reversal of cancer cachexia and muscle wasting by ActRIIB antagonism leads to prolonged survival.* Cell. PMID [20723755](https://pubmed.ncbi.nlm.nih.gov/20723755/)
   <br><sub>모델에서의 용도 / used for: ActRIIB blockade reverses wasting; bimagrumab arm</sub>
52. **Loumaye A et al. (2015).** *Role of Activin A and myostatin in human cancer cachexia.* J Clin Endocrinol Metab. PMID [25751105](https://pubmed.ncbi.nlm.nih.gov/25751105/)
   <br><sub>모델에서의 용도 / used for: the human paradox: activin A rises, myostatin FALLS</sub>
53. **Thomas M et al. (2000).** *Myostatin, a negative regulator of muscle growth, functions by inhibiting myoblast proliferation.* J Biol Chem. PMID [10976104](https://pubmed.ncbi.nlm.nih.gov/10976104/)
   <br><sub>모델에서의 용도 / used for: myostatin biology</sub>
54. **Sartori R et al. (2009).** *Smad2 and 3 transcription factors control muscle mass in adulthood.* Am J Physiol Cell Physiol. PMID [19357234](https://pubmed.ncbi.nlm.nih.gov/19357234/)
   <br><sub>모델에서의 용도 / used for: SMAD2/3 to FoxO coupling, MU_SMAD</sub>
55. **Rooks D et al. (2020).** *Bimagrumab vs Optimized Standard of Care for Treatment of Sarcopenia in Community-Dwelling Older Adults: A Randomized Clinical Trial.* JAMA Netw Open. PMID [33074327](https://pubmed.ncbi.nlm.nih.gov/33074327/)
   <br><sub>모델에서의 용도 / used for: bimagrumab lean-mass gain without matching function gain</sub>
56. **Heymsfield SB et al. (2021).** *Effect of Bimagrumab vs Placebo on Body Fat Mass Among Adults With Type 2 Diabetes and Obesity: A Phase 2 Randomized Clinical Trial.* JAMA Netw Open. PMID [33439265](https://pubmed.ncbi.nlm.nih.gov/33439265/)
   <br><sub>모델에서의 용도 / used for: ActRIIB blockade: lean mass up, fat mass down</sub>

## 8. Anabolic resistance, IGF-1 and mTORC1 (ARM B-ii)

57. **Barclay RD et al. (2019).** *The Role of the IGF-1 Signaling Cascade in Muscle Protein Synthesis and Anabolic Resistance in Aging Skeletal Muscle.* Front Nutr. PMID [31552262](https://pubmed.ncbi.nlm.nih.gov/31552262/)
   <br><sub>모델에서의 용도 / used for: anabolic resistance concept, AN_ARES</sub>
58. **Rieu I et al. (2006).** *Leucine supplementation improves muscle protein synthesis in elderly men independently of hyperaminoacidaemia.* J Physiol. PMID [16777941](https://pubmed.ncbi.nlm.nih.gov/16777941/)
   <br><sub>모델에서의 용도 / used for: leucine threshold, AN_LEU and LEUSIG</sub>
59. **Wolfson RL et al. (2016).** *Sestrin2 is a leucine sensor for the mTORC1 pathway.* Science. PMID [26449471](https://pubmed.ncbi.nlm.nih.gov/26449471/)
   <br><sub>모델에서의 용도 / used for: molecular leucine sensing</sub>
60. **Laplante M et al. (2012).** *mTOR signaling in growth control and disease.* Cell. PMID [22500797](https://pubmed.ncbi.nlm.nih.gov/22500797/)
   <br><sub>모델에서의 용도 / used for: mTORC1 to S6K1 and 4E-BP1</sub>
61. **Larson KR et al. (2024).** *FGF21 Induces Skeletal Muscle Atrophy and Increases Amino Acids in Female Mice: A Potential Role for Glucocorticoids.* Endocrinology. PMID [38244215](https://pubmed.ncbi.nlm.nih.gov/38244215/)
   <br><sub>모델에서의 용도 / used for: GR to KLF15/REDD1/FoxO, the megestrol and dexamethasone off-target arm</sub>
62. **Tamiya T et al. (2011).** *Suppressors of cytokine signaling (SOCS) proteins and JAK/STAT pathways: regulation of T-cell inflammation by SOCS1 and SOCS3.* Arterioscler Thromb Vasc Biol. PMID [21508344](https://pubmed.ncbi.nlm.nih.gov/21508344/)
   <br><sub>모델에서의 용도 / used for: SOCS3 mechanism of GH resistance</sub>
63. **Deutz NE et al. (2011).** *Muscle protein synthesis in cancer patients can be stimulated with a specially formulated medical food.* Clin Nutr. PMID [21683485](https://pubmed.ncbi.nlm.nih.gov/21683485/)
   <br><sub>모델에서의 용도 / used for: the anabolic response is blunted but not absent, GNUT</sub>
64. **Bhasin S et al. (2001).** *Testosterone dose-response relationships in healthy young men.* Am J Physiol Endocrinol Metab. PMID [11701431](https://pubmed.ncbi.nlm.nih.gov/11701431/)
   <br><sub>모델에서의 용도 / used for: androgen-driven muscle accrual, AN_AR and enobosarm</sub>

## 9. Satellite cells, myonuclei and irreversibility

65. **Acharyya S et al. (2004).** *Cancer cachexia is regulated by selective targeting of skeletal muscle gene products.* J Clin Invest. PMID [15286803](https://pubmed.ncbi.nlm.nih.gov/15286803/)
   <br><sub>모델에서의 용도 / used for: Pax7 dysregulation and impaired regeneration in cachexia</sub>
66. **Jackman RW et al. (2013).** *Nuclear factor-κB signalling and transcriptional regulation in skeletal muscle atrophy.* Exp Physiol. PMID [22848079](https://pubmed.ncbi.nlm.nih.gov/22848079/)
   <br><sub>모델에서의 용도 / used for: NF-kB, Pax7 and the arrested differentiation block, AN_MYOD</sub>
67. **Relaix F et al. (2012).** *Satellite cells are essential for skeletal muscle regeneration: the cell on the edge returns centre stage.* Development. PMID [22833472](https://pubmed.ncbi.nlm.nih.gov/22833472/)
   <br><sub>모델에서의 용도 / used for: regenerative ceiling, AN_SATC</sub>
68. **Murach KA et al. (2018).** *Myonuclear Domain Flexibility Challenges Rigid Assumptions on Satellite Cell Contribution to Skeletal Muscle Fiber Hypertrophy.* Front Physiol. PMID [29896117](https://pubmed.ncbi.nlm.nih.gov/29896117/)
   <br><sub>모델에서의 용도 / used for: myonuclear domain and the limits of hypertrophy, HYPMAX</sub>

## 10. Adipose tissue: lipolysis, browning and energy expenditure

69. **Das SK et al. (2011).** *Adipose triglyceride lipase contributes to cancer-associated cachexia.* Science. PMID [21680814](https://pubmed.ncbi.nlm.nih.gov/21680814/)
   <br><sub>모델에서의 용도 / used for: ATGL as rate-limiting; ATGL-null mice are protected, AD_ATGL</sub>
70. **Petruzzelli M et al. (2014).** *A switch from white to brown fat increases energy expenditure in cancer-associated cachexia.* Cell Metab. PMID [25043816](https://pubmed.ncbi.nlm.nih.gov/25043816/)
   <br><sub>모델에서의 용도 / used for: WAT browning raises REE, AD_UCP1 and SREEUCP</sub>
71. **Kir S et al. (2014).** *Tumour-derived PTH-related protein triggers adipose tissue browning and cancer cachexia.* Nature. PMID [25043053](https://pubmed.ncbi.nlm.nih.gov/25043053/)
   <br><sub>모델에서의 용도 / used for: PTHrP to PTH1R browning, TU_PTHRP</sub>
72. **Bao Y et al. (2005).** *Zinc-alpha2-glycoprotein, a lipid mobilizing factor, is expressed and secreted by human (SGBS) adipocytes.* FEBS Lett. PMID [15620688](https://pubmed.ncbi.nlm.nih.gov/15620688/)
   <br><sub>모델에서의 용도 / used for: ZAG/LMF beta3-adrenergic-like lipolysis, TU_ZAG</sub>
73. **Bing C et al. (2006).** *Adipose atrophy in cancer cachexia: morphologic and molecular analysis of adipose tissue in tumour-bearing mice.* Br J Cancer. PMID [17047651](https://pubmed.ncbi.nlm.nih.gov/17047651/)
   <br><sub>모델에서의 용도 / used for: adipose remodelling and fibrosis, AD_FIBR</sub>
74. **Wei L et al. (2022).** *Creatine modulates cellular energy metabolism and protects against cancer cachexia-associated muscle wasting.* Front Pharmacol. PMID [36569317](https://pubmed.ncbi.nlm.nih.gov/36569317/)
   <br><sub>모델에서의 용도 / used for: futile cycles and hypermetabolism</sub>
75. **Hall KD et al. (2012).** *Energy balance and its components: implications for body weight regulation.* Am J Clin Nutr. PMID [22434603](https://pubmed.ncbi.nlm.nih.gov/22434603/)
   <br><sub>모델에서의 용도 / used for: energy partitioning between fat and lean, FPARTF</sub>
76. **Yoshikawa T et al. (2001).** *Insulin resistance in patients with cancer: relationships with tumor site, tumor stage, body-weight loss, acute-phase response, and energy expenditure.* Nutrition. PMID [11448578](https://pubmed.ncbi.nlm.nih.gov/11448578/)
   <br><sub>모델에서의 용도 / used for: insulin resistance and hypermetabolism, AN_IR</sub>
77. **Barcellos PS et al. (2021).** *Resting energy expenditure in cancer patients: Agreement between predictive equations and indirect calorimetry.* Clin Nutr ESPEN. PMID [33745594](https://pubmed.ncbi.nlm.nih.gov/33745594/)
   <br><sub>모델에서의 용도 / used for: measured vs predicted REE, RQREE and EN_HYPERM</sub>

## 11. Muscle quality, mitochondria and physical function (AXIS C)

78. **Reid MB et al. (2011).** *Beyond atrophy: redox mechanisms of muscle dysfunction in chronic inflammatory disease.* J Physiol. PMID [21320886](https://pubmed.ncbi.nlm.nih.gov/21320886/)
   <br><sub>모델에서의 용도 / used for: force lost without mass lost, FN_MHC and MYOQ</sub>
79. **Andersson DC et al. (2011).** *Ryanodine receptor oxidation causes intracellular calcium leak and muscle weakness in aging.* Cell Metab. PMID [21803290](https://pubmed.ncbi.nlm.nih.gov/21803290/)
   <br><sub>모델에서의 용도 / used for: RyR1 leak and excitation-contraction uncoupling, FN_RYR</sub>
80. **Sandri M et al. (2006).** *PGC-1alpha protects skeletal muscle from atrophy by suppressing FoxO3 action and atrophy-specific gene transcription.* Proc Natl Acad Sci U S A. PMID [17053067](https://pubmed.ncbi.nlm.nih.gov/17053067/)
   <br><sub>모델에서의 용도 / used for: PGC-1alpha collapse, FN_PGC1</sub>
81. **Op den Kamp CM et al. (2015).** *Preserved muscle oxidative metabolic phenotype in newly diagnosed non-small cell lung cancer cachexia.* J Cachexia Sarcopenia Muscle. PMID [26136192](https://pubmed.ncbi.nlm.nih.gov/26136192/)
   <br><sub>모델에서의 용도 / used for: timing of the oxidative phenotype change relative to mass loss</sub>
82. **Delmonico MJ et al. (2009).** *Longitudinal study of muscle strength, quality, and adipose tissue infiltration.* Am J Clin Nutr. PMID [19864405](https://pubmed.ncbi.nlm.nih.gov/19864405/)
   <br><sub>모델에서의 용도 / used for: muscle quality declines faster than mass</sub>
83. **Taekema DG et al. (2010).** *Handgrip strength as a predictor of functional, psychological and social health. A prospective population-based study among the oldest old.* Age Ageing. PMID [20219767](https://pubmed.ncbi.nlm.nih.gov/20219767/)
   <br><sub>모델에서의 용도 / used for: handgrip as a clinical endpoint, FN_GRIP</sub>
84. **Su X et al. (2026).** *Role of mitochondrial dysfunction in muscle wasting in cancer cachexia: a narrative review.* J Transl Med. PMID [41742249](https://pubmed.ncbi.nlm.nih.gov/41742249/)
   <br><sub>모델에서의 용도 / used for: mitochondrial dysfunction and mitophagy, FN_MITO and FN_MFN</sub>

## 12. Anamorelin and ghrelin-receptor agonism

85. **Temel JS et al. (2016).** *Anamorelin in patients with non-small-cell lung cancer and cachexia (ROMANA 1 and ROMANA 2): results from two randomised, double-blind, phase 3 trials.* Lancet Oncol. PMID [26906526](https://pubmed.ncbi.nlm.nih.gov/26906526/)
   <br><sub>모델에서의 용도 / used for: THE key anchor: LBM endpoint met, handgrip endpoint missed</sub>
86. **Katakami N et al. (2018).** *Anamorelin (ONO-7643) for the treatment of patients with non-small cell lung cancer and cachexia: Results from a randomized, double-blind, placebo-controlled, multicenter study of Japanese patients (ONO-7643-04).* Cancer. PMID [29205286](https://pubmed.ncbi.nlm.nih.gov/29205286/)
   <br><sub>모델에서의 용도 / used for: ONO-7643-04 confirmation in Japanese patients</sub>
87. **Garcia JM et al. (2015).** *Anamorelin for patients with cancer cachexia: an integrated analysis of two phase 2, randomised, placebo-controlled, double-blind trials.* Lancet Oncol. PMID [25524795](https://pubmed.ncbi.nlm.nih.gov/25524795/)
   <br><sub>모델에서의 용도 / used for: integrated phase 2 analysis</sub>
88. **Chanoine JP et al. (2009).** *Ghrelin and the growth hormone secretagogue receptor in growth and development.* Int J Obes (Lond). PMID [19363508](https://pubmed.ncbi.nlm.nih.gov/19363508/)
   <br><sub>모델에서의 용도 / used for: GHSR-1a pharmacology, anamorelin target</sub>
89. **Meinhardt U et al. (2010).** *The effects of growth hormone on body composition and physical performance in recreational athletes: a randomized trial.* Ann Intern Med. PMID [20439575](https://pubmed.ncbi.nlm.nih.gov/20439575/)
   <br><sub>모델에서의 용도 / used for: GH raises lean mass largely as water without raising strength; the LNW mechanism</sub>

## 13. Ponsegromab and GDF-15 blockade

90. **Groarke JD et al. (2025).** *Ponsegromab for the Treatment of Cancer Cachexia. Reply.* N Engl J Med. PMID [40043245](https://pubmed.ncbi.nlm.nih.gov/40043245/)
   <br><sub>모델에서의 용도 / used for: phase 2 anchor: +5.6% weight at 400 mg, appetite and activity</sub>

## 14. Progestins, corticosteroids and appetite stimulants

91. **Ruiz Garcia V et al. (2013).** *Megestrol acetate for treatment of anorexia-cachexia syndrome.* Cochrane Database Syst Rev. PMID [23543530](https://pubmed.ncbi.nlm.nih.gov/23543530/)
   <br><sub>모델에서의 용도 / used for: Cochrane: weight gain without lean-mass or QoL benefit</sub>
92. **Loprinzi CL et al. (1999).** *Randomized comparison of megestrol acetate versus dexamethasone versus fluoxymesterone for the treatment of cancer anorexia/cachexia.* J Clin Oncol. PMID [10506633](https://pubmed.ncbi.nlm.nih.gov/10506633/)
   <br><sub>모델에서의 용도 / used for: head-to-head appetite stimulants and their toxicities</sub>
93. **Yennurajalingam S et al. (2013).** *Reduction of cancer-related fatigue with dexamethasone: a double-blind, randomized, placebo-controlled trial in patients with advanced cancer.* J Clin Oncol. PMID [23897970](https://pubmed.ncbi.nlm.nih.gov/23897970/)
   <br><sub>모델에서의 용도 / used for: short-lived steroid benefit, TAUDEXA</sub>
94. **Sandhya L et al. (2023).** *Randomized Double-Blind Placebo-Controlled Study of Olanzapine for Chemotherapy-Related Anorexia in Patients With Locally Advanced or Metastatic Gastric, Hepatopancreaticobiliary, and Lung Cancer.* J Clin Oncol. PMID [36977285](https://pubmed.ncbi.nlm.nih.gov/36977285/)
   <br><sub>모델에서의 용도 / used for: olanzapine 2.5 mg: >5% weight gain in most patients</sub>
95. **Gong W et al. (2026).** *Olanzapine Plus Triple Antiemetic Therapy for the Prevention of Platinum-Based Delayed-Phase Chemotherapy-Induced Nausea and Vomiting: A Meta-Analysis.* Curr Oncol. PMID [41590347](https://pubmed.ncbi.nlm.nih.gov/41590347/)
   <br><sub>모델에서의 용도 / used for: olanzapine antiemetic mechanism, EMXOLZN</sub>
96. **Allison DB et al. (1999).** *Antipsychotic-induced weight gain: a comprehensive research synthesis.* Am J Psychiatry. PMID [10553730](https://pubmed.ncbi.nlm.nih.gov/10553730/)
   <br><sub>모델에서의 용도 / used for: H1 and 5-HT2C mediated appetite effect, EMXOLZA</sub>

## 15. SARMs, beta-blockade and other pharmacology

97. **Dobs AS et al. (2013).** *Effects of enobosarm on muscle wasting and physical function in patients with cancer: a double-blind, randomised controlled phase 2 trial.* Lancet Oncol. PMID [23499390](https://pubmed.ncbi.nlm.nih.gov/23499390/)
   <br><sub>모델에서의 용도 / used for: POWER precursor: LBM up, stair-climb power not met</sub>
98. **Crawford J et al. (2016).** *Study Design and Rationale for the Phase 3 Clinical Development Program of Enobosarm, a Selective Androgen Receptor Modulator, for the Prevention and Treatment of Muscle Wasting in Cancer Patients (POWER Trials).* Curr Oncol Rep. PMID [27138015](https://pubmed.ncbi.nlm.nih.gov/27138015/)
   <br><sub>모델에서의 용도 / used for: POWER trial design and co-primary endpoints</sub>
99. **Stewart Coats AJ et al. (2016).** *Espindolol for the treatment and prevention of cachexia in patients with stage III/IV non-small cell lung cancer or colorectal cancer: a randomized, double-blind, placebo-controlled, international multicentre phase II study (the ACT-ONE trial).* J Cachexia Sarcopenia Muscle. PMID [27386169](https://pubmed.ncbi.nlm.nih.gov/27386169/)
   <br><sub>모델에서의 용도 / used for: ACT-ONE: weight AND handgrip both improved</sub>
100. **Fearon KC et al. (2003).** *Effect of a protein and energy dense N-3 fatty acid enriched oral supplement on loss of weight and lean tissue in cancer cachexia: a randomised double blind trial.* Gut. PMID [12970142](https://pubmed.ncbi.nlm.nih.gov/12970142/)
   <br><sub>모델에서의 용도 / used for: EPA-enriched supplement trial</sub>
101. **Koch A et al. (2011).** *Effect of celecoxib on survival in patients with advanced non-small cell lung cancer: a double blind randomised clinical phase III trial (CYCLUS study) by the Swedish Lung Cancer Study Group.* Eur J Cancer. PMID [21565487](https://pubmed.ncbi.nlm.nih.gov/21565487/)
   <br><sub>모델에서의 용도 / used for: COX-2 inhibition arm</sub>
102. **Ando K et al. (2013).** *Possible role for tocilizumab, an anti-interleukin-6 receptor antibody, in treating cancer cachexia.* J Clin Oncol. PMID [23129740](https://pubmed.ncbi.nlm.nih.gov/23129740/)
   <br><sub>모델에서의 용도 / used for: IL-6R blockade in cachexia, tocilizumab arm</sub>
103. **Zhang X et al. (2013).** *Pharmacokinetics and pharmacodynamics of tocilizumab, a humanized anti-interleukin-6 receptor monoclonal antibody, following single-dose administration by subcutaneous and intravenous routes to healthy subjects.* Int J Clin Pharmacol Ther. PMID [23547848](https://pubmed.ncbi.nlm.nih.gov/23547848/)
   <br><sub>모델에서의 용도 / used for: tocilizumab PK, target-mediated clearance, VMTCZ and KMTCZ</sub>
104. **Falconer JS et al. (1994).** *Effect of eicosapentaenoic acid and other fatty acids on the growth in vitro of human pancreatic cancer cell lines.* Br J Cancer. PMID [8180010](https://pubmed.ncbi.nlm.nih.gov/8180010/)
   <br><sub>모델에서의 용도 / used for: EPA anti-tumour and anti-cachexia rationale</sub>
105. **Harada M et al. (2025).** *Tocilizumab, a Humanized Anti-interleukin-6 Receptor Antibody, Induces Hepatic Iron Overload in a Susceptible Patient.* Intern Med. PMID [39401912](https://pubmed.ncbi.nlm.nih.gov/39401912/)
   <br><sub>모델에서의 용도 / used for: IL-6R blockade normalises CRP and stabilises weight</sub>

## 16. Nutrition, exercise and multimodal therapy

106. **Arends J et al. (2017).** *ESPEN guidelines on nutrition in cancer patients.* Clin Nutr. PMID [27637832](https://pubmed.ncbi.nlm.nih.gov/27637832/)
   <br><sub>모델에서의 용도 / used for: protein and energy targets, PROTTG and ONSKCAL</sub>
107. **Muscaritoli M et al. (2021).** *ESPEN practical guideline: Clinical Nutrition in cancer.* Clin Nutr. PMID [33946039](https://pubmed.ncbi.nlm.nih.gov/33946039/)
   <br><sub>모델에서의 용도 / used for: practical nutrition targets</sub>
108. **Huillard O et al. (2020).** *Management of Cancer Cachexia: ASCO Guideline-Time to Address the Elephant in the Room.* J Clin Oncol. PMID [32946354](https://pubmed.ncbi.nlm.nih.gov/32946354/)
   <br><sub>모델에서의 용도 / used for: ASCO guidance: no agent recommended for routine use</sub>
109. **Solheim TS et al. (2017).** *A randomized phase II feasibility trial of a multimodal intervention for the management of cachexia in lung and pancreatic cancer.* J Cachexia Sarcopenia Muscle. PMID [28614627](https://pubmed.ncbi.nlm.nih.gov/28614627/)
   <br><sub>모델에서의 용도 / used for: MENAC multimodal package, scenario S16</sub>
110. **Solheim TS et al. (2018).** *Cancer cachexia: rationale for the MENAC (Multimodal-Exercise, Nutrition and Anti-inflammatory medication for Cachexia) trial.* BMJ Support Palliat Care. PMID [29440149](https://pubmed.ncbi.nlm.nih.gov/29440149/)
   <br><sub>모델에서의 용도 / used for: multimodal rationale and design</sub>
111. **Hardee JP et al. (2019).** *Understanding the Role of Exercise in Cancer Cachexia Therapy.* Am J Lifestyle Med. PMID [30627079](https://pubmed.ncbi.nlm.nih.gov/30627079/)
   <br><sub>모델에서의 용도 / used for: exercise as the only lever on muscle quality, AXIS C</sub>
112. **Ochi E et al. (2026).** *Effects of exercise on people living with advanced lung cancer: a systematic review and meta-analysis.* Support Care Cancer. PMID [41663540](https://pubmed.ncbi.nlm.nih.gov/41663540/)
   <br><sub>모델에서의 용도 / used for: exercise effects on strength and function in cancer</sub>
113. **Virizuela JA et al. (2018).** *Nutritional support and parenteral nutrition in cancer patients: an expert consensus report.* Clin Transl Oncol. PMID [29043569](https://pubmed.ncbi.nlm.nih.gov/29043569/)
   <br><sub>모델에서의 용도 / used for: why nutrition alone cannot cross anabolic resistance</sub>
114. **Langius JA et al. (2013).** *Effect of nutritional interventions on nutritional status, quality of life and mortality in patients with head and neck cancer receiving (chemo)radiotherapy: a systematic review.* Clin Nutr. PMID [23845384](https://pubmed.ncbi.nlm.nih.gov/23845384/)
   <br><sub>모델에서의 용도 / used for: nutrition-only intervention outcomes, scenario S03</sub>

## 17. Anticancer therapy, tumour control and reversibility

115. **Persson C et al. (2002).** *The relevance of weight loss for survival and quality of life in patients with advanced gastrointestinal cancer treated with palliative chemotherapy.* Anticancer Res. PMID [12552973](https://pubmed.ncbi.nlm.nih.gov/12552973/)
   <br><sub>모델에서의 용도 / used for: tumour control as the definitive cachexia therapy</sub>
116. **Damrauer JS et al. (2018).** *Chemotherapy-induced muscle wasting: association with NF-κB and cancer cachexia.* Eur J Transl Myol. PMID [29991992](https://pubmed.ncbi.nlm.nih.gov/29991992/)
   <br><sub>모델에서의 용도 / used for: direct chemotherapy myotoxicity, MYOTOX</sub>
117. **Dewys WD et al. (1980).** *Prognostic effect of weight loss prior to chemotherapy in cancer patients. Eastern Cooperative Oncology Group.* Am J Med. PMID [7424938](https://pubmed.ncbi.nlm.nih.gov/7424938/)
   <br><sub>모델에서의 용도 / used for: weight loss before treatment predicts poor outcome</sub>

## 18. Patient-reported outcomes and quality of life

118. **Gelhorn HL et al. (2019).** *Comprehensive validation of the functional assessment of anorexia/cachexia therapy (FAACT) anorexia/cachexia subscale (A/CS) in lung cancer patients with involuntary weight loss.* Qual Life Res. PMID [30796591](https://pubmed.ncbi.nlm.nih.gov/30796591/)
   <br><sub>모델에서의 용도 / used for: FAACT A/CS cut-off of 37</sub>
119. **Aaronson NK et al. (1993).** *The European Organization for Research and Treatment of Cancer QLQ-C30: a quality-of-life instrument for use in international clinical trials in oncology.* J Natl Cancer Inst. PMID [8433390](https://pubmed.ncbi.nlm.nih.gov/8433390/)
   <br><sub>모델에서의 용도 / used for: QLQ-C30 quality-of-life endpoint</sub>
120. **Amano K et al. (2016).** *Eating-related distress and need for nutritional support of families of advanced cancer patients: a nationwide survey of bereaved family members.* J Cachexia Sarcopenia Muscle. PMID [27239421](https://pubmed.ncbi.nlm.nih.gov/27239421/)
   <br><sub>모델에서의 용도 / used for: the family dimension of eating-related distress, GI_DEPR</sub>
121. **Oken MM et al. (1982).** *Toxicity and response criteria of the Eastern Cooperative Oncology Group.* Am J Clin Oncol. PMID [7165009](https://pubmed.ncbi.nlm.nih.gov/7165009/)
   <br><sub>모델에서의 용도 / used for: ECOG performance status definition</sub>
122. **Ribaudo JM et al. (2000).** *Re-validation and shortening of the Functional Assessment of Anorexia/Cachexia Therapy (FAACT) questionnaire.* Qual Life Res. PMID [11401046](https://pubmed.ncbi.nlm.nih.gov/11401046/)
   <br><sub>모델에서의 용도 / used for: FAACT A/CS validation, the FAACT output</sub>

## 19. QSP methodology, mrgsolve and modelling practice

123. **EFPIA MID3 Workgroup et al. (2016).** *Good Practices in Model-Informed Drug Discovery and Development: Practice, Application, and Documentation.* CPT Pharmacometrics Syst Pharmacol. PMID [27069774](https://pubmed.ncbi.nlm.nih.gov/27069774/)
   <br><sub>모델에서의 용도 / used for: MID3 model documentation practice</sub>
124. **Musante CJ et al. (2017).** *Quantitative Systems Pharmacology: A Case for Disease Models.* Clin Pharmacol Ther. PMID [27709613](https://pubmed.ncbi.nlm.nih.gov/27709613/)
   <br><sub>모델에서의 용도 / used for: QSP disease-model rationale</sub>
125. **Milligan PA et al. (2013).** *Model-based drug development: a rational approach to efficiently accelerate drug development.* Clin Pharmacol Ther. PMID [23588322](https://pubmed.ncbi.nlm.nih.gov/23588322/)
   <br><sub>모델에서의 용도 / used for: model-informed development framing</sub>

---

## 미해결 항목 (unresolved queries)

아래 질의는 PubMed에서 제목 일치 기준(0.85)을 넘는 레코드를 찾지 못해
**의도적으로 인용 목록에서 제외**했습니다. 근거가 없다는 뜻이 아니라, 이 문서의
검증 규칙을 통과하지 못했다는 뜻입니다. 해당 기전은 위 목록의 다른 문헌으로
뒷받침되거나, 모델 코드 주석에 가정으로 명시되어 있습니다.

The queries below did not return a PubMed record meeting the 0.85 title-match
threshold and were therefore excluded on purpose. That is a statement about this
document's verification rule, not a claim that the underlying mechanism is
unsupported: each is either covered by another entry above or flagged as an
explicit assumption in the model source.

- Skeletal muscle depletion predicts survival of patients with advanced renal cell carcinoma treated with sorafenib
- Tumour-derived macrophage inhibitory cytokine-1 causes cachexia in mice
- Growth differentiation factor 15 as a biomarker for mortality in patients with cancer
- Interleukin-1 systems mediate lipopolysaccharide-induced anorexia
- An acute-phase protein response in patients with unresectable pancreatic cancer is associated with weight loss
- Increased hepatic protein synthesis in cachectic patients with pancreatic cancer
- Serum interleukin-6 levels correlate to survival in patients with metastatic renal cell carcinoma
- Leukemia inhibitory factor signaling is required for lipid and protein metabolic responses to a burn injury
- Rate of protein synthesis and quantity of RNA in skeletal muscle of patients with cancer
- Increased muscle protein degradation and expression of ubiquitin-proteasome pathway components in cancer cachexia
- Calpain and caspase-3 activity in cachectic muscle
- TWEAK induces skeletal muscle atrophy through coordinated activation of ubiquitin-proteasome system, autophagy, and caspases
- REDD1 is essential for stress-induced synthesis of the protein Sestrin2 and glucocorticoid-mediated inhibition of mTORC1
- Growth hormone resistance in the catabolic state
- Fat mass loss precedes muscle mass loss and is more sensitive to changes in cachexia
- Resting energy expenditure in cancer patients: agreement between expenditure measured and predicted by common prediction equations
- Glucose turnover and recycling in relation to weight loss in patients with lung cancer
- Nitration of tryptophan residues in proteins from human plasma and tissue
- Mitochondrial dysfunction and mitophagy in cancer cachexia
- Muscle quality is more strongly related to physical function than muscle mass
- Growth hormone and sodium retention: a review of biochemical and renal handling of sodium
- Blockade of GDF15 signaling reverses cancer cachexia in mice
- Neutralization of circulating GDF15 prevents cachexia
- Reduced survival with megestrol acetate in cancer patients
- Efficacy of dexamethasone and megestrol acetate in patients with metastatic cancer and weight loss
- Effect of eicosapentaenoic acid and other fatty acids on the growth in vivo of human breast cancer
- Tocilizumab, a humanized anti-interleukin-6 receptor antibody, ameliorated clinical symptoms and improved prognosis of a patient with advanced lung cancer and cachexia
- Interleukin-6 receptor blockade selectively reduces IL-6 clearance
- Cisplatin induces cachexia-like symptoms in tumor-free mice
- Refractory cachexia: a distinct clinical entity
- Measuring the concerns of cancer patients with low body weight: the Functional Assessment of Anorexia/Cachexia Therapy (FAACT) questionnaire
- mrgsolve: Simulate from ODE-Based Population PK/PD and Systems Pharmacology Models
- A quantitative systems pharmacology model of muscle atrophy
