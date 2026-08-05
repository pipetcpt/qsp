# Moyamoya Disease (MMD) — QSP Model References

**Every PMID in this file was resolved and its title, journal, year and first
authors were read back from the NCBI E-utilities `esummary` endpoint at the time
of writing.** Nothing here is quoted from memory. One candidate reference was
dropped during that check because the PMID resolved to an unrelated paper — which
is the failure mode this procedure exists to catch.

99 references, grouped by the part of the model each one supports.

---

## How to read this list

The model's thesis is that **reserve, not flow, is the state variable**: cerebral
blood flow is defended by a cortical arteriole that is a variable resistor with a
floor, and almost everything clinically distinctive about moyamoya follows from
what happens when that floor is reached. The sections below are ordered so that
each one supplies a specific term in that argument, and each section header note
says which term and how the reference was used — in particular whether a number
was **fitted** to it or **predicted** against it.

Two calibration choices are worth flagging up front:

- **CBF_crit is derived, not fitted.** CMRO2/(CaO2 x OEF_max) = 19.71 mL/100g/min
  at Hb 15 g/dL. The textbook penumbral threshold of ~20 falls out of the oxygen
  arithmetic; the references in section 5 are the warrant for CMRO2, OEF_max and
  CaO2 individually, not for the threshold itself.
- **The JAM Trial is used asymmetrically.** Only the conservative arm's 5-year
  rebleeding rate set `HEM_HAZ0`. The bypass hazard ratio is therefore a model
  prediction that can be checked against the published trial result.

---

## 1. Definition, epidemiology and natural history

Suzuki & Takaku's original description supplies the angiographic staging the
model's `STEN` state is read against; Kobayashi 2000 and Morioka 2003 supply the
haemorrhagic natural history that `HEM_HAZ0` is calibrated to; Kuroda 2007 and
the AMORE work supply the asymptomatic arm that the model has to predict is a
*different decision*.

1. Suzuki J, Takaku A Cerebrovascular "moyamoya" disease. Disease showing abnormal net-like vessels in base of brain. *Arch Neurol* 1969. PMID [5775283](https://pubmed.ncbi.nlm.nih.gov/5775283/)
2. Kuroda S, Houkin K Moyamoya disease: current concepts and future perspectives. *Lancet Neurol* 2008. PMID [18940695](https://pubmed.ncbi.nlm.nih.gov/18940695/)
3. Scott RM, Smith ER Moyamoya disease and moyamoya syndrome. *N Engl J Med* 2009. PMID [19297575](https://pubmed.ncbi.nlm.nih.gov/19297575/)
4. Fujimura M, Tominaga T, Kuroda S et al. 2021 Japanese Guidelines for the Management of Moyamoya Disease: Guidelines from the Research Committee on Moyamoya Disease and Japan Stroke Society. *Neurol Med Chir (Tokyo)* 2022. PMID [35197402](https://pubmed.ncbi.nlm.nih.gov/35197402/)
5. Rifino N, Hervè D, Acerbi F et al. Diagnosis and management of adult Moyamoya angiopathy: An overview of guideline recommendations and identification of future research directions. *Int J Stroke* 2025. PMID [39425621](https://pubmed.ncbi.nlm.nih.gov/39425621/)
6. Baba T, Houkin K, Kuroda S Novel epidemiological features of moyamoya disease. *J Neurol Neurosurg Psychiatry* 2008. PMID [18077479](https://pubmed.ncbi.nlm.nih.gov/18077479/)
7. Hayashi K, Horie N, Suyama K et al. An epidemiological survey of moyamoya disease, unilateral moyamoya disease and quasi-moyamoya disease in Japan. *Clin Neurol Neurosurg* 2013. PMID [23041378](https://pubmed.ncbi.nlm.nih.gov/23041378/)
8. Goyal MS, Hallemeier CL, Zipfel GJ et al. Clinical features and outcome in North American adults with idiopathic basal arterial occlusive disease without moyamoya collaterals. *Neurosurgery* 2010. PMID [20562658](https://pubmed.ncbi.nlm.nih.gov/20562658/)
9. Khan N, Achrol AS, Guzman R et al. Sex differences in clinical presentation and treatment outcomes in Moyamoya disease. *Neurosurgery* 2012. PMID [22718024](https://pubmed.ncbi.nlm.nih.gov/22718024/)
10. Kuroda S, Hashimoto N, Yoshimoto T et al. Radiological findings, clinical course, and outcome in asymptomatic moyamoya disease: results of multicenter survey in Japan. *Stroke* 2007. PMID [17395863](https://pubmed.ncbi.nlm.nih.gov/17395863/)
11. Kuroda S [Pathophysiology and Outcome of Asymptomatic Moyamoya Disease]. *Brain Nerve* 2026. PMID [41804510](https://pubmed.ncbi.nlm.nih.gov/41804510/)
12. Kobayashi E, Saeki N, Oishi H et al. Long-term natural history of hemorrhagic moyamoya disease in 42 patients. *J Neurosurg* 2000. PMID [11117870](https://pubmed.ncbi.nlm.nih.gov/11117870/)
13. Morioka M, Hamada J, Todaka T et al. High-risk age for rebleeding in patients with hemorrhagic moyamoya disease: long-term follow-up study. *Neurosurgery* 2003. PMID [12699546](https://pubmed.ncbi.nlm.nih.gov/12699546/)
14. Hishikawa T, Murai S, Ito YM et al. The AMORE Score for Predicting the 5-Year Risk of Hemorrhagic Stroke in Asymptomatic Moyamoya Disease. *Cerebrovasc Dis* 2026. PMID [41615863](https://pubmed.ncbi.nlm.nih.gov/41615863/)
15. Yamamoto S, Funaki T, Fujimura M et al. Development of Hemorrhage-prone Anastomoses in Asymptomatic Moyamoya Disease-A Comparative Study with Japan Adult Moyamoya Trial. *J Stroke Cerebrovasc Dis* 2019. PMID [31471213](https://pubmed.ncbi.nlm.nih.gov/31471213/)
16. Ogasawara K, Misaki T, Miyoshi K et al. Long-term outcomes of nonsurgical management alone for adult patients with ischemic moyamoya disease without cerebral misery perfusion: combined analyses of two prospective cohorts. *Cerebrovasc Dis* 2026. PMID [42497101](https://pubmed.ncbi.nlm.nih.gov/42497101/)

## 2. RNF213 and the genetics of the intimal lesion

The model gives RNF213 **two entries of opposite sign**: it raises `K_SMC`
(the intimal lesion) and it lowers the collateral ceiling through `RNF_CAP`.
Kamada 2011 is the GWAS; Miyatake 2012 supplies the homozygote/heterozygote
gene-dose contrast the model reproduces; Sonobe 2014 and Li 2024 supply the
MMP-9 link; Choi 2025 supplies the caveolin-1 arm; Kawabori 2023 reports that
the variant impairs *post-operative* collateral development, which is the
clinical form of the ceiling.

17. Kamada F, Aoki Y, Narisawa A et al. A genome-wide association study identifies RNF213 as the first Moyamoya disease gene. *J Hum Genet* 2011. PMID [21048783](https://pubmed.ncbi.nlm.nih.gov/21048783/)
18. Miyatake S, Touho H, Miyake N et al. Sibling cases of moyamoya disease having homozygous and heterozygous c.14576G>A variant in RNF213 showed varying clinical course and severity. *J Hum Genet* 2012. PMID [22931863](https://pubmed.ncbi.nlm.nih.gov/22931863/)
19. Koizumi A, Kobayashi H, Hitomi T et al. A new horizon of moyamoya disease and associated health risks explored through RNF213. *Environ Health Prev Med* 2016. PMID [26662949](https://pubmed.ncbi.nlm.nih.gov/26662949/)
20. Morito D Molecular structure and function of mysterin/RNF213. *J Biochem* 2024. PMID [38378744](https://pubmed.ncbi.nlm.nih.gov/38378744/)
21. Choi J, Inoue R, Masuo Y et al. RNF213 Acts as a Molecular Switch for Cav-1 Ubiquitination and Phosphorylation in Human Cells. *Cells* 2025. PMID [40497951](https://pubmed.ncbi.nlm.nih.gov/40497951/)
22. Mineharu Y, Kamata T, Tomoto M et al. Peripheral blood GATA2 expression impacts RNF213 mutation penetrance and clinical severity in moyamoya disease. *Stroke Vasc Neurol* 2025. PMID [40268337](https://pubmed.ncbi.nlm.nih.gov/40268337/)
23. Akagawa H [Association of Rare RNF213 Variants and Moyamoya Disease]. *No Shinkei Geka* 2025. PMID [40438012](https://pubmed.ncbi.nlm.nih.gov/40438012/)
24. Kawabori M, Ito M, Kazumata K et al. Impact of RNF213 c.14576G>A Variant on the Development of Direct and Indirect Revascularization in Pediatric Moyamoya Disease. *Cerebrovasc Dis* 2023. PMID [36063804](https://pubmed.ncbi.nlm.nih.gov/36063804/)
25. Wu M, Li S, Liu W et al. RNF213 deficiency in human iPSC-derived vascular organoids captures key feature of moyamoya disease vasculopathy. *Microvasc Res* 2026. PMID [42289232](https://pubmed.ncbi.nlm.nih.gov/42289232/)
26. Sonobe S, Fujimura M, Niizuma K et al. Increased vascular MMP-9 in mice lacking RNF213: moyamoya disease susceptibility gene. *Neuroreport* 2014. PMID [25383461](https://pubmed.ncbi.nlm.nih.gov/25383461/)
27. Li Z, Liu Y, Li X et al. Knockdown the moyamoya disease susceptibility gene, RNF213, upregulates the expression of basic fibroblast growth factor and matrix metalloproteinase-9 in bone marrow derived mesenchymal stem cells. *Neurosurg Rev* 2024. PMID [38811382](https://pubmed.ncbi.nlm.nih.gov/38811382/)
28. Wang C, Sun C, Zhao Y et al. RNF213 gene silencing upregulates transforming growth factor β1 expression in bone marrow-derived mesenchymal stem cells and is involved in the onset of Moyamoya disease. *Exp Ther Med* 2021. PMID [34373710](https://pubmed.ncbi.nlm.nih.gov/34373710/)

## 3. Histopathology of the terminal ICA lesion

These establish that the lesion is luminal smooth-muscle proliferation with a
tortuous, duplicated internal elastic lamina and an attenuated media — **not**
atheroma and not vasculitis. That is why the model's lesion is a growth process
(`PDGF -> SMC -> STEN`) with no inflammatory or lipid term, and why the only
drug that touches it is the statin's antiproliferative arm.

29. Masuda J, Ogata J, Yutani C Smooth muscle cell proliferation and localization of macrophages and T cells in the occlusive intracranial major arteries in moyamoya disease. *Stroke* 1993. PMID [7902623](https://pubmed.ncbi.nlm.nih.gov/7902623/)
30. Takekawa Y, Umezawa T, Ueno Y et al. Pathological and immunohistochemical findings of an autopsy case of adult moyamoya disease. *Neuropathology* 2004. PMID [15484702](https://pubmed.ncbi.nlm.nih.gov/15484702/)
31. Takagi Y, Kikuta K, Sadamasa N et al. Caspase-3-dependent apoptosis in middle cerebral arteries in patients with moyamoya disease. *Neurosurgery* 2006. PMID [17038954](https://pubmed.ncbi.nlm.nih.gov/17038954/)
32. Takagi Y, Kikuta K, Nozaki K et al. Expression of hypoxia-inducing factor-1 alpha and endoglin in intimal hyperplasia of the middle cerebral artery of patients with Moyamoya disease. *Neurosurgery* 2007. PMID [17290185](https://pubmed.ncbi.nlm.nih.gov/17290185/)
33. He S, Zhang J, Wang X et al. Organoid Modeling and Single-Cell Profiling Reveal Smooth Muscle Cell Migration in Moyamoya Disease. *Commun Biol* 2026. PMID [41501150](https://pubmed.ncbi.nlm.nih.gov/41501150/)
34. Asselman C, Meersschaut J, Willems P et al. Increased plasma fibronectin mirrors intimal phenotypic switching of vascular smooth muscle cells in moyamoya arteriopathy. *Sci Rep* 2025. PMID [41387995](https://pubmed.ncbi.nlm.nih.gov/41387995/)
35. Yamamoto S, Yamamoto S, Akai T et al. Differentiation of Fibroblasts Into Myofibroblasts in the Arachnoid Membrane of Moyamoya Disease. *Stroke* 2022. PMID [36039752](https://pubmed.ncbi.nlm.nih.gov/36039752/)
36. Kamata I, Terai Y, Ohmoto T Attempt to establish an experimental animal model of moyamoya disease using immuno-embolic material--histological changes of the arterial wall resulting from immunological reaction in cats. *Acta Med Okayama* 2003. PMID [12908012](https://pubmed.ncbi.nlm.nih.gov/12908012/)

## 4. Angiogenesis, collateral biology and the two collateral routes

The model's central anatomical claim is that there are **two** collateral
routes with opposite safety profiles: the leptomeningeal route donated by the
PCA (safe) and the periventricular/choroidal route (the bleeder). Hamano 2025
and Yamamoto 2019 (section 1) tie persistent periventricular anastomoses to
rebleeding; Mugikura 2026 and Tanaka 2026 tie PCA involvement to territorial
mismatch and to cognitive cost. Mukawa 2016 is the autopsy demonstration that an
indirect construct really does grow a new arterial network — the biological
content of `K_BYP_IND * VEGF * ANGIO`.

37. Phi JH, Suzuki N, Moon YJ et al. Chemokine Ligand 5 (CCL5) Derived from Endothelial Colony-Forming Cells (ECFCs) Mediates Recruitment of Smooth Muscle Progenitor Cells (SPCs) toward Critical Vascular Locations in Moyamoya Disease. *PLoS One* 2017. PMID [28072843](https://pubmed.ncbi.nlm.nih.gov/28072843/)
38. Cheng YW, Yang LY, Chen YT et al. Endothelial progenitor cell-derived conditioned medium mitigates chronic cerebral ischemic injury through macrophage migration inhibitory factor-activated AKT pathway. *Stem Cell Res Ther* 2024. PMID [39543689](https://pubmed.ncbi.nlm.nih.gov/39543689/)
39. Gupta T, Bharti R, Kumar M et al. Moyamoya disease and angiogenesis: a quantitative analysis of key angiogenic markers. *Neuroscience* 2026. PMID [42323005](https://pubmed.ncbi.nlm.nih.gov/42323005/)
40. Li Z, Zhao L, Yan P et al. Collateral circulation after revascularization in moyamoya disease: influencing factors and underlying mechanisms. *Front Neurol* 2026. PMID [42518953](https://pubmed.ncbi.nlm.nih.gov/42518953/)
41. Park GH, Shin HS, Choi ES et al. Cranial burr hole with erythropoietin administration induces reverse arteriogenesis from the enriched extracranium. *Neurobiol Dis* 2019. PMID [31344491](https://pubmed.ncbi.nlm.nih.gov/31344491/)
42. Mukawa M, Nariai T, Inaji M et al. First autopsy analysis of a neovascularized arterial network induced by indirect bypass surgery for moyamoya disease: case report. *J Neurosurg* 2016. PMID [26406800](https://pubmed.ncbi.nlm.nih.gov/26406800/)
43. Hamano E, Funaki T, Kataoka H et al. Persistent Periventricular Anastomosis Associated With Rebleeding After Bypass Surgery for Hemorrhagic Moyamoya Disease. *Stroke Vasc Interv Neurol* 2025. PMID [41573187](https://pubmed.ncbi.nlm.nih.gov/41573187/)
44. Mugikura S, Mori N Preoperative ischemic territorial mismatch due to PCA involvement, postoperative STA enlargement, and transient neurological events in moyamoya disease. *J Clin Neurosci* 2026. PMID [42537508](https://pubmed.ncbi.nlm.nih.gov/42537508/)
45. Mugikura S, Mori N Ischemic territorial mismatch due to posterior cerebral artery occlusion illustrated by the FLAIR ivy sign in moyamoya disease. *Jpn J Radiol* 2026. PMID [42159913](https://pubmed.ncbi.nlm.nih.gov/42159913/)
46. Tanaka K, Kusano Y, Funaki T et al. Neurocognitive profile of pediatric moyamoya disease with posterior cerebral artery involvement. *J Neurosurg Pediatr* 2026. PMID [42497452](https://pubmed.ncbi.nlm.nih.gov/42497452/)

## 5. Cerebral haemodynamics: autoregulation, reserve and the oxygen ceiling

This is where the model's arithmetic comes from. Lassen 1959 and Powers 1991
give autoregulation and the staged haemodynamic failure; **Derdeyn 2002** gives
the CBV/OEF staging the model's `OEF -> ceiling -> CMRO2 falls` sequence
reproduces; **Baron 1981** is the demonstration that misery perfusion is
*reversible by EC-IC bypass*, which is the physiological warrant for treating
`g_bypass` as a parallel conductance. Astrup, An 2015 and Guadagno supply the
penumbral threshold that the model does **not** fit but derives as
CMRO2/(CaO2 x OEF_max). Faraci 1987/1992 supply the large-artery share of
cerebrovascular resistance (`FRAC_PROX`); Caldwell 2021 supplies CO2 reactivity
(`K_CO2`); Ito 2006 supplies normal CBF.

47. LASSEN NA Cerebral blood flow and oxygen consumption in man. *Physiol Rev* 1959. PMID [13645234](https://pubmed.ncbi.nlm.nih.gov/13645234/)
48. Powers WJ Cerebral hemodynamics in ischemic cerebrovascular disease. *Ann Neurol* 1991. PMID [2042939](https://pubmed.ncbi.nlm.nih.gov/2042939/)
49. Derdeyn CP, Videen TO, Yundt KD et al. Variability of cerebral blood volume and oxygen extraction: stages of cerebral haemodynamic impairment revisited. *Brain* 2002. PMID [11872616](https://pubmed.ncbi.nlm.nih.gov/11872616/)
50. Baron JC, Bousser MG, Rey A et al. Reversal of focal "misery-perfusion syndrome" by extra-intracranial arterial bypass in hemodynamic cerebral ischemia. A case study with 15O positron emission tomography. *Stroke* 1981. PMID [6976022](https://pubmed.ncbi.nlm.nih.gov/6976022/)
51. Astrup J The Ischemic Penumbra 50 Years: A Personal Vignette. *Stroke* 2025. PMID [41144586](https://pubmed.ncbi.nlm.nih.gov/41144586/)
52. An H, Ford AL, Chen Y et al. Defining the ischemic penumbra using magnetic resonance oxygen metabolic index. *Stroke* 2015. PMID [25721017](https://pubmed.ncbi.nlm.nih.gov/25721017/)
53. Guadagno JV, Warburton EA, Aigbirhio FI et al. Does the acute diffusion-weighted imaging lesion represent penumbra as well as core? A combined quantitative PET/MRI voxel-based study. *J Cereb Blood Flow Metab* 2004. PMID [15545920](https://pubmed.ncbi.nlm.nih.gov/15545920/)
54. Ito H, Inoue K, Goto R et al. Database of normal human cerebral blood flow measured by SPECT: I. Comparison between I-123-IMP, Tc-99m-HMPAO, and Tc-99m-ECD as referred with O-15 labeled water PET and voxel-based morphometry. *Ann Nucl Med* 2006. PMID [16615422](https://pubmed.ncbi.nlm.nih.gov/16615422/)
55. Caldwell HG, Smith KJ, Lewis NCS et al. Regulation of cerebral blood flow by arterial PCO(2) independent of metabolic acidosis at 5050 m. *J Physiol* 2021. PMID [34047356](https://pubmed.ncbi.nlm.nih.gov/34047356/)
56. Faraci FM, Mayhan WG, Heistad DD Segmental vascular responses to acute hypertension in cerebrum and brain stem. *Am J Physiol* 1987. PMID [3565591](https://pubmed.ncbi.nlm.nih.gov/3565591/)
57. Faraci FM, Heistad DD Endothelium-derived relaxing factor inhibits constrictor responses of large cerebral arteries to serotonin. *J Cereb Blood Flow Metab* 1992. PMID [1569143](https://pubmed.ncbi.nlm.nih.gov/1569143/)
58. Webb KL, Mason CE, Prink JD et al. Impaired cerebral autoregulation in patients with moyamoya disease and moyamoya syndrome: A prospective cross-sectional study. *Clin Neurol Neurosurg* 2026. PMID [42155304](https://pubmed.ncbi.nlm.nih.gov/42155304/)
59. Chen S, Yu T, Xiang Y et al. Advances in cerebrovascular reserve function assessment for moyamoya disease. *BMJ Neurol Open* 2026. PMID [42465816](https://pubmed.ncbi.nlm.nih.gov/42465816/)

## 6. Diagnostic probes: acetazolamide, CO2, BOLD, and steal

The model predicts that acetazolamide and PaCO2 are **not interchangeable
probes** in moyamoya, because acetazolamide acts on the arteriole while PaCO2
also acts on the pial collateral conduits. Shi 2026 documents persistent steal
physiology on a graded reactivity challenge; Sebok 2025 shows pre-operative BOLD
reactivity predicting intraoperative flow; Park 2026 documents haemodynamic
change in the *unoperated* hemisphere, which is the cross-territory coupling the
three-node network makes explicit.

60. Uwano I, Kobayashi M, Setta K et al. Assessment of Impaired Cerebrovascular Reactivity in Chronic Cerebral Ischemia using Intravoxel Incoherent Motion Magnetic Resonance Imaging. *J Stroke Cerebrovasc Dis* 2021. PMID [34562793](https://pubmed.ncbi.nlm.nih.gov/34562793/)
61. Sebök M, Stumpo V, Bellomo J et al. Preoperative BOLD cerebrovascular reactivity correlates with intraoperative STA-MCA bypass flow and influences postoperative CVR improvement. *Eur Stroke J* 2025. PMID [40347485](https://pubmed.ncbi.nlm.nih.gov/40347485/)
62. Shi RB, Poublanc J, Sayin ES et al. Persistent Steal Physiology during a Ramp Cerebrovascular Reactivity Assessment Correlates with Severe Cortical Thinning in Steno-Occlusive Disease. *AJNR Am J Neuroradiol* 2026. PMID [41260673](https://pubmed.ncbi.nlm.nih.gov/41260673/)
63. McKetton L, Venkatraghavan L, Poublanc J et al. Importance of Collateralization in Patients With Large Artery Intracranial Occlusive Disease: Long-Term Longitudinal Assessment of Cerebral Hemodynamic Function. *Front Neurol* 2018. PMID [29681886](https://pubmed.ncbi.nlm.nih.gov/29681886/)
64. Shulgina AA, Lukshin VA, Korshunov AA et al. New Classification of the Degree of Cerebrovascular Insufficiency in Patients with Moyamoya Disease Measured According to ASL-MRI Perfusion. *Acta Neurochir Suppl* 2025. PMID [40632262](https://pubmed.ncbi.nlm.nih.gov/40632262/)
65. Yamada I, Matsushima Y, Suzuki S Moyamoya disease: diagnosis with three-dimensional time-of-flight MR angiography. *Radiology* 1992. PMID [1509066](https://pubmed.ncbi.nlm.nih.gov/1509066/)
66. Park TY, Lee SH, Chung Y et al. Hemodynamic Changes in Contralateral Unoperated Hemispheres Following Unilateral Combined Bypass Surgery in Adult Patients With Moyamoya Disease. *Neurosurgery* 2026. PMID [42132407](https://pubmed.ncbi.nlm.nih.gov/42132407/)
67. Kimura K, Akamatsu Y, Fujimoto K et al. Susceptibility Changes on Preoperative Acetazolamide-Loaded 7T MR Quantitative Susceptibility Mapping Predict Post-Carotid Endarterectomy Cerebral Hyperperfusion. *AJNR Am J Neuroradiol* 2025. PMID [39909571](https://pubmed.ncbi.nlm.nih.gov/39909571/)

## 7. Medical therapy, precipitants and quasi-moyamoya

**An evidence gap, stated as such.** There is no randomised trial of
antiplatelet therapy, statins, calcium-channel blockers or minocycline in
moyamoya disease; the literature is registry data, cohorts and case reports.
The model's medical arms are therefore *mechanistic hypotheses*, not
reproductions of trial results, and they are deliberately built so that every
one of them moves a **hazard** (embolic, wall repair, lesion growth rate) and
none of them moves `gS` by more than a few percent. Yang 2026 is included
because it records the EEG response to hyperventilation in moyamoya — the
laboratory form of the crying child. Santoro 2026 and Hamano 2025 cover the
Down-syndrome and autoimmune-thyroid associations that enter as `SEC`.

68. Yang H, Zhu Y, Lu J et al. EEG response to hyperventilation in patients with neurofibromatosis type 1 and Moyamoya syndrome: Two case reports. *Neurophysiol Clin* 2026. PMID [41825152](https://pubmed.ncbi.nlm.nih.gov/41825152/)
69. Rifino N, Aamodt AH, Wiedmann M et al. The Spectrum of Headaches in Moyamoya Angiopathy: From Mechanisms to Management Strategies-A Consensus Review From the NEUROVASC Working Group. *Eur J Neurol* 2025. PMID [41039799](https://pubmed.ncbi.nlm.nih.gov/41039799/)
70. Patterson AA, Vogel SKD, Sharp A et al. Cilostazol in a Child With Moyamoya Disease. *J Child Neurol* 2026. PMID [41505343](https://pubmed.ncbi.nlm.nih.gov/41505343/)
71. Alkhanafsa M, Wafi J, Aloqaily M et al. Adult-onset moyamoya disease presenting with recurrent ischemic stroke: The role of CT Perfusion and digital subtraction angiography in bypass planning. *Radiol Case Rep* 2026. PMID [42518706](https://pubmed.ncbi.nlm.nih.gov/42518706/)
72. Santoro JD, Silverman M, Wang AC et al. Delayed Recognition and Stroke-Predominant Presentation in Down Syndrome-Associated Moyamoya Syndrome. *Stroke* 2026. PMID [42495733](https://pubmed.ncbi.nlm.nih.gov/42495733/)
73. Hamano E, Kataoka H [Association between Moyamoya Disease and Autoimmune Thyroid Disorders]. *No Shinkei Geka* 2025. PMID [40438021](https://pubmed.ncbi.nlm.nih.gov/40438021/)

## 8. Surgical revascularisation: direct, indirect and combined

**Miyamoto 2014 is the JAM Trial**, and it is used asymmetrically: only the
*conservative* arm's 5-year rebleeding rate (31.6%) was used to set `HEM_HAZ0`,
so the bypass hazard ratio the model produces is a prediction. Nguyen 2022 and
Musmar 2026 supply the direct-versus-indirect comparison the model reproduces
from a single parameter (`ANGIO`); El-Hajj 2026 supplies the paediatric-versus-
adult contrast; **Shiino 2025** reports that presenting cerebrovascular
reactivity determines which operation succeeds, which is the model's selection
claim; Wang 2026 reports that EDAS reduces rebleeding in adult haemorrhagic
disease, i.e. that the haemorrhage benefit does not require a direct
anastomosis.

74. Miyamoto S, Yoshimoto T, Hashimoto N et al. Effects of extracranial-intracranial bypass for patients with hemorrhagic moyamoya disease: results of the Japan Adult Moyamoya Trial. *Stroke* 2014. PMID [24668203](https://pubmed.ncbi.nlm.nih.gov/24668203/)
75. Hishikawa T, Date I [Evidence of Efficacy of Superficial Temporal Artery-Middle Cerebral Artery Bypass in Japan]. *No Shinkei Geka* 2022. PMID [35946362](https://pubmed.ncbi.nlm.nih.gov/35946362/)
76. Nguyen VN, Motiwala M, Elarjani T et al. Direct, Indirect, and Combined Extracranial-to-Intracranial Bypass for Adult Moyamoya Disease: An Updated Systematic Review and Meta-Analysis. *Stroke* 2022. PMID [36134563](https://pubmed.ncbi.nlm.nih.gov/36134563/)
77. Musmar B, Roy JM, Abdalrazeq H et al. Direct Versus Indirect Bypass in Early-Stage Moyamoya (Suzuki I-III): A Propensity Score-Weighted Study. *Transl Stroke Res* 2026. PMID [42138779](https://pubmed.ncbi.nlm.nih.gov/42138779/)
78. Mertens R, Efe I, Mrosk F et al. Clinical Efficacy of Revascularization Surgery for Moyamoya Angiopathy: Long-Term Results of a European Cohort. *Eur J Neurol* 2026. PMID [42286416](https://pubmed.ncbi.nlm.nih.gov/42286416/)
79. de Liyis BG, Benet A, Kusdiansah M et al. Risk of post-revascularization stroke in Moyamoya disease: A systematic review, meta-analysis, and meta-regression. *Neurosurg Rev* 2026. PMID [41888499](https://pubmed.ncbi.nlm.nih.gov/41888499/)
80. Fujimura M, Tominaga T Lessons learned from moyamoya disease: outcome of direct/indirect revascularization surgery for 150 affected hemispheres. *Neurol Med Chir (Tokyo)* 2012. PMID [22688070](https://pubmed.ncbi.nlm.nih.gov/22688070/)
81. El-Hajj VG, Roy JM, Musmar B et al. Safety and long-term outcomes following bypass surgery in pediatric versus adult patients with Moyamoya disease: a multicenter cohort study. *Childs Nerv Syst* 2026. PMID [42026373](https://pubmed.ncbi.nlm.nih.gov/42026373/)
82. Shiino S, Han C, Garza M et al. Presenting cerebrovascular reactivity as a determinant of direct and indirect surgical revascularization success in North American patients with moyamoya vasculopathy. *J Neurosurg* 2025. PMID [39951721](https://pubmed.ncbi.nlm.nih.gov/39951721/)
83. Fujimura M, Ito M, Uchino H et al. Efficacy and Safety of Combined Revascularization Surgery for Moyamoya Disease: Standard Procedure and Perioperative Management. *Acta Neurochir Suppl* 2025. PMID [40632259](https://pubmed.ncbi.nlm.nih.gov/40632259/)
84. Kuroda S, Houkin K, Ishikawa T et al. Determinants of intellectual outcome after surgical revascularization in pediatric moyamoya disease: a multivariate analysis. *Childs Nerv Syst* 2004. PMID [15045517](https://pubmed.ncbi.nlm.nih.gov/15045517/)
85. Wang QN, Wang XP, Li XM et al. Encephaloduroarteriosynangiosis effectively reduces rebleeding risk in adults with hemorrhagic moyamoya disease: A long-term ambidirectional cohort study. *J Cereb Blood Flow Metab* 2026. PMID [42219999](https://pubmed.ncbi.nlm.nih.gov/42219999/)
86. Savoldi AM, Demartini Z Jr, Cordeiro ML Indirect Revascularization for Pediatric Moyamoya Angiopathy: Insights from a Brazilian Cohort. *J Clin Med* 2025. PMID [41227135](https://pubmed.ncbi.nlm.nih.gov/41227135/)

## 9. Post-operative cerebral hyperperfusion syndrome

The model's most heavily reworked claim. A lumped territory **cannot**
hyperperfuse, so the syndrome is represented as focal (peri-anastomotic) and
relative (against the flow the barrier had adapted to). Nishizawa 2020 and Jung
2025 supply the predictors and the day 2-7 timing; Kusdiansah 2026 supplies
early post-operative seizure as its clinical face; Chauhan 2026 supplies the
anaesthetic management that keeps PaCO2 out of the harmful range.

87. Nishizawa T, Fujimura M, Katsuki M et al. Prediction of Cerebral Hyperperfusion after Superficial Temporal Artery-Middle Cerebral Artery Anastomosis by Three-Dimensional-Time-of-Flight Magnetic Resonance Angiography in Adult Patients with Moyamoya Disease. *Cerebrovasc Dis* 2020. PMID [32829323](https://pubmed.ncbi.nlm.nih.gov/32829323/)
88. Jung MK, Ha EJ, Kim JH et al. Prediction of Cerebral Hyperperfusion Syndrome After Combined Bypass Surgery in Moyamoya Disease Using Hemodynamic and Clinical Data. *Clin Nucl Med* 2025. PMID [40173304](https://pubmed.ncbi.nlm.nih.gov/40173304/)
89. Kimata J, Tokairin K, Uchino H et al. Symptomatic cerebral hyperperfusion after occipital artery-posterior cerebral artery bypass in a patient with moyamoya disease: illustrative case. *J Neurosurg Case Lessons* 2025. PMID [40889378](https://pubmed.ncbi.nlm.nih.gov/40889378/)
90. Ito M, Uchino H, Fujimura M Reply to the Letter to the Editor: Intraoperative Cortical Indocyanine Green Extravasation as a Predictor of Cerebral Hyperperfusion following Direct Revascularization for Moyamoya Disease - Impact of Prolonged Observations of Indocyanine Green Videoangiography. *Cerebrovasc Dis* 2026. PMID [42030201](https://pubmed.ncbi.nlm.nih.gov/42030201/)
91. Kusdiansah M, de Liyis BG, Alhaq AMG et al. Early postoperative seizure after extracranial-intracranial revascularization for moyamoya disease and moyamoya syndrome: single centre post-revascularization moyamoya cohort study in Indonesia. *Neurosurg Rev* 2026. PMID [42481904](https://pubmed.ncbi.nlm.nih.gov/42481904/)
92. Hori S, Miyata Y, Shimohigoshi W et al. Thalamic anastomosis as a potential marker of perioperative ischemic risk in moyamoya disease. *Neurosurg Rev* 2026. PMID [42493721](https://pubmed.ncbi.nlm.nih.gov/42493721/)
93. Chauhan V Anesthetic Management for Encephaloduroarteriosynangiosis in Moyamoya Disease: A Hemodynamic and Neuromonitoring-Integrated Framework. *J Clin Med* 2026. PMID [42452416](https://pubmed.ncbi.nlm.nih.gov/42452416/)

## 10. Sickle-cell moyamoya and the transfusion threshold

**The cleanest test of the derived threshold.** In sickle-cell moyamoya the
vascular lesion is ordinary but CaO2 is halved, so CBF_crit moves from 19.7 to
37.0 mL/100g/min and the same flow becomes ischaemic. Adams 1998 (STOP) and
DeBaun 2014 (SIT) are the trials showing that **transfusion** — which touches no
vessel — is what works, which is precisely what the model predicts because
transfusion is the only intervention that moves the threshold rather than the
flow.

94. Adams RJ, McKie VC, Hsu L et al. Prevention of a first stroke by transfusions in children with sickle cell anemia and abnormal results on transcranial Doppler ultrasonography. *N Engl J Med* 1998. PMID [9647873](https://pubmed.ncbi.nlm.nih.gov/9647873/)
95. DeBaun MR, Casella JF Transfusions for silent cerebral infarcts in sickle cell anemia. *N Engl J Med* 2014. PMID [25372094](https://pubmed.ncbi.nlm.nih.gov/25372094/)
96. Richerson WT, Aumann MA, Song AK et al. Cerebral Blood Transit in Sickle Cell Anemia. *J Magn Reson Imaging* 2026. PMID [42210666](https://pubmed.ncbi.nlm.nih.gov/42210666/)
97. Papadakis JE, Archer NM, Singh N et al. Understanding stroke risk phenotypes in pediatric patients with sickle cell disease and concurrent moyamoya arteriopathy: insights from 61 cases at a single institution. *J Neurosurg Pediatr* 2026. PMID [41825074](https://pubmed.ncbi.nlm.nih.gov/41825074/)
98. Chen Y, Wang Y, Phuah CL et al. Toward Automated Detection of Silent Cerebral Infarcts in Children and Young Adults With Sickle Cell Anemia. *Stroke* 2023. PMID [37387218](https://pubmed.ncbi.nlm.nih.gov/37387218/)
99. Kwiatkowski JL, Voeks JH, Kanter J et al. Ischemic stroke in children and young adults with sickle cell disease in the post-STOP era. *Am J Hematol* 2019. PMID [31489983](https://pubmed.ncbi.nlm.nih.gov/31489983/)

---

## Deliberate omissions

- **No randomised evidence exists** for antiplatelet therapy, statins,
  calcium-channel blockers, minocycline or edaravone in moyamoya disease. The
  model's medical arms are mechanistic hypotheses and are labelled as such in
  `mmd_mrgsolve_model.R`. Padding this list with case reports would misrepresent
  the strength of that evidence.
- **The JET study** (Japanese EC/IC Bypass Trial) is represented by Hishikawa's
  Japanese-language review of the STA-MCA bypass evidence rather than by a primary
  English-language report, because a primary report indexed under that trial name
  could not be resolved with confidence.
- **Parameter values taken from physiology rather than from moyamoya papers**
  (normal CBF, CMRO2, OEF, CaO2, the large-artery resistance share, CO2
  reactivity) are cited in section 5 and are deliberately *not* moyamoya-specific:
  they are the healthy baseline the disease is a perturbation of.

## Disclaimer

This model is an educational and research artefact. It has not been validated
against patient-level data and must not be used for clinical decision-making.

