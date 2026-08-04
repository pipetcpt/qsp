# 악성 종양 관련 고칼슘혈증 — 참고문헌 (Hypercalcaemia of Malignancy — References)

모든 PMID는 NCBI E-utilities(`esearch` → `esummary`)로 조회하여 저자·저널·연도·제목을
확인했습니다. 각 항목 뒤의 이탤릭 문구는 **그 문헌이 모델의 어느 구조·파라미터를
뒷받침하는지**를 적은 것입니다. 근거가 없는 곳은 "가정(assumption)"이라고 명시했습니다.

All PMIDs were resolved and title-checked through NCBI E-utilities. The italic
note after each entry states which structure or parameter of the model the
reference supports; where the model is asserting rather than citing, the note
says so.

---

## 1. 개요 · 역학 · 진료지침 (Overview, epidemiology, guidelines)

1. Stewart AF. Clinical practice. Hypercalcemia associated with cancer. *N Engl J Med* 2005. [PMID 15673803](https://pubmed.ncbi.nlm.nih.gov/15673803/)
   *The four-mechanism taxonomy the model is organised around (humoral / osteolytic / calcitriol / ectopic PTH) and their approximate frequencies.*
2. El-Hajj Fuleihan G, Clines GA, Hu MI, et al. Treatment of Hypercalcemia of Malignancy in Adults: An Endocrine Society Clinical Practice Guideline. *J Clin Endocrinol Metab* 2023. [PMID 36545746](https://pubmed.ncbi.nlm.nih.gov/36545746/)
   *Therapy sequencing (saline → antiresorptive, calcitonin as a bridge, denosumab in bisphosphonate failure and renal impairment) against which the scenario set is benchmarked.*
3. Zagzag J, Hu MI, Fisher SB, Perrier ND. Hypercalcemia and cancer: differential diagnosis and treatment. *CA Cancer J Clin* 2018. [PMID 30240520](https://pubmed.ncbi.nlm.nih.gov/30240520/)
   *The diagnostic panel encoded in the model outputs: PTH, PTHrP, 1,25D, phosphate.*
4. Rosner MH, Dalkin AC. Onco-nephrology: the pathophysiology and treatment of malignancy-associated hypercalcemia. *Clin J Am Soc Nephrol* 2012. [PMID 22879438](https://pubmed.ncbi.nlm.nih.gov/22879438/)
   *Volume depletion → reduced GFR → reduced calcium clearance, the closed loop the model formalises as a saddle-node.*
5. Mundy GR, Guise TA. Hypercalcemia of malignancy. *Am J Med* 1997. [PMID 9274897](https://pubmed.ncbi.nlm.nih.gov/9274897/)
   *Bone resorption rates in malignancy of roughly 2–5× normal; the model's J_res range of 13–34 mmol/day against a baseline of 8.*
6. Ralston SH, Gallacher SJ, Patel U, Campbell J, Boyle IT. Cancer-associated hypercalcemia: morbidity and mortality. Clinical experience in 126 treated patients. *Ann Intern Med* 1990. [PMID 2138442](https://pubmed.ncbi.nlm.nih.gov/2138442/)
   *Presenting calcium distribution and the prognosis that motivates note J's warning against reading the collapse boundary as survival.*
7. Goldner W. Cancer-Related Hypercalcemia. *J Oncol Pract* 2016. [PMID 27170690](https://pubmed.ncbi.nlm.nih.gov/27170690/)
8. Asonitis N, Angelousi A, Zafeiris C, Lambrou GI, Dontas I, Kassi E. Diagnosis, pathophysiology and management of hypercalcemia in malignancy: a review. *Horm Metab Res* 2019. [PMID 31826272](https://pubmed.ncbi.nlm.nih.gov/31826272/)
9. Lindner G, Felber R, Schwarz C, et al. Hypercalcemia in the ED: prevalence, etiology, and outcome. *Am J Emerg Med* 2013. [PMID 23246111](https://pubmed.ncbi.nlm.nih.gov/23246111/)
10. Legrand SB. Modern management of malignant hypercalcemia. *Am J Hosp Palliat Care* 2011. [PMID 21724679](https://pubmed.ncbi.nlm.nih.gov/21724679/)
11. Sternlicht H, Glezerman IG. Hypercalcemia of malignancy and new treatment options. *Ther Clin Risk Manag* 2015. [PMID 26675713](https://pubmed.ncbi.nlm.nih.gov/26675713/)
12. Bilezikian JP. Management of acute hypercalcemia. *N Engl J Med* 1992. [PMID 1532633](https://pubmed.ncbi.nlm.nih.gov/1532633/)
    *Saline rates of 200–300 mL/h; the model's IV_RATE of 4.8 L/day.*
13. Ariyan CE, Sosa JA. Assessment and management of patients with abnormal calcium. *Crit Care Med* 2004. [PMID 15064673](https://pubmed.ncbi.nlm.nih.gov/15064673/)
14. Bushinsky DA, Monk RD. Electrolyte quintet: Calcium. *Lancet* 1998. [PMID 9690425](https://pubmed.ncbi.nlm.nih.gov/9690425/)
    *Whole-body calcium fluxes: skeletal pool ≈1 kg, ECF pool ≈36 mmol, daily turnover — the compartment sizes in `$MAIN`.*

---

## 2. PTHrP와 체액성 고칼슘혈증 (PTHrP and the humoral mechanism)

15. Suva LJ, Winslow GA, Wettenhall RE, et al. A parathyroid hormone-related protein implicated in malignant hypercalcemia: cloning and expression. *Science* 1987. [PMID 3616618](https://pubmed.ncbi.nlm.nih.gov/3616618/)
16. Burtis WJ, Brady TG, Orloff JJ, et al. Immunochemical characterization of circulating parathyroid hormone-related protein in patients with humoral hypercalcemia of cancer. *N Engl J Med* 1990. [PMID 2320080](https://pubmed.ncbi.nlm.nih.gov/2320080/)
    *Circulating PTHrP concentrations of a few to tens of pmol/L; the model's PTHRP trajectory reaches 14.3 pmol/L at presentation.*
17. Yates AJ, Gutierrez GE, Smolens P, et al. Effects of a synthetic peptide of a parathyroid hormone-related protein on calcium homeostasis, renal tubular calcium reabsorption, and bone metabolism in vivo. *J Clin Invest* 1988. [PMID 3343349](https://pubmed.ncbi.nlm.nih.gov/3343349/)
    *The single most important reference for THESIS 3: PTHrP simultaneously increases bone resorption AND renal tubular calcium reabsorption. This is the experimental basis for the model's PS-dependent distal Tm.*
18. Guise TA, Yin JJ, Taylor SD, et al. Evidence for a causal role of parathyroid hormone-related protein in the pathogenesis of human breast cancer-mediated osteolysis. *J Clin Invest* 1996. [PMID 8833902](https://pubmed.ncbi.nlm.nih.gov/8833902/)
19. Wysolmerski JJ, Broadus AE. Hypercalcemia of malignancy: the central role of parathyroid hormone-related protein. *Annu Rev Med* 1994. [PMID 8198376](https://pubmed.ncbi.nlm.nih.gov/8198376/)
20. Wysolmerski JJ. Parathyroid hormone-related protein: an update. *J Clin Endocrinol Metab* 2012. [PMID 22745236](https://pubmed.ncbi.nlm.nih.gov/22745236/)
21. Horwitz MJ, Tedesco MB, Sereika SM, et al. A 7-day continuous infusion of PTH or PTHrP suppresses bone formation and uncouples bone turnover. *J Bone Miner Res* 2011. [PMID 21544866](https://pubmed.ncbi.nlm.nih.gov/21544866/)
    *Direct human evidence for the model's uncoupling term: continuous (not intermittent) PTH1R signalling raises resorption while SUPPRESSING formation. Parameters `OB_PTH_EMAX` and `O_PTH_EMAX`.*
22. Fukumoto S, Matsumoto T, Yamoto H, et al. Suppression of serum 1,25-dihydroxyvitamin D in humoral hypercalcemia of malignancy. *Endocrinology* 1989. [PMID 2539966](https://pubmed.ncbi.nlm.nih.gov/2539966/)
    *Why the model gives PTHrP only 15% of PTH's potency at renal 1α-hydroxylase (`W1A_PTHRP` = 0.15): 1,25D is LOW in humoral disease despite a high PTH1R signal in bone.*
23. Nakayama K, Fukumoto S, Takeda S, et al. Differences in bone and vitamin D metabolism between primary hyperparathyroidism and malignancy-associated hypercalcemia. *J Clin Endocrinol Metab* 1996. [PMID 8636276](https://pubmed.ncbi.nlm.nih.gov/8636276/)
    *The paired comparison the model reproduces in note D: same PTH1R signal, opposite 1,25D and opposite bone-formation markers.*
24. Kremer R, Shustik C, Tabak T, Papavasiliou V, Goltzman D. Parathyroid-hormone-related peptide in hematologic malignancies. *Am J Med* 1996. [PMID 8610726](https://pubmed.ncbi.nlm.nih.gov/8610726/)
25. Ratcliffe WA, Hutchesson AC, Bundred NJ, Ratcliffe JG. Role of assays for parathyroid-hormone-related protein in investigation of hypercalcaemia. *Lancet* 1992. [PMID 1346019](https://pubmed.ncbi.nlm.nih.gov/1346019/)
26. Pecherstorfer M, Schilling T, Blind E, et al. Parathyroid hormone-related protein and life expectancy in hypercalcemic cancer patients. *J Clin Endocrinol Metab* 1994. [PMID 8175989](https://pubmed.ncbi.nlm.nih.gov/8175989/)
27. Gurney H, Grill V, Martin TJ. Parathyroid hormone-related protein and response to pamidronate in tumour-induced hypercalcaemia. *Lancet* 1993. [PMID 8099988](https://pubmed.ncbi.nlm.nih.gov/8099988/)
    *Clinical confirmation of the model's central therapeutic prediction: a high PTHrP predicts a POORER bisphosphonate response at the same calcium (scenarios 5 vs 11).*
28. Wimalawansa SJ. Significance of plasma PTH-rp in patients with hypercalcemia of malignancy treated with bisphosphonate. *Cancer* 1994. [PMID 8156530](https://pubmed.ncbi.nlm.nih.gov/8156530/)
29. Walls J, Ratcliffe WA, Howell A, Bundred NJ. Response to intravenous bisphosphonate therapy in hypercalcaemic patients with and without bone metastases: the role of parathyroid hormone-related protein. *Br J Cancer* 1994. [PMID 8018531](https://pubmed.ncbi.nlm.nih.gov/8018531/)

---

## 3. 국소 골용해성 기전 (Local osteolytic mechanism)

30. Roodman GD. Mechanisms of bone metastasis. *N Engl J Med* 2004. [PMID 15084698](https://pubmed.ncbi.nlm.nih.gov/15084698/)
31. Roodman GD. MIP-1 alpha and myeloma bone disease. *Cancer Treat Res* 2004. [PMID 15043189](https://pubmed.ncbi.nlm.nih.gov/15043189/)
    *The osteoclast-precursor chemoattractant arm; the model's `OCP_CYT_E` term.*
32. Tian E, Zhan F, Walker R, et al. The role of the Wnt-signaling antagonist DKK1 in the development of osteolytic lesions in multiple myeloma. *N Engl J Med* 2003. [PMID 14695408](https://pubmed.ncbi.nlm.nih.gov/14695408/)
    *The osteoblast-suppression arm of the osteolytic mechanism. NOTE: the current model represents myeloma uncoupling only through the shared PTH1R/TGF-β coupling terms; an explicit DKK1 state is a documented omission.*

---

## 4. 칼시트리올 매개 기전 (Calcitriol-mediated mechanism)

33. Seymour JF, Gagel RF. Calcitriol: the major humoral mediator of hypercalcemia in Hodgkin's disease and non-Hodgkin's lymphomas. *Blood* 1993. [PMID 8364192](https://pubmed.ncbi.nlm.nih.gov/8364192/)
34. Adams JS, Singer FR, Gacad MA, et al. Isolation and structural identification of 1,25-dihydroxyvitamin D3 produced by cultured alveolar macrophages in sarcoidosis. *J Clin Endocrinol Metab* 1985. [PMID 2984238](https://pubmed.ncbi.nlm.nih.gov/2984238/)
35. Barbour GL, Coburn JW, Slatopolsky E, Norman AW, Horst RL. Hypercalcemia in an anephric patient with sarcoidosis: evidence for extrarenal generation of 1,25-dihydroxyvitamin D. *N Engl J Med* 1981. [PMID 6894783](https://pubmed.ncbi.nlm.nih.gov/6894783/)
    *The classic demonstration that the extrarenal enzyme is outside the FGF23/calcium feedback loop entirely — the model's `ext_1a` term has no negative feedback of any kind, deliberately.*
36. Hewison M, Kantorovich V, Liker HR, et al. Vitamin D-mediated hypercalcemia in lymphoma: evidence for hormone production by tumor-adjacent macrophages. *J Bone Miner Res* 2003. [PMID 12619944](https://pubmed.ncbi.nlm.nih.gov/12619944/)
    *Also the basis for the substrate-consumption term `D25_USE`: 25(OH)D runs low because it is being consumed, which is what bounds the calcitriol mechanism in the model.*

---

## 5. 이소성 PTH · 부갑상선암 · 칼슘유사체 (Ectopic PTH, parathyroid carcinoma, calcimimetics)

37. Silverberg SJ, Rubin MR, Faiman C, et al. Cinacalcet hydrochloride reduces the serum calcium concentration in inoperable parathyroid carcinoma. *J Clin Endocrinol Metab* 2007. [PMID 17666472](https://pubmed.ncbi.nlm.nih.gov/17666472/)
    *Scenario 14; the magnitude of the calcium fall on 90 mg bd is the anchor for `CIN_EMAX`.*
38. Peacock M, Bilezikian JP, Klassen PS, Guo MD, Turner SA, Shoback D. Cinacalcet hydrochloride maintains long-term normocalcemia in patients with primary hyperparathyroidism. *J Clin Endocrinol Metab* 2005. [PMID 15522938](https://pubmed.ncbi.nlm.nih.gov/15522938/)
39. Nemeth EF, Shoback D. Calcimimetic and calcilytic drugs for treating bone and mineral-related disorders. *Best Pract Res Clin Endocrinol Metab* 2013. [PMID 23856266](https://pubmed.ncbi.nlm.nih.gov/23856266/)
    *Mechanism of the set-point left shift the model implements as `PTH_SET × (1 − E_cin)`.*

---

## 6. 골 리모델링 · RANKL/OPG (Bone remodelling and the RANKL/OPG axis)

40. Simonet WS, Lacey DL, Dunstan CR, et al. Osteoprotegerin: a novel secreted protein involved in the regulation of bone density. *Cell* 1997. [PMID 9108485](https://pubmed.ncbi.nlm.nih.gov/9108485/)
41. Yasuda H, Shima N, Nakagawa N, et al. A novel molecular mechanism modulating osteoclast differentiation and function. *Bone* 1999. [PMID 10423033](https://pubmed.ncbi.nlm.nih.gov/10423033/)
42. Boyle WJ, Simonet WS, Lacey DL. Osteoclast differentiation and activation. *Nature* 2003. [PMID 12748652](https://pubmed.ncbi.nlm.nih.gov/12748652/)
    *The RANKL/RANK/OPG structure of `$ODE` section 21, including PTH's opposite effects on the two ligands.*
43. Sims NA, Martin TJ. Coupling the activities of bone formation and resorption: a multitude of signals within the basic multicellular unit. *Bonekey Rep* 2014. [PMID 24466412](https://pubmed.ncbi.nlm.nih.gov/24466412/)
    *The TGF-β coupling term (`TGF_K`, `OBP_TGF_EMAX`) that makes antiresorptive therapy also suppress formation.*
44. Everts V, Delaissé JM, Korper W, et al. The bone lining cell: its role in cleaning Howship's lacunae and initiating bone formation. *J Bone Miner Res* 2002. [PMID 11771672](https://pubmed.ncbi.nlm.nih.gov/11771672/)
45. Garnero P, Ferreras M, Karsdal MA, et al. The type I collagen fragments ICTP and CTX reveal distinct enzymatic pathways of bone collagen degradation. *J Bone Miner Res* 2003. [PMID 12733725](https://pubmed.ncbi.nlm.nih.gov/12733725/)
    *CTX as the cathepsin-K-dependent resorption readout; the model's CTX compartment is a first-order lag on J_res.*
46. Vasikaran S, Eastell R, Bruyère O, et al. Markers of bone turnover for the prediction of fracture risk and monitoring of osteoporosis treatment. *Osteoporos Int* 2011. [PMID 21184054](https://pubmed.ncbi.nlm.nih.gov/21184054/)
    *Reference intervals for CTX and P1NP used as the healthy initial conditions (0.35 ng/mL and 45 µg/L).*

---

## 7. CaSR · 신장 칼슘 처리 (Calcium-sensing receptor and renal handling)

47. Riccardi D, Brown EM. Physiology and pathophysiology of the calcium-sensing receptor in the kidney. *Am J Physiol Renal Physiol* 2010. [PMID 19923405](https://pubmed.ncbi.nlm.nih.gov/19923405/)
48. Hebert SC. Extracellular calcium-sensing receptor: implications for calcium and magnesium handling in the kidney. *Kidney Int* 1996. [PMID 8943500](https://pubmed.ncbi.nlm.nih.gov/8943500/)
    *The TAL CaSR term (`CASR_EMAX`, `CASR_EC50`, `CASR_H`) and the coupled magnesium wasting.*
49. Loupy A, Ramakrishnan SK, Wootla B, et al. PTH-independent regulation of blood calcium concentration by the calcium-sensing receptor. *J Clin Invest* 2012. [PMID 22886306](https://pubmed.ncbi.nlm.nih.gov/22886306/)
    *Evidence that the renal CaSR defends plasma calcium independently of PTH — which is why the model's calciuresis survives complete PTH suppression.*
50. Hoenderop JG, Nilius B, Bindels RJ. Calcium absorption across epithelia. *Physiol Rev* 2005. [PMID 15618484](https://pubmed.ncbi.nlm.nih.gov/15618484/)
    *TRPV5/TRPV6, calbindin, NCX1/PMCA1b; the basis for treating the distal segment as TRANSPORT-limited (a Tm) rather than fraction-limited.*
51. Hoenderop JG, van Leeuwen JP, van der Eerden BC, et al. Renal Ca2+ wasting, hyperabsorption, and reduced bone thickness in mice lacking TRPV5. *J Clin Invest* 2003. [PMID 14679186](https://pubmed.ncbi.nlm.nih.gov/14679186/)
52. Alexander RT, Dimke H, Cordat E. Proximal tubular NHEs: sodium, protons and calcium? *Am J Physiol Renal Physiol* 2013. [PMID 23761670](https://pubmed.ncbi.nlm.nih.gov/23761670/)
    *Sodium-coupled paracellular proximal calcium reabsorption — the mechanism behind `AVID_P`, and therefore behind why saline is calciuretic at all.*
53. Blaine J, Chonchol M, Levi M. Renal control of calcium, phosphate, and magnesium homeostasis. *Clin J Am Soc Nephrol* 2015. [PMID 26384363](https://pubmed.ncbi.nlm.nih.gov/26384363/)
    *Segmental fractional reabsorption (≈65% proximal, ≈25% TAL, ≈8% distal, FE ≈1–2%) reproduced exactly by the healthy steady state.*
54. Peacock M, Robertson WG, Nordin BE. Relation between serum and urinary calcium with particular reference to parathyroid activity. *Lancet* 1969. [PMID 4179224](https://pubmed.ncbi.nlm.nih.gov/4179224/)
    *The classic serum-calcium/urine-calcium relation whose slope and PTH dependence the model's `U_Ca(Ca_tot)` curve has to reproduce.*

---

## 8. 신성 요붕증 · 농축능 장애 · 고칼슘 신병증 (Nephrogenic DI, concentrating defect, hypercalcaemic nephropathy)

55. Sands JM, Naruse M, Baum M, et al. Apical extracellular calcium/polyvalent cation-sensing receptor regulates vasopressin-elicited water permeability in rat kidney inner medullary collecting duct. *J Clin Invest* 1997. [PMID 9077550](https://pubmed.ncbi.nlm.nih.gov/9077550/)
56. Earm JH, Christensen BM, Frøkiaer J, et al. Decreased aquaporin-2 expression and apical plasma membrane delivery in kidney collecting ducts of polyuric hypercalcemic rats. *J Am Soc Nephrol* 1998. [PMID 9848772](https://pubmed.ncbi.nlm.nih.gov/9848772/)
    *The AQP2 state variable and its ~2-day time constant in both directions — the "memory" that makes relapse a volume event before it is a bone event.*
57. Levi M, Peterson L, Berl T. Mechanism of concentrating defect in hypercalcemia. Role of polydipsia and prostaglandins. *Kidney Int* 1983. [PMID 6573545](https://pubmed.ncbi.nlm.nih.gov/6573545/)
58. Gill JR Jr, Bartter FC. On the impairment of renal concentrating ability in prolonged hypercalcemia and hypercalciuria in man. *J Clin Invest* 1961. [PMID 13705309](https://pubmed.ncbi.nlm.nih.gov/13705309/)
    *Maximal urine osmolality falling to the 300–400 mosm/kg range; the model's `UOSM_MAX × AQP2` ceiling.*
59. Fenton RA, Knepper MA. Mouse models and the urinary concentrating mechanism in the new millennium. *Physiol Rev* 2007. [PMID 17928581](https://pubmed.ncbi.nlm.nih.gov/17928581/)
60. Benabe JE, Martinez-Maldonado M. Hypercalcemic nephropathy. *Arch Intern Med* 1978. [PMID 646542](https://pubmed.ncbi.nlm.nih.gov/646542/)
61. Moysés-Neto M, Guimarães FM, Ayoub FH, Vieira-Neto OM, Costa JA, Dantas M. Acute renal failure and hypercalcemia. *Ren Fail* 2006. [PMID 16538974](https://pubmed.ncbi.nlm.nih.gov/16538974/)
    *The direct-injury arm (`W_ICA`) alongside the prerenal arm.*

---

## 9. 칼슘 화학종 · 측정 (Speciation and measurement)

62. Payne RB, Little AJ, Williams RB, Milner JR. Interpretation of serum calcium in patients with abnormal serum proteins. *Br Med J* 1973. [PMID 4758544](https://pubmed.ncbi.nlm.nih.gov/4758544/)
    *The correction itself. The model computes it as `Ca_corr` purely so that note G can show where it fails.*
63. Baird GS. Ionized calcium. *Clin Chim Acta* 2011. [PMID 21238441](https://pubmed.ncbi.nlm.nih.gov/21238441/)
    *pH dependence of protein binding (≈0.05 mmol/L per 0.1 pH unit) — the model's `PH_GAMMA` = 0.41 reproduces 0.063 mmol/L per 0.1 unit at a total of 3.0.*
64. Gauci C, Moranne O, Fouqueray B, et al. Pitfalls of measuring total blood calcium in patients with CKD. *J Am Soc Nephrol* 2008. [PMID 18400941](https://pubmed.ncbi.nlm.nih.gov/18400941/)
    *Independent demonstration that albumin-based correction misclassifies; the model's version of this is quantitative and mechanism-based rather than regression-based.*

---

## 10. 수액 · 이뇨제 (Saline and diuretics)

65. Suki WN, Yium JJ, Von Minden M, Saller-Hebert C, Eknoyan G, Martinez-Maldonado M. Acute treatment of hypercalcemia with furosemide. *N Engl J Med* 1970. [PMID 5458033](https://pubmed.ncbi.nlm.nih.gov/5458033/)
    *The original report, in which furosemide was given with meticulous volume replacement — which is exactly the condition under which the model also finds it useful (scenario 15).*
66. LeGrand SB, Leskuski D, Zama I. Narrative review: furosemide for hypercalcemia: an unproven yet common practice. *Ann Intern Med* 2008. [PMID 18711156](https://pubmed.ncbi.nlm.nih.gov/18711156/)
    *The systematic case against. The model resolves the disagreement by making the volume underneath the diuretic the deciding variable (scenarios 15 vs 16 vs 23).*
67. Wermers RA, Kearns AE, Jenkins GD, Melton LJ 3rd. Incidence and clinical spectrum of thiazide-associated hypercalcemia. *Am J Med* 2007. [PMID 17904464](https://pubmed.ncbi.nlm.nih.gov/17904464/)
    *Scenario 20 and the 13% reduction in the fold reported in note B.*

---

## 11. 칼시토닌 (Calcitonin)

68. Wisneski LA. Salmon calcitonin in the acute management of hypercalcemia. *Calcif Tissue Int* 1990. [PMID 2137363](https://pubmed.ncbi.nlm.nih.gov/2137363/)
    *Onset within 2–6 h and a fall of roughly 0.3–0.5 mmol/L; the model gives 0.19 mmol/L at 6 h over and above saline.*
69. Binstock ML, Mundy GR. Effect of calcitonin and glucocorticoids in combination on the hypercalcemia of malignancy. *Ann Intern Med* 1980. [PMID 7406378](https://pubmed.ncbi.nlm.nih.gov/7406378/)
70. Ralston SH, Alzaid AA, Gardner MD, Boyle IT. Treatment of cancer associated hypercalcaemia with combined aminohydroxypropylidene diphosphonate and calcitonin. *Br Med J (Clin Res Ed)* 1986. [PMID 3087513](https://pubmed.ncbi.nlm.nih.gov/3087513/)
    *The bridge combination itself; note C is the model's account of why it is a relay rather than a synergy.*
71. Ralston SH, Gardner MD, Dryburgh FJ, Jenkins AS, Cowan RA, Boyle IT. Comparison of aminohydroxypropylidene diphosphonate, mithramycin, and corticosteroids/calcitonin in treatment of cancer-associated hypercalcaemia. *Lancet* 1985. [PMID 2865417](https://pubmed.ncbi.nlm.nih.gov/2865417/)
72. Takahashi S, Goldring S, Katz M, Hilsenbeck S, Williams R, Roodman GD. Downregulation of calcitonin receptor mRNA expression by calcitonin during human osteoclast-like cell differentiation. *J Clin Invest* 1995. [PMID 7814611](https://pubmed.ncbi.nlm.nih.gov/7814611/)
    *The molecular basis of tachyphylaxis and therefore of the `R_CT` state and `CT_KOUT`/`CT_KIN`.*
73. Wada S, Udagawa N, Akatsu T, Nagata N, Martin TJ, Findlay DM. Regulation by calcitonin and glucocorticoids of calcitonin receptor gene expression in mouse osteoclasts. *Endocrinology* 1997. [PMID 9002981](https://pubmed.ncbi.nlm.nih.gov/9002981/)
74. Chambers TJ, Magnus CJ. Calcitonin alters behaviour of isolated osteoclasts. *J Pathol* 1982. [PMID 7057295](https://pubmed.ncbi.nlm.nih.gov/7057295/)
    *Ruffled-border retraction within minutes — the reason the calcitonin effect is modelled on osteoclast ACTIVITY rather than number.*
75. Hirsch PF, Baruch H. Is calcitonin an important physiological substance? *Endocrine* 2003. [PMID 14515002](https://pubmed.ncbi.nlm.nih.gov/14515002/)
76. Deftos LJ, Neer R. Medical management of the hypercalcemia of malignancy. *Annu Rev Med* 1974. [PMID 4277568](https://pubmed.ncbi.nlm.nih.gov/4277568/)

---

## 12. 비스포스포네이트 (Bisphosphonates)

77. Major P, Lortholary A, Hon J, et al. Zoledronic acid is superior to pamidronate in the treatment of hypercalcemia of malignancy: a pooled analysis of two randomized, controlled clinical trials. *J Clin Oncol* 2001. [PMID 11208851](https://pubmed.ncbi.nlm.nih.gov/11208851/)
    *The primary efficacy anchor: normalisation by day 4–10 in the large majority, and the notably FLAT 4 mg → 8 mg increment that note F reproduces from the delivery step.*
78. Gucalp R, Ritch P, Wiernik PH, et al. Comparative study of pamidronate disodium and etidronate disodium in the treatment of cancer-related hypercalcemia. *J Clin Oncol* 1992. [PMID 1727915](https://pubmed.ncbi.nlm.nih.gov/1727915/)
79. Nussbaum SR, Younger J, Vandepol CJ, et al. Single-dose intravenous therapy with pamidronate for the treatment of hypercalcemia of malignancy: comparison of 30-, 60-, and 90-mg dosages. *Am J Med* 1993. [PMID 8368227](https://pubmed.ncbi.nlm.nih.gov/8368227/)
    *A second, independent flat-topped dose-response.*
80. Rogers MJ, Crockett JC, Coxon FP, Mönkkönen J. Biochemical and molecular mechanisms of action of bisphosphonates. *Bone* 2011. [PMID 21111853](https://pubmed.ncbi.nlm.nih.gov/21111853/)
    *FPP synthase inhibition → loss of prenylation → osteoclast apoptosis; the `ZOL_EMAX`/`ZOL_IMAX` pair.*
81. Sato M, Grasser W, Endo N, et al. Bisphosphonate action. Alendronate localization in rat bone and effects on osteoclast ultrastructure. *J Clin Invest* 1991. [PMID 1661297](https://pubmed.ncbi.nlm.nih.gov/1661297/)
    *The physical basis of the self-delivery term: drug concentrates under the resorbing osteoclast.*
82. Coxon FP, Thompson K, Rogers MJ. Recent advances in understanding the mechanism of action of bisphosphonates. *Curr Opin Pharmacol* 2006. [PMID 16650801](https://pubmed.ncbi.nlm.nih.gov/16650801/)
    *Uptake requires active resorption — the sentence the model turns into `uptake = ZOL_KUP × ZOL_B × (J_res/J_RES0)`.*
83. Cremers S, Papapoulos S. Pharmacology of bisphosphonates. *Bone* 2011. [PMID 21281748](https://pubmed.ncbi.nlm.nih.gov/21281748/)
    *Skeletal retention of roughly 40–60% of an intravenous dose, and the multi-year terminal half-life produced in the model by the buried-in-matrix compartment.*
84. Chen T, Berenson J, Vescio R, et al. Pharmacokinetics and pharmacodynamics of zoledronic acid in cancer patients with bone metastases. *J Clin Pharmacol* 2002. [PMID 12412821](https://pubmed.ncbi.nlm.nih.gov/12412821/)
    *Plasma disposition and renal clearance; `ZOL_KEL` scales with GFR for this reason.*
85. Rosen LS, Gordon D, Kaminski M, et al. Long-term efficacy and safety of zoledronic acid compared with pamidronate disodium in the treatment of skeletal complications in patients with advanced multiple myeloma or breast carcinoma. *Cancer* 2003. [PMID 14534891](https://pubmed.ncbi.nlm.nih.gov/14534891/)
86. Markowitz GS, Fine PL, Stack JI, et al. Toxic acute tubular necrosis following treatment with zoledronate (Zometa). *Kidney Int* 2003. [PMID 12787420](https://pubmed.ncbi.nlm.nih.gov/12787420/)
    *The toxicity that makes the renal-impairment comparison in note H clinically consequential. NOT implemented as an ODE — a documented omission.*
87. Chennuru S, Koduri J, Baumann MA. Risk factors for symptomatic hypocalcaemia complicating treatment with zoledronic acid. *Intern Med J* 2008. [PMID 18284458](https://pubmed.ncbi.nlm.nih.gov/18284458/)

---

## 13. 데노수맙 (Denosumab)

88. Hu MI, Glezerman IG, Leboulleux S, et al. Denosumab for treatment of hypercalcemia of malignancy. *J Clin Endocrinol Metab* 2014. [PMID 24915117](https://pubmed.ncbi.nlm.nih.gov/24915117/)
    *The single-arm study in bisphosphonate-refractory disease; the model reproduces both the deeper nadir and the longer normocalcaemic interval (scenario 8 vs 5).*
89. Thosani S, Hu MI. Denosumab: a new agent in the management of hypercalcemia of malignancy. *Future Oncol* 2015. [PMID 26403973](https://pubmed.ncbi.nlm.nih.gov/26403973/)
90. Diel IJ, Body JJ, Stopeck AT, et al. The role of denosumab in the prevention of hypercalcaemia of malignancy in cancer patients with metastatic bone disease. *Eur J Cancer* 2015. [PMID 25976743](https://pubmed.ncbi.nlm.nih.gov/25976743/)
91. Sutjandra L, Rodriguez RD, Doshi S, et al. Population pharmacokinetic meta-analysis of denosumab in healthy subjects and postmenopausal women with osteopenia or osteoporosis. *Clin Pharmacokinet* 2011. [PMID 22087866](https://pubmed.ncbi.nlm.nih.gov/22087866/)
    *`DMB_KA`, `DMB_F`, `DMB_VD`, `DMB_KEL` and the ~7-day t_max the model reproduces.*
92. Gibiansky L, Sutjandra L, Doshi S, et al. Population pharmacokinetic analysis of denosumab in patients with bone metastases from solid tumours. *Clin Pharmacokinet* 2012. [PMID 22420579](https://pubmed.ncbi.nlm.nih.gov/22420579/)
    *Target-mediated disposition and the absence of a renal clearance pathway — the structural basis of note H.*
93. Fizazi K, Carducci M, Smith M, et al. Denosumab versus zoledronic acid for treatment of bone metastases in men with castration-resistant prostate cancer: a randomised, double-blind study. *Lancet* 2011. [PMID 21353695](https://pubmed.ncbi.nlm.nih.gov/21353695/)
94. Stopeck AT, Lipton A, Body JJ, et al. Denosumab compared with zoledronic acid for the treatment of bone metastases in patients with advanced breast cancer: a randomized, double-blind study. *J Clin Oncol* 2010. [PMID 21060033](https://pubmed.ncbi.nlm.nih.gov/21060033/)
    *Greater and more sustained CTX suppression with denosumab, which the model produces from the TMDD structure (CTX nadir 0.038 vs 0.153).*
95. Body JJ, Bone HG, de Boer RH, et al. Hypocalcaemia in patients with metastatic bone disease treated with denosumab. *Eur J Cancer* 2015. [PMID 26093811](https://pubmed.ncbi.nlm.nih.gov/26093811/)
    *The predicted harm: with no rate-limiting delivery step, RANKL blockade over-suppresses more readily than a bisphosphonate does.*

---

## 14. 체외 순환 요법 · 기타 약제 (Extracorporeal therapy and other agents)

96. Camus C, Charasse C, Jouannic-Montier I, et al. Calcium free hemodialysis: experience in the treatment of 33 patients with severe hypercalcemia. *Intensive Care Med* 1996. [PMID 8857118](https://pubmed.ncbi.nlm.nih.gov/8857118/)
    *Scenario 17: the only intervention in the model that removes calcium faster than bone can add it.*
97. Warrell RP Jr, Israel R, Frisone M, Snyder T, Gaynor JJ, Bockman RS. Gallium nitrate for acute treatment of cancer-related hypercalcemia. A randomized, double-blind comparison to calcitonin. *Ann Intern Med* 1988. [PMID 3282463](https://pubmed.ncbi.nlm.nih.gov/3282463/)
    *Historical comparator; not implemented.*

---

## 15. 부동 · 침상안정 (Immobilisation)

98. Stewart AF, Adler M, Byers CM, Segre GV, Broadus AE. Calcium homeostasis in immobilization: an example of resorptive hypercalciuria. *N Engl J Med* 1982. [PMID 6280047](https://pubmed.ncbi.nlm.nih.gov/6280047/)
    *The magnitude of the resorption increase and formation decrease behind the model's `IMMOB` coefficients (+18% / −25%).*
99. Bergstrom WH. Hypercalciuria and hypercalcemia complicating immobilization. *Am J Dis Child* 1978. [PMID 350037](https://pubmed.ncbi.nlm.nih.gov/350037/)

---

## 16. 정량적 시스템 약리학 · 골-칼슘 모델링 (QSP and bone-calcium modelling)

100. Lemaire V, Tobin FL, Greller LD, Cho CR, Suva LJ. Modeling the interactions between osteoblast and osteoclast activities in bone remodeling. *J Theor Biol* 2004. [PMID 15234198](https://pubmed.ncbi.nlm.nih.gov/15234198/)
     *The ancestor of this model's bone-cell block: RANKL/OPG-driven osteoclast and osteoblast populations with PTH acting in opposite directions on the two ligands.*
101. Peterson MC, Riggs MM. A physiologically based mathematical model of integrated calcium homeostasis and bone remodeling. *Bone* 2010. [PMID 19732857](https://pubmed.ncbi.nlm.nih.gov/19732857/)
     *The reference multiscale calcium model. The present model differs mainly in taking the RENAL and VOLUME arms seriously enough to produce a bifurcation, which a fixed-volume formulation cannot.*
102. Marathe A, Peterson MC, Mager DE. Integrated cellular bone homeostasis model for denosumab pharmacodynamics in multiple myeloma patients. *J Pharmacol Exp Ther* 2008. [PMID 18460643](https://pubmed.ncbi.nlm.nih.gov/18460643/)
     *Precedent for the denosumab TMDD block and for the CTX readout.*
103. Riggs MM, Peterson MC, Gastonguay MR. Multiscale physiology-based modeling of mineral bone disorder in patients with impaired kidney function. *J Clin Pharmacol* 2012. [PMID 22232752](https://pubmed.ncbi.nlm.nih.gov/22232752/)
104. Peterson MC, Riggs MM. Predicting nonlinear changes in bone mineral density over time using a multiscale systems pharmacology model. *CPT Pharmacometrics Syst Pharmacol* 2012. [PMID 23835796](https://pubmed.ncbi.nlm.nih.gov/23835796/)
105. Raposo JF, Sobrinho LG, Ferreira HG. A minimal mathematical model of calcium homeostasis. *J Clin Endocrinol Metab* 2002. [PMID 12213894](https://pubmed.ncbi.nlm.nih.gov/12213894/)
     *The minimal PTH-calcium feedback loop; the sigmoid CaSR secretion curve (`PTH_SET`, `PTH_N`) follows this form.*

---

## 17. 근거가 없는 부분 (Where the model asserts rather than cites)

정직하게 적어 둡니다. 아래 항목은 문헌으로 직접 뒷받침되지 않은 **모델링 선택**입니다.

The following are modelling choices, not findings. They are stated here so that
nobody mistakes an output of the model for a result from the literature.

- **The CNS adaptation fraction (`ADAPT_FRAC` = 0.45).** That chronic
  hypercalcaemia is better tolerated than acute hypercalcaemia is universally
  reported and nowhere quantified. The value 0.45 was chosen so that a patient
  at an ionised calcium of 1.9 mmol/L reached over months is symptomatic but
  not obtunded, and one who arrives there in three days is. It is a placeholder
  with a mechanism attached, and any conclusion that depends on its exact value
  should be treated as untested.
- **The collapse boundary** (ECF deficit > 38% or total calcium > 5.0 mmol/L)
  is the edge of the region over which the constitutive relations were
  calibrated. It is *not* a mortality model and the day on which a trajectory
  crosses it is not a predicted time of death.
- **`VOM_MAX` = 1.4 L/day** and the symptom-to-vomiting Hill function are
  order-of-magnitude estimates. Note B shows the fold is more sensitive to this
  arm (38%) than to the concentrating defect (11%), so this is the parameter
  most worth measuring and the one the model is least sure of.
- **The saddle-node location (23.07 mmol/day)** is an output of the calibrated
  system, not an observed quantity. What is testable is the *ordering* it
  implies — that the renal action of PTHrP costs more of the safety margin than
  its bone action does at the same plasma calcium, and that fluid restores the
  margin faster than any antiresorptive.
- **Tumour dynamics** are a single logistic compartment. Growth rate,
  carrying capacity and the chemotherapy kill rate in scenario 21 are
  illustrative; no regimen is being simulated.
- **`OC_IND` = 0.06**, the RANKL-independent fraction of osteoclastogenesis,
  is inferred from the observation that denosumab suppresses CTX by roughly
  90% rather than completely, and from the existence of TNF-driven
  osteoclastogenesis (ref. 42). It is a fitted, not a measured, number.
