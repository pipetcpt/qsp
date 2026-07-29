# References — Chemotherapy-Induced Neutropenia / Febrile Neutropenia QSP Model

**Every entry below was resolved through the NCBI E-utilities API** (`esearch` →
`esummary`) rather than written from memory: the search returned the PMID, and
the first author, year, journal and title printed here are the API's own record
for that PMID. Nothing in this list is a recalled citation. Where the model
uses a specific *number* from a paper, that number was read out of the abstract
text fetched from the same API, and the section marked **[ANCHOR]** or
**[PREDICTION TARGET]** says exactly which number and what it is used for.

Sections:

1. [Semi-mechanistic myelosuppression modelling](#1-semi-mechanistic-myelosuppression-modelling)
2. [Granulopoiesis and neutrophil kinetics](#2-granulopoiesis-and-neutrophil-kinetics)
3. [G-CSF pharmacology and target-mediated disposition](#3-g-csf-pharmacology-and-target-mediated-disposition)
4. [G-CSF clinical trials — the calibration anchors](#4-g-csf-clinical-trials--the-calibration-anchors)
5. [Trilaciclib and CDK4/6-mediated myelopreservation](#5-trilaciclib-and-cdk46-mediated-myelopreservation)
6. [Cytotoxic pharmacokinetics and pharmacogenomics](#6-cytotoxic-pharmacokinetics-and-pharmacogenomics)
7. [Mucosal barrier injury, the microbiome, and translocation](#7-mucosal-barrier-injury-the-microbiome-and-translocation)
8. [Febrile neutropenia — epidemiology, risk scores, outcomes](#8-febrile-neutropenia--epidemiology-risk-scores-outcomes)
9. [Antimicrobial prophylaxis and empirical therapy](#9-antimicrobial-prophylaxis-and-empirical-therapy)
10. [Guidelines and the 20% decision threshold](#10-guidelines-and-the-20-decision-threshold)
11. [Dose intensity and tumour control](#11-dose-intensity-and-tumour-control)
12. [Platelet and erythroid lineages](#12-platelet-and-erythroid-lineages)
13. [Count versus function — corticosteroid and the ANC that lies](#13-count-versus-function--corticosteroid-and-the-anc-that-lies)
14. [Parameter provenance table](#14-parameter-provenance-table)
15. [Literature the model does not reproduce](#15-literature-the-model-does-not-reproduce)

---

## 1. Semi-mechanistic myelosuppression modelling

The transit-compartment topology (`Prol → Tr1 → Tr2 → Tr3 → Circ`, feedback
`(Circ0/Circ)^γ`, drug effect `slope·C` on the proliferation term) is taken
directly from the first paper below and is not re-derived here. The system
parameters MTT = 89.3 h and γ = 0.170 are its published estimates and are held
**fixed** in this model.

| # | Reference | Used for |
|---|-----------|----------|
| 1 | Friberg LE, et al. *Model of chemotherapy-induced myelosuppression with parameter consistency across drugs.* J Clin Oncol 2002. [PMID 12488418](https://pubmed.ncbi.nlm.nih.gov/12488418/) | The entire neutrophil topology; MTT and γ; the fact that the drug slopes are large enough that (1−E_drug) goes negative |
| 2 | Quartino AL, et al. *A simultaneous analysis of the time-course of leukocytes and neutrophils following docetaxel administration using a semi-mechanistic myelosuppression model.* Invest New Drugs 2012. [PMID 21153753](https://pubmed.ncbi.nlm.nih.gov/21153753/) | Docetaxel-specific application; leukocyte vs neutrophil distinction |
| 3 | Quartino AL, et al. *Characterization of endogenous G-CSF and the inverse correlation to chemotherapy-induced neutropenia in patients with breast cancer using population modeling.* Pharm Res 2014. [PMID 24919931](https://pubmed.ncbi.nlm.nih.gov/24919931/) | Endogenous G-CSF rises as the ANC falls — the `ENDOG` node and the negative-feedback arrow in cluster 7 |
| 4 | Soto E, et al. *Predictive ability of a semi-mechanistic model for neutropenia in the development of novel anti-cancer agents: two case studies.* Invest New Drugs 2011. [PMID 20449627](https://pubmed.ncbi.nlm.nih.gov/20449627/) | Evidence that the topology transfers to new agents, which is what licenses the per-agent slope approach |
| 5 | Wallin JE, et al. *Model-based neutrophil-guided dose adaptation in chemotherapy: evaluation of predicted outcome with different types and amounts of information.* Basic Clin Pharmacol Toxicol 2010. [PMID 20050841](https://pubmed.ncbi.nlm.nih.gov/20050841/) | Precedent for the closed-loop clinician controller in analysis A12 |
| 6 | Sandström M, et al. *Model describing the relationship between pharmacokinetics and hematologic toxicity of the epirubicin-docetaxel regimen in breast cancer patients.* J Clin Oncol 2005. [PMID 15585753](https://pubmed.ncbi.nlm.nih.gov/15585753/) | Additive drug effects in combination regimens |
| 7 | Sandström M, et al. *Population analysis of the pharmacokinetics and the haematological toxicity of the fluorouracil-epirubicin-cyclophosphamide regimen in breast cancer patients.* Cancer Chemother Pharmacol 2006. [PMID 16465545](https://pubmed.ncbi.nlm.nih.gov/16465545/) | Three-agent combination toxicity |
| 8 | Craig M, et al. *A Mathematical Model of Granulopoiesis Incorporating the Negative Feedback Dynamics and Kinetics of G-CSF/Neutrophil Binding and Internalization.* Bull Math Biol 2016. [PMID 27324993](https://pubmed.ncbi.nlm.nih.gov/27324993/) | The G-CSF **receptor-mediated internalisation** term — the model's `elim_tm` |
| 9 | Foley C, Mackey MC. *Dynamic hematological disease: a review.* J Math Biol 2009. [PMID 18317766](https://pubmed.ncbi.nlm.nih.gov/18317766/) | Delay-driven oscillation and why a long transit chain behaves as a pure delay |
| 10 | Shochat E, et al. *Novel strategies for granulocyte colony-stimulating factor treatment of severe prolonged neutropenia suggested by mathematical modeling.* Clin Cancer Res 2008. [PMID 18927273](https://pubmed.ncbi.nlm.nih.gov/18927273/) | Precedent for using a model to choose G-CSF timing rather than dose |
| 11 | Krzyzanski W. *A cell-level model of pharmacodynamics-mediated drug disposition.* J Pharmacokinet Pharmacodyn 2016. [PMID 27612462](https://pubmed.ncbi.nlm.nih.gov/27612462/) | The formal statement of "the drug is cleared by the cells it creates" |
| 12 | Krzyzanski W, et al. *Pharmacodynamic model for chemoradiotherapy-induced thrombocytopenia in mice.* J Pharmacokinet Pharmacodyn 2015. [PMID 26341875](https://pubmed.ncbi.nlm.nih.gov/26341875/) | Transit-model transfer to the platelet lineage |
| 13 | Schmitt A, et al. *Factors for hematopoietic toxicity of carboplatin: refining the targeting of carboplatin systemic exposure.* J Clin Oncol 2010. [PMID 20855828](https://pubmed.ncbi.nlm.nih.gov/20855828/) | Carboplatin exposure → platelet nadir; why the platelet slope for carboplatin is the largest in the agent table |

---

## 2. Granulopoiesis and neutrophil kinetics

These set the compartment sizes and the three-pool structure (marrow storage,
circulating, marginated) that make the "the number you measure is not the pool
that matters" results possible.

| # | Reference | Used for |
|---|-----------|----------|
| 14 | Cartwright GE, et al. *The kinetics of granulopoiesis in normal man.* Blood 1964. [PMID 14235362](https://pubmed.ncbi.nlm.nih.gov/14235362/) | Marrow transit time and the size of the marrow storage pool relative to blood |
| 15 | Athens JW, et al. *Leukokinetic studies. IV. The total blood, circulating and marginal granulocyte pools and the granulocyte turnover rate in normal subjects.* J Clin Invest 1961. [PMID 13684958](https://pubmed.ncbi.nlm.nih.gov/13684958/) | `marg_ratio = 1.0` — the marginated pool is comparable in size to the circulating pool |
| 16 | Dancey JT, et al. *Neutrophil kinetics in man.* J Clin Invest 1976. [PMID 956397](https://pubmed.ncbi.nlm.nih.gov/956397/) | Circulating neutrophil half-life (`tau_circ_N = 7.9 h`) and marrow pool sizes |
| 17 | Lord BI, et al. *The kinetics of human granulopoiesis following treatment with granulocyte colony-stimulating factor in vivo.* Proc Natl Acad Sci U S A 1989. [PMID 2480603](https://pubmed.ncbi.nlm.nih.gov/2480603/) | G-CSF **shortens marrow transit time** — the mechanistic basis of `Emax_ktr` |
| 18 | Lord BI, et al. *Haemopoietic progenitor and myeloid cell kinetics in humans treated with interleukin-3 and granulocyte/macrophage colony-stimulating factor in combination.* Int J Cancer 1994. [PMID 7960217](https://pubmed.ncbi.nlm.nih.gov/7960217/) | Growth-factor effects on progenitor amplification |
| 19 | Summers C, et al. *Neutrophil kinetics in health and disease.* Trends Immunol 2010. [PMID 20620114](https://pubmed.ncbi.nlm.nih.gov/20620114/) | Review of the pool structure; the disputed circulating lifespan |
| 20 | Hidalgo A, et al. *The Neutrophil Life Cycle.* Trends Immunol 2019. [PMID 31153737](https://pubmed.ncbi.nlm.nih.gov/31153737/) | Maturation stages drawn in cluster 4 |
| 21 | Hidalgo A, et al. *Dimensions of neutrophil life and fate.* Semin Immunol 2021. [PMID 34711490](https://pubmed.ncbi.nlm.nih.gov/34711490/) | Heterogeneity of released cells — the `LEFTSHIFT` node |
| 22 | von Vietinghoff S, Ley K. *Homeostatic regulation of blood neutrophil counts.* J Immunol 2008. [PMID 18832668](https://pubmed.ncbi.nlm.nih.gov/18832668/) | The feedback loop that γ represents |
| 23 | Eash KJ, et al. *CXCR2 and CXCR4 antagonistically regulate neutrophil trafficking from murine bone marrow.* J Clin Invest 2010. [PMID 20516641](https://pubmed.ncbi.nlm.nih.gov/20516641/) | The retention/release switch behind `krel` |
| 24 | Semerad CL, et al. *G-CSF potently inhibits osteoblast activity and CXCL12 mRNA expression in the bone marrow.* Blood 2005. [PMID 16037394](https://pubmed.ncbi.nlm.nih.gov/16037394/) | Mechanism of the G-CSF-driven storage release |
| 25 | Layton JE, et al. *Evidence for a novel in vivo control mechanism of granulopoiesis: mature cell-related control of a regulatory growth factor.* Blood 1989. [PMID 2475185](https://pubmed.ncbi.nlm.nih.gov/2475185/) | The original demonstration that mature neutrophils control G-CSF levels — the direct empirical basis for the target-mediated clearance term |
| 26 | Montaldo E, et al. *Cellular and transcriptional dynamics of human neutrophils at steady state and upon stress.* Nat Immunol 2022. [PMID 36138183](https://pubmed.ncbi.nlm.nih.gov/36138183/) | Modern staging of the maturation chain |
| 27 | Ng LG, et al. *From complexity to consensus: A roadmap for neutrophil classification.* Immunity 2025. [PMID 40763729](https://pubmed.ncbi.nlm.nih.gov/40763729/) | Current nomenclature for the compartments |
| 28 | Bendall LJ, Bradstock KF. *G-CSF: From granulopoietic stimulant to bone marrow stem cell mobilizing agent.* Cytokine Growth Factor Rev 2014. [PMID 25131807](https://pubmed.ncbi.nlm.nih.gov/25131807/) | The full range of G-CSF actions, which is why the model splits it into five |

---

## 3. G-CSF pharmacology and target-mediated disposition

The single most consequential piece of pharmacology in the model:
pegfilgrastim's renal route is abolished by PEGylation, so the only clearance
route left is internalisation through G-CSF receptors on neutrophils and their
precursors — the drug is eliminated by the cells it creates.

| # | Reference | Used for |
|---|-----------|----------|
| 29 | Molineux G. *Pegfilgrastim: using pegylation technology to improve neutropenia support in cancer patients.* Anticancer Drugs 2003. [PMID 12679729](https://pubmed.ncbi.nlm.nih.gov/12679729/) | PEGylation removes renal clearance; the rationale for a single fixed dose |
| 30 | Roskos LK, et al. *Pharmacokinetic/pharmacodynamic modeling of pegfilgrastim in healthy subjects.* J Clin Pharmacol 2006. [PMID 16809800](https://pubmed.ncbi.nlm.nih.gov/16809800/) | The reference PK/PD structure for `Vmax_nm`, `Km_nm`, and the flip-flop absorption |
| 31 | Yang BB, Kido A. *Pharmacokinetics and pharmacodynamics of pegfilgrastim.* Clin Pharmacokinet 2011. [PMID 21456630](https://pubmed.ncbi.nlm.nih.gov/21456630/) | Concentration ranges and the ANC-dependent half-life |
| 32 | Arvedson T, et al. *Design Rationale and Development Approach for Pegfilgrastim as a Long-Acting Granulocyte Colony-Stimulating Factor.* BioDrugs 2015. [PMID 25998211](https://pubmed.ncbi.nlm.nih.gov/25998211/) | Why fixed weight-independent dosing is self-titrating |
| 33 | Wiczling P, et al. *Population pharmacokinetic modelling of filgrastim in healthy adults following intravenous and subcutaneous administrations.* Clin Pharmacokinet 2009. [PMID 19902989](https://pubmed.ncbi.nlm.nih.gov/19902989/) | Filgrastim `V`, `ka`, and the linear (renal) clearance that gives its ~3.5 h half-life |
| 34 | Krzyzanski W, et al. *Population modeling of filgrastim PK-PD in healthy adults following intravenous and subcutaneous administrations.* J Clin Pharmacol 2010. [PMID 20881223](https://pubmed.ncbi.nlm.nih.gov/20881223/) | The joint PK-PD structure the effect compartment stands in for |
| 35 | Dale DC, et al. *A systematic literature review of the efficacy, effectiveness, and safety of filgrastim.* Support Care Cancer 2018. [PMID 28939926](https://pubmed.ncbi.nlm.nih.gov/28939926/) | Effect sizes across settings; adverse-effect nodes in cluster 8 |

---

## 4. G-CSF clinical trials — the calibration anchors

**All four numbers the model is fitted to come from two of these trials.** They
are quoted here verbatim from the abstracts retrieved via `efetch`.

| # | Reference | Number used |
|---|-----------|-------------|
| 36 | Crawford J, et al. *Reduction by granulocyte colony-stimulating factor of fever and neutropenia induced by chemotherapy in patients with small-cell lung cancer.* N Engl J Med 1991. [PMID 1711156](https://pubmed.ncbi.nlm.nih.gov/1711156/) | **[ANCHOR ×2]** "the median duration of grade IV neutropenia (absolute neutrophil count, less than 0.5 × 10⁹ per liter) was **six days with placebo as compared with one day with G-CSF**" on cyclophosphamide/doxorubicin/etoposide, G-CSF days 4–17. → fits `sens_global` and `GPD["Emax_amp"]`. Also: FN in 77% vs 40%, IV antibiotic days and hospital days each roughly halved, bone pain in 20% |
| 37 | Green MD, et al. *A randomized double-blind multicenter phase III study of fixed-dose single-administration pegfilgrastim versus daily filgrastim in patients receiving myelosuppressive chemotherapy.* Ann Oncol 2003. [PMID 12488289](https://pubmed.ncbi.nlm.nih.gov/12488289/) | **[ANCHOR]** "The mean duration of grade 4 neutropenia in cycle 1 was **1.8 and 1.6 days for the pegfilgrastim and filgrastim groups**" after doxorubicin 60 + docetaxel 75 mg/m². → the 1.8 d figure fits `sens_taxane`; **[PREDICTION TARGET]** the 1.6 d filgrastim arm is not fitted. Also FN 13% (peg) vs 20% (filgrastim) across all cycles |
| 38 | Vogel CL, et al. *First and subsequent cycle use of pegfilgrastim prevents febrile neutropenia in patients with breast cancer: a multicenter, double-blind, placebo-controlled phase III study.* J Clin Oncol 2005. [PMID 15718314](https://pubmed.ncbi.nlm.nih.gov/15718314/) | **[ANCHOR ×2]** docetaxel 100 mg/m² q3w, pegfilgrastim 6 mg on day 2, n = 928: FN **1% vs 17%**. → fits `hFN_p` and `hFN_k`. **[PREDICTION TARGET]** FN-related hospitalisation 1% vs 14%; IV anti-infectives 2% vs 10%; on-time dosing 80% vs 78% |
| 39 | Holmes FA, et al. *Blinded, randomized, multicenter study to evaluate single administration pegfilgrastim once per cycle versus daily filgrastim as an adjunct to chemotherapy in patients with high-risk stage II or stage III/IV breast cancer.* J Clin Oncol 2002. [PMID 11821454](https://pubmed.ncbi.nlm.nih.gov/11821454/) | Independent replication of the equivalence of one pegfilgrastim dose and daily filgrastim ("the difference in the mean duration of severe neutropenia … was less than 1 day") |
| 40 | Kuderer NM, et al. *Impact of primary prophylaxis with granulocyte colony-stimulating factor on febrile neutropenia and mortality in adult cancer patients receiving chemotherapy: a systematic review.* J Clin Oncol 2007. [PMID 17634496](https://pubmed.ncbi.nlm.nih.gov/17634496/) | Pooled FN and mortality effect sizes used as a plausibility check on A5 |
| 41 | Lyman GH, et al. *The impact of the granulocyte colony-stimulating factor on chemotherapy dose intensity and cancer survival: a systematic review and meta-analysis of randomized controlled trials.* Ann Oncol 2013. [PMID 23788754](https://pubmed.ncbi.nlm.nih.gov/23788754/) | The RDI ↔ survival link used — with explicit caveats — in the tumour layer |
| 42 | Lyman GH, et al. *Long-term outcomes of myeloid growth factor treatment.* J Natl Compr Canc Netw 2011. [PMID 21900223](https://pubmed.ncbi.nlm.nih.gov/21900223/) | The secondary MDS/AML question drawn in cluster 8 |
| 43 | Burris HA, et al. *Pegfilgrastim on the Same Day Versus Next Day of Chemotherapy in Patients With Breast Cancer, Non-Small-Cell Lung Cancer, Ovarian Cancer, and Non-Hodgkin's Lymphoma: Results of Four Multicenter, Double-Blind, Randomized Phase II Studies.* J Oncol Pract 2010. [PMID 20808556](https://pubmed.ncbi.nlm.nih.gov/20808556/) | The same-day-dosing question that analysis A6 answers mechanistically |
| 44 | Crawford J, et al. *Use of prophylactic pegfilgrastim for chemotherapy-induced neutropenia in the US: A review of adherence to present guidelines for usage.* Cancer Treat Res Commun 2021. [PMID 34655862](https://pubmed.ncbi.nlm.nih.gov/34655862/) | Real-world timing and adherence, i.e. how often the day-2 rule is actually followed |
| 45 | Meropol NJ, et al. *Randomized, placebo-controlled, multicenter trial of granulocyte-macrophage colony-stimulating factor as infection prophylaxis in oncologic surgery.* J Clin Oncol 1998. [PMID 9508204](https://pubmed.ncbi.nlm.nih.gov/9508204/) | Contrast case: raising the ANC without a neutropenic deficit to correct |
| 46 | Yau JC, et al. *Randomized placebo-controlled trial of granulocyte-macrophage colony-stimulating-factor support for dose-intensive cyclophosphamide, etoposide, and cisplatin.* Am J Hematol 1996. [PMID 8602629](https://pubmed.ncbi.nlm.nih.gov/8602629/) | Growth-factor support of a dose-intensified regimen |

---

## 5. Trilaciclib and CDK4/6-mediated myelopreservation

| # | Reference | Used for |
|---|-----------|----------|
| 47 | He S, et al. *Transient CDK4/6 inhibition protects hematopoietic stem cells from chemotherapy-induced exhaustion.* Sci Transl Med 2017. [PMID 28446688](https://pubmed.ncbi.nlm.nih.gov/28446688/) | The whole mechanism: G1 arrest removes HSPC from the drug-sensitive cycling fraction. Implemented as `Fcyc = 1 − Imax·C/(IC50+C)` acting only on the cycle-dependent share `f_cyc` of each agent's kill |
| 48 | Weiss JM, et al. *Myelopreservation with the CDK4/6 inhibitor trilaciclib in patients with small-cell lung cancer receiving first-line chemotherapy: a phase Ib/randomized phase II trial.* Ann Oncol 2019. [PMID 31504118](https://pubmed.ncbi.nlm.nih.gov/31504118/) | Direction of effect: "Improvements were seen with trilaciclib in neutrophil, RBC and lymphocyte measures"; ≥G3 AEs 50% vs 83.8%; ORR 66.7% vs 56.8%, PFS 6.2 vs 5.0 mo, OS 10.9 vs 10.6 mo — i.e. **no antitumour penalty**, which is what licenses the model treating trilaciclib as marrow-selective |
| 49 | Weiss J, et al. *Effects of Trilaciclib on Chemotherapy-Induced Myelosuppression and Patient-Reported Outcomes in Patients with Extensive-Stage Small Cell Lung Cancer: Pooled Results from Three Phase II Randomized, Double-Blind, Placebo-Controlled Studies.* Clin Lung Cancer 2021. [PMID 33895103](https://pubmed.ncbi.nlm.nih.gov/33895103/) | Primary endpoints were "duration of severe neutropenia … in cycle 1 and occurrence of severe neutropenia"; multilineage reduction confirmed. **The abstract does not state the numerical DSN values**, so the model's trilaciclib output (A10) is reported as a prediction whose *direction* is validated and whose *magnitude* is not checked against these papers |

---

## 6. Cytotoxic pharmacokinetics and pharmacogenomics

| # | Reference | Used for |
|---|-----------|----------|
| 50 | Bruno R, et al. *Population pharmacokinetics/pharmacodynamics of docetaxel in phase II studies in patients with cancer.* J Clin Oncol 1998. [PMID 9440742](https://pubmed.ncbi.nlm.nih.gov/9440742/) | Docetaxel 3-compartment PK: CL ≈ 21.6 L/h/m², V1 ≈ 7.4 L/m², and the α₁-acid-glycoprotein/albumin covariate on unbound fraction |
| 51 | Baille P, et al. *Optimal sampling strategies for bayesian estimation of docetaxel (Taxotere) clearance.* Clin Cancer Res 1997. [PMID 9815840](https://pubmed.ncbi.nlm.nih.gov/9815840/) | Confirmatory docetaxel clearance estimates |
| 52 | Henningsson A, et al. *Mechanism-based pharmacokinetic model for paclitaxel.* J Clin Oncol 2001. [PMID 11600609](https://pubmed.ncbi.nlm.nih.gov/11600609/) | Paclitaxel's saturable distribution and elimination; the reason its PD is modelled as a **threshold** rather than an AUC effect |
| 53 | Henningsson A, et al. *Population pharmacokinetic modelling of unbound and total plasma concentrations of paclitaxel in cancer patients.* Eur J Cancer 2003. [PMID 12736110](https://pubmed.ncbi.nlm.nih.gov/12736110/) | Unbound paclitaxel, the species the threshold applies to |
| 54 | Joerger M, et al. *Population pharmacokinetics and pharmacodynamics of paclitaxel and carboplatin in ovarian cancer patients …* Clin Cancer Res 2007. [PMID 17975154](https://pubmed.ncbi.nlm.nih.gov/17975154/) | The paclitaxel time-above-threshold → neutropenia relationship, and joint carboplatin PD |
| 55 | Calvert AH, et al. *Carboplatin dosage: prospective evaluation of a simple formula based on renal function.* J Clin Oncol 1989. [PMID 2681557](https://pubmed.ncbi.nlm.nih.gov/2681557/) | `CL_carboplatin (mL/min) = GFR + 25`, implemented literally, which is what makes renal function a covariate on the platelet nadir |
| 56 | Calvert AH, et al. *Carboplatin dosage: prospective evaluation of a simple formula based on renal function.* J Clin Oncol 2023 (classic-paper reissue). [PMID 37757592](https://pubmed.ncbi.nlm.nih.gov/37757592/) | Same formula, reissued |
| 57 | Amstutz U, et al. *Clinical Pharmacogenetics Implementation Consortium (CPIC) Guideline for Dihydropyrimidine Dehydrogenase Genotype and Fluoropyrimidine Dosing: 2017 Update.* Clin Pharmacol Ther 2018. [PMID 29152729](https://pubmed.ncbi.nlm.nih.gov/29152729/) | The `DPYD` node — a genotype that changes exposure, not sensitivity |
| 58 | Karas S, Innocenti F. *All You Need to Know About UGT1A1 Genetic Testing for Patients Treated With Irinotecan: A Practitioner-Friendly Guide.* JCO Oncol Pract 2022. [PMID 34860573](https://pubmed.ncbi.nlm.nih.gov/34860573/) | The `UGT1A1` node and SN-38 exposure |
| 59 | Chen S, et al. *A novel UGT1 marker associated with better tolerance against irinotecan-induced severe neutropenia in metastatic colorectal cancer patients.* Pharmacogenomics J 2015. [PMID 25778466](https://pubmed.ncbi.nlm.nih.gov/25778466/) | Genotype → severe neutropenia, quantitatively |
| 60 | Björn N, et al. *Genes and variants in hematopoiesis-related pathways are associated with gemcitabine/carboplatin-induced thrombocytopenia.* Pharmacogenomics J 2020. [PMID 31616045](https://pubmed.ncbi.nlm.nih.gov/31616045/) | Host determinants of the platelet nadir specifically |

---

## 7. Mucosal barrier injury, the microbiome, and translocation

The model's second pathway: the same cytotoxic that removes the defence also
opens the door. Barrier integrity is a state, translocation is proportional to
the breach, and the FN hazard is multiplied by it.

| # | Reference | Used for |
|---|-----------|----------|
| 61 | Sonis ST. *The pathobiology of mucositis.* Nat Rev Cancer 2004. [PMID 15057287](https://pubmed.ncbi.nlm.nih.gov/15057287/) | The five-phase model behind the `ROS → NF-κB → TNF → MMP → barrier` chain in cluster 11 |
| 62 | Spielberger R, et al. *Palifermin for oral mucositis after intensive therapy for hematologic cancers.* N Engl J Med 2004. [PMID 15602019](https://pubmed.ncbi.nlm.nih.gov/15602019/) | The `PALIF` intervention node — a barrier-directed rather than marrow-directed therapy |
| 63 | Taur Y, et al. *Intestinal domination and the risk of bacteremia in patients undergoing allogeneic hematopoietic stem cell transplantation.* Clin Infect Dis 2012. [PMID 22718773](https://pubmed.ncbi.nlm.nih.gov/22718773/) | Microbiome domination precedes bloodstream infection — the `MICROB → TRANSL` arrow |
| 64 | Ubeda C, et al. *Vancomycin-resistant Enterococcus domination of intestinal microbiota is enabled by antibiotic treatment in mice and precedes bloodstream invasion in humans.* J Clin Invest 2010. [PMID 21099116](https://pubmed.ncbi.nlm.nih.gov/21099116/) | Antibiotic-driven domination — the cost side of prophylaxis |
| 65 | Montassier E, et al. *Chemotherapy-driven dysbiosis in the intestinal microbiome.* Aliment Pharmacol Ther 2015. [PMID 26147207](https://pubmed.ncbi.nlm.nih.gov/26147207/) | Chemotherapy itself as the dysbiosis driver |

---

## 8. Febrile neutropenia — epidemiology, risk scores, outcomes

| # | Reference | Used for |
|---|-----------|----------|
| 66 | Bodey GP, et al. *Quantitative relationships between circulating leukocytes and infection in patients with acute leukemia.* Ann Intern Med 1966. [PMID 5216294](https://pubmed.ncbi.nlm.nih.gov/5216294/) | The founding observation that infection risk depends on **both depth and duration** of neutropenia — the reason the model's hazard is `∫ depth^p dt` and not a function of the nadir |
| 67 | Sickles EA, et al. *Clinical presentation of infection in granulocytopenic patients.* Arch Intern Med 1975. [PMID 1052668](https://pubmed.ncbi.nlm.nih.gov/1052668/) | Localising signs are absent when there are no neutrophils — the `BLUNT` node |
| 68 | Kuderer NM, et al. *Mortality, morbidity, and cost associated with febrile neutropenia in adult cancer patients.* Cancer 2006. [PMID 16575919](https://pubmed.ncbi.nlm.nih.gov/16575919/) | FN case-fatality and cost, used for the endpoint cluster |
| 69 | Lyman GH, et al. *Predicting individual risk of neutropenic complications in patients receiving cancer chemotherapy.* Cancer 2011. [PMID 21509769](https://pubmed.ncbi.nlm.nih.gov/21509769/) | The patient-level covariates in the `PATRISK` node |
| 70 | Crawford J, et al. *Chemotherapy-induced neutropenia: risks, consequences, and new directions for its management.* Cancer 2004. [PMID 14716755](https://pubmed.ncbi.nlm.nih.gov/14716755/) | Overall framing; the dose-delay/reduction consequence chain |
| 71 | Lyman GH, et al. *Risk of febrile neutropenia among patients with intermediate-grade non-Hodgkin's lymphoma receiving CHOP chemotherapy.* Leuk Lymphoma 2003. [PMID 14959849](https://pubmed.ncbi.nlm.nih.gov/14959849/) | CHOP-21 FN risk, used as a plausibility check on the CHOP21 regimen |
| 72 | Lyman GH, et al. *Risk and timing of hospitalization for febrile neutropenia in patients receiving CHOP, CHOP-R, or CNOP chemotherapy for intermediate-grade non-Hodgkin lymphoma.* Cancer 2003. [PMID 14635075](https://pubmed.ncbi.nlm.nih.gov/14635075/) | **Timing** of FN within the cycle — the observable the delay-line structure predicts |
| 73 | Klastersky J, et al. *Bacteraemia in febrile neutropenic cancer patients.* Int J Antimicrob Agents 2007. [PMID 17689933](https://pubmed.ncbi.nlm.nih.gov/17689933/) | Pathogen mix in cluster 12 |
| 74 | Wisplinghoff H, et al. *Current trends in the epidemiology of nosocomial bloodstream infections in patients with hematological malignancies and solid neoplasms in hospitals in the United States.* Clin Infect Dis 2003. [PMID 12715303](https://pubmed.ncbi.nlm.nih.gov/12715303/) | Gram-positive/gram-negative balance and the catheter route |
| 75 | Gustinetti G, Mikulska M. *Bloodstream infections in neutropenic cancer patients: A practical update.* Virulence 2016. [PMID 27002635](https://pubmed.ncbi.nlm.nih.gov/27002635/) | Current epidemiology and resistance |
| 76 | Rolston KV. *Infections in Cancer Patients with Solid Tumors: A Review.* Infect Dis Ther 2017. [PMID 28160269](https://pubmed.ncbi.nlm.nih.gov/28160269/) | Solid-tumour-specific infection profile |
| 77 | Culakova E, et al. *Patterns of chemotherapy-associated toxicity and supportive care in US oncology practice: a nationwide prospective cohort study.* Cancer Med 2014. [PMID 24706592](https://pubmed.ncbi.nlm.nih.gov/24706592/) | Real-world toxicity and support rates |
| 78 | Culakova E, et al. *The impact of chemotherapy dose intensity and supportive care on the risk of febrile neutropenia in patients with early stage breast cancer: a prospective cohort study.* Springerplus 2015. [PMID 26251780](https://pubmed.ncbi.nlm.nih.gov/26251780/) | The joint RDI–FN relationship that analysis A11 computes |
| 79 | Truong J, et al. *Interpreting febrile neutropenia rates from randomized, controlled trials for consideration of primary prophylaxis in the real world: a systematic review and meta-analysis.* Ann Oncol 2016. [PMID 26712901](https://pubmed.ncbi.nlm.nih.gov/26712901/) | Trial FN rates systematically under-state real-world rates — an important caveat on the Vogel anchor |
| 80 | Klastersky J, et al. *The Multinational Association for Supportive Care in Cancer (MASCC) risk index score: 10 years of use for identifying low-risk febrile neutropenic cancer patients.* Support Care Cancer 2013. [PMID 23443617](https://pubmed.ncbi.nlm.nih.gov/23443617/) | The `MASCC` node |
| 81 | Carmona-Bayonas A, et al. *Prediction of serious complications in patients with seemingly stable febrile neutropenia: validation of the Clinical Index of Stable Febrile Neutropenia in a prospective cohort of patients from the FINITE study.* J Clin Oncol 2015. [PMID 25559804](https://pubmed.ncbi.nlm.nih.gov/25559804/) | CISNE |
| 82 | Carmona-Bayonas A, et al. *Performance of the clinical index of stable febrile neutropenia (CISNE) in different types of infections and tumors.* Clin Transl Oncol 2017. [PMID 27525978](https://pubmed.ncbi.nlm.nih.gov/27525978/) | CISNE generalisability |
| 83 | Paesmans M, et al. *Predicting febrile neutropenic patients at low risk using the MASCC score: does bacteremia matter?* Support Care Cancer 2011. [PMID 20596732](https://pubmed.ncbi.nlm.nih.gov/20596732/) | Score behaviour when bacteraemia is present |

---

## 9. Antimicrobial prophylaxis and empirical therapy

| # | Reference | Used for |
|---|-----------|----------|
| 84 | Bucaneve G, et al. *Levofloxacin to prevent bacterial infection in patients with cancer and neutropenia.* N Engl J Med 2005. [PMID 16148283](https://pubmed.ncbi.nlm.nih.gov/16148283/) | The levofloxacin arm: dose, schedule, and effect size for `abx_Imax` |
| 85 | Cullen MH, et al. *Rational selection of patients for antibacterial prophylaxis after chemotherapy.* J Clin Oncol 2007. [PMID 17947731](https://pubmed.ncbi.nlm.nih.gov/17947731/) | Who benefits — the selection question A13 asks quantitatively |
| 86 | Gafter-Gvili A, et al. *Antibiotic prophylaxis for bacterial infections in afebrile neutropenic patients following chemotherapy.* Cochrane Database Syst Rev 2012. [PMID 22258955](https://pubmed.ncbi.nlm.nih.gov/22258955/) | Pooled prophylaxis effect |
| 87 | Gafter-Gvili A, et al. *Effect of quinolone prophylaxis in afebrile neutropenic patients on microbial resistance: systematic review and meta-analysis.* J Antimicrob Chemother 2007. [PMID 17077101](https://pubmed.ncbi.nlm.nih.gov/17077101/) | The resistance cost — the `GNEGRES` node |
| 88 | Leibovici L, et al. *Antibiotic prophylaxis in neutropenic patients: new evidence, practical decisions.* Cancer 2006. [PMID 16977651](https://pubmed.ncbi.nlm.nih.gov/16977651/) | The decision framing |
| 89 | Averbuch D, et al. *European guidelines for empirical antibacterial therapy for febrile neutropenic patients in the era of growing resistance: summary of the 2011 4th European Conference on Infections in Leukemia.* Haematologica 2013. [PMID 24323983](https://pubmed.ncbi.nlm.nih.gov/24323983/) | Empirical therapy structure in cluster 15 |
| 90 | Link H, et al. *Antimicrobial therapy of unexplained fever in neutropenic patients — guidelines of the Infectious Diseases Working Party (AGIHO) …* Ann Hematol 2003. [PMID 13680173](https://pubmed.ncbi.nlm.nih.gov/13680173/) | Alternative guideline structure |

---

## 10. Guidelines and the 20% decision threshold

The model is used in A14 to ask where the 20% threshold comes from. These are
the documents that state it, and the two health-economic papers that show it is
a cost-effectiveness boundary rather than a biological one.

| # | Reference | Used for |
|---|-----------|----------|
| 91 | Smith TJ, et al. *Recommendations for the Use of WBC Growth Factors: American Society of Clinical Oncology Clinical Practice Guideline Update.* J Clin Oncol 2015. [PMID 26169616](https://pubmed.ncbi.nlm.nih.gov/26169616/) | The 20% primary-prophylaxis threshold and the secondary-prophylaxis rule implemented in the A12 controller |
| 92 | Aapro MS, et al. *2010 update of EORTC guidelines for the use of granulocyte-colony stimulating factor …* Eur J Cancer 2011. [PMID 21095116](https://pubmed.ncbi.nlm.nih.gov/21095116/) | European statement of the same threshold and the patient-risk modifiers |
| 93 | Aapro MS, et al. *EORTC guidelines for the use of granulocyte-colony stimulating factor to reduce the incidence of chemotherapy-induced febrile neutropenia in adult patients with lymphomas and solid tumours.* Eur J Cancer 2006. [PMID 16750358](https://pubmed.ncbi.nlm.nih.gov/16750358/) | The earlier threshold formulation |
| 94 | Crawford J, et al. *Myeloid growth factors.* J Natl Compr Canc Netw 2013. [PMID 24142827](https://pubmed.ncbi.nlm.nih.gov/24142827/) | NCCN risk stratification by regimen |
| 95 | Klastersky J, et al. *Management of febrile neutropaenia: ESMO Clinical Practice Guidelines.* Ann Oncol 2016. [PMID 27664247](https://pubmed.ncbi.nlm.nih.gov/27664247/) | FN management pathway |
| 96 | Freifeld AG, et al. *Clinical practice guideline for the use of antimicrobial agents in neutropenic patients with cancer: 2010 update by the Infectious Diseases Society of America.* Clin Infect Dis 2011. [PMID 21258094](https://pubmed.ncbi.nlm.nih.gov/21258094/) | FN definition (ANC < 0.5 with T ≥ 38.3 °C once or ≥ 38.0 °C for 1 h) used in the model |
| 97 | Taplitz RA, et al. *Outpatient Management of Fever and Neutropenia in Adults Treated for Malignancy: American Society of Clinical Oncology and Infectious Diseases Society of America Clinical Practice Guideline Update.* J Clin Oncol 2018. [PMID 29461916](https://pubmed.ncbi.nlm.nih.gov/29461916/) | Risk-stratified management |
| 98 | Flowers CR, et al. *Antimicrobial prophylaxis and outpatient management of fever and neutropenia in adults treated for malignancy: American Society of Clinical Oncology clinical practice guideline.* J Clin Oncol 2013. [PMID 23319691](https://pubmed.ncbi.nlm.nih.gov/23319691/) | Prophylaxis indications |
| 99 | Calhoun EA, et al. *Granulocyte colony-stimulating factor for chemotherapy-induced neutropenia in patients with small cell lung cancer: the 40% rule revisited.* Pharmacoeconomics 2005. [PMID 16097839](https://pubmed.ncbi.nlm.nih.gov/16097839/) | Direct evidence that the prophylaxis threshold is an economic boundary that has moved (40% → 20%) without any change in the biology |
| 100 | Cosler LE, et al. *Effect of outpatient treatment of febrile neutropenia on the risk threshold for the use of CSF in patients with cancer treated with chemotherapy.* Value Health 2005. [PMID 15841893](https://pubmed.ncbi.nlm.nih.gov/15841893/) | The threshold moves when the *cost of an FN episode* moves — used in A14 |

---

## 11. Dose intensity and tumour control

The model's weakest layer, and these are the papers it leans on. The
association between relative dose intensity and survival is largely
observational; the model's tumour module exists only to price dose reduction
against FN risk, and every result that uses it says so.

| # | Reference | Used for |
|---|-----------|----------|
| 101 | Hryniuk W, Levine MN. *Analysis of dose intensity for adjuvant chemotherapy trials in stage II breast cancer.* J Clin Oncol 1986. [PMID 3525765](https://pubmed.ncbi.nlm.nih.gov/3525765/) | The definition of relative dose intensity the model computes |
| 102 | Wood WC, et al. *Dose and dose intensity of adjuvant chemotherapy for stage II, node-positive breast carcinoma.* N Engl J Med 1994. [PMID 8080512](https://pubmed.ncbi.nlm.nih.gov/8080512/) | CALGB 8541: dose–response in the adjuvant setting |
| 103 | Citron ML, et al. *Randomized trial of dose-dense versus conventionally scheduled and sequential versus concurrent combination chemotherapy as postoperative adjuvant treatment of node-positive primary breast cancer: first report of Intergroup Trial C9741/Cancer and Leukemia Group B Trial 9741.* J Clin Oncol 2003. [PMID 12668651](https://pubmed.ncbi.nlm.nih.gov/12668651/) | Dose density made possible *only* by G-CSF — the ddAC regimen in the model |
| 104 | Norton L, Simon R. *Tumor size, sensitivity to therapy, and design of treatment schedules.* Cancer Treat Rep 1977. [PMID 589597](https://pubmed.ncbi.nlm.nih.gov/589597/) | The interval-regrowth term that makes cycle DELAY cost something distinct from dose REDUCTION |
| 105 | Simon R, Norton L. *The Norton-Simon hypothesis: designing more effective and less toxic chemotherapeutic regimens.* Nat Clin Pract Oncol 2006. [PMID 16894366](https://pubmed.ncbi.nlm.nih.gov/16894366/) | Restatement of the same principle |
| 106 | Weycker D, et al. *Incidence of reduced chemotherapy relative dose intensity among women with early stage breast cancer in US clinical practice.* Breast Cancer Res Treat 2012. [PMID 22270932](https://pubmed.ncbi.nlm.nih.gov/22270932/) | How often RDI is actually lost — the base rate the A12 controller should reproduce |
| 107 | Bonadonna G, et al. *Combined modality approach for high-risk breast cancer. The Milan Cancer Institute experience.* Surg Oncol Clin N Am 1995. [PMID 8535906](https://pubmed.ncbi.nlm.nih.gov/8535906/) | Long-term dose-intensity follow-up |
| 108 | Levin L, Hryniuk W. *Importance of multiagent chemotherapy regimens in ovarian carcinoma: dose intensity analysis.* J Natl Cancer Inst 1993. [PMID 8411257](https://pubmed.ncbi.nlm.nih.gov/8411257/) | Dose intensity in a second tumour type |

---

## 12. Platelet and erythroid lineages

| # | Reference | Used for |
|---|-----------|----------|
| 109 | Soff GA, et al. *Romiplostim Treatment of Chemotherapy-Induced Thrombocytopenia.* J Clin Oncol 2019. [PMID 31545663](https://pubmed.ncbi.nlm.nih.gov/31545663/) | The `ROMI` intervention node and the c-Mpl axis |
| 110 | Soff GA, et al. *Romiplostim in chemotherapy-induced thrombocytopenia: A review of the literature.* Cancer Med 2024. [PMID 39135303](https://pubmed.ncbi.nlm.nih.gov/39135303/) | Current evidence base |
| 111 | Song AB, Al-Samkari H. *Chemotherapy-induced thrombocytopenia: modern diagnosis and treatment.* Br J Haematol 2025. [PMID 40040262](https://pubmed.ncbi.nlm.nih.gov/40040262/) | Definitions and thresholds used for the platelet endpoints |
| 112 | Elting LS, et al. *Incidence, cost, and outcomes of bleeding and chemotherapy dose modification among solid tumor patients with chemotherapy-induced thrombocytopenia.* J Clin Oncol 2001. [PMID 11181679](https://pubmed.ncbi.nlm.nih.gov/11181679/) | Platelet-driven dose modification — the second gate in the A12 controller |

---

## 13. Count versus function — corticosteroid and the ANC that lies

| # | Reference | Used for |
|---|-----------|----------|
| 113 | Nakagawa M, et al. *Glucocorticoid-induced granulocytosis: contribution of marrow release and demargination of intravascular granulocytes.* Circulation 1998. [PMID 9826319](https://pubmed.ncbi.nlm.nih.gov/9826319/) | Quantitative split of the steroid-induced ANC rise into marrow release and demargination — the two terms the model implements, and the direct source of the claim that the day-1 CBC gate can be passed for the wrong reason |
| 114 | Dale DC, et al. *Comparison of agents producing a neutrophilic leukocytosis in man. Hydrocortisone, prednisone, endotoxin, and etiocholanolone.* J Clin Invest 1975. [PMID 1159089](https://pubmed.ncbi.nlm.nih.gov/1159089/) | Magnitude and time course of the corticosteroid ANC rise (`dex_demarg`) |

---

## 14. Parameter provenance table

Every numeric parameter in `cin_reference_impl.py` and `cin_mrgsolve_model.R`
falls into exactly one of four classes. This table is the audit trail.

| Class | Parameters | Source |
|-------|-----------|--------|
| **Published, held fixed** | `MTT_N` 89.3 h, `gamma_N` 0.170 | Friberg 2002 [12488418] |
| | `tau_circ_N` 7.9 h, `marg_ratio` 1.0 | Dancey 1976 [956397]; Athens 1961 [13684958] |
| | `MTT_P` 200 h, `life_P` 240 h, `MTT_E` 150 h, `life_E` 2880 h | standard lineage kinetics; Krzyzanski 2015 [26341875] |
| | docetaxel `CL/V1/V2/V3/Q2/Q3` | Bruno 1998 [9440742] |
| | paclitaxel PK + `Cthr` 0.05 µM | Henningsson 2001 [11600609]; Joerger 2007 [17975154] |
| | carboplatin `CL = GFR + 25` | Calvert 1989 [2681557] |
| | filgrastim `V`, `ka`, `CL_lin` | Wiczling 2009 [19902989]; Krzyzanski 2010 [20881223] |
| | pegfilgrastim `ka`, `CL_lin ≈ 0`, `Vmax_nm`, `Km_nm` | Roskos 2006 [16809800]; Yang 2011 [21456630]; Molineux 2003 [12679729] |
| | `Emax_surv` 0.25 (t½ 7.9 → 10.5 h) | anti-apoptotic effect, Bendall 2014 [25131807] |
| | `Emax_demarg` = 0 for G-CSF | Nakagawa 1998 [9826319] attributes demargination to corticosteroid, not G-CSF |
| **Derived from structure** | `phi` = 2⁻⁷ | mitotic amplification 2⁵–2⁹; the exponent is pinned by the requirement that the ANC floor lie below the deepest observed nadirs |
| | `Emax_amp` initial value 0.357 | solves `((1−φ)(1+E))^(1/γ) = 5` for a 5-fold saturated response |
| | `ktr = 4/MTT`, `Prol0 = Circ0·kcirc/ktr` | steady-state algebra of the transit chain |
| **Fitted, one per anchor** | `sens_global` | CAE placebo DSN 6.0 d — Crawford 1991 [1711156] |
| | `sens_taxane` | AT + pegfilgrastim DSN 1.8 d — Green 2003 [12488289] |
| | `GPD["Emax_amp"]` (final) | CAE + filgrastim DSN 1.0 d — Crawford 1991 [1711156] |
| | `hFN_p`, `hFN_k` | FN 17% / 1% — Vogel 2005 [15718314] |
| | per-agent `slope_P` | typical monotherapy platelet nadirs |
| **Assumed, sensitivity-tested** | `khsc_kill`, `kfib` | no direct human anchor; set so cumulative nadir deepening across six cycles is modest, and swept in A9 |
| | `kdam_bar`, `kheal_bar`, `k_translocate`, `kgrow_B`, `kkill_B`, `B50` | mucosal/infection layer; internally consistent, not independently identified. The FN hazard is fitted, so these affect the *shape* of the hazard, not its level |
| | `logkill`, `Tdouble`, `res_kill_frac` | tumour layer — explicitly the weakest assumption in the model |

---

## 15. Literature the model does not reproduce

Recorded rather than tuned away.

| Observation | Source | Model | Comment |
|-------------|--------|-------|---------|
| Filgrastim 5 µg/kg/day raises the ANC 4–6× in an intact marrow | Lord 1989 [2480603]; Dale 2018 [28939926] | under-predicts (see A1) | The model's G-CSF production effect is bounded by γ; to reach a 5× rise the amplification must be large enough to also abolish post-chemotherapy severe neutropenia, which the anchors forbid. The model therefore chooses to fit the *therapeutic* setting and miss the *healthy-volunteer* setting, and says so |
| Carboplatin's ANC nadir is late (≈ day 21) | Schmitt 2010 [20855828] | nadir ≈ day 8 | The model's nadir day is a function of MTT alone, so it cannot produce an agent-specific nadir day. A structural limitation of the Friberg topology, not a parameter problem |
| Trial FN rates understate real-world FN rates | Truong 2016 [26712901] | fitted to the trial rate | The FN hazard is calibrated to Vogel 2005's 17%, so every absolute P(FN) in this model is a trial-population number and will read low against registry data |
| Relative dose intensity ↔ survival | Lyman 2013 [23788754] | implemented as a log-kill | The association is observational and confounded; the tumour layer should be read as a *price list* for dose reduction, not a survival prediction |
| Trilaciclib DSN magnitude | Weiss 2021 [33895103] | prediction only | The pooled abstract states the endpoints but not their values, so the model's trilaciclib numbers are unvalidated in magnitude |

---

*All PMIDs in this file were resolved and verified through the NCBI E-utilities
API during construction of the model. Titles, first authors, journals and years
are the API's records.*
