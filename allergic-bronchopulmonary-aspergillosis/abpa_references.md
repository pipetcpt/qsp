# Allergic Bronchopulmonary Aspergillosis (ABPA) — annotated references

**Every PMID below was resolved against the PubMed E-utilities API on 2026-07-29**, and
the first author, year, journal, volume, pages and title were written into this file
programmatically from the API response. No citation here is asserted from recall.

Each entry states what the paper supplies to the model — a parameter, a calibration
anchor, a structural decision, or a claim the model contradicts. Papers read only for
background are not listed.

Two conventions are used deliberately:

- Where a number was **fitted**, the reference is named as the fit target and the
  analysis that consumed it is named too (see §⑥, Varis 1998 → A4).
- Where the model **predicts** a published number it was never shown, that is stated —
  including when the prediction is wrong (see §⑥, Varis 2000: model 1.12× vs observed
  1.24×).

---

## ① Disease definition, diagnostic criteria and staging — what the model must reproduce

- **Agarwal R (2024)** Revised ISHAM-ABPA working group clinical practice guidelines for diagnosing, classifying and treating allergic bronchopulmonary aspergillosis/mycoses. *Eur Respir J* 63. [PMID 38423624](https://pubmed.ncbi.nlm.nih.gov/38423624/)  
  The current reference standard. Supplies the diagnostic thresholds hard-coded in the model's readouts (Af-specific IgE > 0.35 kUA/L, total IgE > 500–1000 IU/mL, eosinophils > 500/µL) and the ISHAM medium-dose prednisolone taper implemented verbatim in `pred_taper()`.

- **Asano K (2021)** New clinical diagnostic criteria for allergic bronchopulmonary aspergillosis/mycosis and its validation. *J Allergy Clin Immunol* 147:1261-1268.e5. [PMID 32920094](https://pubmed.ncbi.nlm.nih.gov/32920094/)  
  The Japanese ABPM criteria, which weight mucus plugging far more heavily than the ISHAM set. This is the clinical counterpart of the model's decision to make PLUG a state variable rather than a symptom.

- **Agarwal R (2014)** Cut-off values of serum IgE (total and A. fumigatus -specific) and eosinophil count in differentiating allergic bronchopulmonary aspergillosis from asthma. *Mycoses* 57:659-63. [PMID 24963741](https://pubmed.ncbi.nlm.nih.gov/24963741/)  
  Source of the total-IgE and eosinophil cut-offs used to place the untreated endotype at total IgE 2000 IU/mL and blood eosinophils 800/µL.

- **Patterson R (1987)** Prolonged evaluation of patients with corticosteroid-dependent asthma stage of allergic bronchopulmonary aspergillosis. *J Allergy Clin Immunol* 80:663-8. [PMID 3316345](https://pubmed.ncbi.nlm.nih.gov/3316345/)  
  Patterson's original staging, including the corticosteroid-dependent stage IV that the model's `pred_maint` arm represents.

- **Greenberger PA (1987)** Allergic bronchopulmonary aspergillosis. Model of bronchopulmonary disease with defined serologic, radiologic, pathologic and clinical findings from asthma to fatal destructive lung disease. *Chest* 91:165S-171S. [PMID 3556066](https://pubmed.ncbi.nlm.nih.gov/3556066/)  
  The serologic/radiologic/immunologic definition the model's endpoint block is built to output.

- **Lee TM (1987)** Stage V (fibrotic) allergic bronchopulmonary aspergillosis. A review of 17 cases followed from diagnosis. *Arch Intern Med* 147:319-23. [PMID 3545117](https://pubmed.ncbi.nlm.nih.gov/3545117/)  
  Stage V (fibrotic) ABPA. The clinical destination of the model's BRON integrator, and the reason BRON is written so that it can never fall.

- **Sehgal IS (2024)** EQUAL ABPA Score 2024: A Tool to Measure Guideline Adherence for Managing Allergic Bronchopulmonary Aspergillosis. *Mycoses* 67:e13810. [PMID 39462638](https://pubmed.ncbi.nlm.nih.gov/39462638/)  
  EQUAL ABPA score — guideline adherence as a measurable quantity, which is what the model's scenario protocols are attempting to encode.

- **Tillie-Leblond I (2005)** Allergic bronchopulmonary aspergillosis. *Allergy* 60:1004-13. [PMID 15969680](https://pubmed.ncbi.nlm.nih.gov/15969680/)  
  General review; used for the qualitative wiring of clusters ②–⑦ of the mechanistic map.

- **Denning DW (2013)** Global burden of allergic bronchopulmonary aspergillosis with asthma and its complication chronic pulmonary aspergillosis in adults. *Med Mycol* 51:361-70. [PMID 23210682](https://pubmed.ncbi.nlm.nih.gov/23210682/)  
  Global burden: several million ABPA cases among asthmatics worldwide. Motivates the choice of disease.

- **Agarwal R (2023)** Prevalence of Aspergillus sensitization and Allergic Bronchopulmonary Aspergillosis in bronchial asthma: A systematic review of Indian studies. *Lung India* 40:527-536. [PMID 37961961](https://pubmed.ncbi.nlm.nih.gov/37961961/)  
  Prevalence of Aspergillus sensitisation and ABPA in asthma cohorts.

- **Chakrabarti A (2002)** Eight-year study of allergic bronchopulmonary aspergillosis in an Indian teaching hospital. *Mycoses* 45:295-9. [PMID 12572718](https://pubmed.ncbi.nlm.nih.gov/12572718/)  
  Eight-year Indian cohort; the clinical spectrum the untreated natural-history arm (A1) is checked against.

- **Sunman B (2020)** Current Approach in the Diagnosis and Management of Allergic Bronchopulmonary Aspergillosis in Children With Cystic Fibrosis. *Front Pediatr* 8:582964. [PMID 33194914](https://pubmed.ncbi.nlm.nih.gov/33194914/)  
  ABPA in cystic fibrosis in children; the CFTR node in cluster ②.

- **Stevens DA (2003)** Allergic bronchopulmonary aspergillosis in cystic fibrosis--state of the art: Cystic Fibrosis Foundation Consensus Conference. *Clin Infect Dis* 37 Suppl 3:S225-64. [PMID 12975753](https://pubmed.ncbi.nlm.nih.gov/12975753/)  
  Cystic Fibrosis Foundation consensus on ABPA in CF — the alternative predisposing condition.

---

## ② The mucus plug as an object — the evidence the sanctuary compartment rests on

- **Agarwal R (2010)** High attenuation mucoid impaction in allergic bronchopulmonary aspergillosis. *World J Radiol* 2:41-3. [PMID 21160739](https://pubmed.ncbi.nlm.nih.gov/21160739/)  
  High-attenuation mucus in ABPA: mucus denser than paraspinal muscle on CT. The radiological existence proof for a distinct plug compartment, and the reason PLUG is drawn as an integrating compartment rather than a flux.

- **Tang M (2022)** Mucus Plugs Persist in Asthma, and Changes in Mucus Plugs Associate with Changes in Airflow over Time. *Am J Respir Crit Care Med* 205:1036-1045. [PMID 35104436](https://pubmed.ncbi.nlm.nih.gov/35104436/)  
  Mucus plugs PERSIST over years, and changes in plug score track changes in airflow. The single most important observation behind treating the plug as a slow state variable with its own clearance rate rather than as a fast symptom.

- **Bonser LR (2016)** Epithelial tethering of MUC5AC-rich mucus impairs mucociliary transport in asthma. *J Clin Invest* 126:2367-71. [PMID 27183390](https://pubmed.ncbi.nlm.nih.gov/27183390/)  
  Epithelial tethering of MUC5AC-rich mucus impairs mucociliary transport. The mechanistic basis for modelling reduced `kout0` rather than merely raised mucin production.

- **Bonser LR (2017)** Airway Mucus and Asthma: The Role of MUC5AC and MUC5B. *J Clin Med* 6. [PMID 29186064](https://pubmed.ncbi.nlm.nih.gov/29186064/)  
  MUC5AC versus MUC5B division of labour; the two mucin nodes in cluster ⑧.

- **Kuperman DA (2002)** Direct effects of interleukin-13 on epithelial cells cause airway hyperreactivity and mucus overproduction in asthma. *Nat Med* 8:885-9. [PMID 12091879](https://pubmed.ncbi.nlm.nih.gov/12091879/)  
  IL-13 acting directly on epithelium to cause mucus overproduction and hyper-responsiveness. The `smuc` and `ca_fev` terms come from here, and it is why dupilumab (anti-IL-4Rα) is the strongest single agent on PLUG in scenario 14.

- **Persson EK (2019)** Protein crystallization promotes type 2 immunity and is reversible by antibody treatment. *Science* 364. [PMID 31123109](https://pubmed.ncbi.nlm.nih.gov/31123109/)  
  Galectin-10 crystallisation drives type-2 immunity and is REVERSIBLE BY ANTIBODY. Charcot–Leyden crystals as a druggable target, and direct evidence that the plug's physical state can be changed pharmacologically — i.e. that `f_pen` and `k_out` are not constants of nature.

- **Rodríguez-Alcázar JF (2019)** Charcot-Leyden Crystals Activate the NLRP3 Inflammasome and Cause IL-1β Inflammation in Human Macrophages. *J Immunol* 202:550-558. [PMID 30559319](https://pubmed.ncbi.nlm.nih.gov/30559319/)  
  Charcot–Leyden crystals activate NLRP3. The link from crystallised galectin-10 back to innate inflammation (cluster ③).

- **Tang M (2024)** Utility of eosinophil peroxidase as a biomarker of eosinophilic inflammation in asthma. *J Allergy Clin Immunol* 154:580-591.e6. [PMID 38663815](https://pubmed.ncbi.nlm.nih.gov/38663815/)  
  Eosinophil peroxidase as a biomarker of eosinophilic inflammation — the measurable counterpart of the model's EPX state, which is what couples eosinophils to plug clearance through `g_epx`.

- **Everman JL (2024)** A common polymorphism in the Intelectin-1 gene influences mucus plugging in severe asthma. *Nat Commun* 15:3900. [PMID 38724552](https://pubmed.ncbi.nlm.nih.gov/38724552/)  
  Intelectin-1 polymorphism influences mucus plugging: genetic evidence that plug formation is partly separable from type-2 inflammation, which is exactly what giving PLUG its own compartment assumes.

- **Aegerter H (2026)** Effectiveness of biologics for reducing occlusive mucus plugs in patients with severe asthma: a systematic review. *Respir Res* 27:69. [PMID 41559746](https://pubmed.ncbi.nlm.nih.gov/41559746/)  
  Biologics reduce occlusive mucus plugs in severe asthma (systematic review). Clinical evidence that the `BIO_MUC → PLUG_CLEAR` edge is real.

---

## ③ Aspergillus biology, biofilm and the penetration barrier — where f_pen comes from

- **Latgé JP (2019)** Aspergillus fumigatus and Aspergillosis in 2019. *Clin Microbiol Rev* 33. [PMID 31722890](https://pubmed.ncbi.nlm.nih.gov/31722890/)  
  Comprehensive review of A. fumigatus and aspergillosis; source for cluster ①'s germination sequence and secreted products.

- **Latgé JP (2017)** The Cell Wall of the Human Fungal Pathogen Aspergillus fumigatus: Biosynthesis, Organization, Immune Response, and Virulence. *Annu Rev Microbiol* 71:99-116. [PMID 28701066](https://pubmed.ncbi.nlm.nih.gov/28701066/)  
  The cell wall: β-glucan and chitin unmasking on swelling, which is why the map draws `SWELLING → DECTIN1` rather than `CONIDIA → DECTIN1`.

- **Loussert C (2010)** In vivo biofilm composition of Aspergillus fumigatus. *Cell Microbiol* 12:405-10. [PMID 19889082](https://pubmed.ncbi.nlm.nih.gov/19889082/)  
  In vivo biofilm composition of A. fumigatus — galactosaminogalactan-containing extracellular matrix in tissue, not only in vitro.

- **Seidler MJ (2008)** Aspergillus fumigatus forms biofilms with reduced antifungal drug susceptibility on bronchial epithelial cells. *Antimicrob Agents Chemother* 52:4130-6. [PMID 18710910](https://pubmed.ncbi.nlm.nih.gov/18710910/)  
  A. fumigatus forms biofilms on bronchial epithelial cells WITH REDUCED ANTIFUNGAL DRUG SUSCEPTIBILITY. With the tissue-penetration review below, this is the primary justification for f_pen ≈ 0.1 rather than f_pen = 1.

- **Felton T (2014)** Tissue penetration of antifungal agents. *Clin Microbiol Rev* 27:68-88. [PMID 24396137](https://pubmed.ncbi.nlm.nih.gov/24396137/)  
  Tissue penetration of antifungal agents. The quantitative anchor for how much less drug reaches a poorly perfused, matrix-rich compartment.

- **Ashbee HR (2014)** Therapeutic drug monitoring (TDM) of antifungal agents: guidelines from the British Society for Medical Mycology. *J Antimicrob Chemother* 69:1162-76. [PMID 24379304](https://pubmed.ncbi.nlm.nih.gov/24379304/)  
  BSMM therapeutic drug monitoring guidelines: itraconazole trough targets (parent plus hydroxy-metabolite) used to check that the model's 200 mg BID steady state lands in the therapeutic range.

- **Siopi M (2014)** Susceptibility breakpoints and target values for therapeutic drug monitoring of voriconazole and Aspergillus fumigatus in an in vitro pharmacokinetic/pharmacodynamic model. *J Antimicrob Chemother* 69:1611-9. [PMID 24550381](https://pubmed.ncbi.nlm.nih.gov/24550381/)  
  Voriconazole susceptibility breakpoints and TDM target values against A. fumigatus.

- **Howard SJ (2009)** Frequency and evolution of Azole resistance in Aspergillus fumigatus associated with treatment failure. *Emerg Infect Dis* 15:1068-76. [PMID 19624922](https://pubmed.ncbi.nlm.nih.gov/19624922/)  
  Frequency and evolution of azole resistance in A. fumigatus associated with treatment failure. Supports the `RESIST_SEL` edge from sub-MIC intra-plug exposure.

- **Chowdhary A (2014)** Multi-azole-resistant Aspergillus fumigatus in the environment in Tanzania. *J Antimicrob Chemother* 69:2979-83. [PMID 25006238](https://pubmed.ncbi.nlm.nih.gov/25006238/)  
  Environmental multi-azole-resistant A. fumigatus; the `CYP51_MUT` node.

- **Agbetile J (2012)** Isolation of filamentous fungi from sputum in asthma is associated with reduced post-bronchodilator FEV1. *Clin Exp Allergy* 42:782-91. [PMID 22515394](https://pubmed.ncbi.nlm.nih.gov/22515394/)  
  Filamentous fungi in sputum in asthma associated with reduced post-bronchodilator FEV₁ — the clinical link from fungal burden to fixed obstruction.

- **Rick EM (2020)** The airway fungal microbiome in asthma. *Clin Exp Allergy* 50:1325-1341. [PMID 32808353](https://pubmed.ncbi.nlm.nih.gov/32808353/)  
  The airway fungal microbiome in asthma; context for the `seed` term, which makes true sterilisation impossible by construction because conidia keep arriving.

- **Ozyigit LP (2021)** Fungal bronchitis is a distinct clinical entity which is responsive to antifungal therapy. *Chron Respir Dis* 18:1479973120964448. [PMID 33472416](https://pubmed.ncbi.nlm.nih.gov/33472416/)  
  Fungal bronchitis as a distinct antifungal-responsive entity — the phenotype in which the luminal compartment dominates and the sanctuary matters least.

---

## ④ Epithelium, proteases and the type-2 amplifier

- **Kauffman HF (2000)** Protease-dependent activation of epithelial cells by fungal allergens leads to morphologic changes and cytokine production. *J Allergy Clin Immunol* 105:1185-93. [PMID 10856154](https://pubmed.ncbi.nlm.nih.gov/10856154/)  
  Protease-dependent activation of epithelial cells by fungal allergens. The `PROTEASE → PAR2 → EPI` chain.

- **Tomee JF (1997)** Proteases from Aspergillus fumigatus induce release of proinflammatory cytokines and cell detachment in airway epithelial cell lines. *J Infect Dis* 176:300-3. [PMID 9207388](https://pubmed.ncbi.nlm.nih.gov/9207388/)  
  Aspergillus proteases release pro-inflammatory cytokines and cause epithelial cell detachment.

- **Borger P (1999)** Proteases from Aspergillus fumigatus induce interleukin (IL)-6 and IL-8 production in airway epithelial cell lines by transcriptional mechanisms. *J Infect Dis* 180:1267-74. [PMID 10479157](https://pubmed.ncbi.nlm.nih.gov/10479157/)  
  Aspergillus proteases induce IL-6 and IL-8 in airway epithelium.

- **Homma T (2016)** Role of Aspergillus fumigatus in Triggering Protease-Activated Receptor-2 in Airway Epithelial Cells and Skewing the Cells toward a T-helper 2 Bias. *Am J Respir Cell Mol Biol* 54:60-70. [PMID 26072921](https://pubmed.ncbi.nlm.nih.gov/26072921/)  
  A. fumigatus triggering protease-activated receptor 2 in airway epithelial cells — the specific receptor drawn in the map.

- **Redes JL (2019)** Aspergillus fumigatus-Secreted Alkaline Protease 1 Mediates Airways Hyperresponsiveness in Severe Asthma. *Immunohorizons* 3:368-377. [PMID 31603851](https://pubmed.ncbi.nlm.nih.gov/31603851/)  
  Alp1 (Asp f13) mediates airway hyper-responsiveness in severe asthma: a named allergen with a named physiological consequence.

---

## ⑤ Itraconazole and voriconazole pharmacokinetics — the nonlinearity the model implements

- **Poirier JM (1998)** Optimisation of itraconazole therapy using target drug concentrations. *Clin Pharmacokinet* 35:461-73. [PMID 9884817](https://pubmed.ncbi.nlm.nih.gov/9884817/)  
  Optimisation of itraconazole therapy using target drug concentrations; parent plus hydroxy-itraconazole, and the ~1.5–2× metabolite/parent ratio the model reproduces (1.72× at 200 mg BID).

- **Barone JA (1998)** Enhanced bioavailability of itraconazole in hydroxypropyl-beta-cyclodextrin solution versus capsules in healthy volunteers. *Antimicrob Agents Chemother* 42:1862-5. [PMID 9661037](https://pubmed.ncbi.nlm.nih.gov/9661037/)  
  Enhanced bioavailability of the cyclodextrin oral solution over capsules — the `F_itra` / `ITRA_ABS` node.

- **Barone JA (1998)** Food interaction and steady-state pharmacokinetics of itraconazole oral solution in healthy volunteers. *Pharmacotherapy* 18:295-301. [PMID 9545149](https://pubmed.ncbi.nlm.nih.gov/9545149/)  
  Food interaction and steady-state pharmacokinetics of itraconazole oral solution.

- **Purkins L (2003)** Voriconazole, a novel wide-spectrum triazole: oral pharmacokinetics and safety. *Br J Clin Pharmacol* 56 Suppl 1:10-6. [PMID 14616408](https://pubmed.ncbi.nlm.nih.gov/14616408/)  
  Voriconazole oral PK: the saturable, supra-proportional behaviour the model implements as mixed Michaelis–Menten plus linear clearance. Pure Michaelis–Menten was tried first and diverges whenever the input rate exceeds Vmax, which is not what CYP2C19 poor metabolisers actually do — the mixed form was adopted for that reason and the reason is recorded in the code.

- **Agarwal R (2026)** Safety and efficacy of inhaled itraconazole in adults with asthma and allergic bronchopulmonary aspergillosis (PUR1900-ABPA): a randomized, double-blind, parallel group, placebo-controlled, multicenter, phase 2 trial. *Eur Respir J*. [PMID 42463270](https://pubmed.ncbi.nlm.nih.gov/42463270/)  
  Inhaled itraconazole in asthma with ABPA. The route that raises luminal exposure without raising the systemic CYP3A4 interaction — A5 is, in effect, an argument for why that separation matters.

- **Bergagnini-Kolev M (2023)** Evaluation of the Potential for Drug-Drug Interactions with Inhaled Itraconazole Using Physiologically Based Pharmacokinetic Modelling, Based on Phase 1 Clinical Data. *AAPS J* 25:62. [PMID 37344751](https://pubmed.ncbi.nlm.nih.gov/37344751/)  
  PBPK evaluation of drug–drug interaction potential for INHALED itraconazole; directly relevant to whether changing route escapes the confound A5 describes.

---

## ⑥ The CYP3A4 interaction — the numbers A4 is calibrated and validated against

- **Varis T (1998)** Plasma concentrations and effects of oral methylprednisolone are considerably increased by itraconazole. *Clin Pharmacol Ther* 64:363-8. [PMID 9797792](https://pubmed.ncbi.nlm.nih.gov/9797792/)  
  Oral methylprednisolone concentrations and effects are CONSIDERABLY INCREASED by itraconazole (~2.6× AUC). This is the ONE interaction magnitude fitted in A4; everything else in that analysis is prediction.

- **Varis T (1999)** Itraconazole decreases the clearance and enhances the effects of intravenously administered methylprednisolone in healthy volunteers. *Pharmacol Toxicol* 85:29-32. [PMID 10426160](https://pubmed.ncbi.nlm.nih.gov/10426160/)  
  Itraconazole decreases the clearance of intravenously administered methylprednisolone — separates the hepatic component from the gut-wall component, which the model needs in order to split `s3A4h` from `s3A4g`.

- **Varis T (2000)** The effect of itraconazole on the pharmacokinetics and pharmacodynamics of oral prednisolone. *Eur J Clin Pharmacol* 56:57-60. [PMID 10853878](https://pubmed.ncbi.nlm.nih.gov/10853878/)  
  The effect of itraconazole on oral PREDNISOLONE (~1.24× AUC). Not fitted: the model predicts 1.12×, under-predicting by ~10%, and that discrepancy is reported rather than tuned away.

- **Raaska K (2002)** Plasma concentrations of inhaled budesonide and its effects on plasma cortisol are increased by the cytochrome P4503A4 inhibitor itraconazole. *Clin Pharmacol Ther* 72:362-9. [PMID 12386638](https://pubmed.ncbi.nlm.nih.gov/12386638/)  
  Plasma concentrations of INHALED budesonide and its effects on plasma cortisol are increased by itraconazole (~4.2×). Used in A4 to show that a single-site interaction model is arithmetically IMPOSSIBLE here — the ceiling of any single-site model is 1/I = 3.4 at this inhibition level, so two sequential sites (gut wall and liver) are required.

- **Skov M (2002)** Iatrogenic adrenal insufficiency as a side-effect of combined treatment of itraconazole and budesonide. *Eur Respir J* 20:127-33. [PMID 12166560](https://pubmed.ncbi.nlm.nih.gov/12166560/)  
  Iatrogenic adrenal insufficiency from combined itraconazole and budesonide. The clinical event A10 reproduces from the PK block alone, with nothing in the disease model changed.

- **Main KM (2002)** Cushing's syndrome due to pharmacological interaction in a cystic fibrosis patient. *Acta Paediatr* 91:1008-11. [PMID 12412882](https://pubmed.ncbi.nlm.nih.gov/12412882/)  
  Cushing's syndrome from the same pharmacological interaction in a cystic fibrosis patient.

- **Blondin MC (2013)** Iatrogenic Cushing syndrome in patients receiving inhaled budesonide and itraconazole or ritonavir: two cases and literature review. *Endocr Pract* 19:e138-41. [PMID 23807527](https://pubmed.ncbi.nlm.nih.gov/23807527/)  
  Iatrogenic Cushing syndrome in patients receiving inhaled budesonide with itraconazole or ritonavir.

- **Garg VK (2023)** Iatrogenic Cushing's syndrome in a case of allergic bronchopulmonary aspergillosis treated with oral itraconazole and inhaled budesonide. *BMJ Case Rep* 16. [PMID 37813554](https://pubmed.ncbi.nlm.nih.gov/37813554/)  
  Iatrogenic Cushing's syndrome in a case of ABPA treated with itraconazole — A10's exact scenario, in the exact disease.

- **Templeton IE (2024)** Creation of Novel Sensitive Probe Substrate and Moderate Inhibitor Models for a Comprehensive Prediction of CYP2C8 Interactions for Tucatinib. *Clin Pharmacol Ther* 115:299-308. [PMID 37971208](https://pubmed.ncbi.nlm.nih.gov/37971208/)  
  Sensitive probe substrate and moderate inhibitor models for CYP3A interaction prediction; methodological context for expressing the interaction as a single fractional-activity term I.

---

## ⑦ Corticosteroid PK/PD and its toxicity currency

- **Czock D (2005)** Pharmacokinetics and pharmacodynamics of systemically administered glucocorticoids. *Clin Pharmacokinet* 44:61-98. [PMID 15634032](https://pubmed.ncbi.nlm.nih.gov/15634032/)  
  Pharmacokinetics and pharmacodynamics of systemically administered glucocorticoids. Source of the potency ratios (`pot_mpred`, `pot_bud`) and the unbound fractions.

- **Meibohm B (1999)** A pharmacokinetic/pharmacodynamic approach to predict the cumulative cortisol suppression of inhaled corticosteroids. *J Pharmacokinet Biopharm* 27:127-47. [PMID 10567952](https://pubmed.ncbi.nlm.nih.gov/10567952/)  
  A PK/PD approach predicting CUMULATIVE CORTISOL SUPPRESSION by inhaled corticosteroids — the structural template for the model's CORT compartment and `Imax_c`.

- **Boudinot FD (1986)** Receptor-mediated pharmacodynamics of prednisolone in the rat. *J Pharmacokinet Biopharm* 14:469-93. [PMID 2879901](https://pubmed.ncbi.nlm.nih.gov/2879901/)  
  Receptor-mediated pharmacodynamics of prednisolone: the classical GR-occupancy formulation that the model's `CS` term is a reduced form of.

- **Ayyar VS (2019)** Modeling Corticosteroid Pharmacokinetics and Pharmacodynamics, Part II: Sex Differences in Methylprednisolone Pharmacokinetics and Corticosterone Suppression. *J Pharmacol Exp Ther* 370:327-336. [PMID 31197019](https://pubmed.ncbi.nlm.nih.gov/31197019/)  
  Modern corticosteroid PK/PD modelling (methylprednisolone); confirms the receptor-mediated structure at the level of detail this model abstracts away.

- **Van Staa TP (2000)** The use of a large pharmacoepidemiological database to study exposure to oral corticosteroids and risk of fractures: validation of study population and results. *Pharmacoepidemiol Drug Saf* 9:359-66. [PMID 19025840](https://pubmed.ncbi.nlm.nih.gov/19025840/)  
  Oral corticosteroid exposure and fracture risk in a large pharmacoepidemiological database. The `kbmd` term and the fracture endpoint.

- **Waljee AK (2017)** Short term use of oral corticosteroids and related harms among adults in the United States: population based cohort study. *BMJ* 357:j1415. [PMID 28404617](https://pubmed.ncbi.nlm.nih.gov/28404617/)  
  Short-term oral corticosteroid use and related harms — the reason cumulative prednisolone-equivalent dose is carried as a state (CUMO) rather than as an afterthought.

- **Rogers MAM (2018)** Longitudinal study of short-term corticosteroid use by working-age adults with diabetes mellitus: Risks and mitigating factors. *J Diabetes* 10:546-555. [PMID 29193668](https://pubmed.ncbi.nlm.nih.gov/29193668/)  
  Short-term corticosteroid use and glycaemic risk; the HBA1C compartment.

- **Agarwal R (2023)** Long-term follow-up of allergic bronchopulmonary aspergillosis treated with glucocorticoids: A study of 182 subjects. *Mycoses* 66:953-959. [PMID 37555291](https://pubmed.ncbi.nlm.nih.gov/37555291/)  
  Long-term follow-up of ABPA treated with glucocorticoids. The observed relapse-after-taper behaviour that the model's ISHAM-taper arm reproduces, with total IgE back near baseline by week 52.

---

## ⑧ Azole and antifungal trials in ABPA — what the scenario arms are compared against

- **Stevens DA (2000)** A randomized trial of itraconazole in allergic bronchopulmonary aspergillosis. *N Engl J Med* 342:756-62. [PMID 10717010](https://pubmed.ncbi.nlm.nih.gov/10717010/)  
  The randomised trial that established itraconazole in ABPA. Its partial response rate is the number A2's knife-edge threshold result is read against.

- **Wark PA (2003)** Azoles for allergic bronchopulmonary aspergillosis associated with asthma. *Cochrane Database Syst Rev*:CD001108. [PMID 12917898](https://pubmed.ncbi.nlm.nih.gov/12917898/)  
  Cochrane review of azoles for ABPA.

- **Agarwal R (2018)** A Randomized Trial of Itraconazole vs Prednisolone in Acute-Stage Allergic Bronchopulmonary Aspergillosis Complicating Asthma. *Chest* 153:656-664. [PMID 29331473](https://pubmed.ncbi.nlm.nih.gov/29331473/)  
  Itraconazole versus prednisolone in acute-stage ABPA: the head-to-head that makes the model's scenarios 02 and 04 directly comparable.

- **Agarwal R (2022)** A randomised trial of prednisolone versus prednisolone and itraconazole in acute-stage allergic bronchopulmonary aspergillosis complicating asthma. *Eur Respir J* 59. [PMID 34503983](https://pubmed.ncbi.nlm.nih.gov/34503983/)  
  Prednisolone versus prednisolone PLUS itraconazole in acute-stage ABPA. The trial design A5 argues is confounded by the CYP3A4 interaction — and, being prednisolone-based, also the design in which that confound is smallest.

- **Agbetile J (2014)** Effectiveness of voriconazole in the treatment of Aspergillus fumigatus-associated asthma (EVITA3 study). *J Allergy Clin Immunol* 134:33-9. [PMID 24290286](https://pubmed.ncbi.nlm.nih.gov/24290286/)  
  EVITA3: voriconazole in A. fumigatus-associated asthma.

- **Godet C (2022)** Nebulised liposomal amphotericin-B as maintenance therapy in allergic bronchopulmonary aspergillosis: a randomised, multicentre trial. *Eur Respir J* 59. [PMID 34764182](https://pubmed.ncbi.nlm.nih.gov/34764182/)  
  NEBULAMB: nebulised liposomal amphotericin B as maintenance therapy in ABPA. High airway exposure with no systemic CYP3A4 interaction — the cleanest available test of an antifungal effect uncontaminated by a steroid drug interaction.

- **Ram B (2016)** A pilot randomized trial of nebulized amphotericin in patients with allergic bronchopulmonary aspergillosis. *J Asthma* 53:517-24. [PMID 26666774](https://pubmed.ncbi.nlm.nih.gov/26666774/)  
  Pilot randomised trial of nebulised amphotericin in ABPA.

- **Muthu V (2023)** Nebulized amphotericin B for preventing exacerbations in allergic bronchopulmonary aspergillosis: A systematic review and meta-analysis. *Pulm Pharmacol Ther* 81:102226. [PMID 37230237](https://pubmed.ncbi.nlm.nih.gov/37230237/)  
  Systematic review of nebulised amphotericin B for preventing ABPA exacerbations.

- **Agarwal R (2025)** New insights into the treatment of asthma complicated by allergic bronchopulmonary aspergillosis. *Expert Rev Respir Med* 19:967-979. [PMID 40474578](https://pubmed.ncbi.nlm.nih.gov/40474578/)  
  Recent synthesis of the treatment of asthma complicated by ABPA.

---

## ⑨ Biologics in ABPA and in type-2 airway disease

- **Voskamp AL (2015)** Clinical efficacy and immunologic effects of omalizumab in allergic bronchopulmonary aspergillosis. *J Allergy Clin Immunol Pract* 3:192-9. [PMID 25640470](https://pubmed.ncbi.nlm.nih.gov/25640470/)  
  Omalizumab in ABPA: clinical efficacy and immunologic effects. Note that studies in this class routinely report total IgE, which A6 shows moves in the opposite direction to the pharmacology.

- **Ramonell RP (2020)** Dupilumab treatment for allergic bronchopulmonary aspergillosis: A case series. *J Allergy Clin Immunol Pract* 8:742-743. [PMID 31811944](https://pubmed.ncbi.nlm.nih.gov/31811944/)  
  Dupilumab treatment for ABPA — the case series behind scenario 14.

- **Tomomatsu K (2023)** Real-world efficacy of anti-IL-5 treatment in patients with allergic bronchopulmonary aspergillosis. *Sci Rep* 13:5468. [PMID 37015988](https://pubmed.ncbi.nlm.nih.gov/37015988/)  
  Real-world efficacy of anti-IL-5 treatment in ABPA.

- **Tomomatsu K (2020)** Rapid clearance of mepolizumab-resistant bronchial mucus plugs in allergic bronchopulmonary aspergillosis with benralizumab treatment. *Allergol Int* 69:636-638. [PMID 32247541](https://pubmed.ncbi.nlm.nih.gov/32247541/)  
  Rapid clearance of MEPOLIZUMAB-RESISTANT bronchial mucus plugs in ABPA. Direct clinical evidence that the plug compartment can be refractory to one type-2 pathway while responding to another — the observation that most motivates giving PLUG its own clearance term.

- **Chen X (2024)** Efficacy of Biologics in Patients with Allergic Bronchopulmonary Aspergillosis: A Systematic Review and Meta-Analysis. *Lung* 202:367-383. [PMID 38898129](https://pubmed.ncbi.nlm.nih.gov/38898129/)  
  Systematic review and meta-analysis of biologics in ABPA.

- **Asano K (2025)** Treatment of allergic bronchopulmonary aspergillosis with biologics. *Chin Med J Pulm Crit Care Med* 3:6-11. [PMID 40226607](https://pubmed.ncbi.nlm.nih.gov/40226607/)  
  Treatment of ABPA with biologics: current position.

- **Pavord ID (2012)** Mepolizumab for severe eosinophilic asthma (DREAM): a multicentre, double-blind, placebo-controlled trial. *Lancet* 380:651-9. [PMID 22901886](https://pubmed.ncbi.nlm.nih.gov/22901886/)  
  DREAM: mepolizumab in severe eosinophilic asthma; the anti-IL-5 exposure–response anchor.

- **Ortega HG (2014)** Mepolizumab treatment in patients with severe eosinophilic asthma. *N Engl J Med* 371:1198-207. [PMID 25199059](https://pubmed.ncbi.nlm.nih.gov/25199059/)  
  MENSA: mepolizumab; the eosinophil-depletion magnitude behind `Emax_mep`.

- **Busse WW (2019)** Long-term safety and efficacy of benralizumab in patients with severe, uncontrolled asthma: 1-year results from the BORA phase 3 extension trial. *Lancet Respir Med* 7:46-59. [PMID 30416083](https://pubmed.ncbi.nlm.nih.gov/30416083/)  
  Long-term benralizumab safety and efficacy.

- **Dagher R (2022)** Novel mechanisms of action contributing to benralizumab's potent anti-eosinophilic activity. *Eur Respir J* 59. [PMID 34289975](https://pubmed.ncbi.nlm.nih.gov/34289975/)  
  Mechanisms of benralizumab's anti-eosinophilic activity including ADCC — the basis for a separate `kADCC` term rather than treating anti-IL-5Rα as merely a stronger anti-IL-5.

- **Wenzel S (2016)** Dupilumab efficacy and safety in adults with uncontrolled persistent asthma despite use of medium-to-high-dose inhaled corticosteroids plus a long-acting β2 agonist: a randomised double-blind placebo-controlled pivotal phase 2b dose-ranging trial. *Lancet* 388:31-44. [PMID 27130691](https://pubmed.ncbi.nlm.nih.gov/27130691/)  
  Dupilumab in uncontrolled persistent asthma; `Emax_dup` and the transient blood eosinophilia the model reproduces via blocked tissue egress.

- **Castro M (2020)** Dupilumab improves lung function in patients with uncontrolled, moderate-to-severe asthma. *ERJ Open Res* 6. [PMID 32010719](https://pubmed.ncbi.nlm.nih.gov/32010719/)  
  Dupilumab improves lung function; the FEV₁ magnitude the dupilumab arm is checked against.

---

## ⑩ Omalizumab pharmacology — the TMDD that makes total IgE the wrong readout

- **Lowe PJ (2009)** Relationship between omalizumab pharmacokinetics, IgE pharmacodynamics and symptoms in patients with severe persistent allergic (IgE-mediated) asthma. *Br J Clin Pharmacol* 68:61-76. [PMID 19660004](https://pubmed.ncbi.nlm.nih.gov/19660004/)  
  Lowe's mechanistic account of omalizumab PK, IgE PD and symptoms. The rapid-equilibrium binding formulation used here — totals carried as states, complex solved algebraically each step — follows this paper's structure, and it is why the model can be non-stiff without pretending binding is slow.

- **Lowe PJ (2011)** Omalizumab decreases IgE production in patients with allergic (IgE-mediated) asthma; PKPD analysis of a biomarker, total IgE. *Br J Clin Pharmacol* 72:306-20. [PMID 21392073](https://pubmed.ncbi.nlm.nih.gov/21392073/)  
  Omalizumab DECREASES IgE PRODUCTION (PK/PD analysis). This is one of three mutually indistinguishable explanations A7 offers for why the observed total-IgE rise (2–5×) is smaller than a fixed-production model predicts, and it is the one this paper argues for. The model does not adjudicate between them and says so.

- **Lowe PJ (2015)** Revision of omalizumab dosing table for dosing every 4 instead of 2 weeks for specific ranges of bodyweight and baseline IgE. *Regul Toxicol Pharmacol* 71:68-77. [PMID 25497995](https://pubmed.ncbi.nlm.nih.gov/25497995/)  
  Revision of the omalizumab dosing table for q4w dosing. Evidence that the table is a derived object rather than a fixed fact — which is what A8 exploits when it reconstructs the table's functional form from a molar flux balance.

- **Honma W (2016)** Ethnic sensitivity assessment of pharmacokinetics and pharmacodynamics of omalizumab with dosing table expansion. *Drug Metab Pharmacokinet* 31:173-84. [PMID 27238573](https://pubmed.ncbi.nlm.nih.gov/27238573/)  
  Omalizumab dosing table EXTENSION and ethnic sensitivity assessment. Directly relevant to A8: the table's boundaries are a modelling decision, and ABPA sits outside them.

- **Sorkness CA (2013)** Reassessment of omalizumab-dosing strategies and pharmacodynamics in inner-city children and adolescents. *J Allergy Clin Immunol Pract* 1:163-71. [PMID 24565455](https://pubmed.ncbi.nlm.nih.gov/24565455/)  
  Reassessment of omalizumab dosing strategies and pharmacodynamics — free-IgE targets are not reached in a substantial fraction of patients on table doses.

- **Azzano P (2021)** Determinants of omalizumab dose-related efficacy in oral immunotherapy: Evidence from a cohort of 181 patients. *J Allergy Clin Immunol* 147:233-243. [PMID 32980425](https://pubmed.ncbi.nlm.nih.gov/32980425/)  
  Determinants of omalizumab DOSE-RELATED efficacy: response tracking dose per unit IgE burden rather than dose. The clinical shape of A8's prediction.

- **Guo G (2023)** Physiologically-Based Pharmacokinetic Modeling of Omalizumab to Predict the Pharmacokinetics and Pharmacodynamics in Pediatric Patients. *Clin Pharmacol Ther* 113:724-734. [PMID 36495063](https://pubmed.ncbi.nlm.nih.gov/36495063/)  
  PBPK modelling of omalizumab; independent confirmation of the disposition parameters used here.

- **Gauvreau GM (2016)** Efficacy and safety of multiple doses of QGE031 (ligelizumab) versus omalizumab and placebo in inhibiting allergen-induced early asthmatic responses. *J Allergy Clin Immunol* 138:1051-1059. [PMID 27185571](https://pubmed.ncbi.nlm.nih.gov/27185571/)  
  Ligelizumab versus omalizumab: a higher-affinity anti-IgE achieving deeper free-IgE suppression. The cleanest available demonstration that free IgE, not total IgE, is the pharmacologically meaningful species.

- **Shields RL (1995)** Anti-IgE monoclonal antibodies that inhibit allergen-specific histamine release. *Int Arch Allergy Immunol* 107:412-3. [PMID 7542094](https://pubmed.ncbi.nlm.nih.gov/7542094/)  
  Anti-IgE monoclonal antibodies that inhibit allergen-specific histamine release — the original in vitro basis for the FcεRI occupancy term.

---

## ⑪ Bronchiectasis as the irreversible endpoint

- **Sehgal IS (2025)** Impact of Bronchiectasis Severity on Clinical Outcomes in Patients With Allergic Bronchopulmonary Aspergillosis: A Retrospective Cohort Study. *J Allergy Clin Immunol Pract* 13:1103-1109.e2. [PMID 40088971](https://pubmed.ncbi.nlm.nih.gov/40088971/)  
  Impact of bronchiectasis severity on clinical outcomes in ABPA. The justification for BRON entering FEV₁ as a fixed, non-reversible term.

- **Phadnis S (2024)** Bronchiectasis Severity Index and FACED scores in patients with allergic bronchopulmonary aspergillosis complicating asthma: do they correlate with immunological severity or high-attenuation mucus?. *J Asthma* 61:1242-1247. [PMID 38520686](https://pubmed.ncbi.nlm.nih.gov/38520686/)  
  Bronchiectasis Severity Index and FACED scores in ABPA.

- **Tiew PY (2026)** Aspergillus fumigatus Sensitization Is Associated With High-Risk Bronchiectasis. *Chest* 169:932-946. [PMID 41386457](https://pubmed.ncbi.nlm.nih.gov/41386457/)  
  A. fumigatus sensitisation is associated with high-risk bronchiectasis — the outer loop of the model's vicious-vortex cluster.

- **Agarwal R (2026)** Does ABPA Contribute to Bronchiectasis? A Structured Evaluation of Competing Hypotheses. *Clin Exp Allergy* 56:498-505. [PMID 41943086](https://pubmed.ncbi.nlm.nih.gov/41943086/)  
  A structured evaluation of whether ABPA actually causes bronchiectasis. Read this before believing the model's BRON integrator: the causal direction the model assumes is not settled, and this paper is the reason that caveat is stated in the README rather than omitted.

---

## Reference-integrity note

97 references, all resolved through
`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi`. Where a search returned
no usable record the intended citation was **dropped rather than guessed**, which is why
some standard wiring in the mechanistic map (for example the ILC2 and Treg edges in
cluster ④) carries no citation here: that wiring is textbook type-2 immunology and no
specific paper is claimed for it.

Two things this list deliberately does **not** contain, both of them load-bearing:

1. **A reference for `g_p`, the intra-plug fungal growth rate.** There is none, because
   the quantity has never been measured in a human airway. A2(e) shows that the model's
   central conclusion turns on `sign(g_p − k_out)`, so this is the most consequential
   missing measurement in the whole model. Manufacturing a citation for it would have
   concealed exactly the thing most worth knowing.
2. **A trial comparing itraconazole's steroid-sparing effect across steroid backbones.**
   A5 predicts such a trial would find a several-fold larger effect on
   methylprednisolone or inhaled budesonide than on prednisolone, from identical
   antifungal pharmacology. That comparison has not been run, so the prediction stands
   untested rather than supported.
