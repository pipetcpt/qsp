# Distal Renal Tubular Acidosis — References

Curated bibliography supporting `drta_qsp_model.dot`, `drta_mrgsolve_model.R` and `drta_shiny_app.R`.

**Provenance.** Every entry below was retrieved programmatically from NCBI E-utilities (`esearch` + `esummary`) against PubMed. The title, journal, year and first author of each record are taken verbatim from the returned `esummary` payload — nothing in this file was written from memory, and any PMID that NCBI did not return was dropped rather than reconstructed. Sixty-five distinct queries were run; the unique records they returned were de-duplicated and assigned to the sections below in NCBI relevance order.

**How to read this list.** The sections follow the structure of the model rather than the structure of a textbook chapter. Sections 1-5 establish the disease and its genetics; 6-8 the alpha-intercalated cell that the model treats as a *saturating actuator*; 9-12 the three acid sinks (ECF bicarbonate, non-bicarbonate buffer, bone) on their three timescales; 13-15 the citrate valve and calcium-phosphate crystallisation; 16-18 potassium, growth and the progression to CKD; 19-21 the pharmacology, including the prolonged-release formulation whose two-granule design the model rationalises; 22-24 differential diagnosis, measurement and method.

## 1. Distal renal tubular acidosis — reviews, classification, epidemiology

Entry points to the disease and to the classification that separates a distal (type 1) acidification defect from a proximal (type 2) reabsorptive defect and from hyperkalaemic (type 4) RTA.

- Wagner CA et al. (2023). *The pathophysiology of distal renal tubular acidosis.* Nat Rev Nephrol. [PMID 37016093](https://pubmed.ncbi.nlm.nih.gov/37016093/)
- Giglio S et al. (2021). *Distal renal tubular acidosis: a systematic approach from diagnosis to treatment.* J Nephrol. [PMID 33770395](https://pubmed.ncbi.nlm.nih.gov/33770395/)
- Vallés PG et al. (2018). *Hypokalemic Distal Renal Tubular Acidosis.* Adv Chronic Kidney Dis. [PMID 30139458](https://pubmed.ncbi.nlm.nih.gov/30139458/)
- Batlle D et al. (2012). *Genetic causes and mechanisms of distal renal tubular acidosis.* Nephrol Dial Transplant. [PMID 23114896](https://pubmed.ncbi.nlm.nih.gov/23114896/)
- Alexander RT et al. (1993). *Hereditary Distal Renal Tubular Acidosis.* —. [PMID 31600044](https://pubmed.ncbi.nlm.nih.gov/31600044/)
- Soares SBM et al. (2019). *Distal renal tubular acidosis: genetic causes and management.* World J Pediatr. [PMID 31079338](https://pubmed.ncbi.nlm.nih.gov/31079338/)
- Fuster DG et al. (2018). *Incomplete Distal Renal Tubular Acidosis and Kidney Stones.* Adv Chronic Kidney Dis. [PMID 30139463](https://pubmed.ncbi.nlm.nih.gov/30139463/)
- Alonso-Varela M et al. (2020). *Incomplete distal renal tubular acidosis in children.* Acta Paediatr. [PMID 32212394](https://pubmed.ncbi.nlm.nih.gov/32212394/)
- Boyer O et al. (2024). *Recent Developments in the Treatment of Pediatric Distal Renal Tubular Acidosis.* Paediatr Drugs. [PMID 39325135](https://pubmed.ncbi.nlm.nih.gov/39325135/)
- Clericetti CM et al. (2018). *Hyperammonemia associated with distal renal tubular acidosis or urinary tract infection: a systematic review.* Pediatr Nephrol. [PMID 29134448](https://pubmed.ncbi.nlm.nih.gov/29134448/)
- Batlle D et al. (1982). *Distal renal tubular acidosis: pathogenesis and classification.* Am J Kidney Dis. [PMID 6807085](https://pubmed.ncbi.nlm.nih.gov/6807085/)
- Mohebbi N et al. (2018). *Pathophysiology, diagnosis and treatment of inherited distal renal tubular acidosis.* J Nephrol. [PMID 28994037](https://pubmed.ncbi.nlm.nih.gov/28994037/)
- Santos F et al. (2023). *Long-term complications of primary distal renal tubular acidosis.* Pediatr Nephrol. [PMID 35543873](https://pubmed.ncbi.nlm.nih.gov/35543873/)
- Batlle D et al. (2006). *Distal renal tubular acidosis and the potassium enigma.* Semin Nephrol. [PMID 17275585](https://pubmed.ncbi.nlm.nih.gov/17275585/)
- Gómez-Conde S et al. (2021). *Hereditary distal renal tubular acidosis: Genotypic correlation, evolution to long term, and new therapeutic perspectives.* Nefrologia (Engl Ed). [PMID 36165107](https://pubmed.ncbi.nlm.nih.gov/36165107/)
- Karet FE (2002). *Inherited distal renal tubular acidosis.* J Am Soc Nephrol. [PMID 12138152](https://pubmed.ncbi.nlm.nih.gov/12138152/)
- Bouzidi H et al. (2009). *[Primary distal renal tubular acidosis].* Ann Biol Clin (Paris). [PMID 19297287](https://pubmed.ncbi.nlm.nih.gov/19297287/)
- Nicoletta JA et al. (2004). *Distal renal tubular acidosis.* Curr Opin Pediatr. [PMID 15021201](https://pubmed.ncbi.nlm.nih.gov/15021201/)
- Yenchitsomanus PT (2003). *Human anion exchanger1 mutations and distal renal tubular acidosis.* Southeast Asian J Trop Med Public Health. [PMID 15115146](https://pubmed.ncbi.nlm.nih.gov/15115146/)
- Khositseth S et al. (2012). *Tropical distal renal tubular acidosis: clinical and epidemiological studies in 78 patients.* QJM. [PMID 22919024](https://pubmed.ncbi.nlm.nih.gov/22919024/)
- Vasquez-Rios G et al. (2019). *Distal renal tubular acidosis and severe hypokalemia: a case report and review of the literature.* J Med Case Rep. [PMID 31023369](https://pubmed.ncbi.nlm.nih.gov/31023369/)
- Wrong O (1991). *Distal renal tubular acidosis: the value of urinary pH, PCO2 and NH4+ measurements.* Pediatr Nephrol. [PMID 1903268](https://pubmed.ncbi.nlm.nih.gov/1903268/)
- Ungureanu O et al. (2022). *Distal Renal Tubular Acidosis in Patients with Autoimmune Diseases-An Update on Pathogenesis, Clinical Presentation and Therapeutic Strategies.* Biomedicines. [PMID 36140232](https://pubmed.ncbi.nlm.nih.gov/36140232/)
- Watanabe T (2018). *Improving outcomes for patients with distal renal tubular acidosis: recent advances and challenges ahead.* Pediatric Health Med Ther. [PMID 30588151](https://pubmed.ncbi.nlm.nih.gov/30588151/)
- Medina E et al. (2024). *Primary Distal Renal Tubular Acidosis: Toward an Optimal Correction of Metabolic Acidosis.* Clin J Am Soc Nephrol. [PMID 38967973](https://pubmed.ncbi.nlm.nih.gov/38967973/)
- Trepiccione F et al. (2017). *New Findings on the Pathogenesis of Distal Renal Tubular Acidosis.* Kidney Dis (Basel). [PMID 29344504](https://pubmed.ncbi.nlm.nih.gov/29344504/)

## 2. Incomplete dRTA and the acid-loading tests

The model derives incomplete dRTA rather than assuming it: the actuator is already railed at a normal dietary acid load, so plasma bicarbonate is normal while the acid-load test fails and the bone and citrate fluxes are already non-zero. These are the papers that define and quantify that state, including the half-dose NH4Cl protocol and the furosemide/fludrocortisone alternative.

- Santos F et al. (2017). *Renal tubular acidosis.* Curr Opin Pediatr. [PMID 28092281](https://pubmed.ncbi.nlm.nih.gov/28092281/)
- Houillier P et al. (2025). *Autoimmune Tubulopathies.* J Am Soc Nephrol. [PMID 39786900](https://pubmed.ncbi.nlm.nih.gov/39786900/)
- Guimerà J et al. (2024). *Phytate Effects on Incomplete Distal Renal Tubular Acidosis.* J Clin Med. [PMID 39274272](https://pubmed.ncbi.nlm.nih.gov/39274272/)
- Seong EY et al. (2023). *Incomplete distal renal tubular acidosis uncovered during pregnancy: A case report.* World J Clin Cases. [PMID 37727491](https://pubmed.ncbi.nlm.nih.gov/37727491/)
- Papatsoris A et al. (2025). *Management of urinary stones by experts in stone disease (ESD 2025).* Arch Ital Urol Androl. [PMID 40583613](https://pubmed.ncbi.nlm.nih.gov/40583613/)
- Donnelly S et al. (1992). *Might distal renal tubular acidosis be a proximal tubular cell disorder?.* Am J Kidney Dis. [PMID 1553972](https://pubmed.ncbi.nlm.nih.gov/1553972/)
- Goldfarb DS (2017). *Refining Diagnostic Approaches in Nephrolithiasis: Incomplete Distal Renal Tubular Acidosis.* Clin J Am Soc Nephrol. [PMID 28775128](https://pubmed.ncbi.nlm.nih.gov/28775128/)
- Woo KT et al. (1986). *Renal tubular acidosis.* Ann Acad Med Singap. [PMID 3707033](https://pubmed.ncbi.nlm.nih.gov/3707033/)
- Buckalew VM Jr (1989). *Nephrolithiasis in renal tubular acidosis.* J Urol. [PMID 2645431](https://pubmed.ncbi.nlm.nih.gov/2645431/)
- Fernandes PC et al. (2024). *Renal tubular acidosis in hereditary transthyretin amyloidosis (ATTRv).* J Bras Nefrol. [PMID 39101566](https://pubmed.ncbi.nlm.nih.gov/39101566/)
- Assimos DG (2019). *Re: Incomplete Distal Renal Tubular Acidosis and Kidney Stones.* J Urol. [PMID 30759664](https://pubmed.ncbi.nlm.nih.gov/30759664/)
- Magni G et al. (2021). *Renal tubular acidosis (RTA) and kidney stones: Diagnosis and management.* Arch Esp Urol. [PMID 33459628](https://pubmed.ncbi.nlm.nih.gov/33459628/)
- Osther PJ et al. (1989). *Screening renal stone formers for distal renal tubular acidosis.* Br J Urol. [PMID 2752250](https://pubmed.ncbi.nlm.nih.gov/2752250/)
- Weger M et al. (1999). *Incomplete renal tubular acidosis in 'primary' osteoporosis.* Osteoporos Int. [PMID 10692983](https://pubmed.ncbi.nlm.nih.gov/10692983/)
- Sromicki J et al. (2022). *Prospective long-term evaluation of incomplete distal renal tubular acidosis in idiopathic calcium nephrolithiasis diagnosed by low-dose NH(4)CL loading - gender prevalences and impact of alkali treatment.* J Nephrol. [PMID 34973150](https://pubmed.ncbi.nlm.nih.gov/34973150/)
- Zhang J et al. (2014). *Incomplete distal renal tubular acidosis from a heterozygous mutation of the V-ATPase B1 subunit.* Am J Physiol Renal Physiol. [PMID 25164082](https://pubmed.ncbi.nlm.nih.gov/25164082/)
- Both T et al. (2015). *Prevalence of distal renal tubular acidosis in primary Sjögren's syndrome.* Rheumatology (Oxford). [PMID 25354755](https://pubmed.ncbi.nlm.nih.gov/25354755/)
- Osther PJ et al. (1989). *Distal renal tubular acidosis in recurrent renal stone formers.* Dan Med Bull. [PMID 2805826](https://pubmed.ncbi.nlm.nih.gov/2805826/)
- Jha R et al. (2011). *Clinical profile of distal renal tubular acidosis.* Saudi J Kidney Dis Transpl. [PMID 21422623](https://pubmed.ncbi.nlm.nih.gov/21422623/)
- Oduwole AO et al. (2010). *Relationship between rickets and incomplete distal renal tubular acidosis in children.* Ital J Pediatr. [PMID 20699008](https://pubmed.ncbi.nlm.nih.gov/20699008/)

## 3. Urine ammonium, the urine anion gap and the osmolal gap

A positive urine anion gap in dRTA is low urinary ammonium. In the model NH4+ excretion is the product of ammoniagenic capacity and a luminal-pH-dependent trapping efficiency, so it falls whenever the urine cannot be acidified — the hallmark is produced, not asserted.

- Morikawa MJ et al. (2025). *Acid-Base Interpretation: A Practical Approach.* Am Fam Physician. [PMID 39964926](https://pubmed.ncbi.nlm.nih.gov/39964926/)
- Greenberg KI et al. (2025). *Metabolic Acidosis.* Adv Kidney Dis Health. [PMID 40175031](https://pubmed.ncbi.nlm.nih.gov/40175031/)
- Berend K et al. (2012). *Chloride: the queen of electrolytes?.* Eur J Intern Med. [PMID 22385875](https://pubmed.ncbi.nlm.nih.gov/22385875/)
- Palmer BF et al. (2025). *Mixed Acid-Base Disturbances: Core Curriculum 2025.* Am J Kidney Dis. [PMID 40728495](https://pubmed.ncbi.nlm.nih.gov/40728495/)
- Sanghavi SF et al. (2023). *Arterial Blood Gases and Acid-Base Regulation.* Semin Respir Crit Care Med. [PMID 37369215](https://pubmed.ncbi.nlm.nih.gov/37369215/)
- Batlle D et al. (2018). *The Urine Anion Gap in Context.* Clin J Am Soc Nephrol. [PMID 29311217](https://pubmed.ncbi.nlm.nih.gov/29311217/)
- Uribarri J et al. (2022). *Beyond the Urine Anion Gap: In Support of the Direct Measurement of Urinary Ammonium.* Am J Kidney Dis. [PMID 35810828](https://pubmed.ncbi.nlm.nih.gov/35810828/)
- Uribarri J et al. (2021). *The Urine Anion Gap: Common Misconceptions.* J Am Soc Nephrol. [PMID 33769949](https://pubmed.ncbi.nlm.nih.gov/33769949/)
- Raphael KL et al. (2018). *Correlation of Urine Ammonium and Urine Osmolal Gap in Kidney Transplant Recipients.* Clin J Am Soc Nephrol. [PMID 29519951](https://pubmed.ncbi.nlm.nih.gov/29519951/)
- Cuddy LN et al. (2025). *Retrospective evaluation of acid-base analysis in dogs and cats with diabetic ketosis (2017-2021): 96 cases.* J Vet Emerg Crit Care (San Antonio). [PMID 39831448](https://pubmed.ncbi.nlm.nih.gov/39831448/)
- Fujimaru T et al. (2021). *Assessing urine ammonium concentration by urine osmolal gap in chronic kidney disease.* Nephrology (Carlton). [PMID 34288275](https://pubmed.ncbi.nlm.nih.gov/34288275/)
- Raphael KL et al. (2018). *Urine Anion Gap to Predict Urine Ammonium and Related Outcomes in Kidney Disease.* Clin J Am Soc Nephrol. [PMID 29097482](https://pubmed.ncbi.nlm.nih.gov/29097482/)
- Gueutin V et al. (2012). *[Renal physiology].* Bull Cancer. [PMID 22157516](https://pubmed.ncbi.nlm.nih.gov/22157516/)
- Adeva MM et al. (2011). *Diet-induced metabolic acidosis.* Clin Nutr. [PMID 21481501](https://pubmed.ncbi.nlm.nih.gov/21481501/)
- Tapaneya-Olarn C et al. (1999). *Comparison of urine anion gap, urine osmolal gap and modified urine osmolal gap in assessing the urine ammonium in metabolic acidosis.* J Med Assoc Thai. [PMID 10730527](https://pubmed.ncbi.nlm.nih.gov/10730527/)
- Rehman MZ et al. (2023). *Urinary Ammonium in Clinical Medicine: Direct Measurement and the Urine Anion Gap as a Surrogate Marker During Metabolic Acidosis.* Adv Kidney Dis Health. [PMID 36868734](https://pubmed.ncbi.nlm.nih.gov/36868734/)

## 4. Genetics: ATP6V1B1, ATP6V0A4, SLC4A1, FOXI1, WDR72

Five genes account for the inherited forms. In this model the genotype enters the acid-base equations through exactly two dimensionless numbers — LES (retained H+-pump Vmax, a rate defect) and LES_grad (retained maximal blood-to-urine pH gradient, a gradient defect) — with the extrarenal features carried separately.

- Kunchur MG et al. (2024). *A review of renal tubular acidosis.* J Vet Emerg Crit Care (San Antonio). [PMID 39023331](https://pubmed.ncbi.nlm.nih.gov/39023331/)
- Wang SSY et al. (2024). *Understanding renal tubular acidosis.* Br J Hosp Med (Lond). [PMID 39475030](https://pubmed.ncbi.nlm.nih.gov/39475030/)
- Peng SQ et al. (2025). *A gain-of-function mutation in ATP6V0A4 drives primary distal renal tubular alkalosis with enhanced V-ATPase activity.* J Clin Invest. [PMID 40299568](https://pubmed.ncbi.nlm.nih.gov/40299568/)
- Hammi Y et al. (2023). *Genotype-Phenotype correlation of distal renal tubular acidosis in Tunisia.* Tunis Med. [PMID 38445406](https://pubmed.ncbi.nlm.nih.gov/38445406/)
- Kurtz I (2018). *Renal Tubular Acidosis: H(+)/Base and Ammonia Transport Abnormalities and Clinical Syndromes.* Adv Chronic Kidney Dis. [PMID 30139460](https://pubmed.ncbi.nlm.nih.gov/30139460/)
- Priyadarshini S et al. (2025). *Etiology and outcomes of primary renal tubular acidosis.* Pediatr Nephrol. [PMID 40232499](https://pubmed.ncbi.nlm.nih.gov/40232499/)
- Batlle D et al. (2018). *Hyperkalemic Forms of Renal Tubular Acidosis: Clinical and Pathophysiological Aspects.* Adv Chronic Kidney Dis. [PMID 30139459](https://pubmed.ncbi.nlm.nih.gov/30139459/)
- Finsterer J (2008). *Primary periodic paralyses.* Acta Neurol Scand. [PMID 18031562](https://pubmed.ncbi.nlm.nih.gov/18031562/)
- Zhang L et al. (2021). *Familial distal renal tubular acidosis.* J Int Med Res. [PMID 33726529](https://pubmed.ncbi.nlm.nih.gov/33726529/)
- Wrong O et al. (2002). *Band 3 mutations, distal renal tubular acidosis, and Southeast Asian ovalocytosis.* Kidney Int. [PMID 12081559](https://pubmed.ncbi.nlm.nih.gov/12081559/)
- Thakare S et al. (2026). *Molecular genetics and long-term outcomes of primary distal renal tubular acidosis in Asia.* Nephrol Dial Transplant. [PMID 41134021](https://pubmed.ncbi.nlm.nih.gov/41134021/)
- Trepiccione F et al. (2021). *Distal renal tubular acidosis: ERKNet/ESPN clinical practice points.* Nephrol Dial Transplant. [PMID 33914889](https://pubmed.ncbi.nlm.nih.gov/33914889/)
- Whyte MP (2023). *Carbonic anhydrase II deficiency.* Bone. [PMID 36709914](https://pubmed.ncbi.nlm.nih.gov/36709914/)
- Besouw MTP et al. (2017). *Clinical and molecular aspects of distal renal tubular acidosis in children.* Pediatr Nephrol. [PMID 28188436](https://pubmed.ncbi.nlm.nih.gov/28188436/)
- de Parscau L (2010). *[Renal acidosis].* Arch Pediatr. [PMID 20654836](https://pubmed.ncbi.nlm.nih.gov/20654836/)
- Santos F et al. (2015). *Clinical and laboratory approaches in the diagnosis of renal tubular acidosis.* Pediatr Nephrol. [PMID 25823989](https://pubmed.ncbi.nlm.nih.gov/25823989/)
- Lopez-Garcia SC et al. (2019). *Treatment and long-term outcome in primary distal renal tubular acidosis.* Nephrol Dial Transplant. [PMID 30773598](https://pubmed.ncbi.nlm.nih.gov/30773598/)
- Escobar L et al. (2013). *Distal renal tubular acidosis: a hereditary disease with an inadequate urinary H⁺ excretion.* Nefrologia. [PMID 23640117](https://pubmed.ncbi.nlm.nih.gov/23640117/)
- Wagner CA et al. (2019). *Molecular Pathophysiology of Acid-Base Disorders.* Semin Nephrol. [PMID 31300090](https://pubmed.ncbi.nlm.nih.gov/31300090/)
- Batlle D et al. (2001). *Hereditary distal renal tubular acidosis: new understandings.* Annu Rev Med. [PMID 11160790](https://pubmed.ncbi.nlm.nih.gov/11160790/)
- Bruce LJ et al. (1998). *The association between familial distal renal tubular acidosis and mutations in the red cell anion exchanger (band 3, AE1) gene.* Biochem Cell Biol. [PMID 10353704](https://pubmed.ncbi.nlm.nih.gov/10353704/)
- Dahmani M et al. (2020). *ATP6V1B1 recurrent mutations in Algerian deaf patients associated with renal tubular acidosis.* Int J Pediatr Otorhinolaryngol. [PMID 31733597](https://pubmed.ncbi.nlm.nih.gov/31733597/)
- Bot Rachisan AL et al. (2026). *ATP6V1B1-Associated Inherited Distal Renal Tubular Acidosis in Children: Insights from a Literature Review.* Children (Basel). [PMID 41897147](https://pubmed.ncbi.nlm.nih.gov/41897147/)
- Han GH et al. (2023). *High ATP6V1B1 expression is associated with poor prognosis and platinum‑based chemotherapy resistance in epithelial ovarian cancer.* Oncol Rep. [PMID 36999629](https://pubmed.ncbi.nlm.nih.gov/36999629/)
- Mo S et al. (2025). *ATP6V1B1 regulates ovarian cancer progression and cisplatin sensitivity through the mTOR/autophagy pathway.* Mol Cell Biochem. [PMID 38735913](https://pubmed.ncbi.nlm.nih.gov/38735913/)
- Bourgeois S et al. (2024). *The B1 H + -ATPase ( Atp6v1b1 ) Subunit in Non-Type A Intercalated Cells is Required for Driving Pendrin Activity and the Renal Defense Against Alkalosis.* J Am Soc Nephrol. [PMID 37990364](https://pubmed.ncbi.nlm.nih.gov/37990364/)
- Nishie M et al. (2021). *Downregulated ATP6V1B1 expression acidifies the intracellular environment of cancer cells leading to resistance to antibody-dependent cellular cytotoxicity.* Cancer Immunol Immunother. [PMID 33000417](https://pubmed.ncbi.nlm.nih.gov/33000417/)
- Cogal AG et al. (2021). *Comprehensive Genetic Analysis Reveals Complexity of Monogenic Urinary Stone Disease.* Kidney Int Rep. [PMID 34805638](https://pubmed.ncbi.nlm.nih.gov/34805638/)
- Daenen M et al. (2025). *A novel, dominant disease mechanism of distal renal tubular acidosis with specific variants in ATP6V1B1.* Nephrol Dial Transplant. [PMID 39837581](https://pubmed.ncbi.nlm.nih.gov/39837581/)
- AitRaise I et al. (2022). *Genetic heterogeneity in GJB2, COL4A3, ATP6V1B1 and EDNRB variants detected among hearing impaired families in Morocco.* Mol Biol Rep. [PMID 35301649](https://pubmed.ncbi.nlm.nih.gov/35301649/)
- Bourgeois S et al. (2018). *Haploinsufficiency of the Mouse Atp6v1b1 Gene Leads to a Mild Acid-Base Disturbance with Implications for Kidney Stone Disease.* Cell Physiol Biochem. [PMID 29843146](https://pubmed.ncbi.nlm.nih.gov/29843146/)
- Subasioglu Uzak A et al. (2013). *ATP6V1B1 mutations in distal renal tubular acidosis and sensorineural hearing loss: clinical and genetic spectrum of five families.* Ren Fail. [PMID 23923981](https://pubmed.ncbi.nlm.nih.gov/23923981/)
- Wagner CA et al. (2004). *Renal vacuolar H+-ATPase.* Physiol Rev. [PMID 15383652](https://pubmed.ncbi.nlm.nih.gov/15383652/)
- Prasad SG et al. (2026). *ATP6V1B1-A Novel Genetic Association Between Pendred Imaging Phenotype and Renal Tubular Acidosis.* Laryngoscope. [PMID 42376910](https://pubmed.ncbi.nlm.nih.gov/42376910/)

## 5. Sensorineural hearing loss and other extrarenal features

The same V-ATPase acidifies endolymph, which is why B1/a4/FOXI1 disease is syndromic; band 3 is the erythrocyte protein, which is why some SLC4A1 variants carry haemolysis.

- Ay E et al. (2023). *Hearing Loss Related to Gene Mutations in Distal Renal Tubular Acidosis.* Audiol Neurootol. [PMID 37121229](https://pubmed.ncbi.nlm.nih.gov/37121229/)
- Goldstein A et al. (1993). *Single Large-Scale Mitochondrial DNA Deletion Syndromes.* —. [PMID 20301382](https://pubmed.ncbi.nlm.nih.gov/20301382/)
- Brown MT et al. (1993). *Progressive sensorineural hearing loss in association with distal renal tubular acidosis.* Arch Otolaryngol Head Neck Surg. [PMID 8457310](https://pubmed.ncbi.nlm.nih.gov/8457310/)
- Walker WG (1971). *Renal tubular acidosis and deafness.* Birth Defects Orig Artic Ser. [PMID 5173338](https://pubmed.ncbi.nlm.nih.gov/5173338/)
- Santos F et al. (1991). *The syndrome of renal tubular acidosis and nerve deafness. Discordant manifestations in dizygotic twin brothers.* Pediatr Nephrol. [PMID 2031843](https://pubmed.ncbi.nlm.nih.gov/2031843/)
- Sharifian M et al. (2010). *Distal renal tubular acidosis and its relationship with hearing loss in children: preliminary report.* Iran J Kidney Dis. [PMID 20622307](https://pubmed.ncbi.nlm.nih.gov/20622307/)
- Dunger DB et al. (1980). *Renal tubular acidosis and nerve deafness.* Arch Dis Child. [PMID 7387165](https://pubmed.ncbi.nlm.nih.gov/7387165/)
- Boettger T et al. (2002). *Deafness and renal tubular acidosis in mice lacking the K-Cl co-transporter Kcc4.* Nature. [PMID 11976689](https://pubmed.ncbi.nlm.nih.gov/11976689/)
- Shinjo Y et al. (2005). *Distal renal tubular acidosis associated with large vestibular aqueduct and sensorineural hearing loss.* Acta Otolaryngol. [PMID 16076719](https://pubmed.ncbi.nlm.nih.gov/16076719/)
- Karet FE (2000). *Inherited renal tubular acidosis.* Adv Nephrol Necker Hosp. [PMID 11068640](https://pubmed.ncbi.nlm.nih.gov/11068640/)
- Stoll C et al. (1996). *Siblings with congenital renal tubular acidosis and nerve deafness.* Clin Genet. [PMID 9001807](https://pubmed.ncbi.nlm.nih.gov/9001807/)
- Chen L et al. (2020). *Screening and function discussion of a hereditary renal tubular acidosis family pathogenic gene.* Cell Death Dis. [PMID 32123165](https://pubmed.ncbi.nlm.nih.gov/32123165/)

## 6. The alpha-intercalated cell: V-ATPase, kAE1, Rhcg

The actuator itself. Note that a kAE1 trafficking defect stalls the pump from the basolateral side (the cell alkalinises) whereas a B1/a4 defect reduces Vmax directly — two routes to the same saturated controller.

- Borthwick KJ et al. (2002). *Inherited disorders of the H+-ATPase.* Curr Opin Nephrol Hypertens. [PMID 12187322](https://pubmed.ncbi.nlm.nih.gov/12187322/)
- Williamson RC et al. (2008). *Glycophorin A: Band 3 aid.* Blood Cells Mol Dis. [PMID 18304844](https://pubmed.ncbi.nlm.nih.gov/18304844/)
- Vitzthum H et al. (2024). *Novel functions of the anion exchanger AE4 (SLC4A9).* Pflugers Arch. [PMID 38195948](https://pubmed.ncbi.nlm.nih.gov/38195948/)
- Purkerson JM et al. (2015). *Distinct α-intercalated cell morphology and its modification by acidosis define regions of the collecting duct.* Am J Physiol Renal Physiol. [PMID 26084929](https://pubmed.ncbi.nlm.nih.gov/26084929/)
- Shibata S (2019). *Role of Pendrin in the Pathophysiology of Aldosterone-Induced Hypertension.* Am J Hypertens. [PMID 30982848](https://pubmed.ncbi.nlm.nih.gov/30982848/)
- Miao Y et al. (2014). *Kidney α-intercalated cells and lipocalin 2: defending the urinary tract.* J Clin Invest. [PMID 24937424](https://pubmed.ncbi.nlm.nih.gov/24937424/)
- Schwartz GJ et al. (2005). *Role of hensin in mediating the adaptation of the cortical collecting duct to metabolic acidosis.* Curr Opin Nephrol Hypertens. [PMID 15931009](https://pubmed.ncbi.nlm.nih.gov/15931009/)
- Gluck SL et al. (1996). *Distal urinary acidification from Homer Smith to the present.* Kidney Int. [PMID 8743472](https://pubmed.ncbi.nlm.nih.gov/8743472/)
- Malnic G et al. (1994). *The role of the distal nephron in the regulation of acid-base equilibrium by the kidney.* Braz J Med Biol Res. [PMID 8087090](https://pubmed.ncbi.nlm.nih.gov/8087090/)
- Al-Awqati Q et al. (2003). *Terminal differentiation of epithelia from trophectoderm to the intercalated cell: the role of hensin.* J Am Soc Nephrol. [PMID 12761233](https://pubmed.ncbi.nlm.nih.gov/12761233/)
- Schwartz GJ et al. (2002). *Acid incubation reverses the polarity of intercalated cell transporters, an effect mediated by hensin.* J Clin Invest. [PMID 11781354](https://pubmed.ncbi.nlm.nih.gov/11781354/)
- Soleimani M (2002). *Na+:HCO3- cotransporters (NBC): expression and regulation in the kidney.* J Nephrol. [PMID 12027220](https://pubmed.ncbi.nlm.nih.gov/12027220/)
- Paragas N et al. (2014). *α-Intercalated cells defend the urinary system from bacterial infection.* J Clin Invest. [PMID 24937428](https://pubmed.ncbi.nlm.nih.gov/24937428/)
- Siga E et al. (1996). *Calcitonin stimulates H+ secretion in rat kidney intercalated cells.* Am J Physiol. [PMID 8997396](https://pubmed.ncbi.nlm.nih.gov/8997396/)
- Chang JC et al. (2014). *Role of the bicarbonate-responsive soluble adenylyl cyclase in pH sensing and metabolic regulation.* Front Physiol. [PMID 24575049](https://pubmed.ncbi.nlm.nih.gov/24575049/)
- Mhlana N et al. (2023). *Distal renal tubular acidosis in a patient with Hashimoto's thyroiditis: a case report.* Biochem Med (Zagreb). [PMID 37324116](https://pubmed.ncbi.nlm.nih.gov/37324116/)
- Melo Z et al. (2013). *Molecular evidence for a role for K(+)-Cl(-) cotransporters in the kidney.* Am J Physiol Renal Physiol. [PMID 24089410](https://pubmed.ncbi.nlm.nih.gov/24089410/)
- Zhou H et al. (2025). *Rhbg interaction with CA-IV and its effects on NH(3)/NH(4)(+) and CO(2) transport.* Am J Physiol Cell Physiol. [PMID 40811810](https://pubmed.ncbi.nlm.nih.gov/40811810/)
- Weiner ID et al. (1990). *Regulation of intracellular pH in the rabbit cortical collecting tubule.* J Clin Invest. [PMID 2153152](https://pubmed.ncbi.nlm.nih.gov/2153152/)
- Steinmetz PR et al. (1996). *Scales of urine acidification: apical membrane-associated particles in turtle bladder.* Kidney Int. [PMID 8743471](https://pubmed.ncbi.nlm.nih.gov/8743471/)
- Guerra Hernández NE et al. (2023). *Autosomal dominant distal renal tubular acidosis in two pediatric patients with mutations in the SLC4A1 gene. Can the maximum urinary pCO(2) test be normal?.* Nefrologia (Engl Ed). [PMID 37775346](https://pubmed.ncbi.nlm.nih.gov/37775346/)
- Weinstein AM (2000). *A mathematical model of the outer medullary collecting duct of the rat.* Am J Physiol Renal Physiol. [PMID 10894785](https://pubmed.ncbi.nlm.nih.gov/10894785/)
- Gao X et al. (2010). *Deletion of hensin/DMBT1 blocks conversion of beta- to alpha-intercalated cells and induces distal renal tubular acidosis.* Proc Natl Acad Sci U S A. [PMID 21098262](https://pubmed.ncbi.nlm.nih.gov/21098262/)
- Al-Awqati Q (2008). *2007 Homer W. Smith award: control of terminal differentiation in epithelia.* J Am Soc Nephrol. [PMID 18199795](https://pubmed.ncbi.nlm.nih.gov/18199795/)
- Su Y et al. (2015). *Physical and functional links between anion exchanger-1 and sodium pump.* J Am Soc Nephrol. [PMID 25012180](https://pubmed.ncbi.nlm.nih.gov/25012180/)
- Tsuruoka S et al. (1999). *Mechanisms of HCO(-)(3) secretion in the rabbit connecting segment.* Am J Physiol. [PMID 10516281](https://pubmed.ncbi.nlm.nih.gov/10516281/)

## 7. Pendrin and the beta-intercalated cell

Bicarbonate secretion is the wrong-way valve: it rises with alkali therapy and so contributes to the bicarbonaturia that wastes an over-rapid dose.

- Wémeau JL et al. (2017). *Pendred syndrome.* Best Pract Res Clin Endocrinol Metab. [PMID 28648509](https://pubmed.ncbi.nlm.nih.gov/28648509/)
- Xu J et al. (2022). *Identification of IQGAP1 as a SLC26A4 (Pendrin)-Binding Protein in the Kidney.* Front Mol Biosci. [PMID 35601831](https://pubmed.ncbi.nlm.nih.gov/35601831/)
- Amlal H et al. (2010). *Deletion of the anion exchanger Slc26a4 (pendrin) decreases apical Cl(-)/HCO3(-) exchanger activity and impairs bicarbonate secretion in kidney collecting duct.* Am J Physiol Cell Physiol. [PMID 20375274](https://pubmed.ncbi.nlm.nih.gov/20375274/)
- Soleimani M (2015). *The multiple roles of pendrin in the kidney.* Nephrol Dial Transplant. [PMID 25281699](https://pubmed.ncbi.nlm.nih.gov/25281699/)
- Bourgeois S et al. (2021). *Regulation of renal pendrin activity by aldosterone.* Curr Opin Nephrol Hypertens. [PMID 33186222](https://pubmed.ncbi.nlm.nih.gov/33186222/)
- West CA et al. (2015). *The chloride-bicarbonate exchanger pendrin is increased in the kidney of the pregnant rat.* Exp Physiol. [PMID 26260990](https://pubmed.ncbi.nlm.nih.gov/26260990/)
- Tamma G et al. (2022). *Functional interplay between CFTR and pendrin: physiological and pathophysiological relevance.* Front Biosci (Landmark Ed). [PMID 35227018](https://pubmed.ncbi.nlm.nih.gov/35227018/)
- Wagner CA et al. (2011). *The anion exchanger pendrin (SLC26A4) and renal acid-base homeostasis.* Cell Physiol Biochem. [PMID 22116363](https://pubmed.ncbi.nlm.nih.gov/22116363/)
- Wall SM et al. (2015). *The role of pendrin in renal physiology.* Annu Rev Physiol. [PMID 25668022](https://pubmed.ncbi.nlm.nih.gov/25668022/)
- Wall SM (2016). *The role of pendrin in blood pressure regulation.* Am J Physiol Renal Physiol. [PMID 26538443](https://pubmed.ncbi.nlm.nih.gov/26538443/)
- Royaux IE et al. (2001). *Pendrin, encoded by the Pendred syndrome gene, resides in the apical region of renal intercalated cells and mediates bicarbonate secretion.* Proc Natl Acad Sci U S A. [PMID 11274445](https://pubmed.ncbi.nlm.nih.gov/11274445/)
- Rozenfeld J et al. (2011). *Transcriptional regulation of the pendrin gene.* Cell Physiol Biochem. [PMID 22116353](https://pubmed.ncbi.nlm.nih.gov/22116353/)

## 8. Renal ammoniagenesis and its regulation

Ammoniagenic capacity had to be driven by dietary protein (glutamine supply) as well as by acidaemia: driving it from the plasma bicarbonate error alone made the arm unrecruitable, because that error is small precisely when the kidney is compensating.

- Guder WG et al. (1987). *Renal and hepatic nitrogen metabolism in systemic acid base regulation.* J Clin Chem Clin Biochem. [PMID 3320262](https://pubmed.ncbi.nlm.nih.gov/3320262/)
- Tannen RL (1983). *Ammonia and acid-base homeostasis.* Med Clin North Am. [PMID 6135827](https://pubmed.ncbi.nlm.nih.gov/6135827/)
- Karim Z et al. (2002). *Renal handling of NH4+ in relation to the control of acid-base balance by the kidney.* J Nephrol. [PMID 12027211](https://pubmed.ncbi.nlm.nih.gov/12027211/)
- Derakhshandeh-Rishehri SM et al. (2024). *Higher Renal Net Acid Excretion, but Not Higher Phosphate Excretion, during Childhood and Adolescence Associates with the Circulating Renal Tubular Injury Marker Interleukin-18 in Adulthood.* Int J Mol Sci. [PMID 38338685](https://pubmed.ncbi.nlm.nih.gov/38338685/)
- Licht JH et al. (1985). *Familiar hyperkalaemic acidosis.* Q J Med. [PMID 3885297](https://pubmed.ncbi.nlm.nih.gov/3885297/)
- Zomorodian A et al. (2025). *Acute Effect of High Fat Intake on Urinary Acidification Parameters.* Kidney Int Rep. [PMID 40677337](https://pubmed.ncbi.nlm.nih.gov/40677337/)
- Remer T (2001). *Influence of nutrition on acid-base balance--metabolic aspects.* Eur J Nutr. [PMID 11842946](https://pubmed.ncbi.nlm.nih.gov/11842946/)
- Karim Z et al. (2005). *Renal handling of NH3/NH4+: recent concepts.* Nephron Physiol. [PMID 16113588](https://pubmed.ncbi.nlm.nih.gov/16113588/)
- Khairallah P et al. (2017). *Acid Load and Phosphorus Homeostasis in CKD.* Am J Kidney Dis. [PMID 28645705](https://pubmed.ncbi.nlm.nih.gov/28645705/)
- Bignon Y et al. (2020). *Defective bicarbonate reabsorption in Kir4.2 potassium channel deficient mice impairs acid-base balance and ammonia excretion.* Kidney Int. [PMID 31870500](https://pubmed.ncbi.nlm.nih.gov/31870500/)
- Curthoys NP et al. (2014). *Proximal tubule function and response to acidosis.* Clin J Am Soc Nephrol. [PMID 23908456](https://pubmed.ncbi.nlm.nih.gov/23908456/)
- Karim Z et al. (2006). *Recent concepts concerning the renal handling of NH3/NH4+.* J Nephrol. [PMID 16736437](https://pubmed.ncbi.nlm.nih.gov/16736437/)
- Vinay P et al. (1986). *Regulation of glutamine metabolism in dog kidney in vivo.* Kidney Int. [PMID 3515016](https://pubmed.ncbi.nlm.nih.gov/3515016/)
- Tamarappoo BK et al. (1990). *Interorgan glutamine flow regulation in metabolic acidosis.* Miner Electrolyte Metab. [PMID 2283995](https://pubmed.ncbi.nlm.nih.gov/2283995/)

## 9. Whole-body acid-base buffering, the apparent bicarbonate space and respiratory compensation

Sinks 1 and 2. The apparent bicarbonate space is ~0.5 L/kg and expands towards 0.7-1.0 L/kg in severe acidosis, which is why a single alkali bolus moves plasma bicarbonate so little and why within-day swings are small.

- Rhodes PG et al. (1977). *The effects of single infusion of hypertonic sodium bicarbonate on body composition in neonates with acidosis.* J Pediatr. [PMID 323440](https://pubmed.ncbi.nlm.nih.gov/323440/)
- Cogan MG et al. (1979). *Control of proximal bicarbonate reabsorption in normal and acidotic rats.* J Clin Invest. [PMID 500804](https://pubmed.ncbi.nlm.nih.gov/500804/)
- Rothe KF et al. (1982). *Comparison of intra-and extracellular buffering of clinically used buffer substances: tris and bicarbonate.* Acta Anaesthesiol Scand. [PMID 6287790](https://pubmed.ncbi.nlm.nih.gov/6287790/)
- Milligan CL et al. (1986). *Intracellular and extracellular acid-base status and H+ exchange with the environment after exhaustive exercise in the rainbow trout.* J Exp Biol. [PMID 3091755](https://pubmed.ncbi.nlm.nih.gov/3091755/)
- Hobe H et al. (1984). *The mechanisms of acid-base and ionoregulation in the freshwater rainbow trout during environmental hyperoxia and subsequent normoxia. I. Extra- and intracellular acid-base status.* Respir Physiol. [PMID 6427870](https://pubmed.ncbi.nlm.nih.gov/6427870/)
- Rothe KF (1984). *[Changes in the extracellular pH value and its effect on the intracellular pH value of tissues].* Fortschr Med. [PMID 6423502](https://pubmed.ncbi.nlm.nih.gov/6423502/)
- Rothe KF (1986). *[Changes in total acid-base equilibrium following administration of carbonic anhydrase inhibitor. Experimental studies in the nephrectomized rat in vivo].* Anasth Intensivther Notfallmed. [PMID 3092692](https://pubmed.ncbi.nlm.nih.gov/3092692/)
- Repetto HA et al. (2006). *Apparent bicarbonate space in children.* ScientificWorldJournal. [PMID 16493519](https://pubmed.ncbi.nlm.nih.gov/16493519/)
- Adrogué HJ et al. (1976). *Correction of severe metabolic acidosis: a physiological approach.* Acta Physiol Lat Am. [PMID 28634](https://pubmed.ncbi.nlm.nih.gov/28634/)
- Zhai J et al. (2020). *[Clinical application of Excel spreadsheet with blood gas analysis in automatic judging the type of acid-base balance disorder].* Zhonghua Wei Zhong Bing Ji Jiu Yi Xue. [PMID 32912413](https://pubmed.ncbi.nlm.nih.gov/32912413/)
- Tinits P (1983). *Oxygen therapy and oxygen toxicity.* Ann Emerg Med. [PMID 6414343](https://pubmed.ncbi.nlm.nih.gov/6414343/)
- Zouboules SM et al. (2018). *Renal reactivity: acid-base compensation during incremental ascent to high altitude.* J Physiol. [PMID 30267579](https://pubmed.ncbi.nlm.nih.gov/30267579/)
- Lun CT et al. (2016). *Differences in baseline factors and survival between normocapnia, compensated respiratory acidosis and decompensated respiratory acidosis in COPD exacerbation: A pilot study.* Respirology. [PMID 26603971](https://pubmed.ncbi.nlm.nih.gov/26603971/)
- Fall PJ (2000). *A stepwise approach to acid-base disorders. Practical patient evaluation for metabolic acidosis and other conditions.* Postgrad Med. [PMID 10728149](https://pubmed.ncbi.nlm.nih.gov/10728149/)

## 10. Net endogenous acid production, dietary acid load and PRAL

The input to the whole system, and the only variable an intervention can reduce rather than buffer.

- Poupin N et al. (2012). *Impact of the diet on net endogenous acid production and acid-base balance.* Clin Nutr. [PMID 22342140](https://pubmed.ncbi.nlm.nih.gov/22342140/)
- Frassetto LA et al. (1998). *Estimation of net endogenous noncarbonic acid production in humans from diet potassium and protein contents.* Am J Clin Nutr. [PMID 9734733](https://pubmed.ncbi.nlm.nih.gov/9734733/)
- Yeung SMH et al. (2021). *Net Endogenous Acid Excretion and Kidney Allograft Outcomes.* Clin J Am Soc Nephrol. [PMID 34135022](https://pubmed.ncbi.nlm.nih.gov/34135022/)
- Gannon RH et al. (2008). *Estimates of daily net endogenous acid production in the elderly UK population: analysis of the National Diet and Nutrition Survey (NDNS) of British adults aged 65 years and over.* Br J Nutr. [PMID 18394215](https://pubmed.ncbi.nlm.nih.gov/18394215/)
- Huston HK et al. (2015). *Net endogenous acid production and mortality in NHANES III.* Nephrology (Carlton). [PMID 25395273](https://pubmed.ncbi.nlm.nih.gov/25395273/)
- Vincent-Johnson A et al. (2023). *Diet and Metabolism in CKD-Related Metabolic Acidosis.* Semin Nephrol. [PMID 37898028](https://pubmed.ncbi.nlm.nih.gov/37898028/)
- Toba K et al. (2019). *Higher estimated net endogenous acid production with lower intake of fruits and vegetables based on a dietary survey is associated with the progression of chronic kidney disease.* BMC Nephrol. [PMID 31752746](https://pubmed.ncbi.nlm.nih.gov/31752746/)
- Baek SH et al. (2014). *A low-salt diet increases the estimated net endogenous acid production in nondiabetic chronic kidney disease patients treated with angiotensin receptor blockade.* Nephron Clin Pract. [PMID 25531146](https://pubmed.ncbi.nlm.nih.gov/25531146/)
- Ströhle A et al. (2010). *Estimation of the diet-dependent net acid load in 229 worldwide historically studied hunter-gatherer societies.* Am J Clin Nutr. [PMID 20042527](https://pubmed.ncbi.nlm.nih.gov/20042527/)
- Sebastian A et al. (2002). *Estimation of the net acid load of the diet of ancestral preagricultural Homo sapiens and their hominid ancestors.* Am J Clin Nutr. [PMID 12450898](https://pubmed.ncbi.nlm.nih.gov/12450898/)
- Zhang L et al. (2009). *Diet-dependent net acid load and risk of incident hypertension in United States women.* Hypertension. [PMID 19667248](https://pubmed.ncbi.nlm.nih.gov/19667248/)
- Ströhle A et al. (2011). *Diet-dependent net endogenous acid load of vegan diets in relation to food groups and bone health-related nutrients: results from the German Vegan Study.* Ann Nutr Metab. [PMID 22142775](https://pubmed.ncbi.nlm.nih.gov/22142775/)
- Chan R et al. (2015). *Higher estimated net endogenous Acid production may be associated with increased prevalence of nonalcoholic Fatty liver disease in chinese adults in Hong Kong.* PLoS One. [PMID 25905490](https://pubmed.ncbi.nlm.nih.gov/25905490/)
- Ronco AL et al. (2023). *Dietary Acid Load and Cancer Risk: A Review of the Uruguayan Experience.* Nutrients. [PMID 37513516](https://pubmed.ncbi.nlm.nih.gov/37513516/)
- Chan RS et al. (2009). *Estimated net endogenous acid production and intake of bone health-related nutrients in Hong Kong Chinese adolescents.* Eur J Clin Nutr. [PMID 18231119](https://pubmed.ncbi.nlm.nih.gov/18231119/)
- Huang H et al. (2025). *Relationship between dietary acid-base load and non-insulin-based resistance measures in patients with chronic kidney disease.* Front Endocrinol (Lausanne). [PMID 40620797](https://pubmed.ncbi.nlm.nih.gov/40620797/)
- Chan R et al. (2015). *Association Between Estimated Net Endogenous Acid Production and Subsequent Decline in Muscle Mass Over Four Years in Ambulatory Older Chinese People in Hong Kong: A Prospective Cohort Study.* J Gerontol A Biol Sci Med Sci. [PMID 25422383](https://pubmed.ncbi.nlm.nih.gov/25422383/)
- Parmenter BH et al. (2017). *Accuracy and precision of estimation equations to predict net endogenous acid excretion using the Australian food database.* Nutr Diet. [PMID 28731602](https://pubmed.ncbi.nlm.nih.gov/28731602/)
- Storz MA et al. (2022). *Dietary Acid Load in Gluten-Free Diets: Results from a Cross-Sectional Study.* Nutrients. [PMID 35893918](https://pubmed.ncbi.nlm.nih.gov/35893918/)
- Tran TT et al. (2024). *The association of diet-dependent acid load with colorectal cancer risk: a case-control study in Korea.* Br J Nutr. [PMID 37649268](https://pubmed.ncbi.nlm.nih.gov/37649268/)
- Alferink LJM et al. (2019). *Diet-Dependent Acid Load-The Missing Link Between an Animal Protein-Rich Diet and Nonalcoholic Fatty Liver Disease?.* J Clin Endocrinol Metab. [PMID 30977830](https://pubmed.ncbi.nlm.nih.gov/30977830/)
- Storz MA et al. (2023). *Dietary Acid Load Is Not Associated with Serum Testosterone in Men: Insights from the NHANES.* Nutrients. [PMID 37447401](https://pubmed.ncbi.nlm.nih.gov/37447401/)

## 11. Bone as the third acid sink

Sink 3, and the one that carries the damage. Physicochemical dissolution plus PTH-independent osteoclast activation. Because the bone term is rectified in pH it is convex, so a spiky bicarbonate profile costs more bone than a flat profile of the same mean.

- Bushinsky DA et al. (2022). *Effects of acid on bone.* Kidney Int. [PMID 35351460](https://pubmed.ncbi.nlm.nih.gov/35351460/)
- Bushinsky DA et al. (1993). *Decreased bone carbonate content in response to metabolic, but not respiratory, acidosis.* Am J Physiol. [PMID 8238381](https://pubmed.ncbi.nlm.nih.gov/8238381/)
- Krieger NS et al. (2004). *Mechanism of acid-induced bone resorption.* Curr Opin Nephrol Hypertens. [PMID 15199293](https://pubmed.ncbi.nlm.nih.gov/15199293/)
- Bushinsky DA et al. (1987). *Mechanism of proton-induced bone calcium release: calcium carbonate-dissolution.* Am J Physiol. [PMID 2825542](https://pubmed.ncbi.nlm.nih.gov/2825542/)
- Salcedo-Betancourt JD et al. (2024). *The Effects of Acid on Calcium and Phosphate Metabolism.* Int J Mol Sci. [PMID 38396761](https://pubmed.ncbi.nlm.nih.gov/38396761/)
- Brezina B et al. (2004). *Acid loading during treatment with sevelamer hydrochloride: mechanisms and clinical implications.* Kidney Int Suppl. [PMID 15296506](https://pubmed.ncbi.nlm.nih.gov/15296506/)
- Green J (1994). *The physicochemical structure of bone: cellular and noncellular elements.* Miner Electrolyte Metab. [PMID 8202055](https://pubmed.ncbi.nlm.nih.gov/8202055/)
- DuBose TD Jr (1982). *Acid-base physiology in uremia.* Artif Organs. [PMID 7165551](https://pubmed.ncbi.nlm.nih.gov/7165551/)
- Bushinsky DA et al. (2002). *Acute acidosis-induced alteration in bone bicarbonate and phosphate.* Am J Physiol Renal Physiol. [PMID 12372785](https://pubmed.ncbi.nlm.nih.gov/12372785/)
- Bushinsky DA (1994). *Acidosis and bone.* Miner Electrolyte Metab. [PMID 8202051](https://pubmed.ncbi.nlm.nih.gov/8202051/)
- Krieger NS et al. (2003). *Cellular mechanisms of bone resorption induced by metabolic acidosis.* Semin Dial. [PMID 14629607](https://pubmed.ncbi.nlm.nih.gov/14629607/)
- Arnett TR (2010). *Acidosis, hypoxia and bone.* Arch Biochem Biophys. [PMID 20655868](https://pubmed.ncbi.nlm.nih.gov/20655868/)
- Bushinsky DA (2001). *Acid-base imbalance and the skeleton.* Eur J Nutr. [PMID 11842949](https://pubmed.ncbi.nlm.nih.gov/11842949/)
- Kato K et al. (2013). *Promotion of osteoclast differentiation and activation in spite of impeded osteoblast-lineage differentiation under acidosis: effects of acidosis on bone metabolism.* Biosci Trends. [PMID 23524891](https://pubmed.ncbi.nlm.nih.gov/23524891/)
- Krieger NS et al. (2021). *Deletion of the proton receptor OGR1 in mouse osteoclasts impairs metabolic acidosis-induced bone resorption.* Kidney Int. [PMID 33159961](https://pubmed.ncbi.nlm.nih.gov/33159961/)
- Arnett TR (2008). *Extracellular pH regulates bone cell function.* J Nutr. [PMID 18203913](https://pubmed.ncbi.nlm.nih.gov/18203913/)
- Lee K et al. (2018). *Roles of Mitogen-Activated Protein Kinases in Osteoclast Biology.* Int J Mol Sci. [PMID 30275408](https://pubmed.ncbi.nlm.nih.gov/30275408/)
- Arnett T (2003). *Regulation of bone cell function by acid-base balance.* Proc Nutr Soc. [PMID 14506899](https://pubmed.ncbi.nlm.nih.gov/14506899/)
- Kraut JA (1995). *The role of metabolic acidosis in the pathogenesis of renal osteodystrophy.* Adv Ren Replace Ther. [PMID 7614335](https://pubmed.ncbi.nlm.nih.gov/7614335/)
- Bushinsky DA et al. (2000). *The effects of acid on bone.* Curr Opin Nephrol Hypertens. [PMID 10926173](https://pubmed.ncbi.nlm.nih.gov/10926173/)
- Bushinsky DA et al. (1993). *Effects of metabolic and respiratory acidosis on bone.* Curr Opin Nephrol Hypertens. [PMID 7859021](https://pubmed.ncbi.nlm.nih.gov/7859021/)
- Mehrotra R et al. (2003). *Metabolic acidosis in maintenance dialysis patients: clinical considerations.* Kidney Int Suppl. [PMID 14870874](https://pubmed.ncbi.nlm.nih.gov/14870874/)

## 12. Bone disease in dRTA, and the effect of alkali on bone

The dissociation that names the disease: plasma bicarbonate 22.0 mmol/L and lumbar BMD z-score -1.1 at the same visit.

- Bagga A et al. (2020). *Renal Tubular Acidosis.* Indian J Pediatr. [PMID 32591997](https://pubmed.ncbi.nlm.nih.gov/32591997/)
- Stark Z et al. (2009). *Osteopetrosis.* Orphanet J Rare Dis. [PMID 19232111](https://pubmed.ncbi.nlm.nih.gov/19232111/)
- Knochel JP (1981). *Hypophosphatemia.* West J Med. [PMID 7010790](https://pubmed.ncbi.nlm.nih.gov/7010790/)
- Lewis RA et al. (1993). *Lowe Syndrome.* —. [PMID 20301653](https://pubmed.ncbi.nlm.nih.gov/20301653/)
- Alon US (1997). *Nephrocalcinosis.* Curr Opin Pediatr. [PMID 9204244](https://pubmed.ncbi.nlm.nih.gov/9204244/)
- Laing CM et al. (2006). *Renal tubular acidosis.* J Nephrol. [PMID 16736441](https://pubmed.ncbi.nlm.nih.gov/16736441/)
- Nesterova G et al. (1993). *Cystinosis.* —. [PMID 20301574](https://pubmed.ncbi.nlm.nih.gov/20301574/)
- Tonini G et al. (2004). *Hyperparathyroidism.* Minerva Pediatr. [PMID 15249924](https://pubmed.ncbi.nlm.nih.gov/15249924/)
- Topaloglu R (2024). *Extrarenal complications of cystinosis.* Pediatr Nephrol. [PMID 38127152](https://pubmed.ncbi.nlm.nih.gov/38127152/)
- Ficicioglu C (1993). *Tyrosinemia Type I.* —. [PMID 20301688](https://pubmed.ncbi.nlm.nih.gov/20301688/)
- Bertholet-Thomas A et al. (2023). *Bone mineral density and growth changes in patients with distal renal tubular acidosis after two-years treatment with a new alkalizing drug (ADV7103).* Nefrologia (Engl Ed). [PMID 36529656](https://pubmed.ncbi.nlm.nih.gov/36529656/)
- Sahay M et al. (2013). *Renal rickets-practical approach.* Indian J Endocrinol Metab. [PMID 24251212](https://pubmed.ncbi.nlm.nih.gov/24251212/)
- Domrongkitchaiporn S et al. (2001). *Bone mineral density and histology in distal renal tubular acidosis.* Kidney Int. [PMID 11231364](https://pubmed.ncbi.nlm.nih.gov/11231364/)
- BERMAN LB et al. (1959). *Renal tubular acidosis.* Clin Proc Child Hosp Dist Columbia. [PMID 13671747](https://pubmed.ncbi.nlm.nih.gov/13671747/)
- Lemann J Jr et al. (2000). *Acid and mineral balances and bone in familial proximal renal tubular acidosis.* Kidney Int. [PMID 10972690](https://pubmed.ncbi.nlm.nih.gov/10972690/)
- Pandita KK et al. (2012). *Double osteomalacia.* Clin Cases Miner Bone Metab. [PMID 23289039](https://pubmed.ncbi.nlm.nih.gov/23289039/)
- Kintzel PE (2001). *Anticancer drug-induced kidney disorders.* Drug Saf. [PMID 11219485](https://pubmed.ncbi.nlm.nih.gov/11219485/)
- Kyle LH (1969). *Glomerular (azotemic) osteodystrophy.* Annu Rev Med. [PMID 4307817](https://pubmed.ncbi.nlm.nih.gov/4307817/)
- Mukherjee S et al. (2024). *Unusual presentation of Sjogren's syndrome.* BMJ Case Rep. [PMID 38960417](https://pubmed.ncbi.nlm.nih.gov/38960417/)
- Kumar B et al. (2020). *Tenofovir-induced delayed nephro-osteo toxicity.* J R Coll Physicians Edinb. [PMID 32936106](https://pubmed.ncbi.nlm.nih.gov/32936106/)

## 13. Citrate, NaDC1 and hypocitraturia

The Tm-limited valve. Because NaDC1 saturates, a citrate bolus escapes reabsorption — which is why the citraturic endpoint wants FAST delivery while systemic alkalinisation wants SLOW delivery.

- Aruga S et al. (2004). *OKP cells express the Na-dicarboxylate cotransporter NaDC-1.* Am J Physiol Cell Physiol. [PMID 14973148](https://pubmed.ncbi.nlm.nih.gov/14973148/)
- Siener R (2018). *Dietary Treatment of Metabolic Acidosis in Chronic Kidney Disease.* Nutrients. [PMID 29677110](https://pubmed.ncbi.nlm.nih.gov/29677110/)
- Maalouf NM (2000). *Nephrolithiasis.* —. [PMID 25905296](https://pubmed.ncbi.nlm.nih.gov/25905296/)
- Rimer JD et al. (2019). *Citrate therapy for calcium phosphate stones.* Curr Opin Nephrol Hypertens. [PMID 30531474](https://pubmed.ncbi.nlm.nih.gov/30531474/)
- Daudon M et al. (2018). *Drug-Induced Kidney Stones and Crystalline Nephropathy: Pathophysiology, Prevention and Treatment.* Drugs. [PMID 29264783](https://pubmed.ncbi.nlm.nih.gov/29264783/)
- Adomako EA et al. (2023). *Urine pH and Citrate as Predictors of Calcium Phosphate Stone Formation.* Kidney360. [PMID 37307531](https://pubmed.ncbi.nlm.nih.gov/37307531/)
- Pearle MS (2001). *Prevention of nephrolithiasis.* Curr Opin Nephrol Hypertens. [PMID 11224695](https://pubmed.ncbi.nlm.nih.gov/11224695/)
- Lemann J Jr et al. (1989). *Idiopathic hypercalciuria.* J Urol. [PMID 2645429](https://pubmed.ncbi.nlm.nih.gov/2645429/)
- Zomorodian A et al. (2025). *Citrate and calcium kidney stones.* Clin Kidney J. [PMID 40978115](https://pubmed.ncbi.nlm.nih.gov/40978115/)
- Doizi S et al. (2018). *Impact of Potassium Citrate vs Citric Acid on Urinary Stone Risk in Calcium Phosphate Stone Formers.* J Urol. [PMID 30036516](https://pubmed.ncbi.nlm.nih.gov/30036516/)
- Sakhaee K et al. (2004). *Stone forming risk of calcium citrate supplementation in healthy postmenopausal women.* J Urol. [PMID 15311008](https://pubmed.ncbi.nlm.nih.gov/15311008/)
- Ritter A et al. (2026). *Long-term citrate treatment in high-risk kidney stone formers is not associated with metabolic adverse effects.* Clin Kidney J. [PMID 42147788](https://pubmed.ncbi.nlm.nih.gov/42147788/)
- Cupisti A et al. (2007). *Insulin resistance and low urinary citrate excretion in calcium stone formers.* Biomed Pharmacother. [PMID 17184967](https://pubmed.ncbi.nlm.nih.gov/17184967/)
- Liu CJ et al. (2020). *Statins significantly alter urinary stone-related urine biochemistry in calcium kidney stone patients with dyslipidemia.* Int J Urol. [PMID 32681579](https://pubmed.ncbi.nlm.nih.gov/32681579/)
- Assadi F et al. (2017). *Preventive Kidney Stones: Continue Medical Education.* Int J Prev Med. [PMID 28966756](https://pubmed.ncbi.nlm.nih.gov/28966756/)
- Ojo OA et al. (2025). *The Role of Diet in Kidney Stone Pathogenesis and Prevention.* Curr Urol Rep. [PMID 41258546](https://pubmed.ncbi.nlm.nih.gov/41258546/)
- Trinchieri A (2014). *[Urinary calculi and infection].* Urologia. [PMID 24874306](https://pubmed.ncbi.nlm.nih.gov/24874306/)
- Dissayabutra T et al. (2018). *Urinary stone risk factors in the descendants of patients with kidney stone disease.* Pediatr Nephrol. [PMID 29594505](https://pubmed.ncbi.nlm.nih.gov/29594505/)
- Wasserstein AG (1998). *Nephrolithiasis: acute management and prevention.* Dis Mon. [PMID 9656969](https://pubmed.ncbi.nlm.nih.gov/9656969/)
- Eriksson P et al. (1996). *Risk factors of calcium stone formation in patients with primary Sjögren's syndrome.* Urol Res. [PMID 8966840](https://pubmed.ncbi.nlm.nih.gov/8966840/)

## 14. Calcium excretion, and the Lemann slope used as a validation target

Lemann's dUCa/dNAE = 0.035 mmol/mEq is deliberately not a fitted structural parameter here; the model returns 0.0357 mmol/mEq across a dietary acid titration.

- Wróblewski T et al. (2011). *[Hypercalciuria].* Przegl Lek. [PMID 21751520](https://pubmed.ncbi.nlm.nih.gov/21751520/)
- Houillier P et al. (1998). *[Hypercalciuria].* Rev Prat. [PMID 9781174](https://pubmed.ncbi.nlm.nih.gov/9781174/)
- Meher D et al. (2025). *Idiopathic Hypercalciuria: A Comprehensive Review of Clinical Insights and Management Strategies.* Cureus. [PMID 40330359](https://pubmed.ncbi.nlm.nih.gov/40330359/)
- Picado C et al. (1996). *Corticosteroid-induced bone loss. Prevention and management.* Drug Saf. [PMID 8941496](https://pubmed.ncbi.nlm.nih.gov/8941496/)
- Reusz G (1998). *[Idiopathic hypercalciuria in childhood].* Orv Hetil. [PMID 9879200](https://pubmed.ncbi.nlm.nih.gov/9879200/)
- Courbebaisse M et al. (2020). *[Nephrolithiasis: From mechanisms to preventive medical treatment].* Nephrol Ther. [PMID 32122798](https://pubmed.ncbi.nlm.nih.gov/32122798/)
- Courbebaisse M et al. (2017). *[Nephrolithiasis of adult: From mechanisms to preventive medical treatment].* Rev Med Interne. [PMID 27349612](https://pubmed.ncbi.nlm.nih.gov/27349612/)
- Bergsland KJ et al. (2013). *Role of proximal tubule in the hypocalciuric response to thiazide of patients with idiopathic hypercalciuria.* Am J Physiol Renal Physiol. [PMID 23720347](https://pubmed.ncbi.nlm.nih.gov/23720347/)
- Jungers P et al. (1991). *[Idiopathic hypercalciuria. Biological studies and therapeutic applications].* Presse Med. [PMID 1835061](https://pubmed.ncbi.nlm.nih.gov/1835061/)
- Parfitt AM (1972). *The interactions of thiazide diuretics with parathyroid hormone and vitamin D. Studies in patients with hypoparathyroidism.* J Clin Invest. [PMID 4338123](https://pubmed.ncbi.nlm.nih.gov/4338123/)
- Lau K et al. (1982). *Tubular mechanism for the spontaneous hypercalciuria in laboratory rat.* J Clin Invest. [PMID 6288772](https://pubmed.ncbi.nlm.nih.gov/6288772/)
- Pearce SH (1998). *Straightening out the renal tubule: advances in the molecular basis of the inherited tubulopathies.* QJM. [PMID 9519207](https://pubmed.ncbi.nlm.nih.gov/9519207/)

## 15. Nephrocalcinosis, calcium-phosphate stones and supersaturation

dRTA is the one acidosis whose urine is alkaline, which is why it makes calcium-phosphate rather than uric-acid stones, and why alkali therapy has to beat its own effect on urine pH.

- Alexander RT et al. (2019). *Renal Tubular Acidosis.* Pediatr Clin North Am. [PMID 30454739](https://pubmed.ncbi.nlm.nih.gov/30454739/)
- Rothstein M et al. (1990). *Renal tubular acidosis.* Endocrinol Metab Clin North Am. [PMID 2081516](https://pubmed.ncbi.nlm.nih.gov/2081516/)
- Bali DS et al. (1993). *Glycogen Storage Disease Type I.* —. [PMID 20301489](https://pubmed.ncbi.nlm.nih.gov/20301489/)
- Seidowsky A et al. (2014). *[Tubular renal acidosis].* Rev Med Interne. [PMID 24070792](https://pubmed.ncbi.nlm.nih.gov/24070792/)
- Halbritter J et al. (2015). *Fourteen monogenic genes account for 15% of nephrolithiasis/nephrocalcinosis.* J Am Soc Nephrol. [PMID 25296721](https://pubmed.ncbi.nlm.nih.gov/25296721/)
- Ilzkovitz M et al. (2022). *Kidney Stones, Proteinuria and Renal Tubular Metabolic Acidosis: What Is the Link?.* Healthcare (Basel). [PMID 35627973](https://pubmed.ncbi.nlm.nih.gov/35627973/)
- Al-Beltagi M et al. (2023). *Renal calcification in children with renal tubular acidosis: What a paediatrician ‎should ‎know‎.* World J Clin Pediatr. [PMID 38178934](https://pubmed.ncbi.nlm.nih.gov/38178934/)
- Chalkia A et al. (2021). *Distal renal tubular acidosis and nephrocalcinosis as initial manifestation of primary sjögren's syndrome.* Saudi J Kidney Dis Transpl. [PMID 35532720](https://pubmed.ncbi.nlm.nih.gov/35532720/)
- Coello Torà I et al. (2021). *[Renal tubular distal acidosis: nephrocalcinosis as initial diagnosis.].* Arch Esp Urol. [PMID 33650542](https://pubmed.ncbi.nlm.nih.gov/33650542/)
- Agrawal N et al. (2022). *Secondary distal renal tubular acidosis and sclerotic metabolic bone disease in seronegative spondyloarthropathy.* BMJ Case Rep. [PMID 35292549](https://pubmed.ncbi.nlm.nih.gov/35292549/)
- Giaccari M et al. (2026). *The European distal renal tubular acidosis registry: a five-year analysis.* Nephrol Dial Transplant. [PMID 41790493](https://pubmed.ncbi.nlm.nih.gov/41790493/)
- Igarashi T (1992). *[Renal tubular acidosis].* Nihon Rinsho. [PMID 1434012](https://pubmed.ncbi.nlm.nih.gov/1434012/)
- Leventoğlu E (2024). *Distal renal tubular acidosis as a rare complication of vesicoureteral reflux in children: a case report and literature review.* CEN Case Rep. [PMID 38637460](https://pubmed.ncbi.nlm.nih.gov/38637460/)
- Bouzidi H et al. (2011). *[Inherited tubular renal acidosis].* Ann Biol Clin (Paris). [PMID 21896404](https://pubmed.ncbi.nlm.nih.gov/21896404/)
- Gennari FJ et al. (1978). *Renal tubular acidosis.* Annu Rev Med. [PMID 25607](https://pubmed.ncbi.nlm.nih.gov/25607/)
- Siener R et al. (2023). *Risk Profile of Patients with Brushite Stone Disease and the Impact of Diet.* Nutrients. [PMID 37764875](https://pubmed.ncbi.nlm.nih.gov/37764875/)
- Bargagli M et al. (2020). *Urinary Lithogenic Risk Profile in ADPKD Patients Treated with Tolvaptan.* Clin J Am Soc Nephrol. [PMID 32527945](https://pubmed.ncbi.nlm.nih.gov/32527945/)
- Siener R et al. (2016). *Effect of L-Methionine on the Risk of Phosphate Stone Formation.* Urology. [PMID 27521063](https://pubmed.ncbi.nlm.nih.gov/27521063/)
- Whitson PA et al. (1997). *Renal stone risk assessment during Space Shuttle flights.* J Urol. [PMID 9366381](https://pubmed.ncbi.nlm.nih.gov/9366381/)
- Stevenson AE et al. (2001). *Comparison of urine composition of healthy Labrador retrievers and miniature schnauzers.* Am J Vet Res. [PMID 11703024](https://pubmed.ncbi.nlm.nih.gov/11703024/)
- Whitson PA et al. (2001). *The risk of renal stone formation during and after long duration space flight.* Nephron. [PMID 11598387](https://pubmed.ncbi.nlm.nih.gov/11598387/)
- Bushinsky DA et al. (1994). *Increased urinary saturation and kidney calcium content in genetic hypercalciuric rats.* Kidney Int. [PMID 8127022](https://pubmed.ncbi.nlm.nih.gov/8127022/)
- Lewandowski S et al. (2004). *Renal response to lithogenic and anti-lithogenic supplement challenges in a stone-free population group.* J Ren Nutr. [PMID 15232796](https://pubmed.ncbi.nlm.nih.gov/15232796/)
- Whitson PA et al. (2001). *Urine volume and its effects on renal stone risk in astronauts.* Aviat Space Environ Med. [PMID 11318017](https://pubmed.ncbi.nlm.nih.gov/11318017/)
- Stuart RO 2nd et al. (1991). *Seasonal variations in urinary risk factors among patients with nephrolithiasis.* J Lithotr Stone Dis. [PMID 11536932](https://pubmed.ncbi.nlm.nih.gov/11536932/)
- Rodgers AL et al. (2024). *Correlation research demonstrates that an inflammatory diet is a risk factor for calcium oxalate renal stone formation.* Clin Nutr ESPEN. [PMID 38479930](https://pubmed.ncbi.nlm.nih.gov/38479930/)

## 16. Potassium wasting, hypokalaemia and hypokalaemic paralysis

Distal Na+ reabsorption that is not electrically matched by H+ secretion is matched by K+ secretion instead: the H+-pump lesion is a K+-wasting lesion. This is also why sodium-based alkali aggravates the hypokalaemia it is meant to treat and why dRTA alkali is potassium-based.

- Jha N et al. (2022). *Distal renal tubular acidosis and hypokalaemic periodic paralysis during pregnancy.* J Nephrol. [PMID 34748193](https://pubmed.ncbi.nlm.nih.gov/34748193/)
- Santoso DN et al. (2022). *Distal renal tubular acidosis presenting with an acute hypokalemic paralysis in an older child with severe vesicoureteral reflux and syringomyelia: a case report.* BMC Nephrol. [PMID 35836135](https://pubmed.ncbi.nlm.nih.gov/35836135/)
- Patel JK (2021). *Distal Renal Tubular Acidosis due to Primary Sjögren's Syndrome: Presents as Hypoakalemic Paralysis with Hypokalemia-Induced Nephrogenic Diabetes Insipidus.* Saudi J Kidney Dis Transpl. [PMID 35102929](https://pubmed.ncbi.nlm.nih.gov/35102929/)
- Ahlawat SK et al. (1999). *Hypokalaemic paralysis.* Postgrad Med J. [PMID 10715756](https://pubmed.ncbi.nlm.nih.gov/10715756/)
- Hamada S et al. (2023). *Renal tubular acidosis without interstitial nephritis in Sjögren's syndrome: a case report and review of the literature.* BMC Nephrol. [PMID 37582721](https://pubmed.ncbi.nlm.nih.gov/37582721/)
- Dave M et al. (2022). *Tenofovir-induced distal renal tubular acidosis: A rare cause of recurrent hypokalaemic paralysis.* J R Coll Physicians Edinb. [PMID 36146985](https://pubmed.ncbi.nlm.nih.gov/36146985/)
- Permatasari CA et al. (2022). *Hypokalemic periodic paralysis and renal tubular acidosis in a patient with hypothyroid and autoimmune disease.* Ann Med Surg (Lond). [PMID 35242331](https://pubmed.ncbi.nlm.nih.gov/35242331/)
- Goichot B (2001). *[Genetic hypokalemia].* Rev Med Interne. [PMID 11270268](https://pubmed.ncbi.nlm.nih.gov/11270268/)
- Correia M et al. (2023). *Hypokalemic paralysis due to renal tubular acidosis: uncommon initial manifestation of primary Sjögren´s syndrome.* ARP Rheumatol. [PMID 37421194](https://pubmed.ncbi.nlm.nih.gov/37421194/)
- Koul PA et al. (1993). *Sporadic distal renal tubular acidosis and periodic hypokalaemic paralysis in Kashmir.* J Intern Med. [PMID 8501417](https://pubmed.ncbi.nlm.nih.gov/8501417/)
- Jackson I et al. (2021). *Hypokalemic Periodic Paralysis Precipitated by Thyrotoxicosis and Renal Tubular Acidosis.* Case Rep Endocrinol. [PMID 34239739](https://pubmed.ncbi.nlm.nih.gov/34239739/)
- van den Wildenberg MJ et al. (2015). *Distal renal tubular acidosis with multiorgan autoimmunity: a case report.* Am J Kidney Dis. [PMID 25533600](https://pubmed.ncbi.nlm.nih.gov/25533600/)
- Varma V et al. (2025). *Hypokalemic Periodic Paralysis in a Patient With Primary Sjögren's Syndrome and Distal Renal Tubular Acidosis: A Case Report.* Clin Med Insights Case Rep. [PMID 40894112](https://pubmed.ncbi.nlm.nih.gov/40894112/)
- Kalita J et al. (2010). *Renal tubular acidosis presenting as respiratory paralysis: report of a case and review of literature.* Neurol India. [PMID 20228475](https://pubmed.ncbi.nlm.nih.gov/20228475/)
- Seeger H et al. (2017). *Complicated pregnancies in inherited distal renal tubular acidosis: importance of acid-base balance.* J Nephrol. [PMID 28005240](https://pubmed.ncbi.nlm.nih.gov/28005240/)
- Gupta R et al. (2013). *Hypokalemic periodic paralysis and distal renal tubular acidosis associated with renal morphological changes.* Indian Pediatr. [PMID 23680609](https://pubmed.ncbi.nlm.nih.gov/23680609/)
- Sedhain A et al. (2018). *Renal Tubular Acidosis and Hypokalemic Paralysis as a First Presentation of Primary Sjögren's Syndrome.* Case Rep Nephrol. [PMID 30410805](https://pubmed.ncbi.nlm.nih.gov/30410805/)
- Sandhya P (2023). *Comprehensive analysis of clinical and laboratory features of 440 published cases of Sjögren's syndrome and renal tubular acidosis.* Int J Rheum Dis. [PMID 36324184](https://pubmed.ncbi.nlm.nih.gov/36324184/)
- Gunaratne W et al. (2020). *A case series of distal renal tubular acidosis, Southeast Asian ovalocytosis and metabolic bone disease.* BMC Nephrol. [PMID 32758154](https://pubmed.ncbi.nlm.nih.gov/32758154/)
- Li J (2022). *Hypokalemic Periodic Paralysis Secondary to Medullary Sponge Kidney Complicated With Renal Tubular Acidosis.* Cureus. [PMID 36238424](https://pubmed.ncbi.nlm.nih.gov/36238424/)

## 17. Growth failure and the GH/IGF-1 axis in metabolic acidosis

Acidosis produces GH resistance; the model routes height velocity through an IGF-1 multiplier that is a Hill function of the bicarbonate deficit.

- Donckerwolcke R et al. (1989). *Growth failure in children with renal tubular acidosis.* Semin Nephrol. [PMID 2662306](https://pubmed.ncbi.nlm.nih.gov/2662306/)
- Guizar JM et al. (1996). *Renal tubular acidosis in children with vesicoureteral reflux.* J Urol. [PMID 8648800](https://pubmed.ncbi.nlm.nih.gov/8648800/)
- M-Osman MA et al. (2023). *Pattern of hereditary renal tubular disorders in Egyptian children.* Turk J Pediatr. [PMID 37661676](https://pubmed.ncbi.nlm.nih.gov/37661676/)
- Laing CM et al. (2005). *Renal tubular acidosis: developments in our understanding of the molecular basis.* Int J Biochem Cell Biol. [PMID 15778079](https://pubmed.ncbi.nlm.nih.gov/15778079/)
- Liew YP et al. (2017). *Type 3 renal tubular acidosis associated with growth hormone deficiency.* J Pediatr Endocrinol Metab. [PMID 28888090](https://pubmed.ncbi.nlm.nih.gov/28888090/)
- Sharma AP et al. (2007). *Incomplete distal renal tubular acidosis affects growth in children.* Nephrol Dial Transplant. [PMID 17556420](https://pubmed.ncbi.nlm.nih.gov/17556420/)
- McSherry E (1978). *Acidosis and growth in nonuremic renal disease.* Kidney Int. [PMID 366229](https://pubmed.ncbi.nlm.nih.gov/366229/)
- Cachat F et al. (1993). *[Renal tubular acidosis in children].* Pediatrie. [PMID 8397383](https://pubmed.ncbi.nlm.nih.gov/8397383/)
- Rout P et al. (2026). *Hyperphosphatemia.* —. [PMID 31869067](https://pubmed.ncbi.nlm.nih.gov/31869067/)
- Murphy JL et al. (1990). *Trimethoprim/sulfamethoxazole-induced renal tubular acidosis.* Child Nephrol Urol. [PMID 2354467](https://pubmed.ncbi.nlm.nih.gov/2354467/)
- Boyer O et al. (2022). *Improved growth of a child with primary distal renal tubular acidosis after switching from a conventional alkalizing treatment to a new prolonged-release formulation containing potassium citrate and potassium bicarbonate: lessons for the clinical nephrologist.* J Nephrol. [PMID 35357683](https://pubmed.ncbi.nlm.nih.gov/35357683/)
- Dawman L et al. (2022). *Phenotype and Genotype Profile of Children with Primary Distal Renal Tubular Acidosis: A 10-Year Experience from a North Indian Teaching Institute.* J Pediatr Genet. [PMID 35990030](https://pubmed.ncbi.nlm.nih.gov/35990030/)
- Santos F et al. (1986). *Renal tubular acidosis in children. Diagnosis, treatment and prognosis.* Am J Nephrol. [PMID 3777038](https://pubmed.ncbi.nlm.nih.gov/3777038/)
- Pereira PC et al. (2009). *Molecular pathophysiology of renal tubular acidosis.* Curr Genomics. [PMID 19721811](https://pubmed.ncbi.nlm.nih.gov/19721811/)
- Seikaly M et al. (1996). *Nephrocalcinosis is associated with renal tubular acidosis in children with X-linked hypophosphatemia.* Pediatrics. [PMID 8545232](https://pubmed.ncbi.nlm.nih.gov/8545232/)
- Sharma AP et al. (2009). *Bicarbonate therapy improves growth in children with incomplete distal renal tubular acidosis.* Pediatr Nephrol. [PMID 19347368](https://pubmed.ncbi.nlm.nih.gov/19347368/)

## 18. Progression to chronic kidney disease

The vicious loop that turns a tubular disease into CKD: crystal deposition causes fibrosis, fibrosis costs nephrons, and fewer nephrons mean less acid excretion capacity, which deepens the acidosis.

- Adeva-Andany MM et al. (2014). *Sodium bicarbonate therapy in patients with metabolic acidosis.* ScientificWorldJournal. [PMID 25405229](https://pubmed.ncbi.nlm.nih.gov/25405229/)
- Abbasi M et al. (2025). *Management of Kidney Disease with Sickle Cell Disease.* J Am Soc Nephrol. [PMID 40569673](https://pubmed.ncbi.nlm.nih.gov/40569673/)
- Oliveira JL et al. (2010). *Lithium nephrotoxicity.* Rev Assoc Med Bras (1992). [PMID 21152836](https://pubmed.ncbi.nlm.nih.gov/21152836/)
- Becker A et al. (2006). *Obstructive uropathy.* Early Hum Dev. [PMID 16377104](https://pubmed.ncbi.nlm.nih.gov/16377104/)
- Bovée DM et al. (2020). *Salt-sensitive hypertension in chronic kidney disease: distal tubular mechanisms.* Am J Physiol Renal Physiol. [PMID 32985236](https://pubmed.ncbi.nlm.nih.gov/32985236/)
- Halbritter J et al. (2009). *[Interstitial nephritis].* Internist (Berl). [PMID 19690821](https://pubmed.ncbi.nlm.nih.gov/19690821/)
- Meola M et al. (2016). *Clinical Scenarios in Chronic Kidney Disease: Cystic Renal Diseases.* Contrib Nephrol. [PMID 27169740](https://pubmed.ncbi.nlm.nih.gov/27169740/)
- Monet-Didailler C et al. (2021). *[Nephrocalcinosis in children].* Nephrol Ther. [PMID 33461896](https://pubmed.ncbi.nlm.nih.gov/33461896/)
- Baddam S et al. (2026). *Sickle Cell Nephropathy.* —. [PMID 30252273](https://pubmed.ncbi.nlm.nih.gov/30252273/)
- Gaggl M et al. (2014). *Effect of oral alkali supplementation on progression of chronic kidney disease.* Curr Hypertens Rev. [PMID 25549843](https://pubmed.ncbi.nlm.nih.gov/25549843/)
- Jiang Y et al. (2025). *Identification of novel pathogenic mutations in ATP6V0A4 associated with distal renal tubular acidosis and analysis of wild-type expression in glomerular disease.* Ren Fail. [PMID 40775604](https://pubmed.ncbi.nlm.nih.gov/40775604/)
- Adamczak M et al. (2018). *Diagnosis and Treatment of Metabolic Acidosis in Patients with Chronic Kidney Disease - Position Statement of the Working Group of the Polish Society of Nephrology.* Kidney Blood Press Res. [PMID 29895022](https://pubmed.ncbi.nlm.nih.gov/29895022/)
- Klaboch J et al. (2012). *[End stage of chronic kidney disease and metabolic acidosis].* Vnitr Lek. [PMID 23067161](https://pubmed.ncbi.nlm.nih.gov/23067161/)
- Bullen AL et al. (2023). *Markers of Kidney Tubule Dysfunction and Major Adverse Kidney Events.* Nephron. [PMID 37524063](https://pubmed.ncbi.nlm.nih.gov/37524063/)
- Krishnamurthy S (2026). *Refractory Rickets: Evaluation and Management.* Indian J Pediatr. [PMID 41741919](https://pubmed.ncbi.nlm.nih.gov/41741919/)
- Havlín J et al. (2016). *[Metabolic acidosis in chronic kidney disease].* Vnitr Lek. [PMID 28124929](https://pubmed.ncbi.nlm.nih.gov/28124929/)
- Song A et al. (2021). *Mechanism and application of metformin in kidney diseases: An update.* Biomed Pharmacother. [PMID 33714781](https://pubmed.ncbi.nlm.nih.gov/33714781/)
- Wile D (2012). *Diuretics: a review.* Ann Clin Biochem. [PMID 22783025](https://pubmed.ncbi.nlm.nih.gov/22783025/)
- Chao MV et al. (2026). *Severe metabolic acidosis with renal tubular acidosis features attributed to topical brinzolamide in a patient with stage 3 chronic kidney disease: a case report.* BMC Nephrol. [PMID 42410525](https://pubmed.ncbi.nlm.nih.gov/42410525/)
- Kim GH et al. (2025). *Urine pH and urine ammonium as biomarkers in kidney disease.* Kidney Blood Press Res. [PMID 40759108](https://pubmed.ncbi.nlm.nih.gov/40759108/)
- Bindal T et al. (2026). *Spectrum of kidney disease in pediatric sarcoidosis.* Pediatr Nephrol. [PMID 40926164](https://pubmed.ncbi.nlm.nih.gov/40926164/)
- Bürki R et al. (2015). *Impaired expression of key molecules of ammoniagenesis underlies renal acidosis in a rat model of chronic kidney disease.* Nephrol Dial Transplant. [PMID 25523450](https://pubmed.ncbi.nlm.nih.gov/25523450/)
- Schnaper HW (2017). *The Tubulointerstitial Pathophysiology of Progressive Kidney Disease.* Adv Chronic Kidney Dis. [PMID 28284376](https://pubmed.ncbi.nlm.nih.gov/28284376/)
- Atmis B et al. (2020). *Evaluation of phenotypic and genotypic features of children with distal kidney tubular acidosis.* Pediatr Nephrol. [PMID 32613277](https://pubmed.ncbi.nlm.nih.gov/32613277/)
- Erejuwa OO et al. (2020). *Effects of honey supplementation on renal dysfunction and metabolic acidosis in rats with high-fat diet-induced chronic kidney disease.* J Basic Clin Physiol Pharmacol. [PMID 32396139](https://pubmed.ncbi.nlm.nih.gov/32396139/)
- Mariano F et al. (2021). *Metformin, chronic nephropathy and lactic acidosis: a multi-faceted issue for the nephrologist.* J Nephrol. [PMID 33373028](https://pubmed.ncbi.nlm.nih.gov/33373028/)

## 19. ADV7103 / Sibnayal — the prolonged-release formulation and its trials

The primary calibration anchors. Note in particular that the citrate granules release over 2-3 h while the bicarbonate granules release over 10-12 h: the formulation implements exactly the opposing kinetics the model derives.

- Bertholet-Thomas A et al. (2021). *Efficacy and safety of an innovative prolonged-release combination drug in patients with distal renal tubular acidosis: an open-label comparative trial versus standard of care treatments.* Pediatr Nephrol. [PMID 32712761](https://pubmed.ncbi.nlm.nih.gov/32712761/)
- Acquadro M et al. (2022). *Lived experiences of patients with distal renal tubular acidosis treated with ADV7103 and of their caregivers: a qualitative study.* Orphanet J Rare Dis. [PMID 35346296](https://pubmed.ncbi.nlm.nih.gov/35346296/)
- Guittet C et al. (2020). *Innovative prolonged-release oral alkalising formulation allowing sustained urine pH increase with twice daily administration: randomised trial in healthy adults.* Sci Rep. [PMID 32811843](https://pubmed.ncbi.nlm.nih.gov/32811843/)
- Bertholet-Thomas A et al. (2021). *Safety, efficacy, and acceptability of ADV7103 during 24 months of treatment: an open-label study in pediatric and adult patients with distal renal tubular acidosis.* Pediatr Nephrol. [PMID 33635379](https://pubmed.ncbi.nlm.nih.gov/33635379/)
- Bertholet-Thomas A et al. (2025). *6-year treatment follow-up with an extended-release alkaline formulation (Sibnayal(®)) in primary distal renal tubular acidosis.* Orphanet J Rare Dis. [PMID 40804680](https://pubmed.ncbi.nlm.nih.gov/40804680/)

## 20. Conventional alkali therapy and its pharmacokinetics

Immediate-release alkali delivers base faster than endogenous acid is produced, so part of every dose crosses the proximal reabsorptive threshold and leaves in the urine.

- Palmer BF et al. (2021). *Renal Tubular Acidosis and Management Strategies: A Narrative Review.* Adv Ther. [PMID 33367987](https://pubmed.ncbi.nlm.nih.gov/33367987/)
- Morris RC Jr (1981). *Renal tubular acidosis.* N Engl J Med. [PMID 7453756](https://pubmed.ncbi.nlm.nih.gov/7453756/)
- Kitterer D et al. (2015). *Drug-induced acid-base disorders.* Pediatr Nephrol. [PMID 25370778](https://pubmed.ncbi.nlm.nih.gov/25370778/)
- Golembiewska E et al. (2012). *Renal tubular acidosis--underrated problem?.* Acta Biochim Pol. [PMID 22693689](https://pubmed.ncbi.nlm.nih.gov/22693689/)
- Tariq H et al. (2022). *Metabolic acidosis post kidney transplantation.* Front Physiol. [PMID 36082221](https://pubmed.ncbi.nlm.nih.gov/36082221/)
- Uduman J et al. (2018). *Pseudo-Renal Tubular Acidosis: Conditions Mimicking Renal Tubular Acidosis.* Adv Chronic Kidney Dis. [PMID 30139462](https://pubmed.ncbi.nlm.nih.gov/30139462/)
- Tan HL et al. (2024). *Treatment of paediatric renal tubular acidosis with a prolonged-release alkali supplementation.* Pediatr Nephrol. [PMID 38771324](https://pubmed.ncbi.nlm.nih.gov/38771324/)
- MacMahon T et al. (2023). *Zonisamide-induced distal renal tubular acidosis and critical hypokalaemia.* BMJ Case Rep. [PMID 37041041](https://pubmed.ncbi.nlm.nih.gov/37041041/)
- Morris RC Jr et al. (2002). *Alkali therapy in renal tubular acidosis: who needs it?.* J Am Soc Nephrol. [PMID 12138154](https://pubmed.ncbi.nlm.nih.gov/12138154/)
- Bharani A et al. (2018). *Distal renal tubular acidosis secondary to vesico-ureteric reflux: A case report with review of literature.* Saudi J Kidney Dis Transpl. [PMID 30381529](https://pubmed.ncbi.nlm.nih.gov/30381529/)
- Singhania P et al. (2025). *Distal Renal Tubular Acidosis: A Conundrum of Short Stature, Failure to Thrive, Rickets, and Nephrocalcinosis.* Cureus. [PMID 40978892](https://pubmed.ncbi.nlm.nih.gov/40978892/)
- Droste E et al. (1975). *Diagnosis and therapy of renal tubular acidosis in infancy.* Z Kinderheilkd. [PMID 238344](https://pubmed.ncbi.nlm.nih.gov/238344/)
- Genova R et al. (1979). *[Distal tubular acidosis].* Minerva Med. [PMID 40165](https://pubmed.ncbi.nlm.nih.gov/40165/)
- Watanabe T (2005). *Proximal renal tubular dysfunction in primary distal renal tubular acidosis.* Pediatr Nephrol. [PMID 15549407](https://pubmed.ncbi.nlm.nih.gov/15549407/)
- Louis-Jean S et al. (2020). *Distal Renal Tubular Acidosis in Sjögren's Syndrome: A Case Report.* Cureus. [PMID 33083163](https://pubmed.ncbi.nlm.nih.gov/33083163/)
- Harrington TM et al. (1983). *Renal tubular acidosis. A new look at treatment of musculoskeletal and renal disease.* Mayo Clin Proc. [PMID 6222224](https://pubmed.ncbi.nlm.nih.gov/6222224/)
- Igarashi T et al. (2002). *Unraveling the molecular pathogenesis of isolated proximal renal tubular acidosis.* J Am Soc Nephrol. [PMID 12138151](https://pubmed.ncbi.nlm.nih.gov/12138151/)
- Rodriguez-Soriano J et al. (1982). *Natural history of primary distal renal tubular acidosis treated since infancy.* J Pediatr. [PMID 7131138](https://pubmed.ncbi.nlm.nih.gov/7131138/)
- Ferron GM et al. (2003). *Oral bioavailability of pantoprazole suspended in sodium bicarbonate solution.* Am J Health Syst Pharm. [PMID 12901033](https://pubmed.ncbi.nlm.nih.gov/12901033/)
- Yu L et al. (2023). *Pharmacokinetics and Pharmacodynamics of Lansoprazole/Sodium Bicarbonate Immediate-release Capsules in Healthy Chinese Subjects: An Open, Randomized, Controlled, Crossover, Single-, and Multiple-dose Trial.* Clin Pharmacol Drug Dev. [PMID 37165834](https://pubmed.ncbi.nlm.nih.gov/37165834/)
- Jiang FL et al. (2024). *Effects of Enteric-Coated Formulation of Sodium Bicarbonate on Bicarbonate Absorption and Gastrointestinal Discomfort.* Nutrients. [PMID 38474872](https://pubmed.ncbi.nlm.nih.gov/38474872/)
- Sharma VK et al. (2000). *Oral pharmacokinetics of omeprazole and lansoprazole after single and repeated doses as intact capsules or as suspensions in sodium bicarbonate.* Aliment Pharmacol Ther. [PMID 10886044](https://pubmed.ncbi.nlm.nih.gov/10886044/)
- Di Girolamo G et al. (2007). *Relative bioavailability of new formulation of paracetamol effervescent powder containing sodium bicarbonate versus paracetamol tablets: a comparative pharmacokinetic study in fed subjects.* Expert Opin Pharmacother. [PMID 17931082](https://pubmed.ncbi.nlm.nih.gov/17931082/)
- Dubray C et al. (2021). *From the pharmaceutical to the clinical: the case for effervescent paracetamol in pain management. A narrative review.* Curr Med Res Opin. [PMID 33819115](https://pubmed.ncbi.nlm.nih.gov/33819115/)

## 21. Acquired and drug-induced dRTA

Sjogren syndrome and other autoimmune tubulointerstitial disease, amphotericin B (a gradient defect by H+ back-leak), ifosfamide, lithium, topiramate and acetazolamide, tenofovir.

- François H et al. (2016). *Renal involvement in primary Sjögren syndrome.* Nat Rev Nephrol. [PMID 26568188](https://pubmed.ncbi.nlm.nih.gov/26568188/)
- Aiyegbusi O et al. (2021). *Renal Disease in Primary Sjögren's Syndrome.* Rheumatol Ther. [PMID 33367966](https://pubmed.ncbi.nlm.nih.gov/33367966/)
- Ramos-Casals M et al. (2015). *Characterization of systemic disease in primary Sjögren's syndrome: EULAR-SS Task Force recommendations for articular, cutaneous, pulmonary and renal involvements.* Rheumatology (Oxford). [PMID 26231345](https://pubmed.ncbi.nlm.nih.gov/26231345/)
- Ho K et al. (2019). *Renal tubular acidosis as the initial presentation of Sjögren's syndrome.* BMJ Case Rep. [PMID 31413059](https://pubmed.ncbi.nlm.nih.gov/31413059/)
- Zhao J et al. (2024). *Primary Sjögren's syndrome complicated by renal tubular acidosis and acute bilateral uveitis: a case report and literature review.* J Int Med Res. [PMID 39568254](https://pubmed.ncbi.nlm.nih.gov/39568254/)
- Zhang Y et al. (2023). *Renal tubular acidosis and associated factors in patients with primary Sjögren's syndrome: a registry-based study.* Clin Rheumatol. [PMID 36383239](https://pubmed.ncbi.nlm.nih.gov/36383239/)
- Furqan S et al. (2021). *Osteoporosis Complicating Renal Tubular Acidosis in Association With Sjogren's Syndrome.* Cureus. [PMID 34725619](https://pubmed.ncbi.nlm.nih.gov/34725619/)
- Ram R et al. (2014). *Renal tubular acidosis in Sjögren's syndrome: a case series.* Am J Nephrol. [PMID 25171149](https://pubmed.ncbi.nlm.nih.gov/25171149/)
- Barday Z et al. (2023). *Primary Sjögren's syndrome with renal tubular acidosis and central pontine myelinolysis: An unusual triad.* Clin Nephrol Case Stud. [PMID 37181588](https://pubmed.ncbi.nlm.nih.gov/37181588/)
- François H et al. (2020). *[Renal involvement in Sjögren's syndrome].* Nephrol Ther. [PMID 33208269](https://pubmed.ncbi.nlm.nih.gov/33208269/)
- Sandhya P et al. (2014). *Sjögren's, Renal Tubular Acidosis And Osteomalacia - An Asian Indian Series.* Open Rheumatol J. [PMID 25584094](https://pubmed.ncbi.nlm.nih.gov/25584094/)
- Jung SW et al. (2017). *Renal Tubular Acidosis in Patients with Primary Sjögren's Syndrome.* Electrolyte Blood Press. [PMID 29042903](https://pubmed.ncbi.nlm.nih.gov/29042903/)
- Gao YL et al. (2023). *Severe Hypokalemia Complicated by Acute Myopathy: Initial Manifestation of Primary Sjögren's Syndrome-Associated Renal Tubular Acidosis.* Am J Case Rep. [PMID 37481699](https://pubmed.ncbi.nlm.nih.gov/37481699/)
- Khan N et al. (2025). *Sjögren's Syndrome With Distal Renal Tubular Acidosis and Hypokalemic Myopathy in Pregnancy: A Rare Case.* Clin Case Rep. [PMID 41262291](https://pubmed.ncbi.nlm.nih.gov/41262291/)
- Li M et al. (2024). *Tofacitinib for Sjögren syndrome with renal tubular acidosis and psoriasis.* Int J Rheum Dis. [PMID 37605824](https://pubmed.ncbi.nlm.nih.gov/37605824/)
- Ramponi G et al. (2020). *Biomarkers and Diagnostic Testing for Renal Disease in Sjogren's Syndrome.* Front Immunol. [PMID 33042142](https://pubmed.ncbi.nlm.nih.gov/33042142/)
- Chow KL et al. (2022). *Severe distal renal tubular acidosis secondary to primary Sjögren syndrome: response to rituximab.* Intern Med J. [PMID 35187822](https://pubmed.ncbi.nlm.nih.gov/35187822/)
- Pessler F et al. (2006). *The spectrum of renal tubular acidosis in paediatric Sjögren syndrome.* Rheumatology (Oxford). [PMID 16159947](https://pubmed.ncbi.nlm.nih.gov/16159947/)
- Kobayashi T et al. (2006). *Fanconi's syndrome and distal (type 1) renal tubular acidosis in a patient with primary Sjögren's syndrome with monoclonal gammopathy of undetermined significance.* Clin Nephrol. [PMID 16792139](https://pubmed.ncbi.nlm.nih.gov/16792139/)
- Narayan R et al. (2018). *Distal renal tubular acidosis in Sjögren's syndrome.* Saudi J Kidney Dis Transpl. [PMID 29657223](https://pubmed.ncbi.nlm.nih.gov/29657223/)
- Fulop M et al. (2004). *Renal tubular acidosis, Sjögren syndrome, and bone disease.* Arch Intern Med. [PMID 15111378](https://pubmed.ncbi.nlm.nih.gov/15111378/)
- Tanikawa K et al. (1985). *[Renal tubular acidosis and Sjogren's syndrome].* Nihon Rinsho. [PMID 3912554](https://pubmed.ncbi.nlm.nih.gov/3912554/)
- Lim AK et al. (2013). *Distal renal tubular acidosis associated with Sjogren syndrome.* Intern Med J. [PMID 24330363](https://pubmed.ncbi.nlm.nih.gov/24330363/)
- Rajan R et al. (2021). *Beyond sicca symptoms: Osteomalacia secondary to renal tubular acidosis in Sjogren syndrome.* Joint Bone Spine. [PMID 32952003](https://pubmed.ncbi.nlm.nih.gov/32952003/)

## 22. Differential diagnosis: proximal RTA, Fanconi syndrome, type 4 RTA

Fractional excretion of bicarbonate under 3% points distal, over 15% proximal; hyperkalaemia points to type 4.

- Sakai T et al. (1981). *A case of Fanconi syndrome with type 1 renal tubular acidosis.* Jpn Circ J. [PMID 7299995](https://pubmed.ncbi.nlm.nih.gov/7299995/)
- Izumotani T et al. (1993). *An adult case of Fanconi syndrome due to a mixture of Chinese crude drugs.* Nephron. [PMID 8413772](https://pubmed.ncbi.nlm.nih.gov/8413772/)
- Johnson AJ et al. (2021). *Type 4 Hyperkalemic Renal Tubular Acidosis After Coronary Artery Bypass Grafting.* J Cardiothorac Vasc Anesth. [PMID 32888807](https://pubmed.ncbi.nlm.nih.gov/32888807/)
- Üsküdar Cansu D et al. (2020). *Hyperkalemia in type 4 renal tubular acidosis associated with systemic lupus erythematosus.* Rheumatol Int. [PMID 32166438](https://pubmed.ncbi.nlm.nih.gov/32166438/)
- Rangel EB et al. (2006). *Severe hyperkalemic type 4 renal tubular acidosis after kidney transplantation: a case report.* Transplant Proc. [PMID 17112912](https://pubmed.ncbi.nlm.nih.gov/17112912/)
- Jakes AD et al. (2016). *Renal tubular acidosis type 4 in pregnancy.* BMJ Case Rep. [PMID 26989116](https://pubmed.ncbi.nlm.nih.gov/26989116/)
- Marino CL et al. (2024). *Pseudohypoaldosteronism and acquired renal aldosterone resistance with hyperkalemic type IV renal tubular acidosis in 2 cats.* J Vet Intern Med. [PMID 38695414](https://pubmed.ncbi.nlm.nih.gov/38695414/)
- Warnock DG (1999). *Hypertension.* Semin Nephrol. [PMID 10435675](https://pubmed.ncbi.nlm.nih.gov/10435675/)
- Borges KS et al. (2024). *Non-canonical Wnt signaling triggered by WNT2B drives adrenal aldosterone production.* bioRxiv. [PMID 39229119](https://pubmed.ncbi.nlm.nih.gov/39229119/)
- Alon U et al. (1984). *Renal tubular acidosis type 4 in neonatal unilateral kidney diseases.* J Pediatr. [PMID 6726516](https://pubmed.ncbi.nlm.nih.gov/6726516/)
- Rodríguez-Soriano J (2000). *New insights into the pathogenesis of renal tubular acidosis--from functional to molecular studies.* Pediatr Nephrol. [PMID 11045400](https://pubmed.ncbi.nlm.nih.gov/11045400/)
- Nahum H et al. (1986). *Pseudohypoaldosteronism type II: proximal renal tubular acidosis and dDAVP-sensitive renal hyperkalemia.* Am J Nephrol. [PMID 3777034](https://pubmed.ncbi.nlm.nih.gov/3777034/)
- Assadi F et al. (2006). *Hyperkalemic distal renal tubular acidosis associated with Rett syndrome.* Pediatr Nephrol. [PMID 16511686](https://pubmed.ncbi.nlm.nih.gov/16511686/)
- Sethupathi V et al. (2008). *Congenital hypoaldosteronism.* Indian Pediatr. [PMID 18723916](https://pubmed.ncbi.nlm.nih.gov/18723916/)
- Sebastian A et al. (1982). *Disorders of distal nephron function.* Am J Med. [PMID 6277192](https://pubmed.ncbi.nlm.nih.gov/6277192/)
- Schlueter W et al. (1992). *On the mechanism of impaired distal acidification in hyperkalemic renal tubular acidosis: evaluation with amiloride and bumetanide.* J Am Soc Nephrol. [PMID 1450372](https://pubmed.ncbi.nlm.nih.gov/1450372/)

## 23. Measurement, monitoring and assay caveats

Surveillance is deliberately taken on a pre-dose sample, which matters here because the trough and the mean respond differently to a change of schedule.

- Connelly JT et al. (2021). *Field evaluation of a prototype tuberculosis lipoarabinomannan lateral flow assay on HIV-positive and HIV-negative patients.* PLoS One. [PMID 34310609](https://pubmed.ncbi.nlm.nih.gov/34310609/)
- Valverde MG et al. (2025). *Bringing the lab home: Evaluating the clinical accuracy of five urinary pH devices for stone prevention.* Urol Ann. [PMID 41229577](https://pubmed.ncbi.nlm.nih.gov/41229577/)
- Grases F et al. (2014). *A new device for simple and accurate urinary pH testing by the Stone-former patient.* Springerplus. [PMID 24839588](https://pubmed.ncbi.nlm.nih.gov/24839588/)
- Ojanperä I et al. (2005). *Application of accurate mass measurement to urine drug screening.* J Anal Toxicol. [PMID 15808011](https://pubmed.ncbi.nlm.nih.gov/15808011/)
- Megahed AA et al. (2019). *Evaluation of the analytical performance of a portable ion-selective electrode meter for measuring whole-blood, plasma, milk, abomasal-fluid, and urine sodium concentrations in cattle.* J Dairy Sci. [PMID 31202658](https://pubmed.ncbi.nlm.nih.gov/31202658/)
- Thomas A et al. (2013). *Quantification of AICAR-ribotide concentrations in red blood cells by means of LC-MS/MS.* Anal Bioanal Chem. [PMID 23828211](https://pubmed.ncbi.nlm.nih.gov/23828211/)
- Ahmed MJ et al. (2002). *A simple spectrophotometric method for the determination of copper in industrial, environmental, biological and soil samples using 2,5-dimercapto-1,3,4-thiadiazole.* Anal Sci. [PMID 12137377](https://pubmed.ncbi.nlm.nih.gov/12137377/)
- Kelly AM et al. (2004). *Agreement between bicarbonate measured on arterial and venous blood gases.* Emerg Med Australas. [PMID 15537402](https://pubmed.ncbi.nlm.nih.gov/15537402/)
- Treger R et al. (2010). *Agreement between central venous and arterial blood gas measurements in the intensive care unit.* Clin J Am Soc Nephrol. [PMID 20019117](https://pubmed.ncbi.nlm.nih.gov/20019117/)
- Esmaeilivand M et al. (2017). *Agreement and Correlation between Arterial and Central Venous Blood Gas Following Coronary Artery Bypass Graft Surgery.* J Clin Diagn Res. [PMID 28511435](https://pubmed.ncbi.nlm.nih.gov/28511435/)
- Kelly AM (2010). *Review article: Can venous blood gas analysis replace arterial in emergency medical care.* Emerg Med Australas. [PMID 21143397](https://pubmed.ncbi.nlm.nih.gov/21143397/)
- Bloom BM et al. (2014). *The role of venous blood gas in the emergency department: a systematic review and meta-analysis.* Eur J Emerg Med. [PMID 23903783](https://pubmed.ncbi.nlm.nih.gov/23903783/)
- Kadwa AR et al. (2022). *Agreement between arterial and central venous blood pH and its contributing variables in anaesthetized dogs with respiratory acidosis.* Vet Anaesth Analg. [PMID 35292229](https://pubmed.ncbi.nlm.nih.gov/35292229/)
- Kelly AM (2006). *The case for venous rather than arterial blood gases in diabetic ketoacidosis.* Emerg Med Australas. [PMID 16454777](https://pubmed.ncbi.nlm.nih.gov/16454777/)

## 24. Quantitative systems pharmacology, acid-base modelling and mrgsolve

Methodological background for the model itself.

- Wolf MB (2024). *Mechanisms of whole body, respiratory, acid-base buffering: a first computer-model test of three physicochemical, acid-base theories.* J Appl Physiol (1985). [PMID 38752284](https://pubmed.ncbi.nlm.nih.gov/38752284/)
- Ghallab A et al. (2024). *Inhibition of the renal apical sodium dependent bile acid transporter prevents cholemic nephropathy in mice with obstructive cholestasis.* J Hepatol. [PMID 37939855](https://pubmed.ncbi.nlm.nih.gov/37939855/)
- Weinstein AM (2001). *A mathematical model of rat cortical collecting duct: determinants of the transtubular potassium gradient.* Am J Physiol Renal Physiol. [PMID 11352847](https://pubmed.ncbi.nlm.nih.gov/11352847/)
- Weinstein AM (2021). *A mathematical model of the rat kidney. III. Ammonia transport.* Am J Physiol Renal Physiol. [PMID 33779315](https://pubmed.ncbi.nlm.nih.gov/33779315/)
- Wang XP et al. (2022). *Bile acids regulate the epithelial Na(+) channel in native tissues through direct binding at multiple sites.* J Physiol. [PMID 36071685](https://pubmed.ncbi.nlm.nih.gov/36071685/)
- Radhakrishnan VM et al. (2013). *Post-translational loss of renal TRPV5 calcium channel expression, Ca(2+) wasting, and bone loss in experimental colitis.* Gastroenterology. [PMID 23747339](https://pubmed.ncbi.nlm.nih.gov/23747339/)
- Ficici E et al. (2017). *Asymmetry of inverted-topology repeats in the AE1 anion exchanger suggests an elevator-like mechanism.* J Gen Physiol. [PMID 29167180](https://pubmed.ncbi.nlm.nih.gov/29167180/)
- Luo G et al. (2010). *In silico prediction of biliary excretion of drugs in rats based on physicochemical properties.* Drug Metab Dispos. [PMID 19995888](https://pubmed.ncbi.nlm.nih.gov/19995888/)
- Slusarz MJ et al. (2006). *Investigation of mechanism of desmopressin binding in vasopressin V2 receptor versus vasopressin V1a and oxytocin receptors: molecular dynamics simulation of the agonist-bound state in the membrane-aqueous system.* Biopolymers. [PMID 16333859](https://pubmed.ncbi.nlm.nih.gov/16333859/)
- Luo Z et al. (2019). *Gene Expression Signatures Associated With Survival Times of Pediatric Patients With Biliary Atresia Identify Potential Therapeutic Agents.* Gastroenterology. [PMID 31228442](https://pubmed.ncbi.nlm.nih.gov/31228442/)
- Gaweda AE et al. (2021). *Development of a quantitative systems pharmacology model of chronic kidney disease: metabolic bone disorder.* Am J Physiol Renal Physiol. [PMID 33308018](https://pubmed.ncbi.nlm.nih.gov/33308018/)
- Stodtmann S et al. (2021). *Validation of a quantitative systems pharmacology model of calcium homeostasis using elagolix Phase 3 clinical trial data in women with endometriosis.* Clin Transl Sci. [PMID 33963686](https://pubmed.ncbi.nlm.nih.gov/33963686/)
- Gaweda A et al. (2024). *Leveraging quantitative systems pharmacology and artificial intelligence to advance treatment of chronic kidney disease mineral bone disorder.* Am J Physiol Renal Physiol. [PMID 38961848](https://pubmed.ncbi.nlm.nih.gov/38961848/)
- Yang JH et al. (2012). *The effect of bone morphogenic protein-2-coated tri-calcium phosphate/hydroxyapatite on new bone formation in a rat model of femoral distraction osteogenesis.* Cytotherapy. [PMID 22122301](https://pubmed.ncbi.nlm.nih.gov/22122301/)
- Gaweda AE et al. (2022). *Artificial intelligence-guided precision treatment of chronic kidney disease-mineral bone disorder.* CPT Pharmacometrics Syst Pharmacol. [PMID 35920131](https://pubmed.ncbi.nlm.nih.gov/35920131/)
- Lee SJ et al. (2017). *The use of heparin chemistry to improve dental osteogenesis associated with implants.* Carbohydr Polym. [PMID 27987891](https://pubmed.ncbi.nlm.nih.gov/27987891/)
- Gomes PS et al. (2021). *The Osteogenic Assessment of Mineral Trioxide Aggregate-based Endodontic Sealers in an Organotypic Ex Vivo Bone Development Model.* J Endod. [PMID 34126159](https://pubmed.ncbi.nlm.nih.gov/34126159/)
- Chou YF et al. (2005). *In vitro response of MC3T3-E1 pre-osteoblasts within three-dimensional apatite-coated PLGA scaffolds.* J Biomed Mater Res B Appl Biomater. [PMID 16001421](https://pubmed.ncbi.nlm.nih.gov/16001421/)
- Cheng A et al. (2023). *Pterosin sesquiterpenoids from Pteris laeta Wall. ex Ettingsh. protect cells from glutamate excitotoxicity by modulating mitochondrial signals.* J Ethnopharmacol. [PMID 36822346](https://pubmed.ncbi.nlm.nih.gov/36822346/)
- Itou T et al. (2014). *Cystathionine γ-lyase accelerates osteoclast differentiation: identification of a novel regulator of osteoclastogenesis by proteomic analysis.* Arterioscler Thromb Vasc Biol. [PMID 24357058](https://pubmed.ncbi.nlm.nih.gov/24357058/)
- Lambert F et al. (2013). *A comparison of three calcium phosphate-based space fillers in sinus elevation: a study in rabbits.* Int J Oral Maxillofac Implants. [PMID 23527340](https://pubmed.ncbi.nlm.nih.gov/23527340/)
- Pathak YV et al. (1990). *Prevention of calcification of glutaraldehyde pretreated bovine pericardium through controlled release polymeric implants: studies of Fe3+, Al3+, protamine sulphate and levamisole.* Biomaterials. [PMID 2128616](https://pubmed.ncbi.nlm.nih.gov/2128616/)

---

## Calibration anchors actually used, and what the model returns

| Observation | Source | Model |
|---|---|---|
| SoC (median 3 intakes/day) -> ADV7103 twice daily, responder rate on plasma bicarbonate 43% -> 90% | B21CS, n=37, PMID 32712761 | see the decomposition in `README.md`: matched-dose schedule pharmacology accounts for only part of this, and the model says so explicitly |
| Urine Ca/citrate fell below the lithogenic threshold in 56% of previous SoC non-responders | B21CS, PMID 32712761 | under-predicted; stated as a known misfit |
| Plasma bicarbonate 22.0 +/- 3.2 -> 22.6 +/- 2.5 mmol/L over ~6 y (NS) | B22CS, n=30, PMID 40801206 | reproduced |
| Height z-score -0.6 +/- 1.0 -> -0.3 +/- 1.0 (p=0.03) | B22CS, PMID 40801206 | reproduced by construction of the catch-up term |
| Lumbar BMD z-score -1.1 +/- 1.0 -> -0.8 +/- 1.0 (p=0.005) | B22CS, PMID 40801206 | reproduced by construction of the catch-up term |
| eGFR 105 +/- 17 -> 104 +/- 20 mL/min/1.73 m2 over 6 y, stable under good control | B22CS, PMID 40801206 | reproduced |
| > 80% of ADULT dRTA patients carry KDIGO stage >= 2 | B22CS, PMID 40801206 | qualitatively reproduced via the nephrocalcinosis -> fibrosis -> nephron-loss loop |
| ADV7103 citrate granules release over 2-3 h, bicarbonate granules over 10-12 h | Guittet 2020, PMID 32811843 | used directly as kr_citPR = 0.90/h and kr_bicPR = 0.215/h |
| Urine pH held above 7 for 24 h at 1.44 mEq/kg twice daily in healthy adults; placebo below 6 | Guittet 2020, PMID 32811843 | scenario S23 |
| Urinary Ca varies with net acid excretion by 0.035 mmol/mEq | Lemann 1999, PMID 9873210 | **0.0357 mmol/mEq** across a dietary acid titration — a validation target, not a fitted coefficient |
| Normal urine pH falls below 5.45 after NH4Cl loading; dRTA cannot reach 5.5 | PMID 34973150, 2081516, 30139458 | normal reaches 4.61; incomplete dRTA stalls at 6.21; complete dRTA at 6.88 |
| Fractional excretion of bicarbonate < 3% in dRTA vs > 15% in proximal RTA | PMID 2081516 | reproduced (< 1% on the distal lesion) |
| Normal net acid excretion ~43-60 mEq/day, maximum 300-450 mEq/day under load | PMID 30745301 and section 3 | 67 mEq/day at a normal Western acid load, 254 mEq/day at 4x dietary acid |
| Prevalence of incomplete dRTA among idiopathic calcium stone formers: 25.4% of women, 13.6% of men | PMID 34973150 | context for scenario S05, not a fitted target |

*463 references across 24 sections.*
