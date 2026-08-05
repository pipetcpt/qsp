# Pyruvate Kinase Deficiency (PKLR) — QSP Model References

> 267 references, every PMID retrieved by direct query against the NCBI PubMed
> E-utilities API while this model was being built. Nothing here is recalled from
> memory. Sections follow the structure of the mechanistic map
> (`pkd_qsp_model.dot`) and the mrgsolve model (`pkd_mrgsolve_model.R`).

## §0 — The load-bearing references, and exactly what each one carries

The model has four fitted parameters and roughly two hundred that are either
taken from literature or back-calculated from the normal erythrocyte operating
point. The list below is not "key papers"; it is the set of papers a specific
number in the model came from, with the number named. If one of these is wrong,
a stated result moves.

| Quantity in the model | Value used | Source |
|---|---|---|
| Reduced glycolysis / 2,3-BPG kinetic structure and reference metabolite set | ATP 1.70, ADP 0.22, 2,3-BPG 4.80, 3-PG 0.060, PEP 0.017, FBP 0.012 mM; flux 1.60 mmol/L RBC/h | Mulquiney & Kuchel 1999, [PMID 10477269](https://pubmed.ncbi.nlm.nih.gov/10477269/) and [PMID 10477270](https://pubmed.ncbi.nlm.nih.gov/10477270/) |
| Mitapivat haemoglobin response rate and endpoint definition | 16/40 (40%) vs 0/40 placebo, ≥1.5 g/dL sustained at ≥2 of weeks 16/20/24 | ACTIVATE, Al-Samkari 2022, [PMID 35417638](https://pubmed.ncbi.nlm.nih.gov/35417638/) |
| Time-to-response, and the genotype restriction on response | median 10 d to first >1.0 g/dL (range 7–187); responses **only** with ≥1 missense variant, correlated with baseline PK-R protein | DRIVE-PK, Grace 2019, [PMID 31483964](https://pubmed.ncbi.nlm.nih.gov/31483964/) |
| Transfusion-burden endpoint | ≥33% reduction in 10/27 (37%) | ACTIVATE-T, Glenthøj 2022, [PMID 35988546](https://pubmed.ncbi.nlm.nih.gov/35988546/) |
| Splenectomy effect; complication frequencies | median Hb **+1.6 g/dL**; iron overload 48%, gallstones 45%, post-splenectomy thrombosis 11%, neonatal phototherapy 93%, exchange transfusion 46%; 48% of splenectomies without simultaneous cholecystectomy later need one | PK Deficiency Natural History Study, Grace 2018, [PMID 29549173](https://pubmed.ncbi.nlm.nih.gov/29549173/) |
| Iron axis response to mitapivat (direction and magnitude of every analyte) | hepcidin +4770 ng/L, erythroferrone −9835 ng/L, sTfR −56.0 nmol/L, EPO −32.9 IU/L at wk 24; liver iron −2.0 mg Fe/g dw at wk 96 | van Beers 2024, [PMID 38330179](https://pubmed.ncbi.nlm.nih.gov/38330179/) |
| Splenic **sequestration** as distinct from destruction | organ sequestration measured alongside red cell life-span | Nathan 1968, [PMID 5634483](https://pubmed.ncbi.nlm.nih.gov/5634483/) |
| Dose-dependent ↑ATP and ↓2,3-BPG in healthy volunteers — the observation the whole "two opposing trunks" structure rests on | glycolytic intermediates shifted consistently with pathway activation at all multiple-dose levels | Yang 2019, [PMID 30091852](https://pubmed.ncbi.nlm.nih.gov/30091852/) |
| Genotype–phenotype relation and molecular heterogeneity | missense vs non-missense classification driving the "activatable protein" parameter | Zanella 2007 [PMID 17360088](https://pubmed.ncbi.nlm.nih.gov/17360088/); Bianchi 2020 [PMID 32043619](https://pubmed.ncbi.nlm.nih.gov/32043619/), [PMID 33054047](https://pubmed.ncbi.nlm.nih.gov/33054047/) |
| UGT1A1 co-inheritance multiplying gallstone risk in chronic haemolysis | used to set the UGT activity multiplier (1.0 / 0.7 / 0.3) | del Giudice 1999, [PMID 10498597](https://pubmed.ncbi.nlm.nih.gov/10498597/) |
| Statistical structure of red cell survival (age-dependent vs random hazard) | motivates the 14-cohort age axis rather than a single pool | Korell 2011, [PMID 20950630](https://pubmed.ncbi.nlm.nih.gov/20950630/) |
| Gene-therapy correction thresholds (fraction of corrected output needed) | sets the GTFRAC values simulated | Garcia-Gomez 2016 [PMID 27138040](https://pubmed.ncbi.nlm.nih.gov/27138040/); Navarro 2021 [PMID 34514027](https://pubmed.ncbi.nlm.nih.gov/34514027/) |

**The single most load-bearing parameter is NOT in this table**, and that is worth
saying plainly: the 3-PG inhibition constant of the 2,3-BPG phosphatase
(`KIPG3` = 30 µM). It is what makes 2,3-BPG rise 2-fold rather than 9% when
pyruvate kinase activity falls, and the entire oxygen-transport argument scales
with the size of that excursion. It is the least well pinned number in the model.
The break-even table in §3.2 of `pkd_reference_output.txt` exists so the
conclusion can be re-scored against a different 2,3-BPG estimate without
rebuilding anything.

---

### 1. Pyruvate kinase deficiency: disease, genetics, natural history

*31 references*

- Tanaka KR. Pyruvate kinase deficiency. *Semin Hematol 1971*. [PMID 4942588](https://pubmed.ncbi.nlm.nih.gov/4942588/)
- Waller HD. [Enzyme deficiencies in glycolysis and nucleotide metabolism of red blood cells in nonspherocytic hemolytic anemia (author's transl)]. *Klin Wochenschr 1976*. [PMID 184346](https://pubmed.ncbi.nlm.nih.gov/184346/)
- Zanella A. Hereditary pyruvate kinase deficiency: role of the abnormal enzyme in red cell pathophysiology. *Br J Haematol 1978*. [PMID 728372](https://pubmed.ncbi.nlm.nih.gov/728372/)
- Beutler E. Mutations in pyruvate kinase. *Hum Mutat 1996*. [PMID 8664896](https://pubmed.ncbi.nlm.nih.gov/8664896/)
- Orosz F. Triosephosphate isomerase deficiency: predictions and facts. *J Theor Biol 1996*. [PMID 8944178](https://pubmed.ncbi.nlm.nih.gov/8944178/)
- Orosz F. Triosephosphate isomerase deficiency: facts and doubts. *IUBMB Life 2006*. [PMID 17424909](https://pubmed.ncbi.nlm.nih.gov/17424909/)
- Zanella A. Pyruvate kinase deficiency. *Haematologica 2007*. [PMID 17550841](https://pubmed.ncbi.nlm.nih.gov/17550841/)
- Vijayan VK. Parasitic lung infections. *Curr Opin Pulm Med 2009*. [PMID 19276810](https://pubmed.ncbi.nlm.nih.gov/19276810/)
- Bakhramov SM. [Erythrocytic enzymopathy in Uzbekistan]. *Lik Sprava 2011*. [PMID 22768742](https://pubmed.ncbi.nlm.nih.gov/22768742/)
- Jamwal M. Next-generation sequencing unravels homozygous mutation in glucose-6-phosphate isomerase, GPIc.1040G>A (p.Arg347His) causing hemolysis in an Indi. *Clin Chim Acta 2017*. [PMID 28223188](https://pubmed.ncbi.nlm.nih.gov/28223188/)
- Srivastava D. Structural Investigation of a Dimeric Variant of Pyruvate Kinase Muscle Isoform 2. *Biochemistry 2017*. [PMID 29182273](https://pubmed.ncbi.nlm.nih.gov/29182273/)
- Grace RF. Clinical spectrum of pyruvate kinase deficiency: data from the Pyruvate Kinase Deficiency Natural History Study. *Blood 2018*. [PMID 29549173](https://pubmed.ncbi.nlm.nih.gov/29549173/)
- Grace RF. How we manage patients with pyruvate kinase deficiency. *Br J Haematol 2019*. [PMID 30681718](https://pubmed.ncbi.nlm.nih.gov/30681718/)
- Segal J. Low catalytic activity is insufficient to induce disease pathology in triosephosphate isomerase deficiency. *J Inherit Metab Dis 2019*. [PMID 31111503](https://pubmed.ncbi.nlm.nih.gov/31111503/)
- Srivastava D. Mechanistic and Structural Insights into Cysteine-Mediated Inhibition of Pyruvate Kinase Muscle Isoform 2. *Biochemistry 2019*. [PMID 31386812](https://pubmed.ncbi.nlm.nih.gov/31386812/)
- Agarwal AM. Laboratory approach to investigation of anemia with a focus on pyruvate kinase deficiency. *Int J Lab Hematol 2020*. [PMID 32543069](https://pubmed.ncbi.nlm.nih.gov/32543069/)
- Grace RF. Management of pyruvate kinase deficiency in children and adults. *Blood 2020*. [PMID 32702739](https://pubmed.ncbi.nlm.nih.gov/32702739/)
- Nandi S. Structural basis for allosteric regulation of pyruvate kinase M2 by phosphorylation and acetylation. *J Biol Chem 2020*. [PMID 33453989](https://pubmed.ncbi.nlm.nih.gov/33453989/)
- Secrest MH. Prevalence of pyruvate kinase deficiency: A systematic literature review. *Eur J Haematol 2020*. [PMID 32279356](https://pubmed.ncbi.nlm.nih.gov/32279356/)
- Chhipa AS. Targeting pyruvate kinase muscle isoform 2 (PKM2) in cancer: What do we know so far?. *Life Sci 2021*. [PMID 34102192](https://pubmed.ncbi.nlm.nih.gov/34102192/)
- Dongerdiye R. Novel pathogenic variant c.2714C>A (p. Thr905Lys) in the HK1 gene causing severe haemolytic anaemia with developmental delay in an Indian family. *J Clin Pathol 2021*. [PMID 33361148](https://pubmed.ncbi.nlm.nih.gov/33361148/)
- Falak S. Molecular cloning, expression in Escherichia coli and structural-functional analysis of a pyruvate kinase from Pyrobaculum calidifontis. *Int J Biol Macromol 2022*. [PMID 35472364](https://pubmed.ncbi.nlm.nih.gov/35472364/)
- Fattizzo B. Pyruvate Kinase Deficiency: Current Challenges and Future Prospects. *J Blood Med 2022*. [PMID 36072510](https://pubmed.ncbi.nlm.nih.gov/36072510/)
- Hou DY. OGA activated glycopeptide-based nano-activator to activate PKM2 tetramerization for switching catabolic pathways and sensitizing chemotherapy re. *Biomaterials 2022*. [PMID 35462306](https://pubmed.ncbi.nlm.nih.gov/35462306/)
- Johnson S. Diagnosis, monitoring, and management of pyruvate kinase deficiency in children. *Pediatr Blood Cancer 2022*. [PMID 35452178](https://pubmed.ncbi.nlm.nih.gov/35452178/)
- VanDemark AP. Itavastatin and resveratrol increase triosephosphate isomerase protein in a newly identified variant of TPI deficiency. *Dis Model Mech 2022*. [PMID 35315486](https://pubmed.ncbi.nlm.nih.gov/35315486/)
- Grace RF. The Pyruvate Kinase Deficiency Global Longitudinal (Peak) Registry: rationale and study design. *BMJ Open 2023*. [PMID 36958777](https://pubmed.ncbi.nlm.nih.gov/36958777/)
- Al-Samkari H. Diagnosis and management of pyruvate kinase deficiency: international expert guidelines. *Lancet Haematol 2024*. [PMID 38330977](https://pubmed.ncbi.nlm.nih.gov/38330977/)
- Holme S. Glucose phosphate isomerase deficiency demasked by whole-genome sequencing: a case report. *J Med Case Rep 2024*. [PMID 38539245](https://pubmed.ncbi.nlm.nih.gov/38539245/)
- Ahmed SH. Efficacy and safety of pyruvate kinase activator in treating hemolytic anemias: a systematic review. *Expert Rev Hematol 2025*. [PMID 40526104](https://pubmed.ncbi.nlm.nih.gov/40526104/)
- Williams A. TPI deficiency: A case report and review of the literature. *Mol Genet Metab 2025*. [PMID 40897044](https://pubmed.ncbi.nlm.nih.gov/40897044/)

### 2. Erythrocyte glycolysis, the Rapoport-Luebering shunt and kinetic modelling

*8 references*

- Mueggler PA. Postnatal regulation of canine oxygen delivery: control of erythrocyte 2,3-DPG levels. *Am J Physiol 1982*. [PMID 6278964](https://pubmed.ncbi.nlm.nih.gov/6278964/)
- Mulquiney PJ. Model of 2,3-bisphosphoglycerate metabolism in the human erythrocyte based on detailed enzyme kinetic equations: in vivo kinetic characterization. *Biochem J 1999*. [PMID 10477268](https://pubmed.ncbi.nlm.nih.gov/10477268/)
- Kor DJ. Red blood cell storage lesion. *Bosn J Basic Med Sci 2009*. [PMID 19912115](https://pubmed.ncbi.nlm.nih.gov/19912115/)
- Peng Z. Erythrocyte Adenosine A2B Receptor-Mediated AMPK Activation: A Missing Component Counteracting CKD by Promoting Oxygen Delivery. *J Am Soc Nephrol 2019*. [PMID 31278195](https://pubmed.ncbi.nlm.nih.gov/31278195/)
- D'Alessandro A. Erythrocyte adaptive metabolic reprogramming under physiological and pathological hypoxia. *Curr Opin Hematol 2020*. [PMID 32141895](https://pubmed.ncbi.nlm.nih.gov/32141895/)
- Xie T. Erythrocyte Metabolic Reprogramming by Sphingosine 1-Phosphate in Chronic Kidney Disease and Therapies. *Circ Res 2020*. [PMID 32284030](https://pubmed.ncbi.nlm.nih.gov/32284030/)
- Chen C. Erythrocyte ENT1-AMPD3 Axis is an Essential Purinergic Hypoxia Sensor and Energy Regulator Combating CKD in a Mouse Model. *J Am Soc Nephrol 2023*. [PMID 37725437](https://pubmed.ncbi.nlm.nih.gov/37725437/)
- Jaafar LS. 2,3-Diphosphoglycerate: the forgotten metabolic regulator of oxygen affinity. *Br J Nutr 2025*. [PMID 41070558](https://pubmed.ncbi.nlm.nih.gov/41070558/)

### 3. 2,3-BPG, haemoglobin oxygen affinity and oxygen delivery

*22 references*

- Woodson RD. Effect of increased blood oxygen affinity on work performance of rats. *J Clin Invest 1973*. [PMID 4748508](https://pubmed.ncbi.nlm.nih.gov/4748508/)
- Makino N. Oxygen equilibria of hybrid-heme hemoglobins containing proto- and mesoheme groups. On the nonequivalence of alpha and beta chains. *J Biol Chem 1978*. [PMID 24050](https://pubmed.ncbi.nlm.nih.gov/24050/)
- Jaffé ER. Methemoglobin pathophysiology. *Prog Clin Biol Res 1981*. [PMID 7022466](https://pubmed.ncbi.nlm.nih.gov/7022466/)
- Jaffé ER. Methaemoglobinaemia. *Clin Haematol 1981*. [PMID 7011627](https://pubmed.ncbi.nlm.nih.gov/7011627/)
- Buerk DG. An evaluation of Easton's paradigm for the oxyhemoglobin equilibrium curve. *Adv Exp Med Biol 1984*. [PMID 6534109](https://pubmed.ncbi.nlm.nih.gov/6534109/)
- Madsen H. Red cell 2,3-diphosphoglycerate and hemoglobin--oxygen affinity during normal pregnancy. *Acta Obstet Gynecol Scand 1984*. [PMID 6496042](https://pubmed.ncbi.nlm.nih.gov/6496042/)
- Edelstone DI. Effects of maternal anemia on cardiac output, systemic oxygen consumption, and regional blood flow in pregnant sheep. *Am J Obstet Gynecol 1987*. [PMID 3826225](https://pubmed.ncbi.nlm.nih.gov/3826225/)
- Wilkerson DK. Oxygen extraction ratio: a valid indicator of myocardial metabolism in anemia. *J Surg Res 1987*. [PMID 3586630](https://pubmed.ncbi.nlm.nih.gov/3586630/)
- Wilkerson DK. Limits of cardiac compensation in anemic baboons. *Surgery 1988*. [PMID 3375993](https://pubmed.ncbi.nlm.nih.gov/3375993/)
- King LG. Anemia of chronic renal failure in dogs. *J Vet Intern Med 1992*. [PMID 1432900](https://pubmed.ncbi.nlm.nih.gov/1432900/)
- Mansouri A. Concise review: methemoglobinemia. *Am J Hematol 1993*. [PMID 8416301](https://pubmed.ncbi.nlm.nih.gov/8416301/)
- Fox GA. Hematocrit modifies the circulatory control of systemic and myocardial oxygen utilization in septic sheep. *Crit Care Med 1994*. [PMID 8124998](https://pubmed.ncbi.nlm.nih.gov/8124998/)
- Ostgaard G. [Perioperative and postoperative normovolemic anemia. Physiological compensation, monitoring and risk evaluation]. *Tidsskr Nor Laegeforen 1996*. [PMID 8553339](https://pubmed.ncbi.nlm.nih.gov/8553339/)
- Papassotiriou I. Synthesized allosteric effectors of the hemoglobin molecule: a possible mechanism for improved erythrocyte oxygen release capability in hemoglobi. *Exp Hematol 1998*. [PMID 9728926](https://pubmed.ncbi.nlm.nih.gov/9728926/)
- Metivier F. Pathophysiology of anaemia: focus on the heart and blood vessels. *Nephrol Dial Transplant 2000*. [PMID 11032352](https://pubmed.ncbi.nlm.nih.gov/11032352/)
- Eckardt KU. Anaemia in end-stage renal disease: pathophysiological considerations. *Nephrol Dial Transplant 2001*. [PMID 11590249](https://pubmed.ncbi.nlm.nih.gov/11590249/)
- Habler O. Hyperoxia in extreme hemodilution. *Eur Surg Res 2002*. [PMID 11867921](https://pubmed.ncbi.nlm.nih.gov/11867921/)
- Villela NR. Microcirculatory effects of changing blood hemoglobin oxygen affinity during hemorrhagic shock resuscitation in an experimental model. *Shock 2009*. [PMID 18948853](https://pubmed.ncbi.nlm.nih.gov/18948853/)
- Boĭko NV. [Biochemical characteristics of compensation of posthemorrhagic anemia in patients presenting with nasal bleeding]. *Vestn Otorinolaringol 2010*. [PMID 21105337](https://pubmed.ncbi.nlm.nih.gov/21105337/)
- Saleh MC. NADH-dependent cytochrome b5 reductase and NADPH methemoglobin reductase activity in the erythrocytes of Oncorhynchus mykiss. *Fish Physiol Biochem 2012*. [PMID 22733093](https://pubmed.ncbi.nlm.nih.gov/22733093/)
- Maeda S. Methemoglobin reduction by NADH-cytochrome b(5) reductase in Propsilocerus akamusi larvae. *Comp Biochem Physiol B Biochem Mol Biol 2015*. [PMID 25829149](https://pubmed.ncbi.nlm.nih.gov/25829149/)
- Ericsson A. FT-4202, a selective pyruvate kinase R activator for sickle cell disease. *Exp Hematol 2025*. [PMID 39549740](https://pubmed.ncbi.nlm.nih.gov/39549740/)

### 4. Red cell ATP, cation homeostasis, deformability and lifespan

*2 references*

- Glogowska E. Mutations in the Gardos channel (KCNN4) are associated with hereditary xerocytosis. *Blood 2015*. [PMID 26198474](https://pubmed.ncbi.nlm.nih.gov/26198474/)
- Vives-Corrons JL. Hereditary Spherocytosis: Linking Ion Transport Defects to Osmotic Gradient Ektacytometry Profiles-A Review. *Int J Mol Sci 2026*. [PMID 41596371](https://pubmed.ncbi.nlm.nih.gov/41596371/)

### 5. The spleen: filtration, sequestration and splenectomy

*25 references*

- Cokelet GR. Dynamics of erythrocyte motion in filtration tests and in vivo flow. *Scand J Clin Lab Invest Suppl 1981*. [PMID 6948404](https://pubmed.ncbi.nlm.nih.gov/6948404/)
- Weiss L. The red pulp of the spleen: structural basis of blood flow. *Clin Haematol 1983*. [PMID 6352110](https://pubmed.ncbi.nlm.nih.gov/6352110/)
- Shatney CH. Complications of splenectomy. *Acta Anaesthesiol Belg 1987*. [PMID 3327338](https://pubmed.ncbi.nlm.nih.gov/3327338/)
- Jakubovský J. [Functional morphology of the spleen]. *Bratisl Lek Listy 1995*. [PMID 8624745](https://pubmed.ncbi.nlm.nih.gov/8624745/)
- Brew I. Post-splenectomy sepsis--the role of occupational health. *Occup Med (Lond) 1996*. [PMID 8695777](https://pubmed.ncbi.nlm.nih.gov/8695777/)
- Chrobák L. [Splenomegaly (clinical importance, diagnosis and therapy)]. *Vnitr Lek 2002*. [PMID 12061183](https://pubmed.ncbi.nlm.nih.gov/12061183/)
- Cohen AR. Thalassemia. *Hematology Am Soc Hematol Educ Program 2004*. [PMID 15561674](https://pubmed.ncbi.nlm.nih.gov/15561674/)
- Kaplinsky C. Post-splenectomy antibiotic prophylaxis--unfinished story: to treat or not to treat?. *Pediatr Blood Cancer 2006*. [PMID 16933244](https://pubmed.ncbi.nlm.nih.gov/16933244/)
- Soyer T. Portal vein thrombosis after splenectomy in pediatric hematologic disease: risk factors, clinical features, and outcome. *J Pediatr Surg 2006*. [PMID 17101367](https://pubmed.ncbi.nlm.nih.gov/17101367/)
- Buesing KL. Partial splenectomy for hereditary spherocytosis: a multi-institutional review. *J Pediatr Surg 2011*. [PMID 21238662](https://pubmed.ncbi.nlm.nih.gov/21238662/)
- Tuchscherer D. [What is the spleen needed for?]. *Ther Umsch 2013*. [PMID 23454560](https://pubmed.ncbi.nlm.nih.gov/23454560/)
- Taniguchi LU. Overwhelming post-splenectomy infection: narrative review of the literature. *Surg Infect (Larchmt) 2014*. [PMID 25318011](https://pubmed.ncbi.nlm.nih.gov/25318011/)
- Guizzetti L. Total versus partial splenectomy in pediatric hereditary spherocytosis: A systematic review and meta-analysis. *Pediatr Blood Cancer 2016*. [PMID 27300151](https://pubmed.ncbi.nlm.nih.gov/27300151/)
- Rogulski R. Laparoscopic splenectomy for hereditary spherocytosis-preliminary report. *Eur J Haematol 2016*. [PMID 26268883](https://pubmed.ncbi.nlm.nih.gov/26268883/)
- White NJ. Malaria parasite clearance. *Malar J 2017*. [PMID 28231817](https://pubmed.ncbi.nlm.nih.gov/28231817/)
- Tahir F. Post-splenectomy Sepsis: A Review of the Literature. *Cureus 2020*. [PMID 32195065](https://pubmed.ncbi.nlm.nih.gov/32195065/)
- Horvat M. Audit of Post-Splenectomy Prophylaxis in a Single Tertiary Center in Slovenia: Where Are We and What Should Be Done?. *Surg Infect (Larchmt) 2021*. [PMID 32639189](https://pubmed.ncbi.nlm.nih.gov/32639189/)
- Celik SS. Clinical Characteristics and Treatment Outcome of Hereditary Spherocytosis: A Single Center's Experience. *Sisli Etfal Hastan Tip Bul 2023*. [PMID 38268662](https://pubmed.ncbi.nlm.nih.gov/38268662/)
- Saldanha A. The immune thrombocytopenia paradox: Should we be concerned about thrombosis in ITP?. *Thromb Res 2024*. [PMID 39137700](https://pubmed.ncbi.nlm.nih.gov/39137700/)
- Tang X. The efficacy of partial versus total splenectomy in the treatment of hereditary spherocytosis in children: a systematic review and meta-analysis. *Pediatr Surg Int 2024*. [PMID 39470805](https://pubmed.ncbi.nlm.nih.gov/39470805/)
- Cappellini MD. Thalassemia and hypercoagulability. *Hematology Am Soc Hematol Educ Program 2025*. [PMID 41348010](https://pubmed.ncbi.nlm.nih.gov/41348010/)
- Lluís N. Management of splanchnic venous thrombosis after splenectomy in hematologic diseases. systematic review, meta-analysis and consensus guidelines. *Int J Surg 2025*. [PMID 41417989](https://pubmed.ncbi.nlm.nih.gov/41417989/)
- Saffioti NA. Escherichia coli α-hemolysin induces red blood cell retention in a microfluidic spleen-like device. *Biophys J 2025*. [PMID 40077968](https://pubmed.ncbi.nlm.nih.gov/40077968/)
- Tefferi A. Essential Thrombocythemia: A Review. *JAMA 2025*. [PMID 39869325](https://pubmed.ncbi.nlm.nih.gov/39869325/)
- Tessier B. Long-Term Outcome After Partial Splenectomy Compared to Total Splenectomy in Children With Spherocytosis. *J Pediatr Surg 2025*. [PMID 40774475](https://pubmed.ncbi.nlm.nih.gov/40774475/)

### 6. Erythropoiesis, erythropoietin and reticulocyte kinetics

*37 references*

- Anderson MJ. Human parvovirus infections. *J Virol Methods 1987*. [PMID 2822752](https://pubmed.ncbi.nlm.nih.gov/2822752/)
- Rotbart HA. Human parvovirus infections. *Annu Rev Med 1990*. [PMID 2158761](https://pubmed.ncbi.nlm.nih.gov/2158761/)
- Bender MA. Sickle Cell Disease. *1993*. [PMID 20301551](https://pubmed.ncbi.nlm.nih.gov/20301551/)
- Kalfa TA. EPB42-Related Hereditary Spherocytosis. *1993*. [PMID 24624460](https://pubmed.ncbi.nlm.nih.gov/24624460/)
- Brown KE. Parvovirus B19 infection and hematopoiesis. *Blood Rev 1995*. [PMID 8563519](https://pubmed.ncbi.nlm.nih.gov/8563519/)
- Mittman N. Reticulocyte hemoglobin content predicts functional iron deficiency in hemodialysis patients receiving rHuEPO. *Am J Kidney Dis 1997*. [PMID 9398141](https://pubmed.ncbi.nlm.nih.gov/9398141/)
- Choi JW. Change in erythropoiesis with gestational age during pregnancy. *Ann Hematol 2001*. [PMID 11233772](https://pubmed.ncbi.nlm.nih.gov/11233772/)
- Heegaard ED. Human parvovirus B19. *Clin Microbiol Rev 2002*. [PMID 12097253](https://pubmed.ncbi.nlm.nih.gov/12097253/)
- Meyer O. Parvovirus B19 and autoimmune diseases. *Joint Bone Spine 2003*. [PMID 12639611](https://pubmed.ncbi.nlm.nih.gov/12639611/)
- Koury MJ. New insights into erythropoiesis: the roles of folate, vitamin B12, and iron. *Annu Rev Nutr 2004*. [PMID 15189115](https://pubmed.ncbi.nlm.nih.gov/15189115/)
- Eckardt KU. Regulation of erythropoietin production. *Eur J Clin Invest 2005*. [PMID 16281953](https://pubmed.ncbi.nlm.nih.gov/16281953/)
- Johnson DW. Erythropoiesis-stimulating agent hyporesponsiveness. *Nephrology (Carlton) 2007*. [PMID 17635745](https://pubmed.ncbi.nlm.nih.gov/17635745/)
- Tanno T. Growth differentiation factor 15 in erythroid health and disease. *Curr Opin Hematol 2010*. [PMID 20182355](https://pubmed.ncbi.nlm.nih.gov/20182355/)
- Jelkmann W. Regulation of erythropoietin production. *J Physiol 2011*. [PMID 21078592](https://pubmed.ncbi.nlm.nih.gov/21078592/)
- Wenger RH. Erythropoietin. *Compr Physiol 2011*. [PMID 23733688](https://pubmed.ncbi.nlm.nih.gov/23733688/)
- Bunn HF. Erythropoietin. *Cold Spring Harb Perspect Med 2013*. [PMID 23457296](https://pubmed.ncbi.nlm.nih.gov/23457296/)
- Haase VH. Regulation of erythropoiesis by hypoxia-inducible factors. *Blood Rev 2013*. [PMID 23291219](https://pubmed.ncbi.nlm.nih.gov/23291219/)
- Orphanidou-Vlachou E. Extramedullary hemopoiesis. *Semin Ultrasound CT MR 2014*. [PMID 24929265](https://pubmed.ncbi.nlm.nih.gov/24929265/)
- Nemtsas P. Neurological complications of beta-thalassemia. *Ann Hematol 2015*. [PMID 25903043](https://pubmed.ncbi.nlm.nih.gov/25903043/)
- Piva E. Clinical utility of reticulocyte parameters. *Clin Lab Med 2015*. [PMID 25676377](https://pubmed.ncbi.nlm.nih.gov/25676377/)
- Wongtong N. Monocytosis is associated with hemolysis in sickle cell disease. *Hematology 2015*. [PMID 25875078](https://pubmed.ncbi.nlm.nih.gov/25875078/)
- Saliba AN. Morbidities in non-transfusion-dependent thalassemia. *Ann N Y Acad Sci 2016*. [PMID 27186941](https://pubmed.ncbi.nlm.nih.gov/27186941/)
- Yoon S. Comparable pharmacokinetics and pharmacodynamics of two epoetin alfa formulations Eporon(®) and Eprex(®) following a single subcutaneous administ. *Drug Des Devel Ther 2017*. [PMID 29138535](https://pubmed.ncbi.nlm.nih.gov/29138535/)
- Bolfa P. Thoracic and paraspinal extramedullary hematopoiesis in a cat with chronic non-regenerative anemia. *JFMS Open Rep 2018*. [PMID 30245843](https://pubmed.ncbi.nlm.nih.gov/30245843/)
- Cherry-Bukowiec JR. Hepcidin and Anemia in Surgical Critical Care: A Prospective Cohort Study. *Crit Care Med 2018*. [PMID 29517550](https://pubmed.ncbi.nlm.nih.gov/29517550/)
- El Nemer W. Ineffective erythropoiesis in sickle cell disease: new insights and future implications. *Curr Opin Hematol 2021*. [PMID 33631786](https://pubmed.ncbi.nlm.nih.gov/33631786/)
- Hasegawa S. Evaluation of recombinant human erythropoietin responsiveness by measuring erythrocyte creatine content in haemodialysis patients. *BMC Nephrol 2021*. [PMID 34895154](https://pubmed.ncbi.nlm.nih.gov/34895154/)
- Baird DC. Alpha- and Beta-thalassemia: Rapid Evidence Review. *Am Fam Physician 2022*. [PMID 35289581](https://pubmed.ncbi.nlm.nih.gov/35289581/)
- Chaichompoo P. The Roles of Mitophagy and Autophagy in Ineffective Erythropoiesis in β-Thalassemia. *Int J Mol Sci 2022*. [PMID 36142738](https://pubmed.ncbi.nlm.nih.gov/36142738/)
- Taylor CT. The effect of HIF on metabolism and immunity. *Nat Rev Nephrol 2022*. [PMID 35726016](https://pubmed.ncbi.nlm.nih.gov/35726016/)
- Guzzo I. Anemia after kidney transplantation. *Pediatr Nephrol 2023*. [PMID 36282330](https://pubmed.ncbi.nlm.nih.gov/36282330/)
- Liang R. Elevated CDKN1A (P21) mediates β-thalassemia erythroid apoptosis, but its loss does not improve β-thalassemic erythropoiesis. *Blood Adv 2023*. [PMID 37672319](https://pubmed.ncbi.nlm.nih.gov/37672319/)
- Nahm CH. Lipocalin-2, Soluble Transferrin Receptor, and Erythropoietin in Anemia During Mild Renal Dysfunction. *Int J Gen Med 2023*. [PMID 37637706](https://pubmed.ncbi.nlm.nih.gov/37637706/)
- Sato T. The roles of HIF-1α signaling in cardiovascular diseases. *J Cardiol 2023*. [PMID 36127212](https://pubmed.ncbi.nlm.nih.gov/36127212/)
- Bloise S. Parvovirus B19 infection in children: a comprehensive review of clinical manifestations and management. *Ital J Pediatr 2024*. [PMID 39696462](https://pubmed.ncbi.nlm.nih.gov/39696462/)
- Liu FF. Malaria and dyserythropoiesis: a mini review. *Front Cell Infect Microbiol 2025*. [PMID 40980010](https://pubmed.ncbi.nlm.nih.gov/40980010/)
- Peng Z. KDM4B modulates autocrine IL6 in erythroblasts to prevent ineffective erythropoiesis. *Leukemia 2025*. [PMID 40074853](https://pubmed.ncbi.nlm.nih.gov/40074853/)

### 7. Iron: erythroferrone, hepcidin, overload and chelation

*32 references*

- R'zik S. Serum soluble transferrin receptor concentration is an accurate estimate of the mass of tissue receptors. *Exp Hematol 2001*. [PMID 11378262](https://pubmed.ncbi.nlm.nih.gov/11378262/)
- Schumacher YO. Effects of exercise on soluble transferrin receptor and other variables of the iron status. *Br J Sports Med 2002*. [PMID 12055114](https://pubmed.ncbi.nlm.nih.gov/12055114/)
- Beguin Y. Soluble transferrin receptor for the evaluation of erythropoiesis and iron status. *Clin Chim Acta 2003*. [PMID 12589962](https://pubmed.ncbi.nlm.nih.gov/12589962/)
- Bérez V. Soluble transferrin receptor and mutations in hemochromatosis and transferrin genes in a general Catalan population. *Clin Chim Acta 2005*. [PMID 15698609](https://pubmed.ncbi.nlm.nih.gov/15698609/)
- Neufeld EJ. Oral chelators deferasirox and deferiprone for transfusional iron overload in thalassemia major: new data, new questions. *Blood 2006*. [PMID 16627763](https://pubmed.ncbi.nlm.nih.gov/16627763/)
- Cappellini MD. Oral iron chelators. *Annu Rev Med 2009*. [PMID 19630568](https://pubmed.ncbi.nlm.nih.gov/19630568/)
- Malyszko J. Hemojuvelin: the hepcidin story continues. *Kidney Blood Press Res 2009*. [PMID 19287179](https://pubmed.ncbi.nlm.nih.gov/19287179/)
- McLeod C. Deferasirox for the treatment of iron overload associated with regular blood transfusions (transfusional haemosiderosis) in patients suffering wi. *Health Technol Assess 2009*. [PMID 19068191](https://pubmed.ncbi.nlm.nih.gov/19068191/)
- Finberg KE. Down-regulation of Bmp/Smad signaling by Tmprss6 is required for maintenance of systemic iron homeostasis. *Blood 2010*. [PMID 20200349](https://pubmed.ncbi.nlm.nih.gov/20200349/)
- Lee DH. Neogenin inhibits HJV secretion and regulates BMP-induced hepcidin expression and iron homeostasis. *Blood 2010*. [PMID 20065295](https://pubmed.ncbi.nlm.nih.gov/20065295/)
- Lulla RR. Elevated soluble transferrin receptor levels reflect increased erythropoietic drive rather than iron deficiency in pediatric sickle cell disease. *Pediatr Blood Cancer 2010*. [PMID 20486179](https://pubmed.ncbi.nlm.nih.gov/20486179/)
- Meerpohl JJ. Deferasirox for managing transfusional iron overload in people with sickle cell disease. *Cochrane Database Syst Rev 2010*. [PMID 20687088](https://pubmed.ncbi.nlm.nih.gov/20687088/)
- Musallam KM. Iron chelation therapy for transfusional iron overload: a swift evolution. *Hemoglobin 2011*. [PMID 21910602](https://pubmed.ncbi.nlm.nih.gov/21910602/)
- Fuqua BK. Intestinal iron absorption. *J Trace Elem Med Biol 2012*. [PMID 22575541](https://pubmed.ncbi.nlm.nih.gov/22575541/)
- Patel M. Non transferrin bound iron: nature, manifestations and analytical approaches for estimation. *Indian J Clin Biochem 2012*. [PMID 24082455](https://pubmed.ncbi.nlm.nih.gov/24082455/)
- Gkouvatsos K. Iron-dependent regulation of hepcidin in Hjv-/- mice: evidence that hemojuvelin is dispensable for sensing body iron levels. *PLoS One 2014*. [PMID 24409331](https://pubmed.ncbi.nlm.nih.gov/24409331/)
- Lesjak M. Quercetin inhibits intestinal iron absorption and ferroportin transporter expression in vivo and in vitro. *PLoS One 2014*. [PMID 25058155](https://pubmed.ncbi.nlm.nih.gov/25058155/)
- Marsella M. Transfusional iron overload and iron chelation therapy in thalassemia major and sickle cell disease. *Hematol Oncol Clin North Am 2014*. [PMID 25064709](https://pubmed.ncbi.nlm.nih.gov/25064709/)
- Przybyszewska J. The role of hepcidin, ferroportin, HCP1, and DMT1 protein in iron absorption in the human digestive tract. *Prz Gastroenterol 2014*. [PMID 25276251](https://pubmed.ncbi.nlm.nih.gov/25276251/)
- Kent P. Hfe and Hjv exhibit overlapping functions for iron signaling to hepcidin. *J Mol Med (Berl) 2015*. [PMID 25609138](https://pubmed.ncbi.nlm.nih.gov/25609138/)
- Frazer DM. Ferroportin Is Essential for Iron Absorption During Suckling, But Is Hyporesponsive to the Regulatory Hormone Hepcidin. *Cell Mol Gastroenterol Hepatol 2017*. [PMID 28462381](https://pubmed.ncbi.nlm.nih.gov/28462381/)
- Sato M. Increased Duodenal Iron Absorption through Upregulation of Ferroportin 1 due to the Decrement in Serum Hepcidin in Patients with Chronic Hepatiti. *Can J Gastroenterol Hepatol 2018*. [PMID 30186818](https://pubmed.ncbi.nlm.nih.gov/30186818/)
- Knutson MD. Non-transferrin-bound iron transporters. *Free Radic Biol Med 2019*. [PMID 30316781](https://pubmed.ncbi.nlm.nih.gov/30316781/)
- Silvestri L. Hepcidin and the BMP-SMAD pathway: An unexpected liaison. *Vitam Horm 2019*. [PMID 30798817](https://pubmed.ncbi.nlm.nih.gov/30798817/)
- Xiao X. Bone morphogenic proteins in iron homeostasis. *Bone 2020*. [PMID 32585319](https://pubmed.ncbi.nlm.nih.gov/32585319/)
- Sugiura T. Analytical evaluation of serum non-transferrin-bound iron and its relationships with oxidative stress and cardiac load in the general population. *Medicine (Baltimore) 2021*. [PMID 33607814](https://pubmed.ncbi.nlm.nih.gov/33607814/)
- Silva AMN. The (Bio)Chemistry of Non-Transferrin-Bound Iron. *Molecules 2022*. [PMID 35335148](https://pubmed.ncbi.nlm.nih.gov/35335148/)
- —. . *2023*. [PMID 38502749](https://pubmed.ncbi.nlm.nih.gov/38502749/)
- Bruzzese A. Iron chelation therapy. *Eur J Haematol 2023*. [PMID 36708354](https://pubmed.ncbi.nlm.nih.gov/36708354/)
- Girelli D. Diagnosis and management of hereditary hemochromatosis: lifestyle modification, phlebotomy, and blood donation. *Hematology Am Soc Hematol Educ Program 2024*. [PMID 39644049](https://pubmed.ncbi.nlm.nih.gov/39644049/)
- Duca L. The Relationship Between Non-Transferrin-Bound Iron (NTBI), Labile Plasma Iron (LPI), and Iron Toxicity. *Int J Mol Sci 2025*. [PMID 40650208](https://pubmed.ncbi.nlm.nih.gov/40650208/)
- Enko D. Physiology of Iron Metabolism. *Clin Lab 2025*. [PMID 40066539](https://pubmed.ncbi.nlm.nih.gov/40066539/)

### 8. Bilirubin, UGT1A1, gallstones and the neonatal presentation

*31 references*

- Lundh B. Heme catabolism, carbon monoxide production and red cell survival in anemia. *Acta Med Scand 1975*. [PMID 1124665](https://pubmed.ncbi.nlm.nih.gov/1124665/)
- Soloway RD. Pigment gallstones. *Gastroenterology 1977*. [PMID 318581](https://pubmed.ncbi.nlm.nih.gov/318581/)
- Werner B. Endogenous carbon monoxide production after bicycle exercise in healthy subjects and in patients with hereditary spherocytosis. *Scand J Clin Lab Invest 1980*. [PMID 7414250](https://pubmed.ncbi.nlm.nih.gov/7414250/)
- Rothuizen J. Bilirubin metabolism in canine hepatobiliary and haemolytic disease. *Vet Q 1987*. [PMID 3672859](https://pubmed.ncbi.nlm.nih.gov/3672859/)
- Trotman BW. Pigment gallstone disease. *Gastroenterol Clin North Am 1991*. [PMID 2022417](https://pubmed.ncbi.nlm.nih.gov/2022417/)
- Iyer L. Phenotype-genotype correlation of in vitro SN-38 (active metabolite of irinotecan) and bilirubin glucuronidation in human liver tissue with UGT1A. *Clin Pharmacol Ther 1999*. [PMID 10340924](https://pubmed.ncbi.nlm.nih.gov/10340924/)
- Maruo Y. [UDP-glucuronosyltransferase]. *Nihon Eiseigaku Zasshi 2002*. [PMID 11868392](https://pubmed.ncbi.nlm.nih.gov/11868392/)
- American Academy of Pediatrics Subcommittee on Hyperbilirubinemia. Management of hyperbilirubinemia in the newborn infant 35 or more weeks of gestation. *Pediatrics 2004*. [PMID 15231951](https://pubmed.ncbi.nlm.nih.gov/15231951/)
- Amin SB. Clinical assessment of bilirubin-induced neurotoxicity in premature infants. *Semin Perinatol 2004*. [PMID 15686265](https://pubmed.ncbi.nlm.nih.gov/15686265/)
- Frank JE. Diagnosis and management of G6PD deficiency. *Am Fam Physician 2005*. [PMID 16225031](https://pubmed.ncbi.nlm.nih.gov/16225031/)
- Chang JL. UGT1A1 polymorphism is associated with serum bilirubin concentrations in a randomized, controlled, fruit and vegetable feeding trial. *J Nutr 2007*. [PMID 17374650](https://pubmed.ncbi.nlm.nih.gov/17374650/)
- Watson RL. Hyperbilirubinemia. *Crit Care Nurs Clin North Am 2009*. [PMID 19237047](https://pubmed.ncbi.nlm.nih.gov/19237047/)
- Cariati A. Limits and perspective of oral therapy with statins and aspirin for the prevention of symptomatic cholesterol gallstone disease. *Expert Opin Pharmacother 2012*. [PMID 22607008](https://pubmed.ncbi.nlm.nih.gov/22607008/)
- Gil J. Gilbert syndrome: the UGT1A1*28 promoter polymorphism as a biomarker of multifactorial diseases and drug metabolism. *Biomark Med 2012*. [PMID 22448797](https://pubmed.ncbi.nlm.nih.gov/22448797/)
- Stinton LM. Epidemiology of gallbladder disease: cholelithiasis and cancer. *Gut Liver 2012*. [PMID 22570746](https://pubmed.ncbi.nlm.nih.gov/22570746/)
- Bhutani VK. Bilirubin neurotoxicity in preterm infants: risk and prevention. *J Clin Neonatol 2013*. [PMID 24049745](https://pubmed.ncbi.nlm.nih.gov/24049745/)
- Cariati A. Blackberry pigment (whitlockite) gallstones in uremic patient. *Clin Res Hepatol Gastroenterol 2013*. [PMID 22959097](https://pubmed.ncbi.nlm.nih.gov/22959097/)
- Hulzebos CV. Bilirubin-albumin binding, bilirubin/albumin ratios, and free bilirubin levels: where do we stand?. *Semin Perinatol 2014*. [PMID 25304058](https://pubmed.ncbi.nlm.nih.gov/25304058/)
- Morioka I. Disorders of bilirubin binding to albumin and bilirubin-induced neurologic dysfunction. *Semin Fetal Neonatal Med 2015*. [PMID 25432488](https://pubmed.ncbi.nlm.nih.gov/25432488/)
- Ahlfors CE. The Bilirubin Binding Panel: A Henderson-Hasselbalch Approach to Neonatal Hyperbilirubinemia. *Pediatrics 2016*. [PMID 27609825](https://pubmed.ncbi.nlm.nih.gov/27609825/)
- Amin SB. Bilirubin Binding Capacity in the Preterm Neonate. *Clin Perinatol 2016*. [PMID 27235205](https://pubmed.ncbi.nlm.nih.gov/27235205/)
- Mitra S. Neonatal jaundice: aetiology, diagnosis and treatment. *Br J Hosp Med (Lond) 2017*. [PMID 29240507](https://pubmed.ncbi.nlm.nih.gov/29240507/)
- Kaplan M. Hemolysis and Glucose-6-Phosphate Dehydrogenase Deficiency-Related Neonatal Hyperbilirubinemia. *Neonatology 2018*. [PMID 29940590](https://pubmed.ncbi.nlm.nih.gov/29940590/)
- Chávez-Peña T. Prevalence of UGT1A1 (TA)n promoter polymorphism in Panamanians neonates with G6PD deficiency. *J Genet 2020*. [PMID 33622990](https://pubmed.ncbi.nlm.nih.gov/33622990/)
- Stevenson DK. Increased Carbon Monoxide Washout Rates in Newborn Infants. *Neonatology 2020*. [PMID 31634890](https://pubmed.ncbi.nlm.nih.gov/31634890/)
- Par EJ. Neonatal Hyperbilirubinemia: Evaluation and Treatment. *Am Fam Physician 2023*. [PMID 37192079](https://pubmed.ncbi.nlm.nih.gov/37192079/)
- Chastain AP. Managing neonatal hyperbilirubinemia: An updated guideline. *JAAPA 2024*. [PMID 39259272](https://pubmed.ncbi.nlm.nih.gov/39259272/)
- Mobaraki S. Cabozantinib Induces Isolated Hyperbilirubinemia in Renal Cell Carcinoma Patients carrying the UGT1A1*28 Polymorphism. *Clin Genitourin Cancer 2024*. [PMID 39155162](https://pubmed.ncbi.nlm.nih.gov/39155162/)
- Jagroo J. Confronting Cholelithiasis: A Case Series of Patients With Sickle Cell Disease and Gallstones. *Cureus 2025*. [PMID 40225473](https://pubmed.ncbi.nlm.nih.gov/40225473/)
- Wickremasinghe AC. Neonatal Hyperbilirubinemia. *Pediatr Clin North Am 2025*. [PMID 40619190](https://pubmed.ncbi.nlm.nih.gov/40619190/)
- Yang G. End-tidal carbon monoxide levels as a point-of-care biomarker for the early detection of hemolytic disease in Chinese neonates. *Clin Chim Acta 2025*. [PMID 40383361](https://pubmed.ncbi.nlm.nih.gov/40383361/)

### 9. Mitapivat and other PK-R activators: trials and pharmacology

*13 references*

- Vanderschueren D. Skeletal effects of estrogen deficiency as induced by an aromatase inhibitor in an aged male rat model. *Bone 2000*. [PMID 11062346](https://pubmed.ncbi.nlm.nih.gov/11062346/)
- Nathan L. Testosterone inhibits early atherogenesis by conversion to estradiol: critical role of aromatase. *Proc Natl Acad Sci U S A 2001*. [PMID 11248122](https://pubmed.ncbi.nlm.nih.gov/11248122/)
- Kanakis GA. EAA clinical practice guidelines-gynecomastia evaluation and management. *Andrology 2019*. [PMID 31099174](https://pubmed.ncbi.nlm.nih.gov/31099174/)
- Yang C. Clinical application of aromatase inhibitors to treat male infertility. *Hum Reprod Update 2021*. [PMID 34871401](https://pubmed.ncbi.nlm.nih.gov/34871401/)
- Guo B. Efficacy and safety of letrozole or anastrozole in the treatment of male infertility with low testosterone-estradiol ratio: A meta-analysis and s. *Andrology 2022*. [PMID 35438843](https://pubmed.ncbi.nlm.nih.gov/35438843/)
- Idris IM. Epidemiology and treatment of priapism in sickle cell disease. *Hematology Am Soc Hematol Educ Program 2022*. [PMID 36485155](https://pubmed.ncbi.nlm.nih.gov/36485155/)
- van Dijk MJ. One-year safety and efficacy of mitapivat in sickle cell disease: follow-up results of a phase 2, open-label study. *Blood Adv 2023*. [PMID 37934880](https://pubmed.ncbi.nlm.nih.gov/37934880/)
- van Beers EJ. Mitapivat improves ineffective erythropoiesis and iron overload in adult patients with pyruvate kinase deficiency. *Blood Adv 2024*. [PMID 38330179](https://pubmed.ncbi.nlm.nih.gov/38330179/)
- Algeri M. The ENERGIZE trial: Is mitapivat ready to take center stage in NTDT management?. *Med 2025*. [PMID 40945502](https://pubmed.ncbi.nlm.nih.gov/40945502/)
- Conrey A. Long-term mitapivat treatment is safe and efficacious in patients with sickle cell disease. *Blood Red Cells Iron 2025*. [PMID 41809202](https://pubmed.ncbi.nlm.nih.gov/41809202/)
- Idowu M. Safety and efficacy of mitapivat in sickle cell disease (RISE UP): results from the phase 2 portion of a global, double-blind, randomised, placeb. *Lancet Haematol 2025*. [PMID 39644907](https://pubmed.ncbi.nlm.nih.gov/39644907/)
- Kuo KHM. Long-term efficacy and safety of mitapivat in non-transfusion-dependent α- or β-thalassaemia: An open-label phase 2 study. *Br J Haematol 2025*. [PMID 40394935](https://pubmed.ncbi.nlm.nih.gov/40394935/)
- Taher AT. Mitapivat in adults with non-transfusion-dependent α-thalassaemia or β-thalassaemia (ENERGIZE): a phase 3, international, randomised, double-blin. *Lancet 2025*. [PMID 40544857](https://pubmed.ncbi.nlm.nih.gov/40544857/)

### 10. Curative approaches: gene therapy, editing, transplantation

*8 references*

- Trobridge GD. Stem cell selection in vivo using foamy vectors cures canine pyruvate kinase deficiency. *PLoS One 2012*. [PMID 23028826](https://pubmed.ncbi.nlm.nih.gov/23028826/)
- van Straaten S. Worldwide study of hematopoietic allogeneic stem cell transplantation in pyruvate kinase deficiency. *Haematologica 2018*. [PMID 29242305](https://pubmed.ncbi.nlm.nih.gov/29242305/)
- Quintana-Bustamante O. Gene editing of PKLR gene in human hematopoietic progenitors through 5' and 3' UTR modified TALEN mRNA. *PLoS One 2019*. [PMID 31618280](https://pubmed.ncbi.nlm.nih.gov/31618280/)
- Fañanas-Baquero S. Clinically relevant gene editing in hematopoietic stem cells for the treatment of pyruvate kinase deficiency. *Mol Ther Methods Clin Dev 2021*. [PMID 34485608](https://pubmed.ncbi.nlm.nih.gov/34485608/)
- Ma ZY. Allogeneic hematopoietic stem cell transplantation in a 3-year-old boy with congenital pyruvate kinase deficiency: A case report. *World J Clin Cases 2021*. [PMID 33969077](https://pubmed.ncbi.nlm.nih.gov/33969077/)
- Morado M. Consensus document for the diagnosis and treatment of pyruvate kinase deficiency. *Med Clin (Barc) 2021*. [PMID 33431182](https://pubmed.ncbi.nlm.nih.gov/33431182/)
- Fañanas-Baquero S. Specific correction of pyruvate kinase deficiency-causing point mutations by CRISPR/Cas9 and single-stranded oligodeoxynucleotides. *Front Genome Ed 2023*. [PMID 37188156](https://pubmed.ncbi.nlm.nih.gov/37188156/)
- Pang Y. Case report: Modified transplantation for pediatric patients with pyruvate kinase deficiency. *Front Immunol 2024*. [PMID 39635530](https://pubmed.ncbi.nlm.nih.gov/39635530/)

### 11. Systemic complications and patient-reported outcomes

*27 references*

- Ghidini A. Hepatosplenomegaly as the only prenatal finding in a fetus with pyruvate kinase deficiency anemia. *Am J Perinatol 1991*. [PMID 1987968](https://pubmed.ncbi.nlm.nih.gov/1987968/)
- Filosa A. Longitudinal monitoring of bone mineral density in thalassemic patients. Genetic structure and osteoporosis. *Acta Paediatr 1997*. [PMID 9174216](https://pubmed.ncbi.nlm.nih.gov/9174216/)
- Ghidini A. Severe pyruvate kinase deficiency anemia. A case report. *J Reprod Med 1998*. [PMID 9749428](https://pubmed.ncbi.nlm.nih.gov/9749428/)
- Almeida A. Bone involvement in sickle cell disease. *Br J Haematol 2005*. [PMID 15877730](https://pubmed.ncbi.nlm.nih.gov/15877730/)
- Lin EE. Hemolytic anemia-associated pulmonary hypertension in sickle cell disease. *Curr Hematol Rep 2005*. [PMID 15720960](https://pubmed.ncbi.nlm.nih.gov/15720960/)
- Kato GJ. Lactate dehydrogenase as a biomarker of hemolysis-associated nitric oxide resistance, priapism, leg ulceration, pulmonary hypertension, and death. *Blood 2006*. [PMID 16291595](https://pubmed.ncbi.nlm.nih.gov/16291595/)
- Machado RF. Sickle cell anemia-associated pulmonary arterial hypertension. *J Bras Pneumol 2007*. [PMID 18026658](https://pubmed.ncbi.nlm.nih.gov/18026658/)
- Vermylen C. What is new in iron overload?. *Eur J Pediatr 2008*. [PMID 17899187](https://pubmed.ncbi.nlm.nih.gov/17899187/)
- Newaskar M. Asthma in sickle cell disease. *ScientificWorldJournal 2011*. [PMID 21623460](https://pubmed.ncbi.nlm.nih.gov/21623460/)
- Mohamed S. A case of severe pyruvate kinase deficiency in a primigravida: successful outcome. *Obstet Med 2013*. [PMID 27757165](https://pubmed.ncbi.nlm.nih.gov/27757165/)
- Adams-Graves P. Bone mineral density patterns in vitamin D deficient African American men with sickle cell disease. *Am J Med Sci 2014*. [PMID 23538935](https://pubmed.ncbi.nlm.nih.gov/23538935/)
- Irwig MS. Bone health in hypogonadal men. *Curr Opin Urol 2014*. [PMID 25144148](https://pubmed.ncbi.nlm.nih.gov/25144148/)
- Barcellini W. Clinical Applications of Hemolytic Markers in the Differential Diagnosis and Management of Hemolytic Anemia. *Dis Markers 2015*. [PMID 26819490](https://pubmed.ncbi.nlm.nih.gov/26819490/)
- Potoka KP. Vasculopathy and pulmonary hypertension in sickle cell disease. *Am J Physiol Lung Cell Mol Physiol 2015*. [PMID 25398989](https://pubmed.ncbi.nlm.nih.gov/25398989/)
- Gordeuk VR. Pathophysiology and treatment of pulmonary hypertension in sickle cell disease. *Blood 2016*. [PMID 26758918](https://pubmed.ncbi.nlm.nih.gov/26758918/)
- Lal A. Assessment and treatment of pain in thalassemia. *Ann N Y Acad Sci 2016*. [PMID 27124110](https://pubmed.ncbi.nlm.nih.gov/27124110/)
- Kato GJ. Intravascular hemolysis and the pathophysiology of sickle cell disease. *J Clin Invest 2017*. [PMID 28248201](https://pubmed.ncbi.nlm.nih.gov/28248201/)
- Abid S. New Nitric Oxide Donor NCX 1443: Therapeutic Effects on Pulmonary Hypertension in the SAD Mouse Model of Sickle Cell Disease. *J Cardiovasc Pharmacol 2018*. [PMID 29438213](https://pubmed.ncbi.nlm.nih.gov/29438213/)
- Bernardi MH. Hemoadsorption does not Have Influence on Hemolysis During Cardiopulmonary Bypass. *ASAIO J 2019*. [PMID 30325849](https://pubmed.ncbi.nlm.nih.gov/30325849/)
- Niccoli Asabella A. Sickle cell diseases: What can nuclear medicine offer?. *Hell J Nucl Med 2019*. [PMID 30843001](https://pubmed.ncbi.nlm.nih.gov/30843001/)
- Etemad K. Quality of Life and Related Factors in β-Thalassemia Patients. *Hemoglobin 2021*. [PMID 34409903](https://pubmed.ncbi.nlm.nih.gov/34409903/)
- Jang JH. Iptacopan monotherapy in patients with paroxysmal nocturnal hemoglobinuria: a 2-cohort open-label proof-of-concept study. *Blood Adv 2022*. [PMID 35561315](https://pubmed.ncbi.nlm.nih.gov/35561315/)
- Al-Samkari H. Mitapivat for Acquired Pyruvate Kinase Deficiency. *Pediatr Blood Cancer 2025*. [PMID 39538432](https://pubmed.ncbi.nlm.nih.gov/39538432/)
- Jain P. Spectrum of Ophthalmic Manifestations in Patients With Transfusion-Dependent Thalassemia. *Cureus 2025*. [PMID 40370899](https://pubmed.ncbi.nlm.nih.gov/40370899/)
- Kariki O. Biochemical evidence and clinical significance of hemolysis following catheter ablation with a lattice-tip PFA catheter. *J Interv Card Electrophysiol 2025*. [PMID 41129021](https://pubmed.ncbi.nlm.nih.gov/41129021/)
- Koh TWY. Intravascular haemolysis following pulsed field ablation pulmonary vein isolation for atrial fibrillation: a systematic review. *Open Heart 2025*. [PMID 41213824](https://pubmed.ncbi.nlm.nih.gov/41213824/)
- Ahmed K. Adult survivors of sickle cell disease, transfusion-dependent beta-thalassaemia and childhood acute leukaemia in England: protocol for a mixed me. *BMJ Open 2026*. [PMID 41628927](https://pubmed.ncbi.nlm.nih.gov/41628927/)

### 12. Related red cell disorders, transfusion biology and redox

*18 references*

- Siems WG. Erythrocyte free radical and energy metabolism. *Clin Nephrol 2000*. [PMID 10746800](https://pubmed.ncbi.nlm.nih.gov/10746800/)
- Cappellini MD. Glucose-6-phosphate dehydrogenase deficiency. *Lancet 2008*. [PMID 18177777](https://pubmed.ncbi.nlm.nih.gov/18177777/)
- Hess JR. Red cell storage. *J Proteomics 2010*. [PMID 19914410](https://pubmed.ncbi.nlm.nih.gov/19914410/)
- van Zwieten R. Inborn defects in the antioxidant systems of human red blood cells. *Free Radic Biol Med 2014*. [PMID 24316370](https://pubmed.ncbi.nlm.nih.gov/24316370/)
- D'Alessandro A. Red blood cell storage in additive solution-7 preserves energy and redox metabolism: a metabolomics approach. *Transfusion 2015*. [PMID 26271632](https://pubmed.ncbi.nlm.nih.gov/26271632/)
- Al-Abdi SY. Decreased Glutathione S-transferase Level and Neonatal Hyperbilirubinemia Associated with Glucose-6-phosphate Dehydrogenase Deficiency: A Perspec. *Am J Perinatol 2017*. [PMID 27464020](https://pubmed.ncbi.nlm.nih.gov/27464020/)
- Belfield KD. Review and drug therapy implications of glucose-6-phosphate dehydrogenase deficiency. *Am J Health Syst Pharm 2018*. [PMID 29305344](https://pubmed.ncbi.nlm.nih.gov/29305344/)
- Gehrke S. Metabolomics evaluation of early-storage red blood cell rejuvenation at 4°C and 37°C. *Transfusion 2018*. [PMID 29687892](https://pubmed.ncbi.nlm.nih.gov/29687892/)
- Ghashghaeinia M. Proliferating tumor cells mimick glucose metabolism of mature human erythrocytes. *Cell Cycle 2019*. [PMID 31154896](https://pubmed.ncbi.nlm.nih.gov/31154896/)
- La Vieille S. Dietary restrictions for people with glucose-6-phosphate dehydrogenase deficiency. *Nutr Rev 2019*. [PMID 30380124](https://pubmed.ncbi.nlm.nih.gov/30380124/)
- Li R. Exploring the role of glucose‑6‑phosphate dehydrogenase in cancer (Review). *Oncol Rep 2020*. [PMID 33125150](https://pubmed.ncbi.nlm.nih.gov/33125150/)
- Luzzatto L. Glucose-6-phosphate dehydrogenase deficiency. *Blood 2020*. [PMID 32702756](https://pubmed.ncbi.nlm.nih.gov/32702756/)
- Lee HY. Glucose-6-Phosphate Dehydrogenase Deficiency and Neonatal Hyperbilirubinemia: Insights on Pathophysiology, Diagnosis, and Gene Variants in Diseas. *Front Pediatr 2022*. [PMID 35685917](https://pubmed.ncbi.nlm.nih.gov/35685917/)
- D'Alessandro A. Red Blood Cell Storage: From Genome to Exposome Towards Personalized Transfusion Medicine. *Transfus Med Rev 2023*. [PMID 37574398](https://pubmed.ncbi.nlm.nih.gov/37574398/)
- D'Alessandro A. Red Blood Cell Metabolism In Vivo and In Vitro. *Metabolites 2023*. [PMID 37512500](https://pubmed.ncbi.nlm.nih.gov/37512500/)
- Nemkov T. Supercooled storage of red blood cells slows down the metabolic storage lesion. *Sci Rep 2025*. [PMID 41044137](https://pubmed.ncbi.nlm.nih.gov/41044137/)
- Wang X. Obesity-aggravated erythrocyte injury in ischaemia-reperfusion: Interlinked oxidative stress, metabolic reprogramming, and cytoskeletal destabili. *Life Sci 2025*. [PMID 40975374](https://pubmed.ncbi.nlm.nih.gov/40975374/)
- Mak GK. Glucose-6-Phosphate Dehydrogenase Deficiency. *2026*. [PMID 29262208](https://pubmed.ncbi.nlm.nih.gov/29262208/)

### 13. QSP and pharmacometric methodology

*13 references*

- Betts A. A Translational Quantitative Systems Pharmacology Model for CD3 Bispecific Molecules: Application to Quantify T Cell-Mediated Tumor Cell Killing. *AAPS J 2019*. [PMID 31119428](https://pubmed.ncbi.nlm.nih.gov/31119428/)
- El-Khateeb E. Quantitative mass spectrometry-based proteomics in the era of model-informed drug development: Applications in translational pharmacology and rec. *Pharmacol Ther 2019*. [PMID 31376433](https://pubmed.ncbi.nlm.nih.gov/31376433/)
- Elmokadem A. Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial. *CPT Pharmacometrics Syst Pharmacol 2019*. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
- Nguyen LM. A quantitative systems pharmacology model of hyporesponsiveness to erythropoietin in rats. *J Pharmacokinet Pharmacodyn 2021*. [PMID 34100188](https://pubmed.ncbi.nlm.nih.gov/34100188/)
- Ball K. Strategies for clinical dose optimization of T cell-engaging therapies in oncology. *MAbs 2023*. [PMID 36823042](https://pubmed.ncbi.nlm.nih.gov/36823042/)
- Ji Y. Quantitative systems pharmacology model of GITR-mediated T cell dynamics in tumor microenvironment. *CPT Pharmacometrics Syst Pharmacol 2023*. [PMID 36710369](https://pubmed.ncbi.nlm.nih.gov/36710369/)
- Li X. Combining network pharmacology, molecular docking, molecular dynamics simulation, and experimental verification to examine the efficacy and immun. *Front Immunol 2023*. [PMID 37575231](https://pubmed.ncbi.nlm.nih.gov/37575231/)
- Arsène S. In Silico Clinical Trials: Is It Possible?. *Methods Mol Biol 2024*. [PMID 37702936](https://pubmed.ncbi.nlm.nih.gov/37702936/)
- GBD 2021 US Burden of Disease and Forecasting Collaborators. Burden of disease scenarios by state in the USA, 2022-50: a forecasting analysis for the Global Burden of Disease Study 2021. *Lancet 2024*. [PMID 39645377](https://pubmed.ncbi.nlm.nih.gov/39645377/)
- Luo Y. Neoadjuvant PARPi or chemotherapy in ovarian cancer informs targeting effector Treg cells for homologous-recombination-deficient tumors. *Cell 2024*. [PMID 38971151](https://pubmed.ncbi.nlm.nih.gov/38971151/)
- GBD 2023 Disease and Injury and Risk Factor Collaborators. Burden of 375 diseases and injuries, risk-attributable burden of 88 risk factors, and healthy life expectancy in 204 countries and territories, i. *Lancet 2025*. [PMID 41092926](https://pubmed.ncbi.nlm.nih.gov/41092926/)
- GBD 2023 Lower Respiratory Infections and Antimicrobial Resistance Collaborators. Global burden of lower respiratory infections and aetiologies, 1990-2023: a systematic analysis for the Global Burden of Disease Study 2023. *Lancet Infect Dis 2026*. [PMID 41412141](https://pubmed.ncbi.nlm.nih.gov/41412141/)
- GBD 2023 Meningitis & Antimicrobial Resistance Collaborators. Global, regional, and national burden of meningitis, its risk factors, and aetiologies, 1990-2023: a systematic analysis for the Global Burden of. *Lancet Neurol 2026*. [PMID 41911930](https://pubmed.ncbi.nlm.nih.gov/41911930/)
