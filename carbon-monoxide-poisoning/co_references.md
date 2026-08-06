# Carbon Monoxide Poisoning — annotated reference list

Organised by the model component each reference supports.

**Provenance.** Every entry was resolved against PubMed programmatically through the
NCBI E-utilities API, and the title, journal, year, first author and PMID are
reproduced verbatim from the returned record. No citation here was written from
memory. An earlier attempt to look up remembered titles was discarded: of 98
remembered citations, 7 returned nothing at all and dozens more silently matched
unrelated papers — one "Carbon monoxide and the brain" query returned a review of
nutrients in dog and cat cognition. That is precisely the failure mode that puts
confident-looking but wrong citations into a bibliography, so the list was rebuilt
from topic searches plus explicitly verified landmark lookups.

Verification also corrected a specific error the author would otherwise have
committed: Weaver's carboxyhaemoglobin half-life study — the source of the 320 and
74 minute anchors this model is calibrated against — is in *Chest* (2000), not the
journal it was misremembered as.

Landmark papers are marked ★ and are the ones the model's quantitative anchors and
structural commitments actually rest on.

---

## A. Coburn-Forster-Kane kinetics and carboxyhaemoglobin half-life

> These fix the pulmonary-exchange block. The CFK constants used in the model (alveolar ventilation 4.2 L/min, DLco 25 mL/min/mmHg, Haldane ratio M = 245, blood volume 5.5 L) are taken from resting physiology rather than fitted, and reproduce the measured half-lives to within 2-3% at normobaric pressure: 313 min against Weaver's observed 320 on room air, and 72 against 74 on a tight 100% oxygen mask. They FAIL at hyperbaric pressure, and the model reports that failure rather than fitting around it.

- ★ Weaver LK et al. (2002). *Hyperbaric oxygen for acute carbon monoxide poisoning.* N Engl J Med. [PMID 12362006](https://pubmed.ncbi.nlm.nih.gov/12362006/)
- ★ Weaver LK et al. (2000). *Carboxyhemoglobin half-life in carbon monoxide-poisoned patients treated with 100% oxygen at atmospheric pressure.* Chest. [PMID 10713010](https://pubmed.ncbi.nlm.nih.gov/10713010/)
- ★ Coburn RF et al. (1965). *Considerations of the physiological variables that determine the blood carboxyhemoglobin concentration in man.* J Clin Invest. [PMID 5845666](https://pubmed.ncbi.nlm.nih.gov/5845666/)
- ★ Peterson JE et al. (1970). *Absorption and elimination of carbon monoxide by inactive young men.* Arch Environ Health. [PMID 5430002](https://pubmed.ncbi.nlm.nih.gov/5430002/)
- Coburn RF (2013). *Carbon monoxide uptake and excretion: testing assumptions made in deriving the Coburn-Forster-Kane equation.* Respir Physiol Neurobiol. [PMID 23602912](https://pubmed.ncbi.nlm.nih.gov/23602912/)
- Kuo YW et al. (2020). *Endogenous Carbon Monoxide Production and Diffusing Capacity of the Lung for Carbon Monoxide in Sepsis-Induced Acute Respiratory Distress Syndrome.* Crit Care Explor. [PMID 33283194](https://pubmed.ncbi.nlm.nih.gov/33283194/)
- Bleecker ML (2015). *Carbon monoxide intoxication.* Handb Clin Neurol. [PMID 26563790](https://pubmed.ncbi.nlm.nih.gov/26563790/)
- Culnan DM et al. (2018). *Carbon Monoxide and Cyanide Poisoning in the Burned Pregnant Patient: An Indication for Hyperbaric Oxygen Therapy.* Ann Plast Surg. [PMID 29461288](https://pubmed.ncbi.nlm.nih.gov/29461288/)
- Ozturan IU et al. (2019). *Determination of carboxyhemoglobin half-life in patients with carbon monoxide toxicity treated with high flow nasal cannula oxygen therapy.* Clin Toxicol (Phila). [PMID 30689450](https://pubmed.ncbi.nlm.nih.gov/30689450/)
- Bernard TE et al. (1981). *Modeling carbon monoxide uptake during work.* Am Ind Hyg Assoc J. [PMID 7223644](https://pubmed.ncbi.nlm.nih.gov/7223644/)
- Guan L et al. (2011). *Dynamic changes of heme oxygenase-1 in the hippocampus of rats after acute carbon monoxide poisoning.* Arch Environ Contam Toxicol. [PMID 20422170](https://pubmed.ncbi.nlm.nih.gov/20422170/)
- Sakamoto A et al. (2005). *Does carboxy-hemoglobin serve as a stress-induced inflammatory marker reflecting surgical insults?.* J Nippon Med Sch. [PMID 15834204](https://pubmed.ncbi.nlm.nih.gov/15834204/)

## B. Haldane relation, oxygen content and the allosteric left shift

> The allosteric left shift is the second, less appreciated hit. At matched arterial oxygen content, carbon monoxide leaves tissue PO2 about 9 mmHg lower than the anaemia it mimics, because the residual haemoglobin will not release what it holds.

- Halebian P et al. (1986). *Whole body oxygen utilization during acute carbon monoxide poisoning and isocapneic nitrogen hypoxia.* J Trauma. [PMID 3944835](https://pubmed.ncbi.nlm.nih.gov/3944835/)
- Buehler JH et al. (1975). *Lactic acidosis from carboxyhemoglobinemia after smoke inhalation.* Ann Intern Med. [PMID 237451](https://pubmed.ncbi.nlm.nih.gov/237451/)
- Koehler RC et al. (1984). *Comparison of cerebrovascular response to hypoxic and carbon monoxide hypoxia in newborn and adult sheep.* J Cereb Blood Flow Metab. [PMID 6420426](https://pubmed.ncbi.nlm.nih.gov/6420426/)
- Zock JP (1990). *Carbon monoxide binding in a model of hemoglobin differs between the T and the R conformation.* Adv Exp Med Biol. [PMID 2096625](https://pubmed.ncbi.nlm.nih.gov/2096625/)
- Rhodes CE et al. (2026). *Physiology, Oxygen Transport.* —. [PMID 30855920](https://pubmed.ncbi.nlm.nih.gov/30855920/)
- Fischbach A et al. (2022). *Hyperbaric phototherapy augments blood carbon monoxide removal.* Lasers Surg Med. [PMID 34658052](https://pubmed.ncbi.nlm.nih.gov/34658052/)

## C. Cytochrome c oxidase, myoglobin and the intracellular pool

> The intracellular pool. Carbon monoxide binds only the REDUCED a3 haem, so oxygen acts as a competitive protector and the tissue occupancy behaves nothing like COHb: it loads more readily when tissue PO2 is low and unloads over hours rather than minutes. This is the structural basis for the model's central claim of two clocks.

- Cheng Y et al. (2015). *Carbon monoxide modulates cytochrome oxidase activity and oxidative stress in the developing murine brain during isoflurane exposure.* Free Radic Biol Med. [PMID 26032170](https://pubmed.ncbi.nlm.nih.gov/26032170/)
- Lee HM et al. (2010). *Differential inhibition of mitochondrial respiratory complexes by inhalation of combustion smoke and carbon monoxide, in vivo, in the rat brain.* Inhal Toxicol. [PMID 20429857](https://pubmed.ncbi.nlm.nih.gov/20429857/)
- Savolainen H et al. (1980). *Biochemical effects of carbon monoxide poisoning in rat brain with special reference to blood carboxyhemoglobin and cerebral cytochrome oxidase activity.* Neurosci Lett. [PMID 6302602](https://pubmed.ncbi.nlm.nih.gov/6302602/)
- Wattel F et al. (1996). *Carbon monoxide poisoning.* Presse Med. [PMID 8958870](https://pubmed.ncbi.nlm.nih.gov/8958870/)
- McGrath JJ (2000). *Biological plausibility for carbon monoxide as a copollutant in PM epidemiologic studies.* Inhal Toxicol. [PMID 12881888](https://pubmed.ncbi.nlm.nih.gov/12881888/)
- De Sanctis G et al. (1986). *Mini-myoglobin: preparation and reaction with oxygen and carbon monoxide.* J Mol Biol. [PMID 3712445](https://pubmed.ncbi.nlm.nih.gov/3712445/)
- Mahan VL (2020). *Cardiac function dependence on carbon monoxide.* Med Gas Res. [PMID 32189668](https://pubmed.ncbi.nlm.nih.gov/32189668/)
- Rose JJ et al. (2020). *A neuroglobin-based high-affinity ligand trap reverses carbon monoxide-induced mitochondrial poisoning.* J Biol Chem. [PMID 32205448](https://pubmed.ncbi.nlm.nih.gov/32205448/)
- Vos MH et al. (2001). *Dynamics of nitric oxide in the active site of reduced cytochrome c oxidase aa3.* Biochemistry. [PMID 11425307](https://pubmed.ncbi.nlm.nih.gov/11425307/)

## D. Reoxygenation injury, xanthine oxidase and lipid peroxidation

> Reoxygenation injury. The model requires BOTH a prior energy failure, to convert xanthine dehydrogenase to the oxidase, AND restored oxygen as co-substrate, so the oxidative burst is generated by the treatment rather than by the poison.

- Thom SR (1993). *Leukocytes in carbon monoxide-mediated brain oxidative injury.* Toxicol Appl Pharmacol. [PMID 8248931](https://pubmed.ncbi.nlm.nih.gov/8248931/)
- Thom SR (1992). *Dehydrogenase conversion to oxidase and lipid peroxidation in brain after carbon monoxide poisoning.* J Appl Physiol (1985). [PMID 1447108](https://pubmed.ncbi.nlm.nih.gov/1447108/)
- Ischiropoulos H et al. (1996). *Nitric oxide production and perivascular nitration in brain after carbon monoxide poisoning in the rat.* J Clin Invest. [PMID 8636405](https://pubmed.ncbi.nlm.nih.gov/8636405/)
- Koehler RC et al. (2002). *Cerebrovascular effects of carbon monoxide.* Antioxid Redox Signal. [PMID 12006179](https://pubmed.ncbi.nlm.nih.gov/12006179/)
- Angelova PR et al. (2023). *Carbon monoxide neurotoxicity is triggered by oxidative stress induced by ROS production from three distinct cellular sources.* Redox Biol. [PMID 36640724](https://pubmed.ncbi.nlm.nih.gov/36640724/)
- Abramov AY et al. (2024). *Carbon Monoxide: A Pleiotropic Redox Regulator of Life and Death.* Antioxidants (Basel). [PMID 39334780](https://pubmed.ncbi.nlm.nih.gov/39334780/)
- Zhang J et al. (1992). *Mitochondrial oxidative stress after carbon monoxide hypoxia in the rat brain.* J Clin Invest. [PMID 1328293](https://pubmed.ncbi.nlm.nih.gov/1328293/)
- Wu MY et al. (2018). *Current Mechanistic Concepts in Ischemia and Reperfusion Injury.* Cell Physiol Biochem. [PMID 29694958](https://pubmed.ncbi.nlm.nih.gov/29694958/)
- Granger DN et al. (2015). *Reperfusion injury and reactive oxygen species: The evolution of a concept.* Redox Biol. [PMID 26484802](https://pubmed.ncbi.nlm.nih.gov/26484802/)
- Bagheri F et al. (2016). *Reactive oxygen species-mediated cardiac-reperfusion injury: Mechanisms and therapies.* Life Sci. [PMID 27667751](https://pubmed.ncbi.nlm.nih.gov/27667751/)

## E. Neutrophil adhesion, myeloperoxidase and the innate cascade

> The innate amplification step, and the source of the claim that hyperbaric oxygen has an action on beta-2 integrin adhesion that is neither carbon monoxide clearance nor oxygen delivery. Thom's 1993 rat work is the anchor for that term.

- ★ Thom SR (1993). *Functional inhibition of leukocyte B2 integrins by hyperbaric oxygen in carbon monoxide-mediated brain injury in rats.* Toxicol Appl Pharmacol. [PMID 8248932](https://pubmed.ncbi.nlm.nih.gov/8248932/)
- Moon JS et al. (2022). *Prognostic value of the myeloperoxidase index for early prediction of neurologic outcome in acute carbon monoxide poisoning.* Clin Exp Emerg Med. [PMID 36116774](https://pubmed.ncbi.nlm.nih.gov/36116774/)
- Thom SR et al. (1997). *Inhibition of human neutrophil beta2-integrin-dependent adherence by hyperbaric O2.* Am J Physiol. [PMID 9124510](https://pubmed.ncbi.nlm.nih.gov/9124510/)
- Jones SR et al. (2010). *Hyperbaric oxygen inhibits ischemia-reperfusion-induced neutrophil CD18 polarization by a nitric oxide mechanism.* Plast Reconstr Surg. [PMID 20679826](https://pubmed.ncbi.nlm.nih.gov/20679826/)
- Baiula M et al. (2020). *Integrin-mediated adhesive properties of neutrophils are reduced by hyperbaric oxygen therapy in patients with chronic non-healing wound.* PLoS One. [PMID 32810144](https://pubmed.ncbi.nlm.nih.gov/32810144/)

## F. Immune-mediated demyelination and delayed neurological sequelae

> The immune-mediated demyelination hypothesis, implemented here as a bistable switch with a computable threshold on adduct burden. Thom's 2004 PNAS paper, in which the delayed lesion transfers with lymphocytes, is what licenses treating adaptive immunity as the generator of delayed sequelae rather than as a bystander.

- ★ Thom SR et al. (2004). *Delayed neuropathology after carbon monoxide poisoning is immune-mediated.* Proc Natl Acad Sci U S A. [PMID 15342916](https://pubmed.ncbi.nlm.nih.gov/15342916/)
- ★ Choi IS (1983). *Delayed neurologic sequelae in carbon monoxide intoxication.* Arch Neurol. [PMID 6860181](https://pubmed.ncbi.nlm.nih.gov/6860181/)
- ★ Thom SR et al. (2006). *Hyperbaric oxygen reduces delayed immune-mediated neuropathology in experimental carbon monoxide toxicity.* Toxicol Appl Pharmacol. [PMID 16325878](https://pubmed.ncbi.nlm.nih.gov/16325878/)
- Zhang Y et al. (2024). *Association between serum neuron-specific enolase at admission and the risk of delayed neuropsychiatric sequelae in adults with carbon monoxide poisoning: A meta-analysis.* Biomol Biomed. [PMID 38850112](https://pubmed.ncbi.nlm.nih.gov/38850112/)
- Sarı Doğan F et al. (2020). *Demographic characteristics and delayed neurological sequelae risk factors in carbon monoxide poisoning.* Am J Emerg Med. [PMID 31889577](https://pubmed.ncbi.nlm.nih.gov/31889577/)
- Liao SC et al. (2018). *Predictive Role of QTc Prolongation in Carbon Monoxide Poisoning-Related Delayed Neuropsychiatric Sequelae.* Biomed Res Int. [PMID 30356348](https://pubmed.ncbi.nlm.nih.gov/30356348/)
- Han S et al. (2021). *Cox regression model of prognostic factors for delayed neuropsychiatric sequelae in patients with acute carbon monoxide poisoning: A prospective observational study.* Neurotoxicology. [PMID 33232744](https://pubmed.ncbi.nlm.nih.gov/33232744/)
- Kim JH et al. (2003). *Delayed encephalopathy of acute carbon monoxide intoxication: diffusivity of cerebral white matter lesions.* AJNR Am J Neuroradiol. [PMID 13679276](https://pubmed.ncbi.nlm.nih.gov/13679276/)
- Chang KH et al. (1992). *Delayed encephalopathy after acute carbon monoxide intoxication: MR imaging features and distribution of cerebral white matter lesions.* Radiology. [PMID 1609067](https://pubmed.ncbi.nlm.nih.gov/1609067/)

## G. Hyperbaric oxygen: randomised trials, meta-analyses and timing

> The trial literature the model is trying to reconcile. Weaver 2002 and Scheinkestel 1999 are the two arms of the controversy; the model's account is that the therapeutic window is set by adduct formation rather than by carboxyhaemoglobin, so trials differing chiefly in time-to-first-session sit on opposite sides of it.

- ★ Scheinkestel CD et al. (1999). *Hyperbaric or normobaric oxygen for acute carbon monoxide poisoning: a randomised controlled clinical trial.* Med J Aust. [PMID 10092916](https://pubmed.ncbi.nlm.nih.gov/10092916/)
- ★ Buckley NA et al. (2005). *Hyperbaric oxygen for carbon monoxide poisoning : a systematic review and critical analysis of the evidence.* Toxicol Rev. [PMID 16180928](https://pubmed.ncbi.nlm.nih.gov/16180928/)
- ★ Juurlink DN et al. (2005). *Hyperbaric oxygen for carbon monoxide poisoning.* Cochrane Database Syst Rev. [PMID 15674890](https://pubmed.ncbi.nlm.nih.gov/15674890/)
- Weaver LK (2020). *Carbon monoxide poisoning.* Undersea Hyperb Med. [PMID 32176957](https://pubmed.ncbi.nlm.nih.gov/32176957/)
- Weaver LK (2014). *Hyperbaric oxygen therapy for carbon monoxide poisoning.* Undersea Hyperb Med. [PMID 25109087](https://pubmed.ncbi.nlm.nih.gov/25109087/)
- Buckley NA et al. (2011). *Hyperbaric oxygen for carbon monoxide poisoning.* Cochrane Database Syst Rev. [PMID 21491385](https://pubmed.ncbi.nlm.nih.gov/21491385/)
- Wang W et al. (2019). *Effect of Hyperbaric Oxygen on Neurologic Sequelae and All-Cause Mortality in Patients with Carbon Monoxide Poisoning: A Meta-Analysis of Randomized Controlled Trials.* Med Sci Monit. [PMID 31606731](https://pubmed.ncbi.nlm.nih.gov/31606731/)
- Lin CH et al. (2018). *Treatment with normobaric or hyperbaric oxygen and its effect on neuropsychometric dysfunction after carbon monoxide poisoning: A systematic review and meta-analysis of randomized controlled trials.* Medicine (Baltimore). [PMID 30278526](https://pubmed.ncbi.nlm.nih.gov/30278526/)
- Ho YW et al. (2022). *Should We Use Hyperbaric Oxygen for Carbon Monoxide Poisoning Management? A Network Meta-Analysis of Randomized Controlled Trials.* Healthcare (Basel). [PMID 35885838](https://pubmed.ncbi.nlm.nih.gov/35885838/)
- Lee Y et al. (2021). *Effect of Hyperbaric Oxygen Therapy Initiation Time in Acute Carbon Monoxide Poisoning.* Crit Care Med. [PMID 34074856](https://pubmed.ncbi.nlm.nih.gov/34074856/)
- Liao SC et al. (2019). *Targeting optimal time for hyperbaric oxygen therapy following carbon monoxide poisoning for prevention of delayed neuropsychiatric sequelae: A retrospective study.* J Neurol Sci. [PMID 30481656](https://pubmed.ncbi.nlm.nih.gov/30481656/)
- Jia Y et al. (2025). *Hyperbaric oxygen therapy for acute carbon monoxide poisoning patients with coma onset.* Eur J Med Res. [PMID 39987100](https://pubmed.ncbi.nlm.nih.gov/39987100/)

## H. Guidelines, practice recommendations and epidemiology

> Guidelines and epidemiology, including the disaster and portable-generator exposures that motivate the environmental compartment.

- ★ Hampson NB et al. (2012). *Practice recommendations in the diagnosis, management, and prevention of carbon monoxide poisoning.* Am J Respir Crit Care Med. [PMID 23087025](https://pubmed.ncbi.nlm.nih.gov/23087025/)
- Jüttner B et al. (2021). *S2k guideline diagnosis and treatment of carbon monoxide poisoning.* Ger Med Sci. [PMID 34867135](https://pubmed.ncbi.nlm.nih.gov/34867135/)
- Eichhorn L et al. (2018). *The Diagnosis and Treatment of Carbon Monoxide Poisoning.* Dtsch Arztebl Int. [PMID 30765023](https://pubmed.ncbi.nlm.nih.gov/30765023/)
- Fichtner A et al. (2022). *Carbon monoxide intoxication-New aspects and current guideline-based recommendations.* Anaesthesiologie. [PMID 35925170](https://pubmed.ncbi.nlm.nih.gov/35925170/)
- de Pont AC (2006). *The guideline 'Treatment of acute carbon-monoxide poisoning' from doctors in clinics with a tank for hyperbaric ventilation.* Ned Tijdschr Geneeskd. [PMID 16613249](https://pubmed.ncbi.nlm.nih.gov/16613249/)
- Haines D (2016). *Carbon monoxide poisoning.* Med Leg J. [PMID 27130458](https://pubmed.ncbi.nlm.nih.gov/27130458/)
- Kao LW et al. (2004). *Carbon monoxide poisoning.* Emerg Med Clin North Am. [PMID 15474779](https://pubmed.ncbi.nlm.nih.gov/15474779/)
- Meredith T et al. (1988). *Carbon monoxide poisoning.* Br Med J (Clin Res Ed). [PMID 3122961](https://pubmed.ncbi.nlm.nih.gov/3122961/)
- Wu PE et al. (2014). *Carbon monoxide poisoning.* CMAJ. [PMID 24396094](https://pubmed.ncbi.nlm.nih.gov/24396094/)
- Schnall A et al. (2017). *Characterization of Carbon Monoxide Exposure During Hurricane Sandy and Subsequent Nor'easter.* Disaster Med Public Health Prep. [PMID 28438227](https://pubmed.ncbi.nlm.nih.gov/28438227/)
- Van Sickle D et al. (2007). *Carbon monoxide poisoning in Florida during the 2004 hurricane season.* Am J Prev Med. [PMID 17383566](https://pubmed.ncbi.nlm.nih.gov/17383566/)
-  (2005). *Carbon monoxide poisoning from hurricane-associated use of portable generators--Florida, 2004.* MMWR Morb Mortal Wkly Rep. [PMID 16034315](https://pubmed.ncbi.nlm.nih.gov/16034315/)

## I. Neuroimaging, biomarkers and outcome measurement

> Outcome measurement. The globus pallidus lesion motivates the watershed compartment; the biomarker literature motivates the S100B and lactate readouts.

- ★ Hopkins RO et al. (2007). *Apolipoprotein E genotype and response of carbon monoxide poisoning to hyperbaric oxygen treatment.* Am J Respir Crit Care Med. [PMID 17702967](https://pubmed.ncbi.nlm.nih.gov/17702967/)
- Wang T et al. (2022). *Neurological sequelae in acute carbon monoxide poisoning: A prospective observational study with MRI data.* Acta Neurol Scand. [PMID 35102571](https://pubmed.ncbi.nlm.nih.gov/35102571/)
- Park JH et al. (2025). *Bilateral globus pallidus lesions associated with COVID-19: Mimicking acute carbon monoxide poisoning.* Medicine (Baltimore). [PMID 40898492](https://pubmed.ncbi.nlm.nih.gov/40898492/)
- Moon JM et al. (2018). *Initial diffusion-weighted MRI and long-term neurologic outcomes in charcoal-burning carbon monoxide poisoning.* Clin Toxicol (Phila). [PMID 28753048](https://pubmed.ncbi.nlm.nih.gov/28753048/)
- Ozcan N et al. (2016). *Correlation of computed tomography, magnetic resonance imaging and clinical outcome in acute carbon monoxide poisoning.* Braz J Anesthesiol. [PMID 27591467](https://pubmed.ncbi.nlm.nih.gov/27591467/)
- Akdemir HU et al. (2014). *The role of S100B protein, neuron-specific enolase, and glial fibrillary acidic protein in the evaluation of hypoxic brain injury in acute carbon monoxide poisoning.* Hum Exp Toxicol. [PMID 24505052](https://pubmed.ncbi.nlm.nih.gov/24505052/)
- Zhang L et al. (2021). *Serum NSE and S100B protein levels for evaluating the impaired consciousness in patients with acute carbon monoxide poisoning.* Medicine (Baltimore). [PMID 34160445](https://pubmed.ncbi.nlm.nih.gov/34160445/)
- Akelma AZ et al. (2013). *Neuron-specific enolase and S100B protein in children with carbon monoxide poisoning: children are not just small adults.* Am J Emerg Med. [PMID 23380091](https://pubmed.ncbi.nlm.nih.gov/23380091/)
- Hopkins RO et al. (2006). *Neuroimaging, cognitive, and neurobehavioral outcomes following carbon monoxide poisoning.* Behav Cogn Neurosci Rev. [PMID 16891556](https://pubmed.ncbi.nlm.nih.gov/16891556/)
- Parkinson RB et al. (2002). *White matter hyperintensities and neuropsychological outcome following carbon monoxide poisoning.* Neurology. [PMID 12034791](https://pubmed.ncbi.nlm.nih.gov/12034791/)
- Lee JS et al. (2022). *Usefulness of a modified poisoning severity score for predicting prognosis in acute carbon monoxide poisoning.* Am J Emerg Med. [PMID 34739869](https://pubmed.ncbi.nlm.nih.gov/34739869/)
- Wang S et al. (2023). *Clinical application of BIS combined with LCR in assessing brain function and prognosis of patients with severe carbon monoxide poisoning.* Clin Neurol Neurosurg. [PMID 36502652](https://pubmed.ncbi.nlm.nih.gov/36502652/)
- Geng S et al. (2020). *Cardiac injury after acute carbon monoxide poisoning and its clinical treatment scheme.* Exp Ther Med. [PMID 32742349](https://pubmed.ncbi.nlm.nih.gov/32742349/)

## J. Pulse oximetry, co-oximetry and the measurement problem

> The measurement model. A two-wavelength pulse oximeter reports (O2Hb + COHb), so the displayed saturation becomes MORE reassuring as the patient becomes more poisoned.

- Papin M et al. (2023). *Accuracy of pulse CO-oximetry to evaluate blood carboxyhemoglobin level: a systematic review and meta-analysis of diagnostic test accuracy studies.* Eur J Emerg Med. [PMID 37171830](https://pubmed.ncbi.nlm.nih.gov/37171830/)
- Sebbane M et al. (2013). *Emergency department management of suspected carbon monoxide poisoning: role of pulse CO-oximetry.* Respir Care. [PMID 23513247](https://pubmed.ncbi.nlm.nih.gov/23513247/)
- Lee HY et al. (2026). *Agreement and clinical utility of non-invasive SpCO versus arterial COHb in acute carbon monoxide poisoning: a prospective observational study.* Sci Rep. [PMID 42086770](https://pubmed.ncbi.nlm.nih.gov/42086770/)

## K. Myocardial injury and long-term cardiovascular outcome

> Myocardial involvement, which the model generates from cardiac myoglobin occupancy plus terminal-oxidase inhibition in a tissue that already runs at high oxygen extraction. Henry 2006 is the anchor for the late mortality signal.

- ★ Henry CR et al. (2006). *Myocardial injury and long-term mortality following moderate to severe carbon monoxide poisoning.* JAMA. [PMID 16434630](https://pubmed.ncbi.nlm.nih.gov/16434630/)
- Patel B et al. (2023). *The Clinical Association between Carbon Monoxide Poisoning and Myocardial Injury as Measured by Elevated Troponin I Levels.* J Clin Med. [PMID 37685595](https://pubmed.ncbi.nlm.nih.gov/37685595/)
- Cho DH et al. (2024). *Practical Recommendations for the Evaluation and Management of Cardiac Injury Due to Carbon Monoxide Poisoning.* JACC Heart Fail. [PMID 38385937](https://pubmed.ncbi.nlm.nih.gov/38385937/)
- Ling YA et al. (2026). *Myocardial injury predicts mortality in elderly carbon monoxide poisoning: development and internal validation of a simple prognostic model.* Biomarkers. [PMID 42199161](https://pubmed.ncbi.nlm.nih.gov/42199161/)
- Jung YS et al. (2014). *Carbon monoxide-induced cardiomyopathy.* Circ J. [PMID 24705389](https://pubmed.ncbi.nlm.nih.gov/24705389/)
- Gülçiçek H et al. (2022). *Carbon monoxide intoxication in a shisha lounge employee, a possible cause of cardiomyopathy.* Ned Tijdschr Geneeskd. [PMID 35736345](https://pubmed.ncbi.nlm.nih.gov/35736345/)
- Lakhani M et al. (2021). *The Poisoned Heart: A Case of Takotsubo Cardiomyopathy Induced by Carbon Monoxide Poisoning.* J Emerg Med. [PMID 33674139](https://pubmed.ncbi.nlm.nih.gov/33674139/)

## L. Cyanide, fire smoke and hydroxocobalamin

> Fire smoke. Cyanide converges on the same terminal oxidase, so the two toxins compose on one axis while their antidotes share nothing. Baud 1991 is the anchor for cyanide concentrations actually measured in fire victims.

- ★ Baud FJ et al. (1991). *Elevated blood cyanide concentrations in victims of smoke inhalation.* N Engl J Med. [PMID 1944484](https://pubmed.ncbi.nlm.nih.gov/1944484/)
- ★ Borron SW et al. (2007). *Prospective study of hydroxocobalamin for acute cyanide poisoning in smoke inhalation.* Ann Emerg Med. [PMID 17481777](https://pubmed.ncbi.nlm.nih.gov/17481777/)
- Doman G et al. (2022). *Cyanide Poisoning.* J Educ Teach Emerg Med. [PMID 37465777](https://pubmed.ncbi.nlm.nih.gov/37465777/)
- Stoll S et al. (2017). *Concentrations of cyanide in blood samples of corpses after smoke inhalation of varying origin.* Int J Legal Med. [PMID 27470320](https://pubmed.ncbi.nlm.nih.gov/27470320/)
- Tabian D et al. (2021). *Toxic Blood Hydrogen Cyanide Concentration as a Vital Sign of a Deceased Room Fire Victim-Case Report.* Toxics. [PMID 33669200](https://pubmed.ncbi.nlm.nih.gov/33669200/)
- Thompson JP et al. (2012). *Hydroxocobalamin in cyanide poisoning.* Clin Toxicol (Phila). [PMID 23163594](https://pubmed.ncbi.nlm.nih.gov/23163594/)
- Mégarbane B et al. (2003). *Antidotal treatment of cyanide poisoning.* J Chin Med Assoc. [PMID 12854870](https://pubmed.ncbi.nlm.nih.gov/12854870/)

## M. Pregnancy, the fetus and paediatric exposure

> The fetal compartment, from which the bedside rule about treating for several times the maternal clearance time is derived rather than asserted.

- ★ Longo LD (1977). *The biological effects of carbon monoxide on the pregnant woman, fetus, and newborn infant.* Am J Obstet Gynecol. [PMID 561541](https://pubmed.ncbi.nlm.nih.gov/561541/)
- Friedman P et al. (2015). *Carbon Monoxide Exposure During Pregnancy.* Obstet Gynecol Surv. [PMID 26584719](https://pubmed.ncbi.nlm.nih.gov/26584719/)
- Raub JA et al. (2000). *Carbon monoxide poisoning--a public health perspective.* Toxicology. [PMID 10771127](https://pubmed.ncbi.nlm.nih.gov/10771127/)
- Eleftheriou G et al. (2022). *Open issues in management of carbon monoxide poisoning in pregnancy: practical suggestions.* J Obstet Gynaecol. [PMID 35648870](https://pubmed.ncbi.nlm.nih.gov/35648870/)
- Fleta Zaragozano J et al. (2005). *Carbon monoxide poisoning.* An Pediatr (Barc). [PMID 15927126](https://pubmed.ncbi.nlm.nih.gov/15927126/)
- Cho CH et al. (2008). *Carbon monoxide poisoning in children.* Pediatr Neonatol. [PMID 19054917](https://pubmed.ncbi.nlm.nih.gov/19054917/)
- Kaplan O et al. (2026). *Pediatric Carbon Monoxide Poisoning in Southern Israel-Causality and Outcome.* Pediatr Emerg Care. [PMID 41636140](https://pubmed.ncbi.nlm.nih.gov/41636140/)

## N. Cerebral blood flow, autoregulation and CO2 reactivity

> Cerebral haemodynamics, which set both the autoregulatory reserve and the cost of hyperventilating a poisoned brain.

- Claassen JAHR et al. (2021). *Regulation of cerebral blood flow in humans: physiology and clinical implications of autoregulation.* Physiol Rev. [PMID 33769101](https://pubmed.ncbi.nlm.nih.gov/33769101/)
- Aaslid R (2006). *Cerebral autoregulation and vasomotor reactivity.* Front Neurol Neurosci. [PMID 17290140](https://pubmed.ncbi.nlm.nih.gov/17290140/)
- Sinha AK et al. (1991). *Cerebral regional capillary perfusion and blood flow after carbon monoxide exposure.* J Appl Physiol (1985). [PMID 1757341](https://pubmed.ncbi.nlm.nih.gov/1757341/)
- Mendelman A et al. (2002). *Blood flow and ionic responses in the awake brain due to carbon monoxide.* Neurol Res. [PMID 12500698](https://pubmed.ncbi.nlm.nih.gov/12500698/)
- Lu YY et al. (2012). *Regional cerebral blood flow in patients with carbon monoxide intoxication.* Ann Nucl Med. [PMID 22872586](https://pubmed.ncbi.nlm.nih.gov/22872586/)

## O. Rhabdomyolysis, kidney injury and adjunct therapies

> Secondary organ injury and the adjuncts. N-acetylcysteine and allopurinol are in the model because they are the mechanistically obvious targets, not because human outcome evidence supports them; their predicted effects are hypotheses, not claims.

- Ito H et al. (2022). *Rhabdomyolysis secondary to carbon monoxide poisoning: A retrospective cohort study.* Am J Emerg Med. [PMID 35773173](https://pubmed.ncbi.nlm.nih.gov/35773173/)
- Al Khaldi T et al. (2023). *Acute carbon monoxide poisoning as a cause of rhabdomyolysis in a case of flame burn.* BMJ Case Rep. [PMID 37202107](https://pubmed.ncbi.nlm.nih.gov/37202107/)
- Kim SG et al. (2019). *A case report on the acute and late complications associated with carbon monoxide poisoning: Acute kidney injury, rhabdomyolysis, and delayed leukoencephalopathy.* Medicine (Baltimore). [PMID 31083215](https://pubmed.ncbi.nlm.nih.gov/31083215/)
- Howard RJ et al. (1987). *Allopurinol/N-acetylcysteine for carbon monoxide poisoning.* Lancet. [PMID 2887913](https://pubmed.ncbi.nlm.nih.gov/2887913/)
- Kim SJ et al. (2020). *Effects of Adjunctive Therapeutic Hypothermia Combined With Hyperbaric Oxygen Therapy in Acute Severe Carbon Monoxide Poisoning.* Crit Care Med. [PMID 32697512](https://pubmed.ncbi.nlm.nih.gov/32697512/)
- Kamijo Y et al. (2011). *Severe carbon monoxide poisoning complicated by hypothermia: a case report.* Am J Emerg Med. [PMID 20674229](https://pubmed.ncbi.nlm.nih.gov/20674229/)
- Oh BJ et al. (2016). *Treatment of acute carbon monoxide poisoning with induced hypothermia.* Clin Exp Emerg Med. [PMID 27752625](https://pubmed.ncbi.nlm.nih.gov/27752625/)

## P. Carbon monoxide as a signalling molecule and endogenous mediator

> Endogenous carbon monoxide and haem oxygenase-1, which the model includes as a positive feedback: the adaptive response to carbon monoxide makes more carbon monoxide.

- Cheng Y et al. (2017). *Therapeutic Potential of Heme Oxygenase-1/carbon Monoxide System Against Ischemia-Reperfusion Injury.* Curr Pharm Des. [PMID 28412905](https://pubmed.ncbi.nlm.nih.gov/28412905/)
- Otterbein LE et al. (2016). *Heme Oxygenase-1 and Carbon Monoxide in the Heart: The Balancing Act Between Danger Signaling and Pro-Survival.* Circ Res. [PMID 27283533](https://pubmed.ncbi.nlm.nih.gov/27283533/)
- Ryter SW (2020). *Therapeutic Potential of Heme Oxygenase-1 and Carbon Monoxide in Acute Organ Injury, Critical Illness, and Inflammatory Disorders.* Antioxidants (Basel). [PMID 33228260](https://pubmed.ncbi.nlm.nih.gov/33228260/)
- Dennery PA (2014). *Signaling function of heme oxygenase proteins.* Antioxid Redox Signal. [PMID 24180238](https://pubmed.ncbi.nlm.nih.gov/24180238/)
- Choi HI et al. (2022). *Controlled therapeutic delivery of CO from carbon monoxide-releasing molecules (CORMs).* J Control Release. [PMID 36063960](https://pubmed.ncbi.nlm.nih.gov/36063960/)
- Fukuda W et al. (2014). *Anti-inflammatory effects of carbon monoxide-releasing molecule on trinitrobenzene sulfonic acid-induced colitis in mice.* Dig Dis Sci. [PMID 24442266](https://pubmed.ncbi.nlm.nih.gov/24442266/)
- Ruopp M et al. (2023). *Transdermal carbon monoxide delivery.* J Control Release. [PMID 36958403](https://pubmed.ncbi.nlm.nih.gov/36958403/)

## Q. QSP methodology and simulation tooling

> Simulation tooling and the general QSP method.

- Elmokadem A et al. (2019). *Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.* CPT Pharmacometrics Syst Pharmacol. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)

---

## Notes on evidence strength

The model leans on evidence of very unequal quality and it is worth saying which is which.

**Well supported.** The Coburn-Forster-Kane formulation and the normobaric half-lives;
the oxygen-content arithmetic and the allosteric left shift; the failure of
two-wavelength pulse oximetry; the existence of delayed neurological sequelae and loss
of consciousness as its strongest clinical predictor; cyanide co-exposure in
closed-space fire.

**Supported in animals, inferred in humans.** The immune-mediated demyelination
mechanism — lipid peroxidation of myelin basic protein, charge modification, lymphocyte
sensitisation, adoptive transfer of the delayed lesion — rests principally on rodent
work from one group. The model treats it as *the* generator of delayed sequelae, which
is a strong commitment to an animal mechanism and is the single assumption most likely
to invalidate the DNS results if it is wrong.

**Contested.** Whether hyperbaric oxygen improves neurocognitive outcome. The model does
not settle this and cannot. It offers a structural reason why trials differing mainly in
time-to-first-session would reach opposite conclusions — a hypothesis about the trials,
not evidence about the therapy.

**Speculative.** The quantitative form of the tolerance threshold; the epitope-spreading
feedback that makes the switch latch; the clone kinetics that set the latency; and
N-acetylcysteine, allopurinol and carbogen as interventions. These are the parts most
likely to be wrong, and the sensitivity analysis in `co_reference_output.txt` shows the
45-day cognitive score is more sensitive to this arm than to the CO pharmacokinetics.
