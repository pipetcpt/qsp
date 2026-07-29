# Hypophosphatasia (HPP) — References

Literature underpinning the mechanistic map (`hpp_qsp_model.dot`), the ODE model
(`hpp_mrgsolve_model.R`), and the reference implementation
(`hpp_reference_model.py`).

**Verification status.** Every entry in sections 1–10 carries a **PMID that was
resolved against PubMed (E-utilities) while writing this file** — the PMID, the
journal, the year and the title were read back from PubMed rather than recalled.
Section 11 contains items for which a specific record was *not* confirmed; those
are given as PubMed **search links** only, and are explicitly labelled. Where a
model parameter rests on an unpublished or inferred value, that is stated in the
"used for" line rather than attached to a citation that does not support it.

---

## 1. Foundational pathophysiology: alkaline phosphatase and mineralization

| # | Reference | PMID | Used for |
|---|-----------|------|----------|
| 1 | Fleisch H, Bisaz S. **Mechanism of calcification: inhibitory role of pyrophosphate.** *Nature* 1962. | [13893487](https://pubmed.ncbi.nlm.nih.gov/13893487/) | The founding observation of the whole model: PPi inhibits calcification. Sets the existence and order of magnitude of `KIPPI`. |
| 2 | Meyer JL. **Can biological calcification occur in the presence of pyrophosphate?** *Arch Biochem Biophys* 1984. | [6326671](https://pubmed.ncbi.nlm.nih.gov/6326671/) | Quantitative inhibition of hydroxyapatite growth by micromolar PPi — the competitive-inhibition form `1/(1 + PPi/Ki)`. |
| 3 | Whyte MP. **Hypophosphatasia and the role of alkaline phosphatase in skeletal mineralization.** *Endocr Rev* 1994. | [7988481](https://pubmed.ncbi.nlm.nih.gov/7988481/) | The canonical synthesis: TNSALP substrates (PPi, PLP, PEA), the phenotypic spectrum, and the substrate-accumulation logic the model formalises. |
| 4 | Hessle L, et al. **Tissue-nonspecific alkaline phosphatase and plasma cell membrane glycoprotein-1 are central antagonistic regulators of bone mineralization.** *Proc Natl Acad Sci U S A* 2002. | [12082181](https://pubmed.ncbi.nlm.nih.gov/12082181/) | The TNSALP-vs-ENPP1 antagonism that makes perivesicular PPi a *local balance*, i.e. the structural basis of the `JPPI` / `KCATP*ELOC` competition. |
| 5 | Harmey D, et al. **Concerted regulation of inorganic pyrophosphate and osteopontin by akp2, enpp1, and ank: an integrated model of the pathogenesis of osteomalacia/chondromalacia.** *Am J Pathol* 2004. | [15039209](https://pubmed.ncbi.nlm.nih.gov/15039209/) | PPi–osteopontin co-regulation → the `OPN` state and its second inhibitory term; also the ANKH efflux arm of PPi production. |
| 6 | Millán JL, Whyte MP. **Alkaline Phosphatase and Hypophosphatasia.** *Calcif Tissue Int* 2016. | [26590809](https://pubmed.ncbi.nlm.nih.gov/26590809/) | Modern review tying enzyme biology to the clinical spectrum and to bone-targeted enzyme replacement. |
| 7 | Millán JL. **The role of phosphatases in the initiation of skeletal mineralization.** *Calcif Tissue Int* 2013. | [23183786](https://pubmed.ncbi.nlm.nih.gov/23183786/) | Intravesicular PHOSPHO1 vs extravesicular TNSALP division of labour — why the model puts the inhibitor term at the *perivesicular* surface. |
| 8 | Yadav MC, et al. **Loss of skeletal mineralization by the simultaneous ablation of PHOSPHO1 and alkaline phosphatase function: a unified model of the mechanisms of initiation of skeletal calcification.** *J Bone Miner Res* 2011. | [20684022](https://pubmed.ncbi.nlm.nih.gov/20684022/) | Two-step mineralization (Pi generation then propagation), justifying the separate supersaturation (`SSATN`) and inhibition (`MINH`) factors. |
| 9 | Bottini M, et al. **Matrix vesicles from chondrocytes and osteoblasts: their biogenesis, properties, functions and biomimetic models.** *Biochim Biophys Acta Gen Subj* 2018. | [29108957](https://pubmed.ncbi.nlm.nih.gov/29108957/) | The matrix vesicle as the reaction volume — the physical meaning of the `FVOL` perivesicular/plasma volume ratio. |
| 10 | **Hypophosphatasia — pathophysiological understanding, preclinical data looking beyond the skeleton, and upcoming treatments.** *J Bone Miner Res* 2026. | [41055578](https://pubmed.ncbi.nlm.nih.gov/41055578/) | Current framing of extra-skeletal HPP biology and the therapeutic pipeline. |

## 2. Genetics of ALPL

| # | Reference | PMID | Used for |
|---|-----------|------|----------|
| 11 | Del Angel G, et al. **Large-scale in vitro functional testing and novel variant scoring via protein modeling provide insights into alkaline phosphatase activity in hypophosphatasia.** *Hum Mutat* 2020. | [32160374](https://pubmed.ncbi.nlm.nih.gov/32160374/) | Direct empirical support for treating **residual enzyme activity as a continuous variable** (`FRACENZ`) rather than a categorical genotype. |
| 12 | Mornet E, et al. **A molecular-based estimation of the prevalence of hypophosphatasia in the European population.** *Ann Hum Genet* 2011. | [21488855](https://pubmed.ncbi.nlm.nih.gov/21488855/) | Prevalence of severe (~1:100,000–300,000) vs mild/heterozygous disease used in the map's genetics cluster. |
| 13 | Mornet E. **Genetics of hypophosphatasia.** *Arch Pediatr* 2017. | [29405932](https://pubmed.ncbi.nlm.nih.gov/29405932/) | Allelic architecture, dominant-negative variants, genotype–phenotype correlation. |
| 14 | Mornet E, et al. **Identification of fifteen novel mutations in the tissue-nonspecific alkaline phosphatase (TNSALP) gene in European patients with severe hypophosphatasia.** *Eur J Hum Genet* 1998. | [9781036](https://pubmed.ncbi.nlm.nih.gov/9781036/) | Variant spectrum in severe disease. |
| 15 | **ALPL Mutations With Dominant-Negative Effect in Infantile Hypophosphatasia Monozygotic Twins.** *Hum Mutat* 2026. | [41993131](https://pubmed.ncbi.nlm.nih.gov/41993131/) | Dominant-negative dimer poisoning — the `Var_domneg` → `Dimer_stoich` edge in the map. |
| 16 | Mornet E, Nunes ME. **Hypophosphatasia.** *GeneReviews* (updated). | [20301329](https://pubmed.ncbi.nlm.nih.gov/20301329/) | Clinical form definitions (perinatal, infantile, childhood, adult, odonto) used to name the scenarios. |
| 17 | **Intrafamilial phenotypic variability in hypophosphatasia: evidence from two families and the literature.** *Front Endocrinol* 2026. | [42181201](https://pubmed.ncbi.nlm.nih.gov/42181201/) | Same genotype, different severity — the model's `Somatic_mosaic` / modifier edges, and a caution on how far `FRACENZ` alone can go. |

## 3. Animal models — where the causal chain was established

| # | Reference | PMID | Used for |
|---|-----------|------|----------|
| 18 | Waymire KG, et al. **Mice lacking tissue non-specific alkaline phosphatase die from seizures due to defective metabolism of vitamin B-6.** *Nat Genet* 1995. | [7550313](https://pubmed.ncbi.nlm.nih.gov/7550313/) | The B6 axis: TNSALP is required to make pyridoxal available to the CNS. Directly motivates making CNS influx depend on **local** enzyme, not plasma enzyme. |
| 19 | Narisawa S, et al. **Inactivation of two mouse alkaline phosphatase genes and establishment of a model of infantile hypophosphatasia.** *Dev Dyn* 1997. | [9056646](https://pubmed.ncbi.nlm.nih.gov/9056646/) | The Alpl-null phenotype (rickets + seizures) that the model must reproduce jointly. |
| 20 | Millán JL, et al. **Dose response of bone-targeted enzyme replacement for murine hypophosphatasia.** *Bone* 2011. | [21458605](https://pubmed.ncbi.nlm.nih.gov/21458605/) | The only real dose–response data for D10-targeted TNSALP; used to keep HA occupancy at the label dose inside the *linear* range so dose still matters. |
| 21 | Matsumoto T, et al. **Prevention of lethal murine hypophosphatasia by neonatal ex vivo gene therapy using lentivirally transduced bone marrow cells.** *Hum Gene Ther* 2015. | [26467745](https://pubmed.ncbi.nlm.nih.gov/26467745/) | Gene-therapy arm in the map; also evidence that restoring enzyme *supply* is sufficient. |
| 22 | Yadav MC, et al. **Enzyme replacement prevents enamel defects in hypophosphatasia mice.** *J Bone Miner Res* 2012. | [22461224](https://pubmed.ncbi.nlm.nih.gov/22461224/) | Dental compartment responsiveness to ERT. |

## 4. Vitamin B6, seizures and the blood-brain barrier step

| # | Reference | PMID | Used for |
|---|-----------|------|----------|
| 23 | Baumgartner-Sigl S, et al. **Pyridoxine-responsive seizures as the first symptom of infantile hypophosphatasia caused by two novel missense mutations.** *Bone* 2007. | [17395561](https://pubmed.ncbi.nlm.nih.gov/17395561/) | Clinical anchor for the seizure module and for pyridoxine responsiveness (`KPN`, `PNDOSE`). |
| 24 | **Infantile hypophosphatasia secondary to a novel compound heterozygous mutation presenting with pyridoxine-responsive seizures.** *JIMD Rep* 2013. | [23479201](https://pubmed.ncbi.nlm.nih.gov/23479201/) | Second clinical anchor; seizures as a marker of the severe end. |
| 25 | **Vitamin B6 Status in Hypophosphatasia: Association With Clinical Severity, Diagnostic Utility, and Effects on Vitamin B6 Metabolism.** *J Inherit Metab Dis* 2025. | [40387451](https://pubmed.ncbi.nlm.nih.gov/40387451/) | Plasma PLP magnitudes across severity; directly relevant to the model's claim that plasma PLP and CNS cofactor move in *opposite* directions. |
| 26 | **Vitamin B6 challenge as a tool for detecting ALPL mutations and diagnosing hypophosphatasia.** *Osteoporos Int* 2025. | [40579471](https://pubmed.ncbi.nlm.nih.gov/40579471/) | PLP handling as a functional test of enzyme activity — the `KCATPL * EPL` term. |

## 5. Dental / periodontal compartment

| # | Reference | PMID | Used for |
|---|-----------|------|----------|
| 27 | Foster BL, et al. **Central role of pyrophosphate in acellular cementum formation.** *PLoS One* 2012. | [22675556](https://pubmed.ncbi.nlm.nih.gov/22675556/) | Why acellular cementum is the most PPi-sensitive mineralised tissue → the low `KD50` that makes premature tooth loss the earliest sign. |
| 28 | Foster BL, et al. **Overlapping functions of bone sialoprotein and pyrophosphate regulators in directing cementogenesis.** *Bone* 2017. | [28866368](https://pubmed.ncbi.nlm.nih.gov/28866368/) | Redundancy in the cementum inhibitor network. |
| 29 | **Genetic and pharmacologic modulation of cementogenesis via pyrophosphate regulators.** *Bone* 2020. | [32224162](https://pubmed.ncbi.nlm.nih.gov/32224162/) | Pharmacological reversibility of the dental phenotype. |
| 30 | **Alpl ablation in dental epithelium disrupts ameloblasts and incisor enamel mineralization in male mice.** *JBMR Plus* 2026. | [41496794](https://pubmed.ncbi.nlm.nih.gov/41496794/) | Enamel arm of the dental cluster. |

## 6. Asfotase alfa — registrational and long-term data

| # | Reference | PMID | Used for |
|---|-----------|------|----------|
| 31 | Whyte MP, et al. **Enzyme-replacement therapy in life-threatening hypophosphatasia.** *N Engl J Med* 2012. | [22397652](https://pubmed.ncbi.nlm.nih.gov/22397652/) | The registrational infantile study: radiographic and respiratory response, and the timescale (weeks–months) the model must reproduce. |
| 32 | Whyte MP, et al. **Asfotase Alfa Treatment Improves Survival for Perinatal and Infantile Hypophosphatasia.** *J Clin Endocrinol Metab* 2016. | [26529632](https://pubmed.ncbi.nlm.nih.gov/26529632/) | Survival benefit vs historical control — the calibration target for `HAZ0`, `KHZ`, `HAZS`. |
| 33 | Whyte MP, et al. **Asfotase alfa therapy for children with hypophosphatasia.** *JCI Insight* 2016. | [27699270](https://pubmed.ncbi.nlm.nih.gov/27699270/) | Childhood-HPP efficacy: growth, mobility, radiographic healing. |
| 34 | Whyte MP, et al. **Asfotase alfa for infants and young children with hypophosphatasia: 7 year outcomes of a single-arm, open-label, phase 2 extension trial.** *Lancet Diabetes Endocrinol* 2019. | [30558909](https://pubmed.ncbi.nlm.nih.gov/30558909/) | Durability of response; growth trajectories over years. |
| 35 | Kishnani PS, et al. **Five-year efficacy and safety of asfotase alfa therapy for adults and adolescents with hypophosphatasia.** *Bone* 2019. | [30576866](https://pubmed.ncbi.nlm.nih.gov/30576866/) | Adult/adolescent endpoints (6MWT, pain, fractures) used to scale `MUS` and `PAIN`. |
| 36 | Seefried L, et al. **Pharmacodynamics of asfotase alfa in adults with pediatric-onset hypophosphatasia.** *Bone* 2021. | [32987199](https://pubmed.ncbi.nlm.nih.gov/32987199/) | **Central to T2**: how the circulating PD markers behave on therapy relative to clinical response. |
| 37 | **The Effect of Asfotase Alfa on Plasma and Urine Pyrophosphate Levels and Pseudofractures in a Patient With Adult-Onset Hypophosphatasia.** *JBMR Plus* 2023. | [38130758](https://pubmed.ncbi.nlm.nih.gov/38130758/) | Direct observation of plasma/urine PPi on treatment — the measurable end of the compartment mismatch. |
| 38 | Whyte MP, et al. **Validation of a Novel Scoring System for Changes in Skeletal Manifestations of Hypophosphatasia in Newborns, Infants, and Children (RGI-C).** *J Bone Miner Res* 2018. | [29297597](https://pubmed.ncbi.nlm.nih.gov/29297597/) | The radiographic endpoint the `RSS` state is meant to stand in for, including its lag. |
| 39 | **Dual X-ray absorptiometry has limited utility in detecting bone pathology in children with hypophosphatasia: a pooled post hoc analysis.** *Bone* 2020. | [32417537](https://pubmed.ncbi.nlm.nih.gov/32417537/) | Why the model does *not* treat DXA/BMD as the primary skeletal readout — and why `BMIN` is nearly preserved while osteoid is grossly abnormal. |
| 40 | **Safety, pharmacokinetics, and pharmacodynamics of efzimfotase alfa, a second-generation enzyme replacement therapy: phase 1, dose-escalation study.** *J Bone Miner Res* 2024. | [39135540](https://pubmed.ncbi.nlm.nih.gov/39135540/) | Next-generation ERT (higher potency, longer half-life) — the alternative-PK arm in the map. |

## 7. Practical management, monitoring and real-world outcomes

| # | Reference | PMID | Used for |
|---|-----------|------|----------|
| 41 | Whyte MP. **Hypophosphatasia: Enzyme Replacement Therapy Brings New Opportunities and New Challenges.** *J Bone Miner Res* 2017. | [28084648](https://pubmed.ncbi.nlm.nih.gov/28084648/) | Ectopic calcification, injection-site reactions, craniosynostosis on therapy — the model's toxicity/limitation terms. |
| 42 | Kishnani PS, et al. **Monitoring guidance for patients with hypophosphatasia treated with asfotase alfa.** *Mol Genet Metab* 2017. | [28888853](https://pubmed.ncbi.nlm.nih.gov/28888853/) | What is monitored in practice (including vitamin B6 and ectopic calcification) — the context for the model's B6 caveat. |
| 43 | **Enzyme replacement therapy in perinatal hypophosphatasia: case report of a negative outcome and lessons for clinical practice.** *Mol Genet Metab Rep* 2018. | [29159075](https://pubmed.ncbi.nlm.nih.gov/29159075/) | Non-responders exist; a guard against over-reading the treated-arm predictions. |
| 44 | **Case report: suboptimal response to standard-dose asfotase alfa in perinatal hypophosphatasia indicates a need for individualized dosing.** *Front Endocrinol* 2025. | [40607219](https://pubmed.ncbi.nlm.nih.gov/40607219/) | Clinical counterpart of the under-dose scenario (S8) and of dose escalation to 3 mg/kg TIW (S7). |
| 45 | **Six-year clinical outcomes of enzyme replacement therapy for perinatal lethal and infantile hypophosphatasia in Korea: two case reports.** *Medicine (Baltimore)* 2023. | [36820543](https://pubmed.ncbi.nlm.nih.gov/36820543/) | Long-term real-world course including residual deficits. |
| 46 | **Effectiveness of asfotase alfa for treatment of adults with hypophosphatasia: results from a global registry.** *Orphanet J Rare Dis* 2024. | [38459585](https://pubmed.ncbi.nlm.nih.gov/38459585/) | Adult effectiveness outside trials. |
| 47 | **Clinical profiles of treated and untreated adults with hypophosphatasia in the Global HPP Registry.** *Orphanet J Rare Dis* 2022. | [35854311](https://pubmed.ncbi.nlm.nih.gov/35854311/) | Untreated adult natural history (S4). |
| 48 | **Mobility and quality of life in adults with paediatric-onset hypophosphatasia treated with asfotase alfa: results from the UK managed access agreement.** *Adv Ther* 2025. | [40138164](https://pubmed.ncbi.nlm.nih.gov/40138164/) | 6MWT-like mobility outcomes. |
| 49 | Nakamura-Utsunomiya A, et al. **Clinical characteristics of perinatal lethal hypophosphatasia: a report of 6 cases.** *Clin Pediatr Endocrinol* 2010. | [23926372](https://pubmed.ncbi.nlm.nih.gov/23926372/) | Untreated perinatal course (S1), including respiratory failure as the proximate cause of death. |

## 8. Adult disease, fractures, and the bisphosphonate question

| # | Reference | PMID | Used for |
|---|-----------|------|----------|
| 50 | Shapiro JR, Lewiecki EM. **Hypophosphatasia in Adults: Clinical Assessment and Treatment Considerations.** *J Bone Miner Res* 2017. | [28731215](https://pubmed.ncbi.nlm.nih.gov/28731215/) | Adult phenotype; the standard caution that antiresorptives are inappropriate in HPP. |
| 51 | Genest F, Seefried L. **Subtrochanteric and diaphyseal femoral fractures in hypophosphatasia — not atypical at all.** *Osteoporos Int* 2018. | [29774402](https://pubmed.ncbi.nlm.nih.gov/29774402/) | Fracture phenotype the `FX` state is scaled against. |
| 52 | **Pathophysiology of Femoral Fractures in Hypophosphatasia.** *Curr Osteoporos Rep* 2025. | [40906226](https://pubmed.ncbi.nlm.nih.gov/40906226/) | Mechanism linking osteomalacic bone quality to the fracture pattern. |
| 53 | **Atypical Fracture From Bisphosphonate Use in Hypophosphatasia With Improved Bone Response to Teriparatide Therapy.** *JCEM Case Rep* 2026. | [41503045](https://pubmed.ncbi.nlm.nih.gov/41503045/) | The clinical counterpart of **T4** (bisphosphonate harm) *and* of S15 (teriparatide benefit) in one patient. |
| 54 | **Adult hypophosphatasia treated with reduced frequency of teriparatide dosing.** *J Musculoskelet Neuronal Interact* 2021. | [34854399](https://pubmed.ncbi.nlm.nih.gov/34854399/) | Teriparatide in adult HPP — the osteoblast-expansion → endogenous-enzyme route (`ETPTD`). |
| 55 | **Lifetime follow-up of an adult patient with pediatric-onset hypophosphatasia complicated with advanced chronic kidney disease.** *Bone Rep* 2025. | [40894392](https://pubmed.ncbi.nlm.nih.gov/40894392/) | Renal endpoint of the nephrocalcinosis pathway. |

## 9. Contrast conditions that constrain the model

| # | Reference | PMID | Used for |
|---|-----------|------|----------|
| 56 | **Generalized Arterial Calcification of Infancy (GACI): state of the art and clinical perspectives.** *J Cardiovasc Dev Dis* 2026. | [42188070](https://pubmed.ncbi.nlm.nih.gov/42188070/) | The mirror-image disease: ENPP1 loss → PPi *deficiency* → pathological calcification. It is the reason the model carries an ectopic-calcification term that activates when PPi is driven *below* normal. |

## 10. Modelling methodology

| # | Reference | PMID | Used for |
|---|-----------|------|----------|
| 57 | Elmokadem A, Riggs MM, Baron KT. **Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.** *CPT Pharmacometrics Syst Pharmacol* 2019. | [31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/) | mrgsolve implementation idioms used in `hpp_mrgsolve_model.R`. |

---

## 11. Sources used but **not** resolved to a specific PubMed record

These are cited honestly as searches, not as records. Two of them carry real
weight in the model, so the uncertainty matters:

- **Strensiq (asfotase alfa) prescribing information** — the source of
  `FBIO = 0.458`, tmax 24–48 h and terminal t½ = 2.28 d, and of the label
  regimens (2 mg/kg three times weekly; 1 mg/kg six times weekly; escalation to
  3 mg/kg three times weekly). A drug label is not a PubMed record; the values
  should be checked against the current label for the relevant jurisdiction.
  [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=asfotase+alfa+pharmacokinetics+bioavailability+half-life)
- **Whyte MP. Hypophosphatasia — aetiology, nosology, pathogenesis, diagnosis
  and treatment.** *Nat Rev Endocrinol* (review). Used for the clinical-form
  taxonomy; the specific record was not confirmed here.
  [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Whyte+hypophosphatasia+aetiology+nosology+pathogenesis+diagnosis+treatment)
- **Asfotase alfa immunogenicity / anti-drug antibodies.** The `ADA` module's
  structure (most patients seroconvert, a minority show neutralising titres with
  reduced exposure) is taken from label-level safety reporting rather than a
  confirmed primary paper. `KADAON`, `KADAOFF`, `GAMADA` and `ADAMAX` are
  therefore **assumed**, and the high-titre scenario (S11) should be read as a
  sensitivity analysis, not a prediction.
  [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=asfotase+alfa+anti-drug+antibodies+immunogenicity)
- **Consensus management guidance for HPP in adults / children.** Referenced in
  the map's therapeutic cluster (multidisciplinary care, orthopaedic
  load-sharing fixation, calcium/vitamin-D restriction in infantile
  hypercalcaemia) without a confirmed record.
  [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=hypophosphatasia+management+consensus+recommendations)

## 12. Parameters that no citation supports

Stated here so they are not mistaken for literature values. All four are the
model's main identifiability gap and are discussed in README §8:

| Parameter | Meaning | Status |
|-----------|---------|--------|
| `KON`, `KOFF`, `KDEGB` | plasma ↔ hydroxyapatite binding and turnover of the bound enzyme | **Assumed.** Chosen so the bound pool has a ~3.5-day half-life and the label dose sits inside the linear range of occupancy. |
| `BMAXKG` | hydroxyapatite binding capacity | **Assumed.** Absolute occupancy values are therefore a *relative* scale only. |
| `KBONE` | normalised local activity per unit occupancy | **Assumed.** It sets how far above healthy the local activity can be pushed. |
| `KCS`, `TCS` | craniosynostosis progression | **Phenomenological.** Reproduces "common, and not reliably prevented by ERT" by construction, and explains nothing. |

---

### How to read this list against the model

- **T1 (threshold)** rests on refs 1, 2, 4, 5, 7, 8 for structure and ref 11 for
  the continuity of residual activity.
- **T2 (wrong compartment)** rests on refs 20, 36, 37, 39 — the dose–response of
  a *bone-targeted* enzyme, the behaviour of circulating markers on therapy, and
  the demonstrated weakness of BMD as a readout.
- **T3 (two clocks)** rests on refs 31, 32, 34, 41, 44, 45: radiographic and
  respiratory recovery within months, against permanent stature deficit and
  craniosynostosis that ERT does not reliably prevent.
- **T4 (bisphosphonate harm)** rests on refs 1, 2 (PPi analogue chemistry) and
  50, 51, 53 (clinical fracture phenotype).
- **The B6 dissociation** rests on refs 18, 23, 25, 42 — and its extension to
  "ERT may lower the CNS substrate gradient" is the model's own inference, with
  no supporting clinical data.
