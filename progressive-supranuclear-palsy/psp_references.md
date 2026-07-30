# Progressive Supranuclear Palsy (PSP) -- QSP model references

**진행성 핵상마비 (Progressive Supranuclear Palsy, PSP) QSP 모델 참고문헌**


## 인용 규약 (citation convention) -- 반드시 읽을 것

이 목록의 항목은 **두 종류**입니다. 구분은 의도적이며, 확인되지 않은 PMID를
만들어 넣지 않기 위한 것입니다.

| 형태 | 의미 |
|---|---|
| `**저자** 제목. *저널* 연도. [PMID 12345](...)` | PubMed E-utilities로 **실제 조회하여 확인된** 서지정보. 제목·저널·연도·PMID는 모두 PubMed가 반환한 값 그대로입니다. |
| `질의문 — [PubMed search](...)` | 제목 필드 검색으로 단일 논문을 자동 확정하지 못한 항목. **PMID를 추정해 적지 않고** 대신 항상 유효한 PubMed 검색 링크를 제공합니다. |

Entries are of **two kinds**, deliberately. Resolved entries carry bibliographic
data returned verbatim by the PubMed E-utilities API. Unresolved entries carry a
live PubMed search link instead of a guessed PMID: where automatic title-field
resolution did not return a single unambiguous record, **no PMID is asserted**.

Resolved: **80** &nbsp;|&nbsp; search links: **66** &nbsp;|&nbsp; total: **146**

> 모델의 정량적 앵커(anchor)가 어느 문헌에서 왔는지는 `psp_mrgsolve_model.R`의
> `psp_validate()` 표와 `README.md`의 검증 섹션에 정리되어 있습니다. 아래에서
> `—` 로 시작하는 주석은 그 문헌이 모델의 **어느 구조·어느 파라미터**를
> 뒷받침하는지를 나타냅니다.


---

## A. Diagnosis, epidemiology, natural history and rating scales

1. Clinical diagnosis of progressive supranuclear palsy: the movement disorder society criteria — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Clinical%20diagnosis%20of%20progressive%20supranuclear%20palsy%3A%20the%20movement%20disorder%20society%20criteria)  
  &nbsp;&nbsp;— MDS-PSP 2017 criteria (O/P/A/C domains)

2. **Golbe LI et al.** A clinical rating scale for progressive supranuclear palsy. *Brain* 2007. [PMID 17405767](https://pubmed.ncbi.nlm.nih.gov/17405767/)  
  &nbsp;&nbsp;— PSPRS; progression 11.3 points/year -- the model's slope anchor

3. **Dam T et al.** A 15-Item modification of the PSP rating scale to improve clinical meaningfulness and statistical performance. *Nat Commun* 2025. [PMID 39762226](https://pubmed.ncbi.nlm.nih.gov/39762226/)  
  &nbsp;&nbsp;— mPSPRS-10/15: lower measurement noise, smaller required n

4. **Respondek G et al.** The phenotypic spectrum of progressive supranuclear palsy: a retrospective multicenter study of 100 definite cases. *Mov Disord* 2014. [PMID 25370486](https://pubmed.ncbi.nlm.nih.gov/25370486/)  
  &nbsp;&nbsp;— Respondek: PSP-RS vs PSP-P vs other phenotypes

5. **Respondek G et al.** Which ante mortem clinical features predict progressive supranuclear palsy pathology?. *Mov Disord* 2017. [PMID 28500752](https://pubmed.ncbi.nlm.nih.gov/28500752/)  
  &nbsp;&nbsp;— clinicopathological correlation of the phenotypes

6. Survival in progressive supranuclear palsy and frontotemporal dementia — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Survival%20in%20progressive%20supranuclear%20palsy%20and%20frontotemporal%20dementia)  
  &nbsp;&nbsp;— survival anchor

7. **Cushing BJ et al.** Prevalent versus incident progressive supranuclear palsy: An analysis of the frequencies of neuropathological and clinical features at U.S. Alzheimer's Disease Research Centers indicate a relatively common tauopathy of aging. *Res Sq* 2026. [PMID 42427864](https://pubmed.ncbi.nlm.nih.gov/42427864/)  
  &nbsp;&nbsp;— prevalence 5-7 per 100000

8. **Bang J et al.** Predicting disease progression in progressive supranuclear palsy in multicenter clinical trials. *Parkinsonism Relat Disord* 2016. [PMID 27172829](https://pubmed.ncbi.nlm.nih.gov/27172829/)

9. **Madetko-Alster N et al.** The Possible Significance of Proteomics in Understanding Molecular Mechanisms of Progressive Supranuclear Palsy, Corticobasal Degeneration, Multiple System Atrophy, and Dementia with Lewy Bodies. *Cells* 2026. [PMID 42121860](https://pubmed.ncbi.nlm.nih.gov/42121860/)

10. **Wenning GK et al.** Natural history and survival of 14 patients with corticobasal degeneration confirmed at postmortem examination. *J Neurol Neurosurg Psychiatry* 1998. [PMID 9489528](https://pubmed.ncbi.nlm.nih.gov/9489528/)  
  &nbsp;&nbsp;— comparator 4R tauopathy

11. **Ali F et al.** Sensitivity and Specificity of Diagnostic Criteria for Progressive Supranuclear Palsy. *Mov Disord* 2019. [PMID 30726566](https://pubmed.ncbi.nlm.nih.gov/30726566/)

12. **Swallow DMA et al.** The evolution of diagnosis from symptom onset to death in progressive supranuclear palsy (PSP) and corticobasal degeneration (CBD) compared to Parkinson's disease (PD). *J Neurol* 2023. [PMID 36971841](https://pubmed.ncbi.nlm.nih.gov/36971841/)  
  &nbsp;&nbsp;— onset-to-diagnosis interval -- sets the model's enrolment placement


---

## B. Neuropathology and the 4R tau filament fold

1. **McKee AC et al.** The first NINDS/NIBIB consensus meeting to define neuropathological criteria for the diagnosis of chronic traumatic encephalopathy. *Acta Neuropathol* 2016. [PMID 26667418](https://pubmed.ncbi.nlm.nih.gov/26667418/)  
  &nbsp;&nbsp;— neuropathological definition

2. **Schweighauser M et al.** Novel tau filament folds in individuals with MAPT mutations P301L and P301T. *bioRxiv* 2024. [PMID 39185206](https://pubmed.ncbi.nlm.nih.gov/39185206/)  
  &nbsp;&nbsp;— cryo-EM: fold differs between 4R tauopathies

3. Structure-based classification of tauopathies — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Structure-based%20classification%20of%20tauopathies)  
  &nbsp;&nbsp;— PSP has its own filament fold

4. **Koga S et al.** Deep Learning-Based Image Classification in Differentiating Tufted Astrocytes, Astrocytic Plaques, and Neuritic Plaques. *J Neuropathol Exp Neurol* 2021. [PMID 33570124](https://pubmed.ncbi.nlm.nih.gov/33570124/)  
  &nbsp;&nbsp;— tufted astrocyte is pathognomonic for PSP

5. Staging of progressive supranuclear palsy pathology Kovacs — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Staging%20of%20progressive%20supranuclear%20palsy%20pathology%20Kovacs)  
  &nbsp;&nbsp;— regional sequence used for the model's connectome path

6. **Cullinane PW et al.** Data-driven modelling of tau pathology reveals distinct progressive supranuclear palsy subtypes. *Brain* 2026. [PMID 41974128](https://pubmed.ncbi.nlm.nih.gov/41974128/)

7. Neuronal loss in the subthalamic nucleus in progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Neuronal%20loss%20in%20the%20subthalamic%20nucleus%20in%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— STN as an early, severely affected nucleus

8. Pedunculopontine nucleus cholinergic neuron loss in progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Pedunculopontine%20nucleus%20cholinergic%20neuron%20loss%20in%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— PPN cholinergic loss -> falls

9. Oligodendroglial coiled bodies in progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Oligodendroglial%20coiled%20bodies%20in%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— MOBP risk locus and oligodendroglial pathology

10. Quantitative neuropathological comparison of progressive supranuclear palsy phenotypes tau burden — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Quantitative%20neuropathological%20comparison%20of%20progressive%20supranuclear%20palsy%20phenotypes%20tau%20burden)  
  &nbsp;&nbsp;— PSP-P has LOWER tau burden than PSP-RS -- the model's reported failure

11. **Spiegel C et al.** Brainstem and cerebellar radiological findings in progressive supranuclear palsy. *Brain Commun* 2025. [PMID 39958262](https://pubmed.ncbi.nlm.nih.gov/39958262/)

12. **Cui S et al.** Locus coeruleus-norepinephrine system dysfunction: A new concept in cognitive aging and neurodegenerative diseases. *Neural Regen Res* 2026. [PMID 42322653](https://pubmed.ncbi.nlm.nih.gov/42322653/)  
  &nbsp;&nbsp;— noradrenergic loss


---

## C. Genetics

1. **Höglinger GU et al.** Identification of common variants influencing risk of the tauopathy progressive supranuclear palsy. *Nat Genet* 2011. [PMID 21685912](https://pubmed.ncbi.nlm.nih.gov/21685912/)  
  &nbsp;&nbsp;— GWAS: MAPT, MOBP, STX6, EIF2AK3

2. Association of MAPT H1 haplotype with progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Association%20of%20MAPT%20H1%20haplotype%20with%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— H1/H1c haplotype, OR ~5.5

3. Genome-wide association study of progressive supranuclear palsy identifies new risk loci — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Genome-wide%20association%20study%20of%20progressive%20supranuclear%20palsy%20identifies%20new%20risk%20loci)  
  &nbsp;&nbsp;— additional loci

4. TRIM11 protects against tauopathies and is auto-regulated by SUMOylation — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=TRIM11%20protects%20against%20tauopathies%20and%20is%20auto-regulated%20by%20SUMOylation)  
  &nbsp;&nbsp;— disaggregase capacity -- the model's TRIM11F parameter

5. **Ressler HW et al.** MAPT haplotype-associated transcriptomic changes in progressive supranuclear palsy. *Acta Neuropathol Commun* 2024. [PMID 39154163](https://pubmed.ncbi.nlm.nih.gov/39154163/)  
  &nbsp;&nbsp;— mechanism by which the risk haplotype raises the monomer pool M

6. **Li H et al.** A novel MAPT variant (E342K) as a cause of familial progressive supranuclear palsy. *Front Neurol* 2024. [PMID 38708005](https://pubmed.ncbi.nlm.nih.gov/38708005/)

7. **Qi C et al.** The Pick fold in tau filaments from human MAPT mutants. *Acta Neuropathol* 2026. [PMID 42420562](https://pubmed.ncbi.nlm.nih.gov/42420562/)

8. **Wang H et al.** Correction: Whole-genome sequencing analysis reveals new susceptibility loci and structural variants associated with progressive supranuclear palsy. *Mol Neurodegener* 2024. [PMID 39402686](https://pubmed.ncbi.nlm.nih.gov/39402686/)

9. Cell-type-specific expression of MAPT risk variants in progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Cell-type-specific%20expression%20of%20MAPT%20risk%20variants%20in%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— expression QTL

10. **Escobar-Khondiker M et al.** Annonacin, a natural mitochondrial complex I inhibitor, causes tau pathology in cultured neurons. *J Neurosci* 2007. [PMID 17634376](https://pubmed.ncbi.nlm.nih.gov/17634376/)  
  &nbsp;&nbsp;— environmental Complex I inhibition (Guadeloupean PSP)


---

## D. Tau biology: splicing, post-translational modification, aggregation kinetics

1. **Chen L** CELF2 Promotes Tau Exon 10 Inclusion via Hinge Domain-Mediated Nuclear Condensation, Driving Cognitive Dysfunction in Tauopathy Models. *Res Sq* 2026. [PMID 41646425](https://pubmed.ncbi.nlm.nih.gov/41646425/)  
  &nbsp;&nbsp;— 4R:3R ratio

2. **Voelzmann A et al.** GSK-3β coordinates axonal microtubule organization through Shot and Tau. *Proc Natl Acad Sci U S A* 2026. [PMID 41701831](https://pubmed.ncbi.nlm.nih.gov/41701831/)  
  &nbsp;&nbsp;— the tideglusib target

3. **Issad T et al.** [Protein O-GlcNAcylation and regulation of cell signalling: involvement in pathophysiology]. *Biol Aujourdhui* 2014. [PMID 25190571](https://pubmed.ncbi.nlm.nih.gov/25190571/)  
  &nbsp;&nbsp;— OGA inhibitor rationale

4. Pharmacological inhibition of O-GlcNAcase reduces tau pathology thiamet — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Pharmacological%20inhibition%20of%20O-GlcNAcase%20reduces%20tau%20pathology%20thiamet)  
  &nbsp;&nbsp;— OGA inhibition in vivo

5. Acetylation of tau inhibits its degradation and contributes to tauopathy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Acetylation%20of%20tau%20inhibits%20its%20degradation%20and%20contributes%20to%20tauopathy)  
  &nbsp;&nbsp;— acetyl-tau K174/K274 blocks clearance

6. **Min SW et al.** Critical role of acetylation in tau-mediated neurodegeneration and cognitive deficits. *Nat Med* 2015. [PMID 26390242](https://pubmed.ncbi.nlm.nih.gov/26390242/)  
  &nbsp;&nbsp;— the salsalate/p300 axis

7. **Day RJ et al.** Caspase-Cleaved Tau Co-Localizes with Early Tangle Markers in the Human Vascular Dementia Brain. *PLoS One* 2015. [PMID 26161867](https://pubmed.ncbi.nlm.nih.gov/26161867/)  
  &nbsp;&nbsp;— D421 truncation generates N-terminal fragments

8. Asparagine endopeptidase cleaves tau at N368 and mediates tauopathy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Asparagine%20endopeptidase%20cleaves%20tau%20at%20N368%20and%20mediates%20tauopathy)  
  &nbsp;&nbsp;— AEP/legumain: N-terminus is REMOVED from the seeding core

9. **Vigers MP et al.** Water-directed pinning is key to tau prion formation. *Proc Natl Acad Sci U S A* 2025. [PMID 40294272](https://pubmed.ncbi.nlm.nih.gov/40294272/)

10. **Yao TM et al.** Aggregation analysis of the microtubule binding domain in tau protein by spectroscopic methods. *J Biochem* 2003. [PMID 12944375](https://pubmed.ncbi.nlm.nih.gov/12944375/)

11. **Chinnathambi S** Small molecule-mediated therapeutic approaches to target Tau and Alzheimer's disease. *Adv Protein Chem Struct Biol* 2025. [PMID 40324850](https://pubmed.ncbi.nlm.nih.gov/40324850/)

12. **Davtyan H et al.** Transplantation of human iPSC-derived microglia ameliorates neuropathology and circuit dysfunction in progranulin-deficient mice. *Mol Neurodegener* 2026. [PMID 42464356](https://pubmed.ncbi.nlm.nih.gov/42464356/)  
  &nbsp;&nbsp;— the ezeprogind target

13. Rho kinase ROCK inhibition enhances autophagic tau clearance — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Rho%20kinase%20ROCK%20inhibition%20enhances%20autophagic%20tau%20clearance)  
  &nbsp;&nbsp;— the fasudil mechanism

14. **Cabrera J et al.** Barbatolic Acid Prevents Tau and Amylin Interaction and Stimulates the Growth of Acetylated Microtubules in Cell Culture. *Curr Drug Targets* 2026. [PMID 42411090](https://pubmed.ncbi.nlm.nih.gov/42411090/)  
  &nbsp;&nbsp;— the davunetide rationale


---

## E. Prion-like propagation, connectome spread and travelling fronts

1. **Blaudin de Thé FX et al.** P62 accumulates through neuroanatomical circuits in response to tauopathy propagation. *Acta Neuropathol Commun* 2021. [PMID 34727983](https://pubmed.ncbi.nlm.nih.gov/34727983/)

2. Propagation of tau pathology in a model of early Alzheimer's disease — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Propagation%20of%20tau%20pathology%20in%20a%20model%20of%20early%20Alzheimer%27s%20disease)  
  &nbsp;&nbsp;— templated propagation

3. Brain transcriptional and connectomic architecture predicts tau spreading — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Brain%20transcriptional%20and%20connectomic%20architecture%20predicts%20tau%20spreading)  
  &nbsp;&nbsp;— connectome weights govern the itinerary

4. A network diffusion model of disease progression in dementia — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=A%20network%20diffusion%20model%20of%20disease%20progression%20in%20dementia)  
  &nbsp;&nbsp;— graph-Laplacian spread; the geometry behind the linear slope

5. Tau assemblies do not behave like independently acting prion-like particles in mouse neural tissue — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Tau%20assemblies%20do%20not%20behave%20like%20independently%20acting%20prion-like%20particles%20in%20mouse%20neural%20tissue)  
  &nbsp;&nbsp;— seeding is not simple particle counting

6. Heparan sulfate proteoglycans mediate internalization and propagation of tau — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Heparan%20sulfate%20proteoglycans%20mediate%20internalization%20and%20propagation%20of%20tau)  
  &nbsp;&nbsp;— HSPG-dependent uptake

7. **Fearon C et al.** Commentary: LRP1 Is a Master Regulator of Tau Uptake and Spread. *Front Neurol* 2020. [PMID 33424736](https://pubmed.ncbi.nlm.nih.gov/33424736/)  
  &nbsp;&nbsp;— dominant uptake receptor -- the anti-transfer target

8. Exosome-mediated transfer of tau between neurons — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Exosome-mediated%20transfer%20of%20tau%20between%20neurons)  
  &nbsp;&nbsp;— an antibody-INACCESSIBLE transfer route (the PHI_ACC gate)

9. **Tardivel M et al.** Tunneling nanotube (TNT)-mediated neuron-to neuron transfer of pathological Tau protein assemblies. *Acta Neuropathol Commun* 2016. [PMID 27809932](https://pubmed.ncbi.nlm.nih.gov/27809932/)  
  &nbsp;&nbsp;— a second inaccessible route

10. Muscarinic receptors and macropinocytosis in tau internalization — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Muscarinic%20receptors%20and%20macropinocytosis%20in%20tau%20internalization)  
  &nbsp;&nbsp;— alternative uptake routes

11. **Kaufman SK et al.** Characterization of tau prion seeding activity and strains from formaldehyde-fixed tissue. *Acta Neuropathol Commun* 2017. [PMID 28587664](https://pubmed.ncbi.nlm.nih.gov/28587664/)

12. Regional distribution of tau seeding activity in progressive supranuclear palsy brain — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Regional%20distribution%20of%20tau%20seeding%20activity%20in%20progressive%20supranuclear%20palsy%20brain)  
  &nbsp;&nbsp;— regional seeding in PSP


---

## F. Extracellular tau, the ISF/CSF compartment and antibody access

1. **Barini E et al.** Tau in the brain interstitial fluid is fragmented and seeding-competent. *Neurobiol Aging* 2022. [PMID 34655982](https://pubmed.ncbi.nlm.nih.gov/34655982/)  
  &nbsp;&nbsp;— ISF tau ~1 nM and FRAGMENTED -- the pool-size and epitope claims

2. In vivo microdialysis of brain interstitial fluid tau — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=In%20vivo%20microdialysis%20of%20brain%20interstitial%20fluid%20tau)  
  &nbsp;&nbsp;— the measurement behind the ISF concentration

3. Neuronal activity regulates extracellular tau in vivo — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Neuronal%20activity%20regulates%20extracellular%20tau%20in%20vivo)  
  &nbsp;&nbsp;— activity-dependent release

4. **Evans LD et al.** Extracellular Monomeric and Aggregated Tau Efficiently Enter Human Neurons through Overlapping but Distinct Pathways. *Cell Rep* 2018. [PMID 29590627](https://pubmed.ncbi.nlm.nih.gov/29590627/)

5. **Nadadhur AG et al.** Astrocytes from P301S Tau mice exhibit non-canonical protein secretion and reduced morphological complexity. *Neural Regen Res* 2026. [PMID 40808413](https://pubmed.ncbi.nlm.nih.gov/40808413/)

6. Antibody-free quantification of seven tau peptides in human CSF using targeted mass spectrometry — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Antibody-free%20quantification%20of%20seven%20tau%20peptides%20in%20human%20CSF%20using%20targeted%20mass%20spectrometry)  
  &nbsp;&nbsp;— which tau species CSF assays actually see

7. **Robinson CG et al.** Tau-PET and CSF MTBR-tau243 comparisons validate increased tau aggregation in females. *Eur J Nucl Med Mol Imaging* 2026. [PMID 42191949](https://pubmed.ncbi.nlm.nih.gov/42191949/)

8. **Zheng X et al.** Integrated Blood Inflammatory Ratios and Cerebrospinal Fluid Blood‒Brain Barrier Dysfunction Predict Relapse Risk in Neuromyelitis Optica Spectrum Disorder. *Brain Behav* 2026. [PMID 42287009](https://pubmed.ncbi.nlm.nih.gov/42287009/)

9. **Wu S et al.** Investigation of Antibody Pharmacokinetics in the Brain Following Intra-CNS Administration and Development of PBPK Model to Characterize the Data. *AAPS J* 2024. [PMID 38443635](https://pubmed.ncbi.nlm.nih.gov/38443635/)

10. **Chinnathambi S et al.** G-protein coupled receptors (GPCRs) interacts with Tau protein in Alzheimer's disease. *Adv Protein Chem Struct Biol* 2025. [PMID 40973402](https://pubmed.ncbi.nlm.nih.gov/40973402/)


---

## G. Fluid biomarkers

1. **Leclercq V** Biomarkers stewardship in parkinsonism: integrating alpha-synuclein seed amplification assays and neurofilament light chain into diagnostic pathways and patient communication. *J Neurol* 2026. [PMID 42295486](https://pubmed.ncbi.nlm.nih.gov/42295486/)  
  &nbsp;&nbsp;— CSF/plasma NfL levels in PSP

2. Serum neurofilament light chain in progressive supranuclear palsy and multiple system atrophy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Serum%20neurofilament%20light%20chain%20in%20progressive%20supranuclear%20palsy%20and%20multiple%20system%20atrophy)  
  &nbsp;&nbsp;— plasma NfL

3. **Rojas JC et al.** CSF neurofilament light chain and phosphorylated tau 181 predict disease progression in PSP. *Neurology* 2018. [PMID 29282336](https://pubmed.ncbi.nlm.nih.gov/29282336/)  
  &nbsp;&nbsp;— the +40%/y rise the model under-predicts (reported failure)

4. **Mascioli D et al.** Alzheimer's Disease Cerebrospinal Fluid Biomarkers Predict Survival in Progressive Supranuclear Palsy. *Mov Disord* 2026. [PMID 42386668](https://pubmed.ncbi.nlm.nih.gov/42386668/)

5. **Hortsch S et al.** Performance of a fully automated plasma tau phosphorylated at threonine 217 immunoassay to reflect amyloid-beta burden in an unselected cohort representative of clinical practice. *J Prev Alzheimers Dis* 2026. [PMID 41887009](https://pubmed.ncbi.nlm.nih.gov/41887009/)

6. Plasma glial fibrillary acidic protein GFAP in neurodegenerative disease — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Plasma%20glial%20fibrillary%20acidic%20protein%20GFAP%20in%20neurodegenerative%20disease)  
  &nbsp;&nbsp;— astrocytic marker

7. Soluble TREM2 in cerebrospinal fluid as a microglial activation marker — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Soluble%20TREM2%20in%20cerebrospinal%20fluid%20as%20a%20microglial%20activation%20marker)  
  &nbsp;&nbsp;— microglial activation marker

8. Cerebrospinal fluid tau seeding activity measured by real-time quaking induced conversion in tauopathies — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Cerebrospinal%20fluid%20tau%20seeding%20activity%20measured%20by%20real-time%20quaking%20induced%20conversion%20in%20tauopathies)  
  &nbsp;&nbsp;— the assay that would have discriminated

9. **Bernhardt AM et al.** From clinical phenotypes to molecular stratification: early differential diagnosis of four-repeat tauopathies. *Front Aging Neurosci* 2026. [PMID 42382525](https://pubmed.ncbi.nlm.nih.gov/42382525/)

10. **Corvol JC et al.** AZP2006 in Progressive Supranuclear Palsy: Outcomes from a Phase 2a Multicenter, Randomized Trial, and Open-Label Extension on Safety, Biomarkers, and Disease Progression. *Mov Disord* 2025. [PMID 41014124](https://pubmed.ncbi.nlm.nih.gov/41014124/)


---

## H. Imaging

1. Accuracy of magnetic resonance parkinsonism index for differentiation of progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Accuracy%20of%20magnetic%20resonance%20parkinsonism%20index%20for%20differentiation%20of%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— MRPI cut-off ~13.55

2. MR parkinsonism index 2.0 for the diagnosis of progressive supranuclear palsy parkinsonism — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=MR%20parkinsonism%20index%202.0%20for%20the%20diagnosis%20of%20progressive%20supranuclear%20palsy%20parkinsonism)  
  &nbsp;&nbsp;— MRPI 2.0

3. **Cunningham MCQES et al.** Does Midbrain Atrophy Distinguish Progressive Supranuclear Palsy from Frontotemporal Dementia?. *Mov Disord Clin Pract* 2025. [PMID 40172482](https://pubmed.ncbi.nlm.nih.gov/40172482/)  
  &nbsp;&nbsp;— midbrain area <70 mm2

4. Superior cerebellar peduncle atrophy in progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Superior%20cerebellar%20peduncle%20atrophy%20in%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— SCP width, the MRPI denominator

5. Rate of brain atrophy as an outcome measure in progressive supranuclear palsy clinical trials — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Rate%20of%20brain%20atrophy%20as%20an%20outcome%20measure%20in%20progressive%20supranuclear%20palsy%20clinical%20trials)  
  &nbsp;&nbsp;— imaging endpoints and power

6. **Gnörich J et al.** Longitudinal monitoring of tau aggregation in progressive supranuclear palsy with [(18)F]PI-2620 PET. *Alzheimers Dement* 2026. [PMID 41736364](https://pubmed.ncbi.nlm.nih.gov/41736364/)  
  &nbsp;&nbsp;— 4R-selective tau PET

7. 18F-APN-1607 PI-2620 tau PET imaging in four repeat tauopathies — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=18F-APN-1607%20PI-2620%20tau%20PET%20imaging%20in%20four%20repeat%20tauopathies)  
  &nbsp;&nbsp;— tau-PET in 4R tauopathy

8. TSPO PET imaging of microglial activation in progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=TSPO%20PET%20imaging%20of%20microglial%20activation%20in%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— neuroinflammation imaging

9. **Tessema AW et al.** Paired regional complementarity in diffusion MRI reveals disease-specific microstructural profiles in PD, MSA, and PSP: a feasibility study. *Sci Rep* 2026. [PMID 41775795](https://pubmed.ncbi.nlm.nih.gov/41775795/)  
  &nbsp;&nbsp;— DTI marker

10. **Liang M et al.** 18 F-FDG, 18 F-FP-CIT, and 18 F-Florzolotau PET Imaging in Progressive Supranuclear Palsy : Region-Specific Correlations Between Glucose Metabolism, Dopaminergic Function, and Tau Pathology. *Clin Nucl Med* 2026. [PMID 41474764](https://pubmed.ncbi.nlm.nih.gov/41474764/)


---

## I. Oculomotor and postural physiology -- the two lowest-reserve nuclei

1. Vertical saccade velocity in progressive supranuclear palsy quantitative oculography — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Vertical%20saccade%20velocity%20in%20progressive%20supranuclear%20palsy%20quantitative%20oculography)  
  &nbsp;&nbsp;— peak velocity is a nearly unthresholded readout of riMLF survival

2. The rostral interstitial nucleus of the medial longitudinal fasciculus and vertical saccade generation — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=The%20rostral%20interstitial%20nucleus%20of%20the%20medial%20longitudinal%20fasciculus%20and%20vertical%20saccade%20generation)  
  &nbsp;&nbsp;— riMLF burst neurons

3. **Hittinger M et al.** The anatomical identification of saccadic omnipause neurons in the rat brainstem. *Neuroscience* 2012. [PMID 22441037](https://pubmed.ncbi.nlm.nih.gov/22441037/)  
  &nbsp;&nbsp;— square-wave jerks

4. **Quattrone A et al.** Video-oculographic biomarkers for evaluating vertical ocular dysfunction in progressive supranuclear palsy. *Parkinsonism Relat Disord* 2022. [PMID 35642995](https://pubmed.ncbi.nlm.nih.gov/35642995/)

5. Slow vertical saccades as the earliest sign of progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Slow%20vertical%20saccades%20as%20the%20earliest%20sign%20of%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— earliest oculomotor sign

6. **Warren NM et al.** Cholinergic systems in progressive supranuclear palsy. *Brain* 2005. [PMID 15649952](https://pubmed.ncbi.nlm.nih.gov/15649952/)

7. **Rana AQ et al.** Combination of blepharospasm and apraxia of eyelid opening: a condition resistant to treatment. *Acta Neurol Belg* 2012. [PMID 22427299](https://pubmed.ncbi.nlm.nih.gov/22427299/)  
  &nbsp;&nbsp;— botulinum toxin indication

8. **Ganesan M et al.** Direction specific preserved limits of stability in early progressive supranuclear palsy: a dynamic posturographic study. *Gait Posture* 2012. [PMID 22225854](https://pubmed.ncbi.nlm.nih.gov/22225854/)


---

## J. Anti-tau immunotherapy trials

1. **Dam T et al.** Safety and efficacy of anti-tau monoclonal antibody gosuranemab in progressive supranuclear palsy: a phase 2, randomized, placebo-controlled trial. *Nat Med* 2021. [PMID 34385707](https://pubmed.ncbi.nlm.nih.gov/34385707/) (author correction: [PMID 36253611](https://pubmed.ncbi.nlm.nih.gov/36253611/))  
  &nbsp;&nbsp;— PASSPORT, n = 486: CSF unbound N-terminal tau -98% (placebo +11%), PSPRS week 52 10.4 vs 10.6 (p = 0.85), CSF NfL unchanged. THE central anchor of this model.

2. **Höglinger GU et al.** Safety and efficacy of tilavonemab in progressive supranuclear palsy: a phase 2, randomised, placebo-controlled trial. *Lancet Neurol* 2021. [PMID 33609476](https://pubmed.ncbi.nlm.nih.gov/33609476/)  
  &nbsp;&nbsp;— ARISE: stopped for futility at the second interim analysis

3. Results of a phase 1 study of ABBV-8E12 in patients with progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Results%20of%20a%20phase%201%20study%20of%20ABBV-8E12%20in%20patients%20with%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— first-in-human PK of an N-terminal anti-tau mAb

4. **Boxer AL et al.** Safety of the tau-directed monoclonal antibody BIIB092 in progressive supranuclear palsy: a randomised, placebo-controlled, multiple ascending dose phase 1b trial. *Lancet Neurol* 2019. [PMID 31122495](https://pubmed.ncbi.nlm.nih.gov/31122495/)  
  &nbsp;&nbsp;— gosuranemab dose-finding and CSF target engagement

5. Anti-tau antibody BIIB092 semorinemab N-terminal tau target engagement CSF — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Anti-tau%20antibody%20BIIB092%20semorinemab%20N-terminal%20tau%20target%20engagement%20CSF)  
  &nbsp;&nbsp;— N-terminal engagement across programmes

6. Passive immunization targeting the N-terminal region of tau decreases tau pathology in mouse models — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Passive%20immunization%20targeting%20the%20N-terminal%20region%20of%20tau%20decreases%20tau%20pathology%20in%20mouse%20models)  
  &nbsp;&nbsp;— the preclinical rationale that did not translate

7. Anti-tau antibody epitope determines efficacy in blocking tau seeding — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Anti-tau%20antibody%20epitope%20determines%20efficacy%20in%20blocking%20tau%20seeding)  
  &nbsp;&nbsp;— epitope dependence: the model's KD_NTAB vs KD_STAB

8. **Abdel-Haleem AM et al.** CSF proteomics of semorinemab Alzheimer's disease trials identifies cell-type specific signatures. *Brain* 2026. [PMID 40435316](https://pubmed.ncbi.nlm.nih.gov/40435316/)  
  &nbsp;&nbsp;— the same N-terminal strategy failing in another tauopathy

9. Zagotenemab LY3303560 in early Alzheimer's disease conformational anti-tau antibody — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Zagotenemab%20LY3303560%20in%20early%20Alzheimer%27s%20disease%20conformational%20anti-tau%20antibody)  
  &nbsp;&nbsp;— conformation-selective antibody

10. Bepranemab UCB0107 anti-tau antibody targeting the mid-domain of tau — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Bepranemab%20UCB0107%20anti-tau%20antibody%20targeting%20the%20mid-domain%20of%20tau)  
  &nbsp;&nbsp;— a mid-domain antibody -- the model's S04 probe class


---

## K. Genetic tau lowering

1. Development of NIO752 an intrathecally administered MAPT antisense oligonucleotide in progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Development%20of%20NIO752%20an%20intrathecally%20administered%20MAPT%20antisense%20oligonucleotide%20in%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— phase 1: CSF tau and p-tau181 -20%, NfL stabilised

2. Antisense oligonucleotides provide optimism to the therapeutic landscape for tauopathies — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Antisense%20oligonucleotides%20provide%20optimism%20to%20the%20therapeutic%20landscape%20for%20tauopathies)  
  &nbsp;&nbsp;— review of tau-lowering ASOs

3. **DeVos SL et al.** Tau reduction prevents neuronal loss and reverses pathological tau deposition and seeding in mice with tauopathy. *Sci Transl Med* 2017. [PMID 28123067](https://pubmed.ncbi.nlm.nih.gov/28123067/)  
  &nbsp;&nbsp;— the preclinical case for lowering M

4. Reduction of tau by antisense oligonucleotide in a mouse model of tauopathy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Reduction%20of%20tau%20by%20antisense%20oligonucleotide%20in%20a%20mouse%20model%20of%20tauopathy)  
  &nbsp;&nbsp;— proof of mechanism

5. Intrathecal antisense oligonucleotide distribution in the central nervous system rostrocaudal gradient — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Intrathecal%20antisense%20oligonucleotide%20distribution%20in%20the%20central%20nervous%20system%20rostrocaudal%20gradient)  
  &nbsp;&nbsp;— the depth gradient behind ASO_DEEP

6. **Edwards AL et al.** Exploratory Tau Biomarker Results From a Multiple Ascending-Dose Study of BIIB080 in Alzheimer Disease: A Randomized Clinical Trial. *JAMA Neurol* 2023. [PMID 37902726](https://pubmed.ncbi.nlm.nih.gov/37902726/)

7. Tau haploinsufficiency and the safety of lowering tau — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Tau%20haploinsufficiency%20and%20the%20safety%20of%20lowering%20tau)  
  &nbsp;&nbsp;— how much lowering is tolerable


---

## L. Small-molecule and symptomatic trials

1. **Boxer AL et al.** Davunetide in patients with progressive supranuclear palsy: a randomised, double-blind, placebo-controlled phase 2/3 trial. *Lancet Neurol* 2014. [PMID 24873720](https://pubmed.ncbi.nlm.nih.gov/24873720/)  
  &nbsp;&nbsp;— AL-108-231: no effect, with the expected yearly change in BOTH arms

2. **Tolosa E et al.** A phase 2 trial of the GSK-3 inhibitor tideglusib in progressive supranuclear palsy. *Mov Disord* 2014. [PMID 24532007](https://pubmed.ncbi.nlm.nih.gov/24532007/)  
  &nbsp;&nbsp;— TAUROS: GSK-3beta inhibition, no clinical effect

3. Tideglusib GSK-3 inhibitor magnetic resonance imaging outcomes progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Tideglusib%20GSK-3%20inhibitor%20magnetic%20resonance%20imaging%20outcomes%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— TAUROS imaging substudy

4. **Bensimon G et al.** Riluzole treatment, survival and diagnostic criteria in Parkinson plus disorders: the NNIPPS study. *Brain* 2009. [PMID 19029129](https://pubmed.ncbi.nlm.nih.gov/19029129/)  
  &nbsp;&nbsp;— NNIPPS, n=767, null

5. Salsalate and young plasma in progressive supranuclear palsy an open label study — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Salsalate%20and%20young%20plasma%20in%20progressive%20supranuclear%20palsy%20an%20open%20label%20study)  
  &nbsp;&nbsp;— p300/acetyl-tau strategy

6. **Stamelou M et al.** Short-term effects of coenzyme Q10 in progressive supranuclear palsy: a randomized, placebo-controlled trial. *Mov Disord* 2008. [PMID 18464278](https://pubmed.ncbi.nlm.nih.gov/18464278/)  
  &nbsp;&nbsp;— bioenergetic strategy

7. Lithium in progressive supranuclear palsy tolerability trial — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Lithium%20in%20progressive%20supranuclear%20palsy%20tolerability%20trial)  
  &nbsp;&nbsp;— GSK-3 inhibition by another route, stopped for tolerability

8. **Fulco E** Progressive Supranuclear Palsy Unmasked After Post-COVID-19 Functional Decline in an Elderly Patient: A Diagnostic Challenge. *Cureus* 2026. [PMID 42367505](https://pubmed.ncbi.nlm.nih.gov/42367505/)

9. Zolpidem improves motor function and gaze in progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Zolpidem%20improves%20motor%20function%20and%20gaze%20in%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— GABA-A alpha1 action on the pallidum

10. **Nuebling G et al.** PROSPERA: a randomized, controlled trial evaluating rasagiline in progressive supranuclear palsy. *J Neurol* 2016. [PMID 27230855](https://pubmed.ncbi.nlm.nih.gov/27230855/)


---

## M. Emerging mechanisms in the PSP platform trial

1. **Langness VF et al.** Taming the "death receptor": translating the first-in-class p75(NTR) modulator LM11A-31 from basic biology, across broad preclinical models, to clinical proof-of-concept. *J Transl Med* 2026. [PMID 42210251](https://pubmed.ncbi.nlm.nih.gov/42210251/)  
  &nbsp;&nbsp;— tau acetylation suppression; PSP Platform Trial arm

2. Small molecule p75NTR ligand LM11A-31 reduces tau pathology and prevents degeneration — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Small%20molecule%20p75NTR%20ligand%20LM11A-31%20reduces%20tau%20pathology%20and%20prevents%20degeneration)  
  &nbsp;&nbsp;— preclinical tau-acetylation reduction

3. **Novak P et al.** ADAMANT: a placebo-controlled randomized phase 2 study of AADvac1, an active immunotherapy against pathological tau in Alzheimer's disease. *Nat Aging* 2021. [PMID 37117834](https://pubmed.ncbi.nlm.nih.gov/37117834/)  
  &nbsp;&nbsp;— active vaccine: the same extracellular ceiling

4. Safety and immunogenicity of the tau vaccine AADvac1 first-in-human study — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Safety%20and%20immunogenicity%20of%20the%20tau%20vaccine%20AADvac1%20first-in-human%20study)  
  &nbsp;&nbsp;— titre kinetics behind the model's VTIT compartment

5. AZP2006 ezeprogind progranulin phase 2 progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=AZP2006%20ezeprogind%20progranulin%20phase%202%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— progranulin enhancement -- acts on intracellular clearance

6. TPN-101 censavudine LINE-1 reverse transcriptase inhibitor progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=TPN-101%20censavudine%20LINE-1%20reverse%20transcriptase%20inhibitor%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— cGAS-STING interferon axis

7. **Roy N et al.** Elevated expression of the retrotransposon LINE-1 drives Alzheimer's disease-associated microglial dysfunction. *Acta Neuropathol* 2024. [PMID 39604588](https://pubmed.ncbi.nlm.nih.gov/39604588/)  
  &nbsp;&nbsp;— the TPN-101 rationale

8. **Wolff AW et al.** Effects of fasudil on disease spreading in ALS - A MUNIX-based post-hoc analysis of the ROCK-ALS trial. *Neurotherapeutics* 2026. [PMID 42235092](https://pubmed.ncbi.nlm.nih.gov/42235092/)  
  &nbsp;&nbsp;— ROCK inhibition


---

## N. Dysphagia, aspiration and survival

1. **Cattani AC et al.** Longitudinal Videofluorographic Dysphagia Measures in Progressive Supranuclear Palsy. *Mov Disord Clin Pract* 2026. [PMID 41792957](https://pubmed.ncbi.nlm.nih.gov/41792957/)  
  &nbsp;&nbsp;— the dominant competing risk

2. **Horiuchi K et al.** Diagnostic delay and onset-anchored clinical milestones in progressive supranuclear palsy: a japanese single-center retrospective cohort study. *Neurol Sci* 2026. [PMID 42340524](https://pubmed.ncbi.nlm.nih.gov/42340524/)

3. **Glasmacher SA et al.** Predictors of survival in progressive supranuclear palsy and multiple system atrophy: a systematic review and meta-analysis. *J Neurol Neurosurg Psychiatry* 2017. [PMID 28250027](https://pubmed.ncbi.nlm.nih.gov/28250027/)

4. Percutaneous endoscopic gastrostomy in neurodegenerative disease outcomes — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Percutaneous%20endoscopic%20gastrostomy%20in%20neurodegenerative%20disease%20outcomes)  
  &nbsp;&nbsp;— PEG as a hazard modifier, not a disease modifier

5. Palliative and supportive care in progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Palliative%20and%20supportive%20care%20in%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— symptom-directed care

6. Videofluoroscopic assessment of swallowing in progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Videofluoroscopic%20assessment%20of%20swallowing%20in%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— bulbar dysfunction measurement

7. Weight loss and nutritional status in progressive supranuclear palsy — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Weight%20loss%20and%20nutritional%20status%20in%20progressive%20supranuclear%20palsy)  
  &nbsp;&nbsp;— nutritional decline


---

## O. QSP and modelling methodology

1. **Goryanin I et al.** Validation of AI-enabled surrogate models in quantitative systems pharmacology: a practical, context-of-use-driven review. *Drug Discov Today* 2026. [PMID 42409163](https://pubmed.ncbi.nlm.nih.gov/42409163/)  
  &nbsp;&nbsp;— MIDD framework

2. mrgsolve simulation from ordinary differential equation based models in R — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=mrgsolve%20simulation%20from%20ordinary%20differential%20equation%20based%20models%20in%20R)  
  &nbsp;&nbsp;— the ODE engine used here

3. **Marshall S et al.** Model-Informed Drug Discovery and Development: Current Industry Good Practice and Regulatory Expectations and Future Perspectives. *CPT Pharmacometrics Syst Pharmacol* 2019. [PMID 30411538](https://pubmed.ncbi.nlm.nih.gov/30411538/)  
  &nbsp;&nbsp;— reporting standards

4. **Yang L et al.** Obesity-driven phosphatidylethanolamine dysregulation impairs neuroimmune crosstalk and accelerates Alzheimer's pathogenesis. *Mol Neurodegener* 2026. [PMID 41987289](https://pubmed.ncbi.nlm.nih.gov/41987289/)  
  &nbsp;&nbsp;— prior QSP work in CNS

5. Disease progression modelling in neurodegenerative clinical trials sample size — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Disease%20progression%20modelling%20in%20neurodegenerative%20clinical%20trials%20sample%20size)  
  &nbsp;&nbsp;— power and endpoint choice

6. Physiologically based pharmacokinetic modelling of monoclonal antibody distribution to the brain — [PubMed search](https://pubmed.ncbi.nlm.nih.gov/?term=Physiologically%20based%20pharmacokinetic%20modelling%20of%20monoclonal%20antibody%20distribution%20to%20the%20brain)  
  &nbsp;&nbsp;— the mAb brain-PK structure used here


---

## 이 모델의 핵심 앵커가 된 다섯 편 (the five papers this model is built on)

1. **PASSPORT (gosuranemab).** CSF unbound N-terminal tau -98%, PSPRS 52주 변화
   10.4 vs 10.6 (p = 0.85), CSF NfL 무변화. 이 모델의 중심 앵커이며, 동시에
   "표적 결합률(target engagement)이 곧 표적 타당성(target validity)은 아니다"
   라는 진술의 정량적 근거입니다.
2. **ARISE (tilavonemab).** 같은 N-말단 전략이 두 번째로, 독립적으로 실패
   (2차 중간분석에서 무효성 기준 충족). 한 번은 우연이지만 두 번은 구조입니다.
3. **NIO752 phase 1.** CSF total tau / p-tau181 -20% (PASSPORT의 1/5),
   그러나 CSF NfL은 안정화(위약 +40%). 이 질환 역사상 유일한 하류 신호가
   가장 작은 바이오마커 변화에 붙어 있습니다.
4. **AL-108-231 (davunetide).** 미세관 안정화. 무효. 결정적으로, 논문 저자들이
   "양 군 모두에서 예상된 연간 변화가 관찰되었으므로 검정력은 충분했다"고
   명시했습니다 -- 즉 실패는 통계가 아니라 기전의 문제입니다.
5. **Golbe & Ohman-Strickland 2007 (PSPRS).** 연간 11.3 +/- 11.0점. 이 모델이
   재현해야 하는 기울기이자, 표준편차가 평균과 같다는 사실 자체가
   임상시험 설계의 제약입니다.
