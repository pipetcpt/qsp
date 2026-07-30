# Paediatric Growth Hormone Deficiency (GHD) — References

79 PubMed-verified references supporting the mechanistic map
([`ghd_qsp_model.dot`](ghd_qsp_model.dot)), the 42-ODE mrgsolve model
([`ghd_mrgsolve_model.R`](ghd_mrgsolve_model.R)) and the Shiny dashboard
([`ghd_shiny_app.R`](ghd_shiny_app.R)).

**Every PMID below was checked against the NCBI E-utilities API and the title,
journal, year, volume and pages are reproduced from the PubMed record itself.**
Papers whose PMID could not be confirmed were dropped rather than guessed.

## How the literature maps onto the model

| Model element | Sections below |
|---|---|
| Hypothalamic GHRH/somatostatin/ghrelin oscillator | 2 |
| Somatotroph mass, aetiology, genetics | 1, 16 |
| GHR → JAK2 → STAT5b, SOCS2 brake, receptor down-regulation | 3 |
| IGF-1 / IGFBP-3 / ALS ternary complex, free vs total IGF-1, PAPP-A2 | 4 |
| Growth plate: resting-zone reserve, senescence, catch-up growth | 5 |
| Oestrogen entering twice with opposite sign (spurt vs fusion) | 6 |
| Bone mineral accrual and SCFE | 7 |
| Glucocorticoid brake on the growth plate | 8 |
| Daily somatropin dose–response, IGF-1-based titration, adherence | 9 |
| Weekly long-acting constructs and the within-week IGF-1 swing | 10 |
| Mecasermin (rhIGF-1) bypassing the GH receptor | 11 |
| Provocation testing, IGF-1/IGFBP-3 biomarkers | 12 |
| Free-T4 fall and adrenal unmasking on starting rhGH | 13 |
| GnRH analogue and aromatase inhibitor as "buy growing time" adjuvants | 14 |
| GH-induced insulin resistance, lipolysis, body composition | 15 |
| Long-term safety, the IGF-1 SDS ceiling, registry cohorts | 17 |
| Stopping criteria, retesting, adult GHD | 18 |

---

### 1. 지침·총설 (Guidelines, Consensus Statements & Major Reviews)

- Growth Hormone Research Society **Consensus guidelines for the diagnosis and treatment of growth hormone (GH) deficiency in childhood and adolescence: summary statement of the GH Research Society. GH Research Society.** *J Clin Endocrinol Metab 2000;85:3990-3.* [PMID 11095419](https://pubmed.ncbi.nlm.nih.gov/11095419/)
- Grimberg A et al. **Guidelines for Growth Hormone and Insulin-Like Growth Factor-I Treatment in Children and Adolescents: Growth Hormone Deficiency, Idiopathic Short Stature, and Primary Insulin-Like Growth Factor-I Deficiency.** *Horm Res Paediatr 2016;86:361-397.* [PMID 27884013](https://pubmed.ncbi.nlm.nih.gov/27884013/)
- Collett-Solberg PF et al. **Diagnosis, Genetics, and Therapy of Short Stature in Children: A Growth Hormone Research Society International Perspective.** *Horm Res Paediatr 2019;92:1-14.* [PMID 31514194](https://pubmed.ncbi.nlm.nih.gov/31514194/)
- Alatzoglou KS et al. **Isolated growth hormone deficiency (GHD) in childhood and adolescence: recent advances.** *Endocr Rev 2014;35:376-432.* [PMID 24450934](https://pubmed.ncbi.nlm.nih.gov/24450934/)
- Alatzoglou KS & Dattani MT **Phenotype-genotype correlations in congenital isolated growth hormone deficiency (IGHD).** *Indian J Pediatr 2012;79:99-106.* [PMID 22139958](https://pubmed.ncbi.nlm.nih.gov/22139958/)
- Allen DB et al. **GH safety workshop position paper: a critical appraisal of recombinant human GH therapy in children and adults.** *Eur J Endocrinol 2016;174:P1-9.* [PMID 26563978](https://pubmed.ncbi.nlm.nih.gov/26563978/)

### 2. GH 분비의 신경내분비 조절 (Neuroendocrine Control of GH Secretion)

- Giustina A & Veldhuis JD **Pathophysiology of the neuroregulation of growth hormone secretion in experimental animals and the human.** *Endocr Rev 1998;19:717-97.* [PMID 9861545](https://pubmed.ncbi.nlm.nih.gov/9861545/)
- Müller EE et al. **Neuroendocrine control of growth hormone secretion.** *Physiol Rev 1999;79:511-607.* [PMID 10221989](https://pubmed.ncbi.nlm.nih.gov/10221989/)

### 3. GHR–JAK2–STAT5b 신호전달과 SOCS2 브레이크 (Transduction & Feedback)

- Kofoed EM et al. **Growth hormone insensitivity associated with a STAT5b mutation.** *N Engl J Med 2003;349:1139-47.* [PMID 13679528](https://pubmed.ncbi.nlm.nih.gov/13679528/)
- Hwa V et al. **Severe growth hormone insensitivity resulting from total absence of signal transducer and activator of transcription 5b.** *J Clin Endocrinol Metab 2005;90:4260-6.* [PMID 15827093](https://pubmed.ncbi.nlm.nih.gov/15827093/)
- Fang P et al. **A mutant signal transducer and activator of transcription 5b, associated with growth hormone insensitivity and insulin-like growth factor-I deficiency, cannot function as a signal transducer or transcription factor.** *J Clin Endocrinol Metab 2006;91:1526-34.* [PMID 16464942](https://pubmed.ncbi.nlm.nih.gov/16464942/)
- Hwa V et al. **Genetic causes of growth hormone insensitivity beyond GHR.** *Rev Endocr Metab Disord 2021;22:43-58.* [PMID 33029712](https://pubmed.ncbi.nlm.nih.gov/33029712/)
- Metcalf D et al. **Gigantism in mice lacking suppressor of cytokine signalling-2.** *Nature 2000;405:1069-73.* [PMID 10890450](https://pubmed.ncbi.nlm.nih.gov/10890450/)
- Greenhalgh CJ & Alexander WS **Suppressors of cytokine signalling and regulation of growth hormone action.** *Growth Horm IGF Res 2004;14:200-6.* [PMID 15125881](https://pubmed.ncbi.nlm.nih.gov/15125881/)
- Li K et al. **SOCS2 regulation of growth hormone signaling requires a canonical interaction with phosphotyrosine.** *Biosci Rep 2022;42.* [PMID 36398696](https://pubmed.ncbi.nlm.nih.gov/36398696/)

### 4. GH 저항성과 IGF 시스템 — 결합 상태가 생물활성을 정한다 (GH Insensitivity & the IGF System)

- Laron Z **Laron syndrome (primary growth hormone resistance or insensitivity): the personal experience 1958-2003.** *J Clin Endocrinol Metab 2004;89:1031-44.* [PMID 15001582](https://pubmed.ncbi.nlm.nih.gov/15001582/)
- Domené HM et al. **Deficiency of the circulating insulin-like growth factor system associated with inactivation of the acid-labile subunit gene.** *N Engl J Med 2004;350:570-7.* [PMID 14762184](https://pubmed.ncbi.nlm.nih.gov/14762184/)
- Boisclair YR et al. **The acid-labile subunit (ALS) of the 150 kDa IGF-binding protein complex: an important but forgotten component of the circulating IGF system.** *J Endocrinol 2001;170:63-70.* [PMID 11431138](https://pubmed.ncbi.nlm.nih.gov/11431138/)
- Juul A **Serum levels of insulin-like growth factor I and its binding proteins in health and disease.** *Growth Horm IGF Res 2003;13:113-70.* [PMID 12914749](https://pubmed.ncbi.nlm.nih.gov/12914749/)
- Dauber A et al. **Mutations in pregnancy-associated plasma protein A2 cause short stature due to low IGF-I availability.** *EMBO Mol Med 2016;8:363-74.* [PMID 26902202](https://pubmed.ncbi.nlm.nih.gov/26902202/)
- Blum WF & Ranke MB **Use of insulin-like growth factor-binding protein 3 for the evaluation of growth disorders.** *Horm Res 1990;33 Suppl 4:31-7.* [PMID 1700965](https://pubmed.ncbi.nlm.nih.gov/1700965/)
- Campos VC et al. **IGF-I bioavailability in congenital isolated growth hormone deficiency.** *Eur J Endocrinol 2026;194:136-145.* [PMID 41528724](https://pubmed.ncbi.nlm.nih.gov/41528724/)
- Domené HM et al. **Normal growth spurt and final height despite low levels of all forms of circulating insulin-like growth factor-I in a patient with acid-labile subunit deficiency.** *Horm Res 2007;67:243-9.* [PMID 17213728](https://pubmed.ncbi.nlm.nih.gov/17213728/)

### 5. 성장판 생물학 — 증식 예비능과 노화 (Growth-Plate Biology: Reserve & Senescence)

- Kronenberg HM **Developmental regulation of the growth plate.** *Nature 2003;423:332-6.* [PMID 12748651](https://pubmed.ncbi.nlm.nih.gov/12748651/)
- Nilsson O et al. **Endocrine regulation of the growth plate.** *Horm Res 2005;64:157-65.* [PMID 16205094](https://pubmed.ncbi.nlm.nih.gov/16205094/)
- Wit JM & Camacho-Hübner C **Endocrine regulation of longitudinal bone growth.** *Endocr Dev 2011;21:30-41.* [PMID 21865752](https://pubmed.ncbi.nlm.nih.gov/21865752/)
- Lui JC et al. **Growth plate senescence and catch-up growth.** *Endocr Dev 2011;21:23-29.* [PMID 21865751](https://pubmed.ncbi.nlm.nih.gov/21865751/)
- Schrier L et al. **Depletion of resting zone chondrocytes during growth plate senescence.** *J Endocrinol 2006;189:27-36.* [PMID 16614378](https://pubmed.ncbi.nlm.nih.gov/16614378/)
- Nilsson O et al. **Growth plate senescence is associated with loss of DNA methylation.** *J Endocrinol 2005;186:241-9.* [PMID 16002553](https://pubmed.ncbi.nlm.nih.gov/16002553/)
- Nilsson O & Baron J **Impact of growth plate senescence on catch-up growth and epiphyseal fusion.** *Pediatr Nephrol 2005;20:319-22.* [PMID 15723267](https://pubmed.ncbi.nlm.nih.gov/15723267/)
- Green H et al. **A dual effector theory of growth-hormone action.** *Differentiation 1985;29:195-8.* [PMID 3908201](https://pubmed.ncbi.nlm.nih.gov/3908201/)
- Wang J et al. **Igf1 promotes longitudinal bone growth by insulin-like actions augmenting chondrocyte hypertrophy.** *FASEB J 1999;13:1985-90.* [PMID 10544181](https://pubmed.ncbi.nlm.nih.gov/10544181/)
- Vortkamp A et al. **Regulation of rate of cartilage differentiation by Indian hedgehog and PTH-related protein.** *Science 1996;273:613-22.* [PMID 8662546](https://pubmed.ncbi.nlm.nih.gov/8662546/)
- Deng C et al. **Fibroblast growth factor receptor 3 is a negative regulator of bone growth.** *Cell 1996;84:911-21.* [PMID 8601314](https://pubmed.ncbi.nlm.nih.gov/8601314/)
- Bartels CF et al. **Mutations in the transmembrane natriuretic peptide receptor NPR-B impair skeletal growth and cause acromesomelic dysplasia, type Maroteaux.** *Am J Hum Genet 2004;75:27-34.* [PMID 15146390](https://pubmed.ncbi.nlm.nih.gov/15146390/)
- Marino R et al. **Catch-up growth after hypothyroidism is caused by delayed growth plate senescence.** *Endocrinology 2008;149:1820-8.* [PMID 18174286](https://pubmed.ncbi.nlm.nih.gov/18174286/)

### 6. 오에스트로겐과 골단 폐쇄 (Oestrogen and Epiphyseal Fusion)

- Smith EP et al. **Estrogen resistance caused by a mutation in the estrogen-receptor gene in a man.** *N Engl J Med 1994;331:1056-61.* [PMID 8090165](https://pubmed.ncbi.nlm.nih.gov/8090165/)
- Morishima A et al. **Aromatase deficiency in male and female siblings caused by a novel mutation and the physiological role of estrogens.** *J Clin Endocrinol Metab 1995;80:3689-98.* [PMID 8530621](https://pubmed.ncbi.nlm.nih.gov/8530621/)
- Simm PJ et al. **Estrogens and growth.** *Pediatr Endocrinol Rev 2008;6:32-41.* [PMID 18806723](https://pubmed.ncbi.nlm.nih.gov/18806723/)

### 7. 뼈·미네랄 (Bone and Mineral)

- Ohlsson C et al. **Growth hormone and bone.** *Endocr Rev 1998;19:55-79.* [PMID 9494780](https://pubmed.ncbi.nlm.nih.gov/9494780/)
- Mittal M et al. **The Effect of Human Growth Hormone Treatment on the Development of Slipped Capital Femoral Epiphysis: A Cohort Analysis With 6 Years of Follow-up.** *J Pediatr Orthop 2024;44:e344-e350.* [PMID 38225906](https://pubmed.ncbi.nlm.nih.gov/38225906/)

### 8. 글루코코르티코이드에 의한 성장 억제 (Glucocorticoid Growth Suppression)

- Hochberg Z **Mechanisms of steroid impairment of growth.** *Horm Res 2002;58 Suppl 1:33-8.* [PMID 12373012](https://pubmed.ncbi.nlm.nih.gov/12373012/)
- Chrysis D et al. **Dexamethasone induces apoptosis in proliferative chondrocytes through activation of caspases and suppression of the Akt-phosphatidylinositol 3'-kinase signaling pathway.** *Endocrinology 2005;146:1391-7.* [PMID 15576458](https://pubmed.ncbi.nlm.nih.gov/15576458/)
- Klaus G et al. **Suppression of growth plate chondrocyte proliferation by corticosteroids.** *Pediatr Nephrol 2000;14:612-5.* [PMID 10912528](https://pubmed.ncbi.nlm.nih.gov/10912528/)
- Lui JC & Baron J **Effects of glucocorticoids on the growth plate.** *Endocr Dev 2011;20:187-193.* [PMID 21164272](https://pubmed.ncbi.nlm.nih.gov/21164272/)

### 9. 일일 소마트로핀 — 반응 예측과 용량 (Daily Somatropin: Response Prediction & Dosing)

- Ranke MB et al. **Derivation and validation of a mathematical model for predicting the response to exogenous recombinant human growth hormone (GH) in prepubertal children with idiopathic GH deficiency. KIGS International Board. Kabi Pharmacia International Growth Study.** *J Clin Endocrinol Metab 1999;84:1174-83.* [PMID 10199749](https://pubmed.ncbi.nlm.nih.gov/10199749/)
- Ranke MB et al. **Observed and predicted growth responses in prepubertal children with growth disorders: guidance of growth hormone treatment by empirical variables.** *J Clin Endocrinol Metab 2010;95:1229-37.* [PMID 20097713](https://pubmed.ncbi.nlm.nih.gov/20097713/)
- Cohen P et al. **Dose-sparing and safety-enhancing effects of an IGF-I-based dosing regimen in short children treated with growth hormone in a 2-year randomized controlled trial: therapeutic and pharmacoeconomic considerations.** *Clin Endocrinol (Oxf) 2014;81:71-6.* [PMID 24428305](https://pubmed.ncbi.nlm.nih.gov/24428305/)
- Cohen P et al. **Variable degree of growth hormone (GH) and insulin-like growth factor (IGF) sensitivity in children with idiopathic short stature compared with GH-deficient patients: evidence from an IGF-based dosing study of short children.** *J Clin Endocrinol Metab 2010;95:2089-98.* [PMID 20207829](https://pubmed.ncbi.nlm.nih.gov/20207829/)
- Fisher BG & Acerini CL **Understanding the growth hormone therapy adherence paradigm: a systematic review.** *Horm Res Paediatr 2013;79:189-96.* [PMID 23635797](https://pubmed.ncbi.nlm.nih.gov/23635797/)

### 10. 주 1회 지속형 GH (Weekly Long-Acting GH)

- Christiansen JS et al. **Growth Hormone Research Society perspective on the development of long-acting growth hormone preparations.** *Eur J Endocrinol 2016;174:C1-8.* [PMID 27009113](https://pubmed.ncbi.nlm.nih.gov/27009113/)
- Thornton PS et al. **Weekly Lonapegsomatropin in Treatment-Naïve Children With Growth Hormone Deficiency: The Phase 3 heiGHt Trial.** *J Clin Endocrinol Metab 2021;106:3184-3195.* [PMID 34272849](https://pubmed.ncbi.nlm.nih.gov/34272849/)
- Deal CL et al. **Efficacy and Safety of Weekly Somatrogon vs Daily Somatropin in Children With Growth Hormone Deficiency: A Phase 3 Study.** *J Clin Endocrinol Metab 2022;107:e2717-e2728.* [PMID 35405011](https://pubmed.ncbi.nlm.nih.gov/35405011/)
- Miller BS et al. **Effective GH Replacement With Somapacitan in Children With GHD: REAL4 2-year Results and After Switch From Daily GH.** *J Clin Endocrinol Metab 2023;108:3090-3099.* [PMID 37406251](https://pubmed.ncbi.nlm.nih.gov/37406251/)
- Maniatis AK et al. **Safety and Efficacy of Lonapegsomatropin in Children With Growth Hormone Deficiency: enliGHten Trial 2-Year Results.** *J Clin Endocrinol Metab 2022;107:e2680-e2689.* [PMID 35428884](https://pubmed.ncbi.nlm.nih.gov/35428884/)
- Maniatis AK et al. **Switching to Weekly Lonapegsomatropin from Daily Somatropin in Children with Growth Hormone Deficiency: The fliGHt Trial.** *Horm Res Paediatr 2022;95:233-243.* [PMID 35263755](https://pubmed.ncbi.nlm.nih.gov/35263755/)
- Backeljauw PF et al. **Growth Response to Weekly Somapacitan Therapy in Children With GH Deficiency Is Related to GH Thresholds in GH Stimulation Testing.** *J Endocr Soc 2025;9:bvaf038.* [PMID 40104568](https://pubmed.ncbi.nlm.nih.gov/40104568/)

### 11. rhIGF-1 (메카세르민) (Recombinant IGF-1 for GH Insensitivity)

- Chernausek SD et al. **Long-term treatment with recombinant insulin-like growth factor (IGF)-I in children with severe IGF-I deficiency due to growth hormone insensitivity.** *J Clin Endocrinol Metab 2007;92:902-10.* [PMID 17192294](https://pubmed.ncbi.nlm.nih.gov/17192294/)
- Bang P et al. **Effectiveness and safety of rhIGF1 therapy in patients with or without Laron syndrome.** *Eur J Endocrinol 2021;184:267-276.* [PMID 33434161](https://pubmed.ncbi.nlm.nih.gov/33434161/)
- Muthuvel G et al. **Recombinant Human Insulin-Like Growth Factor-1 Treatment of Severe Growth Failure in Three Siblings with STAT5B Deficiency.** *Horm Res Paediatr 2024;97:195-202.* [PMID 37586336](https://pubmed.ncbi.nlm.nih.gov/37586336/)

### 12. 진단 (Diagnosis: Provocation Testing & Biomarkers)

- Sizonenko PC et al. **Diagnosis and management of growth hormone deficiency in childhood and adolescence. Part 1: diagnosis of growth hormone deficiency.** *Growth Horm IGF Res 2001;11:137-65.* [PMID 11735230](https://pubmed.ncbi.nlm.nih.gov/11735230/)
- Garcia JM et al. **Macimorelin as a Diagnostic Test for Adult GH Deficiency.** *J Clin Endocrinol Metab 2018;103:3083-3093.* [PMID 29860473](https://pubmed.ncbi.nlm.nih.gov/29860473/)
- Smyczyńska J et al. **Significance of Direct Confirmation of Growth Hormone Insensitivity for the Diagnosis of Primary IGF-I Deficiency.** *J Clin Med 2020;9.* [PMID 31963242](https://pubmed.ncbi.nlm.nih.gov/31963242/)

### 13. 다른 뇌하수체 축과의 상호작용 (Cross-Axis Interactions: Thyroid & Adrenal)

- Porretti S et al. **Recombinant human GH replacement therapy and thyroid function in a large group of adult GH-deficient patients: when does L-T(4) therapy become mandatory?.** *J Clin Endocrinol Metab 2002;87:2042-5.* [PMID 11994338](https://pubmed.ncbi.nlm.nih.gov/11994338/)
- Giavoli C et al. **Recombinant hGH replacement therapy and the hypothalamus-pituitary-thyroid axis in children with GH deficiency: when should we be concerned about the occurrence of central hypothyroidism?.** *Clin Endocrinol (Oxf) 2003;59:806-10.* [PMID 14974926](https://pubmed.ncbi.nlm.nih.gov/14974926/)
- Giavoli C et al. **Effect of recombinant human growth hormone (GH) replacement on the hypothalamic-pituitary-adrenal axis in adult GH-deficient patients.** *J Clin Endocrinol Metab 2004;89:5397-401.* [PMID 15531488](https://pubmed.ncbi.nlm.nih.gov/15531488/)

### 14. 보조 약물 — GnRH 유사체·아로마타제 억제제 (Adjuvants: GnRH Analogues & Aromatase Inhibitors)

- Carel JC et al. **Consensus statement on the use of gonadotropin-releasing hormone analogs in children.** *Pediatrics 2009;123:e752-62.* [PMID 19332438](https://pubmed.ncbi.nlm.nih.gov/19332438/)
- Hero M et al. **Inhibition of estrogen biosynthesis with a potent aromatase inhibitor increases predicted adult height in boys with idiopathic short stature: a randomized controlled trial.** *J Clin Endocrinol Metab 2005;90:6396-402.* [PMID 16189252](https://pubmed.ncbi.nlm.nih.gov/16189252/)
- Mauras N et al. **Anastrozole increases predicted adult height of short adolescent males treated with growth hormone: a randomized, placebo-controlled, multicenter trial for one to three years.** *J Clin Endocrinol Metab 2008;93:823-31.* [PMID 18165285](https://pubmed.ncbi.nlm.nih.gov/18165285/)
- Lem AJ et al. **Bone mineral density and body composition in short children born SGA during growth hormone and gonadotropin releasing hormone analog treatment.** *J Clin Endocrinol Metab 2013;98:77-86.* [PMID 23125290](https://pubmed.ncbi.nlm.nih.gov/23125290/)
- van der Steen M et al. **Insulin Sensitivity and β-Cell Function in SGA Children Treated With GH and GnRHa: Results of a Long-Term Trial.** *J Clin Endocrinol Metab 2016;101:705-13.* [PMID 26653111](https://pubmed.ncbi.nlm.nih.gov/26653111/)

### 15. 대사 효과 (Metabolic Effects of GH)

- Møller N & Jørgensen JO **Effects of growth hormone on glucose, lipid, and protein metabolism in human subjects.** *Endocr Rev 2009;30:152-77.* [PMID 19240267](https://pubmed.ncbi.nlm.nih.gov/19240267/)

### 16. 획득성 GHD (Acquired GHD: Irradiation & Tumours)

- Lövgren I et al. **The late effects of cranial irradiation in childhood on the hypothalamic-pituitary axis: a radiotherapist's perspective.** *Endocr Connect 2022;11.* [PMID 36269600](https://pubmed.ncbi.nlm.nih.gov/36269600/)
- Tselovalnikova T et al. **Prevalence of growth hormone deficiency in brain tumor survivors: a systematic review and meta-analysis.** *Endocr Oncol 2025;5:e250025.* [PMID 40641631](https://pubmed.ncbi.nlm.nih.gov/40641631/)

### 17. 안전성 (Long-Term Safety)

- Carel JC et al. **Long-term mortality after recombinant growth hormone treatment for isolated growth hormone deficiency or childhood short stature: preliminary report of the French SAGhE study.** *J Clin Endocrinol Metab 2012;97:416-25.* [PMID 22238382](https://pubmed.ncbi.nlm.nih.gov/22238382/)
- Swerdlow AJ et al. **Cancer Risks in Patients Treated With Growth Hormone in Childhood: The SAGhE European Cohort Study.** *J Clin Endocrinol Metab 2017;102:1661-1672.* [PMID 28187225](https://pubmed.ncbi.nlm.nih.gov/28187225/)
- Krasnow MD et al. **Growth hormone therapy does not impact the development of intracranial hypertension in children with Chiari malformation.** *J Pediatr Endocrinol Metab 2024;37:630-634.* [PMID 38776636](https://pubmed.ncbi.nlm.nih.gov/38776636/)

### 18. 전환기·성인 GHD (Transition and Adult GHD)

- Clayton PE et al. **Consensus statement on the management of the GH-treated adolescent in the transition to adult care.** *Eur J Endocrinol 2005;152:165-70.* [PMID 15745921](https://pubmed.ncbi.nlm.nih.gov/15745921/)
- Loche S et al. **Growth Hormone Deficiency in the Transition Age.** *Endocr Dev 2018;33:46-56.* [PMID 29886481](https://pubmed.ncbi.nlm.nih.gov/29886481/)

---

## 참고 사항 (Notes on use)

- Sections 5 and 6 are the load-bearing literature for this model's central
  claim: the growth plate has a **finite** proliferative reserve that is spent
  by cycling and by oestrogen exposure, so catch-up growth is a consumable and
  the year-on-year fall in height velocity on an unchanged dose needs no fitted
  decay term. Schrier (PMID 16614378), Nilsson (PMID 16002553, PMID 15723267),
  Lui & Baron (PMID 21865751) and Marino (PMID 18174286) are the direct
  experimental support; Marino in particular shows that catch-up growth after
  hypothyroidism happens **because** senescence was delayed while growth was
  slow, which is exactly the mechanism the model's `RZ` compartment encodes.
- Section 4 supports the "double hit" in GHD: hepatic IGF-1 synthesis falls AND
  the IGFBP-3/ALS carrier pool falls, so the free fraction rises and clearance
  of what little is made accelerates. Domené (PMID 14762184) and Boisclair
  (PMID 11431138) establish the ternary complex's role as the circulating
  reservoir; Campos (PMID 41528724) and Domené (PMID 17213728) document that
  free/bioactive IGF-1 is far better preserved than total IGF-1, which is why
  the model does **not** let circulating free IGF-1 alone drive the growth plate.
- Section 10's three phase-3 programmes (heiGHt, PMID 34272849; somatrogon,
  PMID 35405011; REAL4, PMID 37406251) are the calibration targets for the
  weekly arms. Each product's `POT_W` (bioactive equivalents per mg) was fitted
  so that its own label dose reproduces its own reported annualised height
  velocity — the model does not assume the three are interchangeable per mg.
- The height and IGF-1 SDS reference tables built into the model are smooth
  approximations for internal bookkeeping. They are **not** a validated national
  growth chart and **not** an assay-specific IGF-1 reference; do not read
  absolute SDS values out of this model into a clinical context.
