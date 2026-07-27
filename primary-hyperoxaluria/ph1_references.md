# 원발성 과옥살산뇨증 (Primary Hyperoxaluria) QSP 모델 — 참고문헌

**Primary Hyperoxaluria (PH1 / PH2 / PH3) — Reference List for the QSP Model**

---

## 이 문헌 목록을 읽는 방법 (How to read this list)

이 목록의 **모든 PMID는 NCBI E-utilities(esearch + esummary)를 통해 실제로
조회하여 확인한 것**이며, 제목·저널·연도는 PubMed가 반환한 원본 레코드에서
그대로 가져왔습니다. 기억에 의존해 작성한 인용은 한 건도 포함하지 않았습니다
(과거 세션에서 기억으로 쓴 PMID 다섯 건 중 네 건이 전혀 무관한 논문을
가리켰던 사고가 있었기 때문입니다). 링크는 모두
`https://pubmed.ncbi.nlm.nih.gov/<PMID>/` 형식입니다.

Every PMID below was resolved live through NCBI E-utilities; titles, journals
and years are PubMed's own strings, not recalled ones.

### 파라미터 근거의 정직한 구분 (Honest provenance of parameters)

QSP 모델의 신뢰성은 "어떤 숫자가 문헌에서 왔고 어떤 숫자가 추정인지"를
구분하는 데서 나옵니다. 이 모델에서 각각은 다음과 같습니다.

**문헌에 직접 근거한 값 (literature-anchored)**

| 값 | 모델에서의 위치 | 근거 섹션 |
|---|---|---|
| 정상 24시간 요중 옥살산 < 0.46 mmol/1.73 m²/일 | `D06` 앵커 | 1, 7 |
| PH1 요중 옥살산 1.0–2.5 mmol/1.73 m²/일 | `FAGT` 스캔의 표현형 범위 | 1, 3 |
| 정상 혈장 옥살산 1–6 µmol/L, ESKD-PH1 60–120 µmol/L | `POX_OUT` | 6 |
| 혈장 CaOx 용해도 한계 ≈ 30 µmol/L | `POX_CRIT` | 6 |
| 옥살산 신클리어런스 / 크레아티닌 클리어런스 ≈ 1.0–1.5 (세뇨관 분비) | `F_SEC_OX` | 6, 9 |
| ILLUMINATE-A 루마시란 요중 옥살산 −65.4 % | `EMAX_HAO1` 보정 표적 | 10 |
| 루마시란 용법 3 mg/kg 월 1회 ×3 → 분기 1회 | `lum_ev()` | 10 |
| 루마시란 혈장 t½ ≈ 5 시간, Tmax ≈ 4 시간 | `KA_LUM`, `CL_LUM` | 10, 12 |
| 네도시란 160 mg SC 월 1회 | `ned_ev()` | 11 |
| PHYOX2에서 PH2 무반응 | `D07` (구조적으로 재현) | 11 |
| G170R 동형접합에서 피리독신 반응 −30 ~ −50 % | `B6RESC` | 3 |
| 피리독신 5–20 mg/kg/일, 고용량 감각신경병증 | `PN_MGKG`, `B6_NEURO` | 3 |
| 과수분요법 2–3 L/m²/일 | `U_VOL_TGT` | 14 |
| Tiselius AP(CaOx) index 산식 | `ap_index` | 7 |
| PH1의 약 20–50 %가 25세까지 ESKD | `K_LOSS` 보정 표적 | 1, 15 |
| 하이드록시프롤린이 옥살산의 약 15 % 기여 | `J_HYP`, `FR_M_*` | 5 |
| GalNAc-siRNA의 ASGPR 매개 간세포 선택적 흡수 | `CLUP_LUM`, `ASGPR` 노드 | 12 |

**추정값 — 문헌으로 검증되지 않았으며 그렇게 표시함 (estimates, flagged as such)**

다음 파라미터에는 검증 가능한 단일 인용이 없습니다. 모델 내적 정합성(정상
상태가 정확한 고정점이 되도록)과 임상 앵커 재현을 위해 선택된 값이며,
**문헌값으로 제시하지 않습니다.**

- `PHI_GOOX` (글리옥실산 중 글리콜산 산화효소가 옥살산으로 산화하는 비율,
  0.15) — 진단 `D12`는 PHYOX2 결과를 역산해 0.077을 제시하며, 이 값 자체를
  간세포 실험으로 반증 가능한 예측으로 내놓습니다.
- `FR_LDH` / `FR_GRH` / `FR_GXSP` (세포질 글리옥실산의 분기 비율)
- `GXP_REF`, `KM_AGT` (AGT의 포화 정도 — `D03`의 비선형성 크기를 결정하며
  민감도가 가장 큰 가정)
- `CLDEP_B/S/K`, `POX_CRIT_B/S` (조직 침착 클리어런스와 역치)
- `BMAX_BONE` = 600 mmol (골 옥살산 저장 용량) — 물리적으로 가능한 상한
  근처의 값이며, 문헌의 골 옥살산 정량 범위가 넓어 확정하지 못했습니다.
- `K_RES` / `K_BURY` / `K_RES_D` (골 표면 교환 풀과 심부 결정 결합 풀의 2구획
  동역학) — 칼슘 골동역학의 표준 구조를 차용한 것이며 옥살산에 대해 직접
  측정된 값이 아닙니다.
- `K_NUC` … `K_LOSS` (결정 → NLRP3 → IL-1β → TGF-β → 섬유화 → 네프론 소실
  캐스케이드의 모든 속도상수). 개별 상수는 임의이며, 오직 **합성된 결과**
  (PH1의 ESKD 도달 연령, eGFR 기울기)만 문헌에 맞춰 보정했습니다.
- `OX_XHEP2` = 600 µmol/일 (PH2에서 드러나는 간외 옥살산 생성) — `D07`은 이
  값을 0으로 두어도 결론이 유지됨을 보이므로, 헤드라인 결과는 이 추정값에
  의존하지 않습니다.
- `IL1_BLOCK` 시나리오(S33)의 효과 크기 — 승인된 치료가 없는 기전적 예측이며,
  임상 근거가 전혀 없습니다.

---


## 1. 질환 개요·역학·진단 (Overview, Epidemiology, Diagnosis)

1. Cunneely OP (2025). *Hyperoxaluria: Diagnosis and Treatment*. Urol Clin North Am. [PMID 40610080](https://pubmed.ncbi.nlm.nih.gov/40610080/)
2. An L (2025). *Gut microbiota modulation via fecal microbiota transplantation mitigates hyperoxaluria and calcium oxalate crystal depositions induced by high oxalate diet*. Gut Microbes. [PMID 39873191](https://pubmed.ncbi.nlm.nih.gov/39873191/)
3. Michael M (2024). *Diagnosis and management of primary hyperoxalurias: best practices*. Pediatr Nephrol. [PMID 38753085](https://pubmed.ncbi.nlm.nih.gov/38753085/)
4. Liu Y (2024). *Clinical features and mutational spectrum of Chinese patients with primary hyperoxaluria type 2*. Urolithiasis. [PMID 38727838](https://pubmed.ncbi.nlm.nih.gov/38727838/)
5. Anderegg MA (2024). *Prevalence and characteristics of genetic disease in adult kidney stone formers*. Nephrol Dial Transplant. [PMID 38544324](https://pubmed.ncbi.nlm.nih.gov/38544324/)
6. Zheng Y (2024). *Infant primary hyperoxaluria type 1: A case report and literature review*. Zhong Nan Da Xue Xue Bao Yi Xue Ban. [PMID 39311781](https://pubmed.ncbi.nlm.nih.gov/39311781/)
7. Groothoff JW (2023). *Clinical practice recommendations for primary hyperoxaluria: an expert consensus statement from ERKNet and OxalEurope*. Nat Rev Nephrol. [PMID 36604599](https://pubmed.ncbi.nlm.nih.gov/36604599/)
8. Boussetta A (2023). *Primary hyperoxaluria type 1: Clinical, genetic, and evolutionary characteristics in Tunisian children*. Tunis Med. [PMID 38445424](https://pubmed.ncbi.nlm.nih.gov/38445424/)
9. Ni T (2023). *A case report of invasive infantile primary hyperoxaluria type 1 and literature review*. CEN Case Rep. [PMID 36194362](https://pubmed.ncbi.nlm.nih.gov/36194362/)
10. Demoulin N (2022). *Pathophysiology and Management of Hyperoxaluria and Oxalate Nephropathy: A Review*. Am J Kidney Dis. [PMID 34508834](https://pubmed.ncbi.nlm.nih.gov/34508834/)
11. Xie X (2022). *Primary Hyperoxaluria*. N Engl J Med. [PMID 35245013](https://pubmed.ncbi.nlm.nih.gov/35245013/)
12. Bacchetta J (2022). *Primary hyperoxaluria type 1: time for prime time?*. Clin Kidney J. [PMID 35592621](https://pubmed.ncbi.nlm.nih.gov/35592621/)
13. Fargue S (2022). *Primary hyperoxaluria type 1: pathophysiology and genetics*. Clin Kidney J. [PMID 35592619](https://pubmed.ncbi.nlm.nih.gov/35592619/)
14. Gupta A (2022). *Treatment of primary hyperoxaluria type 1*. Clin Kidney J. [PMID 35592620](https://pubmed.ncbi.nlm.nih.gov/35592620/)
15. Cai Z (2021). *Primary hyperoxaluria diagnosed after kidney transplantation: a case report and literature review*. BMC Nephrol. [PMID 34837989](https://pubmed.ncbi.nlm.nih.gov/34837989/)
16. Monet-Didailler C (2021). *[Nephrocalcinosis in children]*. Nephrol Ther. [PMID 33461896](https://pubmed.ncbi.nlm.nih.gov/33461896/)
17. Du DF (2018). *Updated Genetic Testing of Primary Hyperoxaluria Type 1 in a Chinese Population: Results from a Single Center Study and a Systematic Review*. Curr Med Sci. [PMID 30341509](https://pubmed.ncbi.nlm.nih.gov/30341509/)
18. Rootman MS (2018). *Imaging features of primary hyperoxaluria*. Clin Imaging. [PMID 30253334](https://pubmed.ncbi.nlm.nih.gov/30253334/)
19. Jiang D (2017). *Primary Hyperoxaluria*. N Engl J Med. [PMID 28402768](https://pubmed.ncbi.nlm.nih.gov/28402768/)
20. Strauss SB (2017). *Primary hyperoxaluria: spectrum of clinical and imaging findings*. Pediatr Radiol. [PMID 27844104](https://pubmed.ncbi.nlm.nih.gov/27844104/)
21. Martin-Higueras C (2017). *Molecular therapy of primary hyperoxaluria*. J Inherit Metab Dis. [PMID 28425073](https://pubmed.ncbi.nlm.nih.gov/28425073/)
22. Ben-Shalom E (2015). *Primary hyperoxalurias: diagnosis and treatment*. Pediatr Nephrol. [PMID 25519509](https://pubmed.ncbi.nlm.nih.gov/25519509/)
23. Hopp K (2015). *Phenotype-Genotype Correlations and Estimated Carrier Frequencies of Primary Hyperoxaluria*. J Am Soc Nephrol. [PMID 25644115](https://pubmed.ncbi.nlm.nih.gov/25644115/)
24. Lorenzo V (2014). *Primary hyperoxaluria*. Nefrologia. [PMID 24798559](https://pubmed.ncbi.nlm.nih.gov/24798559/)
25. Cochat P (2013). *Primary hyperoxaluria*. N Engl J Med. [PMID 23944302](https://pubmed.ncbi.nlm.nih.gov/23944302/)
26. Hoppe B (2009). *The primary hyperoxalurias*. Kidney Int. [PMID 19225556](https://pubmed.ncbi.nlm.nih.gov/19225556/)
27. Watts RW (1997). *Primary hyperoxaluria*. Contrib Nephrol. [PMID 9399058](https://pubmed.ncbi.nlm.nih.gov/9399058/)
28. Folz SJ (1991). *The peroxisome and the eye*. Surv Ophthalmol. [PMID 1710072](https://pubmed.ncbi.nlm.nih.gov/1710072/)
29. HOCKADAY TD (1964). *PRIMARY HYPEROXALURIA*. Medicine (Baltimore). [PMID 14170789](https://pubmed.ncbi.nlm.nih.gov/14170789/)

## 2. 글리옥실산 대사와 AGT 효소학 (Glyoxylate Metabolism and AGT Enzymology)

30. Zhang D (2025). *LNP-mediated in vivo base editing corrects Agxt to cure primary hyperoxaluria type 1*. Clin Transl Med. [PMID 41275431](https://pubmed.ncbi.nlm.nih.gov/41275431/)
31. Mackinnon SR (2025). *Luminescence-based complementation assay to assess target engagement and cell permeability of glycolate oxidase (HAO1) inhibitors*. Biochimie. [PMID 39151880](https://pubmed.ncbi.nlm.nih.gov/39151880/)
32. Jiang Y (2025). *Efficient and safe in vivo treatment of primary hyperoxaluria type 1 via LNP-CRISPR-Cas9-mediated glycolate oxidase disruption*. Mol Ther. [PMID 39385468](https://pubmed.ncbi.nlm.nih.gov/39385468/)
33. Das S (2024). *Inhibition of hepatic oxalate overproduction ameliorates metabolic dysfunction-associated steatohepatitis*. Nat Metab. [PMID 39333384](https://pubmed.ncbi.nlm.nih.gov/39333384/)
34. Patel YP (2024). *Glycolate oxidase-1 gene variants influence the risk of hyperoxaluria and renal stone development*. World J Urol. [PMID 38214752](https://pubmed.ncbi.nlm.nih.gov/38214752/)
35. Baltazar P (2023). *Oxalate (dys)Metabolism: Person-to-Person Variability, Kidney and Cardiometabolic Toxicity*. Genes (Basel). [PMID 37761859](https://pubmed.ncbi.nlm.nih.gov/37761859/)
36. Hoppe B (2022). *Improving Treatment Options for Primary Hyperoxaluria*. Drugs. [PMID 35779234](https://pubmed.ncbi.nlm.nih.gov/35779234/)
37. Zeng Z (2022). *HAO1-mediated oxalate metabolism promotes lung pre-metastatic niche formation by inducing neutrophil extracellular traps*. Oncogene. [PMID 35739335](https://pubmed.ncbi.nlm.nih.gov/35739335/)
38. Dejban P (2022). *New therapeutics for primary hyperoxaluria type 1*. Curr Opin Nephrol Hypertens. [PMID 35266883](https://pubmed.ncbi.nlm.nih.gov/35266883/)
39. Chen H (2022). *HAO1 negatively regulates liver macrophage activation via the NF-κB pathway in alcohol-associated liver disease*. Cell Signal. [PMID 35953025](https://pubmed.ncbi.nlm.nih.gov/35953025/)
40. Ming S (2022). *Oxalate-induced apoptosis through ERS-ROS-NF-κB signalling pathway in renal tubular epithelial cell*. Mol Med. [PMID 35922749](https://pubmed.ncbi.nlm.nih.gov/35922749/)
41. Xie J (2022). *Ferrostatin‑1 alleviates oxalate‑induced renal tubular epithelial cell injury, fibrosis and calcium oxalate stone formation by inhibiting ferroptosis*. Mol Med Rep. [PMID 35703358](https://pubmed.ncbi.nlm.nih.gov/35703358/)
42. Garrelfs SF (2021). *Lumasiran, an RNAi Therapeutic for Primary Hyperoxaluria Type 1*. N Engl J Med. [PMID 33789010](https://pubmed.ncbi.nlm.nih.gov/33789010/)
43. Scott LJ (2021). *Lumasiran: First Approval*. Drugs. [PMID 33405070](https://pubmed.ncbi.nlm.nih.gov/33405070/)
44. Wang BJ (2020). *Diet and Adaptive Evolution of Alanine-Glyoxylate Aminotransferase Mitochondrial Targeting in Birds*. Mol Biol Evol. [PMID 31702777](https://pubmed.ncbi.nlm.nih.gov/31702777/)
45. Yang X (2020). *AhR activation attenuates calcium oxalate nephrocalcinosis by diminishing M1 macrophage polarization and promoting M2 macrophage polarization*. Theranostics. [PMID 33204326](https://pubmed.ncbi.nlm.nih.gov/33204326/)
46. Sun Y (2019). *Loss of alanine-glyoxylate and serine-pyruvate aminotransferase expression accelerated the progression of hepatocellular carcinoma and predicted poor prognosis*. J Transl Med. [PMID 31771612](https://pubmed.ncbi.nlm.nih.gov/31771612/)
47. Marangella M (2016). *[The Hyperoxalurias]*. G Ital Nefrol. [PMID 27960020](https://pubmed.ncbi.nlm.nih.gov/27960020/)
48. Oppici E (2016). *Natural and Unnatural Compounds Rescue Folding Defects of Human Alanine: Glyoxylate Aminotransferase Leading to Primary Hyperoxaluria Type I*. Curr Drug Targets. [PMID 26931357](https://pubmed.ncbi.nlm.nih.gov/26931357/)
49. Dellero Y (2016). *Photorespiratory glycolate-glyoxylate metabolism*. J Exp Bot. [PMID 26994478](https://pubmed.ncbi.nlm.nih.gov/26994478/)
50. Oppici E (2015). *Liver peroxisomal alanine:glyoxylate aminotransferase and the effects of mutations associated with Primary Hyperoxaluria Type I: An overview*. Biochim Biophys Acta. [PMID 25620715](https://pubmed.ncbi.nlm.nih.gov/25620715/)
51. Rodionov RN (2014). *AGXT2: a promiscuous aminotransferase*. Trends Pharmacol Sci. [PMID 25294000](https://pubmed.ncbi.nlm.nih.gov/25294000/)
52. Frishberg Y (2014). *Mutations in HAO1 encoding glycolate oxidase cause isolated glycolic aciduria*. J Med Genet. [PMID 24996905](https://pubmed.ncbi.nlm.nih.gov/24996905/)
53. Cellini B (2011). *Human liver peroxisomal alanine:glyoxylate aminotransferase: characterization of the two allelic forms and their pathogenic variants*. Biochim Biophys Acta. [PMID 21176891](https://pubmed.ncbi.nlm.nih.gov/21176891/)
54. Ichiyama A (2011). *Studies on a unique organelle localization of a liver enzyme, serine:pyruvate (or alanine:glyoxylate) aminotransferase*. Proc Jpn Acad Ser B Phys Biol Sci. [PMID 21558762](https://pubmed.ncbi.nlm.nih.gov/21558762/)
55. Xu HW (2006). *Oxalate accumulation and regulation is independent of glycolate oxidase in rice leaves*. J Exp Bot. [PMID 16595582](https://pubmed.ncbi.nlm.nih.gov/16595582/)
56. Wanders RJ (2006). *Biochemistry of mammalian peroxisomes revisited*. Annu Rev Biochem. [PMID 16756494](https://pubmed.ncbi.nlm.nih.gov/16756494/)
57. Leth PM (2005). *Ethylene glycol poisoning*. Forensic Sci Int. [PMID 16226155](https://pubmed.ncbi.nlm.nih.gov/16226155/)
58. Danpure CJ (1997). *Variable peroxisomal and mitochondrial targeting of alanine: glyoxylate aminotransferase in mammalian evolution and disease*. Bioessays. [PMID 9136629](https://pubmed.ncbi.nlm.nih.gov/9136629/)
59. Watts RW (1992). *Alanine glyoxylate aminotransferase deficiency: biochemical and molecular genetic lessons from the study of a human disease*. Adv Enzyme Regul. [PMID 1496924](https://pubmed.ncbi.nlm.nih.gov/1496924/)
60. Kisaki T (1969). *Glycolate and glyoxylate metabolism by isolated peroxisomes or chloroplasts*. Plant Physiol. [PMID 16657053](https://pubmed.ncbi.nlm.nih.gov/16657053/)

## 3. 유전형-표현형 및 피리독신 반응성 (Genotype-Phenotype, Pyridoxine Responsiveness)

61. Zhu X (2024). *Mutation Characteristics of Primary Hyperoxaluria in the Chinese Population and Current International Diagnosis and Treatment Status*. Kidney Dis (Basel). [PMID 39131880](https://pubmed.ncbi.nlm.nih.gov/39131880/)
62. Bhasin B (2015). *Primary and secondary hyperoxaluria: Understanding the enigma*. World J Nephrol. [PMID 25949937](https://pubmed.ncbi.nlm.nih.gov/25949937/)
63. Hoppe B (2012). *An update on primary hyperoxaluria*. Nat Rev Nephrol. [PMID 22688746](https://pubmed.ncbi.nlm.nih.gov/22688746/)
64. Ortiz-Alvarado O (2011). *Pyridoxine and dietary counseling for the management of idiopathic hyperoxaluria in stone-forming patients*. Urology. [PMID 21334732](https://pubmed.ncbi.nlm.nih.gov/21334732/)
65. Bobrowski AE (2008). *The primary hyperoxalurias*. Semin Nephrol. [PMID 18359396](https://pubmed.ncbi.nlm.nih.gov/18359396/)
66. Bobrowski AE (2006). *Hyperoxaluria and systemic oxalosis: current therapy and future directions*. Expert Opin Pharmacother. [PMID 17020415](https://pubmed.ncbi.nlm.nih.gov/17020415/)
67. Monico CG (2005). *Implications of genotype and enzyme phenotype in pyridoxine response of patients with type I primary hyperoxaluria*. Am J Nephrol. [PMID 15849466](https://pubmed.ncbi.nlm.nih.gov/15849466/)
68. Asplin JR (2002). *Hyperoxaluric calcium nephrolithiasis*. Endocrinol Metab Clin North Am. [PMID 12474639](https://pubmed.ncbi.nlm.nih.gov/12474639/)
69. Milliner DS (1994). *Results of long-term treatment with orthophosphate and pyridoxine in patients with primary hyperoxaluria*. N Engl J Med. [PMID 7969325](https://pubmed.ncbi.nlm.nih.gov/7969325/)
70. Edwards P (1991). *Metabolism of pyridoxine in mild metabolic hyperoxaluria and primary hyperoxaluria (type 1)*. Urol Int. [PMID 1771698](https://pubmed.ncbi.nlm.nih.gov/1771698/)
71. Latta K (1990). *Primary hyperoxaluria type I*. Eur J Pediatr. [PMID 2189732](https://pubmed.ncbi.nlm.nih.gov/2189732/)

## 4. PH2 / PH3 — 다른 효소, 다른 약리 (PH2 / PH3)

72. Yan X (2025). *Identification of a novel GRHPR mutation in primary hyperoxaluria type 2 and establishment of patient-derived iPSC line*. Hum Cell. [PMID 39757298](https://pubmed.ncbi.nlm.nih.gov/39757298/)
73. Hoppe B (2025). *Effective Newborn Screening for Type 1 and 3 Primary Hyperoxaluria*. Kidney Int Rep. [PMID 39810772](https://pubmed.ncbi.nlm.nih.gov/39810772/)
74. Birtel J (2023). *The retinal phenotype in primary hyperoxaluria type 2 and 3*. Pediatr Nephrol. [PMID 36260161](https://pubmed.ncbi.nlm.nih.gov/36260161/)
75. Ge Y (2023). *HOGA1 variants in Chinese patients with primary hyperoxaluria type 3: genetic features and genotype-phenotype relationships*. World J Urol. [PMID 37318624](https://pubmed.ncbi.nlm.nih.gov/37318624/)
76. Baum MA (2023). *PHYOX2: a pivotal randomized study of nedosiran in primary hyperoxaluria type 1 or 2*. Kidney Int. [PMID 36007597](https://pubmed.ncbi.nlm.nih.gov/36007597/)
77. Garrelfs SF (2019). *Patients with primary hyperoxaluria type 2 have significant morbidity and require careful follow-up*. Kidney Int. [PMID 31685312](https://pubmed.ncbi.nlm.nih.gov/31685312/)
78. Konkoľová J (2017). *Severe child form of primary hyperoxaluria type 2 - a case report revealing consequence of GRHPR deficiency on metabolism*. BMC Med Genet. [PMID 28569194](https://pubmed.ncbi.nlm.nih.gov/28569194/)
79. Takayama T (2014). *Ethnic differences in GRHPR mutations in patients with primary hyperoxaluria type 2*. Clin Genet. [PMID 24116921](https://pubmed.ncbi.nlm.nih.gov/24116921/)
80. Levin-Iaina N (2009). *Late diagnosis of primary hyperoxaluria type 2 in the adult: effect of a novel mutation in GRHPR gene on enzymatic activity and molecular modeling*. J Urol. [PMID 19296982](https://pubmed.ncbi.nlm.nih.gov/19296982/)
81. Cregeen DP (2003). *Molecular analysis of the glyoxylate reductase (GRHPR) gene and description of mutations underlying primary hyperoxaluria type 2*. Hum Mutat. [PMID 14635115](https://pubmed.ncbi.nlm.nih.gov/14635115/)
82. Kemper MJ (1997). *Primary hyperoxaluria type 2*. Eur J Pediatr. [PMID 9243228](https://pubmed.ncbi.nlm.nih.gov/9243228/)
83. Mansell MA (1995). *Primary hyperoxaluria type 2*. Nephrol Dial Transplant. [PMID 8592629](https://pubmed.ncbi.nlm.nih.gov/8592629/)

## 5. 하이드록시프롤린·아스코르브산·에틸렌글리콜 — 기질 공급 (Substrate Supply)

84. Deng SK (2025). *Ethylene glycol poisoning: A case report and review of the literature*. World J Clin Cases. [PMID 40671749](https://pubmed.ncbi.nlm.nih.gov/40671749/)
85. Yang K (2024). *Isovaleramide attenuates ethylene glycol poisoning-induced acute kidney injury and reduces mortality by inhibiting alcohol dehydrogenase activity in rats*. Basic Clin Pharmacol Toxicol. [PMID 39324373](https://pubmed.ncbi.nlm.nih.gov/39324373/)
86. Ghannoum M (2023). *Extracorporeal treatment for ethylene glycol poisoning: systematic review and recommendations from the EXTRIP workgroup*. Crit Care. [PMID 36765419](https://pubmed.ncbi.nlm.nih.gov/36765419/)
87. Puiguriguer J (2022). *Calcium oxalate monohydrate crystalluria in ethylene glycol poisoning confirmed by scanning electron microscopy*. Clin Chim Acta. [PMID 35283093](https://pubmed.ncbi.nlm.nih.gov/35283093/)
88. Buchalski B (2020). *The effects of the inactivation of Hydroxyproline dehydrogenase on urinary oxalate and glycolate excretion in mouse models of primary hyperoxaluria*. Biochim Biophys Acta Mol Basis Dis. [PMID 31821850](https://pubmed.ncbi.nlm.nih.gov/31821850/)
89. Crivelli JJ (2020). *Contribution of Dietary Oxalate and Oxalate Precursors to Urinary Oxalate Excretion*. Nutrients. [PMID 33379176](https://pubmed.ncbi.nlm.nih.gov/33379176/)
90. Le Dudal M (2019). *Stiripentol protects against calcium oxalate nephrolithiasis and ethylene glycol poisoning*. J Clin Invest. [PMID 30946030](https://pubmed.ncbi.nlm.nih.gov/30946030/)
91. Fargue S (2018). *Hydroxyproline Metabolism and Oxalate Synthesis in Primary Hyperoxaluria*. J Am Soc Nephrol. [PMID 29588429](https://pubmed.ncbi.nlm.nih.gov/29588429/)
92. Dijcker JC (2014). *The effect of dietary hydroxyproline and dietary oxalate on urinary oxalate excretion in cats*. J Anim Sci. [PMID 24664562](https://pubmed.ncbi.nlm.nih.gov/24664562/)
93. Jiang J (2012). *Metabolism of [13C5]hydroxyproline in vitro and in vivo: implications for primary hyperoxaluria*. Am J Physiol Gastrointest Liver Physiol. [PMID 22207577](https://pubmed.ncbi.nlm.nih.gov/22207577/)
94. Salido E (2012). *Primary hyperoxalurias: disorders of glyoxylate detoxification*. Biochim Biophys Acta. [PMID 22446032](https://pubmed.ncbi.nlm.nih.gov/22446032/)
95. Khan SR (2007). *Dietary oxalate and calcium oxalate nephrolithiasis*. J Urol. [PMID 17870111](https://pubmed.ncbi.nlm.nih.gov/17870111/)
96. Lovrić M (2007). *Ethylene glycol poisoning*. Forensic Sci Int. [PMID 17629645](https://pubmed.ncbi.nlm.nih.gov/17629645/)
97. Knight J (2006). *Hydroxyproline ingestion and urinary oxalate and glycolate excretion*. Kidney Int. [PMID 17021603](https://pubmed.ncbi.nlm.nih.gov/17021603/)
98. Pearle MS (2001). *Prevention of nephrolithiasis*. Curr Opin Nephrol Hypertens. [PMID 11224695](https://pubmed.ncbi.nlm.nih.gov/11224695/)
99. Brent J (2001). *Current management of ethylene glycol poisoning*. Drugs. [PMID 11434452](https://pubmed.ncbi.nlm.nih.gov/11434452/)
100. Ogawa Y (2000). *Oxalate and urinary stones*. World J Surg. [PMID 11071450](https://pubmed.ncbi.nlm.nih.gov/11071450/)
101. Marangella M (1999). *Idiopathic calcium nephrolithiasis*. Nephron. [PMID 9873213](https://pubmed.ncbi.nlm.nih.gov/9873213/)
102. Gerster H (1997). *No contribution of ascorbic acid to renal calcium oxalate stones*. Ann Nutr Metab. [PMID 9429689](https://pubmed.ncbi.nlm.nih.gov/9429689/)
103. Fituri N (1983). *Urinary and plasma oxalate during ingestion of pure ascorbic acid: a re-evaluation*. Eur Urol. [PMID 6628476](https://pubmed.ncbi.nlm.nih.gov/6628476/)
104. Schmidt KH (1981). *Urinary oxalate excretion after large intakes of ascorbic acid in man*. Am J Clin Nutr. [PMID 7211731](https://pubmed.ncbi.nlm.nih.gov/7211731/)
105. Hatch M (1980). *Effect of megadoses of ascorbic acid on serum and urinary oxalate*. Eur Urol. [PMID 7371664](https://pubmed.ncbi.nlm.nih.gov/7371664/)
106. Parry MF (1974). *Ethylene glycol poisoning*. Am J Med. [PMID 4834513](https://pubmed.ncbi.nlm.nih.gov/4834513/)
107. Watts RW (1973). *Oxaluria*. J R Coll Physicians Lond. [PMID 4348042](https://pubmed.ncbi.nlm.nih.gov/4348042/)
108. Briggs MH (1973). *Urinary oxalate and vitamin-C supplements*. Lancet. [PMID 4124273](https://pubmed.ncbi.nlm.nih.gov/4124273/)

## 6. 옥살산 분포·혈장 옥살산·전신 옥살로시스 (Distribution, Plasma Oxalate, Systemic Oxalosis)

109. Chen X (2025). *Dysbiosis of the gut microbiota in calcium oxalate nephrolithiasis is associated with impaired short-chain fatty acid production and systemic metabolomic disruptions*. Microbiome. [PMID 41398611](https://pubmed.ncbi.nlm.nih.gov/41398611/)
110. Bargagli M (2025). *Kidney stone disease: risk factors, pathophysiology and management*. Nat Rev Nephrol. [PMID 40790363](https://pubmed.ncbi.nlm.nih.gov/40790363/)
111. Manolis AA (2025). *Nephrolithiasis and Cardiovascular Disease*. Cardiol Rev. [PMID 41398712](https://pubmed.ncbi.nlm.nih.gov/41398712/)
112. Sun H (2025). *Rising phytate and oxalate intake, declining calcium intake, and bone health in United States adults: 1999-2023, a serial cross-sectional analysis*. Am J Clin Nutr. [PMID 40409467](https://pubmed.ncbi.nlm.nih.gov/40409467/)
113. Puurunen M (2024). *Twenty-four-hour urine oxalate and risk of chronic kidney disease*. Nephrol Dial Transplant. [PMID 37804181](https://pubmed.ncbi.nlm.nih.gov/37804181/)
114. Ermer T (2023). *Oxalate homeostasis*. Nat Rev Nephrol. [PMID 36329260](https://pubmed.ncbi.nlm.nih.gov/36329260/)
115. Grocholski C (2023). *[Oxalate: from physiology to pathology]*. Nephrol Ther. [PMID 37166780](https://pubmed.ncbi.nlm.nih.gov/37166780/)
116. Nóbrega DF (2023). *Systemic oxalosis in a free-ranging green turtle (Chelonia mydas)*. J Comp Pathol. [PMID 36646034](https://pubmed.ncbi.nlm.nih.gov/36646034/)
117. Stamatelou K (2023). *Epidemiology of Kidney Stones*. Healthcare (Basel). [PMID 36766999](https://pubmed.ncbi.nlm.nih.gov/36766999/)
118. Stepanova N (2023). *Oxalate Homeostasis in Non-Stone-Forming Chronic Kidney Disease: A Review of Key Findings and Perspectives*. Biomedicines. [PMID 37371749](https://pubmed.ncbi.nlm.nih.gov/37371749/)
119. Oka Y (2022). *Calcium-Based Phosphate Binders and Plasma Oxalate Concentration in Dialysis Patients*. J Am Soc Nephrol. [PMID 35500940](https://pubmed.ncbi.nlm.nih.gov/35500940/)
120. Bargagli M (2021). *Calcium and Vitamin D Supplementation and Their Association with Kidney Stone Disease: A Narrative Review*. Nutrients. [PMID 34959915](https://pubmed.ncbi.nlm.nih.gov/34959915/)
121. Alshaikh AE (2021). *Gut-kidney axis in oxalate homeostasis*. Curr Opin Nephrol Hypertens. [PMID 33427760](https://pubmed.ncbi.nlm.nih.gov/33427760/)
122. Ben-Shalom E (2021). *Long-term complications of systemic oxalosis in children-a retrospective single-center cohort study*. Pediatr Nephrol. [PMID 33651179](https://pubmed.ncbi.nlm.nih.gov/33651179/)
123. El-Saygeh S (2021). *Calciphylaxis or vascular oxalosis?*. Clin Kidney J. [PMID 33564451](https://pubmed.ncbi.nlm.nih.gov/33564451/)
124. Siener R (2021). *Nutrition and Kidney Stone Disease*. Nutrients. [PMID 34204863](https://pubmed.ncbi.nlm.nih.gov/34204863/)
125. Uribarri J (2020). *Chronic kidney disease and kidney stones*. Curr Opin Nephrol Hypertens. [PMID 31972597](https://pubmed.ncbi.nlm.nih.gov/31972597/)
126. Ermer T (2016). *Oxalate, inflammasome, and progression of kidney disease*. Curr Opin Nephrol Hypertens. [PMID 27191349](https://pubmed.ncbi.nlm.nih.gov/27191349/)
127. Khan SR (2016). *Kidney stones*. Nat Rev Dis Primers. [PMID 27188687](https://pubmed.ncbi.nlm.nih.gov/27188687/)
128. Beck BB (2013). *Hyperoxaluria and systemic oxalosis: an update on current therapy and future directions*. Expert Opin Investig Drugs. [PMID 23167815](https://pubmed.ncbi.nlm.nih.gov/23167815/)
129. Maldonado I (2002). *Oxalate crystal deposition disease*. Curr Rheumatol Rep. [PMID 12010612](https://pubmed.ncbi.nlm.nih.gov/12010612/)
130. Adams ND (1992). *Nephrocalcinosis*. Clin Perinatol. [PMID 1576767](https://pubmed.ncbi.nlm.nih.gov/1576767/)
131. Olsen S (1989). *Primary acute renal failure ("acute tubular necrosis") in the transplanted kidney: morphology and pathogenesis*. Medicine (Baltimore). [PMID 2654537](https://pubmed.ncbi.nlm.nih.gov/2654537/)

## 7. 결정화 물리화학과 과포화도 (Crystallisation and Supersaturation)

132. De Mul A (2025). *Snacks and urinary oxalate: Which wins, almonds or chocolate?*. Fr J Urol. [PMID 41083049](https://pubmed.ncbi.nlm.nih.gov/41083049/)
133. Fargue S (2025). *Factors Influencing Oxalate Synthesis in Healthy Volunteers*. J Endourol. [PMID 40663545](https://pubmed.ncbi.nlm.nih.gov/40663545/)
134. Wang S (2025). *ACOT4 and ACOT6 Activate Akt-mTOR Pathway and Inhibit Calcium Oxalate-Induced Renal Tubular Cell Injury*. Kidney Blood Press Res. [PMID 40544826](https://pubmed.ncbi.nlm.nih.gov/40544826/)
135. Ferraro PM (2024). *24-Hour Urinary Chemistries and Kidney Stone Risk*. Am J Kidney Dis. [PMID 38583757](https://pubmed.ncbi.nlm.nih.gov/38583757/)
136. Siener R (2023). *Urinary Risk Profile, Impact of Diet, and Risk of Calcium Oxalate Urolithiasis in Idiopathic Uric Acid Stone Disease*. Nutrients. [PMID 36771279](https://pubmed.ncbi.nlm.nih.gov/36771279/)
137. Zainodini N (2023). *Associations of Oxalate Consumption and Some Individual Habits with the Risk of Kidney Stones*. Chin Med Sci J. [PMID 37643873](https://pubmed.ncbi.nlm.nih.gov/37643873/)
138. Stepanova N (2022). *Synbiotic supplementation and oxalate homeostasis in rats: focus on microbiota oxalate-degrading activity*. Urolithiasis. [PMID 35129638](https://pubmed.ncbi.nlm.nih.gov/35129638/)
139. Kavouras SA (2021). *Urine osmolality predicts calcium-oxalate crystallization risk in patients with recurrent urolithiasis*. Urolithiasis. [PMID 33635363](https://pubmed.ncbi.nlm.nih.gov/33635363/)
140. Queau Y (2019). *Nutritional Management of Urolithiasis*. Vet Clin North Am Small Anim Pract. [PMID 30583809](https://pubmed.ncbi.nlm.nih.gov/30583809/)
141. Marshall DJ (2018). *A combined liquid chromatography tandem mass spectrometry assay for the quantification of urinary oxalate and citrate in patients with nephrolithiasis*. Ann Clin Biochem. [PMID 28990817](https://pubmed.ncbi.nlm.nih.gov/28990817/)
142. Tiselius HG (2017). *Metabolic Work-up of Patients with Urolithiasis: Indications and Diagnostic Algorithm*. Eur Urol Focus. [PMID 28720369](https://pubmed.ncbi.nlm.nih.gov/28720369/)
143. Upala S (2016). *Risk of nephrolithiasis, hyperoxaluria, and calcium oxalate supersaturation increased after Roux-en-Y gastric bypass surgery: a systematic review and meta-analysis*. Surg Obes Relat Dis. [PMID 27396545](https://pubmed.ncbi.nlm.nih.gov/27396545/)
144. Massey LK (2007). *Food oxalate: factors affecting measurement, biological variation, and bioavailability*. J Am Diet Assoc. [PMID 17604750](https://pubmed.ncbi.nlm.nih.gov/17604750/)
145. Massey LK (2003). *Dietary influences on urinary oxalate and risk of kidney stones*. Front Biosci. [PMID 12700096](https://pubmed.ncbi.nlm.nih.gov/12700096/)
146. Holmes RP (1999). *Urinary oxalate and citrate*. Methods Mol Med. [PMID 21374302](https://pubmed.ncbi.nlm.nih.gov/21374302/)
147. Yamaguchi K (1997). *[Determination of urinary glycolate by ion chromatography: clinical and experimental implication]*. Nihon Hinyokika Gakkai Zasshi. [PMID 9465597](https://pubmed.ncbi.nlm.nih.gov/9465597/)
148. Tiselius HG (1991). *Aspects on estimation of the risk of calcium oxalate crystallization in urine*. Urol Int. [PMID 1781112](https://pubmed.ncbi.nlm.nih.gov/1781112/)
149. Ogawa Y (1984). *Determination of urinary oxalate by ion chromatography: some modifications*. Hinyokika Kiyo. [PMID 6377856](https://pubmed.ncbi.nlm.nih.gov/6377856/)
150. Menon M (1982). *Oxalate metabolism and renal calculi*. J Urol. [PMID 7035692](https://pubmed.ncbi.nlm.nih.gov/7035692/)
151. Brinkley L (1981). *Bioavailability of oxalate in foods*. Urology. [PMID 7245443](https://pubmed.ncbi.nlm.nih.gov/7245443/)

## 8. 결정 유발 신손상·염증·섬유화 (Crystal-Induced Injury, Inflammation, Fibrosis)

152. Sun Y (2025). *CREB1/CRTC2 regulated tubular epithelial-derived exosomal miR-93-3p promotes kidney injury induced by calcium oxalate via activating M1 polarization and macrophage extracellular trap formation*. J Nanobiotechnology. [PMID 40069788](https://pubmed.ncbi.nlm.nih.gov/40069788/)
153. Dong C (2025). *CHAC1 Mediates Endoplasmic Reticulum Stress-Dependent Ferroptosis in Calcium Oxalate Kidney Stone Formation*. Adv Sci (Weinh). [PMID 39836526](https://pubmed.ncbi.nlm.nih.gov/39836526/)
154. He Y (2025). *ROS Responsive Cerium Oxide Biomimetic Nanoparticles Alleviates Calcium Oxalate Crystals Induced Kidney Injury via Suppressing Oxidative Stress and M1 Macrophage Polarization*. Small. [PMID 39629501](https://pubmed.ncbi.nlm.nih.gov/39629501/)
155. Yuan T (2025). *PRMT1-mediated methylation of UBE2m promoting calcium oxalate crystal-induced kidney injury by inhibiting fatty acid metabolism*. Cell Death Dis. [PMID 40744915](https://pubmed.ncbi.nlm.nih.gov/40744915/)
156. Boldt AM (2025). *Targeting the NLRP3 inflammasome for calcium oxalate stones: pathophysiology and emerging pharmacological interventions*. Front Physiol. [PMID 40529990](https://pubmed.ncbi.nlm.nih.gov/40529990/)
157. Yang S (2025). *Regulating the balance between GSDMD-mediated pyroptosis and CHMP4B-dependent cell repair attenuates calcium oxalate kidney stone formation*. Int J Biol Sci. [PMID 40384863](https://pubmed.ncbi.nlm.nih.gov/40384863/)
158. Liu L (2025). *Gut microbiota-bile acid crosstalk contributes to calcium oxalate nephropathy through Hsp90α-mediated ferroptosis*. Cell Rep. [PMID 40591459](https://pubmed.ncbi.nlm.nih.gov/40591459/)
159. Ye Z (2025). *Lgals3 Promotes Calcium Oxalate Crystal Formation and Kidney Injury Through Histone Lactylation-Mediated FGFR4 Activation*. Adv Sci (Weinh). [PMID 39903812](https://pubmed.ncbi.nlm.nih.gov/39903812/)
160. Xia Y (2025). *EZH2-mediated macrophage-to-myofibroblast transition contributes to calcium oxalate crystal-induced kidney fibrosis*. Commun Biol. [PMID 39987296](https://pubmed.ncbi.nlm.nih.gov/39987296/)
161. Duan C (2024). *Sirtuin1 Suppresses Calcium Oxalate Nephropathy via Inhibition of Renal Proximal Tubular Cell Ferroptosis Through PGC-1α-mediated Transcriptional Coactivation*. Adv Sci (Weinh). [PMID 39498889](https://pubmed.ncbi.nlm.nih.gov/39498889/)
162. Yan X (2024). *The SOX4/EZH2/SLC7A11 signaling axis mediates ferroptosis in calcium oxalate crystal deposition-induced kidney injury*. J Transl Med. [PMID 38169402](https://pubmed.ncbi.nlm.nih.gov/38169402/)
163. Yan Q (2023). *NEAT1 Regulates Calcium Oxalate Crystal-Induced Renal Tubular Oxidative Injury via miR-130/IRF1*. Antioxid Redox Signal. [PMID 36242511](https://pubmed.ncbi.nlm.nih.gov/36242511/)
164. Duan C (2023). *Sirtuin1 inhibits calcium oxalate crystal-induced kidney injury by regulating TLR4 signaling and macrophage-mediated inflammatory activation*. Cell Signal. [PMID 37717713](https://pubmed.ncbi.nlm.nih.gov/37717713/)
165. Xu Z (2023). *Metabolic changes in kidney stone disease*. Front Immunol. [PMID 37228601](https://pubmed.ncbi.nlm.nih.gov/37228601/)
166. Syed YY (2023). *Nedosiran: First Approval*. Drugs. [PMID 38060091](https://pubmed.ncbi.nlm.nih.gov/38060091/)
167. Khan SR (2021). *Randall's plaque and calcium oxalate stone formation: role for immunity and inflammation*. Nat Rev Nephrol. [PMID 33514941](https://pubmed.ncbi.nlm.nih.gov/33514941/)
168. Wang Z (2021). *Recent advances on the mechanisms of kidney stone formation (Review)*. Int J Mol Med. [PMID 34132361](https://pubmed.ncbi.nlm.nih.gov/34132361/)
169. Milliner DS (2021). *Plasma oxalate and eGFR are correlated in primary hyperoxaluria patients with maintained kidney function-data from three placebo-controlled studies*. Pediatr Nephrol. [PMID 33515281](https://pubmed.ncbi.nlm.nih.gov/33515281/)
170. Daudon M (2018). *Drug-Induced Kidney Stones and Crystalline Nephropathy: Pathophysiology, Prevention and Treatment*. Drugs. [PMID 29264783](https://pubmed.ncbi.nlm.nih.gov/29264783/)
171. Shavit L (2015). *What is nephrocalcinosis?*. Kidney Int. [PMID 25807034](https://pubmed.ncbi.nlm.nih.gov/25807034/)

## 9. 장-신장 옥살산 축과 Oxalobacter (Gut-Kidney Axis)

172. Suryavanshi M (2026). *Predicting probiotic success: lessons from Oxalobacter and oxalate metabolism*. NPJ Biofilms Microbiomes. [PMID 41495071](https://pubmed.ncbi.nlm.nih.gov/41495071/)
173. Fargue S (2025). *Inducing Oxalobacter formigenes Colonization Reduces Urinary Oxalate in Healthy Adults*. Kidney Int Rep. [PMID 40485679](https://pubmed.ncbi.nlm.nih.gov/40485679/)
174. Suryavanshi M (2025). *Baseline abundance of oxalate-degrading bacteria determines response to Oxalobacter formigenes probiotic therapy*. Gut Microbes. [PMID 40984794](https://pubmed.ncbi.nlm.nih.gov/40984794/)
175. Joly PF (2024). *Pathophysiology and management of enteric hyperoxaluria*. Clin Res Hepatol Gastroenterol. [PMID 38734370](https://pubmed.ncbi.nlm.nih.gov/38734370/)
176. Desenclos J (2024). *Pathophysiology and management of enteric hyperoxaluria*. Clin Res Hepatol Gastroenterol. [PMID 38503362](https://pubmed.ncbi.nlm.nih.gov/38503362/)
177. Montoya A (2024). *Bidentate Substrate Binding Mode in Oxalate Decarboxylase*. Molecules. [PMID 39339409](https://pubmed.ncbi.nlm.nih.gov/39339409/)
178. Zan X (2024). *Recent Advances of Oxalate Decarboxylase: Biochemical Characteristics, Catalysis Mechanisms, and Gene Expression and Regulation*. J Agric Food Chem. [PMID 38653191](https://pubmed.ncbi.nlm.nih.gov/38653191/)
179. Arvans D (2023). *Sel1-like proteins and peptides are the major Oxalobacter formigenes-derived factors stimulating oxalate transport by human intestinal epithelial cells*. Am J Physiol Cell Physiol. [PMID 37125773](https://pubmed.ncbi.nlm.nih.gov/37125773/)
180. Wu F (2023). *Zn(2+) regulates human oxalate metabolism by manipulating oxalate decarboxylase to treat calcium oxalate stones*. Int J Biol Macromol. [PMID 36682657](https://pubmed.ncbi.nlm.nih.gov/36682657/)
181. Hiremath S (2022). *Oxalobacter formigenes: A new hope as a live biotherapeutic agent in the management of calcium oxalate renal stones*. Anaerobe. [PMID 35443224](https://pubmed.ncbi.nlm.nih.gov/35443224/)
182. Verhulst A (2022). *Oxalobacter formigenes treatment confers protective effects in a rat model of primary hyperoxaluria by preventing renal calcium oxalate deposition*. Urolithiasis. [PMID 35122487](https://pubmed.ncbi.nlm.nih.gov/35122487/)
183. Rosenstock JL (2022). *Oxalate nephropathy: a review*. Clin Kidney J. [PMID 35145635](https://pubmed.ncbi.nlm.nih.gov/35145635/)
184. Lemoine S (2022). *Reloxaliase in Enteric Hyperoxaluria - The Recent Brake*. NEJM Evid. [PMID 38319262](https://pubmed.ncbi.nlm.nih.gov/38319262/)
185. Alexander RT (2022). *Mechanisms Underlying Calcium Nephrolithiasis*. Annu Rev Physiol. [PMID 34699268](https://pubmed.ncbi.nlm.nih.gov/34699268/)
186. Cornière N (2022). *Dominant negative mutation in oxalate transporter SLC26A6 associated with enteric hyperoxaluria and nephrolithiasis*. J Med Genet. [PMID 35115415](https://pubmed.ncbi.nlm.nih.gov/35115415/)
187. Lieske JC (2022). *Randomized Placebo-Controlled Trial of Reloxaliase in Enteric Hyperoxaluria*. NEJM Evid. [PMID 38319254](https://pubmed.ncbi.nlm.nih.gov/38319254/)
188. Daniel SL (2021). *Forty Years of Oxalobacter formigenes, a Gutsy Oxalate-Degrading Specialist*. Appl Environ Microbiol. [PMID 34190610](https://pubmed.ncbi.nlm.nih.gov/34190610/)
189. Witting C (2021). *Pathophysiology and Treatment of Enteric Hyperoxaluria*. Clin J Am Soc Nephrol. [PMID 32900691](https://pubmed.ncbi.nlm.nih.gov/32900691/)
190. Whittamore JM (2021). *The anion exchanger PAT-1 (Slc26a6) does not participate in oxalate or chloride transport by mouse large intestine*. Pflugers Arch. [PMID 33205229](https://pubmed.ncbi.nlm.nih.gov/33205229/)
191. Chamberlain CA (2020). *Oxalobacter formigenes produces metabolites and lipids undetectable in oxalotrophic Bifidobacterium animalis*. Metabolomics. [PMID 33219444](https://pubmed.ncbi.nlm.nih.gov/33219444/)
192. Liu M (2019). *Enteric hyperoxaluria: role of microbiota and antibiotics*. Curr Opin Nephrol Hypertens. [PMID 31145706](https://pubmed.ncbi.nlm.nih.gov/31145706/)
193. Knauf F (2019). *Characterization of renal NaCl and oxalate transport in Slc26a6(-/-) mice*. Am J Physiol Renal Physiol. [PMID 30427220](https://pubmed.ncbi.nlm.nih.gov/30427220/)
194. Lingeman JE (2019). *ALLN-177, oral enzyme therapy for hyperoxaluria*. Int Urol Nephrol. [PMID 30783888](https://pubmed.ncbi.nlm.nih.gov/30783888/)
195. Weigert A (2018). *Novel therapeutic approaches in primary hyperoxaluria*. Expert Opin Emerg Drugs. [PMID 30540923](https://pubmed.ncbi.nlm.nih.gov/30540923/)
196. Arvans D (2017). *Oxalobacter formigenes-Derived Bioactive Factors Stimulate Oxalate Transport by Intestinal Epithelial Cells*. J Am Soc Nephrol. [PMID 27738124](https://pubmed.ncbi.nlm.nih.gov/27738124/)
197. Assimos DG (2017). *Re: A Double-Blind, Placebo Controlled, Randomized Phase 1 Cross-Over Study with ALLN-177, an Orally Administered Oxalate Degrading Enzyme*. J Urol. [PMID 28208535](https://pubmed.ncbi.nlm.nih.gov/28208535/)
198. Langman CB (2016). *A Double-Blind, Placebo Controlled, Randomized Phase 1 Cross-Over Study with ALLN-177, an Orally Administered Oxalate Degrading Enzyme*. Am J Nephrol. [PMID 27529510](https://pubmed.ncbi.nlm.nih.gov/27529510/)
199. Peck AB (2016). *Oxalate-degrading microorganisms or oxalate-degrading enzymes: which is the future therapy for enzymatic dissolution of calcium-oxalate uroliths in recurrent stone disease?*. Urolithiasis. [PMID 26645869](https://pubmed.ncbi.nlm.nih.gov/26645869/)
200. Iyalomhe O (2015). *The Structure and Function of OxlT, the Oxalate Transporter of Oxalobacter formigenes*. J Membr Biol. [PMID 25224873](https://pubmed.ncbi.nlm.nih.gov/25224873/)
201. Li X (2015). *Oxalobacter formigenes Colonization and Oxalate Dynamics in a Mouse Model*. Appl Environ Microbiol. [PMID 25979889](https://pubmed.ncbi.nlm.nih.gov/25979889/)
202. Knight J (2013). *The genetic composition of Oxalobacter formigenes and its relationship to colonization and calcium oxalate stone disease*. Urolithiasis. [PMID 23632911](https://pubmed.ncbi.nlm.nih.gov/23632911/)
203. Torzewska A (2013). *[Oxalobacter formigenes--characteristics and role in development of calcium oxalate urolithiasis]*. Postepy Hig Med Dosw (Online). [PMID 24379255](https://pubmed.ncbi.nlm.nih.gov/24379255/)
204. Aronson PS (2010). *Role of SLC26A6-mediated Cl⁻-oxalate exchange in renal physiology and pathophysiology*. J Nephrol. [PMID 21170874](https://pubmed.ncbi.nlm.nih.gov/21170874/)
205. Stewart CS (2004). *Oxalobacter formigenes and its role in oxalate metabolism in the human gut*. FEMS Microbiol Lett. [PMID 14734158](https://pubmed.ncbi.nlm.nih.gov/14734158/)
206. Verkoelen CF (1996). *Oxalate transport and calcium oxalate renal stone disease*. Urol Res. [PMID 8873376](https://pubmed.ncbi.nlm.nih.gov/8873376/)

## 10. 루마시란 (Lumasiran, HAO1 siRNA)

207. Frishberg Y (2026). *Final Results of the ILLUMINATE-A Phase 3 Clinical Trial of Lumasiran for Primary Hyperoxaluria 1*. Clin J Am Soc Nephrol. [PMID 41343248](https://pubmed.ncbi.nlm.nih.gov/41343248/)
208. Saland JM (2024). *Efficacy and Safety of Lumasiran in Patients With Primary Hyperoxaluria Type 1: Results from a Phase III Clinical Trial*. Kidney Int Rep. [PMID 39081738](https://pubmed.ncbi.nlm.nih.gov/39081738/)
209. Frishberg Y (2024). *Efficacy and safety of lumasiran for infants and young children with primary hyperoxaluria type 1: 30-month analysis of the phase 3 ILLUMINATE-B trial*. Front Pediatr. [PMID 39355649](https://pubmed.ncbi.nlm.nih.gov/39355649/)
210. Kang C (2024). *Lumasiran: A Review in Primary Hyperoxaluria Type 1*. Drugs. [PMID 38252335](https://pubmed.ncbi.nlm.nih.gov/38252335/)
211. Hayes W (2023). *Efficacy and safety of lumasiran for infants and young children with primary hyperoxaluria type 1: 12-month analysis of the phase 3 ILLUMINATE-B trial*. Pediatr Nephrol. [PMID 35913563](https://pubmed.ncbi.nlm.nih.gov/35913563/)
212. Gang X (2022). *Lumasiran for primary hyperoxaluria type 1: What we have learned?*. Front Pediatr. [PMID 36704142](https://pubmed.ncbi.nlm.nih.gov/36704142/)
213. Hulton SA (2022). *Randomized Clinical Trial on the Long-Term Efficacy and Safety of Lumasiran in Patients With Primary Hyperoxaluria Type 1*. Kidney Int Rep. [PMID 35257062](https://pubmed.ncbi.nlm.nih.gov/35257062/)
214. Sas DJ (2022). *Phase 3 trial of lumasiran for primary hyperoxaluria type 1: A new RNAi therapeutic in infants and young children*. Genet Med. [PMID 34906487](https://pubmed.ncbi.nlm.nih.gov/34906487/)

## 11. 네도시란 및 LDHA 표적 (Nedosiran, LDHA-directed therapy)

215. Lieske JC (2025). *PHYOX3: Nedosiran Long-Term Safety and Efficacy in Patients With Primary Hyperoxaluria Type 1*. Kidney Int Rep. [PMID 40630298](https://pubmed.ncbi.nlm.nih.gov/40630298/)
216. Zhang S (2025). *Population Pharmacokinetic and Pharmacodynamic Modelling and Simulation for Nedosiran Clinical Development and Dose Guidance in Pediatric Patients with Primary Hyperoxaluria Type 1*. Clin Pharmacokinet. [PMID 40601241](https://pubmed.ncbi.nlm.nih.gov/40601241/)
217. Traber GM (2024). *The Growing Class of Novel RNAi Therapeutics*. Mol Pharmacol. [PMID 38719476](https://pubmed.ncbi.nlm.nih.gov/38719476/)
218. Zhang S (2024). *Nedosiran population pharmacokinetic and pharmacodynamic modelling and simulation to guide clinical development and dose selection in patients with primary hyperoxaluria type 1*. Br J Clin Pharmacol. [PMID 39113219](https://pubmed.ncbi.nlm.nih.gov/39113219/)
219. Hoppe B (2022). *Safety, pharmacodynamics, and exposure-response modeling results from a first-in-human phase 1 study of nedosiran (PHYOX1) in primary hyperoxaluria*. Kidney Int. [PMID 34481803](https://pubmed.ncbi.nlm.nih.gov/34481803/)
220. Liu A (2022). *Nedosiran, a Candidate siRNA Drug for the Treatment of Primary Hyperoxaluria: Design, Development, and Clinical Studies*. ACS Pharmacol Transl Sci. [PMID 36407951](https://pubmed.ncbi.nlm.nih.gov/36407951/)
221. Moya-Garzon MD (2022). *New salicylic acid derivatives, double inhibitors of glycolate oxidase and lactate dehydrogenase, as effective agents decreasing oxalate production*. Eur J Med Chem. [PMID 35500475](https://pubmed.ncbi.nlm.nih.gov/35500475/)
222. Bacchetta J (2022). *Primary hyperoxaluria type 1: novel therapies at a glance*. Clin Kidney J. [PMID 35592618](https://pubmed.ncbi.nlm.nih.gov/35592618/)
223. Shee K (2021). *Nedosiran Dramatically Reduces Serum Oxalate in Dialysis-Dependent Primary Hyperoxaluria 1: A Compassionate Use Case Report*. Urology. [PMID 33774044](https://pubmed.ncbi.nlm.nih.gov/33774044/)
224. Letavernier E (2020). *Stiripentol identifies a therapeutic target to reduce oxaluria*. Curr Opin Nephrol Hypertens. [PMID 32452916](https://pubmed.ncbi.nlm.nih.gov/32452916/)
225. Wood KD (2019). *Reduction in urinary oxalate excretion in mouse models of Primary Hyperoxaluria by RNA interference inhibition of liver lactate dehydrogenase activity*. Biochim Biophys Acta Mol Basis Dis. [PMID 31055082](https://pubmed.ncbi.nlm.nih.gov/31055082/)
226. Stevens JS (2019). *Lactate dehydrogenase 5: identification of a druggable target to reduce oxaluria*. J Clin Invest. [PMID 31107247](https://pubmed.ncbi.nlm.nih.gov/31107247/)
227. Lai C (2018). *Specific Inhibition of Hepatic Lactate Dehydrogenase Reduces Oxalate Production in Mouse Models of Primary Hyperoxaluria*. Mol Ther. [PMID 29914758](https://pubmed.ncbi.nlm.nih.gov/29914758/)
228. Nickels KC (2017). *Stiripentol in the Management of Epilepsy*. CNS Drugs. [PMID 28434133](https://pubmed.ncbi.nlm.nih.gov/28434133/)
229. Sharma V (1992). *Oxalate production from glyoxylate by lactate dehydrogenase in vitro: inhibition by reduced glutathione, cysteine, cysteamine*. Biochem Int. [PMID 1417880](https://pubmed.ncbi.nlm.nih.gov/1417880/)

## 12. GalNAc-siRNA 전달과 PK/PD (GalNAc-siRNA Delivery and PK/PD)

230. An G (2024). *Pharmacokinetics and Pharmacodynamics of GalNAc-Conjugated siRNAs*. J Clin Pharmacol. [PMID 37589246](https://pubmed.ncbi.nlm.nih.gov/37589246/)
231. Tang Q (2024). *RNAi-based drug design: considerations and future directions*. Nat Rev Drug Discov. [PMID 38570694](https://pubmed.ncbi.nlm.nih.gov/38570694/)
232. Traber GM (2023). *RNAi-Based Therapeutics and Novel RNA Bioengineering Technologies*. J Pharmacol Exp Ther. [PMID 35680378](https://pubmed.ncbi.nlm.nih.gov/35680378/)
233. Brown KM (2022). *Expanding RNAi therapeutics to extrahepatic tissues with lipophilic conjugates*. Nat Biotechnol. [PMID 35654979](https://pubmed.ncbi.nlm.nih.gov/35654979/)
234. Belostotsky R (2022). *Catabolism of Hydroxyproline in Vertebrates: Physiology, Evolution, Genetic Diseases and New siRNA Approach for Treatment*. Int J Mol Sci. [PMID 35055190](https://pubmed.ncbi.nlm.nih.gov/35055190/)
235. Mackinnon SR (2022). *Novel Starting Points for Human Glycolate Oxidase Inhibitors, Revealed by Crystallography-Based Fragment Screening*. Front Chem. [PMID 35601556](https://pubmed.ncbi.nlm.nih.gov/35601556/)
236. Alshaer W (2021). *siRNA: Mechanism of action, challenges, and therapeutic approaches*. Eur J Pharmacol. [PMID 34044011](https://pubmed.ncbi.nlm.nih.gov/34044011/)
237. Brown CR (2020). *Investigating the pharmacodynamic durability of GalNAc-siRNA conjugates*. Nucleic Acids Res. [PMID 32808038](https://pubmed.ncbi.nlm.nih.gov/32808038/)
238. Debacker AJ (2020). *Delivery of Oligonucleotides to the Liver with GalNAc: From Research to Registered Therapeutic Drug*. Mol Ther. [PMID 32592692](https://pubmed.ncbi.nlm.nih.gov/32592692/)
239. Dindo M (2019). *Molecular basis of primary hyperoxaluria: clues to innovative treatments*. Urolithiasis. [PMID 30430197](https://pubmed.ncbi.nlm.nih.gov/30430197/)
240. Springer AD (2018). *GalNAc-siRNA Conjugates: Leading the Way for Delivery of RNAi Therapeutics*. Nucleic Acid Ther. [PMID 29792572](https://pubmed.ncbi.nlm.nih.gov/29792572/)
241. Foster DJ (2018). *Advanced siRNA Designs Further Improve In Vivo Performance of GalNAc-siRNA Conjugates*. Mol Ther. [PMID 29456020](https://pubmed.ncbi.nlm.nih.gov/29456020/)
242. Janas MM (2018). *The Nonclinical Safety Profile of GalNAc-conjugated RNAi Therapeutics in Subacute Studies*. Toxicol Pathol. [PMID 30139307](https://pubmed.ncbi.nlm.nih.gov/30139307/)
243. Nair JK (2017). *Impact of enhanced metabolic stability on pharmacokinetics and pharmacodynamics of GalNAc-siRNA conjugates*. Nucleic Acids Res. [PMID 28981809](https://pubmed.ncbi.nlm.nih.gov/28981809/)
244. Fitzgerald K (2017). *A Highly Durable RNAi Therapeutic Inhibitor of PCSK9*. N Engl J Med. [PMID 27959715](https://pubmed.ncbi.nlm.nih.gov/27959715/)
245. Parmar R (2016). *5'-(E)-Vinylphosphonate: A Stable Phosphate Mimic Can Improve the RNAi Activity of siRNA-GalNAc Conjugates*. Chembiochem. [PMID 27121751](https://pubmed.ncbi.nlm.nih.gov/27121751/)
246. Martin-Higueras C (2016). *Glycolate Oxidase Is a Safe and Efficient Target for Substrate Reduction Therapy in a Mouse Model of Primary Hyperoxaluria Type I*. Mol Ther. [PMID 26689264](https://pubmed.ncbi.nlm.nih.gov/26689264/)
247. Assimos DG (2016). *Re: Glycolate Oxidase is a Safe and Efficient Target for Substrate Reduction Therapy in a Mouse Model of Primary Hyperoxaluria Type I*. J Urol. [PMID 27321540](https://pubmed.ncbi.nlm.nih.gov/27321540/)
248. Nair JK (2014). *Multivalent N-acetylgalactosamine-conjugated siRNA localizes in hepatocytes and elicits robust RNAi-mediated gene silencing*. J Am Chem Soc. [PMID 25434769](https://pubmed.ncbi.nlm.nih.gov/25434769/)

## 13. 유전자치료·단백질 폴딩 구제 (Gene Therapy and Chaperone Rescue)

249. Ruta L (2025). *A Minor Haplotype Variant Determines the Pathogenicity of the p.Ile279Thr Substitution in the Primary Hyperoxaluria Type 1 Gene, AGXT*. J Inherit Metab Dis. [PMID 40495747](https://pubmed.ncbi.nlm.nih.gov/40495747/)
250. Zhang D (2024). *Lipid nanoparticle-mediated base-editing of the Hao1 gene achieves sustainable primary hyperoxaluria type 1 therapy in rats*. Sci China Life Sci. [PMID 39425833](https://pubmed.ncbi.nlm.nih.gov/39425833/)
251. Friedrich M (2022). *Therapeutic siRNA: State-of-the-Art and Future Perspectives*. BioDrugs. [PMID 35997897](https://pubmed.ncbi.nlm.nih.gov/35997897/)
252. Fernández-Higuero JÁ (2019). *Structural and functional insights on the roles of molecular chaperones in the mistargeting and aggregation phenotypes associated with primary hyperoxaluria type I*. Adv Protein Chem Struct Biol. [PMID 30635080](https://pubmed.ncbi.nlm.nih.gov/30635080/)
253. Oppici E (2018). *Folding Defects Leading to Primary Hyperoxaluria*. Handb Exp Pharmacol. [PMID 29071511](https://pubmed.ncbi.nlm.nih.gov/29071511/)
254. Pey AL (2013). *Protein homeostasis defects of alanine-glyoxylate aminotransferase: new therapeutic strategies in primary hyperoxaluria type I*. Biomed Res Int. [PMID 23956997](https://pubmed.ncbi.nlm.nih.gov/23956997/)
255. Danpure CJ (2003). *Alanine:glyoxylate aminotransferase peroxisome-to-mitochondrion mistargeting in human hereditary kidney stone disease*. Biochim Biophys Acta. [PMID 12686111](https://pubmed.ncbi.nlm.nih.gov/12686111/)
256. Danpure CJ (1993). *Primary hyperoxaluria type 1 and peroxisome-to-mitochondrion mistargeting of alanine:glyoxylate aminotransferase*. Biochimie. [PMID 8507692](https://pubmed.ncbi.nlm.nih.gov/8507692/)

## 14. 지지요법 — 수분·구연산 (Supportive Care: Fluids and Citrate)

257. Zomorodian A (2025). *Citrate and calcium kidney stones*. Clin Kidney J. [PMID 40978115](https://pubmed.ncbi.nlm.nih.gov/40978115/)
258. Balawender K (2024). *The Multidisciplinary Approach in the Management of Patients with Kidney Stone Disease-A State-of-the-Art Review*. Nutrients. [PMID 38931286](https://pubmed.ncbi.nlm.nih.gov/38931286/)
259. Peerapen P (2023). *Kidney Stone Prevention*. Adv Nutr. [PMID 36906146](https://pubmed.ncbi.nlm.nih.gov/36906146/)
260. Osther SS (2023). *Kidney stone disease*. Ugeskr Laeger. [PMID 37057692](https://pubmed.ncbi.nlm.nih.gov/37057692/)
261. Travers S (2023). *How to Monitor Hydration Status and Urine Dilution in Patients with Nephrolithiasis*. Nutrients. [PMID 37049482](https://pubmed.ncbi.nlm.nih.gov/37049482/)
262. Ferraro PM (2020). *Risk of Kidney Stones: Influence of Dietary Factors, Dietary Patterns, and Vegetarian-Vegan Diets*. Nutrients. [PMID 32183500](https://pubmed.ncbi.nlm.nih.gov/32183500/)
263. Malieckal DA (2020). *Occupational kidney stones*. Curr Opin Nephrol Hypertens. [PMID 31895162](https://pubmed.ncbi.nlm.nih.gov/31895162/)
264. Fontenelle LF (2019). *Kidney Stones: Treatment and Prevention*. Am Fam Physician. [PMID 30990297](https://pubmed.ncbi.nlm.nih.gov/30990297/)
265. Goldfarb DS (2019). *Empiric therapy for kidney stones*. Urolithiasis. [PMID 30478476](https://pubmed.ncbi.nlm.nih.gov/30478476/)
266. Moe OW (2006). *Kidney stones: pathophysiology and medical management*. Lancet. [PMID 16443041](https://pubmed.ncbi.nlm.nih.gov/16443041/)

## 15. 투석과 이식 (Dialysis and Transplantation)

267. Arena M (2025). *Simultaneous or sequential kidney-liver transplantation in primary hyperoxaluria*. J Nephrol. [PMID 39382784](https://pubmed.ncbi.nlm.nih.gov/39382784/)
268. Yi NJ (2024). *Combined liver-kidney transplantation in pediatric patients*. Pediatr Transplant. [PMID 38059323](https://pubmed.ncbi.nlm.nih.gov/38059323/)
269. Ranawaka R (2020). *Combined liver and kidney transplantation in children and long-term outcome*. World J Transplant. [PMID 33134116](https://pubmed.ncbi.nlm.nih.gov/33134116/)
270. Knotek M (2020). *Combined liver-kidney transplantation for rare diseases*. World J Hepatol. [PMID 33200012](https://pubmed.ncbi.nlm.nih.gov/33200012/)
271. Kotb MA (2019). *Combined liver-kidney transplantation for primary hyperoxaluria type I in children: Single Center Experience*. Pediatr Transplant. [PMID 30475440](https://pubmed.ncbi.nlm.nih.gov/30475440/)
272. Grenda R (2018). *Combined and sequential liver-kidney transplantation in children*. Pediatr Nephrol. [PMID 29322327](https://pubmed.ncbi.nlm.nih.gov/29322327/)
273. Dhondup T (2018). *Combined Liver-Kidney Transplantation for Primary Hyperoxaluria Type 2: A Case Report*. Am J Transplant. [PMID 28681512](https://pubmed.ncbi.nlm.nih.gov/28681512/)
274. Ermer T (2017). *Impact of Regular or Extended Hemodialysis and Hemodialfiltration on Plasma Oxalate Concentrations in Patients With End-Stage Renal Disease*. Kidney Int Rep. [PMID 29270514](https://pubmed.ncbi.nlm.nih.gov/29270514/)
275. Jalanko H (2014). *Combined liver and kidney transplantation in children*. Pediatr Nephrol. [PMID 23644898](https://pubmed.ncbi.nlm.nih.gov/23644898/)
276. Plumb TJ (2013). *Nocturnal home hemodialysis for a patient with type 1 hyperoxaluria*. Am J Kidney Dis. [PMID 23830800](https://pubmed.ncbi.nlm.nih.gov/23830800/)
277. Cochat P (2011). *[Primary hyperoxaluria]*. Nephrol Ther. [PMID 21636340](https://pubmed.ncbi.nlm.nih.gov/21636340/)
278. Illies F (2006). *Clearance and removal of oxalate in children on intensified dialysis for primary hyperoxaluria type 1*. Kidney Int. [PMID 16955107](https://pubmed.ncbi.nlm.nih.gov/16955107/)
279. Kemper MJ (2005). *Concurrent or sequential liver and kidney transplantation in children with primary hyperoxaluria type 1?*. Pediatr Transplant. [PMID 16269037](https://pubmed.ncbi.nlm.nih.gov/16269037/)
280. Cochat P (1999). *Combined liver-kidney transplantation in primary hyperoxaluria type 1*. Eur J Pediatr. [PMID 10603104](https://pubmed.ncbi.nlm.nih.gov/10603104/)
281. Kemper MJ (1998). *Preemptive liver transplantation in primary hyperoxaluria type 1: timing and preliminary results*. J Nephrol. [PMID 9604810](https://pubmed.ncbi.nlm.nih.gov/9604810/)
282. Costello JF (1992). *Effect of vitamin B6 supplementation on plasma oxalate and oxalate removal rate in hemodialysis patients*. J Am Soc Nephrol. [PMID 1450364](https://pubmed.ncbi.nlm.nih.gov/1450364/)
283. Marangella M (1992). *Plasma profiles and dialysis kinetics of oxalate in patients receiving hemodialysis*. Nephron. [PMID 1738418](https://pubmed.ncbi.nlm.nih.gov/1738418/)
284. Costello JF (1991). *Plasma oxalate levels rise in hemodialysis patients despite increased oxalate removal*. J Am Soc Nephrol. [PMID 1912391](https://pubmed.ncbi.nlm.nih.gov/1912391/)
285. Jacobsen D (1986). *Methanol and ethylene glycol poisonings. Mechanism of toxicity, clinical course, diagnosis and treatment*. Med Toxicol. [PMID 3537623](https://pubmed.ncbi.nlm.nih.gov/3537623/)
286. Ramsay AG (1984). *Oxalate removal by hemodialysis in end-stage renal disease*. Am J Kidney Dis. [PMID 6475942](https://pubmed.ncbi.nlm.nih.gov/6475942/)
287. Watts RW (1984). *Oxalate dynamics and removal rates during haemodialysis and peritoneal dialysis in patients with primary hyperoxaluria and severe renal failure*. Clin Sci (Lond). [PMID 6368103](https://pubmed.ncbi.nlm.nih.gov/6368103/)

## 16. QSP 방법론 (QSP Methodology)

288. Ye Z (2025). *Luteolin alleviated calcium oxalate crystal induced kidney injury by inhibiting Nr4a1-mediated ferroptosis*. Phytomedicine. [PMID 39662099](https://pubmed.ncbi.nlm.nih.gov/39662099/)
289. Li S (2025). *Gut microbiota-regulated unconjugated bilirubin metabolism drives renal calcium oxalate crystal deposition*. Gut Microbes. [PMID 40849919](https://pubmed.ncbi.nlm.nih.gov/40849919/)
290. Yuan T (2023). *STAT6 promoting oxalate crystal deposition-induced renal fibrosis by mediating macrophage-to-myofibroblast transition via inhibiting fatty acid oxidation*. Inflamm Res. [PMID 37924395](https://pubmed.ncbi.nlm.nih.gov/37924395/)
291. Maddah E (2022). *A quantitative systems pharmacology model of plasma potassium regulation by the kidney and aldosterone*. J Pharmacokinet Pharmacodyn. [PMID 35776281](https://pubmed.ncbi.nlm.nih.gov/35776281/)
292. Yasui T (2017). *Pathophysiology-based treatment of urolithiasis*. Int J Urol. [PMID 27539983](https://pubmed.ncbi.nlm.nih.gov/27539983/)
293. Strazzullo P (1994). *Hypertension, calcium metabolism, and nephrolithiasis*. Am J Med Sci. [PMID 8141146](https://pubmed.ncbi.nlm.nih.gov/8141146/)
294. Kleeman CR (1980). *Kidney stones*. West J Med. [PMID 7385835](https://pubmed.ncbi.nlm.nih.gov/7385835/)

---

## 문헌 수 (Reference count)

총 **294편**. 모든 PMID는 NCBI E-utilities로 실제 조회하여 확인했습니다.

## 데이터베이스·도구 (Databases and tools)

- PubMed / NCBI E-utilities — <https://www.ncbi.nlm.nih.gov/books/NBK25501/>
- OMIM: PH1 259900 (AGXT) · PH2 260000 (GRHPR) · PH3 613616 (HOGA1) —
  <https://www.omim.org/>
- OxalEurope / Rare Kidney Stone Consortium 레지스트리 —
  <https://www.oxaleurope.com/>, <https://www.rarekidneystones.org/>
- mrgsolve — <https://mrgsolve.org/>
- Graphviz — <https://graphviz.org/>

