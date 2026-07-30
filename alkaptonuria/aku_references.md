# 알캅톤뇨증 (Alkaptonuria, AKU) — 참고문헌
### References for the AKU QSP model

모든 인용은 PubMed E-utilities(`esummary`)로 저자·연도·저널·제목을 **기계적으로 검증**한
것입니다. 기억에 의존해 쓴 PMID는 하나도 남기지 않았습니다. 본문에서 어떤 인용이
모델의 어떤 항(파라미터·구조·검증 목표)에 쓰였는지 밝혔고, **모델이 문헌과 불일치하는
지점도 그 인용 아래에 그대로 적었습니다.**

총 **131편**. 앵커(보정에 쓰인 값)와 held-out(검증에만 쓴 값)의 구분은
[README.md](README.md)의 검증 표와 대응합니다.

---

## 1. 정의·종합 리뷰 (Definition and comprehensive reviews)

1. Bernardini G et al. (2024). *Alkaptonuria.* Nat Rev Dis Primers. [PMID 38453957](https://pubmed.ncbi.nlm.nih.gov/38453957/)
   > 현 시점의 결정판 리뷰. 본 모델의 병태생리 골격(HGD 결손 → HGA 축적 → BQA → 콜라겐 중합 → 오크로노시스)과 장기별 침범 분포의 1차 근거.

2. Davison AS et al. (2023). *Alkaptonuria - Past, present and future.* Adv Clin Chem. [PMID 37268334](https://pubmed.ncbi.nlm.nih.gov/37268334/)
   > AKU 연구의 역사·현재·미래를 임상화학 관점에서 정리. 측정 가능한 바이오마커 목록의 근거.

3. Zatkova A et al. (2020). *Alkaptonuria: Current Perspectives.* Appl Clin Genet. [PMID 32158253](https://pubmed.ncbi.nlm.nih.gov/32158253/)

4. Ranganath LR et al. (2013). *Recent advances in management of alkaptonuria (invited review; best practice article).* J Clin Pathol. [PMID 23486607](https://pubmed.ncbi.nlm.nih.gov/23486607/)

5. Gao Y (2026). *Alkaptonuria.* J Clin Rheumatol. [PMID 41451484](https://pubmed.ncbi.nlm.nih.gov/41451484/)

6. Efridi W et al. (2026). *Ochronosis.* StatPearls. [PMID 32809369](https://pubmed.ncbi.nlm.nih.gov/32809369/)


## 2. HGD 유전학과 집단유전 (HGD genetics and population genetics)

7. Zatkova A (2011). *An update on molecular genetics of Alkaptonuria (AKU).* J Inherit Metab Dis. [PMID 21720873](https://pubmed.ncbi.nlm.nih.gov/21720873/)
   > HGD 변이 스펙트럼. 모델의 RESACT(잔존 효소활성 0–5%) 파라미터가 미스센스/널 구분에 대응하는 근거.

8. Zatková A et al. (2000). *High frequency of alkaptonuria in Slovakia: evidence for the appearance of multiple mutations in HGO involving different mutational hot spots.* Am J Hum Genet. [PMID 11017803](https://pubmed.ncbi.nlm.nih.gov/11017803/)
   > 슬로바키아 고빈도(약 1:19,000)와 다중 변이의 출현. PREV/FOUNDER 노드의 근거.

9. Srsen S et al. (2002). *Alkaptonuria in Slovakia: thirty-two years of research on phenotype and genotype.* Mol Genet Metab. [PMID 12051967](https://pubmed.ncbi.nlm.nih.gov/12051967/)

10. Srsen S et al. (1978). *Screening for alkaptonuria in the newborn in Slovakia.* Lancet. [PMID 79941](https://pubmed.ncbi.nlm.nih.gov/79941/)

11. Srsen S et al. (1978). *Alkaptonuria in the Trencín District of Czechoslovakia.* Am J Med Genet. [PMID 263435](https://pubmed.ncbi.nlm.nih.gov/263435/)

12. Zatkova A et al. (2012). *Identification of 11 Novel Homogentisate 1,2 Dioxygenase Variants in Alkaptonuria Patients and Establishment of a Novel LOVD-Based HGD Mutation Database.* JIMD Rep. [PMID 23430897](https://pubmed.ncbi.nlm.nih.gov/23430897/)

13. Danda S et al. (2020). *Founder effects of the homogentisate 1,2-dioxygenase (HGD) gene in a gypsy population and mutation spectrum in the gene among alkaptonuria patients from India.* Clin Rheumatol. [PMID 32212000](https://pubmed.ncbi.nlm.nih.gov/32212000/)

14. Soltysova A et al. (2022). *Alkaptonuria in Russia.* Eur J Hum Genet. [PMID 34504318](https://pubmed.ncbi.nlm.nih.gov/34504318/)

15. Müller CR et al. (1999). *Allelic heterogeneity of alkaptonuria in Central Europe.* Eur J Hum Genet. [PMID 10482952](https://pubmed.ncbi.nlm.nih.gov/10482952/)

16. Tao L et al. (2022). *A novel mutation in the homogentisate 1,2 dioxygenase gene identified in Chinese Hani pediatric patients with Alkaptonuria.* Clin Chim Acta. [PMID 35550814](https://pubmed.ncbi.nlm.nih.gov/35550814/)

17. Abraham SSC et al. (2024). *Gene expression & biochemical analysis in alkaptonuria caused by a founder pathogenic variant across different age groups from India.* Indian J Med Res. [PMID 39737503](https://pubmed.ncbi.nlm.nih.gov/39737503/)

18. Cicaloni V et al. (2019). *Interactive alkaptonuria database: investigating clinical data to improve patient care in a rare disease.* FASEB J. [PMID 31462106](https://pubmed.ncbi.nlm.nih.gov/31462106/)


## 3. 자연사·중증도 척도 (Natural history and severity scoring)

19. Phornphutkul C et al. (2002). *Natural history of alkaptonuria.* N Engl J Med. [PMID 12501223](https://pubmed.ncbi.nlm.nih.gov/12501223/)
   > 본 모델의 임상 시점 앵커 전부의 출처: 관절치환 평균 55세, 심장 판막 침범 54세, 관상동맥 칼슘 59세, 신결석 64세, 방사선학적 점수는 30세 이후
   > 상승(남성에서 더 빠름). 모델은 이 다섯 시점을 파라미터로 넣지 않고 검증용으로만 사용.

20. Cox TF et al. (2011). *A quantitative assessment of alkaptonuria: testing the reliability of two disease severity scoring systems.* J Inherit Metab Dis. [PMID 21744089](https://pubmed.ncbi.nlm.nih.gov/21744089/)
   > AKUSSI 원 개발 논문. 임상/관절/척추 3개 도메인 구조가 모델의 AKJ·AKS·AKC 도메인 분해에 대응.

21. Langford B et al. (2018). *Alkaptonuria Severity Score Index Revisited: Analysing the AKUSSI and Its Subcomponent Features.* JIMD Rep. [PMID 29654544](https://pubmed.ncbi.nlm.nih.gov/29654544/)
   > AKUSSI 하위 항목의 주성분 분석 — 어떤 항목이 정보를 담고 어떤 항목이 중복인지.

22. Cant HEO et al. (2022). *Improving the clinical accuracy and flexibility of the Alkaptonuria severity score index.* JIMD Rep. [PMID 35822087](https://pubmed.ncbi.nlm.nih.gov/35822087/)
   > cAKUSSI 2.0. 편향된 저정보 측정치를 제거한 개정판.

23. Khedr M et al. (2023). *First decade anniversary of the United Kingdom National Alkaptonuria Centre.* JIMD Rep. [PMID 36873083](https://pubmed.ncbi.nlm.nih.gov/36873083/)
   > 영국 국립 AKU 센터(NAC) 첫 10년. 2 mg 저용량 실사용 코호트의 근거.

24. Azami A et al. (2014). *Alkaptonuric ochronosis: a clinical study from Ardabil, Iran.* Int J Rheum Dis. [PMID 24447956](https://pubmed.ncbi.nlm.nih.gov/24447956/)

25. Akbaba AI et al. (2020). *Presentation of 14 alkaptonuria patients from Turkey.* J Pediatr Endocrinol Metab. [PMID 31927521](https://pubmed.ncbi.nlm.nih.gov/31927521/)

26. Srsen S (1979). *Alkaptonuria.* Johns Hopkins Med J. [PMID 513428](https://pubmed.ncbi.nlm.nih.gov/513428/)

27. Ranganath LR et al. (2024). *Joint replacement risk is markedly increased in alkaptonuria (AKU) in those with prior arthroplasty.* Mol Genet Metab Rep. [PMID 38846518](https://pubmed.ncbi.nlm.nih.gov/38846518/)
   > 이전 관절치환 이력이 있는 환자에서 추가 치환 위험이 현저히 증가 — 모델의 HZJR 가속(연골 소실이 부하 재분배를 통해 반대측을 가속)과 일치.


## 4. 티로신 경로 플럭스와 질량보존 (Pathway flux and mass balance) — 본 모델의 BALANCE 1

28. Milan AM et al. (2019). *Quantification of the flux of tyrosine pathway metabolites during nitisinone treatment of Alkaptonuria.* Sci Rep. [PMID 31296884](https://pubmed.ncbi.nlm.nih.gov/31296884/)
   > 니티시논을 플럭스 프로브로 사용해 티로신 경로 대사물의 흐름을 정량화. 모델의 '플럭스는 보존되고 출구만 바뀐다'는 구조의 직접 근거.

29. Ranganath LR et al. (2022). *Revisiting Quantification of Phenylalanine/Tyrosine Flux in the Ochronotic Pathway during Long-Term Nitisinone Treatment of Alkaptonuria.* Metabolites. [PMID 36295821](https://pubmed.ncbi.nlm.nih.gov/36295821/)
   > 위 논문의 재검토. 무치료 시 HGA 등가물 7.1 g/일(총 Phe/Tyr 플럭스의 46.4%), 치료 시 총 15.32 g/일이 드러나며 8.21
   > g/일(53.6%)이 미계상. 저자들은 담도·장 배출을 제안. 본 모델은 이 미계상분이 비가역적 배출로 보기에는 과대(70년 누적 200 kg 이상)하며, 확장된
   > 티로신 풀과 가역적 아미노기전이가 더 그럴듯하다는 대안 해석을 정량적으로 제시함(README의 '모델이 중재하는 두 해석' 참조).

30. Norman BP et al. (2022). *Comprehensive Biotransformation Analysis of Phenylalanine-Tyrosine Metabolism Reveals Alternative Routes of Metabolite Clearance in Nitisinone-Treated Alkaptonuria.* Metabolites. [PMID 36295829](https://pubmed.ncbi.nlm.nih.gov/36295829/)
   > 치료 중 우회 배출 경로의 정량: 요중 HPPA 20–23배, HPLA 9.4–10.6배, N-아세틸-티로신 5.7–6.8배, 티로신-글루쿠로나이드 8배, 티로신
   > 4–4.7배. 모델의 티로신 탈출로(J_TYRU·J_CONJ·J_HPPU·J_HPLAU) 용량 배분의 근거.

31. Ranganath LR et al. (2022). *Determinants of tyrosinaemia during nitisinone therapy in alkaptonuria.* Sci Rep. [PMID 36167967](https://pubmed.ncbi.nlm.nih.gov/36167967/)
   > SONIA 2의 307개 채혈점 분석. 티로시네미아의 정도를 좌우하는 것은 용량이 아니라 HPPA→HPLA 전환의 상대적 감소, 즉 탈출로 용량이라는 결론. 모델의
   > '용량은 HGA를, 식이와 탈출로가 티로신을 정한다'는 비대칭과 독립적으로 일치.

32. Norman BP et al. (2022). *Metabolomic studies in the inborn error of metabolism alkaptonuria reveal new biotransformations in tyrosine metabolism.* Genes Dis. [PMID 35685462](https://pubmed.ncbi.nlm.nih.gov/35685462/)

33. Norman BP et al. (2019). *A Comprehensive LC-QTOF-MS Metabolic Phenotyping Strategy: Application to Alkaptonuria.* Clin Chem. [PMID 30782595](https://pubmed.ncbi.nlm.nih.gov/30782595/)

34. Davison AS et al. (2019). *Evaluation of the serum metabolome of patients with alkaptonuria before and after two years of treatment with nitisinone using LC-QTOF-MS.* JIMD Rep. [PMID 31392115](https://pubmed.ncbi.nlm.nih.gov/31392115/)

35. Grasso D et al. (2022). *Untargeted NMR Metabolomics Reveals Alternative Biomarkers and Pathways in Alkaptonuria.* Int J Mol Sci. [PMID 36555443](https://pubmed.ncbi.nlm.nih.gov/36555443/)

36. Serafimov K et al. (2025). *Targeted and untargeted urinary metabolomics of alkaptonuria patients using ultra high-performance liquid chromatography-tandem mass spectrometry.* J Pharm Biomed Anal. [PMID 39842076](https://pubmed.ncbi.nlm.nih.gov/39842076/)

37. Gertsman I et al. (2015). *Perturbations of tyrosine metabolism promote the indolepyruvate pathway via tryptophan in host and microbiome.* Mol Genet Metab. [PMID 25680927](https://pubmed.ncbi.nlm.nih.gov/25680927/)

38. Milan AM et al. (2017). *The effect of nitisinone on homogentisic acid and tyrosine: a two-year survey of patients attending the National Alkaptonuria Centre, Liverpool.* Ann Clin Biochem. [PMID 28081634](https://pubmed.ncbi.nlm.nih.gov/28081634/)
   > NAC 2년 추적에서 HGA와 티로신의 실사용 값 — 2 mg 코호트 앵커의 출처.


## 5. HGA의 신장 처리 (Renal handling of HGA) — 본 모델의 BALANCE 3

39. Ranganath LR et al. (2020). *Homogentisic acid is not only eliminated by glomerular filtration and tubular secretion but also produced in the kidney in alkaptonuria.* J Inherit Metab Dis. [PMID 31609457](https://pubmed.ncbi.nlm.nih.gov/31609457/)
   > 225명 분석. HGA 신클리어런스와 분획배설이 이론적 최대 신혈류량을 초과 → 세뇨관 분비가 사구체여과보다 크고 신장 자체가 HGA를 생산. 또한 순환 HGA는
   > 연령과 함께 증가하며 이는 HGA 클리어런스 감소와 유의하게 연관. 모델의 CL_HGA 상한·FREN(신장내 생산)·RFUN(연령 의존 감소) 세 항의 유일한
   > 근거이자, '무치료 진행이 노년에 가속되는 이유가 색소가 색소를 부르기 때문만은 아니다'라는 주장의 근거.

40. Wolff F et al. (2015). *Renal and prostate stones composition in alkaptonuria: a case report.* Clin Nephrol. [PMID 26396096](https://pubmed.ncbi.nlm.nih.gov/26396096/)
   > AKU 신결석·전립선 결석의 조성 분석.

41. Huledal G et al. (2019). *Non randomized study on the potential of nitisinone to inhibit cytochrome P450 2C9, 2D6, 2E1 and the organic anion transporters OAT1 and OAT3 in healthy volunteers.* Eur J Clin Pharmacol. [PMID 30443705](https://pubmed.ncbi.nlm.nih.gov/30443705/)
   > 니티시논이 CYP2C9/2D6/2E1과 유기음이온수송체(OAT)를 억제할 잠재력에 대한 연구. 모델이 요중/혈청 HGA 괴리를 설명하기 위해 도입한 OAT 경쟁 항의
   > 약리학적 타당성 근거.


## 6. 산화·중합·색소 침착과 전임상 모델 (Oxidation, polymerisation, pigment; preclinical)

42. Preston AJ et al. (2014). *Ochronotic osteoarthropathy in a mouse model of alkaptonuria, and its inhibition by nitisinone.* Ann Rheum Dis. [PMID 23511227](https://pubmed.ncbi.nlm.nih.gov/23511227/)
   > AKU 마우스 모델에서 오크로노시스성 골관절병증과 니티시논에 의한 억제. 색소 침착 속도만 낮출 수 있고 이미 침착된 색소는 되돌릴 수 없다는 모델 구조의 in
   > vivo 근거.

43. Hughes JH et al. (2021). *Anatomical Distribution of Ochronotic Pigment in Alkaptonuric Mice is Associated with Calcified Cartilage Chondrocytes at Osteochondral Interfaces.* Calcif Tissue Int. [PMID 33057760](https://pubmed.ncbi.nlm.nih.gov/33057760/)
   > 마우스에서 색소의 해부학적 분포가 석회화 연골의 연골세포와 연관. 모델이 색소를 연골·디스크·판막에 우선 배분하고(FCART·FDISC·FVALV) 표피성 조직과
   > 구분하는 근거.

44. Ranganath LR et al. (2020). *Reversal of ochronotic pigmentation in alkaptonuria following nitisinone therapy: Analysis of data from the United Kingdom National Alkaptonuria Centre.* JIMD Rep. [PMID 32904992](https://pubmed.ncbi.nlm.nih.gov/32904992/)
   > 니티시논 치료 후 오크로노시스성 색소의 '역전'. 본 모델은 이를 반증이 아니라 예측으로 다룸: 콜라겐이 실제로 교체되는 조직(귀·공막·피부)에서만 소실항이 존재할
   > 수 있고 관절 연골·디스크에는 존재할 수 없으므로, 관찰된 역전이 표피성 조직에 국한된다는 것이 모델의 검증 가능한 귀결(모델에서 PSCL·PEAR·PSKIN에만
   > 1차 소실항을 둠).

45. Geminiani M et al. (2017). *Cytoskeleton Aberrations in Alkaptonuric Chondrocytes.* J Cell Physiol. [PMID 27454006](https://pubmed.ncbi.nlm.nih.gov/27454006/)

46. Tinti L et al. (2011). *A novel ex vivo organotypic culture model of alkaptonuria-ochronosis.* Clin Exp Rheumatol. [PMID 21813063](https://pubmed.ncbi.nlm.nih.gov/21813063/)

47. Tinti L et al. (2010). *Evaluation of antioxidant drugs for the treatment of ochronotic alkaptonuria in an in vitro human cell model.* J Cell Physiol. [PMID 20648626](https://pubmed.ncbi.nlm.nih.gov/20648626/)

48. Braconi D et al. (2010). *Proteomic and redox-proteomic evaluation of homogentisic acid and ascorbic acid effects on human articular chondrocytes.* J Cell Biochem. [PMID 20665660](https://pubmed.ncbi.nlm.nih.gov/20665660/)
   > HGA와 아스코르브산이 사람 관절연골세포에 미치는 프로테옴·레독스 효과. 아스코르브산이 BQA를 되돌리지만 임상적으로 실패하는 이유를 모델이 구조로 설명하는
   > 근거(중합은 병렬 비가역 단계).

49. Mastroeni P et al. (2024). *An in vitro cell model for exploring inflammatory and amyloidogenic events in alkaptonuria.* J Cell Physiol. [PMID 39351877](https://pubmed.ncbi.nlm.nih.gov/39351877/)

50. Mastroeni P et al. (2026). *HGA-Induced Oxidative Stress Impairs Autophagy via Lysosomal Dysfunction in Alkaptonuria.* Antioxidants (Basel). [PMID 42510554](https://pubmed.ncbi.nlm.nih.gov/42510554/)

51. Davison AS et al. (2022). *Impact of Nitisinone on the Cerebrospinal Fluid Metabolome of a Murine Model of Alkaptonuria.* Metabolites. [PMID 35736410](https://pubmed.ncbi.nlm.nih.gov/35736410/)

52. Davison AS et al. (2019). *Assessing the effect of nitisinone induced hypertyrosinaemia on monoamine neurotransmitters in brain tissue from a murine model of alkaptonuria using mass spectrometry imaging.* Metabolomics. [PMID 31037385](https://pubmed.ncbi.nlm.nih.gov/31037385/)
   > 마우스 뇌조직에서 니티시논 유발 고티로신혈증과 모노아민 신경전달물질. 모델의 NEUROCOG 노드에 대응.

53. Gallagher JA et al. (2015). *Lessons from rare diseases of cartilage and bone.* Curr Opin Pharmacol. [PMID 25978274](https://pubmed.ncbi.nlm.nih.gov/25978274/)


## 7. 아밀로이드·염증·산화 바이오마커 (Amyloid, inflammation, oxidative-stress markers)

54. Millucci L et al. (2012). *Alkaptonuria is a novel human secondary amyloidogenic disease.* Biochim Biophys Acta. [PMID 22850426](https://pubmed.ncbi.nlm.nih.gov/22850426/)
   > AKU가 이차성 아밀로이드 생성 질환이라는 제안.

55. Millucci L et al. (2015). *Amyloidosis in alkaptonuria.* J Inherit Metab Dis. [PMID 25868666](https://pubmed.ncbi.nlm.nih.gov/25868666/)

56. Millucci L et al. (2014). *Amyloidosis, inflammation, and oxidative stress in the heart of an alkaptonuric patient.* Mediators Inflamm. [PMID 24876668](https://pubmed.ncbi.nlm.nih.gov/24876668/)

57. Braconi D et al. (2017). *Homogentisic acid induces aggregation and fibrillation of amyloidogenic proteins.* Biochim Biophys Acta Gen Subj. [PMID 27865997](https://pubmed.ncbi.nlm.nih.gov/27865997/)

58. Mastroeni P et al. (2024). *HGA Triggers SAA Aggregation and Accelerates Fibril Formation in the C20/A4 Alkaptonuria Cell Model.* Cells. [PMID 39273071](https://pubmed.ncbi.nlm.nih.gov/39273071/)

59. Braconi D et al. (2022). *Effects of Nitisinone on Oxidative and Inflammatory Markers in Alkaptonuria: Results from SONIA1 and SONIA2 Studies.* Cells. [PMID 36429096](https://pubmed.ncbi.nlm.nih.gov/36429096/)
   > SONIA 1과 SONIA 2에서 니티시논의 산화·염증 지표에 대한 효과. 모델의 활막염 증폭 루프가 전신 염증지표를 거의 올리지 않는(CRP 거의 정상) 구조와
   > 일치.

60. Braconi D et al. (2018). *Inflammatory and oxidative stress biomarkers in alkaptonuria: data from the DevelopAKUre project.* Osteoarthritis Cartilage. [PMID 29852277](https://pubmed.ncbi.nlm.nih.gov/29852277/)

61. Trezza A et al. (2025). *Integrated Clinomics and Molecular Dynamics Simulation Approaches Reveal the SAA1.1 Allele as a Biomarker in Alkaptonuria Disease Severity.* Biomolecules. [PMID 40001497](https://pubmed.ncbi.nlm.nih.gov/40001497/)


## 8. 관절·척추 병변과 정형외과적 치료 (Arthropathy, spine, surgery)

62. Ranganath LR et al. (2021). *Characterising the arthroplasty in spondyloarthropathy in a large cohort of eighty-seven patients with alkaptonuria.* J Inherit Metab Dis. [PMID 33314212](https://pubmed.ncbi.nlm.nih.gov/33314212/)
   > AKU 87명의 관절치환 특성 분석 — 모델의 관절치환 위험도(HZJR) 시점 검증의 주요 근거.

63. Salem KH et al. (2025). *Ochronotic arthropathy: skeletal manifestations and orthopaedic treatment.* EFORT Open Rev. [PMID 40071956](https://pubmed.ncbi.nlm.nih.gov/40071956/)

64. Borman P et al. (2002). *Ochronotic arthropathy.* Rheumatol Int. [PMID 11958438](https://pubmed.ncbi.nlm.nih.gov/11958438/)

65. Wu K et al. (2019). *Musculoskeletal manifestations of alkaptonuria: A case report and literature review.* Eur J Rheumatol. [PMID 30451653](https://pubmed.ncbi.nlm.nih.gov/30451653/)

66. Al-Mahfoudh R et al. (2008). *Alkaptonuria presenting with ochronotic spondyloarthropathy.* Br J Neurosurg. [PMID 19085367](https://pubmed.ncbi.nlm.nih.gov/19085367/)

67. Sang P et al. (2025). *Alkaptonuria presenting as lumbar degenerative disease: A case report and literature review.* Medicine (Baltimore). [PMID 39833049](https://pubmed.ncbi.nlm.nih.gov/39833049/)

68. Pesciallo C et al. (2022). *Total Knee Replacement in Alkaptonuric Ochronosis.* Acta Biomed. [PMID 35671127](https://pubmed.ncbi.nlm.nih.gov/35671127/)

69. Patel VG (2015). *Total knee arthroplasty in ochronosis.* Arthroplast Today. [PMID 28326376](https://pubmed.ncbi.nlm.nih.gov/28326376/)

70. Basanagoudar PL et al. (2025). *Total Joint Arthroplasty in Ochronotic Arthritis of Lower Extremities.* Indian J Orthop. [PMID 40852538](https://pubmed.ncbi.nlm.nih.gov/40852538/)

71. Tanios M et al. (2023). *Spondyloarthropathies That Mimic Ankylosing Spondylitis: A Narrative Review.* Clin Med Insights Arthritis Musculoskelet Disord. [PMID 37533960](https://pubmed.ncbi.nlm.nih.gov/37533960/)

72. Fisher AA et al. (2004). *Alkaptonuric ochronosis with aortic valve and joint replacements and femoral fracture: a case report and literature review.* Clin Med Res. [PMID 15931360](https://pubmed.ncbi.nlm.nih.gov/15931360/)


## 9. 심혈관 침범 (Cardiovascular involvement)

73. Bruce C et al. (2025). *Effect of Nitisinone on Aortic Stenosis Disease Progression in Patients With Alkaptonuria: An Analysis of the Suitability of Nitisinone in Alkaptonuria (SONIA) 2 Study.* Cureus. [PMID 41589138](https://pubmed.ncbi.nlm.nih.gov/41589138/)
   > SONIA 2의 대동맥판 협착 진행 분석. 기저 대동맥판 협착 13.0%(경증 44.4%·중등증 33.3%·중증 22.2%), 대동맥판 경화 18.1%, 4년간 최대
   > 압력차 진행률의 군간 차이 0.0093 mmHg/년(p=0.53, 비유의). 모델의 판막 축(PVALV→VALVCA→PMAXS)이 재현해야 할 held-out 목표.

74. Ather N et al. (2020). *Cardiovascular ochronosis.* Cardiovasc Pathol. [PMID 32473412](https://pubmed.ncbi.nlm.nih.gov/32473412/)

75. Thimmapuram R et al. (2020). *Aortic distensibility in alkaptonuria.* Mol Genet Metab. [PMID 32466960](https://pubmed.ncbi.nlm.nih.gov/32466960/)


## 10. 안과·이과·피부 소견 (Ocular, aural and cutaneous signs)

76. Lindner M et al. (2014). *On the ocular findings in ochronosis: a systematic review of literature.* BMC Ophthalmol. [PMID 24479547](https://pubmed.ncbi.nlm.nih.gov/24479547/)
   > AKU 안과 소견의 체계적 문헌고찰 — 공막 색소(Osler 징후)의 시기와 분포.

77. Carlson DM et al. (1991). *Ocular ochronosis from alkaptonuria.* J Am Optom Assoc. [PMID 1813514](https://pubmed.ncbi.nlm.nih.gov/1813514/)

78. Okutucu M et al. (2019). *Glaucoma With Alkaptonuria as a Result of Pigment Accumulation.* J Glaucoma. [PMID 31274704](https://pubmed.ncbi.nlm.nih.gov/31274704/)

79. de Azevedo Magalhaes O et al. (2022). *Descemet's membrane folds in ochronosis: a case report.* J Med Case Rep. [PMID 36183119](https://pubmed.ncbi.nlm.nih.gov/36183119/)

80. Pau HW (1984). *[Involvement of the tympanic membrane and ear ossicle system in ochronotic alkaptonuria].* Laryngol Rhinol Otol (Stuttg). [PMID 6503572](https://pubmed.ncbi.nlm.nih.gov/6503572/)

81. Yin ES et al. (2017). *A 54-year-old woman with arthritis and discoloration of the hands, ears, and sclerae.* Int J Dermatol. [PMID 27651033](https://pubmed.ncbi.nlm.nih.gov/27651033/)


## 11. 골대사·힘줄 (Bone and tendon)

82. Ranganath LR et al. (2021). *Frequency, diagnosis, pathogenesis and management of osteoporosis in alkaptonuria: data analysis from the UK National Alkaptonuria Centre.* Osteoporos Int. [PMID 33118050](https://pubmed.ncbi.nlm.nih.gov/33118050/)
   > NAC 데이터에서 AKU 골다공증의 빈도·기전·관리. 모델의 BMD 항(만성 활막염 매개 골소실)의 근거.

83. Abate M et al. (2016). *Tendons Involvement in Congenital Metabolic Disorders.* Adv Exp Med Biol. [PMID 27535253](https://pubmed.ncbi.nlm.nih.gov/27535253/)

84. Yokoe T et al. (2024). *Direct repair of the chronic ochronotic Achilles tendon rupture: a case report.* BMC Musculoskelet Disord. [PMID 39448996](https://pubmed.ncbi.nlm.nih.gov/39448996/)


## 12. HPD 효소학과 억제제 분자약리 (HPD enzymology and inhibitor pharmacology)

85. Ellis MK et al. (1995). *Inhibition of 4-hydroxyphenylpyruvate dioxygenase by 2-(2-nitro-4-trifluoromethylbenzoyl)-cyclohexane-1,3-dione and 2-(2-chloro-4-methanesulfonylbenzoyl)-cyclohexane-1,3-dione.* Toxicol Appl Pharmacol. [PMID 7597701](https://pubmed.ncbi.nlm.nih.gov/7597701/)
   > NTBC의 HPD 억제 동역학 원 논문. 모델의 KI_NT(약 16 nmol/L로 적합)와 보고된 IC50 수준의 비교 근거.

86. Kavana M et al. (2003). *Interaction of (4-hydroxyphenyl)pyruvate dioxygenase with the specific inhibitor 2-[2-nitro-4-(trifluoromethyl)benzoyl]-1,3-cyclohexanedione.* Biochemistry. [PMID 12939152](https://pubmed.ncbi.nlm.nih.gov/12939152/)
   > HPD와 NTBC의 상호작용 — 느리고 단단한 결합. 모델이 단순 경쟁적 억제보다 급한 용량-반응(Hill 지수 > 1)을 쓰는 근거.

87. Lin HY et al. (2019). *Molecular insights into the mechanism of 4-hydroxyphenylpyruvate dioxygenase inhibition: enzyme kinetics, X-ray crystallography and computational simulations.* FEBS J. [PMID 30632699](https://pubmed.ncbi.nlm.nih.gov/30632699/)

88. Santucci A et al. (2017). *4-Hydroxyphenylpyruvate Dioxygenase and Its Inhibition in Plants and Animals: Small Molecules as Herbicides and Agents for the Treatment of Human Inherited Diseases.* J Med Chem. [PMID 28128559](https://pubmed.ncbi.nlm.nih.gov/28128559/)

89. Liu YX et al. (2020). *Identification of key residues determining the binding specificity of human 4-hydroxyphenylpyruvate dioxygenase.* Eur J Pharm Sci. [PMID 32750420](https://pubmed.ncbi.nlm.nih.gov/32750420/)

90. Yang DY (2003). *4-Hydroxyphenylpyruvate dioxygenase as a drug discovery target.* Drug News Perspect. [PMID 14668946](https://pubmed.ncbi.nlm.nih.gov/14668946/)

91. Brownlee JM et al. (2010). *Product analysis and inhibition studies of a causative Asn to Ser variant of 4-hydroxyphenylpyruvate dioxygenase suggest a simple route to the treatment of Hawkinsinuria.* Biochemistry. [PMID 20677779](https://pubmed.ncbi.nlm.nih.gov/20677779/)


## 13. 니티시논 PK/PD (Nitisinone pharmacokinetics and pharmacodynamics)

92. Olsson B et al. (2015). *Relationship Between Serum Concentrations of Nitisinone and Its Effect on Homogentisic Acid and Tyrosine in Patients with Alkaptonuria.* JIMD Rep. [PMID 25772318](https://pubmed.ncbi.nlm.nih.gov/25772318/)
   > 혈청 니티시논 농도와 HGA·티로신 효과의 관계. 모델의 PK(2 mg에서 약 1 umol/L, 10 mg에서 약 5 umol/L)와 노출–반응 연결의 직접 근거.

93. Guffon N et al. (2018). *Open-Label Single-Sequence Crossover Study Evaluating Pharmacokinetics, Efficacy, and Safety of Once-Daily Dosing of Nitisinone in Patients with Hereditary Tyrosinemia Type 1.* JIMD Rep. [PMID 28643275](https://pubmed.ncbi.nlm.nih.gov/28643275/)
   > 1일 1회 투여의 PK·효능·안전성 교차연구. 반감기 약 54시간이 투여간격을 무의미하게 만든다는 모델 예측(S03 vs 분할투여 대조)의 근거.

94. Hall MG et al. (2001). *Pharmacokinetics and pharmacodynamics of NTBC (2-(2-nitro-4-fluoromethylbenzoyl)-1,3-cyclohexanedione) and mesotrione, inhibitors of 4-hydroxyphenyl pyruvate dioxygenase (HPPD) following a single dose to healthy male volunteers.* Br J Clin Pharmacol. [PMID 11488774](https://pubmed.ncbi.nlm.nih.gov/11488774/)

95. Anonymous (2002). *Nitisinone. Ntbc, Orfadin.* Drugs R D. [PMID 12001819](https://pubmed.ncbi.nlm.nih.gov/12001819/)

96. Gertsman I et al. (2015). *Metabolic Effects of Increasing Doses of Nitisinone in the Treatment of Alkaptonuria.* JIMD Rep. [PMID 25665838](https://pubmed.ncbi.nlm.nih.gov/25665838/)
   > 니티시논 용량을 올릴 때의 대사 효과 — 용량-반응의 포화 지점을 보여주는 자료.

97. Kienstra NS et al. (2018). *Daily variation of NTBC and its relation to succinylacetone in tyrosinemia type 1 patients comparing a single dose to two doses a day.* J Inherit Metab Dis. [PMID 29170874](https://pubmed.ncbi.nlm.nih.gov/29170874/)

98. Schlune A et al. (2012). *Single dose NTBC-treatment of hereditary tyrosinemia type I.* J Inherit Metab Dis. [PMID 22307209](https://pubmed.ncbi.nlm.nih.gov/22307209/)


## 14. 임상시험과 규제 (Clinical trials and regulatory)

99. Ranganath LR et al. (2016). *Suitability Of Nitisinone In Alkaptonuria 1 (SONIA 1): an international, multicentre, randomised, open-label, no-treatment controlled, parallel-group, dose-response study to investigate the effect of once daily nitisinone on 24-h urinary homogentisic acid excretion in patients with alkaptonuria after 4 weeks of treatment.* Ann Rheum Dis. [PMID 25475116](https://pubmed.ncbi.nlm.nih.gov/25475116/)
   > SONIA 1(무작위, 5군, 4주). 4주 시점 보정 기하평균 24시간 요중 HGA가 무치료 31.53 mmol, 1 mg 3.26, 2 mg 1.44, 4 mg
   > 0.57, 8 mg 0.15 mmol. 8 mg에서 기저 대비 98.8% 감소. 본 모델 보정의 1차 앵커 5개 전부의 출처.

100. Ranganath LR et al. (2020). *Efficacy and safety of once-daily nitisinone for patients with alkaptonuria (SONIA 2): an international, multicentre, open-label, randomised controlled trial.* Lancet Diabetes Endocrinol. [PMID 32822600](https://pubmed.ncbi.nlm.nih.gov/32822600/)
   > SONIA 2(무작위, 4년, 10 mg/일, n=138). 12개월 u-HGA24가 대조군 대비 99.7% 감소(기하평균비 0.003), 48개월 cAKUSSI
   > 증가폭이 대조군보다 8.6점 적음(-16.0 ~ -1.2, p=0.023). 모델의 held-out 임상 목표.

101. Introne WJ et al. (2011). *A 3-year randomized therapeutic trial of nitisinone in alkaptonuria.* Mol Genet Metab. [PMID 21620748](https://pubmed.ncbi.nlm.nih.gov/21620748/)
   > NIH 3년 무작위 시험(2 mg/일). HGA는 크게 감소했으나 1차 평가변수(고관절 회전)는 유의하게 개선되지 않음. 모델은 이 '음성 결과'를 재현해야 하며,
   > 그 이유(46세에 시작하면 예방 가능한 적분이 이미 대부분 소진)를 정량적으로 제시함.

102. Ranganath LR et al. (2022). *Comparing nitisinone 2 mg and 10 mg in the treatment of alkaptonuria-An approach using statistical modelling.* JIMD Rep. [PMID 35028273](https://pubmed.ncbi.nlm.nih.gov/35028273/)
   > 2 mg과 10 mg의 통계적 모델링 비교. 혈청 니티시논 1.26 vs 4.34 umol/L, 요중 HGA 약 94% vs 99.5% 감소, 혈청 HGA
   > 26.2→3.86 vs 30.3→2.23 umol/L, 혈청 티로신 53.9→782 vs 65.3→875 umol/L(용량 간 상승률 차이 유의하지 않음), 각막병증
   > 3/60(5%) vs 10/69(14.5%), AKUSSI 기울기 0.19 vs 0.06 점/월. 본 모델의 중심 주장(용량은 HGA를, 식이가 티로신을 정한다)을
   > 지지하는 결정적 자료이자, 2 mg의 혈청/요중 HGA 쌍이 질량보존과 양립하지 않는다는 모델의 지적 대상.

103. Sloboda N et al. (2019). *Efficacy of low dose nitisinone in the management of alkaptonuria.* Mol Genet Metab. [PMID 31235217](https://pubmed.ncbi.nlm.nih.gov/31235217/)

104. Abbas K et al. (2022). *Adequacy of nitisinone for the management of alkaptonuria.* Ann Med Surg (Lond). [PMID 36045846](https://pubmed.ncbi.nlm.nih.gov/36045846/)

105. Ranganath LR et al. (2023). *Clinical development innovation in rare diseases: overcoming barriers to successful delivery of a randomised clinical trial in alkaptonuria-a mini-review.* Orphanet J Rare Dis. [PMID 36600285](https://pubmed.ncbi.nlm.nih.gov/36600285/)

106. Chandani HK et al. (2025). *Harliku (nitisinone): first FDA-approved disease-modifying therapy for alkaptonuria.* Ann Med Surg (Lond). [PMID 41377225](https://pubmed.ncbi.nlm.nih.gov/41377225/)
   > 니티시논의 AKU 적응증 FDA 승인(Harliku). EMA 승인에 이어 규제적으로 확립된 최초의 질병조절 치료.

107. Ooi N et al. (2023). *Evaluation of Homogentisic Acid, a Prospective Antibacterial Agent Highlighted by the Suitability of Nitisinone in Alkaptonuria 2 (SONIA 2) Clinical Trial.* Cells. [PMID 37443717](https://pubmed.ncbi.nlm.nih.gov/37443717/)


## 15. 티로시네미아·각막병증·식이 관리 (Tyrosinaemia, keratopathy and dietary management)

108. Hughes JH et al. (2020). *Dietary restriction of tyrosine and phenylalanine lowers tyrosinemia associated with nitisinone therapy of alkaptonuria.* J Inherit Metab Dis. [PMID 31503358](https://pubmed.ncbi.nlm.nih.gov/31503358/)
   > 단백 제한만으로도, 그리고 Tyr/Phe-무함유 아미노산 보충을 병용하면 티로신이 유의하게 감소(10명 중 4명이 700 umol/L 미만 달성). 모델의 '식이가
   > 티로신을 정한다'는 예측의 직접 검증 자료이자 S15/S17 시나리오의 근거.

109. Olsson B et al. (2022). *Effects of a protein-restricted diet on body weight and serum tyrosine concentrations in patients with alkaptonuria.* JIMD Rep. [PMID 35028270](https://pubmed.ncbi.nlm.nih.gov/35028270/)
   > 단백 제한식이의 체중과 혈청 티로신에 대한 효과.

110. Ranganath LR et al. (2024). *Anthropometric, Body Composition, and Nutritional Indicators with and without Nutritional Intervention during Nitisinone Therapy in Alkaptonuria.* Nutrients. [PMID 39203858](https://pubmed.ncbi.nlm.nih.gov/39203858/)
   > 니티시논 치료 중 영양 개입 유무에 따른 인체계측·체성분·영양지표.

111. Teke Kisa P et al. (2022). *Efficacy of Phenylalanine- and Tyrosine-Restricted Diet in Alkaptonuria Patients on Nitisinone Treatment: Case Series and Review of Literature.* Ann Nutr Metab. [PMID 34736252](https://pubmed.ncbi.nlm.nih.gov/34736252/)
   > Phe/Tyr 제한식이 사례군과 문헌고찰. 티로신 501–700에서 0.9 g/kg, 701–900에서 0.8 g/kg, 900 초과에서 Phe/Tyr-무함유 교환식
   > 추가라는 실제 프로토콜의 출처 — 모델 S14–S17의 식이 수준 설정 근거.

112. Ranganath LR et al. (2022). *Comparing the Phenylalanine/Tyrosine Pathway and Related Factors between Keratopathy and No-Keratopathy Groups as Well as between Genders in Alkaptonuria during Nitisinone Treatment.* Metabolites. [PMID 36005644](https://pubmed.ncbi.nlm.nih.gov/36005644/)
   > 각막병증군과 비각막병증군의 Phe/Tyr 경로 비교. 모델의 CORTYR(각막 결정 부하) 임계 구조가 겨냥하는 자료.

113. Ranganath L et al. (2023). *Increased prevalence of Parkinson's disease in alkaptonuria.* JIMD Rep. [PMID 37404676](https://pubmed.ncbi.nlm.nih.gov/37404676/)
   > AKU에서 파킨슨병 유병률 증가. 티로신-도파민 축을 통한 잠재적 연결(모델의 NEUROCOG 노드, 불확실로 표시).

114. Rudebeck M et al. (2020). *A patient survey on the impact of alkaptonuria symptoms as perceived by the patients and their experiences of receiving diagnosis and care.* JIMD Rep. [PMID 32395411](https://pubmed.ncbi.nlm.nih.gov/32395411/)

115. Davison AS et al. (2016). *Acute fatal metabolic complications in alkaptonuria.* J Inherit Metab Dis. [PMID 26596578](https://pubmed.ncbi.nlm.nih.gov/26596578/)
   > AKU의 급성 치명적 대사 합병증 — 드물지만 모델이 다루지 않는 경로임을 명시하기 위한 인용.


## 16. 유전성 티로시네미아 1형에서의 니티시논 장기 경험 (외삽 근거)

116. Spiekerkoetter U et al. (2021). *Long-term safety and outcomes in hereditary tyrosinaemia type 1 with nitisinone treatment: a 15-year non-interventional, multicentre study.* Lancet Diabetes Endocrinol. [PMID 34023005](https://pubmed.ncbi.nlm.nih.gov/34023005/)
   > 15년 비개입 장기 안전성 연구. 소아기부터의 평생 노출에서 관찰된 안전성 프로파일 — 모델의 S09(5세 개시) 시나리오가 기대는 유일한 인체 근거이자, 그 외삽이
   > AKU에서 시험되지 않았음을 명시해야 하는 이유.

117. Geppert J et al. (2017). *Evaluation of pre-symptomatic nitisinone treatment on long-term outcomes in Tyrosinemia type 1 patients: a systematic review.* Orphanet J Rare Dis. [PMID 28893311](https://pubmed.ncbi.nlm.nih.gov/28893311/)
   > 증상 발현 전 니티시논 치료의 장기 결과 — 조기 개시의 가치에 대한 HT-1에서의 증거.

118. Chinsky JM et al. (2017). *Diagnosis and treatment of tyrosinemia type I: a US and Canadian consensus group review and recommendations.* Genet Med. [PMID 28771246](https://pubmed.ncbi.nlm.nih.gov/28771246/)

119. van Ginkel WG et al. (2019). *Long-Term Outcomes and Practical Considerations in the Pharmacological Management of Tyrosinemia Type 1.* Paediatr Drugs. [PMID 31667718](https://pubmed.ncbi.nlm.nih.gov/31667718/)

120. Kuypers AM et al. (2025). *Overview of European Practices for Management of Tyrosinemia Type 1: Towards European Guidelines.* J Inherit Metab Dis. [PMID 40965374](https://pubmed.ncbi.nlm.nih.gov/40965374/)

121. Holme E et al. (1995). *Diagnosis and management of tyrosinemia type I.* Curr Opin Pediatr. [PMID 8776026](https://pubmed.ncbi.nlm.nih.gov/8776026/)

122. Aktuglu Zeybek AC et al. (2022). *Evaluation of dynamic thiol/disulfide homeostasis in hereditary tyrosinemia type 1 patients.* Pediatr Res. [PMID 34628487](https://pubmed.ncbi.nlm.nih.gov/34628487/)


## 17. 아스코르브산과 기타 시도된 치료 (Ascorbate and other attempted therapies)

123. Wolff JA et al. (1989). *Effects of ascorbic acid in alkaptonuria: alterations in benzoquinone acetic acid and an ontogenic effect in infancy.* Pediatr Res. [PMID 2771520](https://pubmed.ncbi.nlm.nih.gov/2771520/)
   > 아스코르브산이 BQA를 변화시키지만 임상 경과를 바꾸지 못함, 영아기의 연령 효과. 모델이 아스코르브산 실패를 구조로(비가역 중합이 병렬 경로) 재현하는 근거.

124. Mayatepek E et al. (1998). *Effects of ascorbic acid and low-protein diet in alkaptonuria.* Eur J Pediatr. [PMID 9809834](https://pubmed.ncbi.nlm.nih.gov/9809834/)

125. Morava E et al. (2003). *Reversal of clinical symptoms and radiographic abnormalities with protein restriction and ascorbic acid in alkaptonuria.* Ann Clin Biochem. [PMID 12542920](https://pubmed.ncbi.nlm.nih.gov/12542920/)

126. Forslind K et al. (1988). *Alkaptonuria and ochronosis in three siblings. Ascorbic acid treatment monitored by urinary HGA excretion.* Clin Exp Rheumatol. [PMID 3180550](https://pubmed.ncbi.nlm.nih.gov/3180550/)

127. Kamoun P et al. (1992). *Ascorbic acid and alkaptonuria.* Eur J Pediatr. [PMID 1537362](https://pubmed.ncbi.nlm.nih.gov/1537362/)


## 18. 개발 중 치료 전략 (Investigational strategies)

128. Lequeue S et al. (2025). *A robust bacterial high-throughput screening assay to identify pharmacological chaperones targeting human homogentisate 1,2-dioxygenase missense variants in alkaptonuria.* Eur J Pharmacol. [PMID 40784658](https://pubmed.ncbi.nlm.nih.gov/40784658/)
   > HGD를 표적하는 약리학적 샤페론 발굴을 위한 세균 기반 고속 스크리닝. 모델의 CHAPERONE→RESACT 경로에 대응하며, 미스센스 유전형에서만 의미가 있다는
   > 구조적 제약을 공유.

129. Paulk NK et al. (2012). *In vivo selection of transplanted hepatocytes by pharmacological inhibition of fumarylacetoacetate hydrolase in wild-type mice.* Mol Ther. [PMID 22871666](https://pubmed.ncbi.nlm.nih.gov/22871666/)

130. Vernon HJ et al. (2021). *Milestones in treatments for inborn errors of metabolism: Reflections on Where chemistry and medicine meet.* Am J Med Genet A. [PMID 34165242](https://pubmed.ncbi.nlm.nih.gov/34165242/)


## 19. QSP 방법론 (QSP methodology)

131. Elmokadem A et al. (2019). *Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.* CPT Pharmacometrics Syst Pharmacol. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
   > mrgsolve로 QSP·PBPK 모델을 구축하는 실습형 튜토리얼 — 본 모델 구현의 도구 근거.

---

## 인용을 쓰지 않은 곳 (Where no citation is claimed)

아래 값들은 문헌에서 직접 인용한 것이 **아니고** 모델 내부에서 유도하거나 임의로
설정한 것입니다. 근거가 있는 것처럼 보이지 않도록 여기에 모아 둡니다.

- **색소 분배 비율** (`FCART` 0.30 · `FDISC` 0.22 · `FVALV` 0.08 · `FTEND` 0.12 ·
  `FSCL` 0.03 · `FEAR` 0.05 · `FSKIN` 0.05 · `FOTH` 0.15): 부검·영상 기술 문헌의
  정성적 침범 빈도를 순서만 보존해 배분한 것으로, 정량적 근거는 없습니다.
- **콜라겐 결합 용량** (`COLLC` 3000 · `COLLD` 2200 · `COLLV` 800 · `COLLT` 1500 umol):
  색소 밀도를 0–1 범위로 규격화하기 위한 스케일링 상수입니다. 절대값에 의미는 없고
  `PD50`과 함께 하나의 비율로만 작동합니다.
- **취성화 Hill 지수** `HPD_H` = 6: 3번째 10년대에 증상이 급하게 나타나는 관찰을
  재현할 수 있는 최소한의 급경사로 선택했습니다. 조직학적 근거는 없습니다.
- **취성화 임계값** `PD50` 0.60 / `PD50D` 0.45 / `PD50V` 0.70 / `PD50T` 0.80:
  디스크 석회화가 가장 먼저, 힘줄 파열이 가장 늦게 나타나는 임상 순서만을 근거로
  한 순서 제약이며, 개별 값은 임의입니다.
- **AKUSSI 도메인 가중치** (45 / 40 / 30 / 15): 실제 AKUSSI의 57개 하위 항목 구조를
  4개 도메인으로 축약한 것으로, 원 척도의 항목별 배점을 재현하지 않습니다. 따라서
  모델의 cAKUSSI 절대값은 임상 점수와 직접 비교할 수 없고, **변화율(점/월)만**
  비교 가능합니다. 검증 표가 절대값이 아니라 기울기를 쓰는 이유입니다.
- **통증 모듈 계수** (`WNOCI_*`, `KCS`, `KCSOFF`): AKU에 특화된 통증 정신물리학
  자료가 없어, 만성 골관절염 문헌의 일반적 형태만 차용했습니다.
- **니티시논 말초 구획** (`V2`, `Q`): 2구획 구조를 유지하기 위한 값으로, AKU 환자의
  집단 PK 분석에서 추정된 것이 아닙니다. 정상 상태 농도는 `CLNT` 하나로 결정됩니다.
- **`SRC_INS` = 7.5 umol/day**: 니티시논에 반응하지 않는 HGA 공급원을 가정하고
  적합시킨 결과 **거의 0으로 수렴**했습니다. 즉 자료는 '억제되지 않는 HGA 공급원'
  가설보다 '용량–반응이 단순 경쟁적 억제보다 급하다'(`HNT` 1.45)는 가설을 지지합니다.
  가설을 세우고 자료에 기각당한 항으로 남겨 둡니다.

---

## 문헌과 모델이 어긋나는 지점 (Documented disagreements)

1. **2 mg에서의 혈청 HGA.** 보고값 3.86 umol/L. 모델은 15.1 umol/L을 냅니다(3.9배 과대).
   이는 파라미터 문제가 아니라 질량보존의 문제입니다. HGA는 사실상 전부 소변으로
   나가므로 요중 HGA ≈ 총 HGA 생산량입니다. 2 mg에서 요중 1,200–1,440 umol/day이고
   10 mg에서 181 umol/day라면, 혈청이 3.86과 2.23이 되기 위해서는 HGA 클리어런스가
   2 mg에서 10 mg보다 **4.6배 높아야** 합니다. 그런데 경쟁 유기음이온(HPPA·HPLA)
   농도는 두 용량에서 거의 같습니다. 어떤 단일 클리어런스 모델로도 네 값을 동시에
   만족시킬 수 없습니다. 후보 설명은 (i) 낮은 농도 구간에서의 혈청 HGA 정량 한계,
   (ii) NAC(2 mg, 적극적 식이관리)와 SONIA 2(10 mg, 정보 제공만)의 코호트 차이,
   (iii) 아직 기술되지 않은 포화성 분비 경로입니다. 판별 실험: **동일 코호트에서 같은
   방문에 2 mg과 10 mg의 혈청·요중 HGA를 동시 측정**하는 것 (PMID 35028273 대상).
2. **혈청 HPPA 상승폭.** 보고 14.3–15.0배, 모델 24.9배. `HNT` = 1.45의 급한 억제
   곡선이 HPP를 필요 이상으로 밀어올립니다. 같은 1.70배 과대가 요중 HGA의
   2 mg/10 mg 비(보고 6.63, 모델 11.3)에도 나타나는데, 두 오차가 같은 배수인 것은
   우연이 아니라 둘 다 동일한 억제 급경사에 의해 결정되기 때문입니다.
3. **니티시논의 대동맥판 효과.** SONIA 2는 4년간 최대 압력차 진행률 차이
   0.0093 mmHg/년(p=0.53)으로 효과 없음을 보고했습니다. 모델은 판막 경로를
   가역 채널에서 제외해 -0.057 mmHg/년을 냅니다 — 부호는 맞고 크기는 여전히 6배
   과대하지만, 대조군의 4년 진행폭(약 4 mmHg)에 대비하면 1.4%로 실질적 null입니다.
4. **색소 침착의 역전 속도.** PMID 32904992는 니티시논 치료 후 피부·귀·공막 색소의
   가시적 역전을 보고합니다. 모델은 2년간 0.5%만 감소시킵니다(기대 약 20%). 모델의
   진피 교체 반감기 6년이 2년 만의 가시적 변화를 만들 수 없기 때문입니다. 관찰된
   역전은 **가시 색소의 유효 교체가 수개월 규모**임을 의미하며, 이는 측정 가능한
   파라미터입니다(연속 사진 계측 + 니티시논 중단 후 재침착 속도).
5. **신결석 발생 시점.** 요중 HGA 과포화는 출생 시부터 존재하므로 1차 반응
   모델로는 유아기에 결석이 생깁니다. 신 파필라 색소를 핵으로 gating하여
   59.4세(보고 64세)로 개선했지만, 모델에는 결석 핵형성 유도기(induction time)
   구조가 없어 '임상적으로 드러나는 시점'은 임의의 검출 임계값에 의존합니다.

---

## 도구 (Tools)

- **mrgsolve** — <https://mrgsolve.org> · <https://github.com/metrumresearchgroup/mrgsolve>
- **Graphviz** — <https://graphviz.org>
- **R / Shiny** — <https://shiny.posit.co>
- **PubMed E-utilities** — <https://www.ncbi.nlm.nih.gov/books/NBK25501/>
