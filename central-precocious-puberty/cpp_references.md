# 중추성 성조숙증 QSP 모델 — 참고문헌
# Central Precocious Puberty (CPP) QSP Model — References

이 문서는 `cpp_qsp_model.dot`, `cpp_mrgsolve_model.R`, `cpp_shiny_app.R`,
`cpp_reference_check.py`에 들어간 **모든 구조적 가정과 파라미터의 근거**를 정리한 것이다.

**총 97편** · 모든 PMID는 PubMed E-utilities로 조회하여 저자·저널·연도·제목을 대조 검증했다
(검증 스크립트가 반환한 서지정보를 그대로 사용했으며, 기억에 의존해 적은 PMID는 없다).

각 항목의 주석은 "이 논문이 모델의 **어느 파라미터 또는 어느 구조적 선택**을 지지하는가"를
적은 것이다. 인용만 있고 모델에 쓰이지 않은 문헌은 넣지 않았다.

---

## 1. 사춘기 개시의 유전학 — 브레이크가 풀리는 것이 시작이다 (Genetics of pubertal onset: CPP is a LOST BRAKE)

사람의 사춘기는 능동적 억제 상태로 유지되며, 단일유전자 CPP의 대부분은 구동 획득이 아니라 **억제 소실**이다. 모델에서 MKRN3(0) = `MK0`가 유일한 환자 수준 knob이고 사춘기 시작 연령은 출력이다.

1. **Abreu AP** et al. *Central precocious puberty caused by mutations in the imprinted gene MKRN3.* N Engl J Med. 2013.  
   [PMID 23738509](https://pubmed.ncbi.nlm.nih.gov/23738509/)  
   → MKRN3 기능소실이 가족성 CPP의 최다 단일유전자 원인임을 최초 보고. 모델의 `MK0` 파라미터가 존재하는 이유.

2. **Macedo DB** et al. *Central precocious puberty that appears to be sporadic caused by paternally inherited mutations in the imprinted gene makorin ring finger 3.* J Clin Endocrinol Metab. 2014.  
   [PMID 24628548](https://pubmed.ncbi.nlm.nih.gov/24628548/)  
   → '산발성'으로 보이는 CPP에서도 부계 유래 MKRN3 변이가 발견됨 — 각인 유전자이므로 모계 전달은 표현형을 만들지 않는다.

3. **Macedo DB** et al. *Central Precocious Puberty Caused by a Heterozygous Deletion in the MKRN3 Promoter Region.* Neuroendocrinology. 2018.  
   [PMID 29763903](https://pubmed.ncbi.nlm.nih.gov/29763903/)  
   → MKRN3 프로모터 영역 결손 — 코딩 변이가 아니어도 브레이크를 잃을 수 있다.

4. **Stecchini MF** et al. *Time Course of Central Precocious Puberty Development Caused by an MKRN3 Gene Mutation: A Prismatic Case.* Horm Res Paediatr. 2016.  
   [PMID 27424312](https://pubmed.ncbi.nlm.nih.gov/27424312/)  
   → MKRN3 변이 CPP의 시간 경과를 추적한 증례 — `TAUKND` 60일의 근거.

5. **Abreu AP** et al. *A new pathway in the control of the initiation of puberty: the MKRN3 gene.* J Mol Endocrinol. 2015.  
   [PMID 25957321](https://pubmed.ncbi.nlm.nih.gov/25957321/)  
   → MKRN3 경로 리뷰.

6. **Busch AS** et al. *Circulating MKRN3 Levels Decline During Puberty in Healthy Boys.* J Clin Endocrinol Metab. 2016.  
   [PMID 27057785](https://pubmed.ncbi.nlm.nih.gov/27057785/)  
   → 건강한 남아에서 순환 MKRN3가 사춘기 전에 감소함을 보임 — 모델의 `KMK` 지수 감쇠(반감기 3.0년)가 여기서 나온다.

7. **Roberts SA** et al. *Hypothalamic Overexpression of Makorin Ring Finger Protein 3 Results in Delayed Puberty in Female Mice.* Endocrinology. 2022.  
   [PMID 35974456](https://pubmed.ncbi.nlm.nih.gov/35974456/)  
   → 시상하부 MKRN3 과발현이 암컷 마우스의 사춘기를 지연시킴 — `MK0`↑ → 사춘기 지연이라는 모델 방향성의 실험적 근거.

8. **Dauber A** et al. *Paternally Inherited DLK1 Deletion Associated With Familial Central Precocious Puberty.* J Clin Endocrinol Metab. 2017.  
   [PMID 28324015](https://pubmed.ncbi.nlm.nih.gov/28324015/)  
   → DLK1 부계 결손이 가족성 CPP + 중심성 비만을 유발 — 모델의 `DLK1R` 잔여 억제항.

9. **Teles MG** et al. *A GPR54-activating mutation in a patient with central precocious puberty.* N Engl J Med. 2008.  
   [PMID 18272894](https://pubmed.ncbi.nlm.nih.gov/18272894/)  
   → KISS1R(GPR54) 활성화 변이 CPP — 구동 획득형의 반례.

10. **Silveira LG** et al. *Mutations of the KISS1 gene in disorders of puberty.* J Clin Endocrinol Metab. 2010.  
   [PMID 20237166](https://pubmed.ncbi.nlm.nih.gov/20237166/)  
   → KISS1 유전자 변이와 사춘기 장애.

11. **de Roux N** et al. *Hypogonadotropic hypogonadism due to loss of function of the KiSS1-derived peptide receptor GPR54.* Proc Natl Acad Sci U S A. 2003.  
   [PMID 12944565](https://pubmed.ncbi.nlm.nih.gov/12944565/)  
   → GPR54 기능소실 → 저성선자극호르몬성 저성선증. 동일 축의 반대 극단.

12. **Topaloglu AK** et al. *TAC3 and TACR3 mutations in familial hypogonadotropic hypogonadism reveal a key role for Neurokinin B in the central control of reproduction.* Nat Genet. 2009.  
   [PMID 19079066](https://pubmed.ncbi.nlm.nih.gov/19079066/)  
   → TAC3/TACR3 변이 — NKB가 KNDy 자가흥분(모델의 pulse ignition)에 필수임을 확립.

13. **Canton APM** et al. *Rare variants in the MECP2 gene in girls with central precocious puberty: a translational cohort study.* Lancet Diabetes Endocrinol. 2023.  
   [PMID 37385287](https://pubmed.ncbi.nlm.nih.gov/37385287/)  
   → MECP2 희귀변이와 여아 CPP.

14. **Canton APM** et al. *MECP2 Rare Variants in Boys With Central Precocious Puberty.* J Clin Endocrinol Metab. 2026.  
   [PMID 41139199](https://pubmed.ncbi.nlm.nih.gov/41139199/)  
   → MECP2 희귀변이와 남아 CPP.

15. **Day FR** et al. *Genomic analyses identify hundreds of variants associated with age at menarche and support a role for puberty timing in cancer risk.* Nat Genet. 2017.  
   [PMID 28436984](https://pubmed.ncbi.nlm.nih.gov/28436984/)  
   → 초경 연령 관련 수백 개 변이 GWAS — 사춘기 시점의 유전적 분산 비중.

16. **Abreu AP** et al. *Pubertal development and regulation.* Lancet Diabetes Endocrinol. 2016.  
   [PMID 26852256](https://pubmed.ncbi.nlm.nih.gov/26852256/)  
   → 사춘기 발달과 조절 — 축 전체를 정리한 리뷰.

17. **Roberts SA** et al. *GENETICS IN ENDOCRINOLOGY: Genetic etiologies of central precocious puberty and the role of imprinted genes.* Eur J Endocrinol. 2020.  
   [PMID 32698138](https://pubmed.ncbi.nlm.nih.gov/32698138/)  
   → CPP의 유전적 원인과 각인 유전자 리뷰.

---

## 2. KNDy 박동발생기 — 상태변수는 농도가 아니라 빈도다 (The KNDy pulse generator: the state variable is FREQUENCY)

모델은 GnRH 박동 **빈도**(`PULS`)를 명시적 상태변수로 둔다. NKB는 박동을 점화하고 dynorphin은 종료시키며, 그 결과 박동 간격이 정해진다.

18. **Goodman RL** et al. *Kisspeptin neurons in the arcuate nucleus of the ewe express both dynorphin A and neurokinin B.* Endocrinology. 2007.  
   [PMID 17823266](https://pubmed.ncbi.nlm.nih.gov/17823266/)  
   → 양의 궁상핵 kisspeptin 뉴런이 dynorphin A와 NKB를 함께 발현 — KNDy 개념의 해부학적 근거.

19. **Merkley CM** et al. *KNDy (kisspeptin/neurokinin B/dynorphin) neurons are activated during both pulsatile and surge secretion of LH in the ewe.* Endocrinology. 2012.  
   [PMID 22989631](https://pubmed.ncbi.nlm.nih.gov/22989631/)  
   → KNDy 뉴런이 박동성 분비와 LH surge 모두에서 활성화됨.

20. **Clarkson J** et al. *Definition of the hypothalamic GnRH pulse generator in mice.* Proc Natl Acad Sci U S A. 2017.  
   [PMID 29109258](https://pubmed.ncbi.nlm.nih.gov/29109258/)  
   → 마우스에서 시상하부 GnRH 박동발생기를 광유전학적으로 정의 — KNDy 집단 활동이 곧 박동임을 증명.

21. **McQuillan HJ** et al. *Definition of the estrogen negative feedback pathway controlling the GnRH pulse generator in female mice.* Nat Commun. 2022.  
   [PMID 36460649](https://pubmed.ncbi.nlm.nih.gov/36460649/)  
   → 박동발생기를 제어하는 에스트로겐 음성 피드백 경로의 정의 — 모델의 `KFB` 항.

22. **Garcia JP** et al. *Kisspeptin and Neurokinin B Signaling Network Underlies the Pubertal Increase in GnRH Release in Female Rhesus Monkeys.* Endocrinology. 2017.  
   [PMID 28977601](https://pubmed.ncbi.nlm.nih.gov/28977601/)  
   → 영장류에서 kisspeptin-NKB 신호망이 사춘기 GnRH 분비 증가의 기반임을 보임.

23. **Toro CA** et al. *Trithorax dependent changes in chromatin landscape at enhancer and promoter regions drive female puberty.* Nat Commun. 2018.  
   [PMID 29302059](https://pubmed.ncbi.nlm.nih.gov/29302059/)  
   → Trithorax 의존적 염색질 변화가 여성 사춘기를 구동 — 후성유전 스위치.

24. **Vazquez MJ** et al. *SIRT1 mediates obesity- and nutrient-dependent perturbation of pubertal timing by epigenetically controlling Kiss1 expression.* Nat Commun. 2018.  
   [PMID 30305620](https://pubmed.ncbi.nlm.nih.gov/30305620/)  
   → SIRT1이 영양/비만 의존적 사춘기 시점 교란을 매개 — 모델의 `KLEPB` 렙틴 허용 신호.

25. **Bessa DS** et al. *Methylome profiling of healthy and central precocious puberty girls.* Clin Epigenetics. 2018.  
   [PMID 30466473](https://pubmed.ncbi.nlm.nih.gov/30466473/)  
   → CPP 여아의 메틸롬 프로파일링.

26. **Teilmann G** et al. *Increased risk of precocious puberty in internationally adopted children in Denmark.* Pediatrics. 2006.  
   [PMID 16882780](https://pubmed.ncbi.nlm.nih.gov/16882780/)  
   → 국제 입양아에서 성조숙증 위험 증가 — 모델의 대사/환경 허용 신호를 뒷받침하는 역학.

---

## 3. 박동성 해독과 GnRH 수용체 약리 — 작용제는 막지 않고 '코드를 파괴'한다 (Pulsatility decoding and GnRHR pharmacology)

이 절의 문헌이 모델의 핵심 방정식 `S = RS·(Sendo·ffree + AINT·fa)`를 정당화한다. 작용제는 내인활성 `AINT` = 1.6 때문에 처음에는 자극을 **올리고**(flare), 억제는 감작 수용체 `RS`가 붕괴한 뒤에 온다.

27. **Belchetz PE** et al. *Hypophysial responses to continuous and intermittent delivery of hypopthalamic gonadotropin-releasing hormone.* Science. 1978.  
   [PMID 100883](https://pubmed.ncbi.nlm.nih.gov/100883/)  
   → Knobil 그룹의 고전 실험: 지속 투여는 성선자극호르몬을 억제하고 박동 투여는 유지시킨다. 모델 전체의 치료 원리가 이 한 논문이다.

28. **Wildt L** et al. *Frequency and amplitude of gonadotropin-releasing hormone stimulation and gonadotropin secretion in the rhesus monkey.* Endocrinology. 1981.  
   [PMID 6788538](https://pubmed.ncbi.nlm.nih.gov/6788538/)  
   → 붉은털원숭이에서 GnRH 자극의 빈도와 진폭이 성선자극호르몬 분비를 결정 — `PULSMIN`/`PULSMAX` 범위의 출처.

29. **Conn PM** et al. *Gonadotropin-releasing hormone and its analogues.* N Engl J Med. 1991.  
   [PMID 1984190](https://pubmed.ncbi.nlm.nih.gov/1984190/)  
   → GnRH와 그 유사체 — 임상 약리 총설.

30. **Millar RP** et al. *Gonadotropin-releasing hormone receptors.* Endocr Rev. 2004.  
   [PMID 15082521](https://pubmed.ncbi.nlm.nih.gov/15082521/)  
   → GnRH 수용체 총설. II형 수용체가 C-말단 꼬리를 갖지 않아 재순환이 느리고 불완전한 점이 `KREC` = 0.035/일의 근거.

31. **Sealfon SC** et al. *Molecular mechanisms of ligand interaction with the gonadotropin-releasing hormone receptor.* Endocr Rev. 1997.  
   [PMID 9101136](https://pubmed.ncbi.nlm.nih.gov/9101136/)  
   → GnRH 수용체-리간드 상호작용의 분자기전 — `EC50` 값들의 근거.

32. **Kaiser UB** et al. *Differential effects of gonadotropin-releasing hormone (GnRH) pulse frequency on gonadotropin subunit and GnRH receptor messenger ribonucleic acid levels in vitro.* Endocrinology. 1997.  
   [PMID 9048630](https://pubmed.ncbi.nlm.nih.gov/9048630/)  
   → GnRH 박동 빈도가 성선자극호르몬 아단위와 수용체 mRNA에 미치는 차별적 효과 — 빠른 박동은 LHB, 느린 박동은 FSHB. 모델의 pulsatility decoder.

33. **Kanasaki H** et al. *Gonadotropin-releasing hormone pulse frequency-dependent activation of extracellular signal-regulated kinase pathways in perifused LbetaT2 cells.* Endocrinology. 2005.  
   [PMID 16141398](https://pubmed.ncbi.nlm.nih.gov/16141398/)  
   → GnRH 박동 빈도 의존적 ERK 경로 활성화.

34. **Tornøe CW** et al. *Pharmacokinetic/pharmacodynamic modelling of GnRH antagonist degarelix: a comparison of the non-linear mixed-effects programs NONMEM and NLME.* J Pharmacokinet Pharmacodyn. 2004.  
   [PMID 16222784](https://pubmed.ncbi.nlm.nih.gov/16222784/)  
   → GnRH 길항제 degarelix의 PK/PD 모델링 — 모델의 길항제 팔(`KAX`, `EC50X`)의 구조적 참고.

35. **Tornøe CW** et al. *Population pharmacokinetic/pharmacodynamic (PK/PD) modelling of the hypothalamic-pituitary-gonadal axis following treatment with GnRH analogues.* Br J Clin Pharmacol. 2007.  
   [PMID 17096678](https://pubmed.ncbi.nlm.nih.gov/17096678/)  
   → degarelix 투여 후 시상하부-뇌하수체-성선 축의 집단 PK/PD 모델링 — 억제 동역학의 정량적 참고.

---

## 4. 성장판 — 오직 에스트로겐이 판을 닫고, 그것은 시간이 아니라 '소모된 예비능'이다 (The growth plate: oestrogen closes it, by DEPLETING a reserve)

모델의 (−) 팔. 인간 기능소실 실험(ESR1, 아로마타제 결핍)이 두 성 모두에서 골성숙 종료가 **에스트로겐 의존적**임을 확정했고, Baron 그룹이 그것을 시계가 아닌 예비능 소모로 정식화했다. 모델의 `GPRES`가 이것이다.

36. **Morishima A** et al. *Aromatase deficiency in male and female siblings caused by a novel mutation and the physiological role of estrogens.* J Clin Endocrinol Metab. 1995.  
   [PMID 8530621](https://pubmed.ncbi.nlm.nih.gov/8530621/)  
   → 남녀 형제의 아로마타제 결핍 — 에스트로겐 없이는 성장판이 닫히지 않는다.

37. **Carani C** et al. *Effect of testosterone and estradiol in a man with aromatase deficiency.* N Engl J Med. 1997.  
   [PMID 9211678](https://pubmed.ncbi.nlm.nih.gov/9211678/)  
   → 아로마타제 결핍 남성에서 테스토스테론과 에스트라디올의 효과 비교 — 골성숙을 닫는 것은 E2이며 T가 아니다. 모델이 남아에서 아로마타제 억제제를 '두 부호를 분리하는 도구'로 쓸 수 있는 근거.

38. **Nilsson O** et al. *Fundamental limits on longitudinal bone growth: growth plate senescence and epiphyseal fusion.* Trends Endocrinol Metab. 2004.  
   [PMID 15380808](https://pubmed.ncbi.nlm.nih.gov/15380808/)  
   → 종축 성장의 근본적 한계: 성장판 노화와 골성숙 종료. `GPRES` 1 → 0의 개념적 출처.

39. **Weise M** et al. *Effects of estrogen on growth plate senescence and epiphyseal fusion.* Proc Natl Acad Sci U S A. 2001.  
   [PMID 11381135](https://pubmed.ncbi.nlm.nih.gov/11381135/)  
   → 에스트로겐이 성장판 노화와 골성숙 종료를 가속함을 실험적으로 보임.

40. **Nilsson O** et al. *Evidence that estrogen hastens epiphyseal fusion and cessation of longitudinal bone growth by irreversibly depleting the number of resting zone progenitor cells in female rabbits.* Endocrinology. 2014.  
   [PMID 24708243](https://pubmed.ncbi.nlm.nih.gov/24708243/)  
   → 에스트로겐이 성장판 전구세포 예비능을 **비가역적으로 소모**하여 골성숙을 앞당긴다는 직접 증거 — 모델에서 `GPRES`가 절대 회복되지 않는 이유.

41. **Emons J** et al. *Mechanisms of growth plate maturation and epiphyseal fusion.* Horm Res Paediatr. 2011.  
   [PMID 21540578](https://pubmed.ncbi.nlm.nih.gov/21540578/)  
   → 성장판 성숙과 골성숙 종료의 기전 총설.

42. **Grumbach MM** et al. *Estrogen, bone, growth and sex: a sea change in conventional wisdom.* J Pediatr Endocrinol Metab. 2000.  
   [PMID 11202221](https://pubmed.ncbi.nlm.nih.gov/11202221/)  
   → 에스트로겐·뼈·성장 — 통설의 전환을 정리한 논평.

43. **Emons J** et al. *Expression of vascular endothelial growth factor in the growth plate is stimulated by estradiol and increases during pubertal development.* J Endocrinol. 2010.  
   [PMID 20093283](https://pubmed.ncbi.nlm.nih.gov/20093283/)  
   → 성장판 VEGF 발현이 에스트라디올에 의해 자극되고 사춘기에 증가.

44. **BAYLEY N** et al. *Tables for predicting adult height from skeletal age: revised for use with the Greulich-Pyle hand standards.* J Pediatr. 1952.  
   [PMID 14918032](https://pubmed.ncbi.nlm.nih.gov/14918032/)  
   → Bayley-Pinneau 표 — 임상에서 계산하는 예측 성인키. 모델은 이 표를 별도로 두어 '임상 예측치'와 '모델이 실제 적분한 최종키'를 나란히 비교한다.

---

## 5. 진료 지침과 최종 성인키 결과 (Guidelines and adult-height outcomes)

모델의 사인 플립 결과를 검증하는 기준 문헌들. 합의문은 골연령 12세 이후 시작 시 신장 이득이 없다고 기술하며, 모델은 그 지점을 골연령 11.6년(이득 1 cm 미만)으로 **계산해서** 재현한다.

45. **Carel JC** et al. *Consensus statement on the use of gonadotropin-releasing hormone analogs in children.* Pediatrics. 2009.  
   [PMID 19332438](https://pubmed.ncbi.nlm.nih.gov/19332438/)  
   → 소아에서 GnRH 유사체 사용에 관한 국제 합의문 — 치료 적응증, 모니터링, 중단 시점의 기준.

46. **Carel JC** et al. *Clinical practice. Precocious puberty.* N Engl J Med. 2008.  
   [PMID 18509122](https://pubmed.ncbi.nlm.nih.gov/18509122/)  
   → 성조숙증 임상 총설.

47. **Latronico AC** et al. *Causes, diagnosis, and treatment of central precocious puberty.* Lancet Diabetes Endocrinol. 2016.  
   [PMID 26852255](https://pubmed.ncbi.nlm.nih.gov/26852255/)  
   → 중추성 성조숙증의 원인·진단·치료.

48. **Bangalore Krishna K** et al. *Use of Gonadotropin-Releasing Hormone Analogs in Children: Update by an International Consortium.* Horm Res Paediatr. 2019.  
   [PMID 31319416](https://pubmed.ncbi.nlm.nih.gov/31319416/)  
   → GnRH 유사체 사용에 관한 국제 컨소시엄 업데이트.

49. **Latronico AC** et al. *Central precocious puberty: an Endocrine Society clinical practice guideline.* J Clin Endocrinol Metab. 2026.  
   [PMID 42287186](https://pubmed.ncbi.nlm.nih.gov/42287186/)  
   → 중추성 성조숙증 — Endocrine Society 임상진료지침(최신).

50. **Carel JC** et al. *Final height after long-term treatment with triptorelin slow release for central precocious puberty: importance of statural growth after interruption of treatment. French study group of Decapeptyl in Precocious Puberty.* J Clin Endocrinol Metab. 1999.  
   [PMID 10372696](https://pubmed.ncbi.nlm.nih.gov/10372696/)  
   → 트립토렐린 서방형 장기 치료 후 최종 신장 — 치료군 최종키의 주요 기준값.

51. **Klein KO** et al. *Increased final height in precocious puberty after long-term treatment with LHRH agonists: the National Institutes of Health experience.* J Clin Endocrinol Metab. 2001.  
   [PMID 11600530](https://pubmed.ncbi.nlm.nih.gov/11600530/)  
   → LHRH 작용제 장기 치료 후 최종 신장 증가: NIH 경험.

52. **Heger S** et al. *Long-term outcome after depot gonadotropin-releasing hormone agonist treatment of central precocious puberty: final height, body proportions, body composition, bone mineral density, and reproductive function.* J Clin Endocrinol Metab. 1999.  
   [PMID 10599723](https://pubmed.ncbi.nlm.nih.gov/10599723/)  
   → 데포 GnRH 작용제 치료 후 장기 결과 — 최종 신장, 체형 비율, 체성분, 골밀도.

53. **Lazar L** et al. *Growth pattern and final height after cessation of gonadotropin-suppressive therapy in girls with central sexual precocity.* J Clin Endocrinol Metab. 2007.  
   [PMID 17579199](https://pubmed.ncbi.nlm.nih.gov/17579199/)  
   → 성선억제 치료 중단 후 성장 패턴과 최종 신장 — 모델의 치료 종료 후 회복·초경 시점 검증.

54. **Pasquino AM** et al. *Long-term observation of 87 girls with idiopathic central precocious puberty treated with gonadotropin-releasing hormone analogs: impact on adult height, body mass index, bone mineral content, and reproductive function.* J Clin Endocrinol Metab. 2008.  
   [PMID 17940112](https://pubmed.ncbi.nlm.nih.gov/17940112/)  
   → 특발성 CPP 여아 87명의 장기 관찰 — 최종 신장, BMI, 난소 기능.

55. **Antoniazzi F** et al. *End results in central precocious puberty with GnRH analog treatment: the data of the Italian Study Group for Physiopathology of Puberty.* J Pediatr Endocrinol Metab. 2000.  
   [PMID 10969920](https://pubmed.ncbi.nlm.nih.gov/10969920/)  
   → 이탈리아 연구그룹의 CPP GnRH 유사체 치료 최종 결과.

56. **Arrigo T** et al. *When to stop GnRH analog therapy: the experience of the Italian Study Group for Physiopathology of Puberty.* J Pediatr Endocrinol Metab. 2000.  
   [PMID 10969918](https://pubmed.ncbi.nlm.nih.gov/10969918/)  
   → GnRH 유사체를 언제 중단할 것인가 — 이탈리아 연구그룹.

57. **Bereket A** et al. *A Critical Appraisal of the Effect of Gonadotropin-Releasing Hormon Analog Treatment on Adult Height of Girls with Central Precocious Puberty.* J Clin Res Pediatr Endocrinol. 2017.  
   [PMID 29280737](https://pubmed.ncbi.nlm.nih.gov/29280737/)  
   → GnRH 유사체 치료가 최종 신장에 미치는 효과에 대한 비판적 평가 — 이득 크기가 시작 시점에 강하게 의존한다는 점을 강조.

58. **Guaraldi F** et al. *MANAGEMENT OF ENDOCRINE DISEASE: Long-term outcomes of the treatment of central precocious puberty.* Eur J Endocrinol. 2016.  
   [PMID 26466612](https://pubmed.ncbi.nlm.nih.gov/26466612/)  
   → 중추성 성조숙증 치료의 장기 결과 총설.

59. **Magiakou MA** et al. *The efficacy and safety of gonadotropin-releasing hormone analog treatment in childhood and adolescence: a single center, long-term follow-up study.* J Clin Endocrinol Metab. 2010.  
   [PMID 19897682](https://pubmed.ncbi.nlm.nih.gov/19897682/)  
   → 소아·청소년기 GnRH 유사체 치료의 유효성과 안전성 — 단일기관 장기 추적.

60. **Hwangbo J** et al. *Long-term outcomes of gonadotropin-releasing hormone agonist treatment in girls with central precocious puberty.* Ann Pediatr Endocrinol Metab. 2025.  
   [PMID 40049673](https://pubmed.ncbi.nlm.nih.gov/40049673/)  
   → CPP 여아에서 GnRH 작용제 치료의 장기 결과(최근 코호트).

61. **Cho KW** et al. *Final adult height in male patients with central precocious puberty after gonadotropin-releasing hormone agonist treatment.* Ann Pediatr Endocrinol Metab. 2026.  
   [PMID 41787710](https://pubmed.ncbi.nlm.nih.gov/41787710/)  
   → 남아 CPP에서 GnRH 작용제 치료 후 최종 성인키.

---

## 6. 제형·용량·PK — 결정 변수는 효능이 아니라 최저농도 커버리지다 (Formulation and PK: trough coverage is the design variable)

모델은 정확히 투여된 데포 제형들 사이에 거의 차이가 없고, 실패는 (a) 비강분무, (b) 지연된 주사, (c) 저용량에서 생긴다고 계산한다.

62. **Eugster EA** et al. *Efficacy and safety of histrelin subdermal implant in children with central precocious puberty: a multicenter trial.* J Clin Endocrinol Metab. 2007.  
   [PMID 17327379](https://pubmed.ncbi.nlm.nih.gov/17327379/)  
   → 히스트렐린 피하 임플란트의 유효성·안전성 다기관 시험 — 0차 방출 65 µg/일, 최저농도가 없는 제형.

63. **Silverman LA** et al. *Long-Term Continuous Suppression With Once-Yearly Histrelin Subcutaneous Implants for the Treatment of Central Precocious Puberty: A Final Report of a Phase 3 Multicenter Trial.* J Clin Endocrinol Metab. 2015.  
   [PMID 25803268](https://pubmed.ncbi.nlm.nih.gov/25803268/)  
   → 연 1회 히스트렐린 임플란트의 장기 지속 억제.

64. **Fuld K** et al. *A randomized trial of 1- and 3-month depot leuprolide doses in the treatment of central precocious puberty.* J Pediatr. 2011.  
   [PMID 21798557](https://pubmed.ncbi.nlm.nih.gov/21798557/)  
   → 1개월 vs 3개월 데포 류프롤리드 무작위 비교 — 모델이 두 제형을 사실상 구분하지 못한다는 예측의 검증 대상.

65. **Klein K** et al. *Efficacy and safety of triptorelin 6-month formulation in patients with central precocious puberty.* J Pediatr Endocrinol Metab. 2016.  
   [PMID 26887034](https://pubmed.ncbi.nlm.nih.gov/26887034/)  
   → 트립토렐린 6개월 제형의 유효성과 안전성.

66. **Mostafa NM** et al. *Pharmacokinetic and exposure-response analyses of leuprolide following administration of leuprolide acetate 3-month depot formulations to children with central precocious puberty.* Clin Drug Investig. 2014.  
   [PMID 24756362](https://pubmed.ncbi.nlm.nih.gov/24756362/)  
   → 류프롤리드 3개월 데포의 PK 및 노출-반응 분석 — `FBURSTL`, `KDISL`, `CLL` 보정의 출처.

67. **Lim CN** et al. *Fixed Dosing of Leuprolide Acetate, a GnRH Agonist, in Children with Central Precocious Puberty: A Population Pharmacokinetic Justification.* Paediatr Drugs. 2026.  
   [PMID 41824265](https://pubmed.ncbi.nlm.nih.gov/41824265/)  
   → 소아 CPP에서 류프롤리드 고정용량 집단 PK 분석.

68. **Wang X** et al. *Reflections on GnRH analogues dosage in the context of escape phenomenon: a case report of hypothalamic hamartoma with central precocious puberty and literature review.* J Pediatr Endocrinol Metab. 2026.  
   [PMID 42381421](https://pubmed.ncbi.nlm.nih.gov/42381421/)  
   → 탈출(escape) 현상 맥락에서의 GnRH 유사체 용량 재고 — 시상하부 하마르토마 증례. 모델의 `KNDLES` ectopic drive 항.

---

## 7. 진단과 치료 중 모니터링 (Diagnosis and on-treatment monitoring)

모델의 모니터링 탭이 검증하려는 대상. **성장속도는 효과적 억제에서 오히려 떨어지므로 신뢰할 수 없고**, ΔBA/ΔCA와 기저 LH·자궁 부피가 작동한다.

69. **Neely EK** et al. *Normal ranges for immunochemiluminometric gonadotropin assays.* J Pediatr. 1995.  
   [PMID 7608809](https://pubmed.ncbi.nlm.nih.gov/7608809/)  
   → 면역화학발광법 성선자극호르몬 정상범위 — 모델 LH 보정의 앵커(사춘기전 0.1-0.15 IU/L).

70. **Houk CP** et al. *Adequacy of a single unstimulated luteinizing hormone level to diagnose central precocious puberty in girls.* Pediatrics. 2009.  
   [PMID 19482738](https://pubmed.ncbi.nlm.nih.gov/19482738/)  
   → 단일 무자극 LH로 여아 CPP를 진단할 수 있는가 — 특이도는 높고 민감도는 낮다.

71. **Resende EA** et al. *Assessment of basal and gonadotropin-releasing hormone-stimulated gonadotropins by immunochemiluminometric and immunofluorometric assays in normal children.* J Clin Endocrinol Metab. 2007.  
   [PMID 17284632](https://pubmed.ncbi.nlm.nih.gov/17284632/)  
   → 기저 및 GnRH 자극 성선자극호르몬의 측정법 비교 — 자극 후 최고 LH > 5 IU/L 기준의 근거.

72. **de Vries L** et al. *Role of pelvic ultrasound in girls with precocious puberty.* Horm Res Paediatr. 2011.  
   [PMID 21228561](https://pubmed.ncbi.nlm.nih.gov/21228561/)  
   → 성조숙증 여아에서 골반 초음파의 역할 — 자궁 부피/길이 기준(> 3.4 mL / > 34 mm).

73. **Chalumeau M** et al. *Selecting girls with precocious puberty for brain imaging: validation of European evidence-based diagnosis rule.* J Pediatr. 2003.  
   [PMID 14571217](https://pubmed.ncbi.nlm.nih.gov/14571217/)  
   → 성조숙증 여아에서 뇌 영상 대상 선별 — 유럽 근거기반 진단 규칙의 검증.

74. **Cho MH** et al. *Diagnostic Usefulness of Subcutaneous Triptorelin Stimulation Testing in the Evaluation of Central Precocious Puberty in Girls.* Horm Metab Res. 2026.  
   [PMID 41956123](https://pubmed.ncbi.nlm.nih.gov/41956123/)  
   → 피하 트립토렐린 자극검사의 진단적 유용성.

---

## 8. 골밀도·체성분·장기 안전성 (Bone mass, body composition, long-term safety)

모델은 골밀도 Z-점수의 '하강'이 대부분 **진행된 기저치가 같은 나이 기준으로 회귀하는 것**이며 최대 골질량은 회복된다고 계산한다.

75. **Boot AM** et al. *Bone mineral density and body composition before and during treatment with gonadotropin-releasing hormone agonist in children with central precocious and early puberty.* J Clin Endocrinol Metab. 1998.  
   [PMID 9467543](https://pubmed.ncbi.nlm.nih.gov/9467543/)  
   → CPP 소아에서 GnRH 작용제 치료 전·중의 골밀도와 체성분 — `KZ`, `KCATCH` 보정의 출처.

76. **van der Sluis IM** et al. *Longitudinal follow-up of bone density and body composition in children with precocious or early puberty before, during and after cessation of GnRH agonist therapy.* J Clin Endocrinol Metab. 2002.  
   [PMID 11836277](https://pubmed.ncbi.nlm.nih.gov/11836277/)  
   → 성조숙증/조기 사춘기 소아의 골밀도·체성분 종단 추적.

77. **Antoniazzi F** et al. *Prevention of bone demineralization by calcium supplementation in precocious puberty during gonadotropin-releasing hormone agonist treatment.* J Clin Endocrinol Metab. 1999.  
   [PMID 10372699](https://pubmed.ncbi.nlm.nih.gov/10372699/)  
   → 칼슘 보충이 GnRH 작용제 치료 중 골 탈미네랄화를 예방 — 모델의 `CAVD`/`KCAVD` 항.

78. **Heger S** et al. *Long-term GnRH agonist treatment for female central precocious puberty does not impair reproductive function.* Mol Cell Endocrinol. 2006.  
   [PMID 16757104](https://pubmed.ncbi.nlm.nih.gov/16757104/)  
   → 여아 CPP의 장기 GnRH 작용제 치료가 생식 기능을 손상시키지 않음.

79. **Franceschi R** et al. *Prevalence of polycystic ovary syndrome in young women who had idiopathic central precocious puberty.* Fertil Steril. 2010.  
   [PMID 19135667](https://pubmed.ncbi.nlm.nih.gov/19135667/)  
   → 특발성 CPP였던 젊은 여성의 다낭성 난소 증후군 유병률.

80. **Schoelwer MJ** et al. *Psychological assessment of mothers and their daughters at the time of diagnosis of precocious puberty.* Int J Pediatr Endocrinol. 2015.  
   [PMID 25780366](https://pubmed.ncbi.nlm.nih.gov/25780366/)  
   → 성조숙증 진단 시점의 어머니와 딸에 대한 심리평가 — 모델의 심리사회 지표(`QOL`) 항의 근거.

---

## 9. GnRH-비의존성(말초성) 성조숙증 — 작용제가 구조적으로 도달할 수 없는 영역 (GnRH-INDEPENDENT precocious puberty)

모델의 가장 선명한 결과 중 하나: 동일한 류프롤리드 요법이 중추성 질환에서는 +9.8 cm, McCune-Albright에서는 **+0.1 cm**를 준다. 표적이 기전과 맞지 않으면 약효는 0이다.

81. **Weinstein LS** et al. *Activating mutations of the stimulatory G protein in the McCune-Albright syndrome.* N Engl J Med. 1991.  
   [PMID 1944469](https://pubmed.ncbi.nlm.nih.gov/1944469/)  
   → McCune-Albright 증후군의 Gsα 활성화 변이 — LH 없이 cAMP가 켜지는 기전. 모델의 `AUTSET` 항.

82. **Shenker A** et al. *Severe endocrine and nonendocrine manifestations of the McCune-Albright syndrome associated with activating mutations of stimulatory G protein GS.* J Pediatr. 1993.  
   [PMID 8410501](https://pubmed.ncbi.nlm.nih.gov/8410501/)  
   → 활성화 Gsα 변이와 연관된 McCune-Albright의 중증 내분비·비내분비 표현형.

83. **Eugster EA** et al. *Tamoxifen treatment for precocious puberty in McCune-Albright syndrome: a multicenter trial.* J Pediatr. 2003.  
   [PMID 12915825](https://pubmed.ncbi.nlm.nih.gov/12915825/)  
   → McCune-Albright 성조숙증의 타목시펜 치료 다기관 시험 — 모델의 `TAMEFF` 항.

84. **Feuillan P** et al. *Letrozole treatment of precocious puberty in girls with the McCune-Albright syndrome: a pilot study.* J Clin Endocrinol Metab. 2007.  
   [PMID 17405850](https://pubmed.ncbi.nlm.nih.gov/17405850/)  
   → McCune-Albright 여아 성조숙증의 레트로졸 치료 — 모델에서 고효능 AI가 +8.8 cm를 내는 근거.

85. **Sims EK** et al. *Fulvestrant treatment of precocious puberty in girls with McCune-Albright syndrome.* Int J Pediatr Endocrinol. 2012.  
   [PMID 22999294](https://pubmed.ncbi.nlm.nih.gov/22999294/)  
   → McCune-Albright 여아에서 fulvestrant 치료.

86. **Reiter EO** et al. *Bicalutamide plus anastrozole for the treatment of gonadotropin-independent precocious puberty in boys with testotoxicosis: a phase II, open-label pilot study (BATT).* J Pediatr Endocrinol Metab. 2010.  
   [PMID 21158211](https://pubmed.ncbi.nlm.nih.gov/21158211/)  
   → testotoxicosis 남아의 gonadotropin-비의존성 성조숙증에 대한 bicalutamide + anastrozole.

---

## 10. 아로마타제 억제제와 rhGH — 두 부호를 분리하려는 시도 (Aromatase inhibitors and rhGH: trying to separate the two signs)

AI는 E2를 차단하면서 테스토스테론을 남기므로 원리적으로는 (−) 팔만 제거한다. 그런데 모델은 AI 단독이 남아에서 신장을 개선하지 못한다고 계산하며, 이는 실제 시험 결과와 일치한다 — E2를 없애면 GH/IGF-1 증폭, 즉 (+) 팔도 함께 사라지기 때문이다.

87. **Hero M** et al. *Inhibition of estrogen biosynthesis with a potent aromatase inhibitor increases predicted adult height in boys with idiopathic short stature: a randomized controlled trial.* J Clin Endocrinol Metab. 2005.  
   [PMID 16189252](https://pubmed.ncbi.nlm.nih.gov/16189252/)  
   → 강력한 아로마타제 억제로 특발성 저신장 남아의 예측 성인키가 증가함.

88. **Wickman S** et al. *A specific aromatase inhibitor and potential increase in adult height in boys with delayed puberty: a randomised controlled trial.* Lancet. 2001.  
   [PMID 11403810](https://pubmed.ncbi.nlm.nih.gov/11403810/)  
   → 사춘기 지연 남아에서 아로마타제 억제제와 성인키 — 무작위 시험.

89. **Varimo T** et al. *Letrozole Monotherapy in Pre- and Early-Pubertal Boys Does Not Increase Adult Height.* Front Endocrinol (Lausanne). 2019.  
   [PMID 31024444](https://pubmed.ncbi.nlm.nih.gov/31024444/)  
   → **레트로졸 단독 요법은 사춘기 전·초기 남아의 성인키를 증가시키지 않는다** — 예측 신장 개선이 실제 신장으로 이어지지 않았다. 모델이 AI 단독에서 −0.3 cm를 계산하는 것과 방향이 일치하는, 사후에 확인된 예측.

90. **Mauras N** et al. *Anastrozole increases predicted adult height of short adolescent males treated with growth hormone: a randomized, placebo-controlled, multicenter trial for one to three years.* J Clin Endocrinol Metab. 2008.  
   [PMID 18165285](https://pubmed.ncbi.nlm.nih.gov/18165285/)  
   → 성장호르몬 치료 중인 저신장 남아에서 아나스트로졸이 예측 성인키를 증가시킴.

91. **Pasquino AM** et al. *Adult height in girls with central precocious puberty treated with gonadotropin-releasing hormone analogues and growth hormone.* J Clin Endocrinol Metab. 1999.  
   [PMID 10022399](https://pubmed.ncbi.nlm.nih.gov/10022399/)  
   → CPP 여아에서 GnRH 유사체 + 성장호르몬 병용 시 성인키 — rhGH 추가 이득의 크기(수 cm) 기준.

92. **Pasquino AM** et al. *Adult height in short normal girls treated with gonadotropin-releasing hormone analogs and growth hormone.* J Clin Endocrinol Metab. 2000.  
   [PMID 10690865](https://pubmed.ncbi.nlm.nih.gov/10690865/)  
   → 정상 저신장 여아에서 GnRH 유사체 + 성장호르몬.

---

## 11. 역학과 세속적 추세 (Epidemiology and secular trends)

모델의 정상 기준 아동(유방발육 10.3세, 초경 13.0세)이 어느 인구집단을 대표하는지 정하는 문헌들.

93. **Teilmann G** et al. *Prevalence and incidence of precocious pubertal development in Denmark: an epidemiologic study based on national registries.* Pediatrics. 2005.  
   [PMID 16322154](https://pubmed.ncbi.nlm.nih.gov/16322154/)  
   → 덴마크의 성조숙증 유병률과 발생률 — 국가 등록 기반.

94. **Soriano-Guillén L** et al. *Central precocious puberty in children living in Spain: incidence, prevalence, and influence of adoption and immigration.* J Clin Endocrinol Metab. 2010.  
   [PMID 20554707](https://pubmed.ncbi.nlm.nih.gov/20554707/)  
   → 스페인 거주 소아의 CPP 발생률·유병률과 입양·이민의 영향.

95. **Bräuner EV** et al. *Trends in the Incidence of Central Precocious Puberty and Normal Variant Puberty Among Children in Denmark, 1998 to 2017.* JAMA Netw Open. 2020.  
   [PMID 33044548](https://pubmed.ncbi.nlm.nih.gov/33044548/)  
   → 덴마크 1998-2017년 중추성 성조숙증과 정상변이 사춘기의 발생률 추세.

96. **Eckert-Lind C** et al. *Worldwide Secular Trends in Age at Pubertal Onset Assessed by Breast Development Among Girls: A Systematic Review and Meta-analysis.* JAMA Pediatr. 2020.  
   [PMID 32040143](https://pubmed.ncbi.nlm.nih.gov/32040143/)  
   → 여아 유방발육 개시 연령의 전 세계 세속적 추세 — 체계적 문헌고찰 및 메타분석. 모델의 `KMK` 보정 목표(유방발육 10.3세)의 출처.

97. **Stagi S** et al. *Increased incidence of precocious and accelerated puberty in females during and after the Italian lockdown for the coronavirus 2019 (COVID-19) pandemic.* Ital J Pediatr. 2020.  
   [PMID 33148304](https://pubmed.ncbi.nlm.nih.gov/33148304/)  
   → 이탈리아 COVID-19 봉쇄 기간 중·후 여아 성조숙증 및 가속 사춘기 발생 증가.

---

## 문헌으로 검증되지 않은 부분 (Explicitly NOT evidence-based)

정직성을 위해 명시한다. 아래 항목들은 위 문헌에서 직접 유도되지 않았다.

- **`ENDOBLEED` = 0.58** — 모델 전체에서 **단 하나** 문헌 수치에 맞춰 피팅한 파라미터다.
  첫 데포 후 이탈혈(withdrawal bleeding) 발생률이 보고된 5-10% 범위에 들어오도록 정했다.
- **심리사회 지표(`QOL`)와 열성홍조 지표(`HF`)** — 0-10 척도의 합성 지수이며, 어떤
  검증된 측정도구에도 대응하지 않는다. 방향과 상대적 크기만 해석해야 한다.
- **`KMKSEX` = 0.75 (남아에서 브레이크가 더 늦게 풀림)** — 남아 사춘기가 여아보다
  약 1.8년 늦다는 관찰을 재현하기 위한 현상학적 스케일이며, MKRN3 자체의 성별 차이를
  보고한 문헌에 근거한 것이 아니다.
- **GnRH 길항제 요법(18 mg q28d)** — CPP에 승인된 길항제 제형은 없다. 모델의 길항제
  팔은 "flare 없는 억제"의 반사실(counterfactual)을 계산하기 위한 **가상 제형**이며,
  PK는 degarelix 문헌(PMID 16222784, 17096678)의 구조를 빌려왔다.
- **`XTRA` = 0.10** — 골연령 아틀라스가 읽는 성숙도보다 예비능이 약간 더 빨리 소모된다는
  가정. 방향은 PMID 24708243과 일치하지만 크기는 정해진 값이 없다.

---

## 모델이 문헌과 어긋나는 지점 (Where this model DISAGREES)

`README.md`의 마지막 절에 정량적으로 정리해 두었다. 요약하면 세 가지다.

1. 정상 아동의 최대 성장속도가 인구 자료보다 약 1 cm/년 낮다(성장판 예비능 지수가
   최종키를 맞추면 성장 급증기를 평탄하게 만든다).
2. 조기 시작 치료에서 ΔBA/ΔCA가 약 0.95-1.00으로, 흔히 인용되는 0.4-0.7보다 높다.
   늦게 시작한 경우에는 0.56-0.72로 떨어져 문헌 범위에 들어온다.
3. Bayley-Pinneau 예측 성인키가 모델이 실제로 적분한 최종키보다 5-10 cm 높게 나온다.
   이것은 버그가 아니라 모델의 예측이며, BP식이 골연령이 진행된 아동에서 과대추정한다는
   임상 관찰과 방향이 같다.

