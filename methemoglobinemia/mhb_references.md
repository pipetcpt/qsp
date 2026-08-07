# 메트헤모글로빈혈증 QSP 모델 — 참고문헌
## Methaemoglobinaemia · References

---

## 이 목록이 만들어진 방법 (How this list was produced)

이 저장소의 다른 모델들과 마찬가지로, **PMID를 기억에 의존해 쓰지 않았습니다.**
각 주제별 검색어를 PubMed E-utilities(`esearch` → `esummary`)에 직접 질의해
반환된 실제 레코드에서 PMID·제목·저널·연도를 가져왔습니다. 따라서 아래 링크는
전부 실재하는 PubMed 레코드입니다. (검색어 자체는 각 항목 끝에 `⌕`로 남겨
두었으므로, 왜 그 논문이 이 절에 있는지 추적할 수 있습니다.)

각 절 머리의 문단은 **그 절이 모델의 어느 파라미터 또는 어느 주장을 떠받치는지**를
적은 것입니다. 문헌 목록과 방정식을 따로 읽지 않아도 되게 하기 위한 것입니다.

---

## 근거의 강도에 대한 정직한 진술 (Where the evidence is weak)

이 모델의 헤드라인 결과 — "MetHb 30%는 Hb 10.5가 아니라 Hb 8.3에 해당한다" —
는 **C절**의 알로스테릭 좌측 이동(`ALPHAM`)에 전적으로 의존합니다. 그런데
메트헤모글로빈이 남은 아철(Fe²⁺) 소단위의 P50을 정확히 얼마나 낮추는지에 대한
**정량적** 인체 데이터는 놀라울 만큼 빈약합니다. Darling & Roughton(1942)의
고전적 관찰과 원자가 혼성(valency hybrid) 연구들이 방향과 대략의 크기를 지지하지만,
`ALPHAM = 0.50`은 그 문헌들과 정합적인 값일 뿐 그로부터 유일하게 결정되는 값은
아닙니다.

`mhb_reference_check.py --all`의 민감도 절은 이것을 숨기지 않고 계산해 보여 줍니다:

| 좌측 이동 크기 | MetHb 30%의 등가 Hb |
|---|---|
| 0 (이동 없음) | 10.50 g/dL — 즉 통상적 산술이 옳고, 이 모델의 주장은 사라진다 |
| 절반 | 약 9.4 g/dL — 주장은 약해지지만 방향은 유지된다 |
| 모델값 | 8.31 g/dL |

**따라서 이 모델은 반증 가능합니다.** 필요한 실험은 새로운 것이 아닙니다: 같은
검체에서 동시산소측정 %MetHb와 실측 P50을 함께 재고, 둘의 관계를 보면 됩니다.
그 기울기가 0이면 이 모델의 중심 주장은 틀린 것입니다.

---


## A. 개관 · 역학 · 진단

> 모델의 임상적 좌표. 후향 코호트들이 반복해서 보여주는 것: 보고되는 %MetHb와 실제 중증도의 상관은 생각보다 약하며, 그 이유가 이 모델의 주제다.

1. **Iolascon A.** Recommendations for diagnosis and treatment of methemoglobinemia. *Am J Hematol* 2021. [PMID 34467556](https://pubmed.ncbi.nlm.nih.gov/34467556/)  
   <sub>⌕ `methemoglobinemia review clinical diagnosis management`</sub>
2. **Gangat N.** JAK2 Unmutated Erythrocytosis: 2026 Update on Diagnosis and Management. *Am J Hematol* 2025. [PMID 41123216](https://pubmed.ncbi.nlm.nih.gov/41123216/)  
   <sub>⌕ `methemoglobinemia review clinical diagnosis management`</sub>
3. **Gangat N.** JAK2 unmutated erythrocytosis: 2023 Update on diagnosis and management. *Am J Hematol* 2023. [PMID 36966432](https://pubmed.ncbi.nlm.nih.gov/36966432/)  
   <sub>⌕ `methemoglobinemia review clinical diagnosis management`</sub>
4. **Alagha I.** Methemoglobinemia. *J Educ Teach Emerg Med* 2022. [PMID 37465138](https://pubmed.ncbi.nlm.nih.gov/37465138/)  
   <sub>⌕ `methemoglobinemia review clinical diagnosis management`</sub>
5. **Ash-Bernal R.** Acquired methemoglobinemia: a retrospective series of 138 cases at 2 teaching hospitals. *Medicine (Baltimore)* 2004. [PMID 15342970](https://pubmed.ncbi.nlm.nih.gov/15342970/)  
   <sub>⌕ `acquired methemoglobinemia case series retrospective etiology`</sub>
6. **Friedman N.** Acquired methemoglobinemia presenting to the pediatric emergency department: a clinical challenge. *CJEM* 2020. [PMID 32396060](https://pubmed.ncbi.nlm.nih.gov/32396060/)  
   <sub>⌕ `methemoglobinemia emergency department presentation`</sub>
7. **Özbilen M.** Ferric Carboxymaltose-mediated Methemoglobinemia. *Curr Drug Saf* 2024. [PMID 36779493](https://pubmed.ncbi.nlm.nih.gov/36779493/)  
   <sub>⌕ `methemoglobinemia emergency department presentation`</sub>
8. **Sather JE.** Toxins. *Anesthesiol Clin* 2006. [PMID 17240611](https://pubmed.ncbi.nlm.nih.gov/17240611/)  
   <sub>⌕ `methemoglobinemia emergency department presentation`</sub>

## B. 산화환원 생화학 · 헤모글로빈 자가산화

> KAUTO(자가산화 2.8%/일)와 VMXB5/KMB5(내인성 환원 ~15%/h)의 근거. 이 두 숫자만으로 정상 MetHb 0.8%가 결정된다.

9. **Perrone S.** Oxidative injury in neonatal erythrocytes. *J Matern Fetal Neonatal Med* 2012. [PMID 23025782](https://pubmed.ncbi.nlm.nih.gov/23025782/)  
   <sub>⌕ `hemoglobin autoxidation superoxide mechanism erythrocyte`</sub>
10. **Watkins JA.** Autoxidation reactions of hemoglobin A free from other red cell components: a minimal mechanism. *Biochem Biophys Res Commun* 1985. [PMID 2998382](https://pubmed.ncbi.nlm.nih.gov/2998382/)  
   <sub>⌕ `hemoglobin autoxidation superoxide mechanism erythrocyte`</sub>
11. **Johnson RM.** Hemoglobin autoxidation and regulation of endogenous H2O2 levels in erythrocytes. *Free Radic Biol Med* 2005. [PMID 16274876](https://pubmed.ncbi.nlm.nih.gov/16274876/)  
   <sub>⌕ `hemoglobin autoxidation superoxide mechanism erythrocyte`</sub>
12. **Clemens MR.** Lipid peroxidation in erythrocytes. *Chem Phys Lipids* 1987. [PMID 3319229](https://pubmed.ncbi.nlm.nih.gov/3319229/)  
   <sub>⌕ `hemoglobin autoxidation superoxide mechanism erythrocyte`</sub>
13. **Tomoda A.** Kinetic studies on methemoglobin reduction by human red cell NADH cytochrome b5 reductase. *J Biol Chem* 1979. [PMID 429336](https://pubmed.ncbi.nlm.nih.gov/429336/)  
   <sub>⌕ `methemoglobin reduction cytochrome b5 reductase erythrocyte kinetics`</sub>
14. **Kuma F.** Properties of methemoglobin reductase and kinetic study of methemoglobin reduction. *J Biol Chem* 1981. [PMID 7240153](https://pubmed.ncbi.nlm.nih.gov/7240153/)  
   <sub>⌕ `methemoglobin reduction cytochrome b5 reductase erythrocyte kinetics`</sub>
15. **Sannes LJ.** Effects of hemolysate concentration, ionic strength and cytochrome b5 concentration on the rate of methemoglobin reduction in hemolysates of human erythrocytes. *Biochim Biophys Acta* 1978. [PMID 31928](https://pubmed.ncbi.nlm.nih.gov/31928/)  
   <sub>⌕ `methemoglobin reduction cytochrome b5 reductase erythrocyte kinetics`</sub>
16. **Yubisui T.** Reduction of methemoglobin through flavin at the physiological concentration by NADPH-flavin reductase of human erythrocytes. *J Biochem* 1980. [PMID 7400118](https://pubmed.ncbi.nlm.nih.gov/7400118/)  
   <sub>⌕ `methemoglobin reduction cytochrome b5 reductase erythrocyte kinetics`</sub>
17. **Elahian F.** Human cytochrome b5 reductase: structure, function, and potential applications. *Crit Rev Biotechnol* 2014. [PMID 23113554](https://pubmed.ncbi.nlm.nih.gov/23113554/)  
   <sub>⌕ `NADH cytochrome b5 reductase purification mechanism methemoglobin`</sub>
18. **Yubisui T.** Purification and properties of soluble NADH-cytochrome b5 reductase of rabbit erythrocytes. *J Biochem* 1982. [PMID 7096301](https://pubmed.ncbi.nlm.nih.gov/7096301/)  
   <sub>⌕ `NADH cytochrome b5 reductase purification mechanism methemoglobin`</sub>

## C. 산소해리곡선 · 알로스테릭 효과

> 이 모델의 헤드라인 결과가 서 있는 곳. ALPHAM(좌측 이동)과 NEFF(협동성 소실)의 근거이며, 동시에 이 모델에서 문헌적 근거가 가장 약한 파라미터이기도 하다 — 그 사실을 감추지 않고 민감도 분석으로 제시한다.

19. **Böning D.** The oxygen dissociation curve of blood in COVID-19. *Am J Physiol Lung Cell Mol Physiol* 2021. [PMID 33978488](https://pubmed.ncbi.nlm.nih.gov/33978488/)  
   <sub>⌕ `methemoglobin oxygen dissociation curve left shift affinity`</sub>
20. **Hess W.** Affinity of oxygen for hemoglobin--its significance under physiological and pathological conditions. *Anaesthesist* 1987. [PMID 3318547](https://pubmed.ncbi.nlm.nih.gov/3318547/)  
   <sub>⌕ `methemoglobin oxygen dissociation curve left shift affinity`</sub>
21. **Srikanth MS.** Topical benzocaine (Hurricaine) induced methemoglobinemia during endoscopic procedures in gastric bypass patients. *Obes Surg* 2005. [PMID 15946444](https://pubmed.ncbi.nlm.nih.gov/15946444/)  
   <sub>⌕ `methemoglobin oxygen dissociation curve left shift affinity`</sub>
22. **Aimo G.** Preliminary study on P50 and on other parameters correlated with hemoglobin function in normal subjects chosen at random. *Minerva Med* 1981. [PMID 7290468](https://pubmed.ncbi.nlm.nih.gov/7290468/)  
   <sub>⌕ `methemoglobin oxygen dissociation curve left shift affinity`</sub>
23. **Kiger L.** Recombinant Phe(beta)63hemoglobin shows rapid oxidation of the beta chains and low-affinity, non-cooperative oxygen binding to the alpha subunits. *Eur J Biochem* 1997. [PMID 9030761](https://pubmed.ncbi.nlm.nih.gov/9030761/)  
   <sub>⌕ `valency hybrid hemoglobin ferric subunit cooperativity oxygen affinity`</sub>
24. **Versmold HT.** Oxygen transport in congenital heart disease: influence of fetal hemoglobin, red cell pH, and 2,3-diphosphoglycerate. *Pediatr Res* 1976. [PMID 5699](https://pubmed.ncbi.nlm.nih.gov/5699/)  
   <sub>⌕ `2,3-diphosphoglycerate P50 oxygen affinity erythrocyte adaptation`</sub>
25. **Bureau MA.** Maternal cigarette smoking and fetal oxygen transport: a study of P50, 2,3-diphosphoglycerate, total hemoglobin, hematocrit, and type F hemoglobin in fetal blood. *Pediatrics* 1983. [PMID 6191270](https://pubmed.ncbi.nlm.nih.gov/6191270/)  
   <sub>⌕ `2,3-diphosphoglycerate P50 oxygen affinity erythrocyte adaptation`</sub>
26. **Bursaux E.** Oxygen transport in children on maintenance haemodialysis. *Clin Sci Mol Med* 1978. [PMID 620497](https://pubmed.ncbi.nlm.nih.gov/620497/)  
   <sub>⌕ `2,3-diphosphoglycerate P50 oxygen affinity erythrocyte adaptation`</sub>
27. **Paglia DE.** Atypical hematological response to combined calorie restriction and chronic hypoxia in Biosphere 2 crew: a possible link to latent features of hibernation capacity. *Habitation (Elmsford)* 2005. [PMID 15742531](https://pubmed.ncbi.nlm.nih.gov/15742531/)  
   <sub>⌕ `2,3-diphosphoglycerate P50 oxygen affinity erythrocyte adaptation`</sub>
28. **Lavrinenko IA.** Mathematical models describing oxygen binding by hemoglobin. *Biophys Rev* 2023. [PMID 37974982](https://pubmed.ncbi.nlm.nih.gov/37974982/)  
   <sub>⌕ `Hill coefficient hemoglobin cooperativity ligand binding`</sub>
29. **Lavrinenko IA.** A New Model of Hemoglobin Oxygenation. *Entropy (Basel)* 2022. [PMID 36141103](https://pubmed.ncbi.nlm.nih.gov/36141103/)  
   <sub>⌕ `Hill coefficient hemoglobin cooperativity ligand binding`</sub>
30. **Holt JM.** The Hill coefficient: inadequate resolution of cooperativity in human hemoglobin. *Methods Enzymol* 2009. [PMID 19289207](https://pubmed.ncbi.nlm.nih.gov/19289207/)  
   <sub>⌕ `Hill coefficient hemoglobin cooperativity ligand binding`</sub>
31. **Edelstein SJ.** Contributions of individual molecular species to the Hill coefficient for ligand binding by an oligomeric protein. *J Mol Biol* 1997. [PMID 9096203](https://pubmed.ncbi.nlm.nih.gov/9096203/)  
   <sub>⌕ `Hill coefficient hemoglobin cooperativity ligand binding`</sub>

## D. 맥박산소측정 · 동시산소측정

> E660M·E940M을 맞춘 대상 데이터. 85% 바닥은 파라미터가 아니라 R→1에서의 보정 직선 값으로 계산되어 나온다.

32. **Toffaletti J.** Misconceptions in reporting oxygen saturation. *Anesth Analg* 2007. [PMID 18048899](https://pubmed.ncbi.nlm.nih.gov/18048899/)  
   <sub>⌕ `pulse oximetry dyshemoglobin saturation gap methemoglobin`</sub>
33. **Feiner JR.** Accuracy of methemoglobin detection by pulse CO-oximetry during hypoxia. *Anesth Analg* 2010. [PMID 20007731](https://pubmed.ncbi.nlm.nih.gov/20007731/)  
   <sub>⌕ `co-oximetry multiwavelength methemoglobin measurement`</sub>
34. **Feiner JR.** Improved accuracy of methemoglobin detection by pulse CO-oximetry during hypoxia. *Anesth Analg* 2010. [PMID 20841412](https://pubmed.ncbi.nlm.nih.gov/20841412/)  
   <sub>⌕ `co-oximetry multiwavelength methemoglobin measurement`</sub>
35. **Kuleš J.** Co-oximetry in clinically healthy dogs and effects of time of post sampling on measurements. *J Small Anim Pract* 2011. [PMID 21981712](https://pubmed.ncbi.nlm.nih.gov/21981712/)  
   <sub>⌕ `co-oximetry multiwavelength methemoglobin measurement`</sub>
36. **Suzaki H.** Noninvasive measurement of total hemoglobin and hemoglobin derivatives using multiwavelength pulse spectrophotometry -In vitro study with a mock circulatory system. *Conf Proc IEEE Eng Med Biol Soc* 2006. [PMID 17946423](https://pubmed.ncbi.nlm.nih.gov/17946423/)  
   <sub>⌕ `co-oximetry multiwavelength methemoglobin measurement`</sub>
37. **Michaelis G.** Effect of dyshemoglobinemia (methemoglobinemia and carboxyhemoglobinemia) on accuracy of measurement in pulse oximetry in operations of long duration. *Anasth Intensivther Notfallmed* 1988. [PMID 3394901](https://pubmed.ncbi.nlm.nih.gov/3394901/)  
   <sub>⌕ `pulse oximetry limitations dyshemoglobinemia carboxyhemoglobin accuracy`</sub>
38. **Thrush D.** Accuracy of pulse oximetry during hypoxemia. *South Med J* 1994. [PMID 8153783](https://pubmed.ncbi.nlm.nih.gov/8153783/)  
   <sub>⌕ `methemoglobin 85% pulse oximeter plateau saturation`</sub>
39. **Rausch-Madison S.** Methodologic problems encountered with cooximetry in methemoglobinemia. *Am J Med Sci* 1997. [PMID 9298047](https://pubmed.ncbi.nlm.nih.gov/9298047/)  
   <sub>⌕ `arterial blood gas calculated saturation error methemoglobinemia`</sub>

## E. 벤조카인 · 국소마취제

> 급성 노출 시나리오의 보정 근거(스프레이 1회 → 20–40분 내 15–40%).

40. **Cooper HA.** Methemoglobinemia caused by benzocaine topical spray. *South Med J* 1997. [PMID 9305310](https://pubmed.ncbi.nlm.nih.gov/9305310/)  
   <sub>⌕ `benzocaine induced methemoglobinemia topical spray`</sub>
41. **Linares LA.** Methemoglobinemia induced by topical anesthetic (benzocaine). *Radiother Oncol* 1990. [PMID 2217871](https://pubmed.ncbi.nlm.nih.gov/2217871/)  
   <sub>⌕ `benzocaine induced methemoglobinemia topical spray`</sub>
42. **Harvey JW.** Benzocaine-induced methemoglobinemia in dogs. *J Am Vet Med Assoc* 1979. [PMID 511741](https://pubmed.ncbi.nlm.nih.gov/511741/)  
   <sub>⌕ `benzocaine induced methemoglobinemia topical spray`</sub>
43. **Sewell CR.** A Case Report of Benzocaine-Induced Methemoglobinemia. *J Pharm Pract* 2018. [PMID 28803519](https://pubmed.ncbi.nlm.nih.gov/28803519/)  
   <sub>⌕ `benzocaine induced methemoglobinemia topical spray`</sub>
44. **Gunes H.** Prilocaine-induced Methemoglobinemia. *J Coll Physicians Surg Pak* 2017. [PMID 28903852](https://pubmed.ncbi.nlm.nih.gov/28903852/)  
   <sub>⌕ `prilocaine methemoglobinemia local anesthetic`</sub>
45. **Centers for Disease Control and Prevention (CDC).** Prilocaine-induced methemoglobinemia--Wisconsin, 1993. *MMWR Morb Mortal Wkly Rep* 1994. [PMID 8072476](https://pubmed.ncbi.nlm.nih.gov/8072476/)  
   <sub>⌕ `prilocaine methemoglobinemia local anesthetic`</sub>
46. **Rehman HU.** Methemoglobinemia. *West J Med* 2001. [PMID 11527852](https://pubmed.ncbi.nlm.nih.gov/11527852/)  
   <sub>⌕ `prilocaine methemoglobinemia local anesthetic`</sub>
47. **Yıldız A.** MEP-43 Methemoglobinemia After Central Venous Catheterization Due to Local Anesthesia with Prilocaine: A Case Report. *Turk Gogus Kalp Damar Cerrahisi Derg* 2024. [PMID 40322132](https://pubmed.ncbi.nlm.nih.gov/40322132/)  
   <sub>⌕ `prilocaine methemoglobinemia local anesthetic`</sub>

## F. 답손 · 아릴하이드록실아민 공동산화

> 모델에서 가장 중요한 구조적 주장의 근거: 아릴하이드록실아민은 촉매이고 그 조효소는 글루타티온이다. 시메티딘 실험(생성 항 공격)도 여기에 있다.

48. **Veggi LM.** Dapsone induces oxidative stress and impairs antioxidant defenses in rat liver. *Life Sci* 2008. [PMID 18602405](https://pubmed.ncbi.nlm.nih.gov/18602405/)  
   <sub>⌕ `dapsone induced methemoglobinemia hemolysis glutathione`</sub>
49. **MacDonald RD.** Acute dapsone intoxication: a pediatric case report. *Pediatr Emerg Care* 1997. [PMID 9127424](https://pubmed.ncbi.nlm.nih.gov/9127424/)  
   <sub>⌕ `dapsone overdose methemoglobinemia prolonged methylene blue`</sub>
50. **Zia S.** Intentional dapsone overdose: a rare cause of acquired methaemoglobinaemia in the emergency department. *BMJ Case Rep* 2022. [PMID 36127031](https://pubmed.ncbi.nlm.nih.gov/36127031/)  
   <sub>⌕ `dapsone overdose methemoglobinemia prolonged methylene blue`</sub>
51. **Jollow DJ.** Dapsone-induced hemolytic anemia. *Drug Metab Rev* 1995. [PMID 7641572](https://pubmed.ncbi.nlm.nih.gov/7641572/)  
   <sub>⌕ `dapsone hydroxylamine hemolytic anemia mechanism red cell`</sub>
52. **Paika S.** Molecular Mechanisms of Drug-Induced Hemolysis in G6PD Deficiency: Mechanistic Insights. *Oxid Med Cell Longev* 2025. [PMID 40799291](https://pubmed.ncbi.nlm.nih.gov/40799291/)  
   <sub>⌕ `dapsone hydroxylamine hemolytic anemia mechanism red cell`</sub>
53. **Haas M.** Stimulation of K-C1 cotransport in rat red cells by a hemolytic anemia-producing metabolite of dapsone. *Am J Physiol* 1989. [PMID 2919657](https://pubmed.ncbi.nlm.nih.gov/2919657/)  
   <sub>⌕ `dapsone hydroxylamine hemolytic anemia mechanism red cell`</sub>
54. **Bian Y.** Dapsone Hydroxylamine, an Active Metabolite of Dapsone, Can Promote the Procoagulant Activity of Red Blood Cells and Thrombosis. *Toxicol Sci* 2019. [PMID 31428780](https://pubmed.ncbi.nlm.nih.gov/31428780/)  
   <sub>⌕ `dapsone hydroxylamine hemolytic anemia mechanism red cell`</sub>
55. **Zuidema J.** Clinical pharmacokinetics of dapsone. *Clin Pharmacokinet* 1986. [PMID 3530584](https://pubmed.ncbi.nlm.nih.gov/3530584/)  
   <sub>⌕ `dapsone pharmacokinetics enterohepatic half-life human`</sub>
56. **Di Girolamo F.** Mass spectrometric identification of hemoglobin modifications induced by nitrosobenzene. *Ecotoxicol Environ Saf* 2009. [PMID 18973939](https://pubmed.ncbi.nlm.nih.gov/18973939/)  
   <sub>⌕ `aniline nitrobenzene metabolite methemoglobin formation kinetics`</sub>
57. **Alfirevic A.** Slow acetylator phenotype and genotype in HIV-positive patients with sulphamethoxazole hypersensitivity. *Br J Clin Pharmacol* 2003. [PMID 12580987](https://pubmed.ncbi.nlm.nih.gov/12580987/)  
   <sub>⌕ `NAT2 acetylator status dapsone metabolism`</sub>
58. **Aynacioglu AS.** N-acetyltransferase polymorphism in patients with Behçet's disease. *Eur J Clin Pharmacol* 2001. [PMID 11791896](https://pubmed.ncbi.nlm.nih.gov/11791896/)  
   <sub>⌕ `NAT2 acetylator status dapsone metabolism`</sub>
59. **Brocvielle H.** N-acetyltransferase 2 acetylation polymorphism: prevalence of slow acetylators does not differ between atopic dermatitis patients and healthy subjects. *Skin Pharmacol Appl Skin Physiol* 2003. [PMID 14528063](https://pubmed.ncbi.nlm.nih.gov/14528063/)  
   <sub>⌕ `NAT2 acetylator status dapsone metabolism`</sub>
60. **Furet Y.** Clinical relevance of N-acetyltransferase type 2 (NAT2) genetic polymorphism. *Therapie* 2002. [PMID 12611196](https://pubmed.ncbi.nlm.nih.gov/12611196/)  
   <sub>⌕ `NAT2 acetylator status dapsone metabolism`</sub>

## G. 아질산염 · 질산염

> KNIT2(자가촉매 항)의 근거. 지연 후 급상승이라는 형태는 이 항 하나에서 나온다.

61. **Kosaka H.** Mechanism of autocatalytic oxidation of oxyhemoglobin by nitrite. *Environ Health Perspect* 1987. [PMID 2822381](https://pubmed.ncbi.nlm.nih.gov/2822381/)  
   <sub>⌕ `nitrite oxidation oxyhemoglobin autocatalytic kinetics`</sub>
62. **Kosaka H.** Mechanism of autocatalytic oxidation of oxyhemoglobin by nitrite. *Biomed Biochim Acta* 1983. [PMID 6326765](https://pubmed.ncbi.nlm.nih.gov/6326765/)  
   <sub>⌕ `nitrite oxidation oxyhemoglobin autocatalytic kinetics`</sub>
63. **Doyle MP.** Autocatalytic oxidation of hemoglobin induced by nitrite: activation and chemical inhibition. *J Free Radic Biol Med* 1985. [PMID 3836241](https://pubmed.ncbi.nlm.nih.gov/3836241/)  
   <sub>⌕ `nitrite oxidation oxyhemoglobin autocatalytic kinetics`</sub>
64. **Kosaka H.** Mechanism of autocatalytic oxidation of oxyhemoglobin by nitrite. An intermediate detected by electron spin resonance. *Biochim Biophys Acta* 1982. [PMID 6282334](https://pubmed.ncbi.nlm.nih.gov/6282334/)  
   <sub>⌕ `nitrite oxidation oxyhemoglobin autocatalytic kinetics`</sub>
65. **Dean DE.** Fatal methemoglobinemia in three suicidal sodium nitrite poisonings. *J Forensic Sci* 2021. [PMID 33598944](https://pubmed.ncbi.nlm.nih.gov/33598944/)  
   <sub>⌕ `sodium nitrite poisoning methemoglobinemia fatal`</sub>
66. **Neth MR.** Fatal Sodium Nitrite Poisoning: Key Considerations for Prehospital Providers. *Prehosp Emerg Care* 2021. [PMID 33074043](https://pubmed.ncbi.nlm.nih.gov/33074043/)  
   <sub>⌕ `sodium nitrite poisoning methemoglobinemia fatal`</sub>
67. **Mun SH.** Two cases of fatal methemoglobinemia caused by self-poisoning with sodium nitrite: A case report. *Medicine (Baltimore)* 2022. [PMID 35363170](https://pubmed.ncbi.nlm.nih.gov/35363170/)  
   <sub>⌕ `sodium nitrite poisoning methemoglobinemia fatal`</sub>
68. **Cvetković D.** Sodium nitrite food poisoning in one family. *Forensic Sci Med Pathol* 2019. [PMID 30293223](https://pubmed.ncbi.nlm.nih.gov/30293223/)  
   <sub>⌕ `sodium nitrite poisoning methemoglobinemia fatal`</sub>
69. **Knobeloch L.** Blue babies and nitrate-contaminated well water. *Environ Health Perspect* 2000. [PMID 10903623](https://pubmed.ncbi.nlm.nih.gov/10903623/)  
   <sub>⌕ `infant methemoglobinemia well water nitrate blue baby`</sub>
70. **Fossen Johnson S.** Methemoglobinemia: Infants at risk. *Curr Probl Pediatr Adolesc Health Care* 2019. [PMID 30956100](https://pubmed.ncbi.nlm.nih.gov/30956100/)  
   <sub>⌕ `infant methemoglobinemia well water nitrate blue baby`</sub>
71. **Kross BC.** Methemoglobinemia: nitrate toxicity in rural America. *Am Fam Physician* 1992. [PMID 1621630](https://pubmed.ncbi.nlm.nih.gov/1621630/)  
   <sub>⌕ `infant methemoglobinemia well water nitrate blue baby`</sub>
72. **Miller LW.** Methemoglobinemia associated with well water. *JAMA* 1971. [PMID 5108507](https://pubmed.ncbi.nlm.nih.gov/5108507/)  
   <sub>⌕ `infant methemoglobinemia well water nitrate blue baby`</sub>

## H. 메틸렌블루 약동학 · 약력학

> MB의 PK와, 갈림길(BR1/BR2) 구조 및 고용량 역설의 근거.

73. **Clifton J 2nd.** Methylene blue. *Am J Ther* 2003. [PMID 12845393](https://pubmed.ncbi.nlm.nih.gov/12845393/)  
   <sub>⌕ `methylene blue pharmacokinetics human intravenous`</sub>
74. **Bateman RM.** 36th International Symposium on Intensive Care and Emergency Medicine : Brussels, Belgium. 15-18 March 2016. *Crit Care* 2016. [PMID 27885969](https://pubmed.ncbi.nlm.nih.gov/27885969/)  
   <sub>⌕ `methylene blue pharmacokinetics human intravenous`</sub>
75. **Aeschlimann C.** Comparative pharmacokinetics of oral and intravenous ifosfamide/mesna/methylene blue therapy. *Drug Metab Dispos* 1998. [PMID 9733667](https://pubmed.ncbi.nlm.nih.gov/9733667/)  
   <sub>⌕ `methylene blue pharmacokinetics human intravenous`</sub>
76. **Peter C.** Pharmacokinetics and organ distribution of intravenous and oral methylene blue. *Eur J Clin Pharmacol* 2000. [PMID 10952480](https://pubmed.ncbi.nlm.nih.gov/10952480/)  
   <sub>⌕ `methylene blue pharmacokinetics human intravenous`</sub>
77. **Haouzi P.** Revisiting the physiological effects of methylene blue as a treatment of cyanide intoxication. *Clin Toxicol (Phila)* 2018. [PMID 29451035](https://pubmed.ncbi.nlm.nih.gov/29451035/)  
   <sub>⌕ `methylene blue methemoglobinemia treatment dose efficacy`</sub>
78. **Rodriguez P.** Methylene blue treatment delays progression of perfusion-diffusion mismatch to infarct in permanent ischemic stroke. *Brain Res* 2014. [PMID 25218555](https://pubmed.ncbi.nlm.nih.gov/25218555/)  
   <sub>⌕ `methylene blue methemoglobinemia treatment dose efficacy`</sub>
79. **Rothenberg R.** Effectiveness and tolerability of methylthioninium chloride (methylene blue) for the treatment of methemoglobinemia: twenty-four years of experience at a single poison center. *Clin Toxicol (Phila)* 2025. [PMID 40062661](https://pubmed.ncbi.nlm.nih.gov/40062661/)  
   <sub>⌕ `methylene blue methemoglobinemia treatment dose efficacy`</sub>
80. **Bradberry SM.** Occupational methaemoglobinaemia. Mechanisms of production, features, diagnosis and management including the use of methylene blue. *Toxicol Rev* 2003. [PMID 14579544](https://pubmed.ncbi.nlm.nih.gov/14579544/)  
   <sub>⌕ `methylene blue methemoglobinemia treatment dose efficacy`</sub>
81. **Lo YH.** High dose vitamin C induced methemoglobinemia and hemolytic anemia in glucose-6-phosphate dehydrogenase deficiency. *Am J Emerg Med* 2020. [PMID 32561141](https://pubmed.ncbi.nlm.nih.gov/32561141/)  
   <sub>⌕ `methylene blue paradoxical methemoglobinemia high dose`</sub>
82. **Blank O.** Interactions of the antimalarial drug methylene blue with methemoglobin and heme targets in Plasmodium falciparum: a physico-biochemical study. *Antioxid Redox Signal* 2012. [PMID 22256987](https://pubmed.ncbi.nlm.nih.gov/22256987/)  
   <sub>⌕ `leucomethylene blue NADPH methemoglobin reductase mechanism`</sub>

## I. G6PD 결핍 · 용혈

> NADPH 공유 상한(CAP)과 하인츠소체 이중 구동의 근거. 모델은 '금기'를 규칙으로 넣지 않고 이 절의 생화학에서 유도한다.

83. **Youngster I.** Medications and glucose-6-phosphate dehydrogenase deficiency: an evidence-based review. *Drug Saf* 2010. [PMID 20701405](https://pubmed.ncbi.nlm.nih.gov/20701405/)  
   <sub>⌕ `glucose-6-phosphate dehydrogenase deficiency methylene blue hemolysis`</sub>
84. **Mak GK.** Glucose-6-Phosphate Dehydrogenase Deficiency. ** 2026. [PMID 29262208](https://pubmed.ncbi.nlm.nih.gov/29262208/)  
   <sub>⌕ `glucose-6-phosphate dehydrogenase deficiency methylene blue hemolysis`</sub>
85. **Belfield KD.** Review and drug therapy implications of glucose-6-phosphate dehydrogenase deficiency. *Am J Health Syst Pharm* 2018. [PMID 29305344](https://pubmed.ncbi.nlm.nih.gov/29305344/)  
   <sub>⌕ `glucose-6-phosphate dehydrogenase deficiency methylene blue hemolysis`</sub>
86. **Pannu AK.** Naphthalene Toxicity in Clinical Practice. *Curr Drug Metab* 2020. [PMID 31755382](https://pubmed.ncbi.nlm.nih.gov/31755382/)  
   <sub>⌕ `glucose-6-phosphate dehydrogenase deficiency methylene blue hemolysis`</sub>
87. **Nannelli C.** Genetic variants causing G6PD deficiency: Clinical and biochemical data support new WHO classification. *Br J Haematol* 2023. [PMID 37415281](https://pubmed.ncbi.nlm.nih.gov/37415281/)  
   <sub>⌕ `G6PD deficiency variants classification WHO activity`</sub>
88. **Geck RC.** Functional interpretation, cataloging, and analysis of 1,341 glucose-6-phosphate dehydrogenase variants. *Am J Hum Genet* 2023. [PMID 36681081](https://pubmed.ncbi.nlm.nih.gov/36681081/)  
   <sub>⌕ `G6PD deficiency variants classification WHO activity`</sub>
89. **Watchko JF.** Revised World Health Organization (WHO) classification of G6PD gene variants: Relevance to neonatal hyperbilirubinemia. *Semin Fetal Neonatal Med* 2025. [PMID 40023662](https://pubmed.ncbi.nlm.nih.gov/40023662/)  
   <sub>⌕ `G6PD deficiency variants classification WHO activity`</sub>
90. **Pfeffer DA.** Genetic Variants of Glucose-6-Phosphate Dehydrogenase and Their Associated Enzyme Activity: A Systematic Review and Meta-Analysis. *Pathogens* 2022. [PMID 36145477](https://pubmed.ncbi.nlm.nih.gov/36145477/)  
   <sub>⌕ `G6PD deficiency variants classification WHO activity`</sub>
91. **Christopher MM.** Erythrocyte pathology and mechanisms of Heinz body-mediated hemolysis in cats. *Vet Pathol* 1990. [PMID 2238384](https://pubmed.ncbi.nlm.nih.gov/2238384/)  
   <sub>⌕ `Heinz body formation oxidative hemolysis glutathione erythrocyte`</sub>
92. **Gordon-Smith EC.** Oxidative haemolysis and Heinz body haemolytic anaemia. *Br J Haematol* 1974. [PMID 4367628](https://pubmed.ncbi.nlm.nih.gov/4367628/)  
   <sub>⌕ `Heinz body formation oxidative hemolysis glutathione erythrocyte`</sub>
93. **Hopkins J.** The effects of drugs on erythrocytes in vitro: Heinz body formation, glutathione peroxidase inhibition and changes in mechanical fragility. *Br J Clin Pharmacol* 1974. [PMID 22454946](https://pubmed.ncbi.nlm.nih.gov/22454946/)  
   <sub>⌕ `Heinz body formation oxidative hemolysis glutathione erythrocyte`</sub>
94. **Furdak P.** Effect of Garlic Extract on the Erythrocyte as a Simple Model Cell. *Int J Mol Sci* 2024. [PMID 38791153](https://pubmed.ncbi.nlm.nih.gov/38791153/)  
   <sub>⌕ `Heinz body formation oxidative hemolysis glutathione erythrocyte`</sub>

## J. 선천성 메트헤모글로빈혈증 · HbM

> 만성 대 급성 비대칭(2,3-BPG 우측 이동이 좌측 이동을 상쇄)과 HbM 비반응성의 근거.

95. **Kedar PS.** Novel mutation (R192C) in CYB5R3 gene causing NADH-cytochrome b5 reductase deficiency in eight Indian patients associated with autosomal recessive congenital methemoglobinemia type-I. *Hematology* 2018. [PMID 29482478](https://pubmed.ncbi.nlm.nih.gov/29482478/)  
   <sub>⌕ `congenital methemoglobinemia cytochrome b5 reductase deficiency type`</sub>
96. **Nicolas-Jilwan M.** Recessive congenital methemoglobinemia type II: Hypoplastic basal ganglia in two siblings with a novel mutation of the cytochrome b5 reductase gene. *Neuroradiol J* 2019. [PMID 30614390](https://pubmed.ncbi.nlm.nih.gov/30614390/)  
   <sub>⌕ `congenital methemoglobinemia cytochrome b5 reductase deficiency type`</sub>
97. **Kedar PS.** Congenital methemoglobinemia due to NADH-methemoglobin reductase deficiency in three Indian families. *Haematologia (Budap)* 2002. [PMID 12803131](https://pubmed.ncbi.nlm.nih.gov/12803131/)  
   <sub>⌕ `congenital methemoglobinemia cytochrome b5 reductase deficiency type`</sub>
98. **Kaplan JC.** Prenatal diagnosis of generalized cytochrome b5 reductase deficiency (congenital methemoglobinemia with mental retardation, type II) (author's transl). *Ann Med Interne (Paris)* 1981. [PMID 7235451](https://pubmed.ncbi.nlm.nih.gov/7235451/)  
   <sub>⌕ `congenital methemoglobinemia cytochrome b5 reductase deficiency type`</sub>
99. **Kim DS.** The First Korean Family with Hemoglobin-M Milwaukee-2 Leading to Hereditary Methemoglobinemia. *Yonsei Med J* 2020. [PMID 33251782](https://pubmed.ncbi.nlm.nih.gov/33251782/)  
   <sub>⌕ `hemoglobin M variant methemoglobinemia`</sub>
100. **Schnedl WJ.** Hereditary methemoglobinemia caused by Hb M-Hyde Park (Hb M-Akita) (HBB:c.277C > T; p.His93Tyr). *Wien Klin Wochenschr* 2019. [PMID 31267164](https://pubmed.ncbi.nlm.nih.gov/31267164/)  
   <sub>⌕ `hemoglobin M variant methemoglobinemia`</sub>
101. **Spears F.** Hemoglobin M variant and congenital methemoglobinemia: methylene blue will not be effective in the presence of hemoglobin M. *Can J Anaesth* 2008. [PMID 18566207](https://pubmed.ncbi.nlm.nih.gov/18566207/)  
   <sub>⌕ `hemoglobin M variant methemoglobinemia`</sub>
102. **Rangan A.** Interpreting sulfhemoglobin and methemoglobin in patients with cyanosis: An overview of patients with M-hemoglobin variants. *Int J Lab Hematol* 2021. [PMID 34092029](https://pubmed.ncbi.nlm.nih.gov/34092029/)  
   <sub>⌕ `hemoglobin M variant methemoglobinemia`</sub>
103. **Yamaji F.** A new mutation of congenital methemoglobinemia exacerbated after methylene blue treatment. *Acute Med Surg* 2018. [PMID 29657736](https://pubmed.ncbi.nlm.nih.gov/29657736/)  
   <sub>⌕ `riboflavin treatment congenital methemoglobinemia`</sub>
104. **Neven J.** Recessive congenital methemoglobinemia: a systematic review of reported cases. *Orphanet J Rare Dis* 2026. [PMID 41639689](https://pubmed.ncbi.nlm.nih.gov/41639689/)  
   <sub>⌕ `riboflavin treatment congenital methemoglobinemia`</sub>
105. **Hirano M.** Congenital methaemoglobinaemia due to NADH methaemoglobin reductase deficiency: successful treatment with oral riboflavin. *Br J Haematol* 1981. [PMID 6893937](https://pubmed.ncbi.nlm.nih.gov/6893937/)  
   <sub>⌕ `riboflavin treatment congenital methemoglobinemia`</sub>
106. **Akompong T.** In vitro activity of riboflavin against the human malaria parasite Plasmodium falciparum. *Antimicrob Agents Chemother* 2000. [PMID 10602728](https://pubmed.ncbi.nlm.nih.gov/10602728/)  
   <sub>⌕ `riboflavin treatment congenital methemoglobinemia`</sub>

## K. 대체 · 구제 치료

> 세 번째 항(운반)을 공격하는 치료들. 보에레마 산술의 원 출처 포함.

107. **Menakuru SR.** Phenazopyridine-Induced Methemoglobinemia in a Jehovah's Witness Treated with High-Dose Ascorbic Acid Due to Methylene Blue Contradictions: A Case Report and Review of the Literature. *Hematol Rep* 2023. [PMID 37367083](https://pubmed.ncbi.nlm.nih.gov/37367083/)  
   <sub>⌕ `ascorbic acid treatment methemoglobinemia high dose`</sub>
108. **Lee KW.** High-dose vitamin C as treatment of methemoglobinemia. *Am J Emerg Med* 2014. [PMID 24972962](https://pubmed.ncbi.nlm.nih.gov/24972962/)  
   <sub>⌕ `ascorbic acid treatment methemoglobinemia high dose`</sub>
109. **Visclosky T.** Primaquine overdose in a toddler. *Am J Emerg Med* 2021. [PMID 33279327](https://pubmed.ncbi.nlm.nih.gov/33279327/)  
   <sub>⌕ `ascorbic acid treatment methemoglobinemia high dose`</sub>
110. **Chen HJ.** New-Onset Methemoglobinemia After Receiving Ozone Autohemotherapy and High-Dose Vitamin C Injection. *Transfus Med Rev* 2022. [PMID 34815107](https://pubmed.ncbi.nlm.nih.gov/34815107/)  
   <sub>⌕ `ascorbic acid treatment methemoglobinemia high dose`</sub>
111. **Lavonas EJ.** 2023 American Heart Association Focused Update on the Management of Patients With Cardiac Arrest or Life-Threatening Toxicity Due to Poisoning: An Update to the American Heart Association Guidelines for Cardiopulmonary Resuscitation and Emergency Cardiovascular Care. *Circulation* 2023. [PMID 37721023](https://pubmed.ncbi.nlm.nih.gov/37721023/)  
   <sub>⌕ `hyperbaric oxygen methemoglobinemia treatment`</sub>
112. **Altintop I.** Methemoglobinemia treated with hyperbaric oxygen therapy: A case report. *Turk J Emerg Med* 2018. [PMID 30533564](https://pubmed.ncbi.nlm.nih.gov/30533564/)  
   <sub>⌕ `hyperbaric oxygen methemoglobinemia treatment`</sub>
113. **Malik MJ.** Methemoglobinemia: A Case Report. *Cureus* 2023. [PMID 38021620](https://pubmed.ncbi.nlm.nih.gov/38021620/)  
   <sub>⌕ `hyperbaric oxygen methemoglobinemia treatment`</sub>
114. **Padhiyar V.** Nitrofurantoin-Induced Hemolytic Anemia and Methemoglobinemia in G6PD Deficiency Treated With Hyperbaric Oxygen Treatment. *Undersea Hyperb Med* 2026. [PMID 42365943](https://pubmed.ncbi.nlm.nih.gov/42365943/)  
   <sub>⌕ `hyperbaric oxygen methemoglobinemia treatment`</sub>
115. **Lanas Mendoza ES.** Exchange transfusion as rescue therapy for severe methemoglobinemia in the absence of methylene blue. *Med Clin (Barc)* 2026. [PMID 42208138](https://pubmed.ncbi.nlm.nih.gov/42208138/)  
   <sub>⌕ `exchange transfusion methemoglobinemia severe`</sub>
116. **Ranasinghe P.** Exchange transfusion can be life-saving in severe propanil poisoning: a case report. *BMC Res Notes* 2014. [PMID 25292188](https://pubmed.ncbi.nlm.nih.gov/25292188/)  
   <sub>⌕ `exchange transfusion methemoglobinemia severe`</sub>
117. **Goyal R.** When methylene blue fails: a case of methemoglobinemia induced by azoxystrobin and propiconazole managed by exchange transfusion. *Int J Emerg Med* 2025. [PMID 41184768](https://pubmed.ncbi.nlm.nih.gov/41184768/)  
   <sub>⌕ `exchange transfusion methemoglobinemia severe`</sub>
118. **Patnaik S.** Methylene blue unresponsive methemoglobinemia. *Indian J Crit Care Med* 2014. [PMID 24872659](https://pubmed.ncbi.nlm.nih.gov/24872659/)  
   <sub>⌕ `exchange transfusion methemoglobinemia severe`</sub>
119. **Jones MW.** Hyperbaric Physics. ** 2026. [PMID 28846268](https://pubmed.ncbi.nlm.nih.gov/28846268/)  
   <sub>⌕ `Boerema life without blood hyperbaric dissolved oxygen`</sub>

## L. 메틸렌블루의 다른 약리작용

> MAO-A 억제와 세로토닌 독성 — 해독제 자체의 약리학.

120. **Delport A.** Methylene Blue Analogues with Marginal Monoamine Oxidase Inhibition Retain Antidepressant-like Activity. *ACS Chem Neurosci* 2018. [PMID 29976053](https://pubmed.ncbi.nlm.nih.gov/29976053/)  
   <sub>⌕ `methylene blue monoamine oxidase inhibition serotonin toxicity`</sub>
121. **Idle JR.** Ifosfamide - History, efficacy, toxicity and encephalopathy. *Pharmacol Ther* 2023. [PMID 36842616](https://pubmed.ncbi.nlm.nih.gov/36842616/)  
   <sub>⌕ `methylene blue monoamine oxidase inhibition serotonin toxicity`</sub>
122. **Delport A.** The monoamine oxidase inhibition properties of selected structural analogues of methylene blue. *Toxicol Appl Pharmacol* 2017. [PMID 28377303](https://pubmed.ncbi.nlm.nih.gov/28377303/)  
   <sub>⌕ `methylene blue monoamine oxidase inhibition serotonin toxicity`</sub>
123. **Gillman PK.** CNS toxicity involving methylene blue: the exemplar for understanding and predicting drug interactions that precipitate serotonin toxicity. *J Psychopharmacol* 2011. [PMID 20142303](https://pubmed.ncbi.nlm.nih.gov/20142303/)  
   <sub>⌕ `methylene blue monoamine oxidase inhibition serotonin toxicity`</sub>
124. **Bužga M.** Methylene blue: a controversial diagnostic acid and medication?. *Toxicol Res (Camb)* 2022. [PMID 36337249](https://pubmed.ncbi.nlm.nih.gov/36337249/)  
   <sub>⌕ `methylene blue serotonin syndrome parathyroid surgery`</sub>
125. **Rowley M.** Methylene blue-associated serotonin syndrome: a 'green' encephalopathy after parathyroidectomy. *Neurocrit Care* 2009. [PMID 19263250](https://pubmed.ncbi.nlm.nih.gov/19263250/)  
   <sub>⌕ `methylene blue serotonin syndrome parathyroid surgery`</sub>
126. **Stanford SC.** Risk of severe serotonin toxicity following co-administration of methylene blue and serotonin reuptake inhibitors: an update on a case report of post-operative delirium. *J Psychopharmacol* 2010. [PMID 19423610](https://pubmed.ncbi.nlm.nih.gov/19423610/)  
   <sub>⌕ `methylene blue serotonin syndrome parathyroid surgery`</sub>
127. **Gillman K.** In reference to Parathyroid surgery and methylene blue: a review with guidelines for safe intraoperative use. *Laryngoscope* 2010. [PMID 19950381](https://pubmed.ncbi.nlm.nih.gov/19950381/)  
   <sub>⌕ `methylene blue serotonin syndrome parathyroid surgery`</sub>

## M. 설프헤모글로빈혈증

> MB에 반응하지 않는 감별진단. 모델은 이를 규칙이 아니라 '환원 방정식에 아예 들어가지 않음'으로 표현한다.

128. **Askew SW.** On the dysfunctional hemoglobins and cyanosis connection: practical implications for the clinical detection and differentiation of methemoglobinemia and sulfhemoglobinemia. *Biomed Opt Express* 2018. [PMID 29984098](https://pubmed.ncbi.nlm.nih.gov/29984098/)  
   <sub>⌕ `sulfhemoglobinemia diagnosis differential methemoglobinemia`</sub>
129. **Noor M.** Acquired sulfhemoglobinemia. An underreported diagnosis?. *West J Med* 1998. [PMID 9866446](https://pubmed.ncbi.nlm.nih.gov/9866446/)  
   <sub>⌕ `sulfhemoglobinemia diagnosis differential methemoglobinemia`</sub>
130. **Stepanenko T.** Sulfhemoglobin under the spotlight - Detection and characterization of SHb and HbFe(III)-SH. *Biochim Biophys Acta Mol Cell Res* 2023. [PMID 36220452](https://pubmed.ncbi.nlm.nih.gov/36220452/)  
   <sub>⌕ `sulfhemoglobinemia diagnosis differential methemoglobinemia`</sub>
131. **Baranoski GV.** On the noninvasive optical monitoring and differentiation of methemoglobinemia and sulfhemoglobinemia. *J Biomed Opt* 2012. [PMID 22975680](https://pubmed.ncbi.nlm.nih.gov/22975680/)  
   <sub>⌕ `sulfhemoglobinemia diagnosis differential methemoglobinemia`</sub>

## N. 산소운반 · 조직 저산소 · 임계 DO2

> PVCRIT(임계 모세혈관 PO₂ ~20 mmHg), 청색증 역치, 보상 기전의 근거.

132. **Bassett DR Jr.** Limiting factors for maximum oxygen uptake and determinants of endurance performance. *Med Sci Sports Exerc* 2000. [PMID 10647532](https://pubmed.ncbi.nlm.nih.gov/10647532/)  
   <sub>⌕ `critical oxygen delivery threshold anaerobic metabolism tissue`</sub>
133. **Ronco JJ.** Identification of the critical oxygen delivery for anaerobic metabolism in critically ill septic and nonseptic humans. *JAMA* 1993. [PMID 8411504](https://pubmed.ncbi.nlm.nih.gov/8411504/)  
   <sub>⌕ `critical oxygen delivery threshold anaerobic metabolism tissue`</sub>
134. **Kisaka T.** Mechanisms That Modulate Peripheral Oxygen Delivery during Exercise in Heart Failure. *Ann Am Thorac Soc* 2017. [PMID 28679061](https://pubmed.ncbi.nlm.nih.gov/28679061/)  
   <sub>⌕ `critical oxygen delivery threshold anaerobic metabolism tissue`</sub>
135. **Steltzer H.** Oxygen kinetics during liver transplantation: the relationship between delivery and consumption. *J Crit Care* 1993. [PMID 8343854](https://pubmed.ncbi.nlm.nih.gov/8343854/)  
   <sub>⌕ `critical oxygen delivery threshold anaerobic metabolism tissue`</sub>
136. **Zimmerman R.** Posttransfusion Increase of Hematocrit per se Does Not Improve Circulatory Oxygen Delivery due to Increased Blood Viscosity. *Anesth Analg* 2017. [PMID 28328758](https://pubmed.ncbi.nlm.nih.gov/28328758/)  
   <sub>⌕ `anemia oxygen delivery compensation cardiac output transfusion threshold`</sub>
137. **Bisarya R.** Serum lactate poorly predicts central venous oxygen saturation in critically ill patients: a retrospective cohort study. *J Intensive Care* 2019. [PMID 31516712](https://pubmed.ncbi.nlm.nih.gov/31516712/)  
   <sub>⌕ `mixed venous oxygen saturation critical threshold lactate`</sub>
138. **Réminiac F.** Are central venous lactate and arterial lactate interchangeable? A human retrospective study. *Anesth Analg* 2012. [PMID 22745117](https://pubmed.ncbi.nlm.nih.gov/22745117/)  
   <sub>⌕ `mixed venous oxygen saturation critical threshold lactate`</sub>
139. **Burchert H.** Revisiting cardiac output estimated noninvasively from oxygen uptake during exercise: an exploratory hypothesis-generating replication study. *Am J Physiol Heart Circ Physiol* 2023. [PMID 37505473](https://pubmed.ncbi.nlm.nih.gov/37505473/)  
   <sub>⌕ `mixed venous oxygen saturation critical threshold lactate`</sub>
140. **Schumacker PT.** Oxygen delivery and uptake by peripheral tissues: physiology and pathophysiology. *Crit Care Clin* 1989. [PMID 2650817](https://pubmed.ncbi.nlm.nih.gov/2650817/)  
   <sub>⌕ `oxygen extraction ratio shock tissue hypoxia review`</sub>
141. **Vincent JL.** Septic shock: particular type of acute circulatory failure. *Crit Care Med* 1990. [PMID 2403516](https://pubmed.ncbi.nlm.nih.gov/2403516/)  
   <sub>⌕ `oxygen extraction ratio shock tissue hypoxia review`</sub>
142. **Peng X.** Research progress of diagnostic and therapeutic value of carbon dioxide-derived indicators in patients with sepsis. *Zhonghua Wei Zhong Bing Ji Jiu Yi Xue* 2024. [PMID 38813642](https://pubmed.ncbi.nlm.nih.gov/38813642/)  
   <sub>⌕ `oxygen extraction ratio shock tissue hypoxia review`</sub>
143. **Convertino VA.** Wearable technology for compensatory reserve to sense hypovolemia. *J Appl Physiol (1985)* 2018. [PMID 28751369](https://pubmed.ncbi.nlm.nih.gov/28751369/)  
   <sub>⌕ `oxygen extraction ratio shock tissue hypoxia review`</sub>

## P. 신생아 · 소아

> 영아 취약성(HbF, 미성숙 CYB5R, 알칼리성 장내 세균총).

144. **Centers for Disease Control and Prevention (CDC).** Methemoglobinemia in an infant--Wisconsin, 1992. *MMWR Morb Mortal Wkly Rep* 1993. [PMID 8450825](https://pubmed.ncbi.nlm.nih.gov/8450825/)  
   <sub>⌕ `infant methemoglobinemia diarrhea acidosis`</sub>
145. **Gebara BM.** Life-threatening methemoglobinemia in infants with diarrhea and acidosis. *Clin Pediatr (Phila)* 1994. [PMID 8200173](https://pubmed.ncbi.nlm.nih.gov/8200173/)  
   <sub>⌕ `infant methemoglobinemia diarrhea acidosis`</sub>
146. **May RB.** An infant with sepsis and methemoglobinemia. *J Emerg Med* 1985. [PMID 4093577](https://pubmed.ncbi.nlm.nih.gov/4093577/)  
   <sub>⌕ `infant methemoglobinemia diarrhea acidosis`</sub>
147. **Yano SS.** Transient methemoglobinemia with acidosis in infants. *J Pediatr* 1982. [PMID 7062175](https://pubmed.ncbi.nlm.nih.gov/7062175/)  
   <sub>⌕ `infant methemoglobinemia diarrhea acidosis`</sub>

## Q. 중환자 · 약물유발 사례

> 기타 원인 약물 — 지도의 노출 클러스터를 채우는 근거.

148. **Ahmed M.** Rasburicase-Induced Methemoglobinemia. *Cureus* 2021. [PMID 33987056](https://pubmed.ncbi.nlm.nih.gov/33987056/)  
   <sub>⌕ `rasburicase methemoglobinemia G6PD`</sub>
149. **Peshin S.** Risky Remedy: Rasburicase-Induced Methemoglobinemia in Tumor Lysis Syndrome Complicated by G6PD Deficiency. *Cureus* 2024. [PMID 39677280](https://pubmed.ncbi.nlm.nih.gov/39677280/)  
   <sub>⌕ `rasburicase methemoglobinemia G6PD`</sub>
150. **Vidhyashree BH.** Rasburicase induced methemoglobinemia: A systematic review of descriptive studies. *J Oncol Pharm Pract* 2022. [PMID 35119341](https://pubmed.ncbi.nlm.nih.gov/35119341/)  
   <sub>⌕ `rasburicase methemoglobinemia G6PD`</sub>
151. **Surya N.** Rasburicase-induced methemoglobinemia and haemolytic anaemia in a patient with G6PD deficiency. *BMJ Case Rep* 2025. [PMID 40164481](https://pubmed.ncbi.nlm.nih.gov/40164481/)  
   <sub>⌕ `rasburicase methemoglobinemia G6PD`</sub>
152. **Flume PA.** Inhaled nitric oxide for adults with pulmonary non-tuberculous mycobacterial infection. *Respir Med* 2023. [PMID 36493605](https://pubmed.ncbi.nlm.nih.gov/36493605/)  
   <sub>⌕ `inhaled nitric oxide methemoglobin monitoring`</sub>
153. **Greenough A.** Inhaled nitric oxide in the neonatal period. *Expert Opin Investig Drugs* 2000. [PMID 11060764](https://pubmed.ncbi.nlm.nih.gov/11060764/)  
   <sub>⌕ `inhaled nitric oxide methemoglobin monitoring`</sub>
154. **Wessel DL.** Delivery and monitoring of inhaled nitric oxide in patients with pulmonary hypertension. *Crit Care Med* 1994. [PMID 8205825](https://pubmed.ncbi.nlm.nih.gov/8205825/)  
   <sub>⌕ `inhaled nitric oxide methemoglobin monitoring`</sub>
155. **Bruno G.** Methemoglobin and Nitrogen Dioxide Levels During Inhaled Nitric Oxide Therapy at 300 Parts per Million in Healthy Volunteers. *Respir Care* 2026. [PMID 41823202](https://pubmed.ncbi.nlm.nih.gov/41823202/)  
   <sub>⌕ `inhaled nitric oxide methemoglobin monitoring`</sub>
156. **Joseph J.** Rare cause of dyspnoea: phenazopyridine-induced methemoglobinemia. *BMJ Support Palliat Care* 2024. [PMID 38167589](https://pubmed.ncbi.nlm.nih.gov/38167589/)  
   <sub>⌕ `phenazopyridine methemoglobinemia`</sub>
157. **Koch AA.** Phenazopyridine-Induced Methemoglobinemia: A Case Report. *Cureus* 2023. [PMID 36788851](https://pubmed.ncbi.nlm.nih.gov/36788851/)  
   <sub>⌕ `phenazopyridine methemoglobinemia`</sub>
158. **Murphy T.** Acquired methemoglobinemia from phenazopyridine use. *Int J Emerg Med* 2018. [PMID 31179913](https://pubmed.ncbi.nlm.nih.gov/31179913/)  
   <sub>⌕ `phenazopyridine methemoglobinemia`</sub>
159. **Terrell JR.** Phenazopyridine-induced methemoglobinemia. *Drug Intell Clin Pharm* 1988. [PMID 3234262](https://pubmed.ncbi.nlm.nih.gov/3234262/)  
   <sub>⌕ `phenazopyridine methemoglobinemia`</sub>

## O. QSP · 모델링 방법론

> 모델링 도구와 방법론.

160. **Ball K.** Strategies for clinical dose optimization of T cell-engaging therapies in oncology. *MAbs* 2023. [PMID 36823042](https://pubmed.ncbi.nlm.nih.gov/36823042/)  
   <sub>⌕ `quantitative systems pharmacology model review drug development`</sub>
161. **Bai JP.** Quantitative Systems Pharmacology for Rare Disease Drug Development. *J Pharm Sci* 2023. [PMID 37422281](https://pubmed.ncbi.nlm.nih.gov/37422281/)  
   <sub>⌕ `quantitative systems pharmacology model review drug development`</sub>
162. **Helmlinger G.** Quantitative Systems Pharmacology: An Exemplar Model-Building Workflow With Applications in Cardiovascular, Metabolic, and Oncology Drug Development. *CPT Pharmacometrics Syst Pharmacol* 2019. [PMID 31087533](https://pubmed.ncbi.nlm.nih.gov/31087533/)  
   <sub>⌕ `quantitative systems pharmacology model review drug development`</sub>
163. **GBD 2021 Stroke Risk Factor Collaborators.** Global, regional, and national burden of stroke and its risk factors, 1990-2021: a systematic analysis for the Global Burden of Disease Study 2021. *Lancet Neurol* 2024. [PMID 39304265](https://pubmed.ncbi.nlm.nih.gov/39304265/)  
   <sub>⌕ `quantitative systems pharmacology model review drug development`</sub>
164. **Elmokadem A.** Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial. *CPT Pharmacometrics Syst Pharmacol* 2019. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)  
   <sub>⌕ `mrgsolve simulation pharmacometrics R`</sub>
165. **Lu T.** gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve. *CPT Pharmacometrics Syst Pharmacol* 2024. [PMID 38082557](https://pubmed.ncbi.nlm.nih.gov/38082557/)  
   <sub>⌕ `mrgsolve simulation pharmacometrics R`</sub>
166. **Hooijmaijers R.** Building an adaptive dose simulation framework to aid dose and schedule selection. *CPT Pharmacometrics Syst Pharmacol* 2023. [PMID 37574587](https://pubmed.ncbi.nlm.nih.gov/37574587/)  
   <sub>⌕ `mrgsolve simulation pharmacometrics R`</sub>
167. **Ponthier L.** Application of machine-learning models to predict the ganciclovir and valganciclovir exposure in children using a limited sampling strategy. *Antimicrob Agents Chemother* 2024. [PMID 39194260](https://pubmed.ncbi.nlm.nih.gov/39194260/)  
   <sub>⌕ `mrgsolve simulation pharmacometrics R`</sub>
168. **Patton JN.** Numerical simulation of oxygen delivery to muscle tissue in the presence of hemoglobin-based oxygen carriers. *Biotechnol Prog* 2006. [PMID 16889379](https://pubmed.ncbi.nlm.nih.gov/16889379/)  
   <sub>⌕ `physiologically based model oxygen transport hemoglobin simulation`</sub>
169. **Cuddington C.** Next-generation polymerized human hemoglobins in hepatic bioreactor simulations. *Biotechnol Prog* 2020. [PMID 31922354](https://pubmed.ncbi.nlm.nih.gov/31922354/)  
   <sub>⌕ `physiologically based model oxygen transport hemoglobin simulation`</sub>
170. **Kaesler A.** Computational Modeling of Oxygen Transfer in Artificial Lungs. *Artif Organs* 2018. [PMID 30043394](https://pubmed.ncbi.nlm.nih.gov/30043394/)  
   <sub>⌕ `physiologically based model oxygen transport hemoglobin simulation`</sub>
171. **Minakawa M.** Measuring and imaging of transcutaneous bilirubin, hemoglobin, and melanin based on diffuse reflectance spectroscopy. *J Biomed Opt* 2023. [PMID 37915398](https://pubmed.ncbi.nlm.nih.gov/37915398/)  
   <sub>⌕ `physiologically based model oxygen transport hemoglobin simulation`</sub>

---

## 문헌이 아닌 것 (Not literature)

아래는 문헌이 아니라 이 모델 안에서 **계산된** 것이므로 인용 대신 재현 방법을 적어 둡니다.

| 결과 | 재현 방법 |
|---|---|
| MetHb 30% ≡ 빈혈 8.31 g/dL | `python3 mhb_reference_check.py --all` §2 |
| 85% 바닥과 감도 붕괴 | 같은 파일 §4 |
| 메틸렌블루 용량–반응의 뒤집힌 U | §6 |
| G6PD 구배(작동 → 실패 → 해악) | §7 |
| 답손 반동과 누적 천장 | §8 |
| 만성 25% vs 급성 25% | §10 |
| 보에레마 산술(HBO) | §11 |

---

## 라이선스 · 면책

본 모델과 문헌 목록은 **교육 및 연구 목적**입니다. 임상 의사결정에 직접 사용해서는
안 됩니다. 특히 메틸렌블루 투여, 누적 용량, G6PD 결핍에서의 사용 여부는 반드시
최신 임상 지침과 독성학 자문에 따라야 합니다.
