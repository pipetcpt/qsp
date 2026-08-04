# G6PD 결핍증 QSP 모델 — 참고문헌
# Glucose-6-Phosphate Dehydrogenase Deficiency — Reference List

이 모델의 모든 구조적 주장과 파라미터 중심값의 출처. 섹션 1이 가장 중요하다 —
이 모델 전체가 "효소 활성은 숫자가 아니라 **세포 연령의 함수**이고, 변이란 그
함수의 **기울기**"라는 한 줄에 얹혀 있고, 그 기울기를 실제로 사람에서 측정한
논문이 섹션 1에 있다.

Every structural claim and central parameter value in
`g6pd_mrgsolve_model.R` traces to something here. Section 1 carries the
model's entire thesis: the papers that measured, *in vivo*, the rate at
which G6PD decays inside a red cell that can no longer make any.

---

## 1. 모델의 핵심 주장 — 효소는 세포 연령의 함수다 (THE THESIS)

**이 섹션의 처음 두 논문이 모델의 `E0`와 `TAU` 파라미터 그 자체다.**

1. **Piomelli S, Corash LM, Davenport DD, Miraglia J, Amorosi EL.**
   In vivo lability of glucose-6-phosphate dehydrogenase in GdA− and
   Gd Mediterranean deficiency. *J Clin Invest.* 1968;47(4):940–948.
   <https://pubmed.ncbi.nlm.nih.gov/5641628/>
   — 이 모델의 `TAU`. 밀도 원심분리로 적혈구를 연령별로 분리해 효소 활성을
   측정, 정상 B형 반감기 약 62일 대 A− 약 13일을 직접 보였다. 성숙 적혈구는
   단백질을 새로 만들 수 없으므로 효소는 오직 감소만 하고, 그래서 활성은
   연령의 지수함수가 된다. **모델의 `E(a) = E0·exp(−ln2·a/TAU)` 한 줄이
   전부 이 논문에서 나온다.**

2. **Dern RJ, Beutler E, Alving AS.** The hemolytic effect of primaquine.
   II. The natural course of the hemolytic anemia and the mechanism of its
   self-limited character. *J Lab Clin Med.* 1954;44(2):171–176.
   <https://pubmed.ncbi.nlm.nih.gov/13184228/>
   — 시나리오 1의 원본. 프리마퀸 30 mg/일을 **끊지 않고** 계속 투여했는데도
   혈색소가 7–8일에 바닥을 치고 **약을 먹는 중에 회복**했다. 저자들이 이미
   1954년에 그 이유를 정확히 지목했다: 살아남은 세포가 젊고, 새로 나온
   망상적혈구는 효소가 가득하다. 모델은 이 곡선을 이 현상을 아는 파라미터
   하나 없이 재현한다.

3. **Beutler E.** Glucose-6-phosphate dehydrogenase deficiency: a historical
   perspective. *Blood.* 2008;111(1):16–24.
   <https://pubmed.ncbi.nlm.nih.gov/18156501/>

4. **Morelli A, Benatti U, Gaetani GF, De Flora A.** Biochemical mechanisms
   of glucose-6-phosphate dehydrogenase deficiency. *Proc Natl Acad Sci USA.*
   1978;75(4):1979–1983.
   <https://pubmed.ncbi.nlm.nih.gov/273925/>
   — 변이 효소의 **불안정성**(빠른 TAU)과 **낮은 초기 활성**(낮은 E0)이
   서로 다른 생화학적 결함임을 분리해 보인 논문. 모델이 두 파라미터를
   따로 두는 이유.

5. **Alving AS, Johnson CF, Tarlov AR, Brewer GJ, Kellermeyer RW, Carson PE.**
   Mitigation of the haemolytic effect of primaquine and enhancement of its
   action against exoerythrocytic forms of the Chesson strain of *Plasmodium
   vivax* by intermittent regimens of drug administration.
   *Bull World Health Organ.* 1960;22:621–631.
   <https://pubmed.ncbi.nlm.nih.gov/13793053/>
   — 시나리오 5(주 1회 요법)의 근거. 주간 투여가 단순한 감량이 아니라 **다른
   기전**이라는 것 — 매 펄스가 가장 늙은 얇은 층만 깎고, 7일의 간격이
   망상적혈구가 그 자리를 채울 시간을 준다.

---

## 2. 종설 · 임상 전체 그림 (Comprehensive reviews)

6. **Luzzatto L, Ally M, Notaro R.** Glucose-6-phosphate dehydrogenase
   deficiency. *Blood.* 2020;136(11):1225–1240.
   <https://pubmed.ncbi.nlm.nih.gov/32702756/>
   — 현재 표준 종설.

7. **Luzzatto L, Arese P.** Favism and glucose-6-phosphate dehydrogenase
   deficiency. *N Engl J Med.* 2018;378(1):60–71.
   <https://pubmed.ncbi.nlm.nih.gov/29298156/>
   — 시나리오 6(잠두중독)의 임상 서술과 시간 경과(수시간~48시간).

8. **Cappellini MD, Fiorelli G.** Glucose-6-phosphate dehydrogenase
   deficiency. *Lancet.* 2008;371(9606):64–74.
   <https://pubmed.ncbi.nlm.nih.gov/18177777/>

9. **Beutler E.** G6PD deficiency. *Blood.* 1994;84(11):3613–3636.
   <https://pubmed.ncbi.nlm.nih.gov/7949118/>

10. **Bancone G, Chu CS.** G6PD variants and haemolytic sensitivity to
    primaquine and other drugs. *Front Pharmacol.* 2021;12:638885.
    <https://pubmed.ncbi.nlm.nih.gov/33790795/>
    — 변이별 용혈 감수성이 "활성 %"만으로 예측되지 않는다는 점을 정면으로
    다룬다. 모델의 A− 대 Mediterranean 대조(시나리오 1 vs 2)와 같은 주장.

---

## 3. 유전학 · 변이 분류 · 역학 (Genetics, classification, epidemiology)

11. **WHO Working Group.** Glucose-6-phosphate dehydrogenase deficiency.
    *Bull World Health Organ.* 1989;67(6):601–611.
    <https://pubmed.ncbi.nlm.nih.gov/2633878/>
    — 고전적 Class I–V 분류.

12. **Luzzatto L, Bancone G, Dugué P-A, et al.** New WHO classification of
    genetic variants causing G6PD deficiency.
    *Bull World Health Organ.* 2024;102(8):615–617.
    <https://pubmed.ncbi.nlm.nih.gov/39070596/>
    — 2022 개정 분류(활성 % 기준). 모델의 `VARIANTS` 목록이 따르는 틀.

13. **Howes RE, Piel FB, Patil AP, et al.** G6PD deficiency prevalence and
    estimated affected population size: a geostatistical model-based map and
    population estimates. *PLoS Med.* 2012;9(11):e1001339.
    <https://pubmed.ncbi.nlm.nih.gov/23152723/>
    — 전 세계 약 4억 명.

14. **Nkhoma ET, Poole C, Vannappagari V, Hall SA, Beutler E.** The global
    prevalence of glucose-6-phosphate dehydrogenase deficiency: a systematic
    review and meta-analysis. *Blood Cells Mol Dis.* 2009;42(3):267–278.
    <https://pubmed.ncbi.nlm.nih.gov/19233695/>

15. **Minucci A, Moradkhani K, Hwang MJ, Zuppi C, Giardina B, Capoluongo E.**
    Glucose-6-phosphate dehydrogenase (G6PD) mutations database: review of
    the "old" and update of the new mutations.
    *Blood Cells Mol Dis.* 2012;48(3):154–165.
    <https://pubmed.ncbi.nlm.nih.gov/22293322/>

16. **Beutler E, Vulliamy TJ.** Hematologically important mutations:
    glucose-6-phosphate dehydrogenase.
    *Blood Cells Mol Dis.* 2002;28(2):93–103.
    <https://pubmed.ncbi.nlm.nih.gov/12064901/>

17. **Ruwende C, Khoo SC, Snow RW, et al.** Natural selection of hemi- and
    heterozygotes for G6PD deficiency in Africa by resistance to severe
    malaria. *Nature.* 1995;376(6537):246–249.
    <https://pubmed.ncbi.nlm.nih.gov/7617034/>
    — 왜 이 대립유전자가 그렇게 흔한가.

---

## 4. 여성 이형접합자 · 라이온화 (Heterozygous females, X-inactivation)

**모델의 `mix_mosaic()`가 존재하는 이유.** 이형접합 여성은 하나의 표현형이
아니라 **두 개의 적혈구 집단**이며, 전혈 검사는 그 둘을 평균낸다.

18. **Chu CS, Bancone G, Nosten F, White NJ, Luzzatto L.** Primaquine-induced
    haemolysis in females heterozygous for G6PD deficiency.
    *Malar J.* 2018;17(1):101.
    <https://pubmed.ncbi.nlm.nih.gov/29499733/>

19. **Chu CS, Bancone G, Moore KA, et al.** Haemolysis in G6PD heterozygous
    females treated with primaquine for *Plasmodium vivax* malaria: a nested
    cohort in a trial of radical curative regimens.
    *PLoS Med.* 2017;14(2):e1002224.
    <https://pubmed.ncbi.nlm.nih.gov/28170391/>

20. **Rueangweerayut R, Bancone G, Harrell EJ, et al.** Hemolytic potential of
    tafenoquine in female volunteers heterozygous for glucose-6-phosphate
    dehydrogenase (G6PD) deficiency (G6PD Mahidol variant) versus G6PD-normal
    volunteers. *Am J Trop Med Hyg.* 2017;97(3):702–711.
    <https://pubmed.ncbi.nlm.nih.gov/28749773/>
    — 시나리오 11의 임상 대응물. 활성 40–60%인 이형접합 여성에서 타페노퀸
    300 mg 단회 투여 후 혈색소 강하가 관찰되었다.

21. **Bancone G, Kalnoky M, Chu CS, et al.** The G6PD flow-cytometric assay is
    a reliable tool for diagnosis of G6PD deficiency in women and anaemic
    subjects. *Sci Rep.* 2017;7(1):9822.
    <https://pubmed.ncbi.nlm.nih.gov/28852037/>
    — **분포**를 보는 유일한 검사. 평균만 보는 정량 검사가 왜 여성에서
    실패하는지에 대한 직접적 근거.

---

## 5. 산화 손상의 생화학 — GSH · 헤미크롬 · 밴드3 (Redox biochemistry)

22. **Arese P, Gallo V, Pantaleo A, Turrini F.** Life and death of
    glucose-6-phosphate dehydrogenase (G6PD) deficient erythrocytes — role of
    redox stress and band 3 modifications.
    *Transfus Med Hemother.* 2012;39(5):328–334.
    <https://pubmed.ncbi.nlm.nih.gov/23801924/>
    — 모델 클러스터 6–8의 뼈대: GSH 고갈 → 헤미크롬 → 하인츠소체 →
    밴드3 클러스터링 → 자가 IgG/보체 → 대식세포 제거.

23. **Turrini F, Arese P, Yuan J, Low PS.** Clustering of integral membrane
    proteins of the human erythrocyte membrane stimulates autologous IgG
    binding, complement deposition, and phagocytosis.
    *J Biol Chem.* 1991;266(35):23611–23617.
    <https://pubmed.ncbi.nlm.nih.gov/1748639/>

24. **Low FM, Hampton MB, Winterbourn CC.** Peroxiredoxin 2 and peroxide
    metabolism in the erythrocyte.
    *Antioxid Redox Signal.* 2008;10(9):1621–1630.
    <https://pubmed.ncbi.nlm.nih.gov/18479200/>
    — 모델 `VMAXOX`(정상 젊은 세포의 과산화물 처리 용량)의 크기 규모.

25. **Gaetani GF, Ferraris AM, Rolfo M, Mangerini R, Arena S, Kirkman HN.**
    Predominant role of catalase in the disposal of hydrogen peroxide within
    human erythrocytes. *Blood.* 1996;87(4):1595–1599.
    <https://pubmed.ncbi.nlm.nih.gov/8608252/>

26. **Rifkind JM, Nagababu E.** Hemoglobin redox reactions and red blood cell
    aging. *Antioxid Redox Signal.* 2013;18(17):2274–2283.
    <https://pubmed.ncbi.nlm.nih.gov/23025272/>
    — 모델 `OXBASE`(약이 없어도 존재하는 기저 자가산화 부하).

27. **Lang F, Lang KS, Lang PA, Huber SM, Wieder T.** Mechanisms and
    significance of eryptosis. *Antioxid Redox Signal.* 2006;8(7–8):1183–1192.
    <https://pubmed.ncbi.nlm.nih.gov/16910766/>

28. **Fibach E, Rachmilewitz E.** The role of oxidative stress in hemolytic
    anemia. *Curr Mol Med.* 2008;8(7):609–619.
    <https://pubmed.ncbi.nlm.nih.gov/18991647/>

29. **Tang HY, Ho HY, Wu PR, et al.** Inability to maintain GSH pool in G6PD-
    deficient red cells causes futile AMPK activation and irreversible
    metabolic disturbance. *Antioxid Redox Signal.* 2015;22(9):744–759.
    <https://pubmed.ncbi.nlm.nih.gov/25580850/>

---

## 6. 프리마퀸 · 타페노퀸 (8-aminoquinolines: the half-life IS the toxicology)

30. **Watson J, Taylor WRJ, Menard D, Kheng S, White NJ.** Modelling
    primaquine-induced haemolysis in G6PD deficiency.
    *eLife.* 2017;6:e23061.
    <https://pubmed.ncbi.nlm.nih.gov/28553934/>
    — **선행 정량 모델.** 역시 적혈구 연령 구조를 명시적으로 다루며, 주간
    요법이 안전한 이유를 같은 방식으로 설명한다. 본 모델은 여기에 약물
    PK(다중 약물), 메트헤모글로빈/NADPH 이중 소비, 빌리루빈-핵황달,
    혈관내 용혈-AKI 축을 추가한 것에 해당한다.

31. **Baird JK.** 8-aminoquinoline therapy for latent malaria.
    *Clin Microbiol Rev.* 2019;32(4):e00011-19.
    <https://pubmed.ncbi.nlm.nih.gov/31366609/>

32. **Ashley EA, Recht J, White NJ.** Primaquine: the risks and the benefits.
    *Malar J.* 2014;13:418.
    <https://pubmed.ncbi.nlm.nih.gov/25363455/>

33. **Llanos-Cuentas A, Lacerda MVG, Hien TT, et al.** Tafenoquine versus
    primaquine to prevent relapse of *Plasmodium vivax* malaria.
    *N Engl J Med.* 2019;380(3):229–241.
    <https://pubmed.ncbi.nlm.nih.gov/30650326/>

34. **Lacerda MVG, Llanos-Cuentas A, Krudsood S, et al.** Single-dose
    tafenoquine to prevent relapse of *Plasmodium vivax* malaria.
    *N Engl J Med.* 2019;380(3):215–228.
    <https://pubmed.ncbi.nlm.nih.gov/30650322/>
    — 시나리오 4a의 임상 맥락. 타페노퀸 300 mg 단회, 반감기 약 15일.

35. **Kheng S, Muth S, Taylor WRJ, et al.** Tolerability and safety of weekly
    primaquine against relapse of *Plasmodium vivax* in Cambodians with
    glucose-6-phosphate dehydrogenase deficiency.
    *BMC Med.* 2015;13:203.
    <https://pubmed.ncbi.nlm.nih.gov/26303162/>
    — 시나리오 5의 전향적 임상 근거.

36. **Pybus BS, Marcsisin SR, Jin X, et al.** The metabolism of primaquine to
    its active metabolite is dependent on CYP 2D6.
    *Malar J.* 2013;12:212.
    <https://pubmed.ncbi.nlm.nih.gov/23782898/>
    — 모델 `AS2D6`. 같은 효소가 **치료 효과와 용혈을 동시에** 만들기 때문에
    CYP2D6 불량 대사자는 용혈도 적고 근치도 실패한다.

37. **Bennett JW, Pybus BS, Yadava A, et al.** Primaquine failure and
    cytochrome P-450 2D6 in *Plasmodium vivax* malaria.
    *N Engl J Med.* 2013;369(14):1381–1382.
    <https://pubmed.ncbi.nlm.nih.gov/24088113/>

38. **Mihaly GW, Ward SA, Edwards G, Orme ML, Breckenridge AM.**
    Pharmacokinetics of primaquine in man: identification of the carboxylic
    acid derivative as a major plasma metabolite.
    *Br J Clin Pharmacol.* 1984;17(4):441–446.
    <https://pubmed.ncbi.nlm.nih.gov/6721990/>
    — 모델의 프리마퀸 PK 파라미터(`VPQ`, `CLPQ`, 반감기 약 7시간).

39. **Charles BG, Miller AK, Nasveld PE, Reid MG, Harris IE, Edstein MD.**
    Population pharmacokinetics of tafenoquine during malaria prophylaxis in
    healthy subjects. *Antimicrob Agents Chemother.* 2007;51(8):2709–2715.
    <https://pubmed.ncbi.nlm.nih.gov/17517845/>
    — 모델의 타페노퀸 PK(`VTQ` 약 1600 L, 반감기 약 15일).

40. **Commons RJ, Simpson JA, Thriemer K, et al.** The haematological
    consequences of *Plasmodium vivax* malaria after chloroquine treatment
    with and without primaquine: a WorldWide Antimalarial Resistance Network
    systematic review and individual patient data meta-analysis.
    *BMC Med.* 2019;17(1):151.
    <https://pubmed.ncbi.nlm.nih.gov/31366382/>

---

## 7. 메트헤모글로빈 · 메틸렌블루 — 같은 NADPH 풀의 두 번째 소비자

**모델에서 메틸렌블루 금기는 외워야 할 사실이 아니라 계산 결과다.**

41. **Rosen PJ, Johnson C, McGehee WG, Beutler E.** Failure of methylene blue
    treatment in toxic methemoglobinemia: association with glucose-6-phosphate
    dehydrogenase deficiency. *Ann Intern Med.* 1971;75(1):83–86.
    <https://pubmed.ncbi.nlm.nih.gov/5091568/>
    — 시나리오 7의 원본 임상 관찰.

42. **Wright RO, Lewander WJ, Woolf AD.** Methemoglobinemia: etiology,
    pharmacology, and clinical management.
    *Ann Emerg Med.* 1999;34(5):646–656.
    <https://pubmed.ncbi.nlm.nih.gov/10533013/>

43. **Cortazzo JA, Lichtman AD.** Methemoglobinemia: a review and
    recommendations for management.
    *J Cardiothorac Vasc Anesth.* 2014;28(4):1043–1047.
    <https://pubmed.ncbi.nlm.nih.gov/24075607/>

44. **Coleman MD, Coleman NA.** Drug-induced methaemoglobinaemia: treatment
    issues. *Drug Saf.* 1996;14(6):394–405.
    <https://pubmed.ncbi.nlm.nih.gov/8828017/>

45. **Coleman MD, Jacobus DP.** Reduction of dapsone hydroxylamine to dapsone
    during methaemoglobin formation in human erythrocytes in vitro.
    *Biochem Pharmacol.* 1993;45(5):1027–1033.
    <https://pubmed.ncbi.nlm.nih.gov/8461032/>
    — 모델 `KMETDAP`. 하나의 하이드록실아민 분자가 촉매적으로 여러 헤모글로빈을
    산화시키는 공동산화 회로.

46. **Zuidema J, Hilbers-Modderman ES, Merkus FW.** Clinical pharmacokinetics
    of dapsone. *Clin Pharmacokinet.* 1986;11(4):299–315.
    <https://pubmed.ncbi.nlm.nih.gov/3530584/>
    — 모델의 답손 PK(`VDP`, `CLDP`, 반감기 20–30시간).

---

## 8. 라스부리카제 · 종양용해증후군 — 산화제 용량은 종양이 정한다

47. **Relling MV, McDonagh EM, Chang T, et al.** Clinical Pharmacogenetics
    Implementation Consortium (CPIC) guidelines for rasburicase therapy in the
    context of G6PD deficiency genotype.
    *Clin Pharmacol Ther.* 2014;96(2):169–174.
    <https://pubmed.ncbi.nlm.nih.gov/24787449/>

48. **Browning LA, Kruse JA.** Hemolysis and methemoglobinemia secondary to
    rasburicase administration.
    *Ann Pharmacother.* 2005;39(11):1932–1935.
    <https://pubmed.ncbi.nlm.nih.gov/16204389/>
    — 요산 1몰당 과산화수소 1몰이라는 화학량론이 임상에서 어떻게 나타나는지.

49. **Sherwood GB, Paschal RD, Adamski J.** Rasburicase-induced
    methemoglobinemia: case report, literature review, and proposed treatment
    algorithm. *Clin Case Rep.* 2016;4(4):315–319.
    <https://pubmed.ncbi.nlm.nih.gov/27099721/>

50. **Coiffier B, Altman A, Pui C-H, Younes A, Cairo MS.** Guidelines for the
    management of pediatric and adult tumor lysis syndrome: an evidence-based
    review. *J Clin Oncol.* 2008;26(16):2767–2778.
    <https://pubmed.ncbi.nlm.nih.gov/18509186/>
    — 시나리오 8의 요산 부하 크기(15–25 mg/dL).

51. **Cheuk DKL, Chiang AKS, Chan GCF, Ha SY.** Urate oxidase for the
    prevention and treatment of tumour lysis syndrome in children with cancer.
    *Cochrane Database Syst Rev.* 2017;3:CD006945.
    <https://pubmed.ncbi.nlm.nih.gov/28272834/>

---

## 9. 신생아 황달 · 핵황달 — 곱셈이지 덧셈이 아니다

52. **Kaplan M, Renbaum P, Levy-Lahad E, Hammerman C, Lahad A, Beutler E.**
    Gilbert syndrome and glucose-6-phosphate dehydrogenase deficiency: a
    dose-dependent genetic interaction crucial to neonatal hyperbilirubinemia.
    *Proc Natl Acad Sci USA.* 1997;94(22):12128–12132.
    <https://pubmed.ncbi.nlm.nih.gov/9342374/>
    — **시나리오 10의 근거 논문.** 두 결함 중 하나만 있으면 안전하고 둘이
    함께 있을 때만 위험해진다는 상호작용을 직접 보였다. 모델이 생성과 포합을
    분자·분모로 나눠 쓰는 이유.

53. **Kaplan M, Hammerman C.** Glucose-6-phosphate dehydrogenase deficiency
    and severe neonatal hyperbilirubinemia: a complexity of interactions
    between genes and environment.
    *Semin Fetal Neonatal Med.* 2010;15(3):148–156.
    <https://pubmed.ncbi.nlm.nih.gov/19942984/>

54. **Watchko JF, Tiribelli C.** Bilirubin-induced neurologic damage —
    mechanisms and management approaches.
    *N Engl J Med.* 2013;369(21):2021–2030.
    <https://pubmed.ncbi.nlm.nih.gov/24256380/>
    — 자유 빌리루빈(`BFREE`)이 총 빌리루빈이 아니라 뇌로 넘어가는 종이라는
    점. 모델이 알부민 결합을 명시적으로 푸는 이유.

55. **Bhutani VK, Johnson LH.** Kernicterus in the 21st century: frequently
    asked questions. *J Perinatol.* 2009;29(Suppl 1):S20–S24.
    <https://pubmed.ncbi.nlm.nih.gov/19177056/>

56. **Slusher TM, Zamora TG, Appiah D, et al.** Burden of severe neonatal
    jaundice: a systematic review and meta-analysis.
    *BMJ Paediatr Open.* 2017;1(1):e000105.
    <https://pubmed.ncbi.nlm.nih.gov/29637134/>
    — G6PD 결핍이 전 세계 핵황달의 주요 원인이라는 역학적 근거.

57. **Bosma PJ, Chowdhury JR, Bakker C, et al.** The genetic basis of the
    reduced expression of bilirubin UDP-glucuronosyltransferase 1 in
    Gilbert's syndrome. *N Engl J Med.* 1995;333(18):1171–1175.
    <https://pubmed.ncbi.nlm.nih.gov/7565971/>
    — 모델 `FUGT`의 (TA)n 프로모터 유전형.

---

## 10. 혈관내 용혈 · 유리 헤모글로빈 · 급성 신손상

58. **Rother RP, Bell L, Hillmen P, Gladwin MT.** The clinical sequelae of
    intravascular hemolysis and extracellular plasma hemoglobin: a novel
    mechanism of human disease. *JAMA.* 2005;293(13):1653–1662.
    <https://pubmed.ncbi.nlm.nih.gov/15811985/>

59. **Schaer DJ, Buehler PW, Alayash AI, Belcher JD, Vercellotti GM.**
    Hemolysis and free hemoglobin revisited: exploring hemoglobin and hemin
    scavengers as a novel class of therapeutic proteins.
    *Blood.* 2013;121(8):1276–1284.
    <https://pubmed.ncbi.nlm.nih.gov/23264591/>
    — 모델의 합토글로빈 포화(`HP0`, `STOIHP`)와 그 이후의 신장 여과.

60. **Van Avondt K, Nur E, Zeerleder S.** Mechanisms of haemolysis-induced
    kidney injury. *Nat Rev Nephrol.* 2019;15(11):671–692.
    <https://pubmed.ncbi.nlm.nih.gov/31455889/>
    — 모델의 `TUB` → `GFRF` → `CREA` 축.

---

## 11. 적혈구 수명 · 적혈구 생성 동역학 (모델의 구조적 뼈대)

61. **Franco RS.** Measurement of red cell lifespan and aging.
    *Transfus Med Hemother.* 2012;39(5):302–307.
    <https://pubmed.ncbi.nlm.nih.gov/23801920/>
    — 120일 수명과 그 분포. 모델이 12개 통과구획(Erlang) 사슬을 쓰는 근거와
    그 근사의 한계.

62. **Higgins JM, Mahadevan L.** Physiological and pathological population
    dynamics of circulating human red blood cells.
    *Proc Natl Acad Sci USA.* 2010;107(47):20587–20592.
    <https://pubmed.ncbi.nlm.nih.gov/21059904/>
    — 적혈구를 연령 구조 집단으로 다루는 정량적 틀.

63. **Cohen RM, Franco RS, Khera PK, et al.** Red cell life span heterogeneity
    in hematologically normal people is sufficient to alter HbA1c.
    *Blood.* 2008;112(10):4284–4291.
    <https://pubmed.ncbi.nlm.nih.gov/18694998/>

64. **Krzyzanski W, Perez-Ruixo JJ.** An assessment of recombinant human
    erythropoietin effect on reticulocyte production rate and lifespan
    distribution in healthy subjects.
    *Pharm Res.* 2007;24(4):758–772.
    <https://pubmed.ncbi.nlm.nih.gov/17318416/>
    — 모델의 EPO → 생성률 · 골수 통과시간 단축 · shift 망상적혈구 구조.

65. **Fisher JW.** Erythropoietin: physiology and pharmacology update.
    *Exp Biol Med (Maywood).* 2003;228(1):1–14.
    <https://pubmed.ncbi.nlm.nih.gov/12524467/>
    — 혈색소 결핍에 대한 EPO의 로그-선형 반응(`GEPO`).

---

## 12. 진단 — 검사는 필요할 때 정확히 거짓말한다

66. **Domingo GJ, Satyagraha AW, Anvikar A, et al.** G6PD testing in support
    of treatment and elimination of malaria: recommendations for evaluation
    of G6PD tests. *Malar J.* 2013;12:391.
    <https://pubmed.ncbi.nlm.nih.gov/24188096/>

67. **LaRue N, Kahn M, Murray M, et al.** Comparison of quantitative and
    qualitative tests for glucose-6-phosphate dehydrogenase deficiency.
    *Am J Trop Med Hyg.* 2014;91(4):854–861.
    <https://pubmed.ncbi.nlm.nih.gov/25071003/>

68. **Bancone G, Chu CS, Chowwiwat N, et al.** Suitability of capillary blood
    for quantitative assessment of G6PD activity and performances of G6PD
    point-of-care tests. *Am J Trop Med Hyg.* 2015;92(4):818–824.
    <https://pubmed.ncbi.nlm.nih.gov/25646252/>

69. **Beutler E, Blume KG, Kaplan JC, Löhr GW, Ramot B, Valentine WN.**
    International Committee for Standardization in Haematology: recommended
    methods for red-cell enzyme analysis.
    *Br J Haematol.* 1977;35(2):331–340.
    <https://pubmed.ncbi.nlm.nih.gov/857853/>
    — 표준 정량법이 **전혈 평균 활성**을 보고한다는 사실. 모델의 `G6PDPCT`가
    바로 이 양을 흉내내며, 망상적혈구 증가가 그것을 어떻게 정상 범위로
    끌어올리는지 보여준다(시나리오 9).

70. **Ley B, Winasti Satyagraha A, Rahmat H, et al.** Performance of the
    Access Bio/CareStart rapid diagnostic test for the detection of
    glucose-6-phosphate dehydrogenase deficiency: a systematic review and
    meta-analysis. *PLoS Med.* 2019;16(12):e1002992.
    <https://pubmed.ncbi.nlm.nih.gov/31834890/>

---

## 13. 약물 안전성 목록 · 처방 지침

71. **Gammal RS, Pirmohamed M, Somogyi AA, et al.** Expanded Clinical
    Pharmacogenetics Implementation Consortium guideline for medication use
    in the context of G6PD genotype.
    *Clin Pharmacol Ther.* 2023;113(5):973–985.
    <https://pubmed.ncbi.nlm.nih.gov/36049896/>
    — 현행 CPIC 처방 지침. 모델의 약물 선택이 따르는 목록.

72. **Youngster I, Arcavi L, Schechmaster R, et al.** Medications and
    glucose-6-phosphate dehydrogenase deficiency: an evidence-based review.
    *Drug Saf.* 2010;33(9):713–726.
    <https://pubmed.ncbi.nlm.nih.gov/20701405/>
    — "금기 약물 목록"의 상당 부분이 근거가 약하다는 점을 정면으로 다룬
    드문 논문.

73. **World Health Organization.** *WHO Guidelines for Malaria.* Geneva: WHO.
    <https://www.who.int/publications/i/item/guidelines-for-malaria>
    — 주 1회 프리마퀸 45 mg × 8주 권고, 타페노퀸 전 정량 검사 요건.

74. **US FDA.** KRINTAFEL (tafenoquine) tablets — prescribing information.
    <https://www.accessdata.fda.gov/drugsatfda_docs/label/2018/210795s000lbl.pdf>
    — 투여 전 **정량적** G6PD 활성 검사 요건(정상의 70% 초과)의 근거 문서.

---

## 14. 잠두중독의 화학 (Favism)

75. **Arese P, De Flora A.** Pathophysiology of hemolysis in glucose-6-
    phosphate dehydrogenase deficiency.
    *Semin Hematol.* 1990;27(1):1–40.
    <https://pubmed.ncbi.nlm.nih.gov/2405497/>

76. **Chevion M, Navok T, Glaser G, Mager J.** The chemistry of favism-
    inducing compounds: the properties of isouramil and divicine and their
    reaction with glutathione.
    *Eur J Biochem.* 1982;127(2):405–409.
    <https://pubmed.ncbi.nlm.nih.gov/7140768/>
    — 모델 `SFV`. 디비신/이소우라밀이 GSH와 직접 반응하며 산화환원 순환한다.

77. **Mavelli I, Ciriolo MR, Rotilio G, De Sole P, Castorino M, Stabile A.**
    Favism: a hemolytic disease associated with increased superoxide
    dismutase and decreased glutathione peroxidase activities in red blood
    cells. *Eur J Biochem.* 1984;139(1):13–18.
    <https://pubmed.ncbi.nlm.nih.gov/6698003/>

---

## 15. 도구 (Tools)

- **mrgsolve** — R 기반 ODE/PK-PD 시뮬레이션: <https://mrgsolve.org/>
- **Graphviz** — 기계론적 지도 렌더링: <https://graphviz.org/>
- **Shiny** — 인터랙티브 대시보드: <https://shiny.posit.co/>
- **G6PD Deficiency Association 약물 목록**: <https://www.g6pd.org/>
- **PharmGKB — G6PD**: <https://www.pharmgkb.org/gene/PA28469>

---

## 이 모델의 한계에 대한 정직한 메모 (Honest limitations)

1. **적합(fitting)을 하지 않았다.** 파라미터는 위 문헌의 중심값에 고정했고,
   시나리오 1(Dern 곡선)과 시나리오 10(신생아 빌리루빈)이 보고된 범위 안에
   떨어지도록 산화제 스케일 인자(`SPQ` 등)만 조정했다. 개별 환자 데이터에
   대한 정식 추정·검증은 하지 않았다.
2. **연령 구조는 12개 구획 Erlang 근사**다. 실제 적혈구 수명 분포보다
   퍼져 있어(CV 약 29%) a\*가 칼날이 아니라 좁은 띠가 된다. 10일보다 미세한
   연령 차이는 이 모델로 해상할 수 없다.
3. **글루타티온을 준정상상태로 축약**했다. 분 단위 과정을 일 단위 질환에
   대해 대수식으로 푼 것이므로 노출 직후 수분~수십분의 동역학은 말할 수 없다.
4. **이형접합 여성은 두 번 돌려서 섞는다**(`mix_mosaic()`). 적혈구 집단에
   대해서는 정확하지만, 혈장 지표(빌리루빈·유리 Hb)는 질량 가중 재조합이라
   근사다.
5. **Class I(CNSHA) 표현형은 실제보다 가볍다.** 모델이 내는 기저 빈혈은
   경증 범위이며, 실제 Class I 환자의 만성 수혈 의존을 재현하지 않는다.
6. 교육·연구·가설 생성 목적이며 **임상 의사결정에 사용해서는 안 된다.**
