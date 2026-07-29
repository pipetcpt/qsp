# 만성 D형 간염 (Chronic Hepatitis D, HDV) — QSP 모델 참고문헌

> **인용 규칙 (citation convention).** 아래 링크는 두 종류입니다.
> - **`PMID xxxxxxxx`** — 본 세션에서 PubMed에서 **직접 조회하여 확인한** 레코드입니다.
>   저자·연도·저널·제목이 모두 검증되었습니다.
> - **`PubMed 검색`** — 확인하지 못한 PMID를 추측해서 적지 않기 위해, 대신 해당
>   논문을 찾을 수 있는 **검색 링크**를 제공합니다. 학회 발표 자료(D-LIVR 등)나
>   본문 확인이 되지 않은 항목이 여기에 해당합니다.
>
> 파라미터 값이 아래 문헌에서 온 것인지, 정상상태 조건으로부터 역산(back-solve)된
> 것인지, 아니면 모델의 **가설**인지는 `hdv_model_report.txt`의 A1/A3/A6/A13절과
> `hdv_reference_model.py`의 주석에 항목별로 명시되어 있습니다.

---

## 1. 수용체와 진입 (NTCP entry axis) — 불레비티드의 표적

1. Yan H, Zhong G, Xu G, et al. **Sodium taurocholate cotransporting polypeptide is a functional receptor for human hepatitis B and D virus.** *eLife* 2012;1:e00049. — HDV/HBV 수용체가 NTCP(SLC10A1)임을 규명한 논문. 본 모델의 진입 항 전체가 여기에 근거합니다. [PMID 23150796](https://pubmed.ncbi.nlm.nih.gov/23150796/)
2. Urban S, Neumann-Haefelin C, Lampertico P. **Hepatitis D virus in 2021: virology, immunology and new treatment approaches for a difficult-to-treat disease.** *Gut* 2021;70:1782-1794. — 진입·복제·조립·면역을 통합한 표준 리뷰. [PMID 34103404](https://pubmed.ncbi.nlm.nih.gov/34103404/)
3. Mentha N, Clément S, Negro F, Alfaiate D. **A review on hepatitis D: From virology to new therapies.** *J Adv Res* 2019;17:3-15. [PMID 31193285](https://pubmed.ncbi.nlm.nih.gov/31193285/)
4. Ni Y, Lempp FA, Mehrle S, et al. **Hepatitis B and D viruses exploit sodium taurocholate co-transporting polypeptide for species-specific entry into hepatocytes.** *Gastroenterology* 2014. — myrcludex B(불레비티드 전신)의 preS1 모방 기전. [PubMed 검색](https://pubmed.ncbi.nlm.nih.gov/?term=Ni+Lempp+hepatitis+B+D+sodium+taurocholate+species-specific+entry+hepatocytes)

## 2. 세포내 복제, RNA 편집(ADAR1), L-HDAg

5. Casey JL. **RNA editing in hepatitis delta virus.** *Curr Top Microbiol Immunol* 2006. — amber/W 부위 A→I 편집으로 S-HDAg에서 L-HDAg가 생성되는 기전. [PMID 16903221](https://pubmed.ncbi.nlm.nih.gov/16903221/)
6. Casey JL. **Control of ADAR1 editing of hepatitis delta virus RNAs.** *Curr Top Microbiol Immunol* 2012. [PMID 21732238](https://pubmed.ncbi.nlm.nih.gov/21732238/)
7. Jayan GC, Casey JL. **Increased RNA editing and inhibition of hepatitis delta virus replication by high-level expression of ADAR1 and ADAR2.** *J Virol* 2002;76:3819-3827. — 모델의 `f_edit` ↑ → L-HDAg ↑ → 복제 억제 경로. [PMID 11907222](https://pubmed.ncbi.nlm.nih.gov/11907222/)
8. Jayan GC, Casey JL. **Inhibition of hepatitis delta virus RNA editing by short inhibitory RNA-mediated knockdown of ADAR1 but not ADAR2.** *J Virol* 2002;76:12399-12404. [PMID 12414985](https://pubmed.ncbi.nlm.nih.gov/12414985/)
9. Chen R, Linnstaedt SD, Casey JL. **RNA editing and its control in hepatitis delta virus replication.** *Viruses* 2010;2:131-146. [PMID 21994604](https://pubmed.ncbi.nlm.nih.gov/21994604/)
10. Polson AG, Ley HL 3rd, Bass BL, Casey JL. **Hepatitis delta virus RNA editing is highly specific for the amber/W site and is suppressed by hepatitis delta antigen.** *Mol Cell Biol* 1998;18:1919-1926. — HDAg 자체가 편집을 억제하는 음성 피드백. [PMID 9528763](https://pubmed.ncbi.nlm.nih.gov/9528763/)
11. Jayan GC, Casey JL. **Effects of conserved RNA secondary structures on hepatitis delta virus genotype I RNA editing, replication, and virus production.** *J Virol* 2005;79:11187-11193. [PMID 16103170](https://pubmed.ncbi.nlm.nih.gov/16103170/)
12. Linnstaedt SD, Kasprzak WK, Shapiro BA, Casey JL. **The fraction of RNA that folds into the correct branched secondary structure determines hepatitis delta virus type 3 RNA editing levels.** *RNA* 2009;15:1177-1187. — 유전형 3의 편집 생물학 차이(모델은 HDV-1 기준). [PMID 19383766](https://pubmed.ncbi.nlm.nih.gov/19383766/)
13. Linnstaedt SD, Kasprzak WK, Shapiro BA, Casey JL. **The role of a metastable RNA secondary structure in hepatitis delta virus genotype III RNA editing.** *RNA* 2006;12:1521-1533. [PMID 16790843](https://pubmed.ncbi.nlm.nih.gov/16790843/)

## 3. 프레닐화(farnesylation)와 비리온 조립 — 로나파닙의 표적

14. Otto JC, Casey PJ. **The hepatitis delta virus large antigen is farnesylated both in vitro and in animal cells.** *J Biol Chem* 1996;271:4569-4572. — L-HDAg CXXQ(Cys211) 파네실화의 최초 규명. [PMID 8617711](https://pubmed.ncbi.nlm.nih.gov/8617711/)
15. Glenn JS, Marsters JC Jr, Greenberg HB. **Use of a prenylation inhibitor as a novel antiviral agent.** *J Virol* 1998;72:9303-9306. — 프레닐화 억제 = 항바이러스 전략의 개념 증명. [PMID 9765479](https://pubmed.ncbi.nlm.nih.gov/9765479/)
16. Bordier BB, Ohkanda J, Liu P, et al. **In vivo antiviral efficacy of prenylation inhibitors against hepatitis delta virus.** *J Clin Invest* 2003;112:407-414. — 모델의 `I_FT → 조립 차단` 항의 근거. [PMID 12897208](https://pubmed.ncbi.nlm.nih.gov/12897208/)
17. Heller T, Liang TJ. **Denying the wolf access to sheep's clothing.** *J Clin Invest* 2003;112:319-321. — 위 논문의 해설. "외피 공급을 끊는다"는 본 모델의 중심 은유의 출처. [PMID 12897197](https://pubmed.ncbi.nlm.nih.gov/12897197/)
18. Verrier ER, Weiss A, Bach C, et al. **Loss of hepatitis D virus infectivity upon farnesyl transferase inhibitor treatment associates with increasing RNA editing rates revealed by a new RT-ddPCR method.** *Antiviral Res* 2022;198:105250. — FTI가 편집률을 올린다는 관찰(지도의 점선 화살표). [PMID 35051490](https://pubmed.ncbi.nlm.nih.gov/35051490/)
19. Liang YJ, Teng W, Chen CL, et al. **Statin inhibits large hepatitis delta antigen-Smad3-twist-mediated epithelial-to-mesenchymal transition and hepatitis D virus secretion.** *J Biomed Sci* 2020;27:65. — 메발론산 경로가 조립 플럭스에 미치는 영향. [PMID 32434501](https://pubmed.ncbi.nlm.nih.gov/32434501/)

## 4. 세포분열·재생을 통한 지속 — 모델의 "바닥(floor)"

20. Giersch K, Bhadra OD, Volz T, et al. **Hepatitis delta virus persists during liver regeneration and is amplified through cell division both in vitro and in vivo.** *Gut* 2019;68:150-157. — **본 모델의 가장 중요한 구조적 근거.** HDV는 간세포 분열을 통과하여 딸세포로 전달되며 증폭됩니다. 이것이 진입억제제가 넘을 수 없는 바닥을 만듭니다. [PMID 29217749](https://pubmed.ncbi.nlm.nih.gov/29217749/)
21. Winer BY, Gaska JM, Ding Q, et al. **Preclinical assessment of antiviral combination therapy in a genetically humanized mouse model for hepatitis delta virus infection.** *Sci Transl Med* 2018;10:eaap9328. — 진입억제 + IFN 병용의 전임상 근거. [PMID 29950446](https://pubmed.ncbi.nlm.nih.gov/29950446/)

## 5. 바이러스 동태 수리모델 (본 모델이 직접 계승하는 계열)

22. Guedj J, Rotman Y, Cotler SJ, et al. **Understanding early serum hepatitis D virus and hepatitis B surface antigen kinetics during pegylated interferon-alpha therapy via mathematical modeling.** *Hepatology* 2014;60:1902-1910. — pegIFN의 1상/2상 동태와 HBsAg 커플링. 모델의 `eps_prod`(1상) / `delta_IFN`(2상) 분리의 근거. [PMID 25098971](https://pubmed.ncbi.nlm.nih.gov/25098971/)
23. Goyal A, Murray JM. **Dynamics of in vivo hepatitis D virus infection.** *J Theor Biol* 2016;398:9-19. — HDV 감염세포 풀과 HBsAg 의존성을 다룬 ODE 모델. [PMID 27012516](https://pubmed.ncbi.nlm.nih.gov/27012516/)
24. Canini L, Koh C, Cotler SJ, et al. **Pharmacokinetics and pharmacodynamics modeling of lonafarnib in patients with chronic hepatitis delta virus infection.** *Hepatol Commun* 2017;1:288-292. — 로나파닙 PK/PD 모델링. 본 모델의 `IC50_ft`·조립차단 효과 크기의 비교 기준. [PMID 29404459](https://pubmed.ncbi.nlm.nih.gov/29404459/)
25. Shekhtman L, Cotler SJ, Hershkovich L, et al. **Modelling hepatitis D virus RNA and HBsAg dynamics during nucleic acid polymer monotherapy suggest rapid turnover of HBsAg.** *Sci Rep* 2020;10:7837. — 혈중 HBsAg 반감기가 짧다는 추정(모델의 `c_s`). [PMID 32398799](https://pubmed.ncbi.nlm.nih.gov/32398799/)
26. Zakh R, Churkin A, Bietsch W, et al. **A Mathematical Model for early HBV and -HDV Kinetics during Anti-HDV Treatment.** *Mathematics (Basel)* 2021;9:3323. [PMID 35282153](https://pubmed.ncbi.nlm.nih.gov/35282153/)
27. Maya S, Hershkovich L, Rotman Y, et al. **Hepatitis delta virus RNA decline post-inoculation in human NTCP transgenic mice is biphasic.** *mBio* 2023;14:e0060023. — 이중지수 감쇠 구조. [PMID 37436080](https://pubmed.ncbi.nlm.nih.gov/37436080/)
28. Maya S, et al. **Hepatitis delta virus RNA decline post inoculation in human NTCP transgenic mice is biphasic.** *bioRxiv* 2023 (preprint). [PMID 36824865](https://pubmed.ncbi.nlm.nih.gov/36824865/)
29. Dahari H, Shudo E, Ribeiro RM, Perelson AS. **Modelling hepatitis C virus kinetics: the relationship between the infected cell loss rate and the final slope of viral decay.** *Antivir Ther* 2009;14:459-464. — "최종 기울기 = 감염세포 소실률"이라는 방법론적 근거. 본 모델에서 진입억제제가 1상이 없고 2상만 갖는 이유가 여기에 있습니다. [PMID 19474480](https://pubmed.ncbi.nlm.nih.gov/19474480/)

## 6. 불레비티드 (Bulevirtide) — 임상시험과 PK/담즙산

30. Wedemeyer H, Aleman S, Brunetto MR, et al. **A Phase 3, Randomized Trial of Bulevirtide in Chronic Hepatitis D.** *N Engl J Med* 2023;389:22-32. — **MYR301.** 48주 복합반응 2 mg 45% / 10 mg 48% / 무치료 2%, HDV RNA 반응 71%/76%/4%, ALT 정상화 51%/56%/12%. 본 모델의 앵커 2개(A3, A6)와 주요 예측 대조군. [PMID 37345876](https://pubmed.ncbi.nlm.nih.gov/37345876/)
31. Asselah T, Chulanov V, Lampertico P, et al. **Bulevirtide Combined with Pegylated Interferon for Chronic Hepatitis D.** *N Engl J Med* 2024;390:1935-1946. — **MYR204.** 병용의 치료 종료 후 지속반응 우위. 본 모델 A8의 시너지 예측 대조군. [PMID 38842520](https://pubmed.ncbi.nlm.nih.gov/38842520/)
32. Asselah T, Wedemeyer H, Aleman S, et al. **Bulevirtide Monotherapy Is Safe and Well Tolerated in Chronic Hepatitis Delta: An Integrated Safety Analysis of Bulevirtide Clinical Trials at Week 48.** *Liver Int* 2025. — 총담즙산 상승의 용량의존성과 무증상성. 모델 A1의 담즙산 앵커가 참조하는 안전성 데이터. [PMID 39648559](https://pubmed.ncbi.nlm.nih.gov/39648559/)
33. Lampertico P, Anolli MP, Roulot D, Wedemeyer H. **Antiviral therapy for chronic hepatitis delta: new insights from clinical trials and real-life studies.** *Gut* 2025. — 실제 진료 코호트에서 장기 불레비티드 반응이 시간에 따라 계속 증가한다는 관찰. 모델이 고원(plateau)보다 완만한 지속 감소를 예측하는 근거. [PMID 39663120](https://pubmed.ncbi.nlm.nih.gov/39663120/)
34. Sandmann L, et al. **Treatment response to bulevirtide is linked to amelioration of portal hypertension in patients with chronic hepatitis D.** *JHEP Rep* 2025. — 반응이 문맥압/혈소판으로 이어진다는 임상 근거(모델의 `PLT`·`Fib` 축). [PMID 41674894](https://pubmed.ncbi.nlm.nih.gov/41674894/)
35. Zhu V, et al. **Evaluation of the drug-drug interaction potential of the novel hepatitis B and D virus entry inhibitor bulevirtide at OATP1B in healthy volunteers.** *Front Pharmacol* 2023;14:1128547. — 불레비티드가 OATP1B를 약하게만 억제한다는 데이터. 모델에서 담즙산 보조 흡수 경로 `r_oatp`를 별도 항으로 둔 근거. [PMID 37089922](https://pubmed.ncbi.nlm.nih.gov/37089922/)
36. Bogomolov P, Alexandrov A, Voronkova N, et al. **Treatment of chronic hepatitis D with the entry inhibitor myrcludex B: First results of a phase Ib/IIa study.** *J Hepatol* 2016;65:490-498. — 최초의 개념 증명. [PubMed 검색](https://pubmed.ncbi.nlm.nih.gov/?term=Bogomolov+chronic+hepatitis+D+entry+inhibitor+myrcludex+B+phase+Ib%2FIIa)
37. Blank A, Meier K, Urban S, Haefeli WE, Blank A. **Clinical pharmacology of the HBV/HDV entry inhibitor myrcludex B / bulevirtide — pharmacokinetics and bile acid changes.** — 2 mg vs 10 mg에서 총담즙산 상승 폭. 모델 A1의 두 앵커. [PubMed 검색](https://pubmed.ncbi.nlm.nih.gov/?term=Blank+myrcludex+B+bulevirtide+pharmacokinetics+bile+acids)
38. European Medicines Agency. **Hepcludex (bulevirtide) — EPAR product information.** — 2 mg SC 1일 1회 승인 용법, 용량 초비례적(more than dose-proportional) 노출. [EMA](https://www.ema.europa.eu/en/medicines/human/EPAR/hepcludex)

## 7. 인터페론 알파 / 람다

39. Wedemeyer H, Yurdaydìn C, Dalekos GN, et al. **Peginterferon plus adefovir versus either drug alone for hepatitis delta.** *N Engl J Med* 2011;364:322-331. — **HIDIT-1.** 48주 pegIFN, 치료 종료 24주 후 HDV RNA 음성 약 26-28%. 아데포비어 추가는 무효(= NUC는 HDV에 듣지 않는다는 모델 구조와 일치). [PMID 21268724](https://pubmed.ncbi.nlm.nih.gov/21268724/)
40. Heidrich B, Yurdaydìn C, Kabaçam G, et al. **Late HDV RNA relapse after peginterferon alpha-based therapy of chronic hepatitis delta.** *Hepatology* 2014;60:87-97. — 장기 추적에서의 후기 재발. 모델의 "탈진 회복이 서서히 되돌아간다"는 항의 근거. [PMID 24585488](https://pubmed.ncbi.nlm.nih.gov/24585488/)
41. Etzion O, Hamid S, Lurie Y, et al. **Treatment of chronic hepatitis D with peginterferon lambda—the phase 2 LIMT-1 clinical trial.** *Hepatology* 2023;77:2093-2103. — 180 µg에서 치료 종료 24주 후 반응 36%(5/14), 고빌리루빈혈증 문제. 모델의 람다 arm(간세포 한정 수용체·혈액학적 독성 없음·간 독성 있음). [PMID 36800850](https://pubmed.ncbi.nlm.nih.gov/36800850/)
42. Grabowski J, Yurdaydìn C, Zachou K, et al. **Hepatitis D virus-specific cytokine responses in patients with chronic hepatitis delta before and during interferon alfa-treatment.** *Liver Int* 2011;31:1395-1405. — IFN 치료 중 HDV 특이 면역반응의 변화. 모델의 탈진 회복 항. [PMID 21762356](https://pubmed.ncbi.nlm.nih.gov/21762356/)
43. Mederacke I, Yurdaydin C, Grosshennig A, et al. **Renal function during treatment with adefovir plus peginterferon alfa-2a vs either drug alone in hepatitis B/D co-infection.** *J Viral Hepat* 2012;19:387-395. [PMID 22571900](https://pubmed.ncbi.nlm.nih.gov/22571900/)
44. Mederacke I, Bremer B, Heidrich B, et al. **Anti-HDV immunoglobulin M testing in hepatitis delta revisited: correlations with disease activity and response to pegylated interferon-alpha2a treatment.** *Antivir Ther* 2012;17:305-312. — 모델의 anti-HD IgM 활성 지표. [PMID 22293066](https://pubmed.ncbi.nlm.nih.gov/22293066/)
45. Wedemeyer H, Yurdaydin C, Hardtke S, et al. **Peginterferon alfa-2a plus tenofovir disoproxil fumarate for hepatitis D (HIDIT-II): a randomised, placebo-controlled, phase 2 trial.** *Lancet Infect Dis* 2019;19:275-286. — **HIDIT-2.** 96주로 늘려도 치료 종료 후 반응은 크게 개선되지 않음. [PubMed 검색](https://pubmed.ncbi.nlm.nih.gov/?term=HIDIT-II+peginterferon+alfa-2a+tenofovir+hepatitis+D+randomised)

## 8. 로나파닙 / 리토나비르

46. Yurdaydin C, Keskin O, Yurdcu E, et al. **A phase 2 dose-finding study of lonafarnib and ritonavir with or without interferon alpha for chronic delta hepatitis.** *Hepatology* 2022;75:1551-1565. — **LOWR HDV-2/4 계열.** 리토나비르 부스팅으로 저용량 로나파닙의 노출을 올리고 위장관 독성을 줄이는 전략. [PMID 34860418](https://pubmed.ncbi.nlm.nih.gov/34860418/)
47. Yurdaydin C, Idilman R, Kalkan Ç, et al. **Optimizing lonafarnib treatment for the management of chronic delta hepatitis: The LOWR HDV-1 study.** *Hepatology* 2018;67:1224-1236. — 용량 제한 독성이 위장관계임을 보인 연구. [PMID 29152762](https://pubmed.ncbi.nlm.nih.gov/29152762/)
48. Keskin O, Yurdaydin C. **Emerging drugs for hepatitis D.** *Expert Opin Emerg Drugs* 2023;28:1-11. [PMID 37096555](https://pubmed.ncbi.nlm.nih.gov/37096555/)
49. Etzion O, Asselah T, Yurdaydin C, et al. **D-LIVR: a phase 3 study of lonafarnib boosted with ritonavir with and without peginterferon alfa in chronic hepatitis delta.** — 48주 복합반응이 로나파닙/리토나비르 약 10%, 삼제 병용 약 19%, 위약 약 2%로 보고된 3상. 학회 발표 중심 자료이므로 검색 링크로 제공합니다. [PubMed 검색](https://pubmed.ncbi.nlm.nih.gov/?term=D-LIVR+lonafarnib+ritonavir+peginterferon+hepatitis+delta+phase+3)

## 9. 역학, 자연경과, 임상 결과

50. Stockdale AJ, Kreuels B, Henrion MYR, et al. **The global prevalence of hepatitis D virus infection: Systematic review and meta-analysis.** *J Hepatol* 2020;73:523-532. — HBsAg 양성자 중 항HDV 양성률 추정. [PMID 32335166](https://pubmed.ncbi.nlm.nih.gov/32335166/)
51. Chen HY, Shen DT, Ji DZ, et al. **Prevalence and burden of hepatitis D virus infection in the global population: a systematic review and meta-analysis.** *Gut* 2019;68:512-521. [PMID 30228220](https://pubmed.ncbi.nlm.nih.gov/30228220/)
52. Gish RG, Wong RJ, Di Tanna GL, et al. **Association of hepatitis delta virus with liver morbidity and mortality: A systematic literature review and meta-analysis.** *Hepatology* 2024;79:1129-1140. — HDV 동반 시 간경변·비대상성·HCC·사망 위험 증가폭. 모델 A11의 대조 기준. [PMID 37870278](https://pubmed.ncbi.nlm.nih.gov/37870278/)
53. Stockdale AJ, Chaponda M, Beloukas A, et al. **Prevalence of hepatitis D virus infection in sub-Saharan Africa: a systematic review and meta-analysis.** *Lancet Glob Health* 2017;5:e992-e1003. [PMID 28911765](https://pubmed.ncbi.nlm.nih.gov/28911765/)
54. Wong RJ, Kaufman HW, Niles JK, et al. **Estimating the prevalence of hepatitis delta virus infection among adults in the United States: A meta-analysis.** *Liver Int* 2024;44:1715-1723. [PMID 38563728](https://pubmed.ncbi.nlm.nih.gov/38563728/)
55. Shen DT, Han J, Ji DZ, et al. **Epidemiology estimates of hepatitis D in individuals co-infected with human immunodeficiency virus and hepatitis B virus, 2002-2018: A systematic review and meta-analysis.** *J Viral Hepat* 2021;28:1057-1067. [PMID 33877742](https://pubmed.ncbi.nlm.nih.gov/33877742/)
56. Negro F, Lok AS. **Hepatitis D: A Review.** *JAMA* 2023;330:2376-2387. — 임상 표준 리뷰(진단·치료·자연경과). [PMID 37943548](https://pubmed.ncbi.nlm.nih.gov/37943548/)
57. Degasperi E, Anolli MP, Lampertico P, et al. **Hepatitis D Virus Infection: Pathophysiology, Epidemiology and Treatment. Report From the Third Delta Cure Meeting 2024.** *Liver Int* 2025. — 최신 개발 파이프라인 요약(siRNA, NAP, HBsAg 표적). [PMID 40540405](https://pubmed.ncbi.nlm.nih.gov/40540405/)
58. Rizzetto M. **Chronic Hepatitis D; at a Standstill?** *Dig Dis* 2016;34:303-307. — HDV를 처음 기술한 저자의 관점. [PMID 27170382](https://pubmed.ncbi.nlm.nih.gov/27170382/)
59. Abbas Z, Abbas M. **Management of hepatitis delta: Need for novel therapeutic options.** *World J Gastroenterol* 2015;21:9461-9465. [PMID 26327754](https://pubmed.ncbi.nlm.nih.gov/26327754/)
60. Alfaiate D, Clément S, Gomes D, Goossens N, Negro F. **Chronic hepatitis D and hepatocellular carcinoma: A systematic review and meta-analysis of observational studies.** *J Hepatol* 2020;73:533-539. — HDV의 HCC 위험 증가. 모델 `bF_hcc`·`bV_hcc`의 근거. [PubMed 검색](https://pubmed.ncbi.nlm.nih.gov/?term=Alfaiate+chronic+hepatitis+D+hepatocellular+carcinoma+systematic+review+meta-analysis)

## 10. HBsAg 표적 치료 (외피 공급 차단)

61. Bazinet M, Pântea V, Cebotarescu V, et al. **Safety and efficacy of REP 2139 and pegylated interferon alfa-2a for treatment-naive patients with chronic hepatitis B virus and hepatitis D virus co-infection (REP 301 and REP 301-LTF).** *Lancet Gastroenterol Hepatol* 2017. — 핵산 폴리머(NAP)로 HBsAg를 낮추는 전략. [PubMed 검색](https://pubmed.ncbi.nlm.nih.gov/?term=Bazinet+REP+2139+pegylated+interferon+hepatitis+B+D+co-infection+REP+301)
62. **Elebsiran (VIR-2218) / bepirovirsen 등 HBsAg 표적제제의 HDV 적용.** — siRNA·ASO로 HBsAg를 낮추어 HDV 외피 공급을 끊는 접근. 본 모델 A9의 `siRNA` arm. [PubMed 검색](https://pubmed.ncbi.nlm.nih.gov/?term=elebsiran+OR+VIR-2218+OR+bepirovirsen+hepatitis+delta+HBsAg)

---

## 파라미터 출처 요약 (Where each number comes from)

| 모델 항목 | 출처 유형 | 근거 |
|---|---|---|
| NTCP가 HDV 수용체 | 문헌 | #1, #2 |
| `Kd_ntcp`, `r_oatp` | **적합(2개)** | #32, #37의 총담즙산 상승 배수 2개 앵커 → 점유율은 역산으로 **도출** |
| 불레비티드 용량 초비례 PK | 문헌 | #38 (포화성 표적매개 소실로 구현) |
| `dth_Id_immune` | **적합(1개)** | #30 MYR301 2 mg 48주 바이러스 반응률 71% |
| `ALT_base` | **적합(1개)** | #30 MYR301 2 mg 48주 ALT 정상화 51% |
| 분열매개 지속(바닥) | 문헌(구조) | #20 — 정량값은 정상상태 역산 |
| 재생 국소성 `loc_renew` | 역산 | 관찰되는 HBsAg 안정성으로부터 |
| IFN 1상 `eps_prod` / 2상 `delta_IFN` | 문헌(구조) | #22, #29 |
| `IC50_ft`, 조립 차단 | 문헌 | #14-16, #24 |
| 세포내 RNA 축적(로나파닙 역설) | 모델 예측 | 구조적 귀결, #18과 방향 일치 |
| ADAR1 IFN 유도 → 편집 ↑ | 문헌 | #5-#7 |
| 탈진 회복이 지속반응을 만든다 | **모델 가설** | #40, #42와 방향 일치, 정량 근거 없음 |
| `kappa_fresh` (ALT-RNA 해리) | **모델 가설** | #30의 관찰을 설명하기 위한 최소 구조. 반증 조건은 보고서 A13 |
| 섬유화 진행/회귀 | 문헌 | #52, #34 |
| HCC 위험 | 문헌 | #52, #60 |

> 모델의 **가설(hypothesis)** 로 표시된 항목은 문헌에서 직접 측정된 값이 아니며,
> `hdv_model_report.txt` A13절에 각각의 **반증 조건(falsifier)** 을 명시했습니다.
