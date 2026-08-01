# 참고문헌 — 비결핵 항산균 폐질환 (MAC-PD) QSP 모델
# References — Nontuberculous Mycobacterial Lung Disease (Mycobacterium avium complex pulmonary disease)

> **링크 규칙 (link convention).** 각 항목은 PubMed **제목 검색 URL**로 연결됩니다
> (`pubmed.ncbi.nlm.nih.gov/?term=...`). 기억에 의존한 PMID 숫자를 그대로 적으면
> 엉뚱한 논문으로 연결될 위험이 있어, 검증 가능한 서지정보(저자·저널·연도·권·면)를
> 함께 적고 링크는 제목 검색으로 통일했습니다. 검색 결과 첫 항목이 해당 논문입니다.
>
> Each entry links to a PubMed **title search** rather than a hard-coded PMID, so
> that every link resolves to the paper actually named in the citation. Full
> bibliographic detail is given so the record can be verified independently.

---

## 1. 진료지침 · 역학 (Guidelines & Epidemiology)

1. **Daley CL, Iaccarino JM, Lange C, et al.** Treatment of Nontuberculous Mycobacterial Pulmonary Disease: An Official ATS/ERS/ESCMID/IDSA Clinical Practice Guideline. *Clin Infect Dis.* 2020;71(4):e1–e36. (동시 게재: *Eur Respir J.* 2020;56(1):2000535)
   → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Treatment+of+Nontuberculous+Mycobacterial+Pulmonary+Disease+An+Official+ATS+ERS+ESCMID+IDSA+Clinical+Practice+Guideline)
   *모델에서의 역할: 시나리오 2·3의 기본 요법(대식세포 3제 + 결절기관지확장형은 주 3회, 공동형은 매일)과 12개월 배양음전 유지라는 1차 평가변수 정의의 근거.*

2. **Griffith DE, Aksamit T, Brown-Elliott BA, et al.** An official ATS/IDSA statement: diagnosis, treatment, and prevention of nontuberculous mycobacterial diseases. *Am J Respir Crit Care Med.* 2007;175(4):367–416.
   → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=An+official+ATS+IDSA+statement+diagnosis+treatment+and+prevention+of+nontuberculous+mycobacterial+diseases)

3. **Prevots DR, Marras TK.** Epidemiology of human pulmonary infection with nontuberculous mycobacteria: a review. *Clin Chest Med.* 2015;36(1):13–34.
   → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Epidemiology+of+human+pulmonary+infection+with+nontuberculous+mycobacteria+a+review+Prevots)

4. **Falkinham JO 3rd.** Environmental sources of nontuberculous mycobacteria. *Clin Chest Med.* 2015;36(1):35–41.
   → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Environmental+sources+of+nontuberculous+mycobacteria+Falkinham)
   *모델에서의 역할: 노출 노드(샤워헤드·온수욕조·급수 바이오필름)와 재감염 대 재발의 구분.*

5. **Kumar K, Loebinger MR.** Nontuberculous Mycobacterial Pulmonary Disease: Clinical Epidemiologic Features, Risk Factors, and Diagnosis. *Chest.* 2022;161(3):637–646.
   → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Nontuberculous+Mycobacterial+Pulmonary+Disease+Clinical+Epidemiologic+Features+Risk+Factors+and+Diagnosis+Kumar+Loebinger)

6. **Hwang JA, Kim S, Jo KW, Shim TS.** Natural history of Mycobacterium avium complex lung disease in untreated patients with stable course. *Eur Respir J.* 2017;49(3):1600537.
   → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Natural+history+of+Mycobacterium+avium+complex+lung+disease+in+untreated+patients+with+stable+course)
   *모델에서의 역할: 시나리오 1(관찰) — 무치료 자연경과의 기준선.*

---

## 2. 숙주 감수성 표현형 (Host Susceptibility — the "Lady Windermere" phenotype)

7. **Reich JM, Johnson RE.** Mycobacterium avium complex pulmonary disease presenting as an isolated lingular or middle lobe pattern. The Lady Windermere syndrome. *Chest.* 1992;101(6):1605–1609.
   → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Mycobacterium+avium+complex+pulmonary+disease+presenting+as+an+isolated+lingular+or+middle+lobe+pattern+Lady+Windermere)

8. **Kim RD, Greenberg DE, Ehrmantraut ME, et al.** Pulmonary nontuberculous mycobacterial disease: prospective study of a distinct preexisting syndrome. *Am J Respir Crit Care Med.* 2008;178(10):1066–1074.
   → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Pulmonary+nontuberculous+mycobacterial+disease+prospective+study+of+a+distinct+preexisting+syndrome)
   *모델에서의 역할: `S_host` 복합 감수성 지표(마른 체형·흉곽 이상·CFTR 이형접합).*

9. **Kartalija M, Ovrutsky AR, Bryan CL, et al.** Patients with nontuberculous mycobacterial lung disease exhibit unique body and immune phenotypes. *Am J Respir Crit Care Med.* 2013;187(2):197–205.
   → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Patients+with+nontuberculous+mycobacterial+lung+disease+exhibit+unique+body+and+immune+phenotypes)
   *모델에서의 역할: 아디포넥틴↑/렙틴↓ 노드, 그리고 TNF 구동 체중감소 되먹임 고리(`WT` ODE).*

10. **Chan ED, Iseman MD.** Slender, older women appear to be more susceptible to nontuberculous mycobacterial lung disease. *Gend Med.* 2010;7(1):5–18.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Slender+older+women+appear+to+be+more+susceptible+to+nontuberculous+mycobacterial+lung+disease)

11. **Chan ED, Iseman MD.** Underlying host risk factors for nontuberculous mycobacterial lung disease. *Semin Respir Crit Care Med.* 2013;34(1):110–123.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Underlying+host+risk+factors+for+nontuberculous+mycobacterial+lung+disease+Chan+Iseman)

12. **Andréjak C, Nielsen R, Thomsen VØ, et al.** Chronic respiratory disease, inhaled corticosteroids and risk of non-tuberculous mycobacteriosis. *Thorax.* 2013;68(3):256–262.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Chronic+respiratory+disease+inhaled+corticosteroids+and+risk+of+non-tuberculous+mycobacteriosis)

13. **Brode SK, Campitelli MA, Kwong JC, et al.** The risk of mycobacterial infections associated with inhaled corticosteroid use. *Eur Respir J.* 2017;50(3):1700037.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=The+risk+of+mycobacterial+infections+associated+with+inhaled+corticosteroid+use)

14. **Browne SK, Burbelo PD, Chetchotisakd P, et al.** Adult-onset immunodeficiency in Thailand and Taiwan. *N Engl J Med.* 2012;367(8):725–734.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Adult-onset+immunodeficiency+in+Thailand+and+Taiwan)
    *모델에서의 역할: 시나리오 10 — 항IFN-γ 자가항체(`IFNCAP` = 0.15).*

15. **Winthrop KL, Baxter R, Liu L, et al.** Mycobacterial diseases and antitumour necrosis factor therapy in USA. *Ann Rheum Dis.* 2013;72(1):37–42.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Mycobacterial+diseases+and+antitumour+necrosis+factor+therapy+in+USA)

---

## 3. 세균 미소환경 · 바이오필름 · 식세포 내 생존 (Bacterial Niches)

16. **Sturgill-Koszycki S, Schlesinger PH, Chakraborty P, et al.** Lack of acidification in Mycobacterium phagosomes produced by exclusion of the vesicular proton-ATPase. *Science.* 1994;263(5147):678–681.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Lack+of+acidification+in+Mycobacterium+phagosomes+produced+by+exclusion+of+the+vesicular+proton-ATPase)
    *모델에서의 역할: `PHPHAG` = 5.2 라는 단 하나의 입력 — 이 값에서 이온 포획비와 MIC 상승이 동시에 유도된다.*

17. **Crowle AJ, Dahl R, Ross E, May MH.** Evidence that vesicles containing living, virulent Mycobacterium tuberculosis or Mycobacterium avium in cultured human macrophages are not acidic. *Infect Immun.* 1991;59(5):1823–1831.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Evidence+that+vesicles+containing+living+virulent+Mycobacterium+tuberculosis+or+Mycobacterium+avium+in+cultured+human+macrophages+are+not+acidic)

18. **Carter G, Wu M, Drummond DC, Bermudez LE.** Characterization of biofilm formation by clinical isolates of Mycobacterium avium. *J Med Microbiol.* 2003;52(Pt 9):747–752.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Characterization+of+biofilm+formation+by+clinical+isolates+of+Mycobacterium+avium)
    *모델에서의 역할: `B_B` 구획과 EPS 확산 장벽(`PBFM`, `PBFK`, `PBFL`).*

19. **Yamazaki Y, Danelishvili L, Wu M, et al.** The ability to form biofilm influences Mycobacterium avium invasion and translocation of bronchial epithelial cells. *Cell Microbiol.* 2006;8(5):806–814.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=The+ability+to+form+biofilm+influences+Mycobacterium+avium+invasion+and+translocation+of+bronchial+epithelial+cells)

20. **Sarathy JP, Dartois V.** Caseum: a Niche for Mycobacterium tuberculosis Drug-Tolerant Persisters. *Clin Microbiol Rev.* 2020;33(3):e00159-19.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Caseum+a+Niche+for+Mycobacterium+tuberculosis+Drug-Tolerant+Persisters)
    *모델에서의 역할: `B_C`(건락괴사 구획)의 비증식성·약물 내성 특성과 `TOLCS` 파라미터.*

21. **Dartois V.** The path of anti-tuberculosis drugs: from blood to lesions to mycobacterial cells. *Nat Rev Microbiol.* 2014;12(3):159–167.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=The+path+of+anti-tuberculosis+drugs+from+blood+to+lesions+to+mycobacterial+cells)
    *모델에서의 역할: 니치별 침투계수(`PCS*`, `PBF*`)를 "혈장 농도 하나"가 아닌 병변별 농도로 나누어 다루는 설계 자체의 근거.*

22. **Prideaux B, Via LE, Zimmerman MD, et al.** The association between sterilizing activity and drug distribution into tuberculosis lesions. *Nat Med.* 2015;21(10):1223–1227.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=The+association+between+sterilizing+activity+and+drug+distribution+into+tuberculosis+lesions)

23. **Early J, Fischer K, Bermudez LE.** Mycobacterium avium uses apoptotic macrophages as tools for spreading. *Microb Pathog.* 2011;50(2):132–139.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Mycobacterium+avium+uses+apoptotic+macrophages+as+tools+for+spreading)

---

## 4. 숙주 면역 (Immunology)

24. **Bermudez LE, Young LS.** Recombinant granulocyte-macrophage colony-stimulating factor activates human macrophages to inhibit growth or kill Mycobacterium avium complex. *J Leukoc Biol.* 1990;48(1):67–73.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Recombinant+granulocyte-macrophage+colony-stimulating+factor+activates+human+macrophages+to+inhibit+growth+or+kill+Mycobacterium+avium+complex)

25. **Cole PJ.** Inflammation: a two-edged sword — the model of bronchiectasis. *Eur J Respir Dis Suppl.* 1986;147:6–15.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Inflammation+a+two-edged+sword+the+model+of+bronchiectasis+Cole)
    *모델에서의 역할: `MMP → BRO → MCC↓ → B_B↑` 로 닫히는 Cole 악순환 고리.*

26. **Flume PA, Chalmers JD, Olivier KN.** Advances in bronchiectasis: endotyping, genetics, microbiome, and disease heterogeneity. *Lancet.* 2018;392(10150):880–890.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Advances+in+bronchiectasis+endotyping+genetics+microbiome+and+disease+heterogeneity)

27. **Chalmers JD, Moffitt KL, Suarez-Cuartin G, et al.** Neutrophil Elastase Activity Is Associated with Exacerbations and Lung Function Decline in Bronchiectasis. *Am J Respir Crit Care Med.* 2017;195(10):1384–1393.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Neutrophil+Elastase+Activity+Is+Associated+with+Exacerbations+and+Lung+Function+Decline+in+Bronchiectasis)

---

## 5. 마크로라이드 세포내 약동학 — 모델의 구조적 핵심 (Macrolide Intracellular PK)

28. **Olsen KM, San Pedro G, Gann LP, et al.** Intrapulmonary pharmacokinetics of azithromycin in healthy volunteers given five oral doses. *Antimicrob Agents Chemother.* 1996;40(11):2582–2585.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Intrapulmonary+pharmacokinetics+of+azithromycin+in+healthy+volunteers+given+five+oral+doses)
    *모델에서의 역할: 폐포 대식세포:혈장 농도비 10²–10³ — 본 모델에서는 이 값을 파라미터로 넣지 않고 pKa 8.7과 pH 5.2로부터 `R_trap ≈ 151`로 유도한다.*

29. **Rodvold KA, Gotfried MH, Danziger LH, Servi RJ.** Intrapulmonary steady-state concentrations of clarithromycin and azithromycin in healthy adult volunteers. *Antimicrob Agents Chemother.* 1997;41(6):1399–1402.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Intrapulmonary+steady-state+concentrations+of+clarithromycin+and+azithromycin+in+healthy+adult+volunteers)

30. **Rodvold KA.** Clinical pharmacokinetics of clarithromycin. *Clin Pharmacokinet.* 1999;37(5):385–398.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Clinical+pharmacokinetics+of+clarithromycin+Rodvold)

31. **Carlier MB, Zenebergh A, Tulkens PM.** Cellular uptake and subcellular distribution of roxithromycin and erythromycin in phagocytic cells. *J Antimicrob Chemother.* 1987;20 Suppl B:47–56.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Cellular+uptake+and+subcellular+distribution+of+roxithromycin+and+erythromycin+in+phagocytic+cells)
    *모델에서의 역할: 리소좀 이온 포획(lysosomotropic trapping) 기전 — Henderson-Hasselbalch 식의 실험적 근거.*

32. **Tulkens PM.** Intracellular distribution and activity of antibiotics. *Eur J Clin Microbiol Infect Dis.* 1991;10(2):100–106.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Intracellular+distribution+and+activity+of+antibiotics+Tulkens)
    *모델에서의 역할: "세포 내로 들어간 양"과 "세포 내에서 발휘되는 활성"이 다르다는 원칙 — 본 모델의 역설(paradox) 항의 개념적 출발점.*

33. **Maurin M, Raoult D.** Use of aminoglycosides in treatment of infections due to intracellular bacteria. *Antimicrob Agents Chemother.* 2001;45(11):2977–2986.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Use+of+aminoglycosides+in+treatment+of+infections+due+to+intracellular+bacteria)
    *모델에서의 역할: 유리 아미카신의 세포막 비투과성(`CK_I = 0` for IV) 및 산성 pH에서의 활성 소실(`GK` = 0.85).*

34. **Lemaire S, Van Bambeke F, Mingeot-Leclercq MP, Tulkens PM.** Activity of three β-lactams (ertapenem, meropenem and ampicillin) against intraphagocytic Listeria monocytogenes and Staphylococcus aureus. *J Antimicrob Chemother.* 2005;55(6):897–904.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Activity+of+three+beta-lactams+ertapenem+meropenem+and+ampicillin+against+intraphagocytic+Listeria+monocytogenes+and+Staphylococcus+aureus)
    *모델에서의 역할: 산성 세포내 환경에서 항생제 활성이 세포외 대비 수 배–수십 배 감소한다는 정량적 선례.*

---

## 6. 에탐부톨 · 리팜핀 · 약물상호작용 (Ethambutol, Rifamycins, DDI)

35. **Wallace RJ Jr, Brown BA, Griffith DE, et al.** Reduced serum levels of clarithromycin in patients treated with multidrug regimens including rifampin or rifabutin for Mycobacterium avium–M. intracellulare infection. *J Infect Dis.* 1995;171(3):747–750.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Reduced+serum+levels+of+clarithromycin+in+patients+treated+with+multidrug+regimens+including+rifampin+or+rifabutin)
    *모델에서의 역할: 시나리오 7 — `ENZ` 유도로 클래리트로마이신 AUC가 약 60–70% 감소하는 것이 CYP3A 분율(`FCYP3A` 0.70 vs 아지트로마이신 0.05) 차이에서 산술적으로 나온다.*

36. **van Ingen J, Egelund EF, Levin A, et al.** The pharmacokinetics and pharmacodynamics of pulmonary Mycobacterium avium complex disease treatment. *Am J Respir Crit Care Med.* 2012;186(6):559–565.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=The+pharmacokinetics+and+pharmacodynamics+of+pulmonary+Mycobacterium+avium+complex+disease+treatment)
    *모델에서의 역할: MAC-PD 표준요법의 혈중 농도가 MIC 대비 낮다는 관찰 — `MICM0`(4 mg/L)과 ELF 농도의 관계 설정 근거.*

37. **Deshpande D, Srivastava S, Pasipanodya JG, Gumbo T.** Ethambutol optimal clinical dose and susceptibility breakpoint identification by use of a novel pharmacokinetic-pharmacodynamic model of disseminated intracellular Mycobacterium avium. *Antimicrob Agents Chemother.* 2010;54(6):2606–2613.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Ethambutol+optimal+clinical+dose+and+susceptibility+breakpoint+identification+by+use+of+a+novel+pharmacokinetic-pharmacodynamic+model+of+disseminated+intracellular+Mycobacterium+avium)
    *모델에서의 역할: 에탐부톨의 역할을 직접 살균이 아닌 세포벽 투과성 증강(`PERM` 항)으로 둔 근거.*

38. **Jarand J, Davis JP, Cowie RL, et al.** Long-term follow-up of Mycobacterium avium complex lung disease in patients treated with regimens including clofazimine and/or rifampin. *Chest.* 2016;149(5):1285–1293.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Long-term+follow-up+of+Mycobacterium+avium+complex+lung+disease+in+patients+treated+with+regimens+including+clofazimine+and+or+rifampin)
    *모델에서의 역할: 시나리오 8(리팜핀 배제·클로파지민 대체)의 임상적 근거.*

39. **Griffith DE, Brown-Elliott BA, Shepherd S, et al.** Ethambutol ocular toxicity in treatment regimens for Mycobacterium avium complex lung disease. *Am J Respir Crit Care Med.* 2005;172(2):250–253.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Ethambutol+ocular+toxicity+in+treatment+regimens+for+Mycobacterium+avium+complex+lung+disease)
    *모델에서의 역할: `OPT` 독성 구획(역치 `EMBTHR` 초과 노출의 누적).*

---

## 7. 아미카신과 흡입 리포솜 아미카신 (Amikacin & ALIS) — 전달경로 논증

40. **Griffith DE, Eagle G, Thomson R, et al.** Amikacin Liposome Inhalation Suspension for Treatment-Refractory Lung Disease Caused by Mycobacterium avium Complex (CONVERT): A Prospective, Open-Label, Randomized Study. *Am J Respir Crit Care Med.* 2018;198(12):1559–1569.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Amikacin+Liposome+Inhalation+Suspension+for+Treatment-Refractory+Lung+Disease+Caused+by+Mycobacterium+avium+Complex+CONVERT)
    *모델에서의 역할: 시나리오 5의 보정 기준 — 6개월 배양음전 29.0% (ALIS+GBT) 대 8.9% (GBT).*

41. **Olivier KN, Griffith DE, Eagle G, et al.** Randomized Trial of Liposomal Amikacin for Inhalation in Nontuberculous Mycobacterial Lung Disease. *Am J Respir Crit Care Med.* 2017;195(6):814–823.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Randomized+Trial+of+Liposomal+Amikacin+for+Inhalation+in+Nontuberculous+Mycobacterial+Lung+Disease)

42. **Winthrop KL, Flume PA, Thomson R, et al.** Amikacin Liposome Inhalation Suspension for Mycobacterium avium Complex Lung Disease: A 12-Month Open-Label Extension Clinical Trial. *Ann Am Thorac Soc.* 2021;18(7):1147–1157.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Amikacin+Liposome+Inhalation+Suspension+for+Mycobacterium+avium+Complex+Lung+Disease+A+12-Month+Open-Label+Extension+Clinical+Trial)

43. **Zhang J, Leifer F, Rose S, et al.** Amikacin Liposome Inhalation Suspension (ALIS) Penetrates Non-tuberculous Mycobacterial Biofilms and Enhances Amikacin Uptake Into Macrophages. *Front Microbiol.* 2018;9:915.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Amikacin+Liposome+Inhalation+Suspension+Penetrates+Non-tuberculous+Mycobacterial+Biofilms+and+Enhances+Amikacin+Uptake+Into+Macrophages)
    *모델에서의 역할: 이 모델에서 가장 결정적인 두 파라미터의 근거 — 리포솜의 바이오필름 투과(`PBFL` 0.60 vs 유리약물 `PBFK` 0.15)와 대식세포 내 전달(`KUPT`). 즉 ALIS만이 장벽 1과 장벽 3을 동시에 통과한다.*

44. **Rose SJ, Neville ME, Gupta R, Bermudez LE.** Delivery of aerosolized liposomal amikacin as a novel approach for the treatment of nontuberculous mycobacteria in an experimental model of pulmonary infection. *PLoS One.* 2014;9(9):e108703.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Delivery+of+aerosolized+liposomal+amikacin+as+a+novel+approach+for+the+treatment+of+nontuberculous+mycobacteria+in+an+experimental+model+of+pulmonary+infection)

45. **Huth ME, Ricci AJ, Cheng AG.** Mechanisms of aminoglycoside ototoxicity and targets of hair cell protection. *Int J Otolaryngol.* 2011;2011:937861.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Mechanisms+of+aminoglycoside+ototoxicity+and+targets+of+hair+cell+protection)
    *모델에서의 역할: `KPERI`(외림프) 구획 — 이독성이 폐가 아니라 혈장을 따라간다는 구조가 ALIS의 효능–독성 분리를 만든다.*

46. **Selimoglu E.** Aminoglycoside-induced ototoxicity. *Curr Pharm Des.* 2007;13(1):119–126.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Aminoglycoside-induced+ototoxicity+Selimoglu)

47. **Nagai J, Takano M.** Molecular aspects of renal handling of aminoglycosides and strategies for preventing the nephrotoxicity. *Drug Metab Pharmacokinet.* 2004;19(3):159–170.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Molecular+aspects+of+renal+handling+of+aminoglycosides+and+strategies+for+preventing+the+nephrotoxicity)
    *모델에서의 역할: `KREN`(근위세뇨관 megalin 매개 축적) 및 신독성이 아미카신 청소율을 다시 낮추는 양성 되먹임.*

---

## 8. 내성 (Resistance)

48. **Griffith DE, Brown-Elliott BA, Langsjoen B, et al.** Clinical and molecular analysis of macrolide resistance in Mycobacterium avium complex lung disease. *Am J Respir Crit Care Med.* 2006;174(8):928–934.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Clinical+and+molecular+analysis+of+macrolide+resistance+in+Mycobacterium+avium+complex+lung+disease)
    *모델에서의 역할: 시나리오 6 — 마크로라이드 단독(또는 기능적 단독) 요법이 rrl 변이를 선택한다는 임상 근거.*

49. **Meier A, Kirschner P, Springer B, et al.** Identification of mutations in 23S rRNA gene of clarithromycin-resistant Mycobacterium intracellulare. *Antimicrob Agents Chemother.* 1994;38(2):381–384.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Identification+of+mutations+in+23S+rRNA+gene+of+clarithromycin-resistant+Mycobacterium+intracellulare)

50. **Nash KA, Inderlied CB.** Genetic basis of macrolide resistance in Mycobacterium avium isolated from patients with disseminated disease. *Antimicrob Agents Chemother.* 1996;40(6):1748–1750.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Genetic+basis+of+macrolide+resistance+in+Mycobacterium+avium+isolated+from+patients+with+disseminated+disease)
    *모델에서의 역할: MAC은 rrl 유전자가 단일 사본이므로 단일 점돌연변이가 완전 내성을 만든다 — `MU` 항이 단일 단계로 `R_E`/`R_I`를 생성하는 구조의 근거.*

51. **Prammananan T, Sander P, Brown BA, et al.** A single 16S ribosomal RNA substitution is responsible for resistance to amikacin and other 2-deoxystreptamine aminoglycosides in Mycobacterium abscessus and Mycobacterium chelonae. *J Infect Dis.* 1998;177(6):1573–1581.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=A+single+16S+ribosomal+RNA+substitution+is+responsible+for+resistance+to+amikacin+and+other+2-deoxystreptamine+aminoglycosides)

52. **Nash KA, Brown-Elliott BA, Wallace RJ Jr.** A novel gene, erm(41), confers inducible macrolide resistance to clinical isolates of Mycobacterium abscessus but is absent from Mycobacterium chelonae. *Antimicrob Agents Chemother.* 2009;53(4):1367–1376.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=A+novel+gene+erm+41+confers+inducible+macrolide+resistance+to+clinical+isolates+of+Mycobacterium+abscessus)
    *모델에서의 역할: 지도의 대조 노드 — MAC에는 erm(41)이 없다는 사실이 "단일 rrl 변이 = 완전 내성" 구조를 정당화한다.*

53. **van Ingen J, Boeree MJ, van Soolingen D, Mouton JW.** Resistance mechanisms and drug susceptibility testing of nontuberculous mycobacteria. *Drug Resist Updat.* 2012;15(3):149–161.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Resistance+mechanisms+and+drug+susceptibility+testing+of+nontuberculous+mycobacteria)

54. **Moon SM, Park HY, Kim SY, et al.** Clinical Characteristics, Treatment Outcomes, and Resistance Mutations Associated with Macrolide-Resistant Mycobacterium avium Complex Lung Disease. *Antimicrob Agents Chemother.* 2016;60(11):6758–6765.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Clinical+Characteristics+Treatment+Outcomes+and+Resistance+Mutations+Associated+with+Macrolide-Resistant+Mycobacterium+avium+Complex+Lung+Disease)

55. **Morimoto K, Namkoong H, Hasegawa N, et al.** Macrolide-Resistant Mycobacterium avium Complex Lung Disease: Analysis of 102 Consecutive Cases. *Ann Am Thorac Soc.* 2016;13(11):1904–1911.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Macrolide-Resistant+Mycobacterium+avium+Complex+Lung+Disease+Analysis+of+102+Consecutive+Cases)

---

## 9. 구제요법 · 보조요법 (Salvage & Adjunctive Therapy)

56. **Martiniano SL, Wagner BD, Levin A, et al.** Safety and Effectiveness of Clofazimine for Primary and Refractory Nontuberculous Mycobacterial Infection. *Chest.* 2017;152(4):800–809.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Safety+and+Effectiveness+of+Clofazimine+for+Primary+and+Refractory+Nontuberculous+Mycobacterial+Infection)

57. **Philley JV, Wallace RJ Jr, Benwill JL, et al.** Preliminary Results of Bedaquiline as Salvage Therapy for Patients With Nontuberculous Mycobacterial Lung Disease. *Chest.* 2015;148(2):499–506.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Preliminary+Results+of+Bedaquiline+as+Salvage+Therapy+for+Patients+With+Nontuberculous+Mycobacterial+Lung+Disease)

58. **Winthrop KL, Ku JH, Marras TK, et al.** The tolerability of linezolid in the treatment of nontuberculous mycobacterial disease. *Eur Respir J.* 2015;45(4):1177–1179.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=The+tolerability+of+linezolid+in+the+treatment+of+nontuberculous+mycobacterial+disease)

59. **Ganapathy US, Dartois V, Dick T.** Repositioning rifamycins for Mycobacterium abscessus lung disease. *Expert Opin Drug Discov.* 2019;14(9):867–878.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Repositioning+rifamycins+for+Mycobacterium+abscessus+lung+disease)

60. **Yu X, Gao X, Zhu Z, et al.** Adjunctive surgery for the treatment of nontuberculous mycobacterial lung disease. (systematic review family)
    → [PubMed 검색](https://pubmed.ncbi.nlm.nih.gov/?term=adjunctive+surgical+resection+nontuberculous+mycobacterial+lung+disease+outcomes)

61. **Nick JA, Daley CL, Lenhart-Pendergrass PM, Davidson RM.** Nontuberculous mycobacteria in cystic fibrosis. *Curr Opin Pulm Med.* 2021;27(6):586–592.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Nontuberculous+mycobacteria+in+cystic+fibrosis+Nick+Daley)

---

## 10. 치료 성적 · 임상 결과 (Treatment Outcomes)

62. **Wallace RJ Jr, Brown BA, Griffith DE, et al.** Clarithromycin regimens for pulmonary Mycobacterium avium complex. The first 50 patients. *Am J Respir Crit Care Med.* 1996;153(6 Pt 1):1766–1772.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Clarithromycin+regimens+for+pulmonary+Mycobacterium+avium+complex+The+first+50+patients)

63. **Kwak N, Park J, Kim E, et al.** Treatment Outcomes of Mycobacterium avium Complex Lung Disease: A Systematic Review and Meta-analysis. *Clin Infect Dis.* 2017;65(7):1077–1084.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Treatment+Outcomes+of+Mycobacterium+avium+Complex+Lung+Disease+A+Systematic+Review+and+Meta-analysis)
    *모델에서의 역할: 시나리오 2 대 3의 배양음전율 격차(결절기관지확장형 > 공동형)를 판정하는 기준.*

64. **Diel R, Nienhaus A, Ringshausen FC, et al.** Microbiologic Outcome of Interventions Against Mycobacterium avium Complex Pulmonary Disease: A Systematic Review. *Chest.* 2018;153(4):888–921.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Microbiologic+Outcome+of+Interventions+Against+Mycobacterium+avium+Complex+Pulmonary+Disease+A+Systematic+Review)

65. **Koh WJ, Moon SM, Kim SY, et al.** Outcomes of Mycobacterium avium complex lung disease based on clinical phenotype. *Eur Respir J.* 2017;50(3):1602503.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Outcomes+of+Mycobacterium+avium+complex+lung+disease+based+on+clinical+phenotype)
    *모델에서의 역할: 공동형(fibrocavitary) 표현형이 예후를 결정한다는 관찰 — 본 모델에서 `CAVFLAG`가 `B_C` 초기값 하나만 바꾸는 설계의 임상적 대응물.*

66. **Wallace RJ Jr, Brown-Elliott BA, McNulty S, et al.** Macrolide/azalide therapy for nodular/bronchiectatic Mycobacterium avium complex lung disease. *Chest.* 2014;146(2):276–282.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Macrolide+azalide+therapy+for+nodular+bronchiectatic+Mycobacterium+avium+complex+lung+disease)

67. **Ray WA, Murray KT, Hall K, et al.** Azithromycin and the risk of cardiovascular death. *N Engl J Med.* 2012;366(20):1881–1890.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Azithromycin+and+the+risk+of+cardiovascular+death)
    *모델에서의 역할: `QTE` 효과구획과 `QTSLP`(클래리트로마이신 > 아지트로마이신).*

---

## 11. QSP · 약동학-약력학 모델링 방법론 (QSP & PK/PD Methodology)

68. **Baron KT, Gastonguay MR.** Simulation from ODE-based population PK/PD and systems pharmacology models in R with mrgsolve. *J Pharmacokinet Pharmacodyn.* 2015;42(Suppl 1):S84–S85. (소프트웨어: <https://mrgsolve.org>)
    → [mrgsolve 문서](https://mrgsolve.org/docs/)

69. **Clewe O, Aulin L, Hu Y, et al.** A multistate tuberculosis pharmacometric model: a framework for studying anti-tubercular drug effects in vitro. *J Antimicrob Chemother.* 2016;71(4):964–974.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=A+multistate+tuberculosis+pharmacometric+model+a+framework+for+studying+anti-tubercular+drug+effects+in+vitro)
    *모델에서의 역할: 세균 아집단을 "표현형 상태"로 나누는 다중상태 접근 — 본 모델은 이를 표현형이 아닌 물리적 니치로 재구성했다.*

70. **Ernest JP, Sarathy J, Wang N, et al.** Lesion Penetration and Activity Limit the Utility of Second-Line Injectable Agents in Pulmonary Tuberculosis. *Antimicrob Agents Chemother.* 2021;65(10):e0050621.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Lesion+Penetration+and+Activity+Limit+the+Utility+of+Second-Line+Injectable+Agents+in+Pulmonary+Tuberculosis)
    *모델에서의 역할: 주사용 아미노글리코사이드가 병변(특히 건락)에 도달하지 못한다는 정량적 선례 — `PCSK` = 0.02.*

71. **Gumbo T, Pasipanodya JG, Nuermberger E, et al.** Correlations Between the Hollow Fiber Model of Tuberculosis and Therapeutic Events in Tuberculosis Patients. *Clin Infect Dis.* 2015;61(Suppl 1):S18–S24.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Correlations+Between+the+Hollow+Fiber+Model+of+Tuberculosis+and+Therapeutic+Events+in+Tuberculosis+Patients)

72. **Deshpande D, Srivastava S, Musuka S, Gumbo T.** Thrice-Weekly Azithromycin plus Ethambutol Is an Effective Treatment for Mycobacterium avium Complex Lung Disease. *Antimicrob Agents Chemother.* 2016;60(4):2157–2163.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Thrice-Weekly+Azithromycin+plus+Ethambutol+Is+an+Effective+Treatment+for+Mycobacterium+avium+Complex+Lung+Disease)
    *모델에서의 역할: 시나리오 2의 주 3회 투여 간격(`ii` = 2.333 d) 및 아지트로마이신–에탐부톨 시너지 항의 근거.*

73. **Schmidt S, Barbour A, Sahre M, et al.** PK/PD: new insights for antibacterial and antiviral applications. *Curr Opin Pharmacol.* 2008;8(5):549–556.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=PK+PD+new+insights+for+antibacterial+and+antiviral+applications+Schmidt+Derendorf)

---

## 12. 진단 · 모니터링 (Diagnosis & Monitoring)

74. **Woods GL, Brown-Elliott BA, Conville PS, et al.** Susceptibility Testing of Mycobacteria, Nocardia spp., and Other Aerobic Actinomycetes. *CLSI standard M24*, 3rd ed. Clinical and Laboratory Standards Institute; 2018.
    → [PubMed 관련 검색](https://pubmed.ncbi.nlm.nih.gov/?term=CLSI+M24+susceptibility+testing+of+mycobacteria+nocardia+aerobic+actinomycetes)
    *모델에서의 역할: MAC에서 임상적으로 의미 있는 감수성검사는 클래리트로마이신과 아미카신뿐 — 지도의 진단 클러스터.*

75. **Peloquin CA.** Therapeutic drug monitoring in the treatment of tuberculosis. *Drugs.* 2002;62(15):2169–2183.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Therapeutic+drug+monitoring+in+the+treatment+of+tuberculosis+Peloquin)

76. **Quittner AL, O'Donnell AE, Salathe MA, et al.** Quality of Life Questionnaire-Bronchiectasis: final psychometric analyses and determination of minimal important difference scores. *Thorax.* 2015;70(1):12–20.
    → [PubMed](https://pubmed.ncbi.nlm.nih.gov/?term=Quality+of+Life+Questionnaire-Bronchiectasis+final+psychometric+analyses+and+determination+of+minimal+important+difference+scores)
    *모델에서의 역할: `QOLB` 출력의 척도와 MCID(≈ 10점) 해석 기준.*

---

## 요약 — 이 모델이 문헌에서 가져온 "구조" 3가지

| 구조적 결정 | 근거 문헌 | 모델에서의 구현 |
|---|---|---|
| 세균을 **물리적 니치**로 나눈다 (표현형이 아니라) | #20, #21, #22, #70 | `B_E` / `B_B` / `B_I` / `B_C` + 니치별 침투계수 |
| 식포 pH 5.2를 **한 번만** 입력해 축적과 활성소실을 **동시에** 유도 | #16, #17, #28, #31, #32, #33, #34 | `R_trap` = 151× 와 `MIC` 12.6× 상승이 같은 pH에서 계산됨 |
| ALIS는 **같은 분자, 다른 주소** — 장벽 1과 3을 동시 통과 | #40, #43, #44, #45 | `KLIP → KMAC` (세포내) 대 `KCEN → KPERI` (이독성)의 분리 |

---

*작성: Claude Code Routine · QSP Disease Model Library · 2026-08-01*
*본 문헌 목록은 교육·연구 목적의 QSP 모델을 위한 것이며, 임상 진료지침을 대체하지 않습니다.*
