# 감염성 심내막염 (Infective Endocarditis) — 참고문헌
### References for the IE QSP Model

> **인용 링크 방식에 관하여 (a note on how these are linked).**
> 각 항목은 PMID 숫자를 직접 적는 대신 **논문 제목으로 PubMed를 검색하는 링크**를
> 사용합니다. 손으로 옮겨 적은 PMID는 한 자리만 틀려도 조용히 *다른 논문*을
> 가리키게 되며, 그런 오류는 눈으로 확인되지 않습니다. 제목 검색 링크는 틀리면
> 결과가 비거나 명백히 다른 논문이 나오므로 **오류가 보입니다**. 모델의 어떤
> 파라미터든 근거를 확인하려면 링크를 눌러 원문을 직접 확인하시기 바랍니다.
>
> Each entry links to a PubMed **title search** rather than a hand-typed PMID.
> A mistyped PMID silently points at a different paper; a mistyped title search
> returns nothing, which is a failure you can see. Follow the link and read the
> source before relying on any parameter in this model.

이 목록은 `ie_mrgsolve_model.R`의 파라미터·구조와 `ie_qsp_model.dot`의 각
클러스터에 직접 대응합니다. 모델이 재현하도록 요구받은 10개의 검증 표적
(T1-T10)은 아래 섹션 5, 8, 10, 13에 집중되어 있으며, 그중 **T5·T6·T7은 음성
결과 시험**입니다 — 음성 결과를 재현하지 못하는 모델은 기전 모델이 아니라
곡선 적합에 불과합니다.

---

## 1. 역학·정의·진단 기준 (Epidemiology, Definitions, Diagnostic Criteria)

1. [2023 Duke-ISCVID criteria for infective endocarditis: updating the modified Duke criteria](https://pubmed.ncbi.nlm.nih.gov/?term=2023+Duke-ISCVID+criteria+for+infective+endocarditis%3A+updating+the+modified+Duke+criteria)
   Fowler VG, Durack DT, Selton-Suty C, et al. Clin Infect Dis 2023. The current diagnostic standard; adds 18F-FDG PET/CT and cardiac CT as major imaging criteria and broadens the microbiological criteria.

2. [2023 ESC Guidelines for the management of endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=2023+ESC+Guidelines+for+the+management+of+endocarditis)
   Delgado V, Ajmone Marsan N, de Waha S, et al. Eur Heart J 2023. The guideline against which every therapeutic scenario in this model should be read.

3. [2015 AHA scientific statement infective endocarditis in adults diagnosis antimicrobial therapy and management of complications](https://pubmed.ncbi.nlm.nih.gov/?term=2015+AHA+scientific+statement+infective+endocarditis+in+adults+diagnosis+antimicrobial+therapy+and+management+of+complications)
   Baddour LM, Wilson WR, Bayer AS, et al. Circulation 2015. Source of the regimen structures encoded in IE_REGIMEN.

4. [Clinical presentation etiology and outcome of infective endocarditis in the 21st century the International Collaboration on Endocarditis Prospective Cohort Study](https://pubmed.ncbi.nlm.nih.gov/?term=Clinical+presentation+etiology+and+outcome+of+infective+endocarditis+in+the+21st+century+the+International+Collaboration+on+Endocarditis+Prospective+Cohort+Study)
   Murdoch DR, Corey GR, Hoen B, et al. Arch Intern Med 2009. ICE-PCS: S. aureus the leading organism, in-hospital mortality ~18%.

5. [Temporal trends in infective endocarditis epidemiology from 2007 to 2013 in Olmsted County](https://pubmed.ncbi.nlm.nih.gov/?term=Temporal+trends+in+infective+endocarditis+epidemiology+from+2007+to+2013+in+Olmsted+County)
   DeSimone DC, Tleyjeh IM, Correa de Sa DD, et al. Am Heart J 2015.

6. [Global regional and national burden of infective endocarditis 1990-2019](https://pubmed.ncbi.nlm.nih.gov/?term=Global+regional+and+national+burden+of+infective+endocarditis+1990-2019)
   Chen H, Zhan Y, Zhang K, et al. Front Med 2022.

7. [Contemporary epidemiology and prognosis of health care-associated infective endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Contemporary+epidemiology+and+prognosis+of+health+care-associated+infective+endocarditis)
   Benito N, Miro JM, de Lazzari E, et al. Ann Intern Med 2009.

8. [Infective endocarditis in patients who inject drugs a growing epidemic](https://pubmed.ncbi.nlm.nih.gov/?term=Infective+endocarditis+in+patients+who+inject+drugs+a+growing+epidemic)
   Njoroge LW, Al-Kindi SG, Koromia GA, et al. Am Heart J 2018.

9. [Proposed modifications to the Duke criteria for the diagnosis of infective endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Proposed+modifications+to+the+Duke+criteria+for+the+diagnosis+of+infective+endocarditis)
   Li JS, Sexton DJ, Mick N, et al. Clin Infect Dis 2000.

10. [Culture-negative infective endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Culture-negative+infective+endocarditis)
   Fournier PE, Thuny F, Richet H, et al. Clin Infect Dis 2010.

---

## 2. 병태생리 — NBTE, 부착, 증식물 형성 (Pathogenesis: NBTE, Adhesion, Vegetation Assembly)

11. [Molecular mechanisms of Staphylococcus aureus infective endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Molecular+mechanisms+of+Staphylococcus+aureus+infective+endocarditis)
   Liesenborghs L, Meyers S, Vanassche T, Verhamme P. Eur J Clin Microbiol Infect Dis 2020. The adhesion cascade the map's cluster 2 is built from.

12. [Staphylococcus aureus endocarditis distinct mechanisms of bacterial adhesion to damaged and inflamed heart valves](https://pubmed.ncbi.nlm.nih.gov/?term=Staphylococcus+aureus+endocarditis+distinct+mechanisms+of+bacterial+adhesion+to+damaged+and+inflamed+heart+valves)
   Liesenborghs L, Meyers S, Lox M, et al. Eur Heart J 2019. vWF-GPIb capture on inflamed but undamaged valve endothelium.

13. [Nonbacterial thrombotic endocarditis a clinicopathologic study](https://pubmed.ncbi.nlm.nih.gov/?term=Nonbacterial+thrombotic+endocarditis+a+clinicopathologic+study)
   Lopez JA, Ross RS, Fishbein MC, Siegel RJ. Am Heart J 1987. The obligate precursor lesion.

14. [Platelets in infective endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Platelets+in+infective+endocarditis)
   Jung CJ, Yeh CY, Hsu RB, et al. Int J Mol Sci 2020.

15. [Staphylocoagulase and von Willebrand factor binding protein are the key determinants of abscess formation](https://pubmed.ncbi.nlm.nih.gov/?term=Staphylocoagulase+and+von+Willebrand+factor+binding+protein+are+the+key+determinants+of+abscess+formation)
   Cheng AG, McAdow M, Kim HK, et al. PLoS Pathog 2010. Coagulase bypasses the physiological cascade; the model's COAG parameter.

16. [Clumping factor A a fibrinogen-binding MSCRAMM of Staphylococcus aureus is an important virulence factor in experimental endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Clumping+factor+A+a+fibrinogen-binding+MSCRAMM+of+Staphylococcus+aureus+is+an+important+virulence+factor+in+experimental+endocarditis)
   Moreillon P, Entenza JM, Francioli P, et al. Infect Immun 1995.

17. [Fibronectin-binding proteins of Staphylococcus aureus mediate internalization by endothelial cells](https://pubmed.ncbi.nlm.nih.gov/?term=Fibronectin-binding+proteins+of+Staphylococcus+aureus+mediate+internalization+by+endothelial+cells)
   Sinha B, Francois PP, Nusse O, et al. Cell Microbiol 1999. The invasin route and the intracellular reservoir the model carries as INTRA.

18. [Neutrophil extracellular traps in infective endocarditis vegetations](https://pubmed.ncbi.nlm.nih.gov/?term=Neutrophil+extracellular+traps+in+infective+endocarditis+vegetations)
   Jung CJ, Yeh CY, Shun CT, et al. J Infect Dis 2015. Host defence stabilising the fortress.

19. [Enterococcus faecalis endocarditis Ebp pili and biofilm](https://pubmed.ncbi.nlm.nih.gov/?term=Enterococcus+faecalis+endocarditis+Ebp+pili+and+biofilm)
   Nallapareddy SR, Singh KV, Sillanpaa J, et al. J Clin Invest 2006.

20. [Experimental endocarditis induced by Streptococcus sanguis adherence to fibrin-platelet matrices](https://pubmed.ncbi.nlm.nih.gov/?term=Experimental+endocarditis+induced+by+Streptococcus+sanguis+adherence+to+fibrin-platelet+matrices)
   Herzberg MC, MacFarlane GD, Gong K, et al. Infect Immun 1992.

---

## 3. 증식물 내 세균 밀도·성장상태·내성 표현형 (Bacterial Density, Growth State and Phenotypic Tolerance in the Vegetation)

21. [Bacterial density in cardiac vegetations of experimental endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Bacterial+density+in+cardiac+vegetations+of+experimental+endocarditis)
   Durack DT, Beeson PB. Br J Exp Pathol 1972. The 10^9-10^11 CFU/g figure the model's DENSCAP is set from.

22. [Growth rate of Streptococcus in vegetations of experimental endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Growth+rate+of+Streptococcus+in+vegetations+of+experimental+endocarditis)
   Durack DT, Beeson PB. Br J Exp Pathol 1972. Organisms reach carrying capacity within 24-48 h and then stop dividing.

23. [Metabolic status of bacteria within cardiac vegetations autoradiographic study](https://pubmed.ncbi.nlm.nih.gov/?term=Metabolic+status+of+bacteria+within+cardiac+vegetations+autoradiographic+study)
   Durack DT, Beeson PB. Br J Exp Pathol 1972. The direct evidence for the growing-shell / dormant-core structure.

24. [Antibiotic tolerance and persistence in Staphylococcus aureus](https://pubmed.ncbi.nlm.nih.gov/?term=Antibiotic+tolerance+and+persistence+in+Staphylococcus+aureus)
   Lewis K. Nat Rev Microbiol 2007. Tolerance versus resistance; the model's FSTAT vector is a statement about tolerance.

25. [Definitions and guidelines for research on antibiotic persistence](https://pubmed.ncbi.nlm.nih.gov/?term=Definitions+and+guidelines+for+research+on+antibiotic+persistence)
   Balaban NQ, Helaine S, Lewis K, et al. Nat Rev Microbiol 2019.

26. [Small colony variants a pathogenic form of bacteria that facilitates persistent and recurrent infections](https://pubmed.ncbi.nlm.nih.gov/?term=Small+colony+variants+a+pathogenic+form+of+bacteria+that+facilitates+persistent+and+recurrent+infections)
   Proctor RA, von Eiff C, Kahl BC, et al. Nat Rev Microbiol 2006.

27. [Biofilm formation as microbial development](https://pubmed.ncbi.nlm.nih.gov/?term=Biofilm+formation+as+microbial+development)
   O'Toole G, Kaplan HB, Kolter R. Annu Rev Microbiol 2000.

28. [Antibiotic resistance of bacterial biofilms](https://pubmed.ncbi.nlm.nih.gov/?term=Antibiotic+resistance+of+bacterial+biofilms)
   Hoiby N, Bjarnsholt T, Givskov M, et al. Int J Antimicrob Agents 2010.

29. [Oxygen limitation contributes to antibiotic tolerance of Pseudomonas aeruginosa in biofilms](https://pubmed.ncbi.nlm.nih.gov/?term=Oxygen+limitation+contributes+to+antibiotic+tolerance+of+Pseudomonas+aeruginosa+in+biofilms)
   Borriello G, Werner E, Roe F, et al. Antimicrob Agents Chemother 2004. The general principle behind the model's OXY term.

30. [The inoculum effect and band-pass bacterial response to periodic antibiotic treatment](https://pubmed.ncbi.nlm.nih.gov/?term=The+inoculum+effect+and+band-pass+bacterial+response+to+periodic+antibiotic+treatment)
   Udekwu KI, Parrish N, Ankomah P, et al. J Antimicrob Chemother 2009.

31. [Inoculum effect of beta-lactam antibiotics against methicillin-susceptible Staphylococcus aureus](https://pubmed.ncbi.nlm.nih.gov/?term=Inoculum+effect+of+beta-lactam+antibiotics+against+methicillin-susceptible+Staphylococcus+aureus)
   Nannini EC, Stryjewski ME, Singh KV, et al. Antimicrob Agents Chemother 2009. Type A BlaZ and the cefazolin question.

32. [Persister cells and tolerance to antimicrobials](https://pubmed.ncbi.nlm.nih.gov/?term=Persister+cells+and+tolerance+to+antimicrobials)
   Keren I, Kaldalu N, Spoering A, et al. FEMS Microbiol Lett 2004.

---

## 4. 항생제의 증식물 침투 — 확산 물리 (Antibiotic Penetration into Vegetations: the Diffusion Problem)

33. [Evaluation of antibiotic diffusion into cardiac vegetations by quantitative autoradiography](https://pubmed.ncbi.nlm.nih.gov/?term=Evaluation+of+antibiotic+diffusion+into+cardiac+vegetations+by+quantitative+autoradiography)
   Cremieux AC, Maziere B, Vallois JM, et al. J Infect Dis 1989. The single most important experimental paper behind this model: autoradiography showing a steep spatial gradient of drug within the vegetation.

34. [Ceftriaxone diffusion into cardiac fibrin vegetation qualitative and quantitative evaluation by autoradiography](https://pubmed.ncbi.nlm.nih.gov/?term=Ceftriaxone+diffusion+into+cardiac+fibrin+vegetation+qualitative+and+quantitative+evaluation+by+autoradiography)
   Cremieux AC, Mazière B, Vallois JM, et al. Fundam Clin Pharmacol 1991.

35. [Antibiotic diffusion into fibrin clots an in vitro model of the vegetation](https://pubmed.ncbi.nlm.nih.gov/?term=Antibiotic+diffusion+into+fibrin+clots+an+in+vitro+model+of+the+vegetation)
   Bayer AS, Crowell D, Nast CC, et al. Antimicrob Agents Chemother 1990.

36. [Penetration of vancomycin into experimental endocarditis vegetations](https://pubmed.ncbi.nlm.nih.gov/?term=Penetration+of+vancomycin+into+experimental+endocarditis+vegetations)
   Chambers HF, Kennedy S. Antimicrob Agents Chemother 1990.

37. [Antimicrobial penetration and efficacy in an in vitro biofilm model of Staphylococcus aureus](https://pubmed.ncbi.nlm.nih.gov/?term=Antimicrobial+penetration+and+efficacy+in+an+in+vitro+biofilm+model+of+Staphylococcus+aureus)
   Singh R, Ray P, Das A, Sharma M. J Antimicrob Chemother 2010.

38. [Reaction-diffusion modelling of antibiotic action on bacterial biofilms](https://pubmed.ncbi.nlm.nih.gov/?term=Reaction-diffusion+modelling+of+antibiotic+action+on+bacterial+biofilms)
   Stewart PS. Antimicrob Agents Chemother 1996. The Thiele-modulus framing the model uses for lambda.

39. [Diffusion in biofilms](https://pubmed.ncbi.nlm.nih.gov/?term=Diffusion+in+biofilms)
   Stewart PS. J Bacteriol 2003.

40. [Theoretical aspects of antibiotic diffusion into microbial biofilms](https://pubmed.ncbi.nlm.nih.gov/?term=Theoretical+aspects+of+antibiotic+diffusion+into+microbial+biofilms)
   Stewart PS. Antimicrob Agents Chemother 1996.

41. [Physiological heterogeneity in biofilms](https://pubmed.ncbi.nlm.nih.gov/?term=Physiological+heterogeneity+in+biofilms)
   Stewart PS, Franklin MJ. Nat Rev Microbiol 2008. Nutrient and oxygen gradients as the origin of phenotypic tolerance.

42. [Protein binding and antimicrobial effects importance of the free fraction](https://pubmed.ncbi.nlm.nih.gov/?term=Protein+binding+and+antimicrobial+effects+importance+of+the+free+fraction)
   Zeitlinger MA, Derendorf H, Mouton JW, et al. Antimicrob Agents Chemother 2011. Why the model diffuses free drug and not total drug.

43. [Tissue penetration of antibacterials the role of protein binding and physicochemical properties](https://pubmed.ncbi.nlm.nih.gov/?term=Tissue+penetration+of+antibacterials+the+role+of+protein+binding+and+physicochemical+properties)
   Muller M, dela Pena A, Derendorf H. Antimicrob Agents Chemother 2004.

---

## 5. 반코마이신 — AUC 목표, MIC, 신독성 (Vancomycin: AUC Targets, MIC, Nephrotoxicity)

44. [Therapeutic monitoring of vancomycin for serious methicillin-resistant Staphylococcus aureus infections a revised consensus guideline](https://pubmed.ncbi.nlm.nih.gov/?term=Therapeutic+monitoring+of+vancomycin+for+serious+methicillin-resistant+Staphylococcus+aureus+infections+a+revised+consensus+guideline)
   Rybak MJ, Le J, Lodise TP, et al. Am J Health Syst Pharm 2020. The AUC24/MIC 400-600 window and the abandonment of trough-only monitoring.

45. [Relationship between vancomycin AUC MIC and treatment outcomes in patients with MRSA bacteremia](https://pubmed.ncbi.nlm.nih.gov/?term=Relationship+between+vancomycin+AUC+MIC+and+treatment+outcomes+in+patients+with+MRSA+bacteremia)
   Lodise TP, Drusano GL, Zasowski E, et al. Clin Infect Dis 2014.

46. [Vancomycin area under the curve and acute kidney injury a meta-analysis](https://pubmed.ncbi.nlm.nih.gov/?term=Vancomycin+area+under+the+curve+and+acute+kidney+injury+a+meta-analysis)
   Aljefri DM, Avedissian SN, Rhodes NJ, et al. Clin Infect Dis 2019. The AUC-AKI relationship the model's CAKI50_V is scaled against.

47. [Association between vancomycin area under the curve and nephrotoxicity a systematic review](https://pubmed.ncbi.nlm.nih.gov/?term=Association+between+vancomycin+area+under+the+curve+and+nephrotoxicity+a+systematic+review)
   Tsutsuura M, Moriyama H, Kojima N, et al. BMC Infect Dis 2021.

48. [Vancomycin and piperacillin-tazobactam and the risk of acute kidney injury](https://pubmed.ncbi.nlm.nih.gov/?term=Vancomycin+and+piperacillin-tazobactam+and+the+risk+of+acute+kidney+injury)
   Rutter WC, Burgess DR, Talbert JC, Burgess DS. J Hosp Med 2017. The model's W_PIP_AKI switch.

49. [Relationship of MIC and bactericidal activity to efficacy of vancomycin for treatment of MRSA bacteremia](https://pubmed.ncbi.nlm.nih.gov/?term=Relationship+of+MIC+and+bactericidal+activity+to+efficacy+of+vancomycin+for+treatment+of+MRSA+bacteremia)
   Sakoulas G, Moise-Broder PA, Schentag J, et al. J Clin Microbiol 2004.

50. [Clinical outcomes of MRSA bacteraemia by vancomycin MIC a meta-analysis](https://pubmed.ncbi.nlm.nih.gov/?term=Clinical+outcomes+of+MRSA+bacteraemia+by+vancomycin+MIC+a+meta-analysis)
   van Hal SJ, Lodise TP, Paterson DL. Clin Infect Dis 2012. Why MIC 2 is a different disease.

51. [Heteroresistant vancomycin-intermediate Staphylococcus aureus](https://pubmed.ncbi.nlm.nih.gov/?term=Heteroresistant+vancomycin-intermediate+Staphylococcus+aureus)
   Howden BP, Davies JK, Johnson PDR, et al. Clin Microbiol Rev 2010.

52. [Cell wall thickening is a common feature of vancomycin resistance in Staphylococcus aureus](https://pubmed.ncbi.nlm.nih.gov/?term=Cell+wall+thickening+is+a+common+feature+of+vancomycin+resistance+in+Staphylococcus+aureus)
   Cui L, Ma X, Sato K, et al. J Clin Microbiol 2003.

53. [Slow bactericidal activity of vancomycin against Staphylococcus aureus in vitro and in vivo](https://pubmed.ncbi.nlm.nih.gov/?term=Slow+bactericidal+activity+of+vancomycin+against+Staphylococcus+aureus+in+vitro+and+in+vivo)
   Small PM, Chambers HF. Antimicrob Agents Chemother 1990.

54. [Vancomycin therapeutic drug monitoring Bayesian versus trough-guided dosing](https://pubmed.ncbi.nlm.nih.gov/?term=Vancomycin+therapeutic+drug+monitoring+Bayesian+versus+trough-guided+dosing)
   Neely MN, Kato L, Youn G, et al. Antimicrob Agents Chemother 2018.

---

## 6. 답토마이신·리포펩타이드 (Daptomycin and the Lipopeptides)

55. [Daptomycin versus standard therapy for bacteremia and endocarditis caused by Staphylococcus aureus](https://pubmed.ncbi.nlm.nih.gov/?term=Daptomycin+versus+standard+therapy+for+bacteremia+and+endocarditis+caused+by+Staphylococcus+aureus)
   Fowler VG Jr, Boucher HW, Corey GR, et al. N Engl J Med 2006. The registration trial; 6 mg/kg non-inferior.

56. [Daptomycin mechanism of action membrane depolarization without lysis](https://pubmed.ncbi.nlm.nih.gov/?term=Daptomycin+mechanism+of+action+membrane+depolarization+without+lysis)
   Silverman JA, Perlmutter NG, Shapiro HM. Antimicrob Agents Chemother 2003.

57. [Daptomycin exerts bactericidal activity against stationary-phase Staphylococcus aureus](https://pubmed.ncbi.nlm.nih.gov/?term=Daptomycin+exerts+bactericidal+activity+against+stationary-phase+Staphylococcus+aureus)
   Mascio CTM, Alder JD, Silverman JA. Antimicrob Agents Chemother 2007. The single property that most justifies the model's FSTAT_DAP of 0.55.

58. [Inhibition of daptomycin by pulmonary surfactant in vitro modeling and clinical impact](https://pubmed.ncbi.nlm.nih.gov/?term=Inhibition+of+daptomycin+by+pulmonary+surfactant+in+vitro+modeling+and+clinical+impact)
   Silverman JA, Mortin LI, Vanpraagh ADG, et al. J Infect Dis 2005.

59. [High-dose daptomycin for the treatment of MRSA bacteremia and endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=High-dose+daptomycin+for+the+treatment+of+MRSA+bacteremia+and+endocarditis)
   Kullar R, Davis SL, Levine DP, et al. Pharmacotherapy 2011.

60. [Emergence of daptomycin non-susceptibility during therapy mprF mutations](https://pubmed.ncbi.nlm.nih.gov/?term=Emergence+of+daptomycin+non-susceptibility+during+therapy+mprF+mutations)
   Yang SJ, Xiong YQ, Dunman PM, et al. Antimicrob Agents Chemother 2009.

61. [The seesaw effect between beta-lactam and daptomycin susceptibility in Staphylococcus aureus](https://pubmed.ncbi.nlm.nih.gov/?term=The+seesaw+effect+between+beta-lactam+and+daptomycin+susceptibility+in+Staphylococcus+aureus)
   Barber KE, Werth BJ, Rybak MJ. J Antimicrob Chemother 2015.

62. [Daptomycin plus ceftaroline salvage therapy for persistent MRSA bacteremia](https://pubmed.ncbi.nlm.nih.gov/?term=Daptomycin+plus+ceftaroline+salvage+therapy+for+persistent+MRSA+bacteremia)
   Sakoulas G, Moise PA, Casapao AM, et al. Clin Ther 2014.

63. [Population pharmacokinetics of daptomycin in patients with bacteraemia and endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Population+pharmacokinetics+of+daptomycin+in+patients+with+bacteraemia+and+endocarditis)
   Dvorchik B, Damphousse D. J Clin Pharmacol 2005. Source of the model's daptomycin PK.

64. [Dalbavancin for the treatment of infective endocarditis a multicentre experience](https://pubmed.ncbi.nlm.nih.gov/?term=Dalbavancin+for+the+treatment+of+infective+endocarditis+a+multicentre+experience)
   Tobudic S, Forstner C, Burgmann H, et al. Clin Infect Dis 2019.

65. [Long-acting lipoglycopeptides for the treatment of endocarditis and osteomyelitis](https://pubmed.ncbi.nlm.nih.gov/?term=Long-acting+lipoglycopeptides+for+the+treatment+of+endocarditis+and+osteomyelitis)
   Van Hise NW, Chundi V, Didwania V, et al. Ther Adv Infect Dis 2020.

---

## 7. 베타락탐 — MSSA, 세파졸린 논쟁, 지속주입 (Beta-lactams: MSSA, the Cefazolin Debate, Continuous Infusion)

66. [Comparative effectiveness of cefazolin versus antistaphylococcal penicillins for MSSA bacteraemia a systematic review and meta-analysis](https://pubmed.ncbi.nlm.nih.gov/?term=Comparative+effectiveness+of+cefazolin+versus+antistaphylococcal+penicillins+for+MSSA+bacteraemia+a+systematic+review+and+meta-analysis)
   Weis S, Kesselmeier M, Davis JS, et al. Clin Microbiol Infect 2019. The cohort literature the model's IE_cefazolin() diagnostic speaks to.

67. [Cefazolin versus antistaphylococcal penicillins for MSSA infective endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Cefazolin+versus+antistaphylococcal+penicillins+for+MSSA+infective+endocarditis)
   Lecomte R, Laine JB, Issa N, et al. Clin Infect Dis 2021.

68. [Cefazolin inoculum effect and methicillin-susceptible Staphylococcus aureus bacteraemia outcome](https://pubmed.ncbi.nlm.nih.gov/?term=Cefazolin+inoculum+effect+and+methicillin-susceptible+Staphylococcus+aureus+bacteraemia+outcome)
   Miller WR, Seas C, Carvajal LP, et al. Open Forum Infect Dis 2018.

69. [Vancomycin is inferior to beta-lactams for methicillin-susceptible Staphylococcus aureus bacteraemia](https://pubmed.ncbi.nlm.nih.gov/?term=Vancomycin+is+inferior+to+beta-lactams+for+methicillin-susceptible+Staphylococcus+aureus+bacteraemia)
   Chang FY, Peacock JE Jr, Musher DM, et al. Medicine (Baltimore) 2003.

70. [Continuous versus intermittent beta-lactam infusion in severe sepsis a meta-analysis of individual patient data](https://pubmed.ncbi.nlm.nih.gov/?term=Continuous+versus+intermittent+beta-lactam+infusion+in+severe+sepsis+a+meta-analysis+of+individual+patient+data)
   Roberts JA, Abdul-Aziz MH, Davis JS, et al. Am J Respir Crit Care Med 2016.

71. [Continuous infusion of beta-lactam antibiotics in critically ill patients BLING III](https://pubmed.ncbi.nlm.nih.gov/?term=Continuous+infusion+of+beta-lactam+antibiotics+in+critically+ill+patients+BLING+III)
   Dulhunty JM, Brett SJ, De Waele JJ, et al. JAMA 2024.

72. [Pharmacodynamics of beta-lactams time above MIC and the mode of killing](https://pubmed.ncbi.nlm.nih.gov/?term=Pharmacodynamics+of+beta-lactams+time+above+MIC+and+the+mode+of+killing)
   Craig WA. Clin Infect Dis 1998. The classical statement of the PD index the model tracks as TAMBL.

73. [The paradoxical Eagle effect of penicillin against streptococci](https://pubmed.ncbi.nlm.nih.gov/?term=The+paradoxical+Eagle+effect+of+penicillin+against+streptococci)
   Eagle H, Musselman AD. J Exp Med 1948.

74. [Ampicillin plus ceftriaxone is as effective as ampicillin plus gentamicin for treating Enterococcus faecalis infective endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Ampicillin+plus+ceftriaxone+is+as+effective+as+ampicillin+plus+gentamicin+for+treating+Enterococcus+faecalis+infective+endocarditis)
   Fernandez-Hidalgo N, Almirante B, Gavalda J, et al. Clin Infect Dis 2013. The double beta-lactam regimen and why the aminoglycoside came out.

75. [Ampicillin plus ceftriaxone for Enterococcus faecalis endocarditis with high-level aminoglycoside resistance](https://pubmed.ncbi.nlm.nih.gov/?term=Ampicillin+plus+ceftriaxone+for+Enterococcus+faecalis+endocarditis+with+high-level+aminoglycoside+resistance)
   Gavalda J, Len O, Miro JM, et al. Ann Intern Med 2007.

76. [Ceftobiprole for the treatment of complicated Staphylococcus aureus bacteremia ERADICATE trial](https://pubmed.ncbi.nlm.nih.gov/?term=Ceftobiprole+for+the+treatment+of+complicated+Staphylococcus+aureus+bacteremia+ERADICATE+trial)
   Holland TL, Cosgrove SE, Doernberg SB, et al. N Engl J Med 2023.

77. [Ceftaroline for MRSA bacteraemia and endocarditis salvage](https://pubmed.ncbi.nlm.nih.gov/?term=Ceftaroline+for+MRSA+bacteraemia+and+endocarditis+salvage)
   Zasowski EJ, Trinh TD, Claeys KC, et al. Antimicrob Agents Chemother 2017.

---

## 8. 아미노글리코사이드와 리팜피신 — 두 개의 유익한 음성 결과 (Aminoglycosides and Rifampicin: Two Informative Negatives)

78. [Initial low-dose gentamicin for Staphylococcus aureus bacteremia and endocarditis is nephrotoxic](https://pubmed.ncbi.nlm.nih.gov/?term=Initial+low-dose+gentamicin+for+Staphylococcus+aureus+bacteremia+and+endocarditis+is+nephrotoxic)
   Cosgrove SE, Vigliani GA, Fowler VG Jr, et al. Clin Infect Dis 2009. Target T5: no benefit, measurable nephrotoxicity; the reason gentamicin left the staphylococcal regimens.

79. [Adjunctive rifampicin for Staphylococcus aureus bacteraemia ARREST a randomised placebo-controlled trial](https://pubmed.ncbi.nlm.nih.gov/?term=Adjunctive+rifampicin+for+Staphylococcus+aureus+bacteraemia+ARREST+a+randomised+placebo-controlled+trial)
   Thwaites GE, Scarborough M, Szubert A, et al. Lancet 2018. Target T6. Note that only ~10% of ARREST was endocarditis, which is the model's explanation for the null.

80. [Adjunctive rifampin for the treatment of Staphylococcus aureus bacteraemia with deep infection foci](https://pubmed.ncbi.nlm.nih.gov/?term=Adjunctive+rifampin+for+the+treatment+of+Staphylococcus+aureus+bacteraemia+with+deep+infection+foci)
   Ryder JH, Tong SYC, Gallagher JC, et al. Open Forum Infect Dis 2023.

81. [Aminoglycoside uptake requires the proton-motive force energy-dependent phase II](https://pubmed.ncbi.nlm.nih.gov/?term=Aminoglycoside+uptake+requires+the+proton-motive+force+energy-dependent+phase+II)
   Taber HW, Mueller JP, Miller PF, Arrow AS. Microbiol Rev 1987. The mechanism behind the model's OXY term.

82. [Anaerobiosis abolishes the bactericidal activity of aminoglycosides](https://pubmed.ncbi.nlm.nih.gov/?term=Anaerobiosis+abolishes+the+bactericidal+activity+of+aminoglycosides)
   Bryan LE, Kwan S. Antimicrob Agents Chemother 1983.

83. [Synergism between penicillin and aminoglycosides against enterococci mechanism](https://pubmed.ncbi.nlm.nih.gov/?term=Synergism+between+penicillin+and+aminoglycosides+against+enterococci+mechanism)
   Moellering RC Jr, Weinberg AN. J Clin Invest 1971.

84. [Rifampin activity against intracellular and biofilm-embedded Staphylococcus aureus](https://pubmed.ncbi.nlm.nih.gov/?term=Rifampin+activity+against+intracellular+and+biofilm-embedded+Staphylococcus+aureus)
   Zheng Z, Stewart PS. Antimicrob Agents Chemother 2002.

85. [Role of rifampin against staphylococcal biofilm infections in vitro in animal models and in orthopaedic device infections](https://pubmed.ncbi.nlm.nih.gov/?term=Role+of+rifampin+against+staphylococcal+biofilm+infections+in+vitro+in+animal+models+and+in+orthopaedic+device+infections)
   Zimmerli W, Sendi P. Antimicrob Agents Chemother 2019. The strongest case for rifampicin, and the setting the model agrees with.

86. [Role of rifampin for treatment of orthopedic implant-related staphylococcal infections a randomized controlled trial](https://pubmed.ncbi.nlm.nih.gov/?term=Role+of+rifampin+for+treatment+of+orthopedic+implant-related+staphylococcal+infections+a+randomized+controlled+trial)
   Zimmerli W, Widmer AF, Blatter M, et al. JAMA 1998.

87. [Emergence of rifampin resistance during therapy of Staphylococcus aureus infection](https://pubmed.ncbi.nlm.nih.gov/?term=Emergence+of+rifampin+resistance+during+therapy+of+Staphylococcus+aureus+infection)
   Zinner SH, Lagast H, Klastersky J. J Infect Dis 1981. Monotherapy failure by pre-existing rpoB mutants.

88. [Mutations in rpoB conferring rifampin resistance in Staphylococcus aureus and their fitness cost](https://pubmed.ncbi.nlm.nih.gov/?term=Mutations+in+rpoB+conferring+rifampin+resistance+in+Staphylococcus+aureus+and+their+fitness+cost)
   Wichelhaus TA, Boddinghaus B, Besier S, et al. Antimicrob Agents Chemother 2002. Source of the model's FIT_COST.

89. [Rifampicin drug interactions induction of CYP3A4 and P-glycoprotein](https://pubmed.ncbi.nlm.nih.gov/?term=Rifampicin+drug+interactions+induction+of+CYP3A4+and+P-glycoprotein)
   Niemi M, Backman JT, Fromm MF, et al. Clin Pharmacokinet 2003. The model's RIFIND state and FRIND_* parameters.

---

## 9. 내성·관용의 집단유전학 (Population Genetics of Resistance and Tolerance)

90. [Mutations of bacteria from virus sensitivity to virus resistance](https://pubmed.ncbi.nlm.nih.gov/?term=Mutations+of+bacteria+from+virus+sensitivity+to+virus+resistance)
   Luria SE, Delbruck M. Genetics 1943. Why the mutant count at diagnosis exceeds f_mut x N; the model reproduces the excess without being told to.

91. [The mutant selection window and antimicrobial resistance](https://pubmed.ncbi.nlm.nih.gov/?term=The+mutant+selection+window+and+antimicrobial+resistance)
   Drlica K. J Antimicrob Chemother 2003.

92. [Mutant prevention concentration as a measure of fluoroquinolone potency against mycobacteria and staphylococci](https://pubmed.ncbi.nlm.nih.gov/?term=Mutant+prevention+concentration+as+a+measure+of+fluoroquinolone+potency+against+mycobacteria+and+staphylococci)
   Dong Y, Zhao X, Domagala J, Drlica K. Antimicrob Agents Chemother 1999.

93. [Antibiotic tolerance facilitates the evolution of resistance](https://pubmed.ncbi.nlm.nih.gov/?term=Antibiotic+tolerance+facilitates+the+evolution+of+resistance)
   Levin-Reisman I, Ronin I, Gefen O, et al. Science 2017. Tolerance is the on-ramp to resistance, which is why the model tracks both.

94. [Distinguishing between resistance tolerance and persistence to antibiotic treatment](https://pubmed.ncbi.nlm.nih.gov/?term=Distinguishing+between+resistance+tolerance+and+persistence+to+antibiotic+treatment)
   Brauner A, Fridman O, Gefen O, Balaban NQ. Nat Rev Microbiol 2016.

95. [Heteroresistance a cause of unexplained antibiotic treatment failure](https://pubmed.ncbi.nlm.nih.gov/?term=Heteroresistance+a+cause+of+unexplained+antibiotic+treatment+failure)
   Andersson DI, Nicoloff H, Hjort K. Nat Rev Microbiol 2019.

96. [The biological cost of antibiotic resistance and compensatory evolution](https://pubmed.ncbi.nlm.nih.gov/?term=The+biological+cost+of+antibiotic+resistance+and+compensatory+evolution)
   Andersson DI, Hughes D. Nat Rev Microbiol 2010.

97. [Combination antibiotic therapy and the suppression of resistance emergence](https://pubmed.ncbi.nlm.nih.gov/?term=Combination+antibiotic+therapy+and+the+suppression+of+resistance+emergence)
   Drusano GL, Louie A, MacGowan A, Hope W. Antimicrob Agents Chemother 2016. The condition the model encodes as Comb_suppress.

---

## 10. 색전증과 수술 (Embolic Events and Surgery)

98. [The relationship between the initiation of antimicrobial therapy and the incidence of stroke in infective endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=The+relationship+between+the+initiation+of+antimicrobial+therapy+and+the+incidence+of+stroke+in+infective+endocarditis)
   Dickerman SA, Abrutyn E, Barsic B, et al. Am Heart J 2007. Target T3: embolic hazard is front-loaded and falls steeply once therapy is effective.

99. [Risk of embolism and death in infective endocarditis prognostic value of echocardiography](https://pubmed.ncbi.nlm.nih.gov/?term=Risk+of+embolism+and+death+in+infective+endocarditis+prognostic+value+of+echocardiography)
   Thuny F, Di Salvo G, Belliard O, et al. Circulation 2005. Target T4: the >10 mm and >15 mm thresholds.

100. [Early surgery versus conventional treatment for infective endocarditis EASE trial](https://pubmed.ncbi.nlm.nih.gov/?term=Early+surgery+versus+conventional+treatment+for+infective+endocarditis+EASE+trial)
   Kang DH, Kim YJ, Kim SH, et al. N Engl J Med 2012. Target T9: the composite is driven by embolic events.

101. [Impact of early surgery on embolic events in patients with infective endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Impact+of+early+surgery+on+embolic+events+in+patients+with+infective+endocarditis)
   Kim DH, Kang DH, Lee MZ, et al. Circulation 2010.

102. [Association between surgical indications operative risk and clinical outcome in infective endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Association+between+surgical+indications+operative+risk+and+clinical+outcome+in+infective+endocarditis)
   Chu VH, Park LP, Athan E, et al. Circulation 2015.

103. [Aspirin for the prevention of embolic events in infective endocarditis a randomized controlled trial](https://pubmed.ncbi.nlm.nih.gov/?term=Aspirin+for+the+prevention+of+embolic+events+in+infective+endocarditis+a+randomized+controlled+trial)
   Chan KL, Dumesnil JG, Cujec B, et al. J Am Coll Cardiol 2003. Target T7: no reduction in embolic events, more bleeding.

104. [Antiplatelet and anticoagulant therapy in infective endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Antiplatelet+and+anticoagulant+therapy+in+infective+endocarditis)
   Anavekar NS, Schultz JC, De Sa DDC, et al. Clin Infect Dis 2007.

105. [Neurological complications of infective endocarditis risk factors outcome and impact of cardiac surgery](https://pubmed.ncbi.nlm.nih.gov/?term=Neurological+complications+of+infective+endocarditis+risk+factors+outcome+and+impact+of+cardiac+surgery)
   Garcia-Cabrera E, Fernandez-Hidalgo N, Almirante B, et al. Circulation 2013.

106. [Silent cerebral embolism in infective endocarditis detected by magnetic resonance imaging](https://pubmed.ncbi.nlm.nih.gov/?term=Silent+cerebral+embolism+in+infective+endocarditis+detected+by+magnetic+resonance+imaging)
   Duval X, Iung B, Klein I, et al. Ann Intern Med 2010.

107. [Infectious intracranial aneurysms in infective endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Infectious+intracranial+aneurysms+in+infective+endocarditis)
   Ducruet AF, Hickman ZL, Zacharia BE, et al. Neurosurg Rev 2010.

108. [Timing of surgery in infective endocarditis with neurological complications](https://pubmed.ncbi.nlm.nih.gov/?term=Timing+of+surgery+in+infective+endocarditis+with+neurological+complications)
   Yoshioka D, Toda K, Sakaguchi T, et al. Ann Thorac Surg 2014.

---

## 11. 심부전·판막파괴·판막주위 확장 (Heart Failure, Valve Destruction, Perivalvular Extension)

109. [Heart failure in infective endocarditis prevalence and prognostic implications](https://pubmed.ncbi.nlm.nih.gov/?term=Heart+failure+in+infective+endocarditis+prevalence+and+prognostic+implications)
   Nadji G, Rusinaru D, Remadi JP, et al. Heart 2009.

110. [Perivalvular abscess in infective endocarditis diagnosis and management](https://pubmed.ncbi.nlm.nih.gov/?term=Perivalvular+abscess+in+infective+endocarditis+diagnosis+and+management)
   Graupner C, Vilacosta I, San Roman J, et al. J Am Coll Cardiol 2002.

111. [Prolonged PR interval as a marker of perivalvular extension in aortic valve endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Prolonged+PR+interval+as+a+marker+of+perivalvular+extension+in+aortic+valve+endocarditis)
   DiNubile MJ, Calderwood SB, Steinhaus DM, Karchmer AW. Am J Cardiol 1986. The model's PRMS state.

112. [Natriuretic peptides in infective endocarditis prognostic value](https://pubmed.ncbi.nlm.nih.gov/?term=Natriuretic+peptides+in+infective+endocarditis+prognostic+value)
   Kahveci G, Bayrak F, Mutlu B, et al. Am J Cardiol 2007.

113. [Acute severe mitral regurgitation pathophysiology and haemodynamics](https://pubmed.ncbi.nlm.nih.gov/?term=Acute+severe+mitral+regurgitation+pathophysiology+and+haemodynamics)
   Stout KK, Verrier ED. Circulation 2009. Why an acute regurgitant lesion decompensates at a normal ventricular size.

114. [Staphylococcus aureus alpha-toxin and host tissue destruction](https://pubmed.ncbi.nlm.nih.gov/?term=Staphylococcus+aureus+alpha-toxin+and+host+tissue+destruction)
   Berube BJ, Bubeck Wardenburg J. Toxins 2013.

---

## 12. 지속 균혈증·전이 병소·세포내 저장고 (Persistent Bacteraemia, Metastatic Foci, the Intracellular Reservoir)

115. [Persistence in Staphylococcus aureus bacteremia incidence characteristics of patients and outcome](https://pubmed.ncbi.nlm.nih.gov/?term=Persistence+in+Staphylococcus+aureus+bacteremia+incidence+characteristics+of+patients+and+outcome)
   Khatib R, Johnson LB, Fakih MG, et al. Scand J Infect Dis 2006. Target T1: the time-to-clearance benchmarks.

116. [Slower clearance of methicillin-resistant than methicillin-susceptible Staphylococcus aureus bacteraemia](https://pubmed.ncbi.nlm.nih.gov/?term=Slower+clearance+of+methicillin-resistant+than+methicillin-susceptible+Staphylococcus+aureus+bacteraemia)
   Levine DP, Fromm BS, Reddy BR. Ann Intern Med 1991.

117. [Persistent methicillin-resistant Staphylococcus aureus bacteremia clinical and molecular characteristics](https://pubmed.ncbi.nlm.nih.gov/?term=Persistent+methicillin-resistant+Staphylococcus+aureus+bacteremia+clinical+and+molecular+characteristics)
   Fowler VG Jr, Sakoulas G, McIntyre LM, et al. J Infect Dis 2004.

118. [Clinical management of Staphylococcus aureus bacteraemia](https://pubmed.ncbi.nlm.nih.gov/?term=Clinical+management+of+Staphylococcus+aureus+bacteraemia)
   Holland TL, Arnold C, Fowler VG Jr. JAMA 2014.

119. [Intracellular Staphylococcus aureus and the failure of antibiotic therapy](https://pubmed.ncbi.nlm.nih.gov/?term=Intracellular+Staphylococcus+aureus+and+the+failure+of+antibiotic+therapy)
   Fraunholz M, Sinha B. Front Cell Infect Microbiol 2012. The model's INTRA compartment.

120. [Intracellular activity of antibiotics against Staphylococcus aureus in a model of THP-1 macrophages](https://pubmed.ncbi.nlm.nih.gov/?term=Intracellular+activity+of+antibiotics+against+Staphylococcus+aureus+in+a+model+of+THP-1+macrophages)
   Barcia-Macay M, Seral C, Mingeot-Leclercq MP, et al. Antimicrob Agents Chemother 2006. Source of the model's FCELL_* ordering: rifampicin high, vancomycin and daptomycin very low.

121. [Activity of daptomycin against intracellular Staphylococcus aureus pH dependence](https://pubmed.ncbi.nlm.nih.gov/?term=Activity+of+daptomycin+against+intracellular+Staphylococcus+aureus+pH+dependence)
   Lemaire S, Van Bambeke F, Mingeot-Leclercq MP, Tulkens PM. Antimicrob Agents Chemother 2007.

122. [Blood volume cultured and the yield of blood cultures in adults](https://pubmed.ncbi.nlm.nih.gov/?term=Blood+volume+cultured+and+the+yield+of+blood+cultures+in+adults)
   Lamy B, Dargere S, Arendrup MC, et al. Front Microbiol 2016. Why the model computes culture positivity as a Poisson function of volume.

---

## 13. 정맥에서 경구로 — POET와 치료 기간 (From IV to Oral: POET and Duration of Therapy)

123. [Partial oral versus intravenous antibiotic treatment of endocarditis POET](https://pubmed.ncbi.nlm.nih.gov/?term=Partial+oral+versus+intravenous+antibiotic+treatment+of+endocarditis+POET)
   Iversen K, Ihlemann N, Gill SU, et al. N Engl J Med 2019. Target T8.

124. [Partial oral antibiotic treatment for endocarditis five-year outcomes of the POET trial](https://pubmed.ncbi.nlm.nih.gov/?term=Partial+oral+antibiotic+treatment+for+endocarditis+five-year+outcomes+of+the+POET+trial)
   Pries-Heje MM, Wiingaard C, Ihlemann N, et al. N Engl J Med 2022.

125. [Oral versus intravenous antibiotics for bone and joint infection OVIVA](https://pubmed.ncbi.nlm.nih.gov/?term=Oral+versus+intravenous+antibiotics+for+bone+and+joint+infection+OVIVA)
   Li HK, Rombach I, Zambellas R, et al. N Engl J Med 2019. The parallel result in another deep-seated infection.

126. [Outpatient parenteral antimicrobial therapy for infective endocarditis feasibility and outcome](https://pubmed.ncbi.nlm.nih.gov/?term=Outpatient+parenteral+antimicrobial+therapy+for+infective+endocarditis+feasibility+and+outcome)
   Duncan CJA, Barr DA, Seaton RA. J Antimicrob Chemother 2013.

127. [Two-week combination therapy for right-sided Staphylococcus aureus endocarditis in injection drug users](https://pubmed.ncbi.nlm.nih.gov/?term=Two-week+combination+therapy+for+right-sided+Staphylococcus+aureus+endocarditis+in+injection+drug+users)
   Chambers HF, Miller RT, Newman MD. Ann Intern Med 1988. Why right-sided disease is physically a different problem.

128. [Duration of antibiotic therapy for bacteraemia and endocarditis current evidence](https://pubmed.ncbi.nlm.nih.gov/?term=Duration+of+antibiotic+therapy+for+bacteraemia+and+endocarditis+current+evidence)
   Hanretty AM, Gallagher JC. Pharmacotherapy 2018.

---

## 14. 영상·감시 (Imaging and Monitoring)

129. [18F-FDG PET/CT in the diagnosis of prosthetic valve endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=18F-FDG+PET%2FCT+in+the+diagnosis+of+prosthetic+valve+endocarditis)
   Saby L, Laas O, Habib G, et al. J Am Coll Cardiol 2013.

130. [Role of multi-modality imaging in infective endocarditis EACVI recommendations](https://pubmed.ncbi.nlm.nih.gov/?term=Role+of+multi-modality+imaging+in+infective+endocarditis+EACVI+recommendations)
   Habib G, Badano L, Tribouilloy C, et al. Eur J Echocardiogr 2010.

131. [Value of C-reactive protein kinetics in infective endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Value+of+C-reactive+protein+kinetics+in+infective+endocarditis)
   Verhagen DWM, Hermanides J, Korevaar JC, et al. Clin Infect Dis 2008. The model's CRP trajectory readout.

132. [Procalcitonin in the diagnosis of infective endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Procalcitonin+in+the+diagnosis+of+infective+endocarditis)
   Yu CW, Juan LI, Hsu SC, et al. Am J Emerg Med 2013.

133. [Cardiac computed tomography angiography for perivalvular complications of infective endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Cardiac+computed+tomography+angiography+for+perivalvular+complications+of+infective+endocarditis)
   Feuchtner GM, Stolzmann P, Dichtl W, et al. J Am Coll Cardiol 2009.

---

## 15. QSP·PK/PD 방법론 (QSP and PK/PD Methodology)

134. [mrgsolve simulation from ODE-based population PK/PD and systems pharmacology models](https://pubmed.ncbi.nlm.nih.gov/?term=mrgsolve+simulation+from+ODE-based+population+PK%2FPD+and+systems+pharmacology+models)
   Baron KT. R package documentation and vignettes; the simulation engine used here.

135. [Pharmacokinetic-pharmacodynamic modelling of antibacterial effects and the emergence of resistance](https://pubmed.ncbi.nlm.nih.gov/?term=Pharmacokinetic-pharmacodynamic+modelling+of+antibacterial+effects+and+the+emergence+of+resistance)
   Nielsen EI, Friberg LE. Pharmacol Rev 2013. The canonical review of the model structures used in this file.

136. [Semimechanistic pharmacokinetic-pharmacodynamic model for assessment of activity of antibacterial agents from time-kill curve experiments](https://pubmed.ncbi.nlm.nih.gov/?term=Semimechanistic+pharmacokinetic-pharmacodynamic+model+for+assessment+of+activity+of+antibacterial+agents+from+time-kill+curve+experiments)
   Nielsen EI, Viberg A, Lowdin E, et al. Antimicrob Agents Chemother 2007. The growing / resting compartment structure this model generalises to three shells.

137. [Mechanism-based pharmacokinetic-pharmacodynamic modeling of antimicrobial drug effects](https://pubmed.ncbi.nlm.nih.gov/?term=Mechanism-based+pharmacokinetic-pharmacodynamic+modeling+of+antimicrobial+drug+effects)
   Czock D, Keller F. J Pharmacokinet Pharmacodyn 2007.

138. [Pharmacokinetic-pharmacodynamic modelling of bacterial kill and regrowth including adaptive resistance](https://pubmed.ncbi.nlm.nih.gov/?term=Pharmacokinetic-pharmacodynamic+modelling+of+bacterial+kill+and+regrowth+including+adaptive+resistance)
   Mohamed AF, Nielsen EI, Cars O, Friberg LE. Antimicrob Agents Chemother 2012.

139. [Suppression of emergence of resistance in pathogenic bacteria why do we know so little](https://pubmed.ncbi.nlm.nih.gov/?term=Suppression+of+emergence+of+resistance+in+pathogenic+bacteria+why+do+we+know+so+little)
   Drusano GL, Louie A, Deziel M, Gumbo T. Curr Opin Infect Dis 2006.

140. [Model-informed drug development for anti-infectives state of the art and future](https://pubmed.ncbi.nlm.nih.gov/?term=Model-informed+drug+development+for+anti-infectives+state+of+the+art+and+future)
   Rayner CR, Smith PF, Andes D, et al. Clin Pharmacol Ther 2021.

141. [Quantitative systems pharmacology in the age of model-informed drug development](https://pubmed.ncbi.nlm.nih.gov/?term=Quantitative+systems+pharmacology+in+the+age+of+model-informed+drug+development)
   Bradshaw EL, Spilker ME, Zang R, et al. CPT Pharmacometrics Syst Pharmacol 2019.

142. [Applications of quantitative systems pharmacology in anti-infective drug development](https://pubmed.ncbi.nlm.nih.gov/?term=Applications+of+quantitative+systems+pharmacology+in+anti-infective+drug+development)
   Nikolaou M, Tam VH. Pharm Res 2006.

143. [A new modeling approach to the effect of antimicrobial agents on heterogeneous microbial populations](https://pubmed.ncbi.nlm.nih.gov/?term=A+new+modeling+approach+to+the+effect+of+antimicrobial+agents+on+heterogeneous+microbial+populations)
   Nikolaou M, Tam VH. J Math Biol 2006.

---

## 16. 실험적 심내막염 모델 (Experimental Endocarditis Models — where most of the mechanism is measured)

144. [Experimental bacterial endocarditis in rabbits a model for the study of antibiotic therapy](https://pubmed.ncbi.nlm.nih.gov/?term=Experimental+bacterial+endocarditis+in+rabbits+a+model+for+the+study+of+antibiotic+therapy)
   Perlman BB, Freedman LR. Yale J Biol Med 1971.

145. [Experimental endocarditis models and their relevance to human disease](https://pubmed.ncbi.nlm.nih.gov/?term=Experimental+endocarditis+models+and+their+relevance+to+human+disease)
   Entenza JM, Moreillon P. Curr Opin Infect Dis 2008.

146. [Correlation of antibiotic bactericidal activity in vitro and in experimental endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Correlation+of+antibiotic+bactericidal+activity+in+vitro+and+in+experimental+endocarditis)
   Bayer AS, Lam K. Antimicrob Agents Chemother 1985.

147. [Efficacy of daptomycin in experimental endocarditis due to methicillin-resistant Staphylococcus aureus](https://pubmed.ncbi.nlm.nih.gov/?term=Efficacy+of+daptomycin+in+experimental+endocarditis+due+to+methicillin-resistant+Staphylococcus+aureus)
   Sakoulas G, Eliopoulos GM, Alder J, Eliopoulos CT. Antimicrob Agents Chemother 2003.

148. [Comparative efficacy of vancomycin and nafcillin in experimental MSSA endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=Comparative+efficacy+of+vancomycin+and+nafcillin+in+experimental+MSSA+endocarditis)
   Chambers HF, Miller RT. J Infect Dis 1987.

149. [In vivo emergence of daptomycin resistance in experimental endocarditis](https://pubmed.ncbi.nlm.nih.gov/?term=In+vivo+emergence+of+daptomycin+resistance+in+experimental+endocarditis)
   Silverman JA, Oliver N, Andrew T, Li T. Antimicrob Agents Chemother 2001.

---

## 이 모델의 핵심 주장과 근거 문헌의 대응 (Claim-to-Evidence Map)

| 모델의 주장 | 코드상의 위치 | 주된 근거 |
|---|---|---|
| 증식물 내 균 밀도는 10^9-10^11 CFU/g이고 정지기(stationary phase)이다 | `DENSCAP`, `MUMAX`, `AVL_*` | Durack & Beeson 1972 (3편) |
| 약물은 확산으로만 들어가며 깊이에 따라 지수적으로 감쇠한다 | `LAM_*`, `exp(-d/lambda)` | Cremieux 1989 autoradiography; Stewart 1996 |
| β-락탐의 살균은 성장속도에 비례하므로 코어에서 소실된다 | `FSTAT_BL = 0.02` | Craig 1998; Lewis 2007; Nielsen 2007 |
| 답토마이신은 비분열 세포에도 작용한다 | `FSTAT_DAP = 0.55` | Mascio 2007 |
| 아미노글리코사이드 흡수는 PMF 의존적이므로 무산소 코어에서 실패한다 | `OXY = AVL^1.5` | Taber 1987; Bryan & Kwan 1983 → **T5** |
| rpoB 변이체는 진단 시점에 이미 존재한다 | `MUTRATE * mu * B` | Luria & Delbruck 1943; Wichelhaus 2002 |
| 반코마이신 AKI는 AUC에 급격히 의존하고 양성 되먹임을 만든다 | `CAKI50_V`, `E_AKI` | Rybak 2020; Aljefri 2019 → **T10** |
| 색전 위험은 표면 *세균* 집단을 따르지 기질을 따르지 않는다 | `FRAG = f(NSG)` | Dickerman 2007; Thuny 2005 → **T3, T4, T7** |
| 세포내 저장고는 리팜피신만 따라 들어간다 | `FCELL_*` | Barcia-Macay 2006; Lemaire 2007 |
| 경로(정맥/경구)가 아니라 자유 약물 농도와 lambda가 결정한다 | `IE_poet()` | POET 2019, 2022; OVIVA 2019 → **T8** |

---

## 면책 (Disclaimer)

본 모델과 참고문헌 목록은 **교육 및 연구 목적**으로 작성되었습니다. 실제 임상
의사결정, 처방, 항생제 선택, 수술 시기 결정, 또는 규제 제출에 사용해서는 안
됩니다. 감염성 심내막염의 치료는 반드시 최신 지침(2023 ESC, AHA)과 감염내과·
심장내과·심장외과·임상약학이 참여하는 **엔도카르디티스 팀**의 판단을 따라야
합니다.
