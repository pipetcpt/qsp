# 인공관절 감염 (Prosthetic Joint Infection, PJI) — 참고문헌

> Implant-associated *Staphylococcus aureus* osteomyelitis — QSP 모델 근거 문헌
>
> **PMID 확인 방법**: 아래 68편의 PMID는 모두 NCBI E-utilities(`esearch` + `esummary`)로
> 제목·저널·연도·제1저자를 조회해 대조한 것입니다. 기억에 의존해 적은 번호가 아니며,
> 각 항목의 저널·연도는 PubMed가 반환한 값 그대로입니다.

---

## 0. 이 모델의 네 기둥과 그 근거 (The four pillars and where they come from)

| 기둥 | 정량적 주장 | 1차 근거 |
|------|-------------|----------|
| ① 이물이 ID50를 옮긴다 | 이물 존재 시 감염 성립에 필요한 균량이 10⁴–10⁵배 감소 | Elek & Conen 1957 [13499821] · Zimmerli 1982 [7119479] |
| ② 살균은 C_bone,free / MBEC | 전신 투여로 달성 가능한 유리 골농도는 MBEC의 0.004–0.28 | Ceri 1999 [10325322] · Landersdorfer 2009 [19271782] · Thabit 2019 [30772469] |
| ③ 수술은 변이 공급 조작 | P(사전 rpoB 변이체) = 1 − e^(−μN), μ ≈ 10⁻⁸ | O'Neill 2001 [11328777] · Achermann 2013 [22987291] · Drlica 2007 [17278059] |
| ④ 바이오필름 성숙은 시계 | DAIR 성공률은 증상 지속기간의 함수 | Byren 2009 [19336454] · Lora-Tamayo 2013 [22942204] · Nishitani 2015 [25820925] |

---

## 1. 총설 및 질환 개관 (Reviews & disease overview)

1. Zimmerli W, Trampuz A, Ochsner PE. **Prosthetic-joint infections.** *N Engl J Med.* 2004. — 이 분야의 기준 총설. 병인·진단·치료 알고리즘의 원형.
   <https://pubmed.ncbi.nlm.nih.gov/15483283/>
2. Tande AJ, Patel R. **Prosthetic joint infection.** *Clin Microbiol Rev.* 2014.
   <https://pubmed.ncbi.nlm.nih.gov/24696437/>
3. Kavanagh N, et al. **Staphylococcal Osteomyelitis: Disease Progression, Treatment Challenges, and Future Directions.** *Clin Microbiol Rev.* 2018.
   <https://pubmed.ncbi.nlm.nih.gov/29444953/>
4. Masters EA, et al. **Skeletal infections: microbial pathogenesis, immunity and clinical management.** *Nat Rev Microbiol.* 2022.
   <https://pubmed.ncbi.nlm.nih.gov/35169289/>
5. Masters EA, et al. **Evolving concepts in bone infection: redefining "biofilm", "acute vs. chronic osteomyelitis", "the immune proteome" and "local antibiotic therapy".** *Bone Res.* 2019.
   <https://pubmed.ncbi.nlm.nih.gov/31646012/>
6. Lew DP, Waldvogel FA. **Osteomyelitis.** *Lancet.* 2004.
   <https://pubmed.ncbi.nlm.nih.gov/15276398/>
7. Izakovicova P, Borens O, Trampuz A. **Periprosthetic joint infection: current concepts and outlook.** *EFORT Open Rev.* 2019.
   <https://pubmed.ncbi.nlm.nih.gov/31423332/>
8. Zimmerli W, Sendi P. **Orthopaedic biofilm infections.** *APMIS.* 2017.
   <https://pubmed.ncbi.nlm.nih.gov/28407423/>

---

## 2. 기둥 ① — 이물과 감염 성립 균량 (Foreign body and the infective dose)

9. Elek SD, Conen PE. **The virulence of *Staphylococcus pyogenes* for man; a study of the problems of wound infection.** *Br J Exp Pathol.* 1957. — 사람 피부에 실크 봉합사가 있을 때와 없을 때 병변을 만드는 균량이 10⁴–10⁵배 차이 난다는 고전 실험. 모델의 `IMPL`/`FBFACT`/`KATT` 조합이 재현하려는 대상.
   <https://pubmed.ncbi.nlm.nih.gov/13499821/>
10. Zimmerli W, Waldvogel FA, Vaudaux P, Nydegger UE. **Pathogenesis of foreign body infection: description and characteristics of an animal model.** *J Infect Dis.* 1982. — 조직 케이지 모델. 국소 호중구 기능 결손(모델의 `FBFACT`)의 실험적 근거.
    <https://pubmed.ncbi.nlm.nih.gov/7119479/>
11. Thurlow LR, et al. **Staphylococcus aureus biofilms prevent macrophage phagocytosis and attenuate inflammation in vivo.** *J Immunol.* 2011. — 좌절 식균작용(`FRUST`).
    <https://pubmed.ncbi.nlm.nih.gov/21525381/>
12. Heim CE, et al. **Myeloid-derived suppressor cells contribute to Staphylococcus aureus orthopedic biofilm infection.** *J Immunol.* 2014. — 모델의 `MDSC`·`EMDSC` 항.
    <https://pubmed.ncbi.nlm.nih.gov/24646737/>
13. Heim CE, et al. **Interleukin-10 production by myeloid-derived suppressor cells contributes to bacterial persistence during Staphylococcus aureus orthopedic biofilm infection.** *J Leukoc Biol.* 2015.
    <https://pubmed.ncbi.nlm.nih.gov/26232453/>

---

## 3. 바이오필름 생물학 (Biofilm biology)

14. Costerton JW, Stewart PS, Greenberg EP. **Bacterial biofilms: a common cause of persistent infections.** *Science.* 1999.
    <https://pubmed.ncbi.nlm.nih.gov/10334980/>
15. Stewart PS, Costerton JW. **Antibiotic resistance of bacteria in biofilms.** *Lancet.* 2001. — 확산 지연·대사 정지·지속균의 3중 기전.
    <https://pubmed.ncbi.nlm.nih.gov/11463434/>
16. Hall CW, Mah TF. **Molecular mechanisms of biofilm-based antibiotic resistance and tolerance in pathogenic bacteria.** *FEMS Microbiol Rev.* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/28369412/>
17. Moormeier DE, Bayles KW. **Staphylococcus aureus biofilm: a complex developmental organism.** *Mol Microbiol.* 2017. — 부착 → 증식 → 탈리의 발달 단계(모델의 `KATT`/`EPS`/`KDISP`).
    <https://pubmed.ncbi.nlm.nih.gov/28142193/>
18. McConoughey SJ, et al. **Biofilms in periprosthetic orthopedic infections.** *Future Microbiol.* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/25302955/>
19. Nishitani K, et al. **Quantifying the natural history of biofilm formation in vivo during the establishment of chronic implant-associated Staphylococcus aureus osteomyelitis in mice.** *J Orthop Res.* 2015. — 모델의 EPS 성숙 시상수(~10일)의 실험적 대응.
    <https://pubmed.ncbi.nlm.nih.gov/25820925/>
20. Foster TJ, Geoghegan JA, Ganesh VK, Höök M. **Adhesion, invasion and evasion: the many functions of the surface proteins of Staphylococcus aureus.** *Nat Rev Microbiol.* 2014. — MSCRAMM(FnBPA/ClfA/SdrC-E).
    <https://pubmed.ncbi.nlm.nih.gov/24336184/>
21. Geoghegan JA, Foster TJ. **Cell Wall-Anchored Surface Proteins of Staphylococcus aureus: Many Proteins, Multiple Functions.** *Curr Top Microbiol Immunol.* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/26667044/>
22. Otto M. **Phenol-soluble modulins.** *Int J Med Microbiol.* 2014. — PSM: 기질 아밀로이드·호중구 용해·분산.
    <https://pubmed.ncbi.nlm.nih.gov/24447915/>
23. Cheung GY, Joo HS, Chatterjee SS, Otto M. **Phenol-soluble modulins — critical determinants of staphylococcal virulence.** *FEMS Microbiol Rev.* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/24372362/>

---

## 4. 표현형 아집단 — 지속균·SCV·세포내 저장고 (Phenotypic subpopulations)

24. Lewis K. **Persister cells.** *Annu Rev Microbiol.* 2010. — 모델의 `NPER`·`KPER`·`KWAKE`.
    <https://pubmed.ncbi.nlm.nih.gov/20528688/>
25. Proctor RA, et al. **Small colony variants: a pathogenic form of bacteria that facilitates persistent and recurrent infections.** *Nat Rev Microbiol.* 2006. — `NSCV`.
    <https://pubmed.ncbi.nlm.nih.gov/16541137/>
26. de Mesy Bentley KL, et al. **Evidence of Staphylococcus aureus Deformation, Proliferation, and Migration in Canaliculi of Live Cortical Bone in Murine Models of Osteomyelitis.** *J Bone Miner Res.* 2017. — 골소관 저장고: 모델 `NIC`가 왜 항생제로 잘 죽지 않는지의 근거.
    <https://pubmed.ncbi.nlm.nih.gov/27933662/>
27. Ellington JK, et al. **Intracellular Staphylococcus aureus and antibiotic resistance: implications for treatment of staphylococcal osteomyelitis.** *J Orthop Res.* 2006. — 세포내 균에 대한 항생제별 활성 차이(`CAR`·`PHI*`).
    <https://pubmed.ncbi.nlm.nih.gov/16419973/>
28. Josse J, Velard F, Gangloff SC. **Staphylococcus aureus vs. Osteoblast: Relationship and Consequences in Osteomyelitis.** *Front Cell Infect Microbiol.* 2015.
    <https://pubmed.ncbi.nlm.nih.gov/26636047/>

---

## 5. 기둥 ② — MBEC, 골 침투, PK/PD (The C_bone,free / MBEC ratio)

29. Ceri H, et al. **The Calgary Biofilm Device: new technology for rapid determination of antibiotic susceptibilities of bacterial biofilms.** *J Clin Microbiol.* 1999. — MBEC 개념·측정법의 출처.
    <https://pubmed.ncbi.nlm.nih.gov/10325322/>
30. Landersdorfer CB, Bulitta JB, Kinzig M, Holzgrabe U, Sörgel F. **Penetration of antibacterials into bone: pharmacokinetic, pharmacodynamic and bioanalytical considerations.** *Clin Pharmacokinet.* 2009. — 모델의 `PEN*` 값의 1차 출처.
    <https://pubmed.ncbi.nlm.nih.gov/19271782/>
31. Thabit AK, et al. **Antibiotic penetration into bone and joints: An updated review.** *Int J Infect Dis.* 2019.
    <https://pubmed.ncbi.nlm.nih.gov/30772469/>
32. Drusano GL. **Antimicrobial pharmacodynamics: critical interactions of 'bug and drug'.** *Nat Rev Microbiol.* 2004.
    <https://pubmed.ncbi.nlm.nih.gov/15031728/>
33. Post V, et al. **Vancomycin displays time-dependent eradication of mature Staphylococcus aureus biofilms.** *J Orthop Res.* 2017. — 반코마이신의 바이오필름 제거가 농도가 아니라 노출시간에 의존한다는 관찰(모델의 `EMXVAN` 저값과 `MBCVAN` 고값의 조합).
    <https://pubmed.ncbi.nlm.nih.gov/27175462/>
34. Widmer AF, Frei R, Rajacic Z, Zimmerli W. **Correlation between in vivo and in vitro efficacy of antimicrobial agents against foreign body infections.** *J Infect Dis.* 1990. — 이물 감염에서 in vitro MIC가 in vivo 효능을 예측하지 못한다는 직접 증거.
    <https://pubmed.ncbi.nlm.nih.gov/2355207/>
35. Rybak MJ, et al. **Therapeutic monitoring of vancomycin for serious methicillin-resistant Staphylococcus aureus infections: A revised consensus guideline and review.** *Am J Health Syst Pharm.* 2020. — AUC24 400–600 목표.
    <https://pubmed.ncbi.nlm.nih.gov/32191793/>

---

## 6. 기둥 ③ — 리팜피신, 변이 공급, 내성 (Rifampicin and mutant supply)

36. Zimmerli W, Widmer AF, Blatter M, Frei R, Ochsner PE. **Role of rifampin for treatment of orthopedic implant-related staphylococcal infections: a randomized controlled trial.** *JAMA.* 1998. — 시프로플록사신+리팜피신 12/12(100%) vs 시프로플록사신 단독 7/12(58%). 모델 시나리오 05/06의 보정 기준점.
    <https://pubmed.ncbi.nlm.nih.gov/9605897/>
37. Widmer AF, Gaechter A, Ochsner PE, Zimmerli W. **Antimicrobial treatment of orthopedic implant-related infections with rifampin combinations.** *Clin Infect Dis.* 1992.
    <https://pubmed.ncbi.nlm.nih.gov/1623081/>
38. Zimmerli W, Sendi P. **Role of Rifampin against Staphylococcal Biofilm Infections In Vitro, in Animal Models, and in Orthopedic-Device-Related Infections.** *Antimicrob Agents Chemother.* 2019.
    <https://pubmed.ncbi.nlm.nih.gov/30455229/>
39. Sendi P, Zimmerli W. **The use of rifampin in staphylococcal orthopaedic-device-related infections.** *Clin Microbiol Infect.* 2017.
    <https://pubmed.ncbi.nlm.nih.gov/27746393/>
40. Sendi P, Zimmerli W. **Antimicrobial treatment concepts for orthopaedic device-related infection.** *Clin Microbiol Infect.* 2012.
    <https://pubmed.ncbi.nlm.nih.gov/23046277/>
41. O'Neill AJ, Cove JH, Chopra I. **Mutation frequencies for resistance to fusidic acid and rifampicin in Staphylococcus aureus.** *J Antimicrob Chemother.* 2001. — μ_rpoB ≈ 10⁻⁷–10⁻⁸의 출처. 모델 `MURPOB`.
    <https://pubmed.ncbi.nlm.nih.gov/11328777/>
42. Aubry-Damon H, Soussy CJ, Courvalin P. **Characterization of mutations in the rpoB gene that confer rifampin resistance in Staphylococcus aureus.** *Antimicrob Agents Chemother.* 1998.
    <https://pubmed.ncbi.nlm.nih.gov/9756760/>
43. Achermann Y, et al. **Factors associated with rifampin resistance in staphylococcal periprosthetic joint infections (PJI): a matched case-control study.** *Infection.* 2013. — 임상에서 실제로 리팜피신 내성이 나오는 조건(불충분한 변연절제·단독 노출).
    <https://pubmed.ncbi.nlm.nih.gov/22987291/>
44. Drlica K, Zhao X. **Mutant selection window hypothesis updated.** *Clin Infect Dis.* 2007. — MSW 개념(모델의 `RES_MSW` 노드).
    <https://pubmed.ncbi.nlm.nih.gov/17278059/>

---

## 7. 기둥 ④ — 수술 전략과 시간 (Surgical strategy and the maturity clock)

45. Byren I, et al. **One hundred and twelve infected arthroplasties treated with 'DAIR' (debridement, antibiotics and implant retention): antibiotic duration and outcome.** *J Antimicrob Chemother.* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/19336454/>
46. Lora-Tamayo J, et al. **A large multicenter study of methicillin-susceptible and methicillin-resistant Staphylococcus aureus prosthetic joint infections managed with implant retention.** *Clin Infect Dis.* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/22942204/>
47. Trebse R, Pisot V, Trampuz A. **Treatment of infected retained implants.** *J Bone Joint Surg Br.* 2005.
    <https://pubmed.ncbi.nlm.nih.gov/15736752/>
48. Bejon P, et al. **Two-stage revision for prosthetic joint infection: predictors of outcome and the role of reimplantation microbiology.** *J Antimicrob Chemother.* 2010.
    <https://pubmed.ncbi.nlm.nih.gov/20053693/>
49. Kunutsor SK, Whitehouse MR, Blom AW, Beswick AD. **Re-Infection Outcomes following One- and Two-Stage Surgical Revision of Infected Hip Prosthesis: A Systematic Review and Meta-Analysis.** *PLoS One.* 2015.
    <https://pubmed.ncbi.nlm.nih.gov/26407003/>
50. Bertazzoni Minelli E, Benini A, Magnan B, Bartolozzi P. **Release of gentamicin and vancomycin from temporary human hip spacers in two-stage revision of infected arthroplasty.** *J Antimicrob Chemother.* 2004. — 모델의 스페이서 2-풀 용출 파라미터(`KFAST`/`KSLOW`/`SPCF0`/`SPCS0`)의 근거.
    <https://pubmed.ncbi.nlm.nih.gov/14688051/>

---

## 8. 임상시험 — 항생제 요법 (Clinical trials of antimicrobial regimens)

51. Li HK, et al. (OVIVA). **Oral versus Intravenous Antibiotics for Bone and Joint Infection.** *N Engl J Med.* 2019. — 경구가 정주에 비열등(실패 13.2% vs 14.6%). 모델 시나리오 11/12: 경로가 아니라 골 노출이 결정 인자.
    <https://pubmed.ncbi.nlm.nih.gov/30699315/>
52. Bernard L, et al. (DATIPO). **Antibiotic Therapy for 6 or 12 Weeks for Prosthetic Joint Infection.** *N Engl J Med.* 2021. — 6주는 12주에 비열등하지 **않았다**. 모델 시나리오 13/14가 재현하는 결과.
    <https://pubmed.ncbi.nlm.nih.gov/34042388/>
53. Byren I, et al. **Randomized controlled trial of the safety and efficacy of Daptomycin versus standard-of-care therapy for management of patients with osteomyelitis associated with prosthetic devices undergoing two-stage revision arthroplasty.** *Antimicrob Agents Chemother.* 2012.
    <https://pubmed.ncbi.nlm.nih.gov/22908174/>
54. John AK, et al. **Efficacy of daptomycin in implant-associated infection due to methicillin-resistant Staphylococcus aureus: importance of combination with rifampin.** *Antimicrob Agents Chemother.* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/19364845/>
55. Saleh-Mghir A, et al. **Adjunctive rifampin is crucial to optimizing daptomycin efficacy against rabbit prosthetic joint infection due to methicillin-resistant Staphylococcus aureus.** *Antimicrob Agents Chemother.* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/21825285/>
56. Baldoni D, Haschke M, Rajacic Z, Zimmerli W, Trampuz A. **Linezolid alone or combined with rifampin against methicillin-resistant Staphylococcus aureus in experimental foreign-body infection.** *Antimicrob Agents Chemother.* 2009.
    <https://pubmed.ncbi.nlm.nih.gov/19075065/>

---

## 9. 진단 기준과 바이오마커 (Diagnostic criteria & biomarkers)

57. Parvizi J, Tan TL, Goswami K, et al. **The 2018 Definition of Periprosthetic Hip and Knee Infection: An Evidence-Based and Validated Criteria.** *J Arthroplasty.* 2018. — 모델의 `ICMSC` 산출식.
    <https://pubmed.ncbi.nlm.nih.gov/29551303/>
58. Osmon DR, et al. **Diagnosis and management of prosthetic joint infection: clinical practice guidelines by the Infectious Diseases Society of America.** *Clin Infect Dis.* 2013.
    <https://pubmed.ncbi.nlm.nih.gov/23223583/>
59. Trampuz A, et al. **Sonication of removed hip and knee prostheses for diagnosis of infection.** *N Engl J Med.* 2007. — 바이오필름을 탈리시켜야 배양이 된다는 직접 증거(모델의 `NB` → `DX_SONIC` 경로).
    <https://pubmed.ncbi.nlm.nih.gov/17699815/>
60. Deirmengian C, et al. **Diagnosing periprosthetic joint infection: has the era of the biomarker arrived?** *Clin Orthop Relat Res.* 2014. — α-defensin.
    <https://pubmed.ncbi.nlm.nih.gov/24590839/>
61. Renz N, et al. **Alpha Defensin Lateral Flow Test for Diagnosis of Periprosthetic Joint Infection: Not a Screening but a Confirmatory Test.** *J Bone Joint Surg Am.* 2018.
    <https://pubmed.ncbi.nlm.nih.gov/29715222/>

---

## 10. 골 파괴 기전 (Bone destruction)

62. Nair SP, Meghji S, Wilson M, et al. **Bacterially induced bone destruction: mechanisms and misconceptions.** *Infect Immun.* 1996.
    <https://pubmed.ncbi.nlm.nih.gov/8698454/>
63. Boyce BF, Xing L. **Functions of RANKL/RANK/OPG in bone modeling and remodeling.** *Arch Biochem Biophys.* 2008. — 모델의 `RANKL`/`OPG`/`OCL` 블록.
    <https://pubmed.ncbi.nlm.nih.gov/18395508/>

---

## 11. 약물 독성 및 상호작용 (Toxicity & drug interactions)

64. Lodise TP, et al. **Relationship between initial vancomycin concentration-time profile and nephrotoxicity among hospitalized patients.** *Clin Infect Dis.* 2009. — 모델의 `AUCTHR`/`ENEPH`.
    <https://pubmed.ncbi.nlm.nih.gov/19586413/>
65. Boak LM, et al. **Clinical population pharmacokinetics and toxicodynamics of linezolid.** *Antimicrob Agents Chemother.* 2014. — 리네졸리드 혈소판 감소의 노출-반응(`EPLT`/`KPLT`).
    <https://pubmed.ncbi.nlm.nih.gov/24514086/>
66. Gandelman K, et al. **Unexpected effect of rifampin on the pharmacokinetics of linezolid: in silico and in vitro approaches to explain its mechanism.** *J Clin Pharmacol.* 2011. — 리팜피신에 의한 리네졸리드 AUC ~32% 감소. 모델의 `RIFLZD` 항.
    <https://pubmed.ncbi.nlm.nih.gov/20371736/>

---

## 12. 질병 부담과 비용 (Burden & cost)

67. Kurtz SM, Lau E, Watson H, Schmier JK, Parvizi J. **Economic burden of periprosthetic joint infection in the United States.** *J Arthroplasty.* 2012.
    <https://pubmed.ncbi.nlm.nih.gov/22554729/>
68. Kurtz SM, et al. **Hospital Costs for Unsuccessful Two-Stage Revisions for Periprosthetic Joint Infection.** *J Arthroplasty.* 2022.
    <https://pubmed.ncbi.nlm.nih.gov/34763048/>

---

## 부록: 모델 파라미터가 문헌에서 온 곳 (Parameter provenance)

| 모델 파라미터 | 값 | 근거 |
|---------------|-----|------|
| `PENVAN` / `PENRIF` / `PENLVX` / `PENDAP` / `PENLZD` / `PENCFZ` | 0.20 / 0.50 / 0.50 / 0.15 / 0.45 / 0.20 | [19271782] · [30772469] |
| `MBCVAN` … `MBCRIF` | 512 … 1.0 mg/L | [10325322] · [27175462] · [30455229] |
| `MURPOB` | 1×10⁻⁸ /분열 | [11328777] |
| `FBFACT` (이물 면역결손) | 0.12 | [7119479] · [13499821] |
| `FRUST` (좌절 식균) | 0.03 | [21525381] |
| `KSMDSC` / `EMDSC` | 0.03 /h, 0.70 | [24646737] · [26232453] |
| `KFAST` / `KSLOW` / `SPCF0` / `SPCS0` | 0.05 /h, 0.0015 /h, 1200 mg, 700 mg | [14688051] |
| `AUCTHR` / `ENEPH` | 600 mg·h/L, 0.55 | [19586413] · [32191793] |
| `EPLT` / `KPLT` | 0.60, 9 mg/L | [24514086] |
| `RIFLZD` | 0.50 (LZD AUC −32%) | [20371736] |
| 검증 기준점 (시나리오 05/06) | RIF 병용 100% vs 단독 58% | [9605897] |
| 검증 기준점 (시나리오 11/12) | 경구 비열등 | [30699315] |
| 검증 기준점 (시나리오 13/14) | 6주 ≠ 12주 | [34042388] |

---

*본 문헌 목록은 교육·연구 목적의 QSP 모델 문서의 일부입니다. 임상 의사결정에 직접 사용하지 마십시오.*
