# 인플루엔자 A — 참고문헌 (References)

이 목록의 모든 항목은 PubMed E-utilities로 **조회하여 확인한 서지정보**이며,
저자·연도·저널·권·페이지는 PubMed 레코드에서 그대로 가져왔다. 기억에 의존해
적은 PMID는 포함하지 않았다.

각 절의 머리글은 그 문헌들이 모델의 **어느 부분을 지지하는지**를 밝힌다.
파라미터의 출처 표기는 `flu_mrgsolve_model.R`의 `$PARAM` 블록에 `[LIT]`
(문헌값) / `[CAL]` (본 모델에서 보정) / `[ASM]` (구조적 가정)으로 달려 있고,
`[CAL]`과 `[ASM]`은 문헌값이 아니라는 뜻이다.

> **보정에 실제로 쓰인 자료는 두 건뿐이다.** CAPSTONE-1(성인 3상)의 세 평가
> 변수와 인체 감염 실험의 자연 경과. 나머지 문헌은 구조를 정당화하거나
> 결과를 대조하는 데 쓰였고, 파라미터를 맞추는 데는 쓰이지 않았다.


## 1. 숙주 내 바이러스 동역학과 모델링 · Within-host viral dynamics and modelling

이 모델의 골격 — 표적세포 제한 복제 루프, eclipse 단계, 인터페론에 의한 refractory 전환 — 은 아래 문헌의 구조를 그대로 따른다. β와 p는 역가 자료만으로는 개별 식별이 불가능하고 곱 β·T₀·p 만이 성장속도로 고정된다는 점도 이들 논문의 공통된 지적이다.

- Baccam P et al. Kinetics of influenza A virus infection in humans. *J Virol* 2006;80:7590-9. [PMID 16840338](https://pubmed.ncbi.nlm.nih.gov/16840338/)
- Smith AM et al. Influenza A virus infection kinetics: quantitative data and models. *Wiley Interdiscip Rev Syst Biol Med* 2011;3:429-45. [PMID 21197654](https://pubmed.ncbi.nlm.nih.gov/21197654/)
- Beauchemin CA et al. A review of mathematical models of influenza A infections within a host or cell culture: lessons learned and challenges ahead. *BMC Public Health* 2011;11 Suppl 1:S7. [PMID 21356136](https://pubmed.ncbi.nlm.nih.gov/21356136/)
- Canini L et al. Viral kinetic modeling: state of the art. *J Pharmacokinet Pharmacodyn* 2014;41:431-43. [PMID 24961742](https://pubmed.ncbi.nlm.nih.gov/24961742/)
- Canini L et al. Impact of different oseltamivir regimens on treating influenza A virus infection and resistance emergence: insights from a modelling study. *PLoS Comput Biol* 2014;10:e1003568. [PMID 24743564](https://pubmed.ncbi.nlm.nih.gov/24743564/)
- Sachak-Patwa R et al. A target-cell limited model can reproduce influenza infection dynamics in hosts with differing immune responses. *J Theor Biol* 2023;567:111491. [PMID 37044357](https://pubmed.ncbi.nlm.nih.gov/37044357/)
- Li K et al. Modelling the Effect of MUC1 on Influenza Virus Infection Kinetics and Macrophage Dynamics. *Viruses* 2021;13. [PMID 34066999](https://pubmed.ncbi.nlm.nih.gov/34066999/)
- Handel A et al. Neuraminidase inhibitor resistance in influenza: assessing the danger of its generation and spread. *PLoS Comput Biol* 2007;3:e240. [PMID 18069885](https://pubmed.ncbi.nlm.nih.gov/18069885/)

## 2. 인체 감염 실험과 자연 경과 · Human challenge studies and natural history

A0의 보정 목표(정점 역가·정점 시각·배출 기간·증상 발현 시각·IL-6 최고치)는 전부 이 절의 문헌에서 나왔다. Kaiser 2001과 Gentile 2001은 증상이 바이러스가 아니라 IL-6를 따른다는 이 모델의 SYM 방정식의 근거이다.

- Carrat F et al. Time lines of infection and disease in human influenza: a review of volunteer challenge studies. *Am J Epidemiol* 2008;167:775-85. [PMID 18230677](https://pubmed.ncbi.nlm.nih.gov/18230677/)
- Hayden FG et al. Local and systemic cytokine responses during experimental human influenza A virus infection. Relation to symptom formation and host defense. *J Clin Invest* 1998;101:643-9. [PMID 9449698](https://pubmed.ncbi.nlm.nih.gov/9449698/)
- Fritz RS et al. Nasal cytokine and chemokine responses in experimental influenza A virus infection: results of a placebo-controlled trial of intravenous zanamivir treatment. *J Infect Dis* 1999;180:586-93. [PMID 10438343](https://pubmed.ncbi.nlm.nih.gov/10438343/)
- Hayden FG et al. Use of the oral neuraminidase inhibitor oseltamivir in experimental human influenza: randomized controlled trials for prevention and treatment. *JAMA* 1999;282:1240-6. [PMID 10517426](https://pubmed.ncbi.nlm.nih.gov/10517426/)
- Calfee DP et al. Safety and efficacy of once daily intranasal zanamivir in preventing experimental human influenza A infection. *Antivir Ther* 1999;4:143-9. [PMID 12731753](https://pubmed.ncbi.nlm.nih.gov/12731753/)
- Kaiser L et al. Symptom pathogenesis during acute influenza: interleukin-6 and other cytokine responses. *J Med Virol* 2001;64:262-8. [PMID 11424113](https://pubmed.ncbi.nlm.nih.gov/11424113/)
- Gentile DA et al. Effect of intranasal challenge with interleukin-6 on upper airway symptomatology and physiology in allergic and nonallergic patients. *Ann Allergy Asthma Immunol* 2001;86:531-6. [PMID 11379804](https://pubmed.ncbi.nlm.nih.gov/11379804/)

## 3. 발록사비르 임상시험과 약리 · Baloxavir trials and pharmacology

A3 원장(ledger)의 대조 수치(TTAS 80.2 / 53.8 / 53.7 h, 배출 중단 96 / 72 / 24 h, day-2 역가 −1.3 / −2.8 / −4.8 log₁₀)는 CAPSTONE-1에서 왔다. A10의 예방 투여 논거는 BLOCKSTONE이다.

- Hayden FG et al. Baloxavir Marboxil for Uncomplicated Influenza in Adults and Adolescents. *N Engl J Med* 2018;379:913-923. [PMID 30184455](https://pubmed.ncbi.nlm.nih.gov/30184455/)
- Ison MG et al. Early treatment with baloxavir marboxil in high-risk adolescent and adult outpatients with uncomplicated influenza (CAPSTONE-2): a randomised, placebo-controlled, phase 3 trial. *Lancet Infect Dis* 2020;20:1204-1214. [PMID 32526195](https://pubmed.ncbi.nlm.nih.gov/32526195/)
- Ikematsu H et al. Baloxavir Marboxil for Prophylaxis against Influenza in Household Contacts. *N Engl J Med* 2020;383:309-320. [PMID 32640124](https://pubmed.ncbi.nlm.nih.gov/32640124/)
- Baker J et al. Baloxavir Marboxil Single-dose Treatment in Influenza-infected Children: A Randomized, Double-blind, Active Controlled Phase 3 Safety and Efficacy Trial (miniSTONE-2). *Pediatr Infect Dis J* 2020;39:700-705. [PMID 32516282](https://pubmed.ncbi.nlm.nih.gov/32516282/)
- Baker JB et al. Safety and Efficacy of Baloxavir Marboxil in Influenza-infected Children 5-11 Years of Age: A Post Hoc Analysis of a Phase 3 Study. *Pediatr Infect Dis J* 2023;42:983-989. [PMID 37595103](https://pubmed.ncbi.nlm.nih.gov/37595103/)
- Palmu S et al. A Phase 3 Safety and Efficacy Study of Baloxavir Marboxil in Children Less Than 1 Year Old With Suspected or Confirmed Influenza. *Pediatr Infect Dis J* 2025;44:645-649. [PMID 40279637](https://pubmed.ncbi.nlm.nih.gov/40279637/)
- Noshi T et al. In vitro characterization of baloxavir acid, a first-in-class cap-dependent endonuclease inhibitor of the influenza virus polymerase PA subunit. *Antiviral Res* 2018;160:109-117. [PMID 30316915](https://pubmed.ncbi.nlm.nih.gov/30316915/)
- Liu H et al. Simultaneous quantification of baloxavir marboxil and its active metabolite in human plasma using UHPLC-MS/MS: Application to a human pharmacokinetic study with different anticoagulants. *J Pharm Biomed Anal* 2024;249:116387. [PMID 39083919](https://pubmed.ncbi.nlm.nih.gov/39083919/)
- Ikematsu H et al. Comparative Effectiveness of Baloxavir Marboxil and Oseltamivir Treatment in Reducing Household Transmission of Influenza: A Post Hoc Analysis of the BLOCKSTONE Trial. *Influenza Other Respir Viruses* 2024;18:e13302. [PMID 38706384](https://pubmed.ncbi.nlm.nih.gov/38706384/)
- Ishiguro N et al. Clinical and virological outcomes with baloxavir compared with oseltamivir in pediatric patients aged 6 to < 12 years with influenza: an open-label, randomized, active-controlled trial protocol. *BMC Infect Dis* 2021;21:777. [PMID 34372769](https://pubmed.ncbi.nlm.nih.gov/34372769/)
- Ishiguro N et al. Clinical and Virologic Outcomes of Baloxavir Compared with Oseltamivir in Pediatric Patients with Influenza in Japan. *Infect Dis Ther* 2025;14:833-846. [PMID 40155497](https://pubmed.ncbi.nlm.nih.gov/40155497/)

## 4. 발록사비르 내성 (PA/I38T) · Baloxavir resistance

A6의 경쟁적 해방(competitive release) 분석에서 쓰인 EC50 50배 이동과 복제 적합도 비용 18%는 이 절에서 나왔다. CAPSTONE-1은 성인 발록사비르 투여군의 9.7%에서 PA/I38T·M·F 출현을 보고했다.

- Omoto S et al. Characterization of influenza virus variants induced by treatment with the endonuclease inhibitor baloxavir marboxil. *Sci Rep* 2018;8:9633. [PMID 29941893](https://pubmed.ncbi.nlm.nih.gov/29941893/)
- Uehara T et al. Treatment-Emergent Influenza Variant Viruses With Reduced Baloxavir Susceptibility: Impact on Clinical and Virologic Outcomes in Uncomplicated Influenza. *J Infect Dis* 2020;221:346-355. [PMID 31309975](https://pubmed.ncbi.nlm.nih.gov/31309975/)
- Lee LY et al. Evaluating the fitness of PA/I38T-substituted influenza A viruses with reduced baloxavir susceptibility in a competitive mixtures ferret model. *PLoS Pathog* 2021;17:e1009527. [PMID 33956888](https://pubmed.ncbi.nlm.nih.gov/33956888/)
- Kuroda T et al. In Vivo Antiviral Activity of Baloxavir against PA/I38T-Substituted Influenza A Viruses at Clinically Relevant Doses. *Viruses* 2023;15. [PMID 37243240](https://pubmed.ncbi.nlm.nih.gov/37243240/)
- Case JR et al. Impaired host shutoff is a fitness cost associated with baloxavir marboxil resistance mutations in influenza A virus PA/PA-X nuclease domain. *PLoS Pathog* 2026;22:e1013550. [PMID 41662358](https://pubmed.ncbi.nlm.nih.gov/41662358/)
- Sun B et al. Emergence of the novel PA-D27G mutation conferring reduced baloxavir susceptibility in influenza A viruses circulating in China, 2018-2025. *Emerg Microbes Infect* 2026;15:2620222. [PMID 41555529](https://pubmed.ncbi.nlm.nih.gov/41555529/)

## 5. 뉴라미니다제 억제제 · Neuraminidase inhibitors

오셀타미비르 카복실레이트의 PK(CL/F 18.8 L/h, ELF:plasma ≈ 1)와 임상 효과 크기의 논쟁(Cochrane 대 Dobson IPD 메타분석) 모두 이 절에 있다. 이 모델은 어느 쪽 편도 들지 않고 두 값을 원장에 나란히 적는다.

- Treanor JJ et al. Efficacy and safety of the oral neuraminidase inhibitor oseltamivir in treating acute influenza: a randomized controlled trial. US Oral Neuraminidase Study Group. *JAMA* 2000;283:1016-24. [PMID 10697061](https://pubmed.ncbi.nlm.nih.gov/10697061/)
- Nicholson KG et al. Efficacy and safety of oseltamivir in treatment of acute influenza: a randomised controlled trial. Neuraminidase Inhibitor Flu Treatment Investigator Group. *Lancet* 2000;355:1845-50. [PMID 10866439](https://pubmed.ncbi.nlm.nih.gov/10866439/)
- Dobson J et al. Oseltamivir treatment for influenza in adults: a meta-analysis of randomised controlled trials. *Lancet* 2015;385:1729-1737. [PMID 25640810](https://pubmed.ncbi.nlm.nih.gov/25640810/)
- Jefferson T et al. Oseltamivir for influenza in adults and children: systematic review of clinical study reports and summary of regulatory comments. *BMJ* 2014;348:g2545. [PMID 24811411](https://pubmed.ncbi.nlm.nih.gov/24811411/)
- Kaiser L et al. Impact of oseltamivir treatment on influenza-related lower respiratory tract complications and hospitalizations. *Arch Intern Med* 2003;163:1667-72. [PMID 12885681](https://pubmed.ncbi.nlm.nih.gov/12885681/)
- Kamal MA et al. Identification of new oral dosing regimens for the neuraminidase inhibitor oseltamivir in patients with moderate and severe renal impairment. *Clin Pharmacol Drug Dev* 2015;4:326-36. [PMID 27137141](https://pubmed.ncbi.nlm.nih.gov/27137141/)
- Su CP et al. Inhaled Zanamivir vs Oral Oseltamivir to Prevent Influenza-related Hospitalization or Death: A Nationwide Population-based Quasi-experimental Study. *Clin Infect Dis* 2022;75:1273-1279. [PMID 35299245](https://pubmed.ncbi.nlm.nih.gov/35299245/)
- Saisho Y et al. Pharmacokinetics and safety of intravenous peramivir, neuraminidase inhibitor of influenza virus, in healthy Japanese subjects. *Antivir Ther* 2017;22:313-323. [PMID 27805571](https://pubmed.ncbi.nlm.nih.gov/27805571/)
- Dong L et al. In Vitro and In Vivo Assessment of Pharmacokinetic Profile of Peramivir in the Context of Inhalation Therapy. *Pharmaceuticals (Basel)* 2025;18. [PMID 40005995](https://pubmed.ncbi.nlm.nih.gov/40005995/)
- Li L et al. Comparison of double-dose vs standard-dose oseltamivir in the treatment of influenza: A systematic review and meta-analysis. *J Clin Pharm Ther* 2020;45:918-926. [PMID 32497319](https://pubmed.ncbi.nlm.nih.gov/32497319/)
- Peng Y et al. Bioequivalence and Safety of Two Oseltamivir Phosphate for Oral Suspension in Healthy Chinese Subjects Under Fasting and Fed Conditions: A Randomized, Open‑Label, Single‑Dose, Crossover Study. *Clin Pharmacol Drug Dev* 2026;15:e1612. [PMID 41045035](https://pubmed.ncbi.nlm.nih.gov/41045035/)

## 6. NAI 내성과 감시 · NAI resistance and surveillance

H275Y 프로파일(EC50 300배, 적합도 비용 28%)의 근거이며, Baz 2009는 예방 투여 중 내성이 출현한 사례로서 A6의 '조기 투여가 임상 이익과 내성 선택을 동시에 산다'는 예측과 직접 맞닿아 있다.

- Baz M et al. Emergence of oseltamivir-resistant pandemic H1N1 virus during prophylaxis. *N Engl J Med* 2009;361:2296-7. [PMID 19907034](https://pubmed.ncbi.nlm.nih.gov/19907034/)
- Takashita E et al. Global update on the susceptibilities of human influenza viruses to neuraminidase inhibitors and the cap-dependent endonuclease inhibitor baloxavir, 2017-2018. *Antiviral Res* 2020;175:104718. [PMID 32004620](https://pubmed.ncbi.nlm.nih.gov/32004620/)
- Farrukee R et al. Characterization of Influenza B Virus Variants with Reduced Neuraminidase Inhibitor Susceptibility. *Antimicrob Agents Chemother* 2018;62. [PMID 30201817](https://pubmed.ncbi.nlm.nih.gov/30201817/)
- Treurnicht FK et al. Replacement of neuraminidase inhibitor-susceptible influenza A(H1N1) with resistant phenotype in 2008 and circulation of susceptible influenza A and B viruses during 2009-2013, South Africa. *Influenza Other Respir Viruses* 2019;13:54-63. [PMID 30218485](https://pubmed.ncbi.nlm.nih.gov/30218485/)
- Duwe SC et al. Increase of Synergistic Secondary Antiviral Mutations in the Evolution of A(H1N1)pdm09 Influenza Virus Neuraminidases. *Viruses* 2024;16. [PMID 39066271](https://pubmed.ncbi.nlm.nih.gov/39066271/)
- Oh DY et al. Preparing for the Next Influenza Season: Monitoring the Emergence and Spread of Antiviral Resistance. *Infect Drug Resist* 2023;16:949-959. [PMID 36814825](https://pubmed.ncbi.nlm.nih.gov/36814825/)
- Ait-Aissa A et al. Surveillance for antiviral resistance among influenza viruses circulating in Algeria during five consecutive influenza seasons (2009-2014). *J Med Virol* 2018;90:844-853. [PMID 29315673](https://pubmed.ncbi.nlm.nih.gov/29315673/)
- Ngiam JN et al. Antiviral Resistance in Influenza: Clinical and Public Health Implications. *Influenza Other Respir Viruses* 2026;20:e70269. [PMID 42306964](https://pubmed.ncbi.nlm.nih.gov/42306964/)
- Cochin M et al. Characterization of oseltamivir-resistant A(H5N1) clade 2.3.4.4b, genotype D1.1 variants identified in poultry farms of British Columbia, Canada. *Emerg Microbes Infect* 2026;15:2686474. [PMID 42419257](https://pubmed.ncbi.nlm.nih.gov/42419257/)

## 7. 선천 면역과 인터페론 · Innate immunity and interferon

T → R (refractory) 전환, NS1의 IFN 길항, 그리고 A11에서 세균 이차감염의 두 조건 중 하나인 '인터페론에 의한 항세균 방어 억제'의 근거.

- Shahangian A et al. Type I IFNs mediate development of postinfluenza bacterial pneumonia in mice. *J Clin Invest* 2009;119:1910-20. [PMID 19487810](https://pubmed.ncbi.nlm.nih.gov/19487810/)
- Cheung PH et al. Virus subtype-specific suppression of MAVS aggregation and activation by PB1-F2 protein of influenza A (H7N9) virus. *PLoS Pathog* 2020;16:e1008611. [PMID 32511263](https://pubmed.ncbi.nlm.nih.gov/32511263/)
- Ihazmade H et al. Impact of IFITM3 rs12252 and IFNAR2 rs2236757 SNP polymorphisms on influenza virus infection outcomes in Moroccan patients. *Mol Biol Rep* 2025;53:2. [PMID 41128776](https://pubmed.ncbi.nlm.nih.gov/41128776/)
- An S et al. Initial Influenza Virus Replication Can Be Limited in Allergic Asthma Through Rapid Induction of Type III Interferons in Respiratory Epithelium. *Front Immunol* 2018;9:986. [PMID 29867963](https://pubmed.ncbi.nlm.nih.gov/29867963/)
- Blümke J et al. Interferon lambda in anti-viral defense and cancer: dual roles, mechanism and therapeutic potential. *J Transl Med* 2026;24. [PMID 42092929](https://pubmed.ncbi.nlm.nih.gov/42092929/)

## 8. 적응 면역 · Adaptive immunity

CD8 확장의 시점과 크기(정점 d6–8, 나이브 대비 20–100배), 점막 IgA에 의한 진입 차단, 항체에 의한 비리온 제거 증강의 근거.

- Saidu A et al. Highly focused human CD8+ T-cell response in the lower airways during acute influenza infection. *J Immunol* 2026;215. [PMID 42152613](https://pubmed.ncbi.nlm.nih.gov/42152613/)
- Hertoghs N et al. A group 1 hemagglutinin stem vaccine elicits broad humoral responses against influenza in phase 1/2a study. *Nat Commun* 2026;17. [PMID 41974675](https://pubmed.ncbi.nlm.nih.gov/41974675/)
- Vanderven HA et al. Fc-mediated functions and the treatment of severe respiratory viral infections with passive immunotherapy - a balancing act. *Front Immunol* 2023;14:1307398. [PMID 38077353](https://pubmed.ncbi.nlm.nih.gov/38077353/)
- Papargyris L et al. Innate immune responsiveness predicts enhanced cellular immunity and symptomatic disease after controlled human influenza infection. *Nat Med* 2026;32:2556-2569. [PMID 42387215](https://pubmed.ncbi.nlm.nih.gov/42387215/)

## 9. 이차 세균 감염 · Secondary bacterial infection

A11의 BAC 구획은 두 조건의 곱으로 만들어졌다 — 탈락한 상피의 노출된 부착 부위, 그리고 I형 인터페론에 의한 항세균 방어 억제. 임계값 자체는 보정되지 않은 가정이며 README에 그렇게 명시했다.

- McCullers JA Insights into the interaction between influenza virus and pneumococcus. *Clin Microbiol Rev* 2006;19:571-82. [PMID 16847087](https://pubmed.ncbi.nlm.nih.gov/16847087/)
- Shahangian A et al. Type I IFNs mediate development of postinfluenza bacterial pneumonia in mice. *J Clin Invest* 2009;119:1910-20. [PMID 19487810](https://pubmed.ncbi.nlm.nih.gov/19487810/)
- Kato K et al. Influenza A virus infection induces initial proliferation of commensal Streptococcus pneumoniae in the larynx leading to dissemination into the lower respiratory tract. *J Virol* 2026;100:e0055526. [PMID 42370676](https://pubmed.ncbi.nlm.nih.gov/42370676/)
- Seo JH et al. Current insights into bacterial secondary infection following influenza A virus infection. *Front Microbiol* 2026;17:1851115. [PMID 42267108](https://pubmed.ncbi.nlm.nih.gov/42267108/)
- Malainou C et al. TNF superfamily member 14 drives post-influenza depletion of alveolar macrophages, enabling secondary pneumococcal pneumonia. *J Clin Invest* 2026;136. [PMID 41252214](https://pubmed.ncbi.nlm.nih.gov/41252214/)
- Alshammari AK et al. Understanding the Molecular Interactions Between Influenza A Virus and Streptococcus Proteins in Co-Infection: A Scoping Review. *Pathogens* 2025;14. [PMID 40005491](https://pubmed.ncbi.nlm.nih.gov/40005491/)
- Zhang H Concerns of using sialidase fusion protein as an experimental drug to combat seasonal and pandemic influenza. *J Antimicrob Chemother* 2008;62:219-23. [PMID 18238888](https://pubmed.ncbi.nlm.nih.gov/18238888/)

## 10. 중증 인플루엔자 · 면역저하 숙주 · 보조 요법

A9(코르티코스테로이드)과 A12(면역저하 숙주)의 근거. 스테로이드의 관찰연구 사망률 신호는 이 모델에서 사망률 항 없이 기전만으로 방향이 재현된다.

- Lansbury L et al. Corticosteroids as adjunctive therapy in the treatment of influenza. *Cochrane Database Syst Rev* 2019;2:CD010406. [PMID 30798570](https://pubmed.ncbi.nlm.nih.gov/30798570/)
- Lansbury LE et al. Corticosteroids as Adjunctive Therapy in the Treatment of Influenza: An Updated Cochrane Systematic Review and Meta-analysis. *Crit Care Med* 2020;48:e98-e106. [PMID 31939808](https://pubmed.ncbi.nlm.nih.gov/31939808/)
- Muthuri SG et al. Effectiveness of neuraminidase inhibitors in reducing mortality in patients admitted to hospital with influenza A H1N1pdm09 virus infection: a meta-analysis of individual participant data. *Lancet Respir Med* 2014;2:395-404. [PMID 24815805](https://pubmed.ncbi.nlm.nih.gov/24815805/)
- Heldman MR et al. Influenza Antivirals for Prevention and Treatment in Immunocompromised People. *J Infect Dis* 2025;232:S243-S253. [PMID 41102612](https://pubmed.ncbi.nlm.nih.gov/41102612/)
- Hui DSC et al. Host Immunomodulatory Interventions in Severe Influenza. *J Infect Dis* 2025;232:S262-S272. [PMID 41102617](https://pubmed.ncbi.nlm.nih.gov/41102617/)
- Salvatore M et al. Baloxavir for the treatment of Influenza in allogeneic hematopoietic stem cell transplant recipients previously treated with oseltamivir. *Transpl Infect Dis* 2020;22:e13336. [PMID 32449254](https://pubmed.ncbi.nlm.nih.gov/32449254/)
- Euzen V et al. Zanamivir and baloxavir combination to cure persistent influenza and coronavirus infections after hematopoietic stem cell transplant. *Int J Antimicrob Agents* 2024;64:107281. [PMID 39047913](https://pubmed.ncbi.nlm.nih.gov/39047913/)
- Vetter P et al. Use of baloxavir as adjunctive antiviral therapy to neuraminidase inhibitors in severely immunocompromised individuals infected with influenza. *Antimicrob Agents Chemother* 2026;70:e0165925. [PMID 41891860](https://pubmed.ncbi.nlm.nih.gov/41891860/)
- Margaroli C et al. Spatial mapping of SARS-CoV-2 and H1N1 lung injury identifies differential transcriptional signatures. *Cell Rep Med* 2021;2:100242. [PMID 33778787](https://pubmed.ncbi.nlm.nih.gov/33778787/)

## 11. 그 밖의 항바이러스제 · Other antivirals and biologics

파비피라비르(치사적 돌연변이 유발), 페라미비르, M2 억제제(S31N으로 사실상 폐기), 항-HA 단클론항체의 근거. 이 중 오셀타미비르와 발록사비르만이 임상자료에 보정되었고 나머지는 가정값이다.

- Stevaert A et al. Nucleoside analogs for management of respiratory virus infections: mechanism of action and clinical efficacy. *Curr Opin Virol* 2022;57:101279. [PMID 36403338](https://pubmed.ncbi.nlm.nih.gov/36403338/)
- Geraghty RJ et al. Broad-Spectrum Antiviral Strategies and Nucleoside Analogues. *Viruses* 2021;13. [PMID 33924302](https://pubmed.ncbi.nlm.nih.gov/33924302/)
- Gregor J et al. Structural and Thermodynamic Analysis of the Resistance Development to Pimodivir (VX-787), the Clinical Inhibitor of Cap Binding to PB2 Subunit of Influenza A Polymerase. *Molecules* 2021;26. [PMID 33673017](https://pubmed.ncbi.nlm.nih.gov/33673017/)
- Stampolaki M et al. Adamantane-based inhibitors of the influenza A M2 proton channel: structure-based design, biological evaluation, and synthetic approaches. *RSC Med Chem* 2026;17:1811-1846. [PMID 41908668](https://pubmed.ncbi.nlm.nih.gov/41908668/)
- Dong J et al. Discovery of adamantane-based α-hydroxycarboxylic acid derivatives as potent M2-S31N blockers of influenza A virus. *Bioorg Chem* 2025;166:109124. [PMID 41151327](https://pubmed.ncbi.nlm.nih.gov/41151327/)
- Lim JJ et al. A Phase 2 Randomized, Double-Blind, Placebo-Controlled Trial of the Monoclonal Antibody MHAA4549A in Patients With Acute Uncomplicated Influenza A Infection. *Open Forum Infect Dis* 2022;9:ofab630. [PMID 35106315](https://pubmed.ncbi.nlm.nih.gov/35106315/)
- Mota KG et al. Monoclonal antibodies against influenza viruses: a clinical trials review. *Front Immunol* 2025;16:1669073. [PMID 41142796](https://pubmed.ncbi.nlm.nih.gov/41142796/)
- Wei J et al. Inhibition of cap-dependent endonuclease in influenza virus with ADC189: a pre-clinical analysis and phase I trial. *Front Med* 2025;19:347-358. [PMID 39832023](https://pubmed.ncbi.nlm.nih.gov/39832023/)

## 12. 평가변수 · 역학 · 상피 재생

TTAS의 정의(7개 증상이 모두 경증 이하로 21.5시간 지속)와 그 계측 도구, 그리고 상피 재생(LREG)의 근거.

- Keeley TJH et al. Content validity and psychometric properties of the inFLUenza Patient-Reported Outcome Plus (FLU-PRO Plus(©)) instrument in patients with COVID-19. *Qual Life Res* 2023;32:1645-1657. [PMID 36703019](https://pubmed.ncbi.nlm.nih.gov/36703019/)
- Saito R et al. Duration of fever and symptoms in children after treatment with baloxavir marboxil and oseltamivir during the 2018-2019 season and detection of variant influenza a viruses with polymerase acidic subunit substitutions. *Antiviral Res* 2020;183:104951. [PMID 32987032](https://pubmed.ncbi.nlm.nih.gov/32987032/)
- Chung JR et al. Influenza vaccine effectiveness against outpatient acute respiratory illness with laboratory-confirmed influenza, United States, 2024-25 season. *Clin Infect Dis* 2026. [PMID 42442752](https://pubmed.ncbi.nlm.nih.gov/42442752/)
- Okoli GN et al. A systematic meta-analytic comparative evaluation of seasonal influenza vaccine effectiveness from test-negative design studies in the Northern Hemisphere pre/post COVID-19 pandemic. *Vaccine* 2026;88:128956. [PMID 42486051](https://pubmed.ncbi.nlm.nih.gov/42486051/)
- Chan LYH et al. Understanding spatiotemporal clustering of seasonal influenza in the United States. *BMC Infect Dis* 2026;26. [PMID 41781891](https://pubmed.ncbi.nlm.nih.gov/41781891/)
- Brcko IC et al. Phylodynamic reconstruction of H1N1pdm09 influenza virus transmission in Brazil: a decade of evolutionary dynamics. *Emerg Microbes Infect* 2026;15:2620237. [PMID 41555524](https://pubmed.ncbi.nlm.nih.gov/41555524/)
- Grieco L et al. Exploring the role of mass immunisation in influenza pandemic preparedness: A modelling study for the UK context. *Vaccine* 2020;38:5163-5170. [PMID 32576461](https://pubmed.ncbi.nlm.nih.gov/32576461/)
- Gagneux P et al. Human-specific regulation of alpha 2-6-linked sialic acids. *J Biol Chem* 2003;278:48245-50. [PMID 14500706](https://pubmed.ncbi.nlm.nih.gov/14500706/)
- Zeng H et al. Tropism and Infectivity of a Seasonal A(H1N1) and a Highly Pathogenic Avian A(H5N1) Influenza Virus in Primary Differentiated Ferret Nasal Epithelial Cell Cultures. *J Virol* 2019;93. [PMID 30814288](https://pubmed.ncbi.nlm.nih.gov/30814288/)
- Weiner AI et al. ΔNp63 drives dysplastic alveolar remodeling and restricts epithelial plasticity upon severe lung injury. *Cell Rep* 2022;41:111805. [PMID 36516758](https://pubmed.ncbi.nlm.nih.gov/36516758/)
- Lu T et al. Dysplastic epithelial repair promotes the tissue residence of lymphocytes to inhibit alveolar regeneration post viral infection. *Cell Stem Cell* 2026;33:108-124.e6. [PMID 41443194](https://pubmed.ncbi.nlm.nih.gov/41443194/)
- Yu Y et al. Timing-specific efficacy of antiviral baloxavir and anti-inflammatory oclacitinib monotherapies, and the benefits of their combination in treating influenza in mice. *Front Microbiol* 2026;17:1741128. [PMID 41918532](https://pubmed.ncbi.nlm.nih.gov/41918532/)


---

## 이 모델이 재현하지 못하는 것 (Reported discrepancies)

문헌과 어긋나는 지점은 감추지 않고 분석 안에서 번호를 붙여 보고한다. 전부
`flu_reference_check.py`의 출력에 그대로 나온다.

| # | 어디서 | 내용 |
|---|---|---|
| 1 | A2 | Hill 기울기 1(Michaelis)로는 CAPSTONE-1의 24시간 역가 하강(−4.8 log₁₀)을 **어떤 Emax로도** 만들 수 없다. 잔여 생산 분율의 하한이 EC50/(C+EC50)이기 때문이다. 시험 자료는 효력이 아니라 농도-반응 **곡선의 모양**에 대한 증거다. |
| 2 | A12 | 면역저하 숙주에서 늦은 항바이러스 치료가 가치를 유지한다는 예측 — 임상 관행과는 맞지만 이를 검증한 무작위 시험이 없다. 반증되지 않았을 뿐 검증된 것이 아니다. |
| 3 | A8 | 증상이 바이러스를 조금이라도 따라간다면 발록사비르가 오셀타미비르보다 증상 이익이 커야 한다. CAPSTONE-1은 day-2 역가가 2 log 다른데도 53.7 대 53.8시간으로 구별하지 못했다. |
| 4 | A6 | 결정론적 모델이므로 10⁻⁴ 수준의 내성 계통이 모든 가상 환자에게 '존재'한다. 실제 I38T 출현은 성인의 약 10%에서 일어나는 확률적 정착 사건이다. |
| 5 | A7 | **순서가 거꾸로 나온다.** 이 모델에서는 기존 면역이 야생형을 늦춰 표적세포 밭을 더 많이 남기므로 내성 선택이 오히려 늘어난다. 소아에서 I38T가 더 흔하다는 관찰과 반대 방향이며, 잘 섞인(well-mixed) 상피 가정이 이 질문에는 틀렸다는 가장 분명한 신호다. |

---

## 도구 · 방법론 (Tools)

- mrgsolve — <https://mrgsolve.org/>
- Graphviz — <https://graphviz.org/>
- PubMed E-utilities (이 목록의 검증에 사용) — <https://www.ncbi.nlm.nih.gov/books/NBK25501/>

**총 97편** (전부 PubMed 조회로 확인).

> 교육·연구 목적의 모델입니다. 임상 의사결정에 사용해서는 안 됩니다.
