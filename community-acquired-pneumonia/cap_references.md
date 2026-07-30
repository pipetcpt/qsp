# 지역사회획득 폐렴 (Community-Acquired Pneumonia, CAP) — 참고문헌

QSP 모델(`cap_qsp_model.dot`, `cap_mrgsolve_model.R`, `cap_shiny_app.R`)의 구조·파라미터·
검증 앵커에 사용된 문헌 목록입니다. 모든 PMID는 PubMed E-utilities로 제목을 대조하여
확인했습니다.

각 섹션 끝의 **▶ 모델 연결** 항목은 그 문헌이 모델의 어느 방정식/파라미터를 뒷받침하는지
명시합니다.

---

## 1. 역학 · 질병 부담 · 총론 (Epidemiology, Burden, Overview)

1. Torres A, Cilloniz C, Niederman MS, et al. **Pneumonia.** *Nat Rev Dis Primers.* 2021;7:25.
   <https://pubmed.ncbi.nlm.nih.gov/33833230/>
2. Musher DM, Thorner AR. **Community-acquired pneumonia.** *N Engl J Med.* 2014;371:1619-1628.
   <https://pubmed.ncbi.nlm.nih.gov/25337751/>
3. Wunderink RG, Waterer GW. **Clinical practice. Community-acquired pneumonia.**
   *N Engl J Med.* 2014;370:543-551. <https://pubmed.ncbi.nlm.nih.gov/24499212/>
4. Prina E, Ranzani OT, Torres A. **Community-acquired pneumonia.** *Lancet.* 2015;386:1097-1108.
   <https://pubmed.ncbi.nlm.nih.gov/26277247/>
5. Jain S, Self WH, Wunderink RG, et al. **Community-acquired pneumonia requiring
   hospitalization among U.S. adults (EPIC study).** *N Engl J Med.* 2015;373:415-427.
   <https://pubmed.ncbi.nlm.nih.gov/26172429/>
6. Welte T, Torres A, Nathwani D. **Clinical and economic burden of community-acquired
   pneumonia among adults in Europe.** *Thorax.* 2012;67:71-79.
   <https://pubmed.ncbi.nlm.nih.gov/20729232/>
7. Ramirez JA, Wiemken TL, Peyrani P, et al. **Adults hospitalized with pneumonia in the
   United States: incidence, epidemiology, and mortality.** *Clin Infect Dis.* 2017;65:1806-1812.
   <https://pubmed.ncbi.nlm.nih.gov/29020164/>
8. Mizgerd JP. **Lung infection — a public health priority.** *PLoS Med.* 2006;3:e76.
   <https://pubmed.ncbi.nlm.nih.gov/16401173/>

**▶ 모델 연결**: 초기 접종량 `B0`·중증도 `SEV0`의 분포 범위, 그리고 `HAZ0`/`HAZ_AGE`
기저 위험도 스케일 설정.

---

## 2. 진료지침 · 중증도 평가 (Guidelines & Severity Scoring)

9. Metlay JP, Waterer GW, Long AC, et al. **Diagnosis and treatment of adults with
   community-acquired pneumonia. An official ATS/IDSA clinical practice guideline.**
   *Am J Respir Crit Care Med.* 2019;200:e45-e67.
   <https://pubmed.ncbi.nlm.nih.gov/31573350/>
10. Mandell LA, Wunderink RG, Anzueto A, et al. **IDSA/ATS consensus guidelines on the
    management of community-acquired pneumonia in adults.** *Clin Infect Dis.* 2007;44 Suppl 2:S27-72.
    <https://pubmed.ncbi.nlm.nih.gov/17278083/>
11. Lim WS, van der Eerden MM, Laing R, et al. **Defining community acquired pneumonia
    severity on presentation to hospital (CURB-65).** *Thorax.* 2003;58:377-382.
    <https://pubmed.ncbi.nlm.nih.gov/12728155/>
12. Fine MJ, Auble TE, Yealy DM, et al. **A prediction rule to identify low-risk patients
    with community-acquired pneumonia (PSI).** *N Engl J Med.* 1997;336:243-250.
    <https://pubmed.ncbi.nlm.nih.gov/8995086/>
13. Singer M, Deutschman CS, Seymour CW, et al. **The Third International Consensus
    Definitions for Sepsis and Septic Shock (Sepsis-3).** *JAMA.* 2016;315:801-810.
    <https://pubmed.ncbi.nlm.nih.gov/26903338/>
14. Evans L, Rhodes A, Alhazzani W, et al. **Surviving Sepsis Campaign: international
    guidelines for management of sepsis and septic shock 2021.** *Crit Care Med.* 2021;49:e1063-e1143.
    <https://pubmed.ncbi.nlm.nih.gov/34643578/>
15. Salluh JIF, Souza-Dantas VC, Martin-Loeches I, et al. **Challenges for a broad
    international implementation of the current severe community-acquired pneumonia
    guidelines.** *Intensive Care Med.* 2024;50:526-538.
    <https://pubmed.ncbi.nlm.nih.gov/38546855/>

**▶ 모델 연결**: `[TABLE]`의 CURB-65 구성요소(의식·요소질소 대용 `AKID`·호흡수·혈압·연령)와
SOFA 하위점수 컷오프, Halm 임상안정 기준(`Clinically_stable`)을 그대로 구현.

---

## 3. 병원체 생물학 · 병독인자 (Pathogen Biology & Virulence)

16. Weiser JN, Ferreira DM, Paton JC. ***Streptococcus pneumoniae*: transmission,
    colonization and invasion.** *Nat Rev Microbiol.* 2018;16:355-367.
    <https://pubmed.ncbi.nlm.nih.gov/29599457/>
17. Kadioglu A, Weiser JN, Paton JC, Andrew PW. **The role of *Streptococcus pneumoniae*
    virulence factors in host respiratory colonization and disease.** *Nat Rev Microbiol.*
    2008;6:288-301. <https://pubmed.ncbi.nlm.nih.gov/18340341/>
18. Loughran AJ, Orihuela CJ, Tuomanen EI. ***Streptococcus pneumoniae*: invasion and
    inflammation.** *Microbiol Spectr.* 2019;7. <https://pubmed.ncbi.nlm.nih.gov/30873934/>
19. Marriott HM, Mitchell TJ, Dockrell DH. **Pneumolysin: a double-edged sword during the
    host-pathogen interaction.** *Curr Mol Med.* 2008;8:497-509.
    <https://pubmed.ncbi.nlm.nih.gov/18781957/>
20. Bewley MA, Naughton M, Preston J, et al. **Pneumolysin activates macrophage lysosomal
    membrane permeabilization and executes apoptosis by distinct mechanisms.** *mBio.* 2014;5:e01710-14.
    <https://pubmed.ncbi.nlm.nih.gov/25293758/>
21. Anderson R, Nel JG, Feldman C. **Multifaceted role of pneumolysin in the pathogenesis
    of myocardial injury in community-acquired pneumonia.** *Int J Mol Sci.* 2018;19:1147.
    <https://pubmed.ncbi.nlm.nih.gov/29641429/>
22. Mizgerd JP. **Pathogenesis of severe pneumonia: advances and knowledge gaps.**
    *Curr Opin Pulm Med.* 2017;23:193-197. <https://pubmed.ncbi.nlm.nih.gov/28221171/>
23. McCullers JA. **Insights into the interaction between influenza virus and pneumococcus.**
    *Clin Microbiol Rev.* 2006;19:571-582. <https://pubmed.ncbi.nlm.nih.gov/16847087/>

**▶ 모델 연결**: 협막에 의한 탐식 포화(`KM_PHAG`), 뉴몰라이신의 상피 손상·TLR4/NLRP3 자극
(`PLYFRAC`), 인플루엔자 NanA 상승효과(`EVIR_ADH`, `EVIR_BAR`).

---

## 4. 용균-결합 PAMP 방출 — 이 모델의 핵심 항 (Lysis-Coupled PAMP Liberation)

24. Nau R, Eiffert H. **Modulation of release of proinflammatory bacterial compounds by
    antibacterials: potential impact on course and outcome of bacterial infections.**
    *Clin Microbiol Rev.* 2002;15:95-110. <https://pubmed.ncbi.nlm.nih.gov/11781269/>
25. Spreer A, Kerstan H, Böttcher T, et al. **Reduced release of pneumolysin by
    *Streptococcus pneumoniae* in vitro and in vivo after treatment with nonbacteriolytic
    antibiotics.** *Antimicrob Agents Chemother.* 2003;47:2649-2654.
    <https://pubmed.ncbi.nlm.nih.gov/12878534/>
26. Lepper PM, Held TK, Schneider EM, et al. **Clinical implications of antibiotic-induced
    endotoxin release in septic shock.** *Intensive Care Med.* 2002;28:824-833.
    <https://pubmed.ncbi.nlm.nih.gov/12122518/>

**▶ 모델 연결**: `dxdt_PAMP`의 `lysis_release = (YLYS_BL*kbl + YLYS_FQ*klvx + YLYS_MAC*kazi)*BE`
항 전체. β-락탐의 용균 수율을 1.0, 퀴놀론 0.35, 마크로라이드 0.12로 둔 근거이며,
첫 투여 후 6-24시간의 염증 급등과 일시적 임상 악화를 만드는 유일한 기전입니다.
`YLYS_BL = 0`으로 두면 이 현상이 사라지므로 반증 가능한 구조입니다.

---

## 5. 마크로라이드 면역조절 — 죽이지 않고 이득을 주는 경로 (Macrolide Immunomodulation)

27. Kovaleva A, Remmelts HHF, Rijkers GT, et al. **Immunomodulatory effects of macrolides
    during community-acquired pneumonia: a literature review.** *J Antimicrob Chemother.*
    2012;67:530-540. <https://pubmed.ncbi.nlm.nih.gov/22190609/>
28. Martin-Loeches I, Lisboa T, Rodriguez A, et al. **Combination antibiotic therapy with
    macrolides improves survival in intubated patients with community-acquired pneumonia.**
    *Intensive Care Med.* 2010;36:612-620. <https://pubmed.ncbi.nlm.nih.gov/19953222/>
29. Sligl WI, Asadi L, Eurich DT, et al. **Macrolides and mortality in critically ill
    patients with community-acquired pneumonia: a systematic review and meta-analysis.**
    *Crit Care Med.* 2014;42:420-432. <https://pubmed.ncbi.nlm.nih.gov/24158175/>
30. Restrepo MI, Anzueto A. **Macrolide therapy of pneumonia: is it necessary, and how does
    it help?** *Curr Opin Infect Dis.* 2016;29:212-217.
    <https://pubmed.ncbi.nlm.nih.gov/26836375/>

**▶ 모델 연결**: `MAC_ANTI`(사이토카인 증폭기 억제, `IMAX_MAC`), `MAC_PLY`(뉴몰라이신 억제,
`IMAX_PLY`), `MAC_MIG`(호중구 유주 억제, `IMAX_MIGM`) — 셋 모두 `kazi`(살균 경로)와
독립입니다. 따라서 `MIC_AZI = 64`(ermB 내성)로 두어도 이득이 대부분 남습니다.

---

## 6. 항생제 약동/약력학 (Antibiotic PK/PD)

31. Craig WA. **Pharmacokinetic/pharmacodynamic parameters: rationale for antibacterial
    dosing of mice and men.** *Clin Infect Dis.* 1998;26:1-10.
    <https://pubmed.ncbi.nlm.nih.gov/9455502/>
32. Drusano GL. **Antimicrobial pharmacodynamics: critical interactions of 'bug and drug'.**
    *Nat Rev Microbiol.* 2004;2:289-300. <https://pubmed.ncbi.nlm.nih.gov/15031728/>
33. Craig WA, Andes D. **Pharmacokinetics and pharmacodynamics of antibiotics in otitis
    media.** *Pediatr Infect Dis J.* 1996;15:255-259.
    <https://pubmed.ncbi.nlm.nih.gov/8852915/>
34. Preston SL, Drusano GL, Berman AL, et al. **Pharmacodynamics of levofloxacin: a new
    paradigm for early clinical trials.** *JAMA.* 1998;279:125-129.
    <https://pubmed.ncbi.nlm.nih.gov/9440662/>
35. Drwiega EN, Rodvold KA. **Penetration of antibacterial agents into pulmonary epithelial
    lining fluid: an update.** *Clin Pharmacokinet.* 2022;61:17-46.
    <https://pubmed.ncbi.nlm.nih.gov/34651282/>
36. Kiem S, Schentag JJ. **Interpretation of epithelial lining fluid concentrations of
    antibiotics.** *Infect Chemother.* 2014;46:219-225.
    <https://pubmed.ncbi.nlm.nih.gov/25566401/>
37. Nightingale CH. **Pharmacokinetics and pharmacodynamics of newer macrolides.**
    *Pediatr Infect Dis J.* 1997;16:438-443. <https://pubmed.ncbi.nlm.nih.gov/9109156/>
38. Rapp RP, Kuhn R. **Comparison of five beta-lactam antibiotics against common nosocomial
    pathogens using the time above MIC at different creatinine clearances.**
    *Pharmacotherapy.* 2000;20:184S-190S. <https://pubmed.ncbi.nlm.nih.gov/10809349/>

**▶ 모델 연결**: `EMAX_*`/`FEC50_*`/`H_*` Emax-Hill 구조와 ELF 이행률(`PEN_CEF` 0.40,
`PEN_AMX` 0.30, `PEN_LVX` 1.20, `KP_AZI_ELF` 22), 그리고 `TAM`(fT>MIC)·`AUCF`(fAUC/MIC)
누적 구획. **중요**: 두 지표는 파라미터가 아니라 적분 결과입니다. β-락탐과 퀴놀론의
차이는 Hill 계수(1.5 vs 2.6) 하나에서 나옵니다.

---

## 7. 코르티코스테로이드 — 양방향 효과 (Corticosteroids: The Two-Armed Intervention)

39. Dequin PF, Meziani F, Quenot JP, et al. **Hydrocortisone in severe community-acquired
    pneumonia (CAPE-COD).** *N Engl J Med.* 2023;388:1931-1941.
    <https://pubmed.ncbi.nlm.nih.gov/36942789/>
40. Torres A, Sibila O, Ferrer M, et al. **Effect of corticosteroids on treatment failure
    among hospitalized patients with severe community-acquired pneumonia and high
    inflammatory response.** *JAMA.* 2015;313:677-686.
    <https://pubmed.ncbi.nlm.nih.gov/25688779/>
41. Blum CA, Nigro N, Briel M, et al. **Adjunct prednisone therapy for patients with
    community-acquired pneumonia: a multicentre, double-blind, randomised, placebo-
    controlled trial.** *Lancet.* 2015;385:1511-1518.
    <https://pubmed.ncbi.nlm.nih.gov/25608756/>
42. Confalonieri M, Urbino R, Potena A, et al. **Hydrocortisone infusion for severe
    community-acquired pneumonia: a preliminary randomized study.** *Am J Respir Crit Care Med.*
    2005;171:242-248. <https://pubmed.ncbi.nlm.nih.gov/15557131/>
43. Meduri GU, Shih MC, Bridges L, et al. **Low-dose methylprednisolone treatment in
    critically ill patients with severe community-acquired pneumonia (ESCAPe).**
    *Intensive Care Med.* 2022;48:1009-1023. <https://pubmed.ncbi.nlm.nih.gov/35723686/>
44. Saleem N, Kulkarni A, Snow TAC, et al. **Effect of corticosteroids on mortality and
    clinical cure in community-acquired pneumonia: a systematic review, meta-analysis,
    and trial sequential analysis.** *Chest.* 2023;163:484-497.
    <https://pubmed.ncbi.nlm.nih.gov/36087797/>
45. Stern A, Skalsky K, Avni T, et al. **Corticosteroids for pneumonia.**
    *Cochrane Database Syst Rev.* 2017;12:CD007720.
    <https://pubmed.ncbi.nlm.nih.gov/29236286/>
46. Feldman C, Anderson R. **Corticosteroids in the adjunctive therapy of community-acquired
    pneumonia: an appraisal of recent meta-analyses of clinical trials.**
    *J Thorac Dis.* 2016;8:E162-171. <https://pubmed.ncbi.nlm.nih.gov/27076965/>

**▶ 모델 연결**: 단일 효과구획 `HCE`에서 두 방향이 갈립니다 — 보호 팔 `GC_ANTI`
(NF-κB transrepression, `IMAX_GC` 0.72)와 유해 팔 `GC_PHAG`(옵소닌 탐식 억제,
`IMAX_GCP` 0.45). 시나리오 6과 7은 스테로이드 용량이 동일하고 MIC만 다르며,
이때 순효과의 **부호가 뒤집히는** 것이 이 모델이 주장하는 바입니다.
`IMAX_GCP = 0`으로 두면 그 주장이 무너집니다.

---

## 8. 바이오마커와 치료기간 (Biomarkers & Treatment Duration)

47. Schuetz P, Wirz Y, Sager R, et al. **Effect of procalcitonin-guided antibiotic treatment
    on mortality in acute respiratory infections: a patient level meta-analysis.**
    *Lancet Infect Dis.* 2018;18:95-107. <https://pubmed.ncbi.nlm.nih.gov/29037960/>
48. Schuetz P, Christ-Crain M, Thomann R, et al. **Effect of procalcitonin-based guidelines
    vs standard guidelines on antibiotic use in lower respiratory tract infections
    (ProHOSP).** *JAMA.* 2009;302:1059-1066. <https://pubmed.ncbi.nlm.nih.gov/19738090/>
49. Chalmers JD, Singanayagam A, Hill AT. **C-reactive protein is an independent predictor
    of severity in community-acquired pneumonia.** *Am J Med.* 2008;121:219-225.
    <https://pubmed.ncbi.nlm.nih.gov/18328306/>
50. Uranga A, España PP, Bilbao A, et al. **Duration of antibiotic treatment in
    community-acquired pneumonia: a multicenter randomized clinical trial.**
    *JAMA Intern Med.* 2016;176:1257-1265. <https://pubmed.ncbi.nlm.nih.gov/27455166/>
51. Dinh A, Ropers J, Duran C, et al. **Discontinuing β-lactam treatment after 3 days for
    patients with community-acquired pneumonia in non-critical care wards (PTC): a
    double-blind, randomised, placebo-controlled, non-inferiority trial.**
    *Lancet.* 2021;397:1195-1203. <https://pubmed.ncbi.nlm.nih.gov/33773631/>
52. el Moussaoui R, de Borgie CAJM, van den Broek P, et al. **Effectiveness of discontinuing
    antibiotic treatment after three days versus eight days in mild to moderate-severe
    community acquired pneumonia.** *BMJ.* 2006;332:1355.
    <https://pubmed.ncbi.nlm.nih.gov/16763247/>

**▶ 모델 연결**: CRP(t½ 19 h, IL-6 구동)와 PCT(t½ 24 h, 세균 PAMP 구동 + 인터페론에 의한
억제 `IPCT_VIR`)를 서로 다른 구동원으로 분리한 이유. 짧은 치료기간의 안전성은 `BE`가
아니라 잠복 아집단 `BP`(`KPERS`, `FTOL_BL`)가 결정하도록 구조화했습니다.

---

## 9. 예후 · 치료시점 · 장기 결과 (Outcomes, Timing, Long-Term Sequelae)

53. Halm EA, Fine MJ, Marrie TJ, et al. **Time to clinical stability in patients hospitalized
    with community-acquired pneumonia: implications for practice guidelines.**
    *JAMA.* 1998;279:1452-1457. <https://pubmed.ncbi.nlm.nih.gov/9600479/>
54. Kumar A, Roberts D, Wood KE, et al. **Duration of hypotension before initiation of
    effective antimicrobial therapy is the critical determinant of survival in human septic
    shock.** *Crit Care Med.* 2006;34:1589-1596.
    <https://pubmed.ncbi.nlm.nih.gov/16625125/>
55. Houck PM, Bratzler DW. **Administration of first hospital antibiotics for
    community-acquired pneumonia: does timeliness affect outcomes?**
    *Curr Opin Infect Dis.* 2005;18:151-156. <https://pubmed.ncbi.nlm.nih.gov/15735420/>
56. Corrales-Medina VF, Alvarez KN, Weissfeld LA, et al. **Association between
    hospitalization for pneumonia and subsequent risk of cardiovascular disease.**
    *JAMA.* 2015;313:264-274. <https://pubmed.ncbi.nlm.nih.gov/25602997/>
57. Corrales-Medina VF, Taljaard M, Fine MJ, et al. **Risk stratification for cardiac
    complications in patients hospitalized for community-acquired pneumonia.**
    *Mayo Clin Proc.* 2014;89:60-68. <https://pubmed.ncbi.nlm.nih.gov/24388023/>
58. Rudd KE, Johnson SC, Agesa KM, et al. **Global, regional, and national sepsis incidence
    and mortality, 1990-2017: analysis for the Global Burden of Disease Study.**
    *Lancet.* 2020;395:200-211. <https://pubmed.ncbi.nlm.nih.gov/31954465/>
59. Piqueras M, et al. **Persistent systemic interleukin-6 elevation after community-acquired
    pneumonia is associated with one-year outcomes.** *Respir Res.* 2025.
    <https://pubmed.ncbi.nlm.nih.gov/41239435/>

**▶ 모델 연결**: `dxdt_CHZ` 누적위험 함수(SOFA·연령·균혈증 log CFU 가중)와 시나리오 9의
door-to-antibiotic 지연 기하학 — 지연 동안 누적된 위험은 이후 회복되지 않습니다.

---

## 10. 폐 손상 · 장벽 · 가스교환 · 해소 (Lung Injury, Barrier, Gas Exchange, Resolution)

60. Matthay MA, Zemans RL, Zimmerman GA, et al. **Acute respiratory distress syndrome.**
    *Nat Rev Dis Primers.* 2019;5:18. <https://pubmed.ncbi.nlm.nih.gov/30872586/>
61. Ware LB, Matthay MA. **Alveolar fluid clearance is impaired in the majority of patients
    with acute lung injury and the acute respiratory distress syndrome.**
    *Am J Respir Crit Care Med.* 2001;163:1376-1383.
    <https://pubmed.ncbi.nlm.nih.gov/11371404/>
62. Roux J, McNicholas CM, Carles M, et al. **IL-8 inhibits cAMP-stimulated alveolar
    epithelial fluid transport via a GRK2/PI3K-dependent mechanism.** *FASEB J.*
    2013;27:1095-1106. <https://pubmed.ncbi.nlm.nih.gov/23221335/>
63. Herold S, Mayer K, Lohmeyer J. **Acute lung injury: how macrophages orchestrate
    resolution of inflammation and tissue repair.** *Front Immunol.* 2011;2:65.
    <https://pubmed.ncbi.nlm.nih.gov/22566854/>
64. Serhan CN. **Pro-resolving lipid mediators are leads for resolution physiology.**
    *Nature.* 2014;510:92-101. <https://pubmed.ncbi.nlm.nih.gov/24899309/>
65. Aberdein JD, Cole J, Bewley MA, et al. **Alveolar macrophages in pulmonary host defence
    — the unrecognised role of apoptosis as a mechanism of intracellular bacterial killing.**
    *Clin Exp Immunol.* 2013;174:193-202. <https://pubmed.ncbi.nlm.nih.gov/23841514/>
66. Dockrell DH, Whyte MKB, Mitchell TJ. **Pneumococcal pneumonia: mechanisms of infection
    and resolution.** *Chest.* 2012;142:482-491.
    <https://pubmed.ncbi.nlm.nih.gov/22871758/>
67. Grudzinska F, et al. **Hospitalised older adults with community-acquired pneumonia and
    sepsis have dysregulated neutrophil function.** *Thorax.* 2025.
    <https://pubmed.ncbi.nlm.nih.gov/39689942/>
68. Mohanty T, Fisher J, Bakochi A, et al. **Neutrophil extracellular traps in the central
    nervous system hinder bacterial clearance during pneumococcal meningitis.**
    *Nat Commun.* 2019;10:1667. <https://pubmed.ncbi.nlm.nih.gov/30971685/>
69. Archer SL, Dunham-Snary KJ, et al. **Hypoxic pulmonary vasoconstriction: an important
    component of the homeostatic oxygen sensing system.** *Physiol Res.* 2024;73:S493-S510.
    <https://pubmed.ncbi.nlm.nih.gov/39589299/>

**▶ 모델 연결**: `IAFC_IL8`(IL-8에 의한 폐포액 청소 억제, ref 62), `KS_SPM`/`ESPM`
(효페로사이토시스 → SPM → 증폭기 감쇠, ref 63-64), `KPHAG_AM`가 호중구보다 높은 근거
(ref 65), 그리고 `HPV_MAX`/`IHPV_CYT`(저산소성 폐혈관수축과 그 염증성 둔화, ref 69).

---

## 11. 예방 (Prevention)

70. Bonten MJM, Huijts SM, Bolkenbaas M, et al. **Polysaccharide conjugate vaccine against
    pneumococcal pneumonia in adults (CAPiTA).** *N Engl J Med.* 2015;372:1114-1125.
    <https://pubmed.ncbi.nlm.nih.gov/25785969/>
71. Kobayashi M, Leidner AJ, Gierke R, et al. **Use of 21-valent pneumococcal conjugate
    vaccine among U.S. adults: recommendations of the Advisory Committee on Immunization
    Practices.** *MMWR Morb Mortal Wkly Rep.* 2024;73:793-798.
    <https://pubmed.ncbi.nlm.nih.gov/39264843/>

**▶ 모델 연결**: `VACC` 공변량 → `OPS_VACC` 1.45배 옵소닌 활성 상승.

---

## 부록: 파라미터 출처 요약표

| 파라미터군 | 값 | 출처 |
|---|---|---|
| Ceftriaxone CL 1.0 L/h, V1 8 L, fu 0.10, ELF 0.40 | 2 g q24h → Cmax ~150, trough ~10 mg/L | 35, 36, 38 |
| Amoxicillin ka 1.5/h, F 0.80, CL 20 L/h, fu 0.82 | 1 g q8h | 33, 35 |
| Azithromycin F 0.37, Vc 450 L, 심부 2200 L, ELF 계수 22 | 혈청 ~0.4 → ELF 1-3 mg/L | 37, 35 |
| Levofloxacin CL 8 L/h, V 90 L, fu 0.70, ELF 1.2 | 750 mg → AUC24 ~100, fAUC/MIC ~70 | 34, 35 |
| Hydrocortisone CL 18 L/h, V 35 L, ke0 0.2/h | 200 mg/일 지속주입 | 39, 42 |
| CRP t½ 19 h · PCT t½ 24 h | 바이오마커 감쇠 | 47-49 |
| 용균 수율 β-락탐 1.0 / FQ 0.35 / 마크로라이드 0.12 | 항생제 계열별 PAMP 방출 | 24-26 |
| 스테로이드 보호 팔 0.72 / 유해 팔 0.45 | CAPE-COD 방향성 재현 | 39-46 |
| 임상안정까지 중앙값 ~3일 (유효 치료 시) | 검증 목표 | 53 |
| 무치료 자연경과 "위기" 7-9일 | 검증 목표 | 2, 3 |

---

## 면책 (Disclaimer)

본 모델과 문헌 정리는 **교육 및 연구 목적**입니다. 파라미터는 공개 문헌에 기반한 근사치이며
환자 수준 데이터에 적합(fit)되지 않았습니다. 임상 의사결정, 처방, 규제 제출에 사용할 수 없습니다.
