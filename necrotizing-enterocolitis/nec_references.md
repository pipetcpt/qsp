# 신생아 괴사성 장염 (NEC) — 참고문헌
# Necrotising Enterocolitis — Reference List for the QSP Model

> **링크 방식에 대한 안내 / A note on how these links are built.**
> 아래 항목은 대부분 **제목 기반 PubMed 검색 링크**입니다. 존재하지 않는 PMID를
> 만들어 넣는 것보다, 제목으로 조회되는 링크가 항상 올바른 논문을 가리키기
> 때문입니다. 확신하는 소수의 landmark 논문에는 PMID를 함께 적었습니다.
>
> Most entries below are **title-resolved PubMed search links** rather than
> bare PMIDs. A title query always resolves to the correct paper, whereas a
> mis-remembered accession number silently points at the wrong one. PMIDs are
> given only for the handful of landmark papers where the number is certain.

각 절의 제목 뒤 대괄호는 **이 문헌이 모델의 어느 항을 지지하는지**를
가리킵니다 — 문헌 목록이 아니라 모델의 근거 대장으로 읽으십시오.
The bracket after each section heading names **which term of the model that
section supports**, so this reads as a provenance ledger rather than a booklist.

---

## 1. 정의 · 병기 · 역학 (Definition, staging, epidemiology)
*[모델 항: Bell 병기 판독 규칙, GA별 발생률 보정 목표]*

1. Bell MJ, Ternberg JL, Feigin RD, et al. **Neonatal necrotizing enterocolitis:
   therapeutic decisions based upon clinical staging.** *Ann Surg.* 1978;187(1):1-7.
   (PMID 413500) —
   <https://pubmed.ncbi.nlm.nih.gov/413500/>
   → 모델의 `BellStage` 판독 규칙(의심 / 확진 / 진행)의 원전.
2. Walsh MC, Kliegman RM. **Necrotizing enterocolitis: treatment based on staging
   criteria.** *Pediatr Clin North Am.* 1986. —
   <https://pubmed.ncbi.nlm.nih.gov/?term=Walsh+Kliegman+necrotizing+enterocolitis+treatment+based+on+staging+criteria>
3. Neu J, Walker WA. **Necrotizing enterocolitis.** *N Engl J Med.*
   2011;364(3):255-264. (PMID 21247316) —
   <https://pubmed.ncbi.nlm.nih.gov/21247316/>
4. Battersby C, Santhalingam T, Costeloe K, Modi N. **Incidence of neonatal
   necrotising enterocolitis in high-income countries: a systematic review.**
   *Arch Dis Child Fetal Neonatal Ed.* 2018. —
   <https://pubmed.ncbi.nlm.nih.gov/?term=Battersby+incidence+of+neonatal+necrotising+enterocolitis+in+high-income+countries>
   → GA 구간별 발생률 보정 목표(24-25주 12-15 %, 30-32주 1-3 %)의 근거.
5. Patel RM, Kandefer S, Walsh MC, et al. **Causes and timing of death in
   extremely premature infants from 2000 through 2011.** *N Engl J Med.* 2015. —
   <https://pubmed.ncbi.nlm.nih.gov/?term=Patel+causes+and+timing+of+death+in+extremely+premature+infants+2000+through+2011>
6. Berman L, Moss RL. **Necrotizing enterocolitis: an update.**
   *Semin Fetal Neonatal Med.* 2011. —
   <https://pubmed.ncbi.nlm.nih.gov/?term=Berman+Moss+necrotizing+enterocolitis+an+update>
7. Sharma R, Hudak ML. **A clinical perspective of necrotizing enterocolitis:
   past, present, and future.** *Clin Perinatol.* 2013. —
   <https://pubmed.ncbi.nlm.nih.gov/?term=Sharma+Hudak+clinical+perspective+necrotizing+enterocolitis+past+present+future>
8. Gephart SM, Spitzer AR, Effken JA, et al. **Discrimination of GutCheck-NEC:
   a clinical risk index for necrotizing enterocolitis.** *J Perinatol.* 2014. —
   <https://pubmed.ncbi.nlm.nih.gov/?term=Gephart+GutCheck+NEC+clinical+risk+index+necrotizing+enterocolitis>

---

## 2. 미숙 장 상피의 TLR4 과발현 (Immature-gut TLR4 over-expression)
*[모델 항: `TLR4expr = TLR4max·(1 − phiTLR·MAT)` — 왜 미숙할수록 같은 균에 더 크게 반응하는가]*

9. Nanthakumar N, Meng D, Goldstein AM, et al. **The mechanism of excessive
   intestinal inflammation in necrotizing enterocolitis: an immature innate
   immune response.** *PLoS One.* 2011. —
   <https://pubmed.ncbi.nlm.nih.gov/?term=Nanthakumar+mechanism+of+excessive+intestinal+inflammation+necrotizing+enterocolitis+immature+innate+immune+response>
   → 미숙 장에서 IκB 발현이 낮고 NF-κB 반응이 과도하다는 직접 근거.
10. Sodhi CP, Neal MD, Siggers R, et al. **Intestinal epithelial Toll-like
    receptor 4 regulates goblet cell development and is required for
    necrotizing enterocolitis in mice.** *Gastroenterology.* 2012. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Sodhi+intestinal+epithelial+Toll-like+receptor+4+regulates+goblet+cell+development+necrotizing+enterocolitis>
11. Leaphart CL, Cavallo J, Gribar SC, et al. **A critical role for TLR4 in the
    pathogenesis of necrotizing enterocolitis by modulating intestinal injury
    and repair.** *J Immunol.* 2007. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Leaphart+critical+role+for+TLR4+pathogenesis+necrotizing+enterocolitis+modulating+intestinal+injury+and+repair>
    → **증식 차단항 `block = 1/(1+(INJ/Ki)²)`** 의 근거: TLR4 신호가 손상만
    늘리는 것이 아니라 크립트 재생을 동시에 억제한다.
12. Hackam DJ, Sodhi CP. **Toll-like receptor-mediated intestinal inflammatory
    imbalance in the pathogenesis of necrotizing enterocolitis.**
    *Cell Mol Gastroenterol Hepatol.* 2018. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Hackam+Sodhi+Toll-like+receptor+mediated+intestinal+inflammatory+imbalance+pathogenesis+necrotizing+enterocolitis>
13. Niño DF, Sodhi CP, Hackam DJ. **Necrotizing enterocolitis: new insights into
    pathogenesis and mechanisms.** *Nat Rev Gastroenterol Hepatol.* 2016. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Nino+Sodhi+Hackam+necrotizing+enterocolitis+new+insights+into+pathogenesis+and+mechanisms>
14. Good M, Sodhi CP, Hackam DJ. **Evidence-based feeding strategies before and
    after the development of necrotizing enterocolitis.**
    *Expert Rev Clin Immunol.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Good+Sodhi+Hackam+evidence-based+feeding+strategies+before+and+after+development+necrotizing+enterocolitis>
15. Sampah MES, Hackam DJ. **Dysregulated mucosal immunity and associated
    pathogeneses in preterm neonates.** *Front Immunol.* 2020. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Sampah+Hackam+dysregulated+mucosal+immunity+associated+pathogeneses+preterm+neonates>
16. Lu P, Sodhi CP, Hackam DJ. **Toll-like receptor regulation of intestinal
    development and inflammation in the pathogenesis of necrotizing
    enterocolitis.** *Pathophysiology.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Lu+Sodhi+Hackam+Toll-like+receptor+regulation+intestinal+development+inflammation+necrotizing+enterocolitis>

---

## 3. 장 미생물 군집과 이상균총 (Microbiome and dysbiosis)
*[모델 항: `B`, `C`, `alphaX` 교차경쟁, `fpH`, 두 개의 안정 상태]*

17. Warner BB, Deych E, Zhou Y, et al. **Gut bacteria dysbiosis and necrotising
    enterocolitis in very low birthweight infants: a prospective case-control
    study.** *Lancet.* 2016;387(10031):1928-1936. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Warner+gut+bacteria+dysbiosis+and+necrotising+enterocolitis+very+low+birthweight+prospective+case-control>
    → NEC 이전에 Gammaproteobacteria(Enterobacteriaceae) 상대풍부도가 증가하고
    편성 혐기균이 감소한다 — 모델의 `B` 우세 상태에 해당.
18. Pammi M, Cope J, Tarr PI, et al. **Intestinal dysbiosis in preterm infants
    preceding necrotizing enterocolitis: a systematic review and
    meta-analysis.** *Microbiome.* 2017. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Pammi+intestinal+dysbiosis+preterm+infants+preceding+necrotizing+enterocolitis+systematic+review+meta-analysis>
19. Olm MR, Bhattacharya N, Crits-Christoph A, et al. **Necrotizing
    enterocolitis is preceded by increased gut bacterial replication, Klebsiella,
    and fimbriae-encoding bacteria.** *Sci Adv.* 2019. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Olm+necrotizing+enterocolitis+preceded+by+increased+gut+bacterial+replication+Klebsiella+fimbriae>
    → 모델이 `B`의 *증식률*(`muB`)과 부하를 함께 다루는 이유.
20. Stewart CJ, Ajami NJ, O'Brien JL, et al. **Temporal development of the gut
    microbiome in early childhood from the TEDDY study.** *Nature.* 2018. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Stewart+temporal+development+of+the+gut+microbiome+in+early+childhood+TEDDY>
21. Claud EC, Walker WA. **Hypothesis: inappropriate colonization of the
    premature intestine can cause neonatal necrotizing enterocolitis.**
    *FASEB J.* 2001. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Claud+Walker+hypothesis+inappropriate+colonization+premature+intestine+neonatal+necrotizing+enterocolitis>
22. Ward DV, Scholz M, Zolfo M, et al. **Metagenomic sequencing with
    strain-level resolution implicates uropathogenic E. coli in necrotizing
    enterocolitis and mortality in preterm infants.** *Cell Rep.* 2016. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ward+metagenomic+sequencing+strain-level+resolution+uropathogenic+coli+necrotizing+enterocolitis+preterm>
23. Heida FH, van Zoonen AGJF, Hulscher JBF, et al. **A necrotizing
    enterocolitis-associated gut microbiota is present in the meconium:
    results of a prospective study.** *Clin Infect Dis.* 2016. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Heida+necrotizing+enterocolitis+associated+gut+microbiota+present+in+the+meconium+prospective>
24. Berrington JE, Stewart CJ, Cummings SP, Embleton ND. **The neonatal
    bowel microbiome in health and infection.** *Curr Opin Infect Dis.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Berrington+neonatal+bowel+microbiome+in+health+and+infection>

---

## 4. 모유 · HMO · B. infantis (Human milk, HMOs, and the strain that can eat them)
*[모델 항: `HMO` → `Ktlr`(문턱), `muCH·binfant`(생태), `eEGF`(영양), `SIGA`]*

25. Lucas A, Cole TJ. **Breast milk and neonatal necrotising enterocolitis.**
    *Lancet.* 1990;336(8730):1519-1523. (PMID 1979363) —
    <https://pubmed.ncbi.nlm.nih.gov/1979363/>
    → 모유의 보호효과를 정량화한 최초의 대규모 연구.
26. Meinzen-Derr J, Poindexter B, Wrage L, et al. **Role of human milk intake
    and persistent organ failure on the risk of necrotizing enterocolitis or
    death.** *J Perinatol.* 2009. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Meinzen-Derr+role+of+human+milk+intake+risk+of+necrotizing+enterocolitis+or+death>
27. Quigley M, Embleton ND, McGuire W. **Formula versus donor breast milk for
    feeding preterm or low birth weight infants.** *Cochrane Database Syst Rev.*
    2019. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Quigley+Embleton+McGuire+formula+versus+donor+breast+milk+for+feeding+preterm+or+low+birth+weight+infants>
    → 모델의 `fDM` 아암(저온살균으로 sIgA·EGF는 잃고 HMO는 보존)의 근거.
28. Cristofalo EA, Schanler RJ, Blanco CL, et al. **Randomized trial of
    exclusive human milk versus preterm formula diets in extremely premature
    infants.** *J Pediatr.* 2013. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Cristofalo+randomized+trial+of+exclusive+human+milk+versus+preterm+formula+diets+extremely+premature>
29. Autran CA, Kellman BP, Kim JH, et al. **Human milk oligosaccharide
    composition predicts risk of necrotising enterocolitis in preterm
    infants.** *Gut.* 2018. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Autran+human+milk+oligosaccharide+composition+predicts+risk+of+necrotising+enterocolitis+preterm>
    → 특정 HMO(DSLNT) 농도가 NEC 위험과 연관 — 모델에서 HMO가 **문턱**을
    올리는 항(`hHMO`)으로 들어가는 근거.
30. Masi AC, Embleton ND, Lamb CA, et al. **Human milk oligosaccharide DSLNT
    and gut microbiome in preterm infants predicts necrotising
    enterocolitis.** *Gut.* 2021. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Masi+human+milk+oligosaccharide+DSLNT+and+gut+microbiome+preterm+infants+predicts+necrotising+enterocolitis>
    → HMO와 미생물 군집이 **함께** 위험을 결정한다는 관찰 — 모델이 모유를
    문턱-이동과 상태-이동 **양쪽**으로 넣은 이유.
31. Jantscher-Krenn E, Zherebtsov M, Nissan C, et al. **The human milk
    oligosaccharide disialyllacto-N-tetraose prevents necrotising enterocolitis
    in neonatal rats.** *Gut.* 2012. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Jantscher-Krenn+human+milk+oligosaccharide+disialyllacto-N-tetraose+prevents+necrotising+enterocolitis+neonatal+rats>
32. Good M, Sodhi CP, Yamaguchi Y, et al. **The human milk oligosaccharide
    2'-fucosyllactose attenuates the severity of experimental necrotising
    enterocolitis by enhancing mesenteric perfusion in the neonatal
    intestine.** *Br J Nutr.* 2016. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Good+human+milk+oligosaccharide+2-fucosyllactose+attenuates+severity+experimental+necrotising+enterocolitis+mesenteric+perfusion>
    → 2′-FL이 **장혈류**를 통해서도 작용한다는 근거 — 모델의 `PERF` 경로와 연결.
33. Bode L. **Human milk oligosaccharides: every baby needs a sugar mama.**
    *Glycobiology.* 2012. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Bode+human+milk+oligosaccharides+every+baby+needs+a+sugar+mama>
34. Underwood MA, German JB, Lebrilla CB, Mills DA. **Bifidobacterium longum
    subspecies infantis: champion colonizer of the infant gut.**
    *Pediatr Res.* 2015. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Underwood+Bifidobacterium+longum+subspecies+infantis+champion+colonizer+of+the+infant+gut>
    → 모델의 `binfant` 스위치(HMO를 실제로 대사할 수 있는 균주의 보유 여부)의
    근거. HMO가 있어도 그것을 쓸 수 있는 균이 없으면 그냥 배설된다.
35. Nolan LS, Rimer JM, Good M. **The role of human milk oligosaccharides and
    probiotics on the neonatal microbiome and risk of necrotizing
    enterocolitis.** *Nutrients.* 2020. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Nolan+role+of+human+milk+oligosaccharides+and+probiotics+neonatal+microbiome+risk+of+necrotizing+enterocolitis>
36. Watkins DJ, Besner GE. **The role of the intestinal microcirculation in
    necrotizing enterocolitis.** *Semin Pediatr Surg.* 2013. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Watkins+Besner+role+of+the+intestinal+microcirculation+in+necrotizing+enterocolitis>
37. Feng J, El-Assal ON, Besner GE. **Heparin-binding EGF-like growth factor
    (HB-EGF) and necrotizing enterocolitis.** *Semin Pediatr Surg.* 2005. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Feng+Besner+heparin-binding+EGF-like+growth+factor+HB-EGF+and+necrotizing+enterocolitis>
    → 모델의 `eEGF` 항(모유 성장인자가 크립트 증식을 올린다).

---

## 5. 프로바이오틱스 (Probiotics — the pure state-mover)
*[모델 항: `P`, `prob_dose`, `EmxPamp/gen/mtz` — 왜 항생제와 같이 주면 듣지 않는가]*

38. Sharif S, Meader N, Oddie SJ, Rojas-Reyes MX, McGuire W. **Probiotics to
    prevent necrotising enterocolitis in very preterm or very low birth weight
    infants.** *Cochrane Database Syst Rev.* 2023. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Sharif+probiotics+to+prevent+necrotising+enterocolitis+very+preterm+or+very+low+birth+weight+infants+Cochrane>
    → 모델의 프로바이오틱스 아암 효과 크기 보정 목표.
39. Jacobs SE, Tobin JM, Opie GF, et al. **Probiotic effects on late-onset
    sepsis in very preterm infants: a randomized controlled trial (ProPrems).**
    *Pediatrics.* 2013. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Jacobs+probiotic+effects+on+late-onset+sepsis+in+very+preterm+infants+randomized+controlled+trial>
40. Costeloe K, Hardy P, Juszczak E, Wilks M, Millar MR. **Bifidobacterium
    breve BBG-001 in very preterm infants: a randomised controlled phase 3
    trial (PiPS).** *Lancet.* 2016. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Costeloe+Bifidobacterium+breve+BBG-001+in+very+preterm+infants+randomised+controlled+phase+3+trial>
    → **음성 결과** — 모델에서 정착 실패(`gP` < `kwash` + `killP`)로 재현된다.
41. Dermyshi E, Wang Y, Yan C, et al. **The 'golden age' of probiotics: a
    systematic review and meta-analysis of randomized and observational studies
    in preterm infants.** *Neonatology.* 2017. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Dermyshi+golden+age+of+probiotics+systematic+review+meta-analysis+randomized+observational+studies+preterm+infants>
42. van den Akker CHP, van Goudoever JB, Shamir R, et al. **Probiotics and
    preterm infants: a position paper by the ESPGHAN Committee on Nutrition and
    the ESPGHAN Working Group for Probiotics and Prebiotics.**
    *J Pediatr Gastroenterol Nutr.* 2020. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=van+den+Akker+probiotics+and+preterm+infants+position+paper+ESPGHAN+Committee+on+Nutrition>

---

## 6. 경험적 항생제 · 노출 기간 (Empirical antibiotics and duration)
*[모델 항: `killB` vs `killC` 비대칭, `abx_days` — 부호가 바뀌는 지점]*

43. Cotten CM, Taylor S, Stoll B, et al. **Prolonged duration of initial
    empirical antibiotic treatment is associated with increased rates of
    necrotizing enterocolitis and death for extremely low birth weight
    infants.** *Pediatrics.* 2009;123(1):58-66. (PMID 19117861) —
    <https://pubmed.ncbi.nlm.nih.gov/19117861/>
    → 모델의 **부호 전환** 결과가 재현하려는 관찰. 항생제는 `B`를 즉시 줄이지만
    `C`를 잃는 대가는 지연되어 나타나므로 두 적분이 어딘가에서 교차한다.
44. Alexander VN, Northrup V, Bizzarro MJ. **Antibiotic exposure in the newborn
    intensive care unit and the risk of necrotizing enterocolitis.**
    *J Pediatr.* 2011. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Alexander+antibiotic+exposure+in+the+newborn+intensive+care+unit+and+the+risk+of+necrotizing+enterocolitis>
45. Greenwood C, Morrow AL, Lagomarcino AJ, et al. **Early empiric antibiotic
    use in preterm infants is associated with lower bacterial diversity and
    higher relative abundance of Enterobacter.** *J Pediatr.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Greenwood+early+empiric+antibiotic+use+in+preterm+infants+lower+bacterial+diversity+higher+relative+abundance+Enterobacter>
    → `killC` > `killB` 비대칭의 직접 근거.
46. Ting JY, Synnes A, Roberts A, et al. **Association between antibiotic use
    and neonatal mortality and morbidities in very low-birth-weight infants
    without culture-proven sepsis or necrotizing enterocolitis.**
    *JAMA Pediatr.* 2016. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ting+association+between+antibiotic+use+and+neonatal+mortality+and+morbidities+very+low-birth-weight+without+culture-proven+sepsis>
47. Rao SC, Srinivasjois R, Moon K. **One dose per day compared to multiple
    doses per day of gentamicin for treatment of suspected or proven sepsis in
    neonates.** *Cochrane Database Syst Rev.* 2016. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Rao+one+dose+per+day+compared+to+multiple+doses+per+day+of+gentamicin+neonates+Cochrane>
    → 모델의 젠타마이신 q36h 요법 및 2구획 축적(`GENP` → `KIN`)의 근거.

---

## 7. 수유 전략 · 증량 속도 (Feeding strategy and advance rate)
*[모델 항: `feedrate`, `Ftroph`, `Sin`, `demand` — 양쪽 부호를 가진 입력]*

48. Dorling J, Abbott J, Berrington J, et al. **Controlled trial of two
    incremental milk-feeding rates in preterm infants (SIFT).**
    *N Engl J Med.* 2019;381(15):1434-1443. (PMID 31597020) —
    <https://pubmed.ncbi.nlm.nih.gov/31597020/>
    → 증량 속도(18 vs 30 mL/kg/d)가 NEC에 차이를 만들지 않았다는 대규모 결과.
    모델은 이것을 성숙도가 높은 아이에서 **기울기가 평평해지는 구간**으로 재현한다.
49. Leaf A, Dorling J, Kempley S, et al. **Early or delayed enteral feeding for
    preterm growth-restricted infants: a randomized trial (ADEPT).**
    *Pediatrics.* 2012. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Leaf+early+or+delayed+enteral+feeding+for+preterm+growth-restricted+infants+randomized+trial>
50. Oddie SJ, Young L, McGuire W. **Slow advancement of enteral feed volumes to
    prevent necrotising enterocolitis in very low birth weight infants.**
    *Cochrane Database Syst Rev.* 2021. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Oddie+slow+advancement+of+enteral+feed+volumes+to+prevent+necrotising+enterocolitis+very+low+birth+weight>
51. Morgan J, Young L, McGuire W. **Delayed introduction of progressive enteral
    feeds to prevent necrotising enterocolitis in very low birth weight
    infants.** *Cochrane Database Syst Rev.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Morgan+delayed+introduction+of+progressive+enteral+feeds+to+prevent+necrotising+enterocolitis>
52. Berseth CL, Bisquera JA, Paje VU. **Prolonging small feeding volumes early
    in life decreases the incidence of necrotizing enterocolitis in very low
    birth weight infants.** *Pediatrics.* 2003. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Berseth+prolonging+small+feeding+volumes+early+in+life+decreases+incidence+of+necrotizing+enterocolitis>
    → 모델이 실제로 재현하는 방향(느린 증량 = 낮은 NEC)의 근거.
53. Ellsbury DL, Ursprung R. **Comprehensive Oxygen Management for the
    Prevention of Retinopathy of Prematurity** *(및 NICU 표준화 관련)*. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ellsbury+Ursprung+comprehensive+oxygen+management+prevention+retinopathy+of+prematurity>

---

## 8. 장혈류 · 허혈 · PAF (Splanchnic perfusion, ischaemia, and PAF)
*[모델 항: `PERF`, `demand`, `SDRcrit`, `PAF`, `PAFAH = dPAF0 + dPAFm·MAT`]*

54. Caplan MS, Sun XM, Hsueh W, Hageman JR. **Role of platelet activating factor
    and tumor necrosis factor-alpha in neonatal necrotizing enterocolitis.**
    *J Pediatr.* 1990. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Caplan+role+of+platelet+activating+factor+and+tumor+necrosis+factor+alpha+in+neonatal+necrotizing+enterocolitis>
    → PAF 축의 원전. 모델에서 PAF는 **혈류를 떨어뜨리는 유일한 내인성 항**이다.
55. Caplan MS, Lickerman M, Adler L, et al. **The role of recombinant
    platelet-activating factor acetylhydrolase in a neonatal rat model of
    necrotizing enterocolitis.** *Pediatr Res.* 1997. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Caplan+role+of+recombinant+platelet-activating+factor+acetylhydrolase+neonatal+rat+model+necrotizing+enterocolitis>
    → 미숙아의 PAF-AH 결핍과 그 보충 실험 — `dPAF0`(바닥값)과 `dPAFm·MAT`의 근거.
56. Yazji I, Sodhi CP, Robinson EK, et al. **Endothelial TLR4 activation impairs
    intestinal microcirculatory perfusion in necrotizing enterocolitis via
    eNOS-NO-nitrite signaling.** *Proc Natl Acad Sci U S A.* 2013. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Yazji+endothelial+TLR4+activation+impairs+intestinal+microcirculatory+perfusion+necrotizing+enterocolitis+eNOS>
    → TLR4 → 혈류 저하 경로. 모델에서 염증과 허혈이 **하나의 `INJ`로 합산되는**
    구조적 근거.
57. Nowicki PT. **Ischemia and necrotizing enterocolitis: where, when, and
    how.** *Semin Pediatr Surg.* 2005. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Nowicki+ischemia+and+necrotizing+enterocolitis+where+when+and+how>
    → 수유 후 산소 요구 증가와 자동조절 미숙 — 모델의 `demand = 1 + aDEM·FEED/200`
    및 `SDRcrit` 예비능 개념의 근거.
58. Reber KM, Nankervis CA, Nowicki PT. **Newborn intestinal circulation:
    physiology and pathophysiology.** *Clin Perinatol.* 2002. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Reber+Nankervis+Nowicki+newborn+intestinal+circulation+physiology+and+pathophysiology>
59. Ohlsson A, Walia R, Shah SS. **Ibuprofen for the treatment of patent
    ductus arteriosus in preterm or low birth weight (or both) infants.**
    *Cochrane Database Syst Rev.* 2020. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ohlsson+ibuprofen+for+the+treatment+of+patent+ductus+arteriosus+in+preterm+or+low+birth+weight+infants+Cochrane>
    → 인도메타신 대비 이부프로펜의 장혈류 영향 차이 — 모델의 `EC50IND` vs `EC50IBU`.
60. Coombs RC, Morgan ME, Durbin GM, et al. **Gut blood flow velocities in the
    newborn: effects of patent ductus arteriosus and parenteral indomethacin.**
    *Arch Dis Child.* 1990. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Coombs+gut+blood+flow+velocities+in+the+newborn+effects+of+patent+ductus+arteriosus+and+parenteral+indomethacin>

---

## 9. 빈혈 · 수혈 관련 NEC (Anaemia, transfusion, and TANEC)
*[모델 항: 침습 삽화(`INS*DEP`, `INS*NFK`) — 분수령을 넘기는 방아쇠]*

61. Patel RM, Knezevic A, Shenvi N, et al. **Association of red blood cell
    transfusion, anemia, and necrotizing enterocolitis in very low-birth-weight
    infants.** *JAMA.* 2016;315(9):889-897. (PMID 26954412) —
    <https://pubmed.ncbi.nlm.nih.gov/26954412/>
    → 빈혈 자체가 위험이고 수혈은 그 위험을 즉각 바꾸지 않는다 — 모델에서
    수혈 전 빈혈 구간과 수혈 구간을 **서로 다른 두 boxcar**로 넣은 근거.
62. Mohamed A, Shah PS. **Transfusion associated necrotizing enterocolitis: a
    meta-analysis of observational data.** *Pediatrics.* 2012. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Mohamed+Shah+transfusion+associated+necrotizing+enterocolitis+meta-analysis+of+observational+data>
63. Marin T, Josephson CD, Kosmetatos N, et al. **Feeding preterm infants
    during red blood cell transfusion is associated with a decline in
    postprandial mesenteric oxygenation.** *J Pediatr.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Marin+feeding+preterm+infants+during+red+blood+cell+transfusion+decline+in+postprandial+mesenteric+oxygenation>
    → 수혈 중 수유가 장 산소 수급을 악화시킨다 — 모델의 `demand` × 침습 상호작용.
64. Kirpalani H, Bell EF, Hintz SR, et al. **Higher or lower hemoglobin
    transfusion thresholds for preterm infants (TOP trial).**
    *N Engl J Med.* 2020. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Kirpalani+higher+or+lower+hemoglobin+transfusion+thresholds+for+preterm+infants>

---

## 10. 장벽 기능 · 밀착연접 · 점액 (Barrier, tight junctions, mucus)
*[모델 항: `TJ`, `MUC`, `BI`, `Pb = Pmin + (Pmax−Pmin)(1−BI)^hP`]*

65. Anderson JM, Van Itallie CM. **Physiology and function of the tight
    junction.** *Cold Spring Harb Perspect Biol.* 2009. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Anderson+Van+Itallie+physiology+and+function+of+the+tight+junction>
66. Bein A, Eventov-Friedman S, Arbell D, Schwartz B. **Intestinal tight
    junctions are severely altered in NEC preterm neonates.**
    *Pediatr Neonatol.* 2018. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Bein+intestinal+tight+junctions+are+severely+altered+in+NEC+preterm+neonates>
67. McElroy SJ, Prince LS, Weitkamp JH, et al. **Tumor necrosis factor
    receptor 1-dependent depletion of mucus in immature small intestine: a
    potential role in neonatal necrotizing enterocolitis.**
    *Am J Physiol Gastrointest Liver Physiol.* 2011. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=McElroy+tumor+necrosis+factor+receptor+1+dependent+depletion+of+mucus+in+immature+small+intestine+necrotizing+enterocolitis>
68. Halpern MD, Denning PW. **The role of intestinal epithelial barrier
    function in the development of NEC.** *Tissue Barriers.* 2015. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Halpern+Denning+role+of+intestinal+epithelial+barrier+function+in+the+development+of+NEC>
69. Whitehouse JS, Riggle KM, Purpi DP, et al. **The protective role of
    intestinal alkaline phosphatase in necrotizing enterocolitis.**
    *J Surg Res.* 2010. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Whitehouse+protective+role+of+intestinal+alkaline+phosphatase+in+necrotizing+enterocolitis>
    → 모델의 IAP 노드(LPS 탈인산화로 신호 문턱을 올리는 개입).
70. Underwood MA. **Paneth cells and necrotizing enterocolitis.**
    *Gut Microbes.* 2012. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Underwood+Paneth+cells+and+necrotizing+enterocolitis>
    → `DEF`(디펜신) 항의 근거.

---

## 11. 사이토카인 · 세포 면역 (Cytokines and cellular immunity)
*[모델 항: `TNFA`, `IL1B`, `IL8`, `IL10`, `NEU`(포화 되먹임), `NOX`]*

71. Egan CE, Sodhi CP, Good M, et al. **Toll-like receptor 4-mediated
    lymphocyte influx induces neonatal necrotizing enterocolitis.**
    *J Clin Invest.* 2016. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Egan+Toll-like+receptor+4+mediated+lymphocyte+influx+induces+neonatal+necrotizing+enterocolitis>
72. Weitkamp JH, Koyama T, Rock MT, et al. **Necrotising enterocolitis is
    characterised by disrupted immune regulation and diminished mucosal
    regulatory (FOXP3)/effector (CD4, CD8) T cell ratios.** *Gut.* 2013. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Weitkamp+necrotising+enterocolitis+disrupted+immune+regulation+diminished+mucosal+regulatory+FOXP3+effector+T+cell+ratios>
73. MohanKumar K, Namachivayam K, Song T, et al. **A murine neonatal model of
    necrotizing enterocolitis caused by anemia and red blood cell
    transfusions.** *Nat Commun.* 2019. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=MohanKumar+murine+neonatal+model+of+necrotizing+enterocolitis+caused+by+anemia+and+red+blood+cell+transfusions>
74. Emami CN, Chokshi N, Wang J, et al. **Role of interleukin-10 in the
    pathogenesis of necrotizing enterocolitis.** *Am J Surg.* 2012. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Emami+role+of+interleukin-10+in+the+pathogenesis+of+necrotizing+enterocolitis>
    → `kIL10·MAT` — 미숙할수록 IL-10 대항조절이 약하다는 항의 근거.
75. Chan KL, Ho JC, Chan KW, Tam PK. **A study of gut immunity to enteral
    endotoxin in rats of different ages.** *Pediatr Surg Int.* 2002. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Chan+study+of+gut+immunity+to+enteral+endotoxin+in+rats+of+different+ages>
76. Ford H, Watkins S, Reblock K, Rowe M. **The role of inflammatory cytokines
    and nitric oxide in the pathogenesis of necrotizing enterocolitis.**
    *J Pediatr Surg.* 1997. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ford+role+of+inflammatory+cytokines+and+nitric+oxide+in+the+pathogenesis+of+necrotizing+enterocolitis>
    → `wNO·NOX` 항(질산화 손상)의 근거.

---

## 12. 진단 표지자 (Diagnostic biomarkers)
*[모델 항: `CRP`, `PLT`, `LAC`, `PNEU`, I-FABP · calprotectin 판독 노드]*

77. Thuijls G, Derikx JPM, van Wijck K, et al. **Non-invasive markers for early
    diagnosis and determination of the severity of necrotizing enterocolitis.**
    *Ann Surg.* 2010. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Thuijls+non-invasive+markers+for+early+diagnosis+and+determination+of+the+severity+of+necrotizing+enterocolitis>
78. Ng PC, Ma TP, Lam HS. **The use of laboratory biomarkers for surveillance,
    diagnosis and prediction of clinical outcomes in neonatal sepsis and
    necrotising enterocolitis.** *Arch Dis Child Fetal Neonatal Ed.* 2015. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ng+use+of+laboratory+biomarkers+for+surveillance+diagnosis+prediction+clinical+outcomes+neonatal+sepsis+necrotising+enterocolitis>
79. Ragazzi S, Pierro A, Peters M, et al. **Early full blood count and severity
    of disease in neonates with necrotizing enterocolitis.**
    *Pediatr Surg Int.* 2003. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ragazzi+early+full+blood+count+and+severity+of+disease+in+neonates+with+necrotizing+enterocolitis>
    → 혈소판 감소와 중증도의 관계 — 모델의 Bell III 판정에 `PLT < 100`이 들어간 근거.
80. Christensen RD, Yoder BA, Baer VL, et al. **Early-onset neutropenia in
    small-for-gestational-age infants** *(및 NEC 혈액학적 소견 관련 연구들)*. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Christensen+thrombocytopenia+necrotizing+enterocolitis+neonate>
81. Sylvester KG, Ling XB, Liu GY, et al. **Urine protein biomarkers for the
    diagnosis and prognosis of necrotizing enterocolitis in infants.**
    *J Pediatr.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Sylvester+urine+protein+biomarkers+for+the+diagnosis+and+prognosis+of+necrotizing+enterocolitis+in+infants>

---

## 13. 수술 · 결과 · 후유증 (Surgery, outcomes, sequelae)
*[모델 항: `NECrev`, `NECa_crit = INJth/wNEC`, 협착 · 단장증후군 · NDI]*

82. Moss RL, Dimmitt RA, Barnhart DC, et al. **Laparotomy versus peritoneal
    drainage for necrotizing enterocolitis and perforation.**
    *N Engl J Med.* 2006;354(21):2225-2234. (PMID 16707747) —
    <https://pubmed.ncbi.nlm.nih.gov/16707747/>
83. Blakely ML, Tyson JE, Lally KP, et al. **Initial laparotomy versus
    peritoneal drainage in extremely low birthweight infants with surgical
    necrotizing enterocolitis or isolated intestinal perforation (NEST).**
    *Ann Surg.* 2021. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Blakely+initial+laparotomy+versus+peritoneal+drainage+extremely+low+birthweight+surgical+necrotizing+enterocolitis+isolated+intestinal+perforation>
84. Rees CM, Eaton S, Kiely EM, et al. **Peritoneal drainage or laparotomy for
    neonatal bowel perforation? A randomized controlled trial.**
    *Ann Surg.* 2008. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Rees+peritoneal+drainage+or+laparotomy+for+neonatal+bowel+perforation+randomized+controlled+trial>
85. Hull MA, Fisher JG, Gutierrez IM, et al. **Mortality and management of
    surgical necrotizing enterocolitis in very low birth weight neonates: a
    prospective cohort study.** *J Am Coll Surg.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Hull+mortality+and+management+of+surgical+necrotizing+enterocolitis+very+low+birth+weight+neonates+prospective+cohort>
86. Hintz SR, Kendrick DE, Stoll BJ, et al. **Neurodevelopmental and growth
    outcomes of extremely low birth weight infants after necrotizing
    enterocolitis.** *Pediatrics.* 2005. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Hintz+neurodevelopmental+and+growth+outcomes+of+extremely+low+birth+weight+infants+after+necrotizing+enterocolitis>
    → 모델의 `NIN` → NDI 판독 경로의 근거.
87. Fredriksson F, Engstrand Lilja H. **Survival rates for surgically treated
    necrotising enterocolitis have improved over the last four decades.**
    *Acta Paediatr.* 2019. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Fredriksson+Engstrand+Lilja+survival+rates+for+surgically+treated+necrotising+enterocolitis+improved+over+the+last+four+decades>
88. Niemarkt HJ, De Meij TG, van Ganzewinkel CJ, et al. **Necrotizing
    enterocolitis, gut microbiota, and brain development: role of the
    brain-gut axis.** *Neonatology.* 2019. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Niemarkt+necrotizing+enterocolitis+gut+microbiota+and+brain+development+role+of+the+brain-gut+axis>
89. Robinson JR, Rellinger EJ, Hatch LD, et al. **Surgical necrotizing
    enterocolitis.** *Semin Perinatol.* 2017. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Robinson+surgical+necrotizing+enterocolitis+Semin+Perinatol>
90. Zani A, Pierro A. **Necrotizing enterocolitis: controversies and
    challenges.** *F1000Res.* 2015. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Zani+Pierro+necrotizing+enterocolitis+controversies+and+challenges>

---

## 14. 신생아 약동학 (Neonatal pharmacokinetics)
*[모델 항: `CL·(PMA/40)^1.3`, 각 약물의 V·CL, `GENP` 축적 → `KIN`]*

91. Anderson BJ, Holford NHG. **Mechanism-based concepts of size and maturity
    in pharmacokinetics.** *Annu Rev Pharmacol Toxicol.* 2008. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Anderson+Holford+mechanism-based+concepts+of+size+and+maturity+in+pharmacokinetics>
    → 모델이 청소율을 PMA로 스케일한 방식(성숙 함수)의 근거.
92. Tremoulet A, Le J, Poindexter B, et al. **Characterization of ampicillin
    pharmacokinetics in neonates using a population approach.**
    *Antimicrob Agents Chemother.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Tremoulet+characterization+of+ampicillin+pharmacokinetics+in+neonates+using+a+population+approach>
93. Fuchs A, Guidi M, Giannoni E, et al. **Population pharmacokinetic study of
    gentamicin in a large cohort of premature and term neonates.**
    *Br J Clin Pharmacol.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Fuchs+population+pharmacokinetic+study+of+gentamicin+in+a+large+cohort+of+premature+and+term+neonates>
94. Suyagh M, Collier PS, Millership JS, et al. **Metronidazole population
    pharmacokinetics in preterm neonates using dried blood-spot sampling.**
    *Pediatrics.* 2011. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Suyagh+metronidazole+population+pharmacokinetics+in+preterm+neonates+dried+blood-spot+sampling>
95. Smyth JM, Collier PS, Darwish M, et al. **Intravenous indometacin in
    preterm infants with symptomatic patent ductus arteriosus: a population
    pharmacokinetic study.** *Br J Clin Pharmacol.* 2004. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Smyth+intravenous+indometacin+in+preterm+infants+with+symptomatic+patent+ductus+arteriosus+population+pharmacokinetic>
96. Hirt D, Van Overmeire B, Treluyer JM, et al. **An optimized ibuprofen
    dosing scheme for preterm neonates with patent ductus arteriosus, based on
    a population pharmacokinetic and pharmacodynamic study.**
    *Br J Clin Pharmacol.* 2008. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Hirt+optimized+ibuprofen+dosing+scheme+for+preterm+neonates+with+patent+ductus+arteriosus+population+pharmacokinetic>

---

## 15. QSP 방법론 · 이중안정성 · 임계 전이 (QSP method, bistability, critical transitions)
*[모델 항: `g(E)`의 근 개수, 안장-결절 분기 `B_lo`/`B_hi`, 분수령 `E*`의 로그 전이시간]*

97. Strogatz SH. **Nonlinear Dynamics and Chaos.** Westview Press. —
    안장-결절 분기와 히스테리시스의 표준 교과서적 서술.
    (도서 — PubMed 항목 없음)
98. Scheffer M, Bascompte J, Brock WA, et al. **Early-warning signals for
    critical transitions.** *Nature.* 2009. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Scheffer+early-warning+signals+for+critical+transitions>
    → 임계 전이 직전에 회복 속도가 느려진다(critical slowing down)는 결과.
    모델의 `E*` 근처 고유값 λ ≈ 0.35/d 및 로그 전이시간과 직접 대응한다.
99. Angeli D, Ferrell JE, Sontag ED. **Detection of multistability,
    bifurcations, and hysteresis in a large class of biological positive-feedback
    systems.** *Proc Natl Acad Sci U S A.* 2004. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Angeli+Ferrell+Sontag+detection+of+multistability+bifurcations+and+hysteresis+biological+positive-feedback+systems>
100. Baker RE, Peña JM, Jayamohan J, Jérusalem A. **Mechanistic models versus
     machine learning, a fight worth fighting for the biological community?**
     *Biol Lett.* 2018. —
     <https://pubmed.ncbi.nlm.nih.gov/?term=Baker+mechanistic+models+versus+machine+learning+a+fight+worth+fighting+for+the+biological+community>
101. Elmokadem A, Riggs M, Baron K. **Quantitative systems pharmacology and
     physiologically-based pharmacokinetic modeling with mrgsolve: a hands-on
     tutorial.** *CPT Pharmacometrics Syst Pharmacol.* 2019. —
     <https://pubmed.ncbi.nlm.nih.gov/?term=Elmokadem+quantitative+systems+pharmacology+and+physiologically-based+pharmacokinetic+modeling+with+mrgsolve+hands-on+tutorial>
     → 이 저장소의 mrgsolve 구현 관례.

---

## 16. 종합 리뷰 (Comprehensive reviews used for cross-checking)

102. Frost BL, Modi BP, Jaksic T, Caplan MS. **New medical and surgical insights
     into neonatal necrotizing enterocolitis.** *JAMA Pediatr.* 2017. —
     <https://pubmed.ncbi.nlm.nih.gov/?term=Frost+new+medical+and+surgical+insights+into+neonatal+necrotizing+enterocolitis>
103. Neu J. **Necrotizing enterocolitis: the future.** *Neonatology.* 2020. —
     <https://pubmed.ncbi.nlm.nih.gov/?term=Neu+necrotizing+enterocolitis+the+future+Neonatology>
104. Bazacliu C, Neu J. **Necrotizing enterocolitis: long term complications.**
     *Curr Pediatr Rev.* 2019. —
     <https://pubmed.ncbi.nlm.nih.gov/?term=Bazacliu+Neu+necrotizing+enterocolitis+long+term+complications>
105. Denning NL, Prince JM. **Neonatal intestinal dysbiosis in necrotizing
     enterocolitis.** *Mol Med.* 2018. —
     <https://pubmed.ncbi.nlm.nih.gov/?term=Denning+Prince+neonatal+intestinal+dysbiosis+in+necrotizing+enterocolitis>
106. Tanner SM, Berryhill TF, Ellenburg JL, et al. **Pathogenesis of
     necrotizing enterocolitis: modeling the innate immune response.**
     *Am J Pathol.* 2015. —
     <https://pubmed.ncbi.nlm.nih.gov/?term=Tanner+pathogenesis+of+necrotizing+enterocolitis+modeling+the+innate+immune+response>

---

## 이 모델이 만드는, 문헌으로 아직 검증되지 않은 예측
## Predictions this model makes that the literature has NOT yet tested

정직하게 적자면, 아래는 **위 문헌에서 나온 것이 아니라 모델의 구조에서 나온 것**
입니다. 반증 가능하도록 명시합니다.

1. **모유의 보호효과는 HMO를 대사할 수 있는 균주(`binfant`)의 보유 여부에 따라
   갈린다.** HMO 농도만 재는 연구는 효과를 희석해서 볼 것이다 (참고 #29, #30, #34
   가 이 방향을 시사하지만 상호작용을 직접 검정한 연구는 아직 없다).
2. **경험적 항생제의 해악은 모유 수유아에서만 관측 가능해야 한다.** 잃을 상재균이
   없는 완전 분유 수유아에서는 항생제가 오히려 보호적으로 보일 수 있다 (참고 #43
   의 코호트를 수유 형태로 층화하면 검정 가능).
3. **회복 불가점은 병변 크기의 함수이며 시간의 함수가 아니다.** 모델에서
   `NECa_crit = INJth/wNEC` 이고, 진단-치료 지연(`dx_lag`)이 이 문턱을 넘기느냐가
   내과적 NEC와 수술적 NEC를 가른다.
4. **Bell II에서 울리는 표지자는 원리적으로 시간을 거의 벌어주지 못한다.** 분수령
   근처의 전이시간이 로그 의존이므로, 앞서 움직이는 것은 상태변수(`E`, `Jtr`)이고
   임상 표지자는 이미 스위치가 넘어간 뒤에 움직인다.

---

**면책 (Disclaimer).** 본 문헌 목록과 모델은 교육 및 연구 목적입니다. 파라미터는
공개 문헌의 **집계 지표**에 맞춰 손으로 보정되었으며 개별 환자 데이터에 적합된
것이 아닙니다. 임상 의사결정·처방·규제 제출에 사용해서는 안 됩니다.
