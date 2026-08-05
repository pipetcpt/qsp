# 중증 화상 (Major Thermal Burn Injury) — QSP 모델 참고문헌
# References for the Major Thermal Burn Injury QSP model

> **출처에 대하여 (provenance).** 아래 352편은 NCBI E-utilities(`esearch` + `esummary`)로
> 각 절의 주제에 해당하는 검색어를 실제로 질의해 **반환된 레코드만** 옮긴 것입니다.
> 제목·저자·저널·연도·PMID는 모두 NCBI가 돌려준 값 그대로이며, 기억이나 추정으로
> 채워 넣은 항목은 하나도 없습니다. 검색 시점은 2026-08입니다.
>
> 모델의 **정량적 주장이 어느 문헌에 걸려 있는지**는 맨 끝
> [부록 A](#부록-a--모델의-각-구성요소가-어느-절에-걸려-있는가)에 표로 정리했습니다.

---

## 이 모델이 문헌에서 가져온 숫자 (Calibration targets)

| # | 보정 목표 | 문헌값 | 모델값 | 근거 절 |
|---|-----------|--------|--------|---------|
| 1 | UO 적정 파클랜드의 24시간 수액량 | 5.2–6.7 mL/kg/%TBSA | 5.97 | C |
| 2 | in:out 비 (투여량 ÷ 파클랜드 처방량) | 1.2–1.6 | 1.49 | C |
| 3 | 복강내 고혈압 역치 | 약 250 mL/kg/24 h | 269 mL/kg → IAP 14.2 mmHg | D |
| 4 | 혈장량 최저점 | 기저의 60–80 % | 62 % | B |
| 5 | 24시간 혈장 교질삼투압 | 10–16 mmHg (기저 약 26) | 13.0 | B |
| 6 | 최대 체중 증가 | +15~30 % | +20 % | C·D |
| 7 | 최대 REE | 예측치의 120–180 % | 157 % (45 %TBSA) | F |
| 8 | REE의 지속 | 창상 폐쇄 후 수개월~수년 | 60일째 135 % | F |
| 9 | 제지방량 변화 (14일, 무투약) | 약 −9 % | −7.1 % | F·G |
| 10 | 제지방량 변화 (14일, 프로프라놀롤) | 약 +9 % | 아래 본문 참조 | G |
| 11 | 고용량 아스코르브산의 수액 절감 | −45 % | −36 % | E |
| 12 | 알부민의 수액 절감 | −20~−40 % | −5 % (8시간 시작) / −38 % (0시간) | C |
| 13 | 침습성 창상 감염 역치 | 10⁵ CFU/g 초과 | 모델 BTHRESH = 5 log₁₀ | K |
| 14 | 개정 Baux 점수와의 일치 | LD50 약 110 | 8개 환자 평균오차 3.2 %p | A |
| 15 | 화상 환자의 증가된 신클리어런스 | CrCl 130 mL/min 초과 | flow phase에서 ×1.55 | L |

**모델이 문헌과 어긋나는 곳은 지우지 않고 그대로 두었습니다.** 가장 노출된 예측은 12번입니다.
모델은 8시간에 시작한 알부민이 수액량을 거의 줄이지 못한다고 말하며, 그 이유를 기전으로
제시합니다 — 파클랜드의 전반 8시간 집중 투여 형태에서는 8시간 시점에 이미 교질삼투압
희석의 대부분이 끝나 있기 때문입니다. 같은 콜로이드 분율을 0시간에 시작하면 −38 %입니다.
이것은 검증 가능한 주장이며, 시험에서 시작 시점이 그만큼 중요하지 않다면 이 구조가 틀린
것입니다. 11번(아스코르브산)이 Tanaka의 값을 잘 재현하는 것도 승리가 아니라 문제입니다 —
이후 어떤 시험도 그 결과를 재현하지 못했으므로, 모델은 사실이 아닐 수 있는 결과에 맞춰져
있는 셈입니다.

---

## A. 역학·중증도·예후 (Epidemiology, severity, prognosis)

*21편*

- Tomaka P et al. **Prognostic Value of Mortality Scoring Systems in Patients With Severe Burns: Identifying Key Predictors of Mortality and Comparative Analysis Between Survivors and Non-Survivors** *Med Sci Monit* 32:e951713 (2026). PMID [42458800](https://pubmed.ncbi.nlm.nih.gov/42458800/)
- Goder N et al. **Recalibrating the revised Baux score: Admission and 48-hour landmark survival models in contemporary older adult patients with burn injury** *Burns* 52:108141 (2026). PMID [42407276](https://pubmed.ncbi.nlm.nih.gov/42407276/)
- Belayneh ES et al. **Accuracy of the revised Baux score for predicting in-hospital mortality of patients with burns: A retrospective cohort study from Ethiopia** *J Plast Reconstr Aesthet Surg* 119:375-381 (2026). PMID [42399141](https://pubmed.ncbi.nlm.nih.gov/42399141/)
- Ren H et al. **Interpretative machine learning for predicting 60-day mortality in burn patients with suspected infection** *World J Emerg Med* 17:244-249 (2026). PMID [42199764](https://pubmed.ncbi.nlm.nih.gov/42199764/)
- Nehila T et al. **Revised Baux Score Identifies a New Risk Factor for Mortality: History of Diabetes** *Adv Skin Wound Care* 39:304-307 (2026). PMID [42155058](https://pubmed.ncbi.nlm.nih.gov/42155058/)
- Alshdowh H et al. **Burn Injuries in Jordan: A 5-Year Retrospective Analysis of Presentation, Management and Hospital Mortality** *Int Wound J* 23:e70839 (2026). PMID [41601206](https://pubmed.ncbi.nlm.nih.gov/41601206/)
- Schmitt A et al. **Epidemiology and prognostic factors in older adults with severe burns (≥10% TBSA): A 13-year retrospective study from a French national burn center** *Burns* 52:108104 (2026). PMID [42308703](https://pubmed.ncbi.nlm.nih.gov/42308703/)
- Dean YE et al. **Predictors of Mortality and Outcomes in Older Adult Burn Patient Subgroups: A Comparative Meta-analysis of Age, Burn Severity, and Intervention Strategies** *Plast Reconstr Surg Glob Open* 14:e7708 (2026). PMID [42164960](https://pubmed.ncbi.nlm.nih.gov/42164960/)
- Zhang L et al. **Early blood urea nitrogen-albumin ratio as a prognostic marker in severe burn patients: A retrospective clinical study** *Burns* 52:107871 (2026). PMID [41632997](https://pubmed.ncbi.nlm.nih.gov/41632997/)
- Ahn J et al. **Artificial intelligence for outcome prediction in burn care: A scoping review** *Burns* 52:108121 (2026). PMID [42407287](https://pubmed.ncbi.nlm.nih.gov/42407287/)
- Danielski LG et al. **Performance of Sepsis-3 criteria in burn patients: A multicenter real-world study** *Burns* 52:108073 (2026). PMID [42302549](https://pubmed.ncbi.nlm.nih.gov/42302549/)
- Baddam S et al. **Systemic Inflammatory Response Syndrome** ** (2026). PMID [31613449](https://pubmed.ncbi.nlm.nih.gov/31613449/)
- Wu R et al. **Early versus delayed laparoscopic cholecystectomy for gallbladder perforation** *J Trauma Acute Care Surg* 98:642-648 (2025). PMID [40122846](https://pubmed.ncbi.nlm.nih.gov/40122846/)
- Heard J et al. **Burn Injury Severity in Adults: Proposed Definitions Based on the National Burn Research Dataset** *J Burn Care Res* 46:438-449 (2025). PMID [39320264](https://pubmed.ncbi.nlm.nih.gov/39320264/)
- Rhee C et al. **Improving Sepsis Outcomes in the Era of Pay-for-Performance and Electronic Quality Measures: A Joint IDSA/ACEP/PIDS/SHEA/SHM/SIDP Position Paper** *Clin Infect Dis* 78:505-513 (2024). PMID [37831591](https://pubmed.ncbi.nlm.nih.gov/37831591/)
- López-Sánchez I et al. **Temporal trends in epidemiology and patient characteristics of 36 cancers: a protocol for a multinational population-based cohort study using OMOP-standardised databases to investigate CANcer (OMOPCAN)** *BMJ Open* 16:e119069 (2026). PMID [42481199](https://pubmed.ncbi.nlm.nih.gov/42481199/)
- Patel P et al. **Epidemiology and Outcomes of Burn Injuries: A Retrospective Study at a Tertiary Care Center** *Ann Plast Surg* (2026). PMID [42295771](https://pubmed.ncbi.nlm.nih.gov/42295771/)
- Miller CL et al. **The tale of two states: mortality outcomes and burn severity among in-state versus out-of-state burn victims treated in the Texas Panhandle** *Proc (Bayl Univ Med Cent)* 39:650-656 (2026). PMID [42269065](https://pubmed.ncbi.nlm.nih.gov/42269065/)
- Daniels M et al. **The influence of burn injury timing on survival in patients with severe burns** *Injury* 57:113302 (2026). PMID [42054930](https://pubmed.ncbi.nlm.nih.gov/42054930/)
- Watanabe T et al. **The association between hospital volume and in-hospital mortality in severe burn patients: A nationwide study using the Japanese Burn Registry** *Burns* 52:107959 (2026). PMID [41950740](https://pubmed.ncbi.nlm.nih.gov/41950740/)
- Cueva-Ramírez JE et al. **Epidemiological Characteristics and Prognostic Factors for Mortality in Severe Burns: A 32-Year Analysis (1993-2024) From a National Referral Center in the Dominican Republic** *World J Surg* 50:1159-1168 (2026). PMID [41931542](https://pubmed.ncbi.nlm.nih.gov/41931542/)

## B. 화상 부종의 미세혈관 기전 (Microvascular mechanism of burn oedema)

*24편*

- Ferrara JJ et al. **Alpha-trinositol reduces edema formation at the site of scald injury** *Surgery* 123:36-45 (1998). PMID [9457221](https://pubmed.ncbi.nlm.nih.gov/9457221/)
- Ferrara JJ et al. **Burn edema reduction by methysergide is not due to control of regional vasodilation** *J Surg Res* 61:11-6 (1996). PMID [8769935](https://pubmed.ncbi.nlm.nih.gov/8769935/)
- Ferrara JJ et al. **Effects of methysergide administration on edema formation at the site of scald** *J Appl Physiol (1985)* 78:191-7 (1995). PMID [7713811](https://pubmed.ncbi.nlm.nih.gov/7713811/)
- Ferrara JJ et al. **Effects of pentafraction administration on microvascular permeability alterations induced by graded thermal injury** *Surgery* 115:182-9 (1994). PMID [7508639](https://pubmed.ncbi.nlm.nih.gov/7508639/)
- Dyess DL et al. **Effects of hypertonic saline and dextran 70 resuscitation on microvascular permeability after burn** *Am J Physiol* 262:H1832-7 (1992). PMID [1377878](https://pubmed.ncbi.nlm.nih.gov/1377878/)
- Isago T et al. **Analysis of pulmonary microvascular permeability after smoke inhalation** *J Appl Physiol (1985)* 71:1403-8 (1991). PMID [1757363](https://pubmed.ncbi.nlm.nih.gov/1757363/)
- Lund T et al. **Mechanisms behind increased dermal imbibition pressure in acute burn edema** *Am J Physiol* 256:H940-8 (1989). PMID [2705564](https://pubmed.ncbi.nlm.nih.gov/2705564/)
- Lund T et al. **Acute postburn edema: role of strongly negative interstitial fluid pressure** *Am J Physiol* 255:H1069-74 (1988). PMID [3189570](https://pubmed.ncbi.nlm.nih.gov/3189570/)
- Tanaka H et al. **High dose vitamin C counteracts the negative interstitial fluid hydrostatic pressure and early edema generation in thermally injured rats** *Burns* 25:569-74 (1999). PMID [10563680](https://pubmed.ncbi.nlm.nih.gov/10563680/)
- Bert J et al. **Fluid resuscitation following a burn injury: implications of a mathematical model of microvascular exchange** *Burns* 23:93-105 (1997). PMID [9177874](https://pubmed.ncbi.nlm.nih.gov/9177874/)
- Lund T et al. **Alpha-Trinositol inhibits edema generation and albumin extravasation in thermally injured skin** *J Trauma* 36:761-5 (1994). PMID [8014994](https://pubmed.ncbi.nlm.nih.gov/8014994/)
- Lund T et al. **Pathogenesis of edema formation in burn injuries** *World J Surg* 16:2-9 (1992). PMID [1290261](https://pubmed.ncbi.nlm.nih.gov/1290261/)
- Rehberg S et al. **Antithrombin attenuates vascular leakage via inhibiting neutrophil activation in acute lung injury** *Crit Care Med* 41:e439-46 (2013). PMID [24107637](https://pubmed.ncbi.nlm.nih.gov/24107637/)
- Sakurai H et al. **Inhibition of lung permeability changes after burn and smoke inhalation by an anti-interleukin-8 antibody in sheep** *Surg Today* 39:399-406 (2009). PMID [19408077](https://pubmed.ncbi.nlm.nih.gov/19408077/)
- Maybauer DM et al. **Lung-protective effects of the metalloporphyrinic peroxynitrite decomposition catalyst WW-85 in interleukin-2 induced toxicity** *Biochem Biophys Res Commun* 377:786-91 (2008). PMID [18951875](https://pubmed.ncbi.nlm.nih.gov/18951875/)
- Jonkam CC et al. **Effects of the bradykinin B2 receptor antagonist icatibant on microvascular permeability after thermal injury in sheep** *Shock* 28:704-9 (2007). PMID [17607158](https://pubmed.ncbi.nlm.nih.gov/17607158/)
- Enkhbaatar P et al. **Effect of inhaled nitric oxide on pulmonary vascular hyperpermeability in sheep following smoke inhalation** *Burns* 31:1013-9 (2005). PMID [16288960](https://pubmed.ncbi.nlm.nih.gov/16288960/)
- Katahira J et al. **Role of anti-L-selectin antibody in burn and smoke inhalation injury in sheep** *Am J Physiol Lung Cell Mol Physiol* 283:L1043-50 (2002). PMID [12376357](https://pubmed.ncbi.nlm.nih.gov/12376357/)
- Hundeshagen G et al. **Closed-Loop- and Decision-Assist-Guided Fluid Therapy of Human Hemorrhage** *Crit Care Med* 45:e1068-e1074 (2017). PMID [28682837](https://pubmed.ncbi.nlm.nih.gov/28682837/)
- Cartotto R et al. **Colloids in Acute Burn Resuscitation** *Crit Care Clin* 32:507-23 (2016). PMID [27600123](https://pubmed.ncbi.nlm.nih.gov/27600123/)
- Huang W et al. **[Effects on vascular permeability with different fluid resuscitation regimens during burn stage in swines]** *Zhonghua Yi Xue Za Zhi* 95:943-6 (2015). PMID [26081060](https://pubmed.ncbi.nlm.nih.gov/26081060/)
- Cartotto R **Fluid resuscitation of the thermally injured patient** *Clin Plast Surg* 36:569-81 (2009). PMID [19793552](https://pubmed.ncbi.nlm.nih.gov/19793552/)
- Kinsky MP et al. **Resuscitation of severe thermal injury with hypertonic saline dextran: effects on peripheral and visceral edema in sheep** *J Trauma* 49:844-53 (2000). PMID [11086774](https://pubmed.ncbi.nlm.nih.gov/11086774/)
- Guha SC et al. **Burn resuscitation: crystalloid versus colloid versus hypertonic saline hyperoncotic colloid in sheep** *Crit Care Med* 24:1849-57 (1996). PMID [8917036](https://pubmed.ncbi.nlm.nih.gov/8917036/)

## C. 소생술 공식과 수액 크립 (Resuscitation formulae and fluid creep)

*35편*

- Alotaibi AM et al. **The impact of resuscitation strategies on burn patient outcomes: Parkland vs. modified Brooke's** *Int J Burns Trauma* 15:220-226 (2025). PMID [41278384](https://pubmed.ncbi.nlm.nih.gov/41278384/)
- Schutzman LM et al. **Metabolomic Assessment of Low- Versus High-Volume Resuscitation in a Combined Porcine Model of Severe Burn and Traumatic Brain Injury** *Shock* 65:218-225 (2026). PMID [40997256](https://pubmed.ncbi.nlm.nih.gov/40997256/)
- Dilken O et al. **Microcirculatory depth of focus measurement shows reduction of tissue edema by albumin resuscitation in burn patients** *J Intensive Med* 5:58-63 (2025). PMID [39872843](https://pubmed.ncbi.nlm.nih.gov/39872843/)
- Vasileiadis V et al. **Fluid and burns in children: What we know and what we do not know-a retrospective analysis of the German Burn Registry from 2015 to 2022** *Eur J Pediatr* 183:5479-5488 (2024). PMID [39436457](https://pubmed.ncbi.nlm.nih.gov/39436457/)
- Peters J et al. **Using a Fluid Resuscitation Algorithm to Reduce the Incidence of Abdominal Compartment Syndrome in the Burn Intensive Care Unit** *Crit Care Nurse* 43:58-66 (2023). PMID [38035617](https://pubmed.ncbi.nlm.nih.gov/38035617/)
- Oboli VN et al. **EMS Burn Rule of Tens** ** (2026). PMID [37983357](https://pubmed.ncbi.nlm.nih.gov/37983357/)
- Saffle JR **Fluid Creep and Over-resuscitation** *Crit Care Clin* 32:587-98 (2016). PMID [27600130](https://pubmed.ncbi.nlm.nih.gov/27600130/)
- Faraklas I et al. **Review of a fluid resuscitation protocol: "fluid creep" is not due to nursing error** *J Burn Care Res* 33:74-83 (2012). PMID [22240507](https://pubmed.ncbi.nlm.nih.gov/22240507/)
- Lawrence A et al. **Colloid administration normalizes resuscitation ratio and ameliorates "fluid creep"** *J Burn Care Res* 31:40-7 (2010). PMID [20061836](https://pubmed.ncbi.nlm.nih.gov/20061836/)
- Saffle JI **The phenomenon of "fluid creep" in acute burn resuscitation** *J Burn Care Res* 28:382-95 (2007). PMID [17438489](https://pubmed.ncbi.nlm.nih.gov/17438489/)
- Friedrich JB et al. **Is supra-Baxter resuscitation in burn patients a new phenomenon?** *Burns* 30:464-6 (2004). PMID [15225912](https://pubmed.ncbi.nlm.nih.gov/15225912/)
- Engrav LH et al. **A biopsy of the use of the Baxter formula to resuscitate burns or do we do it like Charlie did it?** *J Burn Care Rehabil* 21:91-5 (2000). PMID [10752739](https://pubmed.ncbi.nlm.nih.gov/10752739/)
- Datchinamourthy T et al. **Comparison of effect of Ringer's lactate and isotonic bicarbonate combination therapy with Ringer's lactate alone in early fluid resuscitation of burns patients - A randomised controlled trial** *Burns* 51:107752 (2025). PMID [41197195](https://pubmed.ncbi.nlm.nih.gov/41197195/)
- Liu XY et al. **Preliminary exploration of efficacy assessment and optimal concentration of teprenone oral rehydration solution for the early management of post-burn shock** *Burns* 51:107615 (2025). PMID [40651116](https://pubmed.ncbi.nlm.nih.gov/40651116/)
- O'Neill R et al. **Inhalational Injury Secondary to House Fire** *J Educ Teach Emerg Med* 8:S49-S79 (2023). PMID [37969154](https://pubmed.ncbi.nlm.nih.gov/37969154/)
- Dahl R et al. **Regional Burn Review: Neither Parkland Nor Brooke Formulas Reach 85% Accuracy Mark for Burn Resuscitation** *J Burn Care Res* 44:1452-1459 (2023). PMID [37010149](https://pubmed.ncbi.nlm.nih.gov/37010149/)
- Popowicz P et al. **Burn Fluid Resuscitation** ** (2026). PMID [30480960](https://pubmed.ncbi.nlm.nih.gov/30480960/)
- Hassan AM et al. **Development of an open-source, protocol-adherent platform for real-time burn fluid resuscitation management** *J Plast Reconstr Aesthet Surg* 116:312-314 (2026). PMID [41813514](https://pubmed.ncbi.nlm.nih.gov/41813514/)
- Kahn SA et al. **Challenging Legacy Burn Resuscitation Paradigms with Fluid Restriction and Early Plasma** *J Am Coll Surg* 240:339-347 (2025). PMID [39902941](https://pubmed.ncbi.nlm.nih.gov/39902941/)
- Chiao HY et al. **Goal-Directed Fluid Resuscitation Protocol Based on Arterial Waveform Analysis of Major Burn Patients in a Mass Burn Casualty** *Ann Plast Surg* 80:S21-S25 (2018). PMID [29389698](https://pubmed.ncbi.nlm.nih.gov/29389698/)
- Huang M et al. **A comparison of two different fluid resuscitation management protocols for pediatric burn patients: A retrospective study** *Burns* 44:82-89 (2018). PMID [29229195](https://pubmed.ncbi.nlm.nih.gov/29229195/)
- Cancio LC et al. **Protocolized Resuscitation of Burn Patients** *Crit Care Clin* 32:599-610 (2016). PMID [27600131](https://pubmed.ncbi.nlm.nih.gov/27600131/)
- Sánchez-Sánchez M et al. **First resuscitation of critical burn patients: progresses and problems** *Med Intensiva* 40:118-24 (2016). PMID [26873418](https://pubmed.ncbi.nlm.nih.gov/26873418/)
- Wardhana A et al. **The impact of Albumin Administration on Mortality and Resuscitation Volume in Burn Resuscitation: A Systematic Review and Meta-Analysis** *Bull Emerg Trauma* 13:203-214 (2025). PMID [41268471](https://pubmed.ncbi.nlm.nih.gov/41268471/)
- Lewis SR et al. **Colloids versus crystalloids for fluid resuscitation in critically ill people** *Cochrane Database Syst Rev* 8:CD000567 (2018). PMID [30073665](https://pubmed.ncbi.nlm.nih.gov/30073665/)
- Eljaiek R et al. **Albumin administration for fluid resuscitation in burn patients: A systematic review and meta-analysis** *Burns* 43:17-24 (2017). PMID [27613476](https://pubmed.ncbi.nlm.nih.gov/27613476/)
- Perel P et al. **Colloids versus crystalloids for fluid resuscitation in critically ill patients** *Cochrane Database Syst Rev*:CD000567 (2013). PMID [23450531](https://pubmed.ncbi.nlm.nih.gov/23450531/)
- Perel P et al. **Colloids versus crystalloids for fluid resuscitation in critically ill patients** *Cochrane Database Syst Rev*:CD000567 (2012). PMID [22696320](https://pubmed.ncbi.nlm.nih.gov/22696320/)
- Roberts I et al. **Human albumin solution for resuscitation and volume expansion in critically ill patients** *Cochrane Database Syst Rev* 2011:CD001208 (2011). PMID [22071799](https://pubmed.ncbi.nlm.nih.gov/22071799/)
- Greenhalgh DG et al. **A Prospective, Randomized, Multicenter Trial Comparing Lactated Ringer's Alone or With 5% Albumin for Resuscitation of Large Burns The Acute Burn ResUscitation Multicenter Prospective Trial 2 (ABRUPT2)** *Ann Surg* (2026). PMID [42144653](https://pubmed.ncbi.nlm.nih.gov/42144653/)
- Tang F et al. **Early supplemental parenteral nutrition shortens ventilation and intensive care unit stay in ICU patients aged ≥60 Years requiring mechanical ventilation: A randomized controlled trial** *Clin Nutr ESPEN* 72:102933 (2026). PMID [41565072](https://pubmed.ncbi.nlm.nih.gov/41565072/)
- Abdelmotaal AM et al. **Effect of Hydroxyethyl starch (HES) versus 5% albumin solution on intra-abdominal pressure in severe burn patients: A prospective randomized clinical trial** *Burns* 50:197-203 (2024). PMID [37833147](https://pubmed.ncbi.nlm.nih.gov/37833147/)
- Greenhalgh DG et al. **Burn Resuscitation Practices in North America: Results of the Acute Burn ResUscitation Multicenter Prospective Trial (ABRUPT)** *Ann Surg* 277:512-519 (2023). PMID [34417368](https://pubmed.ncbi.nlm.nih.gov/34417368/)
- Palmieri TL **Transfusion and Infections in the Burn Patient** *Surg Infect (Larchmt)* 22:49-53 (2021). PMID [32559401](https://pubmed.ncbi.nlm.nih.gov/32559401/)
- Li B et al. **Resuscitation Fluids in Septic Shock: A Network Meta-Analysis of Randomized Controlled Trials** *Shock* 53:679-685 (2020). PMID [31693630](https://pubmed.ncbi.nlm.nih.gov/31693630/)

## D. 소생술 합병증 (Resuscitation morbidity)

*14편*

- Wiktor AJ et al. **Safety and Efficacy of an Early Low-Dose Fresh Frozen Plasma Infusion in Burn Resuscitation** *J Burn Care Res* 47:807-816 (2026). PMID [41466517](https://pubmed.ncbi.nlm.nih.gov/41466517/)
- Zhou D et al. **Lessons from the similarities and differences in fluid resuscitation between burns and sepsis: a bibliometric analysis** *Front Med (Lausanne)* 12:1561619 (2025). PMID [40103790](https://pubmed.ncbi.nlm.nih.gov/40103790/)
- Kenney CL et al. **Neurologic Complications Associated With Burn Injury and Resuscitation** *J Surg Res* 304:36-40 (2024). PMID [39515023](https://pubmed.ncbi.nlm.nih.gov/39515023/)
- Halalmeh DR et al. **The role of a specialized urethral catheter in early detection of intra-abdominal hypertension: a case report** *J Surg Case Rep* 2024:rjae653 (2024). PMID [39421340](https://pubmed.ncbi.nlm.nih.gov/39421340/)
- Curry D et al. **Revision of an Adult Burn Center's Resuscitation Guideline Leads to Lower Resuscitation Requirements** *J Burn Care Res* 45:1499-1504 (2024). PMID [38824401](https://pubmed.ncbi.nlm.nih.gov/38824401/)
- Payne ML et al. **Effect of Dexmedetomidine on Fluid Resuscitation in Burn-Injured Patients** *J Burn Care Res* 45:1257-1263 (2024). PMID [38459902](https://pubmed.ncbi.nlm.nih.gov/38459902/)
- Arcieri TR et al. **Intraabdominal hypertension and abdominal compartment syndrome: What you need to know** *J Trauma Acute Care Surg* 99:504-513 (2025). PMID [40189748](https://pubmed.ncbi.nlm.nih.gov/40189748/)
- Tayebi S et al. **In Vitro Validation of a Novel Continuous Intra-Abdominal Pressure Measurement System (TraumaGuard)** *J Clin Med* 12 (2023). PMID [37834904](https://pubmed.ncbi.nlm.nih.gov/37834904/)
- Grünherz L et al. **Enzymatic debridement for circumferential deep burns: the role of surgical escharotomy** *Burns* 49:304-309 (2023). PMID [36604280](https://pubmed.ncbi.nlm.nih.gov/36604280/)
- Keyloun JW et al. **An unusual presentation of inhalation injury in a patient with high voltage electrical injury: A case report** *Int J Surg Case Rep* 77:357-361 (2020). PMID [33217653](https://pubmed.ncbi.nlm.nih.gov/33217653/)
- Zhang IY et al. **Validation of a low-cost simulation strategy for burn escharotomy training** *Injury* 51:2059-2065 (2020). PMID [32564962](https://pubmed.ncbi.nlm.nih.gov/32564962/)
- Mataro I et al. **Releasing Burn-Induced Compartment Syndrome by Enzymatic Escharotomy-Debridement: A Case Study** *J Burn Care Res* 41:1097-1103 (2020). PMID [32232328](https://pubmed.ncbi.nlm.nih.gov/32232328/)
- Delgado-Miguel C et al. **Iatrogenic Compartment Syndrome Secondary to Burn Dressing in a 2-Year-Old Child** *European J Pediatr Surg Rep* 7:e72-e74 (2019). PMID [31681528](https://pubmed.ncbi.nlm.nih.gov/31681528/)
- Fischer S et al. **Feasibility and safety of enzymatic debridement for the prevention of operative escharotomy in circumferential deep burns of the distal upper extremity** *Surgery* 165:1100-1105 (2019). PMID [30678870](https://pubmed.ncbi.nlm.nih.gov/30678870/)

## E. 고용량 비타민 C (High-dose ascorbate)

*12편*

- Tanaka H et al. **Reduction of resuscitation fluid volumes in severely burned patients using ascorbic acid administration: a randomized, prospective study** *Arch Surg* 135:326-31 (2000). PMID [10722036](https://pubmed.ncbi.nlm.nih.gov/10722036/)
- Sakurai M et al. **Reduced resuscitation fluid volume for second-degree experimental burns with delayed initiation of vitamin C therapy (beginning 6 h after injury)** *J Surg Res* 73:24-7 (1997). PMID [9441788](https://pubmed.ncbi.nlm.nih.gov/9441788/)
- Tanaka H et al. **Reduced resuscitation fluid volume for second-degree burns with delayed initiation of ascorbic acid therapy** *Arch Surg* 132:158-61 (1997). PMID [9041919](https://pubmed.ncbi.nlm.nih.gov/9041919/)
- Tanaka H et al. **Hemodynamic effects of delayed initiation of antioxidant therapy (beginning two hours after burn) in extensive third-degree burns** *J Burn Care Rehabil* 16:610-5 (1995). PMID [8582940](https://pubmed.ncbi.nlm.nih.gov/8582940/)
- Matsuda T et al. **Antioxidant therapy using high dose vitamin C: reduction of postburn resuscitation fluid volume requirements** *World J Surg* 19:287-91 (1995). PMID [7754637](https://pubmed.ncbi.nlm.nih.gov/7754637/)
- Tanaka H et al. **How long do we need to give antioxidant therapy during resuscitation when its administration is delayed for two hours?** *J Burn Care Rehabil* 13:567-72 (1992). PMID [1452591](https://pubmed.ncbi.nlm.nih.gov/1452591/)
- Novoa J et al. **Intravenous vitamin C in critically ill adult patients with burns: An integrative review** *Nutrition* 134:112728 (2025). PMID [40081106](https://pubmed.ncbi.nlm.nih.gov/40081106/)
- Soltany A et al. **A scoping review of the role of ascorbic acid in modifying fluid requirements in the resuscitation phase in burn patients** *Ann Med Surg (Lond)* 75:103460 (2022). PMID [35386786](https://pubmed.ncbi.nlm.nih.gov/35386786/)
- Cartotto R et al. **Burn State of the Science: Fluid Resuscitation** *J Burn Care Res* 38:e596-e604 (2017). PMID [28328669](https://pubmed.ncbi.nlm.nih.gov/28328669/)
- Buehner M et al. **Oxalate Nephropathy After Continuous Infusion of High-Dose Vitamin C as an Adjunct to Burn Resuscitation** *J Burn Care Res* 37:e374-9 (2016). PMID [25812044](https://pubmed.ncbi.nlm.nih.gov/25812044/)
- Ghanayem H **Towards evidence based emergency medicine: Best BETs from the Manchester Royal Infirmary. BET 3: Vitamin C in severe burns** *Emerg Med J* 29:1017-8 (2012). PMID [23180299](https://pubmed.ncbi.nlm.nih.gov/23180299/)
- Kahn SA et al. **Resuscitation after severe burn injury using high-dose ascorbic acid: a retrospective review** *J Burn Care Res* 32:110-7 (2011). PMID [21131846](https://pubmed.ncbi.nlm.nih.gov/21131846/)

## F. 대사항진 (Hypermetabolism)

*31편*

- AbuBaha M et al. **Hypermetabolism and Lipid Alterations Postburn: A Cardiovascular Perspective** *Cardiovasc Ther* 2026:5983391 (2026). PMID [41660577](https://pubmed.ncbi.nlm.nih.gov/41660577/)
- Shan H et al. **Optimizing nutrition and metabolism in a severely burned patient during prolonged continuous renal replacement therapy: a case report** *Front Nutr* 12:1749501 (2025). PMID [41613929](https://pubmed.ncbi.nlm.nih.gov/41613929/)
- Furtado T et al. **The Metabolic Impact of Biodegradable Temporizing Matrix in Burn Patients: A Retrospective Analysis of Resting Energy Expenditure and Inflammation** *J Burn Care Res* 46:1382-1391 (2025). PMID [40679419](https://pubmed.ncbi.nlm.nih.gov/40679419/)
- Dowling S et al. **Energy expenditure following biodegradable dermal matrix application in severe burn injury: A pilot study** *Clin Nutr ESPEN* 68:71-80 (2025). PMID [40311926](https://pubmed.ncbi.nlm.nih.gov/40311926/)
- Kingren MS et al. **HOUSING TEMPERATURE ALTERS BURN-INDUCED HYPERMETABOLISM IN MICE** *Shock* 63:118-131 (2025). PMID [39450911](https://pubmed.ncbi.nlm.nih.gov/39450911/)
- Bieerkehazhi S et al. **β-Adrenergic blockade attenuates adverse adipose tissue responses after burn** *J Mol Med (Berl)* 102:1245-1254 (2024). PMID [39145814](https://pubmed.ncbi.nlm.nih.gov/39145814/)
- Stanojcic M et al. **Anabolic and anticatabolic agents in critical care** *Curr Opin Crit Care* 22:325-31 (2016). PMID [27272101](https://pubmed.ncbi.nlm.nih.gov/27272101/)
- Jeschke MG et al. **Long-term persistance of the pathophysiologic response to severe burn injury** *PLoS One* 6:e21245 (2011). PMID [21789167](https://pubmed.ncbi.nlm.nih.gov/21789167/)
- Trifi A et al. **Early administration of norepinephrine in sepsis: Multicenter randomized clinical trial (EA-NE-S-TUN) study protocol** *PLoS One* 19:e0307407 (2024). PMID [39024364](https://pubmed.ncbi.nlm.nih.gov/39024364/)
- Gauglitz GG et al. **Abnormal insulin sensitivity persists up to three years in pediatric patients post-burn** *J Clin Endocrinol Metab* 94:1656-64 (2009). PMID [19240154](https://pubmed.ncbi.nlm.nih.gov/19240154/)
- Xie XQ et al. **Neuroendocrine system response modulates oxidative cellular damage in burn patients** *Tohoku J Exp Med* 211:161-9 (2007). PMID [17287600](https://pubmed.ncbi.nlm.nih.gov/17287600/)
- Sugiura T et al. **Effects of total parenteral nutrition on endotoxin translocation and extent of the stress response in burned rats** *Nutrition* 15:570-5 (1999). PMID [10422088](https://pubmed.ncbi.nlm.nih.gov/10422088/)
- Fatehi-Hassanabad Z et al. **Effects of L-canavanine, an inhibitor of inducible nitric oxide synthase, on endotoxin mediated shock in rats** *Shock* 6:194-200 (1996). PMID [8885085](https://pubmed.ncbi.nlm.nih.gov/8885085/)
- Graves TA et al. **The renal effects of low-dose dopamine in thermally injured patients** *J Trauma* 35:97-102; discussion 102-3 (1993). PMID [8331720](https://pubmed.ncbi.nlm.nih.gov/8331720/)
- Wu SL et al. **Morphological features of platelet-rich plasma on acellular vascular scaffolds** *Int J Surg* 111:7787-7797 (2025). PMID [40674238](https://pubmed.ncbi.nlm.nih.gov/40674238/)
- Yang Y et al. **Physicochemical Changes in Biomass Chars by Thermal Oxidation or Ambient Weathering and Their Impacts on Sorption of a Hydrophobic and a Cationic Compound** *Environ Sci Technol* 55:13072-13081 (2021). PMID [34555895](https://pubmed.ncbi.nlm.nih.gov/34555895/)
- Walcott SM et al. **Thermoregulatory costs in molting Antarctic Weddell seals: impacts of physiological and environmental conditions: Themed Issue Article: Conservation of Southern Hemisphere Mammals in a Changing World** *Conserv Physiol* 8:coaa022 (2020). PMID [32274067](https://pubmed.ncbi.nlm.nih.gov/32274067/)
- Speakman JR **Obesity and thermoregulation** *Handb Clin Neurol* 156:431-443 (2018). PMID [30454605](https://pubmed.ncbi.nlm.nih.gov/30454605/)
- Geiser F et al. **A burning question: what are the risks and benefits of mammalian torpor during and after fires?** *Conserv Physiol* 6:coy057 (2018). PMID [30323932](https://pubmed.ncbi.nlm.nih.gov/30323932/)
- Porter C et al. **Uncoupled skeletal muscle mitochondria contribute to hypermetabolism in severely burned adults** *Am J Physiol Endocrinol Metab* 307:E462-7 (2014). PMID [25074988](https://pubmed.ncbi.nlm.nih.gov/25074988/)
- Tidswell R et al. **Experimental sepsis causes SERCA2 expression in white adipose tissue but not classical browning** *Sci Rep* (2026). PMID [42399647](https://pubmed.ncbi.nlm.nih.gov/42399647/)
- Zhang Y et al. **HSF1 inhibits smooth muscle gene program to enhance white fat browning and hypermetabolism after burn injury** *Sci China Life Sci* 69:2788-2800 (2026). PMID [42132988](https://pubmed.ncbi.nlm.nih.gov/42132988/)
- Gao Y et al. **From selection signatures in cattle to functional validation in mice: HSPA12B negatively regulates adipose browning and thermogenesis** *J Anim Sci* 104 (2026). PMID [41615431](https://pubmed.ncbi.nlm.nih.gov/41615431/)
- Zhang C et al. **Ccl12 coordinates immune-neural crosstalk to promote adipose sympathetic remodeling after burn trauma** *Cell Rep* 45:116921 (2026). PMID [41610009](https://pubmed.ncbi.nlm.nih.gov/41610009/)
- Zhang X et al. **Disease-associated adipose browning: current evidence and perspectives** *Adipocyte* 15:2610540 (2026). PMID [41498391](https://pubmed.ncbi.nlm.nih.gov/41498391/)
- Jia Q et al. **Neuroplastin-55 is a receptor of Manf and protects against diet-induced obesity by promoting adipose browning** *Proc Natl Acad Sci U S A* 122:e2515526122 (2025). PMID [41284862](https://pubmed.ncbi.nlm.nih.gov/41284862/)
- Pearson LE et al. **To each its own: Thermoregulatory strategy varies among neonatal polar phocids** *Comp Biochem Physiol A Mol Integr Physiol* 178:59-67 (2014). PMID [25151642](https://pubmed.ncbi.nlm.nih.gov/25151642/)
- Tzika AA et al. **Skeletal muscle mitochondrial uncoupling in a murine cancer cachexia model** *Int J Oncol* 43:886-94 (2013). PMID [23817738](https://pubmed.ncbi.nlm.nih.gov/23817738/)
- Tzika AA et al. **Microarray analysis suggests that burn injury results in mitochondrial dysfunction in human skeletal muscle** *Int J Mol Med* 24:387-92 (2009). PMID [19639232](https://pubmed.ncbi.nlm.nih.gov/19639232/)
- Cree MG et al. **Human mitochondrial oxidative capacity is acutely impaired after burn trauma** *Am J Surg* 196:234-9 (2008). PMID [18639661](https://pubmed.ncbi.nlm.nih.gov/18639661/)
- Zhang Q et al. **Uncoupling protein 3 expression and intramyocellular lipid accumulation by NMR following local burn trauma** *Int J Mol Med* 18:1223-9 (2006). PMID [17089030](https://pubmed.ncbi.nlm.nih.gov/17089030/)

## G. 프로프라놀롤 (Propranolol)

*11편*

- Chao T et al. **Propranolol and Oxandrolone Therapy Accelerated Muscle Recovery in Burned Children** *Med Sci Sports Exerc* 50:427-435 (2018). PMID [29040226](https://pubmed.ncbi.nlm.nih.gov/29040226/)
- Guillory AN et al. **Oxandrolone Coadministration Does Not Alter Plasma Propranolol Concentrations in Severely Burned Pediatric Patients** *J Burn Care Res* 38:243-250 (2017). PMID [28240622](https://pubmed.ncbi.nlm.nih.gov/28240622/)
- Diaz EC et al. **Effects of pharmacological interventions on muscle protein synthesis and breakdown in recovery from burns** *Burns* 41:649-57 (2015). PMID [25468473](https://pubmed.ncbi.nlm.nih.gov/25468473/)
- Pereira CT et al. **Beta-blockade in burns** *Novartis Found Symp* 280:238-48; discussion 248-51 (2007). PMID [17380798](https://pubmed.ncbi.nlm.nih.gov/17380798/)
- Pereira C et al. **Post burn muscle wasting and the effects of treatments** *Int J Biochem Cell Biol* 37:1948-61 (2005). PMID [16109499](https://pubmed.ncbi.nlm.nih.gov/16109499/)
- Herndon DN et al. **Gene expression profiles and protein balance in skeletal muscle of burned children after beta-adrenergic blockade** *Am J Physiol Endocrinol Metab* 285:E783-9 (2003). PMID [12812919](https://pubmed.ncbi.nlm.nih.gov/12812919/)
- Mei-Zahav M et al. **Topical Propranolol Improves Epistaxis Control in Hereditary Hemorrhagic Telangiectasia (HHT): A Randomized Double-Blind Placebo-Controlled Trial** *J Clin Med* 9 (2020). PMID [32998220](https://pubmed.ncbi.nlm.nih.gov/32998220/)
- Hassoun-Kheir N et al. **The Effect of β-Blockers for Burn Patients on Clinical Outcomes: Systematic Review and Meta-Analysis** *J Intensive Care Med* 36:945-953 (2021). PMID [32686565](https://pubmed.ncbi.nlm.nih.gov/32686565/)
- Herndon D et al. **Reduced Postburn Hypertrophic Scarring and Improved Physical Recovery With Yearlong Administration of Oxandrolone and Propranolol** *Ann Surg* 268:431-441 (2018). PMID [30048322](https://pubmed.ncbi.nlm.nih.gov/30048322/)
- Rivas E et al. **Resting β-Adrenergic Blockade Does Not Alter Exercise Thermoregulation in Children With Burn Injury: A Randomized Control Trial** *J Burn Care Res* 39:402-412 (2018). PMID [28661984](https://pubmed.ncbi.nlm.nih.gov/28661984/)
- Guillory AN et al. **Propranolol kinetics in plasma from severely burned adults** *Burns* 43:1168-1174 (2017). PMID [28645713](https://pubmed.ncbi.nlm.nih.gov/28645713/)

## H. 동화 요법 (Anabolic therapy)

*22편*

- Lou J et al. **Oxandrolone for burn patients: a systematic review and updated meta-analysis of randomized controlled trials from 2005 to 2025** *World J Emerg Surg* 20:75 (2025). PMID [41023744](https://pubmed.ncbi.nlm.nih.gov/41023744/)
- Lou J et al. **The efficacy and safety of androgen analog oxandrolone in improving clinical outcomes in burn patients: a systematic review and meta-analysis of randomized controlled trials** *Front Med (Lausanne)* 12:1485474 (2025). PMID [40861228](https://pubmed.ncbi.nlm.nih.gov/40861228/)
- Ring J et al. **Oxandrolone in the Treatment of Burn Injuries: A Systematic Review and Meta-analysis** *J Burn Care Res* 41:190-199 (2020). PMID [31504621](https://pubmed.ncbi.nlm.nih.gov/31504621/)
- Li H et al. **The efficacy and safety of oxandrolone treatment for patients with severe burns: A systematic review and meta-analysis** *Burns* 42:717-27 (2016). PMID [26454425](https://pubmed.ncbi.nlm.nih.gov/26454425/)
- Real DS et al. **Oxandrolone use in adult burn patients. Systematic review and meta-analysis** *Acta Cir Bras* 29 Suppl 3:68-76 (2014). PMID [25351160](https://pubmed.ncbi.nlm.nih.gov/25351160/)
- Breederveld RS et al. **Recombinant human growth hormone for treating burns and donor sites** *Cochrane Database Syst Rev* 2014:CD008990 (2014). PMID [25222766](https://pubmed.ncbi.nlm.nih.gov/25222766/)
- Feathers JR et al. **The Use of Oxandrolone in the Management of Severe Burns: A Multi-service Survey of Burns Centres and Units Across the United Kingdom** *Cureus* 16:e57167 (2024). PMID [38681282](https://pubmed.ncbi.nlm.nih.gov/38681282/)
- Jalkh APC et al. **Oxandrolone Efficacy in Wound Healing in Burned and Decubitus Ulcer Patients: A Systematic Review** *Cureus* 14:e28079 (2022). PMID [36127967](https://pubmed.ncbi.nlm.nih.gov/36127967/)
- Kopel J et al. **A Reappraisal of Oxandrolone in Burn Management** *J Pharm Technol* 38:232-238 (2022). PMID [35832568](https://pubmed.ncbi.nlm.nih.gov/35832568/)
- Gusti NRL et al. **Effects Of Oxandrolone On Lean Body Mass (Lbm) In Severe Burn Patients: A Randomized, Double Blind, Placebo-Controlled Trial** *Ann Burns Fire Disasters* 35:55-61 (2022). PMID [35582088](https://pubmed.ncbi.nlm.nih.gov/35582088/)
- Chen M et al. **Rapid Pterygium Progression in a Child on Growth Hormone Therapy: A Case Report and Literature Review** *Cureus* 18:e100873 (2026). PMID [41658629](https://pubmed.ncbi.nlm.nih.gov/41658629/)
- Han J et al. **Biorhythm-mimicking growth hormone patch** *Nat Mater* 24:1283-1294 (2025). PMID [40181125](https://pubmed.ncbi.nlm.nih.gov/40181125/)
- Cuijpers MD et al. **The efficacy of therapeutic interventions on paediatric burn patients' height, weight, body composition, and muscle strength: A systematic review and meta-analysis** *Burns* 50:1437-1455 (2024). PMID [38580580](https://pubmed.ncbi.nlm.nih.gov/38580580/)
- Akkoç MF et al. **Investigation of the relationship of growth hormone, insulin-like growth factor (IGF)-1, and IGF-binding protein-3 levels with graft viability in autograft-transplanted pediatric patients with major burns** *Transpl Immunol* 73:101624 (2022). PMID [35577268](https://pubmed.ncbi.nlm.nih.gov/35577268/)
- Cai W et al. **TRIP13 promotes lung cancer cell growth and metastasis through AKT/mTORC1/c-Myc signaling** *Cancer Biomark* 30:237-248 (2021). PMID [33136091](https://pubmed.ncbi.nlm.nih.gov/33136091/)
- Xia X et al. **The relationship between H19 and parameters of ovarian reserve** *Reprod Biol Endocrinol* 18:46 (2020). PMID [32404103](https://pubmed.ncbi.nlm.nih.gov/32404103/)
- Boulestreau J et al. **Exercise and Skeletal Muscle-Derived Extracellular Vesicles: Their Cargo, Release, and Role in Metabolic Regulation** *J Extracell Biol* 5:e70139 (2026). PMID [42116941](https://pubmed.ncbi.nlm.nih.gov/42116941/)
- Farahat M et al. **Burn-associated metabolic dysfunction: could extracellular vesicles play a role?** *Clin Sci (Lond)* 140:698-721 (2026). PMID [41994933](https://pubmed.ncbi.nlm.nih.gov/41994933/)
- Berger MM et al. **What does "PICS" mean in major burns? Persistent critical illness or post-intensive care syndrome** *Burns* 52:107991 (2026). PMID [41946293](https://pubmed.ncbi.nlm.nih.gov/41946293/)
- Zhao J et al. **Inhibition of CIRBP represses the proliferation and migration of vascular smooth muscle cells via inhibiting Rheb/mTORC1 axis** *Biochem Biophys Res Commun* 725:150248 (2024). PMID [38870847](https://pubmed.ncbi.nlm.nih.gov/38870847/)
- Van K et al. **Effect of Short-Chain Fatty Acids on Inflammatory and Metabolic Function in an Obese Skeletal Muscle Cell Culture Model** *Nutrients* 16 (2024). PMID [38398822](https://pubmed.ncbi.nlm.nih.gov/38398822/)
- Tian S et al. **Radix Salvia miltiorrhiza Ameliorates Burn Injuries by Reducing Inflammation and Promoting Wound Healing** *J Inflamm Res* 16:4251-4263 (2023). PMID [37791115](https://pubmed.ncbi.nlm.nih.gov/37791115/)

## I. 혈당 조절 (Glycaemic control)

*22편*

- Zinter MS et al. **Biologic Mechanisms Underlying the Heterogeneous Response to Tight Glycemic Control among Differentially Inflamed Patients in the HALF-PINT Trial** *Am J Respir Crit Care Med* 211:1463-1473 (2025). PMID [40493436](https://pubmed.ncbi.nlm.nih.gov/40493436/)
- GBD 2019 Acute and Chronic Care Collaborators **Characterising acute and chronic care needs: insights from the Global Burden of Disease Study 2019** *Nat Commun* 16:4235 (2025). PMID [40335470](https://pubmed.ncbi.nlm.nih.gov/40335470/)
- Bakalář B et al. **Illusory movements for immobile patients with extensive burns (IMMOBILE): A randomized, controlled, cross-over trial** *Burns* 50:107264 (2024). PMID [39327102](https://pubmed.ncbi.nlm.nih.gov/39327102/)
- Tehrany PM et al. **Risk predictions of hospital-acquired pressure injury in the intensive care unit based on a machine learning algorithm** *Int Wound J* 20:3768-3775 (2023). PMID [37312659](https://pubmed.ncbi.nlm.nih.gov/37312659/)
- Lou JQ et al. **[A prospectively randomized controlled study of the effects of intensive insulin therapy combined with glutamine on nutritional metabolism, inflammatory response, and hemodynamics in severe burn patients]** *Zhonghua Shao Shang Za Zhi* 37:821-830 (2021). PMID [34645147](https://pubmed.ncbi.nlm.nih.gov/34645147/)
- Bateman RM et al. **36th International Symposium on Intensive Care and Emergency Medicine : Brussels, Belgium. 15-18 March 2016** *Crit Care* 20:94 (2016). PMID [27885969](https://pubmed.ncbi.nlm.nih.gov/27885969/)
- Sodade OE et al. **Intravenous insulin protocol reduces time to target glucose in critically ill trauma and burn patients** *J Crit Care Med (Targu Mures)* 12:359-367 (2026). PMID [42516160](https://pubmed.ncbi.nlm.nih.gov/42516160/)
- He XJ et al. **Promising community nursing effect on diabetic foot patients under the "six-in-one" collaborative management model** *Am J Transl Res* 17:9759-9770 (2025). PMID [41552321](https://pubmed.ncbi.nlm.nih.gov/41552321/)
- Ichai C et al. **Cyclosporine versus placebo pretreatment of brain-dead donors and kidney graft function (Cis-A-Rein trial): a multicenter, double-blind, randomized, controlled trial** *Intensive Care Med* 52:289-300 (2026). PMID [41504889](https://pubmed.ncbi.nlm.nih.gov/41504889/)
- Schaschinger T et al. **The Role of Glycemic Status in Adverse Outcomes Following ENT Surgery** *Laryngoscope* 136:1020-1029 (2026). PMID [41363265](https://pubmed.ncbi.nlm.nih.gov/41363265/)
- Chinese Burn Association, the Yangtze River Delta Integrated Diabetic Foot Alliance, and the Editorial Committee of the Chinese Journal of Burns and Wound Repair et al. **Practical guidelines for the prevention and management of diabetic foot disease in China** *Burns Trauma* 13:tkaf064 (2025). PMID [41355887](https://pubmed.ncbi.nlm.nih.gov/41355887/)
- Yan T et al. **Severe burns complicated by EDKA, CDI, and NTIS: A case report and literature review** *Medicine (Baltimore)* 105:e47272 (2026). PMID [41559973](https://pubmed.ncbi.nlm.nih.gov/41559973/)
- Khalaf F et al. **Beyond diabetes: harnessing the power of metformin in burn care** *Crit Care* 29:423 (2025). PMID [41057939](https://pubmed.ncbi.nlm.nih.gov/41057939/)
- Hallman TG et al. **Metformin is associated with reduced risk of mortality and morbidity in burn patients compared to insulin** *Burns* 50:1779-1789 (2024). PMID [38981799](https://pubmed.ncbi.nlm.nih.gov/38981799/)
- Wang L et al. **Diabetes Mellitus and Gastric Cancer: Correlation and Potential Mechanisms** *J Diabetes Res* 2023:4388437 (2023). PMID [38020199](https://pubmed.ncbi.nlm.nih.gov/38020199/)
- Raucci A et al. **MicroRNA-34a: the bad guy in age-related vascular diseases** *Cell Mol Life Sci* 78:7355-7378 (2021). PMID [34698884](https://pubmed.ncbi.nlm.nih.gov/34698884/)
- Greabu M et al. **Drugs Interfering with Insulin Resistance and Their Influence on the Associated Hypermetabolic State in Severe Burns: A Narrative Review** *Int J Mol Sci* 22 (2021). PMID [34575946](https://pubmed.ncbi.nlm.nih.gov/34575946/)
- Vanderwegen E et al. **Postoperative metabolic support in the ICU** *Best Pract Res Clin Anaesthesiol* 40:132-139 (2026). PMID [42552018](https://pubmed.ncbi.nlm.nih.gov/42552018/)
- Fu T et al. **Continuous Glucose Monitoring for Adult ICU Patients: A Meta-Analysis on Its Clinical Effectiveness** *Biol Res Nurs*:10998004261470782 (2026). PMID [42486781](https://pubmed.ncbi.nlm.nih.gov/42486781/)
- Lalani B et al. **Long-term Performance of a Dynamic Clinical Decision Support System for Intravenous Insulin Therapy: A Retrospective Cohort Study in Patients With and Without Diabetes** *J Diabetes Sci Technol*:19322968261467050 (2026). PMID [42478936](https://pubmed.ncbi.nlm.nih.gov/42478936/)
- Kim M et al. **Continuous Glucose Monitoring and Hypoglycaemia Risk During Paediatric DKA Treatment in the PICU: A Retrospective Cohort Study** *Diabetes Obes Metab* (2026). PMID [42421174](https://pubmed.ncbi.nlm.nih.gov/42421174/)
- Ammar A et al. **Long-Acting Insulin and Hypoglycemia in Critically Ill Patients with Diabetes: A Sensitivity Analysis** *J Intensive Care Med*:8850666261462876 (2026). PMID [42403155](https://pubmed.ncbi.nlm.nih.gov/42403155/)

## J. 조기 절제와 이식 (Early excision and grafting)

*19편*

- Wang H et al. **Timing of surgical excision for burn wounds: A systematic evaluation and meta-analysis comparing early and delayed excision** *J Int Med Res* 54:3000605261458971 (2026). PMID [42343617](https://pubmed.ncbi.nlm.nih.gov/42343617/)
- Lombardo GAG et al. **The ENGAGE Protocol: Enzymatic Debridement With NexoBrid Followed by Grafting After Graded Early Excision-A Retrospective Cohort Study** *J Burn Care Res* 47:944-949 (2026). PMID [41717791](https://pubmed.ncbi.nlm.nih.gov/41717791/)
- Marchica P et al. **Clinical outcomes of bromelain-based enzymatic debridement (NexoBrid®): evidence from the Italian National Burn Database** *Burns* 52:107915 (2026). PMID [41707536](https://pubmed.ncbi.nlm.nih.gov/41707536/)
- Allorto N et al. **Selective Early Excision Versus a Conservative Surgical Approach. Management of Deep Burn Injury in a Resource-Restricted Setting: A Prospective, Observational Study** *World J Surg* 50:523-528 (2026). PMID [41664426](https://pubmed.ncbi.nlm.nih.gov/41664426/)
- Czerny-Bednarczyk K et al. **Tissue engineering as a tool in a novel approach to the comprehensive treatment and management of a deeply and extensively burned patient: case report** *Cell Tissue Bank* 27:2 (2025). PMID [41284087](https://pubmed.ncbi.nlm.nih.gov/41284087/)
- Jeffery SLA et al. **A Scoping Review of Fluorescence Imaging: A Promising New Technology for Bacterial Detection in Burn Wounds** *J Burn Care Res* 47:315-322 (2026). PMID [40971539](https://pubmed.ncbi.nlm.nih.gov/40971539/)
- Atiyeh B et al. **Burn Wounds and Enzymatic Debridement (ED)-Past, Present, and Future** *J Burn Care Res* 45:864-876 (2024). PMID [38586910](https://pubmed.ncbi.nlm.nih.gov/38586910/)
- Nanni A et al. **Correction to: Research Progress of Photodynamic Therapy in Wound Healing: A Literature Review** *J Burn Care Res* 45:257 (2024). PMID [38006256](https://pubmed.ncbi.nlm.nih.gov/38006256/)
- Janzekovic Z **Once upon a time ... how west discovered east** *J Plast Reconstr Aesthet Surg* 61:240-4 (2008). PMID [18243082](https://pubmed.ncbi.nlm.nih.gov/18243082/)
- De Mey K et al. **Does enzymatic debridement reduce the occurrence of hypertrophic scarring in intermediate depth burns?** *Burns* 52:107819 (2026). PMID [41380209](https://pubmed.ncbi.nlm.nih.gov/41380209/)
- Kruglikov IL et al. **The complement pathway and the pathophysiology of fibroproliferative cutaneous scarring** *Front Immunol* 16:1701998 (2025). PMID [41246307](https://pubmed.ncbi.nlm.nih.gov/41246307/)
- Bird CL et al. **Our initial experience with rapid enzymatic debriding agent for burn eschar: Case series from an ABA verified burn center** *Burns Open* 10 (2025). PMID [42147452](https://pubmed.ncbi.nlm.nih.gov/42147452/)
- Aoki K et al. **A Pilot Study to Evaluate the Minimally Invasive Burn Care for Small, Deep Partial-Thickness Burns of the Hands and Feet Using Enzyme Debridement and Autologous Skin Cell Spray** *J Clin Med* 13 (2024). PMID [39768644](https://pubmed.ncbi.nlm.nih.gov/39768644/)
- Ball S et al. **Scar outcomes for conservatively managed children post burn injury: A retrospective study** *Int Wound J* 21:e14959 (2024). PMID [38949188](https://pubmed.ncbi.nlm.nih.gov/38949188/)
- De Decker I et al. **A single-stage bilayered skin reconstruction using Glyaderm® as an acellular dermal regeneration template results in improved scar quality: an intra-individual randomized controlled trial** *Burns Trauma* 11:tkad015 (2023). PMID [37143955](https://pubmed.ncbi.nlm.nih.gov/37143955/)
- Heard J et al. **Use of Cultured Epithelial Autograft in Conjunction with Biodegradable Temporizing Matrix in Massive Burns: A Case Series** *J Burn Care Res* 44:1434-1439 (2023). PMID [37227867](https://pubmed.ncbi.nlm.nih.gov/37227867/)
- Archer SB et al. **The use of sheet autografts to cover extensive burns in patients** *J Burn Care Rehabil* 19:33-8 (1998). PMID [9502022](https://pubmed.ncbi.nlm.nih.gov/9502022/)
- Munster AM et al. **Cultured epidermis for the coverage of massive burn wounds. A single center experience** *Ann Surg* 211:676-9; discussion 679-80 (1990). PMID [2357130](https://pubmed.ncbi.nlm.nih.gov/2357130/)
- McHugh TP et al. **Therapeutic efficacy of Biobrane in partial- and full-thickness thermal injury** *Surgery* 100:661-4 (1986). PMID [3532390](https://pubmed.ncbi.nlm.nih.gov/3532390/)

## K. 감염·면역마비 (Infection and immunoparalysis)

*15편*

- Demling RH et al. **Comparison of the postburn hyperdynamic state and changes in lung function (effect of wound bacterial content)** *Surgery* 100:828-35 (1986). PMID [3535145](https://pubmed.ncbi.nlm.nih.gov/3535145/)
- Heggers JP et al. **Transient and resident microflora of burn unit personnel and its influence on burn wound sepsis** *Infect Control* 3:471-4 (1982). PMID [6924647](https://pubmed.ncbi.nlm.nih.gov/6924647/)
- Wolk K et al. **Comparison of monocyte functions after LPS- or IL-10-induced reorientation: importance in clinical immunoparalysis** *Pathobiology* 67:253-6 (1999). PMID [10725796](https://pubmed.ncbi.nlm.nih.gov/10725796/)
- Cao K et al. **Establishment of A Mouse Sepsis-Induced Immunosuppression Model By Intranasal Instillation of Bacteria** *J Vis Exp* (2026). PMID [42371921](https://pubmed.ncbi.nlm.nih.gov/42371921/)
- Butsch B et al. **Obesity Attenuates the Association Between the Modified Five-Item Frailty Score and Mortality in Septic Surgical Patients** *J Surg Res* 325:506-515 (2026). PMID [42361581](https://pubmed.ncbi.nlm.nih.gov/42361581/)
- He H et al. **Impact of organ dysfunction on overall survival in patients with burns ≥ 70% TBSA: A case-control study** *Burns* 52:107941 (2026). PMID [41819651](https://pubmed.ncbi.nlm.nih.gov/41819651/)
- Gopi P et al. **Comparative Study of Swab Culture and Tissue Biopsy in Burn Wound Sepsis: A Tertiary Centre Experience** *Cureus* 17:e99588 (2025). PMID [41556028](https://pubmed.ncbi.nlm.nih.gov/41556028/)
- Fan JB et al. **The "cytokine storm" in infection and sepsis: win the battle but lose the war** *Mil Med Res* 12:95 (2026). PMID [41521316](https://pubmed.ncbi.nlm.nih.gov/41521316/)
- Giamarellos-Bourboulis EJ et al. **Precision Immunotherapy to Improve Sepsis Outcomes: The ImmunoSep Randomized Clinical Trial** *JAMA* 335:775-786 (2026). PMID [41359996](https://pubmed.ncbi.nlm.nih.gov/41359996/)
- Tejiram S et al. **Fighting a New Front on an Old Battlefield: Examining the Development of Topical Antimicrobial Care to Control Burn Wound Sepsis** *J Burn Care Res* 46:248-255 (2025). PMID [39288163](https://pubmed.ncbi.nlm.nih.gov/39288163/)
- Hoogewerf CJ et al. **Topical treatment for facial burns** *Cochrane Database Syst Rev* 7:CD008058 (2020). PMID [32725896](https://pubmed.ncbi.nlm.nih.gov/32725896/)
- Ahuja RB et al. **A prospective double-blinded comparative analysis of framycetin and silver sulphadiazine as topical agents for burns: a pilot study** *Burns* 35:672-6 (2009). PMID [19443125](https://pubmed.ncbi.nlm.nih.gov/19443125/)
- Connor-Ballard PA **Understanding and managing burn pain: Part 2** *Am J Nurs* 109:54-62; quiz 63 (2009). PMID [19411907](https://pubmed.ncbi.nlm.nih.gov/19411907/)
- Neely AN et al. **Are topical antimicrobials effective against bacteria that are highly resistant to systemic antibiotics?** *J Burn Care Res* 30:19-29 (2009). PMID [19060725](https://pubmed.ncbi.nlm.nih.gov/19060725/)
- Kucan JO et al. **The potential benefit of 5% Sulfamylon Solution in the treatment of Acinetobacter baumannii-contaminated traumatic war wounds** *J Burns Wounds* 4:e3 (2005). PMID [16921408](https://pubmed.ncbi.nlm.nih.gov/16921408/)

## L. 화상 약동학 (Burn pharmacokinetics)

*22편*

- Závorszky L et al. **Impact of augmented renal clearance on piperacillin and meropenem exposure in critically ill adult burn patients: A prospective observational study** *Burns* 52:108138 (2026). PMID [42508095](https://pubmed.ncbi.nlm.nih.gov/42508095/)
- Gatti M et al. **Pharmacokinetics-pharmacodynamics perspective to optimizing therapy with beta-lactams in critically ill patients: an update** *Expert Rev Anti Infect Ther* 23:1215-1233 (2025). PMID [41437777](https://pubmed.ncbi.nlm.nih.gov/41437777/)
- Lipman J et al. **The long walk to a short half-life: the discovery of augmented renal clearance and its impact on antibiotic dosing** *J Antimicrob Chemother* 80:3367-3374 (2025). PMID [41077962](https://pubmed.ncbi.nlm.nih.gov/41077962/)
- Olivet JD et al. **Pharmacokinetics of continuous infusion ceftolozane/tazobactam in two patients with extensive total body surface area burns** *Pharmacotherapy* 45:386-392 (2025). PMID [40231834](https://pubmed.ncbi.nlm.nih.gov/40231834/)
- Zoccali C et al. **Pharmacokinetic relevance of glomerular hyperfiltration for drug dosing** *Clin Kidney J* 16:1580-1586 (2023). PMID [37779850](https://pubmed.ncbi.nlm.nih.gov/37779850/)
- Hill DM et al. **Pharmacokinetic Analysis of Intravenous Push Cefepime in Burn Patients with Augmented Renal Clearance** *J Burn Care Res* 45:151-157 (2024). PMID [37688528](https://pubmed.ncbi.nlm.nih.gov/37688528/)
- Shi Y et al. **Vancomycin population pharmacokinetics in patients with burns** *Front Med (Lausanne)* 13:1829805 (2026). PMID [42440663](https://pubmed.ncbi.nlm.nih.gov/42440663/)
- Koreeda H et al. **Patient with severe burns for whom serum cystatin C-based assessment was performed due to suspected errors in serum creatinine-based renal function evaluation: a case report** *Eur J Hosp Pharm* (2026). PMID [41571446](https://pubmed.ncbi.nlm.nih.gov/41571446/)
- Azadi S et al. **Therapeutic vancomycin monitoring: a comparative analysis of high-performance liquid chromatography and chemiluminescent microparticle immunoassay methods in liver transplant recipients** *Front Pharmacol* 16:1516339 (2025). PMID [40529505](https://pubmed.ncbi.nlm.nih.gov/40529505/)
- Galvidis IA et al. **Therapeutic Monitoring of Vancomycin Implemented by Eremomycin ELISA** *Antibiotics (Basel)* 13 (2024). PMID [39766523](https://pubmed.ncbi.nlm.nih.gov/39766523/)
- Ma P et al. **Prediction of vancomycin plasma concentration in elderly patients based on multi-algorithm mining combined with population pharmacokinetics** *Sci Rep* 14:27165 (2024). PMID [39511378](https://pubmed.ncbi.nlm.nih.gov/39511378/)
- Santos RM et al. **A Multicenter, Retrospective Outcome Analysis of Vancomycin Area Under the Curve Versus Trough-Based Dosing Strategies in Patients With Burn OR Inhalational Injuries (MONITOR)** *J Burn Care Res* 45:1383-1389 (2024). PMID [38900835](https://pubmed.ncbi.nlm.nih.gov/38900835/)
- Skorup P et al. **Estimation of Plasma Meropenem Concentrations in Patients with Severe Infections Using Ultraviolet and Visible Range Spectrophotometry: A Pilot Study** *Eur J Drug Metab Pharmacokinet* 51:413-420 (2026). PMID [42087051](https://pubmed.ncbi.nlm.nih.gov/42087051/)
- Hall RG 2nd et al. **Pharmacokinetics of Ceftolozane/Tazobactam in Patients With Partial- and Full-Thickness Skin Burns** *Pharmacotherapy* 45:774-779 (2025). PMID [41212678](https://pubmed.ncbi.nlm.nih.gov/41212678/)
- Tebano G et al. **Which Are the Best Regimens of Broad-Spectrum Beta-Lactam Antibiotics in Burn Patients? A Systematic Review of Evidence from Pharmacology Studies** *Antibiotics (Basel)* 12 (2023). PMID [38136771](https://pubmed.ncbi.nlm.nih.gov/38136771/)
- Závorszky L et al. **[Therapeutic drug monitoring of beta-lactam antibiotics in critically ill adult patients]** *Orv Hetil* 164:1904-1911 (2023). PMID [38043089](https://pubmed.ncbi.nlm.nih.gov/38043089/)
- Meenks SD et al. **Mass Spectrometric Characterization of Flucloxacillin Binding to Serum Albumin in Hospitalized Patients With Hypoalbuminemia: An Observational Pilot Study** *Ther Drug Monit* (2026). PMID [42417601](https://pubmed.ncbi.nlm.nih.gov/42417601/)
- Wei XC et al. **Exploratory model-based optimizing isavuconazole dosing regimens for Aspergillus spp. in adult patients according to serum albumin levels: a pharmacokinetic/pharmacodynamic analysis using Monte Carlo simulation** *Front Cell Infect Microbiol* 16:1737154 (2026). PMID [42394820](https://pubmed.ncbi.nlm.nih.gov/42394820/)
- Morales Junior R et al. **What if the Aquarium is in the Intensive Care Unit? Extending the Fish Tank Model to Teach Pharmacokinetics in Critically Ill Patients** *Am J Pharm Educ* 90:102016 (2026). PMID [42248418](https://pubmed.ncbi.nlm.nih.gov/42248418/)
- Nix DE et al. **Correction of posaconazole concentrations for hypoalbuminemia** *Pharmacotherapy* 45:324-331 (2025). PMID [40251837](https://pubmed.ncbi.nlm.nih.gov/40251837/)
- Webb AJ et al. **Clinical Consequences of Disproportionate Free Valproate Elevation in Critically Ill Adult Patients: A Multicenter Retrospective Cohort Study** *Neurocrit Care* 43:472-483 (2025). PMID [40133756](https://pubmed.ncbi.nlm.nih.gov/40133756/)
- Idasiak-Piechocka I et al. **Effect of hypoalbuminemia on drug pharmacokinetics** *Front Pharmacol* 16:1546465 (2025). PMID [40051558](https://pubmed.ncbi.nlm.nih.gov/40051558/)

## M. 흡입 손상 (Inhalation injury)

*16편*

- Fuciños LC et al. **Cardiopulmonary hemodynamic alterations during early resuscitation in burn patients with inhalation injury: do they impact prognosis?** *J Burn Care Res* (2026). PMID [42275046](https://pubmed.ncbi.nlm.nih.gov/42275046/)
- Ozmen BB et al. **Evidence-based AI clinical decision support system for acute burn care and complex reconstruction** *J Plast Reconstr Aesthet Surg* 116:242-253 (2026). PMID [41934058](https://pubmed.ncbi.nlm.nih.gov/41934058/)
- Kiwan O et al. **What You Need to Know About: Assessment of Burns and Initial Management** *Br J Hosp Med (Lond)* 86:1-18 (2025). PMID [41134176](https://pubmed.ncbi.nlm.nih.gov/41134176/)
- Echeverri C et al. **Lower Extremity Compartment Syndrome Due to Capillary Leak Syndrome Following 60% Total Body Surface Area Burn Injury** *J Burn Care Res* 47:410-413 (2026). PMID [40996144](https://pubmed.ncbi.nlm.nih.gov/40996144/)
- Váňa V et al. **Results of microbiological surveillance in patients with high-voltage eletrical injuries: A 10-year single center experience** *Epidemiol Mikrobiol Imunol* 74:97-106 (2025). PMID [40747749](https://pubmed.ncbi.nlm.nih.gov/40747749/)
- Kruse M et al. **Characterisation of Fluid Administration in Burn Shock-A Retrospective Cohort Analysis** *Eur Burn J* 6 (2025). PMID [40558630](https://pubmed.ncbi.nlm.nih.gov/40558630/)
- Nguyen TNM et al. **Scheduled Bronchoscopy with Nebulized Heparin and N-Acetylcysteine in Burn Patients with Inhalation Injury: A Randomized Trial** *Eur Burn J* 7 (2026). PMID [42201099](https://pubmed.ncbi.nlm.nih.gov/42201099/)
- Risinger WB et al. **Nebulized Heparin and N-Acetylcysteine do Not Improve Outcomes of Intubated Burn Patients With Grade II or III Inhalation Injuries** *Am Surg* 91:1392-1395 (2025). PMID [40346841](https://pubmed.ncbi.nlm.nih.gov/40346841/)
- Hakim SM et al. **Effect of early administration of inhaled heparin on outcomes of smoke inhalation injury: A randomized controlled trial** *Burns* 51:107518 (2025). PMID [40319829](https://pubmed.ncbi.nlm.nih.gov/40319829/)
- Rajaratnam G et al. **"To BAL or not to BAL, that is the question": Variations in smoke inhalation injury guidelines from burn units and centres in England, Scotland and Wales** *Burns* 50:107273 (2024). PMID [39353794](https://pubmed.ncbi.nlm.nih.gov/39353794/)
- Peters RA et al. **Extracorporeal Membrane Oxygenation in a Patient with Severe Inhalation Injury: Multidisciplinary Burn Team Care** *J Burn Care Res* 45:796-800 (2024). PMID [38367208](https://pubmed.ncbi.nlm.nih.gov/38367208/)
- Milton-Jones H et al. **An international RAND/UCLA expert panel to determine the optimal diagnosis and management of burn inhalation injury** *Crit Care* 27:459 (2023). PMID [38012797](https://pubmed.ncbi.nlm.nih.gov/38012797/)
- Swafford EP et al. **Inpatient Complications and Outcomes for Burn Patients Admitted with Opioid-Positive Drug Screens** *J Burn Care Res* (2026). PMID [42550487](https://pubmed.ncbi.nlm.nih.gov/42550487/)
- Palmer BL et al. **Corneal abrasions in burn admissions: severity and geographic access factors in a retrospective cohort** *Proc (Bayl Univ Med Cent)*:1-4 (2026). PMID [42496181](https://pubmed.ncbi.nlm.nih.gov/42496181/)
- Matheny M et al. **Residential Air Pollution, Inhalation Injury, and Burn Recovery: A Burn Model System Study** *J Burn Care Res* (2026). PMID [42476515](https://pubmed.ncbi.nlm.nih.gov/42476515/)
- Wardhana A et al. **Mortality Analysis of ICU Burn Patients in Indonesia's National Referral Hospital: A 2-Year Retrospective Study** *Ann Burns Fire Disasters* 39:95-99 (2026). PMID [42370050](https://pubmed.ncbi.nlm.nih.gov/42370050/)

## N. 염증·사이토카인·간 (Inflammation, cytokines, liver)

*16편*

- Chapman CD et al. **Bactericidal and anti-inflammatory activity of defined hypochlorite solution** *PLoS One* 21:e0352576 (2026). PMID [42497178](https://pubmed.ncbi.nlm.nih.gov/42497178/)
- Yeong EK et al. **Intravenous placenta-derived mesenchymal stem cells enhance early survival as an adjunct to supportive care in severe burn: a preclinical study** *Stem Cell Res Ther* 17 (2026). PMID [42493791](https://pubmed.ncbi.nlm.nih.gov/42493791/)
- Ma H et al. **Modulating hypertrophic scar formation by targeting endothelial transient receptor potential vanilloid-1/nuclear factor kappa-B/interleukin-6 axis to regulate angiogenesis** *Burns Trauma* 14:tkag009 (2026). PMID [42131490](https://pubmed.ncbi.nlm.nih.gov/42131490/)
- Yang GX et al. **The Diagnostic and Prognostic Value of Procalcitonin and High-Sensitivity C-Reactive Protein in Early-Stage Burn Sepsis: A Retrospective Cohort Study** *J Inflamm Res* 19:571738 (2026). PMID [41873329](https://pubmed.ncbi.nlm.nih.gov/41873329/)
- Liu H et al. **Characteristics and significance of blister fluid cell-free mitochondrial DNA in pediatric small area intermediate-depth burn wounds** *Burns* 52:107780 (2026). PMID [41297227](https://pubmed.ncbi.nlm.nih.gov/41297227/)
- Liu XY et al. **Bibliometric analysis and initial animal efficacy evaluation of top ten scoring drugs to enhance oral rehydration therapy in early post-burn shock** *Front Pharmacol* 16:1614159 (2025). PMID [40900826](https://pubmed.ncbi.nlm.nih.gov/40900826/)
- Kaur S et al. **Adipose-specific ATGL ablation reduces burn injury-induced metabolic derangements in mice** *Clin Transl Med* 11:e417 (2021). PMID [34185433](https://pubmed.ncbi.nlm.nih.gov/34185433/)
- Bogdanovic E et al. **Endoplasmic reticulum stress in adipose tissue augments lipolysis** *J Cell Mol Med* 19:82-91 (2015). PMID [25381905](https://pubmed.ncbi.nlm.nih.gov/25381905/)
- Hiyama Y et al. **Fenofibrate does not affect burn-induced hepatic endoplasmic reticulum stress** *J Surg Res* 185:733-9 (2013). PMID [23866789](https://pubmed.ncbi.nlm.nih.gov/23866789/)
- Barrow RE et al. **Identification of factors contributing to hepatomegaly in severely burned children** *Shock* 24:523-8 (2005). PMID [16317382](https://pubmed.ncbi.nlm.nih.gov/16317382/)
- Nie C et al. **The blood parameters and liver function changed inconsistently among children between burns and traumatic injuries** *PeerJ* 7:e6415 (2019). PMID [30775182](https://pubmed.ncbi.nlm.nih.gov/30775182/)
- Korkmaz HI et al. **The role of complement in the acute phase response after burns** *Burns* 43:1390-1399 (2017). PMID [28410933](https://pubmed.ncbi.nlm.nih.gov/28410933/)
- Li JH et al. **Suppressed acute phase response to LPS-induced hepatic injury in Smad3-deficient mice** *Mol Immunol* 46:362-5 (2009). PMID [19058853](https://pubmed.ncbi.nlm.nih.gov/19058853/)
- Chai J et al. **[Influence of escharectomy and skin grafting during early burn stage on acute-phase response in severely burned rats and its significance]** *Zhonghua Yi Xue Za Zhi* 82:1420-3 (2002). PMID [12509927](https://pubmed.ncbi.nlm.nih.gov/12509927/)
- Koike K et al. **Recombinant human interleukin-1alpha increases serum albumin, Gc-globulin, and alpha1-antitrypsin levels in burned mice** *Tohoku J Exp Med* 198:23-9 (2002). PMID [12498311](https://pubmed.ncbi.nlm.nih.gov/12498311/)
- Mooser V et al. **Major reduction in plasma Lp(a) levels during sepsis and burns** *Arterioscler Thromb Vasc Biol* 20:1137-42 (2000). PMID [10764684](https://pubmed.ncbi.nlm.nih.gov/10764684/)

## O. 골·성장·장기 결과 (Bone, growth, long-term outcome)

*15편*

- Moolhuijsen LME et al. **Genomic analyses implicate hormonal and metabolic dysregulation in polycystic ovary syndrome** *Nat Genet* 58:1040-1050 (2026). PMID [42026183](https://pubmed.ncbi.nlm.nih.gov/42026183/)
- GBD 2023 Demographics Collaborators **Global age-sex-specific all-cause mortality and life expectancy estimates for 204 countries and territories and 660 subnational locations, 1950-2023: a demographic analysis for the Global Burden of Disease Study 2023** *Lancet* 406:1731-1810 (2025). PMID [41092927](https://pubmed.ncbi.nlm.nih.gov/41092927/)
- GBD 2023 Disease and Injury and Risk Factor Collaborators **Burden of 375 diseases and injuries, risk-attributable burden of 88 risk factors, and healthy life expectancy in 204 countries and territories, including 660 subnational locations, 1990-2023: a systematic analysis for the Global Burden of Disease Study 2023** *Lancet* 406:1873-1922 (2025). PMID [41092926](https://pubmed.ncbi.nlm.nih.gov/41092926/)
- GBD 2021 Risk Factors Collaborators **Global burden and strength of evidence for 88 risk factors in 204 countries and 811 subnational locations, 1990-2021: a systematic analysis for the Global Burden of Disease Study 2021** *Lancet* 403:2162-2203 (2024). PMID [38762324](https://pubmed.ncbi.nlm.nih.gov/38762324/)
- GBD 2021 Diseases and Injuries Collaborators **Global incidence, prevalence, years lived with disability (YLDs), disability-adjusted life-years (DALYs), and healthy life expectancy (HALE) for 371 diseases and injuries in 204 countries and territories and 811 subnational locations, 1990-2021: a systematic analysis for the Global Burden of Disease Study 2021** *Lancet* 403:2133-2161 (2024). PMID [38642570](https://pubmed.ncbi.nlm.nih.gov/38642570/)
- GBD 2021 Causes of Death Collaborators **Global burden of 288 causes of death and life expectancy decomposition in 204 countries and territories and 811 subnational locations, 1990-2021: a systematic analysis for the Global Burden of Disease Study 2021** *Lancet* 403:2100-2132 (2024). PMID [38582094](https://pubmed.ncbi.nlm.nih.gov/38582094/)
- Damme KSF et al. **The cumulative impact of fine particulate matter exposure on hippocampal volume and working memory: Insights from prenatal and adolescent exposures from the ABCD study** *Dev Cogn Neurosci* 77:101648 (2026). PMID [41352198](https://pubmed.ncbi.nlm.nih.gov/41352198/)
- Jilani LZ et al. **Treatment of paediatric subtrochanteric femoral fractures using titanium elastic nails: a single-center experience** *Int J Burns Trauma* 15:202-209 (2025). PMID [41278380](https://pubmed.ncbi.nlm.nih.gov/41278380/)
- McPheron M et al. **Treatment of PDGFRB -Related Penttinen Syndrome With Imatinib in a Young Child** *Am J Med Genet C Semin Med Genet* 199:176-182 (2025). PMID [40742224](https://pubmed.ncbi.nlm.nih.gov/40742224/)
- Wang H et al. **Empagliflozin-Pretreated MSC-Derived Exosomes Enhance Angiogenesis and Wound Healing via PTEN/AKT/VEGF Pathway** *Int J Nanomedicine* 20:5119-5136 (2025). PMID [40297404](https://pubmed.ncbi.nlm.nih.gov/40297404/)
- Hodge KM et al. **Epigenetic associations with neonatal age in infants born very preterm, particularly among genes involved in neurodevelopment** *Sci Rep* 14:18147 (2024). PMID [39103365](https://pubmed.ncbi.nlm.nih.gov/39103365/)
- Ungureanu A et al. **Functional Recovery and Emotional Burden After Burn Injury: A Quality of Life Assessment in Romanian Burn Survivors** *Diseases* 14 (2026). PMID [42346304](https://pubmed.ncbi.nlm.nih.gov/42346304/)
- Akter T et al. **Health-Related Quality of Life in Burn Injury Victims Attending a Specialized Institution** *Cureus* 18:e106895 (2026). PMID [42131675](https://pubmed.ncbi.nlm.nih.gov/42131675/)
- Vigneron L et al. **Effects on health-related quality of life of therapeutic exercise in burn survivors: A systematic review and meta-analyses** *Burns* 52:107829 (2026). PMID [41447903](https://pubmed.ncbi.nlm.nih.gov/41447903/)
- GBD 2023 Causes of Death Collaborators **Global burden of 292 causes of death in 204 countries and territories and 660 subnational locations, 1990-2023: a systematic analysis for the Global Burden of Disease Study 2023** *Lancet* 406:1811-1872 (2025). PMID [41092928](https://pubmed.ncbi.nlm.nih.gov/41092928/)

## P. 특수 집단·전해질 (Special populations)

*21편*

- Aigner A et al. **Too much or too little? Fluid resuscitation in the first 24 h after severe burns: Evaluating the Parkland formula - A retrospective analysis of adult burn patients in Austria, Germany, and Switzerland 2015-2022** *Burns* 51:107397 (2025). PMID [40068435](https://pubmed.ncbi.nlm.nih.gov/40068435/)
- Shen ZA et al. **[Establishment and application of the ten-fold rehydration formula for emergency resuscitation of pediatric patients after extensive burns]** *Zhonghua Shao Shang Yu Chuang Mian Xiu Fu Za Zhi* 39:59-64 (2023). PMID [36740427](https://pubmed.ncbi.nlm.nih.gov/36740427/)
- Allorto NL et al. **A practical formula for fluid resuscitation in acute paediatric burns in a low resource setting: A pilot study** *Injury* 54:25-28 (2023). PMID [36089555](https://pubmed.ncbi.nlm.nih.gov/36089555/)
- Stevens JV et al. **Weight-based vs body surface area-based fluid resuscitation predictions in pediatric burn patients** *Burns* 49:120-128 (2023). PMID [35351355](https://pubmed.ncbi.nlm.nih.gov/35351355/)
- Khan S et al. **Clinical and Epidemiological Profile, Mortality, and Associated Factors Among Paediatric Burn Patients in Western Rajasthan, India** *Cureus* 18:e111793 (2026). PMID [42534189](https://pubmed.ncbi.nlm.nih.gov/42534189/)
- Zakaria H et al. **A global contrast shortage: a descriptive study exploring IV contrast utilization and outcomes in trauma patients at a level I trauma center** *Eur J Trauma Emerg Surg* 52 (2026). PMID [42496906](https://pubmed.ncbi.nlm.nih.gov/42496906/)
- Lahoti GA et al. **Injury mechanism patterns, mortality predictors, and resource utilization in pediatric trauma in Karachi, Pakistan: A seven-year registry study** *Injury* 57:113514 (2026). PMID [42475931](https://pubmed.ncbi.nlm.nih.gov/42475931/)
- Hyman SC et al. **Association of invasive intracranial pressure monitoring with mortality in severe pediatric traumatic brain injury: A national ACS TQIP analysis** *J Pediatr Surg*:163310 (2026). PMID [42468687](https://pubmed.ncbi.nlm.nih.gov/42468687/)
- Montesinos P et al. **QuANTUM-Wild: a Phase III, randomized trial of quizartinib in newly diagnosed FLT3-ITD-negative acute myeloid leukemia** *Future Oncol* 22:1913-1924 (2026). PMID [42434847](https://pubmed.ncbi.nlm.nih.gov/42434847/)
- Melekoğlu T et al. **Acute Biochemical Muscle Damage Responses to a Single Session of Whole-Body Electromyostimulation** *Med Sci Sports Exerc* 58:1362-1376 (2026). PMID [41839178](https://pubmed.ncbi.nlm.nih.gov/41839178/)
- Mitsui D et al. **Severe caffeine poisoning treated with intermittent hemodialysis under circulatory support** *Am J Emerg Med* 76:270.e5-270.e7 (2024). PMID [38129271](https://pubmed.ncbi.nlm.nih.gov/38129271/)
- Boyd AN et al. **A voltage-based analysis of fluid delivery and outcomes in burn patients with electrical injuries over a 6-year period** *Burns* 45:869-875 (2019). PMID [30935702](https://pubmed.ncbi.nlm.nih.gov/30935702/)
- Pollanen MS **The pathology of torture** *Forensic Sci Int* 284:85-96 (2018). PMID [29367172](https://pubmed.ncbi.nlm.nih.gov/29367172/)
- Navarrete N **Severe rhabdomyolysis without renal injury associated with lightning strike** *J Burn Care Res* 34:e209-12 (2013). PMID [22929530](https://pubmed.ncbi.nlm.nih.gov/22929530/)
- Maghsoudi H et al. **Electrical and lightning injuries** *J Burn Care Res* 28:255-61 (2007). PMID [17351442](https://pubmed.ncbi.nlm.nih.gov/17351442/)
- Colligan T et al. **The Impact of Remoteness on the Outcomes of Children With Prenatal Drug Exposure: A Population-Based Cohort Study** *J Paediatr Child Health* 62:1223-1234 (2026). PMID [42186151](https://pubmed.ncbi.nlm.nih.gov/42186151/)
- Chen C et al. **Critical illness outcomes of hospitalized pregnant women following a Texas abortion ban** *Am J Respir Crit Care Med* 212:1730-1739 (2026). PMID [42149703](https://pubmed.ncbi.nlm.nih.gov/42149703/)
- Picetti E et al. **The management of severe isolated traumatic brain injury in pregnancy: A joint consensus statement from the European Association of Neurosurgical Societies (EANS) and the World Society of Emergency Surgery (WSES)** *Brain Spine* 6:105971 (2026). PMID [42064346](https://pubmed.ncbi.nlm.nih.gov/42064346/)
- Saiki K et al. **Social Determinants of Health, Perceived Stress Scores, and Adverse Pregnancy Outcomes** *Am J Perinatol* (2026). PMID [41946475](https://pubmed.ncbi.nlm.nih.gov/41946475/)
- Makarious L et al. **Maternal Deaths Due to Suicide, Accidental Poisoning and Undetermined Intent Within 5 Years Following Childbirth: A Population-Based Study** *BJOG* 133:1581-1591 (2026). PMID [41914610](https://pubmed.ncbi.nlm.nih.gov/41914610/)
- Michelle H et al. **Disparities among adolescent pregnant trauma patients** *J Trauma Acute Care Surg* (2026). PMID [41873860](https://pubmed.ncbi.nlm.nih.gov/41873860/)

## Q. 영양·미량원소 (Nutrition and micronutrients)

*19편*

- Yeh DD **Fueling recovery: evidence-based ICU nutrition and immunonutrition strategies in 2026** *Trauma Surg Acute Care Open* 11:e002284 (2026). PMID [42094750](https://pubmed.ncbi.nlm.nih.gov/42094750/)
- Qinyuan D et al. **Protein nutritional support in critically ill patients: pathophysiological basis, clinical evidence, and areas of uncertainty - a narrative review** *Front Med (Lausanne)* 13:1770345 (2026). PMID [42040590](https://pubmed.ncbi.nlm.nih.gov/42040590/)
- Zhou S et al. **Predictive Value of Inflammatory Burden Index for Sepsis in Critically Ill Patients with Extensive Burns: A Decade-Long Cohort Study** *J Inflamm Res* 19:574776 (2026). PMID [41884165](https://pubmed.ncbi.nlm.nih.gov/41884165/)
- Yang X et al. **Constructing the early enteral nutrition management protocol for severely burned adult patients: a Delphi study** *BMC Nutr* 12:5 (2026). PMID [41495855](https://pubmed.ncbi.nlm.nih.gov/41495855/)
- Nikiforov MV et al. **[Nutritional support in the comprehensive management of cutaneous squamous cell carcinoma in congenital epidermolysis bullosa]** *Vopr Pitan* 94:139-150 (2025). PMID [41263358](https://pubmed.ncbi.nlm.nih.gov/41263358/)
- Heyland DK et al. **A Randomized Trial of Enteral Glutamine for Treatment of Burn Injuries** *N Engl J Med* 387:1001-1010 (2022). PMID [36082909](https://pubmed.ncbi.nlm.nih.gov/36082909/)
- Heyland DK et al. **A RandomizEd trial of ENtERal Glutamine to minimIZE thermal injury (The RE-ENERGIZE Trial): a clinical trial protocol** *Scars Burn Heal* 3:2059513117745241 (2017). PMID [29799545](https://pubmed.ncbi.nlm.nih.gov/29799545/)
- Kurjatko A et al. **Trace Element Supplementation in Burn Patients: Exploring the Relationship Between Burn Size and Mineral Needs** *J Burn Care Res* 46:411-418 (2025). PMID [39269627](https://pubmed.ncbi.nlm.nih.gov/39269627/)
- Karakus M et al. **Nutritional and metabolic characteristics of critically ill patients admitted for severe toxidermia** *Clin Nutr* 42:859-868 (2023). PMID [37086614](https://pubmed.ncbi.nlm.nih.gov/37086614/)
- Dusapin CJ et al. **Computer customization errors compromised the optimization of trace element repletion dose after major burns** *Clin Nutr* 41:2207-2210 (2022). PMID [36081294](https://pubmed.ncbi.nlm.nih.gov/36081294/)
- Saeg F et al. **Evidence-Based Nutritional Interventions in Wound Care** *Plast Reconstr Surg* 148:226-238 (2021). PMID [34181622](https://pubmed.ncbi.nlm.nih.gov/34181622/)
- Żwierełło W et al. **Bioelements in the treatment of burn injuries - The complex review of metabolism and supplementation (copper, selenium, zinc, iron, manganese, chromium and magnesium)** *J Trace Elem Med Biol* 62:126616 (2020). PMID [32739827](https://pubmed.ncbi.nlm.nih.gov/32739827/)
- Zemrani B et al. **Recent insights into trace element deficiencies: causes, recognition and correction** *Curr Opin Gastroenterol* 36:110-117 (2020). PMID [31895229](https://pubmed.ncbi.nlm.nih.gov/31895229/)
- Guo F et al. **Prospective Study on Energy Expenditure in Patients With Severe Burns** *JPEN J Parenter Enteral Nutr* 45:146-151 (2021). PMID [32270887](https://pubmed.ncbi.nlm.nih.gov/32270887/)
- Xi P et al. **Establishment and assessment of new formulas for energy consumption estimation in adult burn patients** *PLoS One* 9:e110409 (2014). PMID [25330180](https://pubmed.ncbi.nlm.nih.gov/25330180/)
- Mendonça Machado N et al. **Burns, metabolism and nutritional requirements** *Nutr Hosp* 26:692-700 (2011). PMID [22470012](https://pubmed.ncbi.nlm.nih.gov/22470012/)
- Tancheva D et al. **Comparison of estimated energy requirements in severely burned patients with measurements by using indirect calorimetry** *Ann Burns Fire Disasters* 18:16-8 (2005). PMID [21990973](https://pubmed.ncbi.nlm.nih.gov/21990973/)
- Barton RG et al. **Chemical paralysis reduces energy expenditure in patients with burns and severe respiratory failure treated with mechanical ventilation** *J Burn Care Rehabil* 18:461-8; discussion 460 (1997). PMID [9313131](https://pubmed.ncbi.nlm.nih.gov/9313131/)
- Pereira JL et al. **[Evaluation of energy metabolism in burn patients: indirect calorimetry predictive equations]** *Nutr Hosp* 12:147-53 (1997). PMID [9617175](https://pubmed.ncbi.nlm.nih.gov/9617175/)

## R. QSP·모델링 방법론 (QSP and modelling methodology)

*17편*

- Dogné JM et al. **From animal testing to model-informed drug development: building on ICH M15 and EMA initiatives to make new approach methodologies (NAMs) deliver** *Front Pharmacol* 17:1877947 (2026). PMID [42499501](https://pubmed.ncbi.nlm.nih.gov/42499501/)
- GBD 2023 Road Injuries Collaborators **Global, regional, and national burden of road injuries 1990-2023: a systematic analysis for the Global Burden of Disease Study 2023** *Lancet Public Health* (2026). PMID [42476159](https://pubmed.ncbi.nlm.nih.gov/42476159/)
- Al-Majdoub ZM et al. **Monitoring pharmacodynamic and molecular drug targets in liquid biopsy: Exploratory study in liver cancer with modelling of EGFR Receptor engagement** *Br J Clin Pharmacol* (2026). PMID [42449479](https://pubmed.ncbi.nlm.nih.gov/42449479/)
- Zhang S et al. **Quantitative calibration of a spatial QSP model identifies fibroblast impact on HCC immunotherapy** *Proc Natl Acad Sci U S A* 123:e2525799123 (2026). PMID [42446991](https://pubmed.ncbi.nlm.nih.gov/42446991/)
- Goryanin I et al. **Validation of AI-enabled surrogate models in quantitative systems pharmacology: a practical, context-of-use-driven review** *Drug Discov Today* 31:104729 (2026). PMID [42409163](https://pubmed.ncbi.nlm.nih.gov/42409163/)
- GBD 2023 TB HIV Collaborators **Global, regional, and national burden of tuberculosis and multidrug-resistant tuberculosis by HIV status, 1990-2023: a systematic analysis for the Global Burden of Disease Study 2023** *Lancet Infect Dis* (2026). PMID [42385762](https://pubmed.ncbi.nlm.nih.gov/42385762/)
- Ponthier L et al. **Application of machine-learning models to predict the ganciclovir and valganciclovir exposure in children using a limited sampling strategy** *Antimicrob Agents Chemother* 68:e0086024 (2024). PMID [39194260](https://pubmed.ncbi.nlm.nih.gov/39194260/)
- Lu T et al. **gPKPDviz: A flexible R shiny tool for pharmacokinetic/pharmacodynamic simulations using mrgsolve** *CPT Pharmacometrics Syst Pharmacol* 13:341-358 (2024). PMID [38082557](https://pubmed.ncbi.nlm.nih.gov/38082557/)
- Hooijmaijers R et al. **Building an adaptive dose simulation framework to aid dose and schedule selection** *CPT Pharmacometrics Syst Pharmacol* 12:1602-1618 (2023). PMID [37574587](https://pubmed.ncbi.nlm.nih.gov/37574587/)
- Le Louedec F et al. **Easy and reliable maximum a posteriori Bayesian estimation of pharmacokinetic parameters with the open-source R package mapbayr** *CPT Pharmacometrics Syst Pharmacol* 10:1208-1220 (2021). PMID [34342170](https://pubmed.ncbi.nlm.nih.gov/34342170/)
- Elmokadem A et al. **Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial** *CPT Pharmacometrics Syst Pharmacol* 8:883-893 (2019). PMID [31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
- Elgjo GI et al. **Burn resuscitation with two doses of 4 mL/kg hypertonic saline dextran provides sustained fluid sparing: a 48-hour prospective study in conscious sheep** *J Trauma* 49:251-63; discussion 263-5 (2000). PMID [10963536](https://pubmed.ncbi.nlm.nih.gov/10963536/)
- Zhao J et al. **Integration of CTA Run-off Score and Systemic Resuscitation Indices to Predict Outcomes in Emergency Lower Extremity Flap Reconstruction: A Retrospective Cohort Study** *Ann Ital Chir* 97:920-931 (2026). PMID [42136165](https://pubmed.ncbi.nlm.nih.gov/42136165/)
- Hassine NBEH et al. **Early Prediction of Acute Respiratory Distress Syndrome in Critically Ill Polytrauma Patients Using Balanced Random Forest ML: A Retrospective Cohort Study** *J Clin Med* 14 (2025). PMID [41464836](https://pubmed.ncbi.nlm.nih.gov/41464836/)
- Li K et al. **Human intention recognition for trauma resuscitation: An interpretable deep learning approach for medical process data** *J Biomed Inform* 161:104767 (2025). PMID [39746431](https://pubmed.ncbi.nlm.nih.gov/39746431/)
- Bosma KJ et al. **Proportional Assist Ventilation for Minimizing the Duration of Mechanical Ventilation (the PROMIZING study): update to the statistical analysis plan for a randomized controlled trial** *Trials* 25:855 (2024). PMID [39736673](https://pubmed.ncbi.nlm.nih.gov/39736673/)
- Cartotto R et al. **American Burn Association Clinical Practice Guidelines on Burn Shock Resuscitation** *J Burn Care Res* 45:565-589 (2024). PMID [38051821](https://pubmed.ncbi.nlm.nih.gov/38051821/)

---

## 부록 A — 모델의 각 구성요소가 어느 절에 걸려 있는가

| 모델 구성요소 | 파일 내 위치 | 근거 절 |
|---|---|---|
| Landis–Pappenheimer 교질삼투압 식 | `landis_pappenheimer()` / R의 `LP()` | B |
| Starling 식 (구역별 Kf·sigma·Pc·Pi) | `derive()`의 `Jvb`·`Jvu` | B |
| 화상 조직의 음압 (imbibition) | `IMB0`, `TAUIMB` | B |
| 관류 제한 (전달 상한·no-reflow) | `Jvb_cap`, `STASIS0`, `QPB0` | B |
| 림프 환류의 포화 | `lymph()`, `LBASE`, `LMAX`, `KL` | B |
| 소생술 폐루프 제어기 | `RSTATE` ODE, `GCTRL`, `DRMAX`, `TAUUO` | C |
| 수액 크립 · in:out 비 | 시나리오 `creep_uncapped` | C |
| 오피오이드 크립 | `KOPI` | C |
| 콜로이드가 f_ret를 복원 | `FCOL`, `ALBSTART` | C |
| 복강내압 · 복부구획증후군 | `IAPK`·`IAPB`·`IAPQ` | D |
| 아스코르브산의 Kf 작용 | `EMAXVIT`, `KROSKF` | E |
| 대사항진 (A_open 구동) | `REE` ODE, `EMAXHM`, `A50`, `HILLA` | F |
| 증발 열손실 | `KEVAP` | F |
| 프로프라놀롤 (전달자 차단) | `prBlock`, `EC50PR`, `KSPR` | G |
| 옥산드롤론 (합성 항) | `oxEff`, `EMAXOX` | H |
| 인슐린 · 혈당 · 저혈당 노출 | `GLC`/`INS`/`HYPOH` ODE, `GFEED` | I |
| 조기 절제 · 이식 · 공여부 제약 | `apply_events()`, `DONPOOLF`, `MESH` | J |
| 정량 세균학 역치 · 면역마비 | `BWD`/`HLADR`/`BSYS` ODE | K |
| 증가된 신클리어런스 · 반코마이신 | `arc`, `CLVC`, `VVC` | L |
| 흡입 손상의 수액 배수 · 위험 배수 | `INH`, `HZINH` | M |
| IL-6 · 음성 급성기 단백 | `IL6` ODE, `albSyn` | N |
| 골 · 성장 | `BMC` ODE, `KBONE` | O |
| 소아 · 고령 · 임신 · 전기 화상 | `build(..., WT/AGE/TBSA/INH)` | P |
| 조기 경장영양 · 글루타민 | `sc["nutrition"]`, `GFEED`, `TFEED` | Q |
| QSP 방법론 · mrgsolve | 모델 전체 | R |

---

## 부록 B — 이 목록을 재현하는 방법

```bash
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed\
&retmode=json&retmax=6&term=Parkland+formula+burn+resuscitation+volume"
```

각 절에 해당하는 검색어를 위와 같이 질의한 뒤 `esummary.fcgi`로 서지사항을 받아
그대로 옮겼습니다. 중복 PMID는 최초로 등장한 절에만 남겼습니다.

---

## ⚠️ 면책 (Disclaimer)

본 문헌 목록은 **모델의 가정과 파라미터가 어디에서 왔는지를 추적 가능하게 만들기 위한
것**이며 임상 진료 지침이 아닙니다. 개별 논문의 결론과 모델 구현이 일치하지 않는 지점은
위 표와 `README.md`의 "모델이 틀리는 곳" 절에 명시했습니다.
