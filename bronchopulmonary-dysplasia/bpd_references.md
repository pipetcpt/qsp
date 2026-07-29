# Bronchopulmonary Dysplasia (BPD) — QSP Model References

**Verification status.** Every PMID below was resolved against the PubMed
E-utilities `esummary` endpoint during the session that built this model, and
the author / year / journal / title shown is the record PubMed returned. Where
a citation I intended could not be resolved to a matching record, it was
**removed** rather than guessed. Nothing in this list is reconstructed from
memory. Interpretations and the mapping onto model parameters are mine, not the
authors'; a reference appearing here means the paper is the evidence I used for
a structure or a number, not that its authors would endorse this model.

Each section ends with a short **→ model** note saying exactly which parameter,
equation or scenario in `bpd_mrgsolve_model.R` the section supports.

---

## 1. Definition, grading, epidemiology, and why the endpoint is a MODE of support

| # | Reference | PMID |
|---|-----------|------|
| 1 | Northway WH Jr, Rosan RC, Porter DY. Pulmonary disease following respirator therapy of hyaline-membrane disease. Bronchopulmonary dysplasia. *N Engl J Med* 1967. | [5334613](https://pubmed.ncbi.nlm.nih.gov/5334613/) |
| 2 | Jobe AH, Bancalari E. Bronchopulmonary dysplasia. *Am J Respir Crit Care Med* 2001. (the 2001 NICHD severity definition) | [11401896](https://pubmed.ncbi.nlm.nih.gov/11401896/) |
| 3 | Jobe AH. The new BPD: an arrest of lung development. *Pediatr Res* 1999. | [10590017](https://pubmed.ncbi.nlm.nih.gov/10590017/) |
| 4 | Bancalari E, Claure N, Sosenko IR. Bronchopulmonary dysplasia: changes in pathogenesis, epidemiology and definition. *Semin Neonatol* 2003. | [12667831](https://pubmed.ncbi.nlm.nih.gov/12667831/) |
| 5 | Coalson JJ. Pathology of new bronchopulmonary dysplasia. *Semin Neonatol* 2003. | [12667832](https://pubmed.ncbi.nlm.nih.gov/12667832/) |
| 6 | Jensen EA, Dysart K, Gantz MG, et al. The diagnosis of bronchopulmonary dysplasia in very preterm infants: an evidence-based approach. *Am J Respir Crit Care Med* 2019. (the 2019 grade 1/2/3 definition used as this model's endpoint) | [30995069](https://pubmed.ncbi.nlm.nih.gov/30995069/) |
| 7 | Higgins RD, Jobe AH, Koso-Thomas M, et al. Bronchopulmonary dysplasia: executive summary of a workshop. *J Pediatr* 2018. | [29551318](https://pubmed.ncbi.nlm.nih.gov/29551318/) |
| 8 | Isayama T, Lee SK, Yang J, et al. Revisiting the definition of bronchopulmonary dysplasia: effect of changing panoply of respiratory support. *JAMA Pediatr* 2017. | [28114678](https://pubmed.ncbi.nlm.nih.gov/28114678/) |
| 9 | Thébaud B, Goss KN, Laughon M, et al. Bronchopulmonary dysplasia. *Nat Rev Dis Primers* 2019. | [31727986](https://pubmed.ncbi.nlm.nih.gov/31727986/) |
| 10 | Abman SH, Bancalari E, Jobe A. The evolution of bronchopulmonary dysplasia after 50 years. *Am J Respir Crit Care Med* 2017. | [28199157](https://pubmed.ncbi.nlm.nih.gov/28199157/) |
| 11 | Bancalari E, Jain D. Bronchopulmonary dysplasia: 50 years after the original description. *Neonatology* 2019. | [30974430](https://pubmed.ncbi.nlm.nih.gov/30974430/) |
| 12 | Kalikkot Thekkeveedu R, Guaman MC, Shivanna B. Bronchopulmonary dysplasia: a review of pathogenesis and pathophysiology. *Respir Med* 2017. | [29229093](https://pubmed.ncbi.nlm.nih.gov/29229093/) |
| 13 | Stoll BJ, Hansen NI, Bell EF, et al. Trends in care practices, morbidity, and mortality of extremely preterm neonates, 1993-2012. *JAMA* 2015. | [26348753](https://pubmed.ncbi.nlm.nih.gov/26348753/) |
| 14 | Stoll BJ, Hansen NI, Bell EF, et al. Neonatal outcomes of extremely preterm infants from the NICHD Neonatal Research Network. *Pediatrics* 2010. | [20732945](https://pubmed.ncbi.nlm.nih.gov/20732945/) |
| 15 | Bell EF, Hintz SR, Hansen NI, et al. Mortality, in-hospital morbidity, care practices, and 2-year outcomes for extremely preterm infants in the US, 2013-2018. *JAMA* 2022. | [35040888](https://pubmed.ncbi.nlm.nih.gov/35040888/) |
| 16 | Laughon MM, Langer JC, Bose CL, et al. Prediction of bronchopulmonary dysplasia by postnatal age in extremely premature infants. *Am J Respir Crit Care Med* 2011. | [21471086](https://pubmed.ncbi.nlm.nih.gov/21471086/) |
| 17 | Álvarez-Fuente M, Arruza L, Muro M, et al. The economic impact of prematurity and bronchopulmonary dysplasia. *Eur J Pediatr* 2017. | [28889192](https://pubmed.ncbi.nlm.nih.gov/28889192/) |
| 18 | Kair LR, Leonard DT, Anderson JM. Bronchopulmonary dysplasia. *Pediatr Rev* 2012. | [22659256](https://pubmed.ncbi.nlm.nih.gov/22659256/) |

**→ model.** These fix the endpoint. BPD is scored in `$TABLE` as a mode of
support at 36 weeks PMA (`GRADE` 0/1/2/3 from `MVF`, `SUP` and `FIO2`), not as
a threshold on any biomarker, because that is what the definition actually is
(#2, #6). Refs #13–#15 set the survival and morbidity magnitudes that
`H0`/`H_HYPOX`/`H_PVR`/`H_G3` are scaled against. Ref #8 is the reason the model
scores the mode of support rather than an FiO2 cut-off: as support technology
changed, so did the measured incidence, without the biology changing at all.
Ref #3 is the sentence this entire model is built around.

---

## 2. Lung development and the finite window — the model's central structure

| # | Reference | PMID |
|---|-----------|------|
| 19 | Schittny JC. Development of the lung. *Cell Tissue Res* 2017. | [28144783](https://pubmed.ncbi.nlm.nih.gov/28144783/) |
| 20 | Herring MJ, Putney LF, Wyatt G, et al. Growth of alveoli during postnatal development in humans based on stereological estimation. *Am J Physiol Lung Cell Mol Physiol* 2014. | [24907055](https://pubmed.ncbi.nlm.nih.gov/24907055/) |
| 21 | Narayanan M, Owers-Bradley J, Beardsmore CS, et al. Alveolarization continues during childhood and adolescence: new evidence from helium-3 magnetic resonance. *Am J Respir Crit Care Med* 2012. | [22071328](https://pubmed.ncbi.nlm.nih.gov/22071328/) |
| 22 | Thurlbeck WM. Postnatal human lung growth. *Thorax* 1982. | [7179184](https://pubmed.ncbi.nlm.nih.gov/7179184/) |
| 23 | Schmid L, Tschanz SA, Yamada M, et al. Microvascular maturation of the septal capillary layers takes place in parallel to alveolarization. *Am J Physiol Lung Cell Mol Physiol* 2023. | [37605833](https://pubmed.ncbi.nlm.nih.gov/37605833/) |
| 24 | Bourbon JR, Boucherat O, Boczkowski J, et al. Bronchopulmonary dysplasia and emphysema: in search of common therapeutic targets. *Trends Mol Med* 2009. | [19303361](https://pubmed.ncbi.nlm.nih.gov/19303361/) |
| 25 | Massaro GD, Massaro D. Postnatal treatment with retinoic acid increases the number of pulmonary alveoli in rats. *Am J Physiol* 1996. | [8780001](https://pubmed.ncbi.nlm.nih.gov/8780001/) |
| 26 | Blanco LN, Massaro GD, Massaro D. On increasing the number of pulmonary alveoli in rats by postnatal treatment with retinoic acid. *Am J Physiol* 1997. | [9142956](https://pubmed.ncbi.nlm.nih.gov/9142956/) |
| 27 | Lévy M, Maurey C, Dinh-Xuan AT, et al. Developmental expression of vasoactive and growth factors in human lung: role in pulmonary vascular resistance adaptation at birth. *Pediatr Res* 2005. | [15817500](https://pubmed.ncbi.nlm.nih.gov/15817500/) |
| 28 | Hilgendorff A, Reiss I, Ehrhardt H, et al. Chronic lung disease in the preterm infant: lessons learned from animal models. *Am J Respir Cell Mol Biol* 2014. | [24024524](https://pubmed.ncbi.nlm.nih.gov/24024524/) |

**→ model.** `S_dev(PMA)` — the asymmetric Gaussian septation-drive window with
`PMAPK`, `SDW_L`, `SDW_R` — is a phenomenological summary of #19–#23: septation
and microvascular maturation proceed together (#23, hence `CAPD` shares the same
`S_dev` and `r_prog`), the process is largely postnatal (#20, #22), and it has a
long tail into childhood and adolescence rather than a hard stop at term (#21,
which is why the model can be run past 36 weeks and why `LOSTW` is only
partially recoverable rather than not at all). The retinoic-acid papers (#25,
#26) are why `RASIG` is a first-class member of the growth gate `G` and why the
gate has a `WR` weight at all — retinoate does not merely permit septation, it
increases alveolar number.

**Honest limitation.** The exact onset of human secondary septation is
contested. `SDW_L = 6` weeks (giving `S_dev` ≈ 0.14 at PMA 20 and ≈ 0.61 at
PMA 24) is a modelling choice that places substantial drive in the late
saccular stage. A narrower rising limb makes the gestational gradient steeper
and a wider one flattens it; `BPD_window_remaining()` in the R file and
experiment A1 in `bpd_reference_impl.py` both expose this so the assumption can
be varied rather than trusted.

---

## 3. Angiogenesis, VEGF, NO and IGF-1 — the components of the growth gate G

| # | Reference | PMID |
|---|-----------|------|
| 29 | Bhatt AJ, Pryhuber GS, Huyck H, et al. Disrupted pulmonary vasculature and decreased VEGF, Flt-1 and TIE-2 in human infants dying with bronchopulmonary dysplasia. *Am J Respir Crit Care Med* 2001. | [11734454](https://pubmed.ncbi.nlm.nih.gov/11734454/) |
| 30 | Thébaud B, Abman SH. Bronchopulmonary dysplasia: where have all the vessels gone? Roles of angiogenic growth factors in chronic lung disease. *Am J Respir Crit Care Med* 2007. | [17272782](https://pubmed.ncbi.nlm.nih.gov/17272782/) |
| 31 | Jakkula M, Le Cras TD, Gebb S, et al. Inhibition of angiogenesis decreases alveolarization in the developing rat lung. *Am J Physiol Lung Cell Mol Physiol* 2000. | [10956636](https://pubmed.ncbi.nlm.nih.gov/10956636/) |
| 32 | Alapati D, Rong M, Chen S, et al. Inhibition of β-catenin signaling improves alveolarization and reduces pulmonary hypertension in experimental bronchopulmonary dysplasia. *Am J Respir Cell Mol Biol* 2014. | [24484510](https://pubmed.ncbi.nlm.nih.gov/24484510/) |
| 33 | Benjamin JT, Carver BJ, Plosa EJ, et al. NF-κB activation limits airway branching through inhibition of Sp1-mediated fibroblast growth factor-10 expression. *J Immunol* 2010. | [20861353](https://pubmed.ncbi.nlm.nih.gov/20861353/) |
| 34 | Nakanishi H, Sugiura T, Streisand JB, et al. TGF-β-neutralizing antibodies improve pulmonary alveologenesis and vasculogenesis in the injured newborn lung. *Am J Physiol Lung Cell Mol Physiol* 2007. | [17400601](https://pubmed.ncbi.nlm.nih.gov/17400601/) |
| 35 | Bland RD, Ertsey R, Mokres LM, et al. Mechanical ventilation uncouples synthesis and assembly of elastin and increases apoptosis in lungs of newborn mice. *Am J Physiol Lung Cell Mol Physiol* 2008. | [17934062](https://pubmed.ncbi.nlm.nih.gov/17934062/) |
| 36 | Thibeault DW, Mabry SM, Ekekezie II, Truog WE. Lung elastic tissue maturation and perturbations during the evolution of chronic lung disease. *Pediatrics* 2000. | [11099603](https://pubmed.ncbi.nlm.nih.gov/11099603/) |
| 37 | Löfqvist C, Hellgren G, Niklasson A, et al. Low postnatal serum IGF-I levels are associated with bronchopulmonary dysplasia. *Acta Paediatr* 2012. | [22924869](https://pubmed.ncbi.nlm.nih.gov/22924869/) |
| 38 | Hellström A, Ley D, Hansen-Pupp I, et al. Role of insulin-like growth factor 1 in fetal development and in the early postnatal life of premature infants. *Am J Perinatol* 2016. | [27603537](https://pubmed.ncbi.nlm.nih.gov/27603537/) |
| 39 | Rozance PJ, Seedorf GJ, Brown A, et al. Intrauterine growth restriction decreases pulmonary alveolar and vessel growth and causes pulmonary artery endothelial cell dysfunction in vitro in fetal sheep. *Am J Physiol Lung Cell Mol Physiol* 2011. | [21873446](https://pubmed.ncbi.nlm.nih.gov/21873446/) |

**→ model.** Ref #31 is the load-bearing citation for the whole architecture:
blocking angiogenesis alone is sufficient to reduce alveolarisation, which is
why `G` multiplies the septation term rather than being one additive injury
among many, and why `Capillary_sprout → Secondary_septa` is drawn as a
permissive edge in the map. #29 and #30 set the direction and rough magnitude of
`AV_IL1`, `AV_ROS`, `AV_SFLT`. #34 justifies `AG_TGF` acting on the gate rather
than only on `SEPT`. #35 and #36 are the basis for the model's distinction
between elastin *amount* and elastin *location* (`Elastin_mislocal` in the map)
and for `KDEGALV`. #37–#39 set `IGFB0`, `IGFBSL`, `AI_INFL` and the `IUGR`
multiplier — and #39 is why `IUGR` in this model lowers **both** `IGF1` and
`VEGF` (via `sflt`) at birth rather than acting as a generic risk flag.

---

## 4. Inflammation as the transduction stage, and IL-1β as the specific blocker

| # | Reference | PMID |
|---|-----------|------|
| 40 | Bry K, Whitsett JA, Lappalainen U. IL-1β disrupts postnatal lung morphogenesis in the mouse. *Am J Respir Cell Mol Biol* 2007. | [16888287](https://pubmed.ncbi.nlm.nih.gov/16888287/) |
| 41 | Nold MF, Mangan NE, Rudloff I, et al. Interleukin-1 receptor antagonist prevents murine bronchopulmonary dysplasia induced by perinatal inflammation and hyperoxia. *Proc Natl Acad Sci U S A* 2013. | [23946428](https://pubmed.ncbi.nlm.nih.gov/23946428/) |
| 42 | Toth A, Steinmeyer S, Kannan P, et al. Inflammatory blockade prevents injury to the developing pulmonary gas exchange surface in preterm primates. *Sci Transl Med* 2022. | [35353543](https://pubmed.ncbi.nlm.nih.gov/35353543/) |
| 43 | D'Angio CT, Ambalavanan N, Carlo WA, et al. Blood cytokine profiles associated with distinct patterns of bronchopulmonary dysplasia among extremely low birth weight infants. *J Pediatr* 2016. | [27117196](https://pubmed.ncbi.nlm.nih.gov/27117196/) |
| 44 | Watterberg KL, Demers LM, Scott SM, Murphy S. Chorioamnionitis and early lung inflammation in infants in whom bronchopulmonary dysplasia develops. *Pediatrics* 1996. | [8584379](https://pubmed.ncbi.nlm.nih.gov/8584379/) |
| 45 | Been JV, Rours IG, Kornelisse RF, et al. Histological chorioamnionitis and respiratory outcome in preterm infants. *Arch Dis Child Fetal Neonatal Ed* 2009. | [19131431](https://pubmed.ncbi.nlm.nih.gov/19131431/) |
| 46 | Been JV, Rours IG, Kornelisse RF, et al. Chorioamnionitis alters the response to surfactant in preterm infants. *J Pediatr* 2010. | [19833352](https://pubmed.ncbi.nlm.nih.gov/19833352/) |
| 47 | Van Marter LJ, Dammann O, Allred EN, et al. Chorioamnionitis, mechanical ventilation, and postnatal sepsis as modulators of chronic lung disease in preterm infants. *J Pediatr* 2002. | [11865267](https://pubmed.ncbi.nlm.nih.gov/11865267/) |
| 48 | Bhattacharya S, Go D, Krenitsky DL, et al. Genome-wide transcriptional profiling reveals connective tissue mast cell accumulation in bronchopulmonary dysplasia. *Am J Respir Crit Care Med* 2012. | [22723293](https://pubmed.ncbi.nlm.nih.gov/22723293/) |

**→ model.** #40 is why `IL1B` has a **direct** septation-blocking term
(`AG_IL1`) in addition to its effect through `VEGF`: IL-1β overexpression alone
arrests lung morphogenesis. #41 and #42 are the basis for the `IL1RA` switch —
and #42, in preterm primates, is the closest thing to a validation that
blocking the transduction stage protects the gas-exchange surface itself.
#44–#47 set `KANTE`, `TAU_ANTE` and the `CHORIO` initial conditions, and #47 in
particular is why the model treats chorioamnionitis and mechanical ventilation
as *interacting* rather than additive (both feed one `INJ` term that then
saturates). #43 is the reason the model exposes distinct phenotypes from one
parameter set rather than one average trajectory.

---

## 5. Oxidative stress and the immature antioxidant defence

| # | Reference | PMID |
|---|-----------|------|
| 49 | Saugstad OD. Bronchopulmonary dysplasia — oxidative stress and antioxidants. *Semin Neonatol* 2003. | [12667829](https://pubmed.ncbi.nlm.nih.gov/12667829/) |
| 50 | Frank L, Sosenko IR. Development of lung antioxidant enzyme system in late gestation: possible implications for the prematurely born infant. *J Pediatr* 1987. | [3540251](https://pubmed.ncbi.nlm.nih.gov/3540251/) |

**→ model.** #50 is the direct source of `AOX0`, `AOXPMA = 31` weeks and
`AOXSL`: the antioxidant enzyme system matures in *late* gestation, so the
infant who most needs it has least of it. In the model this appears as the
`1/(1 + aoxm)` divisor on `oxtox` and on `ROS` clearance, which is why the same
FiO2 is more damaging at PMA 26 than at PMA 34 — a PMA-dependence that would
otherwise have to be asserted.

---

## 6. Mechanical injury and ventilation strategy — the exposure term

| # | Reference | PMID |
|---|-----------|------|
| 51 | Hillman NH, Moss TJ, Kallapur SG, et al. Brief, large tidal volume ventilation initiates lung injury and a systemic response in fetal sheep. *Am J Respir Crit Care Med* 2007. | [17641159](https://pubmed.ncbi.nlm.nih.gov/17641159/) |
| 52 | Björklund LJ, Ingimarsson J, Curstedt T, et al. Manual ventilation with a few large breaths at birth compromises the therapeutic effect of subsequent surfactant replacement in immature lambs. *Pediatr Res* 1997. | [9284276](https://pubmed.ncbi.nlm.nih.gov/9284276/) |
| 53 | Brew N, Hooper SB, Allison BJ, et al. Injury and repair in the very immature lung following brief mechanical ventilation. *Am J Physiol Lung Cell Mol Physiol* 2011. | [21890511](https://pubmed.ncbi.nlm.nih.gov/21890511/) |
| 54 | Klingenberg C, Wheeler KI, McCallion N, et al. Volume-targeted versus pressure-limited ventilation in neonates. *Cochrane Database Syst Rev* 2017. | [29039883](https://pubmed.ncbi.nlm.nih.gov/29039883/) |
| 55 | Schmölzer GM, Kumar M, Pichler G, et al. Non-invasive versus invasive respiratory support in preterm infants at birth: systematic review and meta-analysis. *BMJ* 2013. | [24136633](https://pubmed.ncbi.nlm.nih.gov/24136633/) |
| 56 | Abdel-Latif ME, Davis PG, Wheeler KI, et al. Surfactant therapy via thin catheter in preterm infants with or at risk of respiratory distress syndrome. *Cochrane Database Syst Rev* 2021. | [33970483](https://pubmed.ncbi.nlm.nih.gov/33970483/) |
| 57 | Kirpalani H, Ratcliffe SJ, Keszler M, et al. Effect of sustained inflations vs intermittent positive pressure ventilation on bronchopulmonary dysplasia or death among extremely preterm infants (SAIL). *JAMA* 2019. | [30912836](https://pubmed.ncbi.nlm.nih.gov/30912836/) |
| 58 | Silveira RC, Panceri C, Munõz NP, et al. Less invasive surfactant administration versus intubation-surfactant-extubation in preterm neonates. *J Pediatr (Rio J)* 2024. | [37353207](https://pubmed.ncbi.nlm.nih.gov/37353207/) |

**→ model.** `KVILI`, `NVILI`, `KNIV`, `EFF_VTV`, `EFF_PHC` and the `LISA`
switch come from this block. #51–#53 justify a *superlinear* exposure term
(`NVILI > 1`) and a very early onset: injury is initiated within minutes to
hours, not weeks, which in the model means `vili` contributes to `INJ` from the
delivery room onward. #54 sets `EFF_VTV ≈ 0.35` and is also why `HYPOCAP` is
gated off by `VTV` in the NDI accumulator — the meta-analysis found less
hypocarbia as well as less BPD. #55, #56 and #58 set the LISA/nCPAP arm. #57 is
in the list as a **negative** trial and is represented in the map by a dashed
edge: a plausible delivery-room manoeuvre that did not improve the endpoint.

---

## 7. Oxygen saturation targeting — the U-shaped exposure

| # | Reference | PMID |
|---|-----------|------|
| 59 | Askie LM, Darlow BA, Finer N, et al. Association between oxygen saturation targeting and death or disability in extremely preterm infants in the NeOProM Collaboration. *JAMA* 2018. | [29872859](https://pubmed.ncbi.nlm.nih.gov/29872859/) |
| 60 | Carlo WA, Finer NN, Walsh MC, et al. (SUPPORT Study Group). Target ranges of oxygen saturation in extremely preterm infants. *N Engl J Med* 2010. | [20472937](https://pubmed.ncbi.nlm.nih.gov/20472937/) |
| 61 | Stenson BJ, Tarnow-Mordi WO, Darlow BA, et al. Oxygen saturation and outcomes in preterm infants (BOOST-II). *N Engl J Med* 2013. | [23642047](https://pubmed.ncbi.nlm.nih.gov/23642047/) |
| 62 | Poets CF, Roberts RS, Schmidt B, et al. Association between intermittent hypoxemia or bradycardia and late death or disability in extremely preterm infants. *JAMA* 2015. | [26262797](https://pubmed.ncbi.nlm.nih.gov/26262797/) |
| 63 | Di Fiore JM, Bloom JN, Orge F, et al. A higher incidence of intermittent hypoxemic episodes is associated with severe retinopathy of prematurity. *J Pediatr* 2010. | [20304417](https://pubmed.ncbi.nlm.nih.gov/20304417/) |

**→ model.** These are the reason the model refuses to produce a single optimal
SpO2 target. The low band lowers `FiO2` for a given lung (`BANDF`) and therefore
lowers `oxtox`, while raising `hb` and hence `H_HYPOX`, `KHYPOX` and the HPV
term. Experiment A8 sweeps the band and reports the lung-optimal and
survival-optimal bands separately, because that is the actual finding (#59):
lower targets reduced retinopathy, higher targets improved survival, and the
composite was a wash. #62 and #63 are why intermittent hypoxaemia enters the NDI
and ROP paths and not only the pulmonary-vascular path.

---

## 8. Caffeine — pharmacokinetics and the CAP trial

| # | Reference | PMID |
|---|-----------|------|
| 64 | Schmidt B, Roberts RS, Davis P, et al. Caffeine therapy for apnea of prematurity (CAP). *N Engl J Med* 2006. | [16707748](https://pubmed.ncbi.nlm.nih.gov/16707748/) |
| 65 | Schmidt B, Roberts RS, Davis P, et al. Long-term effects of caffeine therapy for apnea of prematurity. *N Engl J Med* 2007. | [17989382](https://pubmed.ncbi.nlm.nih.gov/17989382/) |
| 66 | Aranda JV, Cook CE, Gorman W, et al. Pharmacokinetic profile of caffeine in the premature newborn infant with apnea. *J Pediatr* 1979. | [430317](https://pubmed.ncbi.nlm.nih.gov/430317/) |
| 67 | Lee TC, Charles B, Steer P, et al. Population pharmacokinetics of intravenous caffeine in neonates with apnea of prematurity. *Clin Pharmacol Ther* 1997. | [9209245](https://pubmed.ncbi.nlm.nih.gov/9209245/) |
| 68 | Charles BG, Townsend SR, Steer PA, et al. Caffeine citrate treatment for extremely premature infants with apnea: population pharmacokinetics. *Ther Drug Monit* 2008. | [19057373](https://pubmed.ncbi.nlm.nih.gov/19057373/) |
| 69 | Bruschettini M, Brattström P, Russo C, et al. Caffeine dosing regimens in preterm infants with or at risk for apnea of prematurity. *Cochrane Database Syst Rev* 2023. | [37040532](https://pubmed.ncbi.nlm.nih.gov/37040532/) |

**→ model.** `CAF_VD = 0.85 L/kg`, `CAF_CLI = 4.6 L/day/kg`,
`CAF_PMA50 = 60 weeks`, `CAF_HILL = 4` were chosen so that the maturation curve
passes through the published anchors from #66–#68: a half-life of roughly 90–100
hours in the very preterm neonate falling towards the adult 5–6 hours by
term-corrected age, and a trough of about 10–20 mg/L on 5 mg/kg/day of caffeine
citrate. Experiment A0 prints the resulting table so the fit can be checked
rather than believed. #64 supplies the effect size the model is asked to
reproduce (BPD 36.3% vs 47.2%), and #69 is the reason `CAF_MD` is a free
parameter — the optimal dose is still open.

**The model makes a claim here.** Caffeine is represented with two separable
arms: a central-drive arm that shortens ventilation (`CAF_DRIVE`, acting through
`apnoea → mvdrive → MVF → vili`) and a small direct adenosine-A2A arm on the
gate (`CAF_EALV`). Experiment A3 turns the first off and reports how much of the
36-week benefit survives. The model attributes most of it to exposure reduction.
That is a mechanistic hypothesis, not a finding of CAP, and it is falsifiable:
a caffeine analogue with no respiratory-drive effect should, on this model,
lose most of the BPD benefit.

---

## 9. Corticosteroids — one receptor, several schedules, opposite conclusions

| # | Reference | PMID |
|---|-----------|------|
| 70 | Baud O, Maury L, Lebail F, et al. Effect of early low-dose hydrocortisone on survival without bronchopulmonary dysplasia in extremely preterm infants (PREMILOC). *Lancet* 2016. | [26916176](https://pubmed.ncbi.nlm.nih.gov/26916176/) |
| 71 | Baud O, Trousson C, Biran V, et al. Association between early low-dose hydrocortisone therapy in extremely preterm neonates and neurodevelopmental outcomes at 2 years of age. *JAMA* 2017. | [28384828](https://pubmed.ncbi.nlm.nih.gov/28384828/) |
| 72 | Doyle LW, Davis PG, Morley CJ, et al. Low-dose dexamethasone facilitates extubation among chronically ventilator-dependent infants (DART). *Pediatrics* 2006. | [16396863](https://pubmed.ncbi.nlm.nih.gov/16396863/) |
| 73 | Doyle LW, Cheong JL, Hay S, et al. Early (<7 days) systemic postnatal corticosteroids for prevention of bronchopulmonary dysplasia. *Cochrane Database Syst Rev* 2021. | [34674229](https://pubmed.ncbi.nlm.nih.gov/34674229/) |
| 74 | Doyle LW, Cheong JL, Hay S, et al. Late (≥7 days) systemic postnatal corticosteroids for prevention of bronchopulmonary dysplasia. *Cochrane Database Syst Rev* 2021. | [34758507](https://pubmed.ncbi.nlm.nih.gov/34758507/) |
| 75 | Doyle LW, Halliday HL, Ehrenkranz RA, et al. An update on the impact of postnatal systemic corticosteroids on mortality and cerebral palsy in preterm infants: effect modification by risk of bronchopulmonary dysplasia. *J Pediatr* 2014. | [25217197](https://pubmed.ncbi.nlm.nih.gov/25217197/) |
| 76 | Barrington KJ. The adverse neuro-developmental effects of postnatal steroids in the preterm infant: a systematic review of RCTs. *BMC Pediatr* 2001. | [11248841](https://pubmed.ncbi.nlm.nih.gov/11248841/) |
| 77 | Yeh TF, Lin YJ, Lin HC, et al. Outcomes at school age after postnatal dexamethasone therapy for lung disease of prematurity. *N Engl J Med* 2004. | [15044641](https://pubmed.ncbi.nlm.nih.gov/15044641/) |
| 78 | Powell K, Kerkering KW, Barker G, Rozycki HJ. Dexamethasone dosing, mechanical ventilation and the risk of cerebral palsy. *J Matern Fetal Neonatal Med* 2006. | [16492591](https://pubmed.ncbi.nlm.nih.gov/16492591/) |
| 79 | Watterberg KL, Gerdes JS, Cole CH, et al. Prophylaxis of early adrenal insufficiency to prevent bronchopulmonary dysplasia: a multicenter trial. *Pediatrics* 2004. | [15574629](https://pubmed.ncbi.nlm.nih.gov/15574629/) |
| 80 | Watterberg KL, Walsh MC, Li L, et al. Hydrocortisone to improve survival without bronchopulmonary dysplasia. *N Engl J Med* 2022. | [35320643](https://pubmed.ncbi.nlm.nih.gov/35320643/) |
| 81 | Onland W, Cools F, Kroon A, et al. Effect of hydrocortisone therapy initiated 7 to 14 days after birth on mortality or bronchopulmonary dysplasia (SToP-BPD). *JAMA* 2019. | [30694322](https://pubmed.ncbi.nlm.nih.gov/30694322/) |
| 82 | Yeh TF, Chen CM, Wu SY, et al. Intratracheal administration of budesonide/surfactant to prevent bronchopulmonary dysplasia. *Am J Respir Crit Care Med* 2016. | [26351971](https://pubmed.ncbi.nlm.nih.gov/26351971/) |
| 83 | Bassler D, Plavka R, Shinwell ES, et al. Early inhaled budesonide for the prevention of bronchopulmonary dysplasia (NEUROSIS). *N Engl J Med* 2015. | [26465983](https://pubmed.ncbi.nlm.nih.gov/26465983/) |
| 84 | Onland W, Offringa M, van Kaam A. Late (≥7 days) inhaled corticosteroids to reduce bronchopulmonary dysplasia. *Cochrane Database Syst Rev* 2022. | [36521169](https://pubmed.ncbi.nlm.nih.gov/36521169/) |
| 85 | Hay S, Ovelman C, Zupancic JA, et al. Systemic corticosteroids for the prevention of bronchopulmonary dysplasia, a network meta-analysis. *Cochrane Database Syst Rev* 2023. | [37650547](https://pubmed.ncbi.nlm.nih.gov/37650547/) |
| 86 | Lok IM, Bourgonje LM, Nijman J, et al. Effects of postnatal corticosteroids on lung development in newborn animals: a systematic review. *Pediatr Res* 2024. | [38493255](https://pubmed.ncbi.nlm.nih.gov/38493255/) |

**→ model.** This is the block the model was built to explain. All three
ligands feed **one** occupancy equation (`GR_EC50`, `GRocc_l`, `GRocc_s`) with
`HC_POT = 0.04` and `BUD_POT = 2.0`; the DART taper (#72) and the PREMILOC
schedule (#70) are hard-coded so the exposure comparison is not fudged. Three
distinct Emax terms come off that one occupancy:

- `GR_EMAXR` — transrepression of NF-κB, the benefit;
- `GR_ECOMP` — the compliance/extubation effect, which is what DART measured;
- `AG_GR` — the **antiproliferative effect on septation itself**, the cost.

`AG_GR` exists because of #86: postnatal corticosteroids impair alveolarisation
in animal models. Without that term a model would predict that more steroid,
earlier, is always better — the prediction the trial record refutes. The NDI
accumulator's `W_DEX` and `W_EARLY` come from #75–#78, and the near-neutral
2-year outcome of low-dose hydrocortisone (#71) is why `HC_POT` scales the NDI
input rather than hydrocortisone carrying a separate penalty.

**Where the model and the literature disagree, stated plainly.** #80
(a large hydrocortisone trial started at 14–28 days) did **not** show improved
survival without BPD, and #81 (SToP-BPD, days 7–14) did not reach its primary
endpoint, whereas #70 (day 1) did. The model reproduces that ordering by timing
alone, which is encouraging. But #85, the network meta-analysis, is less kind to
a clean timing story than this model is: it finds moderate-to-late dexamethasone
among the more effective regimens for the BPD endpoint. A reader should treat
the model's timing claim as a mechanistic hypothesis that the trial record is
*consistent with* rather than one it has established.

---

## 10. Vitamin A

| # | Reference | PMID |
|---|-----------|------|
| 87 | Tyson JE, Wright LL, Oh W, et al. Vitamin A supplementation for extremely-low-birth-weight infants. *N Engl J Med* 1999. | [10379020](https://pubmed.ncbi.nlm.nih.gov/10379020/) |

**→ model.** The `VITA_*` block implements the Tyson regimen exactly (5000 IU
IM three times weekly for four weeks, 12 doses), routed through a hepatic store
(`RETS`) to plasma retinol (`RETP`) to `RASIG`. The IU-to-µmol conversion is
1 IU retinol = 0.3 µg = 0.001047 µmol. `RET_S0 = 3 µmol` encodes the fact that
the preterm infant is born retinol-deficient, which is why supplementation
works at all and why the effect size is modest — it repletes a deficiency
rather than supraphysiologically driving `RAR`.

---

## 11. Ureaplasma and azithromycin — the subgroup drug

| # | Reference | PMID |
|---|-----------|------|
| 88 | Lowe J, Watkins WJ, Edwards MO, et al. Association between pulmonary Ureaplasma colonization and bronchopulmonary dysplasia in preterm infants: updated systematic review and meta-analysis. *Pediatr Infect Dis J* 2014. | [24445836](https://pubmed.ncbi.nlm.nih.gov/24445836/) |
| 89 | Ballard HO, Shook LA, Bernard P, et al. Use of azithromycin for the prevention of bronchopulmonary dysplasia in preterm infants: a randomized, double-blind, placebo-controlled trial. *Pediatr Pulmonol* 2011. | [20963840](https://pubmed.ncbi.nlm.nih.gov/20963840/) |

**→ model.** `UREA` is a logistic-growth state killed by `E_azi`, so
azithromycin can only act where there is something to kill. Experiment A6 runs
the colonised and non-colonised cases separately and shows that an unselected
trial dilutes a real effect with a null one — which is a trial-design
conclusion generated by mechanism, and the main reason the drug is in the model.
`AZI_AI` carries the non-antimicrobial anti-inflammatory effect separately so
the two mechanisms are not conflated.

---

## 12. Inhaled NO, PDE5 inhibition and BPD-associated pulmonary hypertension

| # | Reference | PMID |
|---|-----------|------|
| 90 | Ballard RA, Truog WE, Cnaan A, et al. Inhaled nitric oxide in preterm infants undergoing mechanical ventilation (NO CLD). *N Engl J Med* 2006. | [16870913](https://pubmed.ncbi.nlm.nih.gov/16870913/) |
| 91 | Cole FS, Alleyne C, Barks JD, et al. NIH consensus development conference statement: inhaled nitric-oxide therapy for premature infants. *Pediatrics* 2011. | [21220405](https://pubmed.ncbi.nlm.nih.gov/21220405/) |
| 92 | Krishnan U, Feinstein JA, Adatia I, et al. Evaluation and management of pulmonary hypertension in children with bronchopulmonary dysplasia. *J Pediatr* 2017. | [28645441](https://pubmed.ncbi.nlm.nih.gov/28645441/) |
| 93 | Abman SH, Hansmann G, Archer SL, et al. Pediatric pulmonary hypertension: guidelines from the American Heart Association and American Thoracic Society. *Circulation* 2015. | [26534956](https://pubmed.ncbi.nlm.nih.gov/26534956/) |
| 94 | Mourani PM, Sontag MK, Younoszai A, et al. Early pulmonary vascular disease in preterm infants at risk for bronchopulmonary dysplasia. *Am J Respir Crit Care Med* 2015. | [25389562](https://pubmed.ncbi.nlm.nih.gov/25389562/) |
| 95 | Khemani E, McElhinney DB, Rhein L, et al. Pulmonary artery hypertension in formerly premature infants with bronchopulmonary dysplasia: clinical features and outcomes in the surfactant era. *Pediatrics* 2007. | [18055675](https://pubmed.ncbi.nlm.nih.gov/18055675/) |
| 96 | An HS, Bae EJ, Kim GB, et al. Pulmonary hypertension in preterm infants with bronchopulmonary dysplasia. *Korean Circ J* 2010. | [20339498](https://pubmed.ncbi.nlm.nih.gov/20339498/) |
| 97 | Berkelhamer SK, Mestan KK, Steinhorn R. An update on the diagnosis and management of bronchopulmonary dysplasia (BPD)-associated pulmonary hypertension. *Semin Perinatol* 2018. | [30384985](https://pubmed.ncbi.nlm.nih.gov/30384985/) |
| 98 | Nyp M, Sandritter T, Poppinga N, et al. Sildenafil citrate, bronchopulmonary dysplasia and disordered pulmonary gas exchange. *J Perinatol* 2012. | [21941230](https://pubmed.ncbi.nlm.nih.gov/21941230/) |
| 99 | Mukherjee A, Dombi T, Wittke B, Lalonde R. Population pharmacokinetics of sildenafil in term neonates. *Clin Pharmacol Ther* 2009. | [18800037](https://pubmed.ncbi.nlm.nih.gov/18800037/) |

**→ model.** `PVR` is written as `(1/CAP_rel^PXS) · (1 + KHPV·hpv) · (1 + KVR·VREM)`
— an **area** term, a **tone** term and a **remodelling** term. That structure
is the model's explanation of #94 and #95: BPD-PH is not a separate disease
downstream of BPD, it is the same microvascular arrest read through a different
instrument, which is why it appears early (#94) and why it carries the mortality
it does (#95, #96). A PDE5 inhibitor enters only the tone term; the sildenafil
PK is from #99 with the acknowledged caveat below. #90 and #91 are in the list
as the negative iNO result and are drawn as a dashed edge to `Secondary_septa` —
iNO improves the tone term without improving the 36-week endpoint, exactly as
sildenafil does in experiment A7.

**PK caveat, stated because it matters.** #99 is a **term**-neonate population
PK study. Preterm sildenafil disposition is less well characterised, and the
model's `SIL_F = 0.35`, `SIL_CL = 2.0 L/day/kg`, `SIL_VD = 1.2 L/kg` (giving a
half-life near 10 h and a steady state near 0.5 mg/L on 3 mg/kg/day) are a
plausible extrapolation, not a fitted result. The A7 conclusion depends on the
*structure* of the PVR equation rather than on this PK, but the absolute
sildenafil concentrations should not be used for anything.

---

## 13. Patent ductus arteriosus — exposure-days, not presence

| # | Reference | PMID |
|---|-----------|------|
| 100 | Benitz WE; Committee on Fetus and Newborn. Patent ductus arteriosus in preterm infants. *Pediatrics* 2016. | [26672023](https://pubmed.ncbi.nlm.nih.gov/26672023/) |
| 101 | Mirza H, Garcia J, McKinley G, et al. Duration of significant patent ductus arteriosus and bronchopulmonary dysplasia in extremely preterm infants. *J Perinatol* 2019. | [31554913](https://pubmed.ncbi.nlm.nih.gov/31554913/) |
| 102 | Clyman RI, Hills NK, Liebowitz M, Johng S. Relationship between duration of infant exposure to a moderate-to-large patent ductus arteriosus shunt and the risk of prolonged ventilation. *Am J Perinatol* 2020. | [31600791](https://pubmed.ncbi.nlm.nih.gov/31600791/) |
| 103 | Schena F, Ciarmoli E, Mosca F. Association between hemodynamically significant patent ductus arteriosus and bronchopulmonary dysplasia. *J Pediatr* 2015. | [25882876](https://pubmed.ncbi.nlm.nih.gov/25882876/) |

**→ model.** The PDA is parameterised as an **interval** (`PDA_T0`, `PDA_T1`)
rather than a flag, precisely because #101–#103 report exposure-duration
effects. Scenario S05 closes the duct on day 3 and the benefit comes entirely
from the shortened interval — there is no separate "PDA closed" bonus term.

---

## 14. Diuretics and bronchodilators — the drugs that move a number

| # | Reference | PMID |
|---|-----------|------|
| 104 | Stewart A, Brion LP. Intravenous or enteral loop diuretics for preterm infants with (or developing) chronic lung disease. *Cochrane Database Syst Rev* 2011. | [21901676](https://pubmed.ncbi.nlm.nih.gov/21901676/) |
| 105 | Stewart A, Brion LP, Ambrosio-Perez I. Diuretics acting on the distal renal tubule for preterm infants with (or developing) chronic lung disease. *Cochrane Database Syst Rev* 2011. | [21901679](https://pubmed.ncbi.nlm.nih.gov/21901679/) |
| 106 | Greenberg RG, Gayam S, Savage D, et al. Furosemide exposure and prevention of bronchopulmonary dysplasia in premature infants. *J Pediatr* 2019. | [30579586](https://pubmed.ncbi.nlm.nih.gov/30579586/) |
| 107 | Ng G, Ohlsson A. Inhaled bronchodilators for the prevention and treatment of chronic lung disease in preterm infants. *Cochrane Database Syst Rev* 2024. | [38591664](https://pubmed.ncbi.nlm.nih.gov/38591664/) |

**→ model.** Furosemide acts on exactly one thing, `KOUT_LW`, and therefore
enters the model through the *shunt* term of the oxygen requirement, which has
no reserve. It improves `FIO2` within a day (#104: short-term improvement in
lung mechanics) and leaves `ALV` and `LOSTW` essentially unchanged (#106: no
BPD prevention). Experiment A5 prints both halves of that. The bronchodilator
(#107) is the same shape and acts only on `ASM`.

---

## 15. Nutrition, growth and IGF-1 replacement

| # | Reference | PMID |
|---|-----------|------|
| 108 | Ehrenkranz RA, Dusick AM, Vohr BR, et al. Growth in the neonatal intensive care unit influences neurodevelopmental and growth outcomes of extremely low birth weight infants. *Pediatrics* 2006. | [16585322](https://pubmed.ncbi.nlm.nih.gov/16585322/) |
| 109 | Ley D, Hallberg B, Hansen-Pupp I, et al. rhIGF-1/rhIGFBP-3 in preterm infants: a phase 2 randomized controlled trial. *J Pediatr* 2019. | [30471715](https://pubmed.ncbi.nlm.nih.gov/30471715/) |

**→ model.** Weight is a state (`WT`) rather than a covariate because it does
three jobs at once: it sets the mg/kg dose, it scales every clearance, and poor
growth feeds back on `IGF1` and hence on `G`. That closes a loop that a
fixed-weight model cannot represent — a steroid that buys lung function also
suppresses growth (`AW_DEX`), which suppresses IGF-1, which closes the gate the
steroid was meant to open. #109 supplies the `RHIGF` arm.

---

## 16. Mesenchymal stromal cells

| # | Reference | PMID |
|---|-----------|------|
| 110 | Chang YS, Ahn SY, Yoo HS, et al. Mesenchymal stem cells for bronchopulmonary dysplasia: phase 1 dose-escalation clinical trial. *J Pediatr* 2014. | [24508444](https://pubmed.ncbi.nlm.nih.gov/24508444/) |
| 111 | Ahn SY, Chang YS, Kim JH, et al. Two-year follow-up outcomes of premature infants enrolled in the phase I trial of mesenchymal stem cells transplantation for bronchopulmonary dysplasia. *J Pediatr* 2017. | [28341525](https://pubmed.ncbi.nlm.nih.gov/28341525/) |

**→ model.** MSC therapy is implemented as a decaying paracrine effect
(`MSC_EMAX`, `MSC_TAU`) on `VEGF` and `EPC` rather than as engraftment, which
matches the current mechanistic understanding. It is a **programme-class**
intervention in the model's taxonomy and is deliberately shown as
early-sensitive for the same reason vitamin A is.

---

## 17. Genetic susceptibility

| # | Reference | PMID |
|---|-----------|------|
| 112 | Bhandari V, Bizzarro MJ, Shetty A, et al. Familial and genetic susceptibility to major neonatal morbidities in preterm twins. *Pediatrics* 2006. | [16740829](https://pubmed.ncbi.nlm.nih.gov/16740829/) |
| 113 | Hadchouel A, Durrmeyer X, Bouzigon E, et al. Identification of SPOCK2 as a susceptibility gene for bronchopulmonary dysplasia. *Am J Respir Crit Care Med* 2011. | [21836138](https://pubmed.ncbi.nlm.nih.gov/21836138/) |
| 114 | Wang H, St Julien KR, Stevenson DK, et al. A genome-wide association study (GWAS) for bronchopulmonary dysplasia. *Pediatrics* 2013. | [23897914](https://pubmed.ncbi.nlm.nih.gov/23897914/) |

**→ model.** `GENRISK` shifts the NF-κB set-point only. That is a deliberately
modest representation: #112 reports substantial heritability but #114 found no
genome-wide significant common variant, so the model encodes susceptibility as
a set-point shift rather than as a named pathway lesion.

---

## 18. Antenatal determinants

| # | Reference | PMID |
|---|-----------|------|
| 115 | Hansen AR, Barnés CM, Folkman J, McElrath TF. Maternal preeclampsia predicts the development of bronchopulmonary dysplasia. *J Pediatr* 2010. | [20004912](https://pubmed.ncbi.nlm.nih.gov/20004912/) |
| 116 | Bose C, Van Marter LJ, Laughon M, et al. Fetal growth restriction and chronic lung disease among infants born before the 28th week of gestation. *Pediatrics* 2009. | [19706590](https://pubmed.ncbi.nlm.nih.gov/19706590/) |
| 117 | Morrow LA, Wagner BD, Ingram DA, et al. Antenatal determinants of bronchopulmonary dysplasia and late respiratory disease in preterm infants. *Am J Respir Crit Care Med* 2017. | [28249118](https://pubmed.ncbi.nlm.nih.gov/28249118/) |
| 118 | Keller RL, Feng R, DeMauro SB, et al. Bronchopulmonary dysplasia and perinatal characteristics predict 1-year respiratory outcomes in newborns born at extremely low gestational age. *J Pediatr* 2017. | [28528221](https://pubmed.ncbi.nlm.nih.gov/28528221/) |

**→ model.** These are HIT #1: the antenatal `sflt` term and the `IUGR`
multipliers exist because #115 and #116 report that an anti-angiogenic
intrauterine environment and fetal growth restriction predict BPD independently
of gestational age. In the model this means `G(0) < 1` — the infant arrives with
the gate already partly closed, which is the single biggest determinant of the
irreducible floor computed in experiment A10.

---

## 19. Long-term trajectory — why LOSTW and not FiO2 is the readout

| # | Reference | PMID |
|---|-----------|------|
| 119 | Islam JY, Keller RL, Aschner JL, et al. Understanding the short- and long-term respiratory outcomes of prematurity and bronchopulmonary dysplasia. *Am J Respir Crit Care Med* 2015. | [26038806](https://pubmed.ncbi.nlm.nih.gov/26038806/) |
| 120 | Baraldi E, Filippone M. Chronic lung disease after premature birth. *N Engl J Med* 2007. | [17989387](https://pubmed.ncbi.nlm.nih.gov/17989387/) |
| 121 | Doyle LW, Andersson S, Bush A, et al. Expiratory airflow in late adolescence and early adulthood in individuals born very preterm or with very low birthweight compared with controls. *Lancet Respir Med* 2019. | [31078498](https://pubmed.ncbi.nlm.nih.gov/31078498/) |
| 122 | Simpson SJ, Turkovic L, Wilson AC, et al. Lung function trajectories throughout childhood in survivors of very preterm birth: a longitudinal cohort study. *Lancet Child Adolesc Health* 2018. | [30169268](https://pubmed.ncbi.nlm.nih.gov/30169268/) |
| 123 | Vom Hove M, Prenzel F, Uhlig HH, Robel-Tillig E. Pulmonary outcome in former preterm, very low birth weight children with bronchopulmonary dysplasia: a case-control follow-up at school age. *J Pediatr* 2014. | [24055328](https://pubmed.ncbi.nlm.nih.gov/24055328/) |
| 124 | Wong PM, Lees AN, Louw J, et al. Emphysema in young adult survivors of moderate-to-severe bronchopulmonary dysplasia. *Eur Respir J* 2008. | [18385172](https://pubmed.ncbi.nlm.nih.gov/18385172/) |
| 125 | Doyle LW; Victorian Infant Collaborative Study Group. Bronchopulmonary dysplasia in very low birthweight subjects and lung function in late adolescence. *Pediatrics* 2006. | [16818555](https://pubmed.ncbi.nlm.nih.gov/16818555/) |
| 126 | Boucherat O, Morissette MC, Provencher S, et al. Bridging lung development with chronic obstructive pulmonary disease: relevance of developmental pathways in COPD pathogenesis. *Am J Respir Crit Care Med* 2016. | [26681127](https://pubmed.ncbi.nlm.nih.gov/26681127/) |
| 127 | Davidson LM, Berkelhamer SK. Bronchopulmonary dysplasia: chronic lung disease of infancy and long-term pulmonary outcomes. *J Clin Med* 2017. | [28067830](https://pubmed.ncbi.nlm.nih.gov/28067830/) |

**→ model.** This block is the justification for carrying `LOSTW` as a separate
state. The 36-week support level is the *definition* of BPD but a poor measure
of the developmental deficit, because the lung has large reserve — which is why
the model's oxygen requirement has an explicit reserve threshold (`DTHR`) below
which FiO2 does not move at all. What predicts the adult FEV1 plateau (#121,
#122, #124, #125) is the surface that was never built, and that is `LOSTW`. #126
is the link that makes the deficit matter fifty years later.

---

## 20. How to reproduce every number this model claims

```bash
cd bronchopulmonary-dysplasia
python3 bpd_reference_impl.py          # all experiments A0-A10, ~3 min
python3 bpd_reference_impl.py A3       # just the caffeine decomposition
```

`bpd_reference_impl.py` is a dependency-free (standard library only, fixed-step
RK4) reimplementation of the same equations as `bpd_mrgsolve_model.R`. It is the
authoritative source for every figure quoted in this directory's `README.md`. If
a number in the README and a number printed by that script ever disagree, the
script is right and the README is stale.

Tool references used by this library as a whole:

- mrgsolve: <https://mrgsolve.org/>
- QSP in R: <https://vantage-research.net/qsp-in-r/>
- gPKPDviz (mrgsolve-based Shiny PK/PD tooling):
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/>

---

## ⚠️ Disclaimer

This model is for **education and research only**. It is a mechanistic
simplification calibrated to reproduce the direction and rough magnitude of
published trial results; it has not been fitted to individual patient data,
validated prospectively, or reviewed by a regulator. Neonatal dosing in
particular is quoted here to explain model parameters and **must not** be used
for prescribing. Nothing in this directory should inform the care of a real
infant.
