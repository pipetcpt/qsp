# Mitral Regurgitation — QSP model references

**116 references, every one resolved against PubMed programmatically.**

The bibliography was not written from memory. Each entry began as a query
against the NCBI E-utilities API; the author, title, journal, year and PMID
below are the values NCBI returned. Where a query's top hit turned out to be
a paper other than the one intended, the entry has been re-labelled to
describe what it actually is rather than silently mis-cited, and two junk
hits (a conference-abstract compilation and a single case report) were
dropped. Regenerate or re-verify with the scripts noted at the foot of this
file.

The organising idea of the model is that a regurgitant volume means nothing
until it is divided by something, and that the disease is really a family of
diseases indexed by the denominator. The sections below follow the five
denominators.

| | denominator | what it decides | section |
|---|---|---|---|
| 1 | operating LA compliance `C_op` | congestion; acute vs chronic | §3 |
| 2 | total stroke volume | regurgitant fraction, the grade | §2 |
| 3 | LV end-diastolic volume | valve or ventricle; COAPT vs MITRA-FR | §1 |
| 4 | the afterload the leak removes | hidden contractility, the EF 60% rule | §4 |
| 5 | `k_PISA` | the measurement's own bias | §2 |

---

## 1. The two pivotal trials, and the proportionality argument they created

*The model is calibrated on the COAPT control arm alone; every other trial number in this list is an out-of-sample prediction target.*

- Stone GW, Lindenfeld J, Abraham WT, et al. **Transcatheter Mitral-Valve Repair in Patients with Heart Failure.** *N Engl J Med* 2018. [PMID 30280640](https://pubmed.ncbi.nlm.nih.gov/30280640/)  
  COAPT. The trial the model is calibrated against: HF hospitalisation 35.8 vs 67.9 per patient-year, death 29.1% vs 46.1% at 24 months.
- Obadia JF, Messika-Zeitoun D, Leurent G, et al. **Percutaneous Repair or Medical Treatment for Secondary Mitral Regurgitation.** *N Engl J Med* 2018. [PMID 30145927](https://pubmed.ncbi.nlm.nih.gov/30145927/)  
  MITRA-FR. Same device, opposite result: death or HF hospitalisation 54.6% vs 51.3% at 12 months.
- Grayburn PA, Sannino A, Packer M. **Proportionate and Disproportionate Functional Mitral Regurgitation: A New Conceptual Framework That Reconciles the Results of the MITRA-FR and COAPT Trials.** *JACC Cardiovasc Imaging* 2019. [PMID 30553663](https://pubmed.ncbi.nlm.nih.gov/30553663/)  
  The proportionality framework. This is the paper the model turns into arithmetic.
- Packer M, Grayburn PA. **Contrasting Effects of Pharmacological, Procedural, and Surgical Interventions on Proportionate and Disproportionate Functional Mitral Regurgitation in Chronic Heart Failure.** *Circulation* 2019. [PMID 31070944](https://pubmed.ncbi.nlm.nih.gov/31070944/)  
  Why drug and device effects on functional MR are not interchangeable.
- Packer M, Grayburn PA. **New Evidence Supporting a Novel Conceptual Framework for Distinguishing Proportionate and Disproportionate Functional Mitral Regurgitation.** *JAMA Cardiol* 2020. [PMID 32074243](https://pubmed.ncbi.nlm.nih.gov/32074243/)  
  The framework stated as a testable hypothesis.
- Mack MJ, Lindenfeld J, Abraham WT, et al. **3-Year Outcomes of Transcatheter Mitral Valve Repair in Patients With Heart Failure.** *J Am Coll Cardiol* 2021. [PMID 33632476](https://pubmed.ncbi.nlm.nih.gov/33632476/)  
  Durability of the COAPT effect.
- Stone GW, Abraham WT, Lindenfeld J, et al. **Five-Year Follow-up after Transcatheter Repair of Secondary Mitral Regurgitation.** *N Engl J Med* 2023. [PMID 36876756](https://pubmed.ncbi.nlm.nih.gov/36876756/)  
  Five-year COAPT follow-up.
- Anker SD, Friede T, von Bardeleben RS, et al. **Percutaneous repair of moderate-to-severe or severe functional mitral regurgitation in patients with symptomatic heart failure: Baseline characteristics of patients in the RESHAPE-HF2 trial and comparison to COAPT and MITRA-FR trials.** *Eur J Heart Fail* 2024. [PMID 38847420](https://pubmed.ncbi.nlm.nih.gov/38847420/)  
  A third randomised trial in the same space.
- Pibarot P, Delgado V, Bax JJ. **MITRA-FR vs. COAPT: lessons from two trials with diametrically opposed results.** *Eur Heart J Cardiovasc Imaging* 2019. [PMID 31115470](https://pubmed.ncbi.nlm.nih.gov/31115470/)  
  Why the two trials differed: contemporaneous editorial analysis.

## 2. Quantifying the lesion: EROA, PISA, and the geometry PISA assumes

*Denominator 5. The severity threshold everyone uses was calibrated on a round orifice and is then applied to a crescentic one.*

- Enriquez-Sarano M, Seward JB, Bailey KR, et al. **Effective regurgitant orifice area: a noninvasive Doppler development of an old hemodynamic concept.** *J Am Coll Cardiol* 1994. [PMID 8294699](https://pubmed.ncbi.nlm.nih.gov/8294699/)  
  The effective regurgitant orifice area as a clinical quantity.
- Enriquez-Sarano M, Avierinos JF, Messika-Zeitoun D, et al. **Quantitative determinants of the outcome of asymptomatic mitral regurgitation.** *N Engl J Med* 2005. [PMID 15745978](https://pubmed.ncbi.nlm.nih.gov/15745978/)  
  EROA predicts outcome in primary MR; the origin of the 0.4 cm2 threshold.
- Bargiggia GS, Tronconi L, Sahn DJ, et al. **A new method for quantitation of mitral regurgitation based on color flow Doppler imaging of flow convergence proximal to regurgitant orifice.** *Circulation* 1991. [PMID 1914090](https://pubmed.ncbi.nlm.nih.gov/1914090/)  
  PISA: the proximal isovelocity surface area method.
- Grayburn PA, Weissman NJ, Zamorano JL. **Quantitation of mitral regurgitation.** *Circulation* 2012. [PMID 23071176](https://pubmed.ncbi.nlm.nih.gov/23071176/)  
  Review of what each method actually measures.
- Biner S, Rafique A, Rafii F, et al. **Reproducibility of proximal isovelocity surface area, vena contracta, and regurgitant jet area for assessment of mitral regurgitation severity.** *JACC Cardiovasc Imaging* 2010. [PMID 20223419](https://pubmed.ncbi.nlm.nih.gov/20223419/)  
  Reproducibility of the measurements the trials relied on.
- Utsunomiya T, Ogawa T, Tang HA, et al. **Doppler color flow mapping of the proximal isovelocity surface area: a new method for measuring volume flow rate across a narrowed orifice.** *J Am Soc Echocardiogr* 1991. [PMID 1910832](https://pubmed.ncbi.nlm.nih.gov/1910832/)  
  PISA validated in vitro, where the orifice really is circular.
- Zoghbi WA, Adams D, Bonow RO, et al. **Recommendations for Noninvasive Evaluation of Native Valvular Regurgitation: A Report from the American Society of Echocardiography Developed in Collaboration with the Society for Cardiovascular Magnetic Resonance.** *J Am Soc Echocardiogr* 2017. [PMID 28314623](https://pubmed.ncbi.nlm.nih.gov/28314623/)  
  The guideline grading scheme and its own caveats.
- Uretsky S, Gillam L, Lang R, et al. **Discordance between echocardiography and MRI in the assessment of mitral regurgitation severity: a prospective multicenter trial.** *J Am Coll Cardiol* 2015. [PMID 25790878](https://pubmed.ncbi.nlm.nih.gov/25790878/)  
  Echo and CMR disagree; the model treats this as information.
- Dujardin KS, Enriquez-Sarano M, Bailey KR, et al. **Grading of mitral regurgitation by quantitative Doppler echocardiography: calibration by left ventricular angiography in routine clinical practice.** *Circulation* 1997. [PMID 9396435](https://pubmed.ncbi.nlm.nih.gov/9396435/)  
  Quantitative Doppler grading and its cut-points.
- Militaru S, Bonnefous O, Hami K, et al. **Validation of Semiautomated Quantification of Mitral Valve Regurgitation by Three-Dimensional Color Doppler Transesophageal Echocardiography.** *J Am Soc Echocardiogr* 2020. [PMID 32143780](https://pubmed.ncbi.nlm.nih.gov/32143780/)  
  Cross-modality validation.
- Sugeng L, Chandra S, Lang RM. **Three-dimensional echocardiography for assessment of mitral valve regurgitation.** *Curr Opin Cardiol* 2009. [PMID 19593121](https://pubmed.ncbi.nlm.nih.gov/19593121/)  
  Three-dimensional echocardiography of the regurgitant orifice, which is how its non-circular shape became visible.
- Ashikhmina E, Shook D, Cobey F, et al. **Three-dimensional versus two-dimensional echocardiographic assessment of functional mitral regurgitation proximal isovelocity surface area.** *Anesth Analg* 2015. [PMID 25166465](https://pubmed.ncbi.nlm.nih.gov/25166465/)  
  Three- versus two-dimensional assessment of functional MR: the measurement disagreement, quantified.

## 3. Acute versus chronic, and left atrial compliance

*Denominator 1. The same regurgitant volume floods one lung and is silent in another, and the difference is a property of the atrium.*

- BRAUNWALD E, AWE WC. **The syndrome of severe mitral regurgitation with normal left atrial pressure.** *Circulation* 1963. [PMID 14015085](https://pubmed.ncbi.nlm.nih.gov/14015085/)  
  The observation the model reproduces: severe MR can coexist with a normal atrial pressure.
- Braunwald E. **Mitral regurgitation: physiologic, clinical and surgical considerations.** *N Engl J Med* 1969. [PMID 5797187](https://pubmed.ncbi.nlm.nih.gov/5797187/)  
  The physiology of the acute-versus-chronic distinction, stated early and clearly.
- Nishimura RA, Schaff HV, Shub C, et al. **Papillary muscle rupture complicating acute myocardial infarction: analysis of 17 patients.** *Am J Cardiol* 1983. [PMID 6823851](https://pubmed.ncbi.nlm.nih.gov/6823851/)  
  Acute severe MR from papillary rupture: analysis of the clinical syndrome.
- Dernellis JM, Stefanadis CI, Zacharoulis AA, et al. **Left atrial mechanical adaptation to long-standing hemodynamic loads based on pressure-volume relations.** *Am J Cardiol* 1998. [PMID 9605056](https://pubmed.ncbi.nlm.nih.gov/9605056/)  
  The left atrial pressure-volume relation: denominator 1.
- Estévez-Loureiro R, Lorusso R, Taramasso M, et al. **Management of Severe Mitral Regurgitation in Patients With Acute Myocardial Infarction: JACC Focus Seminar 2/5.** *J Am Coll Cardiol* 2024. [PMID 38692830](https://pubmed.ncbi.nlm.nih.gov/38692830/)  
  Management of acute severe MR.
- Welch TD, Coylewright M, Powell BD, et al. **Symptomatic pulmonary hypertension with giant left atrial v waves after surgical maze procedures: evaluation by comprehensive hemodynamic catheterization.** *Heart Rhythm* 2013. [PMID 24050987](https://pubmed.ncbi.nlm.nih.gov/24050987/)  
  Giant left atrial v waves and their pulmonary consequences -- denominator 1 seen directly in a patient.
- Yellin EL, Yoran C, Sonnenblick EH, et al. **Dynamic changes in the canine mitral regurgitant orifice area during ventricular ejection.** *Circ Res* 1979. [PMID 487530](https://pubmed.ncbi.nlm.nih.gov/487530/)  
  The regurgitant orifice area changes DURING ejection -- the dynamic orifice, measured in 1979.

## 4. Ejection fraction, hidden contractility, and the operative thresholds

*Denominator 4. EF is computed against a stroke volume that contains the leak, which is why the surgical trigger sits at 60% and not 50%.*

- Tribouilloy C, Grigioni F, Avierinos JF, et al. **Survival implication of left ventricular end-systolic diameter in mitral regurgitation due to flail leaflets a long-term follow-up multicenter study.** *J Am Coll Cardiol* 2009. [PMID 19909877](https://pubmed.ncbi.nlm.nih.gov/19909877/)  
  The LV end-systolic diameter trigger.
- Wisenbaugh T, Skudicky D, Sareli P. **Prediction of outcome after valve replacement for rheumatic mitral regurgitation in the era of chordal preservation.** *Circulation* 1994. [PMID 8281646](https://pubmed.ncbi.nlm.nih.gov/8281646/)  
  Post-operative outcome and the ventricle's hidden state.
- Carabello BA. **The current therapy for mitral regurgitation.** *J Am Coll Cardiol* 2008. [PMID 18652937](https://pubmed.ncbi.nlm.nih.gov/18652937/)  
  Carabello's synthesis of the loading problem.
- Zile MR, Gaasch WH, Carroll JD, et al. **Chronic mitral regurgitation: predictive value of preoperative echocardiographic indexes of left ventricular function and wall stress.** *J Am Coll Cardiol* 1984. [PMID 6693615](https://pubmed.ncbi.nlm.nih.gov/6693615/)  
  Preoperative indices that survive the afterload problem.
- Suga H, Sagawa K, Shoukas AA. **Load independence of the instantaneous pressure-volume ratio of the canine left ventricle and effects of epinephrine and heart rate on the ratio.** *Circ Res* 1973. [PMID 4691336](https://pubmed.ncbi.nlm.nih.gov/4691336/)  
  The elastance concept the model's ventricle is built on.
- Sunagawa K, Maughan WL, Burkhoff D, et al. **Left ventricular interaction with arterial load studied in isolated canine ventricle.** *Am J Physiol* 1983. [PMID 6638199](https://pubmed.ncbi.nlm.nih.gov/6638199/)  
  Ventricular-arterial coupling: where E_a comes from.
- Senzaki H, Chen CH, Kass DA. **Single-beat estimation of end-systolic pressure-volume relation in humans. A new method with the potential for noninvasive application.** *Circulation* 1996. [PMID 8921794](https://pubmed.ncbi.nlm.nih.gov/8921794/)  
  Single-beat estimation of E_es, and its normal range.
- Chen CH, Fetics B, Nevo E, et al. **Noninvasive single-beat determination of left ventricular end-systolic elastance in humans.** *J Am Coll Cardiol* 2001. [PMID 11738311](https://pubmed.ncbi.nlm.nih.gov/11738311/)  
  Non-invasive E_es, used for the model's baseline value.
- Enriquez-Sarano M, Tajik AJ, Schaff HV, et al. **Echocardiographic prediction of survival after surgical correction of organic mitral regurgitation.** *Circulation* 1994. [PMID 8044955](https://pubmed.ncbi.nlm.nih.gov/8044955/)  
  The origin of the EF 60% operative threshold.
- Gaynor JW, Feneley MP, Gall SA Jr, et al. **Left ventricular adaptation to aortic regurgitation in conscious dogs.** *J Thorac Cardiovasc Surg* 1997. [PMID 9011684](https://pubmed.ncbi.nlm.nih.gov/9011684/)  
  Volume-overload adaptation measured in a conscious large-animal model (aortic regurgitation); the ventricular response to a regurgitant load.
- Carabello BA, Nolan SP, McGuire LB. **Assessment of preoperative left ventricular function in patients with mitral regurgitation: value of the end-systolic wall stress-end-systolic volume ratio.** *Circulation* 1981. [PMID 7296794](https://pubmed.ncbi.nlm.nih.gov/7296794/)  
  Why EF misleads in MR: the wall-stress correction.
- Klotz S, Hay I, Dickstein ML, et al. **Single-beat estimation of end-diastolic pressure-volume relationship: a novel method with potential for noninvasive application.** *Am J Physiol Heart Circ Physiol* 2006. [PMID 16428349](https://pubmed.ncbi.nlm.nih.gov/16428349/)  
  The EDPVR form the model uses.
- Burkhoff D, Mirsky I, Suga H. **Assessment of systolic and diastolic ventricular properties via pressure-volume analysis: a guide for clinical, translational, and basic researchers.** *Am J Physiol Heart Circ Physiol* 2005. [PMID 16014610](https://pubmed.ncbi.nlm.nih.gov/16014610/)  
  Pressure-volume analysis as used here.

## 5. Secondary MR: tethering, annular geometry, and the dynamic orifice

*Where the functional orifice comes from, and the negative results that constrain any mechanism -- annular dilation alone is often not enough.*

- Otsuji Y, Handschumacher MD, Schwammenthal E, et al. **Insights from three-dimensional echocardiography into the mechanism of functional mitral regurgitation: direct in vivo demonstration of altered leaflet tethering geometry.** *Circulation* 1997. [PMID 9323092](https://pubmed.ncbi.nlm.nih.gov/9323092/)  
  The tethering mechanism demonstrated in vivo.
- Yiu SF, Enriquez-Sarano M, Tribouilloy C, et al. **Determinants of the degree of functional mitral regurgitation in patients with systolic left ventricular dysfunction: A quantitative clinical study.** *Circulation* 2000. [PMID 10993859](https://pubmed.ncbi.nlm.nih.gov/10993859/)  
  What actually determines functional MR severity.
- Kwan J, Shiota T, Agler DA, et al. **Geometric differences of the mitral apparatus between ischemic and dilated cardiomyopathy with significant mitral regurgitation: real-time three-dimensional echocardiography study.** *Circulation* 2003. [PMID 12615791](https://pubmed.ncbi.nlm.nih.gov/12615791/)  
  Ischaemic versus dilated geometry.
- Levine RA, Schwammenthal E. **Ischemic mitral regurgitation on the threshold of a solution: from paradoxes to unifying concepts.** *Circulation* 2005. [PMID 16061756](https://pubmed.ncbi.nlm.nih.gov/16061756/)  
  The unifying geometric account of ischaemic MR.
- Komeda M, Glasson JR, Bolger AF, et al. **Geometric determinants of ischemic mitral regurgitation.** *Circulation* 1997. [PMID 9386087](https://pubmed.ncbi.nlm.nih.gov/9386087/)  
  Geometric determinants, quantified.
- Hung J, Otsuji Y, Handschumacher MD, et al. **Mechanism of dynamic regurgitant orifice area variation in functional mitral regurgitation: physiologic insights from the proximal flow convergence technique.** *J Am Coll Cardiol* 1999. [PMID 9973036](https://pubmed.ncbi.nlm.nih.gov/9973036/)  
  Why the functional orifice changes with load.
- Silbiger JJ. **Mechanistic insights into atrial functional mitral regurgitation: Far more complicated than just left atrial remodeling.** *Echocardiography* 2019. [PMID 30620100](https://pubmed.ncbi.nlm.nih.gov/30620100/)  
  Atrial functional MR: the ventricle need not be involved at all.
- Deferm S, Bertrand PB, Verbrugge FH, et al. **Atrial Functional Mitral Regurgitation: JACC Review Topic of the Week.** *J Am Coll Cardiol* 2019. [PMID 31097168](https://pubmed.ncbi.nlm.nih.gov/31097168/)  
  Atrial functional MR as a distinct entity.
- Kaplan SR, Bashein G, Sheehan FH, et al. **Three-dimensional echocardiographic assessment of annular shape changes in the normal and regurgitant mitral valve.** *Am Heart J* 2000. [PMID 10689248](https://pubmed.ncbi.nlm.nih.gov/10689248/)  
  Annular shape change with regurgitation.
- Gorman JH 3rd, Gupta KB, Streicher JT, et al. **Dynamic three-dimensional imaging of the mitral valve and left ventricle by rapid sonomicrometry array localization.** *J Thorac Cardiovasc Surg* 1996. [PMID 8800160](https://pubmed.ncbi.nlm.nih.gov/8800160/)  
  Dynamic mitral apparatus imaging.
- Schwammenthal E, Chen C, Benning F, et al. **Dynamics of mitral regurgitant flow and orifice area. Physiologic application of the proximal flow convergence method: clinical data and experimental testing.** *Circulation* 1994. [PMID 8026013](https://pubmed.ncbi.nlm.nih.gov/8026013/)  
  The orifice is dynamic within systole.
- Ormiston JA, Shah PM, Tei C, et al. **Size and motion of the mitral valve annulus in man. I. A two-dimensional echocardiographic method and findings in normal subjects.** *Circulation* 1981. [PMID 7237707](https://pubmed.ncbi.nlm.nih.gov/7237707/)  
  Normal annular size and motion.
- Sanfilippo AJ, Weyman AE, Levine RA. **The problem of echocardiographic detection of mitral valve prolapse and determination of its true prevalence.** *Herz* 1988. [PMID 3053381](https://pubmed.ncbi.nlm.nih.gov/3053381/)  
  Echocardiographic detection of mitral valve prolapse and the geometry of the annulus.
- Otsuji Y, Kumanohoso T, Yoshifuku S, et al. **Isolated annular dilation does not usually cause important functional mitral regurgitation: comparison between patients with lone atrial fibrillation and those with idiopathic or ischemic cardiomyopathy.** *J Am Coll Cardiol* 2002. [PMID 12020493](https://pubmed.ncbi.nlm.nih.gov/12020493/)  
  A key negative result: annular dilation alone is often not enough.
- He Z, Ritchie J, Grashow JS, et al. **In vitro dynamic strain behavior of the mitral valve posterior leaflet.** *J Biomech Eng* 2005. [PMID 16060357](https://pubmed.ncbi.nlm.nih.gov/16060357/)  
  Leaflet mechanics.

## 6. Leaflet plasticity: the valve is not a passive gasket

*The mitral leaflets grow under tethering, but only partially, which is why the SPEED of ventricular dilation matters and not only its size.*

- Dal-Bianco JP, Aikawa E, Bischoff J, et al. **Active adaptation of the tethered mitral valve: insights into a compensatory mechanism for functional mitral regurgitation.** *Circulation* 2009. [PMID 19597052](https://pubmed.ncbi.nlm.nih.gov/19597052/)  
  Leaflets actively grow under tethering: the model's leaflet-supply state.
- Dal-Bianco JP, Aikawa E, Bischoff J, et al. **Myocardial Infarction Alters Adaptation of the Tethered Mitral Valve.** *J Am Coll Cardiol* 2016. [PMID 26796392](https://pubmed.ncbi.nlm.nih.gov/26796392/)  
  Infarction alters the adaptation, which is why fast dilation is worse.
- Beaudoin J, Dal-Bianco JP, Aikawa E, et al. **Mitral Leaflet Changes Following Myocardial Infarction: Clinical Evidence for Maladaptive Valvular Remodeling.** *Circ Cardiovasc Imaging* 2017. [PMID 29042413](https://pubmed.ncbi.nlm.nih.gov/29042413/)  
  Clinical evidence of maladaptive leaflet remodelling.
- Chaput M, Handschumacher MD, Tournoux F, et al. **Mitral leaflet adaptation to ventricular remodeling: occurrence and adequacy in patients with functional mitral regurgitation.** *Circulation* 2008. [PMID 18678770](https://pubmed.ncbi.nlm.nih.gov/18678770/)  
  Adaptation occurs but is inadequate: the phi<1 assumption.
- Rausch MK, Tibayan FA, Ingels NB Jr, et al. **Mechanics of the mitral annulus in chronic ischemic cardiomyopathy.** *Ann Biomed Eng* 2013. [PMID 23636575](https://pubmed.ncbi.nlm.nih.gov/23636575/)  
  Annular mechanics in chronic ischaemic disease.
- Bischoff J, Casanovas G, Wylie-Sears J, et al. **CD45 Expression in Mitral Valve Endothelial Cells After Myocardial Infarction.** *Circ Res* 2016. [PMID 27750208](https://pubmed.ncbi.nlm.nih.gov/27750208/)  
  Endothelial-to-mesenchymal transition in the adapting leaflet.
- Grande-Allen KJ, Borowski AG, Troughton RW, et al. **Apparently normal mitral valves in patients with heart failure demonstrate biochemical and structural derangements: an extracellular matrix and echocardiographic study.** *J Am Coll Cardiol* 2005. [PMID 15629373](https://pubmed.ncbi.nlm.nih.gov/15629373/)  
  The 'normal' valve in heart failure is not normal.

## 7. Remodelling, growth laws, fibrosis and contractile loss

*The growth laws the model implements, and the reason the hypertrophy of mitral regurgitation does not protect the ventricle.*

- Grossman W, Jones D, McLaurin LP. **Wall stress and patterns of hypertrophy in the human left ventricle.** *J Clin Invest* 1975. [PMID 124746](https://pubmed.ncbi.nlm.nih.gov/124746/)  
  The stress-driven growth law the model implements, and the concentric/eccentric split.
- Nagatsu M, Zile MR, Tsutsui H, et al. **Native beta-adrenergic support for left ventricular dysfunction in experimental mitral regurgitation normalizes indexes of pump and contractile function.** *Circulation* 1994. [PMID 8313571](https://pubmed.ncbi.nlm.nih.gov/8313571/)  
  Contractile dysfunction in experimental MR.
- Spinale FG, Ishihra K, Zile M, et al. **Structural basis for changes in left ventricular function and geometry because of chronic mitral regurgitation and after correction of volume overload.** *J Thorac Cardiovasc Surg* 1993. [PMID 8246553](https://pubmed.ncbi.nlm.nih.gov/8246553/)  
  Structural basis of the ventricular change.
- Weber KT, Brilla CG. **Pathological hypertrophy and cardiac interstitium. Fibrosis and renin-angiotensin-aldosterone system.** *Circulation* 1991. [PMID 1828192](https://pubmed.ncbi.nlm.nih.gov/1828192/)  
  The RAAS-fibrosis link the model uses.
- Gaasch WH, Meyer TE. **Left ventricular response to mitral regurgitation: implications for management.** *Circulation* 2008. [PMID 19029478](https://pubmed.ncbi.nlm.nih.gov/19029478/)  
  The ventricular response and its clinical implications.
- Enriquez-Sarano M, Sundt TM 3rd. **Early surgery is recommended for mitral regurgitation.** *Circulation* 2010. [PMID 20159841](https://pubmed.ncbi.nlm.nih.gov/20159841/)  
  The case for early operation, i.e. against watchful waiting.
- Carabello BA, Nakano K, Corin W, et al. **Left ventricular function in experimental volume overload hypertrophy.** *Am J Physiol* 1989. [PMID 2523200](https://pubmed.ncbi.nlm.nih.gov/2523200/)  
  Volume-overload hypertrophy: thin-walled, not protective.
- Corporan D, Saadeh M, Yoldas A, et al. **Passive mechanical properties of the left ventricular myocardium and extracellular matrix in hearts with chronic volume overload from mitral regurgitation.** *Physiol Rep* 2022. [PMID 35871778](https://pubmed.ncbi.nlm.nih.gov/35871778/)  
  Matrix remodelling in chronic MR.
- Paulus WJ, Zile MR. **From Systemic Inflammation to Myocardial Fibrosis: The Heart Failure With Preserved Ejection Fraction Paradigm Revisited.** *Circ Res* 2021. [PMID 33983831](https://pubmed.ncbi.nlm.nih.gov/33983831/)  
  Inflammation-to-fibrosis signalling in the myocardium; the pathway behind the model's collagen state.
- Park SJ, Kim M, Son J, et al. **Long-Term Outcomes of Early Surgery Versus Conventional Treatment for Asymptomatic Severe Mitral Regurgitation: A Propensity Analysis.** *Circulation* 2025. [PMID 40799133](https://pubmed.ncbi.nlm.nih.gov/40799133/)  
  Long-term outcomes of early surgery versus watchful waiting in asymptomatic severe MR.
- Rosenhek R, Rader F, Klaar U, et al. **Outcome of watchful waiting in asymptomatic severe mitral regurgitation.** *Circulation* 2006. [PMID 16651470](https://pubmed.ncbi.nlm.nih.gov/16651470/)  
  What watchful waiting actually delivers.

## 8. Pulmonary vascular and right heart consequences

*The second barrier: fix the valve late enough and the lung has become the disease.*

- Barbieri A, Bursi F, Grigioni F, et al. **Prognostic and therapeutic implications of pulmonary hypertension complicating degenerative mitral regurgitation due to flail leaflet: a multicenter long-term international study.** *Eur Heart J* 2011. [PMID 20829213](https://pubmed.ncbi.nlm.nih.gov/20829213/)  
  Pulmonary hypertension as an independent barrier in MR.
- Gorter TM, van Veldhuisen DJ, Bauersachs J, et al. **Right heart dysfunction and failure in heart failure with preserved ejection fraction: mechanisms and management. Position statement on behalf of the Heart Failure Association of the European Society of Cardiology.** *Eur J Heart Fail* 2018. [PMID 29044932](https://pubmed.ncbi.nlm.nih.gov/29044932/)  
  Right heart dysfunction and RV-PA coupling in left-sided heart failure.
- Naeije R, Manes A. **The right ventricle in pulmonary arterial hypertension.** *Eur Respir Rev* 2014. [PMID 25445946](https://pubmed.ncbi.nlm.nih.gov/25445946/)  
  RV-PA coupling, which the model solves in closed form.
- Tello K, Wan J, Dalmer A, et al. **Validation of the Tricuspid Annular Plane Systolic Excursion/Systolic Pulmonary Artery Pressure Ratio for the Assessment of Right Ventricular-Arterial Coupling in Severe Pulmonary Hypertension.** *Circ Cardiovasc Imaging* 2019. [PMID 31500448](https://pubmed.ncbi.nlm.nih.gov/31500448/)  
  TAPSE/PASP as a clinical read-out of the RV coupling term.
- Ghoreishi M, Evans CF, DeFilippi CR, et al. **Pulmonary hypertension adversely affects short- and long-term survival after mitral valve operation for mitral regurgitation: implications for timing of surgery.** *J Thorac Cardiovasc Surg* 2011. [PMID 21962906](https://pubmed.ncbi.nlm.nih.gov/21962906/)  
  PH after operation: the second barrier, measured.
- Humbert M, Kovacs G, Hoeper MM, et al. **2022 ESC/ERS Guidelines for the diagnosis and treatment of pulmonary hypertension.** *Eur Heart J* 2022. [PMID 36017548](https://pubmed.ncbi.nlm.nih.gov/36017548/)  
  Pulmonary hypertension definitions, including the post-capillary and combined categories the model's PVR state maps onto.
- Bartko PE, Heitzinger G, Pavo N, et al. **Burden, treatment use, and outcome of secondary mitral regurgitation across the spectrum of heart failure: observational cohort study.** *BMJ* 2021. [PMID 34193442](https://pubmed.ncbi.nlm.nih.gov/34193442/)  
  Outcome burden of secondary MR.

## 9. Medical therapy

*Effect sizes the drug block is scaled against, and the vasodilator argument that follows from the leak and the aorta being in parallel.*

- McMurray JJ, Packer M, Desai AS, et al. **Angiotensin-neprilysin inhibition versus enalapril in heart failure.** *N Engl J Med* 2014. [PMID 25176015](https://pubmed.ncbi.nlm.nih.gov/25176015/)  
  PARADIGM-HF: the ARNI effect size the model's drug block is scaled against.
- Zannad F, McMurray JJ, Krum H, et al. **Eplerenone in patients with systolic heart failure and mild symptoms.** *N Engl J Med* 2011. [PMID 21073363](https://pubmed.ncbi.nlm.nih.gov/21073363/)  
  MRA effect size.
- McMurray JJV, Solomon SD, Inzucchi SE, et al. **Dapagliflozin in Patients with Heart Failure and Reduced Ejection Fraction.** *N Engl J Med* 2019. [PMID 31535829](https://pubmed.ncbi.nlm.nih.gov/31535829/)  
  SGLT2 inhibitor effect size and plasma volume effect.
- Packer M, Fowler MB, Roecker EB, et al. **Effect of carvedilol on the morbidity of patients with severe chronic heart failure: results of the carvedilol prospective randomized cumulative survival (COPERNICUS) study.** *Circulation* 2002. [PMID 12390947](https://pubmed.ncbi.nlm.nih.gov/12390947/)  
  Beta-blocker effect in severe HF.
- Ellison DH, Felker GM. **Diuretic Treatment in Heart Failure.** *N Engl J Med* 2017. [PMID 29141174](https://pubmed.ncbi.nlm.nih.gov/29141174/)  
  Diuretic pharmacology and its limits.
- Levine HJ, Gaasch WH. **Vasoactive drugs in chronic regurgitant lesions of the mitral and aortic valves.** *J Am Coll Cardiol* 1996. [PMID 8890799](https://pubmed.ncbi.nlm.nih.gov/8890799/)  
  The rationale, and the limits, of vasodilators in regurgitant lesions.
- Waagstein F, Bristow MR, Swedberg K, et al. **Beneficial effects of metoprolol in idiopathic dilated cardiomyopathy. Metoprolol in Dilated Cardiomyopathy (MDC) Trial Study Group.** *Lancet* 1993. [PMID 7902479](https://pubmed.ncbi.nlm.nih.gov/7902479/)  
  Beta-blockade recovers contractility.
- Sasayama S, Ohyagi A, Lee JD, et al. **Effect of the vasodilator therapy in regurgitant valvular disease.** *Jpn Circ J* 1982. [PMID 7087161](https://pubmed.ncbi.nlm.nih.gov/7087161/)  
  Vasodilators reduce regurgitant fraction acutely.
- Murphy SP, Ibrahim NE, Januzzi JL Jr. **Heart Failure With Reduced Ejection Fraction: A Review.** *JAMA* 2020. [PMID 32749493](https://pubmed.ncbi.nlm.nih.gov/32749493/)  
  Contemporary HFrEF pharmacotherapy, including rate control, as a reference for the drug block.

## 10. Devices, surgery, recurrence and guidelines

*Which arm of the feedback loop each procedure actually cuts.*

- Feldman T, Kar S, Elmariah S, et al. **Randomized Comparison of Percutaneous Repair and Surgery for Mitral Regurgitation: 5-Year Results of EVEREST II.** *J Am Coll Cardiol* 2015. [PMID 26718672](https://pubmed.ncbi.nlm.nih.gov/26718672/)  
  EVEREST II: the device compared with surgery in primary MR.
- Acker MA, Parides MK, Perrault LP, et al. **Mitral-valve repair versus replacement for severe ischemic mitral regurgitation.** *N Engl J Med* 2014. [PMID 24245543](https://pubmed.ncbi.nlm.nih.gov/24245543/)  
  Repair versus replacement in ischaemic MR: recurrence is the story.
- Goldstein D, Moskowitz AJ, Gelijns AC, et al. **Two-Year Outcomes of Surgical Treatment of Severe Ischemic Mitral Regurgitation.** *N Engl J Med* 2016. [PMID 26550689](https://pubmed.ncbi.nlm.nih.gov/26550689/)  
  Two-year recurrence after annuloplasty.
- Kron IL, Hung J, Overbey JR, et al. **Predicting recurrent mitral regurgitation after mitral valve repair for severe ischemic mitral regurgitation.** *J Thorac Cardiovasc Surg* 2015. [PMID 25500293](https://pubmed.ncbi.nlm.nih.gov/25500293/)  
  What predicts recurrence: the tethering arm the ring does not cut.
- Cleland JG, Daubert JC, Erdmann E, et al. **Longer-term effects of cardiac resynchronization therapy on mortality in heart failure [the CArdiac REsynchronization-Heart Failure (CARE-HF) trial extension phase].** *Eur Heart J* 2006. [PMID 16782715](https://pubmed.ncbi.nlm.nih.gov/16782715/)  
  CRT effect size.
- Otto CM, Nishimura RA, Bonow RO, et al. **2020 ACC/AHA Guideline for the Management of Patients With Valvular Heart Disease: Executive Summary: A Report of the American College of Cardiology/American Heart Association Joint Committee on Clinical Practice Guidelines.** *Circulation* 2021. [PMID 33332149](https://pubmed.ncbi.nlm.nih.gov/33332149/)  
  The guideline thresholds the model's EF/LVESD results speak to.
- Vahanian A, Beyersdorf F, Praz F, et al. **2021 ESC/EACTS Guidelines for the management of valvular heart disease.** *Eur Heart J* 2022. [PMID 34453165](https://pubmed.ncbi.nlm.nih.gov/34453165/)  
  European thresholds and the timing question.
- Hung J, Papakostas L, Tahta SA, et al. **Mechanism of recurrent ischemic mitral regurgitation after annuloplasty: continued LV remodeling as a moving target.** *Circulation* 2004. [PMID 15364844](https://pubmed.ncbi.nlm.nih.gov/15364844/)  
  Recurrence as continued ventricular remodelling.
- Michler RE, Smith PK, Parides MK, et al. **Two-Year Outcomes of Surgical Treatment of Moderate Ischemic Mitral Regurgitation.** *N Engl J Med* 2016. [PMID 27040451](https://pubmed.ncbi.nlm.nih.gov/27040451/)  
  Adding an annuloplasty ring to revascularisation in moderate ischaemic MR.
- Breithardt OA, Sinha AM, Schwammenthal E, et al. **Acute effects of cardiac resynchronization therapy on functional mitral regurgitation in advanced systolic heart failure.** *J Am Coll Cardiol* 2003. [PMID 12628720](https://pubmed.ncbi.nlm.nih.gov/12628720/)  
  CRT reduces functional MR by restoring closing force.
- Du Y, Han H, Zhang T, et al. **Prognosis of Elevated Mitral Valve Pressure Gradient After Transcatheter Edge-to-Edge Repair: Systematic Review and Meta-Analysis.** *Curr Probl Cardiol* 2024. [PMID 37778430](https://pubmed.ncbi.nlm.nih.gov/37778430/)  
  The iatrogenic gradient as a device trade-off.
- Fang JX, Giustino G, So KC, et al. **Management Approach for Residual and Recurrent Mitral Regurgitation After Transcatheter Edge-to-Edge Repair.** *JACC Cardiovasc Interv* 2026. [PMID 42120115](https://pubmed.ncbi.nlm.nih.gov/42120115/)  
  Residual MR after TEER.
- Westaby S. **Preservation of left ventricular function in mitral valve surgery.** *Heart* 1996. [PMID 8705754](https://pubmed.ncbi.nlm.nih.gov/8705754/)  
  Chordal preservation and the ventricular scaffold.

## 11. Atrial fibrillation and atrial myopathy

*The second arm of the vortex, and the reason a long-compensated patient decompensates with an unchanged valve.*

- Marsan NA, Maffessanti F, Tamborini G, et al. **Left atrial reverse remodeling and functional improvement after mitral valve repair in degenerative mitral regurgitation: a real-time 3-dimensional echocardiography study.** *Am Heart J* 2011. [PMID 21315214](https://pubmed.ncbi.nlm.nih.gov/21315214/)  
  Atrial reverse remodelling after repair.
- Bisbal F, Baranchuk A, Braunwald E, et al. **Atrial Failure as a Clinical Entity: JACC Review Topic of the Week.** *J Am Coll Cardiol* 2020. [PMID 31948652](https://pubmed.ncbi.nlm.nih.gov/31948652/)  
  Atrial myopathy as an entity in its own right.
- Tribouilloy C, Bohbot Y, Essayagh B, et al. **Prognostic implications of functional tricuspid regurgitation in asymptomatic degenerative mitral regurgitation.** *ESC Heart Fail* 2025. [PMID 40122615](https://pubmed.ncbi.nlm.nih.gov/40122615/)  
  Right-sided consequences of degenerative mitral regurgitation in asymptomatic patients.
- Stassen J, Namazi F, van der Bijl P, et al. **Left Atrial Reservoir Function and Outcomes in Secondary Mitral Regurgitation.** *J Am Soc Echocardiogr* 2022. [PMID 35074443](https://pubmed.ncbi.nlm.nih.gov/35074443/)  
  Atrial reservoir function.
- Borger MA, Mansour MC, Levine RA. **Atrial Fibrillation and Mitral Valve Prolapse: Time to Intervene?.** *J Am Coll Cardiol* 2019. [PMID 30678756](https://pubmed.ncbi.nlm.nih.gov/30678756/)  
  AF and device outcomes.

## 12. Modelling methodology

*Closed-loop lumped-parameter circulation, elastance mechanics, stress-driven growth, and the simulation engine.*

- Elmokadem A, Riggs MM, Baron KT. **Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.** *CPT Pharmacometrics Syst Pharmacol* 2019. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)  
  mrgsolve in QSP practice.
- Burkhoff D, Tyberg JV. **Why does pulmonary venous pressure rise after onset of LV dysfunction: a theoretical analysis.** *Am J Physiol* 1993. [PMID 8238596](https://pubmed.ncbi.nlm.nih.gov/8238596/)  
  A closed-loop lumped-parameter analysis of exactly the question this model asks.
- Santamore WP, Burkhoff D. **Hemodynamic consequences of ventricular interaction as assessed by model analysis.** *Am J Physiol* 1991. [PMID 1992793](https://pubmed.ncbi.nlm.nih.gov/1992793/)  
  Ventricular interaction in a lumped model.
- Guyton AC, Coleman TG, Granger HJ. **Circulation: overall regulation.** *Annu Rev Physiol* 1972. [PMID 4334846](https://pubmed.ncbi.nlm.nih.gov/4334846/)  
  The volume/pressure-natriuresis closure the model uses.
- Witzenburg CM, Holmes JW. **Predicting the Time Course of Ventricular Dilation and Thickening Using a Rapid Compartmental Model.** *J Cardiovasc Transl Res* 2018. [PMID 29550925](https://pubmed.ncbi.nlm.nih.gov/29550925/)  
  A compartmental growth model of exactly this kind.
- Arts T, Delhaas T, Bovendeerd P, et al. **Adaptation to mechanical load determines shape and properties of heart and circulation: the CircAdapt model.** *Am J Physiol Heart Circ Physiol* 2005. [PMID 15550528](https://pubmed.ncbi.nlm.nih.gov/15550528/)  
  Load-adaptation modelling of chamber geometry.
- Smith BW, Chase JG, Nokes RI, et al. **Minimal haemodynamic system model including ventricular interaction and valve dynamics.** *Med Eng Phys* 2004. [PMID 15036180](https://pubmed.ncbi.nlm.nih.gov/15036180/)  
  A minimal closed-loop model with valve dynamics.
- Sharifi H, Lee LC, Campbell KS, et al. **A multiscale finite element model of left ventricular mechanics incorporating baroreflex regulation.** *Comput Biol Med* 2024. [PMID 37984204](https://pubmed.ncbi.nlm.nih.gov/37984204/)  
  Multiscale finite-element ventricular mechanics with baroreflex control; a higher-resolution counterpart to the lumped model used here.

---

## Reproducing this bibliography

Each entry was resolved with `esearch` (relevance-sorted, top hit) and
described with `esummary` against `eutils.ncbi.nlm.nih.gov`. The harvest
and rendering scripts are `refs.py`, `refs2.py`, `refs3.py` and
`mkrefs.py`; the raw NCBI payloads are kept in `refs_raw*.json`. Any PMID
in this file can be checked directly:

```bash
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"\
     "?db=pubmed&retmode=json&id=30280640" | python3 -m json.tool
```

## A note on what is calibrated and what is predicted

Only two numbers in the whole model are fitted to outcome data, and both
come from the **COAPT control arm**: the baseline hazard for heart-failure
hospitalisation and the baseline hazard for death. Every other trial
quantity referenced above — the COAPT device arm, both MITRA-FR arms, and
the baseline echocardiographic and haemodynamic values the virtual patients
reproduce — is a prediction, and the ones the model gets wrong are reported
as misses in `README.md` rather than absorbed by refitting.

