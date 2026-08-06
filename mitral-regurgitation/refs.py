#!/usr/bin/env python3
"""Harvest the reference list for the mitral regurgitation QSP model.

Every entry is resolved against PubMed: the query pins the intended paper, and
the title / journal / year / PMID written into the markdown come from NCBI, not
from memory.  Anything that fails to resolve is reported rather than invented.
"""
import json
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def get(url):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=40) as r:
                return r.read().decode()
        except Exception:
            time.sleep(1.5 * (attempt + 1))
    return None


def esearch(q):
    u = EUTILS + "esearch.fcgi?db=pubmed&retmax=1&sort=relevance&retmode=json&term=" \
        + urllib.parse.quote(q)
    t = get(u)
    if not t:
        return None
    try:
        ids = json.loads(t)["esearchresult"]["idlist"]
        return ids[0] if ids else None
    except Exception:
        return None


def esummary(pmids):
    out = {}
    for i in range(0, len(pmids), 150):
        chunk = pmids[i:i + 150]
        u = EUTILS + "esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(chunk)
        t = get(u)
        if not t:
            continue
        try:
            res = json.loads(t)["result"]
        except Exception:
            continue
        for p in chunk:
            if p in res:
                out[p] = res[p]
        time.sleep(0.34)
    return out


# (section, query, editorial note)
Q = [
 # ---- 1. the two pivotal trials and the proportionality argument -------------
 ("1", "Stone Lindenfeld Abraham Transcatheter mitral-valve repair patients heart failure COAPT 2018 N Engl J Med", "COAPT. The trial the model is calibrated against: HF hospitalisation 35.8 vs 67.9 per patient-year, death 29.1% vs 46.1% at 24 months."),
 ("1", "Obadia Messika-Zeitoun Percutaneous repair medical treatment secondary mitral regurgitation MITRA-FR 2018", "MITRA-FR. Same device, opposite result: death or HF hospitalisation 54.6% vs 51.3% at 12 months."),
 ("1", "Grayburn Sannino Packer proportionate disproportionate functional mitral regurgitation new conceptual framework", "The proportionality framework. This is the paper the model turns into arithmetic."),
 ("1", "Packer Grayburn Contrasting effects of pharmacological transcatheter interventions functional mitral regurgitation", "Why drug and device effects on functional MR are not interchangeable."),
 ("1", "Packer Grayburn New evidence supporting a novel conceptual framework for distinguishing proportionate disproportionate functional mitral regurgitation", "The framework stated as a testable hypothesis."),
 ("1", "Mack Lindenfeld Abraham COAPT 3-year outcomes transcatheter mitral valve repair heart failure", "Durability of the COAPT effect."),
 ("1", "Iung Armoiry Vahanian percutaneous repair mitral regurgitation MITRA-FR 24 months outcomes", "MITRA-FR at 24 months: still no benefit."),
 ("1", "Pibarot Delgado Bax MITRA-FR versus COAPT lessons echocardiography", "Echocardiographic reconciliation attempt between the two trials."),
 ("1", "Stone Abraham Lindenfeld five-year follow-up COAPT transcatheter mitral valve repair", "Five-year COAPT follow-up."),
 ("1", "Anker Friede von Bardeleben RESHAPE-HF2 percutaneous repair functional mitral regurgitation heart failure", "A third randomised trial in the same space."),

 # ---- 2. quantifying the lesion: EROA, PISA and its geometry -----------------
 ("2", "Enriquez-Sarano Seward Bailey Tajik effective regurgitant orifice area noninvasive Doppler quantitation", "The effective regurgitant orifice area as a clinical quantity."),
 ("2", "Enriquez-Sarano Avierinos Messika-Zeitoun quantitative determinants clinical outcome mitral regurgitation", "EROA predicts outcome in primary MR; the origin of the 0.4 cm2 threshold."),
 ("2", "Bargiggia Tronconi Sahn new method quantitation mitral regurgitation based on proximal flow convergence", "PISA: the proximal isovelocity surface area method."),
 ("2", "Utsunomiya Ogawa Doshi proximal isovelocity surface area method mitral regurgitation in vitro validation", "PISA validated in vitro, where the orifice really is circular."),
 ("2", "Hopmeyer He Yoganathan Levine three-dimensional geometry proximal flow convergence non-circular orifice", "Why a non-hemispheric convergence zone breaks the PISA assumption."),
 ("2", "Yosefy Levine Solis Vaturi Handschumacher proximal flow convergence region elongated shape functional mitral regurgitation PISA underestimates", "The functional orifice is elongated, and PISA mis-measures it."),
 ("2", "Buck Plicht Kahlert Erbel Hunold Corrective single-plane PISA method quantification eccentric mitral regurgitation", "Correcting PISA for a non-circular orifice."),
 ("2", "Zoghbi Adams Bonow recommendations noninvasive evaluation native valvular regurgitation ASE 2017", "The guideline grading scheme, including its own caveats about method disagreement."),
 ("2", "Uretsky Gillam Lang discordance echocardiography cardiac magnetic resonance quantification mitral regurgitation", "Echo and CMR disagree; the model treats this as information, not noise."),
 ("2", "Grayburn Weissman Zamorano quantitation mitral regurgitation Circulation review", "Review of what each method actually measures."),
 ("2", "Thomas Foster Hagege Vandervoort Levine proximal isovelocity surface area method quantification mitral regurgitation limitations", "Limitations of PISA in practice."),
 ("2", "Dujardin Enriquez-Sarano Bailey Seward Grading of mitral regurgitation quantitative Doppler comparison", "Quantitative Doppler grading, and where the cut-points came from."),
 ("2", "Biner Rafique Rafii Tolstrup Noorani Shiota Reproducibility proximal isovelocity surface area vena contracta mitral regurgitation", "Reproducibility of the measurements the trials relied on."),
 ("2", "Militaru Kinnaird Lancellotti Choure Piazza validation Doppler quantification mitral regurgitation cardiac magnetic resonance", "Cross-modality validation."),

 # ---- 3. the acute / chronic distinction and atrial compliance ---------------
 ("3", "Braunwald Awe syndrome severe mitral regurgitation with normal left atrial pressure", "The observation the model reproduces: severe MR can coexist with a normal atrial pressure."),
 ("3", "Braunwald mitral regurgitation physiologic clinical considerations", "The physiology of the acute-versus-chronic distinction, stated early and clearly."),
 ("3", "Pomerantz Nunez Bhatt giant v wave pulmonary capillary wedge pressure mitral regurgitation", "The giant v wave: the atrial denominator made visible."),
 ("3", "Braunwald Frahm studies on Starling law of the heart left atrial left ventricular pressures man", "Atrial pressure-volume behaviour in man."),
 ("3", "Barbier Solomon Schiller Glantz left atrial relaxation chamber stiffness in humans pressure volume", "Left atrial chamber stiffness measured directly."),
 ("3", "Stefanadis Dernellis Toutouzas clinical appraisal left atrial function pressure volume relation", "The left atrial pressure-volume relation, which is the model's denominator 1."),
 ("3", "Nishimura Schaff Shub Gersh Edwards Tajik papillary muscle rupture complicating acute myocardial infarction", "Acute severe MR from papillary rupture: analysis of the clinical syndrome."),
 ("3", "Thompson Balser Wolfe acute mitral regurgitation haemodynamics management cardiogenic shock", "Management of acute severe MR."),
 ("3", "Rosario Kern Vasodilator therapy nitroprusside acute mitral regurgitation forward stroke volume", "Afterload reduction in acute MR: the impedance argument."),
 ("3", "Yoran Yellin Becker Gabbay Frater Sonnenblick dynamic aspects of acute mitral regurgitation effects of ventricular volume pressure atrial compliance", "Atrial compliance as an explicit determinant of the presentation."),
 ("3", "Hodges Kirkorian Sohmer left atrial compliance functional mitral regurgitation pulmonary hypertension", "Atrial compliance and the pulmonary consequences."),
 ("3", "Grose Strain Yellin Cohen Adaptation of left ventricle to chronic mitral regurgitation left ventricular volume pressure", "Chronic adaptation of the ventricle."),

 # ---- 4. ejection fraction, hidden contractility, operative thresholds -------
 ("4", "Enriquez-Sarano Tajik Schaff Orszulak Bailey Frye echocardiographic prediction survival after surgical correction organic mitral regurgitation", "The origin of the EF 60% operative threshold."),
 ("4", "Tribouilloy Grigioni Avierinos Barbieri Survival implication of left ventricular end-systolic diameter mitral regurgitation surgery", "The LV end-systolic diameter trigger."),
 ("4", "Starling Kirklin Bolling Effects of valve surgery on left ventricular contractile function patients with long-term mitral regurgitation", "Contractile function measured before and after operation."),
 ("4", "Wisenbaugh Skudicky Sareli Prediction of outcome after valve replacement for rheumatic mitral regurgitation era of chordal preservation", "Post-operative outcome and the ventricle's hidden state."),
 ("4", "Carabello Nolan McGuire Assessment of preoperative left ventricular function in patients with mitral regurgitation value of the end-systolic wall stress-ejection fraction relation", "Why EF alone misleads in MR: the wall-stress correction."),
 ("4", "Carabello The current therapy for mitral regurgitation JACC", "Carabello's synthesis of the loading problem."),
 ("4", "Zile Gaasch Carroll Levine Chronic mitral regurgitation predictive value of preoperative echocardiographic indexes of left ventricular function wall stress", "Preoperative indices that survive the afterload problem."),
 ("4", "Suga Sagawa Shoukas Load independence of the instantaneous pressure-volume ratio", "The elastance concept the model's ventricle is built on."),
 ("4", "Sunagawa Maughan Burkhoff Sagawa Left ventricular interaction with arterial load studied in isolated canine ventricle", "Ventricular-arterial coupling: where E_a comes from."),
 ("4", "Kelly Ting Yang Liu Maughan Chang Effective arterial elastance as index of arterial vascular load in humans", "Effective arterial elastance in humans."),
 ("4", "Senzaki Chen Kass Single-beat estimation of end-systolic pressure-volume relation in humans", "Single-beat estimation of E_es, and its normal range."),
 ("4", "Chen Fetics Nevo Rochitte Chiou Noninvasive single-beat determination of left ventricular end-systolic elastance in humans", "Non-invasive E_es, used for the model's baseline value."),
 ("4", "Klotz Hay Dickstein Yi Wang Maurer Development and validation of empirical formula to estimate normalized end-diastolic pressure volume relationship", "The EDPVR form the model uses."),
 ("4", "Burkhoff Mirsky Suga Assessment of systolic and diastolic ventricular properties via pressure-volume analysis guide for clinical translational research", "Pressure-volume analysis as used here."),

 # ---- 5. secondary MR mechanism: tethering, annulus, geometry ---------------
 ("5", "Otsuji Handschumacher Schwammenthal Levine Insights from three-dimensional echocardiography into the mechanism of functional mitral regurgitation direct in vivo demonstration of altered leaflet tethering geometry", "The tethering mechanism demonstrated in vivo."),
 ("5", "Yiu Enriquez-Sarano Tribouilloy Seward Tajik Determinants of the degree of functional mitral regurgitation in patients with systolic left ventricular dysfunction quantitative clinical study", "What actually determines functional MR severity."),
 ("5", "Kwan Shiota Agler Popovic Qin Geometric differences of the mitral apparatus between ischemic and dilated cardiomyopathy with significant mitral regurgitation real-time three-dimensional", "Ischaemic versus dilated geometry."),
 ("5", "Levine Schwammenthal Ischemic mitral regurgitation on the threshold of a solution from paradoxes to unifying concepts", "The unifying geometric account of ischaemic MR."),
 ("5", "Komeda Glasson Bolger Daughters Ingels Miller Geometric determinants of ischemic mitral regurgitation", "Geometric determinants, quantified."),
 ("5", "Schwammenthal Chen Ramires Guerrero Levine Dynamics of mitral regurgitant flow and orifice area physiologic application of the proximal flow convergence method", "The orifice is dynamic within systole."),
 ("5", "Hung Otsuji Handschumacher Schwammenthal Levine Mechanism of dynamic regurgitant orifice area variation in functional mitral regurgitation", "Why the functional orifice changes with load."),
 ("5", "Silbiger Mechanistic insights into atrial functional mitral regurgitation far more complicated than just left atrial remodeling", "Atrial functional MR: the ventricle need not be involved at all."),
 ("5", "Deferm Bertrand Verbrugge Verhaert Atrial functional mitral regurgitation JACC state of the art review", "Atrial functional MR as a distinct entity."),
 ("5", "Ormiston Shah Tyberg Bell Size and motion of the mitral valve annulus in man normal", "Normal annular size and motion, used for the model's annular area."),
 ("5", "Levine Handschumacher Sanfilippo Hagege Harrigan Three-dimensional echocardiographic reconstruction of the mitral valve with implications for the diagnosis of mitral valve prolapse saddle", "The annular saddle shape."),
 ("5", "Kaplan Bashein Sheehan Legget Munt Three-dimensional echocardiographic assessment of annular shape changes in the normal and regurgitant mitral valve", "Annular shape change with regurgitation."),
 ("5", "Otsuji Kumanohoso Yoshifuku Matsukida Koriyama Isolated annular dilation does not usually cause important functional mitral regurgitation comparison between patients with lone atrial fibrillation and those with cardiomyopathy", "A key negative result the model must reproduce: annular dilation alone is often not enough."),
 ("5", "He Fontaine Sauren Grimm Yoganathan In vitro dynamic strain behavior of the mitral valve posterior leaflet", "Leaflet mechanics."),
 ("5", "Gorman Gupta Streicher Gorman Jackson Dynamic three-dimensional imaging of the mitral valve and left ventricle by rapid sonomicrometry array localization", "Dynamic mitral apparatus imaging."),

 # ---- 6. leaflet plasticity -------------------------------------------------
 ("6", "Dal-Bianco Aikawa Bischoff Guerrero Handschumacher Active adaptation of the tethered mitral valve insights into a compensatory mechanism for functional mitral regurgitation", "Leaflets actively grow under tethering: the model's leaflet-supply state."),
 ("6", "Dal-Bianco Aikawa Bischoff Guerrero Sullivan Myocardial infarction alters adaptation of the tethered mitral valve", "Infarction alters the adaptation, which is why fast dilation is worse."),
 ("6", "Beaudoin Dal-Bianco Aikawa Bischoff Guerrero Mitral leaflet changes following myocardial infarction clinical evidence for maladaptive valvular remodeling", "Clinical evidence of maladaptive leaflet remodelling."),
 ("6", "Chaput Handschumacher Tournoux Hua Guerrero Vlahakes Mitral leaflet adaptation to ventricular remodeling occurrence and adequacy in patients with functional mitral regurgitation", "Adaptation occurs but is inadequate: the phi<1 assumption."),
 ("6", "Rausch Tibayan Ingels Miller Kuhl Mechanics of the mitral annulus in chronic ischemic cardiomyopathy", "Annular mechanics in chronic ischaemic disease."),
 ("6", "Grande-Allen Borowski Troughton Thomas Vesely Apparently normal mitral valves in patients with heart failure demonstrate biochemical and structural derangement", "The 'normal' valve in heart failure is not normal."),
 ("6", "Bischoff Casanovas Wylie-Sears Kim Bartko Guerrero CD45 expression in mitral valve endothelial cells after myocardial infarction", "Endothelial-to-mesenchymal transition in the adapting leaflet."),

 # ---- 7. remodelling, growth, fibrosis, contractility ----------------------
 ("7", "Grossman Jones McLaurin Wall stress and patterns of hypertrophy in the human left ventricle", "The stress-driven growth law the model implements, and the concentric/eccentric split."),
 ("7", "Carabello Zile Tanaka Cooper Volume overload hypertrophy in a closed-chest model of mitral regurgitation", "Volume-overload hypertrophy: thin-walled, not protective."),
 ("7", "Nagatsu Zile Tsutsui Schmid Carabello Native beta-adrenergic support for left ventricular dysfunction in experimental mitral regurgitation normalizes systolic function", "Contractile dysfunction in experimental MR."),
 ("7", "Spinale Ishihra Zile DeFryte Crawford Carabello Structural basis for changes in left ventricular function and geometry because of chronic mitral regurgitation", "Structural basis of the ventricular change."),
 ("7", "Corporan Kono Onohara Padala Chronic volume overload mitral regurgitation extracellular matrix remodeling", "Matrix remodelling in chronic MR."),
 ("7", "Kramann Hasenfuss Seidler B-type natriuretic peptide and left ventricular remodeling chronic mitral regurgitation", "Biomarkers of remodelling in MR."),
 ("7", "Weber Brilla Pathological hypertrophy and cardiac interstitium fibrosis and renin-angiotensin-aldosterone system", "The RAAS-fibrosis link the model uses."),
 ("7", "Zile Baicu Ikonomidis Stroud Nietert Bradshaw Myocardial stiffness in patients with heart failure and a preserved ejection fraction contributions of collagen and titin", "Collagen and titin as the two stiffness terms."),
 ("7", "Gaasch Meyer Left ventricular response to mitral regurgitation implications for management", "The ventricular response and its clinical implications."),
 ("7", "Enriquez-Sarano Sundt Early surgery is recommended for mitral regurgitation", "The case for early operation, i.e. against watchful waiting."),
 ("7", "Kang Kim Yun Song Choi Chung Comparison of early surgery versus conventional treatment in asymptomatic severe mitral regurgitation", "Early surgery versus watchful waiting, randomised-adjacent evidence."),
 ("7", "Rosenhek Rader Klaar Gabriel Krejc Binder Outcome of watchful waiting in asymptomatic severe mitral regurgitation", "What watchful waiting actually delivers."),

 # ---- 8. pulmonary vascular and right heart consequences --------------------
 ("8", "Barbieri Bursi Grigioni Tribouilloy Avierinos Prognostic and therapeutic implications of pulmonary hypertension complicating degenerative mitral regurgitation", "Pulmonary hypertension as an independent barrier in MR."),
 ("8", "Ghoreishi Evans DeFilippi Hobbs Young Griffith Pulmonary hypertension adversely affects short- and long-term survival after mitral valve operation", "PH after operation: the second barrier, measured."),
 ("8", "Guazzi Borlaug Pulmonary hypertension due to left heart disease", "Post-capillary pulmonary hypertension and its reversible/fixed components."),
 ("8", "Vachiery Tedford Rosenkranz Palazzini Lang Guazzi Pulmonary hypertension due to left heart disease definition classification", "Definitions used for the combined post-capillary component."),
 ("8", "Naeije Manes The right ventricle in pulmonary arterial hypertension", "RV-PA coupling, which the model solves in closed form."),
 ("8", "Tello Wan Dalmer Vanderpool Ghofrani Naeije Validation of the tricuspid annular plane systolic excursion pulmonary artery systolic pressure ratio for the assessment of right ventricular arterial coupling", "TAPSE/PASP as a clinical read-out of the RV coupling term."),
 ("8", "Bartko Hulsmann Hung Frangieh Levine Praz Burden of secondary mitral regurgitation clinical outcome", "Outcome burden of secondary MR including right heart involvement."),

 # ---- 9. medical therapy ---------------------------------------------------
 ("9", "McMurray Packer Desai Gong Lefkowitz Rizkala Angiotensin-neprilysin inhibition versus enalapril in heart failure PARADIGM-HF", "PARADIGM-HF: the ARNI effect size the model's drug block is scaled against."),
 ("9", "Kang Park Sun Lee Kim Sacubitril valsartan versus valsartan in secondary functional mitral regurgitation PRIME", "PRIME: ARNI reduces functional MR, the model's afterload-plus-remodelling prediction."),
 ("9", "Zannad McMurray Krum van Veldhuisen Swedberg Shi Eplerenone in patients with systolic heart failure and mild symptoms EMPHASIS-HF", "MRA effect size."),
 ("9", "McMurray Solomon Inzucchi Kober Kosiborod Martinez Dapagliflozin in patients with heart failure and reduced ejection fraction DAPA-HF", "SGLT2 inhibitor effect size and plasma volume effect."),
 ("9", "Packer Coats Fowler Katus Krum Mohacsi Effect of carvedilol on survival in severe chronic heart failure COPERNICUS", "Beta-blocker effect in severe HF."),
 ("9", "Waagstein Bristow Swedberg Camerini Fowler Beneficial effects of metoprolol in idiopathic dilated cardiomyopathy reverse remodeling", "Beta-blockade recovers contractility: the model's slow Ees recovery term."),
 ("9", "Hallberg Brater Wilkin Loop diuretic pharmacokinetics pharmacodynamics furosemide braking phenomenon", "The braking phenomenon the model implements."),
 ("9", "Ellison Felker Diuretic treatment in heart failure New England Journal of Medicine review", "Diuretic pharmacology and its limits."),
 ("9", "Ayoub Prabhu Chatterjee Nitroprusside afterload reduction mitral regurgitation regurgitant fraction reduction haemodynamic", "Vasodilators reduce regurgitant fraction acutely."),
 ("9", "Levine Gaasch Vasoactive drugs in chronic regurgitant lesions of the mitral and aortic valves", "The rationale, and the limits, of vasodilators in regurgitant lesions."),
 ("9", "Comin-Colet Sanchez-Corral Manito Gomez-Otero Ivabradine heart rate reduction heart failure mitral regurgitation", "Heart-rate reduction, which the model shows cuts both ways."),

 # ---- 10. devices, surgery and recurrence ---------------------------------
 ("10", "Feldman Foster Glower Kaul Rinaldi Percutaneous repair or surgery for mitral regurgitation EVEREST II", "EVEREST II: the device compared with surgery in primary MR."),
 ("10", "Acker Parides Perrault Argenziano Ascheim Gelijns Mitral-valve repair versus replacement for severe ischemic mitral regurgitation", "Repair versus replacement in ischaemic MR: recurrence is the story."),
 ("10", "Goldstein Moskowitz Gelijns Ailawadi Parides Perrault Two-year outcomes of surgical treatment of severe ischemic mitral regurgitation", "Two-year recurrence after annuloplasty."),
 ("10", "Kron Hung Overbey Bouchard Gelijns Ailawadi Predicting recurrent mitral regurgitation after mitral valve repair for severe ischemic mitral regurgitation", "What predicts recurrence: the tethering arm the ring does not cut."),
 ("10", "Hung Papakostas Tahta Parranto Hardy Duran Mechanism of recurrent ischemic mitral regurgitation after annuloplasty continued LV remodeling as a moving target", "Recurrence as continued ventricular remodelling."),
 ("10", "Michler Smith Ailawadi Gammie Chikwe Two-year outcomes of surgical treatment of moderate ischemic mitral regurgitation CABG", "Adding a ring to revascularisation."),
 ("10", "Cleland Daubert Erdmann Freemantle Gras Kappenberger The effect of cardiac resynchronization on morbidity and mortality in heart failure CARE-HF", "CRT effect size."),
 ("10", "Breithardt Sinha Schwammenthal Bidaoui Markiewicz Levine Acute effects of cardiac resynchronization therapy on functional mitral regurgitation in advanced systolic heart failure", "CRT reduces functional MR by restoring closing force."),
 ("10", "Neuss Schueler Bauer Senges Nickenig Mitral valve gradient after transcatheter edge-to-edge repair outcome iatrogenic stenosis", "The iatrogenic gradient the model tracks as a device trade-off."),
 ("10", "Yoon Bhatt Koseki Kaewkes Patel Predictors of clinical outcomes transcatheter edge-to-edge repair residual mitral regurgitation", "Residual MR after TEER."),
 ("10", "Yun Sintek Fletcher Pfeffer Kerendi Miller Time course of left ventricular remodeling after mitral valve repair chordal preservation", "Chordal preservation and the ventricular scaffold."),
 ("10", "Rogers Bolling The tricuspid valve mitral annuloplasty ring undersizing recurrent regurgitation", "Undersized annuloplasty and its consequences."),
 ("10", "Nishimura Otto Bonow Carabello Erwin Guidelines management patients valvular heart disease ACC AHA 2020", "The guideline thresholds the model's EF/LVESD results speak to."),
 ("10", "Vahanian Beyersdorf Praz Milojevic Baldus 2021 ESC EACTS Guidelines for the management of valvular heart disease", "European thresholds and the timing question."),

 # ---- 11. atrial fibrillation, atrial myopathy and its consequences --------
 ("11", "Grigioni Avierinos Ling Scott Bailey Tajik Atrial fibrillation complicating the course of degenerative mitral regurgitation determinants and long-term outcome", "AF in MR: determinants and prognosis."),
 ("11", "Marsan Maffessanti Tamborini Gripari Caiani Left atrial reverse remodeling and functional improvement after mitral valve repair", "Atrial reverse remodelling after repair."),
 ("11", "Kihara Sasayama Suzuki Left atrial function mitral regurgitation reservoir conduit", "Atrial reservoir function."),
 ("11", "Bisbal Baranchuk Braunwald Bayes de Luna Bayes-Genis Atrial failure as a clinical entity JACC review", "Atrial myopathy as an entity in its own right."),
 ("11", "Rottlaender Motloch Schmidt Reda Larbig Wolny Clinical impact of atrial fibrillation in patients undergoing percutaneous mitral valve repair", "AF and device outcomes."),

 # ---- 12. modelling methodology ------------------------------------------
 ("12", "Baraldi Kaimal mrgsolve simulation from ordinary differential equation based models R", "The simulation engine the R model targets."),
 ("12", "Elmokadem Riggs Baron Quantitative systems pharmacology and physiologically-based pharmacokinetic modeling with mrgsolve an open-source R package", "mrgsolve in QSP practice."),
 ("12", "Burkhoff Tyberg Why does pulmonary venous pressure rise after onset of LV dysfunction a theoretical analysis", "A closed-loop lumped-parameter analysis of exactly the question this model asks."),
 ("12", "Santamore Burkhoff Hemodynamic consequences of ventricular interaction as assessed by model analysis", "Ventricular interaction in a lumped model."),
 ("12", "Smith Chase Nokes Shaw Wake Minimal haemodynamic system model including ventricular interaction and valve dynamics", "A minimal closed-loop cardiovascular model with valve dynamics."),
 ("12", "Guyton Coleman Granger Circulation overall regulation Annual Review of Physiology", "The volume/pressure-natriuresis closure the model uses."),
 ("12", "Witzenburg Holmes Predicting the time course of ventricular dilation and thickening using a rapid compartmental model", "A compartmental growth model of exactly this kind."),
 ("12", "Arts Delhaas Bovendeerd Verbeek Prinzen Adaptation to mechanical load determines shape and properties of heart and circulation the CircAdapt model", "Load-adaptation modelling of chamber geometry."),
 ("12", "Lee Wall Genet Hinson Guccione Bioengineering analysis of mitral regurgitation finite element", "A finite-element counterpart to the lumped approach used here."),
]


def main():
    print("resolving %d queries..." % len(Q))
    rows = []
    for sec, q, note in Q:
        pmid = esearch(q)
        rows.append((sec, q, note, pmid))
        time.sleep(0.34)
    got = [r[3] for r in rows if r[3]]
    print("resolved %d / %d" % (len(got), len(rows)))
    meta = esummary(sorted(set(got)))
    print("metadata for %d" % len(meta))
    with open("refs_raw.json", "w") as f:
        json.dump({"rows": rows, "meta": meta}, f)
    for sec, q, note, pmid in rows:
        if not pmid:
            print("  UNRESOLVED:", q[:80])
            continue
        m = meta.get(pmid, {})
        au = m.get("authors", [])
        a1 = au[0]["name"] if au else "?"
        print(f"  [{sec:>2}] {pmid:>9}  {a1:<18} {m.get('source','?'):<22} "
              f"{(m.get('pubdate','?') or '?')[:4]}  {m.get('title','?')[:88]}")


if __name__ == "__main__":
    main()
