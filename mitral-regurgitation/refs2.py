import json, time
from refs import esearch, esummary
R2 = [
 ("1","MITRA-FR COAPT lessons from two trials mitral regurgitation","Why the two trials differed: contemporaneous editorial analysis."),
 ("2","Utsunomiya proximal isovelocity surface area in vitro","PISA validated in vitro, where the orifice really is circular."),
 ("2","Hopmeyer flow convergence non-planar orifice geometry","Non-hemispheric convergence breaks the PISA assumption."),
 ("2","Yosefy proximal flow convergence functional mitral regurgitation elongated","The functional orifice is elongated, and PISA mis-measures it."),
 ("2","Buck corrective PISA eccentric mitral regurgitation","Correcting PISA for a non-circular orifice."),
 ("2","Zoghbi recommendations noninvasive evaluation valvular regurgitation 2017","The guideline grading scheme and its own caveats."),
 ("2","Uretsky discordance echocardiography magnetic resonance mitral regurgitation","Echo and CMR disagree; the model treats this as information."),
 ("2","Thomas physical and physiological determinants proximal isovelocity","Determinants and limits of the PISA method."),
 ("2","Dujardin grading of mitral regurgitation quantitative Doppler","Quantitative Doppler grading and its cut-points."),
 ("2","Militaru validation mitral regurgitation quantification cardiac magnetic resonance","Cross-modality validation."),
 ("3","giant v wave mitral regurgitation wedge pressure","The giant v wave: the atrial denominator made visible."),
 ("3","Barbier left atrial relaxation and chamber stiffness in humans","Left atrial chamber stiffness measured directly."),
 ("3","Stefanadis left atrial function pressure volume relation","The left atrial pressure-volume relation: denominator 1."),
 ("3","acute mitral regurgitation management cardiogenic shock papillary","Management of acute severe MR."),
 ("3","Yoran dynamic aspects of acute mitral regurgitation atrial compliance","Atrial compliance as an explicit determinant of the presentation."),
 ("3","left atrial compliance pulmonary hypertension mitral regurgitation","Atrial compliance and the pulmonary consequences."),
 ("3","Grose left ventricular volume and function chronic mitral regurgitation","Chronic ventricular adaptation to MR."),
 ("4","Enriquez-Sarano echocardiographic prediction survival surgical correction mitral regurgitation","The origin of the EF 60% operative threshold."),
 ("4","Starling left ventricular contractile function after valve surgery mitral regurgitation","Contractile function measured before and after operation."),
 ("4","Carabello assessment preoperative left ventricular function mitral regurgitation wall stress","Why EF misleads in MR: the wall-stress correction."),
 ("4","Kelly effective arterial elastance humans","Effective arterial elastance in humans."),
 ("4","Klotz single-beat estimation end-diastolic pressure volume relationship","The EDPVR form the model uses."),
 ("4","Burkhoff assessment of systolic and diastolic ventricular properties pressure-volume","Pressure-volume analysis as used here."),
 ("5","Schwammenthal dynamics of mitral regurgitant flow and orifice area","The orifice is dynamic within systole."),
 ("5","Ormiston size and motion of the mitral annulus in man","Normal annular size and motion."),
 ("5","Levine three-dimensional echocardiographic reconstruction mitral valve saddle prolapse","The annular saddle shape."),
 ("5","Otsuji isolated annular dilation functional mitral regurgitation atrial fibrillation","A key negative result: annular dilation alone is often not enough."),
 ("5","He dynamic strain mitral valve leaflet in vitro","Leaflet mechanics."),
 ("6","Grande-Allen apparently normal mitral valves heart failure derangement","The 'normal' valve in heart failure is not normal."),
 ("7","Carabello volume overload hypertrophy mitral regurgitation experimental","Volume-overload hypertrophy: thin-walled, not protective."),
 ("7","Corporan extracellular matrix remodeling mitral regurgitation volume overload","Matrix remodelling in chronic MR."),
 ("7","Zile myocardial stiffness collagen titin heart failure preserved ejection fraction","Collagen and titin as the two stiffness terms."),
 ("7","Kang early surgery asymptomatic severe mitral regurgitation","Early surgery versus watchful waiting."),
 ("7","Rosenhek outcome of watchful waiting asymptomatic severe mitral regurgitation","What watchful waiting actually delivers."),
 ("8","Ghoreishi pulmonary hypertension survival mitral valve surgery","PH after operation: the second barrier, measured."),
 ("8","Vachiery pulmonary hypertension due to left heart disease","Post-capillary PH: reversible and fixed components."),
 ("8","Bartko burden of secondary mitral regurgitation outcome","Outcome burden of secondary MR."),
 ("9","Kang sacubitril valsartan functional mitral regurgitation PRIME","PRIME: ARNI reduces functional MR."),
 ("9","Waagstein metoprolol idiopathic dilated cardiomyopathy","Beta-blockade recovers contractility."),
 ("9","Brater diuretic therapy braking phenomenon furosemide","The braking phenomenon."),
 ("9","nitroprusside mitral regurgitation regurgitant fraction afterload","Vasodilators reduce regurgitant fraction acutely."),
 ("9","ivabradine heart rate mitral regurgitation","Heart-rate reduction, which cuts both ways."),
 ("10","Hung mechanism recurrent ischemic mitral regurgitation annuloplasty","Recurrence as continued ventricular remodelling."),
 ("10","Michler two-year outcomes moderate ischemic mitral regurgitation coronary bypass","Adding a ring to revascularisation."),
 ("10","Breithardt cardiac resynchronization functional mitral regurgitation","CRT reduces functional MR by restoring closing force."),
 ("10","mitral valve gradient after transcatheter edge-to-edge repair outcome","The iatrogenic gradient as a device trade-off."),
 ("10","residual mitral regurgitation after transcatheter edge-to-edge repair outcomes","Residual MR after TEER."),
 ("10","chordal preservation left ventricular function mitral valve replacement","Chordal preservation and the ventricular scaffold."),
 ("11","Grigioni atrial fibrillation degenerative mitral regurgitation outcome","AF in MR: determinants and prognosis."),
 ("11","left atrial reservoir function mitral regurgitation","Atrial reservoir function."),
 ("11","atrial fibrillation outcome percutaneous mitral valve repair","AF and device outcomes."),
 ("12","mrgsolve R package simulation ordinary differential equation pharmacometrics","The simulation engine the R model targets."),
 ("12","Smith minimal cardiovascular system model ventricular interaction valve dynamics","A minimal closed-loop model with valve dynamics."),
 ("12","finite element model mitral regurgitation left ventricle simulation","A finite-element counterpart to the lumped approach."),
]
rows=[]
for sec,q,note in R2:
    pid=esearch(q); rows.append((sec,q,note,pid)); time.sleep(0.34)
got=[r[3] for r in rows if r[3]]
meta=esummary(sorted(set(got)))
json.dump({"rows":rows,"meta":meta},open("refs_raw2.json","w"))
print("resolved %d/%d"%(len(got),len(rows)))
for sec,q,note,pid in rows:
    if not pid: print("  STILL UNRESOLVED:",q); continue
    m=meta.get(pid,{}); au=m.get("authors",[])
    print(f"  [{sec:>2}] {pid:>9} {(au[0]['name'] if au else '?'):<16} {m.get('source','?'):<20} {(m.get('pubdate','?') or '?')[:4]} {m.get('title','?')[:86]}")
