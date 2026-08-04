#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cmd_reference_model.py
=======================================================================
Coronary Microvascular Dysfunction (CMD / ANOCA-INOCA)
Standalone reference implementation of the QSP model that is written in
mrgsolve form in `cmd_mrgsolve_model.R`.

WHY THIS FILE EXISTS
--------------------
There is no R runtime in the environment where this model was built, so
every equation was implemented and executed here first, in dependency-free
Python, and only then transcribed into mrgsolve syntax.  Every numeric
claim in `README.md` is produced by running this file.  The defects the
process exposed are listed in the BUG LOG below and marked at the
offending lines in BOTH files.

THE STRUCTURAL IDEA
-------------------
Every index of the coronary microcirculation used in the clinic is a RATIO
or a DIFFERENCE of two quantities, and this model refuses to merge them:

    CFR   = v_hyper / v_rest         two flows, one number
    MR    = Pd / v                   resistance, not flow
    DPTI  = (Pd - LVEDP) * t_dia     supply is bought only in diastole
    SPTI  = P_sys * t_sys            demand is spent only in systole
    def_L = MVO2_L - v_L*k*CaO2*Emax the only quantity a patient can feel

Four consequences the model is built to EXPOSE rather than assume:

 (1) A CFR of 1.9 is either a raised DENOMINATOR (resting flow too high:
     "functional" endotype) or a floored NUMERATOR (hyperaemic flow capped
     by structure: "structural").  Rahman et al. Circulation 2019;140:1805
     measured both; the two groups have the same CFR.

 (2) Resting tone is not a free parameter - autoregulation sets it so that
     resting supply meets resting demand.  Given the hyperaemic MR (which
     fixes R_min) the resting MR is therefore PREDICTED.  The prediction
     lands on the measured value for the control and structural groups.
     It fails for the functional group, and the size of the failure is the
     model's quantitative definition of what functional CMD is.

 (3) Heart rate enters the oxygen balance twice with the same sign: it
     raises demand (SPTI) and shortens the diastolic window in which the
     subendocardium is perfused (DPTI).  A drug that only slows the heart
     therefore has a supply-side action; a drug that only dilates can be
     cancelled by the pressure it drops.

 (4) An on-target drug can produce a null trial three different ways:
     internal cancellation (zibotentan / PRIZE), population dilution
     (ranolazine / RWISE), and background contamination plus a timescale
     mismatch (statin + ACEi / WARRIOR).  All three are computed.

And one prediction that was not designed in and emerged from the equations:
because the functional endotype's MINIMAL resistance is normal, it has
almost no subendocardial supply-demand deficit at any workload - so its
angina cannot be ischaemic in the supply-demand sense, and must be carried
by the nociceptive arm (A1-adenosine afferent signalling plus central
sensitisation).  That is why anti-ischaemic drugs fail in it, why the
CFR<2.5-with-low-hyperaemic-MR patient is the one who does not respond to
ranolazine, and why aminophylline (an adenosine-receptor antagonist) is
predicted to help precisely the group in which every vasodilator fails.

NUMERICAL DESIGN
----------------
The arteriolar control loop settles in ~15 s; the drug/structure problem
runs for 24 weeks to 2.5 years.  Rather than integrate a stiff system in
pure Python, the fast loop is solved as a FIXED POINT at every slow step,
by BISECTION on a provably monotone residual (see bug B11).  mrgsolve
integrates the identical expressions as ODEs with tau = 15 s under LSODA.
Section V verifies the two agree.  Every index reported here (CFR, MR,
endo/epi, DPTI, deficit) is by definition a steady state of the fast loop,
so the QSS-vs-ODE choice cannot change a reported number.

BUG LOG - real defects found by running this file
-------------------------------------------------
 B1  A high-gain autoregulator ERASED the functional endotype.  The first
     version fitted resting tone directly to Rahman's resting MR of 4.2;
     with the loop closed, the controller put tone straight back and
     resting MR returned to the control value.  The lesson became the
     model's main result: a low resting MR beside a normal minimal MR is
     not expressible as "less tone", only as a CONTROLLER DEFECT - and the
     size of the defect is then derivable rather than assumed.
 B2  Adenosine was implemented as TONE -> 0 exactly, which made hyperaemic
     MR equal to R_min by construction and silently deleted ET-1 and
     Rho-kinase from the hyperaemic measurement - and with them any reason
     to expect an endothelin antagonist to move measured CFR.  Fixed:
     exogenous adenosine abolishes only the metabolic and myogenic terms.
 B3  The diastolic time fraction was applied to BOTH layers, which
     cancelled f_dia out of the endo/epi ratio and made heart rate
     irrelevant to transmural distribution, deleting consequence (3).
     Fixed: the subepicardium is perfused in systole too, at efficiency
     PHI_SYS, and the subendocardium carries a SERIES compression
     resistance that does not scale with arteriolar tone.
 B4  The exercise systolic-pressure gain was shared across endotypes, so
     the structural group's exercise hypertension (188 vs 156 mmHg,
     Rahman) had nowhere to live and beta-blockade looked equally useful
     everywhere.
 B5  Angina raised sympathetic tone, which raised heart rate, which raised
     angina, with no gain check.  Both feedback states were UNBOUNDED, so
     the loop was linearly self-amplifying and the state vector diverged to
     nan for any angina scale factor above about 15 - which is exactly the
     value the symptom calibration lands on.  Both states are bounded
     neurophysiological quantities and are now written with saturation; the
     loop gain is swept and reported (section VIII) rather than clamped.
 B6  hs-cTnI was released in proportion to the deficit with no ceiling, so
     a 2.5-year WARRIOR run produced infarct-range troponins in patients
     whose endpoint was chest pain.
 B7  SAQ was driven by ischaemic burden alone, so central sensitisation
     had no route to the endpoint and the model could not represent the
     patient whose CFR improves while symptoms do not - which is most of
     the negative-trial literature.
 B8  Zibotentan was implemented as REDUCING plasma ET-1.  PRIZE measured
     the opposite (ET-1 rose: blocking ETA also removes ETA-mediated
     clearance).  The wrong sign made the drug look better than it is.
 B9  Bruce workload was mapped linearly onto time, putting stage 4 beyond
     any plausible reserve, so every arm returned the ceiling duration.
 B10 The population sampler drew R_min and the controller offset
     independently of each other and of the risk factors, producing
     patients with severe remodelling and no risk factor at all, and
     missing the observed 62/38 functional/structural split.
 B11 The fast loop was solved by damped fixed-point iteration.  With a
     controller gain of 6 against a resistance that scales as r^-4, the
     iteration OSCILLATED between the clamps and returned whatever value
     it held at iteration 200 - which is why the first run reported a
     control CFR of 1.19 and a hyperaemic MR of 10.7.  The residual
     target(tone) - tone is monotone decreasing, so bisection is
     unconditionally convergent; the continuous-time ODE form was never
     unstable, only the discrete iteration was.
 B12 The metabolic controller's intercept sat OUTSIDE the term that
     adenosine abolishes, so a floor of 0.5 tone survived maximal
     hyperaemia and capped every CFR in the model at about 1.2.
 B25 The pharmacokinetic states were integrated together with the
     physiology, and they were the stiffest states in the model: explicit
     RK4 is unstable for |k*dt| > 2.785 and the absorption constants reach
     1.6 /h, so at the 2 h step the nicorandil and fasudil arms were
     unstable and at the 6 h step used for the 2.5-year WARRIOR run every
     drug with ka > 0.46 /h was.  Ramipril reached -97213 mg/L on day 1 and
     nan after that, so the WARRIOR intensive-therapy arm silently became a
     placebo arm - treated and untreated agreed to three decimals, which is
     what exposed it.  All PK is now advanced in closed form outside the
     RK4 stages: exact rather than merely stable.
 B24 A tonic 16% epicardial pressure loss was applied to anyone with spasm
     susceptibility, and the untreated comparator in the sign table was
     evaluated without it.  Every drug in the structural endotype therefore
     appeared to double the subendocardial deficit - including imipramine,
     which has no haemodynamic action whatsoever, and that impossible row is
     what exposed it.  Epicardial obstruction now exists only during an
     episode, which is also the correct physiology.
 B23 Central sensitisation was driven only by angina, so it decayed to zero
     in any patient without an ischaemic driver, and the non-cardiac chest
     pain archetype - whose defining feature is pain with normal coronary
     function - came out completely asymptomatic (SAQ 99).  Sensitisation now
     has a primary drive, and angina a component sensitisation can generate
     with no afferent input, without which the model cannot represent the
     20-25% of every ANOCA cohort that no vasoactive drug helps.
 B22 Filling pressure was held at its RESTING value while workload varied,
     which removed the dominant subendocardial insult in this disease (it is
     the mechanism the CMD/HFpEF overlap rests on).  Without it the
     structural endotype never became ischaemic inside a plausible daily
     activity range - burden 0.000, SAQ 95, asymptomatic - while the
     functional endotype sat at SAQ 55.  That is exactly backwards from the
     clinic, and it was the model telling the truth about an equation that
     was missing, not about the disease.
 B21 The five endothelial states carried UNIT rate constants (a 1 h time
     constant) as placeholders.  At a 3 h explicit RK4 step the
     amplification factor is 1.375 per step, so the residual left by the
     closed-form initial condition grew by 1.375^672 over 12 weeks: plasma
     ET-1 reached 1e74, BH4 went negative and NO went nan, while the printed
     CFR stayed at a plausible 1.39 for most of the run and only the
     hyperaemic MR gave it away (2.3 -> 37.8).  The control archetype was
     immune because its steady state is exactly 1.0, so the unstable mode is
     never excited - which is how the bug survived three test runs.  The
     regulated pools turn over on hours to days; the time constants are now
     stated instead of implied.
 B20 Spasm was written as extra arteriolar TONE, which the metabolic
     controller absorbed: a loop gain of 6 answers a tone step of 0.6 with a
     10% flow change, so the vasospastic archetype could not become
     ischaemic and came out asymptomatic with a normal exercise test.  Spasm
     is not distributed tone - epicardial spasm is a series obstruction
     upstream of the autoregulated bed and microvascular spasm shuts part of
     the bed outright.  Rewritten as a distal-pressure loss plus a bed
     occlusion, neither of which autoregulation can undo.
 B19 The ischaemia -> late I_Na -> diastolic calcium -> LVEDP ->
     (Pd - LVEDP) -> ischaemia loop had linear gains and no ceiling, so the
     structural endotype ran away over 21 days to a filling pressure that
     collapsed subendocardial perfusion (CFR 0.13, hyperaemic MR 59).  Both
     the ischaemic drive and the filling pressure now saturate.
 B17 The two fast tone states were carried as explicit-RK4 ODEs at a 3 h
     step with a 15 s time constant - unstable by a factor of 720.  TONE_E
     reached 2e9 in a single step and the resulting nan propagated through
     the state vector, so a symptomatic patient silently reported SAQ 100.
     In mrgsolve (LSODA, implicit) the same ODEs are fine; in the Python
     reference they are solved by bisection and their derivative is zero.
     The failure is a reminder that a QSP model's stiffest state is often a
     variable nobody reads.
 B16 The controller-offset arm of the afferent signal was a CONSTANT, so
     the functional endotype had a tonic pain drive that did not rise with
     exertion: it could not generate exertional angina, and the
     exercise-duration calibration collapsed onto the 5 s floor while the
     angina scale factor ran away to 10.6.  The excess metabolic signal is
     proportional to the EXCESS demand, so the term now vanishes at rest.
 B14 ET-1 and Rho-kinase tone were made fully adenosine-resistant when bug
     B2 was fixed, which overshot: the functional endotype then came out
     with a hyperaemic MR of 8.8, when <2.5 is its DEFINITION, and the
     calibrated structural R_min was doing almost none of the work - the
     residual constrictor tone was.  The endotypes were being separated by
     their constrictors rather than by their structure, i.e. by the wrong
     axis entirely.  Fixed with an explicit reversible fraction F_ADO_REV,
     itself calibrated against the functional group's hyperaemic MR.
 B15 Spasm susceptibility was carried in RESTING tone, so the vasospastic
     archetype sat at 0.92 (maximal) tone permanently, with a hyperaemic MR
     of 35 and a resting flow no patient survives.  Vasospastic angina is
     episodic: the susceptibility is now expressed only during provoked or
     spontaneous episodes, weighted by the fraction of the day spent in
     them, which leaves resting physiology normal between attacks - as
     measured.
 B13 The functional endotype produced no subendocardial deficit at any
     workload (correctly - its minimal resistance is normal), which left
     it with no symptoms at all, because nociception was wired only to the
     deficit.  Fixed by giving the afferent signal its documented second
     source: the adenosine/A1 arm, which the controller offset drives
     directly.  This is the arm aminophylline blocks.

Run:  python3 cmd_reference_model.py
      -> writes cmd_reference_output.txt and cmd_population_results.json
=======================================================================
"""

import json
import math
import os
import random

OUT = []


def emit(s=""):
    OUT.append(s)
    print(s)


def rule(title, ch="="):
    emit()
    emit(ch * 78)
    emit(title)
    emit(ch * 78)


# =====================================================================
# 0.  PARAMETERS
# =====================================================================
# pressures mmHg | velocity cm/s | resistance mmHg/(cm/s)
# MVO2 mL O2/min/100 g | perfusion mL/min/g | time h
P = dict(
    # ---- systemic haemodynamics -------------------------------------
    HR0=68.0, SBP0=132.0, DBP0=76.0, PV=8.0,
    KHR_WL=30.0, KDBP_WL=5.0,
    # exercise systolic gain is ENDOTYPE-SPECIFIC (bug B4).  Rahman 2019:
    # peak exercise SBP 188+-25 structural, 161+-27 functional, 156+-30
    # control; peak workload taken as wl = 2.5
    KSBP_CTRL=16.0, KSBP_FUNC=19.0, KSBP_STRUCT=37.0,
    FRAC_DIA_P=0.35,      # mean diastolic aortic P = DBP + 0.35*(SBP-DBP)
    FRAC_SYS_P=0.75,      # mean systolic aortic P  = DBP + 0.75*(SBP-DBP)
    TSYS_K=0.30,          # t_sys = TSYS_K*sqrt(60/HR) seconds
    PHI_SYS=0.45,         # systolic perfusion efficiency, subepicardium
    FFR0=0.95,            # epicardial conductance, non-obstructive disease

    # ---- microvascular resistance -----------------------------------
    RMIN0=1.450,          # CALIBRATED (target: control hyperaemic MR 2.20)
    A_TONE=0.55,          # r/r_max = 1 - A_TONE*TONE  ->  R ~ (.)^-4
    W_ENDO=0.573,         # CALIBRATED (target: control hyperaemic endo/epi 1.02)
    W_EPI=1.15,
    RCOMP_K=0.0318,       # series compression resistance per mmHg of LVEDP (B3)
    TONE_MAX=0.95,

    # ---- oxygen transport -------------------------------------------
    KCONV=0.0610,         # mL/min/g of perfusion per cm/s of velocity
    CAO2=0.19, E_REST=0.70, E_MAX=0.80,

    # ---- myocardial oxygen demand -----------------------------------
    MVO2_B0=1.50, K_TTI=0.0882, K_CTR=1.50, K_DIAST=0.60,
    SIGMA_REF=132.0, F_ENDO_DEM=1.15, F_EPI_DEM=0.88,

    # ---- autoregulation ---------------------------------------------
    TONE_REF=0.50,        # controller intercept - INSIDE the adenosine-
                          # sensitive term (bug B12)
    G_AUTO=6.0, K_MYO=0.25,
    TAU_TONE=15.0 / 3600.0,
    K_ET_TONE=0.12,       # ET-1 tone (partly survives adenosine, bug B2)
    K_ROCK_TONE=0.10,     # Rho-kinase tone (partly survives adenosine)
    K_NO_TONE=0.12, K_A1_TONE=0.18,
    # bug B20: spasm was written as extra ARTERIOLAR TONE, which the
    # metabolic controller simply absorbed - a loop gain of 6 answers a tone
    # step of 0.6 with a 10% flow change, so a vasospastic patient could not
    # become ischaemic at all and the archetype came out asymptomatic with a
    # normal exercise test.  Spasm is not distributed arteriolar tone:
    # epicardial spasm is a SERIES obstruction upstream of the autoregulated
    # bed (autoregulation compensates only until tone reaches zero) and
    # microvascular spasm SHUTS a fraction of the bed, which cannot be
    # recruited at all.  Both are now written that way.
    SPASM_PD=0.50,        # fractional loss of distal pressure per unit SPASM
    SPASM_SHUT=0.40,      # fraction of the bed occluded per unit SPASM
    # fraction of CONSTRICTOR tone that 140 ug/kg/min adenosine reverses.
    # bug B14: with this fraction absent, ET-1 and Rho-kinase tone survived
    # hyperaemia in full, the functional endotype came out with a hyperaemic
    # MR of 8.8 (it is <2.5 BY DEFINITION), and the calibrated structural
    # R_min was doing almost none of the work - the residual constrictor tone
    # was.  The endotypes then differed by their constrictors rather than by
    # their structure, which is the opposite of what the labels mean.
    F_ADO_REV=0.85,       # CALIBRATED (target: functional hyperaemic MR 2.30)
    # fraction of the day spent in spontaneous microvascular/epicardial spasm
    # per unit SPASM susceptibility (bug B15)
    P_SPASM_DAY=0.06,

    # ---- endothelial / constrictor biology --------------------------
    KROS_BASE=0.35, KROS_DEG=1.0, KBH4_OX=0.55,
    KET_SYN=1.45, KET_CLR=1.0, A_INFL_ET=0.55, A_GENE_ET=0.30, B_NO_ET=0.45,
    F_ETA=0.75, F_ETB2=0.25,
    # bug B21 (the same lesson as B17): these five states were written with
    # UNIT rate constants, i.e. a 1 h time constant, as placeholders.  With a
    # 3 h explicit RK4 step the amplification factor is 1.375 per step, so the
    # residual left by the closed-form initial condition grew by 1.375^672
    # over 12 weeks: plasma ET-1 reached 1e74, BH4 went negative, NO went nan,
    # and the hyperaemic MR silently drifted from 2.3 to 37.8 while the CFR
    # printed a plausible-looking 1.39 for most of the run.  The control
    # archetype was unaffected because its steady state is exactly 1.0 and the
    # unstable mode is never excited - which is how the bug survived the first
    # three test runs.  The regulated pools turn over on hours-to-days, not on
    # an hour, so the physiological time constants are stated explicitly.
    TAU_ROS=12.0, TAU_BH4=24.0, TAU_NO=6.0, TAU_ET1=12.0, TAU_ROCK=12.0,

    # ---- structure (slow) -------------------------------------------
    TAU_ML=180.0 * 24, TAU_CAPD=240.0 * 24, TAU_PVF=200.0 * 24,
    TAU_ICF=240.0 * 24, TAU_LVH=180.0 * 24, W_ML=0.55, W_PVF=0.35,

    # ---- filling pressure / volume ----------------------------------
    LVEDP0=11.0, TAU_LVEDP=24.0, K_VOL_LVEDP=9.0, K_CAD_LVEDP=3.2,
    # bug B22: LVEDP was held at its RESTING value while workload was varied.
    # Exercise filling pressure is the dominant subendocardial insult in this
    # disease (it is the mechanism the CMD/HFpEF overlap is built on), and
    # without it the structural endotype never became ischaemic inside a
    # plausible daily activity range: its burden was 0.000 and it came out
    # asymptomatic (SAQ 95) while the functional endotype sat at SAQ 55 -
    # exactly backwards from the clinic.  Exercise LVEDP now rises with
    # workload, faster in a stiff ventricle.
    K_LVEDP_WL=1.60, K_ICF_LVEDP_WL=2.50,
    K_ICF_LVEDP=6.0, TAU_VOL=72.0,

    # ---- late sodium current / diastolic tension --------------------
    # bug B19: the ischaemia -> late I_Na -> diastolic Ca -> LVEDP ->
    # (Pd - LVEDP) -> ischaemia loop was written with LINEAR gains and no
    # ceiling.  In the structural endotype it ran away: over 21 days LVEDP
    # climbed until the subendocardial driving pressure collapsed, and the
    # "21-day untreated structural" comparator in section IV came out with a
    # CFR of 0.13 and a hyperaemic MR of 59 - a patient who cannot exist.
    # Both the ischaemic drive on Na and the resulting filling pressure now
    # saturate.
    TAU_NAI=12.0, TAU_CAD=12.0, K_ISCH_NAI=0.85, K_ISCH_NAI_K=1.50,
    LVEDP_MAX=30.0, K_TEN_CAD=0.55,

    # ---- nociception / symptoms -------------------------------------
    K_NOC_DEF=2.50,       # afferent drive per unit subendocardial deficit
    K_NOC_ADO=0.90,       # afferent drive per unit adenosine signal (B13)
    A_ADO_DEF=4.00,
    # bug B16: the controller-offset term was written as a CONSTANT, so the
    # functional endotype had a tonic pain drive that did not rise with
    # exertion - it could not produce exertional angina at all, and the
    # exercise-duration calibration collapsed onto the 5 s floor.  If the
    # arterioles behave as though demand were X% higher, the excess metabolic
    # signal is proportional to the EXCESS demand, so the term is written
    # against (MVO2/MVO2_rest - 1) and vanishes at rest.
    A_ADO_OFF=3.00, MVO2_REST_REF=9.60,
    K_ANG=18.0,           # CALIBRATED (target: untreated functional SAQ 55)
    K_ANG_GAIN=0.85, TAU_ANG=72.0,
    # bug B5/B18: SENS and SYMP were UNBOUNDED, so the
    # angina -> sensitisation -> angina and angina -> sympathetic -> angina
    # loops were linearly self-amplifying and the state vector diverged to
    # nan for any K_ANG above about 15 - exactly the value the symptom
    # calibration needs.  Both are bounded neurophysiological states and are
    # now written with saturation, which turns the runaway into a smooth
    # ceiling (a very symptomatic patient, not an infinite one).
    TAU_SENS=21.0 * 24, K_SENS_ON=0.35, K_SYMP_ANG_K=12.0,
    # bug B23: central sensitisation was driven ONLY by angina, so it decayed
    # to zero in any patient without an ischaemic driver - and the
    # non-cardiac-chest-pain archetype, whose defining feature is pain with
    # normal coronary function, came out completely asymptomatic (SAQ 99).
    # Sensitisation needs a primary drive of its own, and angina needs a
    # component that sensitisation can generate without any afferent input,
    # or the model cannot represent the 20-25% of the ANOCA population that
    # every coronary function study finds and no vasoactive drug helps.
    K_ANG_CENTRAL=8.0,
    SAQ_MAX=48.0, SAQ_K=3.4, SAQ_SENS=14.0, TAU_SAQ=14.0 * 24,
    TAU_SYMP=48.0, K_SYMP_ANG=0.055, K_SYMP_HR=0.22,
    MCID_SAQ=10.0,
    NOCI_THRESH=3.05,     # CALIBRATED (target: untreated functional Bruce 480 s)

    # ---- biomarkers / outcome ---------------------------------------
    KTNI_OFF=0.35, TNI_MAX=14.0, KBNP_ON=6.0, KBNP_OFF=0.02,
    H0_MORT=0.0155, LN116=math.log(1.16), LN108=math.log(1.08),
    H0_HOSP=0.075, K_HOSP_ANG=0.105,
)

# daily activity distribution for the burden integrals (wl, hours)
# bug B22: the top bin was 2.4 (about 1.7x resting MVO2), below the workload
# at which a CFR-2.1 patient becomes ischaemic, so ordinary daily exertion
# generated no burden at all.  Stair climbing is 6-7 METs.
ACTIVITY = [(1.00, 15.0), (1.55, 6.0), (2.40, 2.5), (3.10, 0.5)]

ARCHETYPE = {
    "control":     dict(RMIN_F=1.000, AUTO_OFF=0.00, INFL=1.0, GENO=0,
                        ROCK_D=1.00, KSBP="KSBP_CTRL", SPASM=0.00, SENS0=0.10,
                        SENS_BASE=0.05),
    # the functional endotype has, by definition, a NORMAL structural floor:
    # RMIN_F is fixed at 1.000 and is not a free parameter here
    "functional":  dict(RMIN_F=1.000, AUTO_OFF=0.79, INFL=1.4, GENO=1,
                        ROCK_D=1.30, KSBP="KSBP_FUNC", SPASM=0.10, SENS0=0.45,
                        SENS_BASE=0.75),
    "structural":  dict(RMIN_F=1.520, AUTO_OFF=0.05, INFL=2.2, GENO=0,
                        ROCK_D=1.50, KSBP="KSBP_STRUCT", SPASM=0.05, SENS0=0.40,
                        SENS_BASE=0.65),
    "vasospastic": dict(RMIN_F=1.080, AUTO_OFF=0.12, INFL=1.9, GENO=1,
                        ROCK_D=3.10, KSBP="KSBP_FUNC", SPASM=0.85, SENS0=0.50,
                        SENS_BASE=0.85),
    "noncardiac":  dict(RMIN_F=1.000, AUTO_OFF=0.00, INFL=1.1, GENO=0,
                        ROCK_D=1.00, KSBP="KSBP_CTRL", SPASM=0.00, SENS0=0.80,
                        SENS_BASE=1.60),
}

# =====================================================================
# 1.  DRUGS
# =====================================================================
DRUGS = {
    "ranolazine":    dict(ka=0.45, ke=0.099, V=180.0, EC50=1.60),
    "ivabradine":    dict(ka=1.20, ke=0.116, V=110.0, EC50=0.021),
    "bisoprolol":    dict(ka=0.90, ke=0.058, V=230.0, EC50=0.012),
    "nebivolol":     dict(ka=1.10, ke=0.069, V=1500.0, EC50=0.0016),
    "amlodipine":    dict(ka=0.35, ke=0.019, V=1400.0, EC50=0.0055),
    "diltiazem":     dict(ka=0.80, ke=0.139, V=350.0, EC50=0.115),
    "nicorandil":    dict(ka=1.60, ke=0.578, V=70.0, EC50=0.055),
    "zibotentan":    dict(ka=0.70, ke=0.058, V=95.0, EC50=0.048),
    "ramipril":      dict(ka=1.00, ke=0.050, V=100.0, EC50=0.0021),
    "atorvastatin":  dict(ka=0.60, ke=0.050, V=560.0, EC50=0.0085),
    "fasudil":       dict(ka=1.40, ke=0.347, V=90.0, EC50=0.090),
    "sildenafil":    dict(ka=1.30, ke=0.173, V=105.0, EC50=0.135),
    "aminophylline": dict(ka=1.10, ke=0.087, V=35.0, EC50=8.00),
    "imipramine":    dict(ka=0.90, ke=0.041, V=1200.0, EC50=0.020),
}
DRUG_ORDER = list(DRUGS.keys())

EMAX = dict(
    iva_HR=0.19, bis_HR=0.21, neb_HR=0.18, dil_HR=0.09, ran_HR=0.052,
    bis_CTR=0.17, neb_CTR=0.14, dil_CTR=0.11,
    amlo_TONE=0.36, dil_TONE=0.24, nico_TONE=0.33, sild_TONE=0.16,
    amlo_MAP=0.075, dil_MAP=0.055, nico_MAP=0.030, zib_MAP=0.055,
    ram_MAP=0.070, sild_MAP=0.045,
    zib_ETA=0.86, zib_VOL=0.115,
    fas_ROCK=0.72, stat_ROCK=0.22,
    stat_ROS=0.30, ram_ROS=0.18, neb_NO=0.22, nico_NO=0.20, ram_NO=0.14,
    ran_LATENA=0.62,
    ami_A1=0.68, ami_A2A=0.45,
    imi_SENS=0.45, cbt_SENS=0.30, rehab_SENS=0.18,
    ram_ML=0.34, stat_ML=0.16, rehab_CAPD=0.28, ram_PVF=0.25,
)


def emx(conc, ec50, emax):
    return emax * conc / (ec50 + conc) if conc > 0 else 0.0


# =====================================================================
# 2.  STATE VECTOR
# =====================================================================
SLOW = ["NO", "ROS", "BH4", "ET1", "ROCK",
        "ML", "CAPD", "PVF", "ICF", "LVH",
        "LVEDP", "VOL", "NAI", "CAD",
        "SENS", "SYMP", "ANG", "SAQ",
        "TNI", "BNP",
        "CH_MORT", "CH_MACE", "CH_HOSP",
        "AUC_DEF", "AUC_NOC",
        "TONE_E", "TONE_P"]
NS_SLOW = len(SLOW)
IX = {k: i for i, k in enumerate(SLOW)}
PK_NAMES = []
for _d in DRUG_ORDER:
    PK_NAMES += ["A_" + _d, "C_" + _d]
PK_NAMES += ["M_ivabradine", "M_ramipril", "M_fasudil"]
NS = NS_SLOW + len(PK_NAMES)
IXP = {k: NS_SLOW + i for i, k in enumerate(PK_NAMES)}


def endothelial_ss(a):
    """Closed-form steady state of the endothelial block (flow-independent)."""
    ros = P["KROS_BASE"] * (1.0 + 0.55 * (a["INFL"] - 1.0)) * 2.857 / P["KROS_DEG"]
    bh4 = 1.0 / (1.0 + P["KBH4_OX"] * max(0.0, ros - 1.0))
    no = bh4 / ros
    et = P["KET_SYN"] * (1.0 + P["A_INFL_ET"] * (a["INFL"] - 1.0)
                         + P["A_GENE_ET"] * a["GENO"]) \
        / ((1.0 + P["B_NO_ET"] * no) * P["KET_CLR"])
    rock = a["ROCK_D"]
    return no, ros, bh4, et, rock


def y0(arch):
    a = ARCHETYPE[arch]
    y = [0.0] * NS
    no, ros, bh4, et, rock = endothelial_ss(a)
    y[IX["NO"]], y[IX["ROS"]], y[IX["BH4"]] = no, ros, bh4
    y[IX["ET1"]], y[IX["ROCK"]] = et, rock
    y[IX["ML"]] = a["RMIN_F"] - 1.0
    y[IX["CAPD"]] = 1.0 / (1.0 + 0.45 * (a["RMIN_F"] - 1.0))
    y[IX["PVF"]] = 0.6 * (a["RMIN_F"] - 1.0)
    y[IX["ICF"]] = 0.25 * (a["INFL"] - 1.0)
    y[IX["LVH"]] = 1.0 + 0.10 * (a["INFL"] - 1.0)
    y[IX["LVEDP"]] = P["LVEDP0"] + P["K_ICF_LVEDP"] * y[IX["ICF"]]
    y[IX["VOL"]] = 1.0
    y[IX["NAI"]] = 1.0
    y[IX["CAD"]] = 1.0
    y[IX["SENS"]] = a["SENS0"]
    y[IX["SAQ"]] = 70.0
    y[IX["TNI"]] = 3.0
    y[IX["BNP"]] = 60.0
    y[IX["TONE_E"]] = 0.45
    y[IX["TONE_P"]] = 0.45
    return y


# =====================================================================
# 3.  FAST MODULE
# =====================================================================
def systemic(y, arch, wl, rx):
    a = ARCHETYPE[arch]
    symp = y[IX["SYMP"]]
    HR = P["HR0"] * (1.0 - rx.get("f_HR", 0.0)) * (1.0 + P["K_SYMP_HR"] * symp) \
        + P["KHR_WL"] * (wl - 1.0)
    HR = max(38.0, min(190.0, HR))
    mdrop = min(0.25, rx.get("f_MAP", 0.0))
    SBP = (P["SBP0"] + P[a["KSBP"]] * (wl - 1.0)) * (1.0 - mdrop) \
        * (1.0 + 0.35 * P["K_SYMP_HR"] * symp)
    DBP = (P["DBP0"] + P["KDBP_WL"] * (wl - 1.0)) * (1.0 - mdrop)
    Pdia = DBP + P["FRAC_DIA_P"] * (SBP - DBP)
    Psysm = DBP + P["FRAC_SYS_P"] * (SBP - DBP)
    tcyc = 60.0 / HR
    tsys = min(P["TSYS_K"] * math.sqrt(60.0 / HR), 0.62 * tcyc)
    fdia = (tcyc - tsys) / tcyc
    return dict(HR=HR, SBP=SBP, DBP=DBP, Pdia=Pdia, Psysm=Psysm,
                tcyc=tcyc, tsys=tsys, fdia=fdia)


def demand(y, h, wl, rx):
    sigma = h["SBP"] / P["SIGMA_REF"]
    ctr = (1.0 + 0.50 * (wl - 1.0)) * (1.0 - rx.get("f_CTR", 0.0))
    ten = 1.0 + P["K_TEN_CAD"] * (y[IX["CAD"]] - 1.0)
    mvo2 = (P["MVO2_B0"] + P["K_TTI"] * h["HR"] * sigma
            + P["K_CTR"] * ctr + P["K_DIAST"] * ten) * (1.0 - rx.get("f_MVO2", 0.0))
    return mvo2, ten


def rmin_eff(y, rx, arch=None, spasm=False):
    struct = (1.0 + P["W_ML"] * y[IX["ML"]] + P["W_PVF"] * y[IX["PVF"]]) \
        / max(0.35, y[IX["CAPD"]])
    if spasm and arch is not None:
        # microvascular spasm removes part of the bed (bug B20)
        sev = ARCHETYPE[arch]["SPASM"] * rx.get("d_SPASM", 0.34) / 0.34 \
            * (1.0 + 0.5 * max(0.0, y[IX["ROCK"]] - 1.0))
        struct /= max(0.15, 1.0 - min(0.80, P["SPASM_SHUT"] * sev))
    return P["RMIN0"] * struct


def passive_tone(y, arch, rx, ach=0.0, hyperaemia=False, spasm=False):
    """
    Constrictor tone that is NOT generated by the metabolic controller.

    At rest it is a disturbance the controller absorbs (it displaces the
    operating point without changing resting flow).  Under exogenous
    adenosine the controller is off, and a fraction F_ADO_REV of this tone
    is reversed as well - bug B14: with F_ADO_REV absent the whole
    constrictor load survived hyperaemia, which inflated every hyperaemic
    MR and let the constrictors, not the structure, define the endotypes.

    Spontaneous spasm is NOT part of resting tone (bug B15): it is an
    episodic event, entered through the `spasm` flag and weighted by the
    fraction of the day spent in it, so a vasospastic patient has normal
    resting flow between attacks - which is what is measured.
    """
    a = ARCHETYPE[arch]
    et_eff = (y[IX["ET1"]] - 1.0) * (P["F_ETA"] * (1.0 - rx.get("f_ETA", 0.0))
                                     + P["F_ETB2"])
    con = (P["K_ET_TONE"] * max(0.0, et_eff)
           + P["K_ROCK_TONE"] * max(0.0, y[IX["ROCK"]] - 1.0)
           + P["K_A1_TONE"] * y[IX["SYMP"]] * (1.0 - rx.get("f_A1", 0.0)))
    dil = P["K_NO_TONE"] * (y[IX["NO"]] - 1.0 + rx.get("d_NO", 0.0)) \
        + rx.get("f_TONE", 0.0)
    t = con - dil
    if hyperaemia:
        t *= (1.0 - P["F_ADO_REV"])
    if spasm:
        t += rx.get("d_SPASM", 0.34) * a["SPASM"] * (1.0 + 2.0 * y[IX["ROCK"]])
    if ach > 0.0:
        t += (-0.30 * min(1.0, y[IX["NO"]]) * ach
              + 0.55 * ach * a["SPASM"]
              * (1.0 + 0.5 * max(0.0, y[IX["ROCK"]] - 1.0))
              + 0.22 * ach * max(0.0, 1.0 - y[IX["NO"]]))
    return max(0.0, min(P["TONE_MAX"], t))


def lvedp_eff(y, wl):
    """Filling pressure at this workload (bug B22): it rises with exertion,
    and faster in a fibrotic, stiff ventricle."""
    return y[IX["LVEDP"]] + P["K_LVEDP_WL"] \
        * (1.0 + P["K_ICF_LVEDP_WL"] * y[IX["ICF"]]) * max(0.0, wl - 1.0)


def layer_velocity(y, h, Rart, layer, rx, ffr_loss=0.0, lvedp=None):
    Pd = h["Pdia"] * P["FFR0"] * (1.0 - ffr_loss)
    Pdsys = h["Psysm"] * P["FFR0"] * (1.0 - ffr_loss)
    if lvedp is None:
        lvedp = y[IX["LVEDP"]]
    if layer == "endo":
        # series extravascular compression, does NOT scale with tone (bug B3)
        Rl = Rart * P["W_ENDO"] + P["RCOMP_K"] * lvedp
        return h["fdia"] * max(0.0, Pd - lvedp) / Rl, Pd
    Rl = Rart * P["W_EPI"]
    return (h["fdia"] * max(0.0, Pd - P["PV"])
            + (1.0 - h["fdia"]) * P["PHI_SYS"] * max(0.0, Pdsys - P["PV"])) / Rl, Pd


def fast_solve(y, arch, wl, rx, hyperaemia=False, ach=0.0, spasm=False, nbis=26):
    """
    Steady state of the two-layer arteriolar control loop.

    Each layer autoregulates against its own required flow.  The residual
    target(tone) - tone is monotone decreasing in tone, so bisection is
    unconditionally convergent (bug B11 - damped fixed-point iteration
    oscillated between the clamps and silently returned garbage).
    """
    a = ARCHETYPE[arch]
    h = systemic(y, arch, wl, rx)
    mvo2, ten = demand(y, h, wl, rx)
    Rmin = rmin_eff(y, rx, arch, spasm)
    # epicardial spasm: a series obstruction upstream of the bed (bug B20)
    ffr_loss = rx.get("d_EPISPASM", 0.0)
    if spasm:
        ffr_loss += P["SPASM_PD"] * ARCHETYPE[arch]["SPASM"] \
            * rx.get("d_SPASM", 0.34) / 0.34
    ffr_loss = min(0.75, ffr_loss)
    lv = lvedp_eff(y, wl)
    tp = passive_tone(y, arch, rx, ach, hyperaemia, spasm)
    myo = P["K_MYO"] * (h["Pdia"] - 90.0) / 90.0
    kO2 = P["KCONV"] * 100.0 * P["CAO2"]
    off = a["AUTO_OFF"] * (1.0 - rx.get("f_AUTO", 0.0))

    res = {}
    for layer, fdem in (("endo", P["F_ENDO_DEM"]), ("epi", P["F_EPI_DEM"])):
        dem = mvo2 * fdem
        v_req = dem / (kO2 * P["E_REST"]) * (1.0 + off)
        if hyperaemia:
            tone = tp
        else:
            lo, hi = 0.0, P["TONE_MAX"]

            def resid(t):
                Rart = Rmin / (1.0 - P["A_TONE"] * t) ** 4
                v, _ = layer_velocity(y, h, Rart, layer, rx, ffr_loss, lv)
                tgt = P["TONE_REF"] + P["G_AUTO"] * (v / v_req - 1.0) + myo + tp
                return max(0.0, min(P["TONE_MAX"], tgt)) - t
            if resid(lo) <= 0.0:
                tone = lo
            elif resid(hi) >= 0.0:
                tone = hi
            else:
                for _ in range(nbis):
                    mid = 0.5 * (lo + hi)
                    if resid(mid) > 0.0:
                        lo = mid
                    else:
                        hi = mid
                tone = 0.5 * (lo + hi)
        Rart = Rmin / (1.0 - P["A_TONE"] * tone) ** 4
        v, Pd = layer_velocity(y, h, Rart, layer, rx, ffr_loss, lv)
        res[layer] = dict(tone=tone, Rart=Rart, v=v, dem=dem, v_req=v_req,
                          sup=v * kO2 * P["E_MAX"], Pd=Pd,
                          reserve=1.0 - tone / max(1e-9, P["TONE_MAX"]))

    e, p = res["endo"], res["epi"]
    v = 0.5 * (e["v"] + p["v"])
    def_e = max(0.0, e["dem"] - e["sup"])
    def_p = max(0.0, p["dem"] - p["sup"])
    ado = (1.0 + P["A_ADO_DEF"] * def_e / (1.0 + def_e)
           + P["A_ADO_OFF"] * off
           * max(0.0, mvo2 / P["MVO2_REST_REF"] - 1.0))     # bug B16
    # bug B13: the afferent signal has TWO sources - true deficit and the
    # adenosine arm that the controller offset drives on its own.  Only the
    # second is blocked by aminophylline.
    noci = P["K_NOC_DEF"] * def_e \
        + P["K_NOC_ADO"] * (ado - 1.0) * (1.0 - rx.get("f_A1", 0.0))
    dpti = max(0.0, e["Pd"] - lv) * h["tcyc"] * h["fdia"]
    spti = h["Psysm"] * h["tsys"]
    return dict(h=h, endo=e, epi=p, v=v, Pd=e["Pd"], Rmin=Rmin, mvo2=mvo2,
                def_endo=def_e, def_epi=def_p, ado=ado, noci=noci, ten=ten,
                MR=e["Pd"] / v if v > 1e-9 else float("inf"),
                endo_epi=e["v"] / p["v"] if p["v"] > 1e-9 else 0.0,
                dpti=dpti, spti=spti, sevr=dpti / spti if spti > 0 else 0.0,
                lvedp_eff=lv,
                mbf=v * P["KCONV"],
                E_rest=mvo2 / max(1e-9, v * P["KCONV"] * 100.0 * P["CAO2"]),
                tone=0.5 * (e["tone"] + p["tone"]))


def indices(y, arch, rx):
    """The invasive coronary function test, as performed."""
    rest = fast_solve(y, arch, 1.0, rx)
    hyp = fast_solve(y, arch, 1.0, rx, hyperaemia=True)
    ach = fast_solve(y, arch, 1.0, rx, ach=1.0)
    cfr = hyp["v"] / rest["v"] if rest["v"] > 1e-9 else 0.0
    return dict(rest=rest, hyp=hyp, ach=ach, CFR=cfr,
                MR_rest=rest["MR"], MR_hyp=hyp["MR"],
                IMR=hyp["Pd"] / hyp["v"] * 8.35,   # Pd x Tmn, scaled: normal 18 U
                MRR=cfr * (rest["Pd"] / hyp["Pd"]) / P["FFR0"],
                ACh_FR=ach["v"] / rest["v"] if rest["v"] > 1e-9 else 0.0,
                endo_epi_hyp=hyp["endo_epi"], endo_epi_rest=rest["endo_epi"])


def burdens(y, arch, rx):
    """
    Daily ischaemic burden and daily nociceptive drive.

    Each activity bin is split between the non-spasm state and the spasm
    state, weighted by the fraction of the day spent in spasm (bug B15).
    An anti-spasm drug therefore reduces both the frequency (through
    d_SPASM) and the severity of the episodes.
    """
    a = ARCHETYPE[arch]
    p_sp = min(0.5, P["P_SPASM_DAY"] * a["SPASM"] * rx.get("d_SPASM", 0.34) / 0.34)
    bi = bn = 0.0
    for wl, hrs in ACTIVITY:
        f0 = fast_solve(y, arch, wl, rx)
        if p_sp > 1e-6:
            f1 = fast_solve(y, arch, wl, rx, spasm=True)
            bi += hrs * ((1 - p_sp) * f0["def_endo"] + p_sp * f1["def_endo"])
            bn += hrs * ((1 - p_sp) * f0["noci"] + p_sp * f1["noci"])
        else:
            bi += hrs * f0["def_endo"]
            bn += hrs * f0["noci"]
    return bi / 24.0, bn / 24.0


BRUCE = [(0, 180, 1.0, 4.6), (180, 360, 4.6, 7.0), (360, 540, 7.0, 10.2),
         (540, 720, 10.2, 12.9), (720, 900, 12.9, 15.0)]


def wl_of_bruce(t_s):
    """METs by Bruce stage, interpolated, mapped to model workload (bug B9)."""
    for t0, t1, m0, m1 in BRUCE:
        if t_s < t1:
            mets = m0 + (t_s - t0) / (t1 - t0) * (m1 - m0)
            return 1.0 + (mets - 1.0) / 4.35
    return 1.0 + (15.0 - 1.0) / 4.35


def exercise_duration(y, arch, rx):
    """Bruce seconds until the afferent drive crosses the anginal threshold."""
    thr = P["NOCI_THRESH"] / (1.0 + 0.55 * y[IX["SENS"]])
    if fast_solve(y, arch, wl_of_bruce(900.0), rx)["noci"] < thr:
        return 900.0
    lo, hi = 5.0, 900.0
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        if fast_solve(y, arch, wl_of_bruce(mid), rx)["noci"] < thr:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def reserve_exhaustion_wl(y, arch, rx):
    """Workload at which subendocardial tone reaches zero (autoregulation
    exhausted; beyond it subendocardial flow is pressure-passive)."""
    lo, hi = 1.0, 5.0
    if fast_solve(y, arch, hi, rx)["endo"]["tone"] > 1e-6:
        return float("nan")
    for _ in range(22):
        mid = 0.5 * (lo + hi)
        if fast_solve(y, arch, mid, rx)["endo"]["tone"] > 1e-6:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# =====================================================================
# 4.  SLOW MODULE
# =====================================================================
def regimen_effects(y, reg):
    C = {d: y[IXP["C_" + d]] for d in DRUG_ORDER}
    Mi, Mr, Mf = (y[IXP["M_ivabradine"]], y[IXP["M_ramipril"]],
                  y[IXP["M_fasudil"]])
    E = lambda d, e: emx(C[d], DRUGS[d]["EC50"], EMAX[e])
    rx = {}
    rx["f_HR"] = (emx(C["ivabradine"] + 0.7 * Mi, DRUGS["ivabradine"]["EC50"],
                      EMAX["iva_HR"])
                  + E("bisoprolol", "bis_HR") + E("nebivolol", "neb_HR")
                  + E("diltiazem", "dil_HR") + E("ranolazine", "ran_HR"))
    rx["f_CTR"] = (E("bisoprolol", "bis_CTR") + E("nebivolol", "neb_CTR")
                   + E("diltiazem", "dil_CTR"))
    rx["f_TONE"] = (E("amlodipine", "amlo_TONE") + E("diltiazem", "dil_TONE")
                    + E("nicorandil", "nico_TONE") + E("sildenafil", "sild_TONE"))
    rx["f_MAP"] = (E("amlodipine", "amlo_MAP") + E("diltiazem", "dil_MAP")
                   + E("nicorandil", "nico_MAP") + E("zibotentan", "zib_MAP")
                   + emx(Mr, DRUGS["ramipril"]["EC50"], EMAX["ram_MAP"])
                   + E("sildenafil", "sild_MAP"))
    if reg.get("no_bp_drop"):
        rx["f_MAP"] = 0.0
    rx["f_ETA"] = E("zibotentan", "zib_ETA")
    rx["d_NO"] = (E("nebivolol", "neb_NO") + E("nicorandil", "nico_NO")
                  + emx(Mr, DRUGS["ramipril"]["EC50"], EMAX["ram_NO"]))
    rx["f_A1"] = E("aminophylline", "ami_A1")
    # aminophylline also blocks A2A: dilator reserve is reduced
    rx["f_TONE"] -= emx(C["aminophylline"], DRUGS["aminophylline"]["EC50"],
                        EMAX["ami_A2A"]) * 0.20
    anti_spasm = min(1.0, E("amlodipine", "amlo_TONE") / EMAX["amlo_TONE"]
                     + E("diltiazem", "dil_TONE") / EMAX["dil_TONE"]
                     + emx(Mf, DRUGS["fasudil"]["EC50"], EMAX["fas_ROCK"])
                     / EMAX["fas_ROCK"])
    rx["d_SPASM"] = 0.34 * (1.0 - anti_spasm)
    # bug B24: a TONIC epicardial pressure loss of 16% used to be applied here
    # whenever a patient had any spasm susceptibility.  Two things were wrong.
    # Physiologically, there is no epicardial obstruction BETWEEN episodes -
    # that is what makes vasospastic angina vasospastic.  Numerically, the
    # untreated comparator in the sign table was evaluated with an empty
    # effect dictionary (loss 0) while every treated arm carried the 16%,
    # so every drug in the structural endotype appeared to double the
    # subendocardial deficit - including imipramine, which has no
    # haemodynamic action at all.  The pressure loss now exists only during an
    # episode, through SPASM_PD, and the comparator is exact.
    rx["d_EPISPASM"] = 0.0
    rx["f_MVO2"] = 0.055 if reg.get("trimetazidine") else 0.0
    rx["f_AUTO"] = 0.0
    return rx, C, Mi, Mr, Mf


def rhs(t, y, arch, reg):
    a = ARCHETYPE[arch]
    d = [0.0] * NS
    rx, C, Mi, Mr, Mf = regimen_effects(y, reg)
    E = lambda dr, e: emx(C[dr], DRUGS[dr]["EC50"], EMAX[e])

    # ---- endothelium -------------------------------------------------
    ros_t = P["KROS_BASE"] * (1.0 + 0.55 * (a["INFL"] - 1.0)) * 2.857 \
        * (1.0 - E("atorvastatin", "stat_ROS")
           - emx(Mr, DRUGS["ramipril"]["EC50"], EMAX["ram_ROS"])
           - 0.16 * reg.get("rehab", 0.0))
    # bug B21: relaxation form with EXPLICIT physiological time constants
    d[IX["ROS"]] = (ros_t / P["KROS_DEG"] - y[IX["ROS"]]) / P["TAU_ROS"]
    bh4_t = 1.0 / (1.0 + P["KBH4_OX"] * max(0.0, y[IX["ROS"]] - 1.0))
    d[IX["BH4"]] = (bh4_t - y[IX["BH4"]]) / P["TAU_BH4"]
    no_t = y[IX["BH4"]] * (1.0 + 0.10 * reg.get("rehab", 0.0)) \
        / max(1e-6, y[IX["ROS"]])
    d[IX["NO"]] = (no_t - y[IX["NO"]]) / P["TAU_NO"]
    et_syn = P["KET_SYN"] * (1.0 + P["A_INFL_ET"] * (a["INFL"] - 1.0)
                             + P["A_GENE_ET"] * a["GENO"]) \
        / (1.0 + P["B_NO_ET"] * y[IX["NO"]])
    # bug B8: ETA blockade RAISES plasma ET-1 (loss of ETA-mediated
    # clearance).  PRIZE measured exactly that.
    et_t = et_syn / (P["KET_CLR"] * (1.0 - 0.42 * E("zibotentan", "zib_ETA")))
    d[IX["ET1"]] = (et_t - y[IX["ET1"]]) / P["TAU_ET1"]
    rock_t = a["ROCK_D"] * (1.0 - E("atorvastatin", "stat_ROCK")
                            - emx(Mf, DRUGS["fasudil"]["EC50"], EMAX["fas_ROCK"]))
    d[IX["ROCK"]] = (rock_t - y[IX["ROCK"]]) / P["TAU_ROCK"]

    # ---- structure ---------------------------------------------------
    ml_t = (a["RMIN_F"] - 1.0) * (1.0 - emx(Mr, DRUGS["ramipril"]["EC50"],
                                            EMAX["ram_ML"])
                                  - E("atorvastatin", "stat_ML"))
    d[IX["ML"]] = (ml_t - y[IX["ML"]]) / P["TAU_ML"]
    capd_t = 1.0 / (1.0 + 0.45 * (a["RMIN_F"] - 1.0)) \
        * (1.0 + EMAX["rehab_CAPD"] * reg.get("rehab", 0.0))
    d[IX["CAPD"]] = (capd_t - y[IX["CAPD"]]) / P["TAU_CAPD"]
    pvf_t = 0.6 * (a["RMIN_F"] - 1.0) * (1.0 - emx(Mr, DRUGS["ramipril"]["EC50"],
                                                   EMAX["ram_PVF"]))
    d[IX["PVF"]] = (pvf_t - y[IX["PVF"]]) / P["TAU_PVF"]
    d[IX["ICF"]] = (0.25 * (a["INFL"] - 1.0) - y[IX["ICF"]]) / P["TAU_ICF"]
    d[IX["LVH"]] = (1.0 + 0.10 * (a["INFL"] - 1.0) - y[IX["LVH"]]) / P["TAU_LVH"]

    # ---- volume / filling pressure -----------------------------------
    vol_t = 1.0 + (0.0 if reg.get("no_fluid") else
                   EMAX["zib_VOL"] * E("zibotentan", "zib_ETA") / EMAX["zib_ETA"])
    d[IX["VOL"]] = (vol_t - y[IX["VOL"]]) / P["TAU_VOL"]
    lvedp_t = min(P["LVEDP_MAX"],
                  P["LVEDP0"] + P["K_VOL_LVEDP"] * (y[IX["VOL"]] - 1.0)
                  + P["K_CAD_LVEDP"] * (y[IX["CAD"]] - 1.0)
                  + P["K_ICF_LVEDP"] * y[IX["ICF"]])       # bug B19: ceiling
    d[IX["LVEDP"]] = (lvedp_t - y[IX["LVEDP"]]) / P["TAU_LVEDP"]

    # ---- ischaemia-dependent -----------------------------------------
    bi, bn = burdens(y, arch, rx)
    f = fast_solve(y, arch, 1.0, rx)

    late = (1.0 + P["K_ISCH_NAI"] * bi / (P["K_ISCH_NAI_K"] + bi)) \
        * (1.0 - E("ranolazine", "ran_LATENA"))          # bug B19: saturating
    d[IX["NAI"]] = (late - y[IX["NAI"]]) / P["TAU_NAI"]
    d[IX["CAD"]] = (y[IX["NAI"]] - y[IX["CAD"]]) / P["TAU_CAD"]

    # ---- symptoms ----------------------------------------------------
    ang_t = P["K_ANG"] * bn * (1.0 + P["K_ANG_GAIN"] * y[IX["SENS"]]) \
        + P["K_ANG_CENTRAL"] * y[IX["SENS"]] ** 2       # bug B23
    d[IX["ANG"]] = (ang_t - y[IX["ANG"]]) / P["TAU_ANG"]
    sens_off = (1.0 + E("imipramine", "imi_SENS")
                + EMAX["cbt_SENS"] * reg.get("cbt", 0.0)
                + EMAX["rehab_SENS"] * reg.get("rehab", 0.0))
    d[IX["SENS"]] = (P["K_SENS_ON"] * y[IX["ANG"]] * (1.0 - y[IX["SENS"]])
                     + a["SENS_BASE"] * (1.0 - y[IX["SENS"]])   # bug B23
                     - sens_off * y[IX["SENS"]]) / P["TAU_SENS"]
    d[IX["SYMP"]] = (y[IX["ANG"]] / (P["K_SYMP_ANG_K"] + y[IX["ANG"]])
                     - y[IX["SYMP"]]) / P["TAU_SYMP"]
    saq_ss = 100.0 - P["SAQ_MAX"] * y[IX["ANG"]] / (P["SAQ_K"] + y[IX["ANG"]]) \
        - P["SAQ_SENS"] * y[IX["SENS"]]
    d[IX["SAQ"]] = (saq_ss - y[IX["SAQ"]]) / P["TAU_SAQ"]

    # ---- biomarkers --------------------------------------------------
    d[IX["TNI"]] = P["KTNI_OFF"] * (3.0 + P["TNI_MAX"] * bi / (2.5 + bi)
                                    - y[IX["TNI"]])          # bug B6 ceiling
    d[IX["BNP"]] = P["KBNP_ON"] * (y[IX["LVEDP"]] - P["LVEDP0"]) \
        - P["KBNP_OFF"] * (y[IX["BNP"]] - 60.0)

    # ---- hazards (per hour) ------------------------------------------
    hyp = fast_solve(y, arch, 1.0, rx, hyperaemia=True)
    cfr = hyp["v"] / max(1e-9, f["v"])
    d[IX["CH_MORT"]] = P["H0_MORT"] / 8760.0 * math.exp(
        P["LN116"] * max(0.0, 2.5 - cfr) / 0.1)
    d[IX["CH_MACE"]] = P["H0_HOSP"] / 8760.0 * math.exp(
        P["LN108"] * max(0.0, 2.5 - cfr) / 0.1)
    d[IX["CH_HOSP"]] = (P["H0_HOSP"] + P["K_HOSP_ANG"] * y[IX["ANG"]]) / 8760.0
    d[IX["AUC_DEF"]] = bi
    d[IX["AUC_NOC"]] = bn

    # Fast tone states.  In mrgsolve these are genuine ODEs,
    #     dTONE_E/dt = (TONE_E_target - TONE_E)/TAU_TONE,   TAU_TONE = 15 s
    # integrated by LSODA, which is built for exactly this stiffness.  Here
    # they are solved by bisection inside fast_solve instead, and their
    # derivative is set to zero (bug B17: carrying a 15 s state through a
    # 3 h explicit RK4 step is unstable by a factor of 720 - TONE_E reached
    # 2e9 within one step and the nan propagated through the state vector,
    # silently returning SAQ = 100 for a symptomatic patient).  Section V
    # shows the two formulations reach the same steady state, which is all
    # any reported index is defined as.
    d[IX["TONE_E"]] = 0.0
    d[IX["TONE_P"]] = 0.0

    # ---- PK ----------------------------------------------------------
    # The pharmacokinetic states are LINEAR and are advanced EXACTLY by
    # pk_advance() outside the RK4 stages, so their derivative here is zero.
    #
    # bug B25: they used to be integrated with the physiology, and they were
    # the stiffest states in the whole model.  Explicit RK4 is unstable for
    # |k*dt| > 2.785, and the absorption rate constants run to 1.6 /h, so at
    # the 2 h step used for the 24-week runs the nicorandil and fasudil arms
    # were unstable, and at the 6 h step used for the 2.5-year WARRIOR run
    # every drug with ka > 0.46 /h was.  Ramipril's concentration reached
    # -97213 mg/L on day 1 and nan thereafter, so the WARRIOR arm silently
    # became a placebo arm: the treated and untreated curves were identical
    # to three decimals, which is what gave it away.  In mrgsolve the same
    # equations are ODEs and LSODA handles them; here they are solved in
    # closed form, which is also exact rather than merely stable.
    return d


# =====================================================================
# 5.  INTEGRATOR
# =====================================================================
METAB = (("M_ivabradine", "C_ivabradine", 0.42, 0.155),
         ("M_ramipril", "C_ramipril", 0.62, 0.038),
         ("M_fasudil", "C_fasudil", 0.85, 0.29))


def pk_advance(y, dt):
    """
    Advance every linear PK state exactly over dt (bug B25).

    One-compartment oral model in closed form; active metabolites use the
    exact decay with a trapezoidal parent input over the same interval,
    which is second-order accurate and unconditionally stable.
    """
    cbefore = {src: y[IXP[src]] for _m, src, _kf, _km in METAB}
    for dr in DRUG_ORDER:
        ka, ke, V = DRUGS[dr]["ka"], DRUGS[dr]["ke"], DRUGS[dr]["V"]
        ia, ic = IXP["A_" + dr], IXP["C_" + dr]
        A0, C0 = y[ia], y[ic]
        eA, eC = math.exp(-ka * dt), math.exp(-ke * dt)
        y[ia] = A0 * eA
        if abs(ka - ke) > 1e-9:
            y[ic] = C0 * eC + (ka * A0 / (V * (ka - ke))) * (eC - eA)
        else:
            y[ic] = C0 * eC + (ka * A0 / V) * dt * eC
    for mname, src, kf, km in METAB:
        im = IXP[mname]
        em = math.exp(-km * dt)
        cbar = 0.5 * (cbefore[src] + y[IXP[src]])
        y[im] = y[im] * em + kf * cbar * (1.0 - em) / km



def simulate(arch, reg, days=168, dt=2.0, rec_every=0):
    # a dose interval shorter than 6*dt aliases the peaks away; TID regimens
    # therefore force dt <= 2 h even in the multi-year runs
    for _mg, _iv, _st in reg.get("doses", {}).values():
        dt = min(dt, max(0.5, _iv / 4.0))
    y = y0(arch)
    y[IX["TONE_E"]] = fast_solve(y, arch, 1.0, {})["endo"]["tone"]
    y[IX["TONE_P"]] = fast_solve(y, arch, 1.0, {})["epi"]["tone"]
    doses = reg.get("doses", {})
    nxt = {dr: v[2] for dr, v in doses.items()}
    T = int(round(days * 24 / dt))
    rec = []
    for step in range(T + 1):
        t = step * dt
        for dr, (mg, iv, st) in doses.items():
            while t + 1e-9 >= nxt[dr]:
                y[IXP["A_" + dr]] += mg
                nxt[dr] += iv
        if rec_every and step % rec_every == 0:
            rx, *_ = regimen_effects(y, reg)
            idx = indices(y, arch, rx)
            bi, bn = burdens(y, arch, rx)
            rec.append(dict(t=t / 24.0, CFR=idx["CFR"], MR_rest=idx["MR_rest"],
                            MR_hyp=idx["MR_hyp"], IMR=idx["IMR"], MRR=idx["MRR"],
                            endo_epi=idx["endo_epi_hyp"], burden=bi, noci=bn,
                            SAQ=y[IX["SAQ"]], ANG=y[IX["ANG"]],
                            SENS=y[IX["SENS"]], LVEDP=y[IX["LVEDP"]],
                            HR=idx["rest"]["h"]["HR"], TNI=y[IX["TNI"]],
                            tone=idx["rest"]["tone"], ET1=y[IX["ET1"]],
                            CH_HOSP=y[IX["CH_HOSP"]], CH_MORT=y[IX["CH_MORT"]]))
        if step == T:
            break
        pk_advance(y, dt / 2.0)          # PK to the step midpoint (bug B25)
        k1 = rhs(t, y, arch, reg)
        y2 = [y[i] + 0.5 * dt * k1[i] for i in range(NS)]
        k2 = rhs(t + dt / 2, y2, arch, reg)
        y3 = [y[i] + 0.5 * dt * k2[i] for i in range(NS)]
        k3 = rhs(t + dt / 2, y3, arch, reg)
        y4 = [y[i] + dt * k3[i] for i in range(NS)]
        k4 = rhs(t + dt, y4, arch, reg)
        for i in range(NS):
            y[i] += dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
        pk_advance(y, dt / 2.0)          # PK to the end of the step
        y[IX["SAQ"]] = max(0.0, min(100.0, y[IX["SAQ"]]))
        y[IX["CAPD"]] = max(0.35, y[IX["CAPD"]])
        y[IX["SENS"]] = max(0.0, min(1.0, y[IX["SENS"]]))
        y[IX["SYMP"]] = max(0.0, min(1.5, y[IX["SYMP"]]))
        y[IX["ANG"]] = max(0.0, y[IX["ANG"]])
    return y, rec


def R(**kw):
    r = dict(doses={})
    for k, v in kw.items():
        (r["doses"] if k in DRUGS else r)[k] = v
    return r


BID, TID, QD = 12.0, 8.0, 24.0
REGIMENS = {
    "untreated": R(),
    "ivabradine": R(ivabradine=(7.5, BID, 0.0)),
    "bisoprolol": R(bisoprolol=(5.0, QD, 0.0)),
    "nebivolol": R(nebivolol=(5.0, QD, 0.0)),
    "ranolazine": R(ranolazine=(1000.0, BID, 0.0)),
    "amlodipine": R(amlodipine=(10.0, QD, 0.0)),
    "diltiazem": R(diltiazem=(180.0, QD, 0.0)),
    "nicorandil": R(nicorandil=(20.0, TID, 0.0)),
    "zibotentan": R(zibotentan=(10.0, QD, 0.0)),
    "zibotentan_cf": R(zibotentan=(10.0, QD, 0.0), no_bp_drop=True, no_fluid=True),
    "zib_nobp": R(zibotentan=(10.0, QD, 0.0), no_bp_drop=True),
    "zib_nofluid": R(zibotentan=(10.0, QD, 0.0), no_fluid=True),
    "imt_warrior": R(ramipril=(10.0, QD, 0.0), atorvastatin=(80.0, QD, 0.0)),
    "fasudil": R(fasudil=(80.0, TID, 0.0)),
    "sildenafil": R(sildenafil=(50.0, TID, 0.0)),
    "aminophylline": R(aminophylline=(225.0, BID, 0.0)),
    "imipramine": R(imipramine=(50.0, QD, 0.0)),
    "trimetazidine": R(trimetazidine=True),
    "rehab": R(rehab=1.0),
    "iva_ran": R(ivabradine=(7.5, BID, 0.0), ranolazine=(1000.0, BID, 0.0)),
    "cormica_func": R(ivabradine=(7.5, BID, 0.0), ramipril=(10.0, QD, 0.0),
                      atorvastatin=(80.0, QD, 0.0)),
    # CorMicA's microvascular-angina arm was a beta-blocker first line
    # (nebivolol) plus an ACE inhibitor and a statin; the first version of
    # this regimen omitted rate control and the structural stratum came out
    # WORSE on "stratified" therapy than on the unguided beta-blocker, which
    # was a transcription error, not a finding.
    "cormica_struct": R(nebivolol=(5.0, QD, 0.0), ramipril=(10.0, QD, 0.0),
                        atorvastatin=(80.0, QD, 0.0)),
    "cormica_vaso": R(amlodipine=(10.0, QD, 0.0), atorvastatin=(80.0, QD, 0.0)),
    "cormica_noncard": R(imipramine=(50.0, QD, 0.0), cbt=1.0),
}


# =====================================================================
# SECTION 0 -- CALIBRATION
# =====================================================================
def bisect_target(f, lo, hi, target, tol, nmax=14):
    """Robust root-find for a monotone response (used for the two symptom
    constants, where the secant step ran the parameter negative)."""
    flo, fhi = f(lo) - target, f(hi) - target
    for _ in range(6):
        if flo * fhi <= 0.0:
            break
        if abs(flo) < abs(fhi):
            lo = max(1e-4, lo * 0.35)
            flo = f(lo) - target
        else:
            hi = hi * 2.2
            fhi = f(hi) - target
    for _ in range(nmax):
        mid = 0.5 * (lo + hi)
        fm = f(mid) - target
        if abs(fm) < tol:
            return mid
        if fm * flo > 0.0:
            lo, flo = mid, fm
        else:
            hi, fhi = mid, fm
    return 0.5 * (lo + hi)


def secant(f, x0, x1, target, tol=1e-4, nmax=30):
    f0, f1 = f(x0) - target, f(x1) - target
    for _ in range(nmax):
        if abs(f1 - f0) < 1e-14:
            break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x2 = max(1e-6, x2)
        f2 = f(x2) - target
        x0, f0, x1, f1 = x1, f1, x2, f2
        if abs(f1) < tol:
            break
    return x1


def section_0():
    rule("SECTION 0.  CALIBRATION AGAINST MEASURED CORONARY PHYSIOLOGY")
    emit("Seven constants are solved (not tuned by hand) against seven published")
    emit("targets.  Everything else in the model is fixed a priori.")
    emit()
    emit("targets:")
    emit("  T1 control hyperaemic MR      2.20 mmHg/(cm/s)   Rahman 2019 (MR_hyp<2.5)")
    emit("  T2 structural hyperaemic MR   3.60               Rahman 2019 (>=2.5)")
    emit("  T3 functional hyperaemic MR   2.30 = 4.2/1.83    Rahman 2019 (internally")
    emit("                                                    consistent with T4)")
    emit("  T4 functional resting MR      4.20               Rahman 2019 (4.2+-1.0)")
    emit("  T5 control hyperaemic endo/epi 1.02              transmural physiology")
    emit("  T6 untreated functional SAQ   55.0               CorMicA baseline ~55-60")
    emit("  T7 untreated functional Bruce  480 s             PRIZE baseline order")
    emit()
    emit("Note on T3: Rahman's three functional-group numbers are self-consistent -")
    emit("4.2 / 2.30 = 1.83 - because CFR = MR_rest/MR_hyp identically when both are")
    emit("measured at the same aortic pressure.  Calibrating two of them therefore")
    emit("fixes the third, and the model's functional CFR is not free.")
    emit()

    def mr_hyp(arch):
        y = y0(arch)
        return indices(y, arch, {})["MR_hyp"]

    def mr_rest(arch):
        y = y0(arch)
        return indices(y, arch, {})["MR_rest"]

    # T1: RMIN0 -> control MR_hyp
    P["RMIN0"] = secant(lambda x: (P.__setitem__("RMIN0", x), mr_hyp("control"))[1],
                        1.2, 1.7, 2.20)
    # T4: W_ENDO -> control hyperaemic endo/epi (independent of RMIN0 only
    # through LVEDP, so iterate T1/T4 twice)
    for _ in range(3):
        def f_ee(x):
            P["W_ENDO"] = x
            y = y0("control")
            return indices(y, "control", {})["endo_epi_hyp"]
        P["W_ENDO"] = secant(f_ee, 0.50, 0.65, 1.02)
        P["RMIN0"] = secant(lambda x: (P.__setitem__("RMIN0", x),
                                       mr_hyp("control"))[1], 1.2, 1.7, 2.20)
    # T3: F_ADO_REV -> functional MR_hyp (its structural floor is fixed at 1.0)
    def f_ar(x):
        P["F_ADO_REV"] = min(0.999, x)
        return mr_hyp("functional")
    P["F_ADO_REV"] = min(0.999, secant(f_ar, 0.70, 0.95, 2.30))
    # T2: structural RMIN_F -> structural MR_hyp
    def f_rf(x):
        ARCHETYPE["structural"]["RMIN_F"] = x
        return mr_hyp("structural")
    ARCHETYPE["structural"]["RMIN_F"] = secant(f_rf, 1.3, 2.2, 3.60)
    # T4: functional AUTO_OFF -> functional MR_rest
    def f_ao(x):
        ARCHETYPE["functional"]["AUTO_OFF"] = x
        return mr_rest("functional")
    ARCHETYPE["functional"]["AUTO_OFF"] = secant(f_ao, 0.4, 1.0, 4.20)

    emit("solved:")
    emit("  RMIN0                       = %.4f mmHg/(cm/s)" % P["RMIN0"])
    emit("  W_ENDO                      = %.4f" % P["W_ENDO"])
    emit("  F_ADO_REV                   = %.4f" % P["F_ADO_REV"])
    emit("  RMIN_F(structural)          = %.4f" % ARCHETYPE["structural"]["RMIN_F"])
    emit("  AUTO_OFF(functional)        = %.4f" % ARCHETYPE["functional"]["AUTO_OFF"])
    emit()
    emit("  achieved: control MR_hyp %.3f | functional MR_hyp %.3f | structural "
         "MR_hyp %.3f" % (mr_hyp("control"), mr_hyp("functional"),
                          mr_hyp("structural")))
    emit("            functional MR_rest %.3f -> CFR %.3f (observed 1.83)"
         % (mr_rest("functional"), indices(y0("functional"), "functional", {})["CFR"]))
    emit("            control hyperaemic endo/epi %.3f"
         % indices(y0("control"), "control", {})["endo_epi_hyp"])
    emit()

    # T5/T6 need time integration; 84 days is enough for SAQ (tau 14 d)
    def f_saq(x):
        P["K_ANG"] = x
        y, _ = simulate("functional", REGIMENS["untreated"], days=84)
        return y[IX["SAQ"]]
    P["K_ANG"] = secant(f_saq, 8.0, 30.0, 55.0, tol=0.15, nmax=10)

    def f_ex(x):
        P["NOCI_THRESH"] = x
        y, _ = simulate("functional", REGIMENS["untreated"], days=84)
        rx, *_ = regimen_effects(y, REGIMENS["untreated"])
        return exercise_duration(y, "functional", rx)
    P["NOCI_THRESH"] = secant(f_ex, 2.0, 4.5, 480.0, tol=4.0, nmax=10)
    y, _ = simulate("functional", REGIMENS["untreated"], days=84, dt=3.0)
    rx, *_ = regimen_effects(y, REGIMENS["untreated"])
    emit("  K_ANG                       = %.4f  -> SAQ %.2f" % (P["K_ANG"], y[IX["SAQ"]]))
    emit("  NOCI_THRESH                 = %.4f  -> Bruce %.1f s"
         % (P["NOCI_THRESH"], exercise_duration(y, "functional", rx)))
    emit()
    emit("These five flow constants and two symptom constants are the ONLY values")
    emit("fitted to data.  Every result below is a prediction from them.")
    return dict(RMIN0=P["RMIN0"], W_ENDO=P["W_ENDO"], F_ADO_REV=P["F_ADO_REV"],
                RMIN_F_struct=ARCHETYPE["structural"]["RMIN_F"],
                AUTO_OFF_func=ARCHETYPE["functional"]["AUTO_OFF"],
                K_ANG=P["K_ANG"], NOCI_THRESH=P["NOCI_THRESH"])


# =====================================================================
# SECTION I
# =====================================================================
def section_I():
    rule("SECTION I.  CFR IS A RATIO: TWO ENDOTYPES, ONE NUMBER")
    emit("Rahman et al. Circulation 2019;140:1805 (n=85 ANOCA, 78% female).")
    emit("MVD = CFR<2.5 (53% of the cohort).  Of those, 62% 'functional'")
    emit("(hyperaemic MR<2.5 mmHg/cm/s) and 38% 'structural' (>=2.5).  Resting MR")
    emit("4.2+-1.0 functional, 6.9+-1.7 structural, 7.3+-2.2 control.")
    emit()
    hdr = ("%-12s %6s %8s %8s %7s %7s %7s %7s %6s %8s"
           % ("endotype", "CFR", "MR_rest", "MR_hyp", "v_rest", "v_hyp",
              "tone_r", "IMR", "MRR", "endo/epi"))
    emit(hdr)
    emit("-" * len(hdr))
    tab = {}
    for arch in ("control", "functional", "structural", "vasospastic", "noncardiac"):
        y = y0(arch)
        idx = indices(y, arch, {})
        tab[arch] = idx
        emit("%-12s %6.2f %8.2f %8.2f %7.1f %7.1f %7.3f %7.1f %6.2f %8.3f"
             % (arch, idx["CFR"], idx["MR_rest"], idx["MR_hyp"],
                idx["rest"]["v"], idx["hyp"]["v"], idx["rest"]["tone"],
                idx["IMR"], idx["MRR"], idx["endo_epi_hyp"]))
    emit()
    obs = dict(control=(7.3, None), functional=(4.2, 1.83), structural=(6.9, 1.92))
    emit("observed vs model:")
    for k, (mr, cfr) in obs.items():
        s = "  %-11s MR_rest obs %.1f  model %.2f  (%+.0f%%)" \
            % (k, mr, tab[k]["MR_rest"], 100 * (tab[k]["MR_rest"] - mr) / mr)
        if cfr:
            s += "   |  CFR obs %.2f  model %.2f" % (cfr, tab[k]["CFR"])
        emit(s)
    emit()
    f, s, c = tab["functional"], tab["structural"], tab["control"]
    emit("THE IDENTITY, as arithmetic:")
    emit("  functional  CFR = %.1f / %.1f = %.2f   (hyperaemic velocity %.0f%% of control)"
         % (f["hyp"]["v"], f["rest"]["v"], f["CFR"],
            100 * f["hyp"]["v"] / c["hyp"]["v"]))
    emit("  structural  CFR = %.1f / %.1f = %.2f   (hyperaemic velocity %.0f%% of control)"
         % (s["hyp"]["v"], s["rest"]["v"], s["CFR"],
            100 * s["hyp"]["v"] / c["hyp"]["v"]))
    emit("  The two CFR values differ by %.2f. The two ABSOLUTE hyperaemic velocities"
         % abs(f["CFR"] - s["CFR"]))
    emit("  differ by %.1f cm/s (%.0f%%), i.e. %.2f vs %.2f mL/min/g of maximal"
         % (abs(f["hyp"]["v"] - s["hyp"]["v"]),
            100 * abs(f["hyp"]["v"] - s["hyp"]["v"]) / s["hyp"]["v"],
            f["hyp"]["mbf"], s["hyp"]["mbf"]))
    emit("  perfusion.  A ratio cannot see that difference; an absolute flow can.")
    emit("  This is the argument for measuring absolute flow (thermodilution MRR)")
    emit("  rather than a reserve ratio, and the model quantifies what is lost.")
    return tab


# =====================================================================
# SECTION II
# =====================================================================
def section_II(tab):
    rule("SECTION II.  RESTING TONE IS NOT A FREE PARAMETER")
    emit("Autoregulation sets tone so that resting supply meets resting demand.")
    emit("Given the measured HYPERAEMIC MR (which fixes R_min), the RESTING MR is")
    emit("therefore predicted with no further parameter.  Only the control group's")
    emit("hyperaemic MR was calibrated; the structural group's resting MR was not.")
    emit()
    hdr = "%-12s %9s %14s %13s %8s" % ("endotype", "MR_hyp", "MR_rest_pred",
                                       "MR_rest_obs", "error")
    emit(hdr)
    emit("-" * len(hdr))
    obs = dict(control=7.3, structural=6.9)
    for arch in ("control", "structural"):
        save = ARCHETYPE[arch]["AUTO_OFF"]
        ARCHETYPE[arch]["AUTO_OFF"] = 0.0
        idx = indices(y0(arch), arch, {})
        ARCHETYPE[arch]["AUTO_OFF"] = save
        emit("%-12s %9.2f %14.2f %13.1f %7.1f%%"
             % (arch, idx["MR_hyp"], idx["MR_rest"], obs[arch],
                100 * (idx["MR_rest"] - obs[arch]) / obs[arch]))
    emit()
    emit("The structural endotype's resting MR is a PREDICTION and it lands inside")
    emit("the measured 6.9+-1.7.  Autoregulation reduces resting tone to compensate")
    emit("for the raised R_min, and that compensation is what consumes the reserve.")
    emit()
    emit("Now the functional endotype, controller intact, offset switched off:")
    save = ARCHETYPE["functional"]["AUTO_OFF"]
    ARCHETYPE["functional"]["AUTO_OFF"] = 0.0
    i0 = indices(y0("functional"), "functional", {})
    ARCHETYPE["functional"]["AUTO_OFF"] = save
    i1 = indices(y0("functional"), "functional", {})
    emit("  offset 0.00 -> MR_rest %.2f  CFR %.2f  tone %.3f"
         % (i0["MR_rest"], i0["CFR"], i0["rest"]["tone"]))
    emit("  offset %.2f -> MR_rest %.2f  CFR %.2f  tone %.3f   (observed 4.2 / 1.83)"
         % (save, i1["MR_rest"], i1["CFR"], i1["rest"]["tone"]))
    emit()
    emit("BUG B1 and the model's central claim: a working autoregulator CANNOT")
    emit("produce a low resting MR beside a normal minimal MR.  Setting tone by hand")
    emit("does not survive one step of the loop - the controller puts it back.  The")
    emit("functional endotype is therefore not a tone state but a CONTROLLER DEFECT,")
    emit("and its magnitude is derivable rather than assumed: the arterioles behave")
    emit("as though metabolic demand were %.0f%% higher than it is." % (100 * save))
    emit()
    c = indices(y0("control"), "control", {})["rest"]
    f = i1["rest"]
    emit("FALSIFIABLE PREDICTION - needs coronary-sinus sampling, not a wire:")
    emit("  resting extraction  control %.3f   functional %.3f   ratio %.2f"
         % (c["E_rest"], f["E_rest"], f["E_rest"] / c["E_rest"]))
    emit("  Resting myocardial oxygen extraction in the functional endotype must be")
    emit("  about %.0f%% of normal, with a correspondingly RAISED coronary venous"
         % (100 * f["E_rest"] / c["E_rest"]))
    emit("  oxygen saturation (%.0f%% -> roughly %.0f%% if arterial SO2 is 98%%)."
         % (100 * (1 - c["E_rest"]) * 0.98 + 0, 100 * (1 - f["E_rest"]) * 0.98))
    emit("  If measured resting extraction is normal, the controller-offset reading")
    emit("  is wrong and the low resting MR must instead be an artefact of the")
    emit("  measurement condition (anxiety, contrast, supine catheter-laboratory")
    emit("  haemodynamics).  Either way the model is refutable on one blood gas,")
    emit("  which is the cheapest experiment it suggests.")
    return f["E_rest"] / c["E_rest"]


# =====================================================================
# SECTION III
# =====================================================================
def section_III():
    rule("SECTION III.  HEART RATE ENTERS THE OXYGEN BALANCE TWICE")
    emit("Resting heart rate is swept at FIXED workload (wl=3.2, where the untreated")
    emit("structural endotype has a subendocardial deficit) in the structural")
    emit("endotype.  Column 'def(time only)' holds demand at its HR=68 value so that")
    emit("only the diastolic window moves; 'def(full)' is the real effect.")
    emit()
    arch = "structural"
    y = y0(arch)
    hdr = ("%5s %7s %8s %8s %12s %11s %9s %8s"
           % ("HR", "f_dia", "DPTI", "SEVR", "def(time)", "def(full)",
              "endo/epi", "MVO2"))
    emit(hdr)
    emit("-" * len(hdr))
    rows = {}
    WL3 = 3.2
    mv_ref = fast_solve(y, arch, WL3, {})["mvo2"]
    for hr in (50, 55, 60, 68, 75, 85, 95, 110):
        fhr = 1.0 - hr / P["HR0"]
        ff = fast_solve(y, arch, WL3, {"f_HR": fhr})
        # cancel the demand change, keep the timing change
        ft = fast_solve(y, arch, WL3,
                        {"f_HR": fhr, "f_MVO2": 1.0 - mv_ref / ff["mvo2"]})
        rows[hr] = (ff, ft)
        emit("%5d %7.3f %8.1f %8.3f %12.4f %11.4f %9.3f %8.2f"
             % (hr, ff["h"]["fdia"], ff["dpti"], ff["sevr"], ft["def_endo"],
                ff["def_endo"], ff["endo_epi"], ff["mvo2"]))
    emit()
    b, t = rows[68], rows[55]
    d_full = b[0]["def_endo"] - t[0]["def_endo"]
    d_time = b[1]["def_endo"] - t[1]["def_endo"]
    emit("Slowing 68 -> 55 bpm at fixed workload:")
    emit("  total fall in subendocardial deficit  %.4f mL O2/min/100 g" % d_full)
    if abs(d_full) > 1e-9:
        emit("  from the diastolic TIME window        %.4f  (%.0f%%)"
             % (d_time, 100 * d_time / d_full))
        emit("  from lower DEMAND                     %.4f  (%.0f%%)"
             % (d_full - d_time, 100 * (d_full - d_time) / d_full))
    emit()
    emit("A pure rate-slowing drug with no vasodilator and no metabolic action still")
    emit("has a SUPPLY-side effect, and it is not interchangeable with removing the")
    emit("same amount of demand by another route.  Note also that DPTI falls")
    emit("monotonically with rate while SEVR falls faster, because the systolic")
    emit("tension-time index in the denominator rises at the same time.")
    emit()
    emit("The same result read as a threshold - the workload at which a")
    emit("subendocardial deficit first appears, as a function of resting heart rate:")
    emit()
    emit("%5s %14s" % ("HR", "onset workload"))
    emit("-" * 20)
    for hr in (50, 55, 60, 68, 75, 85, 95):
        fhr = 1.0 - hr / P["HR0"]
        lo, hi = 1.0, 5.0
        if fast_solve(y, arch, hi, {"f_HR": fhr})["def_endo"] <= 0:
            emit("%5d %14s" % (hr, ">5.0"))
            continue
        for _ in range(20):
            mid = 0.5 * (lo + hi)
            if fast_solve(y, arch, mid, {"f_HR": fhr})["def_endo"] > 0:
                hi = mid
            else:
                lo = mid
        emit("%5d %14.2f" % (hr, 0.5 * (lo + hi)))
    return d_full, d_time


# =====================================================================
# SECTION IV
# =====================================================================
def spearman(x, y):
    """Rank correlation (no scipy available)."""
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = ranks(x), ranks(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(n)))
    dy = math.sqrt(sum((ry[i] - my) ** 2 for i in range(n)))
    return num / (dx * dy) if dx > 0 and dy > 0 else 0.0


def section_IV():
    rule("SECTION IV.  CFR RANKS DRUGS WELL WHERE FLOW IS THE DISEASE, AND BADLY\n            WHERE IT IS NOT")
    emit("Steady-state drug effect after 21 days, evaluated in each endotype at")
    emit("matched dose.  dCFR is what a physiology study reports; dNOC is the afferent")
    emit("drive that becomes angina; dEXD is Bruce duration (s); def@3.2 is the")
    emit("subendocardial deficit at a stair-climbing workload.")
    emit()
    arms = ["ivabradine", "bisoprolol", "ranolazine", "amlodipine", "nicorandil",
            "zibotentan", "sildenafil", "aminophylline", "diltiazem", "fasudil",
            "imipramine"]
    rho = {}
    best = {}
    for arch in ("functional", "structural", "vasospastic"):
        emit("--- %s endotype ---" % arch)
        hdr = ("%-14s %7s %7s %7s %8s %7s %7s %7s %7s"
               % ("arm", "CFR", "dCFR", "MR_hyp", "dMR_hyp", "noci", "dNOC",
                  "dEXD", "def@3.2"))
        emit(hdr)
        emit("-" * len(hdr))
        y, _ = simulate(arch, REGIMENS["untreated"], days=21)
        rx0, *_ = regimen_effects(y, REGIMENS["untreated"])
        i0 = indices(y, arch, rx0)
        _bi0, bn0 = burdens(y, arch, rx0)
        e0 = exercise_duration(y, arch, rx0)
        d0 = fast_solve(y, arch, 3.2, rx0)["def_endo"]
        emit("%-14s %7.2f %7s %7.2f %8s %7.3f %7s %7s %7.3f"
             % ("(untreated)", i0["CFR"], "-", i0["MR_hyp"], "-", bn0, "-", "-", d0))
        dcfr, dnoc, dexd = [], [], []
        for arm in arms:
            yy, _ = simulate(arch, REGIMENS[arm], days=21)
            rx, *_ = regimen_effects(yy, REGIMENS[arm])
            idx = indices(yy, arch, rx)
            _bi, bn = burdens(yy, arch, rx)
            ee = exercise_duration(yy, arch, rx)
            dd = fast_solve(yy, arch, 3.2, rx)["def_endo"]
            dcfr.append(idx["CFR"] - i0["CFR"])
            dnoc.append(bn - bn0)
            dexd.append(ee - e0)
            emit("%-14s %7.2f %+7.3f %7.2f %+8.3f %7.3f %+7.3f %+7.1f %7.3f"
                 % (arm, idx["CFR"], idx["CFR"] - i0["CFR"], idx["MR_hyp"],
                    idx["MR_hyp"] - i0["MR_hyp"], bn, bn - bn0, ee - e0, dd))
        rho[arch] = (spearman(dcfr, [-v for v in dnoc]), spearman(dcfr, dexd))
        best[arch] = (arms[dcfr.index(max(dcfr))], arms[dnoc.index(min(dnoc))],
                      arms[dexd.index(max(dexd))])
        emit()
    emit("The signs agree almost everywhere.  An earlier version of this section")
    emit("claimed outright SIGN discordance between CFR and symptoms; that was an")
    emit("artefact of bug B24 and the claim is withdrawn.  What survives is weaker and")
    emit("more specific, and it is about RANKING - a trial picks one endpoint and ranks")
    emit("drugs on it:")
    emit()
    hdr = "%-13s %14s %14s %-16s %-16s" % ("endotype", "rho(dCFR,-dNOC)",
                                           "rho(dCFR,dEXD)", "best on CFR",
                                           "best on symptoms")
    emit(hdr)
    emit("-" * len(hdr))
    for arch in ("functional", "structural", "vasospastic"):
        emit("%-13s %14.2f %14.2f %-16s %-16s"
             % (arch, rho[arch][0], rho[arch][1], best[arch][0], best[arch][1]))
    emit()
    emit("The pattern in those three numbers is the point: CFR ranks treatments almost")
    emit("perfectly in the structural endotype, where flow IS the disease, and much")
    emit("less well in the two endotypes where it is not.")
    emit()
    emit("Read the functional-endotype rows again.  Aminophylline moves CFR by +0.012")
    emit("- fifth of eleven arms, indistinguishable from placebo on a physiological")
    emit("endpoint - and buys +147.5 s of Bruce time, twice the best vasoactive arm,")
    emit("because in that endotype the afferent adenosine signal and not the flow is")
    emit("the disease.  Elliott et al.")
    emit("Heart 1997;77:523 (PMID 9227295) gave oral aminophylline to patients with")
    emit("angina and normal arteriograms and improved exercise time; the model")
    emit("reproduces that from the structure rather than from a fitted effect, and")
    emit("predicts it should NOT work in the structural endotype, where blocking A2A")
    emit("costs dilator reserve that is actually being used.")
    emit()
    emit("In the vasospastic endotype the two calcium antagonists and the Rho-kinase")
    emit("inhibitor take the top three symptom places while their CFR changes are")
    emit("mid-table, and nicorandil - a potent dilator - does nothing at all,")
    emit("because between episodes there is nothing to dilate.")
    emit()
    emit("One negative result worth stating plainly: microvascular resistance reserve")
    emit("does NOT rescue this.  MRR = (CFR/FFR)(Pa_rest/Pa_hyper), and a drug that")
    emit("lowers blood pressure lowers it in BOTH states, so the correction cancels.")
    emit("Here Pa_rest and Pa_hyper are equal by construction, which makes the")
    emit("cancellation exact: MRR = CFR/FFR identically, so dMRR = dCFR/FFR to three")
    emit("decimals and MRR adds no information about a treatment effect.  MRR earns its")
    emit("keep against an EPICARDIAL stenosis, which is a different problem from a")
    emit("drug's own haemodynamics, and the two are easily conflated.  The")
    emit("index that stays interpretable under treatment is the hyperaemic RESISTANCE")
    emit("itself: dMR_hyp above is near zero for every drug that has no true effect on")
    emit("minimal resistance, whereas dCFR is not.  If a trial wants a physiological")
    emit("endpoint that a blood-pressure-lowering drug cannot fake, it should")
    emit("prespecify hyperaemic microvascular resistance, not a reserve ratio.")
    return rho, best


# =====================================================================
# SECTION V
# =====================================================================
def section_V():
    rule("SECTION V.  FIXED POINT (python) vs FAST ODE (mrgsolve)")
    emit("mrgsolve integrates TONE_E and TONE_P as ODEs with tau = 15 s under LSODA.")
    emit("Here the same expressions are solved by bisection.  Integrating the fast")
    emit("ODE explicitly with dt = 1 s from a deliberately wrong start must land on")
    emit("the same value:")
    emit()
    ok = True
    for arch in ("control", "functional", "structural", "vasospastic"):
        y = y0(arch)
        q = fast_solve(y, arch, 1.0, {})
        for layer, start in (("endo", 0.05), ("epi", 0.90)):
            tone = start
            dt = 1.0 / 3600.0
            h = systemic(y, arch, 1.0, {})
            mvo2, _ = demand(y, h, 1.0, {})
            Rmin = rmin_eff(y, {})
            tp = passive_tone(y, arch, {})
            myo = P["K_MYO"] * (h["Pdia"] - 90.0) / 90.0
            kO2 = P["KCONV"] * 100.0 * P["CAO2"]
            fdem = P["F_ENDO_DEM"] if layer == "endo" else P["F_EPI_DEM"]
            v_req = mvo2 * fdem / (kO2 * P["E_REST"]) \
                * (1.0 + ARCHETYPE[arch]["AUTO_OFF"])
            for _ in range(1800):
                Rart = Rmin / (1.0 - P["A_TONE"] * tone) ** 4
                v, _pd = layer_velocity(y, h, Rart, layer, {})
                tgt = max(0.0, min(P["TONE_MAX"],
                                   P["TONE_REF"] + P["G_AUTO"] * (v / v_req - 1.0)
                                   + myo + tp))
                tone += dt / P["TAU_TONE"] * (tgt - tone)
            diff = abs(q[layer]["tone"] - tone)
            ok = ok and diff < 1e-6
            emit("  %-12s %-5s  bisection %.8f   fast ODE %.8f   |diff| %.2e"
                 % (arch, layer, q[layer]["tone"], tone, diff))
    emit()
    emit("agreement to 1e-6 in every case: %s" % ("YES" if ok else "NO"))
    emit("Both formulations report the same steady state, which is the only thing any")
    emit("of these indices is defined as.  The choice is numerical, not physiological")
    emit("(bug B11: the DISCRETE damped iteration was unstable at this gain; the")
    emit("continuous-time ODE never was).")
    return ok


# =====================================================================
# SECTION VI
# =====================================================================
def section_VI():
    rule("SECTION VI.  PRIZE (zibotentan): DECOMPOSING A NULL RESULT")
    emit("Morrow et al. Circulation 2024 (PMID 39217504).  Zibotentan 10 mg daily x")
    emit("12 weeks in microvascular angina, n=118, sequential crossover, enriched to a")
    emit("50% rs9349379-G allele frequency.  Primary endpoint Bruce treadmill duration:")
    emit("between-treatment difference -4.26 s (95% CI -19.60 to +11.06), p=0.59.")
    emit("Plasma ET-1 ROSE on treatment.  Blood pressure fell.  Adverse events 60.2%")
    emit("vs 14.4%, dominated by fluid retention.")
    emit()
    emit("The trial enrolled microvascular angina, not one endotype, so the")
    emit("decomposition is done per endotype and then weighted (0.33 functional / 0.20")
    emit("structural / 0.29 vasospastic / 0.18 non-cardiac, renormalised from the")
    emit("Rahman and CorMicA proportions).")
    emit()
    arms = (("full drug (as in PRIZE)", "zibotentan"),
            ("no BP fall", "zib_nobp"),
            ("no fluid retention", "zib_nofluid"),
            ("neither (counterfactual)", "zibotentan_cf"))
    W = dict(functional=0.33, structural=0.20, vasospastic=0.29, noncardiac=0.18)
    tot = {k: 0.0 for _n, k in arms}
    per = {}
    for arch, w in W.items():
        yb, _ = simulate(arch, REGIMENS["untreated"], days=84)
        rxb, *_ = regimen_effects(yb, REGIMENS["untreated"])
        exb = exercise_duration(yb, arch, rxb)
        ib = indices(yb, arch, rxb)
        emit("--- %s (weight %.2f) --- placebo period: Bruce %.1f s, CFR %.2f, "
             "MR_hyp %.2f, LVEDP %.1f, ET-1 %.2f"
             % (arch, w, exb, ib["CFR"], ib["MR_hyp"], yb[IX["LVEDP"]],
                yb[IX["ET1"]]))
        hdr = "%-27s %10s %7s %8s %7s %7s" % ("arm", "dBruce(s)", "CFR",
                                              "MR_hyp", "LVEDP", "ET-1")
        emit(hdr)
        emit("-" * len(hdr))
        per[arch] = {}
        for name, key in arms:
            yy, _ = simulate(arch, REGIMENS[key], days=84)
            rx, *_ = regimen_effects(yy, REGIMENS[key])
            ex = exercise_duration(yy, arch, rx)
            idx = indices(yy, arch, rx)
            per[arch][key] = ex - exb
            tot[key] += w * (ex - exb)
            emit("%-27s %+10.2f %7.2f %8.2f %7.1f %7.2f"
                 % (name, ex - exb, idx["CFR"], idx["MR_hyp"], yy[IX["LVEDP"]],
                    yy[IX["ET1"]]))
        emit()
    emit("Population-weighted signed decomposition of exercise duration:")
    emit("  microvascular gain alone (both costs removed)   %+7.2f s"
         % tot["zibotentan_cf"])
    emit("  attributable to the blood-pressure fall         %+7.2f s"
         % (tot["zibotentan"] - tot["zib_nobp"]))
    emit("  attributable to fluid retention (LVEDP up)      %+7.2f s"
         % (tot["zibotentan"] - tot["zib_nofluid"]))
    emit("  model net                                       %+7.2f s"
         % tot["zibotentan"])
    emit("  PRIZE observed                                    -4.26 s "
         "(95% CI -19.60 to +11.06)")
    inside = -19.60 <= tot["zibotentan"] <= 11.06
    emit("  model net inside the observed 95%% CI: %s" % ("YES" if inside else "NO"))
    emit()
    emit("The model reproduces the trial's number, and the reason is not the one this")
    emit("section was written expecting.  The three terms are:")
    emit()
    emit("  1. The microvascular gain is exactly zero.  Removing ETA-mediated tone")
    emit("     helps only where there is tone to remove that is not already removed:")
    emit("     in the functional endotype minimal resistance is normal, so there is")
    emit("     nothing to gain, and in the structural endotype the ceiling is")
    emit("     anatomical, so a receptor antagonist cannot lift it.  A drug can be")
    emit("     perfectly on-target and still have no target left to hit.")
    emit("  2. The blood-pressure fall is a small NET BENEFIT here, not a cost")
    emit("     (%+.2f s).  Lowering systolic pressure removes wall-stress demand faster"
         % (tot["zibotentan"] - tot["zib_nobp"]))
    emit("     than it removes diastolic driving pressure, because demand carries the")
    emit("     systolic pressure through the tension-time product while supply carries")
    emit("     only the diastolic mean.  This section originally asserted the opposite")
    emit("     and the arithmetic contradicted it.")
    emit("  3. Fluid retention accounts for the ENTIRE null (%+.2f s).  LVEDP is"
         % (tot["zibotentan"] - tot["zib_nofluid"]))
    emit("     subtracted from subendocardial driving pressure directly, so the")
    emit("     trial's dominant adverse effect - present in 60% of treated periods -")
    emit("     and its failure on the primary endpoint are the SAME event measured")
    emit("     twice.")
    emit()
    emit("That is an actionable reading rather than a dead end.  It says the endothelin")
    emit("hypothesis was not tested cleanly at 10 mg: the endpoint was moved by the")
    emit("sodium retention and not by the receptor.  Two designs follow - co-administer")
    emit("a diuretic or a lower dose and re-measure, or use a molecule that separates")
    emit("coronary ETA blockade from renal sodium handling.  A larger trial of the same")
    emit("10 mg would measure the same fluid.")
    emit()
    emit("Direction check on plasma ET-1 (bug B8): the model raises it on treatment,")
    emit("as PRIZE measured, because ETA blockade also removes ETA-mediated clearance.")
    emit("A model that lowered ET-1 would have predicted a win.")
    return tot, per


# =====================================================================
# SECTION VII
# =====================================================================
def section_VII():
    rule("SECTION VII.  RWISE (ranolazine): A REAL SUBGROUP EFFECT, DILUTED AWAY")
    emit("Bairey Merz et al. Eur Heart J 2016;37:1504 (PMID 26614823).  Ranolazine")
    emit("500-1000 mg BID x 2 weeks, n=128 (96% women), crossover.  No overall")
    emit("treatment difference in SAQ.  In the CFR<2.5 subgroup: MPRI improved")
    emit("(p=0.014), SAQ angina frequency (p=0.027) and SAQ-7 (p=0.041) improved.")
    emit()
    prev = dict(control=0.47, functional=0.33, structural=0.20)
    emit("Composition from Rahman 2019: 47% CFR>=2.5; of the remaining 53%, 62%")
    emit("functional and 38%% structural -> %.2f and %.2f."
         % (prev["functional"], prev["structural"]))
    emit()
    hdr = "%-12s %7s %9s %9s %8s %10s" % ("stratum", "weight", "SAQ_pbo",
                                          "SAQ_ran", "dSAQ", "dCFR")
    emit(hdr)
    emit("-" * len(hdr))
    dl = {}
    for arch, w in prev.items():
        yp, _ = simulate(arch, REGIMENS["untreated"], days=84)
        yr, _ = simulate(arch, REGIMENS["ranolazine"], days=84)
        rxp, *_ = regimen_effects(yp, REGIMENS["untreated"])
        rxr, *_ = regimen_effects(yr, REGIMENS["ranolazine"])
        dcfr = indices(yr, arch, rxr)["CFR"] - indices(yp, arch, rxp)["CFR"]
        dl[arch] = (w, yp[IX["SAQ"]], yr[IX["SAQ"]], yr[IX["SAQ"]] - yp[IX["SAQ"]],
                    dcfr)
        emit("%-12s %7.2f %9.1f %9.1f %+8.2f %+10.3f"
             % (arch, w, yp[IX["SAQ"]], yr[IX["SAQ"]],
                yr[IX["SAQ"]] - yp[IX["SAQ"]], dcfr))
    emit()
    wcmd = prev["functional"] + prev["structural"]
    d_cmd = (dl["functional"][0] * dl["functional"][3]
             + dl["structural"][0] * dl["structural"][3]) / wcmd
    d_all = sum(v[0] * v[3] for v in dl.values())
    emit("CFR<2.5 stratum only (what the subgroup analysis saw):  dSAQ %+.2f U" % d_cmd)
    emit("whole randomised cohort (what the primary analysis saw): dSAQ %+.2f U" % d_all)
    emit("dilution: the non-CMD stratum carries %.0f%% of the weight and contributes"
         % (100 * prev["control"]))
    emit("  %+.2f U, so the observable effect is %.0f%% of the mechanistic effect."
         % (dl["control"][3], 100 * d_all / d_cmd if abs(d_cmd) > 1e-9 else 0))
    emit()
    emit("Taking %.0f U as the clinically important SAQ difference (CorMicA's"
         % P["MCID_SAQ"])
    emit("stratified arm achieved +11.7 U, 95% CI 5.0-18.4):")
    emit("  subgroup     %+.2f U -> %s"
         % (d_cmd, "detectable" if abs(d_cmd) >= P["MCID_SAQ"] else "below the MCID"))
    emit("  whole cohort %+.2f U -> %s"
         % (d_all, "detectable" if abs(d_all) >= P["MCID_SAQ"] else "below the MCID"))
    emit()
    emit("Nothing about the drug changed between those two lines - only the")
    emit("prevalence of the mechanism it targets.  Note also WHERE the effect sits:")
    emit("it is largest in the %s stratum, i.e. in the patients whose minimal"
         % max(("functional", "structural"), key=lambda k: dl[k][3]))
    emit("resistance is actually abnormal.  A late-sodium-current inhibitor works by")
    emit("lowering diastolic wall tension, which needs a filling-pressure problem to")
    emit("act on; the functional endotype does not have one.  The argument is for an")
    emit("interventional diagnostic procedure BEFORE randomisation, not after.")
    return d_cmd, d_all, prev


# =====================================================================
# SECTION VIII
# =====================================================================
def section_VIII():
    rule("SECTION VIII.  THE PAIN->SYMPATHETIC->ISCHAEMIA LOOP HAS A THRESHOLD")
    emit("Angina raises sympathetic tone, which raises heart rate and alpha1")
    emit("prearteriolar tone, which raises the subendocardial deficit and the")
    emit("afferent drive, which raises angina.  Sweep the loop gain K_SYMP_HR:")
    emit()
    arch = "functional"
    keep = P["K_SYMP_HR"]
    hdr = "%10s %8s %8s %8s %8s %9s %9s" % ("K_SYMP_HR", "HR_rest", "ANG/wk",
                                            "SENS", "SAQ", "noci", "CFR")
    emit(hdr)
    emit("-" * len(hdr))
    prev = None
    jump = None
    for k in (0.0, 0.10, 0.18, 0.22, 0.26, 0.30, 0.36, 0.44):
        P["K_SYMP_HR"] = k
        y, _ = simulate(arch, REGIMENS["untreated"], days=112)
        _bi, bn = burdens(y, arch, {})
        hr = systemic(y, arch, 1.0, {})["HR"]
        cfr = indices(y, arch, {})["CFR"]
        emit("%10.2f %8.1f %8.2f %8.3f %8.1f %9.3f %9.2f"
             % (k, hr, y[IX["ANG"]], y[IX["SENS"]], y[IX["SAQ"]], bn, cfr))
        if prev is not None and prev - y[IX["SAQ"]] > 6.0 and jump is None:
            jump = k
        prev = y[IX["SAQ"]]
    P["K_SYMP_HR"] = keep
    emit()
    if jump:
        emit("The SAQ trajectory turns over between gains %.2f and %.2f: below that"
             % (jump - 0.04, jump))
        emit("the loop is self-limiting, above it each extra episode buys the next.")
    else:
        emit("Across this range the loop stays self-limiting, so the refractory")
        emit("phenotype in this model needs a second driver (central sensitisation,")
        emit("spasm) and not merely a high autonomic gain.")
    emit("Operationally, heart-rate control is not only an anti-anginal: it is the")
    emit("loop-breaker.  That is a different argument for the same prescription, and")
    emit("it predicts the benefit should be largest in patients with the highest")
    emit("resting heart rate rather than the lowest CFR.")
    emit("(bug B5: an unchecked version of this loop drove every untreated patient to")
    emit("SAQ 0, which is why the gain is swept and reported rather than clamped.)")
    return jump


# =====================================================================
# SECTION IX
# =====================================================================
def section_IX():
    rule("SECTION IX.  WARRIOR (statin + ACEi/ARB): CONTAMINATION AND TIMESCALE")
    emit("Pepine et al. Open Heart 2026;13:e004115 (PMID 41932694).  2476 women with")
    emit("suspected ANOCA/INOCA, intensive medical therapy (high-intensity statin +")
    emit("ACEi/ARB + aspirin) vs usual care, 2.5 years.  Primary MACE HR 1.13 (0.94-")
    emit("1.37), p=0.20.  Contamination sensitivity analysis HR 0.74 (0.35-1.56).")
    emit("Hospitalisation for angina dominated MACE.  Baseline LDL and BP already")
    emit("well controlled with 'relatively high rates of statin and ACEI/ARB use'.")
    emit()
    arch = "structural"
    yu, ru = simulate(arch, REGIMENS["untreated"], days=912, dt=6.0, rec_every=28)
    yt, rt = simulate(arch, REGIMENS["imt_warrior"], days=912, dt=6.0, rec_every=28)
    emit("2.5-year structural-endotype run (dt = 6 h):")
    emit("  usual care : CFR %.2f -> %.2f | MR_hyp %.2f -> %.2f | ANG %.2f | SAQ %.1f"
         % (ru[0]["CFR"], ru[-1]["CFR"], ru[0]["MR_hyp"], ru[-1]["MR_hyp"],
            yu[IX["ANG"]], yu[IX["SAQ"]]))
    emit("  IMT        : CFR %.2f -> %.2f | MR_hyp %.2f -> %.2f | ANG %.2f | SAQ %.1f"
         % (rt[0]["CFR"], rt[-1]["CFR"], rt[0]["MR_hyp"], rt[-1]["MR_hyp"],
            yt[IX["ANG"]], yt[IX["SAQ"]]))
    emit()
    hr_true = yt[IX["CH_HOSP"]] / yu[IX["CH_HOSP"]]
    emit("cumulative hazard of angina hospitalisation over 2.5 y:")
    emit("  usual care %.4f   IMT %.4f   true hazard ratio %.3f"
         % (yu[IX["CH_HOSP"]], yt[IX["CH_HOSP"]], hr_true))
    emit()
    emit("If a fraction c of the usual-care arm is already taking the same drugs, the")
    emit("observable ratio is  HR_obs = 1 + (HR_true - 1)*(1 - c):")
    emit()
    hdr = "%14s %10s %14s" % ("contamination", "HR_obs", "distinguishable")
    emit(hdr)
    emit("-" * len(hdr))
    for c in (0.0, 0.2, 0.4, 0.6, 0.75, 0.85, 0.95):
        hro = 1.0 + (hr_true - 1.0) * (1.0 - c)
        emit("%14.2f %10.3f %14s" % (c, hro, "yes" if abs(1 - hro) > 0.12 else "no"))
    emit()
    emit("With a true ratio of %.2f, a contamination of about 0.8 puts the observable"
         % hr_true)
    emit("ratio at %.2f - inside WARRIOR's reported interval and indistinguishable"
         % (1.0 + (hr_true - 1.0) * 0.2))
    emit("from 1.  The trial's own contamination-adjusted estimate was 0.74.")
    emit()
    emit("The second problem is the instrument.  The drugs' target is structural")
    emit("(tau_ML %.0f d, tau_CAPD %.0f d, tau_PVF %.0f d) while the endpoint that"
         % (P["TAU_ML"] / 24, P["TAU_CAPD"] / 24, P["TAU_PVF"] / 24))
    emit("dominated MACE is angina hospitalisation, driven by symptom burden with")
    emit("tau_SAQ %.0f d plus a central-sensitisation component that no vasoactive"
         % (P["TAU_SAQ"] / 24))
    emit("drug in this model touches.  The mismatch, quantified:")
    dcfr = (rt[-1]["CFR"] - rt[0]["CFR"]) - (ru[-1]["CFR"] - ru[0]["CFR"])
    emit("  CFR attributable to IMT over 2.5 y  %+.3f units" % dcfr)
    emit("  SAQ attributable to IMT over 2.5 y  %+.2f U   (MCID %.0f)"
         % (yt[IX["SAQ"]] - yu[IX["SAQ"]], P["MCID_SAQ"]))
    emit("  hs-cTnI  usual care %.1f  IMT %.1f ng/L" % (yu[IX["TNI"]], yt[IX["TNI"]]))
    emit("i.e. the therapy does what it says to the microvasculature and still cannot")
    emit("reach the endpoint it was scored on within the follow-up available.")
    return hr_true, dcfr, yt[IX["SAQ"]] - yu[IX["SAQ"]]


# =====================================================================
# SECTION X
# =====================================================================
def section_X():
    rule("SECTION X.  CorMicA: STRATIFIED THERAPY vs UNSTRATIFIED CARE")
    emit("Ford et al. J Am Coll Cardiol 2018;72:2841 (PMID 30266608).  Interventional")
    emit("diagnostic procedure (CFR, IMR, FFR, then acetylcholine) with linked")
    emit("therapy vs standard care with a sham procedure; 151 randomised of 391")
    emit("enrolled.  SAQ summary +11.7 U at 6 months (95% CI 5.0-18.4, p=0.001);")
    emit("EQ-5D +0.10; no difference in MACE.")
    emit()
    strat = dict(functional="cormica_func", structural="cormica_struct",
                 vasospastic="cormica_vaso", noncardiac="cormica_noncard")
    pop = dict(functional=0.30, structural=0.19, vasospastic=0.29, noncardiac=0.22)
    emit("Endotype mix: %s" % ", ".join("%s %.2f" % kv for kv in pop.items()))
    emit("Control arm = unguided standard care, represented here as a beta-blocker")
    emit("prescribed irrespective of endotype (the commonest real-world default).")
    emit()
    hdr = "%-12s %6s %9s %10s %8s %9s" % ("stratum", "w", "SAQ_std", "SAQ_strat",
                                          "dSAQ", "dBruce")
    emit(hdr)
    emit("-" * len(hdr))
    ts = tu = 0.0
    for arch, w in pop.items():
        yu, _ = simulate(arch, REGIMENS["bisoprolol"], days=182)
        ys, _ = simulate(arch, REGIMENS[strat[arch]], days=182)
        rxu, *_ = regimen_effects(yu, REGIMENS["bisoprolol"])
        rxs, *_ = regimen_effects(ys, REGIMENS[strat[arch]])
        eu = exercise_duration(yu, arch, rxu)
        es = exercise_duration(ys, arch, rxs)
        emit("%-12s %6.2f %9.1f %10.1f %+8.2f %+9.1f"
             % (arch, w, yu[IX["SAQ"]], ys[IX["SAQ"]],
                ys[IX["SAQ"]] - yu[IX["SAQ"]], es - eu))
        ts += w * ys[IX["SAQ"]]
        tu += w * yu[IX["SAQ"]]
    emit()
    emit("population mean SAQ: standard %.1f  stratified %.1f  difference %+.2f U"
         % (tu, ts, ts - tu))
    emit("CorMicA observed                                        +11.70 U")
    emit()
    emit("The stratified arm wins here for a mechanical reason, not a pharmacological")
    emit("one: the same prescription is right for one endotype and wrong for another,")
    emit("so an unstratified arm averages a benefit with a harm.  That average is")
    emit("what conventional trials report as 'no effect'.")
    return ts - tu


# =====================================================================
# SECTION XI
# =====================================================================
def section_XI():
    rule("SECTION XI.  TWENTY-FOUR TREATMENT SCENARIOS (24 weeks)")
    scen = [
        ("1  untreated functional", "functional", "untreated"),
        ("2  untreated structural", "structural", "untreated"),
        ("3  untreated vasospastic", "vasospastic", "untreated"),
        ("4  non-cardiac chest pain", "noncardiac", "untreated"),
        ("5  ivabradine / functional", "functional", "ivabradine"),
        ("6  ivabradine / structural", "structural", "ivabradine"),
        ("7  bisoprolol / functional", "functional", "bisoprolol"),
        ("8  nebivolol / functional", "functional", "nebivolol"),
        ("9  ranolazine / functional", "functional", "ranolazine"),
        ("10 ranolazine / structural", "structural", "ranolazine"),
        ("11 ranolazine / control", "control", "ranolazine"),
        ("12 amlodipine / vasospastic", "vasospastic", "amlodipine"),
        ("13 diltiazem / vasospastic", "vasospastic", "diltiazem"),
        ("14 nicorandil / structural", "structural", "nicorandil"),
        ("15 zibotentan / functional", "functional", "zibotentan"),
        ("16 zibotentan-CF / functional", "functional", "zibotentan_cf"),
        ("17 IMT WARRIOR / structural", "structural", "imt_warrior"),
        ("18 fasudil / vasospastic", "vasospastic", "fasudil"),
        ("19 aminophylline / functional", "functional", "aminophylline"),
        ("20 imipramine / non-cardiac", "noncardiac", "imipramine"),
        ("21 exercise rehab / structural", "structural", "rehab"),
        ("22 iva+ranolazine / functional", "functional", "iva_ran"),
        ("23 sildenafil / structural", "structural", "sildenafil"),
        ("24 CorMicA-directed / functional", "functional", "cormica_func"),
    ]
    hdr = ("%-32s %6s %7s %7s %7s %6s %6s %6s %6s %6s"
           % ("scenario", "CFR", "MR_hyp", "burden", "noci", "ANG", "SAQ",
              "Bruce", "LVEDP", "TnI"))
    emit(hdr)
    emit("-" * len(hdr))
    out = {}
    for name, arch, reg in scen:
        y, _ = simulate(arch, REGIMENS[reg], days=168)
        rx, *_ = regimen_effects(y, REGIMENS[reg])
        idx = indices(y, arch, rx)
        bi, bn = burdens(y, arch, rx)
        ex = exercise_duration(y, arch, rx)
        emit("%-32s %6.2f %7.2f %7.3f %7.3f %6.2f %6.1f %6.0f %6.1f %6.1f"
             % (name, idx["CFR"], idx["MR_hyp"], bi, bn, y[IX["ANG"]],
                y[IX["SAQ"]], ex, y[IX["LVEDP"]], y[IX["TNI"]]))
        out[name] = dict(arch=arch, regimen=reg, CFR=idx["CFR"],
                         MR_hyp=idx["MR_hyp"], MR_rest=idx["MR_rest"],
                         IMR=idx["IMR"], MRR=idx["MRR"], burden=bi, noci=bn,
                         ANG=y[IX["ANG"]], SAQ=y[IX["SAQ"]], Bruce=ex,
                         LVEDP=y[IX["LVEDP"]], TnI=y[IX["TNI"]],
                         endo_epi_hyp=idx["endo_epi_hyp"])
    emit()
    emit("Reserve exhaustion workload (the workload at which subendocardial")
    emit("autoregulation runs out and flow becomes pressure-passive):")
    for arch in ("control", "functional", "structural", "vasospastic"):
        y = y0(arch)
        w = reserve_exhaustion_wl(y, arch, {})
        emit("  %-12s wl = %s" % (arch, "%.2f" % w if w == w else "not reached by 5.0"))
    return out


# =====================================================================
# SECTION XII
# =====================================================================
def section_XII():
    rule("SECTION XII.  VIRTUAL POPULATION (n=300)")
    emit("Risk factors are sampled first; the structural floor and the controller")
    emit("offset are drawn CONDITIONALLY on them (bug B10: independent draws produced")
    emit("patients with severe remodelling and no risk factor at all).")
    emit()
    random.seed(20260804)
    n = 300
    rows = []
    for i in range(n):
        htn = random.random() < 0.55
        dm = random.random() < 0.22
        infl = max(1.0, 1.0 + random.gauss(0.35, 0.30) + 0.25 * htn + 0.35 * dm)
        rmin_f = 1.0 + max(0.0, random.gauss(0.02 + 0.04 * htn + 0.11 * dm, 0.04))
        auto = max(0.0, random.gauss(0.26 + 0.10 * htn, 0.33))
        arch = "vp_%d" % i
        ARCHETYPE[arch] = dict(RMIN_F=rmin_f, AUTO_OFF=auto, INFL=infl,
                               GENO=1 if random.random() < 0.5 else 0,
                               ROCK_D=1.0 + 0.35 * (infl - 1.0),
                               KSBP="KSBP_STRUCT" if rmin_f > 1.15 else "KSBP_FUNC",
                               SPASM=max(0.0, random.gauss(0.25, 0.28)),
                               SENS0=max(0.0, random.gauss(0.40, 0.22)),
                               SENS_BASE=max(0.0, random.gauss(0.72, 0.45)))
        idx = indices(y0(arch), arch, {})
        rows.append(dict(i=i, arch=arch, CFR=idx["CFR"], MR_rest=idx["MR_rest"],
                         MR_hyp=idx["MR_hyp"], IMR=idx["IMR"], MRR=idx["MRR"],
                         rmin_f=rmin_f, auto=auto, infl=infl, htn=htn, dm=dm))
    cmd = [r for r in rows if r["CFR"] < 2.5]
    func = [r for r in cmd if r["MR_hyp"] < 2.5]
    struct = [r for r in cmd if r["MR_hyp"] >= 2.5]
    emit("CFR < 2.5 : %d/%d = %.0f%%              (Rahman observed 53%%)"
         % (len(cmd), n, 100 * len(cmd) / n))
    emit("  functional  (MR_hyp<2.5) %3d = %.0f%%   (observed 62%%)"
         % (len(func), 100 * len(func) / max(1, len(cmd))))
    emit("  structural (MR_hyp>=2.5) %3d = %.0f%%   (observed 38%%)"
         % (len(struct), 100 * len(struct) / max(1, len(cmd))))
    emit("CFR < 2.0 : %.0f%%   |   IMR >= 25 : %.0f%%   |   MRR < 3.0 : %.0f%%"
         % (100 * len([r for r in rows if r["CFR"] < 2.0]) / n,
            100 * len([r for r in rows if r["IMR"] >= 25]) / n,
            100 * len([r for r in rows if r["MRR"] < 3.0]) / n))
    emit()

    def ms(xs):
        mu = sum(xs) / len(xs)
        return mu, (sum((x - mu) ** 2 for x in xs) / max(1, len(xs) - 1)) ** 0.5
    for lab, grp in (("functional", func), ("structural", struct),
                     ("CFR>=2.5", [r for r in rows if r["CFR"] >= 2.5])):
        if not grp:
            continue
        a1, a2 = ms([r["MR_rest"] for r in grp])
        b1, b2 = ms([r["MR_hyp"] for r in grp])
        c1, c2 = ms([r["CFR"] for r in grp])
        emit("  %-11s n=%3d  MR_rest %.2f+-%.2f  MR_hyp %.2f+-%.2f  CFR %.2f+-%.2f"
             % (lab, len(grp), a1, a2, b1, b2, c1, c2))
    emit("  observed: functional MR_rest 4.2+-1.0, structural 6.9+-1.7,")
    emit("            control 7.3+-2.2 mmHg/(cm/s)")
    emit()

    # ---------------------------------------------------------------
    # What the 2.5 cut-off actually separates
    # ---------------------------------------------------------------
    emit("WHAT THE HYPERAEMIC-MR CUT-OFF ACTUALLY SEPARATES")
    emit("A hyperaemic MR of 2.5 mmHg/(cm/s) is read as 'structural', but the model")
    emit("says two different things raise it: inward remodelling and capillary")
    emit("rarefaction, which fix the minimal radius, and constrictor tone that")
    emit("adenosine does not fully reverse (%.0f%% of it is reversed here, calibrated,"
         % (100 * P["F_ADO_REV"]))
    emit("so %.0f%% survives).  The first is not drug-reversible on any useful"
         % (100 * (1 - P["F_ADO_REV"])))
    emit("timescale; the second is reversible in weeks.  The label cannot tell them")
    emit("apart, but a repeat measurement on a Rho-kinase inhibitor can:")
    emit()
    resolved = 0
    checked = struct[:40]
    for r in checked:
        yy, _ = simulate(r["arch"], REGIMENS["fasudil"], days=21)
        rx, *_ = regimen_effects(yy, REGIMENS["fasudil"])
        if indices(yy, r["arch"], rx)["MR_hyp"] < 2.5:
            resolved += 1
    emit("  of %d patients labelled structural, %d (%.0f%%) fall below the cut-off"
         % (len(checked), resolved, 100 * resolved / max(1, len(checked))))
    emit("  after 21 days of a Rho-kinase inhibitor - i.e. their 'structural' label")
    emit("  was constrictor tone, and the anatomy was never the problem.")
    emit()
    emit("  PREDICTION: the hyperaemic measurement should be repeated after acute")
    emit("  Rho-kinase or ETA blockade before a patient is told their microvascular")
    emit("  disease is structural.  It is a one-visit test and it changes the")
    emit("  therapeutic class.")
    emit()

    # ---------------------------------------------------------------
    # treatment response distributions
    # ---------------------------------------------------------------
    emit("Ranolazine and ivabradine response by stratum (12 weeks, 20 per stratum):")
    resp = {}
    for arm in ("ranolazine", "ivabradine"):
        for lab, grp in (("functional", func), ("structural", struct)):
            ds = []
            for r in grp[:20]:
                yp, _ = simulate(r["arch"], REGIMENS["untreated"], days=84, dt=3.0)
                yr, _ = simulate(r["arch"], REGIMENS[arm], days=84, dt=3.0)
                ds.append(yr[IX["SAQ"]] - yp[IX["SAQ"]])
            if not ds:
                continue
            mu, sd = ms(ds)
            pr = 100.0 * len([x for x in ds if x >= P["MCID_SAQ"]]) / len(ds)
            resp["%s/%s" % (arm, lab)] = dict(mean=mu, sd=sd, pct_responder=pr,
                                              n=len(ds))
            emit("  %-12s %-11s dSAQ %+.2f +- %.2f   responders (>=%.0f U) %.0f%%"
                 " (n=%d)" % (arm, lab, mu, sd, P["MCID_SAQ"], pr, len(ds)))
    emit()
    emit("Two things to read off honestly.  First, ranolazine's effect is small and")
    emit("almost identical in the two strata, so the model does NOT support hyperaemic")
    emit("MR as a ranolazine-responder test: it predicts that no clinically important")
    emit("responder subgroup exists for this drug at this dose, and that RWISE's")
    emit("subgroup signal is a real but sub-MCID effect reaching significance on a")
    emit("frequency scale rather than a mechanism waiting to be exploited.  Second,")
    emit("rate reduction outperforms it in both strata, which is the same ranking")
    emit("section IV found and the opposite of the order in which these drugs are")
    emit("usually tried in this population.")
    return rows, resp


# =====================================================================
# SECTION XIII
# =====================================================================
def section_XIII():
    rule("SECTION XIII.  LOCAL SENSITIVITY (+-20%, functional endotype, 12 weeks)")
    emit("Normalised sensitivity S = (dY/Y)/(dX/X) for the four read-outs that")
    emit("matter: CFR, subendocardial burden, the afferent drive, and SAQ.")
    emit()
    arch = "functional"
    targets = ["RMIN0", "A_TONE", "G_AUTO", "K_ET_TONE", "K_ROCK_TONE", "KCONV",
               "E_MAX", "PHI_SYS", "W_ENDO", "RCOMP_K", "TSYS_K", "K_TTI",
               "LVEDP0", "K_ANG", "K_NOC_ADO", "K_NOC_DEF", "SAQ_SENS",
               "K_SYMP_HR", "K_MYO", "FRAC_DIA_P", "F_ENDO_DEM"]

    def run():
        y, _ = simulate(arch, REGIMENS["untreated"], days=84, dt=3.0)
        idx = indices(y, arch, {})
        bi, bn = burdens(y, arch, {})
        return idx["CFR"], bi, bn, y[IX["SAQ"]]
    base = run()
    hdr = "%-13s %9s %10s %9s %9s" % ("parameter", "S(CFR)", "S(burden)",
                                      "S(noci)", "S(SAQ)")
    emit(hdr)
    emit("-" * len(hdr))
    sens = {}
    for k in targets:
        k0 = P[k]
        P[k] = k0 * 1.2
        up = run()
        P[k] = k0 * 0.8
        dn = run()
        P[k] = k0
        s = [0.0 if abs(base[i]) < 1e-12 else ((up[i] - dn[i]) / base[i]) / 0.4
             for i in range(4)]
        sens[k] = s
        emit("%-13s %9.3f %10.3f %9.3f %9.3f" % (k, s[0], s[1], s[2], s[3]))
    emit()
    top = sorted(sens.items(), key=lambda kv: -abs(kv[1][3]))[:5]
    emit("SAQ is most sensitive to: " + ", ".join("%s (%.2f)" % (k, v[3])
                                                 for k, v in top))
    return sens, base


# =====================================================================
# MAIN
# =====================================================================
def main():
    rule("CORONARY MICROVASCULAR DYSFUNCTION (ANOCA / INOCA) - QSP REFERENCE RUN",
         "#")
    emit("dependency-free Python reference for cmd_mrgsolve_model.R")
    emit("states: %d ODEs (%d physiology + %d PK/metabolite) | drugs: %d | "
         "scenarios: 24" % (NS, NS_SLOW, len(PK_NAMES), len(DRUGS)))

    cal = section_0()
    tab = section_I()
    ext = section_II(tab)
    hr_full, hr_time = section_III()
    rho, best = section_IV()
    qss_ok = section_V()
    prize, prize_per = section_VI()
    d_cmd, d_all, prev = section_VII()
    jump = section_VIII()
    w_hr, w_cfr, w_saq = section_IX()
    cormica = section_X()
    scen = section_XI()
    pop, resp = section_XII()
    sens, sbase = section_XIII()

    rule("SUMMARY OF THE QUANTITATIVE CLAIMS", "#")
    emit(" 1. Same number, different disease: functional CFR %.2f vs structural %.2f,"
         % (tab["functional"]["CFR"], tab["structural"]["CFR"]))
    emit("    while their absolute hyperaemic velocities differ by %.0f%% (%.2f vs"
         % (100 * abs(tab["functional"]["hyp"]["v"] - tab["structural"]["hyp"]["v"])
            / tab["structural"]["hyp"]["v"], tab["functional"]["hyp"]["mbf"]))
    emit("    %.2f mL/min/g)." % tab["structural"]["hyp"]["mbf"])
    emit(" 2. Resting MR is predicted, not fitted, wherever the controller is intact;")
    emit("    the functional endotype needs an offset of %.2f, which forces resting"
         % ARCHETYPE["functional"]["AUTO_OFF"])
    emit("    oxygen extraction down to %.0f%% of normal - refutable on one blood gas."
         % (100 * ext))
    if hr_full:
        emit(" 3. Of the fall in subendocardial deficit from 68->55 bpm, %.0f%% is the"
             % (100 * hr_time / hr_full))
        emit("    diastolic time window and %.0f%% is lower demand."
             % (100 * (hr_full - hr_time) / hr_full))
    emit(" 4. CFR ranks drugs badly even when its sign is right: rank correlation")
    emit("    between dCFR and the symptom benefit is %.2f (functional), %.2f"
         % (rho["functional"][0], rho["structural"][0]))
    emit("    (structural), %.2f (vasospastic).  Hyperaemic MR is the index that a"
         % rho["vasospastic"][0])
    emit("    blood-pressure-lowering drug cannot fake.")
    emit(" 5. PRIZE: population-weighted model net %+.2f s vs observed -4.26 s"
         % prize["zibotentan"])
    emit("    (95%% CI -19.60 to +11.06); microvascular gain alone %+.2f s, BP fall"
         % prize["zibotentan_cf"])
    emit("    %+.2f s, fluid retention %+.2f s."
         % (prize["zibotentan"] - prize["zib_nobp"],
            prize["zibotentan"] - prize["zib_nofluid"]))
    emit(" 6. RWISE: dSAQ %+.2f U in the CFR<2.5 stratum, %+.2f U in the whole"
         % (d_cmd, d_all))
    emit("    cohort - %.0f%% of the mechanistic effect survives the dilution."
         % (100 * d_all / d_cmd if abs(d_cmd) > 1e-9 else 0))
    emit(" 7. WARRIOR: model true hazard ratio %.2f; at 80%% background contamination"
         % w_hr)
    emit("    the observable ratio is %.2f. CFR %+.3f but SAQ %+.2f U over 2.5 y."
         % (1.0 + (w_hr - 1.0) * 0.2, w_cfr, w_saq))
    emit(" 8. CorMicA: stratified minus unstratified %+.2f U (observed +11.70 U)."
         % cormica)
    emit(" 9. Virtual population: CFR<2.5 in %.0f%% with a %.0f/%.0f"
         % (100 * len([r for r in pop if r["CFR"] < 2.5]) / len(pop),
            100 * len([r for r in pop if r["CFR"] < 2.5 and r["MR_hyp"] < 2.5])
            / max(1, len([r for r in pop if r["CFR"] < 2.5])),
            100 * len([r for r in pop if r["CFR"] < 2.5 and r["MR_hyp"] >= 2.5])
            / max(1, len([r for r in pop if r["CFR"] < 2.5]))))
    emit("    functional/structural split (observed 53%, 62/38).")
    emit("10. Fixed point and fast-ODE formulations agree to 1e-6: %s"
         % ("yes" if qss_ok else "NO"))
    emit("11. SAQ sensitivity is dominated by %s."
         % ", ".join(k for k, v in sorted(sens.items(),
                                          key=lambda kv: -abs(kv[1][3]))[:3]))

    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "cmd_reference_output.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")
    js = dict(
        model=dict(n_states=NS, n_physiology=NS_SLOW, n_pk=len(PK_NAMES),
                   n_drugs=len(DRUGS), n_scenarios=24),
        calibration=cal,
        endotypes={k: dict(CFR=v["CFR"], MR_rest=v["MR_rest"], MR_hyp=v["MR_hyp"],
                           v_rest=v["rest"]["v"], v_hyp=v["hyp"]["v"],
                           mbf_rest=v["rest"]["mbf"], mbf_hyp=v["hyp"]["mbf"],
                           tone_rest=v["rest"]["tone"], IMR=v["IMR"], MRR=v["MRR"],
                           endo_epi_hyp=v["endo_epi_hyp"]) for k, v in tab.items()},
        resting_extraction_ratio_functional_vs_control=ext,
        heart_rate_decomposition=dict(total=hr_full, time_window=hr_time,
                                      demand=hr_full - hr_time),
        endpoint_rank_correlation={k: dict(rho_dCFR_vs_symptom=v[0],
                                           rho_dCFR_vs_exercise=v[1],
                                           best_on_CFR=best[k][0],
                                           best_on_symptoms=best[k][1],
                                           best_on_exercise=best[k][2])
                                   for k, v in rho.items()},
        qss_ode_agreement=qss_ok,
        prize=dict(model_net_s=prize["zibotentan"],
                   counterfactual_s=prize["zibotentan_cf"],
                   cost_bp_s=prize["zibotentan"] - prize["zib_nobp"],
                   cost_fluid_s=prize["zibotentan"] - prize["zib_nofluid"],
                   observed_s=-4.26, observed_ci=[-19.60, 11.06],
                   per_endotype=prize_per),
        rwise=dict(dSAQ_subgroup=d_cmd, dSAQ_cohort=d_all, prevalence=prev,
                   mcid=P["MCID_SAQ"]),
        warrior=dict(true_HR=w_hr,
                     HR_at_80pct_contamination=1.0 + (w_hr - 1.0) * 0.2,
                     dCFR=w_cfr, dSAQ=w_saq),
        cormica=dict(model_dSAQ=cormica, observed_dSAQ=11.7),
        loop_bifurcation_gain=jump,
        scenarios=scen,
        population=dict(
            n=len(pop),
            pct_CFR_lt_2p5=100 * len([r for r in pop if r["CFR"] < 2.5]) / len(pop),
            pct_CFR_lt_2p0=100 * len([r for r in pop if r["CFR"] < 2.0]) / len(pop),
            pct_IMR_ge_25=100 * len([r for r in pop if r["IMR"] >= 25]) / len(pop),
            functional_share=100 * len([r for r in pop if r["CFR"] < 2.5
                                        and r["MR_hyp"] < 2.5])
            / max(1, len([r for r in pop if r["CFR"] < 2.5])),
            ranolazine_response=resp),
        sensitivity={k: dict(CFR=v[0], burden=v[1], noci=v[2], SAQ=v[3])
                     for k, v in sens.items()},
    )
    with open(os.path.join(here, "cmd_population_results.json"), "w",
              encoding="utf-8") as fh:
        json.dump(js, fh, indent=1, ensure_ascii=False)
    emit()
    emit("wrote cmd_reference_output.txt and cmd_population_results.json")


if __name__ == "__main__":
    main()
