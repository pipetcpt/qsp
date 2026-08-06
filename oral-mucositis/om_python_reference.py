#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
om_python_reference.py -- executable reference implementation of the oral
mucositis (OM) QSP model.

THE ONE IDEA
------------
Oral mucositis is not a toxicity that a drug "has".  It is a BALANCE-SHEET
FAILURE OF A RENEWING EPITHELIUM.  The buccal mucosa is a conveyor belt: a
clonogenic basal pool S feeds three transit-amplifying generations P1..P3,
which feed a post-mitotic barrier D, which is shed.  An ulcer exists exactly
while D has fallen below the threshold at which the barrier is continuous.

Every cytotoxic insult -- alkylator, antimetabolite, photon -- acts on the
PROLIFERATIVE end of that belt (S, P) and spares the post-mitotic end (D).
So the insult is invisible for one epithelial transit time, and then the
barrier collapses because nothing arrived to replace what was shed.

That single structure forces the two clocks that organise the whole disease:

    ONSET    is set by  (transit time of the belt)   ~ 9-14 d, and by the
             SIZE of the kill, i.e. by the insult term;
    DURATION is set by  (regrowth rate of S)         ~ days, i.e. by the
             REGENERATION term.

The insult term is over in hours-to-days (melphalan t1/2 ~ 1.2 h).  The
regeneration term runs for weeks.  Therefore:

    * every intervention that scales the INSULT (cryotherapy, dose
      reduction, fractionation) moves INCIDENCE and ONSET;
    * every intervention that scales REGENERATION (palifermin/KGF,
      photobiomodulation, glutamine) moves DURATION;

and an agent aimed at the wrong term cannot buy the other endpoint no matter
how large its effect size.  om_analysis.py computes this and its
consequences, including one that refuted this file's own first framing:
because KGF ENLARGES the cycling pool it hands a CYCLE-ACTIVE cytotoxic a
bigger target -- but that penalty is structurally ABSENT for alkylators and
photons, whose kill is first-order in the pool, so the fraction killed does
not change.  Section 5 reports the flat sweep and says why it is flat.

STATE VECTOR (50 ODEs)
----------------------
  epithelium (8)     S Sq Sd P1 P2 P3 D Dh
  inflammation (8)   ROS NFkB TNF IL1b IL6 CER MB NS
  regeneration (3)   KGFe PBMe GLNe
  PK melphalan (2)   A_mel_c A_mel_p
  PK 5-FU (2)        A_5fu_c A_5fu_gut
  PK methotrexate(3) A_mtx_c A_mtx_p1 A_mtx_p2
  PK cisplatin (2)   A_cis_c A_cis_p
  PK palifermin (2)  A_pal_c A_pal_p
  PK morphine (3)    A_mor_c A_mor_p A_mor_e
  mucosal tissue (4) Cm_mel Cm_5fu Cm_mtx Cm_cis
  myelosuppression(5) Prol Tr1 Tr2 Tr3 Circ
  accumulators (8)   cAUCmuc cUlcD cOpiD cBEDt cBEDm cInf cSev cPain

All rates are per DAY.  Concentrations are ug/mL (melphalan, 5-FU, palifermin
as ng/mL, morphine as ng/mL) -- each PK block documents its own units.
"""
import numpy as np
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------------
# state index
# ----------------------------------------------------------------------------
NAMES = [
    "S", "Sq", "Sd", "P1", "P2", "P3", "D", "Dh",
    "ROS", "NFkB", "TNF", "IL1b", "IL6", "CER", "MB", "NS",
    "KGFe", "PBMe", "GLNe",
    "A_mel_c", "A_mel_p",
    "A_5fu_c", "A_5fu_gut",
    "A_mtx_c", "A_mtx_p1", "A_mtx_p2",
    "A_cis_c", "A_cis_p",
    "A_pal_c", "A_pal_p",
    "A_mor_c", "A_mor_p", "A_mor_e",
    "Cm_mel", "Cm_5fu", "Cm_mtx", "Cm_cis",
    "Prol", "Tr1", "Tr2", "Tr3", "Circ",
    "cAUCmuc", "cUlcD", "cOpiD", "cBEDt", "cBEDm", "cInf", "cSev", "cPain",
]
IX = {n: i for i, n in enumerate(NAMES)}
NST = len(NAMES)

# ----------------------------------------------------------------------------
# parameters
# ----------------------------------------------------------------------------
P = dict(
    # --- epithelial conveyor belt -------------------------------------------
    # Human buccal epithelium turnover 9-14 d (Squier & Kremer 2001).  The belt
    # is S -> P1 -> P2 -> P3 -> D -> shed.  k_p sets the transit-amplifying
    # residence; k_shed sets barrier lifetime.  Their sum must land the total
    # transit at ~11 d or the ONSET LAG of the model is wrong, and the onset
    # lag is a calibration target, not a free choice.
    lamS=0.14,        # /d  net basal clonogen renewal at homeostasis
    S0=1.0,           #     normalised homeostatic clonogen density
    aS=0.55,          # /d  flux S -> P1 (per unit S)
    k_p=0.62,         # /d  P1->P2->P3->D maturation (3 stages -> mean 4.8 d)
    k_shed=0.30,      # /d  desquamation of the barrier layer
    fmax=6.0,         #     ceiling on the regeneration multiplier
    greg=7.0,         #     gain of the barrier-deficit -> proliferation loop
    Kreg=0.35,        #     half-saturation of that loop (deficit fraction)
    k_ab=0.55,        # /d  abortive-division death of lethally hit basal cells
    ndiv=0.35,        #     residual flux a doomed cell still contributes
    # --- the QUIESCENT CLONOGEN RESERVE ------------------------------------
    # WITHOUT THIS THE MODEL CANNOT HEAL, AND THE FAILURE IS SILENT.
    # dS/dt = lamS*freg*S*(1-S/S0) - k*S is proportional to S, so S = 0 is an
    # ABSORBING state: any insult large enough to drive the cycling pool to
    # numerical zero sterilises the tissue permanently.  The first calibration
    # run hit exactly that -- to reach a 5.5 d onset the fit chose a potency
    # that annihilated S, and severe mucositis then lasted 39 days and never
    # resolved, because there was nothing left to repopulate from.
    #
    # The biology that is missing is the label-retaining, out-of-cycle stem
    # cell.  It is spared by cycle-active agents, damaged much less by the
    # rest, and it is what actually repopulates an ablated epithelium.  Adding
    # it removes the absorbing state AND puts a second, slower time constant
    # into the REGENERATION arm -- which is a strengthening of the two-clock
    # claim, not a rescue of it: kact is now part of clock 2 and can be
    # perturbed independently of lamS.
    Sq0=0.16,         #     size of the quiescent reserve, relative to S0
    q_res=0.18,       #     quiescent kill rate, as a fraction of the cycling
    kact=1.35,        # /d  activation of reserve -> cycling pool on deficit
    krest=0.10,       # /d  restoration of the reserve once S has recovered
    # --- ulceration link ----------------------------------------------------
    Dcrit=0.34,       #     barrier fraction at which the epithelium ulcerates
    wD=0.055,         #     steepness of that transition
    # --- inflammation -------------------------------------------------------
    kROSd=6.0,        # /d
    kNFb=0.5,         # /d  basal NF-kB turnover
    gROS=1.10,        #     ROS  -> NF-kB
    gTNF=0.42,        #     TNF  -> NF-kB  (POSITIVE FEEDBACK: keep gain<1)
    gMB=0.75,         #     PAMP -> NF-kB
    KTNF=1.0, KMB=1.0, KROS=1.0,
    kTNF=2.4, kIL1=2.2, kIL6=1.8,   # /d cytokine turnover
    pTNF=2.4, pIL1=2.2, pIL6=1.8,   # /d production per unit NF-kB
    kCER=1.6, pCER=1.9,             # /d ceramide (aSMase) turnover/production
    aCER=0.22,        #     ceramide -> direct apoptosis of D (post-mitotic)
    aTNFs=0.30,       #     TNF -> extra basal apoptosis
    # --- microbial colonisation of the ulcer bed ----------------------------
    kMBg=1.5, kMBd=0.9, MBmax=1.0, KANC=0.5,
    # --- nociception --------------------------------------------------------
    kNS=0.45, gNS=1.5, KNSc=0.8,
    # --- regeneration modifiers ---------------------------------------------
    kKGFe=0.35,       # /d  turnover of the KGF pharmacodynamic effect pool
    eKGF=2.6,         #     max fold-increase of aS/lamS at saturating KGF
    KpalPD=8.0,       # ng/mL  EC50 of palifermin on KGFR
    fcyc_KGF=0.85,    #     fraction of the KGF proliferative gain that also
                      #     enlarges the CYTOTOXIC TARGET (the paradox knob)
    kPBMe=0.5, ePBM=0.55, KPBM=4.0,   # photobiomodulation (J/cm2 per session)
    kGLNe=0.7, eGLN=0.30,             # glutamine
    # --- benzydamine / topical anti-inflammatory ----------------------------
    eBZD=0.45,        #     max fractional suppression of NF-kB production
    # --- melphalan PK (ug/mL; V in mL, CL in mL/d) --------------------------
    # Nath 2010, Mougenot 2004: t1/2 ~1.2 h, V ~0.5 L/kg, AUC ~12 ug.h/mL at
    # 200 mg/m2.  The model reproduces Cmax ~9.6 ug/mL and AUC ~16 ug.h/mL.
    V_mel=37500.0, CL_mel=520000.0, Q_mel=150000.0, Vp_mel=30000.0,
    # --- 5-FU PK (ug/mL) with saturable DPD ---------------------------------
    # bolus t1/2 ~12 min; 1000 mg/m2/d CI gives Css ~0.6 ug/mL
    V_5fu=25000.0, Vmax_5fu=2.2e6, Km_5fu=2.0, CL_5fu_lin=2.0e6,
    # --- methotrexate PK (uM; V in L, CL in L/d) ----------------------------
    V_mtx=20.0, CL_mtx=140.0, Q1_mtx=90.0, Vp1_mtx=12.0,
    Q2_mtx=8.0, Vp2_mtx=30.0,
    # --- cisplatin, free (ultrafilterable) platinum (ug/mL) -----------------
    V_cis=45000.0, CL_cis=1.83e6, Q_cis=2.0e5, Vp_cis=60000.0,
    # --- palifermin PK (ng/mL) ---------------------------------------------
    # label / Spielberger 2004: t1/2 ~4.5 h, Vss ~5 L/kg
    V_pal=120000.0, CL_pal=1.386e6, Q_pal=4.0e5, Vp_pal=255000.0,
    # --- morphine PK (ng/mL) -----------------------------------------------
    V_mor=1.0e5, CL_mor=2.16e6, Q_mor=3.0e5, Vp_mor=1.5e5, ke0_mor=6.0,
    EC50_mor=25.0, Emax_mor=6.5, k_titr=1.35e7,
    # --- mucosal (perfusion-limited) tissue compartment ---------------------
    # THIS IS THE STRUCTURE THAT MAKES CRYOTHERAPY A FLOW PROBLEM.
    # dCm/dt = (Q/Vm)*(Cart - Cm/Kp).  Q is scaled by the cryotherapy factor.
    # Oral mucosal blood flow is ~0.3-0.6 mL/min/g, so with a tissue:plasma
    # partition near unity the equilibration rate constant is ~500/d (t_eq a
    # couple of minutes).  THIS NUMBER IS NOT FREE -- and it is what makes the
    # naive "cryotherapy works by cutting delivery" story fail: a tissue that
    # re-equilibrates in minutes tracks plasma no matter what the flow is.
    keq_muc=500.0,    # /d
    Kp_mel=0.75, Kp_5fu=0.90, Kp_mtx=0.55, Kp_cis=0.45,
    fQ_cryo=0.22,     #     mucosal blood flow during ice, fraction of basal
    # The second, temperature arm of cryotherapy: cooling the mucosa to ~24 C
    # slows alkylation chemistry and cell-cycle progression with a Q10 of
    # ~2.5.  Section 5 decomposes the observed cryotherapy effect into this
    # arm and the flow arm and shows the flow arm is nearly empty.
    Q10_cryo=2.5, dT_cryo=13.0,
    # --- potency of each insult on the clonogen pool ------------------------
    # per (ug/mL * d) of MUCOSAL exposure
    pot_mel=0.115, pot_5fu=0.030, pot_mtx=0.085, pot_cis=0.020,
    wP1=1.0, wP2=0.6, wP3=0.3,   # cytotoxic weight down the amplifying belt
    acc_kp=0.35,      #     regenerative shortening of the transit time
    sens=1.0,         #     individual sensitivity multiplier (IIV handle)
    fS_cyc=1.0,       #     cycle-specificity weight for S
    fP_cyc=1.0,
    # 5-FU and MTX are S-phase specific: their kill scales with the cycling
    # fraction, which the regeneration loop RAISES.  Alkylators/photons do not.
    cycspec_5fu=1.0, cycspec_mtx=1.0, cycspec_mel=0.0, cycspec_cis=0.0,
    # --- radiation ----------------------------------------------------------
    alpha_m=0.30, beta_m=0.030,     # mucosa alpha/beta = 10 Gy
    alpha_t=0.30, beta_t=0.030,     # HNSCC tumour alpha/beta = 10 Gy
    rad_pot=0.62,     #     conversion of mucosal log-cell-kill to clonogen loss
    Tk_tum=28.0, lam_tum=0.693 / 3.5,   # tumour accelerated repopulation
    # --- myelosuppression (Friberg) ----------------------------------------
    MTT=125.0 / 24.0, gam_anc=0.161, Circ0=5.0, ktr_scale=4.0,
    slope_mel=0.055, slope_5fu=0.0016, slope_cis=0.030, slope_mtx=0.020,
    # --- clinical linkage ---------------------------------------------------
    pain_max=10.0, pain_n=1.0, pain_k=8.5,
    inf_base=0.004, inf_gain=0.55,
    opi_target=4.0,
)

# radiation and infusions are delivered through these time-varying inputs
ZERO_SCHED = dict(rt=[], inf=[], bolus=[], cryo=[], pbm=[], bzd=[], gln=[],
                  opioid=True)


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def _hill(x, K, n=1.0):
    x = max(x, 0.0)
    return x ** n / (K ** n + x ** n)


def _pos(x):
    return x if x > 0.0 else 0.0


def window(t, wins):
    """total instantaneous rate from a list of (t0, t1, rate) windows"""
    r = 0.0
    for t0, t1, v in wins:
        if t0 <= t < t1:
            r += v
    return r


def in_window(t, wins):
    for t0, t1 in wins:
        if t0 <= t < t1:
            return True
    return False


def rt_dose_rate(t, rt):
    """Gy/d delivered at time t.  rt = list of (t_start, dose_Gy, dur_d)."""
    r = 0.0
    for t0, d, dur in rt:
        if t0 <= t < t0 + dur:
            r += d / dur
    return r


# ----------------------------------------------------------------------------
# right-hand side
# ----------------------------------------------------------------------------
def rhs(t, y, p, sch):
    y = np.maximum(y, 0.0)
    (S, Sq, Sd, P1, P2, P3, D, Dh,
     ROS, NFkB, TNF, IL1b, IL6, CER, MB, NS,
     KGFe, PBMe, GLNe,
     Amelc, Amelp, A5fuc, A5fug,
     Amtxc, Amtxp1, Amtxp2, Acisc, Acisp,
     Apalc, Apalp, Amorc, Amorp, Amore,
     Cmmel, Cm5fu, Cmmtx, Cmcis,
     Prol, Tr1, Tr2, Tr3, Circ,
     _a, _b, _c, _d, _e, _f, _g, _h) = y

    dy = np.zeros(NST)

    # ---------------- pharmacokinetics --------------------------------------
    Cmel = Amelc / p["V_mel"]
    Cmelp = Amelp / p["Vp_mel"]
    dy[IX["A_mel_c"]] = (window(t, sch.get("mel", [])) - p["CL_mel"] * Cmel
                         - p["Q_mel"] * (Cmel - Cmelp))
    dy[IX["A_mel_p"]] = p["Q_mel"] * (Cmel - Cmelp)

    C5fu = A5fuc / p["V_5fu"]
    el5 = p["Vmax_5fu"] * C5fu / (p["Km_5fu"] + C5fu) + p["CL_5fu_lin"] * C5fu
    dy[IX["A_5fu_c"]] = window(t, sch.get("fu", [])) - el5
    dy[IX["A_5fu_gut"]] = 0.0

    Cmtx = Amtxc / p["V_mtx"]
    Cm1 = Amtxp1 / p["Vp1_mtx"]
    Cm2 = Amtxp2 / p["Vp2_mtx"]
    dy[IX["A_mtx_c"]] = (window(t, sch.get("mtx", [])) - p["CL_mtx"] * Cmtx
                         - p["Q1_mtx"] * (Cmtx - Cm1)
                         - p["Q2_mtx"] * (Cmtx - Cm2))
    dy[IX["A_mtx_p1"]] = p["Q1_mtx"] * (Cmtx - Cm1)
    dy[IX["A_mtx_p2"]] = p["Q2_mtx"] * (Cmtx - Cm2)

    Ccis = Acisc / p["V_cis"]
    Ccisp = Acisp / p["Vp_cis"]
    dy[IX["A_cis_c"]] = (window(t, sch.get("cis", [])) - p["CL_cis"] * Ccis
                         - p["Q_cis"] * (Ccis - Ccisp))
    dy[IX["A_cis_p"]] = p["Q_cis"] * (Ccis - Ccisp)

    Cpal = Apalc / p["V_pal"]
    Cpalp = Apalp / p["Vp_pal"]
    dy[IX["A_pal_c"]] = (window(t, sch.get("pal", [])) - p["CL_pal"] * Cpal
                         - p["Q_pal"] * (Cpal - Cpalp))
    dy[IX["A_pal_p"]] = p["Q_pal"] * (Cpal - Cpalp)

    # ---------------- mucosal perfusion-limited compartment -----------------
    # Oral cryotherapy is modelled ONLY here: it multiplies Q.  It touches no
    # potency, no repair, no cytokine.  Every cryotherapy result in the
    # analysis is therefore a consequence of flow and drug half-life alone.
    iced = in_window(t, sch.get("cryo", []))
    fQ = p["fQ_cryo"] if iced else 1.0
    Qd = p["keq_muc"] * fQ
    # temperature arm: Q10 on the alkylation / cell-cycle rate
    ftemp = (1.0 / p["Q10_cryo"] ** (p["dT_cryo"] / 10.0)) if iced else 1.0
    dy[IX["Cm_mel"]] = Qd * (Cmel - Cmmel / p["Kp_mel"])
    dy[IX["Cm_5fu"]] = Qd * (C5fu - Cm5fu / p["Kp_5fu"])
    dy[IX["Cm_mtx"]] = Qd * (Cmtx - Cmmtx / p["Kp_mtx"])
    dy[IX["Cm_cis"]] = Qd * (Ccis - Cmcis / p["Kp_cis"])

    # ---------------- regeneration signal -----------------------------------
    deficit = _pos(1.0 - D / 1.0)
    fdef = 1.0 + p["greg"] * deficit / (p["Kreg"] + deficit)

    kgf_drive = _hill(Cpal, p["KpalPD"])
    dy[IX["KGFe"]] = p["kKGFe"] * (kgf_drive - KGFe)
    fkgf = 1.0 + p["eKGF"] * KGFe

    pbm_in = window(t, sch.get("pbm", []))
    dy[IX["PBMe"]] = pbm_in - p["kPBMe"] * PBMe
    fpbm = 1.0 + p["ePBM"] * _hill(PBMe, p["KPBM"])

    gln_in = window(t, sch.get("gln", []))
    dy[IX["GLNe"]] = gln_in - p["kGLNe"] * GLNe
    fgln = 1.0 + p["eGLN"] * _hill(GLNe, 1.0)

    freg = min(fdef * fkgf * fpbm * fgln, p["fmax"])

    # cycling fraction seen by an S-phase-specific cytotoxic.  KGF raises it;
    # that is the whole palifermin-scheduling paradox, and fcyc_KGF is the
    # single parameter that carries it.
    fcyc = (1.0 + p["fcyc_KGF"] * (fkgf - 1.0)) * (fdef ** 0.5)

    # ---------------- cytotoxic insult on the proliferative pool ------------
    def kill_term(Cm, pot, cycspec):
        return ftemp * pot * Cm * (fcyc if cycspec > 0.5 else 1.0)

    k_chem = p["sens"] * (kill_term(Cmmel, p["pot_mel"], p["cycspec_mel"])
              + kill_term(Cm5fu, p["pot_5fu"], p["cycspec_5fu"])
              + kill_term(Cmmtx, p["pot_mtx"], p["cycspec_mtx"])
              + kill_term(Cmcis, p["pot_cis"], p["cycspec_cis"]))

    # radiation: linear-quadratic, applied as an instantaneous log-kill rate.
    # Within a fraction delivered over dur days, the quadratic term is scaled
    # by the instantaneous dose already given in that fraction.
    Rd = rt_dose_rate(t, sch.get("rt", []))
    dpf = sch.get("dose_per_fx", 2.0)
    k_rad = (p["sens"] * ftemp * p["rad_pot"] * Rd
             * (p["alpha_m"] + p["beta_m"] * dpf))

    k_tot = k_chem + k_rad + p["aTNFs"] * _hill(TNF, 1.5)

    # ---------------- epithelium --------------------------------------------
    # S  : viable CYCLING clonogens.  Logistic renewal, first-order lethal
    #      hit, plus activation from the quiescent reserve.
    # Sq : the quiescent reserve.  Hit at q_res x the cycling rate, drained
    #      into S in proportion to the cycling pool's deficit, and slowly
    #      restored once S is back.  Sq > 0 is what makes S = 0 non-absorbing.
    sdef = _pos(1.0 - S / p["S0"])
    act = p["kact"] * Sq * sdef
    dy[IX["S"]] = (p["lamS"] * freg * S * (1.0 - S / p["S0"]) + act
                   - k_tot * S)
    dy[IX["Sq"]] = (-p["q_res"] * k_tot * Sq - act
                    + p["krest"] * _pos(p["Sq0"] - Sq) * (S / p["S0"]))
    # Sd : lethally hit but still present; contributes abortive divisions
    dy[IX["Sd"]] = k_tot * S + p["q_res"] * k_tot * Sq - p["k_ab"] * Sd

    flux_in = p["aS"] * freg * (S + p["ndiv"] * Sd)
    # Regenerating epithelium does not merely divide faster, it also SHORTENS
    # the transit time (Dorr 1997).  Without this term the healing limb is
    # floored at 3/k_p + 1/k_shed ~ 8 d and no drug can shorten it -- which
    # would make the "duration is regeneration" claim untestable by
    # construction.
    kpe = p["k_p"] * (1.0 + p["acc_kp"] * (freg - 1.0))
    dy[IX["P1"]] = flux_in - kpe * P1 - p["wP1"] * k_tot * P1
    dy[IX["P2"]] = kpe * P1 - kpe * P2 - p["wP2"] * k_tot * P2
    dy[IX["P3"]] = kpe * P2 - kpe * P3 - p["wP3"] * k_tot * P3
    # D : post-mitotic barrier.  NO cytotoxic term -- this is the reason the
    # whole disease has a latent period.  Only ceramide/TNF apoptosis reaches
    # it, and that term is small.
    dy[IX["D"]] = (kpe * P3 - p["k_shed"] * D
                   - p["aCER"] * _hill(CER, 1.4) * D)
    dy[IX["Dh"]] = 0.0

    # ---------------- ulcer -------------------------------------------------
    A_ulc = 1.0 / (1.0 + np.exp((D - p["Dcrit"]) / p["wD"]))

    # ---------------- inflammation ------------------------------------------
    ros_in = 2.2 * Rd + 1.6 * (Cmmel * p["pot_mel"] + Cm5fu * p["pot_5fu"]
                               + Cmmtx * p["pot_mtx"] + Cmcis * p["pot_cis"])
    dy[IX["ROS"]] = ros_in - p["kROSd"] * ROS

    bzd = p["eBZD"] * _hill(window(t, sch.get("bzd", [])), 0.5)
    drive = (p["gROS"] * _hill(ROS, p["KROS"])
             + p["gTNF"] * _hill(TNF, p["KTNF"])
             + p["gMB"] * _hill(MB, p["KMB"]))
    dy[IX["NFkB"]] = p["kNFb"] * ((1.0 - bzd) * drive - NFkB)

    dy[IX["TNF"]] = p["pTNF"] * NFkB - p["kTNF"] * TNF
    dy[IX["IL1b"]] = p["pIL1"] * NFkB - p["kIL1"] * IL1b
    dy[IX["IL6"]] = p["pIL6"] * NFkB - p["kIL6"] * IL6
    dy[IX["CER"]] = p["pCER"] * (0.6 * NFkB + 0.4 * _hill(ROS, 0.8)) \
        - p["kCER"] * CER

    anc = max(Circ, 1e-3)
    dy[IX["MB"]] = (p["kMBg"] * A_ulc * (1.0 - MB / p["MBmax"])
                    - p["kMBd"] * MB * anc / (p["KANC"] + anc))

    # ---------------- nociception & analgesia -------------------------------
    Cmor_e = Amore
    opi = p["Emax_mor"] * _hill(Cmor_e, p["EC50_mor"])
    dy[IX["NS"]] = p["kNS"] * (p["gNS"] * A_ulc * (1.0 + 0.5 * _hill(IL1b, 1.0))
                               - NS)
    vas_raw = p["pain_max"] * _hill(p["pain_k"] * A_ulc * (0.4 + 0.6 * NS), 1.0)
    vas = max(vas_raw - opi, 0.0)

    # morphine PK, dosed by a simple titration rule when opioid=True
    Cmor = Amorc / p["V_mor"]
    Cmorp = Amorp / p["Vp_mor"]
    need = 0.0
    if sch.get("opioid", True) and vas_raw > p["opi_target"]:
        need = p["k_titr"] * (vas_raw - p["opi_target"])  # ng/d
    dy[IX["A_mor_c"]] = (need - p["CL_mor"] * Cmor
                         - p["Q_mor"] * (Cmor - Cmorp))
    dy[IX["A_mor_p"]] = p["Q_mor"] * (Cmor - Cmorp)
    dy[IX["A_mor_e"]] = p["ke0_mor"] * (Cmor - Amore)

    # ---------------- myelosuppression (Friberg) ----------------------------
    ktr = p["ktr_scale"] / p["MTT"]
    edrug = (p["slope_mel"] * Cmel + p["slope_5fu"] * C5fu
             + p["slope_cis"] * Ccis + p["slope_mtx"] * Cmtx)
    edrug = min(edrug, 0.98)
    fb = (p["Circ0"] / max(Circ, 1e-3)) ** p["gam_anc"]
    fb = min(fb, 6.0)
    dy[IX["Prol"]] = ktr * Prol * (1.0 - edrug) * fb - ktr * Prol
    dy[IX["Tr1"]] = ktr * (Prol - Tr1)
    dy[IX["Tr2"]] = ktr * (Tr1 - Tr2)
    dy[IX["Tr3"]] = ktr * (Tr2 - Tr3)
    dy[IX["Circ"]] = ktr * (Tr3 - Circ)

    # ---------------- accumulators ------------------------------------------
    dy[IX["cAUCmuc"]] = Cmmel + Cm5fu + Cmmtx + Cmcis
    dy[IX["cUlcD"]] = 1.0 if A_ulc > 0.5 else 0.0
    dy[IX["cOpiD"]] = 1.0 if need > 0 else 0.0
    # tumour BED with accelerated repopulation after Tk
    rt_on = 1.0 if (sch.get("rt") and sch["rt"][0][0] <= t
                    <= sch["rt"][-1][0] + 1.0) else 0.0
    rep = p["lam_tum"] * rt_on if t > p["Tk_tum"] else 0.0
    dy[IX["cBEDt"]] = Rd * (1.0 + dpf / (p["alpha_t"] / p["beta_t"])) \
        - rep / p["alpha_t"]
    dy[IX["cBEDm"]] = Rd * (1.0 + dpf / (p["alpha_m"] / p["beta_m"]))
    dy[IX["cInf"]] = p["inf_base"] + p["inf_gain"] * A_ulc * MB \
        / (p["KANC"] + anc)
    dy[IX["cSev"]] = 1.0 if A_ulc > 0.5 else 0.0
    dy[IX["cPain"]] = vas
    return dy


# ----------------------------------------------------------------------------
# steady state / initial condition
# ----------------------------------------------------------------------------
def y0(p=P):
    y = np.zeros(NST)
    y[IX["S"]] = p["S0"]
    y[IX["Sq"]] = p["Sq0"]
    # homeostatic belt: flux in = aS*S at freg=1 -> P1=P2=P3=aS*S/k_p,
    # D = k_p*P3/k_shed
    Pss = p["aS"] * p["S0"] / p["k_p"]
    y[IX["P1"]] = y[IX["P2"]] = y[IX["P3"]] = Pss
    y[IX["D"]] = p["k_p"] * Pss / p["k_shed"]
    y[IX["Circ"]] = p["Circ0"]
    y[IX["Prol"]] = y[IX["Tr1"]] = y[IX["Tr2"]] = y[IX["Tr3"]] = p["Circ0"]
    return y


def equilibrate(p=P, tmax=200.0):
    """run the drug-free system to its attractor; D is normalised to 1 there"""
    y = y0(p)
    sol = solve_ivp(rhs, (0, tmax), y, args=(p, dict(ZERO_SCHED)),
                    method="LSODA", rtol=1e-8, atol=1e-10,
                    dense_output=False)
    return sol.y[:, -1]


# the barrier scale is set ONCE, at the drug-free attractor, so that D=1 means
# "intact".  Dcrit is then a fraction of a real steady state, not of a guess.
def normalise(p=P):
    ss = equilibrate(p)
    return ss[IX["D"]]


# ----------------------------------------------------------------------------
# simulation driver
# ----------------------------------------------------------------------------
def breakpoints(sched, t_end):
    """
    Every discontinuity in the input functions, as an explicit integration
    boundary.

    THIS IS NOT COSMETIC.  A 2 Gy fraction is delivered over a 10-minute
    window (0.0069 d) and a palifermin bolus over 5 minutes.  An adaptive
    stepper handed a 105-day horizon steps straight over those windows and
    silently delivers a FRACTION of the intended dose -- the answer then
    depends on the integrator's step placement rather than on the regimen.
    An earlier build of this model did exactly that: it delivered 20 Gy of a
    prescribed 70 Gy (mucosal BED 24 instead of 84) and every radiotherapy
    number in the analysis was wrong by a factor of three.  Integrating
    segment by segment between breakpoints removes the failure mode entirely,
    and is also faster, because LSODA is then free to take large steps in
    between events.
    """
    bp = {0.0, t_end}
    for key in ("mel", "fu", "mtx", "cis", "pal", "pbm", "bzd", "gln",
                "cryo"):
        for w in sched.get(key, []):
            bp.add(float(w[0]))
            bp.add(float(w[1]))
    for t0, d, dur in sched.get("rt", []):
        bp.add(float(t0))
        bp.add(float(t0 + dur))
    return sorted(x for x in bp if 0.0 <= x <= t_end)


def simulate(sched, p=P, t_end=60.0, n=1201, y_init=None):
    p = dict(p)
    if y_init is None:
        y_init = equilibrate(p)
    # rescale D so that the homeostatic barrier is exactly 1.0
    Dss = y_init[IX["D"]]
    scale = 1.0 / Dss
    yi = y_init.copy()
    for k in ("P1", "P2", "P3", "D"):
        yi[IX[k]] *= scale
    pp = dict(p)
    pp["aS"] = p["aS"] * scale
    yi[IX["cAUCmuc"]:] = 0.0

    ts = np.linspace(0.0, t_end, n)
    bps = breakpoints(sched, t_end)
    Y = np.zeros((NST, n))
    y = yi.copy()
    filled = 0
    for a, b in zip(bps[:-1], bps[1:]):
        if b - a < 1e-12:
            continue
        last = b >= t_end - 1e-12
        if last:
            mask = (ts >= a - 1e-12) & (ts <= b + 1e-12)
        else:
            mask = (ts >= a - 1e-12) & (ts < b - 1e-12)
        te = ts[mask]
        # ALWAYS ask for the segment END as well.  With t_eval supplied,
        # sol.y[:, -1] is the last REQUESTED point, not the end of the
        # integration span -- carrying that forward silently froze the state
        # at the start of every short window, which is how a 70 Gy course
        # once delivered 0 Gy.
        teq = np.append(te, b) if (not len(te) or te[-1] < b - 1e-12) \
            else te
        sol = solve_ivp(rhs, (a, b), y, args=(pp, sched), t_eval=teq,
                        method="LSODA", rtol=1e-7, atol=1e-9)
        if not sol.success:
            raise RuntimeError("integration failed on [%g,%g]: %s"
                               % (a, b, sol.message))
        if len(te):
            Y[:, filled:filled + len(te)] = sol.y[:, :len(te)]
            filled += len(te)
        y = np.maximum(sol.y[:, -1], 0.0)
    if filled != n:
        raise RuntimeError("t_eval coverage %d/%d" % (filled, n))
    return post(ts, Y, pp, sched)


def post(t, Y, p, sch):
    """derive every reported quantity from the raw trajectory"""
    D = Y[IX["D"]]
    A = 1.0 / (1.0 + np.exp((D - p["Dcrit"]) / p["wD"]))
    NS = Y[IX["NS"]]
    Cmore = Y[IX["A_mor_e"]]
    opi = p["Emax_mor"] * (Cmore / (p["EC50_mor"] + np.maximum(Cmore, 0)))
    x = p["pain_k"] * A * (0.4 + 0.6 * NS)
    vas_raw = p["pain_max"] * x / (1.0 + x)
    vas = np.maximum(vas_raw - opi, 0.0)

    # WHO grade from ulcer area (0-4).  The scale is ORDINAL and SATURATES --
    # section 7 of the analysis uses exactly this to show why WHO-graded
    # trials lose assay sensitivity at the top of the range.
    who = np.zeros_like(A)
    ery = 1.0 / (1.0 + np.exp((D - 0.72) / 0.06))
    who = np.where(ery > 0.4, 1.0, 0.0)
    who = np.where(A > 0.08, 2.0, who)
    who = np.where(A > 0.35, 3.0, who)
    who = np.where(A > 0.72, 4.0, who)
    # OMAS: continuous 0-5 ulceration score, the sensitive instrument
    omas = 5.0 * A ** 0.8

    return dict(t=t, Y=Y, D=D, A=A, who=who, omas=omas, vas=vas,
                vas_raw=vas_raw, S=Y[IX["S"]], Sd=Y[IX["Sd"]],
                Sq=Y[IX["Sq"]],
                P=Y[IX["P1"]] + Y[IX["P2"]] + Y[IX["P3"]],
                TNF=Y[IX["TNF"]], IL1b=Y[IX["IL1b"]], IL6=Y[IX["IL6"]],
                NFkB=Y[IX["NFkB"]], CER=Y[IX["CER"]], MB=Y[IX["MB"]],
                ANC=Y[IX["Circ"]], KGFe=Y[IX["KGFe"]],
                Cm_mel=Y[IX["Cm_mel"]], Cm_5fu=Y[IX["Cm_5fu"]],
                Cm_mtx=Y[IX["Cm_mtx"]], Cm_cis=Y[IX["Cm_cis"]],
                Cpal=Y[IX["A_pal_c"]] / p["V_pal"],
                AUCmuc=Y[IX["cAUCmuc"]], ulcdays=Y[IX["cUlcD"]],
                opidays=Y[IX["cOpiD"]], BEDt=Y[IX["cBEDt"]],
                BEDm=Y[IX["cBEDm"]], infhaz=Y[IX["cInf"]],
                painAUC=Y[IX["cPain"]], p=p)


# ----------------------------------------------------------------------------
# summary metrics
# ----------------------------------------------------------------------------
def metrics(r, sev=3.0):
    t, who, A = r["t"], r["who"], r["A"]
    sevmask = who >= sev
    dt = t[1] - t[0]
    dur = float(sevmask.sum() * dt)
    onset = float(t[sevmask][0]) if sevmask.any() else float("nan")
    end = float(t[sevmask][-1]) if sevmask.any() else float("nan")
    return dict(
        peak_who=float(who.max()),
        peak_area=float(A.max()),
        peak_omas=float(r["omas"].max()),
        dur_sev=dur,
        onset_sev=onset,
        end_sev=end,
        incidence=bool(sevmask.any()),
        ulcdays=float(r["ulcdays"][-1]),
        opidays=float(r["opidays"][-1]),
        painAUC=float(r["painAUC"][-1]),
        nadir_anc=float(r["ANC"].min()),
        infhaz=float(r["infhaz"][-1]),
        auc_muc=float(r["AUCmuc"][-1]),
        BEDt=float(r["BEDt"][-1]),
        area_auc=float(np.trapezoid(A, t)),
    )


# ----------------------------------------------------------------------------
# canonical regimens
# ----------------------------------------------------------------------------
def sched_HDM(dose_mgm2=200.0, bsa=1.8, dur_h=0.5, t0=0.0):
    """high-dose melphalan, single 30-min infusion, day 0 (autologous HSCT)"""
    amt = dose_mgm2 * bsa * 1000.0            # ug
    dur = dur_h / 24.0
    s = dict(ZERO_SCHED)
    s = {k: (list(v) if isinstance(v, list) else v) for k, v in s.items()}
    s["mel"] = [(t0, t0 + dur, amt / dur)]
    return s


def sched_BEAM(bsa=1.8):
    """melphalan 140 mg/m2 within a BEAM-like conditioning block"""
    return sched_HDM(140.0, bsa)


def sched_5FU_bolus(dose_mgm2=425.0, bsa=1.8, days=(0, 1, 2, 3, 4)):
    amt = dose_mgm2 * bsa * 1000.0
    dur = 5.0 / (60.0 * 24.0)                 # 5-min push
    s = {k: (list(v) if isinstance(v, list) else v)
         for k, v in ZERO_SCHED.items()}
    s["fu"] = [(d, d + dur, amt / dur) for d in days]
    return s


def sched_5FU_CI(dose_mgm2_total=4000.0, bsa=1.8, dur_d=4.0, t0=0.0):
    amt = dose_mgm2_total * bsa * 1000.0
    s = {k: (list(v) if isinstance(v, list) else v)
         for k, v in ZERO_SCHED.items()}
    s["fu"] = [(t0, t0 + dur_d, amt / dur_d)]
    return s


def sched_chemoRT(total_Gy=70.0, nfx=35, fx_per_week=5, dose_per_fx=2.0,
                  cis_mgm2=100.0, cis_days=(0, 21, 42), bsa=1.8, t0=0.0,
                  gap=None, fx_per_day=1, bid_gap_h=6.0):
    """
    Definitive head-and-neck chemoradiation.

    fx_per_day > 1 gives genuine b.i.d. delivery with a bid_gap_h interfraction
    interval.  THIS IS NOT COSMETIC: the earlier version expressed
    "hyperfractionation" as fx_per_week = 10 with one fraction per day, which
    silently stretched 68 fractions over 68 DAYS instead of 34.  The overall
    treatment time is the whole point of an altered-fractionation comparison,
    so that arm was answering a different question from the one it was
    labelled with -- and it showed up as a LOWER tumour BED for the
    hyperfractionated arm, which is the opposite of the intent.
    """
    s = {k: (list(v) if isinstance(v, list) else v)
         for k, v in ZERO_SCHED.items()}
    rt = []
    day = t0
    given = 0
    days_per_week = min(fx_per_week, 7) if fx_per_day == 1 else fx_per_week
    while given < nfx:
        dow = int(round(day - t0)) % 7
        on = dow < days_per_week
        if gap is not None and gap[0] <= day - t0 < gap[1]:
            on = False
        if on:
            for j in range(fx_per_day):
                if given >= nfx:
                    break
                rt.append((day + j * bid_gap_h / 24.0, dose_per_fx,
                           10.0 / (60.0 * 24.0)))
                given += 1
        day += 1.0
        if day - t0 > 300:
            break
    s["rt"] = rt
    s["dose_per_fx"] = dose_per_fx
    if cis_mgm2:
        amt = cis_mgm2 * bsa * 1000.0
        dur = 2.0 / 24.0
        s["cis"] = [(t0 + d, t0 + d + dur, amt / dur) for d in cis_days]
    return s


def add_palifermin(s, days, dose_ug_kg=60.0, wt=75.0):
    """palifermin 60 ug/kg/d IV bolus on the listed days"""
    amt = dose_ug_kg * wt * 1000.0            # ng
    dur = 5.0 / (60.0 * 24.0)
    s = dict(s)
    s["pal"] = list(s.get("pal", [])) + [(d, d + dur, amt / dur) for d in days]
    return s


def add_cryo(s, windows):
    s = dict(s)
    s["cryo"] = list(s.get("cryo", [])) + list(windows)
    return s


def add_pbm(s, days, fluence=6.0):
    s = dict(s)
    dur = 5.0 / (60.0 * 24.0)
    s["pbm"] = list(s.get("pbm", [])) + [(d, d + dur, fluence / dur)
                                         for d in days]
    return s


def add_bzd(s, t0, t1, level=1.0):
    s = dict(s)
    s["bzd"] = list(s.get("bzd", [])) + [(t0, t1, level)]
    return s


def add_gln(s, t0, t1, level=1.0):
    s = dict(s)
    s["gln"] = list(s.get("gln", [])) + [(t0, t1, level)]
    return s


if __name__ == "__main__":
    r = simulate(sched_HDM(), t_end=40.0)
    m = metrics(r)
    print("HDM 200 mg/m2 reference run")
    for k, v in m.items():
        print("  %-12s %s" % (k, v))
