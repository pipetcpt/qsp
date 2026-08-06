#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
onfh_python_reference.py
========================
Independent Python/scipy re-implementation of the ONFH (osteonecrosis of the
femoral head) QSP model that is written in mrgsolve form in
`onfh_mrgsolve_model.R`.  Every number quoted in `README.md` is produced by
this file; run it and it writes `onfh_reference_output.txt`.

The purpose of the duplicate implementation is verification: the two files were
written from the same specification but not from each other, and disagreements
between them are model defects.  Defects found this way are listed at the end
of README.md.

------------------------------------------------------------------------------
THE THESIS THE MODEL ENCODES
------------------------------------------------------------------------------
The femoral head is a closed compartment in which three clocks run at once:

  clock 1  PERFUSION  (hours-days)  the ischaemic insult.  Over before the
                                    patient has symptoms, let alone a diagnosis.
  clock 2  REPAIR     (months)      creeping substitution.  Resorption precedes
                                    formation, so the interface between dead and
                                    living bone is WEAKER while it heals than it
                                    was while it was merely dead.
  clock 3  FATIGUE    (months-years) microdamage.  Living bone erases it by
                                    targeted remodelling; necrotic bone cannot,
                                    because targeted remodelling is signalled by
                                    osteocytes and there are none.

Collapse is a clock-3 event whose rate is set by clock 2 through one number:
the stiffness of the reparative interface.  The necrotic bone itself never
fails -- it is fully mineralised and as stiff as it ever was.

------------------------------------------------------------------------------
THE MECHANICAL CORE (why size AND location both matter, from geometry alone)
------------------------------------------------------------------------------
The lesion is a cone with its apex at the centre of the head, half-angle alpha,
axis tilted by phi from the hip joint resultant force.  Two areas follow:

    loaded cap area   A_cap  = 2 pi R^2 (1 - cos alpha)      grows like alpha^2
    conical interface A_int  =    pi R^2  sin alpha           grows like alpha

Every newton that lands on the necrotic cap must leave through the conical
interface, because the apex is a point.  So

    interface stress  ~  A_cap / A_int  =  2(1-cos a)/sin a  =  2 tan(alpha/2)

rises steeply and without bound as the lesion grows, while the material at the
interface -- and ONLY at the interface -- is the reparative tissue whose
strength dips during creeping substitution.  Location enters because only the
part of the cap inside the acetabular contact patch collects any load at all.

That single quotient is the whole of ONFH staging: Kerboull's necrotic angle
measures alpha, the JIC types measure where the cap sits relative to the
contact patch, and Steinberg's percent-volume measures a mixture of the two,
which is why it predicts worse than either.

------------------------------------------------------------------------------
WHAT IS FITTED
------------------------------------------------------------------------------
Three numbers, all stated in PARAMS with a `# FITTED` comment:
   k_dmg    microdamage rate constant   - one anchor: untreated JIC C1 hips
                                          reach 50% collapse near 5 years
   h0_tha   baseline THA hazard         - untreated post-collapse THA rate
   k_pain   pain scale factor           - baseline symptomatic VAS
Everything else is taken from published physiology / geometry / materials.
"""

import json
import math
import sys
import textwrap
from collections import OrderedDict

import numpy as np
from scipy.integrate import solve_ivp

RNG = np.random.default_rng(20260806)
DEG = math.pi / 180.0

# =============================================================================
# 1.  STATE VECTOR
# =============================================================================
STATES = [
    # --- drug / exposure PK (12) -------------------------------------------
    "GC_gut",    # 0   prednisolone-equivalent, gut depot            (mg)
    "GC_cen",    # 1   prednisolone-equivalent, central              (mg)
    "GC_eff",    # 2   marrow transcriptional exposure index         (mg/L-eq)
    "ALN_c",     # 3   alendronate, plasma                           (ug/L)
    "ALN_s",     # 4   alendronate, bone SURFACE (pharmacologic)     (umol/kg)
    "ALN_d",     # 5   alendronate, buried in bone (inactive)        (umol/kg)
    "STA_c",     # 6   statin exposure index                         (norm)
    "ENX",       # 7   enoxaparin anti-Xa                            (IU/mL)
    "ILO",       # 8   iloprost effect-site                          (norm)
    "TPT_sc",    # 9   teriparatide sc depot                         (ug)
    "TPT_c",     # 10  teriparatide plasma                           (pg/mL)
    "DMB_c",     # 11  denosumab plasma                              (ug/mL)
    # --- marrow lineage (5) -------------------------------------------------
    "PPARG",     # 12  adipogenic transcriptional drive              (norm, 1)
    "RUNX2",     # 13  osteogenic transcriptional drive              (norm, 1)
    "ADIPO",     # 14  marrow adipocyte volume fraction              (0-1)
    "MSC",       # 15  marrow stromal progenitor pool                (norm, 1)
    "OPRO",      # 16  osteoprogenitors available at the front       (norm, 1)
    # --- vascular / coagulation (5) ----------------------------------------
    "PAI1",      # 17  plasminogen activator inhibitor-1             (norm, 1)
    "THROMB",    # 18  microvascular occlusion fraction              (0-1)
    "ENDO",      # 19  endothelial/NO competence                     (0-1)
    "EDEMA",     # 20  marrow interstitial oedema fraction           (0-1)
    "ANGIO",     # 21  neovascular density at reparative front       (0-1)
    # --- compartment (2) ----------------------------------------------------
    "PIO",       # 22  intraosseous pressure of the head             (mmHg)
    "CDCH",      # 23  core-decompression channel patency            (0-1)
    # --- ischaemia / viability (3) -----------------------------------------
    "HEMV",      # 24  haematopoietic marrow viability               (0-1)
    "OCYV",      # 25  osteocyte viability in the at-risk zone       (0-1)
    "NECF",      # 26  necrotic volume fraction of the head          (0-1)
    # --- reparative front (6) ----------------------------------------------
    "XF",        # 27  creeping-substitution penetration depth       (mm)
    "CAV",       # 28  OPEN resorption cavity fraction at interface  (0-1)
    "NB",        # 29  new bone that has refilled cavities           (0-1)
    "MINZ",      # 30  mineralisation degree of that new bone        (0-1)
    "RIMD",      # 31  reactive sclerotic rim density index          (norm, 0)
    "FIBR",      # 32  fibrous/granulation fraction at the interface (0-1)
    # --- remodelling regulators (4) ----------------------------------------
    "RANKL",     # 33  RANKL                                         (norm, 1)
    "OPG",       # 34  osteoprotegerin                               (norm, 1)
    "OCL",       # 35  osteoclast activity at the interface          (norm, 1)
    "OBL",       # 36  osteoblast activity at the interface          (norm, 1)
    # --- structure and damage (8) ------------------------------------------
    "SUBPL",     # 37  subchondral plate integrity                   (0-1)
    "D_int",     # 38  microdamage, reparative interface             (0-1+)
    "D_nec",     # 39  microdamage, untouched necrotic bone          (0-1+)
    "D_plate",   # 40  microdamage, subchondral plate at lesion rim  (0-1+)
    "D_rim",     # 41  microdamage under the acetabular rim          (0-1+)
    "D_liv",     # 42  microdamage, adjacent living bone (control)   (0-1+)
    "CRESC",     # 43  subchondral fracture (crescent) extent        (0-1)
    "DEPR",      # 44  femoral head depression                       (mm)
    # --- joint / clinical (4) ----------------------------------------------
    "CART",      # 45  articular cartilage integrity                 (0-1)
    "SYNV",      # 46  synovitis / effusion                          (0-1)
    "PAINS",     # 47  slow pain state                               (0-10)
    "H_THA",     # 48  cumulative hazard of arthroplasty             (-)
]
NS = len(STATES)
IDX = {s: i for i, s in enumerate(STATES)}


# =============================================================================
# 2.  PARAMETERS
# =============================================================================
def base_params():
    p = OrderedDict()

    # ---- anatomy & materials (literature, not fitted) ----------------------
    p["R_head"] = 22.5        # mm     femoral head radius (44-48 mm diameter)
    p["t_sp"] = 1.10          # mm     subchondral bone plate thickness
    p["E_sp0"] = 5000.0       # MPa    apparent modulus, subchondral plate
    p["nu"] = 0.30            # -      Poisson ratio
    p["E_tr0"] = 620.0        # MPa    apparent modulus, femoral head trabeculae
    p["S_tr0"] = 9.5          # MPa    compressive strength, same, at BV/TV = 1
    p["S_pl0"] = 130.0        # MPa    bending strength of the subchondral plate
    # Contact pressure is applied to the JOINT SURFACE; the subchondral plate
    # and the trabecular architecture spread it before it reaches trabecular
    # tissue, so the apparent trabecular stress is a fraction of p(theta).  FE
    # studies of the normal hip put peak trabecular stresses near 1-1.5 MPa
    # against a strength near 9.5 MPa, i.e. a safety factor of 6-8.  The
    # interface traction sigma_int is NOT scaled by this: it is a global
    # equilibrium quantity (all of the plug's load must cross the cone), not a
    # local contact stress.
    p["f_spread"] = 0.30
    p["h_found"] = 6.0        # mm     Winkler foundation depth
    p["theta_c"] = 50.0       # deg    acetabular contact half-angle
    p["BW"] = 70.0            # kg     body weight
    p["f_hip"] = 2.80         # xBW    peak hip contact force, level walking
    p["Ncyc"] = 5000.0        # /day   gait cycles per hip
    # STRUCTURAL (not material) fatigue exponent.  Trabecular bone SPECIMENS
    # give m = 12-16.  A structure sheds load from failing elements to their
    # neighbours, which flattens the structural S-N curve, and the "stress"
    # here is a lumped equilibrium traction rather than a local tissue stress,
    # which flattens it further.  m is not fitted to any collapse rate; it was
    # chosen, together with f_spread, so that (a) normal bone does not fatigue-
    # fail in a heavy, active person and (b) latency to collapse spans months
    # rather than weeks.  Both are qualitative requirements, not data fits.  m
    # is carried in the sensitivity analysis: it is the least certain constant
    # in the model.
    p["m_fat"] = 4.0

    # ---- lesion geometry (per patient) ------------------------------------
    p["alpha_deg"] = 45.0     # deg    lesion cone half-angle (CNA = 4*alpha)
    p["phi_deg"] = 0.0        # deg    lesion axis tilt, +lateral
    p["lesion_preset"] = 1.0  # 1 = lesion present at t=0, 0 = incident model

    # ---- glucocorticoid PK (prednisolone equivalents) ---------------------
    p["GC_F"] = 0.80
    p["GC_ka"] = 48.0         # /d     ka 2 /h
    p["GC_CL"] = 210.0        # L/d    ~8.75 L/h
    p["GC_V"] = 35.0          # L
    p["GC_ke0"] = 0.15        # /d     marrow transcriptional integration
    p["GC_EC50"] = 0.220      # mg/L   effect-site EC50 for adipogenic drive;
                          #        set so that the half-maximal adipogenic
                          #        drive falls near 40 mg/d prednisolone,
                          #        the dose above which the epidemiology
                          #        finds the risk
    p["GC_gam"] = 2.5         # -      Hill; >1 is what makes PEAK dose matter
    p["GC_Emax"] = 3.2        # -      max fold-increase of adipogenic drive

    # ---- alendronate PK/PD -------------------------------------------------
    p["ALN_F"] = 0.0064       # -      oral bioavailability 0.6%
    p["ALN_ka"] = 24.0        # /d
    p["ALN_CL"] = 110.0       # L/d    renal
    p["ALN_V"] = 28.0         # L
    p["ALN_kupt"] = 3.4       # /d     plasma -> bone surface
    p["ALN_kbur"] = 0.0035    # /d     surface -> buried (t1/2 ~ 200 d)
    p["ALN_koff"] = 0.010     # /d     surface -> plasma (desorption)
    p["ALN_IC50"] = 0.55      # umol/kg on osteoclast activity
    p["ALN_Imax"] = 0.85

    # ---- statin, enoxaparin, iloprost, teriparatide, denosumab ------------
    p["STA_kel"] = 0.9        # /d     exposure-index turnover
    p["STA_IC50"] = 0.45      # on PPARG drive
    p["STA_Imax"] = 0.55
    p["STA_Erunx"] = 0.45     # BMP-2 mediated RUNX2 induction
    p["ENX_kel"] = 1.8        # /d
    p["ENX_EC50"] = 0.35      # IU/mL on thrombus dissolution
    p["ENX_Emax"] = 2.2
    p["ILO_kel"] = 0.55       # /d     effect-site (infusion is 5 d, effect weeks)
    p["ILO_Eendo"] = 0.55
    p["ILO_Eedema"] = 1.30
    p["TPT_ka"] = 26.0        # /d
    p["TPT_CL"] = 1900.0      # L/d
    p["TPT_V"] = 100.0        # L
    p["TPT_EC50"] = 12.0      # pg/mL x day-average -> anabolic index
    p["TPT_Eob"] = 1.35       # osteoblast activity
    p["TPT_Eocl"] = 0.35      # coupled osteoclast rise
    p["TPT_Emin"] = 0.55      # mineralisation rate
    p["DMB_kel"] = 0.028      # /d     t1/2 ~ 25 d
    p["DMB_IC50"] = 0.55      # ug/mL  on RANKL
    p["DMB_Imax"] = 0.92

    # ---- marrow lineage ----------------------------------------------------
    p["kin_pparg"] = 0.10
    p["kout_pparg"] = 0.10
    p["kin_runx2"] = 0.12
    p["kout_runx2"] = 0.12
    p["pparg_runx2_I"] = 0.60   # reciprocal lineage inhibition
    p["k_adipo"] = 0.020        # /d
    p["ADIPO0"] = 0.55          # adult femoral head marrow is largely fatty
    p["ADIPO_max"] = 0.90
    p["k_lipo"] = 0.020         # /d  lipolysis back toward ADIPO0
    p["ETOH"] = 0.0             # alcohol adipogenic drive (0-2)
    p["k_msc"] = 0.05
    p["k_opro"] = 0.06

    # ---- vascular / coagulation -------------------------------------------
    p["kin_pai1"] = 0.25
    p["kout_pai1"] = 0.25
    p["PAI1_gc"] = 1.40         # GC induction of PAI-1 (hypofibrinolysis)
    p["THROMBOPHILIA"] = 0.0    # inherited hypofibrinolysis multiplier (0-2)
    p["k_thr_on"] = 0.045       # /d
    p["k_thr_off"] = 0.030      # /d
    p["kin_endo"] = 0.12
    p["ENDO_gc"] = 0.45         # GC impairment of eNOS
    p["k_edema_on"] = 0.10
    p["k_edema_off"] = 0.035
    p["k_angio"] = 0.012        # /d  neovascularisation of the front
    p["angio_lag"] = 21.0       # d   ischaemia -> angiogenic response

    # ---- compartment / perfusion ------------------------------------------
    p["PIO0"] = 16.0            # mmHg  normal femoral head marrow pressure
    p["k_pio_fat"] = 55.0       # mmHg per unit normalised fat excess
    p["k_pio_ed"] = 25.0        # mmHg per unit oedema fraction (above ~31
                            # the oedema-pressure loop has gain > 1)
    p["k_pio_thr"] = 30.0       # mmHg per unit venous occlusion
    p["tau_pio"] = 3.0          # d
    p["k_cd"] = 0.55            # /d    decompression through a patent channel
    p["PIO_cd"] = 9.0           # mmHg  pressure a patent channel vents to
    p["k_close"] = 0.0115       # /d    channel closure, t1/2 ~ 60 d
    p["P_art"] = 65.0           # mmHg  subchondral arteriolar pressure
    p["P_ven"] = 12.0           # mmHg
    p["R_vasc"] = 1.0           # normalised
    p["thr_R"] = 4.0            # resistance rise per unit occlusion
    p["pO2_0"] = 35.0           # mmHg  normal marrow tissue pO2
    p["pO2_exp"] = 0.75

    # ---- viability ---------------------------------------------------------
    p["pO2_crit_hem"] = 15.0
    p["n_hem"] = 4.0
    p["k_hem"] = 3.0            # /d
    p["pO2_crit_ocy"] = 9.0
    p["n_ocy"] = 4.0
    p["k_ocy"] = 0.35           # /d
    p["k_ocy_rec"] = 0.05       # /d  recovery of stunned (not dead) osteocytes
    p["k_nec"] = 0.28           # /d  irreversible commitment to necrosis
    p["NEC_MAX"] = 0.42         # largest infarct the head geometry allows
                            # (alpha = 80.8 deg, Kerboull CNA 323)

    # ---- reparative front --------------------------------------------------
    p["v_front"] = 0.0247       # mm/d  creeping substitution ~0.75 mm/month
    p["L_eng"] = 3.0            # mm  front depth at which the interface is
                                #     fully converted to reparative tissue
    p["W_front"] = 2.0          # mm  thickness of the active reparative zone
    p["k_resp"] = 0.055         # /d  osteoclastic excavation at the interface
    p["RESP_max"] = 0.55        # max total turnover of interface bone volume
    p["k_fill"] = 0.030         # /d  refill of cavities by osteoblasts
    p["k_minz"] = 0.020         # /d  primary mineralisation (t1/2 ~ 35 d)
    p["k_fibr_on"] = 0.045
    p["k_fibr_off"] = 0.018
    p["FIBR_max"] = 0.30
    p["k_rim"] = 0.010          # /d  reactive sclerosis outside the lesion

    # ---- remodelling regulators -------------------------------------------
    p["kin_rankl"] = 0.20
    p["kout_rankl"] = 0.20
    p["RANKL_nec"] = 2.4        # necrotic debris drives RANKL at the interface
    p["kin_opg"] = 0.15
    p["kout_opg"] = 0.15
    p["OPG_runx2"] = 0.60
    p["k_ocl"] = 0.30           # /d
    p["k_obl"] = 0.22           # /d

    # ---- structure, damage, collapse --------------------------------------
    p["k_dmg"] = 0.90           # FITTED (one anchor: JIC C1 ~50% collapse @5y)
    p["k_dmg_rep"] = 0.0015     # /d  targeted removal of microdamage; a BMU
                                #     takes 3-4 months, so t1/2 ~ 460 d
    # Crack coalescence.  Kept deliberately weak.  A strong self term makes the
    # damage law explosive, so that every hip above the endurance limit fails
    # within weeks and time-to-collapse carries no information.  At 0.8 the
    # accumulation is close to linear, so WHETHER a hip collapses is set by the
    # integral of the damage rate over the vulnerable window and WHEN by its
    # magnitude -- two questions with two answers instead of one.
    p["c_Dself"] = 0.8
    p["zeta_max"] = 0.25        # lateral buttress load shedding
    p["dtheta_sh"] = 20.0       # deg over which the buttress engages
    p["chi_rim0"] = 1.35        # acetabular edge stress concentration
    p["c_rim_soft"] = 0.55      # extra edge loading when the rim sits on a plug
    p["c_depr"] = 0.55          # /mm  incongruity pressure rise after collapse
    p["k_cresc"] = 0.030        # /d
    p["k_depr"] = 0.055         # /d
    p["DEPR_max"] = 8.0         # mm
    p["COLLAPSE_MM"] = 2.0      # mm  radiographic collapse (ARCO IIIA/IIIB)
    p["k_subpl"] = 0.010        # /d  plate degradation once damaged

    # ---- joint / clinical --------------------------------------------------
    p["k_cart"] = 0.0022        # /d  (baseline term is now ~0: an
                            #      intact hip does not lose cartilage)
    p["k_synv_on"] = 0.030
    p["k_synv_off"] = 0.020
    p["k_pain"] = 0.16          # FITTED (baseline symptomatic VAS)
    p["tau_pain"] = 12.0        # d
    p["h0_tha"] = 4.20e-4       # FITTED (/d, arthroplasty hazard scale)
    p["a_pain"] = 1.5
    p["a_depr"] = 2.6
    p["a_cart"] = 2.2

    # ---- interventions (event flags filled by scenarios) -------------------
    p["CD_time"] = -1.0         # d   core decompression time (<0 = never)
    p["CD_ntrack"] = 1.0        # 1 = single 8-10 mm track, 3 = multiple small
    p["BMAC_dose"] = 0.0        # x10^4 CFU-F delivered
    p["SUPPORT_GRAFT"] = 0.0    # 0-1 structural support (tantalum/fibula)
    p["PWB"] = 1.0              # protected weight bearing multiplier on Ncyc
    return p


# =============================================================================
# 3.  STATIC GEOMETRY  (computed once per patient, not inside the ODE)
# =============================================================================
def hip_pressure_scale(p):
    """Peak contact pressure p0 (MPa) from body weight and contact geometry.

    Pressure field on the head surface: p(theta) = p0 cos(theta), theta < theta_c.
    Resultant along the load axis:
        F = 2 pi R^2 p0 (1 - cos^3 theta_c)/3
    """
    R = p["R_head"] * 1e-3                     # m
    tc = p["theta_c"] * DEG
    F = p["f_hip"] * p["BW"] * 9.80665         # N
    geom = 2.0 * math.pi * R * R * (1.0 - math.cos(tc) ** 3) / 3.0
    p0 = F / geom                              # Pa
    return p0 * 1e-6                           # MPa


def lesion_geometry(p, ngrid=241):
    """Static geometric quantities for the necrotic cone.

    Returns dict with
        A_cap   mm^2   total lesion cap area (spherical cap, half-angle alpha)
        A_int   mm^2   conical interface area = pi R^2 sin(alpha)
        L_nec   N      load collected by the loaded part of the lesion cap
        L_tot   N      total joint load (check: equals f_hip * BW * g)
        theta_L deg    angular position of the LATERAL lesion boundary
        f_load  -      fraction of the joint load landing on the lesion
        rim_on_lesion  0-1  is the acetabular rim (theta_c) over the lesion?
    """
    R = p["R_head"]
    a = p["alpha_deg"] * DEG
    phi = p["phi_deg"] * DEG
    tc = p["theta_c"] * DEG
    p0 = hip_pressure_scale(p)                 # MPa = N/mm^2

    A_cap = 2.0 * math.pi * R * R * (1.0 - math.cos(a))
    A_int = math.pi * R * R * math.sin(a)

    # Integrate p(theta) * cos(theta) over the sphere, restricted to
    # (i) the contact patch  theta < theta_c  and (ii) the lesion cone
    # (angular distance from the lesion axis < alpha).  Lesion axis lies in the
    # coronal plane at polar angle phi from the load axis.
    th = np.linspace(0.0, tc, ngrid)                       # polar from load axis
    ps = np.linspace(0.0, 2.0 * math.pi, 2 * ngrid, endpoint=False)
    TH, PS = np.meshgrid(th, ps, indexing="ij")
    # angular distance from the lesion axis
    cosd = (np.cos(TH) * math.cos(phi)
            + np.sin(TH) * math.sin(phi) * np.cos(PS))
    inside = cosd >= math.cos(a)
    pr = p0 * np.cos(TH)                                   # MPa
    dA = R * R * np.sin(TH)                                # * dth * dps
    integ = pr * np.cos(TH) * dA                           # axial component
    dth = th[1] - th[0]
    dps = ps[1] - ps[0]
    L_tot = float(np.sum(integ) * dth * dps)
    L_nec = float(np.sum(np.where(inside, integ, 0.0)) * dth * dps)

    theta_L = (p["phi_deg"] + p["alpha_deg"])              # lateral boundary
    # does the acetabular rim sit over necrotic bone?
    rim_on = 1.0 if (theta_L >= p["theta_c"]) else 0.0
    # graded version: how far past the rim
    rim_soft = float(np.clip((theta_L - p["theta_c"]) / 15.0, 0.0, 1.0))

    return dict(A_cap=A_cap, A_int=A_int, L_nec=L_nec, L_tot=L_tot,
                theta_L=theta_L, f_load=(L_nec / L_tot if L_tot > 0 else 0.0),
                p0=p0, rim_on=rim_on, rim_soft=rim_soft)


def jic_type(theta_L, theta_c=50.0):
    """Japanese Investigation Committee type from the lateral lesion boundary.

    The weight-bearing surface spans [-theta_c, +theta_c]; its medial third and
    medial two-thirds boundaries are therefore at -theta_c/3 and +theta_c/3.
    """
    if theta_L <= -theta_c / 3.0:
        return "A"
    if theta_L <= theta_c / 3.0:
        return "B"
    if theta_L <= theta_c:
        return "C1"
    return "C2"


# =============================================================================
# 4.  INITIAL CONDITIONS
# =============================================================================
def initial_state(p):
    y = np.zeros(NS)
    y[IDX["PPARG"]] = 1.0
    y[IDX["RUNX2"]] = 1.0
    y[IDX["ADIPO"]] = p["ADIPO0"]
    y[IDX["MSC"]] = 1.0
    y[IDX["OPRO"]] = 1.0
    y[IDX["PAI1"]] = 1.0
    y[IDX["ENDO"]] = 1.0
    y[IDX["PIO"]] = p["PIO0"]
    y[IDX["HEMV"]] = 1.0
    y[IDX["OCYV"]] = 1.0
    y[IDX["RANKL"]] = 1.0
    y[IDX["OPG"]] = 1.0
    y[IDX["OCL"]] = 1.0
    y[IDX["OBL"]] = 1.0
    y[IDX["SUBPL"]] = 1.0
    y[IDX["CART"]] = 1.0
    if p["lesion_preset"] > 0.5:
        # established lesion: necrosis complete, repair just beginning
        y[IDX["NECF"]] = (1.0 - math.cos(p["alpha_deg"] * DEG)) / 2.0
        y[IDX["OCYV"]] = 0.0
        y[IDX["HEMV"]] = 0.0
        y[IDX["ANGIO"]] = 0.25
        y[IDX["XF"]] = float(p.get("XF0", 0.5))
    return y


# =============================================================================
# 5.  ALGEBRAIC LAYER  (evaluated at every derivative call)
# =============================================================================
def algebra(t, y, p, geo):
    A = {}
    g = lambda k: y[IDX[k]]

    # ---------- drug effects ------------------------------------------------
    aln = max(g("ALN_s"), 0.0)
    A["I_aln"] = p["ALN_Imax"] * aln / (p["ALN_IC50"] + aln)
    sta = max(g("STA_c"), 0.0)
    A["I_sta"] = p["STA_Imax"] * sta / (p["STA_IC50"] + sta)
    A["E_sta_runx"] = p["STA_Erunx"] * sta / (p["STA_IC50"] + sta)
    enx = max(g("ENX"), 0.0)
    A["E_enx"] = p["ENX_Emax"] * enx / (p["ENX_EC50"] + enx)
    ilo = max(g("ILO"), 0.0)
    A["E_ilo_endo"] = p["ILO_Eendo"] * ilo / (0.5 + ilo)
    A["E_ilo_ed"] = p["ILO_Eedema"] * ilo / (0.5 + ilo)
    tpt = max(g("TPT_c"), 0.0)
    A["E_tpt"] = tpt / (p["TPT_EC50"] + tpt)
    dmb = max(g("DMB_c"), 0.0)
    A["I_dmb"] = p["DMB_Imax"] * dmb / (p["DMB_IC50"] + dmb)
    gce = max(g("GC_eff"), 0.0)
    A["E_gc"] = (p["GC_Emax"] * gce ** p["GC_gam"]
                 / (p["GC_EC50"] ** p["GC_gam"] + gce ** p["GC_gam"]))

    # ---------- perfusion (Starling resistor / vascular waterfall) ----------
    pio = max(g("PIO"), 0.0)
    p_out = max(p["P_ven"], pio)
    p_in = p["P_art"] * (0.55 + 0.45 * np.clip(g("ENDO"), 0.0, 1.5))
    R = p["R_vasc"] * (1.0 + p["thr_R"] * np.clip(g("THROMB"), 0.0, 1.0))
    Q = max(0.0, (p_in - p_out) / R)
    # Normalise to the BASELINE operating point, not to venous pressure.  With
    # the venous normalisation a healthy hip sat at Qrel = 0.925, which fed the
    # oedema term, which raised P_io, which lowered Qrel further: the loop gain
    # exceeded one and a normal femoral head slowly infarcted itself.
    Q0 = (p["P_art"] - max(p["P_ven"], p["PIO0"])) / p["R_vasc"]
    A["Qrel"] = Q / Q0
    A["pO2"] = p["pO2_0"] * A["Qrel"] ** p["pO2_exp"]

    # ---------- reparative-interface material ------------------------------
    # Volume bookkeeping inside the reparative interface tissue.  Three phases
    # sum to one:
    #     ORIG = 1 - CAV - NB   original (dead, fully mineralised) trabeculae
    #     CAV                   open resorption cavity (no load-bearing value)
    #     NB                    new bone, relative stiffness 0.35 + 0.65*MINZ
    # so the normalised load-bearing bone volume fraction is
    #     beta = 1 - CAV - 0.65*NB*(1 - MINZ)
    # which is 1 when nothing has happened and returns to 1 once every cavity
    # has been refilled and the new bone has mineralised.  It dips in between:
    # that dip is the whole disease.
    cav = np.clip(g("CAV"), 0.0, 0.95)
    nb = np.clip(g("NB"), 0.0, 0.95)
    fibr = np.clip(g("FIBR"), 0.0, 0.95)
    minz = np.clip(g("MINZ"), 0.0, 1.0)
    beta_rep = 1.0 - cav - 0.65 * nb * (1.0 - minz) - 0.25 * fibr * cav
    beta_rep = float(np.clip(beta_rep, 0.05, 1.25))
    A["beta_rep"] = beta_rep

    # how much of the conical interface has actually been converted from
    # "dead bone still welded to living bone" into "reparative tissue".
    # Before creeping substitution arrives the interface is intact and stiff.
    f_eng = float(np.clip(g("XF") / p["L_eng"], 0.0, 1.0))
    A["f_eng"] = f_eng
    beta_eff = 1.0 * (1.0 - f_eng) + beta_rep * f_eng
    A["beta_eff"] = beta_eff
    A["E_int"] = p["E_tr0"] * beta_eff ** 2.5
    A["S_int"] = p["S_tr0"] * beta_eff ** 2.0

    # ---------- load transfer ----------------------------------------------
    depr = max(g("DEPR"), 0.0)
    chi_depr = 1.0 + p["c_depr"] * depr           # incongruity concentrates load
    A["chi_depr"] = chi_depr

    # lateral buttress: living bone lateral to the lesion, still inside the
    # contact patch, sheds part of the load around the necrotic plug.
    dtheta = p["theta_c"] - geo["theta_L"]
    zeta = p["zeta_max"] * float(np.clip(dtheta / p["dtheta_sh"], 0.0, 1.0))
    A["zeta"] = zeta

    # a structural graft (tantalum rod, vascularised fibula) adds a parallel
    # load path straight to the neck
    supp = float(np.clip(p["SUPPORT_GRAFT"], 0.0, 1.0))

    L_eff = geo["L_nec"] * (1.0 - zeta) * chi_depr * (1.0 - 0.35 * supp)

    # core decompression removes bone from the interior of the head
    A["cd_defect"] = 0.0
    if p["CD_time"] >= 0.0 and t >= p["CD_time"]:
        # a single 8-10 mm track removes ~1.9% of the head cross-section on the
        # load path; multiple 3 mm drillings remove ~5x less per unit relief
        A["cd_defect"] = 0.075 / max(p["CD_ntrack"], 1.0)

    A_int_eff = geo["A_int"] * (1.0 - A["cd_defect"])
    A["sigma_int"] = L_eff / max(A_int_eff, 1.0)             # MPa
    A["sr_int"] = A["sigma_int"] / max(A["S_int"], 1e-6)

    # ---------- untouched necrotic bone (stiff, but cannot repair damage) ---
    # apparent compressive stress inside the plug: load / cap area
    A["sigma_nec"] = (geo["L_nec"] * chi_depr * p["f_spread"]) / max(geo["A_cap"], 1.0)
    A["sr_nec"] = A["sigma_nec"] / p["S_tr0"]

    # ---------- subchondral plate bending at the lesion rim ----------------
    # the plug sinks by w; the plate must accommodate that step over the
    # Winkler decay length of the surrounding healthy bone
    k_seg = A["E_int"] / max(min(max(g("XF"), 1.0), 6.0), 1.0)   # N/mm^3
    A["k_seg"] = k_seg
    w = A["sigma_int"] / max(k_seg, 1e-6) * 1.0                  # mm (order)
    A["w_sink"] = w
    D_pl = p["E_sp0"] * p["t_sp"] ** 3 / (12.0 * (1 - p["nu"] ** 2))
    k_h = p["E_tr0"] / p["h_found"]
    lam_h = (4.0 * D_pl / k_h) ** 0.25                            # healthy side
    lam_s = (4.0 * D_pl / max(k_seg, 1e-3)) ** 0.25               # lesion side
    lam = math.sqrt(lam_h * lam_s)     # the step is spanned across both
    A["lambda_pl"] = lam
    # the plate must accommodate a step w over the decay length lambda:
    # curvature ~ 2w/lambda^2, outer-fibre stress = E * (t/2) * curvature
    sig_pl = p["E_sp0"] * p["t_sp"] * w / (lam ** 2)
    A["sigma_pl"] = sig_pl
    A["sr_pl"] = sig_pl / max(p["S_pl0"] * np.clip(g("SUBPL"), 0.02, 1.0), 1e-6)

    # ---------- acetabular rim ---------------------------------------------
    chi_rim = p["chi_rim0"] + p["c_rim_soft"] * geo["rim_soft"]
    p0 = geo["p0"]
    A["sigma_rim"] = (p0 * math.cos(p["theta_c"] * DEG) * chi_rim * chi_depr
                      * p["f_spread"])
    # what material is under the rim?
    beta_rim = beta_eff if geo["rim_on"] > 0.5 else 1.0
    A["sr_rim"] = A["sigma_rim"] / max(p["S_tr0"] * beta_rim ** 2, 1e-6)
    A["rim_living"] = 0.0 if geo["rim_on"] > 0.5 else 1.0

    # ---------- living-bone control probe ----------------------------------
    # Negative control: bone OUTSIDE the lesion and outside the collapsing
    # segment, so it does not see the incongruity concentration chi_depr.
    A["sr_liv"] = (p0 * p["f_spread"]) / p["S_tr0"]

    # ---------- remodelling activity ---------------------------------------
    A["REM_int"] = min(np.clip(g("OCL"), 0.0, 5.0) / 2.5, 1.2) * A["f_eng"]
    A["REM_liv"] = 1.0
    return A


# =============================================================================
# 6.  RIGHT-HAND SIDE
# =============================================================================
def rhs(t, y, p, geo, dose):
    d = np.zeros(NS)
    A = algebra(t, y, p, geo)
    gv = lambda k: y[IDX[k]]

    # ---------------- PK ----------------------------------------------------
    d[IDX["GC_gut"]] = -p["GC_ka"] * gv("GC_gut")
    d[IDX["GC_cen"]] = (p["GC_F"] * p["GC_ka"] * gv("GC_gut")
                        - p["GC_CL"] / p["GC_V"] * gv("GC_cen"))
    d[IDX["GC_eff"]] = p["GC_ke0"] * (gv("GC_cen") / p["GC_V"] - gv("GC_eff"))

    d[IDX["ALN_c"]] = (-(p["ALN_CL"] / p["ALN_V"]) * gv("ALN_c")
                       - p["ALN_kupt"] * gv("ALN_c")
                       + p["ALN_koff"] * gv("ALN_s") * 10.0)
    d[IDX["ALN_s"]] = (p["ALN_kupt"] * gv("ALN_c") / 10.0
                       - p["ALN_kbur"] * gv("ALN_s")
                       - p["ALN_koff"] * gv("ALN_s"))
    d[IDX["ALN_d"]] = p["ALN_kbur"] * gv("ALN_s")

    d[IDX["STA_c"]] = -p["STA_kel"] * gv("STA_c")
    d[IDX["ENX"]] = -p["ENX_kel"] * gv("ENX")
    d[IDX["ILO"]] = -p["ILO_kel"] * gv("ILO")
    d[IDX["TPT_sc"]] = -p["TPT_ka"] * gv("TPT_sc")
    d[IDX["TPT_c"]] = (p["TPT_ka"] * gv("TPT_sc") * 1e6 / p["TPT_V"] / 1000.0
                       - p["TPT_CL"] / p["TPT_V"] * gv("TPT_c"))
    d[IDX["DMB_c"]] = -p["DMB_kel"] * gv("DMB_c")

    # ---------------- marrow lineage ---------------------------------------
    drive_ad = 1.0 + A["E_gc"] + p["ETOH"]
    d[IDX["PPARG"]] = (p["kin_pparg"] * drive_ad * (1.0 - A["I_sta"])
                       - p["kout_pparg"] * gv("PPARG"))
    d[IDX["RUNX2"]] = (p["kin_runx2"] * (1.0 + A["E_sta_runx"]
                                         + p["TPT_Eob"] * A["E_tpt"])
                       / (1.0 + p["pparg_runx2_I"] * max(gv("PPARG") - 1.0, 0.0))
                       - p["kout_runx2"] * gv("RUNX2"))
    d[IDX["ADIPO"]] = (p["k_adipo"] * max(gv("PPARG") - 1.0, 0.0)
                       * (p["ADIPO_max"] - gv("ADIPO"))
                       - p["k_lipo"] * (gv("ADIPO") - p["ADIPO0"]))
    d[IDX["MSC"]] = p["k_msc"] * (1.0 - gv("MSC")) - 0.02 * max(gv("PPARG") - 1, 0)
    d[IDX["OPRO"]] = (p["k_opro"] * (gv("MSC") * gv("RUNX2") - gv("OPRO")))

    # ---------------- vascular / coagulation -------------------------------
    d[IDX["PAI1"]] = (p["kin_pai1"] * (1.0 + p["PAI1_gc"] * A["E_gc"]
                                       + p["THROMBOPHILIA"])
                      - p["kout_pai1"] * gv("PAI1"))
    thr = np.clip(gv("THROMB"), 0.0, 1.0)
    d[IDX["THROMB"]] = (p["k_thr_on"] * max(gv("PAI1") - 1.0, 0.0)
                        * (1.0 - thr) * (2.0 - np.clip(gv("ENDO"), 0, 1))
                        - p["k_thr_off"] * (1.0 + A["E_enx"]) * thr)
    d[IDX["ENDO"]] = (p["kin_endo"]
                      * (1.0 + A["E_ilo_endo"] + 0.35 * A["E_sta_runx"]
                         - p["ENDO_gc"] * A["E_gc"] - gv("ENDO")))
    ed = np.clip(gv("EDEMA"), 0.0, 1.0)
    # oedema needs a real perfusion deficit, not any deficit at all
    isch_ed = max(1.0 - np.clip(A["Qrel"], 0.0, 1.0) - 0.10, 0.0) / 0.90
    d[IDX["EDEMA"]] = (p["k_edema_on"] * isch_ed * (0.55 - ed)
                       - p["k_edema_off"] * (1.0 + A["E_ilo_ed"]) * ed)
    isch = 1.0 - np.clip(A["Qrel"], 0.0, 1.0)
    d[IDX["ANGIO"]] = (p["k_angio"] * float(np.clip(gv("NECF") / 0.01, 0.0, 1.0))
                       * (1.0 + 0.8 * gv("OPRO"))
                       * (1.0 - np.clip(gv("ANGIO"), 0, 1))
                       - 0.004 * np.clip(gv("ANGIO"), 0, 1) * isch)

    # ---------------- compartment ------------------------------------------
    fat_excess = (gv("ADIPO") - p["ADIPO0"]) / max(1.0 - p["ADIPO0"], 1e-3)
    pio_t = (p["PIO0"] + p["k_pio_fat"] * max(fat_excess, 0.0)
             + p["k_pio_ed"] * ed + p["k_pio_thr"] * thr)
    d[IDX["PIO"]] = ((pio_t - gv("PIO")) / p["tau_pio"]
                     - p["k_cd"] * np.clip(gv("CDCH"), 0, 1)
                     * (gv("PIO") - p["PIO_cd"]))
    d[IDX["CDCH"]] = -p["k_close"] * gv("CDCH")

    # ---------------- viability --------------------------------------------
    h_hem = 1.0 / (1.0 + (A["pO2"] / p["pO2_crit_hem"]) ** p["n_hem"])
    h_ocy = 1.0 / (1.0 + (A["pO2"] / p["pO2_crit_ocy"]) ** p["n_ocy"])
    d[IDX["HEMV"]] = (-p["k_hem"] * h_hem * gv("HEMV")
                      + 0.08 * (1.0 - h_hem) * (1.0 - gv("HEMV")))
    d[IDX["OCYV"]] = (-p["k_ocy"] * h_ocy * gv("OCYV")
                      + p["k_ocy_rec"] * (1.0 - h_ocy)
                      * (1.0 - gv("NECF") / 0.5) * (1.0 - gv("OCYV")))
    if p["lesion_preset"] < 0.5:
        # h_ocy is the fraction of osteocyte territory sitting below the
        # critical oxygen tension, so NEC_MAX * h_ocy is the volume fraction
        # that the current depth of ischaemia can commit to necrosis.  The
        # ratchet (max(., 0)) makes the commitment irreversible: reperfusion
        # does not resurrect an infarct.  Making the CEILING depend on the
        # severity, rather than the RATE, is what lets a mild insult produce a
        # small lesion instead of the same maximal one.
        d[IDX["NECF"]] = p["k_nec"] * max(
            p["NEC_MAX"] * h_ocy - np.clip(gv("NECF"), 0.0, p["NEC_MAX"]), 0.0)
    else:
        d[IDX["NECF"]] = 0.0

    # ---------------- reparative front -------------------------------------
    ocl = np.clip(gv("OCL"), 0.0, 6.0)
    obl = np.clip(gv("OBL"), 0.0, 6.0)
    ang = np.clip(gv("ANGIO"), 0.0, 1.0)
    # the front has to cross the base of the necrotic cone: radius R sin(alpha)
    XFmax = max(2.0, p["R_head"] * math.sin(p["alpha_deg"] * DEG))
    adv = float(np.clip(1.0 - gv("XF") / XFmax, 0.0, 1.0))   # still advancing?
    d[IDX["XF"]] = p["v_front"] * ang * (0.25 + 0.75 * ocl) * adv
    cav = np.clip(gv("CAV"), 0.0, 0.95)
    nb = np.clip(gv("NB"), 0.0, 0.95)
    # osteoclasts excavate; osteoblasts refill.  Total turnover is capped at
    # RESP_max, and excavation stops once the front has finished crossing.
    exc = (p["k_resp"] * ocl * ang * adv
           * max(p["RESP_max"] - cav - nb, 0.0))
    fill = p["k_fill"] * obl * np.clip(gv("OPRO"), 0, 4) \
        * (1.0 - 0.8 * np.clip(gv("FIBR"), 0, 1)) * cav
    # RENEWAL.  Creeping substitution is a MOVING front: as it advances it
    # leaves finished bone behind and invades fresh necrotic bone ahead, so the
    # active zone is continuously reset to "not yet excavated".  Without this
    # term the interface is a single synchronised remodelling cohort that heals
    # itself in weeks and the disease disappears.  nu is the rate at which the
    # active zone of thickness W_front is replaced.
    nu = d[IDX["XF"]] / p["W_front"]
    d[IDX["CAV"]] = exc - fill - nu * cav
    d[IDX["NB"]] = fill - nu * nb
    d[IDX["MINZ"]] = (p["k_minz"] * (1.0 + p["TPT_Emin"] * A["E_tpt"])
                      * (1.0 - np.clip(gv("MINZ"), 0, 1))
                      - nu * np.clip(gv("MINZ"), 0, 1))
    d[IDX["RIMD"]] = p["k_rim"] * ang * (1.0 - np.clip(gv("RIMD"), 0, 1))
    d[IDX["FIBR"]] = (p["k_fibr_on"] * ang * (1.0 - 0.5 * np.clip(gv("OPRO"), 0, 2))
                      * (p["FIBR_max"] - np.clip(gv("FIBR"), 0, 1))
                      - p["k_fibr_off"] * (0.2 + nb) * np.clip(gv("FIBR"), 0, 1))

    # ---------------- remodelling regulators -------------------------------
    nec_debris = np.clip(gv("NECF") / 0.20, 0.0, 1.5) * A["f_eng"]
    d[IDX["RANKL"]] = (p["kin_rankl"] * (1.0 + p["RANKL_nec"] * nec_debris)
                       * (1.0 - A["I_dmb"])
                       - p["kout_rankl"] * gv("RANKL"))
    d[IDX["OPG"]] = (p["kin_opg"] * (1.0 + p["OPG_runx2"] * max(gv("RUNX2") - 1, 0))
                     - p["kout_opg"] * gv("OPG"))
    ratio = np.clip(gv("RANKL"), 0.01, 20.0) / np.clip(gv("OPG"), 0.05, 20.0)
    ocl_t = ratio * (1.0 - A["I_aln"]) * (1.0 + p["TPT_Eocl"] * A["E_tpt"])
    d[IDX["OCL"]] = p["k_ocl"] * (ocl_t - gv("OCL"))
    # coupling: osteoblast recruitment follows resorbed surface.  This is why
    # an antiresorptive cannot be a pure benefit - it removes the stimulus for
    # the formation it is trying to protect.
    coup = 0.30 + 0.70 * np.clip(gv("CAV"), 0, 1) / max(p["RESP_max"], 1e-3)
    obl_t = (np.clip(gv("RUNX2"), 0, 5) * np.clip(gv("OPRO"), 0, 5) * coup
             * (1.0 + p["TPT_Eob"] * A["E_tpt"]))
    d[IDX["OBL"]] = p["k_obl"] * (obl_t - gv("OBL"))

    # ---------------- microdamage ------------------------------------------
    m = p["m_fat"]
    N = p["Ncyc"] * p["PWB"]
    kd = p["k_dmg"]

    def dmg(state, sr, rem):
        """Fatigue microdamage.  Production follows a power-law S-N curve in
        the stress/strength ratio; removal is targeted remodelling, which
        requires living osteocytes to signal and is therefore ZERO inside
        necrotic bone.  N is normalised to 5000 gait cycles per hip per day."""
        D = max(y[IDX[state]], 0.0)
        prod = (kd * (N / 5000.0) * min(sr, 1.6) ** m
                * (1.0 + p["c_Dself"] * min(D, 1.0)))
        rem_t = p["k_dmg_rep"] * rem * D
        return prod - rem_t

    d[IDX["D_int"]] = dmg("D_int", A["sr_int"], A["REM_int"])
    d[IDX["D_nec"]] = dmg("D_nec", A["sr_nec"], 0.0)          # no living cells
    d[IDX["D_plate"]] = dmg("D_plate", A["sr_pl"], 0.25 * A["REM_int"])
    d[IDX["D_rim"]] = dmg("D_rim", A["sr_rim"], A["rim_living"])
    d[IDX["D_liv"]] = dmg("D_liv", A["sr_liv"], A["REM_liv"])

    # ---------------- collapse ---------------------------------------------
    Dmax = max(y[IDX["D_int"]], y[IDX["D_nec"]],
               y[IDX["D_plate"]], y[IDX["D_rim"]], 0.0)
    # failed fraction of the load path; exactly zero until damage approaches 1
    phi_fail = Dmax ** 8 / (1.0 + Dmax ** 8) if Dmax > 0 else 0.0
    cr = np.clip(gv("CRESC"), 0.0, 1.0)
    d[IDX["CRESC"]] = p["k_cresc"] * max(phi_fail - cr, 0.0) * 4.0
    d[IDX["DEPR"]] = (p["k_depr"] * cr * (0.35 + A["sr_int"])
                      * (p["DEPR_max"] - max(gv("DEPR"), 0.0)))
    d[IDX["SUBPL"]] = -p["k_subpl"] * cr * np.clip(gv("SUBPL"), 0, 1)

    # ---------------- joint / clinical -------------------------------------
    depr = max(gv("DEPR"), 0.0)
    # Cartilage is lost because the head has become incongruent, not because
    # time passes.  A baseline term here silently destroyed healthy hips.
    d[IDX["CART"]] = -p["k_cart"] * cr * (1.0 + depr / 2.0) \
        * np.clip(gv("CART"), 0, 1)
    d[IDX["SYNV"]] = (p["k_synv_on"] * (0.5 * cr + 0.5 * (1.0 - gv("CART")))
                      * (1.0 - np.clip(gv("SYNV"), 0, 1))
                      - p["k_synv_off"] * np.clip(gv("SYNV"), 0, 1))
    pain_t = 10.0 * (1.0 - math.exp(-p["k_pain"] * (
        0.09 * max(gv("PIO") - p["PIO0"], 0.0)
        + 3.0 * ed
        + 12.0 * cr
        + 2.2 * depr
        + 6.0 * (1.0 - np.clip(gv("CART"), 0, 1))
        + 4.0 * np.clip(gv("SYNV"), 0, 1))))
    d[IDX["PAINS"]] = (pain_t - gv("PAINS")) / p["tau_pain"]
    # Arthroplasty hazard is PROPORTIONAL to disease.  With an additive
    # constant an asymptomatic, spherical, painless hip accrued a 34% five-year
    # probability of being replaced.
    d[IDX["H_THA"]] = p["h0_tha"] * (
        p["a_pain"] * gv("PAINS") / 10.0
        + p["a_depr"] * depr / 4.0
        + p["a_cart"] * (1.0 - np.clip(gv("CART"), 0, 1)))

    return d


# =============================================================================
# 7.  DOSING EVENTS AND THE INTEGRATOR DRIVER
# =============================================================================
class Regimen:
    """A list of (time, state, amount) bolus events plus repeating schedules."""

    def __init__(self):
        self.events = []

    def add(self, time, state, amt):
        self.events.append((float(time), state, float(amt)))
        return self

    def repeat(self, start, stop, every, state, amt):
        t = float(start)
        while t < stop - 1e-9:
            self.add(t, state, amt)
            t += every
        return self

    def times(self):
        return sorted(set(e[0] for e in self.events))


def simulate(p, reg, tend, nout=400, geo=None):
    if geo is None:
        geo = lesion_geometry(p)
    y = initial_state(p)
    tout = np.linspace(0.0, tend, nout)
    Y = np.zeros((nout, NS))
    breaks = sorted(set([0.0] + [t for t in reg.times() if 0 <= t <= tend]
                        + ([p["CD_time"]] if 0 <= p["CD_time"] <= tend else [])
                        + [tend]))
    ptr = 0
    Y[0] = y
    for i in range(len(breaks) - 1):
        t0, t1 = breaks[i], breaks[i + 1]
        # apply bolus events at t0
        for (te, st, amt) in reg.events:
            if abs(te - t0) < 1e-9:
                y[IDX[st]] += amt
        if p["CD_time"] >= 0 and abs(p["CD_time"] - t0) < 1e-9:
            y[IDX["CDCH"]] = 1.0
            y[IDX["OPRO"]] += 0.9 * p["BMAC_dose"] / (1.0 + p["BMAC_dose"])
            y[IDX["ANGIO"]] = min(1.0, y[IDX["ANGIO"]] + 0.25)
        if t1 <= t0:
            continue
        sol = solve_ivp(rhs, (t0, t1), y, args=(p, geo, reg), method="LSODA",
                        rtol=1e-7, atol=1e-9, max_step=(t1 - t0),
                        dense_output=True)
        if not sol.success:
            raise RuntimeError(f"integration failed at t={t0}: {sol.message}")
        sel = (tout > t0) & (tout <= t1)
        if sel.any():
            Y[sel] = sol.sol(tout[sel]).T
        y = sol.y[:, -1].copy()
    return tout, Y, geo


def derived(tout, Y, p, geo):
    """Post-hoc observations, including the ones the model is judged on."""
    out = {}
    out["t"] = tout
    for s in STATES:
        out[s] = Y[:, IDX[s]]
    # re-evaluate the algebraic layer along the trajectory
    keys = ["sr_int", "sr_nec", "sr_pl", "sr_rim", "sr_liv", "sigma_int",
            "S_int", "beta_eff", "E_int", "Qrel", "pO2", "w_sink", "chi_depr",
            "zeta", "f_eng", "k_seg"]
    for k in keys:
        out[k] = np.array([algebra(tout[i], Y[i], p, geo)[k]
                           for i in range(len(tout))])
    out["S_collapse_free"] = np.exp(0.0)  # deterministic; see collapse_time
    out["HHS"] = (44.0 * (1.0 - out["PAINS"] / 10.0)
                  + 47.0 * np.clip(1.0 - out["DEPR"] / 6.0, 0, 1)
                  * np.clip(out["CART"], 0, 1)
                  + 9.0 * np.clip(1.0 - out["DEPR"] / 8.0, 0, 1))
    out["P_THA"] = 1.0 - np.exp(-out["H_THA"])
    return out


def simulate_incident(p, reg, tend, phi=0.0, t_stage1=540.0):
    """Two stages, because the lesion the insult creates is not known in
    advance: (1) integrate perfusion and cell death under the drug regimen to
    find how much of the head infarcts, (2) convert that volume fraction into
    a cone half-angle and integrate the mechanical consequences of THAT cone.

    Without this the incident scenarios all inherit whatever alpha_deg happened
    to be in the parameter list, and the steroid schedule cannot influence the
    outcome at all -- which is exactly the defect this function repairs.
    """
    p1 = OrderedDict(p)
    p1["lesion_preset"] = 0.0
    t1, Y1, g1 = simulate(p1, reg, t_stage1, nout=200)
    nf = float(Y1[-1, IDX["NECF"]])
    if nf < 5e-3:
        return dict(alpha=0.0, nf=nf, t1=t1, Y1=Y1, p1=p1, geo1=g1,
                    t2=None, Y2=None, p2=None, geo2=None, tcol=math.inf)
    alpha = math.degrees(math.acos(max(-1.0, 1.0 - 2.0 * nf)))
    p2 = OrderedDict(p)
    p2["lesion_preset"] = 1.0
    p2["alpha_deg"] = alpha
    p2["phi_deg"] = phi
    t2, Y2, g2 = simulate(p2, Regimen(), tend, nout=250)
    return dict(alpha=alpha, nf=nf, t1=t1, Y1=Y1, p1=p1, geo1=g1,
                t2=t2, Y2=Y2, p2=p2, geo2=g2,
                tcol=collapse_time(t2, Y2, p2))


def collapse_time(tout, Y, p):
    """First time femoral head depression reaches the radiographic threshold."""
    depr = Y[:, IDX["DEPR"]]
    idx = np.nonzero(depr >= p["COLLAPSE_MM"])[0]
    if idx.size == 0:
        return math.inf
    i = idx[0]
    if i == 0:
        return 0.0
    t0, t1 = tout[i - 1], tout[i]
    d0, d1 = depr[i - 1], depr[i]
    return t0 + (p["COLLAPSE_MM"] - d0) / max(d1 - d0, 1e-12) * (t1 - t0)


# =============================================================================
# 8.  REGIMEN BUILDERS
# =============================================================================
def reg_none():
    return Regimen()


def reg_steroid(daily_mg, days, start=0.0, taper_to=0.0, reg=None):
    """Oral prednisolone-equivalent, once daily, with an optional linear taper."""
    reg = reg or Regimen()
    n = int(round(days))
    for i in range(n):
        frac = 1.0 if taper_to <= 0 else (1.0 - (1.0 - taper_to) * i / max(n - 1, 1))
        reg.add(start + i, "GC_gut", daily_mg * frac)
    return reg


def reg_alendronate(start, days, reg=None):
    """70 mg once weekly oral."""
    reg = reg or Regimen()
    # 70 mg = 70/272.1 mmol = 0.257 mmol -> ug/L scale after F and V
    amt = 70.0 * 1e3 * 1.0  # ug
    reg.repeat(start, start + days, 7.0, "ALN_c",
               amt * 0.0064 / 28.0)  # ug/L into central
    return reg


def reg_statin(start, days, reg=None):
    reg = reg or Regimen()
    for i in range(int(days)):
        reg.add(start + i, "STA_c", 0.9)
    return reg


def reg_enoxaparin(start, days, reg=None):
    reg = reg or Regimen()
    for i in range(int(days)):
        reg.add(start + i, "ENX", 1.0)
    return reg


def reg_iloprost(start, days=5, reg=None):
    reg = reg or Regimen()
    for i in range(int(days)):
        reg.add(start + i, "ILO", 1.2)
    return reg


def reg_teriparatide(start, days, reg=None):
    reg = reg or Regimen()
    for i in range(int(days)):
        reg.add(start + i, "TPT_sc", 20.0)
    return reg


def reg_denosumab(start, days, reg=None):
    reg = reg or Regimen()
    reg.repeat(start, start + days, 182.5, "DMB_c", 60.0 * 1e3 / 3200.0)
    return reg


# =============================================================================
# 9.  SCENARIOS
# =============================================================================
def lesion_params(alpha=45.0, phi=0.0, **kw):
    p = base_params()
    p["alpha_deg"] = alpha
    p["phi_deg"] = phi
    p["lesion_preset"] = 1.0
    p.update(kw)
    return p


def jic_preset(name):
    """Representative geometry for each JIC class (mid-class by construction)."""
    presets = {
        # (alpha, phi) chosen so that theta_L = phi + alpha lands mid-class
        "A":  (28.0, -55.0),
        "B":  (34.0, -34.0),
        "C1": (45.0, -12.0),
        "C2": (58.0, 4.0),
    }
    a, ph = presets[name]
    return a, ph


SCENARIOS = OrderedDict()


def scenario(name, desc):
    def wrap(fn):
        SCENARIOS[name] = (desc, fn)
        return fn
    return wrap


@scenario("S01_healthy", "Healthy hip, no insult, 5 years")
def s01():
    p = base_params()
    p["lesion_preset"] = 0.0
    p["alpha_deg"] = 0.001
    return p, reg_none(), 1826.0


@scenario("S02_steroid_pulse",
          "High-dose steroid: 60 mg/d x 30 d then taper over 90 d (SLE/SARS-like)")
def s02():
    p = base_params()
    p["lesion_preset"] = 0.0
    r = reg_steroid(60.0, 30)
    r = reg_steroid(60.0, 90, start=30.0, taper_to=0.10, reg=r)
    return p, r, 1826.0


@scenario("S03_steroid_low_long",
          "Same cumulative dose, low peak: 12.4 mg/d x 145 d")
def s03():
    p = base_params()
    p["lesion_preset"] = 0.0
    # match cumulative dose of S02
    cum = 60.0 * 30 + sum(60.0 * (1.0 - 0.9 * i / 89) for i in range(90))
    r = reg_steroid(cum / 145.0, 145)
    return p, r, 1826.0


for _nm in ("A", "B", "C1", "C2"):
    def _mk(nm=_nm):
        @scenario(f"S0{4 + ['A', 'B', 'C1', 'C2'].index(nm)}_jic_{nm}",
                  f"Established JIC type {nm} lesion, untreated natural history")
        def _f():
            a, ph = jic_preset(nm)
            return lesion_params(a, ph), reg_none(), 1826.0
        return _f
    _mk()


@scenario("S08_C1_alendronate", "JIC C1 + alendronate 70 mg weekly x 24 months")
def s08():
    a, ph = jic_preset("C1")
    return lesion_params(a, ph), reg_alendronate(0.0, 730.0), 1826.0


@scenario("S09_C1_core_decompression", "JIC C1 + core decompression (single track)")
def s09():
    a, ph = jic_preset("C1")
    p = lesion_params(a, ph, CD_time=30.0, CD_ntrack=1.0, PWB=0.55)
    return p, reg_none(), 1826.0


@scenario("S10_C1_CD_BMAC",
          "JIC C1 + multiple drilling + bone marrow aspirate concentrate")
def s10():
    a, ph = jic_preset("C1")
    p = lesion_params(a, ph, CD_time=30.0, CD_ntrack=3.0, BMAC_dose=2.0, PWB=0.55)
    return p, reg_none(), 1826.0


@scenario("S11_C1_teriparatide", "JIC C1 + teriparatide 20 ug/d x 18 months")
def s11():
    a, ph = jic_preset("C1")
    return lesion_params(a, ph), reg_teriparatide(0.0, 548.0), 1826.0


@scenario("S12_C1_ALN_then_TPTD",
          "JIC C1 + alendronate 12 mo then teriparatide 18 mo")
def s12():
    a, ph = jic_preset("C1")
    r = reg_alendronate(0.0, 365.0)
    r = reg_teriparatide(365.0, 548.0, reg=r)
    return lesion_params(a, ph), r, 1826.0


@scenario("S13_steroid_plus_statin",
          "Steroid pulse + rosuvastatin prophylaxis from day 0 (Pritchett-like)")
def s13():
    p = base_params()
    p["lesion_preset"] = 0.0
    r = reg_steroid(60.0, 30)
    r = reg_steroid(60.0, 90, start=30.0, taper_to=0.10, reg=r)
    r = reg_statin(0.0, 365.0, reg=r)
    return p, r, 1826.0


@scenario("S14_steroid_plus_enoxaparin",
          "Steroid pulse + enoxaparin prophylaxis x 12 weeks (Glueck-like)")
def s14():
    p = base_params()
    p["lesion_preset"] = 0.0
    p["THROMBOPHILIA"] = 0.8
    r = reg_steroid(60.0, 30)
    r = reg_steroid(60.0, 90, start=30.0, taper_to=0.10, reg=r)
    r = reg_enoxaparin(0.0, 84.0, reg=r)
    return p, r, 1826.0


@scenario("S15_C2_core_decompression",
          "JIC C2 + core decompression (the sign flip)")
def s15():
    a, ph = jic_preset("C2")
    p = lesion_params(a, ph, CD_time=30.0, CD_ntrack=1.0, PWB=0.55)
    return p, reg_none(), 1826.0


@scenario("S16_B_core_decompression", "JIC B + core decompression")
def s16():
    a, ph = jic_preset("B")
    p = lesion_params(a, ph, CD_time=30.0, CD_ntrack=1.0, PWB=0.55)
    return p, reg_none(), 1826.0


@scenario("S17_C1_iloprost", "JIC C1 + iloprost 5-day infusion for marrow oedema")
def s17():
    a, ph = jic_preset("C1")
    return lesion_params(a, ph), reg_iloprost(14.0, 5), 1826.0


@scenario("S18_C1_denosumab", "JIC C1 + denosumab 60 mg q6m x 24 months")
def s18():
    a, ph = jic_preset("C1")
    return lesion_params(a, ph), reg_denosumab(0.0, 730.0), 1826.0


@scenario("S19_C1_graft", "JIC C1 + structural (tantalum/fibular) support")
def s19():
    a, ph = jic_preset("C1")
    p = lesion_params(a, ph, CD_time=30.0, CD_ntrack=3.0, SUPPORT_GRAFT=1.0,
                      BMAC_dose=2.0, PWB=0.50)
    return p, reg_none(), 1826.0


@scenario("S20_C1_PWB", "JIC C1 + protected weight bearing alone (crutches)")
def s20():
    a, ph = jic_preset("C1")
    return lesion_params(a, ph, PWB=0.45), reg_none(), 1826.0


# =============================================================================
# 10.  ANALYSES
# =============================================================================
def run_all_scenarios():
    res = OrderedDict()
    for name, (desc, fn) in SCENARIOS.items():
        p, reg, tend = fn()
        tout, Y, geo = simulate(p, reg, tend)
        obs = derived(tout, Y, p, geo)
        res[name] = dict(desc=desc, p=p, geo=geo, t=tout, Y=Y, obs=obs,
                         tcol=collapse_time(tout, Y, p))
    return res


def sample_patient(rng):
    """A virtual patient.

    Geometry is drawn wide enough to populate all four JIC classes; biological
    variability is deliberately SMALLER than geometric variability, because the
    clinical fact the model has to respect is that lesion size and position
    dominate every biological predictor that has ever been measured.

    XF0 is the front penetration already achieved when the patient presents.
    Real cohorts are enrolled at different points along the lesion's evolution;
    a hip diagnosed late with an intact head has already spent part of its
    vulnerable window, and has correspondingly less of it left.
    """
    b = base_params()
    alpha = float(np.clip(rng.normal(38.0, 14.0), 8.0, 78.0))
    phi = float(np.clip(rng.normal(-14.0, 20.0), -70.0, 34.0))
    p = lesion_params(alpha, phi)
    p["BW"] = float(np.clip(rng.normal(70.0, 11.0), 44.0, 112.0))
    p["Ncyc"] = float(np.clip(rng.lognormal(math.log(5000.0), 0.28), 900, 18000))
    p["S_tr0"] = b["S_tr0"] * float(np.clip(rng.normal(1.0, 0.10), 0.65, 1.4))
    p["v_front"] = b["v_front"] * float(np.clip(rng.normal(1.0, 0.20), 0.4, 1.9))
    p["RESP_max"] = float(np.clip(rng.normal(0.55, 0.045), 0.35, 0.72))
    p["k_fill"] = b["k_fill"] * float(np.clip(rng.normal(1.0, 0.15), 0.5, 1.7))
    # anatomical constraint: the cone has to fit inside the epiphysis, so its
    # medial boundary cannot swing past the head-neck junction.  Without this
    # the sampler produces "medial giant" lesions that do not exist.
    if (alpha + phi) - 2.0 * alpha < -105.0:
        alpha = min(alpha, (105.0 + (alpha + phi)) / 2.0)
        alpha = max(alpha, 8.0)
        p = lesion_params(alpha, phi, BW=p["BW"], Ncyc=p["Ncyc"],
                          S_tr0=p["S_tr0"], v_front=p["v_front"],
                          RESP_max=p["RESP_max"], k_fill=p["k_fill"])
    XFmax = max(2.0, p["R_head"] * math.sin(alpha * DEG))
    p["XF0"] = float(rng.uniform(0.3, 0.45 * XFmax))
    return p


def virtual_cohort(n, rng, mod=None, tend=1826.0, nout=120):
    """Simulate n virtual patients; return per-patient records."""
    recs = []
    for _ in range(n):
        p = sample_patient(rng)
        reg = Regimen()
        if mod is not None:
            p, reg = mod(p, reg)
        geo = lesion_geometry(p, ngrid=97)
        try:
            tout, Y, _ = simulate(p, reg, tend, nout=nout, geo=geo)
        except RuntimeError:
            continue
        tc = collapse_time(tout, Y, p)
        recs.append(dict(alpha=p["alpha_deg"], phi=p["phi_deg"],
                         theta_L=geo["theta_L"], jic=jic_type(geo["theta_L"]),
                         CNA=4.0 * p["alpha_deg"], tcol=tc,
                         f_load=geo["f_load"],
                         depr5=Y[-1, IDX["DEPR"]],
                         hhs5=(44.0 * (1 - Y[-1, IDX["PAINS"]] / 10.0)
                               + 47.0 * max(0.0, 1 - Y[-1, IDX["DEPR"]] / 6.0)
                               * Y[-1, IDX["CART"]]
                               + 9.0 * max(0.0, 1 - Y[-1, IDX["DEPR"]] / 8.0)),
                         ptha5=1 - math.exp(-Y[-1, IDX["H_THA"]])))
    return recs


def collapse_rate(recs, horizon):
    if not recs:
        return float("nan")
    return sum(1 for r in recs if r["tcol"] <= horizon) / len(recs)


# =============================================================================
# 11.  REPORT
# =============================================================================
def pre_peak(obs, tout, tcol, key, cap=550.0):
    """Peak of an observable BEFORE the head collapses.

    Post-collapse the incongruity factor chi_depr multiplies every stress by up
    to 5.4, so an unrestricted maximum compares a collapsed hip with an intact
    one and reports nonsense (it made alendronate look 85% protective).
    """
    lim = min(cap, tcol if not math.isinf(tcol) else cap)
    sel = tout <= max(lim, 1.0)
    return float(np.max(obs[key][sel]))


def fmt(x, n=3):
    if x is None:
        return "  n/a"
    if isinstance(x, float) and (math.isinf(x) or math.isnan(x)):
        return "  inf" if math.isinf(x) else "  nan"
    return f"{x:.{n}f}"


def main():
    L = []
    w = L.append
    w("=" * 79)
    w("ONFH QSP MODEL - PYTHON REFERENCE OUTPUT")
    w("osteonecrosis of the femoral head; independent re-implementation")
    w("=" * 79)
    w("")

    # ---------------------------------------------------------------- 0
    p = base_params()
    w("[0] DERIVED MECHANICAL CONSTANTS (no fitting)")
    w("-" * 79)
    p0 = hip_pressure_scale(p)
    w(f"  body weight                          {p['BW']:.0f} kg")
    w(f"  peak hip contact force  {p['f_hip']:.2f} x BW  = "
      f"{p['f_hip'] * p['BW'] * 9.80665:.0f} N")
    w(f"  contact half-angle                   {p['theta_c']:.0f} deg")
    w(f"  => peak contact pressure p0          {p0:.3f} MPa")
    w("     (in-vivo instrumented-endoprosthesis gait values 2-5 MPa)")
    D_pl = p["E_sp0"] * p["t_sp"] ** 3 / (12 * (1 - p["nu"] ** 2))
    k_h = p["E_tr0"] / p["h_found"]
    lam = (4 * D_pl / k_h) ** 0.25
    w(f"  subchondral plate flexural rigidity  {D_pl:.1f} N.mm")
    w(f"  healthy Winkler foundation modulus   {k_h:.1f} N/mm^3")
    w(f"  Winkler decay length lambda          {lam:.3f} mm")
    w("")
    w("  The geometric amplification A_cap/A_int = 2 tan(alpha/2):")
    w("     alpha    CNA     A_cap      A_int      A_cap/A_int   theta_L   JIC")
    for a in (15, 25, 35, 45, 55, 65, 75):
        pp = base_params(); pp["alpha_deg"] = a; pp["phi_deg"] = 0.0
        g = lesion_geometry(pp, ngrid=121)
        w(f"     {a:3d}     {4*a:4d}   {g['A_cap']:8.1f}   {g['A_int']:8.1f}"
          f"      {g['A_cap']/g['A_int']:6.3f}      {g['theta_L']:5.1f}   "
          f"{jic_type(g['theta_L'])}")
    w("")

    # ---------------------------------------------------------------- 1
    w("[1] SCENARIOS")
    w("-" * 79)
    res = run_all_scenarios()
    w(f"{'scenario':<26}{'JIC':>4}{'CNA':>6}{'collapse(mo)':>13}"
      f"{'depr5y':>8}{'HHS5y':>7}{'P(THA)5y':>10}")
    for name, r in res.items():
        g = r["geo"]
        incident = r["p"]["lesion_preset"] < 0.5 and r["p"]["alpha_deg"] > 1
        jt = jic_type(g["theta_L"]) if r["p"]["alpha_deg"] > 1 else "-"
        tc = r["tcol"]
        tcs = "never" if math.isinf(tc) else f"{tc/30.44:.1f}"
        if incident:
            jt, tcs = "->[10]", "see [10]"
        w(f"{name:<26}{jt:>4}{4*r['p']['alpha_deg']:6.0f}{tcs:>13}"
          f"{r['obs']['DEPR'][-1]:8.2f}{r['obs']['HHS'][-1]:7.1f}"
          f"{r['obs']['P_THA'][-1]:10.3f}")
    w("")

    # ---------------------------------------------------------------- 2
    w("[2] THE VULNERABLE WINDOW  (untreated JIC C1)")
    w("-" * 79)
    r = res["S06_jic_C1"]
    w("  month   XF(mm)   CAV     NB   MINZ   beta_eff  E_int(MPa)  S_int(MPa)"
      "  sigma_int  sr_int")
    w(f"  (this hip collapses at month {r['tcol']/30.44:.1f}; rows marked * are"
      f" AFTER collapse, where the incongruity factor inflates every stress)")
    for mo in (0, 3, 6, 9, 12, 18, 24, 36, 48, 60):
        i = int(np.argmin(np.abs(r["t"] - mo * 30.44)))
        o = r["obs"]
        mk = "*" if r["t"][i] > r["tcol"] else " "
        w(f" {mk}{mo:5d}  {o['XF'][i]:6.2f} {o['CAV'][i]:6.3f} {o['NB'][i]:6.3f}"
          f" {o['MINZ'][i]:6.3f}   {o['beta_eff'][i]:7.3f}  {o['E_int'][i]:9.1f}"
          f"  {o['S_int'][i]:9.3f}  {o['sigma_int'][i]:9.3f} {o['sr_int'][i]:7.3f}")
    o = r["obs"]
    itr = int(np.argmin(o["beta_eff"]))
    w("")
    w(f"  interface stiffness minimum at month {r['t'][itr]/30.44:.1f} "
      f"(beta_eff {o['beta_eff'][itr]:.3f}, E_int {o['E_int'][itr]:.1f} MPa, "
      f"{100*o['E_int'][itr]/p['E_tr0']:.1f}% of normal)")
    w("  Before creeping substitution arrives the necrotic bone is fully")
    w("  mineralised and mechanically intact.  The hip becomes vulnerable")
    w("  BECAUSE it started to heal.")
    w("")

    # ---------------------------------------------------------------- 3
    w("[3] WHICH STRUCTURE FAILS FIRST (untreated JIC C1, stress/strength)")
    w("-" * 79)
    w("  month   interface  necrotic-core  subchondral-plate  acetabular-rim"
      "  living-bone")
    for mo in (0, 6, 12, 18, 24, 36, 60):
        i = int(np.argmin(np.abs(r["t"] - mo * 30.44)))
        w(f"  {mo:5d}   {o['sr_int'][i]:9.3f}  {o['sr_nec'][i]:13.3f}"
          f"  {o['sr_pl'][i]:17.3f}  {o['sr_rim'][i]:14.3f}"
          f"  {o['sr_liv'][i]:11.3f}")
    w("")
    w("  Damage at 5 years (1.0 = failure):")
    for s, lab in (("D_int", "reparative interface"), ("D_nec", "necrotic core"),
                   ("D_plate", "subchondral plate"), ("D_rim", "acetabular rim"),
                   ("D_liv", "living bone")):
        w(f"     {lab:<24}{r['obs'][s][-1]:12.4g}")
    w("  The interface fails first, by orders of magnitude.  Collapse in this")
    w("  model is a failure of the REPAIR, not of the necrosis.")
    w("")

    # ---------------------------------------------------------------- 4
    w("[4] VIRTUAL COHORT: JIC CLASS COLLAPSE RATES")
    w("-" * 79)
    rng = np.random.default_rng(7)
    cohort = virtual_cohort(360, rng)
    w(f"  n = {len(cohort)} virtual hips")
    w(f"{'JIC':>5}{'n':>6}{'2y':>9}{'3y':>9}{'5y':>9}"
      f"{'median t_col|collapsed':>24}")
    for jt in ("A", "B", "C1", "C2"):
        sub = [x for x in cohort if x["jic"] == jt]
        if not sub:
            w(f"{jt:>5}{0:>6}")
            continue
        cols = sorted(x["tcol"] for x in sub if not math.isinf(x["tcol"]))
        meds = f"{cols[len(cols)//2]/30.44:.1f}" if cols else "-"
        w(f"{jt:>5}{len(sub):>6}{collapse_rate(sub, 730):>9.3f}"
          f"{collapse_rate(sub, 1096):>9.3f}{collapse_rate(sub, 1826):>9.3f}"
          f"{meds:>24}")
    w("")
    w("  Reported natural history (Japanese Investigation Committee /")
    w("  Nishii, Sugano; Mont 2010 systematic review): type A rarely collapses,")
    w("  type B ~15-25%, C1 ~40-60%, C2 ~70-85% within 3-5 years.")
    w("  ONLY the C1 rate was used to fit k_dmg; A, B and C2 are predictions.")
    w("")
    w("  Timing of collapse among hips that collapse within 5 years:")
    allc = [x["tcol"] for x in cohort if x["tcol"] <= 1826]
    if allc:
        allc = np.array(allc)
        w(f"     within 12 months  {np.mean(allc <= 365):.3f}")
        w(f"     within 24 months  {np.mean(allc <= 730):.3f}")
        w(f"     within 36 months  {np.mean(allc <= 1096):.3f}")
        w("     (clinically, most collapses occur in the first two years and a")
        w("      hip intact at 3-4 years rarely collapses later)")
    w("")

    # ---------------------------------------------------------------- 5
    w("[5] THE NECROTIC-ANGLE THRESHOLD, DERIVED NOT ASSUMED")
    w("-" * 79)
    w("  Kerboull combined necrotic angle (CNA) = AP arc + lateral arc = 4*alpha")
    w("  Sweep CNA with the lesion axis on the load axis; report 5-year collapse.")
    w(f"{'CNA(deg)':>10}{'alpha':>8}{'theta_L':>9}{'JIC':>5}"
      f"{'t_collapse(mo)':>16}{'depr5y(mm)':>12}")
    sweep = []
    for cna in range(80, 301, 15):
        a = cna / 4.0
        pp = lesion_params(a, 0.0)
        tt, YY, gg = simulate(pp, reg_none(), 1826.0, nout=200)
        tc = collapse_time(tt, YY, pp)
        sweep.append((cna, tc))
        tcs = "never" if math.isinf(tc) else f"{tc/30.44:.1f}"
        w(f"{cna:>10}{a:>8.1f}{gg['theta_L']:>9.1f}"
          f"{jic_type(gg['theta_L']):>5}{tcs:>16}{YY[-1, IDX['DEPR']]:>12.2f}")
    thr = None
    for i in range(1, len(sweep)):
        if math.isinf(sweep[i - 1][1]) and not math.isinf(sweep[i][1]):
            thr = (sweep[i - 1][0], sweep[i][0])
    w("")
    if thr:
        w(f"  deterministic collapse threshold between CNA {thr[0]} and {thr[1]} deg")
    w("  Kerboull's published risk bands: <190 low, 190-240 intermediate,")
    w("  >240 high.  The model was never shown those numbers.")
    w("")
    # population version of the same sweep
    w("  Population version (5-year collapse probability by CNA band):")
    for lo, hi in ((0, 150), (150, 190), (190, 240), (240, 400)):
        sub = [x for x in cohort if lo <= x["CNA"] < hi]
        if sub:
            w(f"     CNA {lo:>3}-{hi:<3}  n={len(sub):<4} "
              f"collapse@5y = {collapse_rate(sub, 1826):.3f}")
    w("")

    # ---------------------------------------------------------------- 6
    w("[6] SIZE VERSUS LOCATION: WHY VOLUME-BASED STAGING PREDICTS WORSE")
    w("-" * 79)
    w("  Two lesions with the SAME cone half-angle (same Steinberg volume),")
    w("  differing only in where the cone points.")
    w(f"{'phi(deg)':>10}{'theta_L':>9}{'JIC':>5}{'f_load':>9}"
      f"{'sigma_int@12mo':>16}{'t_collapse(mo)':>16}")
    for phi in (-55, -45, -35, -25, -15, -5, 5, 15):
        pp = lesion_params(45.0, phi)
        tt, YY, gg = simulate(pp, reg_none(), 1826.0, nout=200)
        oo = derived(tt, YY, pp, gg)
        i12 = int(np.argmin(np.abs(tt - 365)))
        tc = collapse_time(tt, YY, pp)
        tcs = "never" if math.isinf(tc) else f"{tc/30.44:.1f}"
        w(f"{phi:>10}{gg['theta_L']:>9.1f}{jic_type(gg['theta_L']):>5}"
          f"{gg['f_load']:>9.3f}{oo['sigma_int'][i12]:>16.3f}{tcs:>16}")
    w("")

    # ---------------------------------------------------------------- 7
    w("[7] BISPHOSPHONATE: WHAT IS THE MECHANISTIC CEILING?")
    w("-" * 79)
    w("  Alendronate acts on ONE limb of the interface: it suppresses the")
    w("  osteoclastic excavation that makes the interface porous.  It cannot")
    w("  make bone.  And because osteoblast recruitment in this model follows")
    w("  the resorbed surface (the coupling term), suppressing resorption also")
    w("  removes part of the stimulus for the formation it is protecting.")
    ra = res["S08_C1_alendronate"]
    w("   month    CAV(ctrl)   CAV(ALN)   NB(ctrl)   NB(ALN)  beta(ctrl)  beta(ALN)")
    for mo in (3, 6, 12, 18, 24, 36):
        i = int(np.argmin(np.abs(r["t"] - mo * 30.44)))
        j = int(np.argmin(np.abs(ra["t"] - mo * 30.44)))
        w(f"   {mo:5d}   {r['obs']['CAV'][i]:9.3f}  {ra['obs']['CAV'][j]:9.3f}"
          f"  {r['obs']['NB'][i]:9.3f} {ra['obs']['NB'][j]:9.3f}"
          f"  {r['obs']['beta_eff'][i]:10.3f} {ra['obs']['beta_eff'][j]:10.3f}")
    b_c = r["obs"]["beta_eff"].min()
    b_a = ra["obs"]["beta_eff"].min()
    sr_c = pre_peak(r["obs"], r["t"], r["tcol"], "sr_int")
    sr_a = pre_peak(ra["obs"], ra["t"], ra["tcol"], "sr_int")
    m = p["m_fat"]
    w("")
    w(f"  minimum interface beta   control {b_c:.3f}   alendronate {b_a:.3f}"
      f"   ({100*(b_a-b_c)/b_c:+.1f}%)")
    w(f"  peak stress/strength     control {sr_c:.3f}   alendronate {sr_a:.3f}"
      f"   ({100*(sr_a-sr_c)/sr_c:+.1f}%)")
    w(f"  implied damage-rate ratio (sr_ALN/sr_ctrl)^m, m = {m:.1f}:"
      f"  {(sr_a/sr_c)**m:.3f}")
    w("")
    w("  The CEILING.  Even at complete osteoclast blockade the excavation")
    w("  term cannot go below zero, so there is a hardest floor on the")
    w("  interface porosity, and therefore a hardest ceiling on the benefit:")
    for imax in (0.0, 0.5, 0.85, 0.999):
        pp = lesion_params(*jic_preset("C1"))
        pp["ALN_Imax"] = imax
        tt, YY, gg = simulate(pp, reg_alendronate(0.0, 730.0), 1100.0, nout=160)
        oo = derived(tt, YY, pp, gg)
        srm = pre_peak(oo, tt, collapse_time(tt, YY, pp), "sr_int")
        w(f"     osteoclast inhibition Imax = {imax:.3f}   beta_min "
          f"{oo['beta_eff'].min():.3f}   peak sr {srm:.3f}   "
          f"rate ratio {(srm/sr_c)**m:.3f}")
    w("")
    w("  Trial-scale simulation, paired virtual patients, enrolment")
    w("  restricted to large lesions (CNA >= 160 deg) as both trials did:")
    rng2 = np.random.default_rng(11)

    def mod_aln(pp, rr):
        return pp, reg_alendronate(0.0, 730.0)

    def big_only(recs):
        return [x for x in recs if x["CNA"] >= 160]
    ctrl = virtual_cohort(180, np.random.default_rng(11), tend=760, nout=60)
    trt = virtual_cohort(180, np.random.default_rng(11), mod=mod_aln,
                         tend=760, nout=60)
    cb, tb = big_only(ctrl), big_only(trt)
    rc = collapse_rate(cb, 730)
    rt = collapse_rate(tb, 730)
    w(f"     control      n={len(cb):<4} collapse @24 mo = {rc:.3f}")
    w(f"     alendronate  n={len(tb):<4} collapse @24 mo = {rt:.3f}")
    w(f"     relative risk = {rt/max(rc,1e-9):.3f}"
      f"   absolute risk reduction {rc-rt:+.3f}")
    w("")
    w("     Chen 2012 (RCT, double-blind, ARCO IIC/IIIC, 2 y): 16/32 vs 17/33")
    w("                                       = 0.500 vs 0.515,  RR 0.97")
    w("     Lai  2005 (Steinberg IIC/IIIC, 2 y):                  2/29 vs 19/25")
    w("                                       = 0.069 vs 0.760,  RR 0.09")
    w("     The two trials disagree by a factor of ten in RR.  The model")
    w("     cannot produce Lai's effect size at ANY degree of osteoclast")
    w("     blockade, because blockade cannot lift beta above the value it")
    w("     would have had with no resorption at all.  It sits close to Chen.")
    w("")

    # ---------------------------------------------------------------- 8
    w("[8] CORE DECOMPRESSION CHANGES SIGN WITH LESION SIZE")
    w("-" * 79)
    w("  Decompression does three things at once: it vents the compartment")
    w("  (good, but the ischaemic insult is over), it opens a channel for")
    w("  repair (good), and it removes load-bearing bone (bad).  The third")
    w("  scales with the lesion; the first two do not.")
    w(f"{'CNA':>6}{'JIC':>5}{'no CD (mo)':>13}{'CD (mo)':>11}"
      f"{'delta (mo)':>12}{'depr5y noCD':>13}{'depr5y CD':>11}{'verdict':>9}")
    for cna in (110, 125, 135, 145, 155, 170, 185, 200, 225, 250):
        a = cna / 4.0
        p1 = lesion_params(a, 0.0)
        t1, Y1, g1 = simulate(p1, reg_none(), 1826.0, nout=200)
        p2 = lesion_params(a, 0.0, CD_time=30.0, CD_ntrack=1.0, PWB=0.55)
        t2, Y2, g2 = simulate(p2, reg_none(), 1826.0, nout=200)
        c1 = collapse_time(t1, Y1, p1)
        c2 = collapse_time(t2, Y2, p2)
        s1 = "never" if math.isinf(c1) else f"{c1/30.44:.1f}"
        s2 = "never" if math.isinf(c2) else f"{c2/30.44:.1f}"
        if math.isinf(c1) and math.isinf(c2):
            dl, verdict = "-", "both OK"
        elif math.isinf(c2):
            dl, verdict = "rescued", "HELP"
        elif math.isinf(c1):
            dl, verdict = "caused", "HARM"
        else:
            dl = f"{(c2 - c1)/30.44:+.1f}"
            verdict = "help" if c2 > c1 + 0.5 else ("harm" if c2 < c1 - 0.5 else "=")
        w(f"{cna:>6}{jic_type(g1['theta_L']):>5}{s1:>13}{s2:>11}{dl:>12}"
          f"{Y1[-1, IDX['DEPR']]:>13.2f}{Y2[-1, IDX['DEPR']]:>11.2f}{verdict:>9}")
    w("")

    # ---------------------------------------------------------------- 9
    w("[9] STEROID: PEAK DOSE VERSUS CUMULATIVE DOSE")
    w("-" * 79)
    w("  Cumulative prednisolone-equivalent held constant at 4650 mg; only the")
    w("  schedule changes.  Each arm is simulated in TWO stages: the insult")
    w("  decides how much of the head infarcts, and the cone that produces")
    w("  then has to carry the load.")
    w(f"{'daily mg':>9}{'days':>6}{'peakGCeff':>11}{'ADIPOmax':>10}"
      f"{'PIOmax':>8}{'Qrelmin':>9}{'NECF':>8}{'alpha':>7}{'CNA':>6}"
      f"{'JIC':>5}{'collapse(mo)':>13}")
    CUM = 4650.0
    for daily in (100, 60, 40, 25, 15, 10, 6):
        days = int(round(CUM / daily))
        pp = base_params()
        rr = reg_steroid(daily, days)
        inc = simulate_incident(pp, rr, 1826.0, phi=0.0,
                                t_stage1=max(540.0, days + 200.0))
        o1 = derived(inc["t1"], inc["Y1"], inc["p1"], inc["geo1"])
        tl = inc["alpha"]
        jt = jic_type(tl) if inc["alpha"] > 1 else "-"
        tcs = ("no lesion" if inc["alpha"] < 1 else
               ("never" if math.isinf(inc["tcol"]) else f"{inc['tcol']/30.44:.1f}"))
        w(f"{daily:>9}{days:>6}{o1['GC_eff'].max():>11.4f}"
          f"{o1['ADIPO'].max():>10.3f}{o1['PIO'].max():>8.1f}"
          f"{o1['Qrel'].min():>9.3f}{inc['nf']:>8.4f}{inc['alpha']:>7.1f}"
          f"{4*inc['alpha']:>6.0f}{jt:>5}{tcs:>13}")
    w("")
    w("  Epidemiology: peak daily dose (above about 40 mg/d prednisolone) and")
    w("  the dose given in the first 2-3 months dominate; low-dose maintenance")
    w("  of the SAME total is far safer.  The model was given one Hill")
    w("  coefficient of 2.5 on the adipogenic drive, with an EC50 placed at")
    w("  the exposure a 40 mg/d schedule produces.  The dominance of peak over")
    w("  cumulative is then arithmetic, not an assumption.")
    w("")

    # ---------------------------------------------------------------- 10
    w("[10] PROPHYLAXIS DURING STEROID EXPOSURE")
    w("-" * 79)
    w("  Statin, anticoagulant and prostacyclin all act on CLOCK 1.  The model")
    w("  therefore predicts they can only work BEFORE or DURING the insult.")
    w(f"{'arm':<40}{'PIOmax':>8}{'Qrelmin':>9}{'NECF':>8}"
      f"{'alpha':>7}{'JIC':>5}{'collapse(mo)':>13}")
    arms = [
        ("60 mg/d x 30 d then 90 d taper", lambda: reg_steroid(
            60.0, 90, start=30.0, taper_to=0.10, reg=reg_steroid(60.0, 30))),
        ("  + rosuvastatin from day 0", lambda: reg_statin(0.0, 365.0,
            reg=reg_steroid(60.0, 90, start=30.0, taper_to=0.10,
                            reg=reg_steroid(60.0, 30)))),
        ("  + enoxaparin x 12 weeks", lambda: reg_enoxaparin(0.0, 84.0,
            reg=reg_steroid(60.0, 90, start=30.0, taper_to=0.10,
                            reg=reg_steroid(60.0, 30)))),
        ("  + iloprost x 5 d at day 14", lambda: reg_iloprost(14.0, 5,
            reg=reg_steroid(60.0, 90, start=30.0, taper_to=0.10,
                            reg=reg_steroid(60.0, 30)))),
        ("same total dose at 12.4 mg/d x 375 d",
            lambda: reg_steroid(12.4, 375)),
        ("statin started 6 months AFTER the pulse", lambda: reg_statin(
            183.0, 365.0, reg=reg_steroid(60.0, 90, start=30.0, taper_to=0.10,
                                          reg=reg_steroid(60.0, 30)))),
    ]
    for lab, mk in arms:
        pp = base_params()
        if "thrombophil" in lab or "enoxaparin" in lab:
            pp["THROMBOPHILIA"] = 0.8
        inc = simulate_incident(pp, mk(), 1826.0, phi=0.0)
        o1 = derived(inc["t1"], inc["Y1"], inc["p1"], inc["geo1"])
        jt = jic_type(inc["alpha"]) if inc["alpha"] > 1 else "-"
        tcs = ("no lesion" if inc["alpha"] < 1 else
               ("never" if math.isinf(inc["tcol"]) else
                f"{inc['tcol']/30.44:.1f}"))
        w(f"{lab:<40}{o1['PIO'].max():>8.1f}{o1['Qrel'].min():>9.3f}"
          f"{inc['nf']:>8.4f}{inc['alpha']:>7.1f}{jt:>5}{tcs:>13}")
    w("")

    # ---------------------------------------------------------------- 11
    w("[11] TREATMENT COMPARISON ON A MATCHED JIC C1 HIP")
    w("-" * 79)
    w(f"{'arm':<30}{'collapse(mo)':>13}{'depr5y':>9}{'HHS5y':>8}"
      f"{'P(THA)5y':>10}{'beta_min':>10}")
    for nm in ("S06_jic_C1", "S08_C1_alendronate", "S09_C1_core_decompression",
               "S10_C1_CD_BMAC", "S11_C1_teriparatide", "S12_C1_ALN_then_TPTD",
               "S17_C1_iloprost", "S18_C1_denosumab", "S19_C1_graft",
               "S20_C1_PWB"):
        rr = res[nm]
        tc = rr["tcol"]
        tcs = "never" if math.isinf(tc) else f"{tc/30.44:.1f}"
        w(f"{nm:<30}{tcs:>13}{rr['obs']['DEPR'][-1]:>9.2f}"
          f"{rr['obs']['HHS'][-1]:>8.1f}{rr['obs']['P_THA'][-1]:>10.3f}"
          f"{rr['obs']['beta_eff'].min():>10.3f}")
    w("")

    # ---------------------------------------------------------------- 12
    w("[12] LOCAL SENSITIVITY  (+-20%, JIC C1)")
    w("-" * 79)
    w("  Output is months to collapse, censored at 60.  Five-year depression is")
    w("  useless as a sensitivity output because it saturates at DEPR_max for")
    w("  every hip that collapses at all.")

    def months_to_collapse(pp):
        tt, YY, gg = simulate(pp, reg_none(), 1826.0, nout=200)
        tc = collapse_time(tt, YY, pp)
        return 60.0 if math.isinf(tc) else tc / 30.44

    base = lesion_params(*jic_preset("C1"))
    ref = months_to_collapse(base)
    rows = []
    for key in ("k_dmg", "m_fat", "RESP_max", "v_front", "k_fill", "k_minz",
                "S_tr0", "E_tr0", "Ncyc", "f_hip", "BW", "k_resp",
                "theta_c", "alpha_deg", "zeta_max", "f_spread", "L_eng",
                "W_front"):
        vals = []
        for fac in (0.8, 1.2):
            pp = lesion_params(*jic_preset("C1"))
            pp[key] = pp[key] * fac
            vals.append(months_to_collapse(pp))
        sens = (vals[1] - vals[0]) / (0.4 * max(ref, 1e-6))
        rows.append((key, vals[0], ref, vals[1], sens))
    rows.sort(key=lambda z: -abs(z[4]))
    w(f"{'parameter':<12}{'-20%':>9}{'base':>9}{'+20%':>9}{'norm. sens.':>13}")
    for k, lo, b, hi, sn in rows:
        w(f"{k:<12}{lo:>9.2f}{b:>9.2f}{hi:>9.2f}{sn:>13.3f}")
    w("")

    # ---------------------------------------------------------------- 13
    w("[13] THE HAZARD IS NOT MONOTONE  (a prediction of clock 2, not a fit)")
    w("-" * 79)
    w("  Because the interface loses strength only while the front is crossing")
    w("  and regains it afterwards, the collapse hazard has to rise, peak and")
    w("  then vanish.  Conditional probability of collapsing in the NEXT year,")
    w("  given the hip is still intact now:")
    tcs = np.array([x["tcol"] for x in cohort])
    for t0 in (0, 365, 730, 1096, 1461):
        alive = tcs > t0
        if alive.sum() == 0:
            continue
        nxt = ((tcs > t0) & (tcs <= t0 + 365)).sum()
        w(f"     intact at {t0/30.44:5.0f} months (n={int(alive.sum()):4d})"
          f"  ->  collapse within the next 12 months: {nxt/alive.sum():.3f}")
    w("  Clinically: the great majority of collapses occur inside two years,")
    w("  and a head still spherical at four years is very unlikely to fail.")
    w("  Nothing in the model was fitted to that; it follows from the fact")
    w("  that creeping substitution finishes.")
    w("")

    w("=" * 79)
    w("END")
    w("=" * 79)

    txt = "\n".join(L)
    with open("onfh_reference_output.txt", "w") as f:
        f.write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
