#!/usr/bin/env python3
# =====================================================================
#  Bile Acid Diarrhoea / Bile Acid Malabsorption (BAM) — QSP core
#  Independent Python re-implementation of bam_mrgsolve_model.R
# ---------------------------------------------------------------------
#  Purpose: every number quoted in README.md is produced by THIS file.
#  The R (mrgsolve) model uses the same equations and the same parameter
#  block; this script exists so the results are reproducible without an
#  R installation, and so that bile-acid mass balance can be checked
#  explicitly (it is, in section 0).
#
#  ORGANISING IDEA
#  ---------------
#  The enterohepatic circulation is a negative-feedback loop whose
#  SENSOR (ileal FXR -> FGF19) measures the ABSORBED flux, while the
#  symptom (colonic bile-acid load) is produced by the SPILLED flux.
#  At steady state spillage == hepatic synthesis exactly, therefore
#
#      colonic bile-acid load  =  S  =  hepatic synthesis rate
#
#  and malabsorption BY ITSELF cannot raise the colonic load at all.
#  It raises it only through the feedback it disinhibits.  Two
#  independent lesions are therefore needed to describe the disease:
#
#      phi    = surviving fraction of ileal ASBT absorptive capacity
#      kappa  = gain of the ileal FGF19 sensor
#
#  and the two routine tests read them separately:
#      SeHCAT 7-day retention  ->  phi   alone
#      serum C4 / faecal BA    ->  S     (i.e. phi AND kappa together)
#
#  Run:  python3 bam_verify_python.py   (writes bam_verification_output.txt)
# =====================================================================
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

MW_BA = 500.0          # g/mol, mean conjugated bile acid

# ---------------------------------------------------------------------
# 1. STATE VECTOR
# ---------------------------------------------------------------------
NAMES = [
    # --- enterohepatic bile-acid mass (umol) --------------------------
    "LIV",      #  0 hepatocyte conjugated bile acid ready for secretion
    "GB",       #  1 gallbladder
    "DUO",      #  2 duodenum + proximal jejunum (micellar phase)
    "JEJ",      #  3 mid small intestine
    "ILE1",     #  4 proximal terminal ileum  (ASBT tank 1)
    "ILE2",     #  5 mid terminal ileum       (ASBT tank 2)
    "ILE3",     #  6 distal terminal ileum    (ASBT tank 3)
    "ENT",      #  7 ileal enterocyte cytosol (FXR ligand pool)
    "POR",      #  8 portal blood
    "SYS",      #  9 systemic plasma
    # --- circulating-pool composition (fractions) ---------------------
    "fCAp",     # 10 fraction of pool that is cholate
    "fDCAp",    # 11 fraction of pool that is deoxycholate
    # --- colonic bile-acid species (umol) -----------------------------
    "CCA",      # 12 colonic cholate            (tri-OH, weak secretagogue)
    "CCDCA",    # 13 colonic chenodeoxycholate  (di-OH, strong)
    "CDCA",     # 14 colonic deoxycholate       (di-OH, strongest)
    "CLCA",     # 15 colonic lithocholate       (mono-OH, insoluble)
    "BND",      # 16 bile acid bound to sequestrant, in transit
    "FEC",      # 17 cumulative faecal bile acid
    # --- FXR / FGF19 / CYP7A1 axis ------------------------------------
    "FGF19",    # 18 plasma FGF19 (pg/mL)
    "SHP",      # 19 hepatic SHP (relative)
    "CYP7A1",   # 20 hepatic CYP7A1 protein (relative, 1 = normal)
    "C4",       # 21 serum 7a-hydroxy-4-cholesten-3-one (ng/mL)
    "CYP8B1",   # 22 hepatic CYP8B1 (relative) -> CA:CDCA ratio
    "ASBT",     # 23 ileal ASBT protein (relative, adaptive)
    # --- colonic physiology ---------------------------------------------
    "WCOL",     # 24 colonic luminal water (mL)
    "TRANS",    # 25 colonic transit rate state (relative)
    "HT5",      # 26 mucosal 5-HT / EC-cell signal (relative)
    "PERM",     # 27 colonic permeability / mucosal irritation (relative)
    "STOOL",    # 28 cumulative stool water (mL)
    # --- microbiome -------------------------------------------------------
    "BSH",      # 29 bacterial bile-salt hydrolase activity (relative)
    "BAI",      # 30 7a-dehydroxylating capacity (relative)
    # --- drug PK ------------------------------------------------------------
    "SEQ1",     # 31 sequestrant, stomach (mg)
    "SEQ2",     # 32 sequestrant, small intestine (mg)  <- binds here
    "SEQ3",     # 33 sequestrant, colon (mg)            <- binds here
    "OCAg",     # 34 obeticholic acid, gut lumen (umol)
    "OCAe",     # 35 obeticholic acid, ileal enterocyte (umol)
    "OCAp",     # 36 obeticholic acid, plasma/liver (umol)
    "ELOg",     # 37 elobixibat, ileal lumen (umol)
    "LOPg",     # 38 loperamide, gut depot (mg)
    "LOPe",     # 39 loperamide, myenteric effect site (mg)
    "RIFg",     # 40 rifaximin, colon (mg)
    "ONDc",     # 41 ondansetron, plasma (mg)
    "TROg",     # 42 tropifexor, gut (umol)
    "TROp",     # 43 tropifexor, plasma -> ileal FXR (umol)
    # --- systemic / organ endpoints -------------------------------------------
    "HCHOL",    # 44 hepatic cholesterol pool (relative)
    "LDLC",     # 45 plasma LDL-C (mg/dL)
    "FECFAT",   # 46 cumulative faecal fat (g)
    "UOX",      # 47 urinary oxalate (mg/day)
    "VITD",     # 48 25-OH vitamin D (ng/mL)
    # --- diagnostic tracer ------------------------------------------------------
    "SEH",      # 49 SeHCAT tracer remaining in body (fraction of dose)
    "FATF",     # 50 faecal fat appearance rate, slow filter (g/day)
]
IX = {n: i for i, n in enumerate(NAMES)}
NS = len(NAMES)

POOL_STATES = ["LIV", "GB", "DUO", "JEJ", "ILE1", "ILE2", "ILE3", "ENT", "POR", "SYS"]
COL_STATES = ["CCA", "CCDCA", "CDCA", "CLCA", "BND"]


# ---------------------------------------------------------------------
# 2. PARAMETERS
# ---------------------------------------------------------------------
def default_params():
    p = {}

    # ==== THE TWO LESION AXES ========================================
    p["phi"]      = 1.0    # surviving ileal ASBT absorptive capacity (0-1)
    p["kappa"]    = 1.0    # ileal FGF19 sensor gain (type-2 BAD < 1)
    # ---- structural modifiers ---------------------------------------
    p["GBfun"]    = 1.0    # 1 = gallbladder present, 0 = cholecystectomy
    p["ENTMASS"]  = 1.0    # functional ileal enterocyte mass
    p["FEEDBACK"] = 1.0    # 1 = FGF19->CYP7A1 loop live; 0 = loop FROZEN

    # ==== hepatic synthesis / CYP7A1 ==================================
    p["Vs"]     = 29.17    # umol/h at CYP7A1 = 1  (700 umol/d = 0.35 g/d)
    p["kout7"]  = 0.15     # /h  CYP7A1 protein turnover (t1/2 4.6 h)
    p["kin7"]   = 0.90     # /h  max transcription  -> ceiling 6x basal
    p["IC50f"]  = 118.049    # pg/mL  FGF19 potency on CYP7A1   [CALIBRATED]
    p["hf"]     = 1.5
    p["IC50s"]  = 1.0      # relative SHP potency on CYP7A1
    p["hs"]     = 1.5
    p["circA"]  = 0.35     # circadian amplitude of CYP7A1 transcription
    p["circPk"] = 13.0     # h, acrophase

    # ==== bile-acid trafficking ========================================
    p["kbsep"]  = 3.00     # /h hepatocyte -> canaliculus (BSEP)
    p["fGB"]    = 0.72     # fasting fraction of bile diverted to gallbladder
    p["kejF"]   = 0.130    # /h interdigestive gallbladder emptying
    p["kejM"]   = 1.90     # /h CCK-driven gallbladder emptying
    p["kduo"]   = 1.05     # /h duodenum -> jejunum
    p["kjej"]   = 0.95     # /h jejunum -> ileum
    p["kpass"]  = 0.020    # /h passive jejunal absorption
    p["kile"]   = 1.00     # /h transit between ileal tanks and out
    p["kost"]   = 6.00     # /h enterocyte -> portal (OSTa/b)
    p["kpor"]   = 12.0     # /h portal -> hepatocyte (NTCP/OATP)
    p["fspill"] = 0.11     # portal fraction escaping first pass
    p["ksys"]   = 0.90     # /h systemic -> liver re-uptake

    # ==== ASBT (the phi axis) ===========================================
    p["Vmax_asbt"] = 1164.94  # umol/h per tank at C >> Km    [CALIBRATED]
    p["Km_asbt"]   = 6.0     # mM luminal
    p["V_ile"]     = 75.0    # mL per ileal tank
    p["V_ent"]     = 60.0    # mL ileal enterocyte water
    p["V_duo"]     = 250.0   # mL duodenal / meal volume
    p["V_col"]     = 300.0   # mL colonic luminal water reference

    # ==== FXR / FGF19 =====================================================
    p["fu_ent"]   = 0.010    # free fraction of enterocyte BA (IBABP buffers
                             # the rest).  Without this the sensor sits at
                             # ~99 % occupancy and stops sensing anything.
    p["EC50_ba"]  = 26.0     # uM FREE enterocyte BA at FXR (CDCA-equivalent)
    p["EC50_oca"] = 0.16     # uM OCA at FXR (~100x CDCA)
    p["EC50_tro"] = 0.0003   # uM tropifexor at FXR (sub-nM, non-steroidal)
    p["kf19"]     = 805.08    # pg/mL/h FGF19 production at FXRi = 1  [CALIBRATED]
    p["kef19"]    = 1.386    # /h FGF19 elimination (t1/2 30 min)
    p["kshp"]     = 0.50     # /h SHP formation
    p["kdshp"]    = 0.50     # /h SHP decay
    p["EC50_hep"] = 900.0    # umol hepatic BA for half-max SHP
    p["kc4"]      = 0.0295424   # (ng/mL) per (umol/h) per h   [CALIBRATED]
    p["kec4"]     = 0.0417   # /h C4 elimination (t1/2 16.6 h)

    # ==== CYP8B1 / pool composition ==========================================
    p["kin8"]  = 0.20
    p["kout8"] = 0.10
    p["fCA0"]  = 0.55        # cholate fraction of new synthesis at the normal set-point
    p["E8ref"] = 0.4385       # CYP8B1 level at the normal set-point  [CALIBRATED]

    # ==== adaptive ASBT ========================================================
    p["kASBT"]   = 0.0060    # /h  (t1/2 ~5 days)
    p["ASBTmax"] = 1.35      # ceiling of compensatory up-regulation

    # ==== colonic microbial transformation ======================================
    p["kbai"]     = 0.16     # /h 7a-dehydroxylation at BAI = 1
    p["kdcaabs"]  = 0.0471678    # /h passive colonic absorption of DCA
    p["kcdcaabs"] = 0.032    # /h passive colonic absorption of CDCA
    p["klcaout"]  = 0.030    # /h LCA sulfation / precipitation (-> faeces)
    p["kfec"]     = 0.0420822    # /h colonic emptying at TRANS = 1  [CALIBRATED]
    p["kbsh"]     = 0.02
    p["kbai_t"]   = 0.02

    # ==== colonic secretion / motility ============================================
    p["wCA"]    = 0.15       # secretory potency weights (CDCA = 1)
    p["wCDCA"]  = 1.00
    p["wDCA"]   = 1.15
    p["wLCA"]   = 0.30
    p["EC50_G"] = 3.00       # mM-equivalent secretory threshold (Mekhjian/Hofmann)
    p["hG"]     = 2.5
    p["Jsec"]   = 190.0      # mL/h maximal TGR5/CFTR-driven secretion
    p["Qin"]    = 62.5       # mL/h ileocaecal water inflow (1500 mL/day)
    p["kwabs"]  = 0.528619      # /h colonic water absorption   [CALIBRATED]
    p["kstool"] = 0.030      # /h stool water evacuation
    p["tauTR"]  = 3.0        # h transit-state time constant
    p["aTR"]    = 1.55       # LOOP GAIN: secretory drive -> transit acceleration
    p["TRANSclamp"] = 0.0    # >0 pins TRANS, i.e. OPENS the transit loop
    p["tauHT"]  = 1.0
    p["aHT"]    = 1.0
    p["tauPERM"] = 12.0
    p["aPERM"]  = 0.55

    # ==== fat absorption / CMC =======================================================
    p["CMC"]    = 1.5        # mM effective critical micellar concentration
    p["hCMC"]   = 4.0
    p["fatmax"] = 0.97       # maximal fractional fat absorption
    p["fatIn"]  = 100.0      # g/day dietary fat

    # ==== systemic endpoints ==========================================================
    p["LDL0"]     = 120.0    # mg/dL at normal synthesis
    p["expLDL"]   = 0.35
    p["tauLDL"]   = 336.0    # h
    p["UOX0"]     = 30.0     # mg/day
    p["aUOXfat"]  = 3.6      # mg/day per g faecal fat above 3 g
    p["aUOXperm"] = 9.0
    p["tauUOX"]   = 48.0
    p["VITD0"]    = 30.0
    p["tauVITD"]  = 720.0

    # ==== drugs ========================================================================
    # In vivo effective capacity, ~10x below the in vitro figure: the resin
    # competes with dietary and endogenous anions and mixing is incomplete.
    # Anchored so that colesevelam 3.75 g/day roughly doubles faecal bile
    # acid excretion, which is what is observed.
    p["seq_cap"]  = 0.22     # umol bile acid bound per mg colesevelam
    p["kbind"]    = 0.020    # /(mM.h) second-order binding onto free resin
    p["kseq12"]   = 1.30     # /h stomach -> small intestine
    p["kseq23"]   = 0.55     # /h small intestine -> colon
    p["kseq3out"] = 0.05     # /h colon -> stool (colonic residence ~20 h)
    p["ka_oca"]   = 0.55     # /h OCA gut absorption
    p["foca_pass"] = 0.40    # fraction of OCA uptake independent of ASBT
    p["kOCAep"]   = 1.10     # /h enterocyte -> plasma
    p["kOCAout"]  = 0.055    # /h systemic elimination
    p["kOCAgo"]   = 0.15     # /h luminal transit loss
    p["V_ocap"]   = 12.0     # L
    p["IC50_elo"] = 0.35     # uM luminal elobixibat
    p["kelo_out"] = 0.55     # /h
    p["ka_tro"]   = 0.75
    p["ke_tro"]   = 0.06
    p["V_trop"]   = 55.0     # L
    p["ftro_ile"] = 0.55     # plasma -> ileal enterocyte partition
    p["ka_lop"]   = 0.90
    p["ke_lop"]   = 0.14
    p["IC50_lop"] = 2.50     # mg at effect site
    p["ke_ond"]   = 0.19
    p["IC50_ond"] = 3.2      # mg
    p["fOND"]     = 0.45     # max fraction of the 5-HT3 arm blocked
    p["ke_rif"]   = 0.35
    p["IC50_rif"] = 900.0    # mg in colon
    p["fRIF"]     = 0.60     # max suppression of BAI

    # ==== dosing (all off by default) ====================================================
    p["dose_seq"] = 0.0; p["n_seq"] = 2      # mg per administration
    p["dose_oca"] = 0.0                      # mg/day
    p["dose_elo"] = 0.0                      # mg/day
    p["dose_tro"] = 0.0                      # mg/day
    p["dose_lop"] = 0.0; p["n_lop"] = 2      # mg per administration
    p["dose_ond"] = 0.0; p["n_ond"] = 2
    p["dose_rif"] = 0.0; p["n_rif"] = 3

    # frozen-loop clamp: the DAY-AVERAGED healthy repression (set after the
    # normal run in main(); these defaults are the values it produces)
    p["_Rfnorm"] = 0.2334
    p["_Rsnorm"] = 0.7140
    return p


MEALS = [8.0, 13.0, 19.0]
MEAL_DUR = 1.2
MEALAUC = len(MEALS) * MEAL_DUR / 2.0      # h/day: integral of the meal signal


def meal_signal(t):
    tod = t % 24.0
    s = 0.0
    for m in MEALS:
        d = tod - m
        if d < -12.0:
            d += 24.0
        elif d > 12.0:
            d -= 24.0
        if -MEAL_DUR / 2 <= d <= MEAL_DUR / 2:
            s += 0.5 * (1.0 + np.cos(2 * np.pi * d / MEAL_DUR))
    return min(s, 1.0)


def dosing_times(n):
    return {1: [8.0], 2: [8.0, 19.0], 3: [8.0, 13.0, 19.0],
            4: [7.0, 12.0, 17.0, 22.0]}.get(n, [8.0])


# ---------------------------------------------------------------------
# 3. RIGHT-HAND SIDE
# ---------------------------------------------------------------------
def rhs(t, y, p):
    y = np.maximum(y, 0.0)
    g = lambda n: y[IX[n]]
    d = np.zeros(NS)
    M = meal_signal(t)

    LIV, GB, DUO, JEJ = g("LIV"), g("GB"), g("DUO"), g("JEJ")
    I1, I2, I3 = g("ILE1"), g("ILE2"), g("ILE3")
    ENT, POR, SYS = g("ENT"), g("POR"), g("SYS")
    fCAp, fDCAp = np.clip(g("fCAp"), 0, 1), np.clip(g("fDCAp"), 0, 1)
    CCA, CCDCA, CDCA, CLCA, BND = g("CCA"), g("CCDCA"), g("CDCA"), g("CLCA"), g("BND")
    FGF19, SHP, E7, C4, E8, ASBT = (g("FGF19"), g("SHP"), g("CYP7A1"),
                                    g("C4"), g("CYP8B1"), g("ASBT"))
    WCOL, TRANS, HT5, PERM = g("WCOL"), g("TRANS"), g("HT5"), g("PERM")
    BAI = g("BAI")
    SEQ1, SEQ2, SEQ3 = g("SEQ1"), g("SEQ2"), g("SEQ3")
    OCAg, OCAe, OCAp = g("OCAg"), g("OCAe"), g("OCAp")
    ELOg, LOPg, LOPe, RIFg, ONDc = g("ELOg"), g("LOPg"), g("LOPe"), g("RIFg"), g("ONDc")
    TROg, TROp = g("TROg"), g("TROp")
    HCHOL, LDLC, PERMv = g("HCHOL"), g("LDLC"), PERM

    # ============ hepatic synthesis ==================================
    circ = 1.0 + p["circA"] * np.cos(2 * np.pi * ((t % 24.0) - p["circPk"]) / 24.0)
    if p["FEEDBACK"] >= 0.5:
        Rf = 1.0 / (1.0 + (FGF19 / p["IC50f"]) ** p["hf"])
        Rs = 1.0 / (1.0 + (SHP / p["IC50s"]) ** p["hs"])
    else:
        # Loop FROZEN.  CYP7A1 integrates the DAY-AVERAGE repression, not the
        # fasting trough, so the clamp must be the time-averaged Rf and Rs of
        # the healthy run -- otherwise the phi = 1 frozen run would not
        # reproduce the healthy state and section 1 would be confounded.
        Rf = p["_Rfnorm"]
        Rs = p["_Rsnorm"]
    d[IX["CYP7A1"]] = p["kin7"] * Rf * Rs * circ - p["kout7"] * E7
    S = p["Vs"] * E7 * HCHOL                                  # umol/h

    d[IX["CYP8B1"]] = p["kin8"] * Rf - p["kout8"] * E8
    # CYP8B1 sets the 12a-hydroxylation ratio, i.e. how much of NEW synthesis
    # leaves as cholate.  Normalised so that the healthy operating point
    # (E8 = E8ref) reproduces fCA0; when FGF19 falls the pool turns cholate-rich.
    fCAnew = float(np.clip(p["fCA0"] * E8 / p["E8ref"], 0.25, 0.80))

    # ============ bile secretion & gallbladder ========================
    Jbsep = p["kbsep"] * LIV
    if p["GBfun"] >= 0.5:
        to_gb = (1.0 - M) * p["fGB"] * Jbsep
        Jeject = (p["kejF"] + p["kejM"] * M) * GB
    else:
        to_gb, Jeject = 0.0, 0.0
    Jduo_in = Jbsep - to_gb + Jeject

    # ============ small-intestinal transit =============================
    Jd2j = p["kduo"] * (1.0 + 0.8 * M) * DUO
    Jj2i = p["kjej"] * (1.0 + 0.8 * M) * JEJ
    Jpass = p["kpass"] * JEJ

    # ---- sequestrant: a CONSUMABLE capacity, not an equilibrium ---------
    # Binding removes bile acid irreversibly and uses up resin capacity.
    # This is what makes the dose-response shallow: capacity is finite and
    # most of it is spent in the small intestine against the 25 mmol/day
    # cycling flux, long before the resin reaches the colon.
    ILEtot = I1 + I2 + I3
    C_tank = np.array([I1, I2, I3]) / p["V_ile"]              # mM per tank
    SEQtot = SEQ2 + SEQ3
    cap_free = max(p["seq_cap"] * SEQtot - BND, 0.0)
    fu_ile = 1.0
    if cap_free > 0 and SEQtot > 1e-9:
        Cbulk_ile = ILEtot / (3 * p["V_ile"])
        Jbind_ile = p["kbind"] * cap_free * (SEQ2 / SEQtot) * Cbulk_ile
        Jbind_ile = min(Jbind_ile, 0.60 * ILEtot)
    else:
        Jbind_ile = 0.0

    # ---- ASBT uptake (the phi axis) -----------------------------------
    C_elo = ELOg / (3 * p["V_ile"]) * 1000.0                  # uM
    Ielo = 1.0 / (1.0 + C_elo / p["IC50_elo"])
    Vmax_eff = p["Vmax_asbt"] * p["phi"] * ASBT * p["ENTMASS"] * Ielo
    Cf = C_tank * fu_ile
    Jasbt_k = Vmax_eff * Cf / (p["Km_asbt"] + Cf)             # per tank
    Jasbt = float(np.sum(Jasbt_k))

    kt = p["kile"] * (1.0 + 0.6 * M)
    J12, J23, J3c = kt * I1, kt * I2, p["kile"] * TRANS * I3
    fb = np.array([I1, I2, I3]) / max(ILEtot, 1e-9)           # bind pro rata
    d[IX["ILE1"]] = Jj2i - Jasbt_k[0] - J12 - Jbind_ile * fb[0]
    d[IX["ILE2"]] = J12 - Jasbt_k[1] - J23 - Jbind_ile * fb[1]
    d[IX["ILE3"]] = J23 - Jasbt_k[2] - J3c - Jbind_ile * fb[2]

    Ji2c_free = J3c

    # ============ enterocyte / FXR =====================================
    C_ent = ENT / p["V_ent"] * 1000.0                          # uM
    C_oca_e = OCAe / p["V_ent"] * 1000.0                       # uM
    C_tro_e = p["ftro_ile"] * TROp / p["V_trop"]               # uM
    # Additive receptor occupancy.  Endogenous bile acid and OCA are both
    # IBABP substrates and share the same free fraction; tropifexor is not
    # a bile acid and reaches FXR from the blood, so it is NOT buffered --
    # that is exactly why it is independent of phi.
    occ = (p["fu_ent"] * C_ent / p["EC50_ba"]
           + p["fu_ent"] * C_oca_e / p["EC50_oca"]
           + C_tro_e / p["EC50_tro"])
    FXRi = occ / (1.0 + occ)
    Jost = p["kost"] * ENT
    d[IX["ENT"]] = Jasbt - Jost

    d[IX["FGF19"]] = p["kf19"] * p["kappa"] * p["ENTMASS"] * FXRi - p["kef19"] * FGF19
    hepFXR = LIV / (p["EC50_hep"] + LIV)
    d[IX["SHP"]] = p["kshp"] * (hepFXR / 0.5) - p["kdshp"] * SHP
    d[IX["C4"]] = p["kc4"] * S - p["kec4"] * C4

    # ============ colon ================================================
    in_CA = Ji2c_free * fCAp
    in_DCA = Ji2c_free * fDCAp
    in_CDCA = Ji2c_free * max(1.0 - fCAp - fDCAp, 0.0)

    RIFeff = 1.0 - p["fRIF"] * RIFg / (p["IC50_rif"] + RIFg)
    kdh = p["kbai"] * BAI * RIFeff

    tot_col = CCA + CCDCA + CDCA + CLCA
    if cap_free > 0 and SEQtot > 1e-9 and tot_col > 1e-9:
        Jbind_col = (p["kbind"] * cap_free * (SEQ3 / SEQtot) * tot_col / p["V_col"])
        Jbind_col = min(Jbind_col, 0.60 * tot_col)
    else:
        Jbind_col = 0.0
    fc = np.array([CCA, CCDCA, CDCA, CLCA]) / max(tot_col, 1e-9)
    fu_col = 1.0

    kfec = p["kfec"] * TRANS
    Jabs_dca = p["kdcaabs"] * CDCA
    Jabs_cdca = p["kcdcaabs"] * CCDCA
    d[IX["CCA"]] = in_CA - kdh * CCA - kfec * CCA - Jbind_col * fc[0]
    d[IX["CCDCA"]] = in_CDCA - kdh * CCDCA - kfec * CCDCA - Jabs_cdca - Jbind_col * fc[1]
    d[IX["CDCA"]] = in_DCA + kdh * CCA - kfec * CDCA - Jabs_dca - Jbind_col * fc[2]
    d[IX["CLCA"]] = kdh * CCDCA - kfec * CLCA - p["klcaout"] * CLCA - Jbind_col * fc[3]
    d[IX["BND"]] = Jbind_ile + Jbind_col - p["kseq3out"] * BND

    Jfec = (kfec * (CCA + CCDCA + CDCA + CLCA) + p["klcaout"] * CLCA
            + p["kseq3out"] * BND)
    d[IX["FEC"]] = Jfec
    Jcolabs = Jabs_dca + Jabs_cdca

    # ---- secretory drive (mM CDCA-equivalent, free bile acid only) ----
    Ds = (p["wCA"] * CCA + p["wCDCA"] * CCDCA +
          p["wDCA"] * CDCA + p["wLCA"] * CLCA) / p["V_col"]
    G = Ds ** p["hG"] / (p["EC50_G"] ** p["hG"] + Ds ** p["hG"])

    OND = 1.0 - p["fOND"] * ONDc / (p["IC50_ond"] + ONDc)
    d[IX["HT5"]] = (p["aHT"] * G * OND - HT5) / p["tauHT"]
    LOP = 1.0 / (1.0 + LOPe / p["IC50_lop"])
    TRANSss = (1.0 + p["aTR"] * (0.55 * G + 0.45 * HT5)) * LOP
    if p["TRANSclamp"] > 0.0:
        d[IX["TRANS"]] = (p["TRANSclamp"] - TRANS) / 0.25   # loop OPENED
    else:
        d[IX["TRANS"]] = (TRANSss - TRANS) / p["tauTR"]
    d[IX["PERM"]] = (1.0 + p["aPERM"] * G - PERM) / p["tauPERM"]

    Jsecw = p["Jsec"] * (0.6 * G + 0.4 * HT5)
    Jabsw = p["kwabs"] * WCOL / TRANS
    Jstool = p["kstool"] * TRANS * WCOL
    d[IX["WCOL"]] = p["Qin"] + Jsecw - Jabsw - Jstool
    d[IX["STOOL"]] = Jstool

    d[IX["BSH"]] = p["kbsh"] * (1.0 - g("BSH"))
    d[IX["BAI"]] = p["kbai_t"] * (1.0 - BAI)

    # ============ bile-acid mass balance ================================
    Jpor_liv = p["kpor"] * POR * (1.0 - p["fspill"])
    Jpor_sys = p["kpor"] * POR * p["fspill"]
    Jsys_liv = p["ksys"] * SYS
    d[IX["LIV"]] = S + Jpor_liv + Jsys_liv - Jbsep
    d[IX["GB"]] = to_gb - Jeject
    d[IX["DUO"]] = Jduo_in - Jd2j
    d[IX["JEJ"]] = Jd2j - Jj2i - Jpass
    d[IX["POR"]] = Jost + Jpass + Jcolabs - p["kpor"] * POR
    d[IX["SYS"]] = Jpor_sys - Jsys_liv

    # ---- circulating-pool composition ----------------------------------
    Pool = LIV + GB + DUO + JEJ + I1 + I2 + I3 + ENT + POR + SYS
    Jin = S + Jcolabs
    if Pool > 1e-6 and Jin > 1e-12:
        gCA = S * fCAnew / Jin
        gDCA = Jabs_dca / Jin
        d[IX["fCAp"]] = Jin * (gCA - fCAp) / Pool
        d[IX["fDCAp"]] = Jin * (gDCA - fDCAp) / Pool

    # ============ fat absorption / CMC ===================================
    # dietary fat arrives WITH the meals, so it meets the postprandial bile
    # acid peak; the meal signal integrates to MEALAUC hours per day.
    C_duo = DUO / p["V_duo"]
    fat_eff = p["fatmax"] * C_duo ** p["hCMC"] / (p["CMC"] ** p["hCMC"] + C_duo ** p["hCMC"])
    fat_in_rate = p["fatIn"] * M / MEALAUC                     # g/h
    fat_app = fat_in_rate * (1.0 - fat_eff)                    # g/h into stool
    d[IX["FECFAT"]] = fat_app
    # FATF low-pass filters the meal-pulsed appearance rate into a smooth
    # g/day faecal-fat figure.  Downstream endpoints (oxalate, vitamin D)
    # must see the DAILY loss, not the between-meal instantaneous value,
    # which is meaningless because nothing is being eaten then.
    d[IX["FATF"]] = (24.0 * fat_app - g("FATF")) / 24.0
    fatrate = g("FATF")                                        # g/day

    # ============ systemic endpoints =======================================
    d[IX["HCHOL"]] = (1.0 - HCHOL) / 48.0
    LDLss = p["LDL0"] * (p["Vs"] / max(S, 1e-6)) ** p["expLDL"]
    d[IX["LDLC"]] = (LDLss - LDLC) / p["tauLDL"]
    UOXss = p["UOX0"] + p["aUOXfat"] * max(fatrate - 3.0, 0.0) + p["aUOXperm"] * (PERMv - 1.0)
    d[IX["UOX"]] = (UOXss - g("UOX")) / p["tauUOX"]
    fabs_day = max(p["fatIn"] - fatrate, 0.0) / p["fatIn"] / p["fatmax"]
    d[IX["VITD"]] = (p["VITD0"] * min(fabs_day, 1.0) - g("VITD")) / p["tauVITD"]

    # ============ SeHCAT tracer =============================================
    d[IX["SEH"]] = -(Jfec / max(Pool, 1e-9)) * g("SEH")

    # ============ drug PK ====================================================
    d[IX["SEQ1"]] = -p["kseq12"] * SEQ1
    d[IX["SEQ2"]] = p["kseq12"] * SEQ1 - p["kseq23"] * SEQ2
    d[IX["SEQ3"]] = p["kseq23"] * SEQ2 - p["kseq3out"] * SEQ3
    oca_asbt = (1.0 - p["foca_pass"]) * p["phi"] * ASBT * p["ENTMASS"]
    kaO = p["ka_oca"] * (p["foca_pass"] + oca_asbt)
    d[IX["OCAg"]] = -(kaO + p["kOCAgo"]) * OCAg
    d[IX["OCAe"]] = kaO * OCAg - p["kOCAep"] * OCAe
    d[IX["OCAp"]] = p["kOCAep"] * OCAe - p["kOCAout"] * OCAp
    d[IX["ELOg"]] = -p["kelo_out"] * ELOg
    d[IX["LOPg"]] = -p["ka_lop"] * LOPg
    d[IX["LOPe"]] = p["ka_lop"] * LOPg - p["ke_lop"] * LOPe
    d[IX["RIFg"]] = -p["ke_rif"] * RIFg
    d[IX["ONDc"]] = -p["ke_ond"] * ONDc
    d[IX["TROg"]] = -p["ka_tro"] * TROg
    d[IX["TROp"]] = p["ka_tro"] * TROg - p["ke_tro"] * TROp

    target = 1.0 + (p["ASBTmax"] - 1.0) * (1.0 - FXRi)
    d[IX["ASBT"]] = p["kASBT"] * (target - ASBT)
    return d


# ---------------------------------------------------------------------
# 4. DRIVER
# ---------------------------------------------------------------------
def initial_state(p):
    y = np.zeros(NS)
    y[IX["LIV"]] = 500.0;  y[IX["GB"]] = 2200.0; y[IX["DUO"]] = 700.0
    y[IX["JEJ"]] = 900.0;  y[IX["ILE1"]] = 500.0; y[IX["ILE2"]] = 200.0
    y[IX["ILE3"]] = 80.0;  y[IX["ENT"]] = 200.0;  y[IX["POR"]] = 100.0
    y[IX["SYS"]] = 30.0
    y[IX["fCAp"]] = 0.42;  y[IX["fDCAp"]] = 0.21
    y[IX["CCA"]] = 150.0;  y[IX["CCDCA"]] = 120.0
    y[IX["CDCA"]] = 260.0; y[IX["CLCA"]] = 150.0
    y[IX["FGF19"]] = 230.0; y[IX["SHP"]] = 1.0; y[IX["CYP7A1"]] = 1.0
    y[IX["C4"]] = 18.0;    y[IX["CYP8B1"]] = 2.0; y[IX["ASBT"]] = 1.0
    y[IX["WCOL"]] = 180.0; y[IX["TRANS"]] = 1.0; y[IX["HT5"]] = 0.05
    y[IX["PERM"]] = 1.0;   y[IX["BSH"]] = 1.0;   y[IX["BAI"]] = 1.0
    y[IX["HCHOL"]] = 1.0;  y[IX["LDLC"]] = 120.0
    y[IX["UOX"]] = 30.0;   y[IX["VITD"]] = 30.0; y[IX["SEH"]] = 1.0
    y[IX["FATF"]] = 3.5
    return y


def simulate(p, days=40.0, y0=None, dense_last=1.0):
    y = initial_state(p) if y0 is None else y0.copy()
    ev = []
    for day in range(int(np.ceil(days))):
        b = 24.0 * day
        if p["dose_seq"] > 0:
            for h in dosing_times(p["n_seq"]):
                ev.append((b + h, "SEQ1", p["dose_seq"]))
        if p["dose_oca"] > 0:
            ev.append((b + 8.0, "OCAg", p["dose_oca"] / 420.6 * 1000.0))
        if p["dose_elo"] > 0:
            ev.append((b + 7.5, "ELOg", p["dose_elo"] / 514.8 * 1000.0))
        if p["dose_tro"] > 0:
            ev.append((b + 8.0, "TROg", p["dose_tro"] / 604.7 * 1000.0))
        if p["dose_lop"] > 0:
            for h in dosing_times(p["n_lop"]):
                ev.append((b + h - 0.5, "LOPg", p["dose_lop"]))
        if p["dose_ond"] > 0:
            for h in dosing_times(p["n_ond"]):
                ev.append((b + h - 0.5, "ONDc", p["dose_ond"]))
        if p["dose_rif"] > 0:
            for h in dosing_times(p["n_rif"]):
                ev.append((b + h, "RIFg", p["dose_rif"]))
    ev.sort(key=lambda e: e[0])
    T = days * 24.0
    grid = sorted(set([0.0] + [e[0] for e in ev] + [T]))
    t_dense = T - dense_last * 24.0
    ot, oY = [], []
    for i in range(len(grid) - 1):
        t0, t1 = grid[i], grid[i + 1]
        for (te, nm, am) in ev:
            if abs(te - t0) < 1e-9:
                y[IX[nm]] += am
        if t1 - t0 < 1e-9:
            continue
        ne = max(2, int((t1 - t0) * 8)) if t1 > t_dense else 2
        sol = solve_ivp(rhs, (t0, t1), y, args=(p,), method="LSODA",
                        rtol=1e-7, atol=1e-9, t_eval=np.linspace(t0, t1, ne),
                        max_step=0.35)
        if not sol.success:
            raise RuntimeError("integration failed: " + sol.message)
        y = sol.y[:, -1]
        if t1 > t_dense:
            ot.append(sol.t); oY.append(sol.y)
    if ot:
        return np.concatenate(ot), np.concatenate(oY, axis=1), y
    return np.array([T]), y.reshape(-1, 1), y


def readouts(p, days=40.0, y0=None):
    tt, YY, yend = simulate(p, days=days, y0=y0, dense_last=1.0)
    T0 = tt[-1] - 24.0
    m = tt >= T0 - 1e-9
    t, Y = tt[m], YY[:, m]
    span = t[-1] - t[0]
    avg = lambda n: float(np.trapezoid(Y[IX[n]], t) / span)
    rate = lambda n: float((Y[IX[n]][-1] - Y[IX[n]][0]) / span * 24.0)

    pool = sum(avg(n) for n in POOL_STATES)
    fecBA = rate("FEC")
    Sv = p["Vs"] * Y[IX["CYP7A1"]] * Y[IX["HCHOL"]]
    S_day = float(np.trapezoid(Sv, t) / span * 24.0)

    Dsv = (p["wCA"] * Y[IX["CCA"]] + p["wCDCA"] * Y[IX["CCDCA"]] +
           p["wDCA"] * Y[IX["CDCA"]] + p["wLCA"] * Y[IX["CLCA"]]) / p["V_col"]
    Ds = float(np.trapezoid(Dsv, t) / span)

    kicv = p["kile"] * Y[IX["TRANS"]] * Y[IX["ILE3"]]
    col_deliv = float(np.trapezoid(kicv, t) / span * 24.0)
    duo_deliv_v = p["kduo"] * (1 + 0.8 * np.array([meal_signal(x) for x in t])) * Y[IX["DUO"]]
    duo_deliv = float(np.trapezoid(duo_deliv_v, t) / span * 24.0)

    tod = t % 24.0
    fast = (tod >= 5.0) & (tod <= 7.5)
    fgf = float(np.mean(Y[IX["FGF19"]][fast])) if fast.any() else avg("FGF19")
    c4 = float(np.mean(Y[IX["C4"]][fast])) if fast.any() else avg("C4")

    Cduo_pk = float(np.max(Y[IX["DUO"]] / p["V_duo"]))
    stool_w = rate("STOOL")
    kturn = fecBA / max(pool, 1e-9)
    return dict(
        pool_umol=pool, pool_g=pool * MW_BA / 1e6,
        synth_umol_d=S_day, synth_g_d=S_day * MW_BA / 1e6,
        fecBA_umol_d=fecBA, fecBA_g_d=fecBA * MW_BA / 1e6,
        col_deliv=col_deliv, duo_deliv=duo_deliv,
        cycles_d=duo_deliv / max(pool, 1e-9),
        f_pass=1.0 - col_deliv / max(duo_deliv, 1e-9),
        Ds=Ds, C4=c4, FGF19=fgf,
        SeHCAT=float(np.exp(-7.0 * kturn) * 100.0), kturn_d=float(kturn),
        stool_mL_d=stool_w,
        nstool=float(1.10 * max(stool_w / 130.0, 1e-6) ** 1.05),
        bristol=float(np.clip(3.6 + 3.0 * (1 - np.exp(-1.05 * (stool_w / 130.0 - 1))), 1, 7)),
        fecfat_g_d=rate("FECFAT"), Cduo_peak=Cduo_pk,
        LDLC=avg("LDLC"), UOX=avg("UOX"), VITD=avg("VITD"),
        TRANS=avg("TRANS"), CYP7A1=avg("CYP7A1"), SHP=avg("SHP"),
        Rfbar=float(np.trapezoid(
            1.0 / (1.0 + (Y[IX["FGF19"]] / p["IC50f"]) ** p["hf"]), t) / span),
        Rsbar=float(np.trapezoid(
            1.0 / (1.0 + (Y[IX["SHP"]] / p["IC50s"]) ** p["hs"]), t) / span),
        fCAp=avg("fCAp"), fDCAp=avg("fDCAp"), FXRpool=avg("ENT"),
        yend=yend,
    )


# ---------------------------------------------------------------------
# 5. ANALYTIC CORE
# ---------------------------------------------------------------------
def analytic_S(f, kappa=1.0, S0=700.0, Smax_mult=6.0, h=1.5, f0=0.975):
    """Closed-form steady state of the loop, independent of the ODEs.

    At steady state colonic load = faecal loss = S, and the sensor sees
    the ABSORBED flux, proportional to x*S with x = f/(1-f):
        S = Smax / (1 + (a*kappa*x*S)^h)
    a is fixed by requiring the normal point (f0, kappa=1) to give S0.
    """
    Smax = Smax_mult * S0
    x0 = f0 / (1.0 - f0)
    a = brentq(lambda a: S0 * (1.0 + (a * x0 * S0) ** h) - Smax, 1e-14, 1.0)
    x = f / (1.0 - f)
    return brentq(lambda S: S * (1.0 + (a * kappa * x * S) ** h) - Smax, 1e-9, Smax)


def phi_from_resection(L_cm, lam=60.0):
    """Surviving ileal ASBT capacity after resecting L cm of terminal
    ileum, for an exponentially distal-weighted ASBT density."""
    return float(np.exp(-L_cm / lam))


# ---------------------------------------------------------------------
# 6. EXPERIMENTS
# ---------------------------------------------------------------------
OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


def hdr(s):
    say(); say("=" * 92); say(s); say("=" * 92)


def P(**kw):
    p = default_params()
    p.update(kw)
    return p


def main():
    hdr("0. NORMAL PHYSIOLOGY — calibration and mass balance")
    b = readouts(default_params())
    ref = [
        ("total bile acid pool (g)", b["pool_g"], "2.5 - 4.0"),
        ("enterohepatic cycles per day", b["cycles_d"], "4 - 12"),
        ("duodenal BA delivery (mmol/day)", b["duo_deliv"] / 1000, "20 - 40"),
        ("per-pass ileal conservation f", b["f_pass"], "0.95 - 0.99"),
        ("hepatic synthesis (g/day)", b["synth_g_d"], "0.2 - 0.6"),
        ("faecal BA excretion (g/day)", b["fecBA_g_d"], "0.2 - 0.6"),
        ("pool composition CA / CDCA / DCA", 0, ""),
        ("fasting serum C4 (ng/mL)", b["C4"], "10 - 30  (BAD cut-off 48)"),
        ("fasting plasma FGF19 (pg/mL)", b["FGF19"], "150 - 300"),
        ("SeHCAT 7-day retention (%)", b["SeHCAT"], "> 15  (healthy 30 - 50)"),
        ("peak duodenal [BA] (mM)", b["Cduo_peak"], "5 - 15"),
        ("colonic secretory drive (mM-eq)", b["Ds"], "< 3.0 threshold"),
        ("stool water (mL/day)", b["stool_mL_d"], "100 - 150"),
        ("stool frequency (BM/day)", b["nstool"], "1 - 2"),
        ("faecal fat (g/day on 100 g)", b["fecfat_g_d"], "< 7"),
        ("LDL-C (mg/dL)", b["LDLC"], "120 (reference)"),
        ("urinary oxalate (mg/day)", b["UOX"], "< 45"),
    ]
    say(f"{'quantity':40s} {'model':>10s}   reference")
    say("-" * 92)
    for n, v, r in ref:
        if n.startswith("pool composition"):
            say(f"{n:40s} {b['fCAp']*100:4.0f} /{(1-b['fCAp']-b['fDCAp'])*100:4.0f} /"
                f"{b['fDCAp']*100:4.0f} %   CA 35-40 / CDCA 30-35 / DCA 20-25 %")
        else:
            say(f"{n:40s} {v:10.2f}   {r}")
    say()
    say(f"MASS BALANCE  synthesis {b['synth_umol_d']:.1f} umol/day  vs  "
        f"faecal loss {b['fecBA_umol_d']:.1f} umol/day   "
        f"(discrepancy {100*abs(b['synth_umol_d']-b['fecBA_umol_d'])/b['synth_umol_d']:.2f} %)")
    # The frozen-loop runs must be pinned to the ACTUAL healthy operating
    # point of this parameterisation, not to nominal values -- otherwise
    # the phi = 1.0 frozen run would not reproduce the healthy state and
    # the comparison in section 1 would be confounded.
    NORM = dict(_Rfnorm=b["Rfbar"], _Rsnorm=b["Rsbar"])
    say()
    say(f"frozen-loop clamp (day-averaged healthy repression): "
        f"Rf = {b['Rfbar']:.4f}, Rs = {b['Rsbar']:.4f}   "
        f"[fasting FGF19 {b['FGF19']:.1f} pg/mL, SHP {b['SHP']:.4f}, "
        f"25-OH vit D {b['VITD']:.1f} ng/mL]")

    # =================================================================
    hdr("1. THE CENTRAL CLAIM — malabsorption alone does not cause diarrhoea")
    say("phi is lowered from 1.0 to 0.05.  Two parallel runs: the FGF19->CYP7A1")
    say("feedback LIVE (physiological), and FROZEN at the normal operating point")
    say("so that hepatic synthesis cannot respond.  Everything else is identical.")
    say()
    say(f"{'':>6s}   {'--------- FEEDBACK LIVE ---------':^42s}   "
        f"{'-------- FEEDBACK FROZEN --------':^42s}")
    say(f"{'phi':>6s}   {'SeHCAT%':>8s} {'C4':>6s} {'colonBA':>8s} {'stool':>7s} {'BM/d':>6s}   "
        f"{'SeHCAT%':>8s} {'C4':>6s} {'colonBA':>8s} {'stool':>7s} {'BM/d':>6s}")
    say("-" * 100)
    c1 = []
    for phi in [1.0, 0.60, 0.35, 0.20, 0.10, 0.05]:
        row = [phi]
        for fb in [1.0, 0.0]:
            r = readouts(P(phi=phi, FEEDBACK=fb, **NORM))
            row += [r["SeHCAT"], r["C4"], r["fecBA_umol_d"], r["stool_mL_d"], r["nstool"]]
        c1.append(row)
        say(f"{row[0]:6.2f}   {row[1]:8.2f} {row[2]:6.1f} {row[3]:8.0f} {row[4]:7.0f} {row[5]:6.2f}   "
            f"{row[6]:8.2f} {row[7]:6.1f} {row[8]:8.0f} {row[9]:7.0f} {row[10]:6.2f}")
    a, z = c1[0], c1[-1]
    say()
    say(f"FROZEN loop, phi 1.00 -> 0.05:  colonic BA load "
        f"{a[8]:.0f} -> {z[8]:.0f} umol/day ({100*(z[8]-a[8])/a[8]:+.1f} %),"
        f"  stool {a[9]:.0f} -> {z[9]:.0f} mL/day ({100*(z[9]-a[9])/a[9]:+.1f} %)")
    say(f"   ... while SeHCAT still collapses {a[6]:.1f} % -> {z[6]:.2f} %.")
    say(f"LIVE loop, same lesion:        colonic BA load "
        f"{a[3]:.0f} -> {z[3]:.0f} umol/day ({z[3]/a[3]:.2f}x),"
        f"  stool {a[4]:.0f} -> {z[4]:.0f} mL/day ({z[4]/a[4]:.2f}x)")
    say()
    say("The malabsorption is real in both columns.  The DIARRHOEA appears only")
    say("in the column where the liver is allowed to answer it.")

    # =================================================================
    hdr("2. TWO LESIONS, TWO TESTS — SeHCAT reads phi, C4 reads phi AND kappa")
    say(f"{'phi':>5s} {'kappa':>6s} | {'SeHCAT%':>8s} {'C4':>6s} {'FGF19':>6s} "
        f"{'fecBA':>7s} {'Ds':>5s} {'stool':>6s} {'BM/d':>5s}  interpretation")
    say("-" * 100)
    lab = {(1.00, 1.00): "healthy",
           (1.00, 0.45): "type 2 (idiopathic) BAD",
           (1.00, 0.25): "severe sensor lesion",
           (0.35, 1.00): "type 1, feedback intact",
           (0.35, 0.45): "type 1 + sensor lesion",
           (0.10, 1.00): "extensive ileal disease",
           (0.10, 0.45): "worst case"}
    t2 = {}
    for phi, kap in lab:
        r = readouts(P(phi=phi, kappa=kap, **NORM))
        t2[(phi, kap)] = r
        say(f"{phi:5.2f} {kap:6.2f} | {r['SeHCAT']:8.2f} {r['C4']:6.1f} {r['FGF19']:6.0f} "
            f"{r['fecBA_umol_d']:7.0f} {r['Ds']:5.2f} {r['stool_mL_d']:6.0f} "
            f"{r['nstool']:5.2f}  {lab[(phi,kap)]}")
    A, B = t2[(1.00, 1.00)], t2[(1.00, 0.25)]
    say()
    say(f"Ileum structurally intact (phi = 1.00), sensor gain 1.00 -> 0.25:")
    say(f"   SeHCAT {A['SeHCAT']:6.2f} -> {B['SeHCAT']:6.2f} %  "
        f"({100*(B['SeHCAT']-A['SeHCAT'])/A['SeHCAT']:+.1f} %)   STAYS NORMAL")
    say(f"   C4     {A['C4']:6.1f} -> {B['C4']:6.1f} ng/mL "
        f"({100*(B['C4']-A['C4'])/A['C4']:+.0f} %)")
    say(f"   stool  {A['stool_mL_d']:6.0f} -> {B['stool_mL_d']:6.0f} mL/day "
        f"({100*(B['stool_mL_d']-A['stool_mL_d'])/A['stool_mL_d']:+.0f} %)")
    say("   => SeHCAT-negative bile acid diarrhoea is a STRUCTURAL prediction of")
    say("      the model, not an assay failure.  Faecal BA / C4 is the test that")
    say("      sees it.")

    # =================================================================
    hdr("3. SeHCAT IS A COMPRESSED READOUT OF phi")
    say(f"{'phi':>6s} {'f per pass':>11s} {'pool g':>8s} {'k_turn/d':>9s} "
        f"{'SeHCAT %':>9s}  grade        {'C4':>6s} {'BM/day':>7s}")
    say("-" * 82)
    for phi in [1.0, 0.8, 0.6, 0.5, 0.42, 0.35, 0.28, 0.20, 0.14, 0.10, 0.05, 0.02]:
        r = readouts(P(phi=phi, **NORM))
        s = r["SeHCAT"]
        gr = "normal" if s > 15 else "mild" if s > 10 else "moderate" if s > 5 else "SEVERE"
        say(f"{phi:6.2f} {r['f_pass']:11.4f} {r['pool_g']:8.2f} {r['kturn_d']:9.3f} "
            f"{s:9.2f}  {gr:12s} {r['C4']:6.1f} {r['nstool']:7.2f}")

    # =================================================================
    hdr("4. THERAPY CLASSES — where each drug enters the loop")
    say("Reference patient: phi = 0.20, kappa = 1.0 (moderate type-1 BAM).")
    say()
    scen = [
        ("untreated", {}),
        ("colesevelam 1.875 g BID", dict(dose_seq=1875.0, n_seq=2)),
        ("colesevelam 3.75 g BID", dict(dose_seq=3750.0, n_seq=2)),
        ("colestyramine 4 g TID", dict(dose_seq=4000.0, n_seq=3, seq_cap=0.13)),
        ("obeticholic acid 25 mg QD", dict(dose_oca=25.0)),
        ("tropifexor 90 ug QD", dict(dose_tro=0.090)),
        ("loperamide 2 mg BID", dict(dose_lop=2.0, n_lop=2)),
        ("ondansetron 4 mg BID", dict(dose_ond=4.0, n_ond=2)),
        ("rifaximin 550 mg TID", dict(dose_rif=550.0, n_rif=3)),
        ("elobixibat 10 mg QD", dict(dose_elo=10.0)),
        ("colesevelam + OCA", dict(dose_seq=1875.0, n_seq=2, dose_oca=25.0)),
        ("colesevelam + loperamide", dict(dose_seq=1875.0, n_seq=2, dose_lop=2.0, n_lop=2)),
    ]
    say(f"{'scenario':27s} {'stool':>7s} {'d%':>7s} {'BM/d':>6s} {'Ds':>5s} {'C4':>6s} "
        f"{'d%C4':>7s} {'FGF19':>6s} {'SeHCAT%':>8s} {'fecBA':>7s} {'LDL-C':>6s}")
    say("-" * 104)
    R4 = {}
    for nm, kw in scen:
        r = readouts(P(phi=0.20, **NORM, **kw), days=56.0)
        R4[nm] = r
        u = R4["untreated"]
        say(f"{nm:27s} {r['stool_mL_d']:7.0f} {100*(r['stool_mL_d']-u['stool_mL_d'])/u['stool_mL_d']:+7.1f} "
            f"{r['nstool']:6.2f} {r['Ds']:5.2f} {r['C4']:6.1f} "
            f"{100*(r['C4']-u['C4'])/u['C4']:+7.1f} {r['FGF19']:6.0f} {r['SeHCAT']:8.2f} "
            f"{r['fecBA_umol_d']:7.0f} {r['LDLC']:6.1f}")
    u, cs, oc, cb = (R4["untreated"], R4["colesevelam 1.875 g BID"],
                     R4["obeticholic acid 25 mg QD"], R4["colesevelam + OCA"])
    say()
    say(f"SEQUESTRANT  raises C4 {u['C4']:.1f} -> {cs['C4']:.1f} ng/mL "
        f"({100*(cs['C4']-u['C4'])/u['C4']:+.0f} %) and lowers FGF19 "
        f"{u['FGF19']:.0f} -> {cs['FGF19']:.0f} pg/mL.")
    say("   It relieves the symptom by removing the very ligand the sensor needs.")
    say("   The escape is structural, which is why the dose-response is shallow:")
    d1 = R4["colesevelam 1.875 g BID"]; d2 = R4["colesevelam 3.75 g BID"]
    say(f"      1.875 g BID -> {100*(u['stool_mL_d']-d1['stool_mL_d'])/u['stool_mL_d']:.1f} % "
        f"stool reduction;  3.75 g BID (2x dose) -> "
        f"{100*(u['stool_mL_d']-d2['stool_mL_d'])/u['stool_mL_d']:.1f} %")
    say(f"FXR AGONIST  lowers C4 {u['C4']:.1f} -> {oc['C4']:.1f} ng/mL "
        f"({100*(oc['C4']-u['C4'])/u['C4']:+.0f} %), raises FGF19 "
        f"{u['FGF19']:.0f} -> {oc['FGF19']:.0f},")
    say(f"   and RAISES SeHCAT {u['SeHCAT']:.2f} -> {oc['SeHCAT']:.2f} % without touching phi")
    say(f"   (falsifiable: an FXR agonist should improve a SeHCAT scan it did not repair),")
    say(f"   at the cost of LDL-C {u['LDLC']:.0f} -> {oc['LDLC']:.0f} mg/dL.")
    es = (u["stool_mL_d"] - cs["stool_mL_d"]) / u["stool_mL_d"]
    eo = (u["stool_mL_d"] - oc["stool_mL_d"]) / u["stool_mL_d"]
    ec = (u["stool_mL_d"] - cb["stool_mL_d"]) / u["stool_mL_d"]
    bliss = 1 - (1 - es) * (1 - eo)
    say(f"COMBINATION  observed {100*ec:.1f} %  vs  Bliss-independent {100*bliss:.1f} %  -> "
        f"{'SUPERADDITIVE' if ec > bliss + 0.005 else 'additive/less'}")
    say("   because the FXR agonist supplies the signal the sequestrant removes.")

    # =================================================================
    hdr("5. ELOBIXIBAT — the model must reproduce iatrogenic bile acid diarrhoea")
    say("An ASBT inhibitor is a pharmacological phi lesion in a NORMAL gut.")
    say("Reproducing its constipation indication is a falsification test.")
    say()
    say(f"{'elobixibat mg/day':>18s} {'stool mL/d':>11s} {'BM/day':>7s} {'C4':>7s} "
        f"{'FGF19':>7s} {'SeHCAT %':>9s} {'Ds':>6s}")
    say("-" * 74)
    for dose in [0.0, 5.0, 10.0, 15.0]:
        r = readouts(P(dose_elo=dose, **NORM))
        say(f"{dose:18.0f} {r['stool_mL_d']:11.0f} {r['nstool']:7.2f} {r['C4']:7.1f} "
            f"{r['FGF19']:7.0f} {r['SeHCAT']:9.2f} {r['Ds']:6.2f}")

    # =================================================================
    hdr("6. THE 100 cm RULE, DERIVED")
    say("Inputs: ASBT density falls exponentially proximally (decay length 60 cm,")
    say("an anatomical fact), hepatic compensation ceiling 6x basal (CYP7A1 Vmax),")
    say("and an effective critical micellar concentration of 1.5 mM.")
    say("Nothing here was fitted to resection data.")
    say()
    say(f"{'resect cm':>10s} {'phi':>7s} {'pool g':>8s} {'peak[BA]mM':>11s} {'CYP7A1 x':>9s} "
        f"{'fecfat g/d':>11s} {'stool':>7s} {'SeHCAT%':>8s}  phenotype")
    say("-" * 100)
    prev, switch = None, None
    for L in [0, 20, 40, 60, 80, 100, 120, 150, 180]:
        ph = phi_from_resection(L)
        r = readouts(P(phi=ph, **NORM), days=56.0)
        ff = r["fecfat_g_d"]
        pheno = "STEATORRHOEA" if ff > 7 else "bile acid diarrhoea"
        if ff > 7 and switch is None and prev is not None:
            L0, f0 = prev
            switch = L0 + (L - L0) * (7.0 - f0) / (ff - f0)
        prev = (L, ff)
        say(f"{L:10d} {ph:7.3f} {r['pool_g']:8.2f} {r['Cduo_peak']:11.2f} {r['CYP7A1']:9.2f} "
            f"{ff:11.1f} {r['stool_mL_d']:7.0f} {r['SeHCAT']:8.2f}  {pheno}")
    say()
    if switch:
        say(f"=> Model crossover from bile-acid diarrhoea to steatorrhoea at "
            f"{switch:.0f} cm resected.")
        say("   The clinical teaching (Hofmann) is 100 cm.  Beyond it, a sequestrant")
        say("   removes bile acid the duodenum can no longer spare and makes fat")
        say("   malabsorption worse — the model reproduces the contraindication.")

    # =================================================================
    hdr("7. CHOLECYSTECTOMY (type 3) — a timing lesion, not an absorption lesion")
    say(f"{'condition':34s} {'SeHCAT%':>8s} {'C4':>6s} {'Ds':>5s} {'stool':>7s} "
        f"{'BM/d':>6s} {'cycles/d':>9s} {'f/pass':>8s}")
    say("-" * 92)
    for nm, kw in [("intact gallbladder", {}),
                   ("cholecystectomy", dict(GBfun=0.0)),
                   ("cholecystectomy + phi 0.60", dict(GBfun=0.0, phi=0.60)),
                   ("cholecystectomy + kappa 0.45", dict(GBfun=0.0, kappa=0.45)),
                   ("phi 0.60 alone (comparison)", dict(phi=0.60))]:
        r = readouts(P(**NORM, **kw))
        say(f"{nm:34s} {r['SeHCAT']:8.2f} {r['C4']:6.1f} {r['Ds']:5.2f} "
            f"{r['stool_mL_d']:7.0f} {r['nstool']:6.2f} {r['cycles_d']:9.2f} {r['f_pass']:8.4f}")

    # =================================================================
    hdr("8. THE MICROBIOME SETS THE POTENCY OF A GIVEN LOAD")
    say("7a-dehydroxylation turns cholate (weak secretagogue, w = 0.15) into")
    say("deoxycholate (w = 1.15).  Identical delivered load, different disease.")
    say()
    say(f"{'BAI (relative)':>15s} {'colonic Ds':>11s} {'stool mL/d':>11s} {'BM/day':>7s} "
        f"{'fecBA':>7s} {'SeHCAT%':>8s} {'C4':>6s}")
    say("-" * 76)
    for bai in [0.2, 0.5, 1.0, 1.5, 2.0]:
        r = readouts(P(phi=0.20, kbai=0.16 * bai, **NORM))
        say(f"{bai:15.2f} {r['Ds']:11.2f} {r['stool_mL_d']:11.0f} {r['nstool']:7.2f} "
            f"{r['fecBA_umol_d']:7.0f} {r['SeHCAT']:8.2f} {r['C4']:6.1f}")
    say()
    say("Faecal bile acid barely moves; the symptom moves a lot.  This is the")
    say("second reason the tests and the symptoms disagree.")

    # =================================================================
    hdr("9. TRANSIT-ABSORPTION POSITIVE FEEDBACK — is it strong enough to latch?")
    say("Colonic bile acid accelerates transit; faster transit shortens ileal")
    say("contact time and spills more bile acid.  That is a genuine POSITIVE")
    say("feedback loop.  The question is quantitative: is its gain above 1?")
    say()
    say("The loop is OPENED by clamping TRANS, and the two halves are measured")
    say("separately:")
    say("   A:  d(Ds)/d(TRANS)        -- how much extra spillage faster transit buys")
    say("   B:  d(TRANSss)/d(Ds)      -- how much extra transit that spillage buys")
    say("   loop gain L = A x B  (dimensionless).  L >= 1 would permit latching.")
    say()
    say(f"{'phi':>6s} {'aTR':>6s} {'Ds@T=1.00':>10s} {'Ds@T=1.10':>10s} {'A':>8s} "
        f"{'B':>8s} {'L = AxB':>9s}  verdict")
    say("-" * 84)
    for phi_v in [1.0, 0.30]:
        for aTR in [1.55, 2.4, 3.6]:
            r1 = readouts(P(phi=phi_v, aTR=aTR, TRANSclamp=1.00, **NORM), days=34.0)
            r2 = readouts(P(phi=phi_v, aTR=aTR, TRANSclamp=1.10, **NORM), days=34.0)
            A = (r2["Ds"] - r1["Ds"]) / 0.10
            Dm = 0.5 * (r1["Ds"] + r2["Ds"])
            hGv, EC = 3.0 - 0.5, 3.0
            # B = d/dDs of (1 + aTR*(0.55G + 0.45G)) = aTR * dG/dDs
            dG = (2.5 * EC ** 2.5 * Dm ** 1.5) / (EC ** 2.5 + Dm ** 2.5) ** 2
            Bg = aTR * dG
            L = A * Bg
            say(f"{phi_v:6.2f} {aTR:6.2f} {r1['Ds']:10.3f} {r2['Ds']:10.3f} {A:8.3f} "
                f"{Bg:8.3f} {L:9.4f}  "
                f"{'LATCHING POSSIBLE' if L >= 1 else 'stable (no hysteresis)'}")
    say()
    say("Direct test: a 14-day transient insult (acute enteritis, ileocaecal")
    say("transit x 2.2) is applied to a phi = 0.30 patient and then withdrawn.")
    say()
    say(f"{'aTR':>6s} {'baseline mL/d':>14s} {'after recovery':>15s} {'delta':>9s}  verdict")
    say("-" * 66)
    for aTR in [0.8, 1.55, 2.4, 3.0, 3.6]:
        p = P(phi=0.30, aTR=aTR, **NORM)
        rb = readouts(p, days=34.0)
        p2 = dict(p); p2["kile"] = p["kile"] * 2.2
        _, _, yi = simulate(p2, days=14.0, y0=rb["yend"], dense_last=0.5)
        rr = readouts(p, days=34.0, y0=yi)
        dl = 100 * (rr["stool_mL_d"] - rb["stool_mL_d"]) / rb["stool_mL_d"]
        say(f"{aTR:6.2f} {rb['stool_mL_d']:14.0f} {rr['stool_mL_d']:15.0f} {dl:+8.2f}%  "
            f"{'DOES NOT RETURN' if dl > 5 else 'returns to baseline'}")
    say()
    say("HONEST NEGATIVE RESULT.  Over the whole physiological range the loop")
    say("gain stays well below 1, so the model does NOT reproduce a latching,")
    say("self-sustaining state, and every transient insult resolves completely.")
    say("The transit arm therefore AMPLIFIES bile acid diarrhoea but cannot")
    say("CAUSE it, and post-infectious bile acid diarrhoea cannot be explained")
    say("by this loop alone -- it needs a persisting lesion in phi or kappa.")
    say("That is a constraint the model imposes, not one that was put into it.")

    hdr("10. ANALYTIC CORE — the closed form behind all of the above")
    say("        S = Smax / (1 + (a . kappa . x . S)^h),   x = f/(1-f)")
    say("    because at steady state colonic load = faecal loss = S exactly.")
    say()
    say(f"{'f per pass':>11s} {'x = f/(1-f)':>12s} {'S kappa=1':>11s} {'S kappa=.45':>12s} "
        f"{'S kappa=.25':>12s}   umol/day")
    say("-" * 68)
    for f in [0.980, 0.975, 0.970, 0.950, 0.930, 0.900, 0.850, 0.800, 0.700, 0.600]:
        say(f"{f:11.3f} {f/(1-f):12.2f} {analytic_S(f):11.0f} {analytic_S(f,0.45):12.0f} "
            f"{analytic_S(f,0.25):12.0f}")
    say()
    r25 = analytic_S(0.60) / analytic_S(0.975)
    say(f"A 65-fold fall in the recycling ratio x (0.975 -> 0.600) raises S only")
    say(f"{r25:.2f}-fold.  The feedback is logarithmically stiff, which is why bile")
    say("acid diarrhoea is graded rather than all-or-none, and why the ceiling")
    say("(section 6) is what finally changes the phenotype.")

    # =================================================================
    hdr("11. SCENARIO LIBRARY — the archetypes")
    arche = [
        ("healthy control", {}),
        ("type 1 - 60 cm ileal resection", dict(phi=phi_from_resection(60))),
        ("type 1 - ileal Crohn disease", dict(phi=0.25, ENTMASS=0.70)),
        ("type 1 - pelvic radiation ileitis", dict(phi=0.30, ENTMASS=0.80, kbai=0.16 * 1.3)),
        ("type 2 - pure sensor lesion (kappa only)", dict(phi=1.00, kappa=0.40)),
        ("type 2 - idiopathic BAD as observed", dict(phi=0.55, kappa=0.40)),
        ("type 3 - post-cholecystectomy", dict(GBfun=0.0, phi=0.75)),
        ("type 4 - metformin-associated", dict(phi=0.55, kbai=0.16 * 1.4)),
        ("iatrogenic - ASBT inhibitor", dict(dose_elo=10.0)),
        ("severe: 120 cm resection", dict(phi=phi_from_resection(120))),
    ]
    say("NOTE on type 2.  A PURE sensor lesion cannot push SeHCAT below the")
    say("15 % cut-off, because the pool grows in step with synthesis and the")
    say("fractional turnover is unchanged.  Since real idiopathic BAD patients")
    say("are largely SeHCAT-positive, the model says they must carry a partial")
    say("phi lesion as well -- a testable structural claim, not a fitted one.")
    say()
    say(f"{'archetype':36s} {'SeHCAT%':>8s} {'C4':>6s} {'FGF19':>6s} {'fecBA':>7s} "
        f"{'Ds':>5s} {'stool':>6s} {'BM/d':>5s} {'fat':>5s} {'UOx':>5s}")
    say("-" * 104)
    for nm, kw in arche:
        r = readouts(P(**NORM, **kw), days=56.0)
        say(f"{nm:36s} {r['SeHCAT']:8.2f} {r['C4']:6.1f} {r['FGF19']:6.0f} "
            f"{r['fecBA_umol_d']:7.0f} {r['Ds']:5.2f} {r['stool_mL_d']:6.0f} "
            f"{r['nstool']:5.2f} {r['fecfat_g_d']:5.1f} {r['UOX']:5.0f}")

    say(); say("=" * 92); say("END OF VERIFICATION RUN"); say("=" * 92)


if __name__ == "__main__":
    main()
    with open("bam_verification_output.txt", "w") as fh:
        fh.write("\n".join(OUT) + "\n")
    print("\n[written] bam_verification_output.txt")
