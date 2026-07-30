## =============================================================================
##  hnscc_mrgsolve_model.R
##  Head and Neck Squamous Cell Carcinoma (HNSCC)
##  Quantitative Systems Pharmacology model for mrgsolve
##
##  72 ODE compartments.  Time unit = DAYS (so that a 25-minute free-platinum
##  half-life, a 2-minute radiation fraction and a 5-year survival curve live
##  in the same system).
##
##  ---------------------------------------------------------------------------
##  THE ONE STRUCTURAL COMMITMENT
##  ---------------------------------------------------------------------------
##
##  The benefit of a resistance-directed agent is NOT a property of the agent.
##  It is the product of what the agent does and how much of the resistance it
##  targets the tumour still has left to lose:
##
##      Delta(log kill) = (agent effect on factor i) x (HEADROOM in factor i)
##
##  This file contains no equation that says "cetuximab is inferior to
##  cisplatin in HPV-positive disease", no equation that says "a treatment gap
##  costs about 1 % of locoregional control per day", and no equation that says
##  "hypoxia modifiers only help hypoxic tumours".  It contains FOUR resistance
##  mechanisms that all act on the same log-kill product:
##
##      R_hypoxia       the oxygen enhancement ratio, applied only to the
##                      clonogens that are actually hypoxic, because CSCO and
##                      CSCH are separate pools that exchange at KOXTR -- which
##                      makes reoxygenation a race rather than a parameter
##
##      R_repair        alpha and the sublethal-damage repair rate MUREP, both
##                      moving with NHEJ capacity, HR capacity and platinum
##                      adduct load
##
##      R_repopulation  clonogen birth, which accelerates once the damage
##                      signal DAMS has accumulated for about three weeks
##
##      R_immune        CD8 kill suppressed by PD-1 engagement, Treg, MDSC and
##                      HLA class-I loss
##
##  and it routes each agent into the factor it actually reaches:
##
##      nimorazole                -> R_hypoxia
##      cisplatin                 -> R_repair (adducts raise alpha, block NHEJ)
##      cetuximab                 -> R_repair AND R_repopulation
##      pembrolizumab / nivolumab -> R_immune
##      shortening treatment time -> R_repopulation
##
##  EVERY HEADROOM IS MEASURED, NOT ASSUMED.  Seven counterfactual integrators
##  (LKACT / LKOXI / LKREP / LKNRP / LKIMM / LKIMX / LKCHM) accumulate the log
##  kill the clonogen pool actually received alongside the log kill it would
##  have received with one resistance factor switched off:
##
##      HDHYP  = 1 - LKACT/LKOXI     fraction of RT log-kill lost to hypoxia
##      HDREPR = 1 - LKACT/LKREP     fraction lost to intact repair
##      HDPOP  = LKNRP/LKACT         fraction handed back by repopulation
##      HDIMM  = 1 - LKIMM/LKIMX     fraction of immune kill lost to escape
##
##  ---------------------------------------------------------------------------
##  WHY FOUR FACTORS AND NOT THREE -- A CLAIM THAT DID NOT SURVIVE
##  ---------------------------------------------------------------------------
##
##  The first draft grouped repair and repopulation into one factor, so that
##  cisplatin and cetuximab would share a single headroom and their combination
##  would be sub-additive -- which is what RTOG 0522 found.  Writing the
##  equations showed that grouping is not supportable: cisplatin acts on alpha
##  and cetuximab acts mostly on clonogen birth, and those are separable.  The
##  model therefore OVER-PREDICTS the cisplatin + cetuximab combination
##  (+0.93 log10 where the trial saw no benefit).  Scenario 12b applies the
##  compliance the trial actually reported -- a 5-day interruption and only two
##  of three cisplatin cycles -- and that removes about half of the excess
##  (+0.45 log10).  The residual is a genuine failure and is reported as one in
##  README.md rather than repaired by re-grouping the factors.
##
##  ---------------------------------------------------------------------------
##  WHAT IS AN INPUT AND WHAT IS AN OUTPUT
##  ---------------------------------------------------------------------------
##
##  INPUTS.  The HPV-positive phenotype is imposed, not derived.  HPV = 1
##  changes exactly six things, each a measured property of HPV-positive
##  oropharyngeal cancer and none of them a treatment outcome:
##
##      HRCAP    0.45 vs 1.00   E7-driven RAD51 mislocalisation / TIP60 loss
##      VASCQ    4.00 vs 1.00   better perfusion, lower 15-gene hypoxia score
##      LAMS0    0.070 vs 0.100 longer T_pot, slower clonogen birth
##      EGAMP    1.00 vs 2.30   EGFR amplification is a HPV-negative event
##      E7BYP    0.55 vs 0.00   Rb bypass makes cycle entry EGFR-independent
##      AGVIR    2.60 vs 1.00   E6/E7 are foreign antigens; TIL-hot tumours
##
##  FITTED.  Six parameters were fitted by Nelder-Mead to eight published
##  locoregional-control values (Bonner 2006 both arms; MACH-NC/RTOG 0129;
##  RTOG 1016 both arms; DAHANCA 5-85 both arms; and an HPV-positive
##  RT-alone reference).  Final objective 0.0089, every anchor within 0.06:
##
##      ALPHA0 0.2484 /Gy   FCSC0 2.58e-5   PHIPT 1.005
##      KHRA   0.4798       KOXTR 0.1625/d  VASCQH 4.00 (at its search bound)
##
##  VASCQH reaching its bound is worth stating plainly: the anchor set wants
##  HPV-positive tumours to be as well perfused as the model permits, which
##  means one perfusion parameter is carrying more of the HPV contrast than it
##  should.
##
##  OUTPUTS.  Everything else, including every number in the validation table
##  in README.md.  The model is not told that treatment gaps cost about 1 % of
##  locoregional control per day, that repopulation starts in week three or
##  four, that hyperfractionation helps, that IMRT preserves saliva, that
##  checkpoint benefit rises with CPS, that the hypoxic fraction of a HNSCC is
##  10-25 %, or that T_pot is 4-5 days.
##
##  ---------------------------------------------------------------------------
##  RADIOTHERAPY IS NOT A DISCRETE MAP IN THIS FILE
##  ---------------------------------------------------------------------------
##
##  Applying survival fractions as instantaneous state jumps between fractions
##  is the usual shortcut and it cannot represent incomplete repair, so this
##  model uses the exact continuous form of the linear-quadratic model instead.
##  Each fraction is a BOLUS of d Gy into a delivery buffer DBUF which empties
##  at KDEL = 2000/d, so the instantaneous dose rate is
##
##      DRATE = KDEL * DBUF        (an exponential pulse of area d Gy and mean
##                                  duration 43 s -- a real linac)
##
##  and the sublethal-damage state obeys
##
##      dxdt_SLD = DRATE - MUREP * SLD
##
##  so the lethal-lesion yield rate
##
##      HZ = ALPHAE*DRATE + 2*BETAE*SLD*DRATE
##
##  integrates to exactly ALPHAE*D + BETAE*G*D^2, with G the Lea-Catcheside
##  protraction factor.  Because MUREP/KDEL = 0.012 the within-fraction G is
##  0.99, and because MUREP is finite the BETWEEN-fraction G is < 1 whenever
##  fractions are closer than about six hours -- which is how
##  hyperfractionation and its incomplete-repair penalty are both outputs.
##  Setting NOREP = 1 abolishes sublethal repair and is used as a numerical
##  check: it raises the log kill of a 70 Gy course from 4.65 to 9.40, which is
##  the analytic value of alpha*D + beta*D^2 with D = 70 taken in one piece.
##
##  ---------------------------------------------------------------------------
##  UNITS
##  ---------------------------------------------------------------------------
##    time                    days
##    cell numbers            cells (1 cm^3 = 1e9 cells)
##    radiation dose          Gy;  dose rate Gy/day
##    cisplatin, carboplatin  mg (amounts), mg/L (concentrations)
##    5-FU, docetaxel         mg, mg/L
##    cetuximab, pembrolizumab, nivolumab   mg, mg/L
##    platinum adducts        arbitrary units, normalised so that the peak
##                            after 100 mg/m^2 cisplatin is close to 1
##    signalling nodes        fold of the HPV-negative untreated reference
##    immune cell pools       arbitrary units, 1.0 = untreated reference
##    ANC                     10^9 cells/L
##    salivary flow           mL/min
##    magnesium               mmol/L
##    TSH                     mIU/L
##    body weight             kg
##
##  Author: QSP Disease Model Library (Claude Code Routine)
##  Licence: see repository LICENSE.  EDUCATIONAL / RESEARCH USE ONLY --
##  not validated for clinical or regulatory use.
## =============================================================================

library(mrgsolve)

hnscc_code <- r"---(
$PROB
# Head and Neck Squamous Cell Carcinoma QSP model (72 ODEs)

$PARAM @annotated
// ======================================================================
//  PATIENT / TUMOUR PHENOTYPE  (the only things a "phenotype" may set)
// ======================================================================
HPV     :  0    : HPV16/p16 status (0 = negative, 1 = positive)
PY      : 35    : Tobacco pack-years (modifies perfusion and antigenicity)
SMOKRT  :  0    : Continued smoking during radiotherapy (0/1)
CPS     : 10    : PD-L1 combined positive score (0-100), sets PD-L1 set-point
HLALOSS :  0.15 : Fraction of clones with HLA-I / B2M loss (immune escape floor)
ERCC1   :  1.0  : NER capacity for platinum adduct removal (1 = reference)
NRF2ON  :  0    : NFE2L2/KEAP1 pathway activation (0/1) -> GSH detoxification
PIK3M   :  0    : PIK3CA activating mutation (0/1)
PTENF   :  1.0  : PTEN function (1 = intact, 0.4 = loss)
METON   :  0    : HGF/MET bypass activity (0-1)
V0      : 30    : Baseline gross tumour volume (cm^3)
FCSC0   : 2.58e-5: Clonogenic-stem-cell fraction of viable cells at baseline
AGE     : 60    : Age (years)
BSA     :  1.80 : Body surface area (m^2)
BW0     : 70    : Baseline body weight (kg)
CRCL0   : 100   : Baseline creatinine clearance (mL/min)
HBG     : 13.5  : Haemoglobin (g/dL)
ECOG    :  1    : ECOG performance status (0-4)
STG     :  3    : AJCC-8 stage group (used only in the survival read-out)

// ======================================================================
//  RADIOBIOLOGY -- TUMOUR
// ======================================================================
ALPHA0  :  0.2484: Tumour alpha for an oxic, repair-competent clonogen (/Gy)
ABTUM   : 10.0   : Tumour alpha/beta ratio (Gy)
KHRA    :  0.4798: Alpha gain per unit loss of HR capacity
PHIPT   :  1.0050: Maximal fractional alpha gain from platinum adducts
KPHIPT  :  0.35  : Adduct level at half-maximal alpha gain
MUREP0  : 24.0   : Sublethal-damage repair rate at full NHEJ (/d; T1/2 = 0.69 h)
NOREP   :  0     : Diagnostic switch: 1 = abolish sublethal repair
KDEL    : 2000   : Radiation delivery-buffer emptying rate (/d)
OERMAX  :  2.80  : Oxygen enhancement ratio in anoxia
FCSCR   :  0.85  : Clonogen radioresistance factor on alpha (CSC vs bulk)
KEGNH   :  0.1998: Nuclear-EGFR stimulation of NHEJ capacity
KPTNH   :  0.80  : Platinum-adduct inhibition of NHEJ capacity
NIMOF   :  0     : Nimorazole / carbogen given (0/1)
NIMOE   :  0.42  : Fractional collapse of (OER-1) achieved by nimorazole

// ======================================================================
//  TUMOUR POPULATION KINETICS
// ======================================================================
LAMS0   :  0.100 : Clonogen (CSC) division rate (/d)
LAMP0   :  0.150 : Bulk proliferating division rate (/d; T_pot = 4.6 d)
RHO0    :  0.615 : Baseline symmetric self-renewal probability
RHOMX   :  0.880 : Self-renewal probability at full repopulation drive
KLOSS0  :  0.096 : Baseline cell-loss rate from the bulk pool (/d)
FREPL   :  1.50  : Division-rate boost at full repopulation drive (loss of Q)
FNEC    :  0.15  : Fraction of lost bulk cells that persist as necrotic mass
FLOSSR  :  0.55  : Fractional collapse of cell loss at full repopulation drive
KDAMIN  :  0.90  : Damage-signal production per unit normalised kill rate
KDAMOUT :  0.025 : Damage-signal decay (/d) -> sets the repopulation lag
KDAM50  :  5.64  : Damage signal at half-maximal repopulation drive
HDAM    :  3.0   : Hill coefficient of the repopulation switch (sets kick-off)
NCAPF   : 25.0   : Carrying capacity as a multiple of baseline cell number
KLOSSH  :  0.030 : Cell-loss rate from the hypoxic quiescent pool (/d)
KNEC    :  0.045 : Anoxic necrosis rate from the hypoxic pool (/d)
KOXTR   :  0.1625: Oxic <-> hypoxic interconversion rate (/d, fitted)
KDOOM   :  0.35  : Clearance of lethally damaged cells (/d)
KCLRN   :  0.060 : Clearance of necrotic debris (/d)
FCSCCH  :  0.45  : Clonogen chemosensitivity relative to bulk
FHYCH   :  0.25  : Hypoxic-pool chemosensitivity relative to oxic
FCSCIM  :  0.70  : Clonogen susceptibility to CD8 kill relative to bulk
FHYIM   :  0.40  : Hypoxic-pool susceptibility to CD8 kill relative to oxic

// ======================================================================
//  PERFUSION, HYPOXIA, ANGIOGENESIS
// ======================================================================
PO2MX   : 64.0   : Scaling constant for the pO2 read-out (mmHg)
PO2H    :  8.0   : pO2 at which half the cells are hypoxic (mmHg)
HHYP    :  2.0   : Hill coefficient of the hypoxia switch
KDEM    :  1.0   : Oxygen demand per unit normalised viable cell number
VASCQ0  :  1.0   : Intrinsic perfusion quality (HPV-negative reference)
VASCQH  :  4.00  : Perfusion-quality multiplier when HPV-positive (at bound)
KPYVQ   :  0.0040: Perfusion-quality loss per pack-year
SMKRTE  :  0.22  : Extra perfusion loss if smoking continues during RT
KHBV    :  0.045 : Perfusion gain per g/dL of haemoglobin above 10
KVGS    :  0.30  : VEGF-A production (/d)
KVGH    :  2.20  : VEGF-A induction by the hypoxic fraction
KVGD    :  0.35  : VEGF-A degradation (/d)
KVIN    :  0.28  : Vessel formation rate per unit VEGF (/d)
KVOUT   :  0.14  : Vessel regression rate (/d)
KVRT    :  0.030 : Endothelial loss per Gy of delivered dose (/Gy)

// ======================================================================
//  SIGNALLING
// ======================================================================
EGAMP0  :  2.30  : EGFR surface density multiplier when HPV-negative
LIGA    :  1.0   : Autocrine TGF-alpha / amphiregulin tone
KPEG    :  6.0   : p-EGFR relaxation rate (/d)
KERK    :  8.0   : p-ERK relaxation rate (/d)
KAKT    :  8.0   : p-AKT relaxation rate (/d)
KSTA    :  4.0   : p-STAT3 relaxation rate (/d)
KCC     :  1.2   : Cyclin D1-CDK4/6 relaxation rate (/d)
KPI3KM  :  0.85  : PI3K gain from a PIK3CA mutation
KMETE   :  0.60  : ERK/AKT gain from MET bypass
KIL6S   :  0.35  : STAT3 gain per unit IL-6
KERKC   :  0.45  : Cyclin D1 drive from p-ERK
KAKTC   :  0.35  : Cyclin D1 drive from p-AKT
E7BYP0  :  0.55  : Fraction of cycle entry that is E7-driven when HPV-positive
PBASE   :  0.55  : EGFR-independent floor on cycle entry
CCNDREF :  0.98  : Cyclin D1-CDK4/6 activity in the untreated HPV-negative case
CDK2AL  :  0.30  : Extra Rb release from CDKN2A deletion when HPV-negative

// ======================================================================
//  IMMUNE MODULE
// ======================================================================
ANTGC   :  0.25  : Constitutive antigen release (normalised)
AGVIR0  :  2.60  : Antigen multiplier from E6/E7 when HPV-positive
KTMB    :  0.020 : Antigen gain per pack-year (neoantigen surrogate)
FICD    :  1.60  : Immunogenic-cell-death boost while RT or platinum on board
KDCIN   :  1.10  : Dendritic-cell activation rate (/d)
KDCOUT  :  0.55  : Dendritic-cell decay (/d)
KTPRIM  :  0.90  : Node priming rate per unit activated DC (/d)
KTNOUT  :  0.30  : Primed-precursor loss (/d)
KTNEMI  :  0.45  : Emigration of primed CD8 to tumour (/d)
FHOME   :  0.60  : Fraction of emigrating CD8 reaching the tumour
KT8OUT  :  0.22  : Intratumoural effector CD8 loss (/d)
KEXH    :  0.55  : Exhaustion rate per unit PD-1 signal (/d)
KTEXO   :  0.10  : Exhausted-cell loss (/d)
KREINV  :  0.16  : Re-invigoration of exhausted CD8 per unit PD-1 blockade (/d)
KIFN    :  0.70  : IFN-gamma production per unit effector CD8 (/d)
KIFND   :  1.40  : IFN-gamma clearance (/d)
KPDLIN  :  1.00  : PD-L1 turnover in (/d)
KPDLO   :  1.00  : PD-L1 turnover out (/d)
KIFNP   :  0.55  : Adaptive PD-L1 induction per unit IFN-gamma
CPSREF  : 20.0   : CPS that corresponds to a PD-L1 set-point of 1.0
KTRIN   :  0.30  : Treg influx (/d)
KTROUT  :  0.30  : Treg loss (/d)
KIL6T   :  0.45  : Treg induction per unit IL-6
KLACTT  :  0.80  : Treg induction per unit hypoxic fraction (lactate proxy)
KMDIN   :  0.30  : MDSC influx (/d)
KMDOUT  :  0.30  : MDSC loss (/d)
KIL6M   :  0.60  : MDSC induction per unit IL-6
KST3M   :  0.40  : MDSC induction per unit p-STAT3
KKILL   :  0.770 : Maximal CD8-mediated kill rate (/d, capacity not realised)
KMT     :  1.00  : Effector CD8 level at half-maximal kill
KTSUP   :  3.43  : Suppression strength of (Treg + MDSC) on CD8 kill
KADCC   :  0.045 : Maximal cetuximab ADCC kill rate (/d)
FCGR3A  :  1.0   : FcgammaRIIIa activity multiplier (V/V 1.35, F/F 0.75)

// ======================================================================
//  CISPLATIN / CARBOPLATIN PK  (amounts in mg)
// ======================================================================
V1CIS   : 12.0   : Free-platinum central volume (L)
V2CIS   : 25.0   : Free-platinum peripheral volume (L)
QCIS    : 120    : Free-platinum intercompartmental clearance (L/d)
CLRCIS  : 288    : Free-platinum renal clearance at CrCl 100 (L/d)
CLBCIS  : 192    : Irreversible protein-binding clearance (L/d)
CLBEL   :  1.20  : Elimination of protein-bound platinum (L/d)
V1BND   : 14.0   : Volume for the protein-bound platinum pool (L)
KGSH    :  0.35  : Fractional platinum inactivation when NRF2 is on
KPCIS   :  0.55  : Tumour-to-plasma partition for free platinum
KADD    :  5.00  : Adduct formation rate per mg/L of tumour platinum (/d)
KNER    :  0.25  : Adduct removal rate at ERCC1 = 1 (/d; T1/2 = 2.8 d)
KUPK    :  0.55  : Renal cortical platinum uptake (/d)
KOUTK   :  0.11  : Renal cortical platinum efflux (/d)
KUPCO   :  0.075 : Cochlear platinum uptake (/d)
KOUTCO  :  0.020 : Cochlear platinum efflux (/d)
CARBOF  :  0     : Carboplatin used instead of cisplatin (0/1)
CARBEQ  :  0.24  : Adduct-forming potency of carboplatin relative to cisplatin
CARBOTO :  0.10  : Ototoxic/nephrotoxic potency of carboplatin relative to cis

// ======================================================================
//  5-FU, DOCETAXEL PK/PD
// ======================================================================
V1FU    : 20.0   : 5-FU central volume (L)
V2FU    : 18.0   : 5-FU peripheral volume (L)
QFU     : 300    : 5-FU intercompartmental clearance (L/d)
CLFU    : 4320   : 5-FU clearance (L/d)
DPDF    :  1.0   : DPD activity multiplier (0.5 = heterozygous DPYD variant)
V1DOC   :  8.0   : Docetaxel central volume (L)
V2DOC   : 80.0   : Docetaxel peripheral volume (L)
QDOC    : 720    : Docetaxel intercompartmental clearance (L/d)
CLDOC   : 500    : Docetaxel clearance (L/d)
EMXCIS  :  0.155 : Maximal platinum direct cytotoxic kill rate (/d)
EC50PT  :  0.55  : Adduct level at half-maximal platinum cytotoxicity
EMXFU   :  0.075 : Maximal 5-FU kill rate (/d)
EC50FU  :  0.32  : 5-FU concentration at half-maximal kill (mg/L)
EMXDOC  :  0.115 : Maximal docetaxel kill rate (/d)
EC50DOC :  0.55  : Docetaxel concentration at half-maximal kill (mg/L)

// ======================================================================
//  CETUXIMAB PK/PD  (amounts in mg)
// ======================================================================
V1CTX   :  3.0   : Cetuximab central volume (L)
V2CTX   :  2.6   : Cetuximab peripheral volume (L)
QCTX    :  0.55  : Cetuximab intercompartmental clearance (L/d)
CLCTX   :  0.52  : Cetuximab linear clearance (L/d)
VMCTX   : 22.0   : Cetuximab target-mediated elimination Vmax (mg/d)
KMCTX   : 12.0   : Cetuximab target-mediated elimination Km (mg/L)
KUPCT   :  0.55  : Tumour-interstitium equilibration rate for cetuximab (/d)
KPCTX   :  0.16  : Tumour-to-plasma partition for cetuximab
KDCTX   :  1.10  : Cetuximab concentration for 50 % EGFR occupancy (mg/L)

// ======================================================================
//  PEMBROLIZUMAB / NIVOLUMAB PK/PD  (amounts in mg)
// ======================================================================
V1PEM   :  3.5   : Anti-PD-1 central volume (L)
V2PEM   :  4.0   : Anti-PD-1 peripheral volume (L)
QPEM    :  0.60  : Anti-PD-1 intercompartmental clearance (L/d)
CLPEM   :  0.22  : Anti-PD-1 clearance (L/d)
KUPPM   :  0.50  : Tumour-interstitium equilibration rate for anti-PD-1 (/d)
KPPEM   :  0.20  : Tumour-to-plasma partition for anti-PD-1
KDPD1   :  0.35  : Anti-PD-1 concentration for 50 % receptor occupancy (mg/L)
CTLARX  :  0     : Anti-CTLA-4 co-administration flag (0/1, priming boost)
KCTLA   :  0.35  : Priming gain from CTLA-4 blockade

// ======================================================================
//  NORMAL TISSUE -- MUCOSA
// ======================================================================
FMUC    :  1.00  : Fraction of the prescribed dose received by oral mucosa
ALPHAM  :  0.062 : Mucosal basal-cell alpha (/Gy)
ABMUC   : 10.0   : Mucosal alpha/beta ratio (Gy)
MUREPM  : 24.0   : Mucosal sublethal-damage repair rate (/d)
KMUCP   :  0.42  : Mucosal basal-cell proliferation rate (/d)
KMUCS   :  0.004 : Mucosal basal-cell influx floor (/d)
KTAIN   :  0.55  : Basal -> transit-amplifying flux (/d)
KTAOUT  :  0.55  : Transit-amplifying maturation rate (/d)
FTAR    :  0.35  : Radiosensitivity of transit-amplifying cells vs basal
KEIN    :  0.50  : Epithelial-integrity formation rate (/d)
KEOUT   :  0.50  : Epithelial-integrity loss rate (/d)
KPTMUC  :  0.24  : Platinum kill rate on mucosal basal cells (/adduct unit /d)
KIIN    :  1.20  : Mucosal inflammation production per unit ulceration (/d)
KAMP    :  0.85  : Amplitude of the NF-kappaB self-amplification loop
KAMP50  :  0.60  : Inflammation level at half-maximal self-amplification
KIOUT   :  0.55  : Mucosal inflammation resolution (/d)
KINFMU  :  0.080 : Inflammation-driven extra basal-cell loss (/d)
MUCGRX  :  4.0   : Maximum reported mucositis grade
MUCGRP  :  0.80  : Exponent mapping ulcerated fraction to WHO grade

// ======================================================================
//  NORMAL TISSUE -- SALIVARY, SWALLOWING, BONE, THYROID
// ======================================================================
FPAR    :  0.35  : Fraction of the prescribed dose received by parotid (IMRT)
ALPHAP  :  0.027 : Parotid acinar alpha (/Gy)
ABPAR   :  3.0   : Parotid alpha/beta ratio (Gy)
MUREPP  :  8.0   : Late-tissue sublethal-damage repair rate (/d; T1/2 = 2.1 h)
KACREG  :  0.0015: Acinar regeneration rate (/d)
SALF0   :  1.50  : Baseline stimulated salivary flow (mL/min)
SALEXP  :  1.35  : Exponent from acinar survival to salivary flow
FCON    :  0.75  : Fraction of the prescribed dose to pharyngeal constrictors
KFIB    :  0.011 : Constrictor fibrosis formation per Gy per unit TGF-beta
KFIBR   :  0.0015: Constrictor fibrosis resolution (/d)
KTGFB   :  0.60  : TGF-beta tone per unit (Treg + MDSC + damage)
KTSHIN  :  0.0045: TSH drift per normalised thyroid dose (/d)
KTSHO   :  0.010 : TSH return rate (/d)
TSH0    :  1.80  : Baseline TSH (mIU/L)
KTHIRA  :  0.90  : TSH gain from immune-related thyroiditis

// ======================================================================
//  SYSTEMIC TOXICITY
// ======================================================================
KHAIR   :  0.0135: Outer-hair-cell loss rate per mg^HHAIR of cochlear platinum
HHAIR   :  1.35  : Hill exponent for cochlear hair-cell loss
HLMAX   : 55.0   : Maximal high-frequency hearing threshold shift (dB)
KKID    :  0.0042: Tubular injury rate per mg of cortical platinum (/d)
KKREP   :  0.014 : Tubular repair rate (/d)
KMGIN   :  1.00  : Magnesium homeostatic influx (/d)
KMGOUT  :  1.00  : Magnesium homeostatic efflux (/d)
MG0     :  0.85  : Baseline serum magnesium (mmol/L)
KMGCTX  :  0.40  : Fractional Mg influx loss at full EGFR occupancy
KMGCIS  :  0.35  : Fractional Mg influx loss at full tubular injury
KRIN    :  0.28  : Rash formation rate at full EGFR occupancy (/d)
KROUT   :  0.09  : Rash resolution rate (/d)
MTT     :  4.20  : Mean neutrophil maturation transit time (d)
GAM     :  0.16  : Feedback exponent of the Friberg neutrophil model
ANC0    :  4.50  : Baseline ANC (10^9/L)
SLPCIS  : 14.00  : ANC slope for free platinum (L/mg)
SLPDOC  :  0.115 : ANC slope for docetaxel (L/mg)
SLPFU   :  0.220 : ANC slope for 5-FU (L/mg)
SARCOF  :  1.0   : Sarcopenia multiplier on myelotoxic slopes
KIRIN   :  0.011 : Immune-related AE accrual at full PD-1 occupancy (/d)
KIROUT  :  0.020 : Immune-related AE resolution (/d)

// ======================================================================
//  HOST STATE
// ======================================================================
KIL6    :  0.85  : IL-6 production per unit normalised tumour burden (/d)
KMUCIL  :  0.75  : IL-6 production per unit mucosal inflammation (/d)
KIL6D   :  1.10  : IL-6 clearance (/d)
IL6B    :  0.52  : IL-6 level below which no cachexia signal is generated
KWL     :  0.0022: Maximal fractional weight-loss rate (/d)
KWL50   :  1.00  : Excess IL-6 at half-maximal weight loss
KWNPO   :  0.0035: Extra fractional weight loss when oral intake fails (/d)
KWREC   :  0.030 : Weight/nutrition recovery rate (/d)
NPOTHR  :  2.60  : Mucositis grade above which oral intake fails

// ======================================================================
//  CLINICAL READ-OUT MAPPING
// ======================================================================
CELLD   : 1.0e9  : Cells per cm^3 of viable tumour
KCTD    :  4.5e3 : HPV ctDNA release per normalised kill flux (copies/mL/d)
KCTDCL  :  6.4   : HPV ctDNA clearance (/d; T1/2 = 2.6 h)
OSA     :  0.145 : Survival mapping: hazard weight on locoregional failure
OSB     :  0.058 : Survival mapping: hazard weight on distant failure
OSSTG   :  0.115 : Survival mapping: hazard per stage group above 1
OSECOG  :  0.145 : Survival mapping: hazard per ECOG point
OSAGE   :  0.011 : Survival mapping: hazard per year above 60
OSCOMP  :  0.020 : Survival mapping: competing hazard per unit organ injury
DMET0   :  0.16  : Baseline 3-year distant-failure probability at stage 3
SIGLK   :  3.6222: Lumped log-scale between-patient heterogeneity (fitted)

$INIT @annotated
CISC  : 0     : Free platinum, central (mg)
CISP  : 0     : Free platinum, peripheral (mg)
CISB  : 0     : Protein-bound platinum (mg)
CISK  : 0     : Renal cortical platinum (mg)
CISCO : 0     : Cochlear platinum (mg)
FUC   : 0     : 5-FU, central (mg)
FUP   : 0     : 5-FU, peripheral (mg)
DOCC  : 0     : Docetaxel, central (mg)
DOCP  : 0     : Docetaxel, peripheral (mg)
CTXC  : 0     : Cetuximab, central (mg)
CTXP  : 0     : Cetuximab, peripheral (mg)
CTXT  : 0     : Cetuximab, tumour interstitium (mg/L)
PEMC  : 0     : Anti-PD-1, central (mg)
PEMP  : 0     : Anti-PD-1, peripheral (mg)
PEMT  : 0     : Anti-PD-1, tumour interstitium (mg/L)
DBUF  : 0     : Radiation delivery buffer (Gy)
SLD   : 0     : Tumour sublethal damage (Gy-equivalent)
SLDM  : 0     : Mucosal sublethal damage (Gy-equivalent)
SLDP  : 0     : Late-tissue sublethal damage (Gy-equivalent)
CUMD  : 0     : Cumulative prescribed dose (Gy)
CSCO  : 5.37e5: Clonogenic stem cells, OXIC (cells)
CSCH  : 1.07e5: Clonogenic stem cells, HYPOXIC / quiescent (cells)
TUMOX : 1.95e10: Proliferating oxic tumour cells (cells)
TUMHY : 3.9e9 : Quiescent hypoxic tumour cells (cells)
TUMDM : 0     : Lethally damaged tumour cells (cells)
TUMNC : 6.6e9 : Necrotic tumour mass (cell-equivalents)
VASC  : 0.50  : Functional vascular density (fraction)
VEGF  : 0.90  : VEGF-A (normalised)
PEGF  : 1.0   : p-EGFR (normalised)
PERK  : 1.0   : p-ERK1/2 (normalised)
PAKT  : 1.0   : p-AKT (normalised)
STA3  : 1.0   : p-STAT3 (normalised)
CCND  : 1.0   : Cyclin D1-CDK4/6 activity (normalised)
NHEJ  : 1.0   : NHEJ repair capacity (normalised)
PTDNA : 0     : Platinum-DNA adduct burden (normalised)
DCA   : 0.45  : Activated dendritic cells (normalised)
TNLN  : 0.60  : Primed CD8 precursors, draining node (normalised)
CD8T  : 0.80  : Intratumoural effector CD8 (normalised)
TEX   : 0.40  : Exhausted intratumoural CD8 (normalised)
TREG  : 1.00  : FOXP3+ Treg (normalised)
MDSC  : 1.00  : MDSC (normalised)
IFNG  : 0.40  : IFN-gamma (normalised)
PDL1  : 0.50  : PD-L1 expression (normalised)
MUCB  : 1.00  : Mucosal basal stem-cell pool (fraction of normal)
MUCTA : 1.00  : Mucosal transit-amplifying pool (fraction of normal)
MUCE  : 1.00  : Mucosal epithelial integrity (fraction of normal)
MUCI  : 0     : Mucosal inflammation (normalised)
ACIN  : 1.00  : Parotid serous acinar survival (fraction)
CONF  : 0     : Pharyngeal-constrictor fibrosis (normalised)
HAIR  : 1.00  : Cochlear outer-hair-cell survival (fraction)
KIDT  : 1.00  : Renal tubular integrity (fraction)
MGSER : 0.85  : Serum magnesium (mmol/L)
RASH  : 0     : Acneiform rash severity (latent, 0-4)
PROLN : 4.50  : Neutrophil proliferating precursors (10^9/L equivalent)
TRN1  : 4.50  : Neutrophil transit compartment 1
TRN2  : 4.50  : Neutrophil transit compartment 2
TRN3  : 4.50  : Neutrophil transit compartment 3
ANC   : 4.50  : Circulating absolute neutrophil count (10^9/L)
IL6   : 0.52  : Interleukin-6 (normalised)
BWT   : 70    : Body weight (kg)
TSHX  : 1.80  : Serum TSH (mIU/L)
IRAE  : 0     : Immune-related adverse-event burden (normalised)
CTD   : 0     : HPV ctDNA (TTMV copies/mL)
DAMS  : 0     : Repopulation damage signal (normalised)
CUMPT : 0     : Cumulative platinum eliminated (mg)
LKACT : 0     : ACTUAL cumulative RT log-kill on clonogens (nats)
LKOXI : 0     : COUNTERFACTUAL RT log-kill if fully oxic (nats)
LKREP : 0     : COUNTERFACTUAL RT log-kill if repair abrogated (nats)
LKNRP : 0     : Cumulative clonogen birth handed back (nats)
LKIMM : 0     : ACTUAL cumulative immune log-kill on clonogens (nats)
LKIMX : 0     : COUNTERFACTUAL immune log-kill if unsuppressed (nats)
LKCHM : 0     : Cumulative chemotherapy log-kill on clonogens (nats)

$GLOBAL
#define CSCT  (CSCO + CSCH)
#define NVIAB (CSCO + CSCH + TUMOX + TUMHY)
#define NTOT  (CSCO + CSCH + TUMOX + TUMHY + TUMDM + TUMNC)

// Baseline reference cell number, used only to normalise the oxygen-demand
// and IL-6 terms.  Fixed at the default initial condition on purpose: it is a
// SCALE, not a state, and letting it move with the tumour would make the
// hypoxia switch scale-free and therefore unable to reoxygenate.
#define NREF  3.5e10

double posf(double x) { return x > 0.0 ? x : 0.0; }
double hillf(double x, double k, double n) {
  double xp = pow(posf(x), n);
  return xp / (pow(k, n) + xp);
}

// The six phenotype quantities, plus renal clearance, are resolved once per
// record in $MAIN and then read by $ODE and $TABLE.  Nothing that varies
// WITHIN a step is stored here: every quantity $TABLE reports is recomputed
// from the state vector, so no output can be a stale derivative-evaluation
// leftover.
double gVQ, gRBFREE, gHRCAP, gLAMSP, gEGAMP, gE7BYP, gAGVIR, gCLR, gALPHA;

$MAIN
// ---------------------------------------------------------------- PHENOTYPE
// Six things, and only six things, are switched by HPV status.  Each is a
// measured property of the tumour, not a treatment outcome.
double HRCAP = HPV > 0.5 ? 0.45 : 1.00;   // homologous-recombination capacity
double VASCQ = HPV > 0.5 ? VASCQH : VASCQ0;
double LAMSP = HPV > 0.5 ? 0.070 : LAMS0;
double EGAMP = HPV > 0.5 ? 1.00 : EGAMP0;
double E7BYP = HPV > 0.5 ? E7BYP0 : 0.0;
double AGVIR = HPV > 0.5 ? AGVIR0 : 1.0;

// Perfusion quality is degraded by tobacco exposure and further degraded if
// the patient keeps smoking through the course.
double VQ = VASCQ * (1.0 - KPYVQ * (PY < 60 ? PY : 60))
                  * (1.0 - SMKRTE * SMOKRT)
                  * (1.0 + KHBV * (HBG - 10.0));
if (VQ < 0.15) VQ = 0.15;
gVQ = VQ;

// Rb release: E7 in HPV-positive disease, CDKN2A deletion in HPV-negative.
gRBFREE = HPV > 0.5 ? 1.0 : CDK2AL;
gHRCAP  = HRCAP;
gLAMSP  = LAMSP;
gEGAMP  = EGAMP;
gE7BYP  = E7BYP;

// Renal function scales free-platinum clearance; tubular injury feeds back.
gCLR = CLRCIS * (CRCL0 / 100.0);

// ---------------------------------------------------- between-patient effects
// Four random effects, and they are USED (see responder_rate() below), not
// declared for decoration.  They sit on the four quantities that dominate
// between-patient variability in outcome: how many clonogens have to be
// sterilised, intrinsic radiosensitivity, perfusion, and how immunologically
// visible the tumour is.  Every deterministic scenario in this file calls
// zero_re(), so the contrasts are unaffected by them.
double FCSC   = FCSC0  * exp(ETA(1));
double ALPHAI = ALPHA0 * exp(ETA(2));
double AGVI   = AGVIR  * exp(ETA(4));
gALPHA = ALPHAI;
gAGVIR = AGVI;
// gVQ was already assigned above, so the perfusion random effect has to be
// applied to gVQ and not to the local VQ -- writing `VQ = VQ * exp(ETA(3))`
// here was dead code and ETA(3) had no effect on anything.
gVQ = gVQ * exp(ETA(3));

// Initial conditions that must follow the covariates
// V0 is the TOTAL radiological volume; 78 % of it is viable at presentation.
double VIAB0 = V0 * CELLD * 0.78;
// The clonogen pool is split the same way the bulk is, because a hypoxic
// clonogen is the one cell in the tumour that both survives radiation and can
// regrow it.  Keeping them in one mean-field pool with an averaged OER makes
// hypoxia far too powerful a lever and abolishes reoxygenation as a mechanism.
double HF0 = HPV > 0.5 ? 0.036 : 0.167;
CSCO_0  = VIAB0 * FCSC * (1.0 - HF0);
CSCH_0  = VIAB0 * FCSC * HF0;
TUMOX_0 = VIAB0 * (1.0 - FCSC) * (1.0 - HF0);
TUMHY_0 = VIAB0 * (1.0 - FCSC) * HF0;
TUMNC_0 = V0 * CELLD * 0.22;
BWT_0   = BW0;
MGSER_0 = MG0;
TSHX_0  = TSH0;
PROLN_0 = ANC0;  TRN1_0 = ANC0;  TRN2_0 = ANC0;  TRN3_0 = ANC0;  ANC_0 = ANC0;

$ODE
// =====================================================================
//  1 · RADIATION DELIVERY AND THE CONTINUOUS LINEAR-QUADRATIC MODEL
// =====================================================================
double DRATE = KDEL * DBUF;              // Gy/day, an exponential pulse
dxdt_DBUF = -KDEL * DBUF;
dxdt_CUMD =  DRATE;

// =====================================================================
//  2 · PLATINUM PK, ADDUCTS AND TISSUE LOADING
// =====================================================================
double CFREE = CISC / V1CIS;                     // mg/L free platinum
double GSHF  = 1.0 - KGSH * NRF2ON;              // NRF2-driven detoxification
double CLREN = gCLR * (KIDT < 0.15 ? 0.15 : KIDT);

dxdt_CISC = -CLREN * CFREE - CLBCIS * CFREE
            - QCIS * (CFREE - CISP / V2CIS)
            - KUPK * CISC - KUPCO * CISC
            + KOUTK * CISK + KOUTCO * CISCO;
dxdt_CISP =  QCIS * (CFREE - CISP / V2CIS);
dxdt_CISB =  CLBCIS * CFREE - CLBEL * CISB / V1BND;
dxdt_CISK =  KUPK * CISC - KOUTK * CISK;
dxdt_CISCO=  KUPCO * CISC - KOUTCO * CISCO;
dxdt_CUMPT=  CLREN * CFREE + CLBEL * CISB / V1BND;

// Carboplatin, when used, enters the same compartments but forms adducts far
// less efficiently and is far less oto-/nephrotoxic per mg.
double PTPOT = CARBOF > 0.5 ? CARBEQ : 1.0;
double PTTOX = CARBOF > 0.5 ? CARBOTO : 1.0;
double CTUM  = KPCIS * CFREE * GSHF * PTPOT;
dxdt_PTDNA = KADD * CTUM - KNER * ERCC1 * PTDNA;

// =====================================================================
//  3 · 5-FU AND DOCETAXEL PK
// =====================================================================
double CFU  = FUC / V1FU;
double CDOC = DOCC / V1DOC;
dxdt_FUC  = -(CLFU * DPDF) * CFU - QFU * (CFU - FUP / V2FU);
dxdt_FUP  =  QFU * (CFU - FUP / V2FU);
dxdt_DOCC = -CLDOC * CDOC - QDOC * (CDOC - DOCP / V2DOC);
dxdt_DOCP =  QDOC * (CDOC - DOCP / V2DOC);

// =====================================================================
//  4 · MONOCLONAL ANTIBODY PK  (cetuximab has a saturable target sink)
// =====================================================================
double CPCT = CTXC / V1CTX;
dxdt_CTXC = -CLCTX * CPCT - QCTX * (CPCT - CTXP / V2CTX)
            - VMCTX * CPCT / (KMCTX + CPCT) * gEGAMP;
dxdt_CTXP =  QCTX * (CPCT - CTXP / V2CTX);
dxdt_CTXT =  KUPCT * (KPCTX * CPCT - CTXT);

double CPPM = PEMC / V1PEM;
dxdt_PEMC = -CLPEM * CPPM - QPEM * (CPPM - PEMP / V2PEM);
dxdt_PEMP =  QPEM * (CPPM - PEMP / V2PEM);
dxdt_PEMT =  KUPPM * (KPPEM * CPPM - PEMT);

double OCCEG = CTXT / (KDCTX + CTXT);            // EGFR occupancy in tumour
double OCCPD = PEMT / (KDPD1 + PEMT);            // PD-1 occupancy in tumour

// =====================================================================
//  5 · SIGNALLING
// =====================================================================
double PEGTG = gEGAMP * (1.0 - OCCEG) * LIGA / (1.0 + LIGA);
dxdt_PEGF = KPEG * (PEGTG - PEGF);
dxdt_PERK = KERK * (PEGF * (1.0 + KMETE * METON) - PERK);
dxdt_PAKT = KAKT * (PEGF * (1.0 + KPI3KM * PIK3M) * (2.0 - PTENF)
                    * (1.0 + KMETE * METON) - PAKT);
dxdt_STA3 = KSTA * (PEGF + KIL6S * IL6 - STA3);
double CCTG = KERKC * PERK + KAKTC * PAKT + 0.20 * gRBFREE;
dxdt_CCND = KCC * (CCTG - CCND);

// Cycle entry has three parts and only ONE of them is reachable from EGFR:
//   PBASE   a ligand- and EGFR-independent floor (no drug abolishes division)
//   CCND    the cyclin-D1 route, which IS reachable from EGFR
//   E7BYP   the Rb bypass, which in HPV-positive disease makes a further
//           fraction of cycle entry indifferent to everything upstream of Rb
// The consequence -- that anti-EGFR therapy has far less proliferative
// headroom to remove in HPV-positive disease -- is arithmetic, not a rule.
double CCNDR = posf(CCND) / CCNDREF;
double PROLD = (1.0 - gE7BYP) * (PBASE + (1.0 - PBASE) * CCNDR) + gE7BYP;

// NHEJ capacity: raised by nuclear EGFR signalling, lowered by adducts.
double NHEJT = (1.0 + KEGNH * (PEGF - 1.0)) / (1.0 + KPTNH * PTDNA);
if (NHEJT < 0.12) NHEJT = 0.12;
dxdt_NHEJ = 3.0 * (NHEJT - NHEJ);

// =====================================================================
//  6 · PERFUSION, HYPOXIA, REOXYGENATION
// =====================================================================
double HFNOW = NVIAB > 1.0 ? TUMHY / NVIAB : 0.0;
dxdt_VEGF = KVGS * (1.0 + KVGH * HFNOW) - KVGD * VEGF;
dxdt_VASC = KVIN * VEGF * (1.0 - VASC) - KVOUT * VASC - KVRT * DRATE * VASC;

double PO2 = PO2MX * posf(VASC) * gVQ / (1.0 + KDEM * NVIAB / NREF);
double HFTGT = 1.0 - hillf(PO2, PO2H, HHYP);     // hypoxic-fraction target

// Oxygen enhancement ratio, collapsed by a hypoxic radiosensitiser
double OERE = 1.0 + (OERMAX - 1.0) * (1.0 - NIMOE * NIMOF);
double OERB = 1.0 + (OERE - 1.0) * HFTGT;        // population-average OER

// =====================================================================
//  7 · EFFECTIVE RADIOSENSITIVITY AND SUBLETHAL REPAIR
// =====================================================================
double ALPHAE = gALPHA * (1.0 + KHRA * (1.0 - gHRCAP))
                      * (1.0 + PHIPT * hillf(PTDNA, KPHIPT, 1.0));
double BETAE  = ALPHAE / ABTUM;

double MUREP = NOREP > 0.5 ? 0.0 : MUREP0 * NHEJ;
dxdt_SLD  = DRATE - MUREP * SLD;
dxdt_SLDM = FMUC * DRATE - MUREPM * SLDM;
dxdt_SLDP = DRATE - MUREPP * SLDP;   // scaled per organ at point of use

// Lethal-lesion yield rates (per day) for each tumour pool
double HZOX  = ALPHAE * DRATE + 2.0 * BETAE * SLD * DRATE;
double HZHY  = (ALPHAE / OERE) * DRATE
             + 2.0 * (BETAE / (OERE * OERE)) * SLD * DRATE;
// Clonogen hazards.  OXIC and HYPOXIC clonogens are separate pools, so the
// oxygen enhancement ratio divides the dose only for the cells that are
// actually hypoxic -- and the KOXTR exchange between the pools between
// fractions IS reoxygenation.  How much protection hypoxia actually buys is
// therefore set by the race between KOXTR and the fraction interval, not by a
// parameter.
double ACO   = FCSCR * ALPHAE;
double BCO   = FCSCR * BETAE;
double ACH   = ACO / OERE;
double BCH   = BCO / (OERE * OERE);
double HZCO  = ACO * DRATE + 2.0 * BCO * SLD * DRATE;
double HZCH  = ACH * DRATE + 2.0 * BCH * SLD * DRATE;
double HZCSC = CSCT > 1e-12 ? (HZCO * CSCO + HZCH * CSCH) / CSCT : 0.0;

// =====================================================================
//  8 · CHEMOTHERAPY AND ADCC KILL
// =====================================================================
double KCIS = EMXCIS * hillf(PTDNA, EC50PT, 1.0);
double KFU  = EMXFU  * hillf(CFU,  EC50FU,  1.0);
double KDOC = EMXDOC * hillf(CDOC, EC50DOC, 1.0);
double KCHEM = KCIS + KFU + KDOC;
double KADCCR = KADCC * OCCEG * FCGR3A;

// =====================================================================
//  9 · REPOPULATION DRIVE AND CELL LOSS
//      Resolved BEFORE the immune module, because natural cell turnover is
//      the dominant source of antigen in an untreated tumour and the immune
//      module needs it as an input.
// =====================================================================
// DAMS is a slow accumulator of recent kill, so the switch from asymmetric to
// symmetric division, the loss of quiescence and the collapse of the cell-loss
// factor all appear about three weeks into a course rather than on day one.
// The kick-off time is therefore an OUTPUT of KDAMOUT, not a parameter.
double KILLRT = HZOX * TUMOX + HZHY * TUMHY;               // cells/day from RT
double KILLCH = KCHEM * TUMOX + KCHEM * FHYCH * TUMHY;     // cells/day, chemo
// DAMS is driven by the SPECIFIC (per-cell) injury rate, not by the absolute
// number of cells dying.  Driving it from the flux was wrong and produced the
// opposite of the observed phenomenon: as the tumour shrank the flux collapsed,
// so the repopulation drive PEAKED in week one and decayed thereafter.  The
// signal a surviving clonogen responds to is the hazard it is experiencing.
dxdt_DAMS = KDAMIN * posf(HZOX + KCHEM) - KDAMOUT * DAMS;
// A Hill-1 saturation cannot give both a LATE kick-off and a HIGH plateau; the
// measured kick-off at 3-4 weeks with a 0.6-0.9 Gy/day cost requires a
// sigmoid, so HDAM > 1.  KDAM50 then sets the kick-off time directly.
double REPD = hillf(DAMS, KDAM50, HDAM);
double RHO  = RHO0 + (RHOMX - RHO0) * REPD;
double KLOSS = KLOSS0 * (1.0 - FLOSSR * REPD);

// =====================================================================
//  10 · IMMUNE MODULE
// =====================================================================
// Antigen supply.  Natural turnover, radiation kill and chemotherapy kill all
// release antigen; radiation and platinum additionally make that death
// immunogenic (calreticulin exposure, ATP and HMGB1 release).  Immune kill is
// deliberately NOT fed back into antigen supply: that loop would be
// self-amplifying and the model has no data to constrain its gain.
double KILLFX = KILLRT + KILLCH + KLOSS * TUMOX;           // cells/day
double KILLN  = KILLFX / NREF;                             // normalised
double ICDON  = (DRATE > 1.0 || PTDNA > 0.10) ? 1.0 : 0.0;
double ANTG   = ANTGC * gAGVIR * (1.0 + KTMB * (PY < 60 ? PY : 60))
              + posf(KILLN) * (1.0 + FICD * ICDON);

dxdt_DCA  = KDCIN * ANTG - KDCOUT * DCA;
dxdt_TNLN = KTPRIM * DCA * (1.0 + KCTLA * CTLARX)
            - KTNOUT * TNLN - KTNEMI * TNLN;

double PDLSET = CPS / CPSREF;
dxdt_PDL1 = KPDLIN * (PDLSET + KIFNP * IFNG) - KPDLO * PDL1;
double PD1S = posf(PDL1) * (1.0 - OCCPD);        // residual PD-1 signal

dxdt_CD8T = KTNEMI * TNLN * FHOME - KT8OUT * CD8T
            - KEXH * PD1S * CD8T + KREINV * OCCPD * TEX;
dxdt_TEX  = KEXH * PD1S * CD8T - KTEXO * TEX - KREINV * OCCPD * TEX;
dxdt_IFNG = KIFN * CD8T - KIFND * IFNG;
dxdt_TREG = KTRIN * (1.0 + KIL6T * IL6 + KLACTT * HFNOW) - KTROUT * TREG;
dxdt_MDSC = KMDIN * (1.0 + KIL6M * IL6 + KST3M * posf(STA3)) - KMDOUT * MDSC;

// Suppression multiplier: HLA loss is a hard ceiling, Treg/MDSC a soft one.
double SUPP = (1.0 - HLALOSS) / (1.0 + KTSUP * (posf(TREG) + posf(MDSC)));
// KIMMB is the capacity the SAME T-cell clones would have if none of them were
// exhausted and nothing suppressed them: the exhausted pool TEX is counted back
// in, HLA loss is removed and Treg/MDSC are removed.  KIMM is what is actually
// realised.  Without counting TEX back in, HDIMM is blind to the PD-1 axis and
// no value of CPS can move it -- which is exactly the bug this replaced.
double KIMMB = KKILL * hillf(CD8T + TEX, KMT, 1.0);   // unsuppressed capacity
double KIMM  = KKILL * hillf(CD8T, KMT, 1.0) * SUPP;  // realised, bulk oxic
double KIMMS = KIMM * FCSCIM;                         // realised on clonogens

// =====================================================================
//  11 · TUMOUR POPULATION DYNAMICS
// =====================================================================
double NCAP  = NCAPF * NREF;
double CROWD = 1.0 - NTOT / NCAP;
if (CROWD < 0.02) CROWD = 0.02;

// Accelerated repopulation is TWO things, not one: cells cycle faster
// (FREPL, loss of the quiescent fraction) and divide more symmetrically
// (RHO).  Either alone is too weak to reproduce the measured 0.6-0.9 Gy/day
// that has to be spent to hold locoregional control constant.
double LAMS = gLAMSP * PROLD * CROWD * (1.0 + FREPL * REPD);
double LAMP = LAMP0  * PROLD * CROWD * (1.0 + FREPL * REPD);

double KILLSO = HZCO + KCHEM * FCSCCH + KIMMS + KADCCR * FCSCCH;
double KILLSH = HZCH + KCHEM * FCSCCH * FHYCH + KIMMS * FHYIM
                     + KADCCR * FCSCCH * FHYCH;
double KILLP = HZOX  + KCHEM          + KIMM  + KADCCR;
double KILLH = HZHY  + KCHEM * FHYCH  + KIMM * FHYIM + KADCCR * FHYCH;

// Only the oxic clonogen pool cycles; the hypoxic pool is quiescent, which is
// why it is chemo-resistant as well as radioresistant.
dxdt_CSCO  = (2.0 * RHO - 1.0) * LAMS * CSCO - KILLSO * CSCO
             - KOXTR * HFTGT * CSCO + KOXTR * (1.0 - HFTGT) * CSCH;
dxdt_CSCH  = KOXTR * HFTGT * CSCO - KOXTR * (1.0 - HFTGT) * CSCH
             - KILLSH * CSCH;
dxdt_TUMOX = 2.0 * (1.0 - RHO) * LAMS * CSCO + LAMP * TUMOX
             - KLOSS * TUMOX - KILLP * TUMOX
             - KOXTR * HFTGT * TUMOX + KOXTR * (1.0 - HFTGT) * TUMHY;
dxdt_TUMHY = KOXTR * HFTGT * TUMOX - KOXTR * (1.0 - HFTGT) * TUMHY
             - KLOSSH * TUMHY - KNEC * TUMHY - KILLH * TUMHY;
dxdt_TUMDM = KILLSO * CSCO + KILLSH * CSCH + KILLP * TUMOX
             + KILLH * TUMHY - KDOOM * TUMDM;
// Only a fraction FNEC of naturally lost cells persists as radiological
// necrotic mass; the rest is cleared by macrophages faster than the model
// resolves.  Without this the necrotic pool integrates the whole cell-loss
// flux and a 30 cm3 tumour reaches 200 cm3 of debris in a month.
dxdt_TUMNC = KDOOM * TUMDM + KNEC * TUMHY + FNEC * KLOSS * TUMOX
             - KCLRN * TUMNC;

// =====================================================================
//  11 · COUNTERFACTUAL LOG-KILL INTEGRATORS  (the headroom measurement)
// =====================================================================
// LKACT vs LKOXI isolates hypoxia: the same course delivered to a clonogen
// pool that never sees an OER.  LKACT vs LKREP isolates repair: the same
// course delivered to a clonogen with maximal HR deficit and maximal platinum
// sensitisation.  LKNRP is what repopulation hands back.  LKIMM vs LKIMX
// isolates immune escape at the CD8 level the tumour actually has.
// All three integrators are evaluated on the SAME clonogen mix that exists at
// that instant, so the only thing that differs between them is the resistance
// factor being switched off.
double FOX  = CSCT > 1e-12 ? CSCO / CSCT : 1.0;
double FHYC = 1.0 - FOX;
double ALMAX = gALPHA * (1.0 + KHRA) * (1.0 + PHIPT);
double BEMAX = ALMAX / ABTUM;
// counterfactual 1: no oxygen effect at all (hypoxic clonogens behave as oxic)
double HZOXI = ACO * DRATE + 2.0 * BCO * SLD * DRATE;
// counterfactual 2: maximal HR deficit and maximal platinum sensitisation,
// hypoxia left exactly as it is
double AMO = FCSCR * ALMAX;
double BMO = FCSCR * BEMAX;
double AMH = AMO / OERE;
double BMH = BMO / (OERE * OERE);
double HZREP = FOX  * (AMO * DRATE + 2.0 * BMO * SLD * DRATE)
             + FHYC * (AMH * DRATE + 2.0 * BMH * SLD * DRATE);

dxdt_LKACT = HZCSC;
dxdt_LKOXI = HZOXI;
dxdt_LKREP = HZREP;
dxdt_LKNRP = posf((2.0 * RHO - 1.0) * LAMS);
dxdt_LKIMM = KIMMS;
dxdt_LKIMX = KIMMB * FCSCIM;
dxdt_LKCHM = KCHEM * FCSCCH + KADCCR * FCSCCH;

// =====================================================================
//  12 · ORAL MUCOSITIS  (the competing clock)
// =====================================================================
double DMUC  = FMUC * DRATE;
double BETAM = ALPHAM / ABMUC;
double HZMUC = ALPHAM * DMUC + 2.0 * BETAM * SLDM * DMUC;

dxdt_MUCB  = KMUCP * MUCB * (1.0 - MUCB) + KMUCS * (1.0 - MUCB)
             - HZMUC * MUCB - KPTMUC * PTDNA * PTTOX * MUCB
             - KINFMU * MUCI * MUCB;
dxdt_MUCTA = KTAIN * MUCB - KTAOUT * MUCTA - HZMUC * FTAR * MUCTA;
dxdt_MUCE  = KEIN * MUCTA - KEOUT * MUCE;
double ULC = posf(1.0 - MUCE);
dxdt_MUCI  = KIIN * ULC * (1.0 + KAMP * hillf(MUCI, KAMP50, 1.0))
             - KIOUT * MUCI;

// =====================================================================
//  13 · SALIVARY GLAND, SWALLOWING, THYROID
// =====================================================================
double DPAR  = FPAR * DRATE;
double BETAP = ALPHAP / ABPAR;
double SLDPP = FPAR * SLDP;                      // organ-scaled sublethal pool
double HZPAR = ALPHAP * DPAR + 2.0 * BETAP * SLDPP * DPAR;
dxdt_ACIN = -HZPAR * ACIN + KACREG * (1.0 - ACIN) * posf(ACIN);

double TGFB = KTGFB * (posf(TREG) + posf(MDSC) + MUCI) / 3.0;
dxdt_CONF = KFIB * FCON * DRATE * (1.0 + TGFB) - KFIBR * CONF;
dxdt_TSHX = KTSHIN * (CUMD / 70.0) * (1.0 + KTHIRA * IRAE)
            - KTSHO * (TSHX - TSH0);

// =====================================================================
//  14 · SYSTEMIC TOXICITY
// =====================================================================
double CCOCH = CISCO / 1.0;                      // cochlear platinum (mg)
dxdt_HAIR = -KHAIR * pow(posf(CCOCH), HHAIR) * PTTOX * HAIR;
dxdt_KIDT = -KKID * (CISK / 1.0) * PTTOX * KIDT + KKREP * (1.0 - KIDT);
dxdt_MGSER= KMGIN * MG0 * (1.0 - KMGCTX * OCCEG - KMGCIS * (1.0 - KIDT))
            - KMGOUT * MGSER;
dxdt_RASH = KRIN * OCCEG * (4.0 - RASH) - KROUT * RASH;
dxdt_IRAE = KIRIN * OCCPD * (1.0 + KCTLA * CTLARX) - KIROUT * IRAE;

// Friberg semi-mechanistic myelosuppression
double KTR   = 4.0 / MTT;
double EDRUG = SARCOF * (SLPCIS * CFREE * PTTOX + SLPDOC * CDOC + SLPFU * CFU);
if (EDRUG > 0.95) EDRUG = 0.95;
double FBK = pow(ANC0 / (ANC > 0.02 ? ANC : 0.02), GAM);
dxdt_PROLN = KTR * PROLN * ((1.0 - EDRUG) * FBK - 1.0);
dxdt_TRN1  = KTR * (PROLN - TRN1);
dxdt_TRN2  = KTR * (TRN1 - TRN2);
dxdt_TRN3  = KTR * (TRN2 - TRN3);
dxdt_ANC   = KTR * TRN3 - KTR * ANC;

// =====================================================================
//  15 · HOST STATE AND ctDNA
// =====================================================================
dxdt_IL6 = KIL6 * (NVIAB / NREF) + KMUCIL * MUCI - KIL6D * IL6;
double MGR  = MUCGRX * pow(ULC, MUCGRP);
double NPOF = MGR > NPOTHR ? 1.0 : 0.0;
dxdt_BWT = -KWL * hillf(IL6 - IL6B, KWL50, 1.0) * BWT - KWNPO * NPOF * BWT
           + KWREC * (BW0 - BWT);
dxdt_CTD = KCTD * HPV * posf(KILLN) - KCTDCL * CTD;

$TABLE
// Every quantity below is recomputed from the state vector and the parameter
// list.  Nothing is carried over from the last derivative evaluation, so an
// output can never be a leftover from an intermediate solver step.

// ---------------------------------------------------------------- geometry
double VOL   = NTOT / CELLD;                            // cm^3
double SLDM_ = 2.0 * pow(posf(3.0 * VOL / (4.0 * 3.14159265)), 1.0/3.0) * 10.0;
double VOLF  = VOL / (V0 > 0.01 ? V0 : 0.01);
double RECST = 100.0 * (pow(VOLF, 1.0/3.0) - 1.0);      // % change in diameter

// -------------------------------------- oxygenation, recomputed from state
double tPO2   = PO2MX * posf(VASC) * gVQ / (1.0 + KDEM * NVIAB / NREF);
double tHF    = 1.0 - hillf(tPO2, PO2H, HHYP);
double tOERE  = 1.0 + (OERMAX - 1.0) * (1.0 - NIMOE * NIMOF);

// ---------------------------------------------- resistance and its headroom
double OERBAR = 1.0 + (tOERE - 1.0) * tHF;
double HDHYP  = LKOXI > 1e-9 ? 1.0 - LKACT / LKOXI : 0.0;
double HDREPR = LKREP > 1e-9 ? 1.0 - LKACT / LKREP : 0.0;
double HDPOP  = LKACT > 1e-9 ? LKNRP / LKACT       : 0.0;
double HDIMM  = LKIMX > 1e-9 ? 1.0 - LKIMM / LKIMX : 0.0;
double RESIST = OERBAR
              * (LKACT > 1e-9 ? LKREP / LKACT : 1.0)
              * (1.0 / (1.0 - (HDIMM > 0.98 ? 0.98 : HDIMM)));

// ------------------------------------------------------------- kill ledger
double TOTLK  = LKACT + LKIMM + LKCHM - LKNRP;          // nats on clonogens
double LOGK10 = (LKACT + LKIMM + LKCHM - LKNRP) / 2.302585;
double CSCLOG = log10(CSCT > 1e-6 ? CSCT : 1e-6);

// ------------------------------------------------- tumour control and stage
// TCP is the Poisson probability that THIS tumour has no surviving clonogen.
// It is a step function of log kill and therefore useless as a population
// read-out: a 0.3-log difference takes it from 0.05 to 0.95.  The population
// curve is obtained by marginalising over a lognormal clonogen number with
// SD = SIGLK, which is what makes clinical dose-response curves shallow.  Both
// are reported, because they answer different questions.
double TCP    = exp(-(CSCT > 0.0 ? CSCT : 0.0));         // Poisson, this tumour
// Marginalising exp(-N) over N ~ lognormal(median = CSC, log-SD = SIGLK):
//     P(N < 1) = Phi(-ln CSC / SIGLK)  ~=  1 / (1 + exp(ln CSC / SIGLK))
// This reads the clonogen STATE, so it keeps working after the course ends --
// a path integral would keep drifting as the survivors regrow.
double LNCSC  = log(CSCT > 1e-12 ? CSCT : 1e-12);
double LRC    = 1.0 / (1.0 + exp(LNCSC / SIGLK));
double PDMET  = DMET0 * (1.0 + 0.25 * (STG - 3.0)) * (HPV > 0.5 ? 0.62 : 1.0);
if (PDMET < 0.02) PDMET = 0.02;
double ORGINJ = (1.0 - HAIR) + (1.0 - KIDT) + posf(1.0 - ACIN) + IRAE;
double HAZ3   = OSA * 10.0 * (1.0 - LRC) + OSB * 10.0 * PDMET
              + OSSTG * (STG - 1.0) + OSECOG * ECOG
              + OSAGE * posf(AGE - 60.0) + OSCOMP * ORGINJ;
double OS3    = exp(-HAZ3);                             // 3-year OS surrogate

// ----------------------------------------------------------- normal tissue
double tULC   = posf(1.0 - MUCE);
double MUCGR  = MUCGRX * pow(tULC, MUCGRP);
double SALFLO = SALF0 * pow(posf(ACIN), SALEXP);
double XEROG  = 4.0 * (1.0 - pow(posf(ACIN), SALEXP));
double DYSPHG = 4.0 * hillf(CONF + 0.45 * tULC, 1.05, 1.6);
double HLOSS  = HLMAX * (1.0 - posf(HAIR));
double EGFRR  = CRCL0 * posf(KIDT);
double PCTWL  = 100.0 * (BW0 - BWT) / BW0;
double RASHG  = RASH;
double ANCNAD = ANC;

// ------------------------------------------------------------- exposures
double CPTOT  = (CUMPT + CISC + CISP + CISB + CISK + CISCO) / BSA;  // mg/m^2
double CFREEO = CISC / V1CIS;
double CCTXP  = CTXC / V1CTX;
double CPEMP  = PEMC / V1PEM;
double PO2OUT = tPO2;
double HFOUT  = tHF;
double OCCEGO = CTXT / (KDCTX + CTXT);
double OCCPDO = PEMT / (KDPD1 + PEMT);
double ALPHAO = gALPHA * (1.0 + KHRA * (1.0 - gHRCAP))
                      * (1.0 + PHIPT * hillf(PTDNA, KPHIPT, 1.0));
double MUREPO = NOREP > 0.5 ? 0.0 : MUREP0 * NHEJ;
double RHOO   = RHO0 + (RHOMX - RHO0) * hillf(DAMS, KDAM50, 1.0);
double NHEJO  = NHEJ;
double CSCN   = CSCT;
double CSCHF  = CSCT > 1e-12 ? CSCH / CSCT : 0.0;

$CAPTURE @annotated
VOL    : Gross tumour volume (cm^3)
SLDM_  : Longest tumour diameter surrogate (mm)
RECST  : Change in longest diameter from baseline (%)
CSCN   : Surviving clonogens (cells)
CSCHF  : Hypoxic fraction of the surviving clonogen pool
CSCLOG : log10 surviving clonogens
LOGK10 : Net log10 kill delivered to the clonogen pool
TOTLK  : Net log kill in nats (RT + immune + chemo - repopulation)
TCP    : Tumour control probability (Poisson on clonogens)
LRC    : Locoregional control probability
OS3    : 3-year overall survival surrogate
OERBAR : Population-average oxygen enhancement ratio
HDHYP  : Headroom lost to hypoxia (fraction of achievable RT log-kill)
HDREPR : Headroom lost to intact repair (fraction)
HDPOP  : Log-kill handed back by repopulation (fraction of RT log-kill)
HDIMM  : Headroom lost to immune escape (fraction)
RESIST : Composite resistance factor (product of the three)
PO2OUT : Tumour pO2 (mmHg)
HFOUT  : Hypoxic fraction
ALPHAO : Effective clonogen alpha (/Gy)
MUREPO : Effective sublethal-damage repair rate (/d)
RHOO   : Symmetric self-renewal probability
NHEJO  : NHEJ repair capacity (normalised)
MUCGR  : Oral mucositis grade (WHO 0-4)
SALFLO : Stimulated salivary flow (mL/min)
XEROG  : Xerostomia grade (0-4)
DYSPHG : Dysphagia grade (0-4)
HLOSS  : High-frequency hearing threshold shift (dB)
EGFRR  : Estimated GFR (mL/min)
ANCNAD : Absolute neutrophil count (10^9/L)
RASHG  : Acneiform rash grade (0-4)
PCTWL  : Weight loss from baseline (%)
CPTOT  : Cumulative platinum dose (mg/m^2)
CFREEO : Free platinum concentration (mg/L)
CCTXP  : Cetuximab plasma concentration (mg/L)
CPEMP  : Anti-PD-1 plasma concentration (mg/L)
OCCEGO : EGFR occupancy in tumour (fraction)
OCCPDO : PD-1 occupancy in tumour (fraction)

$OMEGA @annotated @block
ECLON : 1.40                     : ln-scale variability in clonogen number
EALPH : 0.10 0.090               : ln-scale variability in intrinsic alpha
EPERF : 0.02 0.015 0.120         : ln-scale variability in perfusion quality
EANTG : 0.00 0.000 0.000 0.250   : ln-scale variability in tumour antigenicity

$SIGMA 0

$SET end = 1825, delta = 1
)---"

## -----------------------------------------------------------------------------
##  BUILD
## -----------------------------------------------------------------------------
mod <- mcode_cache("hnscc_qsp", hnscc_code,
                   atol = 1e-8, rtol = 1e-8, maxsteps = 500000)

## =============================================================================
##  DOSING HELPERS
## =============================================================================

## Radiation fraction times.  Fractions are delivered Monday-Friday, so a
## 35-fraction course occupies 46-48 calendar days, not 35.  `gapafter` and
## `gaplen` insert an unplanned treatment break, which is the only way to
## interrogate the overall-treatment-time penalty honestly.
fx_times <- function(n, start = 0, per_day = 1, gapafter = NA, gaplen = 0) {
  out <- numeric(0); day <- start; k <- 0; pending <- gapafter
  while (k < n) {
    if ((day - start) %% 7 < 5) {
      for (j in seq_len(per_day)) {
        if (k >= n) break
        out <- c(out, day + (j - 1) * (6 / 24))   # 6 h between same-day fx
        k <- k + 1
      }
    }
    day <- day + 1
    if (!is.na(pending) && k >= pending) { day <- day + gaplen; pending <- NA }
  }
  out
}

## Radiotherapy course as an mrgsolve event object
rt_course <- function(gy_per_fx = 2, n = 35, start = 0, per_day = 1,
                      gapafter = NA, gaplen = 0) {
  tt <- fx_times(n, start, per_day, gapafter, gaplen)
  as.ev(data.frame(time = tt, amt = gy_per_fx, cmt = "DBUF", evid = 1))
}

## Cisplatin 100 mg/m2 over 2 h, days 1/22/43 of the RT course
cis_q3w <- function(mgm2 = 100, bsa = 1.8, n = 3, start = 0, ii = 21) {
  as.ev(data.frame(time = start + ii * (seq_len(n) - 1),
                   amt = mgm2 * bsa, cmt = "CISC", evid = 1,
                   rate = mgm2 * bsa / (2 / 24)))
}

## Cisplatin 40 mg/m2 weekly over 1 h
cis_wk <- function(mgm2 = 40, bsa = 1.8, n = 7, start = 0) {
  as.ev(data.frame(time = start + 7 * (seq_len(n) - 1),
                   amt = mgm2 * bsa, cmt = "CISC", evid = 1,
                   rate = mgm2 * bsa / (1 / 24)))
}

## Carboplatin AUC 5 (Calvert): dose = AUC x (GFR + 25); modelled through the
## platinum compartments with CARBOF = 1 set at the scenario level.
carbo_q3w <- function(auc = 5, gfr = 90, n = 3, start = 0, ii = 21) {
  as.ev(data.frame(time = start + ii * (seq_len(n) - 1),
                   amt = auc * (gfr + 25), cmt = "CISC", evid = 1,
                   rate = auc * (gfr + 25) / (1 / 24)))
}

## 5-FU 1000 mg/m2/day continuous infusion, days 1-4 of each cycle
fu_ci <- function(mgm2 = 1000, bsa = 1.8, days = 4, n = 6, start = 0, ii = 21) {
  as.ev(data.frame(time = start + ii * (seq_len(n) - 1),
                   amt = mgm2 * bsa * days, cmt = "FUC", evid = 1,
                   rate = mgm2 * bsa))
}

## Docetaxel 75 mg/m2 over 1 h
doce_q3w <- function(mgm2 = 75, bsa = 1.8, n = 3, start = 0, ii = 21) {
  as.ev(data.frame(time = start + ii * (seq_len(n) - 1),
                   amt = mgm2 * bsa, cmt = "DOCC", evid = 1,
                   rate = mgm2 * bsa / (1 / 24)))
}

## Cetuximab: 400 mg/m2 loading one week before RT, then 250 mg/m2 weekly
cetux <- function(bsa = 1.8, n = 8, start = -7) {
  load <- data.frame(time = start, amt = 400 * bsa, cmt = "CTXC", evid = 1,
                     rate = 400 * bsa / (2 / 24))
  wk   <- data.frame(time = start + 7 * seq_len(n), amt = 250 * bsa,
                     cmt = "CTXC", evid = 1, rate = 250 * bsa / (1 / 24))
  as.ev(rbind(load, wk))
}

## Pembrolizumab 200 mg q3w over 30 min
pembro <- function(mg = 200, n = 35, start = 0, ii = 21) {
  as.ev(data.frame(time = start + ii * (seq_len(n) - 1), amt = mg,
                   cmt = "PEMC", evid = 1, rate = mg / (0.5 / 24)))
}

## Nivolumab 3 mg/kg q2w over 30 min
nivo <- function(mgkg = 3, bw = 70, n = 52, start = 0, ii = 14) {
  as.ev(data.frame(time = start + ii * (seq_len(n) - 1), amt = mgkg * bw,
                   cmt = "PEMC", evid = 1, rate = mgkg * bw / (0.5 / 24)))
}

## =============================================================================
##  TWENTY-FOUR TREATMENT SCENARIOS
##
##  The suite is built so that every claim in the header has a matched control.
##  In particular:
##
##    * 04 vs 05 vs 06        identical 70 Gy; cisplatin vs cetuximab vs neither
##    * 07 vs 08 vs 09        the SAME three arms with HPV = 1 -- this is the
##                            RTOG 1016 / De-ESCALaTE contrast, and the only
##                            thing that changed is the six phenotype inputs
##    * 12                    cisplatin + cetuximab + RT (RTOG 0522) -- both
##                            agents on the SAME resistance factor
##    * 13 vs 14              nimorazole in a hypoxic and a well-perfused
##                            tumour at identical dose (DAHANCA 5 geometry)
##    * 15 vs 16              identical prescription, 10-day treatment gap
##    * 17-20                 identical pembrolizumab exposure, CPS 1/10/20/50
##    * 21                    pembrolizumab at fixed CPS but HPV switched, to
##                            show the immune axis is orthogonal to the others
##    * 22                    IMRT vs 3D-CRT parotid dose at identical tumour
##                            dose (toxicity contrast only)
##    * 23                    weekly vs three-weekly cisplatin (same total)
##    * 24                    hyperfractionation, 1.2 Gy BID -- the incomplete
##                            repair penalty is an output of MUREP
## =============================================================================

## NOTE.  zero_re() is used everywhere.  These are DETERMINISTIC contrasts; the
## whole point of 04-vs-07 is that nothing differs except the phenotype, and a
## random draw would let luck masquerade as mechanism.
run_scn <- function(label, dose = NULL, param = list(), end = 1825, delta = 1) {
  m <- zero_re(mod)
  if (length(param)) m <- param(m, param)
  out <- if (is.null(dose)) mrgsim_df(m, end = end, delta = delta)
         else              mrgsim_df(m, events = dose, end = end, delta = delta)
  out$scenario <- label
  out
}

## The two phenotypes.  NOTE: always extend these with modifyList(), never with
## c().  c(HPVNEG, list(CPS = 50)) produces a list with TWO elements named CPS
## and param() silently keeps the first, so the override is discarded without
## any error -- which is exactly how the CPS-gradient scenarios came to be
## running four identical simulations during development.
HPVPOS <- list(HPV = 1, PY = 5, CPS = 20, V0 = 22, HLALOSS = 0.08)
HPVNEG <- list(HPV = 0, PY = 40, CPS = 10, V0 = 32, HLALOSS = 0.18)

scenarios <- list(

  ## ------------------------------------------------- natural history arms
  s01 = function() run_scn("01 Untreated natural history (HPV-negative)",
                           param = HPVNEG, end = 540),
  s02 = function() run_scn("02 Untreated natural history (HPV-positive)",
                           param = HPVPOS, end = 540),

  ## -------------------------------------------- radiation alone, HPV-neg
  s03 = function() run_scn("03 RT alone 70 Gy/35 fx (HPV-negative)",
                           rt_course(2, 35), param = HPVNEG),

  ## ------------------ the three definitive arms, HPV-negative disease
  s04 = function() run_scn("04 CRT cisplatin 100 q3w x3 + 70 Gy (HPV-negative)",
                           c(rt_course(2, 35), cis_q3w(100, 1.8, 3, 0)),
                           param = HPVNEG),
  s05 = function() run_scn("05 Cetuximab-RT 70 Gy (HPV-negative, Bonner)",
                           c(rt_course(2, 35), cetux(1.8, 8, -7)),
                           param = HPVNEG),
  s06 = function() run_scn("06 RT alone 70 Gy (HPV-negative, comparator)",
                           rt_course(2, 35), param = HPVNEG),

  ## ------------------ the SAME three arms, HPV-positive disease
  ##  Nothing changes except the six phenotype inputs.  The de-escalation
  ##  result is read off as (07 - 09) versus (08 - 09).
  s07 = function() run_scn("07 CRT cisplatin 100 q3w x3 + 70 Gy (HPV-positive)",
                           c(rt_course(2, 35), cis_q3w(100, 1.8, 3, 0)),
                           param = HPVPOS),
  s08 = function() run_scn("08 Cetuximab-RT 70 Gy (HPV-positive, RTOG 1016)",
                           c(rt_course(2, 35), cetux(1.8, 8, -7)),
                           param = HPVPOS),
  s09 = function() run_scn("09 RT alone 70 Gy (HPV-positive, comparator)",
                           rt_course(2, 35), param = HPVPOS),

  ## ------------------------------------------------- de-escalated dose
  s10 = function() run_scn("10 De-escalated CRT 60 Gy + cisplatin (HPV-positive)",
                           c(rt_course(2, 30), cis_q3w(100, 1.8, 2, 0)),
                           param = HPVPOS),

  ## --------------------------------- induction chemotherapy then RT (TPF)
  s11 = function() run_scn("11 Induction TPF x3 then 70 Gy (HPV-negative)",
                           c(doce_q3w(75, 1.8, 3, 0), cis_q3w(75, 1.8, 3, 0),
                             fu_ci(1000, 1.8, 4, 3, 0), rt_course(2, 35, 70)),
                           param = HPVNEG),

  ## -------------- two agents on the SAME resistance factor (RTOG 0522)
  s12 = function() run_scn("12 Cisplatin + cetuximab + 70 Gy (RTOG 0522)",
                           c(rt_course(2, 35), cis_q3w(100, 1.8, 3, 0),
                             cetux(1.8, 8, -7)),
                           param = HPVNEG),

  ## ------------------ an agent on an ORTHOGONAL factor: hypoxia
  s13 = function() run_scn("13 RT 70 Gy + nimorazole (hypoxic HPV-negative)",
                           rt_course(2, 35),
                           param = modifyList(HPVNEG, list(NIMOF = 1))),
  s14 = function() run_scn("14 RT 70 Gy + nimorazole (well-perfused HPV-pos)",
                           rt_course(2, 35),
                           param = modifyList(HPVPOS, list(NIMOF = 1))),

  ## ---------------------------------- overall treatment time (10-day gap)
  s15 = function() run_scn("15 CRT 70 Gy, no interruption (reference)",
                           c(rt_course(2, 35), cis_q3w(100, 1.8, 3, 0)),
                           param = HPVNEG),
  s16 = function() run_scn("16 CRT 70 Gy with a 10-day treatment gap",
                           c(rt_course(2, 35, gapafter = 20, gaplen = 10),
                             cis_q3w(100, 1.8, 3, 0)),
                           param = HPVNEG),

  ## -------------------- checkpoint blockade across CPS (KEYNOTE-048/040)
  ##  Identical pembrolizumab exposure in all four; only CPS differs.
  s17 = function() run_scn("17 Pembrolizumab monotherapy, CPS 1 (R/M)",
                           pembro(200, 35, 0),
                           param = modifyList(HPVNEG, list(CPS = 1, V0 = 60)),
                           end = 1095),
  s18 = function() run_scn("18 Pembrolizumab monotherapy, CPS 10 (R/M)",
                           pembro(200, 35, 0),
                           param = modifyList(HPVNEG, list(CPS = 10, V0 = 60)),
                           end = 1095),
  s19 = function() run_scn("19 Pembrolizumab monotherapy, CPS 20 (R/M)",
                           pembro(200, 35, 0),
                           param = modifyList(HPVNEG, list(CPS = 20, V0 = 60)),
                           end = 1095),
  s20 = function() run_scn("20 Pembrolizumab monotherapy, CPS 50 (R/M)",
                           pembro(200, 35, 0),
                           param = modifyList(HPVNEG, list(CPS = 50, V0 = 60)),
                           end = 1095),

  ## ------------- the orthogonality test: same CPS, different HPV status
  s21 = function() run_scn("21 Pembrolizumab, CPS 20, HPV-positive (R/M)",
                           pembro(200, 35, 0),
                           param = modifyList(HPVPOS, list(CPS = 20, V0 = 60)),
                           end = 1095),

  ## ------------------- EXTREME regimen: platinum/5-FU + cetuximab (R/M)
  s22 = function() run_scn("22 EXTREME: cisplatin/5-FU + cetuximab (R/M)",
                           c(cis_q3w(100, 1.8, 6, 0), fu_ci(1000, 1.8, 4, 6, 0),
                             cetux(1.8, 52, 0)),
                           param = modifyList(HPVNEG, list(V0 = 60)), end = 1095),

  ## ----------------- toxicity-only contrast at identical tumour dose
  s23 = function() run_scn("23 CRT 70 Gy with 3D-CRT parotid dose (FPAR 0.80)",
                           c(rt_course(2, 35), cis_q3w(100, 1.8, 3, 0)),
                           param = modifyList(HPVNEG, list(FPAR = 0.80))),

  ## ----------------- schedule contrasts at (nearly) matched total dose
  s24 = function() run_scn("24 CRT weekly cisplatin 40 mg/m2 x7 + 70 Gy",
                           c(rt_course(2, 35), cis_wk(40, 1.8, 7, 0)),
                           param = HPVNEG),
  s25 = function() run_scn("25 Hyperfractionation 1.2 Gy BID to 81.6 Gy",
                           c(rt_course(1.2, 68, per_day = 2),
                             cis_q3w(100, 1.8, 3, 0)),
                           param = HPVNEG),

  ## --------------------------------- renal impairment / carboplatin swap
  s26 = function() run_scn("26 Carboplatin AUC5 + 70 Gy (CrCl 55)",
                           c(rt_course(2, 35), carbo_q3w(5, 55, 3, 0)),
                           param = modifyList(HPVNEG, list(CARBOF = 1, CRCL0 = 55)))
)

run_all <- function() do.call(rbind, lapply(scenarios, function(f) f()))

## =============================================================================
##  READ-OUTS
## =============================================================================

## Nearest-row index helper
i_at <- function(x, t) which.min(abs(x$time - t))

## Index of the clonogen nadir WITHIN A FIXED WINDOW.
##
## Two things go wrong if this window is left open.  Read the control
## probability at a fixed late time and it is dragged down by the regrowth of
## the fraction that was never controlled.  Read it at the global nadir and the
## horizon itself becomes a covariate: in the deeply-responding arms the
## residual immune kill keeps grinding the clonogen pool down for years, the
## nadir slides past the end of any clinically meaningful follow-up, and arms
## that a trial separates cleanly converge.  Clinical locoregional failure is
## declared within about two years, and essentially all of the separation is
## established in the first, so the window is one year.
NADIR_WINDOW <- 365
nadir_idx <- function(x, tmax = NADIR_WINDOW) {
  k <- which(x$time <= tmax)
  if (!length(k)) k <- seq_len(nrow(x))
  k[which.min(x$CSCN[k])]
}
peak_lk <- function(x, tmax = NADIR_WINDOW) {
  k <- which(x$time <= tmax)
  if (!length(k)) k <- seq_len(nrow(x))
  max(x$LOGK10[k])
}

summarise_scn <- function(all) {
  do.call(rbind, lapply(split(all, all$scenario), function(x) {
    ## Read everything at the CLONOGEN NADIR, not at a fixed late time.
    ## Locoregional control is decided by how few clonogens survive the course;
    ## after the nadir the deterministic trajectory shows the regrowth of the
    ## fraction that was NOT controlled, so a fixed-time read-out reports
    ## something far below the actuarial control rate the trials measure.
    j    <- nadir_idx(x)
    jend <- nrow(x)
    data.frame(
      scenario = x$scenario[1],
      logkill  = round(peak_lk(x), 2),
      clonogen = signif(x$CSCN[j], 3),
      LRC      = round(x$LRC[j], 3),
      OS3      = round(x$OS3[j], 3),
      bestRECIST = round(min(x$RECST[x$time <= NADIR_WINDOW]), 1),
      HD_hyp   = round(x$HDHYP[j], 3),
      HD_rep   = round(x$HDREPR[j], 3),
      HD_pop   = round(x$HDPOP[j], 3),
      HD_imm   = round(x$HDIMM[j], 3),
      mucGRmax = round(max(x$MUCGR), 2),
      salivend = round(x$SALFLO[jend], 3),
      heardB   = round(max(x$HLOSS), 1),
      eGFRmin  = round(min(x$EGFRR), 1),
      ANCmin   = round(min(x$ANCNAD), 2),
      wtloss   = round(max(x$PCTWL), 1),
      cumPt    = round(max(x$CPTOT), 0),
      row.names = NULL)
  }))
}

## The central claim, stated as an arithmetic identity rather than a sentence:
## the log-kill a resistance-directed agent buys should scale with the headroom
## its factor had in the arm WITHOUT that agent.
headroom_test <- function(all) {
  g <- function(lab, col, t = 50) {
    x <- all[all$scenario == lab, ]
    x[[col]][nadir_idx(x)]
  }
  lk <- function(lab) peak_lk(all[all$scenario == lab, ])
  data.frame(
    contrast = c("cisplatin added to RT (HPV-neg)",
                 "cetuximab added to RT (HPV-neg)",
                 "cisplatin added to RT (HPV-pos)",
                 "cetuximab added to RT (HPV-pos)",
                 "cetuximab added to cisplatin-RT (HPV-neg)",
                 "nimorazole added to RT (HPV-neg, hypoxic)",
                 "nimorazole added to RT (HPV-pos, oxic)"),
    factor_targeted = c("repair", "repair", "repair", "repair",
                        "repair (already spent)", "hypoxia", "hypoxia"),
    headroom_in_control = round(c(
      g("06 RT alone 70 Gy (HPV-negative, comparator)", "HDREPR"),
      g("06 RT alone 70 Gy (HPV-negative, comparator)", "HDREPR"),
      g("09 RT alone 70 Gy (HPV-positive, comparator)", "HDREPR"),
      g("09 RT alone 70 Gy (HPV-positive, comparator)", "HDREPR"),
      g("04 CRT cisplatin 100 q3w x3 + 70 Gy (HPV-negative)", "HDREPR"),
      g("06 RT alone 70 Gy (HPV-negative, comparator)", "HDHYP"),
      g("09 RT alone 70 Gy (HPV-positive, comparator)", "HDHYP")), 3),
    delta_log10_kill = round(c(
      lk("04 CRT cisplatin 100 q3w x3 + 70 Gy (HPV-negative)") -
        lk("06 RT alone 70 Gy (HPV-negative, comparator)"),
      lk("05 Cetuximab-RT 70 Gy (HPV-negative, Bonner)") -
        lk("06 RT alone 70 Gy (HPV-negative, comparator)"),
      lk("07 CRT cisplatin 100 q3w x3 + 70 Gy (HPV-positive)") -
        lk("09 RT alone 70 Gy (HPV-positive, comparator)"),
      lk("08 Cetuximab-RT 70 Gy (HPV-positive, RTOG 1016)") -
        lk("09 RT alone 70 Gy (HPV-positive, comparator)"),
      lk("12 Cisplatin + cetuximab + 70 Gy (RTOG 0522)") -
        lk("04 CRT cisplatin 100 q3w x3 + 70 Gy (HPV-negative)"),
      lk("13 RT 70 Gy + nimorazole (hypoxic HPV-negative)") -
        lk("06 RT alone 70 Gy (HPV-negative, comparator)"),
      lk("14 RT 70 Gy + nimorazole (well-perfused HPV-pos)") -
        lk("09 RT alone 70 Gy (HPV-positive, comparator)")), 3),
    row.names = NULL)
}

## Overall-treatment-time penalty: the classic "1 % of locoregional control
## per day of prolongation" figure, recovered as a difference rather than
## imposed as a parameter.
ott_penalty <- function(all) {
  a <- all[all$scenario == "15 CRT 70 Gy, no interruption (reference)", ]
  b <- all[all$scenario == "16 CRT 70 Gy with a 10-day treatment gap", ]
  ja <- nadir_idx(a); jb <- nadir_idx(b)
  data.frame(
    LRC_no_gap  = round(a$LRC[ja], 3),
    LRC_10d_gap = round(b$LRC[jb], 3),
    dLRC_per_day = round((a$LRC[ja] - b$LRC[jb]) / 10 * 100, 2),
    dlogkill_per_day = round((peak_lk(a) - peak_lk(b)) / 10, 3),
    Gy_equivalent_per_day = round((peak_lk(a) - peak_lk(b)) /
                                    10 * 2.302585 / 0.2484, 3),
    row.names = NULL)
}

## Virtual-population responder rates.
##
##  IMPORTANT: use TCP here, NOT LRC.  LRC already contains between-patient
##  heterogeneity, because it is exp(-N) marginalised over a lognormal clonogen
##  number with SD = SIGLK.  The $OMEGA block puts that heterogeneity into the
##  simulation explicitly instead.  Using both at once counts it twice and
##  produces a dose-response curve shallower than anything ever measured.
##  So: deterministic contrasts -> zero_re() + LRC;
##      virtual populations     -> $OMEGA + TCP.
responder_rate <- function(n = 200, pp = HPVNEG, dose = NULL, end = 400,
                           seed = 20260730) {
  set.seed(seed)
  m <- param(mod, pp)
  out <- if (is.null(dose)) mrgsim_df(m, end = end, delta = 4, nid = n)
         else mrgsim_df(m, events = dose, end = end, delta = 4, nid = n)
  per <- do.call(rbind, lapply(split(out, out$ID), function(x) data.frame(
    ID = x$ID[1], logkill = max(x$LOGK10), tcp = max(x$TCP),
    best_recist = min(x$RECST), mucGR = max(x$MUCGR))))
  data.frame(
    n              = n,
    logkill_median = round(median(per$logkill), 2),
    logkill_iqr    = paste(round(quantile(per$logkill, c(.25, .75)), 2), collapse = "-"),
    control_rate   = round(mean(per$tcp > 0.5), 3),
    ORR_RECIST     = round(mean(per$best_recist <= -30), 3),
    mucositis_g3p  = round(mean(per$mucGR >= 3), 3),
    row.names = NULL)
}

## Checkpoint benefit versus CPS, at identical drug exposure.
cps_gradient <- function(all) {
  labs <- c("17 Pembrolizumab monotherapy, CPS 1 (R/M)",
            "18 Pembrolizumab monotherapy, CPS 10 (R/M)",
            "19 Pembrolizumab monotherapy, CPS 20 (R/M)",
            "20 Pembrolizumab monotherapy, CPS 50 (R/M)",
            "21 Pembrolizumab, CPS 20, HPV-positive (R/M)")
  do.call(rbind, lapply(labs, function(l) {
    x <- all[all$scenario == l, ]
    j <- i_at(x, 365)
    data.frame(scenario = l,
               CPS_arm = c(1, 10, 20, 50, 20)[match(l, labs)],
               HD_imm = round(x$HDIMM[nadir_idx(x)], 3),
               PD1occ = round(max(x$OCCPDO), 3),
               immune_logkill = round(peak_lk(x), 3),
               volume_1yr = round(x$VOL[j], 1),
               RECIST_1yr = round(x$RECST[j], 1),
               row.names = NULL)
  }))
}

if (identical(environment(), globalenv()) &&
    !is.null(getOption("hnscc.run.scenarios"))) {
  all <- run_all()
  print(summarise_scn(all))
  print(headroom_test(all))
  print(ott_penalty(all))
  print(cps_gradient(all))
}
