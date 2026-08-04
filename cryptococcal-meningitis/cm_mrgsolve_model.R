## =====================================================================
##  HIV-associated Cryptococcal Meningitis (CM)
##  Quantitative Systems Pharmacology model -- mrgsolve implementation
## ---------------------------------------------------------------------
##  53 ODE compartments · 173 parameters · 16 treatment scenarios
##
##  WHAT THIS MODEL IS FOR
##  ----------------------
##  Every randomised trial in this disease has used a fungal endpoint
##  (early fungicidal activity, EFA) as its pharmacodynamic read-out and
##  10-week all-cause mortality as its primary endpoint.  This model
##  exists to make the gap between those two visible, because the disease
##  runs on TWO CLOCKS with a ~20-fold difference in time constant:
##
##    CLOCK 1  viable yeast in CSF                 t1/2 ~ 0.5-2 days on
##             (Fe, Fres, Ft, Fi, Fp)              effective therapy
##    CLOCK 2  capsular polysaccharide already     t1/2 ~ 13 days,
##             in the CSF (GXM)                    NO drug acts on it
##
##  Clock 2 sets CSF outflow resistance -- hence intracranial pressure --
##  and it is the antigen stock that immune reconstitution reacts to.
##  So a patient can be culture-negative and still be generating pressure
##  and IRIS drive.  EFA cannot see any of that.
##
##  TWO SHARED, OPPOSITELY-SIGNED NODES do most of the work:
##
##    ERG  (membrane ergosterol).  Fluconazole depletes it via ERG11.
##         That is how the azole works (growth ∝ ERG^gERG) and it is also
##         how the azole removes amphotericin's binding target
##         (EC50_AmB ∝ 1/ERG^aERG).  This is why flucytosine and
##         fluconazole are NOT interchangeable partners for amphotericin,
##         and the model derives ACTA's partner-drug result from it
##         without any mortality term for "partner drug".
##
##    GFR. Amphotericin lowers it; flucytosine is cleared by it.  The
##         nephrotoxicity of one drug is therefore the marrow toxicity of
##         the other -- which is why the same 5FC dose is more myelotoxic
##         beside amphotericin than beside fluconazole.
##
##  UNITS:  time = days; drug amounts = mg; concentrations = ug/mL;
##          fungal burden = CFU/mL of CSF; GXM = ug/mL; CSF pressure =
##          cmH2O internally and mmH2O in the read-outs (1 cmH2O = 10
##          mmH2O); CSF volumes mL, flows mL/day; Rout = cmH2O per
##          (mL/day) -- multiply by 1440 for clinical cmH2O/(mL/min).
##
##  VERIFICATION -- TWO INDEPENDENT IMPLEMENTATIONS
##  ------------------------------------------------
##  Every equation below is also implemented, independently and with no
##  shared code, in `cm_reference_model.py` (pure-Python fixed-step RK4,
##  standard library only, no dependencies).  Both were executed and
##  compared.  Representative agreement, mrgsolve (LSODA) versus Python
##  (RK4, h = 0.0025 d):
##
##                                   mrgsolve      python
##      presenting CSF GXM             91.5        91.78   ug/mL
##      presenting opening pressure   250          250.1   mmH2O
##      EFA, AmB-d 1 mg/kg alone       -0.3133     -0.313  log10/day
##      EFA, AMBITION regimen          -0.4016     -0.402  log10/day
##      EFA, AmB-d 1 mg/kg + 5FC       -0.4482     -0.448  log10/day
##      10-week mortality, AMBITION     26.09%      26.1%
##      10-week mortality, 2-wk AmB     31.73%      31.8%
##      day of first negative culture   13.50       13.5   days
##      peak plasma 5FC, 2-wk AmB+5FC   53.28       53.3   ug/mL
##      free brain sertraline           33.85       33.8   ng/mL
##      26-week mortality, ART day 7    61.7%       61.6%
##      26-week mortality, ART day 35   45.6%       45.6%
##      high-burden phenotype, GXM     955.9       958     ug/mL
##
##  The cross-check was not decorative.  It exposed three real defects:
##    (1) two `for` loops in $ODE declaring the same index name, which
##        mrgsolve hoists into one scope and g++ then rejects;
##    (2) LSODA taking flucytosine very slightly negative after the last
##        dose, so that pow(Cfc, 1.8) in the marrow-suppression term
##        returned NaN and silently destroyed the mortality output of
##        every arm that sterilised.  All fractional powers are now
##        guarded.  The fixed-step Python integrator never hit this, so
##        the bug was only visible from the other side.
##    (3) the Python burn-in was 70 days and had NOT converged; the R
##        implementation, started from the Python state and run again,
##        drifted to exactly Python's 140-day values.  Rather than paper
##        over it, the 70-day clamp is now documented as the modelling
##        choice it is, and the 6 mmH2O of extra pressure that full
##        convergence adds is reported as a result (see $MAIN).
##
##  One intentional difference remains: the therapeutic lumbar puncture is
##  a proportional drain over a 30-minute window here and an instantaneous
##  pressure-target reset in Python.  This is the only place the two
##  implementations disagree materially -- 138.1 mL drained over 7
##  punctures versus 140.0, and 19.7% versus 19.4% 10-week mortality.
##
##  CALIBRATION -- what was fitted to what
##  --------------------------------------
##  PK: literature values (AmB-d Vss ~4 L/kg; L-AmB Cmax ~120 ug/mL after
##      10 mg/kg with t1/2 ~10 d; 5FC t1/2 ~3-5 h, CSF:plasma 0.74;
##      fluconazole t1/2 ~30 h, CSF:plasma 0.8).
##  CNS delivery + azole antagonism: THREE parameters (CLbrD, aERG,
##      CLbrL) solved by sequential 1-D bisection against three EFA
##      values that each isolate one of them:
##        AmB-d 1 mg/kg alone      EFA -0.31  (Day 2013 group 1) -> CLbrD
##        AmB-d 1 mg/kg + FLU 800  EFA -0.32  (Day 2013 group 3) -> aERG
##        L-AmB 10 mg/kg single    EFA -0.40  (AMBITION)         -> CLbrL
##      The remaining six EFA values are then PREDICTIONS:
##        model / published:  -0.10/-0.07 (FLU 800), -0.13/-0.11 (FLU
##        1200), -0.27/-0.28 (FLU+5FC), -0.38/-0.40 (AmB 0.7+5FC),
##        -0.45/-0.42..-0.56 (AmB 1.0+5FC), -0.36/-0.42 (1-wk AmB+5FC).
##      The last is the model's worst EFA residual and is flagged.
##  Hazard: two-stage.  The seven coefficients no induction trial can
##      identify (pressure, perfusion, anaemia, neutropenia, potassium,
##      GFR, IRIS) are FIXED at mechanistically anchored values; only
##      hB, hAMS and hSTER are fitted, by non-negative least squares, to
##      17 mortality endpoints from 5 randomised trials.  Mean absolute
##      error at the primary 10-week endpoint: 2.1 percentage points.
##      The model over-predicts 2-week mortality by 3-4 points; the cause
##      is identified (the burden hazard is instantaneous rather than
##      integrated) and is not fixed.
## =====================================================================

library(mrgsolve)
library(dplyr)

cm_code <- '
$PARAM @annotated
// ---------------- body ----------------
WT     :  60.0 : body weight (kg)
MAP    :  85.0 : mean arterial pressure (mmHg)

// ---------------- amphotericin B deoxycholate PK (2-cpt) ----------------
VcD    :  30.0 : central volume, AmB-d (L)
VpD    : 210.0 : peripheral volume, AmB-d (L)
CLD    :  34.6 : clearance, AmB-d (L/day)
QD     :  40.0 : intercompartmental clearance, AmB-d (L/day)

// ---------------- liposomal amphotericin B PK (2-cpt) ----------------
VcL    :   5.0 : central volume, L-AmB (L)
VpL    :  12.0 : peripheral volume, L-AmB (L)
CLL    :   1.20 : clearance, L-AmB (L/day)
QL     :   3.00 : intercompartmental clearance, L-AmB (L/day)

// ---------------- CNS and renal distribution of amphotericin ----------------
Vbr    :   1.30 : CNS distribution volume for AmB (L)
CLbrD  :   0.11250 : plasma clearance into CNS, deoxycholate (L/day) [SOLVED]
CLbrL  :   0.003090 : plasma clearance into CNS, liposomal (L/day) [SOLVED]
kbrOut :   0.160 : CNS efflux of AmB (1/day)
Vkid   :   0.30 : renal cortical volume (L)
CLkdD  :   0.900 : plasma clearance into renal cortex, deoxycholate (L/day)
CLkdL  :   0.0125 : plasma clearance into renal cortex, liposomal (L/day)
kkdOut :   0.055 : renal cortical efflux of AmB (1/day)

// ---------------- flucytosine PK ----------------
kaFC   :  30.0 : absorption rate, 5FC (1/day)
FFC    :   0.87 : oral bioavailability, 5FC
VFC    :  45.0 : volume of distribution, 5FC (L)
CLFC0  : 156.0 : clearance at normal GFR, 5FC (L/day)
fcCSF  :   0.74 : CSF:plasma ratio, 5FC
keqFC  :  12.0 : CSF equilibration rate, 5FC (1/day)

// ---------------- fluconazole PK ----------------
kaFL   :  18.0 : absorption rate, fluconazole (1/day)
FFL    :   0.90 : oral bioavailability, fluconazole
VFL    :  45.0 : volume of distribution, fluconazole (L)
CLFL   :  25.0 : clearance, fluconazole (L/day)
flCSF  :   0.80 : CSF:plasma ratio, fluconazole
keqFL  :   3.00 : CSF equilibration rate, fluconazole (1/day)

// ---------------- dexamethasone and sertraline PK ----------------
kaDX   :  25.0 : absorption rate, dexamethasone (1/day)
FDX    :   0.80 : oral bioavailability, dexamethasone
VDX    :  60.0 : volume of distribution, dexamethasone (L)
CLDX   :  17.0 : clearance, dexamethasone (L/day)
KDX    :   0.006 : concentration for 50% glucocorticoid effect (ug/mL)
kaSR   :  15.0 : absorption rate, sertraline (1/day)
FSR    :   0.44 : oral bioavailability, sertraline
VSR    : 1200.0 : volume of distribution, sertraline (L)
CLSR   : 1200.0 : clearance, sertraline (L/day)
kSRbr  :  60.0 : brain uptake, sertraline (L/day)
kSRbrOut : 3.00 : brain efflux, sertraline (1/day)
fuSRbr :   0.015 : free fraction of sertraline in brain tissue

// ---------------- exogenous interferon gamma ----------------
kaIFN  :   1.20 : SC absorption, IFN-gamma (1/day)
FIFN   :   0.60 : SC bioavailability, IFN-gamma
kIFNexo:  48.0 : CSF IFN-gamma per ug absorbed per day (pg/mL per ug/day)

// ---------------- fungal growth and compartment exchange ----------------
g0     :   0.350 : intrinsic net growth rate in CSF (1/day)
Fmax   :   1.0e7 : CSF carrying capacity (CFU/mL)
gp     :   0.160 : parenchymal reservoir growth rate (1/day)
Fpmax  :   3.0e6 : parenchymal carrying capacity (CFU/mL equivalent)
kSeed  :   0.0060 : CSF to parenchyma seeding (1/day)
kShed  :   0.0006 : parenchyma to CSF shedding (1/day)
kPer   :   0.0035 : conversion to persister phenotype (1/day)
kRev   :   0.045 : reversion from persister phenotype (1/day)
kPhag  :   0.075 : phagocytic uptake per unit MAC (1/day)
kEsc   :   0.085 : escape from macrophage (1/day)
kMut   :   2.2e-8 : mutation rate to 5FC resistance (1/day)

// ---------------- amphotericin B pharmacodynamics ----------------
// In the achievable CNS range the kill rate is nearly FIRST ORDER in
// concentration.  That is not an aesthetic choice: it is what the
// observed proportionality between AmB dose and EFA (0.7 vs 1.0 mg/kg,
// Bicanic 2008: -0.45 vs -0.56) forces.  A saturating fit cannot make a
// 1.43x dose step produce a 1.24x rate step.  This is the most
// exposed PD assumption in the model.
KmaxA  :   9.80 : maximum AmB kill rate (1/day, ~2.7 log10/day)
EC50A  :   2.40 : AmB EC50 at ergosterol-replete membrane (ug/mL)
hA     :   1.00 : Hill coefficient, AmB

// ---------------- flucytosine pharmacodynamics ----------------
KmaxF  :   0.72 : maximum 5FC kill rate (1/day)
EC50F  :  22.0 : 5FC EC50 in CSF (ug/mL)
RESFAC :  45.0 : EC50 multiplier in the 5FC-resistant subclone

// ---------------- fluconazole / ergosterol pharmacodynamics ----------------
kERG   :   0.85 : membrane ergosterol turnover (1/day)
IC50erg:  16.0 : CSF fluconazole for 50% ergosterol depletion (ug/mL)
KmaxL  :   0.344 : kill rate attributable to ergosterol depletion (1/day)
gERG   :   1.00 : exponent, growth rate scales as ERG^gERG
aERG   :   0.5437 : exponent, AmB EC50 scales as 1/ERG^aERG [SOLVED]

// ---------------- sertraline pharmacodynamics (from in-vitro MIC) ----------------
KmaxSR :   0.55 : maximum sertraline kill rate (1/day)
EC50SR :   6.00 : sertraline EC50, free concentration (ug/mL)

// ---------------- capsule-mediated drug tolerance ----------------
kCap   :   1.10 : maximum EC50 multiplier from capsule
KCap   : 300.0 : CSF GXM for half-maximal capsule effect (ug/mL)
tolPers:   0.40 : fraction of drug effect retained against persisters
fIntra :   0.22 : fraction of AmB effect against intracellular yeast
fParen :   0.40 : fraction of drug effect in the parenchymal reservoir

// ---------------- immune killing ----------------
kImm   :   0.0345 : extracellular killing per unit MAC (1/day)
KIFN   : 180.0 : CSF IFN-gamma for half-maximal enhancement (pg/mL)
kImmI  :   0.34 : intracellular killing per unit MAC (1/day)
kImmP  :   0.100 : parenchymal killing per unit MAC (1/day)

// ---------------- capsular antigen (GXM) ----------------
sLive  :   5.00e-6 : shedding by live yeast (ug/mL/day per CFU/mL)
sLysis :   4.50e-4 : capsule released per killed cell (ug per CFU)
kGXM   :   0.055 : GXM clearance (1/day, t1/2 = 12.6 d)
KGsat  : 800.0 : saturation of capsule production (ug/mL)

// ---------------- host immunity ----------------
KAg    : 400.0 : CSF GXM for half-maximal antigen drive (ug/mL)
KCD4   : 110.0 : CD4 for half-maximal immune competence (cells/uL)
kMACon :   1.55 : macrophage activation by antigen (1/day)
kMACifn:   0.30 : macrophage activation by IFN-gamma (1/day)
KIFN2  : 250.0 : IFN-gamma for half-maximal macrophage activation (pg/mL)
kMACoff:   0.16 : macrophage deactivation (1/day)
kTHon  :   2.30 : Th1 induction (1/day)
kTHoff :   0.22 : Th1 decay (1/day)
kIFNp  : 210.0 : IFN-gamma production per unit Th1 (pg/mL/day)
kIFNd  :   1.60 : IFN-gamma degradation (1/day)
kPRp   :  95.0 : proinflammatory production (pg/mL/day)
kPRd   :   0.95 : proinflammatory degradation (1/day)
dexPR  :   0.80 : additional proinflammatory clearance on dexamethasone
kILp   :  26.0 : IL-10 production (pg/mL/day)
kILd   :   0.85 : IL-10 degradation (1/day)
kWBCp  :   0.1330 : CSF leucocyte influx per unit PROIN (1/day)
KWCD4  :  55.0 : CD4 for half-maximal leucocyte trafficking (cells/uL)
kWBCd  :   0.20 : CSF leucocyte clearance (1/day)
dexWBC :   0.85 : fractional block of leucocyte trafficking by dexamethasone
klag   :   0.055 : lag of the immune-competence tracker (1/day)
kIRIS  :   2.00 : IRIS induction (per unit recovery rate per day)
KIR    : 200.0 : CSF GXM for half-maximal IRIS antigen drive (ug/mL)
kIRISoff:  0.085 : IRIS resolution (1/day)

// ---------------- antiretroviral therapy ----------------
kCD4   :   0.0062 : CD4 reconstitution rate constant (1/day)
CD4t   : 290.0 : CD4 target on suppressive ART (cells/uL)
kVL    :   0.34 : HIV RNA log10 decline on ART (1/day)

// ---------------- CSF hydrodynamics ----------------
Iform0 : 504.0 : CSF formation rate (mL/day = 0.35 mL/min)
fSupp  :   0.45 : maximum fractional suppression of formation by pressure
KICPf  : 260.0 : ICP for half-maximal suppression of formation (mmH2O)
Pss    :   8.00 : sagittal sinus pressure (cmH2O)
Pel0   :   2.80 : elastic pressure at zero excess volume (cmH2O)
Eel    :   0.0347 : elastance coefficient (1/mL)
Rout0  :   0.005556 : normal outflow resistance (cmH2O per mL/day)
aG     :  11.20 : GXM contribution to outflow resistance (multiplier)
KG     :  60.0 : CSF GXM for half-maximal obstruction (ug/mL)
aW     :   0.65 : CSF leucocyte contribution to outflow resistance
KW2    : 120.0 : CSF leucocytes for half-maximal contribution (cells/uL)
aI     :   1.55 : IRIS contribution to outflow resistance
kRon   :   0.55 : rate of rise of outflow resistance (1/day)
kRoff  :   0.10 : rate of resolution of outflow resistance (1/day)
kEon   :   1.35 : oedema formation (mL/day)
KE     : 340.0 : PROIN for half-maximal oedema (pg/mL)
KE2    :   0.55 : IRIS for half-maximal oedema
kEoff  :   0.17 : oedema resolution (1/day)
Vcsf   : 150.0 : total CSF volume (mL)

// ---------------- safety ----------------
GFR0   :  92.0 : baseline GFR (mL/min)
kGFRrec:   0.115 : GFR recovery rate (1/day)
fGFRtox:   0.62 : maximum fractional GFR loss from renal AmB
KGFRtox:   2.30 : renal cortical AmB for half-maximal toxicity (ug/mL)
Kser0  :   4.05 : baseline serum potassium (mmol/L)
kKrec  :   0.42 : potassium equilibration (1/day)
aKamb  :   1.35 : maximum potassium fall from renal AmB (mmol/L)
KKamb  :   2.10 : renal cortical AmB for half-maximal potassium loss (ug/mL)
Hb0    :   9.60 : baseline haemoglobin (g/dL)
kHbrec :   0.048 : haemoglobin recovery (1/day)
sHbD   :   0.1150 : haemoglobin fall per unit plasma AmB-d (g/dL per ug/mL/day)
sHbL   :   0.00105 : haemoglobin fall per unit plasma L-AmB
sHbF   :   0.130 : haemoglobin fall from 5FC marrow suppression
ANC0   :   3.10 : baseline neutrophil count (1e9/L)
kANCrec:   0.135 : neutrophil recovery (1/day)
sANC   :   2.50 : neutrophil fall from 5FC marrow suppression
KFCmar : 110.0 : plasma 5FC for half-maximal marrow suppression (ug/mL)
hFCmar :   1.80 : Hill coefficient, 5FC marrow suppression
ALT0   :  28.0 : baseline ALT (U/L)
kALTrec:   0.09 : ALT recovery (1/day)
sALTfl :   0.75 : ALT rise per unit plasma fluconazole
sALTfc :  22.0 : ALT rise from 5FC exposure

// ---------------- injury, disability and hazard ----------------
aNicp  :   0.0130 : injury from raised ICP
aNcpp  :   0.0180 : injury from low perfusion pressure
aNpro  :   0.0100 : injury from CSF inflammation
aNbur  :   0.0062 : injury from fungal burden
aNiris :   0.0290 : injury from IRIS
aNdex  :   0.0020 : injury from corticosteroid (myopathy/encephalopathy)
kNrep  :   0.0550 : repair of reversible injury (1/day)
kDIS   :   0.0300 : conversion of injury to permanent disability (1/day)
ANCthr :   1.50 : neutrophil hazard threshold (1e9/L)
Kthr   :   3.00 : potassium hazard threshold (mmol/L)
GFRthr :  60.0 : GFR hazard threshold (mL/min)
h0     :   0.000500 : background hazard of advanced HIV (1/day)
hB     :   0.013906 : hazard per log10 CFU/mL above 3 [FITTED]
hICP   :   0.030000 : hazard per unit of (ICP-250)/250 [fixed]
hCPP   :   0.030000 : hazard per unit of (60-CPP)/20 [fixed, unidentified]
hAN    :   0.000400 : hazard per g/dL fall in haemoglobin [fixed]
hNEUT  :   0.004000 : hazard per 1e9/L below ANCthr [fixed]
hK     :   0.001200 : hazard per mmol/L below Kthr [fixed]
hGFR   :   0.000400 : hazard per unit of (GFRthr-GFR)/30 [fixed]
hIRIS  :   0.008000 : hazard per unit IRIS activity [fixed]
hSTER  :   0.002900 : hazard per day the CSF culture remains positive [FITTED]
hAMS   :   0.004649 : hazard per unit NEUR x 10 [FITTED]
hDEX   :   0.000400 : steroid-attributable infection/metabolic hazard [fixed]

// ---------------- regimen switches (set per scenario) ----------------
AMBD_MGKG :  0.0 : AmB deoxycholate daily dose (mg/kg/day; 0 = none)
AMBD_DAYS :  0.0 : number of days of AmB deoxycholate
AMBD_TINF :  0.1667 : infusion duration (days; 4 h)
AMBL_MGKG :  0.0 : liposomal AmB single dose (mg/kg; 0 = none)
AMBL_TINF :  0.0833 : liposomal infusion duration (days; 2 h)
FC_MGKG   :  0.0 : flucytosine daily dose (mg/kg/day)
FC_DAYS   :  0.0 : days of flucytosine
FLU1_MG   :  0.0 : fluconazole induction dose (mg/day)
FLU1_T0   :  0.0 : fluconazole induction start (day)
FLU1_DAYS :  0.0 : days of fluconazole induction
FLU2_MG   :  0.0 : fluconazole consolidation dose (mg/day)
FLU2_T0   : 14.0 : fluconazole consolidation start (day)
ART_DAY   : 35.0 : day of ART initiation (>1e5 = never)
DEX_MGKG  :  0.0 : dexamethasone dose (mg/kg/day, tapered internally)
DEX_DAYS  :  0.0 : days of dexamethasone
SERT_MG   :  0.0 : sertraline dose (mg/day)
SERT_DAYS :  0.0 : days of sertraline
IFN_UG    :  0.0 : IFN-gamma dose per injection (ug)
IFN_N     :  0.0 : number of IFN-gamma injections (q48h from day 1)
KSUPP     :  0.0 : potassium replacement (mmol/L of target shift)
KSUPP_DAYS: 21.0 : days of potassium replacement
LP_N      :  0.0 : number of therapeutic lumbar punctures
LP_EVERY  :  2.0 : interval between lumbar punctures (days)
LP_T0     :  1.0 : day of the first therapeutic lumbar puncture
LP_WIN    :  0.0208 : duration of the drainage window (days; 30 min)
LP_PELTGT :  8.50 : target elastic pressure at end of drainage (cmH2O)
LP_GAIN   : 900.0 : drainage gain (mL/day per cmH2O of excess)
LP_FRAC   :   6.96 : first-order removal of CSF solutes during drainage (1/day)
CLAMPF    :  0.0 : 1 = hold the fungal states fixed (used by cm_burnin only)
LOGCFU0   :  5.0 : presenting CSF burden (log10 CFU/mL), used when INITMODE=1
CD40      : 25.0 : presenting CD4 count (cells/uL), used when INITMODE=1
INITMODE  :  1.0 : 1 = seed the unrelaxed presenting state in MAIN, 0 = use init()

$CMT @annotated
Ad    : AmB deoxycholate, central (mg)
Ad2   : AmB deoxycholate, peripheral (mg)
Al    : liposomal AmB, central (mg)
Al2   : liposomal AmB, peripheral (mg)
Abr   : AmB at the CNS effect site (mg)
Akid  : AmB in renal cortex (mg)
FCg   : flucytosine, gut depot (mg)
FCc   : flucytosine, central (mg)
FCcsf : flucytosine, CSF concentration (ug/mL)
FLg   : fluconazole, gut depot (mg)
FLc   : fluconazole, central (mg)
FLcsf : fluconazole, CSF concentration (ug/mL)
DXg   : dexamethasone, gut depot (mg)
DXc   : dexamethasone, central (mg)
SRg   : sertraline, gut depot (mg)
SRc   : sertraline, central (mg)
SRbr  : sertraline, brain (mg)
IFNsc : IFN-gamma, SC depot (ug)
Fe    : extracellular 5FC-susceptible yeast (CFU/mL)
Fres  : extracellular 5FC-resistant yeast (CFU/mL)
Ft    : phenotypically tolerant persisters (CFU/mL)
Fi    : intracellular yeast (CFU/mL equivalent)
Fp    : parenchymal reservoir (CFU/mL equivalent)
GXM   : CSF glucuronoxylomannan (ug/mL)
ERG   : membrane ergosterol, relative to untreated
CD4   : peripheral CD4 count (cells/uL)
VL    : HIV RNA (log10 copies/mL)
MAC   : activated CNS macrophage index
TH1   : CNS Th1 effector index
IFNG  : CSF interferon gamma (pg/mL)
PROIN : CSF proinflammatory index (pg/mL)
IL10  : CSF interleukin 10 (pg/mL)
WBC   : CSF leucocyte count (cells/uL)
IMML  : lagged immune competence
IRISa : IRIS activity
Vex   : excess CSF volume (mL)
Rout  : CSF outflow resistance (cmH2O per mL/day)
EDEMA : cerebral oedema, volume-equivalent (mL)
LEAK  : post-puncture dural leak flux (mL/day)
GFR   : glomerular filtration rate (mL/min)
Kser  : serum potassium (mmol/L)
Hb    : haemoglobin (g/dL)
ANC   : absolute neutrophil count (1e9/L)
ALT   : alanine aminotransferase (U/L)
NEUR  : reversible neuronal injury index
DIS   : permanent disability index
HAZ   : cumulative death hazard
AICP  : integral of (ICP - 250)+ (mmH2O.day)
ABR   : integral of CNS AmB concentration (ug.day/mL)
AFC   : integral of CSF 5FC concentration (ug.day/mL)
CLPV  : cumulative CSF volume drained (mL)
CLPG  : cumulative GXM removed (ug)
KILLC : cumulative yeast killed (CFU/mL equivalent)

$MAIN
// --- The presenting state is NOT hard-coded.  With INITMODE = 1 the model
// --- is seeded with an UNRELAXED state at the chosen burden and CD4 count;
// --- cm_burnin() then holds the fungal states fixed (CLAMPF = 1) for 70
// --- days while antigen, resistance, immunity, oedema and pressure relax
// --- to their quasi-steady values, and cm_run() starts every simulation
// --- from that relaxed state with INITMODE = 0.
// --- The 70-day clamp is a modelling choice, not an equilibrium: running
// --- it to convergence raises CSF GXM from 91.8 to 101.6 ug/mL and the
// --- opening pressure from 250 to 256 mmH2O.  That 6 mmH2O is how much the model says
// --- presenting pressure depends on how long the patient has been ill
// --- before reaching a hospital.
if (INITMODE > 0.5) {
  double F0 = pow(10.0, LOGCFU0);
  Fe_0    = 0.9930 * F0;
  Fres_0  = 1.0e-6 * F0;
  Ft_0    = 0.0020 * F0;
  Fi_0    = 0.0050 * F0;
  Fp_0    = 0.1500 * F0;
  GXM_0   = 40.0;
  ERG_0   = 1.0;
  CD4_0   = CD40;
  VL_0    = 5.4;
  MAC_0   = 1.0;
  TH1_0   = 0.05;
  IFNG_0  = 5.0;
  PROIN_0 = 40.0;
  IL10_0  = 10.0;
  WBC_0   = 2.0;
  IMML_0  = CD40 / (CD40 + KCD4);
  Vex_0   = 0.0;
  Rout_0  = Rout0;
  EDEMA_0 = 0.0;
  LEAK_0  = 0.0;
  GFR_0   = GFR0;
  Kser_0  = Kser0;
  Hb_0    = Hb0;
  ANC_0   = ANC0;
  ALT_0   = ALT0;
}

$ODE
// =====================================================================
// 0.  Exogenous input rates (piecewise-constant, regimen-driven)
// =====================================================================
double tt = SOLVERTIME;

// amphotericin B deoxycholate: once-daily 4-hour infusion
double rAmBd = 0.0;
if (AMBD_MGKG > 0.0 && tt < AMBD_DAYS) {
  double frac = tt - floor(tt);
  if (frac < AMBD_TINF) rAmBd = AMBD_MGKG * WT / AMBD_TINF;
}
// liposomal amphotericin B: single dose on day 0
double rAmBl = 0.0;
if (AMBL_MGKG > 0.0 && tt < AMBL_TINF) rAmBl = AMBL_MGKG * WT / AMBL_TINF;

// flucytosine, fluconazole, dexamethasone, sertraline: daily-average rates
double rFCg = (FC_MGKG > 0.0 && tt < FC_DAYS) ? FC_MGKG * WT : 0.0;
double rFLg = 0.0;
if (FLU1_MG > 0.0 && tt >= FLU1_T0 && tt < FLU1_T0 + FLU1_DAYS) rFLg = FLU1_MG;
else if (FLU2_MG > 0.0 && tt >= FLU2_T0) rFLg = FLU2_MG;

// CryptoDex taper: 0.3 mg/kg/day for one week, then 0.2, 0.1, 0.05, 0.03, 0.015
double rDXg = 0.0;
if (DEX_MGKG > 0.0 && tt < DEX_DAYS) {
  int wk = (int)floor(tt / 7.0);
  double sc = 1.0;
  if      (wk <= 0) sc = 1.0;
  else if (wk == 1) sc = 0.6667;
  else if (wk == 2) sc = 0.3333;
  else if (wk == 3) sc = 0.1667;
  else if (wk == 4) sc = 0.1000;
  else              sc = 0.0500;
  rDXg = DEX_MGKG * WT * sc;
}
double rSRg = (SERT_MG > 0.0 && tt < SERT_DAYS) ? SERT_MG :
              ((SERT_MG > 0.0) ? 0.5 * SERT_MG : 0.0);

// IFN-gamma: 100 ug SC every 48 h from day 1, given over 6 h
double rIFN = 0.0;
if (IFN_UG > 0.0) {
  for (int k = 0; k < (int)IFN_N; k++) {
    double td = 1.0 + 2.0 * k;
    if (tt >= td && tt < td + 0.25) rIFN += IFN_UG / 0.25;
  }
}
double ksupp = (KSUPP > 0.0 && tt < KSUPP_DAYS) ? KSUPP : 0.0;
double art   = (tt >= ART_DAY) ? 1.0 : 0.0;

// therapeutic lumbar puncture: proportional drain inside a 30-min window
int inLP = 0;
if (LP_N > 0.0) {
  for (int j = 0; j < (int)LP_N; j++) {
    double tlp = LP_T0 + LP_EVERY * j;
    if (tt >= tlp && tt < tlp + LP_WIN) inLP = 1;
  }
}

// =====================================================================
// 1.  Algebraic quantities.  ICP is a RESIDUAL, never a state variable.
// =====================================================================
double Cd   = Ad / VcD;
double Cd2  = Ad2 / VpD;
double Cl   = Al / VcL;
double Cl2  = Al2 / VpL;
double Cbr  = (Abr > 0.0) ? Abr / Vbr : 0.0;
double Ckid = Akid / Vkid;
double Cfc  = (FCc > 0.0) ? FCc / VFC : 0.0;
double Cfl  = (FLc > 0.0) ? FLc / VFL : 0.0;
double Cdx  = DXc / VDX;
double Csr  = SRc / VSR;
double CsrBrFree = fuSRbr * SRbr / Vbr;

double Vload  = Vex + EDEMA;
double Pel    = Pel0 * exp(Eel * Vload);
double ICPcm  = Pss + Pel;
double ICP    = 10.0 * ICPcm;                       // mmH2O
double ICPmmHg= ICP / 13.6;
double CPP    = MAP - ICPmmHg;
double Iform  = Iform0 * (1.0 - fSupp * ICP / (ICP + KICPf));
double Iabs   = Pel / Rout;

double CFU    = Fe + Fres + Ft + Fi;
if (CFU < 1e-9) CFU = 1e-9;
double logCFU = log10(CFU);

// --- drug effects.  The ergosterol node carries opposite signs. -------
double erg   = (ERG > 0.12) ? ERG : 0.12;
double cap   = 1.0 + kCap * GXM / (GXM + KCap);
double ec50a = EC50A * cap / pow(erg, aERG);
double killA = KmaxA * pow(Cbr, hA) / (pow(ec50a, hA) + pow(Cbr, hA));
double fcs   = (FCcsf > 0.0) ? FCcsf : 0.0;
double killF = KmaxF * fcs / (EC50F + fcs);
double killFr= KmaxF * fcs / (EC50F * RESFAC + fcs);
double killL = KmaxL * (1.0 - erg);
double gsupp = pow(erg, gERG);
double killSR= KmaxSR * CsrBrFree / (EC50SR + CsrBrFree);

double dex   = Cdx / (Cdx + KDX);
double killI = kImm * MAC * (1.0 + IFNG / (IFNG + KIFN)) * (1.0 - 0.75 * dex);
double immcomp = CD4 / (CD4 + KCD4);
double antigen = GXM / (GXM + KAg);
double fcMarrow= pow(Cfc, hFCmar) / (pow(Cfc, hFCmar) + pow(KFCmar, hFCmar));

// =====================================================================
// 2.  Amphotericin B
// =====================================================================
dxdt_Ad  = rAmBd - CLD * Cd - QD * (Cd - Cd2) - CLbrD * Cd - CLkdD * Cd;
dxdt_Ad2 = QD * (Cd - Cd2);
dxdt_Al  = rAmBl - CLL * Cl - QL * (Cl - Cl2) - CLbrL * Cl - CLkdL * Cl;
dxdt_Al2 = QL * (Cl - Cl2);
dxdt_Abr = CLbrD * Cd + CLbrL * Cl - kbrOut * Abr;
dxdt_Akid= CLkdD * Cd + CLkdL * Cl - kkdOut * Akid;

// =====================================================================
// 3.  Flucytosine -- renally cleared, so amphotericin nephrotoxicity
//     raises its exposure.  This is the clearance trap.
// =====================================================================
double clfc = CLFC0 * (0.12 + 0.88 * GFR / GFR0);
dxdt_FCg  = rFCg - kaFC * FCg;
dxdt_FCc  = FFC * kaFC * FCg - clfc * Cfc;
dxdt_FCcsf= keqFC * (fcCSF * Cfc - FCcsf);

// =====================================================================
// 4.  Fluconazole, dexamethasone, sertraline, interferon gamma
// =====================================================================
dxdt_FLg  = rFLg - kaFL * FLg;
dxdt_FLc  = FFL * kaFL * FLg - CLFL * Cfl;
dxdt_FLcsf= keqFL * (flCSF * Cfl - FLcsf);
dxdt_DXg  = rDXg - kaDX * DXg;
dxdt_DXc  = FDX * kaDX * DXg - CLDX * Cdx;
dxdt_SRg  = rSRg - kaSR * SRg;
dxdt_SRc  = FSR * kaSR * SRg - CLSR * Csr - kSRbr * Csr + kSRbrOut * SRbr;
dxdt_SRbr = kSRbr * Csr - kSRbrOut * SRbr;
dxdt_IFNsc= rIFN - kaIFN * IFNsc;

// =====================================================================
// 5.  Ergosterol -- one node, two opposite signs
// =====================================================================
double ergtar = 1.0 / (1.0 + FLcsf / IC50erg);
dxdt_ERG = kERG * (ergtar - ERG);

// =====================================================================
// 6.  Fungal compartments (clock 1)
// =====================================================================
double gr = g0 * gsupp * (1.0 - CFU / Fmax);
if (gr < 0.0) gr = 0.0;
double killE  = killA + killF + killL + killI + killSR;
double killR  = killA + killFr + killL + killI + killSR;
double killT  = tolPers * (killA + killF + killL + killSR) + killI;
double killIn = fIntra * killA + killF + kImmI * MAC * (1.0 - 0.75 * dex);
double killP  = fParen * (killA + killF + killL)
                + kImmP * MAC * (1.0 - 0.75 * dex);
double phag   = kPhag * MAC;
double lpout  = inLP ? LP_FRAC : 0.0;

dxdt_Fe   = gr * Fe - killE * Fe - phag * Fe + kEsc * Fi
            - kPer * Fe + kRev * Ft - kMut * Fe - kSeed * Fe + kShed * Fp
            - lpout * Fe;
dxdt_Fres = gr * Fres - killR * Fres - phag * Fres + kMut * Fe
            - kSeed * Fres - lpout * Fres;
dxdt_Ft   = 0.10 * gr * Ft + kPer * Fe - kRev * Ft - killT * Ft - lpout * Ft;
dxdt_Fi   = 0.30 * gr * Fi + phag * (Fe + Fres) - kEsc * Fi - killIn * Fi
            - lpout * Fi;
double gpp = gp * (1.0 - Fp / Fpmax);
if (gpp < 0.0) gpp = 0.0;
dxdt_Fp   = gpp * Fp + kSeed * (Fe + Fres) - kShed * Fp - killP * Fp;

double killFlux = killE * Fe + killR * Fres + killT * Ft + killIn * Fi
                  + killP * Fp;
// cm_burnin() sets CLAMPF = 1 to hold the fungal states at the presenting
// burden while the SLOW variables (antigen, resistance, immunity, oedema,
// pressure) relax to their quasi-steady values.  The killing FLUX is left
// intact, because it is a real source of capsular antigen even at
// equilibrium -- only the state derivatives are suppressed.
double cf = 1.0 - CLAMPF;
dxdt_Fe   *= cf;  dxdt_Fres *= cf;  dxdt_Ft *= cf;
dxdt_Fi   *= cf;  dxdt_Fp   *= cf;
dxdt_KILLC = killFlux;

// =====================================================================
// 7.  The antigen pool (clock 2).  Killing FEEDS it; no drug clears it.
// =====================================================================
double gsat = 1.0 / (1.0 + GXM / KGsat);
dxdt_GXM = (sLive * (Fe + Fres + Ft + Fi + Fp) + sLysis * killFlux) * gsat
           - kGXM * GXM - lpout * GXM;
dxdt_CLPG = lpout * GXM * Vcsf;

// =====================================================================
// 8.  Host immunity, HIV and ART
// =====================================================================
dxdt_CD4 = art * kCD4 * (CD4t - CD4);
dxdt_VL  = -art * kVL * VL;
dxdt_MAC = kMACon * antigen * (0.30 + 0.70 * immcomp) * (1.0 - 0.80 * dex)
           + kMACifn * IFNG / (IFNG + KIFN2) - kMACoff * (MAC - 1.0);
dxdt_TH1 = kTHon * antigen * immcomp * (1.0 - 0.90 * dex) - kTHoff * TH1;
dxdt_IFNG= kIFNp * TH1 + kIFNexo * FIFN * kaIFN * IFNsc - kIFNd * IFNG;
dxdt_PROIN = kPRp * (MAC * antigen + 3.0 * IRISa) - kPRd * PROIN
             - dexPR * dex * PROIN;
dxdt_IL10= kILp * MAC * antigen - kILd * IL10;
dxdt_WBC = kWBCp * PROIN * (CD4 / (CD4 + KWCD4)) * (1.0 - dexWBC * dex)
           - kWBCd * WBC - lpout * WBC;

// IRIS is the PRODUCT of a rate and a stock.  (immcomp - IMML) is
// proportional to d(immcomp)/dt, so no derivative of a state is needed.
dxdt_IMML = klag * (immcomp - IMML);
double recov = immcomp - IMML;
if (recov < 0.0) recov = 0.0;
dxdt_IRISa = kIRIS * recov * GXM / (GXM + KIR) * (1.0 - 0.70 * dex)
             - kIRISoff * IRISa;

// =====================================================================
// 9.  CSF hydrodynamics
// =====================================================================
double drain = 0.0;
if (inLP && Pel > LP_PELTGT) drain = LP_GAIN * (Pel - LP_PELTGT);
dxdt_Vex  = Iform - Iabs - LEAK - drain;
dxdt_CLPV = drain;
// post-puncture dural leak: the same total impulse (266 mL/day of
// flux, decaying with a 12-hour half-life) that the Python version adds
// instantaneously at each puncture.
dxdt_LEAK = inLP ? (266.0 / LP_WIN - 1.40 * LEAK) : (-1.40 * LEAK);

double Rtar = Rout0 * (1.0 + aG * GXM / (GXM + KG)
                       + aW * WBC / (WBC + KW2) + aI * IRISa);
double kR = (Rtar > Rout) ? kRon : kRoff;
dxdt_Rout = kR * (Rtar - Rout);
dxdt_EDEMA = kEon * (PROIN / (PROIN + KE) + 2.0 * IRISa / (IRISa + KE2))
             - kEoff * EDEMA;

// =====================================================================
// 10.  Safety
// =====================================================================
double gtar = GFR0 * (1.0 - fGFRtox * Ckid / (Ckid + KGFRtox));
dxdt_GFR = kGFRrec * (gtar - GFR);
double ktar = Kser0 - aKamb * Ckid / (Ckid + KKamb) + ksupp;
dxdt_Kser = kKrec * (ktar - Kser);
dxdt_Hb   = kHbrec * (Hb0 - Hb) - sHbD * Cd - sHbL * Cl - sHbF * fcMarrow;
// marrow suppression cannot remove precursors that are already gone
dxdt_ANC  = kANCrec * (ANC0 - ANC) - sANC * fcMarrow * ANC / (ANC + 0.30);
dxdt_ALT  = kALTrec * (ALT0 - ALT) + sALTfl * Cfl + sALTfc * fcMarrow;

// =====================================================================
// 11.  Injury, disability and the separable hazard
// =====================================================================
double icpX = (ICP > 250.0) ? (ICP - 250.0) / 100.0 : 0.0;
double cppX = (CPP < 60.0)  ? (60.0 - CPP) / 10.0 : 0.0;
double burX = (logCFU > 3.0) ? (logCFU - 3.0) : 0.0;

dxdt_DIS  = kDIS * NEUR * (1.0 - DIS);
dxdt_NEUR = (aNicp * icpX + aNcpp * cppX
             + aNpro * PROIN / (PROIN + 400.0)
             + aNbur * burX + aNiris * IRISa + aNdex * dex)
            * (1.0 - NEUR) - kNrep * NEUR;

double haz = h0
  + hB    * burX
  + hICP  * ((ICP > 250.0) ? (ICP - 250.0) / 250.0 : 0.0)
  + hCPP  * ((CPP < 60.0)  ? (60.0 - CPP) / 20.0 : 0.0)
  + hAN   * ((Hb < Hb0)    ? (Hb0 - Hb) : 0.0)
  + hNEUT * ((ANC < ANCthr)? (ANCthr - ANC) : 0.0)
  + hK    * ((Kser < Kthr) ? (Kthr - Kser) : 0.0)
  + hGFR  * ((GFR < GFRthr)? (GFRthr - GFR) / 30.0 : 0.0)
  + hIRIS * IRISa
  + ((CFU >= 1.0) ? hSTER : 0.0)
  + hDEX  * dex
  + hAMS  * NEUR * 10.0;
dxdt_HAZ = haz;

dxdt_AICP = (ICP > 250.0) ? (ICP - 250.0) : 0.0;
dxdt_ABR  = Cbr;
dxdt_AFC  = FCcsf;

$TABLE
double Vload_o = Vex + EDEMA;
double Pel_o   = Pel0 * exp(Eel * Vload_o);
double ICP_o   = 10.0 * (Pss + Pel_o);
double CFU_o   = Fe + Fres + Ft + Fi;
if (CFU_o < 1e-9) CFU_o = 1e-9;
double logCFU_o= log10(CFU_o);
double CPP_o   = MAP - ICP_o / 13.6;
double erg_o   = (ERG > 0.12) ? ERG : 0.12;
double cap_o   = 1.0 + kCap * GXM / (GXM + KCap);
double EC50Aeff= EC50A * cap_o / pow(erg_o, aERG);
double Cbr_o   = (Abr > 0.0) ? Abr / Vbr : 0.0;
double killA_o = KmaxA * pow(Cbr_o, hA) / (pow(EC50Aeff, hA) + pow(Cbr_o, hA));
double CrAg    = 160.0 * ((GXM > 1e-6) ? GXM : 1e-6);
double Mortality = 1.0 - exp(-HAZ);
double Rout_clin = Rout * 1440.0;
double Cfc_o   = (FCc > 0.0) ? FCc / VFC : 0.0;
double Cfl_o   = (FLc > 0.0) ? FLc / VFL : 0.0;
double SrFree  = fuSRbr * SRbr / Vbr;
double Sterile = (CFU_o < 1.0) ? 1.0 : 0.0;

$CAPTURE @annotated
ICP_o     : CSF pressure (mmH2O)
CPP_o     : cerebral perfusion pressure (mmHg)
logCFU_o  : CSF fungal burden (log10 CFU/mL)
CrAg      : CSF cryptococcal antigen titre (1:x)
EC50Aeff  : effective amphotericin EC50 (ug/mL)
killA_o   : amphotericin kill rate (1/day)
Cbr_o     : CNS amphotericin concentration (ug/mL)
Cfc_o     : plasma flucytosine (ug/mL)
Cfl_o     : plasma fluconazole (ug/mL)
SrFree    : free brain sertraline (ug/mL)
Rout_clin : CSF outflow resistance (cmH2O per mL/min)
Mortality : cumulative probability of death
Sterile   : CSF culture negative (1 = yes)
'

mod <- mcode("cm_qsp", cm_code, soloc = tempdir())

## =====================================================================
##  SCENARIOS -- 16 regimens spanning every randomised comparison
##  that has been made in this disease
## =====================================================================

WT <- 60

cm_scenarios <- list(

  ## ---- 1. natural history ------------------------------------------
  untreated = list(
    label = "No antifungal therapy (natural history)",
    par   = list(ART_DAY = 1e6)),

  ## ---- 2-3. fluconazole monotherapy (the resource-limited reality) --
  flu800 = list(
    label = "Fluconazole 800 mg/d x14d, then 400 mg/d",
    par   = list(FLU1_MG = 800, FLU1_DAYS = 14, FLU2_MG = 400)),

  flu1200 = list(
    label = "Fluconazole 1200 mg/d x14d, then 400 mg/d (Longley/Nussbaum)",
    par   = list(FLU1_MG = 1200, FLU1_DAYS = 14, FLU2_MG = 400)),

  ## ---- 4. ACTA all-oral arm ----------------------------------------
  oral_acta = list(
    label = "ACTA oral arm: fluconazole 1200 + flucytosine 100 mg/kg x14d",
    par   = list(FLU1_MG = 1200, FLU1_DAYS = 14, FLU2_MG = 800,
                 FC_MGKG = 100, FC_DAYS = 14)),

  ## ---- 5-8. Day 2013 / Bicanic 2008 amphotericin arms --------------
  ambd_alone = list(
    label = "AmB-d 1 mg/kg x28d, no partner (Day 2013 group 1)",
    par   = list(AMBD_MGKG = 1.0, AMBD_DAYS = 28, FLU2_MG = 400,
                 FLU2_T0 = 28, KSUPP = 0.55)),

  ambd_flu = list(
    label = "AmB-d 1 mg/kg + fluconazole 800 x14d (Day 2013 group 3)",
    par   = list(AMBD_MGKG = 1.0, AMBD_DAYS = 14, FLU1_MG = 800,
                 FLU1_DAYS = 14, FLU2_MG = 400, KSUPP = 0.55)),

  ambd_5fc07 = list(
    label = "AmB-d 0.7 mg/kg + 5FC x14d (Bicanic 2008 group 1)",
    par   = list(AMBD_MGKG = 0.7, AMBD_DAYS = 14, FC_MGKG = 100,
                 FC_DAYS = 14, FLU2_MG = 400, KSUPP = 0.55)),

  ambd_5fc = list(
    label = "AmB-d 1 mg/kg + 5FC x14d (Day 2013 group 2 / Bicanic group 2)",
    par   = list(AMBD_MGKG = 1.0, AMBD_DAYS = 14, FC_MGKG = 100,
                 FC_DAYS = 14, FLU2_MG = 400, KSUPP = 0.55)),

  ## ---- 9. ACTA / AMBITION control: one week of amphotericin --------
  ambd_5fc_1wk = list(
    label = "AmB-d 1 mg/kg + 5FC x7d, then FLU 1200 x7d (ACTA / AMBITION control)",
    par   = list(AMBD_MGKG = 1.0, AMBD_DAYS = 7, FC_MGKG = 100, FC_DAYS = 7,
                 FLU1_MG = 1200, FLU1_T0 = 7, FLU1_DAYS = 7,
                 FLU2_MG = 800, KSUPP = 0.55)),

  ## ---- 10. AMBITION: a single 10 mg/kg liposomal dose --------------
  ambition = list(
    label = "AMBITION: single L-AmB 10 mg/kg + 5FC + fluconazole 1200 x14d",
    par   = list(AMBL_MGKG = 10.0, FC_MGKG = 100, FC_DAYS = 14,
                 FLU1_MG = 1200, FLU1_DAYS = 14, FLU2_MG = 800,
                 KSUPP = 0.55)),

  ## ---- 11. add therapeutic lumbar punctures ------------------------
  ambition_lp = list(
    label = "AMBITION + 7 therapeutic lumbar punctures (days 1-13)",
    par   = list(AMBL_MGKG = 10.0, FC_MGKG = 100, FC_DAYS = 14,
                 FLU1_MG = 1200, FLU1_DAYS = 14, FLU2_MG = 800,
                 KSUPP = 0.55, LP_N = 7, LP_T0 = 1, LP_EVERY = 2)),

  ## ---- 12. CryptoDex ------------------------------------------------
  ambition_dex = list(
    label = "AMBITION + adjunctive dexamethasone (CryptoDex taper, 6 weeks)",
    par   = list(AMBL_MGKG = 10.0, FC_MGKG = 100, FC_DAYS = 14,
                 FLU1_MG = 1200, FLU1_DAYS = 14, FLU2_MG = 800,
                 KSUPP = 0.55, DEX_MGKG = 0.30, DEX_DAYS = 42)),

  ## ---- 13. ASTRO-CM -------------------------------------------------
  ambition_sert = list(
    label = "AMBITION + sertraline 400 mg/d (ASTRO-CM)",
    par   = list(AMBL_MGKG = 10.0, FC_MGKG = 100, FC_DAYS = 14,
                 FLU1_MG = 1200, FLU1_DAYS = 14, FLU2_MG = 800,
                 KSUPP = 0.55, SERT_MG = 400, SERT_DAYS = 14)),

  ## ---- 14. adjunctive interferon gamma ------------------------------
  ambition_ifng = list(
    label = "AMBITION + interferon gamma 100 ug SC x6 (Jarvis 2012)",
    par   = list(AMBL_MGKG = 10.0, FC_MGKG = 100, FC_DAYS = 14,
                 FLU1_MG = 1200, FLU1_DAYS = 14, FLU2_MG = 800,
                 KSUPP = 0.55, IFN_UG = 100, IFN_N = 6)),

  ## ---- 15-16. COAT: the ART-timing question -------------------------
  coat_early = list(
    label = "COAT backbone (AmB-d 0.9 + FLU 800 x14d), ART on day 10",
    par   = list(AMBD_MGKG = 0.9, AMBD_DAYS = 14, FLU1_MG = 800,
                 FLU1_DAYS = 14, FLU2_MG = 400, KSUPP = 0.55,
                 ART_DAY = 10)),

  coat_deferred = list(
    label = "COAT backbone (AmB-d 0.9 + FLU 800 x14d), ART on day 35",
    par   = list(AMBD_MGKG = 0.9, AMBD_DAYS = 14, FLU1_MG = 800,
                 FLU1_DAYS = 14, FLU2_MG = 400, KSUPP = 0.55,
                 ART_DAY = 35))
)

## =====================================================================
##  DRIVERS
## =====================================================================

#' Relax the SLOW variables to quasi-steady state with the fungal burden
#' held fixed at the presenting value, then release.  This is how the
#' presenting state is produced -- it is the model's own construction, not a
#' hand-set vector -- and it is the exact analogue of baseline() in
#' cm_reference_model.py.
cm_burnin <- function(logCFU = 5.0, CD4 = 25.0, burn = 70, mod. = mod) {
  s <- mod. %>%
    param(LOGCFU0 = logCFU, CD40 = CD4, INITMODE = 1, CLAMPF = 1,
          ART_DAY = 1e6, AMBD_MGKG = 0, AMBL_MGKG = 0, FC_MGKG = 0,
          FLU1_MG = 0, FLU2_MG = 0, DEX_MGKG = 0, SERT_MG = 0, IFN_UG = 0,
          LP_N = 0, KSUPP = 0) %>%
    mrgsim(end = burn, delta = burn, hmax = 0.02) %>%
    as_tibble()
  out <- s[nrow(s), ]
  ## the bookkeeping integrals are meaningless during the burn-in
  for (v in c("NEUR", "DIS", "HAZ", "AICP", "ABR", "AFC", "CLPV", "CLPG",
              "KILLC")) out[[v]] <- 0
  out
}

#' Run one scenario, starting from the relaxed presenting state.
cm_run <- function(name, days = 70, delta = 0.25, extra = list(),
                   logCFU = 5.0, CD4 = 25.0, mod. = mod) {
  sc <- cm_scenarios[[name]]
  if (is.null(sc)) stop("unknown scenario: ", name)
  y0   <- cm_burnin(logCFU = logCFU, CD4 = CD4, mod. = mod.)
  keep <- intersect(names(y0), names(init(mod.)))
  p    <- modifyList(sc$par, extra)
  p$INITMODE <- 0
  p$CLAMPF   <- 0
  out <- mod. %>%
    init(as.list(y0[keep])) %>%
    param(p) %>%
    mrgsim(end = days, delta = delta, hmax = 0.02, atol = 1e-10, rtol = 1e-8)
  as_tibble(out) %>% mutate(scenario = name, label = sc$label)
}

#' Run every scenario and return the long data frame.
cm_run_all <- function(days = 70, delta = 0.25) {
  bind_rows(lapply(names(cm_scenarios), cm_run, days = days, delta = delta))
}

#' Early fungicidal activity: least-squares slope of log10 CSF CFU/mL.
cm_efa <- function(sim, t0 = 0, t1 = 14) {
  d <- dplyr::filter(sim, time >= t0, time <= t1)
  unname(coef(lm(logCFU_o ~ time, data = d))[2])
}

#' Compact calibration table: the comparison this model exists to make.
cm_calibration <- function(days = 70) {
  targets <- tibble::tribble(
    ~scenario,       ~efa_pub, ~m10_pub, ~source,
    "flu800",          -0.07,    0.600,  "Longley 2008 CID",
    "flu1200",         -0.11,    0.560,  "Nussbaum 2010 CID",
    "oral_acta",       -0.28,    0.351,  "Nussbaum 2010 / ACTA 2018",
    "ambd_alone",      -0.31,    0.436,  "Day 2013 NEJM group 1",
    "ambd_flu",        -0.32,    0.326,  "Day 2013 NEJM group 3",
    "ambd_5fc07",      -0.40,    NA,     "Bicanic 2008 CID group 1",
    "ambd_5fc",        -0.49,    0.300,  "Day 2013 group 2 / Bicanic group 2",
    "ambd_5fc_1wk",    -0.42,    0.270,  "ACTA 2018 / AMBITION control",
    "ambition",        -0.40,    0.248,  "AMBITION 2022 NEJM")
  res <- lapply(targets$scenario, function(s) {
    sim <- cm_run(s, days = days)
    tibble::tibble(scenario = s,
                   efa_model = cm_efa(sim),
                   m2_model  = sim$Mortality[which.min(abs(sim$time - 14))],
                   m10_model = sim$Mortality[which.min(abs(sim$time - days))],
                   icp_peak  = max(sim$ICP_o),
                   sterile_day = {
                     i <- which(sim$Sterile == 1)
                     if (length(i)) sim$time[i[1]] else NA_real_
                   })
  })
  dplyr::left_join(dplyr::bind_rows(res), targets, by = "scenario")
}

## =====================================================================
##  WHAT THE MODEL SAYS THAT NO TRIAL ENDPOINT DOES
## ---------------------------------------------------------------------
##  1. TWO CLOCKS.  On the AMBITION regimen the CSF culture turns
##     negative on day 13.5 while CSF GXM is still 78 ug/mL -- 85% of
##     its presenting value.  Intracranial pressure stays above
##     250 mmH2O until day 19 and antigen is still 3.5 ug/mL at week 10.
##     EFA measures the clock that finishes first.
##
##  2. AN ANTIGEN SURGE CANNOT EXPLAIN POST-TREATMENT PRESSURE, and the
##     model refutes the hypothesis it was built to test.  The capsule
##     carried by every living yeast in the CSF is only ~50% of the
##     antigen already dissolved there, so complete sterilisation raises
##     GXM by at most ~55% (92 -> 132 ug/mL) and ICP by ~16 mmH2O.
##     Persistently raised pressure is the SLOW DECAY of a large standing
##     pool, not a burst from killing.
##
##  3. THE PARTNER-DRUG QUESTION IS AN ERGOSTEROL QUESTION.  Adding
##     fluconazole to amphotericin raises the effective amphotericin
##     EC50 from 3.16 to 5.05 ug/mL at day 7.  That single interaction
##     reproduces Day 2013's finding that fluconazole adds nothing to
##     amphotericin's fungicidal rate (-0.32 vs -0.31) while fluconazole
##     alone is clearly active (-0.13).  It recovers the DIRECTION of
##     ACTA's partner-drug result but only ~half its magnitude (model
##     hazard ratio 0.85 vs trial 0.62), so antagonism is a partial
##     explanation and the rest must lie in the reservoir dynamics.
##
##  4. THE SINGLE LIPOSOMAL DOSE IS AN EXPOSURE-MATCHING TRICK, NOT A
##     POTENCY ONE.  10 mg/kg once delivers a CNS AmB AUC0-14 of
##     3.4 ug.d/mL against 3.5 for seven daily 1 mg/kg doses -- and a
##     renal cortical AUC of 95 against 189.  Equal efficacy, half the
##     nephrotoxicity, and a LONGER time above the effect threshold
##     (20.2 vs 17.2 days) from one infusion.
##
##  5. AMPHOTERICIN NEPHROTOXICITY IS FLUCYTOSINE MYELOTOXICITY.  The
##     same 100 mg/kg/day of flucytosine gives a plasma peak of
##     33.5 ug/mL beside fluconazole (GFR 92) and 53.3 ug/mL beside two
##     weeks of amphotericin (GFR 42), with neutrophil nadirs of 1.68
##     and 0.90.  The interaction is a shared clearance organ, not
##     additive marrow toxicity.
##
##  6. ART TIMING IS AN ANTIGEN QUESTION, NOT A CD4 QUESTION.  IRIS
##     drive is the PRODUCT of the rate of immune recovery and the
##     antigen stock present when it starts.  Deferring ART from day 7
##     to day 35 lets GXM fall from 129 to 28 ug/mL and cuts 26-week
##     mortality from 61.6% to 45.6% -- a 16-point gap against COAT's
##     15 points, from the same mechanism, with no term for "ART timing".
##
##  7. SERTRALINE HAD TO FAIL.  At 400 mg/day the free brain
##     concentration peaks at 34 ng/mL, 0.56% of the in-vitro EC50.
##     The model predicts an EFA change of 0.000 log10/day.  ASTRO-CM
##     did not need to be run to know this; it needed the free
##     concentration to be compared with the MIC.
##
##  8. THE ONE THERAPY THAT TOUCHES CLOCK 2 IS A NEEDLE.  Seven
##     therapeutic lumbar punctures remove 140 mL of CSF and 11 mg of
##     GXM, hold the pressure-time integral above 250 mmH2O at 1
##     mmH2O.day instead of 173, and cut modelled 10-week mortality from
##     26.1% to 19.4% -- and in the high-burden phenotype from 70.0% to
##     49.4%.  This is the model's largest single intervention effect
##     and it involves no new drug.
## =====================================================================

## Print the calibration table only when this file is run on purpose:
##   Rscript cm_mrgsolve_model.R --calibrate
## Sourcing it (as cm_shiny_app.R does) stays silent and fast.
if ("--calibrate" %in% commandArgs(trailingOnly = TRUE)) {
  print(as.data.frame(cm_calibration()), digits = 4)
}
