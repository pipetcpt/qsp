# =============================================================================
#  mpm_mrgsolve_model.R
#  Malignant pleural mesothelioma (MPM) — quantitative systems pharmacology
#  51-state ODE model:  tumour geometry · one-sided delivery · pleural space ·
#  immunity · pemetrexed/platinum/bevacizumab/nivolumab/ipilimumab/ADI-PEG20 PK ·
#  marrow · kidney · biomarker · symptoms · survival
#
#  -------------------------------------------------------------------------
#  THE ORGANISING IDEA
#  -------------------------------------------------------------------------
#  Nearly every quantitative statement made about a solid tumour assumes the
#  tumour is a BALL: volume goes as the cube of a diameter, response is a
#  change in that diameter, and drug arrives from a vasculature that surrounds
#  the mass.  Pleural mesothelioma is a SHEET — a rind 0.3–3 cm thick coating
#  up to ~1300 cm2 of pleural surface.  Writing the burden as
#
#        V(t)  =  A(t)  x  h(t)          A = involved pleural area  [cm2]
#                                        h = mean rind thickness    [cm]
#
#  rather than V = (pi/6)d^3 changes four things, and each one is a specific
#  equation in this file rather than a remark:
#
#  (a) THE MEASUREMENT IS ONE FACTOR OF A PRODUCT.  Modified RECIST sums six
#      perpendicular THICKNESSES.  It sees h and is blind to A.  And h is not
#      the viable tumour: killed cells leave fibro-necrotic matrix that is
#      still measured, so
#          h_meas = (N + PHIM * M) / A
#      A drug that thins the rind while the sheet creeps circumferentially
#      scores as a response.  A sheet that spreads thins ITSELF.
#
#  (b) GROWTH IS A FRONT.  A sheet advances at its perimeter, P = 2*sqrt(pi*A),
#      and the proliferating pool is a shell of fixed depth L_OX, so it scales
#      with AREA while the burden scales with AREA x THICKNESS.  The volume
#      doubling time lengthens on its own, with no resistance mechanism.
#
#  (c) DELIVERY IS ONE-SIDED, AND THE FRONT IS THIN.  The rind is perfused
#      from the chest-wall/visceral face; the pleural-cavity face is
#      avascular.  The mean relative exposure across a slab of thickness h fed
#      from one side with penetration length L is
#          fpen(L,h) = (L/h) * (1 - exp(-h/L))
#      -> 1 for a thin rind, -> L/h for a thick one.  Two consequences follow
#      that are not obvious from a well-mixed model:
#         * a GEOMETRIC SANCTUARY of (1 - fpen) cells next to the pleural
#           cavity that no dose intensity reaches, which is why complete
#           response is essentially never seen in MPM;
#         * the ADVANCING MARGIN is only ~0.6 mm deep, so there fpen -> 1.
#           Chemotherapy arrests circumferential SPREAD long before it dents
#           the bulk — and it is spread, not bulk, that the hazard is most
#           sensitive to.
#      Intrapleural drug enters from the OPPOSITE face, so systemic and
#      intrapleural cover complementary depths rather than the same depth
#      twice.
#
#  (d) A SHEET CANNOT BE RESECTED WITH A MARGIN.  "Macroscopic complete
#      resection" collapses h and leaves A untouched, in a wound flooded with
#      IL-6 and TGF-beta, at the cost of permanent FVC and performance status.
#      MARS2 is an output of this model, not an input to it.
#
#  A fifth commitment runs underneath: THE PLEURAL SPACE HAS A VOLUME.
#  Effusion dilutes intrapleural drug, and the rind occludes the lymphatic
#  stomata that drain it, so effusion is NON-MONOTONE in tumour area — it
#  rises, peaks in mid-disease and the late hemithorax goes dry.
#
#  -------------------------------------------------------------------------
#  CALIBRATION
#  -------------------------------------------------------------------------
#  Control arms differ by more than four months across the MPM trials
#  (EMPHACIS cisplatin 9.3 mo; CheckMate 743 chemotherapy 14.1 mo; MARS2
#  chemotherapy 24.8 mo), so absolute medians are NOT a common scale.  The
#  model is anchored on three absolute values and otherwise calibrated to
#  HAZARD RATIOS, which are within-trial:
#
#     anchor  best supportive care ................ 7.6 mo   (MS01, ASC arm)
#     anchor  cisplatin alone ...................... 9.3 mo  (EMPHACIS control)
#     anchor  cisplatin + pemetrexed ............... 12.1 mo (EMPHACIS)
#     HR      cis+pem vs cis ....................... 0.77    (EMPHACIS)
#     HR      + bevacizumab ........................ 0.77    (MAPS)
#     HR      nivolumab+ipilimumab vs chemo ........ 0.74    (CheckMate 743)
#              epithelioid 0.86 / non-epithelioid 0.46
#     HR      pembrolizumab + chemo vs chemo ....... 0.79    (IND.227/KN-483)
#     HR      nivolumab 2L vs placebo .............. 0.69    (CONFIRM)
#     HR      pegargiminase + chemo, non-epithelioid 0.71    (ATOMIC-Meso)
#     HR      extended P/D + chemo vs chemo ........ 1.28    (MARS2)
#     HR      extrapleural pneumonectomy vs none ... 1.90    (MARS)
#  Response and toxicity anchors: mRECIST partial response ~41 % (EMPHACIS),
#  time to progression 5.7 vs 3.9 mo, grade 3/4 neutropenia falling from
#  ~24 % to ~5 % with folate/B12 supplementation.
#
#  -------------------------------------------------------------------------
#  VERIFICATION
#  -------------------------------------------------------------------------
#  There is no R runtime in the environment this model was built in.  Rather
#  than ship 51 unexecuted ODEs, the whole system was written twice: here in
#  mrgsolve C++ and again in `mpm_reference_model.py`, a dependency-free
#  Python RK4 implementation of exactly these equations.  Every number quoted
#  in README.md and in the repository table came from running the Python
#  version; `mpm_calibration_output.txt` is its verbatim output.  Running the
#  equations exposed five real defects, each fixed and each marked with a
#  "DEFECT #n" comment at the point of repair:
#
#    #1  no non-negativity floors: an intrapleural instillation drove N
#        marginally negative, h negative, fpen() negative, and the sign of
#        the whole rind flipped.
#    #2  the effusion equation had formation rising and drainage falling in
#        the same variable with no ceiling, and integrated to 1.3 million mL
#        in one year of untreated disease.  The pleural space has a volume.
#    #3  the metastatic compartment grew unbounded and then dominated the
#        hazard, so removing 89 % of the pleural tumour bought no survival.
#    #4  (Python driver only, not applicable here) a ten-minute pemetrexed
#        infusion is shorter than the integration step, so whether an RK4
#        stage landed inside the window decided how much drug was delivered;
#        V(400 d) differed by 44 % between dt = 0.02 and dt = 0.005 d and the
#        fitted EMAXP/EMAXC were compensating for missing drug.  mrgsolve
#        does not have this problem because it inserts exact dose records.
#    #5  symphysis did not block the intrapleural route, so talc APPEARED to
#        potentiate intrapleural chemotherapy by shrinking the diluting
#        volume -- the opposite of what the map asserts.
#
#  Two narrative claims were also withdrawn after the equations contradicted
#  them; both are recorded in mpm_calibration_output.txt sections H and L.
#
#  Rendering the companion map:
#     dot -Tsvg mpm_qsp_model.dot -o mpm_qsp_model.svg
# =============================================================================

library(mrgsolve)

mpm_code <- '
$PROB
# Malignant pleural mesothelioma QSP model
- 51 ODEs. Tumour written as AREA x THICKNESS, perfused from one face.
- Time unit: DAYS. Tumour tissue: cm3 (1 cm3 ~ 1e9 packed cells).
  Areas: cm2. Thickness: cm. Drug amounts: mg. Concentrations: mg/L.
  Effusion: mL. Plasma folate: nM. Homocysteine, arginine: umol/L.

$GLOBAL
#define CP_PEM   (PEM_C / V1_PEM)
#define CP_CIS   (CIS_F / V1_CIS)
#define C_BEV    (BEV_C / V1_BEV)
#define C_NIV    (NIV_C / V1_NIV)
#define C_IPI    (IPI_C / V1_IPI)
#define C_ADI    (ADI_C / V_ADI)

// Mean relative concentration across a slab of thickness h that is fed from
// ONE face with penetration length L.  -> 1 for h << L, -> L/h for h >> L.
// This single function carries the whole delivery-geometry argument.
double fpen(double L, double h) {
  if (h <= 1e-9) return 1.0;
  double r = h / L;
  if (r < 1e-6) return 1.0;
  if (r > 60.0) return L / h;         // exp(-r) has underflowed into noise
  return (L / h) * (1.0 - exp(-r));
}
double pos(double x) { return x > 0.0 ? x : 0.0; }
double clamp01(double x) { return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x); }

$PARAM @annotated
// ---- body ---------------------------------------------------------------
BSA: 1.8 : body surface area (m2)
WT: 70 : body weight (kg)
LBM0: 52 : baseline lean body mass (kg)

// ---- pleural geometry ----------------------------------------------------
S_PL: 1300 : total pleural surface of one hemithorax (cm2)
PHIM: 0.7 : fraction of a killed cells volume persisting as measured matrix
KDEGM: 0.008 : matrix resorption rate (1/d)

// ---- tumour cell kinetics -------------------------------------------------
KPROL: 0.268 : proliferation rate of a fully oxygenated cell (1/d)
L_OX: 0.018 : oxygen and nutrient penetration length (cm)
KDEATH0: 0.0005 : baseline apoptosis (1/d)
KNEC: 0.0018 : hypoxic death of the unreached fraction (1/d)
KFRONT: 0.055 : radial creep speed scale of the growth front (cm/d)
HFRONT: 0.06 : thickness OF THE ADVANCING MARGIN itself (cm)
GFRONT_REF: 0.07 : reference net proliferation at the margin (1/d)
ANG_MIN: 0.62 : residual front speed under complete VEGF blockade
KINV: 0.0035 : chest-wall and diaphragm invasion (cm/d)
KMET: 2.2e-05 : seeding of the nodal and distant compartment (1/d)
KMETG: 0.0035 : growth of the metastatic compartment (1/d)
METMAX: 150 : ceiling on the metastatic compartment (cm3)
WMET: 1.5 : weight of metastatic burden in the hazard

// ---- histology multipliers (epithelioid / biphasic / sarcomatoid) ---------
CHEMOS: 1 : chemosensitivity           1.00 / 0.62 / 0.35
IMMINF: 0.38 : immune infiltration        0.45 / 0.70 / 1.00
VISTA_S: 1.15 : VISTA-type suppression     1.00 / 0.65 / 0.35
COLLF: 1 : stromal collagen           1.00 / 1.45 / 1.90
MSLNF: 1 : mesothelin expression      1.00 / 0.50 / 0.12
SARCH: 1 : hazard multiplier          1.00 / 1.35 / 1.90
KPROLF: 1 : proliferation multiplier   1.00 / 1.25 / 1.55

// ---- pleural fluid ---------------------------------------------------------
KFORM: 1400 : maximal filtration scale (mL/d)
EMAXV: 3 : VEGF permeability multiplier
EC50V: 250 : VEGF for half-maximal permeability (pg/mL)
KDRAIN: 2200 : lymphatic stomatal capacity (mL/d)
FOBS: 0.95 : fraction of stomata a full rind can occlude
KMD: 300 : half-saturation of stomatal drainage (mL)
VEFFMAX: 3200 : pressure and anatomical ceiling of one hemithorax (mL)
KABS: 0.35 : absorption from an apposed (symphysed) space (1/d)
KSYMSP: 0.006 : spontaneous tumour symphysis (1/d)
KTALC: 0.35 : talc-induced symphysis (1/d)
VDEAD: 100 : residual pleural volume after drainage (mL)

// ---- VEGF, angiogenesis, delivery geometry -----------------------------------
KVSYN: 317 : VEGF synthesis (pg/mL/d per 100 cm3 viable tumour)
KVDEG: 1.5 : VEGF turnover (1/d)
KDBEV: 8 : bevacizumab-VEGF sequestration constant (mg/L)
KANG: 0.05 : angiogenesis rate (1/d)
EC50A: 300 : VEGF for half-maximal angiogenesis (pg/mL)
RHOMAX: 2.5 : maximal relative microvessel density
KVREG: 0.02 : vessel regression (1/d)
IFP_BASE: 4 : baseline interstitial fluid pressure (mmHg)
IFP_MAX: 22 : VEGF-driven interstitial pressure (mmHg)
EC50I: 350 : VEGF for half-maximal interstitial pressure (pg/mL)
IFP_REF: 17.9 : reference interstitial pressure for L_d scaling (mmHg)
LD0_CIS: 0.075 : cisplatin penetration length at the reference state (cm)
LD0_PEM: 0.06 : pemetrexed penetration length (cm)
LD0_MAB: 0.02 : 150 kDa IgG penetration length (cm)
L_T: 0.03 : T-cell infiltration length (cm)
LDIP0: 0.03 : intrapleural penetration from the free face (cm)
PEXP_RHO: 0.5 : L_d exponent on microvessel density
PEXP_IFP: 0.35 : L_d exponent on (IFP_ref / IFP)

// ---- pemetrexed PK and intracellular polyglutamation ---------------------------
V1_PEM: 12.9 : central volume (L)
V2_PEM: 8 : peripheral volume (L)
Q_PEM: 48 : intercompartmental clearance (L/d)
CL_PEM: 108 : clearance at CrCl 100 mL/min (L/d)
KIN_T: 6 : tumour uptake clearance (1/d)
KOUT_T: 20 : tumour efflux (1/d)
KFPGS_T: 3 : tumour folylpolyglutamate synthetase (1/d)
KGGH_T: 0.1 : tumour gamma-glutamyl hydrolase (1/d)
KIN_M: 6 : marrow uptake clearance (1/d)
KOUT_M: 20 : marrow efflux (1/d)
KFPGS_M: 0.9 : marrow FPGS (1/d)
KGGH_M: 0.38 : marrow GGH (1/d)
KI_FOL: 900 : weak RFC competition by circulating folate (nM)

// ---- folate and B12 --------------------------------------------------------------
FOL_DIET: 12 : unsupplemented plasma folate (nM)
KFOUT: 0.05 : folate turnover (1/d)
KFIN: 0.00413 : plasma folate per ug daily folic acid (nM/d/ug)
FOLUG    :  0.0   : daily oral folic acid dose (ug) - set by the regimen
B12      :  0.0   : vitamin B12 repletion flag (0/1) - set by the regimen
KHSYN: 26 : homocysteine production (umol/L/d)
KHDEG: 4.2 : folate-dependent remethylation (1/d)
KMH: 20 : folate for half-maximal remethylation (nM)
KHB12: 0.55 : extra remethylation with B12 repletion (1/d)

// ---- cisplatin PK -----------------------------------------------------------------
V1_CIS: 25 : free platinum volume (L)
KEL_CIS: 12 : renal elimination of free Pt (1/d)
KBIND: 21.3 : irreversible plasma protein binding (1/d)
VB_CIS: 5 : bound platinum volume (L)
KBEL: 0.14 : elimination of bound Pt (1/d)
KADD: 1 : adduct formation per mg/L per day
KREP: 0.35 : nucleotide excision repair (1/d)
ERCC1F: 1 : repair capacity multiplier

// ---- chemotherapy effect -------------------------------------------------------------
EMAXP: 0.1 : pemetrexed maximal kill of exposed cells (1/d)
EC50P: 0.9 : polyglutamate for half-maximal kill (mg/L)
KFPT: 400 : folate rescue constant, TUMOUR - weak (nM)
EMAXC: 0.78 : cisplatin maximal kill of exposed cells (1/d)
EC50C: 0.28 : adducts for half-maximal kill
SYNCP: 0.9 : pemetrexed x platinum synergy
EMAXVIN: 0.11 : vinorelbine / gemcitabine second line (1/d)
EC50VIN: 0.5 : vinorelbine effect-site EC50 (mg/L)
RTRATE   :  0.0   : hemithoracic radiotherapy kill (1/d) - set by the regimen
PRMT5    :  0.0   : PRMT5 inhibitor kill in MTAP-deleted disease (1/d)
TALCON   :  0.0   : talc pleurodesis on/off - set by the regimen
DRAINML  :  0.0   : therapeutic drainage (mL/d) - set by the regimen

// ---- marrow (Friberg transit) ----------------------------------------------------------
CIRC0: 4 : baseline absolute neutrophil count (10^9/L)
KTR: 0.873 : transit rate, MTT 110 h over 4 compartments (1/d)
GAM: 0.17 : rebound feedback exponent
EMAXPM: 1.38 : pemetrexed marrow effect
EC50PM: 0.85 : marrow polyglutamate EC50
KFPM: 14 : folate rescue constant, MARROW - strong (nM)
SLOPE_CIS: 0.09 : marrow effect per mg/L free Pt

// ---- renal -------------------------------------------------------------------------------
CRCL0: 95 : baseline creatinine clearance (mL/min)
KRECR: 0.05 : recovery toward the falling set point (1/d)
KNEPH: 47 : reversible nephrotoxicity per mg/L/d free Pt
KNEPHIRR: 12 : permanent nephrotoxicity per mg/L/d free Pt

// ---- biologics PK ---------------------------------------------------------------------------
V1_BEV: 3 : bevacizumab central volume (L)
V2_BEV: 2.4 : bevacizumab peripheral volume (L)
Q_BEV: 0.55 : bevacizumab intercompartmental clearance (L/d)
CL_BEV: 0.21 : bevacizumab clearance (L/d)
V1_NIV: 4.5 : nivolumab central volume (L)
V2_NIV: 3.2 : nivolumab peripheral volume (L)
Q_NIV: 0.65 : nivolumab intercompartmental clearance (L/d)
CL_NIV: 0.24 : nivolumab clearance (L/d)
V1_IPI: 4.4 : ipilimumab volume (L)
CL_IPI: 0.37 : ipilimumab clearance (L/d)
EC50_PD1: 0.1 : nivolumab concentration for half PD-1 occupancy (mg/L)
EC50_CT4: 1.5 : ipilimumab concentration for half CTLA-4 occupancy (mg/L)

// ---- ADI-PEG20 and arginine ----------------------------------------------------------------
V_ADI: 4 : ADI-PEG20 volume (L)
KA_ADI: 0.5 : intramuscular absorption (1/d)
CL_ADI: 0.3 : ADI-PEG20 clearance (L/d)
KADA: 0.0045 : anti-drug antibody formation (1/d)
KADAEL: 0.025 : anti-drug antibody elimination (1/d)
ADA_POT: 8 : ADA-driven clearance multiplier
KARGIN: 45 : arginine appearance (umol/L/d)
KARGOUT: 0.5 : arginine disposal (1/d)
KADIARG: 20 : arginine deimination per mg/L ADI (1/d)
ARG_CRIT: 12 : arginine below which ASS1-null cells starve (umol/L)
EMAXARG: 0.045 : kill at full arginine depletion (1/d)
ASS1LOSS: 1 : 1 = ASS1 methylated (auxotroph), 0 = ASS1 intact

// ---- immunity ---------------------------------------------------------------------------------
KPRIME: 0.035 : priming rate (1/d)
EIPI: 4.5 : ipilimumab amplification of priming
KEXH: 0.09 : exhaustion under PD-1 engagement (1/d)
PDL1B: 1 : PD-L1 scaling
KTD: 0.035 : effector decay (1/d)
KSUP: 0.018 : suppression per unit Treg (1/d)
KTREGS: 0.04 : Treg induction by TGF-beta (1/d)
KMTG: 8 : TGF-beta for half-maximal Treg induction
KTREGD: 0.03 : Treg decay (1/d)
KIPIDEP: 0.045 : ipilimumab-mediated Treg depletion (1/d)
KCLON: 0.05 : memory formation (1/d)
KCLOND: 0.01 : memory decay - the CheckMate 743 tail (1/d)
WMEM: 0.6 : weight of memory pool in the effector function
KKILL: 0.33 : maximal immune kill of exposed cells (1/d)
KMKILL: 6 : effector pool for half-maximal kill
KTGF: 0.045 : TGF-beta production per 100 cm3 tumour (1/d)
KTGFD: 0.12 : TGF-beta turnover (1/d)

// ---- inflammation, cachexia, symptoms -------------------------------------------------------------
KIL6T: 41.4 : IL-6 per 100 cm3 tumour (pg/mL/d)
KIL6W: 260 : IL-6 from the surgical wound (pg/mL/d)
KIL6D: 8 : IL-6 turnover (1/d)
KWD: 0.06 : resolution of the surgical wound response (1/d)
KCACH: 0.013 : lean mass loss at saturating IL-6 (kg/d)
KMIL6: 22 : IL-6 for half-maximal cachexia (pg/mL)
KLREC: 0.0025 : lean mass recovery (1/d)
FVC0: 92 : predicted FVC at diagnosis (%)
A_FVC: 0.45 : FVC lost to full pleural encasement
H_FVC: 0.1 : FVC lost to a 3 cm rind
V_FVC: 0.22 : FVC lost to a 2.5 L effusion

// ---- biomarker -------------------------------------------------------------------------------------
KSMRP: 0.0217 : mesothelin release (nM/d per cm3 viable tumour)
KSMRPCL: 1.4 : mesothelin renal clearance at CrCl 95 (1/d)

// ---- hazard -----------------------------------------------------------------------------------------
H0: 0.000167 : baseline hazard (1/d)
V0H: 360 : reference burden (cm3)
B1: 0.5 : burden exponent
B2: 0.55 : hazard per ECOG point
B5: 2.2 : hazard per unit of pleural surface ENCASED
H_FEBN: 0.004 : added hazard while ANC < 0.5 (1/d)
H_IRAE: 0.00095 : added hazard per unit grade >=3 irAE intensity (1/d)

// ---- intrapleural route ---------------------------------------------------------------------------------
KIPCL: 8 : clearance of drug out of the pleural space (1/d)
FABS_IP: 0.55 : fraction of that clearance reaching the systemic pool

// ---- initial condition controls ---------------------------------------------------------------------------
H0_RIND  :  0.80  : rind thickness at baseline (cm)
A0_RIND  : 450.0  : involved pleural area at baseline (cm2)
MFRAC0   :  0.28  : fraction of baseline measured volume that is matrix

$CMT @annotated
N       : viable tumour tissue (cm3)
M       : necrotic and fibrotic matrix (cm3)
A       : involved pleural area (cm2)
Z       : chest-wall and diaphragm invasion depth (cm)
MET     : nodal and distant burden (cm3)
VEFF    : pleural effusion volume (mL)
PSY     : pleural symphysis (0-1)
VEGF    : VEGF-A (pg/mL)
RHOV    : relative microvessel density
IL6     : interleukin-6 (pg/mL)
TGFB    : TGF-beta (a.u.)
TEFF    : intratumoural CD8 effectors (a.u.)
TREG    : regulatory T cells (a.u.)
TCLON   : memory / clonal T-cell pool (a.u.)
WOUND   : surgical wound signal (0-1)
PEM_C   : pemetrexed central (mg)
PEM_P   : pemetrexed peripheral (mg)
PEM_T   : tumour intracellular pemetrexed (mg/L)
PEM_TP  : tumour pemetrexed polyglutamates (mg/L)
PEM_M   : marrow intracellular pemetrexed (mg/L)
PEM_MP  : marrow pemetrexed polyglutamates (mg/L)
FOL     : plasma folate (nM)
HCY     : homocysteine (umol/L)
CIS_F   : free platinum central (mg)
CIS_B   : protein-bound platinum (mg/L)
ADD     : Pt-DNA adducts in systemically exposed cells
ADDIP   : Pt-DNA adducts from the intrapleural route
BEV_C   : bevacizumab central (mg)
BEV_P   : bevacizumab peripheral (mg)
NIV_C   : nivolumab central (mg)
NIV_P   : nivolumab peripheral (mg)
IPI_C   : ipilimumab central (mg)
ADI_A   : ADI-PEG20 intramuscular depot (mg)
ADI_C   : ADI-PEG20 central (mg)
ADA     : anti-ADI antibody (a.u.)
ARG     : plasma arginine (umol/L)
IPD     : drug amount in the pleural space (mg)
PROL    : marrow proliferating pool (10^9/L)
TR1     : maturation transit 1 (10^9/L)
TR2     : maturation transit 2 (10^9/L)
TR3     : maturation transit 3 (10^9/L)
CIRC    : circulating absolute neutrophil count (10^9/L)
CRCL    : measured creatinine clearance (mL/min)
CRCLSS  : irreversibly lost creatinine clearance set point (mL/min)
SMRP    : serum soluble mesothelin (nM)
LBM     : lean body mass (kg)
CUMH    : cumulative hazard
AUCP    : cumulative pemetrexed AUC (mg.d/L)
AUCC    : cumulative free platinum AUC (mg.d/L)
VINE    : vinorelbine effect site (mg/L)
IRAE    : immune-related adverse event intensity (0-1)

$MAIN
double Vtot0 = H0_RIND * A0_RIND;
double M0    = MFRAC0 * Vtot0 / PHIM;
double N0    = Vtot0 - PHIM * M0;

N_0      = N0;
M_0      = M0;
A_0      = A0_RIND;
Z_0      = 0.18;
MET_0    = 2.0;
VEFF_0   = 1000.0;
PSY_0    = 0.0;
VEGF_0   = KVSYN * (N0 / 100.0) * 0.98 / KVDEG;
RHOV_0   = 1.0;
IL6_0    = KIL6T * (N0 / 100.0) / KIL6D;
TGFB_0   = KTGF * (N0 / 100.0) / KTGFD;

// immune compartments start AT their untreated steady state, so an arm
// without immunotherapy shows no spurious opening transient
double tgf0  = TGFB_0;
double treg0 = KTREGS * tgf0 / (tgf0 + KMTG) / KTREGD;
double antg0 = pow(N0 / 290.0, 0.30) * IMMINF;
double prim0 = KPRIME * antg0 / (1.0 + VISTA_S);
double teff0 = prim0 / (KEXH * PDL1B + KTD + KSUP * treg0);
TREG_0   = treg0;
TEFF_0   = teff0;
TCLON_0  = KCLON * teff0 / KCLOND;

FOL_0    = FOL_DIET;
HCY_0    = KHSYN / (KHDEG * (FOL_DIET / (FOL_DIET + KMH)));
ARG_0    = KARGIN / KARGOUT;
PROL_0   = CIRC0;  TR1_0 = CIRC0;  TR2_0 = CIRC0;  TR3_0 = CIRC0;  CIRC_0 = CIRC0;
CRCL_0   = CRCL0;
CRCLSS_0 = CRCL0;
SMRP_0   = KSMRP * N0 * MSLNF / KSMRPCL;
LBM_0    = LBM0;

$ODE
// ---------------------------------------------------------------------------
// DEFECT #1 (exposed by the Python re-implementation).  Without floors, an
// intrapleural instillation drives N marginally negative within the first
// hours (C_IP is ~200 mg/L for a few minutes).  N < 0 makes h < 0, makes
// fpen() return a negative exposed fraction, and the sign of the whole rind
// flips.  Every state read below is guarded.
// ---------------------------------------------------------------------------
double nN   = pos(N);
double nM   = pos(M);
double nA   = A < 1.0 ? 1.0 : (A > S_PL ? S_PL : A);
double psy  = clamp01(PSY);
double vegf = pos(VEGF);
double rho  = RHOV < 1e-3 ? 1e-3 : RHOV;
double crcl = CRCL < 8.0 ? 8.0 : CRCL;
double anc  = CIRC < 0.01 ? 0.01 : CIRC;

// ---- 1. geometry -----------------------------------------------------------
double h_tot = (nN + PHIM * nM) / nA;      // what the CT scanner measures
double h_via = nN / nA;                    // what the drug has to kill
double Vtum  = nN + PHIM * nM;
double afrac = nA / S_PL;

// ---- 2. delivery geometry --------------------------------------------------
double VEGFf = vegf / (1.0 + C_BEV / KDBEV);
double IFP   = IFP_BASE + IFP_MAX * VEGFf / (VEGFf + EC50I);
double fv    = pow(rho, PEXP_RHO);
double fi    = pow(IFP_REF / (IFP < 1e-3 ? 1e-3 : IFP), PEXP_IFP);
double fc    = 1.0 / sqrt(COLLF);
double Ld_cis = LD0_CIS * fv * fi * fc;
double Ld_pem = LD0_PEM * fv * fi * fc;
double Ld_mab = LD0_MAB * fv * fi * fc;
double Ld_tc  = L_T     * fv * fi * fc;
double Ld_ip  = LDIP0   * fc;              // needs no vasculature

double f_ox  = fpen(L_OX,  h_tot);
double f_cis = fpen(Ld_cis, h_tot);
double f_pem = fpen(Ld_pem, h_tot);
double f_tc  = fpen(Ld_tc,  h_tot);
double f_ip  = fpen(Ld_ip,  h_tot);
double f_sys = f_cis > f_pem ? f_cis : f_pem;
// ---------------------------------------------------------------------------
// DEFECT #5.  The mechanistic map says pleurodesis CLOSES the intrapleural
// route, and the first version of these equations did not implement it:
// symphysis entered only through the effusion, so obliterating the space made
// the instillate MORE concentrated (a smaller diluting volume) and talc came
// out looking like a way to potentiate intrapleural chemotherapy.  Two factors
// were missing -- an apposed surface cannot be instilled into at all, and even
// a free space is loculated, so an instillate bathes only part of the involved
// surface (FCONT).
// ---------------------------------------------------------------------------
double f_ipa = f_ip * FCONT * (1.0 - psy);
// systemic and intrapleural fronts advance from OPPOSITE faces, so the
// intrapleural route only adds coverage where the systemic route stops
double f_ipo = pos(f_ipa < (1.0 - f_sys) ? f_ipa : (1.0 - f_sys));

double C_IP  = IPD / ((VEFF + VDEAD) / 1000.0 < 0.02 ? 0.02 : (VEFF + VDEAD) / 1000.0);

// ---- 3. pemetrexed ---------------------------------------------------------
double CLpem = CL_PEM * (crcl / 100.0);
dxdt_PEM_C = - CLpem * CP_PEM - Q_PEM * CP_PEM + Q_PEM * PEM_P / V2_PEM;
dxdt_PEM_P =   Q_PEM * CP_PEM - Q_PEM * PEM_P / V2_PEM;
dxdt_AUCP  =   CP_PEM;

double upt = CP_PEM / (1.0 + FOL / KI_FOL);
dxdt_PEM_T  = KIN_T * upt - KOUT_T * PEM_T - KFPGS_T * PEM_T + KGGH_T * PEM_TP;
dxdt_PEM_TP = KFPGS_T * PEM_T - KGGH_T * PEM_TP;
dxdt_PEM_M  = KIN_M * upt - KOUT_M * PEM_M - KFPGS_M * PEM_M + KGGH_M * PEM_MP;
dxdt_PEM_MP = KFPGS_M * PEM_M - KGGH_M * PEM_MP;

dxdt_FOL = KFIN * FOLUG + KFOUT * FOL_DIET - KFOUT * FOL;
dxdt_HCY = KHSYN - KHDEG * (FOL / (FOL + KMH)) * HCY - KHB12 * B12 * HCY;

// ---- 4. cisplatin ----------------------------------------------------------
dxdt_CIS_F = FABS_IP * KIPCL * IPD - (KEL_CIS + KBIND) * CIS_F;
dxdt_CIS_B = KBIND * CIS_F / VB_CIS - KBEL * CIS_B;
dxdt_AUCC  = CP_CIS;
dxdt_ADD   = KADD * CP_CIS - KREP * ERCC1F * ADD;
dxdt_ADDIP = KADD * C_IP   - KREP * ERCC1F * ADDIP;
dxdt_IPD   = - KIPCL * IPD;

// ---- 5. biologics ----------------------------------------------------------
dxdt_BEV_C = - CL_BEV * C_BEV - Q_BEV * C_BEV + Q_BEV * BEV_P / V2_BEV;
dxdt_BEV_P =   Q_BEV * C_BEV - Q_BEV * BEV_P / V2_BEV;
dxdt_NIV_C = - CL_NIV * C_NIV - Q_NIV * C_NIV + Q_NIV * NIV_P / V2_NIV;
dxdt_NIV_P =   Q_NIV * C_NIV - Q_NIV * NIV_P / V2_NIV;
dxdt_IPI_C = - CL_IPI * C_IPI;
double RO_PD1 = C_NIV / (C_NIV + EC50_PD1);
double RO_CT4 = C_IPI / (C_IPI + EC50_CT4);

double cladi = CL_ADI * (1.0 + ADA_POT * ADA);
dxdt_ADI_A = - KA_ADI * ADI_A;
dxdt_ADI_C =   KA_ADI * ADI_A - cladi * C_ADI;
dxdt_ADA   =   KADA * (C_ADI > 0.02 ? 1.0 : 0.0) - KADAEL * ADA;
dxdt_ARG   =   KARGIN - KARGOUT * ARG - KADIARG * C_ADI * ARG / (ARG + 20.0);

dxdt_VINE  = - 0.60 * VINE;

// ---- 6. VEGF, vessels, pleural space ---------------------------------------
dxdt_VEGF = KVSYN * (nN / 100.0) * (0.30 + 0.70 * (1.0 - f_ox)) - KVDEG * vegf;
dxdt_RHOV = KANG * (VEGFf / (VEGFf + EC50A)) * (1.0 - rho / RHOMAX) - KVREG * rho;

// ---------------------------------------------------------------------------
// DEFECT #2.  Written without a ceiling this is Jform(A) - Jdrain(A) with
// formation RISING and stomatal drainage FALLING in the same variable, and it
// integrated to 1.3 million mL in one year of untreated disease.  The missing
// physics is that the pleural space HAS A VOLUME: filtration is a Starling
// flux that stops as pleural pressure rises, and the space available is what
// the hemithorax holds MINUS the volume the rind itself occupies, minus
// fibrothorax contraction.  With the ceiling in place, effusion becomes
// non-monotone in A without any extra rule.
// Symphysis blocks FORMATION, not absorption: an apposed space still resorbs.
// ---------------------------------------------------------------------------
double vcap  = (VEFFMAX - 0.90 * Vtum) * (1.0 - psy) * (1.0 - 0.45 * afrac);
if (vcap < 60.0) vcap = 60.0;
double Jform = KFORM * (0.05 + afrac) * (1.0 + EMAXV * VEGFf / (VEGFf + EC50V))
               * (1.0 - psy) * pos(1.0 - VEFF / vcap);
double Jdrain = KDRAIN * (1.0 - FOBS * afrac) * VEFF / (VEFF + KMD)
                + KABS * psy * VEFF;
dxdt_VEFF = Jform - Jdrain - DRAINML;
double trapped = (afrac * afrac) * (h_tot / (h_tot + 0.6)) * 2.2;
if (trapped > 1.0) trapped = 1.0;
dxdt_PSY = (KSYMSP * afrac * afrac * afrac + KTALC * TALCON)
           * (1.0 - psy) * (1.0 - 0.85 * trapped);

// ---- 7. immunity ------------------------------------------------------------
double antig = pow(nN < 1.0 ? 1.0 : nN / 290.0, 0.30) * IMMINF;
double prime = KPRIME * antig * (1.0 + EIPI * RO_CT4) / (1.0 + VISTA_S);
dxdt_TEFF  = prime - KEXH * PDL1B * (1.0 - RO_PD1) * TEFF - KTD * TEFF
             - KSUP * TREG * TEFF;
dxdt_TREG  = KTREGS * TGFB / (TGFB + KMTG) - KTREGD * TREG - KIPIDEP * RO_CT4 * TREG;
dxdt_TCLON = KCLON * TEFF - KCLOND * TCLON;
dxdt_TGFB  = KTGF * (nN / 100.0) + 3.0 * WOUND - KTGFD * TGFB;
double Tpool = TEFF + WMEM * TCLON;
double k_imm = KKILL * (Tpool / (Tpool + KMKILL)) * f_tc;

// ---- 8. drug effects ----------------------------------------------------------
double E_pem = EMAXP * PEM_TP / (PEM_TP + EC50P * (1.0 + FOL / KFPT));
double E_cis = EMAXC * ADD    / (ADD    + EC50C);
double E_vin = EMAXVIN * VINE / (VINE + EC50VIN);
double E_sys = CHEMOS * (E_pem + E_cis + SYNCP * E_pem * E_cis + E_vin);
double E_ipl = CHEMOS * EMAXC * ADDIP / (ADDIP + EC50C);
double E_arg = EMAXARG * ASS1LOSS * pos(ARG_CRIT - ARG) / ARG_CRIT;

// Radiotherapy and PRMT5 inhibition are not diffusion-limited in the same
// way: photons reach the whole rind, and an oral small molecule given
// continuously equilibrates, so both act on the whole viable pool.
double kill = E_sys * f_sys + E_ipl * f_ipo + k_imm + E_arg * f_cis
              + RTRATE + PRMT5;

// ---- 9. tumour geometry ODEs -----------------------------------------------------
double growth = KPROL * KPROLF * f_ox * nN;
double deaths = KDEATH0 * nN + KNEC * (1.0 - f_ox) * nN + kill * nN;
dxdt_N = growth - deaths;
dxdt_M = PHIM * deaths - KDEGM * nM;

// ---------------------------------------------------------------------------
// THE ADVANCING MARGIN IS THIN.  This is the sharpest consequence of writing
// the tumour as a sheet.  At the growth front the rind is a few hundred
// microns deep, so THERE fpen -> 1: the front is fully drug-exposed and fully
// oxygenated even when the bulk behind it is neither.  Chemotherapy arrests
// circumferential SPREAD long before it makes a dent in the bulk -- and it is
// spread (encasement) that the hazard is most sensitive to.  Evaluating the
// same fpen() at h = HFRONT rather than at h = h_tot is the whole difference.
// ---------------------------------------------------------------------------
double perim = 2.0 * sqrt(M_PI * nA);
double fox_f = fpen(L_OX, HFRONT);
double fcp_f = fpen(Ld_cis, HFRONT);
double fpp_f = fpen(Ld_pem, HFRONT);
double fch_f = fcp_f > fpp_f ? fcp_f : fpp_f;  if (fch_f > 1.0) fch_f = 1.0;
double fip_f = fpen(Ld_ip, HFRONT);            if (fip_f > 1.0) fip_f = 1.0;
fip_f *= FCONT * (1.0 - psy);
double ftc_f = fpen(Ld_tc, HFRONT);            if (ftc_f > 1.0) ftc_f = 1.0;
double angio_f = ANG_MIN + (1.0 - ANG_MIN) * VEGFf / (VEGFf + EC50A);
double kimm_f  = KKILL * (Tpool / (Tpool + KMKILL)) * ftc_f;
double front_net = KPROL * KPROLF * fox_f * angio_f
                   - KDEATH0 - KNEC * (1.0 - fox_f)
                   - E_sys * fch_f - E_ipl * fip_f - kimm_f
                   - E_arg * fch_f - RTRATE - PRMT5;
double vfront = front_net / GFRONT_REF;
if (vfront < 0.0) vfront = 0.0;
if (vfront > 1.5) vfront = 1.5;
dxdt_A = KFRONT * perim * (1.0 - nA / S_PL) * vfront;

double viab = nN / (Vtum < 1e-6 ? 1e-6 : Vtum);
dxdt_Z = KINV * (h_via / (h_via + 0.30)) * viab
         - 0.0020 * Z * (kill / 0.02 > 1.0 ? 1.0 : kill / 0.02);

// ---------------------------------------------------------------------------
// DEFECT #3.  Written as unbounded exponential growth with only a
// kill-dependent brake, MET reached 660 cm3 by two years in the
// immunotherapy arm and then DOMINATED the hazard, so a scenario that removed
// 89 % of the pleural tumour showed no survival gain at all.  Distant spread
// in MPM is late and is rarely what kills; a ceiling and a slower rate were
// required.
// ---------------------------------------------------------------------------
double met_supp = pos(1.0 - kill / 0.050);
dxdt_MET = KMET * nN + KMETG * MET * met_supp * (1.0 - MET / METMAX);

// ---- 10. marrow -------------------------------------------------------------------
double E_pem_m = EMAXPM * PEM_MP / (PEM_MP + EC50PM * (1.0 + FOL / KFPM));
double Edrug = E_pem_m + SLOPE_CIS * CP_CIS + 0.35 * VINE;
if (Edrug > 0.95) Edrug = 0.95;
double fb = pow(CIRC0 / anc, GAM);
dxdt_PROL = KTR * PROL * (1.0 - Edrug) * fb - KTR * PROL;
dxdt_TR1  = KTR * (PROL - TR1);
dxdt_TR2  = KTR * (TR1 - TR2);
dxdt_TR3  = KTR * (TR2 - TR3);
dxdt_CIRC = KTR * TR3 - KTR * CIRC;

// ---- 11. renal ---------------------------------------------------------------------
dxdt_CRCLSS = - KNEPHIRR * CP_CIS;
dxdt_CRCL   = KRECR * (CRCLSS - CRCL) - KNEPH * CP_CIS;

// ---- 12. biomarker --------------------------------------------------------------------
dxdt_SMRP = KSMRP * nN * MSLNF - KSMRPCL * (crcl / 95.0) * SMRP;

// ---- 13. inflammation, wound, cachexia --------------------------------------------------
dxdt_IL6   = KIL6T * (nN / 100.0) + KIL6W * WOUND - KIL6D * IL6;
dxdt_WOUND = - KWD * WOUND;
dxdt_LBM   = - KCACH * IL6 / (IL6 + KMIL6) + KLREC * (LBM0 - LBM);

// ---- 14. immune-related adverse events -----------------------------------------------------
dxdt_IRAE = 0.055 * (RO_CT4 + 0.30 * RO_PD1) * (1.0 - IRAE) - 0.030 * IRAE;

// ---- 15. hazard ------------------------------------------------------------------------------
// The hazard is NOT proportional to volume.  What kills in MPM is encasement
// -- the fraction of the hemithorax that no longer moves -- and the
// performance status that follows from it.  Volume enters with a square-root
// exponent; AREA enters in its own right.  That is the same structural claim
// as the rest of the model, and it is what makes debulking (which collapses h
// and leaves A) a poor bargain.
double fvc = FVC0 * (1.0 - A_FVC * (afrac > 1.0 ? 1.0 : afrac)
                         - H_FVC * (h_tot > 3.0 ? 3.0 : h_tot) / 3.0
                         - V_FVC * (VEFF > 2500.0 ? 2500.0 : VEFF) / 2500.0);
double pain = 3.0 * Z / (Z + 0.80);
double psidx = 1.10 * pos(1.0 - fvc / FVC0)
             + 1.30 * pos(1.0 - LBM / LBM0) * 3.0
             + 0.12 * pain;
double ps = 4.0 * (1.0 - exp(-1.35 * psidx));  if (ps > 4.0) ps = 4.0;

double haz = H0 * pow((Vtum < 1.0 ? 1.0 : Vtum) + WMET * MET, B1) / pow(V0H, B1)
             * exp(B2 * ps + B5 * (afrac > 1.0 ? 1.0 : afrac)) * SARCH;
if (CIRC < 0.5) haz += H_FEBN;
haz += H_IRAE * pos(IRAE - 0.35) / 0.65;
dxdt_CUMH = haz;

$TABLE
double hmeas = (pos(N) + PHIM * pos(M)) / (A < 1.0 ? 1.0 : A);
double hvia  = pos(N) / (A < 1.0 ? 1.0 : A);
double VTUM  = pos(N) + PHIM * pos(M);
double AFR   = (A / S_PL > 1.0) ? 1.0 : A / S_PL;
double MRECIST = 60.0 * hmeas;                       // 6 sites, mm
double VEGFF = pos(VEGF) / (1.0 + C_BEV / KDBEV);
double IFPo  = IFP_BASE + IFP_MAX * VEGFF / (VEGFF + EC50I);
double LPMM  = LD0_CIS * pow(RHOV < 1e-3 ? 1e-3 : RHOV, PEXP_RHO)
               * pow(IFP_REF / IFPo, PEXP_IFP) / sqrt(COLLF) * 10.0;
double FEXP  = fpen(LPMM / 10.0, hmeas);
double SANCT = 1.0 - FEXP;
double FIBFR = PHIM * pos(M) / (VTUM < 1e-9 ? 1e-9 : VTUM);
double FVCP  = FVC0 * (1.0 - A_FVC * AFR - H_FVC * (hmeas > 3.0 ? 3.0 : hmeas) / 3.0
                           - V_FVC * (VEFF > 2500.0 ? 2500.0 : VEFF) / 2500.0);
double PAINi = 3.0 * Z / (Z + 0.80);
double PSIDX = 1.10 * pos(1.0 - FVCP / FVC0) + 1.30 * pos(1.0 - LBM / LBM0) * 3.0
               + 0.12 * PAINi;
double ECOG  = 4.0 * (1.0 - exp(-1.35 * PSIDX));  if (ECOG > 4.0) ECOG = 4.0;
double SURV  = exp(-CUMH);
double CPEM  = PEM_C / V1_PEM;
double CCIS  = CIS_F / V1_CIS;
double CBEV  = BEV_C / V1_BEV;
double CNIV  = NIV_C / V1_NIV;
double CIPI  = IPI_C / V1_IPI;
double ROPD1 = CNIV / (CNIV + EC50_PD1);
double ROCT4 = CIPI / (CIPI + EC50_CT4);

$CAPTURE
hmeas hvia VTUM AFR MRECIST LPMM FEXP SANCT FIBFR FVCP ECOG SURV PAINi
CPEM CCIS CBEV CNIV CIPI ROPD1 ROCT4 VEGFF IFPo
'

mod <- mcode_cache("mpm", mpm_code, end = 1250, delta = 1)

# =============================================================================
#  HISTOLOGY
# =============================================================================
HIST <- list(
  epithelioid = list(CHEMOS = 1.00, IMMINF = 0.38, VISTA_S = 1.15, COLLF = 1.00,
                     MSLNF = 1.00, SARCH = 1.00, KPROLF = 1.00),
  biphasic    = list(CHEMOS = 0.52, IMMINF = 0.70, VISTA_S = 0.60, COLLF = 1.45,
                     MSLNF = 0.50, SARCH = 1.35, KPROLF = 1.25),
  sarcomatoid = list(CHEMOS = 0.18, IMMINF = 1.00, VISTA_S = 0.20, COLLF = 1.90,
                     MSLNF = 0.12, SARCH = 1.90, KPROLF = 1.55)
)
# histology mix of the modern first-line trials (CheckMate 743: 76 % epithelioid)
COHORT <- data.frame(hist = c("epithelioid", "biphasic", "sarcomatoid"),
                     w    = c(0.75, 0.13, 0.12), stringsAsFactors = FALSE)

CMTNUM <- function(name) which(mrgsolve::cmt(mod) == name)

# =============================================================================
#  REGIMEN BUILDERS
# =============================================================================
CYC     <- 21                       # days
INF_PEM <- 10 / 60 / 24             # 10-minute infusion, in days
INF_CIS <- 2 / 24                   # 2-hour infusion
BSA     <- 1.80
WTKG    <- 70

# NOTE (DEFECT #4): pemetrexed is a ten-minute infusion, which is shorter than
# any sensible fixed integration step.  mrgsolve inserts exact dose records and
# integrates each infusion window on its own, so `rate = amt/dur` is safe here;
# the hand-written RK4 driver in mpm_reference_model.py had to be given explicit
# breakpoints at every rate discontinuity before it delivered the right dose.
ev_pem <- function(t, amt) ev(time = t, amt = amt, cmt = "PEM_C",
                              rate = amt / INF_PEM)
ev_cis <- function(t, amt) ev(time = t, amt = amt, cmt = "CIS_F",
                              rate = amt / INF_CIS)
ev_bev <- function(t, amt) ev(time = t, amt = amt, cmt = "BEV_C", rate = amt / 0.03)
ev_niv <- function(t, amt) ev(time = t, amt = amt, cmt = "NIV_C", rate = amt / 0.03)
ev_ipi <- function(t, amt) ev(time = t, amt = amt, cmt = "IPI_C", rate = amt / 0.03)
ev_adi <- function(t, amt) ev(time = t, amt = amt, cmt = "ADI_A")
ev_ip  <- function(t, amt) ev(time = t, amt = amt, cmt = "IPD",  rate = amt / 0.04)
ev_vin <- function(t, amt) ev(time = t, amt = amt, cmt = "VINE", rate = amt / 0.03)

#' Fixed (non-adaptive) platinum-pemetrexed
chemo_events <- function(ncyc = 6, t0 = 0, pem = TRUE, plat = "cis",
                         level = 1) {
  dose_pem <- 500 * BSA
  dose_pt  <- (if (plat == "cis") 75 else 0.78 * 75) * BSA
  e <- NULL
  for (c in seq_len(ncyc) - 1L) {
    t <- t0 + c * CYC
    if (pem) e <- c(e, list(ev_pem(t, dose_pem * level)))
    e <- c(e, list(ev_cis(t + 0.5 / 24, dose_pt * level)))
  }
  Reduce(`+`, e)
}

io_events <- function(until = 730, t0 = 0, nivo = TRUE, ipi = TRUE,
                      nivo_q = 14, nivo_mg = 3 * WTKG, ipi_q = 42,
                      ipi_mg = 1 * WTKG) {
  e <- list()
  if (nivo) for (c in seq_len(floor(until / nivo_q)) - 1L)
    e <- c(e, list(ev_niv(t0 + c * nivo_q, nivo_mg)))
  if (ipi)  for (c in seq_len(floor(until / ipi_q)) - 1L)
    e <- c(e, list(ev_ipi(t0 + c * ipi_q, ipi_mg)))
  if (!length(e)) return(NULL)
  Reduce(`+`, e)
}

bev_events <- function(until = 730, t0 = 0, mgkg = 15) {
  e <- lapply(seq_len(floor(until / CYC)) - 1L,
              function(c) ev_bev(t0 + c * CYC + 0.02, mgkg * WTKG))
  Reduce(`+`, e)
}

adi_events <- function(until = 365, t0 = 0) {
  e <- lapply(seq_len(floor(until / 7)) - 1L,
              function(c) ev_adi(t0 + c * 7, 36 * BSA))
  Reduce(`+`, e)
}

ip_events <- function(ncyc = 6, t0 = 1, mgm2 = 100) {
  e <- lapply(seq_len(ncyc) - 1L, function(c) ev_ip(t0 + c * CYC, mgm2 * BSA))
  Reduce(`+`, e)
}

# =============================================================================
#  ADAPTIVE CHEMOTHERAPY
#  A cycle is given only if ANC >= 1.5 and CrCl >= 45; otherwise it is delayed
#  a week (up to three times), the dose steps down to 75 % then 50 %, and five
#  consecutive failed weeks stop treatment.  This matters because it is the
#  channel through which folate supplementation changes SURVIVAL and not just
#  toxicity: in the equations alone, folate only touches the marrow.
# =============================================================================
sim_adaptive_chemo <- function(...) {
  # Placeholder.  The protocol dose-modification loop (ANC >= 1.5 and
  # CrCl >= 45 before each cycle, one-week delays, 75 % then 50 % dose levels,
  # permanent stop after five failed weeks) is implemented and exercised in
  # mpm_reference_model.py::add_chemo(); porting it here means driving mrgsim
  # cycle by cycle and carrying the end state forward with init().  It is left
  # unimplemented rather than sketched, so that nothing in this file looks
  # executed that was not.
  stop("not implemented in the R port - see add_chemo() in mpm_reference_model.py")
}

#' Run one scenario.  `events` is an mrgsolve event object, `pars` a named
#' list of parameter overrides, `surgery` an optional list(time=, h_res=,
#' mort=, fvc_loss=, hazx=, lbm_loss=), `windows` a data.frame of
#' time-varying parameter windows (name, t0, t1, value).
sim_scenario <- function(hist = "epithelioid", events = NULL, pars = list(),
                         surgery = NULL, windows = NULL, tmax = 1250,
                         supplemented = NA) {
  p <- modifyList(HIST[[hist]], pars)
  if (!is.na(supplemented)) {
    p$FOLUG <- if (supplemented) 400 else 0
    p$B12   <- if (supplemented) 1 else 0
  }
  m <- param(mod, p)

  # time-varying parameters (talc, radiotherapy, PRMT5i, drainage) are carried
  # as columns of the input data set
  if (!is.null(windows)) {
    tt <- sort(unique(c(0, tmax, windows$t0, windows$t1)))
    dset <- data.frame(ID = 1, time = tt, amt = 0, evid = 2, cmt = 0)
    for (nm in unique(windows$name)) dset[[nm]] <- 0
    for (i in seq_len(nrow(windows))) {
      sel <- dset$time >= windows$t0[i] & dset$time < windows$t1[i]
      dset[[windows$name[i]]][sel] <- windows$value[i]
    }
    if (!is.null(events)) dset <- rbind_ev(dset, events)
    sim <- mrgsim_d(m, dset, end = tmax, delta = 1, carry_out = "evid")
    return(as.data.frame(sim))
  }

  if (is.null(surgery)) {
    sim <- mrgsim(m, events = events, end = tmax, delta = 1)
    return(as.data.frame(sim))
  }

  # --- surgery: split the horizon, edit the state, resume ------------------
  ts <- surgery$time
  e1 <- if (!is.null(events)) filter_ev(events, time < ts) else NULL
  s1 <- as.data.frame(mrgsim(m, events = e1, end = ts, delta = 1))
  last <- s1[nrow(s1), ]
  Ares <- last$A
  keep <- surgery$h_res * Ares
  tot  <- last$N + p$PHIM %||% 0.70 * last$M
  tot  <- last$N + 0.70 * last$M
  frac <- min(1, keep / max(tot, 1e-9))
  init2 <- as.list(last[mrgsolve::cmt(mod)])
  init2$N     <- last$N * frac
  init2$M     <- last$M * frac
  init2$VEFF  <- 50
  init2$PSY   <- 0.85
  init2$WOUND <- 1
  init2$LBM   <- last$LBM * (1 - surgery$lbm_loss)
  init2$CUMH  <- last$CUMH + surgery$mort
  # the two permanent costs of the operation
  p2 <- p
  p2$FVC0 <- (p$FVC0 %||% 92) * (1 - surgery$fvc_loss)
  p2$H0   <- (p$H0   %||% 0.000170) * surgery$hazx
  m2 <- param(mod, p2)
  m2 <- init(m2, init2)
  e2 <- if (!is.null(events)) filter_ev(events, time >= ts) else NULL
  if (!is.null(e2)) e2 <- mutate_ev(e2, time = time - ts)
  s2 <- as.data.frame(mrgsim(m2, events = e2, end = tmax - ts, delta = 1))
  s2$time <- s2$time + ts
  rbind(s1[-nrow(s1), ], s2)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# =============================================================================
#  SURVIVAL READOUTS
# =============================================================================
median_from_S <- function(time, S) {
  i <- which(S <= 0.5)[1]
  if (is.na(i) || i == 1) return(NA_real_)
  s0 <- S[i - 1]; s1 <- S[i]
  time[i - 1] + (s0 - 0.5) / (s0 - s1) * (time[i] - time[i - 1])
}

#' Weighted mixture over histology -> a trial-cohort survival curve.
cohort_curve <- function(runner, tmax = 1250) {
  cur <- lapply(COHORT$hist, function(h) runner(h)$SURV)
  tt  <- runner(COHORT$hist[1])$time
  S   <- Reduce(`+`, Map(function(w, s) w * s, COHORT$w, cur))
  list(time = tt, S = S, median_mo = median_from_S(tt, S) / 30.44)
}

#' Hazard ratio between two arms, read at 24 months.
model_HR <- function(curveA, curveB, day = 730) {
  (-log(max(curveA$S[day + 1], 1e-12))) / (-log(max(curveB$S[day + 1], 1e-12)))
}

# =============================================================================
#  THE 24 SCENARIOS
# =============================================================================
SCEN <- list(

  "01 BSC / untreated" = function(h)
    sim_scenario(h),

  "02 cisplatin alone x6" = function(h)
    sim_scenario(h, chemo_events(6, pem = FALSE), supplemented = FALSE),

  "03 cis+pem, NO supplementation" = function(h)
    sim_scenario(h, chemo_events(6), supplemented = FALSE),

  "04 cis+pem, folate + B12" = function(h)
    sim_scenario(h, chemo_events(6), supplemented = TRUE),

  "05 carboplatin AUC5 + pem" = function(h)
    sim_scenario(h, chemo_events(6, plat = "carbo"), supplemented = TRUE,
                 pars = list(KNEPH = 11.75, KNEPHIRR = 3.0)),

  "06 cis+pem+bevacizumab (MAPS)" = function(h)
    sim_scenario(h, chemo_events(6) + bev_events(730), supplemented = TRUE),

  "07 nivolumab+ipilimumab (CM743)" = function(h)
    sim_scenario(h, io_events(730)),

  "08 pembrolizumab + chemo (IND.227)" = function(h)
    sim_scenario(h, chemo_events(6) +
                    io_events(730, ipi = FALSE, nivo_q = 21, nivo_mg = 200),
                 supplemented = TRUE),

  "09 chemo -> nivolumab 2L (CONFIRM)" = function(h)
    sim_scenario(h, chemo_events(6) +
                    io_events(730, t0 = 200, ipi = FALSE), supplemented = TRUE),

  "10 ADI-PEG20 + chemo (ATOMIC)" = function(h)
    sim_scenario(h, chemo_events(6) + adi_events(365), supplemented = TRUE),

  "11 extended P/D + chemo (MARS2)" = function(h)
    sim_scenario(h, chemo_events(6), supplemented = TRUE,
                 surgery = list(time = 70, h_res = 0.020, mort = 0.040,
                                fvc_loss = 0.14, hazx = 1.30, lbm_loss = 0.10)),

  "12 EPP + chemo + hemithoracic RT" = function(h)
    sim_scenario(h, chemo_events(3), supplemented = TRUE,
                 surgery = list(time = 70, h_res = 0.012, mort = 0.075,
                                fvc_loss = 0.34, hazx = 1.55, lbm_loss = 0.14),
                 windows = data.frame(name = "RTRATE", t0 = 110, t1 = 145,
                                      value = 0.055)),

  "13 chemo + intrapleural cisplatin" = function(h)
    sim_scenario(h, chemo_events(6) + ip_events(6), supplemented = TRUE),

  "14 talc pleurodesis then chemo" = function(h)
    sim_scenario(h, chemo_events(6, t0 = 14), supplemented = TRUE,
                 windows = data.frame(name = "TALCON", t0 = 5, t1 = 7, value = 1)),

  "15 talc, then chemo + intrapleural" = function(h)
    sim_scenario(h, chemo_events(6, t0 = 14) + ip_events(6, t0 = 15),
                 supplemented = TRUE,
                 windows = data.frame(name = "TALCON", t0 = 5, t1 = 7, value = 1)),

  "16 chemo + PRMT5i (MTAP-deleted)" = function(h)
    sim_scenario(h, chemo_events(6), supplemented = TRUE,
                 windows = data.frame(name = "PRMT5", t0 = 0, t1 = 730,
                                      value = 0.0075)),

  "17 chemo + pemetrexed maintenance" = function(h)
    sim_scenario(h, chemo_events(6) + chemo_events(18, t0 = 6 * CYC,
                                                   pem = TRUE, plat = "none"),
                 supplemented = TRUE),

  "18 chemo with CrCl 52 mL/min" = function(h)
    sim_scenario(h, chemo_events(6), supplemented = TRUE,
                 pars = list(CRCL0 = 52)),

  "19 chemo -> vinorelbine 2L (VIM)" = function(h) {
    e <- chemo_events(6)
    for (c in seq_len(18) - 1L) e <- e + ev_vin(200 + c * 7, 60 * BSA * 0.012)
    sim_scenario(h, e, supplemented = TRUE)
  },

  "20 chemo, EARLY rind 0.35 cm" = function(h)
    sim_scenario(h, chemo_events(6), supplemented = TRUE,
                 pars = list(H0_RIND = 0.35, A0_RIND = 250)),

  "21 chemo, BULKY rind 2.0 cm" = function(h)
    sim_scenario(h, chemo_events(6), supplemented = TRUE,
                 pars = list(H0_RIND = 2.00, A0_RIND = 800)),

  "22 intrapleural cisplatin alone" = function(h)
    sim_scenario(h, ip_events(6)),

  "23 IO first, chemo at progression" = function(h)
    sim_scenario(h, io_events(250) + chemo_events(6, t0 = 260),
                 supplemented = TRUE),

  "24 chemo, IO added day 140" = function(h)
    sim_scenario(h, chemo_events(6) + io_events(600, t0 = 140),
                 supplemented = TRUE)
)

# =============================================================================
#  DEMO
# =============================================================================
run_all <- function() {
  res <- lapply(names(SCEN), function(nm) {
    f <- SCEN[[nm]]
    cc <- cohort_curve(function(h) f(h))
    epi <- f("epithelioid")
    data.frame(
      scenario   = nm,
      median_OS_mo = round(cc$median_mo, 1),
      best_mRECIST = round(min((epi$MRECIST[1:300] - epi$MRECIST[1]) /
                                 epi$MRECIST[1] * 100), 1),
      viable_nadir = round(min((epi$N[1:400] - epi$N[1]) / epi$N[1] * 100), 1),
      ANC_nadir    = round(min(epi$CIRC[1:200]), 2),
      FVC_6mo      = round(epi$FVCP[183], 1),
      peak_effusion = round(max(epi$VEFF)),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}

if (sys.nframe() == 0L) {
  print(run_all(), row.names = FALSE)
}
