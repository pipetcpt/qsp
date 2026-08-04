## =====================================================================
##  PATENT DUCTUS ARTERIOSUS OF PREMATURITY  --  QSP / PK-PD MODEL
##  mrgsolve implementation
## =====================================================================
##
##  The ductus arteriosus is not a lesion.  It is a normal fetal vessel
##  held open by prostaglandin E2 acting through EP4 -> Gs -> cAMP in its
##  own smooth muscle.  In a term infant, losing the placental PGE2
##  source and raising PaO2 closes it in 24-72 h and then remodels it
##  permanently.  In a 25-week infant it often does neither.
##
##  ------------------------------------------------------------------
##  FOUR STRUCTURAL CLAIMS THIS MODEL MAKES
##  ------------------------------------------------------------------
##  1. TWO ENZYME SITES, NOT ONE DRUG CLASS.  Indomethacin and ibuprofen
##     compete with arachidonate in the cyclooxygenase CHANNEL.
##     Acetaminophen reduces the ferryl-protoporphyrin radical at the
##     physically separate PEROXIDASE site, so it competes with
##     PEROXIDE, not with substrate.  Written that way -- I_chan additive
##     in occupancy, I_perox multiplicative with an IC50 that scales with
##     peroxide tone -- the model PREDICTS that acetaminophen loses
##     potency under chorioamnionitis while ibuprofen does not.
##
##  2. CONSTRICTION IS NOT CLOSURE.  Permanent closure needs ductal WALL
##     hypoxia: the constricted lumen stops supplying O2 by diffusion,
##     the media goes hypoxic, and HIF-1a/VEGF/TGF-b1 build a neointimal
##     cushion.  A preterm ductal wall is THIN enough to stay oxygenated
##     by diffusion even while constricted, so the remodelling signal
##     never fires and the duct REOPENS on drug washout.  WALLO2 is the
##     bistability switch and it is why ~1 in 4 preterm ducts reopens
##     while a term duct essentially never does.
##
##  3. SPONTANEOUS CLOSURE IS THE COMPETING PROCESS.  Local ductal
##     prostanoid synthesis involutes postnatally (FSYN), faster at
##     higher gestational age, and this closes almost every preterm duct
##     eventually -- median ~71 days below 26 weeks (Semberova 2017).
##     Any drug effect on a fixed-time closure endpoint is therefore
##     bounded above by the untreated non-closure fraction at that time.
##     That ceiling, not drug potency, is what the negative trials
##     measured.
##
##  4. TONE IS A SIGMOID OF NET DRIVE, NOT A PRODUCT OF GAINS.  An
##     earlier product form (tone = Tmax x O2gain x (1-relax_PGE2) x
##     (1-relax_NO)) made achievable tone the arithmetic product of every
##     immaturity factor, so a 26-week duct could never close at any
##     exposure.  That is wrong -- it closes in roughly half of treated
##     infants.  The gestational limit belongs in TMAXGA, the maximal
##     occlusion a thin immature media can generate: at 24 wk the ceiling
##     leaves a 0.44 mm residual lumen, above the 0.30 mm closure
##     threshold, so no exposure of any drug closes that duct; at 26 wk
##     it leaves 0.24 mm and closure becomes achievable.  The steep
##     clinical GA gradient in treatment success falls OUT of the model
##     rather than being fitted INTO it.
##
##  ------------------------------------------------------------------
##  FITTED VS PREDICTED  (read before trusting a number)
##  ------------------------------------------------------------------
##  FITTED:
##    NET50, NETW, TAUSYN0, KTAUSYN, KINVGA
##        by Nelder-Mead on 8 gestational spontaneous-closure times
##        (Semberova 2017; Nemerofsky 2008; Rolland 2015)
##    KI_IBU, KI_IND, IC50_APAP
##        to pooled day-7 closure rates by drug (Cochrane reviews; the
##        JAMA 2018 network meta-analysis)
##    CLMAT_IBU
##        to the observed loss of efficacy of a fixed regimen with
##        advancing postnatal age
##    B_GA, B_BUR, B_VENT   (3 coefficients, 2 proportions + the GA
##        gradient) to the Baby-OSCAR primary endpoint in BOTH arms
##  PREDICTED -- nothing touched afterwards:
##    reopening rate and its GA dependence; loss of acetaminophen effect
##    under high peroxide tone; high-dose vs standard late rescue;
##    the BeNeDuctus composite in both arms; the TIPP IVH reduction
##    without outcome benefit; renal/cerebral/mesenteric separation of
##    the two NSAIDs; the entire targeted-treatment scenario.
##
##  ------------------------------------------------------------------
##  IMPLEMENTATION NOTES
##  ------------------------------------------------------------------
##  * 51 ODE compartments.  Time in HOURS of postnatal age.  Drug amounts
##    are per kg (mg/kg), so concentrations are mg/L directly.
##  * All shared algebra lives in ONE C++ macro, PDA_ALG(TT), defined in
##    $GLOBAL and invoked from both $ODE and $TABLE.  Duplicating that
##    algebra in two blocks is the classic way these models silently
##    drift, where the reported output stops matching the equations that
##    were integrated.  The macro makes drift impossible.
##  * A Python re-implementation (pda_reference_model.py) integrates the
##    same equations with LSODA, and a script diffs the two parameter sets
##    entry by entry (176 shared parameters, all agreeing exactly; the only
##    Python-only entries are the two infusion rates, which mrgsolve
##    expresses through the event object's `rate` column instead).
##    Development found SIX real defects, each documented in the source at
##    the point where it was fixed: (1) a ductal-resistance gain 56x too
##    large, which capped Qp:Qs at 1.44 and silently zeroed every
##    shunt-driven output; (2) a lung-water -> PVR feedback strong enough
##    for the duct to throttle its own flow; (3) the product-form tone
##    equation above; (4) a threshold fitted to untreated data alone, which
##    put closure beyond the reach of any drug; (5) a wall-O2 form with the
##    right sign and the wrong magnitude, giving no reopening at any
##    gestation; and (6) delta-creatinine reported as a difference of window
##    maxima, which is identically zero for a monotonically falling
##    neonatal creatinine however nephrotoxic the drug is.
##
##  DISCLAIMER: educational / research model.  Not validated for clinical
##  use, prescribing, or regulatory submission.
## =====================================================================

library(mrgsolve)
library(dplyr)

pda_code <- '
$PARAM @annotated
// ================= patient =================
GA       :  26.0   : gestational age at birth (weeks)
BW       :   0.80  : birth weight (kg)
PAO2     :  55.0   : arterial PO2 (mmHg)
SEPSIS   :   0.0   : chorioamnionitis / early-onset sepsis flag (0/1)
HCORT    :   0.0   : concomitant early hydrocortisone (0/1)
ANTESTER :   1.0   : complete antenatal corticosteroid course (0/1)

// ================= ibuprofen PK (per kg) =================
// Reported preterm values are wide AND mutually inconsistent: Aranda 1997
// gives CL 2.06 mL/h/kg with Vd 0.062 L/kg, which implies t1/2 21 h (not
// the 30.5 h reported alongside) and a Cmax of 160 mg/L after 10 mg/kg,
// which no assay report supports.  We adopt a SELF-CONSISTENT set inside
// the reported envelope (t1/2 20-43 h, CL 1.4-4.5 mL/h/kg, Vss 0.06-0.25
// L/kg).  Only unbound drug reaches COX, so the absolute total-
// concentration scale is absorbed into the fitted Ki regardless.
V1_IBU   :   0.140   : central volume, ibuprofen (L/kg)
V2_IBU   :   0.080   : peripheral volume, ibuprofen (L/kg)
Q_IBU    :   0.012   : intercompartmental clearance, ibuprofen (L/h/kg)
CL_IBU0  :   0.00420 : clearance at postnatal age 0, ibuprofen (L/h/kg)
CLMAT_IBU:   0.20    : fractional increase in ibuprofen CL per day [FITTED]
FU_IBU   :   0.0120  : unbound fraction, ibuprofen (preterm)
MW_IBU   : 206.3     : molecular weight, ibuprofen (g/mol)
KE0_IBU  :   0.35    : ductal effect-site equilibration, ibuprofen (1/h)

// ================= indomethacin PK (per kg) =================
V1_IND   :   0.250   : central volume, indomethacin (L/kg)
V2_IND   :   0.120   : peripheral volume, indomethacin (L/kg)
Q_IND    :   0.025   : intercompartmental clearance, indomethacin (L/h/kg)
CL_IND0  :   0.01160 : clearance at postnatal age 0, indomethacin (L/h/kg)
CLMAT_IND:   0.12    : fractional increase in indomethacin CL per day
FU_IND   :   0.010   : unbound fraction, indomethacin
MW_IND   : 357.8     : molecular weight, indomethacin (g/mol)
KE0_IND  :   0.12    : effect-site equilibration, indomethacin (1/h) - slow,
                       representing slowly-reversible tight COX binding

// ================= acetaminophen PK (per kg) =================
V1_APAP  :   0.60    : central volume, acetaminophen (L/kg)
V2_APAP  :   0.40    : peripheral volume, acetaminophen (L/kg)
Q_APAP   :   0.30    : intercompartmental clearance, acetaminophen (L/h/kg)
CL_APAP0 :   0.120   : clearance at postnatal age 0, acetaminophen (L/h/kg)
CLMAT_APAP:  0.06    : fractional increase in acetaminophen CL per day
KA_APAP  :   0.60    : enteral absorption rate constant (1/h)
F_APAP   :   0.85    : enteral bioavailability
FU_APAP  :   0.85    : unbound fraction, acetaminophen
MW_APAP  : 151.2     : molecular weight, acetaminophen (g/mol)

// ================= acetaminophen hepatic safety =================
FNAPQI   :   0.030   : fraction of acetaminophen clearance via CYP oxidation
KGSH     :   0.020   : hepatic glutathione resynthesis (1/h)
GSHCAP   :   1.0     : glutathione capacity (fraction)
KALT     :   0.60    : ALT production gain per unit unquenched NAPQI
KALTOUT  :   0.0173  : ALT elimination (1/h)
ALT0     :  20.0     : baseline ALT (U/L)

// ================= COX enzymology at the ductus =================
KI_IBU   :   0.220   : apparent unbound Ki, ibuprofen, arachidonate channel (uM) [FITTED]
KI_IND   :   0.00083 : apparent unbound Ki, indomethacin (uM) [FITTED]
IC50_APAP:  6.0     : unbound IC50, acetaminophen at peroxidase site, basal peroxide (uM) [FITTED]
KPEROX   :   1.00    : peroxide tone that doubles the acetaminophen IC50

// ================= prostanoid dynamics =================
KSYN_PGE2D:  0.90    : ductal PGE2 synthesis rate (1/h, normalised)
KDEG_PGE2D:  0.90    : ductal PGE2 degradation (1/h)
FSYN_INF :   0.18    : asymptotic fraction of fetal ductal synthesis capacity
TAUSYN0  : 200.0     : postnatal involution time constant at 26 wk (h) [FITTED]
KTAUSYN  :   0.140   : per-week acceleration of synthesis involution [FITTED]
KIN_PGE2D:   0.055   : ductal uptake gain from plasma PGE2
PGE2P0   : 700.0     : fetal circulating PGE2 (pg/mL)
RPROD_PGE2P: 22.0    : postnatal extra-placental PGE2 production (pg/mL/h)
KPUL_PGE2:   0.055   : pulmonary PGE2 clearance gain (1/h, x 15-PGDH x flow)
KOTH_PGE2:   0.010   : non-pulmonary PGE2 clearance (1/h)
PGDH0_GA :  22.0     : gestational age at which pulmonary 15-PGDH ~ 0 (wk)
PGDH_SL  :   0.070   : per-week gestational gain in 15-PGDH
KPGDH    :   0.0060  : postnatal 15-PGDH induction (1/h)
KCAMP    :   2.0     : ductal cAMP turnover (1/h)
KNO      :   0.10    : ductal NO tone turnover (1/h)
NOGA     :   0.055   : per-week decline of ductal NO dependence with GA
KNODEC   :   0.0300  : postnatal loss of ductal NO tone (1/h)
KEP4     :   0.00600 : postnatal EP4 down-regulation at 26 wk (1/h)
EP4GA    :   0.075   : per-week lower EP4 density with advancing GA
EP4INF   :   0.30    : asymptotic EP4 fraction
KINVGA   :   0.140   : per-week acceleration of the whole involution programme [FITTED]
KPEROXIN :   0.030   : ductal peroxide tone turnover (1/h)
PEROXSEP :   3.20    : peroxide tone multiplier with sepsis / chorioamnionitis

// ================= ductal contraction =================
TMAX_GA50:  21.5     : gestational age at half-maximal achievable occlusion (wk)
TMAX_SL  :   0.55    : steepness of the occlusion ceiling in GA (per wk)
// These weights matter STRUCTURALLY, not just numerically.  An earlier version
// had GMAT 0.42 with slow NO decay, so most of the postnatal rise in net drive
// came from contractile maturation and loss of NO -- two terms NO DRUG TOUCHES.
// Fitting NET50 to untreated closure times then placed the threshold above
// anything PGE2 suppression could reach, and every drug became inert: peak
// ductal tone 0.02 at 92% COX inhibition.  The fix is not a bigger dose, it is
// recognising that the dominant driver of postnatal constriction IS the decline
// of the PGE2/EP4 drive -- exactly the axis these drugs act on.
GO2      :   1.00    : weight of O2-driven constriction in net drive
GMAT     :   0.10    : weight of postnatal contractile maturation
GPGE2    :   1.00    : weight of PGE2/EP4 relaxation
GNO      :   0.35    : weight of NO relaxation
NET50    :   0.350   : net drive at half-maximal tone [FITTED]
NETW     :   0.0600   : steepness of the net-drive to tone map [FITTED]
TAU_MAT  : 300.0     : postnatal contractile maturation time constant at 26 wk (h)
KTONE    :   0.28    : ductal tone equilibration (1/h)
P50_O2   :  30.0     : PaO2 for half-maximal O2-driven tone at term (mmHg)
P50_GA   :   0.130   : per-week rise in P50 below term (less O2-sensitive)
EMAX_PGE2:   0.97    : maximal relaxation attributable to PGE2/EP4
EC50_PGE2:   0.300   : normalised PGE2 x EP4 drive for half-maximal relaxation
EMAX_NO  :   0.55    : maximal relaxation attributable to NO
EC50_NO  :   0.45    : normalised NO tone for half-maximal relaxation
ALPHA_D  :   0.965   : fraction of diameter abolished at full tone
DMAX_GA  :   0.135   : unconstricted ductal diameter, per week of GA (mm)
DMAX_INT :  -1.35    : unconstricted ductal diameter intercept (mm)
DCLOSE   :   0.30    : diameter below which the shunt is trivial (mm)

// ================= wall O2 and remodelling =================
WALLTH_GA:   0.075   : ductal wall thickness index per week of GA
WALLTH_INT: -1.20    : ductal wall thickness index intercept
KWALL    :   0.0025  : postnatal wall thickening (1/h)
FLUM_O2  :   0.80    : weight of luminal O2 supply to the ductal wall
FVASA_O2 :   0.20    : weight of vasa vasorum supply
KDIFF    :   2.60    : inverse O2 diffusion length (per unit wall index)
O2CRIT   :  16.0     : wall PO2 below which hypoxic signalling fires (mmHg)
KVEGF    :   0.045   : VEGF induction per mmHg of wall hypoxia (1/h)
KVEGFOUT :   0.030   : VEGF elimination (1/h)
KTGFB    :   0.020   : TGF-beta1 induction by VEGF (1/h)
KTGFBOUT :   0.015   : TGF-beta1 elimination (1/h)
KREMOD   :   0.00060  : neointimal remodelling gain (1/h)
// REMODELLING COMPETENCE is itself gestationally immature, and this term is what
// reproduces the sharpest feature in the data.  the Semberova medians jump 5.5-fold
// -- 71 d below 26+0 wk versus 13 d at 26+0-27+6 -- across about one week of
// gestation.  Above the occlusion threshold (~25.5 wk) tone alone closes the duct
// and that is fast; below it, closure REQUIRES a neointimal cushion.  Without
// this factor the model built the cushion almost as fast in a 24-week duct as in
// a 30-week one and closed it in ~33 d instead of ~85, and neither KREMOD nor
// KWALL fixed it: scanning KREMOD over 6.5-fold moved closure only 30->38 d,
// because remodelling is not rate-limiting once hypoxia is on.  The missing
// physiology was the CAPACITY of the duct to remodel -- the very preterm ductus has
// sparse vasa vasorum and immature endothelial/SMC migration and hyaluronan
// machinery, so it cannot build the cushion.
KREMGA50 :  25.4     : gestational age at half-maximal remodelling competence (wk)
KREMGA_SL:   1.20    : steepness of remodelling competence in GA (per wk)
KUNREMOD :   0.00004 : neointimal regression (1/h) - near-irreversible

// ================= haemodynamics =================
PVR0     :   0.115   : pulmonary vascular resistance at birth (mmHg/(mL/min/kg))
PVRINF   :   0.015   : asymptotic pulmonary vascular resistance
KPVR     :   0.022   : postnatal PVR fall (1/h)
PVR_EVLW :   0.00120 : PVR rise per mL/kg of lung water
SVR0     :   0.165   : systemic vascular resistance at birth
SVRINF   :   0.205   : asymptotic systemic vascular resistance
KSVR     :   0.004   : postnatal SVR rise (1/h)
QS_TARGET: 180.0     : desired systemic blood flow (mL/min/kg)
QMAX_LV  : 400.0     : maximal LV output at 26 wk (mL/min/kg)
QMAX_GA  :  20.0     : added LV reserve per week of GA (mL/min/kg)
PLA0     :   4.0     : baseline left atrial pressure (mmHg)
KPLA     :   0.35    : left atrial pressure equilibration (1/h)
PLA_Q    :   0.022   : LA pressure rise per mL/min/kg of shunt (mmHg)
KRDUCT   :   1.10    : ductal resistance gain (mmHg/(mL/min/kg) x mm^4)
RDUCT_MIN:   0.0060  : non-restrictive ductal resistance floor
KLVDIL   :   0.0035  : LV dilation gain (1/h)
KLVDILOUT:   0.010   : LV dilation regression (1/h)
PDIA0    :  26.0     : baseline diastolic blood pressure (mmHg)
KDIA     :   0.040   : diastolic pressure fall per mL/min/kg of runoff (mmHg)
QSIG     : 100.0     : shunt defining haemodynamic significance (mL/min/kg)

// ================= lung water and mechanics =================
KEVLW    :   0.0035  : lung water accrual per mL/min/kg of excess pulmonary flow
KEVLWOUT :   0.025   : lymphatic lung-water clearance (1/h)
EVLW0    :   6.0     : baseline extravascular lung water (mL/kg)
CRS0     :   0.75    : respiratory compliance at 26 wk (mL/cmH2O/kg)
CRS_GA   :   0.10    : compliance gain per week of GA
KCRS     :   0.030   : compliance equilibration (1/h)
CRS_EVLW :   0.030   : fractional compliance loss per mL/kg of lung water

// ================= kidney =================
GFR0     :   0.62    : GFR at 26 wk (mL/min/kg)
GFR_GA   :   0.115   : GFR gain per week of GA
KGFRMAT  :   0.0045  : postnatal GFR maturation (1/h)
GFR_PGE2 :   0.52    : fraction of GFR that is renal-prostanoid dependent
KPGE2K   :   0.25    : renal PGE2 turnover (1/h)
KI_IBU_K :   2.20    : unbound Ki, ibuprofen at renal COX (uM)
KI_IND_K :   0.0150  : unbound Ki, indomethacin at renal COX (uM)
IC50_APAP_K: 260.0   : unbound IC50, acetaminophen at renal COX (uM) - renal
                       peroxide tone is high, so the peroxidase mechanism is weak
CRPROD   :   0.00360 : creatinine production (mg/dL per h)
KCREL    :   0.01290 : creatinine elimination gain (per h per mL/min/kg)
SCR0     :   0.95    : serum creatinine at birth, maternal (mg/dL)
KUO      :   3.80    : urine output per mL/min/kg of GFR (mL/kg/h)
FLUIDIN  :   3.50    : total fluid intake (mL/kg/h)

// ================= gut, platelet, bilirubin =================
KPGE2G   :   0.20    : gut mucosal PGE2 turnover (1/h)
KI_IBU_G :   1.80    : unbound Ki, ibuprofen at gut COX (uM)
KI_IND_G :   0.0130  : unbound Ki, indomethacin at gut COX (uM)
IC50_APAP_G: 200.0   : unbound IC50, acetaminophen at gut COX (uM)
KTXA2    :   0.010   : platelet pool turnover (1/h)
KI_IBU_PLT:  1.10    : unbound Ki, ibuprofen at platelet COX-1 (uM)
KI_IND_PLT:  0.0060  : unbound Ki, indomethacin at platelet COX-1 (uM)
KBTIME   :   0.15    : bleeding-time index equilibration (1/h)
TBILI0   :   4.0     : baseline total bilirubin (mg/dL)
KBILIIN  :   0.16    : bilirubin production (mg/dL per h)
KBILIOUT :   0.030   : bilirubin conjugation + phototherapy (1/h)
BILIDISP_IBU: 0.024  : free-bilirubin rise per mg/L of total ibuprofen

// ================= non-COX vasoconstriction =================
// The single parameter that separates the two NSAIDs clinically.
// Indomethacin acutely reduces cerebral, mesenteric and renal flow
// velocity; ibuprofen essentially does not.
VC_IND   :   0.30    : maximal extra fractional flow reduction, indomethacin
VC_IBU   :   0.030   : maximal extra fractional flow reduction, ibuprofen
VC50_IND :   0.35    : total indomethacin for half-maximal effect (mg/L)
VC50_IBU :  30.0     : total ibuprofen for half-maximal effect (mg/L)

// ================= lung development =================
ALV0     :   1.0     : initial alveolarisation index
KALV     :   0.00060 : alveolarisation rate (1/h)
ALV_VENT :   0.55    : alveolar arrest per unit ventilator severity
ALV_EVLW :   0.018   : alveolar arrest per mL/kg of lung water

// ================= outcome hazards =================
// THE DRUG-HARM TERM, and why it has to exist.  The burden mechanism ALONE makes
// early ibuprofen beneficial: less shunt -> less burden -> lower BPD hazard.
// Baby-OSCAR found the opposite (69.2% vs 63.5%; death 13.6% vs 10.3%).  So
// either burden barely matters, or the drug harms something the burden term does
// not see.  This model takes the second option and names it: COX-2-derived
// prostaglandins are required for normal alveolar septation, and transient renal
// impairment adds fluid to an already wet lung.  B_BPD_COX is fitted to the arm
// SEPARATION while B_BUR is fixed a priori, so the model answers a question the
// trial cannot: given the burden reduction the drug actually achieves, how large
// must a competing harm be to cancel it?
B_BPD_COX:   0.5996      : BPD hazard, per unit cumulative COX exposure [FITTED]
H_BPD0   :   0.000130608 : baseline hazard, death or moderate/severe BPD (1/h)
B_GA     :   0.335    : BPD hazard, per week below 28 [FITTED]
B_BUR    :   0.00055  : BPD hazard, per unit PDA burden [FITTED]
B_VENT   :   0.052    : BPD hazard, per unit cumulative ventilator severity [FITTED]
H_NEC0   :   0.0000090: baseline hazard, NEC stage >= 2 (1/h)
B_NEC_MES:   3.60     : NEC hazard, per unit cumulative mesenteric deficit
B_NEC_PGE2G: 1.50     : NEC hazard, per unit gut PGE2 suppression
H_IVH0   :   0.000320 : baseline hazard, severe IVH (1/h)
B_IVH_GA :   0.40     : IVH hazard, per week below 28
B_IVH_BT :   0.55     : IVH hazard, per unit bleeding-time prolongation
B_IVH_CER:   1.10     : IVH hazard, per unit cerebral flow instability
IVH_IND_PROT: 0.42    : germinal-matrix maturation benefit of indomethacin
IVH_WINDOW: 96.0      : hours over which severe IVH essentially occurs
H_SIP0   :   0.0000075: baseline hazard, spontaneous intestinal perforation (1/h)
B_SIP_IND:   2.60     : SIP hazard, indomethacin effect
B_SIP_HC :   2.20     : SIP hazard, hydrocortisone effect
H_DTH0   :   0.0000260: baseline hazard, death before 36 wk PMA (1/h)
B_DTH_GA :   0.45     : death hazard, per week below 28
B_DTH_QS :   1.30     : death hazard, per unit systemic flow deficit

$CMT @annotated
IBU1    : ibuprofen central (mg/kg)
IBU2    : ibuprofen peripheral (mg/kg)
IBUE    : ibuprofen ductal effect site (mg/L)
IND1    : indomethacin central (mg/kg)
IND2    : indomethacin peripheral (mg/kg)
INDE    : indomethacin ductal effect site (mg/L)
APAPG   : acetaminophen enteral depot (mg/kg)
APAP1   : acetaminophen central (mg/kg)
APAP2   : acetaminophen peripheral (mg/kg)
GSH     : hepatic glutathione (fraction of normal)
ALT     : alanine aminotransferase (U/L)
PGE2D   : ductal tissue PGE2 (fraction of fetal baseline)
PGE2P   : circulating PGE2 (pg/mL)
PGDH    : pulmonary 15-PGDH capacity (fraction of term)
CAMP    : ductal smooth-muscle cAMP (fraction of fetal baseline)
NOD     : ductal NO tone (fraction of fetal baseline)
EP4     : ductal EP4 receptor density (fraction of fetal baseline)
PEROX   : ductal wall peroxide tone (fraction of normal)
CASMC   : ductal smooth-muscle cytosolic calcium (fraction of maximal)
TONE    : ductal smooth-muscle tone (0-1)
REMOD   : neointimal cushion / anatomic remodelling (0-1)
WALLTH  : ductal wall thickness index (fraction of term)
VEGFD   : ductal wall VEGF (fraction of normal)
TGFBD   : ductal wall TGF-beta1 (fraction of normal)
PVR     : pulmonary vascular resistance (mmHg/(mL/min/kg))
SVR     : systemic vascular resistance (mmHg/(mL/min/kg))
PLA     : left atrial pressure (mmHg)
LVDIL   : LV dilation / volume-load index
EVLW    : extravascular lung water (mL/kg)
CRS     : respiratory system compliance (mL/cmH2O/kg)
PDABUR  : PDA burden, integral of significant shunt ((mL/min/kg) x day)
MESDEF  : cumulative mesenteric flow deficit (fraction x day)
CERDEF  : cumulative cerebral flow deficit (fraction x day)
GFR     : glomerular filtration rate (mL/min/kg)
SCR     : serum creatinine (mg/dL)
UOCUM   : cumulative urine output (mL/kg)
FLUID   : net fluid balance (mL/kg)
PGE2K   : renal prostanoid tone (fraction of normal)
PGE2G   : gut mucosal PGE2 (fraction of normal)
TXA2    : platelet thromboxane A2 capacity (fraction of normal)
BTIME   : bleeding time index (fraction of normal)
TBILI   : total serum bilirubin (mg/dL)
BFREE   : free bilirubin index (fraction of normal)
ALV     : alveolarisation index
COXIX   : cumulative systemic COX-inhibition exposure (fraction x day)
VENTIX  : cumulative respiratory severity (index x day)
HBPD    : cumulative hazard, death or moderate/severe BPD
HNEC    : cumulative hazard, NEC stage >= 2
HIVH    : cumulative hazard, severe IVH
HSIP    : cumulative hazard, spontaneous intestinal perforation
HDEATH  : cumulative hazard, death before 36 wk PMA

$GLOBAL
// ---------------------------------------------------------------------
// Every algebraic quantity of the model lives in this ONE macro, which is
// invoked from both $ODE and $TABLE.  Writing the algebra twice -- once
// for the derivatives and once for the reported outputs -- is the classic
// way a model of this size silently drifts, so that what is reported is
// no longer what was integrated.  A macro cannot drift.
//
// WALL OXYGEN, the bistability switch (Kajino 2002, PMID 11809919).  The media
// has two O2 supplies with OPPOSITE wall-thickness dependence:
//   diffsup  DIFFUSION from lumen and adventitia -- adequate only across a THIN
//            wall, falling with the square of diffusion distance;
//   vasasup  VASA VASORUM -- present only in a THICK wall, and abolished by the
//            constriction itself.
// So a thin preterm wall stays oxygenated while constricted (no hypoxia -> no
// neointimal cushion -> the duct REOPENS when drug clears), whereas a thick term
// wall loses its vasa and goes profoundly hypoxic (permanent closure).  Because
// the wall also thickens postnatally, the same duct constricted at day 20
// remodels where it would not have at day 2 -- so the model PREDICTS that EARLY
// treatment reopens more often than late treatment, which is what is observed.
// An earlier form divided a SINGLE supply by thickness; it gave only a 1.8-fold
// preterm-to-term hypoxia ratio and produced no reopening at any gestation.
// ---------------------------------------------------------------------
#define PDA_ALG(TT)                                                        \\
  double PNAd = (TT) / 24.0;                                               \\
  /* ---- drug concentrations ---- */                                      \\
  double CIBU  = fmax(0.0, IBU1)  / V1_IBU;                                \\
  double CIND  = fmax(0.0, IND1)  / V1_IND;                                \\
  double CAPAP = fmax(0.0, APAP1) / V1_APAP;                               \\
  double UIBU  = CIBU  * FU_IBU  * 1000.0 / MW_IBU;                        \\
  double UIND  = CIND  * FU_IND  * 1000.0 / MW_IND;                        \\
  double UAPAP = CAPAP * FU_APAP * 1000.0 / MW_APAP;                       \\
  double UIBUE = fmax(0.0, IBUE) * FU_IBU * 1000.0 / MW_IBU;               \\
  double UINDE = fmax(0.0, INDE) * FU_IND * 1000.0 / MW_IND;               \\
  double perox = fmax(1e-6, PEROX);                                        \\
  /* ---- COX inhibition: channel (additive occupancy) x peroxidase ---- */\\
  double occD   = UIBUE / KI_IBU + UINDE / KI_IND;                         \\
  double ICHAN  = occD / (1.0 + occD);                                     \\
  double ic50D  = IC50_APAP * (1.0 + perox / KPEROX);                      \\
  double IPEROX = UAPAP / (UAPAP + ic50D);                                 \\
  double ICOXD  = 1.0 - (1.0 - ICHAN) * (1.0 - IPEROX);                    \\
  double occK   = UIBU / KI_IBU_K + UIND / KI_IND_K;                       \\
  double ic50K  = IC50_APAP_K * (1.0 + perox / KPEROX);                    \\
  double ICOXK  = 1.0 - (1.0 - occK / (1.0 + occK))                        \\
                      * (1.0 - UAPAP / (UAPAP + ic50K));                   \\
  double occG   = UIBU / KI_IBU_G + UIND / KI_IND_G;                       \\
  double ic50G  = IC50_APAP_G * (1.0 + perox / KPEROX);                    \\
  double ICOXG  = 1.0 - (1.0 - occG / (1.0 + occG))                        \\
                      * (1.0 - UAPAP / (UAPAP + ic50G));                   \\
  double occP   = UIBU / KI_IBU_PLT + UIND / KI_IND_PLT;                   \\
  double ICOXPLT = occP / (1.0 + occP);                                    \\
  /* ---- ductal relaxant drive ---- */                                    \\
  double drive = fmax(0.0, PGE2D) * fmax(0.0, EP4);                        \\
  double RELAXP = EMAX_PGE2 * drive / (EC50_PGE2 + drive);                 \\
  double nod = fmax(0.0, NOD);                                             \\
  double RELAXN = EMAX_NO * nod / (EC50_NO + nod);                         \\
  /* ---- O2-driven constriction ---- */                                   \\
  double P50 = P50_O2 * exp(P50_GA * (28.0 - GA));                         \\
  double O2GAIN = PAO2 * PAO2 / (PAO2 * PAO2 + P50 * P50);                 \\
  /* ---- maximal achievable occlusion, and the net-drive sigmoid ---- */  \\
  double TMAXGA = 1.0 / (1.0 + exp(-(GA - TMAX_GA50) * TMAX_SL));          \\
  double FINV = exp(KINVGA * (GA - 26.0));                                 \\
  double MATC = 1.0 - exp(-(TT) * FINV / TAU_MAT);                         \\
  double NETD = GO2 * O2GAIN + GMAT * MATC                                 \\
              - GPGE2 * RELAXP - GNO * RELAXN;                             \\
  double TONETGT = TMAXGA / (1.0 + exp(-(NETD - NET50) / NETW));           \\
  TONETGT = fmin(1.0, fmax(0.0, TONETGT));                                 \\
  /* ---- ductal geometry ---- */                                          \\
  double DMAX = fmax(0.6, DMAX_INT + DMAX_GA * GA);                        \\
  double tone = fmin(1.0, fmax(0.0, TONE));                                \\
  double remod = fmin(1.0, fmax(0.0, REMOD));                              \\
  double DDUCT = DMAX * (1.0 - remod) * (1.0 - ALPHA_D * tone);            \\
  double CLOSED = (DDUCT < DCLOSE) ? 1.0 : 0.0;                            \\
  /* ---- wall O2: the bistability switch ---- */                          \\
  double wall = fmax(0.05, WALLTH);                                        \\
  double wallref = fmax(0.05, WALLTH_INT + WALLTH_GA * GA);                \\
  double lumfrac = DDUCT / DMAX;                                          \\
  double wr = wallref * wall;                                             \\
  double diffsup = FLUM_O2 * (0.35 + 0.65 * lumfrac)                      \\
                 / (1.0 + (KDIFF * wr) * (KDIFF * wr));                   \\
  double vasasup = FVASA_O2 * wr * (1.0 - fmin(1.0, TONE));               \\
  double WALLO2 = PAO2 * (diffsup + vasasup);                             \\
  double HYPOX = fmax(0.0, O2CRIT - WALLO2);                               \\
  /* ---- shunt hydraulics and the transitional circulation ---- */        \\
  double dd = fmax(0.05, DDUCT);                                           \\
  double RDUCT = KRDUCT / (dd * dd * dd * dd) + RDUCT_MIN;                 \\
  double pvr = fmax(0.004, PVR);                                           \\
  double svr = fmax(0.05, SVR);                                            \\
  double pla = fmax(0.0, PLA);                                             \\
  double qmax = QMAX_LV + QMAX_GA * (GA - 26.0);                           \\
  double den = RDUCT + pvr;                                                \\
  double QS = QS_TARGET;                                                   \\
  double QSH = fmax(0.0, (QS * (svr - pvr) - pla) / den);                  \\
  double QP = QS + QSH;                                                    \\
  if (QP > qmax) {   /* LV reserve exhausted -> systemic steal */          \\
    QS = (qmax + pla / den) / (1.0 + (svr - pvr) / den);                   \\
    QS = fmax(40.0, QS);                                                   \\
    QSH = fmax(0.0, (QS * (svr - pvr) - pla) / den);                       \\
    QP = QS + QSH;                                                         \\
  }                                                                        \\
  double PAOM = QS * svr;                                                  \\
  double PPA = QP * pvr + pla;                                             \\
  double QPQS = QP / fmax(1.0, QS);                                        \\
  double SHFRAC = QSH / fmax(1.0, QP);                                     \\
  double QSDEF = fmax(0.0, (QS_TARGET - QS) / QS_TARGET);                  \\
  /* ---- diastolic runoff and organ perfusion ---- */                     \\
  double PDIA = fmax(8.0, PDIA0 + 0.20 * (GA - 26.0) - KDIA * QSH);        \\
  double vcind = VC_IND * CIND / (CIND + VC50_IND);                        \\
  double vcibu = VC_IBU * CIBU / (CIBU + VC50_IBU);                        \\
  double VCDRUG = vcind + vcibu;                                           \\
  double perf = (PDIA / PDIA0) * (1.0 - VCDRUG);                           \\
  double QCERREL = fmax(0.05, perf);                                       \\
  double QMESREL = fmax(0.05, perf * (1.0 - 0.25 * SHFRAC));               \\
  double QRENREL = fmax(0.05, perf);                                       \\
  /* ---- lung mechanics and ventilator exposure ---- */                   \\
  double evlw = fmax(0.0, EVLW);                                           \\
  double crsref = CRS0 + CRS_GA * (GA - 26.0);                             \\
  double CRSREL = fmax(0.05, CRS) / crsref;                                \\
  double VENTSEV = fmax(0.0, (1.0 / fmax(0.15, CRSREL) - 1.0)              \\
                     + 0.030 * fmax(0.0, evlw - EVLW0));                   \\
  double FIO2 = fmin(1.0, 0.21 + 0.42 * VENTSEV);                          \\
  double UO = KUO * fmax(0.02, GFR);                                       \\
  double NAPQI = FNAPQI * CL_APAP0 * CAPAP;

$MAIN
GSH_0    = 1.0;
ALT_0    = ALT0;
PGE2D_0  = 1.0;
PGE2P_0  = PGE2P0;
PGDH_0   = fmax(0.02, fmin(1.0, PGDH_SL * (GA - PGDH0_GA)));
CAMP_0   = 1.0;
NOD_0    = exp(-NOGA * (GA - 24.0));
EP4_0    = exp(-EP4GA * (GA - 24.0));
PEROX_0  = 1.0 + (PEROXSEP - 1.0) * SEPSIS;
CASMC_0  = 0.05;
TONE_0   = 0.02;
REMOD_0  = 0.0;
WALLTH_0 = 0.35;
VEGFD_0  = 0.0;
TGFBD_0  = 0.0;
PVR_0    = PVR0;
SVR_0    = SVR0;
PLA_0    = PLA0;
LVDIL_0  = 0.0;
EVLW_0   = EVLW0;
CRS_0    = (CRS0 + CRS_GA * (GA - 26.0)) * (1.0 + 0.10 * ANTESTER);
GFR_0    = GFR0 + GFR_GA * (GA - 26.0);
SCR_0    = SCR0;
PGE2K_0  = 1.0;
PGE2G_0  = 1.0;
TXA2_0   = 1.0;
BTIME_0  = 1.0;
TBILI_0  = TBILI0;
BFREE_0  = 1.0;
ALV_0    = ALV0;

$ODE
PDA_ALG(SOLVERTIME)

// ---------------- drug PK ----------------
double clibu  = CL_IBU0  * (1.0 + CLMAT_IBU  * PNAd);
double clind  = CL_IND0  * (1.0 + CLMAT_IND  * PNAd);
double clapap = CL_APAP0 * (1.0 + CLMAT_APAP * PNAd);

double c1i = IBU1 / V1_IBU, c2i = IBU2 / V2_IBU;
dxdt_IBU1 = -clibu * c1i - Q_IBU * (c1i - c2i);
dxdt_IBU2 =  Q_IBU * (c1i - c2i);
dxdt_IBUE =  KE0_IBU * (CIBU - IBUE);

double c1n = IND1 / V1_IND, c2n = IND2 / V2_IND;
dxdt_IND1 = -clind * c1n - Q_IND * (c1n - c2n);
dxdt_IND2 =  Q_IND * (c1n - c2n);
dxdt_INDE =  KE0_IND * (CIND - INDE);

double c1a = APAP1 / V1_APAP, c2a = APAP2 / V2_APAP;
dxdt_APAPG = -KA_APAP * APAPG;
dxdt_APAP1 =  F_APAP * KA_APAP * APAPG - clapap * c1a - Q_APAP * (c1a - c2a);
dxdt_APAP2 =  Q_APAP * (c1a - c2a);

// ---------------- acetaminophen hepatic safety ----------------
double gsh = fmax(0.0, GSH);
double quench = NAPQI * gsh;
dxdt_GSH = KGSH * (GSHCAP - gsh) - 0.55 * quench;
dxdt_ALT = KALT * NAPQI * (1.0 - gsh) * 100.0 - KALTOUT * (ALT - ALT0);

// ---------------- prostanoid system ----------------
double pgdh = fmax(0.0, PGDH);
double qpfrac = QP / 200.0;
dxdt_PGE2P = RPROD_PGE2P * (1.0 - ICOXG)
           - (KPUL_PGE2 * pgdh * qpfrac + KOTH_PGE2) * PGE2P;
dxdt_PGDH = KPGDH * (1.0 - pgdh);

// FSYN is the postnatal involution of local ductal prostanoid synthesis,
// faster at higher GA.  This term is what makes spontaneous closure happen,
// and hence what bounds every drug effect on a fixed-time endpoint.
double tausyn = TAUSYN0 / exp(KTAUSYN * (GA - 26.0));
double fsyn = FSYN_INF + (1.0 - FSYN_INF) * exp(-SOLVERTIME / tausyn);
dxdt_PGE2D = KSYN_PGE2D * fsyn * (1.0 - ICOXD) - KDEG_PGE2D * PGE2D
           + KIN_PGE2D * (PGE2P / PGE2P0) * (1.0 - ICOXD);

dxdt_CAMP  = KCAMP * (RELAXP / fmax(1e-6, EMAX_PGE2) - CAMP);
dxdt_CASMC = 2.0 * (TONE - CASMC);

double notgt = exp(-NOGA * (GA - 24.0)) * exp(-KNODEC * FINV * SOLVERTIME);
dxdt_NOD = KNO * (notgt - NOD);

double ep4tgt = exp(-EP4GA * (GA - 24.0));
dxdt_EP4 = -KEP4 * FINV * (EP4 - EP4INF * ep4tgt);

double peroxtgt = 1.0 + (PEROXSEP - 1.0) * SEPSIS;
dxdt_PEROX = KPEROXIN * (peroxtgt - PEROX);

// ---------------- ductal contraction ----------------
dxdt_TONE = KTONE * (TONETGT - TONE);

// ---------------- wall O2 and remodelling ----------------
dxdt_WALLTH = KWALL * (1.0 - wall);
dxdt_VEGFD  = KVEGF * HYPOX - KVEGFOUT * VEGFD;
dxdt_TGFBD  = KTGFB * VEGFD - KTGFBOUT * TGFBD;
double frem = 1.0 / (1.0 + exp(-(GA - KREMGA50) * KREMGA_SL));
dxdt_REMOD  = KREMOD * frem * TGFBD * (1.0 - remod) - KUNREMOD * remod;

// ---------------- haemodynamics ----------------
double pvrtgt = PVRINF + PVR_EVLW * fmax(0.0, EVLW - EVLW0);
dxdt_PVR = -KPVR * (PVR - pvrtgt);
dxdt_SVR = KSVR * (SVRINF - SVR);
dxdt_PLA = KPLA * (PLA0 + PLA_Q * QSH - PLA);
dxdt_LVDIL = KLVDIL * SHFRAC * 10.0 - KLVDILOUT * LVDIL;

// ---------------- lung water and mechanics ----------------
double excessqp = fmax(0.0, QP - 200.0);
dxdt_EVLW = KEVLW * excessqp * (1.0 + 0.5 * PLA / 4.0)
          - KEVLWOUT * (EVLW - EVLW0);
double crstgt = crsref * (1.0 - fmin(0.80, CRS_EVLW * fmax(0.0, EVLW - EVLW0)))
              * (1.0 + 0.10 * ANTESTER);
dxdt_CRS = KCRS * (crstgt - CRS);

// ---------------- burden: the variable the outcomes actually see ----------
dxdt_PDABUR = fmax(0.0, QSH - QSIG) / 24.0;
dxdt_MESDEF = fmax(0.0, 1.0 - QMESREL) / 24.0;
dxdt_CERDEF = fmax(0.0, 1.0 - QCERREL) / 24.0;

// ---------------- kidney ----------------
dxdt_PGE2K = KPGE2K * ((1.0 - ICOXK) - PGE2K);
double gfrref = (GFR0 + GFR_GA * (GA - 26.0)) * (1.0 + KGFRMAT * SOLVERTIME);
double gfrtgt = gfrref * ((1.0 - GFR_PGE2) + GFR_PGE2 * fmax(0.0, PGE2K))
              * (0.55 + 0.45 * PDIA / PDIA0) * (1.0 - 0.5 * VCDRUG);
dxdt_GFR = 0.25 * (gfrtgt - GFR);
dxdt_SCR = CRPROD - KCREL * GFR * SCR;
dxdt_UOCUM = UO;
dxdt_FLUID = FLUIDIN - UO - 2.2;

// ---------------- gut, platelet, bilirubin ----------------
dxdt_PGE2G = KPGE2G * ((1.0 - ICOXG) - PGE2G);
dxdt_TXA2  = KTXA2 * ((1.0 - ICOXPLT) - TXA2);
double bttgt = 1.0 + 1.9 * (1.0 - fmax(0.0, TXA2));
dxdt_BTIME = KBTIME * (bttgt - BTIME);
dxdt_TBILI = KBILIIN - KBILIOUT * TBILI;
double bftgt = (TBILI / TBILI0) * (1.0 + BILIDISP_IBU * CIBU);
dxdt_BFREE = 0.5 * (bftgt - BFREE);

// ---------------- lung development ----------------
double arrest = ALV_VENT * VENTSEV + ALV_EVLW * fmax(0.0, EVLW - EVLW0);
dxdt_ALV = KALV * fmax(0.0, 1.0 - arrest);
dxdt_VENTIX = VENTSEV / 24.0;

// ---------------- outcome hazards ----------------
double gapen = exp(B_GA * fmax(0.0, 28.0 - GA));
dxdt_COXIX = ICOXG / 24.0;
dxdt_HBPD = H_BPD0 * gapen * (1.0 + B_BUR * PDABUR)
          * (1.0 + B_VENT * VENTIX) * (1.0 + B_BPD_COX * COXIX)
          * (1.0 - 0.12 * ANTESTER);

dxdt_HNEC = H_NEC0 * exp(0.30 * fmax(0.0, 28.0 - GA))
          * (1.0 + B_NEC_MES * MESDEF)
          * (1.0 + B_NEC_PGE2G * (1.0 - fmax(0.0, PGE2G)));

// Severe IVH is a first-96-hour phenomenon with TWO OPPOSING drug terms:
// platelet COX-1 inhibition prolongs bleeding time (harm) while
// indomethacin matures germinal-matrix vessels (benefit, TIPP).  They
// nearly cancel, which is the model reason prophylactic indomethacin
// reduces IVH without improving long-term outcome.
double ivhgate = (SOLVERTIME < IVH_WINDOW) ? 1.0
               : exp(-(SOLVERTIME - IVH_WINDOW) / 48.0);
double indprot = 1.0 - IVH_IND_PROT * (INDE / (INDE + 0.25));
dxdt_HIVH = H_IVH0 * ivhgate * exp(B_IVH_GA * fmax(0.0, 28.0 - GA))
          * (1.0 + B_IVH_BT * fmax(0.0, BTIME - 1.0))
          * (1.0 + B_IVH_CER * QSDEF) * fmax(0.15, indprot);

dxdt_HSIP = H_SIP0 * exp(0.30 * fmax(0.0, 28.0 - GA))
          * (1.0 + B_SIP_IND * (INDE / (INDE + 0.3)))
          * (1.0 + B_SIP_HC * HCORT);

dxdt_HDEATH = H_DTH0 * exp(B_DTH_GA * fmax(0.0, 28.0 - GA))
            * (1.0 + B_DTH_QS * QSDEF) * (1.0 + 0.6 * VENTSEV);

$TABLE
PDA_ALG(TIME)
double DAY   = TIME / 24.0;
double PBPD  = 1.0 - exp(-HBPD);
double PNEC  = 1.0 - exp(-HNEC);
double PIVH  = 1.0 - exp(-HIVH);
double PSIP  = 1.0 - exp(-HSIP);
double PDTH  = 1.0 - exp(-HDEATH);
// death OR moderate/severe BPD, combining independent hazards -- the
// primary endpoint of Baby-OSCAR and the composite of BeNeDuctus
double PCOMP  = 1.0 - (1.0 - PBPD) * (1.0 - PDTH);
double PCOMP3 = 1.0 - (1.0 - PBPD) * (1.0 - PDTH) * (1.0 - PNEC);

$CAPTURE @annotated
DAY     : postnatal age (days)
CIBU    : total ibuprofen concentration (mg/L)
CIND    : total indomethacin concentration (mg/L)
CAPAP   : total acetaminophen concentration (mg/L)
UIBU    : unbound ibuprofen (uM)
UIND    : unbound indomethacin (uM)
UAPAP   : unbound acetaminophen (uM)
ICOXD   : net ductal COX inhibition (fraction)
ICHAN   : arachidonate-channel component of ductal COX inhibition
IPEROX  : peroxidase-site component of ductal COX inhibition
ICOXK   : renal COX inhibition (fraction)
ICOXG   : gut mucosal COX inhibition (fraction)
ICOXPLT : platelet COX-1 inhibition (fraction)
RELAXP  : PGE2/EP4-mediated relaxation (fraction)
RELAXN  : NO-mediated relaxation (fraction)
O2GAIN  : O2-driven constrictor gain (fraction)
NETD    : net constrictor drive
TONETGT : target ductal tone
TMAXGA  : maximal achievable occlusion at this gestational age
DMAX    : unconstricted ductal diameter (mm)
DDUCT   : ductal diameter (mm)
CLOSED  : ductal closure indicator (diameter < DCLOSE)
WALLO2  : ductal wall PO2 (mmHg) - the remodelling switch
HYPOX   : wall hypoxic signal (mmHg below threshold)
RDUCT   : ductal resistance (mmHg/(mL/min/kg))
QSH     : left-to-right ductal shunt (mL/min/kg)
QP      : pulmonary blood flow (mL/min/kg)
QS      : systemic blood flow (mL/min/kg)
QPQS    : pulmonary-to-systemic flow ratio
SHFRAC  : shunt fraction of pulmonary flow
QSDEF   : systemic flow deficit (fraction)
PAOM    : mean arterial pressure (mmHg)
PPA     : mean pulmonary artery pressure (mmHg)
PDIA    : diastolic blood pressure (mmHg)
QCERREL : cerebral blood flow, fraction of baseline
QMESREL : mesenteric blood flow, fraction of baseline
QRENREL : renal blood flow, fraction of baseline
VCDRUG  : non-COX drug vasoconstriction (fraction)
CRSREL  : respiratory compliance, fraction of gestation-matched
VENTSEV : ventilator severity index
FIO2    : inspired oxygen fraction
UO      : urine output (mL/kg/h)
PBPD    : probability of moderate/severe BPD
PNEC    : probability of NEC stage >= 2
PIVH    : probability of severe IVH
PSIP    : probability of spontaneous intestinal perforation
PDTH    : probability of death before 36 wk PMA
PCOMP   : probability of death or moderate/severe BPD
PCOMP3  : probability of death, moderate/severe BPD or NEC
'

mod <- mcode("pda_qsp", pda_code, soloc = tempdir())

## =====================================================================
##  DOSING REGIMENS
## =====================================================================

## Combine event objects at ABSOLUTE times.  Deliberately not using `+` on ev
## objects: `ev_seq` offsets sequentially and the `+` semantics are easy to get
## wrong for co-administration, where both regimens must keep their own clock.
## Binding the underlying data frames makes the intent unambiguous.
ev_bind <- function(...) {
  ls <- Filter(Negate(is.null), list(...))
  if (!length(ls)) return(NULL)
  d <- do.call(rbind, lapply(ls, as.data.frame))
  d <- d[order(d$time), , drop = FALSE]
  as.ev(d)
}

rx_none <- function() NULL

## IV ibuprofen lysine: 10-5-5 mg/kg q24h (standard) or 20-10-10 (high)
rx_ibuprofen <- function(start_h = 48, high = FALSE) {
  d <- if (high) c(20, 10, 10) else c(10, 5, 5)
  ev(time = start_h + c(0, 24, 48), amt = d, cmt = "IBU1")
}

## Continuous infusion: 20 mg/kg over 24 h, then 10 mg/kg/24 h x 2.
## Maintaining the trough is the point -- see the scenario table.
rx_ibuprofen_inf <- function(start_h = 48) {
  ev(time  = start_h + c(0, 24, 48),
     amt   = c(20, 10, 10),
     rate  = c(20/24, 10/24, 10/24),
     cmt   = "IBU1")
}

rx_indomethacin <- function(start_h = 48, dose = 0.2, maint = 0.1, n = 3,
                            interval = 24) {
  ev(time = start_h + interval * (0:(n - 1)),
     amt  = c(dose, rep(maint, n - 1)), cmt = "IND1")
}

rx_acetaminophen <- function(start_h = 48, dose = 15, interval = 6,
                             days = 3, iv = TRUE) {
  n <- days * 24 / interval
  ev(time = start_h + interval * (0:(n - 1)), amt = dose,
     cmt = if (iv) "APAP1" else "APAPG")
}

## =====================================================================
##  SIMULATION HELPERS
## =====================================================================

sim_pda <- function(doses = NULL, days = 90, delta = 1, ...) {
  p <- list(...)
  m <- mod
  if (length(p)) m <- do.call(param, c(list(m), p))
  if (is.null(doses)) {
    m %>% mrgsim(end = days * 24, delta = delta, hmax = 1.0,
                 atol = 1e-9, rtol = 1e-7)
  } else {
    m %>% ev(doses) %>% mrgsim(end = days * 24, delta = delta, hmax = 1.0,
                               atol = 1e-9, rtol = 1e-7)
  }
}

## first postnatal day on which the duct is closed
closure_day <- function(out) {
  d <- as.data.frame(out)
  i <- which(d$DDUCT < 0.30)
  if (!length(i)) NA_real_ else d$DAY[i[1]]
}

## reached closure, then the diameter recovered -- the reopening phenotype
reopened <- function(out, factor = 1.6) {
  d <- as.data.frame(out)
  i <- which(d$DDUCT < 0.30)
  if (!length(i)) return(FALSE)
  any(d$DDUCT[i[1]:nrow(d)] > 0.30 * factor)
}

val_at <- function(out, name, day) {
  d <- as.data.frame(out)
  d[[name]][which.min(abs(d$DAY - day))]
}

## =====================================================================
##  SIXTEEN TREATMENT SCENARIOS
## =====================================================================

SCENARIOS <- list(
  list(id = "S1",  lab = "Expectant management, 26 wk",
       par = list(GA = 26), rx = NULL),
  list(id = "S2",  lab = "Early ibuprofen d2 (Baby-OSCAR design), 26 wk",
       par = list(GA = 26), rx = rx_ibuprofen(48)),
  list(id = "S3",  lab = "Late ibuprofen d7, standard dose, 26 wk",
       par = list(GA = 26), rx = rx_ibuprofen(168)),
  list(id = "S4",  lab = "Late ibuprofen d7, HIGH dose 20-10-10, 26 wk",
       par = list(GA = 26), rx = rx_ibuprofen(168, high = TRUE)),
  list(id = "S5",  lab = "Continuous ibuprofen infusion d2, 26 wk",
       par = list(GA = 26), rx = rx_ibuprofen_inf(48)),
  list(id = "S6",  lab = "Indomethacin 0.2-0.1-0.1 d2, 26 wk",
       par = list(GA = 26), rx = rx_indomethacin(48)),
  list(id = "S7",  lab = "Prophylactic indomethacin 0.1 q24h x3 from 8 h (TIPP)",
       par = list(GA = 26), rx = rx_indomethacin(8, 0.1, 0.1, 3)),
  list(id = "S8",  lab = "IV acetaminophen 15 mg/kg q6h x 3 d from d2, 26 wk",
       par = list(GA = 26), rx = rx_acetaminophen(48)),
  list(id = "S9",  lab = "Acetaminophen d2 WITH chorioamnionitis (high peroxide)",
       par = list(GA = 26, SEPSIS = 1), rx = rx_acetaminophen(48)),
  list(id = "S10", lab = "Ibuprofen d2 WITH chorioamnionitis (control for S9)",
       par = list(GA = 26, SEPSIS = 1), rx = rx_ibuprofen(48)),
  list(id = "S11", lab = "Ibuprofen + acetaminophen combination d2",
       par = list(GA = 26),
       rx = ev_bind(rx_ibuprofen(48), rx_acetaminophen(48))),
  list(id = "S12", lab = "Ibuprofen d2, second course d7 after reopening",
       par = list(GA = 26),
       rx = ev_bind(rx_ibuprofen(48), rx_ibuprofen(168))),
  list(id = "S13", lab = "Indomethacin d2 + early hydrocortisone (SIP risk)",
       par = list(GA = 26, HCORT = 1), rx = rx_indomethacin(48)),
  list(id = "S14", lab = "Ibuprofen d2, extreme prematurity 24 wk",
       par = list(GA = 24, BW = 0.62), rx = rx_ibuprofen(48)),
  list(id = "S15", lab = "Ibuprofen d2, 29 wk",
       par = list(GA = 29, BW = 1.20), rx = rx_ibuprofen(48)),
  list(id = "S16", lab = "Targeted: high-dose ibuprofen at d10 only (echo-guided)",
       par = list(GA = 26), rx = rx_ibuprofen(240, high = TRUE))
)

run_scenarios <- function(days = 90) {
  rows <- lapply(SCENARIOS, function(s) {
    out <- do.call(sim_pda, c(list(doses = s$rx, days = days), s$par))
    d <- as.data.frame(out)
    early <- d[d$DAY <= 10, ]
    data.frame(
      id        = s$id,
      scenario  = s$lab,
      close_day = closure_day(out),
      reopened  = reopened(out),
      burden    = round(val_at(out, "PDABUR", 36), 0),
      peak_ICOX = round(100 * max(d$ICOXD), 1),
      dmin      = round(min(d$DDUCT), 2),
      SCr_d5    = round(val_at(out, "SCR", 5), 2),
      UO_nadir  = round(min(early$UO), 2),
      QCER_nadir= round(100 * min(early$QCERREL), 0),
      P_BPD     = round(100 * val_at(out, "PBPD", 70), 1),
      P_NEC     = round(100 * val_at(out, "PNEC", 70), 1),
      P_IVH     = round(100 * val_at(out, "PIVH", 70), 1),
      P_comp    = round(100 * val_at(out, "PCOMP", 70), 1),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

## ---------------------------------------------------------------------
##  Spontaneous closure by gestational age -- the calibration anchor.
##  Target (Semberova 2017; Nemerofsky 2008; Rolland 2015):
##  ~85 d at 24 wk, ~71 d at 25 wk, ~48 d at 26 wk, ~21 d at 28 wk,
##  ~10 d at 30 wk, ~5 d at 32 wk, ~1.5 d at term.
## ---------------------------------------------------------------------
spontaneous_closure_table <- function(gas = c(24, 25, 26, 27, 28, 30, 32, 38)) {
  do.call(rbind, lapply(gas, function(g) {
    out <- sim_pda(NULL, days = 150, delta = 3, GA = g,
                   BW = 0.5 + 0.11 * (g - 23))
    data.frame(GA = g, closure_day = round(closure_day(out), 1),
               d_day3 = round(val_at(out, "DDUCT", 3), 2),
               QSH_day3 = round(val_at(out, "QSH", 3), 0),
               QpQs_day3 = round(val_at(out, "QPQS", 3), 2))
  }))
}

## ---------------------------------------------------------------------
##  Clearance maturation: the SAME regimen is a smaller exposure at d7
##  than at d2, because ibuprofen clearance nearly doubles over the first
##  week.  This is the model's account of why late rescue needs a higher
##  dose.
## ---------------------------------------------------------------------
maturation_table <- function() {
  a <- sim_pda(rx_ibuprofen(48),  days = 14, delta = 0.25)
  b <- sim_pda(rx_ibuprofen(168), days = 18, delta = 0.25)
  c <- sim_pda(rx_ibuprofen(168, high = TRUE), days = 18, delta = 0.25)
  auc <- function(out, d0, d1) {
    d <- as.data.frame(out); m <- d$DAY >= d0 & d$DAY <= d1
    sum(diff(d$DAY[m]) * (head(d$CIBU[m], -1) + tail(d$CIBU[m], -1)) / 2)
  }
  data.frame(
    regimen = c("10-5-5 from d2", "10-5-5 from d7", "20-10-10 from d7"),
    AUC_6d  = round(c(auc(a, 2, 8), auc(b, 7, 13), auc(c, 7, 13)), 0),
    peak_ICOXD = round(100 * c(max(as.data.frame(a)$ICOXD),
                              max(as.data.frame(b)$ICOXD),
                              max(as.data.frame(c)$ICOXD)), 1))
}

## ---------------------------------------------------------------------
##  Wall-O2 bistability: why preterm ducts reopen and term ducts do not.
##  Nothing here was fitted to reopening rates.
## ---------------------------------------------------------------------
reopening_table <- function(gas = c(24, 26, 28, 30, 34, 38)) {
  do.call(rbind, lapply(gas, function(g) {
    out <- sim_pda(rx_ibuprofen(48), days = 40, GA = g,
                   BW = 0.5 + 0.11 * (g - 23))
    d <- as.data.frame(out)
    i <- which.max(d$TONE)
    data.frame(GA = g,
               wall_index = round(-1.20 + 0.075 * g, 2),
               wallPO2_at_peak_tone = round(d$WALLO2[i], 1),
               neointima_d20 = round(val_at(out, "REMOD", 20), 2),
               reopened = reopened(out))
  }))
}

## Run the batch report ONLY when this file is executed directly (Rscript).
## sys.nframe() > 0 when the file is sourced -- e.g. by pda_shiny_app.R -- so the
## app does not silently trigger 16 scenario simulations at start-up.
if (sys.nframe() == 0L) {
  cat("\n=== spontaneous closure by gestational age ===\n")
  print(spontaneous_closure_table())
  cat("\n=== ibuprofen clearance maturation ===\n")
  print(maturation_table())
  cat("\n=== 16 treatment scenarios ===\n")
  print(run_scenarios(), row.names = FALSE)
  cat("\n=== wall-O2 bistability / reopening (PREDICTION) ===\n")
  print(reopening_table())
}
