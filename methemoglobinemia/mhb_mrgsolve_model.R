# =============================================================================
#  mhb_mrgsolve_model.R
#  Methaemoglobinaemia — QSP model (mrgsolve)
#  메트헤모글로빈혈증 — 정량적 시스템 약리학 모델
# =============================================================================
#
#  THE CLAIM THIS MODEL IS BUILT TO MAKE
#  -------------------------------------
#  Methaemoglobinaemia is a DELIVERY disease, not a SATURATION disease, and the
#  two numbers that reach the bedside are each wrong in an opposite direction.
#
#    * %MetHb UNDERSTATES the injury.  Ferric subunits do not merely subtract
#      carrying capacity; they LEFT-SHIFT the subunits beside them (Darling &
#      Roughton 1942) and abolish cooperativity.  Oxygen that is still carried
#      is no longer released.  The model computes that "MetHb 30%" is
#      haemodynamically equivalent to an anaemia of 8.3 g/dL, not the 10.5 g/dL
#      that 15 x 0.70 suggests -- and the error GROWS with severity.
#
#    * SpO2 understates it too, but for a reason that has nothing to do with
#      the patient.  The oximeter forms ONE ratio R = A660/A940 and reads it
#      off ONE calibration line.  Any pigment that absorbs equally at both
#      wavelengths drags R towards 1, and the line evaluated at R = 1 is
#      110 - 25 = 85.  The 85% floor is a property of the DEVICE.
#
#  And the antidote is the same molecule as the poison, separated only by an
#  electron supply (NADPH) that the disease's most important comorbidity
#  (G6PD deficiency) destroys.
#
#  WHAT IS *NOT* IN THE PARAMETER LIST
#  -----------------------------------
#  There is no severity scale, no "7 mg/kg" ceiling parameter, no G6PD
#  contraindication switch, no rebound term, no rule that sulfhaemoglobin or
#  HbM fail to respond to methylene blue, and no 85% floor constant.  All of
#  those are OUTPUTS.  What is in the parameter list is:
#
#      - a catalytic co-oxidation cycle whose co-substrate is glutathione
#      - one capped NADPH supply with two named consumers
#      - a leucomethylene-blue branch point with a saturable good branch and a
#        non-saturable futile branch
#      - an allosteric term coupling ferric fraction to P50 and to Hill n
#      - two absorbances per pigment and one straight calibration line
#
#  VERIFICATION
#  ------------
#  R was not available in the container in which this file was written, so it
#  is EQUATION-VERIFIED rather than compile-verified.  Every ODE, parameter and
#  output below is independently re-implemented in Python/scipy in
#  `mhb_reference_check.py`; run
#
#      python3 mhb_reference_check.py --params mhb_mrgsolve_model.R
#
#  to cross-check every parameter in this file against that one, and
#
#      python3 mhb_reference_check.py --all
#
#  to reproduce every number quoted in README.md.  Where the two disagree, the
#  Python file is the one that was actually executed.
#
#  UNITS
#  -----
#    time                  h
#    haemoglobin species   g/dL of WHOLE BLOOD
#    red-cell solutes      uM referred to RED CELL WATER
#    plasma solutes        uM, or mg/L for dapsone / benzocaine / cimetidine
#    methylene blue        umol (amounts), converted to concentrations in $ODE
#
#    1 g/dL Hb == 2482.3 uM haem per litre of blood      (HEMEC)
#    a flux of 1 uM(RBC)/h of haem == 1.8129e-4 g/dL/h   (CONV = HCT/HEMEC)
#
#  USAGE
#  -----
#    library(mrgsolve); library(dplyr); library(ggplot2)
#    mod <- mread("mhb_mrgsolve_model", "mhb_mrgsolve_model.R")
#    out <- mod %>% ev(sc_benzocaine_mb()) %>% mrgsim(end = 48, delta = 0.05)
#    plot(out, MetPct + SPO2 + PVO2 + HBTOT ~ time)
#
# =============================================================================

$PROB
# Methaemoglobinaemia QSP model — 42 ODE states, 139 parameters, 23 scenarios

# ---------------------------------------------------------------------------
# This block is generated from, and cross-checked against, the parameter
# dictionary `P` in mhb_reference_check.py.  Do not edit one without the other.
# ---------------------------------------------------------------------------
$PARAM
  WT = 70, HB0 = 15, FMET0 = 0.008, PAO2 = 95,
  FIO2FLAG = 0, VO2 = 250, CO0 = 5, P50 = 26.8,
  NHILL = 2.7, ALPHAM = 0.5, BBPG = 0.55, PVCRIT = 20,
  GCO = 1.6, TAUCO = 0.3, COMAX = 2.6, GBPG = 1.1,
  TAUBPG = 30, KAUTO = 0.00118, KNOH = 0.042, KTOL = 0.011,
  KNIT = 0.0012, KNIT2 = 0.0006, KROSOX = 0.0022, KMBOXH = 0.0016,
  KSULF = 4e-05, VMXB5 = 0.675, KMB5 = 4.5, EB5SET = 1,
  TAUEB5 = 48, GRIBO = 0.35, FHBM = 0, KASC = 1.05e-05,
  VPPP = 40000, G6PD = 1, VFRMAX = 45000, KMFR = 12,
  VGRMAX = 60000, KMGR = 60, VC = 40, VP = 1360,
  CLMB = 260, QMB = 300, KUP = 25, KOUT = 5,
  MWMB = 319.85, KMBR = 1700, KMHB = 0.8, KO2 = 95,
  KMBLK = 0.6, GSH0 = 2000, KSYN = 0.055, KGSHOX = 30,
  KG3 = 350, KROSB = 0.55, KCAT = 9, WGPX = 1.6,
  KRECYC = 520, KNOARD = 0.62, KTOLND = 0.9, KGR2 = 600,
  KRCROS = 0.02, KHZ = 1.5, KHZG = 0.7, FHZG = 0.5,
  KHZR = 110, NHZ = 4, WHZ = 9, KHZC = 0.03,
  KHEM0 = 0.0005, KHEM1 = 0.012, SELMET = 2.5, TAUEPO = 60,
  KEPO = 0.0042, GHBH = 0.55, KFHB = 0.65, KFHBEL = 0.85,
  HP0 = 12, KHP = 1.1, KHPSYN = 0.019, KBILI = 0.55,
  KBILEL = 0.14, KLDH = 250, KLDHEL = 0.02, KADAP = 1,
  FDAP = 0.9, V1DAP = 40, V2DAP = 30, QDAP = 5,
  CLDAP = 1.6, FNOH = 0.1, KBILDAP = 0.55, KEHC = 0.16,
  MWDAP = 248.3, KNOHEL = 3, KNOHIN = 9, KNOHOUT = 2.2,
  IC50CIM = 2.6, KACIM = 1.4, KELCIM = 0.29, VCIM = 90,
  KABZC = 6, KELBZC = 5.5, FTOL = 0.3, VBZC = 60,
  MWBZC = 165.2, KTOLEL = 5, KTOLIN = 8, KTOLOUT = 2,
  KANIT = 1.2, KNITEL = 1.1, KNITIN = 12, KNITOUT = 3,
  MWNIT = 69, KELASC = 0.35, KASCIN = 0.55, KASCOUT = 0.3,
  VASC = 18, KPBR = 90, KIMAO = 5.5, KSSRI = 1,
  KREL5 = 1, KUP5 = 0.62, KMAO5 = 0.38, E660O = 0.32,
  E940O = 1.21, E660D = 3.23, E940D = 0.69, E660M = 15,
  E940M = 15, E660S = 20, E940S = 4, SPOA = 110,
  SPOB = 25, CYANDEOXY = 5, CYANMET = 1.5, CYANSULF = 0.5,
  KLAC = 0.03, KLACEL = 0.45, RIBO = 0

$GLOBAL
// ---------------------------------------------------------------------------
//  Fixed physical constants.  HEMEC is the only unit conversion in the model
//  and it appears twice: once to turn a haemoglobin mass into a haem molarity,
//  and once (as CONV) to turn a red-cell flux back into a haemoglobin mass.
// ---------------------------------------------------------------------------
#define HEMEC 2482.3            // uM haem per (g/dL) haemoglobin, per L blood
#define HCTF  0.45              // haematocrit
#define CONV  (HCTF/HEMEC)      // = 1.8129e-4  g/dL per uM(RBC)
#define VBL   5.0               // L, blood volume
#define VRBCL (VBL*HCTF)        // L, red cell water
#define PVREF 38.8              // mmHg, the normal mixed-venous PO2

$CMT @annotated
  DAPG    : dapsone, gut (mg)
  DAPC    : dapsone, central (mg)
  DAPP    : dapsone, peripheral (mg)
  DAPB    : dapsone, biliary reservoir (mg)
  NOHC    : dapsone hydroxylamine, plasma (umol)
  NOHR    : dapsone hydroxylamine, red cell (uM RBC)
  NOAR    : dapsone nitrosoarene, red cell (uM RBC)
  BZCD    : benzocaine, mucosal depot (mg)
  BZCC    : benzocaine, plasma (mg)
  TOLR    : o-toluidine N-hydroxy metabolite, red cell (uM RBC)
  TOLN    : o-toluidine nitroso species, red cell (uM RBC)
  NITG    : nitrite, gut (umol)
  NITC    : nitrite, plasma (umol)
  NITR    : nitrite, red cell (uM RBC)
  MBC     : methylene blue, central (umol)
  MBP     : methylene blue, peripheral (umol)
  MBOX    : methylene blue, oxidised, in red cell (umol)
  LMB     : leucomethylene blue, in red cell (umol)
  MBEL    : methylene blue eliminated, cumulative (umol)
  ASCC    : ascorbate, plasma (umol)
  ASCR    : ascorbate, red cell (uM RBC)
  CIMG    : cimetidine, gut (mg)
  CIMC    : cimetidine, plasma (mg)
  SSRI    : SSRI exposure, normalised
  HBF     : FUNCTIONAL (ferrous) haemoglobin (g/dL)
  MHB     : METHAEMOGLOBIN (g/dL)
  SHB     : sulfhaemoglobin (g/dL)
  GSH     : reduced glutathione (uM RBC)
  GSSG    : oxidised glutathione (uM RBC)
  ROS     : reactive oxygen species (uM RBC)
  HEINZ   : Heinz-body burden (0-1)
  FHB     : plasma free haemoglobin (uM haem)
  HAPT    : haptoglobin (uM)
  BILI    : indirect bilirubin (mg/dL)
  LDH     : lactate dehydrogenase (U/L)
  EPOD    : marrow output (g/dL/h)
  BPG     : 2,3-BPG, relative to normal
  CO      : cardiac output (L/min)
  LAC     : lactate (mmol/L)
  INJ     : cumulative hypoxic dose (mmHg.h below critical)
  O2DEBT  : cumulative oxygen debt (mL/kg)
  HT5     : synaptic serotonin, relative to normal
  EB5     : CYB5R3 activity fraction
  MBCUM   : cumulative methylene blue dose (mg/kg)

$MAIN
  // The constitutive HbM floor is carried in the initial condition: those
  // haems are held ferric by the globin itself and no reductase reaches them.
  HBF_0   = HB0 * (1.0 - FMET0 - FHBM);
  MHB_0   = HB0 * (FMET0 + FHBM);
  SHB_0   = 0.0;
  GSH_0   = GSH0;
  GSSG_0  = 1.0;
  ROS_0   = KROSB / (KCAT * (1.0 + WGPX));
  HAPT_0  = HP0;
  BILI_0  = 0.6;
  LDH_0   = 180.0;
  BPG_0   = 1.0;
  CO_0    = CO0;
  LAC_0   = 1.0;
  HT5_0   = 1.0;
  EB5_0   = EB5SET;
  EPOD_0  = KHEM0 * HB0;

$ODE
  // =========================================================================
  //  0.  Guarded state readings
  // =========================================================================
  double hbf   = (HBF > 1e-9) ? HBF : 1e-9;
  double mhb   = (MHB > 0.0)  ? MHB : 0.0;
  double shb   = (SHB > 0.0)  ? SHB : 0.0;
  double gsh   = (GSH > 0.0)  ? GSH : 0.0;
  double gssg  = (GSSG > 0.0) ? GSSG : 0.0;
  double ros   = (ROS > 0.0)  ? ROS : 0.0;
  double hbtot = hbf + mhb + shb;

  // methaemoglobin ABOVE the constitutive HbM floor is the only reducible part
  double mhbred = mhb - FHBM * hbtot;
  if (mhbred < 0.0) mhbred = 0.0;

  // =========================================================================
  //  1.  OXYGEN TRANSPORT  (recomputed here because the ODEs feed back on it)
  // =========================================================================
  double fred = mhb / (hbf + mhb);                 // ferric fraction of haems
  double neff = 1.0 + (NHILL - 1.0) * (1.0 - fred);  // cooperativity is lost
  double p50e = P50 * (1.0 + BBPG * (BPG - 1.0)) * (1.0 - ALPHAM * fred);
  if (p50e < 1.0) p50e = 1.0;

  double sao2 = pow(PAO2, neff) / (pow(p50e, neff) + pow(PAO2, neff));
  double cao2 = 1.34 * hbf * sao2 + 0.003 * PAO2;
  double co   = (CO > 0.5) ? CO : 0.5;
  double cvo2 = cao2 - VO2 / (co * 10.0);
  double caph = 1.34 * hbf;
  double svo2, pvo2;
  if (cvo2 >= caph) {                    // dissolved oxygen alone suffices:
    svo2 = 1.0;                          // this is the hyperbaric-oxygen branch
    pvo2 = (cvo2 - caph) / 0.003;
    if (pvo2 < 0.0) pvo2 = 0.0;
  } else {
    svo2 = cvo2 / caph;
    if (svo2 < 1e-6)       svo2 = 1e-6;
    if (svo2 > 1.0 - 1e-6) svo2 = 1.0 - 1e-6;
    pvo2 = p50e * pow(svo2 / (1.0 - svo2), 1.0 / neff);
  }
  if (pvo2 < 0.0) pvo2 = 0.0;

  // =========================================================================
  //  2.  OXIDANT DRUG PK AND BIOACTIVATION
  // =========================================================================
  double cdap = DAPC / V1DAP;                          // mg/L
  double ccim = CIMC / VCIM;                           // mg/L
  double fnoh = FNOH / (1.0 + ccim / IC50CIM);         // cimetidine acts HERE

  double absdap = KADAP * DAPG;
  double cldap  = CLDAP * cdap;
  double distd  = QDAP * (cdap - DAPP / V2DAP);
  double bild   = KBILDAP * cdap;
  dxdt_DAPG = -absdap + KEHC * DAPB;
  dxdt_DAPC = FDAP * absdap - cldap - distd - bild;
  dxdt_DAPP = distd;
  dxdt_DAPB = bild - KEHC * DAPB;

  double gennoh = cldap * fnoh * 1000.0 / MWDAP;       // umol/h
  double uptnoh = KNOHIN * NOHC - KNOHOUT * NOHR * VRBCL;
  // The nitrosoarene is reduced BACK to the hydroxylamine by glutathione.
  // This single edge is what makes the poison catalytic; without it, a whole
  // therapeutic dose of dapsone could not oxidise 1% of the haemoglobin mass.
  double vrecyc = KRECYC * NOAR * gsh / (KGR2 + gsh);
  double vnoh   = KNOH * NOHR * hbf;                   // g/dL/h of Hb oxidised
  dxdt_NOHC = gennoh - KNOHEL * NOHC - uptnoh;
  dxdt_NOHR = uptnoh / VRBCL - vnoh / (2.0 * CONV) + vrecyc;
  dxdt_NOAR = vnoh / (2.0 * CONV) - vrecyc - KNOARD * NOAR;

  double absbzc = KABZC * BZCD;
  dxdt_BZCD = -absbzc;
  dxdt_BZCC = absbzc - KELBZC * BZCC;
  double gentol = FTOL * KELBZC * BZCC * 1000.0 / MWBZC;   // umol/h
  double vtol   = KTOL * TOLR * hbf;
  double vrcyt  = KRECYC * TOLN * gsh / (KGR2 + gsh);
  dxdt_TOLR = gentol / VRBCL - KTOLEL * TOLR - vtol / (2.0 * CONV) + vrcyt;
  dxdt_TOLN = vtol / (2.0 * CONV) - vrcyt - KTOLND * TOLN;

  double uptnit = KNITIN * NITC - KNITOUT * NITR * VRBCL;
  dxdt_NITG = -KANIT * NITG;
  dxdt_NITC = KANIT * NITG - KNITEL * NITC - uptnit;
  // Nitrite oxidation of oxyhaemoglobin is AUTOCATALYTIC (KNIT2).  This term
  // was added expecting it to produce the clinical lag-then-wall.  IT DOES
  // NOT: switching it off moves the 5 g peak by 0.5 percentage points and the
  // rise time not at all, because at poisoning doses the reaction is limited
  // by absorption on the way up and by stoichiometry at the top (nitrite,
  // unlike the arylhydroxylamines, is CONSUMED as it goes).  The term is kept
  // and the refutation is reported in mhb_reference_output.txt section 9.
  double vnit = (KNIT + KNIT2 * mhb) * NITR * hbf;
  dxdt_NITR = uptnit / VRBCL - vnit / (2.0 * CONV);

  dxdt_CIMG = -KACIM * CIMG;
  dxdt_CIMC = KACIM * CIMG - KELCIM * CIMC;

  double uptasc = KASCIN * ASCC - KASCOUT * ASCR * VRBCL;
  double vasc   = KASC * ASCR * mhbred;
  dxdt_ASCC = -KELASC * ASCC - uptasc;
  dxdt_ASCR = uptasc / VRBCL - 0.10 * ASCR - vasc / CONV;

  dxdt_SSRI = -0.029 * SSRI;

  // =========================================================================
  //  3.  METHYLENE BLUE PK
  // =========================================================================
  double cmbc = MBC / VC;                              // uM in plasma
  dxdt_MBC  = -CLMB * cmbc - QMB * (cmbc - MBP / VP) - KUP * MBC + KOUT * MBOX;
  dxdt_MBP  = QMB * (cmbc - MBP / VP);
  dxdt_MBEL = CLMB * cmbc;

  double cmbox = MBOX / VRBCL;                         // uM in red cell water
  double clmbr = LMB  / VRBCL;

  // =========================================================================
  //  4.  NADPH — ONE CAPPED SUPPLY, TWO CONSUMERS
  //      This block is the whole model.  Nothing downstream knows that G6PD
  //      deficiency is a contraindication, that 7 mg/kg is a ceiling, or that
  //      methylene blue can worsen methaemoglobinaemia.  All three are
  //      consequences of these five lines.
  // =========================================================================
  double vfrdes = VFRMAX * cmbox / (KMFR + cmbox);     // the antidote's demand
  double vgrdes = VGRMAX * gssg  / (KMGR + gssg);      // the membrane's demand
  double capn   = VPPP * G6PD * (1.0 + GRIBO * RIBO);
  double share  = capn / (capn + vfrdes + vgrdes + 1e-9);
  double vfr    = vfrdes * share;                      // uM(RBC)/h leucoMB made
  double vgr    = vgrdes * share;                      // uM(RBC)/h GSSG reduced

  // Oxidised methylene blue also takes electrons straight from haemoglobin.
  // That is not a loss of drug -- it is a LOAN, and it makes leucoMB too.
  double vmboxh = KMBOXH * cmbox * hbf;                          // g/dL/h
  double vmred  = KMBR * clmbr * mhbred / (KMHB + mhbred);       // uM(RBC)/h
  double vo2    = KO2 * clmbr;                    // the futile branch: O2 wins

  dxdt_LMB  = (vfr + vmboxh / (2.0 * CONV) - vmred - vo2 - KMBLK * clmbr) * VRBCL;
  dxdt_MBOX = KUP * MBC - KOUT * MBOX
            + (-vfr - vmboxh / (2.0 * CONV) + vmred + vo2) * VRBCL;

  // =========================================================================
  //  5.  GLUTATHIONE (CONSERVED), ROS, HEINZ BODIES
  // =========================================================================
  double vrcyctot = vrecyc + vrcyt;
  double vgshox   = KGSHOX * ros * gsh / (KG3 + gsh);
  double gt       = gsh + 2.0 * gssg;
  dxdt_GSH  = KSYN * (GSH0 - gt) + 2.0 * vgr - vgshox - 2.0 * vrcyctot;
  dxdt_GSSG = 0.5 * vgshox + vrcyctot - vgr;
  dxdt_ROS  = KROSB + vo2 + KRCROS * vrcyctot
            - KCAT * ros * (1.0 + WGPX * gsh / GSH0);

  // Heinz bodies have TWO drivers and they are not the same driver.  One is
  // frank radical flux (what high-dose methylene blue does to a normal cell).
  // The other is collapse of the glutathione pool itself (what happens in G6PD
  // deficiency, where there is no radical excess at all -- there is simply
  // nothing left holding the beta-93 thiols reduced).  With only the first
  // driver the model predicts, wrongly, that the G6PD-deficient patient
  // haemolyses LESS than the normal one.
  double rr   = pow(ros, NHZ) / (pow(KHZR, NHZ) + pow(ros, NHZ));
  double gcol = 1.0 - gsh / (FHZG * GSH0);
  if (gcol < 0.0) gcol = 0.0;
  gcol = gcol * gcol;
  dxdt_HEINZ = (KHZ * rr + KHZG * gcol) * (1.0 - HEINZ)
             / (1.0 + WHZ * gsh / GSH0) - KHZC * HEINZ;

  double vrosox = KROSOX * ros * hbf;
  double vsulf  = KSULF  * ros * hbf;      // irreversible: never comes back

  // =========================================================================
  //  6.  HAEMOGLOBIN BOOK-KEEPING
  // =========================================================================
  double khem  = KHEM0 + KHEM1 * HEINZ;
  double lysf  = khem * hbf;
  double lysm  = khem * SELMET * mhb;      // oxidised cells are cleared FIRST,
  double lyss  = khem * SELMET * shb;      // which is why %MetHb can FALL while
  double lystot = lysf + lysm + lyss;      // the patient gets worse

  double vb5   = VMXB5 * EB5 * mhbred / (KMB5 + mhbred);
  double vauto = KAUTO * hbf;

  double oxtot  = vauto + vnoh + vtol + vnit + vrosox + vmboxh;
  double redtot = vb5 + vasc + 2.0 * CONV * vmred;

  dxdt_HBF = -oxtot + redtot - vsulf - lysf + EPOD;
  dxdt_MHB =  oxtot - redtot - lysm;
  dxdt_SHB =  vsulf - lyss;

  // =========================================================================
  //  7.  CONSEQUENCES OF HAEMOLYSIS
  // =========================================================================
  dxdt_FHB  = KFHB * lystot * HEMEC / 100.0 - KFHBEL * FHB;
  dxdt_HAPT = KHPSYN * (HP0 - HAPT) - KHP * FHB * HAPT / (HAPT + 2.0);
  dxdt_BILI = KBILI * lystot - KBILEL * BILI;
  dxdt_LDH  = KLDH * lystot - KLDHEL * (LDH - 180.0);

  double hyp0 = 1.0 - pvo2 / PVREF;
  if (hyp0 < 0.0) hyp0 = 0.0;
  double prod = KEPO * (HB0 * (1.0 + GHBH * hyp0) - hbtot);
  if (prod < 0.0) prod = 0.0;
  prod = prod + KHEM0 * HB0;
  dxdt_EPOD = (prod - EPOD) / TAUEPO;

  // =========================================================================
  //  8.  WHOLE-BODY COMPENSATION
  //      2,3-BPG is the ONLY thing in this model that right-shifts the curve,
  //      i.e. the exact opposite allosteric move to the one methaemoglobin
  //      makes.  That is why chronic 25% is well and acute 25% is not.
  // =========================================================================
  double cotgt = CO0 * (1.0 + GCO * hyp0);
  if (cotgt > CO0 * COMAX) cotgt = CO0 * COMAX;
  dxdt_CO  = (cotgt - CO) / TAUCO;
  dxdt_BPG = ((1.0 + GBPG * hyp0) - BPG) / TAUBPG;

  double below = PVCRIT - pvo2;
  if (below < 0.0) below = 0.0;
  dxdt_LAC    = KLAC * below - KLACEL * (LAC - 1.0);
  dxdt_INJ    = below;
  dxdt_O2DEBT = below * VO2 / 1000.0 / WT;

  // =========================================================================
  //  9.  THE ANTIDOTE'S OTHER PHARMACOLOGY
  //      MAO-A sits in tissue, not plasma.  Plasma methylene blue is gone in
  //      minutes; the SSRI interaction lasts hours.  The effect site is
  //      therefore the peripheral compartment.
  // =========================================================================
  double cmbeff = KPBR * MBP / VP;
  double imao   = cmbeff / (KIMAO + cmbeff);
  double isert  = SSRI / (KSSRI + SSRI);
  dxdt_HT5 = KREL5 - KUP5 * (1.0 - isert) * HT5 - KMAO5 * (1.0 - imao) * HT5;

  dxdt_EB5   = (EB5SET * (1.0 + GRIBO * RIBO) - EB5) / TAUEB5;
  dxdt_MBCUM = 0.0;

$TABLE
  // =========================================================================
  //  OUTPUTS.  Everything the model reports but does not integrate.
  // =========================================================================
  double hbf_o   = (HBF > 1e-9) ? HBF : 1e-9;
  double mhb_o   = (MHB > 0.0)  ? MHB : 0.0;
  double shb_o   = (SHB > 0.0)  ? SHB : 0.0;
  double HBTOT   = hbf_o + mhb_o + shb_o;

  double FMET    = mhb_o / HBTOT;
  double MetPct  = 100.0 * FMET;
  double FSULF   = shb_o / HBTOT;
  double FRED    = mhb_o / (hbf_o + mhb_o);

  double NEFF    = 1.0 + (NHILL - 1.0) * (1.0 - FRED);
  double P50E    = P50 * (1.0 + BBPG * (BPG - 1.0)) * (1.0 - ALPHAM * FRED);
  if (P50E < 1.0) P50E = 1.0;

  double SAO2 = pow(PAO2, NEFF) / (pow(P50E, NEFF) + pow(PAO2, NEFF));
  double CAO2 = 1.34 * hbf_o * SAO2 + 0.003 * PAO2;
  double COt  = (CO > 0.5) ? CO : 0.5;
  double CVO2 = CAO2 - VO2 / (COt * 10.0);
  double caph_o = 1.34 * hbf_o;
  double SVO2, PVO2;
  if (CVO2 >= caph_o) { SVO2 = 1.0; PVO2 = (CVO2 - caph_o) / 0.003; }
  else {
    SVO2 = CVO2 / caph_o;
    if (SVO2 < 1e-6) SVO2 = 1e-6;
    if (SVO2 > 1.0 - 1e-6) SVO2 = 1.0 - 1e-6;
    PVO2 = P50E * pow(SVO2 / (1.0 - SVO2), 1.0 / NEFF);
  }
  if (PVO2 < 0.0) PVO2 = 0.0;

  double DO2  = CAO2 * COt * 10.0;
  double ERO2 = (CAO2 > 0.0) ? (CAO2 - CVO2) / CAO2 : 0.0;

  // ---- the pulse oximeter, modelled as the DEVICE rather than as the blood -
  // Note what is and is not fitted here.  E660M and E940M are the only optical
  // parameters that were tuned, and they were tuned to Barker & Tremper's
  // measured SpO2-vs-MetHb curve.  The 85% floor is NOT a parameter: it is
  // SPOA - SPOB evaluated at R = 1, and R -> 1 for any pigment with equal
  // absorbance at the two wavelengths.
  double fO2 = hbf_o * SAO2 / HBTOT;
  double fdx = hbf_o * (1.0 - SAO2) / HBTOT;
  double num = fO2 * E660O + fdx * E660D + FMET * E660M + FSULF * E660S;
  double den = fO2 * E940O + fdx * E940D + FMET * E940M + FSULF * E940S;
  if (den < 1e-9) den = 1e-9;
  double R    = num / den;
  double SPO2 = SPOA - SPOB * R;
  if (SPO2 > 100.0) SPO2 = 100.0;
  if (SPO2 < 0.0)   SPO2 = 0.0;
  double SAO2CO = 100.0 * fO2;            // what a co-oximeter reports
  double GAP    = SPO2 - SAO2CO;          // the saturation gap

  // ---- cyanosis is a SUM of pigments, each with its own threshold ---------
  double CYAN = (hbf_o * (1.0 - SAO2)) / CYANDEOXY
              + mhb_o / CYANMET + shb_o / CYANSULF;

  double CMBPL  = MBC / VC;               // plasma methylene blue (uM)
  double CMBRBC = (MBOX + LMB) / VRBCL;   // total red-cell methylene blue (uM)
  double LEUCOF = (MBOX + LMB > 1e-9) ? LMB / (MBOX + LMB) : 0.0;
  double CDAPL  = DAPC / V1DAP;           // dapsone (mg/L)

$CAPTURE
  HBTOT FMET MetPct FSULF FRED NEFF P50E SAO2 CAO2 CVO2 SVO2 PVO2 DO2 ERO2
  R SPO2 SAO2CO GAP CYAN CMBPL CMBRBC LEUCOF CDAPL

# =============================================================================
#  SCENARIOS
#  ---------
#  23 scenarios.  Each is an mrgsolve event object; the comment above each one
#  states what it is FOR, i.e. which claim it is testing.  Reference values in
#  the comments are the ones printed by `python3 mhb_reference_check.py --all`
#  and reproduced in mhb_reference_output.txt.
# =============================================================================
#
#  --- helper -----------------------------------------------------------------
#  mb_dose(mg_per_kg, time_h, WT = 70) converts a mg/kg dose of methylene blue
#  chloride into umol into the central compartment.
#
#    mb_dose <- function(mgkg, t = 0, wt = 70)
#      ev(time = t, amt = mgkg * wt * 1000 / 319.85, cmt = "MBC")
#
#  S01  HEALTHY BASELINE — the model must first be able to do nothing.
#         mod %>% mrgsim(end = 1440)
#         expect: MetHb 0.83%, PvO2 38.8 mmHg, SpO2 99.6%, Hb drift < 0.03 g/dL,
#                 cumulative hypoxic dose EXACTLY 0.
#
#  S02  BENZOCAINE SPRAY, UNTREATED — the natural history.
#         ev(amt = 250, cmt = "BZCD")
#         expect: peak MetHb ~31% at ~1.7 h, SpO2 ~87%, co-oximetry SaO2 ~66%,
#                 i.e. a saturation gap of ~21 points.
#
#  S03  BENZOCAINE + METHYLENE BLUE 1 mg/kg at 45 min — the intended use.
#         ev(amt = 250, cmt = "BZCD") + mb_dose(1, 0.75)
#
#  S04  BENZOCAINE IN AN ANAEMIC PATIENT — the same percentage, a different
#       disease.  Set HB0 = 8 and repeat S02.  The oximeter reads the SAME
#       number; PvO2 does not.
#
#  S05  DAPSONE 100 mg/d for 21 d — chronic, formation-limited disease.
#         ev(amt = 100, cmt = "DAPG", ii = 24, addl = 20)
#         expect: steady-state MetHb ~5.7%.
#
#  S06  DAPSONE 300 mg/d for 21 d — dose-dependence of the same mechanism.
#
#  S07  DAPSONE 100 mg/d + CIMETIDINE 400 mg tid — attacking the FORMATION
#       term instead of the clearance term.
#
#  S08  DAPSONE 2 g OVERDOSE, untreated — the source term lasts days.
#         ev(amt = 2000, cmt = "DAPG")
#
#  S09  DAPSONE 2 g + a SINGLE methylene blue 2 mg/kg at 6 h — produces
#       REBOUND, which the model was never told about.  Two clocks, different
#       speeds.
#
#  S10  DAPSONE 2 g + FIVE methylene blue doses (7 mg/kg cumulative, i.e. the
#       clinical ceiling) vs S11.
#
#  S11  DAPSONE 2 g + TWO methylene blue doses + CIMETIDINE — the therapeutic
#       argument, made as an argument about which TERM is rate-limiting.
#
#  S12  METHYLENE BLUE DOSE-RESPONSE SWEEP 0.5 -> 10 mg/kg on a fixed
#       benzocaine exposure — the inverted U.  Nothing in the parameter list
#       contains a ceiling; the curve has one anyway.
#
#  S13  G6PD SWEEP 1.00 -> 0.02 with methylene blue 2 mg/kg — the antidote
#       failing, and then doing harm.  Set G6PD.
#
#  S14  G6PD 0.02 GIVEN ASCORBATE INSTEAD — the alternative that is not a
#       slower version of the same thing.
#         ev(amt = 56786, cmt = "ASCC")     # 10 g ascorbic acid, umol
#
#  S15  CONGENITAL CYB5R3 DEFICIENCY — EB5SET = 0.035, run 250 d.  Blue and
#       well: 2,3-BPG and compensatory erythrocytosis restore tissue PO2 while
#       the colour and the oximeter stay abnormal.
#
#  S16  CONGENITAL CYB5R3 + METHYLENE BLUE — a COSMETIC intervention: the
#       colour changes, the tissue PO2 barely moves.
#
#  S17  CONGENITAL CYB5R3 + RIBOFLAVIN — RIBO = 1.
#
#  S18  HbM VARIANT — FHBM = 0.25 + methylene blue 2 mg/kg.  No response, and
#       the model contains no rule saying so.
#
#  S19  SULFHAEMOGLOBINAEMIA — initialise SHB and compare the oximeter, the
#       cyanosis threshold and the (absent) methylene blue response.
#
#  S20  SODIUM NITRITE 2 / 5 / 10 g ingestion — a CONSUMED oxidant, in
#       contrast to the catalytic ones, so the peak tracks the dose.
#         ev(amt = 5000 * 1000 / 69, cmt = "NITG")
#
#  S21  HYPERBARIC OXYGEN at MetHb 70-95% — set PAO2 = 1800.  Dissolved oxygen
#       covers resting extraction with essentially no functional haemoglobin
#       (Boerema's arithmetic: 0.003 x 1800 = 5.4 mL/dL vs VO2/(CO x 10) = 5.0).
#
#  S22  EXCHANGE TRANSFUSION — replaces the carrier rather than repairing it.
#
#  S23  METHYLENE BLUE ON TOP OF AN SSRI — the antidote's other pharmacology.
#         ev(amt = 2.5, cmt = "SSRI") + mb_dose(2, 0.75)
#
# =============================================================================
#  KNOWN LIMITATIONS, STATED RATHER THAN BURIED
#  --------------------------------------------
#  1. ALPHAM (the size of the Darling & Roughton left shift) is the parameter
#     the headline result rests on, and it has the weakest quantitative
#     literature of anything in the file.  The sensitivity analysis in
#     mhb_reference_check.py shows: at half the assumed shift the "MetHb 30%"
#     equivalence moves from 8.3 to ~9.4 g/dL (the argument weakens but does
#     not reverse); at zero shift it disappears entirely and %MetHb becomes an
#     adequate bedside variable.  That is a falsifiable prediction and the
#     experiment is a co-oximeter plus a measured P50 on the same sample.
#
#  2. E660M and E940M are set EQUAL by construction.  The underlying spectra
#     support "methaemoglobin absorbs strongly at both wavelengths", but the
#     exact 940 nm value is the weakest optical number here.  What the model
#     genuinely predicts is the SHAPE (a floor, approached asymptotically, with
#     collapsing sensitivity), not the precise floor height.
#
#  3. EB5 is an IN VIVO FLUX fraction, not the assay activity fraction.
#     Congenital type I patients are usually reported as having 10-20% residual
#     activity on a ferricyanide assay, but reproducing their observed 15-30%
#     methaemoglobin requires an in vivo flux nearer 3.5%.  The model is
#     calibrated to the phenotype and this discrepancy is real, not hidden.
#
#  4. The haemolysis arm is a two-driver population approximation, not a
#     cell-age-structured model.  It reproduces the direction and rough
#     magnitude of oxidant haemolysis but should not be used to predict
#     transfusion thresholds.
#
#  5. R was unavailable in the build container, so this file is
#     EQUATION-VERIFIED against the Python implementation, not compile-verified.
# =============================================================================
