###############################################################################
##  Alcohol Use Disorder (AUD) — Quantitative Systems Pharmacology model
##  ---------------------------------------------------------------------------
##  aud_mrgsolve_model.R          companion to aud_qsp_model.dot
##
##  WHY THIS MODEL IS BUILT AS A CLOSED LOOP
##  ----------------------------------------
##  In every other model in this library the dose is an input. In AUD the dose
##  is an OUTPUT: the patient chooses it, and the pharmacology of what they
##  already drank is what does the choosing. So `drinks per day` is not a
##  covariate here — it is a model variable, generated each simulated evening
##  from the state of the reward, stress and control systems:
##
##      DRIVE ──► drinks/day ──► ethanol PK ──► MOR/dopamine reinforcement
##        ▲                                            │
##        │                                            ▼
##        └── craving ◄── negative affect ◄── allostatic neuroadaptation
##                                                     │
##                          withdrawal (unmasked when BAC falls) ──┘
##
##  Consequences of building it this way:
##    * the disease is GENERATED, not assumed. A healthy drinker and a severe
##      AUD patient differ only in environment (ENVDRIVE) and vulnerability
##      (VULN); the escalation trajectory is emergent.
##    * a medication does not "reduce the endpoint" — it perturbs one arm of
##      the loop and the loop finds a new operating point. That is why
##      naltrexone and acamprosate move DIFFERENT endpoints (see §Validation).
##    * withdrawal is not a separate model. It is what happens to the same
##      equations when ethanol input stops while the adaptation states are
##      still elevated.
##
##  TWO PHARMACOLOGICALLY SEPARATE ARMS ARE MODELLED EXPLICITLY
##  ------------------------------------------------------------
##    (A) POSITIVE-REINFORCEMENT ARM  ethanol → β-endorphin → MOR → disinhibits
##        VTA GABA → phasic dopamine → within-episode escalation (priming).
##        Cut by naltrexone / 6-β-naltrexol / nalmefene / semaglutide /
##        ondansetron / baclofen / topiramate(AMPA).
##        Clinical signature: DRINKS PER DRINKING DAY, %HDD.
##    (B) NEGATIVE-REINFORCEMENT ARM  chronic ethanol → GABA-A subunit
##        reconfiguration + NR2B up + GLT-1 down + dynorphin/CRF up →
##        negative affect, insomnia, protracted withdrawal → relief craving.
##        Cut by acamprosate / gabapentin / topiramate(GABA) / benzodiazepines
##        (detox only).
##        Clinical signature: PERCENT DAYS ABSTINENT, time to first drink.
##
##  A third, upstream arm is aversive (disulfiram → ALDH2 inactivation →
##  acetaldehyde → flush), and it is gated by supervision, not by exposure.
##
##  72 ODE compartments · 268 parameters · 20 scenarios.
##
##  Run:   Rscript aud_mrgsolve_model.R
###############################################################################

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
})

options(stringsAsFactors = FALSE)

###############################################################################
## MODEL CODE
###############################################################################

code <- '
$PROB
# Alcohol Use Disorder QSP model (closed-loop drinking)

$PARAM @annotated
// ---------------- patient / demography -----------------------------------
BW      :  75   : body weight (kg)
SEXF    :   0   : female flag (0/1) - sets Widmark r and gastric ADH
RWIDM   : 0.68  : Widmark distribution ratio, male (L/kg)
RWIDF   : 0.55  : Widmark distribution ratio, female (L/kg)
FED     : 0.5   : fed state 0-1 (slows gastric emptying)
GDRINK  :  14   : grams ethanol per US standard drink (WHO unit = 10 g)

// ---------------- pharmacogenetics ---------------------------------------
ADHF    : 1.0   : ADH1B activity multiplier (ADH1B*2 ~ 2.5)
ALDHF   : 1.0   : ALDH2 low-Km activity multiplier (*2 het ~ 0.2, *2*2 ~ 0.03)
OPRMF   : 1.0   : OPRM1 118G MOR-blockade potency multiplier (G carrier ~ 1.6)
HTTF    : 1.0   : 5-HTTLPR LL ondansetron responsiveness multiplier
VULN    : 1.0   : composite neuro-vulnerability (family history, low LR)

// ---------------- drinking behaviour --------------------------------------
ENVDRIVE : 2.5  : environmental / social drinking pressure (drive units)
DPDMAX   : 11   : maximum sustainable drinks per day
DR50     : 4.3  : drive producing half-maximal drinks per day
TDMID    : 20.0 : centre of the daily drinking window (h of day)
TDSIG    : 0.9  : SD of the daily drinking window (h)
FORCEDPD : -1   : if >=0 override endogenous drinks/day (history/burn-in)
DRINKSW  : 1    : master drinking switch (0 = enforced abstinence/inpatient)
ABSTG    : 0    : abstinence goal strength 0-1 (treatment contract)
WGOAL    : 0.55 : maximum fractional drive suppression by abstinence goal
WPFC     : 0.9  : weight of prefrontal control on drive
WPRIME   : 0.55 : opioid priming gain on within-episode consumption
ESEMINT  : 0.16 : semaglutide brake on within-episode consumption
ETPMINT  : 0.28 : topiramate brake on within-episode consumption
EBACINT  : 0.12 : baclofen brake on within-episode consumption
EONDINT  : 0.12 : ondansetron brake on within-episode consumption
MOR0     : 0.25 : baseline (drug-free, ethanol-free) MOR activation
KDRIVE   : 0.35 : drive equilibration rate (1/h)

// ---------------- ethanol ADME --------------------------------------------
KEMPT0  : 3.2   : gastric emptying rate constant (1/h), fasted
FEDK    : 0.65  : fractional slowing of emptying when fed
KAETH   : 5.0   : intestinal absorption rate constant (1/h)
KGAST   : 0.28  : gastric first-pass ADH rate constant (1/h)
GASTFM  : 1.0   : gastric ADH factor male
GASTFF  : 0.45  : gastric ADH factor female
KRENETH : 0.006 : renal + pulmonary ethanol loss (1/h)
VMADH0  : 6.40  : hepatic ADH Vmax (g/h at 70 kg)
KMADH   : 0.05  : ADH Km (g/L)
VMCYP0  : 1.20  : CYP2E1/MEOS Vmax (g/h) at CYP2E1 = 1
KMCYP   : 0.40  : MEOS Km (g/L) ~ 8-10 mM
VMCAT   : 0.30  : catalase Vmax (g/h)
KMCAT   : 0.50  : catalase Km (g/L)

// ---------------- CYP2E1 induction ----------------------------------------
KOUTCYP : 0.0072 : CYP2E1 turnover (1/h), t1/2 ~ 4 d
EMAXCYP : 1.60   : maximal fold induction above baseline
EC50CYP : 0.45   : ethanol conc for half-maximal induction (g/L)

// ---------------- acetaldehyde --------------------------------------------
VMALD2  : 5000  : ALDH2 (low Km) Vmax (uM/h)
KMALD2  : 1.0   : ALDH2 Km (uM)
VMALD1  : 26000 : ALDH1A1 (high Km) Vmax (uM/h)
KMALD1  : 150   : ALDH1A1 Km (uM)
DSFALD1 : 0.25  : fraction of ALDH1A1 inhibited by disulfiram

// ---------------- redox / lactate -----------------------------------------
KNADH   : 0.9   : NADH equilibration rate (1/h)
ENADH   : 1.8   : maximal NADH/NAD ratio rise
EC50NAD : 3.0   : ethanol oxidation flux for half-max redox shift (g/h)
LACT0   : 1.0   : baseline lactate (mmol/L)
KLACT   : 0.35  : lactate turnover (1/h)
ELACT   : 1.9   : maximal lactate rise from redox shift

// ---------------- chronic-exposure integrator ------------------------------
TAUAVG  : 168   : time constant of the chronic ethanol exposure average (h)
EC50CHR : 0.45  : chronic ethanol average giving half-maximal neuroadaptation (g/L)

// ---------------- GABA-A --------------------------------------------------
KGT     : 4.0   : GABA tone equilibration (1/h)
EMAXGA  : 1.10  : maximal GABA-A potentiation by ethanol
EC50GA  : 25.0  : ethanol concentration for half-max potentiation (mM)
WSUBTOL : 1.30  : subunit-reconfiguration blunting of acute potentiation
WGTSUB  : 0.65  : loss of resting GABA tone per unit subunit index
KSUBON  : 0.0360 : GABA-A reconfiguration on-rate (per CHRON unit per h)
KSUBOFF : 0.0048 : reconfiguration reversal (1/h), t1/2 ~ 6 d
SUBMAX  : 1.0   : maximum subunit reconfiguration index

// ---------------- glutamate / NMDA ----------------------------------------
KNMON   : 0.0220 : NR2B upregulation on-rate (per CHRON unit per h)
KNMOFF  : 0.0041 : NR2B reversal (1/h), t1/2 ~ 7 d
NMMAX   : 1.0    : maximum NMDA upregulation index
KGLTON  : 0.0047 : GLT-1/xCT downregulation rate (per CHRON unit per h)
KGLTOFF : 0.0029 : GLT-1 recovery rate (1/h), t1/2 ~ 10 d
KGLU    : 1.2    : extracellular glutamate equilibration (1/h)
WGLUT   : 1.35   : glutamate rise per unit loss of GLT-1 capacity
WGLUWD   : 0.80  : withdrawal hyperglutamatergia weight
EACP    : 0.70   : maximal acamprosate normalisation of glutamate homeostasis
EC50ACP : 260    : acamprosate plasma conc for half-max effect (ng/mL)

// ---------------- opioid / dopamine ---------------------------------------
KBEND   : 2.5   : beta-endorphin equilibration (1/h)
EBEND   : 2.2   : maximal ethanol-evoked beta-endorphin rise
EC50BE  : 16.0  : ethanol conc for half-max endorphin release (mM)
KDBEND  : 3.0   : beta-endorphin KD at MOR (relative units)
KINTX   : 0.20  : naltrexone Ki at MOR (nM)
KINOL   : 5.0   : 6-beta-naltrexol Ki at MOR (nM)
KINMF   : 0.08  : nalmefene Ki at MOR (nM)
KDA     : 2.0   : accumbal dopamine equilibration (1/h)
EDAMOR  : 2.6   : dopamine gain from MOR activation
EDADIR  : 0.85  : direct (MOR-independent) ethanol dopamine component
EC50DA  : 18.0  : ethanol conc for half-max direct DA effect (mM)
EDAKOR  : 0.22  : dopamine suppression per unit dynorphin tone above 1
FTONIC  : 0.25  : fraction of the MOR gain that applies to TONIC (sub-baseline) blockade

// ---------------- incentive sensitisation / habit --------------------------
KCUEON  : 0.0055 : cue-reactivity acquisition rate (per squared reinforcement unit/h)
REINFTH : 0.35  : phasic-dopamine threshold below which no sensitisation occurs
KCUEOFF : 0.00035 : cue-reactivity extinction (1/h), t1/2 ~ 82 d
CUEMAX  : 12.0  : maximum incentive-salience index
KHABON  : 0.00130 : habit acquisition rate (1/h)
PFCH    : 0.85  : executive-control level below which habit capture begins
KHABOFF : 0.00016 : habit decay (1/h), t1/2 ~ 180 d
KALLON  : 0.00380 : allostatic set-point shift rate (per CHRON unit per h)
KALLOFF : 0.00096 : allostatic recovery (1/h), t1/2 ~ 30 d
ALLOMAX : 6.0   : maximum reward set-point shift

// ---------------- dynorphin / CRF / NPY ------------------------------------
KDYNON  : 0.0074 : dynorphin induction (per CHRON unit per h)
KDYNOFF : 0.0021 : dynorphin reversal (1/h), t1/2 ~ 14 d
DYNMAX  : 3.2   : maximum dynorphin tone
KCRFON  : 0.0062 : CeA CRF recruitment (per CHRON unit per h)
KCRFOFF : 0.0021 : CRF reversal (1/h)
CRFMAX  : 3.6   : maximum CRF tone
KNPYON  : 0.0039 : NPY depletion rate (per CHRON unit per h)
KNPYOFF : 0.0024 : NPY recovery (1/h)
KNE     : 0.8   : locus coeruleus NE equilibration (1/h)
WNECRF  : 0.55  : NE drive from CRF
WNEEXC  : 0.35  : NE drive from CNS excitation

// ---------------- HPA axis --------------------------------------------------
KCRH    : 0.9   : hypothalamic CRH turnover (1/h)
KACTH   : 1.4   : ACTH turnover (1/h)
KCORT   : 0.35  : cortisol turnover (1/h)
WCRHAM  : 0.45  : amygdalar CRF drive on PVN CRH
WCRHETH : 0.55  : acute ethanol drive on PVN CRH
EC50CRH : 20.0  : ethanol conc for half-max acute HPA activation (mM)
KGRON   : 0.0032 : GR desensitisation rate (per CHRON unit per h)
KGROFF  : 0.0016 : GR resensitisation (1/h)
WFEED   : 1.1   : strength of cortisol negative feedback

// ---------------- prefrontal control ----------------------------------------
KPFCDAM : 0.00077 : PFC erosion rate (per CHRON unit per h)
KPFCREC : 0.00048 : PFC recovery rate (1/h), t1/2 ~ 60 d
WGLUPFC : 0.22  : additional PFC impairment from extracellular glutamate

// ---------------- negative affect / sleep -----------------------------------
KNEG    : 0.25  : negative-affect equilibration (1/h)
WNCRF   : 1.35  : negative affect weight on CRF
WNDYN   : 0.95  : negative affect weight on dynorphin
WNGLU   : 0.85  : negative affect weight on glutamate
WNALLO  : 0.22  : negative affect weight on allostatic shift
WNSLP   : 0.35  : negative affect weight on sleep disruption
WNNPY   : 0.60  : NPY buffering of negative affect
KSLP    : 0.30  : sleep-disruption equilibration (1/h)
WSLPEXC : 0.85  : sleep disruption from CNS excitation
WSLPNEG : 0.45  : sleep disruption from negative affect
EGBP    : 0.70  : maximal gabapentin effect on negative affect / sleep
EC50GBP : 3.6   : gabapentin plasma conc for half-max effect (mg/L)

// ---------------- craving ----------------------------------------------------
KCRAVE  : 0.30  : craving equilibration (1/h)
WINC    : 0.42  : incentive-salience weight in craving
WREL    : 0.70  : relief (negative-affect) weight in craving
WHAB    : 1.10  : habit weight in craving
CUEEXP  : 1.0   : environmental cue-exposure multiplier

// ---------------- withdrawal / kindling ---------------------------------------
KCIWA   : 0.060 : CIWA equilibration (1/h) - sets the 24-48 h peak
CIWAMAX : 30.0  : maximum attainable CIWA-Ar score
EC50EXC : 3.40  : excitation index producing half-maximal CIWA
HILLEXC : 2.6   : Hill coefficient of the CIWA response
WNMEXC  : 1.30  : NMDA upregulation weight in the excitation index
WGAEXC  : 1.45  : GABA deficit weight in the excitation index
WNEEXC2 : 0.40  : NE weight in the excitation index
WCRFEXC : 0.45  : CRF weight in the excitation index
WKINDLE : 0.55  : kindling amplification of the excitation index
KKINON  : 0.000012 : kindling accrual per CIWA unit above threshold per h (across episodes)
CIWAKTH : 15.0  : CIWA threshold above which kindling accrues
KKINOFF : 0.000012 : kindling decay (1/h), t1/2 ~ 6.6 y
KINDMAX : 2.5   : maximum kindling index

// ---------------- aversive expectancy (disulfiram) ---------------------------
KAVON   : 0.35  : aversive-expectancy acquisition (1/h)
KAVOFF  : 0.012 : aversive-expectancy decay (1/h), t1/2 ~ 2.4 d
AVMAX   : 3.6   : maximal aversion from an acetaldehyde flush
AVAC50  : 12.0  : acetaldehyde conc for half-max flush aversion (uM)
AVHILL  : 2.0   : Hill coefficient of the flush response
WKNOW   : 2.4   : anticipatory aversion from knowing ALDH2 is blocked
SUPERV  : 0     : supervised administration flag (0/1)

// ---------------- objective biomarkers ----------------------------------------
KPETHF  : 3.00  : PEth formation (ng/mL per (g/L) per h)
KPETHE  : 0.00481 : PEth elimination (1/h), t1/2 ~ 6 d
CDT0    : 1.05  : baseline CDT (%)
KCDT    : 0.00206 : CDT turnover (1/h), t1/2 ~ 14 d
ECDT    : 2.10  : maximal CDT rise (%)
EC50CDT : 0.20  : chronic ethanol average above threshold for half-max CDT rise (g/L)
CDTHR   : 0.15  : chronic ethanol average threshold for CDT (~ 60 g/day)
GGT0    : 25.0  : baseline GGT (U/L)
KGGT    : 0.00144 : GGT turnover (1/h), t1/2 ~ 20 d
EGGT    : 120   : maximal GGT rise (U/L)
EC50GGT : 0.40  : chronic ethanol average for half-max GGT rise (g/L)
MCV0    : 90.0  : baseline MCV (fL)
KMCV    : 0.00072 : MCV turnover (1/h), t1/2 ~ 40 d
EMCV    : 13.0  : maximal MCV rise (fL)
EC50MCV : 0.30  : chronic ethanol average for half-max MCV rise (g/L)

// ---------------- liver ---------------------------------------------------------
STEAT0  : 0.03  : baseline hepatic fat fraction
KSTEAT  : 0.0035 : steatosis turnover (1/h), t1/2 ~ 8 d
ESTEAT  : 0.42  : maximal steatosis fraction
EC50ST  : 0.32  : chronic ethanol average for half-max steatosis (g/L)
LPS0    : 1.0   : baseline portal endotoxin (relative)
KLPS    : 0.09  : endotoxin turnover (1/h)
ELPS    : 3.4   : maximal endotoxin rise
EC50LPS : 0.60  : chronic ethanol average for half-max leak (g/L)
TNF0    : 1.0   : baseline TNF-alpha (relative)
KTNF    : 0.06  : TNF turnover (1/h)
ETNF    : 3.0   : maximal TNF rise
ALT0    : 22.0  : baseline ALT (U/L)
KALT    : 0.0072 : ALT turnover (1/h), t1/2 ~ 4 d
EALT    : 70.0  : maximal ALT rise (U/L)
AST0    : 20.0  : baseline AST (U/L)
KAST    : 0.0139 : AST turnover (1/h), t1/2 ~ 2 d
EAST    : 185   : maximal AST rise (U/L) - larger than ALT (AST/ALT > 2)
KFIBON  : 0.0000090 : fibrosis accrual per TNF unit per h
KFIBOFF : 0.0000042 : fibrosis regression (1/h)
FIBMAX  : 4.0   : maximum fibrosis stage

// ---------------- nutrition / cardiovascular -------------------------------------
THIA0   : 1.0   : baseline thiamine status (relative)
KTHIA   : 0.0035 : thiamine turnover (1/h)
ETHIA   : 0.80  : maximal fractional thiamine depletion
EC50TH  : 0.26  : chronic ethanol average for half-max depletion (g/L)
THIADOSE: 0     : exogenous thiamine input (relative units per h)
SBP0    : 125   : baseline systolic BP (mmHg)
KSBP    : 0.0062 : SBP turnover (1/h), t1/2 ~ 4.7 d
ESBP    : 16.0  : maximal SBP rise (mmHg)
EC50SBP : 0.55  : chronic ethanol average for half-max BP rise (g/L)
KCMYO   : 0.0000060 : cardiomyopathy accrual (per g/L per h)
KCMYOR  : 0.0000030 : cardiomyopathy recovery (1/h)

// ---------------- naltrexone PK/PD --------------------------------------------
KANTX   : 1.5   : naltrexone absorption (1/h)
FNTX    : 0.10  : naltrexone systemic bioavailability (high first pass)
FMNOL   : 0.55  : PRESYSTEMIC fraction of an oral dose converted to 6-beta-naltrexol
FSYSNOL : 0.12  : fraction of SYSTEMIC naltrexone clearance yielding 6-beta-naltrexol
VNTX    : 1100  : naltrexone V (L)
CLNTX   : 190   : naltrexone CL (L/h), t1/2 ~ 4 h
VNOL    : 300   : 6-beta-naltrexol V (L)
CLNOL   : 16.0  : 6-beta-naltrexol CL (L/h), t1/2 ~ 13 h
KRELDEP : 0.0045 : XR-naltrexone depot release (1/h)
FDEPBRST: 0.05  : XR-naltrexone initial burst fraction

// ---------------- acamprosate / disulfiram --------------------------------------
KAACP   : 0.25  : acamprosate absorption (1/h)
FACP    : 0.11  : acamprosate bioavailability
VACP    : 1100  : acamprosate apparent V (L)
CLACP   : 30.5  : acamprosate CL (L/h), t1/2 ~ 25 h
KADSF   : 0.9   : disulfiram absorption (1/h)
FDSF    : 0.85  : disulfiram bioavailability
VDSF    : 200   : disulfiram (+ active metabolite) V (L)
CLDSF   : 16.0  : disulfiram CL (L/h), t1/2 ~ 8.7 h
KINACT  : 0.055 : ALDH2 inactivation rate (per mg/L per h)
KSYNALD : 0.0058 : ALDH2 resynthesis (1/h), t1/2 ~ 5 d

// ---------------- topiramate / gabapentin / baclofen -----------------------------
KATPM   : 0.7   : topiramate absorption (1/h)
FTPM    : 0.80  : topiramate bioavailability
VTPM    : 49.0  : topiramate V (L)
CLTPM   : 1.62  : topiramate CL (L/h), t1/2 ~ 21 h
ETPMDA  : 0.45  : topiramate blunting of ethanol-evoked dopamine (AMPA/GluK1)
ETPMGA  : 0.38  : topiramate GABA-A potentiation
EC50TPM : 5.5   : topiramate conc for half-max effect (mg/L)
KAGBP   : 1.1   : gabapentin absorption (1/h)
FGBP    : 0.35  : gabapentin bioavailability at 600 mg (saturable LAT1)
VGBP    : 58.0  : gabapentin V (L)
CLGBP   : 6.7   : gabapentin CL (L/h)
KABAC   : 1.3   : baclofen absorption (1/h)
FBAC    : 0.70  : baclofen bioavailability
VBAC    : 59.0  : baclofen V (L)
CLBAC   : 10.2  : baclofen CL (L/h)
EBAC    : 0.40  : maximal baclofen blunting of ethanol-evoked dopamine
EC50BAC : 300   : baclofen conc for half-max effect (ng/mL)

// ---------------- nalmefene / ondansetron / semaglutide / benzodiazepine ----------
KANMF   : 1.1   : nalmefene absorption (1/h)
FNMF    : 0.41  : nalmefene bioavailability
VNMF    : 602   : nalmefene V (L)
CLNMF   : 33.4  : nalmefene CL (L/h), t1/2 ~ 12.5 h
EKORNMF : 0.45  : nalmefene KOR partial-agonist capping of dynorphin tone
EC50KOR : 8.0   : nalmefene conc for half-max KOR effect (nM)
KAOND   : 1.4   : ondansetron absorption (1/h)
FOND    : 0.60  : ondansetron bioavailability
VOND    : 160   : ondansetron V (L)
CLOND   : 27.7  : ondansetron CL (L/h)
EOND    : 0.30  : maximal 5-HT3 blockade of ethanol-evoked dopamine
EC50OND : 0.45  : ondansetron conc for half-max effect (ug/L)
KASEM   : 0.020 : semaglutide SC absorption (1/h)
FSEM    : 0.89  : semaglutide bioavailability
VSEM    : 12.5  : semaglutide V (L)
CLSEM   : 0.0525 : semaglutide CL (L/h), t1/2 ~ 7 d
ESEMDA  : 0.30  : semaglutide blunting of ethanol-evoked dopamine
ESEMCUE : 0.28  : semaglutide blunting of incentive salience
EC50SEM : 12.0  : semaglutide conc for half-max effect (nM)
VBZD    : 80.0  : benzodiazepine (diazepam-equivalent) V (L)
CLBZD   : 1.85  : benzodiazepine CL (L/h), t1/2 ~ 30 h
BZDON   : 0     : symptom-triggered benzodiazepine protocol flag (0/1)
BZDGAIN : 0.35  : mg diazepam-equivalent per CIWA unit above threshold per h
BZDTHR  : 8.0   : CIWA-Ar threshold for symptom-triggered dosing
BZDDUR  : 168   : duration of the benzodiazepine protocol (h)
BZDMAXR : 4.0   : maximum benzodiazepine infusion rate (mg/h)
EBZD    : 1.35  : maximal benzodiazepine GABA-A substitution
EC50BZD : 0.30  : benzodiazepine conc for half-max GABA effect (mg/L)
PRAZON  : 0     : prazosin flag (0/1)
EPRAZ   : 0.45  : fractional alpha1 blockade with prazosin

$INIT @annotated
ETH_ST : 0 : stomach ethanol (g)
ETH_GUT : 0 : intestinal ethanol (g)
ETH_C : 0 : central ethanol amount (g)
ACD : 0 : blood acetaldehyde (uM)
ACT : 0 : acetate (mmol/L)
CYP2E1 : 1 : CYP2E1 relative content (fold of baseline)
NADH : 1 : hepatic NADH/NAD+ ratio (relative)
LACT : 1.0 : lactate (mmol/L)
ETHAVG : 0 : chronic ethanol exposure average (g/L)
ETOHCUM : 0 : cumulative ethanol consumed (g)
PETH : 0 : phosphatidylethanol 16:0/18:1 (ng/mL)
CDT : 1.05 : carbohydrate-deficient transferrin (%)
GGT : 25 : gamma-glutamyl transferase (U/L)
MCV : 90 : mean corpuscular volume (fL)
GABA_TONE : 1 : net GABAergic inhibitory tone (relative, 1 = normal)
GLU_NAC : 1 : extracellular glutamate NAc/mPFC (relative)
DA_NAC : 1 : accumbal dopamine (relative)
BEND : 1 : beta-endorphin tone (relative)
NE_LC : 1 : locus coeruleus noradrenergic output (relative)
CRAVE : 0 : composite craving (PACS-like 0-30)
DRIVE : 0 : net drinking drive (drive units)
GABAA_SUB : 0 : GABA-A subunit reconfiguration index (0-1)
NMDA_UP : 0 : NR2B/NMDA upregulation index (0-1)
GLT1 : 1 : astrocytic GLT-1/xCT capacity (relative, 1 = normal)
DYN : 1 : dynorphin/KOR tone (relative)
CRF_CEA : 1 : extended-amygdala CRF tone (relative)
NPY_CEA : 1 : extended-amygdala NPY (relative)
ALLO : 0 : allostatic reward set-point shift
CUE : 0 : cue-conditioned incentive salience
HABIT : 0 : dorsolateral striatal habit strength (0-1)
PFC : 1 : prefrontal executive control capacity (relative)
CRH_HYP : 1 : hypothalamic CRH (relative)
ACTH : 1 : ACTH (relative)
CORT : 1 : plasma cortisol (relative)
GR_SENS : 1 : glucocorticoid receptor sensitivity (relative)
NEGAFF : 0 : negative affect index
SLEEPD : 0 : sleep disruption index
CIWA : 0 : CIWA-Ar withdrawal score
KINDLE : 0 : withdrawal kindling index
AVEXP : 0 : aversive expectancy (disulfiram)
STEAT : 0.03 : hepatic steatosis fraction
LPS : 1 : portal endotoxin (relative)
TNFA : 1 : hepatic TNF-alpha (relative)
ALT : 22 : alanine aminotransferase (U/L)
AST : 20 : aspartate aminotransferase (U/L)
FIB : 0 : liver fibrosis stage (0-4)
THIA : 1 : thiamine status (relative)
SBP : 125 : systolic blood pressure (mmHg)
CMYO : 0 : alcoholic cardiomyopathy index (LVEF decrement units)
NTX_GUT : 0 : naltrexone gut (mg)
NTX_C : 0 : naltrexone central (mg)
NTXOL : 0 : 6-beta-naltrexol central (mg)
NTX_DEP : 0 : XR-naltrexone depot (mg)
ACP_GUT : 0 : acamprosate gut (mg)
ACP_C : 0 : acamprosate central (mg)
DSF_GUT : 0 : disulfiram gut (mg)
DSF_C : 0 : disulfiram + active metabolites central (mg)
ALDH2A : 1 : fraction of ALDH2 catalytically active (0-1)
TPM_GUT : 0 : topiramate gut (mg)
TPM_C : 0 : topiramate central (mg)
GBP_GUT : 0 : gabapentin gut (mg)
GBP_C : 0 : gabapentin central (mg)
BACL_GUT : 0 : baclofen gut (mg)
BACL_C : 0 : baclofen central (mg)
NMF_GUT : 0 : nalmefene gut (mg)
NMF_C : 0 : nalmefene central (mg)
OND_GUT : 0 : ondansetron gut (mg)
OND_C : 0 : ondansetron central (mg)
SEM_SC : 0 : semaglutide SC depot (mg)
SEM_C : 0 : semaglutide central (mg)
BZD_C : 0 : benzodiazepine central (mg diazepam-equivalent)
BZDCUM : 0 : cumulative benzodiazepine dose (mg diazepam-equivalent)

$GLOBAL
#define POSPART(x) ((x) > 0 ? (x) : 0.0)
#define MWETH 46.07
#define MWNTX 341.4
#define MWNOL 343.4
#define MWNMF 339.4
#define MWSEM 4113.6

$MAIN
double RWID  = SEXF > 0.5 ? RWIDF : RWIDM;
double VD    = RWID * BW;
double GASTF = SEXF > 0.5 ? GASTFF : GASTFM;
double VMADH = VMADH0 * pow(BW/70.0, 0.75);
double VMCYP = VMCYP0 * pow(BW/70.0, 0.75);
double KEMPT = KEMPT0 * (1.0 - FEDK*FED);

$ODE
// =====================================================================
// 0. concentrations
// =====================================================================
double AETH  = POSPART(ETH_C);
double CETH  = AETH / VD;                 // g/L
double BAC   = CETH / 10.0;               // g/dL
double ETHMM = CETH / 0.04607;            // mM (brain == plasma)

double CNTX  = POSPART(NTX_C)  / VNTX * 1e6 / 1000.0;    // ng/mL
double CNOL  = POSPART(NTXOL)  / VNOL * 1e6 / 1000.0;    // ng/mL
double CACP  = POSPART(ACP_C)  / VACP * 1e6 / 1000.0;    // ng/mL
double CDSF  = POSPART(DSF_C)  / VDSF;                   // mg/L
double CTPM  = POSPART(TPM_C)  / VTPM;                   // mg/L
double CGBP  = POSPART(GBP_C)  / VGBP;                   // mg/L
double CBAC  = POSPART(BACL_C) / VBAC * 1e6 / 1000.0;    // ng/mL
double CNMF  = POSPART(NMF_C)  / VNMF * 1e6 / 1000.0;    // ng/mL
double COND  = POSPART(OND_C)  / VOND * 1e6 / 1000.0;    // ug/L
double CSEM  = POSPART(SEM_C)  / VSEM * 1e6 / 1000.0;    // ug/L
double CBZD  = POSPART(BZD_C)  / VBZD;                   // mg/L

double NTXnM = CNTX / MWNTX * 1000.0;     // nM
double NOLnM = CNOL / MWNOL * 1000.0;
double NMFnM = CNMF / MWNMF * 1000.0;
double SEMnM = CSEM / MWSEM * 1000.0;

// =====================================================================
// 1. drug pharmacodynamic effect fractions
// =====================================================================
double ACPEFF  = EACP   * CACP  / (EC50ACP + CACP);
double TPMDA   = ETPMDA * CTPM  / (EC50TPM + CTPM);
double TPMGA   = ETPMGA * CTPM  / (EC50TPM + CTPM);
double GBPEFF  = EGBP   * CGBP  / (EC50GBP + CGBP);
double BACEFF  = EBAC   * CBAC  / (EC50BAC + CBAC);
double ONDEFF  = EOND * HTTF * COND / (EC50OND + COND);
double SEMDA   = ESEMDA  * SEMnM / (EC50SEM + SEMnM);
double SEMCUE  = ESEMCUE * SEMnM / (EC50SEM + SEMnM);
double BZDEFF  = EBZD   * CBZD  / (EC50BZD + CBZD);
double KOREFF  = EKORNMF* NMFnM / (EC50KOR + NMFnM);
double A1BLK   = PRAZON > 0.5 ? EPRAZ : 0.0;

// competitive MOR antagonism (naltrexone + 6-beta-naltrexol + nalmefene)
double ANTAG = OPRMF * (NTXnM/KINTX + NOLnM/KINOL + NMFnM/KINMF);
double MORACT = (BEND/KDBEND) / (1.0 + BEND/KDBEND + ANTAG);

// =====================================================================
// 2. drinking behaviour — the dose is generated here
// =====================================================================
double DRIVE_T = (CRAVE + ENVDRIVE)
                 * (1.0 / (1.0 + WPFC*PFC))
                 * (1.0 - WGOAL*ABSTG)
                 / (1.0 + AVEXP);
double DPD_END = DPDMAX * POSPART(DRIVE) / (DR50 + POSPART(DRIVE));
double DPD     = (FORCEDPD >= 0.0) ? FORCEDPD : DPD_END;

double SEMFR = SEMnM/(EC50SEM + SEMnM);
double TPMFR = CTPM /(EC50TPM + CTPM);
double BACFR = CBAC /(EC50BAC + CBAC);
double ONDFR = HTTF*COND/(EC50OND + COND);
double PRIMEF = 1.0 + WPRIME*(MORACT - MOR0)
                - ESEMINT*SEMFR - ETPMINT*TPMFR - EBACINT*BACFR - EONDINT*ONDFR;
if(PRIMEF < 0.08) PRIMEF = 0.08;

double tod  = SOLVERTIME - 24.0*floor(SOLVERTIME/24.0);
double wgt  = exp(-0.5*pow((tod - TDMID)/TDSIG, 2.0)) / (TDSIG*2.5066283);
double DPD_ACT   = DPD * PRIMEF * DRINKSW;
double RATE_ETH  = DPD_ACT * GDRINK * wgt;         // g/h into stomach

// =====================================================================
// 3. ethanol ADME
// =====================================================================
double VADH  = VMADH*ADHF   * CETH / (KMADH + CETH);
double VMEOS = VMCYP*CYP2E1 * CETH / (KMCYP + CETH);
double VCAT  = VMCAT        * CETH / (KMCAT + CETH);
double VOX   = VADH + VMEOS + VCAT;                // g/h

dxdt_ETH_ST  = RATE_ETH - KEMPT*ETH_ST - KGAST*GASTF*ETH_ST;
dxdt_ETH_GUT = KEMPT*ETH_ST - KAETH*ETH_GUT;
dxdt_ETH_C   = KAETH*ETH_GUT - VOX - KRENETH*AETH;
dxdt_ETOHCUM = RATE_ETH;

// CYP2E1 induction (protein stabilisation)
double CYPTGT = 1.0 + EMAXCYP*CETH/(EC50CYP + CETH);
dxdt_CYP2E1  = KOUTCYP*(CYPTGT - CYP2E1);

// acetaldehyde: production from every oxidative route
double PRODACD = VOX / MWETH * 1e6 / VD;                      // uM/h
double DSFINH  = 1.0 - ALDH2A;                                // fraction inactivated
double VMA2    = VMALD2 * ALDHF * ALDH2A;
double VMA1    = VMALD1 * (1.0 - DSFALD1*DSFINH);
double CACD    = POSPART(ACD);
double ELIMACD = VMA2*CACD/(KMALD2 + CACD) + VMA1*CACD/(KMALD1 + CACD);
dxdt_ACD = PRODACD - ELIMACD;
dxdt_ACT = ELIMACD/1000.0 - 3.0*ACT;

// redox
double NADTGT = 1.0 + ENADH*VOX/(EC50NAD + VOX);
dxdt_NADH = KNADH*(NADTGT - NADH);
double LACTGT = LACT0 + ELACT*(NADH - 1.0)/ENADH;
dxdt_LACT = KLACT*(LACTGT - LACT);

// chronic exposure integrator
dxdt_ETHAVG = (CETH - ETHAVG)/TAUAVG;

// =====================================================================
// 4. objective biomarkers
// =====================================================================
dxdt_PETH = KPETHF*CETH - KPETHE*PETH;
double CDTDRV = POSPART(ETHAVG - CDTHR);
dxdt_CDT  = KCDT*(CDT0 + ECDT*CDTDRV/(EC50CDT + CDTDRV) - CDT);
dxdt_GGT  = KGGT*(GGT0 + EGGT*ETHAVG/(EC50GGT + ETHAVG) - GGT);
dxdt_MCV  = KMCV*(MCV0 + EMCV*ETHAVG/(EC50MCV + ETHAVG) - MCV);

// chronic-exposure driver of every slow neuroadaptation: adaptation
// responds to SUSTAINED exposure, not to a single evening. This is what
// separates a stable social drinker from an escalating one.
double CHRON = ETHAVG*ETHAVG/(EC50CHR*EC50CHR + ETHAVG*ETHAVG);

// =====================================================================
// 5. GABA-A signalling and tolerance
// =====================================================================
double GAPOT = EMAXGA*ETHMM/(EC50GA + ETHMM) / (1.0 + WSUBTOL*GABAA_SUB);
double GTTGT = 1.0 + GAPOT + BZDEFF + TPMGA - WGTSUB*GABAA_SUB;
dxdt_GABA_TONE = KGT*(GTTGT - GABA_TONE);
dxdt_GABAA_SUB = KSUBON*VULN*CHRON*(SUBMAX - GABAA_SUB) - KSUBOFF*GABAA_SUB;

// =====================================================================
// 6. glutamate / NMDA
// =====================================================================
dxdt_NMDA_UP = KNMON*VULN*CHRON*(NMMAX - NMDA_UP) - KNMOFF*NMDA_UP;
dxdt_GLT1    = -KGLTON*CHRON*GLT1 + KGLTOFF*(1.0 + 2.0*ACPEFF)*(1.0 - GLT1);
double GLTCAP = GLT1 < 0.05 ? 0.05 : GLT1;
double WDGATE = NMDA_UP*(1.0 - ETHMM/(8.0 + ETHMM));
double GLUTGT = 1.0 + WGLUT*(1.0/GLTCAP - 1.0) + WGLUWD*WDGATE;
GLUTGT = 1.0 + (GLUTGT - 1.0)*(1.0 - ACPEFF);
dxdt_GLU_NAC = KGLU*(GLUTGT - GLU_NAC);

// =====================================================================
// 7. opioid -> dopamine reinforcement arm
// =====================================================================
double BETGT = 1.0 + EBEND*ETHMM/(EC50BE + ETHMM);
dxdt_BEND = KBEND*(BETGT - BEND);

double DIRDA = EDADIR*ETHMM/(EC50DA + ETHMM)
               * (1.0 - ONDEFF) * (1.0 - SEMDA) * (1.0 - BACEFF) * (1.0 - TPMDA);
double MORPH = MORACT - MOR0;
double DAMOR = MORPH >= 0.0 ? EDAMOR*MORPH : EDAMOR*FTONIC*MORPH;
double DATGT = 1.0 + DAMOR + DIRDA - EDAKOR*POSPART(DYN - 1.0);
if(DATGT < 0.05) DATGT = 0.05;
dxdt_DA_NAC = KDA*(DATGT - DA_NAC);

// =====================================================================
// 8. incentive sensitisation, habit, allostasis
// =====================================================================
double REINF = POSPART(DA_NAC - 1.0);
double SENS  = POSPART(REINF - REINFTH);
dxdt_CUE   = KCUEON*VULN*SENS*SENS*(CUEMAX - CUE)*(1.0 - SEMCUE) - KCUEOFF*CUE;
double DRINKFLAG = ETHMM/(25.0 + ETHMM);
double PFCGATE   = POSPART(1.0 - PFC/PFCH);
dxdt_HABIT = KHABON*DRINKFLAG*PFCGATE*(1.0 - HABIT) - KHABOFF*HABIT;
dxdt_ALLO  = KALLON*VULN*CHRON*(ALLOMAX - ALLO) - KALLOFF*ALLO;

// =====================================================================
// 9. dynorphin / CRF / NPY / noradrenaline
// =====================================================================
dxdt_DYN     = KDYNON*VULN*CHRON*(DYNMAX - DYN)*(1.0 - KOREFF) - KDYNOFF*(DYN - 1.0);
dxdt_CRF_CEA = KCRFON*VULN*CHRON*(CRFMAX - CRF_CEA) - KCRFOFF*(CRF_CEA - 1.0);
dxdt_NPY_CEA = -KNPYON*CHRON*NPY_CEA + KNPYOFF*(1.0 - NPY_CEA);

// CNS excitation index (what benzodiazepines titrate)
double GADEF = POSPART(1.0 - GABA_TONE);
double EXCRAW = WNMEXC*NMDA_UP + WGAEXC*GADEF
                + WNEEXC2*POSPART(NE_LC - 1.0)*(1.0 - A1BLK)
                + WCRFEXC*POSPART(CRF_CEA - 1.0);
double EXC = POSPART(EXCRAW)*(1.0 + WKINDLE*KINDLE);

double NETGT = 1.0 + WNECRF*POSPART(CRF_CEA - 1.0) + WNEEXC*EXCRAW;
dxdt_NE_LC = KNE*(NETGT - NE_LC);

// =====================================================================
// 10. HPA axis
// =====================================================================
double CRHTGT = 1.0 + WCRHAM*POSPART(CRF_CEA - 1.0)
                + WCRHETH*ETHMM/(EC50CRH + ETHMM)
                - WFEED*GR_SENS*POSPART(CORT - 1.0);
if(CRHTGT < 0.05) CRHTGT = 0.05;
dxdt_CRH_HYP = KCRH*(CRHTGT - CRH_HYP);
dxdt_ACTH    = KACTH*(CRH_HYP - ACTH);
dxdt_CORT    = KCORT*(ACTH - CORT);
dxdt_GR_SENS = -KGRON*CHRON*GR_SENS + KGROFF*(1.0 - GR_SENS);

// =====================================================================
// 11. prefrontal control
// =====================================================================
dxdt_PFC = KPFCREC*(1.0 - PFC) - KPFCDAM*CHRON*PFC - WGLUPFC*KPFCREC*POSPART(GLU_NAC - 1.0)*PFC;

// =====================================================================
// 12. negative affect, sleep, craving
// =====================================================================
double NEGTGT = WNCRF*POSPART(CRF_CEA - 1.0)
              + WNDYN*POSPART(DYN - 1.0)
              + WNGLU*POSPART(GLU_NAC - 1.0)
              + WNALLO*ALLO
              + WNSLP*SLEEPD
              - WNNPY*POSPART(1.0 - NPY_CEA);
NEGTGT = POSPART(NEGTGT)*(1.0 - GBPEFF)*(1.0 - 0.85*ACPEFF)*(1.0 - 1.40*TPMGA);
dxdt_NEGAFF = KNEG*(NEGTGT - NEGAFF);

double SLPTGT = POSPART(WSLPEXC*EXC + WSLPNEG*NEGAFF - 1.6*BZDEFF)*(1.0 - GBPEFF);
dxdt_SLEEPD = KSLP*(SLPTGT - SLEEPD);

double CRTGT = WINC*CUE*CUEEXP*(1.0 - SEMCUE)
             + WREL*NEGAFF
             + WHAB*HABIT;
dxdt_CRAVE = KCRAVE*(CRTGT - CRAVE);
dxdt_DRIVE = KDRIVE*(DRIVE_T - DRIVE);

// =====================================================================
// 13. withdrawal and kindling
// =====================================================================
double CIWATGT = CIWAMAX*pow(EXC, HILLEXC)/(pow(EC50EXC, HILLEXC) + pow(EXC, HILLEXC));
dxdt_CIWA   = KCIWA*(CIWATGT - CIWA);
dxdt_KINDLE = KKINON*POSPART(CIWA - CIWAKTH)*(KINDMAX - KINDLE) - KKINOFF*KINDLE;

// =====================================================================
// 14. aversive expectancy (disulfiram arm)
// =====================================================================
double FLUSH = AVMAX*pow(CACD, AVHILL)/(pow(AVAC50, AVHILL) + pow(CACD, AVHILL));
double KNOWN = WKNOW*DSFINH*(0.05 + 0.95*SUPERV);
dxdt_AVEXP = KAVON*POSPART(FLUSH + KNOWN - AVEXP) - KAVOFF*AVEXP;

// =====================================================================
// 15. liver, nutrition, cardiovascular
// =====================================================================
dxdt_STEAT = KSTEAT*(STEAT0 + ESTEAT*ETHAVG/(EC50ST + ETHAVG) - STEAT);
dxdt_LPS   = KLPS*(LPS0 + ELPS*ETHAVG/(EC50LPS + ETHAVG) - LPS);
dxdt_TNFA  = KTNF*(TNF0 + ETNF*POSPART(LPS - 1.0)/(1.5 + POSPART(LPS - 1.0)) - TNFA);
double INJ = POSPART(TNFA - 1.0)/(1.5 + POSPART(TNFA - 1.0));
double FATX = POSPART(STEAT - STEAT0)/ESTEAT;
dxdt_ALT   = KALT*(ALT0 + EALT*(0.45*FATX + 0.55*INJ) - ALT);
dxdt_AST   = KAST*(AST0 + EAST*(0.35*FATX + 0.65*INJ) - AST);
dxdt_FIB   = KFIBON*POSPART(TNFA - 1.0)*(FIBMAX - FIB) - KFIBOFF*FIB;
dxdt_THIA  = KTHIA*(THIA0*(1.0 - ETHIA*ETHAVG/(EC50TH + ETHAVG)) - THIA) + THIADOSE;
dxdt_SBP   = KSBP*(SBP0 + ESBP*ETHAVG/(EC50SBP + ETHAVG) - SBP);
dxdt_CMYO  = KCMYO*ETHAVG - KCMYOR*CMYO;

// =====================================================================
// 16. drug pharmacokinetics
// =====================================================================
dxdt_NTX_GUT = -KANTX*NTX_GUT;
dxdt_NTX_C   =  KANTX*NTX_GUT*FNTX + KRELDEP*NTX_DEP - CLNTX/VNTX*NTX_C;
dxdt_NTXOL   =  KANTX*NTX_GUT*FMNOL
                + FSYSNOL*CLNTX/VNTX*NTX_C - CLNOL/VNOL*NTXOL;
dxdt_NTX_DEP = -KRELDEP*NTX_DEP;

dxdt_ACP_GUT = -KAACP*ACP_GUT;
dxdt_ACP_C   =  KAACP*ACP_GUT*FACP - CLACP/VACP*ACP_C;

dxdt_DSF_GUT = -KADSF*DSF_GUT;
dxdt_DSF_C   =  KADSF*DSF_GUT*FDSF - CLDSF/VDSF*DSF_C;
dxdt_ALDH2A  =  KSYNALD*(1.0 - ALDH2A) - KINACT*CDSF*ALDH2A;

dxdt_TPM_GUT = -KATPM*TPM_GUT;
dxdt_TPM_C   =  KATPM*TPM_GUT*FTPM - CLTPM/VTPM*TPM_C;
dxdt_GBP_GUT = -KAGBP*GBP_GUT;
dxdt_GBP_C   =  KAGBP*GBP_GUT*FGBP - CLGBP/VGBP*GBP_C;
dxdt_BACL_GUT= -KABAC*BACL_GUT;
dxdt_BACL_C  =  KABAC*BACL_GUT*FBAC - CLBAC/VBAC*BACL_C;
dxdt_NMF_GUT = -KANMF*NMF_GUT;
dxdt_NMF_C   =  KANMF*NMF_GUT*FNMF - CLNMF/VNMF*NMF_C;
dxdt_OND_GUT = -KAOND*OND_GUT;
dxdt_OND_C   =  KAOND*OND_GUT*FOND - CLOND/VOND*OND_C;
dxdt_SEM_SC  = -KASEM*SEM_SC;
dxdt_SEM_C   =  KASEM*SEM_SC*FSEM - CLSEM/VSEM*SEM_C;

// symptom-triggered benzodiazepine (continuous analogue of a CIWA protocol)
double BZDRATE = (BZDON > 0.5 && SOLVERTIME <= BZDDUR) ? BZDGAIN*POSPART(CIWA - BZDTHR) : 0.0;
if(BZDRATE > BZDMAXR) BZDRATE = BZDMAXR;
dxdt_BZD_C  = BZDRATE - CLBZD/VBZD*BZD_C;
dxdt_BZDCUM = BZDRATE;

$TABLE
double RWIDT  = SEXF > 0.5 ? RWIDF : RWIDM;
double VDT    = RWIDT * BW;
double CETHO  = POSPART(ETH_C)/VDT;
double BACO   = CETHO/10.0;
double ETHMMO = CETHO/0.04607;
double GPD    = DPD_ACT * GDRINK;              // grams ethanol per day
double MOROCC = 100.0*ANTAG/(1.0 + ANTAG);     // % MOR occupied by antagonist
double ASTALT = AST/ALT;
double SEIZP  = 100.0/(1.0 + exp(-(0.72*EXC - 5.9)));
double WHOLVL = GPD < 1.0  ? 0.0 :
                (GPD < 40.0 ? 1.0 : (GPD < 60.0 ? 2.0 : (GPD < 100.0 ? 3.0 : 4.0)));

$CAPTURE @annotated
BACO   : blood alcohol concentration (g/dL)
CETHO  : plasma ethanol (g/L)
ETHMMO : brain/plasma ethanol (mM)
DPD_ACT: drinks per day actually consumed
GPD    : grams ethanol per day
DPD_END: endogenous drinks-per-day set point
DRIVE_T: instantaneous drinking-drive target
MORACT : mu-opioid receptor activation (fraction)
MOROCC : MOR occupancy by antagonist (%)
EXC    : CNS excitation index
CIWATGT: CIWA target
ASTALT : AST/ALT ratio
SEIZP  : modelled withdrawal-seizure probability (%)
WHOLVL : WHO risk drinking level (0-4)
CACD   : blood acetaldehyde (uM)
FLUSH  : acetaldehyde flush aversion signal
CNTX   : naltrexone (ng/mL)
CNOL   : 6-beta-naltrexol (ng/mL)
CACP   : acamprosate (ng/mL)
CDSF   : disulfiram + metabolites (mg/L)
CTPM   : topiramate (mg/L)
CGBP   : gabapentin (mg/L)
CBAC   : baclofen (ng/mL)
CNMF   : nalmefene (ng/mL)
CBZD   : benzodiazepine (mg/L diazepam-equivalent)
SEMnM  : semaglutide (nM)
VOX    : total ethanol oxidation flux (g/h)
'

mod <- mcode_cache("aud_qsp", code, soloc = tempdir())

###############################################################################
## HELPERS
###############################################################################

DAY <- 24

## regimen builders ---------------------------------------------------------
ev_daily <- function(cmt, amt, start_d, dur_d, ii = 24, addl_extra = 0) {
  ev(time = start_d*DAY, amt = amt, cmt = cmt, ii = ii,
     addl = max(0, floor(dur_d*DAY/ii) - 1 + addl_extra))
}

## run a scenario -----------------------------------------------------------
run_scn <- function(param = list(), init = NULL, events = NULL,
                    end_d = 180, delta = 0.5) {
  m <- mod
  if (length(param)) m <- param(m, param)
  if (!is.null(init)) m <- init(m, init)
  if (is.null(events)) {
    m %>% mrgsim(end = end_d*DAY, delta = delta, hmax = 0.5,
                 maxsteps = 500000) %>% as_tibble()
  } else {
    m %>% mrgsim(events = events, end = end_d*DAY, delta = delta,
                 hmax = 0.5, maxsteps = 500000) %>% as_tibble()
  }
}

## daily summary ------------------------------------------------------------
daily <- function(d) {
  d %>%
    mutate(day = floor(time/DAY)) %>%
    group_by(day) %>%
    summarise(
      GPD     = mean(GPD),
      DPD     = mean(DPD_ACT),
      BACmax  = max(BACO),
      CIWAmax = max(CIWA),
      CRAVE   = mean(CRAVE),
      NEGAFF  = mean(NEGAFF),
      PETH    = mean(PETH),
      GGT     = mean(GGT),
      CDT     = mean(CDT),
      MCV     = mean(MCV),
      ALT     = mean(ALT),
      AST     = mean(AST),
      SBP     = mean(SBP),
      STEAT   = mean(STEAT),
      GABAA   = mean(GABAA_SUB),
      NMDAUP  = mean(NMDA_UP),
      GLT1    = mean(GLT1),
      DYN     = mean(DYN),
      CRF     = mean(CRF_CEA),
      CUE     = mean(CUE),
      HABIT   = mean(HABIT),
      PFC     = mean(PFC),
      ALLO    = mean(ALLO),
      KINDLE  = mean(KINDLE),
      MOROCC  = mean(MOROCC),
      SEIZP   = max(SEIZP),
      BZDCUM  = max(BZDCUM),
      THIA    = mean(THIA),
      .groups = "drop")
}

###############################################################################
## STATISTICAL LAYER: deterministic mean consumption -> trial endpoints
## ---------------------------------------------------------------------------
## The ODE system returns a single deterministic mean consumption trajectory.
## Real cohorts do not drink the mean every day: they alternate abstinent days
## and drinking days, and the trial endpoints (%HDD, PDA, DPDD) are properties
## of that DISTRIBUTION, not of the mean. We therefore make the day-level
## distribution explicit and auditable rather than hiding it:
##
##   day i is abstinent with probability  pi
##   otherwise  drinks_i ~ zero-truncated negative binomial(mu = DPDD, size = k)
##
## with
##   pi  = plogis(Z0 + Z1*ABSTG + Z2*dNEG + Z3*dCRAVE)
##   DPDD = mu_model / (1 - pi)
##
## where dNEG and dCRAVE are the FRACTIONAL reductions in negative affect and
## craving relative to the matched untreated/placebo arm. This is the whole
## point of the two-arm architecture: agents that cut the negative-reinforcement
## arm raise pi (percent days abstinent) while agents that cut the
## positive-reinforcement arm lower mu_model (drinks per drinking day and %HDD).
## Nothing else in the mapping is drug-specific.
###############################################################################

Z0 <- -1.35; Z1 <- 1.60; Z2 <- 0.60; Z3 <- 0.12; NBK <- 2.4

ztnb_ge <- function(thresh, mu, size) {
  ## P(X >= thresh | X >= 1) for a negative binomial
  p0 <- dnbinom(0, mu = mu, size = size)
  pl <- pnbinom(thresh - 1, mu = mu, size = size)
  max(0, (1 - pl)/(1 - p0))
}

endpoints <- function(dsum, ref = NULL, abstg = 0, window = NULL, hdd_thr = 5) {
  if (!is.null(window)) dsum <- dsum %>% filter(day >= window[1], day <= window[2])
  mu     <- mean(dsum$DPD)
  neg    <- mean(dsum$NEGAFF)
  crv    <- mean(dsum$CRAVE)
  if (is.null(ref)) { dneg <- 0; dcrv <- 0 } else {
    dneg <- (ref$neg - neg)/max(ref$neg, 1e-6)
    dcrv <- (ref$crv - crv)/max(ref$crv, 1e-6)
  }
  pi_ab <- plogis(Z0 + Z1*abstg + Z2*dneg + Z3*dcrv)
  pi_ab <- min(max(pi_ab, 0.001), 0.999)
  dpdd  <- mu/(1 - pi_ab)
  phdd  <- 100*(1 - pi_ab)*ztnb_ge(hdd_thr, dpdd, NBK)
  list(mu = mu, neg = neg, crv = crv, pi = pi_ab,
       PDA = 100*pi_ab, DPDD = dpdd, HDD = phdd,
       PSNHD = 100*exp(-phdd/100*30),   # P(no heavy day in a 30-day month)
       TAC = mu*14)
}

###############################################################################
## STEP 1 — HEALTHY BASELINE: does the drug-free model actually sit still?
###############################################################################

cat("\n", strrep("=", 78), "\n", sep = "")
cat("ALCOHOL USE DISORDER QSP MODEL — VALIDATION REPORT\n")
cat(strrep("=", 78), "\n", sep = "")

cat("\n[1] DRY BASELINE — no ethanol, no drug, 180 days\n")
dry <- run_scn(param = list(ENVDRIVE = 0, DRINKSW = 0), end_d = 180, delta = 6)
b_states <- c("GABA_TONE","GLU_NAC","DA_NAC","BEND","NE_LC","CRAVE","GABAA_SUB",
              "NMDA_UP","GLT1","DYN","CRF_CEA","NPY_CEA","PFC","CORT","NEGAFF",
              "CIWA","GGT","CDT","MCV","ALT","AST","SBP","STEAT","THIA")
drift <- sapply(b_states, function(s) {
  v <- dry[[s]]; d0 <- v[1]; d1 <- v[length(v)]
  if (abs(d0) < 1e-9) return(abs(d1))
  100*(d1 - d0)/d0
})
cat(sprintf("  max |drift| over 180 d = %.4f %%   (n = %d states)\n",
            max(abs(drift)), length(drift)))
print(round(drift[order(-abs(drift))][1:6], 5))

###############################################################################
## STEP 2 — ETHANOL PK: single 0.8 g/kg challenge (Widmark validation)
###############################################################################

cat("\n[2] ETHANOL PK — single 0.8 g/kg oral challenge in a 75 kg male\n")
chal <- run_scn(param = list(ENVDRIVE = 0, DRINKSW = 0, FED = 0.2),
                events = ev(time = 12, amt = 0.8*75, cmt = "ETH_ST"),
                end_d = 1, delta = 0.05)
pk <- chal %>% filter(time >= 12)
cmax  <- max(pk$BACO)
tmax  <- pk$time[which.max(pk$BACO)] - 12
## descending-limb slope between 40% and 80% of the decline
desc  <- pk %>% filter(time > pk$time[which.max(pk$BACO)], BACO > 0.02)
slope <- if (nrow(desc) > 10) -coef(lm(BACO ~ time, data = desc))[2] else NA
acdmax <- max(pk$CACD)
cat(sprintf("  Cmax BAC       = %.4f g/dL      (expected 0.075-0.105)\n", cmax))
cat(sprintf("  Tmax           = %.2f h          (expected 0.5-1.2)\n", tmax))
cat(sprintf("  beta-slope     = %.4f g/dL/h    (expected 0.012-0.020)\n", slope))
cat(sprintf("  peak acetaldehyde = %.2f uM      (expected 1-3 in ALDH2 wild-type)\n", acdmax))

## ALDH2*2 heterozygote and disulfiram comparison on the same challenge
chal_a2 <- run_scn(param = list(ENVDRIVE = 0, DRINKSW = 0, FED = 0.2, ALDHF = 0.20),
                   events = ev(time = 12, amt = 0.8*75, cmt = "ETH_ST"),
                   end_d = 1, delta = 0.05)
cat(sprintf("  ALDH2*2 het peak acetaldehyde = %.1f uM   (expected 5-20x wild-type)\n",
            max(chal_a2$CACD)))

###############################################################################
## STEP 3 — DISEASE GENERATION: 3-year escalation from social drinking
###############################################################################

cat("\n[3] DISEASE GENERATION — 3 years, cue-rich environment, vulnerable host\n")

esc <- run_scn(param = list(ENVDRIVE = 4.6, VULN = 1.35),
               end_d = 3*365 + 21/24, delta = 3)
escd <- daily(esc)
mk <- function(d) escd %>% filter(day == d)
for (d in c(30, 180, 365, 730, 1095)) {
  r <- mk(d)
  cat(sprintf("  day %4d : %5.1f drinks/d  %5.0f g/d | GABAAsub %.2f  NMDAup %.2f  GLT1 %.2f | CUE %5.2f HABIT %.2f PFC %.2f | GGT %5.0f PEth %5.0f\n",
              d, r$DPD, r$GPD, r$GABAA, r$NMDAUP, r$GLT1, r$CUE, r$HABIT, r$PFC, r$GGT, r$PETH))
}

## the severe-AUD state (detox and natural-course scenarios start here)
AUD_INIT <- as.list(esc[nrow(esc), names(init(mod))])
AUD_INIT$ETOHCUM <- 0; AUD_INIT$BZDCUM <- 0   # patient presents still intoxicated

cat(sprintf("\n  SEVERE-AUD entry state: %.1f drinks/day (%.0f g/d), PEth %.0f ng/mL, GGT %.0f U/L, CDT %.2f%%, MCV %.1f fL\n",
            tail(escd$DPD,1), tail(escd$GPD,1), tail(escd$PETH,1),
            tail(escd$GGT,1), tail(escd$CDT,1), tail(escd$MCV,1)))
cat(sprintf("  AST/ALT %.2f (%.0f/%.0f U/L), steatosis %.0f%%, SBP %.0f mmHg, thiamine %.2f\n",
            tail(escd$AST,1)/tail(escd$ALT,1), tail(escd$AST,1), tail(escd$ALT,1),
            100*tail(escd$STEAT,1), tail(escd$SBP,1), tail(escd$THIA,1)))

## a second, less extreme phenotype matched to what medication trials enrol
## (COMBINE / Jonas 2014 entrants average ~5-6 drinks/day, ~65% drinking days)
esct  <- run_scn(param = list(ENVDRIVE = 3.40, VULN = 1.15),
                 end_d = 3*365 + 21/24, delta = 3)
esctd <- daily(esct)
TRIAL_INIT <- as.list(esct[nrow(esct), names(init(mod))])
TRIAL_INIT$ETOHCUM <- 0; TRIAL_INIT$BZDCUM <- 0
cat(sprintf("  TRIAL-POPULATION entry state: %.1f drinks/day (%.0f g/d), PEth %.0f, GGT %.0f, CDT %.2f%%, GABAAsub %.2f\n",
            tail(esctd$DPD,1), tail(esctd$GPD,1), tail(esctd$PETH,1),
            tail(esctd$GGT,1), tail(esctd$CDT,1), tail(esctd$GABAA,1)))

###############################################################################
## STEP 4 — HANDOVER TEST: is the closed loop self-consistent?
###############################################################################

cat("\n[4] HANDOVER TEST — release the endogenous loop, no treatment, 90 days\n")
hand <- run_scn(init = AUD_INIT, param = list(ENVDRIVE = 4.6, VULN = 1.35),
                end_d = 90, delta = 1)
handd <- daily(hand)
cat(sprintf("  day 0 %.2f -> day 90 %.2f drinks/day  (drift %+.1f %%)\n",
            handd$DPD[1], tail(handd$DPD,1),
            100*(tail(handd$DPD,1) - handd$DPD[1])/handd$DPD[1]))

###############################################################################
## STEP 5 — TREATMENT SCENARIOS
###############################################################################

AUDP  <- list(ENVDRIVE = 4.6, VULN = 1.35)   # severe phenotype
TRLP  <- list(ENVDRIVE = 3.40, VULN = 1.15)  # medication-trial phenotype
TRT_D <- 180                                   # 6-month treatment period
WIN   <- c(30, 180)                            # endpoint window (weeks 5-26)

## ---- detox skeleton: enforced abstinence for the first 5 inpatient days ----
detox_par <- function(extra = list(), bzd = 0) {
  modifyList(c(AUDP, list(ABSTG = 0.55, BZDON = bzd)), extra)
}

scenarios <- list()

## S01 healthy social drinker (self-calibration reference) -------------------
scenarios$S01 <- list(
  name = "S01 stable moderate drinker (no vulnerability, low-cue environment)",
  par  = list(ENVDRIVE = 1.6, VULN = 1.0), init = NULL, ev = NULL, abstg = 0)

## S02 untreated severe AUD --------------------------------------------------
scenarios$S02 <- list(
  name = "S02 untreated severe AUD (natural course)",
  par  = AUDP, init = AUD_INIT, ev = NULL, abstg = 0)

## S03 medical management only (psychosocial platform) -----------------------
scenarios$S03 <- list(
  name = "S03 medical management / CBT only (placebo arm of COMBINE)",
  par  = c(TRLP, list(ABSTG = 0.55)), init = TRIAL_INIT, ev = NULL, abstg = 0.55)

## S04 unmedicated abrupt cessation ------------------------------------------
scenarios$S04 <- list(
  name = "S04 abrupt cessation, NO withdrawal treatment",
  par  = c(AUDP, list(ABSTG = 0.55, DRINKSW = 0)), init = AUD_INIT, ev = NULL,
  abstg = 0.55, detox = TRUE)

## S05 symptom-triggered benzodiazepine detox --------------------------------
scenarios$S05 <- list(
  name = "S05 symptom-triggered benzodiazepine detox (CIWA-Ar >= 8)",
  par  = c(AUDP, list(ABSTG = 0.55, DRINKSW = 0, BZDON = 1, THIADOSE = 0.02)),
  init = AUD_INIT, ev = NULL, abstg = 0.55, detox = TRUE)

## S06 naltrexone PO 50 mg ---------------------------------------------------
scenarios$S06 <- list(
  name = "S06 naltrexone 50 mg PO daily + MM",
  par  = c(TRLP, list(ABSTG = 0.55)), init = TRIAL_INIT,
  ev   = ev_daily("NTX_GUT", 50, 0, TRT_D), abstg = 0.55)

## S07 XR-naltrexone 380 mg IM q4wk -----------------------------------------
scenarios$S07 <- list(
  name = "S07 XR-naltrexone 380 mg IM q4wk + MM",
  par  = c(TRLP, list(ABSTG = 0.55)), init = TRIAL_INIT,
  ev   = ev(time = 0, amt = 380, cmt = "NTX_DEP", ii = 28*DAY, addl = 5),
  abstg = 0.55)

## S08 acamprosate 666 mg tid ------------------------------------------------
scenarios$S08 <- list(
  name = "S08 acamprosate 666 mg tid + MM",
  par  = c(TRLP, list(ABSTG = 0.55)), init = TRIAL_INIT,
  ev   = ev(time = 0, amt = 666, cmt = "ACP_GUT", ii = 8, addl = TRT_D*3 - 1),
  abstg = 0.55)

## S09 naltrexone + acamprosate (COMBINE cell) -------------------------------
scenarios$S09 <- list(
  name = "S09 naltrexone + acamprosate + MM (COMBINE combined cell)",
  par  = c(TRLP, list(ABSTG = 0.55)), init = TRIAL_INIT,
  ev   = ev_daily("NTX_GUT", 50, 0, TRT_D) +
         ev(time = 0, amt = 666, cmt = "ACP_GUT", ii = 8, addl = TRT_D*3 - 1),
  abstg = 0.55)

## S10 disulfiram supervised -------------------------------------------------
scenarios$S10 <- list(
  name = "S10 disulfiram 250 mg daily, SUPERVISED",
  par  = c(TRLP, list(ABSTG = 0.55, SUPERV = 1)), init = TRIAL_INIT,
  ev   = ev_daily("DSF_GUT", 250, 0, TRT_D), abstg = 0.55)

## S11 disulfiram unsupervised (50% adherence) -------------------------------
scenarios$S11 <- list(
  name = "S11 disulfiram 250 mg, UNSUPERVISED (~25% adherence)",
  par  = c(TRLP, list(ABSTG = 0.55, SUPERV = 0)), init = TRIAL_INIT,
  ev   = ev(time = 0, amt = 250, cmt = "DSF_GUT", ii = 96, addl = TRT_D/4 - 1),
  abstg = 0.55)

## S12 topiramate 300 mg/d ---------------------------------------------------
scenarios$S12 <- list(
  name = "S12 topiramate titrated to 300 mg/day",
  par  = c(TRLP, list(ABSTG = 0.55)), init = TRIAL_INIT,
  ev   = ev(time = 0, amt = 50, cmt = "TPM_GUT", ii = 24, addl = 13) +
         ev(time = 14*DAY, amt = 150, cmt = "TPM_GUT", ii = 24, addl = 13) +
         ev(time = 28*DAY, amt = 300, cmt = "TPM_GUT", ii = 24, addl = TRT_D - 29),
  abstg = 0.55)

## S13 gabapentin 1800 mg/d --------------------------------------------------
scenarios$S13 <- list(
  name = "S13 gabapentin 600 mg tid (1800 mg/day)",
  par  = c(TRLP, list(ABSTG = 0.55)), init = TRIAL_INIT,
  ev   = ev(time = 0, amt = 600, cmt = "GBP_GUT", ii = 8, addl = TRT_D*3 - 1),
  abstg = 0.55)

## S14 nalmefene as-needed ---------------------------------------------------
scenarios$S14 <- list(
  name = "S14 nalmefene 18 mg as-needed (targeted reduction)",
  par  = c(TRLP, list(ABSTG = 0.20)), init = TRIAL_INIT,
  ev   = ev(time = 14, amt = 18, cmt = "NMF_GUT", ii = 24, addl = TRT_D - 1),
  abstg = 0.20)

## S14R matched reduction-goal control for the nalmefene comparison ----------
scenarios$S14R <- list(
  name = "S14R reduction-goal control (no medication, ABSTG 0.20)",
  par  = c(TRLP, list(ABSTG = 0.20)), init = TRIAL_INIT, ev = NULL, abstg = 0.20)

## S15 baclofen high dose ----------------------------------------------------
scenarios$S15 <- list(
  name = "S15 baclofen 60 mg tid (180 mg/day, high-dose protocol)",
  par  = c(TRLP, list(ABSTG = 0.55)), init = TRIAL_INIT,
  ev   = ev(time = 0, amt = 60, cmt = "BACL_GUT", ii = 8, addl = TRT_D*3 - 1),
  abstg = 0.55)

## S16 semaglutide -----------------------------------------------------------
scenarios$S16 <- list(
  name = "S16 semaglutide 0.25 -> 1.0 mg SC weekly",
  par  = c(TRLP, list(ABSTG = 0.55)), init = TRIAL_INIT,
  ev   = ev(time = 0, amt = 0.25, cmt = "SEM_SC", ii = 168, addl = 3) +
         ev(time = 4*168, amt = 0.5, cmt = "SEM_SC", ii = 168, addl = 3) +
         ev(time = 8*168, amt = 1.0, cmt = "SEM_SC", ii = 168,
            addl = floor(TRT_D*24/168) - 9),
  abstg = 0.55)

## S17 ondansetron in the 5-HTTLPR LL / early-onset subtype ------------------
scenarios$S17 <- list(
  name = "S17 ondansetron 4 ug/kg bid, 5-HTTLPR LL early-onset subtype",
  par  = c(TRLP, list(ABSTG = 0.55, HTTF = 1.8)), init = TRIAL_INIT,
  ev   = ev(time = 0, amt = 0.004*75, cmt = "OND_GUT", ii = 12,
            addl = TRT_D*2 - 1), abstg = 0.55)

## S18 naltrexone in an OPRM1 118G carrier -----------------------------------
scenarios$S18 <- list(
  name = "S18 naltrexone 50 mg PO in an OPRM1 118G carrier",
  par  = c(TRLP, list(ABSTG = 0.55, OPRMF = 1.6, WPRIME = 0.85)),
  init = TRIAL_INIT, ev = ev_daily("NTX_GUT", 50, 0, TRT_D), abstg = 0.55)

## S19 ALDH2*2 heterozygote (the endogenous disulfiram) ----------------------
scenarios$S19 <- list(
  name = "S19 ALDH2*2 heterozygote, same environment, no medication",
  par  = c(AUDP, list(ALDHF = 0.20)), init = NULL, ev = NULL, abstg = 0)

## S20 kindled patient: 5 prior detoxifications ------------------------------
KIND_INIT <- AUD_INIT; KIND_INIT$KINDLE <- 1.6
scenarios$S20 <- list(
  name = "S20 kindled patient (5 prior detoxifications), unmedicated cessation",
  par  = c(AUDP, list(ABSTG = 0.55, DRINKSW = 0)), init = KIND_INIT, ev = NULL,
  abstg = 0.55, detox = TRUE)

###############################################################################
## RUN
###############################################################################

cat("\n[5] RUNNING ", length(scenarios), " SCENARIOS\n", sep = "")
res <- list()
for (nm in names(scenarios)) {
  s <- scenarios[[nm]]
  end_d <- if (isTRUE(s$detox)) 21 else TRT_D
  dd <- if (isTRUE(s$detox)) 0.25 else 0.5
  out <- run_scn(param = s$par, init = s$init, events = s$ev,
                 end_d = end_d, delta = dd)
  res[[nm]] <- list(name = s$name, raw = out, day = daily(out), abstg = s$abstg,
                    detox = isTRUE(s$detox))
  cat(sprintf("  %s ... ok\n", nm))
}

###############################################################################
## STEP 6 — WITHDRAWAL VALIDATION
###############################################################################

cat("\n[6] WITHDRAWAL — 21-day inpatient course\n")
cat(sprintf("  %-58s %7s %7s %7s %9s\n", "scenario", "peakCIWA", "t_peak", "seiz%", "diazepam"))
for (nm in c("S04","S05","S20")) {
  r <- res[[nm]]$raw
  pc <- max(r$CIWA); tp <- r$time[which.max(r$CIWA)]
  cat(sprintf("  %-58s %7.1f %6.1fh %6.1f%% %8.0f mg\n",
              res[[nm]]$name, pc, tp, max(r$SEIZP), max(r$BZDCUM)))
}

###############################################################################
## STEP 7 — TRIAL ENDPOINTS
###############################################################################

cat("\n[7] SIX-MONTH TRIAL ENDPOINTS (weeks 5-26)\n")

trt <- setdiff(names(scenarios), c("S04","S05","S20"))
refE <- endpoints(res$S03$day, ref = NULL, abstg = 0.55, window = WIN)
ref  <- list(neg = refE$neg, crv = refE$crv)

EP <- lapply(trt, function(nm) {
  e <- endpoints(res[[nm]]$day, ref = ref, abstg = res[[nm]]$abstg, window = WIN)
  data.frame(scn = nm, name = res[[nm]]$name,
             gpd = e$mu*14, DPDD = e$DPDD, HDD = e$HDD,
             PDA = e$PDA, PSNHD = e$PSNHD,
             CRAVE = e$crv, NEGAFF = e$neg)
})
EP <- bind_rows(EP)

cat(sprintf("  %-5s %28s %7s %7s %7s %7s %7s %7s\n",
            "scn", "arm", "g/day", "DPDD", "%HDD", "%PDA", "PSNHD", "craving"))
for (i in seq_len(nrow(EP))) {
  lbl <- substr(sub("^S[0-9]+ ", "", EP$name[i]), 1, 28)
  cat(sprintf("  %-5s %28s %7.1f %7.2f %7.1f %7.1f %7.1f %7.2f\n",
              EP$scn[i], lbl, EP$gpd[i], EP$DPDD[i], EP$HDD[i],
              EP$PDA[i], EP$PSNHD[i], EP$CRAVE[i]))
}

## deltas versus the medical-management-only arm ----------------------------
base <- EP %>% filter(scn == "S03")
cat("\n[8] DELTA VERSUS MEDICAL MANAGEMENT ALONE (S03) — the clinical signature\n")
cat(sprintf("  %-5s %28s %9s %9s %9s\n", "scn", "arm", "d%HDD", "dDPDD", "d%PDA"))
for (i in seq_len(nrow(EP))) {
  if (EP$scn[i] %in% c("S01","S02","S03","S19","S14R")) next
  ## nalmefene is licensed for REDUCTION, not abstinence: compare it with the
  ## matched reduction-goal control (S14R), not with an abstinence-goal arm
  cmp <- if (EP$scn[i] == "S14") EP %>% filter(scn == "S14R") else base
  tag <- if (EP$scn[i] == "S14") " (vs S14R)" else ""
  cat(sprintf("  %-5s %28s %+9.1f %+9.2f %+9.1f%s\n",
              EP$scn[i], substr(sub("^S[0-9]+ ", "", EP$name[i]), 1, 28),
              EP$HDD[i] - cmp$HDD, EP$DPDD[i] - cmp$DPDD,
              EP$PDA[i] - cmp$PDA, tag))
}

cat("\n  Reference effect sizes from the literature (vs placebo, 12-24 weeks):\n")
cat("    naltrexone 50 mg    %HDD  -4 to -6   DPDD -0.5    PDA  +4    Jonas 2014 JAMA\n")
cat("    acamprosate         %HDD   ~ 0       DPDD  ~ 0    PDA  +9    Jonas 2014 JAMA\n")
cat("    topiramate 300 mg   %HDD  -8 to -9   DPDD -1.0    PDA +10    Johnson 2007 JAMA\n")
cat("    gabapentin 1800 mg  %HDD  -9         DPDD -1.0    PDA +13    Mason 2014 JAMA IM\n")
cat("    nalmefene 18 mg prn HDD  -1.7 to -2.7 d/mo        TAC -11 g/d ESENSE1/2\n")
cat("    XR-naltrexone       heavy-drinking events -25%              Garbutt 2005 JAMA\n")

###############################################################################
## STEP 9 — BIOMARKER AND ORGAN RECOVERY
###############################################################################

cat("\n[9] BIOMARKER / ORGAN TRAJECTORY AT 6 MONTHS\n")
cat(sprintf("  %-5s %26s %7s %7s %6s %6s %6s %6s %6s\n",
            "scn", "arm", "PEth", "GGT", "CDT", "MCV", "AST", "ALT", "SBP"))
for (nm in trt) {
  d <- res[[nm]]$day; last <- tail(d, 1)
  cat(sprintf("  %-5s %26s %7.0f %7.0f %6.2f %6.1f %6.0f %6.0f %6.0f\n",
              nm, substr(sub("^S[0-9]+ ", "", res[[nm]]$name), 1, 26),
              last$PETH, last$GGT, last$CDT, last$MCV, last$AST, last$ALT, last$SBP))
}

###############################################################################
## STEP 10 — WHERE EACH DRUG CUTS THE LOOP
###############################################################################

cat("\n[10] MECHANISTIC READ-OUT AT 6 MONTHS (why the endpoints move)\n")
cat(sprintf("  %-5s %24s %7s %7s %7s %7s %7s %7s\n",
            "scn", "arm", "MOROCC", "DA_NAC", "CRAVE", "NEGAFF", "CUE", "GLT1"))
for (nm in trt) {
  r <- res[[nm]]$raw; w <- r %>% filter(time >= WIN[1]*24)
  cat(sprintf("  %-5s %24s %7.1f %7.3f %7.2f %7.2f %7.2f %7.2f\n",
              nm, substr(sub("^S[0-9]+ ", "", res[[nm]]$name), 1, 24),
              mean(w$MOROCC), mean(w$DA_NAC), mean(w$CRAVE),
              mean(w$NEGAFF), mean(w$CUE), mean(w$GLT1)))
}

###############################################################################
## STEP 11 — DRUG EXPOSURE CHECK
###############################################################################

cat("\n[11] DRUG EXPOSURE vs PUBLISHED VALUES (steady state, week 8+)\n")
expo <- function(nm, col, lab, expect) {
  r <- res[[nm]]$raw %>% filter(time >= 56*24)
  cat(sprintf("  %-34s %10.2f   %s\n", lab, mean(r[[col]]), expect))
}
expo("S06","CNTX","naltrexone Cavg (ng/mL)",        "1-10 (Cmax 5-10)")
expo("S06","CNOL","6-beta-naltrexol Cavg (ng/mL)",  "40-100")
expo("S06","MOROCC","MOR occupancy (%)",            ">= 90 at 24 h post-dose")
expo("S07","CNOL","XR-NTX 6-beta-naltrexol (ng/mL)","2-8 sustained")
expo("S07","MOROCC","XR-NTX MOR occupancy (%)",     ">= 80 through the month")
expo("S08","CACP","acamprosate Cavg (ng/mL)",       "350-500 at 2 g/day")
expo("S10","CDSF","disulfiram+metab Cavg (mg/L)",   "0.5-2")
expo("S12","CTPM","topiramate Cavg (mg/L)",         "5-8 at 300 mg/day")
expo("S13","CGBP","gabapentin Cavg (mg/L)",         "3-6 at 1800 mg/day")
expo("S15","CBAC","baclofen Cavg (ng/mL)",          "200-600 at high dose")
expo("S14","CNMF","nalmefene Cavg (ng/mL)",         "Cmax ~ 16 after 18 mg")
expo("S16","SEMnM","semaglutide Cavg (nM)",         "20-30 at 1.0 mg weekly")

cat("\n[12] DISULFIRAM ARM — ALDH2 blockade and the flush\n")
for (nm in c("S10","S11")) {
  r <- res[[nm]]$raw %>% filter(time >= 56*24)
  cat(sprintf("  %-45s ALDH2 active %.2f  peak ACD %5.1f uM  aversion %.2f  g/day %5.1f\n",
              substr(res[[nm]]$name, 1, 45), mean(r$ALDH2A), max(r$CACD),
              mean(r$AVEXP), mean(r$GPD)))
}

cat("\n", strrep("=", 78), "\n", sep = "")
cat("END OF VALIDATION REPORT\n")
cat(strrep("=", 78), "\n\n", sep = "")

invisible(list(mod = mod, res = res, EP = EP, AUD_INIT = AUD_INIT))
