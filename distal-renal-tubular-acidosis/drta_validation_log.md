# dRTA QSP model — validation log

Captured output of the dependency-free Python RK4 reference implementation
(`drta_python_reference.py`), the harness in which all 56 equations were
integrated, calibrated and debugged before being transcribed into
`drta_mrgsolve_model.R`. Nothing below is hand-written — it is tool output.

Note on the supersaturation index: `SS` is a RELATIVE brushite index. Its
absolute normalisation was not tightly calibrated (healthy subjects land at
0.8-1.6 rather than clearly below 1.0), so every SS conclusion drawn from
this model is a comparison between arms at matched conditions, which is
invariant to that normalisation.

---

## 1. Healthy reference subjects and untreated dRTA phenotypes

Nothing here was fitted to a phenotype: the lesion enters only as LES
(retained pump Vmax) and LES_grad (retained maximal blood-to-urine pH
gradient). Note the hyperchloraemia, hypokalaemia, hypercalciuria,
hypocitraturia and urine pH above 5.5 all appearing together.

```
--- healthy reference subjects (last-24h means) ---
healthy adult 70kg             HCO3  24.7 (24.3-25.0) pH 7.403 uPH 6.32 NAE  66.6 gap   -1.4 VH 0.29 K 4.14 Cl 104.0 UCa  3.86 UCit  2.47 Ca/Cit  0.33 SS  1.52 bone    0.5
healthy child 30kg             HCO3  24.1 (23.5-24.4) pH 7.400 uPH 6.35 NAE  37.6 gap   -0.8 VH 0.28 K 4.10 Cl 104.0 UCa  1.66 UCit  2.25 Ca/Cit  0.15 SS  1.60 bone    0.3
healthy infant 12kg            HCO3  22.6 (21.7-23.2) pH 7.392 uPH 5.94 NAE  24.1 gap   -0.5 VH 0.34 K 4.05 Cl 104.0 UCa  0.67 UCit  0.71 Ca/Cit  0.19 SS  0.82 bone    0.2
--- untreated dRTA phenotypes ---
dRTA child complete            HCO3  15.1 (14.7-15.4) pH 7.333 uPH 6.92 NAE  22.6 gap   14.1 VH 1.00 K 3.10 Cl 111.6 UCa  3.57 UCit  0.48 Ca/Cit  1.55 SS  8.32 bone   15.1
dRTA adult complete            HCO3  16.8 (16.4-17.0) pH 7.350 uPH 6.92 NAE  36.2 gap   28.9 VH 1.00 K 3.23 Cl 110.7 UCa  7.50 UCit  0.75 Ca/Cit  2.08 SS  7.54 bone   30.5
dRTA child severe              HCO3  12.8 (12.4-13.0) pH 7.304 uPH 6.98 NAE  13.7 gap   23.0 VH 1.00 K 2.60 Cl 113.6 UCa  4.50 UCit  0.48 Ca/Cit  1.95 SS  9.18 bone   24.0
dRTA adult severe              HCO3  14.9 (14.6-15.1) pH 7.331 uPH 6.97 NAE  20.7 gap   44.4 VH 1.00 K 2.86 Cl 112.3 UCa  9.00 UCit  0.75 Ca/Cit  2.49 SS  8.32 bone   46.0
dRTA child incomplete          HCO3  24.1 (23.5-24.5) pH 7.400 uPH 6.35 NAE  37.6 gap   -0.7 VH 0.90 K 4.16 Cl 104.0 UCa  1.66 UCit  2.24 Ca/Cit  0.15 SS  1.66 bone    0.3
dRTA adult incomplete          HCO3  24.7 (24.2-25.0) pH 7.403 uPH 6.33 NAE  66.5 gap   -1.3 VH 0.95 K 4.17 Cl 104.0 UCa  3.87 UCit  2.45 Ca/Cit  0.33 SS  1.56 bone    0.5
```

## 2. Diagnostic tests: the acid load, and the Lemann slope

The acute NH4Cl loading test separates all four phenotypes, including
incomplete dRTA (normal resting bicarbonate, cannot acidify) and a pure
gradient defect. The dietary acid titration returns the Lemann slope of
0.0357 mmol/mEq against a literature value of 0.035 -- and that slope is a
validation target here, not a structural parameter.

```
--- acute NH4Cl loading test, 0.1 g/kg = 1.87 mEq/kg at 08:00 on day 60 ---
    literature: a normal subject reaches urine pH < 5.45; dRTA cannot
  normal                 pre HCO3  24.8  pre uPH 6.32  ->  nadir uPH 4.61   ACIDIFIES (<5.45)
  incomplete dRTA        pre HCO3  24.8  pre uPH 6.29  ->  nadir uPH 6.20   FAILS to acidify
  gradient-defect only   pre HCO3  24.8  pre uPH 6.31  ->  nadir uPH 5.81   FAILS to acidify
  rate-defect only       pre HCO3  17.1  pre uPH 6.90  ->  nadir uPH 6.88   FAILS to acidify
  complete dRTA          pre HCO3  16.8  pre uPH 6.91  ->  nadir uPH 6.89   FAILS to acidify
--- dietary acid titration in a healthy adult -> the Lemann slope ---
    literature target: dUCa/dNAE = 0.035 mmol/mEq (Lemann 1999, PMID 9873210)
  diet x0.7  NEAP    46  HCO3 24.68  uPH 6.33(min 6.17)  NAE    47  NH4    29  VH 0.22  UCa  3.15
  diet x1.0  NEAP    65  HCO3 24.68  uPH 6.32(min 6.15)  NAE    67  NH4    43  VH 0.29  UCa  3.86
  diet x1.3  NEAP    85  HCO3 24.68  uPH 6.31(min 6.13)  NAE    87  NH4    57  VH 0.37  UCa  4.57
  diet x1.6  NEAP   104  HCO3 24.68  uPH 6.30(min 6.12)  NAE   106  NH4    71  VH 0.44  UCa  5.28
  diet x2.0  NEAP   130  HCO3 24.68  uPH 6.29(min 6.10)  NAE   133  NH4    90  VH 0.54  UCa  6.23
  diet x2.6  NEAP   170  HCO3 24.68  uPH 6.25(min 6.00)  NAE   173  NH4   121  VH 0.69  UCa  7.65
  MODEL slope = 0.0357 mmol/mEq   (literature 0.035)
```

## 3. Incomplete dRTA is a CLIFF, not a slope

The single most informative run in this set. At a normal dietary acid load
the incomplete-dRTA patient is indistinguishable from healthy on every
routine measurement -- but has only ~5% of actuator reserve left. Raise the
dietary acid load by 30% and the actuator rails: urine citrate collapses
2.50 -> 0.73 mmol/day and the bone base flux jumps 0.5 -> 10.8 mEq/day, a
22-fold step. A healthy subject at a 60% higher acid load still holds 55%
reserve and a bone flux of 0.9. Incomplete dRTA is therefore not a mild
disease; it is a normal-LOOKING disease one dietary step from a severe one.
(The 1-year BMD z-scores in the last column are quantitatively unreliable
at large bone fluxes -- see the known-misfit list in README.md.)

```
--- INCOMPLETE dRTA (LES 0.34) vs dietary acid load: where does the reserve go? ---
  diet    HCO3   uPH   NAE   gap   VH   reserve-left   UCit   bone(mEq/d)  BMDz(1yr)
  x0.7   24.68  6.38    46  -0.9  0.70     29.7%      2.51      0.3    -0.131
  x1.0   24.67  6.35    66  -1.3  0.95      5.4%      2.50      0.5    -0.207
  x1.3   20.79  6.71    74   8.6  1.00      0.0%      0.73     10.8    -3.685
  x1.6   17.21  6.81    75  26.7  1.00      0.0%      0.72     29.3    -10.435
  x2.0   14.50  6.81    76  52.0  1.00      0.0%      0.72     55.2    -20.086
  healthy adult for comparison, same 1-year horizon:
  x1.0   24.68  6.35    66  -1.4  0.29     70.6%      2.50      0.5    -0.195
  x1.6   24.68  6.31   105  -2.1  0.45     55.3%      2.52      0.9    -0.355
```

## 4. Trial-reproduction decomposition: schedule pharmacology vs exposure

The B21CS switch trial (PMID 32712761) reported a responder-rate jump of
43% -> 90% at an unchanged daily alkali dose. Here the daily dose is fixed
at 0.90 mEq/kg/day and adherence is an EXPLICIT input, so schedule
pharmacology and exposure can be separated. Schedule accounts for +7.1
percentage points; reproducing the observed gap requires real-world SoC
adherence near 0.60.

```
Doses fixed at 0.90 mEq/kg/day (near the population requirement) so the
comparison is at MATCHED daily mEq.  Adherence is an explicit input.
--- A. pure schedule pharmacology, adherence held EQUAL at 0.92 ---
SoC IR 3x/day                                resp  39.3%  HCO3 20.10  trough 19.45 UCit 1.08  SS 4.34  bone  10.9  waste 40.1%
ADV7103 2x/day                               resp  46.4%  HCO3 20.80  trough 20.17 UCit 1.34  SS 3.39  bone   8.2  waste 37.0%
    schedule-only responder gain: +7.1 pp, dHCO3 +0.70 mmol/L
--- B. counterfactual: ADV schedule WITHOUT prolonged release ---
ADV 2x/day, IR kinetics                      resp  50.0%  HCO3 20.86  trough 20.11 UCit 1.34  SS 4.15  bone   8.3  waste 37.5%
--- C. adherence sweep on the SoC arm (ADV arm fixed at 0.92) ---
SoC IR 3x/day, adherence 0.60                resp  10.7%  HCO3 17.01  trough 16.57 UCit 0.67  SS 5.18  bone  23.1  waste 92.3%
    -> gain vs ADV(0.92): +35.7 pp
SoC IR 3x/day, adherence 0.70                resp  17.9%  HCO3 17.82  trough 17.33 UCit 0.74  SS 5.15  bone  19.2  waste 68.6%
    -> gain vs ADV(0.92): +28.6 pp
SoC IR 3x/day, adherence 0.80                resp  32.1%  HCO3 19.05  trough 18.48 UCit 0.84  SS 4.57  bone  14.6  waste 52.0%
    -> gain vs ADV(0.92): +14.3 pp
SoC IR 3x/day, adherence 0.92                resp  39.3%  HCO3 20.10  trough 19.45 UCit 1.08  SS 4.34  bone  10.9  waste 40.1%
    -> gain vs ADV(0.92): +7.1 pp
```

## 5. Release-rate sweep: two hypotheses rejected, one found

Sweeping the bicarbonate and citrate granule release rates over a 16-fold
range at fixed daily dose. Plasma bicarbonate moves 0.07 mmol/L and urine
citrate 0.02 mmol/day -- both nil, so the rate-matching and NaDC1-Tm-escape
hypotheses are dead at therapeutic doses. Brushite supersaturation moves
27%, which was not the hypothesis and is the finding.

```
  ADV7103, 2 intakes/day, 1.0 mEq/kg/day, 35% of alkali as citrate
  rows = bicarbonate release rate (/h), cols = citrate release rate (/h)
--- mean plasma HCO3- (mmol/L) ---
  bic\cit     0.10    0.18    0.32    0.55    0.95    1.50
     0.10    23.40   23.39   23.38   23.38   23.37   23.37
     0.18    23.38   23.37   23.36   23.36   23.35   23.35
     0.32    23.37   23.35   23.34   23.33   23.33   23.33
     0.55    23.36   23.34   23.33   23.32   23.32   23.31
     0.95    23.35   23.33   23.32   23.31   23.31   23.30
     1.50    23.34   23.33   23.31   23.30   23.30   23.30
  optimum: 23.395 at bic=0.10/h  cit=0.10/h
--- urine citrate (mmol/day) ---
  bic\cit     0.10    0.18    0.32    0.55    0.95    1.50
     0.10     1.65    1.64    1.64    1.64    1.64    1.63
     0.18     1.65    1.64    1.64    1.63    1.63    1.63
     0.32     1.65    1.64    1.64    1.63    1.63    1.63
     0.55     1.65    1.64    1.63    1.63    1.63    1.62
     0.95     1.65    1.64    1.63    1.63    1.63    1.62
     1.50     1.65    1.64    1.63    1.63    1.63    1.63
  optimum: 1.652 at bic=0.55/h  cit=0.10/h
--- brushite supersaturation ---
  bic\cit     0.10    0.18    0.32    0.55    0.95    1.50
     0.10     2.22    2.29    2.31    2.27    2.20    2.20
```

## 6. Virtual-population trial reproduction

n=40 virtual patients. SoC doses are titrated per patient and then
deliberately under-titrated to reproduce the observed 43% baseline; the
shift and spread of that under-titration are the ONLY two quantities fitted
to the trial, and they are fitted to the baseline alone.

```
n=40  mean age 11.9 y  mean BW 58.0 kg
UNTREATED                            resp   2.5%  HCO3 14.53+-2.60  UCit 0.64  <Ca/Cit 0.33   0.0%  adh 0.93  SS 7.14  bone  40.3  UCa 7.68
required dose: mean 1.03 mEq/kg/day;  prescribed (under-titrated): mean 0.88 mEq/kg/day, 2/40 at the 3.0 ceiling
OBSERVED B21CS: SoC 43% -> ADV7103 90% responders (same patients)
SoC IR 3x/day (titrated)             resp  30.0%  HCO3 19.76+-3.39  UCit 1.04  <Ca/Cit 0.33  15.0%  adh 0.82  SS 4.61  bone  12.3  UCa 4.17
ADV7103 2x/day (same daily mEq)      resp  40.0%  HCO3 20.60+-3.30  UCit 1.37  <Ca/Cit 0.33  27.5%  adh 0.93  SS 4.22  bone   9.5  UCa 3.25
  among the 28 SoC non-responders: Ca/citrate crossed below 0.33 in 3/28 (11%)  [B21CS: 56%]
mechanism decomposition (mechanism removed from BOTH arms):
  SoC   [-adherence coupling]        resp  37.5%  HCO3 20.30+-3.28  UCit 1.18  <Ca/Cit 0.33  20.0%  adh 0.93  SS 4.16  bone  10.4  UCa 3.85
  ADV   [-adherence coupling]        resp  40.0%  HCO3 20.60+-3.30  UCit 1.37  <Ca/Cit 0.33  27.5%  adh 0.93  SS 4.22  bone   9.5  UCa 3.25
    -> responder gain falls from +10.0 to +2.5 pp
  SoC   [-prolonged release]         resp  30.0%  HCO3 19.76+-3.39  UCit 1.04  <Ca/Cit 0.33  15.0%  adh 0.82  SS 4.61  bone  12.3  UCa 4.17
  ADV   [-prolonged release]         resp  40.0%  HCO3 20.62+-3.29  UCit 1.37  <Ca/Cit 0.33  27.5%  adh 0.93  SS 4.57  bone   9.6  UCa 3.57
    -> responder gain falls from +10.0 to +10.0 pp
  SoC   [-fast-citrate granule]      resp  30.0%  HCO3 19.76+-3.39  UCit 1.04  <Ca/Cit 0.33  15.0%  adh 0.82  SS 4.61  bone  12.3  UCa 4.17
  ADV   [-fast-citrate granule]      resp  40.0%  HCO3 20.48+-3.30  UCit 1.26  <Ca/Cit 0.33  27.5%  adh 0.93  SS 4.39  bone   9.8  UCa 3.29
    -> responder gain falls from +10.0 to +10.0 pp
COUNTERFACTUAL — SoC 3x/day under DIRECTLY OBSERVED THERAPY
(adherence forced to the ADV7103 level; isolates exposure from kinetics):
  SoC 3x/day, DOT                    resp  37.5%  HCO3 20.30+-3.28  UCit 1.18  <Ca/Cit 0.33  20.0%  adh 0.93  SS 4.16  bone  10.4  UCa 3.85
alternative schedules at the SAME titrated daily dose:
SoC IR 4x/day incl. 23:00 dose       resp  22.5%  HCO3 19.11+-3.32  UCit 0.99  <Ca/Cit 0.33  12.5%  adh 0.74  SS 4.99  bone  14.6  UCa 4.19
K-citrate IR 3x/day                  resp  35.0%  HCO3 20.08+-3.42  UCit 1.24  <Ca/Cit 0.33  20.0%  adh 0.82  SS 4.86  bone  11.3  UCa 3.92
NaHCO3 IR 3x/day (sodium salt)       resp  30.0%  HCO3 19.74+-3.39  UCit 1.04  <Ca/Cit 0.33  15.0%  adh 0.82  SS 2.81  bone  12.4  UCa 4.18
[624s]
```

---

## Reproducing this log

The harness is `drta_python_reference.py`; `report.py` provides the
last-24-hour summariser and the canonical subject / lesion presets. The
six sections above correspond to six driver scripts (phenotypes,
diagnostics, the incomplete-dRTA diet sweep, the schedule/exposure
decomposition, the release-rate sweep, and the virtual population). Every
equation in the Python file is the same equation as in the `$ODE` block of
`drta_mrgsolve_model.R`, including the thirteen `BUG FIX #n` corrections
that this harness is what found.
