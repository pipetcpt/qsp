#!/usr/bin/env python3
"""
pk_atropine.py -- fast ocular/systemic PK subsystem for topical atropine,
solved on its own at fine dt because it is (a) LINEAR and (b) autonomous
(it does not read any disease state).  Operator splitting is therefore exact
up to the nonlinearity of the occupancy readouts, which is precisely why the
occupancy summaries below are computed as time-averages of Occ(C(t)) and NOT
as Occ(mean C).

Outputs handed to the slow disease model, per dose strength, at periodic
steady state under one drop at bedtime (t = 0 within each day):

  OCCP_MEAN  mean over 24 h of the POSTERIOR (choroid/sclera) muscarinic
             occupancy -- the growth-control driver
  OCCA_14H   ANTERIOR (iris / ciliary body) receptor-bound fraction at
             +14 h post-dose -- i.e. what a morning/daytime clinic visit
             measures for pupil size and accommodative amplitude
  CPLA_MAX   peak plasma concentration (pg/mL) -- systemic safety readout
  TAU_OFF    apparent washout time constant (days) of the anterior
             receptor-bound state after the last dose
"""

import math

MW = 289.37          # g/mol, atropine free base
DROP_UL = 30.0       # microlitre per drop

# ---- volumes (mL) of a ~9-year-old eye ------------------------------------
V_TEAR, V_COR, V_AQ, V_IRIS, V_VIT, V_CHOR, V_SCL = 7e-3, 0.07, 0.25, 0.30, 4.0, 0.10, 0.35
V_PLA, V_PER = 30.0e3, 40.0e3          # mL (30 L, 40 L) -- 30 kg child

# ---- rate constants (1/day) ----------------------------------------------
K_DRAIN = 500.0      # tear turnover, t1/2 = 2.0 min
K_COR   = 12.0       # tear -> cornea  (=> ~2.3% transcorneal)
K_TSC   = 1.5        # tear -> conjunctiva -> sclera, direct periocular route
K_CA    = 30.0       # cornea -> aqueous
K_AO    = 14.4       # aqueous outflow (2.5 uL/min / 250 uL)
K_AI    = 20.0       # aqueous -> iris / ciliary body
K_IA    = 8.0        # iris -> aqueous
K_AV    = 1.2        # aqueous -> vitreous (poor anterior->posterior transfer)
K_VA    = 0.6        # vitreous -> aqueous
K_VC    = 2.0        # vitreous -> retina / choroid
K_CS    = 6.0        # choroid -> sclera
K_SC    = 3.0        # sclera  -> choroid
K_CP    = 24.0       # choroid -> plasma (well perfused)
K_SL    = 4.0        # sclera  -> systemic / lymphatic
K_EL    = 5.5        # plasma elimination, t1/2 = 3.0 h
K_CPER  = 2.0        # plasma -> peripheral
K_PERC  = 1.5        # peripheral -> plasma
F_SYS   = 0.08       # nasolacrimal + conjunctival systemic bioavailability
                     # calibrated so 1% gives Cmax ~800 pg/mL (measured range
                     # 300-900) and 0.01% gives ~8 pg/mL (below LLOQ)

# ---- binding ------------------------------------------------------------
# Anterior site: apparent EC50 is ~3 log units above the isolated-receptor Ki
# (0.5-2 nM).  The model attributes the gap to iridial MELANIN binding (free
# fraction FU_IRIS) plus the fact that the clinical readout is taken many
# half-lives after the dose.  FU_IRIS is the iris-colour lever.
FU_IRIS  = 0.02      # free fraction in pigmented iris
KD_ANT   = 4.0       # nM, free-drug Kd at the anterior muscarinic site
KOFF_ANT = 0.35      # 1/day -- slowly reversible receptor complex (t1/2 2.0 d)
KON_ANT  = KOFF_ANT / KD_ANT

# Posterior growth-control site: LOW affinity.  Fitted (one constant) to the
# 0.01 / 0.05 / 1% axial-efficacy triplet; see fit_ec50_post() in myp_slow.py.
EC50_POST = 3000.0   # nM


def dose_nmol(pct):
    """nmol of atropine per drop for a w/v percentage strength."""
    mg_per_ml = pct * 10.0
    ug = mg_per_ml * DROP_UL          # mg/mL * uL = ug
    return ug * 1e-3 / MW * 1e9 / 1e3 # ug -> nmol :  ug/1000 = mg; mg/MW = mmol -> *1e6 = nmol


def _dose_nmol_check():
    # 0.01% -> 0.1 mg/mL * 30 uL = 3 ug = 3e-6 g / 289.37 = 1.037e-8 mol = 10.37 nmol
    return 3e-6 / MW * 1e9


def rhs(y, occ_state_only=False):
    (TEAR, COR, AQ, IRIS, VIT, CHOR, SCL, PLA, PER, RB) = y
    d = [0.0] * 10
    d[0] = -(K_DRAIN + K_COR + K_TSC) * TEAR
    d[1] = K_COR * TEAR - K_CA * COR
    d[2] = K_CA * COR - (K_AO + K_AI + K_AV) * AQ + K_IA * IRIS + K_VA * VIT
    d[3] = K_AI * AQ - K_IA * IRIS
    d[4] = K_AV * AQ - (K_VA + K_VC) * VIT
    d[5] = K_VC * VIT + K_SC * SCL - (K_CP + K_CS) * CHOR
    d[6] = K_TSC * TEAR + K_CS * CHOR - (K_SC + K_SL) * SCL
    d[7] = (F_SYS * K_DRAIN * TEAR + K_CP * CHOR + K_SL * SCL
            - (K_EL + K_CPER) * PLA + K_PERC * PER)
    d[8] = K_CPER * PLA - K_PERC * PER
    # slowly reversible anterior receptor complex, driven by FREE iris drug
    c_iris_free = FU_IRIS * IRIS / V_IRIS * 1e3      # nmol/mL -> uM ; *1e3 -> nM
    d[9] = KON_ANT * c_iris_free * (1.0 - RB) - KOFF_ANT * RB
    return d


def simulate(pct, n_days=80, dt=2.0e-4, stop_day=None, adherence=1.0):
    """RK4 with a drop instilled at the start of every dosing day."""
    y = [0.0] * 10
    D = dose_nmol(pct) * adherence
    steps_per_day = int(round(1.0 / dt))
    # accumulators for the LAST full day (periodic steady state)
    occp_sum, n_occp = 0.0, 0
    occa_14h = None
    cpla_max = 0.0
    rb_at_stop = None
    tau_off = None

    for day in range(n_days):
        dosing = (stop_day is None) or (day < stop_day)
        if dosing:
            y[0] += D
        last = (day == n_days - 1)
        for s in range(steps_per_day):
            t_in_day = s * dt
            k1 = rhs(y)
            y2 = [y[i] + 0.5 * dt * k1[i] for i in range(10)]
            k2 = rhs(y2)
            y3 = [y[i] + 0.5 * dt * k2[i] for i in range(10)]
            k3 = rhs(y3)
            y4 = [y[i] + dt * k3[i] for i in range(10)]
            k4 = rhs(y4)
            for i in range(10):
                y[i] += dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
                if y[i] < 0.0:
                    y[i] = 0.0
            # readouts
            c_pla = y[7] / V_PLA * 1e3            # nmol/mL*1e3 = nM... see note
            cpla_max = max(cpla_max, c_pla)
            if last:
                # posterior driver: occupancy of the LOW-affinity site, read
                # from the mean of choroid and sclera free concentration
                c_post = 0.5 * (y[5] / V_CHOR + y[6] / V_SCL) * 1e3   # nM
                occp_sum += c_post / (c_post + EC50_POST)
                n_occp += 1
                if abs(t_in_day - 14.0 / 24.0) < dt:
                    occa_14h = y[9]
    # plasma nM -> pg/mL
    cpla_pg = cpla_max * MW / 1e3 * 1e3
    return dict(occp_mean=occp_sum / max(n_occp, 1),
                occa_14h=occa_14h if occa_14h is not None else y[9],
                occa_end=y[9],
                cpla_pg=cpla_pg,
                y=y)


if __name__ == "__main__":
    print("dose_nmol(0.01%%) = %.3f   (analytic check %.3f)"
          % (dose_nmol(0.01), _dose_nmol_check()))
    print()
    hdr = "%-8s %12s %12s %12s" % ("dose", "OCCP_MEAN", "OCCA_14H", "Cmax_pg/mL")
    print(hdr); print("-" * len(hdr))
    res = {}
    for pct in (0.01, 0.025, 0.05, 0.1, 0.5, 1.0):
        r = simulate(pct, n_days=60)
        res[pct] = r
        print("%-8s %12.5f %12.5f %12.1f"
              % ("%.3f%%" % pct, r["occp_mean"], r["occa_14h"], r["cpla_pg"]))

    print()
    print("apparent Hill slope of each readout over the 0.01%% -> 1%% range")
    for key, name in (("occp_mean", "posterior (efficacy)"),
                      ("occa_14h", "anterior (side effect)")):
        lo, hi = res[0.01][key], res[1.0][key]
        # odds ratio -> Hill slope over a 100-fold dose range
        o_lo, o_hi = lo / (1 - lo), hi / (1 - hi)
        n = math.log(o_hi / o_lo) / math.log(100.0)
        print("   %-24s  %.4f -> %.4f   apparent n = %.3f" % (name, lo, hi, n))
