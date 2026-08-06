#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Hereditary spherocytosis -- the geometry kernel.

Everything in this model funnels through TWO numbers per red cell: its
membrane surface area S (um^2) and its volume V (um^3).  This file contains
the closed-form maps from (S, V) to every geometric read-out the disease is
diagnosed by, and checks each map against a published number.

The central object is the MINIMUM CYLINDRICAL DIAMETER D_c (Canham & Burton
1968): the diameter of the narrowest cylinder a cell of area S and volume V
can be squeezed into WITHOUT stretching its membrane (red cell membrane is
area-incompressible to within ~4%, so this is a hard constraint, not a soft
one).  Model a deformed cell as a right circular cylinder of diameter D and
length L capped by two hemispheres:

    S = pi*D*L
    V = (pi*D^2/4)*(L - D) + pi*D^3/6

Eliminating L = S/(pi*D):

    V = S*D/4 - pi*D^3/12                                            (*)

For fixed S this is non-monotone in D: it rises to a maximum at D = sqrt(S/pi)
(the sphere, L = D) and falls again.  Given the cell's actual V, (*) is a
cubic in D whose SMALLEST positive root is D_c.  Writing (*) as the depressed
cubic D^3 + pD + q = 0 with p = -3S/pi, q = 12V/pi, the trigonometric solution
gives

    cos(theta) = (3q/2p)*sqrt(-3/p) = -6*V*sqrt(pi)/S^(3/2) = -V/V_sph
    D_k        = 2*sqrt(S/pi)*cos(theta/3 - 2*pi*k/3),  k = 0,1,2

and k = 1 is the physical (smallest positive) root.  V_sph = S^(3/2)/(6*sqrt(pi))
is the volume of a sphere with the SAME membrane area -- i.e. the critical
haemolytic volume, the largest volume the cell can hold at all.

Two consequences that the rest of the model is built on:

  (a) The osmotic fragility curve is not an assay property, it is V_sph.  A
      cell lyses when osmotic swelling drives V up to V_sph, and the tonicity
      at which that happens follows from the van 't Hoff behaviour of the
      osmotically active water.  Nothing is fitted.

  (b) Splenic passage is a function of D_c alone.  A cell whose D_c exceeds
      the effective slit calibre is retained in the red pulp cords, and time
      in the cords is itself what strips membrane -- so D_c feeds back on
      itself.  That loop is the disease.
"""
import math

SQRT_PI = math.sqrt(math.pi)

# ----------------------------------------------------------------- geometry


def v_sphere(S):
    """Volume of the sphere with membrane area S -- the critical haemolytic
    volume (um^3).  A cell cannot exceed this without stretching membrane."""
    return S ** 1.5 / (6.0 * SQRT_PI)


def sphericity(S, V):
    """s = V / V_sph in [0, 1].  1 = sphere (no deformability reserve)."""
    return V / v_sphere(S)


def sphericity_index(S, V):
    """The classical sphericity index SI = A_sphere(V) / S, also in [0,1].
    Related to s by SI = s^(2/3)."""
    return (36.0 * math.pi * V * V) ** (1.0 / 3.0) / S


def d_crit(S, V):
    """Minimum cylindrical diameter (um).  Smallest positive root of
    V = S*D/4 - pi*D^3/12."""
    s = sphericity(S, V)
    if s >= 1.0:
        # already a sphere: only one diameter is possible
        return math.sqrt(S / math.pi)
    if s <= 0.0:
        return 0.0
    theta = math.acos(-s)
    return 2.0 * math.sqrt(S / math.pi) * math.cos(theta / 3.0 - 2.0 * math.pi / 3.0)


def d_crit_long(S, V):
    """The OTHER positive root (k = 0): the diameter of the long thin sausage.
    Physically irrelevant for passage but useful as a check that the cubic is
    being solved correctly (the two roots must bracket the sphere diameter)."""
    s = sphericity(S, V)
    theta = math.acos(-max(min(s, 1.0), -1.0))
    return 2.0 * math.sqrt(S / math.pi) * math.cos(theta / 3.0)


def cyl_volume(S, D):
    """Forward map, for round-tripping the cubic."""
    return S * D / 4.0 - math.pi * D ** 3 / 12.0


# ------------------------------------------------------- osmotic behaviour
#
# Cell water obeys van 't Hoff about the isotonic point:
#     V(C) = V_solid + V_water_iso * (C_iso / C)
# V_solid is the osmotically inactive volume: haemoglobin and other solids.
# Its value is not free -- it is fixed by the haemoglobin mass the cell
# carries and the partial specific volume of haemoglobin.

C_ISO = 290.0            # mOsm/kg, isotonic plasma
NACL_PER_MOSM = 0.9 / 308.0   # %NaCl equivalent to 1 mOsm (0.9% = 308 mOsm)


def v_solid(S, V, mchc):
    """Osmotically inactive volume (um^3) from the cell's own MCHC.
    Hb mass = MCHC(g/dL)*V; partial specific volume of Hb 0.75 mL/g; the
    remaining non-Hb solids add ~4% of cell volume."""
    hb_mass_pg = mchc * V * 1e-15 * 10.0 * 1e12   # g/dL * um^3 -> pg
    return hb_mass_pg * 0.75 + 0.04 * V


def osmotic_lysis_point(S, V, mchc):
    """Tonicity (mOsm/kg and %NaCl) at which this cell reaches V_sph and
    lyses.  This IS the osmotic fragility curve; the population curve is the
    distribution of this number over the red cell population."""
    vs = v_solid(S, V, mchc)
    vw = V - vs
    vmax = v_sphere(S)
    if vmax <= vs:
        return C_ISO, C_ISO * NACL_PER_MOSM
    c = C_ISO * vw / (vmax - vs)
    return c, c * NACL_PER_MOSM


def ektacytometry_omin(S, V, mchc):
    """Osmotic-gradient ektacytometry O_min: the hypotonic osmolality of the
    deformability minimum.  It is the point at which swelling has just
    abolished the shape reserve, i.e. where the cell becomes a sphere -- the
    same event as 50% osmotic lysis.  Reported to track surface-to-volume."""
    return osmotic_lysis_point(S, V, mchc)[0]


def ema_mfi_ratio(S, S_ref=140.0):
    """Eosin-5-maleimide binding is to band 3 (and Rh-associated proteins) in
    the membrane.  Per-cell fluorescence is proportional to the number of
    membrane copies, which for a membrane lost as intact vesicles is
    proportional to area.  So the EMA test reads out S, nothing else."""
    return S / S_ref


# ------------------------------------------------------------- self-checks

def _report(name, got, want, tol, unit=""):
    ok = abs(got - want) <= tol
    print("  %-56s %8.3f %s  (published %s%s)  %s"
          % (name, got, unit, want, unit, "OK" if ok else "<<< MISS"))
    return ok


def main():
    print(__doc__.split("\n")[0])
    print()
    print("1. CUBIC ROOT-FINDING IS EXACT (round-trip through eq. *)")
    for (S, V) in [(140.0, 90.0), (127.6, 83.6), (118.0, 88.0), (105.0, 80.0)]:
        D = d_crit(S, V)
        print("   S=%6.1f V=%6.1f  ->  D_c=%6.4f  ->  V(D_c)=%9.6f "
              "(err %.2e)  D_long=%6.3f  D_sph=%6.3f"
              % (S, V, D, cyl_volume(S, D), abs(cyl_volume(S, D) - V),
                 d_crit_long(S, V), math.sqrt(S / math.pi)))

    print()
    print("2. THE NORMAL RED CELL, WITH NOTHING FITTED")
    # Canham & Burton 1968 measured S = 138 um^2, V = 98 um^3 on 50 cells and
    # reported a mean minimum cylindrical diameter of 2.82 um.
    _report("D_c of a normal red cell (S=138, V=98)",
            d_crit(138.0, 98.0), 2.82, 0.10, " um")
    _report("critical haemolytic volume / actual volume (S=140,V=90)",
            v_sphere(140.0) / 90.0, 1.73, 0.10, " x")
    c, nacl = osmotic_lysis_point(140.0, 90.0, 33.0)
    _report("50% osmotic lysis of a normal cell", nacl, 0.43, 0.06, " %NaCl")

    print()
    print("3. THE SENESCENT NORMAL CELL IS NOT A SPHEROCYTE")
    # Waugh et al. Blood 1992;79:1351: young cells S=138.4 V=95.9,
    # old cells S=127.6 V=83.6.  Area falls 7.8%, volume falls 12.8%.
    dy = d_crit(138.4, 95.9)
    do = d_crit(127.6, 83.6)
    print("   young (S=138.4 V=95.9): D_c = %.3f um,  s = %.4f" % (dy, sphericity(138.4, 95.9)))
    print("   old   (S=127.6 V=83.6): D_c = %.3f um,  s = %.4f" % (do, sphericity(127.6, 83.6)))
    print("   -> normal ageing moves D_c by %+.3f um.  Volume is lost FASTER" % (do - dy))
    print("      than area, so an old normal cell is MORE passable, not less.")
    print("      Normal red cells therefore cannot die of a geometric failure;")
    print("      they must be removed by molecular labelling.  HS cells can.")

    print()
    print("4. THE HS SPECTRUM IS AN AREA-DEFICIT SPECTRUM")
    print("   %-22s %7s %7s %7s %8s %9s" %
          ("phenotype (S,V)", "s", "SI", "D_c", "V_sph/V", "OF %NaCl"))
    for lab, S, V, mchc in [("normal    140, 90", 140.0, 90.0, 33.0),
                            ("mild HS   128, 90", 128.0, 90.0, 34.5),
                            ("moderate  118, 88", 118.0, 88.0, 35.5),
                            ("severe    105, 80", 105.0, 80.0, 36.5),
                            ("cord cell  98, 78", 98.0, 78.0, 37.0)]:
        print("   %-22s %7.4f %7.4f %7.3f %8.3f %9.3f"
              % (lab, sphericity(S, V), sphericity_index(S, V), d_crit(S, V),
                 v_sphere(S) / V, osmotic_lysis_point(S, V, mchc)[1]))
    print("   Published HS median corpuscular fragility is 0.55-0.65 %NaCl;")
    print("   the moderate row lands there from geometry alone.")

    print()
    print("5. THE EMA TEST IS A MEMBRANE-AREA ASSAY")
    for lab, S in [("mild HS", 128.0), ("moderate HS", 118.0), ("severe HS", 105.0)]:
        r = ema_mfi_ratio(S)
        print("   %-12s S=%5.1f um^2 -> EMA MFI %.3f of normal (%.1f%% reduction)"
              % (lab, S, r, 100 * (1 - r)))
    print("   The diagnostic cut-off in use is a 16-21%% reduction, which is")
    print("   exactly the area deficit of a mild-to-moderate HS cell.")


if __name__ == "__main__":
    main()
