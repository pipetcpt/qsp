#!/usr/bin/env python3
# Runs every scenario quoted in README.md / chi_mrgsolve_model.R.
# Usage: python3 chi_verify_run.py
from chi_verify_python import *
import math
DT=0.004
FEED=lambda n=500,q=3.0,amt=845.0:[(i*q,"GGUT",amt,"add") for i in range(n)]
ENT=845.0*8/1440.0
def dzx(d,n=180): return [(i*8.0,"DZXg",d/3.0,"add") for i in range(n)]
def oct_(d,n=260): return [(i*6.0,"OCTs",d/4.0,"add") for i in range(n)]
def loop_(p,hours=48.0,y_init=None,extra=None,target=70.0,tail=12.0):
    q=dict(p); q["loop"]=1.0; q["Gtarget"]=target
    t,y,out=run(q,hours,dt=DT,y_init=y_init,events=FEED()+(extra or []),record=25)
    tl=[o for o in out if o[0]>hours-tail]; m=lambda f:sum(f(o) for o in tl)/len(tl)
    return dict(GIR=m(lambda o:o[2]["GIR"]),G=m(lambda o:o[1][IX["GLU"]]),
                INS=m(lambda o:o[1][IX["INS"]]),BOHB=m(lambda o:o[1][IX["BOHB"]]),
                NH3=m(lambda o:o[1][IX["NH3"]]),GLY=m(lambda o:o[1][IX["GLY"]]),
                CP=m(lambda o:o[1][IX["CPEP"]]),y=y)
def openrun(p,y,h,extra=None,gir=0.0):
    q=dict(p); q["loop"]=0.0; q["GIR_fix"]=gir
    return run(q,h,dt=DT,y_init=y,events=extra or [],record=5)
L=[]
def P_(*a): 
    s=" ".join(str(x) for x in a); print(s); L.append(s)

P_("="*100); P_("CHI QSP MODEL -- VERIFICATION RUN (pure-python RK4 reference, dt=0.004 h)"); P_("="*100)

P_("\n[1] NORMAL NEONATE (3.5 kg, q3h enteral feeds 845 mg/kg systemic) -- calibration targets")
p=P2(); t,yF,out=run(p,60.0,dt=DT,events=FEED(),record=50)
seg=[o for o in out if o[0]>48]; G=[o[1][IX["GLU"]] for o in seg]; I=[o[1][IX["INS"]] for o in seg]
P_(f"    glucose {min(G):.0f}-{max(G):.0f} mg/dL   insulin {min(I):.1f}-{max(I):.1f} uU/mL   "
   f"hepatic glycogen {yF[IX['GLY']]:.0f} mg/kg   ammonia {yF[IX['NH3']]:.0f} umol/L")
t,yfast,outf=openrun(p,yF,24.0)
tgt=[(h,[o for o in outf if o[0]>=h][0]) for h in (4,8,12,18,24) if [o for o in outf if o[0]>=h]]
P_("    FASTING (feeds withheld):")
for h,o in tgt:
    P_(f"      {h:2d} h: G={o[1][IX['GLU']]:5.1f}  INS={o[1][IX['INS']]:5.2f}  glycogen={o[1][IX['GLY']]:5.0f}"
       f"  BOHB={o[1][IX['BOHB']]:5.2f}  FFA={o[1][IX['FFA']]:4.2f}  fuel={o[2]['fuelr']:5.3f}")

P_("\n[2] THE DISEASE AXIS: residual beta-cell K_ATP conductance g")
P_(f"    {'g':>5s} {'Vm(G=60)':>9s} {'IV dextrose':>12s} {'total glc':>10s} {'insulin':>8s} {'BOHB':>6s} {'C-pep':>6s}")
for g in (1.0,0.8,0.6,0.4,0.2,0.1,0.05,0.02):
    r=loop_(P2(g_ab=g,w_ab=1.0))
    d=derived(r["y"],P2(g_ab=g,w_ab=1.0))
    P_(f"    {g:5.2f} {d['Va']:9.1f} {r['GIR']:12.2f} {r['GIR']+ENT:10.2f} {r['INS']:8.1f} {r['BOHB']:6.2f} {r['CP']:6.2f}")

P_("\n[3] DIAZOXIDE vs OCTREOTIDE as a function of g  (effect on IV dextrose requirement, %)")
P_(f"    {'g':>5s} {'GIR0':>7s} {'+DZX 15':>9s} {'red%':>7s} {'+OCT 30':>9s} {'red%':>7s}")
for g in (0.6,0.4,0.3,0.2,0.1,0.05,0.02):
    p0=P2(g_ab=g,w_ab=1.0); b=loop_(p0)
    a1=loop_(p0,72.0,extra=dzx(15)); a2=loop_(p0,72.0,extra=oct_(30))
    f=lambda r: (100*(b["GIR"]-r["GIR"])/b["GIR"]) if b["GIR"]>0.05 else float('nan')
    P_(f"    {g:5.2f} {b['GIR']:7.2f} {a1['GIR']:9.2f} {f(a1):7.1f} {a2['GIR']:9.2f} {f(a2):7.1f}")

P_("\n[4] DIAZOXIDE DOSE-TITRATION BY GENOTYPE (mean glucose, mg/dL; closed loop off)")
def openmean(p,extra=None,hours=54.0,tail=12.0):
    t,y,out=run(p,hours,dt=DT,events=FEED()+(extra or []),record=25)
    tl=[o for o in out if o[0]>hours-tail]; return sum(o[1][IX["GLU"]] for o in tl)/len(tl)
GEN=[("normal",{}),("KATP dominant g=0.60",dict(w_ab=1,g_ab=0.60)),("HNF4A g=0.70",dict(w_ab=1,g_ab=0.70)),
     ("GDH-HI (GLUD1)",dict(w_ab=1,sGTP=0.0)),("GCK activating",dict(w_ab=1,KGshift=0.45)),
     ("KATP recessive g=0.02",dict(w_ab=1,g_ab=0.02))]
P_(f"    {'genotype':24s}"+"".join(f"{d:>8s}" for d in ["0","2.5","5","10","15"]))
for tag,kw in GEN:
    P_(f"    {tag:24s}"+"".join(f"{openmean(P2(**kw),dzx(d) if d>0 else None):8.1f}" for d in (0,2.5,5,10,15)))

P_("\n[5] CEREBRAL FUEL: glucose at which brain supply = demand, vs plasma ketones")
q=P2()
def fuel(G,B,LA=1.3):
    U=q["UbrMax"]*G/(q["Kbr"]+G)*(1+q["acbf"]*max(0,(q["Gcbf"]-G))/q["Gcbf"])
    return U+q["VketMax"]*B/(q["Kket"]+B)+q["VlacMax"]*LA/(q["Klac"]+LA)
P_(f"    {'BOHB mM':>8s} {'threshold G':>12s} {'ketone share':>13s}")
for B in (0.0,0.2,0.5,1.0,2.0,4.0):
    lo,hi=1.0,200.0
    for _ in range(60):
        m=(lo+hi)/2
        if fuel(m,B)<q["CMRreq"]: lo=m
        else: hi=m
    P_(f"    {B:8.2f} {(lo+hi)/2:12.1f} {100*(q['VketMax']*B/(q['Kket']+B))/q['CMRreq']:12.1f}%")

P_("\n[6] GLUCAGON STIMULATION TEST (0.03 mg/kg IV, dextrose held fixed)")
pC=P2(g_ab=0.02,w_ab=1.0); r=loop_(pC,60.0); gir=r["GIR"]
t,y50,_=openrun(pC,r["y"],1.2,gir=gir*0.45)
t,y2,o2=openrun(pC,y50,2.0,extra=[(0.0,"GCGs",30.0,"add")],gir=gir*0.45)
P_(f"    CHI g=0.02   G={y50[IX['GLU']]:5.1f} INS={y50[IX['INS']]:5.0f} BOHB={y50[IX['BOHB']]:4.2f} "
   f"glycogen={y50[IX['GLY']]:5.0f} -> rise {max(o[1][IX['GLU']] for o in o2)-y50[IX['GLU']]:+6.1f} mg/dL")
t,yN,_=run(P2(),60.0,dt=DT,events=FEED(),record=0); t,yNf,_=openrun(P2(),yN,20.0)
t,y3,o3=openrun(P2(),yNf,2.0,extra=[(0.0,"GCGs",30.0,"add")])
P_(f"    normal 20 h fast G={yNf[IX['GLU']]:5.1f} INS={yNf[IX['INS']]:5.1f} BOHB={yNf[IX['BOHB']]:4.2f} "
   f"glycogen={yNf[IX['GLY']]:5.0f} -> rise {max(o[1][IX['GLU']] for o in o3)-yNf[IX['GLU']]:+6.1f} mg/dL")

P_("\n[7] AGENT PANEL, severe recessive diffuse CHI (g=0.02)")
pS=P2(g_ab=0.02,w_ab=1.0); b=loop_(pS)
P_(f"    {'no drug':34s} GIR={b['GIR']:6.2f}  INS={b['INS']:6.1f}")
arms=[("diazoxide 15 mg/kg/d",pS,dzx(15)),("octreotide 30 ug/kg/d",pS,oct_(30)),
      ("glucagon 15 ug/kg/h IV",P2(g_ab=0.02,w_ab=1.0,GCG_inf=15.0),None),
      ("ersodetug 9 mg/kg IV",pS,[(0.0,"ERS",9.0/0.055,"add")]),
      ("sirolimus trough ~16 ng/mL",pS,[(i*24.0,"SIRg",0.06,"add") for i in range(40)]),
      ("nifedipine 0.5 mg/kg q8h",pS,[(i*8.0,"NIF",0.5/1.2,"add") for i in range(180)]),
      ("exendin(9-39) 100 ug/kg/h",P2(g_ab=0.02,w_ab=1.0,EX9inf=100.0),None),
      ("octreotide + diazoxide",pS,oct_(30)+dzx(15)),
      ("octreotide + glucagon",P2(g_ab=0.02,w_ab=1.0,GCG_inf=15.0),oct_(30))]
for n_,p_,e_ in arms:
    r=loop_(p_,extra=e_)
    P_(f"    {n_:34s} GIR={r['GIR']:6.2f}  INS={r['INS']:6.1f}  reduction={100*(b['GIR']-r['GIR'])/b['GIR']:+6.1f}%")

P_("\n[8] RESECTION EXTENT (g stays 0.02 in the remnant)")
P_(f"    {'remnant':>8s} {'GIR':>7s} {'INS':>7s} {'mean G':>8s}  interpretation")
for bm in (1.00,0.50,0.30,0.20,0.10,0.05,0.02):
    r=loop_(P2(g_ab=0.02,w_ab=1.0,BMASS0=bm))
    it = "still dextrose-dependent" if r["GIR"]>0.3 else ("euglycaemic" if r["G"]<120 else "DIABETIC")
    P_(f"    {bm:8.2f} {r['GIR']:7.2f} {r['INS']:7.1f} {r['G']:8.1f}  {it}")

P_("\n[9] LEUCINE LOAD 150 mg/kg PO")
for tag,kw in (("normal",{}),("GDH-HI",dict(w_ab=1,sGTP=0.0)),("KATP g=0.60",dict(w_ab=1,g_ab=0.60))):
    p_=P2(**kw); t,y_,_=run(p_,60.0,dt=DT,events=FEED(),record=0)
    t,y2_,o2_=openrun(p_,y_,4.0,extra=[(0.0,"AA",1.10,"add")])
    P_(f"    {tag:12s} G_pre={y_[IX['GLU']]:6.1f} nadir={min(o[1][IX['GLU']] for o in o2_):6.1f} "
       f"drop={y_[IX['GLU']]-min(o[1][IX['GLU']] for o in o2_):6.1f} peakINS={max(o[1][IX['INS']] for o in o2_):6.1f}")

P_("\n[10] AMMONIA in GDH-HI: diazoxide dependence")
pG=P2(w_ab=1,sGTP=0.0)
P_(f"    GDH-HI            NH3={loop_(pG)['NH3']:6.1f}   GDH-HI+DZX NH3={loop_(pG,72.0,extra=dzx(10))['NH3']:6.1f}"
   f"   normal NH3={loop_(P2())['NH3']:6.1f}  umol/L")

P_("\n[11] FOCAL LESION: burden x density product needed for severe disease")
P_(f"    {'w_ab':>6s} {'density':>8s} {'product':>8s} {'GIR':>7s} {'INS':>7s}")
for w,dn in ((0.02,1),(0.05,3),(0.08,3),(0.10,5),(0.05,10),(0.10,10),(1.0,1)):
    r=loop_(P2(g_ab=0.02,w_ab=w,dens_a=dn))
    P_(f"    {w:6.2f} {dn:8.0f} {w*dn:8.2f} {r['GIR']:7.2f} {r['INS']:7.1f}")
open("chi_verification_output.txt","w").write("\n".join(L)+"\n")
