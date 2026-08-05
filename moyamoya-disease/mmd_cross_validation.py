#!/usr/bin/env python3
"""
mmd_cross_validation.py -- step 1 of the two-implementation check.

Extracts reference values from the Python/scipy implementation into
/tmp/py_ref.json.  Then run mmd_cross_validation.R, which re-computes the
same quantities in mrgsolve/C++ and prints the paired comparison saved as
mmd_cross_validation.txt.

    python3 mmd_cross_validation.py && Rscript mmd_cross_validation.R
"""
import numpy as np, mmd_reference_model as M, json
o = M.simulate(dict(M.ADULT_I), 3650, n=1500)
h = M.simulate(dict(M.ADULT_H), 3650, n=1500)
sc = dict(M.ADULT_I); sc.update(surg_kind="direct", surg_t=3285.0, G_LEAK=0.35)
b = M.simulate(sc, 3285+90, n=2400)
out={}
for lbl,run,days,keys in (("adult_isch",o,[365,730,1825,3650],["CBFA","CBFF","CBFWS_","OEFA","CVR_INTR","CVR_MEAS","INFA","ISCH","STEN","gS","GMOYA","GPVA","ANEU","HEMH","REMOD"]),
                          ("adult_haem",h,[1825,3650],["CBFA","CBFWS_","GPVA","SIG_PVA","ANEU","HEMH","ISCH","INFA"])):
    for d in days:
        for k in keys:
            kk = "CBF_WS" if k=="CBFWS_" else k
            out[f"{lbl}|{d}|{k}"]=M.at(run,d,kk)
m=b["t"]>=3285
out["bypass|peak_HYPER_REL"]=float(b["HYPER_REL"][m].max())
out["bypass|peak_CBFF"]=float(b["CBFF"][m].max())
out["bypass|d90_CBFWS"]=M.at(b,3285+90,"CBF_WS")
out["bypass|d90_GBYP"]=M.at(b,3285+90,"GBYP")
json.dump(out,open("/tmp/py_ref.json","w"),indent=0)
print(len(out),"values")
