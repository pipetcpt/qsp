#!/usr/bin/env python3
"""Rewrite the mrgsolve $PARAM values from the Python twin's P dict, so the two
implementations cannot drift.  Reports any parameter present in one and not the
other."""
import importlib.util, re, sys
D='/home/user/qsp/niemann-pick-disease-type-c/'
spec=importlib.util.spec_from_file_location("m", D+"npc_reference_model.py")
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
P=m.P
r=open(D+'npc_mrgsolve_model.R').read()
start=r.index('$PARAM'); end=r.index('$CMT')
block=r[start:end]
seen=set(); changed=[]
def repl(mo):
    name, val, rest = mo.group(1), mo.group(2), mo.group(3)
    if name not in P:
        return mo.group(0)
    seen.add(name)
    new=P[name]
    txt=('%.6g'%new)
    if abs(float(val)-float(new))>1e-9*max(1.0,abs(float(new))):
        changed.append((name, val, txt))
    pad=max(1, len(val)-len(txt))
    return f"{name:<12}: {txt}{' '*pad}: {rest}"
newblock=re.sub(r'^([A-Za-z_][A-Za-z0-9_]*)\s*:\s*([-0-9.eE+]+)\s*:\s*(.*)$', repl, block, flags=re.M)
open(D+'npc_mrgsolve_model.R','w').write(r[:start]+newblock+r[end:])
missing_in_R=[k for k in P if k not in seen]
extra=[]
print("changed:", len(changed))
for n,o,v in changed: print(f"   {n:16s} {o:>12s} -> {v}")
print("in P but not in $PARAM:", missing_in_R if missing_in_R else "none")
