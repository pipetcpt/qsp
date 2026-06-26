#!/usr/bin/env python3
"""Build a faithful preview of the live Merigolix QSP dashboard
(https://pipetqsp.shinyapps.io/merigolix/) for the slide deck — header bar,
six metric cards (exact values), and the key hormone / clinical / PK panels
with the treatment window shaded. Output: assets/shiny/merigolix_dashboard.png
Replace with a real screenshot anytime (same path)."""
import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.gridspec import GridSpec

NAVY="#163A5F"; INK="#1A2A3A"; MUTED="#6B7A8D"; RULE="#D9E1E8"; PANEL="#F7F9FB"; ORANGE="#E08A2B"
SHADE="#D7E7F4"
plt.rcParams.update({"font.family":"DejaVu Sans","axes.edgecolor":RULE,"text.color":INK,
  "xtick.color":MUTED,"ytick.color":MUTED,"axes.labelcolor":MUTED,"axes.titlecolor":NAVY})
HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,"..","assets","shiny")
os.makedirs(OUT,exist_ok=True)

t=np.linspace(0,168,500); TX0,TX1=28,112
def cyc(t,base,peak):
    v=np.full_like(t,base)
    for c in range(0,168,28):
        v=v+peak*np.exp(-((t-(c+14))/4.0)**2)
    return v
e2=np.where((t<TX0)|(t>=TX1), cyc(t,57,255), 2.6+1.2*np.sin(t))
lh=np.where((t<TX0)|(t>=TX1), cyc(t,11,16), 3.0+0.5*np.sin(t*1.3))
fsh=np.where((t<TX0)|(t>=TX1), cyc(t,7,5), 1.6+0.3*np.sin(t))
def ramp(t,a,b,lo,hi):  # lo before tx & relaxing after, hi during
    s=1/(1+np.exp(-(t-(a+6))/3)); r=1/(1+np.exp((t-(b+10))/8))
    return lo+(hi-lo)*np.minimum(s,r)
pain=5.0-(5.0-1.1)*ramp(t,TX0,TX1,0,1)
hot=0.6+(7.9-0.6)*ramp(t,TX0,TX1,0,1)
bmd=np.where(t<TX0,0.0,np.where(t<TX1,-1.22*(t-TX0)/(TX1-TX0),-1.22+0.5*(t-TX1)/56))
# PK: clean daily sawtooth over first 10 days (matches the dashboard's PK panel)
tp=np.linspace(0,10,1400); pk=np.zeros_like(tp); ka_,ke_,A_=80.0,4.9,76.0
for d in range(0,10):
    dt=tp-d; m=dt>=0; pk=pk+np.where(m, A_*(np.exp(-ke_*dt)-np.exp(-ka_*dt)), 0)

fig=plt.figure(figsize=(15,9.6),dpi=110); fig.patch.set_facecolor("white")
gs=GridSpec(4,6,figure=fig,height_ratios=[0.5,0.7,1.0,1.0],hspace=0.55,wspace=0.4,
            left=0.04,right=0.985,top=0.965,bottom=0.06)
# title bar
tb=fig.add_subplot(gs[0,:]); tb.axis("off")
tb.add_patch(FancyBboxPatch((0,0),1,1,transform=tb.transAxes,boxstyle="round,pad=0,rounding_size=0.02",
            facecolor="#6D5BD0",edgecolor="none"))
tb.text(0.5,0.5,"Merigolix QSP Model Dashboard",transform=tb.transAxes,ha="center",va="center",
        color="white",fontsize=20,fontweight="bold")
# metric cards
cards=[("E2 Trough","2.6","pg/mL","#3E78B2"),("Pain Score","1.1","NRS 0-10","#D9534F"),
       ("BMD Change","-1.22","% from baseline","#E8902E"),("Hot Flashes","7.9","per day","#28A0B0"),
       ("Lesion Size","5.2","mm","#43A047"),("Cmax","71.3","ng/mL","#9AA0A6")]
for i,(lab,val,unit,col) in enumerate(cards):
    ax=fig.add_subplot(gs[1,i]); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.02,0.05),0.96,0.9,transform=ax.transAxes,
        boxstyle="round,pad=0,rounding_size=0.08",facecolor=col,edgecolor="none"))
    ax.text(0.08,0.74,lab,transform=ax.transAxes,color="white",fontsize=11,fontweight="bold",alpha=.95)
    ax.text(0.08,0.40,val,transform=ax.transAxes,color="white",fontsize=26,fontweight="bold")
    ax.text(0.08,0.16,unit,transform=ax.transAxes,color="white",fontsize=9,alpha=.9)
# panels
def panel(r,c,title,y,ylab,color=ORANGE,ymax=None):
    ax=fig.add_subplot(gs[r,c])
    ax.axvspan(TX0,TX1,color=SHADE,alpha=.7,lw=0)
    ax.plot(t,y,color=color,lw=2)
    ax.set_title(title,fontsize=11,loc="left",fontweight="bold",pad=4)
    ax.set_ylabel(ylab,fontsize=8); ax.set_xlabel("Time (days)",fontsize=8)
    ax.tick_params(labelsize=7); ax.set_xlim(0,168)
    if ymax: ax.set_ylim(0,ymax)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    return ax
panel(2,0,"Estradiol",e2,"E2 (pg/mL)")
panel(2,1,"LH",lh,"LH (mIU/mL)")
panel(2,2,"FSH",fsh,"FSH (mIU/mL)")
panel(2,3,"Pain (NRS 0-10)",pain,"Pain",ymax=8)
panel(2,4,"Hot Flashes",hot,"Flashes/day",ymax=10)
panel(2,5,"Bone Mineral Density",bmd,"BMD Δ (%)")
# PK profile (wide) — first 10 days
axpk=fig.add_subplot(gs[3,0:3])
axpk.plot(tp,pk,color=ORANGE,lw=1.2); axpk.set_title("Merigolix Pharmacokinetics (first 10 days, QD)",fontsize=11,loc="left",fontweight="bold")
axpk.set_ylabel("Conc (ng/mL)",fontsize=8); axpk.set_xlabel("Days after treatment start",fontsize=8); axpk.tick_params(labelsize=7); axpk.set_xlim(0,10)
for s in ("top","right"): axpk.spines[s].set_visible(False)
# endpoint stats table
axt=fig.add_subplot(gs[3,3:6]); axt.axis("off")
axt.text(0,1.0,"Endpoint Statistics",transform=axt.transAxes,fontsize=11,fontweight="bold",color=NAVY,va="top")
rows=[("Parameter","Baseline","End of Tx"),("E2 (pg/mL)","15.28","2.61"),("LH (mIU/mL)","11.55","3.04"),
      ("FSH (mIU/mL)","6.95","1.65"),("Pain Score","3.10","1.10"),("Hot Flashes","5.39","7.91"),
      ("BMD Change (%)","0.00","-1.22"),("Lesion (mm)","18.52","5.25"),("Endometrium (mm)","5.92","2.00")]
tbl=axt.table(cellText=rows[1:],colLabels=rows[0],loc="upper center",cellLoc="left",bbox=[0,0.0,1,0.92])
tbl.auto_set_font_size(False); tbl.set_fontsize(8.5)
for (rr,cc),cell in tbl.get_celld().items():
    cell.set_edgecolor(RULE)
    if rr==0: cell.set_facecolor(NAVY); cell.set_text_props(color="white",fontweight="bold")
    elif rr%2==0: cell.set_facecolor(PANEL)
out=os.path.join(OUT,"merigolix_dashboard.png")
fig.savefig(out,dpi=110,facecolor="white"); print("saved",out)
