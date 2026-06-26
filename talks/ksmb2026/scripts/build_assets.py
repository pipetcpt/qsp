#!/usr/bin/env python3
"""
build_assets.py — generate every figure the slides reference, using only
Python libs (cairosvg, Pillow, matplotlib, numpy). No R / Graphviz / ImageMagick
required. Outputs to talks/ksmb2026/assets/{figs,plots,shiny,logos}.

Maps are rasterized from the repo's small SVGs (not the multi-MB PNGs).
Simulation plots are *illustrative reduced-model* dynamics (Euler integration of
a simplified subsystem) — clearly labeled as such; the runnable models are the
mrgsolve files themselves.
"""
import os, glob, io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from PIL import Image, ImageDraw, ImageFont
import cairosvg

# ---- palette ----
NAVY="#163A5F"; TEAL="#0E7C86"; GOLD="#C8961E"; INK="#1A2A3A"; MUTED="#6B7A8D"
PANEL="#F5F7FA"; RULE="#D9E1E8"
plt.rcParams.update({
    "font.size":12, "axes.edgecolor":RULE, "axes.labelcolor":INK,
    "text.color":INK, "xtick.color":MUTED, "ytick.color":MUTED,
    "axes.titlecolor":NAVY, "figure.facecolor":"white", "axes.facecolor":"white",
    "axes.grid":True, "grid.color":RULE, "grid.linewidth":0.6, "font.family":"sans-serif",
})

SCRIPT=os.path.dirname(os.path.abspath(__file__))
ROOT=os.path.abspath(os.path.join(SCRIPT,"..","..",".."))
A=os.path.abspath(os.path.join(SCRIPT,"..","assets"))
FIG=os.path.join(A,"figs"); PLOT=os.path.join(A,"plots")
SHINY=os.path.join(A,"shiny"); LOGO=os.path.join(A,"logos")
for d in (FIG,PLOT,SHINY,LOGO): os.makedirs(d,exist_ok=True)

def find_svg(disease):
    g=sorted(glob.glob(os.path.join(ROOT,disease,"*_qsp*.svg")))
    return g[0] if g else None

def svg_to_png(disease,out,width=1600):
    svg=find_svg(disease)
    if not svg:
        print("  ! no SVG for",disease); return False
    cairosvg.svg2png(url=svg,write_to=os.path.join(FIG,out),output_width=width)
    print("  map",out,"<-",os.path.relpath(svg,ROOT))
    return True

# ---------- 1) hero maps ----------
print("[1] hero maps")
svg_to_png("rheumatoid-arthritis","ra_map.png")
svg_to_png("iga-nephropathy","igan_map.png")
svg_to_png("sickle-cell-disease","scd_map.png")
svg_to_png("multiple-myeloma","mm_map.png")

# ---------- 2) montage ----------
print("[2] montage")
MONTAGE=["rheumatoid-arthritis","systemic-lupus-erythematosus","multiple-sclerosis",
 "type2-diabetes","heart-failure-hfref","pulmonary-arterial-hypertension","copd",
 "breast-cancer","multiple-myeloma","chronic-myeloid-leukemia","glioblastoma",
 "sickle-cell-disease","iga-nephropathy","chronic-kidney-disease","parkinsons-disease",
 "alzheimers-disease","epilepsy","schizophrenia","psoriasis","atopic-dermatitis",
 "crohn-disease","ulcerative-colitis","hiv-aids","cystic-fibrosis"]
cols,rows=6,4; cw,ch=440,300; pad=10
sheet=Image.new("RGB",(cols*cw+pad*(cols+1),rows*ch+pad*(rows+1)),"white")
for i,dis in enumerate(MONTAGE):
    svg=find_svg(dis)
    if not svg: continue
    png=cairosvg.svg2png(url=svg,output_width=560)
    im=Image.open(io.BytesIO(png)).convert("RGB")
    im.thumbnail((cw-12,ch-12))
    r,c=divmod(i,cols)
    x=pad+c*(cw+pad)+(cw-im.width)//2; y=pad+r*(ch+pad)+(ch-im.height)//2
    cell=Image.new("RGB",(cw,ch),PANEL)
    cell.paste(im,((cw-im.width)//2,(ch-im.height)//2))
    sheet.paste(cell,(pad+c*(cw+pad),pad+r*(ch+pad)))
sheet.save(os.path.join(FIG,"montage.png")); print("  montage.png",sheet.size)

# ---------- 3) pipeline diagram ----------
print("[3] pipeline")
def box(ax,x,y,w,h,text,fc="white",ec=NAVY,tc=NAVY,fs=12,bold=True):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.02,rounding_size=0.04",
        linewidth=1.6,edgecolor=ec,facecolor=fc))
    ax.text(x+w/2,y+h/2,text,ha="center",va="center",color=tc,fontsize=fs,
            fontweight="bold" if bold else "normal",wrap=True)
def arrow(ax,x1,y1,x2,y2,color=TEAL):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=18,
        lw=2,color=color,shrinkA=2,shrinkB=2))
fig,ax=plt.subplots(figsize=(11,4.4)); ax.set_xlim(0,12); ax.set_ylim(0,5); ax.axis("off")
box(ax,0.2,2.0,2.0,1.0,"Biomedical\nliterature",fc=PANEL,ec=MUTED,tc=INK,fs=11)
box(ax,2.8,1.8,2.6,1.4,"Claude Code\nRoutine\n(LLM agent)",fc=NAVY,ec=NAVY,tc="white",fs=13)
arrow(ax,2.2,2.5,2.8,2.5)
deliv=[("Mechanistic map\n(Graphviz)",4.0),("mrgsolve\nODE model",2.75),
       ("Shiny\ndashboard",1.5),("Curated\nreferences",0.25)]
for t,y in deliv:
    box(ax,6.3,y,2.6,1.05,t,fc="white",ec=TEAL,tc=NAVY,fs=10.5)
    arrow(ax,5.4,2.5,6.3,y+0.52)
box(ax,9.6,2.0,2.2,1.0,"Git\nrepository",fc=PANEL,ec=GOLD,tc=NAVY,fs=12)
for _,y in deliv: arrow(ax,8.9,y+0.52,9.6,2.5,color=GOLD)
fig.tight_layout(); fig.savefig(os.path.join(PLOT,"pipeline.png"),dpi=150,bbox_inches="tight"); plt.close(fig)
print("  pipeline.png")

# ---------- 4) CCR loop ----------
print("[4] ccr loop")
steps=["Select\ndisease","Research\nliterature","Build\nmap","Author\nODEs",
       "Build\nShiny","Ground\nrefs","Commit\n& push"]
fig,ax=plt.subplots(figsize=(6,6)); ax.set_xlim(-1.5,1.5); ax.set_ylim(-1.5,1.5); ax.axis("off"); ax.set_aspect("equal")
n=len(steps); R=1.05
ang=[np.pi/2-2*np.pi*i/n for i in range(n)]
pts=[(R*np.cos(a),R*np.sin(a)) for a in ang]
for i,(x,y) in enumerate(pts):
    ax.add_patch(plt.Circle((x,y),0.30,facecolor=PANEL,edgecolor=TEAL,lw=1.8,zorder=3))
    ax.text(x,y,steps[i],ha="center",va="center",fontsize=8.5,color=NAVY,fontweight="bold",zorder=4)
for i in range(n):
    x1,y1=pts[i]; x2,y2=pts[(i+1)%n]
    a1=np.arctan2(y2-y1,x2-x1)
    sx,sy=x1+0.31*np.cos(a1),y1+0.31*np.sin(a1)
    ex,ey=x2-0.31*np.cos(a1),y2-0.31*np.sin(a1)
    ax.add_patch(FancyArrowPatch((sx,sy),(ex,ey),arrowstyle="-|>",mutation_scale=14,lw=1.8,color=GOLD,zorder=2))
ax.text(0,0.12,"Claude Code\nRoutine",ha="center",va="center",fontsize=13,color=NAVY,fontweight="bold")
ax.text(0,-0.28,"≈ daily",ha="center",va="center",fontsize=10,color=MUTED,style="italic")
fig.savefig(os.path.join(PLOT,"ccr_loop.png"),dpi=150,bbox_inches="tight"); plt.close(fig)
print("  ccr_loop.png")

# ---------- 5) illustrative reduced-model ODE plots ----------
print("[5] reduced-model plots")
def euler(f,x0,T,dt):
    n=int(T/dt); xs=np.zeros((n+1,len(x0))); xs[0]=x0; t=np.linspace(0,T,n+1)
    for k in range(n): xs[k+1]=xs[k]+dt*np.array(f(t[k],xs[k]))
    return t,xs

def hill(C,emax,ec50,h=1.0): return emax*C**h/(ec50**h+C**h)

# IgAN: g=GdIgA1, a=AutoAb, c=IC, u=UPCR, e=eGFR
def igan_run(E):
    def f(t,s):
        g,a,c,u,e=s
        return [1.0*(1-E)-1.0*g, 0.8*g-0.8*a, 0.9*g*a-0.7*c,
                0.6*c-0.5*u, -0.10*c+0.02*(1-e)]
    return euler(f,[1,1,1,1.0,1.0],24,0.05)
t,unt=igan_run(0.0); _,trt=igan_run(0.7)
fig,ax=plt.subplots(1,2,figsize=(7.2,3.0))
ax[0].plot(t,unt[:,3],color=MUTED,lw=2,label="untreated"); ax[0].plot(t,trt[:,3],color=TEAL,lw=2.4,label="triple Rx")
ax[0].set_title("Proteinuria (UPCR)"); ax[0].set_xlabel("months"); ax[0].set_ylabel("rel. UPCR"); ax[0].legend(fontsize=8)
ax[1].plot(t,unt[:,4],color=MUTED,lw=2,label="untreated"); ax[1].plot(t,trt[:,4],color=NAVY,lw=2.4,label="triple Rx")
ax[1].set_title("eGFR (preserved)"); ax[1].set_xlabel("months"); ax[1].set_ylabel("rel. eGFR"); ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(PLOT,"plot_igan.png"),dpi=150); plt.close(fig)

# SCD: HbF frac, Hb, VOC
def scd_run(HU):
    def f(t,s):
        hbf,hb,voc=s
        ind=hill(HU,0.30,1.0,1.5)
        return [ (0.10+ind)-0.5*hbf, 0.4*(hbf-0.0)+0.2*(1.0-hb), -0.6*voc*(1+1.5*hbf)+0.3 ]
    return euler(f,[0.05,0.8,1.0],24,0.05)
t,u0=scd_run(0.0); _,u1=scd_run(3.0)
fig,ax=plt.subplots(1,2,figsize=(7.2,3.0))
ax[0].plot(t,u0[:,1],color=MUTED,lw=2,label="untreated"); ax[0].plot(t,u1[:,1],color=TEAL,lw=2.4,label="hydroxyurea")
ax[0].set_title("Hemoglobin"); ax[0].set_xlabel("months"); ax[0].set_ylabel("rel. Hb"); ax[0].legend(fontsize=8)
ax[1].plot(t,u0[:,2],color=MUTED,lw=2,label="untreated"); ax[1].plot(t,u1[:,2],color=GOLD,lw=2.4,label="hydroxyurea")
ax[1].set_title("Vaso-occlusion rate"); ax[1].set_xlabel("months"); ax[1].set_ylabel("rel. VOC"); ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(PLOT,"plot_scd.png"),dpi=150); plt.close(fig)

# MM: tumor N (logistic + kill), M-protein ~ N
def mm_run(kill):
    def f(t,s):
        N,=s
        return [0.18*N*(1-N/1.0)-kill*N]
    return euler(f,[0.9],24,0.05)
t,n0=mm_run(0.0); _,n1=mm_run(0.45)
fig,ax=plt.subplots(figsize=(4.0,3.0))
ax.plot(t,n0[:,0],color=MUTED,lw=2,label="untreated")
ax.plot(t,np.clip(n1[:,0],1e-4,None),color=NAVY,lw=2.4,label="VRd regimen")
ax.set_yscale("log"); ax.set_title("Tumor burden (M-protein)"); ax.set_xlabel("months"); ax.set_ylabel("rel. burden (log)"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(PLOT,"plot_mm.png"),dpi=150); plt.close(fig)
print("  plot_igan/scd/mm.png")

# ---------- 6) meta-analysis distributions ----------
print("[6] meta distributions")
csv=os.path.join(PLOT,"library_stats.csv")
cl=[]; rf=[]; nd=[]
if os.path.exists(csv):
    import csv as _csv
    with open(csv) as fh:
        for row in _csv.DictReader(fh):
            try:
                cl.append(int(row["clusters"])); rf.append(int(row["refs"])); nd.append(int(row["nodes"]))
            except: pass
if cl:
    fig,ax=plt.subplots(1,3,figsize=(11,3.1))
    for a,data,title,col in [(ax[0],cl,"Pathway clusters / model",TEAL),
                             (ax[1],rf,"PubMed references / model",NAVY),
                             (ax[2],nd,"Graph elements / model",GOLD)]:
        a.hist(data,bins=16,color=col,edgecolor="white",alpha=0.9)
        a.set_title(title); a.set_ylabel("models")
        a.axvline(np.mean(data),color=INK,ls="--",lw=1.2)
        a.text(0.97,0.92,f"mean {np.mean(data):.0f}",transform=a.transAxes,ha="right",va="top",fontsize=9,color=INK)
    fig.tight_layout(); fig.savefig(os.path.join(PLOT,"meta_dist.png"),dpi=150); plt.close(fig)
    print("  meta_dist.png (n=%d)"%len(cl))
else:
    print("  ! no CSV; skip meta_dist")

# ---------- 7) traditional vs LLM ----------
print("[7] traditional vs llm")
fig,ax=plt.subplots(1,2,figsize=(8.2,3.2))
ax[0].bar(["Traditional\n(typical group)","This library"],[8,191],color=[MUTED,TEAL],edgecolor="white")
ax[0].set_title("QSP models available"); ax[0].set_ylabel("# disease models")
for i,v in enumerate([8,191]): ax[0].text(i,v+3,str(v),ha="center",color=NAVY,fontweight="bold")
ax[1].bar(["Traditional","LLM-augmented"],[180*24,6],color=[MUTED,GOLD],edgecolor="white")
ax[1].set_yscale("log"); ax[1].set_title("Effort per model"); ax[1].set_ylabel("person-hours (log)")
ax[1].text(0,180*24*1.1,"months",ha="center",color=NAVY,fontweight="bold")
ax[1].text(1,8,"hours",ha="center",color=NAVY,fontweight="bold")
fig.tight_layout(); fig.savefig(os.path.join(PLOT,"traditional_vs_llm.png"),dpi=150); plt.close(fig)
print("  traditional_vs_llm.png")

# ---------- 8) placeholders (shiny + logos) ----------
print("[8] placeholders")
def placeholder(path,size,title,sub):
    w,h=size; im=Image.new("RGB",(w,h),"white"); d=ImageDraw.Draw(im)
    # dashed border
    for x in range(0,w,24): d.line([(x,2),(x+12,2)],fill=TEAL,width=3); d.line([(x,h-3),(x+12,h-3)],fill=TEAL,width=3)
    for y in range(0,h,24): d.line([(2,y),(2,y+12)],fill=TEAL,width=3); d.line([(w-3,y),(w-3,y+12)],fill=TEAL,width=3)
    try:
        ft=ImageFont.truetype("DejaVuSans-Bold.ttf",40); fs=ImageFont.truetype("DejaVuSans.ttf",24)
    except Exception:
        ft=ImageFont.load_default(); fs=ImageFont.load_default()
    def ctext(y,txt,font,fill):
        bb=d.textbbox((0,0),txt,font=font); d.text(((w-(bb[2]-bb[0]))/2,y),txt,font=font,fill=fill)
    ctext(h*0.34,title,ft,NAVY); ctext(h*0.54,sub,fs,MUTED)
    im.save(path)
placeholder(os.path.join(SHINY,"igan_dashboard.png"),(1000,620),
            "Shiny dashboard","[ screenshot placeholder — drop a capture here ]")
# logos row: three labeled boxes in one image
w,h=1200,320; im=Image.new("RGB",(w,h),"white"); d=ImageDraw.Draw(im)
try: ft=ImageFont.truetype("DejaVuSans-Bold.ttf",26); fs=ImageFont.truetype("DejaVuSans.ttf",18)
except Exception: ft=ImageFont.load_default(); fs=ImageFont.load_default()
labels=["The Catholic University\nof Korea — Pharmacology","PIPET\nPharmacometrics Institute"]
n=len(labels); bw=420; gap=(w-n*bw)//(n+1)
for i,lab in enumerate(labels):
    x=gap+i*(bw+gap)
    d.rounded_rectangle([x,70,x+bw,250],radius=16,outline=TEAL,width=3,fill=PANEL)
    for j,line in enumerate(lab.split("\n")):
        bb=d.textbbox((0,0),line,font=ft); d.text((x+(bw-(bb[2]-bb[0]))/2,120+j*34),line,font=ft,fill=NAVY)
    bb=d.textbbox((0,0),"logo placeholder",font=fs); d.text((x+(bw-(bb[2]-bb[0]))/2,210),"logo placeholder",font=fs,fill=MUTED)
im.save(os.path.join(LOGO,"logos_row.png"))
print("  placeholders done")

# ---------- 9) repo QR code ----------
print("[9] repo QR")
try:
    import qrcode
    qr=qrcode.QRCode(border=1,box_size=10,error_correction=qrcode.constants.ERROR_CORRECT_M)
    qr.add_data("https://github.com/pipetcpt/qsp"); qr.make(fit=True)
    img=qr.make_image(fill_color=NAVY,back_color="white").convert("RGB")
    img.save(os.path.join(FIG,"repo_qr.png")); print("  repo_qr.png",img.size)
except Exception as e:
    print("  ! QR skipped:",e)

print("ALL ASSETS BUILT ->",A)
