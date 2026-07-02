#!/usr/bin/env python3
"""Build supporting figures for the KBCS 2026 talk (pure Python: cairosvg +
Pillow + matplotlib + numpy + qrcode — no R/Graphviz/ImageMagick needed)."""
import os, glob, io
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
from PIL import Image, ImageDraw, ImageFont
import cairosvg, qrcode

NAVY="#163A5F"; TEAL="#0E7C86"; GOLD="#C8961E"; INK="#1A2A3A"; MUTED="#6B7A8D"
PANEL="#F5F7FA"; RULE="#D9E1E8"
plt.rcParams.update({
    "font.size":12, "axes.edgecolor":RULE, "axes.labelcolor":INK,
    "text.color":INK, "xtick.color":MUTED, "ytick.color":MUTED,
    "axes.titlecolor":NAVY, "figure.facecolor":"white", "axes.facecolor":"white",
    "axes.grid":True, "grid.color":RULE, "grid.linewidth":0.6, "font.family":"sans-serif",
})

HERE=os.path.dirname(os.path.abspath(__file__))
ROOT=os.path.abspath(os.path.join(HERE,"..","..",".."))
A=os.path.abspath(os.path.join(HERE,"..","assets"))
FIG=os.path.join(A,"figs")
os.makedirs(FIG,exist_ok=True)

def find_svg(disease):
    g=sorted(glob.glob(os.path.join(ROOT,disease,"*_qsp*.svg")))
    return g[0] if g else None

# ---------- 1) QSP library montage ----------
print("[1] montage")
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
    cell=Image.new("RGB",(cw,ch),PANEL)
    cell.paste(im,((cw-im.width)//2,(ch-im.height)//2))
    sheet.paste(cell,(pad+c*(cw+pad),pad+r*(ch+pad)))
sheet.save(os.path.join(FIG,"montage.png")); print("  montage.png",sheet.size)

# ---------- 2) translational bridge diagram ----------
print("[2] bridge diagram")
def box(ax,x,y,w,h,text,fc="white",ec=NAVY,tc=NAVY,fs=11,bold=True):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.02,rounding_size=0.04",
        linewidth=1.6,edgecolor=ec,facecolor=fc))
    ax.text(x+w/2,y+h/2,text,ha="center",va="center",color=tc,fontsize=fs,
            fontweight="bold" if bold else "normal")
def arrow(ax,x1,y1,x2,y2,color=TEAL,style="-|>"):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle=style,mutation_scale=18,
        lw=2,color=color,shrinkA=2,shrinkB=2))

fig,ax=plt.subplots(figsize=(12,4.2)); ax.set_xlim(0,13); ax.set_ylim(0,5); ax.axis("off")
box(ax,0.2,1.8,2.1,1.4,"Organ-on-chip\n/ MPS\n(gut·liver·kidney...)",fc=PANEL,ec=MUTED,tc=INK,fs=10.5)
arrow(ax,2.3,2.5,3.0,2.5)
box(ax,3.0,1.8,2.1,1.4,"Chip ODE\nmodel\n(kinetic params)",fc="white",ec=TEAL,tc=NAVY,fs=10.5)
arrow(ax,5.1,2.5,5.8,2.5)
box(ax,5.8,1.8,2.1,1.4,"PBPK\n(whole-body\nscale-up)",fc="white",ec=TEAL,tc=NAVY,fs=10.5)
arrow(ax,7.9,2.5,8.6,2.5)
box(ax,8.6,1.8,2.1,1.4,"QSP\n(disease biology\n+ PD)",fc="white",ec=TEAL,tc=NAVY,fs=10.5)
arrow(ax,10.7,2.5,11.4,2.5,color=GOLD)
box(ax,11.4,1.55,1.5,1.9,"Human\nprediction",fc=NAVY,ec=NAVY,tc="white",fs=11)
ax.text(6.6,4.4,"Human relevance increases left → right",ha="center",color=MUTED,fontsize=10,style="italic")
ax.text(6.6,0.9,"IVIVE  ·  MIDD  ·  dose selection  ·  DDI  ·  special populations",ha="center",color=MUTED,fontsize=9.5)
fig.tight_layout(); fig.savefig(os.path.join(FIG,"bridge.png"),dpi=170,bbox_inches="tight"); plt.close(fig)
print("  bridge.png")

# ---------- 3) future web platform wireframe mock ----------
print("[3] webapp mock")
w,h=1500,900
im=Image.new("RGB",(w,h),"white"); d=ImageDraw.Draw(im)
try:
    ft=ImageFont.truetype("DejaVuSans-Bold.ttf",30); fs=ImageFont.truetype("DejaVuSans.ttf",20)
    fss=ImageFont.truetype("DejaVuSans.ttf",16)
except Exception:
    ft=fs=fss=ImageFont.load_default()
# top bar
d.rectangle([0,0,w,70], fill=NAVY)
d.text((30,18),"QSP Explorer  (concept)",font=ft,fill="white")
# search bar
d.rounded_rectangle([30,95,w-30,150],radius=10,outline=TEAL,width=3,fill=PANEL)
d.text((50,112),"Search a disease, target, or drug...  e.g. 'IL-6', 'sickle cell', 'CD38'",font=fs,fill=MUTED)
# filter chips
chips=["Oncology","Autoimmune","Renal","Cardiovascular","Rare disease","Hepatic"]
x=30
for c in chips:
    tw=d.textlength(c,font=fss)+30
    d.rounded_rectangle([x,165,x+tw,205],radius=18,outline=RULE,width=2,fill="white")
    d.text((x+15,175),c,font=fss,fill=NAVY)
    x+=tw+14
# model cards grid
card_w,card_h=340,190; gap=24
names=[("IgA Nephropathy","20 ODEs · 7 scenarios"),("Sickle Cell Disease","24 ODEs · 7 scenarios"),
       ("Multiple Myeloma","22 ODEs · 6 regimens"),("Rheumatoid Arthritis","20 ODEs · 7 scenarios"),
       ("Heart Failure (HFrEF)","24 ODEs · 5 scenarios"),("NSCLC","18 ODEs · 6 regimens")]
for i,(nm,sub) in enumerate(names):
    r,c=divmod(i,3)
    x0=30+c*(card_w+gap); y0=230+r*(card_h+gap)
    d.rounded_rectangle([x0,y0,x0+card_w,y0+card_h],radius=12,outline=RULE,width=2,fill="white")
    d.rectangle([x0,y0,x0+card_w,y0+8],fill=TEAL)
    d.text((x0+18,y0+26),nm,font=fs,fill=NAVY)
    d.text((x0+18,y0+58),sub,font=fss,fill=MUTED)
    # tiny sparkline mock
    pts=[(x0+18+k*8, y0+140-20*np.sin(k/2.0)-10) for k in range(35)]
    d.line(pts,fill=GOLD,width=3)
    d.rounded_rectangle([x0+18,y0+card_h-40,x0+150,y0+card_h-14],radius=8,outline=TEAL,width=2)
    d.text((x0+30,y0+card_h-36),"Run simulation ▸",font=fss,fill=TEAL)
d.text((30,h-40),"Search → browse mechanistic map & ODEs → run a quick simulation in the browser → export parameters",font=fss,fill=MUTED)
im.save(os.path.join(FIG,"webapp_mock.png"))
print("  webapp_mock.png")

# ---------- 4) repo QR ----------
print("[4] QR")
qr=qrcode.QRCode(border=1,box_size=10,error_correction=qrcode.constants.ERROR_CORRECT_M)
qr.add_data("https://github.com/pipetcpt/qsp"); qr.make(fit=True)
qr.make_image(fill_color=NAVY,back_color="white").convert("RGB").save(os.path.join(FIG,"repo_qr.png"))
print("  repo_qr.png")

print("ALL ASSETS BUILT ->",A)
