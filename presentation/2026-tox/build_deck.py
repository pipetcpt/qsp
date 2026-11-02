#!/usr/bin/env python3
"""Assemble presentation/2026-tox.pptx on the uploaded KSOT template.

Keeps the template's title slide (slide 0) and builds all content slides on the
template's own layouts, so every slide inherits the master chrome: the QR code
(top-right) and the PIPET / PIPET MPS logos (bottom-left). Content is composed
slide-by-slide in the house palette (navy/teal/gold) matching the figure style.
"""
import copy
from pathlib import Path
from PIL import Image
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from pptx.oxml.ns import qn

ROOT = Path("/home/user/qsp/presentation/2026-tox")
FIG = ROOT / "figures"
TEMPLATE = ROOT.parent / "2026-tox-template.pptx"   # pristine uploaded template
OUT = ROOT.parent / "2026-tox.pptx"                 # built deck

# ---- palette --------------------------------------------------------------
NAVY = RGBColor(0x16, 0x3A, 0x5F)
TEAL = RGBColor(0x1F, 0x8A, 0x8C)
GOLD = RGBColor(0xC8, 0x92, 0x2A)
RED  = RGBColor(0xC0, 0x39, 0x2B)
GREEN= RGBColor(0x2E, 0x8B, 0x6F)
BLUE = RGBColor(0x2B, 0x6C, 0xB0)
INK  = RGBColor(0x1A, 0x22, 0x2C)
GREY = RGBColor(0x5F, 0x6B, 0x76)
LGREY= RGBColor(0xD7, 0xDE, 0xE5)
PALE = RGBColor(0xEE, 0xF3, 0xF8)
PALETEAL = RGBColor(0xE4, 0xF1, 0xEB)
WHITE= RGBColor(0xFF, 0xFF, 0xFF)
BODY = "Calibri"      # safe, ships with Office; true-width in QA

SW, SH = Inches(13.333), Inches(7.5)

prs = Presentation(str(TEMPLATE))

# ---- layouts by name ------------------------------------------------------
LAY = {}
for L in prs.slide_layouts:
    LAY[L.name] = L
TITLE_ONLY = LAY.get("제목만") or LAY.get("Title Only")
SECTION    = LAY.get("구역 머리글") or TITLE_ONLY
BLANK      = LAY.get("빈 화면") or TITLE_ONLY

# ---- delete the template's placeholder content slide (slide index 1) ------
# Proper deletion: drop the relationship (so the part is not re-serialized)
# AND remove the sldId entry, otherwise the orphaned part collides on save.
xml_slides = prs.slides._sldIdLst
slide_ids = list(xml_slides)
if len(slide_ids) > 1:
    sldId = slide_ids[1]
    rId = sldId.get(qn("r:id"))
    prs.part.drop_rel(rId)
    xml_slides.remove(sldId)

# ========================================================================
# helpers
# ========================================================================
def _set_runs(tf, text, size, color, bold=False, font=BODY, align=None, italic=False):
    tf.word_wrap = True
    p = tf.paragraphs[0]
    if align is not None: p.alignment = align
    r = p.add_run(); r.text = text
    f = r.font; f.size = Pt(size); f.bold = bold; f.italic = italic
    f.name = font; f.color.rgb = color
    return p

def add_slide(layout, title=None, title_color=NAVY, title_size=None):
    s = prs.slides.add_slide(layout)
    if title is not None:
        ph = None
        for p in s.placeholders:
            if p.placeholder_format.idx == 0:
                ph = p; break
        if ph is not None:
            ph.text = title
            for para in ph.text_frame.paragraphs:
                for r in para.runs:
                    r.font.color.rgb = title_color
                    if title_size: r.font.size = Pt(title_size)
                    r.font.bold = True
    return s

def textbox(s, text, x, y, w, h, size=16, color=INK, bold=False, font=BODY,
            align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, italic=False, line_spacing=1.04):
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame; tf.word_wrap = True
    tf.vertical_anchor = anchor
    for m in ("margin_left","margin_right","margin_top","margin_bottom"):
        setattr(tf, m, 0)
    first = True
    for line in (text if isinstance(text, list) else [text]):
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        p.alignment = align; p.line_spacing = line_spacing
        if isinstance(line, tuple):
            txt, opts = line
        else:
            txt, opts = line, {}
        r = p.add_run(); r.text = txt
        f = r.font
        f.size = Pt(opts.get("size", size)); f.bold = opts.get("bold", bold)
        f.italic = opts.get("italic", italic); f.name = opts.get("font", font)
        f.color.rgb = opts.get("color", color)
    return tb

def bullets(s, items, x, y, w, h, size=17, color=INK, gap=8, bold_lead=True):
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame; tf.word_wrap = True
    for m in ("margin_left","margin_right","margin_top","margin_bottom"):
        setattr(tf, m, 0)
    first = True
    for it in items:
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        p.space_after = Pt(gap); p.line_spacing = 1.03
        # bullet glyph
        lead, rest = (it if isinstance(it, tuple) else (None, it))
        r0 = p.add_run(); r0.text = "•  "
        r0.font.size = Pt(size); r0.font.color.rgb = TEAL; r0.font.bold = True; r0.font.name = BODY
        if lead:
            rl = p.add_run(); rl.text = lead + "  "
            rl.font.size = Pt(size); rl.font.bold = True; rl.font.color.rgb = NAVY; rl.font.name = BODY
        rr = p.add_run(); rr.text = rest
        rr.font.size = Pt(size); rr.font.color.rgb = color; rr.font.name = BODY
    return tb

def eyebrow(s, text, x=0.73, y=1.35, color=GOLD):
    return textbox(s, text.upper(), x, y, 11.0, 0.35, size=13, color=color, bold=True)

def img_box(path, box_x, box_y, box_w, box_h):
    """return (x,y,w,h) in inches to fit image into box, centered, aspect kept."""
    im = Image.open(path); iw, ih = im.size; ar = iw/ih
    bw, bh = box_w, box_h
    if bw/bh > ar:
        h = bh; w = bh*ar
    else:
        w = bw; h = bw/ar
    x = box_x + (bw-w)/2; y = box_y + (bh-h)/2
    return x, y, w, h

def add_image(s, path, box_x, box_y, box_w, box_h, shadow=False):
    x,y,w,h = img_box(str(path), box_x, box_y, box_w, box_h)
    pic = s.shapes.add_picture(str(path), Inches(x), Inches(y), Inches(w), Inches(h))
    return pic

def rrect(s, x, y, w, h, fill, line=None, line_w=1.0, radius=0.08):
    sp = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    sp.fill.solid(); sp.fill.fore_color.rgb = fill
    if line is None:
        sp.line.fill.background()
    else:
        sp.line.color.rgb = line; sp.line.width = Pt(line_w)
    sp.shadow.inherit = False
    try:
        sp.adjustments[0] = radius
    except Exception:
        pass
    return sp

def source(s, text, y=7.04):
    textbox(s, text, 3.1, y, 8.4, 0.3, size=9.5, color=GREY, align=PP_ALIGN.RIGHT, italic=True)

def stat_row(s, stats, y=2.2, x0=0.9, total_w=11.5, h=2.0, num_size=46, lab_size=14):
    n = len(stats); gap = 0.25
    cw = (total_w - gap*(n-1))/n
    for i,(num, lab, col) in enumerate(stats):
        cx = x0 + i*(cw+gap)
        card = rrect(s, cx, y, cw, h, PALE, radius=0.07)
        textbox(s, num, cx+0.1, y+0.28, cw-0.2, 0.95, size=num_size, color=col, bold=True,
                align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
        textbox(s, lab, cx+0.12, y+h-0.72, cw-0.24, 0.6, size=lab_size, color=INK,
                align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.TOP, line_spacing=1.0)

def table(s, rows, x, y, w, col_w=None, header=True, font_size=13, row_h=0.42,
          header_fill=NAVY, zebra=PALE):
    nr = len(rows); nc = len(rows[0])
    gf = s.shapes.add_table(nr, nc, Inches(x), Inches(y), Inches(w), Inches(row_h*nr))
    tbl = gf.table
    # disable default banding style -> we control fills
    tbl.first_row = header; tbl.horz_banding = False
    if col_w:
        tot = sum(col_w)
        for j,cw in enumerate(col_w):
            tbl.columns[j].width = Inches(w*cw/tot)
    for i,row in enumerate(rows):
        tbl.rows[i].height = Inches(row_h)
        for j,val in enumerate(row):
            c = tbl.cell(i,j)
            c.margin_left = Inches(0.08); c.margin_right = Inches(0.06)
            c.margin_top = Inches(0.02); c.margin_bottom = Inches(0.02)
            c.vertical_anchor = MSO_ANCHOR.MIDDLE
            if header and i==0:
                c.fill.solid(); c.fill.fore_color.rgb = header_fill
                col, bold, sz = WHITE, True, font_size
            else:
                c.fill.solid(); c.fill.fore_color.rgb = WHITE if (i%2==1) else zebra
                col, bold, sz = INK, False, font_size
            tf = c.text_frame; tf.word_wrap = True
            p = tf.paragraphs[0]; p.alignment = PP_ALIGN.LEFT if j==0 else PP_ALIGN.CENTER
            r = p.add_run(); r.text = str(val)
            r.font.size = Pt(sz); r.font.bold = bold; r.font.color.rgb = col; r.font.name = BODY
    return gf

def divider(title, num, subtitle):
    s = add_slide(SECTION, None)
    # big navy band feel via text only (no stripes per design guidance)
    textbox(s, f"PART {num}", 0.9, 2.35, 6, 0.5, size=18, color=GOLD, bold=True)
    textbox(s, title, 0.9, 2.9, 11.2, 1.8, size=40, color=NAVY, bold=True, line_spacing=1.02)
    textbox(s, subtitle, 0.92, 4.75, 10.8, 1.0, size=18, color=GREY, line_spacing=1.1)
    return s

CONTENT_TOP = 1.48   # below title
def big_fig_slide(title, fig, caption=None, eyebrow_txt=None, src=None, figbox=None):
    s = add_slide(TITLE_ONLY, title, title_size=26)
    if eyebrow_txt: eyebrow(s, eyebrow_txt)
    bx = figbox or (0.9, CONTENT_TOP+0.25, 11.55, 5.05)
    add_image(s, fig, *bx)
    if caption:
        textbox(s, caption, 0.95, 6.62, 11.4, 0.5, size=14, color=NAVY, bold=True,
                align=PP_ALIGN.CENTER)
    if src: source(s)
    return s

# ========================================================================
# SLIDE SEQUENCE
# ========================================================================

# ---- Slide 0 is the kept template title slide --------------------------

# ---------------- PART 1: the safety problem & regulatory shift ----------
s = add_slide(TITLE_ONLY, "The problem: 90% of trial drugs fail, and toxicity leads", title_size=25)
eyebrow(s, "Why we are here")
stat_row(s, [("~90%", "candidates entering trials never reach patients", RED),
             ("10–15 yr", "per approved drug", NAVY),
             ("$1–2.6 B", "cost per approved drug", NAVY)], y=2.0, h=1.95, num_size=40)
bullets(s, [
    ("Unexpected toxicity", "a leading cause of late-stage attrition — and of post-marketing withdrawal"),
    ("Healthy-volunteer Phase 1", "characterises exposure in the young and fit, not in the patients who take the drug"),
    ("The costliest failures", "are the safety signals that appear only in real, vulnerable patients"),
], 0.95, 4.3, 11.5, 2.3, size=18)
source(s, "DiMasi J Health Econ 2016 · Kola & Landis Nat Rev Drug Discov 2004")

s = add_slide(TITLE_ONLY, "Conventional safety assessment misses the patients who matter", title_size=23)
eyebrow(s, "The translation gap")
cards = [
    ("2D cell cultures", "No tissue architecture, no organ-level function", GREY),
    ("Animal models", "Persistent interspecies differences in drug disposition and response", GOLD),
    ("Phase 1 trials", "Healthy volunteers only — exclude hepatic & renal impairment, pregnancy, the elderly, rare disease", RED),
]
cw = 3.72; x0 = 0.95; y0 = 2.05; ch = 2.5
for i,(h,b,c) in enumerate(cards):
    cx = x0 + i*(cw+0.2)
    rrect(s, cx, y0, cw, ch, PALE, radius=0.06)
    textbox(s, h, cx+0.25, y0+0.25, cw-0.5, 0.7, size=18, color=c, bold=True)
    textbox(s, b, cx+0.25, y0+1.05, cw-0.5, ch-1.2, size=14.5, color=INK, line_spacing=1.08)
textbox(s, "The result: in-vitro signals and animal data rarely translate into human exposure and response for the people actually at risk.",
        0.95, 4.95, 11.5, 1.0, size=18, color=NAVY, bold=True, line_spacing=1.12)
source(s, "Van Norman, JACC Basic Transl Sci 2019")

s = add_slide(TITLE_ONLY, "The regulatory shift: from an animal mandate to NAMs", title_size=24)
eyebrow(s, "What changed")
bullets(s, [
    ("1938 → 2022", "U.S. law required animal testing before human use for 80+ years"),
    ("FDA Modernization Act 2.0 (2022)", "replaced “animal tests” with “nonclinical tests” — explicitly opening the door to in vitro / in silico methods"),
    ("FDA Roadmap to Reducing Animal Testing (2025)", "a phased rollout beginning with monoclonal antibodies; goal: animal testing becomes the exception"),
    ("New approach methodologies (NAMs)", "microphysiological systems (MPS / organ-chips), organoids, and in silico models"),
], 0.95, 2.0, 11.6, 3.2, size=18, gap=11)
textbox(s, "But an in-vitro signal becomes decision-grade evidence only when a mechanistic model translates it into human exposure and response.",
        0.95, 5.5, 11.5, 1.0, size=18, color=NAVY, bold=True, line_spacing=1.12)
source(s, "FDA Modernization Act 2022, Pub. L. 117-328 · FDA Roadmap 2025")

big_fig_slide("NAMs need a translational bridge — and clinical pharmacology is it",
              FIG/"ext/nam_loop.png",
              caption="PBPK and QSP are the qualification & translation step that pattern-matching alone cannot supply",
              eyebrow_txt="The NAM translation loop",
              figbox=(1.3, 1.95, 10.7, 4.3))
source(prs.slides[-1], "Cao Y, Clin Pharmacol Ther 2025")

# ---------------- PART 2: what QSP / PBPK is -----------------------------
divider("QSP & PBPK — the mechanistic bridge", "2",
        "Linking drug exposure to target, pathway, organ, and patient outcomes — mathematically")

s = add_slide(TITLE_ONLY, "QSP and PBPK describe the causal chain, mathematically", title_size=23)
eyebrow(s, "What it is")
rrect(s, 1.4, 2.0, 10.5, 0.95, NAVY, radius=0.12)
textbox(s, "drug  →  target  →  pathway  →  disease  →  patient", 1.5, 2.0, 10.3, 0.95,
        size=26, color=WHITE, bold=True, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
twocol = [
    ("PBPK", "Where does the drug go, and at what concentration, in a whole human body?", BLUE),
    ("QSP", "Given that concentration, what happens to the disease — and to safety?", TEAL),
]
for i,(h,b,c) in enumerate(twocol):
    cx = 1.4 + i*5.35
    rrect(s, cx, 3.35, 5.1, 1.7, PALE, radius=0.07)
    textbox(s, h, cx+0.3, 3.55, 4.5, 0.6, size=22, color=c, bold=True)
    textbox(s, b, cx+0.3, 4.2, 4.5, 0.8, size=15.5, color=INK, line_spacing=1.1)
textbox(s, "Together: first-in-human dose, safety margins, DDI risk, special-population response, and IVIVE of chip-derived parameters.",
        1.4, 5.35, 10.5, 0.9, size=17, color=NAVY, bold=True, line_spacing=1.12)

s = add_slide(TITLE_ONLY, "Anatomy of a QSP model: a coupled nonlinear system", title_size=23)
eyebrow(s, "Under the hood")
rrect(s, 2.6, 1.95, 8.1, 0.95, PALETEAL, radius=0.12)
textbox(s, "dx/dt  =  f ( x(t),  θ,  u(t) ),      x(0) = x₀", 2.6, 1.95, 8.1, 0.95,
        size=24, color=NAVY, bold=True, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE, font="Cambria")
rows = [
    ("States  x", "concentrations · receptor occupancy · signalling species · cell populations · biomarkers · endpoints"),
    ("Parameters  θ", "rate constants · EC₅₀ / Emax · Hill n · clearances · binding affinities"),
    ("Inputs  u(t)", "the dosing regimen — the one controllable handle of therapy"),
]
y = 3.25
for lead, body in rows:
    textbox(s, lead, 1.0, y, 3.2, 0.6, size=18, color=TEAL, bold=True)
    textbox(s, body, 4.2, y, 8.1, 0.7, size=16, color=INK, line_spacing=1.05)
    y += 0.95
textbox(s, "Our library's models: 14–82 coupled ODEs and 70–150+ parameters per disease.",
        1.0, 6.1, 11.3, 0.6, size=17, color=NAVY, bold=True)

s = add_slide(TITLE_ONLY, "The same nonlinear motifs recur across every disease", title_size=23)
eyebrow(s, "A shared grammar")
table(s, [
    ["Motif", "Governing idea", "Typical use"],
    ["Hill / Emax", "E(C) = Emax·Cⁿ / (EC₅₀ⁿ + Cⁿ)", "any drug dose–response"],
    ["Indirect response", "dR/dt = k_in·f(·) − k_out·R", "biomarker turnover"],
    ["Mass-action / TMDD", "drug–target binding, target-mediated clearance", "mAbs, saturating drugs"],
    ["Logistic growth", "dN/dt = rN(1−N/K) − kill(C)·N", "tumour burden under therapy"],
    ["Power-law feedback", "e.g. EPO ∝ (Hb₀/Hb)^γ", "homeostatic loops"],
], 0.95, 2.0, 11.5, col_w=[2.0,3.6,2.4], row_h=0.62, font_size=14.5)
textbox(s, "One grammar across 416 models — and across every chip-to-human bridge.",
        0.95, 6.1, 11.4, 0.5, size=16, color=NAVY, bold=True)

s = add_slide(TITLE_ONLY, "Why it matters: Model-Informed Drug Development", title_size=23)
eyebrow(s, "MIDD — endorsed by FDA / EMA")
table(s, [
    ["Stage", "What QSP/PBPK contributes", "Safety payoff"],
    ["Target validation", "does modulation reach the phenotype?", "kills implausible targets early"],
    ["First-in-human dose", "exposure–response & safety margins from mechanism", "safer, leaner Phase 1"],
    ["Trial design", "virtual patients → enrichment, endpoints", "fewer patients exposed to risk"],
    ["Special populations", "hepatic/renal impairment, age, pregnancy", "dose the excluded patient"],
    ["Combinations", "synergy / antagonism in silico", "rational, safer polypharmacy"],
], 0.95, 2.0, 11.5, col_w=[2.3,4.2,3.0], row_h=0.66, font_size=14)
textbox(s, "QSP turns efficacy–safety trade-offs into tunable optimization problems — the thread through today's case studies.",
        0.95, 6.35, 11.4, 0.6, size=16, color=NAVY, bold=True, line_spacing=1.1)

# ---------------- PART 3: the QSP library (engine) ----------------------
divider("A living QSP model library", "3",
        "An LLM-augmented workflow that builds one literature-grounded, quality-gated disease model per day")

s = add_slide(TITLE_ONLY, "The bottleneck: one QSP model takes expert-months to years", title_size=22)
eyebrow(s, "Why there were so few")
bullets(s, [
    ("Deep synthesis", "across pathophysiology and pharmacology for a single disease"),
    ("Careful formulation", "ODE structure, parameterization, implementation, calibration"),
    ("The cost", "expert-months to years per disease — so only a handful ever get a model"),
    ("The consequence", "models are siloed, rarely open, rarely reproducible; the long tail is never modelled"),
], 0.95, 2.1, 11.5, 2.8, size=18, gap=12)
rrect(s, 0.95, 5.2, 11.5, 1.0, PALE, radius=0.08)
textbox(s, "One disease ≈ one expert-year.    Scarce · slow · closed.", 1.1, 5.2, 11.2, 1.0,
        size=22, color=NAVY, bold=True, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)

s = add_slide(TITLE_ONLY, "Thesis: an autonomous LLM agent builds the models", title_size=23)
eyebrow(s, "The engine")
textbox(s, "An autonomous LLM coding agent builds complete, quality-gated, referenced QSP models — one disease per day — every artifact under version control.",
        0.95, 1.95, 11.5, 1.0, size=18, color=INK, line_spacing=1.15)
arts = [("🗺️  Mechanistic map", "Graphviz network ≥100 nodes", TEAL),
        ("⚙️  mrgsolve ODE model", "PK ⊕ disease PD, ≥15 compartments", NAVY),
        ("📊  Shiny dashboard", "≥6 interactive tabs", GOLD),
        ("📚  Curated references", "≥30 PubMed, calibration memos", GREEN)]
for i,(h,b,c) in enumerate(arts):
    cx = 0.95 + (i%2)*5.85; cy = 3.15 + (i//2)*1.55
    rrect(s, cx, cy, 5.6, 1.35, PALE, radius=0.07)
    textbox(s, h, cx+0.3, cy+0.2, 5.1, 0.6, size=18, color=c, bold=True)
    textbox(s, b, cx+0.3, cy+0.78, 5.1, 0.5, size=14.5, color=INK)
source(s, "github.com/pipetcpt/qsp — fully open")

big_fig_slide("The Claude Code Routine: an autonomous daily loop",
              FIG/"ext/architecture.png",
              caption="Select → research (~50 PubMed) → map → ODE model → dashboard → ground → commit, under a git stop-hook",
              eyebrow_txt="One disease, end-to-end, every run",
              figbox=(1.5, 2.0, 10.3, 4.1))

s = add_slide(TITLE_ONLY, "Guardrails: why you can trust code an LLM wrote", title_size=23)
eyebrow(s, "Structure + grounding + verification")
cards = [
    ("Quality gates", "map ≥100 nodes / ≥8 clusters · ODE ≥15 compartments / ≥5 scenarios · ≥30 references", TEAL),
    ("Structural templates", "reusable motifs (Hill, turnover, TMDD) shrink the failure space", NAVY),
    ("Literature grounding", "every parameter set cites a calibration trial — the main anti-hallucination mechanism", GOLD),
    ("Reproducibility", "diffable, auditable, re-runnable code in git; human review before use", GREEN),
]
for i,(h,b,c) in enumerate(cards):
    cx = 0.95 + (i%2)*5.85; cy = 2.1 + (i//2)*2.0
    rrect(s, cx, cy, 5.6, 1.8, PALE, radius=0.07)
    textbox(s, h, cx+0.3, cy+0.22, 5.0, 0.6, size=18, color=c, bold=True)
    textbox(s, b, cx+0.3, cy+0.85, 5.0, 0.9, size=14.5, color=INK, line_spacing=1.1)

# S1 figures
big_fig_slide("416 models and counting — a library that compounds",
              FIG/"s1_growth.png", eyebrow_txt="Scale",
              caption="Documented milestone trajectory: KSMB 191 → KBCS 258 → KSOT 416 models",
              figbox=(1.8, 2.0, 9.7, 4.2))
big_fig_slide("22 therapeutic areas — three-quarters carry a safety axis",
              FIG/"s1_areas.png", eyebrow_txt="Coverage · the safety atlas",
              caption="74% of 416 models encode an explicit toxicity / adverse-event axis",
              figbox=(1.3, 1.95, 10.7, 4.5))
big_fig_slide("Every model clears hard quality gates",
              FIG/"s1_quality.png", eyebrow_txt="Size & grounding",
              figbox=(1.3, 2.1, 10.7, 4.1))
big_fig_slide("416 models, one consistent schema — the atlas",
              FIG/"s1_atlas.png", eyebrow_txt="One framework",
              caption="Each point is a model: ODE size × map size, marker area ∝ references",
              figbox=(1.7, 2.0, 9.9, 4.3))

# ---------------- PART 4: five case studies -----------------------------
divider("Five safety case studies", "4",
        "Autoimmune · oncology · obesity · aging — each an efficacy–safety trade-off turned into optimization")

# --- Study overview table
s = add_slide(TITLE_ONLY, "Five studies, one pattern", title_size=25)
eyebrow(s, "Roadmap")
table(s, [
    ["Study", "Disease area", "The tunable handle", "Safety axis"],
    ["① Multiple sclerosis", "autoimmune / neuro", "natalizumab dosing interval", "PML risk"],
    ["② Checkpoint inhibitors", "oncology (melanoma)", "nivolumab : ipilimumab ratio", "grade 3–4 irAE"],
    ["③ Incretin therapy", "obesity + aging", "target dose & escalation", "GI AE · lean-mass loss"],
    ["④ DOAC in the elderly", "aging / organ function", "age- & renal-based dose", "bleeding"],
    ["★ Merigolix", "gynaecologic", "GnRH-antagonist dose", "bone loss · hot flashes"],
], 0.9, 2.0, 11.6, col_w=[2.6,2.5,3.4,2.4], row_h=0.66, font_size=13.5)
textbox(s, "Plus Study ⓪ — the library itself as a safety atlas across all 416 models.",
        0.95, 6.35, 11.4, 0.5, size=15.5, color=NAVY, bold=True)

# ===== Study 1 (MS) =====
s = add_slide(TITLE_ONLY, "Study ① — Multiple sclerosis: the sharpest trade-off in neurology", title_size=21)
eyebrow(s, "Autoimmune / neurologic")
bullets(s, [
    ("Natalizumab (anti-α4-integrin)", "among the most effective relapse suppressors — but it lifts CNS immune surveillance"),
    ("The risk", "progressive multifocal leukoencephalopathy (PML) from JC-virus reactivation"),
    ("The question", "can the dosing interval be tuned to keep efficacy while cutting PML risk?"),
], 0.95, 2.0, 11.5, 2.0, size=18, gap=11)
rrect(s, 0.95, 4.4, 11.5, 1.75, PALETEAL, radius=0.07)
textbox(s, "QSP layer:  α4 receptor saturation over the dosing interval  →  relapse-rate suppression  +  a PML-hazard multiplier on a stratified base risk (anti-JCV index · prior immunosuppressant · duration).",
        1.2, 4.6, 11.0, 1.4, size=16.5, color=NAVY, line_spacing=1.15, anchor=MSO_ANCHOR.MIDDLE)
source(s, "Polman NEJM 2006 · Ho Lancet Neurol 2017 · Foley Lancet Neurol 2022")

big_fig_slide("The dosing interval is the controllable handle",
              FIG/"s2_saturation.png", eyebrow_txt="Study ① · MS",
              caption="q6w extended-interval dosing retains ~62% relapse reduction vs 68% at q4w",
              figbox=(1.1, 2.0, 11.1, 4.3))
big_fig_slide("PML risk is stratified — and steep",
              FIG/"s2_pml_table.png", eyebrow_txt="Study ① · MS",
              caption="Risk rises with anti-JCV index, prior immunosuppression, and time on drug",
              figbox=(1.9, 2.0, 9.5, 4.3))
big_fig_slide("Extended-interval dosing: −76% PML risk for 0.3 relapse / 6 yr",
              FIG/"s2_benefit_risk.png", eyebrow_txt="Study ① · MS · the optimization",
              figbox=(1.5, 2.0, 10.3, 4.4))
big_fig_slide("Where natalizumab-q6w sits among MS therapies",
              FIG/"s2_dmt_map.png", eyebrow_txt="Study ① · MS",
              caption="EID moves natalizumab toward the efficacy of ocrelizumab at a fraction of the risk index",
              figbox=(1.6, 2.0, 10.1, 4.3))

# ===== Study 2 (ICI) =====
s = add_slide(TITLE_ONLY, "Study ② — Checkpoint inhibitors: two antibodies, one frontier", title_size=21)
eyebrow(s, "Oncology — metastatic melanoma")
bullets(s, [
    ("Nivolumab + ipilimumab", "the most effective and most toxic melanoma regimen"),
    ("The key biology", "ipilimumab's anti-tumour priming saturates at low dose, but its toxicity keeps climbing with dose"),
    ("The question", "why does lowering ipilimumab (CheckMate 511) preserve efficacy while halving toxicity?"),
], 0.95, 2.0, 11.5, 2.0, size=18, gap=11)
rrect(s, 0.95, 4.4, 11.5, 1.75, PALETEAL, radius=0.07)
textbox(s, "QSP layer:  PD-1 & CTLA-4 occupancy  →  anti-tumour and self-reactive effector pools (Bliss synergy)  →  ORR and grade 3–4 irAE, read through two different CTLA-4 occupancy curves.",
        1.2, 4.6, 11.0, 1.4, size=16.5, color=NAVY, line_spacing=1.15, anchor=MSO_ANCHOR.MIDDLE)
source(s, "Larkin NEJM 2015/2019 (CM067) · Lebbé JCO 2019 (CM511) · Weber JCO 2017")

big_fig_slide("Dosing checkpoint inhibitors: the efficacy–safety frontier",
              FIG/"s3_pareto.png", eyebrow_txt="Study ② · oncology · the optimization",
              caption="Shifting ipilimumab down drops grade 3–4 irAE ~12 points at held efficacy",
              figbox=(1.6, 2.0, 10.1, 4.3))
s = big_fig_slide("Response and toxicity run on two different clocks",
              FIG/"s3_timecourse.png", eyebrow_txt="Study ② · oncology",
              figbox=(1.2, 2.0, 7.0, 4.4))
add_image(s, FIG/"s3_organ_timing.png", 7.9, 2.2, 5.3, 3.9)
textbox(s, "irAEs arrive organ by organ — a kinetic signature for safety monitoring",
        0.95, 6.5, 11.4, 0.4, size=13.5, color=NAVY, bold=True, align=PP_ALIGN.CENTER)

# ===== Study 3 (obesity) =====
s = add_slide(TITLE_ONLY, "Study ③ — Incretin therapy: a three-way optimization", title_size=22)
eyebrow(s, "Obesity / metabolic — with an aging subgroup")
bullets(s, [
    ("Semaglutide & tirzepatide", "the levers are target dose and escalation schedule"),
    ("Three competing outcomes", "weight loss (benefit) vs GI adverse events (tolerability) vs lean-mass loss (safety)"),
    ("The aging angle", "lean-mass loss is a sarcopenia/frailty risk in older adults"),
], 0.95, 2.0, 11.5, 2.0, size=18, gap=11)
rrect(s, 0.95, 4.4, 11.5, 1.75, PALETEAL, radius=0.07)
textbox(s, "QSP layer:  exposure → appetite suppression → energy-balance ODE; weight change partitioned into fat and fat-free mass (age-sensitive Forbes rule); GI AE = nausea driven by escalation, blunted by tolerance.",
        1.2, 4.6, 11.0, 1.4, size=16.5, color=NAVY, line_spacing=1.15, anchor=MSO_ANCHOR.MIDDLE)
source(s, "Wilding NEJM 2021 (STEP 1) · Jastreboff NEJM 2022 (SURMOUNT-1)")

big_fig_slide("Titration is the handle: faster escalation buys only nausea",
              FIG/"s4_titration.png", eyebrow_txt="Study ③ · obesity",
              figbox=(1.2, 2.0, 11.0, 4.3))
big_fig_slide("More weight loss costs more GI burden — the frontier bends",
              FIG/"s4_frontier.png", eyebrow_txt="Study ③ · obesity · the optimization",
              caption="Tirzepatide dominates — more weight loss per unit of GI burden",
              figbox=(1.9, 2.0, 9.5, 4.3))
s = big_fig_slide("The aging safety axis: lean-mass loss",
              FIG/"s4_composition.png", eyebrow_txt="Study ③ · obesity + aging",
              figbox=(1.1, 2.0, 6.6, 4.4))
add_image(s, FIG/"s4_aging.png", 7.7, 2.1, 5.5, 4.0)
textbox(s, "Lean fraction of weight lost rises 36% → 47% from age 45 to 70 — a dose-reduction signal",
        0.95, 6.55, 11.4, 0.4, size=13, color=NAVY, bold=True, align=PP_ALIGN.CENTER)

# ===== Study 4 (DOAC) =====
s = add_slide(TITLE_ONLY, "Study ④ — DOACs in the elderly: dosing the excluded patient", title_size=21)
eyebrow(s, "Aging / renal–hepatic impairment")
bullets(s, [
    ("The gap", "Phase 1 characterises DOAC exposure in young, normal-renal volunteers"),
    ("The reality", "the patients who take DOACs are old and renally impaired — exactly those excluded"),
    ("The question", "can PBPK predict exposure and individualize dosing for patients never studied?"),
], 0.95, 2.0, 11.5, 2.0, size=18, gap=11)
rrect(s, 0.95, 4.4, 11.5, 1.75, PALETEAL, radius=0.07)
textbox(s, "PBPK/PD layer:  age → creatinine clearance → DOAC exposure (renal-fraction scaling) → bleeding & stroke risk; fixed vs individualized dosing in a 4,000-patient virtual elderly cohort.",
        1.2, 4.6, 11.0, 1.4, size=16.5, color=NAVY, line_spacing=1.15, anchor=MSO_ANCHOR.MIDDLE)
source(s, "US product labels (dabigatran, apixaban, rivaroxaban, edoxaban) · Rowland–Tozer")

big_fig_slide("Age is a renal-function covariate — exposure rises silently",
              FIG/"s5_age_crcl.png", eyebrow_txt="Study ④ · aging",
              figbox=(1.5, 2.0, 10.3, 4.4))
big_fig_slide("Validating exposure — and finding where mechanism is needed",
              FIG/"s5_validation.png", eyebrow_txt="Study ④ · aging",
              caption="Renal scaling GMFE 1.18×; dabigatran under-predicted → P-gp term → 1.14× (mechanism, not curve-fitting)",
              figbox=(3.2, 1.78, 6.9, 4.85))
s = big_fig_slide("Two risks oppose; individualized dosing dominates",
              FIG/"s5_exposure_response.png", eyebrow_txt="Study ④ · aging · the optimization",
              figbox=(1.1, 2.0, 6.7, 4.3))
add_image(s, FIG/"s5_individualization.png", 7.7, 2.3, 5.5, 3.6)
textbox(s, "Model-informed dose reduction: over-exposure 68% → 24%, modelled bleeding 6.3 → 4.1 / 100 PY",
        0.95, 6.5, 11.4, 0.4, size=12.5, color=NAVY, bold=True, align=PP_ALIGN.CENTER)

# ---------------- PART 5: flagship + NAMs bridge ------------------------
divider("From optimization to the clinic — and to NAMs", "5",
        "The flagship Merigolix case, then the gut–liver–kidney MPS + PBPK bridge")

s = add_slide(TITLE_ONLY, "Flagship — Merigolix: the estrogen-threshold window", title_size=22)
eyebrow(s, "Gynaecologic · GnRH antagonist for endometriosis")
bullets(s, [
    ("Mechanism", "oral GnRH-receptor antagonist → ↓LH/FSH → ↓estradiol → lesion & pain regression"),
    ("The trade-off", "too little E2 suppression → no relief; too much → bone loss and hot flashes"),
    ("The target", "a partial E2 band (~20–40 pg/mL) — benefit and risk cross inside it"),
], 0.95, 2.0, 11.5, 2.0, size=18, gap=11)
rrect(s, 0.95, 4.45, 11.5, 1.5, PALETEAL, radius=0.07)
textbox(s, "A 2-compartment PK + HPG-axis PD model (E2/LH/FSH turnover, Hill inhibition) makes the estrogen-threshold hypothesis a quantitative optimization — live at pipetqsp.shinyapps.io/merigolix.",
        1.2, 4.6, 11.0, 1.2, size=16.5, color=NAVY, line_spacing=1.15, anchor=MSO_ANCHOR.MIDDLE)

big_fig_slide("Merigolix: the trade-off becomes a tunable optimization",
              FIG/"s6_merigolix.png", eyebrow_txt="Flagship · the named example",
              caption="Dose sets E2 suppression; benefit (pain, lesion) and risk (BMD, flashes) cross in the target band",
              figbox=(1.1, 2.0, 11.1, 4.3))

big_fig_slide("Where this connects to NAMs: the gut–liver–kidney MPS",
              FIG/"ext/glk_validation.png", eyebrow_txt="NAMs bridge · KBCS 2026",
              caption="MPS-fitted parameters → whole-body PBPK predicted human AUC within 4-fold (GMFE 1.85×), no clinical PK input",
              figbox=(1.5, 2.0, 10.3, 4.2))
source(prs.slides[-1], "Han et al., KBCS 2026 — gut–liver–kidney MPS + whole-body PBPK")

s = add_slide(TITLE_ONLY, "Disease-state chips: distinguishing gut- from liver-driven change", title_size=21)
eyebrow(s, "NAMs bridge · KBCS 2026")
bullets(s, [
    ("Same pipeline, disease chips", "DSS colitis and cisplatin-induced kidney injury, unaffected-organ parameters locked"),
    ("Gefitinib on a colitis chip", "AUC 0.49× — gut-localised change (Fg shifts, Fa/Fh unchanged)"),
    ("Simvastatin on a cisplatin-injury chip", "AUC 2.19× — hepatic-localised (Kp,uu halved, OATP1B1 suppression)"),
    ("The payoff", "hypotheses for IBD / AKI patients usually excluded from Phase 1 — the same patients as Study ④"),
], 0.95, 2.05, 11.5, 3.2, size=17.5, gap=12)
rrect(s, 0.95, 5.4, 11.5, 1.0, PALETEAL, radius=0.08)
textbox(s, "The five new studies and the MPS/NAM work are one program: mechanism translates an in-vitro or virtual signal into a human safety decision.",
        1.2, 5.4, 11.0, 1.0, size=16.5, color=NAVY, bold=True, anchor=MSO_ANCHOR.MIDDLE, line_spacing=1.12)

s = add_slide(TITLE_ONLY, "Prerequisites for trusting a chip-derived number", title_size=22)
eyebrow(s, "Discipline for regulatory acceptance")
cards = [
    ("Chip-specific behaviour", "adsorption to plastics / membranes", TEAL),
    ("Dimensional scaling", "correct per-organ cell & fluidic factors", NAVY),
    ("Held-out clinical validation", "verify against data not used to fit", GOLD),
    ("Falsifiable correction rules", "mechanism-based, not curve-fitting", RED),
]
for i,(h,b,c) in enumerate(cards):
    cx = 0.95 + (i%2)*5.85; cy = 2.1 + (i//2)*2.0
    rrect(s, cx, cy, 5.6, 1.8, PALE, radius=0.07)
    textbox(s, f"{i+1}", cx+0.3, cy+0.25, 0.8, 0.8, size=30, color=c, bold=True)
    textbox(s, h, cx+1.15, cy+0.3, 4.2, 0.6, size=17, color=c, bold=True)
    textbox(s, b, cx+1.15, cy+1.0, 4.2, 0.7, size=14, color=INK, line_spacing=1.08)
textbox(s, "Study ④'s DOAC validation was built in exactly this style — a mechanistic (P-gp) correction, not a fitted fudge factor.",
        0.95, 6.3, 11.4, 0.6, size=15.5, color=NAVY, bold=True, line_spacing=1.1)

# ---------------- PART 6: regulatory & conclusion ----------------------
divider("Prerequisites, people, and the path forward", "6",
        "Context of use, validation against clinical data — and the workforce to translate it")

s = add_slide(TITLE_ONLY, "What regulatory acceptance will require", title_size=23)
eyebrow(s, "Context of use")
bullets(s, [
    ("Context of use", "a precise statement of the question each model/chip is qualified to answer"),
    ("Biological relevance", "the human mechanism the NAM reproduces"),
    ("Technical characterization", "reproducibility, dynamic range, failure modes"),
    ("Fit-for-purpose validation", "against held-out clinical data — the step QSP/PBPK supplies"),
], 0.95, 2.1, 11.5, 3.0, size=18, gap=13)
textbox(s, "QSP/PBPK is already an established MIDD pathway — folding NAM/MPS data into it is the credible route to decision-grade safety evidence.",
        0.95, 5.5, 11.5, 0.9, size=17, color=NAVY, bold=True, line_spacing=1.12)

s = add_slide(TITLE_ONLY, "The real bottleneck is people, not technology", title_size=23)
eyebrow(s, "Conclusion 1 — unmet needs in our own field")
bullets(s, [
    ("A shortage", "of clinical pharmacologists trained to read NAM signal through PBPK/QSP"),
    ("No standard workflow", "yet for the four NAM validation principles"),
    ("Few head-to-head validations", "of chip/organoid models against human clinical data"),
], 0.95, 2.1, 11.5, 2.2, size=18, gap=12)
rrect(s, 0.95, 4.7, 11.5, 1.1, NAVY, radius=0.08)
textbox(s, "Not a data problem, not a technology problem — a workforce problem.",
        1.1, 4.7, 11.2, 1.1, size=21, color=WHITE, bold=True, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)

s = add_slide(TITLE_ONLY, "Conclusion", title_size=26)
eyebrow(s, "What to take away")
bullets(s, [
    ("QSP + PBPK are the bridge", "they turn in-vitro and virtual signals into human exposure, safety margins, and benefit–risk"),
    ("Trade-offs become optimizations", "MS dosing interval · ICI dose ratio · incretin titration · DOAC individualization · merigolix window"),
    ("At scale", "an LLM-built, open library of 416 models — 74% already carrying a safety axis — one more each day"),
    ("Validated honestly", "the gut–liver–kidney MPS predicts human AUC within 4-fold with no clinical PK; mechanism, not fudge factors"),
    ("The gap is people", "scientists who can translate NAM-derived signals into clinical safety decisions"),
], 0.95, 2.0, 11.6, 3.9, size=17, gap=11)

# ---- Thank-you / contact
s = add_slide(TITLE_ONLY, None)
textbox(s, "Thank you", 0.95, 2.1, 11.0, 1.1, size=44, color=NAVY, bold=True)
textbox(s, "Questions & collaborators welcome", 0.98, 3.25, 11.0, 0.6, size=20, color=TEAL, bold=True)
textbox(s, [
    ("Sungpil Han, MD PhD", {"size":18,"bold":True,"color":INK}),
    ("Dept. of Pharmacology, College of Medicine, The Catholic University of Korea", {"size":14,"color":GREY}),
    ("Dept. of Clinical Pharmacology & Therapeutics, Seoul St. Mary's Hospital · PIPET", {"size":14,"color":GREY}),
    ("shan@catholic.ac.kr", {"size":15,"bold":True,"color":NAVY}),
], 0.98, 4.2, 11.0, 1.7, line_spacing=1.25)
rrect(s, 7.2, 4.35, 5.2, 0.72, PALE, radius=0.12)
textbox(s, "Fully open:  github.com/pipetcpt/qsp", 7.4, 4.35, 5.0, 0.72, size=16, color=NAVY, bold=True, anchor=MSO_ANCHOR.MIDDLE)

prs.save(str(OUT))
print("saved", OUT, "with", len(prs.slides.__iter__.__self__._sldIdLst), "slides")
print("slides:", len(list(prs.slides._sldIdLst)))
