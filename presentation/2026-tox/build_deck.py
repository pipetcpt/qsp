#!/usr/bin/env python3
"""Assemble presentation/2026-tox.pptx on the uploaded template.

Language rule for this deck
  * cover, slide titles, technical terms, tables and figures: English only
  * body text: Korean-English mix (technical terms stay in English)

Line-break rule
  The template ships with latinLnBrk="1" ("allow Latin text to wrap in the
  middle of a word" / 영어 단어 잘림 허용). Every paragraph is written with
  latinLnBrk="0" and every remaining latinLnBrk="1" in the package (master,
  layouts, presentation defaults) is rewritten to "0" after saving, so English
  words always wrap at word boundaries.

Figures are produced by the R/ggplot2 scripts in R/ (Rscript R/run_all.R);
headline numbers are read from data/*.json so text and figures stay in sync.
"""
import json
import re
import shutil
import tempfile
import zipfile
from pathlib import Path

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Emu, Inches, Pt
from lxml import etree

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
DATA = ROOT / "data"
TEMPLATE = ROOT.parent / "2026-tox-template.pptx"
OUT = ROOT.parent / "2026-tox.pptx"

LIB = json.loads((DATA / "library_summary.json").read_text())
S2 = json.loads((DATA / "s2_results.json").read_text())
S3 = json.loads((DATA / "s3_results.json").read_text())
S4 = json.loads((DATA / "s4_results.json").read_text())
S5 = json.loads((DATA / "s5_results.json").read_text())
S6 = json.loads((DATA / "s6_results.json").read_text())
N_MODELS = LIB["n_gallery"]
PCT_SAFE = LIB["pct_safety"]

# ---- palette ----------------------------------------------------------------
NAVY = RGBColor(0x16, 0x3A, 0x5F)
TEAL = RGBColor(0x1F, 0x8A, 0x8C)
GOLD = RGBColor(0xC8, 0x92, 0x2A)
RED = RGBColor(0xC0, 0x39, 0x2B)
GREEN = RGBColor(0x2E, 0x8B, 0x6F)
BLUE = RGBColor(0x2B, 0x6C, 0xB0)
INK = RGBColor(0x1A, 0x22, 0x2C)
GREY = RGBColor(0x5F, 0x6B, 0x76)
LGREY = RGBColor(0xD7, 0xDE, 0xE5)
PALE = RGBColor(0xEE, 0xF3, 0xF8)
PALETEAL = RGBColor(0xE4, 0xF1, 0xEB)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

# ---- typefaces (S-Core Dream: the template's own family; Hangul + Latin) ----
F_REG = "S-Core Dream 4 Regular"
F_MED = "S-Core Dream 5 Medium"
F_BOLD = "S-Core Dream 6 Bold"
F_MATH = "Cambria Math"

HANGUL = re.compile(r"[가-힣]")

prs = Presentation(str(TEMPLATE))
LAY = {L.name: L for L in prs.slide_layouts}
TITLE_ONLY = LAY.get("제목만") or LAY.get("Title Only")
BLANK = LAY.get("빈 화면") or TITLE_ONLY

# ---- drop the template's placeholder content slide (index 1) -----------------
sld_ids = list(prs.slides._sldIdLst)
if len(sld_ids) > 1:
    prs.part.drop_rel(sld_ids[1].get(qn("r:id")))
    prs.slides._sldIdLst.remove(sld_ids[1])


# =============================================================================
# low-level text helpers
# =============================================================================
def set_font(run, family, size=None, color=None, italic=False):
    f = run.font
    if size is not None:
        f.size = Pt(size)
    f.bold = False            # weight comes from the family, never faux-bold
    f.italic = italic
    if color is not None:
        f.color.rgb = color
    f.name = family           # <a:latin>
    rPr = run._r.get_or_add_rPr()
    for tag in ("a:ea", "a:cs"):
        for el in rPr.findall(qn(tag)):
            rPr.remove(el)
    latin = rPr.find(qn("a:latin"))
    ea = etree.SubElement(rPr, qn("a:ea")); ea.set("typeface", family)
    cs = etree.SubElement(rPr, qn("a:cs")); cs.set("typeface", family)
    latin.addnext(ea); ea.addnext(cs)
    has_ko = bool(HANGUL.search(run.text or ""))
    rPr.set("lang", "ko-KR" if has_ko else "en-US")
    rPr.set("altLang", "en-US" if has_ko else "ko-KR")


def no_word_break(p):
    """latinLnBrk=0 on this paragraph: never split an English word."""
    pPr = p._p.get_or_add_pPr()
    pPr.set("latinLnBrk", "0")
    pPr.set("eaLnBrk", "1")
    pPr.set("hangingPunct", "1")


def add_run(p, text, family=F_REG, size=17, color=INK, italic=False):
    r = p.add_run()
    r.text = text
    set_font(r, family, size, color, italic)
    return r


def _frame(tb, anchor=MSO_ANCHOR.TOP, margins=0.0):
    tf = tb.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    for m in ("margin_left", "margin_right", "margin_top", "margin_bottom"):
        setattr(tf, m, Inches(margins))
    return tf


def textbox(s, lines, x, y, w, h, size=17, color=INK, family=F_REG, align=PP_ALIGN.LEFT,
            anchor=MSO_ANCHOR.TOP, italic=False, spacing=1.15, gap=0):
    """lines: str | list of (str | [(text, {opts}), ...] | (text, {opts}))."""
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = _frame(tb, anchor)
    if not isinstance(lines, list):
        lines = [lines]
    for i, line in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.line_spacing = spacing
        if gap:
            p.space_after = Pt(gap)
        no_word_break(p)
        segs = line if isinstance(line, list) else [line]
        for seg in segs:
            txt, o = (seg if isinstance(seg, tuple) else (seg, {}))
            add_run(p, txt, o.get("family", family), o.get("size", size),
                    o.get("color", color), o.get("italic", italic))
    return tb


def bullets(s, items, x, y, w, h, size=17, gap=10, lead_color=NAVY, spacing=1.12):
    """items: (lead, body) or plain str. Real hanging-indent bullets."""
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = _frame(tb)
    for i, it in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.line_spacing = spacing
        p.space_after = Pt(gap)
        no_word_break(p)
        pPr = p._p.get_or_add_pPr()
        ind = Emu(Inches(0.28))
        pPr.set("marL", str(ind)); pPr.set("indent", str(-ind))
        buClr = etree.SubElement(pPr, qn("a:buClr"))
        etree.SubElement(buClr, qn("a:srgbClr")).set("val", "1F8A8C")
        etree.SubElement(pPr, qn("a:buFont")).set("typeface", "Arial")
        etree.SubElement(pPr, qn("a:buChar")).set("char", "•")
        lead, body = it if isinstance(it, tuple) else (None, it)
        if lead:
            add_run(p, lead + "  ", F_BOLD, size, lead_color)
        add_run(p, body, F_REG, size, INK)
    return tb


def rrect(s, x, y, w, h, fill, line=None, radius=0.08, shape=MSO_SHAPE.ROUNDED_RECTANGLE):
    sp = s.shapes.add_shape(shape, Inches(x), Inches(y), Inches(w), Inches(h))
    sp.fill.solid(); sp.fill.fore_color.rgb = fill
    if line is None:
        sp.line.fill.background()
    else:
        sp.line.color.rgb = line; sp.line.width = Pt(1.2)
    sp.shadow.inherit = False
    if shape == MSO_SHAPE.ROUNDED_RECTANGLE:
        sp.adjustments[0] = radius
    if sp.has_text_frame:
        for p in sp.text_frame.paragraphs:
            no_word_break(p)
    return sp


def add_image(s, path, bx, by, bw, bh, align="center"):
    iw, ih = Image.open(path).size
    ar = iw / ih
    if bw / bh > ar:
        h = bh; w = bh * ar
    else:
        w = bw; h = bw / ar
    x = bx + (bw - w) / 2 if align == "center" else bx
    y = by + (bh - h) / 2
    s.shapes.add_picture(str(path), Inches(x), Inches(y), Inches(w), Inches(h))
    return x, y, w, h


def source(s, text):
    textbox(s, text, 4.9, 7.02, 7.75, 0.3, size=9.5, color=GREY, align=PP_ALIGN.RIGHT, italic=True,
            spacing=1.0)


def eyebrow(s, text, y=1.3):
    textbox(s, text.upper(), 0.75, y, 11.5, 0.32, size=12, color=GOLD, family=F_BOLD, spacing=1.0)


def new_slide(title, eyebrow_txt=None, size=26):
    s = prs.slides.add_slide(TITLE_ONLY)
    ph = next(p for p in s.placeholders if p.placeholder_format.idx == 0)
    tf = ph.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.LEFT
    no_word_break(p)
    add_run(p, title, F_MED, size, NAVY)
    if eyebrow_txt:
        eyebrow(s, eyebrow_txt)
    return s


def callout(s, text, x=0.75, y=5.55, w=11.85, h=0.95, fill=PALETEAL, color=NAVY, size=17):
    rrect(s, x, y, w, h, fill, radius=0.1)
    textbox(s, text, x + 0.3, y, w - 0.6, h, size=size, color=color, family=F_MED,
            anchor=MSO_ANCHOR.MIDDLE, spacing=1.15)


def caption(s, text, y=6.15, size=15):
    textbox(s, text, 0.75, y, 11.85, 0.5, size=size, color=NAVY, family=F_MED,
            align=PP_ALIGN.CENTER, spacing=1.1)


def table(s, rows, x, y, w, col_w, row_h=0.5, size=14, header_fill=NAVY, first_col_bold=True):
    nr, nc = len(rows), len(rows[0])
    gf = s.shapes.add_table(nr, nc, Inches(x), Inches(y), Inches(w), Inches(row_h * nr))
    tbl = gf.table
    tbl.first_row = True; tbl.horz_banding = False
    tot = sum(col_w)
    for j, cw in enumerate(col_w):
        tbl.columns[j].width = Inches(w * cw / tot)
    for i, row in enumerate(rows):
        tbl.rows[i].height = Inches(row_h)
        for j, val in enumerate(row):
            c = tbl.cell(i, j)
            c.margin_left = Inches(0.12); c.margin_right = Inches(0.1)
            c.margin_top = Inches(0.03); c.margin_bottom = Inches(0.03)
            c.vertical_anchor = MSO_ANCHOR.MIDDLE
            c.fill.solid()
            if i == 0:
                c.fill.fore_color.rgb = header_fill
                fam, col = F_BOLD, WHITE
            else:
                c.fill.fore_color.rgb = WHITE if i % 2 else PALE
                fam = F_BOLD if (j == 0 and first_col_bold) else F_REG
                col = NAVY if (j == 0 and first_col_bold) else INK
            tf = c.text_frame; tf.word_wrap = True
            p = tf.paragraphs[0]; p.alignment = PP_ALIGN.LEFT
            no_word_break(p)
            add_run(p, str(val), fam, size, col)
    return gf


def stat_cards(s, stats, x0=0.75, y=1.85, total_w=11.85, h=1.85, num_size=40, lab_size=14.5):
    n = len(stats); gap = 0.3
    cw = (total_w - gap * (n - 1)) / n
    for i, (num, lab, col) in enumerate(stats):
        cx = x0 + i * (cw + gap)
        rrect(s, cx, y, cw, h, PALE, radius=0.08)
        textbox(s, num, cx + 0.15, y + 0.18, cw - 0.3, 0.85, size=num_size, color=col, family=F_BOLD,
                align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE, spacing=1.0)
        textbox(s, lab, cx + 0.2, y + 1.08, cw - 0.4, h - 1.15, size=lab_size, color=INK,
                align=PP_ALIGN.CENTER, spacing=1.12)


def cards(s, items, cols=2, x0=0.75, y0=1.85, total_w=11.85, h=1.75, gap=0.3, body_size=15.5,
          numbered=False):
    cw = (total_w - gap * (cols - 1)) / cols
    for i, (head, body, col) in enumerate(items):
        cx = x0 + (i % cols) * (cw + gap)
        cy = y0 + (i // cols) * (h + gap)
        rrect(s, cx, cy, cw, h, PALE, radius=0.07)
        tx = cx + 0.3
        if numbered:
            textbox(s, str(i + 1), cx + 0.28, cy + 0.2, 0.6, 0.7, size=32, color=col, family=F_BOLD,
                    spacing=1.0)
            tx = cx + 1.0
        textbox(s, head, tx, cy + 0.24, cw - (tx - cx) - 0.25, 0.5, size=18, color=col,
                family=F_BOLD, spacing=1.0)
        textbox(s, body, tx, cy + 0.8, cw - (tx - cx) - 0.25, h - 0.9, size=body_size, color=INK,
                spacing=1.15)


def divider(num, title, subtitle):
    s = prs.slides.add_slide(BLANK)
    rrect(s, 0, 2.05, 13.333, 2.55, NAVY, shape=MSO_SHAPE.RECTANGLE)
    textbox(s, f"PART {num}", 0.95, 2.35, 6.0, 0.45, size=16, color=GOLD, family=F_BOLD, spacing=1.0)
    textbox(s, title, 0.95, 2.85, 11.4, 1.5, size=36, color=WHITE, family=F_MED, spacing=1.05,
            anchor=MSO_ANCHOR.TOP)
    textbox(s, subtitle, 0.97, 4.9, 11.4, 1.0, size=18, color=GREY, spacing=1.2)
    return s


def fig_slide(title, fig, eyebrow_txt=None, cap=None, box=(0.75, 1.75, 11.85, 4.3), src=None):
    s = new_slide(title, eyebrow_txt)
    add_image(s, fig, *box)
    if cap:
        caption(s, cap)
    if src:
        source(s, src)
    return s


def intro_slide(title, eyebrow_txt, items, layer, src):
    s = new_slide(title, eyebrow_txt)
    bullets(s, items, 0.75, 1.95, 11.85, 2.5, size=18, gap=12)
    rrect(s, 0.75, 4.1, 11.85, 1.6, PALETEAL, radius=0.08)
    textbox(s, [[("QSP layer   ", {"family": F_BOLD, "color": TEAL}), (layer, {})]],
            1.05, 4.1, 11.25, 1.6, size=17, color=NAVY, anchor=MSO_ANCHOR.MIDDLE, spacing=1.2)
    source(s, src)
    return s


# =============================================================================
# COVER (template slide 0) -- English only
# =============================================================================
cover = prs.slides[0]
for sh in cover.shapes:
    if not sh.has_text_frame:
        continue
    tf = sh.text_frame
    if "한국독성학회" in tf.text:
        # paragraph 0 = society line <br> session line <br> title line 1:
        # replace only the runs before the first line break
        p0 = tf.paragraphs[0]._p
        society = []
        for ch in list(p0):
            if ch.tag == qn("a:br"):
                break
            if ch.tag == qn("a:r"):
                society.append(ch)
        society[0].find(qn("a:t")).text = "Korean Society of Toxicology (KSOT) · 2026 Fall Meeting (42nd)"
        society[0].find(qn("a:rPr")).set("lang", "en-US")
        for r in society[1:]:
            p0.remove(r)
    if tf.text.strip() == "발표자료 다운로드":
        r = tf.paragraphs[0].runs
        r[0].text = "Download slides"
        for extra in r[1:]:
            extra._r.getparent().remove(extra._r)
        r[0]._r.get_or_add_rPr().set("lang", "en-US")
# layout 1 carries the same Korean caption under its QR code
for L in prs.slide_layouts:
    for sh in L.shapes:
        if sh.has_text_frame and sh.text_frame.text.strip() == "발표자료 다운로드":
            rs = sh.text_frame.paragraphs[0].runs
            rs[0].text = "Download slides"
            for extra in rs[1:]:
                extra._r.getparent().remove(extra._r)

# =============================================================================
# PART 1 -- the safety problem and the regulatory shift
# =============================================================================
s = new_slide("Drug development fails late, and toxicity leads", "Why we are here")
stat_cards(s, [("~90%", "임상시험에 진입한 후보물질 중 승인에 실패하는 비율", RED),
               ("10–15 yr", "신약 하나당 개발 기간", NAVY),
               ("$1–2.6 B", "신약 하나당 개발 비용", NAVY)])
bullets(s, [
    ("Unexpected toxicity", "late-stage attrition과 post-marketing withdrawal의 주요 원인"),
    ("Healthy-volunteer Phase 1", "젊고 건강한 사람의 exposure만 확인 — 실제 약을 복용할 환자와 다름"),
    ("The costliest failures", "취약한 실제 환자에서 뒤늦게 드러나는 safety signal"),
], 0.75, 4.15, 11.85, 2.4, size=18, gap=12)
source(s, "DiMasi J Health Econ 2016 · Kola & Landis Nat Rev Drug Discov 2004")

s = new_slide("Conventional safety assessment misses vulnerable patients", "The translation gap", size=25)
cw = (11.85 - 0.6) / 3
for i, (h, b, c) in enumerate([
        ("2D cell culture", "Tissue architecture와 organ-level function이 없음", GREY),
        ("Animal models", "Drug disposition과 response의 interspecies difference가 지속됨", GOLD),
        ("Phase 1 trials", "Healthy volunteer 위주 — hepatic/renal impairment, pregnancy, elderly, rare disease 환자는 제외", RED)]):
    cx = 0.75 + i * (cw + 0.3)
    rrect(s, cx, 1.85, cw, 2.75, PALE, radius=0.07)
    textbox(s, h, cx + 0.3, 2.1, cw - 0.6, 0.5, size=19, color=c, family=F_BOLD, spacing=1.0)
    textbox(s, b, cx + 0.3, 2.75, cw - 0.6, 1.8, size=16, spacing=1.2)
callout(s, "결과적으로 in vitro signal과 animal data가 실제 위험 환자의 human exposure·response로 잘 translate되지 않음",
        y=5.0, h=1.05)
source(s, "Van Norman, JACC Basic Transl Sci 2019")

s = new_slide("The regulatory shift: from animal mandate to NAMs", "What changed")
table(s, [
    ["Year", "Milestone", "What it means"],
    ["1938", "Federal Food, Drug, and Cosmetic Act", "Animal testing de facto required before human use"],
    ["2022", "FDA Modernization Act 2.0", "“Animal tests” replaced by “nonclinical tests”: in vitro and in silico allowed"],
    ["2025", "FDA Roadmap to Reducing Animal Testing", "Phased rollout, starting with monoclonal antibodies"],
    ["Now", "New approach methodologies (NAMs)", "MPS / organ-chips, organoids, in silico models"],
], 0.75, 1.85, 11.85, col_w=[1.0, 4.0, 6.0], row_h=0.62, size=14.5)
callout(s, "다만 in vitro signal이 decision-grade evidence가 되려면, mechanistic model이 이를 human exposure와 response로 translate해야 함",
        y=5.25, h=1.05)
source(s, "FDA Modernization Act 2.0, Pub. L. 117-328 (2022) · FDA Roadmap (2025)")

s = fig_slide("NAMs need a translational bridge", FIG / "ext/nam_loop.png", "The NAM translation loop",
              cap="PBPK와 QSP가 NAM data를 qualify하고 human으로 translate하는 핵심 단계 — pattern-matching만으로는 대체 불가",
              box=(1.4, 1.85, 10.55, 4.2), src="Cao Y, Clin Pharmacol Ther 2025")

# =============================================================================
# PART 2 -- what QSP / PBPK is
# =============================================================================
divider(2, "QSP & PBPK: the mechanistic bridge",
        "Drug exposure를 target, pathway, organ, patient outcome까지 수학적으로 연결")

s = new_slide("QSP and PBPK describe the causal chain", "What it is")
rrect(s, 0.75, 1.85, 11.85, 0.95, NAVY, radius=0.14)
textbox(s, "Drug  →  Target  →  Pathway  →  Disease  →  Patient", 0.75, 1.85, 11.85, 0.95,
        size=26, color=WHITE, family=F_MED, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE, spacing=1.0)
for i, (h, b, c) in enumerate([
        ("PBPK", "약물이 human body의 어디에, 어떤 concentration으로 분포하는가?", BLUE),
        ("QSP", "그 concentration에서 disease와 safety는 어떻게 변하는가?", TEAL)]):
    cx = 0.75 + i * 6.075
    rrect(s, cx, 3.15, 5.775, 1.75, PALE, radius=0.07)
    textbox(s, h, cx + 0.35, 3.35, 5.0, 0.55, size=22, color=c, family=F_BOLD, spacing=1.0)
    textbox(s, b, cx + 0.35, 4.0, 5.1, 0.8, size=16.5, spacing=1.15)
callout(s, "둘을 결합하면 first-in-human dose, safety margin, DDI risk, special population response, chip-derived parameter의 IVIVE까지 예측 가능",
        y=5.25, h=1.05)

s = new_slide("Anatomy of a QSP model: a coupled nonlinear system", "Under the hood", size=25)
rrect(s, 2.4, 1.85, 8.55, 0.95, PALETEAL, radius=0.14)
textbox(s, "dx/dt = f( x(t), θ, u(t) ),     x(0) = x₀", 2.4, 1.85, 8.55, 0.95, size=24, color=NAVY,
        family=F_MATH, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE, spacing=1.0)
for y, (lead, body) in zip([3.05, 3.95, 4.55], [
        ("States  x", "drug concentration · receptor occupancy · signaling species · cell population · biomarker · clinical endpoint"),
        ("Parameters  θ", "rate constant · EC₅₀ / Emax · Hill coefficient · clearance · binding affinity"),
        ("Inputs  u(t)", "dosing regimen — 치료에서 우리가 조절할 수 있는 유일한 handle")]):
    textbox(s, lead, 0.95, y, 2.9, 0.6, size=19, color=TEAL, family=F_BOLD, spacing=1.0)
    textbox(s, body, 3.9, y, 8.6, 0.8, size=16.5, spacing=1.15)
callout(s, f"Library의 각 질환 모델은 {LIB['min_odes']}–{LIB['max_odes']}개의 coupled ODE (median {LIB['median_odes']})와 수십–수백 개 parameter로 구성",
        y=5.35, h=0.7, size=16)

s = new_slide("The same nonlinear motifs recur across diseases", "A shared grammar")
table(s, [
    ["Motif", "Governing idea", "Typical use"],
    ["Hill / Emax", "E(C) = Emax·Cⁿ / (EC₅₀ⁿ + Cⁿ)", "Any drug dose–response"],
    ["Indirect response", "dR/dt = k_in·f(·) − k_out·R", "Biomarker turnover"],
    ["Mass action / TMDD", "Drug–target binding, target-mediated clearance", "mAbs, saturable drugs"],
    ["Logistic growth", "dN/dt = rN(1 − N/K) − kill(C)·N", "Tumor burden under therapy"],
    ["Power-law feedback", "e.g. EPO ∝ (Hb₀/Hb)^γ", "Homeostatic loops"],
], 0.75, 1.85, 11.85, col_w=[2.3, 4.6, 3.0], row_h=0.56, size=15)
callout(s, f"{N_MODELS}개 모델 전체, 그리고 모든 chip-to-human bridge가 같은 grammar를 공유", y=5.6, h=0.8)

s = new_slide("Why it matters: model-informed drug development", "MIDD, endorsed by FDA and EMA")
table(s, [
    ["Stage", "What QSP / PBPK contributes", "Safety payoff"],
    ["Target validation", "Does target modulation reach the phenotype?", "Drops implausible targets early"],
    ["First-in-human dose", "Exposure–response and safety margins from mechanism", "Safer, leaner Phase 1"],
    ["Trial design", "Virtual patients for enrichment and endpoints", "Fewer patients exposed to risk"],
    ["Special populations", "Hepatic / renal impairment, age, pregnancy", "Doses the excluded patient"],
    ["Combinations", "Synergy and antagonism in silico", "Rational, safer polypharmacy"],
], 0.75, 1.85, 11.85, col_w=[2.4, 4.6, 3.2], row_h=0.56, size=14.5)
callout(s, "QSP는 efficacy–safety trade-off를 tunable optimization problem으로 바꿈 — 오늘 case study를 관통하는 주제",
        y=5.6, h=0.8)

# =============================================================================
# PART 3 -- the living QSP library
# =============================================================================
divider(3, "A living QSP model library",
        "LLM-augmented workflow가 literature-grounded, quality-gated disease model을 하루 하나씩 구축")

s = new_slide("The bottleneck: one model takes expert-months", "Why there were so few")
bullets(s, [
    ("Deep synthesis", "한 질환의 pathophysiology와 pharmacology를 깊이 통합해야 함"),
    ("Careful formulation", "ODE structure, parameterization, implementation, calibration"),
    ("The cost", "질환당 expert-months to years — 소수의 질환만 QSP model을 가짐"),
    ("The consequence", "모델이 silo화·비공개·재현 불가; long-tail disease는 모델링되지 않음"),
], 0.75, 2.1, 11.85, 2.7, size=19, gap=14)
rrect(s, 0.75, 4.9, 11.85, 1.0, NAVY, radius=0.1)
textbox(s, "One disease ≈ one expert-year   ·   scarce, slow, closed", 0.75, 4.9, 11.85, 1.0, size=21,
        color=WHITE, family=F_MED, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE, spacing=1.0)

s = new_slide("An autonomous LLM agent builds the models", "The engine")
textbox(s, "Autonomous LLM coding agent가 quality-gated, referenced QSP model을 하루 한 질환씩 완성 — 모든 산출물은 version control로 관리",
        0.75, 1.8, 11.85, 0.85, size=17, spacing=1.2)
for i, (h, b, c) in enumerate([("Mechanistic map", "Graphviz network, ≥100 nodes", TEAL),
                               ("mrgsolve ODE model", "PK ⊕ disease PD, ≥15 states", NAVY),
                               ("Shiny dashboard", "≥6 interactive tabs", GOLD),
                               ("Curated references", "≥30 PubMed, calibration notes", GREEN)]):
    cx = 0.75 + (i % 2) * 3.2
    cy = 2.95 + (i // 2) * 1.7
    rrect(s, cx, cy, 3.0, 1.5, PALE, radius=0.08)
    textbox(s, h, cx + 0.25, cy + 0.22, 2.6, 0.5, size=16.5, color=c, family=F_BOLD, spacing=1.0)
    textbox(s, b, cx + 0.25, cy + 0.75, 2.6, 0.7, size=14, spacing=1.12)
ix, iy, iw, ih = add_image(s, FIG / "ext/map_igan.png", 7.25, 2.95, 5.35, 3.2)
textbox(s, "Example: IgA nephropathy mechanistic map", ix, iy + ih + 0.08, iw, 0.3, size=11.5,
        color=GREY, align=PP_ALIGN.CENTER, italic=True, spacing=1.0)
source(s, "github.com/pipetcpt/qsp (fully open)")

fig_slide("The daily build loop", FIG / "architecture.png", "One disease, end-to-end, every run",
          cap="Select → research → map → ODE model → dashboard → references → quality gate → commit (git stop-hook로 강제)",
          box=(0.9, 1.75, 11.55, 4.3))

s = new_slide("Guardrails for LLM-written models", "Structure + grounding + verification")
cards(s, [
    ("Quality gates", "Map ≥100 nodes / ≥8 clusters · ODE ≥15 compartments / ≥5 scenarios · ≥30 references", TEAL),
    ("Structural templates", "Hill, turnover, TMDD 같은 reusable motif로 failure space를 줄임", NAVY),
    ("Literature grounding", "모든 parameter set에 calibration trial을 인용 — 핵심 anti-hallucination 장치", GOLD),
    ("Reproducibility", "Diffable·auditable·re-runnable code in git; 사용 전 human review", GREEN),
], h=1.95)

s = new_slide(f"Coverage: {LIB['n_areas']} therapeutic areas, one safety atlas", "The library as a safety atlas")
for i, (num, lab, col) in enumerate([
        (str(N_MODELS), "disease QSP models", NAVY),
        (str(LIB["n_areas"]), "therapeutic areas", NAVY),
        (f"{PCT_SAFE}%", "explicit safety / toxicity axis", TEAL)]):
    cy = 1.85 + i * 1.47
    rrect(s, 0.75, cy, 3.1, 1.3, PALE, radius=0.08)
    textbox(s, num, 0.75, cy + 0.1, 3.1, 0.7, size=34, color=col, family=F_BOLD, align=PP_ALIGN.CENTER,
            anchor=MSO_ANCHOR.MIDDLE, spacing=1.0)
    textbox(s, lab, 0.9, cy + 0.82, 2.8, 0.4, size=13, color=INK, align=PP_ALIGN.CENTER, spacing=1.0)
add_image(s, FIG / "s1_areas.png", 4.15, 1.8, 8.45, 4.3)
caption(s, "Safety axis: model code·map·README·references에 서로 다른 toxicology term(toxicity, QT, DILI, PML 등)이 3개 이상 명시된 모델",
        y=6.22, size=13)

fig_slide("Quality gates hold across the library", FIG / "s1_quality.png", "Size & grounding",
          cap=f"{LIB['n_odes_gate']}/{LIB['n_mined']} 모델이 ≥15 ODE gate (median {LIB['median_odes']}), {LIB['n_refs_gate']}/{LIB['n_mined']}가 ≥30 reference gate 통과 (median {LIB['median_refs']})",
          box=(0.9, 1.75, 11.55, 4.3))

fig_slide("One consistent schema across the library", FIG / "s1_atlas.png", "One framework",
          cap="모든 모델이 같은 4-artifact schema — ODE 수 × map 크기 × reference 수로 한눈에 비교 가능",
          box=(1.4, 1.75, 10.55, 4.3))

# =============================================================================
# PART 4 -- five case studies
# =============================================================================
divider(4, "Five safety case studies",
        "Autoimmune · oncology · obesity · aging · gynecology — efficacy–safety trade-off를 optimization으로")

s = new_slide("Five case studies, one pattern", "Roadmap")
table(s, [
    ["Study", "Disease area", "Tunable handle", "Safety axis"],
    ["① Multiple sclerosis", "Autoimmune / neuro", "Natalizumab dosing interval", "PML risk"],
    ["② Checkpoint inhibitors", "Oncology (melanoma)", "Nivolumab : ipilimumab ratio", "Grade 3–4 irAE"],
    ["③ Incretin therapy", "Obesity + aging", "Target dose and escalation", "GI AE · lean-mass loss"],
    ["④ DOACs in the elderly", "Aging / organ function", "Age- and renal-based dose", "Bleeding"],
    ["★ Merigolix", "Gynecology", "GnRH antagonist dose", "Bone loss · hot flashes"],
], 0.75, 1.85, 11.85, col_w=[2.9, 2.6, 3.4, 2.6], row_h=0.56, size=14.5)
callout(s, "공통 질문: dose·schedule이라는 handle을 어떻게 조절해야 efficacy를 지키면서 safety risk를 줄일 수 있는가?",
        y=5.6, h=0.8)

# ---- Study 1: MS -------------------------------------------------------------
intro_slide("Study ①  Multiple sclerosis: natalizumab vs PML", "Autoimmune / neurology", [
    ("Natalizumab (anti-α4-integrin)", "relapse 억제 효과가 가장 강력한 DMT 중 하나 — 그러나 CNS immune surveillance를 낮춤"),
    ("The risk", "JC virus reactivation에 의한 progressive multifocal leukoencephalopathy (PML)"),
    ("The question", "dosing interval 조절로 efficacy는 유지하고 PML risk만 줄일 수 있는가?"),
], "dosing interval 동안의 α4 receptor saturation → relapse-rate suppression, 그리고 stratified PML hazard (anti-JCV index · prior immunosuppressant · duration)에 extended-interval hazard ratio 적용",
    "Polman NEJM 2006 · Ho Lancet Neurol 2017 · Foley Lancet Neurol 2022 · Ryerson Neurology 2019")

fig_slide("The dosing interval is the controllable handle", FIG / "s2_saturation.png", "Study ① · MS",
          cap=f"q6w extended-interval dosing (EID)에서도 relapse reduction {S2['red_q6w_pct']:.0f}% 유지 (q4w {S2['red_q4w_pct']:.0f}%)")
fig_slide("PML risk is stratified, and steep", FIG / "s2_pml_table.png", "Study ① · MS",
          cap="anti-JCV index, prior immunosuppressant, 치료 기간에 따라 PML risk가 급격히 증가",
          box=(1.4, 1.75, 10.55, 4.3), src="Ho PR et al., Lancet Neurol 2017")
fig_slide(f"Extended-interval dosing cuts PML risk by {S2['pml_drop_pct']}%", FIG / "s2_benefit_risk.png",
          "Study ① · the optimization",
          cap=f"6년간 relapse {S2['relapse_cost_6y']:.1f}회 증가를 대가로 high-risk stratum PML {S2['pml_6y_highrisk_sid']:.1f} → {S2['pml_6y_highrisk_eid']:.1f} /1000",
          box=(1.3, 1.75, 10.75, 4.3))
fig_slide("Where natalizumab q6w sits among MS therapies", FIG / "s2_dmt_map.png", "Study ① · MS",
          cap="EID는 natalizumab을 ocrelizumab에 가까운 efficacy로, 훨씬 낮은 risk 쪽으로 이동 (safety index는 illustrative)",
          box=(1.3, 1.75, 10.75, 4.3),
          src="Polman NEJM 2006 · Kappos NEJM 2010 · Hauser NEJM 2017 · Ryerson Neurology 2019")

# ---- Study 2: ICI ------------------------------------------------------------
intro_slide("Study ②  Checkpoint inhibitors: two antibodies", "Oncology · metastatic melanoma", [
    ("Nivolumab + ipilimumab", "melanoma에서 가장 효과적이면서 가장 독성이 큰 regimen"),
    ("Key biology", "ipilimumab의 anti-tumor priming 효과는 low dose에서 포화되지만, toxicity는 dose와 함께 계속 증가"),
    ("The question", "왜 ipilimumab을 낮추면 (CheckMate 511) efficacy는 유지되고 toxicity는 줄어드는가?"),
], "PD-1·CTLA-4 occupancy → anti-tumor / self-reactive effector pool (Bliss synergy) → ORR과 grade 3–4 irAE — 같은 CTLA-4 occupancy를 서로 다른 두 curve로 읽음",
    "Larkin NEJM 2015, 2019 (CheckMate 067) · Lebbé JCO 2019 (CheckMate 511) · Weber JCO 2017")

r13 = S3["regimens"]["NIVO1 / IPI3"]; r31 = S3["regimens"]["NIVO3 / IPI1"]
fig_slide("The efficacy–safety frontier of combination dosing", FIG / "s3_pareto.png",
          "Study ② · the optimization",
          cap=f"NIVO1/IPI3 → NIVO3/IPI1: grade 3–4 irAE 약 {S3['cm511_tox_drop_pts']:.0f}%p 감소, ORR은 유지 ({r13['orr']:.0f}% → {r31['orr']:.0f}%)",
          box=(1.1, 1.75, 11.15, 4.3))
s = new_slide("Response and toxicity run on two clocks", "Study ② · oncology")
add_image(s, FIG / "s3_timecourse.png", 0.75, 1.8, 5.75, 4.25)
add_image(s, FIG / "s3_organ_timing.png", 6.8, 1.8, 5.8, 4.25)
caption(s, "Tumor response가 먼저, irAE는 뒤이어 organ별로 다른 시점에 발생 — safety monitoring을 위한 kinetic signature")
source(s, "Weber JS et al., J Clin Oncol 2017")

# ---- Study 3: obesity --------------------------------------------------------
intro_slide("Study ③  Incretin therapy: a three-way optimization", "Obesity / metabolic · aging subgroup", [
    ("Semaglutide & tirzepatide", "조절 가능한 handle은 target dose와 escalation schedule"),
    ("Three competing outcomes", "weight loss (benefit) vs GI adverse event (tolerability) vs lean-mass loss (safety)"),
    ("The aging angle", "고령자에서 lean-mass loss는 sarcopenia·frailty risk로 직결"),
], "exposure → appetite suppression → energy-balance ODE; 체중 변화를 fat / fat-free mass로 분배 (age-sensitive Forbes rule); GI AE = escalation이 유발하고 tolerance가 완화하는 nausea state",
    "Wilding NEJM 2021 (STEP 1) · Jastreboff NEJM 2022 (SURMOUNT-1)")

fig_slide("Titration is the handle: fast escalation buys nausea", FIG / "s4_titration.png", "Study ③ · obesity",
          cap="2-weekly escalation은 weight loss 이득 없이 초기 nausea peak만 키움")
fig_slide("More weight loss costs more GI burden", FIG / "s4_frontier.png", "Study ③ · the optimization",
          cap=f"72주 weight loss: semaglutide 2.4 mg {S4['sema24_pct']}% (STEP 1: 14.9%), tirzepatide 15 mg {S4['tirze15_pct']}% (SURMOUNT-1: 20.9%)",
          box=(1.2, 1.75, 10.95, 4.3))
s = new_slide("The aging safety axis: lean-mass loss", "Study ③ · obesity + aging")
add_image(s, FIG / "s4_composition.png", 0.75, 1.8, 5.8, 4.25)
add_image(s, FIG / "s4_aging.png", 6.8, 1.8, 5.8, 4.25)
caption(s, f"Weight loss 중 lean mass 비율이 45세 {S4['lean_frac_young']}% → 70세 {S4['lean_frac_old']}%로 증가 — 고령자 dose 조절 signal")

# ---- Study 4: DOAC -----------------------------------------------------------
intro_slide("Study ④  DOACs in the elderly: dosing the excluded", "Aging / renal and hepatic impairment", [
    ("The gap", "Phase 1은 젊고 renal function이 정상인 volunteer에서 DOAC exposure를 확인"),
    ("The reality", "실제 DOAC 복용 환자는 고령·renal impairment — 바로 임상시험에서 제외된 집단"),
    ("The question", "PBPK로 연구되지 않은 환자의 exposure를 예측하고 dose를 individualize할 수 있는가?"),
], f"age → creatinine clearance → DOAC exposure (renal-fraction scaling) → bleeding·stroke risk; {S5['n_cohort']:,}명 virtual elderly cohort에서 fixed vs individualized dosing 비교",
    "US product labels (dabigatran, apixaban, rivaroxaban, edoxaban) · Rowland & Tozer")

fig_slide("Age is a renal-function covariate", FIG / "s5_age_crcl.png", "Study ④ · aging",
          cap=f"Population mean CrCl은 {S5['age_cross_crcl50']}세에 dose-reduction threshold (50 mL/min) 아래로 — 측정 전에 exposure는 이미 상승",
          box=(1.2, 1.75, 10.95, 4.3))

s = new_slide("Validating exposure, and finding the missing mechanism", "Study ④ · aging", size=25)
add_image(s, FIG / "s5_validation.png", 0.6, 1.8, 7.0, 4.55)
bullets(s, [
    ("Renal scaling alone", f"label의 renal-impairment PK를 GMFE {S5['gmfe_renal_only']:.2f}×로 예측"),
    ("Dabigatran under-predicted", f"moderate RI에서 observed {S5['dabi_moderate_obs']}× vs predicted {S5['dabi_moderate_pred_renal']:.1f}× — GFR만으로 설명 불가"),
    ("Mechanism, not curve-fitting", f"P-gp / absorption term 추가 → GMFE {S5['gmfe_corrected']:.2f}×"),
], 7.85, 2.4, 4.75, 3.8, size=16, gap=14)
source(s, "US product labels: dabigatran, apixaban, rivaroxaban, edoxaban (renal-impairment PK)")

s = new_slide("Two risks oppose; individualized dosing wins", "Study ④ · the optimization")
add_image(s, FIG / "s5_exposure_response.png", 0.75, 1.8, 5.75, 4.25)
add_image(s, FIG / "s5_individualization.png", 6.8, 1.8, 5.8, 4.25)
caption(s, f"Model-informed dose reduction: over-exposure {S5['overexposed_fixed_pct']}% → {S5['overexposed_adj_pct']}%, modeled major bleeding {S5['bleed_fixed']} → {S5['bleed_adj']} /100 PY")

# =============================================================================
# PART 5 -- flagship + NAMs bridge
# =============================================================================
divider(5, "From optimization to the clinic, and to NAMs",
        "Flagship merigolix case, 그리고 gut–liver–kidney MPS + PBPK bridge")

s = new_slide("Flagship: merigolix and the estrogen-threshold window", "Gynecology · GnRH antagonist", size=25)
bullets(s, [
    ("Mechanism", "oral GnRH receptor antagonist → LH/FSH ↓ → estradiol ↓ → lesion·pain 감소"),
    ("The trade-off", "E2 suppression이 부족하면 효과가 없고, 과하면 BMD loss와 hot flash"),
    ("The target", "partial E2 window (20–40 pg/mL) — estrogen threshold hypothesis"),
], 0.75, 1.9, 5.9, 3.0, size=16.5, gap=12)
rrect(s, 0.75, 4.75, 5.9, 1.55, PALETEAL, radius=0.08)
textbox(s, "Phase 1 PK/PD model (SAD/MAD, n = 51) 위에 HPG axis·menstrual cycle·lesion·BMD를 결합한 24-state QSP model",
        1.0, 4.75, 5.4, 1.55, size=15.5, color=NAVY, family=F_MED, anchor=MSO_ANCHOR.MIDDLE, spacing=1.2)
ix, iy, iw, ih = add_image(s, FIG / "ext/merigolix_app.png", 6.95, 1.95, 5.65, 4.0)
textbox(s, "Interactive QSP dashboard · pipetqsp.shinyapps.io/merigolix", ix, iy + ih + 0.1, iw, 0.3,
        size=11.5, color=GREY, align=PP_ALIGN.CENTER, italic=True, spacing=1.0)
source(s, "Bae, Kang et al., Clin Pharmacokinet (in revision) · github.com/pipetcpt/merigolix-qsp")

fig_slide("Merigolix: the trade-off becomes dose optimization", FIG / "s6_merigolix.png",
          "Flagship · the named example",
          cap=f"{S6['best_dose']} mg QD에서 composite therapeutic index 최고 ({S6['best_ti']}) — 120 mg은 tolerability 우선 대안",
          box=(0.75, 1.75, 11.85, 4.3), src="Merigolix QSP model, virtual population at week 12 (github.com/pipetcpt/merigolix-qsp)")

s = new_slide("Connecting to NAMs: gut–liver–kidney MPS + PBPK", "NAMs bridge", size=25)
for i, (num, lab) in enumerate([("6", "drugs"), ("24", "clinical observations"),
                                 ("0", "clinical PK used to fit"), ("24/24", "within 4-fold")]):
    cx = 0.75 + i * 3.0
    textbox(s, [[(num + "  ", {"family": F_BOLD, "size": 24, "color": TEAL}),
                 (lab, {"size": 14.5, "color": INK})]], cx, 1.78, 2.9, 0.5, spacing=1.0)
add_image(s, FIG / "glk_validation.png", 0.75, 2.3, 11.85, 3.75)
caption(s, "Chip-fitted parameter만으로 human AUC 예측 (GMFE 1.85×) — industry PBPK benchmark (2–5×) 범위",
        y=6.15, size=14.5)
source(s, "hiEC–HepaRG–RPTEC gut–liver–kidney MPS + whole-body PBPK (PIPET MPS)")

s = new_slide("Disease-state chips: gut- vs liver-driven change", "NAMs bridge")
for i, (drug, chip, num, where, detail, col) in enumerate([
        ("Gefitinib", "DSS colitis chip", "0.49×", "gut-localized", "Fg 변화, Fa / Fh 불변", TEAL),
        ("Simvastatin", "Cisplatin-injury chip", "2.19×", "hepatic-localized", "Kp,uu,liver 절반, OATP1B1 억제", RED)]):
    cx = 0.75 + i * 6.075
    rrect(s, cx, 1.85, 5.775, 2.3, PALE, radius=0.07)
    textbox(s, f"{drug} · {chip}", cx + 0.35, 2.05, 5.1, 0.45, size=17, color=NAVY, family=F_BOLD, spacing=1.0)
    textbox(s, [[(num, {"family": F_BOLD, "size": 40, "color": col}),
                 ("  AUC, " + where, {"size": 16, "color": INK})]], cx + 0.35, 2.6, 5.2, 0.8,
            anchor=MSO_ANCHOR.MIDDLE, spacing=1.0)
    textbox(s, detail, cx + 0.35, 3.5, 5.1, 0.45, size=15, color=GREY, spacing=1.0)
bullets(s, [
    ("Same pipeline", "disease chip으로 re-fit하되 영향받지 않는 organ의 parameter는 고정"),
    ("The payoff", "Phase 1에서 제외되는 IBD·AKI 환자를 위한 testable hypothesis — Study ④와 같은 문제의식"),
], 0.75, 4.5, 11.85, 1.3, size=17, gap=10)
callout(s, "새 case study와 MPS/NAM 연구는 하나의 program: mechanism이 in vitro·virtual signal을 human safety decision으로 translate",
        y=5.6, h=0.8, size=16)

s = new_slide("Prerequisites for trusting a chip-derived number", "Discipline for regulatory acceptance", size=25)
cards(s, [
    ("Chip-specific behavior", "plastic·membrane adsorption 보정", TEAL),
    ("Dimensional scaling", "organ별 cell·fluidic scaling factor를 정확히", NAVY),
    ("Held-out clinical validation", "fitting에 쓰지 않은 clinical data로 검증", GOLD),
    ("Falsifiable correction rules", "curve-fitting이 아닌 mechanism-based rule", RED),
], h=1.5, numbered=True, body_size=15.5)
callout(s, "Study ④의 DOAC validation도 같은 원칙 — fitted fudge factor가 아닌 mechanistic (P-gp) correction",
        y=5.6, h=0.8, size=16)

# =============================================================================
# PART 6 -- regulatory acceptance, people, conclusion
# =============================================================================
divider(6, "Prerequisites, people, and the path forward",
        "Context of use, clinical data 기반 validation — 그리고 이를 translate할 전문 인력")

s = new_slide("What regulatory acceptance will require", "Context of use")
table(s, [
    ["Principle", "What it asks for"],
    ["Context of use", "The precise question each model or chip is qualified to answer"],
    ["Biological relevance", "The human mechanism the NAM reproduces"],
    ["Technical characterization", "Reproducibility, dynamic range, failure modes"],
    ["Fit-for-purpose validation", "Agreement with held-out clinical data (the step QSP / PBPK supplies)"],
], 0.75, 1.85, 11.85, col_w=[3.3, 8.5], row_h=0.62, size=15)
callout(s, "QSP/PBPK는 이미 확립된 MIDD pathway — NAM/MPS data를 여기에 통합하는 것이 decision-grade safety evidence로 가는 가장 현실적인 길",
        y=5.25, h=1.05)

s = new_slide("The real bottleneck is people, not technology", "Unmet needs in our own field")
bullets(s, [
    ("Workforce shortage", "NAM signal을 PBPK/QSP로 해석할 clinical pharmacologist 부족"),
    ("No standard workflow", "FDA의 네 가지 NAM validation principle을 위한 표준 workflow 부재"),
    ("Few head-to-head validations", "chip·organoid model을 human clinical data와 직접 비교한 연구가 드묾"),
], 0.75, 2.1, 11.85, 2.3, size=19, gap=14)
rrect(s, 0.75, 4.75, 11.85, 1.1, NAVY, radius=0.1)
textbox(s, "Not a data problem, not a technology problem: a workforce problem.", 0.75, 4.75, 11.85, 1.1,
        size=21, color=WHITE, family=F_MED, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE, spacing=1.0)

s = new_slide("Conclusion", "What to take away")
for i, (lead, body, col) in enumerate([
        ("QSP + PBPK are the bridge", "in vitro·virtual signal을 human exposure, safety margin, benefit–risk로 translate", NAVY),
        ("Trade-offs become optimizations", "MS dosing interval · ICI dose ratio · incretin titration · DOAC individualization · merigolix window", TEAL),
        ("At scale", f"LLM이 구축한 {N_MODELS}개 open QSP model library — {PCT_SAFE}%가 이미 explicit safety axis 보유", BLUE),
        ("Validated honestly", "gut–liver–kidney MPS가 clinical PK 없이 human AUC를 4-fold 이내로 예측 — fudge factor가 아닌 mechanism", GOLD),
        ("The gap is people", "NAM-derived signal을 clinical safety decision으로 연결할 전문 인력이 필요", RED)]):
    cy = 1.8 + i * 0.93
    circ = s.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.8), Inches(cy + 0.1), Inches(0.58), Inches(0.58))
    circ.fill.solid(); circ.fill.fore_color.rgb = col; circ.line.fill.background(); circ.shadow.inherit = False
    textbox(s, str(i + 1), 0.8, cy + 0.1, 0.58, 0.58, size=18, color=WHITE, family=F_BOLD,
            align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE, spacing=1.0)
    textbox(s, lead, 1.65, cy + 0.02, 10.9, 0.4, size=18, color=col, family=F_BOLD, spacing=1.0)
    textbox(s, body, 1.65, cy + 0.42, 10.9, 0.45, size=15.5, color=INK, spacing=1.05)

s = prs.slides.add_slide(BLANK)
textbox(s, "Thank you", 0.95, 1.7, 8.0, 1.0, size=44, color=NAVY, family=F_BOLD, spacing=1.0)
textbox(s, "Questions and collaboration welcome", 0.98, 2.8, 8.0, 0.5, size=20, color=TEAL, family=F_MED,
        spacing=1.0)
textbox(s, [
    [("Sungpil Han, MD, PhD", {"family": F_BOLD, "size": 18})],
    [("Dept of Pharmacology, College of Medicine, The Catholic Univ. of Korea", {"size": 14, "color": GREY})],
    [("Dept of Clinical Pharmacology & Therapeutics, Seoul St. Mary’s Hospital", {"size": 14, "color": GREY})],
    [("shan@catholic.ac.kr", {"family": F_BOLD, "size": 15, "color": NAVY})],
], 0.98, 3.75, 7.6, 1.9, spacing=1.35)
rrect(s, 8.65, 3.8, 3.95, 1.35, PALE, radius=0.1)
textbox(s, [[("Fully open", {"family": F_BOLD, "size": 14, "color": GOLD})],
            [("github.com/pipetcpt/qsp", {"family": F_MED, "size": 16, "color": NAVY})]],
        8.9, 3.8, 3.6, 1.35, anchor=MSO_ANCHOR.MIDDLE, spacing=1.3)

# =============================================================================
# enforce latinLnBrk="0" on every paragraph we did not write ourselves
# =============================================================================
for sl in prs.slides:
    for p in sl._element.iter(qn("a:p")):
        pPr = p.find(qn("a:pPr"))
        if pPr is None:
            pPr = etree.Element(qn("a:pPr"))
            p.insert(0, pPr)
        pPr.set("latinLnBrk", "0")

prs.save(str(OUT))


def rewrite_line_breaks(path):
    """Second pass over every XML part: no latinLnBrk="1" may survive."""
    tmp = Path(tempfile.mkstemp(suffix=".pptx")[1])
    n_fixed = 0
    with zipfile.ZipFile(path) as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename.endswith(".xml"):
                txt = data.decode("utf-8")
                n_fixed += txt.count('latinLnBrk="1"')
                data = txt.replace('latinLnBrk="1"', 'latinLnBrk="0"').encode("utf-8")
            zout.writestr(item, data)
    shutil.move(str(tmp), str(path))
    Path(path).chmod(0o644)
    return n_fixed


fixed = rewrite_line_breaks(OUT)
print(f"saved {OUT.name}: {len(prs.slides._sldIdLst)} slides; rewrote {fixed} inherited latinLnBrk=\"1\"")
