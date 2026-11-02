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

Animation and notes
  Every slide except the cover has one click: the shapes added after reveal(s)
  (the conclusion caption/callout, or the second half of the content) fade in
  together. Every slide carries a speaker script (>= 5 sentences) in its notes;
  "[클릭]" in the script marks where to click. The build fails if a script is
  shorter than 5 sentences or if a click marker and an animation disagree.
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
    tb = textbox(s, text, 4.9, 7.02, 7.75, 0.3, size=9.5, color=GREY, align=PP_ALIGN.RIGHT, italic=True,
                 spacing=1.0)
    tb.name = "Source"          # never animated: citations stay visible
    return tb


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
          numbered=False, start=0):
    cw = (total_w - gap * (cols - 1)) / cols
    for i, (head, body, col) in enumerate(items):
        cx = x0 + (i % cols) * (cw + gap)
        cy = y0 + (i // cols) * (h + gap)
        rrect(s, cx, cy, cw, h, PALE, radius=0.07)
        tx = cx + 0.3
        if numbered:
            textbox(s, str(start + i + 1), cx + 0.28, cy + 0.2, 0.6, 0.7, size=32, color=col, family=F_BOLD,
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
    reveal(s)
    textbox(s, subtitle, 0.97, 4.9, 11.4, 1.0, size=18, color=GREY, spacing=1.2)
    return s


def fig_slide(title, fig, eyebrow_txt=None, cap=None, box=(0.75, 1.75, 11.85, 4.3), src=None):
    s = new_slide(title, eyebrow_txt)
    add_image(s, fig, *box)
    if cap:
        reveal(s)
        caption(s, cap)
    if src:
        source(s, src)
    return s


def intro_slide(title, eyebrow_txt, items, layer, src):
    s = new_slide(title, eyebrow_txt)
    bullets(s, items, 0.75, 1.95, 11.85, 2.5, size=18, gap=12)
    reveal(s)
    rrect(s, 0.75, 4.1, 11.85, 1.6, PALETEAL, radius=0.08)
    textbox(s, [[("QSP layer   ", {"family": F_BOLD, "color": TEAL}), (layer, {})]],
            1.05, 4.1, 11.25, 1.6, size=17, color=NAVY, anchor=MSO_ANCHOR.MIDDLE, spacing=1.2)
    source(s, src)
    return s


# =============================================================================
# one-click reveal + speaker notes
# =============================================================================
REVEAL = {}          # slide element id -> (slide, index of first shape to reveal)
NOTES = {}           # slide element id -> (slide, script)


def reveal(s):
    """Everything added to slide s after this call fades in on the first click."""
    REVEAL[id(s._element)] = (s, len(s.shapes))


def say(s, script):
    """Speaker script for slide s; '[클릭]' marks where the one click happens."""
    NOTES[id(s._element)] = (s, script)


def _timing_xml(targets):
    """PowerPoint main sequence: one click, all targets fade in together (0.5 s)."""
    ctn = iter(range(5, 100000))
    effects = []
    for k, (spid, _) in enumerate(targets):
        a, b, c = next(ctn), next(ctn), next(ctn)
        node = "clickEffect" if k == 0 else "withEffect"
        effects.append(
            f'<p:par><p:cTn id="{a}" presetID="10" presetClass="entr" presetSubtype="0" fill="hold" '
            f'grpId="0" nodeType="{node}"><p:stCondLst><p:cond delay="0"/></p:stCondLst><p:childTnLst>'
            f'<p:set><p:cBhvr><p:cTn id="{b}" dur="1" fill="hold"><p:stCondLst><p:cond delay="0"/>'
            f'</p:stCondLst></p:cTn><p:tgtEl><p:spTgt spid="{spid}"/></p:tgtEl><p:attrNameLst>'
            f'<p:attrName>style.visibility</p:attrName></p:attrNameLst></p:cBhvr><p:to>'
            f'<p:strVal val="visible"/></p:to></p:set>'
            f'<p:animEffect transition="in" filter="fade"><p:cBhvr><p:cTn id="{c}" dur="500"/>'
            f'<p:tgtEl><p:spTgt spid="{spid}"/></p:tgtEl></p:cBhvr></p:animEffect>'
            f'</p:childTnLst></p:cTn></p:par>')
    builds = []
    for spid, kind in targets:
        if kind == "textbox":
            builds.append(f'<p:bldP spid="{spid}" grpId="0"/>')
        elif kind == "shape":
            builds.append(f'<p:bldP spid="{spid}" grpId="0" animBg="1"/>')
        elif kind == "table":
            builds.append(f'<p:bldGraphic spid="{spid}" grpId="0"><p:bldAsOne/></p:bldGraphic>')
    bld = f'<p:bldLst>{"".join(builds)}</p:bldLst>' if builds else ""
    return (
        '<p:timing xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main" '
        'xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"><p:tnLst><p:par>'
        '<p:cTn id="1" dur="indefinite" restart="never" nodeType="tmRoot"><p:childTnLst>'
        '<p:seq concurrent="1" nextAc="seek"><p:cTn id="2" dur="indefinite" nodeType="mainSeq">'
        '<p:childTnLst><p:par><p:cTn id="3" fill="hold"><p:stCondLst><p:cond delay="indefinite"/>'
        '</p:stCondLst><p:childTnLst><p:par><p:cTn id="4" fill="hold"><p:stCondLst>'
        '<p:cond delay="0"/></p:stCondLst><p:childTnLst>' + "".join(effects) +
        '</p:childTnLst></p:cTn></p:par></p:childTnLst></p:cTn></p:par></p:childTnLst></p:cTn>'
        '<p:prevCondLst><p:cond evt="onPrev" delay="0"><p:tgtEl><p:sldTgt/></p:tgtEl></p:cond>'
        '</p:prevCondLst><p:nextCondLst><p:cond evt="onNext" delay="0"><p:tgtEl><p:sldTgt/>'
        '</p:tgtEl></p:cond></p:nextCondLst></p:seq></p:childTnLst></p:cTn></p:par></p:tnLst>'
        + bld + '</p:timing>')


def apply_reveal(s, start):
    targets = []
    for sh in list(s.shapes)[start:]:
        if sh.name == "Source":
            continue
        el = sh._element
        if el.tag == qn("p:pic"):
            kind = "picture"
        elif el.tag == qn("p:graphicFrame"):
            kind = "table"
        else:
            cnv = el.find(".//" + qn("p:cNvSpPr"))
            kind = "textbox" if (cnv is not None and cnv.get("txBox") == "1") else "shape"
        targets.append((sh.shape_id, kind))
    if not targets:
        return 0
    sld = s._element
    for old in sld.findall(qn("p:timing")):
        sld.remove(old)
    timing = etree.fromstring(_timing_xml(targets))
    anchor = sld.find(qn("p:clrMapOvr"))
    if anchor is None:
        anchor = sld.find(qn("p:cSld"))
    anchor.addnext(timing)
    return len(targets)


def write_notes(s, script):
    tf = s.notes_slide.notes_text_frame
    parts = [x.strip() for x in script.split("[클릭]")]
    lines = [parts[0]] + [y for part in parts[1:] for y in ("[클릭]", part)]
    tf.text = "\n".join(l for l in lines if l)
    for p in tf.paragraphs:
        no_word_break(p)
        for r in p.runs:
            r.font.size = Pt(14)
            rPr = r._r.get_or_add_rPr(); rPr.set("lang", "ko-KR"); rPr.set("altLang", "en-US")


def n_sentences(text):
    return len(re.findall(r"[.?!](?=\s|$)", text.replace("[클릭]", " ")))


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
say(cover, "안녕하세요, 가톨릭대학교 의과대학 약리학교실 한성필입니다. "
    "오늘은 'Systems-Based Quantitative Modeling for Clinical Safety Assessment', 즉 QSP가 임상 안전성 평가에서 "
    "어떤 역할을 할 수 있는지 말씀드리겠습니다. "
    "먼저 문제 제기와 규제 환경의 변화, 이어서 QSP와 PBPK의 기본 개념을 짚겠습니다. "
    "그다음 저희가 LLM으로 구축하고 있는 QSP model library와 다섯 가지 safety case study를 보여드리겠습니다. "
    "마지막으로 NAM 연계 사례와 앞으로의 과제로 마무리하겠습니다. "
    "오른쪽 위 QR code로 발표자료를 내려받으실 수 있습니다.")

s = new_slide("Drug development fails late, and toxicity leads", "Why we are here")
stat_cards(s, [("~90%", "임상시험에 진입한 후보물질 중 승인에 실패하는 비율", RED),
               ("10–15 yr", "신약 하나당 개발 기간", NAVY),
               ("$1–2.6 B", "신약 하나당 개발 비용", NAVY)])
reveal(s)
bullets(s, [
    ("Unexpected toxicity", "late-stage attrition과 post-marketing withdrawal의 주요 원인"),
    ("Healthy-volunteer Phase 1", "젊고 건강한 사람의 exposure만 확인 — 실제 약을 복용할 환자와 다름"),
    ("The costliest failures", "취약한 실제 환자에서 뒤늦게 드러나는 safety signal"),
], 0.75, 4.15, 11.85, 2.4, size=18, gap=12)
source(s, "DiMasi J Health Econ 2016 · Kola & Landis Nat Rev Drug Discov 2004")
say(s, "신약 개발의 현실부터 보겠습니다. "
    "임상시험에 들어간 후보물질의 약 90%는 결국 승인을 받지 못합니다. "
    "신약 하나를 개발하는 데 10년에서 15년, 비용은 10억에서 26억 달러가 듭니다. "
    "[클릭] 실패 원인 가운데 큰 비중을 차지하는 것이 예상하지 못한 toxicity이고, "
    "이는 late-stage attrition뿐 아니라 시판 후 withdrawal로도 이어집니다. "
    "문제는 Phase 1이 젊고 건강한 volunteer에서 exposure를 확인한다는 점입니다. "
    "실제로 약을 복용하는 환자는 고령이거나 장기 기능이 떨어진 경우가 많습니다. "
    "그래서 가장 비싼 실패는 이런 취약한 환자에서 뒤늦게 드러나는 safety signal입니다.")

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
reveal(s)
callout(s, "결과적으로 in vitro signal과 animal data가 실제 위험 환자의 human exposure·response로 잘 translate되지 않음",
        y=5.0, h=1.05)
source(s, "Van Norman, JACC Basic Transl Sci 2019")
say(s, "기존 safety assessment 도구를 하나씩 보면 이 간극이 왜 생기는지 알 수 있습니다. "
    "2D cell culture는 tissue architecture와 organ-level function을 재현하지 못합니다. "
    "Animal model은 drug disposition과 response에서 종간 차이가 여전히 큽니다. "
    "Phase 1 trial은 healthy volunteer 중심이라 간·신장 기능 저하, 임신, 고령, 희귀질환 환자가 빠져 있습니다. "
    "[클릭] 결국 in vitro signal과 animal data가 실제 위험군 환자의 human exposure와 response로 잘 이어지지 않는다는 것이 핵심 문제입니다. "
    "오늘 발표는 이 translation gap을 mechanistic model로 메우는 방법에 관한 것입니다.")

s = new_slide("The regulatory shift: from animal mandate to NAMs", "What changed")
table(s, [
    ["Year", "Milestone", "What it means"],
    ["1938", "Federal Food, Drug, and Cosmetic Act", "Animal testing de facto required before human use"],
    ["2022", "FDA Modernization Act 2.0", "“Animal tests” replaced by “nonclinical tests”: in vitro and in silico allowed"],
    ["2025", "FDA Roadmap to Reducing Animal Testing", "Phased rollout, starting with monoclonal antibodies"],
    ["Now", "New approach methodologies (NAMs)", "MPS / organ-chips, organoids, in silico models"],
], 0.75, 1.85, 11.85, col_w=[1.0, 4.0, 6.0], row_h=0.62, size=14.5)
reveal(s)
callout(s, "다만 in vitro signal이 decision-grade evidence가 되려면, mechanistic model이 이를 human exposure와 response로 translate해야 함",
        y=5.25, h=1.05)
source(s, "FDA Modernization Act 2.0, Pub. L. 117-328 (2022) · FDA Roadmap (2025)")
say(s, "이 문제에 대해 규제 환경도 빠르게 바뀌고 있습니다. "
    "1938년 FD&C Act 이후 80년 넘게 인체 투여 전 동물시험이 사실상 의무였습니다. "
    "2022년 FDA Modernization Act 2.0은 법 조문의 'animal tests'를 'nonclinical tests'로 바꿔 in vitro와 in silico 방법을 공식적으로 허용했습니다. "
    "2025년에는 FDA가 Roadmap to Reducing Animal Testing을 발표해 monoclonal antibody부터 단계적으로 동물시험을 줄이겠다고 밝혔습니다. "
    "그 대안이 MPS·organ-chip, organoid, in silico model을 아우르는 NAM, 즉 new approach methodologies입니다. "
    "[클릭] 다만 in vitro signal 자체가 곧바로 규제 결정의 근거가 되지는 않습니다. "
    "Mechanistic model이 그 signal을 human exposure와 response로 translate해 줄 때 비로소 decision-grade evidence가 됩니다.")

s = fig_slide("NAMs need a translational bridge", FIG / "ext/nam_loop.png", "The NAM translation loop",
              cap="PBPK와 QSP가 NAM data를 qualify하고 human으로 translate하는 핵심 단계 — pattern-matching만으로는 대체 불가",
              box=(1.4, 1.85, 10.55, 4.2), src="Cao Y, Clin Pharmacol Ther 2025")
say(s, "이 그림은 Cao 등이 2025년 CPT에 제시한 NAM translation loop입니다. "
    "왼쪽은 organoid, organ-on-chip, 2D·3D culture, in silico로 구성된 NAM이고, context-of-use에 따라 선택됩니다. "
    "오른쪽은 biological data와 clinical data가 쌓인 dataset입니다. "
    "가운데가 둘 사이를 오가며 NAM을 qualify하고 human으로 translate하는 단계이며, 여기에 AI/ML 도구와 함께 PBPK·QSP 같은 mechanistic model이 들어갑니다. "
    "[클릭] 저는 바로 이 연결 고리가 임상약리학이 기여할 지점이라고 봅니다. "
    "단순한 pattern-matching만으로는 이 translation을 대신할 수 없습니다.")

# =============================================================================
# PART 2 -- what QSP / PBPK is
# =============================================================================
s = divider(2, "QSP & PBPK: the mechanistic bridge",
            "Drug exposure를 target, pathway, organ, patient outcome까지 수학적으로 연결")
say(s, "이제 QSP와 PBPK가 무엇인지 간단히 정리하겠습니다. "
    "[클릭] 핵심은 drug exposure를 target, pathway, organ을 거쳐 patient outcome까지 수학적으로 연결한다는 점입니다. "
    "이 분야를 잘 아시는 분도 계시겠지만, 뒤의 case study를 이해하는 데 필요한 개념만 짚고 넘어가겠습니다. "
    "특히 dose라는 handle을 돌리면 결과가 어떻게 바뀌는가라는 관점에 주목해 주십시오. "
    "이 관점이 오늘 발표 전체를 관통합니다.")

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
reveal(s)
callout(s, "둘을 결합하면 first-in-human dose, safety margin, DDI risk, special population response, chip-derived parameter의 IVIVE까지 예측 가능",
        y=5.25, h=1.05)
say(s, "QSP와 PBPK는 drug에서 target, pathway, disease를 거쳐 patient에 이르는 인과 사슬을 수식으로 표현합니다. "
    "PBPK는 약물이 우리 몸 어디에 어떤 농도로 분포하는지를 다룹니다. "
    "QSP는 그 농도에서 질환과 safety가 어떻게 변하는지를 다룹니다. "
    "즉 PBPK가 exposure를, QSP가 response를 담당한다고 보시면 됩니다. "
    "[클릭] 둘을 결합하면 first-in-human dose, safety margin, DDI risk, special population response를 예측할 수 있습니다. "
    "나아가 chip에서 얻은 parameter를 사람에게 외삽하는 IVIVE에도 같은 틀을 씁니다.")

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
reveal(s)
callout(s, f"Library의 각 질환 모델은 {LIB['min_odes']}–{LIB['max_odes']}개의 coupled ODE (median {LIB['median_odes']})와 수십–수백 개 parameter로 구성",
        y=5.35, h=0.7, size=16)
say(s, "QSP model을 열어 보면 결국 coupled nonlinear ODE system입니다. "
    "State x는 약물 농도, receptor occupancy, signaling 물질, cell population, biomarker, clinical endpoint 같은 변수입니다. "
    "Parameter θ는 rate constant, EC50과 Emax, Hill coefficient, clearance, binding affinity 같은 값입니다. "
    "Input u(t)는 dosing regimen이고, 치료에서 우리가 직접 조절할 수 있는 유일한 handle입니다. "
    f"[클릭] 저희 library의 질환 모델은 질환당 {LIB['min_odes']}개에서 {LIB['max_odes']}개, 중앙값 {LIB['median_odes']}개의 coupled ODE로 구성되어 있습니다. "
    "이 정도 규모라도 구조가 명확하면 dosing을 바꿨을 때 어떤 일이 생기는지 체계적으로 탐색할 수 있습니다.")

s = new_slide("The same nonlinear motifs recur across diseases", "A shared grammar")
table(s, [
    ["Motif", "Governing idea", "Typical use"],
    ["Hill / Emax", "E(C) = Emax·Cⁿ / (EC₅₀ⁿ + Cⁿ)", "Any drug dose–response"],
    ["Indirect response", "dR/dt = k_in·f(·) − k_out·R", "Biomarker turnover"],
    ["Mass action / TMDD", "Drug–target binding, target-mediated clearance", "mAbs, saturable drugs"],
    ["Logistic growth", "dN/dt = rN(1 − N/K) − kill(C)·N", "Tumor burden under therapy"],
    ["Power-law feedback", "e.g. EPO ∝ (Hb₀/Hb)^γ", "Homeostatic loops"],
], 0.75, 1.85, 11.85, col_w=[2.3, 4.6, 3.0], row_h=0.56, size=15)
reveal(s)
callout(s, f"{N_MODELS}개 모델 전체, 그리고 모든 chip-to-human bridge가 같은 grammar를 공유", y=5.6, h=0.8)
say(s, "질환은 달라도 모델 안에서 반복되는 nonlinear motif는 의외로 몇 가지 되지 않습니다. "
    "Hill 또는 Emax 식은 거의 모든 dose–response에 쓰입니다. "
    "Indirect response model은 biomarker turnover를, mass action과 TMDD는 항체처럼 표적에 결합해 제거되는 약물을 설명합니다. "
    "Logistic growth는 치료 중 tumor burden을, power-law feedback은 EPO와 hemoglobin 같은 항상성 loop를 표현합니다. "
    f"[클릭] {N_MODELS}개 모델 전체와 chip-to-human bridge가 모두 이 같은 문법을 공유합니다. "
    "그래서 LLM이 모델을 만들 때도 이 motif를 template으로 쓰면 실수의 여지가 크게 줄어듭니다.")

s = new_slide("Why it matters: model-informed drug development", "MIDD, endorsed by FDA and EMA")
table(s, [
    ["Stage", "What QSP / PBPK contributes", "Safety payoff"],
    ["Target validation", "Does target modulation reach the phenotype?", "Drops implausible targets early"],
    ["First-in-human dose", "Exposure–response and safety margins from mechanism", "Safer, leaner Phase 1"],
    ["Trial design", "Virtual patients for enrichment and endpoints", "Fewer patients exposed to risk"],
    ["Special populations", "Hepatic / renal impairment, age, pregnancy", "Doses the excluded patient"],
    ["Combinations", "Synergy and antagonism in silico", "Rational, safer polypharmacy"],
], 0.75, 1.85, 11.85, col_w=[2.4, 4.6, 3.2], row_h=0.56, size=14.5)
reveal(s)
callout(s, "QSP는 efficacy–safety trade-off를 tunable optimization problem으로 바꿈 — 오늘 case study를 관통하는 주제",
        y=5.6, h=0.8)
say(s, "이런 모델이 실제 신약 개발 어디에 쓰이는지 MIDD 관점에서 정리했습니다. "
    "Target validation 단계에서는 표적을 조절했을 때 실제 phenotype까지 도달하는지 확인해 가능성이 낮은 target을 일찍 걸러냅니다. "
    "First-in-human dose 단계에서는 mechanism 기반의 exposure–response와 safety margin으로 Phase 1을 더 안전하게 설계합니다. "
    "Trial design에서는 virtual patient로 enrichment와 endpoint를 정해 위험에 노출되는 환자 수를 줄입니다. "
    "Special population과 combination therapy도 모델로 미리 따져볼 수 있습니다. "
    "[클릭] 정리하면 QSP는 efficacy–safety trade-off를 조절 가능한 optimization 문제로 바꿉니다. "
    "이것이 오늘 case study 전체를 관통하는 주제입니다.")

# =============================================================================
# PART 3 -- the living QSP library
# =============================================================================
s = divider(3, "A living QSP model library",
            "LLM-augmented workflow가 literature-grounded, quality-gated disease model을 하루 하나씩 구축")
say(s, "다음은 저희가 구축하고 있는 living QSP model library입니다. "
    "[클릭] LLM-augmented workflow가 문헌에 근거하고 quality gate를 통과한 질환 모델을 하루에 하나씩 만들어 갑니다. "
    "이 library가 왜 필요한지, 어떻게 만들어지는지 순서대로 말씀드리겠습니다. "
    "그리고 LLM이 만든 결과물을 믿을 수 있는 근거가 무엇인지도 함께 보여드리겠습니다. "
    "이 부분은 뒤에 나올 다섯 case study의 기반이기도 합니다.")

s = new_slide("The bottleneck: one model takes expert-months", "Why there were so few")
bullets(s, [
    ("Deep synthesis", "한 질환의 pathophysiology와 pharmacology를 깊이 통합해야 함"),
    ("Careful formulation", "ODE structure, parameterization, implementation, calibration"),
    ("The cost", "질환당 expert-months to years — 소수의 질환만 QSP model을 가짐"),
    ("The consequence", "모델이 silo화·비공개·재현 불가; long-tail disease는 모델링되지 않음"),
], 0.75, 2.1, 11.85, 2.7, size=19, gap=14)
reveal(s)
rrect(s, 0.75, 4.9, 11.85, 1.0, NAVY, radius=0.1)
textbox(s, "One disease ≈ one expert-year   ·   scarce, slow, closed", 0.75, 4.9, 11.85, 1.0, size=21,
        color=WHITE, family=F_MED, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE, spacing=1.0)
say(s, "그동안 QSP model이 많지 않았던 이유는 만드는 데 드는 비용 때문입니다. "
    "한 질환의 pathophysiology와 pharmacology를 깊이 통합해야 합니다. "
    "ODE 구조를 정하고, parameter를 정하고, 코드로 구현하고, calibration까지 해야 합니다. "
    "그래서 질환 하나에 expert-months에서 years가 걸리고, 결국 소수 질환만 모델을 갖게 됩니다. "
    "만들어진 모델도 silo에 갇혀 공개되지 않거나 재현되지 않는 경우가 많습니다. "
    "[클릭] 한 문장으로 요약하면 one disease, one expert-year입니다. "
    "희소하고, 느리고, 닫혀 있습니다.")

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
reveal(s)
ix, iy, iw, ih = add_image(s, FIG / "ext/map_igan.png", 7.25, 2.95, 5.35, 3.2)
textbox(s, "Example: IgA nephropathy mechanistic map", ix, iy + ih + 0.08, iw, 0.3, size=11.5,
        color=GREY, align=PP_ALIGN.CENTER, italic=True, spacing=1.0)
source(s, "github.com/pipetcpt/qsp (fully open)")
say(s, "저희의 접근은 autonomous LLM coding agent가 모델을 직접 만들게 하는 것입니다. "
    "Agent는 하루에 한 질환씩 네 가지 산출물을 만듭니다. "
    "100개 이상 node의 mechanistic map, PK와 disease PD를 결합한 15개 이상 state의 mrgsolve ODE model, "
    "6개 이상 tab의 Shiny dashboard, 그리고 30개 이상의 PubMed reference와 calibration note입니다. "
    "모든 산출물은 git으로 version control되어 누구나 확인할 수 있습니다. "
    "[클릭] 오른쪽은 실제로 만들어진 IgA nephropathy mechanistic map의 예시입니다. "
    "전체 library는 github.com/pipetcpt/qsp에 완전히 공개되어 있습니다.")

s = fig_slide("The daily build loop", FIG / "architecture.png", "One disease, end-to-end, every run",
              cap="Select → research → map → ODE model → dashboard → references → quality gate → commit (git stop-hook로 강제)",
              box=(0.9, 1.75, 11.55, 4.3))
say(s, "하루 한 번 돌아가는 build loop를 그림으로 보면 이렇습니다. "
    "질환 목록과 PubMed 문헌을 입력으로, Claude Code Routine이라는 LLM agent가 작업을 시작합니다. "
    "Agent는 mechanistic map, ODE model, Shiny app, reference 네 가지를 만들고, quality gate를 통과해야만 git에 commit하고 push합니다. "
    "Stop-hook이 걸려 있어 모든 변경이 commit과 push되지 않으면 작업이 끝나지 않습니다. "
    "다음 날에는 다음 질환으로 같은 과정이 반복됩니다. "
    "[클릭] 순서로 보면 select, research, map, ODE model, dashboard, references, quality gate, commit이고, 이 순서를 사람이 아니라 hook이 강제합니다.")

s = new_slide("Guardrails for LLM-written models", "Structure + grounding + verification")
cards(s, [
    ("Quality gates", "Map ≥100 nodes / ≥8 clusters · ODE ≥15 compartments / ≥5 scenarios · ≥30 references", TEAL),
    ("Structural templates", "Hill, turnover, TMDD 같은 reusable motif로 failure space를 줄임", NAVY),
], h=1.95)
reveal(s)
cards(s, [
    ("Literature grounding", "모든 parameter set에 calibration trial을 인용 — 핵심 anti-hallucination 장치", GOLD),
    ("Reproducibility", "Diffable·auditable·re-runnable code in git; 사용 전 human review", GREEN),
], h=1.95, y0=1.85 + 1.95 + 0.3)
say(s, "LLM이 쓴 코드를 어떻게 믿을 수 있느냐는 질문을 가장 많이 받습니다. "
    "첫째는 quality gate입니다. "
    "Map은 node 100개와 cluster 8개 이상, ODE는 compartment 15개와 scenario 5개 이상, reference는 30개 이상이어야 합니다. "
    "둘째는 structural template으로, 앞서 본 Hill, turnover, TMDD 같은 motif를 재사용해 오류가 생길 공간 자체를 줄입니다. "
    "[클릭] 셋째는 literature grounding입니다. "
    "모든 parameter set에 calibration trial을 인용하게 해 hallucination을 막는 핵심 장치로 씁니다. "
    "넷째는 reproducibility로, 모든 코드는 diff와 audit, 재실행이 가능하고 실제로 쓰기 전에는 반드시 사람이 review합니다.")

s = new_slide(f"Coverage: {LIB['n_areas']} therapeutic areas, one safety atlas", "The library as a safety atlas")
add_image(s, FIG / "s1_areas.png", 4.15, 1.8, 8.45, 4.3)
caption(s, "Safety axis: model code·map·README·references에 서로 다른 toxicology term(toxicity, QT, DILI, PML 등)이 3개 이상 명시된 모델",
        y=6.22, size=13)
reveal(s)
for i, (num, lab, col) in enumerate([
        (str(N_MODELS), "disease QSP models", NAVY),
        (str(LIB["n_areas"]), "therapeutic areas", NAVY),
        (f"{PCT_SAFE}%", "explicit safety / toxicity axis", TEAL)]):
    cy = 1.85 + i * 1.47
    rrect(s, 0.75, cy, 3.1, 1.3, PALE, radius=0.08)
    textbox(s, num, 0.75, cy + 0.1, 3.1, 0.7, size=34, color=col, family=F_BOLD, align=PP_ALIGN.CENTER,
            anchor=MSO_ANCHOR.MIDDLE, spacing=1.0)
    textbox(s, lab, 0.9, cy + 0.82, 2.8, 0.4, size=13, color=INK, align=PP_ALIGN.CENTER, spacing=1.0)
say(s, "현재 library가 어디까지 왔는지 보겠습니다. "
    f"그래프는 {LIB['n_areas']}개 therapeutic area별 모델 수이고, 진한 청록색이 explicit safety 또는 toxicity axis를 가진 모델입니다. "
    "Safety axis는 model code, map, README, reference에 toxicity, QT, DILI, PML 같은 서로 다른 독성 용어가 3개 이상 명시된 경우로 정의했습니다. "
    "Oncology, neurology, endocrine과 metabolic 순으로 모델이 많습니다. "
    f"[클릭] 숫자로 정리하면 {N_MODELS}개 질환 모델, {LIB['n_areas']}개 therapeutic area이고, 그중 {PCT_SAFE}%가 이미 safety axis를 갖고 있습니다. "
    "즉 이 library는 efficacy만 보는 모델 모음이 아니라 safety atlas로도 쓸 수 있습니다.")

s = fig_slide("Quality gates hold across the library", FIG / "s1_quality.png", "Size & grounding",
              cap=f"{LIB['n_odes_gate']}/{LIB['n_mined']} 모델이 ≥15 ODE gate (median {LIB['median_odes']}), {LIB['n_refs_gate']}/{LIB['n_mined']}가 ≥30 reference gate 통과 (median {LIB['median_refs']})",
              box=(0.9, 1.75, 11.55, 4.3))
say(s, "Quality gate가 실제로 지켜지는지 분포로 확인했습니다. "
    f"왼쪽은 모델당 ODE compartment 수로, 중앙값은 {LIB['median_odes']}개이고 빨간 점선이 15개 gate입니다. "
    f"오른쪽은 모델당 PubMed reference 수로, 중앙값은 {LIB['median_refs']}개이고 gate는 30개입니다. "
    "Gate 아래에 있는 모델도 일부 있다는 점을 숨기지 않고 그대로 보여드립니다. "
    f"[클릭] {LIB['n_mined']}개 중 {LIB['n_odes_gate']}개가 ODE gate를, {LIB['n_refs_gate']}개가 reference gate를 통과했습니다. "
    "기준에 못 미치는 모델이 숫자로 바로 드러나기 때문에 어디를 보완해야 할지도 분명합니다.")

s = fig_slide("One consistent schema across the library", FIG / "s1_atlas.png", "One framework",
              cap="모든 모델이 같은 4-artifact schema — ODE 수 × map 크기 × reference 수로 한눈에 비교 가능",
              box=(1.4, 1.75, 10.55, 4.3))
say(s, "이 그림은 library 전체를 한 장에 놓은 atlas입니다. "
    "점 하나가 질환 모델 하나이고, x축은 ODE compartment 수, y축은 mechanistic map의 node 수, 점의 크기는 reference 수입니다. "
    "색은 모델이 많은 상위 7개 therapeutic area를 나타냅니다. "
    "질환마다 모델 크기는 다르지만 모두 같은 네 가지 산출물 구조를 따릅니다. "
    "[클릭] 같은 schema를 공유하기 때문에 질환 사이를 한눈에 비교할 수 있습니다. "
    "한 질환에서 만든 safety module을 다른 질환에 옮겨 쓰기도 쉬워집니다.")

# =============================================================================
# PART 4 -- five case studies
# =============================================================================
s = divider(4, "Five safety case studies",
            "Autoimmune · oncology · obesity · aging · gynecology — efficacy–safety trade-off를 optimization으로")
say(s, "이제 이 library를 바탕으로 한 다섯 가지 safety case study를 보여드리겠습니다. "
    "[클릭] Autoimmune, oncology, obesity, aging, gynecology에 걸쳐 있고, 각각 efficacy–safety trade-off를 optimization 문제로 다시 씁니다. "
    "각 study에서 tunable handle이 무엇인지, safety axis가 무엇인지, 모델이 무엇을 계산하는지를 같은 순서로 보시면 됩니다. "
    "수치는 모두 공개된 임상시험 결과에 calibration했습니다. "
    "모델 코드와 그림은 모두 다시 실행해 재현할 수 있습니다.")

s = new_slide("Five case studies, one pattern", "Roadmap")
table(s, [
    ["Study", "Disease area", "Tunable handle", "Safety axis"],
    ["① Multiple sclerosis", "Autoimmune / neuro", "Natalizumab dosing interval", "PML risk"],
    ["② Checkpoint inhibitors", "Oncology (melanoma)", "Nivolumab : ipilimumab ratio", "Grade 3–4 irAE"],
    ["③ Incretin therapy", "Obesity + aging", "Target dose and escalation", "GI AE · lean-mass loss"],
    ["④ DOACs in the elderly", "Aging / organ function", "Age- and renal-based dose", "Bleeding"],
    ["★ Merigolix", "Gynecology", "GnRH antagonist dose", "Bone loss · hot flashes"],
], 0.75, 1.85, 11.85, col_w=[2.9, 2.6, 3.4, 2.6], row_h=0.56, size=14.5)
reveal(s)
callout(s, "공통 질문: dose·schedule이라는 handle을 어떻게 조절해야 efficacy를 지키면서 safety risk를 줄일 수 있는가?",
        y=5.6, h=0.8)
say(s, "다섯 study를 한 표로 정리했습니다. "
    "Multiple sclerosis에서는 natalizumab dosing interval과 PML risk를 봅니다. "
    "Checkpoint inhibitor에서는 nivolumab과 ipilimumab의 dose 비율과 grade 3–4 irAE를 봅니다. "
    "Incretin therapy에서는 target dose와 escalation 속도에 따른 GI AE와 lean-mass loss를, DOAC에서는 나이와 신기능에 따른 dose와 bleeding을 봅니다. "
    "마지막 merigolix는 GnRH antagonist dose와 bone loss, hot flash 사이의 균형입니다. "
    "[클릭] 다섯 study의 공통 질문은 하나입니다. "
    "Dose와 schedule이라는 handle을 어떻게 조절해야 efficacy를 지키면서 safety risk를 줄일 수 있는가입니다.")

# ---- Study 1: MS -------------------------------------------------------------
s = intro_slide("Study ①  Multiple sclerosis: natalizumab vs PML", "Autoimmune / neurology", [
    ("Natalizumab (anti-α4-integrin)", "relapse 억제 효과가 가장 강력한 DMT 중 하나 — 그러나 CNS immune surveillance를 낮춤"),
    ("The risk", "JC virus reactivation에 의한 progressive multifocal leukoencephalopathy (PML)"),
    ("The question", "dosing interval 조절로 efficacy는 유지하고 PML risk만 줄일 수 있는가?"),
], "dosing interval 동안의 α4 receptor saturation → relapse-rate suppression, 그리고 stratified PML hazard (anti-JCV index · prior immunosuppressant · duration)에 extended-interval hazard ratio 적용",
    "Polman NEJM 2006 · Ho Lancet Neurol 2017 · Foley Lancet Neurol 2022 · Ryerson Neurology 2019")
say(s, "첫 번째는 multiple sclerosis입니다. "
    "Natalizumab은 anti-α4-integrin 항체로, relapse 억제 효과가 가장 강한 DMT 중 하나입니다. "
    "하지만 림프구의 CNS 유입을 막아 immune surveillance를 떨어뜨립니다. "
    "그 결과 JC virus가 재활성화되어 progressive multifocal leukoencephalopathy, 즉 PML이 생길 수 있습니다. "
    "질문은 dosing interval만 조절해서 efficacy는 유지하고 PML risk만 줄일 수 있느냐입니다. "
    "[클릭] 모델은 dosing interval 동안 α4 receptor saturation이 어떻게 변하는지 계산하고, 이를 relapse-rate suppression으로 연결합니다. "
    "동시에 anti-JCV index, 이전 면역억제제 사용, 투여 기간으로 층화한 PML hazard에 extended-interval dosing의 hazard ratio를 적용합니다.")

s = fig_slide("The dosing interval is the controllable handle", FIG / "s2_saturation.png", "Study ① · MS",
              cap=f"q6w extended-interval dosing (EID)에서도 relapse reduction {S2['red_q6w_pct']:.0f}% 유지 (q4w {S2['red_q4w_pct']:.0f}%)")
say(s, "왼쪽은 투여 후 다음 투여까지 α4-integrin saturation의 변화입니다. "
    "표준 q4w에서는 trough에서도 saturation이 75% 이상 유지되고, q6w에서는 다음 투여 직전에 50% 근처까지 내려옵니다. "
    "q8w처럼 더 늘리면 desaturation이 훨씬 깊어집니다. "
    f"오른쪽은 dosing interval에 따른 relapse-rate reduction이고, AFFIRM의 q4w 결과인 약 {S2['red_q4w_pct']:.0f}%에 맞춰 calibration했습니다. "
    f"[클릭] q6w extended-interval dosing에서도 relapse reduction은 {S2['red_q6w_pct']:.0f}%로 유지됩니다. "
    "즉 interval을 2주 늘려도 efficacy 손실은 크지 않습니다.")

s = fig_slide("PML risk is stratified, and steep", FIG / "s2_pml_table.png", "Study ① · MS",
              cap="anti-JCV index, prior immunosuppressant, 치료 기간에 따라 PML risk가 급격히 증가",
              box=(1.4, 1.75, 10.55, 4.3), src="Ho PR et al., Lancet Neurol 2017")
say(s, "PML risk는 환자에 따라 매우 다릅니다. "
    "이 표는 Ho 등이 2017년에 발표한 risk stratification으로, 1000명당 추정 발생률입니다. "
    "Anti-JCV 항체가 음성이면 거의 0에 가깝지만, 양성이면서 index가 1.5를 넘으면 투여 49개월에서 72개월 구간에 1000명당 10명까지 올라갑니다. "
    "이전에 면역억제제를 사용한 환자도 고위험군입니다. "
    "[클릭] 즉 PML risk는 anti-JCV index, 면역억제제 병력, 치료 기간에 따라 가파르게 증가합니다. "
    "그래서 고위험군일수록 dosing 전략을 바꾸는 효과가 커집니다.")

s = fig_slide(f"Extended-interval dosing cuts PML risk by {S2['pml_drop_pct']}%", FIG / "s2_benefit_risk.png",
              "Study ① · the optimization",
              cap=f"6년간 relapse {S2['relapse_cost_6y']:.1f}회 증가를 대가로 high-risk stratum PML {S2['pml_6y_highrisk_sid']:.1f} → {S2['pml_6y_highrisk_eid']:.1f} /1000",
              box=(1.3, 1.75, 10.75, 4.3))
say(s, f"이 risk table에 extended-interval dosing의 hazard ratio {S2['eid_pml_hr']}를 적용해 6년 누적 PML risk를 계산했습니다. "
    "빨간 점은 표준 q4w, 청록색 점은 q6w입니다. "
    f"Index 1.5 초과군에서는 1000명당 {S2['pml_6y_highrisk_sid']:.1f}명에서 {S2['pml_6y_highrisk_eid']:.1f}명으로, "
    f"이전 면역억제제 사용군에서는 {S2['pml_6y_priorIS_sid']:.1f}명에서 {S2['pml_6y_priorIS_eid']:.1f}명으로 줄어듭니다. "
    f"모든 stratum에서 약 {S2['pml_drop_pct']}%가 감소합니다. "
    f"[클릭] 대가는 6년 동안 relapse가 환자당 {S2['relapse_cost_6y']:.1f}회 늘어나는 정도입니다. "
    "Efficacy를 거의 잃지 않으면서 가장 심각한 safety risk를 4분의 1 수준으로 낮추는 것, 이것이 trade-off를 optimization으로 바꾼다는 뜻입니다.")

s = fig_slide("Where natalizumab q6w sits among MS therapies", FIG / "s2_dmt_map.png", "Study ① · MS",
              cap="EID는 natalizumab을 ocrelizumab에 가까운 efficacy로, 훨씬 낮은 risk 쪽으로 이동 (safety index는 illustrative)",
              box=(1.3, 1.75, 10.75, 4.3),
              src="Polman NEJM 2006 · Kappos NEJM 2010 · Hauser NEJM 2017 · Ryerson Neurology 2019")
say(s, "다른 MS 치료제와 비교하면 이 optimization의 의미가 더 잘 보입니다. "
    "x축은 serious safety event의 상대적 지표이고, y축은 placebo 대비 relapse-rate reduction입니다. "
    "이 safety index는 비교를 위한 illustrative 값이라는 점을 미리 말씀드립니다. "
    "Natalizumab q4w는 효과가 가장 높지만 risk도 가장 오른쪽에 있습니다. "
    "[클릭] Extended-interval dosing은 natalizumab을 ocrelizumab에 가까운 efficacy를 유지한 채 왼쪽, 즉 훨씬 낮은 risk 쪽으로 옮깁니다. "
    "같은 약이라도 dosing을 바꾸면 benefit–risk 지도 위의 위치가 달라집니다.")

# ---- Study 2: ICI ------------------------------------------------------------
s = intro_slide("Study ②  Checkpoint inhibitors: two antibodies", "Oncology · metastatic melanoma", [
    ("Nivolumab + ipilimumab", "melanoma에서 가장 효과적이면서 가장 독성이 큰 regimen"),
    ("Key biology", "ipilimumab의 anti-tumor priming 효과는 low dose에서 포화되지만, toxicity는 dose와 함께 계속 증가"),
    ("The question", "왜 ipilimumab을 낮추면 (CheckMate 511) efficacy는 유지되고 toxicity는 줄어드는가?"),
], "PD-1·CTLA-4 occupancy → anti-tumor / self-reactive effector pool (Bliss synergy) → ORR과 grade 3–4 irAE — 같은 CTLA-4 occupancy를 서로 다른 두 curve로 읽음",
    "Larkin NEJM 2015, 2019 (CheckMate 067) · Lebbé JCO 2019 (CheckMate 511) · Weber JCO 2017")
say(s, "두 번째는 oncology, 전이성 melanoma의 checkpoint inhibitor입니다. "
    "Nivolumab과 ipilimumab 병용은 melanoma에서 가장 효과적인 regimen이지만 독성도 가장 큽니다. "
    "핵심 biology는 ipilimumab의 anti-tumor priming 효과가 낮은 dose에서 이미 포화되는 반면, toxicity는 dose와 함께 계속 증가한다는 점입니다. "
    "그렇다면 ipilimumab을 낮췄을 때 효과는 유지되고 독성만 줄어야 하며, 실제로 CheckMate 511이 이를 보여주었습니다. "
    "[클릭] 모델은 PD-1과 CTLA-4 occupancy를 anti-tumor effector pool과 self-reactive effector pool로 연결하고, 이를 각각 ORR과 grade 3–4 irAE로 바꿉니다. "
    "핵심은 같은 CTLA-4 occupancy를 efficacy와 toxicity가 서로 다른 curve로 읽는다는 것입니다.")

r13 = S3["regimens"]["NIVO1 / IPI3"]; r31 = S3["regimens"]["NIVO3 / IPI1"]
s = fig_slide("The efficacy–safety frontier of combination dosing", FIG / "s3_pareto.png",
              "Study ② · the optimization",
              cap=f"NIVO1/IPI3 → NIVO3/IPI1: grade 3–4 irAE 약 {S3['cm511_tox_drop_pts']:.0f}%p 감소, ORR은 유지 ({r13['orr']:.0f}% → {r31['orr']:.0f}%)",
              box=(1.1, 1.75, 11.15, 4.3))
say(s, "이 그림은 두 약의 dose 조합을 바꿔 가며 계산한 efficacy–safety frontier입니다. "
    "x축은 grade 3–4 irAE, y축은 objective response rate이고, 작은 점의 색은 ipilimumab dose입니다. "
    "진한 선이 frontier, 즉 같은 toxicity에서 얻을 수 있는 최대 효과입니다. "
    "Monotherapy 결과는 CheckMate 067에, 병용 두 arm은 CheckMate 511에 calibration했습니다. "
    f"[클릭] NIVO1/IPI3에서 NIVO3/IPI1로 바꾸면 grade 3–4 irAE가 약 {S3['cm511_tox_drop_pts']:.0f}%p 줄고, ORR은 {r13['orr']:.0f}%에서 {r31['orr']:.0f}%로 유지됩니다. "
    "NIVO3/IPI1은 frontier 위에 있고 NIVO1/IPI3는 frontier 아래에 있다는 점도 눈여겨보십시오.")

s = new_slide("Response and toxicity run on two clocks", "Study ② · oncology")
add_image(s, FIG / "s3_timecourse.png", 0.75, 1.8, 5.75, 4.25)
reveal(s)
add_image(s, FIG / "s3_organ_timing.png", 6.8, 1.8, 5.8, 4.25)
caption(s, "Tumor response가 먼저, irAE는 뒤이어 organ별로 다른 시점에 발생 — safety monitoring을 위한 kinetic signature")
source(s, "Weber JS et al., J Clin Oncol 2017")
say(s, "같은 regimen 안에서도 response와 toxicity는 서로 다른 시계로 움직입니다. "
    "왼쪽 위는 tumor burden으로, 치료 20주 전후에 partial response 기준인 30% 감소를 넘습니다. "
    "아래는 effector pool로, anti-tumor effector가 먼저 올라가고 irAE를 일으키는 self-reactive effector는 더 천천히 쌓입니다. "
    "[클릭] 오른쪽은 Weber 등이 보고한 organ별 irAE 발생 시점입니다. "
    "피부는 3주에서 4주, 위장관과 간은 6주에서 8주, 내분비와 신장은 10주에서 12주 무렵에 주로 나타납니다. "
    "이 kinetic signature를 알면 언제 어떤 장기를 monitoring해야 하는지 미리 계획할 수 있습니다.")

# ---- Study 3: obesity --------------------------------------------------------
s = intro_slide("Study ③  Incretin therapy: a three-way optimization", "Obesity / metabolic · aging subgroup", [
    ("Semaglutide & tirzepatide", "조절 가능한 handle은 target dose와 escalation schedule"),
    ("Three competing outcomes", "weight loss (benefit) vs GI adverse event (tolerability) vs lean-mass loss (safety)"),
    ("The aging angle", "고령자에서 lean-mass loss는 sarcopenia·frailty risk로 직결"),
], "exposure → appetite suppression → energy-balance ODE; 체중 변화를 fat / fat-free mass로 분배 (age-sensitive Forbes rule); GI AE = escalation이 유발하고 tolerance가 완화하는 nausea state",
    "Wilding NEJM 2021 (STEP 1) · Jastreboff NEJM 2022 (SURMOUNT-1)")
say(s, "세 번째는 obesity의 incretin therapy입니다. "
    "Semaglutide와 tirzepatide에서 조절할 수 있는 handle은 target dose와 escalation schedule입니다. "
    "그런데 결과는 세 가지가 서로 경쟁합니다. "
    "Weight loss라는 benefit, GI adverse event라는 tolerability, 그리고 lean-mass loss라는 safety입니다. "
    "특히 고령자에서 lean-mass loss는 sarcopenia와 frailty risk로 직결됩니다. "
    "[클릭] 모델은 exposure가 appetite suppression을 거쳐 energy-balance ODE로 이어지게 했고, 체중 변화를 fat mass와 fat-free mass로 나누는 Forbes rule에 나이 효과를 넣었습니다. "
    "GI AE는 escalation이 유발하고 tolerance가 완화하는 nausea state로 표현했습니다.")

s = fig_slide("Titration is the handle: fast escalation buys nausea", FIG / "s4_titration.png", "Study ③ · obesity",
              cap="2-weekly escalation은 weight loss 이득 없이 초기 nausea peak만 키움")
say(s, "먼저 titration 속도의 영향입니다. "
    "왼쪽은 semaglutide를 4주마다 증량할 때의 dose 단계와 appetite suppression으로, dose가 오를 때마다 exposure가 계단식으로 올라갑니다. "
    "오른쪽은 nausea intensity로, 4주 escalation과 2주 escalation을 비교했습니다. "
    "2주마다 빠르게 올리면 초기 nausea peak가 훨씬 높아지고, 이후 tolerance가 생기면서 두 schedule 모두 같은 수준으로 수렴합니다. "
    "[클릭] 결국 빠른 escalation은 최종 weight loss에는 이득이 거의 없고 초기 nausea만 키웁니다. "
    "Titration schedule 자체가 tolerability를 조절하는 handle입니다.")

s = fig_slide("More weight loss costs more GI burden", FIG / "s4_frontier.png", "Study ③ · the optimization",
              cap=f"72주 weight loss: semaglutide 2.4 mg {S4['sema24_pct']}% (STEP 1: 14.9%), tirzepatide 15 mg {S4['tirze15_pct']}% (SURMOUNT-1: 20.9%)",
              box=(1.2, 1.75, 10.95, 4.3))
say(s, "다음은 target dose를 바꿨을 때의 efficacy–tolerability frontier입니다. "
    "x축은 escalation 기간의 누적 GI AE burden, y축은 72주 체중 감소율입니다. "
    f"모델은 semaglutide 2.4 mg에서 {S4['sema24_pct']}%를 예측했고, 이는 STEP 1의 14.9%와 비슷합니다. "
    f"Tirzepatide 15 mg은 {S4['tirze15_pct']}%로 SURMOUNT-1의 20.9%와 가깝습니다. "
    "[클릭] 두 곡선 모두 체중을 더 줄일수록 GI burden이 커지는 trade-off를 보입니다. "
    "같은 GI burden에서 tirzepatide 곡선이 더 위에 있어, tolerability 한 단위당 얻는 weight loss가 더 큽니다.")

s = new_slide("The aging safety axis: lean-mass loss", "Study ③ · obesity + aging")
add_image(s, FIG / "s4_composition.png", 0.75, 1.8, 5.8, 4.25)
reveal(s)
add_image(s, FIG / "s4_aging.png", 6.8, 1.8, 5.8, 4.25)
caption(s, f"Weight loss 중 lean mass 비율이 45세 {S4['lean_frac_young']}% → 70세 {S4['lean_frac_old']}%로 증가 — 고령자 dose 조절 signal")
say(s, "세 번째 축인 lean-mass loss는 나이에 따라 달라집니다. "
    "왼쪽은 tirzepatide 15 mg을 72주 투여한 뒤 잃은 fat mass와 lean mass입니다. "
    f"45세에서는 lean mass가 {S4['tirze15_lean_young']} kg, 70세에서는 {S4['tirze15_lean_old']} kg으로 늘어납니다. "
    "Muscle-preserving intervention을 함께 쓴다고 가정하면 lean mass 손실을 크게 줄일 수 있습니다. "
    f"[클릭] 오른쪽은 나이에 따른 lean share로, 체중 감소분 가운데 lean mass 비율이 45세 {S4['lean_frac_young']}%에서 70세 {S4['lean_frac_old']}%로 올라갑니다. "
    "두 약 모두 같은 경향을 보이므로, 고령자에서는 dose를 낮추거나 근육 보존 전략을 함께 쓰라는 signal로 읽을 수 있습니다.")

# ---- Study 4: DOAC -----------------------------------------------------------
s = intro_slide("Study ④  DOACs in the elderly: dosing the excluded", "Aging / renal and hepatic impairment", [
    ("The gap", "Phase 1은 젊고 renal function이 정상인 volunteer에서 DOAC exposure를 확인"),
    ("The reality", "실제 DOAC 복용 환자는 고령·renal impairment — 바로 임상시험에서 제외된 집단"),
    ("The question", "PBPK로 연구되지 않은 환자의 exposure를 예측하고 dose를 individualize할 수 있는가?"),
], f"age → creatinine clearance → DOAC exposure (renal-fraction scaling) → bleeding·stroke risk; {S5['n_cohort']:,}명 virtual elderly cohort에서 fixed vs individualized dosing 비교",
    "US product labels (dabigatran, apixaban, rivaroxaban, edoxaban) · Rowland & Tozer")
say(s, "네 번째는 고령자의 DOAC입니다. "
    "Phase 1은 젊고 신기능이 정상인 volunteer에서 DOAC exposure를 확인합니다. "
    "하지만 실제로 DOAC를 복용하는 환자는 대부분 고령이고 신기능이 떨어져 있어, 바로 임상시험에서 빠진 집단입니다. "
    "질문은 PBPK로 연구되지 않은 환자의 exposure를 예측하고 dose를 개별화할 수 있느냐입니다. "
    "[클릭] 모델은 나이에서 creatinine clearance, exposure, bleeding과 stroke risk로 이어지는 사슬을 만들었습니다. "
    f"그리고 {S5['n_cohort']:,}명의 virtual elderly cohort에서 fixed dose와 개별화된 dose를 비교했습니다.")

s = fig_slide("Age is a renal-function covariate", FIG / "s5_age_crcl.png", "Study ④ · aging",
              cap=f"Population mean CrCl은 {S5['age_cross_crcl50']}세에 dose-reduction threshold (50 mL/min) 아래로 — 측정 전에 exposure는 이미 상승",
              box=(1.2, 1.75, 10.95, 4.3))
say(s, "먼저 나이 자체가 renal function의 covariate라는 점입니다. "
    "회색 점은 virtual population이고, 진한 선은 나이에 따른 평균 creatinine clearance입니다. "
    f"평균 CrCl은 나이가 들수록 꾸준히 떨어져 {S5['age_cross_crcl50']}세 무렵 dose-reduction 기준인 50 mL/min 아래로 내려갑니다. "
    "개인 간 편차도 커서 그보다 젊은 나이에 이미 기준 아래인 환자도 많습니다. "
    "[클릭] 즉 신기능을 측정하기 전에 exposure는 이미 올라가 있을 수 있습니다. "
    "그래서 고령자 dosing에는 나이만이 아니라 신기능에 기반한 예측이 필요합니다.")

s = new_slide("Validating exposure, and finding the missing mechanism", "Study ④ · aging", size=25)
add_image(s, FIG / "s5_validation.png", 0.6, 1.8, 7.0, 4.55)
reveal(s)
bullets(s, [
    ("Renal scaling alone", f"label의 renal-impairment PK를 GMFE {S5['gmfe_renal_only']:.2f}×로 예측"),
    ("Dabigatran under-predicted", f"moderate RI에서 observed {S5['dabi_moderate_obs']}× vs predicted {S5['dabi_moderate_pred_renal']:.1f}× — GFR만으로 설명 불가"),
    ("Mechanism, not curve-fitting", f"P-gp / absorption term 추가 → GMFE {S5['gmfe_corrected']:.2f}×"),
], 7.85, 2.4, 4.75, 3.8, size=16, gap=14)
source(s, "US product labels: dabigatran, apixaban, rivaroxaban, edoxaban (renal-impairment PK)")
say(s, "Exposure 예측이 맞는지 drug label의 renal-impairment PK로 검증했습니다. "
    "네 가지 DOAC의 경증, 중등도, 중증 신장애에서 AUC 증가 배수를 renal-fraction scaling으로 예측하고 관측값과 비교했습니다. "
    f"[클릭] Renal scaling만으로 GMFE는 {S5['gmfe_renal_only']:.2f}배였습니다. "
    f"그런데 dabigatran은 중등도 신장애에서 관측값이 {S5['dabi_moderate_obs']}배인데 예측은 {S5['dabi_moderate_pred_renal']:.1f}배로, GFR만으로는 설명되지 않았습니다. "
    f"여기에 P-gp와 흡수 관련 mechanistic term을 더하자 GMFE가 {S5['gmfe_corrected']:.2f}배로 좋아졌습니다. "
    "중요한 것은 curve-fitting이 아니라 mechanism에 근거한 보정이라는 점입니다.")

s = new_slide("Two risks oppose; individualized dosing wins", "Study ④ · the optimization")
add_image(s, FIG / "s5_exposure_response.png", 0.75, 1.8, 5.75, 4.25)
reveal(s)
add_image(s, FIG / "s5_individualization.png", 6.8, 1.8, 5.8, 4.25)
caption(s, f"Model-informed dose reduction: over-exposure {S5['overexposed_fixed_pct']}% → {S5['overexposed_adj_pct']}%, modeled major bleeding {S5['bleed_fixed']} → {S5['bleed_adj']} /100 PY")
say(s, "Exposure가 정해지면 위험은 두 방향으로 움직입니다. "
    "왼쪽에서 exposure가 높아지면 major bleeding은 증가하고 stroke와 systemic embolism은 감소하며, 80세에서는 두 곡선 모두 위로 올라갑니다. "
    "녹색 구간이 두 위험의 균형을 맞추는 target exposure입니다. "
    "[클릭] 오른쪽은 virtual elderly cohort에서 dabigatran 150 mg bid를 고정 투여한 경우와 label의 dose-reduction rule을 적용한 경우입니다. "
    f"Exposure ceiling을 넘는 환자 비율이 {S5['overexposed_fixed_pct']}%에서 {S5['overexposed_adj_pct']}%로 줄었고, "
    f"모델로 계산한 major bleeding은 100 patient-years당 {S5['bleed_fixed']}건에서 {S5['bleed_adj']}건으로 감소했습니다. "
    "임상시험에 포함되지 않았던 환자 집단의 dosing을 모델로 평가한 사례입니다.")

# =============================================================================
# PART 5 -- flagship + NAMs bridge
# =============================================================================
s = divider(5, "From optimization to the clinic, and to NAMs",
            "Flagship merigolix case, 그리고 gut–liver–kidney MPS + PBPK bridge")
say(s, "다음은 optimization을 실제 임상 개발과 NAM으로 연결하는 사례입니다. "
    "[클릭] 먼저 저희가 직접 참여한 flagship case인 merigolix를 보겠습니다. "
    "이어서 gut–liver–kidney MPS와 PBPK를 결합한 NAM bridge를 소개하겠습니다. "
    "두 사례 모두 앞에서 본 원칙, 즉 mechanism으로 signal을 human decision으로 translate한다는 원칙을 따릅니다. "
    "하나는 dose 선택, 다른 하나는 in vitro data의 human 외삽이라는 점만 다릅니다.")

s = new_slide("Flagship: merigolix and the estrogen-threshold window", "Gynecology · GnRH antagonist", size=25)
bullets(s, [
    ("Mechanism", "oral GnRH receptor antagonist → LH/FSH ↓ → estradiol ↓ → lesion·pain 감소"),
    ("The trade-off", "E2 suppression이 부족하면 효과가 없고, 과하면 BMD loss와 hot flash"),
    ("The target", "partial E2 window (20–40 pg/mL) — estrogen threshold hypothesis"),
], 0.75, 1.9, 5.9, 3.0, size=16.5, gap=12)
reveal(s)
rrect(s, 0.75, 4.75, 5.9, 1.55, PALETEAL, radius=0.08)
textbox(s, "Phase 1 PK/PD model (SAD/MAD, n = 51) 위에 HPG axis·menstrual cycle·lesion·BMD를 결합한 24-state QSP model",
        1.0, 4.75, 5.4, 1.55, size=15.5, color=NAVY, family=F_MED, anchor=MSO_ANCHOR.MIDDLE, spacing=1.2)
ix, iy, iw, ih = add_image(s, FIG / "ext/merigolix_app.png", 6.95, 1.95, 5.65, 4.0)
textbox(s, "Interactive QSP dashboard · pipetqsp.shinyapps.io/merigolix", ix, iy + ih + 0.1, iw, 0.3,
        size=11.5, color=GREY, align=PP_ALIGN.CENTER, italic=True, spacing=1.0)
source(s, "Bae, Kang et al., Clin Pharmacokinet (in revision) · github.com/pipetcpt/merigolix-qsp")
say(s, "Merigolix는 endometriosis 치료를 위한 경구 GnRH receptor antagonist입니다. "
    "LH와 FSH를 낮춰 estradiol을 억제하고, 그 결과 lesion과 pain이 줄어듭니다. "
    "그런데 E2를 너무 적게 억제하면 효과가 없고, 너무 많이 억제하면 BMD loss와 hot flash가 생깁니다. "
    "그래서 E2를 20에서 40 pg/mL 정도로 부분 억제하자는 것이 estrogen threshold hypothesis입니다. "
    "[클릭] 저희는 Phase 1 SAD/MAD 51명의 PK/PD model 위에 HPG axis, menstrual cycle, lesion, BMD를 결합한 24-state QSP model을 만들었습니다. "
    "오른쪽은 누구나 dose를 바꿔 볼 수 있는 interactive dashboard입니다.")

s = fig_slide("Merigolix: the trade-off becomes dose optimization", FIG / "s6_merigolix.png",
              "Flagship · the named example",
              cap=f"{S6['best_dose']} mg QD에서 composite therapeutic index 최고 ({S6['best_ti']}) — 120 mg은 tolerability 우선 대안",
              box=(0.75, 1.75, 11.85, 4.3), src="Merigolix QSP model, virtual population at week 12 (github.com/pipetcpt/merigolix-qsp)")
say(s, "Virtual population에서 12주 결과를 dose별로 보겠습니다. "
    "왼쪽은 12주 estradiol 분포로, dose가 오를수록 중앙값이 target window 안으로 들어오고 240 mg에서는 window 아래로 내려갑니다. "
    "가운데는 window 안에 있는 환자 비율과 20 pg/mL 미만으로 과도하게 억제된 환자 비율입니다. "
    f"Window 비율은 120 mg에서 {S6['in_window_120']}%로 가장 높고, 과억제 비율은 dose와 함께 계속 늘어 240 mg에서 {S6['below20_240']}%에 이릅니다. "
    f"[클릭] 오른쪽 composite therapeutic index는 {S6['best_dose']} mg QD에서 {S6['best_ti']}로 가장 높았습니다. "
    "Tolerability를 우선한다면 120 mg이 대안이 됩니다.")

s = new_slide("Connecting to NAMs: gut–liver–kidney MPS + PBPK", "NAMs bridge", size=25)
for i, (num, lab) in enumerate([("6", "drugs"), ("24", "clinical observations"),
                                 ("0", "clinical PK used to fit"), ("24/24", "within 4-fold")]):
    cx = 0.75 + i * 3.0
    textbox(s, [[(num + "  ", {"family": F_BOLD, "size": 24, "color": TEAL}),
                 (lab, {"size": 14.5, "color": INK})]], cx, 1.78, 2.9, 0.5, spacing=1.0)
add_image(s, FIG / "glk_validation.png", 0.75, 2.3, 11.85, 3.75)
reveal(s)
caption(s, "Chip-fitted parameter만으로 human AUC 예측 (GMFE 1.85×) — industry PBPK benchmark (2–5×)보다 작은 오차",
        y=6.15, size=14.5)
source(s, "hiEC–HepaRG–RPTEC gut–liver–kidney MPS + whole-body PBPK (PIPET MPS)")
say(s, "이제 NAM bridge입니다. "
    "저희 PIPET MPS 팀은 hiEC, HepaRG, RPTEC으로 구성된 gut–liver–kidney MPS를 whole-body PBPK와 결합했습니다. "
    "6개 약물, 24개 clinical observation을 예측했고, 모델을 맞추는 데 clinical PK는 전혀 쓰지 않았습니다. "
    "왼쪽은 약물별 AUC GMFE로, 6개 약물 중 4개가 2배 이내이고 24개 관측 모두 4배 이내입니다. "
    "오른쪽은 scaling rule의 영향으로, 모든 약물에 같은 보정을 적용하면 GMFE가 9.7배로 나빠지지만 mechanism에 따른 selective rule을 쓰면 1.85배가 됩니다. "
    "[클릭] 즉 chip에서 얻은 parameter만으로 human AUC를 GMFE 1.85배로 예측했습니다. "
    "Industry-standard PBPK platform에서 보고되는 2배에서 5배 수준과 비교해도 오차가 작습니다.")

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
reveal(s)
bullets(s, [
    ("Same pipeline", "disease chip으로 re-fit하되 영향받지 않는 organ의 parameter는 고정"),
    ("The payoff", "Phase 1에서 제외되는 IBD·AKI 환자를 위한 testable hypothesis — Study ④와 같은 문제의식"),
], 0.75, 4.5, 11.85, 1.3, size=17, gap=10)
callout(s, "새 case study와 MPS/NAM 연구는 하나의 program: mechanism이 in vitro·virtual signal을 human safety decision으로 translate",
        y=5.6, h=0.8, size=16)
say(s, "같은 pipeline을 disease-state chip으로 확장했습니다. "
    "Gefitinib을 DSS colitis chip에서 평가하면 AUC가 0.49배로 줄었고, Fa와 Fh는 그대로인 채 Fg만 변한 gut-localized change였습니다. "
    "Simvastatin을 cisplatin-injury chip에서 평가하면 AUC가 2.19배로 늘었고, 간의 Kp,uu가 절반으로 줄고 OATP1B1이 억제된 hepatic-localized change였습니다. "
    "[클릭] 이때 영향받지 않는 장기의 parameter는 고정하고 같은 pipeline으로 다시 fitting했습니다. "
    "이 결과는 Phase 1에서 제외되는 IBD나 AKI 환자를 위한 testable hypothesis가 되고, 앞의 DOAC study와 같은 문제의식입니다. "
    "결국 새 case study와 MPS/NAM 연구는 하나의 program입니다.")

s = new_slide("Prerequisites for trusting a chip-derived number", "Discipline for regulatory acceptance", size=25)
cards(s, [
    ("Chip-specific behavior", "plastic·membrane adsorption 보정", TEAL),
    ("Dimensional scaling", "organ별 cell·fluidic scaling factor를 정확히", NAVY),
], h=1.5, numbered=True, body_size=15.5)
reveal(s)
cards(s, [
    ("Held-out clinical validation", "fitting에 쓰지 않은 clinical data로 검증", GOLD),
    ("Falsifiable correction rules", "curve-fitting이 아닌 mechanism-based rule", RED),
], h=1.5, numbered=True, body_size=15.5, y0=1.85 + 1.5 + 0.3, start=2)
callout(s, "Study ④의 DOAC validation도 같은 원칙 — fitted fudge factor가 아닌 mechanistic (P-gp) correction",
        y=5.6, h=0.8, size=16)
say(s, "Chip에서 나온 숫자를 PBPK나 QSP model에 넣으려면 네 가지를 먼저 확인해야 합니다. "
    "첫째, plastic이나 membrane에 대한 흡착처럼 chip 특유의 약물 거동을 보정해야 합니다. "
    "둘째, 장기별 cell 수와 fluidic dimension의 scaling factor가 정확해야 합니다. "
    "[클릭] 셋째, fitting에 쓰지 않은 clinical data로 검증해야 합니다. "
    "넷째, 보정 규칙은 curve-fitting이 아니라 반증 가능한 mechanism 기반 규칙이어야 합니다. "
    "앞의 DOAC validation도 같은 원칙으로, fudge factor가 아니라 P-gp라는 mechanism으로 보정했습니다.")

# =============================================================================
# PART 6 -- regulatory acceptance, people, conclusion
# =============================================================================
s = divider(6, "Prerequisites, people, and the path forward",
            "Context of use, clinical data 기반 validation — 그리고 이를 translate할 전문 인력")
say(s, "마지막으로 규제 수용의 조건과 사람의 문제, 그리고 앞으로의 방향을 말씀드리겠습니다. "
    "[클릭] 첫째는 context of use를 정하고 clinical data로 검증하는 것입니다. "
    "둘째는 이를 실제로 translate할 전문 인력을 키우는 것입니다. "
    "저는 기술 자체보다 이 두 가지가 실제 적용 속도를 결정한다고 생각합니다. "
    "두 가지를 차례로 보겠습니다.")

s = new_slide("What regulatory acceptance will require", "Context of use")
table(s, [
    ["Principle", "What it asks for"],
    ["Context of use", "The precise question each model or chip is qualified to answer"],
    ["Biological relevance", "The human mechanism the NAM reproduces"],
    ["Technical characterization", "Reproducibility, dynamic range, failure modes"],
    ["Fit-for-purpose validation", "Agreement with held-out clinical data (the step QSP / PBPK supplies)"],
], 0.75, 1.85, 11.85, col_w=[3.3, 8.5], row_h=0.62, size=15)
reveal(s)
callout(s, "QSP/PBPK는 이미 확립된 MIDD pathway — NAM/MPS data를 여기에 통합하는 것이 decision-grade safety evidence로 가는 가장 현실적인 길",
        y=5.25, h=1.05)
say(s, "NAM 기반 evidence가 규제에서 받아들여지려면 FDA가 제시한 네 가지 원칙을 충족해야 합니다. "
    "Context of use는 각 model이나 chip이 답할 수 있는 질문을 정확히 정의하는 것입니다. "
    "Biological relevance는 NAM이 재현하는 human mechanism이 무엇인지 밝히는 것이고, technical characterization은 재현성과 dynamic range, failure mode를 확인하는 것입니다. "
    "마지막 fit-for-purpose validation은 held-out clinical data와의 일치를 보이는 단계이며, 바로 QSP와 PBPK가 맡는 부분입니다. "
    "[클릭] QSP와 PBPK는 이미 확립된 MIDD pathway입니다. "
    "NAM과 MPS data를 이 경로에 통합하는 것이 decision-grade safety evidence로 가는 가장 현실적인 길이라고 생각합니다.")

s = new_slide("The real bottleneck is people, not technology", "Unmet needs in our own field")
bullets(s, [
    ("Workforce shortage", "NAM signal을 PBPK/QSP로 해석할 clinical pharmacologist 부족"),
    ("No standard workflow", "FDA의 네 가지 NAM validation principle을 위한 표준 workflow 부재"),
    ("Few head-to-head validations", "chip·organoid model을 human clinical data와 직접 비교한 연구가 드묾"),
], 0.75, 2.1, 11.85, 2.3, size=19, gap=14)
reveal(s)
rrect(s, 0.75, 4.75, 11.85, 1.1, NAVY, radius=0.1)
textbox(s, "Not a data problem, not a technology problem: a workforce problem.", 0.75, 4.75, 11.85, 1.1,
        size=21, color=WHITE, family=F_MED, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE, spacing=1.0)
say(s, "그렇다면 무엇이 병목일까요? "
    "첫째, NAM signal을 PBPK와 QSP로 해석할 수 있는 clinical pharmacologist가 부족합니다. "
    "둘째, FDA의 네 가지 NAM validation principle을 실제로 수행하는 표준 workflow가 아직 없습니다. "
    "셋째, chip이나 organoid model을 human clinical data와 직접 비교한 연구가 드뭅니다. "
    "[클릭] 정리하면 이것은 data의 문제도, technology의 문제도 아닌 workforce의 문제입니다. "
    "저희 분야인 임상약리학이 이 인력을 길러내야 한다고 생각합니다.")

s = new_slide("Conclusion", "What to take away")
for i, (lead, body, col) in enumerate([
        ("QSP + PBPK are the bridge", "in vitro·virtual signal을 human exposure, safety margin, benefit–risk로 translate", NAVY),
        ("Trade-offs become optimizations", "MS dosing interval · ICI dose ratio · incretin titration · DOAC individualization · merigolix window", TEAL),
        ("At scale", f"LLM이 구축한 {N_MODELS}개 open QSP model library — {PCT_SAFE}%가 이미 explicit safety axis 보유", BLUE),
        ("Validated honestly", "gut–liver–kidney MPS가 clinical PK 없이 human AUC를 4-fold 이내로 예측 — fudge factor가 아닌 mechanism", GOLD),
        ("The gap is people", "NAM-derived signal을 clinical safety decision으로 연결할 전문 인력이 필요", RED)]):
    if i == 3:
        reveal(s)
    cy = 1.8 + i * 0.93
    circ = s.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.8), Inches(cy + 0.1), Inches(0.58), Inches(0.58))
    circ.fill.solid(); circ.fill.fore_color.rgb = col; circ.line.fill.background(); circ.shadow.inherit = False
    textbox(s, str(i + 1), 0.8, cy + 0.1, 0.58, 0.58, size=18, color=WHITE, family=F_BOLD,
            align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE, spacing=1.0)
    textbox(s, lead, 1.65, cy + 0.02, 10.9, 0.4, size=18, color=col, family=F_BOLD, spacing=1.0)
    textbox(s, body, 1.65, cy + 0.42, 10.9, 0.45, size=15.5, color=INK, spacing=1.05)
say(s, "오늘 말씀드린 내용을 다섯 가지로 정리하겠습니다. "
    "첫째, QSP와 PBPK는 in vitro와 virtual signal을 human exposure, safety margin, benefit–risk로 translate하는 다리입니다. "
    "둘째, MS dosing interval, ICI dose ratio, incretin titration, DOAC 개별화, merigolix window처럼 trade-off를 optimization 문제로 바꿀 수 있습니다. "
    f"셋째, LLM이 구축한 {N_MODELS}개의 공개 QSP model library가 있고, 그중 {PCT_SAFE}%가 이미 safety axis를 갖고 있습니다. "
    "[클릭] 넷째, gut–liver–kidney MPS는 clinical PK 없이 human AUC를 4배 이내로 예측했고, 보정은 fudge factor가 아닌 mechanism에 근거했습니다. "
    "다섯째, 남은 과제는 사람이며, NAM에서 나온 signal을 clinical safety decision으로 연결할 전문 인력이 필요합니다.")

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
reveal(s)
rrect(s, 8.65, 3.8, 3.95, 1.35, PALE, radius=0.1)
textbox(s, [[("Fully open", {"family": F_BOLD, "size": 14, "color": GOLD})],
            [("github.com/pipetcpt/qsp", {"family": F_MED, "size": 16, "color": NAVY})]],
        8.9, 3.8, 3.6, 1.35, anchor=MSO_ANCHOR.MIDDLE, spacing=1.3)
say(s, "들어주셔서 감사합니다. "
    "오늘 보여드린 모델과 코드, 그림은 모두 공개되어 있습니다. "
    "[클릭] github.com/pipetcpt/qsp에서 library 전체를 보실 수 있고, 각 질환 모델의 R code를 바로 내려받아 실행해 볼 수 있습니다. "
    "공동연구 제안이나 질문은 언제든 환영합니다. "
    "이제 질문 받겠습니다.")

# =============================================================================
# one-click reveal timing + speaker notes
# =============================================================================
n_anim = 0
for s, start in REVEAL.values():
    n_anim += apply_reveal(s, start) > 0
short = []
for s, script in NOTES.values():
    write_notes(s, script)
    if n_sentences(script) < 5:
        short.append(script[:40])
    if (id(s._element) in REVEAL) != ("[클릭]" in script):
        raise SystemExit(f"click marker / animation mismatch: {script[:50]}")
missing = [i + 1 for i, sl in enumerate(prs.slides) if id(sl._element) not in NOTES]
if short or missing:
    raise SystemExit(f"notes check failed: short={short} missing={missing}")
print(f"animated slides: {n_anim}; slides with notes: {len(NOTES)}")

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
