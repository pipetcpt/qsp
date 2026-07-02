#!/usr/bin/env python3
"""
scripts/fix_readme_table.py — normalize the root README.md model gallery table.

WHY THIS EXISTS
----------------
On 2026-07-02 the gallery table broke partway through rendering on GitHub.
Root cause: several rows had picked up a stray, UNMATCHED '**' inside their
one-line summary (left over from citation-heavy prose the daily routine had
written, e.g. "...급성 발작**." with no opening '**'). An odd number of '**'
inside one Markdown table cell corrupts GFM's inline parsing for the rest of
that row and visually breaks every row that follows it. Separately, the
per-disease "상세 섹션" that used to get appended below the table had grown to
thousands of lines and the trailing Disclaimer/References/License footer had
been pushed out / lost entirely.

WHAT THIS SCRIPT DOES
----------------------
1. Scans every disease directory on disk (anything with a *_qsp*.dot or
   *_qsp*.svg file) — this is the authoritative list of models.
2. Parses the EXISTING gallery table rows to recover each directory's
   category / Korean title / English subtitle / summary text (so no manual
   content is lost — only the FORMATTING is normalized).
3. Cleans every summary cell:
     - strips ALL '*'/'**' (titles keep their own separately-generated bold
       pair, so cells can never end up with an odd count again)
     - strips citation clutter such as "(Author 2020 Journal[PMID 123])"
     - truncates to one clean sentence (~190 chars, cut at a sentence break)
4. Caps the English sub-title to its intended "Name · ABBR" form.
5. Normalizes the many ad-hoc category labels down to a fixed canonical set.
6. Rebuilds ONE contiguous table (rows numbered 1..N, always exactly 6
   columns), recomputes the per-category counts, and refreshes every stale
   "N개 모델" mention in the intro.
7. Restores the Disclaimer / References / License footer if it's missing.

USAGE
-----
    # Rebuild + overwrite README.md (do this after adding a new model row):
    python3 scripts/fix_readme_table.py

    # Lint only — check the CURRENT README.md for the exact defects above,
    # print PASS/FAIL, change nothing, exit 1 on failure (useful in CI or
    # right before `git commit`):
    python3 scripts/fix_readme_table.py --check

RULES FOR ADDING A NEW MODEL'S ROW (see CLAUDE.md "README 업데이트")
--------------------------------------------------------------------
When you add a new disease, append its row to the gallery table yourself
first (so the script has something to normalize), following the existing
row format exactly, then run this script with no flags before committing.
The summary cell must be: ONE plain sentence, no bold ('**'), no citations
or PMIDs, roughly under 190 characters. The English sub-title must stay
"English Name · ABBR" — do not append mechanism/drug keyword lists there.
"""
import argparse
import glob
import os
import re
import sys
from collections import Counter

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RD = os.path.join(ROOT, "README.md")

CITATION_PAT = re.compile(
    r"[\(\[][^()\[\]]*(?:PMID|20\d\d\s*[A-Za-z]|NEJM|Lancet|JAMA|Circulation|Blood|Nature)[^()\[\]]*[\)\]]"
)

CATEGORY_RULES = [
    ("혈관염", "혈관염"),
    ("종양", "종양"), ("암", "종양"),
    ("혈액", "혈액"),
    ("정신", "정신·신경"),
    ("신경근육", "근골격·신경근육"), ("근골격", "근골격·신경근육"),
    ("신경", "신경"),
    ("심혈관", "심혈관"), ("자율신경", "심혈관"),
    ("호흡기", "호흡기"), ("이비인후과", "호흡기"),
    ("신장", "신장·비뇨"), ("비뇨", "신장·비뇨"),
    ("소화기", "소화기·간담도"), ("간담도", "소화기·간담도"),
    ("내분비", "내분비·대사"), ("대사", "내분비·대사"),
    ("피부", "피부"),
    ("안과", "안과"), ("망막", "안과"),
    ("감염", "감염"),
    ("자가면역", "자가면역·류마티스"), ("류마티스", "자가면역·류마티스"),
    ("부인", "부인·생식"), ("산과", "부인·생식"), ("생식", "부인·생식"),
    ("희귀", "희귀·유전"), ("유전", "희귀·유전"),
]


def norm_cat(c):
    for k, v in CATEGORY_RULES:
        if k in c:
            return v
    return "기타"


def clean_summary(raw, cap=190):
    s = raw.strip()
    for _ in range(3):
        s = CITATION_PAT.sub("", s)
    s = s.replace("**", "")
    s = re.sub(r"(?<!\w)\*(?!\*)", "", s)
    s = re.sub(r"\s{2,}", " ", s)
    s = re.sub(r"\(\s*\)", "", s)
    s = re.sub(r"\s+([,.;)])", r"\1", s)
    s = s.strip(" -–—")
    if len(s) > cap:
        cut = s[:cap]
        m = list(re.finditer(r"[.。]\s", cut))
        s = cut[: m[-1].end()].strip() if m else cut.rstrip() + "…"
    return s.strip()


def cap_subtitle(sub, max_parts=2, max_len=80):
    parts = [p.strip() for p in sub.split("·")]
    if len(parts) > max_parts or len(sub) > max_len:
        sub = " · ".join(parts[:max_parts])
    return sub


def find_disk_dirs():
    dirs = []
    for d in sorted(os.listdir(ROOT)):
        p = os.path.join(ROOT, d)
        if not os.path.isdir(p) or d in ("talks", "docs", "scripts") or d.startswith("."):
            continue
        if glob.glob(os.path.join(p, "*_qsp*.dot")) or glob.glob(os.path.join(p, "*_qsp*.svg")):
            dirs.append(d)
    return dirs


def pick(d, *pats):
    for pat in pats:
        g = sorted(glob.glob(os.path.join(ROOT, d, pat)))
        g = [x for x in g if "shiny" not in os.path.basename(x).lower() or "shiny_app" in pat]
        if g:
            return os.path.basename(g[0])
    return None


def parse_existing_rows(lines, gidx):
    meta = {}
    for l in lines[gidx:]:
        m = re.match(r"^\|\s*(\d+)\s*\|", l)
        if not m:
            continue
        parts = l.split("|")
        if len(parts) < 7:
            continue
        cat = parts[2].strip()
        model_cell = parts[3]
        last_cell = parts[5]
        md = re.search(r"\]\(([a-z0-9][a-z0-9-]+)/\)", model_cell)
        mk = re.search(r"\[\*\*(.+?)\*\*", model_cell)
        ms = re.search(r"<sub>(.+?)</sub>", model_cell)
        if not (md and mk):
            continue
        d = md.group(1)
        raw_summary = last_cell.split("<br>")[0]
        if d not in meta:
            meta[d] = (cat, mk.group(1).strip(), (ms.group(1).strip() if ms else ""), raw_summary)
    return meta


def build_table(lines, gidx):
    meta = parse_existing_rows(lines, gidx)
    dirs = find_disk_dirs()

    rows, row_cats, missing = [], [], []
    for i, d in enumerate(dirs, 1):
        png = pick(d, "*_qsp_model.png", "*_qsp.png")
        svg = pick(d, "*_qsp_model.svg", "*_qsp.svg")
        R = pick(d, "*_mrgsolve_model.R", "*_model.R")
        refs = pick(d, "*_references.md", "*references.md")
        if not (png and svg):
            continue
        abbr = re.split(r"_qsp", png)[0]
        if d in meta:
            cat, ko, sub, raw_summary = meta[d]
        else:
            missing.append(d)
            cat, ko, sub, raw_summary = ("기타", d.replace("-", " ").title(), "", "")
        cat = norm_cat(cat)
        sub = cap_subtitle(sub)
        summ = clean_summary(raw_summary)
        row_cats.append(cat)
        title = f"[**{ko}**<br><sub>{sub}</sub>]({d}/)"
        img = f'<a href="{d}/{svg}"><img src="{d}/{png}" width="190" alt="{abbr}"></a>'
        links = f"[🗺️ 지도]({d}/{svg}) · [⚙️ mrgsolve]({d}/{R}) · [📚 문헌]({d}/{refs}) · [📄 README]({d}/README.md)"
        rows.append(f"| {i} | {cat} | {title} | {img} | {summ}<br>{links} |")
    return rows, row_cats, missing


def rebuild(check_only=False):
    src = open(RD, encoding="utf-8").read()
    lines = src.splitlines()
    gidx = next(i for i, l in enumerate(lines) if l.startswith("## 📚 모델 갤러리"))
    head = "\n".join(lines[:gidx]).rstrip()

    rows, row_cats, missing = build_table(lines, gidx)
    cats = Counter(row_cats)
    catline = " · ".join(f"{c} {n}" for c, n in cats.most_common())
    n = len(rows)

    if missing:
        print(f"WARNING: {len(missing)} director{'y has' if len(missing)==1 else 'ies have'} no existing "
              f"table row — add a row for them first, then re-run: {missing}", file=sys.stderr)

    head = re.sub(r"\d+개 질환", f"{n}개 질환", head)
    head = re.sub(r"총 \d+개", f"총 {n}개", head)
    head = re.sub(r"models-\d+", f"models-{n}", head)
    head = re.sub(r"매일 \d+개 질환", "매일 1개 질환", head)

    gallery = [
        "## 📚 모델 갤러리 (Model Gallery)",
        "",
        f"전체 **{n}개** QSP 모델입니다. 모델명을 클릭하면 해당 디렉토리로, 그림을 클릭하면 확대 가능한 SVG 지도로 이동합니다. "
        f"각 행의 링크에서 기계론적 지도(🗺️), mrgsolve 모델(⚙️), 참고문헌(📚), 상세 README(📄)에 바로 접근할 수 있습니다.",
        "",
        f"**분류별 모델 수**: {catline}",
        "",
        "| # | 분류 | 모델 | 미리보기 | 요약 및 링크 |",
        "|---|------|------|----------|--------------|",
        *rows,
    ]

    footer = """
---

## ⚠️ 면책 조항 (Disclaimer)

본 라이브러리의 모든 모델은 **교육 및 연구 목적의 정성적·반정량적 QSP 모델**입니다. 공개 문헌과 임상시험 데이터를 바탕으로 구성되었으나 독립적으로 검증·인증되지 않았으며, **실제 임상 의사결정, 처방, 또는 규제 제출에 직접 사용해서는 안 됩니다.** 파라미터와 가정은 설명을 위한 근사치이며, 실제 환자 데이터에 대한 적합·검증이 별도로 필요합니다.

## 📖 참고 자료 (References & Tools)

- mrgsolve를 이용한 R 기반 QSP: <https://vantage-research.net/qsp-in-r/>
- gPKPDviz — mrgsolve 기반 PK/PD 시뮬레이션 Shiny 도구
  - 논문: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/>
  - 코드: <https://github.com/Genentech/gPKPDviz/>

## 📄 라이선스 (License)

본 저장소의 라이선스는 [LICENSE](LICENSE) 파일을 참조하세요.
"""

    new_src = head + "\n\n" + "\n".join(gallery) + "\n" + footer

    if check_only:
        return lint(src, RD)

    open(RD, "w", encoding="utf-8").write(new_src)
    print(f"OK: rewrote README.md — {n} models, {len(cats)} categories, footer restored.")
    return 0


def lint(src, path_label):
    """Validate an EXISTING README.md without rewriting it. Returns exit code."""
    lines = src.splitlines()
    problems = []

    row_lines = [(i, l) for i, l in enumerate(lines, 1) if re.match(r"^\| \d+ \|", l)]
    if not row_lines:
        problems.append("no gallery rows found at all")
    else:
        nums = [int(re.match(r"^\| (\d+) \|", l).group(1)) for _, l in row_lines]
        if nums != list(range(1, len(nums) + 1)):
            problems.append(f"row numbers not contiguous 1..{len(nums)} (found {len(nums)} rows)")
        for i, l in row_lines:
            if l.count("|") != 6:
                problems.append(f"line {i}: expected 6 columns, found {l.count('|')}")
            if l.count("**") % 2 != 0:
                problems.append(f"line {i}: unbalanced '**' (count={l.count('**')}) — will break table rendering")

    if not re.search(r"면책 조항|Disclaimer", src):
        problems.append("footer (Disclaimer/References/License) is missing")

    # broken relative links
    for m in re.finditer(r'src="([^"]+)"', src):
        p = os.path.join(ROOT, m.group(1))
        if not os.path.isfile(p):
            problems.append(f"broken image link: {m.group(1)}")
    for m in re.finditer(r"\]\(([a-z0-9][^)]+\.(?:svg|R|md))\)", src):
        p = os.path.join(ROOT, m.group(1))
        if not os.path.isfile(p):
            problems.append(f"broken file link: {m.group(1)}")

    if problems:
        print(f"FAIL ({path_label}): {len(problems)} problem(s):", file=sys.stderr)
        for p in problems[:50]:
            print(" -", p, file=sys.stderr)
        return 1
    print(f"PASS ({path_label}): gallery table looks structurally sound.")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="lint only, don't rewrite README.md")
    args = ap.parse_args()
    sys.exit(rebuild(check_only=args.check))
