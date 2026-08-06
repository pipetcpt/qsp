#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Re-run ONE section of hdfn_analysis.py alone, verbatim.

WHY THIS FILE EXISTS.  The first version of the IVIG mechanism decomposition
isolated the "FcgammaR only" arm by zeroing the IVIG POOL (ivig_f = 1e-9).
That removes both mechanisms at once -- there is no IVIG-derived IgG left to
block FcgammaR with -- so the arm looked inert and the prose drew the opposite
conclusion from the truth.  The model now carries an explicit `ivig_compete`
switch, and this script re-executes the corrected section so that the block in
hdfn_reference_output.txt is exactly what a full re-run of the analysis would
print.  The section-8 code is EXTRACTED FROM hdfn_analysis.py at run time, not
copied, so the two cannot drift apart.

    python3 hdfn_section_rerun.py 4 --splice
    python3 hdfn_section_rerun.py 16 --with=1,3 --splice

Two sections have needed this.  Section 8 (the IVIG decomposition, whose
"FcgammaR only" arm was switched off by zeroing the IVIG pool, which removes
both mechanisms) was fixed before the final run and is correct in place.
Section 4 (the MCA-PSV threshold) carried a hard-coded viscosity figure and a
mis-escaped format string, and is re-run here.
"""
import io
import re
import sys

SRC = open("hdfn_analysis.py").read()
SEC = int([a for a in sys.argv[1:] if a.isdigit()][0]) if any(
    a.isdigit() for a in sys.argv[1:]) else 8
BEG = SRC.index("hdr(%d, " % SEC)
try:
    END = SRC.index("hdr(%d, " % (SEC + 1))
except ValueError:          # the last section has no successor
    END = len(SRC)
BLOCK = SRC[BEG:END].rstrip()
BLOCK = BLOCK[:BLOCK.rindex("\n#")] if "\n#" in BLOCK else BLOCK

PRE = SRC[:SRC.index("# " + "=" * 78 + "\nhdr(0")]


def block_of(sec):
    a = SRC.index("hdr(%d, " % sec)
    try:
        b = SRC.index("hdr(%d, " % (sec + 1))
    except ValueError:
        b = len(SRC)
    return SRC[a:b]


def render():
    """Execute the prerequisites (whose output is discarded) and then the target
    section.  Section 16 quotes numbers computed in sections 1 and 3, so those
    have to be run first -- which is cheap, and much cheaper than re-running the
    whole analysis for one block."""
    pre = [int(x) for a in sys.argv[1:] if a.startswith("--with=")
           for x in a[7:].split(",")]
    buf = io.StringIO()
    old = sys.stdout
    sys.stdout = buf
    try:
        g = {"__name__": "__main__"}
        exec(compile(PRE, "pre", "exec"), g)
        for sec in pre:
            exec(compile(block_of(sec), "section%d" % sec, "exec"), g)
        buf.truncate(0)
        buf.seek(0)
        exec(compile(BLOCK, "section%d" % SEC, "exec"), g)
    finally:
        sys.stdout = old
    return buf.getvalue()


def main():
    text = render()
    if "--splice" not in sys.argv:
        sys.stdout.write(text)
        return
    out = open("hdfn_reference_output.txt").read()
    pat = (r"\n=+\n%d\.  .*?(?=\n=+\n%d\.  )" % (SEC, SEC + 1)
           if ("\n%d.  " % (SEC + 1)) in out else r"\n=+\n%d\.  .*\Z" % SEC)
    m = re.search(pat, out, re.S)
    if m is None:
        print("section 8 not found in hdfn_reference_output.txt; nothing done")
        return
    open("hdfn_reference_output.txt", "w").write(
        out[:m.start()] + "\n" + text.rstrip("\n") + out[m.end():])
    print("spliced the corrected section %d into hdfn_reference_output.txt" % SEC)


if __name__ == "__main__":
    main()
