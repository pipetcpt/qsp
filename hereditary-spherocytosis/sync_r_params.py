#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rewrite the numeric values in hsph_mrgsolve_model.R's $PARAM block from the
Python reference (P0 overridden by hsph_calibration.json), so the two files
cannot silently drift apart.  Names and annotations in the R file are left
alone; only the numbers are replaced, and any parameter present in Python but
missing from the R block is reported."""
import json
import re
import sys

from hsph_python_reference import P0

FIT = json.load(open("hsph_calibration.json"))
VAL = dict(P0)
VAL.update(FIT)

# Python name -> R name where they differ
ALIAS = {}

src = open("hsph_mrgsolve_model.R").read()
head, rest = src.split("$PARAM", 1)
pblock, tail = rest.split("$CMT", 1)

seen = set()
out_lines = []
for line in pblock.split("\n"):
    m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)(\s*):(\s*)(-?[0-9][0-9.eE+-]*)(\s*):(.*)$",
                 line)
    if not m:
        out_lines.append(line)
        continue
    name, s1, s2, old, s3, ann = m.groups()
    key = ALIAS.get(name, name)
    if key not in VAL:
        out_lines.append(line)
        print("  R parameter not in Python: %s" % name)
        continue
    seen.add(key)
    v = VAL[key]
    new = ("%.6g" % v) if (abs(v) >= 1e-4 or v == 0) else ("%.6g" % v)
    if new != old:
        print("  %-12s %-12s -> %s" % (name, old, new))
    pad = max(1, len(old) + len(s2) - len(new))
    out_lines.append("%s%s:%s%s%s:%s" % (name, s1, " " * pad, new, s3, ann))

missing = sorted(set(VAL) - seen)
if missing:
    print("  Python parameters absent from the R $PARAM block: %s"
          % ", ".join(missing))

open("hsph_mrgsolve_model.R", "w").write(
    head + "$PARAM" + "\n".join(out_lines) + "$CMT" + tail)
print("synced %d parameters into hsph_mrgsolve_model.R" % len(seen))
