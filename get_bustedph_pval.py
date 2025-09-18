#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Extract p-values from HyPhy BUSTED-PH outputs (JSON preferred, OUT as fallback).

Usage:
    python get_bustedph_pval.py <folder> <output.csv>

Output CSV columns:
    ogg_id,test_p,background_p,distribution_p
"""

from __future__ import print_function
import os
import sys
import json
import csv
import re

# -------- Regexes & constants --------

# matches "p = 0.012", "p-value: 1.2e-4", "p: < 0.001"
RE_PVAL = re.compile(r'(?:p[-\s]*value|p)\s*[:=]\s*(<\s*)?([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)')

def _normalize_text(s):
    s = s.replace('**', '')
    s = re.sub(r'\*+', '', s)
    return s

HEADINGS_OUT = {
    "test": re.compile(r'Likelihood\s+ratio\s+test\s+for\s+episodic\s+diversifying\s+positive\s+selection\s+on\s+test\s+branches', re.I),
    "background": re.compile(r'Likelihood\s+ratio\s+test\s+for\s+episodic\s+diversifying\s+positive\s+selection\s+on\s+background\s+branches', re.I),
    "distribution": re.compile(r'Likelihood\s+ratio\s+test\s+for\s+differences\s+in\s+distributions\s+between\s+test\s+and\s+background', re.I),
}

HEADINGS_JSON = {
    "test": re.compile(r'Likelihood\s+ratio\s+test\s+for\s+episodic\s+diversifying\s+positive\s+selection\s+on\s+test\s+branches', re.I | re.S),
    "background": re.compile(r'Likelihood\s+ratio\s+test\s+for\s+episodic\s+diversifying\s+positive\s+selection\s+on\s+background\s+branches', re.I | re.S),
    "distribution": re.compile(r'Likelihood\s+ratio\s+test\s+for\s+differences\s+in\s+distributions\s+between\s+\*?test\*?\s+and\s+\*?background\*?', re.I | re.S),
}

SECTION_KEYS = ["test results", "Test results", "BUSTED-PH test results", "BUSTED test results"]

# -------- Parsers --------

def parse_json(json_path):
    out = {"test": None, "background": None, "distribution": None}
    try:
        with open(json_path, 'r') as fh:
            data = json.load(fh)
    except Exception:
        return out

    try:
        dump_txt = json.dumps(data)
    except:
        dump_txt = str(data)

    for key, pat in HEADINGS_JSON.items():
        m = pat.search(dump_txt)
        if m and out[key] is None:
            tail = dump_txt[m.start(): m.end()+600]
            pm = RE_PVAL.search(tail)
            if pm:
                lt, val = pm.groups()
                out[key] = ("<" + val) if lt else val

    for sect in SECTION_KEYS:
        sect_data = data.get(sect)
        if isinstance(sect_data, dict):
            for k, v in sect_data.items():
                lk = k.lower() if isinstance(k, basestring) else ""
                if 'episodic' in lk and 'test' in lk and 'branches' in lk and 'background' not in lk:
                    if isinstance(v, dict):
                        p = v.get('p') or v.get('p-value') or v.get('p value')
                        if p is not None:
                            out['test'] = str(p)
                if 'episodic' in lk and 'background' in lk and 'branches' in lk:
                    if isinstance(v, dict):
                        p = v.get('p') or v.get('p-value') or v.get('p value')
                        if p is not None:
                            out['background'] = str(p)
                if 'differences' in lk and 'distributions' in lk and 'test' in lk and 'background' in lk:
                    if isinstance(v, dict):
                        p = v.get('p') or v.get('p-value') or v.get('p value')
                        if p is not None:
                            out['distribution'] = str(p)
    return out

def parse_out(out_path):
    out = {"test": None, "background": None, "distribution": None}
    try:
        with open(out_path, 'r') as fh:
            text = fh.read()
    except Exception:
        return out

    text = _normalize_text(text)

    for key, hpat in HEADINGS_OUT.items():
        m = hpat.search(text)
        if m:
            window = text[m.start(): m.start()+800]
            pm = RE_PVAL.search(window)
            if pm:
                lt, val = pm.groups()
                out[key] = ("<" + val) if lt else val

    if (out['test'] is None) or (out['background'] is None) or (out['distribution'] is None):
        chunks = re.split(r'\n\s*\n', text)
        token_map = {
            'test': 'test branches',
            'background': 'background branches',
            'distribution': 'differences in distributions',
        }
        for ch in chunks:
            lch = ch.lower()
            for key, token in token_map.items():
                if out[key] is None and token in lch:
                    pm = RE_PVAL.search(ch)
                    if pm:
                        lt, val = pm.groups()
                        out[key] = ("<" + val) if lt else val
    return out

# -------- Utilities --------

def get_ogg_id_from_filename(fname):
    base = os.path.basename(fname)
    m = re.search(r'(OG[\-_]?\d+)', base, re.I)
    if m:
        return m.group(1)
    parts = base.split('.')
    if parts:
        return parts[0]
    return base

def collect_pairs(folder):
    pairs = {}
    for root, dirs, files in os.walk(folder):
        for fn in files:
            low = fn.lower()
            if low.endswith('.json') or low.endswith('.out') or low.endswith('.txt'):
                full = os.path.join(root, fn)
                ogg = get_ogg_id_from_filename(fn)
                if ogg not in pairs:
                    pairs[ogg] = {'json': None, 'out': None}
                if low.endswith('.json'):
                    pairs[ogg]['json'] = full
                elif low.endswith('.out') or low.endswith('.txt'):
                    pairs[ogg]['out'] = full
    return pairs

def main(argv):
    if len(argv) != 3:
        print("Usage: python get_bustedph_pval.py <folder> <output.csv>")
        return 2

    folder = argv[1]
    out_csv = argv[2]

    if not os.path.isdir(folder):
        print("Error: folder not found: %s" % folder)
        return 2

    pairs = collect_pairs(folder)

    with open(out_csv, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ogg_id', 'test_p', 'background_p', 'distribution_p'])
        for ogg in sorted(pairs.keys()):
            jpath = pairs[ogg].get('json')
            opath = pairs[ogg].get('out')

            vals = {'test': None, 'background': None, 'distribution': None}
            if jpath:
                jvals = parse_json(jpath)
                if isinstance(jvals, dict):
                    vals.update({k: jvals.get(k, None) for k in vals.keys()})
            if opath:
                ovals = parse_out(opath)
                for k in vals.keys():
                    if vals[k] is None and ovals.get(k) is not None:
                        vals[k] = ovals.get(k)

            row = [ogg,
                   '' if vals['test'] is None else str(vals['test']),
                   '' if vals['background'] is None else str(vals['background']),
                   '' if vals['distribution'] is None else str(vals['distribution'])]
            writer.writerow(row)

    print("Done. Wrote: %s" % out_csv)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
