#!/usr/bin/env python3

# Simple markdown mini-report aggregator for multiple BAMs 
# by Diego De Panis, 2026
# This script is part of the GAME pipeline
# note: AI tools may have been used to improve, clean and/or comment this version of the code

import argparse, pathlib

KEYS = [
    ("SAMPLE",                "Sample"),
    ("READ_TYPE",             "Type"),
    ("PRIMARY_MAPPED",        "Primary mapped"),
    ("PRIMARY_MAPPED_NO_DUP", "Primary mapped (no dup)"),
    ("GENOME_LENGTH",         "Genome length"),
    ("BREADTH_PCT",           "Breadth (%)"),
    ("MEAN_DEPTH",            "Mean depth"),
    ("WMD",                   "WMD"),
    ("DEPTH_DISPERSION_95",   "Disp (95%)"),
    ("DEPTH_DISPERSION_100",  "Disp (100%)"),
    ("HIGH_DEPTH_OUTLIERS",   "High-depth outliers"),
]

def parse_kv(path):
    d = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                d[parts[0].strip()] = parts[1].strip()
    return d

def maybe_fmt(val, nd=None, as_int=False):
    if val is None or val == "":
        return "NA"
    try:
        val_f = float(val)
        if as_int:
            return f"{int(val_f):,}"  # Format as integer with commas
        if nd is None:
            return val
        return f"{val_f:,.{nd}f}"     # Format as float with commas and specific decimals
    except Exception:
        return val

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("summaries", nargs="+")
    args = ap.parse_args()

    rows = []
    for p in args.summaries:
        kv = parse_kv(p)
        # build row with light numeric formatting
        row = {
            "Sample":                  kv.get("SAMPLE") or pathlib.Path(p).name.split(".")[0],
            "Type":                    kv.get("READ_TYPE", "NA"),
            "Primary mapped":          maybe_fmt(kv.get("PRIMARY_MAPPED"), as_int=True),
            "Primary mapped (no dup)": maybe_fmt(kv.get("PRIMARY_MAPPED_NO_DUP"), as_int=True),
            "Genome length":           maybe_fmt(kv.get("GENOME_LENGTH"), as_int=True),
            "Breadth (%)":             maybe_fmt(kv.get("BREADTH_PCT"), nd=3),
            "Mean depth":              maybe_fmt(kv.get("MEAN_DEPTH"), nd=2),
            "WMD":                     maybe_fmt(kv.get("WMD"), nd=2),
            "Disp (95%)":              maybe_fmt(kv.get("DEPTH_DISPERSION_95"), nd=3),
            "Disp (100%)":             maybe_fmt(kv.get("DEPTH_DISPERSION_100"), nd=3),
            "High-depth outliers":     maybe_fmt(kv.get("HIGH_DEPTH_OUTLIERS"), as_int=True),
        }
        rows.append(row)

    # stable sort
    rows.sort(key=lambda r: (r["Sample"], r["Type"]))

    headers = [
        "Sample","Type","Primary mapped","Primary mapped (no dup)","Genome length",
        "Breadth (%)","Mean depth","WMD","Disp (95%)","Disp (100%)","High-depth outliers"
    ]
    
    # Left-align text columns, right-align numeric columns
    align = {h: '<' if h in ("Sample", "Type") else '>' for h in headers}
    
    widths = {h: max(len(h), max((len(str(r[h])) for r in rows), default=0)) for h in headers}

    def line(cells):
        return "| " + " | ".join(f"{str(c):{align[h]}{widths[h]}}" for c, h in zip(cells, headers)) + " |\n"

    out = []
    out.append(line(headers))
    
    # Build markdown separator row with colons for right-alignment
    seps = []
    for h in headers:
        w = widths[h]
        if align[h] == '>':
            seps.append("-" * (w - 1) + ":")
        else:
            seps.append("-" * w)
            
    out.append("| " + " | ".join(seps) + " |\n")
    
    for r in rows:
        out.append(line([r[h] for h in headers]))


    # Append the explanatory legend
    legend = """
<br/>

### Column Descriptions
- **Primary mapped**: Total number of reads mapped to the reference as primary alignments.
- **Primary mapped (no dup)**: Primary mapped reads remaining after PCR and optical duplicates are removed.
- **Genome length**: Total length of the reference assembly.
- **Breadth (%)**: Percentage of the reference genome covered by at least one mapped read.
- **Mean depth**: Average sequencing coverage calculated across the entire reference assembly.
- **WMD (Weighted Median Depth)**: The median sequencing depth across all regions, weighted by region length. A much more robust metric than mean depth for assessing baseline coverage.
- **Disp (95%)**: 95% Depth Dispersion. The 95% interquantile range of depths normalized by the median depth. Measures coverage uniformity while ignoring the top and bottom 2.5% extremes.
- **Disp (100%)**: 100% Depth Dispersion. The full range of regional depths (max minus min) normalized by the median depth.
- **High-depth outliers**: Number of regions with a mean depth exceeding 5x the median depth (often suggestive of collapsed repeats, organelles or mapping artifacts).
"""
    out.append(legend)

    pathlib.Path(args.out).write_text("".join(out))

if __name__ == "__main__":
    main()