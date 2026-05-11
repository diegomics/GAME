#!/usr/bin/env python3

# Simple markdown mini-report for BAMs
# by Diego De Panis, 2025
# This script is part of the GAME-pipeline
# See https://github.com/diegomics/GAME/tree/main/scripts/misc

# Please see the README.md for more details

import argparse
from pathlib import Path
import math

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--read-type", required=True)
    ap.add_argument("--coverage", required=True, help="samtools coverage table")
    ap.add_argument("--flagstat", required=True, help="samtools flagstat (not parsed yet, kept for future)")
    ap.add_argument("--primary", required=True)
    ap.add_argument("--primary-nodup", required=True)
    ap.add_argument("--disp-pct", type=float, default=0.95, help="central interval for dispersion (e.g. 0.95)")
    ap.add_argument("--hi-mult", type=float, default=5.0, help="high-depth outlier threshold multiplier of median depth")
    ap.add_argument("--out", required=True, help="summary .txt")
    ap.add_argument("--md-out", help="optional per-sample Markdown table")
    return ap.parse_args()

def _weighted_quantile(values, weights, p):
    """
    Weighted quantile in [0,1]. Returns the first value where cum_w >= p * total_w.
    values, weights: lists of equal length. Assumes non-negative weights.
    """
    pairs = sorted(zip(values, weights), key=lambda x: x[0])
    totw = sum(w for _, w in pairs)
    if totw <= 0:
        return float("nan")
    target = p * totw
    cum = 0.0
    for v, w in pairs:
        cum += w
        if cum >= target:
            return v
    return pairs[-1][0]

def parse_coverage_table(path):
    """
    Reads samtools coverage output.
    Expects header containing columns:
      rname, startpos, endpos, numreads, covbases, coverage, meandepth, meanbaseq, meanmapq
    Returns:
      lengths:   list of region lengths (endpos - startpos + 1)
      covbases:  list of covered bases per region (int)
      depths:    list of mean depth per region (float)
    """
    lengths, covb_list, depths = [], [], []
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if header is None:
                # header can start with '#rname' or 'rname'
                header = [c.lstrip("#") for c in line.split()]
                # sanity check minimal columns
                continue
            parts = line.split()
            if len(parts) < len(header):
                # some samtools builds use tabs; if split by spaces still fine; fallback
                parts = line.split("\t")
            row = dict(zip(header, parts))
            try:
                start = int(row["startpos"])
                end   = int(row["endpos"])
                L = end - start + 1
                covb = int(row["covbases"])
                md   = float(row["meandepth"])
            except Exception:
                # skip malformed lines
                continue
            if L <= 0:
                continue
            lengths.append(L)
            covb_list.append(covb)
            depths.append(md)
    return lengths, covb_list, depths

def safe_div(a, b):
    return (a / b) if b not in (0, 0.0) else float("nan")

# ---- Ratings used only in per-sample Markdown ----
def stars_for_breadth(b):
    if math.isnan(b): return "····"
    if b >= 99.0: return "****"
    if b >= 95.0: return "***-"
    if b >= 90.0: return "**--"
    return "*---"

def stars_for_disp95(x):
    if math.isnan(x): return "····"
    # smaller is better
    if x < 0.10: return "****"
    if x < 0.20: return "***-"
    if x < 0.30: return "**--"
    return "*---"

def stars_for_outliers(n):
    try:
        n = int(n)
    except Exception:
        return "····"
    if n == 0: return "****"
    if n <= 3: return "***-"
    if n <= 10: return "**--"
    return "*---"

def write_text_summary(p, rows):
    Path(p).parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as out:
        for k, v in rows:
            out.write(f"{k}\t{v}\n")

def write_markdown(p, d):
    """
    d is the dict keyed by the text keys below:
    SAMPLE, READ_TYPE, PRIMARY_MAPPED, PRIMARY_MAPPED_NO_DUP, GENOME_LENGTH,
    BREADTH_PCT, MEAN_DEPTH, WMD, DEPTH_DISPERSION_95, DEPTH_DISPERSION_100, HIGH_DEPTH_OUTLIERS
    """
    Path(p).parent.mkdir(parents=True, exist_ok=True)

    # parse numbers safely
    def ffloat(x, nd=3):
        try:
            x = float(x)
            if math.isnan(x): return "NA"
            fmt = f"{{:.{nd}f}}"
            return fmt.format(x)
        except Exception:
            return "NA"

    rows = []
    rows.append(("Sample", d.get("SAMPLE","NA"), "····"))
    rows.append(("Read type", d.get("READ_TYPE","NA"), "····"))
    rows.append(("Primary mapped", d.get("PRIMARY_MAPPED","NA"), "····"))
    rows.append(("Primary mapped (no dup)", d.get("PRIMARY_MAPPED_NO_DUP","NA"), "····"))
    rows.append(("Genome length", d.get("GENOME_LENGTH","NA"), "····"))

    breadth = d.get("BREADTH_PCT","NA")
    rows.append(("Breadth (%)", ffloat(breadth, 3), stars_for_breadth(float(breadth)) if breadth not in ("NA", None, "") else "····"))

    mean_dp = d.get("MEAN_DEPTH","NA")
    rows.append(("Mean depth", ffloat(mean_dp, 2), "····"))

    # WMD = weighted median depth; present as "Median depth" in the human table to avoid confusion
    wmd = d.get("WMD","NA")
    if wmd not in ("NA", None, ""):
        rows.append(("Median depth", ffloat(wmd, 2), "····"))

    disp95 = d.get("DEPTH_DISPERSION_95","NA")
    rows.append(("Depth dispersion (95%)", ffloat(disp95, 3),
                 stars_for_disp95(float(disp95)) if disp95 not in ("NA", None, "") else "····"))

    disp100 = d.get("DEPTH_DISPERSION_100","NA")
    rows.append(("Depth dispersion (100%)", ffloat(disp100, 3), "····"))

    outl = d.get("HIGH_DEPTH_OUTLIERS","0")
    rows.append(("High-depth outliers", str(outl),
                 stars_for_outliers(outl)))

    # column widths
    m_w = max(len("Metric"), max(len(r[0]) for r in rows))
    v_w = max(len("Value"),  max(len(str(r[1])) for r in rows))

    lines = []
    lines.append(f"| {'Metric':<{m_w}} | {'Value':<{v_w}} | Rating |")
    lines.append(f"| {'-'*m_w} | {'-'*v_w} | {'-'*6} |")
    for m,v,r in rows:
        lines.append(f"| {m:<{m_w}} | {v:<{v_w}} | `{r}` |")

    with open(p, "w") as out:
        out.write("\n".join(lines) + "\n")

def main():
    args = parse_args()

    # parse coverage
    lengths, covb_list, depths = parse_coverage_table(args.coverage)
    totL = sum(lengths)
    covb_total = sum(covb_list)

    breadth_pct = safe_div(covb_total, totL) * 100.0
    mean_depth  = safe_div(sum(d * L for d, L in zip(depths, lengths)), totL)

    # weighted median depth (base-weighted)
    wmd = _weighted_quantile(depths, lengths, 0.5)

    # dispersion
    alpha = (1.0 - float(args.disp_pct)) / 2.0
    q_lo = _weighted_quantile(depths, lengths, max(0.0, alpha))
    q_hi = _weighted_quantile(depths, lengths, min(1.0, 1.0 - alpha))
    disp95 = (q_hi - q_lo) / wmd if (wmd and not math.isnan(wmd) and wmd != 0.0) else float("nan")

    # full range normalized by median
    if lengths:
        dmin = min(depths)
        dmax = max(depths)
    else:
        dmin = dmax = float("nan")
    disp100 = (dmax - dmin) / wmd if (wmd and not math.isnan(wmd) and wmd != 0.0 and not math.isnan(dmax) and not math.isnan(dmin)) else float("nan")

    # high-depth outliers
    hi_thr = (args.hi_mult * wmd) if not math.isnan(wmd) else float("inf")
    hi_outliers = sum(1 for d in depths if not math.isnan(d) and d >= hi_thr)

    # prepare text rows (keep keys as before for compatibility)
    rows = [
        ("SAMPLE", args.sample),
        ("READ_TYPE", args.read_type),
        ("PRIMARY_MAPPED", str(args.primary)),
        ("PRIMARY_MAPPED_NO_DUP", str(args.primary_nodup)),
        ("GENOME_LENGTH", str(totL)),
        ("BREADTH_PCT", f"{breadth_pct:.6f}" if not math.isnan(breadth_pct) else "NA"),
        ("MEAN_DEPTH", f"{mean_depth:.6f}" if not math.isnan(mean_depth) else "NA"),
        ("WMD", f"{wmd:.6f}" if not math.isnan(wmd) else "NA"),
        ("DEPTH_DISPERSION_95", f"{disp95:.6f}" if not math.isnan(disp95) else "NA"),
        ("DEPTH_DISPERSION_100", f"{disp100:.6f}" if not math.isnan(disp100) else "NA"),
        ("HIGH_DEPTH_OUTLIERS", str(hi_outliers)),
    ]

    write_text_summary(args.out, rows)

    if args.md_out:
        d = {k: v for k, v in rows}
        write_markdown(args.md_out, d)

if __name__ == "__main__":
    main()
