#!/usr/bin/env python3

# sex-ploidy inference
# by Diego De Panis, 2026
# This script is part of the GAME pipeline
# note: AI tools may have been used to improve, clean and/or comment this version of the code

"""
infer_sex_ploidy.py

Per-contig coverage-based ploidy inference for sex chromosomes when
sample_sex is "unknown".

Reads a samtools `coverage` TSV for one sample (one technology), computes
an autosomal depth baseline (median across non-sex, non-tiny contigs), and
classifies each declared sex contig as haploid / diploid / ambiguous based
on its depth-ratio to that baseline.

Output is a TSV for the user to inspect and to be consumed by the DeepVariant
calling rule to build --haploid_contigs.
"""

import argparse
import csv
import sys
from statistics import median


# -------------------------------------------------------------------------------
# THRESHOLDS
# -------------------------------------------------------------------------------
# The rational would be like this:
# A haploid contig in a heterogametic sample sits at ~0.5x of autosomal
# baseline; a diploid contig at ~1.0x. Midpoint is 0.75. The 0.65 / 0.85
# thresholds sit ±0.10 around that midpoint, giving a 0.20-wide ambiguous
# band that absorbs mapping noise or artifacts without committing a call.
HAPLOID_MAX = 0.65   # ratio ≤ this  -> haploid
DIPLOID_MIN = 0.85   # ratio ≥ this  -> diploid
# Anything in (HAPLOID_MAX, DIPLOID_MIN) is ambiguous -> not flagged.

# Below this autosomal depth, the per-contig depth distribution is too
# noisy to reliably distinguish ratio 0.5 from 1.0. We refuse to infer
# and mark everything ambiguous. 5x is the conventional floor in
# low-coverage WGS literature for sex inference.
MIN_BASELINE_DP = 5.0

# Below this breadth of coverage (% of contig with any reads), the contig is
# treated as effectively absent from this sample rather than haploid. This is
# the typical signature of a W/Y contig (assembled from the heterogametic sex)
# when the sample is the homogametic sex (e.g. a ZZ male mapped to a ZW
# reference): almost nothing maps there. Such contigs are NOT passed to
# DeepVariant's --haploid_contigs, since "one copy" is the wrong model for
# "zero copies"; the few reads that map are usually mismapping artifacts.
MIN_BREADTH_PCT = 10.0


# -------------------------------------------------------------------------------
# COVERAGE TSV PARSING
# -------------------------------------------------------------------------------
# samtools coverage output:
#   #rname  startpos  endpos  numreads  covbases  coverage  meandepth  meanbaseq  meanmapq

def read_coverage_tsv(path):
    """Read samtools coverage output. Returns list of dicts, one per contig."""
    rows = []
    with open(path) as f:
        # First line is "#rname\t..." - strip the leading '#' so csv sees it.
        header_line = f.readline().lstrip("#").strip()
        if not header_line:
            raise ValueError(f"Empty coverage file: {path}")
        fields = header_line.split("\t")
        reader = csv.DictReader(f, fieldnames=fields, delimiter="\t")
        for r in reader:
            try:
                rows.append({
                    "contig":    r["rname"],
                    "length":    int(r["endpos"]) - int(r["startpos"]) + 1,
                    "covbases":  int(r["covbases"]),
                    "coverage":  float(r["coverage"]),   # percent of contig covered
                    "meandepth": float(r["meandepth"]),
                })
            except (KeyError, ValueError) as e:
                # Skip malformed rows but keep going; surface in stderr.
                print(f"[GAME] WARN: skipping malformed coverage row "
                      f"({r}): {e}", file=sys.stderr)
    return rows


# -------------------------------------------------------------------------------
# BASELINE COMPUTATION
# -------------------------------------------------------------------------------

def compute_autosomal_baseline(rows, sex_chr_set, tiny_contig_bp):
    """
    Median meandepth across contigs that are:
      - not in the declared sex_chr set
      - at least `tiny_contig_bp` long
    Returns (baseline, n_contigs_used). baseline is 0.0 if no eligible contigs.
    """
    eligible = [
        r["meandepth"] for r in rows
        if r["contig"] not in sex_chr_set and r["length"] >= tiny_contig_bp
    ]
    if not eligible:
        return 0.0, 0
    return float(median(eligible)), len(eligible)


# -------------------------------------------------------------------------------
# CLASSIFICATION
# -------------------------------------------------------------------------------

def classify(ratio):
    """Map a depth ratio to one of: haploid, diploid, ambiguous."""
    if ratio <= HAPLOID_MAX:
        return "haploid"
    if ratio >= DIPLOID_MIN:
        return "diploid"
    return "ambiguous"


def build_inference_rows(coverage_rows, sex_chr, baseline, tiny_contig_bp,
                        baseline_too_low):
    """
    For each declared sex contig, produce one output row.
    Handles: missing-from-coverage, tiny-contig caveat, low-baseline forcing.
    """
    by_contig = {r["contig"]: r for r in coverage_rows}
    out = []

    for contig in sex_chr:
        row = by_contig.get(contig)

        if row is None:
            # Declared but not seen in the coverage file at all.
            out.append({
                "contig":      contig,
                "length":      "",
                "meandepth":   "",
                "baseline_dp": f"{baseline:.4f}" if baseline else "",
                "ratio":       "",
                "call":        "ambiguous",
                "note":        "contig not found in coverage.tsv",
            })
            continue

        if baseline_too_low:
            # Force ambiguous regardless of ratio (autosomal depth too low)
            # to distinguish 0.5x from 1.0x reliably.
            ratio = row["meandepth"] / baseline if baseline > 0 else float("nan")
            out.append({
                "contig":      contig,
                "length":      row["length"],
                "meandepth":   f"{row['meandepth']:.4f}",
                "baseline_dp": f"{baseline:.4f}",
                "ratio":       f"{ratio:.4f}" if baseline > 0 else "",
                "call":        "ambiguous",
                "note":        f"autosomal baseline {baseline:.2f}x "
                               f"below MIN_BASELINE_DP={MIN_BASELINE_DP}x; "
                               f"inference suppressed",
            })
            continue

        ratio = row["meandepth"] / baseline

        # Effectively-absent contig: almost nothing maps here. Typical of a
        # W/Y contig (from the heterogametic-sex reference) when the sample is
        # the homogametic sex. Call it "likely_absent" and it must NOT be
        # passed to --haploid_contigs (the consumer filters on call=="haploid",
        # so "likely_absent" is excluded automatically).
        if row["coverage"] < MIN_BREADTH_PCT:
            out.append({
                "contig":      contig,
                "length":      row["length"],
                "meandepth":   f"{row['meandepth']:.4f}",
                "baseline_dp": f"{baseline:.4f}",
                "ratio":       f"{ratio:.4f}",
                "call":        "likely_absent",
                "note":        f"only {row['coverage']:.1f}% of contig covered "
                               f"(< {MIN_BREADTH_PCT}%); contig appears absent "
                               f"from this sample, not flagged haploid",
            })
            continue

        call = classify(ratio)

        # Build a note for any reliability caveats. Multiple can stack.
        notes = []
        if row["length"] < tiny_contig_bp:
            notes.append(f"contig shorter than {tiny_contig_bp} bp; "
                         f"ratio less reliable")
        if row["coverage"] < 50.0:
            notes.append(f"only {row['coverage']:.1f}% of contig covered; "
                         f"ratio less reliable")

        out.append({
            "contig":      contig,
            "length":      row["length"],
            "meandepth":   f"{row['meandepth']:.4f}",
            "baseline_dp": f"{baseline:.4f}",
            "ratio":       f"{ratio:.4f}",
            "call":        call,
            "note":        "; ".join(notes),
        })

    return out


# -------------------------------------------------------------------------------
# OUTPUT
# -------------------------------------------------------------------------------

OUTPUT_FIELDS = ["contig", "length", "meandepth", "baseline_dp",
                 "ratio", "call", "note"]


def write_output(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


# -------------------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Infer sex-chromosome ploidy from coverage depth.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--coverage", required=True,
                    help="samtools coverage TSV for the chosen tech")
    ap.add_argument("--sex-chr", required=True,
                    help="Comma-separated declared sex contigs (may be empty)")
    ap.add_argument("--tiny-contig-bp", type=int, required=True,
                    help="Contigs shorter than this are excluded from baseline")
    ap.add_argument("--output", required=True,
                    help="Output sex_inference.tsv path")
    ap.add_argument("--sample-id", default="?",
                    help="Sample identifier (for logging only)")
    args = ap.parse_args()

    sex_chr = [s.strip() for s in args.sex_chr.split(",") if s.strip()]
    sex_chr_set = set(sex_chr)

    print(f"[GAME] infer_sex_ploidy: sample={args.sample_id}, "
          f"declared sex_chr={sex_chr if sex_chr else '(none)'}")

    if not sex_chr:
        # Nothing to classify. Write an empty (header-only) TSV so the
        # downstream DAG step has its input.
        write_output(args.output, [])
        print(f"[GAME] No sex contigs declared; wrote empty {args.output}")
        return 0

    coverage_rows = read_coverage_tsv(args.coverage)
    if not coverage_rows:
        print(f"[GAME] ERROR: no usable rows in {args.coverage}",
              file=sys.stderr)
        return 1

    baseline, n_used = compute_autosomal_baseline(
        coverage_rows, sex_chr_set, args.tiny_contig_bp
    )

    if n_used == 0:
        print(f"[GAME] WARN: no contigs eligible for baseline "
              f"(≥{args.tiny_contig_bp} bp and not in sex_chr). "
              f"Cannot infer ploidy.", file=sys.stderr)
        # Still emit the TSV - every declared contig will be ambiguous
        # with an appropriate note.
        baseline_too_low = True
    else:
        baseline_too_low = baseline < MIN_BASELINE_DP
        print(f"[GAME] Autosomal baseline (median over {n_used} contigs "
              f"≥{args.tiny_contig_bp} bp, excluding sex_chr): "
              f"{baseline:.3f}x")
        if baseline_too_low:
            print(f"[GAME] Baseline below MIN_BASELINE_DP="
                  f"{MIN_BASELINE_DP}x; all calls forced to ambiguous.")

    rows = build_inference_rows(
        coverage_rows, sex_chr, baseline,
        args.tiny_contig_bp, baseline_too_low
    )
    write_output(args.output, rows)

    # Brief summary to stdout for the user.
    counts = {}
    for r in rows:
        counts[r["call"]] = counts.get(r["call"], 0) + 1
    summary = ", ".join(f"{k}={v}" for k, v in sorted(counts.items()))
    print(f"[GAME] Wrote {args.output} ({summary})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
