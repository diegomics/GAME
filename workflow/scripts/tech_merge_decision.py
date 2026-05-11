#!/usr/bin/env python3
"""
Analyze coverage per technology and decide merge strategy for multi-tech samples.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path



def get_coverage_from_summary(summary_path):
    """Extract coverage (mean depth) from QC summary file."""
    try:
        with open(summary_path) as f:
            for line in f:
                if line.startswith("MEAN_DEPTH"):
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        return float(parts[1])
    except Exception as e:
        print(f"WARNING: Could not read summary file {summary_path}: {e}", file=sys.stderr)
    return 0.0


def calculate_coverage(bam_paths, sample_id):
    """Get coverage for each technology from QC summary files."""
    coverage = {}
    
    for bam_path in bam_paths:
        # Extract tech from filename: {sample_id}.{tech}.merged.bam
        bam_name = Path(bam_path).name
        tech = bam_name.replace(f"{sample_id}.", "").replace(".merged.bam", "")
        
        # Build path to summary file
        bam_dir = Path(bam_path).parent
        summary_path = bam_dir / "qc" / f"{sample_id}.{tech}.summary.txt"
        
        # Get coverage from summary
        cov = get_coverage_from_summary(summary_path)
        coverage[tech] = round(cov, 2)
    
    return coverage


def make_decision(coverage, params, log_file):
    """Apply merge policy to make decision."""
    # Calculate fractions
    total_cov = sum(coverage.values())
    fractions = {}
    for tech, cov in coverage.items():
        fractions[tech] = round(cov / total_cov, 3) if total_cov > 0 else 0

    # Log coverage info
    cov_str = ", ".join([f"{t.upper()}={c:.1f}× ({fractions[t]*100:.0f}%)" 
                        for t, c in coverage.items()])
    log_file.write(f"[merge-tech] {params['sample_id']}: {cov_str}\n")
    
    # Initialize decision
    decision = {"coverage": coverage, "fractions": fractions}

    # Handle single-tech case upfront
    if len(coverage) == 1:
        tech = list(coverage.keys())[0]
        decision.update({
            "mode": "single",
            "chosen": tech,
            "reason": "Single technology sample"
        })
        log_file.write(f"[merge-tech] {params['sample_id']}: SINGLE TECH → {tech}\n")
        return decision

    # Parse priority for multi-tech cases
    priority_list = params['priority'].lower().replace(" ", "").split(">")
    mode = params['merge_mode'].lower()
    
    if mode == "off":
        # Never merge - pick by priority and coverage
        chosen = None
        for tech in priority_list:
            if tech in coverage and coverage[tech] >= params['min_cov']:
                chosen = tech
                break
        
        if not chosen and coverage:
            chosen = max(coverage, key=coverage.get)
            log_file.write(f"[merge-tech] WARNING: {chosen} coverage {coverage[chosen]}× "
                          f"< MIN_TECH_COV {params['min_cov']}\n")
        
        decision.update({
            "mode": "single",
            "chosen": chosen,
            "reason": "MERGE_TECH=off"
        })
        log_file.write(f"[merge-tech] {params['sample_id']}: mode=off → SINGLE {chosen}\n")
    
    elif mode == "auto":
        # Check if any tech has good coverage
        good_techs = [t for t in coverage if coverage[t] >= params['good_cov']]
        
        if good_techs:
            # Pick highest priority among good coverage techs
            chosen = None
            for tech in priority_list:
                if tech in good_techs:
                    chosen = tech
                    break
            if not chosen:
                chosen = good_techs[0]
            
            decision.update({
                "mode": "single",
                "chosen": chosen,
                "reason": f"Coverage {coverage[chosen]}× >= GOOD_TECH_COV"
            })
            log_file.write(f"[merge-tech] {params['sample_id']}: mode=auto → SINGLE {chosen} (good coverage)\n")
        else:
            # Merge techs that pass thresholds
            techs_to_merge = []
            dropped = []
            for tech in coverage:
                if coverage[tech] >= params['min_cov'] and fractions[tech] >= params['min_frac']:
                    techs_to_merge.append(tech)
                else:
                    dropped.append(tech)
            
            if len(techs_to_merge) > 1:
                decision.update({
                    "mode": "merge",
                    "techs_to_merge": techs_to_merge,
                    "reason": "No single tech with good coverage"
                })
                merge_str = "{" + ", ".join(t.upper() for t in techs_to_merge) + "}"
                log_file.write(f"[merge-tech] {params['sample_id']}: mode=auto → MERGE {merge_str}")
                if dropped:
                    drop_str = ", ".join([f"{t}(<{params['min_frac']*100:.0f}%)" for t in dropped])
                    log_file.write(f"; dropped {drop_str}")
                log_file.write("\n")
            elif len(techs_to_merge) == 1:
                decision.update({
                    "mode": "single",
                    "chosen": techs_to_merge[0],
                    "reason": "Only one tech passes thresholds"
                })
                log_file.write(f"[merge-tech] {params['sample_id']}: mode=auto → SINGLE {techs_to_merge[0]}\n")
            else:
                # Nothing passes - pick best available
                if coverage:
                    chosen = max(coverage, key=coverage.get)
                    decision.update({
                        "mode": "single",
                        "chosen": chosen,
                        "reason": "No techs pass thresholds, using best available"
                    })
                    log_file.write(f"[merge-tech] WARNING: using {chosen} despite low metrics\n")
    
    elif mode == "on":
        # Always merge techs that pass thresholds
        techs_to_merge = []
        for tech in coverage:
            if coverage[tech] >= params['min_cov'] and fractions[tech] >= params['min_frac']:
                techs_to_merge.append(tech)
        
        if len(techs_to_merge) > 1:
            decision.update({
                "mode": "merge",
                "techs_to_merge": techs_to_merge,
                "reason": "MERGE_TECH=on"
            })
            merge_str = "{" + ", ".join(t.upper() for t in techs_to_merge) + "}"
            log_file.write(f"[merge-tech] {params['sample_id']}: mode=on → MERGE {merge_str}\n")
        elif len(techs_to_merge) == 1:
            decision.update({
                "mode": "single",
                "chosen": techs_to_merge[0],
                "reason": "Only one tech passes thresholds"
            })
            log_file.write(f"[merge-tech] {params['sample_id']}: mode=on → SINGLE {techs_to_merge[0]}\n")
    
    return decision


def main():
    parser = argparse.ArgumentParser(description="Tech merge decision for multi-tech samples")
    parser.add_argument("--bams", nargs="*", default=[], help="Input BAM files")
    parser.add_argument("--sample-id", required=True, help="Sample ID")
    parser.add_argument("--merge-mode", default="auto", help="Merge mode: auto/on/off")
    parser.add_argument("--good-cov", type=float, default=15, help="Good coverage threshold")
    parser.add_argument("--min-cov", type=float, default=1, help="Minimum coverage threshold")
    parser.add_argument("--min-frac", type=float, default=0.05, help="Minimum fraction threshold")
    parser.add_argument("--priority", default="hifi>illumina>ont", help="Tech priority")
    parser.add_argument("--output", required=True, help="Output JSON file")
    parser.add_argument("--log", required=True, help="Log file")
    
    args = parser.parse_args()
    
    # Open log file early to capture any errors
    try:
        with open(args.log, 'w') as log_file:
            try:
                # Check if we have BAMs
                if not args.bams:
                    log_file.write(f"ERROR: No BAM files provided for {args.sample_id}\n")
                    sys.exit(1)
                
                # Get coverage from QC summary files
                coverage = calculate_coverage(args.bams, args.sample_id)
                
                # Make decision
                params = {
                    'sample_id': args.sample_id,
                    'merge_mode': args.merge_mode,
                    'good_cov': args.good_cov,
                    'min_cov': args.min_cov,
                    'min_frac': args.min_frac,
                    'priority': args.priority
                }
                
                decision = make_decision(coverage, params, log_file)
                
                # Write decision
                with open(args.output, 'w') as f:
                    json.dump(decision, f, indent=2)
                
                log_file.write(f"[merge-tech] Decision written to {args.output}\n")
                print(f"Tech merge decision written to {args.output}")
                
            except Exception as e:
                log_file.write(f"ERROR: {str(e)}\n")
                log_file.write(traceback.format_exc())
                raise
                
    except Exception as e:
        print(f"FATAL ERROR: {str(e)}", file=sys.stderr)
        print(traceback.format_exc(), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

