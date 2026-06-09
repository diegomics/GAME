# ===============================================================================
# GAME - Annotate, Soft-Filter (tagging), and Cleaning Rules
# by Diego De Panis, 2026
# note: AI tools may have been used to improve, clean and/or comment this version of the code
# ===============================================================================
#
#  Caller-aware tagging:
#    The filter chain adapts automatically based on which caller produced the
#    BCF (GATK or DeepVariant).  The caller is determined per-sample at DAG
#    time via _caller_for_sample() from E_variant_call.smk.
#
#  Resolution-aware tagging:
#    basepair mode:  every position is an individual record -> all tags apply,
#                    including MASKED (precise per-position annotation)
#    block mode:     ref-blocks are multi-position records -> REF_BLK is
#                    skipped, and MASKED is not applied (bcftools annotate
#                    matches by POS only, making it inaccurate on ref-blocks;
#                    use basepair resolution if you need masked-site info)
#
#  Tag inventory:
#  -------------------------------------------------------------------------
#   Tag         Callers     Description
#  -------------------------------------------------------------------------
#   NO_DATA     both        No depth info (FMT/DP and FMT/MIN_DP missing)
#   REF_N       both        Reference base is N
#   REF_BLK     both        Residual ref block (basepair mode only)
#   REF_LEAK    both*       Alt reads in hom-ref (AD for GATK, VAF for DV)
#   LOW_DP      both        FMT/DP below threshold
#   HIGH_DP     both        FMT/DP above threshold
#   LOW_GQ      both        FMT/GQ below threshold
#   LOW_MQ      gatk        INFO/MQ below threshold
#   LOW_QUAL    deep        Site QUAL below threshold
#   LOW_RGQ     both*       Hom-ref quality below threshold
#                           (RGQ for GATK, GQ for DV)
#   NON_SNP     both        Not a SNP and not a ref site
#   AB_HET      both*       Het allele balance out of range
#                           (AD ratio for GATK, VAF for DV)
#   AB_HOM      both*       Hom-alt with ref leakage
#                           (AD ratio for GATK, VAF for DV)
#   MULTI_ALT   both        Multiple alternate alleles
#   MISS_GT     both        Missing genotype with adequate depth
#   NO_GQ       both        GQ missing at a called position
#  ─────────────────────────────────────────────────────────────────────────
#   INFO tag (not a FILTER):
#   MASKED      both        Overlaps repeat/low-complexity (basepair mode only)
#  ─────────────────────────────────────────────────────────────────────────
#   * Expression differs between callers; threshold is the same.
#
#  Design principle: NO_DATA is the first-line tag for positions that
#  weren't assessed.  Other quality filters only apply to records where
#  depth info is present (FMT/DP for variant records, FMT/MIN_DP for
#  ref-block-expanded records), so uncalled positions get exactly one
#  tag (NO_DATA) instead of a pile-up of 5+ redundant tags.
#
#  Pre-existing FILTER values from the caller (e.g., DeepVariant's
#  RefCall/NoCall) are cleared at the start of tagging.  If you need
#  the caller's raw filter calls, read the gVCF.
# ===============================================================================

# -------------------------------------------------------------------------------
#  HELPER FUNCTIONS
# -------------------------------------------------------------------------------

#def _as_bool(x, default=False):
#    if isinstance(x, bool): return x
#    if x is None: return default
#    s = str(x).strip().lower()
#    return s in {"1", "true", "t", "yes", "y", "on"}

# Parse config
_RUN_TAGGING = _as_bool(config.get("RUN_TAGGING", True))
_KEEP_RAW_BCF = _as_bool(config.get("KEEP_RAW_BCF", False))
_TINY_BP = int(config.get("TINY_CONTIG_BP", 1_000_000))

# Filter thresholds (only used if RUN_TAGGING is on)
_MIN_DEPTH_CFG = str(config.get("MIN_DEPTH", "8")).lower()
_MAX_DEPTH_CFG = str(config.get("MAX_DEPTH", "auto")).lower()
_MIN_GQ = int(config.get("MIN_GQ", 10))
_MIN_MQ = int(config.get("MIN_MQ", 40))
_MIN_QUAL = int(config.get("MIN_QUAL", 1))
_AB_HET_MIN = float(config.get("AB_HET_MIN", 0.25))
_AB_HET_MAX = float(config.get("AB_HET_MAX", 0.75))
_AB_HOM_REF_MAX = float(config.get("AB_HOM_REF_MAX", 0.10))

# Resolution (from E_variant_call.smk, but re-read here for safety)
_RESOLUTION = str(config.get("RESOLUTION", "block")).strip().lower()
_BP_MODE = (_RESOLUTION == "basepair")
_BP_TAG = ".bp" if _BP_MODE else ""

# Joint genotyping tagging thresholds (only used if JOINT_GENO + RUN_TAGGING are on)
_MIN_AQ = int(config.get("MIN_AQ", 10))
_MIN_JOINT_QUAL = int(config.get("MIN_JOINT_QUAL", 1))
_MAX_MISSING = float(config.get("MAX_MISSING", 0.5))
_MIN_AC_CFG = str(config.get("MIN_AC", "auto")).strip().lower()
_MIN_AC_FRAC = float(config.get("MIN_AC_FRAC", 0.01))
_KEEP_RAW_JOINT = _as_bool(config.get("KEEP_RAW_JOINT", False))

# Relatedness check (F03) - KING-robust kinship via plink2.
# Only used when joint genotyping is on (F03's input is the joint BCF).
_KING_CUTOFF = float(config.get("KING_CUTOFF", 0.0884))
_KING_MAF = float(config.get("KING_MAF", 0.05))
_KING_LD_WINDOW = int(config.get("KING_LD_WINDOW", 50))
_KING_LD_STEP = int(config.get("KING_LD_STEP", 5))
_KING_LD_R2 = float(config.get("KING_LD_R2", 0.2))
# Minimum sample size to perform LD pruning.  Below this, plink2's r²
# estimates are statistically unreliable (plink2's own --bad-ld threshold
# is 50).  When N<KING_LD_MIN_SAMPLES, F03 runs KING without LD pruning
# (KING-robust is unaffected by LD per its large-sample theory) and F05
# skips LD-prune annotation entirely.
_KING_LD_MIN_SAMPLES = int(config.get("KING_LD_MIN_SAMPLES", 50))


def _per_sample_stats_for_assembly(w):
    """Collect all per-sample F01 tagging_stats.md files for an assembly.
    Mirrors _joint_gvcfs_for_assembly: only includes samples that have reads.
    Used by F04 to look up PASS_SNPS_PCT for related samples when deciding
    which one of a related pair to drop.
    """
    sp, asm = w.species, w.assembly
    paths = []
    try:
        asm_data = samples_config["sp_name"][sp]["asm_id"][asm]
        if "sample_id" not in asm_data:
            return paths
        for sid in asm_data["sample_id"]:
            if sid in ("None", None):
                continue
            sample_data = asm_data["sample_id"][sid]
            has_reads = False
            for rt_key in sample_data.get("read_type", {}):
                if rt_key not in ("None", None):
                    _, _, reads = _get_sample_node(sp, asm, sid, rt_key)
                    if reads:
                        has_reads = True
                        break
            if has_reads:
                paths.append(os.path.join(
                    config["OUT_FOLDER"], "GAME_results", sp, asm,
                    "samples", sid, "VCFs", f"{sid}.tagging_stats.md"
                ))
    except (KeyError, TypeError, AttributeError):
        pass
    return paths


def fv_get_mask_bed(sp, asm):
    """Get masking BED file path based on the MASKING control-panel variable.

    Returns the expected path of the masking BED for this assembly, or an
    empty string if masking is disabled. The path returned matches whichever
    rule will produce it:
      - MASKING=on   -> {asm}.mask_from_rm.bed.gz   (from RepeatMasker)
      - MASKING=auto -> {asm}.mask_from_file.bed.gz (sibling BED found)
                     or {asm}.mask_from_ref.bed.gz  (extracted from softmask)
      - MASKING=off  -> ""
    """
    mode = str(config.get("MASKING", "off")).strip().lower()

    if mode in ("", "off", "false", "no", "0", "none"):
        return ""

    base = os.path.join(config["OUT_FOLDER"], "GAME_results", sp, asm, "masking")

    if mode in ("on", "true", "yes", "1"):
        return os.path.join(base, f"{asm}.mask_from_rm.bed.gz")

    if mode == "auto":
        if get_sibling_bed_path(sp, asm) is not None:
            return os.path.join(base, f"{asm}.mask_from_file.bed.gz")
        return os.path.join(base, f"{asm}.mask_from_ref.bed.gz")

    # Unknown value — match the driver's behavior and fail loudly
    raise ValueError(
        f"[GAME] Invalid MASKING value: {config.get('MASKING')!r}. "
        f"Expected one of: off, on, auto (or their aliases)."
    )


def fv_qc_summaries_for_sample(w):
    """Per-tech QC summary files needed when MIN/MAX_DEPTH=auto.

    Returns [] when neither threshold is auto (no dependency needed).
    Otherwise reads merge_decision.json to know which tech summaries to
    require. If merge_decision.json doesn't exist yet at DAG time, returns
    [] and the rule's shell will fail loudly when it can't find them — this
    is fine because merge_decision.json is produced upstream and will exist
    by the time F01 actually runs.
    """
    if _MIN_DEPTH_CFG != "auto" and _MAX_DEPTH_CFG != "auto":
        return []

    merge_decision = os.path.join(
        config["OUT_FOLDER"], "GAME_results", w.species, w.assembly,
        "samples", w.sample_id, "BAMs", "merge_decision.json"
    )
    if not os.path.exists(merge_decision):
        return []

    try:
        with open(merge_decision) as fh:
            d = json.load(fh)
    except Exception:
        return []

    mode = d.get("mode", "")
    if mode == "merge":
        techs = d.get("techs_to_merge", [])
    elif mode == "single":
        techs = [d.get("chosen", "")]
    else:
        techs = list(d.get("coverage", {}).keys())

    qc_dir = os.path.join(
        config["OUT_FOLDER"], "GAME_results", w.species, w.assembly,
        "samples", w.sample_id, "BAMs", "qc"
    )
    return [
        os.path.join(qc_dir, f"{w.sample_id}.{t}.summary.txt")
        for t in techs if t
    ]


# -------------------------------------------------------------------------------
#  ANNOTATE WITH INFO TAGS, THEN APPLY FILTER TAGS
# ===============================================================================

rule F01_tag_variants:
    """
    Add QC tags to variants (caller-aware and resolution-aware).
    Step 1: INFO tags (TYPE; MASKED only in basepair mode)
    Step 2: FILTER tags - shared filters applied first, then caller-specific
    """
    input:
        bcf=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "VCFs", "{sample_id}" + _BP_TAG + ".raw.bcf"
        ),
        csi=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "VCFs", "{sample_id}" + _BP_TAG + ".raw.bcf.csi"
        ),
        qc_summaries=fv_qc_summaries_for_sample,
        mask=lambda w: fv_get_mask_bed(w.species, w.assembly) or []
    output:
        bcf=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "VCFs", "{sample_id}" + _BP_TAG + ".tagged.bcf"
        ),
        csi=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "VCFs", "{sample_id}" + _BP_TAG + ".tagged.bcf.csi"
        ),
        stats=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "VCFs", "{sample_id}.tagging_stats.md"
        )
    params:
        mask_bed=lambda w: fv_get_mask_bed(w.species, w.assembly) or "none",
        keep_raw=_KEEP_RAW_BCF,
        min_dp=_MIN_DEPTH_CFG,
        max_dp=_MAX_DEPTH_CFG,
        min_gq=_MIN_GQ,
        min_mq=_MIN_MQ,
        min_qual=_MIN_QUAL,
        ab_het_min=_AB_HET_MIN,
        ab_het_max=_AB_HET_MAX,
        ab_hom_ref_max=_AB_HOM_REF_MAX,
        bp_mode=_BP_MODE,
        caller=lambda w: _caller_for_sample(w.species, w.assembly, w.sample_id),
        qc_dir=lambda w: os.path.join(
            config["OUT_FOLDER"], "GAME_results", w.species, w.assembly,
            "samples", w.sample_id, "BAMs", "qc"
        ),
        merge_decision=lambda w: os.path.join(
            config["OUT_FOLDER"], "GAME_results", w.species, w.assembly,
            "samples", w.sample_id, "BAMs", "merge_decision.json"
        )
    threads: cpu_func("bcftools_concat")
    resources:
        mem_mb=mem_func("bcftools_concat"),
        runtime=time_func("bcftools_concat")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "F01_tag_variants.{sample_id}.log"
        )
    shell:
        r'''
        set -euo pipefail
        export LC_NUMERIC=en_US.UTF-8
        mkdir -p $(dirname {log})
        
        exec > "{log}" 2>&1

        CALLER="{params.caller}"
        BP_MODE="{params.bp_mode}"
        echo "[GAME] Starting QC tagging for {wildcards.sample_id}"
        echo "[GAME] Caller: $CALLER  |  Resolution: {_RESOLUTION}"
        

        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 20)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_tag_{wildcards.species}_{wildcards.assembly}_{wildcards.sample_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        
        CURRENT_BCF="{input.bcf}"
        
        # ============================================
        # STEP 1: Clear pre-existing FILTER values
        # ============================================
        # DeepVariant writes its own FILTER values (RefCall, NoCall) that
        # would clutter our tag namespace.  We own the FILTER column after
        # F01 - any user who needs DV's raw calls can read the gVCF.
        echo "[GAME] Clearing pre-existing FILTER values"
        bcftools annotate --threads {threads} -x FILTER \
            "$CURRENT_BCF" -Ob -o cleared.bcf
        CURRENT_BCF="cleared.bcf"
        
        # ============================================
        # STEP 2: Add INFO tags (TYPE, and MASKED if bp mode)
        # ============================================
        
        # --- MASKED flag (basepair mode only) ---
        # In basepair mode every position is its own record, so bcftools
        # annotate matches precisely.  In block mode, POS-only matching is
        # inaccurate on ref-blocks → we skip masking entirely.  Users who
        # need masked-site information should use RESOLUTION: basepair.
        if [[ "$BP_MODE" == "True" ]] && [[ "{params.mask_bed}" != "none" ]] && [[ -f "{params.mask_bed}" ]]; then
            echo "[GAME] Adding MASKED flag (basepair mode)"
            
            ( [[ "{params.mask_bed}" =~ \.gz$ ]] && gzip -dc "{params.mask_bed}" || cat "{params.mask_bed}" ) \
              | awk 'BEGIN{{OFS="\t"}} NF>=3{{print $1,$2,$3,1}}' \
              | LC_ALL=C sort -k1,1 -k2,2n \
              | bgzip -c > mask4.bed.gz && tabix -p bed mask4.bed.gz

            echo '##INFO=<ID=MASKED,Number=0,Type=Flag,Description="Overlaps repeat/low-complexity regions">' > mask.hdr
            
            bcftools annotate \
                -a mask4.bed.gz \
                -c CHROM,FROM,TO,INFO/MASKED \
                -h mask.hdr \
                -Ob -o masked.bcf \
                "$CURRENT_BCF"
            
            CURRENT_BCF="masked.bcf"
        elif [[ "$BP_MODE" != "True" ]] && [[ "{params.mask_bed}" != "none" ]]; then
            echo "[GAME] Skipping MASKED flag (block mode - use basepair resolution for masking)"
        fi
        
        # --- TYPE tag ---
        echo "[GAME] Adding TYPE tag"
        bcftools +fill-tags --threads {threads} \
            "$CURRENT_BCF" \
            -Ob -o info_tagged.bcf \
            -- -t TYPE
        
        CURRENT_BCF="info_tagged.bcf"
        
        # ============================================
        # STEP 3: Resolve depth thresholds
        # ============================================

        MIN_DP_CFG="{params.min_dp}"
        MAX_DP_CFG="{params.max_dp}"
        BASELINE_DP="0.0"

        if [[ "$MIN_DP_CFG" == "auto" || "$MAX_DP_CFG" == "auto" ]]; then
            if [[ ! -f "{params.merge_decision}" ]]; then
                echo "[GAME] ERROR: auto depth thresholds requested but merge_decision.json not found:" >&2
                echo "[GAME]        {params.merge_decision}" >&2
                exit 1
            fi

            TECHS=$(python - <<'PY'
import json
with open(r"""{params.merge_decision}""") as f:
    d = json.load(f)
mode = d.get("mode", "")
if mode == "merge":
    techs = d.get("techs_to_merge", [])
elif mode == "single":
    techs = [d.get("chosen", "")]
else:
    techs = list(d.get("coverage", {{}}).keys())
print(" ".join(t for t in techs if t))
PY
)

            if [[ -z "$TECHS" ]]; then
                echo "[GAME] ERROR: no techs resolved from merge_decision.json" >&2
                exit 1
            fi

            echo "[GAME] Reading WMD (length-weighted median depth) from techs: $TECHS"

            # note on summing medians across techs:
            # For a merged BAM, the true baseline we want is the median of
            # (tech_A_depth + tech_B_depth + ...) at each position. We
            # approximate that as sum(median(tech_i)). This is exact only
            # for symmetric distributions, but our per-tech depth
            # distributions are tight (DEPTH_DISPERSION_95 typically < 0.5),
            # so the approximation error is small relative to the 3x buffer
            # we apply when deriving MAX_DP. Crucially, sum(medians) is
            # still strictly more robust than sum(means) when one tech has
            # a heavy upper tail (collapsed repeats, organelle contamination),
            # which is the failure mode we're protecting against.
            
            for t in $TECHS; do
                f="{params.qc_dir}/{wildcards.sample_id}.$t.summary.txt"
                if [[ ! -s "$f" ]]; then
                    echo "[GAME] ERROR: expected QC summary not found or empty: $f" >&2
                    exit 1
                fi
                v=$(awk -F'\t' '$1=="WMD"{{print $2}}' "$f")
                if [[ -z "$v" ]]; then
                    echo "[GAME] ERROR: WMD field missing in $f" >&2
                    exit 1
                fi
                BASELINE_DP=$(awk -v a="$BASELINE_DP" -v b="$v" 'BEGIN{{printf "%.4f", a+b}}')
            done

            echo "[GAME] Baseline depth (sum of per-tech WMD): $BASELINE_DP"
        fi

        # Convert auto to concrete thresholds
        if [[ "$MIN_DP_CFG" == "auto" ]]; then
            MIN_DP=$(awk -v m="$BASELINE_DP" 'BEGIN{{v=int(m/3); if(v<1)v=1; print v}}')
        else
            MIN_DP="$MIN_DP_CFG"
        fi
        if [[ "$MAX_DP_CFG" == "auto" ]]; then
            MAX_DP=$(awk -v m="$BASELINE_DP" 'BEGIN{{v=3*m; c=(int(v)==v)?v:int(v)+1; if(c<1)c=1; print c}}')
        else
            MAX_DP="$MAX_DP_CFG"
        fi

        echo "[GAME] Thresholds: MIN_DP=$MIN_DP, MAX_DP=$MAX_DP, MIN_GQ={params.min_gq}"

        
        # ============================================
        # STEP 4: Shared filters (both callers)
        # ============================================
        # Philosophy: NO_DATA is for positions that weren't assessed
        # (no DP and no MIN_DP - i.e. no depth info at all). Other quality
        # filters only apply to records where data IS present. This prevents
        # every uncovered position in basepair mode from piling up 5+ tags.
        #
        # DP vs MIN_DP: DV's gvcf2vcf preserves the gVCF convention where
        # ref-block-expanded records carry MIN_DP (block minimum) instead
        # of DP. Both signal valid coverage - we treat them equivalently.
        # GATK's GenotypeGVCFs output uses DP everywhere, so the MIN_DP
        # clauses simply never fire for GATK records (no behaviour change).
        # We define HAS_DEPTH = (DP present OR MIN_DP present) to keep the
        # filter logic compact.
        
        # ── Ensure FMT/MIN_DP is declared in the header ─────────────────
        # The shared filters below reference FMT/MIN_DP so DeepVariant's
        # ref-block-expanded records (which carry MIN_DP instead of DP) flow
        # through the same compact chain as GATK. bcftools validates EVERY
        # tag in a filter expression against the VCF header at parse time
        # (filters_init), before it reads a single record — so a missing
        # *header declaration* is fatal even when no record would ever match
        # the clause. GATK's GenotypeGVCFs output never declares MIN_DP,
        # which is why this rule failed with:
        #   [filter.c] Error: the tag "MIN_DP" is not defined in the VCF header
        # Declaring the FORMAT line (without adding the field to any record)
        # makes the expression parseable; for records lacking the field the
        # MIN_DP clauses evaluate to missing and never fire — exactly the
        # no-op behaviour the comment above assumes.
        bcftools view -h "$CURRENT_BCF" > pre_filter_header.txt
        if ! grep -q '^##FORMAT=<ID=MIN_DP,' pre_filter_header.txt; then
            echo "[GAME] FMT/MIN_DP absent from header (expected for GATK) - declaring it for filter compatibility"
            echo '##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP within a gVCF reference block; header declared by GAME so the shared depth filters parse on callers that omit this field">' > min_dp.hdr
            bcftools annotate --threads {threads} -h min_dp.hdr \
                "$CURRENT_BCF" -Ob -o min_dp_hdr.bcf
            CURRENT_BCF="min_dp_hdr.bcf"
        fi
        
        echo "[GAME] Applying shared filters"
        
        bcftools filter --threads {threads} -m + -s REF_N    -e 'REF="N"' -Ou "$CURRENT_BCF" | \
        bcftools filter -m + -s NO_DATA  -e 'FMT/DP="." && FMT/MIN_DP="."' -Ou | \
        bcftools filter -m + -s LOW_DP   -e '(FMT/DP!="." && FMT/DP<'$MIN_DP') || (FMT/MIN_DP!="." && FMT/MIN_DP<'$MIN_DP')' -Ou | \
        bcftools filter -m + -s HIGH_DP  -e '(FMT/DP!="." && FMT/DP>'$MAX_DP') || (FMT/MIN_DP!="." && FMT/MIN_DP>'$MAX_DP')' -Ou | \
        bcftools filter -m + -s NO_GQ    -e '(FMT/DP!="." || FMT/MIN_DP!=".") && FMT/GQ="." && INFO/TYPE!="REF"' -Ou | \
        bcftools filter -m + -s LOW_GQ   -e '(FMT/DP!="." || FMT/MIN_DP!=".") && FMT/GQ!="." && FMT/GQ<{params.min_gq}' -Ou | \
        bcftools filter -m + -s NON_SNP  -e '(FMT/DP!="." || FMT/MIN_DP!=".") && INFO/TYPE!="SNP" && INFO/TYPE!="REF"' -Ou | \
        bcftools filter -m + -s MULTI_ALT -e '(FMT/DP!="." || FMT/MIN_DP!=".") && N_ALT>1' -Ou | \
        bcftools filter --threads {threads} -m + -s MISS_GT \
            -e 'GT="mis" && ((FMT/DP!="." && FMT/DP>='$MIN_DP') || (FMT/MIN_DP!="." && FMT/MIN_DP>='$MIN_DP'))' \
            -Ob -o "shared.bcf"
        
        # ============================================
        # STEP 5: Caller-specific filters
        # ============================================
        # All caller-specific filters also require FMT/DP to be present,
        # so no-data records don't accumulate extra tags.
        
        if [[ "$CALLER" == "deep" ]]; then
            # ── DeepVariant filters ──────────────────────────────────
            echo "[GAME] Applying DeepVariant-specific filters"
            echo "[GAME]   LOW_QUAL threshold: {params.min_qual}"
            
            # LOW_QUAL: DV's QUAL is a calibrated NN confidence (not useful for GATK)
            # The QUAL>=0 guard catches missing QUAL (which bcftools treats as negative/missing)
            bcftools filter --threads {threads} -m + -s LOW_QUAL \
                -e 'FMT/DP!="." && QUAL>=0 && QUAL<{params.min_qual} && INFO/TYPE!="REF"' -Ou "shared.bcf" | \
            \
            # LOW_RGQ: DV has no RGQ field - use GQ on hom-ref sites.
            # Hom-ref calls come from both variant records (FMT/DP) and
            # ref-block expanded records (FMT/MIN_DP) → check both.
            bcftools filter -m + -s LOW_RGQ \
                -e '(FMT/DP!="." || FMT/MIN_DP!=".") && (GT="0/0" || GT="0|0") && FMT/GQ!="." && FMT/GQ<{params.min_gq}' -Ou | \
            \
            # REF_LEAK: DV emits FMT/VAF; on hom-ref it should be ~0
            bcftools filter -m + -s REF_LEAK \
                -e '(GT="0/0" || GT="0|0") && FMT/DP!="." && FMT/DP>=1 && FMT/VAF!="." && FMT/VAF>{params.ab_hom_ref_max}' -Ou | \
            \
            # AB_HET: use VAF directly (het should be in [min, max] range)
            bcftools filter -m + -s AB_HET \
                -e 'N_ALT==1 && (GT="0/1" || GT="1/0" || GT="0|1" || GT="1|0") && FMT/DP!="." && FMT/DP>=1 && FMT/VAF!="." && (FMT/VAF<{params.ab_het_min} || FMT/VAF>{params.ab_het_max})' -Ou | \
            \
            # AB_HOM: hom-alt should have VAF close to 1.0
            # Flag if ref fraction > threshold, i.e. VAF < (1 - threshold)
            bcftools filter --threads {threads} -m + -s AB_HOM \
                -e 'N_ALT==1 && (GT="1/1" || GT="1|1") && FMT/DP!="." && FMT/DP>=1 && FMT/VAF!="." && FMT/VAF<(1-{params.ab_hom_ref_max})' \
                -Ob -o "caller.bcf"
        
        else
            # ── GATK filters ────────────────────────────────────────
            echo "[GAME] Applying GATK-specific filters"
            echo "[GAME]   MIN_MQ threshold: {params.min_mq}"
            
            # LOW_MQ: GATK emits INFO/MQ (average mapping quality)
            # Note: numeric INFO fields can't be tested with !="."; bcftools
            # handles missing values correctly in arithmetic expressions.
            bcftools filter --threads {threads} -m + -s LOW_MQ \
                -e 'FMT/DP!="." && INFO/MQ<{params.min_mq} && INFO/TYPE!="REF"' -Ou "shared.bcf" | \
            \
            # LOW_RGQ: GATK emits FMT/RGQ for hom-ref sites
            bcftools filter -m + -s LOW_RGQ \
                -e 'FMT/DP!="." && (GT="0/0" || GT="0|0") && (FMT/RGQ<{params.min_gq} || FMT/RGQ=".")' -Ou | \
            \
            # REF_LEAK: GATK AD array - AD[0:1] is alt depth on hom-ref
            bcftools filter -m + -s REF_LEAK \
                -e '(GT="0/0" || GT="0|0") && FMT/DP!="." && FMT/DP>=1 && FMT/AD[0:1]!="." && FMT/AD[0:1]>=2 && FMT/AD[0:1]/FMT/DP>{params.ab_hom_ref_max}' -Ou | \
            \
            # AB_HET: GATK AD ratio - AD[0:1]/DP is alt fraction
            bcftools filter -m + -s AB_HET \
                -e 'N_ALT==1 && (GT="0/1" || GT="1/0" || GT="0|1" || GT="1|0") && FMT/DP!="." && FMT/DP>=1 && FMT/AD[0:1]!="." && (FMT/AD[0:1]/FMT/DP<{params.ab_het_min} || FMT/AD[0:1]/FMT/DP>{params.ab_het_max})' -Ou | \
            \
            # AB_HOM: hom-alt with ref allele leakage - AD[0:0]/DP is ref fraction
            bcftools filter --threads {threads} -m + -s AB_HOM \
                -e 'N_ALT==1 && (GT="1/1" || GT="1|1") && FMT/DP!="." && FMT/DP>=1 && FMT/AD[0:0]!="." && FMT/AD[0:0]/FMT/DP>{params.ab_hom_ref_max}' \
                -Ob -o "caller.bcf"
        fi
        
        # ============================================
        # STEP 6: Resolution-specific filters
        # ============================================
        
        if [[ "$BP_MODE" == "True" ]]; then
            # In basepair mode, flag residual ref-blocks (gVCF artefacts
            # that survived expansion - shouldn't happen, but safety net)
            echo "[GAME] Applying basepair-mode filter: REF_BLK"
            bcftools filter --threads {threads} -m + -s REF_BLK \
                -e 'INFO/TYPE="REF" && strlen(REF)>1' \
                -Ob -o "output.bcf" "caller.bcf"
        else
            echo "[GAME] Block mode - skipping REF_BLK (ref-blocks are expected)"
            mv "caller.bcf" "output.bcf"
        fi
        
        bcftools index -f --csi --threads {threads} "output.bcf"
        
        # Copy BCF to final location
        echo "[GAME] Copying BCF to final location..."
        mkdir -p "$(dirname {output.bcf})"
        cp "output.bcf" "{output.bcf}"
        cp "output.bcf.csi" "{output.csi}"
        
        # ============================================
        # STEP 7: Generate tagging statistics
        # ============================================
        {{
            echo "# Variant Tagging Statistics"
            echo ""
            echo "**Sample:** {wildcards.sample_id}  "
            echo "**Generated:** $(date)  "
            echo "**Species:** {wildcards.species}  "
            echo "**Reference:** {wildcards.assembly}  "
            if [[ "$CALLER" == "deep" ]]; then
                CALLER_DISPLAY="DeepVariant"
            else
                CALLER_DISPLAY="GATK HaplotypeCaller"
            fi
            echo "**Caller:** $CALLER_DISPLAY  "
            echo "**Resolution:** {_RESOLUTION}"
            echo ""
            echo "---"
            echo ""
            
            if [[ "$BP_MODE" != "True" ]]; then
                echo "> **Note:** Resolution is *block*. Reference sites are stored as"
                echo "> multi-position blocks. MASKED annotation is not applied in block"
                echo "> mode (use basepair resolution for masked-site information)."
                echo ""
                echo "---"
                echo ""
            fi


            echo "## Quality Thresholds"
            echo ""
            
            # Write raw data to a temp file
            tmp_stats_1=$(mktemp)
            {{
                printf "Parameter\tValue\n"
                printf "MIN_DP\t%s\n" "$MIN_DP"
                printf "MAX_DP\t%s\n" "$MAX_DP"
                printf "MIN_GQ\t%s\n" "{params.min_gq}"
                if [[ "$CALLER" == "gatk" ]]; then
                    printf "MIN_MQ\t%s\n" "{params.min_mq}"
                fi
                if [[ "$CALLER" == "deep" ]]; then
                    printf "MIN_QUAL\t%s\n" "{params.min_qual}"
                fi
                printf "AB_HET_MIN\t%s\n" "{params.ab_het_min}"
                printf "AB_HET_MAX\t%s\n" "{params.ab_het_max}"
                printf "AB_HOM_REF_MAX\t%s\n" "{params.ab_hom_ref_max}"
            }} > "$tmp_stats_1"

            # Dynamically size and print the Markdown table
            awk -F'\t' '
                BEGIN {{ max_c1=0; max_c2=0 }}
                NR==FNR {{
                    if(length($1)>max_c1) max_c1=length($1)
                    if(length($2)>max_c2) max_c2=length($2)
                    next
                }}
                FNR==1 {{
                    printf "| %-*s | %-*s |\n", max_c1, $1, max_c2, $2
                    s1=sprintf("%*s", max_c1, ""); gsub(/ /, "-", s1)
                    s2=sprintf("%*s", max_c2, ""); gsub(/ /, "-", s2)
                    printf "|-%s-|-%s-|\n", s1, s2
                }}
                FNR>1 {{
                    printf "| %-*s | %-*s |\n", max_c1, $1, max_c2, $2
                }}
            ' "$tmp_stats_1" "$tmp_stats_1"
            rm -f "$tmp_stats_1"
            
            echo ""

            echo "## Summary Counts"
            echo ""            

            TOTAL=$(bcftools view -H "{output.bcf}" | wc -l)
            PASS_SITES=$(bcftools view -H -f 'PASS,.' "{output.bcf}" | wc -l)
            PASS_SNPS=$(bcftools view -H -f 'PASS,.' -i 'INFO/TYPE="SNP"' "{output.bcf}" | wc -l)
            PASS_PCT=$(awk -v a=$PASS_SITES -v t=$TOTAL 'BEGIN{{if(t>0) printf "%.2f", a*100/t; else print "0"}}')
            PASS_SNPS_PCT=$(awk -v a=$PASS_SNPS -v t=$TOTAL 'BEGIN{{if(t>0) printf "%.4f", a*100/t; else print "0"}}')

            # Count MASKED sites (only present in basepair mode with mask BED)
            if [[ "$BP_MODE" == "True" ]] && [[ "{params.mask_bed}" != "none" ]] && [[ -f "{params.mask_bed}" ]]; then
                MASKED=$(bcftools query -f '%INFO/MASKED\n' "{output.bcf}" | grep -c "1" || true)
                MASKED_PCT=$(awk -v m=$MASKED -v t=$TOTAL 'BEGIN{{if(t>0) printf "%.2f", m*100/t; else print "0"}}')
                PASS_MASKED=$(bcftools view -H -f 'PASS,.' -i 'INFO/MASKED=1' "{output.bcf}" | wc -l)
                PASS_MASKED_PCT=$(awk -v m=$PASS_MASKED -v t=$TOTAL 'BEGIN{{if(t>0) printf "%.2f", m*100/t; else print "0"}}')
                PASS_SNPS_MASKED=$(bcftools view -H -f 'PASS,.' -i 'INFO/TYPE="SNP" && INFO/MASKED=1' "{output.bcf}" | wc -l)
                PASS_SNPS_MASKED_PCT=$(awk -v m=$PASS_SNPS_MASKED -v t=$TOTAL 'BEGIN{{if(t>0) printf "%.4f", m*100/t; else print "0"}}')
            fi

            # Write raw data to a temp file
            tmp_stats_2=$(mktemp)
            {{
                printf "Category\tCount\tPercentage\n"
                printf "**Total Sites**\t%s\t100.00%%\n" "$(printf "%'d" $TOTAL)"
                if [[ "$BP_MODE" == "True" ]] && [[ "{params.mask_bed}" != "none" ]] && [[ -f "{params.mask_bed}" ]]; then
                    printf "**Total Masked**\t%s\t%s%%\n" "$(printf "%'d" $MASKED)" "$MASKED_PCT"
                elif [[ "$BP_MODE" != "True" ]]; then
                    printf "**Total Masked**\tn/a\tn/a (block mode)\n"
                else
                    printf "**Total Masked**\tn/a\tn/a (no mask BED)\n"
                fi
                printf "**PASS Sites**\t%s\t%s%%\n" "$(printf "%'d" $PASS_SITES)" "$PASS_PCT"
                if [[ "$BP_MODE" == "True" ]] && [[ "{params.mask_bed}" != "none" ]] && [[ -f "{params.mask_bed}" ]]; then
                    printf "**PASS Masked**\t%s\t%s%%\n" "$(printf "%'d" $PASS_MASKED)" "$PASS_MASKED_PCT"
                fi
                printf "**PASS SNPs**\t%s\t%s%%\n" "$(printf "%'d" $PASS_SNPS)" "$PASS_SNPS_PCT"
                if [[ "$BP_MODE" == "True" ]] && [[ "{params.mask_bed}" != "none" ]] && [[ -f "{params.mask_bed}" ]]; then
                    printf "**PASS SNPs Masked**\t%s\t%s%%\n" "$(printf "%'d" $PASS_SNPS_MASKED)" "$PASS_SNPS_MASKED_PCT"
                fi
            }} > "$tmp_stats_2"

            # Dynamically size and print the Markdown table
            awk -F'\t' '
                BEGIN {{ max_c1=0; max_c2=0; max_c3=0 }}
                NR==FNR {{
                    if(length($1)>max_c1) max_c1=length($1)
                    if(length($2)>max_c2) max_c2=length($2)
                    if(length($3)>max_c3) max_c3=length($3)
                    next
                }}
                FNR==1 {{
                    # Print header (left, right, right)
                    printf "| %-*s | %*s | %*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                    s1=sprintf("%*s", max_c1, ""); gsub(/ /, "-", s1)
                    s2=sprintf("%*s", max_c2, ""); gsub(/ /, "-", s2)
                    s3=sprintf("%*s", max_c3, ""); gsub(/ /, "-", s3)
                    printf "|-%s-|-%s:|-%s:|\n", s1, s2, s3
                }}
                FNR>1 {{
                    # Print body (left, right, right)
                    printf "| %-*s | %*s | %*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                }}
            ' "$tmp_stats_2" "$tmp_stats_2"
            rm -f "$tmp_stats_2"
            
            echo ""
            
            echo "## Filter Tag Counts"
            echo ""
            
            # Create associative array for descriptions
            declare -A DESC
            DESC[PASS]="Site passed all quality filters"
            DESC[.]="No filter applied (reference sites)"
            DESC[REF_N]="Reference base is N"
            DESC[REF_BLK]="Reference block (gVCF artefact, basepair mode only)"
            DESC[REF_LEAK]="Alt-allele signal in hom-ref call"
            DESC[LOW_DP]="Depth below minimum threshold"
            DESC[HIGH_DP]="Depth above maximum threshold"
            DESC[LOW_GQ]="Genotype quality below threshold"
            DESC[LOW_MQ]="Mapping quality below threshold (GATK only)"
            DESC[LOW_QUAL]="Site QUAL below threshold (DeepVariant only)"
            DESC[LOW_RGQ]="Reference genotype quality below threshold"
            DESC[NON_SNP]="Not a SNP (indel or other variant type)"
            DESC[AB_HET]="Allele balance out of range for heterozygote"
            DESC[AB_HOM]="Ref-allele leakage in homozygous alternate"
            DESC[MULTI_ALT]="Multiple alternate alleles"
            DESC[MISS_GT]="Missing genotype despite adequate depth"
            DESC[NO_DATA]="No depth info at this position (no DP and no MIN_DP)"
            DESC[NO_GQ]="Genotype quality information missing"
            
            # Get filter counts and format as table
            # Build raw rows: tag \t description \t formatted_count \t formatted_pct
            tmp_stats_3=$(mktemp)
            bcftools query -f '%FILTER\n' "{output.bcf}" | \
                awk '{{gsub(/;/,"\n"); print}}' | \
                sort | uniq -c | \
                while read count tag; do
                    desc="${{DESC[$tag]:-}}"
                    [[ -z "$desc" ]] && desc="Unknown filter"
                    pct=$(awk -v c="$count" -v t="$TOTAL" 'BEGIN{{if(t>0) printf "%.2f", c*100/t; else print "0.00"}}')
                    printf "%s\t%s\t%'d\t%s%%\n" "$tag" "$desc" "$count" "$pct"
                done | sort -k1,1 > "$tmp_stats_3"
            
            # Dynamically size and print the Markdown table
            awk -F'\t' '
                BEGIN {{
                    max_t = length("Filter Tag")
                    max_d = length("Description")
                    max_c = length("Count")
                    max_p = length("Percentage")
                }}
                NR == FNR {{
                    if (length($1) > max_t) max_t = length($1)
                    if (length($2) > max_d) max_d = length($2)
                    if (length($3) > max_c) max_c = length($3)
                    if (length($4) > max_p) max_p = length($4)
                    next
                }}
                FNR == 1 {{
                    printf "| %-*s | %-*s | %*s | %*s |\n", max_t, "Filter Tag", max_d, "Description", max_c, "Count", max_p, "Percentage"
                    t_sep = sprintf("%*s", max_t, ""); gsub(/ /, "-", t_sep)
                    d_sep = sprintf("%*s", max_d, ""); gsub(/ /, "-", d_sep)
                    c_sep = sprintf("%*s", max_c, ""); gsub(/ /, "-", c_sep)
                    p_sep = sprintf("%*s", max_p, ""); gsub(/ /, "-", p_sep)
                    printf "|-%s-|-%s-|-%s:|-%s:|\n", t_sep, d_sep, c_sep, p_sep
                }}
                {{
                    printf "| %-*s | %-*s | %*s | %*s |\n", max_t, $1, max_d, $2, max_c, $3, max_p, $4
                }}
            ' "$tmp_stats_3" "$tmp_stats_3"
            rm -f "$tmp_stats_3"
            
            echo ""
            echo "---"
            echo ""
            echo "*Note: The sum of filter counts and percentages may exceed the total because sites can carry multiple tags simultaneously.*"
        }} > "{output.stats}"
        
        # Cleanup raw BCF if configured
        if [[ "{params.keep_raw}" == "False" ]]; then
            echo "[GAME] Removing raw BCF"
            rm -f "{input.bcf}" "{input.csi}" 2>/dev/null || true
        fi
        
        echo "[GAME] ✅ Done"
        
        '''


# ===============================================================================
#  JOINT BCF TAGGING (GLnexus output - variant-sites only, multi-sample)
# ===============================================================================
#
#  Joint genotyping tag inventory (GLnexus output):
#  ─────────────────────────────────────────────────────────────────────────
#   Tag             Type        Description
#  ─────────────────────────────────────────────────────────────────────────
#   MASKED          INFO flag   Overlaps repeat/low-complexity regions
#   REF_N           FILTER      Reference base is N
#   LOW_AQ          FILTER      INFO/AQ below threshold (GLnexus allele quality)
#   LOW_QUAL        FILTER      Site QUAL below threshold
#   NON_SNP         FILTER      Not a SNP (indel, MNP, etc.)
#   MULTI_ALT       FILTER      Multiple alternate alleles (not biallelic)
#   LOW_CALL_RATE   FILTER      Fraction of missing genotypes above threshold
#   LOW_AC          FILTER      Allele count below threshold (auto-scales)
#  ─────────────────────────────────────────────────────────────────────────
#   The joint BCF is variant-sites only (no ref-blocks), so every record
#   has a precise POS → masking annotation is always accurate.
#
#   Per-sample quality (DP, GQ, AB) is not tagged at site level - these
#   are per-genotype concerns.  Use FMT fields for per-sample filtering
#   in downstream analyses.
#
#   Pre-existing FILTER values from GLnexus (e.g., MONOALLELIC) are
#   cleared at the start of tagging.  If you need GLnexus's raw filter
#   calls, keep KEEP_RAW_JOINT: on.
# ===============================================================================

rule F02_tag_joint:
    """
    Add QC tags to the joint-genotyped BCF (GLnexus output).
    Step 1: Clear pre-existing FILTER values from GLnexus
    Step 2: INFO tags (TYPE, F_MISSING, AC, AN; MASKED if BED provided)
    Step 3: Resolve LOW_AC threshold (auto-scales with cohort size)
    Step 4: Apply site-level FILTER tags
    Step 5: Generate tagging statistics
    """
    input:
        bcf=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "jointVCFs", "{assembly}.joint.bcf"
        ),
        csi=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "jointVCFs", "{assembly}.joint.bcf.csi"
        ),
        mask=lambda w: fv_get_mask_bed(w.species, w.assembly) or []
    output:
        bcf=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "jointVCFs", "{assembly}.joint.tagged.bcf"
        ),
        csi=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "jointVCFs", "{assembly}.joint.tagged.bcf.csi"
        ),
        stats=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "jointVCFs", "{assembly}.joint.tagging_stats.md"
        )
    params:
        mask_bed=lambda w: fv_get_mask_bed(w.species, w.assembly) or "none",
        keep_raw=_KEEP_RAW_JOINT,
        min_aq=_MIN_AQ,
        min_qual=_MIN_JOINT_QUAL,
        max_missing=_MAX_MISSING,
        min_ac_cfg=_MIN_AC_CFG,
        min_ac_frac=_MIN_AC_FRAC,
    threads: cpu_func("bcftools_concat")
    resources:
        mem_mb=mem_func("bcftools_concat"),
        runtime=time_func("bcftools_concat")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "logs", "F02_tag_joint.{assembly}.log"
        )
    shell:
        r'''
        set -euo pipefail
        export LC_NUMERIC=en_US.UTF-8
        mkdir -p "$(dirname {output.bcf})" "$(dirname {log})"
        
        exec > "{log}" 2>&1

        echo "[GAME] Starting joint BCF tagging for {wildcards.assembly}"
        

        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 20)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_tag_joint_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        
        CURRENT_BCF="{input.bcf}"
        
        # ============================================
        # STEP 1: Clear pre-existing FILTER values
        # ============================================
        # GLnexus writes its own FILTER values (e.g., MONOALLELIC for sites
        # where no sample carries the alt after joint re-genotyping).  We
        # own the FILTER column after F02 - the raw joint BCF (if KEEP_RAW_JOINT
        # is on) still has GLnexus's native filter calls.
        echo "[GAME] Clearing pre-existing FILTER values"
        bcftools annotate --threads {threads} -x FILTER \
            "$CURRENT_BCF" -Ob -o cleared.bcf
        CURRENT_BCF="cleared.bcf"
        
        # ============================================
        # STEP 2: Add INFO tags (MASKED, TYPE, F_MISSING, AC, AN)
        # ============================================
        
        # --- MASKED flag ---
        # Joint BCF is variant-sites only → every record has a precise POS,
        # so masking annotation is always accurate (no ref-block ambiguity).
        if [[ "{params.mask_bed}" != "none" ]] && [[ -f "{params.mask_bed}" ]]; then
            echo "[GAME] Adding MASKED flag"
            
            ( [[ "{params.mask_bed}" =~ \.gz$ ]] && gzip -dc "{params.mask_bed}" || cat "{params.mask_bed}" ) \
              | awk 'BEGIN{{OFS="\t"}} NF>=3{{print $1,$2,$3,1}}' \
              | LC_ALL=C sort -k1,1 -k2,2n \
              | bgzip -c > mask4.bed.gz && tabix -p bed mask4.bed.gz

            echo '##INFO=<ID=MASKED,Number=0,Type=Flag,Description="Overlaps repeat/low-complexity regions">' > mask.hdr
            
            bcftools annotate \
                -a mask4.bed.gz \
                -c CHROM,FROM,TO,INFO/MASKED \
                -h mask.hdr \
                -Ob -o masked.bcf \
                "$CURRENT_BCF"
            
            CURRENT_BCF="masked.bcf"
        fi
        
        # --- Computed tags: TYPE, F_MISSING, AC, AN ---
        echo "[GAME] Adding computed tags (TYPE, F_MISSING, AC, AN)"
        bcftools +fill-tags --threads {threads} \
            "$CURRENT_BCF" \
            -Ob -o info_tagged.bcf \
            -- -t TYPE,F_MISSING,AC,AN
        
        CURRENT_BCF="info_tagged.bcf"
        
        # ============================================
        # STEP 3: Resolve LOW_AC threshold
        # ============================================
        
        MIN_AC_CFG="{params.min_ac_cfg}"
        
        if [[ "$MIN_AC_CFG" == "auto" ]]; then
            # AN varies per site (depending on missingness), so compute the
            # theoretical max AN from sample count × ploidy (assume diploid).
            # This gives a stable cohort-size-based threshold.
            N_SAMPLES=$(bcftools query -l "$CURRENT_BCF" | wc -l)
            MAX_AN=$(( N_SAMPLES * 2 ))
            MIN_AC=$(awk -v an="$MAX_AN" -v f="{params.min_ac_frac}" \
                'BEGIN{{v=int(an*f+0.999); if(v<1)v=1; print v}}')
            echo "[GAME] AUTO MIN_AC: N_SAMPLES=$N_SAMPLES, MAX_AN=$MAX_AN, fraction={params.min_ac_frac} → MIN_AC=$MIN_AC"
        else
            MIN_AC="$MIN_AC_CFG"
            echo "[GAME] Fixed MIN_AC=$MIN_AC"
        fi
        
        # ============================================
        # STEP 4: Apply site-level FILTER tags
        # ============================================
        
        echo "[GAME] Applying joint filter tags"
        echo "[GAME]   MIN_AQ={params.min_aq}, MIN_QUAL={params.min_qual}"
        echo "[GAME]   MAX_MISSING={params.max_missing}, MIN_AC=$MIN_AC"
        
        bcftools filter --threads {threads} -m + -s REF_N \
            -e 'REF="N"' -Ou "$CURRENT_BCF" | \
        bcftools filter -m + -s LOW_AQ \
            -e 'INFO/AQ<{params.min_aq}' -Ou | \
        bcftools filter -m + -s LOW_QUAL \
            -e 'QUAL>=0 && QUAL<{params.min_qual}' -Ou | \
        bcftools filter -m + -s NON_SNP \
            -e 'INFO/TYPE!="SNP"' -Ou | \
        bcftools filter -m + -s MULTI_ALT \
            -e 'N_ALT>1' -Ou | \
        bcftools filter -m + -s LOW_CALL_RATE \
            -e 'F_MISSING>{params.max_missing}' -Ou | \
        bcftools filter --threads {threads} -m + -s LOW_AC \
            -e 'INFO/AC<'$MIN_AC \
            -Ob -o "output.bcf"
        
        bcftools index -f --csi --threads {threads} "output.bcf"
        
        # Copy BCF to final location
        echo "[GAME] Copying BCF to final location..."
        cp "output.bcf" "{output.bcf}"
        cp "output.bcf.csi" "{output.csi}"
        
        # ============================================
        # STEP 5: Generate tagging statistics
        # ============================================
        {{
            echo "# Joint Genotyping Tagging Statistics"
            echo ""
            echo "**Reference:** {wildcards.assembly}  "
            echo "**Species:** {wildcards.species}  "
            echo "**Generated:** $(date)  "
            echo "**Caller:** DeepVariant + GLnexus"
            echo ""
            
            # Count samples from the BCF header
            N_SAMPLES=$(bcftools query -l "{output.bcf}" | wc -l)
            echo "**Samples:** $N_SAMPLES"
            echo ""
            echo "---"
            echo ""
            
            echo "## Quality Thresholds"
            echo ""
            
            # Write raw data to a temp file
            tmp_stats_1=$(mktemp)
            {{
                printf "Parameter\tValue\tNotes\n"
                printf "MIN_AQ\t%s\tGLnexus allele quality\n" "{params.min_aq}"
                printf "MIN_QUAL\t%s\tSite-level QUAL\n" "{params.min_qual}"
                printf "MAX_MISSING\t%s\tMax fraction of missing genotypes\n" "{params.max_missing}"
                printf "MIN_AC\t%s\tMin allele count (%s)\n" "$MIN_AC" "$MIN_AC_CFG"
            }} > "$tmp_stats_1"

            # Dynamically size and print the Markdown table (3 columns: Left, Left, Left)
            awk -F'\t' '
                BEGIN {{ max_c1=0; max_c2=0; max_c3=0 }}
                NR==FNR {{
                    if(length($1)>max_c1) max_c1=length($1)
                    if(length($2)>max_c2) max_c2=length($2)
                    if(length($3)>max_c3) max_c3=length($3)
                    next
                }}
                FNR==1 {{
                    printf "| %-*s | %-*s | %-*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                    s1=sprintf("%*s", max_c1, ""); gsub(/ /, "-", s1)
                    s2=sprintf("%*s", max_c2, ""); gsub(/ /, "-", s2)
                    s3=sprintf("%*s", max_c3, ""); gsub(/ /, "-", s3)
                    printf "|-%s-|-%s-|-%s-|\n", s1, s2, s3
                }}
                FNR>1 {{
                    printf "| %-*s | %-*s | %-*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                }}
            ' "$tmp_stats_1" "$tmp_stats_1"
            rm -f "$tmp_stats_1"
            
            echo ""
            

            echo "## Summary Counts"
            echo ""
            
            TOTAL=$(bcftools view -H "{output.bcf}" | wc -l)
            PASS_SITES=$(bcftools view -H -f 'PASS,.' "{output.bcf}" | wc -l)
            PASS_SNPS=$(bcftools view -H -f 'PASS,.' -i 'INFO/TYPE="SNP"' "{output.bcf}" | wc -l)
            PASS_PCT=$(awk -v a=$PASS_SITES -v t=$TOTAL 'BEGIN{{if(t>0) printf "%.2f", a*100/t; else print "0"}}')
            PASS_SNPS_PCT=$(awk -v a=$PASS_SNPS -v t=$TOTAL 'BEGIN{{if(t>0) printf "%.2f", a*100/t; else print "0"}}')
            
            if [[ "{params.mask_bed}" != "none" ]] && [[ -f "{params.mask_bed}" ]]; then
                MASKED=$(bcftools query -f '%INFO/MASKED\n' "{output.bcf}" | grep -c "1" || true)
                MASKED_PCT=$(awk -v m=$MASKED -v t=$TOTAL 'BEGIN{{if(t>0) printf "%.2f", m*100/t; else print "0"}}')
                PASS_MASKED=$(bcftools view -H -f 'PASS,.' -i 'INFO/MASKED=1' "{output.bcf}" | wc -l)
                PASS_MASKED_PCT=$(awk -v m=$PASS_MASKED -v t=$TOTAL 'BEGIN{{if(t>0) printf "%.2f", m*100/t; else print "0"}}')
                PASS_SNPS_MASKED=$(bcftools view -H -f 'PASS,.' -i 'INFO/TYPE="SNP" && INFO/MASKED=1' "{output.bcf}" | wc -l)
                PASS_SNPS_MASKED_PCT=$(awk -v m=$PASS_SNPS_MASKED -v t=$TOTAL 'BEGIN{{if(t>0) printf "%.2f", m*100/t; else print "0"}}')
            fi
            
            # Write raw data to a temp file
            tmp_stats_2=$(mktemp)
            {{
                printf "Category\tCount\tPercentage\n"
                printf "**Total Sites**\t%s\t100.00%%\n" "$(printf "%'d" $TOTAL)"
                if [[ "{params.mask_bed}" != "none" ]] && [[ -f "{params.mask_bed}" ]]; then
                    printf "**Total Masked**\t%s\t%s%%\n" "$(printf "%'d" $MASKED)" "$MASKED_PCT"
                else
                    printf "**Total Masked**\tn/a\tn/a (no mask BED)\n"
                fi
                printf "**PASS Sites**\t%s\t%s%%\n" "$(printf "%'d" $PASS_SITES)" "$PASS_PCT"
                if [[ "{params.mask_bed}" != "none" ]] && [[ -f "{params.mask_bed}" ]]; then
                    printf "**PASS Masked**\t%s\t%s%%\n" "$(printf "%'d" $PASS_MASKED)" "$PASS_MASKED_PCT"
                fi
                printf "**PASS SNPs**\t%s\t%s%%\n" "$(printf "%'d" $PASS_SNPS)" "$PASS_SNPS_PCT"
                if [[ "{params.mask_bed}" != "none" ]] && [[ -f "{params.mask_bed}" ]]; then
                    printf "**PASS SNPs Masked**\t%s\t%s%%\n" "$(printf "%'d" $PASS_SNPS_MASKED)" "$PASS_SNPS_MASKED_PCT"
                fi
            }} > "$tmp_stats_2"

            # Dynamically size and print the Markdown table (3 columns: Left, Right, Right)
            awk -F'\t' '
                BEGIN {{ max_c1=0; max_c2=0; max_c3=0 }}
                NR==FNR {{
                    if(length($1)>max_c1) max_c1=length($1)
                    if(length($2)>max_c2) max_c2=length($2)
                    if(length($3)>max_c3) max_c3=length($3)
                    next
                }}
                FNR==1 {{
                    printf "| %-*s | %*s | %*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                    s1=sprintf("%*s", max_c1, ""); gsub(/ /, "-", s1)
                    s2=sprintf("%*s", max_c2, ""); gsub(/ /, "-", s2)
                    s3=sprintf("%*s", max_c3, ""); gsub(/ /, "-", s3)
                    printf "|-%s-|-%s:|-%s:|\n", s1, s2, s3
                }}
                FNR>1 {{
                    printf "| %-*s | %*s | %*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                }}
            ' "$tmp_stats_2" "$tmp_stats_2"
            rm -f "$tmp_stats_2"
            
            echo ""


            echo "## Filter Tag Counts"
            echo ""
            
            declare -A DESC
            DESC[PASS]="Site passed all quality filters"
            DESC[.]="No filter applied"
            DESC[REF_N]="Reference base is N"
            DESC[LOW_AQ]="GLnexus allele quality below threshold"
            DESC[LOW_QUAL]="Site QUAL below threshold"
            DESC[NON_SNP]="Not a SNP (indel, MNP, or other)"
            DESC[MULTI_ALT]="Multiple alternate alleles (not biallelic)"
            DESC[LOW_CALL_RATE]="Too many samples with missing genotype"
            DESC[LOW_AC]="Allele count below threshold"
            



            # Build raw rows: tag \t description \t formatted_count \t formatted_pct
            tmp_stats_3=$(mktemp)
            bcftools query -f '%FILTER\n' "{output.bcf}" | \
                awk '{{gsub(/;/,"\n"); print}}' | \
                sort | uniq -c | \
                while read count tag; do
                    desc="${{DESC[$tag]:-}}"
                    [[ -z "$desc" ]] && desc="Unknown filter"
                    pct=$(awk -v c="$count" -v t="$TOTAL" 'BEGIN{{if(t>0) printf "%.2f", c*100/t; else print "0.00"}}')
                    printf "%s\t%s\t%'d\t%s%%\n" "$tag" "$desc" "$count" "$pct"
                done | sort -k1,1 > "$tmp_stats_3"
            
            # Dynamically size and print the Markdown table
            awk -F'\t' '
                BEGIN {{
                    max_t = length("Filter Tag")
                    max_d = length("Description")
                    max_c = length("Count")
                    max_p = length("Percentage")
                }}
                NR == FNR {{
                    if (length($1) > max_t) max_t = length($1)
                    if (length($2) > max_d) max_d = length($2)
                    if (length($3) > max_c) max_c = length($3)
                    if (length($4) > max_p) max_p = length($4)
                    next
                }}
                FNR == 1 {{
                    printf "| %-*s | %-*s | %*s | %*s |\n", max_t, "Filter Tag", max_d, "Description", max_c, "Count", max_p, "Percentage"
                    t_sep = sprintf("%*s", max_t, ""); gsub(/ /, "-", t_sep)
                    d_sep = sprintf("%*s", max_d, ""); gsub(/ /, "-", d_sep)
                    c_sep = sprintf("%*s", max_c, ""); gsub(/ /, "-", c_sep)
                    p_sep = sprintf("%*s", max_p, ""); gsub(/ /, "-", p_sep)
                    printf "|-%s-|-%s-|-%s:|-%s:|\n", t_sep, d_sep, c_sep, p_sep
                }}
                {{
                    printf "| %-*s | %-*s | %*s | %*s |\n", max_t, $1, max_d, $2, max_c, $3, max_p, $4
                }}
            ' "$tmp_stats_3" "$tmp_stats_3"
            rm -f "$tmp_stats_3"
            
            echo ""
            echo "---"
            echo ""
            echo "*Note: The sum of filter counts and percentages may exceed the total because sites can carry multiple tags simultaneously.*"
        }} > "{output.stats}"
        
        # Cleanup raw joint BCF if configured
        if [[ "{params.keep_raw}" == "False" ]]; then
            echo "[GAME] Removing raw joint BCF"
            rm -f "{input.bcf}" "{input.csi}" 2>/dev/null || true
        fi
        
        echo "[GAME] ✅ Done"
        
        '''


# ===============================================================================
#  RELATEDNESS CHECK - KING-robust kinship via plink2
# ===============================================================================
#
#  Detects pairs of samples with cryptic relatedness using KING-robust
#  estimation (Manichaikul et al. 2010).  Population-structure tools
#  (PCA, ADMIXTURE) and most population-genetic estimators assume
#  unrelated samples - this rule produces a "relatives" list that F04
#  uses to drop one sample from each related pair.
#
#  Pipeline:
#    1. Subset F02-tagged joint BCF → PASS biallelic SNPs (VCF.gz)
#    2. plink2 --maf + --indep-pairwise → LD-pruned marker set
#    3. plink2 --make-king-table → all-pairs KING coefficients
#    4. awk filter on KINSHIP ≥ KING_CUTOFF → relatives table
#    5. Generate stats report (marker counts, degree distribution)
#
#  KING-robust kinship interpretation (standard thresholds):
#  ─────────────────────────────────────────────────────────────────────────
#    Kinship range       Class
#  ─────────────────────────────────────────────────────────────────────────
#    > 0.354             MZ twin / duplicate
#    0.177 - 0.354       1st degree (parent-child or full sibling)
#    0.0884 - 0.177      2nd degree (half-sib, grandparent, avuncular)
#    0.0442 - 0.0884     3rd degree (1st cousin)
#    < 0.0442            Unrelated
#  ─────────────────────────────────────────────────────────────────────────
#
#  Default KING_CUTOFF = 0.0884 → flag pairs as related if 2nd-degree
#  or closer.  Tune via config.
#
#  Note: this rule only runs when joint genotyping is enabled (its
#  input is the joint BCF).  Single-sample runs skip relatedness
#  entirely - there are no pairs to assess.
# ===============================================================================

rule F03_related_check:
    """
    Detect cryptically related sample pairs using KING-robust kinship.
    Step 1: Subset F02-tagged joint BCF to PASS biallelic SNPs
    Step 2: LD prune with plink2 (--maf + --indep-pairwise)
    Step 3: Compute KING-robust all-pairs kinship table
    Step 4: Filter pairs above KING_CUTOFF → relatives list
    Step 5: Generate relatedness statistics report
    """
    input:
        bcf=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "jointVCFs", "{assembly}.joint.tagged.bcf"
        ),
        csi=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "jointVCFs", "{assembly}.joint.tagged.bcf.csi"
        ),
    output:
        kin0=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "relatedness", "{assembly}.king.kin0"
        ),
        full_kin0=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "relatedness", "{assembly}.king.full.kin0"
        ),
        relatives=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "relatedness", "{assembly}.relatives.txt"
        ),
        stats=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "relatedness", "{assembly}.relatedness_stats.md"
        )
    params:
        king_cutoff=_KING_CUTOFF,
        king_maf=_KING_MAF,
        ld_window=_KING_LD_WINDOW,
        ld_step=_KING_LD_STEP,
        ld_r2=_KING_LD_R2,
        ld_min_samples=_KING_LD_MIN_SAMPLES,
        min_snps_warn=1000,
    threads: cpu_func("plink2_king")
    resources:
        mem_mb=mem_func("plink2_king"),
        runtime=time_func("plink2_king")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "logs", "F03_related_check.{assembly}.log"
        )
    shell:
        r'''
        set -euo pipefail
        export LC_NUMERIC=en_US.UTF-8
        mkdir -p "$(dirname {output.kin0})" "$(dirname {log})"
        
        exec > "{log}" 2>&1

        echo "[GAME] Starting relatedness check for {wildcards.assembly}"

        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 20)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_king_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"

        # ============================================
        # STEP 1: Subset to PASS biallelic SNPs
        # ============================================
        # KING wants high-quality biallelic SNPs.  We use FILTER='PASS,.'
        # to match F02's stats counts (records that passed all filters
        # have FILTER='.', not 'PASS', because F02 cleared FILTER first
        # and bcftools filter only ADDS tags - it never sets PASS).
        # -m2 -M2 enforces strictly biallelic.
        echo "[GAME] Step 1: subsetting to PASS biallelic SNPs"
        bcftools view --threads {threads} \
            -f 'PASS,.' \
            -m2 -M2 \
            -v snps \
            -Oz -o hq_snps.vcf.gz \
            "{input.bcf}"
        bcftools index --threads {threads} -t hq_snps.vcf.gz

        N_INPUT_SNPS=$(bcftools view -H hq_snps.vcf.gz | wc -l)
        N_SAMPLES=$(bcftools query -l hq_snps.vcf.gz | wc -l)
        echo "[GAME] Input: $N_INPUT_SNPS PASS biallelic SNPs across $N_SAMPLES samples"

        # ============================================
        # STEP 2: MAF filter (+ optional LD prune)
        # ============================================
        # plink2 refuses to estimate r² with fewer than ~50 samples
        # (--bad-ld safety check), and the underlying issue is real:
        # at small N, r² is statistically meaningless.  We therefore
        # SKIP LD pruning when N<{params.ld_min_samples} and feed all
        # MAF-filtered SNPs directly to KING.  This is supported by
        # the KING-robust paper: kinship inference is unaffected by
        # LD per its large-sample theory, so the LD prune is a
        # convenience (smaller marker set, faster KING) rather than
        # a statistical requirement.
        # --double-id: VCF has only sample IDs (no FID); plink2 needs both,
        #              so it duplicates IID into FID.
        # --allow-extra-chr: this pipeline supports any species, so chrom
        #                    names won't always match plink2's known set.

        DO_LD_PRUNE="True"
        if [[ "$N_SAMPLES" -lt {params.ld_min_samples} ]]; then
            DO_LD_PRUNE="False"
            echo "[GAME] Step 2: MAF filter only (N=$N_SAMPLES < {params.ld_min_samples}, LD prune skipped)"
            echo "[GAME]   MAF >= {params.king_maf}"
        else
            echo "[GAME] Step 2: MAF filter + LD pruning"
            echo "[GAME]   MAF >= {params.king_maf}"
            echo "[GAME]   indep-pairwise: window={params.ld_window} step={params.ld_step} r2<{params.ld_r2}"
        fi

        if [[ "$DO_LD_PRUNE" == "True" ]]; then
            plink2 --threads {threads} \
                --vcf hq_snps.vcf.gz \
                --double-id \
                --allow-extra-chr \
                --vcf-half-call missing \
                --set-missing-var-ids '@:#:$r:$a' \
                --maf {params.king_maf} \
                --indep-pairwise {params.ld_window} {params.ld_step} {params.ld_r2} \
                --out pruning

            N_PRUNED=0
            if [[ -s pruning.prune.in ]]; then
                N_PRUNED=$(wc -l < pruning.prune.in)
            fi
        else
            # MAF-only path: plink2 --write-snplist gives us the IDs of
            # variants that pass MAF; we'll feed that list to KING via
            # --extract, mirroring the LD-prune flow.
            plink2 --threads {threads} \
                --vcf hq_snps.vcf.gz \
                --double-id \
                --allow-extra-chr \
                --vcf-half-call missing \
                --set-missing-var-ids '@:#:$r:$a' \
                --maf {params.king_maf} \
                --write-snplist \
                --out pruning

            N_PRUNED=0
            if [[ -s pruning.snplist ]]; then
                N_PRUNED=$(wc -l < pruning.snplist)
                # Reuse pruning.prune.in name so STEP 3 doesn't need to branch
                cp pruning.snplist pruning.prune.in
            fi
        fi
        echo "[GAME] Markers retained for KING: $N_PRUNED"

        WARNING_NOTE=""
        if [[ "$N_PRUNED" -eq 0 ]]; then
            echo "[GAME] WARNING: no SNPs survived marker filtering."
            if [[ "$DO_LD_PRUNE" == "True" ]]; then
                echo "[GAME] WARNING: try lowering KING_MAF or relaxing KING_LD_R2."
            else
                echo "[GAME] WARNING: try lowering KING_MAF - no MAF-passing markers exist."
            fi
            echo "[GAME] WARNING: Skipping KING - writing empty outputs."
            WARNING_NOTE="NO_MARKERS"
        elif [[ "$N_PRUNED" -lt {params.min_snps_warn} ]]; then
            echo "[GAME] WARNING: only $N_PRUNED informative SNPs available."
            echo "[GAME] WARNING: KING-robust requires several thousand markers for"
            echo "[GAME] WARNING: stable estimates - kinship values may be unreliable."
            WARNING_NOTE="LOW_SNP_COUNT"
        fi

        # ============================================
        # STEP 3: Compute KING-robust all-pairs kinship
        # ============================================
        # --make-king-table writes ALL N*(N-1)/2 pairs to .kin0.
        # Default columns: #FID1 IID1 FID2 IID2 NSNP HETHET IBS0 KINSHIP
        # If marker filtering produced no markers, we write an empty
        # kin0 with just the header and skip the plink2 KING pass -
        # F04 will see "no relatives" and pass all samples through.
        if [[ "$N_PRUNED" -eq 0 ]]; then
            echo "[GAME] Step 3: skipped (no markers); writing empty kin0"
            printf '#FID1\tIID1\tFID2\tIID2\tNSNP\tHETHET\tIBS0\tKINSHIP\n' > king_result.kin0
        else
            echo "[GAME] Step 3: computing KING-robust kinship"
            plink2 --threads {threads} \
                --vcf hq_snps.vcf.gz \
                --double-id \
                --allow-extra-chr \
                --vcf-half-call missing \
                --set-missing-var-ids '@:#:$r:$a' \
                --extract pruning.prune.in \
                --make-king-table \
                --out king_result
        fi

        cp king_result.kin0 "{output.full_kin0}"

        # ============================================
        # STEP 4: Filter to relatives (KINSHIP >= KING_CUTOFF)
        # ============================================
        echo "[GAME] Step 4: filtering to relatives (KINSHIP >= {params.king_cutoff})"

        # KINSHIP is the last column.  Using $NF makes us robust to
        # column-order changes in future plink2 versions.
        awk -F'\t' -v cutoff={params.king_cutoff} '
            NR==1 {{ print; next }}
            $NF >= cutoff
        ' king_result.kin0 > "{output.kin0}"

        # Human-readable relatives list with degree classification.
        # If user lowers KING_CUTOFF below 0.0442, "below_3rd" pairs may
        # appear (still listed as related, just classified honestly).
        {{
            printf "sample1\tsample2\tkinship\tdegree\n"
            awk -F'\t' -v cutoff={params.king_cutoff} '
                NR==1 {{ next }}
                $NF >= cutoff {{
                    kin = $NF + 0.0
                    if (kin > 0.354)        deg = "MZ_or_dup"
                    else if (kin > 0.177)   deg = "1st"
                    else if (kin > 0.0884)  deg = "2nd"
                    else if (kin > 0.0442)  deg = "3rd"
                    else                    deg = "below_3rd"
                    # IID1 = col 2, IID2 = col 4
                    printf "%s\t%s\t%.6f\t%s\n", $2, $4, kin, deg
                }}
            ' king_result.kin0
        }} > "{output.relatives}"

        N_RELATIVES=$(awk 'NR>1' "{output.relatives}" | wc -l)
        echo "[GAME] Found $N_RELATIVES related pairs at KINSHIP >= {params.king_cutoff}"

        # ============================================
        # STEP 5: Generate relatedness statistics report
        # ============================================
        {{
            echo "# Relatedness Analysis Report"
            echo ""
            echo "**Reference:** {wildcards.assembly}  "
            echo "**Species:** {wildcards.species}  "
            echo "**Generated:** $(date)  "
            echo "**Caller:** DeepVariant + GLnexus  "
            echo "**Samples:** $N_SAMPLES"
            echo ""
            echo "---"
            echo ""

            if [[ "$WARNING_NOTE" == "NO_MARKERS" ]]; then
                echo "> ⚠️  **Warning:** No SNPs survived MAF + LD pruning."
                echo "> This typically happens with very small cohorts (N<3) where LD"
                echo "> pruning is too aggressive for the available data. KING was not"
                echo "> run - no relatedness inference is possible from this dataset."
                echo ""
                echo "---"
                echo ""
            elif [[ "$WARNING_NOTE" == "LOW_SNP_COUNT" ]]; then
                echo "> ⚠️  **Warning:** Only $N_PRUNED informative SNPs after MAF+LD pruning."
                echo "> KING-robust kinship estimates require several thousand independent"
                echo "> markers for stability - values below should be interpreted with caution."
                echo ""
                echo "---"
                echo ""
            fi

            echo "## Parameters"
            echo ""
            
            # Write raw data to a temp file
            tmp_t1=$(mktemp)
            {{
                printf "Parameter\tValue\tNotes\n"
                printf "KING_CUTOFF\t%s\tKinship threshold for 'related' (default 0.0884 = 2nd degree)\n" "{params.king_cutoff}"
                printf "KING_MAF\t%s\tMinor allele frequency cutoff\n" "{params.king_maf}"
                printf "KING_LD_MIN_SAMPLES\t%s\tLD prune skipped when N<this\n" "{params.ld_min_samples}"
                if [[ "$DO_LD_PRUNE" == "True" ]]; then
                    printf "KING_LD_WINDOW\t%s\tLD-prune window (SNPs)\n" "{params.ld_window}"
                    printf "KING_LD_STEP\t%s\tLD-prune step (SNPs)\n" "{params.ld_step}"
                    printf "KING_LD_R2\t%s\tLD-prune r² threshold\n" "{params.ld_r2}"
                else
                    printf "KING_LD_WINDOW\t(not used)\tN=$N_SAMPLES below threshold - LD prune skipped\n"
                    printf "KING_LD_STEP\t(not used)\t-\n"
                    printf "KING_LD_R2\t(not used)\t-\n"
                fi
            }} > "$tmp_t1"

            # Dynamically size and print the Markdown table (3 columns: Left, Left, Left)
            awk -F'\t' '
                BEGIN {{ max_c1=0; max_c2=0; max_c3=0 }}
                NR==FNR {{
                    if(length($1)>max_c1) max_c1=length($1)
                    if(length($2)>max_c2) max_c2=length($2)
                    if(length($3)>max_c3) max_c3=length($3)
                    next
                }}
                FNR==1 {{
                    printf "| %-*s | %-*s | %-*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                    s1=sprintf("%*s", max_c1, ""); gsub(/ /, "-", s1)
                    s2=sprintf("%*s", max_c2, ""); gsub(/ /, "-", s2)
                    s3=sprintf("%*s", max_c3, ""); gsub(/ /, "-", s3)
                    printf "|-%s-|-%s-|-%s-|\n", s1, s2, s3
                }}
                FNR>1 {{
                    printf "| %-*s | %-*s | %-*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                }}
            ' "$tmp_t1" "$tmp_t1"
            rm -f "$tmp_t1"
            echo ""

            echo "## Marker Counts"
            echo ""
            
            # Write raw data to a temp file
            tmp_t2=$(mktemp)
            {{
                printf "Stage\tSNPs\n"
                printf "PASS biallelic SNPs\t%s\n" "$(printf "%'d" $N_INPUT_SNPS)"
                if [[ "$DO_LD_PRUNE" == "True" ]]; then
                    printf "After MAF + LD prune\t%s\n" "$(printf "%'d" $N_PRUNED)"
                else
                    printf "After MAF (no LD prune)\t%s\n" "$(printf "%'d" $N_PRUNED)"
                fi
            }} > "$tmp_t2"

            # Dynamically size and print the Markdown table (2 columns: Left, Right)
            awk -F'\t' '
                BEGIN {{ max_c1=0; max_c2=0 }}
                NR==FNR {{
                    if(length($1)>max_c1) max_c1=length($1)
                    if(length($2)>max_c2) max_c2=length($2)
                    next
                }}
                FNR==1 {{
                    printf "| %-*s | %*s |\n", max_c1, $1, max_c2, $2
                    s1=sprintf("%*s", max_c1, ""); gsub(/ /, "-", s1)
                    s2=sprintf("%*s", max_c2, ""); gsub(/ /, "-", s2)
                    printf "|-%s-|-%s:|\n", s1, s2
                }}
                FNR>1 {{
                    printf "| %-*s | %*s |\n", max_c1, $1, max_c2, $2
                }}
            ' "$tmp_t2" "$tmp_t2"
            rm -f "$tmp_t2"
            echo ""

            echo "## Relationship Distribution"
            echo ""
            echo "Counts derived from the full all-pairs KING table."
            echo ""
            
            # Compute counts and write to temp file
            tmp_dist=$(mktemp)
            printf "Class\tKinship range\tPairs\n" > "$tmp_dist"
            awk -F'\t' '
                NR==1 {{ next }}
                {{
                    kin = $NF + 0.0
                    total++
                    if (kin > 0.354)        mz++
                    else if (kin > 0.177)   d1++
                    else if (kin > 0.0884)  d2++
                    else if (kin > 0.0442)  d3++
                    else                    unrel++
                }}
                END {{
                    printf "MZ or duplicate\t> 0.354\t%'\''d\n", mz+0
                    printf "1st degree\t0.177 - 0.354\t%'\''d\n", d1+0
                    printf "2nd degree\t0.0884 - 0.177\t%'\''d\n", d2+0
                    printf "3rd degree\t0.0442 - 0.0884\t%'\''d\n", d3+0
                    printf "Unrelated\t< 0.0442\t%'\''d\n", unrel+0
                    printf "TOTAL\t-\t%'\''d\n", total+0
                }}
            ' "{output.full_kin0}" >> "$tmp_dist"

            # Dynamically size and print the Markdown table (3 columns: Left, Left, Right)
            awk -F'\t' '
                BEGIN {{ max_c1=0; max_c2=0; max_c3=0 }}
                NR==FNR {{
                    if(length($1)>max_c1) max_c1=length($1)
                    if(length($2)>max_c2) max_c2=length($2)
                    if(length($3)>max_c3) max_c3=length($3)
                    next
                }}
                FNR==1 {{
                    printf "| %-*s | %-*s | %*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                    s1=sprintf("%*s", max_c1, ""); gsub(/ /, "-", s1)
                    s2=sprintf("%*s", max_c2, ""); gsub(/ /, "-", s2)
                    s3=sprintf("%*s", max_c3, ""); gsub(/ /, "-", s3)
                    printf "|-%s-|-%s-|-%s:|\n", s1, s2, s3
                }}
                FNR>1 {{
                    printf "| %-*s | %-*s | %*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                }}
            ' "$tmp_dist" "$tmp_dist"
            rm -f "$tmp_dist"
            echo ""


            # Relatives list (those at or above the cutoff)
            echo "## Relatives Detected"
            echo ""
            echo "Sample pairs with KINSHIP ≥ {params.king_cutoff} (default cutoff: 2nd degree)."
            echo ""

            if [[ "$N_RELATIVES" -eq 0 ]]; then
                echo "*No related pairs detected at this threshold.*"
            else
                echo "| Sample 1 | Sample 2 | Kinship | Degree |"
                echo "|----------|----------|---------|--------|"
                awk -F'\t' 'NR>1 {{
                    printf "| %s | %s | %s | %s |\n", $1, $2, $3, $4
                }}' "{output.relatives}"
            fi
            echo ""
            echo "---"
            echo ""
            echo "*Note: KING-robust thresholds (Manichaikul et al. 2010): kinship > 0.354 = MZ/duplicate,*"
            echo "*0.177-0.354 = 1st degree, 0.0884-0.177 = 2nd degree, 0.0442-0.0884 = 3rd degree.*"
        }} > "{output.stats}"

        echo "[GAME] ✅ Done"
        '''


# ===============================================================================
#  FINAL JOINT BCF - drop one sample per related pair, then re-tag
# ===============================================================================
#
#  Reads the relatives list from F03 and resolves it into an exclusion
#  list using a greedy minimum vertex cover heuristic ranked by sample
#  quality (PASS_SNPS_PCT, parsed from F01's per-sample tagging_stats.md).
#  For each related pair, the sample with the LOWER PASS_SNPS_PCT is
#  dropped.  For trios/clusters, the algorithm iteratively drops the
#  worst-quality sample until no related pairs remain.
#
#  Then subsets the joint BCF (with --trim-alt-alleles to clean up
#  alts that become unobserved after sample removal), clears FILTER,
#  refreshes sample-dependent INFO tags (AC, AN, F_MISSING) and re-runs
#  the F02 filter chain.  The auto MIN_AC threshold is recomputed from
#  the new sample count.
#
#  Always runs (even with zero relatives) for pipeline consistency -
#  produces a freshly-tagged "clean" BCF that downstream rules can
#  depend on unconditionally.
# ===============================================================================

rule F04_final_joint:
    """
    Drop one sample per related pair (lowest PASS_SNPS_PCT) and re-tag.
    Step 1: Compute exclusion list from F03 relatives + F01 PASS_SNPS_PCT
    Step 2: Subset BCF (--trim-alt-alleles) and clear FILTER
    Step 3: Refresh INFO tags (AC, AN, F_MISSING, TYPE)
    Step 4: Resolve auto MIN_AC from new sample count
    Step 5: Apply the F02 filter chain
    Step 6: Generate final joint stats report
    """
    input:
        bcf=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "jointVCFs", "{assembly}.joint.tagged.bcf"
        ),
        csi=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "jointVCFs", "{assembly}.joint.tagged.bcf.csi"
        ),
        relatives=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "relatedness", "{assembly}.relatives.txt"
        ),
        sample_stats=_per_sample_stats_for_assembly,
    output:
        bcf=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "finalVCFs", "{assembly}.joint.clean.bcf"
        ),
        csi=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "finalVCFs", "{assembly}.joint.clean.bcf.csi"
        ),
        exclusion=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "finalVCFs", "{assembly}.exclusion.txt"
        ),
        stats=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "finalVCFs", "{assembly}.joint.clean_stats.md"
        )
    params:
        mask_bed=lambda w: fv_get_mask_bed(w.species, w.assembly) or "none",
        min_aq=_MIN_AQ,
        min_qual=_MIN_JOINT_QUAL,
        max_missing=_MAX_MISSING,
        min_ac_cfg=_MIN_AC_CFG,
        min_ac_frac=_MIN_AC_FRAC,
    threads: cpu_func("bcftools_concat")
    resources:
        mem_mb=mem_func("bcftools_concat"),
        runtime=time_func("bcftools_concat")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "logs", "F04_final_joint.{assembly}.log"
        )
    shell:
        r'''
        set -euo pipefail
        export LC_NUMERIC=en_US.UTF-8
        mkdir -p "$(dirname {output.bcf})" "$(dirname {log})"
        
        exec > "{log}" 2>&1

        echo "[GAME] Starting final joint BCF for {wildcards.assembly}"

        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 20)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_final_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"

        # ============================================
        # STEP 1: Compute exclusion list (Python)
        # ============================================
        # Reads relatives.txt → builds adjacency graph → for each connected
        # component, greedily drops the lowest-PASS_SNPS_PCT node touching
        # any edge until no edges remain.  Tiebreak: alphabetical sample ID
        # (for reproducibility across runs).
        echo "[GAME] Step 1: computing exclusion list"

        python3 - "{output.exclusion}" "{input.relatives}" {input.sample_stats} <<'PYEOF'
import sys
import os
import re
from collections import defaultdict

excl_path = sys.argv[1]
relatives_path = sys.argv[2]
stats_paths = sys.argv[3:]

# Build map: sample_id → path to its tagging_stats.md
stats_by_id = {{}}
SUFFIX = ".tagging_stats.md"
for p in stats_paths:
    bn = os.path.basename(p)
    if bn.endswith(SUFFIX):
        sid = bn[:-len(SUFFIX)]
        stats_by_id[sid] = p

# Read F03 relatives.txt: header + tab-sep rows (sample1, sample2, kinship, degree)
pairs = []
try:
    with open(relatives_path) as f:
        next(f, None)  # skip header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] and parts[1]:
                pairs.append((parts[0], parts[1]))
except FileNotFoundError:
    pairs = []

print(f"[GAME-PY] Loaded {{len(pairs)}} related pairs from {{relatives_path}}", flush=True)

# Build adjacency
adj = defaultdict(set)
for s1, s2 in pairs:
    adj[s1].add(s2)
    adj[s2].add(s1)

# Parse PASS_SNPS_PCT from F01 stats markdown.
# Expected line format (from F01 line ~523):
#   | **PASS SNPs** | 12,345,678 | 89.4523% |
PASS_LINE = re.compile(
    r"\|\s*\*\*PASS SNPs\*\*\s*\|\s*[\d,]+\s*\|\s*([\d.]+)\s*%?\s*\|"
)

def get_pass_pct(path):
    """Returns float pct, or None if not parseable."""
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            for line in f:
                m = PASS_LINE.search(line)
                if m:
                    return float(m.group(1))
    except Exception:
        return None
    return None

# Score each sample in the relatedness graph
samples_in_graph = sorted(adj.keys())
qual = {{}}
for s in samples_in_graph:
    pct = get_pass_pct(stats_by_id.get(s))
    if pct is None:
        # Sample missing from stats or unparseable - treat as worst quality
        # so it gets dropped first (safe default: prefer to keep samples
        # whose quality we KNOW is good).
        qual[s] = -1.0
        print(f"[GAME-PY] WARNING: could not parse PASS_SNPS_PCT for {{s}}; "
              f"treating as worst quality", flush=True)
    else:
        qual[s] = pct

# Greedy vertex cover: while edges remain, drop the worst-quality node
# touching any edge.  This is exact for trees and near-optimal for sparse
# graphs (which is what real relatedness graphs look like).
excluded = []
while any(adj.values()):
    active = [s for s, n in adj.items() if n]
    # Sort by (quality_pct ascending, sample_id ascending) → worst first
    active.sort(key=lambda s: (qual[s], s))
    worst = active[0]
    related_to = sorted(adj[worst])
    excluded.append((worst, qual[worst], related_to))
    # Remove this node from the graph
    for nbr in list(adj[worst]):
        adj[nbr].discard(worst)
    adj[worst] = set()

# Write exclusion list: sample_id<TAB>pass_snps_pct<TAB>related_to(comma-sep)
with open(excl_path, "w") as f:
    for sid, pct, rel_to in excluded:
        pct_str = f"{{pct:.4f}}" if pct >= 0 else "NA"
        f.write(f"{{sid}}\t{{pct_str}}\t{{','.join(rel_to)}}\n")

print(f"[GAME-PY] Excluded {{len(excluded)}} samples from "
      f"{{len(samples_in_graph)}} in the relatedness graph", flush=True)
PYEOF

        # Build the bcftools-friendly sample-IDs-only list
        awk -F'\t' '{{ print $1 }}' "{output.exclusion}" > exclusion_samples.txt
        N_EXCLUDED=$(wc -l < exclusion_samples.txt)
        echo "[GAME] Exclusion list: $N_EXCLUDED samples"

        # ============================================
        # STEP 2: Subset BCF + clear FILTER
        # ============================================
        # --trim-alt-alleles: removes alts that no remaining sample carries.
        #   Critical because dropping a sample can leave a multi-allelic site
        #   with one of the alts at AC=0; without trimming, MULTI_ALT would
        #   still fire on that site even though it's now effectively biallelic.
        # If exclusion list is empty, the -S option is skipped (a -S^ with
        # an empty file works in modern bcftools, but the conditional is
        # cheap insurance against version drift).
        echo "[GAME] Step 2: subset + clear FILTER"
        if [[ -s exclusion_samples.txt ]]; then
            bcftools view --threads {threads} \
                -S ^exclusion_samples.txt \
                --trim-alt-alleles \
                -Ou "{input.bcf}" | \
            bcftools annotate --threads {threads} -x FILTER \
                -Ob -o cleared.bcf
        else
            bcftools view --threads {threads} \
                --trim-alt-alleles \
                -Ou "{input.bcf}" | \
            bcftools annotate --threads {threads} -x FILTER \
                -Ob -o cleared.bcf
        fi
        CURRENT_BCF=cleared.bcf

        N_SAMPLES_KEPT=$(bcftools query -l "$CURRENT_BCF" | wc -l)
        echo "[GAME] Samples retained: $N_SAMPLES_KEPT"

        # ============================================
        # STEP 3: Refresh sample-dependent INFO tags
        # ============================================
        # AC, AN, F_MISSING change when samples drop.  TYPE doesn't (it's
        # determined by REF/ALT, not by genotypes), but recomputing is
        # harmless.  MASKED (added in F02) is preserved through subsetting
        # - fill-tags doesn't touch flag tags it didn't create.
        echo "[GAME] Step 3: refresh INFO tags"
        bcftools +fill-tags --threads {threads} "$CURRENT_BCF" \
            -Ob -o info_tagged.bcf -- -t TYPE,F_MISSING,AC,AN
        CURRENT_BCF=info_tagged.bcf

        # ============================================
        # STEP 4: Resolve auto MIN_AC threshold
        # ============================================
        # Recompute from new sample count - a 1% AC fraction means a
        # different absolute count when the cohort is smaller.
        MIN_AC_CFG="{params.min_ac_cfg}"
        if [[ "$MIN_AC_CFG" == "auto" ]]; then
            MAX_AN=$(( N_SAMPLES_KEPT * 2 ))
            MIN_AC=$(awk -v an="$MAX_AN" -v f="{params.min_ac_frac}" \
                'BEGIN{{v=int(an*f+0.999); if(v<1)v=1; print v}}')
            echo "[GAME] AUTO MIN_AC: N_SAMPLES=$N_SAMPLES_KEPT, MAX_AN=$MAX_AN, fraction={params.min_ac_frac} -> MIN_AC=$MIN_AC"
        else
            MIN_AC="$MIN_AC_CFG"
            echo "[GAME] Fixed MIN_AC=$MIN_AC"
        fi

        # ============================================
        # STEP 5: Apply F02 filter chain
        # ============================================
        echo "[GAME] Step 5: applying filter chain"
        echo "[GAME]   MIN_AQ={params.min_aq}, MIN_QUAL={params.min_qual}"
        echo "[GAME]   MAX_MISSING={params.max_missing}, MIN_AC=$MIN_AC"
        bcftools filter --threads {threads} -m + -s REF_N \
            -e 'REF="N"' -Ou "$CURRENT_BCF" | \
        bcftools filter -m + -s LOW_AQ \
            -e 'INFO/AQ<{params.min_aq}' -Ou | \
        bcftools filter -m + -s LOW_QUAL \
            -e 'QUAL>=0 && QUAL<{params.min_qual}' -Ou | \
        bcftools filter -m + -s NON_SNP \
            -e 'INFO/TYPE!="SNP"' -Ou | \
        bcftools filter -m + -s MULTI_ALT \
            -e 'N_ALT>1' -Ou | \
        bcftools filter -m + -s LOW_CALL_RATE \
            -e 'F_MISSING>{params.max_missing}' -Ou | \
        bcftools filter --threads {threads} -m + -s LOW_AC \
            -e 'INFO/AC<'$MIN_AC \
            -Ob -o "{output.bcf}"

        bcftools index --threads {threads} -c "{output.bcf}"

        # ============================================
        # STEP 6: Stats report
        # ============================================
        echo "[GAME] Step 6: generating stats report"

        TOTAL=$(bcftools view -H "{output.bcf}" | wc -l)
        PASS_SITES=$(bcftools view -H -f 'PASS,.' "{output.bcf}" | wc -l)
        PASS_SNPS=$(bcftools view -H -f 'PASS,.' -i 'INFO/TYPE="SNP"' "{output.bcf}" | wc -l)
        PASS_PCT=$(awk -v a=$PASS_SITES -v t=$TOTAL 'BEGIN{{if(t>0) printf "%.2f", a*100/t; else print "0"}}')
        PASS_SNPS_PCT=$(awk -v a=$PASS_SNPS -v t=$TOTAL 'BEGIN{{if(t>0) printf "%.2f", a*100/t; else print "0"}}')

        # Per-FILTER tag counts (one tag per output line, even when ; combined)
        tmp_tags=$(mktemp)
        bcftools view -H "{output.bcf}" | awk -F'\t' '{{
            n = split($7, a, ";")
            for (i=1;i<=n;i++) print a[i]
        }}' | sort | uniq -c | awk '{{printf "%s\t%d\n", $2, $1}}' > "$tmp_tags"

        {{
            echo "# Final Joint BCF - Tagging Stats Report"
            echo ""
            echo "**Reference:** {wildcards.assembly}  "
            echo "**Species:** {wildcards.species}  "
            echo "**Generated:** $(date)  "
            echo "**Caller:** DeepVariant + GLnexus (post-relatedness)  "
            echo "**Samples retained:** $N_SAMPLES_KEPT"
            echo ""
            echo "---"
            echo ""


            # ---- Excluded samples section ----
            echo "## Samples Excluded by Relatedness"
            echo ""
            if [[ "$N_EXCLUDED" -eq 0 ]]; then
                echo "*No samples excluded - no related pairs detected at the configured KING_CUTOFF.*"
            else
                echo "Samples dropped to break related pairs (greedy: lowest PASS_SNPS_PCT first)."
                echo ""
                
                tmp_t1=$(mktemp)
                printf "Sample\tPASS_SNPS_PCT\tRelated to\n" > "$tmp_t1"
                awk -F'\t' '{{
                    pct = ($2 == "NA") ? "NA" : $2 "%"
                    printf "%s\t%s\t%s\n", $1, pct, $3
                }}' "{output.exclusion}" >> "$tmp_t1"

                # 3 columns: Left, Right, Left
                awk -F'\t' '
                    BEGIN {{ max_c1=0; max_c2=0; max_c3=0 }}
                    NR==FNR {{
                        if(length($1)>max_c1) max_c1=length($1)
                        if(length($2)>max_c2) max_c2=length($2)
                        if(length($3)>max_c3) max_c3=length($3)
                        next
                    }}
                    FNR==1 {{
                        printf "| %-*s | %*s | %-*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                        s1=sprintf("%*s", max_c1, ""); gsub(/ /, "-", s1)
                        s2=sprintf("%*s", max_c2, ""); gsub(/ /, "-", s2)
                        s3=sprintf("%*s", max_c3, ""); gsub(/ /, "-", s3)
                        printf "|-%s-|-%s:|-%s-|\n", s1, s2, s3
                    }}
                    FNR>1 {{
                        printf "| %-*s | %*s | %-*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                    }}
                ' "$tmp_t1" "$tmp_t1"
                rm -f "$tmp_t1"
            fi
            echo ""

            # ---- Quality thresholds ----
            echo "## Quality Thresholds"
            echo ""
            tmp_t2=$(mktemp)
            {{
                printf "Filter\tThreshold\tNotes\n"
                printf "LOW_AQ\t< %s\tGLnexus allele quality\n" "{params.min_aq}"
                printf "LOW_QUAL\t< %s\tSite QUAL field\n" "{params.min_qual}"
                printf "LOW_CALL_RATE\tF_MISSING > %s\tFraction of missing genotypes\n" "{params.max_missing}"
                printf "LOW_AC\t< %s\tAuto-scaled from N_SAMPLES=%s, fraction=%s\n" "$MIN_AC" "$N_SAMPLES_KEPT" "{params.min_ac_frac}"
            }} > "$tmp_t2"

            # 3 columns: Left, Left, Left
            awk -F'\t' '
                BEGIN {{ max_c1=0; max_c2=0; max_c3=0 }}
                NR==FNR {{
                    if(length($1)>max_c1) max_c1=length($1)
                    if(length($2)>max_c2) max_c2=length($2)
                    if(length($3)>max_c3) max_c3=length($3)
                    next
                }}
                FNR==1 {{
                    printf "| %-*s | %-*s | %-*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                    s1=sprintf("%*s", max_c1, ""); gsub(/ /, "-", s1)
                    s2=sprintf("%*s", max_c2, ""); gsub(/ /, "-", s2)
                    s3=sprintf("%*s", max_c3, ""); gsub(/ /, "-", s3)
                    printf "|-%s-|-%s-|-%s-|\n", s1, s2, s3
                }}
                FNR>1 {{
                    printf "| %-*s | %-*s | %-*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                }}
            ' "$tmp_t2" "$tmp_t2"
            rm -f "$tmp_t2"
            echo ""

            # ---- Summary counts ----
            echo "## Summary Counts"
            echo ""
            tmp_t3=$(mktemp)
            {{
                printf "Metric\tCount\t%%\n"
                printf "Total Sites\t%s\t100.00%%\n" "$(printf "%'d" $TOTAL)"
                printf "**PASS Sites**\t%s\t%s%%\n" "$(printf "%'d" $PASS_SITES)" "$PASS_PCT"
                printf "**PASS SNPs**\t%s\t%s%%\n" "$(printf "%'d" $PASS_SNPS)" "$PASS_SNPS_PCT"
            }} > "$tmp_t3"

            # 3 columns: Left, Right, Right
            awk -F'\t' '
                BEGIN {{ max_c1=0; max_c2=0; max_c3=0 }}
                NR==FNR {{
                    if(length($1)>max_c1) max_c1=length($1)
                    if(length($2)>max_c2) max_c2=length($2)
                    if(length($3)>max_c3) max_c3=length($3)
                    next
                }}
                FNR==1 {{
                    printf "| %-*s | %*s | %*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                    s1=sprintf("%*s", max_c1, ""); gsub(/ /, "-", s1)
                    s2=sprintf("%*s", max_c2, ""); gsub(/ /, "-", s2)
                    s3=sprintf("%*s", max_c3, ""); gsub(/ /, "-", s3)
                    printf "|-%s-|-%s:|-%s:|\n", s1, s2, s3
                }}
                FNR>1 {{
                    printf "| %-*s | %*s | %*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                }}
            ' "$tmp_t3" "$tmp_t3"
            rm -f "$tmp_t3"
            echo ""

            # ---- Filter tag counts ----
            echo "## Filter Tag Counts"
            echo ""
            declare -A DESC
            DESC[PASS]="Site passed all quality filters"
            DESC[.]="No filter applied"
            DESC[REF_N]="Reference base is N"
            DESC[LOW_AQ]="GLnexus allele quality below threshold"
            DESC[LOW_QUAL]="Site QUAL below threshold"
            DESC[NON_SNP]="Not a SNP (indel, MNP, or other)"
            DESC[MULTI_ALT]="Multiple alternate alleles (not biallelic)"
            DESC[LOW_CALL_RATE]="Too many samples with missing genotype"
            DESC[LOW_AC]="Allele count below threshold"

            tmp_t4=$(mktemp)
            printf "Tag\tCount\tDescription\n" > "$tmp_t4"
            while IFS=$'\t' read -r tag cnt; do
                desc="${{DESC[$tag]:-}}"
                if [[ -z "$desc" ]]; then desc="(no description)"; fi
                printf "%s\t%s\t%s\n" "$tag" "$(printf "%'d" $cnt)" "$desc" >> "$tmp_t4"
            done < "$tmp_tags"

            # 3 columns: Left, Right, Left
            awk -F'\t' '
                BEGIN {{ max_c1=0; max_c2=0; max_c3=0 }}
                NR==FNR {{
                    if(length($1)>max_c1) max_c1=length($1)
                    if(length($2)>max_c2) max_c2=length($2)
                    if(length($3)>max_c3) max_c3=length($3)
                    next
                }}
                FNR==1 {{
                    printf "| %-*s | %*s | %-*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                    s1=sprintf("%*s", max_c1, ""); gsub(/ /, "-", s1)
                    s2=sprintf("%*s", max_c2, ""); gsub(/ /, "-", s2)
                    s3=sprintf("%*s", max_c3, ""); gsub(/ /, "-", s3)
                    printf "|-%s-|-%s:|-%s-|\n", s1, s2, s3
                }}
                FNR>1 {{
                    printf "| %-*s | %*s | %-*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                }}
            ' "$tmp_t4" "$tmp_t4"
            rm -f "$tmp_t4"
            echo ""

            echo "---"
            echo ""
            echo "*Generated by F04_final_joint.  This BCF is the relatedness-cleaned*"
            echo "*counterpart to {wildcards.assembly}.joint.tagged.bcf.*"
        }} > "{output.stats}"

        rm -f "$tmp_tags"
        echo "[GAME] ✅ Done"
        '''


# ===============================================================================
#  LD PRUNING TAG - flag LD-redundant sites in the final joint BCF
# ===============================================================================
#
#  Adds an INFO/LD_REDUNDANT flag to sites that are in linkage
#  disequilibrium with another retained site, as identified by
#  plink2 --indep-pairwise.  Following the philosophy of MASKED:
#  the flag is informational (INFO, not FILTER) because LD-redundancy
#  is an analysis-context concern, not a quality concern.  PCA and
#  ADMIXTURE want the flag respected; coalescent or demographic
#  inference often does not.
#
#  Reuses the same MAF and LD pruning parameters as F03 (KING_MAF,
#  KING_LD_WINDOW, KING_LD_STEP, KING_LD_R2) for now.  These can be
#  split into dedicated LDPRUNE_* config keys later if needed.
#
#  Output BCF is the canonical "final" cohort variant set: relatedness-
#  cleaned (F04) + LD-pruning annotation (F05).
# ===============================================================================

rule F05_ld_prune:
    """
    Add INFO/LD_REDUNDANT flag to LD-redundant sites in the clean BCF.
    Step 1: Subset clean BCF to PASS biallelic SNPs
    Step 2: plink2 --maf + --indep-pairwise → prune.out (LD-redundant sites)
    Step 3: Convert prune.out to bgzipped BED
    Step 4: bcftools annotate → final BCF with LD_REDUNDANT flag
    Step 5: Generate stats report
    """
    input:
        bcf=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "finalVCFs", "{assembly}.joint.clean.bcf"
        ),
        csi=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "finalVCFs", "{assembly}.joint.clean.bcf.csi"
        ),
    output:
        bcf=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "finalVCFs", "{assembly}.joint.final.bcf"
        ),
        csi=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "finalVCFs", "{assembly}.joint.final.bcf.csi"
        ),
        stats=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "finalVCFs", "{assembly}.ld_prune_stats.md"
        )
    params:
        ld_maf=_KING_MAF,
        ld_window=_KING_LD_WINDOW,
        ld_step=_KING_LD_STEP,
        ld_r2=_KING_LD_R2,
        ld_min_samples=_KING_LD_MIN_SAMPLES,
    threads: cpu_func("plink2_king")
    resources:
        mem_mb=mem_func("plink2_king"),
        runtime=time_func("plink2_king")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "logs", "F05_ld_prune.{assembly}.log"
        )
    shell:
        r'''
        set -euo pipefail
        export LC_NUMERIC=en_US.UTF-8
        mkdir -p "$(dirname {output.bcf})" "$(dirname {log})"
        
        exec > "{log}" 2>&1

        echo "[GAME] Starting LD-prune annotation for {wildcards.assembly}"

        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 20)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_ldprune_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"

        # ============================================
        # STEP 0: Sample-count check
        # ============================================
        # plink2's r² estimates are unreliable below ~50 samples; LD
        # pruning at small N would produce a misleading LD_REDUNDANT
        # flag that downstream analyses might trust.  When N is below
        # the configured threshold, we skip plink2 entirely, copy the
        # F04 clean BCF unchanged as the "final" BCF, and document this
        # in the stats report.  This applies to N AFTER F04's
        # relatedness exclusions, since the threshold is about the
        # actual sample count being analysed.
        N_SAMPLES=$(bcftools query -l "{input.bcf}" | wc -l)
        echo "[GAME] N samples in F04 clean BCF: $N_SAMPLES"

        if [[ "$N_SAMPLES" -lt {params.ld_min_samples} ]]; then
            echo "[GAME] N=$N_SAMPLES < KING_LD_MIN_SAMPLES={params.ld_min_samples}"
            echo "[GAME] Skipping LD pruning - copying F04 clean BCF as final."
            
            # Just copy the BCF and re-index.  We use bcftools view as a
            # no-op pass instead of cp because it ensures the output has
            # a fresh header timestamp and is regenerated atomically by
            # Snakemake (avoids hard-link / inode confusion on shared FS).
            bcftools view --threads {threads} -Ob -o "{output.bcf}" "{input.bcf}"
            bcftools index --threads {threads} -c "{output.bcf}"

            # Compute counts for the stats report
            TOTAL=$(bcftools view -H "{output.bcf}" | wc -l)
            PASS_SITES=$(bcftools view -H -f 'PASS,.' "{output.bcf}" | wc -l)
            N_PASS_SNPS=$(bcftools view -H -f 'PASS,.' -m2 -M2 -v snps "{output.bcf}" | wc -l)

            {{
                echo "# LD-Prune Annotation Report"
                echo ""
                echo "**Reference:** {wildcards.assembly}  "
                echo "**Species:** {wildcards.species}  "
                echo "**Generated:** $(date)  "
                echo "**Caller:** DeepVariant + GLnexus (post-relatedness)  "
                echo "**Samples:** $N_SAMPLES"
                echo ""
                echo "---"
                echo ""
                echo "> ⚠️  **LD pruning skipped.**"
                echo ">"
                echo "> Cohort size (N=$N_SAMPLES) is below KING_LD_MIN_SAMPLES={params.ld_min_samples}."
                echo "> plink2's r² estimates are unreliable at small sample sizes, so the"
                echo "> LD_REDUNDANT flag would be statistically meaningless at this N."
                echo ">"
                echo "> The final BCF is therefore identical to the F04 clean BCF, with"
                echo "> NO INFO/LD_REDUNDANT annotation.  If your downstream analysis"
                echo "> requires LD pruning, do it in-analysis using the relevant tool's"
                echo "> built-in handling (e.g., PCA tools that account for small N)."
                echo ""
                echo "## Marker Counts"
                echo ""
                tmp_t1=$(mktemp)
                {{
                    printf "Stage\tCount\n"
                    printf "Total sites in BCF\t%s\n" "$(printf "%'d" $TOTAL)"
                    printf "PASS sites\t%s\n" "$(printf "%'d" $PASS_SITES)"
                    printf "PASS biallelic SNPs\t%s\n" "$(printf "%'d" $N_PASS_SNPS)"
                }} > "$tmp_t1"

                # 2 columns: Left, Right
                awk -F'\t' '
                    BEGIN {{ max_c1=0; max_c2=0 }}
                    NR==FNR {{
                        if(length($1)>max_c1) max_c1=length($1)
                        if(length($2)>max_c2) max_c2=length($2)
                        next
                    }}
                    FNR==1 {{
                        printf "| %-*s | %*s |\n", max_c1, $1, max_c2, $2
                        s1=sprintf("%*s", max_c1, ""); gsub(/ /, "-", s1)
                        s2=sprintf("%*s", max_c2, ""); gsub(/ /, "-", s2)
                        printf "|-%s-|-%s:|\n", s1, s2
                    }}
                    FNR>1 {{
                        printf "| %-*s | %*s |\n", max_c1, $1, max_c2, $2
                    }}
                ' "$tmp_t1" "$tmp_t1"
                rm -f "$tmp_t1"
                echo ""
                echo "---"
                echo ""
                echo "*Generated by F05_ld_prune (skipped path).  This is the canonical"
                echo "final cohort BCF: relatedness-cleaned (F04), no LD annotation*"
                echo "*because cohort size was too small.*"
            }} > "{output.stats}"

            echo "[GAME] ✅ Done (LD prune skipped)"
            exit 0
        fi

        echo "[GAME] N=$N_SAMPLES >= {params.ld_min_samples} - proceeding with LD prune"

        # ============================================
        # STEP 1: Subset to PASS biallelic SNPs
        # ============================================
        # Same gate as F03 - LD pruning is meaningful only on the
        # high-quality biallelic SNP subset.  Note this BCF is from F04
        # (post-relatedness), so AC/AF/F_MISSING already reflect the
        # cleaned cohort.
        echo "[GAME] Step 1: subsetting to PASS biallelic SNPs"
        bcftools view --threads {threads} \
            -f 'PASS,.' \
            -m2 -M2 \
            -v snps \
            -Oz -o hq_snps.vcf.gz \
            "{input.bcf}"
        bcftools index --threads {threads} -t hq_snps.vcf.gz

        N_PASS_SNPS=$(bcftools view -H hq_snps.vcf.gz | wc -l)
        echo "[GAME] PASS biallelic SNPs in input: $N_PASS_SNPS"

        # ============================================
        # STEP 2: plink2 LD prune
        # ============================================
        echo "[GAME] Step 2: plink2 --maf + --indep-pairwise"
        echo "[GAME]   MAF >= {params.ld_maf}"
        echo "[GAME]   indep-pairwise: window={params.ld_window} step={params.ld_step} r2<{params.ld_r2}"

        plink2 --threads {threads} \
            --vcf hq_snps.vcf.gz \
            --double-id \
            --allow-extra-chr \
            --vcf-half-call missing \
            --set-missing-var-ids '@:#:$r:$a' \
            --maf {params.ld_maf} \
            --indep-pairwise {params.ld_window} {params.ld_step} {params.ld_r2} \
            --out pruning

        N_PRUNE_OUT=0
        if [[ -s pruning.prune.out ]]; then
            N_PRUNE_OUT=$(wc -l < pruning.prune.out)
        fi
        N_PRUNE_IN=0
        if [[ -s pruning.prune.in ]]; then
            N_PRUNE_IN=$(wc -l < pruning.prune.in)
        fi
        echo "[GAME] Independent (kept) SNPs: $N_PRUNE_IN"
        echo "[GAME] LD-redundant SNPs:       $N_PRUNE_OUT"

        # ============================================
        # STEP 3: Convert prune.out variant IDs to bgzipped BED
        # ============================================
        # plink2 variant IDs (with --set-missing-var-ids '@:#:$r:$a') are
        # CHROM:POS:REF:ALT.  Awk on ':' to extract CHROM (field 1) and
        # POS (field 2).  BED is half-open 0-based, so FROM=POS-1, TO=POS.
        # Some chromosome names contain ':' (rare, but possible) - guard
        # by reconstructing CHROM as everything except the last 3 colon-
        # separated fields.
        echo "[GAME] Step 3: building LD-redundant BED"

        if [[ "$N_PRUNE_OUT" -eq 0 ]]; then
            echo "[GAME] No LD-redundant sites - writing empty BED."
            # Empty bgzipped BED with tabix index - annotate is a no-op
            : | bgzip -c > prune_out.bed.gz
            tabix -f -p bed prune_out.bed.gz || true
        else
            awk 'BEGIN{{OFS="\t"}} {{
                # ID = CHROM:POS:REF:ALT, where CHROM may contain colons
                # (e.g. "chr1_GL000123:1234"). Strip last 3 colon fields.
                n = split($1, a, ":")
                if (n < 4) next
                chrom = a[1]
                for (i=2; i<=n-3; i++) chrom = chrom ":" a[i]
                pos = a[n-2] + 0
                if (pos < 1) next
                print chrom, pos-1, pos
            }}' pruning.prune.out \
              | LC_ALL=C sort -k1,1 -k2,2n \
              | bgzip -c > prune_out.bed.gz
            tabix -f -p bed prune_out.bed.gz
        fi

        # ============================================
        # STEP 4: Annotate input BCF with LD_REDUNDANT flag
        # ============================================
        echo "[GAME] Step 4: annotating final BCF"
        echo '##INFO=<ID=LD_REDUNDANT,Number=0,Type=Flag,Description="Site in LD with another retained site (plink2 --indep-pairwise '"{params.ld_window}"' '"{params.ld_step}"' '"{params.ld_r2}"', MAF>='"{params.ld_maf}"')">' > ld.hdr

        bcftools annotate --threads {threads} \
            -a prune_out.bed.gz \
            -c CHROM,FROM,TO,INFO/LD_REDUNDANT \
            -h ld.hdr \
            -Ob -o "{output.bcf}" \
            "{input.bcf}"

        bcftools index --threads {threads} -c "{output.bcf}"

        # ============================================
        # STEP 5: Stats report
        # ============================================
        echo "[GAME] Step 5: generating stats report"

        TOTAL=$(bcftools view -H "{output.bcf}" | wc -l)
        PASS_SITES=$(bcftools view -H -f 'PASS,.' "{output.bcf}" | wc -l)
        N_TAGGED=$(bcftools view -H -i 'INFO/LD_REDUNDANT=1' "{output.bcf}" | wc -l)
        N_INDEP=$(bcftools view -H -f 'PASS,.' -m2 -M2 -v snps -e 'INFO/LD_REDUNDANT=1' "{output.bcf}" | wc -l)

        N_SAMPLES=$(bcftools query -l "{output.bcf}" | wc -l)

        {{
            echo "# LD-Prune Annotation Report"
            echo ""
            echo "**Reference:** {wildcards.assembly}  "
            echo "**Species:** {wildcards.species}  "
            echo "**Generated:** $(date)  "
            echo "**Caller:** DeepVariant + GLnexus (post-relatedness, LD-tagged)  "
            echo "**Samples:** $N_SAMPLES"
            echo ""
            echo "---"
            echo ""

            echo "## Parameters"
            echo ""
            tmp_p1=$(mktemp)
            {{
                printf "Parameter\tValue\tSource\n"
                printf "MAF cutoff\t%s\tKING_MAF\n" "{params.ld_maf}"
                printf "LD window (SNPs)\t%s\tKING_LD_WINDOW\n" "{params.ld_window}"
                printf "LD step (SNPs)\t%s\tKING_LD_STEP\n" "{params.ld_step}"
                printf "LD r² threshold\t%s\tKING_LD_R2\n" "{params.ld_r2}"
            }} > "$tmp_p1"

            # 3 columns: Left, Left, Left
            awk -F'\t' '
                BEGIN {{ max_c1=0; max_c2=0; max_c3=0 }}
                NR==FNR {{
                    if(length($1)>max_c1) max_c1=length($1)
                    if(length($2)>max_c2) max_c2=length($2)
                    if(length($3)>max_c3) max_c3=length($3)
                    next
                }}
                FNR==1 {{
                    printf "| %-*s | %-*s | %-*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                    s1=sprintf("%*s", max_c1, ""); gsub(/ /, "-", s1)
                    s2=sprintf("%*s", max_c2, ""); gsub(/ /, "-", s2)
                    s3=sprintf("%*s", max_c3, ""); gsub(/ /, "-", s3)
                    printf "|-%s-|-%s-|-%s-|\n", s1, s2, s3
                }}
                FNR>1 {{
                    printf "| %-*s | %-*s | %-*s |\n", max_c1, $1, max_c2, $2, max_c3, $3
                }}
            ' "$tmp_p1" "$tmp_p1"
            rm -f "$tmp_p1"
            echo ""

            echo "## Marker Counts"
            echo ""
            tmp_p2=$(mktemp)
            {{
                printf "Stage\tCount\n"
                printf "Total sites in BCF\t%s\n" "$(printf "%'d" $TOTAL)"
                printf "PASS sites\t%s\n" "$(printf "%'d" $PASS_SITES)"
                printf "PASS biallelic SNPs (input to plink2)\t%s\n" "$(printf "%'d" $N_PASS_SNPS)"
                printf "LD-redundant (LD_REDUNDANT=1)\t%s\n" "$(printf "%'d" $N_TAGGED)"
                printf "**Independent set (PASS biallelic SNPs, not LD-redundant)**\t%s\n" "$(printf "%'d" $N_INDEP)"
            }} > "$tmp_p2"

            # 2 columns: Left, Right
            awk -F'\t' '
                BEGIN {{ max_c1=0; max_c2=0 }}
                NR==FNR {{
                    if(length($1)>max_c1) max_c1=length($1)
                    if(length($2)>max_c2) max_c2=length($2)
                    next
                }}
                FNR==1 {{
                    printf "| %-*s | %*s |\n", max_c1, $1, max_c2, $2
                    s1=sprintf("%*s", max_c1, ""); gsub(/ /, "-", s1)
                    s2=sprintf("%*s", max_c2, ""); gsub(/ /, "-", s2)
                    printf "|-%s-|-%s:|\n", s1, s2
                }}
                FNR>1 {{
                    printf "| %-*s | %*s |\n", max_c1, $1, max_c2, $2
                }}
            ' "$tmp_p2" "$tmp_p2"
            rm -f "$tmp_p2"
            echo ""

            echo "---"
            echo ""
            echo "*Generated by F05_ld_prune.  This is the canonical final cohort"
            echo "BCF: relatedness-cleaned (F04) + LD-redundancy annotation (F05).*"
        }} > "{output.stats}"

        echo "[GAME] ✅ Done"
        '''
