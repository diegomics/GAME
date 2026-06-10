# ===============================================================================
# GAME - Variant Calling Rules
# by Diego De Panis, 2026
# note: AI tools may have been used to improve, clean and/or comment this version of the code
# ===============================================================================

# -------------------------------------------------------------------------------
#  OUTPUTH PATH TEMPLATES (literal patterns only) ---> add deepv hifi/illu? -> add gatk joint?
# -------------------------------------------------------------------------------

IDX_DIR_TMPL        = os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "idx")
REF_FA_TMPL         = os.path.join(IDX_DIR_TMPL, "{assembly}.fa")
REF_FAI_TMPL        = os.path.join(IDX_DIR_TMPL, "{assembly}.fa.fai")
REF_DICT_TMPL       = os.path.join(IDX_DIR_TMPL, "{assembly}.dict")
INTERVALS_DIR_TMPL  = os.path.join(IDX_DIR_TMPL, "intervals")

SAMPLE_BAMS_TMPL    = os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs")
SAMPLE_VCFS_TMPL    = os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}", "VCFs")

INTERVAL_LIST_TMPL  = os.path.join(INTERVALS_DIR_TMPL, "interval_{shard}.list")
GVCF_PART_TMPL      = os.path.join(SAMPLE_VCFS_TMPL, "tmp", "{sample_id}.bp.{shard}.g.vcf.gz")
GVCF_PART_TBI_TMPL  = GVCF_PART_TMPL + ".tbi"
GENO_PART_TMPL      = os.path.join(SAMPLE_VCFS_TMPL, "tmp", "{sample_id}.bp.{shard}.genotyped.vcf.gz")
GENO_PART_TBI_TMPL  = GENO_PART_TMPL + ".tbi"

# FINAL_* templates are defined after _BP_MODE is set (see "Resolution-aware final paths" below)

# Caller-specific intermediates (routed to FINAL_* by E20/E21)
GATK_GVCF_TMPL     = os.path.join(SAMPLE_VCFS_TMPL, "{sample_id}.gatk.g.vcf.gz")
GATK_GVCF_TBI_TMPL = GATK_GVCF_TMPL + ".tbi"
GATK_BCF_TMPL      = os.path.join(SAMPLE_VCFS_TMPL, "{sample_id}.gatk.raw.bcf")
GATK_BCF_CSI_TMPL  = GATK_BCF_TMPL + ".csi"

DV_GVCF_TMPL       = os.path.join(SAMPLE_VCFS_TMPL, "{sample_id}.dv.g.vcf.gz")
DV_GVCF_TBI_TMPL   = DV_GVCF_TMPL + ".tbi"
DV_BCF_TMPL        = os.path.join(SAMPLE_VCFS_TMPL, "{sample_id}.dv.raw.bcf")
DV_BCF_CSI_TMPL    = DV_BCF_TMPL + ".csi"

# Joint genotyping output (per species/assembly, not per sample)
JOINT_VCFS_DIR_TMPL = os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "jointVCFs")
JOINT_BCF_TMPL      = os.path.join(JOINT_VCFS_DIR_TMPL, "{assembly}.joint.bcf")
JOINT_BCF_CSI_TMPL  = JOINT_BCF_TMPL + ".csi"


# -------------------------------------------------------------------------------
#  HELPER FUNCTIONS (prefixed to avoid collisions) 
# -------------------------------------------------------------------------------

def vc_intervals_dir(sp, asm):
    return os.path.join(config["OUT_FOLDER"], "GAME_results", sp, asm, "idx", "intervals")

def vc_priority_list():
    raw = str(config.get("DATA_PRIORITY", "hifi>illumina>ont"))
    return [x.strip().lower() for x in re.split(r"[>\s,]+", raw) if x.strip()]

def vc_get_calling_bam(sp, asm, sid):
    """Get the appropriate BAM for variant calling based on merge decision"""
    
    # Check if a merged BAM exists (from merge_tech_bams)
    merged_bam = os.path.join(
        config["OUT_FOLDER"], "GAME_results", sp, asm,
        "samples", sid, "BAMs", f"{sid}.merged.bam"
    )
    
    # For multi-tech samples, always use the merged BAM (symlink or actual merge)
    try:
        sample_data = samples_config["sp_name"][sp]["asm_id"][asm]["sample_id"][sid]
        tech_count = len([rt for rt in sample_data.get("read_type", {}) 
                         if rt not in ["None", None]])
        
        if tech_count > 1:
            # Multi-tech sample - use the merged BAM (which might be a symlink)
            return merged_bam
        
        # Single tech sample - use the tech-specific BAM
        for rt_key in sample_data.get("read_type", {}):
            if rt_key not in ["None", None]:
                rt_norm = normalize_read_type(rt_key)
                return os.path.join(
                    config["OUT_FOLDER"], "GAME_results", sp, asm,
                    "samples", sid, "BAMs", f"{sid}.{rt_norm}.merged.bam"
                )
    except Exception as e:
        print(f"Warning: Error getting BAM for {sp}/{asm}/{sid}: {e}")
    
    raise ValueError(f"Could not determine BAM file for {sp}/{asm}/{sid}")


def vc_interval_path(sp, asm, shard):
    return os.path.join(vc_intervals_dir(sp, asm), f"interval_{shard}.list")

def vc_discover_shards(sp, asm):
    # force interval checkpoint to run
    _ = checkpoints.E00_make_intervals.get(species=sp, assembly=asm)
    files = sorted(glob.glob(os.path.join(vc_intervals_dir(sp, asm), "interval_*.list")))
    return [os.path.splitext(os.path.basename(x))[0].split("_", 1)[1] for x in files]

# config knobs
_SAMPLE_PLOIDY   = int(config.get("SAMPLE_PLOIDY", 2))
_KEEP_BP_GVCF    = bool(config.get("KEEP_gVCF", True))

# VCF resolution: "block" (GVCF with ref-blocks, default) or "basepair" (every site)
##_RESOLUTION = str(config.get("RESOLUTION", "block")).strip().lower()
##if _RESOLUTION not in ("block", "basepair"):
##    raise ValueError(
##        f"RESOLUTION must be 'block' or 'basepair', got: '{_RESOLUTION}'"
##    )
##_BP_MODE = (_RESOLUTION == "basepair")

# Resolution-aware final paths
# basepair -> {sid}.bp.g.vcf.gz / {sid}.bp.raw.bcf  (name reflects per-site content)
# block    -> {sid}.g.vcf.gz    / {sid}.raw.bcf      (no prefix needed)
##_BP_TAG = ".bp" if _BP_MODE else ""

FINAL_GVCF_TMPL     = os.path.join(SAMPLE_VCFS_TMPL, "{sample_id}" + _BP_TAG + ".g.vcf.gz")
FINAL_GVCF_TBI_TMPL = FINAL_GVCF_TMPL + ".tbi"
FINAL_BCF_TMPL      = os.path.join(SAMPLE_VCFS_TMPL, "{sample_id}" + _BP_TAG + ".raw.bcf")
FINAL_BCF_CSI_TMPL  = FINAL_BCF_TMPL + ".csi"

# Variant caller: "deep" (DeepVariant, default) or "gatk" (GATK HaplotypeCaller)
# Multi-tech samples always fall back to GATK regardless of this setting
_CALLER = str(config.get("CALLER", "deep")).strip().lower()
if _CALLER not in ("deep", "gatk"):
    raise ValueError(
        f"CALLER must be 'deep' or 'gatk', got: '{_CALLER}'"
    )
if _CALLER == "deep":
    print("[GAME] Variant caller: DeepVariant (multi-tech samples will fall back to GATK)")
else:
    print("[GAME] Variant caller: GATK HaplotypeCaller")

# Joint genotyping via GLnexus (requires CALLER: deep + MERGE_TECH: off)
# YAML parses on/off as True/False — handle both booleans and strings
_joint_raw = config.get("JOINT_GENO", False)
_JOINT_GENO = _joint_raw is True or str(_joint_raw).strip().lower() in ("true", "yes", "on", "1")

if _JOINT_GENO:
    if _CALLER != "deep":
        raise ValueError(
            "JOINT_GENO requires CALLER: deep (joint genotyping uses GLnexus, "
            "which is designed for DeepVariant gVCFs)"
        )
    # MERGE_TECH not parsed yet here — we'll validate after _MERGE_TECH is set (below)

    # Auto-force KEEP_gVCF on (gVCFs are needed as input to GLnexus)
    _keep_raw = config.get("KEEP_gVCF", True)
    _keep_on = _keep_raw is True or str(_keep_raw).strip().lower() in ("true", "yes", "on", "1")
    if not _keep_on:
        print("[GAME] JOINT_GENO=on → setting KEEP_gVCF=on (required for joint genotyping)")
        config["KEEP_gVCF"] = True
    _KEEP_BP_GVCF = True
    print("[GAME] Joint genotyping: enabled (GLnexus)")
else:
    print("[GAME] Joint genotyping: disabled")


def _get_sample_read_category(species, assembly, sample_id):
    """
    Determine if sample uses short or long reads.
    Returns 'short' for illumina/10x, 'long' for hifi/ont.
    If mixed, returns 'long' (more conservative settings).
    """
    try:
        sample_data = samples_config["sp_name"][species]["asm_id"][assembly]["sample_id"][sample_id]
        read_types = sample_data.get("read_type", {})
        
        has_short = False
        has_long = False
        
        for rt_key in read_types:
            if rt_key in ["None", None]:
                continue
            if _is_long(rt_key):
                has_long = True
            elif _is_short(rt_key):
                has_short = True
        
        # If mixed or long reads present, use long read settings (more conservative)
        if has_long:
            return "long"
        elif has_short:
            return "short"
        else:
            return "short"  # default
    except (KeyError, TypeError, AttributeError):
        return "short"


def _vc_pairhmm_impl(read_category):
    """
    Determine PairHMM implementation based on read category and USE_AVX config.
    
    USE_AVX options:
      auto: FASTEST_AVAILABLE for short reads, LOGLESS_CACHING for long reads
      on:   FASTEST_AVAILABLE for all (may crash on some CPUs or with long reads)
      off:  LOGLESS_CACHING for all (safest, slowest)
    """
    mode = str(config.get("USE_AVX", "auto")).lower()
    
    if mode in {"off", "false", "no"}:
        # Disabled: always use safe implementation
        return "LOGLESS_CACHING"
    
    if mode in {"on", "true", "yes"}:
        # Forced on: use AVX for everything (user takes responsibility)
        return "FASTEST_AVAILABLE"
    
    # Auto mode (default): AVX for short reads only
    if read_category == "long":
        return "LOGLESS_CACHING"
    else:
        return "FASTEST_AVAILABLE"


def vc_hc_params(wildcards):
    """
    Return HaplotypeCaller parameters based on read type.
    
    Short reads: optimized for high coverage, PCR artifacts present
    Long reads: optimized for longer reads, no PCR, safer pairhmm
    """
    read_category = _get_sample_read_category(
        wildcards.species, wildcards.assembly, wildcards.sample_id
    )
    
    pairhmm = _vc_pairhmm_impl(read_category)
    
    if read_category == "long":
        return {
            "read_category": "long",
            "pairhmm": pairhmm,
            "pcr_indel_model": "NONE",            # No PCR artifacts in long reads
            "max_haplotypes": 64,                 # Fewer haplotypes needed
            "max_reads_start": 40,                # Each read covers more
            "max_asm_region": 200,                # Smaller assembly regions
        }
    else:  # short reads
        return {
            "read_category": "short",
            "pairhmm": pairhmm,
            "pcr_indel_model": "CONSERVATIVE",    # Handle PCR artifacts
            "max_haplotypes": 128,                # More haplotypes for high coverage
            "max_reads_start": 50,                # Standard for short reads
            "max_asm_region": 300,                # Larger assembly regions
        }



# -------------------------------------------------------------------------------
#  DeepVariant helpers
# -------------------------------------------------------------------------------

# YAML parses unquoted on/off as True/False booleans — normalise back to strings
_merge_raw = config.get("MERGE_TECH", "auto")
if _merge_raw is False:
    _MERGE_TECH = "off"
elif _merge_raw is True:
    _MERGE_TECH = "on"
else:
    _MERGE_TECH = str(_merge_raw).strip().lower()

# Deferred JOINT_GENO validation: MERGE_TECH must be "off"
if _JOINT_GENO and _MERGE_TECH != "off":
    raise ValueError(
        f"JOINT_GENO requires MERGE_TECH: off (got '{_MERGE_TECH}'). "
        "Joint genotyping uses GLnexus with DeepVariant gVCFs, which requires all "
        "samples to be called with DeepVariant. Set MERGE_TECH: off to ensure no "
        "multi-tech merging triggers a GATK fallback."
    )

def _sample_techs(sp, asm, sid):
    """Return list of normalized read types that have actual reads for a sample."""
    try:
        sample_data = samples_config["sp_name"][sp]["asm_id"][asm]["sample_id"][sid]
        techs = []
        for rt_key in sample_data.get("read_type", {}):
            if rt_key not in ["None", None]:
                rt_norm = normalize_read_type(rt_key)
                _, _, reads = _get_sample_node(sp, asm, sid, rt_key)
                if reads:
                    techs.append(rt_norm)
        return techs
    except (KeyError, TypeError, AttributeError):
        return []


def _pick_priority_tech(sp, asm, sid):
    """Pick the highest-priority tech for a sample based on DATA_PRIORITY."""
    techs = _sample_techs(sp, asm, sid)
    if not techs:
        return None
    if len(techs) == 1:
        return techs[0]
    # Walk priority list, return first match
    for prio in vc_priority_list():
        if prio in techs:
            return prio
    # No match in priority list → return first available
    return techs[0]


def _caller_for_sample(sp, asm, sid):
    """Decide caller at DAG time per sample.

    - 1 tech                        → honour CALLER setting
    - >1 tech + MERGE_TECH: off     → honour CALLER (runs on priority tech only)
    - >1 tech + MERGE_TECH: auto    → always GATK (merge outcome unknown at DAG time)
    - >1 tech + MERGE_TECH: on      → always GATK (techs will be merged)
    """
    if _CALLER != "deep":
        return "gatk"
    techs = _sample_techs(sp, asm, sid)
    if len(techs) <= 1:
        return "deep"
    # Multi-tech sample
    if _MERGE_TECH == "off":
        prio = _pick_priority_tech(sp, asm, sid)
        print(f"[GAME] Sample {sid} has {len(techs)} techs, MERGE_TECH=off "
              f"→ DeepVariant on priority tech '{prio}'")
        return "deep"
    else:
        print(f"[GAME] ⚠️  Sample {sid} has {len(techs)} techs, MERGE_TECH={_MERGE_TECH} "
              f"→ falling back to GATK")
        return "gatk"


def _dv_model_type(sp, asm, sid):
    """Map the chosen tech to a DeepVariant --model_type."""
    tech = _pick_priority_tech(sp, asm, sid)
    if tech is None:
        return "WGS"
    if tech in ("hifi",):
        return "PACBIO"
    elif tech in ("ont",):
        return "ONT_R104"
    else:                    # illumina, 10x, etc.
        return "WGS"


def _dv_get_bam(sp, asm, sid):
    """Get the BAM path for DeepVariant (uses priority tech for multi-tech + off)."""
    tech = _pick_priority_tech(sp, asm, sid)
    if tech is None:
        raise ValueError(f"No reads found for {sp}/{asm}/{sid}")
    return os.path.join(
        config["OUT_FOLDER"], "GAME_results", sp, asm,
        "samples", sid, "BAMs", f"{sid}.{tech}.merged.bam"
    )


# -------------------------------------------------------------------------------
#  Caller routing (select caller-specific intermediate → final output)
# -------------------------------------------------------------------------------

def _routed_caller_file(w, tmpl_deep, tmpl_gatk):
    """Route to the correct caller-specific intermediate file."""
    caller = _caller_for_sample(w.species, w.assembly, w.sample_id)
    tmpl = tmpl_deep if caller == "deep" else tmpl_gatk
    return tmpl.format(
        species=w.species, assembly=w.assembly, sample_id=w.sample_id
    )


# -------------------------------------------------------------------------------
#  REF READINESS (depend on map_reads rule)
# ===============================================================================

rule E00_ref_faidx_for_calling:
    input:
        fa=rules.D00_ref_fa_canonical.output.fa
    output:
        fai=REF_FAI_TMPL
    threads: cpu_func("samtools_index")
    resources:
        mem_mb=mem_func("samtools_index"),
        runtime=time_func("samtools_index")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "logs", "E00_ref_faidx.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {log})
        
        exec > {log} 2>&1
        
        echo "[GAME] Creating FASTA index for {input.fa}"
        

        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 25)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_faidx_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        
        echo "[GAME] Copying assembly to temp directory..."
        cp "{input.fa}" local.fa
        
        echo "[GAME] Running samtools faidx..."
        samtools faidx local.fa
        
        echo "[GAME] Copying index back to output..."
        cp local.fa.fai "{output.fai}"
        
        echo "[GAME] Done"
        
        '''

rule E00_ref_dict_for_calling:
    input:
        fa=rules.D00_ref_fa_canonical.output.fa
    output:
        dict=REF_DICT_TMPL
    threads: cpu_func("samtools_index")
    resources:
        mem_mb=mem_func("samtools_index"), 
        runtime=time_func("samtools_index")
    container: CONTAINERS["gatk"]
    log:
        os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "logs", "E00_ref_dict.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {log})
        
        exec > {log} 2>&1
        
        echo "[GAME] Creating sequence dictionary for {input.fa}"
        

        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 25)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_dict_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        
        echo "[GAME] Copying assembly to temp directory..."
        cp "{input.fa}" local.fa
        
        echo "[GAME] Running gatk CreateSequenceDictionary..."
        gatk CreateSequenceDictionary -R local.fa -O local.dict
        
        echo "[GAME] Copying dictionary back to output..."
        cp local.dict "{output.dict}"
        
        echo "[GAME] Done"
        
        '''


# -------------------------------------------------------------------------------
#  PER ASM INTERVALS (no contig splitting)
# ===============================================================================

checkpoint E00_make_intervals:
    input:
        fa=rules.D00_ref_fa_canonical.output.fa,
        fai=rules.E00_ref_faidx_for_calling.output.fai,
        dict=rules.E00_ref_dict_for_calling.output.dict
    output:
        idir=directory(INTERVALS_DIR_TMPL)
    params:
        target_bp=int(50_000_000)
    threads: cpu_func("split_intervals")
    resources:
        mem_mb=mem_func("split_intervals"), runtime=time_func("split_intervals")
    container: CONTAINERS["game_base"]
    run:
        from pathlib import Path
        idir = Path(output.idir); idir.mkdir(parents=True, exist_ok=True)
        contigs = []
        with open(input.fai) as fh:
            for ln in fh:
                toks = ln.rstrip("\n").split("\t")
                if len(toks) >= 2: contigs.append((toks[0], int(toks[1])))

        shard, current, acc, tgt = 1, [], 0, params.target_bp
        def _flush():
            nonlocal shard, current, acc
            if not current: return
            with open(idir / f"interval_{shard:03d}.list", "w") as f:
                for n, L in current: f.write(f"{n}:1-{L}\n")
            shard += 1; current.clear(); acc = 0

        for name, L in contigs:
            if L > tgt and not current:
                with open(idir / f"interval_{shard:03d}.list", "w") as f:
                    f.write(f"{name}:1-{L}\n")
                shard += 1
            elif acc + L <= tgt:
                current.append((name, L)); acc += L
            else:
                _flush(); current.append((name, L)); acc = L
        _flush()

        with open(idir / "MANIFEST.json", "w") as mf:
            json.dump({"intervals": sorted(str(p) for p in idir.glob("interval_*.list"))}, mf, indent=2)


# -------------------------------------------------------------------------------
#  HaplotypeCaller (BP_RESOLUTION), scattered
# ===============================================================================

rule E01_haplotypecaller_scattered:
    input:
        bam=lambda w: vc_get_calling_bam(w.species, w.assembly, w.sample_id),
        # Add merge decision as dependency if it exists
        decision=lambda w: os.path.join(
            config["OUT_FOLDER"], "GAME_results", w.species, w.assembly,
            "samples", w.sample_id, "BAMs", "merge_decision.json"
        ) if len([rt for rt in samples_config["sp_name"][w.species]["asm_id"][w.assembly]["sample_id"][w.sample_id].get("read_type", {}) 
                 if rt not in ["None", None]]) > 1 else [],
        ref=rules.D00_ref_fa_canonical.output.fa,
        fai=rules.E00_ref_faidx_for_calling.output.fai,
        dict=rules.E00_ref_dict_for_calling.output.dict,
        interval=lambda w: vc_interval_path(w.species, w.assembly, w.shard)
    output:
        gvcf=temp(GVCF_PART_TMPL),
        tbi=temp(GVCF_PART_TBI_TMPL)
    wildcard_constraints:
        shard=r"\d{3}"
    params:
        ploidy=_SAMPLE_PLOIDY,
        hc_params=vc_hc_params,
        bp_mode=_BP_MODE,
    threads: cpu_func("gatk_haplotypecaller")
    resources:
        mem_mb=mem_func("gatk_haplotypecaller"), runtime=time_func("gatk_haplotypecaller")
    container: CONTAINERS["gatk"]
    log:
        os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}", "logs",
                     "E01_hc.{sample_id}.{shard}.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.gvcf})" "$(dirname {log})"
        
        exec > {log} 2>&1
        
        # Ensure BAM has a valid index (bai or csi); repair if empty/missing
        if [ -f "{input.bam}.bai" ] && [ ! -s "{input.bam}.bai" ]; then rm -f "{input.bam}.bai"; fi
        if [ ! -s "{input.bam}.bai" ] && [ ! -s "{input.bam}.csi" ]; then
            echo "[GAME] Re-indexing BAM before HC"
            samtools index -@ {threads} "{input.bam}" || true
        fi    


        # TEMP DIRECTORY SETUP
        # Write output to local fast storage, then copy to final location
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 25)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_hc_{wildcards.species}_{wildcards.assembly}_{wildcards.sample_id}_{wildcards.shard}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        echo "[GAME] Working directory: $TEMP_DIR"


        # JAVA/GATK MEMORY SETTINGS
        # -------------------------------------------------------------------
        # samjdk.buffer_size (bytes) based on available RAM
        #  < 8 GB  → 256 KiB
        #  8–16 GB → 512 KiB
        # 16–64 GB →   1 MiB
        # ≥64 GB   →   2 MiB
        MEM_MB={resources.mem_mb}
        XMX=$(( MEM_MB - 2048 )); [ "$XMX" -lt 1024 ] && XMX=1024
        if   [ "$MEM_MB" -ge 65536 ]; then BUF=2097152
        elif [ "$MEM_MB" -ge 16384 ]; then BUF=1048576
        elif [ "$MEM_MB" -ge  8192 ]; then BUF=524288
        else                               BUF=262144
        fi


        # READ-TYPE AWARE PARAMETERS
        # -------------------------------------------------------------------
        READ_CAT="{params.hc_params[read_category]}"
        PAIRHMM="{params.hc_params[pairhmm]}"
        PCR_MODEL="{params.hc_params[pcr_indel_model]}"
        MAX_HAP="{params.hc_params[max_haplotypes]}"
        MAX_READS="{params.hc_params[max_reads_start]}"
        MAX_ASM="{params.hc_params[max_asm_region]}"

        # Use native PairHMM threads only when a native implementation is selected
        NTH={threads}
        if [ "$PAIRHMM" = "LOGLESS_CACHING" ]; then NTH=1; fi

        echo "[GAME] HaplotypeCaller for $READ_CAT reads (resolution={params.bp_mode})"
        echo "[GAME] Memory: XMX=${{XMX}}m, buffer=${{BUF}}, threads={threads}"
        echo "[GAME] Parameters: pairhmm=$PAIRHMM, pcr=$PCR_MODEL, haplotypes=$MAX_HAP, reads_start=$MAX_READS, asm_region=$MAX_ASM"

        # RESOLUTION-DEPENDENT FLAGS
        # -------------------------------------------------------------------
        # basepair: -ERC BP_RESOLUTION + --output-mode EMIT_ALL_CONFIDENT_SITES
        #           (emits a record for every single base pair)
        # block:    -ERC GVCF  (collapses non-variant regions into ref-blocks)
        #           --output-mode EMIT_ALL_CONFIDENT_SITES is incompatible with
        #           GVCF mode and must NOT be used
        if [ "{params.bp_mode}" = "True" ]; then
            ERC_FLAG="-ERC BP_RESOLUTION --output-mode EMIT_ALL_CONFIDENT_SITES"
            echo "[GAME] Resolution: basepair (BP_RESOLUTION)"
        else
            ERC_FLAG="-ERC GVCF"
            echo "[GAME] Resolution: block (GVCF)"
        fi

        # Run HaplotypeCaller - output to local temp
        gatk --java-options "-Xmx${{XMX}}m -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=1 -Dsamjdk.buffer_size=${{BUF}} \
          -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true" \
          HaplotypeCaller \
          -R "{input.ref}" \
          -I "{input.bam}" \
          -L "{input.interval}" \
          -O "output.g.vcf.gz" \
          $ERC_FLAG \
          --sample-ploidy {params.ploidy} \
          --native-pair-hmm-threads $NTH \
          --pair-hmm-implementation $PAIRHMM \
          --pcr-indel-model $PCR_MODEL \
          --max-num-haplotypes-in-population $MAX_HAP \
          --max-reads-per-alignment-start $MAX_READS \
          --max-assembly-region-size $MAX_ASM \
          --interval-padding 0 \
          --tmp-dir "$TEMP_DIR"

        gatk IndexFeatureFile -I "output.g.vcf.gz"
        
        # Copy results to final location
        echo "[GAME] Copying output to final location..."
        cp "output.g.vcf.gz" "{output.gvcf}"
        cp "output.g.vcf.gz.tbi" "{output.tbi}"
        
        echo "[GAME] Done"
        
        '''


# -------------------------------------------------------------------------------
#  GATHER BP gVCF per sample
# ===============================================================================

rule E02_gather_bp_gvcf:
    input:
        gvcfs=lambda w: [
            GVCF_PART_TMPL.format(
                species=w.species, assembly=w.assembly,
                sample_id=w.sample_id, shard=s
            )
            for s in vc_discover_shards(w.species, w.assembly)
        ],
        tbis=lambda w: [
            GVCF_PART_TBI_TMPL.format(
                species=w.species, assembly=w.assembly,
                sample_id=w.sample_id, shard=s
            )
            for s in vc_discover_shards(w.species, w.assembly)
        ]
    output:
        gvcf=temp(GATK_GVCF_TMPL),
        tbi=temp(GATK_GVCF_TBI_TMPL)
    threads: cpu_func("gatk_gathervcfs")
    resources:
        mem_mb=mem_func("gatk_gathervcfs"), runtime=time_func("gatk_gathervcfs")
    container: CONTAINERS["gatk"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "E02_gather_gvcf.{sample_id}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.gvcf})" "$(dirname {log})"
        
        exec > {log} 2>&1


        # TEMP DIRECTORY SETUP
        # Write output to local fast storage, then copy to final location
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 50)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_gather_gvcf_{wildcards.species}_{wildcards.assembly}_{wildcards.sample_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        echo "[GAME] Working directory: $TEMP_DIR"

        MEM_MB={resources.mem_mb}
        XMX=$(( MEM_MB - 2048 )); [ "$XMX" -lt 1024 ] && XMX=1024
        if   [ "$MEM_MB" -ge 65536 ]; then BUF=2097152
        elif [ "$MEM_MB" -ge 16384 ]; then BUF=1048576
        elif [ "$MEM_MB" -ge  8192 ]; then BUF=524288
        else                               BUF=262144
        fi

        # Count the actual number of gVCF files (not including .tbi)
        set -- {input.gvcfs}
        N=$#
        echo "[GAME] gather_bp_gvcf: N=$N gVCFs to gather"

        if [ "$N" -eq 1 ]; then
            # Fast path: single shard → just copy + reindex
            SRC="$1"
            echo "[GAME] Single gVCF, using fast path (copy)"
            cp "$SRC" "output.g.vcf.gz"
            gatk IndexFeatureFile -I "output.g.vcf.gz"
        else
            echo "[GAME] Multiple gVCFs, using GatherVcfs"
            INLIST=""
            for f in {input.gvcfs}; do 
                INLIST="$INLIST -I $f"
            done

            gatk --java-options "-Xmx${{XMX}}m -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=${{BUF}} \
                                 -Dsamjdk.compression_level=2 \
                                 -Dsamjdk.use_async_io_read_samtools=true \
                                 -Dsamjdk.use_async_io_write_samtools=true" \
                 GatherVcfs $INLIST -O "output.g.vcf.gz" --TMP_DIR "$TEMP_DIR"

            gatk IndexFeatureFile -I "output.g.vcf.gz"
        fi
        
        # Copy results to final location
        echo "[GAME] Copying output to final location..."
        cp "output.g.vcf.gz" "{output.gvcf}"
        cp "output.g.vcf.gz.tbi" "{output.tbi}"
        
        echo "[GAME] Done"
        
        '''

# -------------------------------------------------------------------------------
#  SCATTERED SINGLE-SAMPLE GENOTYPING + FINAL BCF
# ===============================================================================

rule E03_genotype_scattered:
    input:
        ref=rules.D00_ref_fa_canonical.output.fa,
        fai=rules.E00_ref_faidx_for_calling.output.fai,
        dict=rules.E00_ref_dict_for_calling.output.dict,
        gvcf=GVCF_PART_TMPL,
        tbi=GVCF_PART_TBI_TMPL,
        interval=lambda w: vc_interval_path(w.species, w.assembly, w.shard)
    output:
        vcf=temp(GENO_PART_TMPL),
        tbi=temp(GENO_PART_TBI_TMPL)
    wildcard_constraints:
        shard=r"\d{3}"
    threads: cpu_func("gatk_genotypegvcfs")
    resources:
        mem_mb=mem_func("gatk_genotypegvcfs"), runtime=time_func("gatk_genotypegvcfs")
    container: CONTAINERS["gatk"]
    params:
        bp_mode=_BP_MODE,
    log:
        os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}", "logs",
                     "E03_genotype.{sample_id}.{shard}.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.vcf})" "$(dirname {log})"
        
        exec > {log} 2>&1


        # TEMP DIRECTORY SETUP
        # Write output to local fast storage, then copy to final location
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 20)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_genotype_{wildcards.species}_{wildcards.assembly}_{wildcards.sample_id}_{wildcards.shard}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        echo "[GAME] Working directory: $TEMP_DIR"

        MEM_MB={resources.mem_mb}
        XMX=$(( MEM_MB - 2048 )); [ "$XMX" -lt 1024 ] && XMX=1024
        if   [ "$MEM_MB" -ge 65536 ]; then BUF=2097152
        elif [ "$MEM_MB" -ge 16384 ]; then BUF=1048576
        elif [ "$MEM_MB" -ge  8192 ]; then BUF=524288
        else                               BUF=262144
        fi

        # RESOLUTION-DEPENDENT FLAGS
        # -------------------------------------------------------------------
        # basepair: emit every site individually (threshold 0 + non-variant sites)
        # block:    emit variant sites + reference blocks (default QUAL threshold)
        if [ "{params.bp_mode}" = "True" ]; then
            EXTRA_FLAGS="--include-non-variant-sites --standard-min-confidence-threshold-for-calling 0"
            echo "[GAME] Resolution: basepair (all sites emitted individually)"
        else
            EXTRA_FLAGS="--include-non-variant-sites"
            echo "[GAME] Resolution: block (variant sites + reference blocks)"
        fi

        # Run GenotypeGVCFs - output to local temp
        gatk --java-options "-Xmx${{XMX}}m -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=${{BUF}} \
          -Dsamjdk.compression_level=1 -Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true" \
          GenotypeGVCFs -R "{input.ref}" -V "{input.gvcf}" -L "{input.interval}" \
          --interval-padding 0 \
          -O "output.vcf.gz" \
          $EXTRA_FLAGS \
          --tmp-dir "$TEMP_DIR"

        gatk IndexFeatureFile -I "output.vcf.gz"
        
        # Copy results to final location
        echo "[GAME] Copying output to final location..."
        cp "output.vcf.gz" "{output.vcf}"
        cp "output.vcf.gz.tbi" "{output.tbi}"
        
        echo "[GAME] Done"
        
        '''

rule E04_gather_genotyped_to_bcf:
    input:
        vcfs=lambda w: [ GENO_PART_TMPL.format(species=w.species, assembly=w.assembly, sample_id=w.sample_id, shard=s)
                    for s in vc_discover_shards(w.species, w.assembly) ],
        tbis=lambda w: [ GENO_PART_TBI_TMPL.format(species=w.species, assembly=w.assembly, sample_id=w.sample_id, shard=s)
                    for s in vc_discover_shards(w.species, w.assembly) ]
    output:
        bcf=temp(GATK_BCF_TMPL),
        csi=temp(GATK_BCF_CSI_TMPL)
    threads: cpu_func("bcftools_concat")
    resources:
        mem_mb=mem_func("bcftools_concat"), runtime=time_func("bcftools_concat")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}", "logs",
                     "E04_gather_bcf.{sample_id}.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.bcf})" "$(dirname {log})"
        
        exec > {log} 2>&1


        # TEMP DIRECTORY SETUP
        # Write output to local fast storage, then copy to final location
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 50)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_gather_bcf_{wildcards.species}_{wildcards.assembly}_{wildcards.sample_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        echo "[GAME] Working directory: $TEMP_DIR"

        # How many genotyped VCFs did we get?
        set -- {input.vcfs}
        N=$#
        echo "[GAME] gather_genotyped_to_bcf: inputs=$N, threads={threads}"
        
        if [ "$N" -eq 1 ]; then
            SRC="$1"
            echo "[GAME] Single input, converting to BCF"
            # single input → just convert to BCF
            export HTS_COMPRESSION_LEVEL=1
            bcftools view -O b -o "output.bcf" "$SRC" --threads {threads}
        else
            echo "[GAME] Multiple inputs, indexing and concatenating"
            # Index all input files first (in case they're not indexed)
            for vcf in {input.vcfs}; do
                if [ ! -f "${{vcf}}.tbi" ] && [ ! -f "${{vcf}}.csi" ]; then
                    echo "[GAME] Indexing $vcf"
                    bcftools index -t "$vcf"
                fi
            done
            
            # multiple inputs → concatenate to BCF
            export HTS_COMPRESSION_LEVEL=1
            bcftools concat -a -D -O b -o "output.bcf" {input.vcfs} --threads {threads}
        fi
        
        echo "[GAME] Creating CSI index for BCF"
        bcftools index -f --csi --threads {threads} "output.bcf"
        
        # Copy results to final location
        echo "[GAME] Copying output to final location..."
        cp "output.bcf" "{output.bcf}"
        cp "output.bcf.csi" "{output.csi}"
        
        echo "[GAME] Done"

        # Clean up empty tmp directory if all temp files are gone
        TMP_DIR="$(dirname {input.vcfs[0]})"
        if [ -d "$TMP_DIR" ] && [ -z "$(ls -A "$TMP_DIR" 2>/dev/null)" ]; then
            rmdir "$TMP_DIR" 2>/dev/null || true
            echo "[GAME] Removed empty tmp directory"
        fi
        
        '''

# ---------

# -------------------------------------------------------------------------------
#  Path template for D06's output (already defined in D, but re-state here
#  for clarity since E references it). Keep in sync with D06_infer_sex_ploidy.
# -------------------------------------------------------------------------------
SEX_INFER_TSV_TMPL = os.path.join(
    config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
    "samples", "{sample_id}", "BAMs", "qc", "{sample_id}.sex_inference.tsv"
)


# -------------------------------------------------------------------------------
#  Reader: extract haploid contigs from D06's sex_inference.tsv
# -------------------------------------------------------------------------------
def _read_inferred_haploid_contigs(tsv_path):
    """
    Parse sex_inference.tsv and return the list of contigs whose call is
    "haploid". Ambiguous and diploid rows are excluded — they are NOT passed
    to --haploid_contigs (safe default per design).

    Returns [] if the file is empty (header only) or unreadable.
    """
    contigs = []
    try:
        with open(tsv_path) as f:
            header = f.readline().strip().split("\t")
            if "contig" not in header or "call" not in header:
                return []
            i_contig = header.index("contig")
            i_call = header.index("call")
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(i_contig, i_call):
                    continue
                if parts[i_call].strip() == "haploid":
                    contigs.append(parts[i_contig].strip())
    except (FileNotFoundError, IOError):
        # Should not happen at rule-execution time (Snakemake guarantees the
        # input file exists), but degrade safely if it does.
        return []
    return contigs


# -------------------------------------------------------------------------------
#  Input function: conditionally depend on sex_inference.tsv
# -------------------------------------------------------------------------------
def _dv_sex_inference_input(w):
    """
    Return the path to sex_inference.tsv as a list — empty list if this
    sample does not need inference (declared sample_sex, missing sex_chr,
    or "none" case).

    Returning [] tells Snakemake "no extra dependency" — D06 will not be
    invoked for this sample.
    """
    contigs, source = resolve_haploid_contigs(
        samples_config, w.species, w.assembly, w.sample_id
    )
    if source == "inferred":
        return [SEX_INFER_TSV_TMPL.format(
            species=w.species, assembly=w.assembly, sample_id=w.sample_id
        )]
    # "declared", "missing_sex_chr", or "none" — no TSV needed.
    return []


# -------------------------------------------------------------------------------
#  Params function: build the --haploid_contigs flag string for make_examples
# -------------------------------------------------------------------------------
def _dv_haploid_contigs_flag(w):
    """
    Build the --haploid_contigs flag string for DeepVariant's make_examples.

    Returns either:
      - "--haploid_contigs=contigA,contigB"  (when contigs to flag exist)
      - ""                                    (no flag, default diploid)

    Selection logic:
      - source == "declared": use the contigs returned by the resolver.
      - source == "inferred": read sex_inference.tsv (from D06) and use
                              only rows where call == "haploid". Ambiguous
                              and diploid rows are dropped (safe default).
      - source == "missing_sex_chr" or "none": no flag.
    """
    contigs, source = resolve_haploid_contigs(
        samples_config, w.species, w.assembly, w.sample_id
    )

    if source == "inferred":
        tsv = SEX_INFER_TSV_TMPL.format(
            species=w.species, assembly=w.assembly, sample_id=w.sample_id
        )
        contigs = _read_inferred_haploid_contigs(tsv)

    if not contigs:
        return ""
    return f"--haploid_contigs={','.join(contigs)}"


# ===============================================================================
#  DeepVariant PATH (single rule replaces the GATK scatter-gather chain)
# ===============================================================================

rule E10_deepvariant:
    """Run DeepVariant as three separate binaries with per-step env.

    Why split from run_deepvariant: make_examples (N parallel processes) needs
    thread caps to avoid RLIMIT_NPROC / cgroup pids.max blowup on clusters.
    call_variants (1 process) needs unrestricted XLA/TF threading for multi-core
    matmul. run_deepvariant sets env once, so we split to control each step.
    """
    input:
        ref=rules.D00_ref_fa_canonical.output.fa,
        fai=rules.E00_ref_faidx_for_calling.output.fai,
        bam=lambda w: _dv_get_bam(w.species, w.assembly, w.sample_id),
        bai=lambda w: _dv_get_bam(w.species, w.assembly, w.sample_id) + ".bai",
        sex_inference=_dv_sex_inference_input,
    output:
        gvcf=temp(DV_GVCF_TMPL),
        tbi=temp(DV_GVCF_TBI_TMPL),
    params:
        model_type=lambda w: _dv_model_type(w.species, w.assembly, w.sample_id),
        model_dir=lambda w: {
            "WGS": "/opt/models/wgs",
            "WES": "/opt/models/wes",
            "PACBIO": "/opt/models/pacbio",
            "ONT_R104": "/opt/models/ont_r104",
            "HYBRID_PACBIO_ILLUMINA": "/opt/models/hybrid_pacbio_illumina",
        }[_dv_model_type(w.species, w.assembly, w.sample_id)],
        haploid_contigs=_dv_haploid_contigs_flag,
    threads: cpu_func("deepvariant")
    resources:
        mem_mb=mem_func("deepvariant"),
        runtime=time_func("deepvariant")
    container: CONTAINERS["deepvariant"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "E10_deepvariant.{sample_id}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.gvcf})" "$(dirname {log})"
        exec > {log} 2>&1

        # Raise thread limit (best-effort, applies to all child processes)
        ulimit -u unlimited 2>/dev/null || ulimit -u "$(ulimit -Hu)" 2>/dev/null || true

        echo "[GAME] DeepVariant (split-steps) for {wildcards.sample_id}"
        echo "[GAME] Model type: {params.model_type} ({params.model_dir})"
        # Validate model type — this rule has been tested only for these.
        case "{params.model_type}" in
            WGS|PACBIO|ONT_R104)
                echo "[GAME] Model type validated: {params.model_type}"
                ;;
            *)
                echo "[GAME] ERROR: Unsupported model type '{params.model_type}'"
                echo "[GAME] This rule supports: WGS, PACBIO, ONT_R104"
                echo "[GAME] WES, HYBRID_PACBIO_ILLUMINA, and others have not been validated."
                exit 1
                ;;
        esac
        echo "[GAME] Threads/shards: {threads}"
        echo "[GAME] ulimit -u (soft/hard): $(ulimit -Su) / $(ulimit -Hu)"

        WORK_DIR="$(game_get_workdir 150)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_deepvariant_{wildcards.species}_{wildcards.assembly}_{wildcards.sample_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        echo "[GAME] Working directory: $TEMP_DIR"

        DV_TMP="$TEMP_DIR/dv_tmp"
        mkdir -p "$DV_TMP"

        EXAMPLES="$DV_TMP/make_examples.tfrecord@{threads}.gz"
        GVCF_RECORDS="$DV_TMP/gvcf.tfrecord@{threads}.gz"
        CV_OUTPUT="$DV_TMP/call_variants_output.tfrecord.gz"
        MODEL_JSON="{params.model_dir}/model.example_info.json"

        # Shared JAX compile cache (helps small-model in make_examples
        # AND reused in call_variants).
        export JAX_COMPILATION_CACHE_DIR="$TEMP_DIR/jax_cache"
        mkdir -p "$JAX_COMPILATION_CACHE_DIR"

        # ==========================================================
        # Mirror run_deepvariant.py decisions that depend on the
        # model.example_info.json and on the model type.
        # ==========================================================

        # Small model: only enable if the JSON declares one
        # (matches run_deepvariant.py:_use_small_model).
        SMALL_MODEL_ME_FLAG=""
        SMALL_MODEL_PP_FLAG=""
        if python3 - "$MODEL_JSON" <<'PYEOF' 2>/dev/null
import json, sys
try:
    with open(sys.argv[1]) as f:
        info = json.load(f)
    flags = info.get('flags_for_calling') or dict()
    sys.exit(0 if flags.get('trained_small_model_path') else 1)
except Exception:
    sys.exit(1)
PYEOF
        then
            SMALL_MODEL_ME_FLAG="--call_small_model_examples"
            SMALL_MODEL_PP_FLAG="--small_model_cvo_records $DV_TMP/make_examples_call_variant_outputs.tfrecord@{threads}.gz"
            echo "[GAME] Small model: ENABLED (declared in model.example_info.json)"
        else
            echo "[GAME] Small model: disabled (not declared in model.example_info.json)"
        fi

        # Phasing: enable for long-read models
        # (matches run_deepvariant.py:_should_phase_vcf).
        PHASE_ME_FLAGS=""
        PHASE_PP_FLAGS=""
        case "{params.model_type}" in
            PACBIO|ONT_R104)
                PHASE_TSV="$DV_TMP/read-phasing_debug@{threads}.tsv"
                PHASE_ME_FLAGS="--output_local_read_phasing $PHASE_TSV --output_phase_info"
                PHASE_PP_FLAGS="--phased_reads_input_path $PHASE_TSV"
                echo "[GAME] Phasing: ENABLED for {params.model_type}"
                ;;
            *)
                echo "[GAME] Phasing: disabled for {params.model_type}"
                ;;
        esac

        export HAPLOID_FLAG="{params.haploid_contigs}"
        if [ -n "$HAPLOID_FLAG" ]; then
            echo "[GAME] Haploid contigs: $HAPLOID_FLAG"
        else
            echo "[GAME] Haploid contigs: (none — all sites called diploid)"
        fi

        # ==========================================================
        # STEP 1: make_examples  --  THREAD-PINNED, PARALLEL
        # ==========================================================
        echo "[GAME] ===== make_examples starting ====="
        (
            export OPENBLAS_NUM_THREADS=1
            export MKL_NUM_THREADS=1
            export OMP_NUM_THREADS=1
            export NUMEXPR_NUM_THREADS=1
            export TF_NUM_INTEROP_THREADS=1
            export TF_NUM_INTRAOP_THREADS=1
            export TF_ENABLE_ONEDNN_OPTS=0
            export TF_CPP_MIN_LOG_LEVEL=2

            time seq 0 $(( {threads} - 1 )) | parallel -q -j {threads} --halt now,fail=1 --line-buffer \
                /opt/deepvariant/bin/make_examples \
                    --mode calling \
                    --ref "{input.ref}" \
                    --reads "{input.bam}" \
                    --examples "$EXAMPLES" \
                    --gvcf "$GVCF_RECORDS" \
                    --checkpoint "{params.model_dir}" \
                    --checkpoint_json "$MODEL_JSON" \
                    $SMALL_MODEL_ME_FLAG \
                    $PHASE_ME_FLAGS \
                    ${{HAPLOID_FLAG:-}} \
                    --sample_name "{wildcards.sample_id}" \
                    --task {{}}
        )
        echo "[GAME] ===== make_examples done ====="

        # ==========================================================
        # STEP 2: call_variants  --  UNPINNED, SINGLE PROCESS
        # ==========================================================
        echo "[GAME] ===== call_variants starting ====="
        (
            export TF_CPP_MIN_LOG_LEVEL=2
            time /opt/deepvariant/bin/call_variants \
                --outfile "$CV_OUTPUT" \
                --examples "$EXAMPLES" \
                --checkpoint "{params.model_dir}"
        )
        echo "[GAME] ===== call_variants done ====="

        # ==========================================================
        # STEP 3: postprocess_variants  --  MULTI-CPU, LIGHTWEIGHT
        # ==========================================================
        echo "[GAME] ===== postprocess_variants starting ====="
        (
            export TF_CPP_MIN_LOG_LEVEL=2
            time /opt/deepvariant/bin/postprocess_variants \
                --ref "{input.ref}" \
                --infile "$CV_OUTPUT" \
                --outfile "$TEMP_DIR/output.vcf.gz" \
                --nonvariant_site_tfrecord_path "$GVCF_RECORDS" \
                --gvcf_outfile "$TEMP_DIR/output.g.vcf.gz" \
                --checkpoint_json "$MODEL_JSON" \
                $SMALL_MODEL_PP_FLAG \
                $PHASE_PP_FLAGS \
                ${{HAPLOID_FLAG:-}} \
                --cpus 4 \
                --sample_name "{wildcards.sample_id}"
        )
        echo "[GAME] ===== postprocess_variants done ====="

        echo "[GAME] Copying gVCF to output..."
        cp "$TEMP_DIR/output.g.vcf.gz"     "{output.gvcf}"
        cp "$TEMP_DIR/output.g.vcf.gz.tbi" "{output.tbi}"
        echo "[GAME] Done"
        '''


rule E11_dv_to_bcf:
    """Convert DeepVariant gVCF to BCF.
    block mode:    gVCF → BCF  (reference blocks preserved)
    basepair mode: gVCF → expand ref-blocks to per-site records → BCF
    """
    input:
        gvcf=DV_GVCF_TMPL,
        tbi=DV_GVCF_TBI_TMPL,
        ref=rules.D00_ref_fa_canonical.output.fa,
        fai=rules.E00_ref_faidx_for_calling.output.fai,
    output:
        bcf=temp(DV_BCF_TMPL),
        csi=temp(DV_BCF_CSI_TMPL),
    params:
        bp_mode=_BP_MODE,
    threads: cpu_func("bcftools_concat")
    resources:
        mem_mb=mem_func("bcftools_concat"),
        runtime=time_func("bcftools_concat")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "E11_dv_to_bcf.{sample_id}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.bcf})" "$(dirname {log})"

        exec > {log} 2>&1

        echo "[GAME] Converting DeepVariant gVCF to BCF for {wildcards.sample_id}"


        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 50)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_dv2bcf_{wildcards.species}_{wildcards.assembly}_{wildcards.sample_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        echo "[GAME] Working directory: $TEMP_DIR"

        export HTS_COMPRESSION_LEVEL=1

        if [ "{params.bp_mode}" = "True" ]; then
            echo "[GAME] Resolution: basepair (expanding reference blocks)"
            bcftools convert --gvcf2vcf "{input.gvcf}" \
                --fasta-ref "{input.ref}" | \
                bcftools view -A -O b -o "output.bcf" --threads {threads}
        else
            echo "[GAME] Resolution: block (preserving reference blocks)"
            bcftools view -A -O b -o "output.bcf" "{input.gvcf}" --threads {threads}
        fi

        echo "[GAME] Creating CSI index for BCF"
        bcftools index -f --csi --threads {threads} "output.bcf"

        # Copy results to final location
        echo "[GAME] Copying output to final location..."
        cp "output.bcf"     "{output.bcf}"
        cp "output.bcf.csi" "{output.csi}"

        echo "[GAME] Done"

        '''


# ===============================================================================
#  ROUTER RULES — caller-specific intermediates → final output paths
# ===============================================================================
#  These rules exist so that downstream modules (F, G, H) always see the same
#  filenames ({sid}[.bp].g.vcf.gz, {sid}[.bp].raw.bcf) regardless of which caller
#  produced them.  The input functions select the GATK or DV intermediate based
#  on _caller_for_sample(), which inspects the per-sample tech count at DAG time.

rule E20_route_gvcf:
    """Route caller-specific gVCF to final output path."""
    input:
        gvcf=lambda w: _routed_caller_file(w, DV_GVCF_TMPL, GATK_GVCF_TMPL),
        tbi=lambda w: _routed_caller_file(w, DV_GVCF_TBI_TMPL, GATK_GVCF_TBI_TMPL),
    output:
        gvcf=(FINAL_GVCF_TMPL if _KEEP_BP_GVCF else temp(FINAL_GVCF_TMPL)),
        tbi=(FINAL_GVCF_TBI_TMPL if _KEEP_BP_GVCF else temp(FINAL_GVCF_TBI_TMPL))
    threads: 1
    resources:
        mem_mb=mem_func("split_intervals"),
        runtime=time_func("split_intervals")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "E20_route_gvcf.{sample_id}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        echo "[GAME] Routing gVCF for {wildcards.sample_id}" > {log}
        cp "{input.gvcf}" "{output.gvcf}"
        cp "{input.tbi}"  "{output.tbi}"
        echo "[GAME] Done" >> {log}
        '''


rule E21_route_bcf:
    """Route caller-specific BCF to final output path."""
    input:
        bcf=lambda w: _routed_caller_file(w, DV_BCF_TMPL, GATK_BCF_TMPL),
        csi=lambda w: _routed_caller_file(w, DV_BCF_CSI_TMPL, GATK_BCF_CSI_TMPL),
    output:
        bcf=FINAL_BCF_TMPL,
        csi=FINAL_BCF_CSI_TMPL
    threads: 1
    resources:
        mem_mb=mem_func("split_intervals"),
        runtime=time_func("split_intervals")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "E21_route_bcf.{sample_id}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        echo "[GAME] Routing BCF for {wildcards.sample_id}" > {log}
        cp "{input.bcf}" "{output.bcf}"
        cp "{input.csi}" "{output.csi}"
        echo "[GAME] Done" >> {log}
        '''


# ===============================================================================
#  JOINT GENOTYPING via GLnexus (CALLER: deep + MERGE_TECH: off + JOINT_GENO: on)
# ===============================================================================

def _joint_gvcfs_for_assembly(w):
    """Collect all per-sample gVCFs under a species/assembly for joint genotyping."""
    sp, asm = w.species, w.assembly
    gvcfs = []
    try:
        asm_data = samples_config["sp_name"][sp]["asm_id"][asm]
        if "sample_id" not in asm_data:
            return gvcfs
        for sid in asm_data["sample_id"]:
            if sid in ("None", None):
                continue
            sample_data = asm_data["sample_id"][sid]
            # Check sample has reads
            has_reads = False
            for rt_key in sample_data.get("read_type", {}):
                if rt_key not in ("None", None):
                    _, _, reads = _get_sample_node(sp, asm, sid, rt_key)
                    if reads:
                        has_reads = True
                        break
            if has_reads:
                gvcfs.append(
                    FINAL_GVCF_TMPL.format(
                        species=sp, assembly=asm, sample_id=sid
                    )
                )
    except (KeyError, TypeError, AttributeError):
        pass
    return gvcfs


rule EJ10_glnexus_joint:
    """Joint genotyping of all samples under one assembly using GLnexus.
    Only runs when JOINT_GENO: on, CALLER: deep, MERGE_TECH: off,
    and there are ≥2 samples for this assembly.
    """
    input:
        gvcfs=_joint_gvcfs_for_assembly,
        ref=rules.D00_ref_fa_canonical.output.fa,
        fai=rules.E00_ref_faidx_for_calling.output.fai,
    output:
        bcf=JOINT_BCF_TMPL,
        csi=JOINT_BCF_CSI_TMPL,
    threads: cpu_func("glnexus")
    resources:
        mem_mb=mem_func("glnexus"),
        runtime=time_func("glnexus")
    container: CONTAINERS["glnexus"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "logs", "EJ10_glnexus_joint.{assembly}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.bcf})" "$(dirname {log})"

        exec > {log} 2>&1

        N_SAMPLES=$(echo {input.gvcfs} | wc -w)
        echo "[GAME] GLnexus joint genotyping for {wildcards.assembly}"
        echo "[GAME] Samples: $N_SAMPLES"
        echo "[GAME] Threads: {threads}"


        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 100)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_glnexus_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        echo "[GAME] Working directory: $TEMP_DIR"

        # Convert mem_mb to GB for GLnexus (floor, minimum 4)
        MEM_GB=$(( {resources.mem_mb} / 1024 ))
        [ "$MEM_GB" -lt 4 ] && MEM_GB=4

        # Run GLnexus
        glnexus_cli \
            --config DeepVariant \
            --dir "$TEMP_DIR/glnexus_db" \
            --threads {threads} \
            --mem-gbytes "$MEM_GB" \
            {input.gvcfs} \
        > "$TEMP_DIR/joint.bcf"

        echo "[GAME] Creating CSI index"
        bcftools index -f --csi --threads {threads} "$TEMP_DIR/joint.bcf"

        # Copy results to final location
        echo "[GAME] Copying output to final location..."
        cp "$TEMP_DIR/joint.bcf"     "{output.bcf}"
        cp "$TEMP_DIR/joint.bcf.csi" "{output.csi}"

        echo "[GAME] Joint genotyping complete ($N_SAMPLES samples)"

        '''
