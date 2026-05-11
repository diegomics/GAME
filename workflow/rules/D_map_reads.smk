# ===============================================================================
#  GAME - Mapping Rules
# ===============================================================================

# consistent idx/base separator for BAM filenames
_SEP = "_"

# -------------------------------------------------------------------------------
#  HELPER FUNCTIONS
# -------------------------------------------------------------------------------

def _norm_rt(rt):
    return normalize_read_type(rt)

def _maybe_temp(path, condition):
    """Wrap path in Snankemake's directive temp() if condition is True to automatically delete intermediate files"""
    return temp(path) if condition else path

def _idx_dir(sp, asm):
    return os.path.join(config["OUT_FOLDER"], "GAME_results", sp, asm, "idx")

def _bams_dir(sp, asm, sid):
    return os.path.join(config["OUT_FOLDER"], "GAME_results", sp, asm, "samples", sid, "BAMs")

def _ref_fa(sp, asm):
    return samples_config["sp_name"][sp]["asm_id"][asm]["asm_file"]

def _platform_from_rt(rt):
    rt = _norm_rt(rt)
    if rt in {"illumina","pe","sr","10x","tenx","10xgenomics","hic","hi-c","linked","linkedreads"}:
        return "ILLUMINA"
    if rt in {"hifi","pacbio","pacbiohifi","pb_hifi","pbhifi"}:
        return "PACBIO"
    if rt in {"ont","nanopore","oxfordnanopore"}:
        return "ONT"
    return "ILLUMINA"

def _rg_string(w, idx, base):
    # read group: ID + SM + PL (+ LB, PU for disambiguation)
    sm = w.sample_id
    rt = _norm_rt(w.read_type)
    pl = _platform_from_rt(rt)
    lb = f"{sm}.{rt}"
    pu = f"{sm}.P{idx}"
    rgid = f"{sm}.{rt}.P{idx}.{base}"
    return f'@RG\\tID:{rgid}\\tSM:{sm}\\tPL:{pl}\\tLB:{lb}\\tPU:{pu}'

def _minimap_preset(rt):
    rt = _norm_rt(rt)
    if rt in {"hifi","pacbio","pacbiohifi","pb_hifi","pbhifi"}:
        return "map-hifi"
    if rt in {"ont","nanopore","oxfordnanopore"}:
        return "map-ont"
    return "map-hifi"

def _reads_for_group(w):
    """
    Return (kind, paths) for mapping inputs of the requested group,
    using the standardized files in reads/.
    """
    sp, asm, sid = w.species, w.assembly, w.sample_id
    # Pass the read_type from wildcards to get the correct data
    rt, trim, read_dict = _get_sample_node(sp, asm, sid, w.read_type)
    if rt is None:
        raise ValueError(f"No valid read data for {sp}/{asm}/{sid}/{w.read_type}")
    
    # Verify the read type matches
    if rt != normalize_read_type(w.read_type):
        raise ValueError(f"read_type mismatch: got {rt}, expected {normalize_read_type(w.read_type)}")
    
    for grp in _enumerate_groups(rt, read_dict):
        if str(grp["idx"]) == str(w.idx) and grp["base"] == w.base:
            outs = _outputs_for_group(sp, asm, sid, grp, trim)
            if grp["kind"] == "long":
                return "long", (outs[0],)
            else:
                return "pe", (outs[0], outs[1])
    raise ValueError(f"Could not find group Path{w.idx}_{w.base} for {sp}/{asm}/{sid}/{w.read_type}")

def _bwa_prefix(sp, asm):
    return os.path.join(_idx_dir(sp, asm), asm)  # prefix

def _minimap_index_path(sp, asm):
    return os.path.join(_idx_dir(sp, asm), f"{asm}.mmi")

def _final_group_bam(sp, asm, sid, grp):
    rt = grp["rt"]; idx = grp["idx"]; base = grp["base"]
    bdir = _bams_dir(sp, asm, sid)
    stem = f"{rt}_Path{idx}{_SEP}{base}"
    if grp["kind"] == "pe":
        return os.path.join(
            bdir,
            f"{stem}.markdup.bam" if config.get("DEDUP_PE", True) else f"{stem}.sorted.bam",
        )
    else:
        return os.path.join(
            bdir,
            f"{stem}.markdup.bam" if config.get("DEDUP_LONG", False) else f"{stem}.sorted.bam",
        )

def _qc_dir(sp, asm, sid):
    return os.path.join(_bams_dir(sp, asm, sid), "qc")

def _reads_for_pe_group(w):
    kind, paths = _reads_for_group(w)
    if kind != "pe":
        raise ValueError(f"Requested PE group but found {kind} for {w.species}/{w.assembly}/{w.sample_id}")
    return paths  # (r1, r2)

def _reads_for_long_group(w):
    kind, paths = _reads_for_group(w)
    if kind != "long":
        raise ValueError(f"Requested long group but found {kind} for {w.species}/{w.assembly}/{w.sample_id}")
    return paths  # (r1,)

def _qc_summaries_for_assembly(w):
    sp, asm = w.species, w.assembly
    ins = []
    for _sp, _asm, sid in all_samples_iter():
        if _sp != sp or _asm != asm:
            continue
        
        # Get all read types for this sample
        try:
            sample_data = samples_config["sp_name"][_sp]["asm_id"][_asm]["sample_id"][sid]
            if not sample_data.get("read_type"):
                continue
            
            # Process each read type
            for rt_key in sample_data["read_type"]:
                if rt_key == "None" or rt_key is None:
                    continue
                
                rt_norm = normalize_read_type(rt_key)
                _, _, reads = _get_sample_node(_sp, _asm, sid, rt_key)
                
                if reads:  # Only if there are actual read files
                    qcdir = os.path.join(config["OUT_FOLDER"], "GAME_results", sp, asm, 
                                       "samples", sid, "BAMs", "qc")
                    ins.append(os.path.join(qcdir, f"{sid}.{rt_norm}.summary.txt"))
        except (KeyError, TypeError, AttributeError):
            continue
    
    return ins


# Tech merge decision logic
def _get_tech_bams_for_sample(w):
    """Get all tech-specific merged BAMs for a sample"""
    sp, asm, sid = w.species, w.assembly, w.sample_id
    bams = []
    
    try:
        sample_data = samples_config["sp_name"][sp]["asm_id"][asm]["sample_id"][sid]
        if not sample_data.get("read_type"):
            return bams
        
        for rt_key in sample_data["read_type"]:
            if rt_key == "None" or rt_key is None:
                continue
            
            rt_norm = normalize_read_type(rt_key)
            _, _, reads = _get_sample_node(sp, asm, sid, rt_key)
            
            if reads:  # Only if there are actual read files
                bam_path = os.path.join(
                    config["OUT_FOLDER"], "GAME_results", sp, asm, 
                    "samples", sid, "BAMs", f"{sid}.{rt_norm}.merged.bam"
                )
                # Don't check if file exists - just add the expected path
                # Snakemake will ensure these are created before running this rule
                bams.append(bam_path)
    except (KeyError, TypeError, AttributeError):
        pass
    
    return bams


# -------------------------------------------------------------------------------
#  INDEXING
# ===============================================================================
# one canonical reference FASTA in idx/: symlink if uncompressed, else decompress once
rule D00_ref_fa_canonical:
    input:
        ref=lambda w: _ref_fa(w.species, w.assembly)
    output:
        fa=os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "idx", "{assembly}.fa")
    container: CONTAINERS["game_base"]
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.fa})"
        if [[ "{input.ref}" == *.gz ]]; then
            gunzip -c "{input.ref}" > "{output.fa}"
        else
            ln -sf "$(realpath "{input.ref}")" "{output.fa}"
        fi
        '''

rule D01_minimap2_index:
    """Build minimap2 index once per assembly."""
    input:
        fa=rules.D00_ref_fa_canonical.output.fa
    output:
        mmi=os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "idx", "{assembly}.mmi")
    threads: cpu_func("minimap2")
    resources:
        mem_mb=mem_func("minimap2"),
        runtime=time_func("minimap2")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "logs", "D01_minimap2_index.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {log})
        
        exec > {log} 2>&1
        
        echo "[GAME] Creating minimap2 index for {wildcards.assembly}"
        

        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 50)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_mm2index_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        
        echo "[GAME] Copying assembly to temp directory..."
        cp "{input.fa}" local.fa
        
        echo "[GAME] Running minimap2 -d..."
        minimap2 -d local.mmi local.fa
        
        echo "[GAME] Copying index back to output..."
        mkdir -p "$(dirname {output.mmi})"
        cp local.mmi "{output.mmi}"
        
        echo "[GAME] Done"
        
        '''

rule D01_bwa_mem2_index:
    """Create bwa-mem2 index once per assembly (prefix-based)."""
    input:
        fa=rules.D00_ref_fa_canonical.output.fa
    output:
        flag=os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "idx", "bwa.done")
    params:
        prefix=lambda w: os.path.join(_idx_dir(w.species, w.assembly), w.assembly),
        idx_dir=lambda w: _idx_dir(w.species, w.assembly)
    threads: cpu_func("bwa_index")
    resources:
        mem_mb=mem_func("bwa_index"),
        runtime=time_func("bwa_index")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "logs", "D01_bwa_mem2_index.log")
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {log})
        
        exec > {log} 2>&1
        
        echo "[GAME] Creating bwa-mem2 index for {wildcards.assembly}"
        

        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 50)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_bwaindex_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        cd "$TEMP_DIR"
        
        echo "[GAME] Copying assembly to temp directory..."
        cp "{input.fa}" local.fa
        
        echo "[GAME] Running bwa-mem2 index..."
        bwa-mem2 index -p local local.fa
        
        echo "[GAME] Copying index files back to output..."
        mkdir -p "{params.idx_dir}"
        for f in local.*; do
            # Skip the .fa file itself
            if [[ "$f" != "local.fa" ]]; then
                ext="${{f#local}}"
                cp "$f" "{params.prefix}$ext"
                echo "[GAME] Copied: {params.prefix}$ext"
            fi
        done
        
        touch "{output.flag}"
        
        echo "[GAME] Done"
        
        '''


# -------------------------------------------------------------------------------
#  LONG READS MAPPING
# ===============================================================================

rule D02_map_long_reads:
    """minimap2 mapping of one long-read group -> sorted BAM + index + flagstat."""
    input:
        r1=lambda w: _reads_for_long_group(w)[0],
        idx=rules.D01_minimap2_index.output.mmi
    output:
        bam=_maybe_temp(
            os.path.join(
                config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
                "samples", "{sample_id}", "BAMs",
                "{read_type}_Path{idx}" + _SEP + "{base}.sorted.bam"
            ),
            config.get("DEDUP_LONG", False)
        ),
        bai=_maybe_temp(
            os.path.join(
                config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
                "samples", "{sample_id}", "BAMs",
                "{read_type}_Path{idx}" + _SEP + "{base}.sorted.bam.bai"
            ),
            config.get("DEDUP_LONG", False)
        )
    wildcard_constraints:
        read_type=r"(hifi|ont)",
        idx=r"\d+"
    threads: cpu_func("minimap2")
    resources:
        mem_mb=mem_func("minimap2"),
        runtime=time_func("minimap2")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}",
            "logs", "D02_map_long_reads.{read_type}.Path{idx}" + _SEP + "{base}.log"
        )
    run:
        r1 = input.r1
        preset = _minimap_preset(wildcards.read_type)
        rg = _rg_string(wildcards, wildcards.idx, wildcards.base)
        ref_mmi = input.idx
        os.makedirs(os.path.dirname(output.bam), exist_ok=True)
        os.makedirs(os.path.dirname(str(log)), exist_ok=True)
        # Streaming pipe - no temp dir needed
        shell(
            'set -euo pipefail\n'
            'minimap2 -t {t} -a -x {preset} -R "{rg}" "{mmi}" "{r1}" | '
            'samtools sort -@ {t} -o "{bam}"\n'
            'samtools index -@ {t} "{bam}"\n'
            'samtools flagstat -@ {t} "{bam}" >> "{log}"\n'
            .format(t=threads, preset=preset, rg=rg, mmi=ref_mmi, r1=r1, bam=output.bam, log=str(log))
        )

rule D03_markdup_long_reads:
    """Optional: mark duplicates for long reads group (enabled via DEDUP_LONG=true)."""
    input:
        bam=rules.D02_map_long_reads.output.bam
    output:
        bam=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{read_type}_Path{idx}" + _SEP + "{base}.markdup.bam",
        ),
        bai=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{read_type}_Path{idx}" + _SEP + "{base}.markdup.bam.bai",
        ),
        metrics=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{read_type}_Path{idx}" + _SEP + "{base}.markdup_metrics.txt",
        )
    wildcard_constraints:
        read_type=r"(hifi|ont)",
        idx=r"\d+"
    threads: cpu_func("markdup")
    resources:
        mem_mb=mem_func("markdup"),
        runtime=time_func("markdup")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}",
            "logs", "D03_markdup_long_reads.{read_type}.Path{idx}" + _SEP + "{base}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.bam})" "$(dirname {log})"
        
        exec > "{log}" 2>&1
        
        echo "[GAME] Starting duplicate marking for {input.bam}"
        

        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 100)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_markduplong_{wildcards.species}_{wildcards.sample_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        
        echo "[GAME] Using temp directory: $TEMP_DIR"
        
        echo "[GAME] Running sambamba markdup"
        sambamba markdup -p -t {threads} --tmpdir "$TEMP_DIR" \
            "{input.bam}" "{output.bam}" 2>&1 | tee "{output.metrics}"
        
        echo "[GAME] Indexing BAM"
        samtools index -@ {threads} "{output.bam}"
        
        echo "[GAME] Done"
        
        '''


# -------------------------------------------------------------------------------
#  SHORT READS MAPPING
# ===============================================================================

rule D02_map_pe_reads:
    """bwa-mem2 mapping of one PE group -> sorted BAM + index + flagstat."""
    input:
        r1=lambda w: _reads_for_pe_group(w)[0],
        r2=lambda w: _reads_for_pe_group(w)[1],
        idx=rules.D01_bwa_mem2_index.output.flag
    output:
        bam=_maybe_temp(
            os.path.join(
                config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
                "samples", "{sample_id}", "BAMs",
                "{read_type}_Path{idx}" + _SEP + "{base}.sorted.bam"
            ),
            config.get("DEDUP_PE", True)
        ),
        bai=_maybe_temp(
            os.path.join(
                config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
                "samples", "{sample_id}", "BAMs",
                "{read_type}_Path{idx}" + _SEP + "{base}.sorted.bam.bai"
            ),
            config.get("DEDUP_PE", True)
        )
    wildcard_constraints:
        read_type=r"(illumina|10x)",
        idx=r"\d+"
    threads: cpu_func("bwa_mem")
    resources:
        mem_mb=mem_func("bwa_mem"),
        runtime=time_func("bwa_mem")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}",
            "logs", "D02_map_pe_reads.{read_type}.Path{idx}" + _SEP + "{base}.log"
        )
    run:
        r1, r2 = input.r1, input.r2
        ref_prefix = _bwa_prefix(wildcards.species, wildcards.assembly)
        rg = _rg_string(wildcards, wildcards.idx, wildcards.base)
        os.makedirs(os.path.dirname(output.bam), exist_ok=True)
        os.makedirs(os.path.dirname(str(log)), exist_ok=True)
        # Streaming pipe - no temp dir needed
        shell(
            'set -euo pipefail\n'
            '(\n'
            '  bwa-mem2 mem -t {t} -M -R "{rg}" "{ref}" "{r1}" "{r2}" | \n'
            '  samtools sort -@ {t} -o "{bam}"\n'
            ') 2>> "{log}"\n'
            'samtools index -@ {t} "{bam}" 2>> "{log}"\n'
            'samtools flagstat -@ {t} "{bam}" >> "{log}"\n'
            .format(t=threads, rg=rg, ref=ref_prefix, r1=r1, r2=r2,
                    bam=output.bam, log=str(log))
        )

rule D03_markdup_pe_reads:
    """Mark duplicates for PE group (enabled by DEDUP_PE)."""
    input:
        bam=rules.D02_map_pe_reads.output.bam
    output:
        bam=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{read_type}_Path{idx}" + _SEP + "{base}.markdup.bam",
        ),
        bai=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{read_type}_Path{idx}" + _SEP + "{base}.markdup.bam.bai",
        ),
        metrics=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{read_type}_Path{idx}" + _SEP + "{base}.markdup_metrics.txt",
        )
    wildcard_constraints:
        read_type=r"(illumina|10x)",
        idx=r"\d+"
    threads: cpu_func("markdup")
    resources:
        mem_mb=mem_func("markdup"),
        runtime=time_func("markdup")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}",
            "logs", "D03_markdup_pe_reads.{read_type}.Path{idx}" + _SEP + "{base}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.bam})" "$(dirname {log})"
        
        exec > "{log}" 2>&1
        
        echo "[GAME] Starting duplicate marking for {input.bam}"
        

        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 100)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_markdup_{wildcards.species}_{wildcards.sample_id}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT
        
        echo "[GAME] Using temp directory: $TEMP_DIR"
        
        echo "[GAME] Running sambamba markdup"
        sambamba markdup -p -t {threads} --tmpdir "$TEMP_DIR" \
            "{input.bam}" "{output.bam}" 2>&1 | tee "{output.metrics}"
        
        echo "[GAME] Indexing BAM"
        samtools index -@ {threads} "{output.bam}"
        
        echo "[GAME] Done"
        
        '''



# ----------------------------- merge per sample & read type -----------------------------

def _sample_group_bams_by_type(w):
    """All final group BAMs for this sample&read_type."""
    sp, asm, sid = w.species, w.assembly, w.sample_id
    rt_req = normalize_read_type(w.read_type)
    
    # Get data for the specific read type requested
    rt, trim, reads = _get_sample_node(sp, asm, sid, w.read_type)
    if rt is None or rt != rt_req:
        return []
    
    outs = []
    for grp in _enumerate_groups(rt, reads):
        outs.append(_final_group_bam(sp, asm, sid, grp))
    return outs

rule D04_merge_sample_readtype_bams:
    input:
        _sample_group_bams_by_type
    output:
        bam=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{sample_id}.{read_type}.merged.bam",
        ),
        bai=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{sample_id}.{read_type}.merged.bam.bai",
        ),
        stat=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}", "BAMs",
            "{sample_id}.{read_type}.merged.flagstat.txt",
        )
    threads: cpu_func("bwa_mem")
    resources:
        mem_mb=mem_func("bwa_mem"),
        runtime=time_func("bwa_mem")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "samples", "{sample_id}",
            "logs", "D04_merge_bams.{sample_id}.{read_type}.log"
        )
    run:
        os.makedirs(os.path.dirname(output.bam), exist_ok=True)
        if len(input) == 0:
            raise ValueError(f"No BAMs to merge for {wildcards.sample_id}/{wildcards.read_type}")
        if len(input) == 1:
            import os as _os
            rel = _os.path.relpath(input[0], _os.path.dirname(output.bam))
            shell('ln -sf "{src}" "{dst}"'.format(src=rel, dst=output.bam))
        else:
            shell('samtools merge -@ {t} -o "{out}" {ins}'.format(
                t=threads, out=output.bam, ins=" ".join(f'"{x}"' for x in input)
            ))
        shell(
            'set -euo pipefail\n'
            'samtools index -@ {t} "{bam}"\n'
            'samtools flagstat -@ {t} "{bam}" > "{stat}"\n'
            .format(t=threads, bam=output.bam, stat=output.stat)
        )

# ----------------------------- QC of merged BAM -----------------------------
rule D05_qc_merged_bam:
    """QC of merged BAM: flagstat, idxstats, coverage, concise summary (incl. dispersion & outliers)."""
    input:
        bam=rules.D04_merge_sample_readtype_bams.output.bam,
        bai=rules.D04_merge_sample_readtype_bams.output.bai
    output:
        flagstat=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "qc",
            "{sample_id}.{read_type}.flagstat.txt",
        ),
        idxstats=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "qc",
            "{sample_id}.{read_type}.idxstats.tsv",
        ),
        coverage=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "qc",
            "{sample_id}.{read_type}.coverage.tsv",
        ),
        summary=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "qc",
            "{sample_id}.{read_type}.summary.txt",
        ),
        summary_md=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "qc",
            "{sample_id}.{read_type}.summary.md",
        )
    params:
        qc_bam_summary = str(scripts_dir / "qc_bam_summary.py"),
        rt=lambda w: _norm_rt(w.read_type),   # canonical read type for summary
        disp_pct=0.95,                        # for Depth dispersion (95%)
        hi_mult=5.0                           # High-depth outliers threshold: depth >= hi_mult * WMD
    threads: cpu_func("bam_qc")
    resources:
        mem_mb=mem_func("bam_qc"),
        runtime=time_func("bam_qc")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs",
            "D05_qc_merged_bam.{sample_id}.{read_type}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.flagstat})" "$(dirname {log})"
        
        exec > "{log}" 2>&1
        
        echo "[GAME] QC for sample={wildcards.sample_id} read_type={params.rt}"
        
        # samtools outputs
        samtools flagstat -@ {threads} "{input.bam}" > "{output.flagstat}"
        samtools idxstats "{input.bam}" > "{output.idxstats}"
        samtools coverage "{input.bam}" > "{output.coverage}"

        # Primary mapped (exclude secondary+supplementary+QCfail+unmapped) = -F 2308
        PRIMARY=$(samtools view -@ {threads} -c -F 2308 "{input.bam}")
        # Primary mapped, excluding duplicates too = -F 3332
        PRIMARY_NODUP=$(samtools view -@ {threads} -c -F 3332 "{input.bam}")

        # Compute summary & dispersion metrics
        python {params.qc_bam_summary:q} \
            --sample "{wildcards.sample_id}" \
            --read-type "{params.rt}" \
            --coverage "{output.coverage}" \
            --flagstat "{output.flagstat}" \
            --primary "$PRIMARY" \
            --primary-nodup "$PRIMARY_NODUP" \
            --disp-pct {params.disp_pct} \
            --hi-mult {params.hi_mult} \
            --out "{output.summary}" \
            --md-out "{output.summary_md}"
        
        echo "[GAME] QC complete"
        '''


rule D06_aggregate_bam_qc:
    """Build a Markdown table combining per-sample BAM QC summaries."""
    input:
        _qc_summaries_for_assembly
    output:
        md=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "bam_qc.md"
        )
    threads: cpu_func("light_task")
    resources:
        mem_mb=mem_func("light_task"),
        runtime=time_func("light_task")
    container: CONTAINERS["game_base"]
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "logs", "D06_aggregate_bam_qc.log"
        )
    params:
        aggregate_bam_qc = str(scripts_dir / "aggregate_bam_qc.py")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "$(dirname {output.md})" "$(dirname {log})"
        python {params.aggregate_bam_qc:q} --out "{output.md}" {input} > "{log}" 2>&1
        '''


rule D07_tech_merge_decision:
    """Analyze coverage per tech and decide merge strategy"""
    input:
        bams=_get_tech_bams_for_sample,
        # Add QC summaries as explicit dependencies
        summaries=lambda w: [
            os.path.join(
                config["OUT_FOLDER"], "GAME_results", w.species, w.assembly,
                "samples", w.sample_id, "BAMs", "qc",
                f"{w.sample_id}.{normalize_read_type(rt_key)}.summary.txt"
            )
            for rt_key in samples_config["sp_name"][w.species]["asm_id"][w.assembly]["sample_id"][w.sample_id].get("read_type", {})
            if rt_key not in ["None", None]
        ]
    output:
        decision=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "merge_decision.json"
        )
    params:
        tech_merge_decision = str(scripts_dir / "tech_merge_decision.py"),
        # YAML parses unquoted on/off as True/False — normalise back to strings
        merge_mode={True: "on", False: "off"}.get(config.get("MERGE_TECH", "auto"), config.get("MERGE_TECH", "auto")),
        good_cov=config.get("GOOD_TECH_COV", 15),
        min_cov=config.get("MIN_TECH_COV", 1),
        min_frac=config.get("MIN_TECH_FRAC", 0.05),
        priority=config.get("DATA_PRIORITY", "hifi>illumina>ont")
    container: CONTAINERS["game_base"]
    threads: cpu_func("light_task")
    resources:
        mem_mb=mem_func("light_task"),
        runtime=time_func("light_task")
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "D07_merge_tech.{sample_id}.log"
        )
    shell:
        """
        python {params.tech_merge_decision:q} \
            --bams {input.bams} \
            --sample-id {wildcards.sample_id} \
            --merge-mode {params.merge_mode} \
            --good-cov {params.good_cov} \
            --min-cov {params.min_cov} \
            --min-frac {params.min_frac} \
            --priority "{params.priority}" \
            --output {output.decision} \
            --log {log}
        """

rule D08_merge_tech_bams:
    """Merge multiple technology BAMs into one, or symlink if single tech chosen"""
    input:
        decision=rules.D07_tech_merge_decision.output.decision,
        bams=_get_tech_bams_for_sample
    output:
        bam=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "{sample_id}.merged.bam"
        ),
        bai=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "BAMs", "{sample_id}.merged.bam.bai"
        )
    container: CONTAINERS["game_base"]
    threads: cpu_func("bwa_mem")
    resources:
        mem_mb=mem_func("bwa_mem"),
        runtime=time_func("bwa_mem")
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "samples", "{sample_id}", "logs", "D08_merge_tech_bams.{sample_id}.log"
        )
    shell:
        r'''
        set -euo pipefail
        mkdir -p $(dirname {log})
        
        exec > "{log}" 2>&1
        
        echo "[GAME] Starting tech BAM merge for {wildcards.sample_id}"
        
        # Read the decision file
        MODE=$(python -c "import json; d=json.load(open('{input.decision}')); print(d.get('mode', ''))")
        
        if [[ "${{MODE}}" == "merge" ]]; then
            echo "[GAME] Mode is merge, combining BAMs"
            # Get techs to merge
            TECHS=$(python -c "import json; d=json.load(open('{input.decision}')); print(' '.join(d.get('techs_to_merge', [])))")
            
            # Build list of BAMs to merge
            BAMS_TO_MERGE=""
            for BAM in {input.bams}; do
                for TECH in ${{TECHS}}; do
                    if [[ "${{BAM}}" == *".${{TECH}}.merged.bam" ]]; then
                        BAMS_TO_MERGE="${{BAMS_TO_MERGE}} ${{BAM}}"
                        break
                    fi
                done
            done
            
            # Count BAMs
            NUM_BAMS=$(echo "${{BAMS_TO_MERGE}}" | wc -w)
            
            if [[ ${{NUM_BAMS}} -gt 1 ]]; then
                echo "[GAME] Merging ${{NUM_BAMS}} BAMs: ${{BAMS_TO_MERGE}}"
                sambamba merge -t {threads} {output.bam} ${{BAMS_TO_MERGE}}
                samtools index -@ {threads} {output.bam}
            else
                echo "[GAME] Warning: Expected multiple BAMs but found ${{NUM_BAMS}}"
                ln -sf ${{BAMS_TO_MERGE}} {output.bam}
                samtools index -@ {threads} {output.bam}
            fi
        elif [[ "${{MODE}}" == "single" ]]; then
            echo "[GAME] Mode is single, symlinking chosen tech BAM"
            # Get the chosen tech
            CHOSEN=$(python -c "import json; d=json.load(open('{input.decision}')); print(d.get('chosen', ''))")
            
            if [[ -z "${{CHOSEN}}" ]]; then
                echo "[GAME] ERROR: No chosen tech in decision file"
                exit 1
            fi
            
            # Find the corresponding BAM
            CHOSEN_BAM=""
            for BAM in {input.bams}; do
                if [[ "${{BAM}}" == *".${{CHOSEN}}.merged.bam" ]]; then
                    CHOSEN_BAM="${{BAM}}"
                    break
                fi
            done
            
            if [[ -z "${{CHOSEN_BAM}}" ]]; then
                echo "[GAME] ERROR: Could not find BAM for chosen tech: ${{CHOSEN}}"
                exit 1
            fi
            
            echo "[GAME] Symlinking ${{CHOSEN_BAM}} to {output.bam}"
            ln -sf "$(basename "${{CHOSEN_BAM}}")" {output.bam}
            
            if [[ -f "${{CHOSEN_BAM}}.bai" ]]; then
                ln -sf "$(basename "${{CHOSEN_BAM}}.bai")" {output.bai}
            else
                samtools index -@ {threads} {output.bam}
            fi
        else
            echo "[GAME] ERROR: Unknown mode: ${{MODE}}"
            exit 1
        fi
        
        echo "[GAME] Done"
        '''
