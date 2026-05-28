# ===============================================================================
# GAME - Assembly Stats Rules
# by Diego De Panis, 2026
# note: AI tools may have been used to improve, clean and/or comment this version of the code
# ===============================================================================

# -------------------------------------------------------------------------------
#  HELPER FUNCTIONS
# -------------------------------------------------------------------------------

# get assembly path
#def get_assembly_path(species, assembly):
#    return samples_config["sp_name"][species]["asm_id"][assembly]["asm_file"]


# Report input collectors (GEP2-style: each returns [] if disabled) ------------

def get_report_gfastats_input(wildcards):
    """Gfastats output path for this assembly (always present)."""
    return os.path.join(
        config["OUT_FOLDER"], "GAME_results", wildcards.species,
        wildcards.assembly, "assembly_stats", "gfastats", "stats.txt"
    )


def get_report_compleasm_input(wildcards):
    """Compleasm tar.gz path (only if RUN_COMPL is enabled).
    The tar.gz contains {lineage}_odb{version}/full_table.tsv which is
    extracted at runtime for --compleasm-full."""
    if not _as_bool(config.get("RUN_COMPL", True)):
        return []
    return [os.path.join(
        config["OUT_FOLDER"], "GAME_results", wildcards.species,
        wildcards.assembly, "assembly_stats", "compleasm",
        "compleasm_results.tar.gz"
    )]


def get_all_report_inputs(wildcards):
    """Collect all inputs for the mini-report rule."""
    inputs = [get_report_gfastats_input(wildcards)]
    inputs.extend(get_report_compleasm_input(wildcards))
    return inputs


# -------------------------------------------------------------------------------
#  COMPLEASM EUKA DB AND PLACEMENT FILES (will be run only once!)
# ===============================================================================

rule A00_download_compleasm_db:
    output:
        flag="busco_lineages/eukaryota_odb12.done",
        placement_flag="busco_lineages/placement_files.done"
    params:
        db_dir="busco_lineages"
    container: CONTAINERS["compleasm"]
    benchmark:
        "logs/A00_download_compleasm_db_benchmark.txt"
    log:
        "logs/A00_download_db.log"
    threads: cpu_func("download_db")
    resources:
        mem_mb=mem_func("download_db"),
        runtime=time_func("download_db")
    shell:
        """
        mkdir -p {params.db_dir}
        compleasm download eukaryota --odb odb12 --library_path {params.db_dir}
        touch {output.flag}
        touch {output.placement_flag}
        """


# -------------------------------------------------------------------------------
#  ASM STATS
# ===============================================================================

rule A01_run_gfastats:
    input:
        assembly=lambda wildcards: get_assembly_path(wildcards.species, wildcards.assembly)
    output:
        stats=os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "assembly_stats", "gfastats", "stats.txt")
    container: CONTAINERS["game_base"]
    benchmark:
        os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "logs", "A01_gfastats_benchmark.txt")
    log:
        os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "logs", "A01_gfastats.log")
    threads: cpu_func("gfastats")
    resources:
        mem_mb=mem_func("gfastats"),
        runtime=time_func("gfastats")
    shell:
        """
        mkdir -p $(dirname {output.stats})
        mkdir -p $(dirname {log})
        gfastats -f {input.assembly} -j {threads} --nstar-report > {output.stats} 2> {log}
        """


rule A02_run_compleasm:
    input:
        assembly=lambda w: get_assembly_path(w.species, w.assembly),
        db_flag="busco_lineages/eukaryota_odb12.done",
        placement_flag="busco_lineages/placement_files.done"
    output:
        summary=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "assembly_stats", "compleasm", "summary.txt"
        ),
        targz=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "assembly_stats", "compleasm", "compleasm_results.tar.gz"
        )
    container: CONTAINERS["compleasm"]
    params:
        outdir=os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "assembly_stats", "compleasm"
        ),
        shared_db="busco_lineages",
        lineage=config.get("LINEAGE", "eukaryota")
    benchmark:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "logs", "A02_compleasm_benchmark.txt"
        )
    log:
        os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "logs", "A02_compleasm.log"
        )
    threads: cpu_func("compleasm")
    resources:
        mem_mb=mem_func("compleasm"),
        runtime=time_func("compleasm")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "{params.outdir}" "$(dirname {log})"

        exec > "{log}" 2>&1

        LINEAGE="{params.lineage}"
        echo "[compleasm] Configured LINEAGE: $LINEAGE"

        # TEMP DIRECTORY SETUP
        # -------------------------------------------------------------------
        WORK_DIR="$(game_get_workdir 50)"
        TEMP_DIR="$(mktemp -d "$WORK_DIR/GAME_compleasm_{wildcards.species}_{wildcards.assembly}_XXXXXX")"
        trap 'rm -rf "$TEMP_DIR"' EXIT

        TMPDB="$TEMP_DIR/compleasm_db"
        mkdir -p "$TMPDB"

        echo "[compleasm] Preparing isolated lineage DB in $TMPDB"
        # Using cp -rL instead of tar piping to avoid "short read" errors on some filesystems
        # This brings in eukaryota_odb12 + placement files (needed for --autolineage)
        cp -rL "{params.shared_db}"/* "$TMPDB/"

        # FETCH ADDITIONAL LINEAGE IF NEEDED (lives only in $TMPDB, dies with trap)
        # -------------------------------------------------------------------
        if [ "$LINEAGE" != "eukaryota" ] && [ "$LINEAGE" != "auto" ]; then
            echo "[compleasm] Downloading additional lineage '$LINEAGE' into temp DB"
            compleasm download "$LINEAGE" --odb odb12 --library_path "$TMPDB"
        fi

        # BUILD LINEAGE FLAG
        # -------------------------------------------------------------------
        if [ "$LINEAGE" = "auto" ]; then
            LINEAGE_FLAG="--autolineage"
            echo "[compleasm] Using --autolineage (placement files required)"
        else
            LINEAGE_FLAG="--lineage $LINEAGE"
            echo "[compleasm] Using --lineage $LINEAGE"
        fi

        # RUN COMPLEASM
        # -------------------------------------------------------------------
        compleasm run \
          -a "{input.assembly}" \
          -o "{params.outdir}" \
          -t {threads} \
          --library_path "$TMPDB" \
          $LINEAGE_FLAG

        # PACKAGE OUTPUTS
        # -------------------------------------------------------------------
        echo "[compleasm] Packaging outputs and cleaning directory"
        TB="{output.targz}"

        # Sync filesystem before packaging (avoids "short read" on network FS)
        sync

        # Package from local scratch to avoid parallel FS tar issues
        PACK_DIR="$TEMP_DIR/pack"
        mkdir -p "$PACK_DIR"
        cp -rL "{params.outdir}"/* "$PACK_DIR/" 2>/dev/null || true
        rm -f "$PACK_DIR"/compleasm_results.tar.gz

        cd "$PACK_DIR"
        tar -czf "$PACK_DIR/compleasm_results.tar.gz" --exclude='summary.txt' .

        # Verify archive integrity
        gzip -t "$PACK_DIR/compleasm_results.tar.gz"
        tar -tzf "$PACK_DIR/compleasm_results.tar.gz" >/dev/null

        # Move verified archive back to output directory
        mv "$PACK_DIR/compleasm_results.tar.gz" "$TB"
        echo "[compleasm] Archive verified and moved to $TB"

        # Handle summary.txt
        cd "{params.outdir}"
        if [ -f "summary.txt" ]; then
            echo "[compleasm] Found summary.txt in output directory"
        else
            echo "[compleasm] Warning: summary.txt not found in expected location, searching..."
            FOUND_SUMMARY=$(find . -name "summary.txt" -type f | head -1)
            if [ -n "$FOUND_SUMMARY" ]; then
                echo "[compleasm] Found summary at: $FOUND_SUMMARY"
                cp "$FOUND_SUMMARY" "{output.summary}"
            else
                echo "[compleasm] ERROR: No summary.txt found anywhere in output"
                echo "[compleasm] Contents of output directory:"
                ls -la
                exit 1
            fi
        fi

        # Clean up intermediate files, keeping only summary.txt and archive
        find . -mindepth 1 -maxdepth 1 ! -name 'summary.txt' ! -name '*_results.tar.gz' -exec rm -rf {{}} \;
        echo "[compleasm] Done: archive created at $TB"
        '''


rule A03_create_mini_report:
    input:
        deps = get_all_report_inputs
    output:
        report = os.path.join(
            config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}",
            "assembly_stats", "{assembly}_mini_report.md"
        )
    params:
        species   = lambda w: w.species,
        assembly  = lambda w: w.assembly,
        gfastats  = lambda w: get_report_gfastats_input(w),
        compleasm = lambda w: get_report_compleasm_input(w),
        lineage   = config.get("LINEAGE", "eukaryota"),
        script_path = str(scripts_dir / "make_gep2_report.py")
    container: CONTAINERS["game_base"]
    benchmark:
        os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "logs", "A03_mini_report_benchmark.txt")
    log:
        os.path.join(config["OUT_FOLDER"], "GAME_results", "{species}", "{assembly}", "logs", "A03_mini_report.log")
    threads: cpu_func("light_task")
    resources:
        mem_mb=mem_func("light_task"),
        runtime=time_func("light_task")
    shell:
        """
        mkdir -p $(dirname {log})
        exec > {log} 2>&1

        cmd="python {params.script_path} -s {params.species} -a {params.assembly} -g {params.gfastats}"

        # Compleasm: extract full_table.tsv from tar.gz files (if any)
        COMPLEASM_FULLS=""
        COMPLEASM_CLEANUP=""
        LINEAGE="{params.lineage}"

        for targz in {params.compleasm}; do
            if [ -n "$targz" ]; then
                extract_dir=$(dirname "$targz")
                tar -xzf "$targz" -C "$extract_dir"

                # Pick the right lineage's full_table.tsv:
                #   - specific lineage (e.g., aves): find that one
                #   - auto: prefer non-eukaryota (the auto-detected result)
                #   - eukaryota: find eukaryota
                if [ "$LINEAGE" = "auto" ]; then
                    # Prefer the non-eukaryota lineage (the interesting auto-detected one)
                    full_table=$(find "$extract_dir" -name "full_table.tsv" -path "*_odb*" \
                        ! -path "*eukaryota_odb*" | head -1)
                    # Fallback to eukaryota if no other lineage found
                    if [ -z "$full_table" ]; then
                        full_table=$(find "$extract_dir" -name "full_table.tsv" -path "*eukaryota_odb*" | head -1)
                    fi
                else
                    # Specific lineage: find exactly that one
                    full_table=$(find "$extract_dir" -name "full_table.tsv" -path "*${{LINEAGE}}_odb*" | head -1)
                fi

                if [ -n "$full_table" ]; then
                    COMPLEASM_FULLS="$COMPLEASM_FULLS $full_table"
                    echo "[GAME] Using compleasm full_table.tsv from: $(dirname "$full_table")"
                else
                    echo "[GAME] WARNING: No full_table.tsv found for lineage '$LINEAGE'"
                fi

                # Track ALL extracted odb directories for cleanup
                for odb_dir in $(find "$extract_dir" -maxdepth 1 -type d -name "*_odb*"); do
                    COMPLEASM_CLEANUP="$COMPLEASM_CLEANUP $odb_dir"
                done
            fi
        done

        if [ -n "$COMPLEASM_FULLS" ]; then
            cmd="$cmd --compleasm-full $COMPLEASM_FULLS"
        fi

        cmd="$cmd -o {output.report}"

        echo "[GAME] Command: $cmd"
        $cmd

        # Cleanup ALL extracted compleasm directories (eukaryota + target lineage)
        for cleanup_dir in $COMPLEASM_CLEANUP; do
            if [ -d "$cleanup_dir" ]; then
                rm -rf "$cleanup_dir"
                echo "[GAME] Cleaned up extracted compleasm: $cleanup_dir"
            fi
        done
        """