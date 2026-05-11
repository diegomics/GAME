# ===============================================================================
#  GAME - Download Data Rules
# ===============================================================================

# -------------------------------------------------------------------------------
#  LOAD DOWNLOAD MANIFEST
# -------------------------------------------------------------------------------

manifest_file = os.path.join(config["OUT_FOLDER"], "GAME_results", "download_manifest.json")

if os.path.exists(manifest_file):
    with open(manifest_file) as f:
        DOWNLOAD_MANIFEST = json.load(f)
else:
    DOWNLOAD_MANIFEST = []


# -------------------------------------------------------------------------------
#  HELPER FUNCTIONS
# -------------------------------------------------------------------------------

def get_download_targets():
    """Get list of files that need to be downloaded"""
    if not os.path.exists(manifest_file):
        return []
    
    with open(manifest_file) as f:
        manifest = json.load(f)
    
    targets = []
    for item in manifest:
        dest = item["destination"]
        
        # For paired enaDataGet, we need both R1 and R2
        if item.get("method") == "enaDataGet" and item.get("paired"):
            targets.append(f"{dest}_1.fastq.gz")
            targets.append(f"{dest}_2.fastq.gz")
        else:
            targets.append(dest)
    
    return targets


def get_manifest_entry(destination):
    """Get download manifest entry for a given destination path"""
    for entry in DOWNLOAD_MANIFEST:
        if entry["destination"] == destination:
            return entry
    return None


# -------------------------------------------------------------------------------
#  WATCHDOG SCRIPT
# -------------------------------------------------------------------------------
# External helper script (workflow/scripts/watchdog.sh) provides:
#   - game_download_with_timeout: wraps a command with timeout(1)
#   - game_ena_download_paired/single: ENA portal API + HTTPS last-resort fallback
#   - game_http_fetch / game_http_fetch_stdout: curl/wget/python urllib chain

WATCHDOG_SCRIPT = str(scripts_dir / "watchdog.sh")


# -------------------------------------------------------------------------------
#  RULE ORDER
# -------------------------------------------------------------------------------
# Prefer paired-end download over single-end when both patterns match
# Prefer SRA downloads over URL downloads for accession-like filenames

ruleorder: _01_download_reads_sra > _01_download_reads_sra_single
ruleorder: _01_download_reads_sra_single > _01_download_reads_url
ruleorder: _01_download_reads_sra > _01_download_reads_url


# -------------------------------------------------------------------------------
#  WILDCARD CONSTRAINTS
# -------------------------------------------------------------------------------

wildcard_constraints:
    species=r"[^/]+",
    sample=r"[^/]+",
    read_type=r"[^/]+",
    accession=r"[SED]RR[0-9]+"


# -------------------------------------------------------------------------------
#  ASSEMBLY DOWNLOAD
# ===============================================================================

rule _00_download_assembly:
    """Download assemblies from URLs or NCBI accessions"""
    output:
        asm = "{outdir}/downloaded_data/assemblies/{species}/{filename}"
    wildcard_constraints:
        filename=r".+\.(fna|fa|fasta)(\.gz)?$"
    params:
        manifest = manifest_file
    threads: cpu_func("download_data")
    resources:
        mem_mb = mem_func("download_data"),
        runtime = time_func("download_data")
    container: CONTAINERS.get("game_base")
    shell:
        """
        #  GET MANIFEST INFO
        # -------------------------------------------------------------------
        MANIFEST_INFO=$(python3 -c "
import json, sys
with open('{params.manifest}') as f:
    manifest = json.load(f)
for item in manifest:
    if item.get('type') == 'assembly' and item['destination'] == '{output.asm}':
        print(item['source'])
        print(item['method'])
        sys.exit(0)
sys.exit(0)
")
        
        if [ -z "$MANIFEST_INFO" ]; then
            echo "[GAME] ❌ Error: No manifest entry found for {output.asm}"
            exit 1
        fi
        
        SOURCE=$(echo "$MANIFEST_INFO" | sed -n '1p')
        METHOD=$(echo "$MANIFEST_INFO" | sed -n '2p')
        
        mkdir -p $(dirname {output.asm})
        

        #  DOWNLOAD BY METHOD
        # -------------------------------------------------------------------
        if [ "$METHOD" = "curl" ]; then
            echo "[GAME] Downloading assembly from URL: $SOURCE"
            curl -L -C - --retry 3 --retry-delay 5 -o {output.asm}.tmp "$SOURCE"
            
            # Validate download
            if [ ! -s {output.asm}.tmp ]; then
                echo "[GAME] ❌ Error: Downloaded file is empty"
                exit 1
            fi
            
            # Check minimum file size (10KB for assemblies)
            FILE_SIZE=$(stat -c%s "{output.asm}.tmp" 2>/dev/null || echo "0")
            if [ "$FILE_SIZE" -lt 10240 ]; then
                echo "[GAME] ❌ Error: Downloaded file is suspiciously small ($FILE_SIZE bytes)"
                rm -f {output.asm}.tmp
                exit 1
            fi
            
            # Validate gzip integrity if compressed
            if [[ "{output.asm}" == *.gz ]]; then
                echo "[GAME] Validating gzip integrity..."
                if ! gzip -t {output.asm}.tmp 2>/dev/null; then
                    echo "[GAME] ❌ Error: Downloaded file failed gzip integrity check"
                    rm -f {output.asm}.tmp
                    exit 1
                fi
            fi
            
            mv {output.asm}.tmp {output.asm}
            echo "[GAME] Downloaded: {output.asm}"
            
        elif [ "$METHOD" = "ncbi_assembly" ]; then
            echo "[GAME] Downloading NCBI assembly: $SOURCE"
            
            # Parse accession using bash regex (e.g., GCA_963854735.1)
            if [[ $SOURCE =~ ^(GC[AF])_([0-9]{{3}})([0-9]{{3}})([0-9]{{3}})\\.([0-9]+)$ ]]; then
                PREFIX=${{BASH_REMATCH[1]}}
                P1=${{BASH_REMATCH[2]}}
                P2=${{BASH_REMATCH[3]}}
                P3=${{BASH_REMATCH[4]}}
                VERSION=${{BASH_REMATCH[5]}}
            else
                echo "[GAME] ❌ Error: Invalid NCBI accession format: $SOURCE"
                exit 1
            fi
            
            # Build base FTP directory URL
            BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/${{PREFIX}}/${{P1}}/${{P2}}/${{P3}}"
            echo "[GAME] Looking in: $BASE_URL"
            
            # Find the assembly directory (contains accession + assembly name)
            ASM_DIR=$(curl -sL "$BASE_URL/" | grep -oP "href=\\"${{SOURCE}}_[^/\\"]+" | head -1 | sed 's/href="//')
            
            if [ -z "$ASM_DIR" ]; then
                echo "[GAME] ❌ Error: Could not find assembly directory for $SOURCE"
                exit 1
            fi
            
            # Construct full URL to genomic.fna.gz
            FULL_URL="${{BASE_URL}}/${{ASM_DIR}}/${{ASM_DIR}}_genomic.fna.gz"
            echo "[GAME] Downloading from: $FULL_URL"
            
            curl -L -C - --retry 5 --retry-delay 10 -o {output.asm}.tmp "$FULL_URL"
            
            # Validate download
            if [ ! -s {output.asm}.tmp ]; then
                echo "[GAME] ❌ Error: Downloaded file is empty"
                exit 1
            fi
            
            # Check minimum file size
            FILE_SIZE=$(stat -c%s "{output.asm}.tmp" 2>/dev/null || echo "0")
            if [ "$FILE_SIZE" -lt 10240 ]; then
                echo "[GAME] ❌ Error: Downloaded file is suspiciously small ($FILE_SIZE bytes)"
                rm -f {output.asm}.tmp
                exit 1
            fi
            
            # Validate gzip file
            echo "[GAME] Validating gzip integrity..."
            if ! gzip -t {output.asm}.tmp 2>/dev/null; then
                echo "[GAME] ❌ Error: Downloaded file is not a valid gzip file"
                rm -f {output.asm}.tmp
                exit 1
            fi
            
            mv {output.asm}.tmp {output.asm}
            echo "[GAME] Downloaded NCBI assembly: {output.asm}"
        else
            echo "[GAME] ❌ Error: Unknown download method: $METHOD"
            exit 1
        fi
        """


# -------------------------------------------------------------------------------
#  READS DOWNLOAD - URL
# ===============================================================================

rule _01_download_reads_url:
    """Download reads from direct URLs (only matches non-SRA filenames)"""
    output:
        reads = "{outdir}/downloaded_data/reads/{sample}/{read_type}/{filename}"
    wildcard_constraints:
        filename=r"(?![SED]RR[0-9]+(_[12])?\.fastq(\.gz)?$).+\.(fastq|fq)(\.gz)?$"
    params:
        manifest = manifest_file
    threads: cpu_func("download_data")
    resources:
        mem_mb = mem_func("download_data"),
        runtime = time_func("download_data")
    container: CONTAINERS.get("game_base")
    shell:
        """
        SOURCE=$(python3 -c "
import json, sys
with open('{params.manifest}') as f:
    manifest = json.load(f)
for item in manifest:
    if item.get('type') == 'reads' and item.get('method') == 'curl' and item['destination'] == '{output.reads}':
        print(item['source'])
        sys.exit(0)
print('')
sys.exit(0)
")
        
        if [ -z "$SOURCE" ]; then
            echo "[GAME] ❌ Error: No URL source found in manifest for {output.reads}"
            exit 1
        fi
        
        mkdir -p $(dirname {output.reads})
        
        echo "[GAME] Downloading reads from URL: $SOURCE"
        curl -L -C - --retry 3 --retry-delay 5 -o {output.reads}.tmp "$SOURCE"
        
        # Validate download
        if [ ! -s {output.reads}.tmp ]; then
            echo "[GAME] ❌ Error: Downloaded file is empty"
            exit 1
        fi
        
        # Check minimum file size (1KB)
        FILE_SIZE=$(stat -c%s "{output.reads}.tmp" 2>/dev/null || echo "0")
        if [ "$FILE_SIZE" -lt 1024 ]; then
            echo "[GAME] ❌ Error: Downloaded file is suspiciously small ($FILE_SIZE bytes)"
            rm -f {output.reads}.tmp
            exit 1
        fi
        
        # Validate gzip integrity for compressed files
        if [[ "{output.reads}" == *.gz ]]; then
            echo "[GAME] Validating gzip integrity..."
            if ! gzip -t {output.reads}.tmp 2>/dev/null; then
                echo "[GAME] ❌ Error: Downloaded file failed gzip integrity check"
                rm -f {output.reads}.tmp
                exit 1
            fi
        fi
        
        mv {output.reads}.tmp {output.reads}
        echo "[GAME] Downloaded reads: {output.reads}"
        """


# -------------------------------------------------------------------------------
#  READS DOWNLOAD - SRA SINGLE-END / LONG READS
# ===============================================================================

rule _01_download_reads_sra_single:
    """Download single-end/long reads from SRA/ENA with aria2c fast path"""
    output:
        reads = "{outdir}/downloaded_data/reads/{sample}/{read_type}/{accession}.fastq.gz"
    log:
        "{outdir}/downloaded_data/reads/{sample}/{read_type}/{accession}.download.log"
    params:
        outdir = lambda w: os.path.join(w.outdir, "downloaded_data", "reads", w.sample, w.read_type),
        manifest = manifest_file,
        watchdog = WATCHDOG_SCRIPT
    threads: cpu_func("download_data")
    resources:
        mem_mb = mem_func("download_data"),
        runtime = time_func("download_data")
    container: CONTAINERS.get("game_base")
    shell:
        """
        exec > {log} 2>&1
        

        # VERIFY MANIFEST ENTRY
        # -------------------------------------------------------------------
        python3 -c "
import json, sys
with open('{params.manifest}') as f:
    manifest = json.load(f)
found = False
for item in manifest:
    if (item.get('type') == 'reads' and 
        item.get('method') == 'enaDataGet' and 
        item.get('source') == '{wildcards.accession}'):
        if not item.get('paired', False):
            found = True
            break
if not found:
    print('[GAME] ❌ Error: Accession {wildcards.accession} not found in manifest as single-end reads')
    sys.exit(1)
"
        
        source {params.watchdog}
        

        # CHECK IF DOWNLOAD PRODUCED FILES
        # -------------------------------------------------------------------
        check_single_files() {{
            local ACC=$1
            local DIR="."
            [ -d "$ACC" ] && DIR="$ACC"
            local SIZE
            
            # Compressed: must be >= 1 KB AND pass gzip integrity
            for f in "$DIR/${{ACC}}.fastq.gz" "$DIR/${{ACC}}_1.fastq.gz"; do
                if [ -f "$f" ]; then
                    SIZE=$(stat -c%s "$f" 2>/dev/null || echo "0")
                    if [ "$SIZE" -ge 1024 ] && gzip -t "$f" 2>/dev/null; then
                        return 0
                    fi
                fi
            done
            
            # Uncompressed: must be >= 1 KB (gzip check N/A)
            for f in "$DIR/${{ACC}}.fastq" "$DIR/${{ACC}}_1.fastq"; do
                if [ -f "$f" ]; then
                    SIZE=$(stat -c%s "$f" 2>/dev/null || echo "0")
                    if [ "$SIZE" -ge 1024 ]; then
                        return 0
                    fi
                fi
            done
            
            return 1
        }}
        

        # MAIN DOWNLOAD LOGIC
        # -------------------------------------------------------------------
        echo "[GAME] Downloading single-end/long reads: {wildcards.accession}"
        
        mkdir -p {params.outdir}
        cd {params.outdir}
        
        MAX_RETRIES=3
        RETRY_DELAY=60
        USE_ARIA2=true  # Start with aria2c, switch to enaDataGet HTTP if it fails
        
        for ATTEMPT in $(seq 1 $MAX_RETRIES); do
            echo ""
            echo "[GAME] ------------------------------------------------------------"
            echo "[GAME] Download attempt $ATTEMPT of $MAX_RETRIES"
            echo "[GAME] ============================================================"
            
            # Clean slate for this attempt
            rm -rf {wildcards.accession}/ {wildcards.accession}.fastq* {wildcards.accession}_*.fastq* 2>/dev/null || true


            # TRY ARIA2C FAST DOWNLOAD (parallel HTTPS via ENA portal API)
            # -----------------------------------------------------------------
            if [ "$USE_ARIA2" = "true" ]; then
                echo "[GAME] Trying aria2c (parallel HTTPS) download..."
                
                URLS=$(game_ena_get_urls {wildcards.accession} single)
                URLS_EXIT=$?
                
                if [ "$URLS_EXIT" -eq 0 ] && [ -n "$URLS" ]; then
                    URL=$(echo "$URLS" | sed -n '1p')
                    echo "[GAME]   URL: $URL"
                    
                    # NOTE: `|| ARIA_EXIT=$?` keeps the command off set -e's kill list.
                    # Failure here is expected (fallback drives the retry loop).
                    ARIA_EXIT=0
                    game_download_with_timeout 14400 aria2c \
                        -x 16 -s 16 -c \
                        --max-tries=3 --retry-wait=10 \
                        --allow-overwrite=true \
                        --auto-file-renaming=false \
                        --file-allocation=none \
                        --console-log-level=warn --summary-interval=30 \
                        -d . "$URL" || ARIA_EXIT=$?
                    
                    # Check if aria2c actually produced valid files (not just exit code!)
                    if check_single_files {wildcards.accession}; then
                        echo "[GAME] aria2c download produced files"
                    else
                        echo "[GAME] ⚠️  aria2c failed or produced no valid files (exit $ARIA_EXIT)"
                        echo "[GAME] Falling back to enaDataGet HTTP..."
                        USE_ARIA2=false
                        rm -rf {wildcards.accession}/ {wildcards.accession}.fastq* {wildcards.accession}_*.fastq* 2>/dev/null || true
                    fi
                else
                    echo "[GAME] ⚠️  Could not get URLs from ENA portal API (exit $URLS_EXIT)"
                    echo "[GAME] Falling back to enaDataGet HTTP..."
                    USE_ARIA2=false
                fi
            fi
            

            # TRY ENADATAGET HTTP (if aria2c failed or was skipped)
            # -----------------------------------------------------------------
            if [ "$USE_ARIA2" = "false" ]; then
                echo "[GAME] Using enaDataGet HTTP download..."
                
                HTTP_EXIT=0
                game_download_with_timeout 14400 enaDataGet.py -f fastq -d . {wildcards.accession} || HTTP_EXIT=$?
                
                echo "[GAME] enaDataGet HTTP finished (exit code: $HTTP_EXIT)"
            fi

            
            # TRY SUBMITTED FORMAT (if fastq format produced nothing)
            # -----------------------------------------------------------------
            if ! check_single_files {wildcards.accession}; then
                echo "[GAME] No FASTQ files found, trying submitted format..."
                rm -rf {wildcards.accession}/ {wildcards.accession}.fastq* {wildcards.accession}_*.fastq* 2>/dev/null || true
                
                SUBMITTED_EXIT=0
                game_download_with_timeout 14400 enaDataGet.py -f submitted -d . {wildcards.accession} || SUBMITTED_EXIT=$?
                
                echo "[GAME] Submitted files download finished (exit code: $SUBMITTED_EXIT)"
                
                # Rename unpredictable submitted file to standard naming
                DIR="."
                [ -d "{wildcards.accession}" ] && DIR="{wildcards.accession}"
                
                SUBMITTED=$(find "$DIR" -maxdepth 1 '(' -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" ')' 2>/dev/null | sort | head -1)
                
                if [ -n "$SUBMITTED" ]; then
                    if [[ "$SUBMITTED" == *.gz ]]; then
                        TARGET="$DIR/{wildcards.accession}.fastq.gz"
                    else
                        TARGET="$DIR/{wildcards.accession}.fastq"
                    fi
                    
                    echo "[GAME] Renaming submitted file to standard naming..."
                    echo "[GAME]   $SUBMITTED -> $TARGET"
                    mv "$SUBMITTED" "$TARGET"
                else
                    echo "[GAME] ⚠️  Warning: No submitted fastq files found"
                fi
            fi
            

            # TRY ENA PORTAL API + HTTPS (last-resort fallback)
            # -----------------------------------------------------------------
            # Handles submissions where enaDataGet.py returns "no files" or
            # crashes, e.g. submissions with FASTQ only in submitted_ftp, not fastq_ftp.
            if ! check_single_files {wildcards.accession}; then
                echo "[GAME] Trying ENA portal API + HTTPS as last resort..."
                rm -rf {wildcards.accession}/ {wildcards.accession}.fastq* {wildcards.accession}_*.fastq* 2>/dev/null || true
                
                API_EXIT=0
                game_ena_download_single {wildcards.accession} || API_EXIT=$?
                
                if [ "$API_EXIT" -eq 0 ] && check_single_files {wildcards.accession}; then
                    echo "[GAME] Portal-API download produced files"
                else
                    echo "[GAME] ⚠️  Portal-API fallback failed (exit $API_EXIT)"
                fi
            fi


            # PROCESS AND VALIDATE FILES
            # -----------------------------------------------------------------
            
            # Move files from subdirectory if created
            if [ -d "{wildcards.accession}" ]; then
                echo "[GAME] Moving files from subdirectory..."
                mv {wildcards.accession}/* . 2>/dev/null || true
                rmdir {wildcards.accession} 2>/dev/null || true
            fi
            
            # Compress any uncompressed fastq files
            for f in {wildcards.accession}*.fastq; do
                if [ -f "$f" ]; then
                    echo "[GAME] Compressing $f..."
                    pigz -p {threads} "$f"
                fi
            done
            
            # Rename variant files (_1, _subreads, etc.) to standard name
            for f in {wildcards.accession}_*.fastq.gz; do
                if [ -f "$f" ] && [ "$f" != "{wildcards.accession}.fastq.gz" ]; then
                    echo "[GAME] Renaming $f to {wildcards.accession}.fastq.gz"
                    mv "$f" "{wildcards.accession}.fastq.gz"
                    break
                fi
            done
            
            # Clean up any leftover paired-style file (shouldn't exist for single-end)
            rm -f "{wildcards.accession}_2.fastq.gz" "{wildcards.accession}_2.fastq" 2>/dev/null || true
            
            # Check if we got the file
            if [ -f "{output.reads}" ]; then
                echo "[GAME] Validating downloaded file..."
                
                FILE_SIZE=$(stat -c%s "{output.reads}" 2>/dev/null || echo "0")
                
                echo "[GAME] File size: $FILE_SIZE bytes"
                
                if [ "$FILE_SIZE" -lt 1024 ]; then
                    echo "[GAME] ⚠️  Downloaded file is suspiciously small"
                else
                    if gzip -t "{output.reads}" 2>/dev/null; then
                        echo "[GAME] Downloaded and validated: {wildcards.accession}"
                        exit 0
                    else
                        echo "[GAME] ⚠️  Downloaded file failed gzip integrity check"
                    fi
                fi
            else
                echo "[GAME] ⚠️  Expected file not found after download"
                echo "[GAME] Looking for: {output.reads}"
                echo "[GAME] Directory contents:"
                ls -la {params.outdir}/ 2>/dev/null || echo "(empty)"
            fi
            

            # RETRY LOGIC
            # -----------------------------------------------------------------
            if [ $ATTEMPT -lt $MAX_RETRIES ]; then
                echo ""
                echo "[GAME] ⚠️  Download attempt $ATTEMPT failed"
                echo "[GAME] Retrying in $RETRY_DELAY seconds..."
                sleep $RETRY_DELAY
                rm -rf {wildcards.accession}/ {wildcards.accession}.fastq* {wildcards.accession}_*.fastq* 2>/dev/null || true
            fi
        done
        
        # All retries exhausted
        echo ""
        echo "[GAME] ------------------------------------------------------------"
        echo "[GAME] ❌ ENA DOWNLOAD FAILED AFTER $MAX_RETRIES ATTEMPTS"
        echo "[GAME] ============================================================"
        echo "[GAME]"
        echo "[GAME] This could be due to:"
        echo "[GAME]   - ENA/EBI server issues"
        echo "[GAME]   - Network problems"
        echo "[GAME]   - Download repeatedly stalling (killed by watchdog)"
        echo "[GAME]"
        echo "[GAME] Please verify the accession at:"
        echo "[GAME]   https://www.ebi.ac.uk/ena/browser/view/{wildcards.accession}"
        echo "[GAME]"
        echo "[GAME] Contents of output directory:"
        ls -lh {params.outdir}/ 2>/dev/null || echo "[GAME] (directory empty or not found)"
        echo ""
        exit 1
        """


# -------------------------------------------------------------------------------
#  READS DOWNLOAD - SRA PAIRED-END
# ===============================================================================

rule _01_download_reads_sra:
    """Download paired-end reads from SRA/ENA with aria2c fast path"""
    output:
        r1 = "{outdir}/downloaded_data/reads/{sample}/{read_type}/{accession}_1.fastq.gz",
        r2 = "{outdir}/downloaded_data/reads/{sample}/{read_type}/{accession}_2.fastq.gz"
    log:
        "{outdir}/downloaded_data/reads/{sample}/{read_type}/{accession}.download.log"
    params:
        outdir = lambda w: os.path.join(w.outdir, "downloaded_data", "reads", w.sample, w.read_type),
        manifest = manifest_file,
        watchdog = WATCHDOG_SCRIPT
    threads: cpu_func("download_data")
    resources:
        mem_mb = mem_func("download_data"),
        runtime = time_func("download_data")
    container: CONTAINERS.get("game_base")
    shell:
        """
        exec > {log} 2>&1
        

        # VERIFY MANIFEST ENTRY
        # -------------------------------------------------------------------
        python3 -c "
import json, sys
with open('{params.manifest}') as f:
    manifest = json.load(f)
found = False
for item in manifest:
    if (item.get('type') == 'reads' and 
        item.get('method') == 'enaDataGet' and 
        item.get('paired') == True and
        item.get('source') == '{wildcards.accession}'):
        found = True
        break
if not found:
    print('[GAME] ❌ Error: Accession {wildcards.accession} not found in manifest as paired-end reads')
    sys.exit(1)
"
        
        source {params.watchdog}
        

        # CHECK IF DOWNLOAD PRODUCED FILES
        # -------------------------------------------------------------------
        check_paired_files() {{
            local ACC=$1
            local DIR="."
            [ -d "$ACC" ] && DIR="$ACC"
            local SIZE1 SIZE2
            
            # Both R1 and R2 must exist and be valid
            if [ -f "$DIR/${{ACC}}_1.fastq.gz" ] && [ -f "$DIR/${{ACC}}_2.fastq.gz" ]; then
                SIZE1=$(stat -c%s "$DIR/${{ACC}}_1.fastq.gz" 2>/dev/null || echo "0")
                SIZE2=$(stat -c%s "$DIR/${{ACC}}_2.fastq.gz" 2>/dev/null || echo "0")
                if [ "$SIZE1" -ge 1024 ] && [ "$SIZE2" -ge 1024 ] && \
                gzip -t "$DIR/${{ACC}}_1.fastq.gz" 2>/dev/null && \
                gzip -t "$DIR/${{ACC}}_2.fastq.gz" 2>/dev/null; then
                    return 0
                fi
            fi
            
            if [ -f "$DIR/${{ACC}}_1.fastq" ] && [ -f "$DIR/${{ACC}}_2.fastq" ]; then
                SIZE1=$(stat -c%s "$DIR/${{ACC}}_1.fastq" 2>/dev/null || echo "0")
                SIZE2=$(stat -c%s "$DIR/${{ACC}}_2.fastq" 2>/dev/null || echo "0")
                if [ "$SIZE1" -ge 1024 ] && [ "$SIZE2" -ge 1024 ]; then
                    return 0
                fi
            fi
            
            return 1
        }}


        # MAIN DOWNLOAD LOGIC
        # -------------------------------------------------------------------
        echo "[GAME] Downloading paired-end reads: {wildcards.accession}"
        
        mkdir -p {params.outdir}
        cd {params.outdir}
        
        MAX_RETRIES=3
        RETRY_DELAY=60
        USE_ARIA2=true  # Start with aria2c, switch to enaDataGet HTTP if it fails
        
        for ATTEMPT in $(seq 1 $MAX_RETRIES); do
            echo ""
            echo "[GAME] ------------------------------------------------------------"
            echo "[GAME] Download attempt $ATTEMPT of $MAX_RETRIES"
            echo "[GAME] ============================================================"
            
            # Clean slate for this attempt
            rm -rf {wildcards.accession}/ {wildcards.accession}_*.fastq* 2>/dev/null || true
            

            # TRY ARIA2C FAST DOWNLOAD (parallel HTTPS via ENA portal API)
            # -----------------------------------------------------------------
            if [ "$USE_ARIA2" = "true" ]; then
                echo "[GAME] Trying aria2c (parallel HTTPS) download..."
                
                URLS=$(game_ena_get_urls {wildcards.accession} paired)
                URLS_EXIT=$?
                
                if [ "$URLS_EXIT" -eq 0 ] && [ -n "$URLS" ]; then
                    URL1=$(echo "$URLS" | sed -n '1p')
                    URL2=$(echo "$URLS" | sed -n '2p')
                    echo "[GAME]   R1: $URL1"
                    echo "[GAME]   R2: $URL2"
                    
                    # NOTE: `|| ARIA_EXIT=$?` keeps the command off set -e's kill list.
                    # Failure here is expected (fallback drives the retry loop).
                    # Two sequential aria2c invocations (16 conns each) — simpler error
                    # tracking than one multi-file call.
                    ARIA_EXIT=0
                    for url in "$URL1" "$URL2"; do
                        game_download_with_timeout 14400 aria2c \
                            -x 16 -s 16 -c \
                            --max-tries=3 --retry-wait=10 \
                            --allow-overwrite=true \
                            --auto-file-renaming=false \
                            --file-allocation=none \
                            --console-log-level=warn --summary-interval=30 \
                            -d . "$url" || ARIA_EXIT=$?
                    done
                    
                    # Check if aria2c actually produced valid files (not just exit code!)
                    if check_paired_files {wildcards.accession}; then
                        echo "[GAME] aria2c download produced files"
                    else
                        echo "[GAME] ⚠️  aria2c failed or produced no valid files (exit $ARIA_EXIT)"
                        echo "[GAME] Falling back to enaDataGet HTTP..."
                        USE_ARIA2=false
                        rm -rf {wildcards.accession}/ {wildcards.accession}_*.fastq* 2>/dev/null || true
                    fi
                else
                    echo "[GAME] ⚠️  Could not get URLs from ENA portal API (exit $URLS_EXIT)"
                    echo "[GAME] Falling back to enaDataGet HTTP..."
                    USE_ARIA2=false
                fi
            fi
            

            # TRY ENADATAGET HTTP (if aria2c failed or was skipped)
            # -----------------------------------------------------------------
            if [ "$USE_ARIA2" = "false" ]; then
                echo "[GAME] Using enaDataGet HTTP download..."
                
                HTTP_EXIT=0
                game_download_with_timeout 14400 enaDataGet.py -f fastq -d . {wildcards.accession} || HTTP_EXIT=$?
                
                echo "[GAME] enaDataGet HTTP finished (exit code: $HTTP_EXIT)"
            fi
            

            # TRY SUBMITTED FORMAT (if fastq format produced nothing)
            # -----------------------------------------------------------------
            if ! check_paired_files {wildcards.accession}; then
                echo "[GAME] No FASTQ files found, trying submitted format..."
                rm -rf {wildcards.accession}/ {wildcards.accession}_*.fastq* 2>/dev/null || true
                
                SUBMITTED_EXIT=0
                game_download_with_timeout 14400 enaDataGet.py -f submitted -d . {wildcards.accession} || SUBMITTED_EXIT=$?
                
                echo "[GAME] Submitted files download finished (exit code: $SUBMITTED_EXIT)"

                # Rename unpredictable submitted files to standard naming
                DIR="."
                [ -d "{wildcards.accession}" ] && DIR="{wildcards.accession}"
                
                FASTQ_FILES=$(find "$DIR" -maxdepth 1 '(' -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" ')' 2>/dev/null | sort)
                NUM_FILES=$(echo "$FASTQ_FILES" | grep -c . || true)
                
                if [ "$NUM_FILES" -ge 2 ]; then
                    R1=$(echo "$FASTQ_FILES" | sed -n '1p')
                    R2=$(echo "$FASTQ_FILES" | sed -n '2p')
                    
                    # Detect extension to avoid renaming .fastq as .fastq.gz
                    if [[ "$R1" == *.gz ]]; then
                        TARGET_R1="$DIR/{wildcards.accession}_1.fastq.gz"
                        TARGET_R2="$DIR/{wildcards.accession}_2.fastq.gz"
                    else
                        TARGET_R1="$DIR/{wildcards.accession}_1.fastq"
                        TARGET_R2="$DIR/{wildcards.accession}_2.fastq"
                    fi
                    
                    echo "[GAME] Renaming submitted files to standard naming..."
                    echo "[GAME]   $R1 -> $TARGET_R1"
                    echo "[GAME]   $R2 -> $TARGET_R2"
                    mv "$R1" "$TARGET_R1"
                    mv "$R2" "$TARGET_R2"
                else
                    echo "[GAME] ⚠️  Warning: Could not find two submitted fastq files to pair"
                    echo "[GAME] Files found:"
                    echo "$FASTQ_FILES"
                fi
            fi
            

            # TRY ENA PORTAL API + HTTPS (last-resort fallback)
            # -----------------------------------------------------------------
            # Handles submissions where enaDataGet.py returns "no files" or
            # crashes, e.g. submissions with FASTQ only in submitted_ftp, not fastq_ftp.
            if ! check_paired_files {wildcards.accession}; then
                echo "[GAME] Trying ENA portal API + HTTPS as last resort..."
                rm -rf {wildcards.accession}/ {wildcards.accession}_*.fastq* 2>/dev/null || true
                
                API_EXIT=0
                game_ena_download_paired {wildcards.accession} || API_EXIT=$?
                
                if [ "$API_EXIT" -eq 0 ] && check_paired_files {wildcards.accession}; then
                    echo "[GAME] Portal-API download produced files"
                else
                    echo "[GAME] ⚠️  Portal-API fallback failed (exit $API_EXIT)"
                fi
            fi


            # PROCESS AND VALIDATE FILES
            # -----------------------------------------------------------------
            
            # Move files from subdirectory if created
            if [ -d "{wildcards.accession}" ]; then
                echo "[GAME] Moving files from subdirectory..."
                mv {wildcards.accession}/* . 2>/dev/null || true
                rmdir {wildcards.accession} 2>/dev/null || true
            fi
            
            # Compress if needed
            if [ -f "{wildcards.accession}_1.fastq" ]; then
                echo "[GAME] Compressing R1..."
                pigz -p {threads} "{wildcards.accession}_1.fastq"
            fi
            if [ -f "{wildcards.accession}_2.fastq" ]; then
                echo "[GAME] Compressing R2..."
                pigz -p {threads} "{wildcards.accession}_2.fastq"
            fi
            
            # Check if we got both files
            if [ -f "{output.r1}" ] && [ -f "{output.r2}" ]; then
                echo "[GAME] Validating downloaded files..."
                
                R1_SIZE=$(stat -c%s "{output.r1}" 2>/dev/null || echo "0")
                R2_SIZE=$(stat -c%s "{output.r2}" 2>/dev/null || echo "0")
                
                echo "[GAME] File sizes: R1=$R1_SIZE bytes, R2=$R2_SIZE bytes"
                
                if [ "$R1_SIZE" -lt 1024 ] || [ "$R2_SIZE" -lt 1024 ]; then
                    echo "[GAME] ⚠️  Downloaded files are suspiciously small"
                else
                    if gzip -t "{output.r1}" 2>/dev/null && gzip -t "{output.r2}" 2>/dev/null; then
                        echo "[GAME] Downloaded and validated paired reads: {wildcards.accession}"
                        exit 0
                    else
                        echo "[GAME] ⚠️  Downloaded files failed gzip integrity check"
                    fi
                fi
            else
                echo "[GAME] ⚠️  Expected files not found after download"
                echo "[GAME] Looking for: {output.r1}"
                echo "[GAME]             {output.r2}"
                echo "[GAME] Directory contents:"
                ls -la {params.outdir}/ 2>/dev/null || echo "(empty)"
            fi
            

            # RETRY LOGIC
            # -----------------------------------------------------------------
            if [ $ATTEMPT -lt $MAX_RETRIES ]; then
                echo ""
                echo "[GAME] ⚠️  Download attempt $ATTEMPT failed"
                echo "[GAME] Retrying in $RETRY_DELAY seconds..."
                sleep $RETRY_DELAY
                rm -rf {wildcards.accession}/ {wildcards.accession}_*.fastq* 2>/dev/null || true
            fi
        done
        
        # All retries exhausted
        echo ""
        echo "[GAME] ------------------------------------------------------------"
        echo "[GAME] ❌ ENA DOWNLOAD FAILED AFTER $MAX_RETRIES ATTEMPTS"
        echo "[GAME] ============================================================"
        echo "[GAME]"
        echo "[GAME] This could be due to:"
        echo "[GAME]   - ENA/EBI server issues"
        echo "[GAME]   - Network problems"
        echo "[GAME]   - Download repeatedly stalling (killed by watchdog)"
        echo "[GAME]"
        echo "[GAME] Please verify the accession at:"
        echo "[GAME]   https://www.ebi.ac.uk/ena/browser/view/{wildcards.accession}"
        echo "[GAME]"
        echo "[GAME] Contents of output directory:"
        ls -lh {params.outdir}/ 2>/dev/null || echo "[GAME] (directory empty or not found)"
        echo ""
        exit 1
        """
