#!/usr/bin/env bash
#
# utils.sh
#
# Shared utility functions for the Illumina IDAT processing pipeline.
# Sourced by stage1, stage2, and other pipeline scripts.
#

# ---------------------------------------------------------------
# detect_sort_memory
#
# Dynamically allocate sort memory: use ~50% of available RAM
# (capped at 64G, floor at 4G). Works on Linux and macOS.
# ---------------------------------------------------------------
detect_sort_memory() {
    local sort_mem="4G"

    if [[ -f /proc/meminfo ]]; then
        # Linux: read available memory from /proc/meminfo
        sort_mem=$(awk '/MemAvailable/{
            mem_gb = int($2 / 1024 / 1024 * 0.5)
            if (mem_gb < 4) mem_gb = 4
            if (mem_gb > 64) mem_gb = 64
            printf "%dG", mem_gb
        }' /proc/meminfo)
    elif command -v sysctl &>/dev/null; then
        # macOS: read total physical memory via sysctl
        local mem_bytes
        mem_bytes=$(sysctl -n hw.memsize 2>/dev/null || true)
        if [[ -n "${mem_bytes}" && "${mem_bytes}" -gt 0 ]]; then
            sort_mem=$(awk -v bytes="${mem_bytes}" 'BEGIN {
                mem_gb = int(bytes / 1024 / 1024 / 1024 * 0.5)
                if (mem_gb < 4) mem_gb = 4
                if (mem_gb > 64) mem_gb = 64
                printf "%dG", mem_gb
            }')
        fi
    fi

    echo "${sort_mem}"
}

# ---------------------------------------------------------------
# get_par_xtr_bed
#
# Write a BED file of PAR1, XTR, and PAR2 regions on chrX and chrY
# for a given genome build.  Used to exclude pseudoautosomal and
# X-transposed regions from sex determination analyses.
#
# Usage: get_par_xtr_bed GENOME OUTPUT_BED
#   GENOME     - CHM13, GRCh38, or GRCh37
#   OUTPUT_BED - Path to write the BED file
# ---------------------------------------------------------------
get_par_xtr_bed() {
    local genome="${1:?Usage: get_par_xtr_bed GENOME OUTPUT_BED}"
    local output_bed="${2:?Usage: get_par_xtr_bed GENOME OUTPUT_BED}"

    case "${genome}" in
        CHM13)
            cat > "${output_bed}" <<'BED'
chrX	0	2781479	PAR1
chrX	2781479	6400875	XTR
chrX	155701382	156040895	PAR2
chrY	0	2458320	PAR1
chrY	2458320	6400875	XTR
chrY	62122809	62460029	PAR2
BED
            ;;
        GRCh38)
            cat > "${output_bed}" <<'BED'
chrX	10001	2781479	PAR1
chrX	2781479	6400000	XTR
chrX	155701383	156030895	PAR2
chrY	10001	2781479	PAR1
chrY	56887903	57217415	PAR2
BED
            ;;
        GRCh37)
            cat > "${output_bed}" <<'BED'
X	60001	2699520	PAR1
X	2699520	6100000	XTR
X	154931044	155260560	PAR2
Y	10001	2649520	PAR1
Y	59034050	59363566	PAR2
BED
            ;;
        *)
            echo "Warning: Unknown genome '${genome}' for PAR/XTR regions; using CHM13 defaults." >&2
            cat > "${output_bed}" <<'BED'
chrX	0	2781479	PAR1
chrX	2781479	6400875	XTR
chrX	155701382	156040895	PAR2
chrY	0	2458320	PAR1
chrY	2458320	6400875	XTR
chrY	62122809	62460029	PAR2
BED
            ;;
    esac
}

