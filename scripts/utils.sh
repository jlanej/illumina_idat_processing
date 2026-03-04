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

