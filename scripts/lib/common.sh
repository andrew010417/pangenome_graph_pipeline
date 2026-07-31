#!/bin/bash
# common.sh — shared helpers for the pangenome graph construction / genotyping pipeline.
# Source this from other scripts:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#   source "${SCRIPT_DIR}/lib/common.sh"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

die() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
    exit 1
}

require_file() {
    local path="$1" label="${2:-file}"
    [[ -f "${path}" ]] || die "${label} not found: ${path}"
}

require_cmd() {
    local cmd="$1"
    command -v "${cmd}" >/dev/null 2>&1 || die "required command not found on PATH: ${cmd}"
}

# Background resource monitor: system memory + RSS of a named process.
# Usage: start_resource_monitor <outdir> <process_name>
#        ... run the job ...
#        stop_resource_monitor
start_resource_monitor() {
    local outdir="$1" proc_name="$2"

    (
        while true; do
            echo "$(date +%s) $(free -g | awk '/Mem:/ {print $3}') $(top -bn1 | grep "${proc_name}" | awk '{print $9}' | head -n1)"
            sleep 30
        done
    ) > "${outdir}/resource.log" &
    _MONITOR_PID=$!

    (
        while true; do
            echo "$(date +%s) $(ps -C "${proc_name}" -o rss= | awk '{sum+=$1} END {print sum/1024/1024}')"
            sleep 30
        done
    ) > "${outdir}/${proc_name}_mem.log" &
    _PROC_MONITOR_PID=$!
}

stop_resource_monitor() {
    [[ -n "${_MONITOR_PID:-}" ]] && kill "${_MONITOR_PID}" 2>/dev/null
    [[ -n "${_PROC_MONITOR_PID:-}" ]] && kill "${_PROC_MONITOR_PID}" 2>/dev/null
}
