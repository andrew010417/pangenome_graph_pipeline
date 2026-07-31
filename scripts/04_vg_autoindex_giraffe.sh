#!/bin/bash
# ============================================================
# 04_vg_autoindex_giraffe.sh
# Purpose : Build the Giraffe-ready index set (.giraffe.gbz + .dist +
#           .min) that 05_vg_giraffe_call.sh needs to map reads.
#
#           Minigraph-Cactus (02a, run with --giraffe) already builds
#           this itself -- you normally only need THIS script for the
#           PGGB path (03's .pg/.vcf.gz output), or if you want to
#           rebuild/refresh an index. Two input modes:
#             -g <graph.gfa>                 build directly from a GFA
#             -r <ref.fa> -v <panel.vcf.gz>   build from ref + VCF panel
#           Give exactly one of these input modes.
#
#           NOTE: `vg autoindex --workflow giraffe` flag names/behavior
#           have shifted across vg releases; run `vg autoindex --help`
#           on your installed version and diff against the invocation
#           below before a production run (see VERIFICATION_TODO.md).
# Usage   : 04_vg_autoindex_giraffe.sh -o <outdir> -p <index_prefix> [-g <graph.gfa> | -r <ref.fa> -v <panel.vcf.gz>] [-t <threads=16>]
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

THREADS=16
OUTDIR=""
PREFIX=""
GFA=""
REF=""
VCF=""

usage() {
    echo "Usage: $(basename "$0") -o <outdir> -p <index_prefix> [-g <graph.gfa> | -r <ref.fa> -v <panel.vcf.gz>] [-t <threads=16>]"
    exit 1
}

while getopts "o:p:g:r:v:t:h" opt; do
    case "${opt}" in
        o) OUTDIR="${OPTARG}" ;;
        p) PREFIX="${OPTARG}" ;;
        g) GFA="${OPTARG}" ;;
        r) REF="${OPTARG}" ;;
        v) VCF="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -n "${OUTDIR}" && -n "${PREFIX}" ]] || usage
if [[ -n "${GFA}" && ( -n "${REF}" || -n "${VCF}" ) ]]; then
    die "give either -g, or -r + -v, not both"
fi
if [[ -z "${GFA}" && ( -z "${REF}" || -z "${VCF}" ) ]]; then
    die "give either -g <graph.gfa>, or both -r <ref.fa> and -v <panel.vcf.gz>"
fi

require_cmd vg
mkdir -p "${OUTDIR}"
OUT_PREFIX="${OUTDIR}/${PREFIX}"

START_TIME=$(date +%s)

if [[ -n "${GFA}" ]]; then
    require_file "${GFA}" "graph GFA"
    log "[1/1] vg autoindex --workflow giraffe (from GFA)..."
    vg autoindex --workflow giraffe -g "${GFA}" -p "${OUT_PREFIX}" -t "${THREADS}" \
        2>&1 | tee "${OUTDIR}/vg_autoindex.log"
else
    require_file "${REF}" "reference FASTA"
    require_file "${VCF}" "panel VCF"
    log "[1/1] vg autoindex --workflow giraffe (from ref + VCF panel)..."
    vg autoindex --workflow giraffe -r "${REF}" -v "${VCF}" -p "${OUT_PREFIX}" -t "${THREADS}" \
        2>&1 | tee "${OUTDIR}/vg_autoindex.log"
fi

END_TIME=$(date +%s)
echo "Runtime: $((END_TIME - START_TIME)) sec" | tee "${OUTDIR}/runtime.log"

log "Done. Expected: ${OUT_PREFIX}.giraffe.gbz, ${OUT_PREFIX}.dist, ${OUT_PREFIX}.min"
