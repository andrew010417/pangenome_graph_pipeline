#!/bin/bash
# ============================================================
# 05_vg_giraffe_call.sh
# Purpose : Per-sample variant discovery by mapping short reads directly
#           onto the pangenome graph (as opposed to PanGenie's approach
#           of checking a known variant panel via k-mers -- see
#           07_pangenie_genotype.sh). Three-stage vg workflow:
#             1. vg giraffe  -- map paired-end short reads to the graph -> GAM
#             2. vg pack     -- collect per-node/edge read-support coverage
#             3. vg call     -- call variants from that coverage -> VCF
#
#           NOTE on `vg pack` GAM input flag: this pipeline's prior
#           hands-on experience (see project background) found GAM input
#           needs `-a`, NOT `-g`, despite `-g` looking like the more
#           obvious flag name for "GAM". This has been carried over as
#           given; flag names/meaning have shifted between vg releases,
#           so confirm against `vg pack --help` on your installed
#           version before trusting it (see VERIFICATION_TODO.md).
# Usage   : 05_vg_giraffe_call.sh -o <outdir> -x <index_prefix (from 04, without .giraffe.gbz)> -s <sample_name> -1 <R1.fastq.gz> -2 <R2.fastq.gz> [-t <threads=16>]
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

THREADS=16
OUTDIR=""
INDEX_PREFIX=""
SAMPLE=""
R1=""
R2=""

usage() {
    echo "Usage: $(basename "$0") -o <outdir> -x <index_prefix (from 04)> -s <sample_name> -1 <R1.fastq.gz> -2 <R2.fastq.gz> [-t <threads=16>]"
    exit 1
}

while getopts "o:x:s:1:2:t:h" opt; do
    case "${opt}" in
        o) OUTDIR="${OPTARG}" ;;
        x) INDEX_PREFIX="${OPTARG}" ;;
        s) SAMPLE="${OPTARG}" ;;
        1) R1="${OPTARG}" ;;
        2) R2="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -n "${OUTDIR}" && -n "${INDEX_PREFIX}" && -n "${SAMPLE}" && -n "${R1}" && -n "${R2}" ]] || usage

require_cmd vg
require_file "${INDEX_PREFIX}.giraffe.gbz" "Giraffe GBZ index"
require_file "${INDEX_PREFIX}.dist" "distance index"
require_file "${INDEX_PREFIX}.min" "minimizer index"
require_file "${R1}" "read 1 FASTQ"
require_file "${R2}" "read 2 FASTQ"

mkdir -p "${OUTDIR}"
GAM="${OUTDIR}/${SAMPLE}.gam"
PACK="${OUTDIR}/${SAMPLE}.pack"
VCF="${OUTDIR}/${SAMPLE}.vcf"

START_TIME=$(date +%s)
start_resource_monitor "${OUTDIR}" "vg"

log "[1/3] vg giraffe: mapping ${SAMPLE} reads to graph..."
vg giraffe \
    -Z "${INDEX_PREFIX}.giraffe.gbz" \
    -m "${INDEX_PREFIX}.min" \
    -d "${INDEX_PREFIX}.dist" \
    -f "${R1}" -f "${R2}" \
    -t "${THREADS}" \
    > "${GAM}" 2> "${OUTDIR}/${SAMPLE}.giraffe.log"

log "[2/3] vg pack: collecting read-support coverage..."
vg pack \
    -x "${INDEX_PREFIX}.giraffe.gbz" \
    -a "${GAM}" \
    -o "${PACK}" \
    -t "${THREADS}" \
    2>&1 | tee "${OUTDIR}/${SAMPLE}.pack.log"

log "[3/3] vg call: calling variants (ploidy 2)..."
vg call \
    "${INDEX_PREFIX}.giraffe.gbz" \
    -k "${PACK}" \
    -s "${SAMPLE}" \
    --ploidy 2 \
    -t "${THREADS}" \
    > "${VCF}" 2> "${OUTDIR}/${SAMPLE}.call.log"

bgzip -f -@ "${THREADS}" "${VCF}"
tabix -p vcf "${VCF}.gz" 2>/dev/null || log "(tabix not on PATH -- index the VCF manually)"

stop_resource_monitor
END_TIME=$(date +%s)
echo "Runtime: $((END_TIME - START_TIME)) sec" | tee "${OUTDIR}/${SAMPLE}.runtime.log"

log "Done. ${VCF}.gz"
