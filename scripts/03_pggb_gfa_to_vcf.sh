#!/bin/bash
# ============================================================
# 03_pggb_gfa_to_vcf.sh
# Purpose : PGGB-only step. Minigraph-Cactus emits a VCF panel directly
#           (see 02a), but PGGB only outputs a graph (GFA) -- this
#           converts it to a variant panel VCF via `vg convert` (GFA ->
#           .vg/.pg) then `vg deconstruct` (graph -> VCF against a
#           chosen reference path).
#
#           -P/--path-prefix selects which sample's paths in the graph
#           to call variants relative to (e.g. the reference genome's
#           PanSN sample id from 01_prepare_input_fasta.sh, such as
#           "grch38"). All other samples' paths become the VCF's ALT
#           alleles/genotype columns.
#
#           If haplotypes were named consistently as
#           <sampleID>#1#<contig> / <sampleID>#2#<contig> in step 01,
#           vg deconstruct should group them into one diploid genotype
#           column per sampleID automatically -- this is a property of
#           vg's PanSN-aware sample grouping, not something this script
#           does explicitly. UNVERIFIED on this pipeline; check
#           `bcftools query -l` on the output VCF and confirm the sample
#           columns match expectations before trusting downstream
#           PanGenie/vg-call results.
# Usage   : 03_pggb_gfa_to_vcf.sh -o <outdir> -g <smooth.final.gfa> -r <reference_sample_id> [-t <threads=16>]
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

THREADS=16
OUTDIR=""
GFA=""
REF=""

usage() {
    echo "Usage: $(basename "$0") -o <outdir> -g <smooth.final.gfa> -r <reference_sample_id> [-t <threads=16>]"
    exit 1
}

while getopts "o:g:r:t:h" opt; do
    case "${opt}" in
        o) OUTDIR="${OPTARG}" ;;
        g) GFA="${OPTARG}" ;;
        r) REF="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -n "${OUTDIR}" && -n "${GFA}" && -n "${REF}" ]] || usage
require_file "${GFA}" "PGGB GFA"
require_cmd vg
require_cmd bgzip
require_cmd tabix

mkdir -p "${OUTDIR}"
PG="${OUTDIR}/pggb_graph.pg"
VCF="${OUTDIR}/pggb.vcf"

START_TIME=$(date +%s)

log "[1/3] vg convert: GFA -> .pg ..."
vg convert -t "${THREADS}" -g "${GFA}" -p > "${PG}"

log "[2/3] vg deconstruct: .pg -> VCF (reference sample: ${REF}) ..."
vg deconstruct -t "${THREADS}" -P "${REF}" -e "${PG}" > "${VCF}"

log "[3/3] bgzip + tabix ..."
bgzip -f -@ "${THREADS}" "${VCF}"
tabix -p vcf "${VCF}.gz"

END_TIME=$(date +%s)
echo "Runtime: $((END_TIME - START_TIME)) sec" | tee "${OUTDIR}/runtime.log"

log "Done. Sample columns in ${VCF}.gz:"
bcftools query -l "${VCF}.gz" 2>/dev/null | tee "${OUTDIR}/sample_columns.txt" || log "(bcftools not on PATH -- skipped sample column listing, check manually)"
