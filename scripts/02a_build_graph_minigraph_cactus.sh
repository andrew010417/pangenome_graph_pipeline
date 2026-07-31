#!/bin/bash
# ============================================================
# 02a_build_graph_minigraph_cactus.sh
# Purpose : Build a pangenome graph with Minigraph-Cactus via the
#           `cactus-pangenome` single-command wrapper, producing a VCF
#           panel + a Giraffe-ready GBZ index directly (no separate
#           vg-deconstruct step needed, unlike the PGGB path — see
#           03_pggb_gfa_to_vcf.sh).
#
#           Benchmarked (see project README/PDF background) as
#           dramatically more practical than PGGB for population-scale
#           inputs: faster, far lower peak memory, and scales to large
#           haplotype counts where PGGB has been observed to fail/stall.
#           This is the recommended default graph for routine use;
#           PGGB (02b) is the secondary/comparison path.
#
#           NOTE: `cactus-pangenome` is the current documented entry
#           point for Minigraph-Cactus as of recent Cactus releases
#           (superseding the older manual
#           minigraph/cactus-graphmap/cactus-align/cactus-graphmap-join
#           step sequence). Flag availability differs across Cactus
#           versions -- run `cactus-pangenome --help` on your installed
#           version and diff against the flags below BEFORE a production
#           run; this has not been executed end-to-end on lab hardware
#           yet (see VERIFICATION_TODO.md).
#
#           Requires a Toil jobstore path (a scratch working directory,
#           NOT the same as --outDir). Reuse the same jobstore path with
#           `--restart` to resume a failed/interrupted run instead of
#           starting over.
# Usage   : 02a_build_graph_minigraph_cactus.sh -o <outdir> -j <jobstore_dir> -s <mc_seqfile.tsv> -r <reference_sample_id> [-t <threads=32>] [-n <outname=pangenome>]
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

THREADS=32
OUTDIR=""
JOBSTORE=""
SEQFILE=""
REF=""
OUTNAME="pangenome"

usage() {
    echo "Usage: $(basename "$0") -o <outdir> -j <jobstore_dir> -s <mc_seqfile.tsv> -r <reference_sample_id> [-t <threads=32>] [-n <outname=pangenome>]"
    exit 1
}

while getopts "o:j:s:r:t:n:h" opt; do
    case "${opt}" in
        o) OUTDIR="${OPTARG}" ;;
        j) JOBSTORE="${OPTARG}" ;;
        s) SEQFILE="${OPTARG}" ;;
        r) REF="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        n) OUTNAME="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -n "${OUTDIR}" && -n "${JOBSTORE}" && -n "${SEQFILE}" && -n "${REF}" ]] || usage

require_cmd cactus-pangenome
require_file "${SEQFILE}" "Minigraph-Cactus seqFile"
grep -q "^${REF}"$'\t' "${SEQFILE}" || die "reference sample id '${REF}' not found as a row in ${SEQFILE}"

mkdir -p "${OUTDIR}"

START_TIME=$(date +%s)
start_resource_monitor "${OUTDIR}" cactus-pangenome

log "[1/1] Running cactus-pangenome (reference=${REF}, outName=${OUTNAME})..."
cactus-pangenome \
    "${JOBSTORE}" \
    "${SEQFILE}" \
    --outDir "${OUTDIR}" \
    --outName "${OUTNAME}" \
    --reference "${REF}" \
    --vcf \
    --gbz \
    --gfa \
    --giraffe \
    --maxCores "${THREADS}" \
    2>&1 | tee "${OUTDIR}/cactus_pangenome.log"

log "Done. Expected outputs under ${OUTDIR}:"
log "  ${OUTNAME}.vcf.gz        (variant panel, directly usable by 06_prepare_pangenie_panel.sh)"
log "  ${OUTNAME}.gbz           (Giraffe-ready graph index, directly usable by 04_vg_autoindex_giraffe.sh)"
log "  ${OUTNAME}.gfa.gz        (graph, GFA format)"
log "NOTE: exact output filenames/suffixes are cactus-version-dependent and UNVERIFIED on this pipeline -- confirm against your run's actual output before wiring into 04/06."

stop_resource_monitor
END_TIME=$(date +%s)
echo "Runtime: $((END_TIME - START_TIME)) sec" | tee "${OUTDIR}/runtime.log"
