#!/bin/bash
# ============================================================
# 02b_build_graph_pggb.sh
# Purpose : Build a pangenome graph with PGGB (reference-free, all-vs-all
#           alignment). Produces a richer/more finely-branched graph than
#           Minigraph-Cactus (edges >> nodes) at substantially higher
#           runtime and memory cost, and has been observed (see project
#           background) to fail outright on large haplotype counts
#           (~100+). Use this as the secondary/comparison path for
#           smaller sample sets or specific loci where PGGB's more
#           complete variation representation is worth the cost; default
#           to 02a_build_graph_minigraph_cactus.sh for routine/larger runs.
#
#           Unlike Minigraph-Cactus, PGGB does NOT emit a VCF directly --
#           run 03_pggb_gfa_to_vcf.sh afterward on its output GFA.
#
#           NOTE: the exact output filename PGGB writes for the final
#           smoothed graph (referred to here as smooth.final.gfa) varies
#           with PGGB version and the -p/-s parameters used, since PGGB
#           embeds those params in its output directory/file names. This
#           script locates it with `find` as a best effort; confirm the
#           match is correct after the first run on your PGGB version
#           (see VERIFICATION_TODO.md).
# Usage   : 02b_build_graph_pggb.sh -o <outdir> -i <pggb_input.fa.gz> -n <num_haplotypes> [-t <threads=32>] [-p <map_pct_id=90>] [-s <segment_length=10000>]
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

THREADS=32
OUTDIR=""
INPUT_FA=""
NHAP=""
PCT_ID=90
SEG_LEN=10000

usage() {
    echo "Usage: $(basename "$0") -o <outdir> -i <pggb_input.fa.gz> -n <num_haplotypes> [-t <threads=32>] [-p <map_pct_id=90>] [-s <segment_length=10000>]"
    exit 1
}

while getopts "o:i:n:t:p:s:h" opt; do
    case "${opt}" in
        o) OUTDIR="${OPTARG}" ;;
        i) INPUT_FA="${OPTARG}" ;;
        n) NHAP="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        p) PCT_ID="${OPTARG}" ;;
        s) SEG_LEN="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -n "${OUTDIR}" && -n "${INPUT_FA}" && -n "${NHAP}" ]] || usage
require_file "${INPUT_FA}" "PGGB input FASTA"
require_file "${INPUT_FA}.fai" "PGGB input FASTA index (run samtools faidx / step 01 first)"

require_cmd pggb

mkdir -p "${OUTDIR}"

START_TIME=$(date +%s)
start_resource_monitor "${OUTDIR}" pggb

log "[1/2] Running pggb (n=${NHAP} haplotypes, -p ${PCT_ID}, -s ${SEG_LEN})..."
pggb \
    -i "${INPUT_FA}" \
    -o "${OUTDIR}" \
    -n "${NHAP}" \
    -t "${THREADS}" \
    -p "${PCT_ID}" \
    -s "${SEG_LEN}" \
    2>&1 | tee "${OUTDIR}/pggb.log"

log "[2/2] Locating final smoothed GFA..."
FINAL_GFA="$(find "${OUTDIR}" -maxdepth 1 -name '*smooth.final.gfa' | head -n1)"
if [[ -z "${FINAL_GFA}" ]]; then
    log "WARNING: no *smooth.final.gfa found under ${OUTDIR} -- inspect ${OUTDIR}/pggb.log and the output directory manually, PGGB's naming convention may differ on your version."
else
    log "Final GFA: ${FINAL_GFA}"
    ln -sf "$(basename "${FINAL_GFA}")" "${OUTDIR}/smooth.final.gfa"
fi

stop_resource_monitor
END_TIME=$(date +%s)
echo "Runtime: $((END_TIME - START_TIME)) sec" | tee "${OUTDIR}/runtime.log"

log "Done. Next: 03_pggb_gfa_to_vcf.sh -g ${OUTDIR}/smooth.final.gfa ..."
