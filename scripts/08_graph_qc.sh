#!/bin/bash
# ============================================================
# 08_graph_qc.sh
# Purpose : Quality-check a pangenome graph BEFORE it's used for
#           genotyping (04/05/06/07), instead of trusting graph
#           construction (02a/02b) blindly. Runs two complementary
#           checks:
#             1. odgi stats   -- graph topology sanity check (node/edge
#                                 counts, connected components, node
#                                 length distribution). Flags structural
#                                 artefacts from bad construction
#                                 parameters (esp. relevant for PGGB --
#                                 see 02b's percent-identity/segment-length
#                                 flags).
#             2. panacus histgrowth -- pangenome growth / core-genome
#                                 curve: how much new sequence each
#                                 additional haplotype contributes, and
#                                 how much sequence is shared by all
#                                 (core) vs most (soft core) haplotypes.
#                                 A curve that never flattens suggests the
#                                 graph hasn't captured enough diversity
#                                 yet (too few/too similar input samples).
#
#           This does not replace VERIFICATION_TODO.md's manual checks
#           on cactus-pangenome/pggb output filenames -- run this AFTER
#           confirming 02a/02b actually produced the GFA you expect.
#
#           NOTE: odgi requires a blunt (non-overlapping-node) graph.
#           Some PGGB/Cactus GFAs may need `odgi build -g <gfa> -o
#           <og>` to succeed without extra flags; if it errors on your
#           graph, check `odgi build --help` for a `--sort`/
#           normalization flag on your installed odgi version -- this
#           has not been run end-to-end yet (see VERIFICATION_TODO.md).
# Usage   : 08_graph_qc.sh -o <outdir> -g <graph.gfa[.gz]> [-t <threads=8>]
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

THREADS=8
OUTDIR=""
GFA=""

usage() {
    echo "Usage: $(basename "$0") -o <outdir> -g <graph.gfa[.gz]> [-t <threads=8>]"
    exit 1
}

while getopts "o:g:t:h" opt; do
    case "${opt}" in
        o) OUTDIR="${OPTARG}" ;;
        g) GFA="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -n "${OUTDIR}" && -n "${GFA}" ]] || usage

require_cmd odgi
require_cmd panacus
require_file "${GFA}" "graph GFA"

mkdir -p "${OUTDIR}"
OG="${OUTDIR}/graph.og"

START_TIME=$(date +%s)
start_resource_monitor "${OUTDIR}" odgi

log "[1/3] odgi build: converting GFA -> odgi format..."
odgi build \
    -g "${GFA}" \
    -o "${OG}" \
    -t "${THREADS}" \
    2>&1 | tee "${OUTDIR}/odgi_build.log"

log "[2/3] odgi stats: topology summary (nodes, edges, length, components)..."
odgi stats \
    -i "${OG}" \
    -S \
    -t "${THREADS}" \
    > "${OUTDIR}/odgi_stats_summary.yaml" 2> "${OUTDIR}/odgi_stats.log"
cat "${OUTDIR}/odgi_stats_summary.yaml"

log "[3/3] panacus histgrowth: pangenome growth / core genome curve..."
panacus histgrowth \
    -t "${THREADS}" \
    "${GFA}" \
    > "${OUTDIR}/panacus_growth.tsv" 2> "${OUTDIR}/panacus.log"
log "  wrote ${OUTDIR}/panacus_growth.tsv -- open in R/pandas or run"
log "  'panacus-visualize ${OUTDIR}/panacus_growth.tsv' (if installed) for a growth-curve plot"
log "  matching Fig. 3d style (growth flattening = healthy; still climbing = add more samples)."

stop_resource_monitor
END_TIME=$(date +%s)
echo "Runtime: $((END_TIME - START_TIME)) sec" | tee "${OUTDIR}/runtime.log"

log "Done. Review ${OUTDIR}/odgi_stats_summary.yaml and ${OUTDIR}/panacus_growth.tsv"
log "BEFORE proceeding to 04_vg_autoindex_giraffe.sh / genotyping."
