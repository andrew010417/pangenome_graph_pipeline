#!/bin/bash
# ============================================================
# 09_augment_graph.sh
# Purpose : Fold a sample's NOVEL variants (found by 05_vg_giraffe_call.sh
#           but absent from the graph's own variant panel) back into the
#           pangenome graph, so the graph keeps growing as new patients
#           are processed instead of staying frozen at whatever samples
#           built it originally.
#
#           This matters specifically because PanGenie (06/07) is a
#           genome-inference genotyper: BY DESIGN it can only report
#           variants that already exist somewhere in the graph/panel. A
#           variant it has never seen cannot be genotyped in the next
#           patient, no matter how real it is, until it's added here.
#
#           Workflow (per the augmentation approach described in the
#           project background reading -- see README):
#             1. bcftools isec -- keep only the sample's VCF (05) records
#                that are NOT already in the graph's panel VCF (02a/02b+03)
#             2. bcftools consensus -- apply those novel variants onto the
#                reference FASTA to materialize an actual haplotype
#                sequence carrying them (a graph can only be extended with
#                real sequence, not bare VCF records)
#             3. seqkit -- rename the new sequence to PanSN
#                (sample#hap#contig), consistent with 01_prepare_input_fasta.sh
#             4. minigraph -- re-align this new haplotype onto the
#                EXISTING graph and emit an augmented GFA, reusing shared
#                structure and adding nodes/edges only for what's new
#
#           NOTE: this produces ONE pseudo-haplotype per sample combining
#           all novel ALT alleles from an unphased vg-call VCF (hap index
#           "3", to avoid colliding with the sample's own 1/2 haplotypes
#           if it was itself an assembly input). If your sample VCF is
#           phased and you want two separate augmented haplotypes, run
#           `bcftools consensus -H 1` and `-H 2` separately instead of the
#           unphased default used below -- NOT done here, unverified
#           which is correct for vg call's output (see
#           VERIFICATION_TODO.md).
#
#           NOTE: `minigraph -cxggs` (construct+extend, graph-to-graph
#           augmentation flags) is carried over from Minigraph's
#           documented augmentation mode but has NOT been run end-to-end
#           on this pipeline. Confirm flags with `minigraph --help` /
#           `man minigraph` on your installed version before trusting
#           this (see VERIFICATION_TODO.md).
#
#           Augmenting after EVERY single sample is wasteful if many new
#           samples arrive close together -- consider accumulating novel
#           haplotype FASTAs (step 1-3 output) across several samples and
#           running the minigraph step (step 4) once in a batch instead.
# Usage   : 09_augment_graph.sh -o <outdir> -g <existing_graph.gfa[.gz]> -p <panel.vcf.gz> -v <sample.vcf.gz (from 05)> -r <reference.fa> -s <sample_name> [-t <threads=16>]
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

THREADS=16
OUTDIR=""
GRAPH=""
PANEL=""
SAMPLE_VCF=""
REF=""
SAMPLE=""

usage() {
    echo "Usage: $(basename "$0") -o <outdir> -g <existing_graph.gfa[.gz]> -p <panel.vcf.gz> -v <sample.vcf.gz (from 05)> -r <reference.fa> -s <sample_name> [-t <threads=16>]"
    exit 1
}

while getopts "o:g:p:v:r:s:t:h" opt; do
    case "${opt}" in
        o) OUTDIR="${OPTARG}" ;;
        g) GRAPH="${OPTARG}" ;;
        p) PANEL="${OPTARG}" ;;
        v) SAMPLE_VCF="${OPTARG}" ;;
        r) REF="${OPTARG}" ;;
        s) SAMPLE="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -n "${OUTDIR}" && -n "${GRAPH}" && -n "${PANEL}" && -n "${SAMPLE_VCF}" && -n "${REF}" && -n "${SAMPLE}" ]] || usage

require_cmd bcftools
require_cmd tabix
require_cmd seqkit
require_cmd minigraph
require_file "${GRAPH}" "existing graph GFA"
require_file "${PANEL}" "graph panel VCF"
require_file "${SAMPLE_VCF}" "sample VCF (from 05_vg_giraffe_call.sh)"
require_file "${REF}" "reference FASTA"

mkdir -p "${OUTDIR}/isec_tmp"
NOVEL_VCF="${OUTDIR}/${SAMPLE}.novel.vcf.gz"
NOVEL_FASTA="${OUTDIR}/${SAMPLE}.novel_hap.fa"
NOVEL_FASTA_RENAMED="${OUTDIR}/${SAMPLE}.novel_hap.pansn.fa"
AUGMENTED_GFA="${OUTDIR}/${SAMPLE}.augmented.gfa"

START_TIME=$(date +%s)
start_resource_monitor "${OUTDIR}" minigraph

log "[1/4] bcftools isec: finding variants in ${SAMPLE} absent from the graph panel..."
[[ -f "${SAMPLE_VCF}.tbi" ]] || tabix -p vcf "${SAMPLE_VCF}"
[[ -f "${PANEL}.tbi" ]] || tabix -p vcf "${PANEL}"
bcftools isec \
    -C \
    -w1 \
    -p "${OUTDIR}/isec_tmp" \
    -Oz \
    "${SAMPLE_VCF}" "${PANEL}" \
    2>&1 | tee "${OUTDIR}/${SAMPLE}.isec.log"
cp "${OUTDIR}/isec_tmp/0000.vcf.gz" "${NOVEL_VCF}"
tabix -p vcf "${NOVEL_VCF}"
log "  $(bcftools view -H "${NOVEL_VCF}" | wc -l) novel variant(s) not already in the panel"

log "[2/4] bcftools consensus: materializing novel alleles onto the reference..."
bcftools consensus \
    -f "${REF}" \
    "${NOVEL_VCF}" \
    > "${NOVEL_FASTA}" \
    2> "${OUTDIR}/${SAMPLE}.consensus.log"

log "[3/4] seqkit: renaming to PanSN (${SAMPLE}#3#contig)..."
seqkit replace -p '^(\S+).*' -r "${SAMPLE}#3#\${1}" "${NOVEL_FASTA}" > "${NOVEL_FASTA_RENAMED}"

log "[4/4] minigraph: augmenting existing graph with ${SAMPLE}'s novel haplotype..."
minigraph \
    -cxggs \
    -t "${THREADS}" \
    "${GRAPH}" \
    "${NOVEL_FASTA_RENAMED}" \
    > "${AUGMENTED_GFA}" \
    2> "${OUTDIR}/${SAMPLE}.minigraph.log"

stop_resource_monitor
END_TIME=$(date +%s)
echo "Runtime: $((END_TIME - START_TIME)) sec" | tee "${OUTDIR}/${SAMPLE}.runtime.log"

log "Done. Augmented graph: ${AUGMENTED_GFA}"
log "NEXT STEPS (not run automatically by this script):"
log "  - re-run 04_vg_autoindex_giraffe.sh -g ${AUGMENTED_GFA} to rebuild the Giraffe index"
log "  - re-extract a VCF panel from ${AUGMENTED_GFA} (vg deconstruct, as in 03) and re-run"
log "    06_prepare_pangenie_panel.sh so future PanGenie genotyping can see ${SAMPLE}'s novel alleles"
