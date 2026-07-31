#!/bin/bash
# ============================================================
# 06_prepare_pangenie_panel.sh
# Purpose : Normalize a graph-derived VCF panel (from 02a's
#           cactus-pangenome --vcf output, or 03's pggb.vcf.gz) into the
#           form PanGenie expects for genotyping:
#             - left-align + split multiallelic sites (`bcftools norm
#               -m -any`), since PanGenie's panel format is biallelic-site
#               based
#             - sort + index
#
#           On diploid genotype columns: if haplotypes were named
#           consistently as <sampleID>#1#<contig> / <sampleID>#2#<contig>
#           back in 01_prepare_input_fasta.sh, `vg deconstruct` (03) or
#           cactus-pangenome (02a) should already emit one diploid
#           genotype column per sampleID, and no further haplotype-merge
#           step should be needed here. This is DIFFERENT from a toy
#           panel built from single-haploid genomes (e.g. two reference
#           assemblies with no #1/#2 pairing), which need an explicit
#           haploid-pair-to-diploid merge that this script does NOT
#           implement. If `bcftools query -l` on your panel shows one
#           column per haplotype instead of per individual, treat that
#           as a sign this step needs the extra merge logic added before
#           trusting the genotyping output -- see VERIFICATION_TODO.md.
# Usage   : 06_prepare_pangenie_panel.sh -o <outdir> -v <panel.vcf.gz> -r <reference.fa> [-n <out_name=panel>]
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

OUTDIR=""
VCF=""
REF=""
OUTNAME="panel"

usage() {
    echo "Usage: $(basename "$0") -o <outdir> -v <panel.vcf.gz> -r <reference.fa> [-n <out_name=panel>]"
    exit 1
}

while getopts "o:v:r:n:h" opt; do
    case "${opt}" in
        o) OUTDIR="${OPTARG}" ;;
        v) VCF="${OPTARG}" ;;
        r) REF="${OPTARG}" ;;
        n) OUTNAME="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -n "${OUTDIR}" && -n "${VCF}" && -n "${REF}" ]] || usage
require_file "${VCF}" "panel VCF"
require_file "${REF}" "reference FASTA"

require_cmd bcftools

mkdir -p "${OUTDIR}"
OUT_VCF="${OUTDIR}/${OUTNAME}.biallelic.vcf.gz"

log "[1/3] bcftools norm: left-align + split multiallelic sites..."
bcftools norm -m -any -f "${REF}" -Oz -o "${OUT_VCF}" "${VCF}"

log "[2/3] Sorting + indexing..."
bcftools sort -Oz -o "${OUT_VCF}.sorted.gz" "${OUT_VCF}"
mv "${OUT_VCF}.sorted.gz" "${OUT_VCF}"
tabix -p vcf "${OUT_VCF}"

log "[3/3] Sample columns (verify per-individual, not per-haplotype -- see header note above):"
bcftools query -l "${OUT_VCF}" | tee "${OUTDIR}/${OUTNAME}_sample_columns.txt"

log "Done. ${OUT_VCF}"
