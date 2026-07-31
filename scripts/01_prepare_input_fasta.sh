#!/bin/bash
# ============================================================
# 01_prepare_input_fasta.sh
# Purpose : Normalize per-sample assembly FASTAs (output of the
#           genome_assembly_ONT(hi-c) pipeline, or any diploid/haploid
#           assembly + reference genomes) into the two input shapes the
#           downstream graph builders need:
#             - PGGB wants ALL sequences in ONE bgzipped multi-FASTA,
#               with headers following the PanSN-spec naming convention
#               (sample#hap#contig), which is how it tells samples/
#               haplotypes apart when building paths in the graph.
#             - Minigraph-Cactus (cactus-pangenome) wants a tab-separated
#               seqFile: one line per SAMPLE with a path to that sample's
#               FASTA (haplotypes of the same sample go in one FASTA, or
#               are given as separate rows with hap-suffixed sample names
#               depending on cactus version — see README).
#
#           Naming samples/haplotypes consistently here (sampleID#1,
#           sampleID#2 for a phased diploid; sampleID#0 for an unphased
#           single assembly) matters beyond cosmetics: `vg deconstruct`
#           (used in 03_pggb_gfa_to_vcf.sh) groups paths into genotype
#           columns by the sample# prefix, so a diploid sample assembled
#           with hap1/hap2 under the SAME sampleID#1 / sampleID#2 prefix
#           should come out of vg deconstruct already genotyped as one
#           diploid VCF column per sample -- no post-hoc haplotype-merge
#           workaround should be needed for that case. This has NOT been
#           empirically confirmed on this pipeline yet; verify the sample
#           columns with `bcftools query -l` after 03/06 and compare
#           against what you expected before trusting it in production.
#
# Input   : a manifest TSV: <sample_id>\t<hap_index>\t<fasta_path>
#           hap_index: 0 for an unphased/haploid assembly (e.g.
#           *_primary.fa from 05_hifiasm_ont_only.sh), 1/2 for a phased
#           diploid assembly (e.g. *_hap1.fa/*_hap2.fa from
#           06_hifiasm_ont_hic.sh). Reference genomes (GRCh38, CHM13,
#           etc.) go in the same manifest with hap_index 0.
# Usage   : 01_prepare_input_fasta.sh -o <outdir> -m <manifest.tsv> [-t <threads=8>]
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

THREADS=8
OUTDIR=""
MANIFEST=""

usage() {
    echo "Usage: $(basename "$0") -o <outdir> -m <manifest.tsv> [-t <threads=8>]"
    echo "  manifest.tsv columns: sample_id<TAB>hap_index<TAB>fasta_path"
    exit 1
}

while getopts "o:m:t:h" opt; do
    case "${opt}" in
        o) OUTDIR="${OPTARG}" ;;
        m) MANIFEST="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -n "${OUTDIR}" && -n "${MANIFEST}" ]] || usage
require_file "${MANIFEST}" "manifest"

require_cmd seqkit
require_cmd bgzip
require_cmd samtools

mkdir -p "${OUTDIR}/renamed" "${OUTDIR}/logs"
MC_SEQFILE="${OUTDIR}/mc_seqfile.tsv"
PGGB_MERGED="${OUTDIR}/pggb_input.fa"
: > "${MC_SEQFILE}"
: > "${PGGB_MERGED}"

log "[1/3] Renaming headers to PanSN-spec (sample#hap#contig) per manifest line..."
while IFS=$'\t' read -r SAMPLE HAP FASTA; do
    [[ -z "${SAMPLE}" || "${SAMPLE}" =~ ^# ]] && continue
    require_file "${FASTA}" "assembly FASTA for ${SAMPLE}#${HAP}"

    RENAMED="${OUTDIR}/renamed/${SAMPLE}.hap${HAP}.fa"
    # sed-based rename: ">chr1 desc" -> ">SAMPLE#HAP#chr1"
    seqkit replace -p '^(\S+).*' -r "${SAMPLE}#${HAP}#\${1}" "${FASTA}" > "${RENAMED}"

    echo -e "${SAMPLE}\t${RENAMED}" >> "${MC_SEQFILE}"
    cat "${RENAMED}" >> "${PGGB_MERGED}"

    log "  ${SAMPLE}#${HAP} <- ${FASTA}"
done < "${MANIFEST}"

log "[2/3] bgzip + faidx the merged PGGB input..."
bgzip -f -@ "${THREADS}" "${PGGB_MERGED}"
samtools faidx "${PGGB_MERGED}.gz"

log "[3/3] Sanity check: sample/haplotype prefixes present in merged FASTA..."
seqkit seq -n "${PGGB_MERGED}.gz" | sed -E 's/#.*//' | sort -u | tee "${OUTDIR}/logs/sample_ids.txt"

log "Done."
log "  Minigraph-Cactus seqFile : ${MC_SEQFILE}"
log "  PGGB merged input        : ${PGGB_MERGED}.gz"
