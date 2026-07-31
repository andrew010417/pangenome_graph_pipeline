#!/bin/bash
# ============================================================
# 07_pangenie_genotype.sh
# Purpose : Genotype a sample against a known pangenome variant panel
#           (from 06) using PanGenie -- k-mer counting from the sample's
#           short reads, then HMM-based genotyping of each panel variant
#           as present/absent. This is the complementary path to
#           05_vg_giraffe_call.sh: vg discovers variants by mapping reads
#           onto the graph directly, PanGenie instead checks a
#           pre-defined variant list against the reads' k-mer content
#           (faster, but can only report what's already in the panel --
#           it does not discover novel variants).
#
#           NOTE: PanGenie CLI flags (-i/-r/-v/-o/-t/-j) match its
#           documented single-sample genotyping mode as of recent
#           releases; confirm against `PanGenie --help` on your
#           installed version before a production run (see
#           VERIFICATION_TODO.md) -- required/optional flags have
#           changed across PanGenie versions (e.g. the newer
#           index+genotype two-binary split in v3+).
# Usage   : 07_pangenie_genotype.sh -o <outdir> -s <sample_name> -q <reads.fastq.gz> -r <reference.fa> -v <panel.biallelic.vcf.gz> [-t <threads=16>]
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

THREADS=16
OUTDIR=""
SAMPLE=""
READS=""
REF=""
PANEL=""

usage() {
    echo "Usage: $(basename "$0") -o <outdir> -s <sample_name> -q <reads.fastq.gz> -r <reference.fa> -v <panel.biallelic.vcf.gz> [-t <threads=16>]"
    exit 1
}

while getopts "o:s:q:r:v:t:h" opt; do
    case "${opt}" in
        o) OUTDIR="${OPTARG}" ;;
        s) SAMPLE="${OPTARG}" ;;
        q) READS="${OPTARG}" ;;
        r) REF="${OPTARG}" ;;
        v) PANEL="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done

[[ -n "${OUTDIR}" && -n "${SAMPLE}" && -n "${READS}" && -n "${REF}" && -n "${PANEL}" ]] || usage

require_cmd PanGenie
require_file "${READS}" "reads FASTQ"
require_file "${REF}" "reference FASTA"
require_file "${PANEL}" "panel VCF"

mkdir -p "${OUTDIR}"
PREFIX="${OUTDIR}/${SAMPLE}"

START_TIME=$(date +%s)
start_resource_monitor "${OUTDIR}" PanGenie

log "[1/1] Running PanGenie (k-mer counting + HMM genotyping) for ${SAMPLE}..."
PanGenie \
    -i "${READS}" \
    -r "${REF}" \
    -v "${PANEL}" \
    -o "${PREFIX}" \
    -s "${SAMPLE}" \
    -t "${THREADS}" \
    -j "${THREADS}" \
    2>&1 | tee "${OUTDIR}/${SAMPLE}.pangenie.log"

stop_resource_monitor
END_TIME=$(date +%s)
echo "Runtime: $((END_TIME - START_TIME)) sec" | tee "${OUTDIR}/${SAMPLE}.runtime.log"

log "Done. Expected: ${PREFIX}_genotyping.vcf(.gz)"
