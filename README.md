# Pangenome Graph Construction & Genotyping Pipeline

문의사항은 Jaehyung Park(JP)에게 문의해주세요.

Builds a pangenome graph from multiple per-sample genome assemblies
(e.g. the output of [`genome_assembly_ONT(hi-c)`](../genome_assembly_ONT(hi-c)))
plus reference genome(s), using **Minigraph-Cactus** and/or **PGGB**,
then genotypes/calls variants for individual short-read samples against
that graph using **vg** (direct read-to-graph mapping) and/or
**PanGenie** (k-mer-based panel genotyping).

Intended as the downstream step after per-sample long-read assembly:
`genome_assembly_ONT(hi-c)` produces individual assemblies; this repo
combines them into a graph and uses it to find variants (including rare
ones a single linear reference would miss) in short-read samples.

## Provenance

Adapted from pangenome graph construction work done previously at KOGIC
(국바빅 범유전체 project) comparing Minigraph-Cactus and PGGB on
CHM13/GRCh38, and their downstream use with `vg` and `PanGenie` for
variant discovery/genotyping on GIAB HG002 data. That project's own
benchmarking (see internal planning notes) found Minigraph-Cactus far
more practical for larger sample counts — faster, much lower peak
memory, scales further — while PGGB produces a more finely-branched
graph (more edges, more bubble structure) from its reference-free
all-vs-all alignment, at the cost of runtime/memory and reliability at
scale (observed failures at ~100+ haplotypes in that prior benchmark).

Unlike `genome_assembly_ONT(hi-c)`, the original scripts from that
project were not available to adapt directly — only summary notes and
benchmark numbers. **The scripts here are freshly written against each
tool's documented CLI**, not copies of a previously-running pipeline.
Treat filenames, exact flags, and version-specific behavior as
unverified until run end-to-end on a small test dataset — see
[VERIFICATION_TODO.md](VERIFICATION_TODO.md).

## Requirements

- `cactus` (>= 9.x with `cactus-pangenome`) — **not** conda-installable,
  see `environment.yml` for install notes
- `pggb`, `vg`, `PanGenie`, `bcftools`, `samtools`, `htslib` (tabix/bgzip), `seqkit`

A version-pinned conda environment for everything except Cactus is
provided in `environment.yml`:

```bash
mamba env create -f environment.yml   # or: conda env create -f environment.yml
conda activate pangenome-graph
```

Scripts are scheduler-agnostic plain bash (no SGE/SLURM headers). Wrap
with `qsub`/`sbatch` yourself if your cluster needs it. Thread counts
(`-t`) default to a conservative value in each script; adjust to your
server's actual core count.

실제 랩 서버에서 처음 실행하기 전에 [VERIFICATION_TODO.md](VERIFICATION_TODO.md)의
검증 항목을 먼저 확인해주세요.

## Pipeline

```
scripts/
├── lib/common.sh                       shared helpers (logging, resource monitor)
├── 01_prepare_input_fasta.sh           PanSN-spec header renaming + MC seqFile / PGGB merged FASTA
├── 02a_build_graph_minigraph_cactus.sh cactus-pangenome  -> VCF panel + Giraffe GBZ (recommended default)
├── 02b_build_graph_pggb.sh             pggb               -> GFA graph (secondary/comparison path)
├── 03_pggb_gfa_to_vcf.sh               vg convert + vg deconstruct  -> VCF panel (PGGB path only)
├── 04_vg_autoindex_giraffe.sh          vg autoindex --workflow giraffe -> GBZ/dist/min (PGGB path, or to rebuild MC's)
├── 05_vg_giraffe_call.sh               vg giraffe -> vg pack -> vg call   (per-sample, direct graph mapping)
├── 06_prepare_pangenie_panel.sh        bcftools norm -m -any  -> biallelic panel VCF
└── 07_pangenie_genotype.sh             PanGenie k-mer counting + HMM genotyping  (per-sample, panel-based)
```

**Which graph builder**: use Minigraph-Cactus (`02a`) by default —
faster, lower memory, handles more samples. Use PGGB (`02b` + `03`) for
smaller sample sets or loci where its more complete/reference-free
variation representation is worth the extra cost, or when you want to
compare both against each other.

**Which variant-calling path**: `vg` (`05`) discovers variants directly
from how reads map onto the graph — it can find things not already in
the graph. `PanGenie` (`06` + `07`) instead checks a fixed panel of
already-known graph variants against a sample's k-mer content — faster,
but limited to what's in the panel. Running both on the same sample is
a reasonable way to cross-check calls.

### Minigraph-Cactus path (recommended default)

```bash
scripts/01_prepare_input_fasta.sh -o prep -m samples_manifest.tsv -t 8

scripts/02a_build_graph_minigraph_cactus.sh \
    -o graph_mc -j /scratch/jobstore_mc \
    -s prep/mc_seqfile.tsv -r grch38 -t 32 -n pangenome

scripts/05_vg_giraffe_call.sh \
    -o calls_vg/sample01 -x graph_mc/pangenome -s sample01 \
    -1 sample01_R1.fastq.gz -2 sample01_R2.fastq.gz -t 16

scripts/06_prepare_pangenie_panel.sh \
    -o panel -v graph_mc/pangenome.vcf.gz -r ref/grch38.fa

scripts/07_pangenie_genotype.sh \
    -o calls_pangenie/sample01 -s sample01 \
    -q sample01_reads.fastq.gz -r ref/grch38.fa -v panel/panel.biallelic.vcf.gz -t 16
```

### PGGB path (secondary / comparison)

```bash
scripts/01_prepare_input_fasta.sh -o prep -m samples_manifest.tsv -t 8

scripts/02b_build_graph_pggb.sh -o graph_pggb -i prep/pggb_input.fa.gz -n <num_haplotypes> -t 32

scripts/03_pggb_gfa_to_vcf.sh -o graph_pggb -g graph_pggb/smooth.final.gfa -r grch38 -t 16

scripts/04_vg_autoindex_giraffe.sh \
    -o graph_pggb -p pangenome_pggb -r ref/grch38.fa -v graph_pggb/pggb.vcf.gz -t 16

scripts/05_vg_giraffe_call.sh \
    -o calls_vg_pggb/sample01 -x graph_pggb/pangenome_pggb -s sample01 \
    -1 sample01_R1.fastq.gz -2 sample01_R2.fastq.gz -t 16
```

`samples_manifest.tsv` format (see `01_prepare_input_fasta.sh` header
comment): `sample_id<TAB>hap_index<TAB>fasta_path`, one row per
haplotype/reference genome, e.g.:

```
sample01	1	assemblies/sample01_hap1.fa
sample01	2	assemblies/sample01_hap2.fa
grch38	0	ref/grch38.fa
```

## Notes / caveats

- **Cactus is not conda-installed here** — see `environment.yml` for
  why and where to get it.
- **PanSN naming matters for diploid genotyping**: `01_prepare_input_fasta.sh`
  renames headers to `sample#hap#contig`. Pairing hap 1/2 of the same
  individual under one `sample_id` is what lets `vg deconstruct` (03)
  and `cactus-pangenome` (02a) emit one diploid genotype column per
  individual automatically — this has not been empirically confirmed on
  this pipeline yet (see VERIFICATION_TODO.md). If your panel VCF ends
  up with one column per haplotype instead of per individual,
  `06_prepare_pangenie_panel.sh` will need an added haplotype-merge step
  before PanGenie genotyping.
- **`vg pack` GAM input flag**: `05_vg_giraffe_call.sh` uses `-a` for
  GAM input (per prior hands-on experience with this exact gotcha, not
  the more obviously-named `-g`). Re-verify against `vg pack --help` on
  your installed vg version.
- **Output filenames are version-dependent and unverified** for
  `cactus-pangenome` (02a) and `pggb`'s final smoothed GFA (02b) — both
  scripts log what's expected but haven't been run end-to-end here.
- **Resource requirements**: prior benchmarking (2 haplotypes, CHM13 +
  GRCh38) saw Minigraph-Cactus finish in ~7.5h at ~108 GB peak memory
  (32 threads) vs. PGGB's ~21h at lower peak memory but far worse
  scaling with more haplotypes/threads. Expect this to scale up
  significantly with real per-sample diploid assemblies and more
  samples — size your job accordingly before submitting to the cluster.
- **This repo does not cover per-sample long-read assembly** — that's
  [`genome_assembly_ONT(hi-c)`](../genome_assembly_ONT(hi-c)), whose
  `_primary.fa`/`_alternate.fa` or `_hap1.fa`/`_hap2.fa` outputs are the
  expected input to `01_prepare_input_fasta.sh` here.

Author: Jaehyung Park (JP)
