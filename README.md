# 판게놈 그래프 구축 및 Genotyping 파이프라인

문의사항은 Jaehyung Park(JP)에게 문의해주세요.

여러 사람의 개인별 genome assembly(예: [`genome_assembly_ONT(hi-c)`](../genome_assembly_ONT(hi-c))의
결과물)와 레퍼런스 게놈을 **Minigraph-Cactus** 그리고/또는 **PGGB**로 합쳐서
하나의 판게놈 그래프를 만들고, 그 그래프를 이용해 개별 short-read 샘플의
변이를 **vg**(read를 그래프에 직접 매핑) 그리고/또는 **PanGenie**(k-mer
기반 패널 genotyping)로 찾아내는 파이프라인입니다.

개인별 long-read assembly 다음 단계로 이어지도록 설계했습니다:
`genome_assembly_ONT(hi-c)`가 개인별 assembly를 만들면, 이 저장소는 그
assembly들을 하나의 그래프로 합치고, 그 그래프로 short-read 샘플에서
변이(단일 레퍼런스로는 놓칠 수 있는 희귀 변이 포함)를 찾는 데 씁니다.

## 배경 (Provenance)

이전에 KOGIC(국바빅 범유전체 프로젝트)에서 진행했던 판게놈 그래프 구축
작업을 참고해서 만들었습니다. 그 프로젝트에서는 CHM13/GRCh38를 대상으로
Minigraph-Cactus와 PGGB를 비교하고, `vg`와 `PanGenie`로 GIAB HG002
데이터에 대한 변이 발굴/genotyping까지 진행한 바 있습니다. 당시 자체
벤치마크(내부 기획 자료 참고)에 따르면, 샘플 수가 많아질수록
Minigraph-Cactus가 훨씬 실용적이었습니다 — 더 빠르고, peak 메모리도
훨씬 낮고, 확장성도 더 좋았습니다. 반면 PGGB는 reference-free
all-vs-all 정렬 방식이라 더 세밀하게 가지 친(branch) 그래프(엣지가 더
많고 bubble 구조가 많음)를 만들지만, 그만큼 시간/메모리 비용이 크고
규모가 커지면 안정성이 떨어졌습니다 (이전 벤치마크에서 haplotype
100개 이상에서 실패 사례 관찰).

`genome_assembly_ONT(hi-c)`와 달리, 이 프로젝트의 원본 스크립트는
확보하지 못했고 요약 노트와 벤치마크 수치만 참고할 수 있었습니다.
**이 저장소의 스크립트들은 예전에 돌아가던 파이프라인을 그대로 복사한
것이 아니라, 각 도구의 공식 CLI 문서를 기준으로 새로 작성**한 것입니다.
파일명, 정확한 플래그, 버전별 동작 방식은 소규모 테스트 데이터로
end-to-end 실행해서 확인하기 전까지는 검증되지 않은 것으로 간주해주세요
— [VERIFICATION_TODO.md](VERIFICATION_TODO.md) 참고.

## 요구 사항

- `cactus` (`cactus-pangenome` 포함, 9.x 이상) — conda로 설치 **불가**,
  설치 방법은 `environment.yml` 참고
- `pggb`, `vg`, `PanGenie`, `bcftools`, `samtools`, `htslib` (tabix/bgzip), `seqkit`

Cactus를 제외한 나머지 도구들은 버전을 고정한 conda 환경을
`environment.yml`로 제공합니다:

```bash
mamba env create -f environment.yml   # 또는: conda env create -f environment.yml
conda activate pangenome-graph
```

스크립트는 스케줄러에 종속되지 않는 순수 bash입니다 (SGE/SLURM 헤더
없음). 클러스터에서 필요하면 `qsub`/`sbatch`로 직접 감싸서 쓰세요.
스레드 수(`-t`)는 각 스크립트마다 보수적인 기본값이 들어있으니 실제
서버 코어 수에 맞게 조정하세요.

실제 랩 서버에서 처음 실행하기 전에 [VERIFICATION_TODO.md](VERIFICATION_TODO.md)의
검증 항목을 먼저 확인해주세요.

## 파이프라인 구성

```
scripts/
├── lib/common.sh                       공통 helper (로깅, 리소스 모니터링)
├── 01_prepare_input_fasta.sh           PanSN 규칙 헤더 정리 + MC seqFile / PGGB 병합 FASTA 생성
├── 02a_build_graph_minigraph_cactus.sh cactus-pangenome  -> VCF 패널 + Giraffe GBZ (기본 권장 경로)
├── 02b_build_graph_pggb.sh             pggb               -> GFA 그래프 (보조/비교 경로)
├── 03_pggb_gfa_to_vcf.sh               vg convert + vg deconstruct  -> VCF 패널 (PGGB 경로 전용)
├── 04_vg_autoindex_giraffe.sh          vg autoindex --workflow giraffe -> GBZ/dist/min (PGGB 경로, 또는 MC 인덱스 재생성용)
├── 05_vg_giraffe_call.sh               vg giraffe -> vg pack -> vg call   (샘플별, 그래프 직접 매핑)
└── 06_prepare_pangenie_panel.sh        bcftools norm -m -any  -> biallelic 패널 VCF
└── 07_pangenie_genotype.sh             PanGenie k-mer counting + HMM genotyping  (샘플별, 패널 기반)
```

**어떤 그래프 구축 도구를 쓸지**: 기본은 Minigraph-Cactus(`02a`)를
사용하세요 — 더 빠르고, 메모리도 적게 쓰고, 샘플 수가 많아도 잘
버팁니다. PGGB(`02b` + `03`)는 샘플 수가 적거나, 더 완전하고
reference-free한 variation 표현이 추가 비용을 감수할 만큼 중요한
locus를 볼 때, 또는 두 방법을 서로 비교하고 싶을 때 사용하세요.

**어떤 변이 발굴 경로를 쓸지**: `vg`(`05`)는 read가 그래프에 매핑되는
방식을 그대로 이용해 변이를 직접 발굴합니다 — 그래프에 아직 없는
변이도 찾을 수 있습니다. `PanGenie`(`06` + `07`)는 대신 이미 알려진
그래프 변이 목록(패널)을 샘플의 k-mer 정보와 대조해서 확인합니다 —
더 빠르지만 패널에 있는 것만 찾을 수 있습니다. 같은 샘플에 두 방법을
모두 돌려서 서로 교차 검증하는 것도 합리적인 방법입니다.

### Minigraph-Cactus 경로 (기본 권장)

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

### PGGB 경로 (보조 / 비교용)

```bash
scripts/01_prepare_input_fasta.sh -o prep -m samples_manifest.tsv -t 8

scripts/02b_build_graph_pggb.sh -o graph_pggb -i prep/pggb_input.fa.gz -n <haplotype 개수> -t 32

scripts/03_pggb_gfa_to_vcf.sh -o graph_pggb -g graph_pggb/smooth.final.gfa -r grch38 -t 16

scripts/04_vg_autoindex_giraffe.sh \
    -o graph_pggb -p pangenome_pggb -r ref/grch38.fa -v graph_pggb/pggb.vcf.gz -t 16

scripts/05_vg_giraffe_call.sh \
    -o calls_vg_pggb/sample01 -x graph_pggb/pangenome_pggb -s sample01 \
    -1 sample01_R1.fastq.gz -2 sample01_R2.fastq.gz -t 16
```

`samples_manifest.tsv` 형식 (`01_prepare_input_fasta.sh` 상단 주석
참고): `sample_id<TAB>hap_index<TAB>fasta_path`, haplotype/레퍼런스
게놈마다 한 줄씩. 예:

```
sample01	1	assemblies/sample01_hap1.fa
sample01	2	assemblies/sample01_hap2.fa
grch38	0	ref/grch38.fa
```

## 참고 사항 / 주의할 점

- **Cactus는 여기서 conda로 설치하지 않습니다** — 이유와 설치 방법은
  `environment.yml` 참고.
- **PanSN 네이밍이 diploid genotyping에 중요합니다**:
  `01_prepare_input_fasta.sh`가 헤더를 `sample#hap#contig` 형식으로
  바꿉니다. 같은 개체의 hap1/hap2를 같은 `sample_id` 아래로 묶어야
  `vg deconstruct`(03)와 `cactus-pangenome`(02a)가 개체당 하나의
  diploid genotype 컬럼을 자동으로 만들어줍니다 — 이 파이프라인에서
  실제로 확인된 사항은 아직 아닙니다 (VERIFICATION_TODO.md 참고).
  만약 패널 VCF에 개체별이 아니라 haplotype별로 컬럼이 따로 나온다면,
  PanGenie genotyping 전에 `06_prepare_pangenie_panel.sh`에 haplotype
  병합 단계를 추가해야 합니다.
- **`vg pack`의 GAM 입력 플래그**: `05_vg_giraffe_call.sh`는 GAM 입력에
  `-g`가 아니라 `-a`를 사용합니다 (이전에 직접 겪었던 사항을 반영한
  것으로, 이름만 보면 `-g`가 더 자연스러워 보이지만 실제로는 `-a`가
  맞았습니다). 설치된 vg 버전의 `vg pack --help`로 다시 확인하세요.
- **출력 파일명은 버전에 따라 달라질 수 있고 아직 검증되지 않았습니다**:
  `cactus-pangenome`(02a)과 `pggb`의 최종 smoothed GFA(02b) 모두
  예상 파일명을 로그로 남기지만, 실제 end-to-end 실행으로 확인된 적은
  없습니다.
- **리소스 요구량**: 이전 벤치마크(haplotype 2개, CHM13 + GRCh38
  기준) 상 Minigraph-Cactus는 약 7.5시간, peak 메모리 약 108GB(32
  threads)로 끝났고, PGGB는 peak 메모리는 더 낮지만 약 21시간이
  걸렸고 haplotype/thread가 늘어날수록 확장성이 훨씬 나빴습니다.
  실제 개인별 diploid assembly와 더 많은 샘플로 돌리면 이보다 훨씬
  커질 수 있으니, 클러스터에 job을 올리기 전에 리소스를 여유 있게
  잡으세요.
- **이 저장소는 개인별 long-read assembly 단계를 다루지 않습니다** —
  그건 [`genome_assembly_ONT(hi-c)`](../genome_assembly_ONT(hi-c))의
  역할이며, 그 결과물인 `_primary.fa`/`_alternate.fa` 또는
  `_hap1.fa`/`_hap2.fa`가 여기 `01_prepare_input_fasta.sh`의 입력으로
  들어갑니다.

Author: Jaehyung Park (JP)
