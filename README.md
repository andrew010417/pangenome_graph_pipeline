# 판게놈 그래프 구축 및 Genotyping 파이프라인

문의사항은 Jaehyung Park(JP)에게 문의해주세요.

## 한 줄 요약

여러 사람의 개인별 genome assembly를 하나의 **판게놈 그래프**로 합치고,
그 그래프를 기준으로 새로운 short-read 샘플의 **변이를 찾아내는**
파이프라인입니다.

## 왜 필요한가

기존의 변이 발굴 방식은 **한 명의 레퍼런스 게놈(예: GRCh38)** 하고만
비교합니다. 이 방식은 그 레퍼런스 한 명과 다른 부분만 "변이"로
잡아내기 때문에, 레퍼런스에는 없는 희귀하거나 구조적으로 복잡한
변이는 애초에 비교 대상이 없어서 놓치기 쉽습니다.

**판게놈 그래프**는 한 명이 아니라 여러 명(수십~수백 명)의 게놈을 한꺼번에
하나의 그래프 자료구조로 합쳐놓은 것입니다. 그래프 안에는 여러 사람이
공유하는 서열은 하나의 경로로 겹쳐 있고, 사람마다 다른 부분은 그래프가
갈라지는(분기하는) 형태로 표현됩니다. 새로운 사람의 DNA를 이 그래프에
대조하면, 단일 레퍼런스로는 안 보이던 변이(이미 그래프 안의 누군가는
갖고 있던 변이)까지 훨씬 잘 잡아낼 수 있습니다.

이 저장소는 개인별 long-read assembly 파이프라인인
[`genome_assembly_ONT(hi-c)`](../genome_assembly_ONT(hi-c))의 다음
단계입니다: 그 파이프라인이 사람 한 명 한 명의 assembly(유전체 지도)를
만들면, 이 저장소는 그 assembly들 + 레퍼런스 게놈을 합쳐서 판게놈
그래프를 만들고, 그 그래프로 다른 short-read 샘플(예: 희귀질환
환자)의 변이를 찾는 데 씁니다.

```
genome_assembly_ONT(hi-c)              pangenome_graph_pipeline
(사람 한 명씩 assembly 생성)     →      (assembly들을 그래프로 합치고,
                                         그래프로 새 샘플의 변이를 발굴)

  sample01_hap1.fa  ─┐
  sample01_hap2.fa  ─┤
  sample02_hap1.fa  ─┼─→  [판게놈 그래프]  ──대조──→  새 샘플의 변이(VCF)
  sample02_hap2.fa  ─┤       ↑
  grch38.fa (레퍼런스)─┘   여러 사람의 정보가
                          다 들어있음
```

## 연구 배경 (관련 국내 연구 동향)

2026년 국내 long-read sequencing 학회에서도 이 방향(long-read assembly
→ 판게놈 → 희귀질환/구조변이 재발굴)이 여러 그룹에서 공통적으로
확인되었습니다:

- **GIST 박지환 교수님**은 발표에서 "이제 ONT를 이용해 이전에
  reference genome에만 의존하던 분석을 판게놈을 통해 population
  genetics 관점에서 진행할 수 있다"고 하시면서, Choi Lab이 보유한
  long-read 데이터로 판게놈을 만들고 이때 assembly 단계가 필요하다는
  점을 직접 언급하셨습니다. `genome_assembly_ONT(hi-c)`가 바로 그
  assembly 단계, 이 저장소가 그 다음의 판게놈 구축 단계에 해당합니다.
- **서울대병원 문장섭 교수님**은 short-read 기반 도구(Expansion
  Hunter 등)로는 놓쳤던 tandem repeat disorder 65종을 ONT
  long-read로 재검출한 사례를 발표하셨습니다 — 이 저장소가 만드는
  판게놈 그래프로 short-read 샘플을 재분석했을 때 기대하는 효과와
  같은 방향입니다.
- **연세대 윤지훈 교수님**은 한국인 56명의 genome을 ONT로 assembly해
  reference-free 방식(예: `asm2d6`)으로 임상적으로 중요한 유전자좌를
  분석한 사례를 보여주셨습니다.

## 핵심 개념 (용어가 낯설다면 먼저 읽어주세요)

- **그래프(Graph)**: 판게놈을 표현하는 자료구조. 서열 조각을 "노드",
  노드 사이의 연결을 "엣지"라고 부릅니다. 한 사람의 전체 게놈은 이
  그래프 위의 하나의 "경로(path)"로 표현됩니다. 사람마다 경로가 갈라지는
  지점이 곧 변이(SNP, indel, 구조변이 등)입니다.
- **PanSN 네이밍 (`sample#hap#contig`)**: 그래프 안에서 각 경로가
  "누구의 몇 번째 haplotype인지" 구분하기 위한 이름 규칙입니다. 예를 들어
  `sample01#1#chr1`은 sample01이라는 사람의 첫 번째 haplotype(엄마 또는
  아빠 쪽)의 chr1 경로를 뜻합니다. 같은 사람의 hap1/hap2를 같은
  `sample_id`로 묶어야 이후 단계에서 이 사람의 변이를 diploid(두 배수체,
  즉 유전자 쌍)로 올바르게 인식합니다.
- **GBZ / dist / min (그래프 인덱스 파일들)**: read를 그래프에 빠르게
  매핑하기 위해 미리 만들어두는 색인 파일들입니다. `.gbz`는 그래프
  자체를 압축해서 저장한 파일, `.dist`는 그래프 위 거리 정보,
  `.min`은 빠른 검색을 위한 minimizer(짧은 서열 조각) 색인입니다.
  일반 참조 게놈 매핑에서 쓰는 `.fai`/BWA 인덱스의 그래프 버전이라고
  생각하면 됩니다.
- **GAM 파일**: read가 그래프의 어느 경로/노드에 매핑됐는지 기록한
  파일입니다. 일반 매핑에서 쓰는 BAM 파일의 "그래프 버전"입니다.
- **Variant calling vs. Genotyping**: 이 두 단어를 혼용하기 쉬운데
  이 파이프라인에서는 구분해서 씁니다.
  - **Variant calling(변이 발굴, `vg` 경로)**: read가 실제로 어떻게
    매핑되는지 보고 "여기에 변이가 있다"를 새로 찾아냅니다. 그래프에
    아직 없던 변이도 발견할 수 있습니다.
  - **Genotyping(패널 기반 확인, `PanGenie` 경로)**: 이미 알려진 변이
    목록(패널)을 놓고 "이 사람이 이 변이를 갖고 있는지 없는지"만
    빠르게 확인합니다. 패널에 없는 변이는 애초에 찾을 수 없습니다.
- **Panel VCF(패널 VCF)**: 그래프 안에 들어있는 모든 변이 목록을 정리한
  VCF 파일입니다. PanGenie가 genotyping할 때 "확인해야 할 목록"으로
  사용합니다.
- **k-mer / HMM (PanGenie가 동작하는 방식)**: PanGenie는 read를 그래프에
  직접 정렬하지 않고, DNA를 일정 길이(k, 보통 31글자)로 잘라서 그
  조각(k-mer)들이 read에서 몇 번 등장하는지 센 다음, 통계 모델(HMM,
  은닉 마르코프 모델)로 "이 변이를 갖고 있을 확률"을 계산합니다. 그래서
  read를 그래프에 올리는 vg 방식보다 훨씬 빠릅니다.

## 파이프라인 출처에 대해

이 저장소의 스크립트는 예전에 어딘가에서 돌아가던 파이프라인을 그대로
가져온 것이 아니라, **Minigraph-Cactus / PGGB / vg / PanGenie 각 도구의
공식 CLI 문서를 기준으로 새로 작성**한 것입니다. 그래서 파일명, 정확한
플래그, 버전별 동작 방식은 소규모 테스트 데이터로 end-to-end 실행해서
확인하기 전까지는 검증되지 않은 것으로 간주해주세요 —
[VERIFICATION_TODO.md](VERIFICATION_TODO.md)에 확인이 필요한 항목을
정리해두었습니다.

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

### 각 스크립트가 정확히 무엇을 입출력하는지

| 스크립트 | 입력 | 출력 | 이 단계가 하는 일 |
|---|---|---|---|
| `01_prepare_input_fasta.sh` | assembly FASTA 여러 개 + manifest.tsv | `mc_seqfile.tsv`, `pggb_input.fa.gz` | 헤더를 `sample#hap#contig`로 통일하고, 두 그래프 도구가 요구하는 입력 형식을 각각 만듦 |
| `02a_build_graph_minigraph_cactus.sh` | `mc_seqfile.tsv` | `<outName>.vcf.gz`, `<outName>.gbz`, `<outName>.gfa.gz` | 그래프 구축 + VCF 패널 + Giraffe 인덱스까지 한 번에 생성 |
| `02b_build_graph_pggb.sh` | `pggb_input.fa.gz` | `smooth.final.gfa` | 그래프 구축만 (VCF는 안 나옴) |
| `03_pggb_gfa_to_vcf.sh` | `smooth.final.gfa` | `pggb.vcf.gz` | GFA 그래프를 VCF 변이 패널로 변환 |
| `04_vg_autoindex_giraffe.sh` | GFA, 또는 (레퍼런스 FASTA + VCF 패널) | `<prefix>.giraffe.gbz`, `.dist`, `.min` | read 매핑에 필요한 인덱스 생성 |
| `05_vg_giraffe_call.sh` | 그래프 인덱스 + 샘플 FASTQ(R1/R2) | `<sample>.vcf.gz` | read를 그래프에 매핑해서 그 샘플의 변이를 새로 발굴 |
| `06_prepare_pangenie_panel.sh` | VCF 패널 + 레퍼런스 FASTA | `<name>.biallelic.vcf.gz` | multiallelic 변이를 분리해서 PanGenie가 읽을 수 있는 형태로 정리 |
| `07_pangenie_genotype.sh` | biallelic 패널 VCF + 샘플 FASTQ | `<sample>_genotyping.vcf(.gz)` | 패널에 있는 변이들을 그 샘플이 갖고 있는지 k-mer 기반으로 판정 |

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
게놈마다 한 줄씩. `hap_index`는 phasing이 안 된(unphased) assembly면
`0`, phased diploid assembly면 `1`/`2`를 씁니다. 예:

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
  만약 패널 VCF에 개체별이 아니라 haplotype별로 컬럼이 따로 나온다면
  (`bcftools query -l`로 확인 가능), PanGenie genotyping 전에
  `06_prepare_pangenie_panel.sh`에 haplotype 병합 단계를 추가해야
  합니다.
- **`vg pack`의 GAM 입력 플래그**: `05_vg_giraffe_call.sh`는 GAM 입력에
  `-g`가 아니라 `-a`를 사용합니다 (이름만 보면 `-g`가 더 자연스러워
  보이지만 실제로는 `-a`가 맞는 경우가 있었습니다). 설치된 vg 버전의
  `vg pack --help`로 다시 확인하세요.
- **출력 파일명은 버전에 따라 달라질 수 있고 아직 검증되지 않았습니다**:
  `cactus-pangenome`(02a)과 `pggb`의 최종 smoothed GFA(02b) 모두
  예상 파일명을 로그로 남기지만, 실제 end-to-end 실행으로 확인된 적은
  없습니다.
- **리소스 요구량**: 사전 벤치마크(haplotype 2개, CHM13 + GRCh38
  기준) 상 Minigraph-Cactus는 약 7.5시간, peak 메모리 약 108GB(32
  threads)로 끝났고, PGGB는 peak 메모리는 더 낮지만 약 21시간이
  걸렸고 haplotype/thread가 늘어날수록 확장성이 훨씬 나빴습니다 (같은
  벤치마크에서 haplotype 100개 이상일 때 PGGB가 실패한 사례도
  있었습니다). 실제 개인별 diploid assembly와 더 많은 샘플로 돌리면
  이보다 훨씬 커질 수 있으니, 클러스터에 job을 올리기 전에 리소스를
  여유 있게 잡으세요.
- **이 저장소는 개인별 long-read assembly 단계를 다루지 않습니다** —
  그건 [`genome_assembly_ONT(hi-c)`](../genome_assembly_ONT(hi-c))의
  역할이며, 그 결과물인 `_primary.fa`/`_alternate.fa` 또는
  `_hap1.fa`/`_hap2.fa`가 여기 `01_prepare_input_fasta.sh`의 입력으로
  들어갑니다.

Author: Jaehyung Park (JP)
