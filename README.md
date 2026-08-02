# 판게놈 그래프 구축 및 Genotyping 파이프라인

문의사항은 Jaehyung Park(JP)에게 문의해주세요.

## 요약

Oxford Nanopore(ONT)로 대표되는 long-read 시퀀싱은 가격과 정확도가
빠르게 개선되면서 희귀질환 진단·인구집단 유전체 연구의 핵심 기술로
자리잡아가고 있습니다. 기존 short-read 기반 분석(Illumina, Q>40 수준의
높은 base-level 정확도)은 read 길이가 150bp 안팎으로 짧다는 근본적인
한계 때문에 다음을 구조적으로 놓칩니다:

- **구조변이(SV)** — 수백~수천 bp 크기의 삽입/결실은 짧은 read로 통째로
  걸쳐 읽을 수 없어 간접적인 신호로만 추정 가능
- **Tandem repeat 확장(expansion) 질환** — 반복서열이 길어지는 질환은
  Expansion Hunter 같은 short-read 전용 도구로도 상당수 놓침
- **Transcript isoform** — short-read RNA-seq은 조각 단위(300~400nt)로만
  읽어서 전체 isoform 구조를 재구성하기 어려움
- **RNA 화학적 변형(modification)** — 170종 이상 존재하는 RNA
  modification은 기존 antibody pulldown/bisulfite 방식으로는 검출
  범위와 비용에 한계가 있었음 (최근 ONT single-molecule 시퀀싱 +
  딥러닝 기반 modification 검출 모델로 이 한계가 풀리고 있음)
- **레퍼런스에 없는 서열** — 단일 레퍼런스 게놈(GRCh38) 하나와만
  비교하는 방식이라, 레퍼런스에 존재하지 않는 서열은 애초에 비교
  대상이 없어 변이로 인식되지 않음 (reference bias)

ONT long-read(HiFi 포함)는 read 하나가 길어서 이런 영역을 통째로 읽고
조립(assembly)할 수 있고, telomere 길이 측정, T2T 수준 assembly,
single-cell 수준 long-read 시퀀싱, 추가 라이브러리 준비 없이 특정
영역만 골라 읽는 adaptive sampling 등으로 활용 범위가 계속
넓어지고 있습니다. 이런 흐름에서, 개별 long-read assembly들을 모아
**판게놈(pangenome) 그래프**로 합치는 것은 단일 레퍼런스 기반 분석의
자연스러운 다음 단계로 여겨집니다 — 한 명이 아니라 여러 명(수십~수백
명)의 게놈을 하나의 그래프 자료구조에 담아두면, 그래프 안에서 여러
사람이 공유하는 서열은 하나의 경로로 겹치고 사람마다 다른 부분은
그래프가 갈라지는(분기하는) 형태로 표현됩니다.

**이 프로젝트의 목표는 바로 이 지점입니다**: short-read 기반 분석으로는
찾지 못했던 변이(희귀질환 관련 변이 포함)를, long-read로 구축한 판게놈
그래프를 이용해 찾아내는 것입니다. 단일 레퍼런스 하나와만 비교하던
기존 방식과 달리, 여러 사람의 long-read assembly가 담긴 그래프에
샘플을 대조하면 레퍼런스에는 없었던(그래서 예전엔 변이로조차 인식되지
못했던) 서열까지 비교 대상에 포함되기 때문에, 그동안 놓쳤던 변이를
훨씬 잘 잡아낼 수 있습니다.

이 저장소는 그 흐름의 후반부, 즉 **판게놈 그래프 구축과 genotyping**
단계를 구현합니다. 개인별 long-read assembly 파이프라인인
`genome_assembly_ONT(hi-c)` (별도 저장소)가 사람 한
명 한 명의 assembly(유전체 지도)를 만들면, 이 저장소는 그
assembly들 + 레퍼런스 게놈을 Minigraph-Cactus/PGGB로 합쳐서 판게놈
그래프를 만들고, vg와 PanGenie로 새로운 short-read 샘플(예: 희귀질환
환자)의 변이를 찾아냅니다.

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

## 왜 "그래프"가 필요한가 (단일 레퍼런스 + long-read만으로는 부족한 이유)

long read를 레퍼런스 하나에 그냥 매핑하는 것만으로도 read 길이가
짧아서 생기는 문제(반복구간 통과 실패, 구조변이를 못 걸쳐서 읽는 문제
등)는 상당 부분 해결됩니다. 하지만 그것만으로는 안 풀리는 문제가
따로 있고, 그게 이 파이프라인이 판게놈 **그래프**를 만드는 이유입니다.

- **레퍼런스에 아예 없는 구조는 애초에 "붙일 자리"가 없습니다.** 어떤
  사람이 GRCh38에 없는 큰 삽입 서열이나 아예 다른 구조(순서가 뒤집힌
  배열 등)를 갖고 있으면, read가 아무리 길어도 "레퍼런스의 이
  위치에 정렬한다"는 방식 자체가 그 사람 구조를 표현하지 못합니다.
  read가 잘려서 붙거나(soft-clip), 억지로 끼워 맞춰져 왜곡된 변이로
  잡히거나, 아예 정렬되지 않습니다. read 길이 문제가 아니라 **비교
  기준(좌표계) 자체가 그 사람의 구조를 표현하지 못하는 문제**입니다.
- **"흔한 변이인지 희귀질환 원인 변이인지" 판단이 안 됩니다.** 희귀질환
  연구에서는 환자에게서 레퍼런스와 다른 부분을 찾는 것만으로는
  부족하고, 그게 건강한 사람들도 흔히 갖는 정상 변이인지 진짜 희귀한
  병인성 변이인지 구분해야 합니다. 레퍼런스 한 명과만 비교해서는
  "다르다"는 것만 알 뿐 "얼마나 흔한지"는 알 수 없습니다. 판게놈
  그래프는 여러 사람의 정보가 이미 담겨 있어서, 그래프 자체가 "이
  부위는 사람들 사이에 이렇게 다양하다"는 인구집단 정보를 갖고
  있습니다 — 그래프에 대조하면 "이건 흔한 경로다 / 이건 아무도
  없던 경로다"까지 함께 알 수 있습니다.
- **재사용성**: 그래프를 한 번 잘 만들어두면, 이후 오는 환자들(특히
  long-read 없이 short-read만 있는 경우도)마다 매번 복잡한 구조를
  새로 분석할 필요 없이 이 그래프에 대조해서 재사용할 수 있습니다.

즉, long read를 단일 레퍼런스에 매핑하는 것이 "더 정확하게 읽는 것"
이라면, 판게놈 그래프는 "비교 기준 자체를 여러 사람 것으로 넓히는
것"입니다. 전자만으로는 "이 환자가 레퍼런스와 다르다"는 것만 알 수
있고, 후자가 있어야 "이게 정상적인 사람들의 다양성 범위 안인지,
진짜 이상한 것인지"까지 판단할 수 있습니다.

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

## 요구 사항

- `cactus` (`cactus-pangenome` 서브커맨드가 있는 버전) — conda로 설치
  **불가**, 설치 방법은 `environment.yml` 참고. 정확한 버전 요구사항은
  미확인 상태이니 설치 후 `cactus-pangenome --help`로 지원 여부를
  먼저 확인하세요.
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
├── 06_prepare_pangenie_panel.sh        bcftools norm -m -any  -> biallelic 패널 VCF
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
  그건 `genome_assembly_ONT(hi-c)`(별도 저장소)의
  역할이며, 그 결과물인 `_primary.fa`/`_alternate.fa` 또는
  `_hap1.fa`/`_hap2.fa`가 여기 `01_prepare_input_fasta.sh`의 입력으로
  들어갑니다.

Author: Jaehyung Park (JP)
