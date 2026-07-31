# 랩 서버 실행 전 검증 필요 항목

이 파이프라인은 이전 프로젝트의 요약 노트/벤치마크 수치만 참고해서
새로 작성한 스크립트입니다 (원본 스크립트는 확보하지 못함). 실제
랩 데이터로 돌리기 전에 아래 항목을 소규모 테스트 데이터로 먼저
확인해주세요.

- [ ] `cactus-pangenome` 설치 버전에서 아래 플래그가 그대로 지원되는지
      확인 (`cactus-pangenome --help`): `--outDir`, `--outName`,
      `--reference`, `--vcf`, `--gbz`, `--gfa`, `--giraffe`, `--maxCores`
      참고: scripts/02a_build_graph_minigraph_cactus.sh
- [ ] cactus-pangenome 실행 후 실제 출력 파일명이 `<outName>.vcf.gz`,
      `<outName>.gbz`, `<outName>.gfa.gz` 형태로 나오는지 확인
      참고: scripts/02a_build_graph_minigraph_cactus.sh
- [ ] `pggb` 설치 버전에서 최종 smoothed GFA 파일명이
      `*smooth.final.gfa` 패턴을 따르는지 확인, 아니면 find 패턴 수정
      참고: scripts/02b_build_graph_pggb.sh
- [ ] `vg convert -g <gfa> -p`, `vg deconstruct -P <ref> -e` 플래그가
      설치된 vg 버전에서 동일하게 동작하는지 확인 (`vg convert --help`,
      `vg deconstruct --help`)
      참고: scripts/03_pggb_gfa_to_vcf.sh
- [ ] `vg autoindex --workflow giraffe` 의 `-g`(gfa) / `-r`+`-v`(ref+vcf)
      입력 모드가 설치된 vg 버전에서 동일하게 지원되는지 확인
      참고: scripts/04_vg_autoindex_giraffe.sh
- [ ] `vg pack` 의 GAM 입력 플래그가 실제로 `-a` 인지 (`-g` 아님) —
      이전 경험 기반으로 넣었지만 vg 버전마다 바뀔 수 있음, `vg pack --help`
      로 재확인
      참고: scripts/05_vg_giraffe_call.sh
- [ ] `01_prepare_input_fasta.sh` 에서 같은 개체의 hap1/hap2를
      `sample#1#contig` / `sample#2#contig` 로 이름 붙였을 때,
      `vg deconstruct`/`cactus-pangenome` 출력 VCF가 개체당 1개의
      diploid genotype 컬럼으로 나오는지 `bcftools query -l` 로 확인.
      만약 haplotype별로 컬럼이 따로 나오면 `06_prepare_pangenie_panel.sh`
      에 haplotype 병합 로직을 추가해야 함
      참고: scripts/01_prepare_input_fasta.sh, scripts/06_prepare_pangenie_panel.sh
- [ ] `PanGenie` 설치 버전에서 `-i/-r/-v/-o/-s/-t/-j` 플래그가 동일하게
      지원되는지 확인 (`PanGenie --help`) — v3+ 부터 index/genotype
      바이너리가 분리되었을 수 있음
      참고: scripts/07_pangenie_genotype.sh
- [ ] `environment.yml` 의 버전 핀(pggb, vg, pangenie, bcftools 등)이
      bioconda에서 실제로 resolve 되는지, 랩 HPC 기존 설치 버전과
      맞는지 확인 후 다르면 실제 버전으로 맞추기
      참고: environment.yml
- [ ] Cactus 설치 (conda 미지원) — 랩 서버에 어떤 방식(바이너리
      tarball / Docker / Singularity)으로 설치할지 결정하고 문서화
      참고: environment.yml
- [ ] `set -euo pipefail` 적용 후 각 스크립트의 파이프(`| tee`) 실패
      시 실제로 스크립트가 멈추는지 확인 (의도적 실패 케이스로 테스트)
      참고: 모든 scripts/*.sh
