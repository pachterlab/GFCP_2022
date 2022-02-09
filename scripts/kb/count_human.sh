#!/bin/bash
# generate count matrices 210517

kb count --verbose \
-i ../ref/refdata-gex-GRCh38-2020-A/kallisto/index.idx \
-g ../ref/refdata-gex-GRCh38-2020-A/t2g_grch38.txt \
-x 10xv3 \
-o ./pbmc_1k_v3/ \
-t 16 -m 30G \
-c1 ../ref/refdata-gex-GRCh38-2020-A/kallisto/cdna_t2c.txt \
-c2 ../ref/refdata-gex-GRCh38-2020-A/kallisto/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../datasets/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz \
../datasets/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz \
../datasets/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz \
../datasets/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz

kb count --verbose \
-i ../ref/refdata-gex-GRCh38-2020-A/kallisto/index.idx \
-g ../ref/refdata-gex-GRCh38-2020-A/t2g_grch38.txt \
-x 10xv3 \
-o ./pbmc_10k_v3/ \
-t 16 -m 30G \
-c1 ../ref/refdata-gex-GRCh38-2020-A/kallisto/cdna_t2c.txt \
-c2 ../ref/refdata-gex-GRCh38-2020-A/kallisto/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../datasets/pbmc_10k_v3_fastqs/pbmc_10k_v3_S1_L001_R1_001.fastq.gz \
../datasets/pbmc_10k_v3_fastqs/pbmc_10k_v3_S1_L001_R2_001.fastq.gz \
../datasets/pbmc_10k_v3_fastqs/pbmc_10k_v3_S1_L002_R1_001.fastq.gz \
../datasets/pbmc_10k_v3_fastqs/pbmc_10k_v3_S1_L002_R2_001.fastq.gz

