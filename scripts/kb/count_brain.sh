#!/bin/bash
# generate count matrices 210622

kb count --verbose \
-i ../ref/refdata-gex-mm10-2020-A/kallisto/index.idx \
-g ../ref/refdata-gex-mm10-2020-A/t2g_mm10.txt \
-x 10xv3 \
-o ./brain_10x_5k/ \
-t 30 -m 30G \
-c1 ../ref/refdata-gex-mm10-2020-A/kallisto/cdna_t2c.txt \
-c2 ../ref/refdata-gex-mm10-2020-A/kallisto/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L001_R1_001.fastq.gz \
../datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L001_R2_001.fastq.gz \
../datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L002_R1_001.fastq.gz \
../datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L002_R2_001.fastq.gz \
../datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L003_R1_001.fastq.gz \
../datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L003_R2_001.fastq.gz \
../datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L004_R1_001.fastq.gz \
../datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L004_R2_001.fastq.gz \
