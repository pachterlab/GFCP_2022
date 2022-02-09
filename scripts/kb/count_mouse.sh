#!/bin/bash
# generate count matrices 210518

kb count --verbose \
-i ../ref/refdata-gex-mm10-2020-A/kallisto/index.idx \
-g ../ref/refdata-gex-mm10-2020-A/t2g_mm10.txt \
-x 10xv3 \
-o ./neuron_1k_v3/ \
-t 30 -m 30G \
-c1 ../ref/refdata-gex-mm10-2020-A/kallisto/cdna_t2c.txt \
-c2 ../ref/refdata-gex-mm10-2020-A/kallisto/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../datasets/neuron_1k_v3_fastqs/neuron_1k_v3_S1_L001_R1_001.fastq.gz \
../datasets/neuron_1k_v3_fastqs/neuron_1k_v3_S1_L001_R2_001.fastq.gz \
../datasets/neuron_1k_v3_fastqs/neuron_1k_v3_S1_L002_R1_001.fastq.gz \
../datasets/neuron_1k_v3_fastqs/neuron_1k_v3_S1_L002_R2_001.fastq.gz

kb count --verbose \
-i ../ref/refdata-gex-mm10-2020-A/kallisto/index.idx \
-g ../ref/refdata-gex-mm10-2020-A/t2g_mm10.txt \
-x 10xv3 \
-o ./neuron_10k_v3/ \
-t 30 -m 30G \
-c1 ../ref/refdata-gex-mm10-2020-A/kallisto/cdna_t2c.txt \
-c2 ../ref/refdata-gex-mm10-2020-A/kallisto/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../datasets/neuron_10k_v3_fastqs/neuron_10k_v3_S1_L001_R1_001.fastq.gz \
../datasets/neuron_10k_v3_fastqs/neuron_10k_v3_S1_L001_R2_001.fastq.gz \
../datasets/neuron_10k_v3_fastqs/neuron_10k_v3_S1_L002_R1_001.fastq.gz \
../datasets/neuron_10k_v3_fastqs/neuron_10k_v3_S1_L002_R2_001.fastq.gz

kb count --verbose \
-i ../ref/refdata-gex-mm10-2020-A/kallisto/index.idx \
-g ../ref/refdata-gex-mm10-2020-A/t2g_mm10.txt \
-x 10xv3 \
-o ./heart_1k_v3/ \
-t 30 -m 30G \
-c1 ../ref/refdata-gex-mm10-2020-A/kallisto/cdna_t2c.txt \
-c2 ../ref/refdata-gex-mm10-2020-A/kallisto/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../datasets/heart_1k_v3_fastqs/heart_1k_v3_S1_L001_R1_001.fastq.gz \
../datasets/heart_1k_v3_fastqs/heart_1k_v3_S1_L001_R2_001.fastq.gz \
../datasets/heart_1k_v3_fastqs/heart_1k_v3_S1_L002_R1_001.fastq.gz \
../datasets/heart_1k_v3_fastqs/heart_1k_v3_S1_L002_R2_001.fastq.gz

kb count --verbose \
-i ../ref/refdata-gex-mm10-2020-A/kallisto/index.idx \
-g ../ref/refdata-gex-mm10-2020-A/t2g_mm10.txt \
-x 10xv3 \
-o ./heart_10k_v3/ \
-t 30 -m 30G \
-c1 ../ref/refdata-gex-mm10-2020-A/kallisto/cdna_t2c.txt \
-c2 ../ref/refdata-gex-mm10-2020-A/kallisto/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../datasets/heart_10k_v3_fastqs/heart_10k_v3_S1_L001_R1_001.fastq.gz \
../datasets/heart_10k_v3_fastqs/heart_10k_v3_S1_L001_R2_001.fastq.gz \
../datasets/heart_10k_v3_fastqs/heart_10k_v3_S1_L002_R1_001.fastq.gz \
../datasets/heart_10k_v3_fastqs/heart_10k_v3_S1_L002_R2_001.fastq.gz




