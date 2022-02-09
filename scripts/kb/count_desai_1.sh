#!/bin/bash

# generate count matrices 211013

kb count --verbose \
-i ../ref/refdata-gex-mm10-2020-A/kallisto/index.idx \
-g ../ref/refdata-gex-mm10-2020-A/t2g_mm10.txt \
-x 10xv2 \
-o ./desai_idu/ \
-t 30 -m 30G \
-c1 ../ref/refdata-gex-mm10-2020-A/kallisto/cdna_t2c.txt \
-c2 ../ref/refdata-gex-mm10-2020-A/kallisto/intron_t2c.txt \
--workflow lamanno --filter bustools --overwrite --loom \
../datasets/desai_idu/idu/SRR14713296.1_2.fastq \
../datasets/desai_idu/idu/SRR14713296.1_3.fastq
