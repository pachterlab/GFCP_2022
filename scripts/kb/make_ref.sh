#!/bin/bash
# construct t2g files
mkdir -p ./refdata-gex-GRCh38-2020-A/kallisto
mkdir -p ./refdata-gex-mm10-2020-A/kallisto
kb ref -i ./refdata-gex-GRCh38-2020-A/kallisto/index.idx -g ./refdata-gex-GRCh38-2020-A/t2g_grch38.txt -f1 ./refdata-gex-GRCh38-2020-A/kallisto/cdna.fa -f2 ./refdata-gex-GRCh38-2020-A/kallisto/intron.fa -c1 ./refdata-gex-GRCh38-2020-A/kallisto/cdna_t2c.txt -c2 ./refdata-gex-GRCh38-2020-A/kallisto/intron_t2c.txt --workflow lamanno ./refdata-gex-GRCh38-2020-A/fasta/genome.fa ./refdata-gex-GRCh38-2020-A/genes/genes.gtf
kb ref -i ./refdata-gex-mm10-2020-A/kallisto/index.idx -g ./refdata-gex-mm10-2020-A/t2g_mm10.txt -f1 ./refdata-gex-mm10-2020-A/kallisto/cdna.fa -f2 ./refdata-gex-mm10-2020-A/kallisto/intron.fa -c1 ./refdata-gex-mm10-2020-A/kallisto/cdna_t2c.txt -c2 ./refdata-gex-mm10-2020-A/kallisto/intron_t2c.txt --workflow lamanno ./refdata-gex-mm10-2020-A/fasta/genome.fa ./refdata-gex-mm10-2020-A/genes/genes.gtf 
