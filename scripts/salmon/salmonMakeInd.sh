#!/bin/bash

Rscript salmonMakeSplici.R

mkdir refdata-gex-mm10-2020-A
mkdir refdata-gex-mm10-2020-A/genes

gffread /home/ggorin/ref/refdata-gex-mm10-2020-A/genes/genes.gtf -o refdata-gex-mm10-2020-A/genes/genes.gff 

grep "gene_name" refdata-gex-mm10-2020-A/genes/genes.gff | cut -f9 | cut -d';' -f2,3 | sed 's/=/ /g' | sed 's/;/ /g' | cut -d' ' -f2,4 | sort | uniq > geneid_to_name.txt 

./salmon-1.6.0_linux_x86_64/bin/salmon index -t transcriptome_mm10_splici_fl146/transcriptome_splici_fl146.fa -i mm10_splici_idx -p 16 

#Make hg38 index
mkdir refdata-gex-GRCh38-2020-A
mkdir refdata-gex-GRCh38-2020-A/genes 

gffread /home/ggorin/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf -o refdata-gex-GRCh38-2020-A/genes/genes.gff 

grep "gene_name" refdata-gex-GRCh38-2020-A/genes/genes.gff | cut -f9 | cut -d';' -f2,3 | sed 's/=/ /g' | sed 's/;/ /g' | cut -d' ' -f2,4$

./salmon-1.6.0_linux_x86_64/bin/salmon index -t transcriptome_hg38_splici_fl146/transcriptome_splici_fl146.fa -i hg38_splici_idx -p 16 


