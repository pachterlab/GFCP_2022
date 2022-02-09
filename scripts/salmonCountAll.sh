#!/bin/bash

./salmon-1.6.0_linux_x86_64/bin/salmon alevin -i mm10_splici_idx -p 20 -l ISR --chromiumV3  --sketch \
-1 /home/ggorin/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L001_R1_001.fastq.gz /home/ggorin/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L002_R1_001.fastq.gz /home/ggorin/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L003_R1_001.fastq.gz /home/ggorin/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L004_R1_001.fastq.gz \
-2 /home/ggorin/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L001_R2_001.fastq.gz /home/ggorin/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L002_R2_001.fastq.gz /home/ggorin/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L003_R2_001.fastq.gz /home/ggorin/datasets/brain_10x_5k_fastqs/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L004_R2_001.fastq.gz \
-o brain_10x_5k

alevin-fry generate-permit-list -d fw -k -i brain_10x_5k -o brain_10x_5k_quant

alevin-fry collate -t 16 -i brain_10x_5k_quant -r brain_10x_5k

alevin-fry quant -t 16 -i brain_10x_5k_quant -o brain_10x_5k_res --tg-map transcriptome_mm10_splici_fl146/transcriptome_splici_fl146_t2g_3col.tsv --resolution cr-like --use-mtx



#brain_nuc_10x_5k

./salmon-1.6.0_linux_x86_64/bin/salmon alevin -i mm10_splici_idx -p 20 -l ISR --chromiumV3  --sketch \
-1 /home/ggorin/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L001_R1_001.fastq.gz /home/ggorin/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L002_R1_001.fastq.gz /home/ggorin/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L003_R1_001.fastq.gz /home/ggorin/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L004_R1_001.fastq.gz \
-2 /home/ggorin/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L001_R2_001.fastq.gz /home/ggorin/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L002_R2_001.fastq.gz /home/ggorin/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L003_R2_001.fastq.gz /home/ggorin/datasets/brain_nuc_10x_5k_fastqs/SC3_v3_NextGem_DI_Nuclei_5K_gex_S6_L004_R2_001.fastq.gz \
-o brain_nuc_10x_5k

alevin-fry generate-permit-list -d fw -k -i brain_nuc_10x_5k -o brain_nuc_10x_5k_quant

alevin-fry collate -t 16 -i brain_nuc_10x_5k_quant -r brain_nuc_10x_5k

alevin-fry quant -t 16 -i brain_nuc_10x_5k_quant -o brain_nuc_10x_5k_res --tg-map transcriptome_mm10_splici_fl146/transcriptome_splici_fl146_t2g_3col.tsv --resolution cr-like --use-mtx



#desai_idu_v2

./salmon-1.6.0_linux_x86_64/bin/salmon alevin -i mm10_splici_idx -p 20 -l ISR --chromium  --sketch \
-1 /home/ggorin/datasets/desai_idu/dmso/SRR14713295_S1_L001_R1_001.fastq \
-2 /home/ggorin/datasets/desai_idu/dmso/SRR14713295_S1_L001_R2_001.fastq \
-o desai_dmso_v2

alevin-fry generate-permit-list -d fw -k -i desai_dmso_v2 -o desai_dmso_v2_quant

alevin-fry collate -t 16 -i desai_dmso_v2_quant -r desai_dmso_v2

alevin-fry quant -t 16 -i desai_dmso_v2_quant -o desai_dmso_v2_res --tg-map transcriptome_mm10_splici_fl146/transcriptome_splici_fl146_t2g_3col.tsv --resolution cr-like --use-mtx



#desai_idu_v2

./salmon-1.6.0_linux_x86_64/bin/salmon alevin -i mm10_splici_idx -p 20 -l ISR --chromium  --sketch \
-1 /home/ggorin/datasets/desai_idu/idu/SRR14713296_S1_L001_R1_001.fastq \
-2 /home/ggorin/datasets/desai_idu/idu/SRR14713296_S1_L001_R2_001.fastq \
-o desai_idu_v2

alevin-fry generate-permit-list -d fw -k -i desai_idu_v2 -o desai_idu_v2_quant

alevin-fry collate -t 16 -i desai_idu_v2_quant -r desai_idu_v2

alevin-fry quant -t 16 -i desai_idu_v2_quant -o desai_idu_v2_res --tg-map transcriptome_mm10_splici_fl146/transcriptome_splici_fl146_t2g_3col.tsv --resolution cr-like --use-mtx



#heart_10k_v3

./salmon-1.6.0_linux_x86_64/bin/salmon alevin -i mm10_splici_idx -p 20 -l ISR --chromiumV3  --sketch \
-1 /home/ggorin/datasets/heart_10k_v3_fastqs/heart_10k_v3_S1_L001_R1_001.fastq.gz /home/ggorin/datasets/heart_10k_v3_fastqs/heart_10k_v3_S1_L002_R1_001.fastq.gz \
-2 /home/ggorin/datasets/heart_10k_v3_fastqs/heart_10k_v3_S1_L001_R2_001.fastq.gz /home/ggorin/datasets/heart_10k_v3_fastqs/heart_10k_v3_S1_L002_R2_001.fastq.gz \
-o heart_10k_v3

alevin-fry generate-permit-list -d fw -k -i heart_10k_v3 -o heart_10k_v3_quant

alevin-fry collate -t 16 -i heart_10k_v3_quant -r heart_10k_v3

alevin-fry quant -t 16 -i heart_10k_v3_quant -o heart_10k_v3_res --tg-map transcriptome_mm10_splici_fl146/transcriptome_splici_fl146_t2g_3col.tsv --resolution cr-like --use-mtx



#heart_1k_v3

./salmon-1.6.0_linux_x86_64/bin/salmon alevin -i mm10_splici_idx -p 20 -l ISR --chromiumV3  --sketch \
-1 /home/ggorin/datasets/heart_1k_v3_fastqs/heart_1k_v3_S1_L001_R1_001.fastq.gz /home/ggorin/datasets/heart_1k_v3_fastqs/heart_1k_v3_S1_L002_R1_001.fastq.gz \
-2 /home/ggorin/datasets/heart_1k_v3_fastqs/heart_1k_v3_S1_L001_R2_001.fastq.gz /home/ggorin/datasets/heart_1k_v3_fastqs/heart_1k_v3_S1_L002_R2_001.fastq.gz \
-o heart_1k_v3

alevin-fry generate-permit-list -d fw -k -i heart_1k_v3 -o heart_1k_v3_quant

alevin-fry collate -t 16 -i heart_1k_v3_quant -r heart_1k_v3

alevin-fry quant -t 16 -i heart_1k_v3_quant -o heart_1k_v3_res --tg-map transcriptome_mm10_splici_fl146/transcriptome_splici_fl146_t2g_3col.tsv --resolution cr-like --use-mtx



#neuron_10k_v3

./salmon-1.6.0_linux_x86_64/bin/salmon alevin -i mm10_splici_idx -p 20 -l ISR --chromiumV3  --sketch \
-1 /home/ggorin/datasets/neuron_10k_v3_fastqs/neuron_10k_v3_S1_L001_R1_001.fastq.gz /home/ggorin/datasets/neuron_10k_v3_fastqs/neuron_10k_v3_S1_L002_R1_001.fastq.gz \
-2 /home/ggorin/datasets/neuron_10k_v3_fastqs/neuron_10k_v3_S1_L001_R2_001.fastq.gz /home/ggorin/datasets/neuron_10k_v3_fastqs/neuron_10k_v3_S1_L002_R2_001.fastq.gz \
-o neuron_10k_v3

alevin-fry generate-permit-list -d fw -k -i neuron_10k_v3 -o neuron_10k_v3_quant

alevin-fry collate -t 16 -i neuron_10k_v3_quant -r neuron_10k_v3

alevin-fry quant -t 16 -i neuron_10k_v3_quant -o neuron_10k_v3_res --tg-map transcriptome_mm10_splici_fl146/transcriptome_splici_fl146_t2g_3col.tsv --resolution cr-like --use-mtx



#neuron_1k_v3

./salmon-1.6.0_linux_x86_64/bin/salmon alevin -i mm10_splici_idx -p 20 -l ISR --chromiumV3  --sketch \
-1 /home/ggorin/datasets/neuron_1k_v3_fastqs/neuron_1k_v3_S1_L001_R1_001.fastq.gz /home/ggorin/datasets/neuron_1k_v3_fastqs/neuron_1k_v3_S1_L002_R1_001.fastq.gz \
-2 /home/ggorin/datasets/neuron_1k_v3_fastqs/neuron_1k_v3_S1_L001_R2_001.fastq.gz /home/ggorin/datasets/neuron_1k_v3_fastqs/neuron_1k_v3_S1_L002_R2_001.fastq.gz \
-o neuron_1k_v3

alevin-fry generate-permit-list -d fw -k -i neuron_1k_v3 -o neuron_1k_v3_quant

alevin-fry collate -t 16 -i neuron_1k_v3_quant -r neuron_1k_v3

alevin-fry quant -t 16 -i neuron_1k_v3_quant -o neuron_1k_v3_res --tg-map transcriptome_mm10_splici_fl146/transcriptome_splici_fl146_t2g_3col.tsv --resolution cr-like --use-mtx



#pbmc_10k_v3

./salmon-1.6.0_linux_x86_64/bin/salmon alevin -i hg38_splici_idx -p 20 -l ISR --chromiumV3  --sketch \
-1 /home/ggorin/datasets/pbmc_10k_v3_fastqs/pbmc_10k_v3_S1_L001_R1_001.fastq.gz /home/ggorin/datasets/pbmc_10k_v3_fastqs/pbmc_10k_v3_S1_L002_R1_001.fastq.gz \
-2 /home/ggorin/datasets/pbmc_10k_v3_fastqs/pbmc_10k_v3_S1_L001_R2_001.fastq.gz /home/ggorin/datasets/pbmc_10k_v3_fastqs/pbmc_10k_v3_S1_L002_R2_001.fastq.gz \
-o pbmc_10k_v3

alevin-fry generate-permit-list -d fw -k -i pbmc_10k_v3 -o pbmc_10k_v3_quant

alevin-fry collate -t 16 -i pbmc_10k_v3_quant -r pbmc_10k_v3

alevin-fry quant -t 16 -i pbmc_10k_v3_quant -o pbmc_10k_v3_res --tg-map transcriptome_hg38_splici_fl146/transcriptome_splici_fl146_t2g_3col.tsv --resolution cr-like --use-mtx



#pbmc_1k_v3

./salmon-1.6.0_linux_x86_64/bin/salmon alevin -i hg38_splici_idx -p 20 -l ISR --chromiumV3  --sketch \
-1 /home/ggorin/datasets/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz /home/ggorin/datasets/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R1_001.fastq.gz \
-2 /home/ggorin/datasets/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz /home/ggorin/datasets/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L002_R2_001.fastq.gz \
-o pbmc_1k_v3

alevin-fry generate-permit-list -d fw -k -i pbmc_1k_v3 -o pbmc_1k_v3_quant

alevin-fry collate -t 16 -i pbmc_1k_v3_quant -r pbmc_1k_v3

alevin-fry quant -t 16 -i pbmc_1k_v3_quant -o pbmc_1k_v3_res --tg-map transcriptome_hg38_splici_fl146/transcriptome_splici_fl146_t2g_3col.tsv --resolution cr-like --use-mtx
