~/cellranger-6.1.1/bin/cellranger count --id=pbmc_1k_v3_cr \
--transcriptome=/home/ggorin/ref/refdata-gex-GRCh38-2020-A \
--fastqs=/home/ggorin/datasets/pbmc_1k_v3_fastqs \
        --localcores=20 \
        --localmem=80 \
--chemistry=SC3Pv3 \
--expect-cells=1000

~/cellranger-6.1.1/bin/cellranger count --id=pbmc_10k_v3_cr \
--transcriptome=/home/ggorin/ref/refdata-gex-GRCh38-2020-A \
--fastqs=/home/ggorin/datasets/pbmc_10k_v3_fastqs \
        --localcores=20 \
        --localmem=80 \
--chemistry=SC3Pv3 \
--expect-cells=10000

~/cellranger-6.1.1/bin/cellranger count --id=heart_1k_v3_cr \
--transcriptome=/home/ggorin/ref/refdata-gex-mm10-2020-A \
--fastqs=/home/ggorin/datasets/heart_1k_v3_fastqs \
        --localcores=20 \
        --localmem=80 \
--chemistry=SC3Pv3 \
--expect-cells=1000

~/cellranger-6.1.1/bin/cellranger count --id=heart_10k_v3_cr \
--transcriptome=/home/ggorin/ref/refdata-gex-mm10-2020-A \
--fastqs=/home/ggorin/datasets/heart_10k_v3_fastqs \
        --localcores=20 \
        --localmem=80 \
--chemistry=SC3Pv3 \
--expect-cells=10000

~/cellranger-6.1.1/bin/cellranger count --id=neuron_1k_v3_cr \
--transcriptome=/home/ggorin/ref/refdata-gex-mm10-2020-A \
--fastqs=/home/ggorin/datasets/neuron_1k_v3_fastqs \
        --localcores=20 \
        --localmem=80 \
--chemistry=SC3Pv3 \
--expect-cells=1000

~/cellranger-6.1.1/bin/cellranger count --id=neuron_10k_v3_cr \
--transcriptome=/home/ggorin/ref/refdata-gex-mm10-2020-A \
--fastqs=/home/ggorin/datasets/neuron_10k_v3_fastqs \
        --localcores=20 \
        --localmem=80 \
--chemistry=SC3Pv3 \
--expect-cells=10000

~/cellranger-6.1.1/bin/cellranger count --id=brain_5k_v3_cr \
--transcriptome=/home/ggorin/ref/refdata-gex-mm10-2020-A \
--fastqs=/home/ggorin/datasets/brain_10x_5k_fastqs \
        --localcores=20 \
        --localmem=80 \
--chemistry=SC3Pv3 \
--expect-cells=5000

~/cellranger-6.1.1/bin/cellranger count --id=brain_nuc_5k_v3_cr \
--transcriptome=/home/ggorin/ref/refdata-gex-mm10-2020-A \
--fastqs=/home/ggorin/datasets/brain_nuc_10x_5k_fastqs \
        --localcores=20 \
        --localmem=80 \
--chemistry=SC3Pv3 \
--expect-cells=5000
