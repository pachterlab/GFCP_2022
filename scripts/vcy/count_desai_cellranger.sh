~/cellranger-6.1.1/bin/cellranger count --id=desai_idu_cr \
--transcriptome=/home/ggorin/ref/refdata-gex-mm10-2020-A \
--fastqs=/home/ggorin/datasets/desai_idu/idu \
        --localcores=20 \
        --localmem=80 \
--chemistry=SC3Pv2 \
--expect-cells=800 

~/cellranger-6.1.1/bin/cellranger count --id=desai_dmso_cr \
--transcriptome=/home/ggorin/ref/refdata-gex-mm10-2020-A \
--fastqs=/home/ggorin/datasets/desai_idu/dmso \
        --localcores=20 \
        --localmem=80 \
--chemistry=SC3Pv2 \
--expect-cells=800
