bwa mem -M -t 8 ~/genomes/Hg38_HIV1_EGFP/indices.bwa/hg38_hiv1_egfp.fa \
  InSA/03_trim_adaptor/$1_R1_val_1.fq.gz \
  InSA/03_trim_adaptor/$1_R2_val_2.fq.gz | samtools view -Shb - > InSA/04_bwa_alignment/$1_bwa.bam

