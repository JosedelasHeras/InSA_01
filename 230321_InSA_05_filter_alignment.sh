samtools view -q 10 -f 66 -F 256 InSA/04_bwa_alignment/$1_bwa.bam > InSA/05_filter_alignments_R1R2/$1_R1.sam
samtools view -q 10 -f 130 -F 256 InSA/04_bwa_alignment/$1_bwa.bam > InSA/05_filter_alignments_R1R2/$1_R2.sam


