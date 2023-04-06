cutadapt -g LTR5=^TAGCAGTGGCGCCCGAACAGGGA -e 0 --overlap 10 --discard-trimmed --report=minimal \
   -o InSA/01_5LTR_filter/$1_R1.fq.gz \
   -p InSA/01_5LTR_filter/$1_R2.fq.gz \
   fastq/$1_L001_R1_001.fastq.gz \
   fastq/$1_L001_R2_001.fastq.gz \
   1> InSA/01_5LTR_filter/report.$1.txt

