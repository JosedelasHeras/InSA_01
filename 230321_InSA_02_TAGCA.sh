cutadapt -g HIV1=^TAGCA -e 0 --overlap 5 --discard-untrimmed --report=minimal \
   -o InSA/02_TAGCA_filter/$1_R1.fq.gz \
   -p InSA/02_TAGCA_filter/$1_R2.fq.gz \
   InSA/01_5LTR_filter/$1_R1.fq.gz \
   InSA/01_5LTR_filter/$1_R2.fq.gz \
   1> InSA/02_TAGCA_filter/report.$1.txt

