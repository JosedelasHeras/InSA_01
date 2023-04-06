read1trim=AGATCGGAAGAGC
read2trim_HIV1=GAGATTTTCCACACTGAC

trim_galore --paired -q 20 --stringency 3 --gzip --no_report_file \
  -o InSA/03_trim_adaptor \
  -a $read1trim -a2 $read2trim_HIV1 \
  InSA/02_TAGCA_filter/$1_R1.fq.gz \
  InSA/02_TAGCA_filter/$1_R2.fq.gz \
  &>> InSA/03_trim_adaptor/trim_report.$1.txt

