#!##############################################################################
#! Process sequences for InSA ##################################################
#!##############################################################################

#! generate "basefilenames2.txt", from fastq file names
# file containing 'meaningful' bits from the sample filenames, excluding R1/R2 as well
ls fastq/*fastq.gz | parallel basename | sed s/_L001_R._001.fastq.gz// | uniq > basefilenames2.txt

cat basefilenames2.txt | head
#NKCD3negN74DS1_S28
#NKCD3negN74DS2_S29
#NKCD3negN74DS3_S30
#NKCD3negWTS1_S25
#NKCD3negWTS2_S26
#NKCD3negWTS3_S27
#NKCD3posN74DS1_S34
#NKCD3posN74DS2_S35
#NKCD3posN74DS3_S36
#NKCD3posWTS1_S31

ls -l fastq | head
#NKCD3negN74DS1_S28_L001_R1_001.fastq.gz -> /datastore/homes3/jdlheras/SEQ_DATA/230201_InSA_Nawal_HIV1/NKCD3negN74DS1_S28_L001_R1_001.fastq.gz
#NKCD3negN74DS1_S28_L001_R2_001.fastq.gz -> /datastore/homes3/jdlheras/SEQ_DATA/230201_InSA_Nawal_HIV1/NKCD3negN74DS1_S28_L001_R2_001.fastq.gz
#NKCD3negN74DS2_S29_L001_R1_001.fastq.gz -> /datastore/homes3/jdlheras/SEQ_DATA/230201_InSA_Nawal_HIV1/NKCD3negN74DS2_S29_L001_R1_001.fastq.gz
#NKCD3negN74DS2_S29_L001_R2_001.fastq.gz -> /datastore/homes3/jdlheras/SEQ_DATA/230201_InSA_Nawal_HIV1/NKCD3negN74DS2_S29_L001_R2_001.fastq.gz
#NKCD3negN74DS3_S30_L001_R1_001.fastq.gz -> /datastore/homes3/jdlheras/SEQ_DATA/230201_InSA_Nawal_HIV1/NKCD3negN74DS3_S30_L001_R1_001.fastq.gz
#NKCD3negN74DS3_S30_L001_R2_001.fastq.gz -> /datastore/homes3/jdlheras/SEQ_DATA/230201_InSA_Nawal_HIV1/NKCD3negN74DS3_S30_L001_R2_001.fastq.gz
#NKCD3negWTS1_S25_L001_R1_001.fastq.gz -> /datastore/homes3/jdlheras/SEQ_DATA/230201_InSA_Nawal_HIV1/NKCD3negWTS1_S25_L001_R1_001.fastq.gz
#NKCD3negWTS1_S25_L001_R2_001.fastq.gz -> /datastore/homes3/jdlheras/SEQ_DATA/230201_InSA_Nawal_HIV1/NKCD3negWTS1_S25_L001_R2_001.fastq.gz
#NKCD3negWTS2_S26_L001_R1_001.fastq.gz -> /datastore/homes3/jdlheras/SEQ_DATA/230201_InSA_Nawal_HIV1/NKCD3negWTS2_S26_L001_R1_001.fastq.gz


#! 01) 5'LTR filtering: discard fragments where R1 starts at the 5' LTR ########
pwd
# /datastore/homes3/jdlheras/230127_InSA_Nawal

mkdir -p InSA/01_5LTR_filter

cat basefilenames2.txt | xargs -n1 -P10 sh 230321_InSA_01_5LTR.sh
# .sh contains:
#cutadapt -g LTR5=^TAGCAGTGGCGCCCGAACAGGGA -e 0 --overlap 10 --discard-trimmed --report=minimal \
#   -o InSA/01_5LTR_filter/$1_R1.fq.gz \
#   -p InSA/01_5LTR_filter/$1_R2.fq.gz \
#   fastq/$1_L001_R1_001.fastq.gz \
#   fastq/$1_L001_R2_001.fastq.gz \
#   1> InSA/01_5LTR_filter/report.$1.txt

## contents of folder = ~1.2Gb => remove later (leave maybe 1 sample's worth, so that there's an example)


#! 02) keep ^TAGCA-R1 fragments only, and remove TAGCA #########################
mkdir InSA/02_TAGCA_filter

cat basefilenames2.txt | xargs -n1 -P10 sh 230321_InSA_02_TAGCA.sh
# .sh contains:
#cutadapt -g HIV1=^TAGCA -e 0 --overlap 5 --discard-untrimmed --report=minimal \
#   -o InSA/02_TAGCA_filter/$1_R1.fq.gz \
#   -p InSA/02_TAGCA_filter/$1_R2.fq.gz \
#   InSA/01_5LTR_filter/$1_R1.fq.gz \
#   InSA/01_5LTR_filter/$1_R2.fq.gz \
#   1> InSA/02_TAGCA_filter/report.$1.txt

## contents of folder = ~1.0Gb => remove later (leave maybe 1 sample's worth, so that there's an example)


#! 03) trim sequencing adaptors and trim on quality ############################
mkdir InSA/03_trim_adaptor

cat basefilenames2.txt | xargs -n1 -P10 sh 230321_InSA_03_trimQadaptors.sh

# .sh contains:
#
#read1trim=AGATCGGAAGAGC
#read2trim_HIV1=GAGATTTTCCACACTGAC
#
#trim_galore --paired -q 20 --stringency 3 --gzip --no_report_file \
#  -o InSA/03_trim_adaptor \
#  -a $read1trim -a2 $read2trim_HIV1 \
#  InSA/02_TAGCA_filter/$1_R1.fq.gz \
#  InSA/02_TAGCA_filter/$1_R2.fq.gz \
#  &>> InSA/03_trim_adaptor/trim_report.$1.txt

## contents of folder = ~680Mb => remove later (leave maybe 1 sample's worth, so that there's an example)


#! 04) align to Hg38+HIV1 combined genomes using BWA MEM #######################
# use the HIV1 genome from the plasmids used, actually HIV1 with a capsid_EGFP fusion
mkdir InSA/04_bwa_alignment

cat basefilenames2.txt | xargs -n1 -P10 sh 230321_InSA_04_bwa_alignment.sh

# .sh contains:
#
#bwa mem -M -t 8 ~/genomes/Hg38_HIV1_EGFP/indices.bwa/hg38_hiv1_egfp.fa \
#  InSA/03_trim_adaptor/$1_R1_val_1.fq.gz \
#  InSA/03_trim_adaptor/$1_R2_val_2.fq.gz | samtools view -Shb - > InSA/04_bwa_alignment/$1_bwa.bam

## contents of folder = ~1 Gb .bam files


#! 05) filter out secondary alignments, keep proper pairs; separate into R1 and R2 sam files
mkdir InSA/05_filter_alignments_R1R2

cat basefilenames2.txt | xargs -n1 -P10 sh 230321_InSA_05_filter_alignment.sh

# .sh contains:
#
#samtools view -q 10 -f 66 -F 256 InSA/04_bwa_alignment/$1_bwa.bam > InSA/05_filter_alignments_R1R2/$1_R1.sam
#samtools view -q 10 -f 130 -F 256 InSA/04_bwa_alignment/$1_bwa.bam > InSA/05_filter_alignments_R1R2/$1_R2.sam

## contents of folder = ~3.3Gb => remove later (leave maybe 1 sample's worth, so that there's an example)



#! 06) extract sites; output .sites and .sisters for each sample ###############
mkdir InSA/06_extract_sites

# created an R script based on Anat's method which takes 4 parameters, and calculates sites.
# parameters required, in order:
#
# 1) basefilename
# 2) R1R2 sam folder (usually 05_filter_alignments_R1R2): filenames will be basefilename_R1.sam (and same for R2)
# 3) output folder name, with path
# 4) report filename
#
# the script is launched from bash:
# e.g.:
# Rscript --vanilla 06_extractsites.R "NKCD3negN74DS1_S28" "InSA/05_filter_alignments_R1R2" "InSA/06_extract_sites" "report_out.txt"

# teh arguments passed to the bash script are 2-4 above, in that order, while basefilename is passed on via xargs, like this:
#
# cat basefilenames2_short.txt | xargs -P1 -n1 sh 230321_InSA_06_extractsites.sh InSA/05_filter_alignments_R1R2 InSA/06_extract_sites out_report.txt

# script reports progress on console; it's a little messy when using multiple threads (-P) as the reporting was made thinking of one sample at a time
# but it doesn't seem worth spending time on this right now.


#.sh contains:
#
#### Script for extracting lists of integration sites (unique read 1) and sisters (unique paires of read1-read2) from aligned data.
#
## 1st: R1R2 sam folder. Filenames will be of the type "basefilename_R1.sam" (and same for R2)
#samfolder=$1
#
## 2nd: outputfolder
#outfolder=$2
#
##3rd: report filename
#reportfilename=$3
#
##4th: basefilename, passed using xargs -n1
#basename=$4
#
#### Check whether files are empty. Skip if they are, process using the R script if they aren't.
#lines1=`wc -l $samfolder/$basename"_R1.sam" | cut -d " " -f 1`
#lines2=`wc -l $samfolder/$basename"_R2.sam" | cut -d " " -f 1`
#
#if [ $lines1 -eq 0 ] || [ $lines2 -eq 0 ]
#then
#	echo -e "\n---\n"$fbname" - no reads remain\n----\n"
#else
#	Rscript --vanilla 06_extractsites.R $basename $samfolder $outfolder $reportfilename
#fi

cat basefilenames2.txt | xargs -P8 -n1 sh 230321_InSA_06_extractsites.sh InSA/05_filter_alignments_R1R2 InSA/06_extract_sites out_report_extractsites.txt


#! 07) report: add nsites and nsisters (per sample) ############################
pwd
# /datastore/homes3/jdlheras/230127_InSA_Nawal

mkdir InSA/07_clean_and_summarise
mkdir InSA/07_clean_and_summarise/bed_replicates
mkdir InSA/07_clean_and_summarise/bed_samples


#! 'targets' file prepared 230403
R
targets<-read.table("220403 - targets.txt", header=T, sep="\t", stringsAsFactors=F)
head(targets)
#  line        sample      i7_seq   i7_id i7_sample_sheet        type virus bad
#1    1   NKH6negWTS1 CTTATGGAAT  UDP0004     ATTCCATAAG    NKH6negWT    wt
#2    2   NKH6negWTS2 TAATCTCGTC  UDP0005     GACGAGATTA    NKH6negWT    wt
#3    3   NKH6negWTS3 GCGCGATGTT  UDP0006     AACATCGCGC    NKH6negWT    wt
#4    4 NKH6negN74DS1  AGAGCACTAG UDP0007     CTAGTGCTCT  NKH6negN74D  n74d
#5    5 NKH6negN74DS2 TGCCTTGATC  UDP0008      GATCAAGGCA NKH6negN74D  n74d
#6    6 NKH6negN74DS3 TTGTAACGGT  UDP0082      ACCGTTACAA NKH6negN74D  n74d
#
# 'bad' simply marks two samples that were prepared with the same i7 index, so they can't be distinguished and will be discarded


#! read report, add colnames & a new column to include number of unique integration sites
report<-read.table("InSA/06_extract_sites/out_report_extractsites.txt", sep="\t", header=F, stringsAsFactors=F)
head(report)
#                  V1     V2     V3     V4    V5    V6
#1 NKCD3negN74DS1_S28  15430  15410  15008 14479 14422
#2 NKCD3posN74DS1_S34  51651  51608  50263 48443 48420
#3 NKCD3posN74DS2_S35  92226  92226  92225 92211 92210
#4   NKCD3negWTS3_S27  60088  59923  58520 55911 55730
#5 NKCD3negN74DS2_S29  69375  69251  68082 65840 65651
#6   NKCD3negWTS2_S26 107957 107888 102403 99426 99212

colnames(report)<-c("sample","R1","R2","ext.R1","ext.R2","n.sites")

n.uis<-rep(0, dim(report)[1])
report<-cbind(report, n.uis)
head(report)
#              sample     R1     R2 ext.R1 ext.R2 n.sites n.uis
#1 NKCD3negN74DS1_S28  15430  15410  15008  14479   14422     0
#2 NKCD3posN74DS1_S34  51651  51608  50263  48443   48420     0
#3 NKCD3posN74DS2_S35  92226  92226  92225  92211   92210     0
#4   NKCD3negWTS3_S27  60088  59923  58520  55911   55730     0
#5 NKCD3negN74DS2_S29  69375  69251  68082  65840   65651     0
#6   NKCD3negWTS2_S26 107957 107888 102403  99426   99212     0

# R1/R2 -> number of reads @
# ext.R1/R2 -> number of reads @
# n.sites -> number of integration sites, full
# n.uis -> number of UNIQUE integration sites


#! loop through samples -> clean site/sister files, assemble adjusted sample reports, make sample BED files
### loop over samples based on 'targets' -> add n.uis to 'report'
# -> clean: filter out a handful of 'non-chr' hits[1], correct n.sites from that (= sum of totaldupl, after filtering non-chr hits)
# -> add n.uis from that, for each sample
# -> output cleaned .sites
#
# [1]: those come from the 5'LTR of HIV1 as far as I can see, but were not filtered out initially as the sequence contained mismatches.
#      it may be better to ignore that initial 5'LTR filter and clean up at teh end, or use a more permissive matching, but we'd have to be very careful to not allow
#      real host genome matches to be thrown out. Because we have thousands of reads, not hundreds of millions, it doesn't seem worth refining at this stage

report2<-report
files<-dir("InSA/06_extract_sites")

for (i in 1:dim(targets)[1])
  {
  sample<-targets[i,"sample"]
  #
  hit<-which(grepl(sample, report[,"sample"]))
  if (i %in% c(49)) next
  if (length(hit)!=1) stop("absent or non-unique hits: check sample names")
  ### sites
  hit2<-which(grepl(sample, files)); hit3<-which(grepl(".sites$", files))
  hit2<-intersect(hit2, hit3)
  sitesfile<-read.table(file.path("InSA/06_extract_sites", files[hit2]), header=T, sep="\t", stringsAsFactors=F)
  hit4<-which(is.na(sitesfile[,1]))
  if (length(hit4)>0) sitesfile<-sitesfile[-hit4,]
  #
  n.uis<-dim(sitesfile)[1]
  report2[hit,"n.uis"]<-n.uis
  #
  report2[hit,"n.sites"]<-sum(sitesfile[,"totaldupl"])
  #
  write.table(sitesfile, file.path("InSA/07_clean_and_summarise", files[hit2]), sep="\t", col.names=T, row.names=F, quote=F)
  #
  bed<-sitesfile[,c(1,2,2,1,4,3)]; bed[,3]<-bed[,2]+1; bed[,4]<-rep(sample, dim(bed)[1])
  #  head(bed)
  #   chr1      pos1    pos1.1      chr1.1 totalshear ort1
  #1 chr16  80995228  80995229 NKH6negWTS1        333    +
  #2 chr11 117246697 117246698 NKH6negWTS1        307    +
  #3  chr9  70295299  70295300 NKH6negWTS1        304    -
  #5  chr2 171821429 171821430 NKH6negWTS1        156    +
  #6  chr9  70295295  70295296 NKH6negWTS1        138    +
  #7 chr12   1765033   1765034 NKH6negWTS1        105    +
  #
  write.table(bed, file.path("InSA/07_clean_and_summarise/bed_replicates", sub("\\.sites$", "\\.bed", files[hit2])), sep="\t", col.names=F, row.names=F, quote=F)
  #
  ### sisters
  hit2<-which(grepl(sample, files)); hit3<-which(grepl(".sisters$", files))
  hit2<-intersect(hit2, hit3)
  sistersfile<-read.table(file.path("InSA/06_extract_sites", files[hit2]), header=T, sep="\t", stringsAsFactors=F)
  hit4<-which(is.na(sistersfile[,1]))
  if (length(hit4)>0) sistersfile<-sistersfile[-hit4,]
  #
  write.table(sistersfile, file.path("InSA/07_clean_and_summarise", files[hit2]), sep="\t", col.names=T, row.names=F, quote=F)
  }

write.table(report2, file.path("InSA/07_clean_and_summarise", "report_extractsites.txt"), sep="\t", row.names=F, col.names=T, quote=F)


#! merge replicates -> bed files, per sample
# call it same if < 5bp
# merge only if same strand

### targets2 (remove two 'bad' samples)
hit<-which(targets[,"bad"]=="x")
targets2<-targets[-hit,-8]
save(targets2, file="230404_targets2.RData")

samples<-unique(targets2[,"type"])      # 20 sample types
files<-dir("InSA/07_clean_and_summarise/bed_replicates")

for (i in 1:length(samples))
  {
  type<-samples[i]
  hit<-which(targets2[,"type"]==type)
  sampleids<-targets2[hit,"sample"]
  #
  hit2<-c()
  for (k in 1:length(sampleids))
    {
    hit2<-c(hit2, which(grepl(sampleids[k], files)))
    }
  filenames<-files[hit2]
  #
  outputfile<-paste0(type,".bed")
  command1<-paste("cat ", paste(file.path("InSA/07_clean_and_summarise/bed_replicates", filenames), collapse=" "), " > tmp", sep="")
  command2<-"sort -k1,1 -k2,2n tmp > tmp2"
  command3<-paste("bedtools merge -s -d 5 -i tmp2 > InSA/07_clean_and_summarise/bed_samples/", outputfile, sep="")
  system(command1, wait=T)
  system(command2, wait=T)
  system(command3, wait=T)
  }


#! report number of UIS per sample
samples<-unique(targets2[,"type"])      # 20 sample types
sample.uis.report<-as.data.frame(cbind(samples, rep(0, length(samples))), stringsAsFactors=F)
colnames(sample.uis.report)[2]<-"n.uis"

for (i in 1:length(samples))
  {
  sample<-samples[i]
  filename<-file.path("InSA/07_clean_and_summarise/bed_samples", paste0(sample,".bed"))
  tmp<-readLines(filename); n=length(tmp)
  hit<-which(sample.uis.report[,1]==sample)
  sample.uis.report[hit, "n.uis"]<-n
  }

write.table(sample.uis.report, "report_n.uis.txt", sep="\t", row.names=F, col.names=T, quote=F)



