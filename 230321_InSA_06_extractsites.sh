### Script for extracting lists of integration sites (unique read 1) and sisters (unique paires of read1-read2) from aligned data. 

# 1st: R1R2 sam folder. Filenames will be of the type "basefilename_R1.sam" (and same for R2)
samfolder=$1

# 2nd: outputfolder
outfolder=$2

#3rd: report filename
reportfilename=$3

#4th: basefilename, passed using xargs -n1
basename=$4

### Check whether files are empty. Skip if they are, process using the R script if they aren't. 
lines1=`wc -l $samfolder/$basename"_R1.sam" | cut -d " " -f 1`
lines2=`wc -l $samfolder/$basename"_R2.sam" | cut -d " " -f 1`

if [ $lines1 -eq 0 ] || [ $lines2 -eq 0 ]
then
	echo -e "\n---\n"$fbname" - no reads remain\n----\n"
else
	Rscript --vanilla 06_extractsites.R $basename $samfolder $outfolder $reportfilename
fi
