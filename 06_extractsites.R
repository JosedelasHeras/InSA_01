# script for extracting and counting integration sites from sam files. Launched from "XdateX_InSA_06_extractsites.sh"
# single sample: meant to run over multiple samples using xargs

# v03_230323 ###################################################################

commandline <- commandArgs(trailingOnly = TRUE)
# retrieves parameters passed onto the script
# in order:
# 1) basefilename
# 2) R1R2 sam folder (usually 05_filter_alignments_R1R2): filenames will be basefilename_R1.sam (and same for R2)
# 3) output folder name, with path
# 4) report filename

basename<-commandline[1]
input.samfolder<-commandline[2]
output.folder<-commandline[3]
report.filename<-commandline[4]

suppressMessages(suppressWarnings(library(Biostrings)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(magrittr)))



### helper functions:
# whichbit            - finds sequence orientation from the sam FLAG field
# parse_cigar_neg     - checks CIGAR string to determine coordinates for start of Hs after HIV 3'LTR (- strand)
# parse_cigar_pos     - same, for (+) strand
# extractread         - extracts read info, using above functions

whichbit <- function(flag, bit) {
  # test imput
  if (!is.vector(flag)) {stop("input must be a vector")}
  # if (length(flag) > 1) {stop("input must be a single flag")}

  if (length(flag) == 1){
    # test flag to find read
    flaginbits <- intToBits(flag)
    res <- ifelse(flaginbits[bit] == 01, TRUE, FALSE)
  }

  # if multiplt flags are given, loop on same function
  if(length(flag) > 1){
    res <- logical(length(flag))
    for(i in 1:length(flag)){
      res[i] <- whichbit(flag[i], bit)
    }
  }
  # report res back
  return(res)
}
# this looks at the sam flag to decide if the match is on the + or - strand
# it looks ike the 5th bit is the one, as in the extractread function below we get: ort = ifelse(whichbit(X2, 5), "-", "+")


parse_cigar_neg <- function(cigar, seqlen){
  if (length(cigar) == 1) {
    # adjust seqlen based on deletions, insertions
    # inspired by script in http://www.bioinformatics.babraham.ac.uk/projects/bismark/deduplicate_bismark
    ops <- unlist(str_split(cigar, "\\d+")) # store the operations
    lens <- unlist(str_split(cigar, "\\D")) # store the length per operation
    ops <- ops[ops != ""]
    lens <- as.integer(lens[lens != ""]) # remove blanks from ends

    longcigar <- tibble(ops, lens)

    correctionfactor = seqlen +
      sum(longcigar$lens[longcigar$ops == "D"]) -
      sum(longcigar$lens[longcigar$ops == "I"]) -
      parse_cigar_pos(cigar)
  } else{
    correctionfactor <- numeric(length(cigar))
    for (i in 1:length(cigar)) {
      correctionfactor[i] <- parse_cigar_neg(cigar[i], seqlen[i])
    }
  }
  return(correctionfactor)
}


parse_cigar_pos <- function(cigar){
  # function for calculating the correction factor required for negative stranded reads.
  # NOTE: This is different from the fusion read version of this function
  # output: the correction factor (numeric) to use.
  if (length(cigar) == 1) {
    if (grepl("^\\d+M", cigar)) {
      correctionfactor <- 0
    } else if (grepl("^\\d+H\\d+M", cigar)) {
      correctionfactor <- 0
    } else if (grepl("^\\d+S\\d+M", cigar)) {
      correctionfactor <- as.numeric(gsub("^(\\d+)S\\d+M.*", "\\1", cigar))
    } else if (grepl("^\\d+I\\d+M", cigar)) {
      correctionfactor <- 0
    } else {
      stop("cigar not planned for ", cigar)
    }
  }

  if (length(cigar) > 1) {
    correctionfactor <- numeric(length(cigar))
    for (i in 1:length(cigar)) {
      correctionfactor[i] <- parse_cigar_pos(cigar[i])
    }
  }

  return(correctionfactor)
}


# corrected 'extractread', - fixed chr names and columns for AS/XS" tags (15/16)
extractread <- function(df){

  if (!require(tidyverse)) {
    stop("tidyverse is required to run this function")
  } else {

    df <- df %>%
      as_tibble() %>%
      filter(grepl("^AS", X15),
             grepl("^XS", X16)) %>%
      mutate(
        seqlen = nchar(X10),
        chr = factor(X3, levels = paste0("chr", c(1:22, "X", "Y"))),
        ort = ifelse(whichbit(X2, 5), "-", "+"),
        cigarcorfact = ifelse(ort == "-",
                              parse_cigar_neg(X6, seqlen),
                              parse_cigar_pos(X6)),
        seq = ifelse(ort == "-",
                     as.character(reverseComplement(DNAStringSet(X10))),
                     X10),
        pos = ifelse(ort == "-",
                     X4 + cigarcorfact - 1,
                     X4 - cigarcorfact),
        AS = as.numeric(gsub(".+:", "", X15)),
        XS = as.numeric(gsub(".+:", "", X16))
      ) %>%
      filter(AS >= 36) %>%
      select(coor = X1, chr, pos, ort, seq, mapq = X5, seqlen,
             cigar = X6, AS, XS)

  }

  return(df)

}
### helper functions END



#basename<-commandline[1]
#input.samfolder<-commandline[2]
#output.folder<-commandline[3]
#report.filename<-commandline[4]

filename1 <- paste0(basename, "_R1.sam")
filename2 <- paste0(basename, "_R2.sam")

outputsisters <- paste0(basename, ".sisters")
outputsites <- paste0(basename, ".sites")

## run starts here ####
cat("\n\n")
timestamp()
cat("\n\n----\nreading in files... basename:",basename,"\n")

#r1 <- read.table(file.path(input.samfolder, filename1),
#                 col.names = paste0("X", 1:16),
#                 fill = TRUE,
#                 quote = "",
#                 stringsAsFactors = FALSE,
#                 comment.char = "")
#
#r2 <- read.table(file.path(input.samfolder, filename2),
#                 col.names = paste0("X", 1:16),
#                 fill = TRUE,
#                 quote = "",
#                 stringsAsFactors = FALSE,
#                 comment.char = "")
                 
r1 <- read.table(file.path(input.samfolder, filename1),
                 #col.names = paste0("X", 1:16),
                 fill = TRUE,
                 quote = "",
                 sep="\t", header=F,
                 stringsAsFactors = FALSE,
                 comment.char = "",
                 flush=T)

r2 <- read.table(file.path(input.samfolder, filename2),
                 #col.names = paste0("X", 1:16),
                 fill = TRUE,
                 quote = "",
                 sep="\t", header=F,
                 stringsAsFactors = FALSE,
                 comment.char = "",
                 flush=T)

colnames(r1)<-paste0("X", 1:dim(r1)[2])
colnames(r2)<-paste0("X", 1:dim(r2)[2])


# Test whether there was a problem with the reading of the quality string
if (length(which(nchar(r1$X11) != nchar(r1$X10))) > 0) {stop("reads not read in correctly!\n"); quit(status = 123)}
if (length(which(nchar(r2$X11) != nchar(r2$X10))) > 0) {stop("reads not read in correctly!\n"); quit(status = 123)}
# ok

orig.r1 <- nrow(r1)
orig.r2 <- nrow(r2)

cat("\n----\nextracting sites...",basename,"\n")
#timestamp()

r1 <- extractread(r1)
r2 <- extractread(r2)

ext.r1 <- nrow(r1)
ext.r2 <- nrow(r2)

### exit if no reads remain (FALSE here)
if (nrow(r1) == 0 | nrow(r2) == 0) {
  tibble(basename, orig.r1, orig.r2, ext.r1, ext.r2, ext.joined = NA) %>%
    write_tsv(file=file.path(output.folder, report.filename), append = TRUE)
  cat("\n", paste("counts:", basename, orig.r1, orig.r2,
                  ext.r1, ext.r2, "NA", sep = "\t"), "\n")
  cat("no reads!\n")
  quit()
  } else {
  cat("\n----\nfiles extracted...\n")
  timestamp() }

joined <- inner_join(r1, r2, by = "coor", suffix = c("1", "2")) %>%
  mutate(length = abs(pos2 - pos1) + 1)

ext.joined <- nrow(joined)

cat("\n----\nsummarising files...",basename,"\n")

joined %>%
  group_by(chr1, pos1, ort1, chr2, pos2, ort2, length) %>%
  arrange(desc(seqlen1)) %>%
  summarise(totaldupl = n(), minseqlen1 = min(seqlen1), minmapq1 = min(mapq1), minmapq2 = min(mapq2),
            minAS1 = min(AS1), minAS2 = min(AS2), seqexample = dplyr::first(seq1)) %>%
  ungroup() %T>%
  write_tsv(file.path(output.folder, outputsisters)) %>%
  group_by(chr1, pos1, ort1) %>%
  arrange(desc(nchar(seqexample))) %>%
  summarise(totalshear = n(), totaldupl = sum(totaldupl), minseqlen1 = min(minseqlen1), minmapq1 = min(minmapq1),
            minmapq2 = min(minmapq2), minAS1 = min(minAS1), minAS2 = min(minAS2), seqexample = dplyr::first(seqexample)) %>%
  ungroup() %>%
  arrange(desc(totalshear), desc(totaldupl)) %>%
  write_tsv(file.path(output.folder, outputsites))

tibble(basename, orig.r1, orig.r2, ext.r1, ext.r2, ext.joined) %>%
  write_tsv(file = file.path(output.folder, report.filename), append = TRUE)

cat("\n", paste("counts:", basename, orig.r1, orig.r2, ext.r1, ext.r2, ext.joined, sep = "\t"), "\n")
cat("\n\nfiles are:", filename1, "and", filename2, "\n")
cat("outputs are:", outputsisters, "and", outputsites, "\n")
cat("-- done!\n")
timestamp()
cat("\n\n--\n")

