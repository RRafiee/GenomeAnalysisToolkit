##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Dr Reza Rafiee, 2017
# Research Associate, Northern Institute for Cancer Research, Newcastle University
# This script creates a text file covering all file names which exist in a folder (".bam", ".fastq", ".gz", etc.) 

setwd("~/ICGC/TargetExome/Fastq")

temp = list.files(path= getwd(), pattern="*.gz")

for (i in 1:length(temp))
 {
  #print(temp[i])
  #write(temp[i], file = "master_list_1.txt", ncolumns = 2, sep = "\t", append = TRUE)
  pos1 = regexpr('TARGET_EXOME', temp[i])
  pos2 = regexpr('targetExtract', temp[i])
  N1 <- substr(temp[i], pos1[1], pos2[1]-2)
  #print(N1)
  print(paste("@RG\tID:TP",i,"\tLB:",N1,"\tSM:",N1,"\tPL:ILLUMINA",sep=""))
 }

# Output:
# [1] "@RG\tID:TP1\tLB:TARGET_EXOME_tumor_MBRep_T71_merged\tSM:TARGET_EXOME_tumor_MBRep_T71_merged\tPL:ILLUMINA"
# [2] ...

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################