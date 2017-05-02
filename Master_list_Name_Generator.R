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
  #i <- 3
  pos11 <- regexpr('EGA', temp[i]) 
  pos12 <- regexpr('MBRep', temp[i]) 
  N1 <- substr(temp[i], pos11[1], pos12[1]-1)
  
  pos13 <- regexpr('MBRep', temp[i]) 
  pos14 <- regexpr('merged', temp[i]) 
  N2 <- substr(temp[i], pos13[1]+5, pos14[1]-1) #you can get like this: _TXX_""#
  #print(N1)
  C1 <- i
  C2 <- paste("@RG\\tID:TP",i,"\\tLB:",N1,"\\tSM:MBRep",N2,"\\tPL:ILLUMINA",sep="")
  C3 <- paste("../FASTQ/",temp[i],sep="")
  df <- data.frame(C1,C2,C3,C3)  # this is a case in which a fastq file includes forwared and reverse read pairs together in a sampe file (i.e., interleaved fastq)
  write.table(df, file = "master_list.txt", 
              append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
 }

# Output:
# [1] "@RG\\tID:TP130\\tLB:EGAR00001031596_TARGET_EXOME_control_\\tSM:MBRep_T49_\\tPL:ILLUMINA"
# [2] ...

# Output:
# a text file: master_list.txt with three columns for WES analysis

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
