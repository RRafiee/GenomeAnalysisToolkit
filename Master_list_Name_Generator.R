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
  df <- data.frame(C1,C2,C3,C3)  # this is a case in which a fastq file includes forwared and reverse read pairs together in a same file (i.e., interleaved fastq)
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
setwd("/data/Nuked")

temp = list.files(path= getwd(), pattern="*.gz")
No_of_sample<- length(temp)/2

for (i in 1:No_of_sample)
{
  #print(temp[i])
  #i <- 1
  j <- 2*i-1
  pos11 <- regexpr('New', temp[j]) 
  pos12 <- regexpr('R', temp[j]) 
  N1 <- substr(temp[j], pos11[1], pos12[1]-1)
  
  pos13 <- regexpr('New', temp[j+1]) 
  pos14 <- regexpr('R', temp[j+1])
  N2 <- substr(temp[j+1], pos13[1], pos14[1]-1) #you can get like this: _TXX_""#
  #print(N1)
  C1 <- i
  C2 <- paste("@RG\\tID:WES",i,"\\tLB:Lib_",N1,"\\tSM:MB_cellline",N1,"\\tPL:ILLUMINA",sep="")
  C3 <- paste("../FASTQ/",temp[j],sep="")
  C4 <- paste("../FASTQ/",temp[j+1],sep="")
  df <- data.frame(C1,C2,C3,C4) 
  write.table(df, file = "master_list.txt", 
              append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

# Output:

# a text file: master_list.txt with four columns for WES analysis

# 1	@RG\tID:IR1\tLB:Lib_01\tSM:ons76_neg\tPL:ILLUMINA	../FASTQ/Newc1_TAAGGCGA_L005_R1_001.fastq.gz    ../FASTQ/Newc1_TAAGGCGA_L005_R2_001.fastq.gz
# 2	@RG\tID:IR2\tLB:Lib_02\tSM:ons76_54gy_1\tPL:ILLUMINA    ../FASTQ/Newc2_CGTACTAG_L005_R1_001.fastq.gz    ../FASTQ/Newc2_CGTACTAG_L005_R2_001.fastq.gz
# 3	@RG\tID:IR3\tLB:Lib_03\tSM:ons76_54gy_2\tPL:ILLUMINA    ../FASTQ/Newc3_AGGCAGAA_L005_R1_001.fastq.gz    ../FASTQ/Newc3_AGGCAGAA_L005_R2_001.fastq.gz
# 4	@RG\tID:IR4\tLB:Lib_04\tSM:ons76_54gy_3\tPL:ILLUMINA    ../FASTQ/Newc4_TCCTGAGC_L005_R1_001.fastq.gz    ../FASTQ/Newc4_TCCTGAGC_L005_R2_001.fastq.gz
# 5	@RG\tID:IR5\tLB:Lib_05\tSM:uw402_36gy_1\tPL:ILLUMINA    ../FASTQ/Newc5_GGACTCCT_L005_R1_001.fastq.gz    ../FASTQ/Newc5_GGACTCCT_L005_R2_001.fastq.gz
# 6	@RG\tID:IR6\tLB:Lib_06\tSM:uw402_neg\tPL:ILLUMINA	../FASTQ/Newc6_CTCTCTAC_L005_R1_001.fastq.gz    ../FASTQ/Newc6_CTCTCTAC_L005_R2_001.fastq.gz

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
