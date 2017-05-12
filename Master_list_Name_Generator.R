##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Dr Reza Rafiee, 2017
# Research Associate, Northern Institute for Cancer Research, Newcastle University
# This script creates a text file covering all file names which exist in a folder (".bam", ".fastq", ".gz", etc.) 

setwd("~/ICGC/TargetExomeInputFiles/Fastq")

temp = list.files(path= getwd(), pattern="*.gz")


# if (as.character(N3) == "tumor_")
# {
#   
# }

for (i in 1:length(temp))
{
  #print(temp[i])
  #i <- 1
  pos11 <- regexpr('EGA', temp[i]) #regexpr('TARGET_EXOME', temp[i])
  pos12 <- regexpr('MBRep', temp[i]) #regexpr('targetExtract', temp[i])
  N1 <- substr(temp[i], pos11[1], pos12[1]-1)
  
  pos13 <- regexpr('MBRep', temp[i]) #regexpr('targetExtract', temp[i])
  pos14 <- regexpr('target', temp[i]) #regexpr('targetExtract', temp[i])
  N2 <- substr(temp[i], pos13[1]+5, pos14[1]-2) #you can get like this: _TXX_""#
  #print(N1)
  
  pos15 <- regexpr('EXOME', temp[i]) #regexpr('TARGET_EXOME', temp[i])
  pos16 <- regexpr('MBRep', temp[i]) #regexpr('targetExtract', temp[i])
  N3 <- substr(temp[i], pos15[1]+6, pos16[1]-2)
  
  C1 <- i
  C2 <- paste("@RG\\tID:TP",i,"\\tLB:Lib_",i,"_",N3,"\\tSM:MBRep",N2,N3,"\\tPL:ILLUMINA",sep="")
  C3 <- paste("../FASTQ/",temp[i],sep="")
  df <- data.frame(C1,C2,C3,C3)  # this is a case in which a fastq file includes forwared and reverse read pairs together in a sampe file (i.e., interleaved fastq)
  print(df)
  write.table(df, file = "master_list.txt", 
              append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

# Output:
# [130] "@RG\\tID:TP130\\tLB:Lib_130_control\\tSM:MBRep_T49_mergedcontrol\\tPL:ILLUMINA"


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
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Dr Reza Rafiee, 2017
# Research Associate, Northern Institute for Cancer Research, Newcastle University
# This script generates a txt file (Mutect pairs list) which we need when running somatic calls in GATK

setwd("~/ICGC/TargetExomeInputFiles/Fastq")

temp = list.files(path= getwd(), pattern="*.gz")

# read master_list.txt from the corresponding folder

Df_master_list_txt <- data.frame(read.table("master_list.txt",header=F, sep="\t") )

for (i in 1:nrow(Df_master_list_txt))
 {
  pos11 <- regexpr('tSM', as.character(Df_master_list_txt$V2[i]))
  pos12 <- regexpr('tPL', as.character(Df_master_list_txt$V2[i])) #regexpr('targetExtract', temp[i])
  N1 <- substr(as.character(Df_master_list_txt$V2[i]), pos11[1]+4, pos12[1]-2)
  Df_master_list_txt$V1[i] <- N1
 }

Df_master_list_txt$V2 <- Df_master_list_txt$V1 
Df_master_list_txt$V1 <- rownames(Df_master_list_txt)
Df_master_list_txt$V3 <- ""
Df_master_list_txt <- Df_master_list_txt[,-4]


Df_master_list_txt$V2 <- sort(Df_master_list_txt$V2, decreasing = FALSE)

df_MuTect_pairs <- data.frame(matrix(nrow=nrow(Df_master_list_txt)/2,ncol = 3,0))

numberofpairs <- nrow(Df_master_list_txt)/2
k <- 1
while (k <= numberofpairs)
{
  l <- (2*k)/2
  df_MuTect_pairs$X1[l] <- (2*k)/2
  df_MuTect_pairs$X2[l] <- Df_master_list_txt$V2[2*k-1]
  df_MuTect_pairs$X3[l] <- Df_master_list_txt$V2[2*k]
  k <- k + 1
}

write.table(df_MuTect_pairs, file = "MuTect_pairs.txt", 
            append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Remove (or replace) everything before or after a specified character in R strings
# > x <- 'aabb.ccdd'
# > sub('.*', '', x)
# [1] ""
# > sub('bb.*', '', x)
# [1] "aa"
# > sub('.*bb', '', x)
# [1] ".ccdd"
# > sub('\\..*', '', x)
# [1] "aabb"
# > sub('.*\\.', '', x)
# [1] "ccdd"

# Output:

# 1	MBRep_T10control	MBRep_T10tumor
# 2	MBRep_T11control	MBRep_T11tumor
# 3	MBRep_T12control	MBRep_T12tumor
# 4 ...
# ...

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

