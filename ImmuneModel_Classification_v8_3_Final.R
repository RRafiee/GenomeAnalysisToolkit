##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Dr Reza Rafiee, June 2018, Version 1.0.0
# Research Fellow, Queen's University Belfast
# Immunogenomics Bioinformatics Platform:
# 1) Generating Immune Signature Scores
# 2) Predicting Immune Subgroup of input RNA-seq data
#
#
# 2.1) Training, validation and test dataset preparation
# Split input RNA-seq gene expression data into training and testing.
# Split the training data into training and validation (again, 80/20 is a fair split).
# Subsample random selections of your training data, train the classifier with this, 
# and record the performance on the validation set.
# Try a series of runs with different amounts of training data: randomly sample 20% of it,
# say, 1000 times and observe performance on the validation data, then do the same with 40%, 60%, 80%. 
# Train on all of your training data, then randomly sample a percentage of your validation data a number of times,
# and observe performance. 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# ------------------------------------------------------------------------------------------------
#                                       Loading libraries
#-------------------------------------------------------------------------------------------------
#install.packages("e1071")
library(shiny)
library(lattice)  # For imputation
#install.packages("mice")
library(mice)     # For mice imputation
#install.packages("NMF")
library(NMF)      # For non-negative matrix factorization technique
#install.packages("e1071")
library(e1071)    # For SVM classifier
library(parallel) # For mclapply which speeds up a process (probability estimation, ...)
#install.packages("caret") # only for precison and recall calcuation of the predicted results
library(caret)
#install.packages("PRROC")
library(PRROC)
#install.packages("dplyr")
library(dplyr)
#install.packages(c("FactoMineR", "factoextra"))
library("FactoMineR") # PCA
library("factoextra") # PCA
#install.packages("randomForest")
#library(randomForest)
#rfNews()
#-------------------------------------------------------------------------------------------------
Number_of_immune_Class <- 6
GoldCohort.Threshold <- 0.7304191
#-------------------------------------------------------------------------------------------------
#                                           Functions
#-------------------------------------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Data_Range_Into_01 <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Spilitting the dataset into a% (or Spilit_Percentage) and (100-a)% (or 100-Spilit_Percentage:Rest). Random subsample selection
index_of_samples_when_spilit <- function(Class_matrix, Spilit_Percentage, Immune_Class)
{
  randomseed <- sample(1:10000)[sample(1:10)[5]] # generating random numbers from 1:10K and then choosing one of the first to ten elements of this vector
  set.seed(randomseed) # The seed in set.seed produces random values which are unique to that seed (and will be same irrespective of the computer you run and hence ensures reproducibility)
  #set.seed(123) # The seed in set.seed produces random values which are unique to that seed (and will be same irrespective of the computer you run and hence ensures reproducibility)
  index_of_RNAseq_Immune_Class <- which(Class_matrix$V1 == Immune_Class)
  Per_spilit <- Spilit_Percentage/100 # Spilit_Percentage of samples in each class
  amount_Immune_Class <- round(Per_spilit*nrow(Class_matrix[index_of_RNAseq_Immune_Class,])) # Spilit_Percentage% of samples in C1 
  index_Spilit_Percentage_Percent_random_sampling_Immune_Class <- sample(as.integer(rownames(Class_matrix[index_of_RNAseq_Immune_Class,])), amount_Immune_Class, replace=F) # Bootstrap resampling, Spilit_Percentage% for each class
  index_Rest_Percent_random_sampling_Immune_Class <- setdiff(index_of_RNAseq_Immune_Class,index_Spilit_Percentage_Percent_random_sampling_Immune_Class) # Extracting the index of samples for the training and validation datasets (100-Spilit_Percentage)% or Rest Percent
  index_of_samples_when_spilit <- list(paste("index of ",Spilit_Percentage,"% of samples",sep = ""), index_Spilit_Percentage_Percent_random_sampling_Immune_Class, paste("index of ",100-Spilit_Percentage,"% of samples",sep = ""),index_Rest_Percent_random_sampling_Immune_Class)  
  return(index_of_samples_when_spilit)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCGA_Cancer_Immune_Subtype <- function(input_string)
{
  TCGA_Sample_id  <- substr(input_string, 1,12) # Sample Id's string starts from the index 1 in the input string
  N12 <- substr(input_string, 14, nchar(input_string)) # cancer subtype starts from the index 14 in the input string
  Nt1 <- regexpr('\\.', N12)
  Cancer_Subtype  <- substr(N12, 1, Nt1-1) 
  Immune_Subtype  <- substr(N12, Nt1+1, nchar(N12)-1) 
  Return_Matrix <- cbind(TCGA_Sample_id, Cancer_Subtype, Immune_Subtype)
  return(Return_Matrix)
  # return variable is a character matrix including three strings: "TCGA.05.4249" "LUAD" "C3"
}
# Test the function, e.g. ===> "TCGA.05.4249.LUAD.C3"  
# Test1 <- TCGA_Cancer_Immune_Subtype("TCGA.05.4249.LUAD.C3_") # ===> Test1[1]: "TCGA.05.4249" "LUAD" "C3"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TCGA_Gene_Signature <- function(input_string)
{
  # input_string: "ACTL6A_S5"
  Nt1 <- regexpr('\\_', input_string)
  Gene_Id <- substr(input_string, 1,Nt1-1) # Sample Id's string starts from the index 1 in the input string
  Signature_Id <- substr(input_string, Nt1+1, nchar(input_string))  
  Return_Matrix <- cbind(Gene_Id, Signature_Id)
  return(Return_Matrix)
}
# Test the function, e.g., ===> "ACTL6A_S5"
# Test1 <- TCGA_Gene_Signature("ACTL6A_S5") # ===> Test1[1]: "ACTL6A", Test1[2]: "S5"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Dataset_Spilitter <- function(Class_matrix, Split_Percentage, RNA_Seq_DataSet)
{
  RNA_Seq_Split_Percentage_percentDataSet_All_Classes <- as.data.frame(matrix(0))
  RNA_Seq_Rest_percentDataSet_All_Classes <- as.data.frame(matrix(0))
  for (t in 1:Number_of_immune_Class)
  {
    Immnue_Class_temp <- paste("C",t, sep = "") # e.g., C1 which is Immune subtype 1
    List_of_index_when_splitting <- index_of_samples_when_spilit(Class_matrix, Split_Percentage, Immnue_Class_temp) # Call the function for splitting
    RNA_Seq_Split_Percentage_percentDataSet_Class_temp <- RNA_Seq_DataSet[,List_of_index_when_splitting[[2]][]]
    RNA_Seq_Rest_percentDataSet_Class_temp <- RNA_Seq_DataSet[,List_of_index_when_splitting[[4]][]]
    # Combine the samples from all immune classes into one dataframe (Split_Percentage%) 
    RNA_Seq_Split_Percentage_percentDataSet_All_Classes <- cbind(RNA_Seq_Split_Percentage_percentDataSet_All_Classes,RNA_Seq_Split_Percentage_percentDataSet_Class_temp) 
    # Combine the samples from all immune classes into one dataframe (Rest%)
    RNA_Seq_Rest_percentDataSet_All_Classes <- cbind(RNA_Seq_Rest_percentDataSet_All_Classes, RNA_Seq_Rest_percentDataSet_Class_temp)
  }
  # Removing the extra column (first column)
  RNA_Seq_Split_Percentage_percentDataSet_All_Classes <- RNA_Seq_Split_Percentage_percentDataSet_All_Classes[,-1] 
  RNA_Seq_Rest_percentDataSet_All_Classes <- RNA_Seq_Rest_percentDataSet_All_Classes[,-1]
  # Removing the temporary matrices from environment (memory)
  #remove(RNA_Seq_Rest_percentDataSet_Class_temp)
  #remove(RNA_Seq_Split_Percentage_percentDataSet_Class_temp)
  #object.size(RNA_Seq_Split_Percentage_percentDataSet_All_Classes) # 5360864 bytes (~5MB)
  return(list(RNA_Seq_Split_Percentage_percentDataSet_All_Classes,RNA_Seq_Rest_percentDataSet_All_Classes))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Immune_subgroup_labels <- function(RNA_Seq_DataSet)
{
  # RNA_Seq_DataSet: feature*sample (row*column)
  Subgroup_labels <- lapply(1:ncol(RNA_Seq_DataSet), function(i)  {
    TCGA_Cancer_Immune_Subtype(colnames(RNA_Seq_DataSet)[i])})
  Subgroup_labels_number <- lapply(1:length(Subgroup_labels), function(j) {
    Subgroup_labels[[j]][3]})
  return(data.frame(unlist(Subgroup_labels_number)))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MultipleImputationModeling <- function(Input_Sample)
{
  # Columns must be features and rows must be samples
  # Show the missing data pattern
  md.pattern(Input_Sample)
  # Column should be variables (features or gene)
  # maxit: a scalar giving the number of iterations. The default is 5
  # m: Number of multiple imputations. The default is m=5
  NumofIteration <- 10 # 5
  imputated_by_mice <- mice(Input_Sample,m=30,maxit=NumofIteration,seed = 123)  # m=20
  #summary(imputated_by_mice)
  ### blue is observed, red is imputed
  #stripplot(imputated_by_mice, PLG_S5, pch = 19, xlab = "Imputation number")  #In general, we would like the imputations to be plausible, i.e., values that could have been observed if they had not been missing.
  #stripplot(imputated_by_mice, MMP3_S2, pch = 19, xlab = "Imputation number")  #In general, we would like the imputations to be plausible, i.e., values that could have been observed if they had not been missing.
  #stripplot(imputated_by_mice, RGS8_S5, pch = 19, xlab = "Imputation number")  #In general, we would like the imputations to be plausible, i.e., values that could have been observed if they had not been missing.
  #imputated_by_mice2 <- mice.mids(imputated_by_mice)
  # verification
  #identical(imputated_by_mice$imp, imputated_by_mice2$imp)
  ### density plot of head circumference per imputation
  ### blue is observed, red is imputed
  #densityplot(imputated_by_mice, ~PLG_S5|.imp)
  #densityplot(imputated_by_mice, ~MMP3_S2|.imp)
  #densityplot(imputated_by_mice, ~RGS8_S5|.imp)
  #xyplot(imputated_by_mice,...)
  #==========================================
  X <- complete(imputated_by_mice, action = "long", include = TRUE)[, -2]
  test <- as.mids(X, where = NULL, .imp = ".imp", .id = ".id")
  is.mids(test)
  Test_dat <- complete(test, action = "long", include = TRUE)
  #================================================
  #output1 <- matrix(unlist(imputated_by_mice2), ncol = 17, byrow = TRUE)
  #imputated_by_mice$data[1:17]
  #firstSample <- imputated_by_mice$data[3,1:17]
  ColumnNoPlus2MoreVaraibles <- ncol(Input_Sample) + 2
  Imputed_mat <- matrix(nrow = nrow(Input_Sample) ,ncol = ColumnNoPlus2MoreVaraibles ,0.0) # ncol = 19
  nrow(Test_dat)
  #for (impx in 0:5) # number of iterations
  for(idx in 1:nrow(Input_Sample)) # number of samples
  {
    for (j in 3:ColumnNoPlus2MoreVaraibles)   # number of columns
    {
      if (is.na(Test_dat[idx,j]))
      {
        sum.imp <- 0
        for (impx in 1:NumofIteration+1) # number of iterations: 6
        {
          sum.imp <- Test_dat[nrow(Input_Sample)*impx+idx,j] + sum.imp
        }
        Imputed_mat[idx,j] <- sum.imp/(NumofIteration+1)
      }
      else
      {
        Imputed_mat[idx,j] <- Test_dat[idx,j]
      }
    }
  }
  Imputed_mat <- Imputed_mat[,-c(1:2)]
  Imputed_mat <- t(Imputed_mat)
  colnames(Imputed_mat) <- rownames(Input_Sample)
  rownames(Imputed_mat) <- colnames(Input_Sample)
  return(Imputed_mat)
}
#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#--------------------------------- The main code starts from here --------------------------------
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
setwd("/home/reza/Documents/Live/")
#RNA_Seq_DataSet <- read.csv("PanCanAtlas_9126RNASeqSamplesWithImmuneSubtypes_440Genes_SampleIdsOrdered_SampleIdWithSubtypes_RR020718_RownamesGenesWithSignature.csv", row.names = 1)
#Number_of_samples_for_each_cancersubtype_Total_9126_RR180618 <- read.csv("Number_of_samples_for_each_cancersubtype_Total_9126_RR180618.csv", row.names = 1)
#RNA_Seq_DataSet <- read.csv("RNA_Seq_80PercentTrainingValidationSet_All_Classes_040718_v4.csv", row.names = 1) # without imputation
RNA_Seq_DataSet_7302 <- read.csv("RNA_Seq_80PercentTrainingValidationSet_All_Classes_040718_v4_ImputedbyMice_110718.csv", row.names = 1)  # with imputation (mice) see the v4 code
#RNA_Seq_DataSet <- read.csv("Final_TrainingSet_ImmuneClassifier_Log2_FPKM_RNASeq_061218.csv", row.names = 1)  # 1910 samples log2 FPKM with imputation (mice) see the v4 code
#RNA_Seq_DataSet <- read.csv("Final_TrainingSet_ImmuneClassifier_1669_440_Log2_FPKM_RNASeq_101218.csv", row.names = 1)  # 1669 samples log2 FPKM with imputation (mice) see the v4 code
#RNA_Seq_DataSet <- read.csv("Final_TrainingSet_ImmuneClassifier_1577_440_Log2_FPKM_RNASeq_101218.csv", row.names = 1)  # 1577 samples log2 FPKM with imputation (mice) see the v4 code
#RNA_Seq_DataSet <- read.csv("Final_TrainingSet_ImmuneClassifier_2024_440_Log2_FPKM_RNASeq_2020205075100_101218.csv", row.names = 1)  # 2024 samples log2 FPKM with imputation (mice) see the v4 code
#RNA_Seq_DataSet <- read.csv("Final_TrainingSet_ImmuneClassifier_2019_440_Log2_FPKM_RNASeq_2020205075100_111218.csv", row.names = 1)  # 2019 samples log2 FPKM with imputation (mice) see the v4 code
RNA_Seq_DataSet <- read.csv("Final_TrainingSet_ImmuneClassifier_2009_440_Log2_FPKM_RNASeq_2020205075100_121218.csv", row.names = 1)  # 2009 samples log2 FPKM with imputation (mice) see the v4 code
RNA_Seq_DataSet_Testset <- read.csv("RNA_Seq_20PercentTestSet_All_Classes_040718_v4_ImputedbyMice_110718.csv",row.names=1)

# convert to double
coltemp <- colnames(RNA_Seq_DataSet_7302)
RNA_Seq_DataSet_7302 <- t(apply(RNA_Seq_DataSet_7302,1,as.numeric))  
colnames(RNA_Seq_DataSet_7302) <- coltemp

# convert to double
coltemp <- colnames(RNA_Seq_DataSet)
RNA_Seq_DataSet <- t(apply(RNA_Seq_DataSet,1,as.numeric))  
colnames(RNA_Seq_DataSet) <- coltemp

# convert to double
coltemp <- colnames(RNA_Seq_DataSet_Testset)
RNA_Seq_DataSet_Testset <- t(apply(RNA_Seq_DataSet_Testset,1,as.numeric))  
colnames(RNA_Seq_DataSet_Testset) <- coltemp

# In RNA_Seq_DataSet, samples are in column and Genes are in row

min(RNA_Seq_DataSet) #[1] -0.7287634
max(RNA_Seq_DataSet) #[1] 1936580

# using log2 transform
RNA_Seq_DataSet_log2 <- log2(RNA_Seq_DataSet_Testset+1)
min(RNA_Seq_DataSet_log2,na.rm = T) # [1] -1.882376
max(RNA_Seq_DataSet_log2, na.rm = T) # [1] 20.88508
RNA_Seq_DataSet_Testset <- RNA_Seq_DataSet_log2
rm(RNA_Seq_DataSet_log2)

RNA_Seq_DataSet_7302 <- log2(RNA_Seq_DataSet_7302+1)
min(RNA_Seq_DataSet_7302,na.rm = T) # [1] -1.882376
max(RNA_Seq_DataSet_7302, na.rm = T) # [1] 20.88508


Total.No.of.Samples <- ncol(RNA_Seq_DataSet_Testset)

#-------------------------------------------------------------------------------------------------
# some further analyses 26/10/18

plot(density(RNA_Seq_DataSet[,1:ncol(RNA_Seq_DataSet)]))
plot(density(RNA_Seq_DataSet[1:440,ncol(RNA_Seq_DataSet)]))

#-------------------------------------------------------------------------------------------------

x <- 500 ## Number of iterations
# amount <- round(0.9*ncol(RNA_Seq_DataSet))
# 
# sel2<- lapply(1:x, function(i) {
#   set.seed(i)
#   sample(1:ncol(RNA_Seq_DataSet), amount, replace=F)
# })

# Optimised Cost and Gamma:
kernel_1 <- "radial"
# for n=2009 samples, 13/12/18
Kernel_optimised_cost  <- 4.4 #4.3 or 4.4 
Kernel_optimised_gamma <- 0.00048828125 #0.000475 

# # for n=7302 samples, 12/12/18
# Kernel_optimised_cost  <- 8.6 
# Kernel_optimised_gamma <- 0.0008 


# #00000000000000000000000000000000000000000000000000000000000000
# # parameter needed for kernel of type polynomial (default: 3)
# kernel_1 <- "radial"
# #kernel_1 <- "polynomial" # "radial"
# #kernel_1 <- "sigmoid"
# 
# #trainH <- t(RNA_Seq_DataSet_log10)
# trainH <- t(RNA_Seq_DataSet)   # 11th July 2018
# 
# # Call the function for extracting labels
# Groups <- Immune_subgroup_labels(RNA_Seq_DataSet)  # feature*sample (row*column)
# colnames(Groups) <- "Class"
# Groups <- as.factor(Groups$Class) # Groups variable needs to be a factor (in tune function)
# 
# # Finding the best parameters of the SVM classifier, columns and rows are features and samples respectively
# Tune_train <- tune(svm, train.x = trainH, train.y = Groups, scale = F,
#                    tolerance = 0.0001, type="C-classification", kernel = kernel_1,
#                    probability = T, seed = 123456, ranges=list(cost = 2^(-1:2), gamma=2^(-1:-2)), # ranges=list(cost = 2^(-2:5), gamma=2^(-2:5)),
#                    tunecontrol = tune.control(sampling = "cross", cross=10))#
# 
# write.csv(Tune_train$performances,"Tune_train_performance_error_dispersion_131218.csv")
# pdf("Tune_train_svm_performance_cv_131218.pdf")
# plot(Tune_train)
# dev.off()
# #00000000000000000000000000000000000000000000000000000000000000  
## This bit causes a delay
# incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 2"))
# Nonlinear.svms <- mclapply(1:x,
#                            mc.cores=4,
#                            function(i)  svm(x = t(RNA_Seq_DataSet)[sel2[[i]],],
#                                             y = Immune_subgroup_labels(RNA_Seq_DataSet)[sel2[[i]],], scale = F,  #
#                                             tolerance = 0.00001, type = "C-classification",
#                                             kernel = "radial",cost = Kernel_optimised_cost,
#                                             gamma=Kernel_optimised_gamma, probability = T,
#                                             seed=i)
# )
# saveRDS(Nonlinear.svms,"Nonlinear_SVMs_500Bootstraps_2009_131218.rds")
# saveRDS(sel2,"All_500bootstraps_indices_training_2009_131218.rds")

Nonlinear.svms <- readRDS("Nonlinear_SVMs_500Bootstraps_2009_131218.rds")
sel2 <- readRDS("All_500bootstraps_indices_training_2009_131218.rds")

#Test_rds <- readRDS("Nonlinear_SVMs_1KBootstrap.rds")
#summary(Nonlinear.svms)
## prediction
#incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 3"))
Radial_test <- mclapply(1:x,
                        mc.cores=2,
                        function(i) predict(Nonlinear.svms[[i]],
                                            newdata=t(RNA_Seq_DataSet_Testset),
                                            decision.values = T,
                                            probability = T)
)

saveRDS(Radial_test,"Radial_test_500Bootstraps_2009_1824_131218.rds")

# Error in mcfork() : 
#   unable to fork, possible reason: Cannot allocate memory

# selection of partial test set
#PartialTestSelection <- sample(1:ncol(RNA_Seq_DataSet_Testset),60,replace = F)
# Radial_test <- mclapply(1:x,
#                         mc.cores=4,
#                         function(i) predict(Nonlinear.svms[[i]],
#                                             newdata=t(RNA_Seq_DataSet_Testset)[PartialTestSelection,],
#                                             decision.values = T,
#                                             probability = T)
# )


# #saveRDS(Radial_test,"Radial_test_1KBootstrap.rds")
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 5"))
prob.test <- (lapply(1:x,
                     function(i) attr(Radial_test[[i]], "probabilities")))
####################################### Creating Pobes2 #############################################
#Number_of_immune_Class

k <- FALSE
for (j in 1:x) # the number of iterations
{
  for (i in 1:Number_of_immune_Class) # 6 subgroups
  {

    predProbTemp <-prob.test[[j]] # j iteration
    predProbTemp <- predProbTemp[,c("C1", "C2", "C3","C4","C5","C6")] # order the matrix based on the subgroup orders

  }
  if (k == FALSE) # Making defult tables
  {

    predProbabilities <- matrix(ncol = Number_of_immune_Class, nrow =nrow(predProbTemp)*x, 0.0)
    predProbabilities <- predProbTemp
    colnames(predProbabilities) <- c("C1","C2", "C3", "C4","C5","C6")

    k <- TRUE
  }
  else
  {
    #Adding other iteration probabilities to the created table in the ordered columns
    predProbabilities <- rbind(predProbabilities,predProbTemp)
  }
}

probs2 <- matrix(ncol=nrow(predProbTemp),nrow=x,0.0)
colnames(probs2) <- rownames(predProbTemp)

for (ttt in 1:nrow(predProbTemp)) # number of samples
{
  mmm <- matrix(ncol = Number_of_immune_Class, nrow =x, 0.0)
  colnames(mmm) <- c("C1","C2", "C3", "C4","C5","C6")
  gg <- 0
  for (fftt in 1:x)
  {
    gg <- gg + 1
    mmm[gg,] <- predProbabilities[ttt+nrow(predProbTemp)*(fftt-1),]
    #predProbabilities[1+3*0,1+3*1,1+3*2,1+nrow(predProbTemp)*(x-1)] # for the first sample, n=3
  }

  ProbSubgroup <- apply(mmm[,1:Number_of_immune_Class],1,max)
  probs2[,ttt] <- ProbSubgroup

}
####################################### creating pobes2 #############################################
# Training stage: using Support Vector Machine (SVM) with RBF kernel. SVM model parameters were optimised based on a grid-based search technique  
i=1234

Nonlinear_SVM_Training_Model_2009 <- svm(t(RNA_Seq_DataSet),Immune_subgroup_labels(RNA_Seq_DataSet), scale = F, tolerance = 0.00001, type = "C-classification", kernel = kernel_1,cost = Kernel_optimised_cost, gamma=Kernel_optimised_gamma, probability = T, seed=i)  
saveRDS(Nonlinear_SVM_Training_Model_2009,"Nonlinear_SVM_Training_Model_2009_131218.RData")

#Nonlinear_SVM_Training_Model_7302 <- svm(t(RNA_Seq_DataSet_7302),Immune_subgroup_labels(RNA_Seq_DataSet_7302), scale = F, tolerance = 0.00001, type = "C-classification", kernel = kernel_1,cost = Kernel_optimised_cost, gamma=Kernel_optimised_gamma, probability = T, seed=i)  
#saveRDS(Nonlinear_SVM_Training_Model_7302,"Nonlinear_SVM_Training_Model_7302_131218.RData")

#Nonlinear_SVM_Training_Model <- Nonlinear_SVM_Training_Model_7302
Nonlinear_SVM_Training_Model <- Nonlinear_SVM_Training_Model_2009
# Prediction stage:                
#test.pred <- predict(object=Nonlinear_SVM_Training_Model, newdata=t(RNA_Seq_DataSet_7302), probability=TRUE)

# test.pred <- predict(object=Nonlinear_SVM_Training_Model, newdata=t(RNA_Seq_DataSet_Testset), probability=TRUE)
# prob.test <- signif(attr(test.pred, "probabilities"), digits=2)
# maxProbs <- apply(prob.test,1,max)

# # 13/12/18
# #install.packages("fitdistrplus")
# library(fitdistrplus)
# krt <- kurtosis(maxProbs)
# descdist(maxProbs, discrete=FALSE, boot=1000)
# set.seed(2017)
# methodname <- "mme"#"mme" # "mle","qme","mge"  Moment matching estimation
# # Moment matching estimation consists in equalizing theoretical and empirical moments.
# # Estimated values of the distribution parameters are computed by a closed-form formula 
# # for the following distributions : "norm", "lnorm", "pois", "exp", "gamma", "nbinom", "geom", "beta", "unif" and "logis". 
# # Otherwise the theoretical and the empirical moments are matched numerically, by minimization of
# # the sum of squared differences between observed and theoretical moments.
# fKS <- fitdist(maxProbs, "beta", method=methodname)#, gof="KS")  #KS AD
# plot(fKS, demp = TRUE,addlegend=FALSE,col="grey",breaks = 64)
# GoldCohort.Threshold1 <- qbeta(p=0.05, shape1=fKS$estimate[1],shape2=fKS$estimate[2])  # 0.7304191, 5 percentile (0.05), 
# summary(fKS)
# fKS$estimate[1] # 5.582643 
# fKS$estimate[2] # 0.4007919

# fit_gamma <- fitdist(maxProbs, distr = "gamma", method = "mle")
# summary(fit_gamma)
# plot(fit_gamma, demp = TRUE,addlegend=FALSE,col="grey",breaks = 32)
# plot(fitdist(maxProbs, "weibull"))
# plot(fitdist(maxProbs, "gamma"))
# plot(fitdist(maxProbs, "lnorm"))
# summary(fitdist(maxProbs, "weibull"))
# summary(fitdist(maxProbs, "gamma"))
# summary(fitdist(maxProbs, "lnorm"))
# summary(fitdist(maxProbs, "beta", method=methodname, gof="KS"))
# anything less than this threshold will be rejected as the least significant samples 

# Prediction stage: test the trained SVM model               
test.pred <- predict(object=Nonlinear_SVM_Training_Model, newdata=t(RNA_Seq_DataSet_Testset), probability=TRUE)

#test.pred <- predict(object=Nonlinear_SVM_Training_Model, newdata=t(RNA_Seq_DataSet_Testset)[PartialTestSelection,], probability=TRUE)
#test.pred <- predict(object=Nonlinear_SVM_Training_Model, newdata=t(RNA_Seq_DataSet), probability=TRUE)

prob.test <- signif(attr(test.pred, "probabilities"), digits=2)
maxProbs <- apply(prob.test,1,max)

Ref_subgroup <- Immune_subgroup_labels(RNA_Seq_DataSet_Testset)
#Ref_subgroup <- Immune_subgroup_labels(RNA_Seq_DataSet_Testset[,PartialTestSelection])
Ref_subgroup <- as.factor(Ref_subgroup$unlist.Subgroup_labels_number.)

#Ref_subgroup <- Immune_subgroup_labels(RNA_Seq_DataSet)
#Ref_subgroup <- as.factor(Ref_subgroup$unlist.Subgroup_labels_number.)

confusionMatrix(test.pred,Ref_subgroup)  # for all 20% of validation

max.col(prob.test)

# 12/12/18, test 2009 samples with SVM trained with same data
# Confusion Matrix and Statistics
# 
# Reference
# Prediction  C1  C2  C3  C4  C5  C6
# C1 379   4   0   4   0   1
# C2   1 407   1   1   0   2
# C3   2   2 375   4   0   3
# C4   0   0   6 448   1   1
# C5   0   0   0   5 230   0
# C6   3   1   1   0   0 127
# 
# Overall Statistics
# 
# Accuracy : 0.9786          
# 95% CI : (0.9713, 0.9845)
# No Information Rate : 0.23            
# P-Value [Acc > NIR] : < 2.2e-16       
# 
# Kappa : 0.9737          
# Mcnemar's Test P-Value : NA              
# 
# Statistics by Class:
# 
#                      Class: C1 Class: C2 Class: C3 Class: C4 Class: C5 Class: C6
# Sensitivity             0.9844    0.9831    0.9791    0.9697    0.9957   0.94776
# Specificity             0.9945    0.9969    0.9932    0.9948    0.9972   0.99733
# Pos Pred Value          0.9768    0.9879    0.9715    0.9825    0.9787   0.96212
# Neg Pred Value          0.9963    0.9956    0.9951    0.9910    0.9994   0.99627
# Prevalence              0.1916    0.2061    0.1906    0.2300    0.1150   0.06670
# Detection Rate          0.1887    0.2026    0.1867    0.2230    0.1145   0.06322
# Detection Prevalence    0.1931    0.2051    0.1921    0.2270    0.1170   0.06570
# Balanced Accuracy       0.9894    0.9900    0.9862    0.9823    0.9964   0.97255

# 12/12/18, Prediction with the test set, n=1824 with the SVM trained with n=2009 samples
# Confusion Matrix and Statistics
# 
# Reference
# Prediction  C1  C2  C3  C4  C5  C6
# C1 410  23  20   7   0   4
# C2  24 462   4   6   0   3
# C3  18  10 406  17   0   6
# C4  14  16  39 190   6   0
# C5   0   0   0   8  71   0
# C6  17   7  10   3   0  23
# 
# Overall Statistics
# 
# Accuracy : 0.8564          
# 95% CI : (0.8394, 0.8721)
# No Information Rate : 0.284           
# P-Value [Acc > NIR] : < 2.2e-16       
# 
# Kappa : 0.813           
# Mcnemar's Test P-Value : NA              
# 
# Statistics by Class:
# 
#                      Class: C1 Class: C2 Class: C3 Class: C4 Class: C5 Class: C6
# Sensitivity             0.8489    0.8919    0.8476    0.8225   0.92208   0.63889
# Specificity             0.9597    0.9717    0.9621    0.9529   0.99542   0.97931
# Pos Pred Value          0.8836    0.9259    0.8884    0.7170   0.89873   0.38333
# Neg Pred Value          0.9463    0.9577    0.9466    0.9737   0.99656   0.99263
# Prevalence              0.2648    0.2840    0.2626    0.1266   0.04221   0.01974
# Detection Rate          0.2248    0.2533    0.2226    0.1042   0.03893   0.01261
# Detection Prevalence    0.2544    0.2736    0.2505    0.1453   0.04331   0.03289
# Balanced Accuracy       0.9043    0.9318    0.9048    0.8877   0.95875   0.80910

CMR <- confusionMatrix(test.pred,Ref_subgroup)
CMR$byClass

# Prediction of the final test set with training set n=2019, 11/12/18
# Confusion Matrix and Statistics
# 
# Reference
# Prediction   C1   C2  C3  C4  C5  C6
#         C1   403  19  22   6   0   1
#         C2   28  467   3   6   0   3
#         C3   18   10 398  15   0   5
#         C4   15   12  43 194   5   0
#         C5    0    0   0   8  72   0
#         C6   19   10  13   2   0  27
# 
# Overall Statistics
# 
# Accuracy : 0.8558          
# 95% CI : (0.8388, 0.8716)
# No Information Rate : 0.284           
# P-Value [Acc > NIR] : < 2.2e-16       
# 
# Kappa : 0.8128          
# Mcnemar's Test P-Value : NA              
# 
# Statistics by Class:
# 
#                      Class: C1 Class: C2 Class: C3 Class: C4 Class: C5 Class: C6
# Sensitivity             0.8344    0.9015    0.8309    0.8398   0.93506   0.75000
# Specificity             0.9642    0.9694    0.9643    0.9529   0.99542   0.97539
# Pos Pred Value          0.8936    0.9211    0.8924    0.7212   0.90000   0.38028
# Neg Pred Value          0.9417    0.9613    0.9412    0.9762   0.99713   0.99487
# Prevalence              0.2648    0.2840    0.2626    0.1266   0.04221   0.01974
# Detection Rate          0.2209    0.2560    0.2182    0.1064   0.03947   0.01480
# Detection Prevalence    0.2473    0.2780    0.2445    0.1475   0.04386   0.03893
# Balanced Accuracy       0.8993    0.9355    0.8976    0.8964   0.96524   0.86270


###########################################################################################################
###########################################################################################################
#maxProbsWhich <- factor(test.pred[1:nrow(prob.test)],levels=c("1", "2", "3", "4"))
#GoldCohort.Threshold <- 0
Threshold.Min <- GoldCohort.Threshold # based on the new threshold for Gold cohort (n=101), 8th of October 2015.
maxProbsWhich <- factor(test.pred[1:nrow(prob.test)],levels=c("C1","C2", "C3", "C4","C5","C6", "NC"))

for (i in 1:nrow(prob.test))
{
  if (maxProbs[i] <= Threshold.Min)
    maxProbsWhich[i] <- "NC"      
}
##
###########################################################################################################

# c("red","yellow", "green","cyan", "darkblue","purple")

maxProbsCol <- ifelse(maxProbsWhich== "C1","red",ifelse(maxProbsWhich=="C2","yellow",
                                                     ifelse(maxProbsWhich=="C3","forestgreen",
                                                            ifelse(maxProbsWhich=="C4","cyan",
                                                                   ifelse(maxProbsWhich=="C5", "dodgerblue3",
                                                                          ifelse(maxProbsWhich=="C6","darkorchid","SlateGray4")))))) # else: "NC"



# "red"
# "yellow"
# "green","forestgreen"
# "cyan"
# "darkblue"
# "purple"
# "SlateGray4"
# https://www.colorhexa.com/0072b2

maxProbsCol2 <- ifelse(maxProbsCol=="red","#FF0000", ifelse(maxProbsCol=="yellow","#FFFF00",
                                                               ifelse(maxProbsCol=="forestgreen","#00b240", 
                                                                      ifelse(maxProbsCol=="cyan","#00ffff",
                                                                             ifelse(maxProbsCol=="dodgerblue3","#5175C9",
                                                                                    ifelse(maxProbsCol=="darkorchid","#da70d6","#6C7B8B"))))))

# crgb <- col2rgb(cc <- colors())
# colnames(crgb) <- cc
# kkkkk <- t(crgb)
# rgb(0, 1, 0)

###########################################################################################################
############################################################
plot.new()
# par(mfrow=c(1,1))
# par(mar=c(6,4,4,5) + 0.1)
# par(cex.axis=0.8)

par(mfrow=c(1,1))
#par(mar=c(6,4,2,1) + 0.1)
par(mar=c(6,4,4,1) + 0.1)
par(cex=1.3)
par(cex.axis=1)

#heading = paste("Immune subgrouping, ","440 genes, classifier ver 3.1.0", sep = "")
heading <- paste("Immune subgroup call confidence intervals for", Total.No.of.Samples, "samples")

# boxplot(yaxt="n",xlab="",main=heading,ylab="Probability",probs2[,order(maxProbsCol, maxProbs)],outpch=NA,ylim=c(0,1),las=2, notch=FALSE,
#         col=maxProbsCol2[order(maxProbsCol,maxProbs)] )

boxplot(yaxt="n",xlab="",main=heading,ylab="Probability",probs2[,order(maxProbsWhich, maxProbs)],outpch=NA,ylim=c(0,1),las=2,
        col=maxProbsCol2[order(maxProbsWhich,maxProbs)] )

abline(col="grey",lty = 1, h = Threshold.Min)
# How many subgroups of each colour are we plotting
tmp <- table(maxProbsCol)
desired_col_order <-c("red","yellow", "forestgreen", "cyan", "dodgerblue3","darkorchid")
to_sort <- names(tmp)
# Re order by correct sub group col order using match on the desired_col_order vector
tmp <- tmp[to_sort[order(match(to_sort,desired_col_order))]]
# Index of where to draw the sub group deviders via cumsum
grp.sum <- cumsum(tmp)
# Add 0.5 to grp.sum for abline
grp.sum <- grp.sum + 0.5
# Index out final element of grp.sum to get rid of unwanted final abline
grp.sum <- grp.sum[1:length(grp.sum)-1]
abline(v=grp.sum)

# #maxProbs[order(maxProbsCol,maxProbs)] >=0.75
# #abline(col="red",lty=1, v=)
# abline(col="grey",lty=1, h=GoldCohort.Threshold)
# #abline(col="grey",lty=1, h=0.70)
# tmp <- table(maxProbsCol)
# 
# # Index of where to draw the sub group deviders via cumsum
# grp.sum <- cumsum(tmp)
# # Add 0.5 to grp.sum for abline
# grp.sum <- grp.sum + 0.5
# # Index out final element of grp.sum to get rid of unwanted final abline
# grp.sum <- grp.sum[1:length(grp.sum)-1]
# # Check
# #grp.sum
# abline(v=grp.sum)

#grp.sum <- cumsum(c(tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[6]))

#abline(v=c(grp.sum[1] + 0.5, grp.sum[2] +0.5, grp.sum[3]+0.5, grp.sum[4]+0.5, grp.sum[5]+0.5, grp.sum[5]+0.5))
#lines(col="black",lwd=2,maxProbs[order(maxProbsCol,maxProbs)])
#points(col=maxProbsCol[order(maxProbsCol,maxProbs)],pch=19, maxProbs[order(maxProbsCol,maxProbs)])
#points(col=maxProbsCol[order(maxProbsCol,maxProbs)],pch=19, maxProbs[order(maxProbsCol,maxProbs)] >= 0.75)

points(col=maxProbsCol[order(maxProbsWhich,maxProbs)],pch=19, maxProbs[order(maxProbsWhich,maxProbs)])

axis(2, las=2)

#legend("bottomleft", legend = c("WNT", "SHH", "Grp3", "Grp4","NC"), col=c("blue", "red", "yellow2", "darkgreen","grey"), pch=19)
#axis(2, las=2)

############################################################
# Evaluate the discordant test samples

confusionMatrix(test.pred,Ref_subgroup)  # for all 20% of validation

CMR <- confusionMatrix(test.pred,Ref_subgroup)
CMR$overall
# Accuracy         Kappa      AccuracyLower   AccuracyUpper
# 0.8558114      0.8127608      0.8388448      0.8716244


# Prediction of the final test set with training set n=2019, 11/12/18
# Confusion Matrix and Statistics
# 
# Reference
# Prediction   C1   C2  C3  C4  C5  C6
#         C1   403  19  22   6   0   1
#         C2   28  467   3   6   0   3
#         C3   18   10 398  15   0   5
#         C4   15   12  43 194   5   0
#         C5    0    0   0   8  72   0
#         C6   19   10  13   2   0  27

Probabilities_testset <- as.data.frame(prob.test)
max(Probabilities_testset[1,])

############################################################
library(randomForest)
# set.seed(71)
# data(iris)
# iris.rf <- randomForest(Species ~ ., data=iris, importance=TRUE,
#                         proximity=TRUE)
# print(iris.rf)

RNA_Seq_DataSet_t <- t(RNA_Seq_DataSet)
RNA_Seq_DataSet_t <- cbind(RNA_Seq_DataSet_t,Immune_subgroup_labels(t(RNA_Seq_DataSet_t)))

#Ref_subgroup_1 <- Immune_subgroup_labels(RNA_Seq_DataSet)
#Ref_subgroup_1 <- as.factor(Ref_subgroup_1$unlist.Subgroup_labels_number.)

colnames(RNA_Seq_DataSet_t)[441] <- "Subgroup"

# wrong feature name: is not allowed to use "-"
colnames(RNA_Seq_DataSet_t)[166] <- "HLA_DMA_S3"
colnames(RNA_Seq_DataSet_t)[167] <- "HLA_DRB1_S3"

# for the 1910 samples
set.seed(17)
TCGA_rf <- randomForest(Subgroup ~ ., data=RNA_Seq_DataSet_t, importance=TRUE, ntree=1000,
                        proximity=TRUE)

#install.packages("party")
library(party)
#cforest(Subgroup ~ ., data=RNA_Seq_DataSet_t, controls=cforest_control(mtry=2, mincriterion=0))

#install.packages("reprtree")
library(devtools)
devtools::install_github('araastat/reprtree')
library(reprtree)

reprtree:::plot.getTree(TCGA_rf)

print(TCGA_rf)
plot(TCGA_rf)
OOB_votes_0 <- predict(TCGA_rf,type="prob")
write.csv(OOB_votes_0,"Probability_of_All_Samples_1910_RF_1000trees_Dec18.csv")

# test.pred_rf <- predict(object=TCGA_rf, newdata=RNA_Seq_DataSet_t, type="response",norm.votes=TRUE)
# 
# #Ref_subgroup <- Immune_subgroup_labels(RNA_Seq_DataSet)
# Ref_subgroup <- Immune_subgroup_labels(RNA_Seq_DataSet_Testset[,PartialTestSelection])
# Ref_subgroup <- as.factor(Ref_subgroup$unlist.Subgroup_labels_number.)
# 
# 
# confusionMatrix(test.pred,Ref_subgroup)  # for all 20% of validation
# 
# CMR <- confusionMatrix(test.pred,Ref_subgroup)
# CMR$byClass


# 07/12/18
# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 20
# 
# OOB estimate of  error rate: 16.28%
# Confusion matrix:
#   C1  C2  C3  C4 C5 C6 class.error
# C1 432  23  17   1  0 10  0.10559006
# C2  36 469  11   2  0  0  0.09459459
# C3  21   1 438  12  0  7  0.08559499
# C4  32  23  26 142  6  2  0.38528139
# C5   0   0   2   6 69  0  0.10389610
# C6  19  21  33   0  0 49  0.59836066

# n=1577
# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 20
# 
# OOB estimate of  error rate: 18.01%
# Confusion matrix:
#   C1  C2  C3  C4 C5 C6 class.error
# C1 331  26  14   1  0 15   0.1447028
# C2  30 376   5   0  0  4   0.0939759
# C3  22   1 326  16  1 18   0.1510417
# C4  26  14  37 102  3  3   0.4486486
# C5   0   0   0   3 59  0   0.0483871
# C6  11  11  23   0  0 99   0.3125000

# 2024
# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 20
# 
# OOB estimate of  error rate: 18.58%
# Confusion matrix:
#   C1  C2  C3  C4  C5 C6 class.error
# C1 317  23  18  17   0 12  0.18087855
# C2  25 371   6  12   0  1  0.10602410
# C3  20   4 303  38   1 18  0.21093750
# C4  29  35  25 357  12  5  0.22894168
# C5   0   0   3  12 216  0  0.06493506
# C6  12  16  25   7   0 84  0.41666667

# n=2019
# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 20
# 
# OOB estimate of  error rate: 17.63%
# Confusion matrix:
#     C1  C2  C3  C4  C5 C6 class.error
# C1 319  29  12  15   0 10  0.17142857
# C2  31 359   8   9   0  7  0.13285024
# C3  16   2 310  42   0 13  0.19060052
# C4  25  24  27 370   9  7  0.19913420
# C5   0   0   1  11 219  0  0.05194805
# C6   9  16  24   9   0 86  0.40277778


# n= 2009, 12/12/2018

# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 20
# 
# OOB estimate of  error rate: 17.52%
# Confusion matrix:
#   C1  C2  C3  C4  C5 C6 class.error
# C1 324  30  13  12   0  6  0.15844156
# C2  30 363   6  10   0  5  0.12318841
# C3  15   2 315  44   0  7  0.17754569
# C4  33  22  26 364   9  8  0.21212121
# C5   0   0   1  11 219  0  0.05194805
# C6  11  13  28  10   0 72  0.46268657


class_center_TCGA <- classCenter(RNA_Seq_DataSet_t,RNA_Seq_DataSet_t$Subgroup, prox = TCGA_rf$proximity,nNbr = 7)

TCGA_rf$err.rate
TCGA_rf$votes
TCGA_rf$mtry


# for all 7302 samples
RNA_Seq_DataSet_t7 <- t(RNA_Seq_DataSet_7302)
RNA_Seq_DataSet_t7 <- cbind(RNA_Seq_DataSet_t7,Immune_subgroup_labels(t(RNA_Seq_DataSet_t7)))

colnames(RNA_Seq_DataSet_t7)[441] <- "Subgroup"

# wrong feature name: is not allowed to use "-"
colnames(RNA_Seq_DataSet_t7)[166] <- "HLA_DMA_S3"
colnames(RNA_Seq_DataSet_t7)[167] <- "HLA_DRB1_S3"

set.seed(17)
TCGA_rf_7 <- randomForest(Subgroup ~ ., data=RNA_Seq_DataSet_t7, importance=TRUE, ntree=1000,
                        proximity=TRUE)
print(TCGA_rf_7)
TCGA_rf_7$votes
plot(TCGA_rf_7)
legend("topright", legend=unique(TCGA_rf_7$classes), col=unique(as.numeric(seq(1:6))), pch=19)
# Type of random forest: classification
# Number of trees: 1000
# No. of variables tried at each split: 20
# 
# OOB estimate of  error rate: 14.04%
# Confusion matrix:
#     C1   C2   C3  C4  C5 C6 class.error
# C1 1739  115   72   6   0  1  0.10036213
# C2  108 1926   31   8   0  0  0.07091172
# C3  107   13 1745  50   2  1  0.09019812
# C4  120   90  115 584  17  0  0.36933045
# C5    0    0   12  16 280  0  0.09090909
# C6   43   33   63   2   0  3  0.97916667

OOB_votes <- predict(TCGA_rf_7,type="prob")
write.csv(OOB_votes,"Probability_of_All_Samples_7302_RF_1000trees_Dec18.csv")
#OOB_pred_C6 <- OOB_votes[,6]

# unsupervised case
TCGA_rf_Unsupervised <- randomForest(RNA_Seq_DataSet_t7[,-441])
MDSplot(TCGA_rf_Unsupervised, RNA_Seq_DataSet_t7$Subgroup)

variableUsedinRF <- varUsed(TCGA_rf_7,by.tree = FALSE, count = TRUE)   # Find out which predictor variables are actually used in the random forest.
varImpPlot(TCGA_rf_7,sort = TRUE,n.var = min(50,nrow(TCGA_rf_7$importance)))
varImpPlot(TCGA_rf,sort = TRUE,n.var = min(50,nrow(TCGA_rf$importance)))  # Dotchart of variable importance as measured by a Random Forest

# #install.packages("smotefamily")
# #install.packages("FNN")
# library(smotefamily)
# library(FNN)
# smoteresults <- SMOTE(RNA_Seq_DataSet_t[,-441], target = RNA_Seq_DataSet_t$Subgroup, K = 5 ,dup_size = 0)  # its just for binary class cases


RNA_Seq_DataSet_t1 <- t(RNA_Seq_DataSet_Testset)
RNA_Seq_DataSet_t1 <- cbind(RNA_Seq_DataSet_t1,Immune_subgroup_labels(t(RNA_Seq_DataSet_t1)))

#Ref_subgroup_1 <- Immune_subgroup_labels(RNA_Seq_DataSet)
#Ref_subgroup_1 <- as.factor(Ref_subgroup_1$unlist.Subgroup_labels_number.)

colnames(RNA_Seq_DataSet_t1)[441] <- "Subgroup"

# wrong feature name: is not allowed to use "-"
colnames(RNA_Seq_DataSet_t1)[166] <- "HLA_DMA_S3"
colnames(RNA_Seq_DataSet_t1)[167] <- "HLA_DRB1_S3"

# for the 1910 samples and unseen test set (n=1824)
set.seed(17)
TCGA_rf1 <- randomForest(Subgroup ~ ., data=RNA_Seq_DataSet_t1, importance=TRUE, ntree=500,
                        proximity=TRUE)

print(TCGA_rf1)
prob_rf_2019 <- TCGA_rf1$votes
plot(TCGA_rf1, col=c("red","yellow","green","cyan","blue","purple"))
importance(TCGA_rf1)
# Type of random forest: classification
# Number of trees: 500
# No. of variables tried at each split: 20
# 
# OOB estimate of  error rate: 16.17%
# Confusion matrix:
#   C1  C2  C3  C4 C5 C6 class.error
# C1 416  40  26   1  0  0  0.13871636
# C2  26 483   7   2  0  0  0.06756757
# C3  21   9 433  16  0  0  0.09603340
# C4  48  24  28 127  4  0  0.45021645
# C5   0   0   2   6 69  0  0.10389610
# C6   3  18  13   1  0  1  0.97222222

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#**************************
#return the rules of a tree
#**************************
getConds<-function(tree) {
  #store all conditions into a list
  conds<-list()
  #start by the terminal nodes and find previous conditions
  id.leafs<-which(tree$status==-1)
  j<-0
  for(i in id.leafs){
    j<-j+1
    prevConds<-prevCond(tree,i)
    conds[[j]]<-prevConds$cond
    while(prevConds$id>1){
      prevConds<-prevCond(tree,prevConds$id)
      conds[[j]]<-paste(conds[[j]]," & ",prevConds$cond)
    }
    if(prevConds$id==1){
      conds[[j]]<-paste(conds[[j]]," => ",tree$prediction[i])
    }
  }
  return(conds)
}

  #return(conds)
#}

#**************************
#find the previous conditions in the tree
#**************************
prevCond<-function(tree,i){
  if(i %in% tree$right_daughter){
    id<-which(tree$right_daughter==i)
    cond<-paste(tree$split_var[id],">",tree$split_point[id])
  }
  if(i %in% tree$left_daughter){
    id<-which(tree$left_daughter==i)
    cond<-paste(tree$split_var[id],"<",tree$split_point[id])
  }
  
  return(list(cond=cond,id=id))
}

#remove spaces in a word
collapse<-function(x){
  x<-sub(" ","_",x)
  
  return(x)
}


# data(iris)
# require(randomForest)
# mod.rf <- randomForest(Species ~ ., data=iris)
# tree<-getTree(mod.rf, k=1, labelVar=TRUE)
# #rename the name of the column
# colnames(tree)<-sapply(colnames(tree),collapse)
# rules<-getConds(tree)
# print(rules)

# perfect is working ... well done REZA!
require(randomForest)
TCGA_rf <- randomForest(Subgroup ~ ., data=RNA_Seq_DataSet_t, importance=TRUE, ntree=1000,proximity=TRUE)
TCGA_rf_tree <- getTree(TCGA_rf, k=1, labelVar=TRUE)
#rename the name of the column
colnames(TCGA_rf_tree)<-sapply(colnames(TCGA_rf_tree),collapse)
TCGA_rf_rules <- getConds(TCGA_rf_tree)
print(TCGA_rf_rules)
 
# #9999999999999999999
# ntree <- 5
# library("party")
# #cf <- cforest(Species~., data=iris,controls=cforest_control(ntree=ntree))
# cf <- cforest(Subgroup ~ ., data=RNA_Seq_DataSet_t, controls=cforest_control(ntree=ntree))
# 
# for(i in 1:ntree){
#   i <- 1
#   pt <- prettytree(cf@ensemble[[i]], names(cf@data@get("input"))) 
#   nt <- new("Random Forest BinaryTree") 
#   nt@tree <- pt 
#   nt@data <- cf@data 
#   nt@responses <- cf@responses 
#   
#   pdf(file=paste0("filex",i,".pdf"))
#   plot(nt, type="simple")
#   dev.off()
# }
# Error in getClass(Class, where = topenv(parent.frame())) : 
#   “Random Forest BinaryTree” is not a defined class
##################################################################################################
###"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
### End
##################################################################################################
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################