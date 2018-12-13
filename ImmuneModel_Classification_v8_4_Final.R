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
#RNA_Seq_DataSet_7302 <- read.csv("RNA_Seq_80PercentTrainingValidationSet_All_Classes_040718_v4_ImputedbyMice_110718.csv", row.names = 1)  # with imputation (mice) see the v4 code
RNA_Seq_DataSet <- read.csv("Final_TrainingSet_ImmuneClassifier_2009_440_Log2_FPKM_RNASeq_2020205075100_121218.csv", row.names = 1)  # 2009 samples log2 FPKM with imputation (mice) see the v4 code
RNA_Seq_DataSet_Testset <- read.csv("RNA_Seq_20PercentTestSet_All_Classes_040718_v4_ImputedbyMice_110718.csv",row.names=1)

# # convert to double
# coltemp <- colnames(RNA_Seq_DataSet_7302)
# RNA_Seq_DataSet_7302 <- t(apply(RNA_Seq_DataSet_7302,1,as.numeric))  
# colnames(RNA_Seq_DataSet_7302) <- coltemp

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

# RNA_Seq_DataSet_7302 <- log2(RNA_Seq_DataSet_7302+1)
# min(RNA_Seq_DataSet_7302,na.rm = T) # [1] -1.882376
# max(RNA_Seq_DataSet_7302, na.rm = T) # [1] 20.88508

#-------------------------------------------------------------------------------------------------
#Total.No.of.Samples <- ncol(RNA_Seq_DataSet_Testset)
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

Nonlinear.svms <- readRDS("Nonlinear_SVMs_500Bootstraps_2009_131218.rds")
sel2 <- readRDS("All_500bootstraps_indices_training_2009_131218.rds")

#summary(Nonlinear.svms)

# Select the number of samples you would like to test
amount <- 50
set.seed(1234)
sel3 <- sample(1:ncol(RNA_Seq_DataSet_Testset),amount, replace = F)
RNA_Seq_DataSet_Testset <- RNA_Seq_DataSet_Testset[,sel3]

Total.No.of.Samples <- ncol(RNA_Seq_DataSet_Testset)

# prediction
#incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 3"))
Radial_test <- mclapply(1:x,
                        mc.cores=4,
                        function(i) predict(Nonlinear.svms[[i]],
                                            newdata=t(RNA_Seq_DataSet_Testset),
                                            decision.values = T,
                                            probability = T)
)

#saveRDS(Radial_test,"Radial_test_500Bootstraps_2009_1824_131218.rds")
#Radial_test <- readRDS("Radial_test_500Bootstraps_2009_1824_131218.rds")


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
#Nonlinear_SVM_Training_Model_2009 <- svm(t(RNA_Seq_DataSet),Immune_subgroup_labels(RNA_Seq_DataSet), scale = F, tolerance = 0.00001, type = "C-classification", kernel = kernel_1,cost = Kernel_optimised_cost, gamma=Kernel_optimised_gamma, probability = T, seed=i)  
#saveRDS(Nonlinear_SVM_Training_Model_2009,"Nonlinear_SVM_Training_Model_2009_131218.RData")

Nonlinear_SVM_Training_Model <- readRDS("Nonlinear_SVM_Training_Model_2009_131218.RData")

# Prediction stage: test the trained SVM model               
test.pred <- predict(object=Nonlinear_SVM_Training_Model, newdata=t(RNA_Seq_DataSet_Testset), probability=TRUE)

prob.test <- signif(attr(test.pred, "probabilities"), digits=2)
maxProbs <- apply(prob.test,1,max)

Ref_subgroup <- Immune_subgroup_labels(RNA_Seq_DataSet_Testset)
#Ref_subgroup <- Immune_subgroup_labels(RNA_Seq_DataSet_Testset[,PartialTestSelection])
Ref_subgroup <- as.factor(Ref_subgroup$unlist.Subgroup_labels_number.)

#Ref_subgroup <- Immune_subgroup_labels(RNA_Seq_DataSet)
#Ref_subgroup <- as.factor(Ref_subgroup$unlist.Subgroup_labels_number.)

confusionMatrix(test.pred,Ref_subgroup)  # for all 20% of validation

#max.col(prob.test)

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

#CMR <- confusionMatrix(test.pred,Ref_subgroup)
#CMR$byClass

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
par(mfrow=c(1,1))
#par(mar=c(6,4,2,1) + 0.1)
par(mar=c(6,4,4,1) + 0.1)
par(cex=1.3)
par(cex.axis=1)

#heading = paste("Immune subgrouping, ","440 genes, classifier ver 3.1.0", sep = "")
heading <- paste("Immune subgroup call confidence intervals for", Total.No.of.Samples, "samples")

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

points(col=maxProbsCol[order(maxProbsWhich,maxProbs)],pch=19, maxProbs[order(maxProbsWhich,maxProbs)])

axis(2, las=2)

#legend("bottomleft", legend = c("WNT", "SHH", "Grp3", "Grp4","NC"), col=c("blue", "red", "yellow2", "darkgreen","grey"), pch=19)
#axis(2, las=2)


##################################################################################################
###"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
### End
##################################################################################################
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################