setwd("~/Desktop/semisupervised-mixture-models-master")
load("saved_images/6C_boxplots_balanced_50")
# ---------------------------------------------------------------------------- #
# Script 6C seq balanced
# Author: Mike Daniels
#
# Uses sequenced-derived features from only positive training for semi-supervised method
# compares to supervised methods (SVM, Decision Tree, LASSO)
# Negatives in training set will be contaminated with various percentages of positives to reflect a real world application
# ---------------------------------------------------------------------------- #

cat('SCRIPT 6C: Train + only (varying n - post 2002), test +/- genes post 2002, seq features\n')

# Get constants and helper functions
source('R/includes.R')
source('R/analysis_helpers.R')
library(randomForest)
library(caret)
library(penalized)
library(e1071)
library(glmnet)
# Load preproc - preprocessed data from Script 0
load('data/Rdata/cerevisiae_preprocessed.RData')

# Format feature.table, date.table, and labels for divideYeastData()
feature.table = as.data.frame(preproc$features)
feature.table = feature.table[,SEQUENCE_DERIVED_FEATURES]
row.names(feature.table) = preproc$ensembl

date.table = as.data.frame(preproc$sgd288)[,c(3,8)]
date.table = date.table[order(date.table$year, decreasing = F),]
date.table = date.table[!duplicated(date.table$ensembl),]
row.names(date.table) = date.table$ensembl

labels = as.numeric(preproc$essential)

# Perform iterative calculations usinf supervised and semisupervised
n1 <- 25  # Number of iterations
p=c(seq(0.01,0.1,0.01)) # percentage training set values
l <- length(p)
n <- p*length(labels)
contam.perc <- c(0.3,0.4,0.5)

AUC.svm <- matrix(0,ncol=n1,nrow=l) # collects AUC values for SVM analysis
AUC.J48 <- matrix(0,ncol=n1,nrow=l) # collects AUC values for decision tree analysis
AUC.sup <- matrix(0,ncol=n1,nrow=l) # collects AUC values for LASSO analysis
AUC.semi <- matrix(0,ncol=n1,nrow=l) # collects AUC values for semi-supervised analysis
label.svm <- matrix(0,ncol=n1,nrow=l) # collects contamination percentages for SVM for contamination confirmation 
label.tree <- matrix(0,ncol=n1,nrow=l) # collects contamination percentages for decision tree for contamination confirmation
label.lasso <- matrix(0,ncol=n1,nrow=l) # collects contamination percentages for LASSO for contamination confirmation
# Contamination does not affect semi-supervised analysis because only positive labels are used in the training set 
rownames(AUC.J48) <- rownames(AUC.svm) <- rownames(AUC.sup) <- rownames(AUC.semi) <- p # row names are training percent levels, columns will be iterations

constant=sum(labels) #769 total number of positive labels in yeast data

## SVM
for(j in 1:n1){
  cat(paste0('Iteration ',j,' of ',n1,'\n'))
  
  # Create lists to be used in analysis
  all.preds.svm = list() # list of predictions
  all.labels.svm = list() # list of labels being predicted
  n.pos = list() # list of number of positive labels in each iteration
  number.pos = vector() # list of number of positive labels in training set
  number.neg = vector() # list of total number of negative labels in training set
  number.contam = vector() # list of number of contaminates in training set
  inTrain.neg.only.length = vector() # list of number of true negative labels in training set
  
  # Within each iteration, perform analysis for each training percent level (l is the number of training sizes being analyzed)
  for(i in 1:l){
    cat(paste0('Iteration ',i,' of ',l,'\n'))
    
    set.seed(100+i*j)
    seed <- 100+i*j
    
    # First identify row numbers for positive and negative labels
    labels.pos <- which(labels==1) # collection of row numbers from labels which are positive (769)
    labels.neg <- which(labels==0) # collection of row numbers from labels which are negative (2731)
    
    # Second determine which positive labels will be used in the training set and how many
    number.pos[i] <- floor(p[i]*length(labels.pos)) # based on p (training %) how many positive labels will be placed in the training set
    idx.pos <- sample(1:length(labels.pos), size = number.pos[i], replace = FALSE) # sample this many (number.pos) from the positive labels
    n.pos$new <- number.pos[i] # collect the number of positive labels in a list for contamination calculations
    names(n.pos)[length(n.pos)] = paste0('SVM (train n = ', p[i]*length(labels), ')') # each assignment must be labeled in order to start a new list entry
    inTrain.pos <- labels.pos[idx.pos] # row numbers of positive labels randomly selected for training set to be concatenated with negative labels for training set
    
    # Third determine number of negative labels and the subset number of positive labels to be the contamination
    number.neg[i] <- floor(p[i]*length(labels.neg)) # calculates the total number of labels to be classified as negative (will contain a percentage contamination)
    number.contam[i] <- floor(contam.perc[3]*number.neg[i]) # among the number of negatives how many will be contamination based on contamination percentages
    contam.idx <- sample(1:length(labels.pos[-idx.pos]), size = number.contam,replace = FALSE) # among the positive labels not utilized in the training set, sample the number of contaminates determined in number.contam
    contam <- labels.pos[-idx.pos][contam.idx] # collection of row numbers for contaminates taken from poistive labels not already allocated to the training set
    
    # Fourth determine true negative labels and assemble full training set
    idx.neg <- sample(1:length(labels.neg), size = number.neg[i]-number.contam[i], replace = FALSE) # sample from negative labels the difference between total negatives and contamination number
    inTrain.neg.only.length[i] <- length(idx.neg) # collect number of true negatives to calculate contamination percentages
    inTrain.neg <- labels.neg[idx.neg] # collection of row numbers for negative labels only that will be included in training set
    inTrain <- c(inTrain.pos,contam,inTrain.neg) # assemble row numbers of full training set
    training <- feature.table[inTrain,] # Training set constructed using row number ids within feature table
    
    # Training Set must not contain a column of zeros or featured variable will be treated as a constant
    repeat{
      if(product(apply(training[,1:7],2,sum))!=0) {
        print(paste0("Seed = ",seed," doesn't contain a column of zeros."))
        break
      } else{ 
        seed <- seed+j*i
        set.seed(seed)
        idx.pos <- sample(1:length(labels.pos), size = number.pos[i], replace = FALSE) # sample this many (number.pos) from the positive labels
        inTrain.pos <- labels.pos[idx.pos] # row numbers of positive labels randomly selected for training set to be concatenated with negative labels for training set
        contam.idx <- sample(1:length(labels.pos[-idx.pos]), size = number.contam,replace = FALSE) # among the positive labels not utilized in the training set, sample the number of contaminates determined in number.contam
        contam <- labels.pos[-idx.pos][contam.idx] # collection of row numbers for contaminates taken from poistive labels not already allocated to the training set
        idx.neg <- sample(1:length(labels.neg), size = number.neg[i]-number.contam[i], replace = FALSE) # sample from negative labels the difference between total negatives and contamination number
        inTrain.neg <- labels.neg[idx.neg] # collection of row numbers for negative labels only that will be included in training set
        inTrain <- c(inTrain.pos,contam,inTrain.neg) # assemble row numbers of full training set
        training <- feature.table[inTrain,] # Training set constructed using row number ids within feature table
      }
    }   
    
    # Construct test set by sampling equal numbers from remaining labels
    # The test set should contain the same number of positive labels as is in the training set
    # and the same number of negative labels as the combined number of true negatives and contaminates in the training set
    
    # Fifth collect all row numbers for remaining positives and negatives
    labels.pos.test <- labels.pos[-c(idx.pos,contam.idx)] # collection of row numbers from labels which are positive but not used in training set (includes true positives and contamination)
    labels.neg.test <- labels.neg[-idx.neg] # collection of row numbers from labels which are negative
    
    # Sixth sample from test labels with corresponding training sizes (balanced approach)
    idx.pos.test <- sample(1:length(labels.pos.test), size=number.pos[i], replace=FALSE) # sample row numbers from remaining positive labels
    idx.neg.test <- sample(1:length(labels.neg.test), size=number.neg[i], replace=FALSE) # sample row numbers from remaining negative labels
    idx.test <- c(labels.pos.test[idx.pos.test],labels.neg.test[idx.neg.test]) # assemble row numbers of full testing set
    testing <- feature.table[idx.test,] # Testing set constructed using row number ids within feature table
    
    # SVM modeling 
    svm.model = svm(training, y=labels[inTrain], kernel="radial") # radial kernel showed the best performance from initial trials
    all.preds.svm$new = predict(svm.model,testing) # predictions of svm model on test set
    names(all.preds.svm)[length(all.preds.svm)] = paste0('SVM (train n = ', p[i]*length(labels), ')') # assign names to predictions in a list
    all.labels.svm$new = labels[idx.test] # collects true labels in corresponding list to be used to produce ROC
    names(all.labels.svm)[length(all.labels.svm)] = paste0('SVM (train n = ', p[i]*length(labels), ')') # assign names to labels in a list
    
  } 
  
  ROC.svm <- multiroc(all.preds.svm,all.labels.svm) # collection of ROC calculations from prediction and label lists
  
  percent.contam = list() # list of contamination percentages
  
  for(m in 1:length(number.pos)){
    percent.contam$new <- number.contam[m]/(inTrain.neg.only.length[m]+number.contam[m]) # contamination percentages determined from number of contaminates/total negative labels
    names(percent.contam)[length(percent.contam)] = p[m] # assigns names to list entries 
  }
  
  # Find AUC and collect contamination percents
  # Each row represents a training percent
  # Each column represents an iteration
  for(k in 1:l){
    AUC.svm[k,j] <- ROC.svm[[k]]$auc # extracts AUC values from ROC and places them in a matrix to be used in boxplots
    label.svm[k,j] <- percent.contam[[k]] # collects contamination percentages in each iteration of analysis
  }
  
}
svm.contam <- apply(label.svm,1,mean) # contamination percent across each row represents average contamination at each training percent (should be close to contam.perc)

## Decision Tree

for(j in 1:n1){
  cat(paste0('Iteration ',j,' of ',n1,'\n'))
  
  # Create lists to be used in analysis
  all.preds.J48 = list() # list of predictions
  all.labels.J48 = list() # list of labels being predicted
  n.pos = list() # list of number of positive labels in each iteration
  number.pos = vector() # list of number of positive labels in training set
  number.neg = vector() # list of total number of negative labels in training set
  number.contam = vector() # list of number of contaminates in training set
  inTrain.neg.only.length = vector() # list of number of true negative labels in training set
  
  # Within each iteration, perform analysis for each training percent level (l is the number of training sizes being analyzed)
  for(i in 1:l){
    cat(paste0('Iteration ',i,' of ',l,'\n'))
    
    set.seed(100+i*j)
    seed <- 100+i*j
    
    # First identify row numbers for positive and negative labels
    labels.pos <- which(labels==1) # collection of row numbers from labels which are positive (769)
    labels.neg <- which(labels==0) # collection of row numbers from labels which are negative (2731)
    
    # Second determine which positive labels will be used in the training set and how many
    number.pos[i] <- floor(p[i]*length(labels.pos)) # based on p (training %) how many positive labels will be placed in the training set
    idx.pos <- sample(1:length(labels.pos), size = number.pos[i], replace = FALSE) # sample this many (number.pos) from the positive labels
    n.pos$new <- number.pos[i] # collect the number of positive labels in a list for contamination calculations
    names(n.pos)[length(n.pos)] = paste0('Tree (train n = ', p[i]*length(labels), ')') # each assignment must be labeled in order to start a new list entry
    inTrain.pos <- labels.pos[idx.pos] # row numbers of positive labels randomly selected for training set to be concatenated with negative labels for training set
    
    # Third determine number of negative labels and the subset number of positive labels to be the contamination
    number.neg[i] <- floor(p[i]*length(labels.neg)) # calculates the total number of labels to be classified as negative (will contain a percentage contamination)
    number.contam[i] <- floor(contam.perc[3]*number.neg[i]) # among the number of negatives how many will be contamination based on contamination percentages
    contam.idx <- sample(1:length(labels.pos[-idx.pos]), size = number.contam,replace = FALSE) # among the positive labels not utilized in the training set, sample the number of contaminates determined in number.contam
    contam <- labels.pos[-idx.pos][contam.idx] # collection of row numbers for contaminates taken from poistive labels not already allocated to the training set
    
    # Fourth determine true negative labels and assemble full training set
    idx.neg <- sample(1:length(labels.neg), size = number.neg[i]-number.contam[i], replace = FALSE) # sample from negative labels the difference between total negatives and contamination number
    inTrain.neg.only.length[i] <- length(idx.neg) # collect number of true negatives to calculate contamination percentages
    inTrain.neg <- labels.neg[idx.neg] # collection of row numbers for negative labels only that will be included in training set
    inTrain <- c(inTrain.pos,contam,inTrain.neg) # assemble row numbers of full training set
    training <- feature.table[inTrain,] # Training set constructed using row number ids within feature table
    
    # Training Set must not contain a column of zeros or featured variable will be treated as a constant
    repeat{
      if(product(apply(training[,1:7],2,sum))!=0) {
        print(paste0("Seed = ",seed," doesn't contain a column of zeros."))
        break
      } else{ 
        seed <- seed+j*i
        set.seed(seed)
        idx.pos <- sample(1:length(labels.pos), size = number.pos[i], replace = FALSE) # sample this many (number.pos) from the positive labels
        inTrain.pos <- labels.pos[idx.pos] # row numbers of positive labels randomly selected for training set to be concatenated with negative labels for training set
        contam.idx <- sample(1:length(labels.pos[-idx.pos]), size = number.contam,replace = FALSE) # among the positive labels not utilized in the training set, sample the number of contaminates determined in number.contam
        contam <- labels.pos[-idx.pos][contam.idx] # collection of row numbers for contaminates taken from poistive labels not already allocated to the training set
        idx.neg <- sample(1:length(labels.neg), size = number.neg[i]-number.contam[i], replace = FALSE) # sample from negative labels the difference between total negatives and contamination number
        inTrain.neg <- labels.neg[idx.neg] # collection of row numbers for negative labels only that will be included in training set
        inTrain <- c(inTrain.pos,contam,inTrain.neg) # assemble row numbers of full training set
        training <- feature.table[inTrain,] # Training set constructed using row number ids within feature table
      }
    }   
    
    # Construct test set by sampling equal numbers from remaining labels
    # The test set should contain the same number of positive labels as is in the training set
    # and the same number of negative labels as the combined number of true negatives and contaminates in the training set
    
    # Fifth collect all row numbers for remaining positives and negatives
    labels.pos.test <- labels.pos[-c(idx.pos,contam.idx)] # collection of row numbers from labels which are positive but not used in training set (includes true positives and contamination)
    labels.neg.test <- labels.neg[-idx.neg] # collection of row numbers from labels which are negative
    
    # Sixth sample from test labels with corresponding training sizes (balanced approach)
    idx.pos.test <- sample(1:length(labels.pos.test), size=number.pos[i], replace=FALSE) # sample row numbers from remaining positive labels
    idx.neg.test <- sample(1:length(labels.neg.test), size=number.neg[i], replace=FALSE) # sample row numbers from remaining negative labels
    idx.test <- c(labels.pos.test[idx.pos.test],labels.neg.test[idx.neg.test]) # assemble row numbers of full testing set
    testing <- feature.table[idx.test,] # Testing set constructed using row number ids within feature table
    
    training.tree <- cbind(as.factor(labels[inTrain]),training)
    colnames(training.tree)[1] <- 'labels'
    J48.model = J48(labels~., data = training.tree)
    
    all.preds.J48$new = predict(J48.model,testing) # predictions of svm model on test set
    names(all.preds.J48)[length(all.preds.J48)] = paste0('Tree (train n = ', p[i]*length(labels), ')') # assign names to predictions in a list
    all.labels.J48$new = labels[idx.test] # collects true labels in corresponding list to be used to produce ROC
    names(all.labels.J48)[length(all.labels.J48)] = paste0('Tree (train n = ', p[i]*length(labels), ')') # assign names to labels in a list
  } 
  
  ROC.J48 <- multiroc(all.preds.J48,all.labels.J48) # collection of ROC calculations from prediction and label lists
  
  percent.contam = list() # list of contamination percentages
  
  for(m in 1:length(number.pos)){
    percent.contam$new <- number.contam[m]/(inTrain.neg.only.length[m]+number.contam[m]) # contamination percentages determined from number of contaminates/total negative labels
    names(percent.contam)[length(percent.contam)] = p[m] # assigns names to list entries 
  }
  
  # Find AUC and collect contamination percents
  # Each row represents a training percent
  # Each column represents an iteration
  for(k in 1:l){
    AUC.J48[k,j] <- ROC.J48[[k]]$auc # extracts AUC values from ROC and places them in a matrix to be used in boxplots
    label.tree[k,j] <- percent.contam[[k]] # collects contamination percentages in each iteration of analysis
  }
  
}
tree.contam <- apply(label.tree,1,mean) # contamination percent across each row represents average contamination at each training percent (should be close to contam.perc)

## supervised method - lasso with L1 penalty 
for(j in 1:n1){
  cat(paste0('Iteration ',j,' of ',n1,'\n'))
  
  # Create lists to be used in analysis
  all.preds.sup = list() # list of predictions
  all.labels.sup = list() # list of labels being predicted
  n.pos = list() # list of number of positive labels in each iteration
  number.pos = vector() # list of number of positive labels in training set
  number.neg = vector() # list of total number of negative labels in training set
  number.contam = vector() # list of number of contaminates in training set
  inTrain.neg.only.length = vector() # list of number of true negative labels in training set
  
  # Within each iteration, perform analysis for each training percent level (l is the number of training sizes being analyzed)
  for(i in 1:l){
    cat(paste0('Iteration ',i,' of ',l,'\n'))
    
    set.seed(100+i*j)
    seed <- 100+i*j
    
    # First identify row numbers for positive and negative labels
    labels.pos <- which(labels==1) # collection of row numbers from labels which are positive (769)
    labels.neg <- which(labels==0) # collection of row numbers from labels which are negative (2731)
    
    # Second determine which positive labels will be used in the training set and how many
    number.pos[i] <- floor(p[i]*length(labels.pos)) # based on p (training %) how many positive labels will be placed in the training set
    idx.pos <- sample(1:length(labels.pos), size = number.pos[i], replace = FALSE) # sample this many (number.pos) from the positive labels
    n.pos$new <- number.pos[i] # collect the number of positive labels in a list for contamination calculations
    names(n.pos)[length(n.pos)] = paste0('sup (train n = ', p[i]*length(labels), ')') # each assignment must be labeled in order to start a new list entry
    inTrain.pos <- labels.pos[idx.pos] # row numbers of positive labels randomly selected for training set to be concatenated with negative labels for training set
    
    # Third determine number of negative labels and the subset number of positive labels to be the contamination
    number.neg[i] <- floor(p[i]*length(labels.neg)) # calculates the total number of labels to be classified as negative (will contain a percentage contamination)
    number.contam[i] <- floor(contam.perc[3]*number.neg[i]) # among the number of negatives how many will be contamination based on contamination percentages
    contam.idx <- sample(1:length(labels.pos[-idx.pos]), size = number.contam,replace = FALSE) # among the positive labels not utilized in the training set, sample the number of contaminates determined in number.contam
    contam <- labels.pos[-idx.pos][contam.idx] # collection of row numbers for contaminates taken from poistive labels not already allocated to the training set
    
    # Fourth determine true negative labels and assemble full training set
    idx.neg <- sample(1:length(labels.neg), size = number.neg[i]-number.contam[i], replace = FALSE) # sample from negative labels the difference between total negatives and contamination number
    inTrain.neg.only.length[i] <- length(idx.neg) # collect number of true negatives to calculate contamination percentages
    inTrain.neg <- labels.neg[idx.neg] # collection of row numbers for negative labels only that will be included in training set
    inTrain <- c(inTrain.pos,contam,inTrain.neg) # assemble row numbers of full training set
    training <- feature.table[inTrain,] # Training set constructed using row number ids within feature table
    
    # Training Set must not contain a column of zeros or featured variable will be treated as a constant
    repeat{
      if(product(apply(training[,1:7],2,sum))!=0) {
        print(paste0("Seed = ",seed," doesn't contain a column of zeros."))
        break
      } else{ 
        seed <- seed+j*i
        set.seed(seed)
        idx.pos <- sample(1:length(labels.pos), size = number.pos[i], replace = FALSE) # sample this many (number.pos) from the positive labels
        inTrain.pos <- labels.pos[idx.pos] # row numbers of positive labels randomly selected for training set to be concatenated with negative labels for training set
        contam.idx <- sample(1:length(labels.pos[-idx.pos]), size = number.contam,replace = FALSE) # among the positive labels not utilized in the training set, sample the number of contaminates determined in number.contam
        contam <- labels.pos[-idx.pos][contam.idx] # collection of row numbers for contaminates taken from poistive labels not already allocated to the training set
        idx.neg <- sample(1:length(labels.neg), size = number.neg[i]-number.contam[i], replace = FALSE) # sample from negative labels the difference between total negatives and contamination number
        inTrain.neg <- labels.neg[idx.neg] # collection of row numbers for negative labels only that will be included in training set
        inTrain <- c(inTrain.pos,contam,inTrain.neg) # assemble row numbers of full training set
        training <- feature.table[inTrain,] # Training set constructed using row number ids within feature table
      }
    }   
    
    # Construct test set by sampling equal numbers from remaining labels
    # The test set should contain the same number of positive labels as is in the training set
    # and the same number of negative labels as the combined number of true negatives and contaminates in the training set
    
    # Fifth collect all row numbers for remaining positives and negatives
    labels.pos.test <- labels.pos[-c(idx.pos,contam.idx)] # collection of row numbers from labels which are positive but not used in training set (includes true positives and contamination)
    labels.neg.test <- labels.neg[-idx.neg] # collection of row numbers from labels which are negative
    
    # Sixth sample from test labels with corresponding training sizes (balanced approach)
    idx.pos.test <- sample(1:length(labels.pos.test), size=number.pos[i], replace=FALSE) # sample row numbers from remaining positive labels
    idx.neg.test <- sample(1:length(labels.neg.test), size=number.neg[i], replace=FALSE) # sample row numbers from remaining negative labels
    idx.test <- c(labels.pos.test[idx.pos.test],labels.neg.test[idx.neg.test]) # assemble row numbers of full testing set
    testing <- feature.table[idx.test,] # Testing set constructed using row number ids within feature table
    
    # LASSO modeling 
    cv.lasso <- cv.glmnet(as.matrix(training), as.double(labels[inTrain]), family='binomial', alpha=1, standardize=TRUE, type.measure='auc')
    
    fit <- glmnet(as.matrix(training), as.double(labels[inTrain]), family='binomial',alpha = 1,lambda=cv.lasso$lambda.min)    
    
    all.preds.sup$new = predict(fit,as.matrix(testing)) # predictions of sup model on test set
    names(all.preds.sup)[length(all.preds.sup)] = paste0('LASSO (train n = ', p[i]*length(labels), ')') # assign names to predictions in a list
    all.labels.sup$new = labels[idx.test] # collects true labels in corresponding list to be used to produce ROC
    names(all.labels.sup)[length(all.labels.sup)] = paste0('LASSO (train n = ', p[i]*length(labels), ')') # assign names to labels in a list
    
  } 
  
  ROC.sup <- multiroc(all.preds.sup,all.labels.sup) # collection of ROC calculations from prediction and label lists
  
  percent.contam = list() # list of contamination percentages
  
  for(m in 1:length(number.pos)){
    percent.contam$new <- number.contam[m]/(inTrain.neg.only.length[m]+number.contam[m]) # contamination percentages determined from number of contaminates/total negative labels
    names(percent.contam)[length(percent.contam)] = p[m] # assigns names to list entries 
  }
  
  # Find AUC and collect contamination percents
  # Each row represents a training percent
  # Each column represents an iteration
  for(k in 1:l){
    AUC.sup[k,j] <- ROC.sup[[k]]$auc # extracts AUC values from ROC and places them in a matrix to be used in boxplots
    label.lasso[k,j] <- percent.contam[[k]] # collects contamination percentages in each iteration of analysis
  }
  
}
sup.contam <- apply(label.lasso,1,mean) # contamination percent across each row represents average contamination at each training percent (should be close to contam.perc)

## Semi-sup method
for(j in 1:n1){
  cat(paste0('Iteration ',j,' of ',n1,'\n'))
  
  all.preds.semi = list()
  all.labels.semi = list()
  
  set.seed(100+j)
  
  for(i in 1:l){
    cat(paste0('Iteration ',i,' of ',l,'\n'))
    
    train.size <- floor(p[i]*constant)
    seed1 <- 100 + j*i
    cat('Modeling with training size =', train.size, '\n')
    divided.data = divideYeastData(feature.table, date.table, labels,
                                   train.class = 'positive',
                                   num.train.examples = train.size, seed = seed1, constant = constant, negatives = 'match')
    repeat{
      if(product(apply(divided.data$test[,c(1,2,3,4,5,6,7)],2,sum))!=0) {
        print(paste0("Seed = ",seed1," doesn't contain a column of zeros."))
        break
      } else{ 
        seed1 <- seed1+j*i
        divided.data = divideYeastData(feature.table, date.table, labels,
                                       train.class = 'positive',
                                       num.train.examples = train.size, seed = seed1, constant = constant, negatives = 'match')
      }
    }   
    
    is.train = divided.data$is.train
    lcmix.preds = getLCMIXmodels.semi.only(feature.table, is.train, labels[!is.train])$preds
    all.preds.semi$new = lcmix.preds$semisup
    names(all.preds.semi)[length(all.preds.semi)] = paste0('Semisupervised (train n = ', 
                                                           train.size, ')')
    all.labels.semi$new = labels[!divided.data$is.train]
    names(all.labels.semi)[length(all.labels.semi)] = paste0('Semisupervised (train n = ', 
                                                             train.size, ')')
  }
  
  ROC.semi <- multiroc(all.preds.semi,all.labels.semi)
  
  for(k in 1:l){
    AUC.semi[k,j] <- ROC.semi[[k]]$auc
  }
  
}

AUC2 <- t(AUC.svm)
AUC1 <- t(AUC.semi)
AUC3 <- t(AUC.sup)
AUC4 <- t(AUC.J48)
AUC <- matrix(0,ncol=4*l,nrow=n1)
for(f in 1:ncol(AUC1)){
  AUC[,4*f-3] <- AUC1[,f]
  AUC[,4*f-2] <- AUC2[,f]
  AUC[,4*f-1] <- AUC3[,f]
  AUC[,4*f]   <- AUC4[,f]
}
AUC <- as.vector(AUC)
n <- rep(p,each=4*n1)
n <- factor(n)
method <- rep(c(rep('semi',n1),
                rep('svm',n1),
                rep('sup',n1),
                rep('tree',n1)),l)

tog <- data.frame(AUC=AUC,n=n,method=method)
box.at <- seq(1:(l*5))
remove <- seq(5,length(box.at),5)
at <- box.at[!(box.at%in%remove)]

png('results/plots/6C_seqFeatures_balanced_50.png',
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC~tog$method*tog$n,at=at,col=c('red','lightblue1','cadetblue3','dodgerblue2'),xlab='% Training Set', ylab='AUC',xaxt="n",
        main='Supervised vs Semi-supervised (14 Sequence-derived Features)\nAUC for various training % with 25 iterations\nTest: balanced (50% Contamination)',ylim=c(0.5,1))
axis(side=1,at=seq(2.5,47.5,5),labels=p*100)
legend("topleft",c("Semi-supervised","Supervised-SVM","Supervised-Lasso","Supervised-Tree"),col=c('red','lightblue1','cadetblue3','dodgerblue2'),pch=15)
dev.off()

# Y axis from range = 0 - 100
png('results/plots/6C_seqFeatures_balanced_0_100_50.png',
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC~tog$method*tog$n,at=at,col=c('red','lightblue1','cadetblue3','dodgerblue2'),xlab='% Training Set', ylab='AUC',xaxt="n",
        main='Supervised vs Semi-supervised (All 22 Features)\nAUC for various training % with 25 iterations\nTest: balanced (50% Contamination)',ylim=c(0,1))
axis(side=1,at=seq(2.5,47.5,5),labels=p*100)
legend("topleft",c("Semi-supervised","Supervised-SVM","Supervised-Lasso","Supervised-Tree"),col=c('red','lightblue1','cadetblue3','dodgerblue2'),pch=15)
dev.off()

# For paper without title
png('results/plots/6C_seqFeatures_balanced_paper_50.png',
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC~tog$method*tog$n,at=at,col=c('red','lightblue1','cadetblue3','dodgerblue2'),xlab='% Training Set', ylab='AUC',xaxt="n",ylim=c(0.50,1))
axis(side=1,at=seq(2.5,47.5,5),labels=p*100)
legend("topleft",c("Semi-supervised","Supervised-SVM","Supervised-Lasso","Supervised-Tree"),col=c('red','lightblue1','cadetblue3','dodgerblue2'),pch=15)
dev.off()

save.image("saved_images/6C_boxplots_balanced_50")
