setwd("~/Desktop/semisupervised-mixture-models-master")
load("saved_images/6D_boxplots_all")
# ---------------------------------------------------------------------------- #
# Script 6D
# Author: Mike Daniels
#
# Uses all features from only positive training for semi-supervised method
# compares to supervised method
# examples to classify positive and negative test examples
# ---------------------------------------------------------------------------- #

cat('SCRIPT 6D: Train + only (varying n - post 2002), test +/- genes post 2002, seq features\n')

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
row.names(feature.table) = preproc$ensembl

date.table = as.data.frame(preproc$sgd288)[,c(3,8)]
date.table = date.table[order(date.table$year, decreasing = F),]
date.table = date.table[!duplicated(date.table$ensembl),]
row.names(date.table) = date.table$ensembl

labels = as.numeric(preproc$essential)

# Perform iterative calculations usinf supervised and semisupervised
n1 <- 25  # Number of iterations
p=seq(0.1,0.9,0.1) # percentage training set values
l <- length(p)

AUC.svm <- matrix(0,ncol=n1,nrow=l)
AUC.J48 <- matrix(0,ncol=n1,nrow=l)
AUC.sup <- matrix(0,ncol=n1,nrow=l)
AUC.semi <- matrix(0,ncol=n1,nrow=l)
rownames(AUC.J48) <- rownames(AUC.svm) <- rownames(AUC.sup) <- rownames(AUC.semi) <- seq(0.1,0.9,0.1)

constant=sum(labels) #769

set.seed(1000)

for(j in 1:n1){
  cat(paste0('Iteration ',j,' of ',n1,'\n'))
  all.preds.svm = list()
  all.labels.svm = list()
  
  set.seed(1000+j)
  for(i in 1:l){
    cat(paste0('Iteration ',i,' of ',l,'\n'))
    
    inTrain <- sample(1:length(labels),size = floor(p[i]*length(labels)), replace = FALSE)
    
    training <- feature.table[inTrain,]
    testing <- feature.table[-inTrain,]
    
    svm.model = svm(training,
                    y=labels[inTrain],
                    kernel="radial")
    all.preds.svm$new = predict(svm.model,testing)
    names(all.preds.svm)[length(all.preds.svm)] = paste0('SVM (train n = ', 
                                                         p[i]*length(labels), ')')
    all.labels.svm$new = labels[-inTrain]
    names(all.labels.svm)[length(all.labels.svm)] = paste0('SVM (train n = ', 
                                                           p[i]*length(labels), ')')
    
  } 
  ROC.svm <- multiroc(all.preds.svm,all.labels.svm)
  
  for(k in 1:l){
    AUC.svm[k,j] <- ROC.svm[[k]]$auc
  }
}

# Decision Tree

for(j in 1:n1){
  cat(paste0('Iteration ',j,' of ',n1,'\n'))
  all.preds.J48 = list()
  all.labels.J48 = list()
  
  set.seed(1000+j)
  for(i in 1:l){
    cat(paste0('Iteration ',i,' of ',l,'\n'))
    
    inTrain <- sample(1:length(labels),size = floor(p[i]*length(labels)), replace = FALSE)
    
    training <- cbind(as.factor(labels[inTrain]),feature.table[inTrain,])
    testing <- feature.table[-inTrain,]
    colnames(training)[1] <- 'labels'
    J48.model = J48(labels~., data = training)
            
    all.preds.J48$new = predict(J48.model,testing)
    names(all.preds.J48)[length(all.preds.J48)] = paste0('J48 (train n = ', 
                                                         p[i]*length(labels), ')')
    all.labels.J48$new = labels[-inTrain]
    names(all.labels.J48)[length(all.labels.J48)] = paste0('J48 (train n = ', 
                                                           p[i]*length(labels), ')')
    
  } 
  ROC.J48 <- multiroc(all.preds.J48,all.labels.J48)
  
  for(k in 1:l){
    AUC.J48[k,j] <- ROC.J48[[k]]$auc
  }
}

# supervised method - lasso with L1 penalty 

for(j in 1:n1){
  cat(paste0('Iteration ',j,' of ',n1,'\n'))
  
  all.preds.sup = list()
  all.labels.sup = list()
  
  set.seed(1000+j)
  
  for(i in 1:l){
    cat(paste0('Iteration ',i,' of ',l,'\n'))
    
    inTrain <- sample(1:length(labels),size = floor(p[i]*length(labels)), replace = FALSE)
    
    training <- feature.table[inTrain,]
    testing <- feature.table[-inTrain,]
    
    cv.lasso <- cv.glmnet(as.matrix(training), as.double(labels[inTrain]), family='binomial', alpha=1, standardize=TRUE, type.measure='auc')
    
    fit <- glmnet(as.matrix(training), as.double(labels[inTrain]), family='binomial',alpha = 1,lambda=cv.lasso$lambda.min)    
    
    'L1 <- optL1(labels[inTrain],training,model="logistic",epsilon = 10)'
    'fit <- penalized(labels[inTrain],training)'
    
    # Accrue prediction values and associated labels for ROC calculations
    all.preds.sup$new = predict(fit,as.matrix(testing))
    names(all.preds.sup)[length(all.preds.sup)] = paste0('Supervised (train n = ', 
                                                         p[i]*length(labels), ')')
    all.labels.sup$new = labels[-inTrain]
    names(all.labels.sup)[length(all.labels.sup)] = paste0('Supervised (train n = ', 
                                                           p[i]*length(labels), ')')
  }
  
  ROC.sup <- multiroc(all.preds.sup,all.labels.sup)
  
  for(k in 1:l){
    AUC.sup[k,j] <- ROC.sup[[k]]$auc
  }
}  


# Semi-sup method

for(j in 1:n1){
  cat(paste0('Iteration ',j,' of ',n1,'\n'))
  
  all.preds.semi = list()
  all.labels.semi = list()
  
  set.seed(1000+j)
  
  for(i in 1:l){
    cat(paste0('Iteration ',i,' of ',l,'\n'))
    
    train.size <- floor(p[i]*constant)
    seed1 <- 1 + j
    cat('Modeling with training size =', train.size, '\n')
    divided.data = divideYeastData(feature.table, date.table, labels,
                                   train.class = 'positive',
                                   num.train.examples = train.size, seed = seed1, constant = constant, negatives = 'all')
    repeat{
      if(product(apply(divided.data$train[,c(2,3,5,7,9,10,13,17,18,20,21,22)],2,sum))!=0){
        print(paste0("Seed = ",seed1," doesn't contain a column of zeros."))
        break
      } else{ 
        seed1 <- seed1+j
        divided.data = divideYeastData(feature.table, date.table, labels,
                                       train.class = 'positive',
                                       num.train.examples = train.size, seed = seed1, constant = constant, negatives = 'all')
      }
    }   
    
    
    is.train = divided.data$is.train
    lcmix.preds = getLCMIXmodels(feature.table, is.train, labels[!is.train])$preds
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

png('results/plots/6D-ROC_allFeatures_pos_only_boxplots_post_post_2002_all.png',
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC~tog$method*tog$n,at=c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24,26,27,28,29,31,32,33,34,36,37,38,39,41,42,43,44),col=c('red','lightblue1','cadetblue3','dodgerblue2'),xlab='% Training Set', ylab='AUC',xaxt="n",
        main='Supervised vs Semi-supervised (All 22 Features)\nAUC for various training % with 25 iterations\nTest: non-training essential genes & all negatives',ylim=c(0.53,0.80))
axis(side=1,at=seq(2.5,42.5,5),labels=p*100)
legend("topleft",c("Semi-supervised","Supervised-Lasso","Supervised-SVM","Supervised-Tree"),col=c('red','lightblue1','cadetblue3','dodgerblue2'),pch=15)
dev.off()

# Y axis from range = 0 - 100
png('results/plots/6D-ROC_allFeatures_pos_only_boxplots_post_post_2002_0_100.png',
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC~tog$method*tog$n,at=c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24,26,27,28,29,31,32,33,34,36,37,38,39,41,42,43,44),col=c('red','lightblue1','cadetblue3','dodgerblue2'),xlab='% Training Set', ylab='AUC',xaxt="n",
        main='Supervised vs Semi-supervised (All 22 Features)\nAUC for various training % with 25 iterations\nTest: non-training essential genes & all negatives',ylim=c(0,1))
axis(side=1,at=seq(2.5,42.5,5),labels=p*100)
legend("topleft",c("Semi-supervised","Supervised-SVM","Supervised-Lasso","Supervised-Tree"),col=c('red','lightblue1','cadetblue3','dodgerblue2'),pch=15)
dev.off()

# For paper without title
png('results/plots/6D-ROC_allFeatures_pos_only_boxplots_post_post_2002_paper.png',
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC~tog$method*tog$n,at=c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24,26,27,28,29,31,32,33,34,36,37,38,39,41,42,43,44),col=c('red','lightblue1','cadetblue3','dodgerblue2'),xlab='% Training Set', ylab='AUC',xaxt="n",ylim=c(0.53,0.80))
axis(side=1,at=seq(2.5,42.5,5),labels=p*100)
legend("topleft",c("Semi-supervised","Supervised-SVM","Supervised-Lasso","Supervised-Tree"),col=c('red','lightblue1','cadetblue3','dodgerblue2'),pch=15)
dev.off()

save.image("saved_images/6D_boxplots_all")
