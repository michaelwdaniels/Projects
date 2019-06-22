setwd("~/Desktop/Projects/Katerina Kechris")
load('data/Rdata/cerevisiae_preprocessed.RData')
load("saved_images/rssl")

source('R/includes.R')
source('R/analysis_helpers.R')
source('R/mikes_functions.R')

library(dplyr,warn.conflicts = FALSE)
library(ggplot2,warn.conflicts = FALSE)
library(AdaSampling)

feature.table <- getFeatureTable(feature='all')
dim(feature.table)
names(feature.table)
str(feature.table)
head(feature.table)
apply(feature.table[,c(2,7,9,10,13,17,20,21,22)],2,sum)/3500
#########################
# remove integer variables with less than 5% information
idx <- which(apply(feature.table[,c(2,7,9,10,13,17,20,21,22)],2,sum)/3500<0.05)
# removed from analysis due to low content -  vacuole, in_how_many_of_5_proks_blast, intron 
idx <- c(2,7,9,10,13,17,20,21,22)[idx]
names(feature.table)[idx]
feature.table <- feature.table[,-idx]
########################
names(feature.table) # 3500 19
date.table = as.data.frame(preproc$sgd288)[,c(3,8)]
date.table = date.table[order(date.table$year, decreasing = F),]
date.table = date.table[!duplicated(date.table$ensembl),]
row.names(date.table) = date.table$ensembl

labels = as.numeric(preproc$essential)

str(feature.table)

# Analysis
p <- c(seq(0.01,0.1,0.005)) # percentage training set values
l <- length(p)
n <- p*length(labels)
contam.perc <- seq(0,0.9,0.1)
ratio <- 1-sum(labels)/length(labels)
iteration <- 100
prec_ada_50 <- prec_ada_md <- prec_ada_80 <- prec_ada_95 <- matrix(0,nrow=iteration,ncol=length(n))
rec_ada_50 <- rec_ada_md <- rec_ada_80 <- rec_ada_95 <- matrix(0,nrow=iteration,ncol=length(n))

labels.pos <- which(labels==1) # collection of row numbers from labels which are positive (769)
labels.neg <- which(labels==0) # collection of row numbers from labels which are negative (2731)
number.in.training.set <- vector() # list of number of positive labels in training set
number.un <- vector() # list of number of positive labels in training set

for(j in 1:iteration){
  cat('Iteration = ',j,' of ',iteration,'\n')
  if(j==41) next
for(i in 1:length(n)){
  set.seed(i*j^2)
  
  # Second determine which positive labels will be used in the training set and how many
  number.in.training.set[i] <- floor(p[i]*length(labels)) # (example n=70) based on p (training %) how many positive labels will be placed in the training set
  idx <- sample(1:length(labels),number.in.training.set[i]) # idx for 70 genes for training set in the 3500 genes
  labels.train <- labels[idx] # labels for training set (70)
  idx.train.neg <- which(labels.train==0) # id (within labels.train) among the 70 rows in the training set how many are 0
  idx.un <- sample(1:length(idx.train.neg),length(idx.train.neg)*ratio) # idxx (47) randomly identifies 78% of negative labels in the 70 training instances
  labels.train[idx.train.neg[idx.un]] <- NA # assign NA to labels.train (70)
  idx.label <- idx[!is.na(labels.train)] # idx of the 3500 in the training set with labels
  idx.unlabel <- idx[is.na(labels.train)] # idx of the 3500 in the training set without labels

#     cat('Training Size = ', number.in.training.set[i],'\n')
X <- as.matrix(feature.table[idx.label,])
y <- as.factor(labels.train[!is.na(labels.train)])
X_un <- as.matrix(feature.table[idx.unlabel,])
X_test <- as.matrix(feature.table[-idx,])
X_train <- rbind(X,X_un)

#ADASampling
t_semi_ada_knn <- adaSample(Ps=rownames(X[which(y==1),]),Ns=rownames(X[which(y==0),]),
                            train.mat=X_train,test.mat=X_test,classifier='knn')
#t_semi_ada_logit <- adaSample(Ps=rownames(X[which(y==1),]),Ns=rownames(X[which(y==0),]),
#                              train.mat=X_train,test.mat=X_test,classifier='logit')

pred_cut_50 <- ifelse(t_semi_ada_knn[,1]<0.5,0,1)
pred_cut_md <- ifelse(t_semi_ada_knn[,1]<median(t_semi_ada_knn[,1]),0,1) # 0.6
pred_cut_80 <- ifelse(t_semi_ada_knn[,1]<0.80,0,1)
pred_cut_95 <- ifelse(t_semi_ada_knn[,1]<0.95,0,1)

confusion_50 <- cbind(predicted=pred_cut_50,truth=labels[-idx])
confusion_md <- cbind(predicted=pred_cut_md,truth=labels[-idx])
confusion_80 <- cbind(predicted=pred_cut_80,truth=labels[-idx])
confusion_95 <- cbind(predicted=pred_cut_95,truth=labels[-idx])

# precision measures proportion of true positives among predicted positives TP/(TP+FP)
# idx predicted positives (predicted==1)
# subset confusion matrix to just true positives
# sum up the predicted column and divide by its length
prec_ada_50[j,i] <- ifelse(length(which(confusion_50[,1]==1))>1,sum(confusion_50[which(confusion_50[,1]==1),][,2])/length(confusion_50[which(confusion_50[,1]==1),][,2]),NA)
prec_ada_md[j,i] <- ifelse(length(which(confusion_md[,1]==1))>1,sum(confusion_md[which(confusion_md[,1]==1),][,2])/length(confusion_md[which(confusion_md[,1]==1),][,2]),NA)
prec_ada_80[j,i] <- ifelse(length(which(confusion_80[,1]==1))>1,sum(confusion_80[which(confusion_80[,1]==1),][,2])/length(confusion_80[which(confusion_80[,1]==1),][,2]),NA)
prec_ada_95[j,i] <- ifelse(length(which(confusion_95[,1]==1))>1,sum(confusion_95[which(confusion_95[,1]==1),][,2])/length(confusion_95[which(confusion_95[,1]==1),][,2]),NA)

# recall measures proportion of true positives predicted TP/P
# idx true positives (truth==1)
# subset confusion matrix to just true positives
# sum up the predicted column and divide by its length
rec_ada_50[j,i] <- ifelse(length(which(confusion_50[,1]==1))>1,sum(confusion_50[which(confusion_50[,1]==1),][,2])/length(confusion_50[which(confusion_50[,1]==1),][,2]),NA)
rec_ada_md[j,i] <- ifelse(length(which(confusion_md[,1]==1))>1,sum(confusion_md[which(confusion_md[,1]==1),][,2])/length(confusion_md[which(confusion_md[,1]==1),][,2]),NA)
rec_ada_80[j,i] <- ifelse(length(which(confusion_80[,1]==1))>1,sum(confusion_80[which(confusion_80[,1]==1),][,2])/length(confusion_80[which(confusion_80[,1]==1),][,2]),NA)
rec_ada_95[j,i] <- ifelse(length(which(confusion_95[,1]==1))>1,sum(confusion_95[which(confusion_95[,1]==1),][,2])/length(confusion_95[which(confusion_95[,1]==1),][,2]),NA)

# f-measure is a geometric mean of precision and recall

}
}

semi_NMC_prec <- replace(semi_NMC_prec,semi_NMC_prec==0,NA)
#sup_NMC_prec <- replace(sup_NMC_prec,sup_NMC_prec==0,NA)
semi_LS_prec <- replace(semi_LS_prec,semi_LS_prec==0,NA)
#semi_LD_prec <- replace(semi_LD_prec,semi_LD_prec==0,NA)
#sup_SVM_prec <- replace(sup_SVM_prec,sup_SVM_prec==0,NA)
#sup_LASSO_prec <- replace(sup_LASSO_prec,sup_LASSO_prec==0,NA)
#semi_LASSO_prec <- replace(semi_LASSO_prec,semi_LASSO_prec==0,NA)

NMC_semi_prec <- apply(semi_NMC_prec,2,function(x) mean(x,na.rm=T))
#NMC_sup_prec <- apply(sup_NMC_prec,2,function(x) mean(x,na.rm=T))
LS_semi_prec <- apply(semi_LS_prec,2,function(x) mean(x,na.rm=T))
#LD_semi_prec <- apply(semi_LD_prec,2,function(x) mean(x,na.rm=T))
#SVM_sup_prec <- apply(sup_SVM_prec,2,function(x) mean(x,na.rm=T))
#LASSO_sup_prec <- apply(sup_LASSO_prec,2,function(x) mean(x,na.rm=T))
#LASSO_semi_prec <- apply(semi_LASSO_prec,2,function(x) mean(x,na.rm=T))

semi_NMC_rec <- replace(semi_NMC_rec,semi_NMC_rec==0,NA)
#sup_NMC_rec <- replace(sup_NMC_rec,sup_NMC_rec==0,NA)
semi_LS_rec <- replace(semi_LS_rec,semi_LS_rec==0,NA)
#semi_LD_rec <- replace(semi_LD_rec,semi_LD_rec==0,NA)
#sup_SVM_rec <- replace(sup_SVM_rec,sup_SVM_rec==0,NA)
#sup_LASSO_rec <- replace(sup_LASSO_rec,sup_LASSO_rec==0,NA)
#semi_LASSO_rec <- replace(semi_LASSO_rec,semi_LASSO_rec==0,NA)

NMC_semi_rec <- apply(semi_NMC_rec,2,function(x) mean(x,na.rm=T))
#NMC_sup_rec <- apply(sup_NMC_rec,2,function(x) mean(x,na.rm=T))
LS_semi_rec <- apply(semi_LS_rec,2,function(x) mean(x,na.rm=T))
#LD_semi_rec <- apply(semi_LD_rec,2,function(x) mean(x,na.rm=T))
#SVM_sup_rec <- apply(sup_SVM_rec,2,function(x) mean(x,na.rm=T))
#LASSO_sup_rec <- apply(sup_LASSO_rec,2,function(x) mean(x,na.rm=T))
#LASSO_semi_rec <- apply(semi_LASSO_rec,2,function(x) mean(x,na.rm=T))

NMC_semi_f <- 2/((1/NMC_semi_prec)+(1/NMC_semi_rec))
#NMC_sup_f <- 2/((1/NMC_sup_prec)+(1/NMC_sup_rec))
LS_semi_f <- 2/((1/LS_semi_prec)+(1/LS_semi_rec))
#LD_semi_f <- 2/((1/LD_semi_prec)+(1/LD_semi_rec))
#SVM_sup_f <- 2/((1/SVM_sup_prec)+(1/SVM_sup_rec))
#LASSO_sup_f <- 2/((1/LASSO_sup_prec)+(1/LASSO_sup_rec))
#LASSO_semi_f <- 2/((1/LASSO_semi_prec)+(1/LASSO_semi_rec))

#prec_rssl <- rbind(NMC_semi_prec,NMC_sup_prec,LS_semi_prec,LD_semi_prec,SVM_sup_prec,LASSO_sup_prec,LASSO_semi_prec)
prec_rssl <- rbind(NMC_semi_prec,LS_semi_prec)
#rec_rssl <- rbind(NMC_semi_rec,NMC_sup_rec,LS_semi_rec,LD_semi_rec,SVM_sup_rec,LASSO_sup_rec,LASSO_semi_rec)
rec_rssl <- rbind(NMC_semi_rec,LS_semi_rec)
f_rssl <- rbind(NMC_semi_f,LS_semi_f)
colnames(prec_rssl) <- colnames(rec_rssl) <- colnames(f_rssl) <- p
prec_rssl
rec_rssl
f_rssl

write.csv(prec_rssl,"results/tables/rsslPrecision.csv")
write.csv(rec_rssl,"results/tables/rsslRecall.csv")
write.csv(f_rssl,"results/tables/rsslf.csv")

#save.image("saved_images/rssl")
