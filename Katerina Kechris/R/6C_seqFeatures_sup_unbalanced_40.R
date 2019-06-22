setwd("~/Desktop/semisupervised-mixture-models-master")
load("saved_images/6C_boxplots_unbalanced")
# ---------------------------------------------------------------------------- #
# Script 6C seq unbalanced
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
source('R/mikes_functions.R')
library(randomForest)
library(caret)
library(penalized)
library(e1071)
library(glmnet)
# Load preproc - preprocessed data from Script 0
load('data/Rdata/cerevisiae_preprocessed.RData')
feature.table <- getFeatureTable(feature='seq')
dim(feature.table)
names(feature.table)
str(feature.table)
head(feature.table)
apply(feature.table[,c(1,2,3,4,5,6,7)],2,sum)/3500
#########################
# remove integer variables with less than 20% information
idx <- which(apply(feature.table[,c(1,2,3,4,5,6,7)],2,sum)/3500<0.2)
idx <- c(c(1,2,3,4,5,6,7))[idx]
names(feature.table)[idx]
feature.table <- feature.table[,-idx]

# Format feature.table, date.table, and labels for divideYeastData()
row.names(feature.table) = preproc$ensembl
date.table = as.data.frame(preproc$sgd288)[,c(3,8)]
date.table = date.table[order(date.table$year, decreasing = F),]
date.table = date.table[!duplicated(date.table$ensembl),]
row.names(date.table) = date.table$ensembl

labels = as.numeric(preproc$essential)

# Perform iterative calculations using supervised and semisupervised
n1 <- 25  # Number of iterations
p <- c(seq(0.01,0.1,0.005)) # percentage training set values
l <- length(p)
n <- p*length(labels)
contam.perc <- seq(0,0.9,0.1)

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
AUC.list <- list()
svm.pred <- list()
svm.perf <- list()


for(g in 1:length(contam.perc)){
  
  for(j in 1:n1){
    cat(paste0('Iteration ',j,' of ',n1,'\n'))
    
    # Create lists to be used in analysis
    all.preds.svm = list() # list of predictions
    all.labels.svm = list() # list of labels being predicted
    all.preds.J48 = list() # list of predictions
    all.labels.J48 = list() # list of labels being predicted
    all.preds.sup = list() # list of predictions
    all.labels.sup = list() # list of labels being predicted
    if(g<2) {
      all.preds.semi = list()
      all.labels.semi = list()
      all.perfs.semi = list()}
    n.pos = list() # list of number of positive labels in each iteration
    number.pos = vector() # list of number of positive labels in training set
    number.neg = vector() # list of total number of negative labels in training set
    number.contam = vector() # list of number of contaminates in training set
    inTrain.neg.only.length = vector() # list of number of true negative labels in training set
    
    # Within each iteration, perform analysis for each training percent level (l is the number of training sizes being analyzed)
    for(i in 1:l){
      cat(paste0('Iteration ',i,' of ',l,'\n'))
      
      seed <- 21+g*i*j
      seed1 <- seed
      set.seed(seed)
      # First identify row numbers for positive and negative labels
      labels.pos <- which(labels==1) # collection of row numbers from labels which are positive (769)
      labels.neg <- which(labels==0) # collection of row numbers from labels which are negative (2731)
      
      # Second determine which positive labels will be used in the training set and how many
      number.pos[i] <- floor(p[i]*length(labels.pos)) # based on p (training %) how many positive labels will be placed in the training set
      idx.pos <- sample(1:length(labels.pos), size = number.pos[i], replace = FALSE) # sample this many (number.pos) from the positive labels
      n.pos$new <- number.pos[i] # collect the number of positive labels in a list for contamination calculations
      names(n.pos)[length(n.pos)] = paste0('train n = ', p[i]*length(labels)) # each assignment must be labeled in order to start a new list entry
      inTrain.pos <- labels.pos[idx.pos] # row numbers of positive labels randomly selected for training set to be concatenated with negative labels for training set
      
      # Third determine number of negative labels and the subset number of positive labels to be the contamination
      number.neg[i] <- floor(p[i]*length(labels.neg)) # calculates the total number of labels to be classified as negative (will contain a percentage contamination)
      number.contam[i] <- floor(contam.perc[g]*number.neg[i]) # among the number of negatives how many will be contamination based on contamination percentages
      contam.idx <- sample(1:length(labels.pos[-idx.pos]), size = number.contam[i],replace = FALSE) # among the positive labels not utilized in the training set, sample the number of contaminates determined in number.contam
      contam <- labels.pos[-idx.pos][contam.idx] # collection of row numbers for contaminates taken from poistive labels not already allocated to the training set
      labels.contam <- labels
      labels.contam[contam] <- 0
      
      # Fourth determine true negative labels and assemble full training set
      idx.neg <- sample(1:length(labels.neg), size = number.neg[i]-number.contam[i], replace = FALSE) # sample from negative labels the difference between total negatives and contamination number
      inTrain.neg.only.length[i] <- length(idx.neg) # collect number of true negatives to calculate contamination percentages
      inTrain.neg <- labels.neg[idx.neg] # collection of row numbers for negative labels only that will be included in training set
      inTrain <- c(inTrain.pos,contam,inTrain.neg) # assemble row numbers of full training set
      training <- feature.table[inTrain,] # Training set constructed using row number ids within feature table
      testing <- feature.table[-inTrain,] # Testing set constructed using row number ids within feature table
      
      repeat{
        if(product(apply(training[,c(2,4,6,7)],2,sum))!=0 & product(apply(testing[,c(2,4,6,7)],2,sum))!=0) {
          print(paste0("Seed = ",seed," doesn't contain a column of zeros."))
          break
        } else{ 
          seed <- 10000+seed+21
          set.seed(seed)
          idx.pos <- sample(1:length(labels.pos), size = number.pos[i], replace = FALSE) # sample this many (number.pos) from the positive labels
          inTrain.pos <- labels.pos[idx.pos] # row numbers of positive labels randomly selected for training set to be concatenated with negative labels for training set
          contam.idx <- sample(1:length(labels.pos[-idx.pos]), size = number.contam,replace = FALSE) # among the positive labels not utilized in the training set, sample the number of contaminates determined in number.contam
          contam <- labels.pos[-idx.pos][contam.idx] # collection of row numbers for contaminates taken from poistive labels not already allocated to the training set
          labels.contam <- labels
          labels.contam[contam] <- 0
          idx.neg <- sample(1:length(labels.neg), size = number.neg[i]-number.contam[i], replace = FALSE) # sample from negative labels the difference between total negatives and contamination number
          inTrain.neg <- labels.neg[idx.neg] # collection of row numbers for negative labels only that will be included in training set
          inTrain <- c(inTrain.pos,contam,inTrain.neg) # assemble row numbers of full training set
          training <- feature.table[inTrain,] # Training set constructed using row number ids within feature table
          testing <- feature.table[-inTrain,] # Testing set constructed using row number ids within feature table
        }
      }  
      
      # SVM modeling 
      svm.model = svm(training, y=labels.contam[inTrain], kernel="radial") # radial kernel showed the best performance from initial trials
      all.preds.svm$new = predict(svm.model,testing) # predictions of svm model on test set
      names(all.preds.svm)[length(all.preds.svm)] = paste0('SVM (train n = ', p[i]*length(labels), ')') # assign names to predictions in a list
      all.labels.svm$new = labels[-inTrain] # collects true labels in corresponding list to be used to produce ROC
      names(all.labels.svm)[length(all.labels.svm)] = paste0('SVM (train n = ', p[i]*length(labels), ')') # assign names to labels in a list
      
      # Decision Tree
      training.tree <- cbind(as.factor(labels.contam[inTrain]),training)
      colnames(training.tree)[1] <- 'labels'
      J48.model = randomForest(labels~., data = training.tree)
      predict.J48 <- predict(J48.model,testing,type = "class") # when type = "class", a factor vector is returned. When type = "prob", a matrix of confidence values is returned (one column per class).
      labels.J48 <- labels[-inTrain]
      
      all.preds.J48$new = predict.J48 # predictions of tree model on test set
      names(all.preds.J48)[length(all.preds.J48)] = paste0('Tree (train n = ', p[i]*length(labels), ')') # assign names to predictions in a list
      all.labels.J48$new = labels.J48 # collects true labels in corresponding list to be used to produce ROC
      names(all.labels.J48)[length(all.labels.J48)] = paste0('Tree (train n = ', p[i]*length(labels), ')') # assign names to labels in a list
      
      
      # LASSO modeling 
      cv.lasso <- cv.glmnet(as.matrix(training), as.double(labels.contam[inTrain]), family='binomial', alpha=1, standardize=TRUE, type.measure='auc')
      
      fit <- glmnet(as.matrix(training), as.double(labels.contam[inTrain]), family='binomial',alpha = 1,lambda=cv.lasso$lambda.min)    
      
      all.preds.sup$new = predict(fit,as.matrix(testing)) # predictions of sup model on test set
      names(all.preds.sup)[length(all.preds.sup)] = paste0('LASSO (train n = ', p[i]*length(labels), ')') # assign names to predictions in a list
      all.labels.sup$new = labels[-inTrain] # collects true labels in corresponding list to be used to produce ROC
      names(all.labels.sup)[length(all.labels.sup)] = paste0('LASSO (train n = ', p[i]*length(labels), ')') # assign names to labels in a list
      
      if(g<2){
        # Semi-sup modeling
        train.size <- floor(p[i]*constant)
        
        cat('Modeling with training size =', train.size, '\n')
        divided.data = divideYeastData(feature.table, date.table, labels,
                                       train.class = 'positive',
                                       num.train.examples = train.size, seed = seed1, constant = constant, negatives = 'all')
        
        repeat{
          if(product(apply(divided.data$train[,c(2,4,6,7)],2,sum))!=0 &
             product(apply(divided.data$test[,c(2,4,6,7)],2,sum))!=0 ) {
            print(paste0("Seed1 = ",seed1," doesn't contain a column of zeros."))
            break
          } else{ 
            seed1 <- 1000 + seed1 
            set.seed(seed1)
            divided.data = divideYeastData(feature.table, date.table, labels,
                                           train.class = 'positive',
                                           num.train.examples = train.size, seed = seed1, constant = constant, negatives = 'all')
          }
        }   
        
        is.train = divided.data$is.train
        
        lcmix.out = getLCMIXmodels.semi.only(feature.table, is.train, labels[!is.train])
        
        all.preds.semi$new = lcmix.out$preds$semisup
        names(all.preds.semi)[length(all.preds.semi)] = paste0('Semisupervised (train n = ', 
                                                               train.size, ')')
        all.labels.semi$new = labels[!divided.data$is.train]
        names(all.labels.semi)[length(all.labels.semi)] = paste0('Semisupervised (train n = ', 
                                                                 train.size, ')')
        all.perfs.semi$new = lcmix.out$performance
        names(all.perfs.semi)[length(all.perfs.semi)] = paste0('Semisupervised (train n = ', 
                                                               train.size, ')')
      }
      
    } 
    if(g<2){
      svm.pred <- prediction(all.preds.svm,all.labels.svm)
      svm.perf <- performance(svm.pred,measure="f")}
    
    'boxplot(unlist((performance(svm.pred,measure="auc"))@y.values))
    svm.perf.roc <- performance(svm.pred,measure="tpr",x.measure="fpr")
    plot(svm.perf,col="blue")
    plot(svm.perf.roc,col="blue")
    plot(performance(svm.pred,measure="prec", x.measure="rec"),col="blue")
    J48.pred <- prediction-class(all.preds.J48,all.labels.J48)
    plot(performance-class(J48.pred,measure="prec", x.measure="rec"),col="red")
    sup.pred <- prediction(all.preds.sup,all.labels.sup)
    sup.perf <- performance(sup.pred,measure="f")
    plot(performance(sup.pred,measure="f"),col="red")'
    
    ROC.svm <- multiroc(all.preds.svm,all.labels.svm) # collection of ROC calculations from prediction and label lists
    ROC.J48 <- multiroc(all.preds.J48,all.labels.J48) # collection of ROC calculations from prediction and label lists
    ROC.sup <- multiroc(all.preds.sup,all.labels.sup) # collection of ROC calculations from prediction and label lists
    if(g<2){ROC.semi <- multiroc(all.preds.semi,all.labels.semi)}
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
      AUC.J48[k,j] <- ROC.J48[[k]]$auc # extracts AUC values from ROC and places them in a matrix to be used in boxplots
      label.tree[k,j] <- percent.contam[[k]] # collects contamination percentages in each iteration of analysis
      AUC.sup[k,j] <- ROC.sup[[k]]$auc # extracts AUC values from ROC and places them in a matrix to be used in boxplots
      label.lasso[k,j] <- percent.contam[[k]] # collects contamination percentages in each iteration of analysis
      if(g<2){AUC.semi[k,j] <- ROC.semi[[k]]$auc}
    }
  }  
  
  
  svm.contam <- apply(label.svm,1,mean) # contamination percent across each row represents average contamination at each training percent (should be close to contam.perc)
  tree.contam <- apply(label.tree,1,mean) # contamination percent across each row represents average contamination at each training percent (should be close to contam.perc)
  sup.contam <- apply(label.lasso,1,mean) # contamination percent across each row represents average contamination at each training percent (should be close to contam.perc)
  
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
  
  
  AUC.list$new <- as.vector(AUC)
  names(AUC.list)[length(AUC.list)] = paste0(contam.perc[g]*100,'% contamination')
  
  length(AUC.list)
  # for(g in 1:length(contam.perc)){  # Use for rerunning plots only
  n <- rep(p,each=4*n1)
  n <- factor(n)
  method <- rep(c(rep('semi',n1),
                  rep('svm',n1),
                  rep('sup',n1),
                  rep('tree',n1)),l)
  
  tog <- data.frame(AUC=AUC.list[[g]],n=n,method=method)
  box.at <- seq(1:(l*5))
  remove <- seq(5,length(box.at),5)
  at <- box.at[!(box.at%in%remove)]
  
  png(paste0('results/plots/6C_seqFeatures_unbalanced_',contam.perc[g]*100,'.png'),
      height = 8, width = 8, res = 300, units = 'in')
  boxplot(tog$AUC~tog$method*tog$n,at=at,col=c('red','lightblue1','cadetblue3','dodgerblue2'),xlab='Training Set', ylab='AUC',xaxt="n",
          main=paste0('Supervised vs Semi-supervised (Seq Features)\nAUC for various training % with 25 iterations\nTest: unbalanced (',contam.perc[g]*100,'% Contamination)'),ylim=c(0.45,1))
  axis(side=1,at=seq(2.5,92.5,5),labels=floor(p*3500),las=2)
  legend("topleft",c("Semi-supervised","Supervised-SVM","Supervised-Lasso","Supervised-Tree"),col=c('red','lightblue1','cadetblue3','dodgerblue2'),pch=15)
  dev.off()
  
  
  # For paper without title
  png(paste0('results/plots/6C_seqFeatures_unbalanced_',contam.perc[g]*100,'_paper.png'),
      height = 8, width = 8, res = 300, units = 'in')
  boxplot(tog$AUC~tog$method*tog$n,at=at,col=c('red','lightblue1','cadetblue3','dodgerblue2'),xlab='Training Set', ylab='AUC',xaxt="n",ylim=c(0.45,1))
  axis(side=1,at=seq(2.5,92.5,5),labels=floor(p*3500),las=2)
  legend("topleft",c("Semi-supervised","Supervised-SVM","Supervised-Lasso","Supervised-Tree"),col=c('red','lightblue1','cadetblue3','dodgerblue2'),pch=15)
  dev.off()
  
}


AUC1.save <- AUC.semi

save.image("saved_images/6C_boxplots_unbalanced")
