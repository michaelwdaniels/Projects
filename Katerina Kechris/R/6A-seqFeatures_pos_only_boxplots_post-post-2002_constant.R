load("saved_images/6A_boxplots_constant")
# ---------------------------------------------------------------------------- #
# Script 6A
# Author: Mike Daniels
#
# Uses sequence derived features from only positive training
# examples to classify positive and negative test examples
# ---------------------------------------------------------------------------- #

cat('SCRIPT 6A: Train + only (varying n), test +/- genes, sequence derived features\n')

# Get constants and helper functions
source('R/includes.R')
source('R/analysis_helpers.R')

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
sum(labels) # 769
# ------------------------ DIVIDE DATA INTO TRAIN/TEST ----------------------- #
# --------------------------- GET LCMIX PREDICTIONS -------------------------- #

# Train on positive genes, test on all 769 positives and a randomly chosen 769 neg genes

n1 <- 25  # Number of iterations

AUC <- matrix(0,ncol=n1,nrow=20)
rownames(AUC) <- c("S700","S600","S500","S400","S300","S200","S100","S50","S25","S10",
                   "U700","U600","U500","U400","U300","U200","U100","U50","U25","U10")
constant=769
for(j in 1:n1){
  seed1 <- 1000 + j
  cat(paste0('Iteration ',j,' of ',n1,'\n'))
  all.preds = list()
  all.labels = list()
  
  for (train.size in c(700,600,500,400,300,200,100,50,25,10)){
    cat('Modeling with training size =', train.size, '\n')
    divided.data = divideYeastData(feature.table, date.table, labels,
                                    train.class = 'positive',
                                    num.train.examples = train.size, seed = seed1, constant = constant, negatives = 'constant')
      
    is.train = divided.data$is.train
    lcmix.preds = getLCMIXmodels(feature.table, is.train, labels[!is.train])$preds
    all.preds$new = lcmix.preds$semisup
    names(all.preds)[length(all.preds)] = paste0('Semisupervised (train n = ', 
                                                 train.size, ')')
    all.preds$new = lcmix.preds$unsup
    names(all.preds)[length(all.preds)] = paste0('Unsupervised (train n = ', 
                                                 train.size, ')')
    all.labels$new = labels[!divided.data$is.train]
    names(all.labels)[length(all.labels)] = paste0('Semisupervised (train n = ', 
                                                   train.size, ')')
    all.labels$new = labels[!divided.data$is.train]
    names(all.labels)[length(all.labels)] = paste0('Unsupervised (train n = ', 
                                                   train.size, ')')
  }
  mroc = multiroc(all.preds, all.labels)
  n <- c(seq(1,19,2),seq(2,20,2))
  for(i in 1:length(n)){
    AUC[i,j] <- mroc[[n[i]]]$auc
  }
}


AUC1 <- as.vector(t(AUC))
n <- rep(rep(c(700,600,500,400,300,200,100,50,25,10),each=n1),2)
n <- factor(n)
method <- c(rep('semi',(length(AUC[,1])/2)*n1),
            rep('unsup',(length(AUC[,1])/2)*n1))

tog <- data.frame(AUC=AUC1,n=n,method=method)

" Jackknife is very time consuming - will hold off until we are sure we will use it
AUC.jack <- matrix(0,ncol=769,nrow=2)

train.size=769
for(k in 1:769){
  seed1 <- 100 + k
  cat(paste0('Iteration ',k,' of 769\n'))
  all.preds.jack = list()
  all.labels.jack = list()
  
  divided.data = divideYeastData.jack(feature.table, date.table, labels, train.class = 'positive',
                                      num.train.examples = train.size-1, seed = seed1,k=k)
  
  is.train = divided.data$is.train
  lcmix.preds = getLCMIXmodels(feature.table, is.train, labels[!is.train])$preds
  all.preds.jack$new = lcmix.preds$semisup
  names(all.preds.jack)[length(all.preds.jack)] = paste0('Semisupervised (train n = ', 
                                                         train.size-1, ')')
  all.preds.jack$new = lcmix.preds$unsup
  names(all.preds.jack)[length(all.preds.jack)] = paste0('Unsupervised (train n = ', 
                                                         train.size-1, ')')
  all.labels.jack$new = labels[!divided.data$is.train]
  names(all.labels.jack)[length(all.labels.jack)] = paste0('Semisupervised (train n = ', 
                                                           train.size-1, ')')
  all.labels.jack$new = labels[!divided.data$is.train]
  names(all.labels.jack)[length(all.labels.jack)] = paste0('Unsupervised (train n = ', 
                                                           train.size-1, ')')
  
  mroc.jack = multiroc(all.preds.jack, all.labels.jack)
  
  n <- c(1,2)
  for(i in 1:length(n)){
    AUC.jack[i,k] <- mroc.jack[[n[i]]]$auc
  }
}

AUC1.jack <- as.vector(t(AUC.jack[,1:64]))
AUC2 <- c(AUC1.jack[1:769],AUC1[1:275],AUC1.jack[770:1538],AUC1[276:550])

n <- c(rep(63,64),
       rep(60,n1),
       rep(55,n1),
       rep(50,n1),
       rep(45,n1),
       rep(40,n1),
       rep(35,n1),
       rep(30,n1),
       rep(25,n1),
       rep(20,n1),
       rep(15,n1),
       rep(10,n1),
       rep(63,64),
       rep(60,n1),
       rep(55,n1),
       rep(50,n1),
       rep(45,n1),
       rep(40,n1),
       rep(35,n1),
       rep(30,n1),
       rep(25,n1),
       rep(20,n1),
       rep(15,n1),
       rep(10,n1))


n <- factor(n)
method <- c(rep('semi',((length(AUC[,1])/2)*n1)+769),
            rep('unsup',((length(AUC[,1])/2)*n1)+769))

tog1 <- data.frame(AUC2,n,method)"

png('results/plots/6A-ROC_seqFeatures_pos_only_boxplots_post_post_2002_constant.png',
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC~tog$method*tog$n,at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29),col=c('aquamarine','pink'),xlab='Training Set = n', ylab='AUC',xaxt="n",
        main='Semi-supervised vs Unsupervised (Sequence Derived only)\nAUC for various n with 25 iterations\nTrain: pre & post/positive only - Test: pre & post/# negatives = 769',ylim=c(0.5,0.7))
axis(side=1,at=seq(1.5,29.5,3),labels=c(10,25,50,100,200,300,400,500,600,700))
legend("topleft",c("semi-supervised","unsupervised"),col=c('aquamarine','pink'),pch=15)
dev.off()

save.image("saved_images/6A_boxplots_constant")
