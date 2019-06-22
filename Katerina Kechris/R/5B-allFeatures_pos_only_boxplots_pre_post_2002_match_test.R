load("saved_images/5B_boxplots_match_test")
# ---------------------------------------------------------------------------- #
# Script 5B
# Author: Mike Daniels
#
# Uses all features from only positive training
# examples to classify positive and negative test examples
# ---------------------------------------------------------------------------- #

cat('SCRIPT 5B: Train + only (varying n - pre 2002), test +/- genes post 2002, all features\n')

# Get constants and helper functions
source('R/includes.R')
source('R/analysis_helpers.R')

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

# ------------------------ DIVIDE DATA INTO TRAIN/TEST ----------------------- #
# --------------------------- GET LCMIX PREDICTIONS -------------------------- #
n1 <- 25  # Number of iterations
AUC <- matrix(0,ncol=n1,nrow=16)
rownames(AUC) <- c("S60","S55","S50","S45","S40","S35","S30","S25",
                   "U60","U55","U50","U45","U40","U35","U30","U25")
for(j in 1:n1){
  seed1 <- 100 + j
  cat(paste0('Iteration ',j,' of ',n1,'\n'))
  all.preds = list()
  all.labels = list()
  
  for (train.size in c(60,55,50,45,40,35,30,25)){
    cat('Modeling with training size =', train.size, ' \n')
    
    divided.data = divideYeastData(feature.table, date.table, labels,
                                    train.date = 2002, train.class = 'positive',
                                    num.train.examples = train.size, seed = seed1, negatives = 'match_test')
    repeat{
      if(product(apply(divided.data$train,2,sum))!=0&rankMatrix(divided.data$train)>=ncol(divided.data$train)){
        print(paste0("Seed = ",seed1," doesn't contain a column of zeros."))
        break
      } else{ 
        seed1 <- seed1+j
        divided.data = divideYeastData(feature.table, date.table, labels,
                                       train.class = 'positive',
                                       num.train.examples = train.size, seed = seed1, constant = constant, negatives = 'match_test')
      }
    }   
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
  n <- c(seq(1,15,2),seq(2,16,2))
  for(i in 1:length(n)){
    AUC[i,j] <- mroc[[n[i]]]$auc
  }
}


AUC1 <- as.vector(t(AUC[,1:n1]))

n <- c(rep(rep(c(60,55,50,45,40,35,30,25),each=n1),2))
n <- factor(n)
method <- c(rep('semi',(length(AUC[,1])/2)*n1),
            rep('unsup',(length(AUC[,1])/2)*n1))

tog <- data.frame(AUC1,n,method)

AUC.jack <- matrix(0,ncol=64,nrow=2)

train.size=64
for(k in 1:64){
  seed1 <- 100 + k
  cat(paste0('Iteration ',k,' of 64\n'))
  all.preds.jack = list()
  all.labels.jack = list()
  
  divided.data = divideYeastData.jack(feature.table, date.table, labels,
                                      train.date = 2002, train.class = 'positive',
                                      num.train.examples = train.size-1, seed = seed1,k=k, negatives = 'match_test')
  repeat{
    if(product(apply(divided.data$train,2,sum))!=0&rankMatrix(divided.data$train)>=ncol(divided.data$train)){
      print(paste0("Seed = ",seed1," doesn't contain a column of zeros."))
      break
    } else{ 
      seed1 <- seed1+j
      divided.data = divideYeastData(feature.table, date.table, labels,
                                     train.class = 'positive',
                                     num.train.examples = train.size, seed = seed1, constant = constant, negatives = 'match_test')
    }
  }   
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
AUC2 <- c(AUC1.jack[1:64],AUC1[1:200],AUC1.jack[65:128],AUC1[201:400])

n <- c(rep(c(rep(63,64),rep(c(60,55,50,45,40,35,30,25),each=n1)),2))

n <- factor(n)
method <- c(rep('semi',((length(AUC[,1])/2)*n1)+64),
            rep('unsup',((length(AUC[,1])/2)*n1)+64))

tog1 <- data.frame(AUC2,n,method)

png('results/plots/5B-ROC_allFeatures_pos_only_boxplots_pre_post_2002_match_test.png',
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog1$AUC2~tog1$method*tog1$n,at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26),col=c('aquamarine','pink'),xlab='Training Set = n', ylab='AUC',xaxt="n",
        main='Semi-supervised vs Unsupervised (All Features)\nAUC for various n with 25 iterations\nTrain: pre-2002/positive - Test: pre & post/equal sizes of positives & negatives',ylim=c(0.5,0.7))
axis(side=1,at=seq(1.5,26.5,3),labels=c(25,30,35,40,45,50,55,60,63))
legend("topleft",c("semi-supervised","unsupervised"),col=c('aquamarine','pink'),pch=15)
dev.off()

save.image("saved_images/5B_boxplots_match_test")
