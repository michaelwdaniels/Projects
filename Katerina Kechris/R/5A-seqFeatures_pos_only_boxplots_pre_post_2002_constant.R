load("saved_images/5A_boxplots_constant")
# ---------------------------------------------------------------------------- #
# Script 5A
# Author: Mike Daniels
#
# Uses sequence derived features from only positive training pre-2002
# examples to classify positive and negative test examples
# ---------------------------------------------------------------------------- #

cat('SCRIPT 5A: Train + only (varying n - pre 2002), test +/- genes post 2002, seq features\n')

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


# ------------------------ DIVIDE DATA INTO TRAIN/TEST ----------------------- #
# --------------------------- GET LCMIX PREDICTIONS -------------------------- #
n1 <- 25  # Number of iterations
AUC <- matrix(0,ncol=n1,nrow=18)
rownames(AUC) <- c("S60","S55","S50","S45","S40","S35","S30","S25","S20",
                   "U60","U55","U50","U45","U40","U35","U30","U25","U20")
for(j in 25:n1){
  seed1 <- 1000 + j
  cat(paste0('Iteration ',j,' of ',n1,'\n'))
  all.preds = list()
  all.labels = list()
  
  for (train.size in c(60,55,50,45,40,35,30,25,20)){
    cat('Modeling with training size =', train.size, ' \n')
    
    divided.data = divideYeastData(feature.table, date.table, labels,
                                   train.date = 2002, train.class = 'positive',
                                   num.train.examples = train.size, seed = seed1, negatives = 'constant')
    repeat{
      if(product(apply(divided.data$train[,1:7],2,sum))!=0&(product(apply(divided.data$train[,1:7],1,sum))!=0|train.size>30)){
        print(paste0("Seed = ",seed1," doesn't contain a column of zeros nor a row of zeros."))
        break
      } else{ 
        seed1 <- seed1+j
        divided.data = divideYeastData(feature.table, date.table, labels,
                                       train.class = 'positive',
                                       num.train.examples = train.size, seed = seed1, constant = constant, negatives = 'constant')
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
  n <- c(seq(1,17,2),seq(2,18,2))
  for(i in 1:length(n)){
    AUC[i,j] <- mroc[[n[i]]]$auc
  }
}

AUC1 <- as.vector(t(AUC))

n <- c(rep(rep(c(60,55,50,45,40,35,30,25,20),each=n1),2))
n <- factor(n)
method <- c(rep('semi',(length(AUC[,1])/2)*n1),
            rep('unsup',(length(AUC[,1])/2)*n1))

tog <- data.frame(AUC1,n,method)

png('results/plots/5A-ROC_seqFeatures_pos_only_boxplots_pre_post_2002_constant.png',
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC1~tog$method*tog$n,at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26),col=c('red','blue'),xlab='Training Set = n', ylab='AUC',xaxt="n",
        main='Semi-supervised vs Unsupervised (Sequence Derived only)\nAUC for various n with 25 iterations\nTrain: pre-2002/positive - Test: pre & post/constant negatives (769)',ylim=c(0.5,0.7))
axis(side=1,at=seq(1.5,26.5,3),labels=c(20,25,30,35,40,45,50,55,60))
legend("topleft",c("Semi-supervised","Unsupervised"),col=c('red','blue'),pch=15)
dev.off()

png('results/plots/5A-ROC_seqFeatures_pos_only_boxplots_pre_post_2002_constant_paper.png',
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC1~tog$method*tog$n,at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26),col=c('red','blue'),xlab='Training Set = n', ylab='AUC',xaxt="n",ylim=c(0.5,0.7))
axis(side=1,at=seq(1.5,26.5,3),labels=c(20,25,30,35,40,45,50,55,60))
legend("topleft",c("Semi-supervised","Unsupervised"),col=c('red','blue'),pch=15)
dev.off()

save.image("saved_images/5A_boxplots_constant")
