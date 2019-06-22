load("saved_images/6A_boxplots_match_test")
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
constant = sum(labels) # 769
# ------------------------ DIVIDE DATA INTO TRAIN/TEST ----------------------- #
# --------------------------- GET LCMIX PREDICTIONS -------------------------- #

# Train on positive genes, test on all 769 positives and a randomly chosen 769 neg genes

n1 <- 25  # Number of iterations

AUC <- matrix(0,ncol=n1,nrow=18)
rownames(AUC) <- c("S700","S600","S500","S400","S300","S200","S100","S50","S25",
                   "U700","U600","U500","U400","U300","U200","U100","U50","U25")

for(j in 1:n1){
  seed1 <- 100 + j
  cat(paste0('Iteration ',j,' of ',n1,'\n'))
  all.preds = list()
  all.labels = list()
  
  for (train.size in c(700,600,500,400,300,200,100,50,25)){
    cat('Modeling with training size =', train.size, '\n')
    divided.data = divideYeastData(feature.table, date.table, labels,
                                    train.class = 'positive',
                                    num.train.examples = train.size, seed = seed1, constant = constant, negatives = 'match_test')
    repeat{
      if(product(apply(divided.data$train[,1:7],2,sum))!=0){
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
  n <- c(seq(1,17,2),seq(2,18,2))
  for(i in 1:length(n)){
    AUC[i,j] <- mroc[[n[i]]]$auc
  }
}


AUC1 <- as.vector(t(AUC))
n <- rep(rep(c(700,600,500,400,300,200,100,50,25),each=n1),2)
n <- factor(n)
method <- c(rep('semi',(length(AUC[,1])/2)*n1),
            rep('unsup',(length(AUC[,1])/2)*n1))

tog <- data.frame(AUC=AUC1,n=n,method=method)

png('results/plots/6A-ROC_seqFeatures_pos_only_boxplots_post_post_2002_match_test.png',
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC~tog$method*tog$n,at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26),col=c('aquamarine','pink'),xlab='Training Set = n', ylab='AUC',xaxt="n",
        main='Semi-supervised vs Unsupervised (Sequence Derived only)\nAUC for various n with 25 iterations\nTrain: pre & post/positive only - Test: pre & post/equal sizes of positives & negatives',ylim=c(0.5,0.7))
axis(side=1,at=seq(1.5,26.5,3),labels=c(25,50,100,200,300,400,500,600,700))
legend("topleft",c("Semi-supervised","Unsupervised"),col=c('aquamarine','pink'),pch=15)
dev.off()

# For paper without title
png('results/plots/6A-ROC_allFeatures_pos_only_boxplots_post_post_2002_match_test_paper.png',
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC~tog$method*tog$n,at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26),col=c('red','blue'),xlab='% Training Set', ylab='AUC',xaxt="n",ylim=c(0.50,0.75))
axis(side=1,at=seq(1.5,26.5,3),labels=c(25,50,100,200,300,400,500,600,700))
legend("topleft",c("Semi-supervised","Unsupervised"),col=c('red','blue'),pch=15)
dev.off()

save.image("saved_images/6A_boxplots_match_test")
