# ---------------------------------------------------------------------------- #
# Script 2B
# Author: Rani Powers
#
# Uses all features from only positive training
# examples to classify positive and negative test examples
# ---------------------------------------------------------------------------- #

cat('SCRIPT 2B: Train + only (varying n), test +/- genes, all features\n')

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

# Train on positive genes, test on pos and neg genes

all.preds = list()
all.labels = list()
set.seed(101)
for (train.size in c(50,45,40,35,30,25,20,15,10)){
  cat('Modeling with training size =', train.size, '\n')
  divided.data = divideYeastData(feature.table, date.table, labels,
                                 train.class = 'positive', 
                                 num.train.examples = train.size)
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


# Plot all lines together
mroc = multiroc(all.preds, all.labels, quasi=F)
cols = brewer.pal(11, 'Spectral')
png('results/plots/2B-ROC_allFeatures_pos_only.png',
    height = 8, width = 8, res = 300, units = 'in')
plot(mroc, grid=TRUE, col=rep(c('red', 'black'), 4),
     cut=c(0.5, 0.8, 0.99), 
     main = 'Train: essential genes (varying n)\nTest: essential and non-essential genes\nData: 22 Seringhaus features')
dev.off()

# BoxPlots various values of n 
AUC <- c(all.preds$`Semisupervised (train n = 50)`,
                   all.preds$`Semisupervised (train n = 45)`,
                   all.preds$`Semisupervised (train n = 40)`,
                   all.preds$`Semisupervised (train n = 35)`,
                   all.preds$`Semisupervised (train n = 30)`,
                   all.preds$`Semisupervised (train n = 25)`,
                   all.preds$`Semisupervised (train n = 20)`,
                   all.preds$`Semisupervised (train n = 15)`,
                   all.preds$`Semisupervised (train n = 10)`,
                   all.preds$`Unsupervised (train n = 50)`,
                   all.preds$`Unsupervised (train n = 45)`,
                   all.preds$`Unsupervised (train n = 40)`,
                   all.preds$`Unsupervised (train n = 35)`,
                   all.preds$`Unsupervised (train n = 30)`,
                   all.preds$`Unsupervised (train n = 25)`,
                   all.preds$`Unsupervised (train n = 20)`,
                   all.preds$`Unsupervised (train n = 15)`,
                   all.preds$`Unsupervised (train n = 10)`)
                 
n <- c(rep(50,length(all.preds$`Semisupervised (train n = 50)`)),
                   rep(45,length(all.preds$`Semisupervised (train n = 45)`)),
                   rep(40,length(all.preds$`Semisupervised (train n = 40)`)),
                   rep(35,length(all.preds$`Semisupervised (train n = 35)`)),
                   rep(30,length(all.preds$`Semisupervised (train n = 30)`)),
                   rep(25,length(all.preds$`Semisupervised (train n = 25)`)),
                   rep(20,length(all.preds$`Semisupervised (train n = 20)`)),
                   rep(15,length(all.preds$`Semisupervised (train n = 15)`)),
                   rep(10,length(all.preds$`Semisupervised (train n = 10)`)),
                   rep(50,length(all.preds$`Unsupervised (train n = 50)`)),
                   rep(45,length(all.preds$`Unsupervised (train n = 45)`)),
                   rep(40,length(all.preds$`Unsupervised (train n = 40)`)),
                   rep(35,length(all.preds$`Unsupervised (train n = 35)`)),
                   rep(30,length(all.preds$`Unsupervised (train n = 30)`)),
                   rep(25,length(all.preds$`Unsupervised (train n = 25)`)),
                   rep(20,length(all.preds$`Unsupervised (train n = 20)`)),
                   rep(15,length(all.preds$`Unsupervised (train n = 15)`)),
                   rep(10,length(all.preds$`Unsupervised (train n = 10)`)))
n <- factor(n)

method <- c(rep("semi",length(c(all.preds$`Semisupervised (train n = 50)`,
                                     all.preds$`Semisupervised (train n = 45)`,
                                     all.preds$`Semisupervised (train n = 40)`,
                                     all.preds$`Semisupervised (train n = 35)`,
                                     all.preds$`Semisupervised (train n = 30)`,
                                     all.preds$`Semisupervised (train n = 25)`,
                                     all.preds$`Semisupervised (train n = 20)`,
                                     all.preds$`Semisupervised (train n = 15)`,
                                     all.preds$`Semisupervised (train n = 10)`))),
            rep("unsup",length(c(all.preds$`Unsupervised (train n = 50)`,
                                   all.preds$`Unsupervised (train n = 45)`,
                                   all.preds$`Unsupervised (train n = 40)`,
                                   all.preds$`Unsupervised (train n = 35)`,
                                   all.preds$`Unsupervised (train n = 30)`,
                                   all.preds$`Unsupervised (train n = 25)`,
                                   all.preds$`Unsupervised (train n = 20)`,
                                   all.preds$`Unsupervised (train n = 15)`,
                                   all.preds$`Unsupervised (train n = 10)`))))

tog <- data.frame(AUC,n,method)

library("ggplot2")
png('results/plots/2B-ROC_allFeatures_pos_only_ggplot.png',
    height = 8, width = 8, res = 300, units = 'in')
p <- ggplot(tog, aes(n, AUC)) 
p.box <- p + geom_boxplot(aes(colour=method),fill = "white") 
p.box
dev.off()

png('results/plots/2B-ROC_allFeatures_pos_only_boxplots.png',
    height = 8, width = 8, res = 300, units = 'in')
boxplot(AUC~method*n,at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26),col=c('aquamarine','pink'),xlab='Training Set = n', ylab='AUC',xaxt="n",
        main='Semi-supervised vs Unsupervised\nAUC for various n')
axis(side=1,at=seq(1.5,25.5,3),labels=c(10,15,20,25,30,35,40,45,50))
legend("topleft",c("semi-supervised","unsupervised"),col=c('aquamarine','pink'),pch=15)
abline(0.672,0)
text(2.5,0.7,'AUC=0.672 (n=64)')
dev.off()
