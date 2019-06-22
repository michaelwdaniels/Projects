# A hub to rerun all plots for contamination



setwd("~/Desktop/semisupervised-mixture-models-master")

###################################################################################################
# All Features Unbalanced
###################################################################################################
source("R/6D_allFeatures_sup_unbalanced.R")

for(g in 1:length(g)){
# For paper without title
png(paste0('results/plots/All Features Unbalanced/6D_allFeatures_unbalanced_',contam.perc[g]*100,'_paper.png'),
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC~tog$method*tog$n,at=at,col=c('red','lightblue1','cadetblue3','dodgerblue2'),xlab='Training Set', ylab='AUC',xaxt="n",ylim=c(0.45,1))
axis(side=1,at=seq(2.5,92.5,5),labels=floor(p*3500),las=2)
legend("topleft",c("Semi-supervised","Supervised-SVM","Supervised-Lasso","Supervised-Tree"),col=c('red','lightblue1','cadetblue3','dodgerblue2'),pch=15)
dev.off()
}

###################################################################################################
# All Features Balanced
###################################################################################################
source("R/6D_allFeatures_sup_balanced.R")

for(g in 1:length(g)){
# For paper without title
png(paste0('results/plots/All Features Balanced/6D_allFeatures_balanced_',contam.perc[g]*100,'_paper.png'),
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC~tog$method*tog$n,at=at,col=c('red','lightblue1','cadetblue3','dodgerblue2'),xlab='Training Set', ylab='AUC',xaxt="n",ylim=c(0.45,1))
axis(side=1,at=seq(2.5,92.5,5),labels=floor(p*3500),las=2)
legend("topleft",c("Semi-supervised","Supervised-SVM","Supervised-Lasso","Supervised-Tree"),col=c('red','lightblue1','cadetblue3','dodgerblue2'),pch=15)
dev.off()
}

###################################################################################################
# Sequenced Derived Features Unbalanced
###################################################################################################
source("R/6C_seqFeatures_sup_unbalanced.R")

for(g in 1:length(g)){
# For paper without title
png(paste0('results/plots/Seq Features Unbalanced/6D_allFeatures_unbalanced_',contam.perc[g]*100,'_paper.png'),
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC~tog$method*tog$n,at=at,col=c('red','lightblue1','cadetblue3','dodgerblue2'),xlab='Training Set', ylab='AUC',xaxt="n",ylim=c(0.45,1))
axis(side=1,at=seq(2.5,92.5,5),labels=floor(p*3500),las=2)
legend("topleft",c("Semi-supervised","Supervised-SVM","Supervised-Lasso","Supervised-Tree"),col=c('red','lightblue1','cadetblue3','dodgerblue2'),pch=15)
dev.off()
}

###################################################################################################
# Sequence Derived Features Balanced
###################################################################################################
source("R/6C_seqFeatures_sup_balanced.R")

for(g in 1:length(g)){
# For paper without title
png(paste0('results/plots/Seq Features Balanced/6D_allFeatures_balanced_',contam.perc[g]*100,'_paper.png'),
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC~tog$method*tog$n,at=at,col=c('red','lightblue1','cadetblue3','dodgerblue2'),xlab='Training Set', ylab='AUC',xaxt="n",ylim=c(0.45,1))
axis(side=1,at=seq(2.5,92.5,5),labels=floor(p*3500),las=2)
legend("topleft",c("Semi-supervised","Supervised-SVM","Supervised-Lasso","Supervised-Tree"),col=c('red','lightblue1','cadetblue3','dodgerblue2'),pch=15)
dev.off()
}
