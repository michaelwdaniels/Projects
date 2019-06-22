setwd("~/Desktop/Projects/Katerina Kechris")
load('data/Rdata/cerevisiae_preprocessed.RData')
load("saved_images/ada")

#source('R/includes.R')
#source('R/analysis_helpers.R')
#source('R/mikes_functions.R')

library(dplyr,warn.conflicts = FALSE)
library(ggplot2,warn.conflicts = FALSE)
library(AdaSampling)
library(grid)
library(gridExtra)

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

# Analysis
p <- c(seq(0.01,0.1,0.005)) # percentage training set values
l <- length(p)
n <- floor(p*length(labels)) # number of positives (negatives) in the training set at each percent
n_un <- c(50,100,300) # number of unlabeled instances in the training set
iteration <- 100

prec_ada_50 <- prec_ada_md <- prec_ada_80 <- prec_ada_95 <- matrix(0,nrow=iteration,ncol=length(n))
rec_ada_50 <- rec_ada_md <- rec_ada_80 <- rec_ada_95 <- matrix(0,nrow=iteration,ncol=length(n))

labels.pos <- which(labels==1) # collection of row numbers from labels which are positive (769)
labels.neg <- which(labels==0) # collection of row numbers from labels which are negative (2731)
number.in.training.set <- vector() # list of number of positive labels in training set
number.un <- vector() # list of number of positive labels in training set

# We need a matrix for allotment of training and test subsets sizes
nn <- matrix(0,ncol=length(p),nrow=4)
rownames(nn) <- c('n_pos','n_neg','n_na','n_test')
colnames(nn) <- paste0(p*100,'%')

AUC <- list(matrix(0,ncol=l,nrow=iteration),matrix(0,ncol=l,nrow=iteration),matrix(0,ncol=l,nrow=iteration))

for(k in 1:length(n_un)){
  cat('Analysis using ',n_un[k],' unlabeled training instances. \n')
  set.seed(1)
nn[1,] <- n
nn[2,] <- n
nn[3,] <- n_un[k]
nn[4,] <- apply(nn,2,function(x) 3500-sum(x[1:3]))
write.csv(nn,paste0('trainTestCount',n_un[k],'.csv'))

for(j in 1:iteration){
  cat('Iteration = ',j,' of ',iteration,'\n')

for(i in 1:l){

  # The training and testing scheme is complex.
  # 1. Randomly chose p[i]*3500=n[i] from positive labels (labels.pos) and the same from negative labels (labels.neg) as row idx
  # 2. Concatenate unused rows c(idxp,idxn)
  # 3. Randomly choose 300 for the unlabeled training set
  # 4. Place remaining genes in test set
  
  idxp <- labels.pos[sample(1:length(labels.pos),n[i])] # idx for positive labels in training set
  idxn <- labels.neg[sample(1:length(labels.neg),n[i])] # idx for negative labels in training set
  idxLabeled <- c(idxp,idxn) # idx for labeled data in training set
  idxAll <- 1:3500 # set of all row numbers
  idxUnlabeled <- idxAll[!idxAll %in% idxLabeled] # removed ids from labeled data
  idxNA <- idxUnlabeled[sample(1:length(idxUnlabeled),n_un[k])] # sampled n_un=300 from remaining subset
  idxTraining <- c(idxp,idxn,idxNA) # concatenate training set ids
  trainingSet <- feature.table[idxTraining,]
  testingSet <- feature.table[-idxTraining,]
  labels.train <- labels[idxTraining] # labels for training set
  labels.train[2*n[i]+1:length(labels.train)] <- NA # remove labels for unlabeling set in training
  
# AdaBoost
# The algorithm works by selecting random subsets of the training set and designs a decision tree
# for each tree to estimate the missing labels in the training set. 
# Learning from positive and unlabeled data frequently occurs in applications
# where only a subset of positive instances is available while the rest of the data
# are unlabeled. In such scenarios, often the goal is to create a discriminant model
# that can accurately classify both positive and negative data by modelling from labeled and 
# unlabeled instances. In this study, we propose an adaptive sampling (AdaSampling) 
# approach that utilises prediction probabilities from a model to iteratively update
# the training data. Starting with equal prior probabilities for all unlabeled data,
# our method “wraps” around a predictive model to iteratively update these probabilities
# to distinguish positive and negative instances in unlabeled data. Subsequently, one or
# more robust negative set(s) can be drawn from unlabeled data, according to the likelihood 
# of each instance being negative, to train a single classification model or ensemble
# of models.

t_semi_ada_knn <- adaSample(Ps=rownames(trainingSet[which(labels.train==1),]),Ns=rownames(trainingSet[which(labels.train==0),]),
                            train.mat=trainingSet,test.mat=testingSet,classifier='knn')
#t_semi_ada_logit <- adaSample(Ps=rownames(X[which(y==1),]),Ns=rownames(X[which(y==0),]),
#                              train.mat=X_train,test.mat=X_test,classifier='logit')

pred_cut_50 <- ifelse(t_semi_ada_knn[,1]<0.5,0,1)
pred_cut_md <- ifelse(t_semi_ada_knn[,1]<median(t_semi_ada_knn[,1]),0,1) # 0.6
pred_cut_80 <- ifelse(t_semi_ada_knn[,1]<0.80,0,1)
pred_cut_95 <- ifelse(t_semi_ada_knn[,1]<0.95,0,1)

confusion_50 <- cbind(predicted=pred_cut_50,truth=labels[-idxTraining])
confusion_md <- cbind(predicted=pred_cut_md,truth=labels[-idxTraining])
confusion_80 <- cbind(predicted=pred_cut_80,truth=labels[-idxTraining])
confusion_95 <- cbind(predicted=pred_cut_95,truth=labels[-idxTraining])

# precision measures proportion of true positives among predicted positives TP/(TP+FP)
# idx predicted positives (predicted==1)
# subset confusion matrix to just predicted positives
# sum up the truth column to get the TP and divide by its length
prec_ada_50[j,i] <- ifelse(length(which(confusion_50[,1]==1))>1,sum(confusion_50[which(confusion_50[,1]==1),][,2])/length(confusion_50[which(confusion_50[,1]==1),][,2]),NA)
prec_ada_md[j,i] <- ifelse(length(which(confusion_md[,1]==1))>1,sum(confusion_md[which(confusion_md[,1]==1),][,2])/length(confusion_md[which(confusion_md[,1]==1),][,2]),NA)
prec_ada_80[j,i] <- ifelse(length(which(confusion_80[,1]==1))>1,sum(confusion_80[which(confusion_80[,1]==1),][,2])/length(confusion_80[which(confusion_80[,1]==1),][,2]),NA)
prec_ada_95[j,i] <- ifelse(length(which(confusion_95[,1]==1))>1,sum(confusion_95[which(confusion_95[,1]==1),][,2])/length(confusion_95[which(confusion_95[,1]==1),][,2]),NA)

# recall measures proportion of true positives predicted TP/P
# idx true positives (truth==1)
# subset confusion matrix to just true positives
# sum up the predicted column and divide by its length
rec_ada_50[j,i] <- ifelse(length(which(confusion_50[,1]==1))>1&length(which(confusion_50[,2]==1))>1,sum(confusion_50[which(confusion_50[,1]==1),][,2])/length(confusion_50[which(confusion_50[,2]==1),][,2]),NA)
rec_ada_md[j,i] <- ifelse(length(which(confusion_md[,1]==1))>1&length(which(confusion_md[,2]==1))>1,sum(confusion_md[which(confusion_md[,1]==1),][,2])/length(confusion_md[which(confusion_md[,2]==1),][,2]),NA)
rec_ada_80[j,i] <- ifelse(length(which(confusion_80[,1]==1))>1&length(which(confusion_80[,2]==1))>1,sum(confusion_80[which(confusion_80[,1]==1),][,2])/length(confusion_80[which(confusion_80[,2]==1),][,2]),NA)
rec_ada_95[j,i] <- ifelse(length(which(confusion_95[,1]==1))>1&length(which(confusion_95[,2]==1))>1,sum(confusion_95[which(confusion_95[,1]==1),][,2])/length(confusion_95[which(confusion_95[,2]==1),][,2]),NA)
                   
# f-measure is a geometric mean of precision and recall
# AUC

true_Y = as.vector(confusion_50[,2])
probs = as.vector(t_semi_ada_knn[,1])

getROC_AUC = function(probs, true_Y){
  probsSort = sort(probs, decreasing = TRUE, index.return = TRUE)
  val = unlist(probsSort$x)
  idx = unlist(probsSort$ix)  
  
  roc_y = true_Y[idx];
  stack_x = cumsum(roc_y == 0)/sum(roc_y == 0)
  stack_y = cumsum(roc_y == 1)/sum(roc_y == 1)    
  
  auc = sum((stack_x[2:length(roc_y)]-stack_x[1:length(roc_y)-1])*stack_y[2:length(roc_y)])
  return(list(stack_x=stack_x, stack_y=stack_y, auc=auc))
}

aList = getROC_AUC(probs, true_Y) 

stack_x = unlist(aList$stack_x)
stack_y = unlist(aList$stack_y)
auc = unlist(aList$auc)

AUC[[k]][j,i] <- auc

#plot(stack_x, stack_y, type = "l", col = "blue", xlab = "False Positive Rate", ylab = "True Positive Rate", main = "ROC")
#axis(1, seq(0.0,1.0,0.1))
#axis(2, seq(0.0,1.0,0.1))
#abline(h=seq(0.0,1.0,0.1), v=seq(0.0,1.0,0.1), col="gray", lty=3)
#legend(0.7, 0.3, sprintf("%3.3f",auc), lty=c(1,1), lwd=c(2.5,2.5), col="blue", title = "AUC")

}
}

write.csv(AUC,paste0('results/tables/AUC',n_un[k],'.csv'))

ada_prec_50 <- replace(prec_ada_50,prec_ada_50==0,NA)
ada_prec_md <- replace(prec_ada_md,prec_ada_md==0,NA)
ada_prec_80 <- replace(prec_ada_80,prec_ada_80==0,NA)
ada_prec_95 <- replace(prec_ada_95,prec_ada_95==0,NA)

ada_rec_50 <- replace(rec_ada_50,rec_ada_50==0,NA)
ada_rec_md <- replace(rec_ada_md,rec_ada_md==0,NA)
ada_rec_80 <- replace(rec_ada_80,rec_ada_80==0,NA)
ada_rec_95 <- replace(rec_ada_95,rec_ada_95==0,NA)

prec_50_avg <- apply(ada_prec_50,2,function(x) mean(x,na.rm=T))
prec_md_avg <- apply(ada_prec_md,2,function(x) mean(x,na.rm=T))
prec_80_avg <- apply(ada_prec_80,2,function(x) mean(x,na.rm=T))
prec_95_avg <- apply(ada_prec_95,2,function(x) mean(x,na.rm=T))

rec_50_avg <- apply(ada_rec_50,2,function(x) mean(x,na.rm=T))
rec_md_avg <- apply(ada_rec_md,2,function(x) mean(x,na.rm=T))
rec_80_avg <- apply(ada_rec_80,2,function(x) mean(x,na.rm=T))
rec_95_avg <- apply(ada_rec_95,2,function(x) mean(x,na.rm=T))

f_50 <- 2/((1/prec_50_avg)+(1/rec_50_avg))
f_md <- 2/((1/prec_md_avg)+(1/rec_md_avg))
f_80 <- 2/((1/prec_80_avg)+(1/rec_80_avg))
f_95 <- 2/((1/prec_95_avg)+(1/rec_95_avg))

prec_ada <- rbind(prec_50_avg,prec_md_avg,prec_80_avg,prec_95_avg)
rec_ada <- rbind(rec_50_avg,rec_md_avg,rec_80_avg,rec_95_avg)
f_ada <- rbind(f_50,f_md,f_80,f_95)
colnames(prec_ada) <- colnames(rec_ada) <- colnames(f_ada) <- p*100
prec_ada
rec_ada
f_ada

write.csv(prec_ada,paste0('results/tables/adaPrecision',n_un[k],'.csv'))
write.csv(rec_ada,paste0('results/tables/adaRecall',n_un[k],'.csv'))
write.csv(f_ada,paste0('results/tables/adaf',n_un[k],'.csv'))

performance.vector.0.ada.50 <- c(prec_ada[1,],rec_ada[1,],f_ada[1,],summary0.med[4,],summary0.med[5,],summary0.med[6,])
performance.vector.0.ada.md <- c(prec_ada[2,],rec_ada[2,],f_ada[2,],summary0.med[4,],summary0.med[5,],summary0.med[6,])
performance.vector.0.ada.80 <- c(prec_ada[3,],rec_ada[3,],f_ada[3,],summary0.med[4,],summary0.med[5,],summary0.med[6,])
performance.vector.0.ada.95 <- c(prec_ada[4,],rec_ada[4,],f_ada[4,],summary0.med[4,],summary0.med[5,],summary0.med[6,])

perc.2 <- rep(as.numeric(rownames(method.means)),6)*100
measure.2 <- rep(c(rep('Precision',l),rep('Recall',l),rep('f-measure',l)),2)
method.2 <- c(rep('Semi-ADA',l*3),rep('Semi-HM',l*3))
group.2 <- interaction(measure.2,method.2)

ada.0.50 <- data.frame(percent=perc.2,performance=performance.vector.0.ada.50,measure=measure.2,method=method.2,group=group.2)
ada.0.md <- data.frame(percent=perc.2,performance=performance.vector.0.ada.md,measure=measure.2,method=method.2,group=group.2)
ada.0.80 <- data.frame(percent=perc.2,performance=performance.vector.0.ada.80,measure=measure.2,method=method.2,group=group.2)
ada.0.95 <- data.frame(percent=perc.2,performance=performance.vector.0.ada.95,measure=measure.2,method=method.2,group=group.2)

levels(ada.0.50$method) <- c("Semi-ADA","Semi-HM")
p7.ada.50 <- ggplot(ada.0.50, aes(x=percent,y=performance,group=group,colour=method.2)) + geom_line(size=1,aes(linetype=measure)) +
  ylim(0,1) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed')) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'))
p7.ada.50 <- p7.ada.50 + xlab('Training Set (%)') + ylab('Performance') + labs(caption = 'cutoff @ 0.50')

levels(ada.0.md$method) <- c("Semi-ADA","Semi-HM")
p7.ada.md <- ggplot(ada.0.md, aes(x=percent,y=performance,group=group,colour=method.2)) + geom_line(size=1,aes(linetype=measure)) +
  ylim(0,1) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed')) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'))
p7.ada.md <- p7.ada.md + xlab('Training Set (%)') + ylab('Performance') + labs(caption = 'cutoff @ median')

levels(ada.0.80$method) <- c("Semi-ADA","Semi-HM")
p7.ada.80 <- ggplot(ada.0.80, aes(x=percent,y=performance,group=group,colour=method.2)) + geom_line(size=1,aes(linetype=measure)) +
  ylim(0,1) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed')) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'))
p7.ada.80 <- p7.ada.80 + xlab('Training Set (%)') + ylab('Performance') + labs(caption = 'cutoff @ 0.80')

levels(ada.0.95$method) <- c("Semi-ADA","Semi-HM")
p7.ada.95 <- ggplot(ada.0.95, aes(x=percent,y=performance,group=group,colour=method.2)) + geom_line(size=1,aes(linetype=measure)) +
  ylim(0,1) + xlim(1,5) + theme(legend.position = 'none') + scale_colour_manual(values=c('dodgerblue2','red')) + scale_linetype_manual(values=c("solid", "dotted",'dashed')) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16,face = 'bold'),
        axis.title=element_text(size=16,face = 'bold'))
p7.ada.95 <- p7.ada.95 + xlab('Training Set (%)') + ylab('Performance') + labs(caption = 'cutoff @ 0.95')

png(paste0('results/plots/measures_ada_',n_un[k],'.png'),height = 8, width = 8, res = 300, units = 'in')
grid.arrange(p7.ada.50,p7.ada.80,p7.ada.95)
dev.off()
}

AUC.list1 <- list()
  
for(k in 1:length(n_un)){
#AUC1 <- t(AUC.semi)
AUC2 <- AUC[[k]]
#AUC3 <- t(AUC.sup)
#AUC4 <- t(AUC.J48)
AUC.ada <- matrix(0,ncol=2*l,nrow=iteration)

for(f in 1:ncol(AUC2)){
  #AUC[,2*f-3] <- AUC1[,f]
  #AUC[,4*f-2] <- AUC2[,f]
  AUC.ada[,2*f-1] <- AUC1[,f]
  AUC.ada[,2*f]   <- AUC2[,f]
}

AUC.list1$new <- as.vector(AUC.ada)
names(AUC.list1)[length(AUC.list1)] = n_un[k]

#for(g in 1:length(contam.perc)){  # Use for rerunning plots only
n <- rep(p,each=2*iteration)
n <- factor(n)
method <- rep(c(rep('semi',iteration),rep('ada',iteration)),l)

tog <- data.frame(AUC=AUC.list1[[k]],n=n,method=method)
levels(tog$method) <- c('semi','ada')
box.at <- seq(1:(l*3))
remove <- seq(3,length(box.at),3)
at <- box.at[!(box.at%in%remove)]

png(paste0('results/plots/6D_allFeatures_unbalanced_ada_',n_un[k],'.png'),
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC~tog$method*tog$n,at=at,col=c('dodgerblue2','red'),xlab='Training Set (%)', ylab='AUC',xaxt="n",
        main=paste0('AdaBoost vs Semi-HM (19 Features)\nAUC for various training % with ',iteration,' iterations'),ylim=c(0.38,0.8),cex.axis=1.5,cex.lab=1.4)
axis(side=1,at=seq(2.5,92.5,5),labels=seq(1,10,0.5),las=2,cex.axis=1.5)
legend("topleft",c("Semi-HM","AdaBoost"),col=c('red','dodgerblue2'),pch=15)
dev.off()


# For paper without title
png(paste0('results/plots/6D_allFeatures_unbalanced_ada_',n_un[k],'_paper.png'),
    height = 8, width = 8, res = 300, units = 'in')
boxplot(tog$AUC~tog$method*tog$n,at=at,col=c('dodgerblue2','red'),xlab='Training Set (%)', ylab='AUC',xaxt="n",ylim=c(0.38,0.8),cex.axis=1.5,cex.lab=1.4)
axis(side=1,at=seq(2.5,92.5,5),labels=seq(1,10,0.5),las=2,cex.axis=1.5)
legend("topleft",c("Semi-HM","AdBoost"),col=c('red','dodgerblue2'),pch=15)
dev.off()
}

# save.image("saved_images/ada")
