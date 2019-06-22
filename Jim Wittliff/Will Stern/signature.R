setwd('~/Desktop/Projects/Jim Wittliff/Will Stern')
load('code/signatureWill')
qPCR_ER_PR <- read.csv("~/Desktop/Projects/Jim Wittliff/Will Stern/data/qPCR_ER_PR.csv")

peptide <- qPCR_ER_PR
head(peptide)
dim(peptide) # 579 88
names <- names(peptide)[-c(1:3)]

# Survival Calculation
library(survival)
erpos <- which(peptide$ER.Status=='POS')
erneg <- which(peptide$ER.Status=='NEG')
prpos <- which(peptide$PR.Status=='POS')
prneg <- which(peptide$PR.Status=='NEG')
length(erpos) # 427
length(erneg) # 152
length(prpos) # 420
length(prneg) # 159

erposPrpos <- which(peptide$ER.Status=='POS'&peptide$PR.Status=='POS')
ernegPrpos <- which(peptide$ER.Status=='NEG'&peptide$PR.Status=='POS')
erposPrneg <- which(peptide$ER.Status=='POS'&peptide$PR.Status=='NEG')
ernegPrneg <- which(peptide$ER.Status=='NEG'&peptide$PR.Status=='NEG')
length(erposPrpos) # 369
length(ernegPrpos) # 51
length(erposPrneg) # 58
length(ernegPrneg) # 101

idxERPR <- list(erpos,erneg,prpos,prneg,erposPrpos,ernegPrpos,erposPrneg,ernegPrneg)
namesERPR <- c('erpos','erneg','prpos','prneg','erposPrpos','ernegPrpos','erposPrneg','ernegPrneg')

library(SGL)

# Retrieve survival data from EIA LBA csv
data.surv <- read.csv("~/Desktop/Projects/Jim Wittliff/Will Stern/data/Total EIA LBA Extraction.csv")
head(data.surv)

names(peptide)[1] <- 'ID'
names(data.surv)[1] <- 'ID'
library(plyr)
xxx <- join(peptide, data.surv, by = c("ID"))
#head(xxx)
#dim(xxx) # 579 98

# First glance of the data reveals a missing data problem
# I can use the mice package to impute the values
# No time for diagnostics - can do before paper but we need this abstract
# we can use a L2 regularization to get a signature

# Multiple imputations by chained equations
library(mice)
set.seed(1)
peptide.mice <- mice(xxx[,-c(1:3)],m=1)
peptide.comp <- complete(peptide.mice)
#str(peptide.comp)
#peptide.mira <- with(peptide.comp,lassoFxPen())
#peptide.pool <- pool(peptide.mira)

#dim(peptide.comp) # 579 95
for(i in 1:length(idxERPR)){
  # DFS
  surv1 <- cbind(xxx$DFS,xxx$DFS.Event,peptide.comp)
  ttt <- complete.cases(surv1[idxERPR[[i]],])
  #optL1 <- optL1(Surv(surv1$DFS[ttt],surv1$DFS.Event[ttt]), penalized = surv1[ttt,-c(1:2,88:97)])
  #optL2 <- optL2(Surv(surv1$DFS[ttt],surv1$DFS.Event[ttt]), penalized = surv1[ttt,-c(1:2,88:97)])
  #fit.pen <- penalized(Surv(surv1$DFS[ttt],surv1$DFS.Event[ttt]), penalized = surv1[ttt,-c(1:2,88:97)],lambda2 = optL2$lambda)

  data.sgl <- list(x=surv1[ttt,-c(1:2,88:97)], time = surv1$DFS[ttt], status = surv1$DFS.Event[ttt])
  sgl.index <- rep(1,ncol(surv1[ttt,-c(1:2,88:97)]))
  #fit.sgl <- SGL(data.sgl, sgl.index, type = "cox")
  cv.sgl <- cvSGL(data.sgl, sgl.index, type = "cox")
  idxx <- which(!cv.sgl$fit$beta[,which.min(cv.sgl$lldiff)]==0)
  dfs.signature <- as.data.frame(cbind(colnames(cv.sgl$fit$X)[idxx],round(as.numeric(cv.sgl$fit$beta[idxx,4]),3)))
  colnames(dfs.signature) <- c('Genes','Beta')
  write.csv(dfs.signature,paste0('code/',namesERPR[i],'.dfs.signature.csv'))

  # OS
  surv2 <- cbind(xxx$OS,xxx$OS.Event,peptide.comp)
  ttt <- complete.cases(surv2)
  #optL1 <- optL1(Surv(surv2$OS[ttt],surv2$OS.Event[ttt]), penalized = surv2[ttt,-c(1:2,88:97)])
  #optL2 <- optL2(Surv(surv2$OS[ttt],surv2$OS.Event[ttt]), penalized = surv2[ttt,-c(1:2,88:97)])
  #fit.pen <- penalized(Surv(surv2$OS[ttt],surv2$OS.Event[ttt]), penalized = surv2[ttt,-c(1:2,88:97)],lambda2 = optL2$lambda)

  data.sgl <- list(x=surv2[ttt,-c(1:2,88:97)], time = surv2$OS[ttt], status = surv2$OS.Event[ttt])
  sgl.index <- rep(1,ncol(surv2[ttt,-c(1:2,88:97)]))
  #fit.sgl <- SGL(data.sgl, sgl.index, type = "cox")
  cv.sgl <- cvSGL(data.sgl, sgl.index, type = "cox")
  idxx <- which(!cv.sgl$fit$beta[,which.min(cv.sgl$lldiff)]==0)
  os.signature <- as.data.frame(cbind(colnames(cv.sgl$fit$X)[idxx],round(as.numeric(cv.sgl$fit$beta[idxx,which.min(cv.sgl$lldiff)]),3)))
  colnames(os.signature) <- c('Genes','Beta')
  write.csv(os.signature,paste0('code/',namesERPR[i],'.os.signature.csv'))
}

save.image('code/signatureWill')

# Survival Calculation
library(survival)
erpos <- which(peptide$ER.Status==1)
erneg <- which(peptide$ER.Status==0)
prpos <- which(peptide$PR.Status==1)
prneg <- which(peptide$PR.Status==0)
length(erpos) #147
length(erneg) #100
length(prpos) #145
length(prneg) #102

names <- colnames(peptide[2:143])

models.DFS <- vector("list", length(names))
models.DFS.pr <- vector("list", length(names))
models.DFS.erpr <- vector("list", length(names))
models.DFS.erpos <- vector("list", length(names))
models.DFS.erneg <- vector("list", length(names))
models.DFS.prpos <- vector("list", length(names))
models.DFS.prneg <- vector("list", length(names))
models.OS <- vector("list", length(names))
models.OS.pr <- vector("list", length(names))
models.OS.erpr <- vector("list", length(names))
models.OS.erpos <- vector("list", length(names))
models.OS.erneg <- vector("list", length(names))
models.OS.prpos <- vector("list", length(names))
models.OS.prneg <- vector("list", length(names))

names(models.OS) <- names(models.OS.pr) <- names(models.OS.erpr) <- names(models.OS.erpos) <- names(models.OS.erneg) <- names(models.OS.prpos) <- names(models.OS.prneg) <- names
names(models.DFS) <- names(models.DFS.pr) <- names(models.DFS.erpr) <- names(models.DFS.erpos) <- names(models.DFS.erneg) <- names(models.DFS.prpos) <- names(models.DFS.prneg) <- names
peptide <- as.data.frame(peptide)

## loops through each gene
for (i in 1:length(names)) {
  models.DFS[[i]] <- coxph(Surv(DFS,DFS.Event) ~ peptide[,names[i]], data=peptide)
  models.DFS.pr[[i]] <- coxph(Surv(DFS,DFS.Event) ~ peptide[,names[i]]+peptide[,"PR.Protein"], data=peptide)
  models.DFS.erpr[[i]] <- coxph(Surv(DFS,DFS.Event) ~ peptide[,names[i]]+peptide[,"PR.Protein"]+peptide[,"ER.Protein"], data=peptide)
  models.DFS.erpos[[i]] <- coxph(Surv(DFS[erpos],DFS.Event[erpos]) ~ peptide[erpos,names[i]], data=peptide)
  models.DFS.erneg[[i]] <- coxph(Surv(DFS[erneg],DFS.Event[erneg]) ~ peptide[erneg,names[i]], data=peptide)
  models.DFS.prpos[[i]] <- coxph(Surv(DFS[prpos],DFS.Event[prpos]) ~ peptide[prpos,names[i]], data=peptide)
  models.DFS.prneg[[i]] <- coxph(Surv(DFS[prneg],DFS.Event[prneg]) ~ peptide[prneg,names[i]], data=peptide)
  models.OS[[i]] <- coxph(Surv(OS,OS.Event) ~ peptide[,names[i]], data=peptide)
  models.OS.pr[[i]] <- coxph(Surv(OS,OS.Event) ~ peptide[,names[i]]+peptide[,"PR.Protein"], data=peptide)
  models.OS.erpr[[i]] <- coxph(Surv(OS,OS.Event) ~ peptide[,names[i]]+peptide[,"PR.Protein"]+peptide[,"ER.Protein"], data=peptide)
  models.OS.erpos[[i]] <- coxph(Surv(OS[erpos],OS.Event[erpos]) ~ peptide[erpos,names[i]], data=peptide)
  models.OS.erneg[[i]] <- coxph(Surv(OS[erneg],OS.Event[erneg]) ~ peptide[erneg,names[i]], data=peptide)
  models.OS.prpos[[i]] <- coxph(Surv(OS[prpos],OS.Event[prpos]) ~ peptide[prpos,names[i]], data=peptide)
  models.OS.prneg[[i]] <- coxph(Surv(OS[prneg],OS.Event[prneg]) ~ peptide[prneg,names[i]], data=peptide)
}

fx <- function(x) {
  sfit <- summary(x)
  return(sfit$coef[1,])
}

res.dfs <- sapply(models.DFS, fx)
res.adj <- p.adjust(round(as.matrix(t(res.dfs)[,5]),4),method= "BH")
res.dfs.1 <- cbind(t(res.dfs),res.adj)

res.dfs.pr <- sapply(models.DFS.pr, fx)
res.adj.pr <- p.adjust(round(as.matrix(t(res.dfs.pr)[,5]),4),method= "BH")
res.dfs.1.pr <- cbind(t(res.dfs.pr),res.adj.pr)

res.dfs.erpr <- sapply(models.DFS.erpr, fx)
res.adj.erpr <- p.adjust(round(as.matrix(t(res.dfs.erpr)[,5]),4),method= "BH")
res.dfs.1.erpr <- cbind(t(res.dfs.erpr),res.adj.erpr)

res.dfs.erpos <- sapply(models.DFS.erpos, fx)
res.adj.erpos <- p.adjust(round(as.matrix(t(res.dfs.erpos)[,5]),4),method= "BH")
res.dfs.1.erpos <- cbind(t(res.dfs.erpos),res.adj.erpos)

res.dfs.erneg <- sapply(models.DFS.erneg, fx)
res.adj.erneg <- p.adjust(round(as.matrix(t(res.dfs.erneg)[,5]),4),method= "BH")
res.dfs.1.erneg <- cbind(t(res.dfs.erneg),res.adj.erneg)

res.dfs.prpos <- sapply(models.DFS.prpos, fx)
res.adj.prpos <- p.adjust(round(as.matrix(t(res.dfs.prpos)[,5]),4),method= "BH")
res.dfs.1.prpos <- cbind(t(res.dfs.prpos),res.adj.prpos)

res.dfs.prneg <- sapply(models.DFS.prneg, fx)
res.adj.prneg <- p.adjust(round(as.matrix(t(res.dfs.prneg)[,5]),4),method= "BH")
res.dfs.1.prneg <- cbind(t(res.dfs.prneg),res.adj.prneg)

res.os <- sapply(models.OS, fx)
res.adj <- p.adjust(round(as.matrix(t(res.os)[,5]),4),method= "BH")
res.os.1 <- cbind(t(res.os),res.adj)

res.os.pr <- sapply(models.OS.pr, fx)
res.adj.pr <- p.adjust(round(as.matrix(t(res.os.pr)[,5]),4),method= "BH")
res.os.1.pr <- cbind(t(res.os.pr),res.adj.pr)

res.os.erpr <- sapply(models.OS.erpr, fx)
res.adj.erpr <- p.adjust(round(as.matrix(t(res.os.erpr)[,5]),4),method= "BH")
res.os.1.erpr <- cbind(t(res.os.erpr),res.adj.erpr)

res.os.erpos <- sapply(models.OS.erpos, fx)
res.adj.erpos <- p.adjust(round(as.matrix(t(res.os.erpos)[,5]),4),method= "BH")
res.os.1.erpos <- cbind(t(res.os.erpos),res.adj.erpos)

res.os.erneg <- sapply(models.OS.erneg, fx)
res.adj.erneg <- p.adjust(round(as.matrix(t(res.os.erneg)[,5]),4),method= "BH")
res.os.1.erneg <- cbind(t(res.os.erneg),res.adj.erneg)

res.os.prpos <- sapply(models.OS.prpos, fx)
res.adj.prpos <- p.adjust(round(as.matrix(t(res.os.prpos)[,5]),4),method= "BH")
res.os.1.prpos <- cbind(t(res.os.prpos),res.adj.prpos)

res.os.prneg <- sapply(models.OS.prneg, fx)
res.adj.prneg <- p.adjust(round(as.matrix(t(res.os.prneg)[,5]),4),method= "BH")
res.os.1.prneg <- cbind(t(res.os.prneg),res.adj.prneg)


'setwd("/Volumes/STORE N GO/MANUSCRIPTS/PR PAPER/PR peptide")'
setwd("MANUSCRIPTS/PR PAPER/PR peptide")
############ List significant univariate cox regression
### Here I identified the genes which are significant at less than 6% pvalues.
a <- res.dfs.1[,6]<1
a1 <- res.dfs.1[a,] 
k <- order(a1[,6])
z <- data.frame(a1[k,])
lowerCI <- round(apply(z,1,function(x) exp(x[1]-1.96*x[3])),2)
upperCI <- round(apply(z,1,function(x) exp(x[1]+1.96*x[3])),2)
CI <- t(t(paste0("(",lowerCI,",",upperCI,")")))
uni.dfs.adj <- cbind(round(z[,1:2],2),CI,round(z[,5:6],6))
write.csv(uni.dfs.adj, "uni.dfs.adj.csv")

a <- res.dfs.1.pr[,6]<1
a1 <- res.dfs.1.pr[a,] 
k <- order(a1[,6])
z <- data.frame(a1[k,])
lowerCI <- round(apply(z,1,function(x) exp(x[1]-1.96*x[3])),2)
upperCI <- round(apply(z,1,function(x) exp(x[1]+1.96*x[3])),2)
CI <- t(t(paste0("(",lowerCI,",",upperCI,")")))
uni.dfs.pr.adj <- cbind(round(z[,1:2],2),CI,round(z[,5:6],6))
write.csv(uni.dfs.pr.adj, "uni.dfs.pr.adj.csv")

a <- res.dfs.1.erpr[,6]<1
a1 <- res.dfs.1.erpr[a,] 
k <- order(a1[,6])
z <- data.frame(a1[k,])
lowerCI <- round(apply(z,1,function(x) exp(x[1]-1.96*x[3])),2)
upperCI <- round(apply(z,1,function(x) exp(x[1]+1.96*x[3])),2)
CI <- t(t(paste0("(",lowerCI,",",upperCI,")")))
uni.dfs.erpr.adj <- cbind(round(z[,1:2],2),CI,round(z[,5:6],6))
write.csv(uni.dfs.erpr.adj, "uni.dfs.erpr.adj.csv")

a <- res.dfs.1.erpos[,6]<1
a1 <- res.dfs.1.erpos[a,] 
k <- order(a1[,6])
z <- data.frame(a1[k,])
lowerCI <- round(apply(z,1,function(x) exp(x[1]-1.96*x[3])),2)
upperCI <- round(apply(z,1,function(x) exp(x[1]+1.96*x[3])),2)
CI <- t(t(paste0("(",lowerCI,",",upperCI,")")))
uni.dfs.adj.erpos <- cbind(round(z[,1:2],2),CI,round(z[,5:6],6))
write.csv(uni.dfs.adj.erpos, "uni.dfs.adj.erpos.csv")

a <- res.dfs.1.erneg[,6]<1
a1 <- res.dfs.1.erneg[a,] 
k <- order(a1[,6])
z <- data.frame(a1[k,])
lowerCI <- round(apply(z,1,function(x) exp(x[1]-1.96*x[3])),2)
upperCI <- round(apply(z,1,function(x) exp(x[1]+1.96*x[3])),2)
CI <- t(t(paste0("(",lowerCI,",",upperCI,")")))
uni.dfs.adj.erneg <- cbind(round(z[,1:2],2),CI,round(z[,5:6],6))
write.csv(uni.dfs.adj.erneg, "uni.dfs.adj.erneg.csv")

a <- res.dfs.1.prpos[,6]<1
a1 <- res.dfs.1.prpos[a,] 
k <- order(a1[,6])
z <- data.frame(a1[k,])
lowerCI <- round(apply(z,1,function(x) exp(x[1]-1.96*x[3])),2)
upperCI <- round(apply(z,1,function(x) exp(x[1]+1.96*x[3])),2)
CI <- t(t(paste0("(",lowerCI,",",upperCI,")")))
uni.dfs.adj.prpos <- cbind(round(z[,1:2],2),CI,round(z[,5:6],6))
write.csv(uni.dfs.adj.prpos, "uni.dfs.adj.prpos.csv")

a <- res.dfs.1.prneg[,6]<1
a1 <- res.dfs.1.prneg[a,] 
k <- order(a1[,6])
z <- data.frame(a1[k,])
lowerCI <- round(apply(z,1,function(x) exp(x[1]-1.96*x[3])),2)
upperCI <- round(apply(z,1,function(x) exp(x[1]+1.96*x[3])),2)
CI <- t(t(paste0("(",lowerCI,",",upperCI,")")))
uni.dfs.adj.prneg <- cbind(round(z[,1:2],2),CI,round(z[,5:6],6))
write.csv(uni.dfs.adj.prneg, "uni.dfs.adj.prneg.csv")

a <- res.os.1[,6]<1
a1 <- res.os.1[a,] 
k <- order(a1[,6])
z <- data.frame(a1[k,])
lowerCI <- round(apply(z,1,function(x) exp(x[1]-1.96*x[3])),2)
upperCI <- round(apply(z,1,function(x) exp(x[1]+1.96*x[3])),2)
CI <- t(t(paste0("(",lowerCI,",",upperCI,")")))
uni.os.adj <- cbind(round(z[,1:2],2),CI,round(z[,5:6],6))
write.csv(uni.os.adj, "uni.os.adj.csv")

a <- res.os.1.pr[,6]<1
a1 <- res.os.1.pr[a,] 
k <- order(a1[,6])
z <- data.frame(a1[k,])
lowerCI <- round(apply(z,1,function(x) exp(x[1]-1.96*x[3])),2)
upperCI <- round(apply(z,1,function(x) exp(x[1]+1.96*x[3])),2)
CI <- t(t(paste0("(",lowerCI,",",upperCI,")")))
uni.os.pr.adj <- cbind(round(z[,1:2],2),CI,round(z[,5:6],6))
write.csv(uni.os.pr.adj, "uni.os.pr.adj.csv")

a <- res.os.1.erpr[,6]<1
a1 <- res.os.1.erpr[a,] 
k <- order(a1[,6])
z <- data.frame(a1[k,])
lowerCI <- round(apply(z,1,function(x) exp(x[1]-1.96*x[3])),2)
upperCI <- round(apply(z,1,function(x) exp(x[1]+1.96*x[3])),2)
CI <- t(t(paste0("(",lowerCI,",",upperCI,")")))
uni.os.erpr.adj <- cbind(round(z[,1:2],2),CI,round(z[,5:6],6))
write.csv(uni.os.erpr.adj, "uni.os.erpr.adj.csv")

a <- res.os.1.erpos[,6]<1
a1 <- res.os.1.erpos[a,] 
k <- order(a1[,6])
z <- data.frame(a1[k,])
lowerCI <- round(apply(z,1,function(x) exp(x[1]-1.96*x[3])),2)
upperCI <- round(apply(z,1,function(x) exp(x[1]+1.96*x[3])),2)
CI <- t(t(paste0("(",lowerCI,",",upperCI,")")))
uni.os.adj.erpos <- cbind(round(z[,1:2],2),CI,round(z[,5:6],6))
write.csv(uni.os.adj.erpos, "uni.os.adj.erpos.csv")

a <- res.os.1.erneg[,6]<1
a1 <- res.os.1.erneg[a,] 
k <- order(a1[,6])
z <- data.frame(a1[k,])
lowerCI <- round(apply(z,1,function(x) exp(x[1]-1.96*x[3])),2)
upperCI <- round(apply(z,1,function(x) exp(x[1]+1.96*x[3])),2)
CI <- t(t(paste0("(",lowerCI,",",upperCI,")")))
uni.os.adj.erneg <- cbind(round(z[,1:2],2),CI,round(z[,5:6],6))
write.csv(uni.os.adj.erneg, "uni.os.adj.erneg.csv")

a <- res.os.1.prpos[,6]<1
a1 <- res.os.1.prpos[a,] 
k <- order(a1[,6])
z <- data.frame(a1[k,])
lowerCI <- round(apply(z,1,function(x) exp(x[1]-1.96*x[3])),2)
upperCI <- round(apply(z,1,function(x) exp(x[1]+1.96*x[3])),2)
CI <- t(t(paste0("(",lowerCI,",",upperCI,")")))
uni.os.adj.prpos <- cbind(round(z[,1:2],2),CI,round(z[,5:6],6))
write.csv(uni.os.adj.prpos, "uni.os.adj.prpos.csv")

a <- res.os.1.prneg[,6]<1
a1 <- res.os.1.prneg[a,] 
k <- order(a1[,6])
z <- data.frame(a1[k,])
lowerCI <- round(apply(z,1,function(x) exp(x[1]-1.96*x[3])),2)
upperCI <- round(apply(z,1,function(x) exp(x[1]+1.96*x[3])),2)
CI <- t(t(paste0("(",lowerCI,",",upperCI,")")))
uni.os.adj.prneg <- cbind(round(z[,1:2],2),CI,round(z[,5:6],6))
write.csv(uni.os.adj.prneg, "uni.os.adj.prneg.csv")

############ KM of significant univariate genes
a2 <- res.dfs.1.prpos[,6]<.3
a3.prpos <- res.dfs.1.prpos[a2,]
which(res.dfs.1.prpos[,6]<.3) # TMSB10
a2 <- res.dfs.1.prneg[,6]<.3
a3.prneg <- res.dfs.1.prneg[a2,]
b2 <- res.os.1.prpos[,6]<.3
b3.prpos <- res.os.1.prpos[b2,]
b2 <- res.os.1.prneg[,6]<.3
b3.prneg <- res.os.1.prneg[b2,]
dim(a3.prpos) # 1 6 There are 3 genes significant for DFS PR+ at the FDR cutoff of 0.5
dim(a3.prneg) # 8 6 There are 3 genes significant for DFS PR- at the FDR cutoff of 0.5
dim(b3.prpos) # 0 6 There are 3 genes significant for DFS PR- at the FDR cutoff of 0.5
dim(b3.prneg) # 3 6 There are 3 genes significant for DFS PR- at the FDR cutoff of 0.5

names.dfs.prpos <- "TMSB10"
names.dfs.prneg <- rownames(a3.prneg)
names.os.prpos <- rownames(b3.prpos) # none
names.os.prneg <- rownames(b3.prneg)

models.km.DFS.prpos <- vector("list", length(names.dfs.prpos))
models.km.DFS.prneg <- vector("list", length(names.dfs.prneg))
models.km.OS.prpos <- vector("list", length(names.os.prpos))
models.km.OS.prneg <- vector("list", length(names.os.prneg))

for (i in 1:length(names.dfs.prpos)){
  models.km.DFS.prpos[[i]] <- ifelse(peptide[prpos,names.dfs.prpos[i]]> median(peptide[prpos,names.dfs.prpos[i]]), 1, 0)
}

for (i in 1:length(names.dfs.prneg)){
  models.km.DFS.prneg[[i]] <- ifelse(peptide[prneg,names.dfs.prneg[i]]> median(peptide[prneg,names.dfs.prneg[i]]), 1, 0)
}

for (i in 1:length(names.os.prpos)){
  models.km.OS.prpos[[i]] <- ifelse(peptide[prpos,names.os.prpos[i]]> median(peptide[prpos,names.os.prpos[i]]), 1, 0)
}

for (i in 1:length(names.os.prneg)){
  models.km.OS.prneg[[i]] <- ifelse(peptide[prneg,names.os.prneg[i]]> median(peptide[prneg,names.os.prneg[i]]), 1, 0)
}

dfs.km.prpos <- vector("list", length(names.dfs.prpos))
a.km.prpos <- vector("list", length(names.dfs.prpos))

dfs.km.prneg <- vector("list", length(names.dfs.prneg))
a.km.prneg <- vector("list", length(names.dfs.prneg))

os.km.prpos <- vector("list", length(names.os.prpos))
b.km.prpos <- vector("list", length(names.os.prpos))

os.km.prneg <- vector("list", length(names.os.prneg))
b.km.prneg <- vector("list", length(names.os.prneg))

for(i in 1:length(names.dfs.prpos)){
  dfs.km.prpos[[i]] <- summary(coxph(Surv(peptide$DFS[prpos], peptide$DFS.Event[prpos]) ~ models.km.DFS.prpos[[i]]))
  a.km.prpos[[i]] <- survdiff(Surv(DFS[prpos], DFS.Event[prpos]) ~ models.km.DFS.prpos[[i]],data=peptide)  
}

dfs.km.p.prpos <- numeric(length(names.dfs.prpos))
for(i in 1:length(names.dfs.prpos)){
  dfs.km.p.prpos[i] <- pchisq(a.km.prpos[[i]]$chisq,1,lower=F)
}
p.dfs.km.prpos <- p.adjust(dfs.km.p.prpos,method="BH")

for(i in 1:length(names.dfs.prneg)){
  dfs.km.prneg[[i]] <- summary(coxph(Surv(peptide$DFS[prneg], peptide$DFS.Event[prneg]) ~ models.km.DFS.prneg[[i]]))
  a.km.prneg[[i]] <- survdiff(Surv(DFS[prneg], DFS.Event[prneg]) ~ models.km.DFS.prneg[[i]],data=peptide)  
}

dfs.km.p.prneg <- numeric(length(names.dfs.prneg))
for(i in 1:length(names.dfs.prneg)){
  dfs.km.p.prneg[i] <- pchisq(a.km.prneg[[i]]$chisq,1,lower=F)
}
p.dfs.km.prneg <- p.adjust(dfs.km.p.prneg,method="BH")


for(i in 1:length(names.os.prpos)){
  os.km.prpos[[i]] <- summary(coxph(Surv(peptide$OS[prpos], peptide$OS.Event[prpos]) ~ models.km.OS.prpos[[i]]))
  b.km.prpos[[i]] <- survdiff(Surv(OS[prpos], OS.Event[prpos]) ~ models.km.OS.prpos[[i]],data=peptide)  
}

os.km.p.prpos <- numeric(length(names.os.prpos))
for(i in 1:length(names.os.prpos)){
  os.km.p.prpos[i] <- pchisq(b.km.prpos[[i]]$chisq,1,lower=F)
}
p.os.km.prpos <- p.adjust(os.km.p.prpos,method="BH")

for(i in 1:length(names.os.prneg)){
  os.km.prneg[[i]] <- summary(coxph(Surv(peptide$OS[prneg], peptide$OS.Event[prneg]) ~ models.km.OS.prneg[[i]]))
  b.km.prneg[[i]] <- survdiff(Surv(OS[prneg], OS.Event[prneg]) ~ models.km.OS.prneg[[i]],data=peptide)  
}

os.km.p.prneg <- numeric(length(names.os.prneg))
for(i in 1:length(names.os.prneg)){
  os.km.p.prneg[i] <- pchisq(b.km.prneg[[i]]$chisq,1,lower=F)
}
p.os.km.prneg <- p.adjust(os.km.p.prneg,method="BH")

setwd("/Volumes/STORE N GO")
source("R/KM plot function.R")
setwd("MANUSCRIPTS/PR PAPER/PR peptide")
# for(i in 1:length(names.dfs.prpos)){
  # if(p.dfs.km.prpos[i]<0.3){
    # jpeg(paste0(names.dfs.prpos[i]," DFS PR+ KM.jpg"))
jpeg("TMSB10 DFS PR+ KM.jpg")    
km <- kmplot(survfit(Surv(DFS[prpos], DFS.Event[prpos]) ~ models.km.DFS.prpos[[1]], data = peptide ), mark='', simple=FALSE,
                 xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
                 lty.surv=1, lwd.surv=2, col.surv=c(1,2), font.lab=2, font.axis=2, # survival.curves
                 col.ci=0, # confidence intervals not plotted
                 group.names=c('Below','Above'),
                 group.order=c(1,2), # order of appearance in the n.risk.at table and legend.
                 extra.left.margin=0, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
                 grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
                 legend=FALSE, loc.legend='bottomleft',
                 cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
    ) 
    text(0,0.05,adj=0,cex=1.5,font=2,paste0(names.dfs.prpos[i],"   adj. p = ",round(p.dfs.km.prpos[i],5)," PR+"))
    mtext("Progression Free Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
    mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
    print(km)
    dev.off()
  


for(i in 1:length(names.dfs.prneg)){
  if(p.dfs.km.prneg[i]<0.3){
    jpeg(paste0(names.dfs.prneg[i]," DFS PR- KM.jpg"))
    km <- kmplot(survfit(Surv(DFS[prneg], DFS.Event[prneg]) ~ models.km.DFS.prneg[[i]], data = peptide ), mark='', simple=FALSE,
                 xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
                 lty.surv=1, lwd.surv=2, col.surv=c(1,2), font.lab=2, font.axis=2, # survival.curves
                 col.ci=0, # confidence intervals not plotted
                 group.names=c('Below','Above'),
                 group.order=c(1,2), # order of appearance in the n.risk.at table and legend.
                 extra.left.margin=0, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
                 grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
                 legend=FALSE, loc.legend='bottomleft',
                 cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
    ) 
    text(0,0.05,adj=0,cex=1.5,font=2,paste0(names.dfs.prneg[i],"   adj. p = ",round(p.dfs.km.prneg[i],5)," PR-"))
    mtext("Progression Free Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
    mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
    print(km)
    dev.off()
  }
}

for(i in 1:length(names.os.prpos)){
  if(p.os.km.prpos[i]<0.3){
    jpeg(paste0(names.os.prpos[i]," OS PR+ KM.jpg"))
    km <- kmplot(survfit(Surv(OS[prpos], OS.Event[prpos]) ~ models.km.OS.prpos[[i]], data = peptide ), mark='', simple=FALSE,
                 xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
                 lty.surv=1, lwd.surv=2, col.surv=c(1,2), font.lab=2, font.axis=2, # survival.curves
                 col.ci=0, # confidence intervals not plotted
                 group.names=c('Below','Above'),
                 group.order=c(1,2), # order of appearance in the n.risk.at table and legend.
                 extra.left.margin=0, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
                 grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
                 legend=FALSE, loc.legend='bottomleft',
                 cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
    ) 
    text(0,0.05,adj=0,cex=1.5,font=2,paste0(names.os.prpos[i],"   adj. p = ",round(p.os.km.prpos[i],5)," PR+"))
    mtext("Overall Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
    mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
    print(km)
    dev.off()
  }
}

for(i in 1:length(names.os.prneg)){
  if(p.os.km.prneg[i]<0.3){
    jpeg(paste0(names.os.prneg[i]," OS PR- KM.jpg"))
    km <- kmplot(survfit(Surv(OS[prneg], OS.Event[prneg]) ~ models.km.OS.prneg[[i]], data = peptide ), mark='', simple=FALSE,
                 xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
                 lty.surv=1, lwd.surv=2, col.surv=c(1,2), font.lab=2, font.axis=2, # survival.curves
                 col.ci=0, # confidence intervals not plotted
                 group.names=c('Below','Above'),
                 group.order=c(1,2), # order of appearance in the n.risk.at table and legend.
                 extra.left.margin=0, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
                 grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
                 legend=FALSE, loc.legend='bottomleft',
                 cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
    ) 
    text(0,0.05,adj=0,cex=1.5,font=2,paste0(names.os.prneg[i],"   adj. p = ",round(p.os.km.prneg[i],5)," PR-"))
    mtext("Overall Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
    mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
    print(km)
    dev.off()
  }
}

##### Molecular Signature
library(penalized)
library(SGL)
data.exp <- peptide[,2:143]

#### Group membership
mydata <- scale(data.exp[prpos,])
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

library(cluster) 
fit2 <- kmeans(mydata, 2) # 2 cluster solution
fit3 <- kmeans(mydata, 3) # 3 cluster solution
fit4 <- kmeans(mydata, 4) # 4 cluster solution
fit5 <- kmeans(mydata, 5) # 5 cluster solution
clusplot(mydata, fit2$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)
clusplot(mydata, fit3$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)
clusplot(mydata, fit4$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)
clusplot(mydata, fit5$cluster, color=TRUE, shade=TRUE, 
         labels=2, lines=0)

# get cluster means 
center2 <- aggregate(mydata,by=list(fit2$cluster),FUN=mean)
center3 <- aggregate(mydata,by=list(fit3$cluster),FUN=mean)
center4 <- aggregate(mydata,by=list(fit4$cluster),FUN=mean)
center5 <- aggregate(mydata,by=list(fit5$cluster),FUN=mean)

# append cluster assignment
mydata2 <- data.frame(mydata, fit2$cluster)
mydata3 <- data.frame(mydata, fit3$cluster)
mydata4 <- data.frame(mydata, fit4$cluster)
mydata5 <- data.frame(mydata, fit5$cluster)
mydata.surv <- data.frame(peptide[prpos,c(1,144:151)], fit2$cluster, fit3$cluster, fit4$cluster, fit5$cluster)

d2 <- dist(mydata2, method = "euclidean") # distance matrix
d3 <- dist(mydata3, method = "euclidean") # distance matrix
d4 <- dist(mydata4, method = "euclidean") # distance matrix
d5 <- dist(mydata5, method = "euclidean") # distance matrix
fit2.ward <- hclust(d2, method="ward.D") 
fit3.ward <- hclust(d3, method="ward.D") 
fit4.ward <- hclust(d4, method="ward.D") 
fit5.ward <- hclust(d5, method="ward.D") 

plot(fit2.ward) # display dendogram
groups2 <- cutree(fit2.ward, k=2) # cut tree into 2 clusters
groups3 <- cutree(fit3.ward, k=3) # cut tree into 2 clusters
groups4 <- cutree(fit4.ward, k=4) # cut tree into 2 clusters
groups5 <- cutree(fit5.ward, k=5) # cut tree into 2 clusters

# draw dendogram with red borders around the 2 clusters 
rect.hclust(fit2.ward, k=2, border="red")

plot(fit3.ward) # display dendogram
groups2 <- cutree(fit3.ward, k=3) # cut tree into 5 clusters
# draw dendogram with red borders around the 3 clusters 
rect.hclust(fit3.ward, k=3, border="red")

plot(fit4.ward) # display dendogram
groups4 <- cutree(fit4.ward, k=4) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters 
rect.hclust(fit4.ward, k=4, border="red")

plot(fit5.ward) # display dendogram
groups5 <- cutree(fit5.ward, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit5.ward, k=5, border="red")

# Centroid Plot against 1st 2 discriminant functions
library(fpc)
plotcluster(mydata, fit2$cluster)
plotcluster(mydata, fit3$cluster)
plotcluster(mydata, fit4$cluster)
plotcluster(mydata, fit5$cluster)

cluster.stats(d2, fit2$cluster)
cluster.stats(d3, fit3$cluster)
cluster.stats(d4, fit4$cluster)
cluster.stats(d5, fit5$cluster)

library(ggfortify)
library(survival)
library(ggplot2)
# DFS
fit2.dfs <- survfit(Surv(DFS, DFS.Event) ~ fit2.cluster, data = mydata.surv)
survdiff(Surv(DFS, DFS.Event) ~ fit2.cluster, data = mydata.surv)
autoplot(fit2.dfs,main='PFS 2-Grouping')

fit3.dfs <- survfit(Surv(DFS, DFS.Event) ~ fit3.cluster, data = mydata.surv)
survdiff(Surv(DFS, DFS.Event) ~ fit3.cluster, data = mydata.surv)
autoplot(fit3.dfs,main='PFS 3-Grouping')

fit4.dfs <- survfit(Surv(DFS, DFS.Event) ~ fit4.cluster, data = mydata.surv)
survdiff(Surv(DFS, DFS.Event) ~ fit4.cluster, data = mydata.surv)
autoplot(fit4.dfs,main='PFS 4-Grouping')

fit5.dfs <- survfit(Surv(DFS, DFS.Event) ~ fit5.cluster, data = mydata.surv)
survdiff(Surv(DFS, DFS.Event) ~ fit5.cluster, data = mydata.surv)
autoplot(fit5.dfs,main='PFS 5-Grouping')
# 3 groups have the most significant survival outcome for DFS

# OS
fit2.os <- survfit(Surv(OS, OS.Event) ~ fit2.cluster, data = mydata.surv)
survdiff(Surv(OS, OS.Event) ~ fit2.cluster, data = mydata.surv)
autoplot(fit2.os,main='OS 2-Grouping')

fit3.os <- survfit(Surv(OS, OS.Event) ~ fit3.cluster, data = mydata.surv)
survdiff(Surv(OS, OS.Event) ~ fit3.cluster, data = mydata.surv)
autoplot(fit3.os,main='OS 3-Grouping')

fit4.os <- survfit(Surv(OS, OS.Event) ~ fit4.cluster, data = mydata.surv)
survdiff(Surv(OS, OS.Event) ~ fit4.cluster, data = mydata.surv)
autoplot(fit4.os,main='OS 4-Grouping')

fit5.os <- survfit(Surv(OS, OS.Event) ~ fit5.cluster, data = mydata.surv)
survdiff(Surv(OS, OS.Event) ~ fit5.cluster, data = mydata.surv)
autoplot(fit5.os,main='OS 5-Grouping')

# 4 groups have the most significant survival outcome for OS

clin <- read.csv("~/Desktop/Publications/PR PAPER/clin.csv")
idx.clin.prpos <- clin[,1]%in%peptide[prpos,1]
clin.prpos <- clin[idx.clin.prpos,]
names(clin.prpos)

library(tidyr)
clin.prpos.iddesc <- clin.prpos %>% arrange(desc(Case))
peptide.prpos.iddesc <- peptide[prpos,] %>% arrange(desc(Case.ID))
cluster.prpos.iddesc <- mydata.surv %>% arrange(desc(Case.ID)) 
head(clin.prpos.iddesc$Case)
head(peptide.prpos.iddesc$Case.ID)
tail(clin.prpos.iddesc$Case)
tail(peptide.prpos.iddesc$Case.ID)
clin.pep <- cbind(peptide.prpos.iddesc,clin.prpos.iddesc,cluster.prpos.iddesc)
names(clin.pep)
idx.3group.1 <- which(clin.pep$fit3.cluster==1)
idx.3group.2 <- which(clin.pep$fit3.cluster==2)
idx.3group.3 <- which(clin.pep$fit3.cluster==3)
idx.4group.1 <- which(clin.pep$fit4.cluster==1)
idx.4group.2 <- which(clin.pep$fit4.cluster==2)
idx.4group.3 <- which(clin.pep$fit4.cluster==3)
idx.4group.4 <- which(clin.pep$fit4.cluster==4)
table.3.1 <- apply(clin.pep[idx.3group.1,c(155,156,157,159,161,163,164,167)],2,table)
means.3.1 <- apply(clin.pep[idx.3group.1,c(153,158,170,172)],2,function(x) mean(x,na.rm=TRUE))
table.3.2 <- apply(clin.pep[idx.3group.2,c(155,156,157,159,161,163,164,167)],2,table)
means.3.2 <- apply(clin.pep[idx.3group.2,c(153,158,170,172)],2,function(x) mean(x,na.rm=TRUE))
table.3.3 <- apply(clin.pep[idx.3group.3,c(155,156,157,159,161,163,164,167)],2,table)
means.3.3 <- apply(clin.pep[idx.3group.3,c(153,158,170,172)],2,function(x) mean(x,na.rm=TRUE))
table.4.1 <- apply(clin.pep[idx.4group.1,c(155,156,157,159,161,163,164,167)],2,table)
means.4.1 <- apply(clin.pep[idx.4group.1,c(153,158,170,172)],2,function(x) mean(x,na.rm=TRUE))
table.4.2 <- apply(clin.pep[idx.4group.2,c(155,156,157,159,161,163,164,167)],2,table)
means.4.2 <- apply(clin.pep[idx.4group.2,c(153,158,170,172)],2,function(x) mean(x,na.rm=TRUE))
table.4.3 <- apply(clin.pep[idx.4group.3,c(155,156,157,159,161,163,164,167)],2,table)
means.4.3 <- apply(clin.pep[idx.4group.3,c(153,158,170,172)],2,function(x) mean(x,na.rm=TRUE))
table.4.4 <- apply(clin.pep[idx.4group.4,c(155,156,157,159,161,163,164,167)],2,table)
means.4.4 <- apply(clin.pep[idx.4group.4,c(153,158,170,172)],2,function(x) mean(x,na.rm=TRUE))

rbind(means.3.1,means.3.2,means.3.3)
rbind(table.3.1$Race,table.3.2$Race,table.3.3$Race)
rbind(table.3.1$Fam.Hx,table.3.2$Fam.Hx,table.3.3$Fam.Hx)
rbind(table.3.1$Tobacco,table.3.2$Tobacco,table.3.3$Tobacco)
rbind(table.3.1$Menop.status,table.3.2$Menop.status,table.3.3$Menop.status)
rbind(table.3.1$Pathology,table.3.2$Pathology,table.3.3$Pathology)
rbind(table.3.1$Node.Status,c(table.3.2$Node.Status,0),c(table.3.3$Node.Status,0))
rbind(c(table.3.1$Grade,'4'=0),c(table.3.2$Grade,'4'=0),table.3.3$Grade)
rbind(c(table.3.1$Stage[1],'0'=0,table.3.1$Stage[2],'2'=0,table.3.1$Stage[3:7],'4'=0),c(table.3.2$Stage[1:6],'3'=0,'3A'=0,table.3.2$Stage[7:8]),table.3.3$Stage)

rbind(means.4.1,means.4.2,means.4.3,means.4.4)
rbind(table.4.1$Race,table.4.2$Race,table.4.3$Race,table.4.4$Race)
rbind(table.4.1$Fam.Hx,table.4.2$Fam.Hx,table.4.3$Fam.Hx,table.4.4$Fam.Hx)
rbind(table.4.1$Tobacco,table.4.2$Tobacco,table.4.3$Tobacco,table.4.4$Tobacco)
rbind(table.4.1$Menop.status,table.4.2$Menop.status,table.4.3$Menop.status,table.4.4$Menop.status)
rbind(table.4.1$Pathology,table.4.2$Pathology,table.4.3$Pathology,table.4.4$Pathology)
rbind(table.3.1$Node.Status,c(table.3.2$Node.Status,0),c(table.3.3$Node.Status,0))
rbind(c('1'=0,table.4.1$Grade,'4'=0),table.4.2$Grade,c(table.4.3$Grade,'4'=0),c(table.4.4$Grade,'4'=0))
rbind(c(table.4.1$Stage[1],'0'=0,table.4.1$Stage[2],'2'=0,table.4.1$Stage[3:5],'3A'=0,'3B'=0,'4'=0),table.4.2$Stage,c(table.4.3$Stage[1:6],'3'=0,'3A'=0,table.4.3$Stage[7:8]),c(table.4.4$Stage[1:3],'2'=0,table.4.4$Stage[4:5],'3'=0,table.4.4$Stage[6:7],'4'=0))

# LASSO OS
data.sgl.os.prpos <- list(x=data.exp[prpos,], time = peptide$OS[prpos], status = peptide$OS.Event[prpos])
data.sgl.os.prneg <- list(x=data.exp[prneg,], time = peptide$OS[prneg], status = peptide$OS.Event[prneg])
sgl.index <- rep(1,142)
fit.sgl.os.prpos <- SGL(data.sgl.os.prpos, sgl.index, type = "cox")
fit.sgl.os.prneg <- SGL(data.sgl.os.prneg, sgl.index, type = "cox")

# try penalized command
optL1.os.prpos <- optL1(Surv(peptide$OS[prpos],peptide$OS.Event[prpos]), penalized = peptide[prpos,2:143]) # 20.672
fit.pen.os.prpos <- penalized(Surv(peptide$OS[prpos],peptide$OS.Event[prpos]), penalized = peptide[prpos,2:143],lambda1 = optL1.os.prpos$lambda)
optL1.os.prneg <- optL1(Surv(peptide$OS[prneg],peptide$OS.Event[prneg]), penalized = peptide[prneg,2:143]) # 22.39418
fit.pen.os.prneg <- penalized(Surv(peptide$OS[prneg],peptide$OS.Event[prneg]), penalized = peptide[prneg,2:143],lambda1 = optL1.os.prpos$lambda)
optL1.os <- optL1(Surv(peptide$OS,peptide$OS.Event), penalized = peptide[,2:143],unpenalized = peptide[,149]) # 22.39418
fit.pen.os <- penalized(Surv(peptide$OS,peptide$OS.Event), penalized = peptide[,2:143],unpenalized = peptide[,149],lambda1 = optL1.os$lambda)

write.csv(coef(fit.pen.os.prpos),"pen.os.prpos.csv")
write.csv(coef(fit.pen.os.prneg),"pen.os.prneg.csv")
write.csv(coef(fit.pen.os),"pen.os.csv")

# LASSO DFS
data.sgl.dfs.prpos <- list(x=data.exp[prpos,], time = peptide$DFS[prpos], status = peptide$DFS.Event[prpos])
data.sgl.dfs.prneg <- list(x=data.exp[prneg,], time = peptide$DFS[prneg], status = peptide$DFS.Event[prneg])
sgl.index <- rep(1,142)
fit.sgl.dfs.prpos <- SGL(data.sgl.dfs.prpos, sgl.index, type = "cox")
fit.sgl.dfs.prneg <- SGL(data.sgl.dfs.prneg, sgl.index, type = "cox")

# try penalized command
optL1.dfs.prpos <- optL1(Surv(peptide$DFS[prpos],peptide$DFS.Event[prpos]), penalized = peptide[prpos,2:143]) # 20.672
fit.pen.dfs.prpos <- penalized(Surv(peptide$DFS[prpos],peptide$DFS.Event[prpos]), penalized = peptide[prpos,2:143],lambda1 = optL1.dfs.prpos$lambda)
optL1.dfs.prneg <- optL1(Surv(peptide$DFS[prneg],peptide$DFS.Event[prneg]), penalized = peptide[prneg,2:143]) # 22.39418
fit.pen.dfs.prneg <- penalized(Surv(peptide$DFS[prneg],peptide$DFS.Event[prneg]), penalized = peptide[prneg,2:143],lambda1 = optL1.dfs.prpos$lambda)
optL1.dfs <- optL1(Surv(peptide$DFS,peptide$DFS.Event), penalized = peptide[,2:143],unpenalized = peptide[,149]) # 22.39418
fit.pen.dfs <- penalized(Surv(peptide$DFS,peptide$DFS.Event), penalized = peptide[,2:143],unpenalized = peptide[,149],lambda1 = optL1.os$lambda)

write.csv(coef(fit.pen.dfs.prpos),"pen.dfs.prpos.csv")
write.csv(coef(fit.pen.dfs.prneg),"pen.dfs.prneg.csv")
write.csv(coef(fit.pen.dfs),"pen.dfs.csv")

# cross validation LASSO OS
cv.sgl.os.prpos <- cvSGL(data.sgl.os.prpos, sgl.index, type = "cox")
cv.sgl.os.prneg <- cvSGL(data.sgl.os.prneg, sgl.index, type = "cox")

which.min(cv.sgl.os.prpos$lldiff) # 2
which.min(cv.sgl.os.prneg$lldiff) # 2

names <- colnames(data.exp)

os.prpos <- which(fit.sgl.os.prpos$beta[,2]!=0)
os.prpos.lasso <- data.frame(cbind(names[os.prpos],round(fit.sgl.os.prpos$beta[,2][os.prpos],5)))
os.prneg <- which(fit.sgl.os.prneg$beta[,2]!=0)
os.prneg.lasso <- data.frame(cbind(names[os.prneg],round(fit.sgl.os.prneg$beta[,2][os.prneg],5)))

colnames(os.erpos.lasso) <- colnames(os.erneg.lasso) <- c("Gene Symbol","Beta Coefficient")
plot(cv.sgl.os.prpos)


# cross validation LASSO DFS
cv.sgl.dfs.prpos <- cvSGL(data.sgl.dfs.prpos, sgl.index, type = "cox")
cv.sgl.dfs.prneg <- cvSGL(data.sgl.dfs.prneg, sgl.index, type = "cox")

which.min(cv.sgl.dfs.prpos$lldiff) # 1
which.min(cv.sgl.dfs.prneg$lldiff) # 3

dfs.prpos <- which(fit.sgl.os.prpos$beta[,2]!=0)
dfs.prpos.lasso <- data.frame(cbind(names[os.prpos],round(fit.sgl.os.prpos$beta[,2][os.prpos],5)))
dfs.prneg <- which(fit.sgl.os.prneg$beta[,3]!=0)
dfs.prneg.lasso <- data.frame(cbind(names[os.prneg],round(fit.sgl.os.prneg$beta[,3][os.prneg],5)))

plot(cv.sgl.dfs.prpos)
plot(cv.sgl.dfs.prneg)

names <- colnames(data.exp)
dfs.prpos <- which(fit.sgl.dfs.prpos$beta[,2]!=0)
dfs.prpos.lasso <- data.frame(cbind(names[dfs.prpos],round(fit.sgl.dfs.prpos$beta[,2][dfs.prpos],5)))
dfs.prneg <- which(fit.sgl.dfs.prneg$beta[,4]!=0)


genes <- os.prpos.lasso[,1]
betas <- os.prpos.lasso[,2]

idx <- match(genes,colnames(peptide))
colnames(peptide)[idx]
lin.pred <- as.numeric(t(as.matrix(betas)))%*%t(as.matrix(peptide[prpos,idx]))
median.split <- ifelse(lin.pred<=median(lin.pred),0,1)
coxph(Surv(OS[prpos], OS.Event[prpos]) ~ as.numeric(median.split), data = peptide)

for(i in 1:length(names.dfs.prpos)){
  jpeg("OS PR+ KM signature.jpg")
    km <- kmplot(survfit(Surv(OS[prpos], OS.Event[prpos]) ~ median.split, data = peptide ), mark='', simple=FALSE,
                 xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
                 lty.surv=1, lwd.surv=2, col.surv=c(1,2), font.lab=2, font.axis=2, # survival.curves
                 col.ci=0, # confidence intervals not plotted
                 group.names=c('Below','Above'),
                 group.order=c(1,2), # order of appearance in the n.risk.at table and legend.
                 extra.left.margin=0, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
                 grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
                 legend=FALSE, loc.legend='bottomleft',
                 cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
    ) 
    text(0,0.05,adj=0,cex=1.5,font=2,paste0(names.dfs.prpos[i],"   adj. p = ",round(p.dfs.km.prpos[i],5)," PR+"))
    mtext("Overall Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
    mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
    print(km)
    dev.off()
  }

for(i in 1:length(names.dfs.prneg)){
  if(p.dfs.km.prneg[i]<0.3){
    jpeg(paste0(names.dfs.prneg[i]," DFS PR- KM.jpg"))
    km <- kmplot(survfit(Surv(DFS[prneg], DFS.Event[prneg]) ~ models.km.DFS.prneg[[i]], data = peptide ), mark='', simple=FALSE,
                 xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
                 lty.surv=1, lwd.surv=2, col.surv=c(1,2), font.lab=2, font.axis=2, # survival.curves
                 col.ci=0, # confidence intervals not plotted
                 group.names=c('Below','Above'),
                 group.order=c(1,2), # order of appearance in the n.risk.at table and legend.
                 extra.left.margin=0, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
                 grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
                 legend=FALSE, loc.legend='bottomleft',
                 cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
    ) 
    text(0,0.05,adj=0,cex=1.5,font=2,paste0(names.dfs.prneg[i],"   adj. p = ",round(p.dfs.km.prneg[i],5)," PR-"))
    mtext("Progression Free Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
    mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
    print(km)
    dev.off()
  }
}

for(i in 1:length(names.os.prpos)){
  if(p.os.km.prpos[i]<0.3){
    jpeg(paste0(names.os.prpos[i]," OS PR+ KM.jpg"))
    km <- kmplot(survfit(Surv(OS[prpos], OS.Event[prpos]) ~ models.km.OS.prpos[[i]], data = peptide ), mark='', simple=FALSE,
                 xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
                 lty.surv=1, lwd.surv=2, col.surv=c(1,2), font.lab=2, font.axis=2, # survival.curves
                 col.ci=0, # confidence intervals not plotted
                 group.names=c('Below','Above'),
                 group.order=c(1,2), # order of appearance in the n.risk.at table and legend.
                 extra.left.margin=0, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
                 grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
                 legend=FALSE, loc.legend='bottomleft',
                 cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
    ) 
    text(0,0.05,adj=0,cex=1.5,font=2,paste0(names.os.prpos[i],"   adj. p = ",round(p.os.km.prpos[i],5)," PR+"))
    mtext("Overall Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
    mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
    print(km)
    dev.off()
  }
}

for(i in 1:length(names.os.prneg)){
  if(p.os.km.prneg[i]<0.3){
    jpeg(paste0(names.os.prneg[i]," OS PR- KM.jpg"))
    km <- kmplot(survfit(Surv(OS[prneg], OS.Event[prneg]) ~ models.km.OS.prneg[[i]], data = peptide ), mark='', simple=FALSE,
                 xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
                 lty.surv=1, lwd.surv=2, col.surv=c(1,2), font.lab=2, font.axis=2, # survival.curves
                 col.ci=0, # confidence intervals not plotted
                 group.names=c('Below','Above'),
                 group.order=c(1,2), # order of appearance in the n.risk.at table and legend.
                 extra.left.margin=0, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
                 grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
                 legend=FALSE, loc.legend='bottomleft',
                 cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
    ) 
    text(0,0.05,adj=0,cex=1.5,font=2,paste0(names.os.prneg[i],"   adj. p = ",round(p.os.km.prneg[i],5)," PR-"))
    mtext("Overall Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
    mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
    print(km)
    dev.off()
  }
}


save.image("prpeptide")
