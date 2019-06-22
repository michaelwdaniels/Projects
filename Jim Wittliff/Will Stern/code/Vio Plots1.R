
install.packages("packages")

par(mfrow=c(2,2)) 
metcaff <-  read.csv("~/Downloads/kmplotproblems/Total EIA LBA Extraction.csv")
head(metcaff)
dim(metcaff) # 1221 11
tail(metcaff)
data <- metcaff
head(data)
names(data)
head(data)

cox <- metcaff
genes <- names(data)[6:9]
l <- length(genes)
head(data)
dim(data)
attach(data)
data

ER.Status.EIA <- ifelse(ER_EIA_RESULT_VALUE<10,0,1) 
er.neg.eia <- which(ER.Status.EIA==0)
er.pos.eia <- which(ER.Status.EIA==1)
length(er.neg.eia) # 232
length(er.pos.eia) # 517

ER.Status.MTA <- ifelse(ER_MTA_SBC_RESULT_VALUE<15,0,1) 
er.neg.mta <- which(ER.Status.MTA==0)
er.pos.mta <- which(ER.Status.MTA==1)
length(er.neg.mta) # 204
length(er.pos.mta) # 343

PR.Status.EIA <- ifelse(PR_EIA_RESULT_VALUE<10,0,1) 
pr.neg.eia <- which(PR.Status.EIA==0)
pr.pos.eia <- which(PR.Status.EIA==1)
length(pr.neg.eia) # 240
length(pr.pos.eia) # 503

PR.Status.MTA <- ifelse(PR_MTA_SBC_RESULT_VALUE<15,0,1) 
pr.neg.mta <- which(PR.Status.MTA==0)
pr.pos.mta <- which(PR.Status.MTA==1)
length(pr.neg.mta) # 163
length(pr.pos.mta) # 370

data$Status.EIA <- interaction(ER.Status.EIA,PR.Status.EIA)
data$Status.MTA <- interaction(ER.Status.MTA,PR.Status.MTA)

model.er.eia <- vector("list",l)
model.pr.eia <- vector("list",l)
model.er.mta <- vector("list",l)
model.pr.mta <- vector("list",l)

for(i in 1:l){
  
  model.er.eia[[i]] <- t.test(data[er.pos.eia,genes[i]], data[er.neg.eia,genes[i]])
  model.pr.eia[[i]] <- t.test(data[pr.pos.eia,genes[i]], data[pr.neg.eia,genes[i]])
  
}  

for(i in 1:l){
  
  if(model.er[[i]]$p.value<0.06){
  if(model.er[[i]]$p.value<0.001){
  boxplot(data[er.pos,genes[i]], data[er.neg,genes[i]],
          main=paste0(" ",genes[i], " p < 0.001 "),
          ylab="Protein Status",
          names=c("ER+","ER-"))
  }
  else  
    boxplot(data[er.pos,genes[i]], data[er.neg,genes[i]],
            main=paste0(" ",genes[i], " p = ",round(model.er[[i]]$p.value,3)),
            ylab="Protein Status",
            names=c("ER+","ER-"))}}


for(i in 1:l){   

  if(model.pr[[i]]$p.value<0.06){
  if(model.pr[[i]]$p.value<0.001){
  boxplot(data[pr.pos,genes[i]], data[pr.neg,genes[i]],
            main=paste0(" ",genes[i], " p < 0.001 "),
            ylab="Protein Status",
            names=c("PR+","PR-")) 
  }
  else
    boxplot(data[pr.pos,genes[i]], data[pr.neg,genes[i]],
            main=paste0(" ",genes[i], " p = ",round(model.pr[[i]]$p.value,3)),
            ylab="Protein Status",
            names=c("PR+","PR-"))}}





########### BW and Violin Plots
par(mfrow=c(1,1))
graphics.off()
cox <- metcaff
head(cox)
dim(cox)
names <- names(cox[2:3])
sig.names.erbox <- vector("list",length(names))
sig.names.prbox <- vector("list",length(names))
model.er <- vector("list",length(names))
model.pr <- vector("list",length(names))
names(model.er) <- names(model.pr) <- names
erpos <- which(cox$ER.Status==1)
erneg <- which(cox$ER.Status==0)
prpos <- which(cox$PR.Status==1)
prneg <- which(cox$PR.Status==0)
length(erpos) # 48
length(erneg) # 48
length(prpos) # 50
length(prneg) # 46
for(i in 1:length(names)){
  gene.er.pos <- cox[,names[i]][erpos]
  gene.er.neg <- cox[,names[i]][erneg]
  model.er[[i]] <- wilcox.test(gene.er.pos,gene.er.neg,exact=FALSE)
  gene.pr.pos <- cox[,names[i]][prpos]
  gene.pr.neg <- cox[,names[i]][prneg]
  model.pr[[i]] <- wilcox.test(gene.pr.pos,gene.pr.neg,exact=FALSE)
}

box.er.p <- sapply(model.er, function(x) x$p.value)
box.pr.p <- sapply(model.pr, function(x) x$p.value)
box.er.padj <- p.adjust(box.er.p,method="BH")
box.pr.padj <- p.adjust(box.pr.p,method="BH")
box.er.padj.1 <- box.er.padj[which(box.er.padj<0.901)]
box.er.names <- names(box.er.padj.1)
box.er.padj.2 <- box.er.padj[which(box.er.padj>0.001&box.er.padj<0.3001)]
box.er.names.1 <- names(box.er.padj.2)
box.pr.padj.1 <- box.pr.padj[which(box.pr.padj<0.901)]
box.pr.names <- names(box.pr.padj.1)
box.pr.padj.2 <- box.pr.padj[which(box.pr.padj>0.001&box.pr.padj<0.3001)]
box.pr.names.1 <- names(box.pr.padj.2)
sort(box.er.padj)
sort(box.pr.padj)
length(box.er.names) # 4
length(box.er.names) # 2
length(box.pr.names.1) # 6
length(box.pr.names.1) # 6
library(vioplot)
#### ER Box and Whisker and violin plots for p < 0.001
par(mai=c(1,1.1,.1,.1))
for(i in 1:length(box.er.names)){
  jpeg(paste0(box.er.names[i]," BW ER.jpeg"))
  bw <- boxplot(cox[erpos,box.er.names[i]],cox[erneg,box.er.names[i]],lwd=2, font.axis=2,xaxt="no")
  legend("topright",adj=0,cex=1.3,paste0(box.er.names[i]," adj. p < 0.001"),text.font=2)
  mtext("Relative Gene Expression",side=2,font=2,line=2.25,cex=1.3)
  mtext(c("ER+","ER-"),side=1,at=c(1,2),line=1,font=2,cex=1.3)
  print(bw)
  dev.off()
}

par(mai=c(0.75,1,.1,.1))
for(i in 1:length(box.er.names)){
  jpeg(paste0(box.er.names[i]," VIOLIN ER.jpeg"))
  vio <- vioplot(cox[erpos,box.er.names[i]],cox[erneg,box.er.names[i]],col=8,lwd=2,names=c("",""))
  mtext(paste0(box.er.names[i]," adj. p < 0.001"),side=3,line=1,font=2,cex=1.9)
  axis(side=1,lwd=4,tck=0,labels = FALSE)
  axis(side=2,cex.axis=1.2,font=2,lwd=4,tck=0)
  axis(3,lwd=4,tck=0,labels = FALSE)
  axis(4,lwd=4,tck=0,labels = FALSE)
  mtext("Relative Gene Expression",side=2,font=2,line=2.5,cex=1.9)
  mtext(c("ER+","ER-"),side=1,at=c(1,2),line=1.8,font=2.9,cex=1.9)
  print(vio)
  dev.off()
}

#### PR Box and Whisker and violin plots for p < 0.001
par(mai=c(1,1.1,.1,.1))
for(i in 1:length(box.pr.names)){
  jpeg(paste0(box.pr.names[i]," BW PR.jpeg"))
  bw <- boxplot(cox[prpos,box.pr.names[i]],cox[prneg,box.pr.names[i]],lwd=2, font.axis=2,xaxt="no")
  legend("topright",adj=0,cex=1.3,paste0(box.pr.names[i]," adj. p < 0.001"),text.font=2)
  mtext("Relative Gene Expression",side=2,font=2,line=2.25,cex=1.3)
  mtext(c("PR+","PR-"),side=1,at=c(1,2),line=1,font=2,cex=1.3)
  print(bw)
  dev.off()
}

par(mai=c(0.75,1,.1,.1))
for(i in 1:length(box.pr.names)){
  jpeg(paste0(box.pr.names[i]," VIOLIN PR.jpeg"))
  vio <- vioplot(cox[prpos,box.pr.names[i]],cox[prneg,box.pr.names[i]],col=8,lwd=2,names=c("",""))
  mtext(paste0(box.pr.names[i]," adj. p = 0.003"),side=3,line=1,font=2,cex=1.9)
  mtext("Relative Gene Expression",side=2,font=2,line=2.5,cex=1.9)
  axis(side=1,lwd=4,tck=0,labels = FALSE)
  axis(side=2,cex.axis=1.2,font=2,lwd=4,tck=0)
  axis(3,lwd=4,tck=0,labels = FALSE)
  axis(4,lwd=4,tck=0,labels = FALSE)
  mtext(c("PR+","PR-"),side=1,at=c(1,2),line=1.8,font=2.9,cex=1.9)
  print(vio)
  dev.off()
}

#### ER Box and Whisker and violin plots for 0.001 < p < 0.01
par(mai=c(1,1.1,.1,.1))
for(i in 1:length(box.er.names.1)){
  jpeg(paste0(box.er.names.1[i]," BW ER.jpeg"))
  bw <- boxplot(cox[erpos,box.er.names.1[i]],cox[erneg,box.er.names.1[i]],lwd=2, font.axis=2,xaxt="no")
  legend("topright",adj=0,cex=1.25,paste0(box.er.names.1[i]," adj. p = ",round(box.er.padj.2[i],3)),text.font=2)
  mtext("Relative Gene Expression",side=2,font=2,line=2.25,cex=1.3)
  mtext(c("ER+","ER-"),side=1,at=c(1,2),line=1,font=2,cex=1.3)
  print(bw)
  dev.off()
}

par(mai=c(0.75,1,.1,.1))
for(i in 1:length(box.er.names.1)){
  jpeg(paste0(box.er.names.1[i]," VIOLIN ER.jpeg"))
  vio <- vioplot(cox[erpos,box.er.names.1[i]],cox[erneg,box.er.names.1[i]],col=8,lwd=2,names=c("",""))
  mtext(paste0(box.er.names[i]," adj. p < 0.001"),side=3,line=1,font=2,cex=1.7)
  axis(side=1,lwd=2,tck=0,labels = FALSE)
  axis(side=2,lwd=2)
  axis(3,lwd=2,tck=0,labels = FALSE)
  axis(4,lwd=2,tck=0,labels = FALSE)
  mtext("Relative Gene Expression",side=2,font=2,line=2.5,cex=1.7)
  mtext(c("Post","Pre"),side=1,at=c(1,2),line=1.8,font=2.9,cex=1.7)
  print(vio)
  dev.off()
}

#### PR Box and Whisker and violin plots for 0.001 < p < 0.01
par(mai=c(1,1.1,.1,.1))
for(i in 1:length(box.pr.names.1)){
  jpeg(paste0(box.pr.names.1[i]," BW PR.jpeg"))
  bw <- boxplot(cox[prpos,box.pr.names.1[i]],cox[prneg,box.pr.names.1[i]],lwd=2, font.axis=2,xaxt="no")
  legend("topright",adj=0,cex=1.3,paste0(box.pr.names.1[i]," adj. p = ",round(box.pr.padj.2[i],3)),text.font=2)
  mtext("Relative Gene Expression",side=2,font=3,line=2.25,cex=1.5)
  mtext(c("Post","Pre"),side=1,at=c(1,2),line=1,font=2,cex=1.5)
  print(bw)
  dev.off()
}

par(mai=c(0.75,1,.1,.1))
for(i in 1:length(box.pr.names.1)){
  jpeg(paste0(box.pr.names.1[i]," VIOLIN PR.jpeg"))
  vio <- vioplot(cox[prpos,box.pr.names.1[i]],cox[prneg,box.pr.names.1[i]],col=8,lwd=2,names=c("",""))
  mtext(paste0(box.pr.names[i]," adj. p < 0.059, n=247"),side=3,line=1,font=2,cex=1.7)
  axis(side=1,lwd=2,tck=0,labels = FALSE)
  axis(side=2,lwd=2)
  axis(3,lwd=2,tck=0,labels = FALSE)
  axis(4,lwd=2,tck=0,labels = FALSE)
  mtext("Relative Gene Expression",side=2,font=2,line=2.5,cex=1.9)
  mtext(c("Post n=143","Pre n=102"),side=1,at=c(1,2),line=1.8,font=2.9,cex=1.7)
  mtext(c("Menopausal Status"),side=1,line=3.5,font=2.9,cex=1.7)
  print(vio)
  dev.off()
}

