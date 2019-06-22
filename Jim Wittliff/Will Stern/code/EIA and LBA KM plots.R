## KM Plots
# Load necessary packages
setwd("~/Desktop/Projects/Will Stern")
library(survival)

# Load necessary functions
source('KM plot function.R')
Total.EIA.LBA.Extraction <- read.csv("Total EIA LBA Extraction.csv")

# EIA and LBA ER/PR PFS
DFS.cut <- ifelse(Total.EIA.LBA.Extraction$DFS>144,144,Total.EIA.LBA.Extraction$DFS)
Total.EIA.LBA.Extraction$EIAcat <- ifelse(Total.EIA.LBA.Extraction$ER_EIA_RESULT_VALUE<=14&Total.EIA.LBA.Extraction$PR_EIA_RESULT_VALUE<=14, 1,
                                           ifelse(Total.EIA.LBA.Extraction$ER_EIA_RESULT_VALUE<=14&Total.EIA.LBA.Extraction$PR_EIA_RESULT_VALUE>14, 2,
                                                  ifelse(Total.EIA.LBA.Extraction$ER_EIA_RESULT_VALUE>14&Total.EIA.LBA.Extraction$PR_EIA_RESULT_VALUE<=14, 3, 4)))
Total.EIA.LBA.Extraction$LBAcat <- ifelse(Total.EIA.LBA.Extraction$ER_MTA_SBC_RESULT_VALUE<=9&Total.EIA.LBA.Extraction$PR_MTA_SBC_RESULT_VALUE<=9, 1,
                                          ifelse(Total.EIA.LBA.Extraction$ER_MTA_SBC_RESULT_VALUE<=9&Total.EIA.LBA.Extraction$PR_MTA_SBC_RESULT_VALUE>9, 2,
                                                 ifelse(Total.EIA.LBA.Extraction$ER_MTA_SBC_RESULT_VALUE>9&Total.EIA.LBA.Extraction$PR_MTA_SBC_RESULT_VALUE<=9, 3, 4)))
table(Total.EIA.LBA.Extraction$EIAcat)
table(Total.EIA.LBA.Extraction$LBAcat)

EIA.LBAcomb <- cbind(Total.EIA.LBA.Extraction$EIAcat, Total.EIA.LBA.Extraction$LBAcat)
table(EIA.LBAcomb)
head(EIA.LBAcomb)
EIA.LBAcomb1 <- cbind(EIA.LBAcomb, ifelse(is.na(EIA.LBAcomb[,1]), EIA.LBAcomb[,2], EIA.LBAcomb[,1]))
table(EIA.LBAcomb1[,3])
sum(table(EIA.LBAcomb1[,3])) #1180
Total.EIA.LBA.Extraction$cat <-EIA.LBAcomb1[,3]

p.dfs <- vector('numeric',1) #change for number of genes
moddiff <- survdiff(Surv(Total.EIA.LBA.Extraction$DFS, Total.EIA.LBA.Extraction$DFS.Event) ~ Total.EIA.LBA.Extraction$cat)

# FOR EXTRA GENES USE BELOW

p.dfs[1] <- pchisq(moddiff$chisq,3,lower.tail = F) # 3 degrees of freemdom due to 4 groups
#moddiff <- survdiff(Surv(DFS, DFS.Event) ~ statusUparEr, data = upa)
#p.dfs[2] <- pchisq(moddiff$chisq,3,lower.tail = F)
#moddiff <- survdiff(Surv(DFS, DFS.Event) ~ statusPai1Er, data = upa)
#p.dfs[3] <- pchisq(moddiff$chisq,3,lower.tail = F)
#p.adj <- p.adjust(p.dfs,method = "fdr")
# summary(EIA.Extraction$DFS[which(EIA.Extraction$ERcat==1)])

jpeg('LBA_ER_DFS_KMplot.jpg')    
km <- kmplot(survfit(Surv(DFS.cut, DFS.Event) ~ EIA.LBAcomb1[,3], data = Total.EIA.LBA.Extraction ), mark='', simple=FALSE,
             xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
             lty.surv=1, lwd.surv=2, col.surv=c(1,2,3,4), font.lab=2, font.axis=2, # survival.curves
             col.ci=0, # confidence intervals not plotted
             group.names=c('ER-/PR-', 'ER-/PR+', 'ER+/PR-', 'ER+/PR+'),
             group.order=c(1,2,3,4), # order of appearance in the n.risk.at table and legend.
             extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
             grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
             legend=FALSE, loc.legend='bottomleft',
             cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
)
text(0,0.05,adj=0,cex=1.5,font=2,paste0('EIA and LBA ER/PR  p = ',signif(p.dfs[1],3)))
mtext('Progression Free Survival (mos)',side=1,line=1.9, font=2, cex=1.5)
mtext('Survival Probability',side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()


# EIA LBA ER/PR (OS)
OS.cut <- ifelse(Total.EIA.LBA.Extraction$OS>144,144,Total.EIA.LBA.Extraction$OS)
Total.EIA.LBA.Extraction$ERPRcat <- ifelse(Total.EIA.LBA.Extraction$ER_EIA_RESULT_VALUE<=14&Total.EIA.LBA.Extraction$PR_EIA_RESULT_VALUE<=14&Total.EIA.LBA.Extraction$ER_MTA_SBC_RESULT_VALUE<=9&Total.EIA.LBA.Extraction$PR_MTA_SBC_RESULT_VALUE<=9, 1,
                                           ifelse(Total.EIA.LBA.Extraction$ER_EIA_RESULT_VALUE<=14&Total.EIA.LBA.Extraction$PR_EIA_RESULT_VALUE>14&Total.EIA.LBA.Extraction$ER_MTA_SBC_RESULT_VALUE<=9&Total.EIA.LBA.Extraction$PR_MTA_SBC_RESULT_VALUE>9, 2,
                                                  ifelse(Total.EIA.LBA.Extraction$ER_EIA_RESULT_VALUE>14&Total.EIA.LBA.Extraction$PR_EIA_RESULT_VALUE<=14&Total.EIA.LBA.Extraction$ER_MTA_SBC_RESULT_VALUE>9&Total.EIA.LBA.Extraction$PR_MTA_SBC_RESULT_VALUE<=9, 3,
                                                         ifelse(Total.EIA.LBA.Extraction$ER_EIA_RESULT_VALUE>14&Total.EIA.LBA.Extraction$PR_EIA_RESULT_VALUE>14&Total.EIA.LBA.Extraction$ER_MTA_SBC_RESULT_VALUE>9&Total.EIA.LBA.Extraction$PR_MTA_SBC_RESULT_VALUE>9, 4, NA))))
p.dfs <- vector('numeric',1) #change for number of genes
moddiff <- survdiff(Surv(Total.EIA.LBA.Extraction$OS, Total.EIA.LBA.Extraction$OS.Event) ~ Total.EIA.LBA.Extraction$cat)

# FOR EXTRA GENES USE BELOW

p.dfs[1] <- pchisq(moddiff$chisq,3,lower.tail = F) # 3 degrees of freemdom due to 4 groups
#moddiff <- survdiff(Surv(DFS, DFS.Event) ~ statusUparEr, data = upa)
#p.dfs[2] <- pchisq(moddiff$chisq,3,lower.tail = F)
#moddiff <- survdiff(Surv(DFS, DFS.Event) ~ statusPai1Er, data = upa)
#p.dfs[3] <- pchisq(moddiff$chisq,3,lower.tail = F)
#p.adj <- p.adjust(p.dfs,method = "fdr")
#summary(EIA.Extraction$DFS[which(EIA.Extraction$ERcat==1)])

jpeg('LBA_ER_OS_KMplot.jpg')    
km <- kmplot(survfit(Surv(OS.cut, OS.Event) ~ EIA.LBAcomb1[,3], data = Total.EIA.LBA.Extraction ), mark='', simple=FALSE,
             xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
             lty.surv=1, lwd.surv=2, col.surv=c(1,2,3,4), font.lab=2, font.axis=2, # survival.curves
             col.ci=0, # confidence intervals not plotted
             group.names=c('ER-/PR-', 'ER-/PR+', 'ER+/PR-', 'ER+/PR+'),
             group.order=c(1,2,3,4), # order of appearance in the n.risk.at table and legend.
             extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
             grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
             legend=FALSE, loc.legend='bottomleft',
             cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
)
text(0,0.05,adj=0,cex=1.5,font=2,paste0('EIA and LBA ER/PR  p = ',signif(p.dfs[1],3)))
mtext('Overall Survival (mos)',side=1,line=1.9, font=2, cex=1.5)
mtext('Survival Probability',side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()





ER_EIA_status <- ifelse(Total.EIA.LBA.Extraction$ER_EIA_RESULT_VALUE<15, 0, 1)
ER_LBA_status <- ifelse(Total.EIA.LBA.Extraction$ER_MTA_SBC_RESULT_VALUE<10, 0, 1)
PR_EIA_status <- ifelse(Total.EIA.LBA.Extraction$PR_EIA_RESULT_VALUE<15, 0, 1)
PR_LBA_status <- ifelse(Total.EIA.LBA.Extraction$PR_MTA_SBC_RESULT_VALUE<10, 0, 1)
table(interaction(ER_EIA_status,ER_LBA_status,PR_EIA_status,PR_LBA_status))
table(interaction(ER_EIA_status, ER_LBA_status))
