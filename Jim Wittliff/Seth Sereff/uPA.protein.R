### uPA protein and ER status interaction on clinical outcomes

# 1) Since ER protein status of breast cancer patients may be divided into positive 
# (good prognosis) and negative (poor prognosis)...does the inclusion of uPA protein 
# content of the biopsies exhibit different progression-free survival times? My thought 
# is that you could calculate the median of the progression free survival and compare 
# these values for the four groups, i.e....ER+/uPA+, ER+/uPA, ER-/uPA+, ER-/uPA-. 
# Repeat for uPAR and PAI-1.  We need to answer whether the addition of the analyses 
# of either uPA, uPAR or PAI-1 protein content to the ER protein content of a biopsy 
# is clinically useful in predicting recurrence free survival.  Your thoughts?

# 2) Repeat the about calculations using PR protein content with uPA, uPAR and PAI-1 
# protein content.

# 3) Repeat the about calculations using HER2 protein content with uPA, uPAR and PAI-1 
# protein content. Addendum: There is not enough HER-2 data to perform analysis

## Analysis Plan

# Does the prescence of ER protein affect the relationship between uPA, uPAR or PAI-1
# protein content and progression free and overall survival?

# Models

# (PFS, PFS Event) ~ B0 + B1*uPA + B2*ER + B3*uPA*ER
# (PFS, PFS Event) ~ uPAR + ER + uPAR*ER
# (PFS, PFS Event) ~ PAI-1 + ER + PAI-1*ER
# (OS, OS Event) ~ uPA + ER + uPA*ER
# (OS, OS Event) ~ uPAR + ER + uPAR*ER
# (OS, OS Event) ~ PAI-1 + ER + PAI-1*ER

# repeat by replacing ER with PR and HER2 for a total of 18 models
# significance of interaction term determines the synergistic/antagonist effects
# of the prescence of ER

# Kaplan Meier plots

# 12 plots could be displayed 
# (3 genes {uPA, uPAR, PAI-1} x 2 outcomes {DFS, OS} x 2 hormone status {ER, PR}) 
# HR Status {ER+,ER-} x gene median split {uPA+,uPA-} (4 groups)

# Set working directory
setwd("~/Desktop/Projects/uPA")

# Read in data
library(readr) 
upa <- read_csv("uPA_protein_final.csv")

# Load necessary packages
library(survival)

# Load necessary functions
source('KM plot function.R')

# Create unique identifier for the four subgroups for ER 
upa$statusUpaEr <- interaction(upa$upa.status,upa$ER.Status)
upa$statusUparEr <- interaction(upa$upar.status,upa$ER.Status)
upa$statusPai1Er <- interaction(upa$pai1.status,upa$ER.Status)

# Create unique identifier for the four subgroups for PR 
upa$statusUpaPr <- interaction(upa$upa.status,upa$PR.Status)
upa$statusUparPr <- interaction(upa$upar.status,upa$PR.Status)
upa$statusPai1Pr <- interaction(upa$pai1.status,upa$PR.Status)

# Create 12 KM plots ({ER PR}x{DFS OS}x{uPA uPAR PAI-1}) 

## DFS

### ER
### Log-rank test for survival curve differences - aggregate p-values from three genes
### and adjust using Holm (default in p.adjust function)
p.dfs <- vector('numeric',3)
moddiff <- survdiff(Surv(DFS, DFS.Event) ~ statusUpaEr, data = upa)
p.dfs[1] <- pchisq(moddiff$chisq,3,lower.tail = F) # 3 degrees of freemdom due to 4 groups
moddiff <- survdiff(Surv(DFS, DFS.Event) ~ statusUparEr, data = upa)
p.dfs[2] <- pchisq(moddiff$chisq,3,lower.tail = F)
moddiff <- survdiff(Surv(DFS, DFS.Event) ~ statusPai1Er, data = upa)
p.dfs[3] <- pchisq(moddiff$chisq,3,lower.tail = F)
p.adj <- p.adjust(p.dfs,method = "fdr")

jpeg("UpaErDfs.jpg")    
km <- kmplot(survfit(Surv(DFS, DFS.Event) ~ statusUpaEr, data = upa ), mark='', simple=FALSE,
             xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
             lty.surv=1, lwd.surv=2, col.surv=c(1,2,3,4), font.lab=2, font.axis=2, # survival.curves
             col.ci=0, # confidence intervals not plotted
             group.names=c('uPA - ER -','uPA + ER -','uPA - ER +','uPA + ER +'),
             group.order=c(1,2,3,4), # order of appearance in the n.risk.at table and legend.
             extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
             grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
             legend=FALSE, loc.legend='bottomleft',
             cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("uPA ER PFS  adj. p = ",round(p.adj[1],3)))
mtext("Progression Free Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

jpeg("UparErDfs.jpg")    
km <- kmplot(survfit(Surv(DFS, DFS.Event) ~ statusUparEr, data = upa ), mark='', simple=FALSE,
             xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
             lty.surv=1, lwd.surv=2, col.surv=c(1,2,3,4), font.lab=2, font.axis=2, # survival.curves
             col.ci=0, # confidence intervals not plotted
             group.names=c('uPAR - ER -','uPAR + ER -','uPAR - ER +','uPAR + ER +'),
             group.order=c(1,2,3,4), # order of appearance in the n.risk.at table and legend.
             extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
             grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
             legend=FALSE, loc.legend='bottomleft',
             cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("uPAR ER PFS  adj. p = ",round(p.adj[2],3)))
mtext("Progression Free Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

jpeg("Pai1ErDfs.jpg")    
km <- kmplot(survfit(Surv(DFS, DFS.Event) ~ statusPai1Er, data = upa ), mark='', simple=FALSE,
             xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
             lty.surv=1, lwd.surv=2, col.surv=c(1,2,3,4), font.lab=2, font.axis=2, # survival.curves
             col.ci=0, # confidence intervals not plotted
             group.names=c('PAI-1 - ER -','PAI-1 + ER -','PAI-1 - ER +','PAI-1 + ER +'),
             group.order=c(1,2,3,4), # order of appearance in the n.risk.at table and legend.
             extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
             grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
             legend=FALSE, loc.legend='bottomleft',
             cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("PAI-1 ER PFS  adj. p = ",round(p.adj[3],3)))
mtext("Progression Free Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

### PR

p.dfs <- vector('numeric',3)
moddiff <- survdiff(Surv(OS, OS.Event) ~ statusUpaPr, data = upa)
p.dfs[1] <- pchisq(moddiff$chisq,3,lower.tail = F)
moddiff <- survdiff(Surv(OS, DFS.Event) ~ statusUparPr, data = upa)
p.dfs[2] <- pchisq(moddiff$chisq,3,lower.tail = F)
moddiff <- survdiff(Surv(OS, DFS.Event) ~ statusPai1Pr, data = upa)
p.dfs[3] <- pchisq(moddiff$chisq,3,lower.tail = F)
p.adj <- p.adjust(p.dfs,method = 'fdr')

jpeg("UpaPrDfs.jpg")    
km <- kmplot(survfit(Surv(DFS, DFS.Event) ~ statusUpaPr, data = upa ), mark='', simple=FALSE,
             xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
             lty.surv=1, lwd.surv=2, col.surv=c(1,2,3,4), font.lab=2, font.axis=2, # survival.curves
             col.ci=0, # confidence intervals not plotted
             group.names=c('uPA - PR -','uPA + PR -','uPA - PR +','uPA + PR +'),
             group.order=c(1,2,3,4), # order of appearance in the n.risk.at table and legend.
             extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
             grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
             legend=FALSE, loc.legend='bottomleft',
             cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("uPA PR PFS  adj. p = ",round(p.adj[1],3)))
mtext("Progression Free Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

jpeg("UparPrDfs.jpg")    
km <- kmplot(survfit(Surv(DFS, DFS.Event) ~ statusUparPr, data = upa ), mark='', simple=FALSE,
             xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
             lty.surv=1, lwd.surv=2, col.surv=c(1,2,3,4), font.lab=2, font.axis=2, # survival.curves
             col.ci=0, # confidence intervals not plotted
             group.names=c('uPAR - PR -','uPAR + PR -','uPAR - PR +','uPAR + PR +'),
             group.order=c(1,2,3,4), # order of appearance in the n.risk.at table and legend.
             extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
             grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
             legend=FALSE, loc.legend='bottomleft',
             cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("uPAR PR PFS  adj. p = ",round(p.adj[2],3)))
mtext("Progression Free Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

jpeg("Pai1PrDfs.jpg")    
km <- kmplot(survfit(Surv(DFS, DFS.Event) ~ statusPai1Pr, data = upa ), mark='', simple=FALSE,
             xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
             lty.surv=1, lwd.surv=2, col.surv=c(1,2,3,4), font.lab=2, font.axis=2, # survival.curves
             col.ci=0, # confidence intervals not plotted
             group.names=c('PAI-1 - PR -','PAI-1 + PR -','PAI-1 - PR +','PAI-1 + PR +'),
             group.order=c(1,2,3,4), # order of appearance in the n.risk.at table and legend.
             extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
             grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
             legend=FALSE, loc.legend='bottomleft',
             cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("PAI-1 PR PFS  adj. p = ",round(p.adj[3],3)))
mtext("Progression Free Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

## OS
### ER
p.os <- vector('numeric',3)
moddiff <- survdiff(Surv(OS, OS.Event) ~ statusUpaEr, data = upa)
p.os[1] <- pchisq(moddiff$chisq,3,lower.tail = F)
moddiff <- survdiff(Surv(OS, OS.Event) ~ statusUparEr, data = upa)
p.os[2] <- pchisq(moddiff$chisq,3,lower.tail = F)
moddiff <- survdiff(Surv(OS, OS.Event) ~ statusPai1Er, data = upa)
p.os[3] <- pchisq(moddiff$chisq,3,lower.tail = F)
p.adj <- p.adjust(p.os,method = 'fdr')

jpeg("UpaErOs.jpg")    
km <- kmplot(survfit(Surv(OS, OS.Event) ~ statusUpaEr, data = upa ), mark='', simple=FALSE,
             xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
             lty.surv=1, lwd.surv=2, col.surv=c(1,2,3,4), font.lab=2, font.axis=2, # survival.curves
             col.ci=0, # confidence intervals not plotted
             group.names=c('uPA - ER -','uPA + ER -','uPA - ER +','uPA + ER +'),
             group.order=c(1,2,3,4), # order of appearance in the n.risk.at table and legend.
             extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
             grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
             legend=FALSE, loc.legend='bottomleft',
             cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("uPA ER OS  adj. p = ",round(p.adj[1],3)))
mtext("Overall Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

jpeg("UparErOs.jpg")    
km <- kmplot(survfit(Surv(OS, OS.Event) ~ statusUparEr, data = upa ), mark='', simple=FALSE,
             xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
             lty.surv=1, lwd.surv=2, col.surv=c(1,2,3,4), font.lab=2, font.axis=2, # survival.curves
             col.ci=0, # confidence intervals not plotted
             group.names=c('uPAR - ER -','uPAR + ER -','uPAR - ER +','uPAR + ER +'),
             group.order=c(1,2,3,4), # order of appearance in the n.risk.at table and legend.
             extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
             grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
             legend=FALSE, loc.legend='bottomleft',
             cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("uPAR ER OS  adj. p = ",round(p.adj[2],3)))
mtext("Overall Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

jpeg("Pai1ErOs.jpg")    
km <- kmplot(survfit(Surv(OS, OS.Event) ~ statusPai1Er, data = upa ), mark='', simple=FALSE,
             xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
             lty.surv=1, lwd.surv=2, col.surv=c(1,2,3,4), font.lab=2, font.axis=2, # survival.curves
             col.ci=0, # confidence intervals not plotted
             group.names=c('PAI-1 - ER -','PAI-1 + ER -','PAI-1 - ER +','PAI-1 + ER +'),
             group.order=c(1,2,3,4), # order of appearance in the n.risk.at table and legend.
             extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
             grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
             legend=FALSE, loc.legend='bottomleft',
             cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("PAI-1 ER OS  adj. p = ",round(p.adj[3],3)))
mtext("Overall Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

### PR

upa$statusUpaPr <- interaction(upa$upa.status,upa$PR.Status)
upa$statusUparPr <- interaction(upa$upar.status,upa$PR.Status)
upa$statusPai1Pr <- interaction(upa$pai1.status,upa$PR.Status)
p.os <- vector('numeric',3)
moddiff <- survdiff(Surv(OS, OS.Event) ~ statusUpaPr, data = upa)
p.os[1] <- pchisq(moddiff$chisq,3,lower.tail = F)
moddiff <- survdiff(Surv(OS, OS.Event) ~ statusUparPr, data = upa)
p.os[2] <- pchisq(moddiff$chisq,3,lower.tail = F)
moddiff <- survdiff(Surv(OS, OS.Event) ~ statusPai1Pr, data = upa)
p.os[3] <- pchisq(moddiff$chisq,3,lower.tail = F)
p.adj <- p.adjust(p.os,method = 'fdr')

jpeg("UpaPrOs.jpg")    
km <- kmplot(survfit(Surv(OS, OS.Event) ~ statusUpaPr, data = upa ), mark='', simple=FALSE,
             xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
             lty.surv=1, lwd.surv=2, col.surv=c(1,2,3,4), font.lab=2, font.axis=2, # survival.curves
             col.ci=0, # confidence intervals not plotted
             group.names=c('uPA - PR -','uPA + PR -','uPA - PR +','uPA + PR +'),
             group.order=c(1,2,3,4), # order of appearance in the n.risk.at table and legend.
             extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
             grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
             legend=FALSE, loc.legend='bottomleft',
             cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("uPA PR OS  adj. p = ",round(p.adj[1],3)))
mtext("Overall Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

jpeg("UparPrOs.jpg")    
km <- kmplot(survfit(Surv(OS, OS.Event) ~ statusUparPr, data = upa ), mark='', simple=FALSE,
             xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
             lty.surv=1, lwd.surv=2, col.surv=c(1,2,3,4), font.lab=2, font.axis=2, # survival.curves
             col.ci=0, # confidence intervals not plotted
             group.names=c('uPAR - PR -','uPAR + PR -','uPAR - PR +','uPAR + PR +'),
             group.order=c(1,2,3,4), # order of appearance in the n.risk.at table and legend.
             extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
             grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
             legend=FALSE, loc.legend='bottomleft',
             cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("uPAR PR OS  adj. p = ",round(p.adj[2],3)))
mtext("Overall Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

jpeg("Pai1PrOs.jpg")    
km <- kmplot(survfit(Surv(OS, OS.Event) ~ statusPai1Pr, data = upa ), mark='', simple=FALSE,
             xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
             lty.surv=1, lwd.surv=2, col.surv=c(1,2,3,4), font.lab=2, font.axis=2, # survival.curves
             col.ci=0, # confidence intervals not plotted
             group.names=c('PAI-1 - PR -','PAI-1 + PR -','PAI-1 - PR +','PAI-1 + PR +'),
             group.order=c(1,2,3,4), # order of appearance in the n.risk.at table and legend.
             extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
             grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
             legend=FALSE, loc.legend='bottomleft',
             cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("PAI-1 PR OS  adj. p = ",round(p.adj[3],3)))
mtext("Overall Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()
library(coxphf)
upa.dfs <- upa[-25,]
mod11 <- coxphf(Surv(DFS, DFS.Event) ~ statusUpaEr,data=upa.dfs)
summary(mod11)
coxphfplot(Surv(DFS, DFS.Event) ~ statusUpaEr,data=upa.dfs)

coxph(Surv(DFS, DFS.Event) ~ statusUparEr, data = upa)
coxph(Surv(DFS, DFS.Event) ~ statusPai1Er, data = upa)
coxph(Surv(DFS, DFS.Event) ~ statusUpaPr, data = upa)
coxph(Surv(DFS, DFS.Event) ~ statusUparPr, data = upa)
coxph(Surv(DFS, DFS.Event) ~ statusPai1Pr, data = upa)

coxph(Surv(OS, OS.Event) ~ statusUpaEr, data = upa)
coxph(Surv(OS, OS.Event) ~ statusUparEr, data = upa)
coxph(Surv(OS, OS.Event) ~ statusPai1Er, data = upa)
coxph(Surv(OS, OS.Event) ~ statusUpaPr, data = upa)
coxph(Surv(OS, OS.Event) ~ statusUparPr, data = upa)
coxph(Surv(OS, OS.Event) ~ statusPai1Pr, data = upa)

coxph(Surv(DFS, DFS.Event) ~ upa + interaction(PR.Status,ER.Status), data = upa)
coxph(Surv(DFS, DFS.Event) ~ interaction(PR.Status,ER.Status), data = upa)
coxph(Surv(DFS, DFS.Event) ~ upa + upar + pai1, data = upa[which(upa$ER.Status==1),])
coxph(Surv(DFS, DFS.Event) ~ upa + upar + pai1, data = upa[which(upa$PR.Status==0),])
coxph(Surv(DFS, DFS.Event) ~ upa + upar + pai1, data = upa[which(upa$PR.Status==1),])

coxph(Surv(OS, OS.Event) ~ upa + upar + pai1, data = upa[which(upa$ER.Status==0),])
coxph(Surv(OS, OS.Event) ~ upa + upar + pai1, data = upa[which(upa$ER.Status==1),])
coxph(Surv(OS, OS.Event) ~ upa + upar + pai1, data = upa[which(upa$PR.Status==0),])
coxph(Surv(OS, OS.Event) ~ upa + upar + pai1, data = upa[which(upa$PR.Status==1),])

# 2 level median gene split KM plots

p.dfs <- vector('numeric',3)
p.mod <- survdiff(Surv(DFS, DFS.Event) ~ upa.status, data = upa)
p.dfs[1] <- pchisq(p.mod$chisq,1,lower.tail = F)
p.mod <- survdiff(Surv(DFS, DFS.Event) ~ upar.status, data = upa)
p.dfs[2] <- pchisq(p.mod$chisq,1,lower.tail = F)
p.mod <- survdiff(Surv(DFS, DFS.Event) ~ pai1.status, data = upa)
p.dfs[3] <- pchisq(p.mod$chisq,1,lower.tail = F)
p.adj <- p.adjust(p.dfs,method = 'fdr')

jpeg("UpaDfs.jpg")
km <- kmplot(survfit(Surv(DFS, DFS.Event) ~ upa.status, data = upa), mark='', simple=FALSE,
             xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
             lty.surv=1, lwd.surv=2, col.surv=c(1,2), font.lab=2, font.axis=2, # survival.curves
             col.ci=0, # confidence intervals not plotted
             group.names=c('uPA -','uPA +'),
             group.order=c(1,2), # order of appearance in the n.risk.at table and legend.
             extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
             grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
             legend=FALSE, loc.legend='bottomleft',
             cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("uPA PFS  adj. p = ",round(p.adj[1],3)))
mtext("Progression Free Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

jpeg("UparDfs.jpg")
km <- kmplot(survfit(Surv(DFS, DFS.Event) ~ upar.status, data = upa), mark='', simple=FALSE,
       xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
       lty.surv=1, lwd.surv=2, col.surv=c(1,2), font.lab=2, font.axis=2, # survival.curves
       col.ci=0, # confidence intervals not plotted
       group.names=c('uPA -','uPA +'),
       group.order=c(1,2), # order of appearance in the n.risk.at table and legend.
       extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
       grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
       legend=FALSE, loc.legend='bottomleft',
       cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("uPAR PFS  adj. p = ",round(p.adj[2],3)))
mtext("Progression Free Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

jpeg("Pai1Dfs.jpg")
km <- kmplot(survfit(Surv(DFS, DFS.Event) ~ pai1.status, data = upa), mark='', simple=FALSE,
       xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
       lty.surv=1, lwd.surv=2, col.surv=c(1,2), font.lab=2, font.axis=2, # survival.curves
       col.ci=0, # confidence intervals not plotted
       group.names=c('uPA -','uPA +'),
       group.order=c(1,2), # order of appearance in the n.risk.at table and legend.
       extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
       grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
       legend=FALSE, loc.legend='bottomleft',
       cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("PAI-1 PFS  adj. p = ",round(p.adj[3],3)))
mtext("Progression Free Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

## OS
p.os <- vector('numeric',3)
p.mod <- survdiff(Surv(OS, OS.Event) ~ upa.status, data = upa)
p.os[1] <- pchisq(p.mod$chisq,1,lower.tail = F)
p.mod <- survdiff(Surv(OS, OS.Event) ~ upar.status, data = upa)
p.os[2] <- pchisq(p.mod$chisq,1,lower.tail = F)
p.mod <- survdiff(Surv(OS, OS.Event) ~ pai1.status, data = upa)
p.os[3] <- pchisq(p.mod$chisq,1,lower.tail = F)
p.adj <- p.adjust(p.os,method = 'fdr')

jpeg("UpaOs.jpg")
km <- kmplot(survfit(Surv(OS, OS.Event) ~ upa.status, data = upa), mark='', simple=FALSE,
       xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
       lty.surv=1, lwd.surv=2, col.surv=c(1,2), font.lab=2, font.axis=2, # survival.curves
       col.ci=0, # confidence intervals not plotted
       group.names=c('uPA -','uPA +'),
       group.order=c(1,2), # order of appearance in the n.risk.at table and legend.
       extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
       grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
       legend=FALSE, loc.legend='bottomleft',
       cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("uPA OS  adj. p = ",round(p.adj[1],3)))
mtext("Overall Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

jpeg("UparOs.jpg")
km <- kmplot(survfit(Surv(OS, OS.Event) ~ upar.status, data = upa), mark='', simple=FALSE,
       xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
       lty.surv=1, lwd.surv=2, col.surv=c(1,2), font.lab=2, font.axis=2, # survival.curves
       col.ci=0, # confidence intervals not plotted
       group.names=c('uPA -','uPA +'),
       group.order=c(1,2), # order of appearance in the n.risk.at table and legend.
       extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
       grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
       legend=FALSE, loc.legend='bottomleft',
       cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("uPAR OS  adj. p = ",round(p.adj[2],3)))
mtext("Overall Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

jpeg("Pai1Os.jpg")
km <- kmplot(survfit(Surv(OS, OS.Event) ~ pai1.status, data = upa), mark='', simple=FALSE,
       xaxis.at=seq(0,144,12), xaxis.lab=seq(0,144,12), # n.risk.at
       lty.surv=1, lwd.surv=2, col.surv=c(1,2), font.lab=2, font.axis=2, # survival.curves
       col.ci=0, # confidence intervals not plotted
       group.names=c('uPA -','uPA +'),
       group.order=c(1,2), # order of appearance in the n.risk.at table and legend.
       extra.left.margin=1.75, label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=.9, # labels
       grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9), 
       legend=FALSE, loc.legend='bottomleft',
       cex.lab=1.5, font.lab=2, xaxs='r', bty='L', las=1, tcl=-.2  # other parameters passed to plot()
) 
text(0,0.05,adj=0,cex=1.5,font=2,paste0("PAI-1 OS  adj. p = ",round(p.adj[3],3)))
mtext("Overall Survival (mos)",side=1,line=1.9, font=2, cex=1.5)
mtext("Survival Probability",side=2,line=2.5, font=2, cex=1.5, las=0)
print(km)
dev.off()

# Are the same patients low in upa, upar, and pai1?
upaERneg.id <- which(upa$upa.status==0&upa$ER.Status==0)
uparERneg.id <- which(upa$upar.status==0&upa$ER.Status==0)
pai1ERneg.id <- which(upa$pai1.status==0&upa$ER.Status==0)
sum(!is.na(match(upaERneg.id,uparERneg.id)))/length(upaERneg.id)
# 75% of patients in low upa & ER are in low upar & ER
sum(!is.na(match(upaERneg.id,pai1ERneg.id)))/length(upaERneg.id)
# 50% of patients in low upa & ER are in low pai1 & ER
sum(!is.na(match(uparERneg.id,pai1ERneg.id)))/length(uparERneg.id)
# 62.5% of patients in low upar & ER are in low pai1 & ER
length(intersect(intersect(upaERneg.id,uparERneg.id),pai1ERneg.id))/length(uparERneg.id)
# 50% of patients in low upa & ER are in low upar & ER and low pai1 & ER 
View(upa[intersect(intersect(upaERneg.id,uparERneg.id),pai1ERneg.id),])
          
# add a, b, c, etc
