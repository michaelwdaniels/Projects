setwd("~/Desktop/Projects/Seth Sereff")
uPA <- read.csv("uPA.csv")
uPA$Nodal_Status[grep('neg',uPA$Nodes)] <- 0
uPA$Nodal_Status[c(grep('pos,1+/',uPA$Nodes,fixed=T),grep('pos, 1+/',uPA$Nodes,fixed=T),
            grep('pos,2+/',uPA$Nodes,fixed=T),grep('pos, 2+/',uPA$Nodes,fixed=T),
            grep('pos,3+/',uPA$Nodes,fixed=T),grep('pos, 3+/',uPA$Nodes,fixed=T))] <- 1
idx <- c(grep('pos,1+/',uPA$Nodes,fixed=T),grep('pos, 1+/',uPA$Nodes,fixed=T),
         grep('pos,2+/',uPA$Nodes,fixed=T),grep('pos, 2+/',uPA$Nodes,fixed=T),
         grep('pos,3+/',uPA$Nodes,fixed=T),grep('pos, 3+/',uPA$Nodes,fixed=T),grep('neg',uPA$Nodes))
uPA$Nodal_Status[-idx] <- 2

table(uPA$Nodal_Status)
# 0 - 125
# 1 - 69
# 2 - 53
uPA$Nodes[idx]

library(vioplot)
par(mfrow=c(2,1))
jpeg('uPA.jpeg')
vioplot(uPA$uPA[uPA$Nodal_Status==0],uPA$uPA[uPA$Nodal_Status==1],uPA$uPA[uPA$Nodal_Status==2],
        names=c('Negative','1-3 Positive','>3 Positive'),col='grey')
title('uPA Expression by Nodal Status',ylab = 'Gene Expression',xlab = 'Nodal Status')
abline(lm(uPA$uPA~uPA$Nodal_Status))
text(3,2,paste0('p = ',round(summary(lm(uPA$uPA~uPA$Nodal_Status))$coef[2,4],3)))
mtext(c('n = 125','n = 69','n = 53'),1,at=1:3,line=2)
dev.off()

vioplot(uPA$uPA,col='grey')


jpeg('uPAR.jpeg')
vioplot(uPA$uPAR[uPA$Nodal_Status==0],uPA$uPAR[uPA$Nodal_Status==1],uPA$uPAR[uPA$Nodal_Status==2],
        names=c('Negative','1-3 Positive','>3 Positive'),col='grey')
title('uPAR Expression by Nodal Status',ylab = 'Gene Expression',xlab = 'Nodal Status')
abline(lm(uPA$uPAR~uPA$Nodal_Status))
text(3,3,paste0('p = ',round(summary(lm(uPA$uPAR~uPA$Nodal_Status))$coef[2,4],3)))
mtext(c('n = 125','n = 69','n = 53'),1,at=1:3,line=2)
dev.off()

jpeg('PAI_1.jpeg')
vioplot(uPA$PAI_1[uPA$Nodal_Status==0],uPA$PAI_1[uPA$Nodal_Status==1],uPA$PAI_1[uPA$Nodal_Status==2],
        names=c('Negative','1-3 Positive','>3 Positive'),col='grey')
title('PAI-1 Expression by Nodal Status',ylab = 'Gene Expression',xlab = 'Nodal Status')
abline(lm(uPA$PAI_1~uPA$Nodal_Status))
text(3,3,paste0('p = ',round(summary(lm(uPA$PAI_1~uPA$Nodal_Status))$coef[2,4],3)))
mtext(c('n = 125','n = 69','n = 53'),1,at=1:3,line=2)
dev.off()

jpeg('PAI_2.jpeg')
vioplot(uPA$PAI_2[uPA$Nodal_Status==0],uPA$PAI_2[uPA$Nodal_Status==1],uPA$PAI_2[uPA$Nodal_Status==2],
        names=c('Negative','1-3 Positive','>3 Positive'),col='grey')
title('PAI-2 Expression by Nodal Status',ylab = 'Gene Expression',xlab = 'Nodal Status')
abline(lm(uPA$PAI_2~uPA$Nodal_Status))
text(3,6,paste0('p = ',round(summary(lm(uPA$PAI_2~uPA$Nodal_Status))$coef[2,4],3)))
mtext(c('n = 125','n = 69','n = 53'),1,at=1:3,line=2)
dev.off()
