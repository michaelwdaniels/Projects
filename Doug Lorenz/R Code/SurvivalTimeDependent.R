# Survival Analysis with Time Dependent Covariates
# Author: Michael Daniels
# Date: June 17, 2019

# Purpose: To produce a general outline for an analysis of survival outcomes with time dependent covariates

library(survival)
set.seed(1953)  # a good year
nvisit <- floor(pmin(lung$time/30.5, 12))
response <- rbinom(nrow(lung), nvisit, .05) > 0
badfit <- survfit(Surv(time/365.25, status) ~ response, data=lung)
plot(badfit, mark.time=FALSE, lty=1:2,
     xlab="Years post diagnosis", ylab="Survival")
legend(1.5, .85, c("Responders", "Non-responders"), lty=2:1, bty='n')
fit <- coxph(Surv(time1, time2, status) ~ age + creatinine,
             data=mydata)
head(cgd0)
newcgd <- tmerge(data1=cgd0[, 1:13], data2=cgd0, id=id, tstop=futime)
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime1))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime2))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime3))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime4))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime5))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime6))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime7))
newcgd <- tmerge(newcgd, newcgd, id, enum=cumtdc(tstart))
dim(newcgd)
attr(newcgd, "tcount")
coxph(Surv(tstart, tstop, infect) ~ treat + inherit + steroids + cluster(id), newcgd)

temp <- subset(pbc, id <= 312, select=c(id:sex, stage)) # baseline
pbc2 <- tmerge(temp, temp, id=id, death = event(time, status)) #set range
pbc2 <- tmerge(pbc2, pbcseq, id=id, ascites = tdc(day, ascites),
                 bili = tdc(day, bili), albumin = tdc(day, albumin),
                 protime = tdc(day, protime), alk.phos = tdc(day, alk.phos))
fit1 <- coxph(Surv(time, status==2) ~ log(bili) + log(protime), pbc)
fit2 <- coxph(Surv(tstart, tstop, death==2) ~ log(bili) + log(protime), pbc2) 
rbind('baseline fit' = coef(fit1),
                                                                                        

setwd("~/Desktop/Projects/Doug Lorenz/SurvivalProject")
HallData <- read.csv("~/Desktop/Projects/Doug Lorenz/SurvivalProject/HallData.csv")



