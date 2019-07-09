setwd("~/Desktop/Projects/Doug Lorenz/Ian Mutchnick")
dat <- read.csv("Mutchnick.csv", header=T)
names(dat)[8:9] <- c('t0','t1') # rename preop temp to t0 and initial temp to t1
View(dat)
# Correlation between prep time and last preop temp
cor(dat$preptime,dat$t0,use='pairwise.complete.obs',method='pearson') # 0.041
which(is.na(dat$preptime)) # rows 67 78
which(is.na(dat$t0)) # rows 19 21
# 116 complete observations have a correlation of 0.041 between prep time and last preop temp

# Correlation between prep time and mininum intraoperative temp
temp_min <- apply(dat[,8:17],1,function(x) min(x,na.rm=T)) # only 3 instances with complete temp data
cor(dat$preptime,temp_min,use='pairwise.complete.obs',method='pearson') # 0.017
which(is.na(dat$preptime)) # rows 67 78
which(is.na(temp_min)) # none
# 118 complete observations have a correlation of 0.017 between prep time and minimum temp (includes last preop temp)

# Correlation between the controls for prep time and last preop temp
cor(dat$preptime[which(dat$group=='control')],dat$t0[which(dat$group=='control')],use='pairwise.complete.obs',method='pearson') # 0.070
which(is.na(dat$preptime[which(dat$group=='control')])) # rows 29 40
which(is.na(dat$t0[which(dat$group=='control')])) # none
# 118 complete observations have a correlation of 0.071 between prep time and last preop temp for controls

# Correlation between the interventions for prep time and last preop temp
cor(dat$preptime[which(dat$group=='intervention')],dat$t0[which(dat$group=='intervention')],use='pairwise.complete.obs',method='pearson') # -0.073
which(is.na(dat$preptime[which(dat$group=='intervention')])) # rows 29 40
which(is.na(dat$t0[which(dat$group=='intervention')])) # rows 19 21
# 118 complete observations have a correlation of -0.073 between prep time and last preop temp for interventions

# T-test between groups for preptime
t.test(dat$preptime[which(dat$group=='control')],dat$preptime[which(dat$group=='intervention')]) 
# t = -0.30115, df = 83.914, p-value = 0.764
# 0.9418750 0.9627193 

# T-test between groups for last preop temp
t.test(dat$t0[which(dat$group=='control')],dat$t0[which(dat$group=='intervention')])
# t = -6.3496, df = 95.798, p-value = 7.2e-09
# 36.25122  36.82222 

# Convert wide format to long
dat.long <- reshape(dat, varying=paste("t",c(0,1,15,30,45,60,75,90,105,120),sep=""), direction="long",
   v.names=c("t"), timevar="time", idvar="id", times=c(0,1,15,30,45,60,75,90,105,120))

# load additional libraries for mixed models and spaghetti plots
library(nlme) # for linear mixed effects models
library(lattice) # for spaghetti plots (xyplot)
library(lme4) # for generalized mixed effects model (logistic regression)

# First let's run the models without imputation for missing data while adjusting for prep time and age and comparing model fits using AIC
mod1.prep <- lme(t~factor(time)*group+preptime, random=~1|id, data=dat.long, na.action=na.omit) # random intercept by id  
mod1.age <- lme(t~factor(time)*group+age, random=~1|id, data=dat.long, na.action=na.omit) # random intercept by id  
mod1.age.prep <- lme(t~factor(time)*group+preptime+age, random=~1|id, data=dat.long, na.action=na.omit) # random intercept by id  
mod1.age.prep.quad <- lme(t~time*group+I(time^2)*group, random=~1|id, data=dat.long, na.action=na.omit) # random slope for time by group
mod1.age.prep.quad.slope <- lme(t~time*group+I(time^2)*group, random=~time|id, data=dat.long, na.action=na.omit) # random slope for time by group
mod1.age.prep.quad.slopesq <- lme(t~time*group+I(time^2)*group, random=~time+I(time^2)|id, data=dat.long, na.action=na.omit) # random slope for time by group
mod1.age.prep.quad.group <- update(mod1.age.prep.quad.slope, random=~time*group|id)
#mod1.age.prep.quad.groupsq <- update(mod1.age.prep.quad.slopesq, random=~time*group+I(time^2)*group|id)

summary(mod1.prep)$AIC                  # 1829.695
summary(mod1.age)$AIC                   # 1840.958
summary(mod1.age.prep)$AIC              # 1826.374 *
summary(mod1.age.prep.quad)$AIC         # 1887.837 
summary(mod1.age.prep.quad.slope)$AIC   # 1851.098 
summary(mod1.age.prep.quad.slopesq)$AIC # 1853.645 
summary(mod1.age.prep.quad.group)$AIC   # 1846.565
#summary(mod1.age.prep.quad.groupsq)$AIC # did not converge

modf <- lme(t~factor(time):group-1, random=~1|id, data=dat.long, na.action=na.omit, subset=t>23)

summary(modf)

f <- fixef(modf)
m <- matrix(f, nrow=10)
x <- c(-15,0,15,30,45,60,75,90,105,120)

plot(x, m[,1], type="o", ylim=c(34,39), las=1, bty="l", ylab='Mean Temperature (\u00B0C)', xlab='Time (min)', pch=19, col="blue")
   lines(x, m[,2], type="o", pch=19, col="red")
   abline(h=36, lty=2)
   legend('bottomright',c('Intervention','Control'), col=c(2,4), lty=1)
   
xyplot(t~time|group, groups=id, data=dat.long, type="l", subset=t>23,
   panel = function(x,y,groups,...) { 
      panel.xyplot(x,y,groups,...);
      panel.lines(x=c(0,120), y=rep(36,2), type="l", lwd=2, col="black")
})

dat.long$resp <- 1*(dat.long$t<36)
mod <- glmer(resp~time*group + (1 | id), data=dat.long, family=binomial)

tmp <- tapply(dat.long$resp, list(dat.long$time, dat.long$group), mean, na.rm=T)
apply(tmp, 1, function(x) (x[1]*(1-x[2]))/(x[2]*(1-x[1])))

# Second, lets impute by carrying forward last value 
# There are two instances with missing data for last preop temp (t0) - impute using mean last preop temp for designated group identification
# This is important due to group difference
dat1 <- dat
t0.mean.int <- mean(dat1$t0[which(dat1$group=='intervention')],na.rm=T)
dat1$t0[which(is.na(dat1$t0))] <- t0.mean.int

# Create a loop that carries forward prior value
for(i in 9:17) {
   dat1[,i] <- ifelse(is.na(dat1[,i]),dat1[,i-1],dat1[,i])
}

# Convert wide format to long
dat.long1 <- reshape(dat1, varying=paste("t",c(0,1,15,30,45,60,75,90,105,120),sep=""), direction="long",
                    v.names=c("t"), timevar="time", idvar="id", times=c(0,1,15,30,45,60,75,90,105,120))
# Repeat above analysis with complete temp data
mod1.prep <- lme(t~factor(time)*group+preptime, random=~1|id, data=dat.long1, na.action=na.omit) # random intercept by id  
mod1.age <- lme(t~factor(time)*group+age, random=~1|id, data=dat.long1, na.action=na.omit) # random intercept by id  
mod1.age.prep <- lme(t~factor(time)*group+preptime+age, random=~1|id, data=dat.long1, na.action=na.omit) # random intercept by id  
mod1.age.prep.quad <- lme(t~time*group+I(time^2)*group, random=~1|id, data=dat.long1, na.action=na.omit) # random slope for time by group
mod1.age.prep.quad.slope <- lme(t~time*group+I(time^2)*group, random=~time|id, data=dat.long1, na.action=na.omit) # random slope for time by group
mod1.age.prep.quad.slopesq <- lme(t~time*group+I(time^2)*group, random=~time+I(time^2)|id, data=dat.long1, na.action=na.omit) # random slope for time by group
mod1.age.prep.quad.group <- update(mod1.age.prep.quad.slope, random=~time*group|id)
#mod1.age.prep.quad.groupsq <- update(mod1.age.prep.quad.slopesq, random=~time*group+I(time^2)*group|id)

summary(mod1.prep)$AIC                  # 2178.055
summary(mod1.age)$AIC                   # 2198.693
summary(mod1.age.prep)$AIC              # 2175.984
summary(mod1.age.prep.quad)$AIC         # 2211.327 
summary(mod1.age.prep.quad.slope)$AIC   # 2112.307 
summary(mod1.age.prep.quad.slopesq)$AIC # 2106.853 *
summary(mod1.age.prep.quad.group)$AIC   # 2112.75 
#summary(mod1.age.prep.quad.groupsq)$AIC # did not converge

modf <- lme(t~factor(time):group-1, random=~1|id, data=dat.long1, na.action=na.omit, subset=t>23)

summary(modf)
   
f <- fixef(modf)
m <- matrix(f, nrow=10)
x <- c(-15,0,15,30,45,60,75,90,105,120)

plot(x, m[,1], type="o", ylim=c(34,39), las=1, bty="l", ylab='Mean Temperature (\u00B0C)', xlab='Time (min)', pch=19, col="blue")
   lines(x, m[,2], type="o", pch=19, col="red")
   abline(h=36, lty=2)
   legend('bottomright',c('Intervention','Control'), col=c(2,4), lty=1)
xyplot(t~time|group, groups=id, data=dat.long1, type="l", subset=t>23,
       panel = function(x,y,groups,...) { 
          panel.xyplot(x,y,groups,...);
          panel.lines(x=c(0,120), y=rep(36,2), type="l", lwd=2, col="black")
       })

dat.long1$resp <- 1*(dat.long1$t<36)

mod <- glmer(resp~time*group + (1 | id), data=dat.long1, family=binomial)
summary(mod)

# Create a plot that compares the proportion of hypothermic patients at each time point
# to the predicted values from the best model
# First find the proportion of patients at each time point

hypo.actual <- tapply(dat.long1$resp, list(dat.long1$time, dat.long1$group), mean, na.rm=T)
apply(tmp, 1, function(x) (x[1]*(1-x[2]))/(x[2]*(1-x[1])))

hypo.pred <- tapply(ifelse(predict(mod1.age.prep.quad.slopesq)<36,1,0), list(dat.long1$time, dat.long1$group), mean, na.rm=T)

plot(x, hypo.actual[,1], type="l", ylim=c(0,1), las=1, bty="l", ylab='Proportion Hypothermic (%)', xlab='Time (min)', pch=19, col="blue",
     main='Comparing Actual to Predicted Proportion\nHypothermic Between Intervention and Controls')
lines(x, hypo.actual[,2], type="l", pch=19, col="red")
lines(x, hypo.pred[,1], type="l", pch=19, col="blue",lty=2)
lines(x, hypo.pred[,2], type="l", pch=19, col="red",lty=2)
legend('topright',c('Intervention','Control','Predicted'), col=c(2,4,1), lty=c(1,1,2))

