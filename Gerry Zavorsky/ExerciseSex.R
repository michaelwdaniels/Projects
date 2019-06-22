# 'Fixed Effect and Random Effects Meta-Analysis'

library(meta)

setwd("~/Desktop/Projects/Gerry Zavorsky")

# 1. Read in the data
data1 <- read.csv("dataset01.csv", as.is=TRUE)
dataGerry <- read.csv("gerry.csv", as.is=TRUE)
dataGerry2 <- read.csv("gerry2.csv", as.is=TRUE)
dataGerry3 <- read.csv("gerry3.csv", as.is=TRUE)

# 2. Calculate mean difference and its standard error for

MD <- with(data1[1,], Me - Mc)
MD1 <- with(dataGerry, Me - Mc)
MD2 <- with(dataGerry2, Me - Mc)
MD3 <- with(dataGerry3, Me - Mc)

seMD <- with(data1[1,], sqrt(Se^2/Ne + Sc^2/Nc))
seMD1 <- with(dataGerry, sqrt(Se^2/Ne + Sc^2/Ne))
seMD2 <- with(dataGerry2, sqrt(Se^2/Ne + Sc^2/Ne))
seMD3 <- with(dataGerry3, sqrt(Se^2/Ne + Sc^2/Ne))

# 3. Print mean difference and limits of 95% confidence
#    interval using round function to show only two digits:

round(c(MD, MD + c(-1,1) * qnorm(1-(0.05/2)) * seMD), 2)
CI1 <- matrix(0,length(MD1),3)
for(i in 1:length(MD1)){CI1[i,] <- round(c(MD1[i], MD1[i] + c(-1,1) * qnorm(1-(0.05/2)) * seMD1[i]), 2)}
CI2 <- matrix(0,length(MD2),3)
for(i in 1:length(MD2)){CI2[i,] <- round(c(MD2[i], MD2[i] + c(-1,1) * qnorm(1-(0.05/2)) * seMD2[i]), 2)}
CI3 <- matrix(0,length(MD3),3)
for(i in 1:length(MD3)){CI3[i,] <- round(c(MD3[i], MD3[i] + c(-1,1) * qnorm(1-(0.05/2)) * seMD3[i]), 2)}





with(data1[1, ],
     print(metacont(Ne, Me, Se, Nc, Mc, Sc),
           digits=2))
with(dataGerry[, ],
     print(metacont(Ne, Me, Se, Ne, Mc, Sc),
           digits=2))

print(metacont(Ne, Me, Se, Nc, Mc, Sc,
               data=data1, subset=1), digits=2)

zscore <- MD/seMD
round(c(zscore, 2*pnorm(abs(zscore), lower.tail=FALSE)), 4)
zscore1 <- MD1/seMD1
round(c(zscore1, 2*pnorm(abs(zscore1), lower.tail=FALSE)), 4)


# 1. Read in the data:
data2 <- read.csv("dataset02.csv")

# 2. As usual, to view an object, type its name:
data2

# 3. Calculate total sample sizes
summary(data2$Ne+data2$Nc)



# 1. Calculate standardised mean difference (SMD) and
#    its standard error (seSMD) for study 1 (Blashki) of
#    dataset data2:
N <- with(dataGerry[,], Ne + Ne)
SMD <- with(dataGerry[,],
            (1 - 3/(4 * N - 9)) * (Me - Mc) /
            sqrt(((Ne - 1) * Se^2 + (Ne - 1) * Sc^2)/(N - 2)))
seSMD <- with(dataGerry[,],
              sqrt(N/(Ne * Ne) + SMD^2/(2 * (N - 3.94))))

# 2. Print standardised mean difference and limits of 95% CI
#    interval using round function to show only two digits:
round(c(SMD, SMD + c(-1,1) * qnorm(1-(0.05/2)) * seSMD), 2)
SMDCI1 <- matrix(0,length(MD1),3)
for(i in 1:length(SMD)){SMDCI1[i,] <- round(c(SMD[i], SMD[i] + c(-1,1) * qnorm(1-(0.05/2)) * seSMD[i]), 2)}

print(metacont(Ne, Me, Se, Ne, Mc, Sc, sm="SMD",
               data=dataGerry), digits=2)

# 1. Calculate mean difference, variance and weights

MD <- with(dataGerry, Me - Mc)
varMD <- with(dataGerry, Se^2/Ne + Sc^2/Ne)
weight <- 1/varMD

# 2. Calculate the inverse variance estimator

round(weighted.mean(MD, weight), 4)

# 3. Calculate the variance

round(1/sum(weight), 4)

data.comb <- rbind(dataGerry,dataGerry2,dataGerry3)
data.comb$Type <- c(rep('Aerobic',4),rep('Strength',4),rep('Endurance',2)) 
data.comb.alt <- rbind(dataGerry,dataGerry2[4,])
data.comb$ID <- c(1,2,3,4,1,2,2,5,1,5) 

mc <- metacont(Ne, Me, Se, Ne, Mc, Sc,
               data=data.comb,comb.random = TRUE,comb.fixed = FALSE,
               label.e = 'Experimental', byvar=data.comb$Type,
               title='Effects on Aerobic, Endurance, and Strength Performance by prior Intercourse',
               studlab=paste(Author))
mc
pdf(file="Gerry-Fig4.pdf", width=13) # uncomment line to generate PDF file
forest(mc, xlab= "Favors Abstinence                Favors sex",digits=1,
       digits.se=1,digits.zval=2,digits.pval = 2,
       digits.pval.Q = 2,digits.Q = 2,digits.tau2 = 2)
invisible(dev.off()) # uncomment line to save PDF file

mc1 <- metacont(Ne, Me, Se, Ne, Mc, Sc,
                data=dataGerry,
                studlab=paste(Author),title='Prior Intercourse Effect on Exercise Performance')
mc1

mc2 <- metacont(Ne, Me, Se, Ne, Mc, Sc,
                data=dataGerry2,
                studlab=paste(Author),title='Prior Intercourse Effect on Exercise Performance')
mc2

mc3 <- metacont(Ne, Me, Se, Ne, Mc, Sc,
                data=dataGerry3,
                studlab=paste(Author),title='Prior Intercourse Effect on Exercise Performance')
mc3

mc4 <- metacont(Ne, Me, Se, Ne, Mc, Sc,
                data=data.comb.alt,
                studlab=paste(Author),title='Prior Intercourse Effect on Exercise Performance')
mc4

round(c(mc1$TE.fixed, mc1$seTE.fixed^2), 4)

print(mc1, digits=2)

mc1$w.fixed[1]
sum(mc1$w.fixed)
round(100*mc1$w.fixed[1] / sum(mc1$w.fixed), 2)

pdf(file="Gerry-Fig1.pdf", width=13) # uncomment line to generate PDF file
forest(mc1, xlab= "Favors Abstinence               Favors sex",digits=1,
       digits.se=1,digits.zval=2,digits.pval = 2,
       digits.pval.Q = 2,digits.Q = 2,digits.tau2 = 2)
invisible(dev.off()) # uncomment line to save PDF file

pdf(file="Gerry-Fig2.pdf", width=13) # uncomment line to generate PDF file
forest(mc2, xlab= "Favors Abstinence               Favors sex",digits=1,
       digits.se=1,digits.zval=2,digits.pval = 2,
       digits.pval.Q = 2,digits.Q = 2,digits.tau2 = 2)
invisible(dev.off()) # uncomment line to save PDF file

pdf(file="Gerry-Fig3.pdf", width=13) # uncomment line to generate PDF file
forest(mc3, xlab= "Favors Abstinence               Favors sex",digits=1,
       digits.se=1,digits.zval=2,digits.pval = 2,
       digits.pval.Q = 2,digits.Q = 2,digits.tau2 = 2)
invisible(dev.off()) # uncomment line to save PDF file

pdf(file="Gerry-Fig4.pdf", width=13) # uncomment line to generate PDF file
forest(mc, xlab= "Favors Abstinence                Favors sex",digits=1,
       digits.se=1,digits.zval=2,digits.pval = 2,
       digits.pval.Q = 2,digits.Q = 2,digits.tau2 = 2)
invisible(dev.off()) # uncomment line to save PDF file

pdf(file="Gerry-Fig5.pdf", width=13) # uncomment line to generate PDF file
forest(mc4, xlab= "Favors Abstinance               Favors sex",digits=1,
       digits.se=1,digits.zval=2,digits.pval = 2,
       digits.pval.Q = 2,digits.Q = 2,digits.tau2 = 2)
invisible(dev.off()) # uncomment line to save PDF file


# 1. Apply generic inverse variance method
mc1.gen <- metagen(mc1$TE, mc1$seTE, sm="MD")
# 2. Same result
mc1.gen <- metagen(TE, seTE, data=mc1, sm="MD")
# 3. Print results for fixed effect and random effects method
c(mc1$TE.fixed, mc1$TE.random)
c(mc1.gen$TE.fixed, mc1.gen$TE.random)





# 1. Calculate standardised mean difference,
#    variance and weights
N <- with(dataGerry, Ne + Ne)
SMD <- with(dataGerry,
            (1 - 3/(4 * N - 9)) * (Me - Mc)/
            sqrt(((Ne - 1) * Se^2 + (Ne - 1) * Sc^2)/(N - 2)))
varSMD <- with(dataGerry,
               N/(Ne * Ne) + SMD^2/(2 * (N - 3.94)))
weight <- 1/varSMD
# 2. Calculate the inverse variance estimator
round(weighted.mean(SMD, weight), 4)
# 3. Calculate the variance
round(1/sum(weight), 4)





mc2 <- metacont(Ne, Me, Se, Ne, Mc, Sc, sm="SMD",
                data=dataGerry)
round(c(mc2$TE.fixed, mc2$seTE.fixed^2), 4)





print(summary(mc2), digits=2)





#
# Web-appendix only
#
# 1. Conduct meta-analyses
# 1a. DerSimonian-Laird estimator (default)
mg1.DL <- metagen(TE, seTE, data=mc1)
# 1b. Paule-Mandel estimator
mg1.PM <- metagen(TE, seTE, data=mc1, method.tau="PM")
# 1c. Restricted maximum-likelihood estimator
mg1.RM <- metagen(TE, seTE, data=mc1, method.tau="REML")
# 1d. Maximum-likelihood estimator
mg1.ML <- metagen(TE, seTE, data=mc1, method.tau="ML")
# 1e. Hunter-Schmidt estimator
mg1.HS <- metagen(TE, seTE, data=mc1, method.tau="HS")
# 1f. Sidik-Jonkman estimator
mg1.SJ <- metagen(TE, seTE, data=mc1, method.tau="SJ")
# 1g. Hedges estimator
mg1.HE <- metagen(TE, seTE, data=mc1, method.tau="HE")
# 1h. Empirical Bayes estimator
mg1.EB <- metagen(TE, seTE, data=mc1, method.tau="EB")
# 2. Extract between-study variance tau-squared
tau2 <- data.frame(tau2=round(c(0,
                     mg1.DL$tau^2, mg1.PM$tau^2,
                     mg1.RM$tau^2, mg1.ML$tau^2,
                     mg1.HS$tau^2, mg1.SJ$tau^2,
                     mg1.HE$tau^2, mg1.EB$tau^2), 2),
                   row.names=c("FE", "DL", "PM", "REML", "ML",
                     "HS", "SJ", "HE", "EB"))
# 3. Print tau-squared values
t(tau2)
# 4. Create dataset with summaries
res <- data.frame(MD=c(mg1.DL$TE.fixed,
                    mg1.DL$TE.random, mg1.PM$TE.random,
                    mg1.RM$TE.random, mg1.ML$TE.random,
                    mg1.HS$TE.random, mg1.SJ$TE.random,
                    mg1.HE$TE.random, mg1.EB$TE.random),
                  seMD=c(mg1.DL$seTE.fixed,
                    mg1.DL$seTE.random, mg1.PM$seTE.random,
                    mg1.RM$seTE.random, mg1.ML$seTE.random,
                    mg1.HS$seTE.random, mg1.SJ$seTE.random,
                    mg1.HE$seTE.random, mg1.EB$seTE.random),
                  method=c("",
                    "DerSimonian-Laird",
                    "Paule-Mandel",
                    "Restricted maximum-likelihood",
                    "Maximum-likelihood",
                    "Hunter-Schmidt",
                    "Sidik-Jonkman",
                    "Hedges",
                    "Empirical Bayes"),
                  tau2=tau2,
                  model=c("Fixed-effect model",
                    rep("Random-effect model", 8)))
# 5. Do meta-analysis
m <- metagen(MD, seMD, data=res,
             studlab=method,
             sm="MD",
             comb.fixed=FALSE, comb.random=FALSE,
             byvar=model)
# 6. Do forest plot
# pdf(file="Schwarzer-Fig2.5.pdf", width=8.1, height=4.0) # uncomment line to generate PDF file
forest(m,
       hetstat=FALSE, smlab="",
       leftcols=c("studlab", "tau2"),
       leftlabs=c("Method", "Between-study\nheterogeneity"),
       print.byvar=FALSE)
# invisible(dev.off()) # uncomment line to save PDF file





mc2.hk <- metacont(Ne, Me, Se, Ne, Mc, Sc, sm="SMD",
                   data=dataGerry, comb.fixed=FALSE,
                   hakn=TRUE)





mc2.hk <- metagen(TE, seTE, data=mc2, comb.fixed=FALSE,
                  hakn=TRUE)





print(summary(mc2.hk), digits=2)





print(summary(mc1, prediction=TRUE), digits=2)





# pdf(file="Schwarzer-Fig2.6.pdf", width=11.0) # uncomment line to generate PDF file
forest(mc1, prediction=TRUE, col.predict="black")
# invisible(dev.off()) # uncomment line to save PDF file





# 1. Read in the data:
data3 <- read.csv("dataset03.csv")
# 2. As usual, to view an object, type its name:
data3





mc3 <- metacont(Ne, Me, Se, Nc, Mc, Sc, data=data3,
                studlab=paste(author, year))





mc3$studlab[mc3$w.fixed==0]





print(summary(mc3), digits=2)





mc3s <- metacont(Ne, Me, Se, Nc, Mc, Sc, data=data3,
                 studlab=paste(author, year),
                 byvar=duration, print.byvar=FALSE)





mc3s <- update(mc3, byvar=duration, print.byvar=FALSE)





#mc3nodata <- metacont(Ne, Me, Se, Nc, Mc, Sc, data=data3,
#                     studlab=paste(author, year),
#                     keepdata=FALSE)
#update(mc3nodata)





print(summary(mc3s), digits=2)





# pdf(file="Schwarzer-Fig2.8.pdf", width=10.35, height=9.25) # uncomment line to generate PDF file
forest(mc3s, xlim=c(-0.5, 0.2),
       xlab="Difference in mean number of acute exacerbations per month")
# invisible(dev.off()) # uncomment line to save PDF file





print(metacont(Ne, Me, Se, Nc, Mc, Sc, data=data3,
               subset=duration=="<= 3 months",
               studlab=paste(author, year)),
      digits=2)





print(update(mc3, subset=duration=="<= 3 months"),
      digits=2)





# 1. Read in the data
data4 <- read.csv("dataset04.csv")
# 2. Print data
data4





mg1 <- metagen(logHR, selogHR,
               studlab=paste(author, year), data=data4,
               sm="HR")





print(mg1, digits=2)





# 1. Read in the data
data5 <- read.csv("dataset05.csv")
# 2. Print data
data5





mg2 <- metagen(mean, SE, studlab=paste(author, year),
               data=data5, sm="MD")





print(summary(mg2), digits=2)





# 1. Read in the data
data6 <- read.csv("dataset06.csv")
# 2. Print data
data6





mg3 <- metagen(b, SE, studlab=paste(author, year),
               data=data6, sm="RR", backtransf=FALSE)





summary(mg3)
