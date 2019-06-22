library(ade4)
setwd("~/Desktop/Projects/Will Stern")

Outcomes.ER.PR.EIA.MTA.Kd <- read.csv("~/Desktop/Projects/Will Stern/Outcomes ER PR EIA MTA Kd.csv")
data <- Outcomes.ER.PR.EIA.MTA.Kd[complete.cases(Outcomes.ER.PR.EIA.MTA.Kd),]
dim(data)
labels.columns <- colnames(data)

par(mfrow=c(1,1))

labels.columns
varidx <- c(2,5,7)
ta <- as.data.frame(data[,varidx])
colnames(ta) <- labels.columns[varidx]
png('DFS ER PR.png')
triangle.plot(ta, label = as.character(1:nrow(ta)), clabel = 0, 
              cpoint = 1, draw.line = TRUE, addaxes = TRUE, addmean = TRUE, 
              labeltriangle = FALSE, sub = "DFS ER PR ", csub = 1, possub = "topright", 
              show.position = FALSE, scale = TRUE, min3 = NULL, max3 = NULL, box = TRUE)
mtext('DFS',side=2,line=-2,las=1)
mtext('ER MTA',side=1,line=2.5)
mtext('PR MTA',side=4,line=-4,las=1,padj=-1)
dev.off()

varidx <- c(3,5,7)
ta <- as.data.frame(data[,varidx])
colnames(ta) <- labels.columns[varidx]
png('OS ER PR.png')
triangle.plot(ta, label = as.character(1:nrow(ta)), clabel = 0, 
              cpoint = 1, draw.line = TRUE, addaxes = TRUE, addmean = TRUE, 
              labeltriangle = FALSE, sub = "OS ER PR ", csub = 1, possub = "topright", 
              show.position = FALSE, scale = TRUE, min3 = NULL, max3 = NULL, box = TRUE)
mtext('OS',side=2,line=-2,las=1,padj=2)
mtext('ER MTA',side=1,line=2.5)
mtext('PR MTA',side=4,line=-4,las=1,padj=-1)
dev.off()

varidx <- c(2,6,8)
ta <- scale(as.data.frame(data[,varidx]), center = FALSE)
apply(ta,2,summary)
colnames(ta) <- labels.columns[varidx]
png('DFS ER PR Kd.png')
triangle.plot(ta, label = as.character(1:nrow(ta)), clabel = 0, 
              cpoint = 1, draw.line = TRUE, addaxes = TRUE, addmean = TRUE, 
              labeltriangle = FALSE, sub = "DFS Kd ER PR ", csub = 1, possub = "topright", 
              show.position = FALSE, scale = TRUE, min3 = NULL, max3 = NULL, box = TRUE)
mtext('DFS',side=2,line=-2,las=1)
mtext('ER MTA Kd',side=1,line=3.5,adj=0.4)
mtext('PR MTA Kd',side=4,line=-4,las=1,padj=-1)
dev.off()

varidx <- c(3,6,8)
ta <- scale(as.data.frame(data[,varidx]), center = FALSE)
colnames(ta) <- labels.columns[varidx]
png('OS ER PR Kd.png')
triangle.plot(ta, label = as.character(1:nrow(ta)), clabel = 0, 
              cpoint = 1, draw.line = TRUE, addaxes = TRUE, addmean = TRUE, 
              labeltriangle = FALSE, sub = "OS Kd ER PR ", csub = 1, possub = "topright", 
              show.position = FALSE, scale = TRUE, min3 = NULL, max3 = NULL, box = TRUE)
mtext('OS',side=2,line=-1.5,las=1,padj=3)
mtext('ER MTA Kd',side=1,line=3.5,adj=0.3)
mtext('PR MTA Kd',side=4,line=-4,las=1,padj=1)
dev.off()


################
# Biplot example with scaling and vignette examples


par(mfrow = c(1, 1))

ta1 <- data[which(data$ER.EIA<median(data$ER.EIA)),]
ta2 <- data[-which(data$ER.EIA<median(data$ER.EIA)),]

triangle.biplot (ta1[,c(2,6,8)], ta2[,c(2,6,8)], label = as.character(1:nrow(ta1)), 
                 draw.line = TRUE, show.position = TRUE, scale = TRUE)

# from vignette
data(euro123)
tot <- rbind.data.frame(euro123$in78, euro123$in86, euro123$in97)
row.names(tot) <- paste(row.names(euro123$in78), rep(c(1, 2, 3), rep(12, 3)), sep = "")
triangle.plot(tot, label = row.names(tot), clab = 1)
par(mfrow = c(2, 2))
triangle.plot(euro123$in78, clab = 0, cpoi = 2, addmean = TRUE, show = FALSE)
triangle.plot(euro123$in86, label = row.names(euro123$in78), clab = 0.8)
triangle.biplot(euro123$in78, euro123$in86)
triangle.plot(rbind.data.frame(euro123$in78, euro123$in86), clab = 1,
              addaxes = TRUE, sub = "Principal axis", csub = 2, possub = "topright")
triangle.plot(euro123[[1]], min3 = c(0, 0.2, 0.3), max3 = c(0.5, 0.7, 0.8),
              clab = 1, label = row.names(euro123[[1]]), addax = TRUE)
triangle.plot(euro123[[2]], min3 = c(0, 0.2, 0.3), max3 = c(0.5, 0.7, 0.8),
              clab = 1, label = row.names(euro123[[1]]), addax = TRUE)
triangle.plot(euro123[[3]], min3 = c(0, 0.2, 0.3), max3 = c(0.5, 0.7, 0.8),
              clab = 1, label = row.names(euro123[[1]]), addax = TRUE)
triangle.plot(rbind.data.frame(euro123[[1]], euro123[[2]], euro123[[3]]))
par(mfrow = c(1, 1))
wtriangleplot <- cbind.data.frame(a = runif(100), b = runif(100), c = runif(100, 4, 5))
wtriangleplot <- triangle.plot(wtriangleplot)
points(wtriangleplot, col = "blue", cex = 2)
wtriangleplot <- colMeans(wtriangleplot)
points(wtriangleplot[1], wtriangleplot[2], pch = 20, cex = 3, col = "red")
rm(wtriangleplot)
