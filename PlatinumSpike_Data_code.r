#-------------------

# Problem 1

#-------------------
# Legacy code from Dr. Gaile
rm(list=ls())
library (knitr)
library (GEOquery)
# originally, I ran this in the HW1 directory:
# gse19439=getGEO("GSE19439",GSEMatrix=T)
# gse19439=gse19439[[1]]
# save(file="HW1/gse19439.RData",gse19439)
# so that now I can save time and just load:
load( le ="HW1/gse19439.RData")
# make a factor object for Control, latent TB and Positive TB
tmp=as.character(pData(phenoData(gse19439))[,1])
J=length(tmp) # J=number of samples
TBgroup=rep("",J)
for(j in 1:J) TBgroup[j]=substring(tmp[j],1,3)
# make a factor for TBgroup
FTB=factor(TBgroup,levels=c("CON","LTB","PTB"))
# get our expression set
X=exprs(gse19439)
# do a quick kruskal-wallis scan
myKrusk=function(i){
  cat(i,"...",fill=F)
  kruskal.test(x=X[i,],g=FTB)$p.value
}
# originally, I ran this in the HW1 directory:
# myPvals=mapply(myKrusk,1:(dim(X)[1])) ;save(file="HW1/myPvals.RData",myPvals)
# so that now I can save time and just load:
load("HW1/myPvals.RData")
# populate vector with last names of the groups in the class.
# note: the code is written this way, becasue it used to be student nameas and not Group names
GroupLabels=c("Group I","Group II","Group III","Group IV")
# pick the best 4 p-values and assign them to the students.
best4=order(myPvals)[1:4]
# print out list of best 4
print(best4)
7
for(i in 1:length(GroupLabels))
  cat("Group Label:",GroupLabels[i],
      "\t\t row assignment:",best4[i],fill=T)
## Group Label: Group I row assignment: 6874
## Group Label: Group II row assignment: 10685
## Group Label: Group III row assignment: 26058
## Group Label: Group IV row assignment: 47526
#-------------------------------------
# Problem 1 (our own part)
#-------------------------------------
?AnnotatedDataFrame
featureNames(gse19439[47526,])
myfeat <???? featureData(gse19439[47526,])
varLabels(myfeat)
myvarmdat=varMetadata(myfeat)
mypdat=pData(myfeat)
mypdat
summary(mypdat$Probe Type)
summary(mypdat$ILMN Gene)
# Preparing for visualization in Problem 2
featVal.g4<???? exprs(gse19439[47526,])
#-------------------------------------
# Problem 2
#-------------------------------------
# Library calls
wants <???? c("beanplot","gplots","PresenceAbsence","plotROC","pROC","ROCR")
has <???? wants %in% rownames(installed.packages())
if (any(!has)) install .packages(wants[!has])
library ("beanplot")
# Create target dataset
dat.g4 <???? data.frame(featVal=as.numeric(featVal.g4), FTB=FTB)
# boxplot
boxplot(featVal~FTB, data=dat.g4, xlab="TB group", ylab="Feature values")
# violin plot
beanplot(featVal~FTB, data=dat.g4, xlab="TB group", ylab="Feature values")
save.image( le ="HW1/RData/AllDataNeeded.RData")
#-------------------------------------
# Problem 3
#-------------------------------------
# Library calls
library (gplots)
best20 <???? order(myPvals)[1:20]
data.matrix <???? X[best20,]
# Labeled each column to the 3 group
colnames(data.matrix) <???? FTB
8
match(data.matrix[4,],featVal.g4) # check if the cols order matched
# Heatmap using actual feature values
jpeg( le ="HW1/Figures/Figure1.jpeg",width=1280,height=1024,pointsize=20)
heatmap.2(data.matrix,trace="none",main="Heatmap for Best20",key=T)
dev.o()
# The above heatmap seems less informative
# Try rank them across all the genes
data.matrix.rk <???? data.matrix
for (i in 1:nrow(data.matrix)) {
  data.matrix.rk[i,] = as.numeric(rank(data.matrix[i,]))
}
jpeg( le ="HW1/Figures/Figure1 rank.jpeg",width=1280,height=1024,pointsize=20)
heatmap.2(data.matrix.rk,trace="none",main="Heatmap for Best20 based on rank",key=T)
dev.o()
#-------------------------------------
# Problem 4
#-------------------------------------
# Library calls
source("http://bioconductor.org/biocLite.R")
biocLite("multtest")
require (multtest)
require (pROC)
# Dataset preparison from Dr. Gaile's lagecy code
# this provides an AffyBatch object
load("HW1/PSpikeData/PSpikeAffyBatch.RData")
spikeDF=read.table( le ="HW1/PSpikeData/AffyProbeSpikeValues.csv",sep="\t")
levels (spikeDF[,2])
summary(spikeDF[,2])
SpikeFC=as.numeric(levels(spikeDF[,2])[spikeDF[,2]])
names(SpikeFC)=spikeDF[,1]
nonZeroDX=which((SpikeFC!=0)&(!is.na(SpikeFC)))3422
nonDEdx = which(SpikeFC[nonZeroDX] == 1)
DEdx = which(SpikeFC[nonZeroDX] != 1)
require (affy)
require (affyPLM)
#------------------------------------------------
# Make our own function to generate the ROC
# Note:
# Here we can set 4 different route-parameter
# and can choose using the maximum statistisc
# or the minimum p-value for testing
#------------------------------------------------
roc.gen <???? function(bgcorr,norm,pmcorr,sum,minP=F) {
  exprVal <???? expresso(affydata, bgcorrect.method = bgcorr, normalize.method = norm,
                       pmcorrect.method = pmcorr, summary.method = sum)
  featVal <???? exprs(exprVal)[nonZeroDX, ]
  nfeat <???? dim(featVal)[1]
  group <???? factor(c(rep("A", 9), rep("B", 9)))
  9
  myresponse <???? rep(NA, nfeat)
  myresponse[nonDEdx] <???? 0
  myresponse[DEdx] <???? 1
  if (!minP) {
    resT <???? mt.maxT(featVal, group, B = 10000)
    testStat <???? resT$teststat[order(resT$index)]
    roc = roc(response = myresponse, predictor = abs(testStat))
  }
  if (minP) {
    resP <???? mt.minP(featVal, group, B = 10000)
    adjpVal <???? resP$adjp[order(resP$index)]
    roc = roc(response = myresponse, predictor = adjpVal)
  }
  return(roc)
}
# You can change parameters below, I'll recommend using minP method to test
roc.1 <???? roc.gen("mas","loess","pmonly","avgdiff",minP=T)
roc.2 <???? roc.gen("mas","constant","mas","avgdiff",minP=T)
roc.3 <???? roc.gen("mas","invariantset","pmonly","medianpolish",minP=T)
roc.4 <???? roc.gen("rma","loess","pmonly","avgdiff",minP=T)
roc.5 <???? roc.gen("rma","constant","mas","medianpolish",minP=T)
roc.6 <???? roc.gen("rma","invariantset","pmonly","avgdiff",minP=T)
roc.7 <???? roc.gen("none","loess","pmonly","avgdiff",minP=T)
roc.8 <???? roc.gen("none","constant","mas","medianpolish",minP=T)
# save(roc.1,roc.2,roc.3,roc.4,roc.5,roc.6,roc.7,roc.8,
# file="HW1/RData/LoadData.RData")
load("HW1/RData/LoadData.RData")
# Plot 8 ROC curves into one figure
jpeg( le ="HW1/Figures/Figure2.jpeg",width=1600,height=1200,pointsize=36)
layout(cbind(1, 2), widths=c(4, 1))
plot(roc.1, main = "ROCs for 8 routes", xlim=c(1,0), ylim=c(0,1))
lines (roc.2, col = "blue")
lines (roc.3, col = "green")
lines (roc.4, col = "red")
lines (roc.5, col = "green4")
lines (roc.6, col = "firebrick")
lines (roc.7, col = "cornflowerblue")
lines (roc.8, col = "darkblue")
par(mar=c(0,0,0,0))
plot.new()
legend("center", lty=c(1,1,1,1,1,1,1,1),lwd=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5),
       c("Route 1","Route 2","Route 3","Route 4","Route 5","Route 6","Route 7","Route 8"),
       col=c("black","blue","green","red","green4","firebrick","cornflowerblue","darkblue"),cex=0.8)
dev.o()
# Extract the AUC info from ROC curves
roc <???? list(roc.1,roc.2,roc.3,roc.4,roc.5,roc.6,roc.7,roc.8)
auc <???? c()
for (i in 1:8) {
  auc[i] <???? as.numeric(roc[[i]]$auc)
}
10
# Make a table and figure to visualize them
auc.tbl <???? as.matrix(auc)
rownames(auc.tbl) <???? c("Route 1","Route 2","Route 3","Route 4","Route 5","Route 6","Route 7","Route 8")
colnames(auc.tbl) <???? "AUC"
print(auc.tbl)
write.csv(auc.tbl, le="HW1/output/auc.csv")
jpeg( le ="HW1/Figures/Figure3.jpeg",width=800,height=600,pointsize=20)
plot(auc.tbl, type="b", ylab="AUC", xlab="Route")
dev.o()
11