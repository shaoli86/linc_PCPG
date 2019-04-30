clin <- read.table("TCGA_PCPG_clinical_molSubtype.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
for(i in 1:nrow(clin)){
alls <- unlist(strsplit(clin[i,1],"-"))
clin[i,1] <- paste(alls[1],alls[2],alls[3],alls[4],sep=".")
}
clinpseudo <- subset(clin, mRNA.Subtype.Clusters=="Pseudohypoxia")
clincortical <- subset(clin, mRNA.Subtype.Clusters=="Cortical admixture")
clinwnt <- subset(clin, mRNA.Subtype.Clusters=="Wnt-altered")
clinkinase <- subset(clin, mRNA.Subtype.Clusters=="Kinase signaling")
clinsdhb <- subset(clin, !is.na(SDHB.Germline.Mutation))
clinsdhd <- subset(clin, !is.na(SDHD.Germline.Mutation))
clinsdhx <- unique(c(clinsdhb[,1],clinsdhd[,1]))
clin <- subset(clin, mRNA.Subtype.Clusters=="Pseudohypoxia" | mRNA.Subtype.Clusters=="Cortical admixture" | mRNA.Subtype.Clusters=="Wnt-altered" | mRNA.Subtype.Clusters=="Kinase signaling" | (!is.na(SDHB.Germline.Mutation)) | (!is.na(SDHD.Germline.Mutation)), select=c("SampleID","mRNA.Subtype.Clusters","Clinically.Aggressive.and.or.Metastatic"))

linc <- read.table("gencode22_lincRNA_GeneID.txt", sep="\t", stringsAsFactors=FALSE)
pcpg <- read.table("PCPG_GDC_CountData_normal_nr.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
pcpg <- na.omit(pcpg)
rownames(pcpg) <- pcpg[,1]
pcpg <- pcpg[,-1]
clinpcpg <- colnames(pcpg)
spseudo <- clinpseudo[,1]
swnt  <- clinwnt[,1]
skinase  <- clinkinase[,1]
scortical <- clincortical[,1]
ssdhb <- clinsdhb[,1]
ssdhd <- clinsdhd[,1]
ssdhx <- unique(c(ssdhb,ssdhd))
snormal <- c("TCGA.P8.A5KD.11A","TCGA.SQ.A6I4.11A","TCGA.P8.A5KC.11A")
spseudo <- setdiff(spseudo,ssdhx)
pcpglinc <- subset(pcpg, rownames(pcpg) %in% linc[,1])
normalexpr <- subset(pcpglinc, select=names(pcpglinc) %in% snormal)
sdhxexpr <- subset(pcpglinc, select=names(pcpglinc) %in% ssdhx)
pseudoexpr <- subset(pcpglinc, select=names(pcpglinc) %in% spseudo)
wntexpr <- subset(pcpglinc, select=names(pcpglinc) %in% swnt)
kinaseexpr <- subset(pcpglinc, select=names(pcpglinc) %in% skinase)
corticalexpr <- subset(pcpglinc, select=names(pcpglinc) %in% scortical)
dataexp <- cbind(sdhxexpr,pseudoexpr, wntexpr, kinaseexpr, corticalexpr,normalexpr)
ms <- numeric(ncol(sdhxexpr))
nms <- numeric(ncol(pseudoexpr))
ws <- numeric(ncol(wntexpr))
ks <- numeric(ncol(kinaseexpr))
cs <- numeric(ncol(corticalexpr))
ns <- numeric(ncol(normalexpr))
for(i in 1:ncol(sdhxexpr))
ms[i] <- 1
for(i in 1:ncol(pseudoexpr))
nms[i] <- 2
for(i in 1:ncol(wntexpr))
ws[i] <- 3
for(i in 1:ncol(kinaseexpr))
ks[i] <- 4
for(i in 1:ncol(corticalexpr))
cs[i] <- 5
for(i in 1:ncol(normalexpr))
ns[i] <- 6
group <- factor(c(ms,nms,ws,ks,cs,ns))
y <- edgeR::DGEList(counts=dataexp,group=group)
y <- edgeR::calcNormFactors(y, na.rm=TRUE)
design <- model.matrix(~group)
y <- edgeR::estimateDisp(y, design)
logcpm <- edgeR::cpm(y, prior.count=2, log=TRUE)
#dataexp <- pcpglinc
#dataexp <- log2(dataexp + 1)
dataexp <- t(scale(t(logcpm)))
dataexp <- as.matrix(t(dataexp))

dataexp <- data.frame(SampleID=rownames(dataexp), as.matrix(dataexp))

clindataexp <- dplyr::inner_join(clin,dataexp)
clindataexp <- within(clindataexp, mRNA.Subtype.Clusters[SampleID %in% ssdhx] <- "SDHX_mutated")
rownames(clindataexp) <- clindataexp[,1]
clindataexp <- clindataexp[,-1]

clindataexp <- within(clindataexp, mRNA.Subtype.Clusters[mRNA.Subtype.Clusters == "SDHX_mutated"] <- 1)
clindataexp <- within(clindataexp, mRNA.Subtype.Clusters[mRNA.Subtype.Clusters == "Pseudohypoxia"] <- 2)
clindataexp <- within(clindataexp, mRNA.Subtype.Clusters[mRNA.Subtype.Clusters == "Wnt-altered"] <- 3)
clindataexp <- within(clindataexp, mRNA.Subtype.Clusters[mRNA.Subtype.Clusters == "Kinase signaling"] <- 4)
clindataexp <- within(clindataexp, mRNA.Subtype.Clusters[mRNA.Subtype.Clusters == "Cortical admixture"] <- 5)
clindataexp <- within(clindataexp, Clinically.Aggressive.and.or.Metastatic[Clinically.Aggressive.and.or.Metastatic == "Yes"] <- 1)
clindataexp <- within(clindataexp, Clinically.Aggressive.and.or.Metastatic[Clinically.Aggressive.and.or.Metastatic == "No"] <- 0)
clindataexp <- within(clindataexp, Clinically.Aggressive.and.or.Metastatic[Clinically.Aggressive.and.or.Metastatic == "no"] <- 0)

clindataexp <- clindataexp[,-1]

lincmeta <- read.table("SelectedlincRNAs_by_ElasticNet_for_PCPG_MolSubtypes.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
#clindataexp <- subset(clindataexp, select=colnames(clindataexp) %in% c("Clinically.Aggressive.and.or.Metastatic","SDHB.Germline.Mutation","ATRX.Somatic.Mutation","Tumor.Location","TERT",lincmeta[,1]))
clindataexp <- subset(clindataexp, select=colnames(clindataexp) %in% c("mRNA.Subtype.Clusters",lincmeta[,2]))
#clindataexp <- subset(clindataexp, select=colnames(clindataexp) %in% c("Clinically.Aggressive.and.or.Metastatic","SDHB.Germline.Mutation","ATRX.Somatic.Mutation","Tumor.Location","TERT")

clindataexp <- na.omit(clindataexp)
clindataexp[is.na(clindataexp)] <- 0

#Random Forest
library(randomForest)
library(dplyr)
library(ROCR)
nt <- 10
a <- 1
#aucr <- numeric(20)
#while(nt < 1000){
classtest <- numeric(nrow(clindataexp))
predicttest <- numeric(nrow(clindataexp))
k=1
clindataridge <- clindataexp
n <- nrow(clindataridge)
for(i in 1:nrow(clindataridge)){
trainexp <- clindataridge[-i,]
testexp <- clindataridge[i,]
cl <- as.numeric(trainexp[,1])
traindata <- as.matrix(trainexp[,-1])
testdata <- as.matrix(testexp[,-1])
model <- randomForest(traindata, cl, , ntree=200)
predicttest[k] <- as.numeric(predict(model, newdata = testdata, type="response"))
classtest[k] <- as.numeric(testexp[,1])
k <- k+1
}
rfroc <- roc(as.numeric(classtest),as.numeric(predicttest))
aucr <- rfroc$auc
#nt <- nt + 50
#a <- a+1
#}

correct <- 0
for(i in 1:length(classtest)){
if(classtest[i]==round(predicttest[i]))
correct=correct+1
}
accuracy <- correct/length(classtest)


#Neural network
if (!require("nnet")) {
   install.packages("nnet", dependencies = TRUE)
   library(nnet)
   }

classtest <- numeric(nrow(clindataexp))
predicttest <- numeric(nrow(clindataexp))
k=1
clindataridge <- clindataexp
n <- nrow(clindataridge)
for(i in 1:nrow(clindataridge)){
trainexp <- clindataridge[-i,]
testexp <- clindataridge[i,]
cl <- as.numeric(trainexp[,1])
traindata <- as.matrix(trainexp[,-1])
testdata <- as.matrix(testexp[,-1])
nmodel <- nnet(traindata, cl, size=4, maxit=10000)
predicttest[k] <- predict(nmodel, newdata=testdata, type="raw")
classtest[k] <- as.numeric(testexp[,1])
k <- k+1
}
nnroc <- roc(as.numeric(classtest),as.numeric(predicttest))
aucn <- nnroc$auc


#Support vector machine
if (!require("e1071")) {
   install.packages("e1071", dependencies = TRUE)
   library(e1071)
   }
classtest <- numeric(nrow(clindataexp))
predicttest <- numeric(nrow(clindataexp))
k=1
clindataridge <- clindataexp
n <- nrow(clindataridge)
for(i in 1:nrow(clindataridge)){
trainexp <- clindataridge[-i,]
testexp <- clindataridge[i,]
cl <- as.numeric(trainexp[,1])
traindata <- as.matrix(trainexp[,-1])
testdata <- as.matrix(testexp[,-1])
#svmmodel <- svm(traindata, cl, kernel="polynomial", degree=3, gamma=0.25)
svmmodel <- svm(traindata, cl, kernel="radial", gamma=0.1)
predicttest[k] <- predict(svmmodel, testdata, type="eps-regression")
classtest[k] <- as.numeric(testexp[,1])
k <- k+1
}
svmroc <- roc(as.numeric(classtest),as.numeric(predicttest))
aucs <- svmroc$auc

#Xtreme gradient boosting
library(xgboost)
classtest <- numeric(nrow(clindataexp))
predicttest <- numeric(nrow(clindataexp))
k=1
clindataridge <- clindataexp
n <- nrow(clindataridge)
for(i in 1:nrow(clindataridge)){
trainexp <- clindataridge[-i,]
testexp <- clindataridge[i,]
dtrain = xgb.DMatrix(as.matrix(trainexp[,-1]),label=as.numeric(trainexp[,1]))
y_train = as.numeric(trainexp[,1])
dtest = xgb.DMatrix(as.matrix(testexp[,-1]),label=as.numeric(testexp[,1]))
y_test = as.numeric(testexp[,1])
#params <- list(max_depth=2, eta=1, silent=1, objective='multi:softprob', num_class=3040)
watchlist <- list(train = dtrain, eval = dtest)
#cv <- xgb.cv(params=param, data = dnorm, nrounds = 2, nfold = 2, metrics = "mlogloss")
param <- list(max_depth = 4, eta = 1, silent = 1, nthread = 2, objective = "multi:softprob", num_class = 6)
bst <- xgb.train(param, dtrain, nrounds = 2, watchlist)
#cv <- xgb.cv(params=param, data = dnorm, nrounds = 100, early_stopping_rounds=50, nfold = 10, prediction=TRUE, metrics = "mlogloss")
predicttest[k] <- as.numeric(predict(bst, dtest, reshape=TRUE))
classtest[k] <- as.numeric(testexp[,1])
k <- k+1
}
xgbroc <- roc(as.numeric(classtest),as.numeric(predicttest))
aucx <- svmroc$auc

correct <- 0
for(i in 1:length(classtest)){
if(classtest[i]==round(predicttest[i]))
correct=correct+1
}
accuracy <- correct/length(classtest)

png("ROC_metasatic_PCPGlinc_multimodel4.png",height=2200, width=2600, res=300)
plot(rfroc, col = "red", cex=2)
plot(xgbroc, add = TRUE, col = "blue", cex=2)
plot(svmroc, add = TRUE, col = "darkgreen", cex=2)
#plot(perfaucsvmradial[[1]], add = TRUE, col = "darkgreen", cex=2)
legend("bottomright", legend=c(paste("RF, AUC:",round(aucr, 2)),paste("XGB, AUC:",round(aucxgb, 2)),paste("SVM, AUC:",round(aucs, 2))), col=c("red","blue","darkgreen"), pch=1)
dev.off()

