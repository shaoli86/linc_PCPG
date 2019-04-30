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

clindataexp <- clindataexp[,-2]

#CART
library(rpart)
set.seed(2017)
n <- nrow(clindataexp)
sample <- sample(seq(n), size = n * 0.6, replace = FALSE)
trainexp <- clindataexp[sample,]
testexp <- clindataexp[-sample,]
fit <- rpart(mRNA.Subtype.Clusters ~ ., trainexp, method="class")
#fit <- rpart(Clinically.Aggressive.and.or.Metastatic ~ ., trainexp, method="class")
tp <- table(predict(fit, testexp, type = "class"), testexp[,1])
cartroc <- multiclass.roc(as.numeric(testexp[,1]), as.numeric(predict(fit, testexp, type = "class")))
cartauc <- cartroc$auc
plot(fit, uniform=TRUE, main="Classification Tree for PCPG Subtypes")
text(fit, use.n=TRUE, all=TRUE, cex=.8)

library(glmnet)
library(pROC)
#clindataexp <- na.omit(clindataexp)
clindataexp[is.na(clindataexp)] <- 0
set.seed(2017)
n <- nrow(clindataexp)
sample <- sample(seq(n), size = n * 0.6, replace = FALSE)
trainexp <- clindataexp[sample,]
testexp <- clindataexp[-sample,]

#Elastic Net 
a <- seq(0.1, 0.9, 0.05)
#Optimize lambda and alpha hyperparameters by 10-fold cross-validation
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(as.matrix(trainexp[,2:ncol(trainexp)]), as.numeric(trainexp[,1]), family = "gaussian", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(as.matrix(trainexp[,2:ncol(trainexp)]), as.numeric(trainexp[,1]), family = "gaussian", lambda = cv3$lambda.1se, alpha = cv3$alpha) #lambda=0.7884815, alpha=0.15
pathselect <- coef(md3)
pathpar <- as.matrix(pathselect)
write.table(pathpar, "C:\\Users\\dass8\\Documents\\NCI\\pacak\\lncRNA_PCPG_ElNet_MolSubtypes_countData.txt", sep="\t", quote=FALSE)
elnetroc <- multiclass.roc(as.numeric(testexp[,1]), as.numeric(predict(md3, as.matrix(testexp[,2:ncol(testexp)]), type = "response")))
elnetauc <- elnetroc$auc #elnetauc=0.9674

lincselect <- rownames(subset(as.data.frame(pathpar), s0!=0))
lincselect <- lincselect[-1]
write.table(lincselect, "SelectedlncRNA_PCPG_ElNet_MolSubtypes_countData.txt", sep="\t", quote=FALSE)

#Ridge 
#Optimize lambda hyperparameter by 10-fold cross-validation
cv <- cv.glmnet(as.matrix(trainexp[,2:ncol(trainexp)]), as.numeric(trainexp[,1]), family = "gaussian", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = 0)
md2 <- glmnet(as.matrix(trainexp[,2:ncol(trainexp)]), as.numeric(trainexp[,1]), family = "gaussian", lambda = cv3$lambda.1se, alpha = 0)
pathselect <- coef(md2)
pathpar <- as.matrix(pathselect)
write.table(pathpar, "lncRNA_PCPG_Ridge_countData.txt", sep="\t", quote=FALSE)
ridgeroc <- multiclass.roc(as.numeric(testexp[,1]), as.numeric(predict(md2, as.matrix(testexp[,2:ncol(testexp)]), type = "response")))
ridgeauc <- ridgeroc$auc
plot(ridgeroc)

#LASSO 
#Optimize lambda hyperparameter by 10-fold cross-validation
cv <- cv.glmnet(as.matrix(trainexp[,2:ncol(trainexp)]), as.numeric(trainexp[,1]), family = "gaussian", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = 1)
md1 <- glmnet(as.matrix(trainexp[,2:ncol(trainexp)]), as.numeric(trainexp[,1]), family = "gaussian", lambda = cv3$lambda.1se, alpha = 1)
pathselect <- coef(md1)
pathpar <- as.matrix(pathselect)
write.table(pathpar, "lncRNA_PCPG_LASSO_countData.txt", sep="\t", quote=FALSE)
lassoroc <- multiclass.roc(as.numeric(testexp[,1]), as.numeric(predict(md1, as.matrix(testexp[,2:ncol(testexp)]), type = "response")))
lassoauc <- lassoroc$auc



