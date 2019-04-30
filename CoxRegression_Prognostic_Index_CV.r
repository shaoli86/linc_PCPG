cancer="PCPG"
linc <- read.table("gencode22_lincRNA_GeneID.txt", sep="\t", stringsAsFactors=FALSE)
pcpg <- read.table("PCPG_GDC_CountData_normal_nr.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
pcpg <- na.omit(pcpg)
rownames(pcpg) <- pcpg[,1]
pcpg <- pcpg[,-1]
clinpcpg <- colnames(pcpg)
dataexp <- pcpg
y <- edgeR::DGEList(counts=dataexp)
y <- edgeR::calcNormFactors(y, na.rm=TRUE)
logcpm <- edgeR::cpm(y, prior.count=2, log=TRUE)
dataexp <- t(scale(t(logcpm)))
mycgds = cgdsr::CGDS("http://www.cbioportal.org/")
mycancerstudy <- paste(tolower(cancer),"tcga",sep="_")
mycaselist = cgdsr::getCaseLists(mycgds,mycancerstudy)[1,1]
#datac <- try(cgdsr::getClinicalData(mycgds,mycaselist))
clin <- read.table("TCGA_PCPG_clinical_molSubtype.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
for(i in 1:nrow(clin)){
alls <- unlist(strsplit(clin[i,1],"-"))
clin[i,1] <- paste(alls[1],alls[2],alls[3],alls[4],sep=".")
}
clin <- within(clin, SDHB.Germline.Mutation[!is.na(SDHB.Germline.Mutation)] <- 1.0)
clin <- within(clin, SDHB.Germline.Mutation[is.na(SDHB.Germline.Mutation)] <- 0.0)
clin <- within(clin, ATRX.Somatic.Mutation[ATRX.Somatic.Mutation != 0] <- 1.0)
clin <- within(clin, ATRX.Somatic.Mutation[ATRX.Somatic.Mutation == 0] <- 0.0)
clin <- within(clin, Tumor.Location[Tumor.Location != "adrenal"] <- 1.0)
clin <- within(clin, Tumor.Location[Tumor.Location == "adrenal"] <- 0.0)
clin <- within(clin, Normetanephrine.secreting[Normetanephrine.secreting == "Yes"] <- 1.0)
clin <- within(clin, Normetanephrine.secreting[Normetanephrine.secreting != 1.0] <- 0.0)
#clin <- within(clin, Normetanephrine.secreting[is.na(Normetanephrine.secreting)] <- 0.0)
clin <- within(clin, Norepinephrine.Secreting[Norepinephrine.Secreting == "Yes"] <- 1.0)
clin <- within(clin, Norepinephrine.Secreting[Norepinephrine.Secreting != 1.0] <- 0.0)
#clin <- within(clin, Norepinephrine.Secreting[is.na(Norepinephrine.Secreting)] <- 0.0)
clin <- within(clin, Epinephrine.Secreting[Epinephrine.Secreting == "Yes"] <- 1.0)
clin <- within(clin, Epinephrine.Secreting[Epinephrine.Secreting != 1.0] <- 0.0)
#clin <- within(clin, Epinephrine.Secreting[is.na(Epinephrine.Secreting)] <- 0.0)
clin <- within(clin, Metanephrine.Secreting[Metanephrine.Secreting == "Yes"] <- 1.0)
clin <- within(clin, Metanephrine.Secreting[Metanephrine.Secreting != 1.0] <- 0.0)
#clin <- within(clin, Metanephrine.Secreting[is.na(Metanephrine.Secreting)] <- 0.0)
clin <- within(clin, Methoxytyramine.Secreting[Methoxytyramine.Secreting == "Yes"] <- 1.0)
clin <- within(clin, Methoxytyramine.Secreting[Methoxytyramine.Secreting != 1.0] <- 0.0)
#clin <- within(clin, Methoxytyramine.Secreting[is.na(Methoxytyramine.Secreting)] <- 0.0)
clin <- within(clin, Dopamine.Secreting[Dopamine.Secreting == "Yes"] <- 1.0)
clin <- within(clin, Dopamine.Secreting[Dopamine.Secreting != 1.0] <- 0.0)
#clin <- within(clin, Dopamine.Secreting[is.na(Dopamine.Secreting)] <- 0.0)

clin$SDHB.Germline.Mutation <- as.numeric(clin$SDHB.Germline.Mutation)
clin$ATRX.Somatic.Mutation <- as.numeric(clin$ATRX.Somatic.Mutation)
clin$Tumor.Location <- as.numeric(clin$Tumor.Location)
clin$Normetanephrine.secreting <- as.numeric(clin$Normetanephrine.secreting)
clin$Norepinephrine.Secreting <- as.numeric(clin$Norepinephrine.Secreting)
clin$Epinephrine.Secreting <- as.numeric(clin$Epinephrine.Secreting)
clin$Metanephrine.Secreting <- as.numeric(clin$Metanephrine.Secreting)
clin$Methoxytyramine.Secreting <- as.numeric(clin$Methoxytyramine.Secreting)
clin$Dopamine.Secreting <- as.numeric(clin$Dopamine.Secreting)
clin$Normetanephrine.secreting <- as.numeric(clin$Normetanephrine.secreting) + as.numeric(clin$Norepinephrine.Secreting)
clin$Normetanephrine.secreting[clin$Normetanephrine.secreting > 1.0] <- 1.0
clin$Normetanephrine.secreting <- as.numeric(clin$Normetanephrine.secreting)
clin$Epinephrine.Secreting <- as.numeric(clin$Epinephrine.Secreting) + as.numeric(clin$Metanephrine.Secreting)
clin$Epinephrine.Secreting[clin$Epinephrine.Secreting > 1.0] <- 1.0
clin$Epinephrine.Secreting <- as.numeric(clin$Epinephrine.Secreting)
clin$Methoxytyramine.Secreting <- as.numeric(clin$Methoxytyramine.Secreting) + as.numeric(clin$Dopamine.Secreting)
clin$Methoxytyramine.Secreting[clin$Methoxytyramine.Secreting > 1.0] <- 1.0
clin$Methoxytyramine.Secreting <- as.numeric(clin$Methoxytyramine.Secreting)
clinpseudo <- subset(clin, mRNA.Subtype.Clusters=="Pseudohypoxia")
clincortical <- subset(clin, mRNA.Subtype.Clusters=="Cortical admixture")
clinwnt <- subset(clin, mRNA.Subtype.Clusters=="Wnt-altered")
clinkinase <- subset(clin, mRNA.Subtype.Clusters=="Kinase signaling")
clinsdhb <- subset(clin, !is.na(SDHB.Germline.Mutation))
clinsdhd <- subset(clin, !is.na(SDHD.Germline.Mutation))
clinret <- subset(clin, !is.na(RET.Germline.Mutation))
clinretsom <- subset(clin, RET.Somatic.Mutation!=0)
clinret <- unique(rbind(clinret,clinretsom))
clinnf1 <- subset(clin, !is.na(NF1.Germline.Mutation))
clinnf1som <- subset(clin, NF1.Somatic.Mutation!=0)
clinnf1 <- unique(rbind(clinnf1,clinnf1som))
clinvhl <- subset(clin, !is.na(VHL.Germline.Mutation))
clinvhlsom <- subset(clin, VHL.Somatic.Mutation!=0)
clinvhl <- unique(rbind(clinvhl,clinvhlsom))
clinepas1 <- subset(clin, EPAS1.Somatic.Mutation!=0)
clinhras <- subset(clin, HRAS.Somatic.Mutation!=0)
colnames(clin)[26] <- "DFS_MONTHS"
clin[,26] <- as.numeric(clin[,26])/30
colnames(clin)[16] <- "DFS_STATUS"
clin <- within(clin, DFS_STATUS[DFS_STATUS == 'yes' | DFS_STATUS == 'Yes'] <- 'Recurred/Progressed')
clin <- within(clin, DFS_STATUS[DFS_STATUS == 'no' | DFS_STATUS == 'No'] <- 'DiseaseFree')
colnames(clin)[6] <- "AGE"
rownames(clin) <- clin[,1]
clin <- clin[,-1]
datac <- data.frame(clin,sampleId=rownames(clin))
mymutprofile = try(paste(tolower(cancer),"_tcga_mutations",sep=""))
survdata <- data.frame(Gene1="",Gene2="",SurvPval=0)
#survdata <- read.table("C:\\Users\\dass8\\Documents\\NCI\\pacak\\prognostic_linc3.txt", sep="\t", stringsAsFactors=FALSE)
survdata <- read.table("Clinical_Aggressive_lincRNA.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
survdata2 <- read.table("SDHX-mut_Pseudohypoxia_PCPG_lincRNA_Survival.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
survdata3 <- read.table("Pseudohypoxia_signaling_PCPG_lincRNA_Survival.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
survdata4 <- read.table("Kinase_signaling_PCPG_lincRNA_Survival.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
survlinc <- unique(c(survdata[,1],survdata2[,2],survdata3[,2],survdata4[,2]))
#survlinc <- unique(c(survdata2[,2],survdata3[,2],survdata4[,2]))
#survlinc <- survdata[,1]
data2 <- subset(dataexp, rownames(dataexp) %in% c("TERT",survlinc))
data2 <- t(data2)
genemut1 <- unique(c(clinsdhb[,1],clinsdhd[,1]))
#gene2mutexp <- subset(data2, row.names(data2) %in% genemut1)
gene2mutexp <- data2
gene2mutexp <- na.omit(gene2mutexp)
gene2mutexp <- data.frame(gene2mutexp, sampleId=rownames(gene2mutexp))
clinical <- subset(datac, select=c(DFS_STATUS,DFS_MONTHS,AGE,sampleId,SDHB.Germline.Mutation,ATRX.Somatic.Mutation,Tumor.Location,Normetanephrine.secreting,Epinephrine.Secreting,Methoxytyramine.Secreting))
clinical <- dplyr::inner_join(clinical,gene2mutexp)
rownames(clinical) <- clinical$sampleId
clinical <- clinical[,-4]
clinical <- subset(clinical, DFS_STATUS!="")
clinical <- within(clinical, DFS_MONTHS[DFS_MONTHS == '[Not Available]'] <- '0')
clinical <- within(clinical, DFS_MONTHS[DFS_MONTHS == 'NA'] <- '0')
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == '[Not Available]'] <- 'NA')
clinical <- subset(clinical, DFS_STATUS!="NA")
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == 'DiseaseFree'] <- 0)
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == 'Recurred/Progressed'] <- 1)
pvalpi <- numeric(100)
samplelist <- list()
for(r in 1:100){
n <- nrow(clinical)
sample <- sample(seq(n), size = n * 0.5, replace = FALSE)
clinical60 <- clinical[-sample,]
clinical40 <- clinical[sample,]
samplelist[[r]] <- sample
f1 <- as.formula(paste("survival::Surv(as.numeric(DFS_MONTHS), as.numeric(DFS_STATUS)) ~", paste(colnames(clinical40)[4:ncol(clinical40)], collapse= "+")))
#f1 <- as.formula(paste("survival::Surv(as.numeric(DFS_MONTHS), as.numeric(DFS_STATUS)) ~", paste(c(colnames(clinical)[4:12],"TERT"), collapse= "+")))
coxdiff <- survival::coxph(f1, data=clinical40)
#survminer::ggforest(coxdiff,data=clinical40)

scox <- survival::survfit(coxdiff)
genecoef <- numeric(length(coxdiff$coef))
pigene <- numeric(nrow(clinical40))
for(k in 1:nrow(clinical40)){
pigene[k]=0
for(i in 1:length(coxdiff$coef)){
genecoef[i] <- as.numeric(coxdiff$coef[i])
j <- i+3
pigene[k] <- pigene[k] + clinical40[k,j]*genecoef[i]
}
}
pigene <- scale(pigene)[,1]
#gene2coef <- as.numeric(coxdiffgenecoef[i]
#gene3coef <- as.numeric(coxdiff$coef[3])
#gene4coef <- as.numeric(coxdiff$coef[4])
#gene5coef <- as.numeric(coxdiff$coef[5])
#gene6coef <- as.numeric(coxdiff$coef[6])
#gene7coef <- as.numeric(coxdiff$coef[7])
#gene8coef <- as.numeric(coxdiff$coef[8])
#gene9coef <- as.numeric(coxdiff$coef[9])
#gene10coef <- as.numeric(coxdiff$coef[10])

##Calculation of prognostic index PI= coef_gene1*exp_gene1 + coef_gene2*exp_gene2
clinical40 <- data.frame(clinical40, PI = pigene)
#clinical <- data.frame(clinical, PI = scale((clinical[,4]*gene1coef)+ (clinical[,5]*gene2coef)+ (clinical[,6]*gene3coef)+ (clinical[,7]*gene4coef)+ (clinical[,8]*gene5coef)+ (clinical[,9]*gene6coef))[,1])
##Stratifying samples on PI>0:Good_prognostic_index and PI<=0:Bad_prognostic_index
piup <- row.names(subset(clinical40, PI> 1))
pidown <- row.names(subset(clinical40, PI< -1))
picon <- character(nrow(clinical40))
for(i in 1:nrow(clinical40)){
if(row.names(clinical40)[i] %in% piup)
picon[i] <- "high> 1"
else if(row.names(clinical40)[i] %in% pidown)
picon[i] <- "low< -1"
else
picon[i] <- "NA"
}
clinical40 <- data.frame(clinical40, prognostic_index=picon, stringsAsFactors=FALSE)
clinical40 <- subset(clinical40, prognostic_index!="NA")
pigood <- nrow(subset(clinical40, prognostic_index=="high> 1"))
pibad <- nrow(subset(clinical40, prognostic_index=="low< -1"))
if(pigood>1 & pibad>1){
diffpi <- try(survival::survdiff(survival::Surv(as.numeric(DFS_MONTHS), as.numeric(DFS_STATUS)) ~ prognostic_index, data=clinical40))
pvalpi[r] <- try(pchisq(diffpi$chisq, length(diffpi$n)-1, lower.tail=FALSE))
#clinical <- within(clinical, prognostic_index[prognostic_index == "high> 1"] <- paste("high> 1:",pigood,sep=""))
#clinical <- within(clinical, prognostic_index[prognostic_index == "low< -1"] <- paste("low< -1:",pibad,sep=""))
mutsurv <- try(survival::survfit(survival::Surv(as.numeric(DFS_MONTHS), as.numeric(DFS_STATUS)) ~ prognostic_index, data=clinical40))
p <- try(survminer::ggsurvplot(mutsurv, pval=round(pvalpi[r],8), risk.table=TRUE, pval.size=12, font.tickslab=c(18,"plain","black")))
#png(paste("C:\\Users\\dass8\\Documents\\NCI\\pacak\\DFS_",gene1,"_PrognosticLinc6_",cancer,".png",sep=""), width=2000, height=800)
png(paste("C:\\Users\\dass8\\Documents\\NCI\\pacak\\pi_cross_validate\\DFS_PrognosticLinc_linc18_clinical_hormonesAll_",cancer,"_",r,".png",sep=""), width=3600, height=1500, res=300)
par(mar=c(4,3,4,6)+0.1)
print(p)
dev.off()
}
clinical <- subset(datac, select=c(DFS_STATUS,DFS_MONTHS,AGE,sampleId,SDHB.Germline.Mutation,ATRX.Somatic.Mutation,Tumor.Location,Normetanephrine.secreting,Epinephrine.Secreting,Methoxytyramine.Secreting))
clinical <- dplyr::inner_join(clinical,gene2mutexp)
rownames(clinical) <- clinical$sampleId
clinical <- clinical[,-4]
clinical <- subset(clinical, DFS_STATUS!="")
clinical <- within(clinical, DFS_MONTHS[DFS_MONTHS == '[Not Available]'] <- '0')
clinical <- within(clinical, DFS_MONTHS[DFS_MONTHS == 'NA'] <- '0')
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == '[Not Available]'] <- 'NA')
clinical <- subset(clinical, DFS_STATUS!="NA")
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == 'DiseaseFree'] <- 0)
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == 'Recurred/Progressed'] <- 1)
}
