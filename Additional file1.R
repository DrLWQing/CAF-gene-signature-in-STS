####GSE21122 data preparation####
library(tidyverse)
library(GEOquery)
chooseBioCmirror()
gset <- getGEO("GSE21122",destdir = ".",AnnotGPL = F,getGPL = F)
pdata <- pData(gset[[1]])
pdata <- pdata[9:18,]
library(stringr)
group_list <- ifelse(str_detect(pdata$source_name_ch1,"sarcoma"),
                     "tumor","normal")
group_list <- factor(group_list,levels = c("tumor","normal"))
colData <- data.frame(row.names = rownames(pdata),group_list)
write.table(colData,file = "GSE21122_conditions.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
exp <- exprs(gset[[1]])
boxplot(exp,outline=F,notch=T,col=group_list,las=2)
library(limma)
exp <- normalizeBetweenArrays(exp)
boxplot(exp,outline=F,notch=T,col=group_list,las=2)
range(exp)
index <- gset[[1]]@annotation
library("hgu133a.db")
ls("package:hgu133a.db")
ids <- toTable(hgu133aSYMBOL)
library(tidyverse)
exp <- as.data.frame(exp)
exp <- exp %>% mutate(probe_id=rownames(exp))
exp <- exp %>% inner_join(ids,by="probe_id")
exp <- exp[!duplicated(exp$symbol),]
rownames(exp) <- exp$symbol
exp <- exp[,-(159:160)]
write.table(exp,file="GSE21122_expr_symbol.txt",sep = "\t",row.names = T)
group <- read.table("GSE21122_conditions.txt",sep = "\t",stringsAsFactors = T,header = T,row.names = 1,check.names = F)
group_list <- group$group_list
exp <- read.table("GSE21122_expr_symbol.txt",sep = "\t",header = T,row.names = 1,stringsAsFactors = F,check.names = F)
library(limma)
design <- model.matrix(~group_list)
fit <- lmFit(exp,design)
fit <- eBayes(fit)
deg <- topTable(fit,coef = 2,number = Inf)
degs <- deg %>% dplyr::filter(abs(logFC)>1,adj.P.Val<0.05)
write.csv(degs,file = "GSE21122_DEGS_tumor_vs_normal.csv")
####enrichment nalysis####
library(tidyverse)
library(clusterProfiler)
library(BiocManager)
library(org.Hs.eg.db)
DEG <- read.table("GSE21122_DEGs_p_tumor_vs_normal.txt",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
gene <- read.table("SARC related CAF genes.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
DEG <- DEG[gene$Genes,]
DEG <- DEG %>% rownames_to_column('Gene')
genelist <- bitr(DEG$Gene,fromType = 'SYMBOL',
                 toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c('Gene'='SYMBOL'))
#GO
ego <- enrichGO(gene=DEG$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = 'all',
                pAdjustMethod = 'BH',
                minGSSize = 1,
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = T)
ego_res <- ego@result
dotplot(ego,showCategory = 10,split='ONTOLOGY')+facet_grid(ONTOLOGY~.,scale='free')
#KEGG
kk <- enrichKEGG(gene = DEG$ENTREZID,
                 organism = 'hsa',
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
kk_res <- kk@result
dotplot(kk,showCategory=20)
####TCGA-SARC data prepartion####
library(tidyverse)
fpkm1 = read.table(file = 'TCGA-SARC.htseq_fpkm.tsv', sep = '\t', header = TRUE) 
rownames(fpkm1) <- fpkm1[,1]  
fpkm1 = fpkm1[,-1]
table(substr(colnames(fpkm1),14,16))
fpkm1 <- fpkm1[,substr(colnames(fpkm1),14,16)%in% c("01A","11A")]
table(substr(colnames(fpkm1),14,16))
rownames(fpkm1) <- substr(rownames(fpkm1),1,15)
fpkm <- fpkm1
Ginfo_0 <- read.table("gene_length_Table.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Ginfo <- Ginfo_0[which(Ginfo_0$genetype == "protein_coding"),]
comgene <- intersect(rownames(fpkm),rownames(Ginfo))
fpkm <- fpkm[comgene,]
Ginfo <- Ginfo[comgene,]
fpkm$Gene <- as.character(Ginfo$genename)  
fpkm <- fpkm[!duplicated(fpkm$Gene),] 
rownames(fpkm) <- fpkm$Gene
fpkm <- fpkm[,-ncol(fpkm)] 
write.table(fpkm, file = "TCGA_SARC_fpkm_mRNA_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
tumor <- colnames(fpkm)[substr(colnames(fpkm),14,16) == "01A"]
fpkm_01A <- fpkm[,tumor]
write.table(fpkm_01A, file = "TCGA_SARC_fpkm_mRNA_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

library(tidyverse)
clini <- read.csv2('sarc_tcga_clinical_data.tsv',header = T,sep = '\t')
pheno <- clini[,c(3,4,6,13,14,21,22,24,32,33,35,36,37,54,55,56,62,69,70,74,76,78,88)]
names(pheno)[4] <- 'DFS.time'
pheno$DFS <- ifelse(substr(pheno$DFS.status,3,13)=='DiseaseFree','0','1')
pheno$OS <- ifelse(substr(pheno$OS.status,3,8)=='LIVING','0','1')
write.table(pheno, file = "TCGA_SARC_phenotype_raw.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
pheno <- read.table('TCGA_phenotype_raw.txt',sep = '')
DFS <- pheno[,c(1,3,24,4)]
DFS <- DFS[complete.cases(DFS),]
OS <- pheno[,c(1,3,25,15)]
OS$sample <- gsub("-",".",OS$sample)
DFS$sample <- gsub("-",".",DFS$sample)
write.table(DFS, file = "TCGA_SARC_DFS.txt",sep = "\t",quote = F)
write.table(OS, file = "TCGA_SARC_OS.txt",sep = "\t",quote = F)

####GSE30929 data preparation####
#similar as GSE21122

####uniCOX in TCGA-SARC####
library(tidyverse)
library(survival)
surv <- read.table(file = "TCGA_SARC_DFS.txt",header = T,row.names = ,sep = "\t")
surv <- surv[!duplicated(surv$sample),]
rownames(surv) <- surv$sample
surv <- surv[,3:4]
expr <-read.table(file = "TCGA_SARC_fpkm_mRNA_01A.txt",row.names= 1,header = T,sep = "\t",check.names = F)
colnames(expr) <- substr(colnames(expr),1,15)
comsample <- intersect(colnames(expr),rownames(surv))
expr <- expr[,comsample]
surv <- surv[comsample,]
mygene <- read.table(file = "SARC related CAF genes.txt",header = T,sep = "\t",check.names = F)
CAF.expr <- expr[mygene$Gene,]行
CAF.expr <- CAF.expr[complete.cases(CAF.expr),]
CAF.expr <- t(CAF.expr)
CAF.expr <- as.data.frame(CAF.expr)
surv.expr <- cbind(surv,CAF.expr)
Coxoutput <- NULL 
for(i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
  cox <- coxph(Surv(DFS.time,DFS) ~ surv.expr[,i], data = surv.expr)
  coxSummary = summary(cox)
  
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}
write.csv(Coxoutput, file = "uniCOX_TCGA_DFS_CAFgene100.csv",quote = F,row.names = T)
write.table(Coxoutput, file = "uniCOX_TCGA_DFS_CAFgene100.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####multiCOX in TCGA-SARC####
library(tidyverse)
library(survival)
library(caret)
library(pacman)
surv <- read.table(file = "TCGA_SARC_DFS.txt",header = T,row.names = ,sep = "\t")
surv <- surv[!duplicated(surv$sample),]
rownames(surv) <- surv$sample
surv <- surv[,3:4]
expr <-read.table(file = "TCGA_SARC_fpkm_mRNA_01A.txt",row.names= 1,header = T,sep = "\t",check.names = F)
colnames(expr) <- substr(colnames(expr),1,15)
comsample <- intersect(colnames(expr),rownames(surv))
expr <- expr[,comsample]
surv <- surv[comsample,]
mygene <- read.table(file = "uniCOX_ressig_DFS.txt",header = T,sep = "\t",check.names = F)
CAF.expr <- expr[mygene$gene,]
CAF.expr <- CAF.expr[complete.cases(CAF.expr),]
CAF.expr <- t(CAF.expr)
CAF.expr <- as.data.frame(CAF.expr)
surv.expr <- cbind(surv,CAF.expr)
#multiCOX
multiCox <- coxph(Surv(DFS.time,DFS)~.,data = surv.expr)
multiCox <- step(multiCox,direction = "both")
multiCoxSum <- summary(multiCox)
outMultiTab <- data.frame()
outMultiTab <- cbind(coef=multiCoxSum$coefficients[,'coef'],
                     HR=multiCoxSum$conf.int[,'exp(coef)'],
                     HR.95L=multiCoxSum$conf.int[,'lower .95'],
                     HR.95H=multiCoxSum$conf.int[,'upper .95'],
                     pvalue=multiCoxSum$coefficients[,'Pr(>|z|)'])
outMultiTab <- cbind(id = row.names(outMultiTab),outMultiTab)
#riskScore
riskScore <- predict(multiCox,type='risk',newdata=surv.expr)
coxGene <- rownames(multiCoxSum$coefficients)
risk <- as.vector(ifelse(riskScore>median(riskScore),'high','low'))
outCol <- c('DFS.time','DFS',coxGene)
dat <- cbind(surv.expr[,outCol],riskScore=as.vector(riskScore),risk)
write.csv(dat, file = "multiCOX3_rickScore_TCGA_DFS_CAFgene5.csv",quote = F,row.names = T)
write.table(dat, file = "multiCOX3_rickScore_TCGA_DFS_CAFgene5.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
####PCA of TCGA-SARC#### 
library(scatterplot3d)      
rt=read.table('multiCOX3_rickScore_TCGA_DFS_CAFgene5.txt',header=T,sep="\t",check.names=F,row.names=1)
data=rt[,3:7]
df_pca <- prcomp(data,scale. = T)
df_pcs <- data.frame(df_pca$x,risk=rt$risk2)
head(df_pcs)
group <- as.vector(rt[,'risk2'])
color <- ifelse(group=='low',"#3182BDFF","#FD8D3CFF")
scatterplot3d(df_pcs[,1:3],pch = 16,color = color,angle = 45)
legend('top',legend = c('low risk','high risk'),pch = 16,
       inset = -0.2,box.col = 'white',xpd=T,horiz = T,col = c("#3182BDFF","#FD8D3CFF"))
####survival analysis of TCGA-SARC####
surv_cluster<- read.table("multiCOX3_rickScore_TCGA_DFS_CAFgene5.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
colnames(surv_cluster) [1]<-"sample" 
surv_cluster <- surv_cluster[,c(1,2,3,9,10)]
surv_cluster$risk <- factor(surv_cluster$risk,levels = c("low",'high'))
class(surv_cluster$risk)
library(survival)
fitd <- survdiff(Surv(DFS.time,DFS)~risk,
                 data = surv_cluster,
                 na.action = na.exclude)
pValue <- 1-pchisq(fitd$chisq,length(fitd$n)-1)
fit <- survfit(Surv(DFS.time, DFS)~ risk, data = surv_cluster)
summary(fit)
library(survminer)
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
ggsurvplot(fit,
           data = surv_cluster,
           pval = p.lab,
           conf.int = T, 
           risk.table = TRUE,
           risk.table.col = "strata",
           palette = c("#3182BDFF", "#E6550DFF"),
           legend.labs = c("Low", "High"), 
           size = 1,
           xlim = c(0,140), 
           break.time.by = 20, 
           legend.title = "",
           surv.median.line = "hv", 
           ylab = "Survival probability (%)", 
           xlab = "Time (Months)",
           ncensor.plot = F,
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
####ROC of TCGA-SARC####
#ROC
library(ROCR)
library(rms)
expr_surv <- read.table("multiCOX3_rickScore_TCGA_DFS_CAFgene5.txt",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
ROC1 <- prediction(expr_surv$riskScore,expr_surv$DFS)
ROC2 <- performance(ROC1,"tpr","fpr")
AUC <- performance(ROC1,"auc")
AUC <- 0.7293115
plot(ROC2,
     col="#E6550DFF",
     xlab="False positive rate", ylab="True positive rate",
     lty=1,lwd=3,
     main="AUC=0.7293115")
abline(0, 1, lty=2, lwd=3)
#time ROC 
library(timeROC)
library(survival)
library(tidyverse)
exp_sur <- read.table("multiCOX3_rickScore_TCGA_DFS_CAFgene5.txt",sep = "\t",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
exp_sur$DFS.time <- exp_sur$DFS.time/12
ROC3 <- timeROC(T=exp_sur$DFS.time,  
                delta=exp_sur$DFS, 
                marker=exp_sur$riskScore,
                cause=1,
                weighting="marginal", 
                times=c(1, 3, 5),
                iid=TRUE)
ROC3
plot(ROC3,
     time=1, col="red")
plot(ROC3,
     time=3, col="green", add=TRUE)
plot(ROC3,
     time=5, col="blue", add=TRUE)
legend("bottomright",
       c("Year-1", "Year-3", "Year-5"),
       col=c("red", "green", "blue"),
       lty=1, lwd=2)
#multi-index ROC
library(pROC)
library(glmnet)
data <- read.table("TCGA_phenotype_data.txt",sep = "\t",header = T,row.names=1)
data$risk <- factor(data$risk,levels = c('low','high'))
data$sex <- factor(data$sex,levels = c('Female','Male'))
data$margin.status <- factor(data$margin.status,levels = c('R1','R0'))
data$DFS <- factor(data$DFS,levels = c('1','0'))
fit1 <- glm(DFS ~ riskScore + age + margin.status + TMB,
            data=data,
            family = binomial())  
summary(fit1)
data$prob <- predict(fit1,newdata = data,type = "response")
ROC1 <- roc(data$DFS,data$prob,levels=c('1','0'))
ROC2 <- roc(data$DFS,data$riskScore)
ROC3 <- roc(data$DFS,data$margin.status)
ROC4 <- roc(data$DFS,data$age)
ROC5 <- roc(data$DFS,data$TMB)
####immune infiltration in TCGA-SARC####
library(immunedeconv)
library(ggplot2)
library(tidyverse)
exprMatrix <- read.table(file = "TCGA_SARC_fpkm_mRNA_01A.txt",row.names=1,sep = "\t",header = T,stringsAsFactors = F)
colnames(exprMatrix) <- substr(colnames(exprMatrix),1,15)
risk_res <- read.table(file = "multiCOX3_rickScore_TCGA_DFS_CAFgene5.txt",row.names = 1, sep = "\t",header = T,stringsAsFactors = F)
comsample <- intersect(colnames(exprMatrix),rownames(risk_res))
exprMatrix <- exprMatrix[,comsample]
write.table(exprMatrix, "TCGA_SARC_FPKM_230.txt", sep="\t", col.names=T, quote=F)
#quantiseq
res1 <- deconvolute(exprMatrix, method="quantiseq")
write.table(res1, "quantiseq_GSE30929.txt", sep="\t", col.names=T, row.names=F, quote=F)
#epic
res2 <- deconvolute(exprMatrix, method="epic")
write.table(res2, "epic_GSE30929.txt", sep="\t", col.names=T, row.names=F, quote=F)
#xcell
library('xCell')
library("ggplot2")
library("ggpubr")
library("tidyverse")
exp <- read.table(file = "TCGA_SARC_FPKM_230.txt",row.names = 1, sep = "\t",header = T,stringsAsFactors = F)
colnames(exp) <- substr(colnames(exp),1,15)
risk_res <- read.table(file = "multiCOX3_rickScore_TCGA_DFS_CAFgene5.txt",row.names = 1, sep = "\t",header = T,stringsAsFactors = F)
comsample <- intersect(colnames(exp),rownames(risk_res))
exp <- exp[,comsample]
celltypeuse <- xCell.data$spill$K
rs <- xCellAnalysis(exp,parallel.sz = 10)
rs <- rs %>% t() %>% as.data.frame()
group <- read.table(file = "group_rickScore_TCGA.txt",sep = "\t",row.names=1 ,header = T,stringsAsFactors = T)
class(group$group)
group <- as.data.frame(group)
identical(rownames(group),rownames(rs))
rs$group <- group$group
class(rs$group)
rs <- rs %>% rownames_to_column("sample")
a <- rs
write.table(rs,"xCell_TCGA.txt",sep = "\t",row.names = T,col.names = T,quote = F)
####chemosensitivity prediction in TCGA-SARC####
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
GDSC2_Expr = readRDS(file=file.path('GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path("GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 
testExpr<- read.table('TCGA_SARC_FPKM_230.txt',sep = "\t",row.names = 1, header=T,stringsAsFactors = F)
testExpr <- as.matrix(testExpr)
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )
####GSEA of TCGA-SARC####
library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)
data <- read.table("DEG_TCGA_high_vs_low.txt",sep = "\t",header = T,row.names=1)
data <- data %>% rownames_to_column('row')
as.numeric(data$logFC)
geneList <- data$logFC
names(geneList) <- data$row
head(geneList)
geneList = sort(geneList, decreasing = TRUE)
genesets <- read.gmt('h.hallmark.all.v2022.1.Hs.symbols.gmt')
y <- GSEA(geneList,TERM2GENE =genesets,nPerm = 10000,pvalueCutoff = 0.25)
yd <- data.frame(y)
write.csv(yd,'GSEA_reactome_TCGA_high_vs_low.csv',quote = F)
library(GseaVis)
gene <- c('ADM')
geneSetID <- c("HALLMARK_HYPOXIA")
gseaNb(object = y,
       geneSetID = geneSetID,
       addGene = gene,
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.8,0.8))
####nomogram and calibrtion####
library(tidyverse)
clinic <- read.table("TCGA_SARC_phenotype_raw.txt",sep = "\t",header = T,row.names=1)
dat <-read.table("multiCOX3_rickScore_TCGA_DFS_CAFgene5.txt",sep = "\t",header = T,row.names=1)
clinic$sample <- gsub('-','.',clinic$sample)
rownames(clinic) <- clinic$sample
clinic <- clinic[,-1]
comsample <- intersect(rownames(clinic),rownames(dat))
clinic <- clinic[comsample,]
dat <- dat[comsample,]
dat <- dat[,-c(3:8)]
clinic <- clinic[,c('age','histologic.type','sex','resection','TMB')]
identical(rownames(clinic),rownames(dat))
rt <- cbind(dat,clinic)
write.csv(rt,'SARC_TCGA_phenotype_clean.csv',quote = F)
rt2 <- read.csv("SARC_TCGA_phenotype_final.csv",sep = ",",header = T)
library(rms)
library(foreign)
library(survival)
library(survminer)
install.packages("riskRegression",destdir = 'D:/load app/R/downloaded_packages')
library(riskRegression)
rt2$DFS.time <- rt2$DFS.time*30
rt2$risk <- factor(rt2$risk,levels = c('low','high'))
rt2$sex <- factor(rt2$sex,levels = c('Female','Male'))
rt2$histologic.type <- factor(rt2$histologic.type,levels = c('Dedifferentiated liposarcoma',
                                                             'Leiomyosarcoma (LMS)',
                                                             'Myxofibrosarcoma',
                                                             'Pleomorphic MFH / Undifferentiated pleomorphic sarcoma',
                                                             'Synovial Sarcoma',
                                                             'Undifferentiated Pleomorphic Sarcoma (UPS)'))
rt2$margin.status <- factor(rt2$margin.status,levels = c('R1','R0'))
ddist <- datadist(rt2)
options(datadist="ddist")
f <- cph(Surv(DFS.time,DFS)~risk+sex+margin.status+age+TMB,x=T,y=T,surv = T,data = rt2,time.inc = 365)
surv <- Survival(f)
nom2 <- nomogram(f,fun=list(function(x) surv(1,x),function(x) surv(2,x),function(x) surv(3,x)),
                 lp=F,funlabel = c('1-year DFS','2-year DFS','3-year DFS'),
                 maxscale = 100,
                 fun.at = c(0.99,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.0))
plot(nom2)
install.packages('PredictABEL',destdir = 'D:/load app/R/downloaded_packages')
library(PredictABEL)
pred.lg <- predict(f,rt2)
rt2$prob <-1/(1+exp(-pred.lg)) 
prob <- rt2$prob
plotCalibration(data = rt2, cOutcome = 3, predRisk = prob, groups = 10, plottitle = "Calibration plot")  #3为结局指标所在列数
f1 <- cph(Surv(DFS.time,DFS)~risk+sex+margin.status+age+TMB,x=T,y=T,surv = T,data = rt2,time.inc = 365)
cal1 <- calibrate(f1, cmethod = "KM", method = "boot", u = 365, m = 64, B = 500)
plot(cal1, lwd = 2, lty = 2, errbar.col = '#3182BDFF',
     xlab = "Nomogram-Predicted Prognosis", ylab = "Actual Prognosis",
     add=F,
     col = '#3182BDFF', subtitles = FALSE, xlim = c(0,1), ylim = c(0, 1), main = "Calibrate plot")
par(new=T)

f3 <- cph(Surv(DFS.time,DFS)~risk+sex+margin.status+age+TMB,x=T,y=T,surv = T,data = rt2,time.inc = 1095)
cal3 <- calibrate(f3, cmethod = "KM", method = "boot", u = 1095, m = 64, B = 500)
plot(cal3, lwd = 2, lty = 2, errbar.col = '#E6550DFF',
     xlab = "Nomogram-Predicted Prognosis", ylab = "Actual Prognosis",
     add=F,
     col = '#E6550DFF', subtitles = FALSE, xlim = c(0,1), ylim = c(0, 1), main = "Calibrate plot")
par(new=T)

f5 <- cph(Surv(DFS.time,DFS)~risk+sex+margin.status+age+TMB,x=T,y=T,surv = T,data = rt2,time.inc = 1825)
cal5 <- calibrate(f5, cmethod = "KM", method = "boot", u = 1825, m = 64, B = 500)
plot(cal5, lwd = 2, lty = 2, errbar.col = '#31A354FF',
     xlab = "Nomogram-Predicted Prognosis", ylab = "Actual Prognosis",
     add=F,
     col = '#31A354FF', subtitles = FALSE, xlim = c(0,1), ylim = c(0, 1), main = "Calibrate plot")
legend("bottomright", legend=c("1-year DFS", "3-year DFS", "5-year DFS"), col=c('#3182BDFF', '#E6550DFF', '#31A354FF'), lwd=2,lty = 2)
lines(cal[, c("mean.predicted", "KM")], type = "l", lwd = 2, col = c(rgb(192, 98,
                                                                         83, maxColorValue = 255)), pch = 16)
####validation in GSE30929####
#similar as code in TCGA
