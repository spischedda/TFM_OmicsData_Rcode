####################TFM_Epigenomics#####################################
PATH = dirname(rstudioapi::getSourceEditorContext()$path) # RStudio
setwd(PATH)

#Loading packages
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(doParallel)
library(limma) 
library(DMRcate) #find for DMRs with limma annotation
library(RColorBrewer) #color palette
library(ggplot2) #Visualizing Data
library(gplots) #Visualizing Data
library(dplyr) #filter function
library(knitr)
library(missMethyl)
library(Gviz)
library(DMRcate)
library(stringr)


################################# Upload Data ################################
Patients<-read.csv2("Patients_characteristics.csv",sep = ",",header = TRUE)
Patients<-Patients[order(Patients$Lab_ID),]
Patients1<-Patients[c(3,4,11,16,24,29,34,39,40),]

idat.folder <- PATH
targets <- read.metharray.sheet("SampleSheet1.csv",base=idat.folder)
View(targets)
targets$Sample<-c(rep("006",6),rep("008",6),rep("022",6),rep("029",3),rep("049",3),rep("058",6),rep("071",6),rep("079",6),rep("080",6))
targets$Sex<-c(rep("F",6),rep("M",6),rep("M",6),rep("F",3),rep("F",3),rep("F",6),rep("M",6),rep("M",6),rep("M",6))
targets$Age<-c(rep(22,6),rep(24,6),rep(32,6),rep(26,3),rep(22,3),rep(29,6),rep(28,6),rep(42,6),rep(23,6))
targets$Dose<-c(rep(1300,6),rep(1000,6),rep(1000,6),rep(825,3),rep(2400,3),rep(700,6),rep(2750,6),rep(1000,6),rep(1300,6))
targets$smoking<-c(rep("Yes",6),rep("No",6),rep("No",6),rep("No",3),rep("Yes",3),rep("Yes",6),rep("Yes",6),rep("No",6),rep("No",6))


CD14<-subset(targets,Pool_ID=="CD14")
CD56<-subset(targets,Pool_ID=="CD56")

#RGChannel object
#empiezo a trabajar con mis datos: creo un objecto rgset
rgsetCD14 <- read.metharray.exp(targets = CD14)
phenoData <- rgsetCD14$Sample_Group
phenoData
manifest <- getManifest(rgsetCD14)
getProbeInfo(manifest, type = "I")
getProbeInfo(manifest, type = "II")

#MethylSet object
MSetCD14 <- preprocessRaw(rgsetCD14) 

#RatioSet object
RSetCD14<- ratioConvert(MSetCD14, what = "both", keepCN = TRUE)
beta <- getBeta(RSetCD14)

#GenomicRatioSet object
GRsetCD14 <- mapToGenome(RSetCD14)
beta <- getBeta(GRsetCD14)

M <- getM(GRsetCD14)
CN <- getCN(GRsetCD14)

#ANNOTATION 
annotation <- getAnnotation(GRsetCD14)
names(annotation)

#Quality control minfi
qc <- getQC(MSetCD14)
plotQC(qc)

controlStripPlot(rgsetCD14, controls="BISULFITE CONVERSION II")

#save PDF qc reports
qcReport(rgsetCD14,pdf="qcReport_GoodQuality.pdf")

#Quality control ChAMP
library(ChAMP)
QC.GUI(beta=betasCD14,
       pheno=CD14$Sample_Group,
       arraytype="EPIC")

#Sex Prediction
predictedSex <- getSex(GRsetCD14, cutoff = -2)$predictedSex
pred<-as.data.frame(predictedSex)
pred$orig<-CD14$Sex
row.names(pred)<-CD14$Sample_Name
colnames(pred)<-c("PredictedSex","OriginalSex")
View(pred)

pdf("sex.pdf")
plotSex(predictedSex)#annotation issue!!
dev.off()

getSex(object = MSet, cutoff = -2)
addSex(GRset, sex = NULL)
plotSex(GRset, id = NULL)

#Detection P-values
detPCD14<-detectionP(rgsetCD14)
keep = colMeans(detPCD14) < 0.01 #or 0.01
rgsetCD14 = rgsetCD14[,keep]
CD14= CD14[keep,]
detPCD14 = detPCD14[,keep]
dim(detPCD14)

barplot(colMeans(detP), las=2, cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")

#Normalization
gRatioSet.quantileCD14 <- preprocessQuantile(rgsetCD14)##SQN

#Filtering
ann850k<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(ann850k)

detPCD14 <- detPCD14[match(featureNames(gRatioSet.quantileCD14),rownames(detPCD14)),] 
keep <- rowSums(detPCD14 < 0.01) == ncol(gRatioSet.quantileCD14) 
table(keep)
gRatioSet.quantileCD14 <- gRatioSet.quantileCD14[keep,]

keep <- !(featureNames(gRatioSet.quantileCD14) %in% ann850k$Name[ann850k$chr %in% c("chrX","chrY")])
table(keep)

gRatioSet.quantileCD14 <- gRatioSet.quantileCD14[keep,]

gRatioSet.quantileCD14<- dropLociWithSnps(gRatioSet.quantileCD14)

xReactiveProbes <- read.csv(file=paste("13059_2016_1066_MOESM1_ESM.csv",sep="/"),stringsAsFactors=FALSE)
keep <- !(featureNames(gRatioSet.quantileCD14) %in% xReactiveProbes$TargetID)
table(keep)
gRatioSet.quantileCD14 <- gRatioSet.quantileCD14[keep,]

betasCD14<-getBeta(gRatioSet.quantileCD14)
MCD14<- getM(gRatioSet.quantileCD14)

colnames(MCD14)<-CD14$Sample_Name
colnames(betasCD14)<-CD14$Sample_Name

ann850k<-as.data.frame(ann850k)
names(ann850k)

#Age estimation
library(wateRmelon)
AGE<-agep(betasCD14, coeff=NULL, method=c('horvath'))
View(AGE)

#Estimation of Cell Composition
library(FlowSorted.Blood.EPIC)
library(ExperimentHub)
hub <- ExperimentHub()
query(hub, "FlowSorted.Blood.EPIC")
FlowSorted.Blood.EPIC <- hub[["EH1136"]]
data(FlowSorted.Blood.EPIC)
cells_CD14<-estimateCellCounts2(rgsetCD14,
                                compositeCellType = "Blood",
                                processMethod = "preprocessQuantile",
                                probeSelect = "any",
                                cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
                                referencePlatform = "IlluminaHumanMethylationEPIC",
                                referenceset = NULL,
                                IDOLOptimizedCpGs = NULL,
                                returnAll = TRUE,
                                meanPlot = TRUE,
                                verbose = TRUE)

cell_counts_CD14 <- as.data.frame(cells_CD14$prop)

cd8t=as.numeric(cell_counts_CD14$CD8T) #Set variables
cd4t=as.numeric(cell_counts_CD14$CD4T)
nk=as.numeric(cell_counts_CD14$NK)
bcell=as.numeric(cell_counts_CD14$Bcell)
mono=as.numeric(cell_counts_CD14$Mono)
neu=as.numeric(cell_counts_CD14$Neu)

shapiro.test(mono) #repeat for each cell type
cell_counts_CD14$pheno<-CD14$Sample_Group
D0<-subset(cell_counts_CD14,Sample_Group=="D0")
D28<-subset(cell_counts_CD14,Sample_Group=="D28")
ED<-subset(cell_counts_CD14,Sample_Group=="ED")
wilcox.test(D0$mono,ED$mono) #repeat for each cell type

#boxplot cell composition
cell_counts_CD14$Names<-CD56$Sample_Name

library(reshape2)
cells<- melt(cell_counts_CD14, id.vars=c("Names","Phenotype"))
library(ggrepel)
library(PupillometryR)
library(ggsignif)
p2 <- ggplot(cells, aes(x=factor(Phenotype,levels = c("D0","ED","D28")),y=value,show.legend = FALSE))+
  facet_wrap(~variable,nrow = 1, ncol =6,labeller=label_context)+
  geom_point(aes(color = Phenotype), 
             position = position_jitter(w = .15), 
             size = 1.5,
             show.legend = F) +
  geom_boxplot(aes(color = Phenotype),
               outlier.shape = NA) +
  #geom_text_repel(aes(x=factor(Phenotype), y=value, label = Names))+
  geom_signif(comparisons = list(c("ED","D0"),c("ED","D28")),y_position = 0.95,step_increase = 0.1,test="wilcox.test",map_signif_level = TRUE,textsize = 6)+ #map_signif_level=c("<0.001"=0.001, "<0.01"=0.01, "<0.05"=0.05)
  scale_fill_manual(values = c("mediumseagreen","lightslateblue","#D55E00"))+ #56B4E9","mediumorchid"
  scale_color_manual(values = c("mediumseagreen","lightslateblue","#D55E00"))  +
  #scale_y_continuous(labels = axis_unit_scaler, breaks = breaks_width(25e3)) +
  labs(x = "",
       y = "Cell proportion",
       fill = "Phenotype: ",
       title = "") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.margin = margin(-5, 5, 5, 5),
        legend.text = element_text(size=25),
        legend.title= element_text(size=20),
        plot.title = element_text(hjust = 0.5,size=20),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=20,color="black"),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 20, b = 0, l = 0)),
        strip.text=element_text(color="black",size = 20))
tiff("cell_estimation_CD14.tiff",units="in", width=20, height=8, res=300)
p2
dev.off()

#check for smoking factor:
D0<-subset(CD14,Sample_Group=="D0")
MD0<-MCD14[,grep("D0",colnames(MCD14))]
group<-factor(D0$smoking)
sex<-factor(D0$Sex)
age<-as.numeric(D0$Age)

design1<-model.matrix(~0+group+sex+age)
fit_D0 <- lmFit(MD0, design)

contMatrix <- makeContrasts(groupyes-groupno,
                            levels=design1)

# fit the contrast
fit1 <- contrasts.fit(fit_D0, contMatrix)
fit2 <- eBayes(fit1)

summary(decideTests(fit2,p.value = 0.05)) 

DMPsD0<- topTable(fit2, num=Inf, coef=NULL, genelist=ann850kSub)

#Manhattan plot
Table<-DMPsD0[,c(1,2,4,52)]
Table$chr<-gsub("chr*","",Table$chr)
Table$chr<-as.numeric(Table$chr)
Table$pos<-as.numeric(Table$pos)
Table$P.Value<-as.numeric(Table$P.Value)

library(qqman)
manhattan(Table, chr = "chr", bp = "pos",  p = "P.Value",  snp = "Name")

#Search for Diferentially Mehtylated Positions (DMPs)
#limma
group<-factor(CD14$Sample_Group)
individ<-CD14$Sample
sex<-factor(CD14$Sex)
age<-as.numeric(CD14$Age)

design<-model.matrix(~0+group+sex+age+cd4t+bcell+nk+neu+cd8t+mono)
dupfit_CD14 <- duplicateCorrelation(MCD14, design, ndups=1, block=individ)
fit_CD14 <- lmFit(MCD14, design, correlation=dupfit_CD14$consensus.correlation, block=individ)

contMatrix <- makeContrasts(groupD0-groupED,
                            groupD28-groupED,
                            groupD0-groupD28,
                            levels=design)

# fit the contrast
fit1 <- contrasts.fit(fit_CD14, contMatrix)
fit2 <- eBayes(fit1)

summary(decideTests(fit2,p.value = 0.05)) 

top.table <- topTable(fit2, coef = NULL, number = Inf, adjust = "fdr")
head(top.table)
hist(top.table$P.Value, breaks = 100, main = "results P")

#QQ-plot
#P.value inflation

library(gaston)
library(RcppParallel)
qqplot.pvalues(top.table$P.Value, col.abline = "red", CB = TRUE, col.CB = "gray80", 
               CB.level = 0.95, thinning = TRUE) 


#Results 
ann850kSub <- ann850k[match(rownames(MCD14),ann850k$Name),
                      c(1:25,26:ncol(ann850k))]

DMPsCD14<- topTable(fit2, num=Inf, coef=NULL, genelist=ann850kSub) #change coefficinet to see the individual contrast results
row.names(DMPsCD14)<-DMPsCD14$Name

PROMSCD14<-merge(DMPsCD14,betasCD14,by="row.names")
rownames(PROMSCD14)<-PROMSCD14$Row.names
PROMSCD14<-PROMSCD14[,-1]

PROMS_001CD14<-subset(PROMSCD14,PROMSCD14$P.Value <0.001)
PROMS_05CD14<-subset(PROMSCD14,PROMSCD14$adj.P.Val<0.05)

write.table(PROMS_05, file="PROMS_05_FDR_CD14.txt", sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

#Volcano plot
Tabla<-PROMSCD14_2[,c(51,81)]

plot(Tabla$DeltaD0_ED,-log10(Tabla$adj.P.Val))
cols <- densCols(Tabla$DeltaD0_ED,-log10(Tabla$adj.P.Val))
plot(Tabla$DeltaD0_ED,-log10(Tabla$adj.P.Val), col=cols, panel.first=grid(),
     main="Volcano plot Hypo/Hypermethylated CpGs", xlab="Delta value", ylab="log10(pvalue)",
     pch=20, cex=0.6)
abline(v=c(-0.1,0.1), col="brown")
abline(h=-log10(alpha), col="brown")

#Data visualization PCA 
l1<-PROMS_05CD14[,c(54:79)]
Data<-t(l1) 

pheno<-CD14$Sample_Group

pca<-prcomp(Data)
summary(pca)

library(factoextra)
library(FactoMineR)
res.pca <- PCA(Data,scale.unit = TRUE, ncp = 5, graph = TRUE)
summary(res.pca)
ind <- get_pca_var(res.pca)

ind.p<-fviz_pca_ind(res.pca,
                    axes.linetype = "blank",
                    geom.ind = c("point"), # show points only (nbut not "text")
                    pointshape = 21,
                    mean.point=FALSE,
                    pointsize = 6,
                    fill.ind = pheno,
                    col.ind = pheno, # color by groups
                    palette=c("mediumseagreen","lightslateblue","#D55E00"),                                        #addEllipses = TRUE, # Concentration ellipses
                    #gradient.cols = c("#E7B800","#e34a33"),
                    legend.title = "Phenotype",labelsize = 8)+
  labs(title=" ")+
  theme(text = element_text(size=20),
        #axis.line=element_blank(),
        axis.text=element_text(color = "black",size=20),
        axis.text.y=element_text(angle = 90, vjust = 0.5, hjust=1,color = "black",margin = margin(15, 0, 0, 0)),
        axis.title = element_text(size=20),
        axis.title.x = element_text(margin = margin(0, 15, 0, 0),color = "black",face = "plain"),
        axis.title.y = element_text(margin = margin(0, 15, 0, 0),color = "black",face = "plain"),
        #legend.position="bottom",
        legend.title=element_text(size=20), 
        panel.background=element_blank(),
        panel.border=element_rect(fill = NA),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background=element_rect(fill="white"),
        panel.spacing = unit(4, "lines"))
tiff("PCA_CD14.tiff",units="in", width=10, height=7, res=300)
ind.p
dev.off()

#PLOT DMPs
PROMS_05CD14<-PROMS_05CD14[c(2,3,6,7,8,9,10,12,13,15),]

CPG_FDR<-PROMS_05CD14[grep("cg05304729|cg053047291|cg17661135|cg15854154|cg13575798|cg04218845|cg25734220|cg11941024",PROMSCD14$Name),]

CPG_FDR<-CPG_FDR[,c(54:79)]
CPG_FDR<-t(CPG_FDR)
CPG_FDR<-as.data.frame(CPG_FDR)
CPG_FDR$Phenotype<-CD14$Sample_Group
CPG_FDR$Names<-row.names(CPG_FDR)

library(reshape2)
Data2<- melt(CPG_FDR, id.vars=c("Phenotype","Names"))
library(PupillometryR)
library(ggsignif)
AA<-ggplot(Data2) +
  aes(x = factor(Phenotype,levels=c("D0","ED","D28"),),y=value,label=Names) +
  geom_flat_violin(aes(color = Phenotype,fill=Phenotype),show.legend = F,position = position_nudge(x = .4), 
                   alpha = .4) +
  facet_wrap(~variable,nrow = 1, ncol =11,labeller=label_context)+
  geom_boxplot(aes(color = Phenotype),width = .45, 
               outlier.shape = NA,
               alpha = 0.5) +
  scale_fill_manual(values = c("mediumseagreen","lightslateblue","#D55E00")) +
  scale_color_manual(values = c("mediumseagreen","lightslateblue","#D55E00"))  +
  labs(x = "",
       y = "Beta Values",
       fill = "Phenotype: ",
       title = "") +
  guides(fill = guide_legend(nrow=1,
                             byrow=TRUE))+
  theme_bw() +
  theme(legend.position = "bottom",
        legend.margin = margin(-5, 5, 5, 5),
        legend.text = element_text(size=22),
        legend.title= element_text(size=22),
        plot.title = element_text(hjust = 0.5,size=22),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=22,color="black"),
        axis.title.y = element_text(size=22,margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(size=22),
        strip.text=element_text(color="black",size = 18.5))
tiff("CD14_Top_CpGannotated_FDR0.05_BoxPlot.tiff",units="in", width=22, height=9, res=300)
AA
dev.off()

# ROC CURVES
library(ROCR)
library(pROC)

roc1 <- roc(D0_ED$Phenotype,D0_ED$cg04218845)
roc2 <- roc(D0_ED$Phenotype,D0_ED$cg05304729)
roc3 <- roc(D0_ED$Phenotype,D0_ED$cg09677297)
roc4 <- roc(D0_ED$Phenotype,D0_ED$cg11941024)
roc5 <- roc(D0_ED$Phenotype,D0_ED$cg13575798)
roc6 <- roc(D0_ED$Phenotype,D0_ED$cg15854154)
roc7 <- roc(D0_ED$Phenotype,D0_ED$cg17661135)
roc8 <- roc(D0_ED$Phenotype,D0_ED$cg17711596)
roc9 <- roc(D0_ED$Phenotype,D0_ED$cg25734220)


Roc1<-ggroc(roc1, legacy.axes = TRUE,size=1)+ 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="red",linetype="dashed",size=1)+
  labs(title="",x="False Positive Rate",y="True Positive Rate")+ 
  geom_text(size=6,aes(0.5, 1,
                       label = paste("AUC =",sprintf("%.3f",roc1$auc)),hjust = 0,vjust=15))+ 
  theme(text = element_text(size=20),
        #axis.line=element_blank(),
        axis.text.x=element_text(color = "black",face = "plain"),
        axis.text.y=element_text(color = "black",face = "plain"),
        axis.title.x = element_text(color = "black",face = "plain",margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.y = element_text(color = "black",face = "plain",margin = margin(t = 20, r = 20, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5),
        #legend.position="bottom",
        #legend.title=element_text(hjust = 0.5), 
        panel.background=element_blank(),
        panel.border=element_rect(fill = NA),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background=element_rect(fill="white"),
        panel.spacing = unit(1, "lines"))
tiff("ROC1_D28-ED.tiff", units="in", width=5, height=5, res=300)
Roc1
dev.off()  

#methylGSA (methylglm)
library(methylGSA)
cpg<-DMPsCD14$P.Value
names(cpg) <- DMPsCD14$Name
head(cpg,20)

resKEGG = methylglm(cpg.pval = cpg, minsize = 50, array.type = "EPIC",
                    maxsize = 500, GS.type = "KEGG", group="all")
View(resKEGG)

res1GO_all = methylglm(cpg.pval = cpg, minsize = 20,group="all",array.type = "EPIC", 
                       maxsize = 500, GS.type = "GO")
View(res1GO_all)

write.table(res1GO_all,file="EnrichPath_all_DMPs_GO_last.txt")

res1Reac = methylglm(cpg.pval = cpg, minsize = 50,group="all",array.type = "EPIC", 
                     maxsize = 500, GS.type = "Reactome")
View(res1Reac)

SignReac<-res1Reac[c(1:20),] #repeat the same with GO and KEGG
SignReac$`FDR P-value`<-SignReac$padj
View(SignReac)
SignReac$Description<-gsub("Homo sapiens: *","",SignReac$Description)

library(forcats) ## for reordering the factor
p1<-ggplot(SignGO, aes(x = Size, y = fct_reorder(Description, Size),fill =`FDR P-value`)) + 
  geom_bar(aes(x = Size),stat="identity") +
  theme_bw(base_size = 10) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  scale_fill_gradient(limits=c(0,0.07),high="blue", low = "red",guide=guide_colourbar(reverse = FALSE)) +
  ylab(NULL) +labs(color="FDR P.value")+
  ggtitle("")
tiff("GO_CD14.tiff",units="in", width=6, height=1.5, res=300)
p1
dev.off()


