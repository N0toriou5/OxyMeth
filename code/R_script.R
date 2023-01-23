### R script to analyze and compare transcriptomics data signatures

homedir <- "D:/Projects/OxyMeth"
#load required environment
#dir.create("R_Analysis/data")
dir.create("results")
dir.create("plots")
options(java.parameters = "-Xmx8G")
require("biomaRt")
require("corto")
require("DESeq2")
require("EnhancedVolcano")
require("fgsea")
require("msigdbr")
require("Rsubread")
require("stringr")
require("vulcan")
require("xlsx")
source("D:/Archive/textplot3.R")

#### count reads
bam.files <- list.files(path = "data/PRJEB33332/trimmed/hisat_bams", pattern = ".bam$", full.names = TRUE)
#props <- propmapped(files=bam.files)
annot <- "D:/genomes/Rat/Rnor6_HIV/Rnor6_hybrid_fixed.gtf"
fc <- featureCounts(bam.files, annot.ext = annot,
                    isGTFAnnotationFile = TRUE,GTF.featureType = "gene",
                    GTF.attrType = "gene_name",strandSpecific = 2,isPairedEnd = TRUE, nthreads = 14)
rawcounts <- fc$counts
# Remove sorted.bam from filenames
colnames(rawcounts) <- gsub("\\.sorted.bam", "", colnames(rawcounts))
colnames(rawcounts) <- gsub("Sample_", "", colnames(rawcounts))
#colnames(rawcounts) <- gsub(".+\\.[0-9]","",colnames(rawcounts)) 
save(rawcounts,file = "results/000_rawcounts_input.rda")

#### Divide the dataset into HIV e WT rats

#hivcounts<-rawcounts[str_detect(colnames(rawcounts), "HIV")]
#wtcounts<-rawcounts[str_detect(colnames(rawcounts), "WT")]

conditions<-c("Type","Treatment")
annotation<-matrix(nrow=ncol(rawcounts),ncol=length(conditions),dimnames = list(colnames(rawcounts),conditions))
annotation[,"Type"] <- ifelse(str_detect(colnames(rawcounts), "H")==TRUE,"HIV","WT")
annotation[,"Treatment"]<-ifelse(str_detect(colnames(rawcounts), "m|M")==TRUE, "Meth", "Naive")
annotation<-as.data.frame(annotation,stringsAsFactors=FALSE)
save(annotation,file = "results/000_annotation_input.rda")  

# Variance Stabilizing Transformation
expmat <- vst(rawcounts, blind = TRUE, nsub = 1000, fitType = "parametric")
save(expmat,file="results/000_expmat.rda")

# PCA

png("plots/000_subreadpca_input.png",w=2000,h=2000,p=25)
pca<-prcomp(t(expmat))
totvar<-sum(pca$sdev^2)
pcavar<-((pca$sdev^2)/totvar)*100
x<-setNames(pca$x[,1],colnames(expmat))
y<-setNames(pca$x[,2],colnames(expmat))
plot(x,y,pch=20, main=paste0("HIV vs. WT Meth"),
     xlim=c(min(x)*1.5,max(x)*1.5),cex=2,col="grey",
     xlab=paste0("PC",1," (",signif(pcavar[1],3),"%)"),
     ylab=paste0("PC",2," (",signif(pcavar[2],3),"%)")
)
#text(x,y,labels=colnames(expmat))
textplot2(x,y,colnames(expmat),new=FALSE)
grid()
dev.off()

#### compare just HIVm vs. HIV and WTm vs. WT
fname<-"results/001_deres.rda"
if(!file.exists(fname)){
  deres<-list()
  ### Dependent
  for(type in c("HIV","WT")){
    subsamples<-rownames(annotation)[annotation$Type==type]
    subraw<-rawcounts[,subsamples]
    subannot<-annotation[subsamples,]
    # DESeq2 block (filter out poorly expressed genes)
    dds <- DESeqDataSetFromMatrix(countData=subraw,colData=subannot,design=~Treatment)
    dds <- dds[rowSums(counts(dds)) >=10,]
    dds$Treatment<-relevel(dds$Treatment, ref = "Naive")
    dea <- DESeq(dds, parallel = TRUE)
    #res <- results(dea, contrast = c("Treatment","Meth","Naive"))
    #res <- as.data.frame(lfcShrink(dds = dea, contrast=c("Treatment","Meth","Naive"), res = res, type = "ashr")) 
    res <- as.data.frame(results(dea,contrast=c("Treatment","Meth","Naive")))
    deres[[type]][["Meth"]]<-res
    #res_A <- lfcShrink(dds=dea, coef=2, res=res)
  }
  save(deres,file=fname)
} else {load(fname)}

#### HIV plots
#Meth
res <- deres$HIV$Meth
res <- na.omit(res)
res <- res[order(res$padj),]
topup <- rownames(res[res$log2FoldChange>=0.5,])[1:10]
topdn <- rownames(res[res$log2FoldChange<= -0.5,])[1:10]
top <- c(topup,topdn)
png("plots/000_volcano_HIVm_vs_HIV_Meth.png", w=2500,h=2500, res=300)
EnhancedVolcano(res, subtitle = "",
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = top,
                xlim = c(-2, 2),
                #ylim = c(0,12),
                title = 'HIVm vs. HIV (Self-administration - Acute)',
                pCutoff = 0.05, #0.05 cutoff
                FCcutoff = 0.5, # 2-fold change
                labFace = "bold",
                labSize = 4,
                col = c('lightgrey', 'lightblue', 'lightpink', 'salmon'),
                colAlpha = 4/5,
                drawConnectors = TRUE,
                widthConnectors = 0.3,colConnectors = 'gray51',maxoverlapsConnectors = Inf,
                legendLabSize = 14,
                legendIconSize = 4.0,
                caption = paste0('Upregulated = ', nrow(res[res$log2FoldChange>0.5&res$padj<=0.05,]), ' genes',"\n",'Downregulated = ',
                                 nrow(res[res$log2FoldChange< -0.5&res$padj<=0.05,]), ' genes'))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
write.xlsx(res,file = "results/001_DEA_HIVm_HIV_Meth.xlsx")


# New corto double GSEA
msigdbr_species()
mdf <- msigdbr(species="Rattus norvegicus", category="C2") # Retrieve all RAT gene sets
length(mdf$gs_name) # Number of associations 504876
length(length(unique(mdf$gs_name))) # Number of pathways
head(mdf)
mlist <- mdf %>% split(x=.$gene_symbol,f=.$gs_name)
mlist <- mlist[grep("(REACTOME|WP|PATHWAY|KEGG)",names(mlist),value="TRUE")]
#grep("(STAT1|STAT2|STAT3)",names(mlist),value="TRUE")
plength <- sapply(mlist,length)
max(plength) # 1975, the biggest pathway
# GSEA
signature <- setNames(res$stat, rownames(res))
set.seed(1)
gseas <- fgseaMultilevel(pathways = mlist, stats = signature, eps = 0, minSize = 15,
                       maxSize = Inf, nproc = 14, nPermSimple = 10000)
gseas <- gseas[order(gseas$pval),]
save(gseas, file="results/001_GSEA_HIVm_HIV_Meth.rda")
write.xlsx2(gseas, file= "results/001_GSEA_HIVm_HIV_Meth.xlsx", row.names = FALSE)

### WT
#Meth
res<-deres$WT$Meth
res<-na.omit(res)
res<-res[order(res$padj),]
# topup<-rownames(res[res$log2FoldChange>=0.5,])[1:10]
# topdn<-rownames(res[res$log2FoldChange<= -0.5,])[1:10]
# top<-c(topup,topdn)
png("plots/001_volcano_WTm_vs_WT_Meth.png", w=2500,h=2500, res=300)
EnhancedVolcano(res, subtitle = "",
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = top,
                xlim = c(-2, 2),
                #ylim = c(0,12),
                title = 'WTm vs. WT (Self-administration - Acute)',
                pCutoff = 0.05, #0.05 cutoff
                FCcutoff = 0.5, # 2-fold change
                labFace = "bold",
                labSize = 4,
                col = c('lightgrey', 'lightblue', 'lightpink', 'salmon'),
                colAlpha = 4/5,
                drawConnectors = TRUE,
                widthConnectors = 0.3,colConnectors = 'gray51',maxoverlapsConnectors = Inf,
                legendLabSize = 14,
                legendIconSize = 4.0,
                caption = paste0('Upregulated = ', nrow(res[res$log2FoldChange>0.5&res$padj<=0.05,]), ' genes',"\n",'Downregulated = ',
                                 nrow(res[res$log2FoldChange< -0.5&res$padj<=0.05,]), ' genes'))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
write.xlsx(res,file = "results/001_DEA_WTm_WT_Meth.xlsx")

# GSEA
signature <- setNames(res$stat, rownames(res))
set.seed(1)
gseas <- fgseaMultilevel(pathways = mlist, stats = signature, eps = 0, minSize = 15,
                         maxSize = Inf, nproc = 14, nPermSimple = 10000)
gseas <- gseas[order(gseas$pval),]
save(gseas, file="results/001_GSEA_WTm_WT_Meth.rda")
write.xlsx2(gseas, file= "results/001_GSEA_WTm_WT_Meth.xlsx", row.names = FALSE)

###### Dataset P35458399
#### count reads
bam.files <- list.files(path = "data/P35458399/trimmed/hisat_bams/", pattern = ".bam$", full.names = TRUE)
#props <- propmapped(files=bam.files)
annot <- "D:/genomes/Rat/Rnor6_HIV/Rnor6_hybrid_fixed.gtf"
fc <- featureCounts(bam.files, annot.ext = annot,
                    isGTFAnnotationFile = TRUE,GTF.featureType = "gene",
                    GTF.attrType = "gene_name",strandSpecific = 2,isPairedEnd = TRUE, nthreads = 14)

#%succesfully assigned reads: 66.5, 64.4, 71.2, 68.4, 63.9, 66.7, 66.9, 61.8, 65.8, 65.9, 69.0, 69.4, 70.0, 70.3, 65.3, 71.4, 67.0, 68.2, 69.0, 67.5, 67.5, 68.8, 
# 62.6, 70.3, 68.6, 69.0, 70.2, 70.5, 66.4, 69.5, 68.8    mean = 67.8
rawcounts <- fc$counts

# Remove sorted.bam from filenames
colnames(rawcounts) <- gsub("\\.sorted.bam", "", colnames(rawcounts))
colnames(rawcounts) <- gsub("_L[0-9]*","",colnames(rawcounts)) 
colnames(rawcounts) <- gsub("_S[0-9]*_*", "",colnames(rawcounts))

#colnames(rawcounts) <- gsub(".+\\.[0-9]","",colnames(rawcounts)) 
save(rawcounts,file = "results/001_rawcounts_Oxy.rda")

conditions<-c("Type","Treatment")
annotation<-matrix(nrow=ncol(rawcounts),ncol=length(conditions),dimnames = list(colnames(rawcounts),conditions))
annotation[,"Type"] <- ifelse(str_detect(colnames(rawcounts), "H")==TRUE,"HIV","WT")
annotation[,"Treatment"]<-ifelse(str_detect(colnames(rawcounts), "D")==TRUE, "Meth", "Naive")
annotation<-as.data.frame(annotation,stringsAsFactors=FALSE)
save(annotation,file = "results/001_annotation_Oxy.rda")  

# Variance Stabilizing Transformation
expmat <- vst(rawcounts, blind = TRUE, nsub = 1000, fitType = "parametric")
save(expmat,file="results/001_expmat_Oxy.rda")

# PCA

png("plots/001_subreadpca_Oxy.png",w=2000,h=2000,p=25)
pca<-prcomp(t(expmat))
totvar<-sum(pca$sdev^2)
pcavar<-((pca$sdev^2)/totvar)*100
x<-setNames(pca$x[,1],colnames(expmat))
y<-setNames(pca$x[,2],colnames(expmat))
plot(x,y,pch=20, main=paste0("HIV vs. WT Meth"),
     xlim=c(min(x)*1.5,max(x)*1.5),cex=2,col="grey",
     xlab=paste0("PC",1," (",signif(pcavar[1],3),"%)"),
     ylab=paste0("PC",2," (",signif(pcavar[2],3),"%)")
)
#text(x,y,labels=colnames(expmat))
textplot2(x,y,colnames(expmat),new=FALSE)
grid()
dev.off()

#### compare just HIVm vs. HIV and WTm vs. WT
fname<-"results/002_deres_Oxy.rda"
if(!file.exists(fname)){
  deres<-list()
  ### Dependent
  for(type in c("HIV","WT")){
    subsamples<-rownames(annotation)[annotation$Type==type]
    subraw<-rawcounts[,subsamples]
    subannot<-annotation[subsamples,]
    # DESeq2 block (filter out poorly expressed genes)
    dds <- DESeqDataSetFromMatrix(countData=subraw,colData=subannot,design=~Treatment)
    dds <- dds[rowSums(counts(dds)) >=10,]
    dds$Treatment<-relevel(dds$Treatment, ref = "Naive")
    dea <- DESeq(dds, parallel = TRUE)
    #res <- results(dea, contrast = c("Treatment","Meth","Naive"))
    #res <- as.data.frame(lfcShrink(dds = dea, contrast=c("Treatment","Meth","Naive"), res = res, type = "ashr")) 
    res <- as.data.frame(results(dea,contrast=c("Treatment","Meth","Naive")))
    deres[[type]][["Meth"]]<-res
    #res_A <- lfcShrink(dds=dea, coef=2, res=res)
  }
  save(deres,file=fname)
} else {load(fname)}

#### HIV plots
#Meth
res<-deres$HIV$Meth
res<-na.omit(res)
res<-res[order(res$padj),]
#topup<-rownames(res[res$log2FoldChange>=0.5,])[1:10]
#topdn<-rownames(res[res$log2FoldChange<= -0.5,])[1:10]
#top<-c(topup,topdn)
png("plots/002_volcano_HIVm_vs_HIV_Oxy.png", w=2500,h=2500, res=300)
EnhancedVolcano(res, subtitle = "",
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = top,
                xlim = c(-2, 2),
                #ylim = c(0,12),
                title = 'HIV + Oxy vs. HIV (Self-administration - Acute)',
                pCutoff = 0.05, #0.05 cutoff
                FCcutoff = 0.5, # 2-fold change
                labFace = "bold",
                labSize = 4,
                col = c('lightgrey', 'lightblue', 'lightpink', 'salmon'),
                colAlpha = 4/5,
                drawConnectors = TRUE,
                widthConnectors = 0.3,colConnectors = 'gray51',maxoverlapsConnectors = Inf,
                legendLabSize = 14,
                legendIconSize = 4.0,
                caption = paste0('Upregulated = ', nrow(res[res$log2FoldChange>0.5&res$padj<=0.05,]), ' genes',"\n",'Downregulated = ',
                                 nrow(res[res$log2FoldChange< -0.5&res$padj<=0.05,]), ' genes'))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
write.xlsx(res, file = "results/001_DEA_HIVm_HIV_Oxy.xlsx")


# New corto double GSEA
msigdbr_species()
mdf<-msigdbr(species="Rattus norvegicus", category="C2") # Retrieve all RAT gene sets
length(mdf$gs_name) # Number of associations 504876
length(length(unique(mdf$gs_name))) # Number of pathways
head(mdf)
mlist<-mdf %>% split(x=.$gene_symbol,f=.$gs_name)
mlist<-mlist[grep("(REACTOME|WP|PATHWAY|KEGG)",names(mlist),value="TRUE")]
#grep("(STAT1|STAT2|STAT3)",names(mlist),value="TRUE")
plength<-sapply(mlist,length)
max(plength) # 1975, the biggest pathway

#GSEA
signature <- setNames(res$stat, rownames(res))
set.seed(1)
gseas<-fgseaMultilevel(pathways = mlist, stats = signature, eps = 0, minSize = 15,
                       maxSize = Inf, nproc = 14, nPermSimple = 10000)
gseas<-gseas[order(gseas$pval),]
save(gseas, file = "results/001_GSEA_HIVm_HIV_Oxy.rda")
write.xlsx2(as.data.frame(gseas), file = "results/001_GSEA_HIVm_HIV_Oxy.xlsx", row.names = FALSE)

### WT
#Meth
res<-deres$WT$Meth
res<-na.omit(res)
res<-res[order(res$padj),]
# topup<-rownames(res[res$log2FoldChange>=0.5,])[1:10]
# topdn<-rownames(res[res$log2FoldChange<= -0.5,])[1:10]
# top<-c(topup,topdn)
png("plots/001_volcano_WTm_vs_WT_Oxy.png", w=2500,h=2500, res=300)
EnhancedVolcano(res, subtitle = "",
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab = top,
                xlim = c(-2, 2),
                #ylim = c(0,12),
                title = 'WT + Oxy vs. WT (Self-administration - Acute)',
                pCutoff = 0.05, #0.05 cutoff
                FCcutoff = 0.5, # 2-fold change
                labFace = "bold",
                labSize = 4,
                col = c('lightgrey', 'lightblue', 'lightpink', 'salmon'),
                colAlpha = 4/5,
                drawConnectors = TRUE,
                widthConnectors = 0.3,colConnectors = 'gray51',maxoverlapsConnectors = Inf,
                legendLabSize = 14,
                legendIconSize = 4.0,
                caption = paste0('Upregulated = ', nrow(res[res$log2FoldChange>0.5&res$padj<=0.05,]), ' genes',"\n",'Downregulated = ',
                                 nrow(res[res$log2FoldChange< -0.5&res$padj<=0.05,]), ' genes'))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
write.xlsx(res,file = "results/001_DEA_WTm_Oxy.xlsx")

# New corto double GSEA
msigdbr_species()
mdf <- msigdbr(species="Rattus norvegicus", category="C2") # Retrieve all RAT gene sets
length(mdf$gs_name) # Number of associations 504876
length(length(unique(mdf$gs_name))) # Number of pathways
head(mdf)
mlist <- mdf %>% split(x=.$gene_symbol,f=.$gs_name)
mlist <- mlist[grep("(REACTOME|WP|PATHWAY|KEGG)",names(mlist),value="TRUE")]
#grep("(STAT1|STAT2|STAT3)",names(mlist),value="TRUE")
plength <- sapply(mlist,length)
max(plength) # 1975, the biggest pathway
# GSEA
signature <- setNames(res$stat, rownames(res))
set.seed(1)
gseas <- fgseaMultilevel(pathways = mlist, stats = signature, eps = 0, minSize = 15,
                         maxSize = Inf, nproc = 14, nPermSimple = 10000)
gseas <- gseas[order(gseas$pval),]
save(gseas, file="results/001_GSEA_HIVm_HIV_Oxy.rda")
write.xlsx2(gseas, file= "results/001_GSEA_HIVm_HIV_Oxy.xlsx", row.names = FALSE)

#### Compare DE
# Meth vs. Oxy in HIV
load("results/001_deres.rda")
meth_DE <- deres
load("results/002_deres.rda")
oxy_DE <- deres

### scatterplot
hiv_meth <- na.omit(meth_DE$HIV$Meth)
hiv_oxy <- na.omit(oxy_DE$HIV$Meth)

x <- setNames(hiv_meth$stat, rownames(hiv_meth))
y <- setNames(hiv_oxy$stat, rownames(hiv_oxy))

common <- intersect(names(x), names(y))
x <- x[common]
y <- y[common]
sig1 <- rownames(hiv_meth[abs(hiv_meth$log2FoldChange) >= 0.5 & hiv_meth$padj <= 0.05, ])
sig2 <- rownames(hiv_oxy[abs(hiv_oxy$log2FoldChange) >= 0.5 & hiv_oxy$padj <= 0.05, ])
length(intersect(sig1,sig2))
png("plots/000_DE_HIV_compare.png",w=2000,h=2000,res=300)
plot(x,y,pch=20,col="grey",xlab="Meth vs. Naive (stat)",ylab="Oxy vs. Naive (stat)",main="HIV")
pcc<-cor.test(x,y)
mtext(paste0("R=",signif(pcc$estimate,2)," p=",signif(pcc$p.value,3)))
lml<-lm(y~x)
abline(lml,lwd=1)
set.seed(1)
points(x[sig1],y[sig1],col="red",pch=20)
points(x[sig2],y[sig2],col="blue",pch=20)
#textplot3(x[top],y[top],words=top,font=2,cex=1,show.lines=F,col="black")
legend("bottomright",pch=20,legend=c(paste0("Significant in Meth: ",length(sig1)),
                                     paste0("Significant in Oxy: ",length(sig2))),col=c("red","blue"),pt.cex=2)
dev.off()

# Meth vs. Oxy in WT
load("results/001_deres.rda")
meth_DE <- deres
load("results/002_deres.rda")
oxy_DE <- deres

### scatterplot
wt_meth <- na.omit(meth_DE$WT$Meth)
wt_oxy <- na.omit(oxy_DE$WT$Meth)

x <- setNames(wt_meth$stat, rownames(wt_meth))
y <- setNames(wt_oxy$stat, rownames(wt_oxy))

common <- intersect(names(x), names(y))
x <- x[common]
y <- y[common]
sig1 <- rownames(wt_meth[abs(wt_meth$log2FoldChange) >= 0.5 & wt_meth$padj <= 0.05, ])
sig2 <- rownames(wt_oxy[abs(wt_oxy$log2FoldChange) >= 0.5 & wt_oxy$padj <= 0.05, ])
length(intersect(sig1,sig2))
png("plots/000_DE_WT_compare.png",w=2000,h=2000,res=300)
plot(x,y,pch=20,col="grey",xlab="Meth vs. Naive (stat)",ylab="Oxy vs. Naive (stat)",main="WT")
pcc<-cor.test(x,y)
mtext(paste0("R=",signif(pcc$estimate,2)," p=",signif(pcc$p.value,3)))
lml<-lm(y~x)
abline(lml,lwd=1)
set.seed(1)
points(x[sig1],y[sig1],col="red",pch=20)
points(x[sig2],y[sig2],col="blue",pch=20)
#textplot3(x[top],y[top],words=top,font=2,cex=1,show.lines=F,col="black")
legend("bottomright",pch=20,legend=c(paste0("Significant in Meth: ",length(sig1)),
                                     paste0("Significant in Oxy: ",length(sig2))),col=c("red","blue"),pt.cex=2)
dev.off()
