

######DIFFERENTIAL GENE EXPRESSION ANALYSIS simple

#using DESeq2 package. Can also use limma/voom or edgeR
#For DESeq2, start with two objects, counts matrix, and column data table (see vignette)

#load package
BiocManager::install("DESeq2")
library(DESeq2)

#construct DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = coldata,
                              design = ~Group)
#summarize dds object stats
dds

#set factor level, specify 'untreated' or it will choose the first group
dds$Group <- relevel(dds$Group, ref = "Non-Neurological Control")

#run DESeq function- may take a while depending on dataset size
dds <- DESeq(dds)

#assign new results to object
res <- results(dds)
res

#summarize results
summary(res)

#make new results object with different FDR cutoff (alpha)
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

#how many genes had an padj less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
#35
#even more strict
sum(res$padj < 0.05, na.rm=TRUE)
#13

#extracting transformed values
#log transforming data
vsd <- vst(dds, blind=F)
rld <- rlog(dds, blind=F)

#set contrasts
resultsNames(dds)

#specify group category name, untreated, then treated samples to compare
res2<-results(dds, contrast = c("Group","Non-Neurological Control", "ALS Spectrum MND"))
vsd2 <- vst(dds, blind=F)

#MA plot, visualize number and significance of DE genes
plotMA(res,main='DE Genes', colSig="red", ylim=c(-4,4))

#identify genes manually
idx <- identify(res$baseMean, res$log2FoldChange)

#assign identified genes to new object
goi <- rownames(res)[idx]

#create PCA plots to visualize sample variance across group
plotPCA(vsd, intgroup="Tissue")
plotPCA(vsd, intgroup="Group")
plotPCA(vsd, intgroup="Instrument")
plotPCA(vsd, intgroup="library_preparation_method")
plotPCA(vsd, intgroup="AvgSpotLen")


###REDOING analysis with Tissue as design choice
dds_tissue <- DESeqDataSetFromMatrix(countData = round(counts),
                                     colData = coldata,
                                     design = ~Tissue)
#summarize dds object stats
dds_tissue

#set factor level
dds_tissue$Tissue <- relevel(dds_tissue$Tissue, ref = "Cortex Temporal")

#run DESeq function
dds_tissue <- DESeq(dds_tissue)

res_tissue <- results(dds_tissue)
res_tissue

summary(res_tissue)
rest0.01 <- results(dds_tissue, alpha = 0.01)
summary(rest0.01)

#extracting transformed values
vsdt <- vst(dds_tissue, blind=F)

#contrasts
resultsNames(dds_tissue)
#have not tried yet
rest2<-results(dds_tissue, contrast = c("Tissue","Cortex Temportal", "Cerebellum"))
vsdt2 <- vst(dds_tissue, blind=F)

#MA plot
plotMA(res,main='DE Genes', colSig="red", ylim=c(-4,4))

idx <- identify(res$baseMean, res$log2FoldChange)

goi <- rownames(res)[idx]

plotPCA(vsd, intgroup="Tissue")
plotPCA(vsd, intgroup="Group")
plotPCA(vsd, intgroup="Instrument")
plotPCA(vsd, intgroup="library_preparation_method")



#######DESeq2 analysis on CORTEX data only


#construct DESeq2 dataset
library(DESeq2)

dds.cor <- DESeqDataSetFromMatrix(countData = round(counts.cor),
                                  colData = col.cor,
                                  design = ~Group)
#summarize dds object stats
dds.cor
#set factor level
dds.cor$Group <- relevel(dds.cor$Group, ref = "Non-Neurological Control")

#run DESeq function
dds.cor <- DESeq(dds.cor)

#assign results to new object
res.cor <- results(dds.cor)
#peek at results
res.cor

#contrasts
res.cor <-results(dds.cor, contrast = c("Group","Non-Neurological Control", "ALS Spectrum MND"))

#extracting transformed values
#log transforming data
vsd.cor <- vst(dds.cor, blind=F)

summary(res.sp)

res.cor <- results(dds.cor, alpha = 0.001)

summary(res.sp0.01)

#how many genes had an padj less than 0.1?
sum(res.sp0.01$padj < 0.1, na.rm=TRUE)

sum(res.cor$padj < 0.001, na.rm=TRUE)

#MA plot
plotMA(res.cor,main='DE Genes in ALS Cortex vs Control\n alpha=0.001', colSig="red", ylim=c(-2,2))

idx <- identify(res$baseMean, res$log2FoldChange)

goi <- rownames(res)[idx]

#PCA plots through tissue and group categories
plotPCA(vsd.cor, intgroup="Tissue")
plotPCA(vsd.cor, intgroup="Group")


#######DESeq2 analysis on just spinal cord data

#construct DESeq2 dataset
library(DESeq2)

dds.sp <- DESeqDataSetFromMatrix(countData = round(counts.sp),
                                 colData = col.sp,
                                 design = ~Group)
#summarize dds object stats
dds.sp
#set factor level
dds.sp$Group <- relevel(dds.sp$Group, ref = "Non-Neurological Control")

#run DESeq function
dds.sp <- DESeq(dds.sp)

res.sp <- results(dds.sp)
res.sp

#contrasts

res.spc <-results(dds.sp, contrast = c("Group","Non-Neurological Control", "ALS Spectrum MND"))

vsdsp <- vst(dds.sp, blind=F)

summary(res.sp)

res.sp0.01 <- results(dds.sp, alpha = 0.01)

summary(res.sp0.01)

#how many genes had an padj less than 0.1?
sum(res.sp0.01$padj < 0.1, na.rm=TRUE)

sum(res.sp$padj < 0.05, na.rm=TRUE)


#extracting transformed values
#log transforming data
vsd.sp <- vst(dds.sp, blind=F)

#contrasts
resultsNames(dds)

res.sp<-results(dds.sp, contrast = c("Group","Non-Neurological Control", "ALS Spectrum MND"))

#MA plot
plotMA(res.sp,main='DE Genes', colSig="red", ylim=c(-4,4))

idx <- identify(res$baseMean, res$log2FoldChange)

goi <- rownames(res)[idx]

#PCA plots through tissue and group categories
plotPCA(vsdsp, intgroup="Tissue")
plotPCA(vsdsp, intgroup="Group")