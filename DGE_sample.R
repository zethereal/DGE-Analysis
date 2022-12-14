

######DIFFERENTIAL GENE EXPRESSION ANALYSIS simple

#using DESeq2 package. Can also use limma/voom or edgeR
#For DESeq2, start with two objects, counts matrix, and column data table (see vignette)


#set working directory
setwd("W:/nl/umw_jonathan_watts/nl_Zach/fastq/raw/htseq/deseq/analysis")

#what files are in this working directory
list.files()

#import file
APO3_11 <- read.delim("W:/nl/umw_jonathan_watts/nl_Zach/fastq/raw/htseq/deseq/analysis/APO3ss01011_count_tab.txt",header=F,sep='\t')
APO3_12 <- read.delim("W:/nl/umw_jonathan_watts/nl_Zach/fastq/raw/htseq/deseq/analysis/APO3ss01012_count_tab.txt",header=F,sep='\t')

NTC31 <- read.delim("W:/nl/umw_jonathan_watts/nl_Zach/fastq/raw/htseq/deseq/analysis/NTC31_count_tab.txt",header=F,sep='\t')
NTC32 <- read.delim("W:/nl/umw_jonathan_watts/nl_Zach/fastq/raw/htseq/deseq/analysis/NTC32_count_tab.txt",header=F,sep='\t')

#view imported file
#make sure format is correct
View(APO3_11)
class(APO3)

#make new object and combine all count tables
count_full <- cbind(APO3_11, APO3_12$V2, NTC31$V2, NTC32$V2)

#rename merged table columns that make sense
colnames(count_full) <- c("gene", "APO3_11", "APO3_12", "NTC31", "NTC32")

##NEED column names of the counts matrix to match rownames of coldata 
#must be in SAME ORDER as well

#create coldata table with sample information and import
coldata <- read.csv("W:/nl/umw_jonathan_watts/nl_Zach/fastq/raw/htseq/deseq/analysis/coldata_sample.csv",header=T,sep=',')

#give coldata rownames that match 1st column named ID
rownames(coldata) <- coldata$ID

#same to count matrix
rownames(count_full) <- count_full$gene

#remove first column of count data
count_full <- count_full[,-1]

#do row names in coldata match column names in counts data
all(colnames(count_full) %in% rownames(coldata)) #TRUE
#are they in the same order?
all(colnames(count_full) == rownames(coldata)) # TRUE

#if both are true, move ahead to importing samples into DESeq

#load package
BiocManager::install("DESeq2")
library(DESeq2)

#construct DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_full,
                              colData = coldata,
                              design = ~treatment)
#summarize dds object stats
dds

#set factor level, specify 'untreated' or it will choose the first group
dds$treatment <- relevel(dds$treatment, ref = "untreated")

#run DESeq function- may take a while depending on dataset size
dds <- DESeq(dds)

#assign new results to object
res <- results(dds)
res

#summarize results
summary(res)

#how many genes had an padj less than 0.1?
sum(res$padj < 0.5, na.rm=TRUE)

#make new results object with different FDR cutoff (alpha)
res <- results(dds, alpha = 0.5)
summary(res0.01)

#extracting transformed values
#log transforming data
vsd <- vst(dds, blind=F)

#set contrasts
#specify group category name, untreated, then treated samples to compare
res<-results(dds, contrast = c("treatment","untreated", "treated"))
vsd2 <- vst(dds, blind=F)

#MA plot, visualize number and significance of DE genes
plotMA(res,main='MA Plot\nDE Genes', colSig="red")

#identify genes manually
idx <- identify(res$baseMean, res$log2FoldChange)

#assign identified genes to new object
goi <- rownames(res)[idx]
goi

#create PCA plots to visualize sample variance across group
plotPCA(vsd, intgroup="treatment")

###volcano plots
#import and install package

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

#making volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 1e-5,
                FCcutoff = 0.5,
                labSize=4,
                boxedLabels = T,
                drawConnectors = T)


