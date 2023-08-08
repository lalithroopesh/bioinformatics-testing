setwd("C:/Users/lalit/Desktop/DEGanalysis")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("ashr")
install.packages("pheatmap")

library(DESeq2)

countsdata <- read.delim("Soybean_counts.txt", row.names="gene_id")
head(countsdata)

coldata <- read.delim("sample_info.txt", row.names = 1)
View(coldata)

all(colnames(countsdata) %in% rownames(coldata))
# -> TRUE, they are in it
all(colnames(countsdata) == rownames(coldata))
# -> TRUE, they are in the same order

dds <- DESeqDataSetFromMatrix(
  countData = countsdata,
  colData = coldata,
  design = ~ condition
)

# sets factor level, compares to untreated
dds$condition <- relevel(dds$condition, ref = "untreated")

# pre-filtering to only keep rows with 10 reads total
keep <- rowSums(counts(dds)) >= 10
keep # -> list of rows
dds[keep,]

# run DESeq
dds <- DESeq(dds)
res <- results(dds)
res

# summarize results
summary(res)
resDS <- results(dds, alpha=0.05)
summary(resDS)

plotMA(res)

# PCA analysis
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("cultivar", "condition"))

# log2(n + 1) transform
ntd <- normTransform(dds)

# Run heatmap
library(pheatmap)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("cultivar", "condition")])

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# export output
write.csv(as.data.frame(res), file="condition_treated_results.csv")

# output genes with p < 0.05
resSig <- subset(res, padj < 0.05)
write.csv(as.data.frame(resSig), 
          file="condition_treated_results05.csv")

# extract the normalized counts for expression levels
dds.countsdata <- estimateSizeFactors(dds)
NormalizedCount <- counts(dds.countsdata, normalized = T)
write.csv(NormalizedCount, file="condition_treated_normcounts.csv", row.names = TRUE)