### Set your working directory 
setwd("/local/workdir/RNA-Seq_workshop/eag/exercise1")

### Read in the gene count file we made in last exercise.
matrixFile <- "/workdir/RNA-Seq_workshop/eag/exercise1/gene_counts.txt"

cts <- as.matrix(read.csv("gene_counts.txt", sep="\t", row.names=1, header=FALSE))

sampleFile <- "/workdir/RNA-Seq_workshop/eag/exercise1/samples.txt"

coldata <- read.csv("samples.txt", sep="\t", row.names=1)

head(coldata)

head(cts) 

### the columns have generic names (V1,V2,...), we can assign the row names of coldata to the column names of cts

colnames(cts) <- rownames(coldata)

head(cts) # looks good now


### Load package
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Treatment)

dds

vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("Treatment"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=Treatment)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text(aes(label=name),vjust=2)


### Differential expression
dds2 <- DESeq(dds)

res <- results(dds2)

summary(res)

### Sort summary list by p-value
res <- res[order(res$padj),]

head(res)









