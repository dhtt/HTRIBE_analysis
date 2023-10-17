# Load DESeq2 library
library(DESeq2)
library(data.table)
library(stringr)

# Set the working directory to where your count data is located
setwd('/home/dhthutrang/TRIBE/mRNA_seq/processed/alignment_salmon_ncbi/analyzed_result')

# # Read count data
# counts <- data.frame(fread("expression_data.csv", sep='\t', header = TRUE))
# countData <- counts[, 2:ncol(counts)]
# countData <- data.frame(lapply(countData, function(y) if(is.numeric(y)) round(y, 0) else y)) 
# rownames(countData) = counts$V1

# # Create a metadata table with sample information

# metadata <- data.frame(
#     sample = rownames(countData),
#     group = factor(sapply(rownames(countData), function(x) strsplit(x, '_')[[1]][1] ))
# )
# rownames(metadata) = rownames(countData)


# # Create a DESeq2 object
# countData <- t(countData)
# print(head(countData))
# print(metadata)
# dds <- DESeqDataSetFromMatrix(countData = countData,
#                               colData = metadata,
#                               design = ~ group)
# saveRDS(dds, file = "dds.RDS")

# # Perform normalization and differential expression analysis
# dds_normed <- DESeq(dds)
# saveRDS(dds_normed, file = "dds_normed.RDS")

# # Get differential expression results
# dds_results <- results(dds_normed)
# saveRDS(dds_results, file = "dds_results.RDS")

# # Filter the results for significantly differentially expressed genes
# significant_results <- subset(dds_results, padj < 0.05)
# saveRDS(significant_results, file = "significant_results.RDS")

dds <- readRDS("dds.RDS")

dds_normed <- readRDS("dds_normed.RDS")
dds_results <- readRDS("dds_results.RDS")
significant_results <- readRDS("significant_results.RDS")
print(significant_results)

# png("PCA_deseq2.png", width = 6, height = 6, unit = 'in')
# plotPCA(dds_rlog, intgroup="group", ntop=500) +
#   theme_bw() + # change theme
#   geom_point(size=5) + # add point size
#   ggtitle(label="Principal Component Analysis (PCA)", subtitle="Top 500 most variable genes")
# dev.off()