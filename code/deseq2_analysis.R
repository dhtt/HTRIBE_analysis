# Load DESeq2 library
library(DESeq2)

# Set the working directory to where your count data is located
setwd('/home/dhthutrang/TRIBE/mRNA_seq/processed/alignment_salmon_ncbi')

# Read count data
counts <- read.table("expression_data.csv", header = TRUE, row.names = 1)
rep <- rep(['WT', 'mCherry', 'IMP2'], each = 3)
# Create a metadata table with sample information

metadata <- data.frame(
  sample = colnames(counts),
  group = factor(rep(['WT', 'mCherry', 'IMP2'], each = 3)),
  type = factor(rep(['control', 'control', 'test'], each = 3))
)
rownames(metadata) = rownames(counts)

# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ group + type)

# Perform normalization and differential expression analysis
dds <- DESeq(dds)

# Get differential expression results
results <- results(dds)

# Filter the results for significantly differentially expressed genes
significant_results <- subset(results, padj < 0.05)

# Print the significant results
print(significant_results)
