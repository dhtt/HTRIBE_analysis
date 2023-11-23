# Load DESeq2 library
library(DESeq2)
library(data.table)
library(stringr)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(rtracklayer)
library(dplyr)
library(biomaRt)
library(pheatmap)

jaccard_index <- function(list1, list2){
    print(paste("Common genes: ", paste(unique(intersect(list1, list2)), collapse=', '), collapse=""))
    print(paste("No common genes: ", length(unique(intersect(list1, list2))), collapse=""))
    print(paste("No all genes: ", length(unique(union(list1, list2))), collapse=""))
    return(length(intersect(list1, list2))/length(union(list1, list2)))
}

# Set the working directory
setwd('/home/dhthutrang/TRIBE/mRNA_seq/processed/alignment_salmon_ncbi/analyzed_result')


################ DESEQ2 ANALYSIS ################ 
# Read count data
counts <- data.frame(fread("expression_data.csv", sep='\t', header = TRUE))
countData <- counts[, 2:ncol(counts)]
countData <- data.frame(lapply(countData, function(y) if(is.numeric(y)) round(y, 0) else y))
rownames(countData) = counts$V1

# Create a metadata table with sample information
metadata <- data.frame(
    sample = rownames(countData),
    group = factor(sapply(rownames(countData), function(x) strsplit(x, '_')[[1]][1] ))
)
rownames(metadata) = rownames(countData)

# Create a DESeq2 object
countData <- t(countData)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = metadata,
                              design = ~ group)
saveRDS(dds, file = "dds.RDS")

# Perform normalization and differential expression analysis
dds_normed <- DESeq(dds)
saveRDS(dds_normed, file = "dds_normed.RDS")


################ RESULTS ANALYSIS ################ 
# Filter the results for significantly differentially expressed genes
dds <- readRDS("analyzed_result/dds.RDS")
dds_normed <- readRDS("analyzed_result/dds_normed.RDS")
normed_data = counts(dds_normed, normalized = TRUE)

dds_result_multifactorial <- results(dds_normed)
dds_result_WT_mCherry <- results(dds_normed, contrast = c('group', 'WT', 'mCherry'))
dds_result_WT_IMP2 <- results(dds_normed, contrast = c('group', 'WT', 'IMP2'))
dds_result_mCherry_IMP2 <- results(dds_normed, contrast = c('group', 'mCherry', 'IMP2'))

vsd = vst(dds)

# PCA plot for transformed counts 
png("analyzed_result/PCA_deseq2.png", width = 6, height = 4, unit = 'in', res = 200)
plotPCA(vsd, intgroup="group", ntop=500) +
    theme_bw() +
    geom_point(size = 4, alpha = 0.8) +
    ggtitle(label="Principal Component Analysis (PCA)", subtitle="Top 500 most variable genes")
dev.off()

# MAplot for pair comparison with transformed counts
png("analyzed_result/MA_deseq2.png", width = 8, height = 8, unit = 'in', res = 200)
par(mfrow=c(2, 2))
plotMA(dds_result_WT_mCherry, main="WT vs mCherry", alpha=0.05)
plotMA(dds_result_WT_IMP2, main="WT vs IMP2", alpha=0.05)
plotMA(dds_result_mCherry_IMP2, main="mCherry vs IMP2", alpha=0.05)
dev.off()

# Heatmap for pair comparison with transformed counts
res_list = list(dds_result_WT_IMP2, dds_result_mCherry_IMP2, dds_result_WT_mCherry, dds_result_multifactorial)
names(res_list) = c("WT_IMP2", "mCherry_IMP2", "WT_mCherry", 'WT/mCherry_IMP2 - WT_mCherry')
sig_gene_list = list() # List of DEGs for 3 comparisons 
for (i in 1:length(res_list)){
    res = res_list[[i]]
    comparison_type = names(res_list)[i]
    comparison_pair = strsplit(comparison_type, "_")[[1]]
    res_lfc <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
    print(paste("Comparison type: ", comparison_type))
    print(paste("Number of DEGs: ", nrow(res_lfc)))
    genes <- rownames(res_lfc)[order(res_lfc$padj, decreasing=TRUE)] 
    sig_gene_list[[i]] = genes
    
    vst_sig <- vsd[rownames(vsd) %in% genes, ]
    plot_title = paste("Differentially Expression Transcripts for ", comparison_pair[1], " (", nrow(res_lfc), " transcripts)", collapse="")
    
    heat <- t(scale(t(assay(vst_sig))))
    heat_rownames = unique(sapply(rownames(heat), function(x) strsplit(x, '.', fixed = T)[[1]][1]))

    png(paste("analyzed_result/heatmap_deseq2_toppadj_", comparison_type, ".png", sep=''), width = 7, height = 8, unit = 'in', res = 200)
    pheatmap(heat, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=FALSE, main=plot_title)
    dev.off()
}

names(sig_gene_list) = names(res_list)

# Map transcript IDs of DEGs to gene IDs using reference genome file 
refgen = import.gff("/Users/trangdo/Documents/BIOINFO/TRIBE/HTRIBE_analysis/refgen/refGene.gtf")
transcript_table = refgen %>%
    as.data.frame() %>%
    dplyr::select(transcript_id, gene_id) %>%
    unique()

map_gene_id <- function(gene_list, refgen_table, strip_char=T, unique=T, reverse=F){
    if (reverse){
        if (strip_char){
            gene_list = sapply(gene_list, function(x) strsplit(x, '.', fixed = T)[[1]][1])
        }
        gene_list = sapply(gene_list, function(x) refgen_table$gene_id[grep(x, refgen_table$transcript_id)])
        if (unique){
            gene_list = unique(gene_list)
        } 
    } else {
        gene_list = unlist(sapply(gene_list, function(x) refgen_table$transcript_id[grep(x, refgen_table$gene_id)]))
    }
    return(gene_list)
}
sig_genes_id = lapply(sig_gene_list, function(x) {
    sig_genes = unique(sapply(x, function(x) strsplit(x, '.', fixed = T)[[1]][1])) 
    sig_genes = map_gene_id(sig_genes, transcript_table)
    return(sig_genes)
})

sig_genes_id[[4]] = setdiff(union(sig_genes_id$mCherry_IMP2, sig_genes_id$WT_IMP2), sig_genes_id$WT_mCherry)
sig_genes_id[[5]] = union(sig_genes_id$mCherry_IMP2, sig_genes_id$WT_IMP2)
names(sig_genes_id) = c(names(sig_gene_list)[2:4], 'WT/mCherry_IMP2 - WT_mCherry', 'WT/mCherry_IMP2')

################ SIMILARITY BETWEEN EACH OTHER AND TO HTRIBE RESULTS ################ 
lapply(sig_genes_id, function(x) lapply(x, function(y) jaccard_index(unlist(x), unlist(y))))
jaccard_index(unlist(sig_genes_id$mCherry_IMP2), unlist(sig_genes_id$WT_mCherry))

sig_genes_HTRIBE <- read.delim("~/Documents/BIOINFO/TRIBE/HTRIBE_analysis/HTRIBE_sig_gene_1%/sig_genesCDS_UTR.txt", header=FALSE, comment.char="#")
ground_truth_lists = list(sig_genes_HTRIBE[1, ], sig_genes_HTRIBE[2, ], sig_genes_HTRIBE[3, ])
ground_truth_lists = lapply(sig_genes_HTRIBE, function(x){return(strsplit(x, ','))})[[1]]
names(ground_truth_lists) = c('WT_IMP2', 'mCherry_IMP2', 'WT_mCherry')
ground_truth_lists[[4]] = setdiff(union(ground_truth_lists$mCherry_IMP2, ground_truth_lists$WT_IMP2), ground_truth_lists$WT_mCherry)
ground_truth_lists[[5]] = union(ground_truth_lists$mCherry_IMP2, ground_truth_lists$WT_IMP2)
names(ground_truth_lists)[4] = 'WT/mCherry_IMP2 - WT_mCherry'
names(ground_truth_lists)[5] = 'WT/mCherry_IMP2'
lapply(ground_truth_lists, length)
lapply(sig_genes_id, length)


round(jaccard_index(sig_genes_id$WT_IMP2, ground_truth_lists$WT_IMP2)[[1]]*1000, 2)
round(jaccard_index(sig_genes_id$mCherry_IMP2, ground_truth_lists$mCherry_IMP2)[[1]]*1000, 2)
round(jaccard_index(sig_genes_id$WT_mCherry, ground_truth_lists$WT_mCherry)[[1]]*1000, 2)
round(jaccard_index(sig_genes_id$`WT/mCherry_IMP2 - WT_mCherry`, ground_truth_lists$`WT/mCherry_IMP2 - WT_mCherry`)[[1]]*1000, 2)
round(jaccard_index(sig_genes_id$`WT/mCherry_IMP2`, ground_truth_lists$`WT/mCherry_IMP2`)[[1]]*1000, 2)



control_genes = unlist(union(intersect(ground_truth_lists$WT_IMP2, sig_genes_id$WT_IMP2), intersect(ground_truth_lists$mCherry_IMP2, sig_genes_id$mCherry_IMP2)))
IMP2_genes = unlist(setdiff(
    unlist(union(intersect(ground_truth_lists$WT_IMP2, sig_genes_id$WT_IMP2), intersect(ground_truth_lists$mCherry_IMP2, sig_genes_id$mCherry_IMP2))),
    unlist(intersect(ground_truth_lists$WT_mCherry, sig_genes_id$WT_mCherry))
    ))
IMP2_genes_method1 = unlist(intersect(sig_genes_id$`WT/mCherry_IMP2 - WT_mCherry`, ground_truth_lists$`WT/mCherry_IMP2 - WT_mCherry`))
HTRIBE_DESeq_intersect_genes_list = list(control_genes, IMP2_genes, IMP2_genes_method1)
names(HTRIBE_DESeq_intersect_genes_list) = c('B: WT_mCherry', 'B: WT/mCherry_IMP2', 'A: WT/mCherry_IMP2')
################ ENRICHMENT ANALYSIS ################ 
library(org.Mm.eg.db)
library(clusterProfiler)
library(viridis)
library(ggplot2)
library(ggpubr)

create_enrichment_plot <- function(gene_list_){
    enriched_terms = list()
    enriched_terms_names = c()
    
    for (i in 1:length(gene_list_)){
        gene_set = unlist(gene_list_[[i]])
        gene_set_name = paste(strsplit(names(gene_list_)[i], '_')[[1]], collapse = ' vs. ')
        enriched_terms_names[i] = gene_set_name
        
        gene_set_Entrez = AnnotationDbi::select(org.Mm.eg.db, 
                                                keys = gene_set,
                                                columns = c("ENTREZID", "SYMBOL"),
                                                keytype = "SYMBOL")
        ego <- enrichGO(gene = gene_set_Entrez$ENTREZID, OrgDb = org.Mm.eg.db, ont="BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05)
        enriched_terms[[i]] = ego
    }
    names(enriched_terms) = enriched_terms_names
    
    enriched_plots = list()
    for (i in 1:length(enriched_terms)){
        data = enriched_terms[[i]]
        gene_ratio = sapply(data$GeneRatio, function(x){
            faction = strsplit(x, "/")[[1]]
            ratio = as.numeric(faction[1])/as.numeric(faction[2])
            return(ratio)
        })
        gene_ratio = head(sort(gene_ratio, decreasing = TRUE), 15)
        plot = enrichplot::dotplot(enriched_terms[[i]], showCategory=15) +
            scale_color_viridis(option = "D") +
            scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 70)) +
            scale_x_continuous(breaks=seq(0, max(gene_ratio), round(max(gene_ratio)/4, 2))) +
            ggtitle(names(enriched_terms)[i]) +
            theme_light() +
            theme(
                aspect.ratio = 1,
                plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
                axis.text.x = element_text(colour="black", size = 12),
                axis.text.y = element_text(colour="black", size = 11),
                plot.margin = unit(c(20,20,20,20), "pt"))
        enriched_plots[[i]] = plot
    }
    names(enriched_plots) = enriched_terms_names
    return(enriched_plots)
}

ground_truth_lists_enriched_plots = create_enrichment_plot(ground_truth_lists)
sig_gene_ids_enriched_plots = create_enrichment_plot(sig_genes_id)
HTRIBE_DESeq_intersect_genes_list_enriched_plots = create_enrichment_plot(HTRIBE_DESeq_intersect_genes_list)

all_plots = list(ground_truth_lists_enriched_plots[[1]], sig_gene_ids_enriched_plots[[1]],
                 ground_truth_lists_enriched_plots[[2]], sig_gene_ids_enriched_plots[[2]],
                 ground_truth_lists_enriched_plots[[3]], sig_gene_ids_enriched_plots[[3]],
                 ground_truth_lists_enriched_plots[[4]], sig_gene_ids_enriched_plots[[4]]
                 )

png(paste("analyzed_result/GO_annots_HTRIBE.png", sep=''), width = 16, height = 8, unit = 'in', res = 200)
ggarrange(plotlist = ground_truth_lists_enriched_plots, ncol=2, nrow=2, common.legend = TRUE)
dev.off()
png(paste("analyzed_result/GO_annots_DEG.png", sep=''), width = 16, height = 8, unit = 'in', res = 200)
ggarrange(plotlist = sig_gene_ids_enriched_plots, ncol=2, nrow=2, common.legend = TRUE)
dev.off()
png(paste("analyzed_result/GO_annots_HTRIBE_DESeq_intersect.png", sep=''), width = 16, height = 8, unit = 'in', res = 200)
ggarrange(plotlist = HTRIBE_DESeq_intersect_genes_list_enriched_plots, ncol=2, nrow=2, common.legend = TRUE)
dev.off()
png(paste("analyzed_result/GO_annots_all_plots.png", sep=''), width = 16, height = 16, unit = 'in', res = 200)
ggarrange(plotlist = all_plots, ncol=2, nrow=4, common.legend = TRUE)
dev.off()

################# CUMULATIVE PLOT ################ 
prepare_plot_ecdf <- function(HTRIBE_gene_set, DESEQ_result, transcript_table){
    HTRIBE_genes = HTRIBE_gene_set
    DESEQ_df = DESEQ_result
    DESEQ_df = DESEQ_df[!is.na(DESEQ_df$padj), ]
    DESEQ_df_rownames = sapply(rownames(DESEQ_df), function(x) {return(strsplit(x, ".", fixed=T)[[1]][1])})
    transcript_list = unlist(sapply(HTRIBE_genes, function(x) transcript_table$transcript_id[grep(x, transcript_table$gene_id)]))
    DESEQ_df$has_A2G = "IMP2-" 
    DESEQ_df$has_A2G[DESEQ_df_rownames %in% transcript_list] = "IMP2+"
    return(as.data.frame(DESEQ_df))
}


ecdf_df_list = list()
comparisons = names(ground_truth_lists)
comparisons_label = sapply(comparisons, function(x){
    return(paste(strsplit(x, split = '_')[[1]], collapse = " vs. "))
})

for (i in 1:length(comparisons)){
    comparison_type = comparisons[i]
    ecdf_df_list[[i]] = prepare_plot_ecdf(ground_truth_lists[[comparison_type]], res_list[[comparison_type]], transcript_table)
}
names(ecdf_df_list) = comparisons_label

plot_ecdf <- function(all_LFC, plot_title){
    ecdf <- ggplot(as.data.frame(all_LFC), aes(log2FoldChange, colour = has_A2G)) + 
        stat_ecdf(aes(colour = has_A2G)) +
        xlab("Log Fold Change") + ylab("Probability") +
        scale_color_discrete(name = "Type") +
        theme_minimal() +
        ggtitle(plot_title) + 
        theme(  
            aspect.ratio = 1,
            plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
            axis.text.x = element_text(colour="black", size = 12),
            axis.text.y = element_text(colour="black", size = 12),
            axis.title = element_text(face = "bold"),
            panel.border = element_rect(colour = "grey", fill = NA, size = 0.5),
            legend.position = "bottom") +
        guides(fill = guide_legend(nrow = 3, byrow = t))
    return(ecdf)
}

ecdf_list = list()
for (i in 1:length(ecdf_df_list)){
    ecdf_list[[i]] = plot_ecdf(ecdf_df_list[[i]], names(ecdf_df_list)[i])
}

png(paste("analyzed_result/ECDF_LFC.png", sep=''), width = 8, height = 8, unit = 'in', res = 200)
ggarrange(plotlist = ecdf_list, ncol=2, nrow=2, common.legend = TRUE)
dev.off()


ks_test = ks.test(ecdf_df_list$`WT/mCherry vs. IMP2 - WT vs. mCherry`[ecdf_df_list$`WT/mCherry vs. IMP2 - WT vs. mCherry`$has_A2G == "IMP2+", "log2FoldChange"], 
                  ecdf_df_list$`WT/mCherry vs. IMP2 - WT vs. mCherry`[ecdf_df_list$`WT/mCherry vs. IMP2 - WT vs. mCherry`$has_A2G ==  "IMP2-", "log2FoldChange"])
ks_test$p.value

plot = ecdf_list[[1]]
plot_ks = plot + 
    geom_text(label = paste('KS test p-value'), x = -20, y = 0.8)
plot_ks
