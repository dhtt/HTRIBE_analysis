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
library(scales)
library(cowplot)
library(org.Mm.eg.db)
library(clusterProfiler)
library(viridis)
library(ggplot2)
library(ggpubr)
library("VennDiagram")

diverging_pallete = colorRampPalette(rev(c("#1a9850", "#a6d96a", "#ffffe4", "#F1C1C1", "#D61818")))(99)
diverging_pallete = c(colorRampPalette(c("#EA6B66", "#EA6B66", "#F1C1C1", "#ffffe4"))(50), 'white', colorRampPalette(c("#ffffe4", "#a6d96a", "#1a9850"))(50))

# Set the working directory
setwd('/home/dhthutrang/TRIBE/mRNA_seq/processed/alignment_salmon_ncbi/analyzed_result')
setwd('~/Documents/BIOINFO/TRIBE/HTRIBE_analysis/analyzed_result/DEGs')


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
#### Filter the results for significantly differentially expressed genes ####
dds <- readRDS("dds.RDS")
dds_normed <- readRDS("dds_normed.RDS")
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

#### MAplot for pair comparison with transformed counts ####
make_MAplot <- function(deseq2_result, plot_title){
    plot = ggmaplot(deseq2_result, main = plot_title, 
             fdr = 0.05, fc = 1, size = 0.5, top = 0,
             xlab = expression(bold("Log"["2"] ~ "mean expression")),
             ylab = expression(bold("Log"["2"] ~ "fold change")),
             palette = c("#64BF52", "#EA6B66", "lightgray"), alpha = 0.6,
             legend = "top", font.legend = c("bold", 12), font.main = "bold", ggtheme = theme_light()) + 
        theme(
            panel.background = element_rect(colour = "grey", linewidth = 1),
            legend.position = c(0.75, 0.85),
            legend.background=element_blank(),
            legend.text = element_text(size=10),
            plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
            panel.grid = element_blank()
        ) + 
        guides(colour = guide_legend(override.aes = list(size=3, alpha=0.8)))
    return(plot)
}

#### Gather DESeq2 result ####
res_list = list(dds_result_WT_IMP2, dds_result_mCherry_IMP2, dds_result_WT_mCherry, dds_result_multifactorial)
names(res_list) = c("WT_IMP2", "mCherry_IMP2", "WT_mCherry", 'WT/mCherry_IMP2 - WT_mCherry')

#### Make MA plot ####
MAplots_label = paste(LETTERS[1:length(res_list)], 
                      sapply(names(res_list), function(x) {return(paste(strsplit(x, "_")[[1]], collapse = " vs. "))}), 
                      sep='. ')
MAplots = list()
for (i in 1:length(res_list)){
    MAplots[[i]] = make_MAplot(res_list[[i]], MAplots_label[i])
}
names(MAplots) = names(res_list)

png("analyzed_result/MA_deseq2.png", width = 10, height = 3.5, unit = 'in', res = 200)
ggarrange(MAplots$WT_IMP2, MAplots$mCherry_IMP2, MAplots$WT_mCherry, ncol = 3, nrow = 1, align = "h")
dev.off()

#### Heatmap for pair comparison with transformed counts ####
sig_gene_list = list() # List of DEGs for 3 comparisons 
heatmaps = list()

for (i in 1:length(res_list)){
    res = res_list[[i]]
    comparison_type = names(res_list)[i]
    comparison_pair = strsplit(comparison_type, "_")[[1]]
    res_lfc <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
    print(paste("Comparison type: ", comparison_type))
    print(paste("Number of DEGs: ", nrow(res_lfc)))
    
    genes <- rownames(res_lfc)[order(res_lfc$padj, decreasing=TRUE)] 
    sig_gene_list[[i]] = genes
    
    vst_sig <- vsd[rownames(vsd) %in% genes, grep(paste(comparison_pair, collapse="|"), colnames(vsd))]
    heat <- t(scale(t(assay(vst_sig))))
    colnames(heat) = gsub("_", " ", colnames(heat))
    
    plot_title = paste(comparison_pair[1], " vs. ", comparison_pair[2], " (", nrow(res_lfc), " transcripts)", sep="")
    if (i == (length(res_list)-1)){ legend = TRUE } else { legend = FALSE}
    
    # png(paste("analyzed_result/heatmap_deseq2_toppadj_", comparison_type, ".png", sep=''), width = 7, height = 8, unit = 'in', res = 200)
    heatmaps[[i]] = pheatmap(heat, cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=FALSE, 
                             breaks=seq(-2.2, 2.2, (4.4/101)), fontsize.x = 14, legend=legend, 
                             color=diverging_pallete, main=plot_title)
    # dev.off()
}

png(paste("heatmap_transcripts/heatmap_deseq2_toppadj_allgroups.png", sep=''), width = 11, height = 6, unit = 'in', res = 200)
ggarrange(plotlist = list(heatmaps[[1]]$gtable, heatmaps[[2]]$gtable, heatmaps[[3]]$gtable), ncol = 3, nrow = 1, labels = c("A", "B", "C"), align="hv", vjust=0.5) +
    theme(plot.margin = margin(1, 0.5, 1, 0.5, "cm")) 
dev.off()
names(sig_gene_list) = names(res_list)

#### Map transcript IDs of DEGs to gene IDs using reference genome file ####
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
    sig_genes = map_gene_id(sig_genes, transcript_table, reverse = T)
    return(sig_genes)
})

sig_genes_id[[4]] = setdiff(union(sig_genes_id$mCherry_IMP2, sig_genes_id$WT_IMP2), sig_genes_id$WT_mCherry)
sig_genes_id[[5]] = union(sig_genes_id$mCherry_IMP2, sig_genes_id$WT_IMP2)
names(sig_genes_id) = c(names(sig_gene_list)[1:3], 'WT/mCherry_IMP2 - WT_mCherry', 'WT/mCherry_IMP2')


################ SIMILARITY BETWEEN EACH OTHER AND TO HTRIBE RESULTS ################ 
#### Prepare datasets ####
sig_genes_HTRIBE <- read.delim("~/Documents/BIOINFO/TRIBE/HTRIBE_analysis/HTRIBE_sig_gene_1%/sig_genesCDS_UTR.txt", header=FALSE, comment.char="#")
ground_truth_lists = lapply(sig_genes_HTRIBE, function(x){return(strsplit(x, ','))})[[1]]
names(ground_truth_lists) = c('WT_IMP2', 'mCherry_IMP2', 'WT_mCherry')
ground_truth_lists[[4]] = setdiff(union(ground_truth_lists$mCherry_IMP2, ground_truth_lists$WT_IMP2), ground_truth_lists$WT_mCherry)
ground_truth_lists[[5]] = union(ground_truth_lists$mCherry_IMP2, ground_truth_lists$WT_IMP2)
names(ground_truth_lists)[4:5] = c('WT/mCherry_IMP2 - WT_mCherry', 'WT/mCherry_IMP2')

#### Compute Jaccard index ####
jaccard_index <- function(list1, list2){
    print(paste("Common genes: ", paste(unique(intersect(list1, list2)), collapse=', '), collapse=""))
    print(paste("No common genes: ", length(unique(intersect(list1, list2))), collapse=""))
    print(paste("No all genes: ", length(unique(union(list1, list2))), collapse=""))
    return(length(intersect(list1, list2))/length(union(list1, list2)))
}

round(jaccard_index(sig_genes_id$WT_IMP2, ground_truth_lists$WT_IMP2)[[1]]*1000, 2)
round(jaccard_index(sig_genes_id$mCherry_IMP2, ground_truth_lists$mCherry_IMP2)[[1]]*1000, 2)
round(jaccard_index(sig_genes_id$WT_mCherry, ground_truth_lists$WT_mCherry)[[1]]*1000, 2)
round(jaccard_index(sig_genes_id$`WT/mCherry_IMP2 - WT_mCherry`, ground_truth_lists$`WT/mCherry_IMP2 - WT_mCherry`)[[1]]*1000, 2)
round(jaccard_index(sig_genes_id$`WT/mCherry_IMP2`, ground_truth_lists$`WT/mCherry_IMP2`)[[1]]*1000, 2)

#### Compare overlapping methods ####
control_genes = unlist(union(intersect(ground_truth_lists$WT_IMP2, sig_genes_id$WT_IMP2), intersect(ground_truth_lists$mCherry_IMP2, sig_genes_id$mCherry_IMP2)))
IMP2_genes = unlist(setdiff(
    unlist(union(intersect(ground_truth_lists$WT_IMP2, sig_genes_id$WT_IMP2), intersect(ground_truth_lists$mCherry_IMP2, sig_genes_id$mCherry_IMP2))),
    unlist(intersect(ground_truth_lists$WT_mCherry, sig_genes_id$WT_mCherry))
))
IMP2_genes_method1 = unlist(intersect(sig_genes_id$`WT/mCherry_IMP2 - WT_mCherry`, ground_truth_lists$`WT/mCherry_IMP2 - WT_mCherry`))
HTRIBE_DESeq_intersect_genes_list = list(control_genes, IMP2_genes, IMP2_genes_method1)
names(HTRIBE_DESeq_intersect_genes_list) = c('B: WT_mCherry', 'B: WT/mCherry_IMP2', 'A: WT/mCherry_IMP2')

#### Compare the overlaps between different windows ####
HTRIBE_DEseq2_overlaps = list(
    c("Crem","Ctnnd1","Taf6","Rabep1","Met","Tmem62","Cplane1","Ccdc9","Ddx58"),
    c("Sult1b1","Ctnnd1","Med16","Tnpo1","Sntg2","Taf6","Rabep1","Tmem62","Macf1","Trim39","Dck","Cplane1","Agfg1","Mug-ps1","Ccdc9"),
    c("Med16","Tardbp","Tnpo1","Mad2l2","Creb1","Taf6","Rabep1","Tmem62","Trim39","Cyp2c38","Dck","Ccdc90b","Agfg1"),
    c("Bdnf","Ctnnd1","Tardbp","Tnpo1","Ap1b1","Sntg2","Mlh3","Rnf128","Mad2l2","Creb1","Trappc13","Rabep1","Tmem62","Serinc2","Mars1","Accs","Trim39","Ugdh","Cyp2c38","Cybc1","Ccdc90b","Foxk1","Agfg1","Ccdc137","Ifi35","Ccdc9")
)
spans = c("Site", "Window", "CDS/UTR", "Transcript")
venn.plot <- venn.diagram(x=HTRIBE_DEseq2_overlaps, 
                          category.names=spans, 
                          filename='analyzed_result/DEGs/HTRIBE_DESeq2overlap.png',
                          output=TRUE, imagetype="png", height=1000, width=1200, resolution=300,
                          compression = "lzw",  lwd = 2, lty = 'blank',
                          fill = c("#a6d96a", "#F1C1C1", "#FFD966", "#C3ABD0"), cat.fontface = "bold", cat.default.pos = "outer",
                          cat.fontfamily = "sans", fontfamily="sans"
                          )

################ ENRICHMENT ANALYSIS ################ 
get_enrichment_data <- function(gene_list_){
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
    return(enriched_terms)
}

create_enrichment_plot <- function(enriched_terms){
    enriched_plots = list()
    for (i in 1:length(enriched_terms)){
        data = enriched_terms[[i]]
        gene_ratio = sapply(data$GeneRatio, function(x){
            faction = strsplit(x, "/")[[1]]
            ratio = as.numeric(faction[1])/as.numeric(faction[2])
            return(ratio)
        })
        gene_ratio = head(sort(gene_ratio, decreasing = TRUE), 15)
        
        if (i %% 2 == 0){
            margin.right = 0.5
            margin.left = -0.5
        } else {
            margin.right = -0.5
            margin.left = 0.5
        }
        plot = enrichplot::dotplot(enriched_terms[[i]], showCategory=8) +
            scale_fill_gradient2(low="#1a9850", high="#FFD966", mid="#a6d96a",
                                 midpoint = 0.000005, guide = "colourbar",
                                 name="Adjusted p-value") +
            scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 70)) +
            scale_x_continuous(breaks = scales::pretty_breaks(n = 2),
                               expand = expand_scale(add = c(0.005, 0.005)), name="Gene Ratio") +
            theme_light() +
            theme(
                legend.text = element_text(size=15),
                legend.title = element_text(size=16),
                legend.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
                legend.key.width = unit(dev.size()[1]/20, "in"),
                axis.title = element_text(size = 15),
                axis.text = element_text(colour="black", size = 17),
                plot.margin = margin(0.5, margin.right, 0.5, margin.left, "cm")
                )
        
        enriched_plots[[i]] = plot
    }
    names(enriched_plots) = names(enriched_terms)
    return(enriched_plots)
}


#### Get data and create plots ####
ground_truth_lists_enriched_data = get_enrichment_data(ground_truth_lists)
sig_gene_ids_enriched_data = get_enrichment_data(sig_genes_id)
HTRIBE_DESeq_intersect_genes_list_enriched_data = get_enrichment_data(HTRIBE_DESeq_intersect_genes_list)

ground_truth_lists_enriched_plots = create_enrichment_plot(ground_truth_lists_enriched_data)
sig_gene_ids_enriched_plots = create_enrichment_plot(sig_gene_ids_enriched_data)
HTRIBE_DESeq_intersect_genes_list_enriched_plots = create_enrichment_plot(HTRIBE_DESeq_intersect_genes_list_enriched_data)

all_plots = list(ground_truth_lists_enriched_plots$`WT vs. IMP2`, sig_gene_ids_enriched_plots$`WT vs. IMP2`,
                 ground_truth_lists_enriched_plots$`mCherry vs. IMP2`, sig_gene_ids_enriched_plots$`mCherry vs. IMP2`,
                 ground_truth_lists_enriched_plots$`WT vs. mCherry`, sig_gene_ids_enriched_plots$`WT vs. mCherry`,
                 ground_truth_lists_enriched_plots$`WT/mCherry vs. IMP2 - WT vs. mCherry`, sig_gene_ids_enriched_plots$`WT/mCherry vs. IMP2 - WT vs. mCherry`)

png(paste("analyzed_result/DEGs/GO_annots_HTRIBE_DESeq_intersect.png", sep=''), width = 17, height = 8, unit = 'in', res = 200)
ggarrange(plotlist = HTRIBE_DESeq_intersect_genes_list_enriched_plots, ncol=2, nrow=2, common.legend = TRUE)
dev.off()
png(paste("analyzed_result/DEGs/GO_annots_all_plots.png", sep=''), width = 18, height = 14, unit = 'in', res = 200)
ggarrange(plotlist = all_plots, ncol=2, nrow=4, common.legend = TRUE, align="hv",
          label.x = c(0, 0, -0.04, -0.03, -0.03, -0.02, -0.18, -0.17) - 0.06,
          font.label = list(size = 17),
          labels = paste(LETTERS[unlist(lapply(seq(1, 4, 1), function(x) {return(c(x, x + 4))}))], 
                rep(c(". HyperTRIBE: ", ". DESeq2: "), times=4), 
                rep(names(ground_truth_lists_enriched_data)[1:4], each=2), 
                sep='')) 
dev.off()

################# CUMULATIVE PLOT ################ 
#### Prepare plot data ####
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

#### Plot cumulative distribution IMP2 ####
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
                  ecdf_df_list$`WT/mCherry vs. IMP2 - WT vs. mCherry`[ecdf_df_list$`WT/mCherry vs. IMP2 - WT vs. mCherry`$has_A2G == "IMP2-", "log2FoldChange"])
ks_test$p.value

plot = ecdf_list[[1]]
plot_ks = plot + 
    geom_text(label = paste('KS test p-value'), x = -20, y = 0.8)
plot_ks
