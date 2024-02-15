library(GEOquery)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(parallel)
library(data.table)
library(limma)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(Orthology.eg.db)
library(ggpubr) 

setwd("~/Documents/BIOINFO/TRIBE/HTRIBE_analysis/code")
phenotypes = c('WT', 'mCherry', "IMP2")
cor_method = "spearman"
comparisons <- c('wt_imp2', 'mcherry_imp2', 'wt_mcherry', 'wt/mcherry-imp2')
comparison_names <- c('WT vs. IMP2', 'mCherry vs. IMP2', 'WT vs. mCherry', 'WT/mCherry vs. IMP2')


#### STEP 1: PREPARE EXPRESSION DATA FROM CURRENT STUDY ####
dataset_id = "current"
dataset_no = 1
count_df = as.data.frame(t(fread("../analyzed_result/expression_data.csv", sep = '\t', quote = F)))
colnames(count_df) = count_df[1, ]
count_df = count_df[2:nrow(count_df), ]

gene_ids = rownames(count_df)
gene_ids = sapply(gene_ids, function(x) strsplit(x, '.', fixed = T)[[1]][1]) 

gene_ids_map <- as.data.table(AnnotationDbi::select(org.Mm.eg.db, keys = gene_ids, columns = c("REFSEQ", "SYMBOL"), keytype = "REFSEQ"))
count_df$gene_id = gene_ids_map$SYMBOL

#### STEP 2: CORRELATION WITH IMP2 ####
samples = seq(1, ncol(count_df)-1, 1)
samples = colnames(count_df)[grep("WT|mCherry", colnames(count_df))]
correlation_result_list = list()
for (i in 1:(length(comparisons)-1)){
    comparison = comparisons[i]
    HTRIBE_list <- fread(paste("../HYPERTRIBE_result_CDS_UTR/", comparison, "_A2G_1%.xls", sep=''))
    HTRIBE_list_all = HTRIBE_list$Gene_name
    
    
    count_df_IMP2 = count_df[grep('^IGF2BP2$|^IMP2$', gene_ids_map$SYMBOL, fixed = F, ignore.case = T), samples]
    count_df_others = count_df[, samples]
    count_df_others_gene_id = count_df[, 'gene_id']

    r_values = apply(count_df_others, 1, function(x){
        return(cor.test(as.numeric(x), as.numeric(count_df_IMP2), method = cor_method)$estimate)
    })
    r_values = reshape2::melt(r_values, value.name = 'R-value', variable.name=rownames(r_values))

    r_values = r_values %>%
        dplyr::mutate(
            gene_id = count_df_others_gene_id,
            # has_A2G_10 = if_else(gene_id %in% HTRIBE_list_all[1:10], '10', '.'),
            has_A2G_50 = if_else(gene_id %in% HTRIBE_list_all[1:50], '50', '.'),
            has_A2G_100 = if_else(gene_id %in% HTRIBE_list_all[1:100], '100', '.'),
            has_A2G_500 = if_else(gene_id %in% HTRIBE_list_all[1:500], '500', '.'),
            has_A2G_all = if_else(gene_id %in% HTRIBE_list_all, 'all', 'other')
        )

    r_values_df_melt = rbindlist(list(
        # r_values[r_values$has_A2G_10 == '10', c('R-value', 'gene_id', 'has_A2G_10')],
        r_values[r_values$has_A2G_50 == '50', c('R-value', 'gene_id', 'has_A2G_50')],
        r_values[r_values$has_A2G_100 == '100', c('R-value', 'gene_id', 'has_A2G_100')],
        r_values[r_values$has_A2G_500 == '500', c('R-value', 'gene_id', 'has_A2G_500')],
        r_values[r_values$has_A2G_all == 'all', c('R-value', 'gene_id', 'has_A2G_all')],
        r_values[r_values$has_A2G_all == 'other', c('R-value', 'gene_id', 'has_A2G_all')]
    ), use.names = FALSE
    )
    colnames(r_values_df_melt) = c('R-value', 'gene_id', 'type')
    correlation_result_list[[i]] = r_values_df_melt
}
names(correlation_result_list) = comparison_names[1:3]

control_vs_IMP2_genes = setdiff(
    union(unique(correlation_result_list$`WT vs. IMP2`$gene_id[correlation_result_list$`WT vs. IMP2`$type != 'other']),
          unique(correlation_result_list$`mCherry vs. IMP2`$gene_id[correlation_result_list$`mCherry vs. IMP2`$type != 'other'])),
    unique(correlation_result_list$`WT vs. mCherry`$gene_id[correlation_result_list$`WT vs. mCherry`$type != 'other'])
)
control_vs_IMP2_result = rbind(correlation_result_list$`WT vs. IMP2`, correlation_result_list$`mCherry vs. IMP2`)
control_vs_IMP2_result = unique(control_vs_IMP2_result[control_vs_IMP2_result$gene_id %in% control_vs_IMP2_genes, ])
correlation_result_list[[4]] = control_vs_IMP2_result
names(correlation_result_list) = comparison_names

saveRDS(correlation_result_list, file = paste("../analyzed_result/expression_correlation/", cor_method, "/result_list_", dataset_id, '_WTmCherry.RDS', sep=''))
# correlation_result_list = readRDS(paste("../analyzed_result/expression_correlation/", cor_method, "/result_list_", dataset_id, '.RDS', sep=''))


### STEP 3: PLOT R-VALUE ECDF ####
htribe_deseq2_genes = list(
    strsplit(c("Calu, Ccdc90b, Cdadc1, Clta, Cped1, Cplane1, Creb1, Crem, Cyp2c38, Dck, Ddx11, Ddx58, Dhx58, Dlg1, Ecm1, Elmod3, Esam, Fktn, Gm20604, Grk6, Hck, Igtp, Macf1, Mad2l2, Me2, Met, Mff, Mmrn2, Mtif3, Mug1, Numb, Parp14, Pecam1, Plec, Ppil3, Rabep1, Reln, Retreg1, Rnf38, Slc66a2, Stx5a, Syt12, Taf6, Tardbp, Tdrd7, Tgtp1, Tmem62, Tnpo1, Trim39, Vps39, Wnk1"), ', ')[[1]],
    strsplit(c('Creb1, Crem, Fam219a, Fktn, Gm20604, Map7d1, Med16, Nfix, Prpf40b, Ptbp3, Rabep1, Smim14, Taf6, Tardbp, Tnpo1, Vps39'), ', ')[[1]],
    strsplit(c('Add3, Cdadc1, Cyp27a1, Dhx58, Dlg1, Dusp12, Fktn, Gm11837, Hsd3b5, Lrp11, Mmrn2, Oasl1, Pigw, Plec, Pnpt1, Prpf40b, Retreg1, Rnf38, Slfn4, Tdrd7, Ttbk2, Ubr5, Wnk1, Zfp708'), ', ')[[1]],
    strsplit(c("Agfg1, Ccdc90b, Creb1, Cyp2c38, Dck, Mad2l2, Med16, Rabep1, Taf6, Tardbp, Tmem62, Tnpo1, Trim39"), ', ')[[1]])

plot_data_list = list()
for (i in 1:length(correlation_result_list)){
    plot_data_list[[i]] = correlation_result_list[[i]] %>%
        dplyr::filter(type != '50' & gene_id != "") %>%
        dplyr::mutate(
                      comparison = paste(LETTERS[i], names(correlation_result_list)[i], sep='. '),
                      is_deseq = ifelse(gene_id %in% htribe_deseq2_genes[[i]], T, F)
                      ) %>%
        dplyr::filter(!(type == 'other' & is_deseq == T))
}
plot_data_list = do.call('rbind', plot_data_list)


pretify_pval <- function(p_val){
    return(ifelse(round(p_val, 3) < 0.05, "< 0.05", round(p_val, 3)))
}

plot_ecdf <- function(result_list, plot_title){
    ecdf <- ggplot(result_list, aes(`R-value`, colour=type)) +
        stat_ecdf(aes(colour=type, linetype=is_deseq), n = 100, alpha = 0.8)  +
        facet_wrap(vars(result_list$comparison)) +
        xlab("R-values") + ylab("Probability") +
        scale_color_discrete(name = "Top A2G genes", label=c('100', '500', 'All', 'Other')) +
        scale_linetype_discrete(name = "A2G genes with altered expression", label=c('False', 'True')) +
        xlim(c(-1, 1)) +
        theme_minimal() +
        theme(
            aspect.ratio = 2/3,
            strip.text = element_text(size = 12, hjust = 0.5, face = "bold"),
            axis.text.x = element_text(colour="black", size = 12),
            axis.text.y = element_text(colour="black", size = 12),
            axis.title = element_text(face = "bold"),
            panel.border = element_rect(colour = "grey", fill = NA, size = 0.5),
            legend.position = "bottom",
            legend.title = element_text(face='bold')) +
        guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
    
    return(ecdf)
}
ecdf = plot_ecdf(plot_data_list, "")

png(paste("../analyzed_result/expression_correlation/", cor_method, "/ECDF_expression_current_WTmCherry_withdeseq2.png", sep=''), width = 7, height = 6, unit = 'in', res = 200)
ecdf
dev.off()



