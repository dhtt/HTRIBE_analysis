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
correlation_result_list = list()
for (i in 1:(length(comparisons)-1)){
    comparison = comparisons[i]
    HTRIBE_list <- fread(paste("../HYPERTRIBE_result_CDS_UTR/", comparison, "_A2G_1%.xls", sep=''))
    HTRIBE_list_all = HTRIBE_list$Gene_name
    
    count_df_IMP2 = count_df[grep('^IGF2BP2$|^IMP2$', gene_ids_map$SYMBOL, fixed = F, ignore.case = T), 1:(ncol(count_df)-1)]
    count_df_others = count_df[, 1:(ncol(count_df)-1)]
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
    ), use.names=FALSE
    )
    colnames(r_values_df_melt) = c('R-value', 'gene_id', 'type')
    correlation_result_list[[i]] = r_values_df_melt
}
names(correlation_result_list) = comparison_names[1:3]
# saveRDS(correlation_result_list, file = paste("../analyzed_result/expression_correlation/", cor_method, "/result_list_", dataset_id, '.RDS', sep=''))
# correlation_result_list = readRDS(paste("../analyzed_result/expression_correlation/", cor_method, "/result_list_", dataset_id, '.RDS', sep=''))

control_vs_IMP2_genes = setdiff(
    union(unique(correlation_result_list$`WT vs. IMP2`$gene_id[correlation_result_list$`WT vs. IMP2`$type != 'other']),
          unique(correlation_result_list$`mCherry vs. IMP2`$gene_id[correlation_result_list$`mCherry vs. IMP2`$type != 'other'])),
    unique(correlation_result_list$`WT vs. mCherry`$gene_id[correlation_result_list$`WT vs. mCherry`$type != 'other'])
)
control_vs_IMP2_result = rbind(correlation_result_list$`WT vs. IMP2`, correlation_result_list$`mCherry vs. IMP2`)
control_vs_IMP2_result = unique(control_vs_IMP2_result[control_vs_IMP2_result$gene_id %in% control_vs_IMP2, ])
correlation_result_list[[4]] = control_vs_IMP2_result
names(correlation_result_list) = comparison_names


#### STEP 3: PLOT R-VALUE ECDF ####
pretify_pval <- function(p_val){
    return(ifelse(round(p_val, 3) < 0.05, ": < 0.05", paste(": ", round(p_val, 3), sep='')))
}

plot_ecdf <- function(result_list, plot_title, show_ks_test = TRUE){
    ecdf_ <- ggplot(result_list, aes(`R-value`, colour=type)) +
        stat_ecdf(aes(colour=type), n = 100, alpha = 0.8) 
    
    if (show_ks_test){
        print(dim(result_list[result_list$type == "100", ]))
        ks_test_all = ks.test(result_list[result_list$type == "all" , "R-value"][[1]], 
                              result_list[result_list$type == "other", "R-value"][[1]])
        p_val_all = pretify_pval(ks_test_all$p.value)
        ks_test_50 = ks.test(result_list[result_list$type == "100" , "R-value"][[1]], 
                             result_list[result_list$type == "other", "R-value"][[1]])
        p_val_50 = pretify_pval(ks_test_50$p.value) 
        
        ecdf_ <- ecdf_ +  
            geom_text(label = paste('All A2G genes ', p_val_all, '\nTop 100 A2G genes', p_val_50, sep=""), 
                      aes(colour=NULL, linetype=NULL), x = 0.5, y = 0.1, show.legend = F, size=3)  
           
    }
    
    ecdf <- ecdf_ + xlab("R-values") + ylab("Probability") +
        scale_color_discrete(name = "Top A2G genes") +
        xlim(c(-1, 1)) +
        theme_minimal() +
        ggtitle(plot_title) +
        theme(
            aspect.ratio = 2/3,
            plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
            axis.text.x = element_text(colour="black", size = 12),
            axis.text.y = element_text(colour="black", size = 12),
            axis.title = element_text(face = "bold"),
            panel.border = element_rect(colour = "grey", fill = NA, size = 0.5),
            legend.position = "bottom") +
        guides(fill = guide_legend(nrow = 3, byrow = t))
    return(ecdf)
}
ecdf_list = c()
for (i in 1:length(correlation_result_list)){
    plot_data = correlation_result_list[[i]]
    plot_data = plot_data[plot_data$type != "50" & plot_data$gene_id != "", ]
    ecdf_list = c(ecdf_list, list(plot_ecdf(plot_data, names(correlation_result_list)[i], show_ks_test = TRUE)))
}
names(ecdf_list) = names(correlation_result_list)

png(paste("../analyzed_result/expression_correlation/", cor_method, "/ECDF_expression_current_100.png", sep=''), width = 7, height = 6, unit = 'in', res = 200)
ggarrange(plotlist = ecdf_list, ncol=2, nrow=2, 
          align = 'hv', common.legend = TRUE)
dev.off()



