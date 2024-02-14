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

all_datasets = data.frame(dataset = c('GSE57957', 'GSE25097', 'GSE14520', 'GSE54236', "current"),
                          gene_symbol = c('Symbol', 'GeneSymbol', 'Gene Symbol', 'GENE_SYMBOL', NULL),
                          conds = c('T$:N$', 'tumor:healthy|non_tumor', 'A$:B$', 'tumor:non-malignant', NULL)
                          )
phenotypes = c('Tumor', 'Normal')
cor_methods = c("pearson", "spearman")
comparisons <- c('wt_imp2', 'mcherry_imp2', 'wt_mcherry', 'wt/mcherry-imp2')
comparison_names <- c('WT vs. IMP2', 'mCherry vs. IMP2', 'WT vs. mCherry', 'WT/mCherry vs. IMP2')

#### STEP 0: Prepare expression data ####
dataset_id = "GSE25097"
dataset_no = 1
dataset <- getGEO(dataset_id)
count_df <- exprs(dataset[[dataset_no]])

# Define gene_synonyms list
gene_ids <- unlist(dataset[[dataset_no]]@featureData@data[all_datasets$gene_symbol[all_datasets$dataset == dataset_id]], use.names = F)
gene_ids[nchar(gene_ids)==max(nchar(gene_ids))]

# Define gene_synonyms list
if (dataset_id == 'GSE57957'){
    gene_alias = mclapply(gene_ids, function(x){
        all_alias = alias2Symbol(x, species = 'Hs')
        all_alias = all_alias[!all_alias %in% x]
        all_alias = paste('_', x, "_", paste(all_alias, collapse='_'), "_", sep='')
        return(all_alias)
        })
    gene_synonyms <- dataset[[dataset_no]]@featureData@data$Synonyms
    gene_synonyms <- lapply(gene_synonyms, function(x) paste(paste(strsplit(x, '; ')[[1]], collapse='_'), '_', sep=''))
    gene_synonyms <- paste(gene_alias, gene_synonyms, sep='')
} else if(dataset_id == "GSE25097" | dataset_id == 'GSE54236'){
    gene_synonyms <- mclapply(unique(gene_ids), function(x){
        all_alias = alias2Symbol(x, species = 'Hs')
        all_alias = all_alias[!all_alias %in% x]
        all_alias = paste('_', paste(all_alias, collapse='_'), '_', sep='')
        return(all_alias)
    })
} else if(dataset_id == "GSE14520"){
    gene_synonyms <- mclapply(unique(gene_ids), function(x){
        syms = strsplit(str_replace_all(x, " ", ""), '/', fixed=T)[[1]]
        syms = syms[syms != ""]
        
        all_alias = alias2Symbol(syms[1], species = 'Hs')
        all_alias = all_alias[!all_alias %in% syms]
        all_alias = paste(all_alias, collapse='_')
        return(paste('_', paste(syms, collapse='_'), '_', all_alias, '_', sep=''))
    })
}


#### STEP 1: PREPARE HTRIBE LIST ####
step_1_result_list = list()
for (i in 1:(length(comparisons)-1)){
    comparison = comparisons[i]
    HTRIBE_list <- fread(paste("../HYPERTRIBE_result_CDS_UTR/", comparison, "_A2G_1%.xls", sep=''))
    HTRIBE_list_all = HTRIBE_list$Gene_name
    
    # Map human-mouse homologs
    human_mouse_map <- mclapply(HTRIBE_list_all, function(x) unique(gene_synonyms[grep(paste('_', x, '_', sep=''), gene_synonyms, ignore.case = T, fixed = FALSE)])[1])
    human_mouse_map <- unlist(lapply(human_mouse_map, function(x) strsplit(toString(x), '_', fixed = F)[[1]][2])) 
    human_mouse_map = data.frame('human'=human_mouse_map, 'mouse'=HTRIBE_list_all)
    exception_genes = human_mouse_map$mouse[is.na(human_mouse_map$human)]
    
    exception_genes_ENTREZ <- mapIds(org.Mm.eg.db, exception_genes, "ENTREZID", "SYMBOL")
    human_mouse_map_ENTREZ <- select(Orthology.eg.db, exception_genes_ENTREZ, "Homo_sapiens","Mus_musculus")
    temp <- AnnotationDbi::select(org.Hs.eg.db, keys=sapply(human_mouse_map_ENTREZ$Homo_sapiens, toString), columns=c("ENTREZID", "SYMBOL"), keytype="ENTREZID", multiVals="first")
    human_mouse_map_ENTREZ$human_homolog <-  temp$SYMBOL
    human_mouse_map_ENTREZ$mouse <- rownames(human_mouse_map_ENTREZ)
    
    human_mouse_map_homolog = human_mouse_map %>% left_join(human_mouse_map_ENTREZ[, c('mouse', 'human_homolog')])
    human_mouse_map_homolog$human = if_else(is.na(human_mouse_map_homolog$human), human_mouse_map_homolog$human_homolog, human_mouse_map_homolog$human)
    print(paste("Dim human_mouse_map_homolog BEFORE filter:", dim(human_mouse_map_homolog)))
    human_mouse_map_homolog = human_mouse_map_homolog[!is.na(human_mouse_map_homolog$human), c('human', 'mouse')]
    print(paste("Dim human_mouse_map_homolog AFTER filter:", dim(human_mouse_map_homolog)))
    
    human_mouse_map_homolog_subset = human_mouse_map_homolog[human_mouse_map_homolog$mouse %in% HTRIBE_list_all, ]
    # Define sample names and tag tumor/normal tissues
    count_df_modified = count_df
    rownames(count_df_modified) <- gene_ids
    conds = strsplit(all_datasets$conds[all_datasets$dataset == dataset_id], ':')[[1]]
    
    if (dataset_id == "GSE57957" ){
        colnames(count_df_modified) <- sample_ids
    } else if (dataset_id == 'GSE54236'){
        patient_id = sapply(dataset[[dataset_no]]$title, function(x) strsplit(x, '_')[[1]][2]) 
        colnames(count_df_modified) = paste(patient_id, colnames(count_df_modified), sample_ids, sep='_')
        paired_samples = names(table(patient_id))[table(patient_id) == 2]
        patient_id_paired = patient_id %in% paired_samples
        count_df_modified = count_df_modified[, patient_id_paired == TRUE]
    } else {
        colnames(count_df_modified) <- paste(colnames(count_df), sample_ids, sep='_')
    }
    
    tumor = colnames(count_df_modified)[grep(conds[1], colnames(count_df_modified), fixed=F)]
    normal = colnames(count_df_modified)[grep(conds[2], colnames(count_df_modified), fixed=F)]
    
    if (dataset_id == "GSE25097") tumor = tumor[grep('non', tumor, invert = T)]
    
    
    step_1_result_list[[i]] = list(human_mouse_map_homolog_subset, count_df_modified, gene_synonyms, tumor, normal)
    names(step_1_result_list[[i]]) = c('human_mouse_map_homolog_subset', 'count_df', 'gene_synonyms', 'tumor', 'normal') 
}
names(step_1_result_list) = comparison_names[1:(length(comparison_names)-1)]
saveRDS(step_1_result_list, file = paste('../analyzed_result/expression_correlation/', cor_method, '/step_1_result_list_', dataset_id, '_', dataset_no, '.RDS', sep=''))

length(step_1_result_list$`WT vs. IMP2`$tumor)
length(step_1_result_list$`WT vs. IMP2`$normal)

#### STEP 2: CORRELATION WITH IMP2 ####
all_datasets$dataset

for (dataset_id in all_datasets$dataset){
    dataset_length = 1
    if (dataset_id == "GSE14520"){
        dataset_length = 2
    }
    for (dataset_no in 1:dataset_length){
        print(paste(dataset_id, dataset_no, sep = ' '))
        step_1_result_list = readRDS(paste('../analyzed_result/expression_correlation/', cor_method, '/step_1_result_list_', dataset_id, '_', dataset_no, '.RDS', sep=''))
        step_2_result_list = list()
        
        for (i in 1:length(step_1_result_list)){
            human_mouse_map_homolog_subset = step_1_result_list[[i]]$human_mouse_map_homolog_subset
            count_df_modified = step_1_result_list[[i]]$count_df
            gene_synonyms = step_1_result_list[[i]]$gene_synonyms
            tumor = step_1_result_list[[i]]$tumor
            normal = step_1_result_list[[i]]$normal
            
            r_values_df = list() 
            for (j in 1:length(phenotypes)){
                phenotype = phenotypes[j]
                if (phenotype == "Tumor") 
                {sample_ids = tumor}
                else
                {sample_ids = normal}
                
                count_df_IMP2 = count_df_modified[grep('^IGF2BP2$|^IMP2$', rownames(count_df_modified), fixed = F), sample_ids]
                if (dataset_id == "GSE54236"){
                    count_df_IMP2 = count_df_IMP2[1, ]
                }
                count_df_others = count_df_modified[grep('^IGF2BP2$|^IMP2$', rownames(count_df_modified), fixed = F, invert = T), sample_ids]
                
                r_values = apply(count_df_others, 1, function(x){
                    return(cor.test(x, count_df_IMP2, method = cor_method)$estimate)
                })
                names(r_values) = make.names(rownames(count_df_others), unique=T)
                r_values = reshape2::melt(r_values, value.name = 'R-value', variable.name=rownames(r_values))
                
                r_values = r_values %>%
                    dplyr::mutate(
                        gene_id = rownames(count_df_others), # Not rownames(r_values) because it is already numbered
                        phenotype = phenotype,
                        # has_A2G_10 = if_else(gene_id %in% human_mouse_map_homolog_subset$human[1:10], '10', '.'),
                        has_A2G_50 = if_else(gene_id %in% human_mouse_map_homolog_subset$human[1:50], '50', '.'),
                        has_A2G_100 = if_else(gene_id %in% human_mouse_map_homolog_subset$human[1:100], '100', '.'),
                        has_A2G_500 = if_else(gene_id %in% human_mouse_map_homolog_subset$human[1:500], '500', '.'),
                        has_A2G_all = if_else(gene_id %in% human_mouse_map_homolog_subset$human, 'all', 'other')
                    )
                r_values_df[[j]] = r_values
            } 
            r_values_df = do.call(rbind, r_values_df)
            dim(r_values_df)
            r_values_df_melt = rbindlist(list( 
                r_values_df[r_values_df$has_A2G_50 == '50', c('R-value', 'gene_id', 'phenotype', 'has_A2G_50')],
                r_values_df[r_values_df$has_A2G_100 == '100', c('R-value', 'gene_id', 'phenotype', 'has_A2G_100')],
                r_values_df[r_values_df$has_A2G_500 == '500', c('R-value', 'gene_id', 'phenotype', 'has_A2G_500')],
                r_values_df[r_values_df$has_A2G_all == 'all', c('R-value', 'gene_id', 'phenotype', 'has_A2G_all')],
                r_values_df[r_values_df$has_A2G_all == 'other', c('R-value', 'gene_id', 'phenotype', 'has_A2G_all')]
            ), use.names=FALSE
            )
            colnames(r_values_df_melt) = c('R-value', 'gene_id', 'phenotype', 'type')
            step_2_result_list[[i]] = r_values_df_melt
        }
        names(step_2_result_list) = names(step_1_result_list)
        
        # Get WT/mCherry vs. IMP2 group
        control_vs_IMP2_genes = setdiff(
            union(step_2_result_list$`WT vs. IMP2`$gene_id[step_2_result_list$`WT vs. IMP2`$type == 'all'], 
                  step_2_result_list$`mCherry vs. IMP2`$gene_id[step_2_result_list$`mCherry vs. IMP2`$type == 'all']), 
            step_2_result_list$`WT vs. mCherry`$gene_id[step_2_result_list$`WT vs. mCherry`$type == 'all'])
        control_vs_IMP2_result = rbind(step_2_result_list$`WT vs. IMP2`, step_2_result_list$`mCherry vs. IMP2`)
        control_vs_IMP2_result = unique(control_vs_IMP2_result[control_vs_IMP2_result$gene_id %in% control_vs_IMP2_genes, ])
        step_2_result_list[[4]] = control_vs_IMP2_result
        names(step_2_result_list) = comparison_names
        
        saveRDS(step_2_result_list, file = paste("../analyzed_result/expression_correlation/", cor_method, "/step_2_result_list_", dataset_id, '_', dataset_no, '.RDS', sep=''))
    }
}
dim(step_2_result_list[[3]])
table(step_2_result_list[[3]]$gene_id[step_2_result_list[[3]]$type == 'other'] == "")
step_2_result_list[[4]]

#### STEP 3: PLOT R-VALUE ####
pretify_pval <- function(p_val){
    return(ifelse(round(p_val, 3) < 0.05, ": < 0.05", paste(": ", round(p_val, 3), sep='')))
}
plot_ecdf <- function(result_list, plot_title, plot_normal_only = TRUE, show_ks_test = TRUE){
    if (plot_normal_only){
        result_list = result_list[result_list$phenotype == "Normal", ]
        ecdf_ <- ggplot(result_list, aes(`R-value`, colour=type)) +
            stat_ecdf(aes(colour=type), n = 500) 
        
        if (show_ks_test){
            ks_test_all = ks.test(result_list[result_list$type == "all" , "R-value"][[1]], 
                              result_list[result_list$type == "other", "R-value"][[1]])
            p_val_all = pretify_pval(ks_test_all$p.value)
            
            ks_test_50 = ks.test(result_list[result_list$type == "100" , "R-value"][[1]], 
                              result_list[result_list$type == "other", "R-value"][[1]])
            p_val_50 = pretify_pval(ks_test_50$p.value) 
            
            ecdf_ <- ecdf_ + 
                geom_text(label = paste('All A2G genes', p_val_all, '\nTop 100 A2G genes', p_val_50, sep=""), 
                          aes(colour=NULL, linetype=NULL), x = 0.5, y = 0.1, show.legend = F, size=3)
            }
        }
    else { 
        ecdf_ <- ggplot(result_list, aes(`R-value`, colour=type, linetype=phenotype)) +
            stat_ecdf(aes(colour=type, linetype=phenotype), n = 500) +
            scale_linetype_discrete(name = "Phenotype")
        
        if (show_ks_test){
            result_list_normal = result_list[result_list$phenotype == "Normal", ]
            result_list_tumor = result_list[result_list$phenotype == "Tumor", ]
            ks_test_normal = ks.test(result_list_normal[result_list_normal$type != "other", "R-value"][[1]], 
                                     result_list_normal[result_list_normal$type == "other", "R-value"][[1]])
            ks_test_tumor = ks.test(result_list_tumor[result_list_tumor$type != "other", "R-value"][[1]], 
                                    result_list_tumor[result_list_tumort$type == "other", "R-value"][[1]])
            p_val_normal = pretify_pval(ks_test_normal$p.value) 
            p_val_tumor = pretify_pval(ks_test_normal$p.value) 
            
            ecdf_ <- ecdf_ +  geom_text(
                label = paste(
                    'p-value (Normal)', p_val_normal, '\np-value (Tumor)', p_val_tumor, 
                    sep=""),  
                aes(colour=NULL, linetype=NULL), x = 0.6, y = 0.08, show.legend = F)
            }
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
ecdf_plotname = c()
for (dataset_id in all_datasets$dataset){
    dataset_length = 1 
    if (dataset_id == "GSE14520"){
        dataset_length = 2
    }
    for (i in 1:dataset_length){
        step_2_result_list = readRDS(paste("../analyzed_result/expression_correlation/", cor_method, "/step_2_result_list_", dataset_id, '_', i, '.RDS', sep=''))
        
        for (j in 1:length(step_2_result_list)){
            plot_data = step_2_result_list[[j]]
            plot_data = plot_data[plot_data$type != "50" & plot_data$gene_id != "", ]
            ecdf_list = c(ecdf_list, list(plot_ecdf(plot_data, comparison_names[j], plot_normal_only = TRUE, show_ks_test = TRUE)))
        }
        plotname <- paste(c(dataset_id, i), collapse="_")
        ecdf_plotname = c(ecdf_plotname, c(plotname))
    }
}

ecdf_plotname <- sapply(ecdf_plotname, function(x) {
    plotname = strsplit(x, "_")[[1]]
    if (plotname[1] == "GSE14520"){
        return(paste(plotname[1], plotname[2], sep=' - Cohort '))
    } else {
        return(plotname[1])
    }
})
names(ecdf_list) = rep(ecdf_plotname, each=4)

png(paste("../analyzed_result/expression_correlation/", cor_method, "/ECDF_expression_combined_100.png", sep=''), width = 14, height = 13, unit = 'in', res = 200)
ggarrange(plotlist = ecdf_list[grep("GSE57957|GSE54236|GSE14520.", names(ecdf_list))], ncol=4, nrow=4, 
          align = 'hv', 
          label.x = rep(c(0, -0.15, -0.15, 0), each=3),
          common.legend = TRUE,
          labels = unlist(lapply(paste(LETTERS[1:4], unique(names(ecdf_list)[grep("GSE57957|GSE54236|GSE14520.", names(ecdf_list))]), sep=". "), function(x) return(c(x, "", "", ""))))
          )
dev.off()

