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

all_datasets = data.frame(dataset = c('GSE57957', 'GSE25097', 'GSE14520', 'GSE54236'),
                          gene_symbol = c('Symbol', 'GeneSymbol', 'Gene Symbol', 'GENE_SYMBOL'),
                          conds = c('T$:N$', 'tumor:healthy|non_tumor', 'A$:B$', 'tumor:non-malignant')
                          )
phenotypes = c('Tumor', 'Normal')
cor_method = "spearman"
comparisons <- c('wt_imp2', 'mcherry_imp2', 'wt_mcherry')
comparison_names <- c('WT vs. IMP2', 'mCherry vs. IMP2', 'WT vs. mCherry')

#### STEP 0: Prepare expression data ####
dataset_id = "GSE25097"
dataset_no = 1
dataset <- getGEO(dataset_id)
count_df <- exprs(dataset[[dataset_no]])

# head( dataset[[dataset_no]]@featureData@GENE_SYMBOL)
sample_ids <- dataset[[dataset_no]]$source_name_ch1
organism <- dataset[[dataset_no]]$taxid_ch1

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
for (i in 1:length(comparisons)){
    
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
names(step_1_result_list) = comparison_names
saveRDS(step_1_result_list, file = paste('../analyzed_result/expression_correlation/', cor_method, '/step_1_result_list_', dataset_id, '_', dataset_no, '.RDS', sep=''))
step_1_result_list = readRDS(paste('../analyzed_result/expression_correlation/', cor_method, '/step_1_result_list_', dataset_id, '_', dataset_no, '.RDS', sep=''))
length(step_1_result_list$`WT vs. IMP2`$tumor)
length(step_1_result_list$`WT vs. IMP2`$normal)

#### STEP 2: CORRELATION WITH IMP2 ####
step_2_result_list = list()
for (i in 1:length(comparisons)){
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
                has_A2G_10 = if_else(gene_id %in% human_mouse_map_homolog_subset$human[1:10], '10', '.'),
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
        # r_values_df[r_values_df$has_A2G_10 == '10', c('R-value', 'gene_id', 'phenotype', 'has_A2G_10')],
        r_values_df[r_values_df$has_A2G_50 == '50', c('R-value', 'gene_id', 'phenotype', 'has_A2G_50')],
        # r_values_df[r_values_df$has_A2G_100 == '100', c('R-value', 'gene_id', 'phenotype', 'has_A2G_100')],
        r_values_df[r_values_df$has_A2G_500 == '500', c('R-value', 'gene_id', 'phenotype', 'has_A2G_500')],
        r_values_df[r_values_df$has_A2G_all == 'all', c('R-value', 'gene_id', 'phenotype', 'has_A2G_all')],
        r_values_df[r_values_df$has_A2G_all == 'other', c('R-value', 'gene_id', 'phenotype', 'has_A2G_all')]
    ), use.names=FALSE
    )
    colnames(r_values_df_melt) = c('R-value', 'gene_id', 'phenotype', 'type')
    step_2_result_list[[i]] = r_values_df_melt
}
names(step_2_result_list) = names(step_1_result_list)
saveRDS(step_2_result_list, file = paste("../analyzed_result/expression_correlation/", cor_method,"/step_2_result_list_", dataset_id, '_', dataset_no, '.RDS', sep=''))

plot_ecdf <- function(result_list, plot_title){
    ecdf <- ggplot(result_list, aes(`R-value`, colour=type, linetype=phenotype)) + 
        stat_ecdf(aes(colour=type, linetype=phenotype), n = 500) +
        xlab("R-values") + ylab("Probability") +
        scale_color_discrete(name = "Top A2G genes") +
        scale_linetype_discrete(name = "Phenotype") +
        xlim(c(-1, 1)) +
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
for (i in 1:length(step_2_result_list)){
    ecdf_list[[i]] = plot_ecdf(step_2_result_list[[i]], names(step_2_result_list)[i])
}

png(paste("../analyzed_result/expression_correlation/", cor_method,"/ECDF_expression_", dataset_id, '_', dataset_no, '.png', sep=''), width = 12, height = 4, unit = 'in', res = 200)
ggarrange(plotlist = ecdf_list, ncol=3, nrow=1, common.legend = TRUE)
dev.off()

