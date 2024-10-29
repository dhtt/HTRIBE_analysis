args<-commandArgs(TRUE)
library(ggplot2)
library(ggpubr)
library(dplyr)

motif_df = read.table(args[1], sep = '\t', header = T)
motif_df = motif_df %>%
    dplyr::mutate(pvalue = as.numeric(pvalue)) %>%
    arrange(pvalue)
print(head(motif_df))

all_region_plots = list()
for (i in 1:length(unique(motif_df$region))){
    region = unique(motif_df$region)[i]

    motif_df_subset = motif_df[motif_df$region == region, ]
    plot <- ggplot(data = motif_df_subset, aes(x = reorder(alt_name, pvalue), y = -log(pvalue))) +
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5)) +
        labs(title = paste(c("Top 10 motifs in", comparison, "binding sites"), collapse = ' '), x = "Motif", y = "-log(pvalue)")
    all_region_plots[[i]] = plot
}

result_path = dirname(args[1])
png(paste(c(result_path, 'motif.png'), collapse='/'), width = 8, height = 6.5, unit = 'in', res = 200)
ggarrange(plotlist = all_region_plots, ncol = 1, nrow = length(all_region_plots))
dev.off()

