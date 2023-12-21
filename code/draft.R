library(ggplot2)
library(dplyr)
library(reshape2)
setwd("~/Documents/BIOINFO/TRIBE/HTRIBE_analysis/1%_bed_motif_meme/_all/jf/")

experiments = c('wt_imp2', 'mcherry_imp2', 'wt_mcherry')
kmers = seq(5, 8, 1)

all_files = list.files(getwd(), pattern = 'histo.txt')
all_histos = lapply(all_files, function(x) read.table(x, quote="\"", comment.char=""))
info = lapply(all_files, function(x) strsplit(x, '_')[[1]])
for (i in 1:12){
    colnames(all_histos[[i]]) = c('freq', 'count')
    all_histos[[i]]$experiment = paste(info[[i]][1], info[[i]][2], sep='_')
    all_histos[[i]]$kmer = info[[i]][3]
}
histo_df = do.call('rbind', all_histos)
histo_df$total = histo_df$freq * histo_df$count
histo_df = histo_df %>%
    dplyr::group_by(experiment, kmer) %>%
    mutate(no_motives = sum(total), 
           no_distinct_motives = sum(count),
           label = paste(kmer, '-mers (', no_distinct_motives, ' distinct / ', no_motives, ' total)', sep=''),
           sublabel = paste(kmer, experiment, sep='_'))

levels = expand.grid(kmers, experiments)
colnames(levels) = c('kmer', 'experiment')
levels$factor = seq(1, nrow(levels), 1)
levels$kmer = as.character(levels$kmer)
histo_df = dplyr::inner_join(histo_df, levels, by=c('kmer', 'experiment'))
levels = unique(histo_df$sublabel[order(histo_df$factor)])
labels = unique(histo_df$label)
names(labels) = unique(histo_df$sublabel)

agact_count = data.frame(
    kmer=5,
    experiment=experiments, 
    freq=c(4332, 4829, 2782), 
    count=c(2.5, 4, 2.75),
    text=rep('AGACT', 3)
    )
agact_count = merge(agact_count, unique(histo_df[, c('experiment', 'kmer', 'label', 'sublabel')]))


png('png', units = 'in', width = 10, height = 7, res=200)
ggplot(data = histo_df, aes(x=freq, y=count, fill=experiment, col=experiment)) +
    geom_col(show.legend = F) +
    facet_wrap(vars(sublabel), labeller = as_labeller(labels),
               ncol=4, nrow=3, scales = 'free', dir='v') +
    # facet_wrap(vars(experiment, kmer), 
    #            ncol=4, nrow=3, scales = 'free') +
    geom_text(data=agact_count, 
              aes(x=freq+1000, y=count, label=paste(text, '(', freq, ')', sep='')), 
              show.legend = F) +
    geom_segment(data=agact_count, 
                 aes(x=freq, xend=freq+1000, y=2, yend=count-0.25), 
                 show.legend = F, linetype='dashed') +
    xlab('Frequency') + ylab('Count') + 
    scale_color_discrete(name='Comparison', labels=c('mCherry vs. IMP2', 'WT vs. IMP2', 'WT vs. mCherry')) +
    # theme_light() +
    theme(
        legend.position = 'bottom'
    )
dev.off()


