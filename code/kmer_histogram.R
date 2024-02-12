library(ggplot2)
library(dplyr)
library(reshape2)
setwd("~/Documents/BIOINFO/TRIBE/HTRIBE_analysis/1%_bed_motif_meme/_all/jf/")

experiments = c('wt_imp2', 'mcherry_imp2', 'wt_mcherry', 'baseline_mm39')
kmers = seq(5, 8, 1)

all_files = list.files(getwd(), pattern = 'histo.txt')
all_histos = lapply(all_files, function(x) read.table(x, quote="\"", comment.char=""))
info = lapply(all_files, function(x) strsplit(x, '_')[[1]])
for (i in 1:length(info)){
    colnames(all_histos[[i]]) = c('freq', 'count')
    all_histos[[i]]$experiment = paste(info[[i]][1], info[[i]][2], sep='_')
    all_histos[[i]]$kmer = info[[i]][3]
}

histo_df = do.call('rbind', all_histos)
histo_df = histo_df[histo_df$kmer %in% kmers, ]
histo_df$total = histo_df$freq * histo_df$count
histo_df = histo_df %>%
    dplyr::group_by(experiment, kmer) %>%
    mutate(no_motives = sum(total), 
           no_distinct_motives = sum(count),
           label = paste(kmer, '-mers (', no_distinct_motives*2, ' distinct / ', no_motives*2, ' total)', sep=''),
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
    freq=c(4332, 4829, 2782, 6352036), 
    count=c(2.5, 4, 2.75, 25),
    text=rep('AGACT', 4)
    )
agact_count = merge(agact_count, unique(histo_df[, c('experiment', 'kmer', 'label', 'sublabel')]))

png('motives_freq_histogram.png', units = 'in', width = 12, height = 9, res=200)
ggplot(data = histo_df, aes(x=freq, y=count, fill=experiment, col=experiment)) +
    geom_col() +
    facet_wrap(vars(sublabel), labeller = as_labeller(labels),
               ncol=4, nrow=4, scales = 'free', dir='v') +
    # facet_wrap(vars(experiment, kmer), 
    #            ncol=4, nrow=3, scales = 'free') +
    geom_text(data=agact_count, 
              aes(x=freq+1000, y=count, label=paste(text, '(', freq, ')', sep='')), 
              show.legend = F) +
    geom_segment(data=agact_count, 
                 aes(x=freq, xend=freq+1000, y=2, yend=count-0.25), 
                 show.legend = F, linetype='dashed') +
    xlab('Frequency') + ylab('Count') + 
    scale_color_discrete(guide = "none") +
    scale_fill_discrete(name='Comparison', labels=c('Reference genome', 'mCherry vs. IMP2', 'WT vs. IMP2', 'WT vs. mCherry')) +
    # theme_light() +
    theme(
        legend.position = 'bottom' 
    )
dev.off()


all_files = list.files(getwd(), pattern = 'dump.txt')
info = lapply(all_files, function(x) strsplit(x, '_')[[1]])
count_df = list()
for (i in 1:length(info)){
    count_file = read.table(all_files[i], quote = '"')
    x_df = data.frame(
        count = count_file[seq(2, nrow(count_file), 2), ],
        motif = count_file[seq(1, nrow(count_file), 2), ])
    
    colnames(x_df) = c('motif', 'count')
    x_df$experiment = paste(info[[i]][1], info[[i]][2], sep='_')
    x_df$kmer = info[[i]][3]
    x_df$count = gsub('>', '', x_df$count, fixed=TRUE)
    count_df[[i]] = x_df
}
count_df = do.call('rbind', count_df)
count_df = count_df[count_df$kmer %in% kmers, ]

exp_group = experiments[1:3]
control_group = experiments[3]
control_group = experiments[3:4]
motives = c('AGACT', head(unique(count_df$motif), 100))
res = list()
res_1 = list()
res_2 = list()

for (i in 1:length(exp_group)){
    for (j in 1:length(control_group)){
        exp_id = exp_group[i]
        ctl_id = control_group[j]
        if (exp_id != ctl_id){
            res = lapply(motives, function(motif){
                exp_true = as.numeric(count_df$count[count_df$motif == motif & count_df$experiment == exp_id])
                ctl_true = as.numeric(count_df$count[count_df$motif == motif & count_df$experiment == ctl_id])
                exp_false = sum(as.numeric(count_df$count[count_df$motif != motif & count_df$experiment == exp_id]))
                ctl_false = sum(as.numeric(count_df$count[count_df$motif != motif & count_df$experiment == ctl_id]))
                prevalence_ratio = (exp_true/exp_false)/(ctl_true/ctl_false)
                
                group_name = paste(c(exp_id, ctl_id, motif), collapse=' + ')
                return(c(exp_id, ctl_id, motif, exp_true, exp_false, ctl_true, ctl_false, prevalence_ratio))
            })
        }
        res_1[[j]] = res
    }
    res_2[[i]] = res_1
}
prevalence_ratio = unlist(unlist(res_2, recursive = F), recursive = F)
prevalence_ratio = do.call(rbind, prevalence_ratio)
colnames(prevalence_ratio) = c('exp_id', 'ctl_id', 'motif', 'exp_true', 'exp_false', 'ctl_true', 'ctl_false', 'prevalence_ratio')
write.table(format(prevalence_ratio, decimal.mark = ','), file='prevalence_ratio.xls', quote = F, sep='\t', row.names = F, col.names = T)
