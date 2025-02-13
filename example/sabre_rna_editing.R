library(ggplot2)
library(dplyr)
library(stringr)
library(reshape2)
library(ggprism)

gtex.edQTL <- read.delim(gzfile('~/Downloads/dataset/edQTL/Whole-Blood.signif_variant_site_pairs.txt.gz'))
gtex.edQTL$variant <- gsub('_[A-Z]*_[A-Z]*_b38$', '', gtex.edQTL$variant_id)
gtex <- read.delim(gzfile('~/Downloads/dataset/edQTL/Whole-Blood.edsite.txt.gz'))
length(table(gtex$gene_id))
length(table(gtex$variant_id))


sample_list <- read.table('~/01_large_files/faser/list_sample.txt')
sample_list$batch <- substr(sample_list$V1, 1, 11)

dirs <- list.dirs('~/01_large_files/faser/output/', full.names = F, recursive = F)
sample_list.dirs <- subset(sample_list, V2 %in% dirs)
sort(table(sample_list.dirs$batch), decreasing = T)



# batch.tmp <- 'SG_HEL_B001'
# sample_list.tmp <- subset(sample_list.dirs, batch == batch.tmp)
sample_list$V2 <- paste0(sample_list$V2, '_', substr(sample_list$V1, 13, 16))
sample_list$V2 <- gsub('_L001$', '', sample_list$V2)

sample_list.tmp <- sample_list
df_all_donors <- data.frame()
for (j in 1:nrow(sample_list.tmp)) {
  if (length(list.files(path = paste0('~/01_large_files/faser/output/', sample_list.tmp$V2[j]), pattern = 'cell_allele_connections_chr' )) == 22) {
    
    barcode <- read.table(paste0('~/01_large_files/faser/meta/meta_', sample_list.tmp$V1[j],'.txt'))
    
    data <- data.frame()
    for (i in seq(22)) {
      data.tmp <- read.delim(paste0('~/01_large_files/faser/output/', sample_list.tmp$V2[j], '/cell_allele_connections_chr', i, '.txt'), colClasses = c(rep('character', 3), rep('double', 1)))
      data <- rbind(data, data.tmp)
    }
    data.ct <- merge(data, barcode, by.x = 'barcode', by.y = 'V1')
    pos1 <- gsub('chr[0-9]*_(.*)_\\._.*', '\\1', str_split_fixed(data.ct$var, ',', 2)[, 1])
    pos2 <- gsub('chr[0-9]*_(.*)_\\._.*', '\\1', str_split_fixed(data.ct$var, ',', 2)[, 2])
    data.ct$var[which(!pos1 < pos2)] <- paste(str_split_fixed(data.ct$var, ',', 2)[, 2], str_split_fixed(data.ct$var, ',', 2)[, 1], sep = ',')[which(!pos1 < pos2)]
    
    data.ct$var1 <- str_split_fixed(data.ct$var, ',', 2)[, 1]
    data.ct$var2 <- str_split_fixed(data.ct$var, ',', 2)[, 2]
    data.ct$g1 <- substr(data.ct$geno, 1, 1)
    data.ct$g2 <- substr(data.ct$geno, 2, 2)
    data.ct$ref1 <- str_split_fixed(str_split_fixed(data.ct$var1, '_\\._', 2)[, 2], '_', 2)[, 1]
    data.ct$alt1 <- str_split_fixed(str_split_fixed(data.ct$var1, '_\\._', 2)[, 2], '_', 2)[, 2]
    data.ct$ref2 <- str_split_fixed(str_split_fixed(data.ct$var2, '_\\._', 2)[, 2], '_', 2)[, 1]
    data.ct$alt2 <- str_split_fixed(str_split_fixed(data.ct$var2, '_\\._', 2)[, 2], '_', 2)[, 2]
    data.ct$v1 <- gsub('_\\._.[A-Z]*_[A-Z]*', '', data.ct$var1)
    data.ct$v2 <- gsub('_\\._.[A-Z]*_[A-Z]*', '', data.ct$var2)
    
    data.ct$v1_editSite <- data.ct$v1 %in% gtex$gene_id
    data.ct$v2_editSite <- data.ct$v2 %in% gtex$gene_id
    
    data.ct$v1_SNP <- data.ct$v1 %in% gtex$variant
    data.ct$v2_SNP <- data.ct$v2 %in% gtex$variant
    
    
    data.ct.edit <- subset(data.ct, v1_editSite == T | v2_editSite == T)
    # data.ct.edQTL <- subset(data.ct.edit, v1_editSite == T & v2_SNP == T | v2_editSite == T & v1_SNP == T)
    
    data.ct.edit$donor <- sample_list.tmp$V2[j]
    # data.ct.edQTL$donor <- sample_list.tmp$V2[j]
    # df_all_donors <- rbind(df_all_donors, data.ct.edQTL)
    df_all_donors <- rbind(df_all_donors, data.ct.edit)
  }
  else{print("not complete")}
}
table(df_all_donors$donor)
length(table(df_all_donors$donor))


# cell types aggregated 
df_all_donors.summary <- df_all_donors %>% group_by(var, geno, V2, v1_editSite, v2_editSite, v1_SNP, v2_SNP) %>% summarise(support = sum(support))
df_all_donors.summary.dcast <- dcast(df_all_donors.summary, var+V2+v1_editSite+v2_editSite+v1_SNP+v2_SNP ~ geno, value.var = 'support', fill = 0)
sum(rowSums(df_all_donors.summary.dcast[, 7:10]))



df_all_donors.summary.dcast.bak <- df_all_donors.summary.dcast
write.table(df_all_donors.summary.dcast.bak, 'df_all_donors.tsv', quote = F, sep = '\t', row.names = F)



### start from 
df_all_donors.summary.dcast.bak <- read.delim('df_all_donors.tsv', check.names = F)

df_all_donors.summary.dcast <- subset(df_all_donors.summary.dcast.bak, v1_editSite == F | v2_editSite == F)
df_all_donors.summary.dcast.sum <- df_all_donors.summary.dcast %>% group_by(var, v1_editSite, v2_editSite, v1_SNP, v2_SNP) %>% summarise(V2 = 'Sum', `00` = sum(`00`), `01` = sum(`01`), `10` = sum(`10`), `11` = sum(`11`))
colnames(df_all_donors.summary.dcast)
colnames(df_all_donors.summary.dcast.sum)
dim(df_all_donors.summary.dcast)
dim(df_all_donors.summary.dcast.sum)
df_all_donors.summary.dcast <- rbind(df_all_donors.summary.dcast, df_all_donors.summary.dcast.sum)



###############



df_all_donors.summary.dcast$variant_ref_edit_freq <- df_all_donors.summary.dcast$`10`/(df_all_donors.summary.dcast$`10`+df_all_donors.summary.dcast$`00`)
df_all_donors.summary.dcast$variant_alt_edit_freq <- df_all_donors.summary.dcast$`11`/(df_all_donors.summary.dcast$`11`+df_all_donors.summary.dcast$`01`)

df_all_donors.summary.dcast$variant_ref_edit_freq2 <- df_all_donors.summary.dcast$`01`/(df_all_donors.summary.dcast$`01`+df_all_donors.summary.dcast$`00`)
df_all_donors.summary.dcast$variant_alt_edit_freq2 <- df_all_donors.summary.dcast$`11`/(df_all_donors.summary.dcast$`11`+df_all_donors.summary.dcast$`10`)



df_all_donors.summary.dcast$sum <- rowSums(df_all_donors.summary.dcast[, 7:10])
df_all_donors.summary.dcast.varSum <- df_all_donors.summary.dcast %>% group_by(var) %>% summarise(varSum = sum(sum))
df_all_donors.summary.dcast.V2Sum <- df_all_donors.summary.dcast %>% group_by(V2) %>% summarise(varSum = sum(sum))
sum(df_all_donors.summary.dcast$sum)

df_all_donors.summary.dcast.melt <- melt(df_all_donors.summary.dcast, id.vars = 1:10, measure.vars = 11:12)
df_all_donors.summary.dcast.melt2 <- melt(df_all_donors.summary.dcast, id.vars = 1:10, measure.vars = 13:14)
df_all_donors.summary.dcast.melt$variant_sum <- df_all_donors.summary.dcast.melt$`10`+df_all_donors.summary.dcast.melt$`00`
df_all_donors.summary.dcast.melt$variant_sum[which(df_all_donors.summary.dcast.melt$variable == 'variant_alt_edit_freq')] <- c(df_all_donors.summary.dcast.melt$`11`+df_all_donors.summary.dcast.melt$`01`)[which(df_all_donors.summary.dcast.melt$variable == 'variant_alt_edit_freq')]
df_all_donors.summary.dcast.melt2$variant_sum <- df_all_donors.summary.dcast.melt2$`01`+df_all_donors.summary.dcast.melt2$`00`
df_all_donors.summary.dcast.melt2$variant_sum[which(df_all_donors.summary.dcast.melt2$variable == 'variant_alt_edit_freq2')] <- c(df_all_donors.summary.dcast.melt2$`11`+df_all_donors.summary.dcast.melt2$`10`)[which(df_all_donors.summary.dcast.melt2$variable == 'variant_alt_edit_freq2')]

df_all_donors.summary.dcast.melt$fisher_p <- apply(df_all_donors.summary.dcast.melt, 1, function(x){
  #typeof(as.numeric(x[7:10]))
  fisher.test(matrix(as.numeric(x[7:10]), nrow = 2))$p
})
df_all_donors.summary.dcast.melt2$fisher_p <- df_all_donors.summary.dcast.melt$fisher_p

df_all_donors.summary.dcast.melt.noB <- subset(df_all_donors.summary.dcast.melt, ! V2 %in% c('B', 'atypical_B', 'Plasma_B', 'IGHMlo_memory_B', 'IGHMhi_memory_B', 'flagged_platelet_sum'))
df_all_donors.summary.dcast.melt.noB.noIG <-df_all_donors.summary.dcast.melt.noB[grep('chr2|chr14|chr22', df_all_donors.summary.dcast.melt.noB$var, invert = T), ]
df_all_donors.summary.dcast.melt.noB.noIG.noAGTC <- df_all_donors.summary.dcast.melt.noB.noIG[grep('.*A_G.*A_G|.*T_C.*T_C', df_all_donors.summary.dcast.melt.noB.noIG$var, invert = T), ]

# at least 100 reads for a variant pair & at least 1000 reads for a cell type
# var_filter <- subset(df_all_donors.summary.dcast.varSum, varSum > 100)$var
# V2_filter <- subset(df_all_donors.summary.dcast.V2Sum, varSum > 1000)$V2
# df_all_donors.summary.dcast.melt.eg <- subset(df_all_donors.summary.dcast.melt, var %in% var_filter & V2 %in% V2_filter)
# ggplot(df_all_donors.summary.dcast.melt.eg, aes(V2, value, fill = variable)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = .1) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# ggsave('plot_edit_level.jitter.pdf', width = 6, height = 5)
# 
# ggplot(df_all_donors.summary.dcast.melt.eg, aes(V2, value, fill = variable)) + geom_boxplot(outlier.size = .1) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# ggsave('plot_edit_level.pdf', width = 6, height = 5)





# 
# # example (Asian specific minor allele)
# df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt, var == 'chr8_11843167_._T_C,chr8_11843196_._C_T')
# #df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt, var == 'chr8_11843167_._T_C,chr8_11843218_._A_C')
# df_all_donors.summary.dcast.melt.oneExample.filter <- df_all_donors.summary.dcast.melt.oneExample[rowSums(df_all_donors.summary.dcast.melt.oneExample[, c('00', '01', '10', '11')]) > 50, ]
# 
# ggplot(df_all_donors.summary.dcast.melt.oneExample.filter, aes(V2, value, fill = variable)) + geom_histogram(stat = 'identity', position = 'dodge') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# ggsave('plot_one_example.both_edit_var.chr8_11843167_._T_C,chr8_11843196_._C_T.pdf', width = 6, height = 5)
# 
# 
# # example (Asian specific edQTL)
# df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt, var == 'chr1_204557582_._A_G,chr1_204557593_._G_A')
# df_all_donors.summary.dcast.melt.oneExample.filter <- df_all_donors.summary.dcast.melt.oneExample[rowSums(df_all_donors.summary.dcast.melt.oneExample[, c('00', '01', '10', '11')]) > 50, ]
# 
# ggplot(df_all_donors.summary.dcast.melt.oneExample.filter, aes(V2, value, fill = variable)) + geom_histogram(stat = 'identity', position = 'dodge') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# ggsave('plot_one_example.edit.chr1_204557582_._A_G,chr1_204557593_._G_A.pdf', width = 6, height = 5)
# 
# 
# # example (exactly what we need, but lack of samples in GTEx)
# df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt, var == 'chr19_17415042_._A_G,chr19_17415200_._C_G')
# df_all_donors.summary.dcast.melt.oneExample.filter <- df_all_donors.summary.dcast.melt.oneExample[rowSums(df_all_donors.summary.dcast.melt.oneExample[, c('00', '01', '10', '11')]) > 50, ]
# 
# ggplot(df_all_donors.summary.dcast.melt.oneExample.filter, aes(V2, value, fill = variable)) + geom_histogram(stat = 'identity', position = 'dodge') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# ggsave('plot_one_example.edit.chr19_17415042_._A_G,chr19_17415200_._C_G.pdf', width = 6, height = 5)
# 
# 
# dat <- data.frame(
#   "smoke_no" = c(940, 980),
#   "smoke_yes" = c(384, 522),
#   row.names = c("Athlete", "Non-athlete"),
#   stringsAsFactors = FALSE
# )
# colnames(dat) <- c("Non-smoker", "Smoker")
# chisq.test(dat)
# fisher.test(dat)
# 
# 
# write.table(df_all_donors.summary.dcast.melt, 'df_all_donors.summary.dcast.melt.tsv', quote = F, sep = '\t', row.names = F)
# 
# df_all_donors.summary.dcast.varSum.noIG <- df_all_donors.summary.dcast.varSum[grep('chr2|chr14|chr22', df_all_donors.summary.dcast.varSum$var, invert = T), ]
# df_all_donors.summary.dcast.varSum.noIG.noAGTC <- df_all_donors.summary.dcast.varSum.noIG[grep('.*A_G.*A_G|.*T_C.*T_C', df_all_donors.summary.dcast.varSum.noIG$var, invert = T), ]
# 
# 
# df_all_donors.summary.dcast.varSum <- df_all_donors.summary.dcast.varSum[order(df_all_donors.summary.dcast.varSum$varSum, decreasing = T), ]
# for (n in 1:100) {
#   now.tmp <- df_all_donors.summary.dcast.varSum$var[n]
#   df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt, var == now.tmp)
#   df_all_donors.summary.dcast.melt.oneExample.filter <- df_all_donors.summary.dcast.melt.oneExample[rowSums(df_all_donors.summary.dcast.melt.oneExample[, c('00', '01', '10', '11')]) > 50, ]
#   ggplot(df_all_donors.summary.dcast.melt.oneExample.filter, aes(V2, value, fill = variable)) + geom_histogram(stat = 'identity', position = 'dodge') +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
#   ggsave(paste0('plot_examples/plot_example.edit.', df_all_donors.summary.dcast.varSum.noIG.noAGTC$varSum[n], '.', now.tmp, '.pdf'), width = 6, height = 5)
# }
# 
# df_all_donors.summary.dcast.varSum.noIG <- df_all_donors.summary.dcast.varSum.noIG[order(df_all_donors.summary.dcast.varSum.noIG$varSum, decreasing = T), ]
# for (n in 1:50) {
#   now.tmp <- df_all_donors.summary.dcast.varSum.noIG$var[n]
#   df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt, var == now.tmp)
#   df_all_donors.summary.dcast.melt.oneExample.filter <- df_all_donors.summary.dcast.melt.oneExample[rowSums(df_all_donors.summary.dcast.melt.oneExample[, c('00', '01', '10', '11')]) > 50, ]
#   ggplot(df_all_donors.summary.dcast.melt.oneExample.filter, aes(V2, value, fill = variable)) + geom_histogram(stat = 'identity', position = 'dodge') +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
#   ggsave(paste0('plot_examples_noIG/plot_example.edit.', df_all_donors.summary.dcast.varSum.noIG.noAGTC$varSum[n], '.', now.tmp, '.pdf'), width = 6, height = 5)
# }
# 
# df_all_donors.summary.dcast.varSum.noIG.noAGTC <- df_all_donors.summary.dcast.varSum.noIG.noAGTC[order(df_all_donors.summary.dcast.varSum.noIG.noAGTC$varSum, decreasing = T), ]
# for (n in 1:100) {
#   now.tmp <- df_all_donors.summary.dcast.varSum.noIG.noAGTC$var[n]
#   df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt, var == now.tmp)
#   df_all_donors.summary.dcast.melt.oneExample.filter <- df_all_donors.summary.dcast.melt.oneExample[rowSums(df_all_donors.summary.dcast.melt.oneExample[, c('00', '01', '10', '11')]) > 50, ]
#   ggplot(df_all_donors.summary.dcast.melt.oneExample.filter, aes(V2, value, fill = variable)) + geom_histogram(stat = 'identity', position = 'dodge') +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
#   ggsave(paste0('plot_examples_noIG_noAGTC/plot_example.edit.', df_all_donors.summary.dcast.varSum.noIG.noAGTC$varSum[n], '.', now.tmp, '.pdf'), width = 6, height = 5)
# }
# 
# 
# 
# # (editing site not found)
# now <- 'chr19_18070688_._T_C,chr19_18070759_._A_G'
# 
# df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt, var == now)
# df_all_donors.summary.dcast.melt.oneExample.filter <- df_all_donors.summary.dcast.melt.oneExample[rowSums(df_all_donors.summary.dcast.melt.oneExample[, c('00', '01', '10', '11')]) > 50, ]
# ggplot(df_all_donors.summary.dcast.melt.oneExample.filter, aes(V2, value, fill = variable)) + geom_histogram(stat = 'identity', position = 'dodge') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# 
# # edQTL slope < 0. CD14+ opposite effect
# now2 <- 'chr19_18560245_._C_T,chr19_18560276_._A_G'
# 
# df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt2, var == now2)
# df_all_donors.summary.dcast.melt.oneExample.filter <- df_all_donors.summary.dcast.melt.oneExample[rowSums(df_all_donors.summary.dcast.melt.oneExample[, c('00', '01', '10', '11')]) > 50, ]
# ggplot(df_all_donors.summary.dcast.melt.oneExample.filter, aes(V2, value, fill = variable)) + geom_histogram(stat = 'identity', position = 'dodge') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# 
# # edQTL slope < 0. 
# now3 <- 'chr16_29670368_._A_G,chr16_29670386_._G_A'
# 
# df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt2, var == now3)
# df_all_donors.summary.dcast.melt.oneExample.filter <- df_all_donors.summary.dcast.melt.oneExample[rowSums(df_all_donors.summary.dcast.melt.oneExample[, c('00', '01', '10', '11')]) > 50, ]
# ggplot(df_all_donors.summary.dcast.melt.oneExample.filter, aes(V2, value, fill = variable)) + geom_histogram(stat = 'identity', position = 'dodge') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# 
# #
# now4 <- 'chr12_7923838_._A_G,chr12_7923842_._T_C'
# 
# df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt2, var == now4)
# df_all_donors.summary.dcast.melt.oneExample.filter <- df_all_donors.summary.dcast.melt.oneExample[rowSums(df_all_donors.summary.dcast.melt.oneExample[, c('00', '01', '10', '11')]) > 50, ]
# ggplot(df_all_donors.summary.dcast.melt.oneExample.filter, aes(V2, value, fill = variable)) + geom_histogram(stat = 'identity', position = 'dodge') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# 
# 
# # final examples 
# source("~/Library/CloudStorage/OneDrive-Personal/workspace/LiuLab/proj_aida_sqtl/reference/cell_type_pretty_names.R")
# 
# #final <- c('chr19_18560245_._C_T,chr19_18560276_._A_G', 'chr12_7923838_._A_G,chr12_7923842_._T_C')
# final <- c('chr22_39018761_._C_T,chr22_39018767_._A_G', 
#            'chr16_29670368_._A_G,chr16_29670413_._C_T',
#            'chr19_17415042_._A_G,chr19_17415200_._C_G',
#            'chr13_19672883_._A_G,chr13_19672915_._A_G',
#            'chr11_108367112_._C_T,chr11_108367230_._A_G')
# 
# for (f in final) {
#   now.tmp <- f
#   if (unique(subset(df_all_donors.summary.dcast.melt, var == now.tmp)$v1_editSite) == T) {
#     df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt, var == now.tmp)
#   }else{
#     df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt2, var == now.tmp)
#   }
#   df_all_donors.summary.dcast.melt.oneExample.filter <- df_all_donors.summary.dcast.melt.oneExample[rowSums(df_all_donors.summary.dcast.melt.oneExample[, c('00', '01', '10', '11')]) > 20, ]
#   df_all_donors.summary.dcast.melt.oneExample.filter <- subset(df_all_donors.summary.dcast.melt.oneExample.filter, V2 %in% gsub('\\\\', '', names(label_dict)))
#   for (i in seq_along(label_dict)) {
#     df_all_donors.summary.dcast.melt.oneExample.filter$V2 <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],df_all_donors.summary.dcast.melt.oneExample.filter$V2)
#   }
#   
#   ggplot(df_all_donors.summary.dcast.melt.oneExample.filter) + geom_histogram(aes(V2, value, fill = variable), stat = 'identity', position = 'dodge') +
#     theme_classic() + scale_y_continuous(expand = c(0, 0)) + geom_label(aes(x = V2, y = 0.1, label = signif(fisher_p, 3))) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('') + ylab('Editing level')
#   ggsave(paste0('plot_examples_final/plot_example.edit.', f, '.pdf'), width = 6, height = 5)
#   
#   # write.table(df_all_donors.summary.dcast.melt.oneExample.filter, paste0('table_examples_final/plot_example.edit.', f, '.tsv'), quote = F, sep = '\t', row.names = F)
# }
# 
# 
# fisher_var_id <- intersect(which(rowSums(df_all_donors.summary.dcast.melt.noB.noIG.noAGTC[, c('00', '01', '10', '11')]) > 50), which(df_all_donors.summary.dcast.melt.noB.noIG.noAGTC$fisher_p < .05))
# final_fisher <- unique(df_all_donors.summary.dcast.melt.noB.noIG.noAGTC$var[fisher_var_id])
# final_baoxingde <- c('chr1_184792052_._G_A,chr1_184792216_._T_C',
#                      'chr1_204557582_._A_G,chr1_204557593_._G_A', 
#                      'chr1_36477832_._T_C,chr1_36478054_._G_A',
#                      'chr11_74455801_._T_C,chr11_74455814_._A_G', 
#                      'chr12_7923089_._A_G,chr12_7923267_._T_C',
#                      'chr13_19672883_._A_G,chr13_19672915_._A_G', 
#                      'chr14_106062251_._T_C,chr14_106062599_._C_T',
#                      'chr14_106062272_._T_A,chr14_106062599_._C_T', 
#                      'chr17_45217304_._T_C,chr17_45217347_._T_C', 
#                      'chr19_17415042_._A_G,chr19_17415200_._C_G',
#                      'chr19_18560245_._C_T,chr19_18560270_._A_G', 
#                      'chr2_201287368_._A_G,chr2_201287439_._T_A', 
#                      'chr22_39018761_._C_T,chr22_39018907_._A_G', 
#                      'chr8_11843077_._T_C,chr8_11843218_._A_C',
#                      'chr2_201287383_._A_G,chr2_201287439_._T_A')
# # for (f in final_fisher) {
# for (f in final_baoxingde) {
#     now.tmp <- f
#   if (unique(subset(df_all_donors.summary.dcast.melt, var == now.tmp)$v1_editSite) == T) {
#     df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt, var == now.tmp)
#     df_all_donors.summary.dcast.melt.oneExample.filter <- df_all_donors.summary.dcast.melt.oneExample[rowSums(df_all_donors.summary.dcast.melt.oneExample[, c('00', '01', '10', '11')]) > 50, ]
#     ref_edit_00 <- sum(df_all_donors.summary.dcast.melt.oneExample.filter[, 7])/2
#     ref_edit_01 <- sum(df_all_donors.summary.dcast.melt.oneExample.filter[, 8])/2
#     ref_edit_10 <- sum(df_all_donors.summary.dcast.melt.oneExample.filter[, 9])/2
#     ref_edit_11 <- sum(df_all_donors.summary.dcast.melt.oneExample.filter[, 10])/2
#     ref_edit_level <- ref_edit_10/(ref_edit_10+ref_edit_00)
#     alt_edit_level <- ref_edit_11/(ref_edit_11+ref_edit_01)
#   }else{
#     df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt2, var == now.tmp)
#     df_all_donors.summary.dcast.melt.oneExample.filter <- df_all_donors.summary.dcast.melt.oneExample[rowSums(df_all_donors.summary.dcast.melt.oneExample[, c('00', '01', '10', '11')]) > 50, ]
#     ref_edit_00 <- sum(df_all_donors.summary.dcast.melt.oneExample.filter[, 7])/2
#     ref_edit_01 <- sum(df_all_donors.summary.dcast.melt.oneExample.filter[, 8])/2
#     ref_edit_10 <- sum(df_all_donors.summary.dcast.melt.oneExample.filter[, 9])/2
#     ref_edit_11 <- sum(df_all_donors.summary.dcast.melt.oneExample.filter[, 10])/2
#     ref_edit_level <- ref_edit_01/(ref_edit_01+ref_edit_00)
#     alt_edit_level <- ref_edit_11/(ref_edit_11+ref_edit_10)
#   }
#   df_all_donors.summary.dcast.melt.oneExample.filter <- subset(df_all_donors.summary.dcast.melt.oneExample.filter, V2 %in% gsub('\\\\', '', names(label_dict)))
#   for (i in seq_along(label_dict)) {
#     df_all_donors.summary.dcast.melt.oneExample.filter$V2 <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],df_all_donors.summary.dcast.melt.oneExample.filter$V2)
#   }
#   
#   df_plot <- df_all_donors.summary.dcast.melt.oneExample.filter[, c('V2', 'variable', 'value', 'fisher_p')]
#   df_plot <- rbind(df_plot, data.frame(V2 = c('Sum', 'Sum'), 
#                                        variable = c(unique(df_plot$variable)),
#                                        value = c(ref_edit_level, alt_edit_level),
#                                        fisher_p = rep(fisher.test(matrix(c(ref_edit_00, ref_edit_01, ref_edit_10, ref_edit_11), nrow = 2))$p, 2)))
#   if (length(unique(df_all_donors.summary.dcast.melt.oneExample.filter$V2)) > 2) {
#   #   ggplot(df_all_donors.summary.dcast.melt.oneExample.filter) + geom_histogram(aes(V2, value, fill = variable), stat = 'identity', position = 'dodge') +
#   #     theme_classic() + scale_y_continuous(expand = c(0, 0)) + geom_label(aes(x = V2, y = 0.1, label = signif(fisher_p, 3))) +
#   #     theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('') + ylab('Editing level') +
#   #     ggtitle(unique(df_all_donors.summary.dcast.melt.oneExample.filter$var))
#   #   ggsave(paste0('plot_examples_final_fisher/plot_example.edit.fisher.', f, '.pdf'), width = 6, height = 5)
#     
#     ggplot(df_plot) + geom_histogram(aes(V2, value, fill = variable, color = ifelse(fisher_p < 0.05, '<005', '>=005')), stat = 'identity', position = 'dodge') +
#       theme_classic() + scale_y_continuous(expand = c(0, 0)) + geom_label(aes(x = V2, y = 0.1, label = signif(fisher_p, 3))) +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('') + ylab('Editing level') + 
#       ggtitle(unique(df_all_donors.summary.dcast.melt.oneExample.filter$var)) + scale_color_manual(values = c('black', 'white'))
#     # ggsave(paste0('plot_examples_final_fisher/plot_example.edit.fisher.colored.', f, '.pdf'), width = 6, height = 5)
#     ggsave(paste0('plot_examples_final_baoxingde/plot_example.edit.fisher.colored.', f, '.pdf'), width = 6, height = 5)
#   }
#   # write.table(df_all_donors.summary.dcast.melt.oneExample.filter, paste0('table_examples_final/plot_example.edit.', f, '.tsv'), quote = F, sep = '\t', row.names = F)
# }





final_baoxingde_v0 <- c('chr1_204557582_._A_G,chr1_204557593_._G_A',
                        'chr2_201287368_._A_G,chr2_201287439_._T_A',
                        'chr2_201287383_._A_G,chr2_201287439_._T_A',
                        'chr8_11843196_._C_T,chr8_11843228_._T_C',
                        'chr1_204557593_._G_A,chr1_204557667_._A_G')

# 
# for (f in final_baoxingde_v0) {
#   now.tmp <- f
#   if (unique(subset(df_all_donors.summary.dcast.melt, var == now.tmp)$v1_editSite) == T) {
#     df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt, var == now.tmp)
#     df_all_donors.summary.dcast.melt.oneExample.filter <- df_all_donors.summary.dcast.melt.oneExample[rowSums(df_all_donors.summary.dcast.melt.oneExample[, c('00', '01', '10', '11')]) > 50, ]
#     ref_edit_00 <- sum(df_all_donors.summary.dcast.melt.oneExample.filter[, 7])/2
#     ref_edit_01 <- sum(df_all_donors.summary.dcast.melt.oneExample.filter[, 8])/2
#     ref_edit_10 <- sum(df_all_donors.summary.dcast.melt.oneExample.filter[, 9])/2
#     ref_edit_11 <- sum(df_all_donors.summary.dcast.melt.oneExample.filter[, 10])/2
#     ref_edit_level <- ref_edit_10/(ref_edit_10+ref_edit_00)
#     alt_edit_level <- ref_edit_11/(ref_edit_11+ref_edit_01)
#   }else{
#     df_all_donors.summary.dcast.melt.oneExample <- subset(df_all_donors.summary.dcast.melt2, var == now.tmp)
#     df_all_donors.summary.dcast.melt.oneExample.filter <- df_all_donors.summary.dcast.melt.oneExample[rowSums(df_all_donors.summary.dcast.melt.oneExample[, c('00', '01', '10', '11')]) > 50, ]
#     ref_edit_00 <- sum(df_all_donors.summary.dcast.melt.oneExample.filter[, 7])/2
#     ref_edit_01 <- sum(df_all_donors.summary.dcast.melt.oneExample.filter[, 8])/2
#     ref_edit_10 <- sum(df_all_donors.summary.dcast.melt.oneExample.filter[, 9])/2
#     ref_edit_11 <- sum(df_all_donors.summary.dcast.melt.oneExample.filter[, 10])/2
#     ref_edit_level <- ref_edit_01/(ref_edit_01+ref_edit_00)
#     alt_edit_level <- ref_edit_11/(ref_edit_11+ref_edit_10)
#   }
#   df_all_donors.summary.dcast.melt.oneExample.filter <- subset(df_all_donors.summary.dcast.melt.oneExample.filter, V2 %in% gsub('\\\\', '', names(label_dict)))
#   for (i in seq_along(label_dict)) {
#     df_all_donors.summary.dcast.melt.oneExample.filter$V2 <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],df_all_donors.summary.dcast.melt.oneExample.filter$V2)
#   }
#   
#   df_plot <- df_all_donors.summary.dcast.melt.oneExample.filter[, c('V2', 'variable', 'value', 'fisher_p', 'variant_sum')]
#   df_plot <- rbind(df_plot, data.frame(V2 = c('Sum', 'Sum'), 
#                                        variable = c(unique(df_plot$variable)),
#                                        value = c(ref_edit_level, alt_edit_level),
#                                        fisher_p = rep(fisher.test(matrix(c(ref_edit_00, ref_edit_01, ref_edit_10, ref_edit_11), nrow = 2))$p, 2),
#                                        variant_sum = c(sum(df_plot[grep('variant_ref_edit', df_plot$variable), 'variant_sum']), sum(df_plot[grep('variant_alt_edit', df_plot$variable), 'variant_sum']))))
#   df_plot$FDR <- signif(p.adjust(df_plot$fisher_p, method = 'fdr'), 3)
#   
#   if (length(unique(df_all_donors.summary.dcast.melt.oneExample.filter$V2)) > 2 & min(df_plot$FDR) < .05) {
#     #   ggplot(df_all_donors.summary.dcast.melt.oneExample.filter) + geom_histogram(aes(V2, value, fill = variable), stat = 'identity', position = 'dodge') +
#     #     theme_classic() + scale_y_continuous(expand = c(0, 0)) + geom_label(aes(x = V2, y = 0.1, label = signif(fisher_p, 3))) +
#     #     theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('') + ylab('Editing level') +
#     #     ggtitle(unique(df_all_donors.summary.dcast.melt.oneExample.filter$var))
#     #   ggsave(paste0('plot_examples_final_fisher/plot_example.edit.fisher.', f, '.pdf'), width = 6, height = 5)
#     df_plot$variable <- gsub('variant_ref_edit_freq|variant_ref_edit_freq2', 'Ref.', df_plot$variable)
#     df_plot$variable <- gsub('variant_alt_edit_freq|variant_alt_edit_freq2', 'Alt.', df_plot$variable)
#     df_plot$variable <- factor(df_plot$variable, levels = c('Ref.', 'Alt.'))
#     df_p_val <- df_plot[, c('V2', 'variable', 'value', 'FDR')] %>% 
#       group_by(V2, FDR) %>%
#       summarise(y.position = max(value) + .15, group1 = 'Ref.', group2 = 'Alt.') %>% arrange(V2, .locale="en") %>%
#       rstatix::add_x_position(x = 'V2', group = 'variable', dodge = 0.5)
#     df_p_val <- subset(df_p_val, FDR < .05)
#     
#     ggplot(df_plot, aes(V2, value)) + 
#       geom_point(aes(color = variable), stat = 'identity', position=position_dodge(width=c(0.5))) +
#       geom_errorbar(aes(color = variable, ymin = value - sqrt(value*(1-value)/variant_sum), ymax = value + sqrt(value*(1-value)/variant_sum)), position=position_dodge(width=c(0.5)), width= .2) + 
#       theme_classic() + scale_y_continuous(expand = c(0, 0)) + 
#       theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('') + ylab('Editing level') + 
#       ggtitle(unique(df_all_donors.summary.dcast.melt.oneExample.filter$var)) + scale_color_manual(values = c('lightblue', 'pink')) +
#       ylim(0,1) + add_pvalue(data = df_p_val, xmin = 'xmin', xmax = 'xmax', label = 'p = {FDR}', tip.length = 0.01, bracket.size = .3)
#     #  geom_label(aes(x = V2, y = 0.1, label = signif(fisher_p, 3))) +
#     # ggsave(paste0('plot_examples_final_fisher/plot_example.edit.fisher.colored.', f, '.pdf'), width = 6, height = 5)
#     ggsave(paste0('plot_examples_final_/plot_example.edit.fisher.colored.width.', f, '.pdf'), width = 2+length(unique(df_all_donors.summary.dcast.melt.oneExample.filter$V2))*0.5, height = 5)
#     write.table(df_plot, paste0('plot_examples_final_/table_example.edit.', gsub(',', '__', f), '.tsv'), quote = F, sep = '\t', row.names = F)
#   }
# }

library("viridis")
df_all_donors.summary.dcast.melt.sum50 <- subset(df_all_donors.summary.dcast.melt, variant_sum > 50)
df_all_donors.summary.dcast.melt.sum50$FDR <- signif(p.adjust(df_all_donors.summary.dcast.melt.sum50$fisher_p, method = 'fdr'), 3)

df_all_donors.summary.dcast.melt2.sum50 <- subset(df_all_donors.summary.dcast.melt2, variant_sum > 50)
df_all_donors.summary.dcast.melt2.sum50$FDR <- signif(p.adjust(df_all_donors.summary.dcast.melt2.sum50$fisher_p, method = 'fdr'), 3)

for (f in final_baoxingde_v0) {
  now.tmp <- f
  if (unique(subset(df_all_donors.summary.dcast.melt.sum50, var == now.tmp)$v1_editSite) == T) {
    df_all_donors.summary.dcast.melt.sum50.oneExample <- subset(df_all_donors.summary.dcast.melt.sum50, var == now.tmp)
    ref_edit_00 <- sum(df_all_donors.summary.dcast.melt.sum50.oneExample[, 7])/2
    ref_edit_01 <- sum(df_all_donors.summary.dcast.melt.sum50.oneExample[, 8])/2
    ref_edit_10 <- sum(df_all_donors.summary.dcast.melt.sum50.oneExample[, 9])/2
    ref_edit_11 <- sum(df_all_donors.summary.dcast.melt.sum50.oneExample[, 10])/2
    ref_edit_level <- ref_edit_10/(ref_edit_10+ref_edit_00)
    alt_edit_level <- ref_edit_11/(ref_edit_11+ref_edit_01)
    edit <- unique(strsplit(df_all_donors.summary.dcast.melt.sum50.oneExample$var, ',')[[1]][1])
    variant <- unique(strsplit(df_all_donors.summary.dcast.melt.sum50.oneExample$var, ',')[[1]][2])
  }else{
    df_all_donors.summary.dcast.melt.sum50.oneExample <- subset(df_all_donors.summary.dcast.melt2.sum50, var == now.tmp)
    ref_edit_00 <- sum(df_all_donors.summary.dcast.melt.sum50.oneExample[, 7])/2
    ref_edit_01 <- sum(df_all_donors.summary.dcast.melt.sum50.oneExample[, 8])/2
    ref_edit_10 <- sum(df_all_donors.summary.dcast.melt.sum50.oneExample[, 9])/2
    ref_edit_11 <- sum(df_all_donors.summary.dcast.melt.sum50.oneExample[, 10])/2
    ref_edit_level <- ref_edit_01/(ref_edit_01+ref_edit_00)
    alt_edit_level <- ref_edit_11/(ref_edit_11+ref_edit_10)
    edit <- unique(strsplit(df_all_donors.summary.dcast.melt.sum50.oneExample$var, ',')[[1]][2])
    variant <- unique(strsplit(df_all_donors.summary.dcast.melt.sum50.oneExample$var, ',')[[1]][1])
  }
  df_all_donors.summary.dcast.melt.sum50.oneExample <- subset(df_all_donors.summary.dcast.melt.sum50.oneExample, V2 %in% c(gsub('\\\\', '', names(label_dict)), 'Sum'))
  for (i in seq_along(label_dict)) {
    df_all_donors.summary.dcast.melt.sum50.oneExample$V2 <- gsub(paste0("^",names(label_dict)[i],"$"),label_dict[[i]],df_all_donors.summary.dcast.melt.sum50.oneExample$V2)
  }
  
  df_plot <- df_all_donors.summary.dcast.melt.sum50.oneExample[, c('V2', 'variable', 'value', 'FDR', 'variant_sum')]
  df_plot <- subset(df_plot, V2 %in% names(table(df_plot$V2)[table(df_plot$V2) == 2]))

  if (length(unique(df_all_donors.summary.dcast.melt.sum50.oneExample$V2)) > 2 & min(df_plot$FDR) < .1) {
    df_plot$variable <- gsub('variant_ref_edit_freq|variant_ref_edit_freq2', 'Ref.', df_plot$variable)
    df_plot$variable <- gsub('variant_alt_edit_freq|variant_alt_edit_freq2', 'Alt.', df_plot$variable)
    df_plot$variable <- factor(df_plot$variable, levels = c('Ref.', 'Alt.'))
    df_p_val <- df_plot[, c('V2', 'variable', 'value', 'FDR')] %>% 
      group_by(V2, FDR) %>%
      summarise(y.position = max(value) + .15, group1 = 'Ref.', group2 = 'Alt.') %>% arrange(V2, .locale="en") %>%
      rstatix::add_x_position(x = 'V2', group = 'variable', dodge = 0.5)
    df_p_val <- subset(df_p_val, FDR < .1)
    
    df_plot$std <- sqrt(df_plot$value*(1-df_plot$value)/df_plot$variant_sum)
    ggplot(df_plot, aes(V2, value)) + 
      geom_point(aes(color = variable), stat = 'identity', position=position_dodge(width=c(0.5))) +
      geom_errorbar(aes(color = variable, ymin = value - std, ymax = value + std), position=position_dodge(width=c(0.5)), width= .2) + 
      theme_classic() + scale_y_continuous(expand = c(0, 0)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab('') + ylab('Editing level') + 
      ggtitle(paste0('Edit site: ', edit, '\n', 'Variant: ', variant)) + 
      scale_color_manual(values = c('blue', 'red')) +
      ylim(0,1) + add_pvalue(data = df_p_val, xmin = 'xmin', xmax = 'xmax', label = '{FDR}', tip.length = 0.01, bracket.size = .3)
    #  geom_label(aes(x = V2, y = 0.1, label = signif(fisher_p, 3))) +
    # ggsave(paste0('plot_examples_final_fisher/plot_example.edit.fisher.colored.', f, '.pdf'), width = 6, height = 5)
    ggsave(paste0('plot_examples_final_FDR/plot_example.edit.fisher.colored.width.', f, '.pdf'), width = 2+length(unique(df_all_donors.summary.dcast.melt.oneExample.filter$V2))*0.5, height = 5)
    write.table(df_plot, paste0('plot_examples_final_FDR/table_example.edit.', gsub(',', '__', f), '.tsv'), quote = F, sep = '\t', row.names = F)
  }
}
# library(bootstrap)
# x <- rnorm(20)
# theta <- function(x){mean(x)}
# results <- jackknife(x,theta)
# 
# mean(x[1:19]) %in% results$jack.values

dim(df_all_donors.summary.dcast.melt)
df_all_donors.summary.dcast.melt.final <- subset(df_all_donors.summary.dcast.melt, v1_editSite == T)
df_all_donors.summary.dcast.melt2.final <- subset(df_all_donors.summary.dcast.melt2, v2_editSite == T)
dim(df_all_donors.summary.dcast.melt.final)
dim(df_all_donors.summary.dcast.melt2.final)

df_all_donors.summary.dcast.melt1and2.final <- rbind(df_all_donors.summary.dcast.melt.final, df_all_donors.summary.dcast.melt2.final)
write.table(df_all_donors.summary.dcast.melt1and2.final[, grep('_SNP', colnames(df_all_donors.summary.dcast.melt1and2.final), invert = T, value = T)], 
            'plot_examples_final_FDR/table_df_all_donors.summary.dcast.melt1and2.final.tsv', sep = '\t', quote = F, row.names = F)


dim(df_all_donors.summary.dcast.melt.sum50)
df_all_donors.summary.dcast.melt.sum50.final <- subset(df_all_donors.summary.dcast.melt.sum50, v1_editSite == T)
df_all_donors.summary.dcast.melt2.sum50.final <- subset(df_all_donors.summary.dcast.melt2.sum50, v2_editSite == T)
dim(df_all_donors.summary.dcast.melt.sum50.final)
dim(df_all_donors.summary.dcast.melt2.sum50.final)

df_all_donors.summary.dcast.melt1and2.sum50.final <- rbind(df_all_donors.summary.dcast.melt.sum50.final, df_all_donors.summary.dcast.melt2.sum50.final)
write.table(df_all_donors.summary.dcast.melt1and2.sum50.final[, grep('_SNP', colnames(df_all_donors.summary.dcast.melt1and2.sum50.final), invert = T, value = T)], 
            'plot_examples_final_FDR/table_df_all_donors.summary.dcast.melt1and2.sum50.final.tsv', sep = '\t', quote = F, row.names = F)

table(df_all_donors.summary.dcast.melt1and2.sum50.final$V2)
df_all_donors.summary.dcast.melt1and2.sum50.final.noB <- subset(df_all_donors.summary.dcast.melt1and2.sum50.final, ! V2 %in% c('B', 'atypical_B', 'naive_B', 'Plasma_B', 'IGHMlo_memory_B', 'IGHMhi_memory_B', 'flagged_platelet_sum'))
df_all_donors.summary.dcast.melt1and2.sum50.final.noB$FDR2 <- p.adjust(df_all_donors.summary.dcast.melt1and2.sum50.final.noB$fisher_p, method = 'fdr')






# cell type editing level
df_all_donors.summary.dcast.melt1and2.sum50.final.summary <- df_all_donors.summary.dcast.melt1and2.sum50.final %>% group_by(V2) %>% summarise(mean = mean(value), sd = sd(value), nrow = n()) %>% arrange(mean)
write.table(df_all_donors.summary.dcast.melt1and2.sum50.final.summary, 'plot_examples_final_FDR/table_df_all_donors.summary.dcast.melt1and2.sum50.final.summary.tsv', quote = F, sep = '\t', row.names = F)













# %
df_all_donors.summary.dcast.melt1and2.final <- read.delim('plot_examples_final_FDR/table_df_all_donors.summary.dcast.melt1and2.final.tsv', check.names = F)
gtex.edQTL <- read.delim(gzfile('~/Downloads/dataset/edQTL/Whole-Blood.signif_variant_site_pairs.txt.gz'))
gtex.edQTL$variant <- gsub('_[A-Z]*_[A-Z]*_b38$', '', gtex.edQTL$variant_id)

length(unique(paste0(gtex.edQTL$variant_id, '__', gtex.edQTL$gene_id)))
length(unique(gtex.edQTL$variant_id))
length(unique(gtex.edQTL$gene_id))

gtex_pairs <- unique(paste0(gtex.edQTL$variant, '__', gtex.edQTL$gene_id))
phased_pairs <- gsub('_\\._.*', '', gsub('_\\._.*,', '__', df_all_donors.summary.dcast.melt1and2.final$var))
phased_pairs_reverse <- paste0(str_split_fixed(phased_pairs, '__', 2)[, 2], '__', str_split_fixed(phased_pairs, '__', 2)[, 1])
phased_pairs_all <- c(phased_pairs, phased_pairs_reverse)
table(gtex_pairs %in% phased_pairs_all)


gtex_pair_distance <- data.frame(distance = as.numeric(gsub('.*_([0-9]*)__.*', '\\1', gtex_pairs)) - as.numeric(gsub('.*_([0-9]*)$*', '\\1', gtex_pairs)))
ggplot(gtex_pair_distance, aes(abs(distance))) + geom_histogram() + scale_x_log10()

gtex_pairs.filt <- gtex_pairs[which(abs(as.numeric(gsub('.*_([0-9]*)__.*', '\\1', gtex_pairs)) - as.numeric(gsub('.*_([0-9]*)$*', '\\1', gtex_pairs))) <= 1000)]
table(gtex_pairs.filt %in% phased_pairs_all)
table(gtex_pairs.filt %in% phased_pairs_all)[2]/table(gtex_pairs.filt %in% phased_pairs_all)[1]
