#### Correlation between transcriptional modules & lung immune responses ####
library(pheatmap)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(knitr)
library(pROC)
library(ggROC)
library(Hmisc)

bcg1pheno <- fread('/labs/khatrilab/yiranliu/BCG/github/data/bcg_routecohort_metadata.csv')
bcg2pheno <- fread('/labs/khatrilab/yiranliu/BCG/github/data/ivbcg_dosecohort_metadata.csv', header = T)
setnames(bcg1pheno, c('description','animal_id','dose_group','protect_outcome'), c('SampleID','id','dose','protect_status'))
setnames(bcg2pheno, c('description','animal_id','dose_group','protect_outcome'), c('SampleID','id','dose','protect_status'))

bcg1pheno$TimeafterBCG <- factor(bcg1pheno$TimeafterBCG, levels = c('pre','d2','wk2','wk12'))
bcg2pheno$TimeafterBCG <- factor(bcg2pheno$TimeafterBCG, levels = c('pre','d2','wk2','wk4','wk12'))
bcg1pheno$dose <- factor(bcg1pheno$dose, levels = c('low','high'))
bcg2pheno$dose <- factor(bcg2pheno$dose, levels = c('low','high'))

vsdmat_bcg1 <- read.csv('/labs/khatrilab/yiranliu/BCG/github/data/bcg_routecohort_processed.csv', row.names = 1, check.names = F)
vsdmat_bcg2 <- read.csv('/labs/khatrilab/yiranliu/BCG/github/data/ivbcg_dosecohort_processed.csv', row.names = 1, check.names = F)

vsdmat_bcg1_norm <- norm_to_baseline(vsdmat_bcg1, pheno_dt = bcg1pheno)
vsdmat_bcg2_norm <- norm_to_baseline(vsdmat_bcg2, pheno_dt = bcg2pheno)

# add imputed baseline for animal with missing baseline sample
nobaseline <- unique(bcg2pheno[TimeafterBCG %in%
                                 c('pre','d2','wk2','wk4','wk12'),id])[which(!unique(bcg2pheno[TimeafterBCG %in% 
                                                                                                 c('pre','d2','wk2','wk4','wk12'),id]) %in% bcg2pheno[TimeafterBCG == 'pre', id])]
vsdmat_bcg2_withimputed <- vsdmat_bcg2
vsdmat_bcg2_withimputed[,paste0(nobaseline,'_pre_imput')] <- rowMedians(as.matrix(vsdmat_bcg2[,bcg2pheno[TimeafterBCG == 'pre', SampleID]]))

fakesamp <- bcg2pheno[id == nobaseline & TimeafterBCG == 'd2']
fakesamp[,`:=`(SampleID=paste0(nobaseline,'_pre_imput'), TimeafterBCG='pre')]
bcg2pheno_withimputed <- rbind(bcg2pheno, fakesamp)

modscores <- fread('/labs/khatrilab/yiranliu/BCG/github//wgcna/modScores.csv')
modscores_bcg1 <- fread('/labs/khatrilab/yiranliu/BCG/github/wgcna/modScores_routecohort.csv')

modscores_norm_cast <- dcast.data.table(modscores[TimeafterBCG %in% c('d2','wk2','wk4')], id + protect_status ~ moduleID + TimeafterBCG, value.var = 'norm_score')

##### Correlation with lung immune responses #####
immune_dt <- fread('/labs/khatrilab/yiranliu/BCG/github/data/immune_data.csv')

### % of CD4 cells in BAL secreting cytokines upon restimulation with WCL
CD4_resp_col <- setdiff(intersect(grep('FoP', colnames(immune_dt), value = T), grep('CD4', colnames(immune_dt), value = T)), 
                        c(grep('Bool|nAUC', colnames(immune_dt), value = T),
                          grep('/Any', colnames(immune_dt), fixed=T, value = T)))
CD4_resp_dt <- immune_dt[,c('id', 'log_Mtb_CFU', 'dose', 'BCG_dose_log10','protect_status', CD4_resp_col), with=F]
CD4_resp_dt <- melt.data.table(CD4_resp_dt, id.vars = c('id', 'log_Mtb_CFU', 'dose', 'BCG_dose_log10','protect_status'))

CD4_resp_dt[,c('cytokine','time') := tstrsplit(variable, '_', fixed=T)]
CD4_resp_dt$time <- factor(CD4_resp_dt$time, levels = sort(unique(as.integer(CD4_resp_dt$time))))

CD4_resp_dt$cytokine <- gsub('FoP|CD4|CD8|Vg9|Marginal','', CD4_resp_dt$cytokine)
CD4_resp_dt$cytokine <- gsub('/','', fixed=T, CD4_resp_dt$cytokine)
CD4_resp_dt[,norm_value := value - value[time == 0], by=.(id,cytokine)]

CD4_cast <- dcast.data.table(CD4_resp_dt[time != 0], id~cytokine+time, value.var = 'norm_value')

mod_CD4_cor <- Hmisc::rcorr(as.matrix(modscores_norm_cast[,grepl('M1|M2|M3|M4|M5|M6|M7|M8',colnames(modscores_norm_cast)),with=F]), 
                            as.matrix(CD4_cast[,!c('id'),with=F]))
mod_CD4_cor <- lapply(mod_CD4_cor, function(mat){
  mod_idx <- which(grepl('M1|M2|M3|M4|M5|M6|M7|M8', colnames(mat)))
  mat <- data.table(mat[mod_idx,-mod_idx], keep.rownames = 'x')
  return(mat)
})

mod_CD4_cor$r_melt <- melt.data.table(mod_CD4_cor$r, variable.name = 'y', value.name = 'r')
mod_CD4_cor$p_melt <- melt.data.table(mod_CD4_cor$P, variable.name = 'y', value.name = 'p')

mod_CD4_cor_dt <- merge.data.table(mod_CD4_cor$r_melt, mod_CD4_cor$p_melt, by=c('x','y'))
mod_CD4_cor_dt[,c('moduleID','mod_time') := tstrsplit(x, '_')]
mod_CD4_cor_dt[,c('cytokine','cyto_time') := tstrsplit(y, '_')]
mod_CD4_cor_dt$cyto_time <- factor(mod_CD4_cor_dt$cyto_time, levels = c(2,4,8,12),
                                   labels = c('wk2','wk4','wk8','wk12'))


### % of CD8 cells in BAL secreting cytokines upon restimulation with WCL
CD8_resp_col <- setdiff(intersect(grep('FoP', colnames(immune_dt), value = T), grep('CD8', colnames(immune_dt), value = T)), 
                        c(grep('Bool|nAUC', colnames(immune_dt), value = T),
                          grep('/Any/', colnames(immune_dt), value = T)))
CD8_resp_dt <- immune_dt[,c('id', 'log_Mtb_CFU', 'dose', 'protect_status', CD8_resp_col), with=F]
CD8_resp_dt <- melt.data.table(CD8_resp_dt, id.vars = c('id', 'log_Mtb_CFU', 'dose', 'protect_status'))

CD8_resp_dt[,c('cytokine','time') := tstrsplit(variable, '_', fixed=T)]
CD8_resp_dt$time <- factor(CD8_resp_dt$time, levels = sort(unique(as.integer(CD8_resp_dt$time))))

CD8_resp_dt$cytokine <- gsub('FoP|CD8|CD8|Vg9|Marginal','', CD8_resp_dt$cytokine)
CD8_resp_dt$cytokine <- gsub('/','', fixed=T, CD8_resp_dt$cytokine)
CD8_resp_dt[,norm_value := value - value[time == 0], by=.(id,cytokine)]

CD8_cast <- dcast.data.table(CD8_resp_dt[time != 0], id~cytokine+time, value.var = 'norm_value')

mod_CD8_cor <- Hmisc::rcorr(as.matrix(modscores_norm_cast[,grepl('M1|M2|M3|M4|M5|M6|M7|M8',colnames(modscores_norm_cast)),with=F]), 
                            as.matrix(CD8_cast[,!c('id'),with=F]))
mod_CD8_cor <- lapply(mod_CD8_cor, function(mat){
  mod_idx <- which(grepl('M1|M2|M3|M4|M5|M6|M7|M8', colnames(mat)))
  mat <- data.table(mat[mod_idx,-mod_idx], keep.rownames = 'x')
  return(mat)
})

mod_CD8_cor$r_melt <- melt.data.table(mod_CD8_cor$r, variable.name = 'y', value.name = 'r')
mod_CD8_cor$p_melt <- melt.data.table(mod_CD8_cor$P, variable.name = 'y', value.name = 'p')

mod_CD8_cor_dt <- merge.data.table(mod_CD8_cor$r_melt, mod_CD8_cor$p_melt, by=c('x','y'))
mod_CD8_cor_dt[,c('moduleID','mod_time') := tstrsplit(x, '_')]
mod_CD8_cor_dt[,c('cytokine','cyto_time') := tstrsplit(y, '_')]
mod_CD8_cor_dt$cyto_time <- factor(mod_CD8_cor_dt$cyto_time, levels = c(2,4,8,12),
                                   labels = c('wk2','wk4','wk8','wk12'))

### Cell counts in BAL 
abs_counts_cols <- setdiff(c(grep('log_B_', colnames(immune_dt), value = T),
                             grep('log_CD4_', colnames(immune_dt), value = T),
                             grep('log_CD8_', colnames(immune_dt), value = T),
                             grep('log_Granu', colnames(immune_dt), value = T),
                             grep('log_iNKT', colnames(immune_dt), value = T),
                             grep('log_Mac', colnames(immune_dt), value = T),
                             grep('log_MAIT', colnames(immune_dt), value = T),
                             grep('log_mDC', colnames(immune_dt), value = T),
                             grep('log_NK', colnames(immune_dt), value = T),
                             grep('log_pDC', colnames(immune_dt), value = T),
                             grep('log_Vg9', colnames(immune_dt), value = T)), grep('nAUC', colnames(immune_dt), value = T))

abs_counts_dt <- immune_dt[,c('id', 'log_Mtb_CFU', 'dose', 'protect_status', abs_counts_cols), with=F]
abs_counts_dt <- melt.data.table(abs_counts_dt, id.vars = c('id', 'log_Mtb_CFU', 'dose','protect_status'))
abs_counts_dt$variable <- gsub('log_', '', abs_counts_dt$variable)
abs_counts_dt[,c('celltype','time') := tstrsplit(variable, '_', fixed=T)]
abs_counts_dt$time <- factor(abs_counts_dt$time, levels = sort(unique(as.integer(abs_counts_dt$time))))
abs_counts_dt[,norm_value := value - value[time == 0], by=c('id','celltype')]

abs_counts_cast <- dcast.data.table(abs_counts_dt[time != 0], id~celltype+time, value.var = 'value')

mod_balcts_cor <- Hmisc::rcorr(as.matrix(modscores_norm_cast[,grepl('M1|M2|M3|M4|M5|M6|M7',colnames(modscores_norm_cast)),with=F]), 
                               as.matrix(abs_counts_cast[,!c('id'),with=F]))
mod_balcts_cor <- lapply(mod_balcts_cor, function(mat){
  mod_idx <- which(grepl('M1|M2|M3|M4|M5|M6|M7', colnames(mat)))
  mat <- data.table(mat[mod_idx,-mod_idx], keep.rownames = 'x')
  return(mat)
})

mod_balcts_cor$r_melt <- melt.data.table(mod_balcts_cor$r, variable.name = 'y', value.name = 'r')
mod_balcts_cor$p_melt <- melt.data.table(mod_balcts_cor$P, variable.name = 'y', value.name = 'p')

mod_balcts_cor_dt <- merge.data.table(mod_balcts_cor$r_melt, mod_balcts_cor$p_melt, by=c('x','y'))
mod_balcts_cor_dt[,c('moduleID','mod_time') := tstrsplit(x, '_')]
mod_balcts_cor_dt[,c('celltype','celltype_time') := tstrsplit(y, '_')]
mod_balcts_cor_dt$celltype_time <- factor(mod_balcts_cor_dt$celltype_time, levels = c(0,2,4,8),
                                          labels = c('pre','wk2','wk4','wk8'))


### Antibodies in BAL and plasma
# IgG or IgA to TB antigens in plasma or BAL
Ab_cols <- setdiff(intersect(grep('IgG|IgA', colnames(immune_dt), value = T),
                             grep('log', colnames(immune_dt), value = T)),
                   grep('nAUC', colnames(immune_dt), value = T))
Ab_dt <- immune_dt[,c('id', 'log_Mtb_CFU', 'dose','protect_status', Ab_cols), with=F]
Ab_dt <- melt.data.table(Ab_dt, id.vars = c('id', 'log_Mtb_CFU', 'dose', 'protect_status'))
Ab_dt$variable <- gsub('log_','',Ab_dt$variable)
Ab_dt[,c('type','time') := tstrsplit(variable, '_', fixed=T)]
Ab_dt$time <- factor(Ab_dt$time, levels = sort(unique(as.integer(Ab_dt$time)))) # 24 wks: time of challenge
Ab_dt_cast <- dcast.data.table(Ab_dt[time %in% c(4)], id~type+time, value.var = 'value')

Ab_dt_cor <- Hmisc::rcorr(as.matrix(modscores_norm_cast[,grepl('M1|M2|M3|M4|M5|M6|M7|M8',colnames(modscores_norm_cast)),with=F]), 
                          as.matrix(Ab_dt_cast[,!c('id'),with=F]))
Ab_dt_cor <- lapply(Ab_dt_cor, function(mat){
  mod_idx <- which(grepl('M1|M2|M3|M4|M5|M6|M7|M8', colnames(mat)))
  mat <- data.table(mat[mod_idx,-mod_idx], keep.rownames = 'x')
  return(mat)
})

Ab_dt_cor$r_melt <- melt.data.table(Ab_dt_cor$r, variable.name = 'y', value.name = 'r')
Ab_dt_cor$p_melt <- melt.data.table(Ab_dt_cor$P, variable.name = 'y', value.name = 'p')

Ab_dt_cor_dt <- merge.data.table(Ab_dt_cor$r_melt, Ab_dt_cor$p_melt, by=c('x','y'))
Ab_dt_cor_dt[,c('moduleID','mod_time') := tstrsplit(x, '_')]
Ab_dt_cor_dt[,c('type','type_time') := tstrsplit(y, '_')]
Ab_dt_cor_dt$type_time <- factor(Ab_dt_cor_dt$type_time, levels = c(4,12),
                                 labels = c('wk4','wk12'))

#### combine to see all correlations
mod_CD4_cor_dt$category <- 'BAL Ag-spec CD4'
mod_CD8_cor_dt$category <- 'BAL Ag-spec CD8'
mod_balcts_cor_dt$category <- 'BAL counts'
Ab_dt_cor_dt[grepl('BAL', type), category := 'BAL titer']
Ab_dt_cor_dt[grepl('Plasma', type), category := 'Plasma titer']

mod_CD4_cor_dt$variable <- mod_CD4_cor_dt$cytokine
mod_CD8_cor_dt$variable <- mod_CD8_cor_dt$cytokine
mod_balcts_cor_dt$variable <- mod_balcts_cor_dt$celltype
Ab_dt_cor_dt$variable <- gsub('BAL.|Plasma.', '', Ab_dt_cor_dt$type)

all_cor <- rbind(mod_CD4_cor_dt[cyto_time == 'wk8'], 
                 mod_CD8_cor_dt[cyto_time == 'wk8'], 
                 mod_balcts_cor_dt[celltype_time == 'wk4' & !celltype %in% c('NK','Granulocytes')], 
                 mod_balcts_cor_dt[celltype_time == 'wk8' & celltype %in% c('NK')],
                 Ab_dt_cor_dt[type_time == 'wk4'], fill=TRUE)

all_cor$fdr <- p.adjust(all_cor$p, method = 'BH')
all_cor$variable <- factor(all_cor$variable, 
                                        levels = c('Anyg2T17','IFNg','IL17','IL2','IL21','TNF',
                                                   'CD4','CD8','MAIT','Vg9','B','NK','iNKT','Mac','mDC','pDC',
                                                   'IgA','IgG'))

ggplot(all_cor, aes(x=mod_time, y=variable, fill=r)) + geom_tile() + 
  scale_x_discrete(expand = c(0,0)) + xlab('') +
  scale_y_discrete(expand = c(0,0), limits=rev) + ylab('') +
  scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0,limits=c(-1,1)) +
  geom_point(data=all_cor[fdr <= 0.01], aes(x=mod_time,y=variable), shape = 8, size = 1) +
  facet_grid(category~moduleID, scales = 'free', space = 'free') + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))

#### scatter plots for notable module score & immune response combinations
CD4_resp_modscores <- merge.data.table(CD4_resp_dt, modscores_norm[,.(id,TimeafterBCG,score,moduleID,norm_score)], by='id', allow.cartesian=T)
CD4_resp_modscores[,mod_time := paste(moduleID, TimeafterBCG, sep='_')]
ggplot(CD4_resp_modscores[time == '8' & mod_time %in% c('M1_d2','M3_wk2','M4_d2')], aes(x=norm_score, y=value)) +
  geom_point(aes(color=mod_time)) + 
  scale_color_manual(values = alpha(get_palette('locuszoom', 8)[c(1,3,4)],0.6)) +
  facet_grid(cytokine~mod_time, scales = 'free') + stat_cor(size=3.5) + theme_classic() + 
  theme(legend.position = 'none', panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.text = element_text(size=12), axis.title = element_text(size=13),
        axis.text.y=element_text(size=10), axis.text.x=element_text(size=10, angle=40, hjust=1)) + 
  xlab('Module summary score (norm. to baseline)') + ylab('% of BAL CD4 responding upon WCL restim.')

CD8_resp_modscores <- merge.data.table(CD8_resp_dt, modscores_norm[,.(id,TimeafterBCG,score,moduleID,norm_score)], by='id', allow.cartesian=T)
CD8_resp_modscores[,mod_time := paste(moduleID, TimeafterBCG, sep='_')]
outliers <- unique(c(CD8_resp_modscores[time == '8' & cytokine %in% c('IL2','IL21') & value > 6, id], # outliers for CD8 responses
                     CD8_resp_modscores[time == '8' & cytokine %in% c('IFNg','TNF') & value > 20, id]))
ggplot(CD8_resp_modscores[time == '8' & mod_time %in% c('M1_d2','M3_wk2','M4_d2') &
                            !id %in% outliers], aes(x=norm_score, y=value)) +
  geom_point(aes(color=mod_time)) + 
  scale_color_manual(values = alpha(get_palette('locuszoom', 8)[c(1,3,4)],0.6)) +
  facet_grid(cytokine~mod_time, scales = 'free') + stat_cor(size=3.5) + theme_classic() + 
  theme(legend.position = 'none', panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.text = element_text(size=12), axis.title = element_text(size=13),
        axis.text.y=element_text(size=10), axis.text.x=element_text(size=10, angle=40, hjust=1)) + 
  xlab('Module summary score (norm. to baseline)') + ylab('% of BAL CD8 responding upon WCL restim.')

abs_counts_modscores <- merge.data.table(abs_counts_dt, modscores_norm[,.(id,TimeafterBCG,score,moduleID,norm_score)], by='id', allow.cartesian=T)
abs_counts_modscores$celltype <- factor(abs_counts_modscores$celltype, 
                                        levels = c('CD4','CD8','MAIT','Vg9','B','NK','iNKT','Granulocytes','Mac','mDC','pDC'))

ggplot(abs_counts_modscores[time == '4' & TimeafterBCG %in% c('d2') & moduleID == 'M1' & celltype %in% c('CD4','CD8','MAIT','Vg9')], aes(x=norm_score, y=value)) + 
  geom_point(aes(color=moduleID)) + scale_color_locuszoom(alpha=0.6) +
  facet_wrap(~celltype, scales = 'free') + stat_cor(size=3.5) + theme_classic() + 
  theme(legend.position = 'none', panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.text = element_text(size=12), axis.title = element_text(size=13),
        axis.text=element_text(size=10)) + 
  xlab('M1 summary score at d2 (norm. to baseline)') + ylab('Number of cells in BAL (log10)')

ggplot(abs_counts_modscores[TimeafterBCG %in% c('d2') & moduleID == 'M1' & 
                              ((time == '4' & celltype %in% c('B','iNKT','Mac','mDC','pDC')) |
                                 (time == '8' & celltype %in% c('NK')))], aes(x=norm_score, y=value)) + 
  geom_point(aes(color=moduleID)) + scale_color_locuszoom(alpha=0.6) +
  facet_wrap(~celltype, scales = 'free', nrow=2) + stat_cor(size=3.5) + theme_classic() + 
  theme(legend.position = 'none', panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.text = element_text(size=12), axis.title = element_text(size=13),
        axis.text=element_text(size=10)) + 
  xlab('M1 summary score at d2 (norm. to baseline)') + ylab('Number of cells in BAL (log10)')

Ab_modscores <- merge.data.table(Ab_dt, modscores_norm[,.(id,TimeafterBCG,score,moduleID,norm_score)], by='id', allow.cartesian=T)
ggplot(Ab_modscores[time == 4 & moduleID == 'M1' & TimeafterBCG %in% c('d2') & grepl('BAL',type)], aes(x=norm_score, y=value)) + 
  geom_point(aes(color=moduleID)) + 
  scale_color_locuszoom(alpha=0.6) +
  facet_wrap(~type, scales = 'free') + stat_cor(size=3.5) + theme_classic() + 
  theme(legend.position = 'none', panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.text = element_text(size=12), axis.title = element_text(size=13),
        axis.text.y=element_text(size=10), axis.text.x=element_text(size=10, angle=40, hjust=1)) + 
  xlab('M1 summary score at d2 (norm. to baseline)') + ylab('Log10 Titer')
