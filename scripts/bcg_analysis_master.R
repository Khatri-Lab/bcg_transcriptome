library(data.table)
library(DESeq2)
library(magrittr)
library(pheatmap)
library(ggsci)
library(apeglm)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(clusterProfiler)
library(WGCNA)
library(enrichplot)
library(org.Mmu.eg.db)
library(org.Hs.eg.db)
library(pROC)
library(ggROC)

source('/labs/khatrilab/yiranliu/BCG/github/scripts/deg_functions.R')

mmulgenes = fread('/labs/khatrilab/yiranliu/BCG/github/mmul_ref_objects/mmul10_genes.txt')
mmulgenes_only <- unique(mmulgenes[,.(`Gene stable ID`, `Gene name`)])

bcg1pheno <- fread('/labs/khatrilab/yiranliu/BCG/github/data/bcg_routecohort_metadata.csv')
bcg2pheno <- fread('/labs/khatrilab/yiranliu/BCG/github/data/ivbcg_dosecohort_metadata.csv', header = T)

setnames(bcg1pheno, c('description','animal_id','dose_group','protect_outcome'), c('SampleID','id','dose','protect_status'))
setnames(bcg2pheno, c('description','animal_id','dose_group','protect_outcome'), c('SampleID','id','dose','protect_status'))

bcg1pheno$TimeafterBCG <- factor(bcg1pheno$TimeafterBCG, levels = c('pre','d2','wk2','wk12'))
bcg2pheno$TimeafterBCG <- factor(bcg2pheno$TimeafterBCG, levels = c('pre','d2','wk2','wk4','wk12'))
bcg1pheno$dose <- factor(bcg1pheno$dose, levels = c('low','high'))
bcg2pheno$dose <- factor(bcg2pheno$dose, levels = c('low','high'))

modGenes_mmul <- readRDS('/labs/khatrilab/yiranliu/BCG/github/mmul_ref_objects/modGenes_mmul.RDS') # BTMs converted to mmul IDs
btm_dt <- data.table(data.frame(gene=unlist(modGenes_mmul)), keep.rownames = 'btm')

bcg1pheno_df <- data.frame(bcg1pheno, row.names = bcg1pheno$SampleID)
bcg2pheno_df <- data.frame(bcg2pheno, row.names = bcg2pheno$SampleID)

# these processed files can be downloaded from GEO GSE218157 and GSE218270
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

# univariate logistic regression for protection outcome by dose
bcg2pheno_withimputed[,class := ifelse(protect_status == 'protected', 1, 0)]
output <- glm(class~BCG_dose_log10, family="binomial", data=bcg2pheno_withimputed[TimeafterBCG == 'wk12'])
summary(output)

# annotation colors for heatmaps
dosecolors <- brewer.pal(8,"Blues")[c(2,5,8)]
names(dosecolors) <- levels(bcg2pheno$dose)

timecolors <- brewer.pal(length(c('pre','d2','wk2','wk4','wk12')),"Set3")
names(timecolors) <- c('pre','d2','wk2','wk4','wk12')

protcolors <- pal_jama()(n=2)
names(protcolors) <- c("not_protected", "protected")

anno_colors <- list(dose = dosecolors, TimeafterBCG = timecolors, protect_status = protcolors)

# Set up annotation df for heatmaps
df <- bcg2pheno[TimeafterBCG %in% c('pre','d2','wk2','wk4','wk12') & dose == 'high' & SampleID %in% colnames(vsdmat_bcg2)]
df <- df[order(TimeafterBCG, BCG_dose_log10, protect_status),c("SampleID","TimeafterBCG")]

breaks <- c()
n_vect <- df[, .N, by='TimeafterBCG']$N
for (i in 1:length(n_vect)-1){
  b <- sum(n_vect[1:i])
  breaks <- c(breaks, b)
}
df <- data.frame(df[,-1], row.names=df$SampleID)
df <- droplevels(df)

#### IV high BCG DEGs ####
# present_txi <- txi_bcg2
# present_design <- ~id + TimeafterBCG
# 
# iv_hi_bcg_dds <- getdds(bcg2pheno[TimeafterBCG %in% c('pre','d2','wk2','wk4','wk12') & dose == 'high', SampleID],
#                         txi=present_txi, pheno=bcg2pheno,
#                         design = present_design)
# resultsNames(iv_hi_bcg_dds)
# 
# iv_hi_bcg_d2_res <- getRes(iv_hi_bcg_dds, "TimeafterBCG_d2_vs_pre")
# iv_hi_bcg_wk2_res <- getRes(iv_hi_bcg_dds, "TimeafterBCG_wk2_vs_pre")
# iv_hi_bcg_wk4_res <- getRes(iv_hi_bcg_dds, "TimeafterBCG_wk4_vs_pre")
# iv_hi_bcg_wk12_res <- getRes(iv_hi_bcg_dds, "TimeafterBCG_wk12_vs_pre")

# read in saved deseq2 objects (dds object too large for github)
iv_hi_bcg_d2_res <- readRDS('/labs/khatrilab/yiranliu/BCG/github/deseq_objects/iv_hi_bcg_d2_res_dosestudyonly.RDS')
iv_hi_bcg_wk2_res <- readRDS('/labs/khatrilab/yiranliu/BCG/github/deseq_objects/iv_hi_bcg_wk2_res_dosestudyonly.RDS')
iv_hi_bcg_wk4_res <- readRDS('/labs/khatrilab/yiranliu/BCG/github/deseq_objects/iv_hi_bcg_wk4_res_dosestudyonly.RDS')
iv_hi_bcg_wk12_res <- readRDS('/labs/khatrilab/yiranliu/BCG/github/deseq_objects/iv_hi_bcg_wk12_res_dosestudyonly.RDS')

# d2
iv_hi_bcg_d2_sig_pos <- getsiggenes(iv_hi_bcg_d2_res, lfc_cutoff = log(1.5, base=2), pos_only = T)
iv_hi_bcg_d2_sig_neg <- getsiggenes(iv_hi_bcg_d2_res, lfc_cutoff = log(1.5, base=2), neg_only = T)
print(paste0('Up genes: ', length(iv_hi_bcg_d2_sig_pos), '; Down genes: ', length(iv_hi_bcg_d2_sig_neg)))

# wk2
iv_hi_bcg_wk2_sig_pos <- getsiggenes(iv_hi_bcg_wk2_res, lfc_cutoff = log(1.5, base=2), pos_only = T)
iv_hi_bcg_wk2_sig_neg <- getsiggenes(iv_hi_bcg_wk2_res, lfc_cutoff = log(1.5, base=2), neg_only = T)
print(paste0('Up genes: ', length(iv_hi_bcg_wk2_sig_pos), '; Down genes: ', length(iv_hi_bcg_wk2_sig_neg)))

# wk4
iv_hi_bcg_wk4_sig_pos <- getsiggenes(iv_hi_bcg_wk4_res, lfc_cutoff = log(1.5, base=2), pos_only = T)
iv_hi_bcg_wk4_sig_neg <- getsiggenes(iv_hi_bcg_wk4_res, lfc_cutoff = log(1.5, base=2), neg_only = T)
print(paste0('Up genes: ', length(iv_hi_bcg_wk4_sig_pos), '; Down genes: ', length(iv_hi_bcg_wk4_sig_neg)))

# wk12
iv_hi_bcg_wk12_sig_pos <- getsiggenes(iv_hi_bcg_wk12_res, lfc_cutoff = log(1.5, base=2), pos_only = T)
iv_hi_bcg_wk12_sig_neg <- getsiggenes(iv_hi_bcg_wk12_res, lfc_cutoff = log(1.5, base=2), neg_only = T)
print(paste0('Up genes: ', length(iv_hi_bcg_wk12_sig_pos), '; Down genes: ', length(iv_hi_bcg_wk12_sig_neg)))

# DEGs in IV high BCG
n_degs <- data.table(time=c(rep('d2',2),rep('wk2',2),rep('wk4',2),rep('wk12',2)),
                     upordown=c(rep(c('upregulated','downregulated'),4)),
                     N=c(length(iv_hi_bcg_d2_sig_pos), length(iv_hi_bcg_d2_sig_neg),
                         length(iv_hi_bcg_wk2_sig_pos), length(iv_hi_bcg_wk2_sig_neg),
                         length(iv_hi_bcg_wk4_sig_pos), length(iv_hi_bcg_wk4_sig_neg),
                         length(iv_hi_bcg_wk12_sig_pos), length(iv_hi_bcg_wk12_sig_neg)))

n_degs$time <- factor(n_degs$time, levels = c('d2','wk2','wk4','wk12'))
ggplot(n_degs, aes(x=time, y=N, fill=upordown)) + geom_col() + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1,size=12),
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=12),
        legend.text = element_text(size=12)) + scale_fill_lancet() + 
  labs(fill='', x='Time after BCG', y='Number of DEGs')

length(unique(c(iv_hi_bcg_d2_sig_pos,iv_hi_bcg_d2_sig_neg,iv_hi_bcg_wk2_sig_pos,iv_hi_bcg_wk2_sig_neg,
                iv_hi_bcg_wk4_sig_pos,iv_hi_bcg_wk4_sig_neg,
                iv_hi_bcg_wk12_sig_pos, iv_hi_bcg_wk12_sig_neg)))

present_genes <- intersect(unique(c(iv_hi_bcg_d2_sig_pos, iv_hi_bcg_d2_sig_neg,
                                    iv_hi_bcg_wk2_sig_pos, iv_hi_bcg_wk2_sig_neg,
                                    iv_hi_bcg_wk4_sig_pos, iv_hi_bcg_wk4_sig_neg,
                                    iv_hi_bcg_wk12_sig_pos, iv_hi_bcg_wk12_sig_neg)), rownames(vsdmat_bcg2))
present_samps <- intersect(bcg2pheno[TimeafterBCG %in% c('pre','d2','wk2','wk4','wk12') & dose == 'high', SampleID],
                           colnames(vsdmat_bcg2_withimputed))

##### WGCNA #####
wgcna_input <- t(vsdmat_bcg2_withimputed[present_genes,present_samps])
check <- goodSamplesGenes(wgcna_input, verbose=3) # check for missing values and/or zero variance in samples and genes
check$allOK

# sample clustering to detect outliers
sampleTree = hclust(dist(wgcna_input), method = "average")
sizeGrWindow(6,12)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2) # none

set.seed(0)
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(wgcna_input, powerVector = powers, verbose = 5)

# Scale-free topology fit index as a function of the soft-thresholding power
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net = blockwiseModules(wgcna_input, power = 9, minModuleSize = 50,
                       mergeCutHeight = 0.2, networkType = 'signed',
                       numericLabels = TRUE, saveTOMs = T, 
                       saveTOMFileBase = '/labs/khatrilab/yiranliu/BCG/github/wgcna/TOM',
                       verbose = 3, randomSeed=2021)
# saveRDS('/labs/khatrilab/yiranliu/BCG/github/wgcna/iv_hi_wgcna.RDS')

table(net$colors)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = data.frame(net$colors)
moduleLabels <- moduleLabels[moduleLabels$net.colors != 0,,drop=F]
moduleLabels$net.colors <- as.character(moduleLabels$net.colors)

# named/ordered based on trends over time
moduleLabels[moduleLabels$net.colors == 1, 'moduleID'] <- 'M1'
moduleLabels[moduleLabels$net.colors == 2, 'moduleID'] <- 'M2'
moduleLabels[moduleLabels$net.colors == 6, 'moduleID'] <- 'M3'
moduleLabels[moduleLabels$net.colors == 4, 'moduleID'] <- 'M4'
moduleLabels[moduleLabels$net.colors == 3, 'moduleID'] <- 'M5'
moduleLabels[moduleLabels$net.colors == 7, 'moduleID'] <- 'M6'
moduleLabels[moduleLabels$net.colors == 5, 'moduleID'] <- 'M7'

moduleLabels$moduleID <- factor(moduleLabels$moduleID, levels = sort(unique(moduleLabels$moduleID)))
moduleLabels <- moduleLabels[order(moduleLabels$moduleID),]

geneTree = net$dendrograms[[1]]
moduleLabels_wnames <- rbind(as.data.table(moduleLabels[,2,drop=F], keep.rownames = 'Gene stable ID'),
                             data.table(`Gene stable ID`=names(net$colors[which(net$colors == 0)]),
                                        moduleID='N/A'))
moduleLabels_wnames <- merge(moduleLabels_wnames, 
                             mmulgenes_only, by='Gene stable ID')

modcolors <- get_palette('locuszoom', length(unique(moduleLabels$moduleID))+1)[1:7]
names(modcolors) <- sort(unique(moduleLabels$moduleID))

anno_colors <- list(dose = dosecolors, dosegroup = dosegroup_colors, 
                    TimeafterBCG = timecolors,
                    protect_status = protcolors, batch=batchcolors,
                    moduleID = modcolors)

pheatmap(vsdmat_bcg2[rownames(moduleLabels),
                     rownames(df)], annotation_col=df, annotation_row = moduleLabels[,2,drop=F], 
         annotation_colors = anno_colors, color = colorRampPalette(c("blue","white","red"))(100), 
         show_rownames = F, border_color = NA, scale = 'row', breaks=seq(-4,4,8/100), cluster_cols = F,
         cluster_rows = F, 
         show_colnames = F, gaps_col = breaks[2:length(breaks)])

#### BTM analysis #####
M1_btms <- get_BTM_dt(iv_hi_bcg_d2_res, sig_genes = rownames(moduleLabels[moduleLabels$moduleID == 'M1',]))
head(M1_btms@enrichmentTable)

M2_btms <- get_BTM_dt(iv_hi_bcg_d2_res, sig_genes = rownames(moduleLabels[moduleLabels$moduleID == 'M2',]))
head(M2_btms@enrichmentTable)

M3_btms <- get_BTM_dt(iv_hi_bcg_d2_res, sig_genes = rownames(moduleLabels[moduleLabels$moduleID == 'M3',]))
head(M3_btms@enrichmentTable)

M4_btms <- get_BTM_dt(iv_hi_bcg_d2_res, sig_genes = rownames(moduleLabels[moduleLabels$moduleID == 'M4',]))
head(M4_btms@enrichmentTable)

M5_btms <- get_BTM_dt(iv_hi_bcg_d2_res, sig_genes = rownames(moduleLabels[moduleLabels$moduleID == 'M5',]))
head(M5_btms@enrichmentTable)

M6_btms <- get_BTM_dt(iv_hi_bcg_d2_res, sig_genes = rownames(moduleLabels[moduleLabels$moduleID == 'M6',]))
head(M6_btms@enrichmentTable)

M7_btms <- get_BTM_dt(iv_hi_bcg_d2_res, sig_genes = rownames(moduleLabels[moduleLabels$moduleID == 'M7',]))
head(M7_btms@enrichmentTable)


iv_hi_btms <- merge_BTMs(list(M1_btms, M2_btms, M3_btms, M4_btms, M5_btms, M6_btms, M7_btms), 
                         suffx = c('M1','M2','M3','M4','M5','M6','M7'), padj_cutoff = 0.01)

# top BTMs by group
setorder(iv_hi_btms, variable, padj)

to_plot_setnames <- as.character(na.omit(iv_hi_btms[!is.na(padj),.SD[1:6], by=variable]$set.name))
iv_hi_btms[,set.name.plot := tstrsplit(set.name, ' (M', fixed=T, keep=1)]
setnames_order <- unique(iv_hi_btms[set.name %in% to_plot_setnames & padj <= 0.01,set.name.plot])
iv_hi_btms$set.name.plot <- factor(iv_hi_btms$set.name.plot, levels = setnames_order)
ggplot(na.omit(iv_hi_btms[padj <= 0.01,.SD[1:6], by=variable]), aes(variable, set.name.plot, color=variable)) + 
  geom_point(aes(size=p_inv)) + ylab("") + xlab("") + 
  theme_bw() + scale_color_manual(values = get_palette('locuszoom', length(unique(full_dt$moduleID))+1)) +
  scale_y_discrete(limits = rev) +
  scale_size_continuous(breaks=c(2.5,5,10,20,50), trans = 'sqrt') +
  labs(size='-log10(p.adj)', color='Module') + 
  theme(axis.text.y = element_text(size=12), axis.text.x=element_blank(),
        axis.ticks.x = element_blank(), legend.text = element_text(size=12),
        legend.title = element_text(size=14))

# Nothing for M7 - try looking at GO terms
iv_hi_geneuniverse <- mmulgenes_only[`Gene stable ID` %in% rownames(iv_hi_bcg_d2_res), `Gene name`]

M7_genenames <- mmulgenes_only[`Gene stable ID` %in% rownames(moduleLabels[moduleLabels$moduleID == 'M7',]), `Gene name`]
M7_GO <- enrichGO(gene=M7_genenames,
                  universe=iv_hi_geneuniverse,
                  OrgDb=org.Hs.eg.db, keyType = 'SYMBOL',
                  ont = 'BP', pvalueCutoff = 0.05)
head(M7_GO@result,20) # nothing with padj < 0.05

##### Mod scores over time #####
full_dt <- data.table()
for (i in unique(moduleLabels$moduleID)){
  genes <- rownames(subset(moduleLabels, moduleID == i))
  score_dt <- get_score_pheno(up_geneIDs = genes, down_geneIDs = c(), genematx = vsdmat_bcg2_withimputed, pheno_dt = bcg2pheno_withimputed)
  score_dt$moduleID <- i
  full_dt <- rbind(full_dt, score_dt)
}
full_dt[,norm_score := score - score[TimeafterBCG == 'pre'], by=c('id','moduleID')]
full_dt$dose <- factor(full_dt$dose, levels = c('high','low'))
full_dt$TimeafterBCG <- factor(full_dt$TimeafterBCG, levels = c(c('pre','d2','wk2','wk4','wk12')))

# line plots with median
medians <- full_dt[,.(med_score = median(score),
                      med_normscore = median(norm_score)), by=c('dose','TimeafterBCG','moduleID')]
# high dose
ggplot(full_dt[TimeafterBCG %in% c('pre','d2','wk2','wk4','wk12') & dose == 'high']) + 
  geom_line(aes(x=TimeafterBCG, y=score, group=id), color='light grey') +
  geom_line(data = medians[TimeafterBCG %in% c('pre','d2','wk2','wk4','wk12') & dose == 'high'], 
            aes(x=TimeafterBCG, y=med_score, color=moduleID, group=dose), size=1.5) + 
  scale_color_manual(values = get_palette('locuszoom', length(unique(full_dt$moduleID))+1)) +
  facet_wrap(~moduleID, scales = 'free', nrow=1) + ylab('Summary Score') + theme_classic() +
  theme(axis.text.x=element_text(angle = 60, hjust=1, size=14), axis.text.y=element_text(size=14),
        axis.title = element_text(size=16), strip.text=element_text(size=16), legend.position = 'none')

# low dose
ggplot(full_dt[TimeafterBCG %in% c('pre','d2','wk2','wk4','wk12') & dose == 'low']) + 
  geom_line(aes(x=TimeafterBCG, y=score, group=id), color='light grey') +
  geom_line(data = medians[TimeafterBCG %in% c('pre','d2','wk2','wk4','wk12') & dose == 'low'], 
            aes(x=TimeafterBCG, y=med_score, color=moduleID, group=dose), size=1.5) + 
  scale_color_manual(values = get_palette('locuszoom', length(unique(full_dt$moduleID))+1)) +
  facet_wrap(~moduleID, scales = 'free', nrow=1) + ylab('Summary Score') + theme_classic() +
  theme(axis.text.x=element_text(angle = 60, hjust=1, size=14), axis.text.y=element_text(size=14),
        axis.title = element_text(size=16), strip.text=element_text(size=16), legend.position = 'none')

# boxplots high and low by protection outcome
ggplot(full_dt[TimeafterBCG %in% c('d2','wk2','wk4','wk12')], aes(x=TimeafterBCG, y=norm_score, color=protect_status)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), size = 0.5) + 
  geom_point(data=full_dt[TimeafterBCG %in% c('d2','wk2','wk4','wk12')], 
             aes(x=TimeafterBCG, y=norm_score+0.2), color='NA') +
  theme_bw() + 
  facet_grid(moduleID~dose, scales = 'free') + ylab('Summary score \n(normalized to baseline)') +
  scale_color_jama() + xlab('Time after BCG') + labs(color='Challenge outcome') +
  stat_compare_means(label.y.npc = 0.94, 
                     label = 'p.signif', method = 'wilcox.test', size = 2.5) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey')

# just d2 boxplots
ggplot(full_dt[TimeafterBCG %in% c('d2') & moduleID == 'M1'], aes(x=TimeafterBCG, y=norm_score, color=protect_status)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_text(data=data.table(TimeafterBCG='d2', norm_score=1.75, dose='high'), label = "", size=4, color='black') +
  geom_text(data=data.table(TimeafterBCG='d2', norm_score=1.15, dose='low'), label = "", size=4, color='black') +
  facet_wrap(~dose, scales = 'free') +
  ylab('M1 summary score \n(norm. to baseline)') + xlab('Time after BCG') + scale_color_jama() + labs(color='') +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0), size = 1) + theme_classic() +
  theme(legend.position = 'bottom', axis.title=element_text(size=13), axis.text = element_text(size=12),
        legend.text = element_text(size=12), strip.text = element_text(size=12)) +
  stat_compare_means(method = 'wilcox', method.args = list(alternative = "less"), label.y.npc = 0.97, size=3)

# M1 individual BTMs on d2
sig_d2_btms <- M1_btms@enrichmentTable[p.adj <= 0.01, set.name]

all_out <- data.table()
for (i in 1:length(M1_btms@sigGenesInSets[sig_d2_btms])){
  btm <- names(M1_btms@sigGenesInSets)[i]
  if (grepl('TBA',btm)) {
    next
  }
  score_dt <- get_score_pheno(up_geneIDs = names(M1_btms@sigGenesInSets[[i]]), down_geneIDs = c(), 
                              genematx = vsdmat_bcg2_withimputed, pheno_dt = bcg2pheno_withimputed[TimeafterBCG %in% c('pre','d2')])
  score_dt[,norm_score := score - score[TimeafterBCG == 'pre'], by=c('id')]
  score_dt$btm <- btm
  all_out <- rbind(all_out, score_dt, fill=T)
}

table(all_out$btm)
M1_score <- get_score_pheno(up_geneIDs = rownames(moduleLabels[moduleLabels$moduleID == 'M1',]), down_geneIDs = c(), 
                            genematx = vsdmat_bcg2_withimputed, pheno_dt = bcg2pheno_withimputed[TimeafterBCG %in% c('pre','d2')])
M1_score[,norm_score := score - score[TimeafterBCG == 'pre'], by=c('id')]

M1_score$btm <- 'All M1 genes'

all_out <- rbind(all_out, M1_score)
all_out[,class := ifelse(protect_status == 'protected', 1, 0)]

hi_auc_vect <- c()
hi_ci_lo_vect <- c()
hi_ci_up_vect <- c()

lo_auc_vect <- c()
lo_ci_lo_vect <- c()
lo_ci_up_vect <- c()

both_auc_vect <- c()
both_ci_lo_vect <- c()
both_ci_up_vect <- c()

for (btm_current in unique(all_out$btm)){
  # print(btm_current)
  hi_roc_current <- roc(all_out[btm == btm_current & TimeafterBCG == 'd2' & dose == 'high', class], 
                        all_out[btm == btm_current & TimeafterBCG == 'd2' & dose == 'high', norm_score])
  hi_auc_vect <- c(hi_auc_vect, auc(hi_roc_current))
  hi_ci_current <- ci.auc(hi_roc_current)
  hi_ci_lo_vect <- c(hi_ci_lo_vect, hi_ci_current[1])
  hi_ci_up_vect <- c(hi_ci_up_vect, hi_ci_current[3])
  
  lo_roc_current <- roc(all_out[btm == btm_current & TimeafterBCG == 'd2' & dose == 'low', class], 
                        all_out[btm == btm_current & TimeafterBCG == 'd2' & dose == 'low', norm_score])
  lo_auc_vect <- c(lo_auc_vect, auc(lo_roc_current))
  lo_ci_current <- ci.auc(lo_roc_current)
  lo_ci_lo_vect <- c(lo_ci_lo_vect, lo_ci_current[1])
  lo_ci_up_vect <- c(lo_ci_up_vect, lo_ci_current[3])
  
  both_roc_current <- roc(all_out[btm == btm_current & TimeafterBCG == 'd2', class], 
                          all_out[btm == btm_current & TimeafterBCG == 'd2', norm_score])
  both_auc_vect <- c(both_auc_vect, auc(both_roc_current))
  both_ci_current <- ci.auc(both_roc_current)
  both_ci_lo_vect <- c(both_ci_lo_vect, both_ci_current[1])
  both_ci_up_vect <- c(both_ci_up_vect, both_ci_current[3])
}

auc_dt <- data.table(btm=rep(unique(all_out$btm),3), auc=c(hi_auc_vect,lo_auc_vect,both_auc_vect), 
                     ci_lo=c(hi_ci_lo_vect,lo_ci_lo_vect,both_ci_lo_vect), 
                     ci_up=c(hi_ci_up_vect, lo_ci_up_vect,both_ci_up_vect),
                     dose=c(rep('high', length(unique(all_out$btm))),rep('low', length(unique(all_out$btm))),
                            rep('combined', length(unique(all_out$btm)))))
auc_dt$dose <- factor(auc_dt$dose, levels = c('high','low','combined'))
setorder(auc_dt, dose, -auc)
auc_dt$btm <- factor(auc_dt$btm, levels = unique(auc_dt$btm))
auc_dt$btm.plot <- auc_dt$btm
auc_dt[grepl('(M', btm, fixed=T),btm.plot := tstrsplit(btm, ' (M', fixed=T, keep=1)]
auc_dt[grepl('(S', btm, fixed=T),btm.plot := tstrsplit(btm, ' (S', fixed=T, keep=1)]
auc_dt$btm.plot <- factor(auc_dt$btm.plot, levels = unique(auc_dt$btm.plot))
ggplot(auc_dt[dose != 'combined'], aes(x=btm.plot, y=auc)) + geom_col(fill='grey') + 
  scale_x_discrete(limits=rev) + theme_bw() + xlab('') + ylab('AUROC') + 
  coord_flip() +
  theme(plot.margin = margin(1,1,1,1,'cm'), axis.text.y = element_text(size=14),
        axis.text.x=element_text(size=12, angle=35, hjust=1),
        strip.text = element_text(size=14), axis.title=element_text(size=14)) + facet_wrap(~dose)

to_plot_btms <- c('All M1 genes', as.character(auc_dt[dose == 'high' & auc >= 0.9, btm]))

all_out$btm <- factor(all_out$btm, levels = c('All M1 genes', setdiff(unique(auc_dt$btm),'All M1 genes')))

# Non-linear dose-response model regressing d2 scores on BCG dose & protection outcome
full_dt[,class := ifelse(protect_status == 'protected', 1, 0)]
set.seed(0)
nls.fit <- nls(norm_score ~ (A + beta_A*class) + (K - (A + beta_A*class))/(1+exp((-(B)*(BCG_dose_log10-M)))),
               data=full_dt[moduleID == 'M1' & TimeafterBCG == 'd2'],
               start = list(B=3,A=0.2,beta_A=0.6,K=1.4,M=6.4),
               control=list(maxiter=1000))
summary(nls.fit)
AIC(nls.fit)

fitdt <- data.table(BCG_dose_log10=c(seq(min(bcg2pheno[dose != 'unvax' & protect_status=='not_protected', BCG_dose_log10])-0.05,
                                  max(bcg2pheno[dose != 'unvax' & protect_status=='not_protected', BCG_dose_log10])+0.05, length.out=100),
                              seq(min(bcg2pheno[dose != 'unvax' & protect_status=='protected', BCG_dose_log10]-0.05),
                                  max(bcg2pheno[dose != 'unvax' & protect_status=='protected', BCG_dose_log10]+0.05), length.out=100)),
                    class=rep(c(0,1), each=100))
fitdt$preds <- predict(nls.fit, newdata = fitdt)

# Correlation plots for d2 btms and mtb CFUs
all_out$btm.plot <- all_out$btm
all_out[grepl('(M', btm, fixed=T),btm.plot := tstrsplit(btm, ' (M', fixed=T, keep=1)]
all_out[grepl('(S', btm, fixed=T),btm.plot := tstrsplit(btm, ' (S', fixed=T, keep=1)]
all_out[grepl('Activated (LPS)', btm, fixed=T), btm.plot := 'Activated (LPS) DC surface sig (9)']
all_out$btm.plot <- factor(all_out$btm.plot, levels = c('All M1 genes', 'RIG-1 like receptor signaling',
                                                        'innate antiviral response', setdiff(unique(all_out$btm.plot), 
                                                                                             c('All M1 genes',
                                                                                               'RIG-1 like receptor signaling',
                                                                                               'innate antiviral response'))),
                           labels = c('All M1 genes (698)', 'RIG-1 like receptor signaling (8)',
                                      'innate antiviral response (9)', setdiff(unique(all_out$btm.plot), 
                                                                               c('All M1 genes',
                                                                                 'RIG-1 like receptor signaling',
                                                                                 'innate antiviral response'))))
ggplot(all_out[btm %in% to_plot_btms & TimeafterBCG == 'd2'], aes(x=norm_score, y=log_Mtb_CFU)) +
  annotate(geom = 'rect', xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0.4, fill=alpha('grey',0.5)) +
  geom_jitter(width=0, height=0.2, size=1.5) + 
  facet_wrap(~btm.plot, ncol=2, scales='free') + xlab('D2 Score \n(Norm. to Baseline)') + 
  ylab(expression(paste("Total ", italic("Mtb "), "CFU (log10)"))) + labs(color='Dose') + theme_classic() + 
  scale_color_lancet() + theme(axis.text.x = element_text(angle = 45, hjust=1, size=12),
                               axis.text.y=element_text(size=12),
                               axis.title=element_text(size=12),
                               strip.text = element_text(size=9)) +
  stat_cor(size=3) + ylim(-0.2,8)

ggplot(all_out[btm %in% to_plot_btms & TimeafterBCG == 'd2'], aes(x=norm_score, y=log_grans_nx)) +
  annotate(geom = 'rect', xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0.12, fill=alpha('grey',0.5)) +
  geom_jitter(width=0, height=0.05, size=1.5) + 
  facet_wrap(~btm.plot, ncol=2, scales='free') + xlab('D2 Score \n(Norm. to Baseline)') + 
  ylab('# of granulomas (log10)') + labs(color='Dose') + theme_classic() + 
  scale_color_lancet() + theme(axis.text.x = element_text(angle = 45, hjust=1, size=12),
                               axis.text.y=element_text(size=12),
                               axis.title=element_text(size=12),
                               strip.text = element_text(size=9)) +
  stat_cor(size=3) + ylim(-.05,2.5)

# ROC curves
to_plot_btms <- c(as.character(auc_dt[dose == 'high' & auc >= 0.9, btm]), "All M1 genes")
btm_names_short <- c('RIG-1 like receptor signaling', 'Innate antiviral response', 'Activated (LPS) DC surface sig', 'All M1 genes')
roc_list <- vector('list', length=length(to_plot_btms) * 2)
names_vect <- c()

for (i in 1:length(to_plot_btms)){
  btm_current <- to_plot_btms[i]
  btm_short <- btm_names_short[i]
  hi_roc <- roc(all_out[btm == btm_current & TimeafterBCG == 'd2' & dose == 'high', class], 
                all_out[btm == btm_current & TimeafterBCG == 'd2' & dose == 'high', norm_score])
  hi_auc <- auc(hi_roc)
  hi_ci <- ci.auc(hi_roc)
  hi_txt <- paste(btm_short, ' (High): AUC=', round(hi_auc,2), ' (95% CI ', round(hi_ci[1],2), '-', round(hi_ci[3],2), ')', sep='')
  
  lo_roc <- roc(all_out[btm == btm_current & TimeafterBCG == 'd2' & dose == 'low', class], 
                all_out[btm == btm_current & TimeafterBCG == 'd2' & dose == 'low', norm_score])
  lo_auc <- auc(lo_roc)
  lo_ci <- ci.auc(lo_roc)
  lo_txt <- paste(btm_short, ' (Low): AUC=', round(lo_auc,2), ' (95% CI ', round(lo_ci[1],2), '-', round(lo_ci[3],2), ')', sep='')
  roc_list[[i*2-1]] <- hi_roc
  roc_list[[i*2]] <- lo_roc
  names_vect <- c(names_vect, hi_txt, lo_txt)
}
names(roc_list) <- names_vect

pROC::ggroc(roc_list, legacy.axes = T, aes=c('color','linetype')) + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dotted") + 
  ggtitle("") + scale_color_nejm() +
  scale_color_manual(values = rep(get_palette('nejm',4),each=2)) + 
  scale_linetype_manual(values = rep(c('dashed','solid'),4)) +
  theme_bw() + theme(legend.title = element_blank(),
                     legend.background = element_blank(), legend.text = element_text(size = 10),
                     axis.title = element_text(size=13),
                     axis.text = element_text(size=10))

### ROUTE cohort ###
full_dt_bcg1 <- data.table()
for (i in unique(moduleLabels$moduleID)){
  genes <- rownames(subset(moduleLabels, moduleID == i))
  score_dt <- get_score_pheno(up_geneIDs = genes, down_geneIDs = c(), genematx = vsdmat_bcg1, pheno_dt = bcg1pheno)
  score_dt$moduleID <- i
  full_dt_bcg1 <- rbind(full_dt_bcg1, score_dt)
}
full_dt_bcg1[,norm_score := score - score[TimeafterBCG == 'pre'], by=c('id','moduleID')]
full_dt_bcg1$TimeafterBCG <- factor(full_dt_bcg1$TimeafterBCG, levels = c(c('pre','d2','wk2','wk12')))
# write.csv(full_dt_bcg1, '/labs/khatrilab/yiranliu/BCG/github/wgcna/modScores_routecohort.csv', row.names = F)

medians_bcg1 <- full_dt_bcg1[,.(med_score = median(score),
                                med_normscore = median(norm_score)), by=c('vax_group','TimeafterBCG','moduleID')]
ggplot(full_dt_bcg1[TimeafterBCG %in% c('pre','d2','wk2','wk12') & vax_group == 'IV']) + 
  geom_line(aes(x=TimeafterBCG, y=score, group=id), color='light grey') +
  geom_line(data = medians_bcg1[TimeafterBCG %in% c('pre','d2','wk2','wk12') & vax_group == 'IV'], 
            aes(x=TimeafterBCG, y=med_score, color=moduleID, group=vax_group), size=1.5) + 
  scale_color_manual(values = get_palette('locuszoom', length(unique(full_dt$moduleID))+1)) +
  facet_wrap(~moduleID, scales = 'free', nrow=1) + ylab('Summary Score') + theme_classic() + xlab('Time after BCG') +
  theme(axis.text.x=element_text(angle = 60, hjust=1, size=14), axis.text.y=element_text(size=14),
        axis.title = element_text(size=16), strip.text=element_text(size=16), legend.position = 'none')

ggplot(full_dt_bcg1[TimeafterBCG %in% c('pre','d2','wk2','wk12')]) + 
  geom_line(aes(x=TimeafterBCG, y=score, group=id), color='light grey') +
  geom_line(data = medians_bcg1[TimeafterBCG %in% c('pre','d2','wk2','wk12')], 
            aes(x=TimeafterBCG, y=med_score, color=moduleID, group=vax_group), size=1.5) + 
  scale_color_manual(values = get_palette('locuszoom', length(unique(full_dt$moduleID))+1)) +
  facet_grid(moduleID~vax_group, scales = 'free') + ylab('Summary Score') + theme_classic() + xlab('Time after BCG') +
  theme(axis.text.x=element_text(angle = 60, hjust=1, size=14), axis.text.y=element_text(size=14),
        axis.title = element_text(size=16), strip.text=element_text(size=16), legend.position = 'none',
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))

all_out_route <- data.table()
for (i in 1:length(M1_btms@sigGenesInSets[1:20])){
  btm <- names(M1_btms@sigGenesInSets)[i]
  if (grepl('TBA',btm)) {
    next
  }
  score_dt <- get_score_pheno(up_geneIDs = names(M1_btms@sigGenesInSets[[i]]), down_geneIDs = c(), 
                              genematx = vsdmat_bcg1, pheno_dt = bcg1pheno[TimeafterBCG %in% c('pre','d2')])
  score_dt[,norm_score := score - score[TimeafterBCG == 'pre'], by=c('id')]
  score_dt$btm <- btm
  all_out_route <- rbind(all_out_route, score_dt, fill=T)
}

M1_score <- get_score_pheno(up_geneIDs = rownames(moduleLabels[moduleLabels$moduleID == 'M1',]), down_geneIDs = c(), 
                            genematx = vsdmat_bcg1, pheno_dt = bcg1pheno[TimeafterBCG %in% c('pre','d2')])
M1_score[,norm_score := score - score[TimeafterBCG == 'pre'], by=c('id')]

M1_score$btm <- 'All M1 genes'

all_out_route <- rbind(all_out_route, M1_score)
all_out_route[,class := ifelse(protect_status == 'protected', 1, 0)]

all_out_route$btm.plot <- all_out_route$btm
all_out_route[grepl('(M', btm, fixed=T),btm.plot := tstrsplit(btm, ' (M', fixed=T, keep=1)]
all_out_route[grepl('(S', btm, fixed=T),btm.plot := tstrsplit(btm, ' (S', fixed=T, keep=1)]
all_out_route[grepl('Activated (LPS)', btm, fixed=T), btm.plot := 'Activated (LPS) DC surface sig (9)']
all_out_route$btm.plot <- factor(all_out_route$btm.plot, levels = c('All M1 genes', 'RIG-1 like receptor signaling',
                                                                    'innate antiviral response', setdiff(unique(all_out_route$btm.plot), 
                                                                                                         c('All M1 genes',
                                                                                                           'RIG-1 like receptor signaling',
                                                                                                           'innate antiviral response'))),
                                 labels = c('All M1 genes (698)', 'RIG-1 like receptor signaling (8)',
                                            'innate antiviral response (9)', setdiff(unique(all_out_route$btm.plot), 
                                                                                     c('All M1 genes',
                                                                                       'RIG-1 like receptor signaling',
                                                                                       'innate antiviral response'))))
ggplot(all_out_route[btm %in% to_plot_btms & TimeafterBCG == 'd2'], aes(x=norm_score, y=log_Mtb_CFU)) +
  annotate(geom = 'rect', xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0.4, fill=alpha('grey',0.5)) +
  geom_jitter(aes(color=vax_group), width=0, height = 0.2, size=1.5) + 
  theme_classic() + facet_wrap(~btm.plot, ncol=2, scales='free') + xlab('D2 Score \n(Norm. to Baseline)') + 
  ylab(expression(paste("Total ", italic("Mtb "), "CFU (log10)")))  + labs(color='Route') +
  scale_color_lancet() + theme(axis.text.x = element_text(angle = 45, hjust=1, size=12),
                               strip.text = element_text(size=9),
                               axis.text.y=element_text(size=12),
                               axis.title =element_text(size=12)) +
  stat_cor(size=3) + ylim(-0.2,7.7)

ggplot(all_out_route[btm %in% to_plot_btms & TimeafterBCG == 'd2'], aes(x=norm_score, y=log10(grans_nx+1))) +
  annotate(geom = 'rect', xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0.1, fill=alpha('grey',0.5)) +
  geom_jitter(aes(color=vax_group), width=0, height = 0.1, size=1.5) + 
  theme_classic() + facet_wrap(~btm.plot, ncol=2, scales='free') + xlab('D2 Score \n(Norm. to Baseline)') + 
  ylab('# of granulomas (log10)')  + labs(color='Route') +
  scale_color_lancet() + theme(axis.text.x = element_text(angle = 45, hjust=1, size=12),
                               strip.text = element_text(size=9),
                               axis.text.y=element_text(size=12),
                               axis.title =element_text(size=12)) +
  stat_cor(size=3) + ylim(-0.1,2.2)

ggplot(all_out_route[btm == "All M1 genes" & TimeafterBCG == 'd2' & vax_group == 'IV'], 
       aes(x=norm_score, y=log10(grans_nx+1))) +
  geom_jitter(aes(color=vax_group), width=0, height = 0.1, size=1.5) + 
  theme_classic() + facet_wrap(~btm, ncol=2, scales='free') + xlab('D2 Score \n(Norm. to Baseline)') + 
  ylab('# of granulomas (log10)')  + labs(color='Route') +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=12),
        strip.text = element_text(size=9),
        axis.text.y=element_text(size=12),
        axis.title =element_text(size=12)) +
  stat_cor(size=3)

