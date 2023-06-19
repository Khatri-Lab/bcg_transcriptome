##############################################################################################################
######################################### FUNCTIONS FOR DEG ANALYSIS #########################################
##############################################################################################################

# function to remove genes from txi object
remove_from_txi <- function(txi, genes_to_remove=c(), samps_to_remove=c()){
  if (is.null(genes_to_remove) & !is.null(samps_to_remove)){
    idx_remove_samps <- which(colnames(txi$counts) %in% samps_to_remove)
    if ('abundance' %in% names(txi)){
      txi$abundance <- txi$abundance[,-idx_remove_samps]
    }
    txi$counts <- txi$counts[,-idx_remove_samps]
    if ('infReps' %in% names(txi)){
      txi$infReps <- lapply(txi$infReps, function(dt) dt[,-idx_remove_samps])
    }
    txi$length <- txi$length[,-idx_remove_samps]
  } else if (!is.null(genes_to_remove) & is.null(samps_to_remove)){
    idx_remove_genes <- which(rownames(txi$counts) %in% genes_to_remove)
    if ('abundance' %in% names(txi)){
      txi$abundance <- txi$abundance[-idx_remove_genes,]
    }
    txi$counts <- txi$counts[-idx_remove_genes,]
    if ('infReps' %in% names(txi)){
      txi$infReps <- lapply(txi$infReps, function(dt) dt[-idx_remove_genes,])
    }
    txi$length <- txi$length[-idx_remove_genes,]
  } else if (!is.null(genes_to_remove) & !is.null(samps_to_remove)){
    idx_remove_genes <- which(rownames(txi$counts) %in% genes_to_remove)
    idx_remove_samps <- which(colnames(txi$counts) %in% samps_to_remove)
    if ('abundance' %in% names(txi)){
      txi$abundance <- txi$abundance[-idx_remove_genes,-idx_remove_samps]
    }
    txi$counts <- txi$counts[-idx_remove_genes,-idx_remove_samps]
    if ('infReps' %in% names(txi)){
      txi$infReps <- lapply(txi$infReps, function(dt) dt[-idx_remove_genes,-idx_remove_samps])
    }
    txi$length <- txi$length[-idx_remove_genes,-idx_remove_samps]
  }
  return(txi)
}

############################################## DESeq2 Functions ##############################################

library(DESeq2)
library(apeglm)

getdds <- function(ids, txi, pheno, design){
  txi_new <- remove_from_txi(txi, samps_to_remove = setdiff(pheno$SampleID, ids))
  pheno <- data.frame(pheno[SampleID %in% ids,], row.names = 'SampleID')
  dds <- DESeqDataSetFromTximport(txi_new,
                                  colData = pheno[colnames(txi_new$counts),],
                                  design = design)
  n <- ncol(counts(dds))
  keep <- rowSums(counts(dds) >= 5) >= round(n / 3) # keep only genes that have >= 5 reads in >= 1/3 of the samples
  dds <- dds[keep,]
  dds <- DESeq(dds)
  return(dds)
}


# Get results object from dds object, feed in resultName to use, optional lfcshrink
getRes <- function(dds, resultName='', resultContrast=c(), padj_cutoff = 0.1, lfc_cutoff=0, lfcShrink=T, lfcShrinkType='apeglm'){
  if (resultName != ''){
    res <- results(dds, name=resultName, alpha = padj_cutoff)
    if (lfcShrink){
      res <- lfcShrink(dds, coef=resultName, res = res, type = lfcShrinkType)
    }
  } else if (!is.null(resultContrast)){
    res <- results(dds, contrast=resultContrast, alpha = padj_cutoff, lfcThreshold = lfc_cutoff)
    if (lfcShrink){
      if (lfcShrinkType == 'apeglm'){
        print("'apeglm' shrinkage only for use with 'coef', not contrast. using normal prior instead")
      }
      res <- lfcShrink(dds, contrast=resultContrast, res = res, type = 'normal')
    }
  }
  
  return(res)
}


# get sig genes from a res object
getsiggenes <- function(res_obj, padj_cutoff=0.1, lfc_cutoff=log(1.5, base=2), pos_only = F, neg_only = F, ids = c(), ctrlgrep = 'pre', trtgrep = ''){
  res <- res_obj %>% as.data.frame %>% as.data.table(keep.rownames = 'GeneID')
  res %<>% subset(baseMean > 0 & !is.na(padj))
  if (pos_only){
    sigGenes = unique(res[log2FoldChange >= lfc_cutoff & padj <= padj_cutoff]$GeneID)
    return(sigGenes)
  } else if (neg_only){
    sigGenes = unique(res[log2FoldChange <= -lfc_cutoff & padj <= padj_cutoff]$GeneID)
    return(sigGenes)
  } else{
    sigGenes = unique(res[abs(log2FoldChange) >= lfc_cutoff & padj <= padj_cutoff]$GeneID)
    return(sigGenes)
  }
}


############################################## BTM Functions ##############################################
source('/labs/khatrilab/yiranliu/functions/setEnrichment.R')

modGenes_mmul <- readRDS('/labs/khatrilab/yiranliu/BCG/mmul/modGenes_mmul.RDS')

get_BTM_dt <- function(res_obj, sig_genes = c(), lfc_cutoff=log(1.5, base = 2), padj_cutoff=0.1, genes_convert=genes_convert){
  res_dt <- as.data.table(data.frame(res_obj), keep.rownames = 'Gene.stable.ID')
  gene_universe = unique(res_dt$Gene.stable.ID)
  if (!is.null(sig_genes)){
    gene_sig_dt = res_dt[Gene.stable.ID %in% sig_genes]
  } else {
    gene_sig_dt = res_dt[abs(log2FoldChange) >= lfc_cutoff & padj <= padj_cutoff,]
  }
  gene_sig_es_vector = gene_sig_dt$log2FoldChange %>% set_names(gene_sig_dt$Gene.stable.ID)
  
  BTMs = setEnrichment(
    gSign = gene_sig_es_vector, # use this if you have an experiment with effect sizes and you want to plot them over the modules
    # gSign = gene_sig, # use this if you do not have effect sizes
    geneSets = modGenes_mmul, 
    sigUniverse = gene_universe, 
    setIDs = names(modGenes_mmul)
  )
  BTMs
}

# merge multiple BTM lists into one dt
merge_BTMs <- function(BTMlist, suffx, merge_by=c('set.name',
                                                  'original.size',
                                                  'filtered.size'), min_relevant_genes=3, subset_cols=c('p.adj',
                                                                                                   'relevant.genes'), padj_cutoff = 0.01){
  if (length(BTMlist) != length(suffx)){
    print('Error: BTM list is not the same length as suffix list')
    break
  }
  
  # filter out BTMs with fewer relevant genes than threshold and subset to relevant columns
  BTM_dtlist <- lapply(BTMlist, function(i) i@enrichmentTable[relevant.genes >= min_relevant_genes, c(merge_by, subset_cols), with=F])
  
  # add suffixes to column names before merging (add suffix to all columns except the one(s) to merge by)
  for (i in 1:length(BTM_dtlist)){
    colnames(BTM_dtlist[[i]])[-which(colnames(BTM_dtlist[[i]]) %in% merge_by)] <- paste(setdiff(colnames(BTM_dtlist[[i]]), merge_by), suffx[i], sep='_')
  }
  BTM_merged <- Reduce(
    function(x, y, ...) merge.data.table(x, y, all = TRUE, by=merge_by),
    BTM_dtlist
  )
  colnames(BTM_merged) <- gsub('p.adj_', '', colnames(BTM_merged))
  
  BTM_merged_melt <- melt.data.table(BTM_merged[,!grepl('relevant.genes',colnames(BTM_merged)),with=F], id.vars = merge_by, value.name = 'padj')
  BTM_merged_melt[,p_inv := log10(1/padj)]
  
  BTM_merged_relevantgenes_melt <- melt(BTM_merged[,c(merge_by, grep('relevant.genes',colnames(BTM_merged), value = T)),
                                                   with=F], id.vars = merge_by, value.name = 'relevant.genes')
  BTM_merged_relevantgenes_melt$variable <- gsub('relevant.genes_', '', BTM_merged_relevantgenes_melt$variable)
  BTM_final_dt <- na.omit(merge(BTM_merged_melt, BTM_merged_relevantgenes_melt, 
                           by=c(merge_by, 'variable')))
  
  BTM_final_dt[,padj_sig := ifelse(padj <= padj_cutoff, 'y','n')]
  BTM_final_dt$set.name <- factor(BTM_final_dt$set.name)
  
  return(BTM_final_dt)
}

spot.theme_2 <- list(
  theme_classic(),
  theme(axis.ticks.y=element_blank(), axis.text=element_text(size = 10),
        axis.text.x=element_text(angle=45, hjust = 1),
        axis.line=element_blank(),
        text = element_text(size = 10),
        plot.margin = unit(c(5,5,5,5), "mm")),
  scale_x_discrete(position = "bottom"))


############################################## Other Functions ##############################################

# get gene name from gene ID for plotting
getgenename_single <- function(geneID){
  if (geneID %in% mmulgenes_only$`Gene stable ID`){
    name <- mmulgenes_only[`Gene stable ID` %in% geneID, 'Gene name']
    if (name != ''){
      geneID <- name
    }
  }
  return(geneID)
}

getgenename <- function(geneIDs){
  geneIDs <- sapply(geneIDs, getgenename_single)
  return(geneIDs)
}

geomMean <- function (x) {
  if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    stop("argument is not numeric or logical: returning NA")
  }
  if (any(x < 0)) {
    stop("'x' contains negative value(s)")
  }
  return(exp(sum(log(x))/length(x)))
}

# function to return melted gene expression matrix, merged with pheno table. add faux sample if needed for later imputation
get_expr_pheno <- function(geneIDs, genematx, pheno_dt, sampleid_colname = 'SampleID', add=NA){
  if ('ENSEMBLE_Gene_IDS' %in% colnames(genematx)){
    expr <- genematx[ENSEMBLE_Gene_IDS %in% geneIDs]
  } else {
    included_genes <- geneIDs[which(geneIDs %in% rownames(genematx))]
    expr <- data.table(genematx[included_genes,], keep.rownames = 'ENSEMBLE_Gene_IDS')
  }
  if (!is.na(add)){
    expr[[add]] <- NA
  }
  expr <- melt(expr, id.vars = 'ENSEMBLE_Gene_IDS')
  colnames(expr) <- c('ENSEMBLE_Gene_IDS','SampleID','value')
  if (sampleid_colname != 'SampleID'){
    pheno_dt <- data.table(pheno_dt)
    setnames(pheno_dt, sampleid_colname, 'SampleID')
  }
  pheno_dt$SampleID <- as.character(pheno_dt$SampleID)
  expr$SampleID <- as.character(expr$SampleID)
  expr_pheno <- data.table(merge(expr, pheno_dt, by = 'SampleID'))
  return(expr_pheno)
}

get_score_pheno <- function(up_geneIDs, down_geneIDs, genematx, pheno_dt){
  if(length(up_geneIDs) == 0){
    scores <- 0 - apply(genematx[down_geneIDs, , drop=F], 2, geomMean)
  } else if (length(down_geneIDs) == 0){
    scores <- apply(genematx[up_geneIDs, , drop=F], 2, geomMean)
  } else {
    scores <- apply(genematx[up_geneIDs, , drop=F], 2, geomMean) - apply(genematx[down_geneIDs, , drop=F], 2, geomMean)
  }
  merged <- data.table(data.frame(scores), keep.rownames = T)
  colnames(merged) <- c('SampleID','score')
  merged <- merge(merged, pheno_dt, by = 'SampleID')
  return(merged)
}

# make counts/VST matrix normalized to baseline
norm_to_baseline <- function(vstmat, pheno_dt=comb_pheno_prechal, impute_baseline=F){
  if (!'ENSEMBLE_Gene_IDS' %in% colnames(vstmat)){
    vstmat <- data.table(vstmat, keep.rownames = 'ENSEMBLE_Gene_IDS')
  }
  vst_melt <- melt.data.table(vstmat, id.vars = 'ENSEMBLE_Gene_IDS', variable.name = 'SampleID')
  vst_melt <- merge.data.table(vst_melt, pheno_dt[,.(SampleID, id, TimeafterBCG)], by='SampleID')
  
  if (impute_baseline){ # if T, add pseudo baseline sample for animals missing baseline sample using median norm expression of all samples
    nobaseline <- unique(pheno_dt$id)[which(!unique(pheno_dt$id) %in% pheno_dt[TimeafterBCG == 'pre', id])]
    temp <- unique(vst_melt[id %in% nobaseline], by=c('id','ENSEMBLE_Gene_IDS'))
    temp$TimeafterBCG <- 'pre'
    temp$SampleID <- paste(temp$id, 'pre_imput', sep='_')
    temp$value <- NA
    
    vst_melt <- rbind(vst_melt, temp)
    vst_melt[, med := median(value[TimeafterBCG == 'pre'], na.rm=T), by=.(ENSEMBLE_Gene_IDS)]
    vst_melt[id %in% nobaseline & TimeafterBCG == 'pre', value := med]

  }
  vst_melt[, norm_expr := value - value[TimeafterBCG == 'pre'], by=.(id,ENSEMBLE_Gene_IDS)]
  vst_melt <- vst_melt[,.(SampleID,ENSEMBLE_Gene_IDS,norm_expr)]
  vst_cast <- dcast(vst_melt, ENSEMBLE_Gene_IDS ~ SampleID)
  vst_cast <- data.frame(vst_cast[,-1], row.names = vst_cast$ENSEMBLE_Gene_IDS, check.names = F)
  vst_cast
}

