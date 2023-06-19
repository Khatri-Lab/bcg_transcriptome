###############
# test if module expression over time significantly differs by protection status in high & low dose recipients
# use general estimating equations
library(lme4)
library(gee)
library(MuMIn)
bcg1pheno <- fread('/labs/khatrilab/yiranliu/BCG/github/data/bcg_routecohort_metadata.csv')
bcg2pheno <- fread('/labs/khatrilab/yiranliu/BCG/github/data/ivbcg_dosecohort_metadata.csv', header = T)

setnames(bcg1pheno, c('description','animal_id','dose_group','protect_outcome'), c('SampleID','id','dose','protect_status'))
setnames(bcg2pheno, c('description','animal_id','dose_group','protect_outcome'), c('SampleID','id','dose','protect_status'))

bcg1pheno$TimeafterBCG <- factor(bcg1pheno$TimeafterBCG, levels = c('pre','d2','wk2','wk12'))
bcg2pheno$TimeafterBCG <- factor(bcg2pheno$TimeafterBCG, levels = c('pre','d2','wk2','wk4','wk12'))
bcg1pheno$dose <- factor(bcg1pheno$dose, levels = c('low','high'))
bcg2pheno$dose <- factor(bcg2pheno$dose, levels = c('low','high'))

modscores <- fread('/labs/khatrilab/yiranliu/BCG/github/wgcna/modScores.csv')

# Transform module scores into long format
modscores_norm <- data.table(modscores)
modscores_norm[,norm_score := score - score[TimeafterBCG == 'pre'], by=.(id,moduleID)]
modscores_cast <- dcast.data.table(modscores[TimeafterBCG %in% c('pre','d2','wk2','wk4')], id + dose + protect_status ~ moduleID + TimeafterBCG, value.var = 'score')
modscores_norm_cast <- dcast(modscores_norm[TimeafterBCG %in% c('d2','wk2','wk4')], id + dose + protect_status ~ moduleID + TimeafterBCG, value.var = 'norm_score')

modscores_cast <- merge.data.table(modscores_cast, bcg2pheno[TimeafterBCG == 'wk2',c('id','log_Mtb_CFU')],by='id')

long <- modscores[TimeafterBCG %in% c('pre','d2','wk2','wk4')]
long[TimeafterBCG == 'pre', time := 0]
long[TimeafterBCG == 'd2', time := 2]
long[TimeafterBCG == 'wk2', time := 14]
long[TimeafterBCG == 'wk4', time := 28]
long$id <- as.factor(long$id)
long$dose <- factor(long$dose, levels = c('low','high'))
long$study <- 'Dose'

# Do the same for route cohort
modscores_bcg1 <- fread('/labs/khatrilab/yiranliu/BCG/github/wgcna/modScores_routecohort.csv')

modscores_bcg1_cast <- dcast.data.table(modscores_bcg1[TimeafterBCG %in% c('pre','d2','wk2')], id + vax_group + protect_status ~ moduleID + TimeafterBCG, value.var = 'score')
modscores_bcg1_cast <- merge.data.table(modscores_bcg1_cast, bcg1pheno[TimeafterBCG == 'wk2',c('id','log_Mtb_CFU','study')],by='id')

long_bcg1 <- modscores_bcg1[TimeafterBCG %in% c('pre','d2','wk2')]
long_bcg1[TimeafterBCG == 'pre', time := 0]
long_bcg1[TimeafterBCG == 'd2', time := 2]
long_bcg1[TimeafterBCG == 'wk2', time := 14]
long_bcg1$id <- as.factor(long_bcg1$id)
long_bcg1$study <- 'Route'

## M1 ##
# high
dt.M1.high <- long[moduleID == 'M1' & dose == 'high']
cor(modscores_cast[dose == 'high',c('M1_pre','M1_d2','M1_wk2','M1_wk4'),with=F], use='pairwise.complete.obs')
gee.fit.M1.high <- gee(score~time + I(time^2) + I(time^3),
                  data = dt.M1.high,
                  id=id,
                  corstr="unstructured")
summary(gee.fit.M1.high)
QIC(gee.fit.M1.high)
M1.high.res <- as.data.table(coef(summary(gee.fit.M1.high)), keep.rownames = 'variable')
M1.high.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M1.high.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M1.high.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M1.high.res)[-1]) set(M1.high.res, j = j, value = round(M1.high.res[[j]],4))
M1.high.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M1.high.res$moduleID <- 'M1'
M1.high.res$dose <- 'high'
M1.high.res$cohort <- 'dose'

# low
dt.M1.low <- long[moduleID == 'M1' & dose == 'low']
cor(modscores_cast[dose == 'low',c('M1_pre','M1_d2','M1_wk2','M1_wk4'),with=F], use='pairwise.complete.obs')
gee.fit.M1.low <- gee(score~time + I(time^2) + I(time^3),
                       data = dt.M1.low,
                       id=id,
                       corstr="unstructured")
summary(gee.fit.M1.low)
QIC(gee.fit.M1.low)
M1.low.res <- as.data.table(coef(summary(gee.fit.M1.low)), keep.rownames = 'variable')
M1.low.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M1.low.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M1.low.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M1.low.res)[-1]) set(M1.low.res, j = j, value = round(M1.low.res[[j]],4))
M1.low.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M1.low.res$moduleID <- 'M1'
M1.low.res$dose <- 'low'
M1.low.res$cohort <- 'dose'

# IV high across dose and route cohorts
dt.M1.bcg1 <- rbind(long_bcg1[moduleID == 'M1' & vax_group == 'IV'], long[moduleID == 'M1' & dose == 'high'], fill=T)
# cor(modscores_bcg1_cast[vax_group == 'IV',c('M1_pre','M1_d2','M1_wk2'),with=F], use='pairwise.complete.obs')
gee.fit.M1.bcg1 <- gee(score~time + I(time^2) + I(time^3) + study,
                       data = dt.M1.bcg1,
                       id=id,
                       corstr="unstructured")
summary(gee.fit.M1.bcg1)
QIC(gee.fit.M1.bcg1)
M1.bcg1.res <- as.data.table(coef(summary(gee.fit.M1.bcg1)), keep.rownames = 'variable')
M1.bcg1.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M1.bcg1.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M1.bcg1.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M1.bcg1.res)[-1]) set(M1.bcg1.res, j = j, value = round(M1.bcg1.res[[j]],4))
M1.bcg1.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M1.bcg1.res$moduleID <- 'M1'
M1.bcg1.res$dose <- 'high'
M1.bcg1.res$cohort <- 'dose & route'

## M2 ##
# high
dt.M2.high <- long[moduleID == 'M2' & dose == 'high']
cor(modscores_cast[dose == 'high',c('M2_pre','M2_d2','M2_wk2','M2_wk4'),with=F], use='pairwise.complete.obs')
gee.fit.M2.high <- gee(score~time + I(time^2) + I(time^3),
                       data = dt.M2.high,
                       id=id,
                       corstr="exchangeable")
summary(gee.fit.M2.high)
QIC(gee.fit.M2.high)
M2.high.res <- as.data.table(coef(summary(gee.fit.M2.high)), keep.rownames = 'variable')
M2.high.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M2.high.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M2.high.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M2.high.res)[-1]) set(M2.high.res, j = j, value = round(M2.high.res[[j]],4))
M2.high.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M2.high.res$moduleID <- 'M2'
M2.high.res$dose <- 'high'
M2.high.res$cohort <- 'dose'

# low
dt.M2.low <- long[moduleID == 'M2' & dose == 'low']
cor(modscores_cast[dose == 'low',c('M2_pre','M2_d2','M2_wk2','M2_wk4'),with=F], use='pairwise.complete.obs')
gee.fit.M2.low <- gee(score~time + I(time^2) + I(time^3),
                      data = dt.M2.low,
                      id=id,
                      corstr="exchangeable")
summary(gee.fit.M2.low)
QIC(gee.fit.M2.low)
M2.low.res <- as.data.table(coef(summary(gee.fit.M2.low)), keep.rownames = 'variable')
M2.low.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M2.low.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M2.low.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M2.low.res)[-1]) set(M2.low.res, j = j, value = round(M2.low.res[[j]],4))
M2.low.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M2.low.res$moduleID <- 'M2'
M2.low.res$dose <- 'low'
M2.low.res$cohort <- 'dose'

# IV high across dose and route cohorts
dt.M2.bcg1 <- rbind(long_bcg1[moduleID == 'M2' & vax_group == 'IV'], long[moduleID == 'M2' & dose == 'high'], fill=T)
# cor(modscores_bcg1_cast[vax_group == 'IV',c('M2_pre','M2_d2','M2_wk2'),with=F], use='pairwise.complete.obs')
gee.fit.M2.bcg1 <- gee(score~time + I(time^2) + I(time^3) + study,
                       data = dt.M2.bcg1,
                       id=id,
                       corstr="unstructured")
summary(gee.fit.M2.bcg1)
QIC(gee.fit.M2.bcg1)
M2.bcg1.res <- as.data.table(coef(summary(gee.fit.M2.bcg1)), keep.rownames = 'variable')
M2.bcg1.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M2.bcg1.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M2.bcg1.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M2.bcg1.res)[-1]) set(M2.bcg1.res, j = j, value = round(M2.bcg1.res[[j]],4))
M2.bcg1.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M2.bcg1.res$moduleID <- 'M2'
M2.bcg1.res$dose <- 'high'
M2.bcg1.res$cohort <- 'dose & route'


## M3 ##
# high
dt.M3.high <- long[moduleID == 'M3' & dose == 'high']
cor(modscores_cast[dose == 'high',c('M3_pre','M3_d2','M3_wk2','M3_wk4'),with=F], use='pairwise.complete.obs')
gee.fit.M3.high <- gee(score~time + I(time^2) + I(time^3),
                       data = dt.M3.high,
                       id=id,
                       corstr="unstructured")
summary(gee.fit.M3.high)
QIC(gee.fit.M3.high)
M3.high.res <- as.data.table(coef(summary(gee.fit.M3.high)), keep.rownames = 'variable')
M3.high.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M3.high.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M3.high.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M3.high.res)[-1]) set(M3.high.res, j = j, value = round(M3.high.res[[j]],4))
M3.high.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M3.high.res$moduleID <- 'M3'
M3.high.res$dose <- 'high'
M3.high.res$cohort <- 'dose'

# low
dt.M3.low <- long[moduleID == 'M3' & dose == 'low']
cor(modscores_cast[dose == 'low',c('M3_pre','M3_d2','M3_wk2','M3_wk4'),with=F], use='pairwise.complete.obs')
gee.fit.M3.low <- gee(score~time + I(time^2) + I(time^3),
                      data = dt.M3.low,
                      id=id,
                      corstr="exchangeable")
summary(gee.fit.M3.low)
QIC(gee.fit.M3.low)
M3.low.res <- as.data.table(coef(summary(gee.fit.M3.low)), keep.rownames = 'variable')
M3.low.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M3.low.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M3.low.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M3.low.res)[-1]) set(M3.low.res, j = j, value = round(M3.low.res[[j]],4))
M3.low.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M3.low.res$moduleID <- 'M3'
M3.low.res$dose <- 'low'
M3.low.res$cohort <- 'dose'

# IV high across dose and route cohorts
dt.M3.bcg1 <- rbind(long_bcg1[moduleID == 'M3' & vax_group == 'IV'], long[moduleID == 'M3' & dose == 'high'], fill=T)
# cor(modscores_bcg1_cast[vax_group == 'IV',c('M3_pre','M3_d2','M3_wk2'),with=F], use='pairwise.complete.obs')
gee.fit.M3.bcg1 <- gee(score~time + I(time^2) + I(time^3) + study,
                       data = dt.M3.bcg1,
                       id=id,
                       corstr="unstructured")
summary(gee.fit.M3.bcg1)
QIC(gee.fit.M3.bcg1)
M3.bcg1.res <- as.data.table(coef(summary(gee.fit.M3.bcg1)), keep.rownames = 'variable')
M3.bcg1.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M3.bcg1.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M3.bcg1.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M3.bcg1.res)[-1]) set(M3.bcg1.res, j = j, value = round(M3.bcg1.res[[j]],4))
M3.bcg1.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M3.bcg1.res$moduleID <- 'M3'
M3.bcg1.res$dose <- 'high'
M3.bcg1.res$cohort <- 'dose & route'

## M4 ##
# high
dt.M4.high <- long[moduleID == 'M4' & dose == 'high']
cor(modscores_cast[dose == 'high',c('M4_pre','M4_d2','M4_wk2','M4_wk4'),with=F], use='pairwise.complete.obs')
gee.fit.M4.high <- gee(score~time + I(time^2) + I(time^3),
                       data = dt.M4.high,
                       id=id,
                       corstr="exchangeable")
summary(gee.fit.M4.high)
QIC(gee.fit.M4.high)
M4.high.res <- as.data.table(coef(summary(gee.fit.M4.high)), keep.rownames = 'variable')
M4.high.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M4.high.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M4.high.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M4.high.res)[-1]) set(M4.high.res, j = j, value = round(M4.high.res[[j]],4))
M4.high.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M4.high.res$moduleID <- 'M4'
M4.high.res$dose <- 'high'
M4.high.res$cohort <- 'dose'

# low
dt.M4.low <- long[moduleID == 'M4' & dose == 'low']
cor(modscores_cast[dose == 'low',c('M4_pre','M4_d2','M4_wk2','M4_wk4'),with=F], use='pairwise.complete.obs')
gee.fit.M4.low <- gee(score~time + I(time^2) + I(time^3),
                      data = dt.M4.low,
                      id=id,
                      corstr="exchangeable")
summary(gee.fit.M4.low)
QIC(gee.fit.M4.low)
M4.low.res <- as.data.table(coef(summary(gee.fit.M4.low)), keep.rownames = 'variable')
M4.low.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M4.low.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M4.low.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M4.low.res)[-1]) set(M4.low.res, j = j, value = round(M4.low.res[[j]],4))
M4.low.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M4.low.res$moduleID <- 'M4'
M4.low.res$dose <- 'low'
M4.low.res$cohort <- 'dose'

# IV high across dose and route cohorts
dt.M4.bcg1 <- rbind(long_bcg1[moduleID == 'M4' & vax_group == 'IV'], long[moduleID == 'M4' & dose == 'high'], fill=T)
# cor(modscores_bcg1_cast[vax_group == 'IV',c('M4_pre','M4_d2','M4_wk2'),with=F], use='pairwise.complete.obs')
gee.fit.M4.bcg1 <- gee(score~time + I(time^2) + I(time^3) + study,
                       data = dt.M4.bcg1,
                       id=id,
                       corstr="unstructured")
summary(gee.fit.M4.bcg1)
QIC(gee.fit.M4.bcg1)
M4.bcg1.res <- as.data.table(coef(summary(gee.fit.M4.bcg1)), keep.rownames = 'variable')
M4.bcg1.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M4.bcg1.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M4.bcg1.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M4.bcg1.res)[-1]) set(M4.bcg1.res, j = j, value = round(M4.bcg1.res[[j]],4))
M4.bcg1.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M4.bcg1.res$moduleID <- 'M4'
M4.bcg1.res$dose <- 'high'
M4.bcg1.res$cohort <- 'dose & route'

## M5 ##
# high
dt.M5.high <- long[moduleID == 'M5' & dose == 'high']
cor(modscores_cast[dose == 'high',c('M5_pre','M5_d2','M5_wk2','M5_wk4'),with=F], use='pairwise.complete.obs')
gee.fit.M5.high <- gee(score~time + I(time^2) + I(time^3),
                       data = dt.M5.high,
                       id=id,
                       corstr="exchangeable")
summary(gee.fit.M5.high)
QIC(gee.fit.M5.high)
M5.high.res <- as.data.table(coef(summary(gee.fit.M5.high)), keep.rownames = 'variable')
M5.high.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M5.high.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M5.high.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M5.high.res)[-1]) set(M5.high.res, j = j, value = round(M5.high.res[[j]],4))
M5.high.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M5.high.res$moduleID <- 'M5'
M5.high.res$dose <- 'high'
M5.high.res$cohort <- 'dose'

# low
dt.M5.low <- long[moduleID == 'M5' & dose == 'low']
cor(modscores_cast[dose == 'low',c('M5_pre','M5_d2','M5_wk2','M5_wk4'),with=F], use='pairwise.complete.obs')
gee.fit.M5.low <- gee(score~time + I(time^2) + I(time^3),
                      data = dt.M5.low,
                      id=id,
                      corstr="exchangeable")
summary(gee.fit.M5.low)
QIC(gee.fit.M5.low)
M5.low.res <- as.data.table(coef(summary(gee.fit.M5.low)), keep.rownames = 'variable')
M5.low.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M5.low.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M5.low.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M5.low.res)[-1]) set(M5.low.res, j = j, value = round(M5.low.res[[j]],4))
M5.low.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M5.low.res$moduleID <- 'M5'
M5.low.res$dose <- 'low'
M5.low.res$cohort <- 'dose'


# IV high across dose and route cohorts
dt.M5.bcg1 <- rbind(long_bcg1[moduleID == 'M5' & vax_group == 'IV'], long[moduleID == 'M5' & dose == 'high'], fill=T)
# cor(modscores_bcg1_cast[vax_group == 'IV',c('M5_pre','M5_d2','M5_wk2'),with=F], use='pairwise.complete.obs')
gee.fit.M5.bcg1 <- gee(score~time + I(time^2) + I(time^3) + study,
                       data = dt.M5.bcg1,
                       id=id,
                       corstr="exchangeable")
summary(gee.fit.M5.bcg1)
QIC(gee.fit.M5.bcg1)
M5.bcg1.res <- as.data.table(coef(summary(gee.fit.M5.bcg1)), keep.rownames = 'variable')
M5.bcg1.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M5.bcg1.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M5.bcg1.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M5.bcg1.res)[-1]) set(M5.bcg1.res, j = j, value = round(M5.bcg1.res[[j]],4))
M5.bcg1.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M5.bcg1.res$moduleID <- 'M5'
M5.bcg1.res$dose <- 'high'
M5.bcg1.res$cohort <- 'dose & route'


## M6 ##
# high
dt.M6.high <- long[moduleID == 'M6' & dose == 'high']
cor(modscores_cast[dose == 'high',c('M6_pre','M6_d2','M6_wk2','M6_wk4'),with=F], use='pairwise.complete.obs')
gee.fit.M6.high <- gee(score~time + I(time^2) + I(time^3),
                       data = dt.M6.high,
                       id=id,
                       corstr="exchangeable")
summary(gee.fit.M6.high)
QIC(gee.fit.M6.high)
M6.high.res <- as.data.table(coef(summary(gee.fit.M6.high)), keep.rownames = 'variable')
M6.high.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M6.high.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M6.high.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M6.high.res)[-1]) set(M6.high.res, j = j, value = round(M6.high.res[[j]],4))
M6.high.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M6.high.res$moduleID <- 'M6'
M6.high.res$dose <- 'high'
M6.high.res$cohort <- 'dose'

# low
dt.M6.low <- long[moduleID == 'M6' & dose == 'low']
cor(modscores_cast[dose == 'low',c('M6_pre','M6_d2','M6_wk2','M6_wk4'),with=F], use='pairwise.complete.obs')
gee.fit.M6.low <- gee(score~time + I(time^2) + I(time^3),
                      data = dt.M6.low,
                      id=id,
                      corstr="exchangeable")
summary(gee.fit.M6.low)
QIC(gee.fit.M6.low)
M6.low.res <- as.data.table(coef(summary(gee.fit.M6.low)), keep.rownames = 'variable')
M6.low.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M6.low.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M6.low.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M6.low.res)[-1]) set(M6.low.res, j = j, value = round(M6.low.res[[j]],4))
M6.low.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M6.low.res$moduleID <- 'M6'
M6.low.res$dose <- 'low'
M6.low.res$cohort <- 'dose'

# IV high across dose and route cohorts
dt.M6.bcg1 <- rbind(long_bcg1[moduleID == 'M6' & vax_group == 'IV'], long[moduleID == 'M6' & dose == 'high'], fill=T)
# cor(modscores_bcg1_cast[vax_group == 'IV',c('M6_pre','M6_d2','M6_wk2'),with=F], use='pairwise.complete.obs')
gee.fit.M6.bcg1 <- gee(score~time + I(time^2) + I(time^3) + study,
                       data = dt.M6.bcg1,
                       id=id,
                       corstr="unstructured")
summary(gee.fit.M6.bcg1)
QIC(gee.fit.M6.bcg1)
M6.bcg1.res <- as.data.table(coef(summary(gee.fit.M6.bcg1)), keep.rownames = 'variable')
M6.bcg1.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M6.bcg1.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M6.bcg1.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M6.bcg1.res)[-1]) set(M6.bcg1.res, j = j, value = round(M6.bcg1.res[[j]],4))
M6.bcg1.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M6.bcg1.res$moduleID <- 'M6'
M6.bcg1.res$dose <- 'high'
M6.bcg1.res$cohort <- 'dose & route'


## M7 ##
# high
dt.M7.high <- long[moduleID == 'M7' & dose == 'high']
cor(modscores_cast[dose == 'high',c('M7_pre','M7_d2','M7_wk2','M7_wk4'),with=F], use='pairwise.complete.obs')
gee.fit.M7.high <- gee(score~time + I(time^2) + I(time^3),
                       data = dt.M7.high,
                       id=id,
                       corstr="exchangeable")
summary(gee.fit.M7.high)
QIC(gee.fit.M7.high)
M7.high.res <- as.data.table(coef(summary(gee.fit.M7.high)), keep.rownames = 'variable')
M7.high.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M7.high.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M7.high.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M7.high.res)[-1]) set(M7.high.res, j = j, value = round(M7.high.res[[j]],4))
M7.high.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M7.high.res$moduleID <- 'M7'
M7.high.res$dose <- 'high'
M7.high.res$cohort <- 'dose'

# low
dt.M7.low <- long[moduleID == 'M7' & dose == 'low']
cor(modscores_cast[dose == 'low',c('M7_pre','M7_d2','M7_wk2','M7_wk4'),with=F], use='pairwise.complete.obs')
gee.fit.M7.low <- gee(score~time + I(time^2) + I(time^3),
                      data = dt.M7.low,
                      id=id,
                      corstr="exchangeable")
summary(gee.fit.M7.low)
QIC(gee.fit.M7.low)
M7.low.res <- as.data.table(coef(summary(gee.fit.M7.low)), keep.rownames = 'variable')
M7.low.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M7.low.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M7.low.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M7.low.res)[-1]) set(M7.low.res, j = j, value = round(M7.low.res[[j]],4))
M7.low.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M7.low.res$moduleID <- 'M7'
M7.low.res$dose <- 'low'
M7.low.res$cohort <- 'dose'

# IV high across dose and route cohorts
dt.M7.bcg1 <- rbind(long_bcg1[moduleID == 'M7' & vax_group == 'IV'], long[moduleID == 'M7' & dose == 'high'], fill=T)
# cor(modscores_bcg1_cast[vax_group == 'IV',c('M7_pre','M7_d2','M7_wk2'),with=F], use='pairwise.complete.obs')
gee.fit.M7.bcg1 <- gee(score~time + I(time^2) + I(time^3) + study,
                       data = dt.M7.bcg1,
                       id=id,
                       corstr="exchangeable")
summary(gee.fit.M7.bcg1)
QIC(gee.fit.M7.bcg1)
M7.bcg1.res <- as.data.table(coef(summary(gee.fit.M7.bcg1)), keep.rownames = 'variable')
M7.bcg1.res[,pvalue := 2*pnorm(abs(`Robust z`), lower.tail = F)]
M7.bcg1.res[,est.lo := Estimate - 1.96*`Robust S.E.`]
M7.bcg1.res[,est.up := Estimate + 1.96*`Robust S.E.`]
for (j in colnames(M7.bcg1.res)[-1]) set(M7.bcg1.res, j = j, value = round(M7.bcg1.res[[j]],4))
M7.bcg1.res[,est_wCI := paste0(Estimate, ' (', est.lo, ', ', est.up, ')')]
M7.bcg1.res$moduleID <- 'M7'
M7.bcg1.res$dose <- 'high'
M7.bcg1.res$cohort <- 'dose & route'

all.res <- rbindlist(list(M1.high.res, M1.low.res, M1.bcg1.res,
                          M2.high.res, M2.low.res, M2.bcg1.res,
                          M3.high.res, M3.low.res, M3.bcg1.res,
                          M4.high.res, M4.low.res, M4.bcg1.res,
                          M5.high.res, M5.low.res, M5.bcg1.res,
                          M6.high.res, M6.low.res, M6.bcg1.res,
                          M7.high.res, M7.low.res, M7.bcg1.res))
