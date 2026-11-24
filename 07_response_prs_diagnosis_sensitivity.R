# PRS sensitivity: same matched participants with vs without 14 baseline diagnoses.
# Uses the saved dsm_14disorder.bl object also used in the original analysis.
library(emmeans)
data_dir <- Sys.getenv('ABCD_DATA_DIR',unset='../ABCD')
prepared_dir <- Sys.getenv('ABCD_PREPARED_DIR',unset='Data_and_Scripts')
output_dir <- Sys.getenv('RESPONSE_OUTPUT_DIR',unset='Revision_Analysis_20261122/results_diagnosis_20261125')
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
e <- new.env();load(file.path(prepared_dir,'abcd_SuicideData.RData'),envir=e)
gen <- read.table(file.path(data_dir,'abcd_prs','abcd_gene_covar.txt'),header=TRUE,stringsAsFactors=FALSE)
names(gen)[names(gen)=='IID'] <- 'subjectkey'
pcs <- paste0('genetic_pc_',1:10)
diags <- c('ADHD','ODD','CD','DEP','GAD','PAN','SEP','SOC','PTSD','PHO','HAL','DEL','BIP','OCD')
dx <- e$dsm_14disorder.bl[,c('subjectkey',diags)]
stopifnot(!anyDuplicated(dx$subjectkey),!anyDuplicated(gen$subjectkey),!anyDuplicated(e$prs_prsice$subjectkey),!anyDuplicated(e$abcd_envir.bl$subjectkey))
names(dx)[-1] <- paste0('diag_',diags)
diagnostic_terms <- names(dx)[-1]
for (term in diagnostic_terms) {
  x <- as.character(dx[[term]])
  if (any(!is.na(x) & !x %in% c('0','1'))) stop('Unexpected diagnostic coding: ',term)
  dx[[term]] <- as.numeric(x)
}
dat <- merge(e$abcd_envir.bl[,c('subjectkey','clust')],e$prs_prsice,by='subjectkey')
dat <- merge(dat,gen[,c('subjectkey','sex',pcs)],by='subjectkey')
dat <- merge(dat,dx,by='subjectkey')
dat$clust <- factor(dat$clust,levels=0:5);dat$sex <- factor(dat$sex)
dat <- dat[!is.na(dat$clust),]
traits <- c('adhd','bip','MDD','OCD','PTSD','scz'); out <- list(); checks <- list()
for(trait in traits){
 vars <- c(trait,'clust','sex',pcs,diagnostic_terms)
 sub <- dat[complete.cases(dat[,vars]),vars,drop=FALSE]
 sub$score <- as.numeric(scale(sub[[trait]]))
 stopifnot(all(table(sub$clust)>0),all(is.finite(sub$score)))
 for(mode in c('sex_pc','sex_pc_diagnoses')){
  adjust <- c('sex',pcs,if(mode=='sex_pc_diagnoses') diagnostic_terms)
  fit <- lm(as.formula(paste('score ~ clust +',paste(adjust,collapse=' + '))),data=sub)
  if(anyNA(coef(fit))) warning('Aliased coefficients: ',trait,' ',mode)
  tab <- as.data.frame(summary(contrast(emmeans(fit,~clust,cov.keep=character(0)),method='trt.vs.ctrl',ref=1,adjust='none'),infer=c(TRUE,TRUE),adjust='none'))
  tab$trait <- trait;tab$model <- mode;tab$n <- nrow(sub);out[[paste(trait,mode)]]<-tab
  checks[[paste(trait,mode)]] <- data.frame(trait=trait,model=mode,n=nrow(sub),rank=fit$rank,parameters=length(coef(fit)),aliased=sum(is.na(coef(fit))))
 }
}
out <- do.call(rbind,out);out$p_holm_30_within_model <- NA_real_
for(mode in unique(out$model)){i<-which(out$model==mode);out$p_holm_30_within_model[i]<-p.adjust(out$p.value[i],method='holm')}
write.csv(out,file.path(output_dir,'prs_sex_pc_diagnosis_comparison.csv'),row.names=FALSE)
write.csv(do.call(rbind,checks),file.path(output_dir,'prs_diagnosis_model_checks.csv'),row.names=FALSE)
writeLines(c('Each model pair uses the same complete cases; six PRSs standardized within their matched sample.',
 'Fourteen diagnostic variables from dsm_14disorder.bl; ascertainment and source mapping require verification.',
 'Holm p-value correction across 30 contrasts separately within each model; 95% intervals pointwise.',
 'No genotyping batch/array variable included; genetic QC and ancestry ascertainment remain to be verified.',
 'Diagnosis-free control selection can cause limited overlap in diagnostic covariates; interpret cautiously.'),file.path(output_dir,'analysis_notes.txt'))
