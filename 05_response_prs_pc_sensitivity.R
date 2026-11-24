# R1.7 / R2.7: exploratory PRS sensitivity analysis with sex and ancestry PCs.
# Run from the project root. Genotyping batch/array and QC provenance are not
# established by these files; do not describe this as a complete QC analysis.
library(emmeans)
data_dir <- Sys.getenv('ABCD_DATA_DIR', unset = '../ABCD')
prepared_dir <- Sys.getenv('ABCD_PREPARED_DIR', unset = 'Data_and_Scripts')
output_dir <- Sys.getenv('RESPONSE_OUTPUT_DIR', unset = 'Revision_Analysis_20261122/results_updated_20261124')
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
load(file.path(prepared_dir, 'abcd_SuicideData.RData'))
genetic_file <- file.path(data_dir, 'abcd_prs', 'abcd_gene_covar.txt')
if (!file.exists(genetic_file)) stop('Genetic covariates missing: ', genetic_file)
gen <- read.table(genetic_file, header=TRUE, stringsAsFactors=FALSE)
names(gen)[names(gen)=='IID'] <- 'subjectkey'
pcs <- paste0('genetic_pc_',1:10)
stopifnot(all(c('subjectkey','sex',pcs) %in% names(gen)))
stopifnot(!anyDuplicated(gen$subjectkey), !anyDuplicated(prs_prsice$subjectkey), !anyDuplicated(abcd_envir.bl$subjectkey))
dat <- merge(abcd_envir.bl[,c('subjectkey','clust')], prs_prsice, by='subjectkey')
dat <- merge(dat, gen[,c('subjectkey','sex',pcs)], by='subjectkey')
dat$clust <- factor(dat$clust, levels=0:5)
dat$sex <- factor(dat$sex)
dat <- dat[!is.na(dat$clust),]
traits <- c('adhd','bip','MDD','OCD','PTSD','scz')
stopifnot(all(traits %in% names(dat)))
results <- list(); details <- list()
for (trait in traits) {
  needed <- c(trait,'clust','sex',pcs)
  sub <- dat[complete.cases(dat[,needed]),needed,drop=FALSE]
  sub$score <- as.numeric(scale(sub[[trait]]))
  if (any(table(sub$clust)==0) || !all(is.finite(sub$score))) stop('Invalid groups or score: ',trait)
  fit <- lm(as.formula(paste('score ~ clust + sex +',paste(pcs,collapse=' + '))),data=sub)
  contrasts <- contrast(emmeans(fit,~clust),method='trt.vs.ctrl',ref=1,adjust='none')
  tab <- as.data.frame(summary(contrasts,infer=c(TRUE,TRUE),adjust='none'))
  tab$trait <- trait; tab$n <- nrow(sub)
  results[[trait]] <- tab
  details[[trait]] <- data.frame(trait=trait,n=nrow(sub),groups=paste(names(table(sub$clust)),table(sub$clust),collapse='; '),
    formula=paste(deparse(formula(fit)),collapse=' '),r_squared=summary(fit)$r.squared)
}
results <- do.call(rbind,results); results$p_holm_30 <- p.adjust(results$p.value,method='holm')
write.csv(results,file.path(output_dir,'prs_pc_sex_sensitivity_contrasts.csv'),row.names=FALSE)
write.csv(do.call(rbind,details),file.path(output_dir,'prs_pc_sex_sensitivity_models.csv'),row.names=FALSE)
writeLines(c('Exploratory sensitivity analysis using saved prs_prsice; sex and ten genetic PCs included.',
  'Genotyping batch/array unavailable in the identified covariates: not adjusted.',
  'PC derivation, ancestry classification, genotype QC and relatedness provenance must be verified.',
  'No psychiatric-diagnosis covariates; Holm p-value adjustment is across 30 profile-vs-comparison tests.',
  'The 95% confidence intervals are pointwise and not multiplicity-adjusted.'),file.path(output_dir,'prs_pc_sex_sensitivity_README.txt'))
