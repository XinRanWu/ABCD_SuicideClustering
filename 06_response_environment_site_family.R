# R1.5 / R2.5: environmental sensitivity models with site and family intercepts.
# Run from the project root. The continuous-score analysis follows Table S22's
# age/sex/race adjustment. Ordinal outcomes and unavailable items are logged.
library(lme4)
library(lmerTest)
library(emmeans)
data_dir <- Sys.getenv('ABCD_DATA_DIR', unset='../ABCD')
prepared_dir <- Sys.getenv('ABCD_PREPARED_DIR', unset='Data_and_Scripts')
out <- Sys.getenv('RESPONSE_OUTPUT_DIR', unset='Revision_Analysis_20261122/results_environment_site_family')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
env <- new.env(); load(file.path(prepared_dir,'abcd_SuicideData.RData'),envir=env)
base <- env$abcd_envir.bl
covars <- read.csv(file.path(data_dir,'abcd_covariables_full.csv'),na.strings=c('','NaN'))
covars <- covars[covars$eventname=='baseline_year_1_arm_1',c('subjectkey','site_id_l','rel_family_id')]
stopifnot(!anyDuplicated(base$subjectkey), !anyDuplicated(covars$subjectkey))
dat <- merge(base,covars,by='subjectkey',all.x=TRUE,suffixes=c('_original',''))
dat$clust <- factor(dat$clust,levels=0:5)
dat$site_id_l <- factor(dat$site_id_l); dat$rel_family_id <- factor(dat$rel_family_id)
dat$sex[trimws(as.character(dat$sex))==''] <- NA
dat$race_6level[trimws(as.character(dat$race_6level))==''] <- NA
dat$sex <- factor(dat$sex); dat$race_6level <- factor(dat$race_6level)
dat <- dat[!is.na(dat$clust),]

# Continuous scales (including multi-item total/mean scores). Income, education,
# binary/ordinal items and raw variables not present are NOT silently treated as
# continuous; they are listed in analysis_scope.csv for further handling.
continuous <- c('nsc_p_ss_mean_3_items','fes_p_ss_fc','fes_y_ss_fc',
  'crpbi_y_ss_parent','crpbi_y_ss_caregiver','pmq_y_ss_mean','psb_y_ss_mean',
  'srpf_y_ss_ses','srpf_y_ss_iiss','srpf_y_ss_dfs','Pubertal',
  'sds_p_ss_total','screen1_y','screen2_y','screen3_y','screen4_y','screen5_y','screen6_y')
ordinal <- c('income','edu','screen13_y','screen14_y')
unavailable <- c('kbi_p_conflict','kbi_p_grades_in_school','kbi_p_how_well_c_school',
  'physical_activity1_y','physical_activity2_y','kbi_p_c_bully','kbi_p_c_gay',
  'kbi_p_c_live_full_time','kbi_p_c_trans','ksads_ptsd_raw_754_p','ksads_ptsd_raw_755_p',
  'ksads_ptsd_raw_756_p','ksads_ptsd_raw_757_p','ksads_ptsd_raw_766_p',
  'ksads_ptsd_raw_769_p','ksads_ptsd_raw_770_p')
scope <- data.frame(outcome=c(continuous,ordinal,unavailable),
  intended=c(rep('continuous mixed model',length(continuous)),rep('ordinal: not fitted',length(ordinal)),rep('not found in prepared data: not fitted',length(unavailable))))
scope$available <- scope$outcome %in% names(dat)
write.csv(scope,file.path(out,'analysis_scope.csv'),row.names=FALSE)

results <- list(); diagnostics <- list()
for (outcome in continuous) {
  if (!(outcome %in% names(dat))) next
  vars <- c(outcome,'clust','sex','interview_age','race_6level','site_id_l','rel_family_id')
  sub <- dat[complete.cases(dat[,vars]),vars,drop=FALSE]
  sub <- droplevels(sub)
  cat('Fitting',outcome,'n=',nrow(sub),'\n')
  fit <- tryCatch(lmer(as.formula(paste(outcome,'~ clust + sex + interview_age + race_6level + (1|site_id_l) + (1|rel_family_id)')),
                       data=sub, REML=FALSE, control=lmerControl(optimizer='bobyqa')),error=function(e)e)
  if (inherits(fit,'error')) {
    diagnostics[[outcome]] <- data.frame(outcome=outcome,n=nrow(sub),status='ERROR',singular=NA,detail=conditionMessage(fit)); next
  }
  if (!('0' %in% levels(sub$clust))) stop('Comparison group missing for ',outcome)
  tab <- as.data.frame(summary(contrast(emmeans(fit,~clust,lmer.df='asymptotic'),method='trt.vs.ctrl',ref=1,adjust='none'),infer=c(TRUE,TRUE),adjust='none'))
  tab$outcome <- outcome; tab$n <- nrow(sub); tab$standardized_estimate <- tab$estimate/sigma(fit)
  results[[outcome]] <- tab
  diagnostics[[outcome]] <- data.frame(outcome=outcome,n=nrow(sub),status='OK',singular=isSingular(fit),
    detail=paste('sites',nlevels(sub$site_id_l),'families',nlevels(sub$rel_family_id), 'convergence',paste(fit@optinfo$conv$lme4$messages,collapse='; ')))
  write.csv(do.call(rbind,results),file.path(out,'continuous_contrasts_unadjusted_interim.csv'),row.names=FALSE)
  write.csv(do.call(rbind,diagnostics),file.path(out,'model_diagnostics_interim.csv'),row.names=FALSE)
}
if (length(results)) {
  result <- do.call(rbind,results); result$p_fdr_all_successful <- p.adjust(result$p.value,method='BH')
  write.csv(result,file.path(out,'continuous_contrasts_FDR.csv'),row.names=FALSE)
}
write.csv(do.call(rbind,diagnostics),file.path(out,'model_diagnostics.csv'),row.names=FALSE)
writeLines(c('Exploratory sensitivity: 18 continuous outcomes x 5 profile-vs-comparison contrasts; BH correction over successfully fitted contrasts.',
 'Raw estimates and pointwise 95% CIs in outcome units; standardized estimates divide by model residual SD.',
 'Models adjust for baseline sex, age and race and include site and family random intercepts.',
 'Income, education and screen13/14 are ordinal and not included; other listed outcomes are absent from the saved environment object.',
 'FDR family is a sensitivity-analysis choice, not a reconstruction of the original correction. Verify singularity, convergence, outcome distributions, missingness and group label alignment.'),file.path(out,'analysis_notes.txt'))
