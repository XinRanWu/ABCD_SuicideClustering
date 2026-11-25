# Refit selected manuscript results using original model specifications.
# No subject-level data are exported. Confidence intervals are pointwise Wald intervals.
library(lmerTest)
e<-new.env();load('Data_and_Scripts/abcd_SuicideData.RData',e)
im<-new.env();load('Data_and_Scripts/abcd_ImagingData.RData',im)
cov<-read.csv('../ABCD/abcd_covariables_full.csv',na.strings='NaN')
base<-merge(cov[,c(1:5,14:23)],e$abcd_envir.bl[,c(1,8:15,27:28,80)],all=TRUE)
out<-'Revision_Analysis_20261122/results_principal_imaging_20261125';dir.create(out,showWarnings=FALSE,recursive=TRUE)
results<-list()
for (modality in c('tbv','rsfc')) {
 brain<-if(modality=='tbv')im$abcd_smrip else im$abcd_rsfc
 dat<-merge(base,brain)
 dat$eventname<-factor(dat$eventname,levels=c('baseline_year_1_arm_1','2_year_follow_up_y_arm_1'),ordered=TRUE)
 variables<-if(modality=='tbv') 'smri_vol_scs_wholeb' else c('rsfmri_c_ngd_fo_ngd_sa','rsfmri_c_ngd_cgc_ngd_vta','rsfmri_c_ngd_n_ngd_sa')
 for(v in variables){
  stopifnot(v%in%names(dat))
  rhs<-if(modality=='tbv') 'clust*eventname+sex+income+edu+race_6level+hisp+Pubertal+anthro_bmi_calc+(1|mri_info_deviceserialnumber)+(1|subjectkey)' else 'clust*eventname+income+edu+sex+race_6level+hisp+Pubertal+anthro_bmi_calc+(1|mri_info_deviceserialnumber)+rsfmri_c_ngd_ntpoints+(1|subjectkey)'
  fit<-lmer(as.formula(paste(v,'~',rhs)),data=dat)
  co<-as.data.frame(coef(summary(fit)));co$term<-rownames(co);co$variable<-v;co$modality<-modality;co$n<-nobs(fit);co$singular<-isSingular(fit)
  co$ci95_low<-co$Estimate-qt(.975,df=co$df)*co$`Std. Error`;co$ci95_high<-co$Estimate+qt(.975,df=co$df)*co$`Std. Error`
  results[[length(results)+1]]<-co
  cat(v,'n',nobs(fit),'singular',isSingular(fit),'\n')
 }
}
x<-do.call(rbind,results);write.csv(x,file.path(out,'principal_imaging_model_coefficients.csv'),row.names=FALSE)
writeLines(c('Original formula includes scanner serial-number random intercept and subject random intercept.',
 'Original RSFC covariate is rsfmri_c_ngd_ntpoints; confirm whether this represents motion before describing it as framewise displacement.',
 'Confidence intervals use t distribution and Satterthwaite degrees of freedom and are pointwise; no multiplicity adjustment.',
 'Original control group and profile label alignment, data versions and recorded numbers should be checked against S7–S9.'),file.path(out,'notes.txt'))
