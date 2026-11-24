# Manuscript neuroimaging: cortical thickness, area, subcortical volume, RSFC, DTI.
# Run after 01. This is extracted legacy code, not yet validated against all manuscript tables.
if (!requireNamespace("lmerTest",quietly=TRUE) || !requireNamespace("stringr",quietly=TRUE) || !requireNamespace("reshape2",quietly=TRUE)) stop("Install lmerTest, stringr and reshape2")
library(lmerTest); library(stringr); library(reshape2)
imaging_file <- file.path(prepared_dir,"abcd_ImagingData.RData")
if (!file.exists(imaging_file)) stop("Missing prepared imaging data")
covariates_file <- file.path(data_dir,"abcd_covariables_full.csv")
if (!file.exists(covariates_file)) stop("Set ABCD_DATA_DIR to directory containing abcd_covariables_full.csv")
imaging_env <- new.env(); load(imaging_file,envir=imaging_env)
abcd_smrip <- imaging_env$abcd_smrip; abcd_rsfc <- imaging_env$abcd_rsfc; abcd_tract <- imaging_env$abcd_tract
abcd_covariables <- read.csv(covariates_file,na.strings="NaN")
# Table S7 includes 14 subcortical ROIs. The previous sensitivity script selected 12;
# two indices have been added here and must be checked against the original table.
abcd_smrip <- abcd_smrip[,-1 * which(stringr::str_detect(names(abcd_smrip), "lesion"))]
abcd_smrip <- abcd_smrip[,-1 * which(stringr::str_detect(names(abcd_smrip), "wmhint"))]


abcd_BrainVar <- merge(abcd_covariables[, c(1:5, 14:23)], abcd_envir.bl[, c(1, 8:15, 27:28, 80)], all = T)
# abcd_BrainVar <- merge(abcd_BrainVar, abcd_tbss01[, c(1:2, 10:12)], all = T)
abcd_BrainVar <- merge(abcd_BrainVar, abcd_smrip)
abcd_BrainVar$eventname <- factor(abcd_BrainVar$eventname,
                                  levels = c("baseline_year_1_arm_1", "2_year_follow_up_y_arm_1"),ordered = T)
var_dic <-  c("eventname.L","clust1","clust2","clust3","clust4","clust5",
              "clust1:eventname.L","clust2:eventname.L","clust3:eventname.L",
              "clust4:eventname.L","clust5:eventname.L")
n <- length(var_dic)
fit.parameter.brain <- data.frame()
brain_dic <- names(abcd_smrip)[c(3:144, 287:321)] # 
aseg_id <-  c(6,23,19,13,3,25,28,8,30,17,27,14,9,26,7,24,10,11,12,35,33,34,32,31)
aseg_id.2nd <- c(1,2,15,16,8,6,13,14,4,11,12,7)
smrip_id <- c(which(str_detect(brain_dic, "_thick_cdk_"))[1:68],
              which(str_detect(brain_dic, "_area_cdk_"))[1:68],
              which(str_detect(brain_dic, "_scs_"))[aseg_id][aseg_id.2nd])
brain_dic2 <- c(brain_dic[smrip_id], "smri_vol_scs_aal", "smri_vol_scs_aar")

for (i in 1:length(brain_dic2)) {
  f <-  as.formula(paste(brain_dic2[i],
                         "~1+clust+eventname+clust*eventname+sex+income+edu+",
                         "race_6level+hisp+Pubertal+anthro_bmi_calc+smri_vol_scs_wholeb+",
                         "(1|mri_info_deviceserialnumber)+(1|subjectkey)",sep = ""))
  fit.lme.brain <- lmerTest::lmer(f, abcd_BrainVar)
  fit.lme.brain.summary <- summary(fit.lme.brain)
  fit.parameter.brain[((i - 1) * n + 1):(i * n), "var"] <- rep(brain_dic2[i], n)
  fit.parameter.brain[((i - 1) * n + 1):(i * n), "contrast"] <- var_dic
  fit.parameter.brain[((i - 1) * n + 1):(i * n), c("t", "df", "p")] <-
    fit.lme.brain.summary$coefficients[var_dic, c("t value", "df", "Pr(>|t|)")]
  print(paste(brain_dic2[i], " linear mixed effect is done. no. ", i, sep = ""))
}

i = length(brain_dic2) + 1
f <- as.formula(paste("smri_vol_scs_wholeb~1+clust+eventname+clust*eventname+sex+income+edu",
                      "+race_6level+hisp+Pubertal+anthro_bmi_calc+(1|mri_info_deviceserialnumber)+",
                      "(1|subjectkey)",sep = ""))
fit.lme.brain <- lmerTest::lmer(f, abcd_BrainVar)
fit.lme.brain.summary <- summary(fit.lme.brain)
fit.parameter.brain[((i - 1) * n + 1):(i * n), "var"] <-rep("smri_vol_scs_wholeb", n)
fit.parameter.brain[((i - 1) * n + 1):(i * n), "contrast"] <- var_dic
fit.parameter.brain[((i - 1) * n + 1):(i * n), c("t", "df", "p")] <-
  fit.lme.brain.summary$coefficients[var_dic, c("t value", "df", "Pr(>|t|)")]

fit.smrip <- fit.parameter.brain
fit.smrip[fit.smrip$var == "smri_vol_scs_wholeb",]

fit.smrip$var<-factor(fit.smrip$var,ordered = TRUE,levels = c(brain_dic2,"smri_vol_scs_wholeb"))
fit.smrip.dcast.t <- dcast(data = fit.smrip,var ~ contrast,value.var = "t")
fit.smrip.dcast.t$item <- rep("t",nrow(fit.smrip.dcast.t))
fit.smrip.dcast.df <- dcast(data = fit.smrip,var ~ contrast,value.var = "df")
fit.smrip.dcast.df$item <- rep("df",nrow(fit.smrip.dcast.df))
fit.smrip.dcast.p <- dcast(data = fit.smrip,var ~ contrast,value.var = "p")
fit.smrip.dcast.p$item <- rep("p",nrow(fit.smrip.dcast.p))
fit.smrip.dcast <- rbind(fit.smrip.dcast.t,fit.smrip.dcast.df,fit.smrip.dcast.p)
fit.smrip.dcast <- fit.smrip.dcast[,c("var", "item", "clust1","clust2","clust3","clust4","clust5",
                                      "clust1:eventname.L","clust2:eventname.L","clust3:eventname.L",
                                      "clust4:eventname.L","clust5:eventname.L","eventname.L")]

write.table(fit.smrip.dcast,file = "fit_smrip_noCryIQ.csv",na = "",sep = ",",row.names = FALSE)

dk_id <- formatC(c(2:4, 6:36), width = 4, flag = '0')
# aseg$data$region[3:26] # aseg$data$hemi[3:26]
# Resting-state functional connectivity: original Part3 lines 416-455.
# RSFC
abcd_BrainVar <- merge(abcd_covariables[, c(1:5, 14:23)], abcd_envir.bl[, c(1, 8:15, 27:28, 80)], all = T)
abcd_BrainVar <- merge(abcd_BrainVar, abcd_rsfc)
abcd_BrainVar$eventname <- factor(
  abcd_BrainVar$eventname,levels = c("baseline_year_1_arm_1","2_year_follow_up_y_arm_1"), ordered = T)
var_dic <-  c("eventname.L","clust1","clust2","clust3","clust4","clust5",
              "clust1:eventname.L","clust2:eventname.L","clust3:eventname.L",
              "clust4:eventname.L","clust5:eventname.L")
n <- length(var_dic)
fit.parameter.brain <- data.frame()
brain_dic <- names(abcd_rsfc)[4:172]# ncol(abcd_rsfc)
for (i in 1:length(brain_dic)) {
  f <- as.formula(paste(brain_dic[i],
        "~1+clust+eventname+clust*eventname+income+edu+",# nihtbx_totalcomp_uncorrected+
        "sex+race_6level+hisp+Pubertal+anthro_bmi_calc+",
        "(1|mri_info_deviceserialnumber)+rsfmri_c_ngd_ntpoints+(1|subjectkey)",sep = ""))
  fit.lme.brain <- lmerTest::lmer(f, abcd_BrainVar)
  fit.lme.brain.summary <- summary(fit.lme.brain)
  fit.parameter.brain[((i - 1) * n + 1):(i * n), "var"] <- rep(brain_dic[i], n)
  fit.parameter.brain[((i - 1) * n + 1):(i * n), "contrast"] <- var_dic
  fit.parameter.brain[((i - 1) * n + 1):(i * n), c("t", "df", "p")] <-
    fit.lme.brain.summary$coefficients[var_dic, c("t value", "df", "Pr(>|t|)")]
  print(paste(brain_dic[i], " linear mixed effect is done. no. ", i, sep = ""))
}
fit.rsfc <- fit.parameter.brain
fit.rsfc.dcast.t <- dcast(data = fit.rsfc,var ~ contrast,value.var = "t")
fit.rsfc.dcast.t$item <- rep("t",nrow(fit.rsfc.dcast.t))
fit.rsfc.dcast.df <- dcast(data = fit.rsfc,var ~ contrast,value.var = "df")
fit.rsfc.dcast.df$item <- rep("df",nrow(fit.rsfc.dcast.df))
fit.rsfc.dcast.p <- dcast(data = fit.rsfc,var ~ contrast,value.var = "p")
fit.rsfc.dcast.p$item <- rep("p",nrow(fit.rsfc.dcast.p))
fit.rsfc.dcast <- rbind(fit.rsfc.dcast.t,fit.rsfc.dcast.df,fit.rsfc.dcast.p)
fit.rsfc.dcast <- fit.rsfc.dcast[,c("var", "item", "clust1","clust2","clust3","clust4","clust5",
                                    "clust1:eventname.L","clust2:eventname.L","clust3:eventname.L",
                                    "clust4:eventname.L","clust5:eventname.L","eventname.L")]
fit.rsfc.dcast$var <- factor(fit.rsfc.dcast$var,ordered = TRUE,levels = brain_dic)

fit.parameter.brain[fit.parameter.brain$p < 0.00025,]
write.table(fit.rsfc.dcast,file = "fit_rsfc.csv",na = "",sep = ",",row.names = FALSE)

# Diffusion tensor imaging: original Part3 lines 499-548.
# Diffusion Data
abcd_BrainVar <-  merge(abcd_covariables[, c(1:5, 14:23)], abcd_envir.bl[, c(1, 8:15, 27:28, 80)], all = T)
# abcd_BrainVar <-  merge(abcd_BrainVar, abcd_tbss01[abcd_tbss01$eventname == "baseline_year_1_arm_1", c(1, 12)], all = T)
abcd_BrainVar <- merge(abcd_BrainVar, abcd_tract)
abcd_BrainVar$eventname <- factor(abcd_BrainVar$eventname,levels = c("baseline_year_1_arm_1","2_year_follow_up_y_arm_1"),ordered = T)
var_dic <-  c("eventname.L","clust1","clust2","clust3","clust4","clust5",
              "clust1:eventname.L","clust2:eventname.L","clust3:eventname.L",
              "clust4:eventname.L","clust5:eventname.L")
n <- length(var_dic)
fit.parameter.brain <- data.frame()
brain_dic <- names(abcd_tract)[4:151]# +nihtbx_totalcomp_uncorrected
for (i in 1:length(brain_dic)) {
  f <- as.formula(paste(brain_dic[i],
        "~1+clust+eventname+clust*eventname+edu+income+",
        "sex+race_6level+hisp+Pubertal+anthro_bmi_calc+",
        "(1|mri_info_deviceserialnumber)+dmdtifp1_1183+(1|subjectkey)",sep = ""))
  fit.lme.brain <- lmerTest::lmer(f, abcd_BrainVar)
  fit.lme.brain.summary <- summary(fit.lme.brain)
  fit.parameter.brain[((i - 1) * n + 1):(i * n), "var"] <- rep(brain_dic[i], n)
  fit.parameter.brain[((i - 1) * n + 1):(i * n), "contrast"] <- var_dic
  fit.parameter.brain[((i - 1) * n + 1):(i * n), c("t", "df", "p")] <-
    fit.lme.brain.summary$coefficients[var_dic, c("t value", "df", "Pr(>|t|)")]
  print(paste(brain_dic[i], " linear mixed effect is done. no. ", i, sep = ""))
}

fit.parameter.brain[fit.parameter.brain$p < .0025,]

fit.tract <- fit.parameter.brain
fit.tract.dcast.t <- dcast(data = fit.tract,var ~ contrast,value.var = "t")
fit.tract.dcast.t$item <- rep("t",nrow(fit.tract.dcast.t))
fit.tract.dcast.df <- dcast(data = fit.tract,var ~ contrast,value.var = "df")
fit.tract.dcast.df$item <- rep("df",nrow(fit.tract.dcast.df))
fit.tract.dcast.p <- dcast(data = fit.tract,var ~ contrast,value.var = "p")
fit.tract.dcast.p$item <- rep("p",nrow(fit.tract.dcast.p))
fit.tract.dcast <- rbind(fit.tract.dcast.t,fit.tract.dcast.df,fit.tract.dcast.p)


dti_dic <-  c('R Fx','L Fx','R CgC','L CgC','R CgH','L CgH','R CST','L CST','R ATR','L ATR',
              'R UNC','L UNC','R ILF','L ILF','R IFO','L IFO','Fmaj','Fmin','CC','R SLF',
              'L SLF','R tSLF','L tSLF','R pSLF','L pSLF','R SCS','L SCS','R fSCS','L fSCS',
              'R pSCS','L pSCS','R SIFC','L SIFC','R IFSFC','L IFSFC','R Fx.no Fim','L Fx.no Fim'  )

dti_labels <- data.frame(dic = rep(dti_dic,4),var = brain_dic,
                         mod = t(matrix(rep(c("FA", "LD", "RD", "TD"),37),4,37))[1:(4*37)])

fit.tract.dcast <- merge(fit.tract.dcast,dti_labels)
fit.tract.dcast <- fit.tract.dcast[,c("var","dic","mod","item","clust1","clust2","clust3","clust4","clust5",
                                    "clust1:eventname.L","clust2:eventname.L","clust3:eventname.L",
                                    "clust4:eventname.L","clust5:eventname.L","eventname.L")]
fit.tract.dcast$var <- factor(fit.tract.dcast$var,ordered = TRUE,levels = brain_dic)
