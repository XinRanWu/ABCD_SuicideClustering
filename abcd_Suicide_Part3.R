########################################################
#Statistic
########################################################

var_id <- c(19, 26:78, 193:200)
Contrast_id = c("1-0","2-0","3-0","4-0","5-0",
                "2-1","3-1","4-1","5-1","3-2",
                "4-2","5-2","4-3","5-3","5-4")
fit.parameter <- data.frame(Var = names(abcd_envir.bl)[var_id],
                            F = rep(0, length(var_id)),p = rep(0, length(var_id)),
                            DF = rep(0, length(var_id)))
fit.HSD <- data.frame(Var = rep(names(abcd_envir.bl)[var_id], length(Contrast_id)),
    Contrast = Contrast_id[gl(length(Contrast_id), length(var_id))],
    diff = rep(0, length(Contrast_id) * length(var_id)),
    p = rep(0, length(Contrast_id) * length(var_id)))

for (i in 1:length(var_id)) {
  f <-
    as.formula(paste(names(abcd_envir.bl)[var_id[i]], "~clust", sep = ""))
  fit <- aov(f, abcd_envir.bl)
  fit.parameter$F[i] <- summary(fit)[[1]]$`F value`[1]
  fit.parameter$p[i] <- summary(fit)[[1]]$`Pr(>F)`[1]
  fit.parameter$DF[i] <- summary(fit)[[1]]$Df[2]
  for (j in 1:length(Contrast_id)) {
    fit.HSD$diff[(j - 1) * length(var_id) + i] <-
      TukeyHSD(fit)[[1]][j, "diff"]
    fit.HSD$p[(j - 1) * length(var_id) + i] <-
      TukeyHSD(fit)[[1]][j, "p adj"]
  }
  
  print(paste(names(abcd_envir.bl)[var_id[i]], " ANOVA is done. no. ", i, sep = ""))
}


fit.HSD[which(fit.HSD$p < .0025 / 15),]
fit.HSD[fit.HSD$Contrast == "5-0",][which(fit.HSD[fit.HSD$Contrast == "5-0",]$p <
                                            .025 / 15),]


fit.HSD.dcast.diff <- dcast(data = fit.HSD,Var ~ Contrast,value.var = "diff")
fit.HSD.dcast.diff$item <- rep("diff",nrow(fit.HSD.dcast.diff))
fit.HSD.dcast.p <- dcast(data = fit.HSD,Var ~ Contrast,value.var = "p")
fit.HSD.dcast.p$item <- rep("p",nrow(fit.HSD.dcast.p))
fit.HSD.dcast <- rbind(fit.HSD.dcast.diff,fit.HSD.dcast.p)
# write.csv(fit.HSD.dcast,file = "fit_HSD_anova.csv")

profile <-
  rbind(fit.HSD[fit.HSD$Contrast == "1-0", "diff"], fit.HSD[fit.HSD$Contrast ==
                                                              "2-0", "diff"],
        fit.HSD[fit.HSD$Contrast == "3-0", "diff"], fit.HSD[fit.HSD$Contrast ==
                                                              "4-0", "diff"],
        fit.HSD[fit.HSD$Contrast == "5-0", "diff"])
profile.p <-
  rbind(fit.HSD[fit.HSD$Contrast == "1-0", "p"], fit.HSD[fit.HSD$Contrast ==
                                                           "2-0", "p"],
        fit.HSD[fit.HSD$Contrast == "3-0", "p"], fit.HSD[fit.HSD$Contrast ==
                                                           "4-0", "p"],
        fit.HSD[fit.HSD$Contrast == "5-0", "p"])
colnames(profile) <- names(abcd_envir.bl)[var_id]
row.names(profile) <-
  c("Subtype 1", "Subtype 2", "Subtype 3", "Subtype 4", "Subtype 5")
colnames(profile.p) <- names(abcd_envir.bl)[var_id]
row.names(profile.p) <-
  c("Subtype 1", "Subtype 2", "Subtype 3", "Subtype 4", "Subtype 5")
# chi test

var_id <- c(4, 14, 15, 24)
fit.xsq.parameter <- data.frame(
  Var = names(abcd_envir.bl)[var_id],
  Xsq = rep(0, length(var_id)),
  p = rep(0, length(var_id)),
  DF = rep(0, length(var_id))
)
fit.xsq.posthoc.parameter <- data.frame()
idx <- 0
for (i in 1:length(var_id)) {
  f <-
    as.formula(paste("~", names(abcd_envir.bl)[var_id[i]], "+clust"))
  fit.xsq <- chisq.test(xtabs(f, data = abcd_envir.bl))
  fit.xsq.parameter$Var[i] <- names(abcd_envir.bl)[var_id[i]]
  fit.xsq.parameter$Xsq[i] <- fit.xsq$statistic
  fit.xsq.parameter$p[i] <- fit.xsq$p.value
  fit.xsq.parameter$DF[i] <- fit.xsq$parameter
  fit.xsq.posthoc <- chisq.posthoc.test(xtabs(f, data = abcd_envir.bl),method = "fdr")
  for (k in 1:(nrow(fit.xsq.posthoc)/2)){
    idx <- idx + 1
    fit.xsq.posthoc.parameter[(1 + (idx - 1)*6):(idx*6),"Name"] <- 
      rep(names(abcd_envir.bl)[var_id[i]],6)
    
    fit.xsq.posthoc.parameter[(1 + (idx - 1)*6):(idx*6),"Cluster"] <- 
      c("0","1","2","3","4","5")
    
    fit.xsq.posthoc.parameter[(1 + (idx - 1)*6):(idx*6),"Dimension"] <- 
      rep(fit.xsq.posthoc[k*2,"Dimension"],6)
    
    fit.xsq.posthoc.parameter[(1 + (idx - 1)*6):(idx*6),"Residuals"] <- 
      as.numeric(fit.xsq.posthoc[which(fit.xsq.posthoc$Value=="Residuals")[k],
                                 c("0","1","2","3","4","5")])
    
    fit.xsq.posthoc.parameter[(1 + (idx - 1)*6):(idx*6),"p"] <- 
      as.character(fit.xsq.posthoc[which(fit.xsq.posthoc$Value=="p values")[k],
                                 c("0","1","2","3","4","5")])
  }
  print(paste(names(abcd_envir.bl)[var_id[i]], " X-square is done. no. ", i, sep = ""))
}
fit.xsq.posthoc.parameter
fit.xsq.parameter[fit.xsq.parameter$p < 0.001,]
# write.csv(fit.xsq.posthoc.parameter,file = "fit_xsq_posthoc_Apart.csv")

var_id <- c(81:192)
fit.xsq.parameter <- data.frame(
  Var = names(abcd_envir.bl)[var_id],
  Xsq = rep(0, length(var_id)),
  p = rep(0, length(var_id)),
  DF = rep(0, length(var_id))
)
abcd_envir.bl.noHC <- abcd_envir.bl[-1*which(abcd_envir.bl$clust == "0"),]
abcd_envir.bl.noHC$clust <- factor(abcd_envir.bl.noHC$clust, 
                                   ordered = TRUE, 
                                   levels = c("1","2","3","4","5"))
fit.xsq.posthoc.parameter <- data.frame()
idx <- 0
for (i in 1:length(var_id)) {
  f <-
    as.formula(paste("~", names(abcd_envir.bl)[var_id[i]], "+clust"))
  x <- xtabs(f, data = abcd_envir.bl.noHC)
  fit.xsq <- chisq.test(x)
  fit.xsq.parameter$Var[i] <- names(abcd_envir.bl.noHC)[var_id[i]]
  fit.xsq.parameter$Xsq[i] <- fit.xsq$statistic
  fit.xsq.parameter$p[i] <- fit.xsq$p.value
  fit.xsq.parameter$DF[i] <- fit.xsq$parameter
  if (nrow(x) != 1) {
    fit.xsq.posthoc <-
      chisq.posthoc.test(x, method = "fdr")
    for (k in 1:(nrow(fit.xsq.posthoc) / 2)) {
      idx <- idx + 1
      fit.xsq.posthoc.parameter[(1 + (idx - 1) * 5):(idx * 5), "Name"] <-
        rep(names(abcd_envir.bl)[var_id[i]], 5)
      
      fit.xsq.posthoc.parameter[(1 + (idx - 1) * 5):(idx * 5), "Cluster"] <-
        c("1", "2", "3", "4", "5")
      
      fit.xsq.posthoc.parameter[(1 + (idx - 1) * 5):(idx * 5), "Dimension"] <-
        rep(fit.xsq.posthoc[k * 2, "Dimension"], 5)
      
      fit.xsq.posthoc.parameter[(1 + (idx - 1) * 5):(idx * 5), "Residuals"] <-
        as.numeric(fit.xsq.posthoc[which(fit.xsq.posthoc$Value == "Residuals")[k],
                                   c("1", "2", "3", "4", "5")])
      
      fit.xsq.posthoc.parameter[(1 + (idx - 1) * 5):(idx * 5), "p"] <-
        as.character(fit.xsq.posthoc[which(fit.xsq.posthoc$Value == "p values")[k],
                                     c("1", "2", "3", "4", "5")])
    }
  }
  print(paste(names(abcd_envir.bl)[var_id[i]], " X-square is done. no. ", i, sep = ""))
}
fit.xsq.posthoc.parameter
View(fit.xsq.posthoc.parameter[
  fit.xsq.posthoc.parameter$Cluster=="5"&
    fit.xsq.posthoc.parameter$Dimension=="1",
  ])
fit.xsq.parameter[fit.xsq.parameter$p < 0.001,]
abcd_ksads.suicide.label <- read.csv2("ksads_suicidal_labels.csv")
fit.xsq.posthoc.parameter <- 
  merge(fit.xsq.posthoc.parameter,
        abcd_ksads.suicide.label,
        by.x = "Name",by.y = "Code")
fit.xsq.posthoc.parameter$Cluster <- 
  paste("Subtype",fit.xsq.posthoc.parameter$Cluster,sep = " ")

# write.csv(fit.xsq.posthoc.parameter,file = "fit_xsq_posthoc_Bpart.csv")

tabone <-
  CreateTableOne(data = abcd_envir.bl[, c(26, 28:79, 81:200)], vars = abcd_envir.bl$clust)
print(tabone)

abcd_envir.bl.factor <- abcd_envir.bl
abcd_envir.bl.factor[, c(4, 14, 15, 24)] <-
  lapply(abcd_envir.bl.factor[, c(4, 14, 15, 24)], factor, exclude =  "")
abcd_envir.bl.factor[, c(81:192)] <-
  lapply(abcd_envir.bl.factor[, c(81:192)],
         factor,
         ordered = TRUE,
         levels = c(0, 1))

tableOne <- CreateTableOne(
  vars = names(abcd_envir.bl.factor)
  [c(3, 4, 14, 15, 19, 24, 25, 27, 28, 179:192)],
  strata = c("clust"),
  data = abcd_envir.bl.factor
)
tableOne <-
  CreateTableOne(
    vars = names(abcd_envir.bl.factor)[c(81:178)],
    strata = c("clust"),
    data = abcd_envir.bl.factor
  )
tableOne
summary(tableOne)

########################################################
# Polygenic Score
########################################################

# abcd_prscs.clust<-merge(abcd_envir.bl[,c("subjectkey","clust")],prs_prscs)
abcd_prscs.clust <-
  merge(abcd_envir.bl[, c("subjectkey", "clust")], prs_prsice)

var_id <- c(3:length(abcd_prscs.clust))
Contrast_id = c("1-0","2-0","3-0","4-0","5-0",
                "2-1","3-1","4-1","5-1","3-2",
                "4-2","5-2","4-3","5-3","5-4")
fit.parameter <- data.frame(
  Var = names(abcd_prscs.clust)[var_id],
  F = rep(0, length(var_id)),
  p = rep(0, length(var_id)),
  DF = rep(0, length(var_id))
)
fit.HSD <-
  data.frame(
    Var = rep(names(abcd_prscs.clust)[var_id], length(Contrast_id)),
    Contrast = Contrast_id[gl(length(Contrast_id), length(var_id))],
    diff = rep(0, length(Contrast_id) * length(var_id)),
    p = rep(0, length(Contrast_id) * length(var_id))
  )

for (i in 1:length(var_id)) {
  f <-
    as.formula(paste(names(abcd_prscs.clust)[var_id[i]], "~clust", sep = ""))
  fit <- aov(f, abcd_prscs.clust)
  fit.parameter$F[i] <- summary(fit)[[1]]$`F value`[1]
  fit.parameter$p[i] <- summary(fit)[[1]]$`Pr(>F)`[1]
  fit.parameter$DF[i] <- summary(fit)[[1]]$Df[2]
  for (j in 1:length(Contrast_id)) {
    fit.HSD$diff[(j - 1) * length(var_id) + i] <-
      TukeyHSD(fit)[[1]][j, "diff"]
    fit.HSD$p[(j - 1) * length(var_id) + i] <-
      TukeyHSD(fit)[[1]][j, "p adj"]
  }
  print(paste(names(abcd_prscs.clust)[var_id[i]], " ANOVA is done. no. ", i, sep = ""))
}

fit.HSD[fit.HSD$Contrast == "1-0",][which(fit.HSD[fit.HSD$Contrast == "1-0",]$p <
                                            0.005),]
fit.HSD[which(fit.HSD$p < 0.025 / 15),]
fit.HSD.dcast.diff <- dcast(data = fit.HSD,Var ~ Contrast,value.var = "diff")
fit.HSD.dcast.diff$item <- rep("diff",nrow(fit.HSD.dcast.diff))
fit.HSD.dcast.p <- dcast(data = fit.HSD,Var ~ Contrast,value.var = "p")
fit.HSD.dcast.p$item <- rep("p",nrow(fit.HSD.dcast.p))
fit.HSD.dcast <- rbind(fit.HSD.dcast.diff,fit.HSD.dcast.p)
# write.csv(fit.HSD.dcast,file = "fit_HSD_anova_PRS.csv")



# Density Ridges Plots of PRS of ADHD
Colormap <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32)
ggplot(abcd_prscs.clust[-1 * which(is.na(abcd_prscs.clust$clust)),],
       aes(y = clust, x = adhd), fill = clust) +
  geom_density_ridges_gradient(aes(fill = ..x..), scale = 3, size = 0.3) +
  scale_fill_gradientn(colours = Colormap, name = "PRS of ADHD") +
  theme_minimal()

ggplot(abcd_prscs.clust[-1 * which(is.na(abcd_prscs.clust$clust)),],
       aes(x = clust,y = adhd,fill = clust,colour = clust,group = clust)) +
  geom_half_violin(side = "r",width = 1.6) + geom_boxplot(width = 0.3,fill = "white",colour = "black",side = "r",width = 2  ) +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) + theme_bw() + coord_flip()  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12))



########################################################
# Brain Statistic
########################################################

abcd_tbss01 <-
  read.abcd.table("F:/Acdamic/Data_ABCD/abcd_tbss01.txt")
abcd_tbss01 <-
  abcd_tbss01[, c(4, 9, which(stringr::str_detect(names(abcd_tbss01), "_uncorrected")))]






abcd_smrip <- abcd_smrip[,-1 * which(stringr::str_detect(names(abcd_smrip), "lesion"))]
abcd_smrip <- abcd_smrip[,-1 * which(stringr::str_detect(names(abcd_smrip), "wmhint"))]

abcd_covariables <- read.csv("abcd_covariables_full.csv", na.strings = "NaN")

abcd_BrainVar <- merge(abcd_covariables[, c(1:5, 14:23)], abcd_envir.bl[, c(1, 8:15, 27:28, 80)], all = T)
abcd_BrainVar <- merge(abcd_BrainVar, abcd_smrip)
abcd_BrainVar$eventname <- factor(abcd_BrainVar$eventname,
                                  levels = c("baseline_year_1_arm_1", "2_year_follow_up_y_arm_1"),ordered = T)
var_dic <-  c("eventname.L","clust1","clust2","clust3","clust4","clust5",
              "clust1:eventname.L","clust2:eventname.L","clust3:eventname.L",
              "clust4:eventname.L","clust5:eventname.L")
n <- length(var_dic)
fit.parameter.brain <- data.frame()
brain_dic <- names(abcd_smrip)[c(3:144, 287:321)] # +nihtbx_totalcomp_uncorrected
aseg_id <-  c(6,23,19,13,3,25,28,8,30,17,27,14,9,26,7,24,10,11,12,35,33,34,32,31)
aseg_id.2nd <- c(1,2,15,16,8,6,13,14,4,11,12,7)
smrip_id <- c(which(str_detect(brain_dic, "_thick_cdk_"))[1:68],
              which(str_detect(brain_dic, "_area_cdk_"))[1:68],
              which(str_detect(brain_dic, "_scs_"))[aseg_id][aseg_id.2nd])

brain_dic <- names(abcd_smrip)[c(3:321)] # +nihtbx_totalcomp_uncorrected
smrip_id <- c(which(str_detect(brain_dic, "_vol_cdk_"))[1:68],
              which(str_detect(brain_dic, "_vol_scs_"))[aseg_id][aseg_id.2nd])
brain_dic2 <- brain_dic[smrip_id]

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

write.table(fit.smrip.dcast,file = "fit_smrip.csv",na = "",sep = ",",row.names = FALSE)

dk_id <- formatC(c(2:4, 6:36), width = 4, flag = '0')
# aseg$data$region[3:26] # aseg$data$hemi[3:26]
# Brain Plots of Thickness & Area
var = "clust2"
t = fit.smrip.dcast[fit.smrip.dcast$item=="t",var]
p = fit.smrip.dcast[fit.smrip.dcast$item=="p",var]

brain_dic2[which(p.adjust(p[1:148], method = "fdr") < 0.05)]

# sig_id <- c(
#   which(p.adjust(fit.smrip.dcast[
#     fit.smrip.dcast$item=="p","clust2"][1:148], method = "fdr") < 0.05),
#   which(p.adjust(fit.smrip.dcast[
#     fit.smrip.dcast$item=="p","clust3"][1:148], method = "fdr") < 0.05),
#   which(p.adjust(fit.smrip.dcast[
#     fit.smrip.dcast$item=="p","clust5"][1:148], method = "fdr") < 0.05))
# fit.smrip.dcast.sig <- rbind(fit.smrip.dcast.t[sig_id,c("var","item",paste("clust",1:5,sep = ""))],
#                              fit.smrip.dcast.df[sig_id,c("var","item",paste("clust",1:5,sep = ""))],
#                              fit.smrip.dcast.p[sig_id,c("var","item",paste("clust",1:5,sep = ""))])

data = data.frame(roi = rep(dk_id, 4),t = t[1:(68*2)] * (p.adjust(p, method = "fdr")[1:(68*2)] < 0.05),
  stringsAsFactors = FALSE,Modal = c(rep("Thickness", 68), rep("Area", 68)),
  hemi = rep(c(rep("left", 34), rep("right", 34)), 2)) # *(p.adjust(p,method = "fdr")<0.05)

data %>%
  group_by(Modal) %>% ggseg(atlas = dk,colour = "black",
                            mapping = aes(fill = t),position = "stacked",size = .5) +
  facet_wrap( ~ Modal, ncol = 1) + theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
    na.value = "lightgrey", breaks = seq(-4.5, 4.5, 0.9),limits = c(-4.5, 4.5)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))


data <- data.frame(roi = aseg$data$region[3:26][aseg_id.2nd],
                  t = t[137:148]*(p.adjust(p[137:148],method = "fdr") < 0.05),
  stringsAsFactors = FALSE,hemi = aseg$data$hemi[3:26][aseg_id.2nd],
  label = aseg$data$label[3:26][aseg_id.2nd],Modal = rep("Volume", 14)) 

ggseg(data, atlas = aseg, colour = "black",mapping = aes(fill = t)) +
  theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],
    na.value = "lightgrey",breaks = seq(-4.5, 4.5, 0.9),limits = c(-4.5, 4.5)) +
  theme_brain() + theme(text = element_text(family = "Arial", size = 15))

# Density Ridges Plots of Whole brain Volume
f <-as.formula(paste("smri_vol_scs_wholeb~1+sex+income+edu",
      "+race_6level+hisp+Pubertal+anthro_bmi_calc+(1|mri_info_deviceserialnumber)",sep = ""))
fit.lme.brain <- lmerTest::lmer(f, abcd_BrainVar)
fit.lme.brain.summary <- summary(fit.lme.brain)
abcd_Clust4WBV<-data.frame(clust<-fit.lme.brain@frame$clust,WBV=fit.lme.brain.summary$residuals)

Colormap <- colorRampPalette(rev(brewer.pal(11, 'RdBu')))(32)
ggplot(abcd_Clust4WBV,aes(x = clust,y = WBV,fill = clust,colour = clust,group = clust)) +
  geom_violin(width = 1.5) + geom_boxplot(width = 0.1,fill = "white",colour = "black") +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) + theme_bw()  +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,size = 12
  )) + scale_y_continuous(limits = c(-2.5,2.5)) + coord_flip()




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

net_dic <-  c('AUD','CO','CP','DMN','DAN','FPN','N','RT','SMH','SMM','SN','VAN','VIS')
net_dic2 <-  c('crcxlh','thplh','cdelh','ptlh','pllh','bs','hplh','aglh','aalh',
               'vtdclh','crcxrh','thprh','cderh','ptrh','plrh','hprh','agrh','aarh','vtdcrh'  )

var = "clust2"
t.mat <- matrix(fit.rsfc[fit.rsfc$contrast == var,]$t[1:169], 13, 13)# [c(1:6,8:13),c(1:6,8:13)]
p.mat <- matrix(fit.rsfc[fit.rsfc$contrast == var,]$p[1:169], 13, 13)# [c(1:6,8:13),c(1:6,8:13)]

which(p.adjust(p.mat[which((1-lower.tri(t.mat)) == 1)],method = "fdr")<0.05)

row.names(t.mat) <- net_dic;colnames(t.mat) <- net_dic
row.names(p.mat) <- net_dic;colnames(p.mat) <- net_dic
corrplot::corrplot(t.mat,  p.mat = p.mat, is.corr = FALSE, method = "circle",  
                   type = "lower",sig.level = c(.001, .01, .05),pch.cex = .9, 
                   insig = "label_sig",  pch.col = "white", tl.col = "black", cl.lim = c(-6.5, 6.5))

p.mat.adjust <- matrix(p.adjust(p.mat,method = "fdr"), 13, 13)
row.names(p.mat.adjust) <- net_dic;colnames(p.mat.adjust) <- net_dic
corrplot::corrplot(t.mat,  p.mat = p.mat.adjust, is.corr = FALSE, method = "circle",  
                   type = "lower",pch.cex = .9, 
                   insig = "blank",  pch.col = "white", tl.col = "black", cl.lim = c(-6.5, 6.5))

corrplot::corrplot(t.mat[c("FPN","SN","CO","VAN"),c("FPN","SN","CO","VAN")],  
                   p.mat = p.mat.adjust[c("FPN","SN","CO","VAN"),c("FPN","SN","CO","VAN")], 
                   is.corr = FALSE, method = "circle",
                   type = "lower",pch.cex = .9, insig = "blank",  
                   pch.col = "white", tl.col = "black", cl.lim = c(-5.5, 5.5))

var = "clust5:eventname.L"
scs_id2 <- c(2,12,3,13,4,14,5,15,7,16,8,17,9,18)
t.mat <-  matrix(fit.rsfc[fit.rsfc$contrast == var,]$t[170:416], 19, 13)[scs_id2,]
p.mat <-  matrix(p.adjust(fit.rsfc[fit.rsfc$contrast == var,]$p[170:416], method = "fdr"), 19, 13)[scs_id2,]
row.names(t.mat) <- net_dic2[scs_id2];colnames(t.mat) <- net_dic
row.names(p.mat) <- net_dic2[scs_id2];colnames(p.mat) <- net_dic
corrplot::corrplot(t.mat,p.mat = p.mat,is.corr = FALSE, method = "circle",sig.level = c(.001, .01, .05),
                   col = brewer.pal(11, "RdBu")[11:1],pch.cex = .9,  insig = "label_sig",  
                   pch.col = "white",  tl.col = "black",cl.lim = c(-6.5, 6.5)) # col = brewer.pal(11, "RdBu")[11:1]

corrplot::corrplot(t.mat,p.mat = p.mat,is.corr = FALSE, method = "circle",
                   col = brewer.pal(11, "RdBu")[11:1],pch.cex = .9,  insig = "blank",  
                   pch.col = "white",  tl.col = "black",cl.lim = c(-6.5, 6.5)) # col = brewer.pal(11, "RdBu")[11:1]


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

write.table(fit.tract.dcast,file = "fit_tract.csv", na = "",sep = ",",row.names = FALSE)

var <- "clust3"
t.mat <- matrix(fit.parameter.brain$t[fit.parameter.brain$contrast == var], 37, 4)
p.mat <- matrix(p.adjust(fit.parameter.brain$p[fit.parameter.brain$contrast == var], method = "fdr"), 37, 4)
row.names(t.mat) <- dti_dic;colnames(t.mat) <- c("FA", "LD", "RD", "TD")
row.names(p.mat) <- dti_dic;colnames(p.mat) <- c("FA", "LD", "RD", "TD")
which(p.mat<0.05)

corrplot::corrplot(t(t.mat),  p.mat = t(p.mat),
                   is.corr = FALSE,  method = "circle",
                   sig.level = c(.001, .01, .05),  pch.cex = .9,
                   col = brewer.pal(11, "RdBu")[11:1], insig = "label_sig", 
                   pch.col = "white",  tl.col = "black",  cl.lim = c(-5.5, 5.5))

corrplot::corrplot(t(t.mat * (p.mat < 0.05)),
                   is.corr = FALSE,  method = "circle",  pch.cex = .9,
                   col = brewer.pal(11, "RdBu")[11:1],
                   pch.col = "white",  tl.col = "black",  cl.lim = c(-5.5, 5.5))

# RSI

abcd_BrainVar <- merge(abcd_covariables[, c(1:5, 14:23)], abcd_envir.bl[, c(1, 8:15, 27:28, 80)], all = T)
abcd_BrainVar <- merge(abcd_BrainVar, abcd_rsi)
abcd_BrainVar$eventname <- factor(abcd_BrainVar$eventname,
                                  levels = c("baseline_year_1_arm_1","2_year_follow_up_y_arm_1"),ordered = T)
var_dic <-  c("eventname.L","clust1","clust2","clust3","clust4","clust5",
              "clust1:eventname.L","clust2:eventname.L","clust3:eventname.L",
              "clust4:eventname.L","clust5:eventname.L")
n <- length(var_dic)
fit.parameter.brain <- data.frame()
brain_dic <- names(abcd_rsi)[str_detect(names(abcd_rsi),"_scs_")] # names(abcd_rsi)[c(6:611)] # +nihtbx_totalcomp_uncorrected
for (i in 1:length(brain_dic)) {
  f <- as.formula(paste(brain_dic[i],
        "~1+clust+eventname+clust*eventname+sex+income+edu+",
        "race_6level+hisp+Pubertal+anthro_bmi_calc+dmri_rsi_meanmotion+",
        "(1|mri_info_deviceserialnumber)+(1|subjectkey)",sep = ""))
  fit.lme.brain <- lmerTest::lmer(f, abcd_BrainVar)
  fit.lme.brain.summary <- summary(fit.lme.brain)
  fit.parameter.brain[((i - 1) * n + 1):(i * n), "var"] <- rep(brain_dic[i], n)
  fit.parameter.brain[((i - 1) * n + 1):(i * n), "contrast"] <- var_dic
  fit.parameter.brain[((i - 1) * n + 1):(i * n), c("t", "df", "p")] <-
    fit.lme.brain.summary$coefficients[var_dic, c("t value", "df", "Pr(>|t|)")]
  print(paste(brain_dic[i], " linear mixed effect is done. no. ", i, sep = ""))
}

aseg_id <-  c(6,23,19,13,3,25,28,8,30,17,27,14,9,26,7,24,10,11,12,35,33,34,32,31)
aseg_id.2nd <- c(1,2,15,16,8,6,13,14,4,11,12,7)


fit.rsi <- fit.parameter.brain
fit.rsi.dcast.t <- dcast(data = fit.rsi,var ~ contrast,value.var = "t")
fit.rsi.dcast.t$item <- rep("t",nrow(fit.rsi.dcast.t))
fit.rsi.dcast.df <- dcast(data = fit.rsi,var ~ contrast,value.var = "df")
fit.rsi.dcast.df$item <- rep("df",nrow(fit.rsi.dcast.df))
fit.rsi.dcast.p <- dcast(data = fit.rsi,var ~ contrast,value.var = "p")
fit.rsi.dcast.p$item <- rep("p",nrow(fit.rsi.dcast.p))
fit.rsi.dcast <- rbind(fit.rsi.dcast.t,fit.rsi.dcast.df,fit.rsi.dcast.p)
fit.rsi.dcast <- fit.rsi.dcast[,c("var","item",var_dic)]
fit.rsi.dcast <- fit.rsi.dcast[,c("var","item","clust1","clust2","clust3","clust4","clust5",
                                  "clust1:eventname.L","clust2:eventname.L","clust3:eventname.L",
                                  "clust4:eventname.L","clust5:eventname.L","eventname.L")]


write.table(fit.rsi.dcast,file = "fit_rsi.csv", na = "",sep = ",",row.names = FALSE)

mri_rsi_mode <- c("rsin0gm","rsin0s2gm","rsindgm","rsinds2gm","rsintgm","rsints2gm")

dk_id <- formatC(c(2:4, 6:36), width = 4, flag = '0')
aseg_id <-  c(6,23,19,13,3,25,28,8,30,17,27,14,9,26)
# Brain Plots of Thickness & Area
var = "clust1"
t = c(fit.rsi[fit.rsi$contrast == var & str_detect(fit.rsi$var, paste(mri_rsi_mode[1], "_cdk_", sep = "")),][1:68,]$t,
      fit.rsi[fit.rsi$contrast == var & str_detect(fit.rsi$var, paste(mri_rsi_mode[2], "_cdk_", sep = "")),][1:68,]$t,
      fit.rsi[fit.rsi$contrast == var & str_detect(fit.rsi$var, paste(mri_rsi_mode[3], "_cdk_", sep = "")),][1:68,]$t,
      fit.rsi[fit.rsi$contrast == var & str_detect(fit.rsi$var, paste(mri_rsi_mode[4], "_cdk_", sep = "")),][1:68,]$t,
      fit.rsi[fit.rsi$contrast == var & str_detect(fit.rsi$var, paste(mri_rsi_mode[5], "_cdk_", sep = "")),][1:68,]$t,
      fit.rsi[fit.rsi$contrast == var & str_detect(fit.rsi$var, paste(mri_rsi_mode[6], "_cdk_", sep = "")),][1:68,]$t)
p = c(fit.rsi[fit.rsi$contrast == var & str_detect(fit.rsi$var, paste(mri_rsi_mode[1], "_cdk_", sep = "")),][1:68,]$p,
      fit.rsi[fit.rsi$contrast == var & str_detect(fit.rsi$var, paste(mri_rsi_mode[2], "_cdk_", sep = "")),][1:68,]$p,
      fit.rsi[fit.rsi$contrast == var & str_detect(fit.rsi$var, paste(mri_rsi_mode[3], "_cdk_", sep = "")),][1:68,]$p,
      fit.rsi[fit.rsi$contrast == var & str_detect(fit.rsi$var, paste(mri_rsi_mode[4], "_cdk_", sep = "")),][1:68,]$p,
      fit.rsi[fit.rsi$contrast == var & str_detect(fit.rsi$var, paste(mri_rsi_mode[5], "_cdk_", sep = "")),][1:68,]$p,
      fit.rsi[fit.rsi$contrast == var & str_detect(fit.rsi$var, paste(mri_rsi_mode[6], "_cdk_", sep = "")),][1:68,]$p)

t.scs = fit.rsi[fit.rsi$contrast == var & str_detect(fit.rsi$var, "_scs_"),]$t
p.scs = fit.rsi[fit.rsi$contrast == var & str_detect(fit.rsi$var, "_scs_"),]$p

brain_dic[which(p.adjust(c(p.scs), method = "fdr") < 0.05)]



data = data.frame(roi = rep(dk_id, 12), t = t,stringsAsFactors = FALSE,
  Modal = c(rep("IN0", 68), rep("IN0_S2", 68), rep("IND0", 68), 
            rep("IND0_S2", 68), rep("INT0", 68), rep("INT0_S2", 68)),
  hemi = rep(c(rep("left", 34), rep("right", 34)), 6)) # *(p.adjust(p,method = "fdr")<0.05)
data %>%
  group_by(Modal) %>% ggseg(atlas = dk,colour = "black",position = "stacked",mapping = aes(fill = t)) +
  facet_wrap( ~ Modal, ncol = 2) + theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],na.value = "lightgrey",
    breaks = seq(-4.5, 4.5, 0.9),limits = c(-4.5, 4.5)) +  theme_brain()

data = data.frame(t = t * (p.adjust(p, method = "fdr") < 0.05),
  roi = rep(dk_id, 12),stringsAsFactors = FALSE,
  Modal = c(rep("IN0", 68), rep("IN0_S2", 68),
            rep("IND0", 68), rep("IND0_S2", 68),
            rep("INT0", 68), rep("INT0_S2", 68)),
  hemi = rep(c(rep("left", 34), rep("right", 34)), 6))

data %>%  group_by(Modal) %>% ggseg(atlas = dk,colour = "black",
                                    mapping = aes(fill = t),position = "stacked",size = .5) +
  facet_wrap( ~ Modal, ncol = 3) + theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],na.value = "lightgrey",
    breaks = seq(-4.5, 4.5, 0.9),limits = c(-4.5, 4.5)) +  theme_brain() +
  theme(text = element_text(family = "Arial", size = 15))


# MID
# anticipation of large reward versus neutral contrast
# anticipation of large loss versus neutral contrast
# reward positive versus negative feedback contrast
# loss positive versus negative feedback contrast

nuclei_id_mid <- c(6, 7, 8, 9, 13, 14, 16, 23:29)
con_id <- c(5, 8, 3, 4);surf_id <- c();scs_id <- c()

for (i in 1:length(con_id)) {
  surf_id <- c(surf_id, which(
    str_detect(names(abcd_MID), paste(mid_dic2[con_id[i]], '_b_cds_', sep = "")    )))
  scs_id <-  c(scs_id, which(str_detect(
    names(abcd_MID), paste(mid_dic[con_id[i]], '_b_scs_', sep = "")    )))
}

abcd_BrainVar <- merge(abcd_covariables[, c(1:5, 14:23)], abcd_envir.bl[, c(1, 8:15, 27:28, 80)], all = T)
abcd_BrainVar <- merge(abcd_BrainVar, abcd_MID[, c(1:5, surf_id, scs_id)])
abcd_BrainVar$eventname <- factor(abcd_BrainVar$eventname,
                                  levels = c("baseline_year_1_arm_1", "2_year_follow_up_y_arm_1"),ordered = T)
var_dic <-  c("eventname.L","clust1","clust2","clust3","clust4","clust5",
              "clust1:eventname.L","clust2:eventname.L","clust3:eventname.L",
              "clust4:eventname.L","clust5:eventname.L")
n <- length(var_dic)
fit.parameter.brain <- data.frame()
brain_dic <- names(abcd_MID)[c(surf_id, scs_id)]
for (i in 1:length(brain_dic)) {
  f <- as.formula(paste(brain_dic[i],
        "~1+clust+eventname+clust*eventname+sex+income+edu+",
        "race_6level+hisp+Pubertal+anthro_bmi_calc+tfmri_mid_all_b_meanmotion+",
        "(1|mri_info_deviceserialnumber)+(1|subjectkey)", sep = ""))
  fit.lme.brain <- lmerTest::lmer(f, abcd_BrainVar)
  fit.lme.brain.summary <- summary(fit.lme.brain)
  fit.parameter.brain[((i - 1) * n + 1):(i * n), "var"] <-  rep(brain_dic[i], n)
  fit.parameter.brain[((i - 1) * n + 1):(i * n), "contrast"] <-  var_dic
  fit.parameter.brain[((i - 1) * n + 1):(i * n), c("t", "df", "p")] <-
    fit.lme.brain.summary$coefficients[var_dic, c("t value", "df", "Pr(>|t|)")]
  print(paste(brain_dic[i], " linear mixed effect is done. no. ", i, sep = ""))
}

fit.parameter.brain[fit.parameter.brain$p < .025,]
fit.mid <- fit.parameter.brain

fit.mid.dcast.t <- dcast(data = fit.mid,var ~ contrast,value.var = "t")
fit.mid.dcast.t$item <- rep("t",nrow(fit.mid.dcast.t))
fit.mid.dcast.df <- dcast(data = fit.mid,var ~ contrast,value.var = "df")
fit.mid.dcast.df$item <- rep("df",nrow(fit.mid.dcast.df))
fit.mid.dcast.p <- dcast(data = fit.mid,var ~ contrast,value.var = "p")
fit.mid.dcast.p$item <- rep("p",nrow(fit.mid.dcast.p))
fit.mid.dcast <- rbind(fit.mid.dcast.t,fit.mid.dcast.df,fit.mid.dcast.p)
fit.mid.dcast <- fit.mid.dcast[,c("var","item","clust1","clust2","clust3","clust4","clust5",
                                  "clust1:eventname.L","clust2:eventname.L","clust3:eventname.L",
                                  "clust4:eventname.L","clust5:eventname.L","eventname.L")]


write.table(fit.mid,file = "fit_mid.csv",na = "",sep = ",",row.names = FALSE)

dk_id <- formatC(c(2:4, 6:36), width = 4, flag = '0')
aseg_id <-  c(6,23,19,13,3,25,28,8,30,17,27,14,9,26,7,24,10,11,12,35,33,34,32,31)

var = "clust4"
con_cds = mid_dic2[con_id[4]]
con_scs = mid_dic[con_id[4]]
t = c(fit.mid[fit.mid$contrast == var & str_detect(fit.mid$var, paste(con_cds, "_b_cds_", sep = "")),][1:68,]$t)
p = c(fit.mid[fit.mid$contrast == var & str_detect(fit.mid$var, paste(con_cds, "_b_cds_", sep = "")),][1:68,]$p)
t.scs = fit.mid[fit.mid$contrast == var & str_detect(fit.mid$var, paste(con_scs, "_b_scs_", sep = "")),]$t[aseg_id]
p.scs = fit.mid[fit.mid$contrast == var & str_detect(fit.mid$var, paste(con_scs, "_b_scs_", sep = "")),]$p[aseg_id]

brain_dic[which(p.adjust(c(p, p.scs[nuclei_id]), method = "fdr") < 0.0025)]

data = data.frame(roi = aseg$data$region[3:26],
                  t = t.scs*(p.adjust(p.scs,"fdr")<0.05),stringsAsFactors = FALSE,
  hemi = aseg$data$hemi[3:26],label = aseg$data$label[3:26],Modal = rep("Volume", 24)) # *(p.adjust(p,method = "fdr")<0.05)
ggseg(data, atlas = aseg, colour = "grey",mapping = aes(fill = t)) +  theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1],na.value = "lightgrey",
    breaks = seq(-4.5, 4.5, 0.9),limits = c(-4.5, 4.5)) +  theme_brain()#theme_void()

data = data.frame(roi = rep(dk_id, 2),t = t*(p.adjust(p,"fdr")<0.05),stringsAsFactors = FALSE,
  Modal = c(rep(con_cds, 68)),hemi = rep(c(rep("left", 34), rep("right", 34)), 1)) # *(p.adjust(p,method = "fdr")<0.05)
data %>%  group_by(Modal) %>% ggseg(atlas = dk, colour = "grey", 
                                    mapping = aes(fill = t),position = "stack") +
  facet_wrap( ~ Modal, ncol = 1) + theme(legend.position = "right") +
  scale_fill_gradientn(colours = brewer.pal(11, "RdBu")[11:1], na.value = "lightgrey",
    breaks = seq(-4.5, 4.5, 0.9),limits = c(-4.5, 4.5)) +  theme_brain()


# nBack
nback_dic = c('sacgvf','sacsvcg','saisvcg','saasvcg','sacsvis','saigvcg','saigvis')
nback_dic2 = c('nBack.0.back.condition','nBack.2.back.condition','nBack.place.condition',
  'nBack.emotion.condition','nBack.2.back.versus.0.back.contrast',
  'nBack.face.versus.place.contrast','nBack.emotion.versus.neutral.face.contrast',
  'nBack.negative.face.versus.neutral.face.contrast','nBack.positive.face.versus.neutral.face.contrast')
abcd_Nback.label = nback_bwroi02.label[c(4, 9, 16, 18, 20, 22:903)]

nuclei_id_nback <- c(6, 7, 8, 9, 13, 14, 16, 23:29)
con_id <- c(5, 8, 9);surf_id <- c();scs_id <- c()

for (i in 1:length(con_id)) {
  surf_id <- c(surf_id, which(str_detect(
      abcd_Nback.label, paste(nback_dic2[con_id[i]], '.in.APARC', sep = ""))))
  scs_id <- c(scs_id, which(str_detect(
      abcd_Nback.label, paste(nback_dic2[con_id[i]], '.in.ASEG', sep = ""))))
}

abcd_BrainVar <- merge(abcd_covariables[, c(1:5, 14:23)], abcd_envir.bl[, c(1, 8:15, 27:28, 80)], all = T)
abcd_BrainVar <- merge(abcd_BrainVar, abcd_nBack[, c(1:5, surf_id, scs_id)])
abcd_BrainVar$eventname <- factor(abcd_BrainVar$eventname,
  levels = c("baseline_year_1_arm_1","2_year_follow_up_y_arm_1"),ordered = T)
var_dic <-  c("eventname.L","clust1","clust2","clust3","clust4","clust5",
              "clust1:eventname.L","clust2:eventname.L","clust3:eventname.L",
              "clust4:eventname.L","clust5:eventname.L")
n <- length(var_dic)
fit.parameter.brain <- data.frame()
brain_dic <- names(abcd_nBack)[c(surf_id, scs_id)]
for (i in 1:length(brain_dic)) {
  f <- as.formula(paste(brain_dic[i],
        "~1+clust+eventname+clust*eventname+sex+income+edu+",
        # nihtbx_totalcomp_uncorrected+
        "race_6level+hisp+Pubertal+anthro_bmi_calc+tfmri_nback_all_beta_mm+",
        "(1|mri_info_deviceserialnumber)+(1|subjectkey)",sep = ""))
  fit.lme.brain <- lmerTest::lmer(f, abcd_BrainVar)
  fit.lme.brain.summary <- summary(fit.lme.brain)
  fit.parameter.brain[((i - 1) * n + 1):(i * n), "var"] <- rep(brain_dic[i], n)
  fit.parameter.brain[((i - 1) * n + 1):(i * n), "contrast"] <- var_dic
  fit.parameter.brain[((i - 1) * n + 1):(i * n), c("t", "df", "p")] <-
    fit.lme.brain.summary$coefficients[var_dic, c("t value", "df", "Pr(>|t|)")]
  print(paste(brain_dic[i], " linear mixed effect is done. no. ", i, sep = ""))
}

fit.parameter.brain[fit.parameter.brain$p < .0025,]
fit.nback <- fit.parameter.brain
fit.nback.dcast.t <- dcast(data = fit.nback,var ~ contrast,value.var = "t")
fit.nback.dcast.t$item <- rep("t",nrow(fit.nback.dcast.t))
fit.nback.dcast.df <- dcast(data = fit.nback,var ~ contrast,value.var = "df")
fit.nback.dcast.df$item <- rep("df",nrow(fit.nback.dcast.df))
fit.nback.dcast.p <- dcast(data = fit.nback,var ~ contrast,value.var = "p")
fit.nback.dcast.p$item <- rep("p",nrow(fit.nback.dcast.p))
fit.nback.dcast <- rbind(fit.nback.dcast.t,fit.nback.dcast.df,fit.nback.dcast.p)
fit.nback.dcast <- fit.nback.dcast[,c("var","item",var_dic)]


write.table(fit.nback,file = "fit_nback.csv", na = "",sep = ",",row.names = FALSE)

dk_id <- formatC(c(2:4, 6:36), width = 4, flag = '0')
aseg_id <-  c(6,23,19,13,3,25,28,8,30,17,27,14,9,26,7,24,10,11,12,35,33,34,32,31)

var = "clust3"
i = 2
t = c(fit.nback[fit.nback$contrast == var &
                  str_detect(abcd_Nback.label[c(surf_id, scs_id)],
                             paste(nback_dic2[con_id[i]], '.in.APARC', sep = "")),][1:68,]$t)
p = c(fit.nback[fit.nback$contrast == var &
                  str_detect(abcd_Nback.label[c(surf_id, scs_id)],
                             paste(nback_dic2[con_id[i]], '.in.APARC', sep = "")),][1:68,]$p)
t.scs = fit.nback[fit.nback$contrast == var &
                    str_detect(abcd_Nback.label[c(surf_id, scs_id)],
                               paste(nback_dic2[con_id[i]], '.in.ASEG', sep = "")),][1:68,]$t[aseg_id]
p.scs = fit.nback[fit.nback$contrast == var &
                    str_detect(abcd_Nback.label[c(surf_id, scs_id)],
                               paste(nback_dic2[con_id[i]], '.in.ASEG', sep = "")),][1:68,]$p[aseg_id]

brain_dic[which(p.adjust(c(p, p.scs[nuclei_id_nback]), method = "fdr") <
                  0.025)]

data = data.frame(
  roi = aseg$data$region[3:26],
  t = t.scs,
  stringsAsFactors = FALSE,
  hemi = aseg$data$hemi[3:26],
  label = aseg$data$label[3:26],
  Modal = rep("Volume", 24)
) # *(p.adjust(p,method = "fdr")<0.05)
ggseg(data,
      atlas = aseg,
      colour = "grey",
      mapping = aes(fill = t)) +
  theme(legend.position = "right") +
  scale_fill_gradientn(
    colours = brewer.pal(11, "RdBu")[11:1],
    na.value = "lightgrey",
    breaks = seq(-4.5, 4.5, 0.9),
    limits = c(-4.5, 4.5)
  ) +
  theme_brain()#theme_void()

data = data.frame(
  roi = rep(dk_id, 2),
  t = t,
  stringsAsFactors = FALSE,
  Modal = c(rep(nback_dic2[con_id[i]], 68)),
  hemi = rep(c(rep("left", 34), rep("right", 34)), 1)
) # *(p.adjust(p,method = "fdr")<0.05)
data %>%
  group_by(Modal) %>% ggseg(atlas = dk,
                            colour = "grey",
                            mapping = aes(fill = t)) +
  facet_wrap( ~ Modal, ncol = 1) + theme(legend.position = "right") +
  scale_fill_gradientn(
    colours = brewer.pal(11, "RdBu")[11:1],
    na.value = "lightgrey",
    breaks = seq(-4.5, 4.5, 0.9),
    limits = c(-4.5, 4.5)
  ) +
  theme_brain()


# SST
con_id <- c(2);surf_id <- c();scs_id <- c()
sst_dic = c('sacgvf','sacsvcg','saisvcg','saasvcg','sacsvis','saigvcg','saigvis')
nuclei_id_sst <- c(6, 7, 8, 9, 13, 14, 16, 23:29)

# SST correct go versus fixation contrast
# SST correct stop versus correct go contrast
# SST incorrect stop versus correct go contrast
# SST any stop versus correct go contrast
# SST correct stop versus incorrect stop contrast
# SST incorrect go versus correct go contrast
# SST incorrect go versus incorrect stop contrast


for (i in 1:length(con_id)) {
  surf_id <-
    c(surf_id, which(str_detect(
      names(abcd_SST), paste(sst_dic[con_id[i]], '_bcdk_', sep = "")
    )))
  scs_id <-
    c(scs_id, which(str_detect(
      names(abcd_SST), paste(sst_dic[con_id[i]], '_bscs_', sep = "")
    )))
}

abcd_BrainVar <- merge(abcd_covariables[, c(1:5, 14:23)], abcd_envir.bl[, c(1, 8:15, 27:28, 80)], all = T)
abcd_BrainVar <- merge(abcd_BrainVar, abcd_SST[, c(1:5, surf_id, scs_id)])
abcd_BrainVar$eventname <- factor(abcd_BrainVar$eventname,
  levels = c("baseline_year_1_arm_1", "2_year_follow_up_y_arm_1"),ordered = T)
var_dic <-  c("eventname.L","clust1","clust2","clust3","clust4","clust5",
              "clust1:eventname.L","clust2:eventname.L","clust3:eventname.L",
              "clust4:eventname.L","clust5:eventname.L")
n <- length(var_dic)
fit.parameter.brain <- data.frame()
brain_dic <- names(abcd_SST)[c(surf_id, scs_id)]
for (i in 1:length(brain_dic)) {
  f <- as.formula(paste(brain_dic[i],
        "~1+clust+eventname+clust*eventname+sex+income+edu+",
        "race_6level+hisp+Pubertal+anthro_bmi_calc+tfmri_sa_beta_mm+",
        "(1|mri_info_deviceserialnumber)+(1|subjectkey)",sep = ""))
  fit.lme.brain <- lmerTest::lmer(f, abcd_BrainVar)
  fit.lme.brain.summary <- summary(fit.lme.brain)
  fit.parameter.brain[((i - 1) * n + 1):(i * n), "var"] <- rep(brain_dic[i], n)
  fit.parameter.brain[((i - 1) * n + 1):(i * n), "contrast"] <- var_dic
  fit.parameter.brain[((i - 1) * n + 1):(i * n), c("t", "df", "p")] <-
    fit.lme.brain.summary$coefficients[var_dic, c("t value", "df", "Pr(>|t|)")]
  print(paste(brain_dic[i], " linear mixed effect is done. no. ", i, sep = ""))
}

fit.parameter.brain[fit.parameter.brain$p < .00025,]
fit.SST <- fit.parameter.brain
write.table(
  fit.parameter.brain,
  file = "abcd_suicidal_SST.csv",
  na = "",
  sep = ",",
  row.names = FALSE
)

dk_id <- formatC(c(2:4, 6:36), width = 4, flag = '0')
aseg_id <-  c(6,23,19,13,3,25,28,8,30,17,27,14,9,26,7,24,10,11,12,35,33,34,32,31)

var = "clust5"
con_cds = mid_dic2[con_id[3]]
con_scs = mid_dic[con_id[3]]
t = c(fit.SST[fit.SST$contrast == var &
                str_detect(fit.SST$var, paste(con_cds, "_b_cds_", sep = "")),][1:68,]$t)
p = c(fit.SST[fit.SST$contrast == var &
                str_detect(fit.SST$var, paste(con_cds, "_b_cds_", sep = "")),][1:68,]$p)
t.scs = fit.SST[fit.SST$contrast == var &
                  str_detect(fit.SST$var, paste(con_scs, "_b_scs_", sep = "")),]$t[aseg_id]
p.scs = fit.SST[fit.SST$contrast == var &
                  str_detect(fit.SST$var, paste(con_scs, "_b_scs_", sep = "")),]$p[aseg_id]
brain_dic[which(p.adjust(c(p, p.scs[nuclei_id_sst]), method = "fdr") < 0.0025)]

data = data.frame(
  roi = aseg$data$region[3:26],
  t = t.scs,
  stringsAsFactors = FALSE,
  hemi = aseg$data$hemi[3:26],
  label = aseg$data$label[3:26],
  Modal = rep("Volume", 24)
) # *(p.adjust(p,method = "fdr")<0.05)
ggseg(data,
      atlas = aseg,
      colour = "grey",
      mapping = aes(fill = t)) +
  theme(legend.position = "right") +
  scale_fill_gradientn(
    colours = brewer.pal(11, "RdBu")[11:1],
    na.value = "lightgrey",
    breaks = seq(-4.5, 4.5, 0.9),
    limits = c(-4.5, 4.5)
  ) +
  theme_brain()#theme_void()

data = data.frame(
  roi = rep(dk_id, 2),
  t = t,
  stringsAsFactors = FALSE,
  Modal = c(rep(nback_dic2[con_id[i]], 68)),
  hemi = rep(c(rep("left", 34), rep("right", 34)), 1)
) # *(p.adjust(p,method = "fdr")<0.05)
data %>%
  group_by(Modal) %>% ggseg(atlas = dk,
                            colour = "grey",
                            mapping = aes(fill = t)) +
  facet_wrap( ~ Modal, ncol = 1) + theme(legend.position = "right") +
  scale_fill_gradientn(
    colours = brewer.pal(11, "RdBu")[11:1],
    na.value = "lightgrey",
    breaks = seq(-4.5, 4.5, 0.9),
    limits = c(-4.5, 4.5)
  ) +
  theme_brain()
