# Run after 01_data_and_clustering.R in the same R session.
# Group comparisons: original code uses saved abcd_envir.bl and positional column indices.
# Verify the variable order, covariates, multiple-testing family and output tables.
if (!exists("abcd_envir.bl")) stop("Run 01_data_and_clustering.R first")
if (!requireNamespace("chisq.posthoc.test", quietly=TRUE)) stop("Install chisq.posthoc.test")
library(chisq.posthoc.test)
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
dir.create("results",showWarnings=FALSE)
write.csv(fit.parameter,"results/behavior_anova.csv",row.names=FALSE)
write.csv(fit.HSD,"results/behavior_pairwise_tukey.csv",row.names=FALSE)
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
# PRS: original executed code used prs_prsice. Its ADHD ANOVA reproduces Table S10.
# The manuscript describes PRS-CS; reconcile this provenance before publication.
# This is a one-way ANOVA; it does not adjust for genetic PCs, batch or diagnoses.
if (!requireNamespace("reshape2",quietly=TRUE)) stop("Install reshape2")
library(reshape2)
# abcd_prscs.clust<-merge(abcd_envir.bl[,c("subjectkey","clust")],prs_prscs)
abcd_prscs.clust <-
  merge(abcd_envir.bl[, c("subjectkey", "clust")], prs_prsice)

var_id <- match(c("adhd", "bip", "MDD", "OCD", "PTSD", "scz"), names(abcd_prscs.clust)); stopifnot(!anyNA(var_id))
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
dir.create("results",showWarnings=FALSE)
write.csv(fit.HSD, file.path("results","prs_pairwise.csv"), row.names=FALSE)
write.csv(fit.parameter, file.path("results","prs_anova.csv"), row.names=FALSE)
