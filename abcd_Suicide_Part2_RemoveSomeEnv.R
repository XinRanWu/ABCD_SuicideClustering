########################################################
#                 Function Definition
########################################################
pacman::p_load(cluster,factoextra,psych,ggradar,pheatmap,ggplot2,
               ggiraph,ggiraphExtra,RColorBrewer,tableone,glmnet,
               rgl,ggseg,dplyr,reshape2,Rmisc,ggpubr,export,Rcpp,
               ggridges,stringr,ggfun,gghalves,chisq.posthoc.test)

# devtools::install_github("ebbertd/chisq.posthoc.test")

########################################################
#                  Behavioral Profile
########################################################
abcd_trait.bl.clust$clust <- abcd_suicide.clust$cluster

x <- abcd_trait.bl
x <- merge(x, abcd_covariables_bl_envir[, c(1, 43, 44)], by = "subjectkey")
x$clust <- NA
x$clust[x$ksad01_nodiagnosis == 1 & x$ksad501_nodiagnosis == 1] <- 0
rownames(x) <- x$subjectkey
x[row.names(abcd_trait.bl.clust), "clust"] <-
  abcd_trait.bl.clust$clust
x$clust <- as.factor(x$clust)
abcd_envir.bl <- merge(abcd_covariables_bl_envir, x)
rm(x)

abcd_ksads_suicide.bl <-
  abcd_ksads_suicide[abcd_ksads_suicide$eventname == "baseline_year_1_arm_1", ]
abcd_stq01.bl <-
  abcd_stq01[abcd_stq01$eventname == "baseline_year_1_arm_1", ]
# names(abcd_ksad501_suicide.bl)[3:51]<-abcd_ksad501.label[which(stringr::str_detect(names(abcd_ksad501),"ksads_23_"))]
abcd_envir.bl <- merge(abcd_envir.bl, abcd_ksads_suicide.bl)
abcd_envir.bl <- merge(abcd_envir.bl, dsm_14disorder.bl)
abcd_envir.bl <- merge(abcd_envir.bl, abcd_stq01.bl)

profile.name <-
  c("picv","flanker","list","cardsort","pattern","pic","reading","anxdep",
    "withdep","somatic","social","tho","attent","ruleb","aggre","depress",
    "anxdis","somaticpr","adhd","opposit","conduct","sct","ocd","stress",
    "pps","neg urge","lack plan","sen seek","pos urge","lack pers","bis",
    "bas rr","bas drive" ,"bas fs")

abcd_envir.bl.longdat <-
  abcd_envir.bl[-1 * which(is.na(abcd_envir.bl$clust)), c(45:78, 80)]
names(abcd_envir.bl.longdat) <- c(profile.name, "Cluster")
abcd_envir.bl.longdat$Cluster <-
  paste("Subtype ", abcd_envir.bl.longdat$Cluster, sep = "")
abcd_envir.bl.longdat$Cluster[abcd_envir.bl.longdat$Cluster == "Subtype 0"] <-
  "HC"
abcd_envir.bl.longdat <-
  abcd_envir.bl.longdat[, c(35, 1:7, 8:11, 16:18, 24, 12:15, 19:23, 26:30, 32:34, 31, 25)]


# Radar Plot
angle <-
  c(90 - seq(
    0 + 360 / ncol(abcd_envir.bl.longdat),
    180,
    360 / ncol(abcd_envir.bl.longdat)
  ),
  90 - seq(
    0 + 360 / ncol(abcd_envir.bl.longdat),
    180,
    360 / ncol(abcd_envir.bl.longdat)
  ))
ggRadar(
  data = abcd_envir.bl.longdat,
  aes(
    fill = Cluster,
    colour = Cluster,
    facet = Cluster
  ),
  size = .8,
  alpha = .8
) +
  theme_minimal() +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  theme(axis.text.x = element_text(angle = angle, hjust = 1))
graph2ppt(
  file = "subtypes.profiles.pptx",
  width = 9,
  aspectr = sqrt(2),
  append = TRUE
)

ggRadar(
  data = abcd_envir.bl.longdat,
  aes(fill = Cluster, colour = Cluster),
  size = .1,
  alpha = .5
) + theme_minimal() +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  theme(text = element_text(size = 16, family = "sans"))
ggDonut(abcd_envir.bl.longdat, aes(donuts = Cluster, count = n))

i = 1
ggRadar(
  data = abcd_envir.bl.longdat[c(
    which(abcd_envir.bl.longdat$Cluster == "HC"),
    which(abcd_envir.bl.longdat$Cluster == paste("Subtype", i, sep = " "))
  ), ],
  aes(fill = Cluster, colour = Cluster),
  size = .8,
  alpha = .4,
  ylim = c(0, .8),
  use.label = FALSE,
  legend.position = "none"
) +
  theme_minimal() +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1")))

# Prepare for Multiple Bar Plot (Err Bar)
behav.longdat = melt(
  data = abcd_envir.bl[, c(1, 45:78)],
  id.vars = c("subjectkey"),
  variable.name = "Feature",
  value.name = "Value"
)
behav.longdat <- merge(behav.longdat, abcd_envir.bl[, c(1, 80)])
behav.longdat <- na.omit(behav.longdat)
behav.summary <-
  summarySE(behav.longdat,
            measurevar = "Value",
            groupvars = c("clust", "Feature"))
behav.summary$Feature <- profile.name
behav.summary$Questionnaire <-
  c(rep(c(
    rep("NIH", 7),
    rep("CBCL", 17),
    "PPS",
    rep("UPPS", 5),
    rep("BIS/BAS", 4)
  ), 6))

x <- abcd_envir.bl[, c(45:78)]
for (i in 1:34) {
  x[, i] <-
    (x[, i] - min(x[, i], na.rm = T)) / (max(x[, i], na.rm = T) - min(x[, i], na.rm = T))
}

behav.profile <-
  aggregate(x, by = list(abcd_envir.bl$clust), mean, na.rm = TRUE)
behav.profile <- behav.profile[, -1]
# names(behav.profile)<-c("pv","fl","li","cs","pa","pic","re",
#                         "ad","wd","sm","soc","tho","att","rb",
#                         "ag","dep" ,"anx","smp","adhd","odd","cd",
#                         "sct","ocd","str","pps","nu","lpl","ss" ,
#                         "pu","lpe","bis","rr","drive" ,"fs")
names(behav.profile) <- profile.name
behav.profile$Cluster <-
  c("HC",
    "Subtype 1",
    "Subtype 2",
    "Subtype 3",
    "Subtype 4",
    "Subtype 5")
behav.profile.longdat = melt(
  data = behav.profile,
  id.vars = c("Cluster"),
  variable.name = "Feature",
  value.name = "Value"
)
behav.profile.longdat$Scale <- c(
  rep("NIH", 7 * 6),
  rep("CBCL", 17 * 6),
  rep("PPS", 6),
  rep("UPPS", 5 * 6),
  rep("BIS/BAS", 4 * 6)
)
angle <-
  c(90 - seq(
    0 + 360 / ncol(abcd_envir.bl.longdat),
    180,
    360 / ncol(abcd_envir.bl.longdat)
  ),
  90 - seq(
    0 + 360 / ncol(abcd_envir.bl.longdat),
    180,
    360 / ncol(abcd_envir.bl.longdat)
  ))
behav.profile.longdat$angle <-
  as.numeric(t(matrix(rep(angle, 6), 34, 6)))

# Multiple Bar Plot (no Err Bar)
ggplot(behav.profile.longdat,
       aes(x = Feature, y = Value, fill = Cluster),
       alpha = 0.5) +
  geom_bar(stat = "identity", position = "dodge") + ylim(-.01, .9) +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  facet_grid( ~ Scale, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    size = 12
  ))

ggplot(
  behav.profile.longdat,
  aes(
    x = Feature,
    y = Value,
    colour = Cluster,
    group = Cluster
  ),
  alpha = 1
) +
  geom_point(size = 2, shape = 16) + geom_line() + ylim(-.01, .9) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  facet_grid( ~ Scale, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    size = 12
  )) # shape=Cluster


graph2ppt(
  file = "subtypes.profiles.bar.pptx",
  width = 9,
  aspectr = sqrt(2),
  append = TRUE
)


# Prepare for Multiple Bar Plot (Err Bar) for Envir Var
envir.name <-
  c("income","edu","neighbor safety","family conflict p","family conflict y",
    "acceptance parent","acceptance caregiver","parental monitor",
    "school environment","school involvement","school disengagement ",
    "pubertal","sleep disturb","physic exercise",
    "screen mature games","screen r-movie","screen tv" ,"screen video","screen games",
    "screen text chat","screen social", "screen video chat"  )
envir.longdat = melt(
  data = abcd_envir.bl[, c(1, 27:29, 30, 38, 40:41, 37,  42:44, 19, 32, 33, 193:200)],
  id.vars = c("subjectkey"),
  variable.name = "Feature",
  value.name = "Value"
)
envir.longdat <- merge(envir.longdat, abcd_envir.bl[, c("subjectkey","clust")])
envir.longdat <- na.omit(envir.longdat)
envir.summary <-
  summarySE(envir.longdat,
            measurevar = "Value",
            groupvars = c("clust", "Feature")
            )
envir.summary$Feature <- rep(envir.name, 6)
envir.summary$Questionnaire <-
  rep(c(
    rep("SES", 3),
    rep("Family", 5),
    rep("School", 3),
    rep("Health", 11)
  ), 6)

x <-
  abcd_envir.bl[, c(27:29, 30, 38, 40:41, 37, 42:44, 19, 32, 33, 193:200)]
for (i in 1:22) {
  x[, i] <-
    (x[, i] - min(x[, i], na.rm = T)) / (max(x[, i], na.rm = T) - min(x[, i], na.rm = T))
}

envir.profile <-
  aggregate(x, by = list(abcd_envir.bl$clust), mean, na.rm = TRUE)
envir.profile <- envir.profile[, -1]
names(envir.profile) <- envir.name
envir.profile$Cluster <-
  c("HC",
    "Subtype 1",
    "Subtype 2",
    "Subtype 3",
    "Subtype 4",
    "Subtype 5")
envir.profile.longdat = melt(
  data = envir.profile,
  id.vars = c("Cluster"),
  variable.name = "Feature",
  value.name = "Value"
)
envir.profile.longdat$Scale <- c(
  rep("SES", 3 * 6),
  rep("Family", 5 * 6),
  rep("School", 3 * 6),
  rep("Health", 3 * 6),
  rep("Screen", 8 * 6)
)

# Multiple Bar Plot (no Err Bar)
ggplot(envir.profile.longdat,
       aes(x = Feature, y = Value, fill = Cluster),
       alpha = 0.5) +
  geom_bar(stat = "identity", position = "dodge") + ylim(-.01, .95) +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  facet_grid( ~ Scale, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    size = 12
  ))
# graph2ppt(file="subtypes.profiles.bar.pptx", width=9, aspectr=sqrt(2), append=TRUE)

# Multiple Line Plot (no Err Bar)
ggplot(
  envir.profile.longdat,
  aes(
    x = Feature,
    y = Value,
    colour = Cluster,
    group = Cluster
  ),
  alpha = 1
) +
  geom_point(size = 2, shape = 16) + geom_line() + ylim(-.01, .95) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  facet_grid( ~ Scale, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    size = 12
  ))



# Multiple Circle Bar Plot (no Err Bar)
grid_data <-
  data.frame(
    start = c(.2, 7.2, 24.2, 25.2, 30.2),
    end = c(6.8, 23.8, 24.8, 29.8, 33.8)
  )
base_data <-
  data.frame(
    group = c("NIH", "CBCL", "PPS", "UPPS", "BIS/BAS"),
    title = rowMeans(grid_data)
  )
g1 <-
  ggplot(behav.profile.longdat,
         aes(x = Feature, y = Value, fill = Cluster),
         alpha = 0.5) +
  geom_bar(stat = "identity") + ylim(-.5, 1) +
  facet_wrap(. ~ Cluster) + coord_polar() + theme_void() +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1")))
g2 <- g1 + geom_segment(
  data = grid_data,
  aes(
    x = start,
    y = -.1,
    xend = end,
    yend = -.1
  ),
  colour = "black",
  alpha = 1,
  size = .5,
  inherit.aes = FALSE
)
g3 <- g2 + geom_text(
  data = behav.profile.longdat,
  aes(
    x = Feature,
    y = Value + .1,
    label = Feature,
    angle = angle
  ),
  color = "black",
  fontface = "bold",
  alpha = 0.6,
  size = 2.5
)
g4 <-
  g3 + geom_text(
    data = base_data,
    aes(x = title, y = -.25, label = group),
    colour = "black",
    alpha = 0.8,
    size = 2,
    fontface = "bold",
    inherit.aes = FALSE
  )
g4
# theme(axis.text.x = element_text(angle = angle, hjust = 1))


# Multiple Bar Plot (Err Bar)
ggplot(behav.summary, aes(x = Feature, y = Value, fill = clust)) +
  geom_bar(stat = "identity",  position = position_dodge()) +
  geom_errorbar(
    aes(ymin = Value - sd, ymax = Value + sd),
    width = 0.05,
    colour = "gray",
    position = position_dodge(0.9)
  ) +
  facet_grid(. ~ Questionnaire, scales = "free", space = "free") +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) # scale_x_discrete(limits=c("2", "0.5", "1"))

# Multiple Point Plot (Err Bar)
ggplot(behav.summary, aes(x = Feature, y = Value, color = clust)) +
  geom_line(position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = Value - sd, ymax = Value + sd),
                width = 0.2,
                position = position_dodge(0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(. ~ Questionnaire, scales = "free", space = "free") +
  scale_fill_brewer(palette = "Spectral", type = "div")


# Multiple Bar Plot of KSADS Diagnosis
abcd_ksads.suicide.label <- read.csv2("ksads_suicidal_labels2.csv")
behav.longdat = melt(
  data = abcd_envir.bl[, c(1, 25, 27:44, 193:200)],
  id.vars = c("subjectkey"),
  variable.name = "Feature",
  value.name = "Value"
)

x <- abcd_envir.bl[, c(1, 81:192)]
# names(x)[2:99]<-abcd_ksads.suicide.label$Names
ksads.longdat.bl = melt(
  data = x,
  id.vars = c("subjectkey"),
  variable.name = "Diagnosis",
  value.name = "Value"
)
ksads.longdat.bl <- na.omit(ksads.longdat.bl)
ksads.longdat.bl <-
  merge(ksads.longdat.bl, abcd_envir.bl[, c(1, 80)], by = "subjectkey", all = T)
ksads.longdat.bl <- ksads.longdat.bl[-1 * which(ksads.longdat.bl$Value == 0), ]
ksads.longdat.bl <-
  merge(ksads.longdat.bl,
        abcd_ksads.suicide.label,
        by.x = "Diagnosis",
        by.y = "Code")
rm(x)

ksads.longdat.Combid.bl <-
  ksads.longdat.bl[which(ksads.longdat.bl$Questionnaire == "Comorbidity"), ]



ggplot(ksads.longdat.Combid.bl, aes(Names, fill = clust)) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values = c(brewer.pal(5, "Set1"),"#B0B0B0")) +
  scale_color_manual(values = c(brewer.pal(5, "Set1"),"#B0B0B0")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(.~Class, scales = "free", space = "free") +
  ylim(0.6,1)

ggplot(ksads.longdat.Combid.bl, aes(Names, fill = clust)) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values = c(brewer.pal(5, "Set1"),"#B0B0B0")) +
  scale_color_manual(values = c(brewer.pal(5, "Set1"),"#B0B0B0")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(Class~., scales = "free", space = "free") + 
  ylim(0.6,1) + 
  coord_flip()

ksads.longdat.noCombid.bl <-
  ksads.longdat.bl[-1 * which(ksads.longdat.bl$Questionnaire == "Comorbidity"), ]
ksads.longdat.noCombid.bl[which(ksads.longdat.noCombid.bl$Class == "Self Injurious Behavior"), "Class"] <-
  "Self Injury"
ggplot(ksads.longdat.noCombid.bl, aes(Names, fill = clust)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(Questionnaire ~ Class, scales = "free", space = "free")

ggplot(data = ksads.longdat.noCombid.bl[
  -1*c(which(is.na(ksads.longdat.noCombid.bl$clust)),
       which(ksads.longdat.noCombid.bl$clust=="0")), ],
  aes(y = Questionnaire,fill = clust)) + 
  geom_bar(position = "fill",width = 0.9) + 
  coord_polar() + 
  scale_fill_manual(values = c(brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c(brewer.pal(5, "Set1"))) +
  theme_void() + scale_y_discrete(limits = c("NA","K-SADS-P","K-SADS-Y"))

tableone::CreateTableOne()

g1 <-
  ggplot(abcd_envir.bl[which(!is.na(abcd_envir.bl$clust) &
                               abcd_envir.bl$clust != "0"), ], aes(sex, fill = clust)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c(brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c(brewer.pal(5, "Set1"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()

g2 <-
  ggplot(abcd_envir.bl[which(!is.na(abcd_envir.bl$clust) &
                               abcd_envir.bl$clust != "0"), ], aes(race_6level, fill = clust)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c(brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c(brewer.pal(5, "Set1"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()

require(gridExtra)
grid.arrange(g1,g2,nrow=2,ncol=1,heights = c(0.9,0.9),widths = c(1,1))

########################################################
#               Longitudinal Analysis
########################################################

abcd_trait.fu1 <-
  abcd_trait[abcd_trait$eventname == "1_year_follow_up_y_arm_1", ]
abcd_trait.fu2 <-
  abcd_trait[abcd_trait$eventname == "2_year_follow_up_y_arm_1", ]
dsm_14disorder.fu2 <-
  dsm_14disorder[dsm_14disorder$eventname == "2_year_follow_up_y_arm_1", ]
dsm_14disorder.fu1 <-
  dsm_14disorder[dsm_14disorder$eventname == "1_year_follow_up_y_arm_1", ]

abcd_trait.bl.clust$clust <- abcd_suicide.clust$cluster

rownames(abcd_trait.fu2) <- abcd_trait.fu2$subjectkey

x <- abcd_trait.fu2
x <- merge(x, abcd_covariables_bl_envir[, c(1, 43, 44)])
x$clust <- NA
x$clust[x$ksad01_nodiagnosis == 1 & x$ksad501_nodiagnosis == 1] <- 0
rownames(x) <- x$subjectkey
x[row.names(abcd_trait.bl.clust), "clust"] <-
  abcd_trait.bl.clust$clust
x$clust <- as.factor(x$clust)
x <- x[, -2]

abcd_envir.fu2 <-
  merge(abcd_covariables_bl_envir[, c(1:26)], x, by = c("subjectkey"))
rm(x)

# abcd_envir.fu2<-
abcd_envir.fu2 <-
  merge(abcd_envir.fu2, abcd_ksads_suicide.fu2[, c(1, 3:100)], by = c("subjectkey"))
abcd_envir.fu2 <-
  merge(abcd_envir.fu2, dsm_14disorder.fu2[, c(1, 5:18)], by = c("subjectkey"))
abcd_envir.fu2[, c(4, 14, 15, 24)] <-
  lapply(abcd_envir.fu2[, c(4, 14, 15, 24)], factor, exclude =  "")
abcd_envir.fu2[, c(61, 62, 64:175)] <-
  lapply(abcd_envir.fu2[, c(61, 62, 64:175)],
         factor,
         ordered = TRUE,
         levels = c(0, 1))

tabone <-
  CreateTableOne(data = abcd_envir.fu2[, c(3, 4, 27:60, 64:175)])
print(tabone)

xt<-xtabs(formula =  ~ clust+Class+Questionnaire,
          data = ksads.longdat.noCombid[
            -1*c(which(is.na(ksads.longdat.noCombid$clust)), 
                 which(ksads.longdat.noCombid$clust=="0")),])

abcd_envir.longitudinal <-
  merge(abcd_envir.bl[, c(1, 3:80, 179:200)], abcd_envir.fu2[, c(1, 64:161)])
# step(fit.glm.longitudinal)


abcd_envir.fu2.longdat <-
  abcd_envir.fu2[-1 * which(is.na(abcd_envir.fu2$clust)), c(27:60, 63)]
names(abcd_envir.fu2.longdat) <- c(profile.name, "Cluster")
abcd_envir.fu2.longdat$Cluster <-
  paste("Subtype ", abcd_envir.fu2.longdat$Cluster, sep = "")
abcd_envir.fu2.longdat$Cluster[abcd_envir.fu2.longdat$Cluster == "Subtype 0"] <-
  "HC"
abcd_envir.fu2.longdat <-
  abcd_envir.fu2.longdat[, c(35, 1:7, 8:11, 16:18, 24, 12:15, 19:23, 26:30, 32:34, 31, 25)]
abcd_envir.fu2.longdat <- abcd_envir.fu2.longdat[, c(-4, -5)]

# Radar Plot
angle <-
  c(90 - seq(
    0 + 360 / ncol(abcd_envir.fu2.longdat),
    180,
    360 / ncol(abcd_envir.fu2.longdat)
  ),
  90 - seq(
    0 + 360 / ncol(abcd_envir.fu2.longdat),
    180,
    360 / ncol(abcd_envir.fu2.longdat)
  ))
ggRadar(
  data = abcd_envir.fu2.longdat,
  aes(
    fill = Cluster,
    colour = Cluster,
    facet = Cluster
  ),
  size = .8,
  alpha = .8
) +
  theme_minimal() +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  theme(axis.text.x = element_text(angle = angle, hjust = 1))

x <- abcd_envir.fu2[, c(27:60)]
for (i in 1:34) {
  x[, i] <-
    (x[, i] - min(x[, i], na.rm = T)) / (max(x[, i], na.rm = T) - min(x[, i], na.rm = T))
}

behav.profile <-
  aggregate(x,
            by = list(abcd_envir.fu2$clust),
            mean,
            na.rm = TRUE)
behav.profile <- behav.profile[, c(-1, -4, -5)]
names(behav.profile) <- c("pv","fl","pa","pic","re","ad","wd","sm","soc",
                          "tho","att","rb","ag","dep" ,"anx","smp","adhd",
                          "odd","cd","sct","ocd","str","pps","nu","lpl",
                          "ss" ,"pu","lpe","bis","rr","drive","fs")

names(behav.profile) <- profile.name[c(1:2, 5:35)]
behav.profile$Cluster <-
  c("HC",
    "Subtype 1",
    "Subtype 2",
    "Subtype 3",
    "Subtype 4",
    "Subtype 5")
behav.profile.longdat = melt(
  data = behav.profile,
  id.vars = c("Cluster"),
  variable.name = "Feature",
  value.name = "Value"
)
behav.profile.longdat$Scale <- c(
  rep("NIH", 5 * 6),
  rep("CBCL", 17 * 6),
  rep("PPS", 6),
  rep("UPPS", 5 * 6),
  rep("BIS/BAS", 4 * 6)
)
angle <-
  c(90 - seq(
    0 + 360 / ncol(abcd_envir.fu2.longdat),
    180,
    360 / ncol(abcd_envir.fu2.longdat)
  ),
  90 - seq(
    0 + 360 / ncol(abcd_envir.fu2.longdat),
    180,
    360 / ncol(abcd_envir.fu2.longdat)
  ))
behav.profile.longdat$angle <-
  as.numeric(t(matrix(rep(angle, 6), 32, 6)))

# Multiple Bar Plot (no Err Bar)
ggplot(behav.profile.longdat,
       aes(x = Feature, y = Value, fill = Cluster),
       alpha = 0.5) +
  geom_bar(stat = "identity", position = "dodge") + ylim(-.01, .9) +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  facet_grid( ~ Scale, scales = "free", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    size = 12
  ))
# graph2ppt(file="subtypes.profiles.bar.pptx", width=9, aspectr=sqrt(2), append=TRUE)


# Multiple Circle Bar Plot (no Err Bar)
grid_data <-
  data.frame(
    start = c(.2, 5.2, 22.2, 23.2, 28.2),
    end = c(4.8, 21.8, 22.8, 27.8, 31.8)
  )
base_data <-
  data.frame(
    group = c("NIH", "CBCL", "PPS", "UPPS", "BIS/BAS"),
    title = rowMeans(grid_data)
  )
g1 <-
  ggplot(behav.profile.longdat,
         aes(x = Feature, y = Value, fill = Cluster),
         alpha = 0.5) +
  geom_bar(stat = "identity") + ylim(-.5, 1) +
  facet_wrap(. ~ Cluster) + coord_polar() + theme_void() +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1")))
g2 <- g1 + geom_segment(
  data = grid_data,
  aes(
    x = start,
    y = -.1,
    xend = end,
    yend = -.1
  ),
  colour = "black",
  alpha = 1,
  size = .5,
  inherit.aes = FALSE
)
g3 <- g2 + geom_text(
  data = behav.profile.longdat,
  aes(
    x = Feature,
    y = Value + .1,
    label = Feature,
    angle = angle
  ),
  color = "black",
  fontface = "bold",
  alpha = 0.6,
  size = 2.5
)
g4 <-
  g3 + geom_text(
    data = base_data,
    aes(x = title, y = -.25, label = group),
    colour = "black",
    alpha = 0.8,
    size = 2,
    fontface = "bold",
    inherit.aes = FALSE
  )
g4



# Multiple Bar Plot of KSADS Diagnosis
abcd_ksads.suicide.label <- read.csv2("ksads_suicidal_labels2.csv")

x <- abcd_envir.bl[, c(1, 81:192)]
# names(x)[2:99]<-abcd_ksads.suicide.label$Names
ksads.longdat.bl = melt(
  data = x,
  id.vars = c("subjectkey"),
  variable.name = "Diagnosis",
  value.name = "Value"
)
ksads.longdat.bl <- na.omit(ksads.longdat.bl)
ksads.longdat.bl <-
  merge(ksads.longdat.bl, abcd_envir.bl[, c(1, 80)], by = "subjectkey", all = T)
ksads.longdat.bl <- ksads.longdat.bl[-1 * which(ksads.longdat.bl$Value == 0), ]
ksads.longdat.bl <-
  merge(ksads.longdat.bl,
        abcd_ksads.suicide.label,
        by.x = "Diagnosis",
        by.y = "Code")
rm(x)

ksads.longdat.Combid.bl <-
  ksads.longdat.bl[which(ksads.longdat.bl$Questionnaire == "Comorbidity"), ]
ksads.longdat.noCombid.bl <-
  ksads.longdat.bl[-1 * which(ksads.longdat.bl$Questionnaire == "Comorbidity"), ]
ksads.longdat.noCombid.bl[which(ksads.longdat.noCombid.bl$Class == "Self Injurious Behavior"), "Class"] <-
  "Self Injury"

x <- abcd_envir.fu2[, c(1, 64:175)]
# names(x)[2:99]<-abcd_ksads.suicide.label$Names
ksads.longdat.fu2 = melt(
  data = x,
  id.vars = c("subjectkey"),
  variable.name = "Diagnosis",
  value.name = "Value"
)
ksads.longdat.fu2 <- na.omit(ksads.longdat.fu2)
ksads.longdat.fu2 <-
  merge(ksads.longdat.fu2, abcd_envir.fu2[, c(1, 63)], by = "subjectkey", all = T)
ksads.longdat.fu2 <- ksads.longdat.fu2[-1 * which(ksads.longdat.fu2$Value == 0), ]
ksads.longdat.fu2 <-
  merge(ksads.longdat.fu2,
        abcd_ksads.suicide.label,
        by.x = "Diagnosis",
        by.y = "Code")
rm(x)

ksads.longdat.Combid.fu2 <-
  ksads.longdat.fu2[which(ksads.longdat.fu2$Questionnaire == "Comorbidity"), ]
ksads.longdat.noCombid.fu2 <-
  ksads.longdat.fu2[-1 * which(ksads.longdat.fu2$Questionnaire == "Comorbidity"), ]
ksads.longdat.noCombid.fu2[which(ksads.longdat.noCombid.fu2$Class == "Self Injurious Behavior"), "Class"] <-
  "Self Injury"


ggplot(ksads.longdat.Combid, aes(Names, fill = clust)) + geom_bar(position = "fill") +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(Wave ~ Class, scales = "free", space = "free") + 
  ylim(0.6,1)
  coord_flip()

ggplot(ksads.longdat.noCombid, aes(Names, fill = clust)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(Questionnaire ~ Class, scales = "free", space = "free")

ggplot(ksads.longdat.noCombid.fu2, aes(Names, fill = clust)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(Questionnaire ~ Class, scales = "free", space = "free")


ksads.longdat.bl$Wave <- rep("BL",length(ksads.longdat.bl),1)
ksads.longdat.fu2$Wave <- rep("FU2",length(ksads.longdat.fu2),1)
ksads.longdat <- rbind(ksads.longdat.bl,ksads.longdat.fu2)
ksads.longdat.Combid <-
  ksads.longdat[which(ksads.longdat$Questionnaire == "Comorbidity"), ]
ksads.longdat.noCombid <-
  ksads.longdat[-1 * which(ksads.longdat$Questionnaire == "Comorbidity"), ]
ksads.longdat.noCombid[which(ksads.longdat.noCombid$Class == "Self Injurious Behavior"), "Class"] <-
  "Self Injury"


xt<-xtabs(formula =  ~ clust+Class+Wave+Questionnaire,
          data = ksads.longdat.noCombid[
            -1*c(which(is.na(ksads.longdat.noCombid$clust)), 
                 which(ksads.longdat.noCombid$clust=="0")),])
xt


abcd_ksads_suicide <- merge(abcd_ksads_suicide,abcd_envir.bl[,c(1,80)],all = TRUE)
names(abcd_ksads_suicide)[101] <- "Cluster"

ksads.longdat.rate <- data.frame()
indx <- 0
for (j in unique(abcd_ksads.suicide.label[1:98,]$Class)[c(2, 1, 3, 5)]) {
  for (k in unique(abcd_ksads_suicide$eventname)[c(2, 1, 3)]) {
    for (l in unique(abcd_ksads.suicide.label[1:98,]$Questionnaire)) {
      x <- abcd_ksads_suicide[which(abcd_ksads_suicide$eventname == k),
                              c(
                                101,
                                which(
                                  abcd_ksads.suicide.label[1:98,]$Item == "Symptom" &
                                    abcd_ksads.suicide.label[1:98,]$Class == j &
                                    abcd_ksads.suicide.label[1:98,]$Questionnaire == l
                                ) + 2
                              )]
      for (g in c("1", "2", "3", "4", "5")) {
        indx <- indx + 1
        ksads.longdat.rate[indx, "Eventname"] <- k
        ksads.longdat.rate[indx, "Cluster"] <- as.character(g)
        ksads.longdat.rate[indx, "Class"] <- j
        ksads.longdat.rate[indx, "Questionnaire"] <- l
        ksads.longdat.rate[indx, "Rate"] <-
          sum(rowSums(x[which(x$Cluster == g), 2:ncol(x)], na.rm = TRUE) >
                0) / nrow(x[which(x$Cluster == g), ])
      }
      indx <- indx + 1
      ksads.longdat.rate[indx, "Eventname"] <- k
      ksads.longdat.rate[indx, "Cluster"] <- "All"
      ksads.longdat.rate[indx, "Class"] <- j
      ksads.longdat.rate[indx, "Questionnaire"] <- l
      ksads.longdat.rate[indx, "Rate"] <-
        sum(rowSums(x[-1 * which(x$Cluster == "0" |
                                   is.na(x$Cluster)), 2:ncol(x)], na.rm = TRUE) > 0) /
        nrow(x[-1 * which(x$Cluster == "0" | is.na(x$Cluster)), ])
    }
  }
}



ksads.longdat.rate[ksads.longdat.rate$Class=="Self Injurious Behavior","Class"] <- 
  rep("Self Injury",sum(ksads.longdat.rate$Class=="Self Injurious Behavior"),1)

ksads.longdat.rate[ksads.longdat.rate$Eventname=="baseline_year_1_arm_1","Eventname"] <- 
  rep("BL",sum(ksads.longdat.rate$Eventname=="baseline_year_1_arm_1"),1)
ksads.longdat.rate[ksads.longdat.rate$Eventname=="1_year_follow_up_y_arm_1","Eventname"] <- 
  rep("FU1",sum(ksads.longdat.rate$Eventname=="1_year_follow_up_y_arm_1"),1)
ksads.longdat.rate[ksads.longdat.rate$Eventname=="2_year_follow_up_y_arm_1","Eventname"] <- 
  rep("FU2",sum(ksads.longdat.rate$Eventname=="2_year_follow_up_y_arm_1"),1)

ksads.longdat.rate$Eventname <-
  factor(
    ksads.longdat.rate$Eventname,
    ordered = TRUE,
    levels = c(
      "BL",
      "FU1",
      "FU2"
    )
  )
ksads.longdat.rate$Class <-
  factor(
    ksads.longdat.rate$Class,
    ordered = TRUE,
    levels = c(
      "Self Injury",
      "Suicidal Ideation",
      "Suicidal Attempt",
      "Suicidal Behavior"
    )
  )

ksads.longdat.rate$Cluster <-
  factor(
    ksads.longdat.rate$Cluster,
    ordered = TRUE,
    levels = c(
      "All","1","2","3","4","5")
  )

ggplot(data = ksads.longdat.rate[-1*which(ksads.longdat.rate$Eventname == "FU1"&
                                            ksads.longdat.rate$Questionnaire == "K-SADS-P"), ],
       aes(
         x = Eventname,
         y = Rate,
         group = Cluster,
         colour = Cluster
       )) +
  geom_line(stat = "identity", size = .5) +
  geom_point(stat = "identity", size = 2) +
  facet_grid(Questionnaire ~ Class, scales = "free", space = "free") +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data = ksads.longdat.rate[-1*which(ksads.longdat.rate$Eventname == "FU1"&
                                            ksads.longdat.rate$Questionnaire == "K-SADS-P"), ],
       aes(
         x = Eventname,
         y = Rate,
         group = Cluster,
         colour = Cluster
       )) +
  geom_line(stat = "identity", size = .5) +
  geom_point(stat = "identity", size = 2) +
  facet_grid( ~Questionnaire +  Class, scales = "free", space = "free") +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(data = ksads.longdat.rate[-1*which(ksads.longdat.rate$Eventname == "FU1"&
                                            ksads.longdat.rate$Questionnaire == "K-SADS-P"), ],
       aes(
         x = Eventname,
         y = Rate,
         group = Cluster,
         colour = Cluster
       )) +
  geom_line(stat = "identity", size = .5) +
  geom_point(stat = "identity", size = 2) +
  facet_grid(Class ~ Questionnaire, scales = "free", space = "free") +
  scale_fill_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c("#B0B0B0", brewer.pal(5, "Set1"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ksads.longdat.rate.dcast <- dcast(ksads.longdat.rate, value.var = "Rate", 
                                  formula = Eventname + Cluster + Questionnaire ~ Class)
ksads.longdat.rate.dcast$Cluster <- paste("Subtype",ksads.longdat.rate.dcast$Cluster,sep = " ")
write.csv(ksads.longdat.rate.dcast,file = "ksads_longdat_rate.csv")



dsm_14disorder <- merge(dsm_14disorder,abcd_envir.bl[,c(1,80)],all = TRUE)
names(dsm_14disorder)[19] <- "Cluster"

ksads.longdat.Combid.rate <- data.frame()
indx <- 0
for (j in unique(abcd_ksads.suicide.label[99:112, ]$Code)) {
  for (k in unique(dsm_14disorder$eventname)[c(2, 1, 3)]) {
    for (l in unique(abcd_ksads.suicide.label[99:112, ]$Questionnaire)) {
      for (g in c("1", "2", "3", "4", "5")) {
        indx <- indx + 1
        ksads.longdat.Combid.rate[indx, "Eventname"] <- k
        ksads.longdat.Combid.rate[indx, "Cluster"] <-
          as.character(g)
        ksads.longdat.Combid.rate[indx, "Code"] <- j
        ksads.longdat.Combid.rate[indx, "Questionnaire"] <- l
        ksads.longdat.Combid.rate[indx, "Rate"] <-
          sum(dsm_14disorder[which(dsm_14disorder$Cluster == g&dsm_14disorder$eventname==k),j])/
          length(dsm_14disorder[which(dsm_14disorder$Cluster == g&dsm_14disorder$eventname==k),j])
      }
    }
  }
}

ksads.longdat.Combid.rate[ksads.longdat.Combid.rate$Eventname=="baseline_year_1_arm_1","Eventname"] <- 
  rep("BL",sum(ksads.longdat.Combid.rate$Eventname=="baseline_year_1_arm_1"),1)
ksads.longdat.Combid.rate[ksads.longdat.Combid.rate$Eventname=="1_year_follow_up_y_arm_1","Eventname"] <- 
  rep("FU1",sum(ksads.longdat.Combid.rate$Eventname=="1_year_follow_up_y_arm_1"),1)
ksads.longdat.Combid.rate[ksads.longdat.Combid.rate$Eventname=="2_year_follow_up_y_arm_1","Eventname"] <- 
  rep("FU2",sum(ksads.longdat.Combid.rate$Eventname=="2_year_follow_up_y_arm_1"),1)

ksads.longdat.Combid.rate$Eventname <-
  factor(
    ksads.longdat.Combid.rate$Eventname,
    ordered = TRUE,
    levels = c(
      "BL",
      "FU1",
      "FU2"
    )
  )
ksads.longdat.Combid.rate$Code <-
  factor(ksads.longdat.Combid.rate$Code,
         ordered = TRUE,
         levels = c(unique(abcd_ksads.suicide.label[99:112,]$Code)))

ksads.longdat.Combid.rate$Cluster <-
  factor(
    ksads.longdat.Combid.rate$Cluster,
    ordered = TRUE,
    levels = c("1", "2", "3", "4", "5")
  )

ksads.longdat.Combid.rate <- merge(ksads.longdat.Combid.rate,abcd_ksads.suicide.label[,c(3,5)])

ggplot(data = ksads.longdat.Combid.rate[-1*which(ksads.longdat.Combid.rate$Eventname=="FU1"),],
       aes(
         x = Eventname,
         y = Rate,
         group = Cluster,
         colour = Cluster
       )) +
  geom_line(stat = "identity", size = .5) +
  geom_point(stat = "identity", size = 2) +
  facet_grid(. ~ Code, scales = "free", space = "free") +
  scale_fill_manual(values = c(brewer.pal(5, "Set1"))) +
  scale_color_manual(values = c(brewer.pal(5, "Set1"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



