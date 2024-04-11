
abcd_ksads_suicide.bl$STB_P_score_BL<-rowSums(abcd_ksads_suicide.bl[,3:51],na.rm = T)
abcd_ksads_suicide.bl$STB_Y_score_BL<-rowSums(abcd_ksads_suicide.bl[,52:100],na.rm = T)

abcd_ksads_suicide.fu2$STB_P_score_FU2<-rowSums(abcd_ksads_suicide.fu2[,3:51],na.rm = T)
abcd_ksads_suicide.fu2$STB_Y_score_FU2<-rowSums(abcd_ksads_suicide.fu2[,52:100],na.rm = T)

abcd_KSADS_STB <- merge(abcd_ksads_suicide.bl[,c(1,101:102)],
                        abcd_ksads_suicide.fu2[,c(1,101:102)])

# names(abcd_KSADS_STB) <- abcd_KSADS_STB$subjectkey
abcd_KSADS_STB <- merge(abcd_trait.bl.clust[c(1,36)],abcd_KSADS_STB)

abcd_KSADS_STB$STB_P_Delta <- abcd_KSADS_STB$STB_P_score_FU2 - abcd_KSADS_STB$STB_P_score_BL
abcd_KSADS_STB$STB_Y_Delta <- abcd_KSADS_STB$STB_Y_score_FU2 - abcd_KSADS_STB$STB_Y_score_BL

abcd_KSADS_STB$STB_Sum_BL <- abcd_KSADS_STB$STB_P_score_BL + abcd_KSADS_STB$STB_Y_score_BL
abcd_KSADS_STB$STB_Sum_FU2 <- abcd_KSADS_STB$STB_P_score_FU2 + abcd_KSADS_STB$STB_Y_score_FU2

STB_Average = data.frame(STB_P_score_BL = NA,STB_Y_score_BL = NA,STB_P_score_FU2 = NA,
                         STB_Y_score_FU2 = NA,STB_P_Delta = NA,STB_Y_Delta = NA)
for (i in 1:5){
  STB_Average[i,] = colMeans(abcd_KSADS_STB[abcd_KSADS_STB$clust==i,c(3:8)])
}



library(ggpubr)
fit.wilcox.KSADS_P = data.frame(wilcox.V = rep(NA,5),wilcox.p = rep(NA,5))
fit.wilcox.KSADS_Y = data.frame(wilcox.V = rep(NA,5),wilcox.p = rep(NA,5))
fit.tstat.KSADS_P = data.frame(tstat = rep(NA,5),tstat.p = rep(NA,5))
fit.tstat.KSADS_Y = data.frame(tstat = rep(NA,5),tstat.p = rep(NA,5))
fit.tstat.KSADS_Sum = data.frame(tstat = rep(NA,5),tstat.p = rep(NA,5))
for (j in 1:5){
  
  fit = wilcox.test(x = abcd_KSADS_STB[abcd_KSADS_STB$clust==j,5], 
                    y = abcd_KSADS_STB[abcd_KSADS_STB$clust==j,3],
                    alternative = "two.sided",paired = T,exact = T)
  fit.wilcox.KSADS_P$wilcox.V[j] = fit$statistic
  fit.wilcox.KSADS_P$wilcox.p[j] = fit$p.value
  
  fit = wilcox.test(x = abcd_KSADS_STB[abcd_KSADS_STB$clust==j,6], 
                    y = abcd_KSADS_STB[abcd_KSADS_STB$clust==j,4], 
                    alternative = "two.sided",paired = T,exact = T)
  fit.wilcox.KSADS_Y$wilcox.V[j] = fit$statistic
  fit.wilcox.KSADS_Y$wilcox.p[j] = fit$p.value
  
  
  fit = t.test(x = abcd_KSADS_STB[abcd_KSADS_STB$clust==j,5], 
               y = abcd_KSADS_STB[abcd_KSADS_STB$clust==j,3],
                        alternative = "two.sided",paired = T)
  fit.tstat.KSADS_P$tstat[j] = fit$statistic
  fit.tstat.KSADS_P$tstat.p[j] = fit$p.value
  
  fit = t.test(x = abcd_KSADS_STB[abcd_KSADS_STB$clust==j,6], 
               y = abcd_KSADS_STB[abcd_KSADS_STB$clust==j,4],
               alternative = "two.sided",paired = T)
  fit.tstat.KSADS_Y$tstat[j] = fit$statistic
  fit.tstat.KSADS_Y$tstat.p[j] = fit$p.value
  
  fit = t.test(x = abcd_KSADS_STB[abcd_KSADS_STB$clust==j,10], 
               y = abcd_KSADS_STB[abcd_KSADS_STB$clust==j,9],
               alternative = "two.sided",paired = T)
  fit.tstat.KSADS_Sum$tstat[j] = fit$statistic
  fit.tstat.KSADS_Sum$tstat.p[j] = fit$p.value

  
}

model <- aov_car(dependent_variable ~ condition + Error(subject), data=data)


# 加载reshape2包
library(reshape2)
library(car)
# library(emmeans)
# 使用melt()函数将宽数据转换为长数据
abcd_KSADS_STB_Long <- melt(abcd_KSADS_STB[,c(1,4,6)], 
                            id.vars = "subjectkey", 
                            variable.name = "Time_point", 
                            value.name = "Measurement")
abcd_KSADS_STB_Long$Time <- 0
abcd_KSADS_STB_Long$Time[abcd_KSADS_STB_Long$Time_point=="STB_P_score_FU2"] <- 1
abcd_KSADS_STB_Long <- merge(abcd_KSADS_STB_Long,abcd_KSADS_STB[,c(1,2)])
abcd_KSADS_STB_Long$clust <- paste("Subtype",abcd_KSADS_STB_Long$clust)
abcd_KSADS_STB_Long$Measurement <- scale(abcd_KSADS_STB_Long$Measurement)
# ANOVA
fit = aov(data = abcd_KSADS_STB_Long, 
              Measurement ~ clust*Time + Error(subjectkey/Time))
summary(fit)

library(rstatix)
anova_test(data = abcd_KSADS_STB_Long,
           dv = Measurement,
           wid = subjectkey,
           within = Time,
           between = clust
)

library(ggplot2)
library(dplyr)
abcd_KSADS_STB_Long %>% 
  group_by(Time,clust) %>% 
  summarise(mm=mean(Measurement)) %>% 
  ggplot(aes(Time,mm))+
  geom_line(aes(group=clust,color=clust),size=1.2)+
  theme_bw()

library(PMCMRplus)

summary(lsdTest(Measurement ~ clust, data = abcd_KSADS_STB_Long))

abcd_KSADS_STB_Long %>% 
  group_by(clust) %>% 
  t_test(Measurement ~ Time, paired = T)
# ref.group = "HC"
abcd_KSADS_STB_Long %>% 
  t_test(Measurement ~ clust, paired = F)

TukeyHSD(fit)
# emmeans(fit, specs = ~clust:Time, type = "marginal")

fit.tstat.KSADS_P$tstat[j] = fit$statistic
fit.tstat.KSADS_P$tstat.p[j] = fit$p.value

coefficients(fit)


ggboxplot(abcd_KSADS_STB,
          y = "STB_P_Delta", 
          x = "clust",
          ylab = "Delta", xlab = FALSE,
          ggtheme = theme_minimal())
ggboxplot(abcd_KSADS_STB,
          y = "STB_Y_Delta", 
          x = "clust",
          ylab = "Delta", xlab = FALSE,
          ggtheme = theme_minimal())
