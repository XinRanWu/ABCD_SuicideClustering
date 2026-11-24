# Run after 01. Reviewer sensitivity analysis; not a validated replacement for original results.
# The original score included NSSI-labelled items and requires scientific verification.
# Do not use estimates for submission until coding and missingness are confirmed.
# Reviewer response: longitudinal analysis using the ORIGINAL score definition.
# Run from the project directory with Rscript. Inspect the item dictionary and
# coding before interpreting the score as a measure of suicidal severity.
library(lme4)
library(lmerTest)
library(emmeans)

load(prepared_file)
out <- 'results/longitudinal'
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# The original LongitudinalRate.R summed columns 3:51 and 52:100 with
# na.rm=TRUE. Report that exact calculation as the primary reproduction;
# an all-missing person-wave is flagged explicitly, not interpreted as zero.
items <- names(abcd_ksads_suicide.bl)[3:100]
stopifnot(identical(items, names(abcd_ksads_suicide.fu2)[3:100]))
labels <- if (exists('abcd_ksads.suicide.label') && length(abcd_ksads.suicide.label) >= 98) as.character(abcd_ksads.suicide.label[1:98]) else items
dictionary <- data.frame(reporter=rep(c('parent','youth'), each=49), variable=items, label=labels)
write.csv(dictionary, file.path(out,'score_item_dictionary.csv'), row.names=FALSE)

bl <- abcd_ksads_suicide.bl
fu <- abcd_ksads_suicide.fu2
bl$parent <- rowSums(bl[,3:51], na.rm=TRUE)
bl$youth <- rowSums(bl[,52:100], na.rm=TRUE)
fu$parent <- rowSums(fu[,3:51], na.rm=TRUE)
fu$youth <- rowSums(fu[,52:100], na.rm=TRUE)
bl$parent[rowSums(!is.na(bl[,3:51])) == 0] <- NA
bl$youth[rowSums(!is.na(bl[,52:100])) == 0] <- NA
fu$parent[rowSums(!is.na(fu[,3:51])) == 0] <- NA
fu$youth[rowSums(!is.na(fu[,52:100])) == 0] <- NA
key <- abcd_trait.bl.clust[,c('subjectkey','clust')]
stopifnot(!anyDuplicated(key$subjectkey), !anyDuplicated(bl$subjectkey), !anyDuplicated(fu$subjectkey))
wide <- merge(key, bl[,c('subjectkey','parent','youth')], by='subjectkey', all.x=TRUE)
wide <- merge(wide, fu[,c('subjectkey','parent','youth')], by='subjectkey', all.x=TRUE, suffixes=c('_bl','_fu'))
wide$clust <- factor(wide$clust, levels=1:5)
flow <- aggregate(cbind(parent_retained=!is.na(parent_fu), youth_retained=!is.na(youth_fu)) ~ clust, wide, sum)
flow$baseline_n <- as.integer(table(wide$clust)[as.character(flow$clust)])
write.csv(flow,file.path(out,'retention_by_profile.csv'),row.names=FALSE)

library(pbkrtest)

for (reporter in c('parent','youth')) {
  dat <- rbind(data.frame(subjectkey=wide$subjectkey, clust=wide$clust, time='baseline', score=wide[[paste0(reporter,'_bl')]]),
               data.frame(subjectkey=wide$subjectkey, clust=wide$clust, time='followup', score=wide[[paste0(reporter,'_fu')]]))
  dat$time <- factor(dat$time, levels=c('baseline','followup'))
  dat <- dat[complete.cases(dat),]
  fit <- lmer(score ~ clust*time + (1|subjectkey), data=dat, REML=TRUE)
  capture.output(summary(fit), anova(fit, type=3), file=file.path(out,paste0(reporter,'_model.txt')))
  # Compare within-profile follow-up minus baseline changes directly.
  changes <- contrast(emmeans(fit, ~ time | clust), method='revpairwise')
  comparisons <- contrast(changes, method='pairwise', by=NULL, adjust='holm')
  result <- as.data.frame(summary(comparisons, infer=c(TRUE,TRUE)))
  write.csv(result,file.path(out,paste0(reporter,'_direct_change_contrasts.csv')),row.names=FALSE)
  write.csv(as.data.frame(summary(changes, infer=c(TRUE,TRUE))),file.path(out,paste0(reporter,'_within_profile_change.csv')),row.names=FALSE)
}
cat('Finished. Verify item coding, missingness and model assumptions before reporting.\n')
