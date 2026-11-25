# Reporter-specific longitudinal analyses of original K-SADS diagnosis items.
# Present/past endorsements are combined at each visit: this is NOT incidence.
library(lme4);library(emmeans)
e<-new.env();load('Data_and_Scripts/abcd_SuicideData.RData',envir=e)
out<-'Revision_Analysis_20261122/results_stb_binary_20261125';dir.create(out,recursive=TRUE,showWarnings=FALSE)
key<-e$abcd_trait.bl.clust[,c('subjectkey','clust')];stopifnot(!anyDuplicated(key$subjectkey))
categories<-list(ideation=c(946:950,957:961),attempt=c(954,965),other_behavior=c(951:953,962:964))
description<-list();summaries<-list();contrast_tables<-list();item_map<-list()
for(reporter in c('p','t'))for(category in names(categories)){
 vars<-paste0('ksads_23_',categories[[category]],'_',reporter)
 stopifnot(all(vars %in% names(e$abcd_ksads_suicide.bl)),all(vars %in% names(e$abcd_ksads_suicide.fu2)))
 item_map[[paste(reporter,category)]]<-data.frame(reporter=reporter,category=category,variable=vars)
 bl<-e$abcd_ksads_suicide.bl[,c('subjectkey',vars)];fu<-e$abcd_ksads_suicide.fu2[,c('subjectkey',vars)]
 bl$wave<-'baseline';fu$wave<-'followup_2y';d<-merge(rbind(bl,fu),key,by='subjectkey')
 x<-as.matrix(d[,vars]);stopifnot(all(na.omit(as.vector(x)) %in% c(0,1)))
 # 1 if any item positive; 0 only when all items observed and none positive.
 d$outcome<-ifelse(rowSums(x==1,na.rm=TRUE)>0,1,ifelse(rowSums(is.na(x))==0,0,NA))
 d$wave<-factor(d$wave,levels=c('baseline','followup_2y'));d$clust<-factor(d$clust,levels=1:5)
 for(g in 1:5)for(w in levels(d$wave)){
   z<-d[d$clust==g & d$wave==w,'outcome'];n<-sum(!is.na(z));positive<-sum(z==1,na.rm=TRUE)
   description[[length(description)+1]]<-data.frame(reporter=reporter,category=category,profile=g,wave=w,observed=n,positive=positive,missing=sum(is.na(z)),proportion=if(n)positive/n else NA)
 }
 use<-d[complete.cases(d[,c('outcome','clust','wave','subjectkey')]),]
 cat('Fitting',reporter,category,'events',sum(use$outcome),'\n')
 fit<-tryCatch(glmer(outcome~clust*wave+(1|subjectkey),data=use,family=binomial,control=glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=200000))),error=function(e)e)
 if(inherits(fit,'error')){summaries[[length(summaries)+1]]<-data.frame(reporter=reporter,category=category,participants=length(unique(use$subjectkey)),observations=nrow(use),events=sum(use$outcome),singular=NA,convergence=conditionMessage(fit),interaction_p=NA);next}
 null<-tryCatch(update(fit,.~.-clust:wave),error=function(e)e)
 convergence<-paste(fit@optinfo$conv$lme4$messages,collapse='; ')
 lrt<-if(inherits(null,'error')) NA_real_ else anova(null,fit,test='Chisq')$`Pr(>Chisq)`[2]
 summaries[[length(summaries)+1]]<-data.frame(reporter=reporter,category=category,participants=length(unique(use$subjectkey)),observations=nrow(use),events=sum(use$outcome),singular=isSingular(fit),convergence=convergence,interaction_p=lrt)
 capture.output(summary(fit),if(!inherits(null,'error'))anova(null,fit,test='Chisq'),file=file.path(out,paste0(reporter,'_',category,'_model.txt')))
 cmp<-tryCatch({em<-emmeans(fit,~wave|clust);delta<-contrast(em,'revpairwise');as.data.frame(summary(contrast(delta,'pairwise',by=NULL,adjust='holm'),infer=c(TRUE,TRUE)))},error=function(e)NULL)
 if(!is.null(cmp)){cmp$reporter<-reporter;cmp$category<-category;contrast_tables[[length(contrast_tables)+1]]<-cmp}
}
write.csv(do.call(rbind,description),file.path(out,'stb_counts_by_profile_wave.csv'),row.names=FALSE)
write.csv(do.call(rbind,summaries),file.path(out,'stb_model_diagnostics.csv'),row.names=FALSE)
if(length(contrast_tables))write.csv(do.call(rbind,contrast_tables),file.path(out,'stb_change_contrasts.csv'),row.names=FALSE)
write.csv(do.call(rbind,item_map),file.path(out,'stb_item_dictionary.csv'),row.names=FALSE)
writeLines(c('Original-version diagnostic items; excluded NSSI 945/956 and negative items 955/966.',
 'Any present OR past endorsement at each visit; other_behavior includes preparations, interrupted and aborted attempts, while attempts are analyzed separately.',
 'A positive item defines a positive outcome even with partial missingness; a negative outcome requires all category items observed.',
 'Do not infer onset or persistence: past and current item windows may overlap across visits.',
 'Models are exploratory unadjusted logistic random-intercept models, on log-odds scale; Holm correction within reporter/outcome across 10 contrasts.',
 'Sparse events, separation, singularity and convergence must be inspected before drawing conclusions.'),file.path(out,'notes.txt'))
