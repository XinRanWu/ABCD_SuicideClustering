########################################################
#                 Function Defination
########################################################
require("pacman")
pacman::p_load(
  cluster,factoextra,psych,ggradar,pheatmap,ggplot2,ggiraph,
  ggiraphExtra,RColorBrewer,tableone,glmnet,rgl,ggseg,dplyr,
  reshape2,Rmisc,ggpubr,export,Rcpp,ggridges,stringr,ggfun,
  gghalves
)
# devtools::install_github("tomwenseleers/export")

read.abcd.table <- function(abcd.table.path) {
  res.table <-
    read.csv(abcd.table.path,
             sep = "\t",
             header = TRUE,
             skip = 1)
  att <-
    attributes(read.csv(abcd.table.path, sep = "\t", header = TRUE))
  res.short_label <- att$name
  colnames(res.table) <- att$name
  return(res.table)
}

read.abcd.table.long_label <- function(abcd.table.path) {
  res.table <-
    read.table(abcd.table.path,
               sep = "\t",
               header = TRUE,
               skip = 1)
  res_att <- attributes(res.table)
  res.long_label <- res_att$names
  return(res.long_label)
}

########################################################
#                      Covaribles
########################################################

abcd_covariables_bl_envir <-
  read.csv("abcd_covariables_bl_envir.csv", na.strings = "NaN")
abcd_covariables <-
  read.csv("abcd_covariables_full.csv", na.strings = "NaN")

abcd_ksad01 <-
  read.abcd.table("F:/Acdamic/Data_ABCD/abcd_ksad01.txt")
abcd_ksad01.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/abcd_ksad01.txt")
abcd_ksad501 <-
  read.abcd.table("F:/Acdamic/Data_ABCD/abcd_ksad501.txt")
abcd_ksad501.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/abcd_ksad501.txt")
abcd_cbcls01 <-
  read.abcd.table("F:/Acdamic/Data_ABCD/abcd_cbcls01.txt")
abcd_tbss01 <-
  read.abcd.table("F:/Acdamic/Data_ABCD/abcd_tbss01.txt")
abcd_mhy02 <- read.abcd.table("F:/Acdamic/Data_ABCD/abcd_mhy02.txt")
abcd_mhy02.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/abcd_mhy02.txt")
# abcd_mhp02<-read.abcd.table("F:/Acdamic/Data_ABCD/abcd_mhp02.txt")
# abcd_mhp02.label<-read.abcd.table.long_label("F:/Acdamic/Data_ABCD/abcd_mhp02.txt")

abcd_tbss01 <-
  abcd_tbss01[, c(4, 9, which(stringr::str_detect(names(abcd_tbss01), "_uncorrected")))]
abcd_mhy02 <-
  abcd_mhy02[, c(
    "subjectkey",
    "eventname",
    "pps_y_ss_severity_score",
    "upps_y_ss_negative_urgency",
    "upps_y_ss_lack_of_planning",
    "upps_y_ss_sensation_seeking",
    "upps_y_ss_positive_urgency",
    "upps_y_ss_lack_of_perseverance" ,
    "bis_y_ss_bis_sum",
    "bis_y_ss_bas_rr",
    "bis_y_ss_bas_drive",
    "bis_y_ss_bas_fs"
  )]
abcd_cbcls01 <- abcd_cbcls01[, c(4, 9, seq(10, 89, 4))]

dsm_14disorder <- read.csv("dsm_14disorder.csv", na.strings = NaN)
dsm_14disorder.bl <-
  dsm_14disorder[dsm_14disorder$eventname == "baseline_year_1_arm_1",]
# screen time
abcd_stq01 <- read.abcd.table("F:/Acdamic/Data_ABCD/abcd_stq01.txt")
abcd_stq01.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/abcd_stq01.txt")
abcd_stq01.trans <- transform(
  abcd_stq01[, c(1:23)],
  screen1_y = screen1_wkdy_y * 5 + screen7_wknd_y *
    2,
  screen2_y = screen2_wkdy_y * 5 + screen8_wknd_y *
    2,
  screen3_y = screen3_wkdy_y * 5 + screen9_wknd_y *
    2,
  screen4_y = screen4_wkdy_y * 5 + screen10_wknd_y *
    2,
  screen5_y = screen5_wkdy_y * 5 + screen11_wknd_y *
    2,
  screen6_y = screen_wkdy_y * 5 + screen12_wknd_y *
    2
)
abcd_stq01 <- abcd_stq01.trans[, c(4, 9, 22:29)]
rm(list = c("abcd_stq01.trans"))

prs_prscs <- read.csv("prs_prscs.csv")
prs_prsice <- read.csv("prs_prsice.csv")
names(prs_prscs)[1] <- "subjectkey"
names(prs_prsice)[1] <- "subjectkey"

abcd_imgincl01 <-
  read.abcd.table("F:/Acdamic/Data_ABCD/abcd_imgincl01.txt")
abcd_smrip101 <-
  read.abcd.table("F:/Acdamic/Data_ABCD/abcd_smrip101.txt")
abcd_smrip201 <-
  read.abcd.table("F:/Acdamic/Data_ABCD/abcd_smrip201.txt")
abcd_betnet02 <-
  read.abcd.table("F:/Acdamic/Data_ABCD/abcd_betnet02.txt")
mrirscor02 <- read.abcd.table("F:/Acdamic/Data_ABCD/mrirscor02.txt")
abcd_dmdtifp101 <-
  read.abcd.table("F:/Acdamic/Data_ABCD/abcd_dmdtifp101.txt")
abcd_dmdtifp202 <-
  read.abcd.table("F:/Acdamic/Data_ABCD/abcd_dmdtifp202.txt")
midaparc03 <- read.abcd.table("F:/Acdamic/Data_ABCD/midaparc03.txt")
midaparcp203 <-
  read.abcd.table("F:/Acdamic/Data_ABCD/midaparcp203.txt")
mrisst02 <- read.abcd.table("F:/Acdamic/Data_ABCD/mrisst02.txt")
nback_bwroi02 <-
  read.abcd.table("F:/Acdamic/Data_ABCD/nback_bwroi02.txt")

abcd_imgincl01.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/abcd_imgincl01.txt")
abcd_smrip101.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/abcd_smrip101.txt")
abcd_smrip201.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/abcd_smrip201.txt")
abcd_betnet02.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/abcd_betnet02.txt")
mrirscor02.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/mrirscor02.txt")
abcd_dmdtifp101.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/abcd_dmdtifp101.txt")
abcd_dmdtifp202.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/abcd_dmdtifp202.txt")
midaparc03.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/midaparc03.txt")
midaparcp203.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/midaparcp203.txt")
mrisst02.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/mrisst02.txt")
nback_bwroi02.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/nback_bwroi02.txt")

abcd_smrip <-
  merge(abcd_smrip101[, c(
    4,
    9,
    which(str_detect(
      colnames(abcd_smrip101), "smri_thick_cdk_"
    )),
    which(str_detect(
      colnames(abcd_smrip101), "smri_area_cdk_"
    )),
    which(str_detect(
      colnames(abcd_smrip101), "smri_vol_cdk_"
    )),
    which(str_detect(
      colnames(abcd_smrip101), "smri_t1wcnt_cdk_"
    ))
  )],
  abcd_smrip201[, c(4, 456, which(str_detect(
    colnames(abcd_smrip201), "smri_vol_scs_"
  )))],
  by = c("subjectkey", "eventname"))
abcd_smrip <-
  merge(abcd_imgincl01[, c(4, 9, 11)], abcd_smrip, by = c("subjectkey", "eventname"))
abcd_smrip <-
  abcd_smrip[-1 * which(abcd_smrip$imgincl_t1w_include != 1),]
abcd_smrip <- abcd_smrip[,-3]

abcd_rsfc <-
  merge(abcd_betnet02[, c(4, 9, 16, 23:191)], mrisst02[, c(4, 9, 23:269)], by =
          c("subjectkey", "eventname"))
abcd_rsfc <-
  merge(abcd_imgincl01[, c(4, 9, 14)], abcd_rsfc, by = c("subjectkey", "eventname"))
abcd_rsfc <-
  abcd_rsfc[-1 * which(abcd_rsfc$imgincl_rsfmri_include != 1),]
abcd_rsfc <- abcd_rsfc[,-3]


abcd_tract <-
  merge(abcd_dmdtifp101[, c(4, 9, which(str_detect(
    abcd_dmdtifp101.label, "within.DTI.atlas.tract"
  )))],
  abcd_dmdtifp202[, c(4, 9, which(str_detect(
    abcd_dmdtifp202.label, "within.DTI.atlas.tract"
  )))],
  by = c("subjectkey", "eventname"))
abcd_tract <-
  merge(abcd_dmdtifp202[, c("subjectkey", "eventname", "dmdtifp1_1183")], abcd_tract, by = c("subjectkey", "eventname"))

abcd_tract.label <-
  c(abcd_dmdtifp101.label[c(which(str_detect(
    abcd_dmdtifp101.label, "within.DTI.atlas.tract"
  )))],
  abcd_dmdtifp202.label[which(str_detect(abcd_dmdtifp101.label, "within.DTI.atlas.tract"))])

abcd_tract <-
  merge(abcd_imgincl01[, c(4, 9, 13)], abcd_tract, by = c("subjectkey", "eventname"))
abcd_tract <-
  abcd_tract[-1 * which(abcd_tract$imgincl_dmri_include != 1),]
abcd_tract <- abcd_tract[,-3]


mri_rsi_p102.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/mri_rsi_p102.txt")
mri_rsi_p202.label <-
  read.abcd.table.long_label("F:/Acdamic/Data_ABCD/mri_rsi_p202.txt")

mri_rsi_p102 <-
  read.abcd.table("F:/Acdamic/Data_ABCD/mri_rsi_p102.txt")
mri_rsi_p202 <-
  read.abcd.table("F:/Acdamic/Data_ABCD/mri_rsi_p202.txt")

mri_rsi_mode <-
  c("rsin0gm",
    "rsin0s2gm",
    "rsindgm",
    "rsinds2gm",
    "rsintgm",
    "rsints2gm")
mri_rsi <-
  merge(mri_rsi_p102[, c(4, 9, 11:913)], mri_rsi_p202[, c(4, 9, 9:861)])
mri_rsi <- 
  mri_rsi[,-1*which(str_detect(colnames(mri_rsi), ".1"))]

abcd_rsi <-
  mri_rsi[, c(
    1:5,
    which(str_detect(colnames(mri_rsi), "gm_cdk_")),
    which(str_detect(colnames(mri_rsi), "_scs_"))
  )]

abcd_rsi <-
  merge(abcd_imgincl01[, c(4, 9, 13)], abcd_rsi, by = c("subjectkey", "eventname"))
abcd_rsi <-
  abcd_rsi[-1 * which(abcd_rsi$imgincl_dmri_include != 1),]
abcd_rsi <- abcd_rsi[,-3]



rm(
  list = c(
    "abcd_smrip101",
    "abcd_smrip201",
    "abcd_betnet02",
    "mrirscor02",
    "abcd_imgincl01"
  )
)
rm(list = c("abcd_dmdtifp101", "abcd_dmdtifp202"))
rm(list = c("mri_rsi_p102", "mri_rsi_p202", "mri_rsi"))

cor_dic = c("tp", "caudate", "putamen", "pallidum", "hpus", "amygdala", "aa")
cdk_dic = c(
  "bankssts",
  "cdacate",
  "cdmdfr",
  "cuneus",
  "ehinal",
  "fusiform",
  "ifpl",
  "iftm",
  "ihcate",
  "locc",
  "lobfr",
  "lingual",
  "mobfr",
  "mdtm",
  "parahpal",
  "paracn",
  "parsopc",
  "parsobis",
  "parstgris",
  "pericc",
  "postcn",
  "ptcate",
  "precn",
  "pc",
  "rracate",
  "rrmdfr",
  "sufr",
  "supl",
  "sutm",
  "sm",
  "frpole",
  "tmpole",
  "trvtm",
  "insula"
)

dk_dic <-
  c(
    "bankssts",
    "caudal anterior cingulate",
    "caudal middle frontal",
    "cuneus",
    "entorhinal",
    "fusiform",
    "inferior parietal",
    "inferior temporal",
    "isthmus cingulate",
    "lateral occipital",
    "lateral orbitofrontal",
    "lingual",
    "medial orbitofrontal",
    "middle temporal",
    "parahippocampal",
    "paracentral",
    "pars opercularis",
    "pars orbitalis",
    "pars triangularis",
    "pericalcarine",
    "postcentral",
    "posterior cingulate",
    "precentral",
    "precuneus",
    "rostral anterior cingulate",
    "rostral middle frontal",
    "superior frontal",
    "superior parietal",
    "superior temporal",
    "supramarginal",
    "frontal pole",
    "temporal pole",
    "transverse temporal",
    "insula"
  )


abcd_smrip <-
  abcd_smrip[,-1 * which(stringr::str_detect(names(abcd_smrip), "lesion"))]
abcd_smrip <-
  abcd_smrip[,-1 * which(stringr::str_detect(names(abcd_smrip), "wmhint"))]

# abcd_smrip.lr<-cbind(abcd_smrip[,c(1:2,251)],(abcd_smrip[,which(stringi::stri_detect_fixed(names(abcd_smrip),"lh"))] +
#                                          abcd_smrip[,which(stringi::stri_detect_fixed(names(abcd_smrip),"rh"))])/2)
# names(abcd_smrip.lr)<-gsub("lh","",names(abcd_smrip.lr))




mid_dic = c(
  'acdn',
  'acvn',
  'rpvnfb',
  'lpvnfb',
  'alrvn',
  'asrvn',
  'alvcr',
  'aclvn',
  'acmvn',
  'acgml'
)
mid_dic2 = c('arvn',
             'acvn',
             'rpvnfb',
             'lvnfb',
             'alrvn',
             'asrvn',
             'alvsr',
             'allvn',
             'asvn',
             'alvsl')
# anticipation (reward (large / small) / neutral / loss)
# feedback (positive / negative)
# arlvn allvn

# anticipation of reward versus neutral contrast
# anticipation of loss versus neutral contrast
# reward positive versus negative feedback contrast
# loss positive versus negative feedback contrast
# anticipation of large reward versus neutral contrast
# anticipation of small reward versus neutral contrast
# anticipation of large reward versus small reward contrast
# anticipation of large loss versus neutral contrast
# anticipation of small loss versus neutral contrast
# anticipation of large versus small loss contrast

abcd_MID <-
  merge(abcd_imgincl01[, c(4, 9, 15)], midaparc03[, c(4, 9, 16, 18, 20, 22:321)], by = c('subjectkey', 'eventname'))

abcd_MID <-
  merge(abcd_MID, midaparcp203[, c(4, 9, 10:689)], by = c('subjectkey', 'eventname'))

abcd_MID <-
  abcd_MID[-1 * which(abcd_MID$imgincl_mid_include != 1),]
abcd_MID <- abcd_MID[,-3]

nuclei_id <- c(6, 7, 8, 9, 13, 14, 16, 23:29)

con_id <- c(5, 8, 3, 4)
surf_id <- c()
scs_id <- c()

for (i in 1:length(con_id)) {
  surf_id <-
    c(surf_id, which(str_detect(
      names(abcd_MID), paste(mid_dic2[con_id[i]], '_b_cds_', sep = "")
    )))
  scs_id <-
    c(scs_id, which(str_detect(
      names(abcd_MID), paste(mid_dic[con_id[i]], '_b_scs_', sep = "")
    ))[nuclei_id])
  
}

abcd_BrainVar.bl <-
  abcd_MID[abcd_MID$eventname == "baseline_year_1_arm_1", c(4, surf_id, scs_id)]


# nBack 0 back condition
# nBack 2 back condition
# nBack place condition
# nBack emotion condition
# nBack 2 back versus 0 back contrast
# nBack face versus place contrast
# nBack emotion versus neutral face contrast
# nBack negative face versus neutral face contrast
# nBack positive face versus neutral face contrast

nback_dic = c('sacgvf',
              'sacsvcg',
              'saisvcg',
              'saasvcg',
              'sacsvis',
              'saigvcg',
              'saigvis')
nback_dic2 = c(
  'nBack 0 back condition',
  'nBack 2 back condition',
  'nBack place condition',
  'nBack emotion condition',
  'nBack 2 back versus 0 back contrast',
  'nBack face versus place contrast',
  'nBack emotion versus neutral face contrast',
  'nBack negative face versus neutral face contrast',
  'nBack positive face versus neutral face contrast'
)

abcd_nBack <-
  merge(abcd_imgincl01[, c(4, 9, 16)],
        nback_bwroi02[, c(4, 9, 16, 18, 20, 22:903)],
        by = c('subjectkey', 'eventname'))

abcd_nBack <-
  abcd_nBack[-1 * which(abcd_nBack$imgincl_nback_include != 1),]
abcd_nBack <- abcd_nBack[,-3]

con_id <- c(5, 8, 9)

# SST correct go versus fixation contrast
# SST correct stop versus correct go contrast
# SST incorrect stop versus correct go contrast
# SST any stop versus correct go contrast
# SST correct stop versus incorrect stop contrast
# SST incorrect go versus correct go contrast
# SST incorrect go versus incorrect stop contrast
sst_dic = c('sacgvf',
            'sacsvcg',
            'saisvcg',
            'saasvcg',
            'sacsvis',
            'saigvcg',
            'saigvis')
nuclei_id_sst <- c(6, 7, 8, 9, 13, 14, 16, 23:29)



abcd_SST <-
  merge(abcd_imgincl01[, c(4, 9, 17)], mrisst02[, c(4, 9, 16, 18, 20, 22:707)], by = c('subjectkey', 'eventname'))

abcd_SST <-
  abcd_SST[-1 * which(abcd_SST$imgincl_sst_include != 1),]
abcd_SST <- abcd_SST[,-3]

con_id <- c(2)

rm(list = c(
  "midaparc03",
  "midaparcp203",
  "mrisst02",
  "nback_bwroi02",
  "abcd_imgincl01"
))

########################################################
#                  Preparing Clustering
########################################################

abcd_ksad01[abcd_ksad01 == 555] <- NA
abcd_ksad01[abcd_ksad01 == 888] <- NA

abcd_ksad01.suicide.label <-
  abcd_ksad01.label[which(stringr::str_detect(names(abcd_ksad01), "ksads_23_"))]
abcd_ksad01_suicide <-
  abcd_ksad01[, c(4, 9, which(stringr::str_detect(names(abcd_ksad01), "ksads_23_")))]

abcd_ksad501[abcd_ksad501 == 555] <- NA
abcd_ksad501[abcd_ksad501 == 888] <- NA

abcd_ksad501.suicide.label <-
  abcd_ksad501.label[which(stringr::str_detect(names(abcd_ksad501), "ksads_23_"))]
abcd_ksad501_suicide <-
  abcd_ksad501[, c(4, 9, which(stringr::str_detect(names(abcd_ksad501), "ksads_23_")))]


abcd_ksads_suicide <-
  merge(abcd_ksad01_suicide, abcd_ksad501_suicide)
abcd_ksads.suicide.label <-
  c(abcd_ksad01.suicide.label, abcd_ksad501.suicide.label)

abcd_ksads_suicide.bl <-
  abcd_ksads_suicide[abcd_ksads_suicide$eventname == "baseline_year_1_arm_1",]
abcd_ksads_suicide.fu1 <-
  abcd_ksads_suicide[abcd_ksads_suicide$eventname == "1_year_follow_up_y_arm_1",]
abcd_ksads_suicide.fu2 <-
  abcd_ksads_suicide[abcd_ksads_suicide$eventname == "2_year_follow_up_y_arm_1",]

sum(rowSums(abcd_ksads_suicide.bl[, which(stringr::str_detect(abcd_ksads.suicide.label, "Diagnosis")) +
                                    2], na.rm = T) > 0)

abcd_trait <-
  merge(
    abcd_tbss01[, c(1:9)],
    abcd_cbcls01[, c(1:10, 14:22)],
    keys = c("subjectkey", "eventname"),
    all = T
  )
abcd_trait <-
  merge(abcd_trait,
        abcd_mhy02,
        keys = c("subjectkey", "eventname"),
        all = T)
abcd_trait[abcd_trait == 888] <- NA
abcd_trait[abcd_trait == 999] <- NA
abcd_trait[abcd_trait == 555] <- NA

diag_id = which(stringr::str_detect(abcd_ksads.suicide.label, "Diagnosis"))
abcd_ksads.suicide.label.diag <- abcd_ksads.suicide.label[diag_id]
del_id <-
  c(
    which(
      stringr::str_detect(abcd_ksads.suicide.label.diag, "SelfInjurious")
    ),
    which(
      stringr::str_detect(
        abcd_ksads.suicide.label.diag,
        "Nopastsuicidalideationorbehavior"
      )
    ),
    which(
      stringr::str_detect(
        abcd_ksads.suicide.label.diag,
        "Nosuicidalideationorbehavior"
      )
    )
  )

diag_id <- diag_id[-1 * del_id]

# abcd_ksads_suicide.bl$ksads01_23_diagnosis<-as.integer(rowSums(abcd_ksads_suicide.bl[,diag_id[1:18]+2],na.rm = T)>=1)
# abcd_ksads_suicide.bl$ksads501_23_diagnosis<-as.integer(rowSums(abcd_ksads_suicide.bl[,diag_id[19:36]+2],na.rm = T)>=1)
abcd_ksads_suicide.bl$ksads_23_diagnosis <-
  as.integer(rowSums(abcd_ksads_suicide.bl[, diag_id[1:18] + 2] +
                       abcd_ksads_suicide.bl[, diag_id[19:36] + 2],
                     na.rm = T) >= 1)

abcd_trait.bl <-
  abcd_trait[abcd_trait$eventname == "baseline_year_1_arm_1",]
abcd_trait.bl <-
  merge(abcd_trait.bl, abcd_ksads_suicide.bl[, c("subjectkey", "ksads_23_diagnosis")])
abcd_trait.bl[which(abcd_trait.bl$ksads_23_diagnosis == 1), 1:36]
rownames(abcd_trait.bl) <- abcd_trait.bl$subjectkey
abcd_trait.bl.clust <-
  abcd_trait.bl[which(abcd_trait.bl$ksads_23_diagnosis == 1), c(1, 3:36)]
abcd_trait.bl.clust <-
  abcd_trait.bl.clust[-1 * which(is.na(rowSums(abcd_trait.bl.clust[, 2:35]))),]

########################################################
#                      Clustering
########################################################

abcd_suicide.clust <-
  eclust(
    scale(abcd_trait.bl.clust[, c(2:35)]),
    nboot = 100,
    FUNcluster = c("kmeans"),
    iter.max = 50
  )#

sigclust(as.matrix(scale(
  abcd_trait.bl.clust[, c(2:35)])))

pvalueSum = rep(0,5)
for (i in 1:5){
  pvalueSum[i] <- sigclust(as.matrix(scale(
    abcd_trait.bl.clust[, c(2:35)])),nsim = i)@pvalnorm
}
pvalue <- sigclust(as.matrix(scale(
  abcd_trait.bl.clust[, c(2:35)])),n_start = 5)


library(CancerSubtypes)
pvalue <- CancerSubtypes::sigclustTest(as.matrix(scale(
  abcd_trait.bl.clust[, c(2:35)])),n_sim = 100,group = sigclustTest)

pvalue$p_norm
x = as.matrix(scale(abcd_trait.bl.clust[, c(2:35)]))
pvalues = matrix(rep(NA,25,1),5,5)
for (i in 1:5){
  for (j in 1:5){
    if (i != j){
      pvalues[i,j] <- sigclust(
        x[which(abcd_suicide.clust$cluster==i|abcd_suicide.clust$cluster==j),],
        nsim = 100,
        label = abcd_suicide.clust$cluster[which(
          abcd_suicide.clust$cluster==i|abcd_suicide.clust$cluster==j)]
        )@pvalnorm
    }
  }
}

clusters <- matrix(rep(NA,nrow(abcd_trait.bl.clust)*100,1),nrow(abcd_trait.bl.clust),100)
for (i in 1:100){
  set.seed(i)
  clust <-
    eclust(scale(abcd_trait.bl.clust[, c(2:35)]),
      nboot = 5,FUNcluster = c("kmeans"),iter.max = 50)#
  clusters[,i] <- clust$cluster
  
}

install.packages("pdfCluster")
library(pdfCluster)
adj.rand.index(cl1, cl2)

ARI = matrix(rep(NA,10000),100,100)
for (i in 1:100){
  for (j in 1:100){
    if (i != j){
      ARI[i,j] <- adj.rand.index(clusters[,i], clusters[,j])
    }
  }
}

fviz_cluster(abcd_suicide.clust,
             axes = c(1, 2),
             geom = c("point")) +
  ggplot2::theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1")
fviz_cluster(abcd_suicide.clust,
             axes = c(2, 3),
             geom = c("point")) +
  ggplot2::theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1")

abcd_trait.bl.clust.principal <-
  principal(abcd_trait.bl.clust[, 2:35], nfactors = 3)
p3d <- plot3d(
  abcd_trait.bl.clust.principal$scores[, 1],
  abcd_trait.bl.clust.principal$scores[, 2],
  abcd_trait.bl.clust.principal$scores[, 3],
  col = brewer.pal(5, 'Set1')[abcd_suicide.clust$cluster],
  xlab = "Dim 1 (29.1%)",
  ylab = "Dim 2 (9.8%)",
  zlab = "Dim 3 (7.1%)",
  size = 1,
  type = "s",
  alpha = 1,
  shininess = 80,
  smooth = FALSE,
  lit = FALSE,
  specular = "white",
  emission = "black"
)#type = c("shade","wire","dots")
bmp(filename = "p3d_cluster.bmp") # Spectral Pastel2 Paired Accent

fviz_gap_stat(abcd_suicide.clust$gap_stat) + ggplot2::theme_bw()
fviz_silhouette(abcd_suicide.clust) # silhouette plot








abcd_suicide.clust.part1 <-
  eclust(
    scale(abcd_trait.bl.clust[1:812, c(2:35)]),
    nboot = 100,
    FUNcluster = c("kmeans"),
    iter.max = 50
  )#

abcd_suicide.clust.part2 <-
  eclust(
    scale(abcd_trait.bl.clust[813:1624, c(2:35)]),
    nboot = 100,
    FUNcluster = c("kmeans"),
    iter.max = 50
  )#


fviz_nbclust(scale(abcd_trait.bl.clust[, c(2:35)]),FUNcluster = kmeans,method = c( "gap_stat"))

corrplot::corrplot(-1*cor(t(abcd_suicide.clust.part1$centers[order(c(4,2,1,5,3)),]),
                          t(abcd_suicide.clust.part2$centers[order(c(4,3,5,2,1)),])))

corrplot::corrplot(-1*cor(t(abcd_suicide.clust$centers),
                          t(abcd_suicide.clust.part1$centers[order(c(4,2,1,5,3)),])))
corrplot::corrplot(-1*cor(t(abcd_suicide.clust$centers),
                          t(abcd_suicide.clust.part2$centers[order(c(4,3,5,2,1)),])))

fviz_cluster(abcd_suicide.clust.part1,
             axes = c(1, 2),
             geom = c("point")) +
  ggplot2::theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1")
fviz_cluster(abcd_suicide.clust.part1,
             axes = c(2, 3),
             geom = c("point")) +
  ggplot2::theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1")

abcd_trait.bl.clust.principal.part1 <-
  principal(abcd_trait.bl.clust[1:812, 2:35], nfactors = 3)
abcd_trait.bl.clust.principal.part1$values[1:3]/sum(abcd_trait.bl.clust.principal.part1$values)
p3d <- plot3d(
  abcd_trait.bl.clust.principal.part1$scores[, 1],
  abcd_trait.bl.clust.principal.part1$scores[, 2],
  abcd_trait.bl.clust.principal.part1$scores[, 3],
  col = brewer.pal(5, 'Set1')[c(4,2,1,5,3)][abcd_suicide.clust.part1$cluster],
  xlab = "Dim 1 (29.7%)",
  ylab = "Dim 2 (9.7%)",
  zlab = "Dim 3 (6.9%)",
  size = 1,
  type = "s",
  alpha = 1,
  shininess = 80,
  smooth = FALSE,
  lit = FALSE,
  specular = "white",
  emission = "black"
)#type = c("shade","wire","dots")

abcd_trait.bl.clust.principal.part2 <-
  principal(abcd_trait.bl.clust[813:1624, 2:35], nfactors = 3)
abcd_trait.bl.clust.principal.part2$values[1:3]/sum(abcd_trait.bl.clust.principal.part2$values)
p3d <- plot3d(
  abcd_trait.bl.clust.principal.part2$scores[, 1],
  abcd_trait.bl.clust.principal.part2$scores[, 2],
  abcd_trait.bl.clust.principal.part2$scores[, 3],
  col = brewer.pal(5, 'Set1')[c(4,3,5,2,1)][abcd_suicide.clust.part2$cluster],
  xlab = "Dim 1 (28.7%)",
  ylab = "Dim 2 (10.1%)",
  zlab = "Dim 3 (7.4%)",
  size = 1,
  type = "s",
  alpha = 1,
  shininess = 80,
  smooth = FALSE,
  lit = FALSE,
  specular = "white",
  emission = "black"
)#type = c("shade","wire","dots")







abcd_ksads_suicide.bl$ksads_p_diagnosis <-
  as.integer(rowSums(abcd_ksads_suicide.bl[, diag_id[1:18] + 2],
                     na.rm = T) >= 1)
abcd_ksads_suicide.bl$ksads_y_diagnosis <-
  as.integer(rowSums(abcd_ksads_suicide.bl[, diag_id[19:36] + 2],
                     na.rm = T) >= 1)

abcd_trait.bl <-
  merge(abcd_trait.bl, abcd_ksads_suicide.bl[, c("subjectkey", "ksads_p_diagnosis","ksads_y_diagnosis")])
abcd_trait.bl[which(abcd_trait.bl$ksads_23_diagnosis == 1), 1:36]
rownames(abcd_trait.bl) <- abcd_trait.bl$subjectkey
abcd_trait.bl.clust <-
  abcd_trait.bl[which(abcd_trait.bl$ksads_23_diagnosis == 1), c(1, c(3:36,38,39))]
abcd_trait.bl.clust <-
  abcd_trait.bl.clust[-1 * which(is.na(rowSums(abcd_trait.bl.clust[, 2:37]))),]



abcd_suicide.clust.ksadsY <-
  eclust(
    scale(abcd_trait.bl.clust[abcd_trait.bl.clust$ksads_y_diagnosis==1, c(2:35)]),
    nboot = 100,
    FUNcluster = c("kmeans"),
    iter.max = 50
  )#

abcd_suicide.clust.ksadsP <-
  eclust(
    scale(abcd_trait.bl.clust[abcd_trait.bl.clust$ksads_p_diagnosis==1, c(2:35)]),
    nboot = 100,
    FUNcluster = c("kmeans"),
    iter.max = 50
  )#


fviz_nbclust(scale(abcd_trait.bl.clust[, c(2:35)]),FUNcluster = kmeans,method = c( "gap_stat"))

corrplot::corrplot(-1*cor(t(abcd_suicide.clust.ksadsP$centers[order(c(4,2,1,5,3)),]),
                          t(abcd_suicide.clust.ksadsY$centers[order(c(4,3,5,2,1)),])))

corrplot::corrplot(-1*cor(t(abcd_suicide.clust$centers),
                          t(abcd_suicide.clust.ksadsP$centers[order(c(4,2,1,5,3)),])))
corrplot::corrplot(-1*cor(t(abcd_suicide.clust$centers),
                          t(abcd_suicide.clust.ksadsY$centers[order(c(4,3,5,2,1)),])))


fviz_nbclust(scale(abcd_trait.bl.clust[abcd_trait.bl.clust$ksads_y_diagnosis==1, c(2:35)]),FUNcluster = kmeans,method = c( "gap_stat"))
fviz_nbclust(scale(abcd_trait.bl.clust[abcd_trait.bl.clust$ksads_p_diagnosis==1, c(2:35)]),FUNcluster = kmeans,method = c( "gap_stat"))


fviz_cluster(abcd_suicide.clust.ksadsP,
             axes = c(1, 2),
             geom = c("point")) +
  ggplot2::theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1")
fviz_cluster(abcd_suicide.clust.ksadsP,
             axes = c(2, 3),
             geom = c("point")) +
  ggplot2::theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1")



fviz_cluster(abcd_suicide.clust.ksadsY,
             axes = c(1, 2),
             geom = c("point")) +
  ggplot2::theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1")
fviz_cluster(abcd_suicide.clust.ksadsY,
             axes = c(2, 3),
             geom = c("point")) +
  ggplot2::theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  scale_colour_brewer(palette = "Set1")



abcd_trait.bl.clust.principal.ksadsP <-
  principal(abcd_trait.bl.clust[abcd_trait.bl.clust$ksads_p_diagnosis==1, 2:35], nfactors = 3)
abcd_trait.bl.clust.principal.ksadsP$values[1:3]/sum(abcd_trait.bl.clust.principal.ksadsP$values)
p3d <- plot3d(
  abcd_trait.bl.clust.principal.ksadsP$scores[, 1],
  abcd_trait.bl.clust.principal.ksadsP$scores[, 2],
  abcd_trait.bl.clust.principal.ksadsP$scores[, 3],
  col = brewer.pal(5, 'Set1')[c(4,2,1,5,3)][abcd_suicide.clust.ksadsP$cluster],
  xlab = "Dim 1 (28.2%)",
  ylab = "Dim 2 (9.9%)",
  zlab = "Dim 3 (7.2%)",
  size = 1,
  type = "s",
  alpha = 1,
  shininess = 80,
  smooth = FALSE,
  lit = FALSE,
  specular = "white",
  emission = "black"
)#type = c("shade","wire","dots")

abcd_trait.bl.clust.principal.ksadsY <-
  principal(abcd_trait.bl.clust[abcd_trait.bl.clust$ksads_y_diagnosis==1, 2:35], nfactors = 3)
abcd_trait.bl.clust.principal.ksadsY$values[1:3]/sum(abcd_trait.bl.clust.principal.ksadsY$values)
p3d <- plot3d(
  abcd_trait.bl.clust.principal.ksadsY$scores[, 1],
  abcd_trait.bl.clust.principal.ksadsY$scores[, 2],
  abcd_trait.bl.clust.principal.ksadsY$scores[, 3],
  col = brewer.pal(6, 'Set1')[c(1:6)][abcd_suicide.clust.ksadsY$cluster],
  xlab = "Dim 1 (28.8%)",
  ylab = "Dim 2 (9.7%)",
  zlab = "Dim 3 (7.2%)",
  size = 1,
  type = "s",
  alpha = 1,
  shininess = 80,
  smooth = FALSE,
  lit = FALSE,
  specular = "white",
  emission = "black"
)#type = c("shade","wire","dots")





fviz_nbclust(scale(abcd_trait.bl.clust[, c(2:35)]),FUNcluster = kmeans,method = c("silhouette"))
fviz_nbclust(scale(abcd_trait.bl.clust[, c(2:35)]),FUNcluster = kmeans,method = c("wss"))
fviz_nbclust(scale(abcd_trait.bl.clust[, c(2:35)]),FUNcluster = kmeans,method = c("gap_stat"))




fviz_nbclust(scale(abcd_trait.bl.clust[1:812, c(2:35)]),FUNcluster = kmeans,method = c("silhouette"))
fviz_nbclust(scale(abcd_trait.bl.clust[813:1624, c(2:35)]),FUNcluster = kmeans,method = c("silhouette"))


