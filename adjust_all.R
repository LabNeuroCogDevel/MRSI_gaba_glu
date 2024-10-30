#!/usr/bin/env Rscript
source('mrsi_funcs.R') # mrsi_add_cols
source('res_with_age.R')
library(tidyr) # nest, unnest
zscore_abs <- function(x) abs(scale(x,center=T,scale=T)[,1])

# 20230531 - intially want GABA and Glu. but easy to add more
# 20230614 - remove MM20. regexp also removes Glu.Gln.Cr and GPC.Cho.Cr
# 20230630 - add ML and GSH
# 20230914 - add Taurine
mets_keep <- c("GABA","Glu","Gln","Cho","Glc", "NAA", "mI","GSH", "Tau")

# select metabolites based on column name using dplyr::match w/ regular expression
mets_patt <- paste0(collapse="|", mets_keep) # GABA|Glu|Gln|Cho|Glc
mets_regex <- glue::glue('^({mets_patt})\\.(Cr|SD)') # Cr ratio and SD

# read in raw data. clean. and reduce columns to just those we care about
mrs <- read.csv('13MP20200207_LCMv2fixidx.csv') %>%
   filter(!failqc) %>%
   mrsi_metqc() %>% mrsi_add_cols() %>%
   select(ld8, region, age, GMrat, dateNumeric,matches(mets_regex))


# longer with new cols 'met', 'Cr', and 'SD' (reshape to remove per metabolite columns)
z_thres <- 3
# added 20230614 (to match previous analysis)
# 20 for Glu and GABA. guess for others
# TODO: VALIDATE Glc and Cho SD thres
sd_thres <- function(x)
   case_when(
             grepl('Glu|GABA',x) ~ 20,
             grepl('Cho',x)      ~ 40, # TODO: check
             grepl('Glc',x)      ~ 120,# TODO: check
             .default=20)              # taurine is 20, just not explicitly
mrs_long <- longer_met_CrSD(mrs) %>%
   # collapse across hemispheres
   # conveniently R,L are always hemisphere. others: MPFC ACC 
   mutate(biregion=gsub('^(R|L)','',region)) %>%
   # work on each region X met separately
   group_by(met, region) %>% 
   #NB. 'Cr' here is Cr ratio for a given metabolite
   #    sd_thres is 20 for gaba and glu
   mutate(Cr=ifelse(SD>sd_thres(met), NA, Cr),
          met_crz=zscore_abs(Cr),
          Cr=ifelse(met_crz>=z_thres, NA, Cr))
write.csv(mrs_long, "out/long_thres.csv", quote=F, row.names=F)

# apply res_with_age to each group
chunked_by_met_region <- mrs_long %>%
   ungroup() %>%
   nest(.by=c("met", "biregion"), .key="metdata")
mrs_long_adj <-
   chunked_by_met_region %>%
   mutate(metdata=lapply(metdata, function(d) res_with_age(d, met='Cr',return_df=TRUE))) %>%
   unnest(cols=c("metdata"))

# 20241023 - added save_model. rerun above but save out model instead
gam_models <- chunked_by_met_region %>%
   mutate(model=lapply(metdata, \(d) res_with_age(d, met='Cr',return_model=TRUE)))
save(list=c("gam_models"), file="mgcv-gam_tibble.Rdata")

# 20241030 - issue with gamViz?
#setwd("/Volumes/Hera/Projects/7TBrainMech/scripts/mri/MRSI_roi/gam_adjust")
#source('res_with_age.R')
#mrs_long <- read.csv("out/long_thres.csv")
#data_gaba_acc <- mrs_long %>% filter(met=="GABA",  biregion =='ACC') %>% na.omit()
#model_gaba_acc <- res_with_age(data_gaba_acc, met='Cr',return_model=TRUE)
#viz <- mgcViz::getViz(model_gaba_acc, allTerms = T)
#plot(viz) # Error ... replacement has 268 rows, data has 256

write.csv(mrs_long_adj, "out/gamadj_long.csv", quote=F, row.names=F)

# NB! In both bilateral average and idv columns,
#     resids are from model that includes both hemis in input
#     values (rows) for both L and R are input to model => ... + region
#     see res_with_age.R

# reshape to look like input again: column per region+met+measure
mrs_wide_adj_bilat <- mrs_long_adj %>%
    select(ld8,region=biregion,met,Cr,gamadj=Cr_gamadj,SD) %>%
    pivot_wider(id_cols = c("ld8"),
                names_from = c("region", "met"),
                values_from = c("gamadj","Cr","SD"),
                names_glue = "{region}_{met}_{.value}",
                # NB. might have value per hemispere. simple mean to collapse
                values_fn = mean)

# keep lateral regions
mrs_wide_adj_lat <- mrs_long_adj %>%
    select(ld8,region,met,Cr,gamadj=Cr_gamadj,SD) %>%
    filter(grepl('^[LR]',region)) %>% 
    pivot_wider(id_cols = c("ld8"),
                names_from = c("region", "met"),
                values_from = c("gamadj","Cr","SD"),
                names_glue = "{region}_{met}_{.value}") %>%
    # 20230614: add Cr and SD, not just gamadj
    select(ld8,matches('gamadj|Cr|SD'))

# GMrat is not per metabolite. so we'll do the collapsing separately
# we want bilateral and averaged GMrats
mrs_wide_gmrat_bi <- mrs_long_adj %>%
    select(ld8,region=biregion, GMrat) %>%
    pivot_wider(id_cols = c("ld8"),
                names_from = c("region"),
                values_from = c("GMrat"),
                # match region_met_.value from above
                names_glue = "{region}_all_GMrat",
                # NB. might have value per hemispere. simple mean to collapse
                values_fn = function(x) mean(x,na.rm=T))

mrs_wide_gmrat_hemi <- mrs_long_adj %>%
    select(ld8,region, GMrat) %>%
    pivot_wider(id_cols = c("ld8"),
                names_from = c("region"),
                values_from = c("GMrat"),
                # match region_met_.value from above
                names_glue = "{region}_all_GMrat",
                # NB. might have value per hemispere. simple mean to collapse
                values_fn = function(x) mean(x,na.rm=T))
bilat_names <- setdiff(names(mrs_wide_gmrat_bi), names(mrs_wide_gmrat_hemi))
mrs_wide_gmrat <- merge(mrs_wide_gmrat_hemi, mrs_wide_gmrat_bi[c("ld8",bilat_names)], all=T)

mrs_wide_adj_gm <- merge(mrs_wide_adj_lat, mrs_wide_gmrat, by="ld8") %>%
   merge(mrs_wide_adj_bilat, by="ld8")

write.csv(mrs_wide_adj_gm, "out/gamadj_wide.csv", quote=F,row.names=F)
