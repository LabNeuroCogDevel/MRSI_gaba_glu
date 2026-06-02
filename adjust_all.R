#!/usr/bin/env Rscript
source('mrsi_funcs.R') # mrsi_add_cols(), mrsi_metqc(), mrs_wide_to_long_cleaned()
source('res_with_age.R')
library(tidyr) # nest, unnest

# 20230531 - intially want GABA and Glu. but easy to add more
# 20230614 - remove MM20. regexp also removes Glu.Gln.Cr and GPC.Cho.Cr
# 20230630 - add ML and GSH
# 20230914 - add Taurine
# 20241108 - move more functions into res_with_age.R
# 20260520 - add GPC
# 20260602 - get raw Cr

mets_keep <- c("GABA","Glu","Gln","Cho","Glc", "NAA", "mI","GSH", "Tau", "NAAG", "GPC")
mets_regex <- colname_sd_or_cr(mets_keep)

# read in raw data. clean. and reduce columns to just those we care about
all_mrs <- read.csv('13MP20200207_LCMv2fixidx.csv')
mrs_qc <- all_mrs %>%
   filter(!failqc) %>%
   mrsi_metqc() %>% mrsi_add_cols()
mrs <- mrs_qc %>%
   select(ld8, region, age, GMrat, dateNumeric,matches(mets_regex))

# 20260602 - raw Cr values (for reviewer)
cr_raw <- mrs_qc |> select(ld8, region, Cr_raw=Cr, Cr_raw.SD=Cr.SD) 
write.csv(cr_raw,"out/cr_raw.csv", quote=F, row.names=F)

# track long format and thresholded; nrow(mrs_long) == 34056;
# will be input to residual modeling
# cols: "ld8", "region", "age", "GMrat", "dateNumeric", "met", "Cr", "SD", "biregion", "met_crz"
mrs_long <- mrs %>% mrs_wide_to_long_cleaned(z_thres=3)
write.csv(mrs_long, "out/long_thres.csv", quote=F, row.names=F)

# save copy of long with Cr appened. Likely don't want this but the modeled version below
mrs_long_wcr <- merge(mrs_long,cr_raw, by=c("ld8","region"), all.x=T)
write.csv(mrs_long_wcr, "out/long_thres_cr.csv", quote=F, row.names=F)

# apply res_with_age to each group
chunked_by_met_region <- mrs_long %>%
   ungroup() %>%
   nest(.by=c("met", "biregion"), .key="metdata")
mrs_long_adj <-
   chunked_by_met_region %>%
   mutate(metdata=lapply(metdata, function(d) res_with_age(d, met='Cr',return_df=TRUE))) %>%
   unnest(cols=c("metdata"))

write.csv(mrs_long_adj, "out/gamadj_long_wfmod.csv", quote=F, row.names=F)

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

# 20260602 - again with exra Cr columns (2 for each region)
# mrs_wide_adj_gm <- read.csv("out/gamadj_wide.csv")
cr_raw_wide <- cr_raw |> pivot_wider(id_cols='ld8',names_from=c('region'), values_from=c('Cr_raw','Cr_raw.SD'))
mrs_wide_adj_gm_wcr <- merge(mrs_wide_adj_gm, cr_raw_wide, by=c("ld8"), all.x=T)
write.csv(mrs_wide_adj_gm_wcr, "out/gamadj_wide_cr.csv", quote=F,row.names=F)
