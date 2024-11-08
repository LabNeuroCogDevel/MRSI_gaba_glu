#!/usr/bin/env Rscript
# save gam models to apply to external datasets

# 20241031WF - init
library(dplyr); library(tidyr);

source('mrsi_funcs.R') # mrsi_add_cols(), mrsi_metqc(), mrs_wide_to_long_cleaned()
source('res_with_age.R') 

mets_keep <- c("GABA","Glu","Gln","Cho","Glc", "NAA", "mI","GSH", "Tau")
# 20241108 - additional metabolites
mets_keep <- c(mets_keep, "NAAG","NAA.NAAG","GPC.Cho","Glu.Gln")
mets_regex <- colname_sd_or_cr(mets_keep)

# from adjust_all.R
# read in raw data. clean. and reduce columns to just those we care about
mrs <- read.csv('13MP20200207_LCMv2fixidx.csv') %>%
   filter(!failqc) %>%
   mrsi_metqc() %>% mrsi_add_cols() %>%
   select(ld8, region, age, GMrat, dateNumeric,matches(mets_regex))

mrs_long <- mrs %>% mrs_wide_to_long_cleaned()
chunked_by_met_region <- mrs_long %>%
   ungroup() %>%
   nest(.by=c("met", "biregion"), .key="metdata")

# all data
gam_models_all <- chunked_by_met_region %>%
   mutate(model=lapply(metdata, \(d) mrsi_model_date(d, met_name='Cr', include.na=FALSE)))
# only with 18yo+
gam_models_18 <- chunked_by_met_region %>%
   mutate(model=lapply(metdata, \(d) mrsi_model_date(d, met_name='Cr', include.na=FALSE, min_age=18)))

save(list=c("gam_models_all","gam_models_18"), file="mgcv-gam_tibble.Rdata")

## compare
data_gaba_thal <- mrs_long %>% filter(met=="GABA",  biregion =='Thalamus') %>% na.omit()

gam_gaba_thal_all <- gam_models_all %>% filter(met=="GABA", biregion=='Thalamus') %>% pull(model) %>% `[[`(1)
gam_gaba_thal_18  <- gam_models_18  %>% filter(met=="GABA", biregion=='Thalamus') %>% pull(model) %>% `[[`(1)

vals <- data.frame(
           gt_gm_all=adj_by_model(gam_gaba_thal_all, data_gaba_thal, center_on=c("dateNumeric","GMrat")),
           gt_nogm_all=adj_by_model(gam_gaba_thal_all, data_gaba_thal, center_on=c("dateNumeric")),
           gt_gm_18=adj_by_model(gam_gaba_thal_18, data_gaba_thal, center_on=c("dateNumeric","GMrat")),
           gt_nogm_18=adj_by_model(gam_gaba_thal_18, data_gaba_thal, center_on=c("dateNumeric"))) %>% 
   cbind(data_gaba_thal %>% select(age,GMrat,dateNumeric,met,biregion))
vals %>% pivot_longer(cols=matches('^gt')) %>% ggplot() + aes(x=age,y=value)

## single model
model_gaba_thal <- mrsi_model_date(data_gaba_thal, met='Cr')
viz <- mgcViz::getViz(model_gaba_thal, allTerms = T)
plot(viz)

