#!/usr/bin/env Rscript
source('mrsi_funcs.R')
source('res_with_age.R')
library(tidyr) # nest, unnest
zscore_abs <- function(x) abs(scale(x,center=T,scale=T)[,1])

# 20230531 - intially want GABA and Glu. but easy to add more
mets_keep <- c("GABA","Glu","Gln","Cho","Glc", "MM20")

# select metabolites based on column name using dplyr::match w/ regular expression
mets_patt <- paste0(collapse="|", mets_keep) # MM20|GABA|Glu|Gln|Cho|Glc
mets_regex <- glue::glue('({mets_patt})\\.(Cr|SD)') # Cr ratio and SD

# read in raw data. clean. and reduce columns to just those we care about
mrs <- read.csv('13MP20200207_LCMv2fixidx.csv') %>%
   mrsi_metqc() %>% mrsi_add_cols() %>%
   select(ld8, region, age, GMrat, dateNumeric,matches(mets_regex))


# longer with new cols 'met', 'Cr', and 'SD' (reshape to remove per metabolite columns)
z_thres <- 3
mrs_long <- longer_met_CrSD(mrs) %>%
   # collapse across hemispheres
   # conviently R,L are always hemisphere. others: MPFC ACC 
   mutate(biregion=gsub('^(R|L)','',region)) %>%
   # work on each region X met separetly
   # TODO: should thresholding be just by met (not also by region?)
    #      OR should group by biregion?
   group_by(met, region) %>% 
   mutate(met_crz=zscore_abs(Cr),
          met_sdz=zscore_abs(SD),
          # TODO: also threshold on SD? or zscore of SD?
          Cr=ifelse(met_crz>=z_thres, NA, Cr))
write.csv(mrs_long, "out/long_thres.csv", quote=F, row.names=F)

# apply res_with_age to each group
mrs_long_adj <- mrs_long %>% ungroup() %>%
   nest(.by=c("met", "biregion"), .key="metdata") %>%
   mutate(metdata=lapply(metdata, function(d) res_with_age(d, met='Cr',return_df=TRUE))) %>%
   unnest(cols=c("metdata"))

write.csv(mrs_long_adj, "out/gamadj_long.csv", quote=F, row.names=F)

# reshape to look like input again: column per region+met+measure
mrs_wide_adj <- mrs_long_adj %>%
    select(ld8,region=biregion,met,Cr,gamadj=Cr_gamadj,SD) %>%
    pivot_wider(id_cols = c("ld8"),
                names_from = c("region", "met"),
                values_from = c("gamadj","Cr","SD"),
                # NB. might have value per hemispere. simple mean to collapse
                values_fn = mean)

write.csv(mrs_wide_adj, "out/gamadj_wide.csv", quote=F,row.names=F)
