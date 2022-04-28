#!/usr/bin/env Rscript
# get MRS (cleaned and with age group and roi labels)
# NB. visually excluded spectrum are not represented in the shared data.
#     nor is bad data from GPC+Cho, NAA+NAAG, Cr
#
# MRS_GabaGlu.R:MRS_gaba_glu also filtered Gaba/glu outliers. so not using that.
# see MRS_GabaGlu.R for keep_rois labels
source("readMRS.R")
keep_rois <- c(1, 2, 7, 8, 9, 10)

gaba_glu_share <- MRS %>%
   filter(roi %in% keep_rois) %>%
   select(id=ld8, age, sex, roi, label, GMrat, Glu, Glu.SD, GABA, GABA.SD) %>%
   # remove visitdate (age + visitdate is too close to DOB)
   mutate(id=gsub("_2.*", "", id))

write.csv(gaba_glu_share, "out/gaba_glu.csv", row.names=FALSE)
