#!/usr/bin/env Rscript
source("./210517_MRSPaper.R")

# https://r-graph-gallery.com/line-chart-dual-Y-axis-ggplot2.html
MRS_glu <- MRS %>% filter(Glu.SD <=20)

MRS_glu <- MRS %>% filter(Glu.SD <=20)
MRS_glu <- MRS_glu %>%
  group_by(roi) %>%
  mutate(zscore=scale(Glu.Cr, center=T, scale=T)) %>%
  filter(abs(zscore) <= z_thres) %>% ungroup 

MRS_corr <- MRS_glu %>%
  filter(GABA.SD <= 20) %>%
  group_by(roi) %>%
  mutate(zscore=scale(GABA.Cr, center=T, scale=T)) %>%
  filter(abs(zscore) <= z_thres) %>% ungroup %>%
  mutate(agegrp = cut(age,
                      breaks=c(0,16,22,Inf),
                      labels=c("10-16","17-22", "23-30")))
