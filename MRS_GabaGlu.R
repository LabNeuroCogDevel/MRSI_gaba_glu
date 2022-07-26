# depends on MRS from readMRS.R
# exports
#  MRS_GabaGlu  -- used by ma_gaba_glu_slope
#  MRS_GabaGlu_r -- correlation with residisualized
source("readMRS.R") # get MRS (cleaned and with age group and roi labels)
select <- dplyr::select
z_thres = 2
keep_rois <- c(1,2,7,8,9,10)

MRS_glu <- MRS %>%
  mutate(agegrp = cut(age, breaks=c(0,16,22,Inf),
                      labels=c("10-16","17-22", "23-30")),
         region=unlist(region_look[as.character(roi)])) %>%
  filter(Glu.SD <=20)  %>%
  group_by(roi) %>%
  mutate(zscore.glu=scale(Glu.Cr, center=T, scale=T)) %>%
  filter(abs(zscore.glu) <= z_thres)

MRS_GabaGlu <- MRS_glu %>%
  filter(GABA.SD <= 20, !is.na(GMrat)) %>%
  group_by(roi) %>%
  mutate(zscore.gaba=scale(GABA.Cr, center=T, scale=T)) %>%
  filter(abs(zscore.gaba) <= z_thres)

# add residuals
glu.resids <- lm(data=MRS_GabaGlu, Glu.Cr ~ sex + GMrat, na.action="na.exclude")
MRS_GabaGlu$glu.resids <- residuals(glu.resids)
gaba.resids <- lm(data=MRS_GabaGlu, GABA.Cr ~ sex + GMrat, na.action="na.exclude")
MRS_GabaGlu$gaba.resids <- residuals(gaba.resids)

# average resids over identicial label/region
# NB. label/region are now all unique within visit.
#     prev might have had region names exclude L/R portion of label
#     to collapse across hemisphere
#     now summarise does nothing: mean of 1 element is identical
# TODO: can probably drop 'label' column from all group_by
#       column is used in plot_GabaGlu_slope_and_bar.R but as an alias for 'roi'
gabaGluRes_region_avg <- MRS_GabaGlu %>%
  filter(roi %in% keep_rois) %>%          # remove any rois not in region_loopup
  select(ld8, agegrp,glu.resids,gaba.resids, label, region, age) %>%
  filter(!is.na(region)) %>%
  group_by(label,region,ld8, agegrp,age) %>%
  summarise_at(vars(glu.resids, gaba.resids), mean, na.rm=T)

# get correlations and confidence intervals
# now (20220726) label==region. holdover from collapsing hemisphers?
# see readMRS.R:region_look
MRS_GabaGlu_r <- gabaGluRes_region_avg %>%
  group_by(label,region,agegrp) %>%
  summarise(GabaGlu_r=cor(glu.resids, gaba.resids), 
            GabaGlu_p=cor.test(glu.resids, gaba.resids, method=c("pearson"), use = "complete.obs")$p.value, 
            GabaGlu_lb=cor.test(glu.resids, gaba.resids, method=c("pearson"), use = "complete.obs", conf.level=.95)$conf.int[1], 
            GabaGlu_ub=cor.test(glu.resids, gaba.resids, method=c("pearson"), use = "complete.obs", conf.level=.95)$conf.int[2], 
            n=n(),
            agegrp_mean=mean(age)) %>%
    mutate(agegrp_mid=c(14,20,26)[agegrp])

