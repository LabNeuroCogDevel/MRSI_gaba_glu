# depends on MRS from readMRS.R
source("readMRS.R")
source("finn.R")
pacman::p_load(dplyr,ggplot2,viridis,RColorBrewer, cowplot)
select <- dplyr::select
z_thres = 2
keep_rois <- c(1,2,7,8,9,10)   # or keep_rois <- names(region_look)

MRS_glu <- MRS %>%
  mutate(agegrp = cut(age, breaks=c(0,16,22,Inf),
                      labels=c("10-16","17-22", "23-30")),
         region=unlist(region_look[as.character(roi)])) %>%
  filter(Glu.SD <=20)  %>%
  group_by(roi) %>%
  mutate(zscore.glu=scale(Glu.Cr, center=T, scale=T)) %>%
  filter(abs(zscore.glu) <= z_thres)

MRS_corr <- MRS_glu %>%
  filter(GABA.SD <= 20, !is.na(GMrat)) %>%
  group_by(roi) %>%
  mutate(zscore.gaba=scale(GABA.Cr, center=T, scale=T)) %>%
  filter(abs(zscore.gaba) <= z_thres)

# moving average model
allSlopes <- ma_gaba_glu_slope(MRS_corr)

# add residuals
glu.resids <- lm(data=MRS_corr, Glu.Cr ~ sex + GMrat, na.action="na.exclude")
MRS_corr$glu.resids <- residuals(glu.resids)
gaba.resids <- lm(data=MRS_corr, GABA.Cr ~ sex + GMrat, na.action="na.exclude")
MRS_corr$gaba.resids <- residuals(gaba.resids)

region_avg <- MRS_corr %>%
  filter(roi %in% keep_rois) %>%          # remove any rois not in region_loopup
  select(ld8, agegrp,glu.resids,gaba.resids, label, region, age) %>%
  filter(!is.na(region)) %>%
  group_by(label,region,ld8, agegrp,age) %>%
  summarise_at(vars(glu.resids, gaba.resids), mean, na.rm=T)


cor_vals <- region_avg %>%
  group_by(label,region,agegrp) %>%
  summarise(GabaGlu_r=cor(glu.resids, gaba.resids), 
            GabaGlu_p=cor.test(glu.resids, gaba.resids, method=c("pearson"), use = "complete.obs")$p.value, 
            GabaGlu_lb=cor.test(glu.resids, gaba.resids, method=c("pearson"), use = "complete.obs", conf.level=.95)$conf.int[1], 
            GabaGlu_ub=cor.test(glu.resids, gaba.resids, method=c("pearson"), use = "complete.obs", conf.level=.95)$conf.int[2], 
            n=n(),
            agegrp_mean=mean(age)) %>%
    mutate(agegrp_mid=c(14,20,26)[agegrp])


    #ggplot()+geom_bar(data=cor_vals, aes(x=agegrp_mid, y=GabaGlu_r, fill=agegrp), stat="identity", width=5.5) + facet_wrap(~region)

corr_bar_and_line <-
 allSlopes %>% mutate(label=roi) %>% filter(label %in% unique(cor_vals$label)) %>%
 ggplot() +
    aes(x = meanAge, y= slope) +
    geom_bar(data=cor_vals,
             aes(x=agegrp_mid, y=GabaGlu_r, fill=agegrp),
             stat="identity", position="dodge",width=5.5) +
    geom_point(alpha=.5, size=.1) +
    geom_errorbarh(aes(xmin = minAge,xmax = maxAge), alpha=0.1) +
    geom_errorbar(aes(ymin = slope - se, ymax = slope+se), alpha=0.1) +
    stat_smooth(method='loess', span=1.5, se=F, size=2) +
    facet_wrap(. ~ label) +
    theme_cowplot() +
    labs(x='Age',
         y='Glu GABA Corr (r)',
         fill = "Age Group")
print(corr_bar_and_line)

