#!/usr/bin/env Rscript
source("MRS_GabaGlu.R") # MRS_GabaGlu MRS_GabaGlu_r 
source("ma_gaba_glu_slope.R") # function: ma_gaba_glu_slope.R
pacman::p_load(dplyr,ggplot2,viridis,RColorBrewer, cowplot)
select <- dplyr::select # might not be necessary, but keep it just in case

# moving average model
allSlopes <- ma_gaba_glu_slope(MRS_GabaGlu)

    #ggplot()+geom_bar(data=cor_vals, aes(x=agegrp_mid, y=GabaGlu_r, fill=agegrp), stat="identity", width=5.5) + facet_wrap(~region)

corr_bar_and_line <-
 allSlopes %>% mutate(label=roi) %>% filter(label %in% unique(MRS_GabaGlu_r$label)) %>%
 ggplot() +
    aes(x = meanAge, y= slope) +
    geom_bar(data=MRS_GabaGlu_r,
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
