#!/usr/bin/env Rscript
library(ggplot2)
library(tidyr) # nest, unnest
library(dplyr)


# show number of excluded by thresholding for 
mrs_long <- read.csv("out/long_thres.csv")
g_cnt <- mrs_long %>% group_by(met,region) %>%
    summarise(n_discard=length(which(is.na(Cr))),
              n_kept=n()-n_discard, 
              #prct_keep=100*(n()-nna)/n())
              ) %>%
    tidyr::pivot_longer(cols=c('n_discard','n_kept'))

p <- ggplot(g_cnt) +
    aes(x=region, fill=name, y=value) +
    facet_wrap(~met) +
    geom_bar(stat='identity', position='stack') +
    see::theme_modern(axis.text.angle = 90)
ggsave(p, file="imgs/thresholding_cnt_met_region.png", width=8.08, height=8.63)


# adjusted compared to raw. relationship looks linear (yay)
mrs_long_adj <- read.csv("out/gamadj_long.csv")
p_adj <- ggplot(mrs_long_adj) + 
  aes(x=Cr, y=Cr_gamadj)+
  facet_wrap(met~biregion, scales = "free") +
  geom_point(aes(size=-pmin(met_crz,3),
                 color=-pmin(SD,10)))+ 
  geom_smooth(method='lm') +
  see::theme_modern() +
  scale_size(range=c(.25,1)) +
  scale_y_continuous(
    limits = ~ c(min(.x), ceiling(max(.x))),
    breaks = ~ .x[2],
    expand = c(0, 0))
ggsave(p_adj, file="imgs/gam_adjusted_Vs_Cr.png", width=16.3, height=8.65)


