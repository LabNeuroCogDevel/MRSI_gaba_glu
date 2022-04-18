#!/usr/bin/Rscript

# also, i can't remember if this is in maria's, but the first thing i do is create MRS_corr--
# MRS_corr <- MRS_glu %>%
#   filter(GABA.SD <= 20) %>%
#   group_by(roi) %>%
#   mutate(zscore=scale(GABA.Cr, center=T, scale=T)) %>%
#   filter(abs(zscore) <= z_thres) %>% ungroup %>%
#   filter(visitnum == 1) %>%
#   mutate(agegrp = cut(age,
#                       breaks=c(0,16,22,Inf),
#                       labels=c("10-16","17-22", "23-30")))

ma_gaba_glu_slope <- function(MRS_corr) {
  # sliding window - approach 2
  #rois <- c('DLPFC', 'L DLPFC', 'R DLPFC',
  #          'MPFC', 'ACC',
  #          'Anterior Insula',
  #          'R Anterior Insula', 'L Anterior Insula')
  rois <- unique(MRS_corr$label)
  
  # this creates a structure called 'allSlopes' that has the model coefs for each interation
  allSlopes <- c()
  for (thisroii in seq(1, length(rois))) {
    thisroi <- rois[thisroii]
    this <- MRS_corr %>%
        ungroup %>%
        filter(grepl(thisroi, label)) %>%
        mutate(subj = as.character(ld8)) %>% 
        select(subj, ld8, GABA.Cr, Glu.Cr,
               label, age, invage, GMrat, sex) %>%
        filter(complete.cases(.)) %>%
        arrange(age)
    ages <- this$age %>% unique %>% sort
    winSize <- 30
    slopes <- c()
    midAges <- c()
    # run models at each age in winSize bins
    for (i in seq(1, length(ages)-winSize)) {
      minAge <- ages[i]
      maxAge <- ages[i+winSize]
      d_agerange <-this %>% filter(age >= minAge & age <= maxAge)
  
      # labels+ld8 combos is the actual number we care about?
  
      #n_labels <- this %>% select(label) %>% distinct() %>% nrow
      n_labels <- this$label %>% unique %>% length
      
      if (n_labels == 1) {
        #print('Using lm')
        lm.model <- lm(GABA.Cr ~ Glu.Cr + GMrat + sex,
                       data = d_agerange)
      } else {
        #print('Using lmer')
        lm.model <- lmer(GABA.Cr ~ Glu.Cr + GMrat + sex + label + (1|ld8),
                         data = d_agerange)
      }
      lm.model.summary <- summary(lm.model)
      
      meanAge <- mean(d_agerange$age)
      
      allSlopes <- rbind(allSlopes, 
                         data.frame(roi = thisroi, 
                                    slope = lm.model.summary$coefficients[2,1], 
                                    se = lm.model.summary$coefficients[2,2], 
                                    r2 = lm.model.summary$r.squared,
                                    adjr2 = lm.model.summary$adj.r.squared,
                                    minAge, maxAge, meanAge, n=nrow(d_agerange)))
    }
  }
 return(allSlopes)
}


# plotting is then just
ma_slop_plot<-function(allSlopes){

  # prev restricted to rois:
  #filter(roi %in% c('L DLPFC','R DLPFC','L Anterior Insula','R Anterior Insula','MPFC','ACC'))) 

  ggplot(allSlopes) +
    aes(x = meanAge, y= slope) +
    geom_point(alpha=.5, size=.1) +
    geom_errorbarh(aes(xmin = minAge,xmax = maxAge), alpha=0.1) +
    geom_errorbar(aes(ymin = slope - se, ymax = slope+se), alpha=0.1) +
    stat_smooth(method='loess', span=1.5, se=F, size=2) +
    facet_wrap(. ~ roi) +
    labs(y='GABA-Glu Association', x='Age')
}
