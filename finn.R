# sliding window - approach 2
rois <- c('DLPFC', 'L DLPFC', 'R DLPFC', 'MPFC', 'ACC', 'Anterior Insula', 'R Anterior Insula', 'L Anterior Insula')

# this creates a structure called 'allSlopes' that has the model coefs for each interation
allSlopes <- c()
for (thisroii in seq(1, length(rois))) {
  thisroi <- rois[thisroii]
  this <- MRS_corr %>% filter(grepl(thisroi, label)) %>% mutate(subj = as.character(ld8)) %>% 
    select(subj, ld8, GABA.Cr, Glu.Cr, label, age, invage, GMrat, sex) %>% filter(complete.cases(.)) %>% arrange(age)
  ages <- this %>% select(age) %>% distinct() %>% arrange(age)
  winSize <- 30
  slopes <- c()
  midAges <- c()
  for (i in seq(1, dim(ages)[1]-winSize)) {
    minAge <- unname(unlist(ages[i,]))
    maxAge <- unname(unlist(ages[i+winSize,]))
    
    if (length(this %>% select(label) %>% distinct()) == 1) {
      #print('Using lm')
      lm.model <- lm(GABA.Cr ~ Glu.Cr + GMrat + sex, data = this %>% filter(age >= minAge & age <= maxAge))
    } else {
      #print('Using lmer')
      lm.model <- lmer(GABA.Cr ~ Glu.Cr + GMrat + sex + label + (1|ld8), data = this %>% filter(age >= minAge & age <= maxAge))
    }
    lm.model.summary <- summary(lm.model)
    
    meanAge <- unname(unlist(this %>% filter(age >= minAge & age <= maxAge) %>% select(age) %>% summarize(midAge = mean(age))))
    
    allSlopes <- rbind(allSlopes, 
                       data.frame(roi = thisroi, 
                                  slope = lm.model.summary$coefficients[2,1], 
                                  se = lm.model.summary$coefficients[2,2], 
                                  r2 = lm.model.summary$r.squared,
                                  adjr2 = lm.model.summary$adj.r.squared,
                                  minAge, maxAge, meanAge))
  }
  
}

# plotting is then just
lunaize(ggplot(data = allSlopes %>% filter(roi %in% c('L DLPFC','R DLPFC','L Anterior Insula','R Anterior Insula','MPFC','ACC')), aes(x = meanAge, y= slope)) +
          geom_point(alpha=.5, size=.1) + geom_errorbarh(aes(xmin = minAge,xmax = maxAge), alpha=0.1) + geom_errorbar(aes(ymin = slope - se, ymax = slope+se), alpha=0.1) +
          stat_smooth(method='loess', span=1.5, se=F, size=2) + facet_wrap(. ~ roi) + labs(y='GABA-Glu Association', x='Age'))

# also, i can't remember if this is in maria's, but the first thing i do is create MRS_corr--
MRS_corr <- MRS_glu %>%
  filter(GABA.SD <= 20) %>%
  group_by(roi) %>%
  mutate(zscore=scale(GABA.Cr, center=T, scale=T)) %>%
  filter(abs(zscore) <= z_thres) %>% ungroup %>%
  filter(visitnum == 1) %>%
  mutate(agegrp = cut(age,
                      breaks=c(0,16,22,Inf),
                      labels=c("10-16","17-22", "23-30")))
