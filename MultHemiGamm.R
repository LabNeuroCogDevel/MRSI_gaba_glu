library(tidyverse)
library(gamm4)
library(readxl)
library(ggplot2)
library(lubridate)
library(mgcViz)
library(tinytex)
library(ggpubr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(gamm4)
library(dplyr)
library(checkmate) # for expect_subset dataframe names test


# multiple hemis
met_res_for_roi <- function(MRSI_input, this_roi, met_name, resid_fml=met ~ s(age, k=3), minage=9, maxage=31) {
  
  # have columns we need. sanity checking for better error message. not needed for analysis
  model_vars <- grep(invert=T,'met',all.vars(resid_fml), value=T) # need whats in model, but we'll make met
  expect_subset(c("roi", "dateNumeric", met_name, model_vars), names(MRSI_input))
  
  roi_met <- MRSI_input %>% 
    filter(roi %in% this_roi, age >=minage, age <=maxage) %>%
    mutate(met = .data[[met_name]],  # map whatever metabolite we're using (met_name) onto the 'met' column
           gamresids  = gam(resid_fml, data=cur_data()) %>% residuals,
           dateAdjust = gam(gamresids ~ s(dateNumeric, k=3), data=cur_data()) %>% predict,
           met_adj = met - dateAdjust)
  
  return(roi_met)
}

plot_met_adjusted <- function(roi_met) {
  # have columns we need
  expect_subset(c("age", "gamresids", "dateAdjust", "met", "met_adj", "visitnum", "subjID", "label"),
                names(roi_met))
  qassert(roi_met$met_adj, "n+") # nonzero numeric vector in met_adj
  
  cowplot::plot_grid(  
    ggplot(roi_met) + aes(x=age, y=met_adj, color=label, group=paste(label, subjID)) + geom_point() + stat_smooth(method='loess', se=F,aes(group=label)) + geom_line() + theme_classic(base_size=20) +  guides(color="none"),
    ggplot(roi_met) + aes(x=age, y=met, color=label, group=paste(label, subjID)) + geom_point() + stat_smooth(method='loess', se=F,aes(group=label)) + geom_line() + theme_classic(base_size=20)
  )
}


met_lm_stats <- function(roi_met, additional_regressors=c()) {
  expect_subset(c("invage", "met_adj", "GMrat", "sex", "subjID", additional_regressors), names(roi_met))
  fml_str <- "met_adj ~ invage + GMrat + sex + (1|subjID)"
  fml_str <- paste(collapse="+",sep="+", fml_str, additional_regressors)
  met_lm <- lmer(data=roi_met, as.formula(fml_str))
  summary(met_lm)
}


adj_df <- met_res_for_roi(MRS_glu, c(9,10), "Glu.Cr")
plot_met_adjusted(adj_df)
met_lm_stats(adj_df, c("label"))
ggsave("DLPFC_Glu.pdf", dpi = 300, height = 6, width = 12)


