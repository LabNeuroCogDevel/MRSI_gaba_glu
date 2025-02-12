library(tidyr)
library(dplyr)
library(ggplot2)
library(mgcv)
library(gratia)
library(mgcViz)
library(patchwork)
library(gamm4)
library(readxl)

# Load models
load('~/scratch/deepak/mgcv-gam_tibble.Rdata')
gam_models <- gam_models_18

# Load SZ MRSI data
szmet_orig <- readxl::read_xlsx('caud-thal_MRSI data_10.11.24.xlsx', sheet = 2) %>%
  rename(scan_date = `scan date`) %>% select(-GM891) %>% mutate(source = 'orig')
szmet <- rbind(szmet_orig,
      readxl::read_xlsx('7T.DARES.baselines.xlsx', sheet = 1) %>%
        rename(region = roi) %>%
        mutate(timepoint = 1) %>%
        mutate(source = 'DARES.baselines') %>%
        select(names(szmet_orig)),
      readxl::read_xlsx('BL.FU.mrsi.xlsx', sheet = 4) %>%
        rename(scan_date = `scan date`) %>%
        mutate(source = 'BL.FU') %>%
        select(names(szmet_orig)))

# Match ROI names, create date value
szmet_new <- szmet %>% 
  separate(region, c('hemi','roi'), ' ', convert=TRUE) %>%
  mutate(roi = ifelse(roi == 'caudate', 'Caudate', roi)) %>%
  mutate(roi = ifelse(roi == 'thalamus', 'Thalamus', roi)) %>%
  mutate(dateNumeric = as.numeric(as.POSIXct(scan_date, format="%Y-%m-%d")),
         GMrat = GM,
         hemi = ifelse(hemi=='left', 'L', ifelse(hemi=='right', 'R', NA))) %>%
  mutate(region = paste0(hemi, roi))

# Grab model based on ROI & metabolite
gam.model <- gam_models %>% filter(met == 'GABA', biregion == 'Thalamus') %>% pull(model)
metdata <- gam_models %>% filter(met == 'GABA', biregion == 'Thalamus') %>% pull(metdata)

# Extract necessary data from sz data frame
this.met <- szmet_new %>%
  select(RECID, timepoint, age, hemi, roi, dateNumeric, GMrat, GABA.Cre = `GABA/Cre`) %>%
  filter(roi == 'Thalamus')

## first, get residual (difference from expectation given date, GMrat, age)
yhat <- unname(predict(gam.model[[1]], this.met))
e <- this.met$GABA.Cre - yhat

## now, predict @ mean date (but with real age & GMrat) and add back in residual
this.met$GABA.Cre.adj <- unname(predict(gam.model[[1]], this.met %>% mutate(dateNumeric = mean(metdata[[1]]$dateNumeric, na.rm=T)))) + e

ggplot(data = this.met, aes(x = GABA.Cre, y = GABA.Cre.adj, color=as.factor(hemi), shape=as.factor(timepoint))) +
  geom_point() +
  geom_path(aes(group = interaction(RECID, hemi))) +
  stat_smooth(method='lm', alpha=0.2) +
  geom_abline(slope=1, intercept = 0, linetype=2) +
  theme(legend.position=c(.1, .8)) +
  theme_bw()


# Loop all
names(szmet_new %>% select(-contains('%') & contains('/Cre')))
gam_models %>% select(met) %>% distinct()

sz_mets <- c('GABA/Cre', 'Glu/Cre', 'Gln/Cre','Glc/Cre', 'GSH/Cre','mI/Cre','NAA/Cre', 'Tau/Cre', 'Glu+Gln', 'NAAG', 'GPC/Cre', 'GPC+Cho/Cre')
mod_mets <- c('GABA',    'Glu',     'Gln',    'Glc',     'GSH',    'mI',    'NAA',     'Tau',     'Glu.Gln', 'NAAG', 'GPC',     'GPC.Cho')

#Glu+Gln, NAAG, GPC, GPC+CHO

adj.df <- c()
for (thisroi in unique(szmet_new$roi)) {
  
  for (meti in seq(1, length(sz_mets))) {
    this_mod_met <- mod_mets[meti]
    this_sz_met <- sz_mets[meti]
    
    # Grab model based on ROI & metabolite
    gam.model <- gam_models %>% filter(met == this_mod_met, biregion == thisroi) %>% pull(model)
    metdata <- gam_models %>% filter(met == this_mod_met, biregion == thisroi) %>% pull(metdata)
    
    # Extract necessary data from sz data frame
    this.met <- szmet_new %>%
      select(RECID, timepoint, age, region, hemi, roi, dateNumeric, GMrat, source, met = !!this_sz_met) %>%
      filter(roi == thisroi)
    
    ## first, get residual (difference from expectation given date, GMrat, age)
    yhat <- unname(predict(gam.model[[1]], this.met))
    e <- this.met$met - yhat
    
    ## now, predict @ mean date (but with real age & GMrat) and add back in residual
    this.met$met.adj <- unname(predict(gam.model[[1]], this.met %>% mutate(dateNumeric = mean(metdata[[1]]$dateNumeric, na.rm=T)))) + e

    adj.df <- rbind(adj.df, this.met %>% mutate(metname = this_mod_met))
  }
}

adj.df.wide <- merge(
  adj.df %>% select(-met.adj) %>% pivot_wider(names_from = metname, values_from = met),
  adj.df %>% select(-met) %>% mutate(metname = paste0(metname, '.adj')) %>% pivot_wider(names_from = metname, values_from = met.adj),
  by = c('RECID','timepoint','age','region','hemi','roi','dateNumeric','GMrat','source'),
  all=T) %>% 
  select(RECID, timepoint, age, hemi, roi, dateNumeric, source, GMrat, GABA, GABA.adj, Glu, Glu.adj, Gln, Gln.adj, Glc, Glc.adj, GSH, GSH.adj, mI, mI.adj, NAA, NAA.adj, Tau, Tau.adj)

dim(adj.df.wide)
View(adj.df.wide)

write.csv(x = adj.df.wide, file = '~/scratch/deepak/sarpal_mrsi_adj_all_20250212.csv')


ggplot(data = adj.df.wide %>% filter(roi == 'Thalamus'), 
       aes(x = Glu, y = Glu.adj, color=as.factor(hemi), shape=as.factor(timepoint))) +
  geom_point() +
  geom_path(aes(group = interaction(RECID, hemi))) +
  stat_smooth(method='lm', alpha=0.2) +
  geom_abline(slope=1, intercept = 0, linetype=2) +
  theme(legend.position=c(.1, .8))


#gam.model <- gam_models %>% filter(met == 'Glu', biregion == 'DLPFC') %>% pull(model)
#print(plot(getViz(gam.model[[1]]), allTerms = T), pages = 1)

head(szmet)
szmet %>% select(region) %>% distinct()

unique(gam_models$biregion)
unique(gam_models$met)

print(plot(getViz(gam.model[[1]]), allTerms = T), pages = 1)


szdata <- adj.df.wide %>% filter(roi == 'Caudate')
gam.model <- gam_models %>% filter(met == 'GABA', biregion == 'Caudate') %>% pull(model)
metdata <- gam_models %>% filter(met == 'GABA', biregion == 'Caudate') %>% pull(metdata)

## adjust original data
luna_yhat <- unname(predict(gam.model[[1]], metdata[[1]]))
luna_e <- metdata[[1]]$Cr - luna_yhat
metdata[[1]]$met.adj <- unname(predict(gam.model[[1]], metdata[[1]] %>% mutate(dateNumeric = mean(metdata[[1]]$dateNumeric, na.rm=T)))) + luna_e


## combine
all.data <- rbind(szdata %>% filter(timepoint == 1) %>% select(age, GABA.adj) %>% mutate(group = 'SZ1'),
      szdata %>% filter(timepoint == 2) %>% select(age, GABA.adj) %>% mutate(group = 'SZ2'),
      metdata[[1]] %>% filter(age >= 18) %>% select(age, GABA.adj = met.adj) %>% mutate(group = 'Cn'))

## plot
ggplot(data = all.data, aes(x = age, y = GABA.adj, color = group)) +
  geom_point() +
  stat_smooth(method = 'gam', alpha = 0.1)

ggplot(data = all.data, aes(x = group, y = GABA.adj, fill=group)) + 
#  geom_violin(draw_quantiles = c(.25, .5, .75))
  geom_boxplot()

# Compute the mean and SEM for each group
data_summary <- all.data %>%
  group_by(group) %>%
  summarise(
    mean_GABA = mean(GABA.adj, na.rm = TRUE),
    sem_GABA = sd(GABA.adj, na.rm = TRUE) / sqrt(n())
  )

# Create the bar plot with error bars
ggplot(data_summary, aes(x = group, y = mean_GABA, fill = group)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = mean_GABA - sem_GABA, ymax = mean_GABA + sem_GABA), 
                width = 0.2) +
  scale_fill_manual(values = c("SZ1" = "blue", "SZ2" = "red", "Cn" = "gray")) +
  labs(x = "Group", y = "Glu (Mean ± SEM)") +
  theme_minimal() +
  theme(legend.position = "none")

