# demographics by region
# would be better as org-babel or RMarkdown document
# initial 20220725. revisited 20220815
source("MRS_GabaGlu.R")

# keep_rois <- c(1,2,7,8,9,10) # from MRS_GabaGlu.R
MRS_GG_keep <- MRS_GabaGlu %>% filter(roi %in% keep_rois)
length(unique(MRS_GG_keep$ld8)) # == 137 

group_demo <- MRS_GG_keep %>%
   group_by(agegrp,sex,region) %>%
   summarise(n=length(unique(ld8)))

ggplot(group_demo) +
   aes(x=agegrp,y=n,fill=sex) +
   geom_col(position='dodge') +
   facet_wrap(~region)

d <- MRS_GabaGlu %>%  filter(roi %in% keep_rois)
plot_grid(ggplot(d) + aes(x=region,fill=agegrp)+ geom_bar(position='dodge'),
          ggplot(d) + aes(x=region,fill=sex)+ geom_bar(position='dodge'),nrow=2)


# 20220815
MRS_GG_keep %>%
    group_by(region,sex) %>%
    summarise(n=length(unique(ld8)), y=min(age), o=max(age)) %>%
    pivot_wider(c("region"), names_from="sex", values_from = c("n","y","o")) %>%
    mutate(n=n_M+n_F, y=min(c_across(c(y_F,y_M))), o=max(c_across(c(o_F,o_M))))

#region       n_F   n_M   y_F   y_M   o_F   o_M     n     y     o
#ACC           57    65  10.6  10.2  29.6  30.4   122  10.2  30.4
#LAntInsula    55    60  10.6  10.2  29.6  30.4   115  10.2  30.4
#LDLPFC        47    54  10.6  10.2  29.6  30.4   101  10.2  30.4
#MPFC          56    57  10.6  11.5  29.6  29.5   113  10.6  29.6
#RAntInsula    62    61  10.6  10.2  29.6  30.4   123  10.2  30.4
#RDLPFC        46    53  10.9  10.2  29.6  29.5    99  10.2  29.6
