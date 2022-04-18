#### Code for Paper ####
# 5/17/21 

library(tidyverse)
library(ggplot2)
library(readxl)
library(Hmisc)
library(lmerTest)
library(corrplot)
library(RColorBrewer)
library(data.table)
library(mediation)

#### Get data and remove bad quality data ####
source("readMRS.R") # makes MRS
source("MRS_glu.R")
z_thres = 2

#### Participant histogram ####
# doing this for the full sample prior to data quality exclusions
# pick one ROI as representative 
ROI1 <- MRS %>% filter(roi==1) # n = 144
agehist <- ggplot(ROI1, aes(x=age,fill=sex)) + geom_histogram(color="black", binwidth = 1) + theme_classic(base_size = 15) + xlab("Age (years)") + ylab("Count") + labs(fill = "Sex")

ggsave(agehist, filename = "0519_agehistogram.pdf", width = 7, height = 6, units = "in", dpi = 300)

#### Glutamate and Age ####

# Create dataframe with good quality Glutamate data 
MRS_glu <- MRS %>% filter(Glu.SD <=20)

# ROI 1 (R Anterior Insula) and 2 (L Anterior Insula)

ROI12_Glu <- MRS_glu %>% filter(roi == 1 | roi == 2)
ROI1_Glu <- MRS_glu %>% filter(roi == 1)
ROI2_Glu <- MRS_glu %>% filter(roi == 2)


ROI12_Glu <- ROI12_Glu%>%
  mutate(zscore_glu=scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
ROI12_Glu <- ROI12_Glu %>% 
  filter(abs(zscore_glu) <= z_thres)

ROI1_Glu <- ROI1_Glu%>%
  mutate(zscore_glu=scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
ROI1_Glu <- ROI1_Glu %>% 
  filter(abs(zscore_glu) <= z_thres)

ROI2_Glu <- ROI2_Glu%>%
  mutate(zscore_glu=scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
ROI2_Glu <- ROI2_Glu %>% 
  filter(abs(zscore_glu) <= z_thres)


#ROI12_Glu %>% select(ld8, label) %>% spread(label, -ld8)

# pick best model
ROI12_Glu_age <- lmer(data=ROI12_Glu, Glu.Cr ~ age + label + sex + GMrat + (1|ld8))
summary(ROI12_Glu_age)
ROI12_Glu_invage <- lmer(data=ROI12_Glu, Glu.Cr ~ invage + label + sex + GMrat + (1|ld8))
summary(ROI12_Glu_invage) 
ROI12_Glu_quadage <- lmer(data=ROI12_Glu, Glu.Cr ~ age + age2 + label + sex + GMrat + (1|ld8))
summary(ROI12_Glu_quadage) 

AIC(ROI12_Glu_age)
AIC(ROI12_Glu_invage) # best fit by a lot
AIC(ROI12_Glu_quadage)

# effect sizes
ROI12_zGlu <- lmer(data=ROI12_Glu, zscore_glu ~ zscore_invage + label + (1|ld8))
summary(ROI12_zGlu) # just glu ~ age effect

ROI12_zGlu_invage <- lmer(data=ROI12_Glu, zscore_glu ~ zscore_invage + label + sex + zscore_gm + (1|ld8))
summary(ROI12_zGlu_invage) # all covariates for table
summary(ROI12_zGlu_invage)$coefficients

# test for interactions w/ age
ROI12_Glu_hemi_int <- lmer(data=ROI12_Glu, Glu.Cr ~ invage * label + sex + GMrat + (1|ld8))
summary(ROI12_Glu_hemi_int) 
ROI12_Glu_sex_int <- lmer(data=ROI12_Glu, Glu.Cr ~ invage * sex + label + GMrat + (1|ld8))
summary(ROI12_Glu_sex_int) 
ROI12_Glu_gmrat_int <- lmer(data=ROI12_Glu, Glu.Cr ~ invage * GMrat + sex + label + (1|ld8))
summary(ROI12_Glu_gmrat_int) 


ggplot(ROI12_Glu, aes(x=age, y=Glu.SD)) + geom_point() + geom_smooth(method="lm") + theme_classic()
aiglusd_age <- lmer(data=ROI12_Glu, scale(Glu.SD) ~ scale(age) +(1|ld8))
summary(aiglusd_age)


# variability by age 
var_Glu1 <- lm(Glu.Cr ~ invage, data=ROI1_Glu)
ROI1_Glu$Glu_ageResids <- abs(var_Glu1$residuals)
gluvar_1 <- lm(Glu_ageResids ~ invage, data=ROI1_Glu)
summary(gluvar_1)
gluvar_1 <- lm(Glu_ageResids ~ invage + GMrat + sex, data=ROI1_Glu)
summary(gluvar_1) # stays sig 

ggplot(ROI1_Glu, aes(x=age, y=Glu_ageResids)) + geom_point() + geom_smooth(method="lm") + theme_classic()
ROI1_Glu <- ROI1_Glu%>%
  mutate(zscore_Glu_ageResids=scale(Glu_ageResids), zscore_age=scale(age), zscore_GM = scale(GMrat))
zscore_gluvar_1 <- lm(zscore_Glu_ageResids ~ zscore_invage, data=ROI1_Glu)
summary(zscore_gluvar_1)
zscore_gluvar_1 <- lm(zscore_Glu_ageResids ~ zscore_invage + zscore_GM + sex, data=ROI1_Glu)
summary(zscore_gluvar_1)

roi1_var <- ggplot(ROI1_Glu) + aes(y=Glu_ageResids, x=age) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Variability") + 
  ggtitle("Right AIns") +
  theme(legend.position = "none")
ggsave(roi1_var, filename = "1203_RAI_var.pdf", width = 5, height = 5, units = "in", dpi = 300)

var_Glu2 <- lm(Glu.Cr ~ invage, data=ROI2_Glu)
ROI2_Glu$Glu_ageResids <- abs(var_Glu2$residuals)
gluvar_2 <- lm(Glu_ageResids ~ invage, data=ROI2_Glu)
summary(gluvar_2)
ggplot(ROI2_Glu, aes(x=age, y=Glu_ageResids)) + geom_point() + geom_smooth(method="lm") + theme_classic()
ROI2_Glu <- ROI2_Glu%>%
  mutate(zscore_Glu_ageResids=scale(Glu_ageResids), zscore_age=scale(age), zscore_GM = scale(GMrat))
zscore_gluvar_2 <- lm(zscore_Glu_ageResids ~ zscore_invage + zscore_GM + sex, data=ROI2_Glu)
summary(zscore_gluvar_2)
roi2_var <- ggplot(ROI2_Glu) + aes(y=Glu_ageResids, x=age) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Variability") + 
  ggtitle("Left AIns") +
  theme(legend.position = "none")
ggsave(roi2_var, filename = "1203_LAI_var.pdf", width = 5, height = 5, units = "in", dpi = 300)


ggplot(ROI1_Glu, aes(x=age, y=Glu.SD)) + geom_point() + geom_smooth(method="lm") + theme_classic()
aiglusd_age <- lm(data=ROI1_Glu, scale(Glu.SD) ~ scale(age))
summary(aiglusd_age)


ggplot(ROI2_Glu, aes(x=age, y=Glu.SD)) + geom_point() + geom_smooth(method="lm") + theme_classic()
aiglusd_age <- lm(data=ROI2_Glu, scale(Glu.SD) ~ scale(age))
summary(aiglusd_age)

# ROI 7 (ACC)
ROI7_Glu <- MRS_glu %>% filter(roi == 7)

ROI7_Glu <- ROI7_Glu%>%
  mutate(zscore_glu=scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))

ROI7_Glu <- ROI7_Glu %>% 
  filter(abs(zscore_glu) <= z_thres)

# ERROR HERE: no label?
#ROI7_Glu %>% select(ld8, label) %>% spread(label, -ld8)

# pick best model
ROI7_Glu_age <- lm(data=ROI7_Glu, Glu.Cr ~ age + sex + GMrat)
summary(ROI7_Glu_age)
ROI7_Glu_invage <- lm(data=ROI7_Glu, Glu.Cr ~ invage + sex + GMrat)
summary(ROI7_Glu_invage)
ROI7_Glu_quadage <- lm(data=ROI7_Glu, Glu.Cr ~ age + age2 + sex + GMrat)
summary(ROI7_Glu_quadage)

AIC(ROI7_Glu_age) 
AIC(ROI7_Glu_invage)
AIC(ROI7_Glu_quadage)

#effect sizes
ROI7_zGlu <- lm(data=ROI7_Glu, zscore_glu ~ zscore_invage)
summary(ROI7_zGlu)
ROI7_zGlu_invage <- lm(data=ROI7_Glu, zscore_glu ~ zscore_invage + sex + zscore_gm)
summary(ROI7_zGlu_invage)

# test for interactions
ROI7_age_sex_int <- lm(data=ROI7_Glu, Glu.Cr ~ invage * sex + GMrat)
summary(ROI7_age_sex_int)
ROI7_gm_sex_int <- lm(data=ROI7_Glu, Glu.Cr ~ invage * GMrat + sex)
summary(ROI7_gm_sex_int)

ggplot(ROI7_Glu, aes(x=age, y=Glu.Cr)) + geom_point() + geom_smooth(method="lm") + theme_classic()

ggplot(ROI7_Glu, aes(x=age, y=Glu.SD)) + geom_point() + geom_smooth(method="lm") + theme_classic()
accglusd_age <- lm(data=ROI7_Glu, scale(Glu.SD) ~ scale(age))
summary(accglusd_age)

#variability and age 
var_Glu7 <- lm(Glu.Cr ~ invage, data=ROI7_Glu)
ROI7_Glu$Glu_ageResids <- abs(var_Glu7$residuals)
gluvar_7 <- lm(Glu_ageResids ~ invage, data=ROI7_Glu)
summary(gluvar_7)
gluvar_7 <- lm(Glu_ageResids ~ invage + GMrat + sex, data=ROI7_Glu)
summary(gluvar_7)

ggplot(ROI7_Glu, aes(x=age, y=Glu_ageResids)) + geom_point() + geom_smooth(method="lm") + theme_classic()
ROI7_Glu <- ROI7_Glu%>%
  mutate(zscore_Glu_ageResids=scale(Glu_ageResids), zscore_age=scale(age), zscore_GM = scale(GMrat))
zscore_gluvar_7 <- lm(zscore_Glu_ageResids ~ zscore_invage + sex + zscore_GM, data=ROI7_Glu)
summary(zscore_gluvar_7)

roi7_var <- ggplot(ROI7_Glu) + aes(y=Glu_ageResids, x=age) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Variability") + 
  ggtitle("ACC") +
  theme(legend.position = "none")
ggsave(roi7_var, filename = "1203_ACC_var.pdf", width = 5, height = 5, units = "in", dpi = 300)




# ROI 8 (MPFC)
ROI8_Glu <- MRS_glu %>% filter( roi == 8)

ROI8_Glu <- ROI8_Glu%>%
  mutate(zscore_glu=scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))

ROI8_Glu <- ROI8_Glu %>% 
  filter(abs(zscore_glu) <= z_thres)

# ERROR HERE
#ROI8_Glu %>% select(ld8, label) %>% spread(label, -ld8)

ROI8_Glu_age <- lm(data=ROI8_Glu, Glu.Cr ~ age + sex + GMrat)
summary(ROI8_Glu_age)
ROI8_Glu_invage <- lm(data=ROI8_Glu, Glu.Cr ~ invage + sex + GMrat)
summary(ROI8_Glu_invage)
ROI8_Glu_quadage <- lm(data=ROI8_Glu, Glu.Cr ~ age + age2 + sex + GMrat)
summary(ROI8_Glu_quadage)

AIC(ROI8_Glu_age)
AIC(ROI8_Glu_invage)
AIC(ROI8_Glu_quadage)

ROI8_zGlu <- lm(data=ROI8_Glu, zscore_glu ~ zscore_invage)
summary(ROI8_zGlu)
ROI8_zGlu_invage <- lm(data=ROI8_Glu, zscore_glu ~ zscore_invage + sex + zscore_gm)
summary(ROI8_zGlu_invage)

# test for interactions
ROI8_age_sex_int <- lm(data=ROI8_Glu, Glu.Cr ~ invage * sex + GMrat)
summary(ROI8_age_sex_int)
ROI8_gm_sex_int <- lm(data=ROI8_Glu, Glu.Cr ~ invage * GMrat + sex)
summary(ROI8_gm_sex_int)

ggplot(ROI8_Glu, aes(x=age, y=Glu.Cr)) + geom_point() + geom_smooth(method="lm") + theme_classic()

ggplot(ROI8_Glu, aes(x=age, y=Glu.SD)) + geom_point() + geom_smooth(method="lm") + theme_classic()
mpfcglusd_age <- lm(data=ROI8_Glu, scale(Glu.SD) ~ scale(age))
summary(mpfcglusd_age)

var_Glu8 <- lm(Glu.Cr ~ invage, data=ROI8_Glu)
ROI8_Glu$Glu_ageResids <- abs(var_Glu8$residuals)
gluvar_8 <- lm(Glu_ageResids ~ invage, data=ROI8_Glu)
summary(gluvar_8)
gluvar_8 <- lm(Glu_ageResids ~ invage + GMrat + sex, data=ROI8_Glu)
summary(gluvar_8)

ggplot(ROI8_Glu, aes(x=age, y=Glu_ageResids)) + geom_point() + geom_smooth(method="lm") + theme_classic()


ROI8_Glu <- ROI8_Glu%>%
  mutate(zscore_Glu_ageResids=scale(Glu_ageResids))
zscore_gluvar_8 <- lm(zscore_Glu_ageResids ~ zscore_invage + zscore_gm + sex, data=ROI8_Glu)
summary(zscore_gluvar_8)

roi8_var <- ggplot(ROI8_Glu) + aes(y=Glu_ageResids, x=age) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Variability") + 
  ggtitle("MPFC") +
  theme(legend.position = "none")
#ggsave(roi8_var, filename = "1203_MPFC_var.pdf", width = 5, height = 5, units = "in", dpi = 300)


# ROI 9 ( R DLPFC) and 10 (L DLPFC)

ROI910_Glu <- MRS_glu %>% filter(roi == 9 | roi == 10)
ROI910_Glu <- ROI910_Glu%>%
  mutate(zscore_glu=scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))


ROI910_Glu <- ROI910_Glu %>% 
  filter(abs(zscore_glu) <= z_thres)

ROI9_Glu <- MRS_glu %>% filter(roi == 9)
ROI9_Glu <- ROI9_Glu%>%
  mutate(zscore_glu=scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))


ROI9_Glu <- ROI9_Glu %>% 
  filter(abs(zscore_glu) <= z_thres)

ROI10_Glu <- MRS_glu %>% filter(roi == 10)
ROI10_Glu <- ROI10_Glu%>%
  mutate(zscore_glu=scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))


ROI10_Glu <- ROI10_Glu %>% 
  filter(abs(zscore_glu) <= z_thres)



ROI910_Glu_age <- lmer(data=ROI910_Glu, Glu.Cr ~ age + label + sex + GMrat + (1|ld8))
summary(ROI910_Glu_age) 

ROI910_Glu_invage <- lmer(data=ROI910_Glu, Glu.Cr ~ invage + label + sex + GMrat + (1|ld8))
summary(ROI910_Glu_invage) 

ROI910_Glu_quadage <- lmer(data=ROI910_Glu, Glu.Cr ~ age + age2+ label + sex + GMrat + (1|ld8))
summary(ROI910_Glu_quadage) 

AIC(ROI910_Glu_age)
AIC(ROI910_Glu_invage)
AIC(ROI910_Glu_quadage)

ROI910_Glu_invagez <- lmer(data=ROI910_Glu, zscore_glu ~ zscore_invage + label + sex + zscore_gm + (1|ld8))
summary(ROI910_Glu_invagez) 
summary(ROI910_Glu_invagez)$coefficients

ROI910_zGlu <- lmer(data=ROI910_Glu, zscore_glu ~ zscore_invage + label + (1|ld8))
summary(ROI910_zGlu)

# test for interactions w/ age
ROI910_age_hemi_int <- lmer(data=ROI910_Glu, Glu.Cr ~ invage * label + sex + GMrat + (1|ld8))
summary(ROI910_age_hemi_int) 
ROI910_age_sex_int <- lmer(data=ROI910_Glu, Glu.Cr ~ invage * sex + label + GMrat + (1|ld8))
summary(ROI910_age_sex_int) 
ROI910_age_gmrat_int <- lmer(data=ROI910_Glu, Glu.Cr ~ invage * GMrat + sex + label + (1|ld8))
summary(ROI910_age_gmrat_int) 


ggplot(ROI910_Glu, aes(x=age, y=Glu.SD, color=label)) + geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + theme_classic()
glusd_age <- lm(data=ROI910_Glu, scale(Glu.SD) ~ scale(age))
summary(glusd_age)

var_Glu9 <- lm(Glu.Cr ~ invage, data=ROI9_Glu)
ROI9_Glu$Glu_ageResids <- abs(var_Glu9$residuals)
gluvar_9 <- lm(Glu_ageResids ~ invage, data=ROI9_Glu)
summary(gluvar_9)
ggplot(ROI9_Glu, aes(x=age, y=Glu_ageResids)) + geom_point() + geom_smooth(method="lm") + theme_classic()


roi9_var <- ggplot(ROI9_Glu) + aes(y=Glu_ageResids, x=age) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Variability") + 
  ggtitle("Right DLPFC") +
  theme(legend.position = "none")
#ggsave(roi9_var, filename = "1203_RDLPFC_var.pdf", width = 5, height = 5, units = "in", dpi = 300)

ROI9_Glu <- ROI9_Glu%>%
  mutate(zscore_Glu_ageResids=scale(Glu_ageResids), zscore_age=scale(age))
zscore_gluvar_9 <- lm(zscore_Glu_ageResids ~ zscore_invage + zscore_gm + sex, data=ROI9_Glu)
summary(zscore_gluvar_9)

var_Glu10 <- lm(Glu.Cr ~ invage, data=ROI10_Glu)
ROI10_Glu$Glu_ageResids <- abs(var_Glu10$residuals)
gluvar_10 <- lm(Glu_ageResids ~ invage, data=ROI10_Glu)
summary(gluvar_10)
ggplot(ROI10_Glu, aes(x=age, y=Glu_ageResids)) + geom_point() + geom_smooth(method="lm") + theme_classic()

ROI10_Glu <- ROI10_Glu%>%
  mutate(zscore_Glu_ageResids=scale(Glu_ageResids), zscore_age=scale(age))
zscore_gluvar_10 <- lm(zscore_Glu_ageResids ~ zscore_invage + zscore_gm + sex, data=ROI10_Glu)
summary(zscore_gluvar_10)

roi10_var <- ggplot(ROI10_Glu) + aes(y=Glu_ageResids, x=age) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Variability") + 
  ggtitle("Left DLPFC") +
  theme(legend.position = "none")
#ggsave(roi10_var, filename = "1203_LDLPFC_var.pdf", width = 5, height = 5, units = "in", dpi = 300)



#### GABA and Age ####

# Create dataframe with good quality GABA data 
MRS_GABA <- MRS %>% filter(GABA.SD <=20)
z_thres = 2

# ROI 1 (R Anterior Insula) and 2 (L Anterior Insula)

ROI12_GABA <- MRS_GABA %>% filter(roi == 1 | roi == 2)
ROI12_GABA <- ROI12_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))

ROI12_GABA <- ROI12_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)

ROI12_GABA %>% select(ld8, label) %>% spread(label, -ld8)

ROI1_GABA <- MRS_GABA %>% filter(roi == 1)
ROI1_GABA <- ROI1_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))

ROI1_GABA <- ROI1_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)


var_GABA1 <- lm(GABA.Cr ~ invage, data=ROI1_GABA)
ROI1_GABA$GABA_ageResids <- abs(var_GABA1$residuals)
gabavar_1 <- lm(GABA_ageResids ~ invage, data=ROI1_GABA)
summary(gabavar_1)
ggplot(ROI1_GABA, aes(x=age, y=GABA_ageResids)) + geom_point() + geom_smooth(method="lm") + theme_classic()

ROI1_GABA <- ROI1_GABA%>%
  mutate(zscore_GABA_ageResids=scale(GABA_ageResids))
zscore_gabavar_1 <- lm(zscore_GABA_ageResids ~ zscore_invage + zscore_gm + sex, data=ROI1_GABA)
summary(zscore_gabavar_1)

roi1_vargab <- ggplot(ROI1_GABA) + aes(y=GABA_ageResids, x=age) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Variability") + 
  ggtitle("Right AIns") +
  theme(legend.position = "none")
#ggsave(roi1_vargab, filename = "1203_RAI_vargab.pdf", width = 5, height = 5, units = "in", dpi = 300)

ROI2_GABA <- MRS_GABA %>% filter(roi == 2)
ROI2_GABA <- ROI2_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))

ROI2_GABA <- ROI2_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)


var_GABA2 <- lm(GABA.Cr ~ invage, data=ROI2_GABA)
ROI2_GABA$GABA_ageResids <- abs(var_GABA2$residuals)
gabavar_2 <- lm(GABA_ageResids ~ invage, data=ROI2_GABA)
summary(gabavar_2)
ggplot(ROI2_GABA, aes(x=age, y=GABA_ageResids)) + geom_point() + geom_smooth(method="lm") + theme_classic()

ROI2_GABA <- ROI2_GABA%>%
  mutate(zscore_GABA_ageResids=scale(GABA_ageResids))
zscore_gabavar_2 <- lm(zscore_GABA_ageResids ~ zscore_invage + zscore_gm + sex, data=ROI2_GABA)
summary(zscore_gabavar_2)

roi2_vargab <- ggplot(ROI2_GABA) + aes(y=GABA_ageResids, x=age) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Variability") + 
  ggtitle("Left AIns") +
  theme(legend.position = "none")
#ggsave(roi2_vargab, filename = "1203_LAI_vargab.pdf", width = 5, height = 5, units = "in", dpi = 300)


# pick best model
ROI12_GABA_age <- lmer(data=ROI12_GABA, GABA.Cr ~ age + label + sex + GMrat + (1|ld8))
summary(ROI12_GABA_age)
ROI12_GABA_invage <- lmer(data=ROI12_GABA, GABA.Cr ~ invage + label + sex + GMrat + (1|ld8))
summary(ROI12_GABA_invage) 
ROI12_GABA_quadage <- lmer(data=ROI12_GABA, GABA.Cr ~ age + age2 + label + sex + GMrat + (1|ld8))
summary(ROI12_GABA_quadage) 

AIC(ROI12_GABA_age)
AIC(ROI12_GABA_invage)
AIC(ROI12_GABA_quadage)

# effect sizes
ROI12_zGABA <- lmer(data=ROI12_GABA, zscore_GABA ~ zscore_invage + label + (1|ld8))
summary(ROI12_zGABA) # just GABA ~ age effect

ROI12_zGABA_invage <- lmer(data=ROI12_GABA, zscore_GABA ~ zscore_invage + label + sex + zscore_gm + (1|ld8))
summary(ROI12_zGABA_invage) # all covariates for table
summary(ROI12_zGABA_invage)$coefficients

# test for interactions w/ age
ROI12_GABA_hemi_int <- lmer(data=ROI12_GABA, GABA.Cr ~ invage * label + sex + GMrat + (1|ld8))
summary(ROI12_GABA_hemi_int) 
ROI12_GABA_hemi_int <- lmer(data=ROI12_GABA, zscore_GABA ~ zscore_invage * label + sex + zscore_gm + (1|ld8))
summary(ROI12_GABA_hemi_int) 
summary(ROI12_GABA_hemi_int)$coefficients
ROI12_GABA_sex_int <- lmer(data=ROI12_GABA, GABA.Cr ~ invage * sex + label + GMrat + (1|ld8))
summary(ROI12_GABA_sex_int) 
ROI12_GABA_gmrat_int <- lmer(data=ROI12_GABA, GABA.Cr ~ invage * GMrat + sex + label + (1|ld8))
summary(ROI12_GABA_gmrat_int) 


ggplot(ROI12_GABA, aes(x=age, y=GABA.SD)) + geom_point() + geom_smooth(method="lm") + theme_classic()
aigabasd_age <- lmer(data=ROI12_GABA, scale(GABA.SD) ~ scale(age) +(1|ld8))
summary(aigabasd_age)

# ROI 7 (ACC)
ROI7_GABA <- MRS_GABA %>% filter(roi == 7)
ROI7_GABA <- ROI7_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))

ROI7_GABA <- ROI7_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)

#variability and age 
var_GABA7 <- lm(GABA.Cr ~ invage, data=ROI7_GABA)
ROI7_GABA$GABA_ageResids <- abs(var_GABA7$residuals)
gabavar_7 <- lm(GABA_ageResids ~ invage, data=ROI7_GABA)
summary(gabavar_7)
ggplot(ROI7_GABA, aes(x=age, y=GABA_ageResids)) + geom_point() + geom_smooth(method="lm") + theme_classic()

ROI7_GABA <- ROI7_GABA%>%
  mutate(zscore_GABA_ageResids=scale(GABA_ageResids))
zscore_gabavar_7 <- lm(zscore_GABA_ageResids ~ zscore_invage + zscore_gm + sex, data=ROI7_GABA)
summary(zscore_gabavar_7)

roi7_vargab <- ggplot(ROI7_GABA) + aes(y=GABA_ageResids, x=age) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Variability") + 
  ggtitle("ACC") +
  theme(legend.position = "none")
#ggsave(roi7_vargab, filename = "1203_ACC_vargab.pdf", width = 5, height = 5, units = "in", dpi = 300)


ROI7_GABA %>% select(ld8, label) %>% spread(label, -ld8)

# pick best model
ROI7_GABA_age <- lm(data=ROI7_GABA, GABA.Cr ~ age + sex + GMrat)
summary(ROI7_GABA_age)
ROI7_GABA_invage <- lm(data=ROI7_GABA, GABA.Cr ~ invage + sex + GMrat)
summary(ROI7_GABA_invage)
ROI7_GABA_quadage <- lm(data=ROI7_GABA, GABA.Cr ~ age + age2 + sex + GMrat)
summary(ROI7_GABA_quadage)

AIC(ROI7_GABA_age) 
AIC(ROI7_GABA_invage) # best fit
AIC(ROI7_GABA_quadage)

#effect sizes
ROI7_zGABA <- lm(data=ROI7_GABA, zscore_GABA ~ zscore_invage)
summary(ROI7_zGABA)
ROI7_zGABA_invage <- lm(data=ROI7_GABA, zscore_GABA ~ zscore_invage + sex + zscore_gm)
summary(ROI7_zGABA_invage)

# test for interactions
ROI7_age_sex_int <- lm(data=ROI7_GABA, GABA.Cr ~ invage * sex + GMrat)
summary(ROI7_age_sex_int)
ROI7_gm_sex_int <- lm(data=ROI7_GABA, GABA.Cr ~ invage * GMrat + sex)
summary(ROI7_gm_sex_int)


ggplot(ROI7_GABA, aes(x=age, y=GABA.Cr)) + geom_point() + geom_smooth(method="lm") + theme_classic()

ggplot(ROI7_GABA, aes(x=age, y=GABA.SD)) + geom_point() + geom_smooth(method="lm") + theme_classic()
accgabasd_age <- lm(data=ROI7_GABA, scale(GABA.SD) ~ scale(age))
summary(accgabasd_age)


# ROI 8 (MPFC)
ROI8_GABA <- MRS_GABA %>% filter( roi == 8)
ROI8_GABA <- ROI8_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))

ROI8_GABA <- ROI8_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)

#variability and age 
var_GABA8 <- lm(GABA.Cr ~ invage, data=ROI8_GABA)
ROI8_GABA$GABA_ageResids <- abs(var_GABA8$residuals)
gabavar_8 <- lm(GABA_ageResids ~ invage, data=ROI8_GABA)
summary(gabavar_8)
ggplot(ROI8_GABA, aes(x=age, y=GABA_ageResids)) + geom_point() + geom_smooth(method="lm") + theme_classic()

ROI8_GABA <- ROI8_GABA%>%
  mutate(zscore_GABA_ageResids=scale(GABA_ageResids))
zscore_gabavar_8 <- lm(zscore_GABA_ageResids ~ zscore_invage + zscore_gm + sex, data=ROI8_GABA)
summary(zscore_gabavar_8)

roi8_vargab <- ggplot(ROI8_GABA) + aes(y=GABA_ageResids, x=age) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Variability") + 
  ggtitle("MPFC") +
  theme(legend.position = "none")
#ggsave(roi8_vargab, filename = "1203_MPFC_vargab.pdf", width = 5, height = 5, units = "in", dpi = 300)




ROI8_GABA %>% select(ld8, label) %>% spread(label, -ld8)

ROI8_GABA_age <- lm(data=ROI8_GABA, GABA.Cr ~ age + sex + GMrat)
summary(ROI8_GABA_age)
ROI8_GABA_invage <- lm(data=ROI8_GABA, GABA.Cr ~ invage + sex + GMrat)
summary(ROI8_GABA_invage)
ROI8_GABA_quadage <- lm(data=ROI8_GABA, GABA.Cr ~ age + age2 + sex + GMrat)
summary(ROI8_GABA_quadage)

AIC(ROI8_GABA_age)
AIC(ROI8_GABA_invage) # best fit
AIC(ROI8_GABA_quadage)

ROI8_zGABA <- lm(data=ROI8_GABA, zscore_GABA ~ zscore_invage)
summary(ROI8_zGABA)
ROI8_zGABA_invage <- lm(data=ROI8_GABA, zscore_GABA ~ zscore_invage + sex + zscore_gm)
summary(ROI8_zGABA_invage)

# test for interactions
ROI8_age_sex_int <- lm(data=ROI8_GABA, GABA.Cr ~ invage * sex + GMrat)
summary(ROI8_age_sex_int)
ROI8_gm_sex_int <- lm(data=ROI8_GABA, GABA.Cr ~ invage * GMrat + sex)
summary(ROI8_gm_sex_int)

ggplot(ROI8_GABA, aes(x=age, y=GABA.Cr)) + geom_point() + geom_smooth(method="lm") + theme_classic()
ggplot(ROI8_GABA, aes(x=age, y=GABA.SD)) + geom_point() + geom_smooth(method="lm") + theme_classic()
mpfcgabasd_age <- lm(data=ROI8_GABA, scale(GABA.SD) ~ scale(age))
summary(mpfcgabasd_age)
ggplot(ROI8_GABA, aes(x=GABA.Cr, y=GABA.SD)) + geom_point() + geom_smooth(method="lm") + theme_classic()
mpfcgabacr_age <- lm(data=ROI8_GABA, scale(GABA.SD) ~ scale(GABA.Cr))
summary(mpfcgabacr_age)




# ROI 9 ( R DLPFC) and 10 (L DLPFC)
ROI910_GABA <- MRS_GABA %>% filter(roi == 9 | roi == 10)
ROI910_GABA <- ROI910_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))


ROI910_GABA <- ROI910_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)

ROI9_GABA <- MRS_GABA %>% filter(roi == 9)
ROI9_GABA <- ROI9_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))


ROI9_GABA <- ROI9_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)

var_GABA9 <- lm(GABA.Cr ~ invage, data=ROI9_GABA)
ROI9_GABA$GABA_ageResids <- abs(var_GABA9$residuals)
gabavar_9 <- lm(GABA_ageResids ~ invage, data=ROI9_GABA)
summary(gabavar_9)
ggplot(ROI9_GABA, aes(x=age, y=GABA_ageResids)) + geom_point() + geom_smooth(method="lm") + theme_classic()

ROI9_GABA <- ROI9_GABA%>%
  mutate(zscore_GABA_ageResids=scale(GABA_ageResids))
zscore_gabavar_9 <- lm(zscore_GABA_ageResids ~ zscore_invage + zscore_gm + sex, data=ROI9_GABA)
summary(zscore_gabavar_9)

roi9_vargab <- ggplot(ROI9_GABA) + aes(y=GABA_ageResids, x=age) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Variability") + 
  ggtitle("Right DLPFC") +
  theme(legend.position = "none")
#ggsave(roi9_vargab, filename = "1203_RDLPFC_vargab.pdf", width = 5, height = 5, units = "in", dpi = 300)


ROI10_GABA <- MRS_GABA %>% filter(roi == 10)
ROI10_GABA <- ROI10_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))


ROI10_GABA <- ROI10_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)

var_GABA10 <- lm(GABA.Cr ~ invage, data=ROI10_GABA)
ROI10_GABA$GABA_ageResids <- abs(var_GABA10$residuals)
gabavar_10 <- lm(GABA_ageResids ~ invage, data=ROI10_GABA)
summary(gabavar_10)
ggplot(ROI10_GABA, aes(x=age, y=GABA_ageResids)) + geom_point() + geom_smooth(method="lm") + theme_classic()

ROI10_GABA <- ROI10_GABA%>%
  mutate(zscore_GABA_ageResids=scale(GABA_ageResids))
zscore_gabavar_10 <- lm(zscore_GABA_ageResids ~ zscore_invage + zscore_gm + sex, data=ROI10_GABA)
summary(zscore_gabavar_10)

roi10_vargab <- ggplot(ROI10_GABA) + aes(y=GABA_ageResids, x=age) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Variability") + 
  ggtitle("Left DLPFC") +
  theme(legend.position = "none")
#ggsave(roi10_vargab, filename = "1203_LDLPFC_vargab.pdf", width = 5, height = 5, units = "in", dpi = 300)



ROI910_GABA_invage <- lmer(data=ROI910_GABA, GABA.Cr ~ invage + label + sex + GMrat + (1|ld8))
summary(ROI910_GABA_invage) 
ROI910_GABA_age <- lmer(data=ROI910_GABA, GABA.Cr ~ age + label + sex + GMrat + (1|ld8))
summary(ROI910_GABA_age) 
ROI910_GABA_quadage <- lmer(data=ROI910_GABA, GABA.Cr ~ age + age2+ label + sex + GMrat + (1|ld8))
summary(ROI910_GABA_quadage) 


AIC(ROI910_GABA_age)
AIC(ROI910_GABA_invage)
AIC(ROI910_GABA_quadage)

ROI910_zGABA <- lmer(data=ROI910_GABA, zscore_GABA ~ zscore_invage + label + (1|ld8))
summary(ROI910_zGABA)

ROI910_GABA_invagez <- lmer(data=ROI910_GABA, zscore_GABA ~ zscore_invage + label + sex + zscore_gm + (1|ld8))
summary(ROI910_GABA_invagez) 
summary(ROI910_GABA_invagez)$coefficients



# test for interactions w/ age
ROI910_age_hemi_int <- lmer(data=ROI910_GABA, GABA.Cr ~ invage * label + sex + GMrat + (1|ld8))
summary(ROI910_age_hemi_int) 
ROI910_age_sex_int <- lmer(data=ROI910_GABA, GABA.Cr ~ invage * sex + label + GMrat + (1|ld8))
summary(ROI910_age_sex_int) 
ROI910_age_gmrat_int <- lmer(data=ROI910_GABA, GABA.Cr ~ invage * GMrat + sex + label + (1|ld8))
summary(ROI910_age_gmrat_int) 


ggplot(ROI910_GABA, aes(x=zscore_invage, y=GABA.Cr, color=label)) + geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + theme_classic()
ggplot(ROI910_GABA, aes(x=age, y=GABA.SD, color=label)) + geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + theme_classic()
gabasd_age <- lm(data=ROI910_GABA, scale(GABA.SD) ~ scale(age))
summary(gabasd_age)

#### Ratio and Age ####
# create dataframe keeping only people who have both good quality Glu and good GABA data to
# make a ratio out of 
#MRS_Ratio <- MRS_glu %>% filter(GABA.SD <=20)


MRS_glu <- MRS_glu %>%
  group_by(roi) %>%
  mutate(zscore=scale(Glu.Cr, center=T, scale=T)) %>%
  filter(abs(zscore) <= z_thres) %>% ungroup 

MRS_Ratio <- MRS_glu %>%
  filter(GABA.SD <= 20) %>%
  group_by(roi) %>%
  mutate(zscore=scale(GABA.Cr, center=T, scale=T)) %>%
  filter(abs(zscore) <= z_thres) %>% ungroup %>%
  mutate(agegrp = cut(age,
                      breaks=c(0,16,22,Inf),
                      labels=c("10-16","17-22", "23-30")))

MRS_Ratio$Ratio <- MRS_Ratio$Glu.Cr/MRS_Ratio$GABA.Cr

# ROI 1 (R Anterior Insula) and 2 (L Anterior Insula)

ROI12_Ratio <- MRS_Ratio %>% filter(roi == 1 | roi == 2)
ROI12_Ratio <- ROI12_Ratio%>%
  mutate(zscore_Ratio=scale(Ratio), zscore_invage=scale(invage), zscore_gm=scale(GMrat))


ROI12_Ratio <- ROI12_Ratio %>% 
 filter(abs(zscore_Ratio) <= z_thres)

ROI12_Ratio %>% select(ld8, label) %>% spread(label, -ld8)

# pick best model
ROI12_Ratio_age <- lmer(data=ROI12_Ratio, Ratio ~ age + label + sex + GMrat + (1|ld8))
summary(ROI12_Ratio_age)
ROI12_Ratio_invage <- lmer(data=ROI12_Ratio, Ratio ~ invage + label + sex + GMrat + (1|ld8))
summary(ROI12_Ratio_invage) 
ROI12_Ratio_quadage <- lmer(data=ROI12_Ratio, Ratio ~ age + age2 + label + sex + GMrat + (1|ld8))
summary(ROI12_Ratio_quadage) 

AIC(ROI12_Ratio_age)
AIC(ROI12_Ratio_invage)
AIC(ROI12_Ratio_quadage)

# effect sizes
ROI12_zRatio <- lmer(data=ROI12_Ratio, zscore_Ratio ~ zscore_invage + label + (1|ld8))
summary(ROI12_zRatio) # just Ratio ~ age effect

ROI12_zRatio_invage <- lmer(data=ROI12_Ratio, zscore_Ratio ~ zscore_invage + label + sex + zscore_gm + (1|ld8))
summary(ROI12_zRatio_invage) # all covariates for table

# test for interactions w/ age
ROI12_Ratio_hemi_int <- lmer(data=ROI12_Ratio, Ratio ~ invage * label + sex + GMrat + (1|ld8))
summary(ROI12_Ratio_hemi_int) 
ROI12_zRatio_invage <- lmer(data=ROI12_Ratio, zscore_Ratio ~ zscore_invage * label + sex + zscore_gm + (1|ld8))
summary(ROI12_zRatio_invage) # all covariates for table

ROI12_Ratio_sex_int <- lmer(data=ROI12_Ratio, Ratio ~ invage * sex + label + GMrat + (1|ld8))
summary(ROI12_Ratio_sex_int) 
ROI12_Ratio_gmrat_int <- lmer(data=ROI12_Ratio, Ratio ~ invage * GMrat + sex + label + (1|ld8))
summary(ROI12_Ratio_gmrat_int) 

ROI1 <- ROI12_Ratio %>% filter(roi==1)
ROI1_zRatio <- lm(data=ROI1, zscore_Ratio ~ zscore_invage + sex + zscore_gm)
summary(ROI1_zRatio) # just Ratio ~ age effect

ROI2 <- ROI12_Ratio %>% filter(roi==2)
ROI2_zRatio <- lm(data=ROI2, zscore_Ratio ~ zscore_invage + sex + zscore_gm)
summary(ROI2_zRatio) # just Ratio ~ age effect

# ROI 7 (ACC)
ROI7_Ratio <- MRS_Ratio %>% filter(roi == 7)
ROI7_Ratio <- ROI7_Ratio%>%
  mutate(zscore_Ratio=scale(Ratio), zscore_invage=scale(invage), zscore_gm=scale(GMrat))


ROI7_Ratio <- ROI7_Ratio %>% 
  filter(abs(zscore_Ratio) <= z_thres)

#ROI7_Ratio %>% select(ld8, label) %>% spread(label, -ld8)

# pick best model
ROI7_Ratio_age <- lm(data=ROI7_Ratio, Ratio ~ age + sex + GMrat)
summary(ROI7_Ratio_age)
ROI7_Ratio_invage <- lm(data=ROI7_Ratio, Ratio ~ invage + sex + GMrat)
summary(ROI7_Ratio_invage)
ROI7_Ratio_quadage <- lm(data=ROI7_Ratio, Ratio ~ age + age2 + sex + GMrat)
summary(ROI7_Ratio_quadage)

AIC(ROI7_Ratio_age) 
AIC(ROI7_Ratio_invage) # best fit
AIC(ROI7_Ratio_quadage)

#effect sizes
ROI7_zRatio <- lm(data=ROI7_Ratio, zscore_Ratio ~ zscore_invage)
summary(ROI7_zRatio)
ROI7_zRatio_invage <- lm(data=ROI7_Ratio, zscore_Ratio ~ zscore_invage + sex + zscore_gm)
summary(ROI7_zRatio_invage)

# test for interactions
ROI7_age_sex_int <- lm(data=ROI7_Ratio, Ratio ~ invage * sex + GMrat)
summary(ROI7_age_sex_int) # sig but wont survive
ROI7_gm_sex_int <- lm(data=ROI7_Ratio, Ratio ~ invage * GMrat + sex)
summary(ROI7_gm_sex_int)

ROI7_zRatio_invage <- lm(data=ROI7_Ratio, zscore_Ratio ~ zscore_invage * sex + zscore_gm)
summary(ROI7_zRatio_invage)

ggplot(ROI7_Ratio, aes(x=age, y=Ratio)) + geom_point() + geom_smooth(method="lm") + theme_classic()


# ROI 8 (MPFC)
ROI8_Ratio <- MRS_Ratio %>% filter( roi == 8)
ROI8_Ratio <- ROI8_Ratio%>%
  mutate(zscore_Ratio=scale(Ratio), zscore_invage=scale(invage), zscore_gm=scale(GMrat))


ROI8_Ratio <- ROI8_Ratio %>% 
  filter(abs(zscore_Ratio) <= z_thres)

#ROI8_Ratio %>% select(ld8, label) %>% spread(label, -ld8)

ROI8_Ratio_age <- lm(data=ROI8_Ratio, Ratio ~ age + sex + GMrat)
summary(ROI8_Ratio_age)
ROI8_Ratio_invage <- lm(data=ROI8_Ratio, Ratio ~ invage + sex + GMrat)
summary(ROI8_Ratio_invage)
ROI8_Ratio_quadage <- lm(data=ROI8_Ratio, Ratio ~ age + age2 + sex + GMrat)
summary(ROI8_Ratio_quadage)

AIC(ROI8_Ratio_age)
AIC(ROI8_Ratio_invage) # best fit
AIC(ROI8_Ratio_quadage)

ROI8_zRatio <- lm(data=ROI8_Ratio, zscore_Ratio ~ zscore_invage)
summary(ROI8_zRatio)
ROI8_zRatio_invage <- lm(data=ROI8_Ratio, zscore_Ratio ~ zscore_invage + sex + zscore_gm)
summary(ROI8_zRatio_invage)

# test for interactions
ROI8_age_sex_int <- lm(data=ROI8_Ratio, Ratio ~ invage * sex + GMrat)
summary(ROI8_age_sex_int)
ROI8_gm_sex_int <- lm(data=ROI8_Ratio, Ratio ~ invage * GMrat + sex)
summary(ROI8_gm_sex_int)

ggplot(ROI8_Ratio, aes(x=age, y=Ratio)) + geom_point() + geom_smooth(method="lm") + theme_classic()


# ROI 9 ( R DLPFC) and 10 (L DLPFC)
ROI910_Ratio <- MRS_Ratio %>% filter(roi == 9 | roi == 10)

#ROI910_Ratio %>% select(ld8,label) %>% gather(label)
ROI910_Ratio <- ROI910_Ratio%>%
  mutate(zscore_Ratio=scale(Ratio), zscore_invage=scale(invage), zscore_gm=scale(GMrat))


ROI910_Ratio <- ROI910_Ratio %>% 
  filter(abs(zscore_Ratio) <= z_thres)

#ROI910_Ratio %>% select(ld8, label) %>% spread(label, -ld8)

ROI910_Ratio_invage <- lmer(data=ROI910_Ratio, Ratio ~ invage + label + sex + GMrat + (1|ld8))
summary(ROI910_Ratio_invage) 
ROI910_Ratio_age <- lmer(data=ROI910_Ratio, Ratio ~ age + label + sex + GMrat + (1|ld8))
summary(ROI910_Ratio_age) 
ROI910_Ratio_quadage <- lmer(data=ROI910_Ratio, Ratio ~ age + age2+ label + sex + GMrat + (1|ld8))
summary(ROI910_Ratio_quadage) 
AIC(ROI910_Ratio_age)
AIC(ROI910_Ratio_invage)
AIC(ROI910_Ratio_quadage)


ROI910_zRatio <- lmer(data=ROI910_Ratio, zscore_Ratio ~ zscore_invage + label + (1|ld8))
summary(ROI910_zRatio)
ROI910_zRatio_invage <- lmer(data=ROI910_Ratio, zscore_Ratio ~ zscore_invage + label + sex + zscore_gm + (1|ld8))
summary(ROI910_zRatio_invage) 

# test for interactions w/ age
ROI910_age_hemi_int <- lmer(data=ROI910_Ratio, Ratio ~ invage * label + sex + GMrat + (1|ld8))
summary(ROI910_age_hemi_int) 
ROI910_age_sex_int <- lmer(data=ROI910_Ratio, Ratio ~ invage * sex + label + GMrat + (1|ld8))
summary(ROI910_age_sex_int) 
ROI910_age_gmrat_int <- lmer(data=ROI910_Ratio, Ratio ~ invage * GMrat + sex + label + (1|ld8))
summary(ROI910_age_gmrat_int) 

ggplot(ROI910_Ratio, aes(x=age, y=Ratio, color=label)) + geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + theme_classic()

#### Plot GABA and Glu on same graph ####
gaba12 <- ROI12_GABA %>% select(concentration=GABA.Cr, age, sex, GMrat, label) %>% mutate(metabolite="GABA")
glu12 <- ROI12_Glu %>% select(concentration=Glu.Cr, age, sex, GMrat, label) %>% mutate(metabolite="Glu")
gabaglu12 <- rbind(gaba12,glu12)
AI_gabaglu <- ggplot(gabaglu12) + aes(y=concentration, x=age, color=metabolite, linetype = label, shape=label) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Concentration (Metabolite/Cr)") + 
  ggtitle("AIns") +
  scale_shape_manual(values=c(16, 1)) + 
  scale_color_manual(values=c("#FF6666", "#66B2FF"), guide = guide_legend(reverse = TRUE)) +
  guides(linetype=FALSE) + 
  theme(legend.position = "none") + 
  coord_cartesian(ylim=c(0, 2.25))
#ggsave(AI_gabaglu, filename = "0715_AI_Glu_GABA.pdf", width = 6, height = 6, units = "in", dpi = 300)

#without legend
AI_gabaglu <- ggplot(gabaglu12) + aes(y=concentration, x=age, color=metabolite, linetype = label, shape=label) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Metabolite Level (/Cr)") + 
  ggtitle("AIns") +
  scale_shape_manual(values=c(16, 1)) + 
  scale_color_manual(values=c("#FF6666", "#66B2FF"), guide = guide_legend(reverse = TRUE)) +
  guides(linetype=FALSE) + 
  theme(legend.position = "none") + coord_cartesian(ylim=c(0, 2.10))
#ggsave(AI_gabaglu, filename = "0908_AI_Glu_GABA.pdf", width = 6, height = 6, units = "in", dpi = 300)


gaba910 <- ROI910_GABA %>% select(concentration=GABA.Cr, age, sex, GMrat, label) %>% mutate(metabolite="GABA")
glu910 <- ROI910_Glu %>% select(concentration=Glu.Cr, age, sex, GMrat, label) %>% mutate(metabolite="Glu")
gabaglu910 <- rbind(gaba910,glu910)
# with legend
DLPFC_gabaglu <- ggplot(gabaglu910) + aes(y=concentration, x=age, color=metabolite, linetype = label, shape=label) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Concentration (Metabolite/Cr)") + 
  labs(color= "Metabolite", shape = "Hemisphere") + 
  scale_shape_manual(values=c(16, 1)) + 
  scale_color_manual(values=c("#FF6666", "#66B2FF"), guide = guide_legend(reverse = TRUE)) +
  guides(linetype=FALSE) + 
  theme(legend.key = element_rect(fill = "white", colour = "black"), legend.title = element_text(size=rel(0.67)), legend.text = element_text(size=rel(0.67))) +
  ggtitle("DLPFC") +
  coord_cartesian(ylim=c(0, 2.25))
#ggsave(DLPFC_gabaglu, filename = "0518_DLPFC_Glu_GABA.pdf", width = 6, height = 6, units = "in", dpi = 300)

# without legend
DLPFC_gabaglu <- ggplot(gabaglu910) + aes(y=concentration, x=age, color=metabolite, linetype = label, shape=label) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Metabolite Level (/Cr)") + 
  ggtitle("DLPFC") +
  scale_shape_manual(values=c(16, 1)) + 
  scale_color_manual(values=c("#FF6666", "#66B2FF"), guide = guide_legend(reverse = TRUE)) +
  guides(linetype=FALSE) + 
  theme(legend.position = "none") + coord_cartesian(ylim=c(0, 2.10))
#ggsave(DLPFC_gabaglu, filename = "0908_DLPFC_Glu_GABA.pdf", width = 6, height = 6, units = "in", dpi = 300)


gaba7 <- ROI7_GABA %>% select(concentration=GABA.Cr, age, sex, GMrat) %>% mutate(metabolite="GABA")
glu7 <- ROI7_Glu %>% select(concentration=Glu.Cr, age, sex, GMrat) %>% mutate(metabolite="Glu")
gabaglu7 <- rbind(gaba7,glu7)
ACC_gabaglu <- ggplot(gabaglu7) + aes(y=concentration, x=age, color=metabolite) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Concentration (Metabolite/Cr)") + 
  labs(color= "Metabolite") + 
  ggtitle("ACC")+ 
  scale_shape_manual(values=c(16, 1)) + 
  scale_color_manual(values=c("#FF6666", "#66B2FF"), guide = guide_legend(reverse = TRUE)) +
  guides(linetype=FALSE) + 
  theme(legend.key = element_rect(fill = "white", colour = "black"), legend.title = element_text(size=rel(0.67)), legend.text = element_text(size=rel(0.67))) +
  coord_cartesian(ylim=c(0, 2.25))
#ggsave(ACC_gabaglu, filename = "0518_ACC_Glu_GABA.pdf", width = 6, height = 6, units = "in", dpi = 300)


# without legend
ACC_gabaglu <- ggplot(gabaglu7) + aes(y=concentration, x=age, color=metabolite) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Metabolite Level (/Cr)") + 
  ggtitle("ACC") +
  scale_shape_manual(values=c(16, 1)) + 
  scale_color_manual(values=c("#FF6666", "#66B2FF"), guide = guide_legend(reverse = TRUE)) +
  guides(linetype=FALSE) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim=c(0, 2.10))
#ggsave(ACC_gabaglu, filename = "0908_ACC_Glu_GABA.pdf", width = 6, height = 6, units = "in", dpi = 300)



gaba8 <- ROI8_GABA %>% select(concentration=GABA.Cr, age, sex, GMrat) %>% mutate(metabolite="GABA")
glu8 <- ROI8_Glu %>% select(concentration=Glu.Cr, age, sex, GMrat) %>% mutate(metabolite="Glu")
gabaglu8 <- rbind(gaba8,glu8)
MPFC_gabaglu <- ggplot(gabaglu8) + aes(y=concentration, x=age, color=metabolite) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Concentration (Metabolite/Cr)") + 
  labs(color= "Metabolite") + 
  ggtitle("MPFC") +
  scale_shape_manual(values=c(16, 1)) + 
  scale_color_manual(values=c("#FF6666", "#66B2FF"), guide = guide_legend(reverse = TRUE)) +
  guides(linetype=FALSE) + 
  theme(legend.key = element_rect(fill = "white", colour = "black"), legend.title = element_text(size=rel(0.67)), legend.text = element_text(size=rel(0.67))) +
  coord_cartesian(ylim=c(0, 2.10))
#ggsave(MPFC_gabaglu, filename = "0518_MPFC_Glu_GABA.png", width = 5, height = 5, units = "in", dpi = 300)

# without legend
MPFC_gabaglu <- ggplot(gabaglu8) + aes(y=concentration, x=age, color=metabolite) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Metabolite Level (/Cr)") + 
  ggtitle("MPFC") +
  scale_shape_manual(values=c(16, 1)) + 
  scale_color_manual(values=c("#FF6666", "#66B2FF"), guide = guide_legend(reverse = TRUE)) +
  guides(linetype=FALSE) + 
  theme(legend.position = "none") + coord_cartesian(ylim=c(0, 2.10))
#ggsave(MPFC_gabaglu, filename = "0908_MPFC_Glu_GABA.png", width = 5, height = 5, units = "in", dpi = 300)


#### Plot Ratio and Age ####
AI_ratio <- ggplot(ROI12_Ratio) + aes(y=Ratio, x=age, linetype = label, shape=label) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Glu/GABA Level") + 
  ggtitle("AIns") +
  theme(legend.position = "none") + 
  scale_shape_manual(values=c(16, 1)) +   coord_cartesian(ylim=c(1.75, 5))
#ggsave(AI_ratio, filename = "0908_AI_ratio.pdf", width = 6, height = 6, units = "in", dpi = 300)

DLPFC_ratio <- ggplot(ROI910_Ratio) + aes(y=Ratio, x=age, linetype = label, shape=label) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Glu/GABA Level") + 
  labs(shape = "Hemisphere") + 
  ggtitle("DLPFC")+
  scale_shape_manual(values=c(16, 1)) + 
  theme(legend.position = "none") + 
  coord_cartesian(ylim=c(1.75, 5))
#ggsave(DLPFC_ratio, filename = "0908_DLPFC_Ratio.pdf", width = 5, height = 5, units = "in", dpi = 300)


ACC_ratio <- ggplot(ROI7_Ratio) + aes(y=Ratio, x=age) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Glu/GABA Level") + 
  ggtitle("ACC") + 
  theme(legend.key = element_rect(fill = "white", colour = "black"), legend.title = element_text(size=rel(0.67)), legend.text = element_text(size=rel(0.67))) +
  coord_cartesian(ylim=c(1.75, 5))
#ggsave(ACC_ratio, filename = "0908_ACC_Ratio.pdf", width = 5, height = 5, units = "in", dpi = 300)

MPFC_ratio <- ggplot(ROI8_Ratio) + aes(y=Ratio, x=age) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Glu/GABA Level") + 
  ggtitle("MPFC")+
  theme(legend.key = element_rect(fill = "white", colour = "black"), legend.title = element_text(size=rel(0.67)), legend.text = element_text(size=rel(0.67))) +
  coord_cartesian(ylim=c(1.75, 5))
#ggsave(MPFC_ratio, filename = "0908_MPFC_Ratio.pdf", width = 5, height = 5, units = "in", dpi = 300)


#### Plot correlations as bar graphs ####
# make dataframe that has only good glutamate and gaba data for all ROIs
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

region_look <- list("1"="RAntInsula",
                    "2"="LAntInsula",
                    "7"='ACC',
                    "8"='MPFC',
                    "9"="RDLPFC",
                    "10"="LDLPFC")
keep_rois <- c(1,2,7,8,9,10)   # or keep_rois <- names(region_look)

region_avg <- MRS_corr %>%
  filter(roi %in% keep_rois) %>%          # remove any rois not in region_loopup
  select(ld8, agegrp,Glu.Cr,GABA.Cr, roi) %>%
  mutate(region=unlist(region_look[as.character(roi)])) %>%
  filter(!is.na(region)) %>%
  group_by(region,ld8, agegrp) %>%
  summarise_at(vars(Glu.Cr, GABA.Cr), mean)

cor_vals <- region_avg %>%
  group_by(region,agegrp) %>%
  summarise(GabaGlu_r=cor(Glu.Cr, GABA.Cr), 
            GabaGlu_p=cor.test(Glu.Cr, GABA.Cr, method=c("pearson"), use = "complete.obs")$p.value, 
            GabaGlu_lb=cor.test(Glu.Cr, GABA.Cr, method=c("pearson"), use = "complete.obs", conf.level=.95)$conf.int[1], 
            GabaGlu_ub=cor.test(Glu.Cr, GABA.Cr, method=c("pearson"), use = "complete.obs", conf.level=.95)$conf.int[2], 
            n=n())

corr_bar <- ggplot(cor_vals) +
  aes(x=agegrp, y=GabaGlu_r, fill=agegrp) +
  geom_bar(stat="identity") +
  facet_wrap(~region) +
  theme_classic(base_size = 15) + labs(x='Age',y='Glu GABA Corr (r)', fill = "Age Group") +
  scale_fill_brewer(palette = "viridis")

#ggsave(corr_bar, filename = "0518_corr_bargraph.pdf", width = 8, height = 6, units = "in", dpi = 300)

# comparing correlations using cocor package 

d_subset <- cor_vals %>% select(agegrp, r=GabaGlu_r, n)
d_pvals <- merge(d_subset,d_subset,by=NULL) %>% 
  mutate(pval=Vectorize(cocor::cocor.indep.groups)(r1.jk=r.x, r2.hm=r.y, n1=n.x, n2=n.y,
                                                   alternative="two.sided", alpha=0.05) %>% 
           sapply(function(x) x@fisher1925$p.value)
  )

p.adjust(d_pvals$pval, method="fdr") < .05
d_pvals %>% mutate(p_adjust = p.adjust(pval, method="fdr")) %>% filter(p_adjust < .05)

#### do it again but residualize out covariates(GMrat, sex) ####
MRS <- read.csv("/Users/mariaperica/Desktop/Lab/Projects/2020_MRSMGS/gmrat_13MP20200207_LCMv2fixidx.csv")

# Step 1 Outlier Detection - visual inspection of LCModel fits/spectra
# create a list of who to remove and remove them
lcm <- read_excel("/Users/mariaperica/Desktop/Lab/Projects/2020_MRSMGS/lcm.xlsx", col_names = FALSE)
lcm <- separate(lcm, "...1", c("ld8", "junk","y","x"),extra="merge", sep = "[-.]")
lcm <- dplyr::select(lcm, -junk)
lcm$bad <- TRUE
MRS <- MRS %>% mutate(x=216+1-x,y=216+1-y)
MRS <- merge(MRS, lcm, by=c("ld8", "x", "y"), all=T) 
MRS <- filter(MRS, is.na(bad))
MRS <- dplyr::select(MRS, -bad)

#keep only visit 1 people
MRS <- MRS %>% filter(visitnum==1)
#keep people's correct coordinates
MRS <- MRS %>% filter(!is.na(roi))
#get rid of people who are actually visit 2 but for some reason aren't filtered out
MRS <- MRS %>% filter(ld8!="10195_20191205")

# get rid of junk data noticed recently 
MRS<- MRS %>% filter(Glu.Cr != 0)

# save out a dataframe to share data after this step

# Step 2 Outlier Detection - get rid of peole who have bad data for 3 major metabolite peaks - GPC+Cho, NAA+NAAG, Cr
MRS<- filter(MRS, GPC.Cho.SD <= 10 | is.na(GPC.Cho.SD))
MRS <- filter(MRS, NAA.NAAG.SD <= 10 | is.na(NAA.NAAG.SD))
MRS <- filter(MRS, Cr.SD <= 10 | is.na(Cr.SD))

# Step 3 Outlier Detection - get rid of people who have lots of macromolecule in their spectra, as that can create distortions
MRS <- filter(MRS, MM20.Cr <= 3 | is.na(MM20.Cr))

#make inverse age column
MRS$invage <- 1/MRS$age
#make age^2 column
MRS$age2 <- (MRS$age - mean(MRS$age))^2

z_thres = 2



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


glu.resids <- lm(data=MRS_corr, Glu.Cr ~ sex + GMrat, na.action="na.exclude")
MRS_corr$glu.resids <- residuals(glu.resids)

gaba.resids <- lm(data=MRS_corr, GABA.Cr ~ sex + GMrat, na.action="na.exclude")
MRS_corr$gaba.resids <- residuals(gaba.resids)


region_look <- list("1"="RAntInsula",
                    "2"="LAntInsula",
                    "7"='ACC',
                    "8"='MPFC',
                    "9"="RDLPFC",
                    "10"="LDLPFC")
keep_rois <- c(1,2,7,8,9,10)   # or keep_rois <- names(region_look)

region_avg <- MRS_corr %>%
  filter(roi %in% keep_rois) %>%          # remove any rois not in region_loopup
  dplyr::select(ld8, agegrp,glu.resids,gaba.resids, roi) %>%
  mutate(region=unlist(region_look[as.character(roi)])) %>%
  filter(!is.na(region)) %>%
  group_by(region,ld8, agegrp) %>%
  summarise_at(vars(glu.resids, gaba.resids), mean)


cor_vals <- region_avg %>%
  group_by(region,agegrp) %>%
  summarise(GabaGlu_r=cor(glu.resids, gaba.resids), 
            GabaGlu_p=cor.test(glu.resids, gaba.resids, method=c("pearson"), use = "complete.obs")$p.value, 
            GabaGlu_lb=cor.test(glu.resids, gaba.resids, method=c("pearson"), use = "complete.obs", conf.level=.95)$conf.int[1], 
            GabaGlu_ub=cor.test(glu.resids, gaba.resids, method=c("pearson"), use = "complete.obs", conf.level=.95)$conf.int[2], 
            n=n())

ggplot(cor_vals) +
  aes(x=agegrp, y=GabaGlu_r, fill=agegrp) +
  geom_bar(stat="identity") +
  facet_wrap(~region) +
  theme_classic(base_size = 15) + labs(x='Age',y='Glu GABA Corr (r)', fill = "Age Group") +
  scale_fill_brewer(palette = "viridis")

#ggsave(corr_bar, filename = "0518_corr_bargraph.pdf", width = 8, height = 6, units = "in", dpi = 300)

# comparing correlations using cocor package 

d_subset <- cor_vals %>% select(agegrp, r=GabaGlu_r, n)
d_pvals <- merge(d_subset,d_subset,by=NULL) %>% 
  mutate(pval=Vectorize(cocor::cocor.indep.groups)(r1.jk=r.x, r2.hm=r.y, n1=n.x, n2=n.y,
                                                   alternative="two.sided", alpha=0.05) %>% 
           sapply(function(x) x@fisher1925$p.value)
  )

d_pvals_within_only <- d_pvals %>% filter(region.x == region.y)


corr_bar <- ggplot(cor_vals) +
  aes(x=agegrp, y=GabaGlu_r, fill=agegrp) +
  geom_bar(stat="identity") +
  facet_wrap(~region) +
  theme_classic(base_size = 15) + labs(x='Age',y='Glu GABA Corr (r)', fill = "Age Group") +
  scale_fill_brewer(palette = "Dark2") + coord_cartesian(ylim=c(0, 1.00))

#ggsave(corr_bar, filename = "1213_corr_bargraph.pdf", width = 8, height = 6, units = "in", dpi = 300)



#### Plot correlations as line graphs ####
MRS1 <- MRS_corr %>% filter(roi==1)
MRS2 <- MRS_corr %>% filter(roi==2)
MRS7 <- MRS_corr %>% filter(roi==7)
MRS8 <- MRS_corr %>% filter(roi==8)
MRS9 <- MRS_corr %>% filter(roi==9)
MRS10 <- MRS_corr %>% filter(roi==10)

MRS1 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% split(.$agegrp) %>% sapply(function(x) cor(x$Glu.Cr, x$GABA.Cr))
MRS1 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% group_by(agegrp) %>% tally
MRS1 <- MRS1 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age))
mrs1<- MRS1 %>% 
  mutate(agegrp = fct_recode(agegrp, "10-16" = "(0,16]","17-22" = "(16,22]", "23-30" = "(22,Inf]")) %>% 
  ggplot() + aes(y=GABA.Cr, x=Glu.Cr) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 18) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) + 
  facet_wrap(~ agegrp) +
  xlab("Glu/Cr") +
  ylab("GABA/Cr")
#ggsave(mrs1, filename = "0518_MRS1corr_scatter.pdf", width = 9, height = 5, units = "in", dpi = 300)


MRS2 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% split(.$agegrp) %>% sapply(function(x) cor(x$Glu.Cr, x$GABA.Cr))
MRS2 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% group_by(agegrp) %>% tally
MRS2 <- MRS2 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age))
mrs2<- MRS2 %>% 
  mutate(agegrp = fct_recode(agegrp, "10-16" = "(0,16]","17-22" = "(16,22]", "23-30" = "(22,Inf]")) %>% 
  ggplot() + aes(y=GABA.Cr, x=Glu.Cr) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 18) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) + 
  facet_wrap(~ agegrp) +
  xlab("Glu/Cr") +
  ylab("GABA/Cr")
#ggsave(mrs2, filename = "0518_MRS2corr_scatter.pdf", width = 9, height = 5, units = "in", dpi = 300)


MRS7 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% split(.$agegrp) %>% sapply(function(x) cor(x$Glu.Cr, x$GABA.Cr))
MRS7 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% group_by(agegrp) %>% tally
MRS7<- MRS7 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age))


mrs7<- MRS7 %>% 
  mutate(agegrp = fct_recode(agegrp, "10-16" = "(0,16]","17-22" = "(16,22]", "23-30" = "(22,Inf]")) %>% 
  ggplot() + aes(y=GABA.Cr, x=Glu.Cr) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 18) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) + 
  facet_wrap(~ agegrp) +
  xlab("Glu/Cr") +
  ylab("GABA/Cr")
#ggsave(mrs7, filename = "0518_MRS7corr_scatter.pdf", width = 9, height = 5, units = "in", dpi = 300)


MRS8 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% split(.$agegrp) %>% sapply(function(x) cor(x$Glu.Cr, x$GABA.Cr))
MRS8 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% group_by(agegrp) %>% tally
MRS8<- MRS8 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age))


ggplot(MRS8) + aes(y=GABA.Cr, x=Glu.Cr) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 18) +xlab("Glu.Cr") + ylab("GABA.Cr") + 
  theme(legend.key = element_rect(fill = "white", colour = "black")) + 
  facet_wrap(~ agegrp) +
  xlab("Glu/Cr") +
  ylab("GABA/Cr")

mrs8<- MRS8 %>% 
  mutate(agegrp = fct_recode(agegrp, "10-16" = "(0,16]","17-22" = "(16,22]", "23-30" = "(22,Inf]")) %>% 
  ggplot() + aes(y=GABA.Cr, x=Glu.Cr) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 18) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) + 
  facet_wrap(~ agegrp) +
  xlab("Glu/Cr") +
  ylab("GABA/Cr")
#ggsave(mrs8, filename = "0518_MRS8corr_scatter.pdf", width = 9, height = 5, units = "in", dpi = 300)


MRS9 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% split(.$agegrp) %>% sapply(function(x) cor(x$Glu.Cr, x$GABA.Cr))
MRS9 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% group_by(agegrp) %>% tally
MRS9<- MRS9 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age))

mrs9<- MRS9 %>% 
  mutate(agegrp = fct_recode(agegrp, "10-16" = "(0,16]","17-22" = "(16,22]", "23-30" = "(22,Inf]")) %>% 
  ggplot() + aes(y=GABA.Cr, x=Glu.Cr) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 18) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) + 
  facet_wrap(~ agegrp) +
  xlab("Glu/Cr") +
  ylab("GABA/Cr")
#ggsave(mrs9, filename = "0518_MRS9corr_scatter.pdf", width = 9, height = 5, units = "in", dpi = 300)

MRS10 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% split(.$agegrp) %>% sapply(function(x) cor(x$Glu.Cr, x$GABA.Cr))
MRS10 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% group_by(agegrp) %>% tally
MRS10<- MRS10 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age))
mrs10<- MRS10 %>% 
  mutate(agegrp = fct_recode(agegrp, "10-16" = "(0,16]","17-22" = "(16,22]", "23-30" = "(22,Inf]")) %>% 
  ggplot() + aes(y=GABA.Cr, x=Glu.Cr) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 18) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) + 
  facet_wrap(~ agegrp) +
  xlab("Glu/Cr") +
  ylab("GABA/Cr")
#ggsave(mrs10, filename = "0518_MRS10corr_scatter.pdf", width = 9, height = 5, units = "in", dpi = 300)


#### Correlation matrices ####

# used cleaned data 
MRS_glu <- MRS_glu %>%
  group_by(roi) %>%
  mutate(zscore=scale(Glu.Cr, center=T, scale=T)) %>%
  filter(abs(zscore) <= z_thres) %>% ungroup %>% mutate(agegrp = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30")))


MRS_GABA <- MRS_GABA %>%
  group_by(roi) %>%
  mutate(zscore=scale(GABA.Cr, center=T, scale=T)) %>%
  filter(abs(zscore) <= z_thres) %>% ungroup %>% mutate(agegrp = cut(age, breaks=c(0,16,22,Inf), labels=c("10-16","17-22", "23-30")))

# How correlated is Glutamate between regions?

# Without age residualized out 
# Step 1 - take good Glu data in MRS_glu dataframe and change the format
MRS_glu_wide <- pivot_wider(MRS_glu, id_cols=ld8, names_from=label, values_from=Glu.Cr)
MRS_glu_wide <- select(MRS_glu_wide, -`R STS`, -`L STS`, -`R Thalamus`, -`L Caudate`, -`R Caudate`, -`L Posterior Insula`, -`R Posterior Insula`)
MRS_glu_wide <- MRS_glu_wide %>% relocate('ld8', 'ACC', 'MPFC', 'R Anterior Insula', 'L Anterior Insula', 'R DLPFC', 'L DLPFC')

#Step 2 - Make correlation matrix
Glu_corr<- cor(MRS_glu_wide %>% select_if(is.numeric), use="pairwise.complete.obs")

#Step 3 - Make significance matrix 
pvals <- cor.mtest(MRS_glu_wide %>% select_if(is.numeric), conf.level = .95)

#Step 4 - plot correlation matrix 
corrplot(Glu_corr, p.mat = pvals$p, sig.level = .05, method = "number", type="lower", tl.col = "black", tl.srt = 45, insig = "blank", col = rev(brewer.pal(n = 8, name = "Spectral")))

#ggsave(corrmat1, filename = "0519_GluCorrMat.pdf", width = 6, height = 6, units = "in", dpi = 300)
# issue with the ggsave; need to resave these out as high res figs; says it cant save a matrix

# With age residualized out 

# Step 1 - Residualize age out of the correlation
mk_age_resid_glu <- function(d) {d[,'Glu.Cr'] <- lm(Glu.Cr ~ invage, d)$residuals; return(d) }
resids_glu <- MRS_glu %>% split(MRS_glu$label) %>% lapply(mk_age_resid_glu)
resids_glu <- resids_glu %>% bind_rows
# Step 1 - take good Glu data in MRS_glu dataframe and change the format
MRS_glu_wide_resids <- pivot_wider(resids_glu, id_cols=ld8, names_from=label, values_from=Glu.Cr)
MRS_glu_wide_resids <- select(MRS_glu_wide_resids, -`R STS`, -`L STS`, -`R Thalamus`, -`L Caudate`, -`R Caudate`, -`L Posterior Insula`, -`R Posterior Insula`)
MRS_glu_wide_resids <- MRS_glu_wide_resids %>% relocate('ld8', 'ACC', 'MPFC', 'R Anterior Insula', 'L Anterior Insula', 'R DLPFC', 'L DLPFC')

#Step 2 - Make correlation matrix
Glu_corr<- cor(MRS_glu_wide_resids %>% select_if(is.numeric), use="pairwise.complete.obs")

#Step 3 - Make significance matrix 
pvals <- cor.mtest(MRS_glu_wide_resids %>% select_if(is.numeric), conf.level = .95)

#Step 4 - plot correlation matrix 
corrplot(Glu_corr, p.mat = pvals$p, sig.level = .05, method = "number", type="lower", tl.col = "black", tl.srt = 45, insig = "blank", col = rev(brewer.pal(n = 8, name = "Spectral")))


# How correlated is GABA between regions?

# without age residualized out 
MRS_GABA_wide <- pivot_wider(MRS_GABA, id_cols=ld8, names_from=label, values_from=GABA.Cr)
MRS_GABA_wide <- select(MRS_GABA_wide, -`R STS`, -`L STS`, -`R Thalamus`, -`L Caudate`, -`R Caudate`, -`L Posterior Insula`, -`R Posterior Insula`)
MRS_GABA_wide <- MRS_GABA_wide %>% relocate('ld8', 'ACC', 'MPFC', 'R Anterior Insula', 'L Anterior Insula', 'R DLPFC', 'L DLPFC')

#Step 2 - Make correlation matrix
GABA_corr<- cor(MRS_GABA_wide %>% select_if(is.numeric), use="pairwise.complete.obs")

#Step 3 - Make significance matrix 
pvals <- cor.mtest(MRS_GABA_wide %>% select_if(is.numeric), conf.level = .95)

#Step 4 - plot correlation matrix 
corrplot(GABA_corr, p.mat = pvals$p, sig.level = .05, method = "number", type="lower", tl.col = "black", tl.srt = 45, insig = "blank", col = rev(brewer.pal(n = 8, name = "Spectral")))

# with age residualized out 
# Step 1 - take good GABA data in MRS_GABA dataframe and change the format
mk_age_resid_GABA <- function(d) {d[,'GABA.Cr'] <- lm(GABA.Cr ~ invage, d)$residuals; return(d) }
resids_GABA <- MRS_GABA %>% split(MRS_GABA$label) %>% lapply(mk_age_resid_GABA)
resids_GABA <- resids_GABA %>% bind_rows

MRS_GABA_wide_resids <- pivot_wider(resids_GABA, id_cols=ld8, names_from=label, values_from=GABA.Cr)
MRS_GABA_wide_resids <- select(MRS_GABA_wide_resids, -`R STS`, -`L STS`, -`R Thalamus`, -`L Caudate`, -`R Caudate`, -`L Posterior Insula`, -`R Posterior Insula`)
MRS_GABA_wide_resids <- MRS_GABA_wide_resids %>% relocate('ld8', 'ACC', 'MPFC', 'R Anterior Insula', 'L Anterior Insula', 'R DLPFC', 'L DLPFC')

#Step 2 - Make correlation matrix
GABA_corr<- cor(MRS_GABA_wide_resids %>% select_if(is.numeric), use="pairwise.complete.obs")

#Step 3 - Make significance matrix 
pvals <- cor.mtest(MRS_GABA_wide_resids %>% select_if(is.numeric), conf.level = .95)

#Step 4 - plot correlation matrix 
corrplot(GABA_corr, p.mat = pvals$p, sig.level = .05, method = "number", type="lower", tl.col = "black", tl.srt = 45, insig = "blank", col = rev(brewer.pal(n = 8, name = "Spectral")))


# How correlated are these between regions and among age groups

MRS_glu_young <- MRS_glu %>% filter(agegrp == "10-16")
MRS_glu_mid <- MRS_glu %>% filter(agegrp == "17-22")
MRS_glu_old <- MRS_glu %>% filter(agegrp == "23-30")

# young glu
MRS_glu_wide <- pivot_wider(MRS_glu_young, id_cols=ld8, names_from=label, values_from=Glu.Cr)
MRS_glu_wide <- select(MRS_glu_wide, -`R STS`, -`L STS`, -`R Thalamus`, -`L Caudate`, -`R Caudate`, -`L Posterior Insula`, -`R Posterior Insula`)
MRS_glu_wide <- MRS_glu_wide %>% relocate('ld8', 'ACC', 'MPFC', 'R Anterior Insula', 'L Anterior Insula', 'R DLPFC', 'L DLPFC')

#Step 2 - Make correlation matrix
Glu_corr<- cor(MRS_glu_wide %>% select_if(is.numeric), use="pairwise.complete.obs")

#Step 3 - Make significance matrix 
pvals <- cor.mtest(MRS_glu_wide %>% select_if(is.numeric), conf.level = .95)

#Step 4 - plot correlation matrix 
corrplot(Glu_corr, p.mat = pvals$p, sig.level = .05, method = "number", type="lower", tl.col = "black", tl.srt = 45, insig = "blank", col = rev(brewer.pal(n = 8, name = "Spectral")))


# mid glu
MRS_glu_wide <- pivot_wider(MRS_glu_mid, id_cols=ld8, names_from=label, values_from=Glu.Cr)
MRS_glu_wide <- select(MRS_glu_wide, -`R STS`, -`L STS`, -`R Thalamus`, -`L Caudate`, -`R Caudate`, -`L Posterior Insula`, -`R Posterior Insula`)
MRS_glu_wide <- MRS_glu_wide %>% relocate('ld8', 'ACC', 'MPFC', 'R Anterior Insula', 'L Anterior Insula', 'R DLPFC', 'L DLPFC')

#Step 2 - Make correlation matrix
Glu_corr<- cor(MRS_glu_wide %>% select_if(is.numeric), use="pairwise.complete.obs")

#Step 3 - Make significance matrix 
pvals <- cor.mtest(MRS_glu_wide %>% select_if(is.numeric), conf.level = .95)

#Step 4 - plot correlation matrix 
corrplot(Glu_corr, p.mat = pvals$p, sig.level = .05, method = "number", type="lower", tl.col = "black", tl.srt = 45, insig = "blank", col = rev(brewer.pal(n = 8, name = "Spectral")))

# old glu
MRS_glu_wide <- pivot_wider(MRS_glu_old, id_cols=ld8, names_from=label, values_from=Glu.Cr)
MRS_glu_wide <- select(MRS_glu_wide, -`R STS`, -`L STS`, -`R Thalamus`, -`L Caudate`, -`R Caudate`, -`L Posterior Insula`, -`R Posterior Insula`)
MRS_glu_wide <- MRS_glu_wide %>% relocate('ld8', 'ACC', 'MPFC', 'R Anterior Insula', 'L Anterior Insula', 'R DLPFC', 'L DLPFC')

#Step 2 - Make correlation matrix
Glu_corr<- cor(MRS_glu_wide %>% select_if(is.numeric), use="pairwise.complete.obs")

#Step 3 - Make significance matrix 
pvals <- cor.mtest(MRS_glu_wide %>% select_if(is.numeric), conf.level = .95)

#Step 4 - plot correlation matrix 
corrplot(Glu_corr, p.mat = pvals$p, sig.level = .05, method = "number", type="lower", tl.col = "black", tl.srt = 45, insig = "blank", col = rev(brewer.pal(n = 8, name = "Spectral")))


MRS_GABA_young <- MRS_GABA %>% filter(agegrp == "10-16")
MRS_GABA_mid <- MRS_GABA %>% filter(agegrp == "17-22")
MRS_GABA_old <- MRS_GABA %>% filter(agegrp == "23-30")

# young GABA
MRS_GABA_wide <- pivot_wider(MRS_GABA_young, id_cols=ld8, names_from=label, values_from=GABA.Cr)
MRS_GABA_wide <- select(MRS_GABA_wide, -`R STS`, -`L STS`, -`R Thalamus`, -`L Caudate`, -`R Caudate`, -`L Posterior Insula`, -`R Posterior Insula`)
MRS_GABA_wide <- MRS_GABA_wide %>% relocate('ld8', 'ACC', 'MPFC', 'R Anterior Insula', 'L Anterior Insula', 'R DLPFC', 'L DLPFC')

#Step 2 - Make correlation matrix
GABA_corr<- cor(MRS_GABA_wide %>% select_if(is.numeric), use="pairwise.complete.obs")

#Step 3 - Make significance matrix 
pvals <- cor.mtest(MRS_GABA_wide %>% select_if(is.numeric), conf.level = .95)

#Step 4 - plot correlation matrix 
corrplot(GABA_corr, p.mat = pvals$p, sig.level = .05, method = "number", type="lower", tl.col = "black", tl.srt = 45, insig = "blank", col = rev(brewer.pal(n = 8, name = "Spectral")))

# mid GABA
MRS_GABA_wide <- pivot_wider(MRS_GABA_mid, id_cols=ld8, names_from=label, values_from=GABA.Cr)
MRS_GABA_wide <- select(MRS_GABA_wide, -`R STS`, -`L STS`, -`R Thalamus`, -`L Caudate`, -`R Caudate`, -`L Posterior Insula`, -`R Posterior Insula`)
MRS_GABA_wide <- MRS_GABA_wide %>% relocate('ld8', 'ACC', 'MPFC', 'R Anterior Insula', 'L Anterior Insula', 'R DLPFC', 'L DLPFC')

#Step 2 - Make correlation matrix
GABA_corr<- cor(MRS_GABA_wide %>% select_if(is.numeric), use="pairwise.complete.obs")

#Step 3 - Make significance matrix 
pvals <- cor.mtest(MRS_GABA_wide %>% select_if(is.numeric), conf.level = .95)

#Step 4 - plot correlation matrix 
corrplot(GABA_corr, p.mat = pvals$p, sig.level = .05, method = "number", type="lower", tl.col = "black", tl.srt = 45, insig = "blank", col = rev(brewer.pal(n = 8, name = "Spectral")))

#old gaba
MRS_GABA_wide <- pivot_wider(MRS_GABA_old, id_cols=ld8, names_from=label, values_from=GABA.Cr)
MRS_GABA_wide <- select(MRS_GABA_wide, -`R STS`, -`L STS`, -`R Thalamus`, -`L Caudate`, -`R Caudate`, -`L Posterior Insula`, -`R Posterior Insula`)
MRS_GABA_wide <- MRS_GABA_wide %>% relocate('ld8', 'ACC', 'MPFC', 'R Anterior Insula', 'L Anterior Insula', 'R DLPFC', 'L DLPFC')

#Step 2 - Make correlation matrix
GABA_corr<- cor(MRS_GABA_wide %>% select_if(is.numeric), use="pairwise.complete.obs")

#Step 3 - Make significance matrix 
pvals <- cor.mtest(MRS_GABA_wide %>% select_if(is.numeric), conf.level = .95)

#Step 4 - plot correlation matrix 
corrplot(GABA_corr, p.mat = pvals$p, sig.level = .05, method = "number", type="lower", tl.col = "black", tl.srt = 45, insig = "blank", col = rev(brewer.pal(n = 8, name = "Spectral")))

#### gaba ~ glu * age ####
MRS1$z_Glu <- scale(MRS1_Glu$Glu.Cr)
MRS1$z_GABA <- scale(MRS1$GABA.Cr)
MRS1$z_invage <- scale(MRS1$invage)
MRS1$z_gm <- scale(MRS1$GMrat)

MRS1_gabaglu <- lm(data=MRS1, z_GABA ~ z_Glu)
summary(MRS1_gabaglu)
MRS1_gabaglu_full <- lm(data=MRS1, z_GABA ~ z_Glu * z_invage + z_gm + sex)
summary(MRS1_gabaglu_full)

MRS1_gabaglu <- lm(data=MRS1, GABA.Cr ~ Glu.Cr * invage)
summary(MRS1_gabaglu)
MRS1_gabaglu <- lm(data=MRS1, z_GABA ~ z_Glu * z_invage)
summary(MRS1_gabaglu)

MRS1_gabaglu <- lm(data=MRS1, GABA.Cr ~ Glu.Cr * invage + GMrat + sex)
summary(MRS1_gabaglu)

ggplot(data=MRS1) + aes(x=Glu.Cr, y=GABA.Cr, color=agegrp) + geom_point() + geom_smooth(method="lm")

RAI_balance <- ggplot(MRS1) + aes(y=Glu.Cr, x=GABA.Cr, color=agegrp) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 15) +xlab("Glu/Cr") + ylab("GABA/Cr") + 
  labs(shape = "Age Group (years)") + 
  ggtitle("Right AIns") +
  scale_shape_manual(values=c(16, 1)) + 
  theme(legend.position = "none")
ggsave(RAI_balance, filename = "0908_RAI_balance.pdf", width = 6, height = 6, units = "in", dpi = 300)

# trying with finn code
thisroi <- 'DLPFC'
this <- MRS_corr %>% filter(grepl(thisroi, label)) %>% mutate(subj = as.character(ld8)) %>% 
  select(subj, ld8, GABA.Cr, Glu.Cr, label, age, invage, GMrat, sex) %>% filter(complete.cases(.))
lm.model.interaction <- lmer(GABA.Cr ~ Glu.Cr*invage + label + GMrat + sex + (1|ld8), data = this)
#lm.model.interaction <- lm(GABA.Cr ~ Glu.Cr*age, data = this)
library(interactions)
interact_plot(lm.model.interaction, pred = Glu.Cr, modx = invage, modx.values = c(12,18,24), main.title = thisroi)
summary(lm.model.interaction)


MRS2$z_Glu <- scale(MRS2$Glu.Cr)
MRS2$z_GABA <- scale(MRS2$GABA.Cr)
MRS2$z_invage <- scale(MRS2$invage)
MRS2$z_gm <- scale(MRS2$GMrat)


mrs2_gabaglu <- lm(data=MRS2, z_GABA ~ z_Glu)
summary(mrs2_gabaglu)

mrs2_gabaglu_full <- lm(data=MRS2, z_GABA ~ z_Glu + z_invage + z_gm + sex)
summary(mrs2_gabaglu_full)
mrs2_gabaglu_full <- lm(data=MRS2, z_GABA ~ z_Glu * z_invage + z_gm + sex)
summary(mrs2_gabaglu_full)


mrs2_gabaglu <- lm(data=MRS2, GABA.Cr ~ Glu.Cr * invage)
summary(mrs2_gabaglu)

mrs2_gabaglu <- lm(data=MRS2, scale(GABA.Cr) ~ scale(Glu.Cr) * scale(invage))
summary(mrs2_gabaglu)

mrs2_gabaglu <- lm(data=MRS2, GABA.Cr ~ Glu.Cr * invage + GMrat + sex)
summary(mrs2_gabaglu)

ggplot(data=MRS2) + aes(x=Glu.Cr, y=GABA.Cr, color=agegrp) + geom_point() + geom_smooth(method="lm")
LAI_balance <- ggplot(MRS2) + aes(y=Glu.Cr, x=GABA.Cr, color=agegrp) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 15) +xlab("Glu/Cr") + ylab("GABA/Cr") + 
  labs(shape = "Age Group (years)") + 
  ggtitle("Left AIns") +
  scale_shape_manual(values=c(16, 1)) + 
  theme(legend.position = "none")
ggsave(LAI_balance, filename = "0908_LAI_balance.pdf", width = 6, height = 6, units = "in", dpi = 300)


MRS7$z_Glu <- scale(MRS7$Glu.Cr)
MRS7$z_GABA <- scale(MRS7$GABA.Cr)
MRS7$z_invage <- scale(MRS7$invage)
MRS7$z_gm <- scale(MRS7$GMrat)


mrs7_gabaglu <- lm(data=MRS7, z_GABA ~ z_Glu)
summary(mrs7_gabaglu)

mrs7_gabaglu_full <- lm(data=MRS7, z_GABA ~ z_Glu + z_invage + z_gm + sex)
summary(mrs7_gabaglu_full)
mrs7_gabaglu_full <- lm(data=MRS7, z_GABA ~ z_Glu * z_invage + z_gm + sex)
summary(mrs7_gabaglu_full)

mrs7_gabaglu <- lm(data=MRS7, GABA.Cr ~ Glu.Cr * invage)
summary(mrs7_gabaglu)
mrs7_gabaglu <- lm(data=MRS7, GABA.Cr ~ Glu.Cr * invage + GMrat + sex)
summary(mrs7_gabaglu)

ggplot(data=MRS7) + aes(x=Glu.Cr, y=GABA.Cr, color=agegrp) + geom_point() + geom_smooth(method="lm")
ACC_balance <- ggplot(MRS7) + aes(y=Glu.Cr, x=GABA.Cr, color=agegrp) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 15) +xlab("Glu/Cr") + ylab("GABA/Cr") + 
  labs(shape = "Age Group (years)") + 
  ggtitle("ACC") +
  scale_shape_manual(values=c(16, 1)) + 
  theme(legend.position = "none")
ggsave(ACC_balance, filename = "0908_ACC_balance.pdf", width = 6, height = 6, units = "in", dpi = 300)


MRS8$z_Glu <- scale(MRS8$Glu.Cr)
MRS8$z_GABA <- scale(MRS8$GABA.Cr)
MRS8$z_invage <- scale(MRS8$invage)
MRS8$z_gm <- scale(MRS8$GMrat)


mrs8_gabaglu <- lm(data=MRS8, z_GABA ~ z_Glu)
summary(mrs8_gabaglu)

mrs8_gabaglu_full <- lm(data=MRS8, z_GABA ~ z_Glu + z_invage + z_gm + sex)
summary(mrs8_gabaglu_full)

mrs8_gabaglu_full <- lm(data=MRS8, z_GABA ~ z_Glu * z_invage + z_gm + sex)
summary(mrs8_gabaglu_full)

mrs8_gabaglu <- lm(data=MRS8, GABA.Cr ~ Glu.Cr * invage + GMrat + sex)
summary(mrs8_gabaglu)
ggplot(data=MRS8) + aes(x=Glu.Cr, y=GABA.Cr, color=agegrp) + geom_point() + geom_smooth(method="lm")


MPFC_balance <- ggplot(MRS8) + aes(y=Glu.Cr, x=GABA.Cr, color=agegrp) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 15) +xlab("Glu/Cr") + ylab("GABA/Cr") + 
  labs(shape = "Age Group (years)") + 
  ggtitle("MPFC") +
  scale_shape_manual(values=c(16, 1)) + 
  theme(legend.position = "none")
ggsave(MPFC_balance, filename = "0908_MPFC_balance.pdf", width = 6, height = 6, units = "in", dpi = 300)

MRS9$z_Glu <- scale(MRS9$Glu.Cr)
MRS9$z_GABA <- scale(MRS9$GABA.Cr)
MRS9$z_invage <- scale(MRS9$invage)
MRS9$z_gm <- scale(MRS9$GMrat)


mrs9_gabaglu <- lm(data=MRS9, z_GABA ~ z_Glu)
summary(mrs9_gabaglu)

mrs9_gabaglu_full <- lm(data=MRS9, z_GABA ~ z_Glu + z_invage + z_gm + sex)
summary(mrs9_gabaglu_full)
mrs9_gabaglu_full <- lm(data=MRS9, z_GABA ~ z_Glu * z_invage + z_gm + sex)
summary(mrs9_gabaglu_full)

mrs9_gabaglu <- lm(data=MRS9, GABA.Cr ~ Glu.Cr * invage + GMrat + sex)
summary(mrs9_gabaglu)
ggplot(data=MRS9) + aes(x=Glu.Cr, y=GABA.Cr, color=agegrp) + geom_point() + geom_smooth(method="lm")

RDLPFC_balance <- ggplot(MRS9) + aes(y=Glu.Cr, x=GABA.Cr, color=agegrp) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 15) +xlab("Glu/Cr") + ylab("GABA/Cr") + 
  labs(shape = "Age Group (years)") + 
  ggtitle("Right DLPFC") +
  scale_shape_manual(values=c(16, 1)) + 
  theme(legend.position = "none")
ggsave(RDLPFC_balance, filename = "0908_RDLPFC_balance.pdf", width = 6, height = 6, units = "in", dpi = 300)


MRS10$z_Glu <- scale(MRS10$Glu.Cr)
MRS10$z_GABA <- scale(MRS10$GABA.Cr)
MRS10$z_invage <- scale(MRS10$invage)
MRS10$z_gm <- scale(MRS10$GMrat)


mrs10_gabaglu <- lm(data=MRS10, z_GABA ~ z_Glu)
summary(mrs10_gabaglu)

mrs10_gabaglu_full <- lm(data=MRS10, z_GABA ~ z_Glu + z_invage + z_gm + sex)
summary(mrs10_gabaglu_full)
mrs10_gabaglu_full <- lm(data=MRS10, z_GABA ~ z_Glu * z_invage + z_gm + sex)
summary(mrs10_gabaglu_full)

mrs10_gabaglu <- lm(data=MRS10, GABA.Cr ~ Glu.Cr * invage + GMrat + sex)
summary(mrs10_gabaglu)
ggplot(data=MRS10) + aes(x=Glu.Cr, y=GABA.Cr, color=agegrp) + geom_point() + geom_smooth(method="lm")

LDLPFC_balance <- ggplot(MRS10) + aes(y=Glu.Cr, x=GABA.Cr, color=agegrp) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 15) +xlab("Glu/Cr") + ylab("GABA/Cr") + 
  labs(shape = "Age Group (years)") + 
  ggtitle("Left DLPFC") +
  scale_shape_manual(values=c(16, 1)) + 
  theme(legend.position = "none")
ggsave(LDLPFC_balance, filename = "0908_LDLPFC_balance.pdf", width = 6, height = 6, units = "in", dpi = 300)




#### Variability residuals and balance residuals ####

# RAIns
# get residual of balance 
MRS1_gabaglu <- lm(data=MRS1, GABA.Cr ~ Glu.Cr * invage + GMrat + sex, na.action = na.exclude)
summary(MRS1_gabaglu)
balance_rsdl <- abs(residuals(MRS1_gabaglu))
MRS1$balance_resids <- balance_rsdl


#variability - get residual of glu onto inverse age 
var_Glu1 <- lm(data=MRS1, Glu.Cr ~ invage + GMrat + sex, na.action = na.exclude) 
summary(var_Glu1)
var_rsdl <- abs(residuals(var_Glu1))
MRS1$var_resids <- var_rsdl


balance_var1 <- lm(data=MRS1, balance_resids ~ var_resids + invage)
summary(balance_var1)

# ACC
# get residual of balance 
MRS7_gabaglu <- lm(data=MRS7, GABA.Cr ~ Glu.Cr * invage + GMrat + sex, na.action = na.exclude)
summary(MRS7_gabaglu)
balance_rsdl <- abs(residuals(MRS7_gabaglu))
MRS7$balance_resids <- balance_rsdl

#variability - get residual of glu onto inverse age 
var_Glu7 <- lm(data=MRS7, Glu.Cr ~ invage + GMrat + sex, na.action = na.exclude) 
summary(var_Glu7)
var_rsdl <- abs(residuals(var_Glu7))
MRS7$var_resids <- var_rsdl


balance_var7 <- lm(data=MRS7, balance_resids ~ var_resids +invage)
summary(balance_var7)

# MPFC

# get residual of balance 
MRS8_gabaglu <- lm(data=MRS8, GABA.Cr ~ Glu.Cr * invage + GMrat + sex, na.action = na.exclude)
summary(MRS8_gabaglu)
balance_rsdl <- abs(residuals(MRS8_gabaglu))
MRS8$balance_resids <- balance_rsdl

#variability - get residual of glu onto inverse age 
var_Glu8 <- lm(data=MRS8, Glu.Cr ~ invage + GMrat + sex, na.action = na.exclude) 
summary(var_Glu8)
var_rsdl <- abs(residuals(var_Glu8))
MRS8$var_resids <- var_rsdl

balance_var8 <- lm(data=MRS8, balance_resids ~ var_resids + invage)
summary(balance_var8)


#### sliding window ####
# sliding window - approach 2
rois <- c('DLPFC', 'L DLPFC', 'R DLPFC', 'MPFC', 'ACC', 'Anterior Insula', 'R Anterior Insula', 'L Anterior Insula')
allSlopes <- c()
for (thisroii in seq(1, length(rois))) {
  thisroi <- rois[thisroii]
  this <- MRS_corr %>% filter(grepl(thisroi, label)) %>% mutate(subj = as.character(ld8)) %>% 
    select(subj, ld8, GABA.Cr, Glu.Cr, label, age, invage, GMrat, sex) %>% filter(complete.cases(.)) %>% arrange(age)
  ages <- this %>% select(age) %>% distinct() %>% arrange(age)
  winSize <- 80
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
    
    thisEffect <- lm.model.summary$coefficients[2,1]
    meanAge <- unname(unlist(this %>% filter(age >= minAge & age <= maxAge) %>% select(age) %>% summarize(midAge = mean(age))))
    
    allSlopes <- rbind(allSlopes, data.frame(roi = thisroi, slope = thisEffect, minAge, maxAge, meanAge))
  }
  
}

LNCDR::lunaize(ggplot(data = allSlopes %>% filter(grepl(' DLPFC', roi)), aes(x = meanAge, y= slope, color=roi)) +
          geom_point() + geom_errorbarh(aes(xmin = minAge,xmax = maxAge), alpha=0.2) +
          stat_smooth(method='loess', span=1.5) + labs(x = 'Age', y = 'GABA-Glu Association', title='')) + 
  theme(legend.position=c(.8,1), legend.title = element_blank())

