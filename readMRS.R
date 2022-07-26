#!/usr/bin/env Rscript
pacman::p_load(readxl, dplyr, tidyr)
select <- dplyr::select

# makes MRS variable from csv and xlsx
csv_file <- "13MP20200207_LCMv2fixidx.csv"
xlsx_file <- "lcm.xlsx"

MRS <- read.csv(csv_file)

# Step 1 Outlier Detection - visual inspection of LCModel fits/spectra
# create a list of who to remove and remove them
lcm <- read_excel(xlsx_file, col_names = FALSE)
lcm <- separate(lcm, "...1", c("ld8", "junk","y","x"),extra="merge", sep = "[-.]")
lcm <- select(lcm, -junk)
lcm$bad <- TRUE
MRS <- MRS %>% mutate(x=216+1-x,y=216+1-y)
MRS <- merge(MRS, lcm, by=c("ld8", "x", "y"), all=T) 
MRS <- filter(MRS, is.na(bad))
MRS <- select(MRS, -bad)

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


## 3. and agegrp and region label
##  NB. could remove R/L from names to later collapse/average across hemisphere
##  13MP20200207_LCMv2fixidx.csv now includes 'label' column
region_look <- list(
  "1" = "RAntInsula",
  "2" = "LAntInsula",
  "3" = "RPostInsula",
  "4" = "LPostInsula",
  "5" = "RCaudate",
  "6" = "LCaudate",
  "7" = "ACC",
  "8" = "MPFC",
  "9" = "RDLPFC",
 "10" = "LDLPFC",
 "11" = "R STS",
 "12" = "L STS",
 "13" = "R Thalamus")

MRS <- MRS %>%
    mutate(agegrp = cut(age, breaks=c(0,16,22,Inf),
                        labels=c("10-16","17-22", "23-30")),
          region=unlist(region_look[as.character(roi)]))
