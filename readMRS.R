#!/usr/bin/env Rscript
pacman::p_load(readxl, dplyr, tidyr)
select <- dplyr::select  # make sure we're using the select we mean too
source('mrsi_funcs.R') # mrsi_metqc, mrsi_add_cols

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
#get rid of people who are actually visit 2 but for some reason aren't filtered out
MRS <- MRS %>% filter(ld8!="10195_20191205")

# remove rows where any of GPC.Cho.SD, NAA.NAAG.SD, Cr.SD, MM20.Cr have very bad values
MRS <- MRS %>% mrsi_metqc()

# add invage, age2 (squared), region (roi labels), dateNumeric, agegrp
MRS <- MRS %>% mrsi_add_columns()
