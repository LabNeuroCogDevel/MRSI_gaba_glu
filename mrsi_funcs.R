require(dplyr)

#' reject rows that fail quality checking
#' 
#' @param MRS input dataframe with SD and Cr columns for macromolicules
#' @return dataframe with noisy rows removed
mrsi_metqc <- function(MRS) {
  checkmate::expect_subset(c("roi", "GPC.Cho.SD", "NAA.NAAG.SD", "Cr.SD", "MM20.Cr", "Glu.Cr"),
                           names(MRS))
  #keep people's correct coordinates
  MRS <- MRS %>% filter(!is.na(roi))
  # get rid of junk data noticed recently 
  MRS<- MRS %>% filter(Glu.Cr != 0)
  
  # save out a dataframe to share data after this step
  
  # Step 2 Outlier Detection - get rid of peole who have bad data for 3 major metabolite peaks - GPC+Cho, NAA+NAAG, Cr
  MRS <- MRS %>%
     filter(GPC.Cho.SD  <= 10 | is.na(GPC.Cho.SD),
            NAA.NAAG.SD <= 10 | is.na(NAA.NAAG.SD),
            Cr.SD       <= 10 | is.na(Cr.SD),
            # Step 3 Outlier Detection - get rid of people who have lots of macromolecule in their spectra, as that can create distortions
            MM20.Cr      <= 3 | is.na(MM20.Cr))
}


#' add calculated columns age,age2,agegrp,region,dateNumeric
#'
#' @param MRS input dataframe with age, roi, and ld8
#' @return dataframe with new columns
mrsi_add_cols <- function(MRS){
  checkmate::expect_subset(c("age","roi", "ld8"), names(MRS))
  #make inverse age column
  MRS$invage <- 1/MRS$age
  #make age^2 column
  MRS$age2 <- (MRS$age - mean(MRS$age))^2

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
   "11" = "RSTS",
   "12" = "LSTS",
   "13" = "RThalamus")
  
  MRS <- MRS %>%
      mutate(agegrp = cut(age, breaks=c(0,16,22,Inf),
                          labels=c("10-16","17-22", "23-30")),
            region=unlist(region_look[as.character(roi)]),
            dateNumeric=as.numeric(as.POSIXct(gsub('.*_','',ld8), format="%Y%m%d")))
}


#' row per metabolite, but keep Cr and SD each in their own column
#'
#' @param mrs  input dataframe with many column names ending in .SD and .Cr
#' @return longer dataframe with new columsn 'met' 'Cr' 'SD'
longer_met_CrSD <- function(mrs){
  mrs %>%
     tidyr::pivot_longer(cols=matches('\\.(SD|Cr)$'),
                         names_pattern="(\\S+)\\.(SD|Cr)",
                         names_to=c("met","mtype")) %>%
     tidyr::spread(mtype,value)
}
