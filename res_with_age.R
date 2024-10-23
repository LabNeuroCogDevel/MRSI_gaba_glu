require(dplyr)
require(mgcv)
#
# 20230530 WF copy of FC
# also see MultHemiGamm.R:met_res_for_roi  

#' adjust metabolte with residuals of gam model
#'
#' 20230531 - removed this_roi option. filter before sending to function
#' 20241023 - added save_model option for Pupi&Sarpal
#' @param MRSI_input dataframe with roi, age, dateNumeric
#' @param met_name  metabolite eg. 'GABA.Cre'. column in MRSI_input
#' @param return_df bool, if NULL (default) return df only when this_roi is not NULL
#' @param return_model bool return gam model instead of dataframe
#' @return vector of adjusted metabolite OR dataframe with new 'met_name'_gamadj column OR gam model
res_with_age <- function(MRSI_input, met_name, return_df=FALSE, return_model=FALSE) {
  checkmate::expect_subset(c("age", "dateNumeric", met_name), names(MRSI_input))

  model <- paste0(met_name, ' ~ s(age, k=3) + s(dateNumeric, k=5)')

  # 20231025 - GMrat not always wanted here. MP wants to report GM model fit in own model
  # but be careful!  from FC: 
  # > gm is correlated to age, so we want to be sure we're removing the effects of %gm and date
  # > and not age itself
  hasGM <- FALSE
  if("GMrat" %in% names(MRSI_input)) {
     model <- paste0(model ,'+ s(GMrat, k=3)')
     hasGM <- TRUE
  } else{
     warning("Not including GM when generating residuals with gam model! no 'GMrat' column")
  }

  # if region and there are multiples, e.g. if have e.g. L and R, model each
  regions <- paste0(collapse=',', unique(MRSI_input[['region']]))
  if("region" %in% names(MRSI_input) && length(unique(MRSI_input$region))>1L) {
     cat(glue::glue("modeling multiple regions: {regions}\n"),"\n")
     MRSI_input$region <- as.factor(MRSI_input$region) # ensure it's a factor for modeling
     model <- paste0(model, '+ region')
  }

  
  mrsi.gam <- tryCatch(mgcv::gam(as.formula(model), data=MRSI_input, na.action = na.exclude),
                       error=function(e){print(e);print(MRSI_input[1,]); return(NULL)})
  if(is.null(mrsi.gam)){
     met_adj <- rep(NA,nrow(MRSI_input))
  }
  else{
     # get residuals at average date and greymatter
     center <- MRSI_input %>%
        mutate(dateNumeric = mean(dateNumeric, na.rm=T),
               GMrat = mean(GMrat, na.rm=T))

     # met_adj is what's left after we've modeled are nuisance regressors
     met_adj <- unname(predict(mrsi.gam,center) + residuals(mrsi.gam))
  }

  if(return_model) return(mrsi.gam)

  # vector or dataframe
  if(!return_df) return(met_adj)

  # otherwise, make a new column
  new_col <- paste0(met_name,'_gamadj')
  MRSI_input[new_col] <- as.numeric(met_adj)
  return(MRSI_input)
}
