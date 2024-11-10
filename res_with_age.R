require(dplyr)
require(mgcv)
#
# 20230530 WF copy of FC
# also see MultHemiGamm.R:met_res_for_roi  

#' adjust metabolte with residuals of gam model
#'
#' 20230531 - removed this_roi option. filter before sending to function
#' 20241023 - added save_model option for Pupi&Sarpal
#' 20241031 - center_on and min_age to model for Sarpal data
#' @param MRSI_input dataframe with roi, age, dateNumeric
#' @param met_name  metabolite eg. 'GABA.Cre'. column in MRSI_input
#' @param return_df bool, if NULL (default) return df only when this_roi is not NULL
#' @param return_model bool return gam model instead of dataframe
#' @param center_on vector of columns to replace with their center before predict-ing on generated model
#' @param include.na bool. model with NA data? default is TRUE but breaks viz packages 
#' @param min_age ages below this value are removed from modeling (20241031 for DS data)
#' @return vector of adjusted metabolite OR dataframe with new 'met_name'_gamadj column OR gam model
res_with_age <- function(MRSI_input, met_name, return_df=FALSE, return_model=FALSE,
                         center_on=c("dateNumeric","GMrat"),
                         include.na=TRUE,
                         min_age=-Inf) {
  checkmate::expect_subset(c("age", "dateNumeric", met_name), names(MRSI_input))

  model <- paste0(met_name, ' ~ s(age, k=3) + s(dateNumeric, k=5)')

  # 20231025 - GMrat not always wanted here. MP wants to report GM model fit in own model
  # but be careful!  from FC: 
  # > gm is correlated to age, so we want to be sure we're removing the effects of %gm and date
  # > and not age itself
  if("GMrat" %in% names(MRSI_input)) {
     model <- paste0(model ,'+ s(GMrat, k=3)')
  } else{
     # 20241031: center_on might inclue GM but it won't be applied if missing from the data
     warning("Not including GM when generating residuals with gam model! no 'GMrat' column")
  }

  # if region and there are multiples, e.g. if have e.g. L and R, model each
  regions <- paste0(collapse=',', unique(MRSI_input[['region']]))
  if("region" %in% names(MRSI_input) && length(unique(MRSI_input$region))>1L) {
     cat(glue::glue("modeling multiple regions: {regions}\n"),"\n")
     MRSI_input$region <- as.factor(MRSI_input$region) # ensure it's a factor for modeling
     model <- paste0(model, '+ region')
  }

  gam_formula <- as.formula(model)
  MRSI_restricted <- MRSI_input %>% filter(age > min_age)
  # 20241030: remove NA columns before going into model to work with viz packages (mgcViz)?
  if(!include.na) {
     have_all_data <- complete.cases(MRSI_restricted[,all.vars(gam_formula)])
     MRSI_restricted <- MRSI_restricted[have_all_data,] 
  }

  # 20241031: optionally remove e.g. less than 18yo
  mrsi.gam <- tryCatch(mgcv::gam(gam_formula, data=MRSI_restricted, na.action = na.exclude),
                       error=function(e){print(e);print(MRSI_input[1,]); return(NULL)})

  if(return_model) return(mrsi.gam)

  if(is.null(mrsi.gam)){
     met_adj <- rep(NA,nrow(MRSI_input))
  }
  else{
     # get residuals at average date and greymatter
     center <- MRSI_input %>%
        mutate(across(any_of(center_on), \(x) mean(x,na.rm=T)))
        # 20241031: maybe we want to keep GMrat in the prediction?
               #dateNumeric = mean(dateNumeric, na.rm=T),
               #GMrat = mean(GMrat, na.rm=T))

     # met_adj is what's left after we've modeled are nuisance regressors
     met_adj <- unname(predict(mrsi.gam,center) + residuals(mrsi.gam))
  }


  # vector or dataframe
  if(!return_df) return(met_adj)

  # otherwise, make a new column
  new_col <- paste0(met_name,'_gamadj')
  MRSI_input[new_col] <- as.numeric(met_adj)
  return(MRSI_input)
}






#' model MRSI with shared variance inputs visit date, age, and gray matter ratio
#'
#' @param MRSI_input dataframe with roi, age, dateNumeric
#' @param met_name  metabolite eg. 'GABA.Cre'. column in MRSI_input
#' @param include.na bool. model with NA data? default is TRUE but breaks viz packages 
#' @param min_age ages below this value are removed from modeling (20241031 for DS data)
#' @return model
mrsi_model_date <- function(MRSI_input, met_name, min_age=-Inf, include.na=TRUE) {
  checkmate::expect_subset(c("age", "dateNumeric", met_name), names(MRSI_input))

  model <- paste0(met_name, ' ~ s(age, k=3) + s(dateNumeric, k=5)')

  # 20231025 - GMrat not always wanted here. MP wants to report GM model fit in own model
  # but be careful!  from FC: 
  # > gm is correlated to age, so we want to be sure we're removing the effects of %gm and date
  # > and not age itself
  if("GMrat" %in% names(MRSI_input)) {
     model <- paste0(model ,'+ s(GMrat, k=3)')
  } else{
     # 20241031: center_on might inclue GM but it won't be applied if missing from the data
     warning("Not including GM when generating residuals with gam model! no 'GMrat' column")
  }

  # if region and there are multiples, e.g. if have e.g. L and R, model each
  regions <- paste0(collapse=',', unique(MRSI_input[['region']]))
  if("region" %in% names(MRSI_input) && length(unique(MRSI_input$region))>1L) {
     cat(glue::glue("modeling multiple regions: {regions}\n"),"\n")
     MRSI_input$region <- as.factor(MRSI_input$region) # ensure it's a factor for modeling
     model <- paste0(model, '+ region')
  }

  # 20241030: remove NA columns before going into model to work with viz packages (mgcViz)
  gam_formula <- as.formula(model)
  MRSI_restricted <- MRSI_input %>% filter(age > min_age)
  # 20241030: remove NA columns before going into model to work with viz packages (mgcViz)?
  if(!include.na) {
     have_all_data <- complete.cases(MRSI_restricted[,all.vars(gam_formula)])
     MRSI_restricted <- MRSI_restricted[have_all_data,] 
  }
  # 20241031: optionally remove e.g. less than 18yo
  mrsi.gam <- tryCatch(mgcv::gam(gam_formula, data=MRSI_restricted, na.action = na.exclude),
                       error=function(e){print(e);print(MRSI_input[1,]); return(NULL)})

  return(mrsi.gam)
}

#' adjust metabolte with residuals of gam model
#'
#' @param mrsi.gam model, likely from 
#' @param MRSI_input dataframe with roi, age, dateNumeric
#' @param return_df bool, if NULL (default) return df only when this_roi is not NULL
#' @param center_on vector of columns to replace with their center before predict-ing
#' @return vector of adjusted metabolite OR dataframe with new 'met_name'_gamadj column
adj_by_model <- function(mrsi.gam, MRSI_input,
                         return_df=FALSE,
                         center_on=c("dateNumeric","GMrat")) {
  checkmate::expect_subset(all.vars(mrsi.gam$formula), names(MRSI_input))
  checkmate::expect_subset(center_on, names(MRSI_input))
  met_name <- formula.tools::lhs(mrsi.gam$formula)

  if(is.null(mrsi.gam)){
     met_adj <- rep(NA,nrow(MRSI_input))
  }
  else{
     # residuals (difference from expectation given date, GMrat, age)
     yhat <- unname(predict(mrsi.gam, MRSI_input))
     e <- MRSI_input[[met_name]] - yhat

     # predict w/all row's date (and maybe GMrat) replaced with sample mean
     # importantly, using real (not-centered) age
     center <- MRSI_input %>%
        mutate(across(any_of(center_on), \(x) mean(x,na.rm=T)))
     center_predict <- unname(predict(mrsi.gam, center))
     # and add back in residual
     met_adj <- center_predict + e
  }

  # vector or dataframe
  if(!return_df) return(met_adj)

  # otherwise, make a new column
  new_col <- paste0(met_name,'_gamadj')
  MRSI_input[new_col] <- as.numeric(met_adj)
  return(MRSI_input)
}
