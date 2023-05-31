require(dplyr)
require(mgcv)
#
# 20230530 WF copy of FC
# also see MultHemiGamm.R:met_res_for_roi  

#' adjust metabolte with residuals of gam model
#'
#' 20230531 - removed this_roi option. filter before sending to function
#' @param MRSI_input dataframe with roi, age, dateNumeric
#' @param met_name  metabolite eg. 'GABA.Cre'. column in MRSI_input
#' @param return_df bool, if NULL (default) return df only when this_roi is not NULL
#' @return vector of adjusted metabolite OR dataframe with new 'met_name'_gamadj column
res_with_age <- function(MRSI_input, met_name, return_df=FALSE) {
  checkmate::expect_subset(c("age", "dateNumeric", "GMrat", met_name), names(MRSI_input))

  model <- paste0(met_name, ' ~ s(age, k=3) + s(dateNumeric, k=5) + s(GMrat, k=3)')

  # if region and there are multiples, e.g. if have e.g. L and R, model each
  if("region" %in% names(MRSI_input) && length(unique(MRSI_input$region))>1L) {
     regions <- paste0(collapse=',', unique(MRSI_input$region))
     cat(glue::glue("modeling multiple regions: {regions}\n"),"\n")
     MRSI_input$region <- as.factor(MRSI_input$region) # ensure it's a factor for modeling
     model <- paste0(model, '+ region')
  }

  mrsi.gam <- mgcv::gam(as.formula(model), data=MRSI_input, na.action = na.exclude)

  # get residuals at average date and greymatter
  center <- MRSI_input %>%
     mutate(dateNumeric = mean(dateNumeric, na.rm=T),
            GMrat = mean(GMrat, na.rm=T))

  # met_adj is what's left after we've modeled are nuisance regressors
  met_adj <- unname(predict(mrsi.gam,center) + residuals(mrsi.gam))

  # vector or dataframe
  if(!return_df) return(met_adj)

  # otherwise, make a new column
  new_col <- paste0(met_name,'_gamadj')
  MRSI_input[new_col] <- as.numeric(met_adj)
  return(MRSI_input)
}
