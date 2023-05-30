res_with_age <- function(MRSI_input, this_roi, met_name) {
  # have columns we need
  print(names(MRSI_input))
  checkmate::expect_subset(c("roi", "age", "dateNumeric", met_name), names(MRiSI_input))
  
  if (length(this_roi) > 1) {
    model <- paste0(met_name, ' ~ s(age, k=3) + s(dateNumeric, k=5) + s(GMrat, k=3) + label')
  } else {
    model <- paste0(met_name, ' ~ s(age, k=3) + s(dateNumeric, k=5) + s(GMrat, k=3)')
  }
  #print(model)
  mrsi.gam <- gam(as.formula(model), 
                  data=MRSI_input %>% 
                    filter(roi %in% this_roi) %>% 
                    mutate(label = as.factor(label)), na.action = na.exclude)
  
  met_out <- MRSI_input %>% filter(roi %in% this_roi)
  met_out$met_adj <- predict(mrsi.gam, met_out %>% 
                               mutate(dateNumeric = mean(met_out$dateNumeric, na.rm=T), 
                                      GMrat = mean(met_out$GMrat, na.rm=T),
                                      label = as.factor(label))) +
    residuals(mrsi.gam)
  
  
  met_out$met_adj <- as.numeric(met_out$met_adj)
  return(met_out)
}
