source('res_with_age.R')      # adj_by_model
source('mrsi_funcs.R')        # mrs_wide_to_long_cleaned
loaded_data <- load('mgcv-gam_tibble.Rdata') # gam_models_all, gam_models_18

## checking against dataframe column names in the model
# this is just to show what we eventaully want to do
# will swap 'metdata' for real data
ex_dlpfc_N <- gam_models_all %>%
   filter(biregion=='DLPFC', met=='NAA.NAAG')
ex_dlpfc_model <- dlpfc_N$model[[1]]
ex_dlpfc_data <- dlpfc_N$metdata[[1]]
ex_dlpfc_out <- adj_by_model(this_model, this_data, return_df=TRUE) %>% head()
####


# add region as used by model
mrs_xlx <- readxl::read_excel('caud-thal_MRSI_20241011.xlsx')

## change met names to same as model
# "Glu+Gln %SD"  "Glu+Gln/Cre"
# "Glu.Gln/SD"  "Glu.Gln.Ce"
names(mrs_xlx) <- names(mrs_xlx) %>%
   gsub('/| |\\+','.',.) %>%
   gsub('^-|%','',.) %>%
   gsub('\\.Cre$','.Cr',.)

## reshape data
mrs <- mrs_xlx |>
   mutate(
    # lunaid_d8 is 5 digit id with date like YYYYmmdd
    # RECID is 4 digits. add 0
    ld8=paste0("0", RECID,"_",format(scan_date,"%Y%m%d")),
    roi=case_when(
       roi == 'right caudate' ~ "5",   # 'RCaudate'
       roi == 'left caudate' ~ "6",    # 'LCaudate'
       roi == 'right thalamus' ~ "13", # 'RThalamus'
       roi == 'left thalamus' ~ "99", # LThalamus not in model!,
       .default = NA)) |> 
   select(ld8, RECID,
          age,
          roi,
          scan_date,
          GMrat=GM,
          matches('Cr$|SD$'))

MRSI_input_sp <- mrs |>
   filter(!is.na(roi), roi!="99") |>
   mrsi_add_cols() |>
   mrs_wide_to_long_cleaned()# |>
   #mrsi_metqc()


## subset model
input_summary <- MRSI_input_sp |> count(biregion,met, name="total_visits")
gam_models <-  inner_join(gam_models_18, input_summary) # nrow(gam_models) == 48, down from 128

# confirm we can work with data we have on a single model/region
ex_model <- gam_models    |> filter(region=='LCaudate', met=="NAA.NAAG") |> pull('model')|> first()
ex_data  <- MRSI_input_sp |> filter(region=='LCaudate', met=="NAA.NAAG")
adj_by_model(ex_model, ex_data)

## function to extract input data given a model (region and met)
gam_on_subset<-function(model, d_region, d_met) {
 d <- MRSI_input_sp|>filter(region %in% d_region, met %in% d_met)
 #cat(d_region,"x", d_met, ":", nrow(d),"\n")
 tryCatch(adj_by_model(model, d, return_df=TRUE),error=\(e) NA)
}

## run each model
corrected_nested <- gam_models %>%
   mutate(adj=mapply(gam_on_subset, model, region, met))

corrected <- corrected_nested %>%
   select(adj) %>%
   tidyr::unnest(cols=c(adj))

