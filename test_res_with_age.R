library(testthat)
source('res_with_age.R')
source('mrsi_funcs.R')
met_example <- read.csv('13MP20200207_LCMv2fixidx.csv') %>% mrsi_metqc() %>% mrsi_add_cols()

test_that("res_with_age return types",{
   res <- res_with_age(met_example, 'Glu.Cr', this_roi=1)
   expect_equal(class(res),'data.frame')
   expect_true(all(res$roi==1))

   res <- res_with_age(met_example, 'Glu.Cr')
   expect_equal(class(res),'array')
   expect_equal(length(res), nrow(met_example))

   res <- res_with_age(met_example, 'Glu.Cr', return_df=TRUE)
   expect_equal(class(res),'data.frame')
})


test_that("res_with_age handle NA",{
   met_example2 <- met_example %>%
      mutate(Glu.Cr = ifelse(Glu.SD>3,NA,Glu.Cr)) %>%
      filter(!is.na(GMrat), !is.na(age))
   nna_glu <- length(which(is.na(met_example2$Glu.Cr)))
   res <- res_with_age(met_example2, 'Glu.Cr')
   nna_glu_m <- length(which(is.na(res)))
   expect_true(nna_glu>10)
   expect_equal(length(res), nrow(met_example2))
   expect_true(all(is.na(res)==is.na(met_example2$Glu.Cr)))
})
