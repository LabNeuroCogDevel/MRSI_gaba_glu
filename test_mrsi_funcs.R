library(testthat)
source('mrsi_funcs.R')
met_example <- read.csv('13MP20200207_LCMv2fixidx.csv')

test_that("mrsi qc filter",{
 mrs <- met_example %>% mrsi_metqc()
 expect_lt(nrow(mrs), nrow(met_example))
 expect_lte(max(mrs[c("GPC.Cho.SD", "NAA.NAAG.SD", "Cr.SD")]), 10)
 expect_lte(max(mrs[c("MM20.Cr")]), 3)
})

test_that("mrsi add cols",{
 mrs_in <- met_example %>% filter(!is.na(roi)) %>% mrsi_add_cols()
 mrs <- mrsi_add_cols(mrs_in)
 expect_equal(nrow(mrs), nrow(mrs_in))
 expect_true(all(c("invage","age2","dateNumeric","region") %in% names(mrs)))
 expect_true(all(mrs$age>mrs$invage))

 #age2 is centered
 #expect_true(all(mrs$age<mrs$age2))

 # time in seconds
 expect_gt(diff(range(mrs$dateNumeric)), 10^7)
 expect_lt(diff(range(mrs$dateNumeric)), 10^9)

})
