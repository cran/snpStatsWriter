library(snpStats)
data(for.exercise, package="snpStats")
test.data <- snps.10[,11:20]

test_that("impute",{
  
