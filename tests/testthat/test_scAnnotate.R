context("test_scAnnotate")
library(testthat)
library(scAnnotate)

data(pbmc1)
data(pbmc2)
predict_label=scAnnotate(train=pbmc1,
                         test=pbmc2[,-1],
                         distribution="normal",
                         correction ="auto",
                         threshold=0)

test_that("output consistency in default setting",{
  expect_is(predict_label,"character")
  expect_equal(length(predict_label),length(pbmc2[,1]))
  expect_true(all(predict_label%in%c(names(table(pbmc2[,1])),"unassigned")))
})
