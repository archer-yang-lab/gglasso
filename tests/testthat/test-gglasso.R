library(gglasso)
context("gglasso model fit")

# load bardet data set
data(bardet)

# load colon data set
data(colon)

# define group index
group1 <- rep(1:20, each = 5)
group2 <- rep(1:20, each = 5)

test_that("no error in fitting gglasso for different loss functions", {
  
  fit_ls <- try(gglasso(x = bardet$x, y = bardet$y, group = group1, loss = "ls"),
            silent = TRUE)
  
  fit_logit <- try(gglasso(x = colon$x, y = colon$y, group = group2, loss = "logit"),
               silent = TRUE)
  
  fit_hsvm <- try(gglasso(x = colon$x, y = colon$y, group = group2, loss = "hsvm"),
                   silent = TRUE)
  
  fit_sqsvm <- try(gglasso(x = colon$x, y = colon$y, group = group2, loss = "sqsvm"),
                  silent = TRUE)
  
  expect_false(inherits(fit_ls, "try-error"))
  expect_false(inherits(fit_logit, "try-error"))
  expect_false(inherits(fit_hsvm, "try-error"))
  expect_false(inherits(fit_sqsvm, "try-error"))
  
})