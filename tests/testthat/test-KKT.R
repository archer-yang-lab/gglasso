library(gglasso)
context("KKT checks")

# load bardet data set
data(bardet)

# load colon data set
data(colon)

# define group index
group1 <- rep(1:20, each = 5)
group2 <- rep(1:20, each = 5)

test_that("KKT checks for loss = 'ls', without observation weights", {

  # fit group lasso penalized least squares
  m1 <- gglasso(x = bardet$x, y = bardet$y, group = group1, loss = "ls")
  violations <- gglasso:::KKT(b0 = m1$b0, beta = m1$beta, y = bardet$y, x = bardet$x, lambda = m1$lambda,
                              pf = rep(sqrt(5),20), group = group1, thr = 1e-3, loss = "ls")
  
  expect_lte(violations, 0)
  
})


test_that("KKT checks for loss = 'logit', without observation weights", {
  
  # fit group lasso penalized logistic regression
  m2 <- gglasso(x = colon$x, y = colon$y, group = group2, loss="logit")

  violations <- gglasso:::KKT(b0 = m2$b0, beta = m2$beta, y = colon$y, x = colon$x, lambda = m2$lambda,
                              pf = rep(sqrt(5),20), group = group2, thr = 1e-3, loss = "logit")
  
  expect_lte(violations, 0)
  
})

# View(gglasso:::KKT)
# 
