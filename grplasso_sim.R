library(grplasso)
data(splice)


contr <- rep(list("contr.sum"), ncol(splice) - 1)

names(contr) <- names(splice)[-1]

fit.splice <- grplasso(y ~ ., data = splice, model = LogReg(), lambda = 20,
contrasts = contr, center = TRUE, standardize = TRUE)


set.seed(79)

n <- 50 ## observations 
p <- 4	## variables

index <- c(NA, 2, 2, 3, 3)

x <- cbind(1, matrix(rnorm(p * n), nrow = n)) 

colnames(x) <- c("Intercept", paste("X", 1:4, sep = ""))

par <- c(0, 2.1, -1.8, 0, 0)

prob <- 1 / (1 + exp(-x %*% par))

mean(pmin(prob, 1 - prob)) ## Bayes risk

y <- rbinom(n, size = 1, prob = prob) ## binary response vector


lambda <- lambdamax(x, y = y, index = index, penscale = sqrt,
model = LogReg()) * 0.5^(0:5)


fit <- grplasso(x, y = y, index = index, lambda = lambda, model = LogReg(),
penscale = sqrt, control = grpl.control(update.hess = "lambda", trace = 0))


plot(fit)