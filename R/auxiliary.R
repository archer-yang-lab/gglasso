dhsvm <- function(v, delta, weight) {
    r <- v[1]
    if (r > 1) 
        dl <- 0 else if (r <= (1 - delta)) 
        dl <- -1 else dl <- (r - 1)/delta
    dl
}


dlogit <- function(r, delta, weight) {
    dl <- -1/(1 + exp(r))
    dl
}

dsqsvm <- function(r, delta, weight) {
    dl <- -2 * ifelse((1 - r) > 0, (1 - r), 0)
    dl
}

dls <- function(r, delta, weight) {
    dl <- -r
}

dwls <- function(r, delta, weight) {
    dl <- -r%*%weight
}

margin <- function(b0, beta, y, x, delta, loss = c("ls", "logit", 
    "sqsvm", "hsvm", "wls"), weight) {
    loss <- match.arg(loss)
    nobs <- nrow(x)
    b0MAT <- matrix(rep(b0, nobs), nrow = nobs, byrow = TRUE)
    link <- x %*% beta + b0MAT
    if (loss %in% c("logit", "sqsvm", "hsvm")) {
        r <- y * link
    } else r <- y - link
    fun <- paste("d", loss, sep = "")
    if (loss %in% c("wls")) {
    	dMat <- apply(r, 2, eval(fun), delta = delta, weight = weight)
	} else dMat <- apply(r, c(1, 2), eval(fun), delta = delta, weight = weight)
    if (loss %in% c("logit", "sqsvm", "hsvm")) {
        yxdMat <- t(x) %*% (dMat * y)/nobs
    } else if(loss=="wls") yxdMat <- t(x) %*% dMat
	else yxdMat <- t(x) %*% dMat/nobs
    yxdMat
}


KKT <- function(b0, beta, y, x, weight = NULL, lambda, pf, group, thr, delta, 
                loss = c("ls", "logit", "sqsvm", "hsvm","wls")) {
    loss <- match.arg(loss)
	y <- drop(y)
    bn <- as.integer(max(group))
    dl <- margin(b0, beta, y, x, delta, loss, weight)
    B <- matrix(NA, ncol = length(lambda))
    ctr <- 0
    for (l in 1:length(lambda)) {
        for (g in 1:bn) {
            ind <- (group == g)
            dl_norm <- sqrt(crossprod(dl[ind, l], dl[ind, l]))
            b_norm <- sqrt(crossprod(beta[ind, l], beta[ind, l]))
            if (b_norm != 0) {
                AA <- dl[ind, l] + beta[ind, l] * lambda[l] * as.vector(pf[g]/b_norm)
                if (sum(abs(AA)) >= thr) {
                  cat("violate at b != 0", sum(abs(AA)), "\n")
                  ctr <- ctr + 1
                }
            } else {
                BB <- dl_norm - pf[g] * lambda[l]
                if (BB > thr) {
                  cat("violate at b = 0", BB, "\n")
                  ctr <- ctr + 1
                }
            }
        }
    }
    # cat("# of violations", ctr/length(lambda), "\n")
    cat("Averge # of violations per lambda", ctr/length(lambda))
    return(ctr)
} 
