ls <- function(bn, bs, ix, iy, nobs, nvars, x, y, pf, dfmax, 
    pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr) {
    #################################################################################
    # call Fortran core
    gamma <- rep(NA, bn)
    for (g in 1:bn) gamma[g] <- max(eigen(crossprod(x[, ix[g]:iy[g]]))$values)
    gamma <- gamma/nobs
    gamma <- as.double(gamma)
    fit <- .Fortran("ls_f", bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
        as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
        b0 = double(nlam), beta = double(nvars * nlam), idx = integer(pmax), 
        nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1))
    #################################################################################
    # output
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr, group = group))
    class(outlist) <- c("ls")
    outlist
}

logit <- function(bn, bs, ix, iy, nobs, nvars, x, y, pf, 
    dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr) {
    #################################################################################
    # call Fortran core
    gamma <- rep(NA, bn)
    for (g in 1:bn) gamma[g] <- max(eigen(crossprod(x[, ix[g]:iy[g]]))$values)
    gamma <- 0.25 * gamma/nobs
    gamma <- as.double(gamma)
    fit <- .Fortran("log_f", bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
        as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
        b0 = double(nlam), beta = double(nvars * nlam), idx = integer(pmax), 
        nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1))
    #################################################################################
    # output
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr, group = group))
    class(outlist) <- c("logit")
    outlist
}


hsvm <- function(delta, bn, bs, ix, iy, nobs, nvars, x, y, 
    pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr) {
    #################################################################################
    # call Fortran core
    gamma <- rep(NA, bn)
    for (g in 1:bn) gamma[g] <- max(eigen(crossprod(x[, ix[g]:iy[g]]))$values)
    gamma <- 2 * gamma/(delta * nobs)
    gamma <- as.double(gamma)
    fit <- .Fortran("hsvm_f", delta, bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
        as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
        b0 = double(nlam), beta = double(nvars * nlam), idx = integer(pmax), 
        nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1))
    #################################################################################
    # output
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr, group = group))
    class(outlist) <- c("hsvm")
    outlist
}

sqsvm <- function(bn, bs, ix, iy, nobs, nvars, x, y, pf, 
    dfmax, pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr) {
    #################################################################################
    # call Fortran core
    gamma <- rep(NA, bn)
    for (g in 1:bn) gamma[g] <- max(eigen(crossprod(x[, ix[g]:iy[g]]))$values)
    gamma <- 4 * gamma/nobs
    gamma <- as.double(gamma)
    fit <- .Fortran("sqsvm_f", bn, bs, ix, iy, gamma, nobs, nvars, as.double(x), 
        as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
        b0 = double(nlam), beta = double(nvars * nlam), idx = integer(pmax), 
        nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1))
    #################################################################################
    # output
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr, group = group))
    class(outlist) <- c("sqsvm")
    outlist
} 


wls <- function(bn, bs, ix, iy, nobs, nvars, x, y, pf, weight, dfmax, 
    pmax, nlam, flmin, ulam, eps, maxit, vnames, group, intr) {
    #################################################################################
    # call Fortran core
    gamma <- rep(NA, bn)
	wx <- weight %*% x
    for (g in 1:bn) gamma[g] <- max(eigen(crossprod(x[, ix[g]:iy[g]], wx[, ix[g]:iy[g]]))$values)
    gamma <- as.double(gamma)
    fit <- .Fortran("wls_f", bn, bs, ix, iy, as.double(weight), gamma, nobs, nvars, as.double(x), 
        as.double(y), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, intr, nalam = integer(1), 
        b0 = double(nlam), beta = double(nvars * nlam), idx = integer(pmax), 
        nbeta = integer(nlam), alam = double(nlam), npass = integer(1), jerr = integer(1))
    #################################################################################
    # output
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr, group = group))
    class(outlist) <- c("ls")
    outlist
}
