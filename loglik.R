loglik <- function(x,y,beta,b0)
  {
    L <- ncol(beta)
    val <- numeric(L)
    for (l in 1:L)
    {
       eta <- x%*%beta[,l]+b0[l]
       pi. <- exp(eta)/(1+exp(eta))
       val[l] <- -sum(y*log(pi.)+(1-y)*log(1-pi.))
    }
    return(val)
  }
