\name{FHT}
\alias{FHT}
\docType{data}
\title{A simulation model by combining the FHT model introduced in Friedman et al. (2010) and the simulation model 3 in Yuan & Lin (2006).}
\description{
The \code{FHT} data set has n = 100 observations and p = 150 predictors. The covariance between predictors Xj and Xj' has the same correlation 0.5. See details in reference.
}
\usage{data(FHT)}
\format{
  This data contains the following objects:
  \describe{
    \item{x}{a matrix with 100 rows and 150 columns}
    \item{y_cls}{class labels}
    \item{y_reg}{regression response}
  }
}
\references{
Yang, Y. and Zou, H. (2012), "A Fast Unified Algorithm for Computing Group-Lasso Penalized Learning Problems," \emph{Journal of Computational and Graphical Statistics}. Under review.\cr
BugReport: \url{http://code.google.com/p/gglasso/}\cr

Yuan, M. & Lin, Y. (2006), "Model selection and estimation in regression with grouped
variables", \emph{Journal of the Royal Statistical Society}, Series B 68, 49--67.\cr

Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\examples{
data(FHT)
}
\keyword{datasets}
