\name{print.gcdnet}
\alias{print.gcdnet}
\title{print a gcdnet object}
\description{
Print a summary of the gcdnet path at each step along the path.
}
\usage{
\method{print}{gcdnet}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{fitted \code{\link{gcdnet}} object}
  \item{digits}{significant digits in printout}
  \item{\dots}{additional print arguments}
}
\details{
The call that produced the \code{\link{gcdnet}} object is printed, followed by a two-column matrix with columns \code{Df} and \code{Lambda}. The \code{Df} column is the number of nonzero coefficients.
}
\value{
a two-column matrix, the first columns is the number of nonzero coefficients and the second column is \code{Lambda}.}

\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "An Efficient Algorithm for Computing The HHSVM and Its Generalizations," \emph{Journal of Computational and Graphical Statistics}. Accepted\cr
BugReport: \url{http://code.google.com/p/gcdnet/}\cr


Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}

\examples{
data(FHT)
m1 <- gcdnet(x=FHT$x,y=FHT$y,delta=1,lambda2=0.1)
print(m1)
}
\keyword{models}
\keyword{regression}

