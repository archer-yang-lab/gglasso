\name{print.gglasso}
\alias{print.gglasso}
\title{print a gglasso object}
\description{
Print the nonzero group counts at each lambda step along the gglasso path.
}
\usage{
\method{print}{gglasso}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{fitted \code{\link{gglasso}} object}
  \item{digits}{significant digits in printout}
  \item{\dots}{additional print arguments}
}
\details{
Print the information about the nonzero group counts at each lambda step in the \code{\link{gglasso}} object. The result is a two-column matrix with columns \code{Df} and \code{Lambda}. The \code{Df} column is the number of the groups that have nonzero within-group coefficients at the corresponding \code{Lambda}.
}
\value{
a two-column matrix, the first columns is the number of nonzero group counts and the second column is \code{Lambda}.}

\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "A Fast Unified Algorithm for Computing Group-Lasso Penalized Learning Problems," \emph{Journal of Computational and Graphical Statistics}. Under review.\cr
BugReport: \url{http://code.google.com/p/gglasso/}\cr
}
\examples{
data(FHT)
group = rep(1:50,each=3)
m1 <- gglasso(x=x,y=y_cls,group=group,loss="logit")
print(m1)
}
\keyword{models}
\keyword{regression}

